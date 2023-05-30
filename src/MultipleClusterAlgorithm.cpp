#include <cassert>
#include <algorithm>
#include "MultipleClusterAlgorithm.h"

using MCA = MultipleClusterAlgorithm;

MCA::MultipleClusterAlgorithm(
  const ClusterAlgorithmFactory & factory,
  const std::vector<std::string> &names,
  const std::vector<std::vector<std::string>> &subset_names,
  const int threads
) :
  ClusterAlgorithm(factory.dconv),
  factory(factory),
  names(names),
  subset_names(subset_names),
  threads(threads),
  subsets(),
  subset_key(names.size()),
  fwd_map(subset_names.size()),
  whichsets(threads)
{
  subsets.reserve(subset_names.size());
  for (auto & ws : whichsets)
    ws.reserve(subset_names.size());
  std::unordered_map<std::string, j_t> namekey;
  for (j_t i = 0; i < names.size(); i++) {
    namekey.emplace(names[i], i);
  }
  for (j_t i = 0; i < subset_names.size(); ++i) {
    subsets.push_back(factory.create(subset_names[i].size()));
    fwd_map[i].reserve(subset_names[i].size());
    for (j_t j = 0; j < subset_names[i].size(); ++j) {
      auto f = namekey.find(subset_names[i][j]);
      assert(f != namekey.end());
      subset_key[f->second].push_back(i);
      fwd_map[i].emplace(f->second, j);
    }
  }
}

MCA::MultipleClusterAlgorithm(MCA * parent) :
  ClusterAlgorithm{parent},
  factory{parent->factory},
  names{parent->names},
  subset_names{parent->subset_names},
  threads{parent->threads},
  subsets{},
  subset_key(parent->subset_key),
  fwd_map(parent->fwd_map),
  whichsets(parent->threads)
  {
    subsets.reserve(subset_names.size());
    for (auto & ss : parent->subsets) {
      subsets.emplace_back(std::move(ss->make_child()));
    }
    for (auto & ws : whichsets)
      ws.reserve(subset_names.size());
}

// does not anything directly!
// relies on locks inside subset algorithms for thread safety.
void MCA::operator()(j_t seq1, j_t seq2, double dist, int thread) {
  d_t i = dconv.convert(dist);
  (*this)(seq1, seq2, i, thread);
}

// does not lock anything directly!
// relies on locks inside subset algorithms for thread safety.
void MCA::operator()(j_t seq1, j_t seq2, int i, int thread) {
  // in practice, we should always be able to rely on the 'whichsets'
  // which was calculated in "max_relevant", if it was called.
  if (ws_keys[i] != std::pair<j_t, j_t>{seq1, seq2}) {
    whichsets[thread].clear();
    std::set_intersection(
      subset_key[seq1].begin(),
      subset_key[seq1].end(),
      subset_key[seq2].begin(),
      subset_key[seq2].end(),
      whichsets[thread].begin()
    );
    ws_keys[thread] = {seq1, seq2};
  }
  for (j_t j : whichsets[thread]) {
    (*subsets[j])(fwd_map[j][seq1], fwd_map[j][seq2], i);
  }
}

// does not lock anything directly!
// relies on locks inside subset algorithms for thread safety.
double MCA::max_relevant(j_t seq1, j_t seq2, int thread) const {
  double max = -1.0;
  whichsets[thread].clear();
  std::set_intersection(
    subset_key[seq1].begin(),
    subset_key[seq1].end(),
    subset_key[seq2].begin(),
    subset_key[seq2].end(),
    whichsets[thread].begin()
  );
  ws_keys[thread] = {seq1, seq2};
  for (j_t j : whichsets[thread]) {
    double maxj = subsets[j]->max_relevant(seq1, seq2);
    max = std::max(max, maxj);
  }
  return max;
}

void MCA::merge_into(DistanceConsumer &consumer) {
  for (auto & ss : this->subsets) {
    ss->merge_into(consumer);
  }
}

void MCA::merge_into(ClusterAlgorithm &consumer) {
  for (auto & ss : this->subsets) {
    ss->merge_into(consumer);
  }
}

// send consumer() pairwise distances to ensure it is up-to-date with this
// clustering
void MCA::merge_into(MCA &consumer) {
  for (size_t i = 0; i < this->subsets.size(); i++) {
    this->subsets[i]->merge_into(*consumer.subsets[i]);
  }
}

void MCA::merge_into_parent() {
  for (auto & ss : this->subsets) {
    ss->merge_into_parent();
  }
}

MultipleClusterAlgorithm * MCA::make_child() {
  tbb::queuing_rw_mutex::scoped_lock lock(this->mutex);
  if (own_child) {
    auto child_ptr = new MultipleClusterAlgorithm(this);
    auto child = std::unique_ptr<ClusterAlgorithm>(
      (ClusterAlgorithm*)child_ptr
    );
    this->children.push_back(std::move(child));
    return child_ptr;
  }
  this->own_child = true;
  return this;
}

void MCA::write_to_matrix(std::vector<internal_matrix_t> &matrix_list) {
  for (size_t i = 0; i < this->subsets.size(); i++) {
    this->subsets[i]->write_to_matrix(matrix_list[i]);
  }
}

#ifdef OPTIMOTU_R
Rcpp::List MCA::as_hclust() const {
  std::vector<Rcpp::List> out;
  for (size_t i = 0; i < this->subsets.size(); i++) {
    out.push_back(this->subsets[i]->as_hclust(Rcpp::wrap(this->subset_names[i])));
  }
  return Rcpp::wrap(out);
}
Rcpp::List MCA::as_hclust(const Rcpp::CharacterVector &seqnames) const {
  return this->as_hclust();
}
#endif //OPTIMOTU_R
