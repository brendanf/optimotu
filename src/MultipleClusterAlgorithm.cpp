#include <cassert>
#include <algorithm>
#include <string>
#include "MultipleClusterAlgorithm.h"
#include "optimotu.h"

using MCA = MultipleClusterAlgorithm;

std::vector<std::vector<j_t>> calculate_subset_indices(
  const std::vector<std::string> &names,
  const std::vector<std::vector<std::string>> &subset_names
) {
  std::unordered_map<std::string, j_t> namekey;
  for (j_t i = 0; i < names.size(); i++) {
    namekey.emplace(names[i], i);
  }
  std::vector<std::vector<j_t>> out;
  out.reserve(subset_names.size());
  for (const auto & my_subset_names : subset_names) {
    std::vector<j_t> subset_indices;
    subset_indices.reserve(my_subset_names.size());
    for (const auto & name : my_subset_names) {
      auto f = namekey.find(name);
      if (f == namekey.end()) {
        OPTIMOTU_STOP("name '" + name +
           "' in subset is not found in the names of the sequences");
      }
      subset_indices.push_back(f->second);
    }
    std::sort(subset_indices.begin(), subset_indices.end());
    out.push_back(subset_indices);
  }
  return out;
}

// Base protected main constructor: initializer list only.
MCA::MultipleClusterAlgorithm(
  const ClusterAlgorithmFactory & factory,
  const std::vector<std::string> &names,
  const std::vector<std::vector<std::string>> &subset_names,
  const int threads
) :
  ClusterAlgorithm(factory.dconv),
  factory(factory),
  names(names),
  subset_indices(calculate_subset_indices(names, subset_names)),
  subset_names(subset_names),
  threads(threads),
  subset_key(names.size()),
  fwd_map(subset_names.size()),
  subsets(),
  whichsets(threads),
  ws_keys(threads)
{}

template<int verbose>
MultipleClusterAlgorithmImpl<verbose>::MultipleClusterAlgorithmImpl(MultipleClusterAlgorithm * parent) :
  MultipleClusterAlgorithm(parent) {}

template<int verbose>
MultipleClusterAlgorithmImpl<verbose>::MultipleClusterAlgorithmImpl(
  MultipleClusterAlgorithm * parent,
  const std::vector<std::string> names,
  const std::vector<std::vector<j_t>> subset_indices,
  const std::vector<std::vector<std::string>> subset_names,
  const std::vector<std::vector<j_t>> subset_key,
  const std::vector<std::unordered_map<j_t, j_t>> fwd_map,
  const std::vector<j_t> child_to_parent_map,
  PairGenerator * pg,
  const int threads
) : MultipleClusterAlgorithm(parent, names, subset_indices, subset_names, subset_key, fwd_map, child_to_parent_map, pg, threads) {}

template<int verbose>
MultipleClusterAlgorithmImpl<verbose>::MultipleClusterAlgorithmImpl(
  const ClusterAlgorithmFactory & factory,
  const std::vector<std::string> &names,
  const std::vector<std::vector<std::string>> &subset_names,
  const int threads,
  int
) :
  MultipleClusterAlgorithm(factory, names, subset_names, threads)
{
  OPTIMOTU_DEBUG(
    1,
    << "  Allocating " << subset_names.size() << " subsets..." << std::flush;
  );
  subsets.reserve(subset_names.size());
  for (auto & ws : whichsets)
    ws.reserve(subset_names.size());

  OPTIMOTU_DEBUG(
    1,
    << "done\n  Generating sequence name key for " << names.size()
    << " sequence names..." << std::flush;
  );
  std::unordered_map<std::string, j_t> namekey;
  for (j_t i = 0; i < names.size(); i++) {
    namekey.emplace(names[i], i);
  }

  OPTIMOTU_DEBUG(
    1,
    << "done\n  Mapping " << subset_names.size()
    << " subsets..." << std::flush;
  );

  OPTIMOTU_DEBUG(
    2,
    << std::endl;
  );

  for (j_t i = 0; i < subset_names.size(); ++i) {
    owned_subsets.emplace_back(factory.create(subset_names[i].size(), verbose));
    subsets.push_back(owned_subsets.back().get());
    fwd_map[i].reserve(subset_names[i].size());
    OPTIMOTU_DEBUG(
      4,
      << "  Subset " << i << ":" << std::endl;
    );
    for (j_t j = 0; j < subset_names[i].size(); ++j) {
      auto f = namekey.find(subset_names[i][j]);
      assert(f != namekey.end());
      OPTIMOTU_DEBUG(
        4,
        << "    adding sequence " << f->second
        << " (" << subset_names[i][j]
        << ") to subset " << i
        << " at position " << fwd_map[i].size()
        << std::endl;
      );
      subset_key[f->second].push_back(i);
      fwd_map[i].emplace(f->second, j);
      OPTIMOTU_DEBUG(
        4,
        << "    sequence " << f->second
        << " now found in " << subset_key[j].size()
        << " subsets\n   precluster " << i
        << " now has " << fwd_map[i].size()
        << " sequences" << std::endl;
      );
    }
  }
}

// // Subset the names of the sequences to only the names of the sequences
// // which can be reached by the PairGenerator.
// std::vector<std::string> child_names(
//   const std::vector<std::string> & names,
//   const PairGenerator * const pg
// ) {
//   std::vector<std::string> out;
//   out.reserve(pg->max_value());
//   for (size_t i = 0; i < pg->max_value(); i++) {
//     out.push_back(names[pg->forward_map(i)]);
//   }
//   return out;
// }

// std::vector<std::vector<j_t>> child_subset_key(
//   const std::vector<std::vector<j_t>> & parent_subset_key,
//   const std::size_t parent_n_subsets,
//   const PairGenerator * const pg
// ) {
//   std::vector<std::vector<j_t>> out(parent_n_subsets);
//   for (std::size_t i = 0; i < pg->max_value(); i++) {
//     std::size_t i0 = pg->forward_map(i);
//   }
//   return out;
// }

// // Subset the indices of the sequences to only those which can be reached by
// // the PairGenerator.  This results in some empty subsets!
// std::vector<std::vector<j_t>> child_subset_indices(
//   const std::vector<std::vector<j_t>> & parent_subset_key,
//   const std::size_t parent_n_subsets,
//   const PairGenerator * const pg
// ) {
//   std::vector<std::vector<j_t>> out(parent_n_subsets);
//   // For each sequence that the PairGenerator can produce,
//   // find the subsets that it belongs to in the parent
//   // and add the mapped index to each of those subsets.
//   for (std::size_t i = 0; i < pg->max_value(); i++) {
//     std::size_t i0 = pg->forward_map(i);
//     for (const j_t j : parent_subset_key[i0]) {
//       out[j].push_back(i);
//     }
//   }
//   return out;
// }

// // Find the names of the sequences in each subset.
// std::vector<std::vector<std::string>> child_subset_names(
//   const std::vector<std::string> & child_names,
//   const std::vector<std::vector<j_t>> & child_subset_indices
// ) {
//   std::vector<std::vector<std::string>> out;
//   out.reserve(child_subset_indices.size());
//   for (const auto & my_subset : child_subset_indices) {
//     std::vector<std::string> my_names;
//     my_names.reserve(my_subset.size());
//     for (const auto & i : my_subset) {
//       my_names.push_back(child_names[i]);
//     }
//     out.push_back(my_names);
//   }
//   return out;
// }

MCA::MultipleClusterAlgorithm(MCA * parent) :
  ClusterAlgorithm(parent),
  factory(parent->factory),
  names(parent->names),
  subset_indices(parent->subset_indices),
  subset_names(parent->subset_names),
  threads(parent->threads),
  subset_key(parent->subset_key),
  fwd_map(parent->fwd_map),
  whichsets(parent->threads),
  ws_keys(parent->threads)
{
  subsets.reserve(subset_names.size());
  for (auto & subset : parent->subsets) {
    subsets.push_back(subset->make_child());
  }
  for (auto & ws : whichsets)
    ws.reserve(subset_names.size());
}

// This constructor is used to create a child algorithm for a
// MultipleClusterAlgorithm.
MCA::MultipleClusterAlgorithm(
  MCA * parent,
  const std::vector<std::string> names,
  const std::vector<std::vector<j_t>> subset_indices,
  const std::vector<std::vector<std::string>> subset_names,
  const std::vector<std::vector<j_t>> subset_key,
  const std::vector<std::unordered_map<j_t, j_t>> fwd_map,
  const std::vector<j_t> child_to_parent_map,
  PairGenerator * pg,
  const int threads
) :
  ClusterAlgorithm{parent},
  factory{parent->factory},
  names{names},
  subset_indices{subset_indices},
  subset_names{subset_names},
  threads{threads},
  subset_key{subset_key},
  fwd_map{fwd_map},
  whichsets{threads},
  ws_keys{threads}
  {
    subsets.reserve(child_to_parent_map.size());
    for (const j_t j0 : child_to_parent_map) {
      subsets.push_back(parent->subsets[j0]->make_child(pg));
    }
    for (auto & ws : whichsets)
      ws.reserve(fwd_map.size());
}

// does not anything directly!
// relies on locks inside subset algorithms for thread safety.
void MCA::operator()(j_t seq1, j_t seq2, double dist, int thread) {
  if (seq1 == seq2) return;
  d_t i = dconv.convert(dist);
  (*this)(seq1, seq2, i, thread);
}

// does not lock anything directly!
// relies on locks inside subset algorithms for thread safety.
void MCA::operator()(j_t seq1, j_t seq2, int i, int thread) {
  // in practice, we should always be able to rely on the 'whichsets'
  // which was calculated in "max_relevant", if it was called.

  std::unique_lock<std::shared_timed_mutex> lock(this->mutex);
  // OPTIMOTU_COUT << "in operator() for MultipleClusterAlgorithm thread " << thread
  //           << std::endl << "length of whichsets is " << whichsets.size()
  //           << std::endl << "length of ws_keys is " << ws_keys.size()
  //           << std::endl << "finding overlap sets for item " << seq1
  //           << " in sets:";
  // for (auto pc : subset_key[seq1]){
  //   OPTIMOTU_COUT << " " << pc;
  // }
  // OPTIMOTU_COUT << std::endl
  //           << "and item " << seq2 << " in sets:";
  // for (auto pc : subset_key[seq2]) OPTIMOTU_COUT << " " << pc;
  // OPTIMOTU_COUT << std::endl
  //           << "cached overlaps are for pair " << ws_keys[thread].first
  //           << ", " << ws_keys[thread].second
  //           << std::endl;
  if (seq1 == seq2) return;
  if (seq1 >= subset_key.size() || seq2 >= subset_key.size()) {
    OPTIMOTU_STOP("MultipleClusterAlgorithm: sequence index out of range (seq1="
      + std::to_string(seq1) + ", seq2=" + std::to_string(seq2)
      + ", names.size()=" + std::to_string(subset_key.size()) + ")");
  }

  if (ws_keys[thread] != std::pair<j_t, j_t>{seq1, seq2}) {
    // OPTIMOTU_COUT << "calculating..." << std::flush;
    whichsets[thread].clear();
    std::set_intersection(
      subset_key[seq1].begin(),
      subset_key[seq1].end(),
      subset_key[seq2].begin(),
      subset_key[seq2].end(),
      std::back_inserter(whichsets[thread])
    );
    ws_keys[thread] = {seq1, seq2};
  } else {
    // OPTIMOTU_COUT << "using cached values..." << std::flush;
  }

  for (j_t j : whichsets[thread]) {
    // if (j == 0) {OPTIMOTU_COUT << "in operator() for MultipleClusterAlgorithm thread " << thread
    //   // << std::endl << "length of whichsets is " << whichsets.size()
    //   // << std::endl << "length of ws_keys is " << ws_keys.size()
    //      << std::endl << "finding overlap sets for item " << seq1
    //      << " in sets:";
    //   for (auto pc : subset_key[seq1]) OPTIMOTU_COUT << " " << pc;
    //   OPTIMOTU_COUT << std::endl
    //             << "and seq " << seq2 << " in sets:";
    //   for (auto pc : subset_key[seq2]) OPTIMOTU_COUT << " " << pc;
    //   OPTIMOTU_COUT << std::endl
    //             << "cached overlaps are for pair " << ws_keys[thread].first
    //             << ", " << ws_keys[thread].second
    //             << std::endl;
    //   OPTIMOTU_COUT << "found " << whichsets[thread].size()
    //             << " overlaps:";
    //   for (auto pc : whichsets[thread]) OPTIMOTU_COUT << " " << pc;
    //   OPTIMOTU_COUT << std::endl;
    // }
    // OPTIMOTU_COUT << "sending seq1=" << fwd_map[j][seq1]
    //               << " seq2=" << fwd_map[j][seq2]
    //               << " i=" << i
    //               << " to subset " << j
    //               << std::endl;
    (*subsets[j])(fwd_map[j].at(seq1), fwd_map[j].at(seq2), i, thread);
  }
}

void MCA::finalize() {
  for (auto s : subsets) s->finalize();
}

// does not lock anything directly!
// relies on locks inside subset algorithms for thread safety.
double MCA::max_relevant(j_t seq1, j_t seq2, int thread) const {
  if (seq1 >= subset_key.size() || seq2 >= subset_key.size()) {
    OPTIMOTU_STOP("MultipleClusterAlgorithm::max_relevant: sequence index out of range");
  }
  double max = -1.0;
  whichsets[thread].clear();
  std::set_intersection(
    subset_key[seq1].begin(),
    subset_key[seq1].end(),
    subset_key[seq2].begin(),
    subset_key[seq2].end(),
    std::back_inserter(whichsets[thread])
  );
  ws_keys[thread] = {seq1, seq2};
  for (j_t j : whichsets[thread]) {
    double maxj = subsets[j]->max_relevant(fwd_map[j].at(seq1), fwd_map[j].at(seq2));
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
  std::unique_lock<std::shared_timed_mutex> lock(this->mutex);
  auto child_ptr = new MultipleClusterAlgorithm(this);
  auto child = std::unique_ptr<ClusterAlgorithm>(
    (ClusterAlgorithm*)child_ptr
  );
  this->children.push_back(std::move(child));
  return child_ptr;
}

template<int verbose>
MultipleClusterAlgorithm * MultipleClusterAlgorithmImpl<verbose>::make_child() {
  std::unique_lock<std::shared_timed_mutex> lock(this->mutex);
  auto child_ptr = new MultipleClusterAlgorithmImpl<verbose>(this);
  auto child = std::unique_ptr<ClusterAlgorithm>(
    (ClusterAlgorithm*)child_ptr
  );
  this->children.push_back(std::move(child));
  return child_ptr;
}

MultipleClusterAlgorithm * MCA::make_child(PairGenerator * pg) {
  std::unique_lock<std::shared_timed_mutex> lock(this->mutex);

  // calculate which subsets have members that the PairGenerator
  // can produce
  std::vector<bool> nonempty(subset_indices.size());
  std::size_t n_nonempty = 0;
  for (std::size_t i = 0; i <= pg->max_value(); i++) {
    std::size_t i0 = pg->forward_map(i);
    for (const j_t j : subset_key[i0]) {
      if (!nonempty[j]) {
        nonempty[j] = true;
        n_nonempty++;
      }
    }
  }

  // map from our subset to the child's subsets
  std::unordered_map<j_t, j_t> parent_to_child_map;
  parent_to_child_map.reserve(n_nonempty);
  // map from the child's subsets to our subsets
  std::vector<j_t> child_to_parent_map;
  child_to_parent_map.reserve(n_nonempty);
  std::size_t child_index = 0;
  for (std::size_t i = 0; i < nonempty.size(); i++) {
    if (nonempty[i]) {
      parent_to_child_map.emplace(i, child_index);
      child_to_parent_map.push_back(i);
      child_index++;
    }
  }

  // Now calculate the child's subset_key, subset_indices, subset_names,
  // and fwd_map
  std::vector<std::vector<j_t>> child_subset_key(pg->max_value() + 1);
  std::vector<std::vector<j_t>> child_subset_indices(n_nonempty);
  std::vector<std::vector<std::string>> child_subset_names(n_nonempty);
  std::vector<std::unordered_map<j_t, j_t>> child_fwd_map(n_nonempty);
  // Loop over each sequence that the PairGenerator can produce.
  for (std::size_t i = 0; i <= pg->max_value(); i++) {
    std::size_t i0 = pg->forward_map(i);
    // Loop over the subsets that the sequence belongs to in the parent.
    for (const j_t j0 : subset_key[i0]) {
      j_t j = parent_to_child_map[j0];
      child_subset_key[i].push_back(j);
      child_subset_indices[j].push_back(i);
      child_subset_names[j].push_back(names[i]);
      child_fwd_map[j].emplace(i0, i);
    }
  }

  auto child_ptr = new MultipleClusterAlgorithm(
    this,
    names,
    child_subset_indices,
    child_subset_names,
    child_subset_key,
    child_fwd_map,
    child_to_parent_map,
    pg,
    threads
  );
  auto child = std::unique_ptr<ClusterAlgorithm>(
    (ClusterAlgorithm*)child_ptr
  );
  this->children.push_back(std::move(child));
  return child_ptr;
}

template<int verbose>
MultipleClusterAlgorithm * MultipleClusterAlgorithmImpl<verbose>::make_child(PairGenerator * pg) {
  std::unique_lock<std::shared_timed_mutex> lock(this->mutex);

  std::vector<bool> nonempty(subset_indices.size());
  std::size_t n_nonempty = 0;
  for (std::size_t i = 0; i <= pg->max_value(); i++) {
    std::size_t i0 = pg->forward_map(i);
    for (const j_t j : subset_key[i0]) {
      if (!nonempty[j]) {
        nonempty[j] = true;
        n_nonempty++;
      }
    }
  }

  std::unordered_map<j_t, j_t> parent_to_child_map;
  parent_to_child_map.reserve(n_nonempty);
  std::vector<j_t> child_to_parent_map;
  child_to_parent_map.reserve(n_nonempty);
  std::size_t child_index = 0;
  for (std::size_t i = 0; i < nonempty.size(); i++) {
    if (nonempty[i]) {
      parent_to_child_map.emplace(i, child_index);
      child_to_parent_map.push_back(i);
      child_index++;
    }
  }

  std::vector<std::vector<j_t>> child_subset_key(pg->max_value() + 1);
  std::vector<std::vector<j_t>> child_subset_indices(n_nonempty);
  std::vector<std::vector<std::string>> child_subset_names(n_nonempty);
  std::vector<std::unordered_map<j_t, j_t>> child_fwd_map(n_nonempty);
  for (std::size_t i = 0; i <= pg->max_value(); i++) {
    std::size_t i0 = pg->forward_map(i);
    for (const j_t j0 : subset_key[i0]) {
      j_t j = parent_to_child_map[j0];
      child_subset_key[i].push_back(j);
      child_subset_indices[j].push_back(i);
      child_subset_names[j].push_back(names[i]);
      child_fwd_map[j].emplace(i0, i);
    }
  }

  auto child_ptr = new MultipleClusterAlgorithmImpl<verbose>(
    this,
    names,
    child_subset_indices,
    child_subset_names,
    child_subset_key,
    child_fwd_map,
    child_to_parent_map,
    pg,
    threads
  );
  auto child = std::unique_ptr<ClusterAlgorithm>(
    (ClusterAlgorithm*)child_ptr
  );
  this->children.push_back(std::move(child));
  return child_ptr;
  // this->own_child = true;
  // for (std::size_t i = 0; i < subsets.size(); i++) {
  //   bool found_match = false;
  //   for (std::size_t j = 0; j < pg->max_value(); j++) {
  //     if (fwd_map[i].count(pg->forward_map(j)) > 0) {
  //       found_match = true;
  //       break;
  //     }
  //   }
  //   if (found_match) {
  //     subsets[i]->make_child();
  //   }
  //     ss->make_child();
  // return this;
}

void MCA::write_to_matrix(std::vector<internal_matrix_t> &matrix_list) {
  for (size_t i = 0; i < this->subsets.size(); i++) {
    // OPTIMOTU_COUT << "writing to matrix " << i << "..." << std::flush;
    //           << "matrix size is " << matrix_list[i].nrow()
    //           << "x" << matrix_list[i].ncol()
    //           << std::endl << "subset size is " << this->m
    //           << "x" << subsets[i]->n
    //           << std::endl;
    this->subsets[i]->write_to_matrix(matrix_list[i]);
    // OPTIMOTU_COUT << "done" << std::endl;
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

template class MultipleClusterAlgorithmImpl<0>;
template class MultipleClusterAlgorithmImpl<1>;
template class MultipleClusterAlgorithmImpl<2>;
template class MultipleClusterAlgorithmImpl<3>;
template class MultipleClusterAlgorithmImpl<4>;
