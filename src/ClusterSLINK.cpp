#include "ClusterSLINK.h"
#include <cstdint>

ClusterSLINK::ClusterSLINK(SingleClusterAlgorithm * parent) :
  SingleClusterAlgorithm(parent), Pi(n), Lambda(n, m), M(n, m), delegate(dconv, n) {}

ClusterSLINK::ClusterSLINK(const DistanceConverter &dconv, const j_t n) :
  SingleClusterAlgorithm(dconv, n),
  Pi(n),
  Lambda(n, m),
  M(n, m),
  delegate(dconv, n) {}

ClusterSLINK::ClusterSLINK(const DistanceConverter &dconv, init_matrix_t im) :
  SingleClusterAlgorithm(dconv, im),
  Pi(n),
  Lambda(n, m),
  M(n, m),
  delegate(dconv, im) {}

void ClusterSLINK::init_iter() {
  // OPTIMOTU_CERR << "### Initializing SLINK iteration " << slink_seq1
  //               << std::endl;
  for (j_t k = 0; k < slink_seq1; k++) {
    M[k] = m;
  }
  slink_seq2 = 0;
}

void ClusterSLINK::update() {
  if (Lambda[slink_seq2] >= M[slink_seq2]) {
    // OPTIMOTU_CERR << "Found closer (or equal) match for sequences " << slink_seq1
    //               << " and " << slink_seq2
    //               << ": " << Lambda[slink_seq2]
    //               << " >= " << M[slink_seq2]
    //               << std::endl;
    // OPTIMOTU_CERR << "Setting M[" << Pi[slink_seq2]
    //               << "] from " << M[Pi[slink_seq2]];
    M[Pi[slink_seq2]] = std::min(M[Pi[slink_seq2]], Lambda[slink_seq2]);
    // OPTIMOTU_CERR <<  " to " << M[Pi[slink_seq2]]
    //               << "\nSetting Lambda[" << slink_seq2
    //               << "] from " << Lambda[slink_seq2];
    Lambda[slink_seq2] = M[slink_seq2];
    // OPTIMOTU_CERR <<  " to " << Lambda[slink_seq2]
    //               << "\nSetting Pi[" << slink_seq2
    //               << "] from " << Pi[slink_seq2];
    Pi[slink_seq2] = slink_seq1;
    // OPTIMOTU_CERR <<  " to " << Pi[slink_seq2] << std::endl;
  } else {
    // OPTIMOTU_CERR << "larger distance for sequences " << slink_seq1
    //               << " and " << slink_seq2
    //               << ": " << Lambda[slink_seq2]
    //               << " < " << M[slink_seq2]
    //               << std::endl;
    // OPTIMOTU_CERR << "Setting M[" << Pi[slink_seq2]
    //               << "] from " << M[Pi[slink_seq2]];
    M[Pi[slink_seq2]] = std::min(M[Pi[slink_seq2]], M[slink_seq2]);
    // OPTIMOTU_CERR <<  " to " << M[Pi[slink_seq2]] << std::endl;
  }
}

void ClusterSLINK::finish_iter() {
  // OPTIMOTU_CERR << "### Finishing SLINK iteration " << slink_seq1
  //               << std::endl;
  while (slink_seq2 < slink_seq1) {
    update();
    slink_seq2++;
  }
  for (j_t k = 0; k < slink_seq1; k++) {
    if (Lambda[k] >= Lambda[Pi[k]]) Pi[k] = slink_seq1;
  }
}

void ClusterSLINK::finalize() {
  finish_iter();
}

ClusterSLINK * ClusterSLINK::make_child() {
  std::unique_lock<std::shared_timed_mutex> lock(this->mutex);
  if (own_child) {
    auto child_ptr = new ClusterSLINK(&delegate);
    auto child = std::unique_ptr<ClusterAlgorithm>(
      (ClusterAlgorithm*)child_ptr
    );
    this->children.push_back(std::move(child));
    return child_ptr;
  }
  this->own_child = true;
  return this;
}

void ClusterSLINK::operator()(j_t seq2, j_t seq1, d_t i, int thread) {
  if (seq1 <= seq2) {
    OPTIMOTU_CERR << "ClusterSLINK requires seq2 < seq1.  seq2=" << seq2
                  << ", seq1=" << seq1
                  << ", thread " << thread << std::endl;
    OPTIMOTU_STOP("ClusterSLINK input error.");
  }
  if (seq1 < slink_seq1) {
    OPTIMOTU_CERR << "ClusterSLINK requires sequences in order."
                  << " Current seq1=" << seq1
                  << ", slink_seq1=" << slink_seq1
                  << ", thread " << thread << std::endl;
    OPTIMOTU_STOP("ClusterSLINK input error.");
  }
  while (seq1 > slink_seq1) {
    finish_iter();
    ++slink_seq1;
    Pi[slink_seq1] = slink_seq1;
    init_iter();
  }
  if (seq2 < slink_seq2) {
    OPTIMOTU_CERR << "ClusterSLINK requires sequences in order."
                  << " Current seq2=" << seq2
                  << ", slink_seq2=" << slink_seq2
                  << ", thread " << thread << std::endl;
    OPTIMOTU_STOP("ClusterSLINK input error.");
  }
  while (slink_seq2 <= seq2) {
    if (slink_seq2 == seq2) {

      // OPTIMOTU_CERR << "Received distance " << i
      //               << " for seq " << slink_seq2
      //               << "; old value was " << M[slink_seq2];
      M[slink_seq2] = std::min(i, M[slink_seq2]);
      // OPTIMOTU_CERR << "; new value is " << M[slink_seq2] << std::endl;
    }
    update();
    slink_seq2++;
  }
}

void ClusterSLINK::write_to_matrix(internal_matrix_t &out) {
  if (own_child && children.size() > 0) return delegate.write_to_matrix(out);
  std::shared_lock<std::shared_timed_mutex> lock(this->mutex);
  // OPTIMOTU_CERR << "preparing to write matrix" << std::endl
  //               << " i Pi Lambda" << std::endl;
  // for (j_t i = 0; i < this->n; i++) {
  //   OPTIMOTU_CERR << std::setw(2) << i
  //                 << std::setw(3) << Pi[i]
  //                 << std::setw(7) << Lambda[i] << std::endl;
  // }
  j_t j, jp;
  std::size_t k = 0;
  for (std::uint32_t i = 0; i < this->n; i++) {
    j = i;
    jp = Pi[j];
    j_t i2 = 0;
    while (i2 < this->m) {
      d_t max = Lambda[j];
      if (this->m < (size_t)max) max = this->m;
      while (i2 < (size_t)max) {
        out[k++] = j;
        i2++;
      }
      if (i2 < this->m) {
        j = Pi[j];
        jp = Pi[j];
      }
    }
  }
}

#ifdef OPTIMOTU_R

struct RevOrderElement {
  RevOrderElement * prev;
  int i;
};

Rcpp::List ClusterSLINK::as_hclust(const Rcpp::CharacterVector &seqnames) const {
  if (own_child && children.size() > 0) return delegate.as_hclust(seqnames);

  std::shared_lock<std::shared_timed_mutex> lock(this->mutex);
  Rcpp::IntegerMatrix merge(this->n - 1, 2);
  Rcpp::NumericVector height(this->n-1);
  Rcpp::IntegerVector order(this->n);

  // keep track of which cluster the active tips are in at each stage of clustering
  std::vector<int> clust_id;
  // ordering of earlier tips in the same cluster
  std::vector<RevOrderElement> ordering;
  // first element of each cluster (reverse_list does not know this)
  std::vector<RevOrderElement*> first;
  clust_id.reserve(this->n);
  ordering.reserve(this->n);
  first.reserve(this->n);

  std::vector<std::tuple<d_t, j_t, j_t>> slink;
  slink.reserve(this->n);
    int last_clust = 0;
  {
    std::shared_lock<std::shared_timed_mutex> lock(this->mutex);
    for (int i = 0; i < (int)n;) {
      slink.emplace_back(Lambda[i], -i, -Pi[i]);
      ++i;
      clust_id.push_back(-i);
      ordering.push_back({NULL, i});
      first.push_back(&ordering.back());
    }
    std::sort(slink.begin(), slink.end());

    for (const auto &c : slink) {
      double d = this->dconv.inverse(std::get<0>(c));
      j_t i = -std::get<1>(c);
      j_t j = -std::get<2>(c);
      merge(last_clust, 0) = clust_id[i];
      merge(last_clust, 1) = clust_id[j];
      height[last_clust] = d;
      first[j]->prev = &ordering[i];
      first[j] = first[i];
      clust_id[j] = ++last_clust;
    }
  }
  int j = 0;
  for (j_t i = 0; i < this->n - 1; ++i) {
    if (1 - clust_id[i] == i) {
      merge(last_clust, 0) = clust_id[i];
      merge(last_clust, 1) = clust_id[this->n - 1];
      height[last_clust] = 1.0;
      clust_id[this->n - 1] = ++last_clust;

      RevOrderElement * e = &ordering[i];
      while (e != NULL) {
        order[j] = e->i;
        ++j;
        e = e->prev;
      }
    }
  }

  Rcpp::List out;
  out["merge"] = merge;
  out["height"] = height;
  out["order"] = order;
  out["labels"] = seqnames;
  out["method"] = "single";
  out.attr("class") = "hclust";
  return out;
}
#endif // OPTIMOTU_R

void ClusterSLINK::merge_into(DistanceConsumer &consumer) {
  std::shared_lock<std::shared_timed_mutex> lock(this->mutex);
  for (j_t j = 0; j < this->n; j++) {
    consumer(Pi[j], j, dconv.inverse(Lambda[j]));
  }
}

void ClusterSLINK::merge_into(ClusterAlgorithm &consumer) {
  std::shared_lock<std::shared_timed_mutex> lock(this->mutex);
  for (j_t j = 0; j < this->n; j++) {
    consumer(Pi[j], j, Lambda[j]);
  }
}

void ClusterSLINK::merge_into_parent() {
  if (parent && !own_child) this->merge_into(*parent);
  if (own_child && children.size() > 0) this->merge_into(delegate);
}

double ClusterSLINK::max_relevant(j_t seq1, j_t seq2, int thread) const {

  if (seq1 <= seq2) {
    OPTIMOTU_CERR << "ClusterSLINK requires seq2 < seq1.  seq2=" << seq2
                  << ", seq1=" << seq1 << std::endl;
    OPTIMOTU_STOP("ClusterSLINK input error.");
  }
  if (seq1 < this->slink_seq1) {
    OPTIMOTU_CERR << "ClusterSLINK requires sequences in order."
                  << " Current seq1=" << seq1
                  << ", slink_seq1=" << this->slink_seq1 << std::endl;
    OPTIMOTU_STOP("ClusterSLINK input error.");
  }
  if (seq1 > slink_seq1) {
    return dconv.inverse(this->m - 1);
  }
  return dconv.inverse(this->M[seq2] - 1);
}
