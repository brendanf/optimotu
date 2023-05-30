#ifndef OPTIMOTU_MULTIPLECLUSTERALGORITHM_H_INCLUDED
#define OPTIMOTU_MULTIPLECLUSTERALGORITHM_H_INCLUDED

#include "single_linkage.h"
#include "DistanceConverter.h"
#include "ClusterAlgorithm.h"
#include <cassert>
#include <algorithm>

#ifdef OPTIMOTU_R
#include <Rcpp.h>
#include "ClusterTree.h"
#endif //OPTIMOTU_R

template <class CLUST_T>
class MultipleClusterAlgorithm : public DistanceConsumer {
protected:
  const DistanceConverter &dconv;
  const d_t m;
  const std::vector<std::string> &names;
  const std::vector<std::vector<std::string>> &subset_names;
  std::vector<CLUST_T> subsets;
  // for each element, which subsets does it belong to? sorted
  std::vector<std::vector<j_t>> subset_key;
  // for each subset, map from universal index to index in the subset
  std::vector<std::unordered_map<j_t, j_t>> fwd_map;
  // temp, declare once and reuse
  std::vector<j_t> whichsets;
public:
  MultipleClusterAlgorithm(
    const DistanceConverter &dconv,
    const std::vector<std::string> &names,
    const std::vector<std::vector<std::string>> &subset_names
  ) : dconv(dconv), m(dconv.m), names(names), subset_names(subset_names), subsets(),
  subset_key(names.size()), fwd_map(subset_names.size()) {
    subsets.reserve(subset_names.size());
    whichsets.reserve(subset_names.size());
    std::unordered_map<std::string, j_t> namekey;
    for (j_t i = 0; i < names.size(); i++) {
      namekey.emplace(names[i], i);
    }
    for (j_t i = 0; i < subset_names.size(); ++i) {
      subsets.emplace_back(dconv, subset_names[i].size(), m);
      fwd_map[i].reserve(subset_names[i].size());
      for (j_t j = 0; j < subset_names[i].size(); ++j) {
        auto f = namekey.find(subset_names[i][j]);
        assert(f != namekey.end());
        subset_key[f->second].push_back(i);
        fwd_map[i].emplace(f->second, j);
      }
    }
  };

  void operator()(j_t seq1, j_t seq2, double dist) {
    d_t i = dconv.convert(dist);
    (*this)(seq1, seq2, i);
  };

  void operator()(j_t seq1, j_t seq2, int i) {
    std::set_intersection(
      subset_key[seq1].begin(),
      subset_key[seq1].end(),
      subset_key[seq2].begin(),
      subset_key[seq2].end(),
      whichsets.begin()
    );
    for (j_t j : whichsets) {
      subsets[j](fwd_map[j][seq1], fwd_map[j][seq2], i);
    }
    whichsets.clear();
  };

  double min_relevant(j_t seq1, j_t seq2) {
    double min = dconv.inverse(m - 1);
    std::set_intersection(
      subset_key[seq1].begin(),
      subset_key[seq1].end(),
      subset_key[seq2].begin(),
      subset_key[seq2].end(),
      whichsets.begin()
    );
    for (j_t j : whichsets) {
      double minj = subsets[j].min_relevant(seq1, seq2);
      min = std::min(min, minj);
    }
    whichsets.clear();
    return min;
  };

  // send consumer() pairwise distances to ensure it is up-to-date with this
  // clustering
  template<class CLUST_U>
  void merge_into(MultipleClusterAlgorithm<CLUST_U> &consumer) {
    for (size_t i = 0; i < this->subsets.size(); i++) {
      this->subsets[i].merge_into(consumer.subsets[i]);
    }
  };

#ifdef OPTIMOTU_R
  void write_to_matrix(std::vector<RcppParallel::RMatrix<int>> &matrix_list) {
    for (size_t i = 0; i < this->subsets.size(); i++) {
      this->subsets[i].template write_to_matrix(matrix_list[i]);
    }
  }
  Rcpp::List as_hclust() = delete;
#endif //OPTIMOTU_R
};

#ifdef OPTIMOTU_R
template<>
Rcpp::List MultipleClusterAlgorithm<ClusterTree>::as_hclust() {
  std::vector<Rcpp::List> out;
  for (size_t i = 0; i < this->subsets.size(); i++) {
    out.push_back(this->subsets[i].as_hclust(Rcpp::wrap(this->subset_names[i])));
  }
  return Rcpp::wrap(out);
}
#endif //OPTIMOTU_R

#endif //OPTIMOTU_MULTIPLECLUSTERALGORITHM_H_INCLUDED
