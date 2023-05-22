#ifndef _MULTIPLE_CLUSTER_ALGORITHM_
#define _MULTIPLE_CLUSTER_ALGORITHM_
#include "single_linkage.h"
#include "DistanceConverter.h"
#include "ClusterAlgorithm.h"

template <class CLUST_T>
class MultipleClusterAlgorithm : DistanceConsumer {
private:
  const DistanceConverter &dconv;
  const d_t m;
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
    const std::vector<std::vector<std::string>> &subset_names,
    const d_t m
  ) : dconv(dconv), m(m), subsets(), subset_key(names.size()), fwd_map(subset_names.size()) {
    subsets.reserve(subset_names.size());
    whichsets.reserve(subset_names.size());
    std::unordered_map<std::string, j_t> namekey;
    for (j_t i = 0; i < names.size(); i++) {
      namekey.emplace(names[i], i);
    }
    for (int i = 0; i < subset_names.size(); ++i) {
      subsets.emplace_back(dconv, subset_names[i].size(), m);
      fwd_map[i].reserve(subset_names[i].size());
      for (int j = 0; j < subset_names[i].size(); ++j) {
        auto f = namekey.find(subset_names[i][j]);
        if (f == namekey.end()) Rcpp::stop("mismatched name");
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

};
#endif
