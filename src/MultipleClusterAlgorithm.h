#ifndef OPTIMOTU_MULTIPLECLUSTERALGORITHM_H_INCLUDED
#define OPTIMOTU_MULTIPLECLUSTERALGORITHM_H_INCLUDED

#include "single_linkage.h"
#include "DistanceConverter.h"
#include "ClusterAlgorithm.h"
#include "ClusterAlgorithmFactory.h"
#include "ClusterTree.h"

class MultipleClusterAlgorithm : public ClusterAlgorithm {
protected:
  const ClusterAlgorithmFactory & factory;
  const std::vector<std::string> &names;
  const std::vector<std::vector<std::string>> &subset_names;
  const int threads;
  std::vector<std::unique_ptr<SingleClusterAlgorithm>> subsets;
  // for each element, which subsets does it belong to? sorted
  std::vector<std::vector<j_t>> subset_key;
  // for each subset, map from universal index to index in the subset
  std::vector<std::unordered_map<j_t, j_t>> fwd_map;
  // temp, declare once (per thread) and reuse
  mutable std::vector<std::vector<j_t>> whichsets;
  // remember which the last sequence pairs were
  mutable std::vector<std::pair<j_t, j_t>> ws_keys;

  MultipleClusterAlgorithm(MultipleClusterAlgorithm * parent);

public:
  MultipleClusterAlgorithm(
    const ClusterAlgorithmFactory & factory,
    const std::vector<std::string> &names,
    const std::vector<std::vector<std::string>> &subset_names,
    const int threads = 1
  );

  void operator()(j_t seq1, j_t seq2, double dist, int thread) override;

  void operator()(j_t seq1, j_t seq2, int i, int thread) override;

  double max_relevant(j_t seq1, j_t seq2, int thread) const override;

  // send consumer() pairwise distances to ensure it is up-to-date with this
  // clustering
  void merge_into(DistanceConsumer &consumer) override;

  // send consumer() pairwise distances to ensure it is up-to-date with this
  // clustering
  void merge_into(ClusterAlgorithm &consumer) override;

  // send consumer() pairwise distances to ensure it is up-to-date with this
  // clustering
  virtual void merge_into(MultipleClusterAlgorithm &consumer);

  void merge_into_parent() override;

  MultipleClusterAlgorithm * make_child() override;

  void write_to_matrix(std::vector<internal_matrix_t> &matrix_list);

#ifdef OPTIMOTU_R
  Rcpp::List as_hclust(const Rcpp::CharacterVector &seqnames) const override;

  Rcpp::List as_hclust() const;
#endif //OPTIMOTU_R
};

#endif //OPTIMOTU_MULTIPLECLUSTERALGORITHM_H_INCLUDED
