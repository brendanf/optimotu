#ifndef OPTIMOTU_MULTIPLECLUSTERALGORITHM_H_INCLUDED
#define OPTIMOTU_MULTIPLECLUSTERALGORITHM_H_INCLUDED

#include "single_linkage.h"
#include "DistanceConverter.h"
#include "ClusterAlgorithm.h"
#include "ClusterAlgorithmFactory.h"
#include "ClusterTree.h"
#include <unordered_map>

class MultipleClusterAlgorithm : public ClusterAlgorithm {
protected:
  const ClusterAlgorithmFactory & factory;
  // the names of the sequences
  const std::vector<std::string> names;
  // the indices of the sequences in each subset
  const std::vector<std::vector<j_t>> subset_indices;
  // the names of the sequences in each subset
  const std::vector<std::vector<std::string>> subset_names;
  // the number of threads
  const int threads;

  // for each element, which subsets does it belong to? sorted
  std::vector<std::vector<j_t>> subset_key;

  // for each subset, map from universal index to index in the subset
  std::vector<std::unordered_map<j_t, j_t>> fwd_map;

  // clustering algorithms which handle subsets
  std::vector<SingleClusterAlgorithm*> subsets;

  // owned subsets are subsets that are owned by this MultipleClusterAlgorithm
  std::vector<std::unique_ptr<SingleClusterAlgorithm>> owned_subsets;

  // temp, declare once (per thread) and reuse
  mutable std::vector<std::vector<j_t>> whichsets;

  // remember which the last sequence pairs were
  mutable std::vector<std::pair<j_t, j_t>> ws_keys;

  // Protected main constructor: initializer list only (no body).
  // Used by MultipleClusterAlgorithmImpl to set up members before filling subsets.
  MultipleClusterAlgorithm(
    const ClusterAlgorithmFactory & factory,
    const std::vector<std::string> &names,
    const std::vector<std::vector<std::string>> &subset_names,
    const int threads
  );

  MultipleClusterAlgorithm(MultipleClusterAlgorithm * parent);

  // Constructor for child objects
  // It makes sense to calculate many const members all together,
  // rather than serially, so these are calculated in make_child().
  MultipleClusterAlgorithm(
    MultipleClusterAlgorithm * parent,
    const std::vector<std::string> names,
    const std::vector<std::vector<j_t>> subset_indices,
    const std::vector<std::vector<std::string>> subset_names,
    const std::vector<std::vector<j_t>> subset_key,
    const std::vector<std::unordered_map<j_t, j_t>> fwd_map,
    const std::vector<j_t> child_to_parent_map,
    PairGenerator * pg,
    const int threads = 1
  );

public:
  template<int verbose>
  friend class MultipleClusterAlgorithmImpl;

  MultipleClusterAlgorithm() = delete;

  void operator()(j_t seq1, j_t seq2, double dist, int thread) override;

  void operator()(j_t seq1, j_t seq2, int i, int thread) override;

  virtual void finalize() override;

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

  MultipleClusterAlgorithm * make_child(PairGenerator * pg) override;

  void write_to_matrix(std::vector<internal_matrix_t> &matrix_list);

#ifdef OPTIMOTU_R
  Rcpp::List as_hclust(const Rcpp::CharacterVector &seqnames) const override;

  Rcpp::List as_hclust() const;
#endif //OPTIMOTU_R
};

template<int verbose>
class MultipleClusterAlgorithmImpl : public MultipleClusterAlgorithm {
protected:
  MultipleClusterAlgorithmImpl(MultipleClusterAlgorithm * parent);

  MultipleClusterAlgorithmImpl(
    MultipleClusterAlgorithm * parent,
    const std::vector<std::string> names,
    const std::vector<std::vector<j_t>> subset_indices,
    const std::vector<std::vector<std::string>> subset_names,
    const std::vector<std::vector<j_t>> subset_key,
    const std::vector<std::unordered_map<j_t, j_t>> fwd_map,
    const std::vector<j_t> child_to_parent_map,
    PairGenerator * pg,
    const int threads = 1
  );

public:
  MultipleClusterAlgorithmImpl(
    const ClusterAlgorithmFactory & factory,
    const std::vector<std::string> &names,
    const std::vector<std::vector<std::string>> &subset_names,
    const int threads,
    int verbose_param
  );

  MultipleClusterAlgorithm * make_child() override;
  MultipleClusterAlgorithm * make_child(PairGenerator * pg) override;
};

#endif //OPTIMOTU_MULTIPLECLUSTERALGORITHM_H_INCLUDED
