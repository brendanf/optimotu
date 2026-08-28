#ifndef OPTIMOTU_MULTIPLECLUSTERALGORITHM_H_INCLUDED
#define OPTIMOTU_MULTIPLECLUSTERALGORITHM_H_INCLUDED

#include "single_linkage.h"
#include "DistanceConverter.h"
#include "ClusterAlgorithm.h"
#include "ClusterAlgorithmFactory.h"
#include "ClusterTree.h"
#include "MemoryBudget.h"
#include <unordered_map>
#include <unordered_set>

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
  std::shared_ptr<MemoryBudgetTracker> memory_budget;

  // for each element, which subsets does it belong to? sorted
  std::vector<std::vector<j_t>> subset_key;

  // for each subset, map from universal index to index in the subset
  std::vector<std::unordered_map<j_t, j_t>> fwd_map;

  // clustering algorithms which handle subsets
  std::vector<SingleClusterAlgorithm*> subsets;

  // owned subsets are subsets that are owned by this MultipleClusterAlgorithm
  std::vector<std::unique_ptr<SingleClusterAlgorithm>> owned_subsets;
  std::vector<std::pair<SingleClusterAlgorithm*, ClusterAlgorithm*>> borrowed_subsets;
  std::vector<std::size_t> tracked_allocations;
  std::size_t tracked_base_allocation = 0;

  // When set, routing tables and names come from the parent (tile merge child).
  const MultipleClusterAlgorithm *tile_routing_parent = nullptr;

  const std::vector<std::vector<j_t>> &routing_subset_key() const;
  const std::vector<std::unordered_map<j_t, j_t>> &routing_fwd_map() const;
  const std::vector<std::vector<std::string>> &routing_subset_names() const;

  // Last-pair overlap cache is thread_local inside ensure_whichsets so
  // concurrent workers do not share a slot. Returns the subsets that
  // contain both sequences.
  const std::vector<j_t> & ensure_whichsets(j_t seq1, j_t seq2) const;
  void apply_pair(
    j_t seq1,
    j_t seq2,
    d_t i,
    double dist,
    int thread,
    bool filter_irrelevant
  );

  // Protected main constructor: initializer list only (no body).
  // Used by MultipleClusterAlgorithmImpl to set up members before filling subsets.
  MultipleClusterAlgorithm(
    const ClusterAlgorithmFactory & factory,
    const std::vector<std::string> &names,
    const std::vector<std::vector<std::string>> &subset_names,
    const int threads,
    std::shared_ptr<MemoryBudgetTracker> memory_budget = nullptr
  );

  MultipleClusterAlgorithm(MultipleClusterAlgorithm * parent);

  // Tile-local merge child: shares parent routing tables.
  MultipleClusterAlgorithm(MultipleClusterAlgorithm *parent, PairGenerator *pg);

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
    const int threads = 1,
    std::shared_ptr<MemoryBudgetTracker> memory_budget = nullptr
  );

public:
  template<int verbose>
  friend class MultipleClusterAlgorithmImpl;

  MultipleClusterAlgorithm() = delete;
  ~MultipleClusterAlgorithm() override;

  void operator()(j_t seq1, j_t seq2, double dist, int thread) override;

  void operator()(j_t seq1, j_t seq2, int i, int thread) override;

  void operator()(PairGenerator & pg, double dist, int thread = 0) override;

  // Distance-matrix workers send every pair; do not per-subset-filter.
  void operator()(DistanceElement d, int thread = 0) override;

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

  bool accepts_unordered_pairs() const override;

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

  MultipleClusterAlgorithmImpl(MultipleClusterAlgorithm *parent, PairGenerator *pg);

  MultipleClusterAlgorithmImpl(
    MultipleClusterAlgorithm * parent,
    const std::vector<std::string> names,
    const std::vector<std::vector<j_t>> subset_indices,
    const std::vector<std::vector<std::string>> subset_names,
    const std::vector<std::vector<j_t>> subset_key,
    const std::vector<std::unordered_map<j_t, j_t>> fwd_map,
    const std::vector<j_t> child_to_parent_map,
    PairGenerator * pg,
    const int threads = 1,
    std::shared_ptr<MemoryBudgetTracker> memory_budget = nullptr
  );

public:
  MultipleClusterAlgorithmImpl(
    const ClusterAlgorithmFactory & factory,
    const std::vector<std::string> &names,
    const std::vector<std::vector<std::string>> &subset_names,
    const int threads,
    int verbose_param,
    std::shared_ptr<MemoryBudgetTracker> memory_budget = nullptr
  );

  MultipleClusterAlgorithm * make_child() override;
  MultipleClusterAlgorithm * make_child(PairGenerator * pg) override;
};

#endif //OPTIMOTU_MULTIPLECLUSTERALGORITHM_H_INCLUDED
