#ifndef OPTIMOTU_MAPPEDCLUSTERALGORITHM_H_INCLUDED
#define OPTIMOTU_MAPPEDCLUSTERALGORITHM_H_INCLUDED

#include "ClusterAlgorithm.h"
#include <memory>

std::unordered_map<std::size_t, std::size_t> invert_map(
    const std::vector<std::size_t> &fwd_map);

// child_local -> parent_subset_local for sequences in subset intersect tile.
std::vector<std::size_t> subset_tile_fwd_map(
    const std::unordered_map<j_t, j_t> &global_to_subset,
    const PairGenerator &pg);

// A MappedClusterAlgorithm is a ClusterAlgorithm that wraps another
// ClusterAlgorithm and maps the indices of the inner algorithm to external
// indices. This is typically used when the inner algorithm only accepts a
// subset of the external indices.
// The forward map is the mapping from internal indices to external indices.
// It is used by functions which output results, such as `merge_into()`,
// `write_to_matrix()`, and `as_hclust()`.
// The reverse map is the mapping from external indices to internal indices.
// It is used by functions which do clustering calculations, such as
// `operator()()` and `max_relevant()`.
class MappedClusterAlgorithm : public SingleClusterAlgorithm {
protected:

  virtual SingleClusterAlgorithm & get_inner() const = 0;

  // Constructors for use by Impl: only initialize SingleClusterAlgorithm.
  MappedClusterAlgorithm(ClusterAlgorithm * parent, j_t n);
  MappedClusterAlgorithm(const DistanceConverter & dconv, j_t n);

  virtual SingleClusterAlgorithm * make_inner_child(ClusterAlgorithm * parent, const j_t n) override = 0;

public:

  void operator()(j_t seq1, j_t seq2, double dist, int thread = 0) override = 0;
  void operator()(j_t seq1, j_t seq2, int i, int thread = 0) override = 0;
  void operator()(PairGenerator & pg, double dist, int thread = 0) override = 0;
  void operator()(PairGenerator & pg, int i, int thread = 0) override = 0;
  void merge_into(DistanceConsumer &consumer) override = 0;
  void merge_into(ClusterAlgorithm &consumer) override = 0;
  void merge_into_parent() override;
  void finalize() override;
  SingleClusterAlgorithm * make_child() override = 0;
  MappedClusterAlgorithm * make_child(PairGenerator * pg) override = 0;
  double max_relevant(j_t seq1, j_t seq2, int thread = 0) const override = 0;
  double max_relevant(PairGenerator & pg, int thread = 0) const override = 0;
  bool accepts_unordered_pairs() const override;
  void write_to_matrix(internal_matrix_t &out) override;
  void prepare_output() override;
  void write_threshold_row(d_t t, int *dest) const override;
  Rcpp::List as_hclust(const Rcpp::CharacterVector &seqnames) const override;
};

template<int verbose = 0>
class MappedClusterAlgorithmImpl : public MappedClusterAlgorithm {
protected:
  // The forward map is dense, so it is implemented as a vector.
  const std::vector<std::size_t> fwd_map;

  // The reverse map is sparse, so it is implemented as an unordered map.
  const std::unordered_map<std::size_t, std::size_t> rev_map;

  // Surrogate used as the parent of the inner algorithm (child constructors only).
  std::unique_ptr<ClusterAlgorithm> surrogate;

  // The inner algorithm is the algorithm that is being wrapped.
  SingleClusterAlgorithm * const inner;

  // Pointer to the PairGenerator used to create this object, if any.
  PairGenerator * const pair_generator = nullptr;

  SingleClusterAlgorithm & get_inner() const override { return *inner; }

  // The surrogate is used as the parent of the inner algorithm.
  // It is used as an interface between the inner algorithm and the parent
  // of the MappedClusterAlgorithm.
  class Surrogate: public ClusterAlgorithm {
    using ClusterAlgorithm::parent;
    const std::vector<std::size_t> &fwd_map;
  public:
    Surrogate(ClusterAlgorithm * parent, const std::vector<std::size_t> &fwd_map);

    void operator()(j_t seq1, j_t seq2, double dist, int thread = 0) override;
    void operator()(j_t seq1, j_t seq2, int i, int thread = 0) override;
    void merge_into(DistanceConsumer &consumer) override;
    void merge_into(ClusterAlgorithm &consumer) override;
    void merge_into_parent() override;
    ClusterAlgorithm * make_child() override;
    ClusterAlgorithm * make_child(PairGenerator * pg) override;
    double max_relevant(j_t seq1, j_t seq2, int thread = 0) const override;
    bool accepts_unordered_pairs() const override;
    #ifdef OPTIMOTU_R
      Rcpp::List as_hclust(const Rcpp::CharacterVector &seqnames) const override;
    #endif // OPTIMOTU_R
  };

  // DistanceForwarder is used to forward operator()() calls from the inner
  // algorithm as a result of merge_into(DistanceConsumer).
  class DistanceForwarder: public DistanceConsumer {
    const std::vector<std::size_t> &fwd_map;
    DistanceConsumer * const target;
  public:
    DistanceForwarder(const std::vector<std::size_t> &fwd_map, DistanceConsumer * target);
    void operator()(j_t seq1, j_t seq2, double dist, int thread = 0) override;
    void operator()(PairGenerator & pg, double dist, int thread = 0) override;
  };

  // IndexForwarder is used to forward operator()() calls from the inner
  // algorithm as a result of merge_into(ClusterAlgorithm).
  class IndexForwarder: public ClusterAlgorithm {
    const std::vector<std::size_t> &fwd_map;
    ClusterAlgorithm * const target;

  public:
    IndexForwarder(const std::vector<std::size_t> &fwd_map, ClusterAlgorithm * target);
    void operator()(j_t seq1, j_t seq2, double dist, int thread = 0) override;
    void operator()(PairGenerator & pg, double dist, int thread = 0) override;
    void operator()(j_t seq1, j_t seq2, int i, int thread = 0) override;
    void operator()(PairGenerator & pg, int i, int thread = 0) override;
    void merge_into(DistanceConsumer &consumer) override;
    void merge_into(ClusterAlgorithm &consumer) override;
    void merge_into_parent() override;
    ClusterAlgorithm * make_child() override;
    ClusterAlgorithm * make_child(PairGenerator * pg) override;
    double max_relevant(j_t seq1, j_t seq2, int thread = 0) const override;
    bool accepts_unordered_pairs() const override;
    #ifdef OPTIMOTU_R
      Rcpp::List as_hclust(const Rcpp::CharacterVector &seqnames) const override;
    #endif // OPTIMOTU_R
  };

  SingleClusterAlgorithm * make_inner_child(ClusterAlgorithm * parent, const j_t n) override;

public:
  // Constructor for child objects (mapping determined by PairGenerator)
  MappedClusterAlgorithmImpl(SingleClusterAlgorithm * parent, PairGenerator * pg);

  // Special constructor for child objects where the parent has a delegate
  MappedClusterAlgorithmImpl(
    SingleClusterAlgorithm * parent,
    SingleClusterAlgorithm * delegate,
    PairGenerator * pg
  );

  // Child with explicit mapping; inner merges into delegate (SLINK tiles).
  MappedClusterAlgorithmImpl(
      SingleClusterAlgorithm *parent,
      SingleClusterAlgorithm *delegate,
      const std::vector<std::size_t> &fwd_map,
      const std::unordered_map<std::size_t, std::size_t> &rev_map);

  // Constructor for parent objects (no surrogate)
  MappedClusterAlgorithmImpl(SingleClusterAlgorithm & inner_ref, const std::vector<std::size_t> &fwd_map);

  // Constructor for nested mappings
  MappedClusterAlgorithmImpl(
    SingleClusterAlgorithm * parent,
    const std::vector<std::size_t> &fwd_map,
    const std::unordered_map<std::size_t, std::size_t> &rev_map
  );

  void operator()(j_t seq1, j_t seq2, double dist, int thread = 0) override;
  void operator()(j_t seq1, j_t seq2, int i, int thread = 0) override;
  void operator()(PairGenerator & pg, double dist, int thread = 0) override;
  void operator()(PairGenerator & pg, int i, int thread = 0) override;
  void merge_into(DistanceConsumer &consumer) override;
  void merge_into(ClusterAlgorithm &consumer) override;
  SingleClusterAlgorithm * make_child() override;
  MappedClusterAlgorithm * make_child(PairGenerator * pg) override;
  double max_relevant(j_t seq1, j_t seq2, int thread = 0) const override;
  double max_relevant(PairGenerator & pg, int thread = 0) const override;
};

#endif // OPTIMOTU_MAPPEDCLUSTERALGORITHM_H_INCLUDED
