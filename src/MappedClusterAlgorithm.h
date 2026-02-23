#ifndef OPTIMOTU_MAPPEDCLUSTERALGORITHM_H_INCLUDED
#define OPTIMOTU_MAPPEDCLUSTERALGORITHM_H_INCLUDED

#include "ClusterAlgorithm.h"

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
    // stub method
    void merge_into(DistanceConsumer &consumer) override;
    // stub method
    void merge_into(ClusterAlgorithm &consumer) override;
    // stub method
    void merge_into_parent() override;
    // stub method
    ClusterAlgorithm * make_child() override;
    // stub method
    ClusterAlgorithm * make_child(PairGenerator * pg) override;
    // stub method
    double max_relevant(j_t seq1, j_t seq2, int thread = 0) const override;
    #ifdef OPTIMOTU_R
      // stub method
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
    // stub method
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
    // stub method
    void operator()(PairGenerator & pg, double dist, int thread = 0) override;
    void operator()(j_t seq1, j_t seq2, int i, int thread = 0) override;
    // stub method
    void operator()(PairGenerator & pg, int i, int thread = 0) override;
    // stub method
    void merge_into(DistanceConsumer &consumer) override;
    // stub method
    void merge_into(ClusterAlgorithm &consumer) override;
    // stub method
    void merge_into_parent() override;
    // stub method
    ClusterAlgorithm * make_child() override;
    // stub method
    ClusterAlgorithm * make_child(PairGenerator * pg) override;
    // stub method
    double max_relevant(j_t seq1, j_t seq2, int thread = 0) const override;
    #ifdef OPTIMOTU_R
      // stub method
      Rcpp::List as_hclust(const Rcpp::CharacterVector &seqnames) const override;
    #endif // OPTIMOTU_R
  };

  // The forward map is dense, so it is implemented as a vector.
  const std::vector<std::size_t> fwd_map;

  // The reverse map is sparse, so it is implemented as an unordered map.
  const std::unordered_map<std::size_t, std::size_t> rev_map;

  const std::unique_ptr<Surrogate> surrogate;

  // The inner algorithm is the algorithm that is being wrapped.
  SingleClusterAlgorithm & inner;

  // Pointer to the PairGenerator used to create this object, if any.
  PairGenerator * const pair_generator = nullptr;

  SingleClusterAlgorithm * make_inner_child(ClusterAlgorithm * parent, const j_t n) override;

  // Explicit constructor for use in creating nested mappings
  MappedClusterAlgorithm(
    SingleClusterAlgorithm * parent,
    const std::vector<std::size_t> &fwd_map,
    const std::unordered_map<std::size_t, std::size_t> &rev_map
  );

public:
  // Constructor for child objects
  // This is called by ClusterAlgorithm::make_child(), which is used by
  // MergeClusterWorker to create child algorithms for each thread.
  // The mapping is determined by the PairGenerator.
  MappedClusterAlgorithm(SingleClusterAlgorithm * parent, PairGenerator * pg);

  // Special constructor for child objects where the parent has a delegate
  // to manage merges from the children (e.g. ClusterSLINK).
  MappedClusterAlgorithm(
    SingleClusterAlgorithm * parent,
    SingleClusterAlgorithm * delegate,
    PairGenerator * pg
  );

  // Constructor for parent objects
  // This is used by MultipleClusterAlgorithm to create algorithms which
  // cover subsets of the original sequences.
  // In this case, there is no outer parent, so we will never need to
  // merge into the parent, and there is no need for a surrogate.
  MappedClusterAlgorithm(SingleClusterAlgorithm & inner, const std::vector<std::size_t> &fwd_map);

  void operator()(j_t seq1, j_t seq2, double dist, int thread = 0) override;
  void operator()(j_t seq1, j_t seq2, int i, int thread = 0) override;
  void operator()(PairGenerator & pg, double dist, int thread = 0) override;
  void operator()(PairGenerator & pg, int i, int thread = 0) override;
  void merge_into(DistanceConsumer &consumer) override;
  void merge_into(ClusterAlgorithm &consumer) override;
  void merge_into_parent() override;
  SingleClusterAlgorithm * make_child() override;
  MappedClusterAlgorithm * make_child(PairGenerator * pg) override;
  double max_relevant(j_t seq1, j_t seq2, int thread = 0) const override;
  double max_relevant(PairGenerator & pg, int thread = 0) const override;
  void write_to_matrix(internal_matrix_t &out) override;
  Rcpp::List as_hclust(const Rcpp::CharacterVector &seqnames) const override;
};

#endif // OPTIMOTU_MAPPEDCLUSTERALGORITHM_H_INCLUDED
