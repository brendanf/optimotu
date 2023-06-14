#ifndef OPTIMOTU_CLUSTERMATRIX_H_INCLUDED
#define OPTIMOTU_CLUSTERMATRIX_H_INCLUDED

#include "ClusterAlgorithm.h"
#ifdef OPTIMOTU_R
#include <RcppParallel.h>
#endif

#define LINEAR_FILL 1
#define BINARY_FILL 2
#define TOPDOWN_FILL 3

// Matrix representation of a hierarchical clustering of n items at m different
// thresholds.
template <bool BINARY_SEARCH=true, //alt would be linear
          int FILL=2,
          class ARRAY_T = std::vector<int>>
class ClusterMatrix : public SingleClusterAlgorithm {
  template<bool, int, typename> friend class ClusterMatrix;
private:
  ARRAY_T clust_array;
  int *const ca;
  std::vector<int> toclust;

protected:
  void initialize();

  ClusterMatrix(SingleClusterAlgorithm * parent) :
    SingleClusterAlgorithm(parent),
    clust_array(m*n), ca(&clust_array[0]), toclust(m, 0) {
    initialize();
  };

public:
  ClusterMatrix(const DistanceConverter &dconv, size_t n) :
  SingleClusterAlgorithm(dconv, n),
  clust_array(m*n), ca(&clust_array[0]), toclust(m, 0) {
    initialize();
  };

  ClusterMatrix(const DistanceConverter &dconv, init_matrix_t im) :
    SingleClusterAlgorithm(dconv, im),
    clust_array(im), ca(&clust_array[0]), toclust(m, 0) {
    initialize();
  };

  using ClusterAlgorithm::operator();

  void operator()(j_t seq1, j_t seq2, d_t i, int thread = 0) override;

  void merge_into(DistanceConsumer &consumer) override final;

  void merge_into(ClusterAlgorithm &consumer) override final;

  SingleClusterAlgorithm * make_child() override;

  double max_relevant(j_t seq1, j_t seq2, int thread = 0) const override;

  void write_to_matrix(internal_matrix_t &out) override;

#ifdef OPTIMOTU_R
  Rcpp::List as_hclust(const Rcpp::CharacterVector &seqnames) const override;
#endif
};

template<>
ClusterMatrix<true, LINEAR_FILL>::ClusterMatrix(
  const DistanceConverter &dconv, init_matrix_t
) = delete;

template<>
ClusterMatrix<true, BINARY_FILL>::ClusterMatrix(
    const DistanceConverter &dconv, init_matrix_t
) = delete;

template<>
ClusterMatrix<true, TOPDOWN_FILL>::ClusterMatrix(
    const DistanceConverter &dconv, init_matrix_t
) = delete;

template<>
ClusterMatrix<false, LINEAR_FILL>::ClusterMatrix(
    const DistanceConverter &dconv, init_matrix_t
) = delete;

template<>
ClusterMatrix<false, BINARY_FILL>::ClusterMatrix(
    const DistanceConverter &dconv, init_matrix_t
) = delete;

template<>
ClusterMatrix<false, TOPDOWN_FILL>::ClusterMatrix(
    const DistanceConverter &dconv, init_matrix_t
) = delete;

#define deleted_funcs template<>                              \
cm::ClusterMatrix(SingleClusterAlgorithm * parent) = delete;         \
template<>                                                     \
cm::ClusterMatrix(const DistanceConverter &dconv, size_t n) = delete

#define cm ClusterMatrix<true, LINEAR_FILL, internal_matrix_t>
deleted_funcs;
#define cm ClusterMatrix<true, BINARY_FILL, internal_matrix_t>
deleted_funcs;
#define cm ClusterMatrix<true, TOPDOWN_FILL, internal_matrix_t>
deleted_funcs;
#define cm ClusterMatrix<false, LINEAR_FILL, internal_matrix_t>
deleted_funcs;
#define cm ClusterMatrix<false, BINARY_FILL, internal_matrix_t>
deleted_funcs;
#define cm ClusterMatrix<false, TOPDOWN_FILL, internal_matrix_t>
deleted_funcs;
#undef cm
#undef deleted_funcs

template class ClusterMatrix<true, LINEAR_FILL, internal_matrix_t>;
template class ClusterMatrix<true, BINARY_FILL, internal_matrix_t>;
template class ClusterMatrix<true, TOPDOWN_FILL, internal_matrix_t>;
template class ClusterMatrix<false, LINEAR_FILL, internal_matrix_t>;
template class ClusterMatrix<false, BINARY_FILL, internal_matrix_t>;
template class ClusterMatrix<false, TOPDOWN_FILL, internal_matrix_t>;

template class ClusterMatrix<true, LINEAR_FILL>;
template class ClusterMatrix<true, BINARY_FILL>;
template class ClusterMatrix<true, TOPDOWN_FILL>;
template class ClusterMatrix<false, LINEAR_FILL>;
template class ClusterMatrix<false, BINARY_FILL>;
template class ClusterMatrix<false, TOPDOWN_FILL>;

#endif //OPTIMOTU_CLUSTERMATRIX_H_INCLUDED
