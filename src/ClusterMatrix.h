#ifndef _ClusterMatrix_
#define _ClusterMatrix_

#include "ClusterAlgorithm.h"
#ifdef OPTIMOTU_R
#include <RcppParallel.h>
#endif

#define LINEAR_FILL 1
#define BINARY_FILL 2
#define TOPDOWN_FILL 3

// Matrix representation of a hierarchical clustering of n items at m different
// thresholds.
template <class ARRAY_T = std::vector<int>,
          bool BINARY_SEARCH=true, //alt would be linear
          int FILL=2>
struct ClusterMatrix : public ClusterAlgorithm {
  template<typename, bool, int> friend class ClusterMatrix;
private:
  ARRAY_T clust_array;
  int *const ca;
  std::vector<int> toclust;

protected:
  void initialize();

  ClusterMatrix(ClusterAlgorithm * parent) :
  ClusterAlgorithm(parent), clust_array(m*n), ca(&clust_array[0]), toclust(m, 0) {
    initialize();
  };

public:
  ClusterMatrix(const DistanceConverter &dconv, size_t n, int m) :
  ClusterAlgorithm(dconv, n, m), clust_array(m*n), ca(&clust_array[0]), toclust(m, 0) {
    initialize();
  };

#ifdef OPTIMOTU_R
  ClusterMatrix(const DistanceConverter &dconv, Rcpp::IntegerMatrix &parent) = delete;
#endif

  using ClusterAlgorithm::operator();

  void operator()(j_t seq1, j_t seq2, d_t i) override;

  void merge_into(DistanceConsumer &consumer) override final;

  void merge_into(ClusterAlgorithm &consumer) override final;

  ClusterAlgorithm * make_child() override;

  double max_relevant(j_t seq1, j_t seq2) const override;

#ifdef OPTIMOTU_R
  void write_to_matrix(RcppParallel::RMatrix<int> &out) override {
    tbb::queuing_rw_mutex::scoped_lock lock(this->mutex, true);
    std::copy(clust_array.begin(), clust_array.end(), out.begin());
  }
#endif
};

#ifdef OPTIMOTU_R
template<>
ClusterMatrix<RcppParallel::RMatrix<int>, true, LINEAR_FILL>::ClusterMatrix(
  const DistanceConverter &dconv, Rcpp::IntegerMatrix &im
);
template<>
ClusterMatrix<RcppParallel::RMatrix<int>, true, LINEAR_FILL>::ClusterMatrix(
    ClusterAlgorithm * parent) = delete;
template<>
ClusterMatrix<RcppParallel::RMatrix<int>, true, LINEAR_FILL>::ClusterMatrix(
    const DistanceConverter &dconv, size_t n, int m) = delete;

template<>
ClusterMatrix<RcppParallel::RMatrix<int>, true, BINARY_FILL>::ClusterMatrix(
    const DistanceConverter &dconv, Rcpp::IntegerMatrix &im
);
template<>
ClusterMatrix<RcppParallel::RMatrix<int>, true, BINARY_FILL>::ClusterMatrix(
    ClusterAlgorithm * parent) = delete;
template<>
ClusterMatrix<RcppParallel::RMatrix<int>, true, BINARY_FILL>::ClusterMatrix(
    const DistanceConverter &dconv, size_t n, int m) = delete;

template<>
ClusterMatrix<RcppParallel::RMatrix<int>, true, TOPDOWN_FILL>::ClusterMatrix(
    const DistanceConverter &dconv, Rcpp::IntegerMatrix &im
);
template<>
ClusterMatrix<RcppParallel::RMatrix<int>, true, TOPDOWN_FILL>::ClusterMatrix(
    ClusterAlgorithm * parent) = delete;
template<>
ClusterMatrix<RcppParallel::RMatrix<int>, true, TOPDOWN_FILL>::ClusterMatrix(
    const DistanceConverter &dconv, size_t n, int m) = delete;

template<>
ClusterMatrix<RcppParallel::RMatrix<int>, false, LINEAR_FILL>::ClusterMatrix(
    const DistanceConverter &dconv, Rcpp::IntegerMatrix &im
);
template<>
ClusterMatrix<RcppParallel::RMatrix<int>, false, LINEAR_FILL>::ClusterMatrix(
    ClusterAlgorithm * parent) = delete;
template<>
ClusterMatrix<RcppParallel::RMatrix<int>, false, LINEAR_FILL>::ClusterMatrix(
    const DistanceConverter &dconv, size_t n, int m) = delete;

template<>
ClusterMatrix<RcppParallel::RMatrix<int>, false, BINARY_FILL>::ClusterMatrix(
    const DistanceConverter &dconv, Rcpp::IntegerMatrix &im
);
template<>
ClusterMatrix<RcppParallel::RMatrix<int>, false, BINARY_FILL>::ClusterMatrix(
    ClusterAlgorithm * parent) = delete;
template<>
ClusterMatrix<RcppParallel::RMatrix<int>, false, BINARY_FILL>::ClusterMatrix(
    const DistanceConverter &dconv, size_t n, int m) = delete;

template<>
ClusterMatrix<RcppParallel::RMatrix<int>, false, TOPDOWN_FILL>::ClusterMatrix(
    const DistanceConverter &dconv, Rcpp::IntegerMatrix &im
);
template<>
ClusterMatrix<RcppParallel::RMatrix<int>, false, TOPDOWN_FILL>::ClusterMatrix(
    ClusterAlgorithm * parent) = delete;
template<>
ClusterMatrix<RcppParallel::RMatrix<int>, false, TOPDOWN_FILL>::ClusterMatrix(
    const DistanceConverter &dconv, size_t n, int m) = delete;

template class ClusterMatrix<RcppParallel::RMatrix<int>, true, LINEAR_FILL>;
template class ClusterMatrix<RcppParallel::RMatrix<int>, true, BINARY_FILL>;
template class ClusterMatrix<RcppParallel::RMatrix<int>, true, TOPDOWN_FILL>;
template class ClusterMatrix<RcppParallel::RMatrix<int>, false, LINEAR_FILL>;
template class ClusterMatrix<RcppParallel::RMatrix<int>, false, BINARY_FILL>;
template class ClusterMatrix<RcppParallel::RMatrix<int>, false, TOPDOWN_FILL>;
#endif

template class ClusterMatrix<std::vector<int>, true, LINEAR_FILL>;
template class ClusterMatrix<std::vector<int>, true, BINARY_FILL>;
template class ClusterMatrix<std::vector<int>, true, TOPDOWN_FILL>;
template class ClusterMatrix<std::vector<int>, false, LINEAR_FILL>;
template class ClusterMatrix<std::vector<int>, false, BINARY_FILL>;
template class ClusterMatrix<std::vector<int>, false, TOPDOWN_FILL>;

#endif
