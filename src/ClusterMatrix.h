#ifndef _ClusterMatrix_
#define _ClusterMatrix_

#include "ClusterAlgorithm.h"
#ifdef OPTIMOTU_R
#include <RcppParallel.h>
#endif

// #define BINARY_MAX  (uint8_t)0b00000001
// #define BINARY_FILL (uint8_t)0b00000010
// #define TOPDOWN_FILL (uint8_t)0b00000100

// Matrix representation of a hierarchical clustering of n items at m different
// thresholds.
template <class ARRAY_T = std::vector<int>,
          bool BINARY_MAX=true,
          bool BINARY_FILL=true,
          bool TOPDOWN_FILL=false>
struct ClusterMatrix : public ClusterAlgorithm {
  static_assert(!(BINARY_FILL && TOPDOWN_FILL),
                "cannot use both BINARY_FILL and TOPDOWN_FILL");
  template<typename, bool, bool, bool> friend class ClusterMatrix;
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

  void merge_into(DistanceConsumer &consumer) const override final;

  void merge_into(ClusterAlgorithm &consumer) const override final;

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
ClusterMatrix<RcppParallel::RMatrix<int>, true, true, false>::ClusterMatrix(
  const DistanceConverter &dconv, Rcpp::IntegerMatrix &im
);
template<>
ClusterMatrix<RcppParallel::RMatrix<int>, true, true, false>::ClusterMatrix(
    ClusterAlgorithm * parent) = delete;
template<>
ClusterMatrix<RcppParallel::RMatrix<int>, true, true, false>::ClusterMatrix(
    const DistanceConverter &dconv, size_t n, int m) = delete;

template<>
ClusterMatrix<RcppParallel::RMatrix<int>, true, false, true>::ClusterMatrix(
    const DistanceConverter &dconv, Rcpp::IntegerMatrix &im
);
template<>
ClusterMatrix<RcppParallel::RMatrix<int>, true, false, true>::ClusterMatrix(
    ClusterAlgorithm * parent) = delete;
template<>
ClusterMatrix<RcppParallel::RMatrix<int>, true, false, true>::ClusterMatrix(
    const DistanceConverter &dconv, size_t n, int m) = delete;

template<>
ClusterMatrix<RcppParallel::RMatrix<int>, true, false, false>::ClusterMatrix(
    const DistanceConverter &dconv, Rcpp::IntegerMatrix &im
);
template<>
ClusterMatrix<RcppParallel::RMatrix<int>, true, false, false>::ClusterMatrix(
    ClusterAlgorithm * parent) = delete;
template<>
ClusterMatrix<RcppParallel::RMatrix<int>, true, false, false>::ClusterMatrix(
    const DistanceConverter &dconv, size_t n, int m) = delete;

template<>
ClusterMatrix<RcppParallel::RMatrix<int>, false, true, false>::ClusterMatrix(
    const DistanceConverter &dconv, Rcpp::IntegerMatrix &im
);
template<>
ClusterMatrix<RcppParallel::RMatrix<int>, false, true, false>::ClusterMatrix(
    ClusterAlgorithm * parent) = delete;
template<>
ClusterMatrix<RcppParallel::RMatrix<int>, false, true, false>::ClusterMatrix(
    const DistanceConverter &dconv, size_t n, int m) = delete;

template<>
ClusterMatrix<RcppParallel::RMatrix<int>, false, false, true>::ClusterMatrix(
    const DistanceConverter &dconv, Rcpp::IntegerMatrix &im
);
template<>
ClusterMatrix<RcppParallel::RMatrix<int>, false, false, true>::ClusterMatrix(
    ClusterAlgorithm * parent) = delete;
template<>
ClusterMatrix<RcppParallel::RMatrix<int>, false, false, true>::ClusterMatrix(
    const DistanceConverter &dconv, size_t n, int m) = delete;

template<>
ClusterMatrix<RcppParallel::RMatrix<int>, false, false, false>::ClusterMatrix(
    const DistanceConverter &dconv, Rcpp::IntegerMatrix &im
);
template<>
ClusterMatrix<RcppParallel::RMatrix<int>, false, false, false>::ClusterMatrix(
    ClusterAlgorithm * parent) = delete;
template<>
ClusterMatrix<RcppParallel::RMatrix<int>, false, false, false>::ClusterMatrix(
    const DistanceConverter &dconv, size_t n, int m) = delete;

template class ClusterMatrix<RcppParallel::RMatrix<int>, true, true, false>;
template class ClusterMatrix<RcppParallel::RMatrix<int>, true, false, true>;
template class ClusterMatrix<RcppParallel::RMatrix<int>, true, false, false>;
template class ClusterMatrix<RcppParallel::RMatrix<int>, false, true, false>;
template class ClusterMatrix<RcppParallel::RMatrix<int>, false, false, true>;
template class ClusterMatrix<RcppParallel::RMatrix<int>, false, false, false>;
#endif

template class ClusterMatrix<std::vector<int>, true, true, false>;
template class ClusterMatrix<std::vector<int>, true, false, true>;
template class ClusterMatrix<std::vector<int>, true, false, false>;
template class ClusterMatrix<std::vector<int>, false, true, false>;
template class ClusterMatrix<std::vector<int>, false, false, true>;
template class ClusterMatrix<std::vector<int>, false, false, false>;

#endif
