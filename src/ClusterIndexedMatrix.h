#ifndef OPTIMOTU_CLUSTERINDEXEDMATRIX_H_INCLUDED
#define OPTIMOTU_CLUSTERINDEXEDMATRIX_H_INCLUDED

#include "ClusterAlgorithm.h"
#ifdef OPTIMOTU_R
#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#endif

// Matrix representation of a hierarchical clustering of n items at m different
// thresholds.
template <class ARRAY_T = std::vector<int>>
struct ClusterIndexedMatrix : public ClusterAlgorithm {
  template<typename> friend class ClusterIndexedMatrix;

private:
  struct tip {
    j_t j = NO_CLUST;
    tip *prev = nullptr, *next = nullptr;
    int *column = nullptr;
    d_t prev_d = NO_DIST, next_d = NO_DIST;
  };

  ARRAY_T clust_array;
  int *const ca;
  int *buffer;
  tip *index;
  tip tfwd, trev;

protected:

  void initialize();

  ClusterIndexedMatrix(ClusterAlgorithm * parent);

public:
  ClusterIndexedMatrix(
    const DistanceConverter &dconv,
    size_t n,
    uint8_t do_binary_search = 3
  ) : ClusterAlgorithm(dconv, n),
  clust_array(m*n),
  ca(&clust_array[0])
  {
    initialize();
  };

#ifdef OPTIMOTU_R
  ClusterIndexedMatrix(const DistanceConverter &dconv, Rcpp::IntegerMatrix &source) = delete;
#endif

  ~ClusterIndexedMatrix();

  using ClusterAlgorithm::operator();

  void dump_index();

  void print_index();

  void verify_index();

  void heal_splice();

  bool index_splice(tip *&t1max, tip *&t2min, tip *&t2max, d_t i);

  void operator()(j_t seq1, j_t seq2, d_t i) override;

  void merge_into(DistanceConsumer &consumer) override;

  void merge_into(ClusterAlgorithm &consumer) override;

  ClusterAlgorithm * make_child() override;

  double max_relevant(j_t seq1, j_t seq2) const override;

#ifdef OPTIMOTU_R
  void write_to_matrix(RcppParallel::RMatrix<int> &out) override;
#endif
};

template class ClusterIndexedMatrix<>;

#ifdef OPTIMOTU_R

template<>
ClusterIndexedMatrix<RcppParallel::RMatrix<int>>::ClusterIndexedMatrix(
    const DistanceConverter &dconv,
    size_t n,
    uint8_t do_binary_search) = delete;

template<>
ClusterIndexedMatrix<RcppParallel::RMatrix<int>>::ClusterIndexedMatrix(
    ClusterAlgorithm * parent
) = delete;

template<>
ClusterIndexedMatrix<RcppParallel::RMatrix<int>>::ClusterIndexedMatrix(
  const DistanceConverter &dconv,
  Rcpp::IntegerMatrix &im
);

template class ClusterIndexedMatrix<RcppParallel::RMatrix<int>>;
#endif

#endif //OPTIMOTU_CLUSTERINDEXEDMATRIX_H_INCLUDED
