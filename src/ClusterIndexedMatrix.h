#ifndef OPTIMOTU_CLUSTERINDEXEDMATRIX_H_INCLUDED
#define OPTIMOTU_CLUSTERINDEXEDMATRIX_H_INCLUDED

#include "ClusterAlgorithm.h"
#ifdef OPTIMOTU_R
#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#endif // OPTIMOTU_R

// Matrix representation of a hierarchical clustering of n items at m different
// thresholds.
template <class ARRAY_T = std::vector<int>>
class ClusterIndexedMatrix : public SingleClusterAlgorithm {
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

  ClusterIndexedMatrix(SingleClusterAlgorithm * parent);

public:
  ClusterIndexedMatrix(const DistanceConverter &dconv, size_t n);

  ClusterIndexedMatrix(const DistanceConverter &dconv, init_matrix_t &im);

  ~ClusterIndexedMatrix();

  using ClusterAlgorithm::operator();

  void dump_index();

  void print_index();

  void verify_index();

  void heal_splice();

  bool index_splice(tip *&t1max, tip *&t2min, tip *&t2max, d_t i);

  void operator()(j_t seq1, j_t seq2, d_t i, int thread = 0) override;

  void merge_into(DistanceConsumer &consumer) override;

  void merge_into(ClusterAlgorithm &consumer) override;

  SingleClusterAlgorithm * make_child() override;

  double max_relevant(j_t seq1, j_t seq2, int thread = 0) const override;

  void write_to_matrix(internal_matrix_t &out) override;

#ifdef OPTIMOTU_R
  Rcpp::List as_hclust(const Rcpp::CharacterVector &seqnames) const override;
#endif // OPTIMOTU_R
};

template<>
ClusterIndexedMatrix<>::ClusterIndexedMatrix(
  const DistanceConverter &dconv,
  init_matrix_t &im
) = delete;

template class ClusterIndexedMatrix<>;

#define cim_internal ClusterIndexedMatrix<internal_matrix_ref_t>

template<>
cim_internal::ClusterIndexedMatrix(
    const DistanceConverter &dconv,
    size_t n
    ) = delete;

template<>
cim_internal::ClusterIndexedMatrix(
    SingleClusterAlgorithm * parent
) = delete;

template class cim_internal;

#undef cim_internal

#endif //OPTIMOTU_CLUSTERINDEXEDMATRIX_H_INCLUDED
