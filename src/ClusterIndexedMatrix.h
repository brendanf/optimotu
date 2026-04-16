// SPDX-FileCopyrightText: 2025 Brendan Furneaux <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_CLUSTERINDEXEDMATRIX_H_INCLUDED
#define OPTIMOTU_CLUSTERINDEXEDMATRIX_H_INCLUDED

#include "ClusterAlgorithm.h"
#include "MappedClusterAlgorithm.h"
#include <memory>
#ifdef OPTIMOTU_R
#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#endif // OPTIMOTU_R

// Matrix representation of a hierarchical clustering of n items at m different
// thresholds.
template <class ARRAY_T = std::vector<int>, int verbose = 0>
class ClusterIndexedMatrix : public SingleClusterAlgorithm {
  template<typename, int> friend class ClusterIndexedMatrix;

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

  ClusterIndexedMatrix(ClusterAlgorithm * parent, const j_t n);

  SingleClusterAlgorithm * make_inner_child(
    ClusterAlgorithm * parent,
    const j_t n
  ) override;

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

  MappedClusterAlgorithm * make_child(PairGenerator * pg) override;

  double max_relevant(j_t seq1, j_t seq2, int thread = 0) const override;

  void write_to_matrix(internal_matrix_t &out) override;

#ifdef OPTIMOTU_R
  Rcpp::List as_hclust(const Rcpp::CharacterVector &seqnames) const override;
#endif // OPTIMOTU_R
};

template<>
ClusterIndexedMatrix<std::vector<int>, 0>::ClusterIndexedMatrix(
  const DistanceConverter &dconv,
  init_matrix_t &im
) = delete;

template<>
ClusterIndexedMatrix<std::vector<int>, 1>::ClusterIndexedMatrix(
  const DistanceConverter &dconv,
  init_matrix_t &im
) = delete;

template<>
ClusterIndexedMatrix<std::vector<int>, 2>::ClusterIndexedMatrix(
  const DistanceConverter &dconv,
  init_matrix_t &im
) = delete;

#define cim_internal ClusterIndexedMatrix<internal_matrix_ref_t, 0>

template<>
cim_internal::ClusterIndexedMatrix(
    const DistanceConverter &dconv,
    size_t n
    ) = delete;

template<>
cim_internal::ClusterIndexedMatrix(
    ClusterAlgorithm * parent,
    const j_t n
) = delete;

#undef cim_internal

std::unique_ptr<SingleClusterAlgorithm> create_cluster_indexed_matrix(
  const DistanceConverter & dconv, j_t n, int verbose = 0);
std::unique_ptr<SingleClusterAlgorithm> create_cluster_indexed_matrix(
  const DistanceConverter & dconv, init_matrix_t & im);

#endif //OPTIMOTU_CLUSTERINDEXEDMATRIX_H_INCLUDED
