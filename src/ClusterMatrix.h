// SPDX-FileCopyrightText: 2025 Brendan Furneaux <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_CLUSTERMATRIX_H_INCLUDED
#define OPTIMOTU_CLUSTERMATRIX_H_INCLUDED

#include "ClusterAlgorithm.h"
#include "MappedClusterAlgorithm.h"
#include <memory>
#ifdef OPTIMOTU_R
#include <RcppParallel.h>
#endif // OPTIMOTU_R

#define LINEAR_FILL 1
#define BINARY_FILL 2
#define TOPDOWN_FILL 3

// Matrix representation of a hierarchical clustering of n items at m different
// thresholds.
template <bool BINARY_SEARCH=true, //alt would be linear
          int FILL=2,
          class ARRAY_T = std::vector<int>,
          int verbose = 0>
class ClusterMatrix : public SingleClusterAlgorithm {
  template<bool, int, typename, int> friend class ClusterMatrix;
  friend class MappedClusterAlgorithm;
private:
  ARRAY_T clust_array;
  int *const ca;
  std::vector<int> toclust;

protected:
  void initialize();

  ClusterMatrix(ClusterAlgorithm * parent, const j_t n);

  SingleClusterAlgorithm * make_inner_child(
    ClusterAlgorithm * parent,
    const j_t n
  ) override;

public:
  ClusterMatrix(const DistanceConverter &dconv, size_t n);

  ClusterMatrix(const DistanceConverter &dconv, init_matrix_t im);

  using ClusterAlgorithm::operator();

  void operator()(j_t seq1, j_t seq2, d_t i, int thread = 0) override;

  void merge_into(DistanceConsumer &consumer) override final;

  void merge_into(ClusterAlgorithm &consumer) override final;

  SingleClusterAlgorithm * make_child() override;
  MappedClusterAlgorithm * make_child(PairGenerator * pg) override;
  ClusterAlgorithm *make_child(const std::vector<std::size_t> &fwd_map) override;

  double max_relevant(j_t seq1, j_t seq2, int thread = 0) const override;

  void write_to_matrix(internal_matrix_t &out) override;

  void write_threshold_row(d_t t, int *dest) const override;

#ifdef OPTIMOTU_R
  Rcpp::List as_hclust(const Rcpp::CharacterVector &seqnames) const override;
#endif // OPTIMOTU_R
};

std::unique_ptr<SingleClusterAlgorithm> create_cluster_matrix(
  const DistanceConverter & dconv, j_t n,
  bool binary_search, int fill_type, int verbose = 0);
std::unique_ptr<SingleClusterAlgorithm> create_cluster_matrix(
  const DistanceConverter & dconv, init_matrix_t & im,
  bool binary_search, int fill_type);

#endif //OPTIMOTU_CLUSTERMATRIX_H_INCLUDED
