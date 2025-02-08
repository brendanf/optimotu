// SPDX-FileCopyrightText: 2025 Brendan Furneaux <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_CLUSTERSLINK_H_INCLUDED
#define OPTIMOTU_CLUSTERSLINK_H_INCLUDED

#include "single_linkage.h"
#include "ClusterAlgorithm.h"
#include "ClusterTree.h"

#ifdef OPTIMOTU_R
#include <Rcpp.h>
#endif

class ClusterSLINK : public SingleClusterAlgorithm {
protected:
  std::vector<j_t> Pi; // first higher-numbered leaf which is joined
  std::vector<d_t> Lambda; // level at which Pi[i] is joined
  std::vector<d_t> M;
  j_t slink_seq1 = 0;
  j_t slink_seq2 = 0;
  ClusterTree delegate;

  void init_iter();
  void update();
  void finish_iter();

  ClusterSLINK(SingleClusterAlgorithm * parent);

public:
  ClusterSLINK(const DistanceConverter &dconv, const j_t n);
  ClusterSLINK(const DistanceConverter &dconv, init_matrix_t im);

  ClusterSLINK * make_child() override;

  virtual void operator()(j_t seq1, j_t seq2, d_t i, int thread = 0) override;

  virtual void finalize() override;

  void write_to_matrix(internal_matrix_t &out) override;

#ifdef OPTIMOTU_R
  Rcpp::List as_hclust(const Rcpp::CharacterVector &seqnames) const override;
#endif // OPTIMOTU_R

  // send consumer() pairwise distances to ensure it is up-to-date with this
  // clustering
  virtual void merge_into(DistanceConsumer &consumer) override;

  // send consumer() pairwise distances to ensure it is up-to-date with this
  // clustering
  virtual void merge_into(ClusterAlgorithm &consumer) override;

  virtual void merge_into_parent() override;

  // calculate the maximum distance between seq1 and seq2 which would actually
  // cause an update
  virtual double max_relevant(j_t seq1, j_t seq2, int thread = 0) const override;
};

#endif
