// SPDX-FileCopyrightText: 2025 Brendan Furneaux <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_CLUSTERSLINK_H_INCLUDED
#define OPTIMOTU_CLUSTERSLINK_H_INCLUDED

#include "single_linkage.h"
#include "ClusterAlgorithm.h"
#include "MappedClusterAlgorithm.h"
#include <memory>
#include <vector>

#ifdef OPTIMOTU_R
#include <Rcpp.h>
#endif

template <int verbose = 0>
class ClusterSLINK : public SingleClusterAlgorithm {
  template<int> friend class ClusterSLINK;
  friend class MappedClusterAlgorithm;
protected:
  struct SlinkState
  {
    std::vector<j_t> Pi;     // first higher-numbered leaf which is joined
    std::vector<d_t> Lambda; // level at which Pi[i] is joined
    std::vector<d_t> M;
    j_t slink_seq1 = 0;
    j_t slink_seq2 = 0;
  };
  struct MergeEdge
  {
    j_t seq1;
    j_t seq2;
    d_t i;
  };
  std::unique_ptr<SlinkState> slink;
  // Unordered MST edges from children; sorted and replayed in finalize().
  // Set when any child has been created. Cannot be derived from `children`
  // because children are erased by release_child() after they merge.
  std::vector<MergeEdge> merge_cache;
  bool uses_cache = false;

  void ensure_slink();
  void apply_ordered(j_t seq1, j_t seq2, d_t i);
  void init_iter();
  void update();
  void finish_iter();

  ClusterSLINK(ClusterAlgorithm * parent, j_t n);

  ClusterSLINK<verbose> * make_inner_child(ClusterAlgorithm * parent, const j_t n) override;

public:
  ClusterSLINK(const DistanceConverter &dconv, const j_t n);
  ClusterSLINK(const DistanceConverter &dconv, init_matrix_t im);

  ClusterSLINK<verbose> * make_child() override;
  MappedClusterAlgorithm * make_child(PairGenerator * pg) override;
  ClusterAlgorithm *make_child(const std::vector<std::size_t> &fwd_map) override;

  virtual void operator()(j_t seq1, j_t seq2, d_t i, int thread = 0) override;
  virtual void operator()(PairGenerator & pg, d_t i, int thread = 0) override;

  virtual void finalize() override;

  void write_to_matrix(internal_matrix_t &out) override;

  void prepare_output() override;

  void write_threshold_row(d_t t, int *dest) const override;

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
  virtual double max_relevant(PairGenerator & pg, int thread = 0) const override;

  bool accepts_unordered_pairs() const override {
    return uses_cache;
  }
};

std::unique_ptr<SingleClusterAlgorithm> create_cluster_slink(
  const DistanceConverter & dconv, j_t n, int verbose = 0);
std::unique_ptr<SingleClusterAlgorithm> create_cluster_slink(
  const DistanceConverter & dconv, init_matrix_t & im, int verbose = 0);

#endif
