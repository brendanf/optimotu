// SPDX-FileCopyrightText: 2026 Brendan Furneaux <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_CLUSTERPRIM_H_INCLUDED
#define OPTIMOTU_CLUSTERPRIM_H_INCLUDED

#include "single_linkage.h"
#include "ClusterAlgorithm.h"
#include "MappedClusterAlgorithm.h"
#include "optimotu.h"

#include <atomic>
#include <condition_variable>
#include <limits>
#include <memory>
#include <mutex>
#include <vector>

#ifdef OPTIMOTU_R
#include <Rcpp.h>
#include <RcppThread.h>
#endif

class ClusterPrim : public SingleClusterAlgorithm {
  friend class MappedClusterAlgorithm;

  struct TournamentNode {
    std::vector<MstEdge> edges;
    std::atomic<int> arrived{0};
  };

  std::vector<MstEdge> mst_;
  bool grown_ = false;
  std::size_t merge_slot_ = 0;

  std::size_t n_leaves_ = 0;
  std::size_t n_pad_ = 0;
  std::vector<TournamentNode> nodes_;
  std::atomic<bool> root_ready_{false};
  std::mutex ready_mu_;
  std::condition_variable ready_cv_;

  ClusterPrim(ClusterAlgorithm *parent, j_t n);

  ClusterPrim *make_inner_child(ClusterAlgorithm *parent, const j_t n) override;

  double key_cutoff(d_t key) const;
  void add_mst_edge(j_t u, j_t v, d_t bin);
  void kruskal_from_edges(std::vector<MstEdge> edges);
  void submit_leaf(std::size_t leaf_id, std::vector<MstEdge> edges);
  void wait_root();
  void union_until(d_t t, std::vector<j_t> &parent) const;
  j_t uf_find(std::vector<j_t> &parent, j_t x) const;

  template <typename Dist>
  void grow_from_start(j_t n_vert, Dist dist, bool bipartite, j_t n_j);

public:
  ClusterPrim(const DistanceConverter &dconv, const j_t n);
  ClusterPrim(const DistanceConverter &dconv, init_matrix_t im);

  ClusterPrim *make_child() override;
  MappedClusterAlgorithm *make_child(PairGenerator *pg) override;
  ClusterAlgorithm *make_child(const std::vector<std::size_t> &fwd_map) override;

  void operator()(j_t seq1, j_t seq2, d_t i, int thread = 0) override;

  void finalize() override;

  void write_to_matrix(internal_matrix_t &out) override;
  void write_threshold_row(d_t t, int *dest) const override;

#ifdef OPTIMOTU_R
  Rcpp::List as_hclust(const Rcpp::CharacterVector &seqnames) const override;
#endif

  void merge_into(DistanceConsumer &consumer) override;
  void merge_into(ClusterAlgorithm &consumer) override;
  void merge_into_parent() override;

  double max_relevant(j_t seq1, j_t seq2, int thread = 0) const override;

  bool accepts_unordered_pairs() const override { return true; }
  bool uses_prim() const override { return true; }
  void prepare_merge_slots(std::size_t n_slots) override;
  void set_merge_slot(std::size_t slot) override;
  void receive_mst(std::size_t slot, std::vector<MstEdge> edges) override;

  template <typename Dist>
  void grow_complete(Dist dist);

  template <typename Dist>
  void grow_bipartite(j_t n_j, j_t n_i, Dist dist);
};

std::unique_ptr<SingleClusterAlgorithm> create_cluster_prim(
    const DistanceConverter &dconv, j_t n, int verbose = 0);
std::unique_ptr<SingleClusterAlgorithm> create_cluster_prim(
    const DistanceConverter &dconv, init_matrix_t &im, int verbose = 0);

template <typename Dist>
void ClusterPrim::grow_from_start(
    j_t n_vert, Dist dist, bool bipartite, j_t n_j)
{
  grown_ = true;
  mst_.clear();
  if (n_vert < 2)
  {
    return;
  }
  const d_t inf = this->m;
  std::vector<d_t> key(n_vert, inf);
  std::vector<j_t> parent(n_vert, 0);
  std::vector<char> in_tree(n_vert, 0);
  std::size_t dist_calls = 0;
  auto maybe_interrupt = [&dist_calls]()
  {
#ifdef OPTIMOTU_R
    if ((++dist_calls & 65535u) == 0u)
    {
      RcppThread::checkUserInterrupt();
    }
#else
    (void)dist_calls;
#endif
  };

  auto update_key = [&](j_t u, j_t v)
  {
    if (in_tree[v])
    {
      return;
    }
    const double cutoff = key_cutoff(key[v]);
    maybe_interrupt();
    const double d = dist(u, v, cutoff);
    const d_t bin = this->dconv.convert(d);
    if (bin < key[v])
    {
      key[v] = bin;
      parent[v] = u;
    }
  };

  auto opposite = [&](j_t u, j_t v)
  {
    if (!bipartite)
    {
      return true;
    }
    const bool u_left = u < n_j;
    const bool v_left = v < n_j;
    return u_left != v_left;
  };

  in_tree[0] = 1;
  key[0] = 0;
  for (j_t v = 1; v < n_vert; ++v)
  {
    if (opposite(0, v))
    {
      update_key(0, v);
    }
  }

  for (j_t added = 1; added < n_vert; ++added)
  {
    j_t u = NO_CLUST;
    d_t best = inf;
    for (j_t v = 0; v < n_vert; ++v)
    {
      if (!in_tree[v] && key[v] < best)
      {
        best = key[v];
        u = v;
      }
    }
    // Distances at or beyond dmax never lower a key below inf, so the
    // graph may be a forest. Seed the next tree; do not add a dummy
    // edge to the previous component.
    if (u == NO_CLUST)
    {
      for (j_t v = 0; v < n_vert; ++v)
      {
        if (!in_tree[v])
        {
          u = v;
          break;
        }
      }
      if (u == NO_CLUST)
      {
        break;
      }
      in_tree[u] = 1;
      key[u] = 0;
      for (j_t v = 0; v < n_vert; ++v)
      {
        if (opposite(u, v))
        {
          update_key(u, v);
        }
      }
      continue;
    }
    in_tree[u] = 1;
    add_mst_edge(parent[u], u, key[u]);
    for (j_t v = 0; v < n_vert; ++v)
    {
      if (opposite(u, v))
      {
        update_key(u, v);
      }
    }
  }
}

template <typename Dist>
void ClusterPrim::grow_complete(Dist dist)
{
  grow_from_start(this->n, dist, false, 0);
}

template <typename Dist>
void ClusterPrim::grow_bipartite(j_t n_j, j_t n_i, Dist dist)
{
  (void)n_i;
  grow_from_start(this->n, dist, true, n_j);
}

#endif
