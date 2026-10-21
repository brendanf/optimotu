// SPDX-FileCopyrightText: 2026 Brendan Furneaux <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_KRUSKAL_TOURNAMENT_H_INCLUDED
#define OPTIMOTU_KRUSKAL_TOURNAMENT_H_INCLUDED

#include "ClusterAlgorithm.h"
#include "optimotu.h"

#include <algorithm>
#include <atomic>
#include <cassert>
#include <condition_variable>
#include <iterator>
#include <mutex>
#include <numeric>
#include <vector>

inline std::vector<MstEdge> kruskal_union(
    const std::vector<MstEdge> &a,
    const std::vector<MstEdge> &b,
    j_t n)
{
  if (a.empty())
  {
    return b;
  }
  if (b.empty())
  {
    return a;
  }
  assert(std::is_sorted(a.begin(), a.end()));
  assert(std::is_sorted(b.begin(), b.end()));
  std::vector<MstEdge> all;
  all.reserve(a.size() + b.size());
  std::merge(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(all));
  std::vector<j_t> parent(n);
  std::vector<j_t> rank(n, 0);
  std::iota(parent.begin(), parent.end(), 0);
  auto find = [&](j_t x)
  {
    while (parent[x] != x)
    {
      parent[x] = parent[parent[x]];
      x = parent[x];
    }
    return x;
  };
  auto unite = [&](j_t x, j_t y)
  {
    x = find(x);
    y = find(y);
    if (x == y)
    {
      return false;
    }
    if (rank[x] < rank[y])
    {
      parent[x] = y;
    }
    else if (rank[y] < rank[x])
    {
      parent[y] = x;
    }
    else
    {
      parent[y] = x;
      ++rank[x];
    }
    return true;
  };
  std::vector<MstEdge> out;
  out.reserve(n > 0 ? n - 1 : 0);
  for (MstEdge e : all)
  {
    if (e.u > e.v)
    {
      std::swap(e.u, e.v);
    }
    if (e.u == e.v || e.u >= n || e.v >= n)
    {
      continue;
    }
    if (unite(e.u, e.v))
    {
      out.push_back(e);
      if (n > 0 && out.size() + 1 == static_cast<std::size_t>(n))
      {
        break;
      }
    }
  }
  return out;
}

inline j_t mst_uf_find(std::vector<j_t> &parent, j_t x)
{
  while (parent[x] != x)
  {
    parent[x] = parent[parent[x]];
    x = parent[x];
  }
  return x;
}

inline void mst_union_until(
    const std::vector<MstEdge> &mst,
    d_t t,
    std::vector<j_t> &parent)
{
  std::iota(parent.begin(), parent.end(), 0);
  for (const MstEdge &e : mst)
  {
    if (e.bin > t)
    {
      continue;
    }
    const j_t a = mst_uf_find(parent, e.u);
    const j_t b = mst_uf_find(parent, e.v);
    if (a != b)
    {
      parent[a] = b;
    }
  }
}

inline void mst_write_to_matrix(
    const std::vector<MstEdge> &mst,
    j_t n,
    d_t m,
    internal_matrix_t &out)
{
  std::vector<std::vector<int>> by_t(m, std::vector<int>(n));
  std::vector<j_t> parent(n);
  for (d_t t = 0; t < m; ++t)
  {
    mst_union_until(mst, t, parent);
    for (j_t i = 0; i < n; ++i)
    {
      by_t[t][i] = static_cast<int>(mst_uf_find(parent, i));
    }
  }
  std::size_t k = 0;
  for (j_t i = 0; i < n; ++i)
  {
    for (d_t t = 0; t < m; ++t)
    {
      out[k++] = by_t[t][i];
    }
  }
}

inline void mst_write_threshold_row(
    const std::vector<MstEdge> &mst,
    j_t n,
    d_t t,
    int *dest)
{
  std::vector<j_t> parent(n);
  mst_union_until(mst, t, parent);
  for (j_t i = 0; i < n; ++i)
  {
    dest[i] = static_cast<int>(mst_uf_find(parent, i));
  }
}

#ifdef OPTIMOTU_R
inline Rcpp::List mst_as_hclust(
    const std::vector<MstEdge> &mst,
    j_t n,
    const DistanceConverter &dconv,
    const Rcpp::CharacterVector &seqnames)
{
  Rcpp::IntegerMatrix merge(std::max(0, static_cast<int>(n) - 1), 2);
  Rcpp::NumericVector height(std::max(0, static_cast<int>(n) - 1));
  Rcpp::IntegerVector order(n);
  if (n == 0)
  {
    Rcpp::List out;
    out["merge"] = merge;
    out["height"] = height;
    out["order"] = order;
    out["labels"] = seqnames;
    out["method"] = "single";
    out.attr("class") = "hclust";
    return out;
  }
  std::vector<j_t> parent(n);
  std::vector<j_t> rank(n, 0);
  std::iota(parent.begin(), parent.end(), 0);
  std::vector<int> label(n);
  for (j_t i = 0; i < n; ++i)
  {
    label[i] = -static_cast<int>(i + 1);
  }
  auto find = [&](j_t x)
  {
    while (parent[x] != x)
    {
      parent[x] = parent[parent[x]];
      x = parent[x];
    }
    return x;
  };
  int row = 0;
  std::vector<MstEdge> edges = mst;
  std::sort(edges.begin(), edges.end());
  for (MstEdge e : edges)
  {
    if (row >= static_cast<int>(n) - 1)
    {
      break;
    }
    j_t a = find(e.u);
    j_t b = find(e.v);
    if (a == b)
    {
      continue;
    }
    merge(row, 0) = label[a];
    merge(row, 1) = label[b];
    height[row] = dconv.inverse(e.bin);
    if (rank[a] < rank[b])
    {
      parent[a] = b;
      label[b] = row + 1;
    }
    else
    {
      parent[b] = a;
      label[a] = row + 1;
      if (rank[a] == rank[b])
      {
        ++rank[a];
      }
    }
    ++row;
  }
  while (row < static_cast<int>(n) - 1)
  {
    j_t ra = find(0);
    j_t other = n;
    for (j_t i = 1; i < n; ++i)
    {
      if (find(i) != ra)
      {
        other = i;
        break;
      }
    }
    if (other >= n)
    {
      break;
    }
    j_t rb = find(other);
    merge(row, 0) = label[ra];
    merge(row, 1) = label[rb];
    height[row] = 1.0;
    parent[rb] = ra;
    label[ra] = row + 1;
    ++row;
  }
  for (j_t i = 0; i < n; ++i)
  {
    order[i] = static_cast<int>(i + 1);
  }
  Rcpp::List out;
  out["merge"] = merge;
  out["height"] = height;
  out["order"] = order;
  out["labels"] = seqnames;
  out["method"] = "single";
  out.attr("class") = "hclust";
  return out;
}
#endif

class KruskalTournament {
  struct Node {
    std::vector<MstEdge> edges;
    std::atomic<int> arrived{0};
  };

  j_t n_ = 0;
  std::size_t n_pad_ = 0;
  std::vector<Node> nodes_;
  std::vector<MstEdge> mst_;
  std::atomic<bool> root_ready_{false};
  std::mutex ready_mu_;
  std::condition_variable ready_cv_;

  void submit_leaf(std::size_t leaf_id, std::vector<MstEdge> edges)
  {
    for (MstEdge &e : edges)
    {
      if (e.u > e.v)
      {
        std::swap(e.u, e.v);
      }
    }
    std::sort(edges.begin(), edges.end());
    std::size_t idx = n_pad_ + leaf_id;
    nodes_[idx].edges = std::move(edges);
    while (idx > 1)
    {
      const std::size_t par = idx / 2;
      const int prev =
          nodes_[par].arrived.fetch_add(1, std::memory_order_acq_rel);
      if (prev == 0)
      {
        return;
      }
      nodes_[par].edges = kruskal_union(
          nodes_[2 * par].edges, nodes_[2 * par + 1].edges, n_);
      nodes_[2 * par].edges.clear();
      nodes_[2 * par].edges.shrink_to_fit();
      nodes_[2 * par + 1].edges.clear();
      nodes_[2 * par + 1].edges.shrink_to_fit();
      idx = par;
    }
    mst_ = nodes_[1].edges;
    root_ready_.store(true, std::memory_order_release);
    ready_cv_.notify_all();
  }

public:
  bool prepared() const { return n_pad_ > 0; }

  const std::vector<MstEdge> &mst() const { return mst_; }

  void prepare(std::size_t n_slots, j_t n)
  {
    n_ = n;
    n_pad_ = 1;
    while (n_pad_ < n_slots)
    {
      n_pad_ *= 2;
    }
    {
      std::vector<Node> fresh(2 * n_pad_);
      nodes_.swap(fresh);
    }
    mst_.clear();
    root_ready_.store(n_slots == 0, std::memory_order_release);
    for (std::size_t i = n_slots; i < n_pad_; ++i)
    {
      submit_leaf(i, {});
    }
  }

  void submit(std::size_t slot, std::vector<MstEdge> edges)
  {
    if (n_pad_ == 0)
    {
      OPTIMOTU_STOP(
          "merge tournament was not prepared before receive_mst()");
    }
    submit_leaf(slot, std::move(edges));
  }

  void wait()
  {
    std::unique_lock<std::mutex> lock(ready_mu_);
    ready_cv_.wait(lock, [this]
                   { return root_ready_.load(std::memory_order_acquire); });
  }
};

// Peak for a parent KruskalTournament. MstEdge is 16B padded; Node is
// vector plus atomic (~32B). Leaves may each hold an n-1 forest while
// last-arrivers run kruskal_union (merge buffer, out, union-find).
inline std::size_t kruskal_tournament_estimate_bytes(
    j_t n,
    std::size_t n_slots)
{
  const std::size_t ns = static_cast<std::size_t>(n);
  constexpr std::size_t k_edge = 16;
  constexpr std::size_t k_node = 32;
  if (n_slots == 0)
  {
    n_slots = 1;
  }
  std::size_t n_pad = 1;
  while (n_pad < n_slots)
  {
    n_pad *= 2;
  }
  const std::size_t skeleton = 2 * n_pad * k_node;
  const std::size_t leaf_edges = n_slots * ns * k_edge;
  const std::size_t conc = std::max<std::size_t>(1, n_slots / 2);
  const std::size_t union_scratch =
      conc * (ns * 8 + 3 * ns * k_edge);
  const std::size_t final_mst = ns * k_edge;
  return skeleton + leaf_edges + union_scratch + final_mst;
}

#endif
