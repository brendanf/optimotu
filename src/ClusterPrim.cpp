// SPDX-FileCopyrightText: 2026 Brendan Furneaux <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#include "ClusterPrim.h"
#include "optimotu.h"

#ifdef OPTIMOTU_R
#include <RcppThread.h>
#endif

#include <algorithm>
#include <numeric>

namespace
{

  std::vector<MstEdge> kruskal_union(
      const std::vector<MstEdge> &a,
      const std::vector<MstEdge> &b,
      j_t n)
  {
    std::vector<MstEdge> all;
    all.reserve(a.size() + b.size());
    all.insert(all.end(), a.begin(), a.end());
    all.insert(all.end(), b.begin(), b.end());
    std::sort(all.begin(), all.end());
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

} // namespace

ClusterPrim::ClusterPrim(ClusterAlgorithm *parent, j_t n)
    : SingleClusterAlgorithm(parent, n) {}

ClusterPrim::ClusterPrim(const DistanceConverter &dconv, const j_t n)
    : SingleClusterAlgorithm(dconv, n) {}

ClusterPrim::ClusterPrim(const DistanceConverter &dconv, init_matrix_t im)
    : SingleClusterAlgorithm(dconv, im) {}

double ClusterPrim::key_cutoff(d_t key) const
{
  if (key <= 0)
  {
    return 0.0;
  }
  return dconv.inverse(key - 1);
}

void ClusterPrim::add_mst_edge(j_t u, j_t v, d_t bin)
{
  if (u > v)
  {
    std::swap(u, v);
  }
  mst_.push_back({u, v, bin});
}

void ClusterPrim::kruskal_from_edges(std::vector<MstEdge> edges)
{
  mst_ = kruskal_union(edges, {}, n);
}

j_t ClusterPrim::uf_find(std::vector<j_t> &parent, j_t x) const
{
  while (parent[x] != x)
  {
    parent[x] = parent[parent[x]];
    x = parent[x];
  }
  return x;
}

void ClusterPrim::union_until(d_t t, std::vector<j_t> &parent) const
{
  std::iota(parent.begin(), parent.end(), 0);
  for (const MstEdge &e : mst_)
  {
    if (e.bin > t)
    {
      continue;
    }
    const j_t a = uf_find(parent, e.u);
    const j_t b = uf_find(parent, e.v);
    if (a != b)
    {
      parent[a] = b;
    }
  }
}

void ClusterPrim::submit_leaf(std::size_t leaf_id, std::vector<MstEdge> edges)
{
  if (n_pad_ == 0)
  {
    kruskal_from_edges(std::move(edges));
    root_ready_.store(true, std::memory_order_release);
    ready_cv_.notify_all();
    return;
  }
  std::size_t idx = n_pad_ + leaf_id;
  nodes_[idx].edges = std::move(edges);
  while (idx > 1)
  {
    const std::size_t par = idx / 2;
    const int prev = nodes_[par].arrived.fetch_add(1, std::memory_order_acq_rel);
    if (prev == 0)
    {
      return;
    }
    nodes_[par].edges = kruskal_union(
        nodes_[2 * par].edges, nodes_[2 * par + 1].edges, n);
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

void ClusterPrim::wait_root()
{
  std::unique_lock<std::mutex> lock(ready_mu_);
  ready_cv_.wait(lock, [this]
                 { return root_ready_.load(std::memory_order_acquire); });
}

void ClusterPrim::prepare_merge_slots(std::size_t n_slots)
{
  n_leaves_ = n_slots;
  n_pad_ = 1;
  while (n_pad_ < n_slots)
  {
    n_pad_ *= 2;
  }
  nodes_ = std::vector<TournamentNode>(2 * n_pad_);
  root_ready_.store(n_slots == 0, std::memory_order_release);
  for (std::size_t i = n_slots; i < n_pad_; ++i)
  {
    submit_leaf(i, {});
  }
}

void ClusterPrim::set_merge_slot(std::size_t slot)
{
  merge_slot_ = slot;
}

void ClusterPrim::receive_mst(std::size_t slot, std::vector<MstEdge> edges)
{
  if (n_pad_ == 0)
  {
    OPTIMOTU_STOP(
        "clust_prim() merge tournament was not prepared before receive_mst()");
  }
  submit_leaf(slot, std::move(edges));
}

ClusterPrim *ClusterPrim::make_inner_child(ClusterAlgorithm *parent, const j_t n)
{
  auto child_ptr = new ClusterPrim(parent, n);
  auto child = std::unique_ptr<ClusterAlgorithm>(
      (ClusterAlgorithm *)child_ptr);
  this->children.push_back(std::move(child));
  return child_ptr;
}

ClusterPrim *ClusterPrim::make_child()
{
  std::unique_lock<std::shared_timed_mutex> lock(this->mutex);
  auto child_ptr = new ClusterPrim(this, n);
  auto child = std::unique_ptr<ClusterAlgorithm>(
      (ClusterAlgorithm *)child_ptr);
  this->children.push_back(std::move(child));
  return child_ptr;
}

MappedClusterAlgorithm *ClusterPrim::make_child(PairGenerator *pg)
{
  std::unique_lock<std::shared_timed_mutex> lock(this->mutex);
  auto child_ptr = new MappedClusterAlgorithmImpl<0>(this, pg);
  auto child = std::unique_ptr<ClusterAlgorithm>(
      (ClusterAlgorithm *)child_ptr);
  this->children.push_back(std::move(child));
  return child_ptr;
}

ClusterAlgorithm *ClusterPrim::make_child(const std::vector<std::size_t> &fwd_map)
{
  if (fwd_map.size() < 2)
  {
    return nullptr;
  }
  std::unique_lock<std::shared_timed_mutex> lock(this->mutex);
  auto rev_map = invert_map(fwd_map);
  auto child_ptr = new MappedClusterAlgorithmImpl<0>(this, fwd_map, rev_map);
  auto child = std::unique_ptr<ClusterAlgorithm>(
      (ClusterAlgorithm *)child_ptr);
  this->children.push_back(std::move(child));
  return child_ptr;
}

void ClusterPrim::operator()(j_t seq1, j_t seq2, d_t i, int thread)
{
  (void)seq1;
  (void)seq2;
  (void)i;
  (void)thread;
  OPTIMOTU_STOP(
      "clust_prim() cannot consume streamed pairwise distances");
}

void ClusterPrim::finalize()
{
  if (n_pad_ > 0 && parent == nullptr)
  {
    wait_root();
    return;
  }
  if (grown_)
  {
    std::sort(mst_.begin(), mst_.end());
  }
}

void ClusterPrim::write_to_matrix(internal_matrix_t &out)
{
  std::shared_lock<std::shared_timed_mutex> lock(this->mutex);
  std::vector<std::vector<int>> by_t(m, std::vector<int>(n));
  std::vector<j_t> parent(n);
  for (d_t t = 0; t < m; ++t)
  {
    union_until(t, parent);
    for (j_t i = 0; i < n; ++i)
    {
      by_t[t][i] = static_cast<int>(uf_find(parent, i));
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

void ClusterPrim::write_threshold_row(d_t t, int *dest) const
{
  std::shared_lock<std::shared_timed_mutex> lock(this->mutex);
  std::vector<j_t> parent(n);
  union_until(t, parent);
  for (j_t i = 0; i < n; ++i)
  {
    dest[i] = static_cast<int>(uf_find(parent, i));
  }
}

#ifdef OPTIMOTU_R
Rcpp::List ClusterPrim::as_hclust(const Rcpp::CharacterVector &seqnames) const
{
  std::shared_lock<std::shared_timed_mutex> lock(this->mutex);
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
  std::vector<MstEdge> edges = mst_;
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

void ClusterPrim::merge_into(DistanceConsumer &consumer)
{
  std::shared_lock<std::shared_timed_mutex> lock(this->mutex);
  for (const MstEdge &e : mst_)
  {
    consumer(e.u, e.v, dconv.inverse(e.bin));
  }
}

void ClusterPrim::merge_into(ClusterAlgorithm &consumer)
{
  std::shared_lock<std::shared_timed_mutex> lock(this->mutex);
  for (const MstEdge &e : mst_)
  {
    consumer(e.u, e.v, e.bin);
  }
}

void ClusterPrim::merge_into_parent()
{
  if (!parent)
  {
    return;
  }
  parent->receive_mst(merge_slot_, mst_);
}

double ClusterPrim::max_relevant(j_t seq1, j_t seq2, int thread) const
{
  (void)seq1;
  (void)seq2;
  (void)thread;
  return dconv.inverse(this->m - 1);
}

std::unique_ptr<SingleClusterAlgorithm> create_cluster_prim(
    const DistanceConverter &dconv, j_t n, int)
{
  return std::make_unique<ClusterPrim>(dconv, n);
}

std::unique_ptr<SingleClusterAlgorithm> create_cluster_prim(
    const DistanceConverter &dconv, init_matrix_t &im, int)
{
  return std::make_unique<ClusterPrim>(dconv, im);
}
