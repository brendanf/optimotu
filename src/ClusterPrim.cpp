// SPDX-FileCopyrightText: 2026 Brendan Furneaux <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#include "ClusterPrim.h"
#include "optimotu.h"

#ifdef OPTIMOTU_R
#include <RcppThread.h>
#endif

#include <algorithm>

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

void ClusterPrim::prepare_merge_slots(std::size_t n_slots)
{
  tourney_.prepare(n_slots, n);
}

void ClusterPrim::set_merge_slot(std::size_t slot)
{
  merge_slot_ = slot;
}

void ClusterPrim::receive_mst(std::size_t slot, std::vector<MstEdge> edges)
{
  tourney_.submit(slot, std::move(edges));
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
  if (tourney_.prepared() && parent == nullptr)
  {
    tourney_.wait();
    mst_ = tourney_.mst();
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
  mst_write_to_matrix(mst_, n, m, out);
}

void ClusterPrim::write_threshold_row(d_t t, int *dest) const
{
  std::shared_lock<std::shared_timed_mutex> lock(this->mutex);
  mst_write_threshold_row(mst_, n, t, dest);
}

#ifdef OPTIMOTU_R
Rcpp::List ClusterPrim::as_hclust(const Rcpp::CharacterVector &seqnames) const
{
  std::shared_lock<std::shared_timed_mutex> lock(this->mutex);
  return mst_as_hclust(mst_, n, dconv, seqnames);
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
