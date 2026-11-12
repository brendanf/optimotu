// SPDX-FileCopyrightText: 2026 Brendan Furneaux <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_PRIM_TILE_H_INCLUDED
#define OPTIMOTU_PRIM_TILE_H_INCLUDED

#include "ClusterPrim.h"
#include "MappedClusterAlgorithm.h"
#include "BipartitePairGenerator.h"

// If algo is ClusterPrim (or a mapped wrapper), grow the tile MST with
// dist_global(global_i, global_j, cutoff) and return true. MCA is not
// ClusterPrim (and multi-subset Prim is rejected at construction).
template <typename DistGlobal>
bool try_grow_prim(
    ClusterAlgorithm *algo,
    PairGenerator &pg,
    DistGlobal dist_global)
{
  if (!algo)
  {
    return false;
  }
  ClusterPrim *prim = dynamic_cast<ClusterPrim *>(algo);
  const std::vector<std::size_t> *fwd = nullptr;
  if (!prim)
  {
    if (auto *mapped = dynamic_cast<MappedClusterAlgorithm *>(algo))
    {
      prim = dynamic_cast<ClusterPrim *>(&mapped->inner_algorithm());
      fwd = mapped->fwd_map_ptr();
    }
  }
  if (!prim)
  {
    return false;
  }
  auto dist = [&](j_t i, j_t j, double cutoff)
  {
    if (fwd)
    {
      return dist_global(
          static_cast<j_t>((*fwd)[i]),
          static_cast<j_t>((*fwd)[j]),
          cutoff);
    }
    return dist_global(i, j, cutoff);
  };
  if (auto *bp = dynamic_cast<BipartitePairGenerator *>(&pg))
  {
    const j_t n_j = static_cast<j_t>(bp->seq_end_j() - bp->seq_begin_j());
    const j_t n_i = static_cast<j_t>(bp->seq_end_i() - bp->seq_begin_i());
    prim->grow_bipartite(n_j, n_i, dist);
  }
  else
  {
    prim->grow_complete(dist);
  }
  return true;
}

// Like try_grow_prim, but score(local_u, dist_by_local) fills distances
// from local vertex u. Local indices stay inside Prim; the scorer maps
// to packed rows when the tile is a MappedClusterAlgorithm.
template <typename Score>
bool try_grow_prim_batch(
    ClusterAlgorithm *algo,
    PairGenerator &pg,
    Score score)
{
  if (!algo)
  {
    return false;
  }
  ClusterPrim *prim = dynamic_cast<ClusterPrim *>(algo);
  if (!prim)
  {
    if (auto *mapped = dynamic_cast<MappedClusterAlgorithm *>(algo))
    {
      prim = dynamic_cast<ClusterPrim *>(&mapped->inner_algorithm());
    }
  }
  if (!prim)
  {
    return false;
  }
  if (auto *bp = dynamic_cast<BipartitePairGenerator *>(&pg))
  {
    const j_t n_j = static_cast<j_t>(bp->seq_end_j() - bp->seq_begin_j());
    const j_t n_i = static_cast<j_t>(bp->seq_end_i() - bp->seq_begin_i());
    prim->grow_bipartite_batch(n_j, n_i, score);
  }
  else
  {
    prim->grow_complete_batch(score);
  }
  return true;
}

#endif
