// SPDX-FileCopyrightText: 2026 Brendan Furneaux <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_TILE_PAIRS_H_INCLUDED
#define OPTIMOTU_TILE_PAIRS_H_INCLUDED

#include "AllPairGenerator.h"
#include "BipartitePairGenerator.h"
#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

// Walk every pair in a tile with nested index ranges for dense
// AllPairGenerator / BipartitePairGenerator tiles. Unknown PairGenerator
// types fall back to operator++ (e.g. sparse KmerPairGenerator).
//
// Callback signature: void(i0, j0, i_local, j_local) where i0/j0 are
// global sequence indices and i_local/j_local use the same tile-local
// numbering as HammingSplitClusterWorker (diagonal: i0 - offset;
// bipartite: i0 - begin_i + nj).
template <typename Fn>
void for_each_tile_pair(PairGenerator &pg, Fn &&fn) {
  if (auto *diag = dynamic_cast<AllPairGenerator *>(&pg)) {
    const std::size_t offset = diag->seq_offset();
    const std::size_t n_tile = diag->seq_count();
    for (std::size_t i0 = offset + 1; i0 < offset + n_tile; ++i0) {
      const std::size_t li = i0 - offset;
      for (std::size_t j0 = offset; j0 < i0; ++j0) {
        fn(i0, j0, li, j0 - offset);
      }
    }
    return;
  }
  if (auto *bip = dynamic_cast<BipartitePairGenerator *>(&pg)) {
    const std::size_t begin_i = bip->seq_begin_i();
    const std::size_t end_i = bip->seq_end_i();
    const std::size_t begin_j = bip->seq_begin_j();
    const std::size_t end_j = bip->seq_end_j();
    const std::size_t nj = end_j - begin_j;
    for (std::size_t i0 = begin_i; i0 < end_i; ++i0) {
      const std::size_t li = i0 - begin_i + nj;
      for (std::size_t j0 = begin_j; j0 < end_j; ++j0) {
        fn(i0, j0, li, j0 - begin_j);
      }
    }
    return;
  }
  // Sparse / unexpected generators: per-pair virtual increment.
  while (pg) {
    fn(pg.i0(), pg.j0(), pg.i(), pg.j());
    ++pg;
  }
}

// Round-robin over tiles assigned to this worker: begin, begin+n_workers, ...
template <typename Fn>
void for_each_assigned_tile(
  const std::vector<std::unique_ptr<PairGenerator>> &generators,
  std::size_t begin,
  std::size_t n_workers,
  Fn &&fn
) {
  if (n_workers == 0) {
    return;
  }
  for (std::size_t pg_index = begin; pg_index < generators.size();
       pg_index += n_workers) {
    fn(pg_index, *generators[pg_index]);
  }
}

#endif // OPTIMOTU_TILE_PAIRS_H_INCLUDED
