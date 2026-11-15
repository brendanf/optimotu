// SPDX-FileCopyrightText: 2026 Brendan Furneaux
// <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#include "optimotu.h"
#include "AaClusterWorker.h"
#include "MappedClusterAlgorithm.h"
#include "MultipleClusterAlgorithm.h"
#include "prim_tile.h"
#include "tile_pairs.h"
#include "AllPairGenerator.h"
#include "BipartitePairGenerator.h"
#ifdef OPTIMOTU_CUDA
#include "evo_gpu.h"
#endif

#include <algorithm>
#include <utility>

namespace {

constexpr std::size_t kAaInterruptBatch = 65536;

inline void aa_maybe_interrupt(std::size_t &pairs_since_check) {
  if (++pairs_since_check >= kAaInterruptBatch) {
    pairs_since_check = 0;
#ifdef OPTIMOTU_R
    RcppThread::checkUserInterrupt();
#endif
  }
}

// Position-major block of Prim vertices still outside the tree.
// n_cols is the original stride; live shrinks by swap-remove.
struct LiveAaBlock {
  int n_cols = 0;
  int live = 0;
  int alen = 0;
  std::vector<uint8_t> pos_major;
  std::vector<int> local_idx;
  std::vector<int> packed_idx;
  std::vector<int> lane_of;

  void init(
    const PackedAaSequenceSet &pss,
    const std::vector<int> &locals,
    const std::vector<int> &packed,
    int n_vert
  ) {
    const int n = static_cast<int>(locals.size());
    n_cols = n;
    live = n;
    alen = pss.alen;
    pss.fill_pos_major_idx(packed.data(), n, pos_major);
    local_idx = locals;
    packed_idx = packed;
    lane_of.assign(static_cast<std::size_t>(n_vert), -1);
    for (int t = 0; t < n; ++t) {
      lane_of[static_cast<std::size_t>(locals[static_cast<std::size_t>(t)])] =
        t;
    }
  }

  void remove(int u) {
    if (u < 0 || static_cast<std::size_t>(u) >= lane_of.size()) {
      return;
    }
    const int lane = lane_of[static_cast<std::size_t>(u)];
    if (lane < 0 || live <= 0) {
      return;
    }
    const int last = live - 1;
    if (lane != last) {
      for (int k = 0; k < alen; ++k) {
        const std::size_t base =
          static_cast<std::size_t>(k) * static_cast<std::size_t>(n_cols);
        std::swap(
          pos_major[base + static_cast<std::size_t>(lane)],
          pos_major[base + static_cast<std::size_t>(last)]
        );
      }
      std::swap(
        local_idx[static_cast<std::size_t>(lane)],
        local_idx[static_cast<std::size_t>(last)]
      );
      std::swap(
        packed_idx[static_cast<std::size_t>(lane)],
        packed_idx[static_cast<std::size_t>(last)]
      );
      lane_of[static_cast<std::size_t>(
        local_idx[static_cast<std::size_t>(lane)]
      )] = lane;
    }
    lane_of[static_cast<std::size_t>(u)] = -1;
    --live;
  }
};

}  // namespace

AaClusterWorker::AaClusterWorker(
  const SequenceSet &seq,
  ClusterAlgorithm &clust_algo,
  DivisiblePairGenerator::Builder &pgb,
  const int min_overlap,
  const bool ignore_gaps,
  const float *weight_lut_in,
  int verbose,
  std::size_t worker_threads,
  const EvoModel *evo_model
) : DistClusterWorker(seq, clust_algo, pgb, verbose, worker_threads),
    pss(seq),
    min_overlap(min_overlap),
    ignore_gaps(ignore_gaps),
    has_evo(evo_model != nullptr),
    evo(evo_model != nullptr ? *evo_model : EvoModel()) {
  if (weight_lut_in == nullptr) {
    weight_lut = nullptr;
  } else {
    lut_storage.assign(
      weight_lut_in,
      weight_lut_in + PackedAaSequenceSet::kLut * PackedAaSequenceSet::kLut
    );
    weight_lut = lut_storage.data();
  }
  tracked_pss_bytes = PackedAaSequenceSet::estimate_bytes(
    seq.size(), pss.alen
  );
  if (tracked_pss_bytes > 0) {
    if (auto *mca = dynamic_cast<MultipleClusterAlgorithm *>(&clust_algo)) {
      mca->acquire_memory(tracked_pss_bytes, "PackedAaSequenceSet");
    }
  }
}

void AaClusterWorker::score_targets(
  int query,
  const uint8_t *pos_major,
  const int *target_idx,
  int n_target,
  int n_cols,
  double *dist_out,
  int *n_ok_out
) const {
  if (has_evo) {
    pss.score_query_evo(
      query, pos_major, target_idx, n_target, n_cols,
      min_overlap, ignore_gaps, evo, dist_out, n_ok_out
    );
    return;
  }
  pss.score_query_vs_pos_major(
    query, pos_major, target_idx, n_target, n_cols,
    min_overlap, ignore_gaps, weight_lut, 1.0, false,
    dist_out, n_ok_out
  );
}

double AaClusterWorker::pair_distance(int i, int j) const {
  if (has_evo) {
    return pss.dist_evo(i, j, min_overlap, ignore_gaps, evo);
  }
  return pss.dist(i, j, min_overlap, ignore_gaps, weight_lut);
}

bool AaClusterWorker::below_unit(double d) const {
  if (has_evo) {
    return d < evo.t_max;
  }
  return d < 1.0;
}

AaClusterWorker::~AaClusterWorker() {
  if (tracked_pss_bytes > 0) {
    if (auto *mca = dynamic_cast<MultipleClusterAlgorithm *>(&clust_algo)) {
      mca->release_memory(tracked_pss_bytes);
    }
  }
}

template<int verbose>
AaSplitClusterWorker<verbose>::AaSplitClusterWorker(
  const SequenceSet &seq,
  ClusterAlgorithm &clust_algo,
  DivisiblePairGenerator::Builder &pgb,
  const int min_overlap,
  const bool ignore_gaps,
  const float *weight_lut,
  std::size_t worker_threads,
  const EvoModel *evo_model
) : AaClusterWorker(
      seq, clust_algo, pgb, min_overlap, ignore_gaps, weight_lut, verbose,
      worker_threads, evo_model
    ) {}

template<int verbose>
void AaSplitClusterWorker<verbose>::operator()(std::size_t begin, std::size_t end) {
  (void)end;
  size_t my_prealigned = 0;
  size_t my_aligned = 0;
  std::size_t pairs_since_check = 0;

  for_each_assigned_tile(
    pair_generators, begin, end,
    [&](std::size_t pg_index, PairGenerator &pg) {
      ClusterAlgorithm *my_algo = tile_algo(&pg);
      if (!my_algo) {
        finish_tile(nullptr, pg_index);
        return;
      }
      my_algo->set_merge_slot(pg_index);
      if (my_algo->uses_prim()) {
        const std::vector<std::size_t> *fwd = nullptr;
        if (auto *mapped = dynamic_cast<MappedClusterAlgorithm *>(my_algo)) {
          fwd = mapped->fwd_map_ptr();
        }
        auto packed_of = [&](j_t local) -> int {
          if (fwd) {
            return static_cast<int>((*fwd)[static_cast<std::size_t>(local)]);
          }
          return static_cast<int>(local);
        };
        auto *bp = dynamic_cast<BipartitePairGenerator *>(&pg);
        const j_t n_j = bp
          ? static_cast<j_t>(bp->seq_end_j() - bp->seq_begin_j())
          : 0;
        const j_t n_i = bp
          ? static_cast<j_t>(bp->seq_end_i() - bp->seq_begin_i())
          : 0;
        const j_t n_vert = bp ? (n_j + n_i) : (
          fwd
            ? static_cast<j_t>(fwd->size())
            : static_cast<j_t>(pss.num_seqs)
        );
        LiveAaBlock left;
        LiveAaBlock right;
        LiveAaBlock all;
        if (bp) {
          std::vector<int> locals_j(static_cast<std::size_t>(n_j));
          std::vector<int> packed_j(static_cast<std::size_t>(n_j));
          for (j_t t = 0; t < n_j; ++t) {
            locals_j[static_cast<std::size_t>(t)] = static_cast<int>(t);
            packed_j[static_cast<std::size_t>(t)] = packed_of(t);
          }
          std::vector<int> locals_i(static_cast<std::size_t>(n_i));
          std::vector<int> packed_i(static_cast<std::size_t>(n_i));
          for (j_t t = 0; t < n_i; ++t) {
            const j_t local = n_j + t;
            locals_i[static_cast<std::size_t>(t)] = static_cast<int>(local);
            packed_i[static_cast<std::size_t>(t)] = packed_of(local);
          }
          left.init(pss, locals_j, packed_j, static_cast<int>(n_vert));
          right.init(pss, locals_i, packed_i, static_cast<int>(n_vert));
        } else {
          std::vector<int> locals(static_cast<std::size_t>(n_vert));
          std::vector<int> packed(static_cast<std::size_t>(n_vert));
          for (j_t t = 0; t < n_vert; ++t) {
            locals[static_cast<std::size_t>(t)] = static_cast<int>(t);
            packed[static_cast<std::size_t>(t)] = packed_of(t);
          }
          all.init(pss, locals, packed, static_cast<int>(n_vert));
        }
        std::vector<double> lane_dist;
        std::vector<int> lane_ok;
        auto score_block = [&](LiveAaBlock &block, int packed_u) {
          if (block.live <= 0) {
            return;
          }
          lane_dist.assign(static_cast<std::size_t>(block.live), 1.0);
          lane_ok.assign(static_cast<std::size_t>(block.live), 0);
          score_targets(
            packed_u,
            block.pos_major.data(),
            block.packed_idx.data(),
            block.live,
            block.n_cols,
            lane_dist.data(),
            lane_ok.data()
          );
          for (int t = 0; t < block.live; ++t) {
            ++my_prealigned;
            if (below_unit(lane_dist[static_cast<std::size_t>(t)])) {
              ++my_aligned;
            }
            aa_maybe_interrupt(pairs_since_check);
          }
        };
        const bool grew = try_grow_prim_batch(
          my_algo, pg,
          [&](j_t u, double *dist_by_local) {
            const int packed_u = packed_of(u);
            if (bp) {
              if (u < n_j) {
                left.remove(static_cast<int>(u));
                score_block(right, packed_u);
                for (int t = 0; t < right.live; ++t) {
                  dist_by_local[right.local_idx[static_cast<std::size_t>(t)]] =
                    lane_dist[static_cast<std::size_t>(t)];
                }
              } else {
                right.remove(static_cast<int>(u));
                score_block(left, packed_u);
                for (int t = 0; t < left.live; ++t) {
                  dist_by_local[left.local_idx[static_cast<std::size_t>(t)]] =
                    lane_dist[static_cast<std::size_t>(t)];
                }
              }
            } else {
              all.remove(static_cast<int>(u));
              score_block(all, packed_u);
              for (int t = 0; t < all.live; ++t) {
                dist_by_local[all.local_idx[static_cast<std::size_t>(t)]] =
                  lane_dist[static_cast<std::size_t>(t)];
              }
            }
          }
        );
        if (!grew) {
          OPTIMOTU_STOP(
            "clust_prim() amino-acid batch grow requires ClusterPrim"
          );
        }
        mutex.lock();
        _aligned += my_aligned;
        _prealigned += my_prealigned;
        mutex.unlock();
        my_aligned = 0;
        my_prealigned = 0;
        finish_tile(my_algo, pg_index);
        return;
      }
      const bool use_local =
        dynamic_cast<MappedClusterAlgorithm *>(my_algo) != nullptr;

      auto apply_pair = [&](std::size_t i0, std::size_t j0, std::size_t li,
                            std::size_t lj, double d) {
        ++my_prealigned;
        if (below_unit(d)) {
          ++my_aligned;
        }
        if (use_local) {
          double threshold = my_algo->max_relevant_local(
            static_cast<j_t>(li), static_cast<j_t>(lj)
          );
          if (d <= threshold) {
            my_algo->apply_local(
              static_cast<j_t>(li), static_cast<j_t>(lj), d
            );
          }
        } else {
          double threshold = my_algo->max_relevant(
            static_cast<j_t>(i0), static_cast<j_t>(j0)
          );
          if (d <= threshold) {
            (*my_algo)(static_cast<j_t>(i0), static_cast<j_t>(j0), d);
          }
        }
        aa_maybe_interrupt(pairs_since_check);
      };

      if (auto *diag = dynamic_cast<AllPairGenerator *>(&pg)) {
        const int offset = static_cast<int>(diag->seq_offset());
        const int n_tile = static_cast<int>(diag->seq_count());
        std::vector<uint8_t> pos_major;
        pss.fill_pos_major(offset, offset + n_tile, pos_major);
        std::vector<int> target_idx(static_cast<std::size_t>(n_tile));
        for (int t = 0; t < n_tile; ++t) {
          target_idx[static_cast<std::size_t>(t)] = offset + t;
        }
        std::vector<double> dists(static_cast<std::size_t>(n_tile));
        std::vector<int> n_ok(static_cast<std::size_t>(n_tile));
        for (int i = offset + 1; i < offset + n_tile; ++i) {
          const int n_j = i - offset;
          score_targets(
            i, pos_major.data(), target_idx.data(), n_j, n_tile,
            dists.data(), n_ok.data()
          );
          const std::size_t li = static_cast<std::size_t>(i - offset);
          for (int t = 0; t < n_j; ++t) {
            apply_pair(
              static_cast<std::size_t>(i),
              static_cast<std::size_t>(offset + t),
              li,
              static_cast<std::size_t>(t),
              dists[static_cast<std::size_t>(t)]
            );
          }
        }
      } else if (auto *bip = dynamic_cast<BipartitePairGenerator *>(&pg)) {
        const int begin_i = static_cast<int>(bip->seq_begin_i());
        const int end_i = static_cast<int>(bip->seq_end_i());
        const int begin_j = static_cast<int>(bip->seq_begin_j());
        const int end_j = static_cast<int>(bip->seq_end_j());
        const int n_j = end_j - begin_j;
        std::vector<uint8_t> pos_major;
        pss.fill_pos_major(begin_j, end_j, pos_major);
        std::vector<int> target_idx(static_cast<std::size_t>(n_j));
        for (int t = 0; t < n_j; ++t) {
          target_idx[static_cast<std::size_t>(t)] = begin_j + t;
        }
        std::vector<double> dists(static_cast<std::size_t>(n_j));
        std::vector<int> n_ok(static_cast<std::size_t>(n_j));
        for (int i = begin_i; i < end_i; ++i) {
          score_targets(
            i, pos_major.data(), target_idx.data(), n_j, n_j,
            dists.data(), n_ok.data()
          );
          const std::size_t li =
            static_cast<std::size_t>(i - begin_i + n_j);
          for (int t = 0; t < n_j; ++t) {
            apply_pair(
              static_cast<std::size_t>(i),
              static_cast<std::size_t>(begin_j + t),
              li,
              static_cast<std::size_t>(t),
              dists[static_cast<std::size_t>(t)]
            );
          }
        }
      } else {
        for_each_tile_pair(pg, [&](
          std::size_t i0, std::size_t j0, std::size_t li, std::size_t lj
        ) {
          double d = pair_distance(static_cast<int>(i0), static_cast<int>(j0));
          apply_pair(i0, j0, li, lj, d);
        });
      }

#ifdef OPTIMOTU_R
      RcppThread::checkUserInterrupt();
#endif
      pairs_since_check = 0;
      mutex.lock();
      _aligned += my_aligned;
      _prealigned += my_prealigned;
      mutex.unlock();
      my_aligned = 0;
      my_prealigned = 0;
      finish_tile(my_algo, pg_index);
    });
}

template<int verbose>
AaConcurrentClusterWorker<verbose>::AaConcurrentClusterWorker(
  const SequenceSet &seq,
  ClusterAlgorithm &clust_algo,
  DivisiblePairGenerator::Builder &pgb,
  const int min_overlap,
  const bool ignore_gaps,
  const float *weight_lut,
  std::size_t worker_threads,
  const EvoModel *evo_model
) : AaClusterWorker(
      seq, clust_algo, pgb, min_overlap, ignore_gaps, weight_lut, verbose,
      worker_threads, evo_model
    ) {}

template<int verbose>
void AaConcurrentClusterWorker<verbose>::operator()(
  std::size_t begin,
  std::size_t end
) {
  (void)end;
  size_t my_prealigned = 0;
  size_t my_aligned = 0;
  std::size_t pairs_since_check = 0;

  for_each_assigned_tile(
    pair_generators, begin, end,
    [&](std::size_t, PairGenerator &pg) {
      for_each_tile_pair(pg, [&](
        std::size_t i0, std::size_t j0, std::size_t, std::size_t
      ) {
        double threshold = clust_algo.max_relevant(
          static_cast<j_t>(i0), static_cast<j_t>(j0)
        );
        ++my_prealigned;
        double d = pair_distance(static_cast<int>(i0), static_cast<int>(j0));
        if (below_unit(d)) {
          ++my_aligned;
        }
        if (d < threshold) {
          clust_algo(static_cast<j_t>(i0), static_cast<j_t>(j0), d);
        }
        aa_maybe_interrupt(pairs_since_check);
      });
#ifdef OPTIMOTU_R
      RcppThread::checkUserInterrupt();
#endif
      pairs_since_check = 0;
      mutex.lock();
      _aligned += my_aligned;
      _prealigned += my_prealigned;
      mutex.unlock();
      my_aligned = 0;
      my_prealigned = 0;
    }
  );
}

bool AaClusterWorker::wants_gpu() const {
  return has_evo;
}

std::unique_ptr<GpuDistanceKernel> AaClusterWorker::make_gpu_kernel() {
#ifdef OPTIMOTU_CUDA
  return make_evo_gpu_kernel(pss, evo);
#else
  return nullptr;
#endif
}

namespace {

double gpu_lookup(
  const GpuTileSpec &spec,
  const GpuTileResult &result,
  int i0,
  int j0
) {
  if (spec.kind == GpuTileSpec::Kind::Triangle) {
    int i = i0 - spec.offset;
    int j = j0 - spec.offset;
    if (i < j) {
      std::swap(i, j);
    }
    const std::size_t p =
      static_cast<std::size_t>(i) * static_cast<std::size_t>(i - 1) / 2 +
      static_cast<std::size_t>(j);
    return result.distances[p];
  }
  const int nj = spec.j1 - spec.j0;
  const int row = i0 - spec.i0;
  const int col = j0 - spec.j0;
  return result.distances[
    static_cast<std::size_t>(row) * static_cast<std::size_t>(nj) +
    static_cast<std::size_t>(col)
  ];
}

}  // namespace

template<int verbose>
void AaSplitClusterWorker<verbose>::consume_gpu_tile(
  std::size_t pg_index,
  PairGenerator &pg,
  const GpuTileResult &result
) {
  ClusterAlgorithm *my_algo = tile_algo(&pg);
  if (!my_algo) {
    finish_tile(nullptr, pg_index);
    return;
  }
  my_algo->set_merge_slot(pg_index);
  const GpuTileSpec spec = gpu_tile_from_generator(pg);
  size_t my_prealigned = 0;
  size_t my_aligned = 0;
  if (my_algo->uses_prim()) {
    const bool grew = try_grow_prim(
      my_algo, pg,
      [&](j_t gi, j_t gj, double) {
        const double d = gpu_lookup(
          spec, result, static_cast<int>(gi), static_cast<int>(gj)
        );
        ++my_prealigned;
        if (below_unit(d)) {
          ++my_aligned;
        }
        return d;
      }
    );
    if (!grew) {
      OPTIMOTU_STOP(
        "clust_prim() amino-acid batch grow requires ClusterPrim"
      );
    }
    mutex.lock();
    _aligned += my_aligned;
    _prealigned += my_prealigned;
    mutex.unlock();
    finish_tile(my_algo, pg_index);
    return;
  }
  const bool use_local =
    dynamic_cast<MappedClusterAlgorithm *>(my_algo) != nullptr;
  for_each_tile_pair(pg, [&](
    std::size_t i0, std::size_t j0, std::size_t li, std::size_t lj
  ) {
    const double d = gpu_lookup(
      spec, result, static_cast<int>(i0), static_cast<int>(j0)
    );
    ++my_prealigned;
    if (below_unit(d)) {
      ++my_aligned;
    }
    if (use_local) {
      const double threshold = my_algo->max_relevant_local(
        static_cast<j_t>(li), static_cast<j_t>(lj)
      );
      if (d <= threshold) {
        my_algo->apply_local(
          static_cast<j_t>(li), static_cast<j_t>(lj), d
        );
      }
    } else {
      const double threshold = my_algo->max_relevant(
        static_cast<j_t>(i0), static_cast<j_t>(j0)
      );
      if (d <= threshold) {
        (*my_algo)(static_cast<j_t>(i0), static_cast<j_t>(j0), d);
      }
    }
  });
  mutex.lock();
  _aligned += my_aligned;
  _prealigned += my_prealigned;
  mutex.unlock();
  finish_tile(my_algo, pg_index);
}

template<int verbose>
void AaConcurrentClusterWorker<verbose>::consume_gpu_tile(
  std::size_t,
  PairGenerator &pg,
  const GpuTileResult &result
) {
  const GpuTileSpec spec = gpu_tile_from_generator(pg);
  size_t my_prealigned = 0;
  size_t my_aligned = 0;
  for_each_tile_pair(pg, [&](
    std::size_t i0, std::size_t j0, std::size_t, std::size_t
  ) {
    const double d = gpu_lookup(
      spec, result, static_cast<int>(i0), static_cast<int>(j0)
    );
    const double threshold = clust_algo.max_relevant(
      static_cast<j_t>(i0), static_cast<j_t>(j0)
    );
    ++my_prealigned;
    if (below_unit(d)) {
      ++my_aligned;
    }
    if (d < threshold) {
      clust_algo(static_cast<j_t>(i0), static_cast<j_t>(j0), d);
    }
  });
  mutex.lock();
  _aligned += my_aligned;
  _prealigned += my_prealigned;
  mutex.unlock();
}

template class AaSplitClusterWorker<0>;
template class AaSplitClusterWorker<1>;
template class AaSplitClusterWorker<2>;
template class AaSplitClusterWorker<3>;
template class AaSplitClusterWorker<4>;
template class AaConcurrentClusterWorker<0>;
template class AaConcurrentClusterWorker<1>;
template class AaConcurrentClusterWorker<2>;
template class AaConcurrentClusterWorker<3>;
template class AaConcurrentClusterWorker<4>;
