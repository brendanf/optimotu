#include "HammingClusterWorker.h"
#include "MappedClusterAlgorithm.h"
#include "MultipleClusterAlgorithm.h"
#include "prim_tile.h"
#include "tile_pairs.h"
#include "AllPairGenerator.h"
#include "BipartitePairGenerator.h"
#include "hamming_kernel.h"

#include <algorithm>
#include <utility>

extern "C" {
#include "defs.h"
}

typedef RcppParallel::RMatrix<int> matrix_t;

namespace {

constexpr std::size_t kHammingInterruptBatch = 65536;

inline void hamming_maybe_interrupt(std::size_t &pairs_since_check) {
  if (++pairs_since_check >= kHammingInterruptBatch) {
    pairs_since_check = 0;
#ifdef OPTIMOTU_R
    RcppThread::checkUserInterrupt();
#endif
  }
}

struct LiveHammingBlock {
  int n_cols = 0;
  int live = 0;
  int ulen = 0;
  int mulen = 0;
  std::vector<uint64_t> seq_major;
  std::vector<uint64_t> mask_major;
  std::vector<int> local_idx;
  std::vector<int> packed_idx;
  std::vector<int> seq_start;
  std::vector<int> seq_end;
  std::vector<int> lane_of;

  void init(
    const PackedSequenceSet &pss,
    const std::vector<int> &locals,
    const std::vector<int> &packed,
    int n_vert
  ) {
    const int n = static_cast<int>(locals.size());
    n_cols = n;
    live = n;
    ulen = pss.ulen;
    mulen = pss.mulen;
    pss.fill_word_major_idx(packed.data(), n, seq_major, mask_major);
    local_idx = locals;
    packed_idx = packed;
    seq_start.resize(static_cast<std::size_t>(n));
    seq_end.resize(static_cast<std::size_t>(n));
    for (int t = 0; t < n; ++t) {
      const int pi = packed[static_cast<std::size_t>(t)];
      seq_start[static_cast<std::size_t>(t)] =
        pss.start[static_cast<std::size_t>(pi)];
      seq_end[static_cast<std::size_t>(t)] =
        pss.end[static_cast<std::size_t>(pi)];
    }
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
      for (int w = 0; w < ulen; ++w) {
        const std::size_t base =
          static_cast<std::size_t>(w) * static_cast<std::size_t>(n_cols);
        std::swap(
          seq_major[base + static_cast<std::size_t>(lane)],
          seq_major[base + static_cast<std::size_t>(last)]
        );
      }
      for (int w = 0; w < mulen; ++w) {
        const std::size_t base =
          static_cast<std::size_t>(w) * static_cast<std::size_t>(n_cols);
        std::swap(
          mask_major[base + static_cast<std::size_t>(lane)],
          mask_major[base + static_cast<std::size_t>(last)]
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
      std::swap(
        seq_start[static_cast<std::size_t>(lane)],
        seq_start[static_cast<std::size_t>(last)]
      );
      std::swap(
        seq_end[static_cast<std::size_t>(lane)],
        seq_end[static_cast<std::size_t>(last)]
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

HammingClusterWorker::HammingClusterWorker(
  const SequenceSet &seq,
  ClusterAlgorithm &clust_algo,
  DivisiblePairGenerator::Builder & pgb,
  const int min_overlap,
  const bool ignore_gaps,
  int verbose,
  std::size_t worker_threads
) : DistClusterWorker(seq, clust_algo, pgb, verbose, worker_threads), pss(seq),
min_overlap(min_overlap), ignore_gaps(ignore_gaps) {
  tracked_pss_bytes = PackedSequenceSet::estimate_bytes(
    seq.size(), pss.alen
  );
  if (tracked_pss_bytes > 0) {
    if (auto *mca = dynamic_cast<MultipleClusterAlgorithm *>(&clust_algo)) {
      mca->acquire_memory(tracked_pss_bytes, "PackedSequenceSet");
    }
  }
}

HammingClusterWorker::~HammingClusterWorker() {
  if (tracked_pss_bytes > 0) {
    if (auto *mca = dynamic_cast<MultipleClusterAlgorithm *>(&clust_algo)) {
      mca->release_memory(tracked_pss_bytes);
    }
  }
}

template<int verbose>
HammingSplitClusterWorker<verbose>::HammingSplitClusterWorker(
  const SequenceSet &seq,
  ClusterAlgorithm &clust_algo,
  DivisiblePairGenerator::Builder & pgb,
  const int min_overlap,
  const bool ignore_gaps,
  std::size_t worker_threads
) : HammingClusterWorker(seq, clust_algo, pgb, min_overlap, ignore_gaps,
   verbose, worker_threads) {};

template<int verbose>
void HammingSplitClusterWorker<verbose>::operator()(std::size_t begin, std::size_t end) {
  (void)end;
  size_t my_prealigned = 0;
  size_t my_aligned = 0;
  std::size_t pairs_since_check = 0;

  for_each_assigned_tile(
      pair_generators, begin, end,
      [&](std::size_t pg_index, PairGenerator &pg)
      {
        ClusterAlgorithm *my_algo = tile_algo(&pg);
        if (!my_algo)
        {
          finish_tile(nullptr, pg_index);
          return;
        }
        my_algo->set_merge_slot(pg_index);
        if (my_algo->uses_prim() &&
            (ignore_gaps || OPTIMOTU_HAMMING_PDISTB2_STRIPE)) {
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
          LiveHammingBlock left;
          LiveHammingBlock right;
          LiveHammingBlock all;
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
          auto score_block = [&](LiveHammingBlock &block, int packed_u) {
            if (block.live <= 0) {
              return;
            }
            lane_dist.assign(static_cast<std::size_t>(block.live), 1.0);
            lane_ok.assign(static_cast<std::size_t>(block.live), 0);
            pss.score_query_vs_word_major(
              packed_u,
              block.seq_major.data(),
              block.mask_major.data(),
              block.live,
              block.n_cols,
              min_overlap,
              lane_dist.data(),
              lane_ok.data(),
              nullptr,
              ignore_gaps ? nullptr : block.seq_start.data(),
              ignore_gaps ? nullptr : block.seq_end.data()
            );
            for (int t = 0; t < block.live; ++t) {
              ++my_prealigned;
              if (lane_dist[static_cast<std::size_t>(t)] < 1.0) {
                ++my_aligned;
              }
              hamming_maybe_interrupt(pairs_since_check);
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
              "clust_prim() Hamming batch grow requires ClusterPrim"
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
        if (try_grow_prim(
                my_algo, pg,
                [&](j_t gi, j_t gj, double cutoff)
                {
                  (void)cutoff;
                  ++my_prealigned;
                  double d = pss.dist(
                      static_cast<int>(gi),
                      static_cast<int>(gj),
                      min_overlap,
                      ignore_gaps);
                  if (d < 1.0)
                  {
                    ++my_aligned;
                  }
                  hamming_maybe_interrupt(pairs_since_check);
                  return d;
                }))
        {
          mutex.lock();
          _aligned += my_aligned;
          _prealigned += my_prealigned;
          mutex.unlock();
          my_aligned = 0;
          my_prealigned = 0;
          finish_tile(my_algo, pg_index);
          return;
        }
        OPTIMOTU_DEBUG(
            2,
            << "HammingSplitClusterWorker thread " << pg_index
            << " entered" << std::endl);

        const bool use_local =
            dynamic_cast<MappedClusterAlgorithm *>(my_algo) != nullptr;

        auto apply_pair = [&](std::size_t i0, std::size_t j0, std::size_t li,
                              std::size_t lj, double d) {
          ++my_prealigned;
          if (d < 1.0) {
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
          hamming_maybe_interrupt(pairs_since_check);
        };

        if ((ignore_gaps || OPTIMOTU_HAMMING_PDISTB2_STRIPE) &&
            dynamic_cast<AllPairGenerator *>(&pg) != nullptr) {
          auto *diag = static_cast<AllPairGenerator *>(&pg);
          const int offset = static_cast<int>(diag->seq_offset());
          const int n_tile = static_cast<int>(diag->seq_count());
          std::vector<uint64_t> seq_major;
          std::vector<uint64_t> mask_major;
          pss.fill_word_major(offset, offset + n_tile, seq_major, mask_major);
          std::vector<double> dists(static_cast<std::size_t>(n_tile));
          std::vector<int> n_ok(static_cast<std::size_t>(n_tile));
          const int *tstart =
            ignore_gaps ? nullptr : pss.start.data() + offset;
          const int *tend = ignore_gaps ? nullptr : pss.end.data() + offset;
          for (int i = offset + 1; i < offset + n_tile; ++i) {
            const int n_j = i - offset;
            pss.score_query_vs_word_major(
              i, seq_major.data(), mask_major.data(), n_j, n_tile,
              min_overlap, dists.data(), n_ok.data(), nullptr, tstart, tend
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
        } else if ((ignore_gaps || OPTIMOTU_HAMMING_PDISTB2_STRIPE) &&
                   dynamic_cast<BipartitePairGenerator *>(&pg) != nullptr) {
          auto *bip = static_cast<BipartitePairGenerator *>(&pg);
          const int begin_i = static_cast<int>(bip->seq_begin_i());
          const int end_i = static_cast<int>(bip->seq_end_i());
          const int begin_j = static_cast<int>(bip->seq_begin_j());
          const int end_j = static_cast<int>(bip->seq_end_j());
          const int n_j = end_j - begin_j;
          std::vector<uint64_t> seq_major;
          std::vector<uint64_t> mask_major;
          pss.fill_word_major(begin_j, end_j, seq_major, mask_major);
          std::vector<double> dists(static_cast<std::size_t>(n_j));
          std::vector<int> n_ok(static_cast<std::size_t>(n_j));
          const int *tstart =
            ignore_gaps ? nullptr : pss.start.data() + begin_j;
          const int *tend = ignore_gaps ? nullptr : pss.end.data() + begin_j;
          for (int i = begin_i; i < end_i; ++i) {
            pss.score_query_vs_word_major(
              i, seq_major.data(), mask_major.data(), n_j, n_j,
              min_overlap, dists.data(), n_ok.data(), nullptr, tstart, tend
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
            double d = pss.dist(
              static_cast<int>(i0), static_cast<int>(j0), min_overlap,
              ignore_gaps
            );
            apply_pair(i0, j0, li, lj, d);
          });
        }

#ifdef OPTIMOTU_R
        RcppThread::checkUserInterrupt();
#endif
        pairs_since_check = 0;

        mutex.lock();
        OPTIMOTU_DEBUG(2, << "thread " << pg_index << " ready to merge" << std::endl);
        _aligned += my_aligned;
        _prealigned += my_prealigned;
        mutex.unlock();
        my_aligned = 0;
        my_prealigned = 0;
        finish_tile(my_algo, pg_index);
        OPTIMOTU_DEBUG(2, << "thread " << pg_index << " done" << std::endl);
      });
}

template<int verbose>
HammingConcurrentClusterWorker<verbose>::HammingConcurrentClusterWorker(
  const SequenceSet &seq,
  ClusterAlgorithm &clust_algo,
  DivisiblePairGenerator::Builder & pgb,
  const int min_overlap,
  const bool ignore_gaps,
  std::size_t worker_threads
) : HammingClusterWorker(seq, clust_algo, pgb, min_overlap, ignore_gaps,
   verbose, worker_threads) {};

template<int verbose>
void HammingConcurrentClusterWorker<verbose>::operator()(std::size_t begin, std::size_t end) {
  (void)end;
  size_t my_prealigned = 0;
  size_t my_aligned = 0;
  std::size_t pairs_since_check = 0;

  for_each_assigned_tile(
    pair_generators, begin, end,
    [&](std::size_t pg_index, PairGenerator &pg) {
      OPTIMOTU_DEBUG(
        2,
        << "HammingConcurrentClusterWorker thread " << pg_index
        << " entered" << std::endl
      );

      auto apply_one = [&](std::size_t i0, std::size_t j0, double d) {
        double threshold = clust_algo.max_relevant(
          static_cast<j_t>(i0), static_cast<j_t>(j0)
        );
        ++my_prealigned;
        if (d < 1.0) {
          ++my_aligned;
        }
        if (d < threshold) {
          clust_algo(static_cast<j_t>(i0), static_cast<j_t>(j0), d);
        }
        hamming_maybe_interrupt(pairs_since_check);
      };

      if ((ignore_gaps || OPTIMOTU_HAMMING_PDISTB2_STRIPE) &&
          dynamic_cast<AllPairGenerator *>(&pg) != nullptr) {
        auto *diag = static_cast<AllPairGenerator *>(&pg);
        const int offset = static_cast<int>(diag->seq_offset());
        const int n_tile = static_cast<int>(diag->seq_count());
        std::vector<uint64_t> seq_major;
        std::vector<uint64_t> mask_major;
        pss.fill_word_major(offset, offset + n_tile, seq_major, mask_major);
        std::vector<double> dists(static_cast<std::size_t>(n_tile));
        std::vector<int> n_ok(static_cast<std::size_t>(n_tile));
        const int *tstart =
          ignore_gaps ? nullptr : pss.start.data() + offset;
        const int *tend = ignore_gaps ? nullptr : pss.end.data() + offset;
        for (int i = offset + 1; i < offset + n_tile; ++i) {
          const int n_j = i - offset;
          pss.score_query_vs_word_major(
            i, seq_major.data(), mask_major.data(), n_j, n_tile,
            min_overlap, dists.data(), n_ok.data(), nullptr, tstart, tend
          );
          for (int t = 0; t < n_j; ++t) {
            apply_one(
              static_cast<std::size_t>(i),
              static_cast<std::size_t>(offset + t),
              dists[static_cast<std::size_t>(t)]
            );
          }
        }
      } else if ((ignore_gaps || OPTIMOTU_HAMMING_PDISTB2_STRIPE) &&
                 dynamic_cast<BipartitePairGenerator *>(&pg) != nullptr) {
        auto *bip = static_cast<BipartitePairGenerator *>(&pg);
        const int begin_i = static_cast<int>(bip->seq_begin_i());
        const int end_i = static_cast<int>(bip->seq_end_i());
        const int begin_j = static_cast<int>(bip->seq_begin_j());
        const int end_j = static_cast<int>(bip->seq_end_j());
        const int n_j = end_j - begin_j;
        std::vector<uint64_t> seq_major;
        std::vector<uint64_t> mask_major;
        pss.fill_word_major(begin_j, end_j, seq_major, mask_major);
        std::vector<double> dists(static_cast<std::size_t>(n_j));
        std::vector<int> n_ok(static_cast<std::size_t>(n_j));
        const int *tstart =
          ignore_gaps ? nullptr : pss.start.data() + begin_j;
        const int *tend = ignore_gaps ? nullptr : pss.end.data() + begin_j;
        for (int i = begin_i; i < end_i; ++i) {
          pss.score_query_vs_word_major(
            i, seq_major.data(), mask_major.data(), n_j, n_j,
            min_overlap, dists.data(), n_ok.data(), nullptr, tstart, tend
          );
          for (int t = 0; t < n_j; ++t) {
            apply_one(
              static_cast<std::size_t>(i),
              static_cast<std::size_t>(begin_j + t),
              dists[static_cast<std::size_t>(t)]
            );
          }
        }
      } else {
        for_each_tile_pair(pg, [&](
          std::size_t i0, std::size_t j0, std::size_t, std::size_t
        ) {
          apply_one(
            i0, j0,
            pss.dist(
              static_cast<int>(i0), static_cast<int>(j0), min_overlap,
              ignore_gaps
            )
          );
        });
      }

#ifdef OPTIMOTU_R
      RcppThread::checkUserInterrupt();
#endif
      pairs_since_check = 0;

      mutex.lock();
      _aligned += my_aligned;
      _prealigned += my_prealigned;
      OPTIMOTU_DEBUG(2, << "thread " << pg_index << " done" << std::endl);
      mutex.unlock();
      my_aligned = 0;
      my_prealigned = 0;
    }
  );
}

template class HammingSplitClusterWorker<0>;
template class HammingSplitClusterWorker<1>;
template class HammingSplitClusterWorker<2>;
template class HammingSplitClusterWorker<3>;
template class HammingSplitClusterWorker<4>;
template class HammingConcurrentClusterWorker<0>;
template class HammingConcurrentClusterWorker<1>;
template class HammingConcurrentClusterWorker<2>;
template class HammingConcurrentClusterWorker<3>;
template class HammingConcurrentClusterWorker<4>;
