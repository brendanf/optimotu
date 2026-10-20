// SPDX-FileCopyrightText: 2026 Brendan Furneaux
// SPDX-License-Identifier: MIT

#include "Ksw2ClusterWorker.h"
#include "Ksw2Aligner.h"
#include "MappedClusterAlgorithm.h"
#include "tile_pairs.h"
#include "pairwise_alignment.h"
#include "wfa_identity_bound.h"

#include <algorithm>

Ksw2ClusterWorker::Ksw2ClusterWorker(
  const SequenceSet &seq,
  ClusterAlgorithm &clust_algo,
  DivisiblePairGenerator::Builder & pgb,
  const int match, const int mismatch,
  const int gap_open, const int gap_extend,
  const int gap_open2, const int gap_extend2,
  int verbose,
  std::size_t worker_threads
) : DistClusterWorker(seq, clust_algo, pgb, verbose, worker_threads),
match(match), mismatch(mismatch), gap_open(gap_open), gap_extend(gap_extend),
gap_open2(gap_open2), gap_extend2(gap_extend2) {};

template<int verbose>
void Ksw2SplitClusterWorker<verbose>::operator()(std::size_t begin, std::size_t end) {
  (void)end;
  Ksw2Aligner aligner{match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2};
  for_each_assigned_tile(
    pair_generators, begin, end,
    [&](std::size_t pg_index, PairGenerator &pg) {
      ClusterAlgorithm *my_algo = tile_algo(&pg);
      if (!my_algo) {
        return;
      }
      OPTIMOTU_DEBUG(2, << "Ksw2SplitClusterWorker thread " << pg_index << " entered" << std::endl);
      size_t my_prealigned = 0;
      size_t my_aligned = 0;
      const bool use_local =
          dynamic_cast<MappedClusterAlgorithm *>(my_algo) != nullptr;

      for_each_tile_pair(pg, [&](
        std::size_t i0, std::size_t j0, std::size_t li, std::size_t lj
      ) {
        double threshold = use_local
            ? my_algo->max_relevant_local(
                  static_cast<j_t>(li), static_cast<j_t>(lj))
            : my_algo->max_relevant(
                  static_cast<j_t>(i0), static_cast<j_t>(j0));
        OPTIMOTU_DEBUG(
          4,
          << "thread" << pg_index
          << ": seqs " << lj << " (j0=" << j0 << ")"
          << " and " << li << " (i0=" << i0 << ")"
          << " max relevant=" << threshold
          << std::endl
        );

        bool is_seqj_longer = seq[j0].size() > seq[i0].size();
        size_t s1 = is_seqj_longer ? i0 : j0;
        size_t s2 = is_seqj_longer ? j0 : i0;
        double l1 = seq[s1].size(), l2 = seq[s2].size();
        OPTIMOTU_DEBUG(
          4,
          << "#### seq " << i0 << " (l1=" << l1 << ") and "
          << j0 << " (l2=" << l2 <<")####" << std::endl
        );

        double sim_threshold = 1.0 - threshold;
        if (l1/l2 >= sim_threshold) {
          ++my_prealigned;
          WfaIdentityBound bound = wfa_identity_bound(
              l1, l2, threshold, match, mismatch,
              gap_open, gap_extend, gap_open2, gap_extend2
          );
          int w = std::max(bound.max_k, -bound.min_k);
          aligner.setBandWidth(w);
          // Longer as query, shorter as target (KSW2 convention in plan)
          double d = distance_ksw2(seq[s2], seq[s1], aligner);
          if (d < 1.0) ++my_aligned;

          OPTIMOTU_DEBUG(
            4,
            << " distance=" << d
            << std::endl
          );
          if (d < threshold) {
            if (use_local) {
              my_algo->apply_local(
                  static_cast<j_t>(li), static_cast<j_t>(lj), d);
            } else {
              (*my_algo)(static_cast<j_t>(i0), static_cast<j_t>(j0), d);
            }
          }
        }
        RcppThread::checkUserInterrupt();
      });
      mutex.lock();
      OPTIMOTU_DEBUG(2, << "thread " << pg_index << " ready to merge" << std::endl);
      _aligned += my_aligned;
      _prealigned += my_prealigned;
      mutex.unlock();
      finish_tile(my_algo);
      OPTIMOTU_DEBUG(2, << "thread " << pg_index << " done" << std::endl);
    }
  );
}

template <int verbose>
void Ksw2ConcurrentClusterWorker<verbose>::operator()(std::size_t begin, std::size_t end) {
  (void)end;
  size_t my_prealigned = 0;
  size_t my_aligned = 0;

  Ksw2Aligner aligner{match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2};
  for_each_assigned_tile(
    pair_generators, begin, end,
    [&](std::size_t pg_index, PairGenerator &pg) {
      OPTIMOTU_DEBUG(
        2,
        << "Ksw2ConcurrentClusterWorker thread " << pg_index
        << " entered" << std::endl
      );
      for_each_tile_pair(pg, [&](
        std::size_t i0, std::size_t j0, std::size_t, std::size_t
      ) {
        double threshold = clust_algo.max_relevant(
            static_cast<j_t>(i0), static_cast<j_t>(j0));
        OPTIMOTU_DEBUG(
          4,
          << "thread" << pg_index
          << ": seqs " << j0 << " (j0=" << j0 << ")"
          << " and " << i0 << " (i0=" << i0 << ")"
          << " max relevant=" << threshold
          << std::endl
        );
        bool is_seqj_longer = seq[j0].size() > seq[i0].size();
        size_t s1 = is_seqj_longer ? i0 : j0;
        size_t s2 = is_seqj_longer ? j0 : i0;
        double l1 = seq[s1].size(), l2 = seq[s2].size();
        OPTIMOTU_DEBUG(
          4,
          << "#### seq " << s1 << " (l1=" << l1 << ") and "
          << s2 << " (l2=" << l2 <<")####" << std::endl
        );

        double sim_threshold = 1.0 - threshold;
        if (l1/l2 >= sim_threshold) {
          ++my_prealigned;
          WfaIdentityBound bound = wfa_identity_bound(
              l1, l2, threshold, match, mismatch,
              gap_open, gap_extend, gap_open2, gap_extend2
          );
          int w = std::max(bound.max_k, -bound.min_k);
          aligner.setBandWidth(w);
          OPTIMOTU_DEBUG(
            4,
            << "ksw2 band width=" << w
            << " min_k=" << bound.min_k
            << " max_k=" << bound.max_k
            << std::endl
          );
          double d = distance_ksw2(seq[s2], seq[s1], aligner);
          if (d < 1.0) ++my_aligned;
          OPTIMOTU_DEBUG(
            4,
            << "Thread " << begin
            << " distance=" << d
            << std::endl
          );
          if (d < threshold) {
            clust_algo(static_cast<j_t>(i0), static_cast<j_t>(j0), d);
          }
        }
        OPTIMOTU_DEBUG(
          4,
          << "Thread " << pg_index
          << ": finished " << j0 << " (j0=" << j0 << ")"
          << " and " << i0 << " (i0=" << i0 << ")"
          << "\n" << std::endl
        );
        RcppThread::checkUserInterrupt();
      });
      OPTIMOTU_DEBUG(2, << "thread " << pg_index << " done" << std::endl);
    }
  );
  mutex.lock();
  _aligned += my_aligned;
  _prealigned += my_prealigned;
  mutex.unlock();
}

template class Ksw2SplitClusterWorker<0>;
template class Ksw2SplitClusterWorker<1>;
template class Ksw2SplitClusterWorker<2>;
template class Ksw2SplitClusterWorker<3>;
template class Ksw2SplitClusterWorker<4>;
template class Ksw2ConcurrentClusterWorker<0>;
template class Ksw2ConcurrentClusterWorker<1>;
template class Ksw2ConcurrentClusterWorker<2>;
template class Ksw2ConcurrentClusterWorker<3>;
template class Ksw2ConcurrentClusterWorker<4>;
