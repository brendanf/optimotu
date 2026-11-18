#include "HybridClusterWorker.h"
#include "MappedClusterAlgorithm.h"
#include "prim_tile.h"
#include "tile_pairs.h"
#include "ClusterAlgorithmFactory.h"
#include "config.h"

#include <bindings/cpp/WFAligner.hpp>
#include <edlib.h>
#include "pairwise_alignment.h"
#include "wfa_identity_bound.h"
#include "score_dissimilarity.h"
#include "wfa_fatal.h"

typedef RcppParallel::RMatrix<int> matrix_t;

HybridClusterWorker::HybridClusterWorker(
  const SequenceSet &seq,
  ClusterAlgorithm &clust_algo,
  DivisiblePairGenerator::Builder & pgb,
  const double breakpoint,
  int verbose,
  std::size_t worker_threads,
  Dissimilarity dissimilarity
) : DistClusterWorker(seq, clust_algo, pgb, verbose, worker_threads, dissimilarity), breakpoint(breakpoint) {};

template<int verbose>
void HybridSplitClusterWorker<verbose>::operator()(std::size_t begin, std::size_t end) {
  (void)end;

  EdlibAlignConfig ed_aligner = edlibNewAlignConfig(
    -1, EdlibAlignMode::EDLIB_MODE_NW,
    dissimilarity == Dissimilarity::SCORE
      ? EdlibAlignTask::EDLIB_TASK_DISTANCE
      : EdlibAlignTask::EDLIB_TASK_PATH,
    0, 0
  );
  wfa::WFAlignerEdit wfa_aligner{
    dissimilarity == Dissimilarity::SCORE
      ? wfa::WFAligner::Score
      : wfa::WFAligner::Alignment
  };
  if (wfa_aligner.fatal())
  {
    optimotu_wfa_note_fatal();
    return;
  }
  for_each_assigned_tile(
    pair_generators, begin, end,
    [&](std::size_t pg_index, PairGenerator &pg) {
      OPTIMOTU_DEBUG(
        2,
        << "HybridSplitClusterWorker thread " << pg_index
        << " entered" << std::endl
      );
      size_t my_prealigned = 0;
      size_t my_aligned = 0;
      ClusterAlgorithm *my_algo = tile_algo(&pg);
      if (!my_algo) {
        finish_tile(nullptr, pg_index);
        return;
      }
      my_algo->set_merge_slot(pg_index);
      if (try_grow_prim(
              my_algo, pg,
              [&](j_t gi, j_t gj, double threshold)
              {
                if (optimotu_wfa_fatal_pending() || wfa_aligner.fatal())
                {
                  optimotu_wfa_note_fatal();
                  return 1.0;
                }
                bool is_seqj_longer = seq[gj].size() > seq[gi].size();
                size_t s1 = is_seqj_longer ? gi : gj;
                size_t s2 = is_seqj_longer ? gj : gi;
                double l1 = seq[s1].size(), l2 = seq[s2].size();
                double sim_threshold = 1.0 - threshold;
                if (l1 / l2 < sim_threshold)
                {
                  return 1.0;
                }
                ++my_prealigned;
                double maxd1 = wfa_identity_max_edits(l1, l2, threshold);
                bool is_close = breakpoint >= 1 ? maxd1 < breakpoint
                                                : threshold < breakpoint;
                double d;
                if (is_close)
                {
                  WfaIdentityBound bound = wfa_identity_bound(
                      l1, l2, threshold, 0, 1, 0, 1, 0, 1);
                  wfa_aligner.setHeuristicBandedStatic(bound.min_k, bound.max_k);
                  wfa_aligner.setMaxAlignmentSteps(bound.max_alignment_steps);
                  d = distance_wfa2<>(seq[s1], seq[s2], wfa_aligner);
                }
                else
                {
                  ed_aligner.k = (int)maxd1 + 1;
                  d = distance_edlib(seq[s1], seq[s2], ed_aligner);
                }
                if (wfa_aligner.fatal() || optimotu_wfa_fatal_pending())
                {
                  optimotu_wfa_note_fatal();
                  return 1.0;
                }
                if (d < 1.0)
                {
                  ++my_aligned;
                }
                return d;
              }))
      {
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
                                 std::size_t i0, std::size_t j0, std::size_t li, std::size_t lj)
                         {
        if (optimotu_wfa_fatal_pending()) return;
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
          << "#### seq " << s1 << " (l1=" << l1 << ") and "
          << s2 << " (l2=" << l2 <<")####" << std::endl
        );

        if (pair_length_ok(
              dissimilarity, AlignmentSpan::GLOBAL, l1, l2, threshold,
              0, 1, 0, 1, 0, 1
            )) {
          ++my_prealigned;
          double maxd1 = dissimilarity == Dissimilarity::SCORE
              ? edlib_score_k(score_ref_length(l1, l2, AlignmentSpan::GLOBAL), threshold)
              : wfa_identity_max_edits(l1, l2, threshold);
          bool is_close = breakpoint >= 1 ? maxd1 < breakpoint : threshold < breakpoint;
          double d;
          if (is_close) {
            WfaIdentityBound bound = choose_alignment_bound(
                dissimilarity, l1, l2, threshold, 0, 1, 0, 1, 0, 1
            );
            wfa_aligner.setHeuristicBandedStatic(bound.min_k, bound.max_k);
            wfa_aligner.setMaxAlignmentSteps(bound.max_alignment_steps);
            d = distance_wfa2<>(seq[s1], seq[s2], wfa_aligner, dissimilarity, 0, 1);
          } else {
            ed_aligner.k = dissimilarity == Dissimilarity::SCORE
                ? (int)maxd1
                : (int)maxd1 + 1;
            d = distance_edlib(seq[s1], seq[s2], ed_aligner, dissimilarity);
          }
          if (wfa_aligner.fatal() || optimotu_wfa_fatal_pending()) {
            optimotu_wfa_note_fatal();
            return;
          }
          if (d < 1.0) ++my_aligned;

          OPTIMOTU_DEBUG(
            2,
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
        OPTIMOTU_DEBUG(
          2,
          << "thread" << pg_index
          << ": finished " << lj << " (j0=" << j0 << ")"
          << " and " << li << " (i0=" << i0 << ")"
          << "\n" << std::endl
        );
        RcppThread::checkUserInterrupt(); });
      mutex.lock();
      OPTIMOTU_DEBUG(1, << "thread" << pg_index << " ready to merge" << std::endl);
      _aligned += my_aligned;
      _prealigned += my_prealigned;
      mutex.unlock();
      finish_tile(my_algo, pg_index);
      OPTIMOTU_DEBUG(1, << "thread" << pg_index << " done" << std::endl);
    }
  );
}

template <int verbose>
void HybridConcurrentClusterWorker<verbose>::operator()(std::size_t begin, std::size_t end) {
  (void)end;
  EdlibAlignConfig ed_aligner = edlibNewAlignConfig(
    -1,
    EdlibAlignMode::EDLIB_MODE_NW,
    dissimilarity == Dissimilarity::SCORE
      ? EdlibAlignTask::EDLIB_TASK_DISTANCE
      : EdlibAlignTask::EDLIB_TASK_PATH,
    0,
    0
  );
  wfa::WFAlignerEdit wfa_aligner{
    dissimilarity == Dissimilarity::SCORE
      ? wfa::WFAligner::Score
      : wfa::WFAligner::Alignment
  };
  if (wfa_aligner.fatal())
  {
    optimotu_wfa_note_fatal();
    return;
  }
  for_each_assigned_tile(
    pair_generators, begin, end,
    [&](std::size_t pg_index, PairGenerator &pg) {
      OPTIMOTU_DEBUG(
        2,
        << "HybridConcurrentClusterWorker thread " << pg_index
        << " entered" << std::endl
      );
      size_t my_prealigned = 0;
      size_t my_aligned = 0;
      for_each_tile_pair(pg, [&](
                                 std::size_t i0, std::size_t j0, std::size_t, std::size_t)
                         {
        if (optimotu_wfa_fatal_pending()) return;
        OPTIMOTU_DEBUG(
          4,
          << "Thread " << pg_index
          << ": seqs " << j0 << " (j0=" << j0 << ")"
          << " and " << i0 << " (i0=" << i0 << ")"
          << std::endl
        );
        double threshold = clust_algo.max_relevant(
            static_cast<j_t>(i0), static_cast<j_t>(j0));
        OPTIMOTU_DEBUG(
          4,
          << "thread" << pg_index
          << ": max relevant=" << threshold
          << std::endl
        );
        bool is_seqj_longer = seq[j0].size() > seq[i0].size();
        size_t s1 = is_seqj_longer ? i0 : j0;
        size_t s2 = is_seqj_longer ? j0 : i0;
        double l1 = seq[s1].size(), l2 = seq[s2].size();
        OPTIMOTU_DEBUG(
          4,
          << "#### seq " << s1
          << " (l1=" << l1
          << ") and "
          << s2
          << " (l2=" << l2
          << ")"
          << "####" << std::endl
        );

        if (pair_length_ok(
              dissimilarity, AlignmentSpan::GLOBAL, l1, l2, threshold,
              0, 1, 0, 1, 0, 1
            )) {
          ++my_prealigned;
          double maxd1 = dissimilarity == Dissimilarity::SCORE
              ? edlib_score_k(score_ref_length(l1, l2, AlignmentSpan::GLOBAL), threshold)
              : wfa_identity_max_edits(l1, l2, threshold);
          bool is_close = breakpoint >= 1 ? maxd1 < breakpoint : threshold < breakpoint;
          double d;
          if (is_close) {
            WfaIdentityBound bound = choose_alignment_bound(
                dissimilarity, l1, l2, threshold, 0, 1, 0, 1, 0, 1
            );
            wfa_aligner.setHeuristicBandedStatic(bound.min_k, bound.max_k);
            wfa_aligner.setMaxAlignmentSteps(bound.max_alignment_steps);
            OPTIMOTU_DEBUG(
              4,
              << "wfa_aligner min_k=" << bound.min_k
              << " max_k=" << bound.max_k
              << " max score=" << bound.max_alignment_steps
              << std::endl
            );
            d = distance_wfa2<>(seq[s1], seq[s2], wfa_aligner, dissimilarity, 0, 1);
          } else {
            ed_aligner.k = dissimilarity == Dissimilarity::SCORE
                ? (int)maxd1
                : (int)maxd1 + 1;
            d = distance_edlib(seq[s1], seq[s2], ed_aligner, dissimilarity);
          }
          if (wfa_aligner.fatal() || optimotu_wfa_fatal_pending()) {
            optimotu_wfa_note_fatal();
            return;
          }
          if (d < 1.0) ++my_aligned;
          OPTIMOTU_DEBUG(
            4,
            << "thread" << pg_index
            << ": distance=" << d
            << std::endl
          );
          if (d < threshold) {
            clust_algo(static_cast<j_t>(i0), static_cast<j_t>(j0), d);
          }
        }
        OPTIMOTU_DEBUG(
          4,
          << "thread" << pg_index
          << ": finished " << j0 << " (j0=" << j0 << ")"
          << " and " << i0 << " (i0=" << i0 << ")"
          << "\n" << std::endl
        );
        RcppThread::checkUserInterrupt(); });
      mutex.lock();
      _aligned += my_aligned;
      _prealigned += my_prealigned;
      OPTIMOTU_DEBUG(2, << "thread" << pg_index << " done" << std::endl);
      mutex.unlock();
    }
  );
}

template class HybridSplitClusterWorker<0>;
template class HybridSplitClusterWorker<1>;
template class HybridSplitClusterWorker<2>;
template class HybridSplitClusterWorker<3>;
template class HybridSplitClusterWorker<4>;
template class HybridConcurrentClusterWorker<0>;
template class HybridConcurrentClusterWorker<1>;
template class HybridConcurrentClusterWorker<2>;
template class HybridConcurrentClusterWorker<3>;
template class HybridConcurrentClusterWorker<4>;
