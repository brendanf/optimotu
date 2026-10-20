#include "EdlibClusterWorker.h"
#include "MappedClusterAlgorithm.h"
#include "tile_pairs.h"
#include <edlib.h>
#include <algorithm>
#include "pairwise_alignment.h"

template <int verbose>
void EdlibSplitClusterWorker<verbose>::operator()(std::size_t begin, std::size_t end) {
  (void)end;

  EdlibAlignConfig ed_aligner = edlibNewAlignConfig(
    -1, EdlibAlignMode::EDLIB_MODE_NW, EdlibAlignTask::EDLIB_TASK_PATH, 0, 0
  );

  for_each_assigned_tile(
    pair_generators, begin, end,
    [&](std::size_t pg_index, PairGenerator &pg) {
      OPTIMOTU_DEBUG(
        2,
        << "EdlibSplit thread " << pg_index
        << " entered" << std::endl
      );

      ClusterAlgorithm *my_algo = tile_algo(&pg);
      if (!my_algo) {
        return;
      }
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
        bool is_seqj_longer = seq[j0].size() > seq[i0].size();
        size_t s1 = is_seqj_longer ? i0 : j0;
        size_t s2 = is_seqj_longer ? j0 : i0;
        double l1 = seq[s1].size(), l2 = seq[s2].size();
        OPTIMOTU_DEBUG(
          4,
          << "thread" << pg_index
          << ": seq " << li
          << " (i0=" << i0
          << ") and " << lj
          << " (j0=" << j0
          << "): s1=" << s1
          << " (l1=" << l1
          << "), s2=" << s2
          << " (l2=" << l2 << ")"
          << " threshold=" << threshold
          << std::endl
        );

        double sim_threshold = 1.0 - threshold;
        if (l1/l2 >= sim_threshold) {
          ++my_prealigned;
          double sim_threshold_plus_1 = 2.0 - threshold;
          double maxd1 = threshold * (l1 + l2) / sim_threshold_plus_1;
          ed_aligner.k = (int)maxd1 + 1;
          double d = distance_edlib(seq[s1], seq[s2], ed_aligner);
          if (d < 1.0) ++my_aligned;

          OPTIMOTU_DEBUG(
            4,
            << "thread" << pg_index
            << ": distance=" << d
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

template<int verbose>
void EdlibConcurrentClusterWorker<verbose>::operator()(std::size_t begin, std::size_t end) {
  (void)end;
  size_t my_prealigned = 0;
  size_t my_aligned = 0;

  EdlibAlignConfig ed_aligner = edlibNewAlignConfig(
    -1,
    EdlibAlignMode::EDLIB_MODE_NW,
    EdlibAlignTask::EDLIB_TASK_PATH,
    0,
    0
  );

  for_each_assigned_tile(
    pair_generators, begin, end,
    [&](std::size_t pg_index, PairGenerator &pg) {
      OPTIMOTU_DEBUG(
        2,
        << "EdlibConcurrent thread " << pg_index
        << " entered" << std::endl
      );
      for_each_tile_pair(pg, [&](
        std::size_t i0, std::size_t j0, std::size_t, std::size_t
      ) {
        OPTIMOTU_DEBUG(
          4,
          << "Thread " << begin
          << ": seqs " << j0
          << " and " << i0
          << std::endl
        );
        double threshold = clust_algo.max_relevant(
            static_cast<j_t>(i0), static_cast<j_t>(j0));
        OPTIMOTU_DEBUG(
          4,
          << "Thread " << begin
          << ": max relevant=" << threshold
          << std::endl
        );
        bool is_seqj_longer = seq[j0].size() > seq[i0].size();
        size_t s1 = is_seqj_longer ? i0 : j0;
        size_t s2 = is_seqj_longer ? j0 : i0;
        double l1 = seq[s1].size(), l2 = seq[s2].size();
        OPTIMOTU_DEBUG(
          4,
          << "#### seq " << i0
          << " (l1=" << l1
          << ") and "<< j0
          << " (l2=" << l2
          <<")####" << std::endl
        );

        double sim_threshold = 1.0 - threshold;
        if (l1/l2 >= sim_threshold) {
          ++my_prealigned;
          double sim_threshold_plus_1 = 2.0 - threshold;
          double maxd1 = threshold * (l1 + l2) / sim_threshold_plus_1;
          ed_aligner.k = (int)maxd1 + 1;
          double d = distance_edlib(seq[s1], seq[s2], ed_aligner);
          if (d < 1.0) ++my_aligned;
          OPTIMOTU_DEBUG(
            4,
            << "Thread " << begin
            << ": distance=" << d
            << std::endl
          );
          if (d < threshold) {
            clust_algo(static_cast<j_t>(i0), static_cast<j_t>(j0), d);
          }
        }
        OPTIMOTU_DEBUG(
          4,
          << "Thread " << begin
          << ": finished " << j0
          << " and " << i0
          << "\n" << std::endl
        );
        RcppThread::checkUserInterrupt();
      });
    }
  );
  mutex.lock();
  _aligned += my_aligned;
  _prealigned += my_prealigned;
  OPTIMOTU_DEBUG(2, << "Exiting thread " << begin << std::endl);
  mutex.unlock();
}

template class EdlibSplitClusterWorker<0>;
template class EdlibSplitClusterWorker<1>;
template class EdlibSplitClusterWorker<2>;
template class EdlibSplitClusterWorker<3>;
template class EdlibSplitClusterWorker<4>;
template class EdlibConcurrentClusterWorker<0>;
template class EdlibConcurrentClusterWorker<1>;
template class EdlibConcurrentClusterWorker<2>;
template class EdlibConcurrentClusterWorker<3>;
template class EdlibConcurrentClusterWorker<4>;
