#include "EdlibClusterWorker.h"
#include <edlib.h>
#include <algorithm>
#include "pairwise_alignment.h"

template <int verbose>
void EdlibSplitClusterWorker<verbose>::operator()(std::size_t begin, std::size_t end) {

  size_t my_prealigned = 0;
  size_t my_aligned = 0;

  EdlibAlignConfig ed_aligner = edlibNewAlignConfig(-1, EdlibAlignMode::EDLIB_MODE_NW, EdlibAlignTask::EDLIB_TASK_PATH, 0, 0);

  for (size_t pg_index = begin; pg_index < pair_generators.size(); pg_index += threads) {
    OPTIMOTU_DEBUG(
      2,
      << "EdlibSplit thread " << pg_index
      << " entered" << std::endl
    );

    auto & pg = pair_generators[pg_index];
    ClusterAlgorithm * my_algo = clust_algo.make_child();
    // iterate over all pairs in the pair generator
    while (*pg) {
      // i, j are the indices of the pair in the cluster algorithm
      size_t i = pg->i();
      size_t j = pg->j();
      // i0, j0 are the indices of the pair in the original sequence list
      size_t i0 = pg->i0();
      size_t j0 = pg->j0();

      double threshold = my_algo->max_relevant(*pg);
      bool is_seqj_longer = seq[j0].size() > seq[i0].size();
      // s1, s2 are the indices of the shorter and longer sequences in
      // the original sequence list
      size_t s1 = is_seqj_longer ? i0 : j0;
      size_t s2 = is_seqj_longer ? j0 : i0;
      double l1 = seq[s1].size(), l2 = seq[s2].size();
      OPTIMOTU_DEBUG(
        4,
        << "thread" << pg_index
        << ": seq " << i
        << " (i0=" << i0
        << ") and " << j
        << " (j0=" << j0
        << "): s1=" << s1
        << " (l1=" << l1
        << "), s2=" << s2
        << " (l2=" << l2 << ")"
        << " threshold=" << threshold
        << std::endl
      );

      double sim_threshold = 1.0 - threshold; // compiler can probably do this?
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
        if (d < threshold) (*my_algo)(*pg, d);
      }
      ++(*pg);
      RcppThread::checkUserInterrupt();
    }
    mutex.lock();
    OPTIMOTU_DEBUG(2, << "thread " << pg_index << " ready to merge" << std::endl);
    _aligned += my_aligned;
    _prealigned += my_prealigned;
    mutex.unlock();
    my_algo->merge_into_parent();
    clust_algo.release_child(my_algo);
    OPTIMOTU_DEBUG(2, << "thread " << pg_index << " done" << std::endl);
  }
}

template<int verbose>
void EdlibConcurrentClusterWorker<verbose>::operator()(std::size_t begin, std::size_t end) {
  size_t my_prealigned = 0;
  size_t my_aligned = 0;

  EdlibAlignConfig ed_aligner = edlibNewAlignConfig(
    -1,
    EdlibAlignMode::EDLIB_MODE_NW,
    EdlibAlignTask::EDLIB_TASK_PATH,
    0,
    0
  );

  for (size_t pg_index = begin; pg_index < pair_generators.size(); pg_index += threads) {
    OPTIMOTU_DEBUG(
      2,
      << "EdlibConcurrent thread " << pg_index
      << " entered" << std::endl
    );
    auto & pg = pair_generators[pg_index];
    while (*pg) {
      size_t i = pg->i0();
      size_t j = pg->j0();
      OPTIMOTU_DEBUG(
        4,
        << "Thread " << begin
        << ": seqs " << j
        << " and " << i
        << std::endl
      );
      double threshold = clust_algo.max_relevant(*pg);
      OPTIMOTU_DEBUG(
        4,
        << "Thread " << begin
        << ": max relevant=" << threshold
        << std::endl
      );
      bool is_seqj_longer = seq[j].size() > seq[i].size();
      size_t s1 = is_seqj_longer ? i : j;
      size_t s2 = is_seqj_longer ? j : i;
      double l1 = seq[s1].size(), l2 = seq[s2].size();
      OPTIMOTU_DEBUG(
        4,
        << "#### seq " << i
        << " (l1=" << l1
        << ") and "<< j
        << " (l2=" << l2
        <<")####" << std::endl
      );

      double sim_threshold = 1.0 - threshold; // compiler can probably do this?
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
        if (d < threshold) clust_algo(*pg, d);
      }
      OPTIMOTU_DEBUG(
        4,
        << "Thread " << begin
        << ": finished " << j
        << " and " << i
        << "\n" << std::endl
      );
      ++(*pg);
      RcppThread::checkUserInterrupt();
    }
  }
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