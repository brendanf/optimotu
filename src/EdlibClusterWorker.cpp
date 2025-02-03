#include "EdlibClusterWorker.h"
#include <edlib.h>
#include "pairwise_alignment.h"

template <int verbose>
void EdlibSplitClusterWorker<verbose>::operator()(std::size_t begin, std::size_t end) {

  double n = seq.size();
  double m = (n*n - 3.0*n + 2.0)/2.0;
  size_t my_prealigned = 0;
  size_t my_aligned = 0;
  size_t begin_i;

  EdlibAlignConfig ed_aligner = edlibNewAlignConfig(-1, EdlibAlignMode::EDLIB_MODE_NW, EdlibAlignTask::EDLIB_TASK_PATH, 0, 0);
  ClusterAlgorithm * my_algo = clust_algo.make_child();

  if (begin == 0) {
    begin_i = 1;
  } else {
    begin_i = round(1.5 + 0.5*sqrt(9.0 + 8.0*((m*begin)/threads - 1.0)));
  }
  size_t end_i   = round(1.5 + 0.5*sqrt(9.0 + 8.0*((m*end)/threads - 1.0)));
  OPTIMOTU_DEBUG(
    1,
    << "EdlibSplit thread " << begin
    << " entered; sequences [" << begin_i
    << ", "<< end_i << ")" << std::endl;
  );
  for (size_t i = begin_i; i < end_i; i++) {
    for (size_t j = 0; j < i; j++) {
      double threshold = my_algo->max_relevant(i, j);
      OPTIMOTU_DEBUG(
        2,
        << "thread" << begin
        << ": seqs " << j
        << " and " << i
        << " max relevant=" << threshold
        << std::endl
      );

      bool is_seqj_longer = seq[j].size() > seq[i].size();
      size_t s1 = is_seqj_longer ? i : j;
      size_t s2 = is_seqj_longer ? j : i;
      double l1 = seq[s1].size(), l2 = seq[s2].size();
      OPTIMOTU_DEBUG(
        2,
        << "thread" << begin
        << ": seq " << i << " (l1=" << l1 << ") and "
        << j << " (l2=" << l2 <<")"
        << std::endl
      );

      double sim_threshold = 1.0 - threshold; // compiler can probably do this?
      if (l1/l2 < sim_threshold) continue;
      ++my_prealigned;
      double sim_threshold_plus_1 = 2.0 - threshold;
      double maxd1 = threshold * (l1 + l2) / sim_threshold_plus_1;
      ed_aligner.k = (int)maxd1 + 1;
      double d = distance_edlib(seq[s1], seq[s2], ed_aligner);
      if (d < 1.0) ++my_aligned;

      OPTIMOTU_DEBUG(
        2,
        << "thread" << begin
        << ": distance=" << d
        << std::endl
      );
      if (d < threshold) (*my_algo)(j, i, d);
      RcppThread::checkUserInterrupt();
    }
  }
  mutex.lock();
  OPTIMOTU_DEBUG(1, << "thread " << begin << " ready to merge" << std::endl);
  _aligned += my_aligned;
  _prealigned += my_prealigned;
  mutex.unlock();
  my_algo->merge_into_parent();
  OPTIMOTU_DEBUG(1, << "thread " << begin << " done" << std::endl);
}

template<int verbose>
void EdlibConcurrentClusterWorker<verbose>::operator()(std::size_t begin, std::size_t end) {
  double n = seq.size();
  double m = (n*n - 3.0*n + 2.0)/2.0;
  size_t my_prealigned = 0;
  size_t my_aligned = 0;
  size_t begin_i;

  EdlibAlignConfig ed_aligner = edlibNewAlignConfig(
    -1,
    EdlibAlignMode::EDLIB_MODE_NW,
    EdlibAlignTask::EDLIB_TASK_PATH,
    0,
    0
  );

  if (begin == 0) {
    begin_i = 1;
  } else {
    begin_i = round(1.5 + 0.5*sqrt(9.0 + 8.0*((m*begin)/threads - 1.0)));
  }
  size_t end_i   = round(1.5 + 0.5*sqrt(9.0 + 8.0*((m*end)/threads - 1.0)));
  OPTIMOTU_DEBUG(
    1,
    << "EdlibConcurrent thread " << begin
    << " entered; sequences [" << begin_i
    << ", "<< end_i << ")" << std::endl;
  );
  for (size_t i = begin_i; i < end_i; i++) {
    for (size_t j = 0; j < i; j++) {
      OPTIMOTU_DEBUG(
        2,
        << "Thread " << begin
        << ": seqs " << j
        << " and " << i
        << std::endl
      );
      double threshold = clust_algo.max_relevant(i, j);
      OPTIMOTU_DEBUG(
        2,
        << "Thread " << begin
        << ": max relevant=" << threshold
        << std::endl
      );
      bool is_seqj_longer = seq[j].size() > seq[i].size();
      size_t s1 = is_seqj_longer ? i : j;
      size_t s2 = is_seqj_longer ? j : i;
      double l1 = seq[s1].size(), l2 = seq[s2].size();
      OPTIMOTU_DEBUG(
        2,
        << "#### seq " << i
        << " (l1=" << l1
        << ") and "<< j
        << " (l2=" << l2
        <<")####" << std::endl
      );

      double sim_threshold = 1.0 - threshold; // compiler can probably do this?
      if (l1/l2 < sim_threshold) continue;
      ++my_prealigned;
      double sim_threshold_plus_1 = 2.0 - threshold;
      double maxd1 = threshold * (l1 + l2) / sim_threshold_plus_1;
      ed_aligner.k = (int)maxd1 + 1;
      double d = distance_edlib(seq[s1], seq[s2], ed_aligner);
      if (d < 1.0) ++my_aligned;
      OPTIMOTU_DEBUG(
        2,
        << "Thread " << begin
        << ": distance=" << d
        << std::endl
      );
      if (d < threshold) clust_algo(j, i, d);
      OPTIMOTU_DEBUG(
        2,
        << "Thread " << begin
        << ": finished " << j
        << " and " << i
        << "\n" << std::endl
      );
      RcppThread::checkUserInterrupt();
    }
  }
  mutex.lock();
  _aligned += my_aligned;
  _prealigned += my_prealigned;
  OPTIMOTU_DEBUG(1, << "Exiting thread " << begin << std::endl);
  mutex.unlock();
}

template class EdlibSplitClusterWorker<0>;
template class EdlibSplitClusterWorker<1>;
template class EdlibSplitClusterWorker<2>;
template class EdlibConcurrentClusterWorker<0>;
template class EdlibConcurrentClusterWorker<1>;
template class EdlibConcurrentClusterWorker<2>;
