#include "HybridClusterWorker.h"
#include "ClusterAlgorithmFactory.h"
#include "config.h"

#include <bindings/cpp/WFAligner.hpp>
#include <edlib.h>
#include "pairwise_alignment.h"

typedef RcppParallel::RMatrix<int> matrix_t;

HybridClusterWorker::HybridClusterWorker(
  const std::vector<std::string> &seq,
  ClusterAlgorithm &clust_algo,
  const std::uint8_t threads,
  const double breakpoint
) : AlignClusterWorker(seq, clust_algo, threads), breakpoint(breakpoint) {};

template<int verbose>
void HybridSplitClusterWorker<verbose>::operator()(std::size_t begin, std::size_t end) {
  double n = seq.size();
  double m = (n*n - 3.0*n + 2.0)/2.0;
  size_t my_prealigned = 0;
  size_t my_aligned = 0;
  size_t begin_i;

  EdlibAlignConfig ed_aligner = edlibNewAlignConfig(-1, EdlibAlignMode::EDLIB_MODE_NW, EdlibAlignTask::EDLIB_TASK_PATH, 0, 0);
  wfa::WFAlignerEdit wfa_aligner{wfa::WFAligner::Alignment};
  ClusterAlgorithm * my_algo = clust_algo.make_child();

  if (begin == 0) {
    begin_i = 1;
  } else {
    begin_i = round(1.5 + 0.5*sqrt(9.0 + 8.0*((m*begin)/threads - 1.0)));
  }
  size_t end_i   = round(1.5 + 0.5*sqrt(9.0 + 8.0*((m*end)/threads - 1.0)));
  OPTIMOTU_DEBUG(
    1,
    << "HybridSplit thread " << begin
    << " entered; sequences [" << begin_i
    << ", "<< end_i << ")" << std::endl
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
        << "#### seq " << i << " (l1=" << l1 << ") and "
        << j << " (l2=" << l2 <<")####" << std::endl
      );

      double sim_threshold = 1.0 - threshold; // compiler can probably do this?
      if (l1/l2 < sim_threshold) continue;
      ++my_prealigned;
      double sim_threshold_plus_1 = 2.0 - threshold;
      double maxd1 = threshold * (l1 + l2) / sim_threshold_plus_1;
      bool is_close = breakpoint >= 1 ? maxd1 < breakpoint : threshold < breakpoint;
      double d;
      if (is_close) {
        int max_k = (int)ceil((l2 - l1 * sim_threshold) / sim_threshold_plus_1);
        int min_k = -(int)ceil((l1 - l2 * sim_threshold) / sim_threshold_plus_1);
        wfa_aligner.setHeuristicBandedStatic(min_k, max_k);
        wfa_aligner.setMaxAlignmentSteps((int)maxd1 + 1);
        d = distance_wfa2(seq[s1], seq[s2], wfa_aligner);
      } else {
        ed_aligner.k = (int)maxd1 + 1;
        d = distance_edlib(seq[s1], seq[s2], ed_aligner);
      }
      if (d < 1.0) ++my_aligned;

      OPTIMOTU_DEBUG(
        2,
        << " distance=" << d
        << std::endl
      );
      if (d < threshold) (*my_algo)(j, i, d);
      OPTIMOTU_DEBUG(
        2,
        << "thread" << begin
        << ": finished " << j
        << " and " << i
        << "\n" << std::endl
      );
      RcppThread::checkUserInterrupt();
    }
  }
  mutex.lock();
  OPTIMOTU_DEBUG(1, << "thread" << begin << " ready to merge" << std::endl);
  _aligned += my_aligned;
  _prealigned += my_prealigned;
  mutex.unlock();
  my_algo->merge_into_parent();
  OPTIMOTU_DEBUG(1, << "thread" << begin << " done" << std::endl);
}

template <int verbose>
void HybridConcurrentClusterWorker<verbose>::operator()(std::size_t begin, std::size_t end) {
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
  wfa::WFAlignerEdit wfa_aligner{wfa::WFAligner::Alignment};

  if (begin == 0) {
    begin_i = 1;
  } else {
    begin_i = round(1.5 + 0.5*sqrt(9.0 + 8.0*((m*begin)/threads - 1.0)));
  }
  size_t end_i   = round(1.5 + 0.5*sqrt(9.0 + 8.0*((m*end)/threads - 1.0)));
  OPTIMOTU_DEBUG(
    1,
    << "HybridConcurrent thread " << begin
    << " entered; sequences [" << begin_i
    << ", "<< end_i << ")" << std::endl
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
        << "#### seq " << i << " (l1=" << l1 << ") and "
        << j << " (l2=" << l2 <<")####" << std::endl
      );

      double sim_threshold = 1.0 - threshold; // compiler can probably do this?
      if (l1/l2 < sim_threshold) continue;
      ++my_prealigned;
      double sim_threshold_plus_1 = 2.0 - threshold;
      double maxd1 = threshold * (l1 + l2) / sim_threshold_plus_1;
      bool is_close = breakpoint >= 1 ? maxd1 < breakpoint : threshold < breakpoint;
      double d;
      if (is_close) {
        int max_k = (int)ceil((l2 - l1 * sim_threshold) / sim_threshold_plus_1);
        int min_k = -(int)ceil((l1 - l2 * sim_threshold) / sim_threshold_plus_1);
        wfa_aligner.setHeuristicBandedStatic(min_k, max_k);
        wfa_aligner.setMaxAlignmentSteps((int)maxd1 + 1);
        OPTIMOTU_DEBUG(
          2,
          << "wfa_aligner min_k=" << min_k
          << " max_k=" << max_k
          << " max score=" << (int)maxd1 + 1
          << std::endl
        );
        d = distance_wfa2(seq[s1], seq[s2], wfa_aligner);
      } else {
        ed_aligner.k = (int)maxd1 + 1;
        d = distance_edlib(seq[s1], seq[s2], ed_aligner);
      }
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

template class HybridSplitClusterWorker<0>;
template class HybridSplitClusterWorker<1>;
template class HybridSplitClusterWorker<2>;
template class HybridConcurrentClusterWorker<0>;
template class HybridConcurrentClusterWorker<1>;
template class HybridConcurrentClusterWorker<2>;
