#include "optimotu.h"
#include "EdlibSearchWorker.h"
#include <edlib.h>
#include "pairwise_alignment.h"

template<int verbose, enum AlignmentSpan span>
void EdlibSearchWorkerImpl<verbose, span>::operator()(std::size_t begin, std::size_t end) {
  std::size_t begin_i = (begin * query.size()) / threads;
  std::size_t end_i = (end * query.size()) / threads;
  OPTIMOTU_DEBUG(
    1,
    << "EdlibSearchWorker thread " << begin
    << " entered; sequences [" << begin_i
    << ", "<< end_i << ")" << std::endl
  );

  EdlibAlignMode mode;
  if constexpr (span == AlignmentSpan::GLOBAL) {
    mode = EdlibAlignMode::EDLIB_MODE_NW;
  } else if constexpr (span == AlignmentSpan::EXTEND) {
    mode = EdlibAlignMode::EDLIB_MODE_SHW;
  } else {
    static_assert(span != span, "Invalid alignment span");
  }
  EdlibAlignConfig ed_aligner = edlibNewAlignConfig(-1, mode, EdlibAlignTask::EDLIB_TASK_PATH, 0, 0);

  std::size_t my_prealigned = 0, my_aligned = 0;
  for (std::size_t i = begin_i; i < end_i; i++) {
    SearchHit & hit = *hits[i];
    for (std::size_t j = 0; j < ref.size(); j++) {
      double max_dist = (hit.best_dist < threshold) ? hit.best_dist : threshold;
      OPTIMOTU_DEBUG(
        2,
        << "thread" << begin
        << ": seqs " << j
        << " and " << i
        << " max relevant=" << max_dist
        << std::endl
      );

      bool is_query_longer = query[i].size() > ref[j].size();
      const std::string & s1 = is_query_longer ? ref[j] : query[i];
      const std::string & s2 = is_query_longer ? query[i] : ref[j];

      double l1 = s1.size(), l2 = s2.size();
      OPTIMOTU_DEBUG(
        2,
        << (is_query_longer ? "#### ref " : "#### query " )
        << (is_query_longer ? j : i)
        << " (l1=" << l1 << ") and "
        << (is_query_longer ? "query " : "ref " )
        << (is_query_longer ? i : j)
        << " (l2=" << l2 <<")####" << std::endl
      );

      double sim_threshold = 1.0 - max_dist; // compiler can probably do this?
      if (l1/l2 < sim_threshold) continue;
      ++my_prealigned;
      double sim_threshold_plus_1 = 2.0 - threshold;
      double maxd1 = threshold * (l1 + l2) / sim_threshold_plus_1;
      ed_aligner.k = (int)maxd1 + 1;
      double d = distance_edlib(s1, s2, ed_aligner);
      OPTIMOTU_DEBUG(
        2,
        << "Thread " << begin
        << ": distance=" << d
        << std::endl
      );
      if (d == 1.0) continue;
      ++my_aligned;
      if (d < hit.best_dist) {
        OPTIMOTU_DEBUG(
          2,
          << "Thread " << begin
          << ": new best distance=" << d
          << std::endl
        );
        hit.best_dist = d;
        hit.best_ref.clear();
        hit.best_ref.push_back(j);
      } else if (d == hit.best_dist) {
        OPTIMOTU_DEBUG(
          2,
          << "Thread " << begin
          << ": tied for best distance=" << d
          << std::endl
        );
        hit.best_ref.push_back(j);
      }
      RcppThread::checkUserInterrupt();
    }
  }
  {
    std::lock_guard<std::mutex> lock(mutex);
    _prealigned += my_prealigned;
    _aligned += my_aligned;
  }
  OPTIMOTU_DEBUG(1, << "Exiting thread " << begin << std::endl);
}

template class EdlibSearchWorkerImpl<0>;
template class EdlibSearchWorkerImpl<1>;
template class EdlibSearchWorkerImpl<2>;
