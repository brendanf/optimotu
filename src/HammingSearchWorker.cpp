#include "optimotu.h"
#include "HammingSearchWorker.h"
extern "C" {
#include "defs.h"
}

template<int verbose>
void HammingSearchWorkerImpl<verbose>::operator()(std::size_t begin, std::size_t end) {
  std::size_t begin_i = (begin * query.size()) / threads;
  std::size_t end_i = (end * query.size()) / threads;
  OPTIMOTU_DEBUG(
    1,
    << "HammingSearchWorker thread " << begin
    << " entered; sequences [" << begin_i
    << ", "<< end_i << ")" << std::endl
  );

  std::size_t my_prealigned = 0, my_aligned = 0;
  for (std::size_t i = begin_i; i < end_i; i++) {
    if (!hits[i]) {
      hits[i] = std::make_unique<SearchHit>();
    }
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
      ++my_prealigned;
      double d = pss.dist(i, j + query.size(), min_overlap, ignore_gaps);

      OPTIMOTU_DEBUG(
        2,
        << "Thread " << begin
        << ": distance=" << d
        << std::endl
      );
      if (d == 1.0) continue;
      ++my_aligned;
      if (d < max_dist) {
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

template class HammingSearchWorkerImpl<0>;
template class HammingSearchWorkerImpl<1>;
template class HammingSearchWorkerImpl<2>;
