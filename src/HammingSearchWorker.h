#ifndef OPTIMOTU_HAMMINGSEARCHWORKER_H
#define OPTIMOTU_HAMMINGSEARCHWORKER_H

#include "SearchWorker.h"
#include "PackedSequenceSet.h"

class HammingSearchWorker : public SearchWorker {
protected:
  const PackedSequenceSet pss;
  const int min_overlap;
  const bool ignore_gaps;
public :
  HammingSearchWorker(
    const std::vector<std::string> &query,
    const std::vector<std::string> &ref,
    const double threshold,
    const std::uint8_t threads,
    const int min_overlap = 0,
    const bool ignore_gaps = TRUE
  ) : SearchWorker(query, ref, threshold, threads),
  pss(query, ref), min_overlap(min_overlap),
  ignore_gaps(ignore_gaps) {};
};

template<int verbose>
class HammingSearchWorkerImpl : public HammingSearchWorker {
  using SearchWorker::query;
  using SearchWorker::ref;
  using SearchWorker::threshold;
  using SearchWorker::threads;
  using SearchWorker::mutex;
  using SearchWorker::_prealigned;
  using SearchWorker::_aligned;
  using SearchWorker::hits;

  using HammingSearchWorker::pss;
  using HammingSearchWorker::min_overlap;
  using HammingSearchWorker::ignore_gaps;

public:
  using HammingSearchWorker::HammingSearchWorker;
  void operator()(std::size_t begin, std::size_t end) override;
};

#endif
