#ifndef OPTIMOTU_HYBRIDSEARCHWORKER_H
#define OPTIMOTU_HYBRIDSEARCHWORKER_H

#include "SearchWorker.h"

class HybridSearchWorker : public SearchWorker {
protected:
  double breakpoint = 0.1;
public :
  HybridSearchWorker(
    const std::vector<std::string> &query,
    const std::vector<std::string> &ref,
    const double threshold,
    const std::uint8_t threads,
    const double breakpoint = 0.1
  ) : SearchWorker(query, ref, threshold, threads),
  breakpoint(breakpoint) {};
};

template<int verbose>
class HybridSearchWorkerImpl : public HybridSearchWorker {
  using SearchWorker::query;
  using SearchWorker::ref;
  using SearchWorker::threshold;
  using SearchWorker::threads;
  using SearchWorker::mutex;
  using SearchWorker::_prealigned;
  using SearchWorker::_aligned;
  using SearchWorker::hits;

  using HybridSearchWorker::breakpoint;

public:
  using HybridSearchWorker::HybridSearchWorker;
  void operator()(std::size_t begin, std::size_t end) override;
};

#endif
