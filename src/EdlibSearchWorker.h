#ifndef OPTIMOTU_EDLIBSEARCHWORKER_H
#define OPTIMOTU_EDLIBSEARCHWORKER_H

#include "SearchWorker.h"

class EdlibSearchWorker : public SearchWorker {
public :
  using SearchWorker::SearchWorker;
};

template<int verbose>
class EdlibSearchWorkerImpl : public EdlibSearchWorker {
  using SearchWorker::query;
  using SearchWorker::ref;
  using SearchWorker::threshold;
  using SearchWorker::threads;
  using SearchWorker::mutex;
  using SearchWorker::_prealigned;
  using SearchWorker::_aligned;
  using SearchWorker::hits;
public:
  using EdlibSearchWorker::EdlibSearchWorker;
  void operator()(std::size_t begin, std::size_t end) override;
};

#endif
