#ifndef OPTIMOTU_WFA2SEARCHWORKER_H
#define OPTIMOTU_WFA2SEARCHWORKER_H

#include "SearchWorker.h"

class Wfa2SearchWorker : public SearchWorker {
protected:
  int match = 0, mismatch = 1,
    gap_open = 0, gap_extend = 1,
    gap_open2 = 0, gap_extend2 = 1;
public:
  Wfa2SearchWorker(
    const std::vector<std::string> &query,
    const std::vector<std::string> &ref,
    const double threshold,
    const std::uint8_t threads,
    const int match = 0, const int mismatch = 1,
    const int gap_open = 0, const int gap_extend = 1,
    const int gap_open2 = 0, const int gap_extend2 = 1
  ) : SearchWorker(query, ref, threshold, threads),
  match(match), mismatch(mismatch), gap_open(gap_open),
  gap_extend(gap_extend), gap_open2(gap_open2), gap_extend2(gap_extend2) {}
};

template<int verbose>
class Wfa2SearchWorkerImpl : public Wfa2SearchWorker {
  using SearchWorker::query;
  using SearchWorker::ref;
  using SearchWorker::threshold;
  using SearchWorker::threads;
  using SearchWorker::mutex;
  using SearchWorker::_prealigned;
  using SearchWorker::_aligned;
  using SearchWorker::hits;

  using Wfa2SearchWorker::match;
  using Wfa2SearchWorker::mismatch;
  using Wfa2SearchWorker::gap_open;
  using Wfa2SearchWorker::gap_extend;
  using Wfa2SearchWorker::gap_open2;
  using Wfa2SearchWorker::gap_extend2;

public:
  using Wfa2SearchWorker::Wfa2SearchWorker;
  void operator()(std::size_t begin, std::size_t end) override;
};

#endif
