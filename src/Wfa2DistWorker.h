#ifndef OPTIMOTU_WFA2DISTWORKER_H
#define OPTIMOTU_WFA2DISTWORKER_H

#ifdef OPTIMOTU_R
#include "DistWorker.h"
#include "alignment_enums.h"

class Wfa2DistWorker : public DistWorker {
protected:
  int match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2;
public:
  Wfa2DistWorker(
    const std::vector<std::string> &seq,
    const double dist_threshold,
    const std::uint8_t threads,
    SparseDistanceMatrix &sdm,
    int match, int mismatch,
    int gap_open, int gap_extend,
    int gap_open2, int gap_extend2
  );
};

template<int verbose, bool is_constrained, enum AlignmentSpan span = AlignmentSpan::GLOBAL, typename SparseDistanceMatrixType = SparseDistanceMatrix>
class Wfa2DistWorkerImpl : public Wfa2DistWorker {
  using DistWorker::seq;
  using DistWorker::dist_threshold;
  using DistWorker::threads;
  using DistWorker::sdm;
  using Wfa2DistWorker::match;
  using Wfa2DistWorker::mismatch;
  using Wfa2DistWorker::gap_open;
  using Wfa2DistWorker::gap_extend;
  using Wfa2DistWorker::gap_open2;
  using Wfa2DistWorker::gap_extend2;
public:
  using Wfa2DistWorker::Wfa2DistWorker;

  virtual void operator()(std::size_t begin, std::size_t end) override;
};

std::unique_ptr<Wfa2DistWorker> create_wfa2_dist_worker(
  const std::vector<std::string> &seq,
  const double dist_threshold,
  const std::uint8_t threads,
  SparseDistanceMatrix &sdm,
  int match = 0,
  int mismatch = 1,
  int gap_open = 0,
  int gap_extend = 1,
  int gap_open2 = 0,
  int gap_extend2 = 1,
  int verbose = 0,
  enum AlignmentSpan span = AlignmentSpan::GLOBAL,
  bool constrain = true
);

#endif //OPTIMOTU_R

#endif //OPTIMOTU_WFA2DISTWORKER_H
