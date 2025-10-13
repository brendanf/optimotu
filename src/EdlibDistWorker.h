#ifndef OPTIMOTU_EDLIBDISTWORKER_H
#define OPTIMOTU_EDLIBDISTWORKER_H

#ifdef OPTIMOTU_R

#include "DistWorker.h"
#include "SparseDistanceMatrix.h"
#include <cstdint>
#include "alignment_enums.h"

struct EdlibDistWorker : public DistWorker {
  using DistWorker::DistWorker;
};

template<int verbose, bool is_constrained, enum AlignmentSpan span = AlignmentSpan::GLOBAL, typename SparseDistanceMatrixType = SparseDistanceMatrix>
struct EdlibDistWorkerImpl : public EdlibDistWorker {
  using DistWorker::seq;
  using DistWorker::dist_threshold;
  using DistWorker::sim_threshold;
  using DistWorker::sim_threshold_plus_1;
  using DistWorker::threads;
  using DistWorker::sdm;
  using DistWorker::_prealigned;
  using DistWorker::_aligned;
  
  public:
  using EdlibDistWorker::EdlibDistWorker;

  void operator()(std::size_t begin, std::size_t end) override; 
};

std::unique_ptr<EdlibDistWorker> create_edlib_dist_worker(
  const std::vector<std::string> &seq,
  const double dist_threshold,
  const std::uint8_t threads,
  SparseDistanceMatrix &sdm,
  int verbose = 0,
  AlignmentSpan span = AlignmentSpan::GLOBAL,
  bool constrain = true
);

#endif
#endif