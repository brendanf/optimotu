#ifndef OPTIMOTU_HAMMINGDISTWORKER_H_INCLUDED
#define OPTIMOTU_HAMMINGDISTWORKER_H_INCLUDED

#ifdef OPTIMOTU_R

#include <memory>

#include "DistWorker.h"
#include "PackedSequenceSet.h"

class HammingDistWorker : public DistWorker {
protected:
  const PackedSequenceSet pss;
  const int min_overlap;
  const bool ignore_gap;

public:
  HammingDistWorker(
    const std::vector<std::string> &seq,
    const double dist_threshold,
    const std::uint8_t threads,
    SparseDistanceMatrix &sdm,
    const int min_overlap,
    const bool ignore_gap
  );
};

template<int verbose, typename SparseDistanceMatrixType = SparseDistanceMatrix>
class HammingDistWorkerImpl : public HammingDistWorker {
  using DistWorker::seq;
  using DistWorker::dist_threshold;
  using DistWorker::threads;
  using DistWorker::sdm;

  using HammingDistWorker::pss;
  using HammingDistWorker::min_overlap;
  using HammingDistWorker::ignore_gap;

public:
  using HammingDistWorker::HammingDistWorker;
  void operator()(std::size_t begin, std::size_t end) override;
};

std::unique_ptr<HammingDistWorker> create_hamming_dist_worker(
  const std::vector<std::string> &seq,
  const double dist_threshold,
  const std::uint8_t threads,
  SparseDistanceMatrix &sdm,
  const int min_overlap,
  const bool ignore_gap,
  int verbose = 0
);

#endif //OPTIMOTU_R

#endif //OPTIMOTU_HAMMINGDISTWORKER_H_INCLUDED
