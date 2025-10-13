#ifndef OPTIMOTU_DISTWORKER_H
#define OPTIMOTU_DISTWORKER_H

#ifdef OPTIMOTU_R

#include <RcppParallel.h>
#include <vector>
#include <string>
#include <memory>
#include <cstdint>

#include "SparseDistanceMatrix.h"

class DistWorker : public RcppParallel::Worker {
protected:
  const std::vector<std::string> &seq;
  const double dist_threshold, sim_threshold, sim_threshold_plus_1;
  const std::uint8_t threads;
  std::size_t _prealigned = 0, _aligned = 0;
  SparseDistanceMatrix &sdm;

public:
  DistWorker(
    const std::vector<std::string> &seq,
    const double dist_threshold,
    const std::uint8_t threads,
    SparseDistanceMatrix &sdm
  ) : seq(seq),
      dist_threshold(dist_threshold),
      sim_threshold(1.0 - dist_threshold),
      sim_threshold_plus_1(1.0 + sim_threshold),
      threads(threads),
      sdm(sdm) {}

  virtual ~DistWorker() = default;

  virtual void operator()(std::size_t begin, std::size_t end) = 0;

  std::uint8_t n_threads() const { return threads; }

  std::size_t prealigned() const { return _prealigned; }

  std::size_t aligned() const { return _aligned; }
};

#endif
#endif