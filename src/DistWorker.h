#ifndef OPTIMOTU_DISTWORKER_H
#define OPTIMOTU_DISTWORKER_H

#ifdef OPTIMOTU_R

#include <RcppParallel.h>
#include <algorithm>
#include <vector>
#include <string>
#include <memory>
#include <cstdint>

#include "SparseDistanceMatrix.h"
#include "PairGenerator.h"

class DistWorker : public RcppParallel::Worker {
protected:
  const std::vector<std::string> &seq;
  const double dist_threshold, sim_threshold, sim_threshold_plus_1;
  std::vector<std::unique_ptr<PairGenerator>> pair_generators;
  const std::size_t threads;
  std::size_t _prealigned = 0, _aligned = 0;
  SparseDistanceMatrix &sdm;

public:
  DistWorker(
    const std::vector<std::string> &seq,
    const double dist_threshold,
    DivisiblePairGenerator::Builder & pair_generator_builder,
    SparseDistanceMatrix &sdm,
    int verbose = 0,
    std::size_t worker_threads = 0
  ) : seq(seq),
      dist_threshold(dist_threshold),
      sim_threshold(1.0 - dist_threshold),
      sim_threshold_plus_1(1.0 + sim_threshold),
      pair_generators(pair_generator_builder.build(verbose)),
      threads(
        worker_threads == 0
          ? std::max<std::size_t>(1, pair_generators.size())
          : std::max<std::size_t>(
              1, std::min(worker_threads, pair_generators.size())
            )
      ),
      sdm(sdm) {}

  virtual ~DistWorker() = default;

  virtual void operator()(std::size_t begin, std::size_t end) = 0;

  std::size_t n_threads() const { return threads; }

  std::size_t prealigned() const { return _prealigned; }

  std::size_t aligned() const { return _aligned; }
};

#endif
#endif
