#ifndef OPTIMOTU_DISTWORKER_H
#define OPTIMOTU_DISTWORKER_H

#ifdef OPTIMOTU_R

#include <RcppParallel.h>
#include <algorithm>
#include <vector>
#include <string>
#include <memory>
#include <cstdint>
#include <functional>

#include "SparseDistanceMatrix.h"
#include "PairGenerator.h"
#include "alignment_enums.h"
#include "gpu_tile.h"

class DistWorker : public RcppParallel::Worker {
protected:
  const std::vector<std::string> &seq;
  const double dist_threshold, sim_threshold, sim_threshold_plus_1;
  std::vector<std::unique_ptr<PairGenerator>> pair_generators;
  const std::size_t threads;
  std::size_t _prealigned = 0, _aligned = 0;
  SparseDistanceMatrix &sdm;
  const Dissimilarity dissimilarity;

public:
  DistWorker(
    const std::vector<std::string> &seq,
    const double dist_threshold,
    DivisiblePairGenerator::Builder & pair_generator_builder,
    SparseDistanceMatrix &sdm,
    int verbose = 0,
    std::size_t worker_threads = 0,
    Dissimilarity dissimilarity = Dissimilarity::IDENTITY
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
      sdm(sdm),
      dissimilarity(dissimilarity) {}

  virtual ~DistWorker() = default;

  virtual void operator()(std::size_t begin, std::size_t end) = 0;

  virtual bool wants_gpu() const { return false; }
  virtual std::unique_ptr<GpuDistanceKernel> make_gpu_kernel() {
    return nullptr;
  }
  virtual void consume_gpu_tile(
    std::size_t pg_index,
    PairGenerator &pg,
    const GpuTileResult &result
  ) {
    (void)pg_index;
    (void)pg;
    (void)result;
  }

  std::size_t n_threads() const { return threads; }

  std::size_t n_tiles() const { return pair_generators.size(); }

  void run() {
    const std::size_t n = pair_generators.size();
    if (n == 0) {
      return;
    }
    if (wants_gpu() && gpu_runtime_available() &&
        !gpu_selected_devices().empty()) {
      gpu_run_tiles(
        pair_generators,
        std::max(threads, gpu_selected_devices().size()),
        [this]() { return make_gpu_kernel(); },
        [this](
          std::size_t idx,
          PairGenerator &pg,
          const GpuTileResult &res
        ) {
          consume_gpu_tile(idx, pg, res);
        }
      );
      return;
    }
    if (threads <= 1) {
      (*this)(0, n);
      return;
    }
    RcppParallel::parallelFor(0, n, *this, 1, threads);
  }

  std::size_t prealigned() const { return _prealigned; }

  std::size_t aligned() const { return _aligned; }
};

#endif
#endif
