// SPDX-FileCopyrightText: 2026 Brendan Furneaux
// <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_AASEARCHWORKER_H
#define OPTIMOTU_AASEARCHWORKER_H

#include "EvoModel.h"
#include "SearchWorker.h"
#include "PackedAaSequenceSet.h"
#include <vector>

class AaSearchWorker : public SearchWorker {
protected:
  const PackedAaSequenceSet pss;
  const int min_overlap;
  const bool ignore_gaps;
  std::vector<float> lut_storage;
  const float *weight_lut;
  bool has_evo = false;
  EvoModel evo;
  std::vector<uint8_t> ref_pos_major;
  std::vector<int> ref_idx;
  void score_refs(
    int query_i,
    double *dist_out,
    int *n_ok_out
  ) const;
public:
  AaSearchWorker(
    const std::vector<std::string> &query,
    const std::vector<std::string> &ref,
    const double threshold,
    const std::uint8_t threads,
    const int min_overlap = 0,
    const bool ignore_gaps = true,
    const float *weight_lut = nullptr,
    const EvoModel *evo_model = nullptr
  );
};

template<int verbose>
class AaSearchWorkerImpl : public AaSearchWorker {
  using SearchWorker::query;
  using SearchWorker::ref;
  using SearchWorker::threshold;
  using SearchWorker::threads;
  using SearchWorker::mutex;
  using SearchWorker::_prealigned;
  using SearchWorker::_aligned;
  using SearchWorker::hits;
  using AaSearchWorker::pss;
  using AaSearchWorker::min_overlap;
  using AaSearchWorker::ignore_gaps;
  using AaSearchWorker::weight_lut;
  using AaSearchWorker::has_evo;
  using AaSearchWorker::evo;
  using AaSearchWorker::score_refs;
  using AaSearchWorker::ref_pos_major;
  using AaSearchWorker::ref_idx;
public:
  using AaSearchWorker::AaSearchWorker;
  void operator()(std::size_t begin, std::size_t end) override;
  bool try_run_gpu() override;
};

#endif
