// SPDX-FileCopyrightText: 2026 Brendan Furneaux
// <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_AACLUSTERWORKER_H_INCLUDED
#define OPTIMOTU_AACLUSTERWORKER_H_INCLUDED

#include "DistClusterWorker.h"
#include "PackedAaSequenceSet.h"
#include "SequenceView.h"
#include <cstddef>
#include <vector>

class AaClusterWorker : public DistClusterWorker {
protected:
  const PackedAaSequenceSet pss;
  const int min_overlap;
  const bool ignore_gaps;
  std::vector<float> lut_storage;
  const float *weight_lut;
  std::size_t tracked_pss_bytes = 0;
public:
  AaClusterWorker(
    const SequenceSet &seq,
    ClusterAlgorithm &clust_algo,
    DivisiblePairGenerator::Builder &pgb,
    const int min_overlap = 0,
    const bool ignore_gaps = true,
    const float *weight_lut = nullptr,
    int verbose = 0,
    std::size_t worker_threads = 0
  );
  ~AaClusterWorker() override;
};

template<int verbose>
class AaSplitClusterWorker : public AaClusterWorker {
  using AaClusterWorker::pss;
  using DistClusterWorker::clust_algo;
  using DistClusterWorker::pair_generators;
  using DistClusterWorker::mutex;
  using AaClusterWorker::min_overlap;
  using AaClusterWorker::ignore_gaps;
  using AaClusterWorker::weight_lut;
  using DistClusterWorker::_prealigned;
  using DistClusterWorker::_aligned;
public:
  AaSplitClusterWorker(
    const SequenceSet &seq,
    ClusterAlgorithm &clust_algo,
    DivisiblePairGenerator::Builder &pgb,
    const int min_overlap = 0,
    const bool ignore_gaps = true,
    const float *weight_lut = nullptr,
    std::size_t worker_threads = 0
  );
  void operator()(std::size_t begin, std::size_t end) override;
};

template<int verbose>
class AaConcurrentClusterWorker : public AaClusterWorker {
  using AaClusterWorker::pss;
  using DistClusterWorker::clust_algo;
  using DistClusterWorker::pair_generators;
  using DistClusterWorker::mutex;
  using AaClusterWorker::min_overlap;
  using AaClusterWorker::ignore_gaps;
  using AaClusterWorker::weight_lut;
  using DistClusterWorker::_prealigned;
  using DistClusterWorker::_aligned;
public:
  AaConcurrentClusterWorker(
    const SequenceSet &seq,
    ClusterAlgorithm &clust_algo,
    DivisiblePairGenerator::Builder &pgb,
    const int min_overlap = 0,
    const bool ignore_gaps = true,
    const float *weight_lut = nullptr,
    std::size_t worker_threads = 0
  );
  void operator()(std::size_t begin, std::size_t end) override;
};

#endif
