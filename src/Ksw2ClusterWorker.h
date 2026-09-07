// SPDX-FileCopyrightText: 2026 Brendan Furneaux
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_KSW2CLUSTERWORKER_H
#define OPTIMOTU_KSW2CLUSTERWORKER_H

#include "DistClusterWorker.h"
#include "SequenceView.h"

class Ksw2ClusterWorker : public DistClusterWorker {
protected:
  int match = 0, mismatch = 1,
    gap_open = 0, gap_extend = 1,
    gap_open2 = 0, gap_extend2 = 1;
public :
  Ksw2ClusterWorker(
    const SequenceSet &seq,
    ClusterAlgorithm &clust_algo,
    DivisiblePairGenerator::Builder & pgb,
    const int match = 0, const int mismatch = 1,
    const int gap_open = 0, const int gap_extend = 1,
    const int gap_open2 = 0, const int gap_extend2 = 1,
    int verbose = 0,
    std::size_t worker_threads = 0
  );
};

template<int verbose>
class Ksw2SplitClusterWorker : public Ksw2ClusterWorker {
  using DistClusterWorker::seq;
  using DistClusterWorker::clust_algo;
  using DistClusterWorker::threads;
  using DistClusterWorker::mutex;
  using DistClusterWorker::_prealigned;
  using DistClusterWorker::_aligned;

  using Ksw2ClusterWorker::match;
  using Ksw2ClusterWorker::mismatch;
  using Ksw2ClusterWorker::gap_open;
  using Ksw2ClusterWorker::gap_extend;
  using Ksw2ClusterWorker::gap_open2;
  using Ksw2ClusterWorker::gap_extend2;
public:
  Ksw2SplitClusterWorker(
    const SequenceSet &seq,
    ClusterAlgorithm &clust_algo,
    DivisiblePairGenerator::Builder & pgb,
    const int match = 0, const int mismatch = 1,
    const int gap_open = 0, const int gap_extend = 1,
    const int gap_open2 = 0, const int gap_extend2 = 1,
    std::size_t worker_threads = 0
  ) : Ksw2ClusterWorker(seq, clust_algo, pgb, match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2, verbose, worker_threads) {};
  void operator()(std::size_t begin, std::size_t end) override;
};

template<int verbose>
class Ksw2ConcurrentClusterWorker : public Ksw2ClusterWorker {
  using DistClusterWorker::seq;
  using DistClusterWorker::clust_algo;
  using DistClusterWorker::threads;
  using DistClusterWorker::mutex;
  using DistClusterWorker::_prealigned;
  using DistClusterWorker::_aligned;

  using Ksw2ClusterWorker::match;
  using Ksw2ClusterWorker::mismatch;
  using Ksw2ClusterWorker::gap_open;
  using Ksw2ClusterWorker::gap_extend;
  using Ksw2ClusterWorker::gap_open2;
  using Ksw2ClusterWorker::gap_extend2;
  public:
  Ksw2ConcurrentClusterWorker(
    const SequenceSet &seq,
    ClusterAlgorithm &clust_algo,
    DivisiblePairGenerator::Builder & pgb,
    const int match = 0, const int mismatch = 1,
    const int gap_open = 0, const int gap_extend = 1,
    const int gap_open2 = 0, const int gap_extend2 = 1,
    std::size_t worker_threads = 0
  ) : Ksw2ClusterWorker(seq, clust_algo, pgb, match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2, verbose, worker_threads) {};
  void operator()(std::size_t begin, std::size_t end) override;
};

#endif
