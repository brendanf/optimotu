#ifndef OPTIMOTU_WFA2CLUSTERWORKER_H
#define OPTIMOTU_WFA2CLUSTERWORKER_H

#include "DistClusterWorker.h"

class Wfa2ClusterWorker : public DistClusterWorker {
protected:
  int match = 0, mismatch = 1,
    gap_open = 0, gap_extend = 1,
    gap_open2 = 0, gap_extend2 = 1;
public :
  Wfa2ClusterWorker(
    const std::vector<std::string> &seq,
    ClusterAlgorithm &clust_algo,
    DivisiblePairGenerator::Builder & pgb,
    const int match = 0, const int mismatch = 1,
    const int gap_open = 0, const int gap_extend = 1,
    const int gap_open2 = 0, const int gap_extend2 = 1,
    int verbose = 0
  );
};

template<int verbose>
class Wfa2SplitClusterWorker : public Wfa2ClusterWorker {
  using DistClusterWorker::seq;
  using DistClusterWorker::clust_algo;
  using DistClusterWorker::threads;
  using DistClusterWorker::mutex;
  using DistClusterWorker::_prealigned;
  using DistClusterWorker::_aligned;

  using Wfa2ClusterWorker::match;
  using Wfa2ClusterWorker::mismatch;
  using Wfa2ClusterWorker::gap_open;
  using Wfa2ClusterWorker::gap_extend;
  using Wfa2ClusterWorker::gap_open2;
  using Wfa2ClusterWorker::gap_extend2;
public:
  Wfa2SplitClusterWorker(
    const std::vector<std::string> &seq,
    ClusterAlgorithm &clust_algo,
    DivisiblePairGenerator::Builder & pgb,
    const int match = 0, const int mismatch = 1,
    const int gap_open = 0, const int gap_extend = 1,
    const int gap_open2 = 0, const int gap_extend2 = 1
  ) : Wfa2ClusterWorker(seq, clust_algo, pgb, match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2, verbose) {};
  void operator()(std::size_t begin, std::size_t end) override;
};

template<int verbose>
class Wfa2ConcurrentClusterWorker : public Wfa2ClusterWorker {
  using DistClusterWorker::seq;
  using DistClusterWorker::clust_algo;
  using DistClusterWorker::threads;
  using DistClusterWorker::mutex;
  using DistClusterWorker::_prealigned;
  using DistClusterWorker::_aligned;

  using Wfa2ClusterWorker::match;
  using Wfa2ClusterWorker::mismatch;
  using Wfa2ClusterWorker::gap_open;
  using Wfa2ClusterWorker::gap_extend;
  using Wfa2ClusterWorker::gap_open2;
  using Wfa2ClusterWorker::gap_extend2;
  public:
  Wfa2ConcurrentClusterWorker(
    const std::vector<std::string> &seq,
    ClusterAlgorithm &clust_algo,
    DivisiblePairGenerator::Builder & pgb,
    const int match = 0, const int mismatch = 1,
    const int gap_open = 0, const int gap_extend = 1,
    const int gap_open2 = 0, const int gap_extend2 = 1
  ) : Wfa2ClusterWorker(seq, clust_algo, pgb, match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2, verbose) {};
  void operator()(std::size_t begin, std::size_t end) override;
};

#endif
