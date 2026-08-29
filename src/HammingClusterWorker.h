#ifndef OPTIMOTU_HAMMINGCLUSTERWORKER_H_INCLUDED
#define OPTIMOTU_HAMMINGCLUSTERWORKER_H_INCLUDED

#include "DistClusterWorker.h"
#include "PackedSequenceSet.h"
#include "SequenceView.h"
#include <cstddef>

class HammingClusterWorker : public DistClusterWorker {
protected:
  const PackedSequenceSet pss;
  const int min_overlap;
  const bool ignore_gaps;
  std::size_t tracked_pss_bytes = 0;
public :
  HammingClusterWorker(
    const SequenceSet &seq,
    ClusterAlgorithm &clust_algo,
    DivisiblePairGenerator::Builder & pgb,
    const int min_overlap = 0,
    const bool ignore_gaps = true,
    int verbose = 0,
    std::size_t worker_threads = 0
  );
  ~HammingClusterWorker() override;
};

template<int verbose>
class HammingSplitClusterWorker : public HammingClusterWorker {
  using HammingClusterWorker::pss;
  using DistClusterWorker::clust_algo;
  using DistClusterWorker::threads;
  using DistClusterWorker::pair_generators;
  using DistClusterWorker::mutex;
  using HammingClusterWorker::min_overlap;
  using HammingClusterWorker::ignore_gaps;
  using DistClusterWorker::_prealigned;
  using DistClusterWorker::_aligned;

public :
  HammingSplitClusterWorker(
    const SequenceSet &seq,
    ClusterAlgorithm &clust_algo,
    DivisiblePairGenerator::Builder & pgb,
    const int min_overlap = 0,
    const bool ignore_gaps = TRUE,
    std::size_t worker_threads = 0
  );
  void operator()(std::size_t begin, std::size_t end);
};

template<int verbose>
class HammingConcurrentClusterWorker : public HammingClusterWorker {
  using HammingClusterWorker::pss;
  using DistClusterWorker::clust_algo;
  using DistClusterWorker::threads;
  using DistClusterWorker::pair_generators;
  using DistClusterWorker::mutex;
  using HammingClusterWorker::min_overlap;
  using HammingClusterWorker::ignore_gaps;
  using DistClusterWorker::_prealigned;
  using DistClusterWorker::_aligned;
public :
  HammingConcurrentClusterWorker(
    const SequenceSet &seq,
    ClusterAlgorithm &clust_algo,
    DivisiblePairGenerator::Builder & pgb,
    const int min_overlap = 0,
    const bool ignore_gaps = true,
    std::size_t worker_threads = 0
  );
  void operator()(std::size_t begin, std::size_t end) override;
};

#endif //OPTIMOTU_HAMMINGCLUSTERWORKER_H_INCLUDED
