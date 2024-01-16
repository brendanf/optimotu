#ifndef OPTIMOTU_HAMMINGCLUSTERWORKER_H_INCLUDED
#define OPTIMOTU_HAMMINGCLUSTERWORKER_H_INCLUDED

#include <Rcpp.h>
#include <RcppParallel.h>
#include <RcppThread.h>
#include "ClusterAlgorithm.h"
extern "C" {
#include "defs.h"
}

// Analogous to SequenceSetB in C code
struct PackedSequenceSet {
  int num_seqs, alen, ulen, mulen;
  std::vector<std::vector<unsigned long int>> packed_seq;
  std::vector<std::vector<unsigned long int>> mask;
  PackedSequenceSet(const std::vector<std::string> &seq);

  double dist(const int i, const int j) const;
};

class HammingClusterWorker : public RcppParallel::Worker {
protected:
  const PackedSequenceSet pss;
  ClusterAlgorithm &clust_algo;
  const uint8_t threads;
  std::mutex mutex;
  size_t _prealigned = 0, _aligned = 0;
  const bool verbose;
public :
  HammingClusterWorker(
    const std::vector<std::string> &seq,
    ClusterAlgorithm &clust_algo,
    const uint8_t threads,
    bool verbose = false
  );

  size_t prealigned() {
    return _prealigned;
  }

  size_t aligned() {
    return _aligned;
  }
};

class HammingSplitClusterWorker : public HammingClusterWorker {

public :
  HammingSplitClusterWorker(
    const std::vector<std::string> &seq,
    ClusterAlgorithm &clust_algo,
    const uint8_t threads,
    bool verbose = false
  );
  void operator()(std::size_t begin, std::size_t end);
};

class HammingConcurrentClusterWorker : public HammingClusterWorker {
public :
  HammingConcurrentClusterWorker(
    const std::vector<std::string> &seq,
    ClusterAlgorithm &clust_algo,
    const uint8_t threads,
    bool verbose = false
  );
  void operator()(std::size_t begin, std::size_t end);
};

#endif //OPTIMOTU_HAMMINGCLUSTERWORKER_H_INCLUDED
