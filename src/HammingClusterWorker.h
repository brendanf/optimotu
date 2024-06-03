#ifndef OPTIMOTU_HAMMINGCLUSTERWORKER_H_INCLUDED
#define OPTIMOTU_HAMMINGCLUSTERWORKER_H_INCLUDED

#include "AlignClusterWorker.h"


// Analogous to SequenceSetB in C code
struct PackedSequenceSet {
  int num_seqs, alen, ulen, mulen;
  std::vector<std::vector<unsigned long int>> packed_seq;
  std::vector<std::vector<unsigned long int>> mask;
  std::vector<int> start;
  std::vector<int> end;
  PackedSequenceSet(const std::vector<std::string> &seq);

  double dist(const int i, const int j, const int min_overlap, const bool ignore_gap) const;
};

class HammingClusterWorker : public AlignClusterWorker {
protected:
  const PackedSequenceSet pss;
  const int min_overlap;
  const bool ignore_gaps;
public :
  HammingClusterWorker(
    const std::vector<std::string> &seq,
    ClusterAlgorithm &clust_algo,
    const uint8_t threads,
    const int min_overlap = 0,
    const bool ignore_gaps = TRUE,
    bool verbose = false
  );
};

class HammingSplitClusterWorker : public HammingClusterWorker {

public :
  HammingSplitClusterWorker(
    const std::vector<std::string> &seq,
    ClusterAlgorithm &clust_algo,
    const uint8_t threads,
    const int min_overlap = 0,
    const bool ignore_gaps = TRUE,
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
    const int min_overlap = 0,
    const bool ignore_gaps = TRUE,
    bool verbose = false
  );
  void operator()(std::size_t begin, std::size_t end);
};

#endif //OPTIMOTU_HAMMINGCLUSTERWORKER_H_INCLUDED
