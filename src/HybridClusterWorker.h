#include "AlignClusterWorker.h"

class HybridClusterWorker : public AlignClusterWorker {
protected:
  double breakpoint = 0.1;
public :
  HybridClusterWorker(
    const std::vector<std::string> &seq,
    ClusterAlgorithm &clust_algo,
    const uint8_t threads,
    const double breakpoint = 0.1,
    bool verbose = false
  );
};

class HybridSplitClusterWorker : public HybridClusterWorker {
public :
  using HybridClusterWorker::HybridClusterWorker;
  void operator()(std::size_t begin, std::size_t end);
};

class HybridConcurrentClusterWorker : public HybridClusterWorker {
public :
  using HybridClusterWorker::HybridClusterWorker;
  void operator()(std::size_t begin, std::size_t end);
};

