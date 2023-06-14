#ifndef OPTIMOTU_CLUSTERWORKER_H_INCLUDED
#define OPTIMOTU_CLUSTERWORKER_H_INCLUDED

#include <atomic>
#include <RcppParallel.h>
#include "ClusterAlgorithm.h"
#include "MultipleClusterAlgorithm.h"

class ClusterWorker : public RcppParallel::Worker {
protected:
  std::istream &file;
  const int threads;
public:
  ClusterWorker(std::istream &file, const int threads);
  virtual ~ClusterWorker() = default;
  virtual void finalize()=0;
  int n_threads();
};

class MergeClusterWorker : public ClusterWorker {
protected:
  std::vector<ClusterAlgorithm*> algo_list;
  std::mutex mutex;
public:
  MergeClusterWorker(ClusterAlgorithm *algo, std::istream &file, const int threads);

  void operator()(size_t begin, size_t end) override;

  void finalize() override;
};

class ConcurrentClusterWorker : public ClusterWorker {
protected:
  ClusterAlgorithm* algo;
  std::mutex mutex;
public:
  ConcurrentClusterWorker(ClusterAlgorithm *algo, std::istream &file, const int threads);

  void operator()(size_t begin, size_t end) override;

  void finalize() override;
};

class HierarchicalClusterWorker : public ClusterWorker {
protected:
  std::vector<ClusterAlgorithm*> algo_list;
  std::mutex mutex;
  std::vector<std::atomic_size_t> thread_count;
  const int shards;
public:
  HierarchicalClusterWorker(
    ClusterAlgorithm *algo,
    std::istream &file,
    const int threads,
    const int shards
  );

  void operator()(size_t begin, size_t end) override;

  void finalize() override;
};

#endif //OPTIMOTU_CLUSTERWORKER_H_INCLUDED
