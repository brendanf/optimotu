// SPDX-FileCopyrightText: 2025 Brendan Furneaux <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_CLUSTERWORKER_H_INCLUDED
#define OPTIMOTU_CLUSTERWORKER_H_INCLUDED

#include <atomic>
#include <RcppParallel.h>
#include "ClusterAlgorithm.h"
#include "MultipleClusterAlgorithm.h"
#include <unordered_map>
#include <memory>

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

template <class id_type>
class MergeClusterWorker : public ClusterWorker {
protected:
  std::vector<ClusterAlgorithm*> algo_list;
  std::mutex mutex;
  typedef std::unordered_map<id_type, int> id_map_type;
  std::shared_ptr<id_map_type> id_map;
public:
  typedef std::vector<id_type> id_list_type;

  MergeClusterWorker(ClusterAlgorithm *algo, std::istream &file,
                                const int threads);

  MergeClusterWorker(ClusterAlgorithm *algo, std::istream &file,
                     const id_list_type & id_list,
                     const int threads);

  void operator()(size_t begin, size_t end) override;

  void finalize() override;
};

template <> MergeClusterWorker<std::string>::MergeClusterWorker(
    ClusterAlgorithm *algo,
    std::istream &file,
    const int threads
) = delete;

template <> MergeClusterWorker<int>::MergeClusterWorker(
    ClusterAlgorithm *algo,
    std::istream &file,
    const id_list_type & id_list,
    const int threads
) = delete;

template <class id_type>
class ConcurrentClusterWorker : public ClusterWorker {
protected:
  ClusterAlgorithm* algo;
  std::mutex mutex;
  typedef std::unordered_map<id_type, int> id_map_type;
  std::shared_ptr<id_map_type> id_map;
public:
  typedef std::vector<id_type> id_list_type;

  ConcurrentClusterWorker(
      ClusterAlgorithm *algo,
      std::istream &file,
      const int threads
  );

  ConcurrentClusterWorker(
      ClusterAlgorithm *algo,
      std::istream &file,
      const id_list_type & id_list,
      const int threads
  );

  void operator()(size_t begin, size_t end) override;

  void finalize() override;
};

template <> ConcurrentClusterWorker<std::string>::ConcurrentClusterWorker(
    ClusterAlgorithm *algo,
    std::istream &file,
    const int threads
) = delete;

template <> ConcurrentClusterWorker<int>::ConcurrentClusterWorker(
    ClusterAlgorithm *algo,
    std::istream &file,
    const id_list_type & id_list,
    const int threads
) = delete;

template <class id_type>
class HierarchicalClusterWorker : public ClusterWorker {
protected:
  std::vector<ClusterAlgorithm*> algo_list;
  std::mutex mutex;
  std::vector<std::atomic_size_t> thread_count;
  const int shards;
  typedef std::unordered_map<id_type, int> id_map_type;
  std::shared_ptr<id_map_type> id_map;
public:
  typedef std::vector<id_type> id_list_type;

  HierarchicalClusterWorker(
    ClusterAlgorithm *algo,
    std::istream &file,
    const int threads,
    const int shards
  );

  HierarchicalClusterWorker(
    ClusterAlgorithm *algo,
    std::istream &file,
    const id_list_type & id_list,
    const int threads,
    const int shards
  );

  void operator()(size_t begin, size_t end) override;

  void finalize() override;
};

template <> HierarchicalClusterWorker<std::string>::HierarchicalClusterWorker(
    ClusterAlgorithm *algo,
    std::istream &file,
    const int threads,
    const int shards
) = delete;

template <> HierarchicalClusterWorker<int>::HierarchicalClusterWorker(
    ClusterAlgorithm *algo,
    std::istream &file,
    const id_list_type & id_list,
    const int threads,
    const int shards
) = delete;

#endif //OPTIMOTU_CLUSTERWORKER_H_INCLUDED
