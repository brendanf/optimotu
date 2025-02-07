#include "ClusterWorker.h"

ClusterWorker::ClusterWorker(std::istream &file, const int threads) :
  file(file), threads(threads) {}

int ClusterWorker::n_threads() {
  return this->threads;
}

template<>
MergeClusterWorker<int>::MergeClusterWorker(
  ClusterAlgorithm *algo,
  std::istream &file,
  const int threads
) :
  ClusterWorker(file, threads),
  algo_list(threads) {
  // OPTIMOTU_COUT << "MergeClusterWorker constructor start...";
  for (int i = 0; i < threads; ++i) {
    algo_list[i] = algo->make_child();
  }
  // OPTIMOTU_COUT << "done" << std::endl;
}

template<>
MergeClusterWorker<std::string>::MergeClusterWorker(
  ClusterAlgorithm *algo,
  std::istream &file,
  const std::vector<std::string> & id_list,
  const int threads
) :
  ClusterWorker(file, threads),
  algo_list(threads) {
  // OPTIMOTU_COUT << "MergeClusterWorker constructor start...";
  id_map = std::make_shared<id_map_type>();
  for (int i = 0; i < id_list.size(); ++i) {
    id_map->insert({id_list[i], i});
  }
  for (int i = 0; i < threads; ++i) {
    algo_list[i] = algo->make_child();
  }
  // OPTIMOTU_COUT << "done" << std::endl;
}

template<typename id_type>
void MergeClusterWorker<id_type>::operator()(size_t begin, size_t end) {
  DistanceElement d;
  // std::vector<DistanceElement> buffer;
  // buffer.reserve(100);
  // OPTIMOTU_COUT << "Starting MergeClusterWorker thread " << begin << std::endl;
  while (true) {
    {
      std::lock_guard<std::mutex> lock(mutex);
    // for (int i = 0; i < 100 && file; ++i) {
    if constexpr (std::is_same<id_list_type, void>::value) {
      if (!(file >> d)) break;
    } else {
      id_type id1, id2;
      double dist;
      if (!(file >> id1 >> id2 >> d)) break;
      auto it1 = id_map->find(id1);
      auto it2 = id_map->find(id2);
      if (it1 == id_map->end() || it2 == id_map->end()) {
        // silently skip the line if either id is not found!
        continue;
      }
      d.seq1 = it1->second;
      d.seq2 = it2->second;
      d.dist = dist;
    }
      // buffer.push_back(d);
    // }
    }
    // if (buffer.size() == 0) break;
    // for (auto d : buffer) {
    algo_list[begin]->operator()(d, begin);
    // }
    // buffer.clear();
  }
  algo_list[begin]->finalize();
  // mutex.lock();
  // OPTIMOTU_COUT << "ClusterWorker thread " << begin << " merging..." << std::endl;
  // mutex.unlock();
  algo_list[begin]->merge_into_parent();

  // mutex.lock();
  // OPTIMOTU_COUT << "ClusterWorker thread " << begin << " exiting" << std::endl;
  // mutex.unlock();
}

template<typename id_type>
void MergeClusterWorker<id_type>::finalize() {}

template<>
ConcurrentClusterWorker<int>::ConcurrentClusterWorker(
  ClusterAlgorithm *algo,
  std::istream &file,
  const int threads
) :
  ClusterWorker(file, threads),
  algo(algo) {}

template<>
ConcurrentClusterWorker<std::string>::ConcurrentClusterWorker(
  ClusterAlgorithm *algo,
  std::istream &file,
  const id_list_type & id_list,
  const int threads
) :
  ClusterWorker(file, threads),
  algo(algo) {
  id_map = std::make_shared<id_map_type>();
  for (int i = 0; i < id_list.size(); ++i) {
    id_map->insert({id_list[i], i});
  }
}

template<typename id_type>
void ConcurrentClusterWorker<id_type>::operator()(size_t begin, size_t end) {
  DistanceElement d;
  // std::vector<DistanceElement> buffer;
  // buffer.reserve(100);
  // OPTIMOTU_COUT << "Starting ConcurrentClusterWorker thread " << begin << std::endl;
  while (true) {
    {
      std::lock_guard<std::mutex> lock(mutex);
    // OPTIMOTU_COUT << "reading line " << ++i << "..." << std::flush;
    // for (int i = 0; i < 100 && file; ++i
    if constexpr (std::is_same<id_list_type, void>::value) {
      if (!(file >> d)) break;
    } else {
      id_type id1, id2;
      double dist;
      if (!(file >> id1 >> id2 >> d)) break;
      auto it1 = id_map->find(id1);
      auto it2 = id_map->find(id2);
      if (it1 == id_map->end() || it2 == id_map->end()) {
        // silently skip the line if either id is not found!
        continue;
      }
      d.seq1 = it1->second;
      d.seq2 = it2->second;
      d.dist = dist;
    }
      // buffer.push_back(d);
    }
    // OPTIMOTU_COUT << "done" << std::endl;
    // if (buffer.size() == 0) break;
    // for (auto d : buffer) {
      algo->operator()(d, begin);
    // }
    // buffer.clear();
  }

  // mutex.lock();
  // OPTIMOTU_COUT << "ConcurrentClusterWorker thread " << begin << " exiting" << std::endl;
  // mutex.unlock();
}

template<typename id_type>
void ConcurrentClusterWorker<id_type>::finalize() {}

template<>
HierarchicalClusterWorker<int>::HierarchicalClusterWorker(
  ClusterAlgorithm *algo,
  std::istream &file,
  const int threads,
  const int shards
) : ClusterWorker(file, threads), algo_list(shards),
thread_count(shards), shards(shards) {
  // OPTIMOTU_COUT << "HieararchicalClusterWorker constructor start...";
  for (int i = 0; i < shards; ++i) {
    algo_list[i] = algo->make_child();
  }
  // OPTIMOTU_COUT << "done" << std::endl;
}

template<>
HierarchicalClusterWorker<std::string>::HierarchicalClusterWorker(
  ClusterAlgorithm *algo,
  std::istream &file,
  const id_list_type & id_list,
  const int threads,
  const int shards
) : ClusterWorker(file, threads), algo_list(shards),
thread_count(shards), shards(shards) {
  // OPTIMOTU_COUT << "HieararchicalClusterWorker constructor start...";
  id_map = std::make_shared<id_map_type>();
  for (int i = 0; i < id_list.size(); ++i) {
    id_map->insert({id_list[i], i});
  }
  for (int i = 0; i < shards; ++i) {
    algo_list[i] = algo->make_child();
  }
  // OPTIMOTU_COUT << "done" << std::endl;
}

template<typename id_type>
void HierarchicalClusterWorker<id_type>::operator()(size_t begin, size_t end) {
  DistanceElement d;
  // std::vector<DistanceElement> buffer;
  // buffer.reserve(100);
  size_t i = begin % shards;
  // mutex.lock();
  // size_t tc =
    ++thread_count[i];
  // OPTIMOTU_COUT << "Starting HierarchicalClusterWorker thread " << begin
  //           << " (shard " << i
  //           << " with " << tc
  //           << " previous threads)"
  //           << std::endl;
  // mutex.unlock();
  while (true) {
    {
      std::lock_guard<std::mutex> lock(mutex);
    // for (int i = 0; i < 100 && file; ++i) {
    if constexpr (std::is_same<id_list_type, void>::value) {
      if (!(file >> d)) break;
    } else {
      id_type id1, id2;
      double dist;
      if (!(file >> id1 >> id2 >> d)) break;
      auto it1 = id_map->find(id1);
      auto it2 = id_map->find(id2);
      if (it1 == id_map->end() || it2 == id_map->end()) {
        // silently skip the line if either id is not found!
        continue;
      }
      d.seq1 = it1->second;
      d.seq2 = it2->second;
      d.dist = dist;
    }
    //   buffer.push_back(d);
    }

    // if (buffer.size() == 0) break;
    // for (auto d : buffer) {
      algo_list[i]->operator()(d, begin);
    // }
    // buffer.clear();
  }
  if (--thread_count[i] == 0) {
    // mutex.lock();
    // OPTIMOTU_COUT << "HierarchicalClusterWorker thread " << begin
    //           << " merging shard " << i
    //           << "..." << std::endl;
    // mutex.unlock();
    algo_list[i]->finalize();
    algo_list[i]->merge_into_parent();
  }

    // mutex.lock();
    // OPTIMOTU_COUT << "HieararchicalClusterWorker thread " << begin
    //           << " exiting" << std::endl;
    // mutex.unlock();
}

template<typename id_type>
void HierarchicalClusterWorker<id_type>::finalize() {}
