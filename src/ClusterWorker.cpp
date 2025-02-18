#include "ClusterWorker.h"

ClusterWorker::ClusterWorker(ClusterAlgorithm * algo, const int threads) :
  threads(threads), algo(algo) {}

int ClusterWorker::n_threads() {
  return this->threads;
}

template <typename distmx_t, typename id_t>
ClusterWorkerImpl<distmx_t, id_t>::ClusterWorkerImpl(
    ClusterAlgorithm * algo,
    distmx_t &distmx,
    const int threads
) :
  ClusterWorker(algo, threads), distmx(distmx), translator(distmx) {
  if constexpr (!std::is_same<id_t, int>::value) {
    OPTIMOTU_STOP("Distance matrix with non-integer id requires id list. ",
                  "Please report this as a bug!");
  }
  if constexpr (std::is_same<distmx_t, Rcpp::DataFrame>::value) {
    if (distmx.ncol() != 3) {
      OPTIMOTU_STOP("Distance matrix must have 3 columns: id1, id2, distance");
    }
    if (!Rcpp::is<Rcpp::IntegerVector>(distmx[0]) ||
        !Rcpp::is<Rcpp::IntegerVector>(distmx[1])) {
        OPTIMOTU_STOP("Distance matrix with integer id must have integer id columns");
    }
    if (!Rcpp::is<Rcpp::NumericVector>(distmx[2])) {
      OPTIMOTU_STOP("Distance matrix must have numeric distance column");
    }
  } else if constexpr (std::is_same<distmx_t, std::istream>::value) {
    if (!distmx.good()) {
      OPTIMOTU_STOP("Failed to read distance matrix");
    }
  } else {
    OPTIMOTU_STOP("Unsupported distmx_t type");
  }
}

template <typename distmx_type, typename id_t>
ClusterWorkerImpl<distmx_type, id_t>::ClusterWorkerImpl(
    ClusterAlgorithm * algo,
    distmx_type &distmx,
    const id_list_type &id_list,
    const int threads
) :
  ClusterWorker(algo, threads), distmx(distmx), translator(distmx) {
  if constexpr (std::is_same<id_t, int>::value) {
    OPTIMOTU_STOP("Distance matrix with integer id requires no id list. ",
                  "Please report this as a bug!");
  }
  if constexpr (std::is_same<distmx_type, Rcpp::DataFrame>::value) {
    if (distmx.ncol() != 3) {
      OPTIMOTU_STOP("Distance matrix must have 3 columns: id1, id2, distance");
    }
    if (!Rcpp::is<Rcpp::CharacterVector>(distmx[0]) ||
        !Rcpp::is<Rcpp::CharacterVector>(distmx[1])) {
        OPTIMOTU_STOP("Distance matrix with name id must have string id columns");
    }
    if (!Rcpp::is<Rcpp::NumericVector>(distmx[2])) {
      OPTIMOTU_STOP("Distance matrix must have numeric distance column");
    }
  } else if constexpr (std::is_same<distmx_type, std::istream>::value) {
    if (!distmx.good()) {
      OPTIMOTU_STOP("Failed to read distance matrix");
    }
  } else {
    OPTIMOTU_STOP("Unsupported distmx_t type");
  }
  id_map = std::make_shared<id_map_type>();
  for (int i = 0; i < id_list.size(); ++i) {
    id_map->insert({id_list[i], i});
  }
}

template <class distmx_type, class id_t>
bool ClusterWorkerImpl<distmx_type, id_t>::next_line(DistanceElement & d) {
  std::lock_guard<std::mutex> lock(mutex);
  if constexpr (std::is_same<distmx_type, std::istream>::value) {
    if constexpr (std::is_same<id_t, int>::value) {
      if (!(distmx >> d)) return false;
    } else {
      bool found = false;
      do {
        id_t id1, id2;
        double dist;
        if (!(distmx >> id1 >> id2 >> dist)) return false;
        auto it1 = id_map->find(id1);
        auto it2 = id_map->find(id2);
        // silently skip the line if either id is not found!
        if (it1 == id_map->end() || it2 == id_map->end()) {
          ++line_number;
          continue;
        }
        d.seq1 = it1->second;
        d.seq2 = it2->second;
        d.dist = dist;
        found = true;
      } while (!found);
    }
  } else if constexpr (std::is_same<distmx_type, Rcpp::DataFrame>::value) {
    if constexpr (std::is_same<id_t, int>::value) {
      if (line_number >= distmx.nrows()) return false;
      d.seq1 = translator.id1(line_number);
      d.seq2 = translator.id2(line_number);
      d.dist = translator.dist(line_number);
    } else {
      bool found = false;
      do {
        if (line_number >= distmx.nrows()) return false;
        auto it1 = id_map->find(translator.id1(line_number));
        auto it2 = id_map->find(translator.id2(line_number));
        if (it1 == id_map->end() || it2 == id_map->end()) {
          // silently skip the line if either id is not found!
          ++line_number;
          continue;
        }
        d.seq1 = it1->second;
        d.seq2 = it2->second;
        d.dist = translator.dist(line_number);
        found = true;
      } while (!found);
    }
  } else {
    OPTIMOTU_STOP("Unsupported distmx_t type");
  }
  ++line_number;
  return true;
}

template<typename distmx_t, typename id_t>
MergeClusterWorker<distmx_t, id_t>::MergeClusterWorker(
    ClusterAlgorithm *algo,
    distmx_t &distmx,
    const int threads
) : ClusterWorkerImpl<distmx_t, id_t>(algo, distmx, threads) {
  algo_list.reserve(threads);
  // OPTIMOTU_COUT << "MergeClusterWorker constructor start...";
  for (int i = 0; i < threads; ++i) {
    algo_list.push_back(algo->make_child());
  }
  // OPTIMOTU_COUT << "done" << std::endl;
}

template<typename distmx_t, typename id_t>
MergeClusterWorker<distmx_t, id_t>::MergeClusterWorker(
    ClusterAlgorithm *algo,
    distmx_t &distmx,
    const id_list_type &id_list,
    const int threads
) : ClusterWorkerImpl<distmx_t, id_t>(algo, distmx, id_list, threads) {
  algo_list.reserve(threads);
  for (int i = 0; i < threads; ++i) {
    algo_list.push_back(algo->make_child());
  }
}

template<typename distmx_t, typename id_type>
void MergeClusterWorker<distmx_t, id_type>::operator()(size_t begin, size_t end) {
  DistanceElement d;
  // OPTIMOTU_COUT << "Starting MergeClusterWorker thread " << begin << std::endl;
  while (next_line(d)) {
    algo_list[begin]->operator()(d, begin);
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

template<typename distmx_t, typename id_type>
void ConcurrentClusterWorker<distmx_t, id_type>::operator()(size_t begin, size_t end) {
  DistanceElement d;
  // OPTIMOTU_COUT << "Starting ConcurrentClusterWorker thread " << begin << std::endl;
  while (next_line(d)) {
    algo->operator()(d, begin);
  }

  // mutex.lock();
  // OPTIMOTU_COUT << "ConcurrentClusterWorker thread " << begin << " exiting" << std::endl;
  // mutex.unlock();
}

template<typename distmx_t, typename id_t>
HierarchicalClusterWorker<distmx_t, id_t>::HierarchicalClusterWorker(
  ClusterAlgorithm *algo,
  distmx_t &distmx,
  const int threads,
  const int shards
) : ClusterWorkerImpl<distmx_t, id_t>(algo, distmx, threads),
thread_count(shards), shards(shards) {
  // OPTIMOTU_COUT << "HieararchicalClusterWorker constructor start...";
  algo_list.reserve(shards);
  for (int i = 0; i < shards; ++i) {
    thread_count[i] = 0;
    algo_list.push_back(algo->make_child());
  }
  // OPTIMOTU_COUT << "done" << std::endl;
}

template<typename distmx_t, typename id_t>
HierarchicalClusterWorker<distmx_t, id_t>::HierarchicalClusterWorker(
    ClusterAlgorithm *algo,
    distmx_t &distmx,
    const id_list_type &id_list,
    const int threads,
    const int shards
) : ClusterWorkerImpl<distmx_t, id_t>(algo, distmx, id_list, threads),
thread_count(shards), shards(shards) {
  // OPTIMOTU_COUT << "HieararchicalClusterWorker constructor start...";
  algo_list.reserve(shards);
  for (int i = 0; i < shards; ++i) {
    thread_count[i] = 0;
    algo_list[i] = algo->make_child();
  }
  // OPTIMOTU_COUT << "done" << std::endl;
}

template<typename distmx_t, typename id_type>
void HierarchicalClusterWorker<distmx_t, id_type>::operator()(size_t begin, size_t end) {
  DistanceElement d;
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
  while (next_line(d)) {
    algo_list[i]->operator()(d, begin);
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

template class ClusterWorkerImpl<Rcpp::DataFrame, int>;
template class ClusterWorkerImpl<std::istream, int>;
template class ClusterWorkerImpl<Rcpp::DataFrame, std::string>;
template class ClusterWorkerImpl<std::istream, std::string>;

template class MergeClusterWorker<Rcpp::DataFrame, int>;
template class MergeClusterWorker<std::istream, int>;
template class MergeClusterWorker<Rcpp::DataFrame, std::string>;
template class MergeClusterWorker<std::istream, std::string>;

template class ConcurrentClusterWorker<Rcpp::DataFrame, int>;
template class ConcurrentClusterWorker<std::istream, int>;
template class ConcurrentClusterWorker<Rcpp::DataFrame, std::string>;
template class ConcurrentClusterWorker<std::istream, std::string>;

template class HierarchicalClusterWorker<Rcpp::DataFrame, int>;
template class HierarchicalClusterWorker<std::istream, int>;
template class HierarchicalClusterWorker<Rcpp::DataFrame, std::string>;
template class HierarchicalClusterWorker<std::istream, std::string>;
