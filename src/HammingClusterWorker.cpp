#include "HammingClusterWorker.h"
extern "C" {
#include "defs.h"
}

typedef RcppParallel::RMatrix<int> matrix_t;

HammingClusterWorker::HammingClusterWorker(
  const std::vector<std::string> &seq,
  ClusterAlgorithm &clust_algo,
  DivisiblePairGenerator::Builder & pgb,
  const int min_overlap,
  const bool ignore_gaps,
  int verbose,
  std::size_t worker_threads
) : DistClusterWorker(seq, clust_algo, pgb, verbose, worker_threads), pss(seq),
min_overlap(min_overlap), ignore_gaps(ignore_gaps) {};

template<int verbose>
HammingSplitClusterWorker<verbose>::HammingSplitClusterWorker(
  const std::vector<std::string> &seq,
  ClusterAlgorithm &clust_algo,
  DivisiblePairGenerator::Builder & pgb,
  const int min_overlap,
  const bool ignore_gaps,
  std::size_t worker_threads
) : HammingClusterWorker(seq, clust_algo, pgb, min_overlap, ignore_gaps,
   verbose, worker_threads) {};

template<int verbose>
void HammingSplitClusterWorker<verbose>::operator()(std::size_t begin, std::size_t end) {
  size_t my_prealigned = 0;
  size_t my_aligned = 0;

  for (size_t pg_index = begin; pg_index < pair_generators.size(); pg_index += threads) {
    auto & pg = pair_generators[pg_index];
    ClusterAlgorithm * my_algo = clust_algo.make_child();
    OPTIMOTU_DEBUG(
      2,
      << "HammingSplitClusterWorker thread " << pg_index
      << " entered" << std::endl
    );
    while (*pg) {
      size_t i = pg->i();
      size_t j = pg->j();
      size_t i0 = pg->i0();
      size_t j0 = pg->j0();
      double threshold = my_algo->max_relevant(*pg);
      OPTIMOTU_DEBUG(
        4,
        << "thread" << pg_index
        << ": seqs " << j << " (j0=" << j0 << ")"
        << " and " << i << " (i0=" << i0 << ")"
        << " max relevant=" << threshold
        << std::endl
      );
      ++my_prealigned;
      double d = pss.dist(i0, j0, min_overlap, ignore_gaps);
      if (d < 1.0) ++my_aligned;

      OPTIMOTU_DEBUG(
        4,
        << (d <= threshold ? "*" : " ")
        << " distance=" << d
        << std::endl;
      );
      if (d <= threshold) (*my_algo)(*pg, d);
      ++(*pg);
      RcppThread::checkUserInterrupt();
    }
    mutex.lock();
    OPTIMOTU_DEBUG(2, << "thread " << pg_index << " ready to merge" << std::endl);
    _aligned += my_aligned;
    _prealigned += my_prealigned;
    mutex.unlock();
    my_algo->merge_into_parent();
    clust_algo.release_child(my_algo);
    OPTIMOTU_DEBUG(2, << "thread " << pg_index << " done" << std::endl);
  }
}

template<int verbose>
HammingConcurrentClusterWorker<verbose>::HammingConcurrentClusterWorker(
  const std::vector<std::string> &seq,
  ClusterAlgorithm &clust_algo,
  DivisiblePairGenerator::Builder & pgb,
  const int min_overlap,
  const bool ignore_gaps,
  std::size_t worker_threads
) : HammingClusterWorker(seq, clust_algo, pgb, min_overlap, ignore_gaps,
   verbose, worker_threads) {};

template<int verbose>
void HammingConcurrentClusterWorker<verbose>::operator()(std::size_t begin, std::size_t end) {
  size_t my_prealigned = 0;
  size_t my_aligned = 0;

  for (size_t pg_index = begin; pg_index < pair_generators.size(); pg_index += threads) {
    auto & pg = pair_generators[pg_index];
    OPTIMOTU_DEBUG(
      2,
      << "HammingConcurrentClusterWorker thread " << pg_index
      << " entered" << std::endl
    );
    while (*pg) {
      size_t i = pg->i();
      size_t j = pg->j();
      size_t i0 = pg->i0();
      size_t j0 = pg->j0();
      OPTIMOTU_DEBUG(
        4,
        << "thread" << pg_index
        << ": seqs " << j << " (j0=" << j0 << ")"
        << " and " << i << " (i0=" << i0 << ")"
        << std::endl
      );
      double threshold = clust_algo.max_relevant(*pg);
      ++my_prealigned;
      double d = pss.dist(i0, j0, min_overlap, ignore_gaps);
      if (d < 1.0) ++my_aligned;
      if (d < threshold) clust_algo(*pg, d);
      ++(*pg);
      RcppThread::checkUserInterrupt();
    }
    mutex.lock();
    _aligned += my_aligned;
    _prealigned += my_prealigned;
    OPTIMOTU_DEBUG(2, << "thread " << pg_index << " done" << std::endl);
    mutex.unlock();
  }
}

template class HammingSplitClusterWorker<0>;
template class HammingSplitClusterWorker<1>;
template class HammingSplitClusterWorker<2>;
template class HammingSplitClusterWorker<3>;
template class HammingSplitClusterWorker<4>;
template class HammingConcurrentClusterWorker<0>;
template class HammingConcurrentClusterWorker<1>;
template class HammingConcurrentClusterWorker<2>;
template class HammingConcurrentClusterWorker<3>;
template class HammingConcurrentClusterWorker<4>;