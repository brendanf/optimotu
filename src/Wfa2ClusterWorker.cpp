#include "Wfa2ClusterWorker.h"
#include <bindings/cpp/WFAligner.hpp>
#include "pairwise_alignment.h"

Wfa2ClusterWorker::Wfa2ClusterWorker(
  const SequenceSet &seq,
  ClusterAlgorithm &clust_algo,
  DivisiblePairGenerator::Builder & pgb,
  const int match, const int mismatch,
  const int gap_open, const int gap_extend,
  const int gap_open2, const int gap_extend2,
  int verbose,
  std::size_t worker_threads
) : DistClusterWorker(seq, clust_algo, pgb, verbose, worker_threads),
match(match), mismatch(mismatch), gap_open(gap_open), gap_extend(gap_extend),
gap_open2(gap_open2), gap_extend2(gap_extend2) {};

template<int verbose>
void Wfa2SplitClusterWorker<verbose>::operator()(std::size_t begin, std::size_t end) {
  wfa::WFAlignerChoose wfa_aligner{match, mismatch, gap_open, gap_extend,
                                   gap_open2, gap_extend2, wfa::WFAligner::Alignment};
  for (size_t pg_index = begin; pg_index < pair_generators.size(); pg_index += threads) {
    auto & pg = pair_generators[pg_index];
    ClusterAlgorithm *my_algo = tile_algo(pg.get());
    if (!my_algo)
    {
      continue;
    }
    OPTIMOTU_DEBUG(2, << "Wfa2SplitClusterWorker thread " << pg_index << " entered" << std::endl);
    size_t my_prealigned = 0;
    size_t my_aligned = 0;
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

      bool is_seqj_longer = seq[j0].size() > seq[i0].size();
      // Use global indices (i0/j0) into seq; tile-local i/j only match
      // when the generator offset is 0 (single full-range tile).
      size_t s1 = is_seqj_longer ? i0 : j0;
      size_t s2 = is_seqj_longer ? j0 : i0;
      double l1 = seq[s1].size(), l2 = seq[s2].size();
      OPTIMOTU_DEBUG(
        4,
        << "#### seq " << i0 << " (l1=" << l1 << ") and "
        << j0 << " (l2=" << l2 <<")####" << std::endl
      );

      double sim_threshold = 1.0 - threshold; // compiler can probably do this?
      if (l1/l2 >= sim_threshold) {
        ++my_prealigned;
        double sim_threshold_plus_1 = 2.0 - threshold;
        double maxd1 = threshold * (l1 + l2) / sim_threshold_plus_1;
        int max_k = (int)ceil((l2 - l1 * sim_threshold) / sim_threshold_plus_1);
        int min_k = -(int)ceil((l1 - l2 * sim_threshold) / sim_threshold_plus_1);
        wfa_aligner.setHeuristicBandedStatic(min_k, max_k);
        wfa_aligner.setMaxAlignmentSteps((int)maxd1 + 1);
        double d = distance_wfa2(seq[s1], seq[s2], wfa_aligner);
        if (d < 1.0) ++my_aligned;

        OPTIMOTU_DEBUG(
          4,
          << " distance=" << d
          << std::endl
        );
        if (d < threshold) (*my_algo)(*pg, d);
      }
      ++(*pg);
      RcppThread::checkUserInterrupt();
    }
    mutex.lock();
    OPTIMOTU_DEBUG(2, << "thread " << pg_index << " ready to merge" << std::endl);
    _aligned += my_aligned;
    _prealigned += my_prealigned;
    mutex.unlock();
    finish_tile(my_algo);
    OPTIMOTU_DEBUG(2, << "thread " << pg_index << " done" << std::endl);
  }
}

template <int verbose>
void Wfa2ConcurrentClusterWorker<verbose>::operator()(std::size_t begin, std::size_t end) {
  size_t my_prealigned = 0;
  size_t my_aligned = 0;

  wfa::WFAlignerChoose wfa_aligner{match, mismatch, gap_open, gap_extend,
                                   gap_open2, gap_extend2, wfa::WFAligner::Alignment};
  for (size_t pg_index = begin; pg_index < pair_generators.size(); pg_index += threads) {
    auto & pg = pair_generators[pg_index];
    OPTIMOTU_DEBUG(
      2,
      << "Wfa2ConcurrentClusterWorker thread " << pg_index
      << " entered" << std::endl
    );
    while (*pg) {
      size_t i = pg->i();
      size_t j = pg->j();
      size_t i0 = pg->i0();
      size_t j0 = pg->j0();
      double threshold = clust_algo.max_relevant(*pg);
      OPTIMOTU_DEBUG(
        4,
        << "thread" << pg_index
        << ": seqs " << j << " (j0=" << j0 << ")"
        << " and " << i << " (i0=" << i0 << ")"
        << " max relevant=" << threshold
        << std::endl
      );
      bool is_seqj_longer = seq[j0].size() > seq[i0].size();
      size_t s1 = is_seqj_longer ? i0 : j0;
      size_t s2 = is_seqj_longer ? j0 : i0;
      double l1 = seq[s1].size(), l2 = seq[s2].size();
      OPTIMOTU_DEBUG(
        4,
        << "#### seq " << s1 << " (l1=" << l1 << ") and "
        << s2 << " (l2=" << l2 <<")####" << std::endl
      );

      double sim_threshold = 1.0 - threshold; // compiler can probably do this?
      if (l1/l2 >= sim_threshold) {
        ++my_prealigned;
        double sim_threshold_plus_1 = 2.0 - threshold;
        double maxd1 = threshold * (l1 + l2) / sim_threshold_plus_1;
        int max_k = (int)ceil((l2 - l1 * sim_threshold) / sim_threshold_plus_1);
        int min_k = -(int)ceil((l1 - l2 * sim_threshold) / sim_threshold_plus_1);
        wfa_aligner.setHeuristicBandedStatic(min_k, max_k);
        wfa_aligner.setMaxAlignmentSteps((int)maxd1 + 1);
        OPTIMOTU_DEBUG(
          4,
          << "wfa_aligner min_k=" << min_k
          << " max_k=" << max_k
          << " max score=" << (int)maxd1 + 1
          << std::endl
        );
        double d = distance_wfa2(seq[s1], seq[s2], wfa_aligner);
        if (d < 1.0) ++my_aligned;
        OPTIMOTU_DEBUG(
          4,
          << "Thread " << begin
          << " distance=" << d
          << std::endl
        );
        if (d < threshold) clust_algo(*pg, d);
      }
      OPTIMOTU_DEBUG(
        4,
        << "Thread " << pg_index
        << ": finished " << j << " (j0=" << j0 << ")"
        << " and " << i << " (i0=" << i0 << ")"
        << "\n" << std::endl
      );
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

template class Wfa2SplitClusterWorker<0>;
template class Wfa2SplitClusterWorker<1>;
template class Wfa2SplitClusterWorker<2>;
template class Wfa2SplitClusterWorker<3>;
template class Wfa2SplitClusterWorker<4>;
template class Wfa2ConcurrentClusterWorker<0>;
template class Wfa2ConcurrentClusterWorker<1>;
template class Wfa2ConcurrentClusterWorker<2>;
template class Wfa2ConcurrentClusterWorker<3>;
template class Wfa2ConcurrentClusterWorker<4>;