#include "Wfa2ClusterWorker.h"
#include <bindings/cpp/WFAligner.hpp>
#include "pairwise_alignment.h"

Wfa2ClusterWorker::Wfa2ClusterWorker(
  const std::vector<std::string> &seq,
  ClusterAlgorithm &clust_algo,
  DivisiblePairGenerator::Builder & pgb,
  const int match, const int mismatch,
  const int gap_open, const int gap_extend,
  const int gap_open2, const int gap_extend2
) : DistClusterWorker(seq, clust_algo, pgb),
match(match), mismatch(mismatch), gap_open(gap_open), gap_extend(gap_extend),
gap_open2(gap_open2), gap_extend2(gap_extend2) {};

template<int verbose>
void Wfa2SplitClusterWorker<verbose>::operator()(std::size_t begin, std::size_t end) {
  wfa::WFAlignerChoose wfa_aligner{match, mismatch, gap_open, gap_extend,
                                   gap_open2, gap_extend2, wfa::WFAligner::Alignment};
  for (size_t pg_index = begin; pg_index < end; pg_index++) {
    auto & pg = pair_generators[pg_index];
    ClusterAlgorithm * my_algo = clust_algo.make_child(pg.get());
    OPTIMOTU_DEBUG(1, << "Wfa2SplitClusterWorker thread " << pg_index << " entered" << std::endl);
    size_t my_prealigned = 0;
    size_t my_aligned = 0;
    while (*pg) {
      size_t i = pg->i();
      size_t j = pg->j();
      size_t i0 = pg->i0();
      size_t j0 = pg->j0();
      double threshold = my_algo->max_relevant(*pg);
      OPTIMOTU_DEBUG(
        2,
        << "thread" << pg_index
        << ": seqs " << j << " (j0=" << j0 << ")"
        << " and " << i << " (i0=" << i0 << ")"
        << " max relevant=" << threshold
        << std::endl
      );

      bool is_seqj_longer = seq[j0].size() > seq[i0].size();
      size_t s1 = is_seqj_longer ? i : j;
      size_t s2 = is_seqj_longer ? j : i;
      double l1 = seq[s1].size(), l2 = seq[s2].size();
      OPTIMOTU_DEBUG(
        2,
        << "#### seq " << i0 << " (l1=" << l1 << ") and "
        << j0 << " (l2=" << l2 <<")####" << std::endl
      );

      double sim_threshold = 1.0 - threshold; // compiler can probably do this?
      if (l1/l2 < sim_threshold) continue;
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
        2,
        << " distance=" << d
        << std::endl
      );
      if (d < threshold) (*my_algo)(*pg, d);
      RcppThread::checkUserInterrupt();
    }
    mutex.lock();
    OPTIMOTU_DEBUG(1, << "thread " << pg_index << " ready to merge" << std::endl);
    _aligned += my_aligned;
    _prealigned += my_prealigned;
    mutex.unlock();
    my_algo->merge_into_parent();
    OPTIMOTU_DEBUG(1, << "thread " << pg_index << " done" << std::endl);
  }
}

template <int verbose>
void Wfa2ConcurrentClusterWorker<verbose>::operator()(std::size_t begin, std::size_t end) {
  size_t my_prealigned = 0;
  size_t my_aligned = 0;

  wfa::WFAlignerChoose wfa_aligner{match, mismatch, gap_open, gap_extend,
                                   gap_open2, gap_extend2, wfa::WFAligner::Alignment};
  for (size_t pg_index = begin; pg_index < end; pg_index++) {
    auto & pg = pair_generators[pg_index];
    OPTIMOTU_DEBUG(
      1,
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
        2,
        << "thread" << pg_index
        << ": seqs " << j << " (j0=" << j0 << ")"
        << " and " << i << " (i0=" << i0 << ")"
        << " max relevant=" << threshold
        << std::endl
      );
      bool is_seqj_longer = seq[j0].size() > seq[i0].size();
      size_t s1 = is_seqj_longer ? i : j;
      size_t s2 = is_seqj_longer ? j : i;
      double l1 = seq[s1].size(), l2 = seq[s2].size();
      OPTIMOTU_DEBUG(
        2,
        << "#### seq " << s1 << " (l1=" << l1 << ") and "
        << s2 << " (l2=" << l2 <<")####" << std::endl
      );

      double sim_threshold = 1.0 - threshold; // compiler can probably do this?
      if (l1/l2 < sim_threshold) continue;
      ++my_prealigned;
      double sim_threshold_plus_1 = 2.0 - threshold;
      double maxd1 = threshold * (l1 + l2) / sim_threshold_plus_1;
      int max_k = (int)ceil((l2 - l1 * sim_threshold) / sim_threshold_plus_1);
      int min_k = -(int)ceil((l1 - l2 * sim_threshold) / sim_threshold_plus_1);
      wfa_aligner.setHeuristicBandedStatic(min_k, max_k);
      wfa_aligner.setMaxAlignmentSteps((int)maxd1 + 1);
      OPTIMOTU_DEBUG(
        2,
        << "wfa_aligner min_k=" << min_k
        << " max_k=" << max_k
        << " max score=" << (int)maxd1 + 1
        << std::endl
      );
      double d = distance_wfa2(seq[s1], seq[s2], wfa_aligner);
      if (d < 1.0) ++my_aligned;
      OPTIMOTU_DEBUG(
        2,
        << "Thread " << begin
        << " distance=" << d
        << std::endl
      );
      if (d < threshold) clust_algo(*pg, d);
      OPTIMOTU_DEBUG(
        2,
        << "Thread " << pg_index
        << ": finished " << j << " (j0=" << j0 << ")"
        << " and " << i << " (i0=" << i0 << ")"
        << "\n" << std::endl
      );
      RcppThread::checkUserInterrupt();
    }
    mutex.lock();
    _aligned += my_aligned;
    _prealigned += my_prealigned;
    OPTIMOTU_DEBUG(1, << "thread " << pg_index << " done" << std::endl);
    mutex.unlock();
  }
}

template class Wfa2SplitClusterWorker<0>;
template class Wfa2SplitClusterWorker<1>;
template class Wfa2SplitClusterWorker<2>;
template class Wfa2ConcurrentClusterWorker<0>;
template class Wfa2ConcurrentClusterWorker<1>;
template class Wfa2ConcurrentClusterWorker<2>;
