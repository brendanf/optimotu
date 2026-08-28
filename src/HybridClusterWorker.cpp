#include "HybridClusterWorker.h"
#include "ClusterAlgorithmFactory.h"
#include "config.h"

#include <bindings/cpp/WFAligner.hpp>
#include <edlib.h>
#include "pairwise_alignment.h"

typedef RcppParallel::RMatrix<int> matrix_t;

HybridClusterWorker::HybridClusterWorker(
  const std::vector<std::string> &seq,
  ClusterAlgorithm &clust_algo,
  DivisiblePairGenerator::Builder & pgb,
  const double breakpoint,
  int verbose,
  std::size_t worker_threads
) : DistClusterWorker(seq, clust_algo, pgb, verbose, worker_threads), breakpoint(breakpoint) {};

template<int verbose>
void HybridSplitClusterWorker<verbose>::operator()(std::size_t begin, std::size_t end) {

  EdlibAlignConfig ed_aligner = edlibNewAlignConfig(-1, EdlibAlignMode::EDLIB_MODE_NW, EdlibAlignTask::EDLIB_TASK_PATH, 0, 0);
  wfa::WFAlignerEdit wfa_aligner{wfa::WFAligner::Alignment};
  for (size_t pg_index = begin; pg_index < pair_generators.size(); pg_index += threads) {
    OPTIMOTU_DEBUG(
      2,
      << "HybridSplitClusterWorker thread " << pg_index
      << " entered" << std::endl
    );
    size_t my_prealigned = 0;
    size_t my_aligned = 0;
    auto & pg = pair_generators[pg_index];
    ClusterAlgorithm *my_algo = clust_algo.make_child(pg.get());
    if (!my_algo)
    {
      continue;
    }

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
        bool is_close = breakpoint >= 1 ? maxd1 < breakpoint : threshold < breakpoint;
        double d;
        if (is_close) {
          int max_k = (int)ceil((l2 - l1 * sim_threshold) / sim_threshold_plus_1);
          int min_k = -(int)ceil((l1 - l2 * sim_threshold) / sim_threshold_plus_1);
          wfa_aligner.setHeuristicBandedStatic(min_k, max_k);
          wfa_aligner.setMaxAlignmentSteps((int)maxd1 + 1);
          d = distance_wfa2<>(seq[s1], seq[s2], wfa_aligner);
        } else {
          ed_aligner.k = (int)maxd1 + 1;
          d = distance_edlib(seq[s1], seq[s2], ed_aligner);
        }
        if (d < 1.0) ++my_aligned;

        OPTIMOTU_DEBUG(
          2,
          << " distance=" << d
          << std::endl
        );
        if (d < threshold) (*my_algo)(*pg, d);
      }
      OPTIMOTU_DEBUG(
        2,
        << "thread" << pg_index
        << ": finished " << j << " (j0=" << j0 << ")"
        << " and " << i << " (i0=" << i0 << ")"
        << "\n" << std::endl
      );
      ++(*pg);
      RcppThread::checkUserInterrupt();
    }
    mutex.lock();
    OPTIMOTU_DEBUG(1, << "thread" << pg_index << " ready to merge" << std::endl);
    _aligned += my_aligned;
    _prealigned += my_prealigned;
    mutex.unlock();
    my_algo->merge_into_parent();
    clust_algo.release_child(my_algo);
    OPTIMOTU_DEBUG(1, << "thread" << pg_index << " done" << std::endl);
  }
}

template <int verbose>
void HybridConcurrentClusterWorker<verbose>::operator()(std::size_t begin, std::size_t end) {
  EdlibAlignConfig ed_aligner = edlibNewAlignConfig(
    -1,
    EdlibAlignMode::EDLIB_MODE_NW,
    EdlibAlignTask::EDLIB_TASK_PATH,
    0,
    0
  );
  wfa::WFAlignerEdit wfa_aligner{wfa::WFAligner::Alignment};
  for (size_t pg_index = begin; pg_index < pair_generators.size(); pg_index += threads) {
    auto & pg = pair_generators[pg_index];
    OPTIMOTU_DEBUG(
      2,
      << "HybridConcurrentClusterWorker thread " << pg_index
      << " entered" << std::endl
    );
    size_t my_prealigned = 0;
    size_t my_aligned = 0;
    while (*pg) {
      size_t i = pg->i();
      size_t j = pg->j();
      size_t i0 = pg->i0();
      size_t j0 = pg->j0();
      OPTIMOTU_DEBUG(
        4,
        << "Thread " << pg_index
        << ": seqs " << j << " (j0=" << j0 << ")"
        << " and " << i << " (i0=" << i0 << ")"
        << std::endl
      );
      double threshold = clust_algo.max_relevant(*pg);
      OPTIMOTU_DEBUG(
        4,
        << "thread" << pg_index
        << ": max relevant=" << threshold
        << std::endl
      );
      bool is_seqj_longer = seq[j0].size() > seq[i0].size();
      size_t s1 = is_seqj_longer ? i0 : j0;
      size_t s2 = is_seqj_longer ? j0 : i0;
      double l1 = seq[s1].size(), l2 = seq[s2].size();
      OPTIMOTU_DEBUG(
        4,
        << "#### seq " << s1
        << " (l1=" << l1
        << ") and "
        << s2
        << " (l2=" << l2
        << ")"
        << "####" << std::endl
      );

      double sim_threshold = 1.0 - threshold; // compiler can probably do this?
      if (l1/l2 >= sim_threshold) {
        ++my_prealigned;
        double sim_threshold_plus_1 = 2.0 - threshold;
        double maxd1 = threshold * (l1 + l2) / sim_threshold_plus_1;
        bool is_close = breakpoint >= 1 ? maxd1 < breakpoint : threshold < breakpoint;
        double d;
        if (is_close) {
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
          d = distance_wfa2<>(seq[s1], seq[s2], wfa_aligner);
        } else {
          ed_aligner.k = (int)maxd1 + 1;
          d = distance_edlib(seq[s1], seq[s2], ed_aligner);
        }
        if (d < 1.0) ++my_aligned;
        OPTIMOTU_DEBUG(
          4,
          << "thread" << pg_index
          << ": distance=" << d
          << std::endl
        );
        if (d < threshold) clust_algo(*pg, d);
      }
      OPTIMOTU_DEBUG(
        4,
        << "thread" << pg_index
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
    OPTIMOTU_DEBUG(2, << "thread" << pg_index << " done" << std::endl);
    mutex.unlock();
  }
}

template class HybridSplitClusterWorker<0>;
template class HybridSplitClusterWorker<1>;
template class HybridSplitClusterWorker<2>;
template class HybridSplitClusterWorker<3>;
template class HybridSplitClusterWorker<4>;
template class HybridConcurrentClusterWorker<0>;
template class HybridConcurrentClusterWorker<1>;
template class HybridConcurrentClusterWorker<2>;
template class HybridConcurrentClusterWorker<3>;
template class HybridConcurrentClusterWorker<4>;