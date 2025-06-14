#include "optimotu.h"
#include "Wfa2SearchWorker.h"
#include <bindings/cpp/WFAligner.hpp>
#include "pairwise_alignment.h"

template<int verbose, bool do_cigar, enum AlignmentSpan span>
void Wfa2SearchWorkerImpl<verbose, do_cigar, span>::operator()(std::size_t begin, std::size_t end) {
  std::size_t begin_i = (begin * query.size()) / threads;
  std::size_t end_i = (end * query.size()) / threads;
  OPTIMOTU_DEBUG(
    1,
    << "Wfa2SearchWorker thread " << begin
    << " entered; sequences [" << begin_i
    << ", "<< end_i << ")" << std::endl
  );

  wfa::WFAlignerChoose wfa_aligner{match, mismatch, gap_open, gap_extend,
                                   gap_open2, gap_extend2,
                                   wfa::WFAligner::Alignment};
  std::size_t my_prealigned = 0, my_aligned = 0;
  for (std::size_t i = begin_i; i < end_i; i++) {
    if (!hits[i]) {
      if constexpr (do_cigar) {
        hits[i] = std::make_unique<SearchCigarHit>();
      } else {
        hits[i] = std::make_unique<SearchHit>();
      }
    }
    SearchHit & hit = *hits[i];
    for (std::size_t j = 0; j < ref.size(); j++) {
      double max_dist = (hit.best_dist < threshold) ? hit.best_dist : threshold;
      OPTIMOTU_DEBUG(
        2,
        << "thread" << begin
        << ": seqs " << j
        << " and " << i
        << " max relevant=" << max_dist
        << std::endl
      );

      // s1 is the shorter sequence, s2 is the longer sequence
      bool is_query_longer = query[i].size() > ref[j].size();
      const std::string & s1 = is_query_longer ? ref[j] : query[i];
      const std::string & s2 = is_query_longer ? query[i] : ref[j];

      double l1 = s1.size(), l2 = s2.size();
      OPTIMOTU_DEBUG(
        2,
        << (is_query_longer ? "#### ref " : "#### query " )
        << (is_query_longer ? j : i)
        << " (l1=" << l1 << ") and "
        << (is_query_longer ? "query " : "ref " )
        << (is_query_longer ? i : j)
        << " (l2=" << l2 <<")####" << std::endl
      );

      double sim_threshold = 1.0 - max_dist; // compiler can probably do this?
      if constexpr (span == AlignmentSpan::GLOBAL) {
        if (l1/l2 < sim_threshold) continue;
      }
      ++my_prealigned;
      double sim_threshold_plus_1 = 2.0 - threshold;
      double maxd1 = threshold * (l1 + l2) / sim_threshold_plus_1;
      int max_k = (int)ceil((l2 - l1 * sim_threshold) / sim_threshold_plus_1);
      int min_k;
      if constexpr (span == AlignmentSpan::GLOBAL) {
        min_k = -(int)ceil((l1 - l2 * sim_threshold) / sim_threshold_plus_1);
      } else if constexpr (span == AlignmentSpan::EXTEND) {
        min_k = -(int)ceil(l1 * sim_threshold);
      } else {
        // static_assert(span != span, "Instatiation of non-implemented AlignmentSpan");
      }
      wfa_aligner.setHeuristicBandedStatic(min_k, max_k);
      wfa_aligner.setMaxAlignmentSteps((int)maxd1 + 1);
      double d;
      std::string cigar;
      if constexpr (do_cigar) {
        std::tie(d, cigar) = distance_and_cigar_wfa2<span>(s1, s2, wfa_aligner);
      } else {
        d = distance_wfa2<span>(s1, s2, wfa_aligner);
      }
      OPTIMOTU_DEBUG(
        2,
        << "Thread " << begin
        << ": distance=" << d
        << std::endl
      );
      if (d == 1.0) continue;
      ++my_aligned;
      if (d < hit.best_dist) {
        OPTIMOTU_DEBUG(
          2,
          << "Thread " << begin
          << ": new best distance=" << d
          << std::endl
        );
        hit.best_dist = d;
        hit.best_ref.clear();
        hit.best_ref.push_back(j);
        if constexpr (do_cigar) {
          hit.best_cigar().clear();
          hit.best_cigar().push_back(cigar);
        }
      } else if (d == hit.best_dist) {
        OPTIMOTU_DEBUG(
          2,
          << "Thread " << begin
          << ": tied for best distance=" << d
          << std::endl
        );
        hit.best_ref.push_back(j);
        if constexpr (do_cigar) {
          hit.best_cigar().push_back(cigar);
        }
      }
      RcppThread::checkUserInterrupt();
    }
  }
  {
    std::lock_guard<std::mutex> lock(mutex);
    _prealigned += my_prealigned;
    _aligned += my_aligned;
  }
  OPTIMOTU_DEBUG(1, << "Exiting thread " << begin << std::endl);
}

template class Wfa2SearchWorkerImpl<0, true, AlignmentSpan::GLOBAL>;
template class Wfa2SearchWorkerImpl<1, true, AlignmentSpan::GLOBAL>;
template class Wfa2SearchWorkerImpl<2, true, AlignmentSpan::GLOBAL>;
template class Wfa2SearchWorkerImpl<0, false, AlignmentSpan::GLOBAL>;
template class Wfa2SearchWorkerImpl<1, false, AlignmentSpan::GLOBAL>;
template class Wfa2SearchWorkerImpl<2, false, AlignmentSpan::GLOBAL>;
template class Wfa2SearchWorkerImpl<0, true, AlignmentSpan::EXTEND>;
template class Wfa2SearchWorkerImpl<1, true, AlignmentSpan::EXTEND>;
template class Wfa2SearchWorkerImpl<2, true, AlignmentSpan::EXTEND>;
template class Wfa2SearchWorkerImpl<0, false, AlignmentSpan::EXTEND>;
template class Wfa2SearchWorkerImpl<1, false, AlignmentSpan::EXTEND>;
template class Wfa2SearchWorkerImpl<2, false, AlignmentSpan::EXTEND>;
