// SPDX-FileCopyrightText: 2026 Brendan Furneaux
// SPDX-License-Identifier: MIT

#include "optimotu.h"
#include "Ksw2SearchWorker.h"
#include "Ksw2Aligner.h"
#include "pairwise_alignment.h"
#include "wfa_identity_bound.h"

#include <algorithm>

template<int verbose, bool do_cigar, enum AlignmentSpan span>
void Ksw2SearchWorkerImpl<verbose, do_cigar, span>::operator()(std::size_t begin, std::size_t end) {
  std::size_t begin_i = (begin * query.size()) / threads;
  std::size_t end_i = (end * query.size()) / threads;
  OPTIMOTU_DEBUG(
    2,
    << "Ksw2SearchWorker thread " << begin
    << " entered; sequences [" << begin_i
    << ", "<< end_i << ")" << std::endl
  );

  Ksw2Aligner aligner{match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2};
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
        4,
        << "thread" << begin
        << ": seqs " << j
        << " and " << i
        << " max relevant=" << max_dist
        << std::endl
      );

      bool is_query_longer = query[i].size() > ref[j].size();
      const std::string & s1 = is_query_longer ? ref[j] : query[i];
      const std::string & s2 = is_query_longer ? query[i] : ref[j];

      double l1 = s1.size(), l2 = s2.size();
      OPTIMOTU_DEBUG(
        4,
        << (is_query_longer ? "#### ref " : "#### query " )
        << (is_query_longer ? j : i)
        << " (l1=" << l1 << ") and "
        << (is_query_longer ? "query " : "ref " )
        << (is_query_longer ? i : j)
        << " (l2=" << l2 <<")####" << std::endl
      );

      double sim_threshold = 1.0 - max_dist;
      if constexpr (span == AlignmentSpan::GLOBAL) {
        if (l1/l2 < sim_threshold) continue;
      }
      ++my_prealigned;
      WfaIdentityBound bound = wfa_identity_bound(
          l1, l2, max_dist, match, mismatch,
          gap_open, gap_extend, gap_open2, gap_extend2,
          span == AlignmentSpan::EXTEND
      );
      aligner.setBandWidth(std::max(bound.max_k, -bound.min_k));
      double d;
      std::string cigar;
      // Longer as query, shorter as target
      if constexpr (do_cigar) {
        std::tie(d, cigar) = distance_and_cigar_ksw2<span>(s2, s1, aligner);
      } else {
        d = distance_ksw2<span>(s2, s1, aligner);
      }
      OPTIMOTU_DEBUG(
        4,
        << "Thread " << begin
        << ": distance=" << d
        << std::endl
      );
      if (d == 1.0) continue;
      ++my_aligned;
      if (d < hit.best_dist) {
        OPTIMOTU_DEBUG(
          4,
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
          4,
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
  OPTIMOTU_DEBUG(2, << "Exiting thread " << begin << std::endl);
}

template class Ksw2SearchWorkerImpl<0, true, AlignmentSpan::GLOBAL>;
template class Ksw2SearchWorkerImpl<1, true, AlignmentSpan::GLOBAL>;
template class Ksw2SearchWorkerImpl<2, true, AlignmentSpan::GLOBAL>;
template class Ksw2SearchWorkerImpl<3, true, AlignmentSpan::GLOBAL>;
template class Ksw2SearchWorkerImpl<4, true, AlignmentSpan::GLOBAL>;
template class Ksw2SearchWorkerImpl<0, false, AlignmentSpan::GLOBAL>;
template class Ksw2SearchWorkerImpl<1, false, AlignmentSpan::GLOBAL>;
template class Ksw2SearchWorkerImpl<2, false, AlignmentSpan::GLOBAL>;
template class Ksw2SearchWorkerImpl<3, false, AlignmentSpan::GLOBAL>;
template class Ksw2SearchWorkerImpl<4, false, AlignmentSpan::GLOBAL>;
template class Ksw2SearchWorkerImpl<0, true, AlignmentSpan::EXTEND>;
template class Ksw2SearchWorkerImpl<1, true, AlignmentSpan::EXTEND>;
template class Ksw2SearchWorkerImpl<2, true, AlignmentSpan::EXTEND>;
template class Ksw2SearchWorkerImpl<3, true, AlignmentSpan::EXTEND>;
template class Ksw2SearchWorkerImpl<4, true, AlignmentSpan::EXTEND>;
template class Ksw2SearchWorkerImpl<0, false, AlignmentSpan::EXTEND>;
template class Ksw2SearchWorkerImpl<1, false, AlignmentSpan::EXTEND>;
template class Ksw2SearchWorkerImpl<2, false, AlignmentSpan::EXTEND>;
template class Ksw2SearchWorkerImpl<3, false, AlignmentSpan::EXTEND>;
template class Ksw2SearchWorkerImpl<4, false, AlignmentSpan::EXTEND>;
