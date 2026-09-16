#ifdef OPTIMOTU_R
#include "EdlibDistWorker.h"
#include "tile_pairs.h"
#include <edlib.h>
#include <algorithm>
#include "pairwise_alignment.h"
#include "optimotu.h"

template<int verbose, bool is_constrained, enum AlignmentSpan span, typename SparseDistanceMatrixType>
void EdlibDistWorkerImpl<verbose, is_constrained, span, SparseDistanceMatrixType>::operator()(std::size_t begin, std::size_t end) {
  (void)end;
  size_t my_prealigned = 0;
  size_t my_aligned = 0;
  std::vector<size_t> my_seq1;
  my_seq1.reserve(1000);
  std::vector<size_t> my_seq2;
  my_seq2.reserve(1000);
  std::vector<int> my_score1;
  my_score1.reserve(1000);
  std::vector<int> my_score2;
  my_score2.reserve(1000);
  std::vector<double> my_dist1;
  my_dist1.reserve(1000);
  std::vector<double> my_dist2;
  my_dist2.reserve(1000);
  std::vector<std::string> my_cigar;
  if constexpr (std::is_same_v<SparseDistanceMatrixType, SparseDistanceMatrixCigar>) {
    my_cigar.reserve(1000);
  }
  std::vector<int> my_align_length;
  std::vector<int> my_n_insert;
  std::vector<int> my_n_delete;
  std::vector<int> my_max_insert;
  std::vector<int> my_max_delete;

  if constexpr (std::is_same_v<SparseDistanceMatrixType, SparseDistanceMatrixGapstats>) {
    my_align_length.reserve(1000);
    my_n_insert.reserve(1000);
    my_n_delete.reserve(1000);
    my_max_insert.reserve(1000);
    my_max_delete.reserve(1000);
  }

  EdlibAlignMode mode;
  if constexpr (span == AlignmentSpan::GLOBAL) {
    mode = EdlibAlignMode::EDLIB_MODE_NW;
  } else if constexpr (span == AlignmentSpan::EXTEND) {
    mode = EdlibAlignMode::EDLIB_MODE_SHW;
  } else {
    static_assert(always_false_v<span>, "Invalid alignment span");
  }

  EdlibAlignConfig aligner =
    edlibNewAlignConfig(-1, mode, EdlibAlignTask::EDLIB_TASK_PATH, 0, 0);

  if constexpr (verbose >= 2) {
    RcppThread::Rcerr << "EdlibDistWorker thread " << begin
                      << " entered" << std::endl;
  }

  auto flush_if_full = [&]() {
    if (my_seq1.size() != 1000) {
      return;
    }
    if constexpr (std::is_same_v<SparseDistanceMatrixType, SparseDistanceMatrixCigar>) {
      sdm.append(
        my_seq1, my_seq2, my_score1, my_score2, my_dist1, my_dist2, my_cigar
      );
    } else if constexpr (std::is_same_v<SparseDistanceMatrixType, SparseDistanceMatrixGapstats>) {
      sdm.append(
        my_seq1, my_seq2, my_score1, my_score2, my_dist1, my_dist2,
        my_align_length, my_n_insert, my_n_delete, my_max_insert, my_max_delete
      );
    } else {
      sdm.append(my_seq1, my_seq2, my_score1, my_score2, my_dist1, my_dist2);
    }
  };

  for_each_assigned_tile(
    pair_generators, begin, threads,
    [&](std::size_t /*pg_index*/, PairGenerator &pg) {
      for_each_tile_pair(pg, [&](
        std::size_t i0, std::size_t j0, std::size_t, std::size_t
      ) {
        const size_t i = i0;
        const size_t j = j0;

        bool is_seqj_longer = seq[j].size() > seq[i].size();
        size_t s1 = is_seqj_longer ? i : j;
        size_t s2 = is_seqj_longer ? j : i;
        double l1 = seq[s1].size(), l2 = seq[s2].size();
        if constexpr (verbose >= 3) {
          RcppThread::Rcerr << "#### seq " << i << " (l1=" << l1 << ") and "
                            << j << " (l2=" << l2 <<") ####" << std::endl;
        }

        if (l1/l2 < sim_threshold) return;

        if constexpr (is_constrained) {
          if constexpr (span == AlignmentSpan::GLOBAL) {
            double maxd1 = dist_threshold * (l1 + l2) / sim_threshold_plus_1;
            aligner.k = (int)maxd1 + 1;
          } else {
            OPTIMOTU_STOP(
              "Alignment span not implemented for constrained banding"
            );
          }
        }
        auto alignResult = edlibAlign(
          seq[s1].c_str(), seq[s1].size(),
          seq[s2].c_str(), seq[s2].size(),
          aligner
        );
        my_prealigned++;
        if (alignResult.status == EDLIB_STATUS_ERROR) {
          OPTIMOTU_DEBUG(4, << "edlibAlign error " << j << std::endl);
          edlibFreeAlignResult(alignResult);
          return;
        }
        if (alignResult.editDistance == -1) {
          OPTIMOTU_DEBUG(
            4,
            << "edlibAlign no alignment within band " << aligner.k
            << std::endl
          );
          edlibFreeAlignResult(alignResult);
          return;
        }
        my_aligned++;
        double d2 = (double)alignResult.editDistance /
          (double) alignResult.alignmentLength;
        if (d2 <= dist_threshold) {
          my_seq1.push_back(j);
          my_seq2.push_back(i);
          my_score1.push_back(0);
          my_score2.push_back(alignResult.editDistance);
          my_dist1.push_back(0);
          my_dist2.push_back(d2);
          if constexpr (std::is_same_v<SparseDistanceMatrixType, SparseDistanceMatrixCigar>) {
            std::string cigar;
            cigar.reserve(alignResult.alignmentLength);
            std::vector<char> key = {'M', 'I', 'D', 'X'};
            for (int k = 0; k < alignResult.alignmentLength; k++) {
              cigar.push_back(key[alignResult.alignment[k]]);
            }
            my_cigar.push_back(cigar);
          }
          if constexpr (std::is_same_v<SparseDistanceMatrixType, SparseDistanceMatrixGapstats>) {
            my_align_length.push_back(alignResult.alignmentLength);
            int ni = 0, nd = 0, mi = 0, md = 0, ii = 0, d = 0;
            for (int k = 0; k < alignResult.alignmentLength; k++) {
              if (alignResult.alignment[k] == 'I') {
                ni++;
                ii++;
                if (ii > mi) mi = ii;
                d = 0;
              } else if (alignResult.alignment[k] == 'D') {
                nd++;
                d++;
                if (d > md) md = d;
                ii = 0;
              } else {
                ii = 0;
                d = 0;
              }
            }
            my_n_insert.push_back(ni);
            my_n_delete.push_back(nd);
            my_max_insert.push_back(mi);
            my_max_delete.push_back(md);
          }
          flush_if_full();
        } else {
          OPTIMOTU_DEBUG(
            4,
            << "edlibAlign distance " << d2 << " above threshold "
            << dist_threshold << std::endl
          );
        }
        edlibFreeAlignResult(alignResult);
        RcppThread::checkUserInterrupt();
      });
    }
  );
  if (my_seq1.size() > 0) {
    if constexpr (std::is_same_v<SparseDistanceMatrixType, SparseDistanceMatrixCigar>) {
      sdm.append(
        my_seq1, my_seq2, my_score1, my_score2, my_dist1, my_dist2, my_cigar
      );
    } else if constexpr (std::is_same_v<SparseDistanceMatrixType, SparseDistanceMatrixGapstats>) {
      sdm.append(
        my_seq1, my_seq2, my_score1, my_score2, my_dist1, my_dist2,
        my_align_length, my_n_insert, my_n_delete, my_max_insert, my_max_delete
      );
    } else {
      sdm.append(my_seq1, my_seq2, my_score1, my_score2, my_dist1, my_dist2);
    }
  }
  sdm.mutex.lock();
  _aligned += my_aligned;
  _prealigned += my_prealigned;
  sdm.mutex.unlock();
}

template<int verbose, bool is_constrained, enum AlignmentSpan span>
std::unique_ptr<EdlibDistWorker> create_edlib_dist_worker(
  const std::vector<std::string> &seq,
  const double dist_threshold,
  DivisiblePairGenerator::Builder & pgb,
  SparseDistanceMatrix &sdm,
  std::size_t worker_threads
) {
  if (dynamic_cast<SparseDistanceMatrixCigar*>(&sdm)) {
    return std::make_unique<EdlibDistWorkerImpl<verbose, is_constrained, span, SparseDistanceMatrixCigar>>(
      seq, dist_threshold, pgb, sdm, worker_threads
    );
  } else if (dynamic_cast<SparseDistanceMatrixGapstats*>(&sdm)) {
    return std::make_unique<EdlibDistWorkerImpl<verbose, is_constrained, span, SparseDistanceMatrixGapstats>>(
      seq, dist_threshold, pgb, sdm, worker_threads
    );
  } else {
    return std::make_unique<EdlibDistWorkerImpl<verbose, is_constrained, span, SparseDistanceMatrix>>(
      seq, dist_threshold, pgb, sdm, worker_threads
    );
  }
}

template<int verbose, bool is_constrained>
std::unique_ptr<EdlibDistWorker> create_edlib_dist_worker(
  const std::vector<std::string> &seq,
  const double dist_threshold,
  DivisiblePairGenerator::Builder & pgb,
  SparseDistanceMatrix &sdm,
  AlignmentSpan span,
  std::size_t worker_threads
) {
  switch (span) {
    case AlignmentSpan::GLOBAL:
      return create_edlib_dist_worker<verbose, is_constrained, AlignmentSpan::GLOBAL>(
        seq, dist_threshold, pgb, sdm, worker_threads
      );
    case AlignmentSpan::EXTEND:
      return create_edlib_dist_worker<verbose, is_constrained, AlignmentSpan::EXTEND>(
        seq, dist_threshold, pgb, sdm, worker_threads
      );
    default:
      OPTIMOTU_STOP("span must be 0 or 1");
  }
}

template<int verbose>
std::unique_ptr<EdlibDistWorker> create_edlib_dist_worker(
  const std::vector<std::string> &seq,
  const double dist_threshold,
  DivisiblePairGenerator::Builder & pgb,
  SparseDistanceMatrix &sdm,
  AlignmentSpan span,
  bool constrain,
  std::size_t worker_threads
) {
  if (constrain) {
    return create_edlib_dist_worker<verbose, true>(
      seq, dist_threshold, pgb, sdm, span, worker_threads
    );
  } else {
    return create_edlib_dist_worker<verbose, false>(
      seq, dist_threshold, pgb, sdm, span, worker_threads
    );
  }
}

std::unique_ptr<EdlibDistWorker> create_edlib_dist_worker(
  const std::vector<std::string> &seq,
  const double dist_threshold,
  DivisiblePairGenerator::Builder & pgb,
  SparseDistanceMatrix &sdm,
  int verbose,
  AlignmentSpan span,
  bool constrain,
  std::size_t worker_threads
) {
  int v = (verbose > 4) ? 4 : verbose;
  switch (v) {
    case 0:
      return create_edlib_dist_worker<0>(
        seq, dist_threshold, pgb, sdm, span, constrain, worker_threads
      );
    case 1:
      return create_edlib_dist_worker<1>(
        seq, dist_threshold, pgb, sdm, span, constrain, worker_threads
      );
    case 2:
      return create_edlib_dist_worker<2>(
        seq, dist_threshold, pgb, sdm, span, constrain, worker_threads
      );
    case 3:
      return create_edlib_dist_worker<3>(
        seq, dist_threshold, pgb, sdm, span, constrain, worker_threads
      );
    case 4:
    default:
      return create_edlib_dist_worker<4>(
        seq, dist_threshold, pgb, sdm, span, constrain, worker_threads
      );
  }
}



template class EdlibDistWorkerImpl<0, false, AlignmentSpan::GLOBAL>;
template class EdlibDistWorkerImpl<1, false, AlignmentSpan::GLOBAL>;
template class EdlibDistWorkerImpl<2, false, AlignmentSpan::GLOBAL>;
template class EdlibDistWorkerImpl<3, false, AlignmentSpan::GLOBAL>;
template class EdlibDistWorkerImpl<4, false, AlignmentSpan::GLOBAL>;
template class EdlibDistWorkerImpl<0, true, AlignmentSpan::GLOBAL>;
template class EdlibDistWorkerImpl<1, true, AlignmentSpan::GLOBAL>;
template class EdlibDistWorkerImpl<2, true, AlignmentSpan::GLOBAL>;
template class EdlibDistWorkerImpl<3, true, AlignmentSpan::GLOBAL>;
template class EdlibDistWorkerImpl<4, true, AlignmentSpan::GLOBAL>;
template class EdlibDistWorkerImpl<0, false, AlignmentSpan::EXTEND>;
template class EdlibDistWorkerImpl<1, false, AlignmentSpan::EXTEND>;
template class EdlibDistWorkerImpl<2, false, AlignmentSpan::EXTEND>;
template class EdlibDistWorkerImpl<3, false, AlignmentSpan::EXTEND>;
template class EdlibDistWorkerImpl<4, false, AlignmentSpan::EXTEND>;
template class EdlibDistWorkerImpl<0, true, AlignmentSpan::EXTEND>;
template class EdlibDistWorkerImpl<1, true, AlignmentSpan::EXTEND>;
template class EdlibDistWorkerImpl<2, true, AlignmentSpan::EXTEND>;
template class EdlibDistWorkerImpl<3, true, AlignmentSpan::EXTEND>;
template class EdlibDistWorkerImpl<4, true, AlignmentSpan::EXTEND>;

#endif
