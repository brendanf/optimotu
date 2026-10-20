#ifdef OPTIMOTU_R
#include "Wfa2DistWorker.h"
#include "tile_pairs.h"
#include <bindings/cpp/WFAligner.hpp>
#include <algorithm>
#include "pairwise_alignment.h"
#include "optimotu.h"
#include "wfa_identity_bound.h"

Wfa2DistWorker::Wfa2DistWorker(
  const std::vector<std::string> &seq,
  const double dist_threshold,
  DivisiblePairGenerator::Builder & pgb,
  SparseDistanceMatrix &sdm,
  int match, int mismatch,
  int gap_open, int gap_extend,
  int gap_open2, int gap_extend2,
  int verbose,
  std::size_t worker_threads
) : DistWorker(seq, dist_threshold, pgb, sdm, verbose, worker_threads),
match(match), mismatch(mismatch),
gap_open(gap_open), gap_extend(gap_extend),
gap_open2(gap_open2), gap_extend2(gap_extend2) {}

template<int verbose, bool is_constrained, enum AlignmentSpan span, typename SparseDistanceMatrixType>
void Wfa2DistWorkerImpl<verbose, is_constrained, span, SparseDistanceMatrixType>::operator()(std::size_t begin, std::size_t end) {
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

  wfa::WFAlignerChoose wfa_aligner{match, mismatch, gap_open, gap_extend,
                                   gap_open2, gap_extend2,
                                   wfa::WFAligner::Alignment};

  OPTIMOTU_DEBUG(
    2,
    << "Wfa2DistWorker thread " << begin
    << " entered" << std::endl
  );

  OPTIMOTU_DEBUG(2,
    << "verbose = " << verbose
    << "\nis_constrained = " << is_constrained
    << "\nspan = " << (span == AlignmentSpan::GLOBAL ? "GLOBAL" : "EXTEND")
    << "\nSparseDistanceMatrixType = " << typeid(SparseDistanceMatrixType).name()
    << std::endl
  );

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
    pair_generators, begin, end,
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

        if constexpr (span == AlignmentSpan::GLOBAL) {
          if (l1/l2 < sim_threshold) return;
        }
        ++my_prealigned;
        if constexpr (is_constrained) {
          WfaIdentityBound bound = wfa_identity_bound(
              l1, l2, dist_threshold, match, mismatch,
              gap_open, gap_extend, gap_open2, gap_extend2,
              span == AlignmentSpan::EXTEND
          );
          wfa_aligner.setHeuristicBandedStatic(bound.min_k, bound.max_k);
          wfa_aligner.setMaxAlignmentSteps(bound.max_alignment_steps);
        }

        std::string cigar = cigar_wfa2<span>(seq[s1], seq[s2], wfa_aligner);
        if (cigar == "") return;
        ++my_aligned;
        int score = wfa_aligner.getAlignmentScore();
        double d = distance_from_cigar(cigar);
        if (d > dist_threshold) return;

        my_seq1.push_back(j);
        my_seq2.push_back(i);
        my_score1.push_back(0);
        my_score2.push_back(score);
        my_dist1.push_back(0);
        my_dist2.push_back(d);

        if constexpr (std::is_same_v<SparseDistanceMatrixType, SparseDistanceMatrixCigar>) {
          my_cigar.push_back(cigar);
        } else if constexpr (std::is_same_v<SparseDistanceMatrixType, SparseDistanceMatrixGapstats>) {
          auto [d2, aln_len, n_ins, n_del, max_ins, max_del] =
            dist_gapstats_from_cigar(cigar);
          (void)d2;
          my_align_length.push_back(aln_len);
          my_n_insert.push_back(n_ins);
          my_n_delete.push_back(n_del);
          my_max_insert.push_back(max_ins);
          my_max_delete.push_back(max_del);
        }
        flush_if_full();
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
std::unique_ptr<Wfa2DistWorker> create_wfa2_dist_worker(
  const std::vector<std::string> &seq,
  const double dist_threshold,
  DivisiblePairGenerator::Builder & pgb,
  SparseDistanceMatrix &sdm,
  int match, int mismatch,
  int gap_open, int gap_extend,
  int gap_open2, int gap_extend2,
  std::size_t worker_threads
) {
  if (dynamic_cast<SparseDistanceMatrixCigar*>(&sdm)) {
    return std::make_unique<Wfa2DistWorkerImpl<verbose, is_constrained, span, SparseDistanceMatrixCigar>>(
      seq, dist_threshold, pgb, sdm, match, mismatch, gap_open, gap_extend,
      gap_open2, gap_extend2, worker_threads
    );
  } else if (dynamic_cast<SparseDistanceMatrixGapstats*>(&sdm)) {
    return std::make_unique<Wfa2DistWorkerImpl<verbose, is_constrained, span, SparseDistanceMatrixGapstats>>(
      seq, dist_threshold, pgb, sdm, match, mismatch, gap_open, gap_extend,
      gap_open2, gap_extend2, worker_threads
    );
  } else {
    return std::make_unique<Wfa2DistWorkerImpl<verbose, is_constrained, span, SparseDistanceMatrix>>(
      seq, dist_threshold, pgb, sdm, match, mismatch, gap_open, gap_extend,
      gap_open2, gap_extend2, worker_threads
    );
  }
}

template<int verbose, bool is_constrained>
std::unique_ptr<Wfa2DistWorker> create_wfa2_dist_worker(
  const std::vector<std::string> &seq,
  const double dist_threshold,
  DivisiblePairGenerator::Builder & pgb,
  SparseDistanceMatrix &sdm,
  int match, int mismatch,
  int gap_open, int gap_extend,
  int gap_open2, int gap_extend2,
  AlignmentSpan span,
  std::size_t worker_threads
) {
  switch (span) {
    case AlignmentSpan::GLOBAL:
      return create_wfa2_dist_worker<verbose, is_constrained, AlignmentSpan::GLOBAL>(
        seq, dist_threshold, pgb, sdm, match, mismatch, gap_open, gap_extend,
        gap_open2, gap_extend2, worker_threads
      );
    case AlignmentSpan::EXTEND:
      return create_wfa2_dist_worker<verbose, is_constrained, AlignmentSpan::EXTEND>(
        seq, dist_threshold, pgb, sdm, match, mismatch, gap_open, gap_extend,
        gap_open2, gap_extend2, worker_threads
      );
    default:
      OPTIMOTU_STOP("span must be 0 or 1");
  }
}

template<int verbose>
std::unique_ptr<Wfa2DistWorker> create_wfa2_dist_worker(
  const std::vector<std::string> &seq,
  const double dist_threshold,
  DivisiblePairGenerator::Builder & pgb,
  SparseDistanceMatrix &sdm,
  int match, int mismatch,
  int gap_open, int gap_extend,
  int gap_open2, int gap_extend2,
  AlignmentSpan span,
  bool is_constrained,
  std::size_t worker_threads
) {
  if (is_constrained) {
    return create_wfa2_dist_worker<verbose, true>(
      seq, dist_threshold, pgb, sdm, match, mismatch, gap_open, gap_extend,
      gap_open2, gap_extend2, span, worker_threads
    );
  } else {
    return create_wfa2_dist_worker<verbose, false>(
      seq, dist_threshold, pgb, sdm, match, mismatch, gap_open, gap_extend,
      gap_open2, gap_extend2, span, worker_threads
    );
  }
}

std::unique_ptr<Wfa2DistWorker> create_wfa2_dist_worker(
  const std::vector<std::string> &seq,
  const double dist_threshold,
  DivisiblePairGenerator::Builder & pgb,
  SparseDistanceMatrix &sdm,
  int match, int mismatch,
  int gap_open, int gap_extend,
  int gap_open2, int gap_extend2,
  int verbose,
  AlignmentSpan span,
  bool constrain,
  std::size_t worker_threads
) {
  int v = (verbose > 4) ? 4 : verbose;
  switch (v) {
    case 0:
      return create_wfa2_dist_worker<0>(
        seq, dist_threshold, pgb, sdm, match, mismatch, gap_open, gap_extend,
        gap_open2, gap_extend2, span, constrain, worker_threads
      );
    case 1:
      return create_wfa2_dist_worker<1>(
        seq, dist_threshold, pgb, sdm, match, mismatch, gap_open, gap_extend,
        gap_open2, gap_extend2, span, constrain, worker_threads
      );
    case 2:
      return create_wfa2_dist_worker<2>(
        seq, dist_threshold, pgb, sdm, match, mismatch, gap_open, gap_extend,
        gap_open2, gap_extend2, span, constrain, worker_threads
      );
    case 3:
      return create_wfa2_dist_worker<3>(
        seq, dist_threshold, pgb, sdm, match, mismatch, gap_open, gap_extend,
        gap_open2, gap_extend2, span, constrain, worker_threads
      );
    case 4:
    default:
      return create_wfa2_dist_worker<4>(
        seq, dist_threshold, pgb, sdm, match, mismatch, gap_open, gap_extend,
        gap_open2, gap_extend2, span, constrain, worker_threads
      );
  }
}



#endif //OPTIMOTU_R
