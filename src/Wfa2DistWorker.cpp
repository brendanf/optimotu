#ifdef OPTIMOTU_R
#include "Wfa2DistWorker.h"
#include <bindings/cpp/WFAligner.hpp>
#include "pairwise_alignment.h"
#include "optimotu.h"

Wfa2DistWorker::Wfa2DistWorker(
  const std::vector<std::string> &seq,
  const double dist_threshold,
  const std::uint8_t threads,
  SparseDistanceMatrix &sdm,
  int match, int mismatch,
  int gap_open, int gap_extend,
  int gap_open2, int gap_extend2
) : DistWorker(seq, dist_threshold, threads, sdm),
match(match), mismatch(mismatch),
gap_open(gap_open), gap_extend(gap_extend),
gap_open2(gap_open2), gap_extend2(gap_extend2) {}

template<int verbose, bool is_constrained, enum AlignmentSpan span, typename SparseDistanceMatrixType>
void Wfa2DistWorkerImpl<verbose, is_constrained, span, SparseDistanceMatrixType>::operator()(std::size_t begin, std::size_t end) {
  double n = seq.size();
  double m = (n*n - 3.0*n + 2.0)/2.0;
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

  size_t begin_i;
  if (begin == 0) {
    begin_i = 1;
  } else {
    begin_i = round(1.5 + 0.5*sqrt(9.0 + 8.0*((m*begin)/threads - 1.0)));
  }
  size_t end_i = round(1.5 + 0.5*sqrt(9.0 + 8.0*((m*end)/threads - 1.0)));

  OPTIMOTU_DEBUG(
    1,
    << "Wfa2DistWorker thread " << begin
    << " entered; sequences [" << begin_i
    << ", "<< end_i << ")" << std::endl
  );

  OPTIMOTU_DEBUG(2,
    << "verbose = " << verbose
    << "\nis_constrained = " << is_constrained
    << "\nspan = " << (span == AlignmentSpan::GLOBAL ? "GLOBAL" : "EXTEND")
    << "\nSparseDistanceMatrixType = " << typeid(SparseDistanceMatrixType).name() 
    << std::endl
  );

  for (size_t i = begin_i; i < end_i; i++) {
    for (size_t j = 0; j < i; j++) {

      bool is_seqj_longer = seq[j].size() > seq[i].size();
      size_t s1 = is_seqj_longer ? i : j;
      size_t s2 = is_seqj_longer ? j : i;
      double l1 = seq[s1].size(), l2 = seq[s2].size();

      if constexpr (span == AlignmentSpan::GLOBAL) {
        if (l1/l2 < sim_threshold) continue;
      }
      ++my_prealigned;
      if constexpr (is_constrained) {
        double max_d1 = dist_threshold * (l1 + l2) / sim_threshold_plus_1;
        int max_k = (int)ceil((l2 - l1 * sim_threshold) / sim_threshold_plus_1);
        int min_k;
        if constexpr (span == AlignmentSpan::GLOBAL) {
          min_k = -(int)ceil((l1 - l2 * sim_threshold) / sim_threshold_plus_1);
        } else if constexpr (span == AlignmentSpan::EXTEND) {
          min_k = -(int)ceil(l1 * sim_threshold);
        } else {
          static_assert(span != span, "Instatiation of non-implemented AlignmentSpan");
        }
        wfa_aligner.setHeuristicBandedStatic(min_k, max_k);
        wfa_aligner.setMaxAlignmentSteps((int)max_d1 + 1);
      }

      std::string cigar = cigar_wfa2<span>(seq[s1], seq[s2], wfa_aligner);
      if (cigar == "") continue;
      ++my_aligned;
      int score = wfa_aligner.getAlignmentScore();
      double d = distance_from_cigar(cigar);
      if (d > dist_threshold) continue;

      my_seq1.push_back(j);
      my_seq2.push_back(i);
      my_score1.push_back(0);
      my_score2.push_back(score);
      my_dist1.push_back(0);
      my_dist2.push_back(d);

      if constexpr (std::is_same_v<SparseDistanceMatrixType, SparseDistanceMatrixCigar>) {
        my_cigar.push_back(cigar);
        if (my_seq1.size() == 1000) {
          sdm.append(my_seq1, my_seq2, my_score1, my_score2, my_dist1, my_dist2, my_cigar);
        }
      } else if constexpr (std::is_same_v<SparseDistanceMatrixType, SparseDistanceMatrixGapstats>) {
        auto [d, aln_len, n_ins, n_del, max_ins, max_del] = dist_gapstats_from_cigar(cigar);
        my_align_length.push_back(aln_len);
        my_n_insert.push_back(n_ins);
        my_n_delete.push_back(n_del);
        my_max_insert.push_back(max_ins);
        my_max_delete.push_back(max_del);
        if (my_seq1.size() == 1000) {
          sdm.append(my_seq1, my_seq2, my_score1, my_score2, my_dist1, my_dist2, my_align_length, my_n_insert, my_n_delete, my_max_insert, my_max_delete);
        }
      } else {
        if (my_seq1.size() == 1000) {
          sdm.append(my_seq1, my_seq2, my_score1, my_score2, my_dist1, my_dist2);
        }
      }
      RcppThread::checkUserInterrupt();
    }
  }
  if (my_seq1.size() > 0) {
    if constexpr (std::is_same_v<SparseDistanceMatrixType, SparseDistanceMatrixCigar>) {
      sdm.append(my_seq1, my_seq2, my_score1, my_score2, my_dist1, my_dist2, my_cigar);
    } else if constexpr (std::is_same_v<SparseDistanceMatrixType, SparseDistanceMatrixGapstats>) {
      sdm.append(my_seq1, my_seq2, my_score1, my_score2, my_dist1, my_dist2, my_align_length, my_n_insert, my_n_delete, my_max_insert, my_max_delete);
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
  const std::uint8_t threads,
  SparseDistanceMatrix &sdm,
  int match, int mismatch,
  int gap_open, int gap_extend,
  int gap_open2, int gap_extend2
) {
  if (dynamic_cast<SparseDistanceMatrixCigar*>(&sdm)) {
    return std::make_unique<Wfa2DistWorkerImpl<verbose, is_constrained, span, SparseDistanceMatrixCigar>>(seq, dist_threshold, threads, sdm, match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2);
  } else if (dynamic_cast<SparseDistanceMatrixGapstats*>(&sdm)) {
    return std::make_unique<Wfa2DistWorkerImpl<verbose, is_constrained, span, SparseDistanceMatrixGapstats>>(seq, dist_threshold, threads, sdm, match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2);
  } else {
    return std::make_unique<Wfa2DistWorkerImpl<verbose, is_constrained, span, SparseDistanceMatrix>>(seq, dist_threshold, threads, sdm, match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2);
  }
}

template<int verbose, bool is_constrained>
std::unique_ptr<Wfa2DistWorker> create_wfa2_dist_worker(
  const std::vector<std::string> &seq,
  const double dist_threshold,
  const std::uint8_t threads,
  SparseDistanceMatrix &sdm,
  int match, int mismatch,
  int gap_open, int gap_extend,
  int gap_open2, int gap_extend2,
  AlignmentSpan span
) {
  switch (span) {
    case AlignmentSpan::GLOBAL:
      return create_wfa2_dist_worker<verbose, is_constrained, AlignmentSpan::GLOBAL>(seq, dist_threshold, threads, sdm, match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2);
    case AlignmentSpan::EXTEND:
      return create_wfa2_dist_worker<verbose, is_constrained, AlignmentSpan::EXTEND>(seq, dist_threshold, threads, sdm, match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2);
    default:
      OPTIMOTU_STOP("span must be 0 or 1");
  }
}

template<int verbose>
std::unique_ptr<Wfa2DistWorker> create_wfa2_dist_worker(
  const std::vector<std::string> &seq,
  const double dist_threshold,
  const std::uint8_t threads,
  SparseDistanceMatrix &sdm,
  int match, int mismatch,
  int gap_open, int gap_extend,
  int gap_open2, int gap_extend2,
  AlignmentSpan span,
  bool is_constrained
) {
  if (is_constrained) {
    return create_wfa2_dist_worker<verbose, true>(seq, dist_threshold, threads, sdm, match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2, span);
  } else {
    return create_wfa2_dist_worker<verbose, false>(seq, dist_threshold, threads, sdm, match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2, span);
  }
}

std::unique_ptr<Wfa2DistWorker> create_wfa2_dist_worker(
  const std::vector<std::string> &seq,
  const double dist_threshold,
  const std::uint8_t threads,
  SparseDistanceMatrix &sdm,
  int match, int mismatch,
  int gap_open, int gap_extend,
  int gap_open2, int gap_extend2,
  int verbose,
  AlignmentSpan span,
  bool constrain
) {
  switch (verbose) {
    case 0:
      return create_wfa2_dist_worker<0>(seq, dist_threshold, threads, sdm, match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2, span, constrain);
    case 1:
      return create_wfa2_dist_worker<1>(seq, dist_threshold, threads, sdm, match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2, span, constrain);
    default:
      return create_wfa2_dist_worker<2>(seq, dist_threshold, threads, sdm, match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2, span, constrain);
  }
}



#endif //OPTIMOTU_R
