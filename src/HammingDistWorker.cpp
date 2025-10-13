#include "HammingDistWorker.h"

#ifdef OPTIMOTU_R

HammingDistWorker::HammingDistWorker(
  const std::vector<std::string> &seq,
  const double dist_threshold,
  const std::uint8_t threads,
  SparseDistanceMatrix &sdm,
  const int min_overlap,
  const bool ignore_gap
) : DistWorker(seq, dist_threshold, threads, sdm),
pss(seq),
min_overlap(min_overlap),
ignore_gap(ignore_gap) {}

template<int verbose, typename SparseDistanceMatrixType>
void HammingDistWorkerImpl<verbose, SparseDistanceMatrixType>::operator()(std::size_t begin, std::size_t end) {
  double n = pss.num_seqs;
  double m = (n*n - 3.0*n + 2.0)/2.0;
  size_t my_aligned = 0;
  size_t my_considered = 0;
  size_t begin_i;
  if (begin == 0) {
    begin_i = 1;
  } else {
    begin_i = round(1.5 + 0.5*sqrt(9.0 + 8.0*((m*begin)/threads - 1.0)));
  }
  size_t end_i   = round(1.5 + 0.5*sqrt(9.0 + 8.0*((m*end)/threads - 1.0)));
  OPTIMOTU_DEBUG(1,
    << "HammingDistWorker thread " << begin
    << " entered; sequences [" << begin_i
    << ", "<< end_i << ")" << std::endl
  );

  OPTIMOTU_DEBUG(2,
    << "verbose = " << verbose
    << "\nSparseDistanceMatrixType = " << typeid(SparseDistanceMatrixType).name()
    << std::endl
  );

  // Buffers for the results
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

  // Buffers for the gap stats (only needed if the distance matrix is a
  // SparseDistanceMatrixGapstats)
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

  for (size_t i = begin_i; i < end_i; i++) {
    for (size_t j = 0; j < i; j++) {
      OPTIMOTU_DEBUG(
        2,
        << "considering " << i
        << " (seq: " << pss.packed_seq[i].size()
        << ", mask: " << pss.mask[i].size()
        << ") and " << j
        << " (seq: " << pss.packed_seq[j].size()
        << ", mask: " << pss.mask[j].size()
        << ")" << std::endl;
      );
      if constexpr (std::is_same_v<SparseDistanceMatrixType, SparseDistanceMatrixGapstats>) {
        auto[success, num_matches, d, num_ins, num_del, max_ins, max_del] =
          pss.success_score_and_dist_gap(i, j, min_overlap, ignore_gap);
        my_considered++;
        if (success) {
          my_aligned++;
          if (d <= dist_threshold) {
            my_seq1.push_back(j);
            my_seq2.push_back(i);
            my_score1.push_back(0);
            my_score2.push_back(num_matches);
            my_dist1.push_back(0);
            my_dist2.push_back(d);
            my_align_length.push_back((double)num_matches/(1.0 - d));
            my_n_insert.push_back(num_ins);
            my_n_delete.push_back(num_del);
            my_max_insert.push_back(max_ins);
            my_max_delete.push_back(max_del);
            if (my_seq1.size() == 1000) {
              sdm.append(my_seq1, my_seq2, my_score1, my_score2, my_dist1, my_dist2, my_align_length, my_n_insert, my_n_delete, my_max_insert, my_max_delete);
            }
          }
        }
      } else {
        auto[success, num_matches, d] =
          pss.success_score_and_dist(i, j, min_overlap, ignore_gap);
        my_considered++;
        if (success) {
          my_aligned++;
          if (d <= dist_threshold) {
            my_seq1.push_back(j);
            my_seq2.push_back(i);
            my_score1.push_back(0);
            my_score2.push_back(num_matches);
            my_dist1.push_back(0);
            my_dist2.push_back(d);
            if (my_seq1.size() == 1000) {
              sdm.append(my_seq1, my_seq2, my_score1, my_score2, my_dist1, my_dist2);
            }
          }
        }
      }
    }
  }
  if (my_seq1.size() > 0) {
    if constexpr (std::is_same_v<SparseDistanceMatrixType, SparseDistanceMatrixGapstats>) {
      sdm.append(my_seq1, my_seq2, my_score1, my_score2, my_dist1, my_dist2, my_align_length, my_n_insert, my_n_delete, my_max_insert, my_max_delete);
    } else {
      sdm.append(my_seq1, my_seq2, my_score1, my_score2, my_dist1, my_dist2);
    }
  }
  sdm.mutex.lock();
  _aligned += my_aligned;
  _prealigned += my_considered;
  sdm.mutex.unlock();
};

template<int verbose>
std::unique_ptr<HammingDistWorker> create_hamming_dist_worker(
  const std::vector<std::string> &seq,
  const double dist_threshold,
  const std::uint8_t threads,
  SparseDistanceMatrix &sdm,
  const int min_overlap,
  const bool ignore_gap
) {
  if (dynamic_cast<SparseDistanceMatrixGapstats*>(&sdm)) {
    return std::make_unique<HammingDistWorkerImpl<verbose, SparseDistanceMatrixGapstats>>(seq, dist_threshold, threads, sdm, min_overlap, ignore_gap);
  } else if (dynamic_cast<SparseDistanceMatrixCigar*>(&sdm)) {
    OPTIMOTU_STOP("Cigar is not implemented for Hamming distance");
  } else {
    return std::make_unique<HammingDistWorkerImpl<verbose, SparseDistanceMatrix>>(seq, dist_threshold, threads, sdm, min_overlap, ignore_gap);
  }
}

std::unique_ptr<HammingDistWorker> create_hamming_dist_worker(
  const std::vector<std::string> &seq,
  const double dist_threshold,
  const std::uint8_t threads,
  SparseDistanceMatrix &sdm,
  const int min_overlap,
  const bool ignore_gap,
  int verbose
) {
  if (verbose == 0) {
    return create_hamming_dist_worker<0>(seq, dist_threshold, threads, sdm, min_overlap, ignore_gap);
  } else if (verbose == 1) {
    return create_hamming_dist_worker<1>(seq, dist_threshold, threads, sdm, min_overlap, ignore_gap);
  } else {
    return create_hamming_dist_worker<2>(seq, dist_threshold, threads, sdm, min_overlap, ignore_gap);
  }
}

#endif //OPTIMOTU_R
