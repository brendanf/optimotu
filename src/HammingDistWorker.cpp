#include "HammingDistWorker.h"
#include "AllPairGenerator.h"
#include "BipartitePairGenerator.h"
#include "tile_pairs.h"
#include "hamming_kernel.h"
#include <algorithm>

#ifdef OPTIMOTU_R

namespace {

// Hamming visits ~1e9 pairs; checking R interrupt every pair dominates.
constexpr std::size_t kHammingInterruptBatch = 65536;

inline void hamming_maybe_interrupt(std::size_t &pairs_since_check) {
  if (++pairs_since_check >= kHammingInterruptBatch) {
    pairs_since_check = 0;
    RcppThread::checkUserInterrupt();
  }
}

}  // namespace

HammingDistWorker::HammingDistWorker(
  const std::vector<std::string> &seq,
  const double dist_threshold,
  DivisiblePairGenerator::Builder & pgb,
  SparseDistanceMatrix &sdm,
  const int min_overlap,
  const bool ignore_gap,
  int verbose,
  std::size_t worker_threads
) : DistWorker(seq, dist_threshold, pgb, sdm, verbose, worker_threads),
pss(seq),
min_overlap(min_overlap),
ignore_gap(ignore_gap) {}

template<int verbose, typename SparseDistanceMatrixType>
void HammingDistWorkerImpl<verbose, SparseDistanceMatrixType>::operator()(
  std::size_t begin, std::size_t end
) {
  (void)end;
  size_t my_aligned = 0;
  size_t my_considered = 0;
  std::size_t pairs_since_check = 0;
  OPTIMOTU_DEBUG(2,
    << "HammingDistWorker thread " << begin
    << " entered" << std::endl
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

  auto consider_plain = [&](
    std::size_t i,
    std::size_t j,
    int num_matches,
    int n_ok,
    double d
  ) {
    my_considered++;
    if (n_ok > 0 && n_ok >= min_overlap) {
      my_aligned++;
      if (d <= dist_threshold) {
        my_seq1.push_back(j);
        my_seq2.push_back(i);
        my_score1.push_back(0);
        my_score2.push_back(num_matches);
        my_dist1.push_back(0);
        my_dist2.push_back(d);
        if (my_seq1.size() == 1000) {
          sdm.append(
            my_seq1, my_seq2, my_score1, my_score2, my_dist1, my_dist2
          );
        }
      }
    }
    hamming_maybe_interrupt(pairs_since_check);
  };

  for_each_assigned_tile(
    pair_generators, begin, end,
    [&](std::size_t /*pg_index*/, PairGenerator &pg) {
      if constexpr (!std::is_same_v<
        SparseDistanceMatrixType, SparseDistanceMatrixGapstats
      >) {
        if ((ignore_gap || OPTIMOTU_HAMMING_PDISTB2_STRIPE) &&
            dynamic_cast<AllPairGenerator *>(&pg) != nullptr) {
          auto *diag = static_cast<AllPairGenerator *>(&pg);
          const int offset = static_cast<int>(diag->seq_offset());
          const int n_tile = static_cast<int>(diag->seq_count());
          std::vector<uint64_t> seq_major;
          std::vector<uint64_t> mask_major;
          pss.fill_word_major(offset, offset + n_tile, seq_major, mask_major);
          std::vector<double> dists(static_cast<std::size_t>(n_tile));
          std::vector<int> n_ok(static_cast<std::size_t>(n_tile));
          std::vector<int> n_match(static_cast<std::size_t>(n_tile));
          const int *tstart =
            ignore_gap ? nullptr : pss.start.data() + offset;
          const int *tend = ignore_gap ? nullptr : pss.end.data() + offset;
          for (int i = offset + 1; i < offset + n_tile; ++i) {
            const int n_j = i - offset;
            pss.score_query_vs_word_major(
              i, seq_major.data(), mask_major.data(), n_j, n_tile,
              min_overlap, dists.data(), n_ok.data(), n_match.data(),
              tstart, tend
            );
            for (int t = 0; t < n_j; ++t) {
              consider_plain(
                static_cast<std::size_t>(i),
                static_cast<std::size_t>(offset + t),
                n_match[static_cast<std::size_t>(t)],
                n_ok[static_cast<std::size_t>(t)],
                dists[static_cast<std::size_t>(t)]
              );
            }
          }
          return;
        }
        if ((ignore_gap || OPTIMOTU_HAMMING_PDISTB2_STRIPE) &&
            dynamic_cast<BipartitePairGenerator *>(&pg) != nullptr) {
          auto *bip = static_cast<BipartitePairGenerator *>(&pg);
          const int begin_i = static_cast<int>(bip->seq_begin_i());
          const int end_i = static_cast<int>(bip->seq_end_i());
          const int begin_j = static_cast<int>(bip->seq_begin_j());
          const int end_j = static_cast<int>(bip->seq_end_j());
          const int n_j = end_j - begin_j;
          std::vector<uint64_t> seq_major;
          std::vector<uint64_t> mask_major;
          pss.fill_word_major(begin_j, end_j, seq_major, mask_major);
          std::vector<double> dists(static_cast<std::size_t>(n_j));
          std::vector<int> n_ok(static_cast<std::size_t>(n_j));
          std::vector<int> n_match(static_cast<std::size_t>(n_j));
          const int *tstart =
            ignore_gap ? nullptr : pss.start.data() + begin_j;
          const int *tend = ignore_gap ? nullptr : pss.end.data() + begin_j;
          for (int i = begin_i; i < end_i; ++i) {
            pss.score_query_vs_word_major(
              i, seq_major.data(), mask_major.data(), n_j, n_j,
              min_overlap, dists.data(), n_ok.data(), n_match.data(),
              tstart, tend
            );
            for (int t = 0; t < n_j; ++t) {
              consider_plain(
                static_cast<std::size_t>(i),
                static_cast<std::size_t>(begin_j + t),
                n_match[static_cast<std::size_t>(t)],
                n_ok[static_cast<std::size_t>(t)],
                dists[static_cast<std::size_t>(t)]
              );
            }
          }
          return;
        }
      }
      for_each_tile_pair(pg, [&](
        std::size_t i0, std::size_t j0, std::size_t, std::size_t
      ) {
        const size_t i = i0;
        const size_t j = j0;
        OPTIMOTU_DEBUG(
          4,
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
            pss.success_score_and_dist_gap(
              static_cast<int>(i), static_cast<int>(j), min_overlap, ignore_gap
            );
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
                sdm.append(
                  my_seq1, my_seq2, my_score1, my_score2, my_dist1, my_dist2,
                  my_align_length, my_n_insert, my_n_delete, my_max_insert,
                  my_max_delete
                );
              }
            }
          }
        } else {
          auto[success, num_matches, d] =
            pss.success_score_and_dist(
              static_cast<int>(i), static_cast<int>(j), min_overlap, ignore_gap
            );
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
                sdm.append(
                  my_seq1, my_seq2, my_score1, my_score2, my_dist1, my_dist2
                );
              }
            }
          }
        }
        hamming_maybe_interrupt(pairs_since_check);
      });
    }
  );
  if (my_seq1.size() > 0) {
    if constexpr (std::is_same_v<SparseDistanceMatrixType, SparseDistanceMatrixGapstats>) {
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
  _prealigned += my_considered;
  sdm.mutex.unlock();
};

template<int verbose>
std::unique_ptr<HammingDistWorker> create_hamming_dist_worker(
  const std::vector<std::string> &seq,
  const double dist_threshold,
  DivisiblePairGenerator::Builder & pgb,
  SparseDistanceMatrix &sdm,
  const int min_overlap,
  const bool ignore_gap,
  std::size_t worker_threads
) {
  if (dynamic_cast<SparseDistanceMatrixGapstats*>(&sdm)) {
    return std::make_unique<HammingDistWorkerImpl<verbose, SparseDistanceMatrixGapstats>>(
      seq, dist_threshold, pgb, sdm, min_overlap, ignore_gap, worker_threads
    );
  } else if (dynamic_cast<SparseDistanceMatrixCigar*>(&sdm)) {
    OPTIMOTU_STOP("Cigar is not implemented for Hamming distance");
  } else {
    return std::make_unique<HammingDistWorkerImpl<verbose, SparseDistanceMatrix>>(
      seq, dist_threshold, pgb, sdm, min_overlap, ignore_gap, worker_threads
    );
  }
}

std::unique_ptr<HammingDistWorker> create_hamming_dist_worker(
  const std::vector<std::string> &seq,
  const double dist_threshold,
  DivisiblePairGenerator::Builder & pgb,
  SparseDistanceMatrix &sdm,
  const int min_overlap,
  const bool ignore_gap,
  int verbose,
  std::size_t worker_threads
) {
  int v = (verbose > 4) ? 4 : verbose;
  if (v == 0) {
    return create_hamming_dist_worker<0>(
      seq, dist_threshold, pgb, sdm, min_overlap, ignore_gap, worker_threads
    );
  } else if (v == 1) {
    return create_hamming_dist_worker<1>(
      seq, dist_threshold, pgb, sdm, min_overlap, ignore_gap, worker_threads
    );
  } else if (v == 2) {
    return create_hamming_dist_worker<2>(
      seq, dist_threshold, pgb, sdm, min_overlap, ignore_gap, worker_threads
    );
  } else if (v == 3) {
    return create_hamming_dist_worker<3>(
      seq, dist_threshold, pgb, sdm, min_overlap, ignore_gap, worker_threads
    );
  } else {
    return create_hamming_dist_worker<4>(
      seq, dist_threshold, pgb, sdm, min_overlap, ignore_gap, worker_threads
    );
  }
}

#endif //OPTIMOTU_R
