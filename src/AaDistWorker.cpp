// SPDX-FileCopyrightText: 2026 Brendan Furneaux
// <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#include "optimotu.h"
#include "AaDistWorker.h"
#include "AllPairGenerator.h"
#include "BipartitePairGenerator.h"
#include "tile_pairs.h"
#include <algorithm>

#ifdef OPTIMOTU_R

namespace {

constexpr std::size_t kAaInterruptBatch = 65536;

inline void aa_maybe_interrupt(std::size_t &pairs_since_check) {
  if (++pairs_since_check >= kAaInterruptBatch) {
    pairs_since_check = 0;
    RcppThread::checkUserInterrupt();
  }
}

void flush_hits(
  SparseDistanceMatrix &sdm,
  std::vector<size_t> &seq1,
  std::vector<size_t> &seq2,
  std::vector<int> &score1,
  std::vector<int> &score2,
  std::vector<double> &dist1,
  std::vector<double> &dist2
) {
  if (!seq1.empty()) {
    sdm.append(seq1, seq2, score1, score2, dist1, dist2);
  }
}

void consider_pair(
  std::size_t i,
  std::size_t j,
  double d,
  int n_ok,
  double dist_threshold,
  bool return_all,
  std::vector<size_t> &seq1,
  std::vector<size_t> &seq2,
  std::vector<int> &score1,
  std::vector<int> &score2,
  std::vector<double> &dist1,
  std::vector<double> &dist2,
  SparseDistanceMatrix &sdm
) {
  if (!return_all && d > dist_threshold) {
    return;
  }
  seq1.push_back(j);
  seq2.push_back(i);
  score1.push_back(0);
  score2.push_back(n_ok);
  dist1.push_back(0);
  dist2.push_back(d);
  if (seq1.size() == 1000) {
    sdm.append(seq1, seq2, score1, score2, dist1, dist2);
  }
}

}  // namespace

AaDistWorker::AaDistWorker(
  const std::vector<std::string> &seq,
  const double dist_threshold,
  DivisiblePairGenerator::Builder &pgb,
  SparseDistanceMatrix &sdm,
  const int min_overlap,
  const bool ignore_gaps,
  const float *weight_lut_in,
  bool return_all,
  int verbose,
  std::size_t worker_threads,
  const EvoModel *evo_model
) : DistWorker(seq, dist_threshold, pgb, sdm, verbose, worker_threads),
    pss(seq),
    min_overlap(min_overlap),
    ignore_gaps(ignore_gaps),
    return_all(return_all),
    has_evo(evo_model != nullptr),
    evo(evo_model != nullptr ? *evo_model : EvoModel()) {
  if (weight_lut_in == nullptr) {
    weight_lut = nullptr;
  } else {
    lut_storage.assign(
      weight_lut_in,
      weight_lut_in + PackedAaSequenceSet::kLut * PackedAaSequenceSet::kLut
    );
    weight_lut = lut_storage.data();
  }
}

void AaDistWorker::score_targets(
  int query,
  const uint8_t *pos_major,
  const int *target_idx,
  int n_target,
  int n_cols,
  double dist_threshold_in,
  bool early_exit,
  double *dist_out,
  int *n_ok_out
) const {
  if (has_evo) {
    pss.score_query_evo(
      query, pos_major, target_idx, n_target, n_cols,
      min_overlap, ignore_gaps, evo, dist_out, n_ok_out
    );
    return;
  }
  pss.score_query_vs_pos_major(
    query, pos_major, target_idx, n_target, n_cols,
    min_overlap, ignore_gaps, weight_lut, dist_threshold_in, early_exit,
    dist_out, n_ok_out
  );
}

std::tuple<bool, int, double> AaDistWorker::score_pair(int i, int j) const {
  if (has_evo) {
    return pss.score_evo_pair(i, j, min_overlap, ignore_gaps, evo);
  }
  return pss.success_score_and_dist(
    i, j, min_overlap, ignore_gaps, weight_lut
  );
}

template<int verbose>
void AaDistWorkerImpl<verbose>::operator()(std::size_t begin, std::size_t end) {
  size_t my_aligned = 0;
  size_t my_considered = 0;
  std::size_t pairs_since_check = 0;
  std::vector<size_t> my_seq1;
  std::vector<size_t> my_seq2;
  std::vector<int> my_score1;
  std::vector<int> my_score2;
  std::vector<double> my_dist1;
  std::vector<double> my_dist2;
  my_seq1.reserve(1000);
  my_seq2.reserve(1000);
  my_score1.reserve(1000);
  my_score2.reserve(1000);
  my_dist1.reserve(1000);
  my_dist2.reserve(1000);

  for_each_assigned_tile(
    pair_generators, begin, end,
    [&](std::size_t, PairGenerator &pg) {
      if (auto *diag = dynamic_cast<AllPairGenerator *>(&pg)) {
        const int offset = static_cast<int>(diag->seq_offset());
        const int n_tile = static_cast<int>(diag->seq_count());
        const int end_i = offset + n_tile;
        std::vector<uint8_t> pos_major;
        pss.fill_pos_major(offset, end_i, pos_major);
        std::vector<int> target_idx(static_cast<std::size_t>(n_tile));
        for (int t = 0; t < n_tile; ++t) {
          target_idx[static_cast<std::size_t>(t)] = offset + t;
        }
        std::vector<double> dists(static_cast<std::size_t>(n_tile));
        std::vector<int> n_ok(static_cast<std::size_t>(n_tile));
        for (int i = offset + 1; i < end_i; ++i) {
          const int n_j = i - offset;
          score_targets(
            i,
            pos_major.data(),
            target_idx.data(),
            n_j,
            n_tile,
            dist_threshold,
            !return_all,
            dists.data(),
            n_ok.data()
          );
          for (int t = 0; t < n_j; ++t) {
            ++my_considered;
            if (n_ok[static_cast<std::size_t>(t)] > 0) {
              ++my_aligned;
            }
            consider_pair(
              static_cast<std::size_t>(i),
              static_cast<std::size_t>(offset + t),
              dists[static_cast<std::size_t>(t)],
              n_ok[static_cast<std::size_t>(t)],
              dist_threshold,
              return_all,
              my_seq1, my_seq2, my_score1, my_score2, my_dist1, my_dist2,
              sdm
            );
            aa_maybe_interrupt(pairs_since_check);
          }
        }
        return;
      }
      if (auto *bip = dynamic_cast<BipartitePairGenerator *>(&pg)) {
        const int begin_i = static_cast<int>(bip->seq_begin_i());
        const int end_i = static_cast<int>(bip->seq_end_i());
        const int begin_j = static_cast<int>(bip->seq_begin_j());
        const int end_j = static_cast<int>(bip->seq_end_j());
        const int n_j = end_j - begin_j;
        std::vector<uint8_t> pos_major;
        pss.fill_pos_major(begin_j, end_j, pos_major);
        std::vector<int> target_idx(static_cast<std::size_t>(n_j));
        for (int t = 0; t < n_j; ++t) {
          target_idx[static_cast<std::size_t>(t)] = begin_j + t;
        }
        std::vector<double> dists(static_cast<std::size_t>(n_j));
        std::vector<int> n_ok(static_cast<std::size_t>(n_j));
        for (int i = begin_i; i < end_i; ++i) {
          score_targets(
            i,
            pos_major.data(),
            target_idx.data(),
            n_j,
            n_j,
            dist_threshold,
            !return_all,
            dists.data(),
            n_ok.data()
          );
          for (int t = 0; t < n_j; ++t) {
            ++my_considered;
            if (n_ok[static_cast<std::size_t>(t)] > 0) {
              ++my_aligned;
            }
            consider_pair(
              static_cast<std::size_t>(i),
              static_cast<std::size_t>(begin_j + t),
              dists[static_cast<std::size_t>(t)],
              n_ok[static_cast<std::size_t>(t)],
              dist_threshold,
              return_all,
              my_seq1, my_seq2, my_score1, my_score2, my_dist1, my_dist2,
              sdm
            );
            aa_maybe_interrupt(pairs_since_check);
          }
        }
        return;
      }
      for_each_tile_pair(pg, [&](
        std::size_t i0, std::size_t j0, std::size_t, std::size_t
      ) {
        auto [success, n_ok, d] = score_pair(
          static_cast<int>(i0),
          static_cast<int>(j0)
        );
        ++my_considered;
        if (success) {
          ++my_aligned;
        }
        consider_pair(
          i0, j0, d, n_ok, dist_threshold, return_all,
          my_seq1, my_seq2, my_score1, my_score2, my_dist1, my_dist2, sdm
        );
        aa_maybe_interrupt(pairs_since_check);
      });
    }
  );
  flush_hits(sdm, my_seq1, my_seq2, my_score1, my_score2, my_dist1, my_dist2);
  sdm.mutex.lock();
  _aligned += my_aligned;
  _prealigned += my_considered;
  sdm.mutex.unlock();
}

template class AaDistWorkerImpl<0>;
template class AaDistWorkerImpl<1>;
template class AaDistWorkerImpl<2>;
template class AaDistWorkerImpl<3>;
template class AaDistWorkerImpl<4>;

std::unique_ptr<AaDistWorker> create_aa_dist_worker(
  const std::vector<std::string> &seq,
  const double dist_threshold,
  DivisiblePairGenerator::Builder &pgb,
  SparseDistanceMatrix &sdm,
  const int min_overlap,
  const bool ignore_gaps,
  const float *weight_lut,
  bool return_all,
  int verbose,
  std::size_t worker_threads,
  const EvoModel *evo_model
) {
  if (dynamic_cast<SparseDistanceMatrixCigar *>(&sdm)) {
    OPTIMOTU_STOP("Cigar is not implemented for amino-acid distance");
  }
  if (dynamic_cast<SparseDistanceMatrixGapstats *>(&sdm)) {
    OPTIMOTU_STOP("gapstats are not implemented for amino-acid distance");
  }
  int v = (verbose > 4) ? 4 : verbose;
  if (v == 0) {
    return std::make_unique<AaDistWorkerImpl<0>>(
      seq, dist_threshold, pgb, sdm, min_overlap, ignore_gaps, weight_lut,
      return_all, 0, worker_threads, evo_model
    );
  }
  if (v == 1) {
    return std::make_unique<AaDistWorkerImpl<1>>(
      seq, dist_threshold, pgb, sdm, min_overlap, ignore_gaps, weight_lut,
      return_all, 1, worker_threads, evo_model
    );
  }
  if (v == 2) {
    return std::make_unique<AaDistWorkerImpl<2>>(
      seq, dist_threshold, pgb, sdm, min_overlap, ignore_gaps, weight_lut,
      return_all, 2, worker_threads, evo_model
    );
  }
  if (v == 3) {
    return std::make_unique<AaDistWorkerImpl<3>>(
      seq, dist_threshold, pgb, sdm, min_overlap, ignore_gaps, weight_lut,
      return_all, 3, worker_threads, evo_model
    );
  }
  return std::make_unique<AaDistWorkerImpl<4>>(
    seq, dist_threshold, pgb, sdm, min_overlap, ignore_gaps, weight_lut,
      return_all, 4, worker_threads, evo_model
  );
}

#endif
