// SPDX-FileCopyrightText: 2026 Brendan Furneaux
// SPDX-License-Identifier: MIT

#ifdef OPTIMOTU_R
#include "Gasal2DistWorker.h"
#include "gasal2_edit.h"
#include "gasal2_trace.h"
#include "optimotu.h"
#include "tile_pairs.h"
#ifdef OPTIMOTU_CUDA
#include "gasal2_gpu.h"
#endif

Gasal2DistWorker::Gasal2DistWorker(
  const std::vector<std::string> &seq,
  const double dist_threshold,
  DivisiblePairGenerator::Builder &pgb,
  SparseDistanceMatrix &sdm,
  int match,
  int mismatch,
  int gap_open,
  int gap_extend,
  int verbose,
  std::size_t worker_threads,
  Dissimilarity dissimilarity
) : DistWorker(
      seq, dist_threshold, pgb, sdm, verbose, worker_threads, dissimilarity
    ),
    match(match),
    mismatch(mismatch),
    gap_open(gap_open),
    gap_extend(gap_extend) {}

bool Gasal2DistWorker::wants_gpu() const {
  return true;
}

std::unique_ptr<GpuDistanceKernel> Gasal2DistWorker::make_gpu_kernel() {
#ifdef OPTIMOTU_CUDA
  return make_gasal2_gpu_kernel(seq, match, mismatch, gap_open, gap_extend);
#else
  OPTIMOTU_STOP("dist_gasal2() requires a CUDA build of optimotu");
  return nullptr;
#endif
}

template<int verbose>
void Gasal2DistWorkerImpl<verbose>::operator()(std::size_t, std::size_t) {
  OPTIMOTU_STOP("dist_gasal2() requires a CUDA device");
}

template<int verbose>
void Gasal2DistWorkerImpl<verbose>::consume_gpu_tile(
  std::size_t,
  PairGenerator &pg,
  const GpuTileResult &result
) {
  const GpuTileSpec spec = gpu_tile_from_generator(pg);
  Gasal2HostAlign host{match, mismatch, gap_open, gap_extend};
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
  size_t my_prealigned = 0;
  size_t my_aligned = 0;
  auto flush_if_full = [&]() {
    if (my_seq1.size() != 1000) {
      return;
    }
    sdm.append(my_seq1, my_seq2, my_score1, my_score2, my_dist1, my_dist2);
  };
  for_each_tile_pair(pg, [&](
    std::size_t i0, std::size_t j0, std::size_t, std::size_t
  ) {
    const int score = gpu_result_score(
      spec, result, static_cast<int>(i0), static_cast<int>(j0)
    );
    const double l1 = static_cast<double>(seq[i0].size());
    const double l2 = static_cast<double>(seq[j0].size());
    if (!host.pair_ok(dissimilarity, l1, l2, dist_threshold)) {
      return;
    }
    ++my_prealigned;
    const int L = score_ref_length(l1, l2, AlignmentSpan::GLOBAL);
    const double d_score = gasal2_score_distance(
      score, L, match, mismatch, gap_open, gap_extend
    );
    double d = d_score;
    if (dissimilarity == Dissimilarity::IDENTITY) {
      if (!(d_score < dist_threshold)) {
        return;
      }
      const bool longer_j = seq[j0].size() > seq[i0].size();
      const std::string &query = longer_j ? seq[j0] : seq[i0];
      const std::string &target = longer_j ? seq[i0] : seq[j0];
      d = host.identity(query, target);
    }
    if (d > dist_threshold) {
      return;
    }
    ++my_aligned;
    my_seq1.push_back(j0);
    my_seq2.push_back(i0);
    my_score1.push_back(0);
    my_score2.push_back(score);
    my_dist1.push_back(0);
    my_dist2.push_back(d);
    flush_if_full();
  });
  if (!my_seq1.empty()) {
    sdm.append(my_seq1, my_seq2, my_score1, my_score2, my_dist1, my_dist2);
  }
  sdm.mutex.lock();
  _aligned += my_aligned;
  _prealigned += my_prealigned;
  sdm.mutex.unlock();
}

template class Gasal2DistWorkerImpl<0>;
template class Gasal2DistWorkerImpl<1>;
template class Gasal2DistWorkerImpl<2>;
template class Gasal2DistWorkerImpl<3>;
template class Gasal2DistWorkerImpl<4>;

std::unique_ptr<Gasal2DistWorker> create_gasal2_dist_worker(
  const std::vector<std::string> &seq,
  const double dist_threshold,
  DivisiblePairGenerator::Builder &pgb,
  SparseDistanceMatrix &sdm,
  int match,
  int mismatch,
  int gap_open,
  int gap_extend,
  int verbose,
  std::size_t worker_threads,
  Dissimilarity dissimilarity
) {
  if (dynamic_cast<SparseDistanceMatrixCigar *>(&sdm) ||
      dynamic_cast<SparseDistanceMatrixGapstats *>(&sdm)) {
    OPTIMOTU_STOP(
      "CIGAR and gapstats are not implemented for dist_gasal2()"
    );
  }
  const int v = (verbose > 4) ? 4 : verbose;
  if (v == 0) {
    return std::make_unique<Gasal2DistWorkerImpl<0>>(
      seq, dist_threshold, pgb, sdm, match, mismatch, gap_open, gap_extend,
      0, worker_threads, dissimilarity
    );
  }
  if (v == 1) {
    return std::make_unique<Gasal2DistWorkerImpl<1>>(
      seq, dist_threshold, pgb, sdm, match, mismatch, gap_open, gap_extend,
      1, worker_threads, dissimilarity
    );
  }
  if (v == 2) {
    return std::make_unique<Gasal2DistWorkerImpl<2>>(
      seq, dist_threshold, pgb, sdm, match, mismatch, gap_open, gap_extend,
      2, worker_threads, dissimilarity
    );
  }
  if (v == 3) {
    return std::make_unique<Gasal2DistWorkerImpl<3>>(
      seq, dist_threshold, pgb, sdm, match, mismatch, gap_open, gap_extend,
      3, worker_threads, dissimilarity
    );
  }
  return std::make_unique<Gasal2DistWorkerImpl<4>>(
    seq, dist_threshold, pgb, sdm, match, mismatch, gap_open, gap_extend,
    4, worker_threads, dissimilarity
  );
}

#endif
