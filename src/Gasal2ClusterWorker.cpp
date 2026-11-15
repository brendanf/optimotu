// SPDX-FileCopyrightText: 2026 Brendan Furneaux
// SPDX-License-Identifier: MIT

#include "Gasal2ClusterWorker.h"
#include "MappedClusterAlgorithm.h"
#include "gasal2_edit.h"
#include "gasal2_trace.h"
#include "optimotu.h"
#include "tile_pairs.h"
#ifdef OPTIMOTU_CUDA
#include "gasal2_gpu.h"
#endif

Gasal2ClusterWorker::Gasal2ClusterWorker(
  const SequenceSet &seq,
  ClusterAlgorithm &clust_algo,
  DivisiblePairGenerator::Builder &pgb,
  int match,
  int mismatch,
  int gap_open,
  int gap_extend,
  int verbose,
  std::size_t worker_threads,
  Dissimilarity dissimilarity
) : DistClusterWorker(
      seq, clust_algo, pgb, verbose, worker_threads, dissimilarity
    ),
    match(match),
    mismatch(mismatch),
    gap_open(gap_open),
    gap_extend(gap_extend) {}

bool Gasal2ClusterWorker::wants_gpu() const {
  return true;
}

std::unique_ptr<GpuDistanceKernel> Gasal2ClusterWorker::make_gpu_kernel() {
#ifdef OPTIMOTU_CUDA
  return make_gasal2_gpu_kernel(seq, match, mismatch, gap_open, gap_extend);
#else
  OPTIMOTU_STOP("dist_gasal2() requires a CUDA build of optimotu");
  return nullptr;
#endif
}

namespace {

void apply_gasal_pair(
  const SequenceSet &seq,
  ClusterAlgorithm *algo,
  bool use_local,
  Dissimilarity dissimilarity,
  Gasal2HostAlign &host,
  std::size_t i0,
  std::size_t j0,
  std::size_t li,
  std::size_t lj,
  int score,
  size_t *my_prealigned,
  size_t *my_aligned
) {
  const double l1 = static_cast<double>(seq[i0].size());
  const double l2 = static_cast<double>(seq[j0].size());
  const double threshold = use_local
    ? algo->max_relevant_local(static_cast<j_t>(li), static_cast<j_t>(lj))
    : algo->max_relevant(static_cast<j_t>(i0), static_cast<j_t>(j0));
  if (!host.pair_ok(dissimilarity, l1, l2, threshold)) {
    return;
  }
  ++(*my_prealigned);
  const int L = score_ref_length(l1, l2, AlignmentSpan::GLOBAL);
  const double d_score = gasal2_score_distance(
    score, L, host.match, host.mismatch, host.gap_open, host.gap_extend
  );
  double d = d_score;
  if (dissimilarity == Dissimilarity::IDENTITY) {
    if (!(d_score < threshold)) {
      return;
    }
    bool longer_j = seq[j0].size() > seq[i0].size();
    const SequenceView &query = seq[longer_j ? j0 : i0];
    const SequenceView &target = seq[longer_j ? i0 : j0];
    d = host.identity(query, target);
  }
  if (d < 1.0) {
    ++(*my_aligned);
  }
  if (d < threshold) {
    if (use_local) {
      algo->apply_local(static_cast<j_t>(li), static_cast<j_t>(lj), d);
    } else {
      (*algo)(static_cast<j_t>(i0), static_cast<j_t>(j0), d);
    }
  }
}

}  // namespace

template<int verbose>
void Gasal2SplitClusterWorker<verbose>::operator()(std::size_t, std::size_t) {
  OPTIMOTU_STOP("dist_gasal2() requires a CUDA device");
}

template<int verbose>
void Gasal2SplitClusterWorker<verbose>::consume_gpu_tile(
  std::size_t pg_index,
  PairGenerator &pg,
  const GpuTileResult &result
) {
  ClusterAlgorithm *my_algo = tile_algo(&pg);
  if (!my_algo) {
    finish_tile(nullptr, pg_index);
    return;
  }
  my_algo->set_merge_slot(pg_index);
  const GpuTileSpec spec = gpu_tile_from_generator(pg);
  const bool use_local =
    dynamic_cast<MappedClusterAlgorithm *>(my_algo) != nullptr;
  Gasal2HostAlign host{match, mismatch, gap_open, gap_extend};
  size_t my_prealigned = 0;
  size_t my_aligned = 0;
  for_each_tile_pair(pg, [&](
    std::size_t i0, std::size_t j0, std::size_t li, std::size_t lj
  ) {
    const int score = gpu_result_score(
      spec, result, static_cast<int>(i0), static_cast<int>(j0)
    );
    apply_gasal_pair(
      seq, my_algo, use_local, dissimilarity, host, i0, j0, li, lj, score,
      &my_prealigned, &my_aligned
    );
  });
  mutex.lock();
  _aligned += my_aligned;
  _prealigned += my_prealigned;
  mutex.unlock();
  finish_tile(my_algo, pg_index);
}

template<int verbose>
void Gasal2ConcurrentClusterWorker<verbose>::operator()(
  std::size_t,
  std::size_t
) {
  OPTIMOTU_STOP("dist_gasal2() requires a CUDA device");
}

template<int verbose>
void Gasal2ConcurrentClusterWorker<verbose>::consume_gpu_tile(
  std::size_t,
  PairGenerator &pg,
  const GpuTileResult &result
) {
  const GpuTileSpec spec = gpu_tile_from_generator(pg);
  Gasal2HostAlign host{match, mismatch, gap_open, gap_extend};
  size_t my_prealigned = 0;
  size_t my_aligned = 0;
  for_each_tile_pair(pg, [&](
    std::size_t i0, std::size_t j0, std::size_t, std::size_t
  ) {
    const int score = gpu_result_score(
      spec, result, static_cast<int>(i0), static_cast<int>(j0)
    );
    apply_gasal_pair(
      seq, &clust_algo, false, dissimilarity, host, i0, j0, 0, 0, score,
      &my_prealigned, &my_aligned
    );
  });
  mutex.lock();
  _aligned += my_aligned;
  _prealigned += my_prealigned;
  mutex.unlock();
}

template class Gasal2SplitClusterWorker<0>;
template class Gasal2SplitClusterWorker<1>;
template class Gasal2SplitClusterWorker<2>;
template class Gasal2SplitClusterWorker<3>;
template class Gasal2SplitClusterWorker<4>;
template class Gasal2ConcurrentClusterWorker<0>;
template class Gasal2ConcurrentClusterWorker<1>;
template class Gasal2ConcurrentClusterWorker<2>;
template class Gasal2ConcurrentClusterWorker<3>;
template class Gasal2ConcurrentClusterWorker<4>;
