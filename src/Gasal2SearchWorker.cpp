// SPDX-FileCopyrightText: 2026 Brendan Furneaux
// SPDX-License-Identifier: MIT

#include "optimotu.h"
#include "Gasal2SearchWorker.h"
#include "gasal2_edit.h"
#include "gasal2_trace.h"
#include "gpu_tile.h"
#ifdef OPTIMOTU_CUDA
#include "gasal2_gpu.h"
#endif

#include <algorithm>
#include <memory>

template<int verbose>
void Gasal2SearchWorkerImpl<verbose>::operator()(std::size_t, std::size_t) {
  OPTIMOTU_STOP("dist_gasal2() requires a CUDA device");
}

template<int verbose>
bool Gasal2SearchWorkerImpl<verbose>::try_run_gpu() {
  if (!gpu_runtime_available()) {
    OPTIMOTU_STOP("dist_gasal2() requires a CUDA device");
  }
  const auto devices = gpu_selected_devices();
  if (devices.empty()) {
    OPTIMOTU_STOP("dist_gasal2() requires a CUDA device");
  }
  const int nq = static_cast<int>(query.size());
  const int nr = static_cast<int>(ref.size());
  if (nq == 0 || nr == 0) {
    return true;
  }
  std::vector<std::string> packed = query;
  packed.insert(packed.end(), ref.begin(), ref.end());
  const int n_chunk = std::max(
    1,
    nq / std::max(1, static_cast<int>(devices.size() * 4))
  );
  std::vector<GpuTileSpec> specs;
  for (int q = 0; q < nq; q += n_chunk) {
    GpuTileSpec spec;
    spec.kind = GpuTileSpec::Kind::Rectangle;
    spec.i0 = q;
    spec.i1 = std::min(q + n_chunk, nq);
    spec.j0 = nq;
    spec.j1 = nq + nr;
    specs.push_back(spec);
  }
  gpu_run_specs(
    specs,
    std::max<std::size_t>(threads, devices.size()),
    [&]() {
#ifdef OPTIMOTU_CUDA
      return make_gasal2_gpu_kernel(
        packed, match, mismatch, gap_open, gap_extend
      );
#else
      OPTIMOTU_STOP("dist_gasal2() requires a CUDA build of optimotu");
      return nullptr;
#endif
    },
    [&](std::size_t, const GpuTileSpec &spec, const GpuTileResult &result) {
      Gasal2HostAlign host{match, mismatch, gap_open, gap_extend};
      std::size_t my_prealigned = 0;
      std::size_t my_aligned = 0;
      for (int i = spec.i0; i < spec.i1; ++i) {
        if (!hits[static_cast<std::size_t>(i)]) {
          hits[static_cast<std::size_t>(i)] = std::make_unique<SearchHit>();
        }
        SearchHit &hit = *hits[static_cast<std::size_t>(i)];
        for (int j = 0; j < nr; ++j) {
          const double max_dist =
            (hit.best_dist < threshold) ? hit.best_dist : threshold;
          const int score = gpu_result_score(spec, result, i, nq + j);
          const double l1 = static_cast<double>(query[static_cast<std::size_t>(i)].size());
          const double l2 = static_cast<double>(ref[static_cast<std::size_t>(j)].size());
          if (!host.pair_ok(dissimilarity, l1, l2, max_dist)) {
            continue;
          }
          ++my_prealigned;
          const int L = score_ref_length(l1, l2, AlignmentSpan::GLOBAL);
          const double d_score = gasal2_score_distance(
            score, L, match, mismatch, gap_open, gap_extend
          );
          double d = d_score;
          if (dissimilarity == Dissimilarity::IDENTITY) {
            if (!(d_score < max_dist)) {
              continue;
            }
            const bool q_longer =
              query[static_cast<std::size_t>(i)].size() >
              ref[static_cast<std::size_t>(j)].size();
            const std::string &s1 = q_longer
              ? ref[static_cast<std::size_t>(j)]
              : query[static_cast<std::size_t>(i)];
            const std::string &s2 = q_longer
              ? query[static_cast<std::size_t>(i)]
              : ref[static_cast<std::size_t>(j)];
            d = host.identity(s2, s1);
          }
          if (d == 1.0) {
            continue;
          }
          ++my_aligned;
          if (d < hit.best_dist) {
            hit.best_dist = d;
            hit.best_ref.clear();
            hit.best_ref.push_back(j);
          } else if (d == hit.best_dist) {
            hit.best_ref.push_back(j);
          }
        }
      }
      std::lock_guard<std::mutex> lock(mutex);
      _prealigned += my_prealigned;
      _aligned += my_aligned;
    }
  );
  return true;
}

template class Gasal2SearchWorkerImpl<0>;
template class Gasal2SearchWorkerImpl<1>;
template class Gasal2SearchWorkerImpl<2>;
template class Gasal2SearchWorkerImpl<3>;
template class Gasal2SearchWorkerImpl<4>;
