// SPDX-FileCopyrightText: 2026 Brendan Furneaux
// <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#include "optimotu.h"
#include "AaSearchWorker.h"
#include "gpu_tile.h"
#ifdef OPTIMOTU_CUDA
#include "evo_gpu.h"
#endif
#include <algorithm>
#include <memory>

AaSearchWorker::AaSearchWorker(
  const std::vector<std::string> &query,
  const std::vector<std::string> &ref,
  const double threshold,
  const std::uint8_t threads,
  const int min_overlap,
  const bool ignore_gaps,
  const float *weight_lut_in,
  const EvoModel *evo_model
) : SearchWorker(query, ref, threshold, threads),
    pss(query, ref),
    min_overlap(min_overlap),
    ignore_gaps(ignore_gaps),
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
  const int nq = static_cast<int>(query.size());
  const int nr = static_cast<int>(ref.size());
  ref_idx.resize(static_cast<std::size_t>(nr));
  for (int j = 0; j < nr; ++j) {
    ref_idx[static_cast<std::size_t>(j)] = nq + j;
  }
  pss.fill_pos_major(nq, nq + nr, ref_pos_major);
}

void AaSearchWorker::score_refs(
  int query_i,
  double *dist_out,
  int *n_ok_out
) const {
  const int nr = static_cast<int>(ref_idx.size());
  if (has_evo) {
    pss.score_query_evo(
      query_i,
      ref_pos_major.data(),
      ref_idx.data(),
      nr,
      nr,
      min_overlap,
      ignore_gaps,
      evo,
      dist_out,
      n_ok_out
    );
    return;
  }
  pss.score_query_vs_pos_major(
    query_i,
    ref_pos_major.data(),
    ref_idx.data(),
    nr,
    nr,
    min_overlap,
    ignore_gaps,
    weight_lut,
    threshold,
    true,
    dist_out,
    n_ok_out
  );
}

template<int verbose>
void AaSearchWorkerImpl<verbose>::operator()(std::size_t begin, std::size_t end) {
  const std::size_t begin_i = (begin * query.size()) / threads;
  const std::size_t end_i = (end * query.size()) / threads;
  const int nr = static_cast<int>(ref.size());
  std::size_t my_prealigned = 0;
  std::size_t my_aligned = 0;
  std::vector<double> dists(static_cast<std::size_t>(nr));
  std::vector<int> n_ok(static_cast<std::size_t>(nr));

  for (std::size_t i = begin_i; i < end_i; ++i) {
    if (!hits[i]) {
      hits[i] = std::make_unique<SearchHit>();
    }
    SearchHit &hit = *hits[i];
    score_refs(static_cast<int>(i), dists.data(), n_ok.data());
    for (int j = 0; j < nr; ++j) {
      ++my_prealigned;
      const double d = dists[static_cast<std::size_t>(j)];
      if (has_evo) {
        if (!(d < evo.t_max)) {
          continue;
        }
      } else if (d == 1.0) {
        continue;
      }
      ++my_aligned;
      const double max_dist =
        (hit.best_dist < threshold) ? hit.best_dist : threshold;
      if (d < max_dist) {
        hit.best_dist = d;
        hit.best_ref.clear();
        hit.best_ref.push_back(j);
      } else if (d == hit.best_dist) {
        hit.best_ref.push_back(j);
      }
    }
#ifdef OPTIMOTU_R
    RcppThread::checkUserInterrupt();
#endif
  }
  {
    std::lock_guard<std::mutex> lock(mutex);
    _prealigned += my_prealigned;
    _aligned += my_aligned;
  }
}

template<int verbose>
bool AaSearchWorkerImpl<verbose>::try_run_gpu() {
  if (!has_evo || !gpu_runtime_available()) {
    return false;
  }
  const auto devices = gpu_selected_devices();
  if (devices.empty()) {
    return false;
  }
  const int nq = static_cast<int>(query.size());
  const int nr = static_cast<int>(ref.size());
  if (nq == 0 || nr == 0) {
    return true;
  }
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
      return make_evo_gpu_kernel(pss, evo);
#else
      return nullptr;
#endif
    },
    [&](std::size_t, const GpuTileSpec &spec, const GpuTileResult &result) {
      std::size_t my_prealigned = 0;
      std::size_t my_aligned = 0;
      for (int i = spec.i0; i < spec.i1; ++i) {
        if (!hits[static_cast<std::size_t>(i)]) {
          hits[static_cast<std::size_t>(i)] = std::make_unique<SearchHit>();
        }
        SearchHit &hit = *hits[static_cast<std::size_t>(i)];
        for (int j = 0; j < nr; ++j) {
          const double d = gpu_result_distance(spec, result, i, nq + j);
          ++my_prealigned;
          if (!(d < evo.t_max)) {
            continue;
          }
          ++my_aligned;
          const double max_dist =
            (hit.best_dist < threshold) ? hit.best_dist : threshold;
          if (d < max_dist) {
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

template class AaSearchWorkerImpl<0>;
template class AaSearchWorkerImpl<1>;
template class AaSearchWorkerImpl<2>;
template class AaSearchWorkerImpl<3>;
template class AaSearchWorkerImpl<4>;
