// SPDX-FileCopyrightText: 2025 Brendan Furneaux <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#include "DistClusterWorker.h"
#include "optimotu.h"
#include "wfa_fatal.h"

size_t DistClusterWorker::prealigned() {
  return _prealigned;
}

size_t DistClusterWorker::aligned() {
  return _aligned;
}

std::size_t DistClusterWorker::n_threads() {
  return threads;
}

std::size_t DistClusterWorker::n_tiles() {
  return pair_generators.size();
}

void DistClusterWorker::consume_gpu_tile(
  std::size_t,
  PairGenerator &,
  const GpuTileResult &
) {
  OPTIMOTU_STOP("GPU tiles are not implemented for this distance method");
}

void DistClusterWorker::run() {
  optimotu_wfa_fatal_clear();
  const std::size_t n = pair_generators.size();
  if (n == 0)
  {
    return;
  }
  if (wants_gpu() && gpu_runtime_available()) {
    const auto devices = gpu_selected_devices();
    if (!devices.empty()) {
      active_workers = std::max(threads, devices.size());
      if (n_active() > 1) {
        clust_algo.prepare_merge_slots(n);
      }
      gpu_run_tiles(
        pair_generators,
        n_active(),
        [this]() { return make_gpu_kernel(); },
        [this](
          std::size_t idx,
          PairGenerator &pg,
          const GpuTileResult &res
        ) {
          consume_gpu_tile(idx, pg, res);
        }
      );
      return;
    }
  }
  active_workers = threads;
  if (threads > 1)
  {
    clust_algo.prepare_merge_slots(n);
  }
  if (threads <= 1)
  {
    (*this)(0, n);
  }
  else
  {
    RcppParallel::parallelFor(0, n, *this, 1, threads);
  }
  optimotu_wfa_stop_if_fatal();
}

ClusterAlgorithm *DistClusterWorker::tile_algo(PairGenerator *pg)
{
  if (n_active() == 1)
  {
    return &clust_algo;
  }
  return clust_algo.make_child(pg);
}

void DistClusterWorker::finish_tile(ClusterAlgorithm *my_algo, std::size_t pg_index)
{
  if (n_active() == 1)
  {
    return;
  }
  if (!my_algo)
  {
    clust_algo.receive_mst(pg_index, {});
    return;
  }
  my_algo->set_merge_slot(pg_index);
  my_algo->finalize();
  my_algo->merge_into_parent();
  clust_algo.release_child(my_algo);
}
