// SPDX-FileCopyrightText: 2026 Brendan Furneaux
// <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_GPU_TILE_H
#define OPTIMOTU_GPU_TILE_H

#include "PairGenerator.h"

#include <cstddef>
#include <functional>
#include <memory>
#include <string>
#include <vector>

// Dense all-pairs tile for a GPU kernel. Triangle is a diagonal
// AllPairGenerator (j < i inside [offset, offset + n)). Rectangle is a
// BipartitePairGenerator.
struct GpuTileSpec {
  enum class Kind { Triangle, Rectangle };
  Kind kind = Kind::Triangle;
  int offset = 0;
  int n = 0;
  int i0 = 0;
  int i1 = 0;
  int j0 = 0;
  int j1 = 0;

  std::size_t n_pairs() const;
};

enum class GpuKernelKind { Distance, Score };

struct GpuTileResult {
  GpuKernelKind kind = GpuKernelKind::Distance;
  std::vector<double> distances;
  std::vector<int> scores;
};

GpuTileSpec gpu_tile_from_generator(const PairGenerator &pg);
double gpu_result_distance(
  const GpuTileSpec &spec,
  const GpuTileResult &result,
  int i0,
  int j0
);
int gpu_result_score(
  const GpuTileSpec &spec,
  const GpuTileResult &result,
  int i0,
  int j0
);

class GpuDistanceKernel {
public:
  virtual ~GpuDistanceKernel() = default;
  virtual GpuKernelKind kind() const = 0;
  virtual void bind_device(int device) = 0;
  virtual void fill_distances(const GpuTileSpec &tile, double *out) = 0;
  virtual void fill_scores(const GpuTileSpec &tile, int *out) = 0;
};

bool gpu_runtime_available();
std::vector<int> gpu_selected_devices();
std::string gpu_device_name(int device);
int gpu_device_count();

// One pump thread per device. Workers pull finished tiles. make_kernel
// is called on the pump after cudaSetDevice so device pointers stay
// local to that GPU.
void gpu_run_specs(
  const std::vector<GpuTileSpec> &specs,
  std::size_t threads,
  const std::function<std::unique_ptr<GpuDistanceKernel>()> &make_kernel,
  const std::function<
    void(std::size_t, const GpuTileSpec &, const GpuTileResult &)> &consume
);

void gpu_run_tiles(
  const std::vector<std::unique_ptr<PairGenerator>> &generators,
  std::size_t threads,
  const std::function<std::unique_ptr<GpuDistanceKernel>()> &make_kernel,
  const std::function<
    void(std::size_t, PairGenerator &, const GpuTileResult &)> &consume
);

#endif
