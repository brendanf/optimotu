// SPDX-FileCopyrightText: 2026 Brendan Furneaux
// <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#include "evo_gpu.h"
#include "optimotu.h"

#include <cuda_runtime.h>

#include <cmath>
#include <cstring>
#include <string>
#include <vector>

namespace {

constexpr int kMaxGrid = 64;

void check_cuda(cudaError_t st, const char *what) {
  if (st != cudaSuccess) {
    OPTIMOTU_STOP(std::string(what) + ": " + cudaGetErrorString(st));
  }
}

struct GpuModel {
  const float *logp;
  const double *t_grid;
  const double *lambda;
  const double *v;
  const double *u;
  const int *code_to_state;
  const int *is_state;
  int n_grid;
  double t_max;
};

__device__ void dot_p(
  const double *v,
  const double *u,
  const double *e,
  const double *de,
  const double *d2e,
  double *p_out,
  double *dp_out,
  double *d2p_out
) {
  double p = 0.0;
  double dp = 0.0;
  double d2p = 0.0;
  for (int k = 0; k < EvoModel::kState; ++k) {
    const double w = v[k] * u[k];
    p += w * e[k];
    dp += w * de[k];
    d2p += w * d2e[k];
  }
  *p_out = p;
  *dp_out = dp;
  *d2p_out = d2p;
}

__device__ void hoist_exp(
  double tt,
  const double *lambda,
  double *e,
  double *de,
  double *d2e
) {
  for (int k = 0; k < EvoModel::kState; ++k) {
    const double lam = lambda[k];
    e[k] = exp(lam * tt);
    de[k] = e[k] * lam;
    d2e[k] = de[k] * lam;
  }
}

__device__ double spectral_loglik(
  const unsigned short *counts,
  const double *v,
  const double *u,
  const double *e
) {
  double sum = 0.0;
  for (int i = 0; i < EvoModel::kState; ++i) {
    for (int j = 0; j < EvoModel::kState; ++j) {
      const int c = counts[i * EvoModel::kState + j];
      if (c == 0) {
        continue;
      }
      double p = 0.0;
      for (int k = 0; k < EvoModel::kState; ++k) {
        p += v[i * EvoModel::kState + k] *
          u[j * EvoModel::kState + k] * e[k];
      }
      if (!(p > 1e-300)) {
        p = 1e-300;
      }
      sum += static_cast<double>(c) * log(p);
    }
  }
  return sum;
}

__device__ double score_pair(
  const uint8_t *qa,
  const uint8_t *tb,
  int alen,
  GpuModel model
) {
  double acc[kMaxGrid];
  for (int g = 0; g < model.n_grid; ++g) {
    acc[g] = 0.0;
  }
  unsigned short counts[EvoModel::kState * EvoModel::kState];
  for (int i = 0; i < EvoModel::kState * EvoModel::kState; ++i) {
    counts[i] = 0;
  }
  int n_ok = 0;
  for (int s = 0; s < alen; ++s) {
    const uint8_t a = qa[s];
    const uint8_t b = tb[s];
    if (a >= EvoModel::kCode || b >= EvoModel::kCode ||
        model.is_state[a] == 0 || model.is_state[b] == 0) {
      continue;
    }
    const int si = model.code_to_state[a];
    const int sj = model.code_to_state[b];
    ++counts[si * EvoModel::kState + sj];
    ++n_ok;
    const float *row = model.logp +
      (static_cast<long long>(a) * EvoModel::kCode + b) * model.n_grid;
    for (int g = 0; g < model.n_grid; ++g) {
      acc[g] += static_cast<double>(row[g]);
    }
  }
  if (n_ok == 0) {
    return model.t_max;
  }
  int best = 0;
  double best_ll = -1.0e300;
  for (int g = 0; g < model.n_grid; ++g) {
    if (acc[g] > best_ll) {
      best_ll = acc[g];
      best = g;
    }
  }
  if (best == 0 || best == model.n_grid - 1) {
    return model.t_grid[best];
  }

  double t = model.t_grid[best];
  const double t_lo = model.t_grid[best - 1];
  const double t_hi = model.t_grid[best + 1];
  double e[EvoModel::kState];
  double de[EvoModel::kState];
  double d2e[EvoModel::kState];
  hoist_exp(t, model.lambda, e, de, d2e);
  const double grid_ll = spectral_loglik(counts, model.v, model.u, e);
  bool moving = true;
  for (int iter = 0; iter < 8 && moving; ++iter) {
    hoist_exp(t, model.lambda, e, de, d2e);
    double g = 0.0;
    double hess = 0.0;
    for (int i = 0; i < EvoModel::kState; ++i) {
      for (int j = 0; j < EvoModel::kState; ++j) {
        const int c = counts[i * EvoModel::kState + j];
        if (c == 0) {
          continue;
        }
        double p = 0.0;
        double dp = 0.0;
        double d2p = 0.0;
        dot_p(
          model.v + i * EvoModel::kState,
          model.u + j * EvoModel::kState,
          e,
          de,
          d2e,
          &p,
          &dp,
          &d2p
        );
        if (!(p > 1e-300)) {
          p = 1e-300;
        }
        const double inv = 1.0 / p;
        const double r = dp * inv;
        const double cd = static_cast<double>(c);
        g += cd * r;
        hess += cd * (d2p * inv - r * r);
      }
    }
    if (!isfinite(g) || !isfinite(hess) || !(hess < -1e-8)) {
      break;
    }
    double t2 = t - g / hess;
    if (t2 <= t_lo) {
      t2 = 0.5 * (t + t_lo);
    }
    if (t2 >= t_hi) {
      t2 = 0.5 * (t + t_hi);
    }
    if (!isfinite(t2)) {
      break;
    }
    const double step = fabs(t2 - t);
    t = t2;
    if (step < 1e-8 * (1.0 + t2)) {
      moving = false;
    }
  }
  hoist_exp(t, model.lambda, e, de, d2e);
  if (spectral_loglik(counts, model.v, model.u, e) + 1e-8 >= grid_ll) {
    return t;
  }
  return model.t_grid[best];
}

__device__ void pair_index_triangle(
  long long pair,
  int n,
  int *i,
  int *j
) {
  long long row = static_cast<long long>(
    (1.0 + sqrt(1.0 + 8.0 * static_cast<double>(pair))) * 0.5
  );
  while (row * (row - 1) / 2 > pair) {
    --row;
  }
  while ((row + 1) * row / 2 <= pair) {
    ++row;
  }
  *i = static_cast<int>(row);
  *j = static_cast<int>(pair - row * (row - 1) / 2);
}

__global__ void evo_triangle(
  const uint8_t *codes,
  int alen,
  int offset,
  int n,
  GpuModel model,
  double *out
) {
  const long long n_pairs =
    static_cast<long long>(n) * (n - 1) / 2;
  const long long pair =
    static_cast<long long>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (pair >= n_pairs) {
    return;
  }
  int i = 0;
  int j = 0;
  pair_index_triangle(pair, n, &i, &j);
  const uint8_t *qa =
    codes + static_cast<long long>(offset + i) * alen;
  const uint8_t *tb =
    codes + static_cast<long long>(offset + j) * alen;
  out[pair] = score_pair(qa, tb, alen, model);
}

__global__ void evo_rectangle(
  const uint8_t *codes,
  int alen,
  int i0,
  int ni,
  int j0,
  int nj,
  GpuModel model,
  double *out
) {
  const long long n_pairs = static_cast<long long>(ni) * nj;
  const long long pair =
    static_cast<long long>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (pair >= n_pairs) {
    return;
  }
  const int qi = i0 + static_cast<int>(pair / nj);
  const int tj = j0 + static_cast<int>(pair % nj);
  const uint8_t *qa = codes + static_cast<long long>(qi) * alen;
  const uint8_t *tb = codes + static_cast<long long>(tj) * alen;
  out[pair] = score_pair(qa, tb, alen, model);
}

class EvoGpuKernel : public GpuDistanceKernel {
  const PackedAaSequenceSet &pss;
  const EvoModel &model;
  uint8_t *d_codes = nullptr;
  float *d_logp = nullptr;
  double *d_t_grid = nullptr;
  double *d_lambda = nullptr;
  double *d_v = nullptr;
  double *d_u = nullptr;
  int *d_code_to_state = nullptr;
  int *d_is_state = nullptr;
  GpuModel gpu{};
  bool bound = false;

public:
  EvoGpuKernel(const PackedAaSequenceSet &pss, const EvoModel &model)
    : pss(pss), model(model) {
    if (model.n_grid() > kMaxGrid) {
      OPTIMOTU_STOP("evo GPU kernel grid is longer than 64");
    }
  }

  ~EvoGpuKernel() override {
    cudaFree(d_codes);
    cudaFree(d_logp);
    cudaFree(d_t_grid);
    cudaFree(d_lambda);
    cudaFree(d_v);
    cudaFree(d_u);
    cudaFree(d_code_to_state);
    cudaFree(d_is_state);
  }

  GpuKernelKind kind() const override {
    return GpuKernelKind::Distance;
  }

  void bind_device(int) override {
    if (bound) {
      return;
    }
    const std::size_t n_code =
      static_cast<std::size_t>(pss.num_seqs) *
      static_cast<std::size_t>(pss.alen);
    check_cuda(cudaMalloc(&d_codes, n_code), "evo codes");
    check_cuda(
      cudaMemcpy(d_codes, pss.codes.data(), n_code, cudaMemcpyHostToDevice),
      "copy evo codes"
    );
    check_cuda(
      cudaMalloc(&d_logp, model.logp.size() * sizeof(float)),
      "evo logp"
    );
    check_cuda(
      cudaMemcpy(
        d_logp,
        model.logp.data(),
        model.logp.size() * sizeof(float),
        cudaMemcpyHostToDevice
      ),
      "copy evo logp"
    );
    check_cuda(
      cudaMalloc(&d_t_grid, model.t_grid.size() * sizeof(double)),
      "evo t_grid"
    );
    check_cuda(
      cudaMemcpy(
        d_t_grid,
        model.t_grid.data(),
        model.t_grid.size() * sizeof(double),
        cudaMemcpyHostToDevice
      ),
      "copy evo t_grid"
    );
    check_cuda(
      cudaMalloc(&d_lambda, model.lambda.size() * sizeof(double)),
      "evo lambda"
    );
    check_cuda(
      cudaMemcpy(
        d_lambda,
        model.lambda.data(),
        model.lambda.size() * sizeof(double),
        cudaMemcpyHostToDevice
      ),
      "copy evo lambda"
    );
    std::vector<double> v(EvoModel::kState * EvoModel::kState);
    std::vector<double> u(EvoModel::kState * EvoModel::kState);
    for (int i = 0; i < EvoModel::kState; ++i) {
      for (int k = 0; k < EvoModel::kState; ++k) {
        v[static_cast<std::size_t>(i) * EvoModel::kState + k] =
          model.V_ik(i, k);
        u[static_cast<std::size_t>(i) * EvoModel::kState + k] =
          model.Vinv_kj(k, i);
      }
    }
    check_cuda(cudaMalloc(&d_v, v.size() * sizeof(double)), "evo v");
    check_cuda(cudaMalloc(&d_u, u.size() * sizeof(double)), "evo u");
    check_cuda(
      cudaMemcpy(
        d_v, v.data(), v.size() * sizeof(double), cudaMemcpyHostToDevice
      ),
      "copy evo v"
    );
    check_cuda(
      cudaMemcpy(
        d_u, u.data(), u.size() * sizeof(double), cudaMemcpyHostToDevice
      ),
      "copy evo u"
    );
    check_cuda(
      cudaMalloc(&d_code_to_state, EvoModel::kCode * sizeof(int)),
      "evo code_to_state"
    );
    check_cuda(
      cudaMalloc(&d_is_state, EvoModel::kCode * sizeof(int)),
      "evo is_state"
    );
    check_cuda(
      cudaMemcpy(
        d_code_to_state,
        model.code_to_state.data(),
        EvoModel::kCode * sizeof(int),
        cudaMemcpyHostToDevice
      ),
      "copy evo code_to_state"
    );
    check_cuda(
      cudaMemcpy(
        d_is_state,
        model.is_state.data(),
        EvoModel::kCode * sizeof(int),
        cudaMemcpyHostToDevice
      ),
      "copy evo is_state"
    );
    gpu.logp = d_logp;
    gpu.t_grid = d_t_grid;
    gpu.lambda = d_lambda;
    gpu.v = d_v;
    gpu.u = d_u;
    gpu.code_to_state = d_code_to_state;
    gpu.is_state = d_is_state;
    gpu.n_grid = model.n_grid();
    gpu.t_max = model.t_max;
    bound = true;
  }

  void fill_distances(const GpuTileSpec &tile, double *out) override {
    const std::size_t np = tile.n_pairs();
    if (np == 0) {
      return;
    }
    double *d_out = nullptr;
    check_cuda(cudaMalloc(&d_out, np * sizeof(double)), "evo out");
    const int threads = 128;
    const int blocks = static_cast<int>((np + threads - 1) / threads);
    if (tile.kind == GpuTileSpec::Kind::Triangle) {
      evo_triangle<<<blocks, threads>>>(
        d_codes, pss.alen, tile.offset, tile.n, gpu, d_out
      );
    } else {
      evo_rectangle<<<blocks, threads>>>(
        d_codes,
        pss.alen,
        tile.i0,
        tile.i1 - tile.i0,
        tile.j0,
        tile.j1 - tile.j0,
        gpu,
        d_out
      );
    }
    check_cuda(cudaDeviceSynchronize(), "evo kernel");
    check_cuda(
      cudaMemcpy(out, d_out, np * sizeof(double), cudaMemcpyDeviceToHost),
      "copy evo distances"
    );
    cudaFree(d_out);
  }

  void fill_scores(const GpuTileSpec &, int *) override {
    OPTIMOTU_STOP("evo GPU kernel does not produce alignment scores");
  }
};

}  // namespace

std::unique_ptr<GpuDistanceKernel> make_evo_gpu_kernel(
  const PackedAaSequenceSet &pss,
  const EvoModel &model
) {
  return std::make_unique<EvoGpuKernel>(pss, model);
}
