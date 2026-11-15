// SPDX-FileCopyrightText: 2026 Brendan Furneaux
// <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

// Packed global affine scores in the GASAL2 style (score only, no
// traceback). Edit penalties (match 0, linear equal sub/indel) keep
// match cost 0 so the score is -edits. Affine / KSW-style penalties
// map match 0 to +1 like Ksw2Aligner.

#include "gasal2_gpu.h"
#include "gasal2_edit.h"
#include "optimotu.h"

#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

namespace {

constexpr int kMaxLen = 1024;

void check_cuda(cudaError_t st, const char *what) {
  if (st != cudaSuccess) {
    OPTIMOTU_STOP(std::string(what) + ": " + cudaGetErrorString(st));
  }
}

__device__ int encode_base(char c) {
  if (c == 'A' || c == 'a') return 0;
  if (c == 'C' || c == 'c') return 1;
  if (c == 'G' || c == 'g') return 2;
  if (c == 'T' || c == 't' || c == 'U' || c == 'u') return 3;
  return 4;
}

__device__ int nw_affine(
  const char *a,
  int na,
  const char *b,
  int nb,
  int match,
  int mismatch,
  int gap_open,
  int gap_extend
) {
  if (na > kMaxLen || nb > kMaxLen) {
    return -0x3fffffff;
  }
  if (na == 0 && nb == 0) {
    return 0;
  }
  int H[kMaxLen + 1];
  int E[kMaxLen + 1];
  const int qe = gap_open + gap_extend;
  H[0] = 0;
  E[0] = -0x3fffffff;
  for (int j = 1; j <= nb; ++j) {
    H[j] = -qe - (j - 1) * gap_extend;
    E[j] = -0x3fffffff;
  }
  for (int i = 1; i <= na; ++i) {
    int h_diag = H[0];
    H[0] = -qe - (i - 1) * gap_extend;
    int f = -0x3fffffff;
    const int ai = encode_base(a[i - 1]);
    for (int j = 1; j <= nb; ++j) {
      const int bj = encode_base(b[j - 1]);
      int s = (ai < 4 && bj < 4 && ai == bj) ? match : mismatch;
      if (ai == 4 || bj == 4) {
        s = mismatch;
      }
      f = max(H[j] - qe, f - gap_extend);
      E[j] = max(H[j - 1] - qe, E[j] - gap_extend);
      const int h = max(h_diag + s, max(E[j], f));
      h_diag = H[j];
      H[j] = h;
    }
  }
  return H[nb];
}

__device__ void pair_index_triangle(
  long long pair,
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

struct SeqDev {
  const char *bytes;
  const int *offset;
  const int *length;
};

__device__ int score_at(
  SeqDev seqs,
  int i,
  int j,
  int match,
  int mismatch,
  int gap_open,
  int gap_extend
) {
  const int na = seqs.length[i];
  const int nb = seqs.length[j];
  return nw_affine(
    seqs.bytes + seqs.offset[i],
    na,
    seqs.bytes + seqs.offset[j],
    nb,
    match,
    mismatch,
    gap_open,
    gap_extend
  );
}

__global__ void gasal_triangle(
  SeqDev seqs,
  int offset,
  int n,
  int match,
  int mismatch,
  int gap_open,
  int gap_extend,
  int *out
) {
  const long long n_pairs = static_cast<long long>(n) * (n - 1) / 2;
  const long long pair =
    static_cast<long long>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (pair >= n_pairs) {
    return;
  }
  int i = 0;
  int j = 0;
  pair_index_triangle(pair, &i, &j);
  out[pair] = score_at(
    seqs, offset + i, offset + j, match, mismatch, gap_open, gap_extend
  );
}

__global__ void gasal_rectangle(
  SeqDev seqs,
  int i0,
  int ni,
  int j0,
  int nj,
  int match,
  int mismatch,
  int gap_open,
  int gap_extend,
  int *out
) {
  const long long n_pairs = static_cast<long long>(ni) * nj;
  const long long pair =
    static_cast<long long>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (pair >= n_pairs) {
    return;
  }
  const int qi = i0 + static_cast<int>(pair / nj);
  const int tj = j0 + static_cast<int>(pair % nj);
  out[pair] = score_at(
    seqs, qi, tj, match, mismatch, gap_open, gap_extend
  );
}

void pack_sequences(
  const SequenceSet &seq,
  std::vector<char> *bytes,
  std::vector<int> *offset,
  std::vector<int> *length
) {
  offset->resize(seq.size());
  length->resize(seq.size());
  int pos = 0;
  for (std::size_t i = 0; i < seq.size(); ++i) {
    const int n = static_cast<int>(seq[i].size());
    if (n > kMaxLen) {
      OPTIMOTU_STOP(
        "dist_gasal2() sequences longer than 1024 bases are not supported"
      );
    }
    (*offset)[i] = pos;
    (*length)[i] = n;
    pos += n;
  }
  bytes->resize(static_cast<std::size_t>(pos));
  for (std::size_t i = 0; i < seq.size(); ++i) {
    std::memcpy(
      bytes->data() + (*offset)[i],
      seq[i].data(),
      static_cast<std::size_t>((*length)[i])
    );
  }
}

class Gasal2GpuKernel : public GpuDistanceKernel {
  std::vector<char> bytes;
  std::vector<int> offset;
  std::vector<int> length;
  int match;
  int mismatch;
  int gap_open;
  int gap_extend;
  char *d_bytes = nullptr;
  int *d_offset = nullptr;
  int *d_length = nullptr;
  SeqDev seqs{};
  bool bound = false;

public:
  Gasal2GpuKernel(
    const SequenceSet &seq,
    int match,
    int mismatch,
    int gap_open,
    int gap_extend
  ) : match(gasal2_gpu_match(match, mismatch, gap_open, gap_extend)),
      mismatch(gasal2_gpu_mismatch(mismatch)),
      gap_open(gap_open),
      gap_extend(gap_extend) {
    pack_sequences(seq, &bytes, &offset, &length);
  }

  ~Gasal2GpuKernel() override {
    cudaFree(d_bytes);
    cudaFree(d_offset);
    cudaFree(d_length);
  }

  GpuKernelKind kind() const override {
    return GpuKernelKind::Score;
  }

  void bind_device(int) override {
    if (bound) {
      return;
    }
    check_cuda(cudaMalloc(&d_bytes, bytes.size()), "gasal2 bytes");
    if (!bytes.empty()) {
      check_cuda(
        cudaMemcpy(
          d_bytes, bytes.data(), bytes.size(), cudaMemcpyHostToDevice
        ),
        "copy gasal2 bytes"
      );
    }
    check_cuda(
      cudaMalloc(&d_offset, offset.size() * sizeof(int)),
      "gasal2 offset"
    );
    check_cuda(
      cudaMemcpy(
        d_offset,
        offset.data(),
        offset.size() * sizeof(int),
        cudaMemcpyHostToDevice
      ),
      "copy gasal2 offset"
    );
    check_cuda(
      cudaMalloc(&d_length, length.size() * sizeof(int)),
      "gasal2 length"
    );
    check_cuda(
      cudaMemcpy(
        d_length,
        length.data(),
        length.size() * sizeof(int),
        cudaMemcpyHostToDevice
      ),
      "copy gasal2 length"
    );
    seqs.bytes = d_bytes;
    seqs.offset = d_offset;
    seqs.length = d_length;
    bound = true;
  }

  void fill_distances(const GpuTileSpec &, double *) override {
    OPTIMOTU_STOP("GASAL2 GPU kernel produces scores, not distances");
  }

  void fill_scores(const GpuTileSpec &tile, int *out) override {
    const std::size_t np = tile.n_pairs();
    if (np == 0) {
      return;
    }
    int *d_out = nullptr;
    check_cuda(cudaMalloc(&d_out, np * sizeof(int)), "gasal2 out");
    const int threads = 128;
    const int blocks = static_cast<int>((np + threads - 1) / threads);
    if (tile.kind == GpuTileSpec::Kind::Triangle) {
      gasal_triangle<<<blocks, threads>>>(
        seqs, tile.offset, tile.n, match, mismatch, gap_open, gap_extend, d_out
      );
    } else {
      gasal_rectangle<<<blocks, threads>>>(
        seqs,
        tile.i0,
        tile.i1 - tile.i0,
        tile.j0,
        tile.j1 - tile.j0,
        match,
        mismatch,
        gap_open,
        gap_extend,
        d_out
      );
    }
    check_cuda(cudaDeviceSynchronize(), "gasal2 kernel");
    check_cuda(
      cudaMemcpy(out, d_out, np * sizeof(int), cudaMemcpyDeviceToHost),
      "copy gasal2 scores"
    );
    cudaFree(d_out);
  }
};

}  // namespace

std::unique_ptr<GpuDistanceKernel> make_gasal2_gpu_kernel(
  const SequenceSet &seq,
  int match,
  int mismatch,
  int gap_open,
  int gap_extend
) {
  return std::make_unique<Gasal2GpuKernel>(
    seq, match, mismatch, gap_open, gap_extend
  );
}

std::unique_ptr<GpuDistanceKernel> make_gasal2_gpu_kernel(
  const std::vector<std::string> &seq,
  int match,
  int mismatch,
  int gap_open,
  int gap_extend
) {
  return make_gasal2_gpu_kernel(
    sequence_views_from_strings(seq), match, mismatch, gap_open, gap_extend
  );
}
