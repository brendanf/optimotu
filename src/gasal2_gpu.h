// SPDX-FileCopyrightText: 2026 Brendan Furneaux
// <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_GASAL2_GPU_H
#define OPTIMOTU_GASAL2_GPU_H

#ifdef OPTIMOTU_CUDA

#include "SequenceView.h"
#include "gpu_tile.h"

#include <memory>
#include <vector>

std::unique_ptr<GpuDistanceKernel> make_gasal2_gpu_kernel(
  const SequenceSet &seq,
  int match,
  int mismatch,
  int gap_open,
  int gap_extend
);

std::unique_ptr<GpuDistanceKernel> make_gasal2_gpu_kernel(
  const std::vector<std::string> &seq,
  int match,
  int mismatch,
  int gap_open,
  int gap_extend
);

#endif

#endif
