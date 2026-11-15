// SPDX-FileCopyrightText: 2026 Brendan Furneaux
// <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_EVO_GPU_H
#define OPTIMOTU_EVO_GPU_H

#ifdef OPTIMOTU_CUDA

#include "EvoModel.h"
#include "PackedAaSequenceSet.h"
#include "gpu_tile.h"

#include <memory>

std::unique_ptr<GpuDistanceKernel> make_evo_gpu_kernel(
  const PackedAaSequenceSet &pss,
  const EvoModel &model
);

#endif

#endif
