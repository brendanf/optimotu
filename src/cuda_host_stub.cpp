// SPDX-FileCopyrightText: 2026 Brendan Furneaux
// SPDX-License-Identifier: MIT

// nvcc emits a host stub for each __device__ function. The stub loads
// 1 into edi and calls exit(), which would terminate R if that stub
// ever ran. Device code itself lives in the fatbin, so the stub is
// not on the real call path. src/Makevars renames exit in the .cu
// objects only; exit in every other object is unchanged and still
// reported by R CMD check.

#include <Rcpp.h>

extern "C" void optimotu_nvcc_device_stub_exit(int status) {
  Rcpp::stop(
    "internal CUDA error: device function invoked on the host (status %d)",
    status
  );
}
