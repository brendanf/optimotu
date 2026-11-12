// SPDX-FileCopyrightText: 2026 Brendan Furneaux
// <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_AACONFIG_H
#define OPTIMOTU_AACONFIG_H

#ifdef OPTIMOTU_R

#include "PackedAaSequenceSet.h"
#include "config.h"
#include <vector>

inline const float *aa_weight_lut_from_config(
  Rcpp::List dist_config,
  std::vector<float> &storage
) {
  if (!dist_config.containsElementNamed("weight_lut")) {
    return nullptr;
  }
  Rcpp::NumericMatrix m = dist_config["weight_lut"];
  if (m.nrow() != PackedAaSequenceSet::kLut ||
      m.ncol() != PackedAaSequenceSet::kLut) {
    OPTIMOTU_STOP(
      "dist_aa weight_lut must be a 32 x 32 numeric matrix"
    );
  }
  storage.resize(
    static_cast<std::size_t>(PackedAaSequenceSet::kLut) *
    static_cast<std::size_t>(PackedAaSequenceSet::kLut)
  );
  for (int a = 0; a < PackedAaSequenceSet::kLut; ++a) {
    for (int b = 0; b < PackedAaSequenceSet::kLut; ++b) {
      storage[
        static_cast<std::size_t>(a) * PackedAaSequenceSet::kLut +
        static_cast<std::size_t>(b)
      ] = static_cast<float>(m(a, b));
    }
  }
  return storage.data();
}

#endif

#endif
