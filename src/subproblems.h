// SPDX-FileCopyrightText: 2026 Brendan Furneaux
// <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_SUBPROBLEMS_H_INCLUDED
#define OPTIMOTU_SUBPROBLEMS_H_INCLUDED

#include <algorithm>
#include <cmath>
#include <string>

// Hamming / USEARCH: square nearest 4 * threads.
// Alignment methods: square in [threads, 4 * threads] with the best
// equal-tile utilization S / (threads * ceil(S / threads)); ties keep
// the smaller square.
inline bool optimotu_dist_method_is_slow(const std::string &method) {
  return method == "wfa2" || method == "edlib" || method == "ksw2" ||
         method == "hybrid";
}

inline int optimotu_square_subproblems(int n) {
  if (n <= 1) {
    return 1;
  }
  const double n_blocks = std::ceil(std::sqrt(static_cast<double>(n)));
  const int nb = static_cast<int>(n_blocks);
  return nb * nb;
}

inline int optimotu_default_subproblems(
    int threads,
    const std::string &method) {
  if (threads <= 1) {
    return 1;
  }
  const int max_s = 4 * threads;
  if (!optimotu_dist_method_is_slow(method)) {
    const double target = 4.0 * static_cast<double>(threads);
    int n_blocks = static_cast<int>(std::lround(std::sqrt(target)));
    if (n_blocks < 1) {
      n_blocks = 1;
    }
    int s = n_blocks * n_blocks;
    if (s < threads) {
      return optimotu_square_subproblems(threads);
    }
    return s;
  }
  int best_s = optimotu_square_subproblems(threads);
  double best_util = -1.0;
  int n_blocks = static_cast<int>(
      std::ceil(std::sqrt(static_cast<double>(threads))));
  for (;;) {
    const int s = n_blocks * n_blocks;
    if (s > max_s) {
      break;
    }
    const int waves = (s + threads - 1) / threads;
    const double util = static_cast<double>(s) /
                        (static_cast<double>(threads) * waves);
    if (util > best_util + 1e-12 ||
        (std::abs(util - best_util) <= 1e-12 && s < best_s)) {
      best_util = util;
      best_s = s;
    }
    ++n_blocks;
  }
  return best_s;
}

inline int optimotu_resolve_subproblems(
    int threads,
    const std::string &method,
    int requested) {
  if (requested > 0) {
    return std::max(threads, requested);
  }
  return optimotu_default_subproblems(threads, method);
}

#endif
