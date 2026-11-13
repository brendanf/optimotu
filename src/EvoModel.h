// SPDX-FileCopyrightText: 2026 Brendan Furneaux
// <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_EVOMODEL_H
#define OPTIMOTU_EVOMODEL_H

#include <array>
#include <vector>

// Reversible amino-acid model prepared for an ML distance scan.
// logp layout is ((query_code * 32 + target_code) * n_grid + g), so the
// grid for one residue pair is contiguous. evec and evec_inv are
// column-major, matching R's as.numeric(matrix).
struct EvoModel {
  static constexpr int kState = 20;
  static constexpr int kCode = 32;

  double t_max = 5.0;
  std::vector<double> t_grid;
  std::vector<float> logp;
  std::vector<double> lambda;
  std::vector<double> evec;
  std::vector<double> evec_inv;
  std::array<int, kCode> code_to_state{};
  std::array<int, kCode> is_state{};

  EvoModel() {
    code_to_state.fill(-1);
    is_state.fill(0);
  }

  int n_grid() const {
    return static_cast<int>(t_grid.size());
  }

  double V_ik(int i, int k) const {
    return evec[static_cast<std::size_t>(k) * kState + static_cast<std::size_t>(i)];
  }

  double Vinv_kj(int k, int j) const {
    return evec_inv[
      static_cast<std::size_t>(j) * kState + static_cast<std::size_t>(k)
    ];
  }

  const float *logp_row(int query_code, int target_code) const {
    return logp.data() +
      ((static_cast<std::size_t>(query_code) * kCode +
        static_cast<std::size_t>(target_code)) *
       static_cast<std::size_t>(n_grid()));
  }

  float logp_at(int g, int query_code, int target_code) const {
    return logp_row(query_code, target_code)[g];
  }
};

#endif
