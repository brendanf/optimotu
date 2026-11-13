// SPDX-FileCopyrightText: 2026 Brendan Furneaux
// <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_EVOCONFIG_H
#define OPTIMOTU_EVOCONFIG_H

#ifdef OPTIMOTU_R

#include "EvoModel.h"
#include "config.h"
#include <string>

inline EvoModel evo_model_from_config(Rcpp::List dist_config) {
  EvoModel m;
  m.t_max = element_as_double(dist_config, "t_max", "dist_config");
  m.t_grid = element_as_double_vector(dist_config, "t_grid", "dist_config");
  const std::vector<double> logp =
    element_as_double_vector(dist_config, "logp", "dist_config");
  m.lambda = element_as_double_vector(dist_config, "lambda", "dist_config");
  m.evec = element_as_double_vector(dist_config, "evec", "dist_config");
  m.evec_inv = element_as_double_vector(dist_config, "evec_inv", "dist_config");
  const std::string states =
    element_as_string(dist_config, "states", "dist_config");

  const int n_grid = m.n_grid();
  const std::size_t n_logp =
    static_cast<std::size_t>(n_grid) * EvoModel::kCode * EvoModel::kCode;
  if (n_grid < 3 || logp.size() != n_logp) {
    OPTIMOTU_STOP("dist_evo logp does not match t_grid");
  }
  if (m.lambda.size() != EvoModel::kState ||
      m.evec.size() != EvoModel::kState * EvoModel::kState ||
      m.evec_inv.size() != EvoModel::kState * EvoModel::kState ||
      states.size() != EvoModel::kState) {
    OPTIMOTU_STOP("dist_evo spectral factors must cover 20 states");
  }
  m.logp.resize(n_logp);
  for (std::size_t i = 0; i < n_logp; ++i) {
    m.logp[i] = static_cast<float>(logp[i]);
  }
  for (int i = 0; i < EvoModel::kState; ++i) {
    const char ch = states[static_cast<std::size_t>(i)];
    if (ch < 'A' || ch > 'Z') {
      OPTIMOTU_STOP("dist_evo state letters must be A-Z");
    }
    const int code = ch - 'A' + 1;
    m.code_to_state[static_cast<std::size_t>(code)] = i;
    m.is_state[static_cast<std::size_t>(code)] = 1;
  }
  return m;
}

#endif

#endif
