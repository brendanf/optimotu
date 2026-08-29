// SPDX-FileCopyrightText: 2025 Brendan Furneaux
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_CLUSTER_MEASURES_H_INCLUDED
#define OPTIMOTU_CLUSTER_MEASURES_H_INCLUDED

#ifdef OPTIMOTU_R

#include "MultipleClusterAlgorithm.h"
#include <Rcpp.h>
#include <string>
#include <vector>

// Score one subset from a live MultipleClusterAlgorithm without allocating
// an n x m cluster matrix. Returns measure / threshold / value as in
// find_best_threshold().
Rcpp::DataFrame find_best_threshold_from_algo(
  MultipleClusterAlgorithm &algo,
  std::size_t subset,
  const Rcpp::IntegerVector &c,
  const std::vector<std::string> &measures,
  const Rcpp::NumericVector &thresholds,
  const Rcpp::IntegerVector &threshold_order,
  int threads
);

// Fill an m x n IntegerMatrix for one subset via write_threshold_row.
Rcpp::IntegerMatrix subset_matrix_via_rows(
  MultipleClusterAlgorithm &algo,
  std::size_t subset
);

#endif // OPTIMOTU_R

#endif // OPTIMOTU_CLUSTER_MEASURES_H_INCLUDED
