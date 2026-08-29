// SPDX-FileCopyrightText: 2025 Brendan Furneaux
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_CLUSTER_RUN_H_INCLUDED
#define OPTIMOTU_CLUSTER_RUN_H_INCLUDED

#ifdef OPTIMOTU_R

#include "config.h"
#include "MemoryBudget.h"
#include <istream>
#include <memory>
#include <string>
#include <vector>

void stop_memory_budget(const MemoryBudgetExceeded &e);

struct MultiClusterJob {
  std::unique_ptr<DistanceConverter> dconv;
  std::unique_ptr<ClusterAlgorithmFactory> factory;
  std::unique_ptr<MultipleClusterAlgorithm> algo;
};

MultiClusterJob run_seq_cluster_multi(
  const std::vector<std::string> &cppseq,
  Rcpp::CharacterVector seqnames,
  Rcpp::ListOf<Rcpp::CharacterVector> which,
  Rcpp::List dist_config,
  Rcpp::List threshold_config,
  Rcpp::List clust_config,
  Rcpp::List parallel_config,
  int verbose,
  double clustering_memory_budget_mb
);

MultiClusterJob run_distmx_cluster_multi(
  std::istream &infile,
  Rcpp::CharacterVector seqnames,
  Rcpp::ListOf<Rcpp::CharacterVector> which,
  Rcpp::List threshold_config,
  Rcpp::List clust_config,
  Rcpp::List parallel_config,
  bool by_name,
  bool verbose
);

#endif // OPTIMOTU_R

#endif // OPTIMOTU_CLUSTER_RUN_H_INCLUDED
