// SPDX-FileCopyrightText: 2025 Brendan Furneaux <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_CONFIG_H_INCLUDED
#define OPTIMOTU_CONFIG_H_INCLUDED

#include "ClusterAlgorithmFactory.h"
#include "MultipleClusterAlgorithm.h"
#include "ClusterWorker.h"
#include "SearchWorker.h"

#include "DistanceConverter.h"

#include "AlignClusterWorker.h"

#ifdef OPTIMOTU_R
#include <Rcpp.h>
#endif //OPTIMOTU_R

std::unique_ptr<DistanceConverter> create_distance_converter(
    const std::string &type,
    const double from,
    const double to,
    const double by,
    std::vector<double> thresholds,
    double precision
);
std::unique_ptr<ClusterAlgorithmFactory> create_cluster_algorithm(
    const std::string &method,
    const bool do_binary_search,
    const int fill_type,
    DistanceConverter * dconv
);
std::unique_ptr<MultipleClusterAlgorithm> create_multiple_cluster_algorithm(
    const int threads,
    ClusterAlgorithmFactory &factory,
    const std::vector<std::string> seqnames,
    const std::vector<std::vector<std::string>> subset_names
);
std::unique_ptr<ClusterWorker> create_cluster_worker(
    const std::string &method,
    const int threads,
    const int shards,
    ClusterAlgorithm * algo,
    std::istream &file,
    const bool by_name,
    const std::vector<std::string> & seqnames
);
std::unique_ptr<AlignClusterWorker> create_align_cluster_worker(
    const std::string &type,
    const std::vector<std::string> &seq,
    const double breakpoint,
    ClusterAlgorithm &cluster,
    const std::uint8_t threads,
    int verbose = 0
);

#ifdef OPTIMOTU_R

double element_as_double(Rcpp::List list, const std::string &e, const std::string &name);
int element_as_int(Rcpp::List list, const std::string &e, const std::string &name);
std::string element_as_string(Rcpp::List list, const std::string &e, const std::string &name);
std::vector<double> element_as_double_vector(Rcpp::List list, const std::string &e, const std::string &name);
bool element_as_bool(Rcpp::List list, const std::string &e, const std::string &name);

std::unique_ptr<DistanceConverter> create_distance_converter(Rcpp::List config);
std::unique_ptr<ClusterAlgorithmFactory> create_cluster_algorithm(
    Rcpp::List config,
    DistanceConverter * dconv
);
std::unique_ptr<MultipleClusterAlgorithm> create_multiple_cluster_algorithm(
    Rcpp::List parallel_config,
    ClusterAlgorithmFactory &factory,
    Rcpp::CharacterVector seqnames,
    Rcpp::ListOf<Rcpp::CharacterVector> subset_names
);
std::unique_ptr<ClusterWorker> create_cluster_worker(
  Rcpp::List config,
  ClusterAlgorithm * algo,
  std::istream &file,
  bool by_name,
  const Rcpp::CharacterVector seqnames
);
std::unique_ptr<AlignClusterWorker> create_align_cluster_worker(
    Rcpp::List dist_config,
    Rcpp::List parallel_config,
    const std::vector<std::string> &seq,
    ClusterAlgorithm &cluster,
    int verbose = 0
);

std::unique_ptr<SearchWorker> create_search_worker(
    Rcpp::List dist_config,
    Rcpp::List parallel_config,
    const std::vector<std::string> & query,
    const std::vector<std::string> & ref,
    double threshold,
    int verbose = 0
);

#endif //OPTIMOTU_R

#endif //OPTIMOTU_CONFIG_H_INCLUDED
