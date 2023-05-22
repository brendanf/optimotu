#ifndef OPTIMOTU_CONFIG_H_INCLUDED
#define OPTIMOTU_CONFIG_H_INCLUDED

#ifdef OPTIMOTU_R

#include <Rcpp.h>

#include "ClusterTree.h"
#include "ClusterMatrix.h"
#include "ClusterIndexedMatrix.h"

#include "DistanceConverter.h"

std::unique_ptr<DistanceConverter> create_distance_converter(Rcpp::List config);
std::unique_ptr<ClusterAlgorithm> create_cluster_algorithm(
    Rcpp::List config,
    DistanceConverter * dconv,
    Rcpp::IntegerMatrix &im
);

#endif //OPTIMOTU_R

#endif //OPTIMOTU_CONFIG_H_INCLUDED
