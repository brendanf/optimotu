#ifndef OPTIMOTU_SINGLE_LINKAGE_H_INCLUDED
#define OPTIMOTU_SINGLE_LINKAGE_H_INCLUDED

#include <cstdint>
#include <limits>

// type of "j", which is the index of a cluster
typedef std::uint32_t j_t;
// type of "d", which is the index corresponding to a distance
typedef std::int32_t d_t;

static const j_t NO_CLUST = std::numeric_limits<std::uint32_t>::max();
static const d_t NO_DIST = std::numeric_limits<std::int32_t>::max();

#ifdef OPTIMOTU_R
#define OPTIMOTU_COUT Rcpp::Rcout
#define OPTIMOTU_CERR Rcpp::Rcerr
#define OPTIMOTU_STOP Rcpp::stop
#else
#include <iostream>
#define OPTIMOTU_COUT std::cout
#define OPTIMOTU_CERR std::cerr
#define OPTIMOTU_STOP(s) do { std::cerr << (s); std::exit(1); } while (false)
#endif // OPTIMOTU_R
#endif //OPTIMOTU_SINGLELINKAGE_H_INCLUDED
