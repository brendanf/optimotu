#ifndef OPTIMOTU_H_INCLUDED

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

#endif // OPTIMOTU_H_INCLUDED
