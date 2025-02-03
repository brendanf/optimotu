#ifndef OPTIMOTU_H_INCLUDED

#ifdef OPTIMOTU_R

#include <RcppThread.h>
#define OPTIMOTU_COUT RcppThread::Rcout
#define OPTIMOTU_CERR RcppThread::Rcerr
#define OPTIMOTU_STOP Rcpp::stop
#define OPTIMOTU_DEBUG(level, s) do {\
  if constexpr (level <= verbose) RcppThread::Rcerr s;\
} while (false)

#else

#include <iostream>
#define OPTIMOTU_COUT std::cout
#define OPTIMOTU_CERR std::cerr
#define OPTIMOTU_STOP(s) do { std::cerr << (s); std::exit(1); } while (false)
#define OPTIMOTU_DEBUG(level, s) do {\
  if constexpr (level <= verbose) std::cerr s;\
} while (false)

#endif // OPTIMOTU_R

#endif // OPTIMOTU_H_INCLUDED
