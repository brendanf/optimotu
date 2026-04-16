#ifndef OPTIMOTU_H_INCLUDED

#ifdef OPTIMOTU_R
#include <Rcpp.h>
#include <RcppThread.h>
#define OPTIMOTU_COUT RcppThread::Rcout
#define OPTIMOTU_CERR RcppThread::Rcerr
#define OPTIMOTU_STOP Rcpp::stop
// Verbosity levels: 0 = off; 1--2 = high-level (once or few per top-level call);
// 3--4 = loop-level debugging (per pair/sequence). Use OPTIMOTU_DEBUG(level, ...)
// and OPTIMOTU_VERBOSE(level, ...) with level in 1--4 accordingly. Values > 4
// are treated as 4 at dispatch.
// "debug" is intended for very low-level debugging messages,
// which are implemented in templates so that there is no runtime overhead when
// not in debug mode.
// "verbose" is intended for more general messages, which will not be called
// inside critical loops.
#define OPTIMOTU_DEBUG(level, s) do {\
  if constexpr (level <= verbose) RcppThread::Rcerr s;\
} while (false)
#define OPTIMOTU_VERBOSE(level, s) do {\
  if (level <= verbose) RcppThread::Rcout s;\
} while (false)

#else

#include <iostream>
#define OPTIMOTU_COUT std::cout
#define OPTIMOTU_CERR std::cerr
#define OPTIMOTU_STOP(s) do { std::cerr << (s); std::exit(1); } while (false)
// Verbosity levels: 0 = off; 1--2 = high-level; 3--4 = loop-level. level in 1--4.
#define OPTIMOTU_DEBUG(level, s) do {\
  if constexpr (level <= verbose) std::cerr s;\
} while (false)

#endif // OPTIMOTU_R

#endif // OPTIMOTU_H_INCLUDED
