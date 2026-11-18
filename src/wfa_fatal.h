#ifndef OPTIMOTU_WFA_FATAL_H
#define OPTIMOTU_WFA_FATAL_H

#include <atomic>

// Set when a WFA alignment runs out of memory or longjmps out of exit().
// Workers only store the flag. The main thread calls
// optimotu_wfa_stop_if_fatal() after the parallel join.
inline std::atomic<int> optimotu_wfa_fatal_flag{0};

inline void optimotu_wfa_note_fatal() {
  optimotu_wfa_fatal_flag.store(1, std::memory_order_relaxed);
}

inline bool optimotu_wfa_fatal_pending() {
  return optimotu_wfa_fatal_flag.load(std::memory_order_relaxed) != 0;
}

inline void optimotu_wfa_fatal_clear() {
  optimotu_wfa_fatal_flag.store(0, std::memory_order_relaxed);
}

#ifdef OPTIMOTU_R
#include <Rcpp.h>

inline void optimotu_wfa_stop_if_fatal() {
  if (optimotu_wfa_fatal_flag.exchange(0, std::memory_order_relaxed) != 0) {
    Rcpp::stop(
        "WFA alignment failed (out of memory or fatal aligner error)");
  }
}
#endif

#endif
