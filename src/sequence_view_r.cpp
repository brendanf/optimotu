// SPDX-FileCopyrightText: 2025 Brendan Furneaux
// SPDX-License-Identifier: MIT

#ifdef OPTIMOTU_R

#include "SequenceView.h"
#include "optimotu.h"
#include <Rcpp.h>

// Extract non-owning views of R character strings on the main thread.
// The CharacterVector must remain alive for the duration of the views.
SequenceSet sequence_views_from_r(const Rcpp::CharacterVector &seq) {
  const R_xlen_t n = seq.length();
  SequenceSet out;
  out.reserve(static_cast<std::size_t>(n));
  for (R_xlen_t i = 0; i < n; ++i) {
    SEXP elt = STRING_ELT(seq, i);
    if (elt == NA_STRING) {
      OPTIMOTU_STOP("sequence vector cannot contain NA");
    }
    out.emplace_back(CHAR(elt), static_cast<std::size_t>(LENGTH(elt)));
  }
  return out;
}

#endif // OPTIMOTU_R
