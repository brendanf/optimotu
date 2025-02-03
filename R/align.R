# SPDX-CopyrightText: (c) 2025, Brendan Furneaux
# SPDX-License-Identifier: MIT

#' Align two sequences and return the pairwise distance or CIGAR string
#'
#' @description While Edlib only calculates the edit distance, WFA2 allows
#' several alignment strategies, and this function selects the least
#' parameterized version possible given the inputs:
#'
#'  - Edit : `gap_extend` == `mismatch` != 0; `match` == `gap_open` ==
#'    `gap_open2` == `gap_extend2` == 0.
#'  - Indel : `gap_extend` != 0; `mismatch` == `match` == `gap_open` ==
#'   `gap_open2` == `gap_extend2` == 0.
#'  - GapLinear : `gap_open` == `gap_open2` == `gap_extend2` == 0; other
#'    parameters do not meet requirements for "Edit" or "Indel".
#'  - GapAffine : `gap_open2` == `gap_extend2` == 0; other parameters do not
#'    meet requirements for "Edit", "Indel", or "GapLinear".
#'  - GapAffine2Pieces : parameters do not meet requirements for `Edit`,
#'    `Indel`, `GapLinear`, or `GapAffine`.
#'
#' @param a (`character` string) first string to align
#' @param b (`character` string) second string to align
#' @param match (`integer` scalar) match score; positive is a bonus.
#' @param mismatch (`integer` scalar) mismatch score; positive is a penalty.
#' @param gap_open (`integer` scalar) per-gap opening score; positive is a penalty.
#' This penalty is applied once per run of consecutive gap characters
#' @param gap_extend (`integer` scalar) gap extension score; positive is a penalty.
#' This penalty is applied for each gap character. This is the appropriate
#' parameter to use for a linear gap penalty.
#' @param gap_open2 (`integer` scalar) alternate gap opening score for two-piece
#' affine gap penalty; positive is penalty.
#' @param gap_extend2 (`integer` scalar) alternate gap extension score for
#' two-piece affine gap penalty; positive is penalty.
#' @param method (`character` scalar) alignment method to use; one of "wfa2"
#' or "edlib".
#' @export
#' @rdname pairwise_alignment
align <- function(a, b, match = 1, mismatch = 1, gap_open = 0, gap_extend = 1,
                  gap_open2 = 0, gap_extend2 = 0, method = c("wfa2", "edlib")) {
  # Check that the inputs are valid
  checkmate::assert_string(a)
  checkmate::assert_string(b)
  checkmate::assert_integerish(match)
  checkmate::assert_integerish(mismatch)
  checkmate::assert_integerish(open)
  checkmate::assert_integerish(extend)
  checkmate::assert_integerish(open2)
  checkmate::assert_integerish(extend2)

  checkmate::assert_string(method)
  method <- match.arg(method)
  if (method == "edlib") {
    if (match != 1 || mismatch != 1 || open != 0 || extend != 1 ||
        open2 != 0 || extend2 != 0) {
      warning("Edlib only supports match = 1, mismatch = 1, open = 0, extend = 1, open2 = 0, extend2 = 0")
    }
    align_edlib(a, b)
  } else {
    align_wfa2(a, b, match, mismatch, open, extend, open2, extend2)
  }
}
