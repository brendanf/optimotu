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
#' @param span (`character` scalar) alignment span; one of "global" or "extend".
#' In "global" mode, the entire sequences are aligned and end gaps are penalized.
#' "extend" mode is slightly between methods: in both cases, left end gaps are
#' penalized. For WFA2, right end gaps in both sequences are not penalized.
#' For Edlib, right end gaps in the *second* sequence are not penalized, but
#' end gaps in the *first* sequence are penalized.
#' @export
#' @rdname pairwise_alignment
align <- function(a, b, match = 0, mismatch = 1, gap_open = 0, gap_extend = 1,
                  gap_open2 = 0, gap_extend2 = 0, method = c("wfa2", "edlib"),
                  span = c("global", "extend")) {
  # Check that the inputs are valid
  checkmate::assert_string(a)
  checkmate::assert_string(b)
  checkmate::assert_integerish(match)
  checkmate::assert_integerish(mismatch)
  checkmate::assert_integerish(gap_open)
  checkmate::assert_integerish(gap_extend)
  checkmate::assert_integerish(gap_open2)
  checkmate::assert_integerish(gap_extend2)
  checkmate::assert_character(span)
  span <- match.arg(span)
  checkmate::assert_choice(span, c("global", "extend"))

  checkmate::assert_character(method)
  method <- match.arg(method)
  if (method == "edlib") {
    if (match != 0 || mismatch != 1 || gap_open != 0 || gap_extend != 1 ||
        gap_open2 != 0 || gap_extend2 != 0) {
      warning("Edlib only supports match = 0, mismatch = 1, open = 0, extend = 1, open2 = 0, extend2 = 0")
    }
    switch(
      span,
      global = align_edlib_global(a, b),
      extend = align_edlib_extend(a, b)
    )
  } else {
    switch(
      span,
      global = align_wfa2_global(a, b, -match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2),
      extend = align_wfa2_extend(a, b, -match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2)
    )
  }
}

#' @return (`character(1)`) CIGAR string
#' @export
#' @keywords internal
#' @describeIn pairwise_alignment Generate alignment CIGAR with WFA2
cigar_wfa2 <- function(a, b, match = 0, mismatch = 1, gap_open = 0,
                                  gap_extend = 1, gap_open2 = 0, gap_extend2 = 0,
                       span = c("global", "extend")) {
  # Check that the inputs are valid
  checkmate::assert_string(a)
  checkmate::assert_string(b)
  checkmate::assert_integerish(match)
  checkmate::assert_integerish(mismatch)
  checkmate::assert_integerish(gap_open)
  checkmate::assert_integerish(gap_extend)
  checkmate::assert_integerish(gap_open2)
  checkmate::assert_integerish(gap_extend2)
  checkmate::assert_character(span)
  span <- match.arg(span)
  checkmate::assert_choice(span, c("global", "extend"))

  switch(
    span,
    global = cigar_wfa2_global(a, b, -match, mismatch, gap_open,
                               gap_extend, gap_open2, gap_extend2),
    extend = cigar_wfa2_extend(a, b, -match, mismatch, gap_open,
                               gap_extend, gap_open2, gap_extend2)
  )
}

cigar_edlib <- function(a, b, span = c("global", "extend")) {
  # Check that the inputs are valid
  checkmate::assert_string(a)
  checkmate::assert_string(b)
  checkmate::assert_character(span)
  span <- match.arg(span)
  checkmate::assert_choice(span, c("global", "extend"))

  # Call edlib to get the CIGAR string
  switch(
    span,
    global = cigar_edlib_global(a, b),
    extend = cigar_edlib_extend(a, b)
  )
}
