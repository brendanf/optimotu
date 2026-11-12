# SPDX-FileCopyrightText: 2026 Brendan Furneaux
# SPDX-License-Identifier: MIT

#' Translate aligned nucleotide sequences
#'
#' Lowercase A2M insertions are stripped. Remaining uppercase letters and
#' gaps must share one length. Codons of `---` become `-`; mixed
#' gap/ambiguous codons become `X`. A trailing remainder of length not
#' divisible by 3 is dropped.
#'
#' @param seq (`character`) aligned nucleotide sequences
#' @param genetic_code NCBI translation table number
#' @param codon_offset (`integer`) 0, 1, or 2
#' @return named `character` amino-acid sequences
#' @export
translate_aligned <- function(seq, genetic_code, codon_offset = 0L) {
  checkmate::assert_character(seq, any.missing = FALSE)
  checkmate::assert_integerish(genetic_code, len = 1L, lower = 1)
  checkmate::assert_integerish(codon_offset, len = 1L, lower = 0, upper = 2)
  codon_offset <- as.integer(codon_offset)
  gc <- Biostrings::getGeneticCode(as.character(as.integer(genetic_code)))
  stripped <- gsub("[a-z]", "", seq)
  widths <- nchar(stripped)
  if (length(unique(widths)) > 1L) {
    stop(
      "Aligned sequences must have equal length after removing ",
      "lowercase insertions"
    )
  }
  names_out <- names(seq)
  if (length(stripped) == 0L) {
    return(stripped)
  }
  width <- widths[[1L]]
  if (width <= codon_offset) {
    out <- rep("", length(stripped))
    names(out) <- names_out
    return(out)
  }
  usable <- width - codon_offset
  n_codon <- usable %/% 3L
  if (n_codon < 1L) {
    out <- rep("", length(stripped))
    names(out) <- names_out
    return(out)
  }
  start <- codon_offset + 1L
  aa <- vapply(
    stripped,
    function(s) {
      pieces <- character(n_codon)
      for (k in seq_len(n_codon)) {
        codon <- substr(s, start + (k - 1L) * 3L, start + k * 3L - 1L)
        codon_u <- toupper(codon)
        if (identical(codon_u, "---")) {
          pieces[[k]] <- "-"
        } else if (
          grepl("-", codon_u, fixed = TRUE) ||
            grepl("[^ACGT]", codon_u)
        ) {
          pieces[[k]] <- "X"
        } else {
          aa_k <- unname(gc[codon_u])
          pieces[[k]] <- if (length(aa_k) != 1L || is.na(aa_k)) {
            "X"
          } else {
            aa_k
          }
        }
      }
      paste(pieces, collapse = "")
    },
    character(1),
    USE.NAMES = FALSE
  )
  names(aa) <- names_out
  aa
}

maybe_translate_aa_seq <- function(seq, dist_config) {
  if (!identical(dist_config$method, "aa")) {
    return(seq)
  }
  code <- dist_config$genetic_code
  if (is.null(code)) {
    return(seq)
  }
  translate_aligned(
    seq,
    genetic_code = code,
    codon_offset = if (is.null(dist_config$codon_offset)) {
      0L
    } else {
      dist_config$codon_offset
    }
  )
}
