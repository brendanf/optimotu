# SPDX-CopyrightText: (c) 2025 Brendan Furneaux
# SPDX-License-Identifier: MIT

#' Convert taxonomic ranks from character vectors and integers to ordered factors
#'
#' The least inclusive rank is given the smallest value in the ordering,
#' so for instance `"species" < "kingdom"`.
#'
#' @param x the taxonomic ranks to convert
#' @return an ordered factor of taxonomic ranks
#' @keywords internal
rank2factor <- function(x, ranks) {
  factor(x, levels = rev(ranks), ordered = TRUE)
}

#' @rdname rank2factor
#' @keywords internal
int2rankfactor <- function(x, ranks) {
  rank2factor(ranks[x])
}

#' Get superordinate or subordinate taxonomic ranks
#'
#' @param x (`character` string) the taxonomic rank to find superordinate or
#' subordinate ranks for
#' @param ranks (`character` vector) the full list of taxonomic ranks, in order
#' from most inclusive to least inclusive
#' @return a `character` vector of superordinate or subordinate taxonomic ranks
#' @keywords internal
superranks <- function(x, ranks) {
  ranks[rank2factor(ranks, ranks) > x]
}

#' @rdname superranks
#' @keywords internal
subranks <- function(x, ranks) {
  ranks[rank2factor(ranks, ranks) < x]
}
