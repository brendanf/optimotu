# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @export
intersect_length <- function(c, k) {
    .Call('_optimotu_intersect_length', PACKAGE = 'optimotu', c, k)
}

#' @export
intersect_length_string <- function(c, k) {
    .Call('_optimotu_intersect_length_string', PACKAGE = 'optimotu', c, k)
}

#' @export
inner_fmeasure <- function(cj, kpartition, nk) {
    .Call('_optimotu_inner_fmeasure', PACKAGE = 'optimotu', cj, kpartition, nk)
}

#' @export
fmeasure <- function(k, c, ncpu = 1L) {
    .Call('_optimotu_fmeasure', PACKAGE = 'optimotu', k, c, ncpu)
}

#' @export
single_linkage_matrix_thread <- function(file, seqnames, dmin, dmax, dstep, threads = 1L, minsplit = 1L) {
    .Call('_optimotu_single_linkage_matrix_thread', PACKAGE = 'optimotu', file, seqnames, dmin, dmax, dstep, threads, minsplit)
}

#' @export
single_linkage_pool <- function(file, seqnames, dmin, dmax, dstep) {
    .Call('_optimotu_single_linkage_pool', PACKAGE = 'optimotu', file, seqnames, dmin, dmax, dstep)
}

#' @export
single_linkage_multi <- function(file, seqnames, dmin, dmax, dstep, preclust, threads = 1L) {
    .Call('_optimotu_single_linkage_multi', PACKAGE = 'optimotu', file, seqnames, dmin, dmax, dstep, preclust, threads)
}

