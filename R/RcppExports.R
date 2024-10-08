# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Confusion matrix for a set of "test" partitions vs. a "true" partition
#'
#' One way to analyze the comparison of two different partitions of the same
#' data is to treat it as a binary classification problem operating on pairs
#' of objects, where pairs should be classified as belonging to the same
#' cluster or not.  Then a "true positive" is a pair which belong to the same
#' cluster in both the "test" partition and the "true" partition; a "false
#' positive" is a pair which belongs to the same cluster in the "test"
#' partition but not the "true" partition; a "false negative" is a pair which
#' belongs to the same cluster in the "true" partition but not in the "test"
#' partition; and a "true" negative is a pair which does not belong to the same
#' cluster in either the "test" partition or the "true" partition.
#'
#' This formulation allows various measures of binary classification
#' performance to be applied to the case of clustering.
#'
#' @return (`data.frame` with `m` rows) for each "test" partition, the number
#' of true positives ("`TP`"), false positives ("`FP`"), false negatives
#' ("`FN`"), and true negatives ("`TN`") relative to the "true" partition.
#'
#' @export
#' @inheritParams mutual_information
confusion_matrix <- function(k, c, threads = 1L) {
    .Call(`_optimotu_confusion_matrix`, k, c, threads)
}

#' @export
#' @describeIn confusion_matrix confusion matrix for different references at
#'  each threshold
confusion_matrix2 <- function(k, c, threads = 1L) {
    .Call(`_optimotu_confusion_matrix2`, k, c, threads)
}

distmx_cluster_single <- function(file, seqnames, threshold_config, clust_config, parallel_config, output_type = "matrix", verbose = FALSE) {
    .Call(`_optimotu_distmx_cluster_single`, file, seqnames, threshold_config, clust_config, parallel_config, output_type, verbose)
}

distmx_cluster_multi <- function(file, seqnames, which, threshold_config, method_config, parallel_config, output_type = "matrix", verbose = FALSE) {
    .Call(`_optimotu_distmx_cluster_multi`, file, seqnames, which, threshold_config, method_config, parallel_config, output_type, verbose)
}

#' Sparse distance matrix between DNA sequences
#'
#' @param seq (`character` vector) DNA sequences to calculate distances for
#' @param dist_threshold (`numeric` scalar) maximum sequence distance (edit
#' distance / alignment length) threshold for reporting
#' @param constrain (`logical` flag) if `TRUE`, the alignment algorithm will
#' use optimizations that will cause it to exit early if the optimal alignment
#' has a distance greater than the distance threshold. This should not change
#' the correctness of distance calculations below the threshold, and results in
#' a large speedup. It is recommended to use `constrain=FALSE` only to verify
#' that the results do not change.
#' @param threads (`integer` count) number of parallel threads to use for
#' computation.
#'
#' @return (`data.frame`) a sparse distance matrix; columns are `seq1` and
#' `seq2` for the 0-based indices of two sequences; `score1` and `score2` are
#' the optimal alignment score for the two sequences in the "prealignment" (if
#' any) and "alignment" stages; `dist1` and `dist2` are the corresponding
#' sequence distances.
#'
#' @export
#' @rdname seq_distmx
seq_distmx_edlib <- function(seq, dist_threshold, constrain = TRUE, threads = 1L) {
    .Call(`_optimotu_seq_distmx_edlib`, seq, dist_threshold, constrain, threads)
}

#' @param breakpoint (`numeric` scalar) threshold for deciding whether to use
#' WFA2 or edlib for edit-distance alignment.  This parameter is interpreted as
#' an edit distance if greater than or equal to `1`, or as a pairwise
#' dissimilarity if less than 1. In either case, WFA2 is used below the
#' breakpoint, and edlib is used above it.
#' @export
#' @rdname seq_distmx
seq_distmx_hybrid <- function(seq, dist_threshold, breakpoint = 0.1, threads = 1L) {
    .Call(`_optimotu_seq_distmx_hybrid`, seq, dist_threshold, breakpoint, threads)
}

#' Size of the intersection between two sorted sets
#'
#' This implementation is much faster that `length(intersect(c, k))`. However
#' it assumes (without checking!) that the sets are sorted.
#'
#' @param c (sorted `integer` or ` character` vector) the first set to compare
#' @param k (sorted `integer` or ` character` vector) the second set to compare
#'
#' @return (`integer` count) the number of elements which occur in both sets
#' @export
intersect_length <- function(c, k) {
    .Call(`_optimotu_intersect_length`, c, k)
}

#' @export
#' @rdname intersect_length
intersect_length_string <- function(c, k) {
    .Call(`_optimotu_intersect_length_string`, c, k)
}

inner_fmeasure <- function(cj, kpartition, nk) {
    .Call(`_optimotu_inner_fmeasure`, cj, kpartition, nk)
}

fmeasure_list <- function(k, c, ncpu = 1L) {
    .Call(`_optimotu_fmeasure_list`, k, c, ncpu)
}

fmeasure_matrix <- function(k, c, ncpu = 1L) {
    .Call(`_optimotu_fmeasure_matrix`, k, c, ncpu)
}

#' @param udist_threshold (`numeric` scalar between 0 and 1) maximum udist
#' (number of shared kmers / number of kmers in the shorter sequence) for full
#' alignment.
#' @export
#' @rdname seq_distmx
seq_distmx_kmer <- function(seq, dist_threshold, udist_threshold, match = 1L, mismatch = 2L, gap_open = 10L, gap_extend = 1L, gap_open2 = 0L, gap_extend2 = 0L, threads = 1L) {
    .Call(`_optimotu_seq_distmx_kmer`, seq, dist_threshold, udist_threshold, match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2, threads)
}

#' Calculate similarity for a set of alternate "test" partitions vs. a "true" partition
#'
#' @param k (`m` x `n` `integer` matrix) `m` alternative "test" partitions; each
#' row gives the cluster assignment for the `n` objects. Objects with the same
#' cluster ID are clustered together. Cluster IDs do not need to be
#' consecutive, and they do not need to correspond between different
#' partitions.
#' @param c (`integer` vector of length `n`) "True" partition of the `n`
#' objects.
#' @param threads (`integer` count) number of parallel threads to use.
#'
#' @return (`numeric` vector of length `m`) The similarity measure between each
#' of the alternative partitions and the "true" partition.
#' @export
mutual_information <- function(k, c, threads = 1L) {
    .Call(`_optimotu_mutual_information`, k, c, threads)
}

#' @rdname mutual_information
#' @export
adjusted_mutual_information <- function(k, c, threads = 1L) {
    .Call(`_optimotu_adjusted_mutual_information`, k, c, threads)
}

cigar_wfa2 <- function(a, b, match = 0L, mismatch = 1L, open1 = 0L, extend1 = 1L, open2 = 0L, extend2 = 1L) {
    .Call(`_optimotu_cigar_wfa2`, a, b, match, mismatch, open1, extend1, open2, extend2)
}

cigar_edlib <- function(a, b) {
    .Call(`_optimotu_cigar_edlib`, a, b)
}

#' Align two sequences using WFA2 and return the pairwise distance
#'
#' @description WFA2 allows several alignment strategies, and this function
#' selects the least parameterized version possible given the inputs:
#'
#'  - Edit : `gap` == `mismatch` != 0; `match` == `extend` == `gap2` == `extend2` == 0.
#'  - Indel : `gap` != 0; `mismatch` == `match` == `extend` == `gap2` == `extend2` == 0.
#'  - GapLinear : `extend` == `gap2` == `extend2` == 0; other parameters do not meet requirements for "Edit" or "Indel".
#'  - GapAffine : `gap2` == `extend2` == 0; other parameters do not meet requirements for "Edit", "Indel", or "GapLinear".
#'  - GapAffine2Pieces : parameters do not meet requirements for `Edit`, `Indel`, `GapLinear`, or `GapAffine`.
#'
#' @param a (`character` string) first string to align
#' @param b (`character` string) second string to align
#' @param match (`integer` scalar) match score; positive is a bonus.
#' @param mismatch (`integer` scalar) mismatch score; positive is a penalty.
#' @param gap (`integer` scalar) gap opening score; positive is a penalty.
#' Alternatively, if `extend`, `gap2`, and `extend2` are all 0, then this is
#' the penalty for all gaps.
#' @param extend (`integer` scalar) gap extension score; positive is a penalty
#' @param gap2 (`integer` scalar) alternate gap opening score for two-piece
#' affine gap penalty; positive is penalty.  Ignored if both `gap2` and
#' `extend2` are 0.
#' @param extend2 (`integer` scalar) alternate gap extension score for
#' two-piece affine gap penalty; positive is penalty. Ignored if both `gap2`
#' and `extend2` are 0.
#' @export
#'
align <- function(a, b, match = 0L, mismatch = 1L, gap = 1L, extend = 0L, gap2 = 0L, extend2 = 0L) {
    .Call(`_optimotu_align`, a, b, match, mismatch, gap, extend, gap2, extend2)
}

#' @param prealign (`logical` flag) if `TRUE`, do a prealignment using
#' edit-distance as alignment score (`match = 0, mismatch = 1, gap_open = 0,
#' gap_extend = 1, gap_open2 = 0, gap_extend2 = 1`) to test feasibility before
#' aligning with alternate alignment scores. Note that, since pairwise distance
#' is defined using edit distance, any other set of scores will always result
#' in a pairwise distance which is equal to or greater than an alignment based
#' on the edit distance score.
#' @export
#' @rdname seq_distmx
seq_distmx_wfa2 <- function(seq, dist_threshold, match = 1L, mismatch = 2L, gap_open = 10L, gap_extend = 1L, gap_open2 = 0L, gap_extend2 = 0L, prealign = TRUE, constrain = TRUE, threads = 1L) {
    .Call(`_optimotu_seq_distmx_wfa2`, seq, dist_threshold, match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2, prealign, constrain, threads)
}

seq_cluster_single <- function(seq, dist_config, threshold_config, clust_config, parallel_config, output_type = "matrix", verbose = FALSE) {
    .Call(`_optimotu_seq_cluster_single`, seq, dist_config, threshold_config, clust_config, parallel_config, output_type, verbose)
}

seq_cluster_multi <- function(seq, which, dist_config, threshold_config, clust_config, parallel_config, output_type = "matrix", verbose = FALSE) {
    .Call(`_optimotu_seq_cluster_multi`, seq, which, dist_config, threshold_config, clust_config, parallel_config, output_type, verbose)
}

#' @param match (non-negative `integer`) alignment score for matching nucleotides
#' @param mismatch (non-negative `integer`) alignment penalty for mismatched
#' nucleotides.
#' @param gap_open (non-negative `integer`) alignment penalty for opening a new
#' gap (i.e., insertion or deletion).
#' @param gap_extend (non-negative `integer`) alignment penalty for each
#' position in a gap.
#' @param gap_open2 (non-negative `integer`) alternate alignment penalty for
#' opening a new gap (i.e., insertion or deletion).
#' @param gap_extend2 (non-negative `integer`) alternate alignment penalty for
#' each position in a gap.
#' @rdname seq_distmx
#' @export
seq_distmx_snsn <- function(seq, dist_threshold, match = 1L, mismatch = 2L, gap_open = 10L, gap_extend = 1L, gap_open2 = 0L, gap_extend2 = 0L, constrain = TRUE, threads = 1L) {
    .Call(`_optimotu_seq_distmx_snsn`, seq, dist_threshold, match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2, constrain, threads)
}

