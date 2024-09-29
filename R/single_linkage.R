verify_which <- function(which, seqnames) {
   checkmate::assert_list(
      which,
      any.missing = FALSE
   )
   if (!all(vapply(lapply(which, `%in%`, table = seqnames), all, TRUE))) {
      stop("Elements of 'which' must all be in 'seqnames'.")
   }
   invisible(TRUE)
}



is_list_of_character <- function(x) {
  all(vapply(x, is.character, TRUE))
}

#' Single linkage clustering at multiple thresholds from a sparse distance
#' matrix
#'
#' @description This function is designed to reduce CPU and memory requirements
#' of clustering by doing multiple clustering "jobs" concurrently from the same
#' sparse distance matrix, which is read once and never kept in memory.  This is
#' especially useful when the algorithm generating the distance matrix can also
#' operate in a streaming fashion, so that distance matrices which are too large
#' to fit in memory, or even to disk, can still be utilized. The motivating
#' application is clustering biological sequences, but the implementation is
#' totally agnostic about the nature of the entities it is clustering; it only
#' needs distances.
#'
#' @param distmx (`character` filename) The name of the file (or, e.g., a named
#' pipe) from which to read the sparse distance matrix. The matrix should be
#' a white-space delimited text file, with three values per line: id1 id2 dist,
#' where id1 and id2 are integer indices of the objects to be clustered,
#' starting with 0, and dist is the distance, typically a real number between 0
#' and 1 (but the algorithm works for any non-negative distance.)
#' @param names (`character` vector) Names of the sequences, used for
#' labeling the columns of the output matrix. Typically, in order to generate
#' the sparse distance matrix with integer indices, an alternate version of the
#' input may be used, where the "real" names are replaced with integers. These
#' are the "real" names.
#' @param threshold_config (`optimotu_threshold_config` object returned by
#' [threshold_config()] or one of its helper functions) Definition of the
#' thresholds to use for clustering.
#' @param clust_config (`optimotu_cluster_config` object returned by
#' [clust_config()] or one of its helpers) The clustering algorithm to use; all
#' algorithms give identical results, but may have different performance
#' characteristics on different problems.
#' @param parallel_config (`optimotu_parallel_config` object returned by
#' [parallel_config()] or one of its helpers) The method to use for
#' parallel clustering. For single-threaded clustering, use the default value
#' ([parallel_concurrent()] with `threads = 1`).
#' @param output_type (`character`) Which type of output to give; one of
#' `"matrix"` or `"hclust"`. `"matrix"` returns an integer matrix giving
#' clustering results, where the element in row `i` and column `j` gives the
#' 0-based index of the first member of the cluster to which sequence `j`
#' belongs when clustered at the `i`th clustering threshold. `"hclust"` returns
#' an object as returned by [stats::hclust()], which requires less memory,
#' especially for large problems. If `which` is given, then both `output_type`s
#' instead return a list whose elements are of the chosen type.
#' @param which (`list` of `character` vectors) Instead of performing clustering
#' on all input sequences, perform independent clustering on subsets of the
#' sequences defined by the elements of `which`. Subsets do not need to be
#' disjoint (and indeed, if they are it is probably faster to calculate the
#' distance matrices separately.)
#'
#' @return An [`integer matrix`][methods::structure-class] if
#' `output_type=="matrix"`, an [`hclust`][stats::hclust] object if
#' `output_type=="hclust"`, or a list of one of these when `which` is a list.
#' @export
distmx_cluster = function(
   distmx,
   names,
   threshold_config,
   clust_config = clust_index(),
   parallel_config = parallel_concurrent(1),
   output_type = c("matrix", "hclust"),
   which = NULL
) {
  output_type = match.arg(output_type)
  out <- if (!is.null(which) && !isTRUE(which)) {
    verify_which(which, names)
    distmx_cluster_multi(
      distmx,
      names,
      which,
      threshold_config,
      clust_config,
      parallel_config,
      output_type
    )
  } else {
    distmx_cluster_single(
      distmx,
      names,
      threshold_config,
      clust_config,
      parallel_config,
      output_type
    )
  }
  out <- reduplicate_thresholds(out, threshold_config)
  out <- rename_thresholds(out, threshold_config)
  out
}
