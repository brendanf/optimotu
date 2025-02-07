#' Single-linkage clustering of nucleotide sequences
#'
#' @inheritParams seq_cluster_usearch
#' @param dist_config (`optimotu_dist_config` object returned by
#' [dist_config()] or one of its helpers) Configuration of the method to
#' calculate distances. If `dist_usearch()`, then this function dispatches to
#' `seq_cluster_usearch()`.
#' @param verbose (`logical(1)` or `integer(1)`) whether to print progress;
#' values greater than 1 (or TRUE) print more
#' @export
seq_cluster <- function(
    seq,
    seq_id = names(seq),
    dist_config,
    threshold_config,
    clust_config = clust_tree(),
    parallel_config = parallel_concurrent(1),
    output_type = c("matrix", "hclust"),
    which = TRUE,
    verbose = FALSE
) {
  UseMethod("seq_cluster", seq)
}

#' @method seq_cluster data.frame
#' @export
seq_cluster.data.frame <- function(
    seq,
    seq_id = seq$seq_id,
    dist_config,
    threshold_config,
    clust_config = clust_tree(),
    parallel_config = parallel_concurrent(1),
    output_type = c("matrix", "hclust"),
    which = TRUE,
    verbose = FALSE
) {
  mycall <- match.call()
  if (missing(seq_id)) {
    newseq_id <- quote(seq$seq_id)
    newseq_id[[2]] <- mycall$seq
    mycall$seq_id <- newseq_id
  }
  if (identical(dist_config$method, "usearch")) {
    mycall[[1]] <- quote(seq_cluster_usearch.DNAStringSet)
    newseq <- quote(Biostrings::DNAStringSet(seq$seq))
    newseq[[2]][[2]] <- mycall$seq
    mycall$usearch <- dist_config$usearch
    mycall$usearch_ncpu <- dist_config$usearch_ncpu
    mycall$dist_config <- NULL
  } else {
    mycall[[1]] <- quote(seq_cluster.character)
    newseq <- quote(seq$seq)
    newseq[[2]] <- mycall$seq
  }
  mycall$seq <- newseq
  eval(mycall, envir = parent.frame())
}

#' @export
seq_cluster.character <- function(
    seq,
    seq_id = names(seq),
    dist_config,
    threshold_config,
    clust_config = clust_tree(),
    parallel_config = parallel_concurrent(1),
    output_type = c("matrix", "hclust"),
    which = TRUE,
    verbose = FALSE
) {
  checkmate::assert_class(dist_config, "optimotu_dist_config")
  if (identical(dist_config$method, "usearch")) {
    mycall <- match.call
    mycall$dist_config <- NULL
    mycall$usearch <- dist_config$usearch
    mycall$usearch_ncpu <- dist_config$usearch_ncpu
    mycall[[1]] <- quote(seq_cluster_usearch.character)
    return(eval(mycall, envir = parent.frame()))
  }
  if (identical(dist_config$method, "file")) {
    mycall <- match.call
    mycall$dist_config <- NULL
    mycall$seq <- NULL
    mycall$names <- mycall$seq_id
    mycall$seq_id <- NULL
    mycall$by_names <- dist_config$by_names
    mycall[[1]] <- quote(distmx_cluster)
    return(eval(mycall, envir = parent.frame()))
  }
  output_type = match.arg(output_type)
  if (length(seq) == 1 && file.exists(seq)) {
    seq <- as.character(Biostrings::readBStringSet(seq))
  }
  if (!missing(seq_id)) names(seq) <- seq_id
  checkmate::assert_class(threshold_config, "optimotu_threshold_config")
  checkmate::assert_class(clust_config, "optimotu_cluster_config")
  checkmate::assert_class(parallel_config, "optimotu_parallel_config")
  checkmate::assert_character(seq_id)
  out <- if (!is.list(which)) {
    seq <- seq[which]
    seq_cluster_single(
      seq,
      dist_config,
      threshold_config,
      clust_config,
      parallel_config,
      output_type,
      as.integer(verbose)
    )
  } else {
    checkmate::assert(
      checkmate::check_list(which, types = "logical", any.missing = FALSE, min.len = 1),
      checkmate::check_list(which, types = "integerish", any.missing = FALSE, min.len = 1),
      checkmate::check_list(which, types = "character", any.missing = FALSE, min.len = 1)
    )
    if (is.list(which) && !is.character(which[[1]])) {
      which <- lapply(which, `[`, x = seq_id)
    }
    verify_which(which, seq_id)
    seq_cluster_multi(
      seq = seq,
      which = which,
      dist_config = dist_config,
      threshold_config = threshold_config,
      clust_config = clust_config,
      parallel_config = parallel_config,
      output_type = output_type,
      verbose = as.integer(verbose)
    )
  }
  out <- reduplicate_thresholds(out, threshold_config)
  out <- rename_thresholds(out, threshold_config)
  out
}

#' @export
seq_cluster.DNAStringSet <- function(
    seq,
    seq_id = names(seq),
    dist_config,
    threshold_config,
    clust_config = clust_index(),
    parallel_config = parallel_concurrent(1),
    output_type = c("matrix", "hclust"),
    which = TRUE,
    verbose = FALSE
) {
  mycall <- match.call
  if (identical(dist_config$method, "usearch")) {
    mycall$dist_config <- NULL
    mycall$usearch <- dist_config$usearch
    mycall$usearch_ncpu <- dist_config$usearch_ncpu
    mycall[[1]] <- quote(seq_cluster_usearch.DNAStringSet)
    return(eval(mycall, envir = parent.frame()))
  }
  mycall[[1]] <- quote(seq_cluster.character)
  newseq <- quote(as.character(seq))
  newseq[[2]] <- mycall$seq
  mycall$seq <- newseq
  eval(mycall, envir = parent.frame())
}
