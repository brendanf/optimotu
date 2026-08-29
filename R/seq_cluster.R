#' Single-linkage clustering of nucleotide sequences
#'
#' @inheritParams seq_cluster_usearch
#' @param seq_idx optional 1-based subset indices; see [seq_as_char()].
#' @param seq_file optional path overrides for index-backed `seq` only; see
#'   [seq_as_char()]. Not used with `dist_file()`.
#' @param dist_config (`optimotu_dist_config` object returned by
#' [dist_config()] or one of its helpers) Configuration of the method to
#' calculate distances. If `dist_usearch()`, then this function dispatches to
#' `seq_cluster_usearch()`.
#' @param verbose (`logical(1)` or `integer(1)`) whether to print progress;
#' values greater than 1 (or TRUE) print more
#' @param clustering_memory_budget_mb (`numeric(1)` or `NULL`) Optional
#' best-effort memory budget in MB for native clustering structures in
#' multi-subset mode.
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
  verbose = FALSE,
  seq_idx = NULL,
  seq_file = NULL,
  clustering_memory_budget_mb = NULL
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
  verbose = FALSE,
  seq_idx = NULL,
  seq_file = NULL,
  clustering_memory_budget_mb = NULL
) {
  if (!is.null(seq_file)) {
    stop(
      "`seq_file` is not supported for `data.frame` `seq` inputs.",
      call. = FALSE
    )
  }
  if (!is.null(seq_idx)) {
    ii <- seq_resolve_linear_seq_idx(nrow(seq), seq_idx)
    seq <- if (length(ii) < 1L) {
      seq[integer(), , drop = FALSE]
    } else {
      seq[ii, , drop = FALSE]
    }
  }
  mycall <- match.call()
  mycall$seq_idx <- NULL
  mycall$seq_file <- NULL
  if (missing(seq_id)) {
    newseq_id <- quote(seq$seq_id)
    newseq_id[[2]] <- mycall$seq
    mycall$seq_id <- newseq_id
  }
  if (identical(dist_config$method, "usearch")) {
    if (!is.null(clustering_memory_budget_mb)) {
      stop(
        "`clustering_memory_budget_mb` is not supported with `dist_usearch()`.",
        call. = FALSE
      )
    }
    mycall[[1]] <- quote(seq_cluster_usearch.DNAStringSet)
    newseq <- quote(Biostrings::DNAStringSet(seq$seq))
    newseq[[2]][[2]] <- mycall$seq
    mycall$usearch <- dist_config$usearch
    mycall$usearch_ncpu <- dist_config$usearch_ncpu
    mycall$dist_config <- NULL
    mycall$clustering_memory_budget_mb <- NULL
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
  verbose = FALSE,
  seq_idx = NULL,
  seq_file = NULL,
  clustering_memory_budget_mb = NULL
) {
  checkmate::assert_class(dist_config, "optimotu_dist_config")
  if (identical(dist_config$method, "usearch")) {
    if (!is.null(clustering_memory_budget_mb)) {
      stop(
        "`clustering_memory_budget_mb` is not supported with `dist_usearch()`.",
        call. = FALSE
      )
    }
    mycall <- match.call()
    mycall$dist_config <- NULL
    mycall$clustering_memory_budget_mb <- NULL
    mycall$usearch <- dist_config$usearch
    mycall$usearch_ncpu <- dist_config$usearch_ncpu
    mycall[[1]] <- quote(optimotu::seq_cluster_usearch)
    return(eval(mycall, envir = parent.frame()))
  }
  if (identical(dist_config$method, "file")) {
    if (!is.null(seq_file)) {
      stop("`seq_file` is not supported with `dist_file()`.", call. = FALSE)
    }
    if (!is.null(clustering_memory_budget_mb)) {
      stop(
        "`clustering_memory_budget_mb` is not supported with `dist_file()`.",
        call. = FALSE
      )
    }
    names_vec <- if (missing(seq_id)) {
      seq_input_seq_ids(seq, seq_idx, NULL)
    } else {
      ii <- seq_resolve_linear_seq_idx(length(seq_id), seq_idx)
      seq_id[ii]
    }
    mycall <- match.call()
    mycall$dist_config <- NULL
    mycall$seq <- NULL
    mycall$names <- names_vec
    mycall$seq_id <- NULL
    mycall$seq_idx <- NULL
    mycall$seq_file <- NULL
    mycall$clustering_memory_budget_mb <- NULL
    mycall$by_names <- dist_config$by_names
    mycall[[1]] <- quote(optimotu::distmx_cluster)
    return(eval(mycall, envir = parent.frame()))
  }
  output_type = match.arg(output_type)
  seq <- seq_as_char(
    seq,
    seq_idx = seq_idx,
    seq_file = seq_file,
    as = "character"
  )
  if (!missing(seq_id)) {
    names(seq) <- seq_id
  }
  checkmate::assert_class(threshold_config, "optimotu_threshold_config")
  checkmate::assert_class(clust_config, "optimotu_cluster_config")
  checkmate::assert_class(parallel_config, "optimotu_parallel_config")
  checkmate::assert_character(seq_id)
  checkmate::assert_number(
    clustering_memory_budget_mb,
    lower = 0,
    null.ok = TRUE,
    finite = TRUE
  )
  out <- if (!is.list(which)) {
    if (!(length(seq) == 0L && isTRUE(which))) {
      seq <- seq[which]
    }
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
      checkmate::check_list(
        which,
        types = "logical",
        any.missing = FALSE,
        min.len = 1
      ),
      checkmate::check_list(
        which,
        types = "integerish",
        any.missing = FALSE,
        min.len = 1
      ),
      checkmate::check_list(
        which,
        types = "character",
        any.missing = FALSE,
        min.len = 1
      ),
    )
    if (is.list(which) && !is.character(which[[1]])) {
      which <- lapply(which, `[`, x = seq_id)
    }
    if (isTRUE(verbose) || verbose >= 1L) {
      cat("Verifying subsets...")
    }
    verify_which(which, seq_id)
    if (isTRUE(verbose) || verbose >= 1L) {
      cat(" done\n")
    }
    # Restrict to the union of subset members so packing and pair
    # generation do not scale with unused sequences.
    keep <- sort(unique(unlist(which, use.names = FALSE)))
    seq <- seq[keep]
    seq_cluster_multi(
      seq = seq,
      which = which,
      dist_config = dist_config,
      threshold_config = threshold_config,
      clust_config = clust_config,
      parallel_config = parallel_config,
      output_type = output_type,
      verbose = as.integer(verbose),
      clustering_memory_budget_mb = if (is.null(clustering_memory_budget_mb)) {
        -1
      } else {
        clustering_memory_budget_mb
      }
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
  verbose = FALSE,
  seq_idx = NULL,
  seq_file = NULL,
  clustering_memory_budget_mb = NULL
) {
  if (!is.null(seq_file)) {
    stop(
      "`seq_file` is not supported for `DNAStringSet` `seq` inputs.",
      call. = FALSE
    )
  }
  if (!is.null(seq_idx)) {
    ii <- seq_resolve_linear_seq_idx(length(seq), seq_idx)
    seq <- if (length(ii) < 1L) {
      seq[integer()]
    } else {
      seq[ii]
    }
  }
  mycall <- match.call()
  mycall$seq_idx <- NULL
  mycall$seq_file <- NULL
  if (identical(dist_config$method, "usearch")) {
    if (!is.null(clustering_memory_budget_mb)) {
      stop(
        "`clustering_memory_budget_mb` is not supported with `dist_usearch()`.",
        call. = FALSE
      )
    }
    mycall$dist_config <- NULL
    mycall$clustering_memory_budget_mb <- NULL
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
