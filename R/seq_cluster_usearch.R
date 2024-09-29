#' Do single-linkage clustering at a series of increasing similarity thresholds
#' using USEARCH to calculate sequence similarities.
#'
#' @description This function uses a unix pipe to direct the output of the
#' USEARCH "`calc_distmx`" command  to [distmx_cluster()]. USEARCH (version 8.0
#' or higher) should be installed separately; it is available with a free
#' license for most users at https://www.drive5.com/usearch/.
#'
#' @inheritParams distmx_cluster
#' @param seq (`character` vector, filename,
#' [DNAStringSet][Biostrings::DNAStringSet()], or `data.frame` with columns
#' "seq_id" (`character`) and "seq" (`character`)) sequences to cluster
#' @param seq_id (`character` vector) names for the sequences.  If they are
#' already named, this will replace the names.  Has no effect if `seq` is a
#' filename.
#' @param usearch_ncpu (`integer` scalar) number of threads to use for
#' calculating the distance matrix.  The number of threads for clustering is
#' specified in the "parallel_config" argument.
#' @param usearch (`character` scalar) path to usearch executable
#'
#' @return An [`integer matrix`][methods::structure-class] if
#' `output_type=="matrix"`, an [`hclust`][stats::hclust] object if
#' `output_type=="hclust"`, or a list of one of these when `which` is a list.
#'
#' @export
seq_cluster_usearch <- function(
    seq,
    seq_id = names(seq),
    threshold_config,
    clust_config = clust_tree(),
    parallel_config = parallel_concurrent(1),
    output_type = c("matrix", "hclust"),
    which = TRUE,
    usearch_ncpu = NULL,
    usearch = Sys.which("usearch")
) {
  UseMethod("seq_cluster_usearch", seq)
}

#' @method seq_cluster_usearch data.frame
#' @export
seq_cluster_usearch.data.frame <- function(
    seq,
    seq_id = seq$seq_id,
    threshold_config,
    clust_config = clust_tree(),
    parallel_config = parallel_concurrent(1),
    output_type = c("matrix", "hclust"),
    which = TRUE,
    usearch_ncpu = NULL,
    usearch = Sys.which("usearch")
) {
  mycall <- match.call()
  mycall[[1]] <- seq_cluster_usearch.DNAStringSet
  if (missing(seq_id)) {
    newseq_id <- quote(seq$seq_id)
    newseq_id[[2]] <- mycall$seq
    mycall$seq_id <- newseq_id
  }
  newseq <- quote(Biostrings::DNAStringSet(seq$seq))
  newseq[[2]][[2]] <- mycall$seq
  mycall$seq <- newseq
  eval(mycall, envir = parent.frame())
}

#' @export
seq_cluster_usearch.character <- function(
    seq,
    seq_id = names(seq),
    threshold_config,
    clust_config = clust_index(),
    parallel_config = parallel_concurrent(1),
    output_type = c("matrix", "hclust"),
    which = TRUE,
    usearch_ncpu = NULL,
    usearch = Sys.which("usearch")
) {
  output_type = match.arg(output_type)
  if (length(seq) == 1 && file.exists(seq)) {
    index <- Biostrings::fasta.seqlengths(seq)
    if (!missing(seq_id))
      warning("'seq_id' has no effect when 'seq' is a file.")
    if (!all(names(index) == as.character(seq_along(index)))) {
      # write a temp version of the file which has headers as integer indices
      tf <- tempfile(pattern = "clust", fileext = ".fasta")
      tc <- file(tf, open = "w")
      on.exit(close(tc))
      sc <- file(seq, open = "r")
      on.exit(close(sc), add = TRUE)
      nseq <- 0L
      while(length(lines <- readLines(sc, n = 10000)) > 0L) {
        headers <- grep("^>", lines)
        lines[headers] <- sprintf(">%d", seq_along(headers) + nseq - 1L)
        writeLines(lines, tc)
        nseq <- nseq + length(headers)
      }
      close(sc)
      close(tc)
      on.exit()
    } else {
      tf <- seq
    }
    do_usearch_singlelink(
      seq_file = tf,
      seq_id = names(index),
      threshold_config = threshold_config,
      clust_config = clust_config,
      parallel_config = parallel_config,
      output_type = output_type,
      which = which,
      usearch_ncpu = usearch_ncpu,
      usearch = usearch
    )
  } else {
    mycall <- match.call()
    mycall[[1]] <- seq_cluster_usearch.DNAStringSet
    newseq <- quote(Biostrings::DNAStringSet(seq))
    newseq[[2]] <- mycall$seq
    mycall$seq <- newseq
    eval(mycall, envir = parent.frame())
  }
}

#' @export
seq_cluster_usearch.DNAStringSet <- function(
    seq,
    seq_id = names(seq),
    threshold_config,
    clust_config = clust_index(),
    parallel_config = parallel_concurrent(1),
    output_type = c("matrix", "hclust"),
    which = TRUE,
    usearch_ncpu = NULL,
    usearch = Sys.which("usearch")
) {
  output_type = match.arg(output_type)
  # rename the sequences if necessary
  if (!isTRUE(all.equal(names(seq), seq_id))) names(seq) <- seq_id
  if (is.list(which)) {
    if (all(vapply(which, is.logical, TRUE))) {
      seq <- seq[Reduce(`|`, which)]
    } else if (all(vapply(which, is.integer, TRUE))) {
      seq <- seq[sort(unique(unlist(which)))]
    } else if (is_list_of_character(which)) {
      seq <- seq[sort(unique(unlist(which)))]
    } else {
      stop("'which' must be a character, integer, or logical vector, or a",
           " list of one of these.")
    }
  } else {
    seq <- seq[which]
  }
  # shortcut if only one sequence
  if (length(seq) == 1) return(names(seq))
  tf <- tempfile(pattern = "clust", fileext = ".fasta")
  Biostrings::writeXStringSet(`names<-`(seq, seq_along(seq) - 1L), tf)
  on.exit(unlink(tf))
  do_usearch_singlelink(
    seq_file = tf,
    seq_id = names(seq),
    threshold_config = threshold_config,
    clust_config = clust_config,
    parallel_config = parallel_config,
    output_type = output_type,
    which = which,
    usearch_ncpu = usearch_ncpu,
    usearch = usearch
  )
}

do_usearch_singlelink <- function(
    seq_file,
    seq_id,
    threshold_config,
    clust_config,
    parallel_config,
    output_type,
    which,
    usearch_ncpu,
    usearch
) {
  checkmate::assert_class(threshold_config, "optimotu_threshold_config")
  usearch_thresh_max <- if (threshold_config$type == "uniform") threshold_config$to else
    threshold_config$thresholds[length(threshold_config$thresholds)]
  checkmate::assert_class(clust_config, "optimotu_cluster_config")
  checkmate::assert_class(parallel_config, "optimotu_parallel_config")
  if (is.list(which)) {
    verify_which(which, seq_id)
  }
  checkmate::assert(
    checkmate::check_integer(which, lower = 1, upper = length(seq_id), any.missing = FALSE),
    checkmate::check_subset(which, seq_id),
    checkmate::check_logical(which, any.missing = FALSE),
    checkmate::check_list(which, types = "logical", any.missing = FALSE, min.len = 1),
    checkmate::check_list(which, types = "integerish", any.missing = FALSE, min.len = 1),
    checkmate::check_list(which, types = "character", any.missing = FALSE, min.len = 1)
  )
  if (is.list(which) && !is.character(which[[1]])) {
    which <- lapply(which, `[`, x = seq_id)
  }
  fifoname <- tempfile("fifo")
  stopifnot(system2("mkfifo", fifoname) == 0)
  on.exit(unlink(fifoname), TRUE)
  args <- c(
    "-calc_distmx", seq_file, # input file
    "-tabbedout", fifoname, # output fifo
    "-maxdist", usearch_thresh_max, # similarity threshold
    "-termdist", min(1, 2*usearch_thresh_max), # threshold for udist
    "-lopen", "1", # gap opening
    "-lext", "1" # gap extend
  )
  if (!is.null(usearch_ncpu)) {
    checkmate::assert_count(usearch_ncpu, positive = TRUE)
    args <- c(args, "-threads", usearch_ncpu)
  } else {
    usearch_ncpu <- 1L
  }
  if (system2(usearch, "-version", stdout = NULL, stderr = NULL) != 0) {
    stop("usearch could not be found at path: ", usearch)
  }
  system2(usearch, args, wait = FALSE)
  if (is.list(which)) {
    distmx_cluster(
      distmx = fifoname,
      names = seq_id,
      threshold_config = threshold_config,
      clust_config = clust_config,
      parallel_config = parallel_config,
      output_type = output_type,
      which = which
    )
  } else {
    distmx_cluster(
      distmx = fifoname,
      names = seq_id,
      threshold_config = threshold_config,
      clust_config = clust_config,
      parallel_config = parallel_config,
      output_type = output_type
    )
  }
}
