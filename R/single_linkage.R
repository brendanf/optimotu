verify_which <- function(which, method, seqnames) {
   checkmate::assert_list(
      which,
      types = "character",
      any.missing = FALSE
   )
   if (!all(vapply(lapply(which, `%in%`, table = seqnames), all, TRUE))) {
      stop("Elements of 'which' must all be in 'seqnames'.")
   }
   if (method=="matrix") {
      warning("Method 'matrix' is not yet implemented for subset clustering.",
              "Falling back to method 'tree'.")
   }
   invisible(TRUE)
}

verify_threshold_steps <- function(thresh_min, thresh_max, thresh_step) {
   if (is.null(thresh_min) || is.null(thresh_max) || is.null(thresh_step)) {
      stop("either 'thresholds' or 'thresh_min', 'thresh_max', and 'thresh_step' must be given")
   }
   if (!is.numeric(thresh_min) || !is.numeric(thresh_max) || !is.numeric(thresh_step)) {
      stop("'thresh_min', 'thresh_max', and 'thresh_step' must all be numbers.")
   }
   if (length(thresh_min) != 1L || length(thresh_max) != 1L || length(thresh_step) != 1L) {
      stop("'thresh_min', 'thresh_max', and 'thresh_step' must all be of length 1.")
   }
   if (is.na(thresh_min) || is.na(thresh_max) || is.na(thresh_step)) {
      stop("'thresh_min', 'thresh_max', and 'thresh_step' may not be NA.")
   }
   if (is.nan(thresh_min) || is.nan(thresh_max) || is.nan(thresh_step)) {
      stop("'thresh_min', 'thresh_max', and 'thresh_step' may not be NaN.")
   }
   if (thresh_step <= 0) {
      stop("'thresh_step' must be positive.")
   }
   if (thresh_max < thresh_min) {
      stop("'thresh_max' must be greater than or equal to 'thresh_min'.")
   }
}

verify_thresholds <- function(thresholds) {
   if (!is.numeric(thresholds)) {
      stop("'thresholds' must be numeric.")
   }
   if (length(thresholds) == 0) {
      stop("At least one threshold must be given in 'thresholds'.")
   }
   if (any(is.na(thresholds))) {
      stop("'thresholds' may not have NA values.")
   }
   if (any(is.nan(thresholds))) {
      stop("'thresholds' may not have NaN values.")
   }
   if (any(thresholds < 0)) {
      stop("'thresholds' may not contain negative values.")
   }
   if (any(utils::head(thresholds, -1) >= utils::tail(thresholds, -1))) {
      stop("'thresholds' must be strictly increasing.")
   }
   invisible(TRUE)
}

verify_precision <- function(precision) {
   checkmate::assert_number(
      precision,
      na.ok = FALSE,
      null.ok = TRUE,
      finite = TRUE,
      lower = .Machine$double.xmin
   )
}

verify_method_output_type <- function(
   method = c("tree", "matrix"),
   output_type = c("matrix", "hclust")
) {
   method = match.arg(method)
   output_type = match.arg(output_type)
   if (output_type == "hclust" && method == "matrix") {
      stop("output_type 'hclust' is not available for method 'matrix'.")
   }
   invisible(TRUE)
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
#' @param method (`character`) The algorithm to use; one of "tree" or "matrix".
#' The "tree" algorithm is faster in at least some large cases, but tends to be
#' slower in smaller cases, and cannot take advantage of parallel computation
#' unless multiple overlapping subsets are specified in `which`. The two
#' algorithms give identical results.
#' @param output_type (`character`) Which type of output to give; one of
#' `"matrix"` or `"hclust"`. `"matrix"` returns an integer matrix giving
#' clustering results, where the element in row `i` and column `j` gives the
#' 0-based index of the first member of the cluster to which sequence `j`
#' belongs when clustered at the `i`th clustering threshold. `"hclust"` returns
#' an object as returned by [stats::hclust()], which requires less memory,
#' especially for large problems, but is only supported for method `"tree"`. If
#' `which` is given, then either `output_type` returns a list whose elements are
#' of the chosen type.
#' @param thresholds (sorted `numeric` vector) An explicit list of clustering
#' thresholds to try.  These do not need to be evenly spaced but must be
#' strictly increasing.
#' @param precision (`numeric` scalar) The precision of the distances in the
#' distance matrix; providing this may give a slight speedup when explicit
#' `thresholds` are provided. If the actual precision of numbers in the distance
#' matrix is smaller than this value, then distances will be rounded to this
#' precision without warning.
#' @param thresh_min (`numeric` scalar) The minimum distance threshold for
#' clustering; should not be given if explicit `thresholds` are specified.
#' @param thresh_max (`numeric` scalar) The maximum distance threshold for
#' clustering; should not be given if explicit `thresholds` are specified.
#' @param thresh_step (`numeric` scalar) The spacing between subsequent distance
#' thresholds for clustering; should not be given if explicit `thresholds` are
#' specified.
#' @param which (`list` of `character` vectors) Instead of performing clustering
#' on all input sequences, perform independent clustering on subsets of the
#' sequences defined by the elements of `which`. Subsets do not need to be
#' disjoint (and indeed, if they are it is probably faster to calculate the
#' distance matrices separately.) Currently `which` is only implemented for the
#' `"tree"` algorithm.
#' @param threads (`integer` scalar) Maximum number of parallel threads.
#' @param minsplit (`integer` scalar) Controls the granularity of parallel
#' processing in the "matrix" algorithm.
#'
#' @return An [`integer matrix`][methods::StructureClasses] if
#' `output_type=="matrix"`, an [`hclust`][stats::hclust] object if
#' `output_type=="hclust"`, or a list of one of these when `which` is a list.
#' @export
#'
#' @examples
single_linkage = function(
   distmx,
   names,
   method = c("tree", "matrix"),
   output_type = c("matrix", "hclust"),
   thresh_min = NULL,
   thresh_max = NULL,
   thresh_step = NULL,
   thresholds = NULL,
   precision = NULL,
   which = NULL,
   threads = 1L,
   minsplit = 1L
) {
   method = match.arg(method)
   output_type = match.arg(output_type)
   verify_method_output_type()
   if (is.null(thresholds)) {
      if (!is.null(precision)) {
         warning("'precision' has no effect when 'thresholds' is not given. Ignoring.\n")
      }
      verify_threshold_steps(thresh_min, thresh_max, thresh_step)
      if (!is.null(which) && !isTRUE(which)) {
         verify_which(which, method, names)
         single_linkage_multi_uniform(
            distmx,
            names,
            thresh_min,
            thresh_max,
            which,
            threads
         )
      } else {
         switch(
            method,
            tree = single_linkage_pool_uniform(
               distmx,
               names,
               thresh_min,
               thresh_max,
               thresh_step,
               output_type
            ),
            matrix = single_linkage_matrix_uniform(
               distmx,
               names,
               thresh_min,
               thresh_max,
               thresh_step,
               threads,
               minsplit
            )
         )
      }
   } else if (!is.null(thresh_min) || !is.null(thresh_max) || !is.null(thresh_step)) {
      stop("If 'thresholds' is given, 'thresh_min', 'thresh_max', and 'thresh_step' must not be.")
   } else {
      verify_thresholds(thresholds)
      verify_precision(precision)
      if (!is.null(which) && !isTRUE(which)) {
         verify_which(which, method, seqnames)
         if (is.null(precision)) {
            single_linkage_multi_array(distmx, names, thresholds, which, threads)
         } else {
            single_linkage_multi_cached(distmx, names, thresholds, precision, which, threads)
         }
      } else {
         if (is.null(precision)) {
            switch(
               method,
               tree = single_linkage_pool_array(distmx, names, thresholds, output_type),
               matrix = single_linkage_matrix_array(distmx, names, thresholds, threads, minsplit)
            )
         } else {
            switch(
               method,
               tree = single_linkage_pool_cached(distmx, names, thresholds, precision, output_type),
               matrix = single_linkage_matrix_cached(distmx, names, thresholds, precision, threads, minsplit)
            )
         }
      }
   }
}

#' Do single-linkage clustering at a series of increasing similarity thresholds
#' using USEARCH to calculate sequence similarities.
#'
#' @description This function uses a unix pipe to direct the output of the
#' USEARCH "`calc_distmx`" command  to [single_linkage()]. USEARCH (version 8.0
#' or higher) should be installed separately; it is available with a free
#' license for most users at [https://www.drive5.com/usearch/].
#'
#' @param seq (`character` vector, filename,
#' [DNAStringSet][Biostrings::DNAStringSet()], or `data.frame` with columns
#' "seq_id" (`character`) and "seq" (`character`)) sequences to cluster
#' @param seq_id (`character` vector) names for the sequences.  If they are
#' already named, this will replace the names.  Has no effect if `seq` is a
#' filename.
#' @param thresh_min (`numeric`) minimum sequence dissimilarity threshold for
#' clustering. Number between 0 and 1.
#' @param thresh_max (`numeric`) maximum sequence dissimilarity threshold for
#' clustering. Number between 0 and 1.
#' @param thresh_step (`numeric`) difference between successive percentage similarity thresholds for
#' clustering. Number <= thresh_max - thresh_min.
#' @param thresh_name (`character` vector) names for the thresholds.
#' @param which (`logical`, `character` or `integer` vector, or a list of these)
#' subset(s) of `seq` to operate on.  If `seq` is a filename, then the distance
#' matrix for the full file will be calculated, but only the selected sequences
#' will be clustered.  Otherwise, distances will only be calculated for
#' sequences occurring in the subset(s)
#' @param ncpu (`integer` scalar) number of threads to use for calculating the
#' distance matrix and clustering
#' @param usearch (`character` scalar) path to usearch executable
#'
#' @return `integer` matrix giving clustering results
#' @export
usearch_single_linkage <- function(
   seq,
   seq_id = names(seq),
   method = c("tree", "matrix"),
   output_type = c("matrix", "hclust"),
   thresh_max = NULL, thresh_min = NULL, thresh_step = NULL,
   thresholds = NULL,
   precision = if (is.null(thresholds)) NULL else 0.001,
   thresh_names = names(thresholds),
   which = TRUE,
   ncpu = local_cpus(),
   usearch = Sys.which("usearch")) {
   UseMethod("usearch_single_linkage", seq)
}

#' @method usearch_single_linkage data.frame
#' @export
usearch_single_linkage.data.frame <- function(
   seq,
   seq_id = seq$seq_id,
   method = c("tree", "matrix"),
   output_type = c("matrix", "hclust"),
   thresh_max = NULL, thresh_min = NULL, thresh_step = NULL,
   thresholds = NULL,
   precision = NULL,
   thresh_names = names(thresholds),
   which = TRUE,
   ncpu = local_cpus(),
   usearch = Sys.which("usearch")
) {
   mycall <- match.call()
   mycall[[1]] <- usearch_single_linkage.DNAStringSet
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
usearch_single_linkage.character <- function(
   seq,
   seq_id = names(seq),
   method = c("tree", "matrix"),
   output_type = c("matrix", "hclust"),
   thresh_max = NULL, thresh_min = NULL, thresh_step = NULL,
   thresholds = NULL,
   precision = NULL,
   thresh_names = names(thresholds),
   which = TRUE,
   ncpu = local_cpus(),
   usearch = Sys.which("usearch")
) {
   verify_method_output_type()
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
      do_usearch_single_linkage(
         seq_file = tf,
         seq_id = names(index),
         method = method,
         output_type = output_type,
         thresh_max = thresh_max,
         thresh_min = thresh_min,
         thresh_step = thresh_step,
         thresholds = thresholds,
         precision = precision,
         thresh_names = thresh_names,
         which = which,
         ncpu = ncpu,
         usearch = usearch
      )
   } else {
      mycall <- match.call()
      mycall[[1]] <- usearch_single_linkage.DNAStringSet
      newseq <- quote(Biostrings::DNAStringSet(seq))
      newseq[[2]] <- mycall$seq
      mycall$seq <- newseq
      eval(mycall, envir = parent.frame())
   }
}

#' @export
usearch_single_linkage.DNAStringSet <- function(
   seq,
   seq_id = names(seq),
   method = c("tree", "matrix"),
   output_type = c("matrix", "hclust"),
   thresh_max = NULL, thresh_min = NULL, thresh_step = NULL,
   thresholds = NULL,
   precision = NULL,
   thresh_names = names(thresholds),
   which = TRUE,
   ncpu = local_cpus(),
   usearch = Sys.which("usearch")
) {
   verify_method_output_type()
   # rename the sequences if necessary
   if (!isTRUE(all.equal(names(seq), seq_id))) names(seq) <- seq_id
   if (is.list(which)) {
      if (all(vapply(which, is.logical, TRUE))) {
         seq <- seq[Reduce(magrittr::or, which)]
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
      method = method,
      output_type = output_type,
      thresh_max = thresh_max,
      thresh_min = thresh_min,
      thresh_step = thresh_step,
      thresholds = thresholds,
      precision = precision,
      thresh_names = thresh_names,
      which = which,
      ncpu = ncpu,
      usearch = usearch
   )
}

do_usearch_singlelink <- function(
   seq_file,
   seq_id,
   method,
   output_type,
   thresh_max, thresh_min, thresh_step,
   thresholds,
   precision,
   thresh_names,
   which,
   ncpu,
   usearch = Sys.which("usearch")
) {

   if (is.list(thresholds)) {
      verify_thresholds(thresholds)
      verify_precision(precision)
      nthresh <- length(thresholds)
   } else {
      verify_threshold_steps(thresh_min, thresh_max, thresh_step)
      nthresh <- length(seq(thresh_min, thresh_max, thresh_step))
   }
   if (is.list(which)) {
      verify_which(which, seq_id)
   }
   checkmate::assert(
      checkmate::check_null(thresh_names),
      checkmate::check_character(thresh_names, len = nthresh)
   )
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
   method <- match.arg(method)
   fifoname <- tempfile("fifo")
   stopifnot(system2("mkfifo", fifoname) == 0)
   on.exit(unlink(fifoname), TRUE)
   f <- fifo(fifoname)
   system2(
      usearch,
      c(
         "-calc_distmx", seq_file, # input file
         "-tabbedout", fifoname, # output fifo
         "-maxdist", thresh_max, # similarity threshold
         "-termdist", min(1, 1.5*thresh_max), # threshold for udist
         "-lopen", "1", # gap opening
         "-lext", "1", # gap extend
         # "-pattern", "111010010111", # pattern gives better result than kmers maybe?
         "-threads", ncpu
      ),
      wait = FALSE
   )
   if (is.list(which)) {
      out <- single_linkage(
         distmx = fifoname,
         seq_id,
         thresh_min = thresh_min,
         thresh_max = thresh_max,
         thresh_step = thresh_step,
         thresholds = thresholds,
         precision = precision,
         method = method,
         output_type = output_type,
         which = which,
         threads = ncpu
      )
      out <- lapply(out, `rownames<-`, thresh_names)
   } else {
      out <- single_linkage(
         distmx = fifoname,
         seq_id,
         thresh_min = thresh_min,
         thresh_max = thresh_max,
         thresh_step = thresh_step,
         thresholds = thresholds,
         precision = precision,
         method = method,
         output_type = output_type,
         threads = ncpu
      )
      row.names(out) <- thresh_names
   }
   out
}
