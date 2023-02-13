verify_which <- function(which, seqnames) {
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

verify_threshold_steps <- function(dmin, dmax, dstep) {
   if (is.null(dmin) || is.null(dmax) || is.null(dstep)) {
      stop("either 'thresholds' or 'dmin', 'dmax', and 'dstep' must be given")
   }
   if (!is.numeric(dmin) || !is.numeric(dmax) || !is.numeric(dstep)) {
      stop("'dmin', 'dmax', and 'dstep' must all be numbers.")
   }
   if (length(dmin) != 1L || length(dmax) != 1L || length(dstep) != 1L) {
      stop("'dmin', 'dmax', and 'dstep' must all be of length 1.")
   }
   if (is.na(dmin) || is.na(dmax) || is.na(dstep)) {
      stop("'dmin', 'dmax', and 'dstep' may not be NA.")
   }
   if (is.nan(dmin) || is.nan(dmax) || is.nan(dstep)) {
      stop("'dmin', 'dmax', and 'dstep' may not be NaN.")
   }
   if (dstep <= 0) {
      stop("'dstep' must be positive.")
   }
   if (dmax < dmin) {
      stop("'dmax' must be greater than or equal to 'dmin'.")
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

#' Title
#'
#' @param file (`character` filename) The name of the file (or, e.g., a named
#' pipe) from which to read the sparse distance matrix. The matrix should be
#' a whitespace delimited text file, with three values per line: seq1 seq2 dist,
#' where seq1 and seq2 are numeric indices of the sequences to be clustered,
#' starting with 0, and dist is the distance, typically a read number between 0
#' and 1 (but the algorithm works for any non-negative distance.)
#' @param seqnames (`character` vector) Names of the sequences, used for
#' labeling the columns of the output matrix. Typically, in order to generate
#' the sparse distance matrix with integer indices, an alternate version of the
#' input file is written, where the sequence names are replaced with integers.
#' These are the "real" names.
#' @param method (`character`) The algorithm to use; one of "tree" or "matrix".
#' The "tree" algorithm is faster in at least some large cases, but tends to be
#' slower in smaller cases, and cannot take advantage of parallel computation
#' unless multiple overlapping subsets are specified in `which`. The two
#' algorithms give identical results.
#' @param output_type (`character`) Which type of output to give; one of
#' "matrix" or "hclust". "matrix" returns an integer matrix giving clustering
#' results, where the element in row `i` and column `j` give the 0-based index
#' of the first member of the cluster to which sequence `j` belongs when
#' clustered at the `i`th clustering threshold. "hclust" returns an object as
#' returned by `stats::hclust()`, which requires less memory, especially for
#' large problems, but is only supported for method "tree". If "which" is given,
#' then either output_type returns a list whose elements are of the chosen type.
#' @param thresholds (sorted `numeric` vector) An explicit list of clustering
#' thresholds to try.  These do not need to be evenly spaced but must be
#' strictly increasing.
#' @param precision (`numeric` scalar) The precision of the distances in the
#' distance matrix; providing this may give a slight speedup when explicit
#' `thresholds` are provided. If the actual precision of numbers in the distance
#' matrix is smaller than this value, then distances will be rounded to this
#' precision without warning.
#' @param which (`list` of `character` vectors) Instead of performing clustering
#' on all input sequences, perform independent clustering on subsets of the
#' sequences defined by the elements of `which`. Subsets do not need to be
#' disjoint (and indeed, if they are it is probably faster to calculate the
#' distance matrices separately.) Currently `which` is only implemented for the
#' "tree" algorithm.
#' @param dmin (`numeric` scalar) The minimum distance threshold for clustering;
#' should not be given if explicit `thresholds` are specified.
#' @param dmax (`numeric` scalar) The maximum distance threshold for clustering;
#' should not be given if explicit `thresholds` are specified.
#' @param dstep (`numeric` scalar) The spacing between subsequent distance
#' thresholds for clustering; should not be given if explicit `thresholds` are
#' specified.
#' @param threads (`integer` scalar) Maximum number of parallel threads.
#' @param minsplit (`integer` scalar) Controls the granularity of parallel
#' processing in the "matrix" algorithm.
#'
#' @return An [`integer matrix`][methods::StructureClasses] if
#' `output_type=="matrix"`, an [`hclust`][stats::hclust] object if
#' `output_type == "hclust"`, or a list of one of these when `which` is a list.
#' @export
#'
#' @examples
single_linkage = function(
   file,
   seqnames,
   method = c("tree", "matrix"),
   output_type = c("matrix", "hclust"),
   thresholds = NULL,
   precision = NULL,
   which = NULL,
   dmin = NULL,
   dmax = NULL,
   dstep = NULL,
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
      verify_threshold_steps(dmin, dmax, dstep)
      if (!is.null(which) && !isTRUE(which)) {
         verify_which(which, method, seqnames)
         single_linkage_multi_uniform(file, seqnames, dmin, dmax, which, threads)
      } else {
         switch(
            method,
            tree = single_linkage_pool_uniform(file, seqnames, dmin, dmax, dstep, output_type),
            matrix = single_linkage_matrix_uniform(file, seqnames, dmin, dmax, dstep, threads, minsplit)
         )
      }
   } else if (!is.null(dmin) || !is.null(dmax) || !is.null(dstep)) {
      stop("If 'thresholds' is given, 'dmin', 'dmax', and 'dstep' must not be.")
   } else {
      verify_thresholds(thresholds)
      verify_precision(precision)
      if (!is.null(which) && !isTRUE(which)) {
         verify_which(which, method, seqnames)
         if (is.null(precision)) {
            single_linkage_multi_array(file, seqnames, thresholds, which, threads)
         } else {
            single_linkage_multi_cached(file, seqnames, thresholds, precision, which, threads)
         }
      } else {
         if (is.null(precision)) {
            switch(
               method,
               tree = single_linkage_pool_array(file, seqnames, thresholds, output_type),
               matrix = single_linkage_matrix_array(file, seqnames, thresholds, threads, minsplit)
            )
         } else {
            switch(
               method,
               tree = single_linkage_pool_cached(file, seqnames, thresholds, precision, output_type),
               matrix = single_linkage_matrix_cached(file, seqnames, thresholds, precision, threads, minsplit)
            )
         }
      }
   }
}

#' Do single-linkage clustering at a series of increasing similarity thresholds.
#'
#' @param seq (`character` vector, filename, or `Biostrings::DNAStringSet`)
#' sequences to cluster
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
#' @param ...
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
      mycall[[1]] <- usearch_singlelink.DNAStringSet
      newseq <- quote(Biostrings::DNAStringSet(seq))
      newseq[[2]] <- mycall$seq
      mycall$seq <- newseq
      eval(mycall, envir = parent.frame())
   }
}

usearch_singlelink.DNAStringSet <- function(
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
         file = fifoname,
         seq_id,
         dmin = thresh_min,
         dmax = thresh_max,
         dstep = thresh_step,
         thresholds = thresholds,
         method = method,
         output_type = output_type,
         preclust = which,
         threads = ncpu
      )
      out <- lapply(out, `rownames<-`, thresh_names)
   } else {
      out <- single_linkage(
         file = fifoname,
         seq_id,
         dmin = thresh_min,
         dmax = thresh_max,
         dstep = thresh_step,
         thresholds = thresholds,
         method = method,
         output_type = output_type,
         threads = ncpu
      )
      row.names(out) <- thresh_names
   }
   out
}
