#' Search for the closest match(es) to sequences in a reference database
#'
#' @param seq (`character`, `data.frame`, or `XStringSet`) DNA sequences to
#' calculate distances for
#' @param seq_id (`character` vector) names for the sequences.  If they
#' are already named, this will replace the names.
#' @param threshold (`numeric` scalar) maximum distance to consider a match, in
#' \[0, 1\] where 0 is identical.
#' @param dist_config (`optimotu_dist_config`) configuration for calculating
#' distances, as returned by `dist_config()` or its helpers.
#' @param parallel_config (`optimotu_parallel_config`) configuration for
#' parallel processing, as returned by `parallel_config()` or its helpers.
#' @param verbose (`logical` or `integer` scalar) print progress messages.
#' @param details (one of "none", "gapstats", or "cigar") if "gapstats", return
#' #' the gap statistics. If "cigar", return the CIGAR string. Otherwise return
#' only the score and distance
#' @param span (`character` string) the span of the alignment; currently
#' accepted values are "global" and "extension".  The default is "global".
#' @param constrain (`logical` flag) if `TRUE`, the alignment algorithm will
#' use optimizations that will cause it to exit early if all possible alignments
#' have a distance greater than the distance threshold. This should not change
#' the correctness of distance calculations below the threshold, and results in
#' a large speedup for WFA2 and edlib. It is recommended to use
#' `constrain=FALSE` only to verify that the results do not change.
#' @param id_is_int (`logical` scalar) if `TRUE`, the sequence IDs are
#' parsed as integers, and the returned IDs are integers.
#' The default is `FALSE`.
#' @param ... passed to methods
#' @return (`data.frame`) with columns "seq_id1" and "seq_id2" (`character`),
#' or "seq_idx1" and "seq_idx2" (`integer`) if `id_is_integer` is `TRUE`,
#' "score1" and "dist1" giving the score and distance for the prealignment
#' stage, if any, and "score2" and "dist2" (`numeric`) giving the final score
#' and distance. If `details` is "gapstats", the gap statistics are returned in
#' columns "align_length", "n_insert", "n_delete", "max_insert", and
#' "max_delete" (`integer`). If `details` is "cigar", the CIGAR string is
#' returned in a column "cigar" (`character`).
#' @export

seq_distmx <- function(
  seq,
  threshold,
  seq_id = NULL,
  dist_config = dist_wfa2(),
  parallel_config = parallel_concurrent(1),
  verbose = FALSE,
  details = c("none", "gapstats", "cigar"),
  span = c("global", "extension"),
  constrain = TRUE,
  id_is_int = is.data.frame(seq) && "seq_idx" %in% names(seq),
  ...
) {
  checkmate::assert_character(seq_id, null.ok = TRUE, len = length(seq),
                              unique = TRUE)
  checkmate::assert_class(dist_config, "optimotu_dist_config")
  checkmate::assert_class(parallel_config, "optimotu_parallel_config")
  if (!missing(details)) checkmate::assert_string(details)
  details <- match.arg(details, several.ok = FALSE)
  details <- match(details, c("none", "gapstats", "cigar")) - 1L
  if (!missing(span)) checkmate::assert_string(span)
  span <- match.arg(span, several.ok = FALSE)
  span <- match(span, c("global", "extension")) - 1L
  checkmate::assert_flag(id_is_int)

  mycall <- match.call()

  if (identical(dist_config$method, "usearch")) {
    mycall[[1]] <- quote(optimotu::seq_distmx_usearch)
    mycall$dist_config <- NULL
    mycall$span <- NULL
    mycall$constrain <- NULL
    mycall$usearch <- dist_config$usearch
    eval(mycall, envir = parent.frame())
  } else {
    seq <- seq_as_char(seq)
    if (!is.null(seq_id)) {
      names(seq) <- seq_id
    }

    if (identical(dist_config$method, "file")) {
      checkmate::assert_file(dist_config$filename)
      warning("seq_distmx called with external distance matrix --",
              " ignoring 'details' and 'span'")
      if (dist_config$by_name == TRUE) {
        out <- utils::read.table(
          file = dist_config$filename,
          header = FALSE,
          col.names = c("seq_id1", "seq_id2", "dist"),
          colClasses = c("character", "character", "numeric")
        )
        out <- out[out$seq_id1 %in% names(seq), ]
        out <- out[out$seq_id2 %in% names(seq), ]
        out <- out[out$dist <= threshold, ]
        out
      } else {
        utils::read.table(
          file = dist_config$filename,
          header = FALSE,
          col.names = c("seq_idx1", "seq_idx2", "dist"),
          colClasses = c("integer", "integer", "numeric")
        )
        out <- out[out$seq_idx1 + 1 <= length(seq), ]
        out <- out[out$seq_idx2 + 1 <= length(seq), ]
        if (id_is_int) {
          out$seq_idx1 <- out$seq_idx1 + 1
          out$seq_idx2 <- out$seq_idx2 + 1
        } else {
          out$seq_id1 <- names(seq)[out$seq_idx1 + 1]
          out$seq_id2 <- names(seq)[out$seq_idx2 + 1]
          out$seq_idx1 <- NULL
          out$seq_idx2 <- NULL
        }
        out <- out[out$dist <= threshold, ]
        out
      }
    } else {
      out <- seq_distmx_internal(
        seq = seq,
        dist_config = dist_config,
        parallel_config = parallel_config,
        threshold = threshold,
        verbose = verbose,
        details = details,
        span = span,
        constrain = constrain
      )
      if (isFALSE(id_is_int)) {
        out$seq_idx1 <- names(seq)[out$seq_idx1 + 1]
        out$seq_idx2 <- names(seq)[out$seq_idx2 + 1]
        names(out)[1:2] <- c("seq_id1", "seq_id2")
      }
      out
    }
  }
}

#' @export
#' @rdname seq_distmx
seq_distmx_edlib <- function(
  seq,
  threshold,
  details = 0L,
  span = 0L,
  constrain = TRUE,
  threads = 1L,
  verbose = 0L,
  ...
) {
  call <- match.call()
  call[[1]] <- quote(optimotu::seq_distmx)
  call$dist_config <- dist_edlib()
  call$parallel_config <- parallel_concurrent(threads)
  call$threads <- NULL
  eval(call, envir = parent.frame())
}

#' @param min_overlap (`integer` scalar) minimum length of overlap between two
#' sequences (i.e., length after subtracting end gaps in either sequence). If
#' the required overlap is not met, the distance is returned as 1.
#' @param ignore_gaps (`logical` flag) if `TRUE`, columns containing an internal
#' gap in one sequence are ignored, i.e. not counted as a mismatch and also not
#' counted towards the total aligned length. If `FALSE`, such columns are
#' counted as mismatches.  End gaps in either sequence, and internal gaps shared
#' by both sequences, are always ignored.
#' @export
#' @rdname seq_distmx
seq_distmx_hamming <- function(
  seq,
  threshold,
  min_overlap = 0L,
  ignore_gaps = FALSE,
  threads = 1L,
  verbose = 0L,
  ...
) {
  call <- match.call()
  call[[1]] <- quote(optimotu::seq_distmx)
  call$dist_config <- dist_hamming(
    min_overlap = min_overlap,
    ignore_gaps = ignore_gaps
  )
  call$parallel_config <- parallel_concurrent(threads)
  call$min_overlap <- NULL
  call$ignore_gaps <- NULL
  call$threads <- NULL
  eval(call, envir = parent.frame())
}

#' @param match (non-negative `integer` scalar) score for a match
#' @param mismatch (positive `integer` scalar) penalty for a mismatch
#' @param gap_open (non-negative `integer` scalar) penalty for gap opening
#' @param gap_extend (positive `integer` scalar) penalty for gap extension
#' @param gap_open2 (non-negative `integer` scalar) alternate penalty for gap
#' opening
#' @param gap_extend2 (positive `integer` scalar) alternate penalty for gap
#' extension
#' @export
#' @rdname seq_distmx
seq_distmx_wfa2 <- function(
  seq,
  threshold,
  match = -1L,
  mismatch = 2L,
  gap_open = 10L,
  gap_extend = 1L,
  gap_open2 = 0L,
  gap_extend2 = 0L,
  constrain = TRUE,
  threads = 1L,
  verbose = 0L,
  ...
) {
  call <- match.call()
  call[[1]] <- quote(optimotu::seq_distmx)
  call$dist_config <- dist_wfa2(
    match = match,
    mismatch = mismatch,
    gap_open = gap_open,
    gap_extend = gap_extend,
    gap_open2 = gap_open2,
    gap_extend2 = gap_extend2
  )
  call$parallel_config <- parallel_concurrent(threads)
  call$match <- NULL
  call$mismatch <- NULL
  call$gap_open <- NULL
  call$gap_extend <- NULL
  call$gap_open2 <- NULL
  call$gap_extend2 <- NULL
  call$threads <- NULL
  eval(call, envir = parent.frame())
}