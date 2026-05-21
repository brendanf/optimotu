#' Search for the closest match(es) to sequences in a reference database
#'
#' @param query (`character`, `data.frame`, or `XStringSet`) sequences to search
#' @param query_id (`character` vector) names for the query sequences.  If they
#' are already named, this will replace the names.
#' @param ref (`character`, `data.frame`, or `XStringSet`) reference sequences to
#' search against
#' @param ref_id (`character` vector) names for the reference sequences.  If they
#' are already named, this will replace the names.
#' @param query_seq_idx,query_seq_file optional subsetting for `query`; see
#'   [seq_as_char()].
#' @param ref_seq_idx,ref_seq_file optional subsetting for `ref`; see
#'   [seq_as_char()].
#' @param threshold (`numeric` scalar) maximum distance to consider a match, in
#' \[0, 1\] where 0 is identical.
#' @param dist_config (`optimotu_dist_config`) configuration for calculating
#' distances, as returned by `dist_config()` or its helpers.
#' @param parallel_config (`optimotu_parallel_config`) configuration for parallel
#' processing, as returned by `parallel_config()` or its helpers.
#' @param verbose (`logical` or `integer` scalar) print progress messages.
#' @param return_cigar (`logical` scalar) if `TRUE`, return the cigar string
#' @param span (`character` string) the span of the alignment; currently
#' accepted values are "global" and "extension".  The default is "global".
#' @param ... passed to methods
#' @return (`data.frame`) with columns "seq_id" (`character`),
#' "ref_id" (`character`), and "dist" (`numeric`) giving the distance between
#' the query and reference. If `return_cigar` is `TRUE`, the CIGAR string is
#' returned in a column "cigar" (`character`).
#' @export

seq_search <- function(
  query,
  ref,
  threshold,
  query_id = NULL,
  ref_id = NULL,
  dist_config = dist_wfa2(),
  parallel_config = parallel_concurrent(1),
  verbose = FALSE,
  return_cigar = FALSE,
  span = c("global", "extension"),
  query_seq_idx = NULL,
  query_seq_file = NULL,
  ref_seq_idx = NULL,
  ref_seq_file = NULL,
  ...
) {
  nq <- seq_input_n_records(query, query_seq_file)
  nr <- seq_input_n_records(ref, ref_seq_file)
  if (!is.null(query_id)) {
    checkmate::assert_character(query_id, len = nq, unique = TRUE)
  }
  if (!is.null(ref_id)) {
    checkmate::assert_character(ref_id, len = nr, unique = TRUE)
  }
  checkmate::assert_class(dist_config, "optimotu_dist_config")
  checkmate::assert_class(parallel_config, "optimotu_parallel_config")
  checkmate::assert_flag(return_cigar)
  checkmate::assert_character(span)
  span <- match.arg(span, c("global", "extension"), several.ok = FALSE)
  checkmate::assert_choice(span, c("global", "extension"))
  span <- match(span, c("global", "extension")) - 1L

  mycall <- match.call()

  if (identical(dist_config$method, "usearch")) {
    mycall[[1]] <- quote(optimotu::seq_search_usearch)
    mycall$dist_config <- NULL
    mycall$usearch <- dist_config$usearch
    mycall$query_seq_idx <- query_seq_idx
    mycall$query_seq_file <- query_seq_file
    mycall$ref_seq_idx <- ref_seq_idx
    mycall$ref_seq_file <- ref_seq_file
    eval(mycall, envir = parent.frame())
  } else {
    query <- seq_as_char(
      query,
      seq_idx = query_seq_idx,
      seq_file = query_seq_file,
      as = "character"
    )
    if (!is.null(query_id)) {
      names(query) <- query_id
    }
    ref <- seq_as_char(
      ref,
      seq_idx = ref_seq_idx,
      seq_file = ref_seq_file,
      as = "character"
    )
    if (!is.null(ref_id)) {
      names(ref) <- ref_id
    }

    if (identical(dist_config$method, "file")) {
      checkmate::assert_file(dist_config$filename)
      if (isTRUE(return_cigar)) {
        stop("CIGAR strings are not supported for external distance files.")
      }
      if (span != 0L) {
        stop("Span is not supported for external distance files.")
      }
      if (dist_config$by_name == TRUE) {
        out <- utils::read.table(
          file = dist_config$filename,
          header = FALSE,
          col.names = c("seq_id", "ref_id", "dist"),
          colClasses = c("character", "character", "numeric")
        )
        # An external distance matrix may not have queries and references
        # distinguished, so search for pairs in both directions
        fwd_matches <- out$seq_id %in% names(query) & out$ref_id %in% names(ref)
        swap_matches <- out$seq_id %in%
          names(ref) &
          out$ref_id %in% names(query)
        data.frame(
          seq_id = c(out$seq_id[fwd_matches], out$ref_id[swap_matches]),
          ref_id = c(out$ref_id[fwd_matches], out$seq_id[swap_matches]),
          dist = c(out$dist[fwd_matches], out$dist[swap_matches]),
          stringsAsFactors = FALSE
        )
      } else {
        utils::read.table(
          file = dist_config$filename,
          header = FALSE,
          col.names = c("seq_id", "ref_id", "dist"),
          colClasses = c("integer", "integer", "numeric")
        )
        out$seq_id <- names(query)[out$seq_id + 1]
        out$ref_id <- names(ref)[out$ref_id + 1]
        out
      }
    } else {
      seq_search_internal(
        query = query,
        ref = ref,
        dist_config = dist_config,
        parallel_config = parallel_config,
        threshold = threshold,
        verbose = verbose,
        return_cigar = return_cigar,
        span = span
      )
    }
  }
}
