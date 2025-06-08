#' Search for the closest match(es) to sequences in a reference database
#'
#' @param query (`character`, `data.frame`, or `XStringSet`) sequences to search
#' @param query_id (`character` vector) names for the query sequences.  If they
#' are already named, this will replace the names.
#' @param ref (`character`, `data.frame`, or `XStringSet`) reference sequences to
#' search against
#' @param ref_id (`character` vector) names for the reference sequences.  If they
#' are already named, this will replace the names.
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
    ...
) {
  checkmate::assert_character(query_id, null.ok = TRUE, len = length(query), unique = TRUE)
  checkmate::assert_character(ref_id, null.ok = TRUE, len = length(ref), unique = TRUE)
  checkmate::assert_class(dist_config, "optimotu_dist_config")
  checkmate::assert_class(parallel_config, "optimotu_parallel_config")
  checkmate::assert_flag(return_cigar)
  checkmate::assert_character(span)
  span <- match.arg(span, c("global", "extension"), several.ok = FALSE)
  checkmate::assert_choice(span, c("global", "extension"))
  span <- match(span, c("global", "extension")) - 1L

  mycall <- match.call()

  if (identical(dist_config$method, "usearch")) {
    mycall[[1]] <- seq_search_usearch
    mycall$dist_config <- NULL
    mycall$usearch <- dist_config$usearch
    eval(mycall, envir = parent.frame())
  } else {
    query <- seq_as_char(query)
    if (!is.null(query_id)) {
      names(query) <- query_id
    }
    ref <- seq_as_char(ref)
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
        swap_matches <- out$seq_id %in% names(ref) & out$ref_id %in% names(query)
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
        out$query_id <- names(ref)[out$ref_id + 1]
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
