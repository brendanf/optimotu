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
#' @param ... passed to methods
#' @return (`data.frame`) with columns "seq_id" (`character`),
#' "ref_id" (`character`), and "dist" (`numeric`) giving the distance between
#' the query and reference.
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
    ...
) {
  checkmate::assert_character(query_id, null.ok = TRUE, len = length(query), unique = TRUE)
  checkmate::assert_character(ref_id, null.ok = TRUE, len = length(ref), unique = TRUE)
  checkmate::assert_class(dist_config, "optimotu_dist_config")
  checkmate::assert_class(parallel_config, "optimotu_parallel_config")

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
    seq_search_internal(
      query = query,
      ref = ref,
      dist_config = dist_config,
      parallel_config = parallel_config,
      threshold = threshold,
      verbose = verbose
    )
  }
}
