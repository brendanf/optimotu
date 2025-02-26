#' Search for the best match(es) to query sequences in a set of reference
#' sequences using USEARCH
#'
#' @description This function uses a unix pipe to read the output of the
#' USEARCH "`usearch_global`" command. USEARCH (version 8.0' or higher) should
#' be installed separately; it is available with a free license for most users
#' at https://www.drive5.com/usearch/.
#'
#' @inheritParams seq_search
#' @param query_id (`character` vector) names for the query sequences.  If they
#' are already named, this will replace the names.  Has no effect if `query` is
#' a filename.
#' @param ref_id (`character` vector) names for the reference sequences.  If
#' they are already named, this will replace the names.  Has no effect if `seq`
#' is a filename.
#' @param usearch (`character` scalar) path to USEARCH executable
#' @param ... passed to methods (currently inactive)
#'
#' @return (`data.frame`) with columns "seq_id" (`character`),
#' "ref_id" (`character`), and "dist" (`numeric`) giving the distance between
#' the query and reference.
#'
#' @export
seq_search_usearch <- function(
    query,
    ref,
    threshold,
    query_id = NULL,
    ref_id = NULL,
    parallel_config = parallel_concurrent(1),
    usearch = find_usearch(),
    verbose = FALSE,
    ...
) {
  if (is.character(query) && length(query) == 1 && file.exists(query)) {
    tquery <- query
  } else {
    if (!is.null(query_id)) {
      checkmate::assert_character(query_id, len = length(query), unique = TRUE)
      seq_names(query) <- query_id
    }
    tquery <-tempfile(pattern = "query", fileext = ".fasta")
    write_sequence(query, tquery)
    on.exit(unlink(tquery))
  }
  if (is.character(ref) && length(ref) == 1 && file.exists(ref)) {
    tref <- ref
  } else {
    if (!is.null(ref_id)) {
      checkmate::assert_character(ref_id, len = length(ref), unique = TRUE)
      seq_names(ref) <- ref_id
    }
    tref <- tempfile(pattern = "ref", fileext = ".fasta")
    write_sequence(ref, tref)
    on.exit(unlink(tref), add = TRUE)
  }
  tout <- tempfile(fileext = ".dat")
  file.create(tout)
  on.exit(unlink(tout), add = TRUE)
  args <- c(
    "-usearch_global", tquery,
    "-db", tref,
    "-id", 1 - threshold,
    "-strand", "plus",
    "-lopen", 0,
    "-lext", 1,
    "-maxaccepts", 100,
    "-top_hits_only",
    "-threads", parallel_config$threads,
    "-userout", tout,
    "-userfields", "query+target+id"
  )
  extra_args <- list(...)
  if (length(extra_args) > 0) {
    checkmate::assert_named(extra_args)
    extra_args <- c(rbind(paste0("-", names(extra_args)), unlist(extra_args)))
    args <- c(args, extra_args)
  }
  if (system2(usearch, stdout = NULL, stderr = NULL) != 0) {
    stop("usearch could not be found at path: ", usearch)
  }
  status <- system2(
    usearch,
    args,
    stdout = if (verbose) "" else NULL,
    stderr = if (verbose) "" else NULL
  )
  stopifnot(status == 0)
  out <- utils::read.table(
    tout,
    header = FALSE,
    col.names = c("seq_id", "ref_id", "dist"),
    colClasses = c("character", "character", "numeric")
  )
  out$dist <- 1 - out$dist/100
  attr(out, "class") <- c("tbl_df", "tbl", "data.frame")
  out
}
