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
#' @param return_cigar (`logical` scalar) if `TRUE`, return the CIGAR string
#' from USEARCH field `caln` in a column named `cigar`.
#' @param ... passed to methods (currently inactive)
#'
#' @return (`data.frame`) with columns "seq_id" (`character`),
#' "ref_id" (`character`), and "dist" (`numeric`) giving the distance between
#' the query and reference. If `return_cigar` is `TRUE`, also includes
#' "cigar" (`character`).
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
  return_cigar = FALSE,
  query_seq_idx = NULL,
  query_seq_file = NULL,
  ref_seq_idx = NULL,
  ref_seq_file = NULL,
  ...
) {
  checkmate::assert_flag(return_cigar)
  tq <- seq_usearch_fasta_file(query, query_id, query_seq_idx, query_seq_file)
  tquery <- tq$path
  if (tq$unlink) {
    on.exit(unlink(tquery), add = TRUE)
  }
  tr <- seq_usearch_fasta_file(ref, ref_id, ref_seq_idx, ref_seq_file)
  tref <- tr$path
  if (tr$unlink) {
    on.exit(unlink(tref), add = TRUE)
  }
  tout <- tempfile(fileext = ".dat")
  file.create(tout)
  on.exit(unlink(tout), add = TRUE)
  userfields <- if (isTRUE(return_cigar)) {
    "query+target+id+caln"
  } else {
    "query+target+id"
  }
  colnames <- if (isTRUE(return_cigar)) {
    c("seq_id", "ref_id", "dist", "cigar")
  } else {
    c("seq_id", "ref_id", "dist")
  }
  coltypes <- if (isTRUE(return_cigar)) {
    c("character", "character", "numeric", "character")
  } else {
    c("character", "character", "numeric")
  }
  args <- c(
    "-usearch_global",
    tquery,
    "-db",
    tref,
    "-id",
    1 - threshold,
    "-strand",
    "plus",
    "-lopen",
    0,
    "-lext",
    1,
    "-maxaccepts",
    100,
    "-top_hits_only",
    "-threads",
    parallel_config$threads,
    "-userout",
    tout,
    "-userfields",
    userfields
  )
  extra_args <- list(...)
  if (length(extra_args) > 0) {
    checkmate::assert_named(extra_args)
    extra_args <- c(rbind(paste0("-", names(extra_args)), unlist(extra_args)))
    args <- c(args, extra_args)
  }
  if (system2(usearch, "-version", stdout = NULL, stderr = NULL) != 0) {
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
    col.names = colnames,
    colClasses = coltypes,
  )
  out$dist <- 1 - out$dist / 100
  attr(out, "class") <- c("tbl_df", "tbl", "data.frame")
  out
}
