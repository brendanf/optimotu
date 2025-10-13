#' Generates a sparse distance matrix using USEARCH
#'
#' @description This function uses a unix pipe to read the output of the
#' USEARCH "`allpairs_global`" command. USEARCH (version 8.0 or later) should
#' be installed separately; it is available with a free license for most users
#' at https://www.drive5.com/usearch/.
#'
#' @inheritParams seq_distmx
#' @param seq_id (`character` vector) names for the query sequences.  If they
#' are already named, this will replace the names.  Has no effect if `query` is
#' a filename.
#' @param usearch (`character` scalar) path to USEARCH executable
#' @param ... passed to methods (currently inactive)
#'
#' @return (`data.frame`) with columns "seq_id" (`character`),
#' "ref_id" (`character`), and "dist" (`numeric`) giving the distance between
#' the query and reference.
#'
#' @export
seq_distmx_usearch <- function(
  seq,
  threshold,
  seq_id = NULL,
  parallel_config = parallel_concurrent(1),
  usearch = find_usearch(),
  verbose = FALSE,
  details = c("none", "gapstats", "cigar"),
  id_is_int = FALSE,
  ...
) {
  if (!missing(details)) checkmate::assert_string(details)
  details <- match.arg(details, several.ok = FALSE)
  checkmate::assert_flag(id_is_int)

  seq_id1 <- "seq_id1"
  seq_id2 <- "seq_id2"
  id_type <- "character"
  if (id_is_int) {
    seq_id1 <- "seq_idx1"
    seq_id2 <- "seq_idx2"
    id_type <- "integer"

  }
  colnames <- if (details == "none") {
    c(seq_id1, seq_id2, "score2", "dist2")
  } else {
    c(seq_id1, seq_id2, "score2", "dist2", "cigar")
  }
  coltypes <- if (details == "none") {
    c(id_type, id_type, "numeric", "numeric")
  } else {
    c(id_type, id_type, "numeric", "numeric", "character")
  }
  userfields <- switch(details,
    "none" = "query+target+raw+id",
    "gapstats" = "query+target+raw+id+caln",
    "cigar" = "query+target+raw+id+caln"
  )
  if (is.character(seq) && length(seq) == 1 && file.exists(seq)) {
    tseq <- seq
  } else {
    if (!is.null(seq_id)) {
      checkmate::assert_character(seq_id, len = length(seq), unique = TRUE)
      seq_names(seq) <- seq_id
    }
    tseq <- tempfile(pattern = "seq", fileext = ".fasta")
    write_sequence(seq, tseq)
    on.exit(unlink(tseq))
  }
  tout <- tempfile(fileext = ".dat")
  file.create(tout)
  on.exit(unlink(tout), add = TRUE)
  args <- c(
    "-allpairs_global", tseq,
    "-id", 1 - threshold,
    "-termdist", min(1, 2 * threshold),
    "-lopen", 0,
    "-lext", 1,
    "-threads", parallel_config$threads,
    "-userout", tout,
    "-userfields", userfields
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
  if (details == "gapstats") {
    out <- add_gapstats(out, "cigar")
    out$cigar <- NULL
  }
  out$dist2 <- 1 - out$dist2 / 100
  out$score1 <- out$score2
  out$dist1 <- out$dist2
  if (details == "none") {
    out <- out[c(seq_id1, seq_id2, "score1", "dist1", "score2", "dist2")]
  } else if (details == "gapstats") {
    out <- out[c(seq_id1, seq_id2, "score1", "dist1", "score2", "dist2",
                 "align_length", "n_insert", "n_delete", "max_insert",
                 "max_delete")]
  } else if (details == "cigar") {
    out <- out[c(seq_id1, seq_id2, "score1", "dist1", "score2", "dist2",
                 "cigar")]
  }
  attr(out, "class") <- c("tbl_df", "tbl", "data.frame")
  out
}
