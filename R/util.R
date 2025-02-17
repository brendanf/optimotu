# SPDX-CopyrightText: (c) 2025, Brendan Furneaux
# SPDX-License-Identifier: MIT

# These functions are modified from the `optimotu.pipeline` package, where
# they are exported. Here they are only for internal use (and thus cannot use
# S3 dispatch)

#### generic sequence helpers ####
# helper functions to work with sequences sets that may be XStringSet,
# (named) character, fasta/fastq file, or data.frame

fasta_regex <- "\\.(fas?|fasta)(\\.gz)?$"
fastq_regex <- "\\.f(ast)?q(\\.gz)?$"

#' Guess the column name in a data frame which refers to the sequence ID
#' @param d (`data.frame`) data frame to search
#' @return (`character`) name of the column
find_name_col <- function(d) {
  stopifnot(is.data.frame(d))
  if ("seq_id" %in% names(d)) return("seq_id")
  if ("name" %in% names(d)) return("name")
  if ("ASV" %in% names(d)) return("ASV")
  if ("OTU" %in% names(d)) return("OTU")
  if (ncol(d) == 2) {
    if ("seq" %in% names(d)) return(names(d)[names(d) != "seq"])
    if ("sequence" %in% names(d)) return(names(d)[names(d) != "sequence"])
    return(names(d)[1])
  }
  stop("unable to determine sequence name column:", names(d))
}

#' Guess the column name in a data frame which refers to the sequence itself
#' @param d (`data.frame`) data frame to search
#' @return (`character`) name of the column
find_seq_col <- function(d) {
  stopifnot(is.data.frame(d))
  if ("seq" %in% names(d)) return("seq")
  if ("sequence" %in% names(d)) return("sequence")
  if (ncol(d) == 2) {
    if ("seq_id" %in% names(d)) return(names(d)[names(d) != "seq_id"])
    if ("name" %in% names(d)) return(names(d)[names(d) != "name"])
    if ("ASV" %in% names(d)) return(names(d)[names(d) != "ASV"])
    if ("OTU" %in% names(d)) return(names(d)[names(d) != "OTU"])
    return(names(d)[2])
  }
  stop("unable to determine sequence column:", names(d))
}

#' Write sequences as a fasta file
#'
#' @param seq (`data.frame`, `character`, or `XStringSet`) sequences to write
#' @param fname (`character` string) file path to write to
#' @param ... additional arguments to pass to the writing function
write_sequence <- function(seq, fname, ...) {
  args <- list(...)
  if (is.data.frame(seq)) {
    name_col <- args$name_col
    args$name_col <- NULL
    if (is.null(name_col)) {
      name_col <- find_name_col(seq)
    }
    seq_col <- args$seq_col
    args$seq_col <- NULL
    if (is.null(seq_col)) {
      seq_col <- find_seq_col(seq)
    }
    seq <- `names<-`(Biostrings::DNAStringSet(seq[[seq_col]]), seq[[name_col]])
  } else if (is.character(seq) && length(seq) == 1 && file.exists(seq)) {
    if (grepl(fasta_regex, seq)) {
      seq <- Biostrings::readDNAStringSet(seq)
    } else if (grepl(fastq_regex, seq)) {
      seq <- Biostrings::readDNAStringSet(seq, format = "fastq")
    } else {
      stop("unknown file type: ", seq)
    }
  } else if (is.character(seq)) {
    seq <- Biostrings::DNAStringSet(seq)
  }
  do.call(Biostrings::writeXStringSet, c(list(seq, fname), args))
  fname
}

#' Choose a subset of sequences
#' @param sequence (`data.frame`, `character`, or `XStringSet`) sequences to
#' subset
#' @param which (`integer`, `logical`, or `character`) which sequences to
#' select
#' @param negate (`logical`) if TRUE, select all sequences except those
#' specified
#' @param name_col (`character`) name of the column in a data frame which
#' contains the sequence names
select_sequence <- function(sequence, which, negate = FALSE, name_col = NULL) {
  checkmate::assert_flag(negate)
  checkmate::assert(
    checkmate::test_integerish(which),
    checkmate::test_logical(which),
    checkmate::test_character(which)
  )
  checkmate::assert_string(name_col, null.ok = TRUE)

  if (isTRUE(negate)) {
    if (is.numeric(which)){
      which <- -which
    } else if (is.logical(which)) {
      which <- !which
    } else if (is.character(which)) {
      which <- setdiff(seq_names(sequence, name_col = name_col), which)
    }
  }

  if (is.character(sequence) && length(sequence) == 1 && file.exists(sequence)) {
    if (grepl(fasta_regex, sequence)) {
      seq <- Biostrings::readBStringSet(sequence)
    } else if (grepl(fastq_regex, sequence)) {
      seq <- Biostrings::readBStringSet(sequence, format = "fastq")
    } else {
      stop("unknown file type: ", sequence)
    }
    seq <- as.character(seq)
  } else if (is.data.frame(sequence)) {
    if (is.null(name_col)) {
      name_col <- find_name_col(sequence)
    }
    if (is.numeric(which)) return(sequence[which,])
    if (is.logical(which)) return(sequence[which,])
    if (is.character(which)) {
      out <- data.frame(which)
      names(out) <- name_col
      return(merge(out, sequence, by = name_col, sort = FALSE))
    }
  }
  seq[which]
}

#' Get the number of sequences in a set
#' @param seq (`data.frame`, `character`, or `XStringSet`) sequences to count
sequence_size <- function(seq) {
  if (is.character(seq) && length(seq) == 1 && file.exists(seq)) {
    if (grepl(fasta_regex, seq)) {
      return(length(Biostrings::fasta.seqlengths(seq)))
    } else if (grepl(fastq_regex, seq)) {
      return(length(Biostrings::fastq.seqlengths(seq))
      )
    } else {
      stop("unknown file type: ", seq)
  }
  } else if (is.data.frame(seq)) {
    return(nrow(seq))
  } else {
    return(length(seq))
  }
}

#' get the names of a set of sequences
#' @param seq (`data.frame`, `character`, or `XStringSet`) sequences to get
#' the names of
#' @param name_col (`character`) name of the column in a data frame which
#' contains the sequence names
#' @return (`character` vector) the names of the sequences
seq_names <- function(seq, name_col = NULL) {
  if (is.character(seq) && length(seq) == 1 && file.exists(seq)) {
    if (grepl(fasta_regex, seq)) {
      names(Biostrings::fasta.seqlengths(seq))
    } else if (grepl(fastq_regex, seq)) {
      fastq_names(seq)
    } else {
      stop("unknown file type: ", seq)
    }
  } else if (is.data.frame(seq)) {
    if (is.null(name_col)) {
      name_col <- find_name_col(seq)
    }
    seq[[name_col]]
  } else {
    names(seq)
  }
}

#' set the names of a set of sequences
#' @param seq (`data.frame`, `character`, or `XStringSet`) sequences to set
#' the names of
#' @param value (`character` vector) the names to set
#' @return the sequences with the names set (invisibly)
`seq_names<-` <- function(seq, value) {
  if (is.character(seq) && length(seq) == 1 && file.exists(seq)) {
    stop("cannot set sequence names on a file")
  } else if (is.data.frame(seq)) {
    name_col <- tryCatch(
      find_name_col(seq),
      error = function(e) "seq_id"
    )
    seq[[name_col]] <- value
  } else {
    names(seq) <- value
  }
  invisible(seq)
}

#' get sequences as a character vector
#' @param seq (`data.frame`, `character`, or `XStringSet`) sequences to get
#' @param seq_col (`character`) name of the column in a data frame which
#' contains the sequences
#' @param name_col (`character`) name of the column in a data frame which
#' contains the sequence names
seq_as_char <- function(seq, seq_col = NULL, name_col = NULL) {
  if (is.character(seq) && length(seq) == 1 && file.exists(seq)) {
    if (grepl(fasta_regex, seq)) {
      Biostrings::readBStringSet(seq)
    } else if (grepl(fastq_regex, seq)) {
      Biostrings::readBStringSet(seq, format = "fastq")
    } else {
      stop("unknown file type: ", seq)
    }
    seq <- as.character(seq)
  } else if (is.data.frame(seq)) {
    if (is.null(seq_col)) {
      seq_col <- find_seq_col(seq)
    }
    if (is.null(name_col)) {
      name_col <- find_name_col(seq)
    }
    seq <- `names<-`(as.character(seq[[seq_col]]), seq[[name_col]])
  } else if (!is.character(seq)) {
    seq <- as.character(seq)
  }
  seq
}

#' Generate names for a set of sequences
#'
#' Names produced are like "ASV0001", "ASV0002", ... The prefix ("ASV" in
#' the example) is customizable, and the number will be zero-padded to fit
#' all values.
#'
#' @param n (`integer`) number of sequences to name
#' @param prefix (`character`) prefix to use for the names
#' @return (`character` vector) names
make_seq_names <- function(n, prefix) {
  sprintf(
    sprintf("%s%%0%dd", prefix, max(floor(log10(n)) + 1L, 0L)),
    seq_len(n)
  )
}

#### helpers for external commands ####

#' Find an executable
#'
#' This function tries to find an executable in the system. It first checks the
#' environment variables, then the system path, and finally the current working
#' directory.
#'
#' @param executable (`character` string) the name of the executable to find.
#' @return (`character` string) the full path to the executable.

find_executable <- function(executable) {
  checkmate::assert_character(executable)
  out <- Sys.getenv(executable)
  if (nchar(out) == 0 || !file.exists(out)) {
    out <- Sys.getenv(toupper(executable))
  }
  if (nchar(out) == 0 || !file.exists(out)) {
    out <- Sys.which(executable)
  }
  if (nchar(out) == 0 || !file.exists(out)) {
    out <- list.files(path = "bin", pattern = executable, recursive = TRUE,
                      full.names = TRUE)
  }
  checkmate::assert_file_exists(out, access = "x", .var.name = executable)
  out
}

#' Try to find the usearch executable
#' @return (`character` string) the full path to the usearch executable.
find_usearch <- function() {
  find_executable("usearch")
}
