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
  if ("seq_id" %in% names(d)) {
    return("seq_id")
  }
  if ("name" %in% names(d)) {
    return("name")
  }
  if ("ASV" %in% names(d)) {
    return("ASV")
  }
  if ("OTU" %in% names(d)) {
    return("OTU")
  }
  if (ncol(d) == 2) {
    if ("seq" %in% names(d)) {
      return(names(d)[names(d) != "seq"])
    }
    if ("sequence" %in% names(d)) {
      return(names(d)[names(d) != "sequence"])
    }
    return(names(d)[1])
  }
  stop("unable to determine sequence name column:", names(d))
}

#' Guess the column name in a data frame which refers to the sequence itself
#' @param d (`data.frame`) data frame to search
#' @return (`character`) name of the column
find_seq_col <- function(d) {
  stopifnot(is.data.frame(d))
  if ("seq" %in% names(d)) {
    return("seq")
  }
  if ("sequence" %in% names(d)) {
    return("sequence")
  }
  if (ncol(d) == 2) {
    if ("seq_id" %in% names(d)) {
      return(names(d)[names(d) != "seq_id"])
    }
    if ("name" %in% names(d)) {
      return(names(d)[names(d) != "name"])
    }
    if ("ASV" %in% names(d)) {
      return(names(d)[names(d) != "ASV"])
    }
    if ("OTU" %in% names(d)) {
      return(names(d)[names(d) != "OTU"])
    }
    return(names(d)[2])
  }
  stop("unable to determine sequence column:", names(d))
}

#' Whether `x` is a character vector of existing `.fqi` index paths
#'
#' @param x object to test
#' @return `TRUE` if `x` is a non-empty character vector of existing paths whose
#'   names end in `.fqi` (case-insensitive).
#' @noRd
seq_is_fqi_path_set <- function(x) {
  if (!is.character(x) || anyNA(x) || length(x) < 1L) {
    return(FALSE)
  }
  if (!all(file.exists(x))) {
    return(FALSE)
  }
  all(grepl("\\.fqi$", x, ignore.case = TRUE))
}

#' Whether sequence input requires random-access extraction via fastqindexr
#'
#' @param seq sequence argument
#' @return `TRUE` if `seq` is a `fastqindexr_index` object or a `.fqi` path set.
#' @noRd
seq_is_index_backed_input <- function(seq) {
  inherits(seq, "fastqindexr_index") || seq_is_fqi_path_set(seq)
}

#' Require suggested package fastqindexr (>= 0.1)
#'
#' @param what short phrase for the error message
#' @noRd
optimotu_require_fastqindexr <- function(what) {
  if (
    !requireNamespace("fastqindexr", quietly = TRUE) ||
      utils::packageVersion("fastqindexr") < "0.1"
  ) {
    stop(
      "Package 'fastqindexr' (>= 0.1) is required ",
      what,
      ".",
      call. = FALSE
    )
  }
}

#' Resolve `seq_idx` against a fixed record count
#'
#' @param n_records total number of sequences (`integer` scalar)
#' @param seq_idx `NULL` or integerish vector of positions in `1:n_records`
#' @return integer vector of selected indices (possibly empty)
#' @noRd
seq_resolve_linear_seq_idx <- function(n_records, seq_idx) {
  checkmate::assert_count(n_records)
  if (is.null(seq_idx)) {
    if (n_records < 1L) {
      return(integer())
    }
    return(seq_len(n_records))
  }
  checkmate::assert_integerish(
    seq_idx,
    lower = 1L,
    upper = n_records,
    any.missing = FALSE
  )
  as.integer(seq_idx)
}

#' Number of sequences represented by an input (before optional `seq_idx`)
#'
#' Used for argument checks when `seq` is a file path or index object.
#'
#' @param seq sequence argument (same kinds as [seq_as_char()]).
#' @param seq_file optional path override for index-backed inputs.
#' @return non-negative integer count
#' @noRd
seq_input_n_records <- function(seq, seq_file = NULL) {
  if (seq_is_index_backed_input(seq)) {
    optimotu_require_fastqindexr("to read index-backed sequence inputs")
    idx <- seq_load_index_object(seq, seq_file)
    return(as.integer(idx$n_records))
  }
  if (
    is.character(seq) && length(seq) == 1L && nzchar(seq) && file.exists(seq)
  ) {
    if (grepl(fasta_regex, seq)) {
      return(as.integer(length(Biostrings::fasta.seqlengths(seq))))
    }
    if (grepl(fastq_regex, seq)) {
      return(as.integer(length(Biostrings::fastq.seqlengths(seq))))
    }
    stop("unknown file type: ", seq, call. = FALSE)
  }
  if (is.data.frame(seq)) {
    return(nrow(seq))
  }
  if (methods::is(seq, "XStringSet")) {
    return(length(seq))
  }
  if (is.character(seq)) {
    return(length(seq))
  }
  stop("Unsupported sequence input type: ", paste(class(seq), collapse = "/"))
}

#' Load a `fastqindexr_index` from an index object or `.fqi` path vector
#'
#' @noRd
seq_load_index_object <- function(seq, seq_file = NULL) {
  if (inherits(seq, "fastqindexr_index")) {
    return(seq)
  }
  if (seq_is_fqi_path_set(seq)) {
    optimotu_require_fastqindexr("to read .fqi index files")
    return(fastqindexr::read_fqi_index(
      fqi_path = seq,
      files = seq_file,
      type = "auto"
    ))
  }
  stop("Internal error: not an index-backed input", call. = FALSE)
}

#' Sequence record names for a subset (without loading full sequence letters
#' when possible for plain FASTA paths)
#'
#' @inheritParams seq_as_char
#' @param name_col optional name column for `data.frame` inputs
#' @return character vector of sequence identifiers in selection order
#' @noRd
seq_input_seq_ids <- function(
  seq,
  seq_idx = NULL,
  seq_file = NULL,
  name_col = NULL
) {
  if (!is.null(seq_idx) && is.data.frame(seq)) {
    ii <- seq_resolve_linear_seq_idx(nrow(seq), seq_idx)
    seq <- if (length(ii) < 1L) {
      seq[integer(), , drop = FALSE]
    } else {
      seq[ii, , drop = FALSE]
    }
    seq_idx <- NULL
  } else if (!is.null(seq_idx) && methods::is(seq, "XStringSet")) {
    ii <- seq_resolve_linear_seq_idx(length(seq), seq_idx)
    seq <- if (length(ii) < 1L) seq[integer()] else seq[ii]
    seq_idx <- NULL
  } else if (
    !is.null(seq_idx) &&
      is.character(seq) &&
      (length(seq) != 1L || !file.exists(seq))
  ) {
    ii <- seq_resolve_linear_seq_idx(length(seq), seq_idx)
    seq <- if (length(ii) < 1L) character() else seq[ii]
    seq_idx <- NULL
  }

  if (seq_is_index_backed_input(seq)) {
    optimotu_require_fastqindexr(
      "to read sequence IDs from index-backed inputs"
    )
    idx <- seq_load_index_object(seq, seq_file)
    nrec <- as.integer(idx$n_records)
    ii <- seq_resolve_linear_seq_idx(nrec, seq_idx)
    if (length(ii) < 1L) {
      return(character())
    }
    if (
      is.null(seq_file) &&
        length(idx$files) == 1L &&
        identical(idx$format, "fasta") &&
        length(ii) == nrec &&
        identical(ii, seq_len(nrec))
    ) {
      p <- as.character(idx$files[1L])
      return(names(Biostrings::fasta.seqlengths(p)))
    }
    df <- fastqindexr::extract_sequences(
      index = idx,
      seq_idx = ii,
      file = seq_file,
      return = "data.frame",
      renumber = "none"
    )
    return(as.character(df$seq_id))
  }

  if (is.character(seq) && length(seq) == 1L && file.exists(seq)) {
    if (grepl(fasta_regex, seq)) {
      nm <- names(Biostrings::fasta.seqlengths(seq))
      ii <- seq_resolve_linear_seq_idx(length(nm), seq_idx)
      return(nm[ii])
    }
    if (grepl(fastq_regex, seq)) {
      nm <- names(Biostrings::fastq.seqlengths(seq))
      ii <- seq_resolve_linear_seq_idx(length(nm), seq_idx)
      return(nm[ii])
    }
    stop("unknown file type: ", seq, call. = FALSE)
  }

  seq_names(seq, name_col = name_col)
}

fastqindexr_ok <- function() {
  requireNamespace("fastqindexr", quietly = TRUE) &&
    utils::packageVersion("fastqindexr") >= "0.1"
}

#' Subset a FASTA file using Biostrings index (no full-file read)
#'
#' @noRd
seq_read_fasta_subset_to_dna <- function(path, seq_idx) {
  fi <- Biostrings::fasta.index(path)
  n <- nrow(fi)
  ii <- seq_resolve_linear_seq_idx(n, seq_idx)
  if (length(ii) < 1L) {
    return(Biostrings::DNAStringSet())
  }
  if (fastqindexr_ok()) {
    idx <- fastqindexr::create_index(path, type = "fasta")
    ex <- fastqindexr::extract_sequences(
      index = idx,
      seq_idx = ii,
      file = NULL,
      return = "seq",
      renumber = "none"
    )
    return(Biostrings::DNAStringSet(ex))
  }
  nm <- names(Biostrings::fasta.seqlengths(path))
  full <- Biostrings::readBStringSet(path)
  methods::as(full[nm[ii]], "DNAStringSet")
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

  if (seq_is_index_backed_input(sequence)) {
    optimotu_require_fastqindexr("for index-backed inputs to select_sequence()")
    sequence <- seq_as_char(sequence, as = "data.frame")
  }

  if (isTRUE(negate)) {
    if (is.numeric(which)) {
      which <- -which
    } else if (is.logical(which)) {
      which <- !which
    } else if (is.character(which)) {
      which <- setdiff(seq_names(sequence, name_col = name_col), which)
    }
  }

  if (
    is.character(sequence) && length(sequence) == 1 && file.exists(sequence)
  ) {
    if (grepl(fasta_regex, sequence)) {
      sequence <- Biostrings::readBStringSet(sequence)
    } else if (grepl(fastq_regex, sequence)) {
      sequence <- Biostrings::readBStringSet(sequence, format = "fastq")
    } else {
      stop("unknown file type: ", sequence)
    }
    sequence <- as.character(sequence)
  } else if (is.data.frame(sequence)) {
    if (is.null(name_col)) {
      name_col <- find_name_col(sequence)
    }
    if (is.numeric(which)) {
      return(sequence[which, ])
    }
    if (is.logical(which)) {
      return(sequence[which, ])
    }
    if (is.character(which)) {
      out <- data.frame(which)
      names(out) <- name_col
      return(merge(out, sequence, by = name_col, sort = FALSE))
    }
  }
  sequence[which]
}

#' Get the number of sequences in a set
#' @param seq (`data.frame`, `character`, or `XStringSet`) sequences to count
sequence_size <- function(seq) {
  if (seq_is_index_backed_input(seq)) {
    return(seq_input_n_records(seq))
  }
  if (is.character(seq) && length(seq) == 1 && file.exists(seq)) {
    if (grepl(fasta_regex, seq)) {
      return(length(Biostrings::fasta.seqlengths(seq)))
    } else if (grepl(fastq_regex, seq)) {
      return(length(Biostrings::fastq.seqlengths(seq)))
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
  if (seq_is_index_backed_input(seq)) {
    return(seq_input_seq_ids(seq, seq_idx = NULL, seq_file = NULL, name_col))
  }
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

#' Harmonize sequence inputs to character, paths, or Biostrings objects
#'
#' Supports named `character` vectors (literal sequences), `data.frame` and
#' [`XStringSet`][Biostrings::XStringSet-class] objects, single FASTA/FASTQ file
#' paths, and (when suggested package `fastqindexr` is installed) a
#' `fastqindexr_index` object or a character vector of existing `.fqi` paths.
#'
#' Optional `seq_idx` selects 1-based records (for files and index-backed
#' inputs) or rows / elements (for `data.frame`, `XStringSet`, and multi-element
#' `character` literals). `NULL` means all records in order.
#'
#' `seq_file` may be supplied only for `fastqindexr_index` or `.fqi` path inputs;
#' it overrides the indexed FASTA/FASTQ paths passed to `fastqindexr`.
#'
#' For index-backed inputs or `.fqi` files, `fastqindexr` (>= 0.1) is required.
#' For plain FASTA files with `seq_idx` not `NULL`, random access uses
#' `fastqindexr` when installed; otherwise the file is read with
#' [Biostrings::readBStringSet()] and subset in memory (avoid for very large
#' inputs). Plain FASTQ subsetting without `fastqindexr` reads the entire file,
#' then subsets (avoid for very large inputs).
#'
#' @param seq (`character`, `data.frame`, `XStringSet`, single FASTA/FASTQ path,
#'   `fastqindexr_index`, or character vector of `.fqi` paths).
#' @param seq_idx optional integerish vector of 1-based indices (`NULL` = all).
#' @param seq_file optional character vector of path overrides for index-backed
#'   inputs only.
#' @param seq_col,name_col optional column names for `data.frame` inputs.
#' @param as desired output: `"character"` (named character vector),
#'   `"dna_string_set"`, `"fasta_path"` (single FASTA path; attribute
#'   `optimotu_unlink` is `TRUE` when the caller should [unlink()] the file),
#'   or `"data.frame"` with columns `seq_id` and `seq`.
#' @return Depends on `as`; see parameter description.
#' @export
seq_as_char <- function(
  seq,
  seq_idx = NULL,
  seq_file = NULL,
  seq_col = NULL,
  name_col = NULL,
  as = c("character", "dna_string_set", "fasta_path", "data.frame")
) {
  as <- match.arg(as)
  if (!is.null(seq_file) && !seq_is_index_backed_input(seq)) {
    stop(
      "`seq_file` is only allowed when `seq` is a fastqindexr_index ",
      "or a character vector of existing .fqi files.",
      call. = FALSE
    )
  }

  # Fast path: preserve legacy behavior for common inputs (no subsetting)
  if (is.null(seq_idx) && is.null(seq_file) && as == "character") {
    if (is.character(seq)) {
      if (length(seq) != 1L || !file.exists(seq)) {
        return(seq)
      }
      if (grepl(fasta_regex, seq)) {
        bs <- Biostrings::readBStringSet(seq)
        return(`names<-`(as.character(bs), names(bs)))
      }
      if (grepl(fastq_regex, seq)) {
        bs <- Biostrings::readBStringSet(seq, format = "fastq")
        return(`names<-`(as.character(bs), names(bs)))
      }
      stop("unknown file type: ", seq, call. = FALSE)
    }
    if (methods::is(seq, "XStringSet")) {
      return(as.character(seq))
    }
  }

  # --- in-memory subsetting (seq_idx as row/element indices) ---
  if (!is.null(seq_idx) && is.data.frame(seq)) {
    ii <- seq_resolve_linear_seq_idx(nrow(seq), seq_idx)
    if (length(ii) < 1L) {
      seq <- seq[integer(), , drop = FALSE]
    } else {
      seq <- seq[ii, , drop = FALSE]
    }
    seq_idx <- NULL
  } else if (!is.null(seq_idx) && methods::is(seq, "XStringSet")) {
    ii <- seq_resolve_linear_seq_idx(length(seq), seq_idx)
    seq <- if (length(ii) < 1L) {
      seq[integer()]
    } else {
      seq[ii]
    }
    seq_idx <- NULL
  } else if (
    !is.null(seq_idx) &&
      is.character(seq) &&
      (length(seq) != 1L || !file.exists(seq))
  ) {
    ii <- seq_resolve_linear_seq_idx(length(seq), seq_idx)
    seq <- if (length(ii) < 1L) {
      character()
    } else {
      seq[ii]
    }
    seq_idx <- NULL
  }

  # --- index-backed (fastqindexr) ---
  if (seq_is_index_backed_input(seq)) {
    optimotu_require_fastqindexr("for index-backed sequence inputs")
    idx <- seq_load_index_object(seq, seq_file)
    nrec <- as.integer(idx$n_records)
    ii <- seq_resolve_linear_seq_idx(nrec, seq_idx)
    if (as == "fasta_path") {
      if (length(ii) < 1L) {
        tf <- tempfile(pattern = "seq", fileext = ".fasta")
        file.create(tf)
        return(structure(tf, optimotu_unlink = TRUE))
      }
      if (
        is.null(seq_file) &&
          length(ii) == nrec &&
          identical(ii, seq_len(nrec)) &&
          length(idx$files) == 1L &&
          identical(idx$format, "fasta") &&
          grepl(fasta_regex, idx$files[1L]) &&
          !grepl("\\.gz$", idx$files[1L], ignore.case = TRUE)
      ) {
        fp <- as.character(idx$files[1L])
        return(structure(fp, optimotu_unlink = FALSE))
      }
      tf <- tempfile(pattern = "seq", fileext = ".fasta")
      fastqindexr::extract_sequences_to_file(
        index = idx,
        seq_idx = ii,
        file = seq_file,
        outfile = tf,
        type = "fasta",
        append = FALSE,
        compress = FALSE,
        collapse_sequence_lines = TRUE
      )
      return(structure(tf, optimotu_unlink = TRUE))
    }
    ex <- fastqindexr::extract_sequences(
      index = idx,
      seq_idx = ii,
      file = seq_file,
      return = "seq",
      renumber = "none"
    )
    if (as == "character") {
      return(ex)
    }
    if (as == "dna_string_set") {
      return(Biostrings::DNAStringSet(ex))
    }
    if (as == "data.frame") {
      return(data.frame(
        seq_id = names(ex),
        seq = as.character(ex),
        stringsAsFactors = FALSE
      ))
    }
  }

  # --- single file path (plain FASTA/FASTQ, not .fqi) ---
  if (is.character(seq) && length(seq) == 1L && file.exists(seq)) {
    if (grepl(fasta_regex, seq)) {
      n <- as.integer(length(Biostrings::fasta.seqlengths(seq)))
      ii <- seq_resolve_linear_seq_idx(n, seq_idx)
      if (as == "fasta_path") {
        if (length(ii) < 1L) {
          tf <- tempfile(pattern = "seq", fileext = ".fasta")
          file.create(tf)
          return(structure(tf, optimotu_unlink = TRUE))
        }
        if (is.null(seq_idx) && !grepl("\\.gz$", seq, ignore.case = TRUE)) {
          return(structure(seq, optimotu_unlink = FALSE))
        }
        dna <- seq_read_fasta_subset_to_dna(seq, ii)
        tf <- tempfile(pattern = "seq", fileext = ".fasta")
        Biostrings::writeXStringSet(dna, tf)
        return(structure(tf, optimotu_unlink = TRUE))
      }
      if (length(ii) < 1L) {
        out <- character()
        if (as == "dna_string_set") {
          return(Biostrings::DNAStringSet())
        }
        if (as == "data.frame") {
          return(data.frame(
            seq_id = character(),
            seq = character(),
            stringsAsFactors = FALSE
          ))
        }
        return(out)
      }
      dna <- seq_read_fasta_subset_to_dna(seq, ii)
      if (as == "character") {
        return(`names<-`(as.character(dna), names(dna)))
      }
      if (as == "dna_string_set") {
        return(dna)
      }
      if (as == "data.frame") {
        return(data.frame(
          seq_id = names(dna),
          seq = as.character(dna),
          stringsAsFactors = FALSE
        ))
      }
    }
    if (grepl(fastq_regex, seq)) {
      n <- as.integer(length(Biostrings::fastq.seqlengths(seq)))
      ii <- seq_resolve_linear_seq_idx(n, seq_idx)
      if (length(ii) < 1L) {
        if (as == "fasta_path") {
          tf <- tempfile(pattern = "seq", fileext = ".fasta")
          file.create(tf)
          return(structure(tf, optimotu_unlink = TRUE))
        }
        if (as == "dna_string_set") {
          return(Biostrings::DNAStringSet())
        }
        if (as == "data.frame") {
          return(data.frame(
            seq_id = character(),
            seq = character(),
            stringsAsFactors = FALSE
          ))
        }
        return(character())
      }
      if (fastqindexr_ok()) {
        idx <- fastqindexr::create_index(seq, type = "fastq")
        ex <- fastqindexr::extract_sequences(
          index = idx,
          seq_idx = ii,
          file = NULL,
          return = "seq",
          renumber = "none"
        )
        if (as == "character") {
          return(ex)
        }
        if (as == "dna_string_set") {
          return(Biostrings::DNAStringSet(ex))
        }
        if (as == "fasta_path") {
          tf <- tempfile(pattern = "seq", fileext = ".fasta")
          fastqindexr::extract_sequences_to_file(
            index = idx,
            seq_idx = ii,
            file = NULL,
            outfile = tf,
            type = "fasta",
            append = FALSE,
            compress = FALSE,
            collapse_sequence_lines = TRUE
          )
          return(structure(tf, optimotu_unlink = TRUE))
        }
        if (as == "data.frame") {
          df <- fastqindexr::extract_sequences(
            index = idx,
            seq_idx = ii,
            file = NULL,
            return = "data.frame",
            renumber = "none"
          )
          return(df[, c("seq_id", "seq"), drop = FALSE])
        }
      }
      # Fallback: load full FASTQ then subset (memory-heavy)
      bs <- Biostrings::readBStringSet(seq, format = "fastq")
      bs <- bs[ii]
      if (as == "character") {
        return(`names<-`(as.character(bs), names(bs)))
      }
      if (as == "dna_string_set") {
        return(methods::as(bs, "DNAStringSet"))
      }
      if (as == "fasta_path") {
        tf <- tempfile(pattern = "seq", fileext = ".fasta")
        Biostrings::writeXStringSet(methods::as(bs, "DNAStringSet"), tf)
        return(structure(tf, optimotu_unlink = TRUE))
      }
      if (as == "data.frame") {
        return(data.frame(
          seq_id = names(bs),
          seq = as.character(bs),
          stringsAsFactors = FALSE
        ))
      }
    }
    stop("unknown file type: ", seq, call. = FALSE)
  }

  # --- data.frame (no prior seq_idx or already row-subset) ---
  if (is.data.frame(seq)) {
    if (is.null(seq_col)) {
      seq_col <- find_seq_col(seq)
    }
    if (is.null(name_col)) {
      name_col <- find_name_col(seq)
    }
    ch <- `names<-`(as.character(seq[[seq_col]]), seq[[name_col]])
    if (as == "character") {
      return(ch)
    }
    if (as == "dna_string_set") {
      return(Biostrings::DNAStringSet(ch))
    }
    if (as == "fasta_path") {
      tf <- tempfile(pattern = "seq", fileext = ".fasta")
      Biostrings::writeXStringSet(Biostrings::DNAStringSet(ch), tf)
      return(structure(tf, optimotu_unlink = TRUE))
    }
    if (as == "data.frame") {
      return(data.frame(
        seq_id = names(ch),
        seq = as.character(ch),
        stringsAsFactors = FALSE
      ))
    }
  }

  # --- XStringSet ---
  if (methods::is(seq, "XStringSet")) {
    if (as == "character") {
      return(as.character(seq))
    }
    if (as == "dna_string_set") {
      return(methods::as(seq, "DNAStringSet"))
    }
    if (as == "fasta_path") {
      tf <- tempfile(pattern = "seq", fileext = ".fasta")
      Biostrings::writeXStringSet(methods::as(seq, "DNAStringSet"), tf)
      return(structure(tf, optimotu_unlink = TRUE))
    }
    if (as == "data.frame") {
      return(data.frame(
        seq_id = names(seq),
        seq = as.character(seq),
        stringsAsFactors = FALSE
      ))
    }
  }

  # --- named character literals ---
  if (is.character(seq)) {
    if (as == "character") {
      return(seq)
    }
    if (as == "dna_string_set") {
      return(Biostrings::DNAStringSet(seq))
    }
    if (as == "fasta_path") {
      tf <- tempfile(pattern = "seq", fileext = ".fasta")
      Biostrings::writeXStringSet(Biostrings::DNAStringSet(seq), tf)
      return(structure(tf, optimotu_unlink = TRUE))
    }
    if (as == "data.frame") {
      return(data.frame(
        seq_id = names(seq),
        seq = as.character(seq),
        stringsAsFactors = FALSE
      ))
    }
  }

  stop("Unsupported sequence input type: ", paste(class(seq), collapse = "/"))
}

#' Build a FASTA path for USEARCH, optionally rewriting headers from `seq_id`
#'
#' @return `list(path, unlink)` where `unlink` is `TRUE` if the caller should
#'   [unlink()] `path` after USEARCH finishes (register with [on.exit()] in the
#'   calling frame).
#' @noRd
seq_usearch_fasta_file <- function(
  seq,
  seq_id = NULL,
  seq_idx = NULL,
  seq_file = NULL
) {
  sp <- seq_as_char(
    seq,
    seq_idx = seq_idx,
    seq_file = seq_file,
    as = "fasta_path"
  )
  tseq <- unname(as.character(sp))
  unlink_prev <- isTRUE(attr(sp, "optimotu_unlink", exact = TRUE))
  if (!is.null(seq_id)) {
    n_seq <- length(Biostrings::fasta.seqlengths(tseq))
    checkmate::assert_character(seq_id, len = n_seq, unique = TRUE)
    dna <- Biostrings::readDNAStringSet(tseq)
    names(dna) <- seq_id
    if (unlink_prev) {
      unlink(tseq)
    }
    tseq <- tempfile(pattern = "sq", fileext = ".fasta")
    Biostrings::writeXStringSet(dna, tseq)
    unlink_prev <- TRUE
  }
  list(path = tseq, unlink = unlink_prev)
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
    out <- list.files(
      path = "bin",
      pattern = executable,
      recursive = TRUE,
      full.names = TRUE
    )
  }
  checkmate::assert_file_exists(out, access = "x", .var.name = executable)
  out
}

#' Try to find the usearch executable
#' @return (`character` string) the full path to the usearch executable.
find_usearch <- function() {
  find_executable("usearch")
}
