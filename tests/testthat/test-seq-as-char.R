test_that("seq_as_char preserves legacy named character output", {
  x <- c(a = "ACGT", b = "ACGA")
  expect_identical(seq_as_char(x), x)
})

test_that("seq_file is rejected for non-index-backed seq", {
  tf <- tempfile(fileext = ".fasta")
  on.exit(unlink(tf), add = TRUE)
  Biostrings::writeXStringSet(
    Biostrings::DNAStringSet(c(s1 = "ACGT", s2 = "AAAA")),
    tf
  )
  expect_error(
    seq_as_char(tf, seq_file = "/tmp/nope"),
    "`seq_file` is only allowed"
  )
})

test_that("seq_idx subsets data.frame rows in order", {
  df <- data.frame(
    seq_id = c("x", "y", "z"),
    seq = c("AAA", "CCC", "GGG"),
    stringsAsFactors = FALSE
  )
  out <- seq_as_char(df, seq_idx = c(3L, 1L))
  expect_identical(names(out), c("z", "x"))
  expect_identical(unname(out), c("GGG", "AAA"))
})

test_that("seq_idx subsets XStringSet", {
  ss <- Biostrings::DNAStringSet(c(a = "ACGT", b = "AAAA", c = "GCGC"))
  out <- seq_as_char(ss, seq_idx = c(2L, 3L))
  expect_identical(names(out), c("b", "c"))
  expect_identical(unname(out), c("AAAA", "GCGC"))
})

test_that("FASTA path with seq_idx uses indexed read (Biostrings)", {
  tf <- tempfile(fileext = ".fasta")
  on.exit(unlink(tf), add = TRUE)
  Biostrings::writeXStringSet(
    Biostrings::DNAStringSet(c(s1 = "ACGT", s2 = "AAAA", s3 = "GCGC")),
    tf
  )
  full <- seq_as_char(tf)
  sub <- seq_as_char(tf, seq_idx = c(3L, 1L))
  expect_identical(sub, full[c("s3", "s1")])
  fp <- seq_as_char(tf, seq_idx = c(2L, 1L), as = "fasta_path")
  on.exit(unlink(fp), add = TRUE)
  expect_true(file.exists(fp))
  expect_true(attr(fp, "optimotu_unlink", exact = TRUE))
  reread <- seq_as_char(fp)
  expect_identical(reread, full[c("s2", "s1")])
})

test_that("fastqindexr_index without fastqindexr errors", {
  skip_if(
    requireNamespace("fastqindexr", quietly = TRUE) &&
      utils::packageVersion("fastqindexr") >= "0.1"
  )
  idx <- structure(list(n_records = 1L), class = "fastqindexr_index")
  expect_error(seq_as_char(idx), "fastqindexr")
})

test_that("fastqindexr index subset matches full-read subset", {
  skip_if_not_installed("fastqindexr")
  skip_if(utils::packageVersion("fastqindexr") < "0.1")
  tf <- tempfile(fileext = ".fastq")
  on.exit(unlink(tf), add = TRUE)
  con <- file(tf, "wt")
  writeLines(
    c(
      "@s1", "ACGT", "+", "!!!!",
      "@s2", "AAAA", "+", "!!!!",
      "@s3", "GCGC", "+", "!!!!"
    ),
    con
  )
  close(con)
  idx <- fastqindexr::create_index(tf, type = "fastq")
  full <- seq_as_char(tf)
  sub <- seq_as_char(idx, seq_idx = c(3L, 1L))
  expect_identical(sub, full[c("s3", "s1")])
})
