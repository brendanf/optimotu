testthat::test_that("optimize_thresholds accepts annotated FASTA via seq_names", {
  seq <- c(
    "ACGTACGTACGT",
    "ACGTACGTTCGT",
    "ACGTTCGTACGT",
    "ACGTTCGTTCGT",
    "ACGTACGGACGT",
    "ACGTACGGTCGT"
  )
  seq_ids <- paste0("s", seq_along(seq))
  # SINTAX-style headers: Biostrings names are the full header line
  headers <- paste0(
    seq_ids,
    ";tax=k:k1,p:",
    rep(c("p1", "p1", "p1", "p2", "p2", "p2"), 1),
    ",c:",
    rep(c("c1", "c1", "c2", "c3", "c3", "c4"), 1)
  )
  names(seq) <- headers

  tf <- tempfile(fileext = ".fasta")
  on.exit(unlink(tf), add = TRUE)
  Biostrings::writeXStringSet(Biostrings::BStringSet(seq), tf)

  taxonomy <- data.frame(
    seq_id = seq_ids,
    kingdom = rep("k1", length(seq)),
    phylum = rep(c("p1", "p1", "p1", "p2", "p2", "p2"), 1),
    class = rep(c("c1", "c1", "c2", "c3", "c3", "c4"), 1),
    stringsAsFactors = FALSE
  )

  common_args <- list(
    taxonomy = taxonomy,
    ranks = c("kingdom", "phylum", "class"),
    dist_config = dist_hamming(),
    threshold_config = threshold_set(c(0.05, 0.10, 0.20)),
    clust_config = clust_slink(),
    parallel_config = parallel_concurrent(threads = 1L),
    min_taxa = 2L,
    min_refseq = 3L,
    measures = "MCC",
    verbose = FALSE
  )

  # Full annotated headers do not match taxonomy$seq_id
  testthat::expect_error(
    do.call(optimize_thresholds, c(list(refseq = tf), common_args)),
    "set equal|Assertion|seq_id|names",
    ignore.case = TRUE
  )

  file_ids <- sub(";.*$", "", names(Biostrings::fasta.seqlengths(tf)))
  seq_idx <- match(taxonomy$seq_id, file_ids)
  testthat::expect_false(anyNA(seq_idx))

  out <- do.call(
    optimize_thresholds,
    c(
      list(
        refseq = tf,
        seq_idx = seq_idx,
        seq_names = taxonomy$seq_id
      ),
      common_args
    )
  )
  testthat::expect_s3_class(out, "data.frame")
  testthat::expect_true(nrow(out) > 0L)
  testthat::expect_true("threshold" %in% names(out))
})
