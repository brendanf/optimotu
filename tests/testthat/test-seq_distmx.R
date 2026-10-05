test_that("seq_distmx_internal works", {
  testseqs <- seqs <-
    c(
      "ACGT",
      "ACGA",
      "ACGG",
      "ACGC",
      "ACGT",
      "AGCA"
    )
  names(testseqs) <- paste0("seq", seq_along(seqs))
  # USEARCH is a proprietary external binary. GitHub Actions does not
  # install it, so include that backend only when it is on PATH.
  methods <- c("wfa2", "hamming", "edlib", "ksw2")
  usearch <- unname(Sys.which("usearch"))
  if (nzchar(usearch) && file.access(usearch, 1L) == 0L) {
    methods <- c("usearch", methods)
  }
  for (method in methods) {
    for (threshold in c(0.1, 0.3, 0.55)) {
      for (threads in c(1, 4)) {
        for (detail in c("cigar", "gapstats", "none")) {
          if (method == "hamming" && detail == "cigar") {
            next
          }
          # cat("\nmethod =", method, ", threshold =", threshold, ", threads =", threads, ", detail =", detail, "\n")

          distmx <- seq_distmx(
            seq = testseqs,
            dist_config = optimotu::dist_config(method = method),
            threshold = threshold,
            parallel_config = optimotu::parallel_concurrent(threads),
            detail = detail,
            # these parameters are for usearch, to make it work in this test
            # normally it cannot be run with such short sequences
            fulldp = "",
            gapopen = "'*E'"
          )
          expect_true(is.data.frame(distmx))
          if (threshold == 0.1) {
            expect_equal(nrow(distmx), 1)
          } else if (threshold == 0.3) {
            expect_equal(nrow(distmx), 10)
          } else if (threshold == 0.5) {
            expect_equal(nrow(distmx), 15)
          }
          if (detail == 2) {
            expect_type(distmx$cigar, "character")
          } else if (detail == 1) {
            expect_type(distmx$align_length, "integer")
            expect_type(distmx$n_insert, "integer")
            expect_type(distmx$n_delete, "integer")
            expect_type(distmx$max_insert, "integer")
            expect_type(distmx$max_delete, "integer")
          }
        }
      }
    }
  }
})

test_that("constrained affine WFA keeps identity-feasible mismatch pairs", {
  # 4/20 mismatches: identity 0.2. Affine score 24 exceeds the old edit-shaped
  # cap (~6) but is under the identity-feasible affine cap.
  a <- paste(rep("A", 20L), collapse = "")
  b <- paste(c(rep("A", 16L), rep("C", 4L)), collapse = "")
  seq <- c(a = a, b = b)
  dist_config <- dist_wfa2(
    match = 0L,
    mismatch = 6L,
    gap_open = 4L,
    gap_extend = 2L
  )
  hits_con <- seq_distmx(
    seq,
    threshold = 0.25,
    dist_config = dist_config,
    constrain = TRUE,
    parallel_config = parallel_concurrent(1L)
  )
  expect_equal(nrow(hits_con), 1L)
  expect_equal(hits_con$dist2, 0.2)
})

test_that("seq_distmx character ids match 1-based indices", {
  # Names are not the index strings, so a shifted lookup cannot pass by
  # accident. "a" and "e" are identical. Threshold 1 keeps every pair,
  # including the last sequence, so an extra +1 becomes NA.
  seqs <- c(
    a = "ACGT",
    b = "ACGA",
    c = "ACGG",
    d = "ACGC",
    e = "ACGT",
    f = "AGCA"
  )
  for (method in c("wfa2", "hamming", "edlib", "ksw2")) {
    by_name <- seq_distmx(
      seqs,
      threshold = 1,
      dist_config = dist_config(method = method),
      parallel_config = parallel_concurrent(1L),
      id_is_int = FALSE
    )
    by_idx <- seq_distmx(
      seqs,
      threshold = 1,
      dist_config = dist_config(method = method),
      parallel_config = parallel_concurrent(1L),
      id_is_int = TRUE
    )
    expect_true(nrow(by_idx) > 0L, info = method)
    expect_false(anyNA(by_name$seq_id1), info = method)
    expect_false(anyNA(by_name$seq_id2), info = method)
    name_ord <- order(by_name$seq_id1, by_name$seq_id2)
    idx_ord <- order(by_idx$seq_idx1, by_idx$seq_idx2)
    expect_equal(
      by_name$seq_id1[name_ord],
      names(seqs)[by_idx$seq_idx1[idx_ord]],
      info = method
    )
    expect_equal(
      by_name$seq_id2[name_ord],
      names(seqs)[by_idx$seq_idx2[idx_ord]],
      info = method
    )
    zero <- which(by_name$dist2 == 0)
    expect_length(zero, 1L)
    expect_equal(
      sort(c(by_name$seq_id1[zero], by_name$seq_id2[zero])),
      c("a", "e"),
      info = method
    )
  }
})
