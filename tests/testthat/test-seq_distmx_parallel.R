make_test_sequences <- function(n, len, max_mutations = 8L) {
  alphabet <- c("A", "C", "G", "T")
  template <- sample(alphabet, size = len, replace = TRUE)
  out <- character(n)
  for (i in seq_len(n)) {
    s <- template
    n_mut <- sample.int(max_mutations + 1L, size = 1L) - 1L
    if (n_mut > 0L) {
      idx <- sample.int(len, size = n_mut, replace = FALSE)
      s[idx] <- sample(alphabet, size = n_mut, replace = TRUE)
    }
    out[[i]] <- paste0(s, collapse = "")
  }
  names(out) <- as.character(seq_len(n))
  out
}

normalize_distmx <- function(df) {
  if (nrow(df) == 0L) {
    return(df)
  }
  ord <- order(df$seq_idx1, df$seq_idx2)
  df[ord, , drop = FALSE]
}

testthat::test_that("seq_distmx handles one-sequence inputs across methods", {
  seq_one <- c("ACGTACGTACGT")
  names(seq_one) <- "1"

  test_cases <- list(
    list(
      method = "wfa2",
      spans = c("global", "extension"),
      details = c("none", "gapstats", "cigar"),
      constrain = TRUE
    ),
    list(
      method = "edlib",
      spans = c("global"),
      details = c("none", "gapstats", "cigar"),
      constrain = TRUE
    ),
    list(
      method = "hamming",
      spans = c("global"),
      details = c("none", "gapstats"),
      constrain = TRUE
    )
  )

  for (tc in test_cases) {
    for (sp in tc$spans) {
      for (det in tc$details) {
        for (thr in c(1L, 4L)) {
          msg <- sprintf(
            "method=%s span=%s details=%s threads=%s",
            tc$method, sp, det, thr
          )
          out <- testthat::expect_no_error(
            seq_distmx(
              seq = seq_one,
              threshold = 0.2,
              dist_config = dist_config(method = tc$method),
              parallel_config = parallel_concurrent(thr),
              details = det,
              span = sp,
              constrain = tc$constrain,
              id_is_int = TRUE
            )
          )
          testthat::expect_true(is.data.frame(out), info = msg)
          testthat::expect_equal(nrow(out), 0L, info = msg)
        }
      }
    }
  }
})

testthat::test_that("seq_distmx handles zero-sequence inputs across methods", {
  seq_none <- character(0)

  test_cases <- list(
    list(
      method = "wfa2",
      spans = c("global", "extension"),
      details = c("none", "gapstats", "cigar"),
      constrain = TRUE
    ),
    list(
      method = "edlib",
      spans = c("global"),
      details = c("none", "gapstats", "cigar"),
      constrain = TRUE
    ),
    list(
      method = "hamming",
      spans = c("global"),
      details = c("none", "gapstats"),
      constrain = TRUE
    )
  )

  for (tc in test_cases) {
    for (sp in tc$spans) {
      for (det in tc$details) {
        for (thr in c(1L, 4L)) {
          msg <- sprintf(
            "method=%s span=%s details=%s threads=%s",
            tc$method, sp, det, thr
          )
          out <- testthat::expect_no_error(
            seq_distmx(
              seq = seq_none,
              threshold = 0.2,
              dist_config = dist_config(method = tc$method),
              parallel_config = parallel_concurrent(thr),
              details = det,
              span = sp,
              constrain = tc$constrain,
              id_is_int = TRUE
            )
          )
          testthat::expect_true(is.data.frame(out), info = msg)
          testthat::expect_equal(nrow(out), 0L, info = msg)
        }
      }
    }
  }
})

testthat::test_that("seq_distmx tiny two-sequence inputs match serial in parallel", {
  seq_two <- c("ACGTACGT", "ACGTTCGT")
  names(seq_two) <- c("1", "2")

  test_cases <- list(
    list(
      method = "wfa2",
      spans = c("global", "extension"),
      details = c("none", "gapstats", "cigar"),
      constrain = TRUE
    ),
    list(
      method = "edlib",
      spans = c("global"),
      details = c("none", "gapstats", "cigar"),
      constrain = TRUE
    ),
    list(
      method = "hamming",
      spans = c("global"),
      details = c("none", "gapstats"),
      constrain = TRUE
    )
  )

  for (tc in test_cases) {
    for (sp in tc$spans) {
      for (det in tc$details) {
        msg <- sprintf("method=%s span=%s details=%s", tc$method, sp, det)
        serial <- testthat::expect_no_error(
          seq_distmx(
            seq = seq_two,
            threshold = 0.5,
            dist_config = dist_config(method = tc$method),
            parallel_config = parallel_concurrent(1L),
            details = det,
            span = sp,
            constrain = tc$constrain,
            id_is_int = TRUE
          )
        )
        parallel <- testthat::expect_no_error(
          seq_distmx(
            seq = seq_two,
            threshold = 0.5,
            dist_config = dist_config(method = tc$method),
            parallel_config = parallel_concurrent(4L),
            details = det,
            span = sp,
            constrain = tc$constrain,
            id_is_int = TRUE
          )
        )

        testthat::expect_equal(normalize_distmx(parallel), normalize_distmx(serial), info = msg)
      }
    }
  }
})

testthat::test_that("seq_distmx edlib gapstats is stable with multithreading", {
  seqs <- make_test_sequences(n = 200L, len = 120L, max_mutations = 10L)
  threshold <- 0.2
  cfg <- dist_config(method = "edlib")

  baseline <- testthat::expect_no_error(
    seq_distmx(
      seq = seqs,
      threshold = threshold,
      dist_config = cfg,
      parallel_config = parallel_concurrent(1L),
      details = "gapstats",
      span = "global",
      constrain = TRUE,
      id_is_int = TRUE
    )
  )
  baseline <- normalize_distmx(baseline)

  for (iter in seq_len(3L)) {
    msg <- sprintf("edlib-gapstats multithread iteration %s", iter)
    out <- testthat::expect_no_error(
      seq_distmx(
        seq = seqs,
        threshold = threshold,
        dist_config = cfg,
        parallel_config = parallel_concurrent(4L),
        details = "gapstats",
        span = "global",
        constrain = TRUE,
        id_is_int = TRUE
      )
    )
    out <- normalize_distmx(out)
    testthat::expect_identical(nrow(out), nrow(baseline), info = msg)
    testthat::expect_equal(out, baseline, tolerance = 0, info = msg)
  }
})

testthat::test_that("seq_distmx edlib cigar is stable with multithreading", {
  seqs <- make_test_sequences(n = 200L, len = 120L, max_mutations = 10L)
  threshold <- 0.2
  cfg <- dist_config(method = "edlib")

  baseline <- testthat::expect_no_error(
    seq_distmx(
      seq = seqs,
      threshold = threshold,
      dist_config = cfg,
      parallel_config = parallel_concurrent(1L),
      details = "cigar",
      span = "global",
      constrain = TRUE,
      id_is_int = TRUE
    )
  )
  baseline <- normalize_distmx(baseline)

  for (iter in seq_len(3L)) {
    msg <- sprintf("edlib-cigar multithread iteration %s", iter)
    out <- testthat::expect_no_error(
      seq_distmx(
        seq = seqs,
        threshold = threshold,
        dist_config = cfg,
        parallel_config = parallel_concurrent(4L),
        details = "cigar",
        span = "global",
        constrain = TRUE,
        id_is_int = TRUE
      )
    )
    out <- normalize_distmx(out)
    testthat::expect_identical(nrow(out), nrow(baseline), info = msg)
    testthat::expect_equal(out, baseline, tolerance = 0, info = msg)
  }
})

testthat::test_that("seq_distmx wfa2 cigar is stable with multithreading", {
  seqs <- make_test_sequences(n = 200L, len = 120L, max_mutations = 10L)
  threshold <- 0.2
  cfg <- dist_config(method = "wfa2")

  baseline <- testthat::expect_no_error(
    seq_distmx(
      seq = seqs,
      threshold = threshold,
      dist_config = cfg,
      parallel_config = parallel_concurrent(1L),
      details = "cigar",
      span = "global",
      constrain = TRUE,
      id_is_int = TRUE
    )
  )
  baseline <- normalize_distmx(baseline)

  for (iter in seq_len(3L)) {
    msg <- sprintf("wfa2-cigar multithread iteration %s", iter)
    out <- testthat::expect_no_error(
      seq_distmx(
        seq = seqs,
        threshold = threshold,
        dist_config = cfg,
        parallel_config = parallel_concurrent(4L),
        details = "cigar",
        span = "global",
        constrain = TRUE,
        id_is_int = TRUE
      )
    )
    out <- normalize_distmx(out)
    testthat::expect_identical(nrow(out), nrow(baseline), info = msg)
    testthat::expect_equal(out, baseline, tolerance = 0, info = msg,
     ignore_attr = TRUE)
  }
})
