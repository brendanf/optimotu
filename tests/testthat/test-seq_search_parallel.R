make_search_sequences <- function(n, len, max_mutations = 8L) {
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
  names(out) <- paste0("s", seq_len(n))
  out
}

normalize_search <- function(df) {
  if (nrow(df) == 0L) {
    return(df)
  }
  ord <- order(df$seq_id, df$ref_id, df$dist)
  df[ord, , drop = FALSE]
}

testthat::test_that("seq_search handles empty inputs across methods", {
  ref_one <- c("ACGTACGT")
  names(ref_one) <- "r1"
  query_one <- c("ACGTTCGT")
  names(query_one) <- "q1"

  input_cases <- list(
    list(
      name = "empty_query",
      query = character(0),
      ref = ref_one,
      query_id = character(0),
      ref_id = names(ref_one)
    ),
    list(
      name = "empty_ref",
      query = query_one,
      ref = character(0),
      query_id = names(query_one),
      ref_id = character(0)
    ),
    list(
      name = "both_empty",
      query = character(0),
      ref = character(0),
      query_id = character(0),
      ref_id = character(0)
    )
  )

  method_cases <- list(
    list(method = "wfa2", spans = c("global", "extension"), return_cigar = c(FALSE, TRUE)),
    list(method = "edlib", spans = c("global"), return_cigar = c(FALSE)),
    list(method = "hamming", spans = c("global"), return_cigar = c(FALSE))
  )

  for (mc in method_cases) {
    for (sp in mc$spans) {
      for (cig in mc$return_cigar) {
        for (threads in c(1L, 4L)) {
          for (ic in input_cases) {
            msg <- sprintf(
              "method=%s span=%s cigar=%s threads=%s case=%s",
              mc$method, sp, cig, threads, ic$name
            )
            out <- testthat::expect_no_error(
              seq_search(
                query = ic$query,
                ref = ic$ref,
                threshold = 0.2,
                query_id = ic$query_id,
                ref_id = ic$ref_id,
                dist_config = dist_config(method = mc$method),
                parallel_config = parallel_concurrent(threads),
                return_cigar = cig,
                span = sp
              )
            )
            testthat::expect_true(is.data.frame(out), info = msg)
            testthat::expect_equal(nrow(out), 0L, info = msg)
          }
        }
      }
    }
  }
})

testthat::test_that("seq_search tiny inputs match serial in parallel", {
  query <- c("ACGTACGT")
  ref <- c("ACGTTCGT")
  names(query) <- "q1"
  names(ref) <- "r1"

  method_cases <- list(
    list(method = "wfa2", spans = c("global", "extension"), return_cigar = c(FALSE, TRUE)),
    list(method = "edlib", spans = c("global"), return_cigar = c(FALSE)),
    list(method = "hamming", spans = c("global"), return_cigar = c(FALSE))
  )

  for (mc in method_cases) {
    for (sp in mc$spans) {
      for (cig in mc$return_cigar) {
        msg <- sprintf("method=%s span=%s cigar=%s", mc$method, sp, cig)
        serial <- testthat::expect_no_error(
          seq_search(
            query = query,
            ref = ref,
            threshold = 0.5,
            dist_config = dist_config(method = mc$method),
            parallel_config = parallel_concurrent(1L),
            return_cigar = cig,
            span = sp
          )
        )
        parallel <- testthat::expect_no_error(
          seq_search(
            query = query,
            ref = ref,
            threshold = 0.5,
            dist_config = dist_config(method = mc$method),
            parallel_config = parallel_concurrent(4L),
            return_cigar = cig,
            span = sp
          )
        )
        testthat::expect_equal(normalize_search(parallel), normalize_search(serial), info = msg)
      }
    }
  }
})

testthat::test_that("seq_search wfa2 cigar is stable with multithreading", {
  set.seed(1)
  query <- make_search_sequences(n = 30L, len = 100L, max_mutations = 10L)
  ref <- make_search_sequences(n = 35L, len = 100L, max_mutations = 10L)

  baseline <- testthat::expect_no_error(
    seq_search(
      query = query,
      ref = ref,
      threshold = 0.2,
      dist_config = dist_config(method = "wfa2"),
      parallel_config = parallel_concurrent(1L),
      return_cigar = TRUE,
      span = "global"
    )
  )
  baseline <- normalize_search(baseline)

  for (iter in seq_len(3L)) {
    msg <- sprintf("wfa2-cigar multithread iteration %s", iter)
    out <- testthat::expect_no_error(
      seq_search(
        query = query,
        ref = ref,
        threshold = 0.2,
        dist_config = dist_config(method = "wfa2"),
        parallel_config = parallel_concurrent(4L),
        return_cigar = TRUE,
        span = "global"
      )
    )
    out <- normalize_search(out)
    testthat::expect_identical(nrow(out), nrow(baseline), info = msg)
    testthat::expect_equal(out, baseline, tolerance = 0, info = msg)
  }
})
