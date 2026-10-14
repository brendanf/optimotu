# Sequences used below are chosen so a non-match is distance 1 from every
# reference (all columns differ). Search drops distance 1, so those queries
# contribute no rows and must not shift the names of later hits.

far_q <- "GATAGATA"
near_q <- "ACGTACGT"
far_r <- "CCCCCCCC"
near_r <- "ACGTACGT"

# Alignment backends that share cigar, extension, and tied-hit behavior.
aligner_methods <- c("wfa2", "edlib", "ksw2", "hybrid")
# Methods that support named exact-match search without cigar/extension.
search_methods <- c(aligner_methods, "hamming")

zero_pairs <- function(out) {
  zero <- out[out$dist == 0, c("seq_id", "ref_id"), drop = FALSE]
  zero <- as.data.frame(zero, stringsAsFactors = FALSE)
  zero <- zero[order(zero$seq_id, zero$ref_id), , drop = FALSE]
  rownames(zero) <- NULL
  zero
}

test_that("seq_search names the best query and reference, not a neighbor", {
  query <- c(miss = far_q, hit_last = near_q, hit_first = far_r)
  ref <- c(first = far_r, last = near_r)
  for (method in search_methods) {
    out <- seq_search(
      query = query,
      ref = ref,
      threshold = 0.2,
      dist_config = dist_config(method = method),
      parallel_config = parallel_concurrent(1L)
    )
    expect_false(anyNA(out$seq_id), info = method)
    expect_false(anyNA(out$ref_id), info = method)
    expect_true(all(out$seq_id %in% names(query)), info = method)
    expect_true(all(out$ref_id %in% names(ref)), info = method)
    expect_equal(
      zero_pairs(out),
      data.frame(
        seq_id = c("hit_first", "hit_last"),
        ref_id = c("first", "last"),
        stringsAsFactors = FALSE
      ),
      info = method
    )
  }
})

test_that("seq_search keeps every reference tied for the best distance", {
  query <- c(q = near_q)
  ref <- c(r_first = near_r, decoy = far_r, r_last = near_r)
  for (method in aligner_methods) {
    out <- seq_search(
      query = query,
      ref = ref,
      threshold = 0.2,
      dist_config = dist_config(method = method),
      parallel_config = parallel_concurrent(1L)
    )
    expect_equal(
      zero_pairs(out),
      data.frame(
        seq_id = c("q", "q"),
        ref_id = c("r_first", "r_last"),
        stringsAsFactors = FALSE
      ),
      info = method
    )
    expect_false("decoy" %in% out$ref_id, info = method)
  }
})

test_that("seq_search query_id and ref_id replace existing names", {
  out <- seq_search(
    query = c(old_q = near_q),
    ref = c(old_r = near_r),
    threshold = 0.2,
    query_id = "new_q",
    ref_id = "new_r",
    dist_config = dist_hamming(),
    parallel_config = parallel_concurrent(1L)
  )
  expect_identical(out$seq_id, "new_q")
  expect_identical(out$ref_id, "new_r")
  expect_equal(out$dist, 0)

  expect_error(
    seq_search(
      query = c(a = near_q, b = far_q),
      ref = c(r = near_r),
      threshold = 0.2,
      query_id = "only_one",
      dist_config = dist_hamming()
    ),
    "query_id"
  )
  expect_error(
    seq_search(
      query = c(q = near_q),
      ref = c(a = near_r, b = far_r),
      threshold = 0.2,
      ref_id = "only_one",
      dist_config = dist_hamming()
    ),
    "ref_id"
  )
})

test_that("seq_search rejects unnamed sequences", {
  expect_error(
    seq_search(
      query = near_q,
      ref = c(r = near_r),
      threshold = 0.2,
      dist_config = dist_hamming()
    ),
    "named"
  )
  expect_error(
    seq_search(
      query = c(q = near_q),
      ref = near_r,
      threshold = 0.2,
      dist_config = dist_hamming()
    ),
    "named"
  )
})

test_that("seq_search seq_idx selects 1-based records", {
  query <- c(far_q = far_q, near_q = near_q)
  ref <- c(far_r = far_r, near_r = near_r)
  for (method in search_methods) {
    kept <- seq_search(
      query = query,
      ref = ref,
      threshold = 0.2,
      dist_config = dist_config(method = method),
      parallel_config = parallel_concurrent(1L),
      query_seq_idx = 2L,
      ref_seq_idx = 2L
    )
    expect_equal(
      zero_pairs(kept),
      data.frame(
        seq_id = "near_q",
        ref_id = "near_r",
        stringsAsFactors = FALSE
      ),
      info = method
    )
    expect_identical(nrow(kept), 1L, info = method)

    # Index 1 is the non-matching sequence. A 0-based reading of seq_idx
    # would select near_q / near_r instead and return a hit.
    dropped <- seq_search(
      query = query,
      ref = ref,
      threshold = 0.2,
      dist_config = dist_config(method = method),
      parallel_config = parallel_concurrent(1L),
      query_seq_idx = 1L,
      ref_seq_idx = 1L
    )
    expect_identical(nrow(dropped), 0L, info = method)
  }
})

test_that("seq_search accepts data.frame and XStringSet inputs", {
  query <- data.frame(
    seq_id = c("far_q", "near_q"),
    seq = c(far_q, near_q),
    stringsAsFactors = FALSE
  )
  ref <- Biostrings::DNAStringSet(c(far_r = far_r, near_r = near_r))
  out <- seq_search(
    query = query,
    ref = ref,
    threshold = 0.2,
    dist_config = dist_hamming(),
    parallel_config = parallel_concurrent(1L),
    query_seq_idx = 2L,
    ref_seq_idx = 2L
  )
  expect_identical(out$seq_id, "near_q")
  expect_identical(out$ref_id, "near_r")
  expect_equal(out$dist, 0)
})

test_that("seq_search reads FASTA paths", {
  fasta <- tempfile(fileext = ".fasta")
  on.exit(unlink(fasta), add = TRUE)
  writeLines(
    c(">near_q", near_q, ">far_q", far_q),
    fasta
  )
  out <- seq_search(
    query = fasta,
    ref = c(near_r = near_r, far_r = far_r),
    threshold = 0.2,
    dist_config = dist_hamming(),
    parallel_config = parallel_concurrent(1L)
  )
  expect_equal(
    zero_pairs(out),
    data.frame(
      seq_id = "near_q",
      ref_id = "near_r",
      stringsAsFactors = FALSE
    )
  )
})

test_that("seq_search file by_name keeps both orientations", {
  path <- tempfile(fileext = ".tsv")
  on.exit(unlink(path), add = TRUE)
  writeLines(
    c(
      "q_first r_first 0.10",
      "q_last r_first 0.15",
      "r_last q_last 0.20",
      "q_first q_last 0.30",
      "nope other 0.40"
    ),
    path
  )
  out <- seq_search(
    query = c(q_first = near_q, q_last = far_q),
    ref = c(r_first = far_r, r_last = near_r),
    threshold = 1,
    dist_config = dist_file(path, by_name = TRUE)
  )
  expect_identical(out$seq_id, c("q_first", "q_last", "q_last"))
  expect_identical(out$ref_id, c("r_first", "r_first", "r_last"))
  expect_equal(out$dist, c(0.10, 0.15, 0.20))
})

test_that("seq_search file indices are 0-based positions", {
  path <- tempfile(fileext = ".tsv")
  on.exit(unlink(path), add = TRUE)
  # 0 is the first sequence and 2 is the last query. An extra +1 would
  # turn the first into the second and the last into NA.
  writeLines(
    c(
      "0 0 0.10",
      "2 1 0.20",
      "1 0 0.05"
    ),
    path
  )
  out <- seq_search(
    query = c(qA = near_q, qB = far_q, qC = "AAAAAAAA"),
    ref = c(rA = far_r, rB = near_r),
    threshold = 1,
    dist_config = dist_file(path, by_name = FALSE)
  )
  expect_identical(out$seq_id, c("qA", "qC", "qB"))
  expect_identical(out$ref_id, c("rA", "rB", "rA"))
  expect_equal(out$dist, c(0.10, 0.20, 0.05))
  expect_false(anyNA(out$seq_id))
  expect_false(anyNA(out$ref_id))
})

test_that("seq_search file method rejects cigar and non-global span", {
  path <- tempfile(fileext = ".tsv")
  on.exit(unlink(path), add = TRUE)
  writeLines("q r 0.1", path)
  cfg <- dist_file(path, by_name = TRUE)
  query <- c(q = near_q)
  ref <- c(r = near_r)
  expect_error(
    seq_search(
      query = query,
      ref = ref,
      threshold = 1,
      dist_config = cfg,
      return_cigar = TRUE
    ),
    "CIGAR"
  )
  expect_error(
    seq_search(
      query = query,
      ref = ref,
      threshold = 1,
      dist_config = cfg,
      span = "extension"
    ),
    "Span"
  )
})

test_that("seq_search aligners return cigar and accept extension span", {
  query <- c(q = near_q)
  ref <- c(r = near_r)
  for (method in aligner_methods) {
    for (sp in c("global", "extension")) {
      for (cig in c(FALSE, TRUE)) {
        msg <- sprintf("method=%s span=%s cigar=%s", method, sp, cig)
        out <- seq_search(
          query = query,
          ref = ref,
          threshold = 0.2,
          dist_config = dist_config(method = method),
          parallel_config = parallel_concurrent(1L),
          return_cigar = cig,
          span = sp
        )
        expect_identical(out$seq_id, "q", info = msg)
        expect_identical(out$ref_id, "r", info = msg)
        expect_equal(out$dist, 0, info = msg)
        if (isTRUE(cig)) {
          expect_type(out$cigar, "character")
          expect_true(nzchar(out$cigar), info = msg)
        } else {
          expect_false("cigar" %in% names(out), info = msg)
        }
      }
    }
  }
  expect_error(
    seq_search(
      query = query,
      ref = ref,
      threshold = 0.2,
      dist_config = dist_hamming(),
      return_cigar = TRUE
    ),
    "cigar"
  )
  expect_error(
    seq_search(
      query = query,
      ref = ref,
      threshold = 0.2,
      dist_config = dist_hamming(),
      span = "extension"
    ),
    "span"
  )
})

test_that("seq_search extension recovers a short query on a longer ref", {
  # Global length screen rejects l1/l2 < 1 - threshold. Extension aligns the
  # query to a prefix of the longer reference, so a leading exact match
  # must score distance 0.
  query <- c(q = "ACGT")
  ref <- c(r = "ACGTAAAAAAAA")
  for (method in aligner_methods) {
    global <- seq_search(
      query = query,
      ref = ref,
      threshold = 0.2,
      dist_config = dist_config(method = method),
      parallel_config = parallel_concurrent(1L),
      span = "global"
    )
    expect_identical(nrow(global), 0L, info = method)

    extend <- seq_search(
      query = query,
      ref = ref,
      threshold = 0.2,
      dist_config = dist_config(method = method),
      parallel_config = parallel_concurrent(1L),
      span = "extension",
      return_cigar = TRUE
    )
    expect_identical(extend$seq_id, "q", info = method)
    expect_identical(extend$ref_id, "r", info = method)
    expect_equal(extend$dist, 0, info = method)
    expect_true(nzchar(extend$cigar), info = method)
  }
})

test_that("seq_search usearch returns sequence names and optional cigar", {
  usearch <- unname(Sys.which("usearch"))
  if (!nzchar(usearch) || file.access(usearch, 1L) != 0L) {
    skip("usearch not available")
  }
  query <- c(q = paste(rep("ACGT", 8L), collapse = ""))
  ref <- c(
    far = paste(rep("GGGG", 8L), collapse = ""),
    hit = paste(rep("ACGT", 8L), collapse = "")
  )
  out <- seq_search(
    query = query,
    ref = ref,
    threshold = 0.2,
    dist_config = dist_usearch(),
    parallel_config = parallel_concurrent(1L),
    fulldp = "",
    gapopen = "'*E'"
  )
  expect_identical(out$seq_id, "q")
  expect_identical(out$ref_id, "hit")
  expect_equal(out$dist, 0)
  expect_false("cigar" %in% names(out))

  out_cig <- seq_search(
    query = query,
    ref = ref,
    threshold = 0.2,
    dist_config = dist_usearch(),
    parallel_config = parallel_concurrent(1L),
    return_cigar = TRUE,
    fulldp = "",
    gapopen = "'*E'"
  )
  expect_identical(out_cig$seq_id, "q")
  expect_identical(out_cig$ref_id, "hit")
  expect_equal(out_cig$dist, 0)
  expect_type(out_cig$cigar, "character")
  expect_true(nzchar(out_cig$cigar))

  expect_error(
    seq_search(
      query = query,
      ref = ref,
      threshold = 0.2,
      dist_config = dist_usearch(),
      span = "extension"
    ),
    "Span"
  )
})

test_that("seq_search dist_wfa2 match bonus does not abort", {
  query <- c(q = "ACGT")
  ref <- c(r = "ACGA")
  out <- seq_search(
    query = query,
    ref = ref,
    threshold = 1,
    dist_config = dist_wfa2(
      match = 2L,
      mismatch = 4L,
      gap_open = 6L,
      gap_extend = 2L
    ),
    parallel_config = parallel_concurrent(1L)
  )
  expect_identical(out$seq_id, "q")
  expect_identical(out$ref_id, "r")
  expect_equal(out$dist, 0.25)
})
