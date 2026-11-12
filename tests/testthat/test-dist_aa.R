test_that("dist_aa Hamming p-distance matches known values", {
  seq <- c(
    A = "ACDEFGHIKL",
    B = "ACDEFGHIKM",
    C = "ACDEFGHIKL"
  )
  d <- seq_distmx(
    seq,
    threshold = 1,
    dist_config = dist_aa(),
    parallel_config = parallel_concurrent(1),
    id_is_int = TRUE
  )
  ab <- d$dist2[
    d$seq_idx1 == 1L & d$seq_idx2 == 2L | d$seq_idx1 == 2L & d$seq_idx2 == 1L
  ]
  ac <- d$dist2[
    d$seq_idx1 == 1L & d$seq_idx2 == 3L | d$seq_idx1 == 3L & d$seq_idx2 == 1L
  ]
  expect_equal(ab, 0.1)
  expect_equal(ac, 0)
})

test_that("dist_aa custom weights and ignore_gaps", {
  seq <- c(
    A = "AA-A",
    B = "AATA"
  )
  w <- matrix(
    1,
    2,
    2,
    dimnames = list(c("A", "T"), c("A", "T"))
  )
  diag(w) <- 0
  w["A", "T"] <- 0.5
  w["T", "A"] <- 0.5
  d_h <- seq_search(
    query = seq[1],
    ref = seq[2],
    threshold = 1,
    dist_config = dist_aa(ignore_gaps = TRUE)
  )
  expect_equal(d_h$dist, 0)
  d_g <- seq_search(
    query = seq[1],
    ref = seq[2],
    threshold = 1,
    dist_config = dist_aa(ignore_gaps = FALSE)
  )
  expect_equal(d_g$dist, 1 / 4)
  d_w <- seq_search(
    query = c(A = "AA"),
    ref = c(B = "AT"),
    threshold = 1,
    dist_config = dist_aa(weights = w)
  )
  expect_equal(d_w$dist, 0.25)
})

test_that("translate_aligned strips A2M lowercase and uses mitochondrial code", {
  # Invertebrate mito (5): AGA is Ser, not Arg (standard).
  nt <- c(s1 = "ATGAGAaaTAA")
  aa <- translate_aligned(nt, genetic_code = 5)
  expect_identical(unname(aa), "MS*")
  std <- translate_aligned("ATGAGATAA", genetic_code = 1)
  expect_identical(unname(std), "MR*")
})

test_that("seq_distmx pairs= returns requested rows for dist_aa", {
  seq <- c(
    a = "AAAA",
    b = "AAAT",
    c = "TTTT"
  )
  pairs <- rbind(c(1L, 2L), c(1L, 3L), c(1L, 2L))
  d <- seq_distmx(
    seq,
    threshold = 0,
    dist_config = dist_aa(),
    pairs = pairs,
    id_is_int = TRUE
  )
  expect_equal(nrow(d), 3L)
  expect_equal(d$seq_idx1, c(1L, 1L, 1L))
  expect_equal(d$seq_idx2, c(2L, 3L, 2L))
  expect_equal(d$dist2[[1]], 0.25)
  expect_equal(d$dist2[[2]], 1)
  expect_equal(d$dist2[[3]], 0.25)
})

test_that("seq_search dist_aa matches seq_distmx on the same pairs", {
  query <- c(q1 = "ACDE", q2 = "ACDF")
  ref <- c(r1 = "ACDE", r2 = "GGGG")
  s <- seq_search(
    query = query,
    ref = ref,
    threshold = 1,
    dist_config = dist_aa()
  )
  d <- seq_distmx(
    c(query, ref),
    threshold = 1,
    dist_config = dist_aa(),
    id_is_int = TRUE
  )
  q1_r1 <- d$dist2[
    (d$seq_idx1 == 1L & d$seq_idx2 == 3L) |
      (d$seq_idx1 == 3L & d$seq_idx2 == 1L)
  ]
  expect_equal(s$dist[s$seq_id == "q1" & s$ref_id == "r1"], q1_r1)
})

test_that("seq_cluster dist_aa agrees with dist_file of AA Hamming", {
  seq <- c(
    s1 = "ACDEFGHIKL",
    s2 = "ACDEFGHIKM",
    s3 = "ACDEFGHIKL",
    s4 = "QQQQQQQQQQ"
  )
  cl <- seq_cluster(
    seq,
    dist_config = dist_aa(),
    threshold_config = threshold_set(c(0, 0.2, 1)),
    clust_config = clust_slink(),
    parallel_config = parallel_merge(1)
  )
  cl2 <- seq_cluster(
    seq,
    dist_config = dist_aa(),
    threshold_config = threshold_set(c(0, 0.2, 1)),
    clust_config = clust_tree(),
    parallel_config = parallel_concurrent(1)
  )
  expect_equal(cl["0"], cl2["0"])
  expect_identical(unname(cl["0", "s1"]), unname(cl["0", "s3"]))
  expect_false(identical(unname(cl["0", "s1"]), unname(cl["0", "s4"])))
})

test_that("estimate_aa_packed_mb is positive for method aa", {
  expect_gt(optimotu:::estimate_aa_packed_mb(100, 200), 0)
})

test_that("optimize_thresholds runs with dist_aa", {
  taxonomy <- data.frame(
    seq_id = paste0("s", 1:6),
    kingdom = "k1",
    phylum = c("p1", "p1", "p1", "p2", "p2", "p2"),
    class = c("c1", "c1", "c2", "c3", "c3", "c4"),
    stringsAsFactors = FALSE
  )
  seq <- c(
    s1 = "ACDEFGHIKL",
    s2 = "ACDEFGHIKM",
    s3 = "ACDEFGHIQL",
    s4 = "QQQQQQQQQQ",
    s5 = "QQQQQQQQQR",
    s6 = "QQQQQQQQQS"
  )
  out <- optimize_thresholds(
    taxonomy = taxonomy,
    refseq = seq,
    ranks = c("kingdom", "phylum", "class"),
    dist_config = dist_aa(),
    threshold_config = threshold_set(c(0, 0.2, 1)),
    clust_config = clust_slink(),
    parallel_config = parallel_merge(1),
    min_taxa = 2L,
    min_refseq = 2L,
    measures = "MCC",
    verbose = FALSE
  )
  expect_gt(nrow(out), 0L)
  expect_true(all(out$measure == "MCC"))
})

ref_aa_pair <- function(a, b, ignore_gaps = TRUE, mismatch = 1) {
  aa <- strsplit(a, "", fixed = TRUE)[[1L]]
  bb <- strsplit(b, "", fixed = TRUE)[[1L]]
  span <- function(x) {
    nz <- which(x != "-")
    if (length(nz) == 0L) {
      return(c(NA_integer_, NA_integer_))
    }
    c(nz[[1L]], nz[[length(nz)]])
  }
  sa <- span(aa)
  sb <- span(bb)
  if (anyNA(sa) || anyNA(sb)) {
    return(1)
  }
  lo <- max(sa[[1L]], sb[[1L]])
  hi <- min(sa[[2L]], sb[[2L]])
  if (lo > hi) {
    return(1)
  }
  ok <- 0
  acc <- 0
  for (k in lo:hi) {
    ga <- aa[[k]] == "-"
    gb <- bb[[k]] == "-"
    if (ga && gb) {
      next
    }
    if (ignore_gaps && (ga || gb)) {
      next
    }
    ok <- ok + 1
    if (aa[[k]] != bb[[k]]) {
      acc <- acc + mismatch
    }
  }
  if (ok == 0) {
    1
  } else {
    acc / ok
  }
}

test_that("vector kernels match scalar distances on long alignments", {
  base <- paste(rep("ACDEFGHIKLMNPQRSTVWY", 2), collapse = "")
  chars0 <- strsplit(base, "", fixed = TRUE)[[1L]]
  seq <- character(20)
  for (i in seq_along(seq)) {
    chars <- chars0
    chars[[i]] <- "Q"
    if (i %% 5L == 0L) {
      chars[[30L]] <- "-"
    }
    if (i == 1L) {
      chars[1:3] <- "-"
    }
    seq[[i]] <- paste(chars, collapse = "")
  }
  names(seq) <- paste0("s", seq_along(seq))
  d <- seq_distmx(
    seq,
    threshold = 1,
    dist_config = dist_aa(),
    parallel_config = parallel_concurrent(1),
    id_is_int = TRUE
  )
  pair_dist <- function(i, j) {
    hit <- d$dist2[
      (d$seq_idx1 == i & d$seq_idx2 == j) |
        (d$seq_idx1 == j & d$seq_idx2 == i)
    ]
    expect_length(hit, 1L)
    hit
  }
  for (i in c(1L, 5L, 19L)) {
    for (j in c(2L, 10L, 20L)) {
      if (i == j) {
        next
      }
      expect_equal(
        pair_dist(i, j),
        ref_aa_pair(seq[[i]], seq[[j]])
      )
    }
  }
  lets <- LETTERS
  w <- matrix(0.25, 26, 26, dimnames = list(lets, lets))
  diag(w) <- 0
  dw <- seq_distmx(
    seq,
    threshold = 1,
    dist_config = dist_aa(weights = w),
    parallel_config = parallel_concurrent(1),
    id_is_int = TRUE
  )
  pair_w <- function(i, j) {
    hit <- dw$dist2[
      (dw$seq_idx1 == i & dw$seq_idx2 == j) |
        (dw$seq_idx1 == j & dw$seq_idx2 == i)
    ]
    expect_length(hit, 1L)
    hit
  }
  for (i in c(1L, 5L, 19L)) {
    for (j in c(2L, 10L, 20L)) {
      if (i == j) {
        next
      }
      expect_equal(
        pair_w(i, j),
        ref_aa_pair(seq[[i]], seq[[j]], mismatch = 0.25)
      )
    }
  }
})
