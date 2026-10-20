testthat::test_that("clust_prim() errors for multi-subset clustering", {
  seq <- c(a = "ACGTACGTACGT", b = "ACGTACGTACGA", c = "ACGTACGTACCT")
  testthat::expect_error(
    seq_cluster(
      seq = seq,
      dist_config = dist_hamming(),
      threshold_config = threshold_set(c(0.1, 0.2)),
      clust_config = clust_prim(),
      parallel_config = parallel_merge(1L),
      output_type = "matrix",
      which = list(all = names(seq), ab = c("a", "b")),
      verbose = FALSE
    ),
    "does not support multi-subset clustering"
  )
  tax <- data.frame(
    seq_id = names(seq),
    kingdom = "k1",
    phylum = c("p1", "p1", "p2"),
    stringsAsFactors = FALSE
  )
  testthat::expect_error(
    optimize_thresholds(
      taxonomy = tax,
      refseq = seq,
      ranks = c("kingdom", "phylum"),
      dist_config = dist_hamming(),
      threshold_config = threshold_set(c(0.1, 0.2)),
      clust_config = clust_prim(),
      parallel_config = parallel_merge(1L),
      min_taxa = 1L,
      min_refseq = 1L,
      measures = "AMI",
      verbose = FALSE
    ),
    "does not support multi-subset clustering"
  )
})

testthat::test_that("parallel_concurrent(threads > 1) errors with clust_prim()", {
  seq <- c(a = "ACGTACGTACGT", b = "ACGTACGTACGA")
  testthat::expect_error(
    seq_cluster(
      seq = seq,
      dist_config = dist_hamming(),
      threshold_config = threshold_set(c(0.1, 0.2)),
      clust_config = clust_prim(),
      parallel_config = parallel_concurrent(4L),
      output_type = "matrix",
      verbose = FALSE
    ),
    "does not support parallel_concurrent"
  )
})

testthat::test_that("clust_prim() partitions match clust_slink() and hclust export", {
  set.seed(2)
  alphabet <- c("A", "C", "G", "T")
  make_clade <- function(base, n) {
    vapply(
      seq_len(n),
      function(i) {
        s <- strsplit(base, "", TRUE)[[1L]]
        pos <- sample.int(length(s), 3L)
        s[pos] <- sample(alphabet, 3L, replace = TRUE)
        paste(s, collapse = "")
      },
      character(1L)
    )
  }
  # Length 81 so Hamming fractions are not exactly 0.02 / 0.05 / 0.1
  # (those equalities are skipped when a worker uses d < threshold).
  bases <- replicate(
    3L,
    paste(sample(alphabet, 81L, replace = TRUE), collapse = ""),
    simplify = TRUE
  )
  seq <- unlist(lapply(bases, make_clade, n = 8L), use.names = FALSE)
  names(seq) <- sprintf("s%02d", seq_along(seq))
  thresholds <- c(0.02, 0.05, 0.1)
  dist_cfgs <- list(
    hamming = dist_hamming(),
    edlib = dist_edlib()
  )
  for (dist_name in names(dist_cfgs)) {
    slink_mat <- seq_cluster(
      seq = seq,
      dist_config = dist_cfgs[[dist_name]],
      threshold_config = threshold_set(thresholds),
      clust_config = clust_slink(),
      parallel_config = parallel_merge(1L),
      output_type = "matrix",
      verbose = FALSE
    )
    thread_set <- if (dist_name == "hamming") c(1L, 4L, 9L) else c(1L, 4L)
    for (threads in thread_set) {
      prim_mat <- seq_cluster(
        seq = seq,
        dist_config = dist_cfgs[[dist_name]],
        threshold_config = threshold_set(thresholds),
        clust_config = clust_prim(),
        parallel_config = parallel_merge(threads),
        output_type = "matrix",
        verbose = FALSE
      )
      conf <- confusion_matrix2(prim_mat, slink_mat)
      testthat::expect_equal(
        conf$FP,
        rep(0L, nrow(conf)),
        info = paste(dist_name, threads)
      )
      testthat::expect_equal(
        conf$FN,
        rep(0L, nrow(conf)),
        info = paste(dist_name, threads)
      )
    }
  }
  slink_mat <- seq_cluster(
    seq = seq,
    dist_config = dist_hamming(),
    threshold_config = threshold_set(thresholds),
    clust_config = clust_slink(),
    parallel_config = parallel_merge(1L),
    output_type = "matrix",
    verbose = FALSE
  )
  hc <- seq_cluster(
    seq = seq,
    dist_config = dist_hamming(),
    threshold_config = threshold_set(thresholds),
    clust_config = clust_prim(),
    parallel_config = parallel_merge(1L),
    output_type = "hclust",
    verbose = FALSE
  )
  testthat::expect_s3_class(hc, "hclust")
  testthat::expect_equal(nrow(hc$merge), length(seq) - 1L)
  cut <- cutree(hc, h = thresholds)
  if (is.vector(cut)) {
    cut <- matrix(cut, nrow = 1L)
  } else {
    cut <- t(cut)
  }
  slink_cut <- slink_mat
  testthat::expect_equal(nrow(cut), nrow(slink_cut))
  for (r in seq_len(nrow(slink_cut))) {
    fa <- as.integer(factor(cut[r, ], levels = unique(cut[r, ])))
    fb <- as.integer(factor(slink_cut[r, ], levels = unique(slink_cut[r, ])))
    testthat::expect_equal(fa, fb, info = rownames(slink_cut)[r])
  }
})
