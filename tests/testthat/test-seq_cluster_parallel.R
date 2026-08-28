simulate_clustered_sequences <- function(
  n_clusters = 5L,
  members_per_cluster = 5L,
  seq_len = 220L,
  centroid_sub_frac = 0.10,
  member_sub_frac = 0.03
) {
  alphabet <- c("A", "C", "G", "T")

  mutate_substitutions <- function(seq, frac) {
    out <- strsplit(seq, split = "", fixed = TRUE)[[1]]
    n_mut <- max(1L, as.integer(round(length(out) * frac)))
    # only mutate every fifth position, and not within 3 positions of either end
    # this helps to endure that the optimal alignment is gap-free even when the
    # alignment algorithm allows gaps.
    idx <- sample.int((length(out) - 6L) %/% 5L, n_mut, replace = FALSE) *
      5L +
      3L
    for (i in idx) {
      out[i] <- sample(alphabet[alphabet != out[i]], size = 1L)
    }
    paste0(out, collapse = "")
  }

  ancestor <- paste0(sample(alphabet, seq_len, replace = TRUE), collapse = "")
  centroids <- vapply(
    seq_len(n_clusters),
    function(i) mutate_substitutions(ancestor, centroid_sub_frac),
    character(1)
  )
  seq <- unlist(
    lapply(seq_len(n_clusters), function(cluster_id) {
      vapply(
        seq_len(members_per_cluster),
        function(member_id) {
          mutate_substitutions(centroids[[cluster_id]], member_sub_frac)
        },
        character(1)
      )
    }),
    use.names = FALSE
  )
  names(seq) <- sprintf(
    "c%02d_s%02d",
    rep(seq_len(n_clusters), each = members_per_cluster),
    rep(seq_len(members_per_cluster), times = n_clusters)
  )
  seq
}

simulate_hierarchical_sequences <- function(
  ancestor_len = 160L,
  branching = 10L,
  substitutions = c(18L, 9L, 5L, 2L),
  seed = 1L
) {
  set.seed(seed)
  alphabet <- c("A", "C", "G", "T")
  mutate_sequence <- function(seq, n_subs) {
    chars <- strsplit(seq, "", fixed = TRUE)[[1]]
    n_subs <- min(n_subs, length(chars))
    if (n_subs <= 0L) {
      return(seq)
    }
    pos <- sample.int(length(chars), size = n_subs, replace = FALSE)
    for (i in pos) {
      chars[i] <- sample(setdiff(alphabet, chars[i]), size = 1L)
    }
    paste(chars, collapse = "")
  }

  current <- paste(
    sample(alphabet, ancestor_len, replace = TRUE),
    collapse = ""
  )
  for (n_subs in substitutions) {
    next_level <- character(length(current) * branching)
    k <- 1L
    for (s in current) {
      for (j in seq_len(branching)) {
        next_level[k] <- mutate_sequence(s, n_subs)
        k <- k + 1L
      }
    }
    current <- next_level
  }
  names(current) <- sprintf("sim_%05d", seq_along(current))
  current
}

normalize_seq_distmx <- function(df) {
  out <- data.frame(
    seq_idx1 = as.integer(df$seq_idx1),
    seq_idx2 = as.integer(df$seq_idx2),
    dist = as.numeric(df$dist2)
  )
  out <- out[order(out$seq_idx1, out$seq_idx2), , drop = FALSE]
  rownames(out) <- NULL
  out
}

distmx_to_dense <- function(df, n) {
  out <- matrix(0, nrow = n, ncol = n)
  if (nrow(df) > 0L) {
    for (k in seq_len(nrow(df))) {
      i <- df$seq_idx1[[k]]
      j <- df$seq_idx2[[k]]
      d <- df$dist[[k]]
      out[i, j] <- d
      out[j, i] <- d
    }
  }
  out
}

hclust2matrix <- function(distmx, thresholds) {
  out <- hclust(pmax(distmx - 1e-8, 0), method = "single") |>
    lapply(thresholds, cutree, k = NULL, tree = _) |>
    lapply(function(x) {
      for (i in seq_along(x)) {
        if (x[i] < i) {
          x[x >= i] <- x[x >= i] + 1L
        }
      }
      x
    }) |>
    do.call(rbind, args = _)
  out - 1L
}

fixture <- local({
  set.seed(11)
  seq <- simulate_clustered_sequences()
  n <- length(seq)
  idx <- seq_len(n)
  subsets <- list(
    idx[idx %% 2L == 1L],
    idx[idx %% 2L == 0L],
    sort(sample.int(n, size = floor(n * 0.7), replace = FALSE))
  )
  thresholds <- 0:10 / 10

  list(
    seq = seq,
    n = n,
    subsets = subsets,
    thresholds = thresholds,
    threshold_cfg = threshold_set(thresholds)
  )
})

testthat::test_that("seq_cluster handles one-sequence inputs across methods", {
  threshold_cfg <- threshold_set(c(0.1, 0.2))
  seq_one <- c("ACGTACGTACGT")
  names(seq_one) <- "seq1"

  for (method in c("wfa2", "edlib", "hamming")) {
    for (threads in c(1L, 4L)) {
      msg <- sprintf("method=%s threads=%s", method, threads)
      out <- testthat::expect_no_error(
        seq_cluster(
          seq = seq_one,
          dist_config = dist_config(method = method),
          threshold_config = threshold_cfg,
          clust_config = clust_tree(),
          parallel_config = parallel_concurrent(threads),
          output_type = "matrix",
          which = TRUE,
          verbose = FALSE
        )
      )
      testthat::expect_true(is.matrix(out), info = msg)
      testthat::expect_equal(nrow(out), 2L, info = msg)
      testthat::expect_equal(ncol(out), 1L, info = msg)
    }
  }
})

testthat::test_that("seq_cluster handles zero-sequence inputs across methods", {
  threshold_cfg <- threshold_set(c(0.1, 0.2))

  for (method in c("wfa2", "edlib", "hamming")) {
    for (threads in c(1L, 4L)) {
      msg <- sprintf("method=%s threads=%s", method, threads)
      out <- testthat::expect_no_error(
        seq_cluster(
          seq = character(0),
          seq_id = character(0),
          dist_config = dist_config(method = method),
          threshold_config = threshold_cfg,
          clust_config = clust_tree(),
          parallel_config = parallel_concurrent(threads),
          output_type = "matrix",
          which = TRUE,
          verbose = FALSE
        )
      )
      testthat::expect_true(is.matrix(out), info = msg)
      testthat::expect_equal(nrow(out), 2L, info = msg)
      testthat::expect_equal(ncol(out), 0L, info = msg)
    }
  }
})

testthat::test_that("seq_distmx backends agree for simulated substitutions", {
  dist_cfg <- list(
    wfa2 = dist_config(method = "wfa2"),
    edlib = dist_config(method = "edlib"),
    hamming = dist_config(method = "hamming")
  )

  out <- lapply(names(dist_cfg), function(method) {
    seq_distmx(
      seq = fixture$seq,
      threshold = 1.0,
      dist_config = dist_cfg[[method]],
      parallel_config = parallel_concurrent(1L),
      span = "global",
      constrain = TRUE,
      id_is_int = TRUE
    )
  })
  names(out) <- names(dist_cfg)
  out <- lapply(out, normalize_seq_distmx)

  testthat::expect_equal(out$wfa2$seq_idx1, out$edlib$seq_idx1)
  testthat::expect_equal(out$wfa2$seq_idx2, out$edlib$seq_idx2)
  testthat::expect_equal(
    out$wfa2$dist,
    out$edlib$dist,
    tolerance = 1e-12
  )

  testthat::expect_equal(out$edlib$seq_idx1, out$hamming$seq_idx1)
  testthat::expect_equal(out$edlib$seq_idx2, out$hamming$seq_idx2)
  testthat::expect_equal(
    out$edlib$dist,
    out$hamming$dist,
    tolerance = 1e-12
  )
})

hclust_matrix <- NULL
distmx_hamming <- NULL
testthat::test_that("hclust truth set can be generated from simulated sequences", {
  distmx_hamming <<- normalize_seq_distmx(
    seq_distmx(
      seq = fixture$seq,
      threshold = 1.0,
      dist_config = dist_config(method = "hamming"),
      parallel_config = parallel_concurrent(1L),
      span = "global",
      constrain = TRUE,
      id_is_int = TRUE
    )
  )

  dense <- distmx_to_dense(distmx_hamming, fixture$n)
  hclust_matrix <<- hclust2matrix(as.dist(dense), fixture$thresholds)
  testthat::expect_true(is.matrix(hclust_matrix))
})

testthat::test_that("parallel_merge does not deadlock for tree+hamming", {
  for (threads in c(1L, 2L)) {
    msg <- sprintf("parallel_merge(%d)", threads)
    out <- testthat::expect_no_error(
      seq_cluster(
        seq = fixture$seq,
        dist_config = dist_hamming(),
        threshold_config = fixture$threshold_cfg,
        clust_config = clust_tree(),
        parallel_config = parallel_merge(threads),
        output_type = "matrix",
        verbose = FALSE
      )
    )
    testthat::expect_true(is.matrix(out), info = msg)
    testthat::expect_equal(ncol(out), fixture$n, info = msg)
  }
})

testthat::test_that("parallel_merge handles slink+hamming on hierarchical sequences", {
  seq_full <- simulate_hierarchical_sequences(seed = 1L)
  set.seed(1)
  seq_subset <- seq_full[sample.int(
    length(seq_full),
    size = 500L,
    replace = FALSE
  )]
  out <- testthat::expect_no_error(
    seq_cluster(
      seq = seq_subset,
      dist_config = dist_hamming(),
      threshold_config = threshold_set(c(0.02, 0.04)),
      clust_config = clust_slink(),
      parallel_config = parallel_merge(1L),
      output_type = "matrix",
      verbose = FALSE
    )
  )
  testthat::expect_true(is.matrix(out))
  testthat::expect_equal(ncol(out), length(seq_subset))
})

testthat::test_that("parallel_merge handles slink+hamming with multiple threads", {
  seq_full <- simulate_hierarchical_sequences(seed = 1L)
  set.seed(1)
  seq_subset <- seq_full[sample.int(
    length(seq_full),
    size = 500L,
    replace = FALSE
  )]
  for (threads in c(1L, 2L)) {
    out <- testthat::expect_no_error(
      seq_cluster(
        seq = seq_subset,
        dist_config = dist_hamming(),
        threshold_config = threshold_set(c(0.02, 0.04)),
        clust_config = clust_slink(),
        parallel_config = parallel_merge(threads),
        output_type = "matrix",
        verbose = FALSE
      )
    )
    testthat::expect_equal(ncol(out), length(seq_subset))
  }
})

algorithms <- list(
  tree = clust_tree(),
  slink = clust_slink(),
  matrix_binary = clust_matrix(binary_search = TRUE, fill_method = "binary"),
  index = clust_index()
)

parallels <- list(
  concurrent1 = parallel_concurrent(1L),
  concurrent4 = parallel_concurrent(4L),
  merge1 = parallel_merge(1L),
  merge2 = parallel_merge(2L)
)

for (algorithm_name in names(algorithms)) {
  for (parallel_name in names(parallels)) {
    # SLINK requires a single ordered pair stream; concurrent(>1) is invalid.
    if (algorithm_name == "slink" && parallel_name == "concurrent4") {
      next
    }
    testthat::test_that(
      sprintf(
        "seq_cluster %s with %s matches hclust truth",
        algorithm_name,
        parallel_name
      ),
      {
        out <- testthat::expect_no_error(
          seq_cluster(
            seq = fixture$seq,
            dist_config = dist_hamming(),
            threshold_config = fixture$threshold_cfg,
            clust_config = algorithms[[algorithm_name]],
            parallel_config = parallels[[parallel_name]],
            output_type = "matrix",
            verbose = FALSE
          )
        )
        conf_mat <- confusion_matrix2(out, hclust_matrix)
        testthat::expect_equal(conf_mat$FP, rep(0L, nrow(conf_mat)))
        testthat::expect_equal(conf_mat$FN, rep(0L, nrow(conf_mat)))
      }
    )
  }
}
