same_cluster_partition <- function(a, b) {
  cols <- intersect(colnames(a), colnames(b))
  testthat::expect_gt(length(cols), 0L)
  a <- a[, cols, drop = FALSE]
  b <- b[, cols, drop = FALSE]
  testthat::expect_equal(dim(a), dim(b))
  mismatches <- integer()
  for (r in seq_len(nrow(a))) {
    fa <- as.integer(factor(a[r, ], levels = unique(a[r, ])))
    fb <- as.integer(factor(b[r, ], levels = unique(b[r, ])))
    if (!identical(fa, fb)) {
      mismatches <- c(mismatches, r)
    }
  }
  testthat::expect_equal(
    mismatches,
    integer(),
    info = paste(
      "differing rows:",
      paste(rownames(a)[mismatches], collapse = ",")
    )
  )
  invisible(TRUE)
}

make_nested_seq_fixture <- function(n_per = 20L, seq_len = 120L, n_mut = 3L) {
  set.seed(1)
  alphabet <- c("A", "C", "G", "T")
  make_seqs <- function(base, n) {
    vapply(
      seq_len(n),
      function(i) {
        s <- strsplit(base, "", fixed = TRUE)[[1]]
        pos <- sample.int(length(s), n_mut)
        s[pos] <- sample(alphabet, n_mut, replace = TRUE)
        paste(s, collapse = "")
      },
      character(1)
    )
  }
  bases <- replicate(
    3,
    paste(sample(alphabet, seq_len, replace = TRUE), collapse = ""),
    simplify = TRUE
  )
  seq <- unlist(lapply(bases, make_seqs, n = n_per), use.names = FALSE)
  names(seq) <- sprintf("s%03d", seq_along(seq))
  taxonomy <- data.frame(
    seq_id = names(seq),
    kingdom = "k1",
    phylum = rep(c("p1", "p2", "p3"), each = n_per),
    class = rep(rep(c("c1", "c2"), each = n_per / 2L), 3),
    stringsAsFactors = FALSE
  )
  list(seq = seq, taxonomy = taxonomy)
}

cluster_as_matrix <- function(
  seq,
  which = TRUE,
  dist_config = dist_hamming(),
  clust_config = clust_tree(),
  parallel_config = parallel_concurrent(1L)
) {
  seq_cluster(
    seq = seq,
    dist_config = dist_config,
    threshold_config = threshold_uniform(0, 0.4, 0.01),
    clust_config = clust_config,
    parallel_config = parallel_config,
    output_type = "matrix",
    which = which,
    verbose = FALSE
  )
}

testthat::test_that("multi-subset clustering matches independent subset clustering", {
  fx <- make_nested_seq_fixture()
  ts <- summarize_by_rank(
    fx$taxonomy,
    c("kingdom", "phylum", "class"),
    "seq_id"
  )
  ts <- ts[
    (ts$n_seq >= 4L & ts$n_taxa >= 2L) |
      ts$supertaxon %in% fx$taxonomy$kingdom,
  ]
  which_sets <- ts$seq_id
  names(which_sets) <- paste(ts$rank, ts$superrank, ts$supertaxon, sep = "/")

  configs <- list(
    list(
      dist = dist_hamming(),
      clust = clust_tree(),
      par = parallel_concurrent(1L)
    ),
    list(
      dist = dist_edlib(),
      clust = clust_tree(),
      par = parallel_concurrent(1L)
    ),
    list(
      dist = dist_hamming(),
      clust = clust_matrix(),
      par = parallel_concurrent(1L)
    ),
    list(
      dist = dist_hamming(),
      clust = clust_tree(),
      par = parallel_merge(1L)
    ),
    list(
      dist = dist_hamming(),
      clust = clust_tree(),
      par = parallel_merge(2L)
    )
  )

  for (cfg in configs) {
    info <- paste(
      cfg$dist$method,
      cfg$clust$method,
      cfg$par$method,
      cfg$par$threads,
      sep = "/"
    )
    multi <- cluster_as_matrix(
      fx$seq,
      which = which_sets,
      dist_config = cfg$dist,
      clust_config = cfg$clust,
      parallel_config = cfg$par
    )
    testthat::expect_length(multi, length(which_sets))
    for (nm in names(which_sets)) {
      independent <- cluster_as_matrix(
        fx$seq[which_sets[[nm]]],
        which = TRUE,
        dist_config = cfg$dist,
        clust_config = cfg$clust,
        parallel_config = parallel_concurrent(1L)
      )
      same_cluster_partition(multi[[nm]], independent)
    }
  }
})

testthat::test_that("crossing subsets do not change clustering of the full set", {
  fx <- make_nested_seq_fixture()
  n <- length(fx$seq)
  which_sets <- list(
    all = names(fx$seq),
    # every other sequence: overlaps all three phyla, not nested
    odd = names(fx$seq)[seq(1L, n, by = 2L)],
    even = names(fx$seq)[seq(2L, n, by = 2L)]
  )
  configs <- list(
    list(dist = dist_hamming(), par = parallel_concurrent(1L)),
    list(dist = dist_hamming(), par = parallel_concurrent(2L)),
    list(dist = dist_hamming(), par = parallel_merge(1L)),
    list(dist = dist_edlib(), par = parallel_concurrent(1L)),
    list(dist = dist_edlib(), par = parallel_merge(1L))
  )
  for (cfg in configs) {
    multi <- cluster_as_matrix(
      fx$seq,
      which = which_sets,
      dist_config = cfg$dist,
      parallel_config = cfg$par
    )
    for (nm in names(which_sets)) {
      independent <- cluster_as_matrix(
        fx$seq[which_sets[[nm]]],
        dist_config = cfg$dist
      )
      same_cluster_partition(multi[[nm]], independent)
    }
  }
})
