ref_evo_ml <- function(cfg, a, b) {
  states <- strsplit(cfg$states, "", fixed = TRUE)[[1L]]
  aa <- strsplit(a, "", fixed = TRUE)[[1L]]
  bb <- strsplit(b, "", fixed = TRUE)[[1L]]
  counts <- matrix(0, 20L, 20L)
  for (k in seq_along(aa)) {
    i <- match(aa[[k]], states)
    j <- match(bb[[k]], states)
    if (is.na(i) || is.na(j)) {
      next
    }
    counts[i, j] <- counts[i, j] + 1L
  }
  evec <- matrix(cfg$evec, 20L, 20L)
  evec_inv <- matrix(cfg$evec_inv, 20L, 20L)
  nll <- function(t) {
    p <- evec %*% diag(exp(cfg$lambda * t), 20L) %*% evec_inv
    p[!is.finite(p) | p < 1e-300] <- 1e-300
    s <- 0
    for (i in seq_len(20L)) {
      for (j in seq_len(20L)) {
        if (counts[i, j] == 0L) {
          next
        }
        s <- s + counts[i, j] * (log(cfg$pi[[i]]) + log(p[i, j]))
      }
    }
    -s
  }
  stats::optimize(nll, c(1e-6, cfg$t_max))$minimum
}

pair_dist <- function(seq, cfg, i, j) {
  d <- seq_distmx(
    seq,
    threshold = 0,
    dist_config = cfg,
    pairs = matrix(c(i, j), nrow = 1L),
    id_is_int = TRUE,
    verbose = FALSE
  )
  d$dist2[[1L]]
}

test_that("dist_evo accepts an abbreviation, a PAML file, and a matrix", {
  named <- dist_evo("mtart")
  path <- system.file("extdata", "evo", "mtArt.dat", package = "optimotu")
  from_file <- dist_evo(path)
  from_matrix <- dist_evo(named$exchangeability, pi = named$pi)
  expect_equal(named$method, "evo")
  expect_equal(named$exchangeability["L", "I"], 515)
  expect_equal(named$logp, from_file$logp)
  expect_equal(named$logp, from_matrix$logp)
  expect_equal(sum(named$pi), 1)
  expect_error(dist_evo("not-a-model"), "abbreviation")
  bad <- named$exchangeability
  bad[1L, 2L] <- bad[1L, 2L] + 10
  expect_error(dist_evo(bad, pi = named$pi), "symmetric")
})

test_that("mtART distances are zero for identical peptides and smaller for Ile/Leu", {
  cfg <- dist_evo("mtart", t_grid = 32L)
  same <- paste(rep("A", 31L), collapse = "")
  il <- paste0(strrep("A", 15L), "I", strrep("A", 15L))
  leu <- paste0(strrep("A", 15L), "L", strrep("A", 15L))
  cys <- paste0(strrep("A", 15L), "C", strrep("A", 15L))
  trp <- paste0(strrep("A", 15L), "W", strrep("A", 15L))
  seq <- c(same, same, il, leu, cys, trp)
  expect_equal(pair_dist(seq, cfg, 1L, 2L), 0)
  d_il <- pair_dist(seq, cfg, 3L, 4L)
  d_cw <- pair_dist(seq, cfg, 5L, 6L)
  expect_lt(d_il, d_cw)
  expect_gt(d_il, 0)
  expect_equal(d_il, ref_evo_ml(cfg, il, leu), tolerance = 1e-3)
  expect_equal(d_cw, ref_evo_ml(cfg, cys, trp), tolerance = 1e-3)
})

test_that("gaps and ambiguous residues are missing data", {
  cfg <- dist_evo("mtart")
  seq <- c("ACDE-FGHI", "ACDEXFGHI", "ACDEFFGHI")
  d_gap <- seq_distmx(
    seq,
    threshold = 0,
    dist_config = dist_evo("mtart", ignore_gaps = TRUE),
    pairs = rbind(c(1L, 3L), c(2L, 3L)),
    id_is_int = TRUE
  )
  d_keep <- seq_distmx(
    seq,
    threshold = 0,
    dist_config = dist_evo("mtart", ignore_gaps = FALSE),
    pairs = rbind(c(1L, 3L), c(2L, 3L)),
    id_is_int = TRUE
  )
  expect_equal(d_gap$dist2, d_keep$dist2)
  expect_equal(d_gap$dist2[[1L]], d_gap$dist2[[2L]])
  expect_equal(d_gap$dist2[[1L]], 0)
})

test_that("seq_distmx threads and seq_cluster agree for dist_evo", {
  seq <- c(
    s1 = "ACDEFGHIKLMNPQRSTVWY",
    s2 = "ACDEFGHIKLMNPQRSTVWY",
    s3 = "ACDEFGHIKLMNPQRSTVWW",
    s4 = "WWWWWWWWWWWWWWWWWWWW"
  )
  cfg <- dist_evo("mtart")
  d1 <- seq_distmx(
    seq,
    threshold = 5,
    dist_config = cfg,
    parallel_config = parallel_concurrent(1),
    id_is_int = TRUE
  )
  d2 <- seq_distmx(
    seq,
    threshold = 5,
    dist_config = cfg,
    parallel_config = parallel_concurrent(2),
    id_is_int = TRUE
  )
  key <- function(x) {
    x[order(pmin(x$seq_idx1, x$seq_idx2), pmax(x$seq_idx1, x$seq_idx2)), ]
  }
  expect_equal(key(d1)$dist2, key(d2)$dist2)
  cl1 <- seq_cluster(
    seq,
    dist_config = cfg,
    threshold_config = threshold_set(c(0, 0.05, 0.5)),
    clust_config = clust_slink(),
    parallel_config = parallel_merge(1)
  )
  cl2 <- seq_cluster(
    seq,
    dist_config = cfg,
    threshold_config = threshold_set(c(0, 0.05, 0.5)),
    clust_config = clust_slink(),
    parallel_config = parallel_merge(2)
  )
  expect_equal(cl1, cl2)
  expect_identical(unname(cl1["0", "s1"]), unname(cl1["0", "s2"]))
  s <- seq_search(
    query = seq[1],
    ref = seq,
    threshold = 0.2,
    dist_config = cfg
  )
  expect_true(any(s$ref_id == "s1"))
})

test_that("optimize_thresholds runs with dist_evo", {
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
    dist_config = dist_evo("mtart"),
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
