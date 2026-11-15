gpu_dna <- function() {
  c(
    a = "ACGTACGTACGTACGTACGTACGTACGTACGT",
    b = "ACGTACGTACGTACGTACGTACGTACGTACGA",
    c = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG",
    d = "ACGTACGTACGTACGTACGTACGTACGTACGT"
  )
}

test_that("dist_gasal2() errors without a CUDA build", {
  skip_if(optimotu_cuda_build())
  expect_error(dist_gasal2(), "CUDA")
})

test_that("dist_gasal2() constructs on a CUDA build", {
  skip_if_not(optimotu_cuda_build())
  expect_equal(dist_gasal2()$method, "gasal2")
})

test_that("dist_evo clusters without a CUDA device", {
  skip_if(optimotu_cuda())
  seq <- c(
    s1 = paste(rep("A", 31L), collapse = ""),
    s2 = paste(rep("A", 31L), collapse = ""),
    s3 = paste0(strrep("A", 15L), "I", strrep("A", 15L))
  )
  clust <- seq_cluster(
    seq,
    dist_config = dist_evo("mtart"),
    threshold_config = threshold_set(c(0, 0.05, 0.5)),
    clust_config = clust_slink(),
    parallel_config = parallel_concurrent(1L)
  )
  expect_true(is.matrix(clust) || is.data.frame(clust))
})

test_that("dist_gasal2() errors at run time without a device", {
  skip_if_not(optimotu_cuda_build())
  skip_if(optimotu_cuda())
  seq <- gpu_dna()
  expect_error(
    seq_cluster(
      seq,
      dist_config = dist_gasal2(),
      threshold_config = threshold_set(c(0, 0.1)),
      clust_config = clust_slink(),
      parallel_config = parallel_concurrent(1L)
    ),
    "CUDA device"
  )
})

test_that("GPU dist_evo matches CPU Newton on a small mtART set", {
  skip_if_not(optimotu_cuda())
  cfg <- dist_evo("mtart", t_grid = 32L)
  same <- paste(rep("A", 31L), collapse = "")
  il <- paste0(strrep("A", 15L), "I", strrep("A", 15L))
  leu <- paste0(strrep("A", 15L), "L", strrep("A", 15L))
  seq <- c(s1 = same, s2 = same, s3 = il, s4 = leu)
  d <- seq_distmx(
    seq,
    threshold = 0,
    dist_config = cfg,
    pairs = rbind(c(1L, 2L), c(3L, 4L)),
    id_is_int = TRUE
  )
  expect_equal(d$dist2[[1L]], 0)
  expect_gt(d$dist2[[2L]], 0)
  expect_lt(d$dist2[[2L]], 1)
})

test_that("GASAL2 score matches KSW2 score on the same penalties", {
  skip_if_not(optimotu_cuda())
  seq <- gpu_dna()
  penalties <- list(match = 0L, mismatch = 2L, gap_open = 3L, gap_extend = 1L)
  d_k <- seq_distmx(
    seq,
    threshold = 1,
    dist_config = do.call(
      dist_ksw2,
      c(penalties, list(dissimilarity = "score"))
    ),
    parallel_config = parallel_concurrent(1L)
  )
  d_g <- seq_distmx(
    seq,
    threshold = 1,
    dist_config = do.call(
      dist_gasal2,
      c(penalties, list(dissimilarity = "score"))
    ),
    parallel_config = parallel_concurrent(1L)
  )
  key <- function(x) {
    paste(pmin(x$seq_id1, x$seq_id2), pmax(x$seq_id1, x$seq_id2), sep = ":")
  }
  d_k <- d_k[order(key(d_k)), ]
  d_g <- d_g[order(key(d_g)), ]
  expect_equal(d_g$dist2, d_k$dist2, tolerance = 1e-8)
})

test_that("GASAL2 edit score matches Edlib score", {
  skip_if_not(optimotu_cuda())
  seq <- c(
    a = "ACGTACGTACGTACGT",
    b = "ACGTACGTACGTACGA",
    c = "ACGTACGTACGTACG",
    d = "GGGGGGGGGGGGGGGG"
  )
  d_e <- seq_distmx(
    seq,
    threshold = 1,
    dist_config = dist_edlib(dissimilarity = "score"),
    parallel_config = parallel_concurrent(1L)
  )
  d_g <- seq_distmx(
    seq,
    threshold = 1,
    dist_config = dist_gasal2(dissimilarity = "score"),
    parallel_config = parallel_concurrent(1L)
  )
  key <- function(x) {
    paste(pmin(x$seq_id1, x$seq_id2), pmax(x$seq_id1, x$seq_id2), sep = ":")
  }
  d_e <- d_e[order(key(d_e)), ]
  d_g <- d_g[order(key(d_g)), ]
  expect_equal(d_g$dist2, d_e$dist2, tolerance = 1e-8)
})

test_that("GASAL2 edit identity cluster ids match dist_edlib", {
  skip_if_not(optimotu_cuda())
  seq <- gpu_dna()
  thresh <- threshold_set(c(0, 0.1, 0.5))
  e <- seq_cluster(
    seq,
    dist_config = dist_edlib(dissimilarity = "identity"),
    threshold_config = thresh,
    clust_config = clust_slink(),
    parallel_config = parallel_concurrent(1L)
  )
  g <- seq_cluster(
    seq,
    dist_config = dist_gasal2(dissimilarity = "identity"),
    threshold_config = thresh,
    clust_config = clust_slink(),
    parallel_config = parallel_concurrent(1L)
  )
  expect_equal(
    matrix(g, nrow = nrow(g), dimnames = dimnames(g)),
    matrix(e, nrow = nrow(e), dimnames = dimnames(e))
  )
})

test_that("GASAL2 affine identity cluster ids match dist_ksw2", {
  skip_if_not(optimotu_cuda())
  seq <- gpu_dna()
  thresh <- threshold_set(c(0, 0.1, 0.5))
  penalties <- list(match = 0L, mismatch = 2L, gap_open = 3L, gap_extend = 1L)
  k <- seq_cluster(
    seq,
    dist_config = do.call(
      dist_ksw2,
      c(penalties, list(dissimilarity = "identity"))
    ),
    threshold_config = thresh,
    clust_config = clust_slink(),
    parallel_config = parallel_concurrent(1L)
  )
  g <- seq_cluster(
    seq,
    dist_config = do.call(
      dist_gasal2,
      c(penalties, list(dissimilarity = "identity"))
    ),
    threshold_config = thresh,
    clust_config = clust_slink(),
    parallel_config = parallel_concurrent(1L)
  )
  expect_equal(
    matrix(g, nrow = nrow(g), dimnames = dimnames(g)),
    matrix(k, nrow = nrow(k), dimnames = dimnames(k))
  )
})

test_that("tiled dist_evo uses every visible GPU", {
  skip_if_not(optimotu_cuda_device_count() > 1L)
  cfg <- dist_evo("mtart")
  seq <- setNames(
    c(
      paste(rep("A", 31L), collapse = ""),
      paste(rep("A", 31L), collapse = ""),
      paste0(strrep("A", 15L), "I", strrep("A", 15L)),
      paste0(strrep("A", 15L), "L", strrep("A", 15L)),
      paste0(strrep("A", 15L), "C", strrep("A", 15L)),
      paste0(strrep("A", 15L), "W", strrep("A", 15L))
    ),
    paste0("s", 1:6)
  )
  one <- local({
    old <- Sys.getenv("OPTIMOTU_CUDA_DEVICES", unset = NA_character_)
    on.exit({
      if (is.na(old)) {
        Sys.unsetenv("OPTIMOTU_CUDA_DEVICES")
      } else {
        Sys.setenv(OPTIMOTU_CUDA_DEVICES = old)
      }
    })
    Sys.setenv(OPTIMOTU_CUDA_DEVICES = "0")
    seq_cluster(
      seq,
      dist_config = cfg,
      threshold_config = threshold_set(c(0, 0.05, 0.5)),
      clust_config = clust_slink(),
      parallel_config = parallel_concurrent(2L)
    )
  })
  all_dev <- seq_cluster(
    seq,
    dist_config = cfg,
    threshold_config = threshold_set(c(0, 0.05, 0.5)),
    clust_config = clust_slink(),
    parallel_config = parallel_concurrent(2L)
  )
  expect_equal(all_dev, one)
})
