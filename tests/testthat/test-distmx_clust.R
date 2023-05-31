set.seed(1)
k <- 5
m <- 5
n <- k*m
# generate centroids for 5 clusters in 2d space
centroids = data.frame(
  x = runif(k),
  y = runif(k)
)
# for each cluster, generate 5 random points around it
points = lapply(centroids, rep, each = m) |>
  lapply(rnorm, n = n, sd = 0.1) |>
  as.data.frame()

# euclidian distance matrix aroudn the points
distmx <- dist(points)

# convert distance matrix to "sparse" format
dist_table <- data.frame(
  seq1 = rep(0:(n-2), (n-1):1),
  seq2 = unlist(lapply(1:(n-1), seq, to = (n-1))),
  dist = unclass(distmx)
)
dist_table <- dist_table[dist_table$dist <= 1,]
# write as a temp file
distmx_file <- tempfile(fileext = ".distmx")
write.table(dist_table, distmx_file, sep = "\t", row.names = FALSE, col.names = FALSE)
withr::defer(unlink(distmx_file), teardown_env())

hclust_matrix <- hclust(distmx, method = "single") |>
  lapply(0:20 / 20, cutree, k = NULL, tree = _) |>
  lapply(function(x) {
    for (i in seq_along(x)) {if (x[i] < i) x[x >= i] = x[x >= i] + 1L}
    x
  }) |>
  do.call(rbind, args = _)
hclust_matrix <- hclust_matrix - 1L

uniform_thresh <- threshold_uniform(0, 1, 0.05)
set_thresh <- threshold_set(0:20 / 20)
cache_thresh <- threshold_cached(0:20/20, 0.05)

test_that("distmx_cluster index method agrees with hclust", {
  expect_equal(
    optimotu:::distmx_cluster_single(
      file = distmx_file,
      seqnames = as.character(1:n),
      threshold_config = uniform_thresh,
      method_config = clust_index(),
      parallel_config = parallel_concurrent(1),
      output_type = "matrix"
    ),
    hclust_matrix
  )
})

test_that("distmx_cluster tree method agrees with hclust", {
  expect_equal(
    optimotu:::distmx_cluster_single(
      file = distmx_file,
      seqnames = as.character(1:n),
      threshold_config = uniform_thresh,
      method_config = clust_tree(),
      parallel_config = parallel_concurrent(1),
      output_type = "matrix"
    ),
    hclust_matrix
  )
})

test_that("distmx_cluster matrix (binary, binary) method agrees with hclust", {
  expect_equal(
    optimotu:::distmx_cluster_single(
      file = distmx_file,
      seqnames = as.character(1:n),
      threshold_config = uniform_thresh,
      method_config = clust_matrix(TRUE, "binary"),
      parallel_config = parallel_concurrent(1),
      output_type = "matrix"
    ),
    hclust_matrix
  )
})

test_that("distmx_cluster matrix (binary, linear) method agrees with hclust", {
  expect_equal(
    optimotu:::distmx_cluster_single(
      file = distmx_file,
      seqnames = as.character(1:n),
      threshold_config = uniform_thresh,
      method_config = clust_matrix(TRUE, "linear"),
      parallel_config = parallel_concurrent(1),
      output_type = "matrix"
    ),
    hclust_matrix
  )
})

test_that("distmx_cluster matrix (binary, topdown) method agrees with hclust", {
  expect_equal(
    optimotu:::distmx_cluster_single(
      file = distmx_file,
      seqnames = as.character(1:n),
      threshold_config = uniform_thresh,
      method_config = clust_matrix(TRUE, "topdown"),
      parallel_config = parallel_concurrent(1),
      output_type = "matrix"
    ),
    hclust_matrix
  )
})

test_that("distmx_cluster matrix (linear, binary) method agrees with hclust", {
  expect_equal(
    optimotu:::distmx_cluster_single(
      file = distmx_file,
      seqnames = as.character(1:n),
      threshold_config = uniform_thresh,
      method_config = clust_matrix(FALSE, "binary"),
      parallel_config = parallel_concurrent(1),
      output_type = "matrix"
    ),
    hclust_matrix
  )
})

test_that("distmx_cluster matrix (linear, linear) method agrees with hclust", {
  expect_equal(
    optimotu:::distmx_cluster_single(
      file = distmx_file,
      seqnames = as.character(1:n),
      threshold_config = uniform_thresh,
      method_config = clust_matrix(FALSE, "linear"),
      parallel_config = parallel_concurrent(1),
      output_type = "matrix"
    ),
    hclust_matrix
  )
})

test_that("distmx_cluster matrix (linear, topdown) method agrees with hclust", {
  expect_equal(
    optimotu:::distmx_cluster_single(
      file = distmx_file,
      seqnames = as.character(1:n),
      threshold_config = uniform_thresh,
      method_config = clust_matrix(FALSE, "topdown"),
      parallel_config = parallel_concurrent(1),
      output_type = "matrix"
    ),
    hclust_matrix
  )
})

