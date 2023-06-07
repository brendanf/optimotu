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

# take some subsets for testing subset clustering
subsets <- lapply(5:20, sample.int, n = n)

# euclidian distance matrix around the points, rounded to precision of 0.05
distmx <- round(dist(points) / 0.05) * 0.05
subsetdistmx <- lapply(subsets, \(i) round(dist(as.matrix(points)[i,]) / 0.05)*0.05)


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
# withr::defer(unlink(distmx_file), teardown_env())

hclust2matrix <- function(distmx, thresholds) {
  out <- hclust(pmax(distmx - 1e-8, 0), method = "single") |>
    lapply(thresholds, cutree, k = NULL, tree = _) |>
    lapply(function(x) {
      for (i in seq_along(x)) {if (x[i] < i) x[x >= i] = x[x >= i] + 1L}
      x
    }) |>
    do.call(rbind, args = _)
  out - 1L
}

hclust_matrix <- hclust2matrix(distmx, 0:20 / 20)
subset_hclust_matrix <- lapply(subsetdistmx, hclust2matrix, 0:20 / 20)

thresholds <- list(
  set_thresh = threshold_set(0:20 / 20),
  uniform_thresh = threshold_uniform(0, 1, 0.05),
  cache_thresh = threshold_cached(0:20 / 20, 0.05)
)

algorithms <- list(
  tree = clust_tree(),
  matrix_binary_binary = clust_matrix(binary_search = TRUE, fill_method = "binary"),
  matrix_binary_linear = clust_matrix(binary_search = TRUE, fill_method = "linear"),
  matrix_binary_topdown = clust_matrix(binary_search = TRUE, fill_method = "topdown"),
  matrix_linear_binary = clust_matrix(binary_search = FALSE, fill_method = "binary"),
  matrix_linear_linear = clust_matrix(binary_search = FALSE, fill_method = "linear"),
  matrix_linear_topdown = clust_matrix(binary_search = FALSE, fill_method = "topdown"),
  index = clust_index()
)

parallels <- list(
  merge = parallel_merge(4),
  serial = parallel_concurrent(1),
  concurrent = parallel_concurrent(4),
  hierarchical = parallel_hierarchical(2, 2)
)

for (p in names(parallels)) {
  for (t in names(thresholds)) {
    for (a in names(algorithms)) {
      # cat(sprintf(
      #   "%s distmx_cluster_single %s method with %s thresholds agrees with hclust\n",
      #   p, a, t
      # ))
      test_that(
        sprintf(
          "%s distmx_cluster_single %s method with %s thresholds agrees with hclust",
          p, a, t
        ),
        {
          expect_equal(
            distmx_cluster(
              distmx = distmx_file,
              names = as.character(1:n),
              threshold_config = thresholds[[t]],
              clust_config = algorithms[[a]],
              parallel_config = parallels[[p]],
              output_type = "matrix"
            ),
            hclust_matrix
          )
        }
      )
      # cat(
      #   sprintf(
      #     "%s distmx_cluster_multi %s method with %s thresholds agrees with hclust\n",
      #     p, a, t
      #   )
      # )
      test_that(
        sprintf(
          "%s distmx_cluster_multi %s method with %s thresholds agrees with hclust",
          p, a, t
        ),
        {
          expect_equal(
            distmx_cluster(
              distmx = distmx_file,
              names = as.character(1:n),
              which = lapply(subsets, as.character),
              threshold_config = thresholds[[t]],
              clust_config = algorithms[[a]],
              parallel_config = parallels[[p]],
              output_type = "matrix"
            ),
            subset_hclust_matrix
          )
        }
      )
    }
  }
}

