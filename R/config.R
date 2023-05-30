#' Configuration for clustering algorithms
#'
#' @details
#' # Tree method
#'
#' The tree method maintains the state of the clustering algorithm using a
#' tree structure. When there are more than a few thresholds, it uses less
#' memory than the other methods, and requires fewer operations to update.
#' However, because its data structure is less cache friendly and it requires
#' multiple checks to determine whether an incoming pairwise distance will
#' lead to a cluster update, it is often slower than the matrix-based methods.
#' It is the only data structure whose result can currently be returned as
#' an `hclust` object.
#'
#' # Matrix method
#'
#' The matrix method maintains the state of the clustering algorithm in a
#' contiguous matrix structure; for matrix output and when using a single thread
#' or concurrent parallelism, this does not require any additional memory beyond
#' what is required to store the output. The matrix method uses the least code,
#' and despite the higher number of operations required for updates, it is
#' competitive due to its \eqn{\mathcal{O}(1)}{O(1)} checks on incoming pairwise
#' distances and good cache locality.
#'
#' # Index(ed matrix) method
#'
#' The indexed matrix method uses a linked-list index to order the columns of
#' the matrix methods such that clusters are contiguous.  This adds some
#' overhead to updates in order to maintain the index, but drastically reduces
#' the number of matrix columns which must be accessed during an update.
#'
#' @param method (`character` string) The clustering algorithm to use. Options
#' are `"tree"`, `"matrix"`, and `"index"`
#' @param ... passed on to variants
#'
#' @return an object describing the clustering algorithm, to pass to
#' `distmx_cluster()` or `seq_cluster()`
#' @export
clust_config <- function(method = c("tree", "matrix", "index"), ...) {
  method = match.arg(method)
  switch(
    method,
    tree = clust_tree(...),
    matrix = clust_matrix(...),
    index = clust_index(...)
  )
}


#' @export
#' @describeIn clust_config helper function for method `"tree"`
clust_tree <- function() {
  structure(
    list(method = "tree"),
    class = "optimotu_cluster_config"
  )
}

#' @param binary_search (`logical` flag) if `TRUE`, use binary search instead
#' of linear search when determining the currently known minimum distance
#' between two sequences.  This may be slightly faster when the number of
#' thresholds is very large.
#' @param fill_method (`character` string) method to use to determine which
#' matrix elements must be updated for a sequence.  `"linear"` and `"binary"`
#' both update the matrix column with a range memory write (using `memcpy`)
#' after using a search (either linear or binary) to determine which range needs
#' to be updated.  `"topdown"` fills element-by element starting at the largest
#' distance, and is primarily included only as a testing option, since it is
#' usually slower.
#' @export
#' @describeIn clust_config helper function for method `"matrix"`
clust_matrix <- function(
    binary_search = TRUE,
    fill_method = c("binary", "linear", "topdown")
) {
  fill_method = match.arg(fill_method)
  structure(
    list(
      method = "matrix",
      binary_search = binary_search,
      fill_method = switch(fill_method(linear = 1L, binary = 2L, topdown = 2L))
    ),
    class = "optimotu_cluster_config"
  )
}

#' @export
#' @describeIn clust_config helper function for method `"index"`
clust_index <- function() {
  structure(
    list(method = "index"),
    class = "optimotu_cluster_config"
  )
}

#' Configuration for threshold representations
#'
#' @param type (`character` string) representation of thresholds.
#' @param ... passed on to variants
#'
#' @return an object representing the thresholds
#' @export
threshold_config <- function(type = c("uniform", "set", "lookup"), ...) {
  type = match.arg(type)
  switch(
    type,
    uniform = threshold_uniform(...),
    set = threshold_array(...),
    lookup = threshold_cached(...)
  )
}

#' @param from (`numeric` scalar) smallest threshold; typically between 0 and 1
#' @param to (`numeric` scalar) largest threshold; typically between 0 and 1
#' @param by (`numeric` scalar) step size
#' @export
#' @describeIn threshold_config helper function for type `"uniform"`
threshold_uniform <- function(from, to, by) {
  verify_threshold_steps(from, to, by)
  structure(
    list(type = "uniform", from = from, to = to, by = by),
    class = "optimotu_threshold_config"
  )
}

gcd <- function(a, b) {
  while (abs(b) > sqrt(.Machine$double.eps)) {
    c <- a %% b
    a <- b
    b <- c
  }
  a
}

#' @param thresholds (`numeric` vector) explicit list of thresholds.  Should be
#' sorted in ascending order.
#' @export
#' @describeIn threshold_config helper function for method `"set"`
threshold_set <- function(thresholds) {
  verify_thresholds(thresholds)
  out <- if (length(thresholds) == 1) {
    list(
      type = "uniform",
      from = thresholds[1],
      to = thresholds[1],
      by = 1
    )
  } else if (length(thresholds) == 2) {
    list(
      type = "uniform",
      from = thresholds[1],
      to = thresholds[2],
      by = thresholds[2] - thresholds[1]
    )
  } else {
    g = thresholds[2] - thresholds[1]
    for (x in thresholds[-(1:2)] - thresholds[1]) g <- gcd(g, x)
    if ((thresholds[length(thresholds)] - thresholds[1]) / g / length(thresholds) > 100) {
      list(type = "set", thresholds = thresholds)
    } else {
      list(type = "lookup", thresholds = thresholds, precision = g)
    }
  }
  structure(out, class = "optimotu_threshold_config")
}

#' @param precision (`numeric` scalar) precision for distances; this is used to
#' generate a look-up table for distances to the smallest encompassing threshold.
#' @export
#' @describeIn threshold_config helper function for method `"lookup"`
threshold_cached <- function(thresholds, precision) {
  verify_thresholds(thresholds)
  verify_precision(precision)
  structure(
    list(type = "lookup", thresholds = thresholds, precision = precision),
    class = "optimotu_threshold_config"
  )
}
