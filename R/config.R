verify_threshold_steps <- function(thresh_min, thresh_max, thresh_step) {
  if (is.null(thresh_min) || is.null(thresh_max) || is.null(thresh_step)) {
    stop("either 'thresholds' or 'thresh_min', 'thresh_max', and 'thresh_step' must be given")
  }
  if (!is.numeric(thresh_min) || !is.numeric(thresh_max) || !is.numeric(thresh_step)) {
    stop("'thresh_min', 'thresh_max', and 'thresh_step' must all be numbers.")
  }
  if (length(thresh_min) != 1L || length(thresh_max) != 1L || length(thresh_step) != 1L) {
    stop("'thresh_min', 'thresh_max', and 'thresh_step' must all be of length 1.")
  }
  if (is.na(thresh_min) || is.na(thresh_max) || is.na(thresh_step)) {
    stop("'thresh_min', 'thresh_max', and 'thresh_step' may not be NA.")
  }
  if (is.nan(thresh_min) || is.nan(thresh_max) || is.nan(thresh_step)) {
    stop("'thresh_min', 'thresh_max', and 'thresh_step' may not be NaN.")
  }
  if (thresh_step <= 0) {
    stop("'thresh_step' must be positive.")
  }
  if (thresh_max < thresh_min) {
    stop("'thresh_max' must be greater than or equal to 'thresh_min'.")
  }
}

verify_thresholds <- function(thresholds) {
  if (!is.numeric(thresholds)) {
    stop("'thresholds' must be numeric.")
  }
  if (length(thresholds) == 0) {
    stop("At least one threshold must be given in 'thresholds'.")
  }
  if (any(is.na(thresholds))) {
    stop("'thresholds' may not have NA values.")
  }
  if (any(is.nan(thresholds))) {
    stop("'thresholds' may not have NaN values.")
  }
  if (any(thresholds < 0)) {
    stop("'thresholds' may not contain negative values.")
  }
  if (any(utils::head(thresholds, -1) > utils::tail(thresholds, -1))) {
    stop("'thresholds' must be non-decreasing.")
  }
  invisible(TRUE)
}

deduplicate_thresholds <- function(thresholds) {
  out <- list()
  if (all(utils::head(thresholds, -1) < utils::tail(thresholds, -1))) {
    out$thresholds <- thresholds
    out$threshold_order <- seq_along(thresholds)
  } else {
    d <- duplicated(thresholds)
    out$thresholds <- thresholds[!d]
    out$threshold_order <- match(thresholds, thresholds[!d])
  }
  out
}

reduplicate_thresholds <- function(out, thresholds) {
  UseMethod("reduplicate_thresholds", out)
}

#' @method reduplicate_thresholds matrix
reduplicate_thresholds.matrix <- function(out, thresholds) {
  checkmate::assert_class(thresholds, "optimotu_threshold_config")
  if (thresholds$type == "uniform") {
    out
  } else if (isTRUE(all.equal(thresholds$threshold_order, seq_len(nrow(out))))) {
    out
  } else {
    out[thresholds$threshold_order,]
  }
}

#' @method reduplicate_thresholds hclust
reduplicate_thresholds.hclust <- function(out, thresholds) {
  out
}

#' @method reduplicate_thresholds list
reduplicate_thresholds.list <- function(out, thresholds) {
  lapply(out, function(x) reduplicate_thresholds(x, thresholds = thresholds))
}

verify_precision <- function(precision) {
  checkmate::assert_number(
    precision,
    na.ok = FALSE,
    null.ok = TRUE,
    finite = TRUE,
    lower = .Machine$double.xmin
  )
}

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
      fill_method = switch(fill_method, linear = 1L, binary = 2L, topdown = 3L)
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
    set = threshold_set(...),
    lookup = threshold_lookup(...)
  )
}

#' @param from (`numeric` scalar) smallest threshold; typically between 0 and 1
#' @param to (`numeric` scalar) largest threshold; typically between 0 and 1
#' @param by (`numeric` scalar) step size
#' @export
#' @describeIn threshold_config helper function for type `"uniform"`
threshold_uniform <- function(from, to, by, thresh_names = NULL) {
  verify_threshold_steps(from, to, by)
  structure(
    list(type = "uniform", from = from, to = to, by = by, thresh_names = thresh_names),
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
threshold_set <- function(thresholds, thresh_names = names(thresholds)) {
  verify_thresholds(thresholds)
  # out <- if (length(thresholds) == 1) {
  #   list(
  #     type = "uniform",
  #     from = thresholds[1],
  #     to = thresholds[1],
  #     by = 1
  #   )
  # } else if (length(thresholds) == 2) {
  #   list(
  #     type = "uniform",
  #     from = thresholds[1],
  #     to = thresholds[2],
  #     by = thresholds[2] - thresholds[1]
  #   )
  # } else {
  #   g = thresholds[2] - thresholds[1]
  #   for (x in thresholds[-(1:2)] - thresholds[1]) g <- gcd(g, x)
  #   if ((thresholds[length(thresholds)] - thresholds[1]) / g / length(thresholds) > 100) {
  #     list(type = "set", thresholds = thresholds)
  #   } else {
  #     list(type = "lookup", thresholds = thresholds, precision = g)
  #   }
  # }
  structure(
    c(
      list(type = "set"),
      deduplicate_thresholds(thresholds),
      list(thresh_names = thresh_names)
    ),
    class = "optimotu_threshold_config"
  )
}

#' @param precision (`numeric` scalar) precision for distances; this is used to
#' generate a look-up table for distances to the smallest encompassing threshold.
#' @export
#' @describeIn threshold_config helper function for method `"lookup"`
threshold_lookup <- function(thresholds, precision, thresh_names = names(thresholds)) {
  verify_thresholds(thresholds)
  verify_precision(precision)
  structure(
    c(
      list(type = "lookup"),
      deduplicate_thresholds(thresholds),
      list(precision = precision, thresh_names = thresh_names)
    ),
    class = "optimotu_threshold_config"
  )
}

#' Configuration for parallelization options
#'
#' @details
#'
#' # Merge method
#'
#' In the merge method, each thread works on an independent clustering of the
#' data, based on its own (disjoint) subset of the distance matrix. When each
#' thread finishes clustering, it merges its results into the master clustering.
#' This avoids concurrency collisions between the threads during the main
#' clustering, although collisions occur if multiple threads are trying to
#' merge at the same time.
#'
#' # Concurrent method
#'
#' In the concurrent method, all threads work jointly on the same clustering of
#' the data. Although many threads can simultaneously read the clustering to
#' determine whether a new pairwise distance will lead to an update (most do
#' not), only one thread can update the clustering at a time, so this method
#' can lead to more concurrency collisions when many threads are in use.
#' However, sharing the state of the clustering between threads leads to fewer
#' total updates than the merge method.
#'
#' # Hierarchical method
#'
#' The hierarchical method is a combination of several "shards", each of which
#' in turn has multiple threads. Each shard has one clustering object, which the
#' threads within the shard update concurrently. When all threads in the shard
#' have finished, then the results from the shard are merged into the master
#' clustering. This is probably the most efficient method when there are very
#' many threads, but the optimal number of shards varies between data sets.
#'
#' @param method (`character` string) parallelization method.
#' @param threads (positive `integer` scalar) total number of threads to use
#' @param ... passed on to variants
#'
#' @return an object representing the thresholds
#' @export
parallel_config <- function(method = c("merge", "concurrent", "hierarchical"), ...) {
  method = match.arg(method)
  switch(
    method,
    merge = parallel_merge(...),
    concurrent = parallel_concurrent(...),
    hierarchical = parallel_hierarchical(...)
  )
}


#' @export
#' @describeIn parallel_config helper function for method `"merge"`
parallel_merge <- function(threads) {
  checkmate::assert_integerish(threads, lower = 1L)
  structure(
    list(method = "merge", threads = as.integer(threads)),
    class = "optimotu_parallel_config"
  )
}


#' @export
#' @describeIn parallel_config helper function for method `"concurrent"`
parallel_concurrent <- function(threads) {
  checkmate::assert_integerish(threads, lower = 1L)
  structure(
    list(method = "concurrent", threads = as.integer(threads)),
    class = "optimotu_parallel_config"
  )
}

#' @param shards (positive `integer` scalar) number of independent working
#' units for clustering. Must be less than or equal to `threads` (and in order
#' for the result to actually be hierarchical, should be at least 2 and at most
#' `threads`/2)
#' @export
#' @describeIn parallel_config helper function for method `"hierarchical"`
parallel_hierarchical <- function(threads, shards) {
  checkmate::assert_integerish(threads, lower = 1L)
  checkmate::assert_integerish(shards, lower = 1L, upper = threads)
  structure(
    list(
      method = "hierarchical",
      threads = as.integer(threads),
      shards = as.integer(shards)
    ),
    class = "optimotu_parallel_config"
  )
}
