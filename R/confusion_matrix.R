#' @export
#'
#' @param x (either a `data.frame` as returned by `confusion_matrix()`, or an
#' `m` x `n` matrix, as described for `k`) data for which to calculate the
#' index.
#' @param ... passed on to methods
#' @rdname confusion_matrix
rand_index <- function(x, ...) {
  UseMethod("rand_index", x)
}

#' @export
#' @rdname confusion_matrix
rand_index.matrix <- function(x, c, threads = 1, ...) {
  cm <- confusion_matrix(x, c, threads)
  rand_index.data.frame(cm)
}

#' @export
#' @rdname confusion_matrix
#' @method rand_index data.frame
rand_index.data.frame <- function(x, ...) {
  with(x, (TP + TN) / (TP + FP + FN + TN))
}

#' @export
#' @rdname confusion_matrix
adjusted_rand_index <- function(x, ...) {
  UseMethod("adjusted_rand_index", x)
}

#' @export
#' @rdname confusion_matrix
adjusted_rand_index.matrix <- function(x, c, threads = 1, ...) {
  cm <- confusion_matrix(x, c, threads)
  adjusted_rand_index.data.frame(cm)
}

#' @export
#' @rdname confusion_matrix
#' @method adjusted_rand_index data.frame
adjusted_rand_index.data.frame <- function(x, ...) {
  with(x, 2*(TP*TN - FP*FN) / ((TP + FP)*(FP + TN) + (TP + FN)*(FN + TN)))
}

#' @export
#' @rdname confusion_matrix
matthews_correlation_coefficient <- function(x, ...) {
  UseMethod("matthews_correlation_coefficient", x)
}

#' @export
#' @rdname confusion_matrix
matthews_correlation_coefficient.matrix <- function(x, c, threads = 1, ...) {
  cm <- confusion_matrix(x, c, threads)
  matthews_correlation_coefficient.data.frame(cm)
}

#' @export
#' @rdname confusion_matrix
#' @method matthews_correlation_coefficient data.frame
matthews_correlation_coefficient.data.frame <- function(x, ...) {
  with(x, (TP*TN - FP*FN) / sqrt(TP + FP) / sqrt(TP + FN) / sqrt(FN + TN) / sqrt(FP + TN))
}

#' @export
#' @rdname confusion_matrix
fowlkes_mallow_index <- function(x, ...) {
  UseMethod("fowlkes_mallow_index", x)
}

#' @export
#' @rdname confusion_matrix
fowlkes_mallow_index.matrix <- function(x, c, threads = 1, ...) {
  cm <- confusion_matrix(x, c, threads)
  fowlkes_mallow_index.data.frame(cm)
}

#' @export
#' @rdname confusion_matrix
#' @method fowlkes_mallow_index data.frame
fowlkes_mallow_index.data.frame <- function(x, ...) {
  with(x, TP / sqrt(TP + FP) / sqrt(TP + FN) )
}
