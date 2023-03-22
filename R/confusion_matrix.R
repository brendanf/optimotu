#' @export
rand_index <- function(x, ...) {
  UseMethod("rand_index", x)
}

#' @export
rand_index.matrix <- function(x, y, ncpu = 1, ...) {
  cm <- confusion_matrix2(x, y, ncpu)
  rand_index.data.frame(cm)
}

#' @export
#' @method rand_index data.frame
rand_index.data.frame <- function(x, ...) {
  with(x, (TP + TN) / (TP + FP + FN + TN))
}

#' @export
adjusted_rand_index <- function(x, ...) {
  UseMethod("adjusted_rand_index", x)
}

#' @export
adjusted_rand_index.matrix <- function(x, y, ncpu = 1, ...) {
  cm <- confusion_matrix2(x, y, ncpu)
  adjusted_rand_index.data.frame(cm)
}

#' @export
#' @method adjusted_rand_index data.frame
adjusted_rand_index.data.frame <- function(x, ...) {
  with(x, 2*(TP*TN - FP*FN) / ((TP + FP)*(FP + TN) + (TP + FN)*(FN + TN)))
}

#' @export
matthews_correlation_coefficient <- function(x, ...) {
  UseMethod("matthews_correlation_coefficient", x)
}

#' @export
matthews_correlation_coefficient.matrix <- function(x, y, ncpu = 1, ...) {
  cm <- confusion_matrix2(x, y, ncpu)
  matthews_correlation_coefficient.data.frame(cm)
}

#' @export
#' @method matthews_correlation_coefficient data.frame
matthews_correlation_coefficient.data.frame <- function(x, ...) {
  with(x, (TP*TN - FP*FN) / sqrt(TP + FP) / sqrt(TP + FN) / sqrt(FN + TN) / sqrt(FP + TN))
}

#' @export
fowlkes_mallow_index <- function(x, ...) {
  UseMethod("fowlkes_mallow_index", x)
}

#' @export
fowlkes_mallow_index.matrix <- function(x, y, ncpu = 1, ...) {
  cm <- confusion_matrix2(x, y, ncpu)
  fowlkes_mallow_index.data.frame(cm)
}

#' @export
#' @method fowlkes_mallow_index data.frame
fowlkes_mallow_index.data.frame <- function(x, ...) {
  with(x, TP / sqrt(TP + FP) / sqrt(TP + FN) )
}
