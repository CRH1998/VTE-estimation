############################################################
#                                                          #
#   The following script includes kernel functions         #
#                                                          #
############################################################


############################################################

# For kernels
library(kerntools)

############################################################

# Private copies of the earlier helpers make this script self-contained.
# A vector is one observation. Matrices and data frames have one observation
# per row, with features in the same order in both inputs.
#' Title
#'
#' @param x Data frame or matrix of covariates (one observation per row).
#' @param y Data frame or matrix of covariates (one observation per row).
#' @param sigma Bandwidth parameter for the Gaussian kernel. If NULL, it is estimated using the median heuristic.
#'
#' @returns A matrix of kernel values, where the (i,j)-th entry is the Gaussian kernel evaluated at the i-th row of x and the j-th row of y.
#' @export
#'
#' @examples
.nf_gaussian_kernel <- function(x, y, sigma = NULL) {
  
  
  as_rows <- function(z) {
    if (is.data.frame(z)) z <- as.matrix(z)
    if (is.null(dim(z))) z <- matrix(z, nrow = 1L)
    if (!is.matrix(z) || !is.numeric(z)) {
      stop("Inputs must be numeric vectors, matrices, or data frames.")
    }
    storage.mode(z) <- "double"
    z
  }
  
  x <- as_rows(x)
  y <- as_rows(y)
  
  if (is.null(sigma)) {
    print("Fitting median heuristic bandwidth using kerTests::med_sigma()")
    sigma <- kerTests::med_sigma(x,y)
  }
  
  if (ncol(x) != ncol(y)) {
    stop("x and y must have the same number of features (columns).")
  }
  if (!is.numeric(sigma) || length(sigma) != 1L ||
      !is.finite(sigma) || sigma <= 0) {
    stop("sigma must be one positive, finite number.")
  }
  sigma <- as.numeric(sigma)
  
  # Accumulate squared distances for every pair of rows, one feature at a time.
  d2 <- matrix(0, nrow = nrow(x), ncol = nrow(y))
  for (j in seq_len(ncol(x))) {
    delta <- outer(x[, j], y[, j], "-") / sigma
    d2 <- d2 + delta^2
  }
  drop(exp(-0.5 * d2))
}




#' Estimate a Gaussian bandwidth from a reference sample of covariates
#'
#' Uses the square root of the median squared Euclidean distance over unordered
#' pairs of different row indices, including zero distances from duplicate rows.
#' With all pairs this matches kerTests::med_sigma(X, Y) for reference = rbind(X, Y).
#' @param reference Numeric matrix or data frame with one observation per row.
#'   A numeric vector denotes a sample of scalar observations.
#' @param max_pairs Inf for the exact calculation. If fewer than choose(n, 2)
#'   pairs are requested, sample pairs uniformly with replacement. Use set.seed()
#'   before this call for a reproducible approximation.
#' @return A positive, finite bandwidth.
#' @export
#' @examples
#' ref <- data.frame(x = c(0, 1, 3, 7))
#' stopifnot(isTRUE(all.equal(median_gaussian_bandwidth(ref), sqrt(12.5))))
#' if (requireNamespace("kerTests", quietly = TRUE)) {
#'   package_sigma <- kerTests::med_sigma(ref[1:2, , drop = FALSE],
#'                                       ref[3:4, , drop = FALSE])
#'   stopifnot(isTRUE(all.equal(median_gaussian_bandwidth(ref), package_sigma)))
#' }
median_gaussian_bandwidth <- function(reference, max_pairs = Inf) {
  if (!(is.matrix(reference) || is.data.frame(reference) ||
        (is.numeric(reference) && is.null(dim(reference))))) {
    stop("Supply a numeric reference sample with observations in rows.")
  }
  X <- as.matrix(reference)
  n <- nrow(X)
  if (!is.numeric(X) || n < 2L || ncol(X) < 1L || any(!is.finite(X))) {
    stop("The reference sample needs at least two rows and finite numeric covariates.")
  }
  if (!is.numeric(max_pairs) || length(max_pairs) != 1L || is.na(max_pairs) ||
      max_pairs < 1 || max_pairs != floor(max_pairs)) {
    stop("max_pairs must be a positive integer or Inf.")
  }
  if (max_pairs >= choose(n, 2L)) {
    d2 <- as.numeric(stats::dist(X))^2
  } else {
    # Draw distinct indices within each pair. Both orientations are equiprobable.
    i <- sample.int(n, max_pairs, replace = TRUE)
    j <- sample.int(n - 1L, max_pairs, replace = TRUE)
    j <- j + (j >= i)
    storage.mode(X) <- "double"
    d2 <- numeric(max_pairs)
    for (p in seq_len(ncol(X))) {
      delta <- X[i, p] - X[j, p]
      d2 <- d2 + delta^2
    }
  }
  sigma <- sqrt(stats::median(d2))
  if (!is.finite(sigma)) stop("The bandwidth is non-finite. Check the covariate scales.")
  if (sigma <= 0) {
    stop("The median pairwise distance is zero. Check for duplicate rows or supply a positive bandwidth.")
  }
  sigma
}

#' Fit a Gaussian kernel once for repeated evaluation
#'
#' @inheritParams median_gaussian_bandwidth
#' @return A two-argument Gaussian kernel with its fitted bandwidth in attr(k, "sigma").
#' @export
#' @examples
#' ref <- data.frame(x = c(0, 1, 3, 7))
#' k <- make_median_gaussian_kernel(ref)
#' stopifnot(isTRUE(all.equal(k(0, 1), exp(-1 / 25))))
#' stopifnot(isTRUE(all.equal(k(0, 3), k(ref[1:2, , drop = FALSE],
#'                                     ref[3:4, , drop = FALSE])[1, 1])))
make_median_gaussian_kernel <- function(reference, max_pairs = Inf) {
  sigma <- median_gaussian_bandwidth(reference, max_pairs = max_pairs)
  rm(reference)  # Only the fitted bandwidth is needed during evaluation.
  k <- function(x, y) .nf_gaussian_kernel(x, y, sigma = sigma)
  attr(k, "sigma") <- sigma
  k
}
