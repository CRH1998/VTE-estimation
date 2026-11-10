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
  
  if (is.null(sigma)) sigma <- kerTests::med_sigma(x,y)
  
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