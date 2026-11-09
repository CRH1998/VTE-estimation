############################################################
#                                                          #
#   The following script includes kernel functions         #
#                                                          #
############################################################


############################################################

# For kernels
library(kerntools)

############################################################


#' Title guassian_kernel
#'
#' @param x Data point 1 (vector or matrix)
#' @param y Data point 2
#' @param sigma Bandwidth parameter for the Gaussian kernel
#'
#' @returns The Gaussian kernel value between the two data points
#' @export
#'
#' @examples
#' 


gaussian_kernel <- function(x, y, sigma = 1) {
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




#' Title
#'
#' @param train_data Training data used to define the kernel function. This should be a matrix or data frame where each row is a data point and each column is a feature.
#' @param sigma Bandwidth parameter for the kernel function. This controls the width of the kernel and affects how much influence each training point has on the kernel evaluation. Default is 1.
#' @param kernel A kernel function that takes two arguments (a data point and the training data) and returns a vector of kernel evaluations. The default is the Gaussian kernel defined above.
#'
#' @returns
#' @export
#'
#' @examples
data_dependent_kernel <- function(train_data, sigma = 1, kernel = gaussian_kernel) {
  
  ddk <- function(x) kernel(x, train_data, sigma)
  
  return(ddk)
}



#' Title
#'
#' @param k kernel function that takes two arguments (a data point and the training data) and returns a vector of kernel evaluations. This function should be defined such that it can handle both single data points and matrices of data points.
#' @param f_hat A function that takes a matrix or data frame of training data and returns a numeric vector of predictions (one for each row of the training data). This function represents the estimated function that we want to embed using the kernel.
#' @param train_data A matrix or data frame of training data where each row is a data point and each column is a feature. This data is used to compute the kernel embedding.
#'
#' @returns
#' @export
#'
#' @examples
make_kernel_embedding <- function(k, f_hat, train_data) {
  stopifnot(is.function(k), is.function(f_hat))
  if (is.null(dim(train_data))) train_data <- matrix(train_data, nrow = 1L)
  if (!(is.matrix(train_data) || is.data.frame(train_data)) || nrow(train_data) == 0L) {
    stop("train_data must be a nonempty sample with one observation per row.")
  }
  n <- nrow(train_data)
  
  f_X <- f_hat(train_data)
  if (!is.numeric(f_X) || length(f_X) != n) {
    stop("f_hat(train_data) must return one numeric prediction per row of train_data.")
  }
  f_X <- as.numeric(f_X)
  
  function(x) {
    k_x <- k(x, train_data)
    if (!is.numeric(k_x) || length(k_x) != n) {
      stop("For a single x, k(x, train_data) must return one numeric value per row of train_data.")
    }
    mean(as.numeric(k_x) * f_X)
  }
}


