############################################################
#                                                          #
#   The following script includes helper functions to      #
#   construct the design matrix, validate covariates,      #
#   and checking inputs for the nuisance_func.R script.    #
#                                                          #
############################################################


############################################################
#                         Packages                         #
############################################################

# For data wrangling
library(tidyverse)

# For median heuristic
library(kerTests)


############################################################
#                         Dependencies                     #
############################################################


############################################################
#       Functions for constructing the design matrix,      #
#       validating covariates, and checking inputs.        #
############################################################


#-------------------------------------------------------------------------------
# 1. Construct the design matrix for the nuisance function.
# ------------------------------------------------------------------------------
.nf_design <- function(a, x) {
  a <- .nf_treatment(a, nrow(x))
  interactions <- as.matrix(x) * a
  d <- data.frame(A = a, x, interactions, check.names = FALSE)
  names(d) <- make.names(c("A", names(x), paste0("A_", names(x))), unique = TRUE)
  d
}


#-------------------------------------------------------------------------------
# 2. Validate X-covariates and return a data frame
#    with the same columns in the same order as x_names.
#-------------------------------------------------------------------------------
.nf_covariates <- function(x, x_names) {
  if (is.null(dim(x))) {
    if (!is.numeric(x) || length(x) != length(x_names)) {
      stop("A covariate vector must describe one observation with all features.")
    }
    x <- matrix(x, nrow = 1L, dimnames = list(NULL, names(x)))
  }
  if (!(is.matrix(x) || is.data.frame(x)) || nrow(x) == 0L) {
    stop("Use a vector, matrix, or data frame with at least one observation.")
  }
  if (is.null(colnames(x))) {
    if (ncol(x) != length(x_names)) stop("Wrong number of unnamed covariates.")
    colnames(x) <- x_names
  }
  if (anyDuplicated(colnames(x)) || !all(x_names %in% colnames(x))) {
    stop("Named inputs must contain each training covariate exactly once.")
  }
  x <- as.data.frame(x, optional = TRUE)[x_names]
  ok <- vapply(x, function(v) is.numeric(v) && is.null(dim(v)) && all(is.finite(v)), logical(1))
  if (!all(ok)) stop("Use finite numeric covariates. Encode categorical features first.")
  x
}


#-------------------------------------------------------------------------------
# 3. Validate that a count is a whole number between lower and upper.
#-------------------------------------------------------------------------------
.nf_count <- function(x, lower, upper, label) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x != floor(x) || x < lower || x > upper) {
    stop(sprintf("%s must be a whole number between %d and %d.", label, lower, upper))
  }
}


#-------------------------------------------------------------------------------
# 4. Validate that a probability is a number between 0 and 1.
#-------------------------------------------------------------------------------
.nf_treatment <- function(a, n) {
  if (!(is.numeric(a) || is.logical(a)) || !(length(a) %in% c(1L, n)) ||
      !all(a %in% 0:1)) {
    stop("Treatment must be 0 or 1, with length one or one value per observation.")
  }
  rep_len(as.numeric(a), n)
}


#-------------------------------------------------------------------------------
# 5. Validate that a vector of values is numeric, finite, and has length n.
#-------------------------------------------------------------------------------
.nf_values <- function(v, n, label) {
  if (!(is.numeric(v) || is.logical(v)) || length(v) != n || any(!is.finite(v))) {
    stop(paste(label, "must supply one finite number per observation."))
  }
  as.numeric(v)
}


#-------------------------------------------------------------------------------
# 6. Construct the empirical measure of a sample of observations.
#-------------------------------------------------------------------------------
# A vector is a sample of scalar observations. Each matrix or data-frame row
# is one observation. Pass vectorized = TRUE to evaluate f on the whole sample.
.nf_empirical_measure <- function(X) {
  n <- NROW(X)
  if (n == 0L) stop("X must contain at least one observation.")
  
  if (is.matrix(X)) {
    observation <- function(i) X[i, , drop = TRUE]
  } else if (is.data.frame(X)) {
    observation <- function(i) X[i, , drop = FALSE]
  } else if (is.null(dim(X)) && (is.atomic(X) || is.list(X))) {
    observation <- function(i) X[[i]]
  } else {
    stop("Supply a vector, a list of observations, a matrix, or a data frame.")
  }
  
  function(f, vectorized = FALSE) {
    if (!is.function(f)) stop("f must be a function.")
    values <- if (vectorized) {
      f(X)
    } else {
      vapply(seq_len(n), function(i) f(observation(i)), numeric(1))
    }
    if (length(values) != n || !(is.numeric(values) || is.logical(values))) {
      stop("f must give one numeric or logical value per observation.")
    }
    mean(values)
  }
}


