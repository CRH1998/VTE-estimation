############################################################
#                                                          #
#   The following script includes estimation functions     #
#   for the test of heterogenous effects simulation study  #
#                                                          #
############################################################


############################################################

# For estimation
library(glmnet)
library(SuperLearner)
library(grf)
library(hal9001)

# For data wrangling
library(tidyverse)


############################################################


############################################################
#       Function for estimating nuisance parameters        #
############################################################


#' Cross-fitted nuisance estimates for the embedded centred CATE
#'
#' Supply only the training data of the current outer split. In inner split l,
#' regressions use rows with fold_id != l and the empirical covariate marginal
#' uses rows with fold_id == l. Each element of folds has its own fitted models.
#'
#' @param train_data Data frame with numeric Y, binary A (0/1), and finite
#'   numeric covariates. All columns except Y and A are treated as covariates.
#' @param splits Number of marginal folds, at least two.
#' @param SL.lib SuperLearner library for the outcome and the unknown propensity.
#' @param family gaussian() for a continuous Y or binomial() for binary Y.
#' @param hal Use hal9001 for the outcome regression. An unknown propensity
#'   is fitted with SuperLearner in either case.
#' @param randomized Use the known constant treatment probability pi_hat.
#' @param pi_hat Known P(A = 1) when randomized is TRUE.
#' @param K Number of internal CV folds for SuperLearner and, if used, HAL.
#' @param kernel Function k(U, x) returning one numeric value per row of U
#'   when x is a one-row data frame. The default is Gaussian with sigma = 1.
#' @param fold_id Optional integer labels from 1 to splits, one per training
#'   row. The default is a balanced random partition independent of the data.
#' @param SL.env Optional environment for custom SuperLearner wrappers.
#'   Defaults to the SuperLearner namespace for its built-in wrappers.
#' @return A list with folds, fold_id, x_names, weights, and average(f).
#'   Each fold contains its row indices, P_X, fitted models, functions mu(a,x),
#'   pi(x), cate(x), centred_cate(x), embedded_cate(x), residual(a,y,x),
#'   and the scalars ate, mean_embedding, and embedding_norm_sq.
#'   A covariate vector is one observation. Rows of a matrix or data frame
#'   are observations. Named inputs are matched to the training covariates.
#'   The two double integrals use the product empirical measure, including
#'   diagonal pairs. average(f) weights f(fold) by the marginal sample size.
#' @export
#' @examples
#' \dontrun{
#' set.seed(12)
#' n <- 120
#' train <- data.frame(A = rbinom(n, 1, 0.5), x1 = rnorm(n), x2 = rnorm(n))
#' train$Y <- with(train, x2 + A * (1 + x1)) + rnorm(n)
#' fit <- nuisance_func(train, SL.lib = "SL.glm", K = 3)
#' x <- data.frame(x1 = c(0, 1), x2 = c(0.5, -0.5))
#' fit$folds[[1]]$cate(x)
#' fit$average(function(s) s$embedded_cate(x))
#'
#' # Check the two product-measure integrals for the default Gaussian kernel.
#' s <- fit$folds[[1]]
#' U <- train[s$marginal_rows, fit$x_names, drop = FALSE]
#' centred <- s$centred_cate(U)
#' gram <- exp(-as.matrix(stats::dist(U))^2 / 2)
#' stopifnot(abs(mean(centred)) < 1e-10,
#'           isTRUE(all.equal(s$mean_embedding, mean(sweep(gram, 1L, centred, "*")))),
#'           isTRUE(all.equal(s$embedding_norm_sq, mean(gram * tcrossprod(centred)))))
#' }
nuisance_func <- function(train_data, splits = 2L,
                          SL.lib = c("SL.glmnet", "SL.randomForest", "SL.xgboost"),
                          family = stats::gaussian(), hal = FALSE,
                          randomized = TRUE, pi_hat = 0.5, K = 5L,
                          kernel = .nf_gaussian_kernel, fold_id = NULL,
                          SL.env = NULL) {
  
  #-----------------------------------------------------------------------------
  # 1.0 Validate inputs
  #-----------------------------------------------------------------------------
  
  # Validate train_data
  if (!is.data.frame(train_data) || !all(c("Y", "A") %in% names(train_data)) ||
      anyNA(names(train_data)) || any(!nzchar(names(train_data))) ||
      anyDuplicated(names(train_data))) {
    stop("train_data needs uniquely named columns including Y and A.")
  }
  
  
  train_data <- as.data.frame(train_data)
  x_names <- setdiff(names(train_data), c("Y", "A"))
  if (!length(x_names)) stop("Supply at least one covariate.")
  
  
  # Validate covariates, treatment, and outcome
  X <- .nf_covariates(train_data, x_names)
  n <- nrow(X)
  A <- .nf_treatment(train_data$A, n)
  Y <- .nf_values(train_data$Y, n, "Y")
  
  
  #-----------------------------------------------------------------------------
  # 1.1 Validate other inputs
  #-----------------------------------------------------------------------------
  
  # Validate hal and randomized
  if (!is.logical(hal) || length(hal) != 1L || is.na(hal) ||
      !is.logical(randomized) || length(randomized) != 1L || is.na(randomized)) {
    stop("hal and randomized must each be TRUE or FALSE.")
  }
  
  # Validate family
  if (!inherits(family, "family") ||
      !(family$family %in% c("gaussian", "binomial")) ||
      !identical(family$link, if (family$family == "gaussian") "identity" else "logit")) {
    stop("Use gaussian() or binomial() for the outcome family.")
  }
  
  if (family$family == "binomial" && !all(Y %in% 0:1)) {
    stop("For a binary outcome, code Y as 0 or 1.")
  }
  
  # Validate pi_hat when randomized
  if (randomized && (!is.numeric(pi_hat) || length(pi_hat) != 1L ||
                     !is.finite(pi_hat) || pi_hat <= 0 || pi_hat >= 1)) {
    stop("The known treatment probability must be strictly between 0 and 1.")
  }
  
  # Validate splits and K
  if (!is.function(kernel)) stop("kernel must be a function.")
  .nf_count(splits, 2L, n, "splits")
  .nf_count(K, if (hal) 3L else 2L, n, "K")

  if (is.null(fold_id)) {
    fold_id <- sample(rep_len(seq_len(splits), n))
  } else if (!is.numeric(fold_id) || !is.null(dim(fold_id)) ||
             length(fold_id) != n || !all(fold_id %in% seq_len(splits)) ||
             !all(seq_len(splits) %in% fold_id)) {
    stop("fold_id must assign every row once and use every label from 1 to splits.")
  }
  
  sizes <- tabulate(as.integer(fold_id), nbins = splits)
  if (any(n - sizes < K)) stop("K exceeds the size of a regression sample.")

  # Validate packages, SL.lib and SL.env
  if (hal && !requireNamespace("hal9001", quietly = TRUE)) {
    stop("Install hal9001 to use hal = TRUE.")
  }
  if (!hal || !randomized) {
    if (!requireNamespace("SuperLearner", quietly = TRUE)) {
      stop("Install SuperLearner and the packages required by SL.lib.")
    }
    if (is.null(SL.env)) SL.env <- asNamespace("SuperLearner")
    if (!is.environment(SL.env)) stop("SL.env must be an environment.")
  }

  
  #-----------------------------------------------------------------------------
  # 2.0 Fit nuisance functions in each split
  #-----------------------------------------------------------------------------
  folds <- lapply(seq_len(splits), function(l) {
    r <- which(fold_id != l)
    m <- which(fold_id == l)
    if (length(unique(A[r])) != 2L) {
      stop(sprintf("Regression sample %d needs observations from both treatment arms.", l))
    }
    s <- .nf_fit_split(
      X_reg = X[r, , drop = FALSE], A_reg = A[r], Y_reg = Y[r],
      U = X[m, , drop = FALSE], kernel = kernel, family = family,
      hal = hal, randomized = randomized, pi_hat = pi_hat, K = K,
      SL.lib = SL.lib, SL.env = SL.env
    )
    s$regression_rows <- r
    s$marginal_rows <- m
    s
  })

  
  #-----------------------------------------------------------------------------
  # 3.0 Average over the marginal folds
  #-----------------------------------------------------------------------------
  weights <- sizes / n
  average <- function(f) {
    values <- lapply(folds, f)
    first <- values[[1L]]
    compatible <- vapply(values, function(v) {
      is.numeric(v) && length(v) == length(first) && identical(dim(v), dim(first))
    }, logical(1))
    if (!all(compatible)) stop("f must return numeric results of the same shape in every split.")
    Reduce(`+`, Map(function(w, v) w * v, weights, values))
  }
  list(folds = folds, fold_id = fold_id, x_names = x_names,
       weights = weights, average = average)
}





# Each call has a separate R environment. Only X is passed for the marginal rows.
#' Title .nf_fit_split
#'
#' @param X_reg covariates for the nuisance regressions, one row per observation
#' @param A_reg treatment for the nuisance regressions, one value per observation
#' @param Y_reg outcome for the nuisance regressions, one value per observation
#' @param U covariates for the marginal empirical measure, one row per observation
#' @param kernel specifies the kernel function k(U, x) for the embedding. 
#' It must return one numeric value per row of U when x is a one-row data frame.
#' @param family Family object for the outcome regression, either gaussian() or binomial().
#' @param hal Boolean indicating whether to use hal9001 for the outcome regression. If FALSE, SuperLearner is used.
#' @param randomized Boolean indicating whether the treatment is randomized. If TRUE, the known treatment probability pi_hat is used; if FALSE, a treatment model is fitted using SuperLearner.
#' @param pi_hat Known treatment probability P(A = 1) when randomized is TRUE. Must be strictly between 0 and 1.
#' @param K Number of internal cross-validation folds for SuperLearner and, if used, HAL. Must be at least 2 (or 3 if hal is TRUE).
#' @param SL.lib SuperLearner library for the outcome and the unknown propensity. It should be a character vector of SuperLearner wrapper names.
#' @param SL.env Optional environment for custom SuperLearner wrappers. Defaults to the SuperLearner namespace for its built-in wrappers.
#'
#' @returns A list containing the fitted nuisance functions and related quantities for the current split. The list includes:
#' - mu(a, x): Function to predict the outcome regression for treatment a and covariates x.
#' - pi(x): Function to predict the treatment probability for covariates x.
#' - cate(x): Function to compute the conditional average treatment effect for covariates x.
#' - centred_cate(x): Function to compute the centred conditional average treatment effect for covariates x.
#' - embedded_cate(x): Function to compute the embedded centred CATE for covariates x.
#' - residual(a, y, x): Function to compute the weighted residuals for treatment a, outcome y, and covariates x.
#' - P_X: Empirical measure function for the marginal covariate distribution.
#' - ate: Scalar representing the average treatment effect (mean of the fitted CATE).
#' - mean_embedding: Scalar representing the mean of the embedded centred CATE.
#' - embedding_norm_sq: Scalar representing the squared norm of the embedded centred CATE.
#' @export
#'
#' @examples
.nf_fit_split <- function(X_reg, A_reg, Y_reg, U, kernel, family, hal,
                          randomized, pi_hat, K, SL.lib, SL.env) {
  
  
  #-----------------------------------------------------------------------------
  # 1.0 Prepare data and helper functions
  #-----------------------------------------------------------------------------
  
  x_names <- names(X_reg)
  as_x <- function(x) .nf_covariates(x, x_names)
  
  # Specify the design matrix with treatment interactions for the outcome regression.
  design <- .nf_design(A_reg, X_reg)
  
  
  # Functions to fit and predict with SuperLearner
  fit_sl <- function(y, x, outcome_family) {
    SuperLearner::SuperLearner(
      Y = y, X = x, family = outcome_family, SL.library = SL.lib,
      cvControl = list(V = K), env = SL.env
    )
  }
  predict_sl <- function(fit, x, train_x, train_y) {
    stats::predict(fit, newdata = x, X = train_x, Y = train_y, onlySL = TRUE)$pred
  }

  
  
  #-----------------------------------------------------------------------------
  # 2.0 Fit the outcome regression and the treatment model using either HAL or SuperLearner.
  #-----------------------------------------------------------------------------
  
  if (hal) {
    
    # Fit the outcome regression with HAL, using the design matrix with treatment interactions.
    fit_mu <- hal9001::fit_hal(
      X = as.matrix(design), Y = Y_reg, family = family$family,
      yolo = FALSE, smoothness_orders = 1L, max_degree = 2L,
      num_knots = c(25L, 10L),
      fit_control = list(cv_select = TRUE, use_min = TRUE, nfolds = K,
                         lambda.min.ratio = 1e-4, prediction_bounds = "default")
      
    )
    
    # Define a prediction function for the outcome regression using the fitted HAL model.
    predict_mu <- function(d) {
      stats::predict(fit_mu, new_data = as.matrix(d), type = "response")
    }
    
  } else {
    
    # Fit the outcome regression with SuperLearner, using the design matrix with treatment interactions.
    fit_mu <- fit_sl(Y_reg, design, family)
    
    # Define a prediction function for the outcome regression using the fitted SuperLearner model.
    predict_mu <- function(d) predict_sl(fit_mu, d, design, Y_reg)
    
  }
  
  
  
  #-----------------------------------------------------------------------------
  # 3.0 Fit mu(a, x)
  #-----------------------------------------------------------------------------
  mu <- function(a, x) {
    d <- .nf_design(a, as_x(x))
    .nf_values(predict_mu(d), nrow(d), "Outcome predictions")
  }
  
  #-----------------------------------------------------------------------------
  # 3.1 Fit cate(x) = mu(1, x) - mu(0, x)
  #-----------------------------------------------------------------------------
  cate <- function(x) mu(1, x) - mu(0, x)
  
  #-----------------------------------------------------------------------------
  # 3.2 Fit pi(x)
  #-----------------------------------------------------------------------------
  # The treatment model uses raw covariates, with no treatment interactions.
  fit_pi <- NULL
  if (randomized) {
    # If the treatment is randomized, use the known constant treatment probability pi_hat.
    pi <- function(x) rep(pi_hat, nrow(as_x(x)))
  } else {
    # If the treatment is not randomized, fit a treatment model using SuperLearner with a binomial family.
    fit_pi <- fit_sl(A_reg, X_reg, stats::binomial())
    pi <- function(x) {
      x <- as_x(x)
      p <- .nf_values(predict_sl(fit_pi, x, X_reg, A_reg), nrow(x), "Propensity predictions")
      if (any(p < 0 | p > 1)) stop("A propensity prediction is outside [0, 1].")
      p
    }
  }
  
  #-----------------------------------------------------------------------------
  # 3.3 Fit the weighted residuals R(z) = (a / pi(x) - (1 - a) / (1 - pi(x))) * (y - mu(a, x))
  #-----------------------------------------------------------------------------
  residual <- function(a, y, x) {
    x <- as_x(x)
    a <- .nf_treatment(a, nrow(x))
    y <- .nf_values(y, nrow(x), "y")
    p <- pi(x)
    if (any(p <= 0 | p >= 1)) stop("The weighted residual requires 0 < pi(x) < 1.")
    ifelse(a == 1, 1 / p, -1 / (1 - p)) * (y - mu(a, x))
  }

  #-----------------------------------------------------------------------------
  # 3.4 Compute the empirical measure, plug-in ATE, and the centred CATE on the marginal sample U.
  #-----------------------------------------------------------------------------
  P_X <- .nf_empirical_measure(U)
  tau_U <- cate(U)
  ate <- mean(tau_U)
  c_U <- tau_U - ate
  centred_cate <- function(x) cate(x) - ate

  
  #-----------------------------------------------------------------------------
  # 3.5 Compute the embedded centred CATE on the marginal sample U.
  #-----------------------------------------------------------------------------
  # One kernel column at a time. The fitted values c_U are reused for every x.
  embedded_cate <- function(x) {
    x <- as_x(x)
    vapply(seq_len(nrow(x)), function(j) {
      column <- .nf_values(kernel(U, x[j, , drop = FALSE]), nrow(U), "Kernel column")
      mean(column * c_U)
    }, numeric(1))
  }
  G_U <- embedded_cate(U)
  
  
  list(mu = mu, 
       pi = pi, 
       cate = cate, 
       centred_cate = centred_cate,
       embedded_cate = embedded_cate, 
       residual = residual, 
       P_X = P_X, 
       ate = ate,
       mean_embedding = mean(G_U), 
       embedding_norm_sq = mean(c_U * G_U),
       fit_mu = fit_mu, 
       fit_pi = fit_pi)
}




# The same design is used for observed and counterfactual treatments.
.nf_design <- function(a, x) {
  a <- .nf_treatment(a, nrow(x))
  interactions <- as.matrix(x) * a
  d <- data.frame(A = a, x, interactions, check.names = FALSE)
  names(d) <- make.names(c("A", names(x), paste0("A_", names(x))), unique = TRUE)
  d
}

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

.nf_count <- function(x, lower, upper, label) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x != floor(x) || x < lower || x > upper) {
    stop(sprintf("%s must be a whole number between %d and %d.", label, lower, upper))
  }
}

.nf_treatment <- function(a, n) {
  if (!(is.numeric(a) || is.logical(a)) || !(length(a) %in% c(1L, n)) ||
      !all(a %in% 0:1)) {
    stop("Treatment must be 0 or 1, with length one or one value per observation.")
  }
  rep_len(as.numeric(a), n)
}

.nf_values <- function(v, n, label) {
  if (!(is.numeric(v) || is.logical(v)) || length(v) != n || any(!is.finite(v))) {
    stop(paste(label, "must supply one finite number per observation."))
  }
  as.numeric(v)
}

# Private copies of the earlier helpers make this script self-contained.
# A vector is one observation. Matrices and data frames have one observation
# per row, with features in the same order in both inputs.
.nf_gaussian_kernel <- function(x, y, sigma = 1) {
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
