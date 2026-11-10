############################################################
#                                                          #
#   The following script includes functions to estimate    #
#   the nuisance parameters in an outer fold of the        #
#   cross-fitting procedure. It includes functions to      #
#   to estimate                                            #
#     - mu(a, x): Outcome regression for treatment a       #
#                 and covariates x.                        #
#     - pi(x): Treatment probability for covariates x.     #
#     - cate(x): Conditional average treatment effect      #
#                for covariates x.                         #
#     - centred_cate(x): Centred conditional average       #
#                        treatment effect for covariates x.#
#     - embedded_cate(x): Embedded centred CATE for        # 
#                         covariates x.                    #
#     - residual(a, y, x): The weighted residuals R(z)     #
#                          for treatment a, outcome y,     #
#                          and covariates x.               #
#     - P_X: The marginal covariate distribution.          #
#     - ate: The average treatment effect (mean of         #
#            the fitted CATE).                             #
#     - mean_embedding: The mean of the embedded           #
#                       centred CATE.                      #
#     - embedding_norm_sq: The squared norm of the         #
#                          embedded centred CATE           #
#                          corresponding to the embedded   #
#                          parameter.                      #
#                                                          #
#   The nuisance functions are estimated using either      #
#   SuperLearner or HAL9001, depending on the user's       #
#   choice.                                                #
#                                                          #
############################################################


############################################################
#                         Packages                         #
############################################################

# For estimation
library(glmnet)
library(SuperLearner)
library(grf)
library(hal9001)

# For data wrangling
library(tidyverse)


############################################################
#                         Dependencies                     #
############################################################

# Load helper functions
source("nf_helpers.R")

# Load kernel functions
source("kernels.R")

############################################################
#       Function for estimating nuisance parameters        #
############################################################






############################################################
#                                                          #
#   Function to estimate nuisance functions for the        #
#   embedded centred CATE with an inner cross-fitting      #
#   procedure. In particular, it splits the training data  #
#   into two folds. On one fold it fits the nuisance       #
#   functions and evaluates them on the other fold which   #
#   corresponds to training the empirical covariate        #
#   marginal on a seperate fold.                           #
#                                                          #
############################################################

#' Cross-fitted nuisance estimates for the embedded centred CATE
#'
#' Supply only the training data of the current outer split. In inner split l,
#' regressions use rows with fold_id != l and the empirical covariate marginal
#' uses rows with fold_id == l. Each element of folds has its own fitted models.
#'
#' Title nuisance_func
#'
#' @param train_data Data frame with numeric Y, binary A (0/1), and finite
#'   numeric covariates. All columns except Y and A are treated as covariates.
#' @param inner_splits Number of marginal folds, at least two.
#' @param SL.lib SuperLearner library for the outcome and the unknown propensity.
#' @param family gaussian() for a continuous Y or binomial() for binary Y.
#' @param hal Use hal9001 for the outcome regression. An unknown propensity
#'   is fitted with SuperLearner in either case.
#' @param randomized Use the known constant treatment probability pi_hat.
#' @param pi_hat Known P(A = 1) when randomized is TRUE.
#' @param K Number of internal CV folds for SuperLearner and, if used, HAL.
#' @param kernel Function k(U, x) returning one numeric value per row of U
#'   when x is a one-row data frame. The default is Gaussian with median heuristic bandwidth.
#' @param bandwidth_pairs Maximum number of pairs to use for estimating the Gaussian bandwidth when kernel = "median". 
#'  If Inf, all pairs are used. If fewer than choose(n, 2) pairs are requested, they are sampled uniformly with replacement. 
#'  Use set.seed() before this call for a reproducible approximation.
#' @param fold_id Optional integer labels from 1 to splits, one per training
#'   row. The default is a balanced random partition independent of the data.
#' @param SL.env Optional environment for custom SuperLearner wrappers.
#'   Defaults to the SuperLearner namespace for its built-in wrappers.
#' @return A list with folds, fold_id, x_names, weights, and average(f).
#'   Each fold contains its row indices, P_X, fitted models, functions mu(a,x),
#'   pi(x), cate(x), centred_cate(x), embedded_cate(x), residual(z),
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
#' }
nuisance_func <- function(train_data, 
                          inner_splits = 2L,
                          SL.lib = c("SL.glmnet", "SL.randomForest", "SL.xgboost"),
                          family = stats::gaussian(), 
                          hal = FALSE,
                          randomized = TRUE, 
                          pi_hat = 0.5, 
                          K = 5L,
                          kernel = .nf_gaussian_kernel, 
                          fold_id = NULL,
                          bandwidth_pairs = Inf,
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
  
  # Validate kernel
  if (!is.function(kernel) && !identical(kernel, "median")) {
    stop('Supply a kernel function or kernel = "median".')
  }
  
  # Validate splits and K
  .nf_count(inner_splits, 2L, n, "inner_splits")
  .nf_count(K, if (hal) 3L else 2L, n, "K")

  # Validate fold_id
  if (is.null(fold_id)) {
    fold_id <- sample(rep_len(seq_len(inner_splits), n))
  } else if (!is.numeric(fold_id) || !is.null(dim(fold_id)) ||
             length(fold_id) != n || !all(fold_id %in% seq_len(inner_splits)) ||
             !all(seq_len(inner_splits) %in% fold_id)) {
    stop("fold_id must assign every row once and use every label from 1 to inner_splits.")
  }
  
  # Validate that each regression sample has at least K observations in each treatment arm.
  sizes <- tabulate(as.integer(fold_id), nbins = inner_splits)
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
  folds <- lapply(seq_len(inner_splits), function(l) {
    
    r <- which(fold_id != l)      # Regression sample rows
    m <- which(fold_id == l)      # Marginal sample rows
    
    if (length(unique(A[r])) != 2L) {
      stop(sprintf("Regression sample %d needs observations from both treatment arms.", l))
    }
    
    # Fit the nuisance functions for the current split using the regression and marginal samples.
    s <- .nf_fit_split(
      X_reg = X[r, , drop = FALSE], 
      A_reg = A[r], 
      Y_reg = Y[r],
      U = X[m, , drop = FALSE], 
      kernel = kernel,
      bandwidth_pairs = bandwidth_pairs,
      family = family,
      hal = hal, 
      randomized = randomized, 
      pi_hat = pi_hat, 
      K = K,
      SL.lib = SL.lib, 
      SL.env = SL.env
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
  
  
  return(list(folds = folds, fold_id = fold_id, x_names = x_names,
              weights = weights, average = average))
}






############################################################
#                                                          #
#   Function to estimate nuisance functions for the        #
#   embedded centred CATE in a single inner fold. In       #
#   particular, it splits the training data into two       #
#   folds. On one fold it fits the nuisance functions and  #
#   evaluates them on the other fold which corresponds     #
#   to training the empirical covariate marginal on a      #
#   seperate fold.                                         #
#                                                          #
############################################################


#' Title .nf_fit_split
#'
#' @param X_reg covariates for the nuisance regressions, one row per observation
#' @param A_reg treatment for the nuisance regressions, one value per observation
#' @param Y_reg outcome for the nuisance regressions, one value per observation
#' @param U covariates for the marginal empirical measure, one row per observation
#' @param kernel specifies the kernel function k(U, x) for the embedding. 
#' It must return one numeric value per row of U when x is a one-row data frame.
#' @param bandwidth_pairs Maximum number of pairs to use for estimating the Gaussian bandwidth when kernel = "median". If Inf, all pairs are used. If fewer than choose(n, 2) pairs are requested, they are sampled uniformly with replacement. Use set.seed() before this call for a reproducible approximation.
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
.nf_fit_split <- function(X_reg, A_reg, Y_reg, U, 
                          kernel, family, bandwidth_pairs = Inf, 
                          randomized, pi_hat, 
                          hal,
                          K, SL.lib, SL.env) {
  
  
  #-----------------------------------------------------------------------------
  # 0.0 Prepare kernel and bandwidth
  #-----------------------------------------------------------------------------
  if (identical(kernel, "median")) {
    kernel <- make_median_gaussian_kernel(X_reg, max_pairs = bandwidth_pairs)
  }
  sigma <- if (identical(kernel, .nf_gaussian_kernel)) 1 else attr(kernel, "sigma", exact = TRUE)
  
  
  #-----------------------------------------------------------------------------
  # 1.0 Prepare data
  #-----------------------------------------------------------------------------
  
  # Prepare and validate covariates
  x_names <- names(X_reg)
  as_x <- function(x) .nf_covariates(x, x_names)
  
  # Specify the design matrix with treatment interactions for the outcome regression.
  design <- .nf_design(A_reg, X_reg)
  
  
  #-----------------------------------------------------------------------------
  # 2.0 Fit the outcome regression and the treatment model using either HAL or SuperLearner.
  #-----------------------------------------------------------------------------
  
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
  # 3.5 Compute embedded parameters on the marginal sample U.
  #-----------------------------------------------------------------------------
  # Reduce columns k(U, point) without storing an m-by-m Gram matrix.
  kernel_column <- function(rows, point) {
    .nf_values(kernel(rows, point), nrow(rows), "Kernel column")
  }
  
  marginal_columns <- function(points, f, prototype = numeric(1)) {
    vapply(seq_len(nrow(points)), function(j) {
      f(kernel_column(U, points[j, , drop = FALSE]))
    }, prototype)
  }
  
  #-----------------------------------------------------------------------------
  # 3.5.1 Compute embedded cate on the marginal sample U.
  #-----------------------------------------------------------------------------
  # Define the embedded centred CATE as the mean of the kernel-weighted centred CATE over the marginal sample U.
  cate_average <- function(column) mean(column * c_U)
  embedded_cate <- function(x) marginal_columns(as_x(x), cate_average)
  
  # Cache both marginal sums in the same pass over all ordered pairs in U.
  means_on_U <- marginal_columns(U, function(column) {
    c(kernel = mean(column), cate = cate_average(column))
  }, c(kernel = 0, cate = 0))
  
  G_U <- means_on_U["cate", ]
  
  
  #-----------------------------------------------------------------------------
  # 3.5.2 Compute mean of the kernel on the marginal sample U.
  #-----------------------------------------------------------------------------
  mean_UU <- mean(means_on_U["kernel", ])
  
  
  #-----------------------------------------------------------------------------
  # 3.6 Compute centred kernel on the marginal sample U.
  #-----------------------------------------------------------------------------
  centred_kernel <- function(x, y) {
    x <- as_x(x)
    y <- as_x(y)
    mean_Ux <- marginal_columns(x, mean)
    mean_Uy <- marginal_columns(y, mean)
    drop(vapply(seq_len(nrow(y)), function(j) {
      kxy_j <- kernel_column(x, y[j, , drop = FALSE])
      kxy_j - mean_Ux - mean_Uy[j] + mean_UU
    }, numeric(nrow(x))))
  }
  
  
  return(list(mu = mu, 
       pi = pi, 
       cate = cate, 
       centred_cate = centred_cate,
       embedded_cate = embedded_cate, 
       centred_kernel = centred_kernel,
       residual = residual, 
       P_X = P_X, 
       ate = ate,
       mean_embedding = mean(G_U), 
       embedding_norm_sq = mean(c_U * G_U),
       fit_mu = fit_mu, 
       fit_pi = fit_pi))
}




