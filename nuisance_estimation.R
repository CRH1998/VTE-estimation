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

# For simulation
library(MASS)

# For kernels
library(kerntools)

############################################################


############################################################
#       Function for estimating nuisance parameters        #
############################################################




#' Title empirical_measure
#'
#' @param X A vector, list, matrix, or data frame of observations. Each observation can be a numeric vector, a list, or a data frame row.
#'
#' @returns A function that takes a function f and an optional boolean vectorized argument. The returned function computes the empirical measure of f over the observations in X, either by applying f to each observation individually (if vectorized is FALSE) or by applying f to the entire dataset at once (if vectorized is TRUE).
#' @export
#'
#' @examples
empirical_measure <- function(X) {
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



#' Title Nuisance parameter estimation for embedded parameter estimation
#'
#' @param train_data A data frame containing the training data, which should include the outcome variable Y, the treatment variable A, and covariates X. The function will use this data to estimate the nuisance functions needed for the pseudo-VCATE estimation.
#' @param eval_data A data frame containing the evaluation data, which should have the same structure as train_data. The function will use this data to predict the nuisance functions (mu0_hat, mu1_hat, and pi_hat) for the evaluation set, which are then used in the pseudo-VCATE estimation.
#' @param SL.lib A character vector of SuperLearner algorithms to use for estimating the nuisance functions. Default includes glmnet, random forest, and xgboost.
#' @param family The family argument for SuperLearner, which specifies the type of outcome. Default is gaussian() for continuous outcomes. For binary outcomes, use binomial().
#' @param randomized A boolean indicating whether the propensity score is known (e.g., from a randomized experiment). If TRUE, the function will use the provided pi_hat as the propensity score. If FALSE, it will estimate the propensity score using SuperLearner. Default is TRUE.
#' @param pi_hat If randomized is TRUE, this should be a numeric value representing the known propensity score (e.g., 0.5 for a balanced randomized experiment). If randomized is FALSE, this argument is ignored. Default is 0.5.
#' @param K Number of folds for cross-validation in SuperLearner. Default is 5.
#'
#' @returns A list containing: S_nk, which is the machine learning proxy for average treatment effect estimates (tau_hat(X) - ate_hat); M_nk, which is the machine learning proxy for outcome under treatment 0 (mu0(X)); and pi_hat, which is either the known propensity score or the estimated propensity score from SuperLearner.
#' @export
#'
#' @examples
nuisance_param_func <- function(train_data, eval_data,
                                kernel = gaussian_kernel,
                                SL.lib = c("SL.glmnet", "SL.randomForest", "SL.xgboost"),
                                family = gaussian(),
                                hal = FALSE, 
                                K = 5) {
  #browser()
  # Training data
  Y_tr <- train_data$Y
  A_tr <- train_data$A
  X_tr <- train_data %>% dplyr::select(-Y, -A)
  
  # Creating the design matrix for SuperLearner, including interactions between treatment and covariates.
  X_design <- cbind(
    A_tr,
    X_tr,
    as.data.frame(A_tr * as.matrix(X_tr))
  )
  names(X_design) <- c("A", names(X_tr), paste0("A_", names(X_tr)))
  
  # Evaluation data
  Y_ev <- eval_data$Y
  A_ev <- eval_data$A
  X_ev <- eval_data %>% dplyr::select(-Y, -A)
  
  newdata_a <- cbind(
    A = A_ev,
    X_ev,
    as.data.frame(A_ev * as.matrix(X_ev))
  )
  
  
  newdata_0 <- cbind(
    A = 0,
    X_ev,
    as.data.frame(0 * as.matrix(X_ev))
  )
  names(newdata_0) <- names(X_design)
  
  newdata_1 <- cbind(
    A = 1,
    X_ev,
    as.data.frame(1 * as.matrix(X_ev))
  )
  names(newdata_1) <- names(X_design)
  
  if (hal) {
    
    fit_Q <- fit_hal(X = X_design, Y = Y_tr, yolo = FALSE, smoothness_orders = 1,
                     max_degree = 2, num_knots = c(25, 10, 5))
    
    mu0_hat <- predict(fit_Q, new_data = newdata_0)
    mu1_hat <- predict(fit_Q, new_data = newdata_1)
    tau_hat <- mu1_hat - mu0_hat
    bar_tau_hat <- mean(tau_hat)
    
  } else {
    
    
    
    # Fitting SuperLearner
    fit_Q <- SuperLearner::SuperLearner(
      Y = Y_tr,
      X = X_design,
      family = family,
      SL.library = SL.lib,
      cvControl = list(V = K)
    )
    
    mua_hat <- predict(fit_Q, newdata = X_design)$pred
    mu0_hat <- predict(fit_Q, newdata = newdata_0)$pred
    mu1_hat <- predict(fit_Q, newdata = newdata_1)$pred
    tau_hat <- mu1_hat - mu0_hat
    bar_tau_hat <- mean(tau_hat)
    
  }
  
  list(
    mua_hat = mua_hat,
    mu0_hat = mu0_hat,
    mu1_hat = mu1_hat,
    tau_hat = tau_hat,
    bar_tau_hat = bar_tau_hat
  )
}






#' Title nuisance parameter estimation for embedded parameter estimation
#'
#' @param train_data A data frame containing the training data, which should include the outcome variable Y, the treatment variable A, and covariates X. The function will use this data to estimate the nuisance functions needed for the pseudo-VCATE estimation.
#' @param SL.lib A character vector of SuperLearner algorithms to use for estimating the nuisance functions. Default includes glmnet, random forest, and xgboost.
#' @param family The family argument for SuperLearner, which specifies the type of outcome. Default is gaussian() for continuous outcomes. For binary outcomes, use binomial().
#' @param randomized A boolean indicating whether the propensity score is known (e.g., from a randomized experiment). If TRUE, the function will use the provided pi_hat as the propensity score. If FALSE, it will estimate the propensity score using SuperLearner. Default is TRUE.
#' @param pi_hat If randomized is TRUE, this should be a numeric value representing the known propensity score (e.g., 0.5 for a balanced randomized experiment). If randomized is FALSE, this argument is ignored. Default is 0.5.
#' @param K Number of folds for cross-validation in SuperLearner. Default is 5.
#' @param splits Number of splits for cross-fitting. The data will be randomly split into this many folds, and the nuisance parameters will be estimated on one fold while the marginal P_X is estimated on the remaining data. Default is 2.
#' @param hal 
#'
#' @returns A list containing:
#' @export
#'
#' @examples
nuisance_func <- function(train_data,
                          splits = 2,
                          SL.lib = c("SL.glmnet", "SL.randomForest", "SL.xgboost"),
                          family = gaussian(),
                          hal = FALSE,
                          randomized = TRUE, pi_hat = 0.5, K = 5) {
  #browser()
  # Training data
  Y_tr <- train_data$Y
  A_tr <- train_data$A
  X_tr <- train_data %>% dplyr::select(-Y, -A)
  
  # Creating the design matrix for SuperLearner, including interactions between treatment and covariates.
  X_design <- cbind(
    A_tr,
    X_tr,
    as.data.frame(A_tr * as.matrix(X_tr))
  )
  names(X_design) <- c("A", names(X_tr), paste0("A_", names(X_tr)))
  
  
  
  
  split_indices <- sample(1:splits, nrow(train_data), replace = TRUE)
  folds <- vector("list", splits)
  
  for (i in 1:splits) {
    
    
    #---------------------------------------------------------------------------
    # 0.0 Split the data into training and evaluation sets for cross-fitting
    #---------------------------------------------------------------------------
    
    marginal_data_split <- X_design[split_indices == i, ]
    nuis_data_split <- X_design[split_indices != i, ]
    
    
    #---------------------------------------------------------------------------
    # 1.0 Estimate marginal distribution P_X using empirical measure on marginal_data_split
    #---------------------------------------------------------------------------
    
    P_X_hat <- empirical_measure(marginal_data_split %>% dplyr::select(-Y, -A))
    
    
    
    
    #-----------------------------------------------------------------------------
    # 2.1 Estimate CATE using HAL or SuperLearner
    #-----------------------------------------------------------------------------
    if (hal) {
      
      fit_CATE <- fit_hal(X = nuis_data_split, Y = Y_tr, yolo = FALSE, smoothness_orders = 1,
                          max_degree = 2, num_knots = c(25, 10, 5))
      
    } else {
      
      # Fitting SuperLearner
      fit_CATE <- SuperLearner::SuperLearner(
        Y = Y_tr,
        X = nuis_data_split,
        family = family,
        SL.library = SL.lib,
        cvControl = list(V = K)
      )
      
    }
    
    
    #-----------------------------------------------------------------------------
    # 2.2 Estimate propensity scores
    #-----------------------------------------------------------------------------
    if (randomized) {
      pi_hat_eval <- rep(pi_hat, nrow(eval_data))
    } else {
      fit_pi <- SuperLearner(
        Y = A_tr,
        X = X_tr,
        family = binomial(),
        SL.library = SL.lib
      )
      pi_hat_eval <- predict(fit_g, newdata = X_ev)$pred
    }
    
    
    #---------------------------------------------------------------------------
    # 2.3 Estimate ATE using nuisance and marginal estimates from the training data
    #---------------------------------------------------------------------------
    fit_ATE <- P_X_hat(function(x) {
      newdata_0 <- cbind(
        A = 0,
        x,
        as.data.frame(0 * as.matrix(x))
      )
      names(newdata_0) <- names(X_design)
      
      newdata_1 <- cbind(
        A = 1,
        x,
        as.data.frame(1 * as.matrix(x))
      )
      names(newdata_1) <- names(X_design)
      
      mu0_hat <- predict(fit_CATE, newdata = newdata_0)$pred
      mu1_hat <- predict(fit_CATE, newdata = newdata_1)$pred
      tau_hat <- mu1_hat - mu0_hat
      
      return(tau_hat)
    }, vectorized = TRUE)
    
    
    
    #-----------------------------------------------------------------------------
    # 2.4 Estimate kernel embedding of the centred CATE function
    #-----------------------------------------------------------------------------
    P_X_hat(function(U) {
      newdata_0 <- cbind(
        A = 0,
        U,
        as.data.frame(0 * as.matrix(U))
      )
      names(newdata_0) <- names(X_design)
      
      newdata_1 <- cbind(
        A = 1,
        U,
        as.data.frame(1 * as.matrix(U))
      )
      names(newdata_1) <- names(X_design)
      
      
      
      gaussian_kernel(x, U) * f_hat(U)
    }, 
    vectorized = TRUE)
    
    
    
    
    centred_cate <- function(x) {
      newdata_0 <- cbind(
        A = 0,
        x,
        as.data.frame(0 * as.matrix(x))
      )
      names(newdata_0) <- names(X_design)
      
      newdata_1 <- cbind(
        A = 1,
        x,
        as.data.frame(1 * as.matrix(x))
      )
      names(newdata_1) <- names(X_design)
      
      mu0_hat <- predict(fit_CATE, newdata = newdata_0)$pred
      mu1_hat <- predict(fit_CATE, newdata = newdata_1)$pred
      tau_hat <- mu1_hat - mu0_hat
      
      return(tau_hat - fit_ATE)
    }
    

  }
  
  list(
    fit_ATE = fit_ATE,
    fit_CATE = fit_CATE,
    fit_pi = fit_pi
  )
}




#' Title Average Treatment Effect (ATE) estimation with Cross-Fitting
#'
#' @param train_data Training data for estimating the nuisance parameters and the average treatment effect (ATE). This should be a data frame containing the outcome variable Y, the treatment variable A, and covariates X.
#' @param splits Number of splits for cross-fitting. The data will be randomly split into this many folds, and the nuisance parameters will be estimated on each fold while evaluating on the remaining data. Default is 2.
#' @param SL.lib A character vector of SuperLearner algorithms to use for estimating the nuisance functions. Default includes random forest.
#' @param family The family argument for SuperLearner, which specifies the type of outcome. Default is gaussian() for continuous outcomes. For binary outcomes, use binomial().
#' @param hal A boolean indicating whether to use the Highly Adaptive Lasso (HAL) for estimating the nuisance functions. If TRUE, HAL will be used; if FALSE, SuperLearner will be used. Default is FALSE.
#' @param K Number of folds for cross-validation in SuperLearner. Default is 5.
#'
#' @returns A cross-fitted estimate of  hat P_X tau_{hat P}(X)
#' @export
#'
#' @examples
ate_est <- function(train_data,
                    splits = 2,
                    SL.lib = c("SL.randomForest"),
                    family = gaussian(),
                    hal = FALSE,
                    K = 5) {
  
  #browser()
  # Split the data into n_splits splits for cross-fitting
  split_indices <- sample(1:splits, nrow(train_data), replace = TRUE)
  
  folds <- vector("list", splits)
  
  for (i in 1:splits) {
    folds[i] <- nuisance_param_func(train_data = train_data[split_indices != i, ], 
                                    eval_data = train_data[split_indices == i, ], 
                                    SL.lib = SL.lib, 
                                    family = family, 
                                    hal = hal,
                                    K = K)$bar_tau_hat
  }
  
  return(mean(unlist(folds)))
}



