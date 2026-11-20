############################################################
#                                                          #
#   The following script includes estimation functions     #
#   for the test of heterogenous effects simulation study  #
#                                                          #
############################################################


############################################################
library(glmnet)
library(SuperLearner)
library(tidyverse)
library(grf)
library(MASS)
library(hal9001)

############################################################


############################################################
#                 Sánchez-Becerra functions                #
############################################################


############################################################
#       Function for estimating nuisance parameters        #
############################################################


#' Title Nuissance parameter estimation for pseudo-VCATE estimation
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
  
  # Testing data
  Y_ev <- eval_data$Y
  A_ev <- eval_data$A
  X_ev <- eval_data %>% dplyr::select(-Y, -A)
  
  
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
    # S_nk <- predict(fit_Q, new_data = cbind(A_ev, X_ev))
    # M_nk <- predict(fit_Q, new_data = cbind(0, X_ev))
    
    M_nk <- predict(fit_Q, new_data = newdata_0)
    mu1_hat <- predict(fit_Q, new_data = newdata_1)
    S_nk <- mu1_hat - M_nk
    
  } else {
    

    
    # Fitting SuperLearner
    fit_Q <- SuperLearner::SuperLearner(
      Y = Y_tr,
      X = X_design,
      family = family,
      SL.library = SL.lib,
      cvControl = list(V = K)
    )
    
    M_nk <- predict(fit_Q, newdata = newdata_0)$pred
    mu1_hat <- predict(fit_Q, newdata = newdata_1)$pred
    S_nk <- mu1_hat - M_nk
  }
  
  # Propensity scores
  if (randomized) {
    pi_hat_eval <- rep(pi_hat, nrow(eval_data))
  } else {
    fit_g <- SuperLearner(
      Y = A_tr,
      X = X_tr,
      family = binomial(),
      SL.library = SL.lib
    )
    pi_hat_eval <- predict(fit_g, newdata = X_ev)$pred
  }
  
  list(
    S_nk = S_nk,
    M_nk = M_nk,
    pi_hat = pi_hat_eval
  )
}



############################################################
#   Functions for estimating pseudo-VCATE and Omega        #
############################################################

#-------------------------------------------------------------------------------

#' Title Sanchéz-Becerra pseudo-VCATE estimate
#'
#' @param S_nk Machine learning proxy for average treatment effect estimates
#' @param M_nk Machine learning proxy for outcome under treatment 0, mu0(X)
#' @param data The entire data frame
#'
#' @returns Estimated pseudo-VCATE as in Sánchez-Becerra (2023)
#' @export
#'
#' @examples
pseudoVCATE_SB <- function(S_nk, M_nk, pi_hat, data){
  
  # Estimating S_hat
  avg_S_nk <- mean(S_nk)
  Snk_hat <- S_nk - avg_S_nk
  
  # Estimating M_hat
  Mnk_hat <- M_nk + pi_hat * S_nk
  
  # Estimating W_hat
  Wnk_hat <- cbind(1, Mnk_hat, (data$A - pi_hat), (data$A - pi_hat) * Snk_hat)
  
  # Estimating theta_hat using weighted least squares regression of Y on W with weights lambda = 1 / (pi_hat * (1 - pi_hat))
  WLS_thetank_hat <- lm.wfit(
    x = Wnk_hat,
    y = data$Y,
    w = as.vector(1 / (pi_hat * (1 - pi_hat)))
  )
  thetank_hat <- unname(coef(WLS_thetank_hat))
  thetank_hat[is.na(thetank_hat)] <- 0
  
  # Estimating V_Snk_hat
  VSnk_hat <- mean(Snk_hat^2)
  
  # For numerical stability, ensure VSnk_hat is not too close to zero
  if (VSnk_hat < 1e-12) {
    return(list(
      V_taunk_star_hat = 0,
      Omegank_hat = matrix(NA_real_, nrow = 2, ncol = 2),
      thetank_hat = thetank_hat
    ))
  }
  
  # Estimating the pseudo-VCATE estimate, V_taunk^*_hat
  V_taunk_star_hat <- VSnk_hat * thetank_hat[4]^2
  
  # Estimating robust sandwich estimator of the asymptotic variance of the pseudo-VCATE estimate, Omega_nk
  Omegank_hat <- RSE_SB(data = data, VSnk_hat = VSnk_hat, Snk_hat = Snk_hat, Wnk_hat = Wnk_hat, thetank_hat = thetank_hat, pi_hat = pi_hat)
  
  return(list(V_taunk_star_hat = V_taunk_star_hat, Omegank_hat = Omegank_hat, thetank_hat = thetank_hat))
}
#-------------------------------------------------------------------------------




#-------------------------------------------------------------------------------
#' Title
#'
#' @param data The entire data frame
#' @param VSnk_hat The estimated VSnk_hat, which is the variance of the S_i = tau_hat(X_i) - ate_hat, as in the SB paper.
#' @param S The vector of estimated S_i = tau_hat(X_i) - ate_hat, as in the SB paper.
#' @param W The matrix of covariates used in the regression of Y on W with weights lambda, as in the SB paper. Note that W includes an intercept, M_hat, (D - pi), and (D - pi) * S.
#' @param thetank_hat The estimated theta vector from the regression of Y on W with weights lambda, as in the SB paper
#' @param pi The propensity score estimates, which are constant in the SB paper but could be estimated in general.
#'
#' @returns An robust estimator of the asymptotic variance of the VCATE estimator as in Sánchez-Becerra (2023)
#' @export
#'
#' @examples
RSE_SB <- function(data, VSnk_hat, Snk_hat, Wnk_hat, thetank_hat, pi_hat){
  # Sample size in fold k
  n <- nrow(data)
  
  # Selection matrix
  SLM <- rbind(c(0,0,0,1,0),
               c(0,0,0,0,1))
  
  # Estimating two auxillary residuals, T_i_hat and U_i_hat
  RT <- 1/VSnk_hat * Snk_hat^2 - 1
  RU <- data$Y - Wnk_hat %*% thetank_hat
  
  
  # Constructing diagonal matrix Pi_nk_hat
  Pi <- diag(c(1,1,1, sqrt(1 / VSnk_hat)))
  
  # Calculating lambda_Xi
  lambda <- c(1/(pi_hat * (1 - pi_hat)))
  
  
  # Constructing matrix J_hat_nk-------------------------------------------#
  J_11 <- (Pi %*% crossprod(Wnk_hat, Wnk_hat * lambda) %*% t(Pi)) / n
  
  J <- matrix(0, nrow = 5, ncol = 5)
  J[1:4, 1:4] <- J_11
  J[5, 5] <- 1
  #------------------------------------------------------------------------#
  
  # Constructing matrix H_hat_nk-------------------------------------------#
  # Top-left block
  wA <- c((lambda * RU)^2)                        # length nk, weights = lambda^2 * U^2
  H11_temp  <- crossprod(Wnk_hat, Wnk_hat * wA)   # 4x4, equals t(W) %*% diag(wA) %*% W
  H11 <- (Pi %*% H11_temp %*% t(Pi)) / n          # 4x4
  
  # Top-right block
  wB  <- c(lambda * RU * RT)                  # length nk
  H12_temp   <- crossprod(Wnk_hat, wB)        # 4x1
  H12 <- (Pi %*% H12_temp) / n                # 4x1
  
  # Bottom-right scalar
  H22 <- mean(RT^2)
  
  # Assemble 5x5 H_nk
  H <- rbind(
    cbind(H11, H12),
    cbind(t(H12), H22)
  )
  #------------------------------------------------------------------------#
  
  
  # Robust sandwich estimator
  # Try solve(J) first, if it fails, try with the generalized inverse
  Jinv <- tryCatch(
    solve(J),
    error = function(e) MASS::ginv(J)
  )
  
  Omega <- SLM %*% Jinv %*% H %*% t(Jinv) %*% t(SLM)
  
  
  # For numerical stability, ensure Omega is symmetric and positive definite
  Omega <- (Omega + t(Omega)) / 2
  eig <- eigen(Omega, symmetric = TRUE, only.values = TRUE)$values
  if (min(eig) < 1e-8) {
    Omega <- as.matrix(Matrix::nearPD(Omega, corr = FALSE)$mat)
  }
  
  
  return(Omega)
}
#-------------------------------------------------------------------------------





############################################################
#   Function constructing adaptive confidence intervals    #
############################################################


# Empirical CDF of G at value 'v', for given value of pseudo-VCATE and zeta
#' Title Empirical CDF of the distribution of the pseudo-VCATE estimator for given parameters
#'
#' @param Omegank_hat The estimated asymptotic variance matrix of the pseudo-VCATE estimator, as returned by the pseudoVCATE_SB function. This should be a 2x2 matrix corresponding to the variance of the leading normal term and the variance of the quadratic term in the asymptotic distribution of the pseudo-VCATE estimator.
#' @param nk The sample size in the fold for which we want to construct the confidence interval. This is used to scale the terms in the distribution of the pseudo-VCATE estimator.
#' @param M The number of Monte Carlo samples to draw for approximating the distribution of the pseudo-VCATE estimator. A larger M will give a more accurate approximation but will take more computational time. Default is 500,000.
#' @param seed Random seed for reproducibility of the Monte Carlo draws. Default is 1.
#'
#' @returns A function that takes as input a value 'v', a pseudo-VCATE estimate 'pVCATE', and a sign 'zeta' (which can be -1 or 1), and returns the empirical CDF of the random variable G at 'v', where G is defined as in the adaptive confidence interval construction in Sánchez-Becerra (2023). This function can be used to compute the critical values for the adaptive confidence interval by evaluating it at different values of 'v' and 'pVCATE'.
#' @export
#'
#' @examples
Fn_hat <- function(Omegank_hat, nk, M = 4000, seed = NA) {
  #browser()
  
  if (!is.na(seed)) set.seed(seed)
  cholO <- chol(Omegank_hat)
  gamma <- matrix(rnorm(M * 2), ncol = 2)
  e1 <- c(1,0)
  e2 <- c(0,1)
  
  e1chol0 <- e1 %*% cholO
  e2chol0 <- e2 %*% cholO
  
  e1chol0gamma <- gamma %*% t(e1chol0)
  e2chol0gamma <- gamma %*% t(e2chol0)
  
  return(function(v, pVCATE, zeta) {
    sqrtpVCATE <- sqrt(pVCATE) 
    G <- (e1chol0gamma^2) / nk + 
      (2 * zeta * sqrtpVCATE * e1chol0gamma) / sqrt(nk) + 
      (pVCATE * e2chol0gamma) / sqrt(nk)
    mean(G <= v)
  }
  )
}


#' Title
#'
#' @param pVCATE The candidate pseudo-VCATE value for which we want to compute the critical values of the adaptive confidence interval. This should be a non-negative numeric value, as the pseudo-VCATE represents a variance-like quantity. The critical values will depend on this candidate value, as well as on the sign 'zeta' which indicates whether we are considering the upper or lower tail of the distribution of the pseudo-VCATE estimator.
#' @param zeta The sign indicating whether we are considering the upper tail (zeta = 1) or the lower tail (zeta = -1) of the distribution of the pseudo-VCATE estimator. This is used in the definition of the random variable G whose distribution we are approximating with the empirical CDF Fn_hat. The critical values for the adaptive confidence interval will differ depending on whether we are looking at the upper or lower tail, which is why zeta is an input to this function.
#' @param alpha The significance level for the confidence interval. This should be a numeric value between 0 and 1, typically something like 0.05 for a 95% confidence interval. The critical values computed by this function will be based on this alpha level, as they will determine the quantiles of the distribution of G that correspond to the desired coverage probability of the confidence interval.
#'
#' @returns A named numeric vector containing the lower and upper critical values (qL and qU) for the adaptive confidence interval corresponding to the given candidate pseudo-VCATE value, sign zeta, and significance level alpha. These critical values are computed based on the empirical CDF of the random variable G, which is defined in terms of the parameters of the pseudo-VCATE estimator's distribution. The lower critical value qL corresponds to the quantile of G that we need to exceed for the lower tail, while the upper critical value qU corresponds to the quantile that we need to be below for the upper tail, in order to include the candidate pseudo-VCATE value in the confidence interval.
#' @export
#'
#' @examples
crit_vals <- function(pVCATE, zeta, alpha = 0.05, Fn_hat_nk) {
  F0 <- Fn_hat_nk(0, pVCATE, zeta)
  qL <- min(alpha/2, F0)
  qU <- 1 - alpha + min(alpha/2, F0)
  c(qL = qL, qU = qU)
}



#' Title
#'
#' @param V_taunk_star_hat The estimated pseudo-VCATE value for the fold of interest, as returned by the pseudoVCATE_SB function. This is the point estimate of the pseudo-VCATE that we want to construct a confidence interval around. The adaptive confidence interval will be constructed by testing whether candidate pseudo-VCATE values are included in the acceptance region defined by the distribution of the pseudo-VCATE estimator, which depends on V_taunk_star_hat as well as on the variance matrix Omegank_hat and the sample size nk.
#' @param Omegank_hat The estimated asymptotic variance matrix of the pseudo-VCATE estimator, as returned by the pseudoVCATE_SB function. This should be a 2x2 matrix corresponding to the variance of the leading normal term and the variance of the quadratic term in the asymptotic distribution of the pseudo-VCATE estimator. This matrix is used to define the distribution of the random variable G, which is used in the construction of the adaptive confidence interval. The critical values for including candidate pseudo-VCATE values in the confidence interval will depend on this variance matrix, as it affects the shape of the distribution of G.
#' @param nk The sample size in the fold for which we want to construct the confidence interval. This is used to scale the terms in the distribution of the pseudo-VCATE estimator, as the distribution of the estimator depends on nk through the scaling of the leading normal term and the quadratic term. The critical values for the adaptive confidence interval will also depend on nk, as it affects the variability of the pseudo-VCATE estimator and therefore the quantiles of the distribution of G.
#' @param alpha The significance level for the confidence interval. This should be a numeric value between 0 and 1, typically something like 0.05 for a 95% confidence interval. The critical values computed by this function will be based on this alpha level, as they will determine the quantiles of the distribution of G that correspond to the desired coverage probability of the confidence interval. A smaller alpha will lead to wider confidence intervals, while a larger alpha will lead to narrower intervals.
#' @param Vmax An optional upper bound for the grid of candidate pseudo-VCATE values to test for inclusion in the confidence interval. If NULL, the function will choose a reasonable upper bound based on V_taunk_star_hat and the variance matrix Omegank_hat. The choice of Vmax can affect the computational time and the accuracy of the confidence interval, as it determines how far out we search for candidate pseudo-VCATE values. If Vmax is too small, we might miss including values that should be in the confidence interval; if it's too large, we might spend unnecessary time testing values that are very unlikely to be included.
#' @param gridN The number of grid points to use when testing candidate pseudo-VCATE values for inclusion in the confidence interval. A larger gridN will give a more accurate approximation of the confidence interval but will take more computational time. The function will create a grid of candidate pseudo-VCATE values between 0 and Vmax, and will test each one for inclusion in the confidence interval based on the distribution of G. The choice of gridN can affect the precision of the confidence interval, as a finer grid allows for a more accurate determination of the interval endpoints.
#' @param M The number of Monte Carlo samples to draw for approximating the distribution of the pseudo-VCATE estimator when computing the empirical CDF Fn_hat. A larger M will give a more accurate approximation of the distribution of G, which in turn will lead to more accurate critical values and confidence intervals, but it will also take more computational time. The choice of M can affect the stability and accuracy of the confidence interval, especially if the distribution of G is complex or if we are looking at extreme quantiles.
#' @param seed Random seed for reproducibility of the Monte Carlo draws used in approximating the distribution of G. Setting a seed allows for consistent results across runs, which can be important for debugging and for comparing results. The choice of seed does not affect the theoretical properties of the confidence interval, but it can affect the specific numerical results due to the randomness in the Monte Carlo approximation.
#'
#' @returns
#' @export
#'
#' @examples
adaptive_vcate_ci <- function(V_taunk_star_hat, Omegank_hat, nk,
                              alpha = 0.05,
                              Vmax = NULL, gridN = 400,
                              M = 4000, seed = 1) {
  
  if (is.na(Omegank_hat[1,1]) || is.na(Omegank_hat[2,2])) {
    return(list(
      ci = c(0, 0)
    ))
  }
  
  stopifnot(is.matrix(Omegank_hat), all(dim(Omegank_hat) == c(2,2)))
  stopifnot(nk > 10, V_taunk_star_hat >= 0, alpha > 0, alpha < 1)
  
  # Choose a reasonable upper grid bound if not provided
  if (is.null(Vmax)) {
    Vmax <- max(1e-8, V_taunk_star_hat + 10 * sqrt(Omegank_hat[2,2])^2) 
    Vmax <- max(Vmax, 5 * max(V_taunk_star_hat, 1e-8))
  }
  
  # Empirical CDF of G at value 'v', for given pVCATE and zeta
  Fn_hat_nk <- Fn_hat(Omegank_hat = Omegank_hat, nk = nk, M = M, seed = seed)
  
  # Membership test: does candidate pVCATE belong to CI for at least one zeta?
  in_CI <- function(pVCATE) {
    delta <- V_taunk_star_hat - pVCATE
    ok <- FALSE
    for (zeta in c(-1, 1)) {
      pval_like <- Fn_hat_nk(v = delta, pVCATE = pVCATE, zeta)
      qs <- crit_vals(pVCATE = pVCATE, zeta = zeta, alpha = alpha, Fn_hat_nk = Fn_hat_nk)
      if (pval_like >= qs["qL"] && pval_like <= qs["qU"]) {
        ok <- TRUE
        break
      }
    }
    ok
  }

  # Test a grid of candidate pVCATE values for inclusion in the CI
  Vmax_local <- max(0.5, 4 * V_taunk_star_hat)                   
  Vgrid1 <- seq(0, Vmax_local, by = max(1e-5, Vmax_local/1000))
  
  # Optional: add a coarse tail if you want to search beyond Vmax_local
  Vgrid2 <- seq(Vmax_local, max(Vmax_local, Vmax), length.out = 200)
  
  Vgrid <- unique(c(Vgrid1, Vgrid2))
  inside <- vapply(Vgrid, in_CI, logical(1))
  
  if (!any(inside)) {
    # If nothing included, expand search once (common when Vmax too small)
    Vmax2 <- 2 * Vmax
    Vgrid <- seq(0, Vmax2, length.out = gridN)
    inside <- vapply(Vgrid, in_CI, logical(1))
  }
  
  if (!any(inside)) {
    return(list(
      ci = c(NA_real_, NA_real_),
      Vgrid = Vgrid,
      inside = inside,
      note = "No grid points satisfied the inversion; increase Vmax/gridN/M."
    ))
  }
  
  # Typically the acceptance region is an interval; take min/max of accepted grid points
  ci_low  <- min(Vgrid[inside])
  ci_high <- max(Vgrid[inside])
  
  list(
    ci = c(ci_low, ci_high)
  )
}








crossfit_pseudoVCATE <- function(data, hal = FALSE,
                                 crossfit_folds = 5,
                                 SL.lib = c("SL.glmnet", "SL.randomForest", "SL.xgboost"),
                                 family = gaussian(),
                                 randomized = TRUE,
                                 pi_hat = 0.5,
                                 SL_cv_folds = 5,
                                 alpha = 0.05,
                                 Vmax = NULL,
                                 gridN = 200,
                                 M = 2000,
                                 seed = NA,
                                 adaptive_ci = TRUE) {
  
  n <- nrow(data)
  
  if (!is.na(seed)) {
    set.seed(seed)
  }
  
  # Random fold allocation
  fold_id <- sample(rep(seq_len(crossfit_folds), length.out = n))
  
  fold_results <- vector("list", crossfit_folds)
  
  for (k in seq_len(crossfit_folds)) {
    eval_idx <- which(fold_id == k)
    train_idx <- setdiff(seq_len(n), eval_idx)
    
    train_data <- data[train_idx, , drop = FALSE]
    eval_data  <- data[eval_idx, , drop = FALSE]
    
    nk <- nrow(eval_data)
    
    nuisance_k <- nuisance_param_func(
      train_data = train_data,
      eval_data = eval_data,
      SL.lib = SL.lib,
      randomized = randomized,
      pi_hat = pi_hat,
      K = SL_cv_folds, 
      family = family, 
      hal = hal
    )
    
    pseudo_k <- pseudoVCATE_SB(
      S_nk = nuisance_k$S_nk,
      M_nk = nuisance_k$M_nk,
      pi_hat = nuisance_k$pi_hat,
      data = eval_data
    )
    
    if (adaptive_ci) {
      ci_k <- adaptive_vcate_ci(
        V_taunk_star_hat = pseudo_k$V_taunk_star_hat,
        Omegank_hat = pseudo_k$Omegank_hat,
        nk = nk,
        alpha = alpha,
        Vmax = Vmax,
        gridN = gridN,
        M = M,
        seed = seed
      )
    } else {
      ci_k <- list(ci = c(NA_real_, NA_real_))
    }

    
    fold_results[[k]] <- as.data.frame(as.list(c(
      fold = k,
      nk = nk,
      V_taunk_star_hat = pseudo_k$V_taunk_star_hat,
      setNames(as.vector(pseudo_k$Omegank_hat), paste0("Omega_", 1:4)),
      setNames(pseudo_k$thetank_hat, paste0("theta_", 1:4)),
      setNames(ci_k$ci, c("ci_low", "ci_high"))
    )))
  }
  
  summary_df <-  dplyr::bind_rows(fold_results)
  return(summary_df)
}








