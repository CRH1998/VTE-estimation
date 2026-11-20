#!/usr/local/bin/Rscript

require(hal9001)
require(mvtnorm)
require(CVXR)
require(kernlab)
require(MASS)
# require(splines2)

#' Constructs basis functions for linear thershold rules with one-dimensional
#' effect modifier.
#'
#' 
#' @param x Effect modifier of interest, an n-dimensional vector.
#' @param no.tresholds Number of evenly-spaced thresholds. Set to sqrt(n) by default.
#'                     Unused if vector of thresholds are supplied by user.
#' @param thresholds User-supplied thresholds for linear threshold rule. 
#' 
#' @return An n x (2 * no.threshold matrix) containing basis functions for linear
#'         threshold rules.
#'                  
#' 
#' @export
ThresholdBasis <- function(x, no.thresholds = ceiling(sqrt(length(x))),
                           thresholds = NULL) {
  n <- length(x)
  
  if(is.null(thresholds)) {
      thresholds <- min(x) + ((max(x) - min(x))/(no.thresholds + 1)) * (1:no.thresholds)
  } else {
    no.thresholds <- length(thresholds)
  }
  
  direc.basis <- matrix(NA, nrow = n, ncol = 2 * no.thresholds)
  for(i in 1:no.thresholds) {
    direc.basis[,i] <- ifelse(x <= thresholds[i], 1, 0)
    direc.basis[,no.thresholds+i] <- ifelse(x > thresholds[i], 1, 0)
  }
  
  return(direc.basis)
}

#' Constructs basis functions for linear threshold rules with multi-dimensional
#' effect modifier.
#'
#' 
#' @param x Effect modifier of interest, an n x d matrix.
#' @param no.thresholds Number of evenly-spaced thresholds b for 
#'                      ifelse(1, 0, sum(weigths * x[,j] > b) )
#' @param no.thresholds.ang Number of evenly-spaced grid points for specifying weights
#'                     
#' 
#' @return An n x (2 * no.threshold matrix) containing basis functions for linear
#'         threshold rules.
#'                  
#' 
#' @export
ThresholdBasisMulti <- function(x, no.thresholds.ang = 10, no.thresholds = 10) {
  n <- nrow(x)
  
  x.scale <- matrix(NA, nrow = n, ncol = ncol(x))
  for(j in 1:ncol(x)) {
    x.scale[,j] <- (x[,j] - min(x[,j]))/(max(x[,j]) - min(x[,j]))
  }
  
  ang <- expand.grid(rep(list(seq(0, pi/2,  length.out = no.thresholds.ang)), ncol(x) - 1))

  weights <- matrix(NA, nrow = nrow(ang), ncol = ncol(x))
  weights[,1] <- cos(ang[,1])
  for(j in 2:ncol(weights)) {
    if(j < ncol(x)) {
      weights[,j] <- apply(as.matrix(sin(ang[,which(1:ncol(ang) < j)])), 1, prod) *
                     cos(ang[, j])
    } else {
      weights[,j] <- apply(as.matrix(sin(ang[,which(1:ncol(ang) < j)])), 1, prod) *
                     sin(ang[,j-1])
    }
  }
  weights <- diag(1/apply(weights, 1, sum)) %*% weights
  
  thresholds <- seq(.05, .95, length.out = no.thresholds) 
  direc.basis <- matrix(NA, nrow = n, ncol = nrow(weights) * length(thresholds))
  
  for(i in 1:nrow(weights)) {
    for(j in 1:length(thresholds)) {
      z <- x.scale %*% c(weights[i,])
      direc.basis[,length(thresholds) * (i - 1) + j] <- ifelse(z <= thresholds[j], 1, 0)
    }
  }
  direc.basis <- cbind(direc.basis, 1 - direc.basis)
  
  return(direc.basis)
}

#' Constructs basis functions for linear thershold rules with one-dimensional
#' effect modifier.
#'
#' 
#' @param x Effect modifier of interest, an n-dimensional vector.
#' @param no.tresholds Number of evenly-spaced thresholds. Set to sqrt(n) by default.
#'                     Unused if vector of thresholds are supplied by user.
#' @param thresholds User-supplied thresholds for interval basis functions.
#' 
#' @return An (n x no.threshold) matrix containing basis functions for total variation class.
#'                  
#' 
#' @export
IntervalBasis <- function(x, no.thresholds = ceiling(sqrt(length(x))),
                          thresholds = NULL) {
  n <- length(x)
  
  if(is.null(thresholds)) {
      thresholds <- quantile(x, (1:no.thresholds)/no.thresholds)
  } else {
    no.thresholds <- length(thresholds)
  }
  
  interval.basis <- matrix(NA, nrow = n, ncol = no.thresholds)
  interval.basis[,1] <- ifelse(x <= thresholds[1], 1, 0)
  for(i in 1:(no.thresholds-1)) {
    interval.basis[,i+1] <- ifelse(x <= thresholds[i+1] & x > thresholds[i], 1, 0)
  }
    
  return(interval.basis)
}

#' Estimates nuisance parameters using the highly adaptive LASSO
#'
#' 
#' @param y Outcome variable, an n-dimesional vector
#' @param a Treatment variable, an n-dimensional vector.
#' @param w Covariates for adjustment, an n x p matrix.
#' @param randomized Boolean variable, indicating whether the treatment probability can be
#'                   estimated using sample proportion  
#' 
#' @return An n x (2 * no.threshold matrix) containing basis functions for linear
#'         threshold rules.
#'                  
#' @return A list containing the following: 
#'         \describe{
#'            \item{outcome.regression} Estimated conditional mean, fitted at observed a and w
#'            \item{outcome.regression.tx0} Estimated conditional mean under a = 0 and observed w
#'            \item{outcome.regression.tx1} Estimated conditional mean under a = 1 and observed w
#'            \item{propensity.score} Estimated conditonal treatment probability at observed w
#'          }
#' 
#' @export
EstimateNuisance <- function(y, a, w, randomized = FALSE, 
                             hal = TRUE, SL.lib = c("SL.glmnet", "SL.randomForest", "SL.xgboost"),
                             family = gaussian(), K = 10) {
  n <- length(y)
  
  if (hal) {
     outcome.reg.fit <- fit_hal(X = cbind(a, w), Y = y, yolo = FALSE)
     
     outcome.regression <- predict(outcome.reg.fit, new_data = cbind(a, w))
     outcome.regression.tx0 <- predict(outcome.reg.fit, new_data = cbind(0, w))
     outcome.regression.tx1 <- predict(outcome.reg.fit, new_data = cbind(1, w))
     
  } else {
    # Creating the design matrix for SuperLearner, including interactions between treatment and covariates
    X_design <- cbind(
      a,
      w
      #,as.data.frame(a * as.matrix(w))
    )
    names(X_design) <- c("A", names(w))# , paste0("A_", names(w)))
    
    # Fitting SuperLearner
    outcome.reg.fit <- SuperLearner::SuperLearner(
      Y = y,
      X = X_design,
      family = family,
      SL.library = SL.lib,
      cvControl = list(V = K)
    )
    
    # For the evaluation data
    outcome.regression <- as.numeric(predict(outcome.reg.fit, new_data = X_design)$pred)
    outcome.regression.tx0 <- as.numeric(predict(outcome.reg.fit, newdata = X_design)$pred)
    outcome.regression.tx1 <- as.numeric(predict(outcome.reg.fit, newdata = X_design)$pred)
  }

  prop.score <- numeric(n)
  if(randomized) {
    prop.score <- rep(mean(a), n)
  } else {
    prop.score.fit <- fit_hal(X = w, Y = a, family = "binomial", yolo = FALSE)
    prop.score <- predict(prop.score.fit, new_data = w)
  }
  
  out <- list(prop.score = prop.score,
              outcome.regression = outcome.regression,
              outcome.regression.tx0 = outcome.regression.tx0,
              outcome.regression.tx1 = outcome.regression.tx1)
  return(out)
}

#' Estimates expected product of centered conditional treatment effect and basis functions
#
#' @param y Outcome variable, an n-dimesional vector
#' @param a Treatment variable, an n-dimensional vector.
#' @param w Covariates for adjustment, an n x p matrix.
#' @param outcome.regression Estimated conditional mean, fitted at observed a and w
#' @param outcome.regression.tx0 Estimated conditional mean under a = 0 and observed w
#' @param outcome.regression.tx1 Estimated conditional mean under a = 1 and observed w
#' @param propensity.score Estimated conditional treatment probability at observed w
#' @param direc.basis (n x no.threshold) matrix containing basis functions
#'    
#' @return A list containing the following: 
#'         \describe{
#'            \item{one.step} one-step estimator of the expected product of centered conditional treatment effect and basis functions
#'            \item{eif} efficient influence function for one-step estimator
#'          }
#'
#' @export
OneStepQuant <- function(y, a, w, 
                         prop.score, outcome.regression,
                         outcome.regression.tx0, outcome.regression.tx1,
                         direc.basis) {
  n <- length(y)
  
  pseudo.outcome <- (outcome.regression.tx1 - outcome.regression.tx0) +
                    (2*a - 1) * (y - outcome.regression)/
                    (prop.score^a * (1-prop.score)^(1-a))
  ate <- mean(pseudo.outcome)
  
  one.step <- t((diag(n) - matrix(1/n, n, n)) %*% (1 - 2 * direc.basis)) %*% 
              (pseudo.outcome - ate) * (1/n)
  
  eif <- diag(pseudo.outcome - ate) %*% ((diag(n) - matrix(1/n, n, n)) %*% (1 - 2 * direc.basis))
  eif <- (diag(n) - matrix(1/n, n, n)) %*% eif
  
  out <- list(one.step = one.step,
              eif = eif)
  
  return(out)
}

fast_one_step <- function(y, a, prop.score, outcome.regression,
                          outcome.regression.tx0, outcome.regression.tx1,
                          direc.basis) {
  n <- length(y)
  
  # Pseudo-outcome
  pseudo.outcome <- (outcome.regression.tx1 - outcome.regression.tx0) +
    (2 * a - 1) * (y - outcome.regression) /
    (prop.score^a * (1 - prop.score)^(1 - a))
  
  # Centered pseudo-outcome
  r <- pseudo.outcome - mean(pseudo.outcome)
  
  # Basis transform
  B <- 1 - 2 * direc.basis
  
  # (I - J) %*% B  == column-center B
  Bc <- sweep(B, 2, colMeans(B), FUN = "-")
  
  # one.step = t(Bc) %*% r / n
  one.step <- crossprod(Bc, r) / n
  
  # eif = (I - J) %*% diag(r) %*% Bc
  # First multiply each row of Bc by r
  eif <- Bc * r
  
  # Then left-multiply by (I - J): center each column
  eif <- sweep(eif, 2, colMeans(eif), FUN = "-")
  
  list(
    one.step = one.step,
    eif = eif
  )
}

#' Estimates expected product of conditional treatment effect and basis functions
#
#' @param y Outcome variable, an n-dimesional vector
#' @param a Treatment variable, an n-dimensional vector.
#' @param w Covariates for adjustment, an n x p matrix.
#' @param outcome.regression Estimated conditional mean, fitted at observed a and w
#' @param outcome.regression.tx0 Estimated conditional mean under a = 0 and observed w
#' @param outcome.regression.tx1 Estimated conditional mean under a = 1 and observed w
#' @param propensity.score Estimated conditional treatment probability at observed w
#' @param direc.basis (n x no.threshold) matrix containing basis functions
#' @param delta threshold for CATE for qualitative heterogeneity
#'    
#' @return A list containing the following: 
#'         \describe{
#'            \item{one.step} one-step estimator of the expected product of centered conditional treatment effect and basis functions
#'            \item{eif} efficient influence function for one-step estimator
#'          }
#'
#' @export
OneStepQual <- function(y, a, w, 
                         prop.score, outcome.regression,
                         outcome.regression.tx0, outcome.regression.tx1,
                         direc.basis, delta = 0) {
  n <- length(y)
  
  pseudo.outcome <- (outcome.regression.tx1 - outcome.regression.tx0) +
                    (2*a - 1) * (y - outcome.regression)/
                    (prop.score^a * (1-prop.score)^(1-a))
  
  one.step.pos <- t(direc.basis) %*% (pseudo.outcome - delta) * (1/n)
  one.step.neg <- t(1 - direc.basis) %*% (pseudo.outcome - delta) * (1/n)
  
  eif.pos <- diag(pseudo.outcome - delta) %*% direc.basis
  eif.pos <- (diag(n) - matrix(1/n, n, n)) %*% eif.pos
  eif.neg <- diag(pseudo.outcome - delta) %*% (1 - direc.basis)
  eif.neg <- (diag(n) - matrix(1/n, n, n)) %*% eif.neg
  
  out <- list(one.step.pos = one.step.pos,
              one.step.neg = one.step.neg,
              eif.pos = eif.pos,
              eif.neg = eif.neg)
  
  return(out)
}

#' Calculates test statistic for quantitative null
#' 
#' @param one.step one-step estimator of the expected product of centered conditional treatment effect and basis functions
#' @param method either "threshold" for linear threshold rule statistic or "interval" for total variation statistic
#' @param no.roots postulated number of roots of centered conditional treatment effect curve
#' @param lambda postulated upper bound on variation norm of optimal decision rule. Overrides no.roots if supplied.
#' 
#' @return scalar test statistic
QuantTestStat <- function(one.step,
                          method = "threshold",
                          no.roots = 2, lambda = NULL) {
  test.stat <- NULL
  if(method == "threshold") {
    test.stat <- max(abs(one.step))
  } else if(method == "interval") {
    if(is.null(lambda)) {
      lambda <- 2^no.roots
    }
    d <- length(one.step)
    D <- matrix(0, nrow = d - 1, ncol = d)
    D[cbind(1:(d-1), 1:(d-1))] <- 1
    D[cbind(1:(d-1), 2:d)] <- -1
    coef <- Variable(d)
    objective <- Minimize(t(one.step) %*% coef)
    prob <- Problem(objective,
                    list(p_norm(coef, "inf") <= 1,
                         p_norm(D %*% coef, 1) <= lambda))
    CVXR.result <- solve(prob, num_iter = 5000, 
                         reltol = 1e-02, feastol = 1e-02, abstol = 1e-02)
    test.stat <- ifelse(CVXR.result$status == "solver_error",
                          NA, -CVXR.result$value)
  }  else {
    stop("direction must be one of...")
  }
  
  out <- test.stat
  return(out)
}

#' Calculates test statistic for quantitative null
#' 
#' @param one.step one-step estimator of the expected product of delta-centered conditional treatment effect and basis functions
#' @param method either "threshold" for linear threshold rule statistic or "interval" for total variation statistic
#' @param no.roots postulated number of roots of centered conditional treatment effect curve
#' @param lambda postulated upper bound on variation norm of optimal decision rule. Overrides no.roots if supplied.
#' 
#' @return scalar test statistic
#' @return A list containing the following: 
#'         \describe{
#'            \item{test.stat.pos} test-statistic for one-sided test of area above the delta-centered conditional treatment curve 
#'            \item{test.stat.neg} test-statistic for one-sided test of area below the delta-centered conditional treatment curve
#'          }
QualTestStat <- function(one.step.pos, one.step.neg,
                          method = "threshold",
                          no.roots = 2, lambda = NULL) {
  test.stat <- NULL
  if(method == "threshold") {
    test.stat.pos <- max(one.step.pos)
    test.stat.neg <- min(one.step.neg)
  } else if(method == "interval") {
    if(is.null(lambda)) {
      lambda <- 2^no.roots
    }
    d <- length(one.step.pos)
    D <- matrix(0, nrow = d - 1, ncol = d)
    D[cbind(1:(d-1), 1:(d-1))] <- 1
    D[cbind(1:(d-1), 2:d)] <- -1
    coef <- Variable(d)
    objective.pos <- Maximize(t(one.step.pos) %*% coef)
    objective.neg <- Minimize(t(one.step.pos) %*% coef)
    prob.pos <- Problem(objective.pos,
                    list(p_norm(coef, "inf") <= 1,
                         p_norm(D %*% coef, 1) <= lambda))
    prob.pos <- Problem(objective.pos,
                        list(max_entries(coef) <= 1,
                             min_entries(coef) >= 0,
                             p_norm(D %*% coef, 1) <= lambda))
    prob.neg <- Problem(objective.neg,
                    list(p_norm(coef, "inf") <= 1,
                         p_norm(D %*% coef, 1) <= lambda))
    prob.neg <- Problem(objective.neg,
                        list(max_entries(coef) <= 1,
                             min_entries(coef) >= 0,
                             p_norm(D %*% coef, 1) <= lambda))
    CVXR.result.pos <- solve(prob.pos, num_iter = 5000, 
                             reltol = 1e-02, feastol = 1e-02, abstol = 1e-02)
    CVXR.result.neg <- solve(prob.neg, num_iter = 5000, 
                             reltol = 1e-02, feastol = 1e-02, abstol = 1e-02)
    test.stat.pos <- ifelse(CVXR.result.pos$status == "solver_error",
                            NA, CVXR.result.pos$value)
    test.stat.neg <- ifelse(CVXR.result.neg$status == "solver_error",
                            NA, CVXR.result.neg$value)
  }  else {
    stop("direction must be one of...")
  }
  
  out <- list(test.stat.pos = test.stat.pos,
              test.stat.neg = test.stat.neg)
}

#' Test for quantitative heterogeneity with one-diemsnional effect modifier
#
#' @param y Outcome variable, an n-dimesional vector
#' @param a Treatment variable, an n-dimensional vector.
#' @param w Covariates for adjustment, an n x p matrix.
#' @param x Effect modifier of interest, an n-dimensional vector
#' @param method either "threshold" for linear threshold rule statistic or "interval" for total variation statistic
#' @param no.tresholds Number of evenly-spaced thresholds. Set to sqrt(n) by default.
#'                     Unused if vector of thresholds are supplied by user.
#' @param thresholds User-supplied thresholds for linear threshold rule or interval basis. 
#' @param no.roots postulated number of roots of centered conditional treatment effect curve
#' @param lambda postulated upper bound on variation norm of optimal decision rule. Overrides no.roots if supplied.
#' @param randomized Boolean variable, indicating whether the treatment probability can be
#'                   estimated using sample proportion  
#' @param no.bs number of multiplier bootstrap samples
#' @param nuisance.est list of nuisance parameter estimators, obtained from EstimateNuisance
#'    
#' @return A list containing the following: 
#'         \describe{
#'            \item{test.stat} observed test statistic
#'            \item{test.stats.bs} multiplier bootstrap test statistics
#'            \item{p.val} p-value
#'          }
#'
#' @export
QuantTest <- function(y, a, w, x,
                      method = "threshold",
                      no.thresholds = ceiling(sqrt(length(x))),
                      thresholds = NULL,
                      no.roots = 2, lambda = NULL,
                      randomized = FALSE, no.bs = 5000,
                      nuisance.est = NULL) {
  n <- length(y)
  # x can be contained within w, though it need not be 
  
  # First, create basis
  if(method == "threshold") {
    direc.basis <- ThresholdBasis(x, no.thresholds = no.thresholds,
                                  thresholds = NULL)
  } else if(method == "interval") {
    direc.basis <- IntervalBasis(x, no.thresholds = no.thresholds,
                                 thresholds = NULL)
  }
  
  # Second, one-step estimator for each basis
  if(is.null(nuisance.est)) {
    nuisance.est <- EstimateNuisance(y = y, a = a, w = w, 
                                     randomized = randomized, K = K, 
                                     hal = hal, SL.lib = SL.lib, family = family)
  }
  one.step <- OneStepQuant(y = y, a = a, w = w, direc.basis = direc.basis,
                           prop.score = nuisance.est$prop.score,
                           outcome.regression = nuisance.est$outcome.regression,
                           outcome.regression.tx0 = nuisance.est$outcome.regression.tx0, 
                           outcome.regression.tx1 = nuisance.est$outcome.regression.tx1)
  
  # Third, calculate test statistic on observed data
  test.stat <- QuantTestStat(sqrt(n) * one.step$one.step,
                             method = method,
                             no.roots = 2, lambda = NULL)
  
  # Fourth, multiplier bootstrap to approximate the limiting distribution
  test.stats.bs <- numeric(no.bs)
  for(i in 1:no.bs) {
    rad <- 2 * rbinom(n, 1, .5) - 1
    # rad <- rnorm(n)
    rad <- rad - mean(rad)
    test.stats.bs[i] <- QuantTestStat(t(one.step$eif) %*% (rad - mean(rad))/sqrt(n),
                                     method = method,
                                     no.roots = no.roots, lambda = lambda)
  }
  
  # Fifth Calculate p-value
  p.val <- mean(test.stats.bs > test.stat, na.rm = TRUE)
  
  # Last, return output
  out <- list(test.stat = test.stat,
              test.stats.bs = test.stats.bs,
              p.val = p.val)
  
  return(out)
}

#' Test for quantitative heterogeneity with multi-diemsnional effect modifier using linear threshold rules
#'
#' @param y Outcome variable, an n-dimesional vector
#' @param a Treatment variable, an n-dimensional vector.
#' @param w Covariates for adjustment, an n x p matrix.
#' @param x Effect modifier of interest, an n x d matrix
#' @param no.thresholds Number of evenly-spaced thresholds b for 
#'                      ifelse(1, 0, sum(weigths * x[,j] > b) )
#' @param no.thresholds.ang Number of evenly-spaced grid points for specifying weights
#' @param randomized Boolean variable, indicating whether the treatment probability can be
#'                   estimated using sample proportion  
#' @param no.bs number of multiplier bootstrap samples
#' @param nuisance.est list of nuisance parameter estimators, obtained from EstimateNuisance
#'    
#' @return A list containing the following: 
#'         \describe{
#'            \item{test.stat} observed test statistic
#'            \item{test.stats.bs} multiplier bootstrap test statistics
#'            \item{p.val} p-value
#'          }
#'
#' @export
QuantTestMulti <- function(y, a, w, x,
                           no.thresholds = 10,
                           no.thresholds.ang = 10,
                           randomized = FALSE, no.bs = 5000,
                           nuisance.est = NULL, 
                           SL.lib = c("SL.glmnet", "SL.randomForest", "SL.xgboost"),
                           family = gaussian(), K = 10,
                           hal = T) {
  n <- length(y)
  # x can be contained within w, though it need not be 
  #browser()
  # First, create basis
  direc.basis <- ThresholdBasisMulti(x, no.thresholds = no.thresholds, 
                                     no.thresholds.ang = no.thresholds.ang)
  
  # Second, one-step estimator for each basis
  if(is.null(nuisance.est)) {
    nuisance.est <- EstimateNuisance(y = y, a = a, w = w, randomized = randomized, hal = hal, SL.lib = SL.lib)
  }
  
  one.step <- fast_one_step(y = y, a = a,
                            prop.score = nuisance.est$prop.score,
                            outcome.regression = nuisance.est$outcome.regression,
                            outcome.regression.tx0 = nuisance.est$outcome.regression.tx0, 
                            outcome.regression.tx1 = nuisance.est$outcome.regression.tx1,
                            direc.basis = direc.basis)
    
    
  # one.step <- OneStepQuant(y = y, a = a, w = w, direc.basis = direc.basis,
  #                          prop.score = nuisance.est$prop.score,
  #                          outcome.regression = nuisance.est$outcome.regression,
  #                          outcome.regression.tx0 = nuisance.est$outcome.regression.tx0,
  #                          outcome.regression.tx1 = nuisance.est$outcome.regression.tx1)
  
  # Third, calculate test statistic on observed data
  test.stat <- QuantTestStat(sqrt(n) * one.step$one.step,
                             method = "threshold")
  
  # Fourth, multiplier bootstrap to approximate the limiting distribution
  test.stats.bs <- numeric(no.bs)
  for(i in 1:no.bs) {
    rad <- 2 * rbinom(n, 1, .5) - 1
    # rad <- rnorm(n)
    rad <- rad - mean(rad)
    test.stats.bs[i] <- QuantTestStat(t(one.step$eif) %*% (rad - mean(rad))/sqrt(n),
                                     method = "threshold")
  }
  
  # Fifth Calculate p-value
  p.val <- mean(test.stats.bs > test.stat, na.rm = TRUE)
  
  # Last, return output
  out <- list(test.stat = test.stat,
              test.stats.bs = test.stats.bs,
              p.val = p.val)
  
  return(out)
}

#' Test for qualitative heterogeneity with one-dimesnional effect modifier
#
#' @param y Outcome variable, an n-dimensional vector
#' @param a Treatment variable, an n-dimensional vector.
#' @param w Covariates for adjustment, an n x p matrix.
#' @param x Effect modifier of interest, an n-dimensional vector
#' @param delta Threshold for qualitative effect heterogeneity
#' @param method either "threshold" for linear threshold rule statistic or "interval" for total variation statistic
#' @param no.tresholds Number of evenly-spaced thresholds. Set to sqrt(n) by default.
#'                     Unused if vector of thresholds are supplied by user.
#' @param thresholds User-supplied thresholds for linear threshold rule or interval basis. 
#' @param no.roots postulated number of roots of centered conditional treatment effect curve
#' @param lambda postulated upper bound on variation norm of optimal decision rule. Overrides no.roots if supplied.
#' @param randomized Boolean variable, indicating whether the treatment probability can be
#'                   estimated using sample proportion  
#' @param no.bs number of multiplier bootstrap samples
#' @param nuisance.est list of nuisance parameter estimates, obtained from EstimateNuisance
#'    
#' @return A list containing the following: 
#'         \describe{
#'            \item{test.stat} observed test statistic for qualitative heterogeneity
#'            \item{test.stats.pos.bs} bootstrap test statitics for one-sided positive test
#'            \item{test.stats.neg.bs} bootstrap test statistic for one-sided negative test
#'            \item{p.val.pos} p-value for one-sided positive test
#'            \item{p.val.neg} p-value for one-sided negative test
#'            \item{p.val} overall p-value
#'          }
#'
#' @export
QualTest <- function(y, a, w, x,
                     delta = 0,
                     method = "threshold",
                     no.thresholds = ceiling(sqrt(length(x))),
                     thresholds = NULL,
                     no.roots = 2, lambda = NULL,
                     randomized = FALSE, no.bs = 5000,
                     nuisance.est = NULL) {
  n <- length(y)
  # x can be contained within w, though it need not be
  
  # First, create basis
  if(ncol(as.matrix(x)) == 1) {
    if(method == "threshold") {
      direc.basis <- ThresholdBasis(x, no.thresholds = no.thresholds,
                                    thresholds = NULL)
    } else if(method == "interval") {
      direc.basis <- IntervalBasis(x, no.thresholds = no.thresholds,
                                   thresholds = NULL)
    }
  } else {
    direc.basis <- ThresholdBasisMulti(x, no.thresholds = no.thresholds)
  }
  
  # Second, one-step estimator for each basis
  if(is.null(nuisance.est)) {
    nuisance.est <- EstimateNuisance(y = y, a = a, w = w, randomized = randomized)
  }
  one.step <- OneStepQual(y = y, a = a, w = w,
                          delta = delta,
                          direc.basis = direc.basis,
                          prop.score = nuisance.est$prop.score,
                          outcome.regression = nuisance.est$outcome.regression,
                          outcome.regression.tx0 = nuisance.est$outcome.regression.tx0, 
                          outcome.regression.tx1 = nuisance.est$outcome.regression.tx1)
  
  # Third, calculate test statistic on observed data
  test.stat <- QualTestStat(one.step.pos = sqrt(n) * one.step$one.step.pos,
                            one.step.neg = sqrt(n) * one.step$one.step.neg,
                            method = method,
                            no.roots = 2, lambda = NULL)
   
  # Fourth, multiplier bootstrap to approximate the limiting distribution
  test.stats.pos.bs <- numeric(no.bs)
  test.stats.neg.bs <- numeric(no.bs)
  for(i in 1:no.bs) {
    rad <- 2 * rbinom(n, 1, .5) - 1
    # rad <- rnorm(n)
    rad <- rad - mean(rad)
    bs.out <- QualTestStat(one.step.pos = t(one.step$eif.pos) %*% (rad - mean(rad))/sqrt(n),
                           one.step.neg = t(one.step$eif.neg) %*% (rad - mean(rad))/sqrt(n),
                           method = method,
                           no.roots = no.roots, lambda = lambda)
    test.stats.pos.bs[i] <- bs.out$test.stat.pos
    test.stats.neg.bs[i] <- bs.out$test.stat.neg
  }
  
  # Fifth Calculate p-value
  p.val.pos <- mean(test.stats.pos.bs > test.stat$test.stat.pos, na.rm = TRUE)
  p.val.neg <- mean(test.stats.neg.bs < test.stat$test.stat.neg, na.rm = TRUE)
  
  # Last, return output
  out <- list(test.stat = test.stat,
              test.stats.pos.bs = test.stats.pos.bs,
              test.stats.neg.bs = test.stats.neg.bs,
              p.val.pos = p.val.pos,
              p.val.neg = p.val.neg,
              p.val = max(c(p.val.pos, p.val.neg)))
  
  return(out)
}

#' Test for qualitative heterogeneity with multi-dimesnional effect modifier using linear threhold rules
#
#' @param y Outcome variable, an n-dimensional vector
#' @param a Treatment variable, an n-dimensional vector.
#' @param w Covariates for adjustment, an n x p matrix.
#' @param x Effect modifier of interest, an n-dimensional vector
#' @param delta Threshold for qualitative effect heterogeneity
#' @param no.thresholds Number of evenly-spaced thresholds b for 
#'                      ifelse(1, 0, sum(weigths * x[,j] > b) )
#' @param no.thresholds.ang Number of evenly-spaced grid points for specifying weights
#' @param randomized Boolean variable, indicating whether the treatment probability can be
#'                   estimated using sample proportion  
#' @param no.bs number of multiplier bootstrap samples
#' @param nuisance.est list of nuisance parameter estimates, obtained from EstimateNuisance
#'    
#' @return A list containing the following: 
#'         \describe{
#'            \item{test.stat} observed test statistic for qualitative heterogeneity
#'            \item{test.stats.pos.bs} bootstrap test statitics for one-sided positive test
#'            \item{test.stats.neg.bs} bootstrap test statistic for one-sided negative test
#'            \item{p.val.pos} p-value for one-sided positive test
#'            \item{p.val.neg} p-value for one-sided negative test
#'            \item{p.val} overall p-value
#'          }
#'
#' @export
QualTestMulti <- function(y, a, w, x,
                          delta = 0,
                          no.thresholds = 10,
                          no.thresholds.ang = 10,
                          randomized = FALSE, no.bs = 5000,
                          nuisance.est = NULL) {
  n <- length(y)
  # x can be contained within w, though it need not be
  
  
  # First, create basis
  direc.basis <- ThresholdBasisMulti(x, no.thresholds = no.thresholds, 
                                     no.thresholds.ang = no.thresholds.ang)
  
  # Second, one-step estimator for each basis
  if(is.null(nuisance.est)) {
    nuisance.est <- EstimateNuisance(y = y, a = a, w = w, randomized = randomized)
  }
  one.step <- OneStepQual(y = y, a = a, w = w,
                          delta = delta,
                          direc.basis = direc.basis,
                          prop.score = nuisance.est$prop.score,
                          outcome.regression = nuisance.est$outcome.regression,
                          outcome.regression.tx0 = nuisance.est$outcome.regression.tx0, 
                          outcome.regression.tx1 = nuisance.est$outcome.regression.tx1)
  
  # Third, calculate test statistic on observed data
  test.stat <- QualTestStat(one.step.pos = sqrt(n) * one.step$one.step.pos,
                            one.step.neg = sqrt(n) * one.step$one.step.neg,
                            method = "threshold")
   
  # Fourth, multiplier bootstrap to approximate the limiting distribution
  test.stats.pos.bs <- numeric(no.bs)
  test.stats.neg.bs <- numeric(no.bs)
  for(i in 1:no.bs) {
    rad <- 2 * rbinom(n, 1, .5) - 1
    # rad <- rnorm(n)
    rad <- rad - mean(rad)
    bs.out <- QualTestStat(one.step.pos = t(one.step$eif.pos) %*% (rad - mean(rad))/sqrt(n),
                           one.step.neg = t(one.step$eif.neg) %*% (rad - mean(rad))/sqrt(n),
                           method = "threshold")
    test.stats.pos.bs[i] <- bs.out$test.stat.pos
    test.stats.neg.bs[i] <- bs.out$test.stat.neg
  }
  
  # Fifth Calculate p-value
  p.val.pos <- mean(test.stats.pos.bs > test.stat$test.stat.pos, na.rm = TRUE)
  p.val.neg <- mean(test.stats.neg.bs < test.stat$test.stat.neg, na.rm = TRUE)
  
  # Last, return output
  out <- list(test.stat = test.stat,
              test.stats.pos.bs = test.stats.pos.bs,
              test.stats.neg.bs = test.stats.neg.bs,
              p.val.pos = p.val.pos,
              p.val.neg = p.val.neg,
              p.val = max(c(p.val.pos, p.val.neg)))
  
  return(out)
}

# quantitative test with categorical effect modifier
QuantTestCategorical <- function(y, a, w, x,
                                 randomized = FALSE, no.bs = 5000,
                                 nuisance.est = NULL) {
  n <- length(y)
  # x can be contained within w, though it need not be
  
  # First, create basis
  x.unique <- unique(x)
  dmat <- matrix(NA, nrow = n, ncol = length(unique(x)))
  for(j in 1:length(unique(x))) {
    dmat[,j] <- ifelse(x == x.unique[j], 1/sum(x == x.unique[j]), 0)
  }
  
  # Second, one-step estimator for each basis
  if(is.null(nuisance.est)) {
    nuisance.est <- EstimateNuisance(y = y, a = a, w = w, randomized = randomized)
  }
  pseudo.outcome <- (nuisance.est$outcome.regression.tx1 - 
                     nuisance.est$outcome.regression.tx0) +
                    (2*a - 1) * (y - nuisance.est$outcome.regression)/
                    (nuisance.est$prop.score^a * (1 - nuisance.est$prop.score)^(1-a))
  
  ate <- mean(pseudo.outcome)
  cate <- c(t(dmat) %*% pseudo.outcome)
  
  eif.ate <- pseudo.outcome - mean(pseudo.outcome)
  eif.cate <- (diag(n) - matrix(1/n, n, n)) %*% diag(pseudo.outcome) %*% (n * dmat)
  
  # Third, calculate test statistic on observed data
  test.stat <- sum(abs(cate - ate))
  
  # Fourth, multiplier bootstrap to approximate the limiting distribution
  test.stats.bs <- numeric(no.bs)
  for(i in 1:no.bs) {
    rad <- 2 * rbinom(n, 1, .5) - 1
    # rad <- rnorm(n)
    rad <- rad - mean(rad)
    
    ate.bs <- mean(rad * eif.ate)
    cate.bs <- apply(diag(rad) %*% eif.cate, 2, mean)
    test.stats.bs[i] <- sum(abs(cate.bs - ate.bs))
  }
  
  # Fifth Calculate p-value
  p.val <- mean(test.stats.bs > test.stat, na.rm = TRUE)
  
  # Last, return output
  out <- list(test.stat = test.stat,
              test.stats.bs = test.stats.bs,
              p.val = p.val)
  
  return(out)
}

# qualitative test with categorical effect modifier, using either Gail-Simon test or range test
QualTestCategorical <- function(y, a, w, x,
                                randomized = FALSE, no.bs = 5000,
                                delta = 0, method = "Gail-Simon",
                                nuisance.est = NULL) {
  n <- length(y)
  # x can be contained within w, though it need not be
  
  # First, create basis
  x.unique <- unique(x)
  dmat <- matrix(NA, nrow = n, ncol = length(unique(x)))
  for(j in 1:length(unique(x))) {
    dmat[,j] <- ifelse(x == x.unique[j], 1/sum(x == x.unique[j]), 0)
  }
  
  # Second, one-step estimator for each basis
  if(is.null(nuisance.est)) {
    nuisance.est <- EstimateNuisance(y = y, a = a, w = w, randomized = randomized)
  }
  pseudo.outcome <- (nuisance.est$outcome.regression.tx1 - 
                     nuisance.est$outcome.regression.tx0) +
                    (2*a - 1) * (y - nuisance.est$outcome.regression)/
                    (nuisance.est$prop.score^a * (1 - nuisance.est$prop.score)^(1-a))
  
  cate <- c(t(dmat) %*% pseudo.outcome)
  eif.cate <- (diag(n) - matrix(1/n, n, n)) %*% diag(pseudo.outcome) %*% (n * dmat)
  se.cate <- apply(eif.cate, 2, sd)/sqrt(n)
  
  # Third, calculate test statistic on observed data
  if(method == "Gail-Simon") {
    test.stat.pos <- sum((cate-delta)^2/se.cate^2 * (cate > delta))
    test.stat.neg <- sum((cate-delta)^2/se.cate^2 * (cate < delta))
  } else {# method == "range"
    test.stat.pos <- max((cate-delta)^2/se.cate^2 * (cate > delta))
    test.stat.neg <- max((cate-delta)^2/se.cate^2 * (cate < delta))
  }
  
  # Fourth, p-value using the Gail-Simon Test
  norm.sample <- matrix(rnorm((length(x.unique) - 1) * no.bs, 0 , 1),
                      nrow = no.bs, ncol = length(x.unique) - 1)
  if(method == "Gail-Simon") {
    test.stat.bs <- apply(norm.sample, 1, function(x){sum(x^2 * ifelse(x > 0, 1, 0))})
  } else {# method == "range"
    test.stat.bs <- apply(norm.sample, 1, function(x){max(x^2 * ifelse(x > 0, 1, 0))})
  }
  p.val.pos <- mean(test.stat.bs > test.stat.pos)
  p.val.neg <- mean(test.stat.bs > test.stat.neg)
  # p.val.pos <- pchisq(test.stat.pos, df = length(unique(x)) - 1, lower.tail = FALSE)
  # p.val.neg <- pchisq(test.stat.neg, df = length(unique(x)) - 1, lower.tail = FALSE)
  p.val <- max(c(p.val.pos, p.val.neg))
  
  # Last, return output
  out <- list(test.stat.pos = test.stat.pos,
              test.stat.neg = test.stat.neg,
              p.val.pos = p.val.pos,
              p.val.neg = p.val.neg,
              p.val = p.val)
  
  return(out)
}

# Quantitative test using kernel ridge regression
QuantTestKRidge <- function(y, a, w, x, dim.basis = min(c(100, round(length(y)/2))),
                            randomized = FALSE, no.mc = 5000, no.folds = 10,
                            nuisance.est = NULL) {
  n <- length(y)
  # x can be contained within w, though it need not be
  
  # First, create basis
  x <- as.matrix(x)
  sd.x <- apply(x, 2, sd)
  if(ncol(x)  == 1) {
    sigma.inv <- 1/sd.x
  } else {
    sigma.inv <- diag(1/sd.x)
  }
  K <- kernelMatrix(rbfdot(sigma = 1), x = x %*% sigma.inv)
  eigen.K <- eigen(K, symmetric = TRUE)
  basis <- eigen.K$vectors[,1:dim.basis]
  pweights <- 1/eigen.K$values[1:dim.basis]
  
  # Second, one-step estimator for each basis
  if(is.null(nuisance.est)) {
    nuisance.est <- EstimateNuisance(y = y, a = a, w = w, randomized = randomized)
  }
  pseudo.outcome <- (nuisance.est$outcome.regression.tx1 - 
                     nuisance.est$outcome.regression.tx0) +
                    (2*a - 1) * (y - nuisance.est$outcome.regression)/
                    (nuisance.est$prop.score^a * (1 - nuisance.est$prop.score)^(1-a))
  
  ate <- mean(pseudo.outcome)
  one.step <- t(basis) %*% (pseudo.outcome - ate) * (1/n)
  eif <- diag(pseudo.outcome - ate) %*% basis
  
  # Third, calculate test statistic on observed data
  folds <- sample(1:no.folds, n, replace = TRUE)
  lam.max <- sum(diag(t(basis) %*% basis) * 1/pweights)/n
  lam.seq <- exp(seq(log(lam.max * .001/log(n)), log(lam.max), length.out = 50))
  pred.error <- matrix(NA, nrow = length(lam.seq), ncol = no.folds)
  for(j in 1:no.folds) {
    for(k in 1:length(lam.seq)) {
      contrast.j <- t(basis[folds != j,]) %*% 
                    (pseudo.outcome[folds != j] - mean(pseudo.outcome[folds != j])) * 
                    (1/sum(folds != j))
      eigen.j <- eigen(t(basis[folds != j,]) %*% basis[folds != j,]/sum(folds != j) + lam.seq[k] * diag(pweights), symmetric = TRUE)
      coef.j <- t(eigen.j$vectors) %*% diag(1/eigen.j$values) %*% eigen.j$vectors %*% contrast.j
      cate.hold.j <- basis[folds == j,] %*% coef.j
      pred.error[k, j] <- mean((pseudo.outcome[folds == j] - cate.hold.j)^2)
    }
  }

  lam.cv <- lam.seq[which.min(apply(pred.error, 2, mean))]
  eigen.full <- eigen(t(basis) %*% basis/n + lam.cv * diag(pweights), symmetric = TRUE)
  coef <- Re(t(eigen.full$vectors) %*% diag(1/eigen.full$values) %*% eigen.full$vectors) %*% one.step 
  fitted <- basis %*% coef
  eif.coef <- eif %*% Re(t(eigen.full$vectors) %*% diag(1/eigen.full$values) %*% eigen.full$vectors)
  eif.fitted <- eif %*% Re(t(eigen.full$vectors) %*% diag(1/eigen.full$values) %*% eigen.full$vectors) %*% t(basis)
  test.stat <- mean(abs(fitted - mean(fitted)))
  
  # Fourth, multiplier bootstrap to approximate the limiting distribution
  eif.fitted.cov <- cov(eif.fitted)/n
  fitted.mc <- rmvnorm(no.mc, mean = rep(0, n), sigma = eif.fitted.cov)
  test.stats.mc <- apply(fitted.mc, 1, function(x){mean(abs(x - mean(x)))})
  
  # Fifth, calculate p-value
  p.val <- mean(test.stats.mc > test.stat, na.rm = TRUE)
  
  # Last, return output
  out <- list(test.stat = test.stat,
              test.stats.mc = test.stats.mc,
              p.val = p.val)
  
  return(out)
  
}

# qualitative test using kernel ridge regression
QualTestKRidge <- function(y, a, w, x, 
                           dim.basis = min(c(100, round(length(y)/2))), no.folds = 10,
                           randomized = FALSE, no.mc = 5000, delta = 0,
                           nuisance.est = NULL) {
  n <- length(y)
  # x can be contained within w, though it need not be
  
  # First, create basis
  x <- as.matrix(x)
  sd.x <- apply(x, 2, sd)
  if(ncol(x)  == 1) {
    sigma.inv <- 1/sd.x
  } else {
    sigma.inv <- diag(1/sd.x)
  }
  K <- kernelMatrix(rbfdot(sigma = 1), x = x %*% sigma.inv)
  eigen.K <- eigen(K, symmetric = TRUE)
  basis <- eigen.K$vectors[,1:dim.basis]
  pweights <- 1/eigen.K$values[1:dim.basis]
  
  # Second, one-step estimator for each basis
  if(is.null(nuisance.est)) {
    nuisance.est <- EstimateNuisance(y = y, a = a, w = w, randomized = randomized)
  }
  pseudo.outcome <- (nuisance.est$outcome.regression.tx1 - 
                     nuisance.est$outcome.regression.tx0) +
                    (2*a - 1) * (y - nuisance.est$outcome.regression)/
                    (nuisance.est$prop.score^a * (1 - nuisance.est$prop.score)^(1-a))
  
  one.step <- t(basis) %*% pseudo.outcome * (1/n)
  eif <- (diag(n) - matrix(1/n, n, n)) %*% (diag(pseudo.outcome) %*% basis)
  
  # Third, calculate test statistic on observed data
  folds <- sample(1:no.folds, n, replace = TRUE)
  lam.max <- sum(diag(t(basis) %*% basis) / pweights)/n
  lam.seq <- exp(seq(log(lam.max * .001/log(n)), log(lam.max), length.out = 50))
  pred.error <- matrix(NA, nrow = length(lam.seq), ncol = no.folds)
  for(j in 1:no.folds) {
    for(k in 1:length(lam.seq)) {
      contrast.j <- t(basis[folds != j,]) %*% pseudo.outcome[folds != j] * (1/sum(folds != j))
      eigen.j <- eigen(t(basis[folds != j,]) %*% basis[folds != j,]/sum(folds != j) + lam.seq[k] * diag(pweights), symmetric = TRUE)
      coef.j <- t(eigen.j$vectors) %*% diag(1/eigen.j$values) %*% eigen.j$vectors %*% contrast.j
      cate.hold.j <- basis[folds == j,] %*% coef.j
      pred.error[k, j] <- mean((pseudo.outcome[folds == j] - cate.hold.j)^2)
    }
  }
  
  lam.cv <- lam.seq[which.min(apply(pred.error, 2, mean))]
  eigen.full <- eigen(t(basis) %*% basis/n + lam.cv * diag(pweights), symmetric = TRUE)
  coef <- t(eigen.full$vectors) %*% diag(1/eigen.full$values) %*% eigen.full$vectors %*% one.step 
  fitted <- basis %*% coef
  eif.fitted <- eif %*% t(eigen.full$vectors) %*% diag(1/eigen.full$values) %*% eigen.full$vectors %*% t(basis)

  test.stat.pos <- mean(abs(fitted * (fitted > delta)))
  test.stat.neg <- mean(abs(fitted * (fitted < delta)))
  
  # Fourth, multiplier bootstrap to approximate the limiting distribution
  eif.fitted.cov <- cov(eif.fitted)/n
  fitted.mc <- rmvnorm(no.mc, mean = rep(0, n), sigma = eif.fitted.cov)
  test.stats.pos.mc <- apply(fitted.mc, 1, function(x){mean(abs(x * (x > delta)))})
  test.stats.neg.mc <- apply(fitted.mc, 1, function(x){mean(abs(x * (x < delta)))})
  
  # Fifth, calculate p-value
  p.val.pos <- mean(test.stat.pos < test.stats.pos.mc, na.rm = TRUE)
  p.val.neg <- mean(test.stat.neg < test.stats.neg.mc, na.rm = TRUE)
  p.val <- max(c(p.val.pos, p.val.neg))
  
  # Last, return output
  out <- list(test.stat.pos = test.stat.pos,
              test.stat.neg = test.stat.neg,
              test.stats.pos.mc = test.stats.pos.mc,
              test.stats.neg.mc = test.stats.neg.mc,
              p.val.pos = p.val.pos,
              p.val.neg = p.val.neg,
              p.val = p.val)
  
  return(out)
  
}

QuantTestSplit <- function(y, a, w, x, split.prob = .5,
                           randomized = FALSE, dim.basis = min(c(100, round(length(y)/2))),
                           no.folds = 10) {
  
  n <- length(y)
  # x can be contained within w, though it need not be
  
  # First, create basis
  x <- as.matrix(x)
  sd.x <- apply(x, 2, sd)
  if(ncol(x)  == 1) {
    sigma.inv <- 1/sd.x
  } else {
    sigma.inv <- diag(1/sd.x)
  }
  K <- kernelMatrix(rbfdot(sigma = 1), x = x %*% sigma.inv)
  eigen.K <- eigen(K, symmetric = TRUE)
  basis <- eigen.K$vectors[,1:dim.basis]
  pweights <- 1/eigen.K$values[1:dim.basis]
  
  # Second, split the data in two partitions
  sample.split <- rbinom(n, size = 1, prob = split.prob)
  y1 <- y[sample.split == 1]
  y0 <- y[sample.split == 0]
  a1 <- a[sample.split == 1]
  a0 <- a[sample.split == 0]
  w1 <- as.matrix(w)[sample.split == 1,]
  w0 <- as.matrix(w)[sample.split == 0,]
  basis1 <- basis[sample.split == 1,]
  basis0 <- basis[sample.split == 0,]
  n1 <- sum(sample.split == 1)
  n0 <- sum(sample.split == 0)
  
  # Third, one-step estimator for each basis
  nuisance.est1 <- EstimateNuisance(y = y1, a = a1, w = w1, randomized = randomized)
  pseudo.outcome1 <- (nuisance.est1$outcome.regression.tx1 - 
                     nuisance.est1$outcome.regression.tx0) +
                     (2*a1 - 1) * (y1 - nuisance.est1$outcome.regression)/
                     (nuisance.est1$prop.score^a1 * (1 - nuisance.est1$prop.score)^(1-a1))
  
  nuisance.est0 <- EstimateNuisance(y = y0, a = a0, w = w0, randomized = randomized)
  pseudo.outcome0 <- (nuisance.est0$outcome.regression.tx1 - 
                      nuisance.est0$outcome.regression.tx0) +
                     (2*a0 - 1) * (y0 - nuisance.est0$outcome.regression)/
                     (nuisance.est0$prop.score^a0 * (1 - nuisance.est0$prop.score)^(1-a0))
  
  ate1 <- mean(pseudo.outcome1)
  ate0 <- mean(pseudo.outcome0)
  one.step1 <- t(basis1) %*% (pseudo.outcome1 - ate1) * (1/n1)
  one.step0 <- t(basis0) %*% (pseudo.outcome0 - ate0) * (1/n0)
  
  # Fourth, calculate test statistic on observed data
  lam.max <- sum(diag(t(basis) %*% basis)/pweights)/n
  lam.seq <- exp(seq(log(lam.max * .001/log(n)), log(lam.max), length.out = 50))
  
  folds1 <- sample(1:no.folds, n1, replace = TRUE)
  pred.error1 <- matrix(NA, nrow = length(lam.seq), ncol = no.folds)
  folds0 <- sample(1:no.folds, n0, replace = TRUE)
  pred.error0 <- matrix(NA, nrow = length(lam.seq), ncol = no.folds)
  for(j in 1:no.folds) {
    for(k in 1:length(lam.seq)) {
      contrast.j <- t(basis0[folds0 != j,]) %*% 
                    (pseudo.outcome0[folds0 != j] - mean(pseudo.outcome0[folds0 != j])) * 
                    (1/sum(folds0 != j))
      eigen.j <- eigen(t(basis0[folds0 != j,]) %*% basis0[folds0 != j,]/sum(folds0 != j) + lam.seq[k] * diag(pweights),
                       symmetric = TRUE)
      coef.j <- t(eigen.j$vectors) %*% diag(1/eigen.j$values) %*% eigen.j$vectors %*% contrast.j
      cate.hold.j <- basis0[folds0 == j,] %*% coef.j
      pred.error0[k, j] <- mean((pseudo.outcome0[folds0 == j] - cate.hold.j)^2)
    }
    for(k in 1:length(lam.seq)) {
      contrast.j <- t(basis1[folds1 != j,]) %*% 
                    (pseudo.outcome1[folds1 != j] - mean(pseudo.outcome1[folds1 != j])) * 
                    (1/sum(folds1 != j))
      eigen.j <- eigen(t(basis1[folds1 != j,]) %*% basis1[folds1 != j,]/sum(folds1 != j) + lam.seq[k] * diag(pweights), symmetric = TRUE)
      coef.j <- t(eigen.j$vectors) %*% diag(1/eigen.j$values) %*% eigen.j$vectors %*% contrast.j
      cate.hold.j <- basis1[folds1 == j,] %*% coef.j
      pred.error1[k, j] <- mean((pseudo.outcome1[folds1 == j] - cate.hold.j)^2)
    }
  }

  lam.cv0 <- lam.seq[which.min(apply(pred.error0, 2, mean))]
  lam.cv1 <- lam.seq[which.min(apply(pred.error1, 2, mean))]
  
  eigen.full1 <- eigen(t(basis1) %*% basis1/n1 + lam.cv1 * diag(pweights), symmetric = TRUE)
  coef1 <- t(eigen.full1$vectors) %*% diag(1/eigen.full1$values) %*% eigen.full1$vectors %*% one.step1 
  fitted0 <- basis0 %*% coef1
  
  eigen.full0 <- eigen(t(basis0) %*% basis0/n0 + lam.cv0 * diag(pweights), symmetric = TRUE)
  coef0 <- t(eigen.full0$vectors) %*% diag(1/eigen.full0$values) %*% eigen.full0$vectors %*% one.step0 
  fitted1 <- basis1 %*% coef0
  
  rule1 <- ifelse(fitted1 > mean(fitted1), 1, 0)
  rule0 <- ifelse(fitted0 > mean(fitted0), 1, 0)
  
  # Fifth, calculate p-value
  quant0 <- OneStepQuant(y = y0, a = a0, w = w0, 
                         prop.score = nuisance.est0$prop.score,
                         outcome.regression = nuisance.est0$outcome.regression,
                         outcome.regression.tx0 = nuisance.est0$outcome.regression.tx0,
                         outcome.regression.tx1 = nuisance.est0$outcome.regression.tx1,
                         direc.basis = rule0)
  p.val0 <- pchisq(n0 * quant0$one.step^2/var(quant0$eif), df = 1, lower.tail = FALSE)
  
  quant1 <- OneStepQuant(y = y1, a = a1, w = w1, 
                         prop.score = nuisance.est1$prop.score,
                         outcome.regression = nuisance.est1$outcome.regression,
                         outcome.regression.tx0 = nuisance.est1$outcome.regression.tx0,
                         outcome.regression.tx1 = nuisance.est1$outcome.regression.tx1,
                         direc.basis = rule1)
  p.val1 <- pchisq(n1 * quant1$one.step^2/var(quant1$eif), df = 1, lower.tail = FALSE)
  
  p.val.adj <- min(c(p.val0 * 2, p.val1 * 2, 1))
  
  # Last, return output
  out <- list(p.val0 = p.val0,
              p.val1 = p.val1,
              p.val.adj = p.val.adj)
  
  return(out)
}

QualTestSplit <- function(y, a, w, x, split.prob = .5,
                          randomized = FALSE, 
                          dim.basis = min(c(100, round(length(y)/2))), no.folds = 10,
                          delta = 0) {
  
  ## Estimate CATE using DR-learner and kernel ridge
  n <- length(y)
  # x can be contained within w, though it need not be
  
  # First, create basis
  x <- as.matrix(x)
  sd.x <- apply(x, 2, sd)
  if(ncol(x)  == 1) {
    sigma.inv <- 1/sd.x
  } else {
    sigma.inv <- diag(1/sd.x)
  }
  K <- kernelMatrix(rbfdot(sigma = 1), x = x %*% sigma.inv)
  eigen.K <- eigen(K, symmetric = TRUE)
  basis <- eigen.K$vectors[,1:dim.basis]
  pweights <- 1/eigen.K$values[1:dim.basis]
  
  # Second, split the data in two partitions
  sample.split <- rbinom(n, size = 1, prob = split.prob)
  y1 <- y[sample.split == 1]
  y0 <- y[sample.split == 0]
  a1 <- a[sample.split == 1]
  a0 <- a[sample.split == 0]
  w1 <- as.matrix(w)[sample.split == 1,]
  w0 <- as.matrix(w)[sample.split == 0,]
  basis1 <- basis[sample.split == 1,]
  basis0 <- basis[sample.split == 0,]
  n1 <- sum(sample.split == 1)
  n0 <- sum(sample.split == 0)
  
  # Second, one-step estimator for each basis
  nuisance.est1 <- EstimateNuisance(y = y1, a = a1, w = w1, randomized = randomized)
  pseudo.outcome1 <- (nuisance.est1$outcome.regression.tx1 - 
                     nuisance.est1$outcome.regression.tx0) +
                    (2*a1 - 1) * (y1 - nuisance.est1$outcome.regression)/
                    (nuisance.est1$prop.score^a1 * (1 - nuisance.est1$prop.score)^(1-a1))
  
  nuisance.est0 <- EstimateNuisance(y = y0, a = a0, w = w0, randomized = randomized)
  pseudo.outcome0 <- (nuisance.est0$outcome.regression.tx1 - 
                     nuisance.est0$outcome.regression.tx0) +
                    (2*a0 - 1) * (y0 - nuisance.est0$outcome.regression)/
                    (nuisance.est0$prop.score^a0 * (1 - nuisance.est0$prop.score)^(1-a0))
  
  one.step1 <- t(basis1) %*% pseudo.outcome1 * (1/n1)
  eif1 <- diag(pseudo.outcome1) %*% basis1
  one.step0 <- t(basis0) %*% pseudo.outcome0 * (1/n0)
  eif0 <- diag(pseudo.outcome0) %*% basis0
  
  # Third, estimate the optimal rule on each data set
  folds1 <- sample(1:no.folds, n1, replace = TRUE)
  folds0 <- sample(1:no.folds, n0, replace = TRUE)
  lam.max <- sum(diag(t(basis) %*% basis) / pweights)/n
  lam.seq <- exp(seq(log(lam.max * .001/log(n)), log(lam.max), length.out = 50))
  pred.error1 <- matrix(NA, nrow = length(lam.seq), ncol = no.folds)
  pred.error0 <- matrix(NA, nrow = length(lam.seq), ncol = no.folds)

  for(j in 1:no.folds) {
    for(k in 1:length(lam.seq)) {
      contrast.j <- t(basis0[folds0 != j,]) %*% 
                    (pseudo.outcome0[folds0 != j] - mean(pseudo.outcome0[folds0 != j])) * 
                    (1/sum(folds0 != j))
      eigen.j <- eigen(t(basis0[folds0 != j,]) %*% basis0[folds0 != j,]/sum(folds0 != j) + lam.seq[k] * diag(pweights), symmetric = TRUE)
      coef.j <- t(eigen.j$vectors) %*% diag(1/eigen.j$values) %*% eigen.j$vectors %*% contrast.j
      cate.hold.j <- basis0[folds0 == j,] %*% coef.j
      pred.error0[k, j] <- sum((pseudo.outcome0[folds0 == j] - cate.hold.j)^2)
    }
    for(k in 1:length(lam.seq)) {
      contrast.j <- t(basis1[folds1 != j,]) %*% 
                    (pseudo.outcome1[folds1 != j] - mean(pseudo.outcome1[folds1 != j])) * 
                    (1/sum(folds1 != j))
      eigen.j <- eigen(t(basis1[folds1 != j,]) %*% basis1[folds1 != j,]/sum(folds1 != j) + lam.seq[k] * diag(pweights), symmetric = TRUE)
      coef.j <- t(eigen.j$vectors) %*% diag(1/eigen.j$values) %*% eigen.j$vectors %*% contrast.j
      cate.hold.j <- basis1[folds1 == j,] %*% coef.j
      pred.error1[k, j] <- sum((pseudo.outcome1[folds1 == j] - cate.hold.j)^2)
    }
  }
  
  lam.cv0 <- lam.seq[which.min(apply(pred.error0, 2, mean))]
  lam.cv1 <- lam.seq[which.min(apply(pred.error1, 2, mean))]
  
  eigen.full1 <- eigen(t(basis1) %*% basis1/n1 + lam.cv1 * diag(pweights), symmetric = TRUE)
  coef1 <- t(eigen.full1$vectors) %*% diag(1/eigen.full1$values) %*% eigen.full1$vectors %*% one.step1 
  fitted0 <- basis0 %*% coef1
  eigen.full0 <- eigen(t(basis0) %*% basis0/n0 + lam.cv0 * diag(pweights), symmetric = TRUE)
  coef0 <- t(eigen.full0$vectors) %*% diag(1/eigen.full0$values) %*% eigen.full0$vectors %*% one.step0 
  fitted1 <- basis1 %*% coef0
  
  rule1 <- ifelse(fitted1 > delta, 1, 0)
  rule0 <- ifelse(fitted0 > delta, 1, 0)
  
  # Fourth, test using optimal rule in each data set
  qual0 <- OneStepQual(y = y0, a = a0, w = w0, 
                       prop.score = nuisance.est0$prop.score,
                       outcome.regression = nuisance.est0$outcome.regression,
                       outcome.regression.tx0 = nuisance.est0$outcome.regression.tx0, 
                       outcome.regression.tx1 = nuisance.est0$outcome.regression.tx1,
                       direc.basis = rule0, delta = delta)
  test.stat.pos0 <- n0 * qual0$one.step.pos^2/var(qual0$eif.pos) * (qual0$one.step.pos > 0)
  test.stat.neg0 <- n0 * qual0$one.step.neg^2/var(qual0$eif.neg) * (qual0$one.step.neg < 0)
  
  qual1 <- OneStepQual(y = y1, a = a1, w = w1, 
                       prop.score = nuisance.est1$prop.score,
                       outcome.regression = nuisance.est1$outcome.regression,
                       outcome.regression.tx0 = nuisance.est1$outcome.regression.tx0, 
                       outcome.regression.tx1 = nuisance.est1$outcome.regression.tx1,
                       direc.basis = rule1, delta = delta)
  test.stat.pos1 <- n1 * qual1$one.step.pos^2/var(qual1$eif.pos) * (qual1$one.step.pos > 0)
  test.stat.neg1 <- n1 * qual1$one.step.neg^2/var(qual1$eif.neg) * (qual1$one.step.neg < 0)
  
  if(all(rule0 == 1) | all(rule0 == 0)) {
    p.val0 <- 1
  } else {
    p.val0 <- .5 * max(c(pchisq(test.stat.pos0, 1, lower.tail = F),
                       pchisq(test.stat.neg0, 1, lower.tail = F)))
  }
  
  if(all(rule1 == 1) | all(rule1 == 0)) {
    p.val1 <- 1
  } else {
    p.val1 <- .5  *max(c(pchisq(test.stat.pos1, 1, lower.tail = F),
                       pchisq(test.stat.neg1, 1, lower.tail = F)))
  }
  p.val.adj <- min(c(p.val0 * 2, p.val1 * 2, 1))
  
  out <- list(p.val0 = p.val0,
              p.val1 = p.val1,
              p.val.adj = p.val.adj)
  return(out)
}

# Example
# set.seed(206)
# n <- 400
# w <- matrix(runif(3 * n, -1, 1), nrow = n, ncol = 3)
# a <- rbinom(n, 1, exp(.25 * (.5 * w[,1] + sin(pi/2 * w[,2])))/(1 + exp(.25 * (w[,1] + sin(pi * w[,2])))))
# h <- function(x) {x[1] + exp(.5 * (x[2] + x[3]))/(1 + exp(.5 * (x[2] + x[3])))}
# g.null <- function(x) {1}
# g.alt <- function(x) x[3]
# g <- function(x) {3 * exp(x)/(1 + exp(x)) * cos(x * pi/2)}
# g <- function(x) {2 * x^2 * sign(x)}
# b <- 1
# y.null <- apply(w, 1, h) - b * a * apply(w, 1, g.null) + rnorm(n, mean = 0, sd = 3)
# y.alt <- apply(w, 1, h) + b * a * g(w[,1]) + rnorm(n, mean = 0, sd = 2)
# 
# x <- w[,1]
# y <- y.alt
# # y <- y.null
# 
# out.thresh <- QuantTest(y = y.null, a = a, w = w, x = x,
#                  method = "threshold",
#                  no.thresholds = ceiling(sqrt(length(x))),
#                  thresholds = NULL,
#                  no.roots = 2, lambda = NULL,
#                  randomized = FALSE, no.bs = 1000)
# 
# out.intvl <- QuantTest(y = y.null, a = a, w = w, x = x,
#                  method = "interval",
#                  no.thresholds = ceiling(sqrt(length(x))),
#                  thresholds = NULL,
#                  no.roots = 2, lambda = NULL,
#                  randomized = FALSE, no.bs = 1000)
# 
# out.thresh.qual <- QualTest(y = y.alt, a = a, w = w, x = x,
#                             method = "threshold",
#                             no.thresholds = ceiling(sqrt(length(x))),
#                             thresholds = NULL,
#                             no.roots = 2, lambda = NULL,
#                             randomized = FALSE, no.bs = 1000)
# 
# out.intvl.qual <- QualTest(y = y.alt, a = a, w = w, x = x,
#                             method = "interval",
#                             no.thresholds = ceiling(sqrt(length(x))),
#                             thresholds = NULL,
#                             no.roots = 3, lambda = NULL,
#                             randomized = FALSE, no.bs = 1000)
