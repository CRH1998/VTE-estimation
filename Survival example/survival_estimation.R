################################################################################
#                          Survival estimation functions                       #
################################################################################
#                                                                              #
#         This script contains functions to estimate the survival type         #
#         estimand and its variance using cross-fitting and U-statistics.      #
#                                                                              #
################################################################################
#                                                                              #
#                               Dependencies                                   #
#                                                                              #
# No scripts are needed for this script                                        #
#                                                                              #
#                                                                              #
################################################################################


#-------------------------------------------------------------------------------
# 1. Defining population estimand to compare with cross-fitted estimate
#-------------------------------------------------------------------------------

population_estimand <- function(beta, tau = 1, p = 0.5, rate = 1) {
  if (beta == 0) return(0)
  stats::integrate(function(u) {
    s0 <- exp(-rate * u)
    s1 <- exp(-rate * exp(beta) * u)
    s <- (1 - p) * s0 + p * s1
    contrast <- p - p * s1 / s
    density <- (1 - p) * rate * s0 + p * rate * exp(beta) * s1
    contrast^2 * density
  }, 0, tau, rel.tol = 1e-10)$value
}




#-------------------------------------------------------------------------------
# 2. Define function to estimate nuisance parameters
#-------------------------------------------------------------------------------

fit_nuisance <- function(train, tau) {
  Pn <- stats::ecdf(train$time)   # The empirical measure for the time random variable
  u <- sort(unique(train$time))   # Sort the time points for integration
  Pn_u <- Pn(u)                   # Estimated CDF at time u using empirical measure
  dPn_u <- diff(c(0, Pn_u))       # Estimated empirical differential
  
  keep <- u <= tau                # Values of time less than tau to keep
  
  list(piP = mean(train$A), grid = u[keep], dPn_u = dPn_u[keep],
       Sn_u = 1 - Pn_u[keep], S_tau = 1 - Pn(tau), method = "ECDF")
}


check_survival_vals <- function(fit, threshold = 0.01){
  S_tau <- fit$S_tau
  if (S_tau < threshold) warning("Survival value at tau is", S_tau, "which is below the threshold of", threshold, ".")
  
  return(fit)
}



#-------------------------------------------------------------------------------
# 3. Define function to estimate estimand
#-------------------------------------------------------------------------------

# Primitive representation of the exact empirical Stieltjes integral:
# h_ij = 2(A_i-piP)(A_j-piP){c0-c1(t_i)-c1(t_j)+c2(min(t_i,t_j))}.
# c_r(t) = sum_{u < t, u <= tau} dPn_u(u) / S(u)^r.
# left.open=TRUE implements the strict inequality in 1(time > u),
# including ties and atoms at the fixed truncation endpoint tau.
.empirical_primitives <- function(time, fit) {
  j <- findInterval(time, fit$grid, left.open = TRUE)
  c1 <- c(0, cumsum(fit$dPn_u / fit$Sn_u))[j + 1L]
  c2 <- c(0, cumsum(fit$dPn_u / fit$Sn_u^2))[j + 1L]
  list(c0 = sum(fit$dPn_u), c1 = c1, c2 = c2)
}




.column_cumsum <- function(x) {
  if (getOption("survival_u_test.use_matrixStats", TRUE) &&
      requireNamespace("matrixStats", quietly = TRUE))
    return(matrixStats::colCumsums(x))
  out <- x
  for (j in seq_len(ncol(x))) out[, j] <- cumsum(x[, j])
  out
}


# U-statistic sum over pairs in a triangle block, with the primitive representation.
.triangle_sums <- function(x, c1, d) {
  if (nrow(x) < 2L) return(rep(0, ncol(x)))
  xd <- x * d
  prefix_x <- .column_cumsum(x) - x
  prefix_d <- .column_cumsum(xd) - xd
  2 * colSums(x * (prefix_d - c1 * prefix_x))
}


quadratic_statistic <- function(object, multipliers = NULL) {
  
  #-----------------------------------------------------------------------------
  # 3.0 Validate inputs
  #-----------------------------------------------------------------------------
  if (is.null(multipliers)) multipliers <- matrix(1, object$n, 1)
  if (is.vector(multipliers)) multipliers <- matrix(multipliers, ncol = 1)
  stopifnot(is.matrix(multipliers), nrow(multipliers) == object$n,
            ncol(multipliers) >= 1L, all(is.finite(multipliers)))
  
  #-----------------------------------------------------------------------------
  # 3.1 Compute quadratic statistic for each multiplier
  #-----------------------------------------------------------------------------
  ans <- rep(0, ncol(multipliers))
  for (b in object$blocks) {
    
    # Inner sum over pairs in the block, weighted by the multipliers and the centered treatment assignment.
    x <- multipliers[b$idx, , drop = FALSE] * b$am
    val <- .triangle_sums(x, b$c1, b$d)
    if (b$type == "rectangle") {
      # All pairs in the union minus within-half pairs = this rectangle.
      for (label in 1:2) {
        keep <- b$label == label
        val <- val - .triangle_sums(x[keep, , drop = FALSE],
                                    b$c1[keep], b$d[keep])
      }
    }
    
    # Summing over blocks (k = 1,...,K) which is l = 1,..., L in the latex document.
    ans <- ans + val
  }
  ans / (object$n - 1)
}






#-------------------------------------------------------------------------------
# 4. Define cross-fit estimator
#-------------------------------------------------------------------------------

crossfit_survival <- function(data, 
                              tau = 1, 
                              K = 2L,
                              fit_nuisance = fit_nuisance_ecdf,
                              s_floor = 0.05, split_seed = NULL) {

  #browser()
  n <- nrow(data)
  folds <- make_pair_folds(n, K, seed = split_seed)
  stopifnot(inherits(folds, "pair_folds"), folds$n == n)
  
  fitted <- vector("list", length(folds$blocks))
  diagnostics <- vector("list", length(fitted))
  
  for (l in seq_along(fitted)) {
    
    #---------------------------------------------------------------------------
    # 4.0 Divide data into train and evaluation sets for the current fold
    #---------------------------------------------------------------------------
    b <- folds$blocks[[l]]
    btrain <- b$train
    beval <- b$eval
    
    #---------------------------------------------------------------------------
    # 4.1 Fit nuisance parameters on training data
    #---------------------------------------------------------------------------
    
    fit <- check_survival_vals(fit_nuisance(data[btrain, , drop = FALSE], tau),
                               threshold = s_floor)
    
    
    
    #---------------------------------------------------------------------------
    # 4.2 Calculate empirical primitives on evaluation data
    #---------------------------------------------------------------------------
    idx <- beval[order(data$time[beval], beval)]
    pr <- .empirical_primitives(data$time[idx], fit)
    
    
    #---------------------------------------------------------------------------
    # 4.3 Store results for the current fold
    #---------------------------------------------------------------------------
    fitted[[l]] <- list(type = b$type, idx = idx,
                        label = ifelse(idx %in% b$left, 1L, 2L),
                        am = data$A[idx] - fit$piP, 
                        c0 = pr$c0, 
                        c1 = pr$c1, 
                        c2 = pr$c2,
                        d = pr$c0 - pr$c1 + pr$c2, 
                        nuisance = fit)
    
    diagnostics[[l]] <- data.frame(block = l, 
                                   type = b$type,
                                   training_n = length(btrain), 
                                   evaluation_n = length(beval),
                                   n_pairs = b$n_pairs, 
                                   piP = fit$piP, 
                                   raw_S_tau = fit$S_tau)
  }
  
  #---------------------------------------------------------------------------
  # 4.4 Create the cross-fitted survival object
  #---------------------------------------------------------------------------
  object <- structure(list(n = n, tau = tau, K = folds$K, data = data,
                           folds = folds, blocks = fitted, s_floor = s_floor,
                           diagnostics = do.call(rbind, diagnostics)), class = "survival_cf")
  object$statistic <- quadratic_statistic(object)[1L]  # n * estimate
  object$estimate <- object$statistic / n
  
  return(object)
}


