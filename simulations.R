############################################################
#                                                          #
#   The following script includes simulation functions     #
#   for the test of heterogenous effects simulation study  #
#                                                          #
############################################################


############################################################
library(MASS)
library(lava)

############################################################


############################################################
#                 Sánchez-Becerra simulations              #
############################################################


#' Title
#'
#' @param n sample size
#' @param p number of covariates in X0 and X1
#' @param rho correlation between X0j and X1j
#' @param mu mean of X0 and X1
#' @param Dp Treatment assignment probability
#' @param tau Average effect size
#' @param c Hyperparameter
#' @param sigma_d Hyperparameter
#' @param tilde_sigma_d Hyperparameter
#' @param lambda Exponential decay of coefficients
#' @param V_mu Baseline variance (variance of X0)
#' @param V_tau VCATE
#'
#' @returns Simulated data set
#' @export
#'
#' @examples
sim_SB <- function(n, J = 5,
                   rho = 0.5, mu = 0,
                   Ap = 0.5,
                   tau = 0.15, c = 1, sigma_d = sqrt(0.7), tilde_sigma_d = 0.21, lambda = 0.7,
                   V_mu = 0.3, V_tau){
  
  #browser()
  # Calculating derived hyperparameters
  ell <- sqrt((1 - lambda)/(1 - lambda^J) * lambda^(1:J-1))     # Geometric sequence for decaying effect
  
  kappa0 <- ell * sqrt(max(sigma_d^2 - tilde_sigma_d^2, 0))
  kappa1 <- ell * sqrt(max(sigma_d^2 - tilde_sigma_d^2, 0))
  beta0 <- ell * sqrt(V_mu)
  beta1 <- ell * sqrt(V_tau)
  
  # Simulating covariates
  Ij <- diag(J)
  Sigma_x <- rbind(cbind(Ij, rho * Ij), 
                   cbind(rho * Ij, Ij))
  mu_x <- rep(mu, 2*J)
  
  X <- MASS::mvrnorm(n = n, mu = rep(0, 2 * J), Sigma = Sigma_x)
  X0 <- X[, 1:J, drop = FALSE]
  X1 <- X[, (J + 1):(2 * J), drop = FALSE]
  
  
  # Simulating treatment
  A <- rbinom(n, 1, Ap)
  
  
  # Simulating errors
  U <- mvrnorm(n = n, mu = rep(0,2), Sigma = diag(2))
  U0 <- U[,1]
  U1 <- U[,2]
  
  
  # Simulating outcome
  Y0 <- c + as.numeric(X0 %*% beta0) +
    U0 * sqrt(tilde_sigma_d^2 + rowSums((X0 %*% diag(kappa0))^2))
  
  Y1 <- (c + tau) + as.numeric(X0 %*% beta0) + as.numeric(X1 %*% beta1) +
    U1 * sqrt(tilde_sigma_d^2 + rowSums((X1 %*% diag(kappa1))^2))
  
  Y <- A * Y1 + (1 - A) * Y0
  # Return output
  return(data.frame(Y, A, X))
}





############################################################
#                 Dukes et al. simulations                 #
############################################################


sim_Dukes <- function(n, heterogeneity = c("None", "QuantHMC", "QuantHNMC"), randomized = FALSE,
                      noise_mean = 0, noise_sd = 3, heterogeneity_degree = 15, pi_01x = 0.5) {
  
  heterogeneity <- match.arg(heterogeneity)
  
  X1 <- runif(n = n, min = -1, max = 1)
  X2 <- runif(n = n, min = -1, max = 1)
  X3 <- runif(n = n, min = -1, max = 1)
  
  if (randomized) {
    pi_01x <- pi_01x
  } else {
    pi_01x <- expit(1/8 * X1 + 1/4*sin(pi*X2))
  }
  A <- rbinom(n = n, size = 1, prob = pi_01x)
  
  hx <- X1 + expit(1/2*(X2 + X3))
  
  if(heterogeneity == "None"){ # No heterogeneity
    gamma <- 3/4
  } else if (heterogeneity == "QuantHMC"){ # Quantitative heterogeneity; monotone CATE
    gamma <- heterogeneity_degree*(X3 - 0.5) * (X3 > 0.5)
  } else if (heterogeneity == "QuantHNMC"){ # Quantitative heterogeneity; non-monotone CATE
    gamma <- 3*(1 - X3^2)
  }
  
  Y <- hx + A * gamma + rnorm(n = n, mean = noise_mean, sd = noise_sd)
  
  return(data.frame(Y = Y, A = A, X1 = X1, X2 = X2, X3 = X3))
}














  
