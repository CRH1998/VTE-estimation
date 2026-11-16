################################################################################
#                           Surival multiplier bootstrap                       #
################################################################################
#                                                                              #
# This script is used to test bootstrap the survival U-statistic               #
#                                                                              #
################################################################################
#                                                                              #
#                               Dependencies                                   #
#                                                                              #
# This script does not depend on any other scripts                             #
################################################################################



# Literal integral implementation for independent numerical verification.
multiplier_test <- function(object, n_boot = 999L, alpha = 0.05,
                            distribution = c("normal", "rademacher"),
                            seed = NULL, multipliers = NULL, return_boot = TRUE) {
  distribution <- match.arg(distribution)
  stopifnot(alpha > 0, alpha < 1)
  
  #browser()
  
  if (is.null(multipliers)) {
    stopifnot(n_boot >= 1L, n_boot == as.integer(n_boot))
    multipliers <- 
      if (distribution == "normal") matrix(stats::rnorm(object$n * n_boot), object$n)
      else matrix(sample(c(-1, 1), object$n * n_boot, replace = TRUE), object$n)
  } else {
    stopifnot(is.matrix(multipliers), nrow(multipliers) == object$n,
              ncol(multipliers) > 0, all(is.finite(multipliers)))
    n_boot <- ncol(multipliers)
    distribution <- "user-supplied"
  }
  
  
  
  # One weight per OBSERVATION, shared across ALL blocks. Do not independently
  # weight pairs, re-center the weights, refit nuisances, or subtract T_n here.
  boot <- quadratic_statistic(object, multipliers)
  
  # Plus-one Monte Carlo p-value; this does not imply exact finite-sample size.
  p_value <- (1 + sum(boot >= object$statistic)) / (n_boot + 1)
  
  # Order-statistic critical value aligned with that p-value (strict rejection).
  k <- ceiling((1 - alpha) * (n_boot + 1))
  critical <- if (k > n_boot) Inf else sort(boot, partial = k)[k]
  list(statistic = object$statistic, estimate = object$estimate,
       p_value = p_value, reject = p_value <= alpha,
       critical_value = critical, alpha = alpha, n_boot = n_boot,
       distribution = distribution, clipped = any(object$diagnostics$clipped),
       bootstrap = if (return_boot) boot else NULL)
}