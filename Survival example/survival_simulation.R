################################################################################
#                           Simulation script                                  #
################################################################################
#                                                                              #
#         This script contains functions to simulate survival data             #
#                                                                              #
################################################################################
#                                                                              #
#                               Dependencies                                   #
#                                                                              #
#             No scripts are need for this script                              #
################################################################################

library(simsurv)

#-------------------------------------------------------------------------------
# 1. Simulate proportional hazards survival data
#-------------------------------------------------------------------------------

#' Title
#'
#' @param n number of observations to simulate
#' @param beta log hazard ratio. beta = 0 gives independence exactly.
#' @param p probability of treatment assignment
#' @param rate baseline hazard rate
#' @param seed random seed for reproducibility
#'
#' @returns A data frame with columns A (treatment assignment) and time (survival time)
#' @export
#'
#' @examples

simulate_survival <- function(n, beta = 0, p = 0.5, rate = 1, seed = NULL) {
  stopifnot(length(n) == 1L, n >= 2, n == as.integer(n),
            length(beta) == 1L, is.finite(beta), p > 0, p < 1, rate > 0)
  
  set.seed(seed)
  
  A <- stats::rbinom(n, 1, p)
  data.frame(A = A, time = stats::rexp(n, rate = rate * exp(beta * A)))

  }

#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
# 2. Population survival function for proportional hazards model
#-------------------------------------------------------------------------------

#' Title
#'
#' @param u time point at which to evaluate the population survival function
#' @param beta log hazard ratio. beta = 0 gives independence exactly.
#' @param p probability of treatment assignment
#' @param rate baseline hazard rate
#'
#' @returns The population survival probability at time u
#' @export
#'
#' @examples
#' 
population_survival <- function(u, beta = 0, p = 0.5, rate = 1) {
  (1 - p) * exp(-rate * u) + p * exp(-rate * exp(beta) * u)
}

#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
# 3. Simulate non-proportional hazards survival data
#-------------------------------------------------------------------------------

simulate_nonproportional_survival <- function(n, 
                                              beta = 0,
                                              time_effect = 0.5,
                                              p = 0.5, 
                                              scale = 1,
                                              shape = 1.5,
                                              seed = NULL) {
  #browser()
  set.seed(seed)
  
  
  # 1. Define a dataset with a binary treatment covariate
  covs <- data.frame(
    id = seq_len(n),
    A = rbinom(n, 1, p)
  )
  
  # 2. Simulate survival times with a time-dependent effect of 'trt'
  # Suppose the log hazard ratio starts negative (protective effect) 
  # but changes linearly with time (non-proportional hazard).
  sim_data <- simsurv(
    dist = "weibull",
    lambdas = scale,
    gammas = shape,
    x = covs,
    betas = c(A = beta),                         # base log HR
    tde = c(A = time_effect),       # time-dependent interaction
    maxt = 5
  )
  
  # Merge survival times back with covariates
  data <- merge(covs, sim_data, by = "id")
  data$id <- NULL      # Remove the id column for clarity
  data$status <- NULL  # Remove the status column since no censoring is simulated here
  names(data)[names(data) == "eventtime"] <- "time"
  return(data)
}
#-------------------------------------------------------------------------------


beta_values <- c(0, 0.1823216, 0.2623643, 0.4054651)   # Baseline log hazard ratios corresponding to HR = 1, 1.2, 1.3, 1.5
time_effect_values <- c(0, -0.1823216, -0.2623643, -0.4054651)               # Time-dependent effect values to test

test_data <- simulate_nonproportional_survival(n = 400, 
                                               beta = 0.4054651, 
                                               time_effect = -0.1823216, 
                                               p = 0.5, 
                                               scale = 1, 
                                               shape = 1, 
                                               seed = NULL)

km_fit <- survfit(Surv(time, rep(1, nrow(test_data))) ~ A, data = test_data)

# Option 1: Publication-ready plot with risk table using survminer
ggsurvplot(
  km_fit,
  data = test_data,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  legend.labs = c("Control (A = 0)", "Treatment (A = 1)"),
  xlab = "Time",
  ylab = "Survival Probability"
)



