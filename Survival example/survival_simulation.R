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

#-------------------------------------------------------------------------------
# 1. Simulate survival data
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
# 2. Population survival function
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




