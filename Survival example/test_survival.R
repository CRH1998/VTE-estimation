################################################################################
#                           Test survival functions                            #
################################################################################
#                                                                              #
# This script is used to test the survival estimation and simulation functions #
#                                                                              #
################################################################################
#                                                                              #
#                               Dependencies                                   #
#                                                                              #
# survival_simluation.R                                                        #
# survival_estimation.R                                                        #
# u_statistic_crossfit.R                                                       #
# survival_multiplier_bootstrap.R                                              #
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("survival_simulation.R")
source("survival_estimation.R")
source("u_statistic_crossfit.R")
source("survival_multiplier_bootstrap.R")

################################################################################
#                              Load required libraries                         #
################################################################################
# For comparison with log-rank test 
library(survival)
library(survminer)
library(Rcpp)
library(ggplot2)

################################################################################


################################################################################
#-------------------------------------------------------------------------------
# 0. Testing functions
#-------------------------------------------------------------------------------
# Simulation specs
n <- 1000           # Number of individuals
beta <- 0.01        # Log hazard ratio
p <- 0.5            # Probability of treatment assignment
rate <- 1           # Baseline hazard rate
tau <- 1            # Upper time limit

#-------------------------------------------------------------------------------
# 0.1 Simulate survival data
#-------------------------------------------------------------------------------
surv_data <- simulate_survival(n = n, beta = beta, p = p, rate = rate)

#-------------------------------------------------------------------------------
# 0.2 Estimate functions
#-------------------------------------------------------------------------------
# Population estimand for reference
oracle_val <- population_estimand(beta = beta, rate = rate, p = p, tau = tau)

# Nuisance parameters
nuisance_fit <- fit_nuisance(train = surv_data, tau = tau)
check_survival_vals(nuisance_fit, threshold = 0.05)

cf_surv <- crossfit_survival(data = surv_data, tau = tau, K = 2L, fit_nuisance = fit_nuisance)
mult_test <- multiplier_test(object = cf_surv, n_boot = 1000L, alpha = 0.05,
                             distribution = "normal")
oracle_val
cf_surv$estimate
mult_test$p_value
################################################################################



################################################################################
# 1. Simulation study
################################################################################

#-------------------------------------------------------------------------------
# 1.1 Define function to run one simulation
#-------------------------------------------------------------------------------
run_one_standard <- function(n, beta, p, rate, tau, K = 2L, n_boot = 1999L, alpha = 0.05) {
  
  # Simulate survival data
  surv_data <- simulate_survival(n = n, beta = beta, p = p, rate = rate)
  
  # Calculate the oracle value for reference
  oracle_val <- population_estimand(beta = beta, rate = rate, p = p, tau = tau)
  
  # Estimate the cross-fit survival function and calculate the multiplier test
  cf_surv <- crossfit_survival(data = surv_data, tau = tau, K = K, fit_nuisance = fit_nuisance)
  mult_test <- multiplier_test(object = cf_surv, 
                               n_boot = n_boot, 
                               alpha = alpha,
                               distribution = "normal")
  mult_test_p_value <- mult_test$p_value
  
  # Estimate survival function and log-rank test with no censoring and no covariates
  surv_obj <- survival::Surv(time = surv_data$time, event = rep(1, n))
  surv_fit <- survival::survfit(surv_obj ~ surv_data$A)
  log_rank_test <- survival::survdiff(surv_obj ~ surv_data$A)
  log_rank_p_value <- log_rank_test$pvalue
  
  
  list(n = n, 
       beta = beta,
       oracle_val = oracle_val,
       estimate = cf_surv$estimate,
       RMSE = sqrt(mean((cf_surv$estimate - oracle_val)^2)),
       mult_test_p_value = mult_test$p_value,
       log_rank_p_value = log_rank_p_value)
}



#-------------------------------------------------------------------------------
# 1.2 Simulation study settings
#-------------------------------------------------------------------------------
seed <- 16112026    # Random seed for reproducibility
p <- 0.5            # Probability of treatment assignment
rate <- 1           # Baseline hazard rate
tau <- 1            # Upper time limit
n_values <- c(400, 800, 1200, 1600, 2000) # Sample sizes to test
beta_values <- c(0, 0.09531018, 0.1823216, 0.2623643, 0.4054651) # Corresponding to HR = 1, 1.1, 1.2, 1.3, 1.5

n_rep <- 250L        # Number of replicates per setting

settings <- expand.grid(
  n = n_values,
  beta = beta_values,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

# Each setting gets n_rep consecutive runs
setting_index <- rep(seq_len(nrow(settings)), each = n_rep)

#-------------------------------------------------------------------------------
# 1.3 Run the simulation study
#-------------------------------------------------------------------------------
start.time <- Sys.time()
old_plan <- future::plan(future::multisession, workers = 4L)

runs <- tryCatch(
  future.apply::future_lapply(
    seq_along(setting_index),
    function(sim_id) {
      setting <- settings[setting_index[sim_id], , drop = FALSE]
      run_one_standard(
        n = setting$n,           # Number of individuals
        beta = setting$beta,     # Log hazard ratio
        p = p,                   # Probability of treatment assignment
        rate = rate,             # Baseline hazard rate
        tau = tau                # Upper time limit
      )
    },
    future.chunk.size = 1L,
    future.seed = seed
  ),
  finally = future::plan(old_plan)
)
end.time <- Sys.time()
end.time - start.time

#-------------------------------------------------------------------------------
# 1.4 Process the results
#-------------------------------------------------------------------------------  
# One list element per setting, with replicates in columns
results <- lapply(seq_len(nrow(settings)), function(j) {
  simplify2array(runs[setting_index == j], higher = FALSE)
})

# Label each element with its parameter combination
names(results) <- with(
  settings,
  paste0(
    "n=", n,
    "_beta=", beta
  )
)


# Process results
results_df <- lapply(results, function(sim_result) {
  matrix(unlist(sim_result), ncol = 7, byrow = T)
})
results_df <- data.frame(do.call(rbind, results_df))
colnames(results_df) <- c("n", 
                          "beta", 
                          "oracle_val",
                          "estimate", 
                          "RMSE", 
                          "mult_test_p_value", 
                          "log_rank_p_value")



# Summarize results
summary_results <- results_df %>% 
  group_by(n, beta) %>%                               
  summarise(LogRankLevel = mean(log_rank_p_value > 0.05 & beta == 0),
            TwoStepLevel = mean(mult_test_p_value > 0.05 & beta == 0),
            LogRankPower = mean(log_rank_p_value < 0.05 & beta > 0),
            TwoStepPower = mean(mult_test_p_value < 0.05 & beta > 0),
            oracle_val = mean(oracle_val),
            est = mean(estimate),
            var = var(estimate),
            RMSE = mean(RMSE),
            .groups = "drop") %>%
  arrange(n, beta)




# Plot empirical levels with n on the x-axis and empirical level on the y-axis, with separate lines for each test (log-rank and multiplier) represented by different line types and different values of beta represented by different colors. Use ggplot2 for visualization.
ggplot_empirical_levels <- function(summary_results) {
  summary_results_long <- summary_results %>%
    filter(beta == 0) %>%
    pivot_longer(cols = c(LogRankLevel, TwoStepLevel), 
                 names_to = "test", 
                 values_to = "empirical_level")
  
  ggplot(summary_results_long, aes(x = n, y = empirical_level, color = factor(beta), linetype = test)) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
    labs(title = "Empirical Levels of Tests by Sample Size and Beta",
         x = "Sample Size (n)",
         y = "Empirical Level",
         color = "Beta",
         linetype = "Test") +
    theme_minimal()
}

ggplot_empirical_levels(summary_results)



# Plot empirical power with n on the x-axis and empirical power on the y-axis, with separate lines for each test (log-rank and multiplier) represented by different line types and different values of beta represented by different colors. Use ggplot2 for visualization.
ggplot_empirical_power <- function(summary_results) {
  summary_results_long <- summary_results %>%
    filter(beta > 0) %>%
    pivot_longer(cols = c(LogRankPower, TwoStepPower), 
                 names_to = "test", 
                 values_to = "empirical_power")
  
  ggplot(summary_results_long, aes(x = n, y = empirical_power, color = factor(beta), linetype = test)) +
    geom_line() +
    geom_point() +
    labs(title = "Empirical Power of Tests by Sample Size and Beta",
         x = "Sample Size (n)",
         y = "Empirical Power",
         color = "Beta",
         linetype = "Test") +
    theme_minimal()
}

ggplot_empirical_power(summary_results)



# Plot empirical RMSE with n on the x-axis and empirical RMSE on the y-axis, with different values of beta represented by different colors. Use ggplot2 for visualization.
ggplot_empirical_RMSE <- function(summary_results) {
  ggplot(summary_results, aes(x = n, y = RMSE, color = factor(beta))) +
    geom_line() +
    geom_point() +
    labs(title = "Empirical RMSE of reduced two-step estimator by Sample Size and Beta",
         x = "Sample Size (n)",
         y = "Empirical RMSE",
         color = "Beta") +
    theme_minimal()
}
ggplot_empirical_RMSE(summary_results)
