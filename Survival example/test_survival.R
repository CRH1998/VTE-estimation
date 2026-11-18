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
source("survival_plots.R")

################################################################################
#                              Load required libraries                         #
################################################################################
# For comparison with log-rank test 
library(survival)
library(survminer)
library(Rcpp)
library(ggplot2)
library(patchwork)

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
# 1. Simulation study with proportional hazards
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
  
  # Truncate the survival times at the upper limit
  time_tau  <- pmin(surv_data$time, tau)
  event_tau <- as.numeric(surv_data$time <= tau)
  
  surv_obj <- survival::Surv(time = time_tau, event = event_tau)
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
tau <- 2            # Upper time limit
n_values <- c(400, 800, 1200, 1600, 2000, 3000) # Sample sizes to test
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
  summarise(
    LogRankLevel = if (first(beta) == 0) {
      mean(log_rank_p_value <= alpha)
    } else NA_real_,
    TwoStepLevel = if (first(beta) == 0) {
      mean(mult_test_p_value <= alpha)
    } else NA_real_,
    LogRankPower = if (first(beta) > 0) {
      mean(log_rank_p_value <= alpha)
    } else NA_real_,
    TwoStepPower = if (first(beta) > 0) {
      mean(mult_test_p_value <= alpha)
    } else NA_real_,
    RMSE = sqrt(mean((estimate - oracle_val)^2)),
    oracle_val = mean(oracle_val),
    est = mean(estimate),
    var = var(estimate),
    .groups = "drop"
  ) %>%
  arrange(n, beta)


p_levels <- ggplot_empirical_levels(summary_results, alpha = 0.05)
p_power <- ggplot_empirical_power(summary_results)
p_rmse <- ggplot_empirical_RMSE(summary_results)
p_combined <- combine_article_plots(p_levels, p_power, p_rmse)

# print(p_levels)
# print(p_power)
# print(p_rmse)
print(p_combined)


ggsave(
  "survival_simulation_combined.pdf",
  plot = p_combined,
  width = 200, height = 100, units = "mm",
  device = grDevices::pdf,
  useDingbats = FALSE, bg = "white"
)





################################################################################
# 2. Simulation study with non-proportional hazards
################################################################################

# Repeat the above simulation study but with simulate_nonproportional_survival instead of simulate_survival, and adjust the run_one_standard function accordingly.

#-------------------------------------------------------------------------------
# 2.1 Define function to run one simulation with non-proportional hazards
#-------------------------------------------------------------------------------
run_one_nonproportional <- function(n, beta, delta, p, rate, tau, K = 2L, n_boot = 1999L, alpha = 0.05) {
  
  # Simulate non-proportional hazards survival data
  surv_data <- simulate_nonproportional_survival(n = n, 
                                                 beta = beta, 
                                                 delta = delta,
                                                 p = p, 
                                                 rate = rate)
  
  
  # Estimate the cross-fit survival function and calculate the multiplier test
  cf_surv <- crossfit_survival(data = surv_data, tau = tau, K = K, fit_nuisance = fit_nuisance)
  mult_test <- multiplier_test(object = cf_surv, 
                               n_boot = n_boot, 
                               alpha = alpha,
                               distribution = "normal")
  mult_test_p_value <- mult_test$p_value
  
  # Estimate survival function and log-rank test with no censoring and no covariates
  time_tau  <- pmin(surv_data$time, tau)
  event_tau <- as.numeric(surv_data$time <= tau)
  
  surv_obj <- survival::Surv(time = time_tau, event = event_tau)
  surv_fit <- survival::survfit(surv_obj ~ surv_data$A)
  log_rank_test <- survival::survdiff(surv_obj ~ surv_data$A)
  log_rank_p_value <- log_rank_test$pvalue
  
  
  list(n = n, 
       beta = beta,
       delta = delta,
       mult_test_p_value = mult_test$p_value,
       log_rank_p_value = log_rank_p_value)
}

#-------------------------------------------------------------------------------
# 2.2 Simulation study settings for non-proportional hazards
#-------------------------------------------------------------------------------

seed <- 17112026    # Random seed for reproducibility
p <- 0.5            # Probability of treatment assignment
rate <- 1           # Rate parameter
tau <- 2            # Upper time limit
n_values <- c(400, 800, 1200, 1600, 2400)                        # Sample sizes to test
beta_values <- c(0, 0.4054651)                                   # Baseline log hazard ratios corresponding to HR = 1, 1.5
delta_values <- c(0, -0.4054651, -0.6931472, -1.098612)          # Time-dependent effect values
alpha_nonprop <- 0.05

n_rep_nonprop <- 250L                                            # Number of replicates per setting

settings_nonprop <- expand.grid(
  n = n_values,
  beta = beta_values,
  delta = delta_values,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

# Each setting gets n_rep consecutive runs
setting_index_nonprop <- rep(seq_len(nrow(settings_nonprop)), each = n_rep_nonprop)
                             
#------------------------------------------------------------------------------
# 2.3 Run the simulation study for non-proportional hazards
#------------------------------------------------------------------------------
start.time <- Sys.time()
old_plan_nonprop <- future::plan(future::multisession, workers = 4L)

runs_nonprop <- tryCatch(
  future.apply::future_lapply(
    seq_along(setting_index_nonprop),
    function(sim_id) {
      setting <- settings_nonprop[setting_index_nonprop[sim_id], , drop = FALSE]
      run_one_nonproportional(
        n = setting$n,                           # Number of individuals
        beta = setting$beta,                     # Log hazard ratio
        delta = setting$delta,                   # Time-dependent effect
        p = p,                                   # Probability of treatment assignment
        rate = rate,                             # Rate parameter for exponential distribution
        tau = tau                                # Upper time limit
      )
    },
    future.chunk.size = 1L,
    future.seed = seed
  ),
  finally = future::plan(old_plan_nonprop)
)
end.time <- Sys.time()
end.time - start.time

#-------------------------------------------------------------------------------
# 2.4 Process the results
#-------------------------------------------------------------------------------  
# One list element per setting, with replicates in columns
results_nonprop <- lapply(seq_len(nrow(settings_nonprop)), function(j) {
  simplify2array(runs_nonprop[setting_index_nonprop == j], higher = FALSE)
})

# Label each element with its parameter combination
names(results_nonprop) <- with(
  settings_nonprop,
  paste0(
    "n=", n,
    "_beta=", beta,
    "_delta=", delta
  )
)


# Process results
results_df_nonprop <- lapply(results_nonprop, function(sim_result) {
  matrix(unlist(sim_result), ncol = 5, byrow = T)
})
results_df_nonprop <- data.frame(do.call(rbind, results_df_nonprop))
colnames(results_df_nonprop) <- c("n", 
                          "beta", 
                          "delta",
                          "mult_test_p_value", 
                          "log_rank_p_value")



# Summarize results
summary_results_nonprop <- results_df_nonprop %>%
  group_by(n, beta, delta) %>%
  summarise(
    LogRankReject = mean(log_rank_p_value <= alpha_nonprop),
    TwoStepReject = mean(mult_test_p_value <= alpha_nonprop),
    n_rep = dplyr::n(),
    .groups = "drop"
  ) %>%
  mutate(
    is_null = beta == 0 & delta == 0,
    LogRankLevel = if_else(is_null, LogRankReject, NA_real_),
    TwoStepLevel = if_else(is_null, TwoStepReject, NA_real_),
    LogRankPower = if_else(!is_null, LogRankReject, NA_real_),
    TwoStepPower = if_else(!is_null, TwoStepReject, NA_real_)
  ) %>%
  arrange(n, beta, delta)


p_levels_nonprop <- ggplot_empirical_levels_nonprop(
  summary_results_nonprop, alpha = alpha_nonprop
)
p_power_nonprop <- ggplot_empirical_power_nonprop(summary_results_nonprop)
p_combined_nonprop <- combine_article_plots_nonprop(
  p_levels_nonprop, p_power_nonprop, alpha = alpha_nonprop
)
print(p_combined_nonprop)


ggsave(
  "survival_non_prop_simulation_combined.pdf",
  plot = p_combined_nonprop,
  width = 200, height = 100, units = "mm",
  device = grDevices::pdf,
  useDingbats = FALSE, bg = "white"
)


