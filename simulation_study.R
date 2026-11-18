#############################
#     Simulation study      #
#############################


############################################################
#                         Packages                         #
############################################################
library(future.apply)
library(future)


############################################################
#                         Dependencies                     #
############################################################
# Testing functions and code
source("/Users/christianhejstvig-larsen/Desktop/Studie/Master/Cambridge/Essay/R/Part-III-Essay-Repo/Estimation functions.R")
source("/Users/christianhejstvig-larsen/Desktop/Studie/Master/Cambridge/Essay/R/Part-III-Essay-Repo/Simulations.R")
source("/Users/christianhejstvig-larsen/Desktop/Studie/Master/Cambridge/Essay/R/code-JASA/np-qualint-JASA.R")
source("nuisance_func.R")



run_one <- function(n, sim_type = c("SB", "Dukes", "Hines", "PM"), randomized = TRUE, pi_hat = 0.5, seed = NA,
                    J, tau = 0.15,                                                                    # SB simulation settings
                    heterogeneity = c("None", "QuantHMC", "QuantHNMC"), heterogeneity_degree,         # Dukes simulation settings
                    sigma = 0.5, alternative_type = "smooth", sim_scenario = "regular",   # PM simulation settings
                    analysis_type = c("SB", "Dukes", "CRH"),                                          # Analysis type
                    family = gaussian(), K = 2, hal = F,                                              # Analysis settings
                    SL.lib = list("SL.glmnet", "SL.xgboost", "SL.ranger", "SL.mean"), 
                    i = NA_integer_) {
  
  
  sim_type <- match.arg(sim_type)
  #analysis_type <- match.arg(analysis_type)
  
  if (!is.na(seed)) {
    set.seed(seed)
  }
  
  
  # Choose simulation type
  if (sim_type == "SB") {
    
    V_tau <- heterogeneity_degree
    data <- sim_SB(n = n, J = J, V_tau = V_tau, tau = tau, Ap = pi_hat)
    print(paste("Running simulation", i, "with n =", n, "and V_tau =", V_tau))
    
  } else if (sim_type == "Dukes") {
    
    data <- sim_Dukes(n = n, heterogeneity = heterogeneity, 
                      heterogeneity_degree = heterogeneity_degree, 
                      randomized = randomized, 
                      pi_01x = pi_hat)
    
    V_tau <- 13/768 * heterogeneity_degree^2
    print(paste("Running simulation", i, "with n =", n, "and V_tau =", V_tau))
    
  } else if (sim_type == "Hines") {
    
    data <- sim_Hines(n = n, sigma = 0.5)
    V_tau <- NA
    print(paste("Running simulation", i, "with n =", n))
    
  } else if (sim_type == "PM") {
    
    beta <- heterogeneity_degree  # Assuming heterogeneity_degree is used as beta for PM simulations
    
    data <- sim_PM(n = n, 
                   sigma = sigma, 
                   beta = beta, 
                   alternative_type = alternative_type, 
                   sim_scenario = sim_scenario,
                   randomized = randomized,
                   pi_hat = pi_hat)
    
    # Calculate V_tau based on simulation settings
    if (sim_scenario == "reduced" && alternative_type == "smooth") {
      V_tau <- beta^2 / 12
    } else if (sim_scenario == "reduced" && alternative_type == "rough") {
      V_tau <- beta^2 * 0.5 * (1 - exp(-50 * pi^2)) / 12
    } else if (sim_scenario == "HD10") {
      V_tau <- 0.5
    } else {
      stop("Invalid simulation scenario.")
    }  
    print(paste("Running simulation", i, "with n =", n))
    
  } else {
    stop("Invalid simulation type.")
  }
  
  ci_low_multifold <- NA
  ci_high_multifold <- NA
  coverage_multifold <- NA
  V_taun <- NA
  DukesP <- NA
  CRHP <- NA
  
  # Choose analysis type
  if ("SB" %in% analysis_type) {
    SB_CF <- crossfit_pseudoVCATE(data = data,
                                  SL.lib = SL.lib,
                                  randomized = randomized,
                                  SL_cv_folds = K,
                                  family = family, 
                                  hal = hal,
                                  crossfit_folds = 2,
                                  pi_hat = pi_hat,
                                  alpha = 0.025,
                                  Vmax = 0.05,
                                  gridN = 400,
                                  M = 2000,
                                  seed = NA)
    ci_low_multifold <- median(SB_CF$ci_low)
    ci_high_multifold <- median(SB_CF$ci_high)
    coverage_multifold <- (ci_low_multifold <= V_tau) & (ci_high_multifold >= V_tau)
    V_taun <- mean(SB_CF$V_taunk_star_hat)
  }
  
  
  if ("Dukes" %in% analysis_type) {
    Y <- data$Y
    A <- data$A
    W <- data[, -c(1, 2)]
    
    DukesP <- QuantTestMulti(y = Y, a = A, w = W, x = W, 
                             no.bs = 100,
                             randomized = randomized, 
                             SL.lib = SL.lib, 
                             K = K, 
                             family = family,
                             hal = hal)$p.val
    
  }

    
  if ("CRH" %in% analysis_type) {
    
    plan <- make_pair_folds(nrow(data), K = 2L)
    
    CRH_CF <- crossfit_embedded_u(
      data = data, 
      pair_folds = plan,
      family = family,
      inner_splits = 2L,
      K = K,
      SL.lib = SL.lib,
      randomized = randomized, 
      pi_hat = pi_hat,
      kernel = "median",
      bandwidth_pairs = 50000L,
      keep_pair_scores = TRUE
    )
    CRHP <- bootstrap_embedded_u(CRH_CF, B = 3999L, alpha = 0.05)$p_value
  }  
  
  
  return(list(n = n,
              V_tau = V_tau,
              sigma = sigma,
              V_taun = V_taun,
              coverage_multifold = coverage_multifold,
              ci_low_multifold = ci_low_multifold,
              ci_high_multifold = ci_high_multifold,
              DukesP = DukesP,
              CRHP = CRHP))
}


#-------------------------------------------------------------------------------
# 1. Set parameters for the simulation study
#-------------------------------------------------------------------------------




#-------------------------------------------------------------------------------
# 1.0 Simulation of data settings
#-------------------------------------------------------------------------------
# Set a random seed
seed <- sample(1:1e6, 1)

# Heterogeneity specific parameters
heterogeneity_degree <- c(0, 0.5, 1.0)    # Corresponding to beta values for PM simulations

# Number of observations for each simulation setting
n_values <- c(400, 800, 1200, 1600)

# Simulation type and settings
sim_type <- "PM"                      # Type of simulation ("SB", "Dukes", "Hines", "PM")
randomized <- TRUE
pi_hat <- 0.5


# PM specific settings
sigma <- c(0.1, 0.4, 0.8)             # Correlation parameter for the copula
alternative_type <- "rough"          # Type of alternative hypothesis ("smooth" or "rough")
sim_scenario <- "reduced"             # Simulation scenario ("reduced", "regular" or "HD10")


#-------------------------------------------------------------------------------
# 1.1 Analysis settings
#-------------------------------------------------------------------------------
# Analysis type and settings
analysis_type <- c("SB", "CRH", "Dukes")
K <- 2                                      # Number of folds for cross-fitting SuperLearners
hal <- FALSE                                # Whether to use HAL instead of SuperLearner for nuisance function estimation
SL.lib <- list("SL.glmnet", "SL.xgboost", "SL.ranger", "SL.mean")


#-------------------------------------------------------------------------------
# 1.2 Simulation study settings
#-------------------------------------------------------------------------------
n_rep <- 20L

settings <- expand.grid(
  n = n_values,
  sigma = sigma,
  heterogeneity_degree = heterogeneity_degree,
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
      
      run_one(
        i = sim_id,
        n = setting$n,
        sigma = setting$sigma,
        heterogeneity_degree = setting$heterogeneity_degree,
        sim_type = sim_type,
        randomized = randomized,
        pi_hat = pi_hat,
        seed = NA,
        J = J,
        alternative_type = alternative_type,
        sim_scenario = sim_scenario,
        K = K,
        family = gaussian(),
        SL.lib = SL.lib,
        hal = hal,
        heterogeneity = heterogeneity_type,
        analysis_type = analysis_type
      )
    },
    future.chunk.size = 1L,
    future.seed = seed,
    future.packages = "SuperLearner"
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
    "_sigma=", sigma,
    "_degree=", heterogeneity_degree
  )
)


# Process results
results_df <- lapply(results, function(sim_result) {
  matrix(unlist(sim_result), ncol = 9, byrow = T)
})
results_df <- data.frame(do.call(rbind, results_df))
colnames(results_df) <- c("n", 
                          "V_tau", 
                          "Sigma",
                          "V_taun", 
                          "coverage_multifold", 
                          "ci_low_multifold", 
                          "ci_high_multifold",
                          "DukesP",
                          "CRHP")



# Summarize results
summary_results <- results_df %>% 
  group_by(n, V_tau, Sigma) %>%                               # Change between group_by(V_tau) and group_by(n)
  summarise(SBLevel = mean(coverage_multifold),
            SBPower = mean(ci_low_multifold > 0),
            DukesLevel = mean(DukesP > 0.05 & V_tau == 0),
            DukesPower = mean(DukesP < 0.05 & V_tau > 0),
            CRHLevel = mean(CRHP > 0.05 & V_tau == 0),
            CRHPower = mean(CRHP < 0.05 & V_tau > 0),
            V_taun = mean(V_taun),
            ci_low_multifold = mean(ci_low_multifold),
            ci_high_multifold = mean(ci_high_multifold))
summary_results$reps <- n_rep



#Save results df as txt
#write.table(results_df, file = "/Users/christianhejstvig-larsen/Desktop/Studie/Master/Cambridge/Essay/R/Part-III-Essay-Repo/DukesSimulations16042026.txt", row.names = FALSE, col.names = TRUE, sep = "\t")









# # Plot results
# ggplot(summary_results, aes(x = n)) +
#   geom_line(aes(y = coverage), color = "blue") +
#   geom_line(aes(y = power), color = "red") +
#   geom_line(aes(y = DukesPower), color = "green") +
#   geom_line(aes(y = CRHPower), color = "purple") +
#   labs(title = "Coverage and Power vs V_tau", y = "Proportion", x = "V_tau") +
#   theme_minimal() +
#   theme(legend.position = "none") +
#   geom_line(aes(y = 0.95), color = "black", linetype = "dashed") +
#   geom_line(aes(y = 0.05), color = "black", linetype = "dashed")
# 
# 
# 


