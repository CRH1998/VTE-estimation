#-------------------------------------------------------------------------------
# 1.0 Simulation of data settings
#-------------------------------------------------------------------------------
# Set a random seed
#seed <- 21112026 #sample(1:1e6, 1)

# Heterogeneity specific parameters
heterogeneity_degree <- c(0, 0.7, 1.4) # (rough)    c(0, 0.5, 1.0) (smooth)    Corresponding to beta values for PM simulations

# Number of observations for each simulation setting
n_values <- c(400, 800, 1200, 1600)

# Simulation type and settings
sim_type <- "PM"                      # Type of simulation ("SB", "Dukes", "Hines", "PM")
randomized <- TRUE
pi_hat <- 0.5


# PM specific settings
sigma <- c(0.8)# c(0.1, 0.4, 0.8)# Correlation parameter for the copula
alternative_type <- "rough"          # Type of alternative hypothesis ("smooth" or "rough")
sim_scenario <- "reduced"            # Simulation scenario ("reduced", "regular" or "HD10")


#-------------------------------------------------------------------------------
# 1.1 Analysis settings
#-------------------------------------------------------------------------------
# Analysis type and settings
analysis_type <- c("SB", "CRH", "PM")
K <- 2                                      # Number of folds for cross-fitting SuperLearners
hal <- FALSE                                # Whether to use HAL instead of SuperLearner for nuisance function estimation
SL.lib <- list("SL.glmnet", "SL.xgboost", "SL.ranger", "SL.mean")


#-------------------------------------------------------------------------------
# 1.2 Simulation study settings
#-------------------------------------------------------------------------------
n_rep <- 30L

settings <- expand.grid(
  n = n_values,
  sigma = sigma,
  heterogeneity_degree = heterogeneity_degree,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

# Each setting gets n_rep consecutive runs
setting_index <- rep(seq_len(nrow(settings)), each = n_rep)







#install.packages("doMC")
#install.packages("mets")
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  BLIS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

library(RhpcBLASctl)
blas_set_num_threads(1)
omp_set_num_threads(1)

library(data.table)
setDTthreads(1)

library(doRNG)
library(doMC)
registerDoMC(4) #registerDoMC(46)




source("SBfunctions.R")     # Sánchez-Becerra (2023)
source("np-qualint-JASA.R") # Dukes et al. (2026)
source("PMfunctions.R")     # Morzywołek et al. (2026)
source("nuisance_func.R")
source("embedded_u_statistic.R")
source("u_statistic_crossfit.R")
source("simulations.R")
source("simulation_study.R")



resl <- foreach(i = seq_along(setting_index)) %dopar% {
  
  # Extract the simulation setting
  setting <- settings[setting_index[i], , drop = FALSE]
  
  # Repetition number within the setting
  sim_id <- (i - 1L) %% n_rep + 1L
  
  run_one(
    i = sim_id,
    n = setting$n,
    sigma = setting$sigma,
    heterogeneity_degree = setting$heterogeneity_degree,
    sim_type = sim_type,
    randomized = randomized,
    pi_hat = pi_hat,
    seed = NA,
    alternative_type = alternative_type,
    sim_scenario = sim_scenario,
    K = K,
    family = gaussian(),
    SL.lib = SL.lib,
    hal = hal,
    analysis_type = analysis_type
  )
}

save(resl, file = "resl.rda")












