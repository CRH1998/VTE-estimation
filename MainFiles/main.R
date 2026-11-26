






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
registerDoMC(60) #registerDoMC(46)

library(ranger)
library(xgboost)


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












