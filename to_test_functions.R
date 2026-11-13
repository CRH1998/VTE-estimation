# Testing simulations and nuisance estimation functions
source("nuisance_func.R")
source("simulations.R")
source("u_statistic_crossfit.R")
source("embedded_u_statistic.R")


test_data <- sim_SB(n = 400, V_tau = 0.05, tau = 0.5, J = 2)
train_data <- sim_SB(n = 400, V_tau = 0.05, tau = 0.5, J = 2)

#ests <- nuisance_param_func(test_data, train_data)
#ate_est(train_data = train_data, splits = 2)



train_data <- sim_PM(n = 4000, alternatives = "rough")
test_data <- sim_PM(n = 4000, alternatives = "rough")
set.seed(42)
fit <- nuisance_func(
  train_data = train_data,  # training complement of the outer evaluation fold
  inner_splits = 2,
  randomized = FALSE,
  pi_hat = 0.5,
  kernel = "median",
  bandwidth_pairs = 50000L,
  K = 5
)

X_ev <- test_data[, fit$x_names, drop = FALSE]
s <- fit$folds[[1]]

s$centred_kernel(X_ev, X_ev)          # centred kernel matrix
s$mu(1, X_ev)                         # outcome prediction under treatment
s$pi(X_ev)                            # treatment probability
s$cate(X_ev)                          # CATE
var(s$cate(X_ev))                     # VTE
s$ate                                 # marginal mean of the fitted CATE
s$centred_cate(X_ev)                  # CATE minus s$ate
s$embedded_cate(X_ev)                 # one embedding value per evaluation row
s$residual(a = test_data$A, y = test_data$Y, x = X_ev)
s$embedding_norm_sq



fit$average(function(s) {
  s$residual(a = test_data$A, y = test_data$Y, x = X_ev) *
    s$embedded_cate(X_ev)
})




source("nuisance_func.R")
source("simulations.R")
source("u_statistic_crossfit.R")
source("embedded_u_statistic.R")


# Simulation settings
randomized <- TRUE
pi_hat <- 0.5
SL.lib <- c("SL.ranger", "SL.xgboost", "SL.glmnet")
K <- 2L
family <- gaussian()

#test_data <- sim_SB(n = 400, V_tau = 0.5, tau = 0.15, J = 2)
#test_data <- sim_PM(n = 400, alternatives = "rough")
test_data <- sim_Dukes(800, heterogeneity = "QuantHMC", randomized = T, heterogeneity_degree = 3, pi_01x = 0.5)

#median_gaussian_kernel <- make_median_gaussian_kernel(test_data2[, c("X1", "X2", "X3")])

plan <- make_pair_folds(nrow(test_data), K = 2L)#, seed = 17L)

#set.seed(18)
out <- crossfit_embedded_u(
  family = family,
  test_data, 
  plan,
  inner_splits = 2L,
  K = K,
  randomized = randomized, 
  pi_hat = pi_hat,
  SL.lib = SL.lib,
  kernel = "median",
  bandwidth_pairs = 50000L,
  keep_pair_scores = TRUE
)
out$estimate

test <- bootstrap_embedded_u(out, B = 1999L, alpha = 0.05)
test[c("statistic", "critical_value", "p_value", "reject")]

head(test$draws)


SB_CF <- crossfit_pseudoVCATE(data = test_data,
                              SL.lib = SL.lib,
                              randomized = randomized,
                              SL_cv_folds = K,
                              family = family, 
                              hal = F,
                              crossfit_folds = 10,
                              pi_hat = pi_hat,
                              alpha = 0.025,
                              Vmax = 0.05,
                              gridN = 400,
                              M = 2000,
                              seed = NA)

ci_low_multifold <- median(SB_CF$ci_low)
ci_high_multifold <- median(SB_CF$ci_high)
c(ci_low_multifold, ci_high_multifold)



