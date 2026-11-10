# Testing simulations and nuisance estimation functions
source("nuisance_func.R")
source("simulations.R")
source("u_statistic_crossfit.R")
source("embedded_u_statistic.R")


test_data <- sim_SB(n = 400, V_tau = 0.05, tau = 0.5, J = 2)
train_data <- sim_SB(n = 400, V_tau = 0.05, tau = 0.5, J = 2)

#ests <- nuisance_param_func(test_data, train_data)
#ate_est(train_data = train_data, splits = 2)




set.seed(42)
fit <- nuisance_func(
  train_data = train_data,  # training complement of the outer evaluation fold
  splits = 2,
  randomized = FALSE,
  pi_hat = 0.5,
  K = 5
)

X_ev <- test_data[, fit$x_names, drop = FALSE]
s <- fit$folds[[1]]

s$centred_kernel(X_ev, X_ev)          # centred kernel matrix
s$mu(1, X_ev)                         # outcome prediction under treatment
s$pi(X_ev)                            # treatment probability
s$cate(X_ev)                          # CATE
s$ate                                 # marginal mean of the fitted CATE
s$centred_cate(X_ev)                  # CATE minus s$ate
s$embedded_cate(X_ev)                 # one embedding value per evaluation row
s$residual(a = test_data$A, y = test_data$Y, x = X_ev)
s$embedding_norm_sq



fit$average(function(s) {
  s$residual(a = test_data$A, y = test_data$Y, x = X_ev) *
    s$embedded_cate(X_ev)
})





# dat contains Y, A, and the numeric covariates.
plan <- make_pair_folds(nrow(test_data), K = 2L, seed = 17L)

out <- crossfit_embedded_u(
  data       = test_data,
  pair_folds = plan,
  splits     = 2L,       # inner splits for regressions and marginal
  K          = 5L,       # learner's internal cross-validation folds
  randomized = TRUE,
  pi_hat     = 0.5,
  SL.lib = c("SL.glmnet", "SL.randomForest", "SL.xgboost")
)

out$estimate




source("nuisance_func.R")
source("simulations.R")
source("u_statistic_crossfit.R")
source("embedded_u_statistic.R")

test_data <- sim_SB(n = 400, V_tau = 0.05, tau = 0.15, J = 4)
#test_data <- sim_Dukes(800, heterogeneity = "QuantHMC", randomized = T, heterogeneity_degree = 7.5)

#median_gaussian_kernel <- make_median_gaussian_kernel(test_data2[, c("X1", "X2", "X3")])

plan <- make_pair_folds(nrow(test_data), K = 2L, seed = 17L)

#set.seed(18)
out <- crossfit_embedded_u(
  test_data, 
  plan,
  inner_splits = 2L,
  K = 2L,
  randomized = TRUE, 
  pi_hat = 0.5,
  SL.lib = c("SL.ranger", "SL.xgboost", "SL.glmnet"),
  kernel = "median",
  bandwidth_pairs = 50000L,
  keep_pair_scores = TRUE
)
out$estimate

test <- bootstrap_embedded_u(out, B = 1999L, alpha = 0.05)

test[c("statistic", "critical_value", "p_value", "reject")]
head(test$draws)




