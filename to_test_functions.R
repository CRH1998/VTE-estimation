# Testing simulations and nuisance estimation functions



test_data <- sim_SB(n = 1000, V_tau = 0.5, tau = 0.5)
train_data <- sim_SB(n = 1000, V_tau = 0.5, tau = 0.5)

#ests <- nuisance_param_func(test_data, train_data)
#ate_est(train_data = train_data, splits = 2)



source("nuisance_func.R")

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

s$mu(1, X_ev)                         # outcome prediction under treatment
s$pi(X_ev)                            # treatment probability
s$cate(X_ev)                          # CATE
s$ate                                 # marginal mean of the fitted CATE
s$centred_cate(X_ev)                  # CATE minus s$ate
s$embedded_cate(X_ev)                 # one embedding value per evaluation row
s$residual(a = eval_data$A, y = eval_data$Y, x = X_ev)



fit$average(function(s) {
  s$residual(a = eval_data$A, y = eval_data$Y, x = X_ev) *
    s$embedded_cate(X_ev)
})




