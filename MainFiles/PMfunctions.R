
### Function to compute nuisance parameters (propensity score and outcome regressions)
comp_nuisance_params = function(Y,A,X,newX,
                                V1,V2,
                                randomized, 
                                pi_hat,
                                nuisance_learner,
                                SL.lib=list("SL.glmnet", "SL.xgboost", "SL.ranger", "SL.mean")){
  
  ### Inputs:
  # Y - outcome
  # A - binary treatment 
  # X - vector of covariates
  # newX - separate vector covariates used to compute predictions
  # V1, V2 - subsets of covariates defining variable importance measure of interest (for LOO: V1=X and V2=X\X1; for KOO: V1=X1 and V2=null)
  # nuisance_learner - chosen model to compute nuisance parameters ("ranger", "HAL" or "SL")
  
  # Size of validation sets
  n_val=nrow(newX)
  
  # Fit the propensity score model
  if (randomized) {
    # If the treatment is randomized, use the known constant treatment probability pi_hat.
    gn_val <- rep(pi_hat, n_val)
  } else if (nuisance_learner=="HAL"){
    # Set the max_degree k for HAL
    k <- 2
    gn_model <- fit_hal(Y=A,X=X, family = "binomial", max_degree = k, fit_control = list(maxit = 10^6))
    gn_val <- predict(gn_model, new_data = newX)
  } else if (nuisance_learner=="ranger"){
    gn_model = ranger(y=A,x=X)
    gn_val = predict(gn_model, newX)$predict
  } else if (nuisance_learner=="SL"){
    SL.lib=SL.lib
    gn_model = SuperLearner(Y=A,X=X,SL.library=SL.lib,family='binomial',cvControl = list(V = 2))
    gn_val = as.vector(predict(gn_model, newdata = newX, onlySL = TRUE)$pred)
  }
  
  # Fit the outcome regression model
  if(nuisance_learner=="HAL"){
    # Set the max_degree k for HAL
    k <- 2
    Qn_model <- fit_hal(Y=Y,X=data.frame(X=X,A=A), max_degree = k, fit_control = list(maxit = 10^6))
    Qn_val <- predict(Qn_model, new_data = data.frame(X=rbind(newX,newX),A=rep(c(0,1),each=n_val)))
  } else if (nuisance_learner=="ranger"){
    Qn_model = ranger(y=Y,x=data.frame(X=X,A=A))
    Qn_val = predict(Qn_model, data.frame(X=rbind(newX,newX),A=rep(c(0,1),each=n_val)))$predict
  } else if (nuisance_learner=="SL"){
    SL.lib=SL.lib
    Qn_model <- SuperLearner(Y=Y,X=data.frame(X=X,A=A),SL.library=SL.lib, cvControl = list(V = 2))
    Qn_val <- as.vector(predict(Qn_model, newdata = data.frame(X=rbind(newX,newX),A=rep(c(0,1),each=n_val)), onlySL = TRUE)$pred)
  }
  
  # Predictions of the outcome under treatment and under no-treatment from the outcome regression model
  Q0n_val = Qn_val[1:n_val]
  Q1n_val = Qn_val[(n_val+1):(2*n_val)]
  
  my_list <- list(gn = gn_val,Q0n = Q0n_val,Q1n = Q1n_val)
  return(my_list)
}


# Function to generate the permutation of the set of variables for the Monte Carlo simulation of the Shapley value
sample_permutations = function(V, n_MC, test_variable){
  
  ### Inputs:
  # V - subset of covariates used for buiding permutations
  # n_MC - number of Monte Carlo simulations to approximate the Shapley value (n_MC = factorial(X_dim) computes the exact values)
  # test_variable - the covariate tested for variable importance (in our simulation study it is X1)
  
  # Vector of indices to permute
  index <- 1:length(V) 
  
  # If n_MC is set to factorial(length(V)) then we take all the permutation to compute Shapley values exactly. If n < factorial(length(V))
  # then we sample a subset of all permutation and these are used to approximate Shapley values. 
  if(n_MC==factorial(length(V))){
    permutations_sample <- t(permutations(length(index),length(index),index))
    V_permutations <- matrix(V[permutations_sample], length(index), n_MC)
  } else {
    permutations_sample <- Rfast::colShuffle(matrix(index, length(index), n_MC))
    V_permutations <- matrix(V[permutations_sample], length(index), n_MC)
  }
  
  V1_list <- list()
  V2_list <- list()
  
  # Find the position of specified test variable in a given permutation and compute 
  # vector V2 (containing all the covariates occurring in a given permutation before 
  # the test_variable) and vector V1 = V2 + X1.
  for (i in 1:n_MC){
    position <- match(test_variable, V_permutations[,i])
    if(position > 1){
      V2 <- V_permutations[1:(position - 1),i]
    } else {
      V2 <- NULL
    }
    
    # Construct V1
    V1 <- c(V2, test_variable)
    
    # Save V1 and V2 to a list
    V1_list[i] <- list(V1)
    V2_list[i] <- list(V2)
  }
  
  my_list <- list(V1_list=V1_list, V2_list=V2_list)
  return(my_list)
}


# Function to compute squared RKHS norm of the one-step estimator (formula S4 in the Appendix)
comp_one_step_RKHSnorm_squared = function(alphaV1, alphaV2, 
                                          betaV1, betaV2, 
                                          W_V1, W_V2, 
                                          K_V1, K_V2, 
                                          kernelMatrix, VIM_method, X_dim, n_MC){
  
  ### Input:
  # alphaV1, alphaV2 - variables as defined in display S1 in the Appendix,
  # betaV1, betaV2 - variables as defined in display S1 in the Appendix,
  # W_V1, W_V2 - matrices as defined in display S3 in the Appendix,
  # K_V1, K_V2 - kernal Gram matrices corresponding to subset of covariates V1 and V2 as defined in display S3 in the Appendix,
  # kernelMatrix - kernal Gram matrix corresponding to the full set of covariates X,
  # VIM_method - variable importance measure that one wishes to compute ("LOCO" or "Shapley")
  # X_dim - dimensionality of the covariates vector
  # n_MC - number of covariates permutations used to approximate Shapley values 
  
  n <- dim(kernelMatrix)[1]
  
  ##############################################################################
  ##### Compute the squared RKHS norm of the one-step estimator (expression S4 in the Appendix)
  ##############################################################################
  
  if (VIM_method=="LOCO"){
  
    # Compute intermediate term for the expression S4 in the Appendix (applied to the case of LOCO - see Section C.2)
    W_V1_K_V1_betaV1 <- W_V1 %*% (K_V1 %*% betaV1)
    W_V2_K_V2_betaV2 <- W_V2 %*% (K_V2 %*% betaV2)
    
    termV1 <- alphaV1 + W_V1_K_V1_betaV1
    termV2 <- alphaV2 + W_V2_K_V2_betaV2
    term <- termV1 - termV2
    
    # Compute expression S4 in the Appendix (applied to the case of LOCO - see Section C.2)
    one_step_RKHSnorm_squared <- ( t(term) %*% kernelMatrix %*% term ) / n^2
    
  } else if (VIM_method=="Shapley"){
    
    W_V1_K_V1_betaV1 <- c()
    W_V2_K_V2_betaV2 <- c()
    
    for (j in 1:n_MC){
      
      # Compute intermediate term for the expression S4 in the Appendix (applied to the case of Shapley - see Section C.4)
      W_V1_K_V1_betaV1[[j]] <- W_V1[[j]] %*% K_V1[[j]] %*% betaV1[[j]]
      W_V2_K_V2_betaV2[[j]] <- W_V2[[j]] %*% K_V2[[j]] %*% betaV2[[j]]
      
      }
  
  # Efficient computation of the expression S4 in the Appendix (to reduce the number of matrix multiplications)
  sum_alphaV1 <- Reduce(`+`, alphaV1) 
  sum_alphaV2 <- Reduce(`+`, alphaV2) 
  
  sum_W_V1_K_V1_betaV1 <- Reduce(`+`, W_V1_K_V1_betaV1) 
  sum_W_V2_K_V2_betaV2 <- Reduce(`+`, W_V2_K_V2_betaV2)
  
  termV1 <- sum_alphaV1 + sum_W_V1_K_V1_betaV1
  termV2 <- sum_alphaV2 + sum_W_V2_K_V2_betaV2
  term <- termV1 - termV2
  
  # Compute expression S4 in the Appendix (applied to the case of LOCO - see Section C.4)
  one_step_RKHSnorm_squared <- ( t(term) %*% kernelMatrix %*% term ) / (factorial(X_dim) * n)^2
  
  }
  return(one_step_RKHSnorm_squared)
}


# Function to compute zeta_hat via bootstrap (see Section D.1 in the Appendix)
comp_zeta_hat = function(alphaV1, alphaV2,
                                  betaV1, betaV2,
                                  W_V1, W_V2,
                                  K_V1, K_V2,
                                  kernelMatrix, n_boot, alpha,
                                  VIM_method, X_dim, n_MC){
  
  ### Input:
  # alphaV1, alphaV2 - variables as defined in display S1 in the Appendix,
  # betaV1, betaV2 - variables as defined in display S1 in the Appendix,
  # W_V1, W_V2 - matrices as defined in display S3 in the Appendix,
  # K_V1, K_V2 - kernal Gram matrices corresponding to subset of covariates V1 and V2 as defined in display S3 in the Appendix,
  # kernelMatrix - kernal Gram matrix corresponding to full set of covariates X,
  # n_boot - number of bootstrap samples,
  # alpha - level of the test,
  # VIM_method - variable importance measure that one wishes to compute ("LOCO" or "Shapley")
  # X_dim - dimensionality of the covariates vector
  # n_MC - number of covariates permutations used to approximate Shapley values 
  
  n <- nrow(kernelMatrix)
  boot <- numeric(n_boot)
  
  # Weights are given as a difference of the empirical distribution of bootstrapped sample and the empirical distribution of the original sample.
  wgt = matrix(c(rmultinom(n_boot,n,rep(1/n,n))-1) / n, nrow = n)
  
  ################################################################################
  ##### Compute the squared RKHS norm of the bootstrap approximation of the Gaussian element - see Section D.1 in the Appendix
  ################################################################################
  
  if (VIM_method == "LOCO"){
    
    ### Compute the squared RKHS norm of the bootstrap approximation of the Gaussian element 
    ### (applied to the case of LOCO - see Section D.2 in the Appendix).
    
    # Compute betas_V tilde (see definition in Section D.1 in the Appendix)
    betaV1tilde <- wgt * as.vector(betaV1)
    betaV2tilde <- wgt * as.vector(betaV2)
    
    # Compute alphas tilde (see definition in Section D.1 in the Appendix)
    alphaV1tilde <- wgt * as.vector(alphaV1)
    alphaV2tilde <- wgt * as.vector(alphaV2)
    
    # Compute intermediate term for squared RKHS norm of the bootstrap approximation of the Gaussian element 
    # (applied to the case of LOCO - see Section D.2 in the Appendix)
    W_V1_K_V1_betaV1tilde <- W_V1 %*% K_V1 %*% betaV1tilde
    W_V2_K_V2_betaV2tilde <- W_V2 %*% K_V2 %*% betaV2tilde
    
    termV1 <- alphaV1tilde + W_V1_K_V1_betaV1tilde
    termV2 <- alphaV2tilde + W_V2_K_V2_betaV2tilde
    term <- termV1 - termV2
    
    # Compute the squared RKHS norm of the bootstrap approximation of the Gaussian element 
    # (applied to the case of LOCO - see Section D.2 in the Appendix)
    result <- t(term) %*% kernelMatrix %*% term
    
  } else if (VIM_method == "Shapley"){
    
    ### Compute the squared RKHS norm of the bootstrap approximation of the Gaussian element 
    ### (applied to the case of Shapley values - see Section D.4 in the Appendix).
    
    ### Compute betas_V tilde (see definition in Section D.1 in the Appendix)
    alphaV1tilde <- c()
    alphaV2tilde <- c()
    betaV1tilde <- c()
    betaV2tilde <- c()
    
    W_V1_K_V1_betaV1tilde <- c()
    W_V2_K_V2_betaV2tilde <- c()
    
    for (k in 1:n_MC){
      
      ### Compute alphas_V tilde (see definition in Section D.1 in the Appendix)
      alphaV1tilde[[k]] <- wgt * as.vector(alphaV1[[k]])
      alphaV2tilde[[k]] <- wgt * as.vector(alphaV2[[k]])
      
      ### Compute betas_V tilde (see definition in Section D.1 in the Appendix)
      betaV1tilde[[k]] <- wgt * as.vector(betaV1[[k]])
      betaV2tilde[[k]] <- wgt * as.vector(betaV2[[k]])
      
      # Compute intermediate term for squared RKHS norm of the bootstrap approximation of the Gaussian element 
      # (applied to the case of Shapley - see Section D.4 in the Appendix)
      W_V1_K_V1_betaV1tilde[[k]] <- W_V1[[k]] %*% K_V1[[k]] %*% betaV1tilde[[k]]
      W_V2_K_V2_betaV2tilde[[k]] <- W_V2[[k]] %*% K_V2[[k]] %*% betaV2tilde[[k]]
      
    }
    
    # Efficient computation of the squared RKHS norm of the bootstrap approximation of the Gaussian element 
    # (to reduce the number of matrix multiplications) - see Section D.1 in the Appendix
    sum_alphaV1tilde <- Reduce(`+`, alphaV1tilde) 
    sum_alphaV2tilde <- Reduce(`+`, alphaV2tilde) 
    
    sum_W_V1_K_V1_betaV1tilde <- Reduce(`+`, W_V1_K_V1_betaV1tilde) 
    sum_W_V2_K_V2_betaV2tilde <- Reduce(`+`, W_V2_K_V2_betaV2tilde)
    
    termV1 <- sum_alphaV1tilde + sum_W_V1_K_V1_betaV1tilde
    termV2 <- sum_alphaV2tilde + sum_W_V2_K_V2_betaV2tilde
    term <- termV1 - termV2
    
    result <- ( t(term) %*% kernelMatrix %*% term ) / (factorial(X_dim))^2
    
  }
  
  boot <- diag(result) * n 
  
  # Compute the threshold estimate
  zeta_hat <- quantile(boot, 1-alpha)
  return(zeta_hat)
  
}


# The main function to compute the test of zero-importance
importance_test = function(V, test_variable, X_dim, n, n_boot, alpha,
                           n_MC, sigma, beta, alternative_type,
                           sim_data = FALSE, 
                           data,
                           VIM_method,
                           randomized = TRUE, pi_hat = 0.5,
                           nuisance_learner, SL.lib=list("SL.glmnet", "SL.xgboost", "SL.ranger", "SL.mean"), 
                           sim_scenario){
  #browser()
  ### Inputs:
  ## V - subset of covariates used for buiding permutations
  ## test_variable - the covariate tested for variable importance (in our simulation study it is X1)
  ## X_dim - dimensionality of the covariates vector
  ## n - sample size
  ## n_sims - number of simulations to compute the variable importance measure
  ## n_boot - number of bootstrap samples,
  ## alpha - level of the test,
  ## n_MC - number of covariates permutations used to approximate Shapley values (n_MC = factorial(X_dim) computes the exact values)
  ## sigma - correlatin between covariates
  ## beta - coefficient of X1 (which is a covariate whose importance we are testing in the simulation study) 
  # in the in the partially linear model generating tau (treatment effect). Parameter beta controls variable importance of X1.
  ## alternative type - specifies form of tau, i.e. treatment effect function, and sets the alternative type (eithr "smooth" or "rough")
  ## VIM_method - variable importance measure that one wishes to compute ("LOCO" or "Shapley")
  ## nuisance_learner - chosen model to compute nuisance parameters ("ranger", "HAL" or "SL")
  ## sim_scenario - sets the dimensionality of the covariates in a given set-up. We considered four set-ups:
  # "reduced" (X_dim = 3), "regular" (X_dim = 5), "HD10" (X_dim = 10) and "HD" (X_dim = 20).
  
  
  # Data frame for the output
  
  if (sim_data) {
  ################################################################################
  ##### Simulate the data
  ################################################################################
  
  # Simulate dataset A (used for cross-fitting)
  dt_setA <- sim.data(n=n/2, sigma=sigma, beta=beta, alternative_type=alternative_type, sim_scenario=sim_scenario)
  X_setA <- dt_setA$X
  A_setA <- dt_setA$A
  Y_setA <- dt_setA$Y
  tau_setA <- dt_setA$tau

  # Simulate dataset B (used for cross-fitting)
  dt_setB <- sim.data(n=n/2, sigma=sigma, beta=beta, alternative_type=alternative_type, sim_scenario=sim_scenario)
  X_setB <- dt_setB$X
  A_setB <- dt_setB$A
  Y_setB <- dt_setB$Y
  tau_setB <- dt_setB$tau

  # Combine datasets
  X <- rbind(X_setA, X_setB)
  A <- c(A_setA, A_setB)
  Y <- c(Y_setA, Y_setB)
  
  } else {
    # Use the provided data
    n <- nrow(data)
    
    Y <- data$Y
    A <- data$A
    X <- data[, V]
    
    # Split the data into two halves for cross-fitting
    n_half <- floor(nrow(data) / 2)
    X_setA <- X[1:n_half, ]
    A_setA <- A[1:n_half]
    Y_setA <- Y[1:n_half]
    
    X_setB <- X[(n_half + 1):nrow(data), ]
    A_setB <- A[(n_half + 1):nrow(data)]
    Y_setB <- Y[(n_half + 1):nrow(data)]
  }
  
  
  # Choose the bandwidth for the kernel based on the "median heuristic"
  bandwidth <- 2 * median(dist(X)^2)
  
  if (VIM_method=="LOCO"){
    
    ### Compute LOCO variable importance measure
    # Define the subsets of covariates V1 and V2
    # V1 <- V
    # if (length(V1)>1){
    #   V2 <- V1[!V1==test_variable]
    # } else {
    #   V2 <- NULL
    # }
    
    V1 <- V
    V2 <- setdiff(V1, test_variable)
    if (length(V2) == 0L) V2 <- NULL
    
    V1_setA <- X_setA[names(X_setA) %in% V1]
    V2_setA <- X_setA[names(X_setA) %in% V2]
    V1_setB <- X_setB[names(X_setB) %in% V1]
    V2_setB <- X_setB[names(X_setB) %in% V2]
    
    V1_vars <- rbind(V1_setA, V1_setB)
    V2_vars <- rbind(V2_setA, V2_setB)
    
    # Compute the Gram matrices for kernels K_V1 and K_V2 (for the subset of covariates V1 and V2)
    # and for the kernelMarnix - Gram matrix for the kernel for the full vector of covariates X
    kernelMatrix <- KRLS::gausskernel(X,sigma=bandwidth)
    kernelMatrixV1 <- KRLS::gausskernel(V1_vars,sigma=bandwidth)
    if(!identical(V2,NULL)){
      kernelMatrixV2 <- KRLS::gausskernel(V2_vars,sigma=bandwidth)
    } else {
      kernelMatrixV2 <- matrix(1,n,n)
    }
    
    ################################################################################
    ##### Compute nuisance parameters and pseudooutcome psi (as defined in Theorem 1 in the paper)
    ################################################################################
    
    nuisance_parameters_setA <- comp_nuisance_params(Y=Y_setB,A=A_setB,X=X_setB,newX=X_setA,V1=V1,V2=V2,nuisance_learner=nuisance_learner,SL.lib=SL.lib,randomized=randomized, pi_hat=pi_hat)
    nuisance_parameters_setB <- comp_nuisance_params(Y=Y_setA,A=A_setA,X=X_setA,newX=X_setB,V1=V1,V2=V2,nuisance_learner=nuisance_learner,SL.lib=SL.lib,randomized=randomized, pi_hat=pi_hat)
    
    # Compute the propensity scores and outcome predictions in the dataset A
    gn_setA <- nuisance_parameters_setA$gn
    Q1n_setA <- nuisance_parameters_setA$Q1n
    Q0n_setA <- nuisance_parameters_setA$Q0n
    
    # Compute the propensity scores and outcome predictions in the dataset B
    gn_setB <- nuisance_parameters_setB$gn
    Q1n_setB <- nuisance_parameters_setB$Q1n
    Q0n_setB <- nuisance_parameters_setB$Q0n
    
    gn <- c(gn_setA, gn_setB)
    Q1n <- c(Q1n_setA, Q1n_setB)
    Q0n <- c(Q0n_setA, Q0n_setB)
    
    ### Compute pseudooutcomes for CATE computation (psi as defined in Theorem 1 in the paper)
    psi <- A/gn * (Y - Q1n) - (1-A)/(1-gn) * (Y - Q0n) + Q1n - Q0n
    
    ################################################################################
    # Compute kappa (kappa as defined in Theorem 1 in the paper)
    ################################################################################
    
    # Construct the outcomes for the kernel ridge regression
    KRR_outcome_V1 <- kernelMatrixV1 %*% rep(1/n,n)
    KRR_outcome_V2 <- kernelMatrixV2 %*% rep(1/n,n)
    
    # Compute matrix W_V1 (defined in display S3 in the Appendix) using kernel ridge regression (KRR) unless V1 is full vector of covariates 
    # (in this case the estimator simplifies and there is no need to perform KRR)
    if(!setequal(V1,names(X))){
      
      # Compute matrix W_V1 using lambda obtained from kernel ridge regression (KRR)
      # It is advised to compute penalization parameter lambda using cross-validation, however for computational reasons we set lambda fixed 
      K_V1 <- kernelMatrixV1
      W_V1 <- chol2inv(chol((K_V1 + sqrt(log(n)/n) * diag(rep(1,n)))))
      
    } else {
      
      K_V1 <- diag(rep(1,n))
      W_V1 <- diag(rep(1,n))
      
    }
    
    # Compute matrix W_V2 (defined in display S3 in the Appendix) using kernel ridge regression (KRR) unless V2 is empty vector of covariates
    # (in this case the estimator simplifies and no need to perform KRR)
    if(!identical(V2,NULL)){
      
      # Compute matrix W_V2 using lambda obtained from kernel ridge regression (KRR)
      # It is advised to compute penalization parameter lambda using cross-validation, however for computational reasons we set lambda fixed
      K_V2 <- kernelMatrixV2
      W_V2 <- chol2inv(chol((K_V2 + sqrt(log(n)/n) * diag(rep(1,n)))))
      
    } else {
      
      K_V2 <- kernelMatrixV2
      W_V2 <- diag(rep(1,n))/n
      
    }
    
    ################################################################################
    ##### Compute the plug-in estimators for the nuisance parameters nuV1 and nuV2 using KRR
    ################################################################################
    
    # Compute matrix CATE function tau_n using kernel ridge regression (KRR)
    # It is advised to compute penalization parameter lambda using cross-validation, however for computational reasons we set lambda fixed 
    M <-  chol2inv(chol((kernelMatrix + sqrt(log(n)/n) * diag(rep(1,n)))))
    
    # Compute CATE function (which will be then regressed on V to obtain nuV1 and nuV2)
    tau_n <- kernelMatrix %*% M %*% psi
    
    ### Recompute the outcome predictions to be compatible with tau_n
    Q0n_new <- Q0n
    Q1n_new <- Q0n + tau_n
    
    ### Recompute pseudooutcome for CATE computation to be compatible with tau_n
    psi_new <- A/gn * (Y - Q1n_new) - (1-A)/(1-gn) * (Y - Q0n_new) + Q1n_new - Q0n_new
    
    # Compute the plug-in estimator for the nuisance paramenter nuV1 using KRR
    nuV1 <- K_V1 %*% W_V1 %*% kernelMatrix %*% M %*% psi_new
    
    # Compute kernel ridge regression unless V2 is empty vector of covariates 
    # (in this case the estimator simplifies and no need to perform KRR)
    if(!identical(V2,NULL)){
      nuV2 <- K_V2 %*% W_V2 %*% kernelMatrix %*% M %*% psi_new
    } else {
      nuV2 <- rep(mean(psi_new), length(psi_new))
    }
    
    # Compute betas (as defined in display S1 in the Appendix)
    betaV1 <- psi_new - nuV1
    betaV2 <- psi_new - nuV2
    
    # Compute alphas (as defined in display S1 in the Appendix)
    alphaV1 <- nuV1
    alphaV2 <- nuV2
    
    ################################################################################
    ##### Compute the squared RKHS norm of the one-step estimator (see display S4 in the Appendix)
    ################################################################################
    
    # Compute the RKHS norm of the one-step estimator
    one_step_RKHSnorm_squared0 <- as.numeric(comp_one_step_RKHSnorm_squared(alphaV1=alphaV1, alphaV2=alphaV2,
                                                                            betaV1=betaV1, betaV2=betaV2,
                                                                            W_V1=W_V1, W_V2=W_V2, 
                                                                            K_V1=K_V1, K_V2=K_V2,
                                                                            kernelMatrix=kernelMatrix,
                                                                            VIM_method=VIM_method, 
                                                                            X_dim=X_dim, n_MC=n_MC))
    
    ################################################################################
    ##### Compute the threshold zeta using bootstrap (see Section D1 in the Appendix)
    ################################################################################
    
    # Compute the threshold estimate using bootstrap
    zeta_hat0 <- as.numeric(comp_zeta_hat(alphaV1=alphaV1, alphaV2=alphaV2,
                                          betaV1=betaV1, betaV2=betaV2,
                                          W_V1=W_V1, W_V2=W_V2,
                                          K_V1=K_V1, K_V2=K_V2,
                                          kernelMatrix=kernelMatrix,
                                          n_boot=n_boot, alpha=alpha,
                                          VIM_method=VIM_method,
                                          X_dim=X_dim, n_MC=n_MC))
    
  } else if (VIM_method=="Shapley"){
    
    ### Compute Shapley variable importance measure
    
    # Compute the permutations of the covariates
    my_list <- sample_permutations(V=V,n_MC=n_MC,test_variable=test_variable)
    
    V1_list <- my_list$V1_list
    V2_list <- my_list$V2_list
    
    alphaV1_list <- c()
    alphaV2_list <- c()
    
    betaV1_list <- c()
    betaV2_list <- c()
    
    K_V1_list <- c()
    K_V2_list <- c()
    
    W_V1_list <- c()
    W_V2_list <- c()
    
    for (j in 1:n_MC){
      
      V1 <- V1_list[[j]]
      V2 <- V2_list[[j]]
      
      V1_setA <- X_setA[names(X_setA) %in% V1]
      V2_setA <- X_setA[names(X_setA) %in% V2]
      V1_setB <- X_setB[names(X_setB) %in% V1]
      V2_setB <- X_setB[names(X_setB) %in% V2]
      
      V1_vars <- rbind(V1_setA, V1_setB)
      V2_vars <- rbind(V2_setA, V2_setB)
      
      # Compute the Gram matrices for kernels K_V1 and K_V2 (for the subset of covariates V1 and V2)
      # and for the kernelMarnix - Gram matrix for the kernel for the full vector of covariates X
      kernelMatrix <- KRLS::gausskernel(X,sigma=bandwidth)
      kernelMatrixV1 <- KRLS::gausskernel(V1_vars,sigma=bandwidth)
      if(!identical(V2,NULL)){
        kernelMatrixV2 <- KRLS::gausskernel(V2_vars,sigma=bandwidth)
      } else {
        kernelMatrixV2 <- matrix(1,n,n)
      }
      
      ################################################################################
      ##### Compute nuisance parameters and pseudooutcome psi (as defined in Theorem 1 in the paper)
      ################################################################################
      
      nuisance_parameters_setA <- comp_nuisance_params(Y=Y_setB,A=A_setB,X=X_setB,newX=X_setA,V1=V1,V2=V2,nuisance_learner=nuisance_learner)
      nuisance_parameters_setB <- comp_nuisance_params(Y=Y_setA,A=A_setA,X=X_setA,newX=X_setB,V1=V1,V2=V2,nuisance_learner=nuisance_learner)
      
      # Compute the propensity scores and outcome predictions in the dataset A
      gn_setA <- nuisance_parameters_setA$gn
      Q1n_setA <- nuisance_parameters_setA$Q1n
      Q0n_setA <- nuisance_parameters_setA$Q0n
      
      # Compute the propensity scores and outcome predictions in the dataset B
      gn_setB <- nuisance_parameters_setB$gn
      Q1n_setB <- nuisance_parameters_setB$Q1n
      Q0n_setB <- nuisance_parameters_setB$Q0n
      
      gn <- c(gn_setA, gn_setB)
      Q1n <- c(Q1n_setA, Q1n_setB)
      Q0n <- c(Q0n_setA, Q0n_setB)
      
      ### Compute pseudooutcomes for CATE computation (psi as defined in Theorem 1 in the paper)
      psi <- A/gn * (Y - Q1n) - (1-A)/(1-gn) * (Y - Q0n) + Q1n - Q0n
      
      ################################################################################
      # Compute kappa (kappa as defined in Theorem 1 in the paper)
      ################################################################################
      
      # Construct the outcomes for the kernel ridge regression
      KRR_outcome_V1 <- kernelMatrixV1 %*% rep(1/n,n)
      KRR_outcome_V2 <- kernelMatrixV2 %*% rep(1/n,n)
      
      # Compute matrix W_V1 (defined in display S3 in the Appendix) using kernel ridge regression (KRR) unless V1 is full vector of covariates 
      # (in this case the estimator simplifies and there is no need to perform KRR)
      if(!setequal(V1,names(X))){
        
        # Compute matrix W_V1 using lambda obtained from kernel ridge regression (KRR)
        # It is advised to compute lambda using cross-validation, however for computational reasons we set lambda fixed  
        K_V1 <- kernelMatrixV1
        W_V1 <- chol2inv(chol((K_V1 + sqrt(log(n)/n) * diag(rep(1,n)))))
        
      } else {
        
        K_V1 <- diag(rep(1,n))
        W_V1 <- diag(rep(1,n))
        
      }
      
      # Compute matrix W_V2 (defined in display S3 in the Appendix) using kernel ridge regression (KRR) unless V2 is empty vector of covariates
      # (in this case the estimator simplifies and no need to perform KRR)
      if(!identical(V2,NULL)){
        
        # Compute matrix W_V2 using lambda obtained from kernel ridge regression (KRR)
        # It is advised to compute lambda using cross-validation, however for computational reasons we set lambda fixed 
        K_V2 <- kernelMatrixV2
        W_V2 <- chol2inv(chol((K_V2 + sqrt(log(n)/n) * diag(rep(1,n)))))
        
      } else {
        
        K_V2 <- kernelMatrixV2
        W_V2 <- diag(rep(1,n))/n
        
      }
      
      ################################################################################
      ##### Compute the plug-in estimators for the nuisance parameters nuV1 and nuV2 using KRR
      ################################################################################
      
      # Compute CATE function tau_n using kernel ridge regression (KRR)
      # It is advised to compute lambda using cross-validation, however for computational reasons we set lambda fixed 
      M <-  chol2inv(chol((kernelMatrix + sqrt(log(n)/n) * diag(rep(1,n)))))
      
      # Compute CATE function (which will be then regressed on V to obtain nuV1 and nuV2)
      tau_n <- kernelMatrix %*% M %*% psi
      
      ### Recompute the outcome predictions to be compatible with tau_n
      Q0n_new <- Q0n
      Q1n_new <- Q0n + tau_n
      
      ### Recompute pseudooutcome for CATE computation to be compatible with tau_n
      psi_new <- A/gn * (Y - Q1n_new) - (1-A)/(1-gn) * (Y - Q0n_new) + Q1n_new - Q0n_new
      
      # Compute the plug-in estimator for the nuisance paramenter nuV1 using KRR
      nuV1 <- K_V1 %*% W_V1 %*% kernelMatrix %*% M %*% psi_new
      
      # Compute kernel ridge regression unless V2 is empty vector of covariates 
      # (in this case the estimator simplifies and no need to perform KRR)
      if(!identical(V2,NULL)){
        nuV2 <- K_V2 %*% W_V2 %*% kernelMatrix %*% M %*% psi_new
      } else {
        nuV2 <- rep(mean(psi_new), length(psi_new))
      }
      
      # Compute betas (as defined in display S1 in the Appendix)
      betaV1 <- psi_new - nuV1
      betaV2 <- psi_new - nuV2
      
      # Compute alphas (as defined in display S1 in the Appendix)
      alphaV1 <- nuV1
      alphaV2 <- nuV2
      
      # Save computed objects in a list
      alphaV1_list[[j]] <- alphaV1
      alphaV2_list[[j]] <- alphaV2
      betaV1_list[[j]] <- betaV1
      betaV2_list[[j]] <- betaV2
      K_V1_list[[j]] <- K_V1
      K_V2_list[[j]] <- K_V2
      W_V1_list[[j]] <- W_V1
      W_V2_list[[j]] <- W_V2
      
    }
    
    ################################################################################
    ##### Compute the squared RKHS norm of the one-step estimator (see display S4 in the Appendix)
    ################################################################################
    
    # Compute the RKHS norm of the one-step estimator
    one_step_RKHSnorm_squared0 <- as.numeric(comp_one_step_RKHSnorm_squared(alphaV1=alphaV1_list, alphaV2=alphaV2_list,
                                                                            betaV1=betaV1_list, betaV2=betaV2_list,
                                                                            W_V1=W_V1_list, W_V2=W_V2_list,
                                                                            K_V1=K_V1_list, K_V2=K_V2_list,
                                                                            kernelMatrix=kernelMatrix,
                                                                            VIM_method=VIM_method,
                                                                            X_dim=X_dim, n_MC=n_MC))
    
    ################################################################################
    ##### Compute the threshold zeta using bootstrap (see Section D1 in the Appendix)
    ################################################################################
    
    # Compute the threshold estimate using bootstrap
    zeta_hat0 <- as.numeric(comp_zeta_hat(alphaV1=alphaV1_list, alphaV2=alphaV2_list,
                                          betaV1=betaV1_list, betaV2=betaV2_list,
                                          W_V1=W_V1_list, W_V2=W_V2_list,
                                          K_V1=K_V1_list, K_V2=K_V2_list,
                                          kernelMatrix=kernelMatrix,
                                          n_boot=n_boot, alpha=alpha,
                                          VIM_method=VIM_method,
                                          X_dim=X_dim, n_MC=n_MC))
    
  }
  
  # Store the estimates 
  res0 <- list("one_step_RKHSnorm_squared0" = one_step_RKHSnorm_squared0, "zeta_hat0" = zeta_hat0)
  
  return(res0)
}
