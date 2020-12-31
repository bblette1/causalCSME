rm(list = ls())
library(geex)
library(matrixStats)
library(rootSolve)

nsims <- 500
n <- 2000
beta1_true <- 0.7
beta3_true <- -0.4
beta5_true <- 0.3
L1_prob <- 0.5
L2_mean <- 1
true_effect <- beta1_true + beta3_true*L1_prob + beta5_true*L2_mean
sigma_me <- 0.09


simulator <- function(trial, sigma_me) {
  
  # Outcome regression wrong, propensity model correct
  
  #V1 <- rbinom(n, 1, 0.5)
  #V2 <- rnorm(n, 0.75, 0.25)
  L1 <- rbinom(n, 1, L1_prob)
  L2 <- rnorm(n, L2_mean, 0.5)
  #X <- rnorm(n, 2 + 1.5*L1 - 1.0*L2, 0.9)
  #Y_logit <- -1 + beta1_true*X + 2*L1 + beta3_true*X*L1 -
    #0.4*L2 + beta5_true*X*L2
  X <- rnorm(n, 2 + 0.4*L1 - 0.5*L2, 0.9)
  Y_logit <- -2 + beta1_true*X - 0.4*L1 + beta3_true*X*L1 -
    0.2*L2 + beta5_true*X*L2 + rnorm(n, 0, 0.2)
  Y_prob <- exp(Y_logit) / (1 + exp(Y_logit))
  Y <- rbinom(n, 1, Y_prob)
  Xstar <- X + rnorm(n, 0, sqrt(sigma_me))
  data <- data.frame("Y" = Y, "Xstar" = Xstar, "L1" = L1, "L2" = L2)
  
  # Fit logistic regression to use for starting values
  mod <- glm(Y ~ Xstar*L2, family = "binomial")
  
  # Model 1: G-formula-CSME
  eefun_csme_gform <- function(data) {
    Y <- data$Y
    Xstar <- data$Xstar
    L2 <- data$L2
    delta <- function(beta1, beta3) {
      Xstar + sigma_me*(beta1 + beta3*L2)*Y
    }
    H <- function(x) {
      1 / (1 + exp(-x))
    }
    condexp <- function(beta0, beta1, beta2, beta3) {
      H(beta0 + (beta1 + beta3*L2)*delta(beta1, beta3) +
          beta2*L2 - ((beta1 + beta3*L2)^2)*sigma_me / 2)
    }
    function(theta) {
      c((Y - condexp(theta[1], theta[2], theta[3], theta[4])),
        (Y - condexp(theta[1], theta[2], theta[3], theta[4]))*L2,
        (Y - condexp(theta[1], theta[2], theta[3], theta[4]))*
          delta(theta[2], theta[4]),
        (Y - condexp(theta[1], theta[2], theta[3], theta[4]))*
          L2*delta(theta[2], theta[4]),
        L2 - theta[5],
        theta[2] + theta[4]*theta[5] - theta[6]
      )
    }
  }
  
  results_csme_gform <- m_estimate(estFUN = eefun_csme_gform, data = data,
                                   compute_roots = TRUE,
                                   root_control = 
                                     setup_root_control(start = c(coef(mod),
                                                                  mean(L2),
                                                                  0.7)))
  
  bias_gform_ps <- coef(results_csme_gform)[6] - true_effect
  se_gform_ps <- sqrt(vcov(results_csme_gform)[6, 6])
  coverage_gform_ps <- 
    1*(coef(results_csme_gform)[6] - 1.96*se_gform_ps < true_effect &
         coef(results_csme_gform)[6] + 1.96*se_gform_ps > true_effect)
  
  # Weighted CSME
  # Estimate weights
  denom_mod <- lm(Xstar ~ L1 + L2)
  p_denom <- predict(denom_mod, type='response')
  dens_denom <- dnorm(Xstar, p_denom, summary(denom_mod)$sigma)
  num_mod <- lm(Xstar ~ 1)
  p_num <- predict(num_mod, type='response')
  dens_num <- dnorm(Xstar, p_num, summary(denom_mod)$sigma)
  data$sw <- dens_num / dens_denom
  
  # Fit weighted regression for starting values
  wmod <- glm(Y ~ Xstar, weights = sw, family = "binomial", data = data)
  
  # Model 2: IPW-CSE
  # Get point estimates and variance using geex
  # First point estimates
  eefun_ipw1 <- function(data) {
    Y <- data$Y
    Xstar <- data$Xstar
    sw <- data$sw
    delta <- function(beta1) {
      Xstar + sigma_me*beta1*Y
    }
    H <- function(x) {
      1 / (1 + exp(-x))
    }
    condexp <- function(beta0, beta1) {
      H(beta0 + beta1*delta(beta1) - (beta1^2)*sigma_me/2)
    }
    function(theta) {
      c(sw*(Y - condexp(theta[1], theta[2])),
        sw*(Y - condexp(theta[1], theta[2]))*
          delta(theta[2])
      )
    }
  }
  
  results_ipw <- 
    m_estimate(estFUN = eefun_ipw1, data = data,
               root_control = setup_root_control(start = c(coef(wmod))))
  bias_ipw_ps <- coef(results_ipw)[2] - true_effect
  
  # Then get variance estimates accounting for weight estimation
  eefun_ipw2 <- function(data, model1, model2) {
    Y <- data$Y
    Xstar <- model.response(model.frame(model1, data = data))
    Lmat <- model.matrix(model1, data = data)
    Vmat <- model.matrix(model2, data = data)
    #sw <- data$sw
    n <- dim(data)[1]
    delta <- function(beta1) {
      Xstar + sigma_me*beta1*Y
    }
    H <- function(x) {
      1 / (1 + exp(-x))
    }
    condexp <- function(beta0, beta1) {
      H(beta0 + beta1*delta(beta1) - (beta1^2)*sigma_me/2)
    }
    function(theta) {
      p  <- length(theta)
      p1 <- length(coef(model1))
      p2 <- length(coef(model2))
      rho1 <- Lmat %*% theta[1:p1]
      rho2 <- Vmat %*% theta[(p1+1):(p1+p2)]
      
      dens_denom <- dnorm(Xstar, rho1, sqrt(theta[p-2]))
      dens_num <- dnorm(Xstar, rho2, sqrt(theta[p-2]))
      sw <- dens_num / dens_denom

      score_eqns1 <- apply(Lmat, 2, function(x) sum((Xstar - rho1) * x))
      score_eqns2 <- apply(Vmat, 2, function(x) sum((Xstar - rho2) * x))
      
      c(score_eqns1,
        score_eqns2,
        (n-p1)/n * theta[p-2] - (Xstar - Lmat %*% theta[1:p1])^2,
        sw*(Y - condexp(theta[p-1], theta[p])),
        sw*(Y - condexp(theta[p-1], theta[p]))*
          delta(theta[p])
      )
    }
  }
  
  results_ipw2 <- 
    m_estimate(estFUN = eefun_ipw2, data = data,
               outer_args = list(model1 = denom_mod, model2 = num_mod),
               compute_roots = FALSE,
               roots = c(coef(denom_mod), coef(num_mod), sigma(denom_mod)^2,
                         coef(results_ipw)))

  se_ipw_ps <- sqrt(vcov(results_ipw2)[7, 7])
  coverage_ipw_ps <- 
    1*(coef(results_ipw2)[7] - 1.96*se_ipw_ps < true_effect &
         coef(results_ipw2)[7] + 1.96*se_ipw_ps > true_effect)
  
  # Model 3: AIPW
  # Weighted version of model 1
  eefun_csme_aipw <- function(data, model) {
    Y <- data$Y
    Xstar <- data$Xstar
    L2 <- data$L2
    Lmat <- model.matrix(model, data = data)
    sw <- data$sw
    delta <- function(beta1, beta3) {
      Xstar + sigma_me*(beta1 + beta3*L2)*Y
    }
    H <- function(x) {
      1 / (1 + exp(-x))
    }
    condexp <- function(beta0, beta1, beta2, beta3) {
      H(beta0 + (beta1 + beta3*L2)*delta(beta1, beta3) +
          beta2*L2 - ((beta1 + beta3*L2)^2)*sigma_me / 2)
    }
    function(theta) {
      p  <- length(theta)
      p1 <- length(coef(model))
      rho <- Lmat %*% theta[1:p1]
      
      score_eqns <- apply(Lmat, 2, function(x) sum((Xstar - rho) * x))
      
      c(score_eqns,
        Xstar - theta[p-6],
        sw*(Y - condexp(theta[p-5], theta[p-4], theta[p-3], theta[p-2])),
        sw*(Y - condexp(theta[p-5], theta[p-4], theta[p-3], theta[p-2]))*L2,
        sw*(Y - condexp(theta[p-5], theta[p-4], theta[p-3], theta[p-2]))*
          delta(theta[p-4], theta[p-2]),
        sw*(Y - condexp(theta[p-5], theta[p-4], theta[p-3], theta[p-2]))*
          L2*delta(theta[p-4], theta[p-2]),
        L2 - theta[p-1],
        theta[p-4] + theta[p-2]*theta[p-1] - theta[p]
      )
    }
  }
  
  failed <- TRUE
  j <- 1
  
  while(failed == TRUE & j < 6) {
    
    failed <- FALSE
    startvec <- c(coef(denom_mod), mean(Xstar), coef(mod), mean(L2), 0.7)*(j == 1) +
      c(coef(denom_mod), mean(Xstar), rnorm(4, 0, j/5), mean(L2), 0.7)*(j > 1) 
    results_csme_aipw <- tryCatch(m_estimate(estFUN = eefun_csme_aipw, data = data,
                                             outer_args = list(denom_mod),
                                       root_control = setup_root_control(start = startvec)),
                            error = function(e) { failed <<- TRUE})
    if (failed == FALSE) {
      if (abs(coef(results_csme_aipw)[6]) > 2) { failed <- TRUE }
    }
    j <- j + 1
    
  }
  
  bias_aipw_ps <- coef(results_csme_aipw)[10] - true_effect
  se_aipw_ps <- sqrt(vcov(results_csme_aipw)[10, 10])
  coverage_aipw_ps <- 
    1*(coef(results_csme_aipw)[10] - 1.96*se_aipw_ps < true_effect &
         coef(results_csme_aipw)[10] + 1.96*se_aipw_ps > true_effect)
  
  if(failed == TRUE) { 
    bias_aipw_ps <- NA
    se_aipw_ps <- NA
    coverage_aipw_ps <- NA
  }


  # Outcome regression correct, propensity model wrong

  
  # Fit logistic regression to use for starting values
  mod <- glm(Y ~ Xstar*L1 + Xstar*L2, family = "binomial")
  
  # Model 1: G-formula-CSME
  eefun_csme_gform <- function(data) {
    Y <- data$Y
    Xstar <- data$Xstar
    L1 <- data$L1
    L2 <- data$L2
    delta <- function(beta1, beta4, beta5) {
      Xstar + sigma_me*(beta1 + beta4*L1 + beta5*L2)*Y
    }
    H <- function(x) {
      1 / (1 + exp(-x))
    }
    condexp <- function(beta0, beta1, beta2, beta3,
                        beta4, beta5) {
      H(beta0 + (beta1 + beta4*L1 + beta5*L2)*
         delta(beta1, beta4, beta5) + beta2*L1 + beta3*L2 -
          ((beta1 + beta4*L1 + beta5*L2)^{2})*sigma_me / 2)
    }
    function(theta) {
      c((Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6])),
        (Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6]))*L1,
        (Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6]))*L2,
        (Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6]))*
          delta(theta[2], theta[5], theta[6]),
        (Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6]))*
          L1*delta(theta[2], theta[5], theta[6]),
        (Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6]))*
          L2*delta(theta[2], theta[5], theta[6]),
        L1 - theta[7],
        L2 - theta[8],
        theta[2] + theta[5]*theta[7] + theta[6]*theta[8] - theta[9]
      )
    }
  }
  
  results_csme_gform <- m_estimate(estFUN = eefun_csme_gform, data = data,
                                   compute_roots = TRUE,
                                   root_control = 
                                     setup_root_control(start = c(coef(mod),
                                                                  mean(L1),
                                                                  mean(L2),
                                                                  0.7)))
  
  bias_gform_or <- coef(results_csme_gform)[9] - true_effect
  se_gform_or <- sqrt(vcov(results_csme_gform)[9, 9])
  coverage_gform_or <- 
    1*(coef(results_csme_gform)[9] - 1.96*se_gform_or < true_effect &
         coef(results_csme_gform)[9] + 1.96*se_gform_or > true_effect)
  
  # Weighted CSME
  # Estimate weights
  denom_mod <- lm(Xstar ~ L2)
  p_denom <- predict(denom_mod, type='response')
  dens_denom <- dnorm(Xstar, p_denom, summary(denom_mod)$sigma)
  num_mod <- lm(Xstar ~ 1)
  p_num <- predict(num_mod, type='response')
  dens_num <- dnorm(Xstar, p_num, summary(denom_mod)$sigma)
  data$sw <- dens_num / dens_denom
  
  # Fit weighted regression for starting values
  wmod <- glm(Y ~ Xstar, weights = data$sw, family = "binomial")
  
  # Model 2: IPW-CSE
  # Get point estimates and variance using geex
  eefun_ipw <- function(data, model) {
    Y <- data$Y
    Xstar <- data$Xstar
    Lmat <- model.matrix(model, data = data)
    sw <- data$sw
    delta <- function(beta1) {
      Xstar + sigma_me*beta1*Y
    }
    H <- function(x) {
      1 / (1 + exp(-x))
    }
    condexp <- function(beta0, beta1) {
      H(beta0 + beta1*delta(beta1) - (beta1^2)*sigma_me/2)
    }
    function(theta) {
      p  <- length(theta)
      p1 <- length(coef(model))
      rho <- Lmat %*% theta[1:p1]

      score_eqns <- apply(Lmat, 2, function(x) sum((Xstar - rho) * x))

      c(score_eqns,
        Xstar - theta[p-2],
        sw*(Y - condexp(theta[p-1], theta[p])),
        sw*(Y - condexp(theta[p-1], theta[p]))*
          delta(theta[p])
      )
    }
  }
  
  failed <- TRUE
  j <- 1
  
  while(failed == TRUE & j < 6) {
    
    failed <- FALSE
    startvec <- c(coef(denom_mod), mean(Xstar), coef(wmod)[1:2])*(j == 1) +
      c(coef(denom_mod), mean(Xstar), coef(mod)[1:2])*(j == 2) +
      c(coef(denom_mod), mean(Xstar), rnorm(2, 0, j/5))*(j > 2) 
    results_ipw <- tryCatch(m_estimate(estFUN = eefun_ipw, data = data,
                                       outer_args = list(denom_mod),
                                       root_control = setup_root_control(start = startvec)),
                            error = function(e) { failed <<- TRUE})
    if (failed == FALSE) {
      if (abs(coef(results_ipw)[2]) > 2) { failed <- TRUE }
    }
    j <- j + 1
    
  }
  
  bias_ipw_or <- coef(results_ipw)[5] - true_effect
  se_ipw_or <- sqrt(vcov(results_ipw)[5, 5])
  coverage_ipw_or <- 
    1*(coef(results_ipw)[5] - 1.96*se_ipw_or < true_effect &
         coef(results_ipw)[5] + 1.96*se_ipw_or > true_effect)
  
  # Model 3: AIPW
  # Include denominator density in g-formula CSME outcome model
  eefun_csme_aipw <- function(data, model) {
    Y <- data$Y
    Xstar <- data$Xstar
    L1 <- data$L1
    L2 <- data$L2
    Lmat <- model.matrix(model, data = data)
    sw <- data$sw
    delta <- function(beta1, beta4, beta5) {
      Xstar + sigma_me*(beta1 + beta4*L1 + beta5*L2)*Y
    }
    H <- function(x) {
      1 / (1 + exp(-x))
    }
    condexp <- function(beta0, beta1, beta2, beta3,
                        beta4, beta5) {
      H(beta0 + (beta1 + beta4*L1 + beta5*L2)*
          delta(beta1, beta4, beta5) + beta2*L1 + beta3*L2 -
          ((beta1 + beta4*L1 + beta5*L2)^{2})*sigma_me / 2)
    }
    function(theta) {
      p  <- length(theta)
      p1 <- length(coef(model))
      rho <- Lmat %*% theta[1:p1]

      score_eqns <- apply(Lmat, 2, function(x) sum((Xstar - rho) * x))

      c(score_eqns,
        Xstar - theta[p-9],
        sw*(Y - condexp(theta[p-8], theta[p-7], theta[p-6], theta[p-5],
                        theta[p-4], theta[p-3])),
        sw*(Y - condexp(theta[p-8], theta[p-7], theta[p-6], theta[p-5],
                        theta[p-4], theta[p-3]))*L1,
        sw*(Y - condexp(theta[p-8], theta[p-7], theta[p-6], theta[p-5],
                        theta[p-4], theta[p-3]))*L2,
        sw*(Y - condexp(theta[p-8], theta[p-7], theta[p-6], theta[p-5],
                        theta[p-4], theta[p-3]))*
          delta(theta[p-7], theta[p-4], theta[p-3]),
        sw*(Y - condexp(theta[p-8], theta[p-7], theta[p-6], theta[p-5],
                        theta[p-4], theta[p-3]))*
          L1*delta(theta[p-7], theta[p-4], theta[p-3]),
        sw*(Y - condexp(theta[p-8], theta[p-7], theta[p-6], theta[p-5],
                     theta[p-4], theta[p-3]))*
          L2*delta(theta[p-7], theta[p-4], theta[p-3]),
        L1 - theta[p-2],
        L2 - theta[p-1],
        theta[p-7] + theta[p-4]*theta[p-2] + theta[p-3]*theta[p-1] - theta[p]
      )
    }
  }
  
  results_csme_aipw <- m_estimate(estFUN = eefun_csme_aipw, data = data,
                                  compute_roots = TRUE,
                                  outer_args = list(denom_mod),
                                  root_control = 
                                    setup_root_control(start = c(coef(denom_mod),
                                                                 mean(Xstar),
                                                                 coef(mod),
                                                                 mean(L1),
                                                                 mean(L2),
                                                                 0.7)))
  
  bias_aipw_or <- coef(results_csme_aipw)[12] - true_effect
  se_aipw_or <- sqrt(vcov(results_csme_aipw)[12, 12])
  coverage_aipw_or <- 
    1*(coef(results_csme_aipw)[12] - 1.96*se_aipw_or < true_effect &
         coef(results_csme_aipw)[12] + 1.96*se_aipw_or > true_effect)
  

  # Both models correct

  
  # Fit linear regression to use for starting values
  mod <- glm(Y ~ Xstar*L1 + Xstar*L2, family = "binomial")
  
  # Model 1: G-formula-CSME
  eefun_csme_gform <- function(data) {
    Y <- data$Y
    Xstar <- data$Xstar
    L1 <- data$L1
    L2 <- data$L2
    delta <- function(beta1, beta4, beta5) {
      Xstar + sigma_me*(beta1 + beta4*L1 + beta5*L2)*Y
    }
    H <- function(x) {
      1 / (1 + exp(-x))
    }
    condexp <- function(beta0, beta1, beta2, beta3,
                        beta4, beta5) {
      H(beta0 + (beta1 + beta4*L1 + beta5*L2)*
          delta(beta1, beta4, beta5) + beta2*L1 + beta3*L2 -
          ((beta1 + beta4*L1 + beta5*L2)^{2})*sigma_me / 2)
    }
    function(theta) {
      c((Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6])),
        (Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6]))*L1,
        (Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6]))*L2,
        (Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6]))*
          delta(theta[2], theta[5], theta[6]),
        (Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6]))*
          L1*delta(theta[2], theta[5], theta[6]),
        (Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6]))*
          L2*delta(theta[2], theta[5], theta[6]),
        L1 - theta[7],
        L2 - theta[8],
        theta[2] + theta[5]*theta[7] + theta[6]*theta[8] - theta[9]
      )
    }
  }
  
  results_csme_gform <- m_estimate(estFUN = eefun_csme_gform, data = data,
                                   compute_roots = TRUE,
                                   root_control = 
                                     setup_root_control(start = c(coef(mod),
                                                                  mean(L1),
                                                                  mean(L2),
                                                                  0.7)))
  
  bias_gform_psor <- coef(results_csme_gform)[9] - true_effect
  se_gform_psor <- sqrt(vcov(results_csme_gform)[9, 9])
  coverage_gform_psor <- 
    1*(coef(results_csme_gform)[9] - 1.96*se_gform_psor < true_effect &
         coef(results_csme_gform)[9] + 1.96*se_gform_psor > true_effect)
  
  # Weighted CSME
  # Estimate weights
  denom_mod <- lm(Xstar ~ L1 + L2)
  p_denom <- predict(denom_mod, type='response')
  dens_denom <- dnorm(Xstar, p_denom, summary(denom_mod)$sigma)
  num_mod <- lm(Xstar ~ 1)
  p_num <- predict(num_mod, type='response')
  dens_num <- dnorm(Xstar, p_num, summary(denom_mod)$sigma)
  data$sw <- dens_num / dens_denom
  
  # Fit weighted regression for starting values
  wmod <- glm(Y ~ Xstar, weights = data$sw, family = "binomial")
  
  # Model 2: IPW-CSE
  # Get point estimates and variance using geex
  eefun_ipw <- function(data, model) {
    Y <- data$Y
    Xstar <- data$Xstar
    Lmat <- model.matrix(model, data = data)
    sw <- data$sw
    delta <- function(beta1) {
      Xstar + sigma_me*beta1*Y
    }
    H <- function(x) {
      1 / (1 + exp(-x))
    }
    condexp <- function(beta0, beta1) {
      H(beta0 + beta1*delta(beta1) - (beta1^2)*sigma_me/2)
    }
    function(theta) {
      p  <- length(theta)
      p1 <- length(coef(model))
      rho <- Lmat %*% theta[1:p1]

      score_eqns <- apply(Lmat, 2, function(x) sum((Xstar - rho) * x))

      c(score_eqns,
        Xstar - theta[p-2],
        sw*(Y - condexp(theta[p-1], theta[p])),
        sw*(Y - condexp(theta[p-1], theta[p]))*
          delta(theta[p])
      )
    }
  }
  
  failed <- TRUE
  j <- 1
  
  while(failed == TRUE & j < 2) {
    
    failed <- FALSE
    startvec <- c(coef(denom_mod), mean(Xstar), coef(wmod)[1:2])*(j == 1) +
      c(coef(denom_mod), mean(Xstar), coef(mod)[1:2])*(j == 2) +
      c(coef(denom_mod), mean(Xstar), rnorm(2, 0, j/5))*(j > 2) 
    results_ipw <- tryCatch(m_estimate(estFUN = eefun_ipw, data = data,
                                       outer_args = list(denom_mod),
                                       root_control = setup_root_control(start = startvec)),
                            error = function(e) { failed <<- TRUE})
    if (failed == FALSE) {
      if (abs(coef(results_ipw)[2]) > 2) { failed <- TRUE }
    }
    j <- j + 1
    
  }
  
  bias_ipw_psor <- coef(results_ipw)[6] - true_effect
  se_ipw_psor <- sqrt(vcov(results_ipw)[6, 6])
  coverage_ipw_psor<- 
    1*(coef(results_ipw)[6] - 1.96*se_ipw_psor < true_effect &
         coef(results_ipw)[6] + 1.96*se_ipw_psor > true_effect)
  
  # Model 3: AIPW
  # Include denominator density in g-formula CSME outcome model
  eefun_csme_aipw <- function(data, model) {
    Y <- data$Y
    Xstar <- data$Xstar
    L1 <- data$L1
    L2 <- data$L2
    Lmat <- model.matrix(model, data = data)
    sw <- data$sw
    delta <- function(beta1, beta4, beta5) {
      Xstar + sigma_me*(beta1 + beta4*L1 + beta5*L2)*Y
    }
    H <- function(x) {
      1 / (1 + exp(-x))
    }
    condexp <- function(beta0, beta1, beta2, beta3,
                        beta4, beta5) {
      H(beta0 + (beta1 + beta4*L1 + beta5*L2)*
          delta(beta1, beta4, beta5) + beta2*L1 + beta3*L2 -
          ((beta1 + beta4*L1 + beta5*L2)^{2})*sigma_me / 2)
    }
    function(theta) {
      p  <- length(theta)
      p1 <- length(coef(model))
      rho <- Lmat %*% theta[1:p1]

      score_eqns <- apply(Lmat, 2, function(x) sum((Xstar - rho) * x))

      c(score_eqns,
        Xstar - theta[p-9],
        sw*(Y - condexp(theta[p-8], theta[p-7], theta[p-6], theta[p-5],
                        theta[p-4], theta[p-3])),
        sw*(Y - condexp(theta[p-8], theta[p-7], theta[p-6], theta[p-5],
                        theta[p-4], theta[p-3]))*L1,
        sw*(Y - condexp(theta[p-8], theta[p-7], theta[p-6], theta[p-5],
                        theta[p-4], theta[p-3]))*L2,
        sw*(Y - condexp(theta[p-8], theta[p-7], theta[p-6], theta[p-5],
                        theta[p-4], theta[p-3]))*
          delta(theta[p-7], theta[p-4], theta[p-3]),
        sw*(Y - condexp(theta[p-8], theta[p-7], theta[p-6], theta[p-5],
                        theta[p-4], theta[p-3]))*
          L1*delta(theta[p-7], theta[p-4], theta[p-3]),
        sw*(Y - condexp(theta[p-8], theta[p-7], theta[p-6], theta[p-5],
                        theta[p-4], theta[p-3]))*
          L2*delta(theta[p-7], theta[p-4], theta[p-3]),
        L1 - theta[p-2],
        L2 - theta[p-1],
        theta[p-7] + theta[p-4]*theta[p-2] + theta[p-3]*theta[p-1] - theta[p]
      )
    }
  }
  
  results_csme_aipw <- m_estimate(estFUN = eefun_csme_aipw, data = data,
                                  compute_roots = TRUE,
                                  outer_args = list(denom_mod),
                                  root_control = 
                                    setup_root_control(start = c(coef(denom_mod),
                                                                 mean(Xstar),
                                                                 coef(mod),
                                                                 mean(L1),
                                                                 mean(L2),
                                                                 0.7)))
  
  bias_aipw_psor <- coef(results_csme_aipw)[13] - true_effect
  se_aipw_psor <- sqrt(vcov(results_csme_aipw)[13, 13])
  coverage_aipw_psor <- 
    1*(coef(results_csme_aipw)[13] - 1.96*se_aipw_psor < true_effect &
         coef(results_csme_aipw)[13] + 1.96*se_aipw_psor > true_effect)
  
  
  return(c(bias_gform_ps, bias_ipw_ps, bias_aipw_ps,
           bias_gform_or, bias_ipw_or, bias_aipw_or,
           bias_gform_psor, bias_ipw_psor, bias_aipw_psor,
           se_gform_ps, se_ipw_ps, se_aipw_ps,
           se_gform_or, se_ipw_or, se_aipw_or,
           se_gform_psor, se_ipw_psor, se_aipw_psor,
           coverage_gform_ps, coverage_ipw_ps, coverage_aipw_ps,
           coverage_gform_or, coverage_ipw_or, coverage_aipw_or,
           coverage_gform_psor, coverage_ipw_psor, coverage_aipw_psor))
  
  
}


trials <- seq(1, nsims)
combos <- data.frame(trials = rep(trials, length(beta1_true)),
                     mes = rep(sigma_me, each = nsims))
i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")) + 1500
combo_i <- combos[(i-1500), ]

set.seed(i*1000)
sim <- with(combo_i, mapply(simulator, trials, mes))

# Output
outfile <- paste("./Results/results_scen2_", i, ".Rdata", sep = "")
save(sim, file = outfile)