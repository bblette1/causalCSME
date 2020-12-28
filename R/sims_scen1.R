rm(list = ls())
library(geex)
library(rootSolve)

nsims <- 500
n <- 800
beta1_true <- 0.3
beta3_true <- -0.4
beta5_true <- -0.05
beta6_true <- -0.15
beta7_true <- -0.2
L1_prob <- 0.5
L2_prob <- 0.35
true_effect1 <- beta1_true + beta6_true*L1_prob
true_effect2 <- beta3_true + beta7_true*L2_prob
true_effect3 <- beta5_true
sigma_me1 <- 0.09
sigma_me3 <- 0.25
sigma_me <- matrix(c(sigma_me1, 0, 0, sigma_me3), byrow = T, ncol = 2)


simulator <- function(trial, beta1_true) {
  
  L1 <- rbinom(n, 1, L1_prob)
  L2 <- rbinom(n, 1, L2_prob)
  X1 <- rnorm(n, 4 + L1, 1)
  X2 <- rnorm(n, 1.25 + 0.5*L2, 0.5)
  X3 <- rnorm(n, 2.5, 0.7)
  Y <- -0.1 + beta1_true*X1 + 0.2*L1 + beta3_true*X2 - 0.3*L2 +
       beta5_true*X3 + beta6_true*X1*L1 + beta7_true*X2*L2 + rnorm(n, 0, 0.5)
  X1star <- X1 + rnorm(n, 0, sqrt(sigma_me1))
  X3star <- X3 + rnorm(n, 0, sqrt(sigma_me3))
  data <- data.frame("Y" = Y, "X1star" = X1star, "X2" = X2, "X3star" = X3star,
                     "L1" = L1, "L2" = L2)
  
  # Model 1: Simple outcome regression
  mod <- lm(Y ~ X1star*L1 + X2*L2 + X3star)
  
  bias1_lm <- coef(mod)[[2]] - true_effect1
  se1_lm <- summary(mod)$coefficients[2, 2]
  coverage1_lm <- 1*(coef(mod)[[2]] - 1.96*se1_lm < true_effect1 &
                       coef(mod)[[2]] + 1.96*se1_lm > true_effect1)
  
  bias2_lm <- coef(mod)[[4]] - true_effect2
  se2_lm <- summary(mod)$coefficients[4, 2]
  coverage2_lm <- 1*(coef(mod)[[4]] - 1.96*se2_lm < true_effect2 &
                       coef(mod)[[4]] + 1.96*se2_lm > true_effect2)
  
  bias3_lm <- coef(mod)[[6]] - true_effect3
  se3_lm <- summary(mod)$coefficients[6, 2]
  coverage3_lm <- 1*(coef(mod)[[6]] - 1.96*se3_lm < true_effect3 &
                       coef(mod)[[6]] + 1.96*se3_lm > true_effect3)
  
  # Model 2: Conditional score ME estimator (CSME)
  eefun_csme <- function(data) {
    Y <- data$Y
    X1star <- data$X1star
    X2 <- data$X2
    X3star <- data$X3star
    L1 <- data$L1
    L2 <- data$L2
    delta1 <- function(beta1, beta6, sigma_ep) {
      X1star + (beta1 + beta6*L1)*sigma_me1*Y / sigma_ep
    }
    delta3 <- function(beta5, sigma_ep) {
      X3star + beta5*sigma_me3*Y / sigma_ep
    }
    condexp <- function(beta0, beta1, beta2, beta3, beta4, beta5, beta6, beta7,
                        sigma_ep) {
      (beta0 + beta2*L1 + beta3*X2 + beta4*L2 +
         (beta1 + beta6*L1)*delta1(beta1, beta6, sigma_ep) +
          beta5*delta3(beta5, sigma_ep) + beta7*X2*L2) /
        (1 + ((beta1 + beta6*L1)^2 * sigma_me1 + beta5^2 * sigma_me3) /
           sigma_ep)[[1]]
    }
    condvar <- function(beta1, beta5, beta6, sigma_ep) {
      sigma_ep /
        (1 + ((beta1 + beta6*L1)^2 * sigma_me1 + beta5^2 * sigma_me3) /
           sigma_ep)[[1]]
    }
    function(theta) {
      c((Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                     theta[6], theta[7], theta[8], theta[9])),
        (Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                     theta[6], theta[7], theta[8], theta[9]))*
          delta1(theta[2], theta[7], theta[9]),
        (Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                     theta[6], theta[7], theta[8], theta[9]))*L1,
        (Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                     theta[6], theta[7], theta[8], theta[9]))*X2,
        (Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                     theta[6], theta[7], theta[8], theta[9]))*L2,
        (Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                     theta[6], theta[7], theta[8], theta[9]))*
          delta3(theta[6], theta[9]),
        (Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                     theta[6], theta[7], theta[8], theta[9]))*L1*
          delta1(theta[2], theta[7], theta[9]),
        (Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                     theta[6], theta[7], theta[8], theta[9]))*L2*X2,
        theta[9] - theta[9]*
          (Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                       theta[6], theta[7], theta[8], theta[9]))^2 / 
          condvar(theta[2], theta[6], theta[7], theta[9])
      )
    }
  }
  
  results_csme <- m_estimate(estFUN = eefun_csme, data = data,
                             compute_roots = TRUE,
                             root_control = 
                               setup_root_control(start = c(coef(mod),
                                                            sigma(mod)^2)))
  
  bias1_csme <- coef(results_csme)[2] - true_effect1
  se1_csme <- sqrt(vcov(results_csme)[2, 2])
  coverage1_csme <- 1*(coef(results_csme)[2] - 1.96*se1_csme < true_effect1 &
                       coef(results_csme)[2] + 1.96*se1_csme > true_effect1)
  bias2_csme <- coef(results_csme)[4] - true_effect2
  se2_csme <- sqrt(vcov(results_csme)[4, 4])
  coverage2_csme <- 1*(coef(results_csme)[4] - 1.96*se2_csme < true_effect2 &
                       coef(results_csme)[4] + 1.96*se2_csme > true_effect2)
  bias3_csme <- coef(results_csme)[6] - true_effect3
  se3_csme <- sqrt(vcov(results_csme)[6, 6])
  coverage3_csme <- 1*(coef(results_csme)[6] - 1.96*se3_csme < true_effect3 &
                         coef(results_csme)[6] + 1.96*se3_csme > true_effect3)
  
  # Model 3: G-formula, no ME correction
  # Find point estimates outside of geex
  EY_X1_5 <- 
    predict(mod, newdata = data.frame(X1star = 5, L1 = 1, X2 = 1.5,
                                      L2 = 0.5, X3star = 2.5))*mean(L1) +
    predict(mod, newdata = data.frame(X1star = 5, L1 = 0, X2 = 1.5,
                                      L2 = 0.5, X3star = 2.5))*(1 - mean(L1))
  EY_X1_4 <- 
    predict(mod, newdata = data.frame(X1star = 4, L1 = 1, X2 = 1.5,
                                      L2 = 0.5, X3star = 2.5))*mean(L1) +
    predict(mod, newdata = data.frame(X1star = 4, L1 = 0, X2 = 1.5,
                                      L2 = 0.5, X3star = 2.5))*(1 - mean(L1))
  bias1_gform <- EY_X1_5 - EY_X1_4 - true_effect1
  
  EY_X2_2 <-
    predict(mod, newdata = data.frame(X1star = 4.5, L1 = 0.5, X2 = 2,
                                      L2 = 1, X3star = 2.5))*mean(L2) +
    predict(mod, newdata = data.frame(X1star = 4.5, L1 = 0.5, X2 = 2,
                                      L2 = 0, X3star = 2.5))*(1 - mean(L2))
  EY_X2_1 <- 
    predict(mod, newdata = data.frame(X1star = 4.5, L1 = 0.5, X2 = 1,
                                      L2 = 1, X3star = 2.5))*mean(L2) +
    predict(mod, newdata = data.frame(X1star = 4.5, L1 = 0.5, X2 = 1,
                                      L2 = 0, X3star = 2.5))*(1 - mean(L2))
  bias2_gform <- EY_X2_2 - EY_X2_1 - true_effect2
  
  
  EY_X3_3 <-
    predict(mod, newdata = data.frame(X1star = 4.5, L1 = 0.5, X2 = 1.5,
                                      L2 = 0.5, X3star = 3))
  EY_X3_2 <- 
    predict(mod, newdata = data.frame(X1star = 4.5, L1 = 0.5, X2 = 1.5,
                                      L2 = 0.5, X3star = 2))
  bias3_gform <- EY_X3_3 - EY_X3_2 - true_effect3
  
  # Use geex for variance estimation
  eefun_gform <- function(data) {
    Y <- data$Y
    X1star <- data$X1star
    X2 <- data$X2
    X3star <- data$X3star
    L1 <- data$L1
    L2 <- data$L2
    mm <- model.matrix(Y ~ X1star*L1 + X2*L2 + X3star, drop = FALSE)
    k <- dim(mm)[2]
    function(theta) {
      mu <- mm %*% theta[3:(k + 2)]
      c(L1 - theta[1],
        L2 - theta[2],
        (Y - (theta[3] + X1star*theta[4] + L1*theta[5] + X2*theta[6] + 
              L2*theta[7] + X3star*theta[8] + X1star*L1*theta[9] +
              X2*L2*theta[10])) *
          c(1, X1star, L1, X2, L2, X3star, X1star*L1, X2*L2),
        t(c(1, 5, 1, 1.5, 0.5, 2.5, 5, 0.75)) %*% theta[3:(k + 2)]*theta[1] +
          t(c(1, 5, 0, 1.5, 0.5, 2.5, 0, 0.75)) %*% theta[3:(k + 2)]*(1 - theta[1]) -
          t(c(1, 4, 1, 1.5, 0.5, 2.5, 4, 0.75)) %*% theta[3:(k + 2)]*theta[1] -
          t(c(1, 4, 0, 1.5, 0.5, 2.5, 0, 0.75)) %*% theta[3:(k + 2)]*(1 - theta[1]) - theta[k + 3],
        t(c(1, 4.5, 0.5, 2, 1, 2.5, 2.25, 2)) %*% theta[3:(k + 2)]*theta[2] +
          t(c(1, 4.5, 0.5, 2, 0, 2.5, 2.25, 0)) %*% theta[3:(k + 2)]*(1 - theta[2]) -
          t(c(1, 4.5, 0.5, 1, 1, 2.5, 2.25, 1)) %*% theta[3:(k + 2)]*theta[2] -
          t(c(1, 4.5, 0.5, 1, 0, 2.5, 2.25, 0)) %*% theta[3:(k + 2)]*(1 - theta[2]) - theta[k + 4],
        t(c(1, 4.5, 0.5, 1.5, 0.5, 3, 2.25, 0.75)) %*% theta[3:(k + 2)] -
          t(c(1, 4.5, 0.5, 1.5, 0.5, 2, 2.25, 0.75)) %*% theta[3:(k + 2)] - theta[k + 5]
      )
    }
  }
  
  results_gform <- m_estimate(estFUN = eefun_gform, data = data,
                              compute_roots = FALSE,
                              roots = c(mean(L1), mean(L2), coef(mod),
                                        EY_X1_5 - EY_X1_4, EY_X2_2 - EY_X2_1,
                                        EY_X3_3 - EY_X3_2))
  
  se1_gform <- sqrt(vcov(results_gform)[11, 11])
  coverage1_gform <- 1*(EY_X1_5 - EY_X1_4 - 1.96*se1_gform < true_effect1 &
                        EY_X1_5 - EY_X1_4 + 1.96*se1_gform > true_effect1)
  
  se2_gform <- sqrt(vcov(results_gform)[12, 12])
  coverage2_gform <- 1*(EY_X2_2 - EY_X2_1 - 1.96*se2_gform < true_effect2 &
                          EY_X2_2 - EY_X2_1 + 1.96*se2_gform > true_effect2)
  
  se3_gform <- sqrt(vcov(results_gform)[13, 13])
  coverage3_gform <- 1*(EY_X3_3 - EY_X3_2 - 1.96*se3_gform < true_effect3 &
                          EY_X3_3 - EY_X3_2 + 1.96*se3_gform > true_effect3)
  
  # Model 4: G-formula, using CSME for outcome prediction
  # Get point estimates outside of geex
  EY_X1_5_EE <- t(coef(results_csme)[1:8]) %*% c(1, 5, 1, 1.5, 0.5, 2.5, 5, 0.75)*mean(L1) +
    t(coef(results_csme)[1:8]) %*% c(1, 5, 0, 1.5, 0.5, 2.5, 0, 0.75)*(1 - mean(L1))
  EY_X1_4_EE <- t(coef(results_csme)[1:8]) %*% c(1, 4, 1, 1.5, 0.5, 2.5, 4, 0.75)*mean(L1) +
    t(coef(results_csme)[1:8]) %*% c(1, 4, 0, 1.5, 0.5, 2.5, 0, 0.75)*(1 - mean(L1))
  bias1_gform_csme <- EY_X1_5_EE - EY_X1_4_EE - true_effect1
  
  EY_X2_2_EE <- t(coef(results_csme)[1:8]) %*% c(1, 4.5, 0.5, 2, 1, 2.5, 2.25, 2)*mean(L2) +
    t(coef(results_csme)[1:8]) %*% c(1, 4.5, 0.5, 2, 0, 2.5, 2.25, 0)*(1 - mean(L2))
  EY_X2_1_EE <- t(coef(results_csme)[1:8]) %*% c(1, 4.5, 0.5, 1, 1, 2.5, 2.25, 1)*mean(L2) +
    t(coef(results_csme)[1:8]) %*% c(1, 4.5, 0.5, 1, 0, 2.5, 2.25, 0)*(1 - mean(L2))
  bias2_gform_csme <- EY_X2_2_EE - EY_X2_1_EE - true_effect2
  
  EY_X3_3_EE <- t(coef(results_csme)[1:8]) %*% c(1, 4.5, 0.5, 1.5, 0.5, 3, 2.25, 0.75)
  EY_X3_2_EE <- t(coef(results_csme)[1:8]) %*% c(1, 4.5, 0.5, 1.5, 0.5, 2, 2.25, 0.75)
  bias3_gform_csme <- EY_X3_3_EE - EY_X3_2_EE - true_effect3
  
  # Estimate variance using geex
  eefun_csme_gform <- function(data) {
    Y <- data$Y
    X1star <- data$X1star
    X2 <- data$X2
    X3star <- data$X3star
    L1 <- data$L1
    L2 <- data$L2
    delta1 <- function(beta1, beta6, sigma_ep) {
      X1star + (beta1 + beta6*L1)*sigma_me1*Y / sigma_ep
    }
    delta3 <- function(beta5, sigma_ep) {
      X3star + beta5*sigma_me3*Y / sigma_ep
    }
    condexp <- function(beta0, beta1, beta2, beta3, beta4, beta5, beta6, beta7,
                        sigma_ep) {
      (beta0 + beta2*L1 + beta3*X2 + beta4*L2 +
         (beta1 + beta6*L1)*delta1(beta1, beta6, sigma_ep) +
         beta5*delta3(beta5, sigma_ep) + beta7*X2*L2) /
        (1 + ((beta1 + beta6*L1)^2 * sigma_me1 + beta5^2 * sigma_me3) /
           sigma_ep)[[1]]
    }
    condvar <- function(beta1, beta5, beta6, sigma_ep) {
      sigma_ep /
        (1 + ((beta1 + beta6*L1)^2 * sigma_me1 + beta5^2 * sigma_me3) /
           sigma_ep[[1]])
    }
    function(theta) {
      c((Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                     theta[6], theta[7], theta[8], theta[9])),
        (Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                     theta[6], theta[7], theta[8], theta[9]))*
          delta1(theta[2], theta[7], theta[9]),
        (Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                     theta[6], theta[7], theta[8], theta[9]))*L1,
        (Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                     theta[6], theta[7], theta[8], theta[9]))*X2,
        (Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                     theta[6], theta[7], theta[8], theta[9]))*L2,
        (Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                     theta[6], theta[7], theta[8], theta[9]))*
          delta3(theta[6], theta[9]),
        (Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                     theta[6], theta[7], theta[8], theta[9]))*L1*
          delta1(theta[2], theta[7], theta[9]),
        (Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                     theta[6], theta[7], theta[8], theta[9]))*L2*X2,
        theta[9] - theta[9]*
          (Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                       theta[6], theta[7], theta[8], theta[9]))^2 / 
          condvar(theta[2], theta[6], theta[7], theta[9]),
        L1 - theta[10],
        L2 - theta[11],
        t(c(1, 5, 1, 1.5, 0.5, 2.5, 5, 0.75)) %*% theta[1:8]*theta[10] +
          t(c(1, 5, 0, 1.5, 0.5, 2.5, 0, 0.75)) %*% theta[1:8]*(1 - theta[10]) -
          t(c(1, 4, 1, 1.5, 0.5, 2.5, 4, 0.75)) %*% theta[1:8]*theta[10] -
          t(c(1, 4, 0, 1.5, 0.5, 2.5, 0, 0.75)) %*% theta[1:8]*(1 - theta[10]) - theta[12],
        t(c(1, 4.5, 0.5, 2, 1, 2.5, 2.25, 2)) %*% theta[1:8]*theta[11] +
          t(c(1, 4.5, 0.5, 2, 0, 2.5, 2.25, 0)) %*% theta[1:8]*(1 - theta[11]) -
          t(c(1, 4.5, 0.5, 1, 1, 2.5, 2.25, 1)) %*% theta[1:8]*theta[11] -
          t(c(1, 4.5, 0.5, 1, 0, 2.5, 2.25, 0)) %*% theta[1:8]*(1 - theta[11]) - theta[13],
        t(c(1, 4.5, 0.5, 1.5, 0.5, 3, 2.25, 0.75)) %*% theta[1:8] -
          t(c(1, 4.5, 0.5, 1.5, 0.5, 2, 2.25, 0.75)) %*% theta[1:8] - theta[14]
      )
    }
  }
  
  results_csme_gform <- m_estimate(estFUN = eefun_csme_gform, data = data,
                                   compute_roots = FALSE,
                                   roots = c(coef(results_csme), mean(L1),
                                             mean(L2), EY_X1_5_EE - EY_X1_4_EE,
                                             EY_X2_2_EE - EY_X2_1_EE,
                                             EY_X3_3_EE - EY_X3_2_EE))
  
  se1_gform_csme <- sqrt(vcov(results_csme_gform)[12, 12])
  coverage1_gform_csme <- 
    1*(EY_X1_5_EE - EY_X1_4_EE - 1.96*se1_gform_csme < true_effect1 &
       EY_X1_5_EE - EY_X1_4_EE + 1.96*se1_gform_csme > true_effect1)
  
  se2_gform_csme <- sqrt(vcov(results_csme_gform)[13, 13])
  coverage2_gform_csme <- 
    1*(EY_X2_2_EE - EY_X2_1_EE - 1.96*se2_gform_csme < true_effect2 &
       EY_X2_2_EE - EY_X2_1_EE + 1.96*se2_gform_csme > true_effect2)
  
  se3_gform_csme <- sqrt(vcov(results_csme_gform)[14, 14])
  coverage3_gform_csme <- 
    1*(EY_X3_3_EE - EY_X3_2_EE - 1.96*se3_gform_csme < true_effect3 &
         EY_X3_3_EE - EY_X3_2_EE + 1.96*se3_gform_csme > true_effect3)
  
  # Model 5: IPW estimator (no ME correction)
  # Estimate weights
  denom_mod1 <- lm(X1star ~ L1)
  p_denom1 <- predict(denom_mod1, type='response')
  dens_denom1 <- dnorm(X1star, p_denom1, summary(denom_mod1)$sigma)
  num_mod1 <- lm(X1star ~ 1)
  p_num1 <- predict(num_mod1, type='response')
  dens_num1 <- dnorm(X1star, p_num1, summary(denom_mod1)$sigma)
  data$sw1 <- dens_num1 / dens_denom1
  
  denom_mod2 <- lm(X2 ~ L2)
  p_denom2 <- predict(denom_mod2, type='response')
  dens_denom2 <- dnorm(X2, p_denom2, summary(denom_mod2)$sigma)
  num_mod2 <- lm(X2 ~ 1)
  p_num2 <- predict(num_mod2, type='response')
  dens_num2 <- dnorm(X2, p_num2, summary(denom_mod2)$sigma)
  data$sw2 <- dens_num2 / dens_denom2
  
  data$sw <- data$sw1*data$sw2
  
  # Fit weighted regression
  wmod <- lm(Y ~ X1star + X2 + X3star, weights = data$sw)
  bias1_ipw <- coef(wmod)[2] - true_effect1
  bias2_ipw <- coef(wmod)[3] - true_effect2
  bias3_ipw <- coef(wmod)[4] - true_effect3
  
  eefun_wr <- function(data, model1, model2, model3, model4) {
    Y <- data$Y
    X3star <- data$X3star
    X1star <- model.response(model.frame(model1, data = data))
    L1 <- model.matrix(model1, data = data)
    V1 <- model.matrix(model2, data = data)
    X2 <- model.response(model.frame(model3, data = data))
    L2 <- model.matrix(model3, data = data)
    V2 <- model.matrix(model4, data = data)
    #sw <- data$sw
    n <- dim(data)[1]
    
    function(theta) {
      p  <- length(theta)
      p1 <- length(coef(model1))
      p2 <- length(coef(model2))
      p3 <- length(coef(model3))
      p4 <- length(coef(model4))
      rho1 <- L1 %*% theta[1:p1]
      rho2 <- V1 %*% theta[(p1+1):(p1+p2)]
      rho3 <- L2 %*% theta[(p1+p2+1):(p1+p2+p3)]
      rho4 <- V2 %*% theta[(p1+p2+p3+1):(p1+p2+p3+p4)]
      
      dens_denom1 <- dnorm(X1star, rho1, sqrt(theta[p-6]))
      dens_num1 <- dnorm(X1star, rho2, sqrt(theta[p-6]))
      dens_denom2 <- dnorm(X2, rho3, sqrt(theta[p-5]))
      dens_num2 <- dnorm(X2, rho4, sqrt(theta[p-5]))
      sw1 <- dens_num1 / dens_denom1
      sw2 <- dens_num2 / dens_denom2
      sw <- sw1*sw2
      
      score_eqns1 <- apply(L1, 2, function(x) sum((X1star - rho1) * x))
      score_eqns2 <- apply(V1, 2, function(x) sum((X1star - rho2) * x))
      score_eqns3 <- apply(L2, 2, function(x) sum((X2 - rho3) * x))
      score_eqns4 <- apply(V2, 2, function(x) sum((X2 - rho4) * x))
      
      c((X1star - rho1)*L1[,1],
        (X1star - rho1)*L1[,2],
        (X1star - rho2)*V1[,1],
        (X2 - rho3)*L2[,1],
        (X2 - rho3)*L2[,2],
        (X2 - rho4)*V2[,1],
        (n-p1)/n * theta[p-6] - (X1star - L1 %*% theta[1:p1])^2,
        (n-p3)/n * theta[p-5] - (X2 - L2 %*% theta[(p1+p2+1):(p1+p2+p3)])^2,
        sw*(Y - (theta[p-4] + theta[p-3]*X1star + theta[p-2]*X2 + theta[p-1]*X3star)),
        sw*(Y - (theta[p-4] + theta[p-3]*X1star + theta[p-2]*X2 + theta[p-1]*X3star))*X1star,
        sw*(Y - (theta[p-4] + theta[p-3]*X1star + theta[p-2]*X2 + theta[p-1]*X3star))*X2,
        sw*(Y - (theta[p-4] + theta[p-3]*X1star + theta[p-2]*X2 + theta[p-1]*X3star))*X3star,
        sw*(theta[p] - (Y - (theta[p-4] + theta[p-3]*X1star + theta[p-2]*X2 + theta[p-1]*X3star))^2)
      )
    }
  }
    
  results_wr <- m_estimate(estFUN = eefun_wr, data = data,
                           outer_args = list(denom_mod1, num_mod1,
                                             denom_mod2, num_mod2),
                           compute_roots = FALSE,
                           roots = c(coef(denom_mod1), coef(num_mod1),
                                     coef(denom_mod2), coef(num_mod2),
                                     sigma(denom_mod1)^2, sigma(denom_mod2)^2,
                                     coef(wmod), sigma(wmod)^2))
  
  se1_ipw <- sqrt(vcov(results_wr)[10, 10])
  coverage1_ipw <- 1*(coef(wmod)[2] - 1.96*se1_ipw < true_effect1 &
                      coef(wmod)[2] + 1.96*se1_ipw > true_effect1)
  
  se2_ipw <- sqrt(vcov(results_wr)[11, 11])
  coverage2_ipw <- 1*(coef(wmod)[3] - 1.96*se2_ipw < true_effect2 &
                      coef(wmod)[3] + 1.96*se2_ipw > true_effect2)
  
  se3_ipw <- sqrt(vcov(results_wr)[12, 12])
  coverage3_ipw <- 1*(coef(wmod)[4] - 1.96*se3_ipw < true_effect3 &
                        coef(wmod)[4] + 1.96*se3_ipw > true_effect3)
  
  # Model 6: IPW-CSME
  # Get point estimates and variance using geex
  # We get point estimates using a geex with only some parameters
  # Then variance using geex with the full set of parameters
  # We don't get point estimates with the second geex because of
  # high level of divergence
  eefun_ipw_csme1 <- function(data) {
    Y <- data$Y
    X1star <- data$X1star
    X2 <- data$X2
    X3star <- data$X3star
    sw <- data$sw
    delta1 <- function(beta1, sigma_ep) {
      X1star + beta1*sigma_me1*Y / sigma_ep
    }
    delta3 <- function(beta3, sigma_ep) {
      X3star + beta3*sigma_me3*Y / sigma_ep
    }
    condexp <- function(beta0, beta1, beta2, beta3, sigma_ep) {
      (beta0 + beta1*delta1(beta1, sigma_ep) + beta2*X2 +
         beta3*delta3(beta3, sigma_ep)) /
        (1 + (beta1^2 * sigma_me1 + beta3^2 * sigma_me3) /
           sigma_ep)
    }
    condvar <- function(beta1, beta3, sigma_ep) {
      sigma_ep /
        (1 + (beta1^2 * sigma_me1 + beta3^2 * sigma_me3) /
           sigma_ep)
    }
    function(theta) {
      c(sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5])),
        sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5]))*
          delta1(theta[2], theta[5]),
        sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5]))*
          X2,
        sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5]))*
          delta3(theta[4], theta[5]),
        sw*(theta[5] - theta[5]*
              (Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5]))^2 / 
              condvar(theta[2], theta[4], theta[5]))
      )
    }
  }
  
  failed <- TRUE
  j <- 1
  
  while(failed == TRUE & j < 6) {
    
    failed <- FALSE
    startvec <- c(coef(wmod)[1:4], sigma(wmod))*(j == 1) +
      c(coef(mod)[1:4], sigma(mod))*(j == 2) +
      c(rnorm(4, 0, j/5), sigma(wmod))*(j > 2) 
    results_ipw_csme1 <- tryCatch(m_estimate(estFUN = eefun_ipw_csme1, data = data,
                                        root_control = setup_root_control(start = startvec)),
                             error = function(e) { failed <<- TRUE})
    if (failed == FALSE) {
      if (abs(coef(results_ipw_csme1)[2]) > 2) { failed <- TRUE }
    }
    j <- j + 1
    
  }
  
  if(failed == TRUE) { next }
  
  bias1_ipw_csme <- coef(results_ipw_csme1)[2] - true_effect1
  bias2_ipw_csme <- coef(results_ipw_csme1)[3] - true_effect2
  bias3_ipw_csme <- coef(results_ipw_csme1)[4] - true_effect3
  
  eefun_ipw_csme2 <- function(data, model1, model2, model3, model4) {
    Y <- data$Y
    X3star <- data$X3star
    X1star <- model.response(model.frame(model1, data = data))
    L1 <- model.matrix(model1, data = data)
    V1 <- model.matrix(model2, data = data)
    X2 <- model.response(model.frame(model3, data = data))
    L2 <- model.matrix(model3, data = data)
    V2 <- model.matrix(model4, data = data)
    #sw <- data$sw
    n <- dim(data)[1]
    delta1 <- function(beta1, sigma_ep) {
      X1star + beta1*sigma_me1*Y / sigma_ep
    }
    delta3 <- function(beta3, sigma_ep) {
      X3star + beta3*sigma_me3*Y / sigma_ep
    }
    condexp <- function(beta0, beta1, beta2, beta3, sigma_ep) {
      (beta0 + beta1*delta1(beta1, sigma_ep) + beta2*X2 +
         beta3*delta3(beta3, sigma_ep)) /
        (1 + (beta1^2 * sigma_me1 + beta3^2 * sigma_me3) /
           sigma_ep)
    }
    condvar <- function(beta1, beta3, sigma_ep) {
      sigma_ep /
        (1 + (beta1^2 * sigma_me1 + beta3^2 * sigma_me3) /
           sigma_ep)
    }
    function(theta) {
      p  <- length(theta)
      p1 <- length(coef(model1))
      p2 <- length(coef(model2))
      p3 <- length(coef(model3))
      p4 <- length(coef(model4))
      rho1 <- L1 %*% theta[1:p1]
      rho2 <- V1 %*% theta[(p1+1):(p1+p2)]
      rho3 <- L2 %*% theta[(p1+p2+1):(p1+p2+p3)]
      rho4 <- V2 %*% theta[(p1+p2+p3+1):(p1+p2+p3+p4)]
      
      dens_denom1 <- dnorm(X1star, rho1, sqrt(theta[p-6]))
      dens_num1 <- dnorm(X1star, rho2, sqrt(theta[p-6]))
      dens_denom2 <- dnorm(X2, rho3, sqrt(theta[p-5]))
      dens_num2 <- dnorm(X2, rho4, sqrt(theta[p-5]))
      sw1 <- dens_num1 / dens_denom1
      sw2 <- dens_num2 / dens_denom2
      sw <- sw1*sw2
      
      score_eqns1 <- apply(L1, 2, function(x) sum((X1star - rho1) * x))
      score_eqns2 <- apply(V1, 2, function(x) sum((X1star - rho2) * x))
      score_eqns3 <- apply(L2, 2, function(x) sum((X2 - rho3) * x))
      score_eqns4 <- apply(V2, 2, function(x) sum((X2 - rho4) * x))
      
      c(score_eqns1,
        score_eqns2,
        score_eqns3,
        score_eqns4,
        (n-p1)/n * theta[p-6] - (X1star - L1 %*% theta[1:p1])^2,
        (n-p3)/n * theta[p-5] - (X2 - L2 %*% theta[(p1+p2+1):(p1+p2+p3)])^2,
        sw*(Y - condexp(theta[p-4], theta[p-3], theta[p-2], theta[p-1], theta[p])),
        sw*(Y - condexp(theta[p-4], theta[p-3], theta[p-2], theta[p-1], theta[p]))*
          delta1(theta[p-3], theta[p]),
        sw*(Y - condexp(theta[p-4], theta[p-3], theta[p-2], theta[p-1], theta[p]))*X2,
        sw*(Y - condexp(theta[p-4], theta[p-3], theta[p-2], theta[p-1], theta[p]))*
          delta3(theta[p-1], theta[p]),
        sw*(theta[p] - theta[p]*
              (Y - condexp(theta[p-4], theta[p-3], theta[p-2], theta[p-1], theta[p]))^2 / 
              condvar(theta[p-3], theta[p-1], theta[p]))
      )
    }
  }
  
  results_ipw_csme2 <- m_estimate(estFUN = eefun_ipw_csme2, data = data,
                             outer_args = list(denom_mod1, num_mod1,
                                               denom_mod2, num_mod2),
                             compute_roots = FALSE,
                             roots = c(coef(denom_mod1), coef(num_mod1),
                                       coef(denom_mod2), coef(num_mod2),
                                       sigma(denom_mod1)^2, sigma(denom_mod2)^2,
                                       coef(results_ipw_csme1)))
  
  se1_ipw_csme <- sqrt(vcov(results_ipw_csme2)[10, 10])
  coverage1_ipw_csme <- 
    1*(coef(results_ipw_csme1)[2] - 1.96*se1_ipw_csme < true_effect1 &
       coef(results_ipw_csme1)[2] + 1.96*se1_ipw_csme > true_effect1)
  
  se2_ipw_csme <- sqrt(vcov(results_ipw_csme2)[11, 11])
  coverage2_ipw_csme <- 
    1*(coef(results_ipw_csme1)[3] - 1.96*se2_ipw_csme < true_effect2 &
       coef(results_ipw_csme1)[3] + 1.96*se2_ipw_csme > true_effect2)
  
  se3_ipw_csme <- sqrt(vcov(results_ipw_csme2)[12, 12])
  coverage3_ipw_csme <- 
    1*(coef(results_ipw_csme1)[4] - 1.96*se3_ipw_csme < true_effect3 &
         coef(results_ipw_csme1)[4] + 1.96*se3_ipw_csme > true_effect3)
  
  return(c(bias1_lm, bias1_csme, bias1_gform, bias1_gform_csme,
           bias1_ipw, bias1_ipw_csme,
           bias2_lm, bias2_csme, bias2_gform, bias2_gform_csme,
           bias2_ipw, bias2_ipw_csme,
           bias3_lm, bias3_csme, bias3_gform, bias3_gform_csme,
           bias3_ipw, bias3_ipw_csme,
           se1_lm, se1_csme, se1_gform, se1_gform_csme, se1_ipw, se1_ipw_csme,
           se2_lm, se2_csme, se2_gform, se2_gform_csme, se2_ipw, se2_ipw_csme,
           se3_lm, se3_csme, se3_gform, se3_gform_csme, se3_ipw, se3_ipw_csme,
           coverage1_lm, coverage1_csme, coverage1_gform, coverage1_gform_csme,
           coverage1_ipw, coverage1_ipw_csme,
           coverage2_lm, coverage2_csme, coverage2_gform, coverage2_gform_csme,
           coverage2_ipw, coverage2_ipw_csme,
           coverage3_lm, coverage3_csme, coverage3_gform, coverage3_gform_csme,
           coverage3_ipw, coverage3_ipw_csme))
  
}


trials <- seq(1, nsims)
combos <- data.frame(trials = rep(trials, length(beta1_true)),
                     betas = rep(beta1_true, each = nsims))
i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
combo_i <- combos[(i), ]

set.seed(i*1000)
sim <- with(combo_i, mapply(simulator, trials, betas))

# Output
outfile <- paste("./Results/results_scen1_", i, ".Rdata", sep = "")
save(sim, file = outfile)
