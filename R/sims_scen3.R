rm(list = ls())
library(geex)
library(rootSolve)

# Set parameter values
nsims <- 1000
n <- 1500
beta1_true <- 0.7
beta3_true <- -0.7
beta5_true <- 0.4
L1_prob <- 0.5
L2_mean <- 1
true_effect <- beta1_true + beta3_true*L1_prob + beta5_true*L2_mean
sigma_me <- 0.16

# Simulator function to be sent to computing cluster
simulator <- function(trial, sigma_me, n) {

  # Generate data
  L1 <- rbinom(n, 1, L1_prob)
  L2 <- rnorm(n, L2_mean, 0.5)
  A <- rnorm(n, 2 + 0.9*L1 - 0.6*L2, 1.1)
  Y_mean <- 1.5 + beta1_true*A + 0.9*L1 + beta3_true*A*L1 - 0.6*L2 +
            beta5_true*A*L2
  Y <- rnorm(n, Y_mean, 0.4)
  Astar <- A + rnorm(n, 0, sqrt(sigma_me))
  data <- data.frame("Y" = Y, "Astar" = Astar, "L1" = L1, "L2" = L2)

  # Initial guess for causal effect (good but variable intuition)
  guess <- true_effect + rnorm(1, 0, 0.1)

  ######################################################################
  # Outcome regression wrong, propensity model correct

  # Fit regression to use for starting values
  mod <- lm(Y ~ Astar*L2, data = data)

  # Model 1: G-formula-CSME
  eefun_csme_gform <- function(data) {
    Y <- data$Y
    Astar <- data$Astar
    L2 <- data$L2
    delta <- function(beta1, beta3, sigma_ep) {
      Astar + sigma_me*(beta1 + beta3*L2)*Y / sigma_ep
    }
    condexp <- function(beta0, beta1, beta2, beta3, sigma_ep) {
      (beta0 + (beta1 + beta3*L2)*delta(beta1, beta3, sigma_ep) +
       beta2*L2) /
        (1 + ((beta1 + beta3*L2)^2)*sigma_me / sigma_ep)[[1]]
    }
    condvar <- function(beta1, beta3, sigma_ep) {
      sigma_ep / (1 + ((beta1 + beta3*L2)^2)*sigma_me / sigma_ep)[[1]]
    }
    function(theta) {
      c((Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5])),
        (Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5]))*
          delta(theta[2], theta[4], theta[5]),
        (Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5]))*L2,
        (Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5]))*
          L2*delta(theta[2], theta[4], theta[5]),
        theta[5] - theta[5]*
          (Y-condexp(theta[1], theta[2], theta[3], theta[4], theta[5]))^2 /
          condvar(theta[2], theta[4], theta[5]),
        L2 - theta[6],
        theta[2] + theta[4]*theta[6] - theta[7]
      )
    }
  }

  # Solve the estimating equations in a try-catch system
  # Almost always works on first try, but gives extra opportunities for
  # new starting values and sets to NA if always divergent
  failed <- TRUE
  j <- 1

  while(failed == TRUE & j < 6) {

    failed <- FALSE
    startvec <- c(coef(mod), sigma(mod)^2, mean(L2), guess)*(j == 1) +
      c(coef(mod) + rnorm(4, 0, j/5), sigma(mod)^2, mean(L2), guess)*(j > 1)
    results_csme_gform <-
      tryCatch(m_estimate(estFUN = eefun_csme_gform, data = data,
                          root_control =
                            setup_root_control(start = startvec)),
               error = function(e) { failed <<- TRUE})
    if (failed == FALSE) {
      if (abs(coef(results_csme_gform)[7]) > 3) { failed <- TRUE }
    }
    j <- j + 1

  }

  bias_gform_ps <- coef(results_csme_gform)[7] - true_effect
  se_gform_ps <- sqrt(vcov(results_csme_gform)[7, 7])
  coverage_gform_ps <-
    1*(coef(results_csme_gform)[7] - 1.96*se_gform_ps < true_effect &
         coef(results_csme_gform)[7] + 1.96*se_gform_ps > true_effect)

  if(failed == TRUE) {
    bias_gform_ps <- NA
    se_gform_ps <- NA
    coverage_gform_ps <- NA
  }

  # Model 2: IPW-CSME
  # Estimate weights
  denom_mod <- lm(Astar ~ L1 + L2)
  p_denom <- predict(denom_mod, type='response')
  dens_denom <- dnorm(Astar, p_denom, summary(denom_mod)$sigma)
  num_mod <- lm(Astar ~ 1)
  p_num <- predict(num_mod, type='response')
  dens_num <- dnorm(Astar, p_num, summary(denom_mod)$sigma)
  data$sw <- dens_num / dens_denom

  # Fit weighted regression for starting values
  wmod <- lm(Y ~ Astar, weights = sw, data = data)

  # Get point estimates and variance using geex
  # First point estimates
  eefun_csme_ipw1 <- function(data) {
    Y <- data$Y
    Astar <- data$Astar
    sw <- data$sw
    delta <- function(beta1, sigma_ep) {
      Astar + sigma_me*beta1*Y / sigma_ep
    }
    condexp <- function(beta0, beta1, sigma_ep) {
      (beta0 + beta1*delta(beta1, sigma_ep)) /
        (1 + (beta1^2)*sigma_me / sigma_ep)
    }
    condvar <- function(beta1, sigma_ep) {
      sigma_ep /
        (1 + (beta1^2)*sigma_me / sigma_ep)
    }
    function(theta) {
      c(sw*(Y - condexp(theta[1], theta[2], theta[3])),
        sw*(Y - condexp(theta[1], theta[2], theta[3]))*
          delta(theta[2], theta[3]),
        sw*(theta[3] - theta[3]*
          (Y - condexp(theta[1], theta[2], theta[3]))^2 /
          condvar(theta[2], theta[3]))
      )
    }
  }

  results_csme_ipw <-
    m_estimate(estFUN = eefun_csme_ipw1, data = data,
               root_control = setup_root_control(start = c(coef(wmod),
                                                           sigma(wmod)^2)))
  bias_ipw_ps <- coef(results_csme_ipw)[2] - true_effect

  # Then get variance estimates accounting for weight estimation
  eefun_ipw2 <- function(data, model1, model2) {
    Y <- data$Y
    Astar <- model.response(model.frame(model1, data = data))
    Lmat <- model.matrix(model1, data = data)
    Vmat <- model.matrix(model2, data = data)
    n <- dim(data)[1]
    delta <- function(beta1, sigma_ep) {
      Astar + sigma_me*beta1*Y / sigma_ep
    }
    condexp <- function(beta0, beta1, sigma_ep) {
      (beta0 + beta1*delta(beta1, sigma_ep)) /
        (1 + (beta1^2)*sigma_me / sigma_ep)
    }
    condvar <- function(beta1, sigma_ep) {
      sigma_ep /
        (1 + (beta1^2)*sigma_me / sigma_ep)
    }
    function(theta) {
      p  <- length(theta)
      p1 <- length(coef(model1))
      p2 <- length(coef(model2))
      rho1 <- Lmat %*% theta[1:p1]
      rho2 <- Vmat %*% theta[(p1+1):(p1+p2)]

      dens_denom <- dnorm(Astar, rho1, sqrt(theta[p-2]))
      dens_num <- dnorm(Astar, rho2, sqrt(theta[p-2]))
      sw <- dens_num / dens_denom

      score_eqns1 <- apply(Lmat, 2, function(x) sum((Astar - rho1) * x))
      score_eqns2 <- apply(Vmat, 2, function(x) sum((Astar - rho2) * x))

      c(score_eqns1,
        score_eqns2,
        (n-p1)/n * theta[p-3] - (Astar - Lmat %*% theta[1:p1])^2,
        sw*(Y - condexp(theta[p-2], theta[p-1], theta[p])),
        sw*(Y - condexp(theta[p-2], theta[p-1], theta[p]))*
          delta(theta[p-1], theta[p]),
        sw*(theta[p] - theta[p]*
              (Y - condexp(theta[p-2], theta[p-1], theta[p]))^2 /
              condvar(theta[p-1], theta[p]))
      )
    }
  }

  results_ipw2 <-
    m_estimate(estFUN = eefun_ipw2, data = data,
               outer_args = list(model1 = denom_mod, model2 = num_mod),
               compute_roots = FALSE,
               roots = c(coef(denom_mod), coef(num_mod),
                         sigma(denom_mod)^2, coef(results_csme_ipw)))

  #se_ipw_ps <- sqrt(vcov(results_ipw2)[7, 7])
  se_ipw_ps <- sqrt(vcov(results_csme_ipw)[2, 2])
  coverage_ipw_ps <-
    1*(coef(results_ipw2)[7] - 1.96*se_ipw_ps < true_effect &
         coef(results_ipw2)[7] + 1.96*se_ipw_ps > true_effect)

  # Model 3: DR-CSME
  eefun_csme_aipw <- function(data, model) {
    Y <- data$Y
    Astar <- data$Astar
    L2 <- data$L2
    Lmat <- model.matrix(model, data = data)
    sw <- data$sw
    delta <- function(beta1, beta3, sigma_ep) {
      Astar + sigma_me*(beta1 + beta3*L2)*Y / sigma_ep
    }
    condexp <- function(beta0, beta1, beta2, beta3, sigma_ep) {
      (beta0 + (beta1 + beta3*L2)*delta(beta1, beta3, sigma_ep) +
         beta2*L2) /
        (1 + ((beta1 + beta3*L2)^2)*sigma_me / sigma_ep)[[1]]
    }
    condvar <- function(beta1, beta3, sigma_ep) {
      sigma_ep / (1 + ((beta1 + beta3*L2)^2)*sigma_me / sigma_ep)[[1]]
    }
    function(theta) {
      p  <- length(theta)
      p1 <- length(coef(model))
      rho <- Lmat %*% theta[1:p1]

      score_eqns <- apply(Lmat, 2, function(x) sum((Astar - rho) * x))

      c(score_eqns,
        Astar - theta[p-7],
        sw*(Y - condexp(theta[p-6], theta[p-5], theta[p-4], theta[p-3],
                        theta[p-2])),
        sw*(Y - condexp(theta[p-6], theta[p-5], theta[p-4], theta[p-3],
                        theta[p-2]))*L2,
        sw*(Y - condexp(theta[p-6], theta[p-5], theta[p-4], theta[p-3],
                        theta[p-2]))*
          delta(theta[p-5], theta[p-3], theta[p-2]),
        sw*(Y - condexp(theta[p-6], theta[p-5], theta[p-4], theta[p-3],
                        theta[p-2]))*
          L2*delta(theta[p-5], theta[p-3], theta[p-2]),
        sw*(theta[p-2] - theta[p-2]*
              (Y - condexp(theta[p-6], theta[p-5], theta[p-4], theta[p-3],
                           theta[p-2]))^2 /
              condvar(theta[p-5], theta[p-3], theta[p-2])),
        L2 - theta[p-1],
        theta[p-5] + theta[p-3]*theta[p-1] - theta[p]
      )
    }
  }

  failed <- TRUE
  j <- 1

  while(failed == TRUE & j < 6) {

    failed <- FALSE
    startvec <- c(coef(denom_mod), mean(Astar), coef(mod), sigma(mod)^2,
                  mean(L2), guess)*(j == 1) +
      c(coef(denom_mod), mean(Astar), rnorm(4, 0, j/5), sigma(mod)^2,
        mean(L2), guess)*(j > 1)
    results_csme_aipw <-
      tryCatch(m_estimate(estFUN = eefun_csme_aipw, data = data,
                          outer_args = list(denom_mod),
                          root_control =
                            setup_root_control(start = startvec)),
               error = function(e) { failed <<- TRUE})
    if (failed == FALSE) {
      if (abs(coef(results_csme_aipw)[11]) > 3) { failed <- TRUE }
    }
    j <- j + 1

  }

  bias_aipw_ps <- coef(results_csme_aipw)[11] - true_effect
  se_aipw_ps <- sqrt(vcov(results_csme_aipw)[11, 11])
  coverage_aipw_ps <-
    1*(coef(results_csme_aipw)[11] - 1.96*se_aipw_ps < true_effect &
       coef(results_csme_aipw)[11] + 1.96*se_aipw_ps > true_effect)

  if(failed == TRUE) {
    bias_aipw_ps <- NA
    se_aipw_ps <- NA
    coverage_aipw_ps <- NA
  }

  #####################################################################
  # Outcome regression correct, propensity model wrong


  # Fit regression to use for starting values
  mod <- lm(Y ~ Astar*L1 + Astar*L2, data = data)

  # Model 1: G-formula-CSME
  eefun_csme_gform <- function(data) {
    Y <- data$Y
    Astar <- data$Astar
    L1 <- data$L1
    L2 <- data$L2
    delta <- function(beta1, beta4, beta5, sigma_ep) {
      Astar + sigma_me*(beta1 + beta4*L1 + beta5*L2)*Y / sigma_ep
    }
    condexp <- function(beta0, beta1, beta2, beta3, beta4, beta5,
                        sigma_ep) {
      (beta0 + beta2*L1 + beta3*L2 + (beta1 + beta4*L1 + beta5*L2)*
          delta(beta1, beta4, beta5, sigma_ep)) /
          (1 + ((beta1 + beta4*L1 + beta5*L2)^2)*sigma_me / sigma_ep)[[1]]
    }
    condvar <- function(beta1, beta4, beta5, sigma_ep) {
      sigma_ep /
        (1 + ((beta1 + beta4*L1 + beta5*L2)^2)*sigma_me / sigma_ep)[[1]]
    }
    function(theta) {
      c((Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6], theta[7])),
        (Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6], theta[7]))*
          delta(theta[2], theta[5], theta[6], theta[7]),
        (Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6], theta[7]))*L1,
        (Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6], theta[7]))*L2,
        (Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6], theta[7]))*
          L1*delta(theta[2], theta[5], theta[6], theta[7]),
        (Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6], theta[7]))*
          L2*delta(theta[2], theta[5], theta[6], theta[7]),
        theta[7] - theta[7]*
          (Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6], theta[7]))^2 /
          condvar(theta[2], theta[5], theta[6], theta[7]),
        L1 - theta[8],
        L2 - theta[9],
        theta[2] + theta[5]*theta[8] + theta[6]*theta[9] - theta[10]
      )
    }
  }

  results_csme_gform <-
    m_estimate(estFUN = eefun_csme_gform, data = data,
               compute_roots = TRUE,
               root_control =
                 setup_root_control(start = c(coef(mod), sigma(mod)^2,
                                              mean(L1), mean(L2), guess)))

  bias_gform_or <- coef(results_csme_gform)[10] - true_effect
  se_gform_or <- sqrt(vcov(results_csme_gform)[10, 10])
  coverage_gform_or <-
    1*(coef(results_csme_gform)[10] - 1.96*se_gform_or < true_effect &
         coef(results_csme_gform)[10] + 1.96*se_gform_or > true_effect)

  # Weighted CSME
  # Estimate weights
  denom_mod <- lm(Astar ~ L2)
  p_denom <- predict(denom_mod, type='response')
  dens_denom <- dnorm(Astar, p_denom, summary(denom_mod)$sigma)
  num_mod <- lm(Astar ~ 1)
  p_num <- predict(num_mod, type='response')
  dens_num <- dnorm(Astar, p_num, summary(denom_mod)$sigma)
  data$sw <- dens_num / dens_denom

  # Fit weighted regression for starting values
  wmod <- lm(Y ~ Astar, weights = sw, data = data)

  # Model 2: IPW-CSME
  # Get point estimates and variance using geex
  eefun_ipw <- function(data, model) {
    Y <- data$Y
    Astar <- data$Astar
    Lmat <- model.matrix(model, data = data)
    sw <- data$sw
    delta <- function(beta1, sigma_ep) {
      Astar + sigma_me*beta1*Y / sigma_ep
    }
    condexp <- function(beta0, beta1, sigma_ep) {
      (beta0 + beta1*delta(beta1, sigma_ep)) /
        (1 + (beta1^2)*sigma_me / sigma_ep)
    }
    condvar <- function(beta1, sigma_ep) {
      sigma_ep /
        (1 + (beta1^2)*sigma_me / sigma_ep)
    }
    function(theta) {
      p  <- length(theta)
      p1 <- length(coef(model))
      rho <- Lmat %*% theta[1:p1]

      score_eqns <- apply(Lmat, 2, function(x) sum((Astar - rho) * x))

      c(score_eqns,
        Astar - theta[p-3],
        sw*(Y - condexp(theta[p-2], theta[p-1], theta[p])),
        sw*(Y - condexp(theta[p-2], theta[p-1], theta[p]))*
          delta(theta[p-1], theta[p]),
        sw*(theta[p] - theta[p]*
              (Y - condexp(theta[p-2], theta[p-1], theta[p]))^2 /
              condvar(theta[p-1], theta[p]))
      )
    }
  }

  failed <- TRUE
  j <- 1

  while(failed == TRUE & j < 6) {

    failed <- FALSE
    startvec <-
      c(coef(denom_mod), mean(Astar), coef(wmod)[1:2], sigma(wmod)^2)*
        (j == 1) +
      c(coef(denom_mod), mean(Astar), coef(mod)[1:2], sigma(mod)^2)*
        (j == 2) +
      c(coef(denom_mod), mean(Astar), rnorm(3, 0, j/5))*(j > 2)
    results_ipw <-
      tryCatch(m_estimate(estFUN = eefun_ipw, data = data,
                          outer_args = list(denom_mod),
                          root_control =
                            setup_root_control(start = startvec)),
               error = function(e) { failed <<- TRUE})
    if (failed == FALSE) {
      if (abs(coef(results_ipw)[5]) > 3) { failed <- TRUE }
    }
    j <- j + 1

  }

  bias_ipw_or <- coef(results_ipw)[5] - true_effect
  se_ipw_or <- sqrt(vcov(results_ipw)[5, 5])
  coverage_ipw_or <-
    1*(coef(results_ipw)[5] - 1.96*se_ipw_or < true_effect &
         coef(results_ipw)[5] + 1.96*se_ipw_or > true_effect)

  # Model 3: AIPW
  eefun_csme_aipw <- function(data, model) {
    Y <- data$Y
    Astar <- data$Astar
    L1 <- data$L1
    L2 <- data$L2
    Lmat <- model.matrix(model, data = data)
    sw <- data$sw
    delta <- function(beta1, beta4, beta5, sigma_ep) {
      Astar + sigma_me*(beta1 + beta4*L1 + beta5*L2)*Y / sigma_ep
    }
    condexp <- function(beta0, beta1, beta2, beta3, beta4, beta5,
                        sigma_ep) {
      (beta0 + beta2*L1 + beta3*L2 +
        (beta1 + beta4*L1 + beta5*L2)*delta(beta1, beta4, beta5,sigma_ep)) /
        (1 + ((beta1 + beta4*L1 + beta5*L2)^2)*sigma_me / sigma_ep)[[1]]
    }
    condvar <- function(beta1, beta4, beta5, sigma_ep) {
      sigma_ep /
        (1 + ((beta1 + beta4*L1 + beta5*L2)^2)*sigma_me / sigma_ep)[[1]]
    }
    function(theta) {
      p  <- length(theta)
      p1 <- length(coef(model))
      rho <- Lmat %*% theta[1:p1]

      score_eqns <- apply(Lmat, 2, function(x) sum((Astar - rho) * x))

      c(score_eqns,
        Astar - theta[p-10],
        sw*(Y - condexp(theta[p-9], theta[p-8], theta[p-7], theta[p-6],
                        theta[p-5], theta[p-4], theta[p-3])),
        sw*(Y - condexp(theta[p-9], theta[p-8], theta[p-7], theta[p-6],
                        theta[p-5], theta[p-4], theta[p-3]))*
          delta(theta[p-8], theta[p-5], theta[p-4], theta[p-3]),
        sw*(Y - condexp(theta[p-9], theta[p-8], theta[p-7], theta[p-6],
                        theta[p-5], theta[p-4], theta[p-3]))*L1,
        sw*(Y - condexp(theta[p-9], theta[p-8], theta[p-7], theta[p-6],
                        theta[p-5], theta[p-4], theta[p-3]))*L2,
        sw*(Y - condexp(theta[p-9], theta[p-8], theta[p-7], theta[p-6],
                        theta[p-5], theta[p-4], theta[p-3]))*
          L1*delta(theta[p-8], theta[p-5], theta[p-4], theta[p-3]),
        sw*(Y - condexp(theta[p-9], theta[p-8], theta[p-7], theta[p-6],
                        theta[p-5], theta[p-4], theta[p-3]))*
          L2*delta(theta[p-8], theta[p-5], theta[p-4], theta[p-3]),
        sw*(theta[p-3] - theta[p-3]*
              (Y - condexp(theta[p-9], theta[p-8], theta[p-7], theta[p-6],
                           theta[p-5], theta[p-4], theta[p-3]))^2 /
              condvar(theta[p-8], theta[p-5], theta[p-4], theta[p-3])),
        L1 - theta[p-2],
        L2 - theta[p-1],
        theta[p-8] + theta[p-5]*theta[p-2] +
          theta[p-4]*theta[p-1] - theta[p]
      )
    }
  }

  results_csme_aipw <-
    m_estimate(estFUN = eefun_csme_aipw, data = data, compute_roots = TRUE,
               outer_args = list(denom_mod),
               root_control =
                 setup_root_control(start = c(coef(denom_mod), mean(Astar),
                                              coef(mod), sigma(mod)^2,
                                              mean(L1), mean(L2), guess)))

  bias_aipw_or <- coef(results_csme_aipw)[13] - true_effect
  se_aipw_or <- sqrt(vcov(results_csme_aipw)[13, 13])
  coverage_aipw_or <-
    1*(coef(results_csme_aipw)[13] - 1.96*se_aipw_or < true_effect &
         coef(results_csme_aipw)[13] + 1.96*se_aipw_or > true_effect)

  #####################################################################
  # Both models correct

  # Fit linear regression to use for starting values
  mod <- lm(Y ~ Astar*L1 + Astar*L2, data = data)

  # Model 1: G-formula-CSME
  bias_gform_psor <- bias_gform_or
  se_gform_psor <- se_gform_or
  coverage_gform_psor <- coverage_gform_or

  # Estimate weights
  denom_mod <- lm(Astar ~ L1 + L2)
  p_denom <- predict(denom_mod, type='response')
  dens_denom <- dnorm(Astar, p_denom, summary(denom_mod)$sigma)
  num_mod <- lm(Astar ~ 1)
  p_num <- predict(num_mod, type='response')
  dens_num <- dnorm(Astar, p_num, summary(denom_mod)$sigma)
  data$sw <- dens_num / dens_denom

  # Model 2: IPW-CSE
  bias_ipw_psor <- bias_ipw_ps
  se_ipw_psor <- se_ipw_ps
  coverage_ipw_psor <- coverage_ipw_ps

  # Model 3: AIPW
  eefun_csme_aipw <- function(data, model) {
    Y <- data$Y
    Astar <- data$Astar
    L1 <- data$L1
    L2 <- data$L2
    Lmat <- model.matrix(model, data = data)
    sw <- data$sw
    delta <- function(beta1, beta4, beta5, sigma_ep) {
      Astar + sigma_me*(beta1 + beta4*L1 + beta5*L2)*Y / sigma_ep
    }
    condexp <- function(beta0, beta1, beta2, beta3, beta4, beta5,
                        sigma_ep) {
      (beta0 + beta2*L1 + beta3*L2 +
        (beta1 + beta4*L1 + beta5*L2)*delta(beta1, beta4, beta5,sigma_ep)) /
        (1 + ((beta1 + beta4*L1 + beta5*L2)^2)*sigma_me / sigma_ep)[[1]]
    }
    condvar <- function(beta1, beta4, beta5, sigma_ep) {
      sigma_ep /
        (1 + ((beta1 + beta4*L1 + beta5*L2)^2)*sigma_me / sigma_ep)[[1]]
    }
    function(theta) {
      p  <- length(theta)
      p1 <- length(coef(model))
      rho <- Lmat %*% theta[1:p1]

      score_eqns <- apply(Lmat, 2, function(x) sum((Astar - rho) * x))

      c(score_eqns,
        Astar - theta[p-10],
        sw*(Y - condexp(theta[p-9], theta[p-8], theta[p-7], theta[p-6],
                        theta[p-5], theta[p-4], theta[p-3])),
        sw*(Y - condexp(theta[p-9], theta[p-8], theta[p-7], theta[p-6],
                        theta[p-5], theta[p-4], theta[p-3]))*
          delta(theta[p-8], theta[p-5], theta[p-4], theta[p-3]),
        sw*(Y - condexp(theta[p-9], theta[p-8], theta[p-7], theta[p-6],
                        theta[p-5], theta[p-4], theta[p-3]))*L1,
        sw*(Y - condexp(theta[p-9], theta[p-8], theta[p-7], theta[p-6],
                        theta[p-5], theta[p-4], theta[p-3]))*L2,
        sw*(Y - condexp(theta[p-9], theta[p-8], theta[p-7], theta[p-6],
                        theta[p-5], theta[p-4], theta[p-3]))*
          L1*delta(theta[p-8], theta[p-5], theta[p-4], theta[p-3]),
        sw*(Y - condexp(theta[p-9], theta[p-8], theta[p-7], theta[p-6],
                        theta[p-5], theta[p-4], theta[p-3]))*
          L2*delta(theta[p-8], theta[p-5], theta[p-4], theta[p-3]),
        sw*(theta[p-3] - theta[p-3]*
              (Y - condexp(theta[p-9], theta[p-8], theta[p-7], theta[p-6],
                           theta[p-5], theta[p-4], theta[p-3]))^2 /
              condvar(theta[p-8], theta[p-5], theta[p-4], theta[p-3])),
        L1 - theta[p-2],
        L2 - theta[p-1],
        theta[p-8] + theta[p-5]*theta[p-2] +
          theta[p-4]*theta[p-1] - theta[p]
      )
    }
  }

  results_csme_aipw <-
    m_estimate(estFUN = eefun_csme_aipw, data = data, compute_roots = TRUE,
               outer_args = list(denom_mod),
               root_control =
                 setup_root_control(start = c(coef(denom_mod), mean(Astar),
                                              coef(mod), sigma(mod)^2,
                                              mean(L1), mean(L2), guess)))

  bias_aipw_psor <- coef(results_csme_aipw)[14] - true_effect
  se_aipw_psor <- sqrt(vcov(results_csme_aipw)[14, 14])
  coverage_aipw_psor <-
    1*(coef(results_csme_aipw)[14] - 1.96*se_aipw_psor < true_effect &
         coef(results_csme_aipw)[14] + 1.96*se_aipw_psor > true_effect)

  ######################################################################
  # Both models wrong

  # Fit regression to use for starting values
  mod <- lm(Y ~ Astar*L2, data = data)

  # Model 1: G-formula-CSME
  bias_gform_bothwrong <- bias_gform_ps
  se_gform_bothwrong <- se_gform_ps
  coverage_gform_bothwrong <- coverage_gform_ps

  # Model 2: IPW-CSME
  # Estimate weights
  denom_mod <- lm(Astar ~ L2)
  p_denom <- predict(denom_mod, type='response')
  dens_denom <- dnorm(Astar, p_denom, summary(denom_mod)$sigma)
  num_mod <- lm(Astar ~ 1)
  p_num <- predict(num_mod, type='response')
  dens_num <- dnorm(Astar, p_num, summary(denom_mod)$sigma)
  data$sw <- dens_num / dens_denom

  # Fit weighted regression for starting values
  bias_ipw_bothwrong <- bias_ipw_or
  se_ipw_bothwrong <- se_ipw_or
  coverage_ipw_bothwrong <- coverage_ipw_or

  # Model 3: DR-CSME
  eefun_csme_aipw <- function(data, model) {
    Y <- data$Y
    Astar <- data$Astar
    L2 <- data$L2
    Lmat <- model.matrix(model, data = data)
    sw <- data$sw
    delta <- function(beta1, beta3, sigma_ep) {
      Astar + sigma_me*(beta1 + beta3*L2)*Y / sigma_ep
    }
    condexp <- function(beta0, beta1, beta2, beta3, sigma_ep) {
      (beta0 + (beta1 + beta3*L2)*delta(beta1, beta3, sigma_ep) +
         beta2*L2) /
        (1 + ((beta1 + beta3*L2)^2)*sigma_me / sigma_ep)[[1]]
    }
    condvar <- function(beta1, beta3, sigma_ep) {
      sigma_ep / (1 + ((beta1 + beta3*L2)^2)*sigma_me / sigma_ep)[[1]]
    }
    function(theta) {
      p  <- length(theta)
      p1 <- length(coef(model))
      rho <- Lmat %*% theta[1:p1]

      score_eqns <- apply(Lmat, 2, function(x) sum((Astar - rho) * x))

      c(score_eqns,
        Astar - theta[p-7],
        sw*(Y - condexp(theta[p-6], theta[p-5], theta[p-4], theta[p-3],
                        theta[p-2])),
        sw*(Y - condexp(theta[p-6], theta[p-5], theta[p-4], theta[p-3],
                        theta[p-2]))*L2,
        sw*(Y - condexp(theta[p-6], theta[p-5], theta[p-4], theta[p-3],
                        theta[p-2]))*
          delta(theta[p-5], theta[p-3], theta[p-2]),
        sw*(Y - condexp(theta[p-6], theta[p-5], theta[p-4], theta[p-3],
                        theta[p-2]))*
          L2*delta(theta[p-5], theta[p-3], theta[p-2]),
        sw*(theta[p-2] - theta[p-2]*
              (Y - condexp(theta[p-6], theta[p-5], theta[p-4], theta[p-3],
                           theta[p-2]))^2 /
              condvar(theta[p-5], theta[p-3], theta[p-2])),
        L2 - theta[p-1],
        theta[p-5] + theta[p-3]*theta[p-1] - theta[p]
      )
    }
  }

  failed <- TRUE
  j <- 1

  while(failed == TRUE & j < 6) {

    failed <- FALSE
    startvec <- c(coef(denom_mod), mean(Astar), coef(mod), sigma(mod)^2,
                  mean(L2), guess)*(j == 1) +
      c(coef(denom_mod), mean(Astar), rnorm(4, 0, j/5), sigma(mod)^2,
        mean(L2), guess)*(j > 1)
    results_csme_aipw <-
      tryCatch(m_estimate(estFUN = eefun_csme_aipw, data = data,
                          outer_args = list(denom_mod),
                          root_control =
                            setup_root_control(start = startvec)),
               error = function(e) { failed <<- TRUE})
    if (failed == FALSE) {
      if (abs(coef(results_csme_aipw)[10]) > 3) { failed <- TRUE }
    }
    j <- j + 1

  }

  bias_aipw_bothwrong <- coef(results_csme_aipw)[10] - true_effect
  se_aipw_bothwrong <- sqrt(vcov(results_csme_aipw)[10, 10])
  coverage_aipw_bothwrong <-
    1*(coef(results_csme_aipw)[10] -
         1.96*se_aipw_bothwrong < true_effect &
       coef(results_csme_aipw)[10] +
         1.96*se_aipw_bothwrong > true_effect)

  if(failed == TRUE) {
    bias_aipw_bothwrong <- NA
    se_aipw_bothwrong <- NA
    coverage_aipw_bothwrong <- NA
  }


  return(c(bias_gform_ps, bias_ipw_ps, bias_aipw_ps,
           bias_gform_or, bias_ipw_or, bias_aipw_or,
           bias_gform_psor, bias_ipw_psor, bias_aipw_psor,
           bias_gform_bothwrong, bias_ipw_bothwrong, bias_aipw_bothwrong,
           se_gform_ps, se_ipw_ps, se_aipw_ps,
           se_gform_or, se_ipw_or, se_aipw_or,
           se_gform_psor, se_ipw_psor, se_aipw_psor,
           se_gform_bothwrong, se_ipw_bothwrong, se_aipw_bothwrong,
           coverage_gform_ps, coverage_ipw_ps, coverage_aipw_ps,
           coverage_gform_or, coverage_ipw_or, coverage_aipw_or,
           coverage_gform_psor, coverage_ipw_psor, coverage_aipw_psor,
           coverage_gform_bothwrong, coverage_ipw_bothwrong,
           coverage_aipw_bothwrong))

}

n <- 2000
low_n <- 400
trials <- seq(1, nsims)
combos <- data.frame(trials = rep(trials, length(beta1_true)),
                     mes = rep(sigma_me, each = nsims),
                     ns = rep(n, each = nsims))
i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
combo_i <- combos[(i), ]

set.seed(i*1000)
sim <- with(combo_i, mapply(simulator, trials, mes, ns))

# Output
outfile <- paste("./Results/results_scen3_", i, ".Rdata", sep = "")
save(sim, file = outfile)
