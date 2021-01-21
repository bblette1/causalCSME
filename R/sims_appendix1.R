rm(list = ls())
library(geex)
library(matrixStats)
library(rootSolve)

# Set parameter values
nsims <- 500
n <- 1250
beta1_true <- 0.7
beta3_true <- -0.4
beta5_true <- 0.3
L1_prob <- 0.5
L2_mean <- 1
true_effect <- beta1_true + beta3_true*L1_prob + beta5_true*L2_mean
sigma_me <- 0.09


simulator <- function(trial, sigma_me) {

  ####################################################################
  # Outcome regression wrong, propensity model correct

  # Generate data
  L1 <- rbinom(n, 1, L1_prob)
  L2 <- rnorm(n, L2_mean, 0.5)
  X <- rnorm(n, 2 + 0.4*L1 - 0.5*L2, 0.9)
  Y_logit <- -3 + beta1_true*X - 0.4*L1 + beta3_true*X*L1 -
    0.2*L2 + beta5_true*X*L2 + rnorm(n, 0, 0.2)
  Y_prob <- exp(Y_logit) / (1 + exp(Y_logit))
  Y <- rbinom(n, 1, Y_prob)
  Xstar <- X + rnorm(n, 0, sqrt(sigma_me))

  # Case-cohort sample, Xstar value when R,Y = 0 is arbitrary
  R <- rbinom(n, 1, 0.1)
  Xstar[R == 0 & Y == 0] <- 0

  data <- data.frame("Y" = Y, "Xstar" = Xstar, "L1" = L1, "L2" = L2,
                     "R" = R)

  # Estimate case-cohort weights
  pi_hat <- mean(R[Y == 0])
  data$ccw <- (1-Y)*R/pi_hat + Y

  # Initial guess for causal effect (good but variable intuition)
  guess <- true_effect + rnorm(1, 0, 0.1)

  # Fit logistic regression to use for starting values
  mod <- glm(Y ~ Xstar*L2, family = "binomial", weights = data$ccw)

  # Model 1: G-formula-CSME
  eefun_csme_gform <- function(data) {
    Y <- data$Y
    Xstar <- data$Xstar
    L2 <- data$L2
    ccw <- data$ccw
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
      c(ccw*(Y - condexp(theta[1], theta[2], theta[3], theta[4])),
        ccw*(Y - condexp(theta[1], theta[2], theta[3], theta[4]))*L2,
        ccw*(Y - condexp(theta[1], theta[2], theta[3], theta[4]))*
          delta(theta[2], theta[4]),
        ccw*(Y - condexp(theta[1], theta[2], theta[3], theta[4]))*
          L2*delta(theta[2], theta[4]),
        L2 - theta[5],
        theta[2] + theta[4]*theta[5] - theta[6]
      )
    }
  }

  # Solve the estimating equations in a try-catch system
  # Almost always works on first try in main sims but Appendix 1 sims
  # have very small effective sample size, leading to more issues
  # Tries new starting values and sets to NA if always divergent
  failed <- TRUE
  j <- 1

  while(failed == TRUE & j < 6) {

    failed <- FALSE
    startvec <- c(coef(mod), mean(L2), guess)*(j == 1) +
      c(coef(mod) + rnorm(4, 0, j/5), mean(L2),
        guess + rnorm(1, 0, j/10))*(j > 1)
    results_csme_gform <-
      tryCatch(m_estimate(estFUN = eefun_csme_gform, data = data,
                          root_control =
                            setup_root_control(start = startvec)),
               error = function(e) { failed <<- TRUE})
    if (failed == FALSE) {
      if (abs(coef(results_csme_gform)[6]) > 3) { failed <- TRUE }
    }
    j <- j + 1

  }

  if (failed == FALSE) {
    bias_gform_ps <- coef(results_csme_gform)[6] - true_effect
    se_gform_ps <- sqrt(vcov(results_csme_gform)[6, 6])
    coverage_gform_ps <-
      1*(coef(results_csme_gform)[6] - 1.96*se_gform_ps < true_effect &
           coef(results_csme_gform)[6] + 1.96*se_gform_ps > true_effect)
  }

  if(failed == TRUE) {
    bias_gform_ps <- NA
    se_gform_ps <- NA
    coverage_gform_ps <- NA
  }

  # Weighted CSME
  # Estimate weights
  denom_mod <- lm(Xstar ~ L1 + L2, weights = data$ccw)
  p_denom <- predict(denom_mod, type='response')
  dens_denom <- dnorm(Xstar, p_denom, summary(denom_mod)$sigma)
  num_mod <- lm(Xstar ~ 1, weights = data$ccw)
  p_num <- predict(num_mod, type='response')
  dens_num <- dnorm(Xstar, p_num, summary(denom_mod)$sigma)
  data$w <- dens_num / dens_denom
  data$sw <- data$w*data$ccw

  # Fit weighted regression for starting values
  wmod <- glm(Y ~ Xstar, weights = sw, family = "binomial", data = data)

  # Model 2: IPW-CSE
  # Get point estimates and variance using geex
  eefun_csme_ipw <- function(data) {
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

  # Solve the estimating equations in a try-catch system
  failed <- TRUE
  j <- 1

  while(failed == TRUE & j < 6) {

    failed <- FALSE
    startvec <- coef(wmod)[1:2]*(j == 1) +
                coef(mod)[1:2]*(j == 2) +
                rnorm(2, c(0, guess), j/5)*(j > 2)
    results_csme_ipw <-
      tryCatch(m_estimate(estFUN = eefun_csme_ipw, data = data,
                          root_control =
                            setup_root_control(start = startvec)),
               error = function(e) { failed <<- TRUE})
    if (failed == FALSE) {
      if (abs(coef(results_csme_ipw)[2]) > 3) { failed <- TRUE }
    }
    j <- j + 1

  }

  if (failed == FALSE) {
    bias_ipw_ps <- coef(results_csme_ipw)[2] - true_effect
    se_ipw_ps <- sqrt(vcov(results_csme_ipw)[2, 2])
    coverage_ipw_ps <-
      1*(coef(results_csme_ipw)[2] - 1.96*se_ipw_ps < true_effect &
           coef(results_csme_ipw)[2] + 1.96*se_ipw_ps > true_effect)
  }

  if(failed == TRUE) {
    bias_ipw_ps <- NA
    se_ipw_ps <- NA
    coverage_ipw_ps <- NA
  }

  # Model 3: AIPW
  # Weighted version of model 1
  eefun_csme_aipw <- function(data) {
    Y <- data$Y
    Xstar <- data$Xstar
    L2 <- data$L2
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
      c(sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4])),
        sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4]))*L2,
        sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4]))*
          delta(theta[2], theta[4]),
        sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4]))*
          L2*delta(theta[2], theta[4]),
        L2 - theta[5],
        theta[2] + theta[4]*theta[5] - theta[6]
      )
    }
  }

  failed <- TRUE
  j <- 1

  while(failed == TRUE & j < 6) {

    failed <- FALSE
    startvec <- c(coef(mod), mean(L2), guess)*(j == 1) +
                c(rnorm(4, 0, j/5), mean(L2), guess)*(j > 1)
    results_csme_aipw <-
      tryCatch(m_estimate(estFUN = eefun_csme_aipw, data = data,
                          root_control =
                            setup_root_control(start = startvec)),
               error = function(e) { failed <<- TRUE})
    if (failed == FALSE) {
      if (abs(coef(results_csme_aipw)[6]) > 2) { failed <- TRUE }
    }
    j <- j + 1

  }

  if (failed == FALSE) {
    bias_aipw_ps <- coef(results_csme_aipw)[6] - true_effect
    se_aipw_ps <- sqrt(vcov(results_csme_aipw)[6, 6])
    coverage_aipw_ps <-
      1*(coef(results_csme_aipw)[6] - 1.96*se_aipw_ps < true_effect &
           coef(results_csme_aipw)[6] + 1.96*se_aipw_ps > true_effect)
  }

  if(failed == TRUE) {
    bias_aipw_ps <- NA
    se_aipw_ps <- NA
    coverage_aipw_ps <- NA
  }

  #####################################################################
  # Outcome regression correct, propensity model wrong

  # Fit logistic regression to use for starting values
  mod <- glm(Y ~ Xstar*L1 + Xstar*L2, family = "binomial",
             weights = data$ccw)

  # Model 1: G-formula-CSME
  eefun_csme_gform <- function(data) {
    Y <- data$Y
    Xstar <- data$Xstar
    L1 <- data$L1
    L2 <- data$L2
    ccw <- data$ccw
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
      c(ccw*(Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6])),
        ccw*(Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6]))*L1,
        ccw*(Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6]))*L2,
        ccw*(Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6]))*
          delta(theta[2], theta[5], theta[6]),
        ccw*(Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6]))*
          L1*delta(theta[2], theta[5], theta[6]),
        ccw*(Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6]))*
          L2*delta(theta[2], theta[5], theta[6]),
        L1 - theta[7],
        L2 - theta[8],
        theta[2] + theta[5]*theta[7] + theta[6]*theta[8] - theta[9]
      )
    }
  }

  failed <- TRUE
  j <- 1

  while(failed == TRUE & j < 6) {

    failed <- FALSE
    startvec <- c(coef(mod), mean(L1), mean(L2), guess)*(j == 1) +
      c(coef(mod) + rnorm(6, 0, j/5), mean(L1), mean(L2),
        guess + rnorm(1, 0, j/10))*(j > 1)
    results_csme_gform <-
      tryCatch(m_estimate(estFUN = eefun_csme_gform, data = data,
                          root_control =
                            setup_root_control(start = startvec)),
               error = function(e) { failed <<- TRUE})
    if (failed == FALSE) {
      if (abs(coef(results_csme_gform)[9]) > 3) { failed <- TRUE }
    }
    j <- j + 1

  }

  if (failed == FALSE) {
    bias_gform_or <- coef(results_csme_gform)[9] - true_effect
    se_gform_or <- sqrt(vcov(results_csme_gform)[9, 9])
    coverage_gform_or <-
      1*(coef(results_csme_gform)[9] - 1.96*se_gform_or < true_effect &
           coef(results_csme_gform)[9] + 1.96*se_gform_or > true_effect)
  }

  if(failed == TRUE) {
    bias_gform_or <- NA
    se_gform_or <- NA
    coverage_gform_or <- NA
  }

  # Weighted CSME
  # Estimate weights
  denom_mod <- lm(Xstar ~ L2)
  p_denom <- predict(denom_mod, type='response')
  dens_denom <- dnorm(Xstar, p_denom, summary(denom_mod)$sigma)
  num_mod <- lm(Xstar ~ 1)
  p_num <- predict(num_mod, type='response')
  dens_num <- dnorm(Xstar, p_num, summary(denom_mod)$sigma)
  data$w <- dens_num / dens_denom
  data$sw <- data$w*data$ccw

  # Fit weighted regression for starting values
  wmod <- glm(Y ~ Xstar, weights = data$sw, family = "binomial")

  # Model 2: IPW-CSE
  # Get point estimates and variance using geex
  eefun_csme_ipw <- function(data) {
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

  failed <- TRUE
  j <- 1

  while(failed == TRUE & j < 6) {

    failed <- FALSE
    startvec <- coef(wmod)[1:2]*(j == 1) +
                coef(mod)[1:2]*(j == 2) +
                rnorm(2, c(0, guess), j/5)*(j > 2)
    results_ipw <-
      tryCatch(m_estimate(estFUN = eefun_csme_ipw, data = data,
                          root_control =
                            setup_root_control(start = startvec)),
               error = function(e) { failed <<- TRUE})
    if (failed == FALSE) {
      if (abs(coef(results_ipw)[2]) > 3) { failed <- TRUE }
    }
    j <- j + 1

  }

  if (failed == FALSE) {
    bias_ipw_or <- coef(results_ipw)[2] - true_effect
    se_ipw_or <- sqrt(vcov(results_ipw)[2, 2])
    coverage_ipw_or <-
      1*(coef(results_ipw)[2] - 1.96*se_ipw_or < true_effect &
           coef(results_ipw)[2] + 1.96*se_ipw_or > true_effect)
  }

  if(failed == TRUE) {
    bias_ipw_or <- NA
    se_ipw_or <- NA
    coverage_ipw_or <- NA
  }

  # Model 3: AIPW
  # Include denominator density in g-formula CSME outcome model
  eefun_csme_aipw <- function(data) {
    Y <- data$Y
    Xstar <- data$Xstar
    L1 <- data$L1
    L2 <- data$L2
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
      c(sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6])),
        sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6]))*L1,
        sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6]))*L2,
        sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6]))*
          delta(theta[2], theta[5], theta[6]),
        sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6]))*
          L1*delta(theta[2], theta[5], theta[6]),
        sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6]))*
          L2*delta(theta[2], theta[5], theta[6]),
        L1 - theta[7],
        L2 - theta[8],
        theta[2] + theta[5]*theta[7] + theta[6]*theta[8] - theta[9]
      )
    }
  }

  failed <- TRUE
  j <- 1

  while(failed == TRUE & j < 6) {

    failed <- FALSE
    startvec <- c(coef(mod), mean(L1), mean(L2), guess)*(j == 1) +
                c(rnorm(6, coef(mod), j/5), mean(L1),
                  mean(L2), rnorm(1, guess, j/10))*(j > 1)
    results_csme_aipw <-
      tryCatch(m_estimate(estFUN = eefun_csme_aipw, data = data,
                          root_control =
                            setup_root_control(start = startvec)),
               error = function(e) { failed <<- TRUE})
    if (failed == FALSE) {
      if (abs(coef(results_csme_aipw)[9]) > 3) { failed <- TRUE }
    }
    j <- j + 1

  }

  if (failed == FALSE) {
    bias_aipw_or <- coef(results_csme_aipw)[9] - true_effect
    se_aipw_or <- sqrt(vcov(results_csme_aipw)[9, 9])
    coverage_aipw_or <-
      1*(coef(results_csme_aipw)[9] - 1.96*se_aipw_or < true_effect &
           coef(results_csme_aipw)[9] + 1.96*se_aipw_or > true_effect)
  }

  if(failed == TRUE) {
    bias_aipw_or <- NA
    se_aipw_or <- NA
    coverage_aipw_or <- NA
  }

  #####################################################################
  # Both models correct

  # Fit linear regression to use for starting values
  mod <- glm(Y ~ Xstar*L1 + Xstar*L2, family = "binomial",
             weights = data$ccw)

  # Model 1: G-formula-CSME
  eefun_csme_gform <- function(data) {
    Y <- data$Y
    Xstar <- data$Xstar
    L1 <- data$L1
    L2 <- data$L2
    ccw <- data$ccw
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
      c(ccw*(Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6])),
        ccw*(Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6]))*L1,
        ccw*(Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6]))*L2,
        ccw*(Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6]))*
          delta(theta[2], theta[5], theta[6]),
        ccw*(Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6]))*
          L1*delta(theta[2], theta[5], theta[6]),
        ccw*(Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6]))*
          L2*delta(theta[2], theta[5], theta[6]),
        L1 - theta[7],
        L2 - theta[8],
        theta[2] + theta[5]*theta[7] + theta[6]*theta[8] - theta[9]
      )
    }
  }

  failed <- TRUE
  j <- 1

  while(failed == TRUE & j < 6) {

    failed <- FALSE
    startvec <- c(coef(mod), mean(L1), mean(L2), guess)*(j == 1) +
                c(coef(mod) + rnorm(6, 0, j/5), mean(L1), mean(L2),
                  guess + rnorm(1, 0, j/10))*(j > 1)
    results_csme_gform <-
      tryCatch(m_estimate(estFUN = eefun_csme_gform, data = data,
                          root_control =
                            setup_root_control(start = startvec)),
               error = function(e) { failed <<- TRUE})
    if (failed == FALSE) {
      if (abs(coef(results_csme_gform)[9]) > 3) { failed <- TRUE }
    }
    j <- j + 1

  }

  if (failed == FALSE) {
    bias_gform_psor <- coef(results_csme_gform)[9] - true_effect
    se_gform_psor <- sqrt(vcov(results_csme_gform)[9, 9])
    coverage_gform_psor <-
      1*(coef(results_csme_gform)[9] - 1.96*se_gform_psor < true_effect &
           coef(results_csme_gform)[9] + 1.96*se_gform_psor > true_effect)
  }

  if(failed == TRUE) {
    bias_gform_psor <- NA
    se_gform_psor <- NA
    coverage_gform_psor <- NA
  }

  # Weighted CSME
  # Estimate weights
  denom_mod <- lm(Xstar ~ L1 + L2)
  p_denom <- predict(denom_mod, type='response')
  dens_denom <- dnorm(Xstar, p_denom, summary(denom_mod)$sigma)
  num_mod <- lm(Xstar ~ 1)
  p_num <- predict(num_mod, type='response')
  dens_num <- dnorm(Xstar, p_num, summary(denom_mod)$sigma)
  data$w <- dens_num / dens_denom
  data$sw <- data$w*data$ccw

  # Fit weighted regression for starting values
  wmod <- glm(Y ~ Xstar, weights = data$sw, family = "binomial")

  # Model 2: IPW-CSE
  # Get point estimates and variance using geex
  eefun_csme_ipw <- function(data) {
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

  failed <- TRUE
  j <- 1

  while(failed == TRUE & j < 6) {

    failed <- FALSE
    startvec <- coef(wmod)[1:2]*(j == 1) +
                coef(mod)[1:2]*(j == 2) +
                rnorm(2, coef(wmod), j/5)*(j > 2)
    results_csme_ipw <-
      tryCatch(m_estimate(estFUN = eefun_csme_ipw, data = data,
                          root_control =
                            setup_root_control(start = startvec)),
               error = function(e) { failed <<- TRUE})
    if (failed == FALSE) {
      if (abs(coef(results_csme_ipw)[2]) > 3) { failed <- TRUE }
    }
    j <- j + 1

  }

  if (failed == FALSE) {
    bias_ipw_psor <- coef(results_csme_ipw)[2] - true_effect
    se_ipw_psor <- sqrt(vcov(results_csme_ipw)[2, 2])
    coverage_ipw_psor<-
      1*(coef(results_csme_ipw)[2] - 1.96*se_ipw_psor < true_effect &
           coef(results_csme_ipw)[2] + 1.96*se_ipw_psor > true_effect)
  }

  if(failed == TRUE) {
    bias_ipw_psor <- NA
    se_ipw_psor <- NA
    coverage_ipw_psor <- NA
  }

  # Model 3: AIPW
  # Include denominator density in g-formula CSME outcome model
  eefun_csme_aipw <- function(data) {
    Y <- data$Y
    Xstar <- data$Xstar
    L1 <- data$L1
    L2 <- data$L2
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
      c(sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4],
                        theta[5], theta[6])),
        sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4],
                        theta[5], theta[6]))*L1,
        sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4],
                        theta[5], theta[6]))*L2,
        sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4],
                        theta[5], theta[6]))*
          delta(theta[2], theta[5], theta[6]),
        sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4],
                        theta[5], theta[6]))*
          L1*delta(theta[2], theta[5], theta[6]),
        sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4],
                        theta[5], theta[6]))*
          L2*delta(theta[2], theta[5], theta[6]),
        L1 - theta[7],
        L2 - theta[8],
        theta[2] + theta[5]*theta[7] + theta[6]*theta[8] - theta[9]
      )
    }
  }

  failed <- TRUE
  j <- 1

  while(failed == TRUE & j < 6) {

    failed <- FALSE
    startvec <- c(coef(mod), mean(L1), mean(L2), guess)*(j == 1) +
                c(rnorm(6, coef(mod), j/5), mean(L1),
                  mean(L2), rnorm(1, guess, j/10))*(j > 1)
    results_csme_aipw <-
      tryCatch(m_estimate(estFUN = eefun_csme_aipw, data = data,
                          root_control =
                            setup_root_control(start = startvec)),
               error = function(e) { failed <<- TRUE})
    if (failed == FALSE) {
      if (abs(coef(results_csme_aipw)[9]) > 3) { failed <- TRUE }
    }
    j <- j + 1

  }

  if (failed == FALSE) {
    bias_aipw_psor <- coef(results_csme_aipw)[9] - true_effect
    se_aipw_psor <- sqrt(vcov(results_csme_aipw)[9, 9])
    coverage_aipw_psor <-
      1*(coef(results_csme_aipw)[9] - 1.96*se_aipw_psor < true_effect &
           coef(results_csme_aipw)[9] + 1.96*se_aipw_psor > true_effect)
  }

  if (failed == TRUE) {
    bias_aipw_psor <- NA
    se_aipw_psor <- NA
    coverage_aipw_psor <- NA
  }


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
i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")) + 500
combo_i <- combos[(i - 500), ]

set.seed(i*1000)
sim <- with(combo_i, mapply(simulator, trials, mes))

# Output
outfile <- paste("./Results/results_app1_1250_", i, ".Rdata", sep = "")
save(sim, file = outfile)
