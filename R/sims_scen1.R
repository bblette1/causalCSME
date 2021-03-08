rm(list = ls())
library(geex)
library(rootSolve)

# Set parameter values
nsims <- 1000
n <- 800
beta1_true <- 0.7
beta4_true <- -0.4
beta5_true <- -0.2
L1_prob <- 0.5
L2_prob <- 0.2
true_effect_a3 <-
  0.4*exp(-2 + 0.7*3) / (1 + exp(-2 + 0.7*3)) +
  0.1*exp(-1.6 + 0.5*3) / (1 + exp(-1.6 + 0.5*3)) +
  0.4*exp(-2.6 + 0.3*3) / (1 + exp(-2.6 + 0.3*3)) +
  0.1*exp(-2.2 + 0.1*3) / (1 + exp(-2.2 + 0.1*3))
sigma_me <- 0.25

# Simulation function to be sent to computing cluster
simulator <- function(trial, beta1_true) {

  L1 <- rbinom(n, 1, L1_prob)
  L2 <- rbinom(n, 1, L2_prob)
  A <- rnorm(n, 2 + 0.3*L1 - 0.5*L2, 0.6)
  Y_logit <- -2 + beta1_true*A - 0.6*L1 + 0.4*L2 + beta4_true*A*L1 +
             beta5_true*A*L2
  Y_prob <- exp(Y_logit) / (1 + exp(Y_logit))
  Y <- rbinom(n, 1, Y_prob)
  Astar <- A + rnorm(n, 0, sqrt(sigma_me))
  data <- data.frame("Y" = Y, "Astar" = Astar, "L1" = L1, "L2" = L2)

  # Model 1: Simple outcome regression
  mod <- glm(Y ~ Astar*L1 + Astar*L2, data = data, family = "binomial")

  # Bias and coverage for a = 3
  bias_lr <- coef(mod)[1] + 3*coef(mod)[2] - true_effect_a3
  se_lr <- sqrt(vcov(mod)[1, 1] + 9*vcov(mod)[2, 2] + 6*vcov(mod)[1, 2])
  coverage_lr <-
    1*(coef(mod)[1] + 3*coef(mod)[2] - 1.96*se_lr < true_effect_a3 &
       coef(mod)[1] + 3*coef(mod)[2] + 1.96*se_lr > true_effect_a3)

  # Model 2: Conditional score ME estimator (CSME)
  eefun_csme <- function(data) {
    Y <- data$Y
    Astar <- data$Astar
    L1 <- data$L1
    L2 <- data$L2
    delta <- function(beta1, beta4, beta5) {
      Astar + sigma_me*(beta1 + beta4*L1 + beta5*L2)*Y
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
                     theta[5], theta[6]))*
          delta(theta[2], theta[5], theta[6]),
        (Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6]))*L1,
        (Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6]))*L2,
        (Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6]))*
          L1*delta(theta[2], theta[5], theta[6]),
        (Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6]))*
          L2*delta(theta[2], theta[5], theta[6])
      )
    }
  }

  results_csme <- m_estimate(estFUN = eefun_csme, data = data,
                             compute_roots = TRUE,
                             root_control =
                               setup_root_control(start = c(coef(mod))))

  # Bias, SE, and coverage for a = 3
  bias_csme <- coef(results_csme)[1] + 3*coef(results_csme)[2] -
               true_effect_a3
  se_csme <- sqrt(vcov(results_csme)[1, 1] + 9*vcov(results_csme)[2, 2] +
                  6*vcov(results_csme)[1, 2])
  coverage_csme <-
    1*(coef(results_csme)[1] + 3*coef(results_csme)[2] - 1.96*se_csme <
         true_effect_a3 &
       coef(results_csme)[1] + 3*coef(results_csme)[2] + 1.96*se_csme >
         true_effect_a3)


  # For the two g-formula methods, evaluate over a grid of values
  # From a=0 to a=4 every 0.1 value

  # Model 3: G-formula, no ME correction and Model 4: G-formula CSME
  gform_ests <- rep(NA, length(seq(0, 4, 0.1)))
  gform_csme_ests <- rep(NA, length(seq(0, 4, 0.1)))

  for (i in seq(0, 4, 0.1)) {

    gform_ests[i*10 + 1] <-
      predict(mod, newdata = data.frame(Astar = i, L1 = 1, L2 = 1),
              type = "response")*mean(L1)*mean(L2) +
      predict(mod, newdata = data.frame(Astar = i, L1 = 1, L2 = 0),
              type = "response")*mean(L1)*(1 - mean(L2)) +
      predict(mod, newdata = data.frame(Astar = i, L1 = 0, L2 = 1),
              type = "response")*(1 - mean(L1))*mean(L2) +
      predict(mod, newdata = data.frame(Astar = i, L1 = 0, L2 = 0),
              type = "response")*(1 - mean(L1))*(1 - mean(L2))

    gform_csme_ests[i*10 + 1] <-
      exp(t(coef(results_csme)[1:6]) %*% c(1, i, 1, 1, i, i)) /
        (1 + exp(t(coef(results_csme)[1:6]) %*% c(1, i, 1, 1, i, i)))*
        mean(L1)*mean(L2) +
      exp(t(coef(results_csme)[1:6]) %*% c(1, i, 1, 0, i, 0)) /
        (1 + exp(t(coef(results_csme)[1:6]) %*% c(1, i, 1, 0, i, 0)))*
        mean(L1)*(1 - mean(L2)) +
      exp(t(coef(results_csme)[1:6]) %*% c(1, i, 0, 1, 0, i)) /
        (1 + exp(t(coef(results_csme)[1:6]) %*% c(1, i, 0, 1, 0, i)))*
        (1 - mean(L1))*mean(L2) +
      exp(t(coef(results_csme)[1:6]) %*% c(1, i, 0, 0, 0, 0)) /
        (1 + exp(t(coef(results_csme)[1:6]) %*% c(1, i, 0, 0, 0, 0)))*
        (1 - mean(L1))*(1 - mean(L2))

  }

  # Get sandwich variance estimates for each when a = 3
  eefun_gform <- function(data) {
    Y <- data$Y
    Astar <- data$Astar
    L1 <- data$L1
    L2 <- data$L2
    H <- function(x) {
      1 / (1 + exp(-x))
    }
    function(theta) {
      c(L1 - theta[1],
        L2 - theta[2],
        (Y - H(theta[3] + Astar*theta[4] + L1*theta[5] + L2*theta[6] +
              Astar*L1*theta[7] + Astar*L2*theta[8])) *
          c(1, Astar, L1, L2, Astar*L1, Astar*L2),
        H(t(c(1, 3, 1, 1, 3, 3)) %*% theta[3:8])*theta[1]*theta[2] +
          H(t(c(1, 3, 1, 0, 3, 0)) %*% theta[3:8])*theta[1]*(1-theta[2]) +
          H(t(c(1, 3, 0, 1, 0, 3)) %*% theta[3:8])*(1-theta[1])*theta[2] +
          H(t(c(1, 3, 0, 0, 0, 0)) %*% theta[3:8])*
          (1-theta[1])*(1-theta[2]) -
          theta[9]
      )
    }
  }

  results_gform <- m_estimate(estFUN = eefun_gform, data = data,
                              compute_roots = FALSE,
                              roots = c(mean(L1), mean(L2), coef(mod),
                                        gform_ests[31]))

  eefun_gform_csme <- function(data) {
    Y <- data$Y
    Astar <- data$Astar
    L1 <- data$L1
    L2 <- data$L2
    delta <- function(beta1, beta4, beta5) {
      Astar + sigma_me*(beta1 + beta4*L1 + beta5*L2)*Y
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
                     theta[5], theta[6]))*
          delta(theta[2], theta[5], theta[6]),
        (Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6]))*L1,
        (Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6]))*L2,
        (Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6]))*
          L1*delta(theta[2], theta[5], theta[6]),
        (Y - condexp(theta[1], theta[2], theta[3], theta[4],
                     theta[5], theta[6]))*
          L2*delta(theta[2], theta[5], theta[6]),
        L1 - theta[7],
        L2 - theta[8],
        H(t(c(1, 3, 1, 1, 3, 3)) %*% theta[1:6])*theta[7]*theta[8] +
          H(t(c(1, 3, 1, 0, 3, 0)) %*% theta[1:6])*theta[7]*(1-theta[8]) +
          H(t(c(1, 3, 0, 1, 0, 3)) %*% theta[1:6])*(1-theta[7])*theta[8] +
          H(t(c(1, 3, 0, 0, 0, 0)) %*% theta[1:6])*
          (1-theta[7])*(1-theta[8]) -
          theta[9]
      )
    }
  }

  results_gform_csme <- m_estimate(estFUN = eefun_gform_csme, data = data,
                                   compute_roots = FALSE,
                                   roots = c(coef(results_csme), mean(L1),
                                             mean(L2), gform_csme_ests[31]))

  # Bias, SE, and coverage for a = 3
  bias_gform <- gform_ests[31] - true_effect_a3
  se_gform <- sqrt(vcov(results_gform)[9, 9])
  coverage_gform <- 1*(gform_ests[31] - 1.96*se_gform < true_effect_a3 &
                       gform_ests[31] + 1.96*se_gform > true_effect_a3)

  bias_gform_csme <- gform_csme_ests[31] - true_effect_a3
  se_gform_csme <- sqrt(vcov(results_gform_csme)[9, 9])
  coverage_gform_csme <-
    1*(gform_csme_ests[31] - 1.96*se_gform_csme < true_effect_a3 &
       gform_csme_ests[31] + 1.96*se_gform_csme > true_effect_a3)

  return(c(coef(mod), coef(results_csme), gform_ests, gform_csme_ests,
           bias_lr, se_lr, coverage_lr, bias_csme, se_csme, coverage_csme,
           bias_gform, se_gform, coverage_gform, bias_gform_csme,
           se_gform_csme, coverage_gform_csme))

}

trials <- seq(1, nsims)
combos <- data.frame(trials = rep(trials, length(beta1_true)),
                     betas = rep(beta1_true, each = nsims))
i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")) + 1000
combo_i <- combos[(i-1000), ]

set.seed(i*1000)
sim <- with(combo_i, mapply(simulator, trials, betas))

# Output
outfile <- paste("./Results/results_scen1_", i, ".Rdata", sep = "")
save(sim, file = outfile)
