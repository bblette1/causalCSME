rm(list = ls())
library(geex)
library(rootSolve)

# Set parameter values
nsims <- 1000
n <- 2000
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

  # Case-cohort sample, Astar value when R,Y = 0 is arbitrary
  R <- rbinom(n, 1, 0.25)
  Astar[R == 0 & Y == 0] <- 0

  data <- data.frame("Y" = Y, "Astar" = Astar, "L1" = L1, "L2" = L2,
                     "R" = R)

  # Estimate case-cohort weights
  pi_hat <- mean(R[Y == 0])
  data$ccw <- (1-Y)*R/pi_hat + Y

  # Estimate weights
  denom_mod <- lm(Astar ~ L1 + L2, weights = data$ccw)
  p_denom <- predict(denom_mod, type='response')
  dens_denom <- dnorm(Astar, p_denom,
                      sd(residuals(denom_mod)[R == 1 | Y == 1]))
  num_mod <- lm(Astar ~ 1, weights = data$ccw)
  p_num <- predict(num_mod, type='response')
  dens_num <- dnorm(Astar, p_num,
                    sd(residuals(denom_mod)[R == 1 | Y == 1]))
  data$w <- rep(0, n)
  data$w[R == 1 | Y == 1] <- dens_num[R == 1 | Y == 1] /
                             dens_denom[R == 1 | Y == 1]
  data$sw <- data$w*data$ccw

  # Starting values
  mod <- glm(Y ~ Astar*L1 + Astar*L2, data = data, family = "binomial",
             weights = sw)

  # DR model estimating equations
  eefun_dr <- function(data) {
    Y <- data$Y
    Astar <- data$Astar
    L1 <- data$L1
    L2 <- data$L2
    sw <- data$sw
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
      c(sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4],
                         theta[5], theta[6])),
        sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4],
                         theta[5], theta[6]))*
          delta(theta[2], theta[5], theta[6]),
        sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4],
                         theta[5], theta[6]))*L1,
        sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4],
                         theta[5], theta[6]))*L2,
        sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4],
                         theta[5], theta[6]))*
          L1*delta(theta[2], theta[5], theta[6]),
        sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4],
                         theta[5], theta[6]))*
          L2*delta(theta[2], theta[5], theta[6])
      )
    }
  }

  results_dr <- m_estimate(estFUN = eefun_dr, data = data,
                           compute_roots = TRUE,
                           root_control =
                             setup_root_control(start = c(coef(mod))))

  # Evaluate over a grid of values from a = 0 to a = 4 every 0.1 value
  dr_ests <- rep(NA, length(seq(0, 4, 0.1)))
  cf <- coef(results_dr)

  for (i in seq(0, 4, 0.1)) {

    dr_ests[i*10 + 1] <-
      (1/n)*sum(exp(cf[1]+cf[2]*i+cf[3]*L1+cf[4]*L2+cf[5]*i*L1+cf[6]*i*L2) /
        (1 + exp(cf[1]+cf[2]*i+cf[3]*L1+cf[4]*L2 +cf[5]*i*L1+cf[6]*i*L2)))

  }

  # Get sandwich variance estimates for when a = 3
  # DR model estimating equations
  eefun_dr <- function(data) {
    Y <- data$Y
    Astar <- data$Astar
    L1 <- data$L1
    L2 <- data$L2
    sw <- data$sw
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
      c(sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4],
                        theta[5], theta[6])),
        sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4],
                        theta[5], theta[6]))*
          delta(theta[2], theta[5], theta[6]),
        sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4],
                        theta[5], theta[6]))*L1,
        sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4],
                        theta[5], theta[6]))*L2,
        sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4],
                        theta[5], theta[6]))*
          L1*delta(theta[2], theta[5], theta[6]),
        sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4],
                        theta[5], theta[6]))*
          L2*delta(theta[2], theta[5], theta[6]),
        theta[7] -
          exp(theta[1] + theta[2]*3 + theta[3]*L1 + theta[4]*L2 +
              theta[5]*3*L1 + theta[6]*3*L2) /
          (1 + exp(theta[1] + theta[2]*3 + theta[3]*L1 + theta[4]*L2 +
                   theta[5]*3*L1 + theta[6]*3*L2))
      )
    }
  }

  results_dr <- m_estimate(estFUN = eefun_dr, data = data,
                           compute_roots = FALSE,
                           roots = c(cf, dr_ests[31]))

  # Bias, SE, and coverage for a = 3
  bias_dr <- dr_ests[31] - true_effect_a3
  se_dr <- sqrt(vcov(results_dr)[7, 7])
  coverage_dr <- 1*(dr_ests[31] - 1.96*se_dr < true_effect_a3 &
                    dr_ests[31] + 1.96*se_dr > true_effect_a3)

  return(c(dr_ests, bias_dr, se_dr, coverage_dr))

}

trials <- seq(1, nsims)
combos <- data.frame(trials = rep(trials, length(beta1_true)),
                     betas = rep(beta1_true, each = nsims))
i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")) + 1000
combo_i <- combos[(i-1000), ]

set.seed(i*1000)
sim <- with(combo_i, mapply(simulator, trials, betas))

# Output
outfile <- paste("./Results/results_app1_3", i, ".Rdata", sep = "")
save(sim, file = outfile)
