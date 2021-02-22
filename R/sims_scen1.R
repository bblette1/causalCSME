rm(list = ls())
library(geex)
library(rootSolve)

# Set parameter values
nsims <- 500
n <- 800
beta1_true <- 0.7
beta4_true <- -0.4
beta5_true <- -0.2
L1_prob <- 0.5
L2_prob <- 0.2
true_effect <- NA
sigma_me <- 0.36

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
          L2*delta(theta[2], theta[5], theta[6])
      )
    }
  }

  results_csme <- m_estimate(estFUN = eefun_csme, data = data,
                             compute_roots = TRUE,
                             root_control =
                               setup_root_control(start = c(coef(mod))))

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

    gform_deriv <- c(1, i, mean(L1), mean(L2), i*mean(L1), i*mean(L2))

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

  return(c(coef(mod), coef(results_csme), gform_ests, gform_csme_ests))

}

outmat <- matrix(rep(NA, 100*94), nrow = 100)
for (j in 1:100) {
  outmat[j, ] <- simulator(1, 0.7)
}
plot(testvec, testout, type = "l", lty = 1, ylim = c(0, 0.7),
     xlab = "a", ylab = "E[Y(a)]")
lines(testvec, colMeans(outmat)[13:53], col = "blue")
lines(testvec, colMeans(outmat)[54:94], col = "green")
lines(testvec, exp(colMeans(outmat)[7] + colMeans(outmat)[8]*testvec) /
        (1 + exp(colMeans(outmat)[7] + colMeans(outmat)[8]*testvec)),
      col = "orange")
lines(testvec, exp(colMeans(outmat)[1] + colMeans(outmat)[2]*testvec) /
        (1 + exp(colMeans(outmat)[1] + colMeans(outmat)[2]*testvec)),
      col = "red")
legend(0, 0.7,
       legend=c("Oracle", "LR", "CSME", "G-formula", "G-formula CSME"),
       col=c("black", "red", "orange", "blue", "green"), cex = 0.6, lty = 1)

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
