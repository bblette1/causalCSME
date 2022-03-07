rm(list = ls())
library(geex)
library(rootSolve)
library(survey)

nsims <- 1000
n <- 800
beta1 <- 0.4
beta2 <- -0.6
beta3 <- -0.4
beta4 <- 0.7
beta5 <- -0.6
beta6 <- -0.9
lambda <- 3
sigma_me1 <- 0.36
sigma_me3 <- 0.25

simulator <- function(trial, beta1, n) {

  L <- rexp(n, lambda)
  # A2 and A3 switch names in the manuscript
  A1 <- rnorm(n, 4 + 0.9*L, 1.1)
  A2 <- rnorm(n, 1.4 + 0.5*L, 0.6)
  A3 <- rnorm(n, 2.5, 0.7)
  Y_prob <- exp(-1.7 + beta1*A1 + beta2*A2 + beta3*A3 + beta4*L +
                  beta5*A1*L + beta6*A2*L -
                  log(lambda / (lambda - beta4 - beta5*A1 - beta6*A2))) /
            (1 + exp(-1.7 + beta1*A1 + beta2*A2 + beta3*A3))
  # Correct Y_prob for very rare instances > 1
  Y_prob[Y_prob > 1] <- 0.999
  Y <- rbinom(n, 1, Y_prob)
  A1star <- A1 + rnorm(n, 0, sqrt(sigma_me1))
  A3star <- A3 + rnorm(n, 0, sqrt(sigma_me3))
  data <- data.frame("Y" = Y, "A1star" = A1star, "A2" = A2,
                     "A3star" = A3star, "L" = L)

  # Model 1: Simple outcome regression
  mod <- glm(Y ~ A1star*L + A2*L + A3star, family = "binomial", data = data)

  bias1_lm <- coef(mod)[[2]] - beta1
  se1_lm <- summary(mod)$coefficients[2, 2]
  coverage1_lm <- 1*(coef(mod)[[2]] - 1.96*se1_lm < (beta1) &
                       coef(mod)[[2]] + 1.96*se1_lm > (beta1))

  bias2_lm <- coef(mod)[[4]] - (beta2)
  se2_lm <- summary(mod)$coefficients[4, 2]
  coverage2_lm <- 1*(coef(mod)[[4]] - 1.96*se2_lm < (beta2) &
                       coef(mod)[[4]] + 1.96*se2_lm > (beta2))

  bias3_lm <- coef(mod)[[5]] - beta3
  se3_lm <- summary(mod)$coefficients[5, 2]
  coverage3_lm <- 1*(coef(mod)[[5]] - 1.96*se3_lm < beta3 &
                       coef(mod)[[5]] + 1.96*se3_lm > beta3)

  # Model 2: Conditional score ME estimator (CSME)
  eefun_csme <- function(data) {
    Y <- data$Y
    A1star <- data$A1star
    A2 <- data$A2
    A3star <- data$A3star
    L <- data$L
    delta1 <- function(beta1, beta5) {
      A1star + (beta1 + beta5*L)*sigma_me1*Y
    }
    delta3 <- function(beta4) {
      A3star + beta4*sigma_me3*Y
    }
    H <- function(x) {
      1 / (1 + exp(-x))
    }
    condexp <- function(beta0, beta1, beta2, beta3, beta4, beta5, beta6) {
      H(beta0 + beta2*L + beta3*A2 + beta6*A2*L +
         (beta1 + beta5*L)*delta1(beta1, beta5) + beta4*delta3(beta4) -
         ((beta1 + beta5*L)^2 * sigma_me1 + beta4^2 * sigma_me3) / 2)
    }
    function(theta) {
      c((Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                     theta[6], theta[7])),
        (Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                     theta[6], theta[7]))*delta1(theta[2], theta[6]),
        (Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                     theta[6], theta[7]))*L,
        (Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                     theta[6], theta[7]))*A2,
        (Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                     theta[6], theta[7]))*delta3(theta[5]),
        (Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                     theta[6], theta[7]))*L*delta1(theta[2], theta[6]),
        (Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                     theta[6], theta[7]))*L*A2
      )
    }
  }

  results_csme <- m_estimate(estFUN = eefun_csme, data = data,
                             compute_roots = TRUE,
                             root_control =
                               setup_root_control(start = coef(mod)))

  bias1_csme <- coef(results_csme)[2] - beta1
  se1_csme <- sqrt(vcov(results_csme)[2, 2])
  coverage1_csme <-
    1*(coef(results_csme)[2] - 1.96*se1_csme < beta1 &
       coef(results_csme)[2] + 1.96*se1_csme > beta1)
  bias2_csme <- coef(results_csme)[4] - beta2
  se2_csme <- sqrt(vcov(results_csme)[4, 4])
  coverage2_csme <-
    1*(coef(results_csme)[4] - 1.96*se2_csme < beta2 &
       coef(results_csme)[4] + 1.96*se2_csme > beta2)
  bias3_csme <- coef(results_csme)[5] - beta3
  se3_csme <- sqrt(vcov(results_csme)[5, 5])
  coverage3_csme <- 1*(coef(results_csme)[5] - 1.96*se3_csme < beta3 &
                       coef(results_csme)[5] + 1.96*se3_csme > beta3)

  # Model 3: IPW estimator (no ME correction)
  # Estimate weights
  denom_mod1 <- lm(A1star ~ L)
  p_denom1 <- predict(denom_mod1, type='response')
  dens_denom1 <- dnorm(A1star, p_denom1, summary(denom_mod1)$sigma)
  num_mod1 <- lm(A1star ~ 1)
  p_num1 <- predict(num_mod1, type='response')
  dens_num1 <- dnorm(A1star, p_num1, summary(denom_mod1)$sigma)
  data$sw1 <- dens_num1 / dens_denom1

  denom_mod2 <- lm(A2 ~ L)
  p_denom2 <- predict(denom_mod2, type='response')
  dens_denom2 <- dnorm(A2, p_denom2, summary(denom_mod2)$sigma)
  num_mod2 <- lm(A2 ~ 1)
  p_num2 <- predict(num_mod2, type='response')
  dens_num2 <- dnorm(A2, p_num2, summary(denom_mod2)$sigma)
  data$sw2 <- dens_num2 / dens_denom2

  data$sw <- data$sw1*data$sw2

  # Fit weighted regression
  #wmod <- glm(Y ~ A1star + A2 + A3star, weights = data$sw,
              #family = "binomial")
  wmod <- svyglm(Y ~ A1star + A2 + A3star,
                 design = svydesign(id = ~1, weights=data$sw, data=data),
                 family = "binomial")

  bias1_ipw <- coef(wmod)[2] - (beta1)
  bias2_ipw <- coef(wmod)[3] - (beta2)
  bias3_ipw <- coef(wmod)[4] - beta3

  run <- F
  if (run == T) {
  eefun_wr <- function(data, model1, model2, model3, model4) {
    Y <- data$Y
    A3star <- data$A3star
    A1star <- model.response(model.frame(model1, data = data))
    L1 <- model.matrix(model1, data = data)
    V1 <- model.matrix(model2, data = data)
    A2 <- model.response(model.frame(model3, data = data))
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

      dens_denom1 <- dnorm(A1star, rho1, sqrt(theta[p-6]))
      dens_num1 <- dnorm(A1star, rho2, sqrt(theta[p-6]))
      dens_denom2 <- dnorm(A2, rho3, sqrt(theta[p-5]))
      dens_num2 <- dnorm(A2, rho4, sqrt(theta[p-5]))
      sw1 <- dens_num1 / dens_denom1
      sw2 <- dens_num2 / dens_denom2
      sw <- sw1*sw2

      score_eqns1 <- apply(L1, 2, function(x) sum((A1star - rho1) * x))
      score_eqns2 <- apply(V1, 2, function(x) sum((A1star - rho2) * x))
      score_eqns3 <- apply(L2, 2, function(x) sum((A2 - rho3) * x))
      score_eqns4 <- apply(V2, 2, function(x) sum((A2 - rho4) * x))

      c((A1star - rho1)*L1[,1],
        (A1star - rho1)*L1[,2],
        (A1star - rho2)*V1[,1],
        (A2 - rho3)*L2[,1],
        (A2 - rho3)*L2[,2],
        (A2 - rho4)*V2[,1],
        (n-p1)/n * theta[p-6] - (A1star - L1 %*% theta[1:p1])^2,
        (n-p3)/n * theta[p-5] - (A2 - L2 %*% theta[(p1+p2+1):(p1+p2+p3)])^2,
        sw*(Y - (theta[p-4] + theta[p-3]*A1star + theta[p-2]*A2 +
                   theta[p-1]*A3star)),
        sw*(Y - (theta[p-4] + theta[p-3]*A1star + theta[p-2]*A2 +
                   theta[p-1]*A3star))*A1star,
        sw*(Y - (theta[p-4] + theta[p-3]*A1star + theta[p-2]*A2 +
                   theta[p-1]*A3star))*A2,
        sw*(Y - (theta[p-4] + theta[p-3]*A1star + theta[p-2]*A2 +
                   theta[p-1]*A3star))*A3star,
        sw*(theta[p] - (Y - (theta[p-4] + theta[p-3]*A1star +
                               theta[p-2]*A2 + theta[p-1]*A3star))^2)
      )
    }
  }

  results_wr <- m_estimate(estFUN = eefun_wr, data = data,
                           outer_args = list(denom_mod1, num_mod1,
                                             denom_mod2, num_mod2),
                           compute_roots = FALSE,
                           roots = c(coef(denom_mod1), coef(num_mod1),
                                     coef(denom_mod2), coef(num_mod2),
                                     sigma(denom_mod1)^2,
                                     sigma(denom_mod2)^2,
                                     coef(wmod), sigma(wmod)^2))


  eefun_wr <- function(data) {
    Y <- data$Y
    A1star <- data$A1star
    A2 <- data$A2
    A3star <- data$A3star
    sw <- data$sw
    n <- dim(data)[1]

    function(theta) {
      c(sw*(Y - (theta[1] + theta[2]*A1star + theta[3]*A2 +
                   theta[4]*A3star)),
        sw*(Y - (theta[1] + theta[2]*A1star + theta[3]*A2 +
                   theta[4]*A3star))*A1star,
        sw*(Y - (theta[1] + theta[2]*A1star + theta[3]*A2 +
                   theta[4]*A3star))*A2,
        sw*(Y - (theta[1] + theta[2]*A1star + theta[3]*A2 +
                   theta[4]*A3star))*A3star,
        sw*(theta[5] - (Y - (theta[1] + theta[2]*A1star +
                               theta[3]*A2 + theta[4]*A3star))^2)
      )
    }
  }

  results_wr <- m_estimate(estFUN = eefun_wr, data = data,
                           compute_roots = FALSE,
                           roots = c(coef(wmod), sigma(wmod)^2))

  eefun_wr <- function(data) {
    Y <- data$Y
    A3star <- data$A3star
    A1star <- data$A1star
    L <- data$L
    A2 <- model.response(model.frame(model3, data = data))
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

      dens_denom1 <- dnorm(A1star, rho1, sqrt(theta[p-6]))
      dens_num1 <- dnorm(A1star, rho2, sqrt(theta[p-6]))
      dens_denom2 <- dnorm(A2, rho3, sqrt(theta[p-5]))
      dens_num2 <- dnorm(A2, rho4, sqrt(theta[p-5]))
      sw1 <- dens_num1 / dens_denom1
      sw2 <- dens_num2 / dens_denom2
      sw <- sw1*sw2

      score_eqns1 <- apply(L1, 2, function(x) sum((A1star - rho1) * x))
      score_eqns2 <- apply(V1, 2, function(x) sum((A1star - rho2) * x))
      score_eqns3 <- apply(L2, 2, function(x) sum((A2 - rho3) * x))
      score_eqns4 <- apply(V2, 2, function(x) sum((A2 - rho4) * x))

      c((A1star - rho1)*L1[,1],
        (A1star - rho1)*L1[,2],
        (A1star - rho2)*V1[,1],
        (A2 - rho3)*L2[,1],
        (A2 - rho3)*L2[,2],
        (A2 - rho4)*V2[,1],
        (n-p1)/n * theta[p-6] - (A1star - theta[1] - theta[2]*L)^2,
        (n-p3)/n * theta[p-5] - (A2 - L2 %*% theta[(p1+p2+1):(p1+p2+p3)])^2,
        sw*(Y - (theta[p-4] + theta[p-3]*A1star + theta[p-2]*A2 +
                   theta[p-1]*A3star)),
        sw*(Y - (theta[p-4] + theta[p-3]*A1star + theta[p-2]*A2 +
                   theta[p-1]*A3star))*A1star,
        sw*(Y - (theta[p-4] + theta[p-3]*A1star + theta[p-2]*A2 +
                   theta[p-1]*A3star))*A2,
        sw*(Y - (theta[p-4] + theta[p-3]*A1star + theta[p-2]*A2 +
                   theta[p-1]*A3star))*A3star,
        sw*(theta[p] - (Y - (theta[p-4] + theta[p-3]*A1star +
                               theta[p-2]*A2 + theta[p-1]*A3star))^2)
      )
    }
  }

  results_wr <- m_estimate(estFUN = eefun_wr, data = data,
                           outer_args = list(denom_mod1, num_mod1,
                                             denom_mod2, num_mod2),
                           compute_roots = FALSE,
                           roots = c(coef(denom_mod1), coef(num_mod1),
                                     coef(denom_mod2), coef(num_mod2),
                                     sigma(denom_mod1)^2,
                                     sigma(denom_mod2)^2,
                                     coef(wmod), sigma(wmod)^2))

  }



  se1_ipw <- sqrt(vcov(wmod)[2, 2])
  coverage1_ipw <- 1*(coef(wmod)[2] - 1.96*se1_ipw < (beta1) &
                      coef(wmod)[2] + 1.96*se1_ipw > (beta1))

  se2_ipw <- sqrt(vcov(wmod)[3, 3])
  coverage2_ipw <- 1*(coef(wmod)[3] - 1.96*se2_ipw < (beta2) &
                      coef(wmod)[3] + 1.96*se2_ipw > (beta2))

  se3_ipw <- sqrt(vcov(wmod)[4, 4])
  coverage3_ipw <- 1*(coef(wmod)[4] - 1.96*se3_ipw < beta3 &
                      coef(wmod)[4] + 1.96*se3_ipw > beta3)

  # Model 4: IPW-CSME
  # Get point estimates and variance using geex
  # We get point estimates using a geex with only some parameters
  # Then variance using geex with the full set of parameters
  # We don't get point estimates with the second geex because of
  # high level of divergence
  eefun_ipw_csme1 <- function(data) {
    Y <- data$Y
    A1star <- data$A1star
    A2 <- data$A2
    A3star <- data$A3star
    sw <- data$sw
    delta1 <- function(beta1) {
      A1star + beta1*sigma_me1*Y
    }
    delta3 <- function(beta3) {
      A3star + beta3*sigma_me3*Y
    }
    H <- function(x) {
      1 / (1 + exp(-x))
    }
    condexp <- function(beta0, beta1, beta2, beta3) {
      H(beta0 + beta1*delta1(beta1) + beta2*A2 + beta3*delta3(beta3) -
        (beta1^2 * sigma_me1 + beta3^2 * sigma_me3) / 2)
    }
    function(theta) {
      c(sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4])),
        sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4]))*
          delta1(theta[2]),
        sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4]))*
          A2,
        sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4]))*
          delta3(theta[4])
      )
    }
  }

  failed <- TRUE
  j <- 1

  while(failed == TRUE & j < 6) {

    failed <- FALSE
    startvec <- c(coef(wmod)[1:4])*(j == 1) +
      c(coef(mod)[1:4])*(j == 2) +
      c(rnorm(4, 0, j/5))*(j > 2)
    results_ipw_csme1 <-
      tryCatch(m_estimate(estFUN = eefun_ipw_csme1, data = data,
                          root_control =
                            setup_root_control(start = startvec)),
               error = function(e) { failed <<- TRUE})
    if (failed == FALSE) {
      if (abs(coef(results_ipw_csme1)[2]) > 3) { failed <- TRUE }
    }
    j <- j + 1

  }

  if(failed == TRUE) { next }

  bias1_ipw_csme <- coef(results_ipw_csme1)[2] - (beta1)
  bias2_ipw_csme <- coef(results_ipw_csme1)[3] - (beta2)
  bias3_ipw_csme <- coef(results_ipw_csme1)[4] - beta3

  eefun_ipw_csme2 <- function(data, model1, model2, model3, model4) {
    Y <- data$Y
    A3star <- data$A3star
    A1star <- model.response(model.frame(model1, data = data))
    L1 <- model.matrix(model1, data = data)
    V1 <- model.matrix(model2, data = data)
    A2 <- model.response(model.frame(model3, data = data))
    L2 <- model.matrix(model3, data = data)
    V2 <- model.matrix(model4, data = data)
    #sw <- data$sw
    n <- dim(data)[1]
    delta1 <- function(beta1, sigma_ep) {
      A1star + beta1*sigma_me1*Y
    }
    delta3 <- function(beta3, sigma_ep) {
      A3star + beta3*sigma_me3*Y
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

  run <- F
  if (run == T) {
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
  }



  se1_ipw_csme <- sqrt(vcov(results_ipw_csme1)[2, 2])
  coverage1_ipw_csme <-
    1*(coef(results_ipw_csme1)[2] - 1.96*se1_ipw_csme < (beta1) &
         coef(results_ipw_csme1)[2] + 1.96*se1_ipw_csme > (beta1))

  se2_ipw_csme <- sqrt(vcov(results_ipw_csme1)[3, 3])
  coverage2_ipw_csme <-
    1*(coef(results_ipw_csme1)[3] - 1.96*se2_ipw_csme < (beta2) &
         coef(results_ipw_csme1)[3] + 1.96*se2_ipw_csme > (beta2))

  se3_ipw_csme <- sqrt(vcov(results_ipw_csme1)[4, 4])
  coverage3_ipw_csme <-
    1*(coef(results_ipw_csme1)[4] - 1.96*se3_ipw_csme < beta3 &
         coef(results_ipw_csme1)[4] + 1.96*se3_ipw_csme > beta3)

  return(c(bias1_lm, bias1_csme, bias1_ipw, bias1_ipw_csme,
           bias2_lm, bias2_csme, bias2_ipw, bias2_ipw_csme,
           bias3_lm, bias3_csme, bias3_ipw, bias3_ipw_csme,
           se1_lm, se1_csme, se1_ipw, se1_ipw_csme,
           se2_lm, se2_csme, se2_ipw, se2_ipw_csme,
           se3_lm, se3_csme, se3_ipw, se3_ipw_csme,
           coverage1_lm, coverage1_csme, coverage1_ipw, coverage1_ipw_csme,
           coverage2_lm, coverage2_csme, coverage2_ipw, coverage2_ipw_csme,
           coverage3_lm, coverage3_csme, coverage3_ipw, coverage3_ipw_csme))

}

n <- 800
low_n <- 400
trials <- seq(1, nsims)
combos <- data.frame(trials = rep(trials, length(beta1)),
                     betas = rep(beta1, each = nsims),
                     ns = rep(low_n, each = nsims))
i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
combo_i <- combos[(i), ]

set.seed(i*1000)
sim <- with(combo_i, mapply(simulator, trials, betas, ns))

# Output
outfile <- paste("./Results/results_scen2_", i, ".Rdata", sep = "")
save(sim, file = outfile)
