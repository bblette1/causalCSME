# Load libraries
rm(list = ls())
library(dplyr)
library(geex)
library(ggplot2)

# Build application data set
load("C:/Users/bblette1/Downloads/HVTN505_2019-08-08/HVTN505/data/dat.505.rda")

assays <- subset(var.505, assay %in% c("fcrR2a", "fcrR3a", "phago"))

primarydat <- read.csv("C:/Users/bblette1/Documents/primary505_for_sharing .csv")
primarydat$ptid <- primarydat$pub_id

fulldat <- merge(dat.505, primarydat, by = "ptid", all = T)
fulldat$marker1 <- log(fulldat$ADCP1)
fulldat$marker2 <- fulldat$R2aConSgp140CFI

# Treat behavior risk as binary variable
fulldat$bhv_bin.y <- 1*(fulldat$bhvrisk.y > 0)

# Simple imputation of missing confounders, marker2 0 value
set.seed(1234)

bmimod <- lm(BMI.y ~ HIVall.y + age.y + bhvrisk.y + race.y, data = fulldat)
bmi_imp_dat <- data.frame(HIVall.y = fulldat$HIVall.y[is.na(fulldat$BMI.y)],
                          age.y = fulldat$age.y[is.na(fulldat$BMI.y)],
                          bhvrisk.y = fulldat$bhvrisk.y[is.na(fulldat$BMI.y)],
                          race.y = fulldat$race.y[is.na(fulldat$BMI.y)])
fulldat$BMI.y[is.na(fulldat$BMI.y)] <- predict(bmimod, bmi_imp_dat)

bhvmod <- glm(bhv_bin.y ~ HIVall.y + age.y + BMI.y + race.y, data = fulldat,
              family = binomial)
bhv_imp_dat <-
  data.frame(HIVall.y = fulldat$HIVall.y[is.na(fulldat$bhv_bin.y)],
             age.y = fulldat$age.y[is.na(fulldat$bhv_bin.y)],
             BMI.y = fulldat$BMI.y[is.na(fulldat$bhv_bin.y)],
             race.y = fulldat$race.y[is.na(fulldat$bhv_bin.y)])
fulldat$bhv_bin.y[is.na(fulldat$bhv_bin.y)] <-
  1*(predict(bhvmod, bhv_imp_dat, type = "response") > 0.5)

imputemax <- min(fulldat$marker2[fulldat$marker2 > 0], na.rm = T)
fulldat$marker2[fulldat$marker2 == 0 & !is.na(fulldat$marker2) &
                  fulldat$trt.y == 1] <- runif(1, 0, imputemax)

fulldat$sampled <- !is.na(fulldat$marker1)
sampmod <- glm(sampled ~ HIVwk28preunbl.y,
               data = fulldat[fulldat$trt.y == 1, ])

# Create analysis dataset
# Restrict to treatment, use all three markers
fulldat$CD4PFS <- fulldat$CD4_ANYVRCENV_PolyfunctionalityScore_score_bin
fulldat$CD8PFS <- fulldat$CD8_ANYVRCENV_PolyfunctionalityScore_score_bin
analysisdat <- fulldat %>%
  filter(trt.y == 1) %>%
  select(HIVall.y, HIVwk28preunbl.y, marker1, marker2,
         race.y, BMI.y, age.y, bhv_bin.y, ADCP1_bin, CD4PFS, CD8PFS)

# Probability selected by treatment
# Weights based on HIV at week 28, not end of study
analysisdat$estweights <-
  (!is.na(analysisdat$marker1))*(analysisdat$HIVwk28preunbl.y == 1) /
  (25 / sum(fulldat$HIVwk28preunbl.y[fulldat$trt.y == 1])) +
  (!is.na(analysisdat$marker1))*(analysisdat$HIVwk28preunbl.y == 0) /
  (125 / sum(1 - fulldat$HIVwk28preunbl.y[fulldat$trt.y == 1]))



#########################################################################
# Analysis
#########################################################################

# Alter the measured 0 value for one individual
analysisdat$marker1[analysisdat$marker1 == 0 &
                      !is.na(analysisdat$marker1)] <- 0.01

# Set NA marker values to 0
analysisdat$marker1[is.na(analysisdat$marker1)] <- 0
analysisdat$marker2[is.na(analysisdat$marker2)] <- 0

# Estimate IPTW
# First marker
denom_mod1 <- lm(marker1 ~ race.y + BMI.y + age.y + bhv_bin.y +
                   CD4PFS + CD8PFS,
                 data = analysisdat, weights = estweights)
p_denom1 <- predict(denom_mod1, type='response')
dens_denom1 <-
  dnorm(analysisdat$marker1[analysisdat$estweights > 0],
        p_denom1,
        sd(residuals(denom_mod1)))

num_mod1 <- lm(marker1 ~ 1, data = analysisdat, weights = estweights)
p_num1 <- predict(num_mod1, type='response')
dens_num1 <-
  dnorm(analysisdat$marker1[analysisdat$estweights > 0],
        p_num1[analysisdat$estweights > 0],
        sd(residuals(num_mod1)[analysisdat$estweights > 0]))

analysisdat$w1 <- rep(0, dim(analysisdat)[1])
analysisdat$w1[analysisdat$estweights > 0] <- dens_num1 / dens_denom1
analysisdat$w1[is.na(analysisdat$w1)] <- 0

analysisdat$sw1 <- analysisdat$w1*analysisdat$estweights

# Second marker
denom_mod2 <- lm(marker2 ~ race.y + BMI.y + age.y + bhv_bin.y +
                   CD4PFS + CD8PFS,
                 data = analysisdat, weights = estweights)
p_denom2 <- predict(denom_mod2, type='response')
dens_denom2 <-
  dnorm(analysisdat$marker2[analysisdat$estweights > 0],
        p_denom2,
        sd(residuals(denom_mod2)))

num_mod2 <- lm(marker2 ~ 1, data = analysisdat, weights = estweights)
p_num2 <- predict(num_mod2, type='response')
dens_num2 <-
  dnorm(analysisdat$marker2[analysisdat$estweights > 0],
        p_num2[analysisdat$estweights > 0],
        sd(residuals(num_mod2)[analysisdat$estweights > 0]))

analysisdat$w2 <- rep(0, dim(analysisdat)[1])
analysisdat$w2[analysisdat$estweights > 0] <- dens_num2 / dens_denom2
analysisdat$w2[is.na(analysisdat$w2)] <- 0

analysisdat$sw2 <- analysisdat$w2*analysisdat$estweights

# Starting values
glmmod1 <-
  lm(HIVwk28preunbl.y ~ marker1 + bhv_bin.y + age.y + race.y + BMI.y +
        CD8PFS + CD4PFS + marker1*bhv_bin.y + marker1*age.y,
      data = analysisdat, weights = sw1)

glmmod2 <-
  lm(HIVwk28preunbl.y ~ marker2 + bhv_bin.y + age.y + race.y + BMI.y +
        CD8PFS + CD4PFS + marker2*bhv_bin.y + marker2*age.y,
      data = analysisdat, weights = sw2)

# DR estimating equations
eefun_dr1 <- function(data, val) {
  Y <- data$HIVwk28preunbl.y
  A1star <- data$marker1
  A1star[is.na(A1star)] <- 0
  L1 <- data$bhv_bin.y
  L2 <- data$age.y
  L3 <- data$race.y
  L4 <- data$BMI.y
  L5 <- data$CD8PFS
  L5[is.na(L5)] <- 0
  L6 <- data$CD4PFS
  L6[is.na(L6)] <- 0
  sw <- data$sw1
  delta1 <- function(beta1, beta8, beta9, sigma_ep) {
    A1star + sigma_me1*(beta1 + beta8*L1 + beta9*L2)*Y / sigma_ep
  }
  condexp <- function(beta0, beta1, beta2, beta3, beta4, beta5, beta6,
                      beta7, beta8, beta9, sigma_ep) {
    (beta0 + beta2*L1 + beta3*L2 + beta4*L3 + beta5*L4 +
       beta6*L5 + beta7*L6 +
       (beta1 + beta8*L1 + beta9*L2)*delta1(beta1, beta8, beta9, sigma_ep)) /
      (1 + ((beta1 + beta8*L1 + beta9*L2)^2)*sigma_me1 / sigma_ep)[[1]]
  }
  condvar <- function(beta1, beta8, beta9, sigma_ep) {
    sigma_ep / (1 + ((beta1 + beta8*L1 + beta9*L2)^2)*sigma_me1 /
                  sigma_ep)[[1]]
  }
  function(theta) {
    c(sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9], theta[10],
                      theta[11])),
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9], theta[10],
                      theta[11]))*
        delta1(theta[2], theta[9], theta[10], theta[11]),
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9], theta[10],
                      theta[11]))*
        L1,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9], theta[10],
                      theta[11]))*
        L2,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9], theta[10],
                      theta[11]))*
        L3,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9], theta[10],
                      theta[11]))*
        L4,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9], theta[10],
                      theta[11]))*
        L5,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9], theta[10],
                      theta[11]))*
        L6,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9], theta[10],
                      theta[11]))*
        delta1(theta[2], theta[9], theta[10], theta[11])*L1,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9], theta[10],
                      theta[11]))*
        delta1(theta[2], theta[9], theta[10], theta[11])*L2,
      sw*(theta[11] - theta[11]*
            (Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                         theta[6], theta[7], theta[8], theta[9], theta[10],
                         theta[11]))^2 /
            condvar(theta[2], theta[9], theta[10], theta[11])),
      theta[12] - (theta[1] + theta[2]*val + theta[3]*L1 + theta[4]*L2 +
              theta[5]*L3 + theta[6]*L4 + theta[7]*L5 + theta[8]*L6 +
              theta[9]*val*L1 + theta[10]*val*L2)
    )
  }
}

eefun_dr2 <- function(data, val) {
  Y <- data$HIVwk28preunbl.y
  A2star <- data$marker2
  A2star[is.na(A2star)] <- 0
  L1 <- data$bhv_bin.y
  L2 <- data$age.y
  L3 <- data$race.y
  L4 <- data$BMI.y
  L5 <- data$CD8PFS
  L5[is.na(L5)] <- 0
  L6 <- data$CD4PFS
  L6[is.na(L6)] <- 0
  sw <- data$sw2
  delta2 <- function(beta1, beta8, beta9, sigma_ep) {
    A2star + sigma_me2*(beta1 + beta8*L1 + beta9*L2)*Y / sigma_ep
  }
  condexp <- function(beta0, beta1, beta2, beta3, beta4, beta5, beta6,
                      beta7, beta8, beta9, sigma_ep) {
    (beta0 + beta2*L1 + beta3*L2 + beta4*L3 + beta5*L4 +
       beta6*L5 + beta7*L6 +
       (beta1 + beta8*L1 + beta9*L2)*delta2(beta1, beta8, beta9, sigma_ep)) /
      (1 + ((beta1 + beta8*L1 + beta9*L2)^2)*sigma_me2 / sigma_ep)[[1]]
  }
  condvar <- function(beta1, beta8, beta9, sigma_ep) {
    sigma_ep / (1 + ((beta1 + beta8*L1 + beta9*L2)^2)*sigma_me2 /
                  sigma_ep)[[1]]
  }
  function(theta) {
    c(sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9], theta[10],
                      theta[11])),
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9], theta[10],
                      theta[11]))*
        delta2(theta[2], theta[9], theta[10], theta[11]),
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9], theta[10],
                      theta[11]))*
        L1,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9], theta[10],
                      theta[11]))*
        L2,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9], theta[10],
                      theta[11]))*
        L3,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9], theta[10],
                      theta[11]))*
        L4,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9], theta[10],
                      theta[11]))*
        L5,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9], theta[10],
                      theta[11]))*
        L6,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9], theta[10],
                      theta[11]))*
        delta2(theta[2], theta[9], theta[10], theta[11])*L1,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9], theta[10],
                      theta[11]))*
        delta2(theta[2], theta[9], theta[10], theta[11])*L2,
      sw*(theta[11] - theta[11]*
            (Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                         theta[6], theta[7], theta[8], theta[9], theta[10],
                         theta[11]))^2 /
            condvar(theta[2], theta[9], theta[10], theta[11])),
      theta[12] - (theta[1] + theta[2]*val + theta[3]*L1 + theta[4]*L2 +
                     theta[5]*L3 + theta[6]*L4 + theta[7]*L5 + theta[8]*L6 +
                     theta[9]*val*L1 + theta[10]*val*L2)
    )
  }
}


# Calculate DR estimator for 4 ME and a grid of values
me_list <- c(0, 1/6, 1/4, 1/3)
m1vals <- seq(0.5, 3, 0.1)
m2vals <- seq(7, 11, 0.2)

dr_ests1 <- array(NA, dim = c(4, length(m1vals)))
dr_ests2 <- array(NA, dim = c(4, length(m2vals)))
dr_se1 <- array(NA, dim = c(4, length(m1vals)))
dr_se2 <- array(NA, dim = c(4, length(m2vals)))

for (me in 1:4) {

  sigma_me1 <- var(analysisdat$marker1[analysisdat$marker1 > 0.5])*me_list[me]
  sigma_me2 <- var(analysisdat$marker2[analysisdat$marker2 > 7])*me_list[me]

  for (i in 1:length(m1vals)) {

    startvec <- c(coef(glmmod1), sigma(glmmod1)^2, 0.05)
    if (i > 1) {
      startvec <- coef(results_dr1)
    }

    results_dr1 <-
      m_estimate(estFUN = eefun_dr1, data = analysisdat,
                 outer_args = list(m1vals[i]), compute_roots = TRUE,
                 root_control =
                   setup_root_control(start = startvec))
    dr_ests1[me, i] <- coef(results_dr1)[12]
    dr_se1[me, i] <- sqrt(vcov(results_dr1)[12, 12])

  }

  for (i in 1:length(m2vals)) {

    startvec <- c(coef(glmmod2), sigma(glmmod2)^2, 0.1)
    if (i > 1) {
      startvec <- coef(results_dr2)
    }

    results_dr2 <-
      m_estimate(estFUN = eefun_dr2, data = analysisdat,
                 outer_args = list(m2vals[i]), compute_roots = TRUE,
                 root_control =
                   setup_root_control(start = startvec))
    dr_ests2[me, i] <- coef(results_dr2)[12]
    dr_se2[me, i] <- sqrt(vcov(results_dr2)[12, 12])

  }

}

# Make plots
par(mfrow = c(2, 4))

# First row is first marker, four plots for four ME
plot(seq(0, 3, 0.1), dr_ests1[1, ], ylim = c(0, 0.35),
     main = expression(paste("No ME: ", sigma^2, "=0")), xlab = "ADCP",
     ylab = "HIV risk at study end", type = "n")
lines(seq(0, 3, 0.1), dr_ests1[1, ], lwd = 2)
polygon(x = c(seq(0, 3, 0.1), rev(seq(0, 3, 0.1))),
        y = c(dr_ests1[1, ] - 1.96*dr_se1[1, ],
              rev(dr_ests1[1, ] + 1.96*dr_se1[1, ])),
        col =  adjustcolor("gray", alpha.f = 0.4), border = NA)

plot(seq(0, 3, 0.1), dr_ests1[2, ], ylim = c(0, 0.35),
     main = expression(paste("Mild ME: ", sigma^2, "=0.1")), xlab = "ADCP",
     ylab = "HIV risk at study end", type = "n")
lines(seq(0, 3, 0.1), dr_ests1[2, ], lwd = 2)
polygon(x = c(seq(0, 3, 0.1), rev(seq(0, 3, 0.1))),
        y = c(dr_ests1[2, ] - 1.96*dr_se1[2, ],
              rev(dr_ests1[2, ] + 1.96*dr_se1[2, ])),
        col =  adjustcolor("gray", alpha.f = 0.4), border = NA)

plot(seq(0, 3, 0.1), dr_ests1[3, ], ylim = c(0, 0.35),
     main = expression(paste("Moderate ME: ", sigma^2, "=0.2")),
     xlab = "ADCP", ylab = "HIV risk at study end", type = "n")
lines(seq(0, 3, 0.1), dr_ests1[3, ], lwd = 2)
polygon(x = c(seq(0, 3, 0.1), rev(seq(0, 3, 0.1))),
        y = c(dr_ests1[3, ] - 1.96*dr_se1[3, ],
              rev(dr_ests1[3, ] + 1.96*dr_se1[3, ])),
        col =  adjustcolor("gray", alpha.f = 0.4), border = NA)

plot(seq(0, 3, 0.1), dr_ests1[4, ], ylim = c(0, 0.35),
     main = expression(paste("High ME: ", sigma^2, "=0.3")), xlab = "ADCP",
     ylab = "HIV risk at study end", type = "n")
lines(seq(0, 3, 0.1), dr_ests1[4, ], lwd = 2)
polygon(x = c(seq(0, 3, 0.1), rev(seq(0, 3, 0.1))),
        y = c(dr_ests1[4, ] - 1.96*dr_se1[4, ],
              rev(dr_ests1[4, ] + 1.96*dr_se1[4, ])),
        col =  adjustcolor("gray", alpha.f = 0.4), border = NA)

# Second row is second marker, four plots for four ME
plot(seq(0, 11, 0.5), dr_ests2[1, ], ylim = c(0, 1),
     main = expression(paste("No ME: ", sigma^2, "=0")), xlab = "RII",
     ylab = "HIV risk at study end", type = "n")
lines(seq(0, 11, 0.5), dr_ests2[1, ], lwd = 2)
polygon(x = c(seq(0, 11, 0.5), rev(seq(0, 11, 0.5))),
        y = c(dr_ests2[1, ] - 1.96*dr_se2[1, ],
              rev(dr_ests2[1, ] + 1.96*dr_se2[1, ])),
        col =  adjustcolor("gray", alpha.f = 0.4), border = NA)

plot(seq(0, 11, 0.5), dr_ests2[2, ], ylim = c(0, 1),
     main = expression(paste("Mild ME: ", sigma^2, "=0.1")), xlab = "RII",
     ylab = "HIV risk at study end", type = "n")
lines(seq(0, 11, 0.5), dr_ests2[2, ], lwd = 2)
polygon(x = c(seq(0, 11, 0.5), rev(seq(0, 11, 0.5))),
        y = c(dr_ests2[2, ] - 1.96*dr_se2[2, ],
              rev(dr_ests2[2, ] + 1.96*dr_se2[2, ])),
        col =  adjustcolor("gray", alpha.f = 0.4), border = NA)

plot(seq(0, 11, 0.5), dr_ests2[3, ], ylim = c(0, 1),
     main = expression(paste("Moderate ME: ", sigma^2, "=0.2")),
     xlab = "RII", ylab = "HIV risk at study end", type = "n")
lines(seq(0, 11, 0.5), dr_ests2[3, ], lwd = 2)
polygon(x = c(seq(0, 11, 0.5), rev(seq(0, 11, 0.5))),
        y = c(dr_ests2[3, ] - 1.96*dr_se2[3, ],
              rev(dr_ests2[3, ] + 1.96*dr_se2[3, ])),
        col =  adjustcolor("gray", alpha.f = 0.4), border = NA)

plot(seq(0, 11, 0.5), dr_ests2[4, ], ylim = c(0, 1),
     main = expression(paste("High ME: ", sigma^2, "=0.3")), xlab = "RII",
     ylab = "HIV risk at study end", type = "n")
lines(seq(0, 11, 0.5), dr_ests2[4, ], lwd = 2)
polygon(x = c(seq(0, 11, 0.5), rev(seq(0, 11, 0.5))),
        y = c(dr_ests2[4, ] - 1.96*dr_se2[4, ],
              rev(dr_ests2[4, ] + 1.96*dr_se2[4, ])),
        col =  adjustcolor("gray", alpha.f = 0.4), border = NA)

# Alternative plot, more of a lattice structure
melist <- c(expression(paste("No ME: ", sigma^2, "=0")),
            expression(paste("Low ME: ", sigma[me], "=", sigma^2, "/6")),
            expression(paste("Moderate ME: ", sigma^2, "=0.2")),
            expression(paste("High ME: ", sigma^2, "=0.3")))

melist <- c("No ME: sigma^2 = 0", "Low ME: sigma^2 = 0.1",
            "Moderate ME: sigma^2 = 0.2", "High ME: sigma^2 = 0.3")

latdat <- data.frame(vals = c(rep(seq(0.5, 3, 0.1), 4),
                              rep(seq(7, 11, 0.2), 4)),
                     Risk = c(dr_ests1[1, ], dr_ests1[2, ],
                              dr_ests1[3, ], dr_ests1[4, ],
                              dr_ests2[1, ], dr_ests2[2, ],
                              dr_ests2[3, ], dr_ests2[4, ]),
                     Risk_low = c(dr_ests1[1, ] - 1.96*dr_se1[1, ],
                                  dr_ests1[2, ] - 1.96*dr_se1[2, ],
                                  dr_ests1[3, ] - 1.96*dr_se1[3, ],
                                  dr_ests1[4, ] - 1.96*dr_se1[4, ],
                                  dr_ests2[1, ] - 1.96*dr_se2[1, ],
                                  dr_ests2[2, ] - 1.96*dr_se2[2, ],
                                  dr_ests2[3, ] - 1.96*dr_se2[3, ],
                                  dr_ests2[4, ] - 1.96*dr_se2[4, ]),
                     Risk_upp = c(dr_ests1[1, ] + 1.96*dr_se1[1, ],
                                  dr_ests1[2, ] + 1.96*dr_se1[2, ],
                                  dr_ests1[3, ] + 1.96*dr_se1[3, ],
                                  dr_ests1[4, ] + 1.96*dr_se1[4, ],
                                  dr_ests2[1, ] + 1.96*dr_se2[1, ],
                                  dr_ests2[2, ] + 1.96*dr_se2[2, ],
                                  dr_ests2[3, ] + 1.96*dr_se2[3, ],
                                  dr_ests2[4, ] + 1.96*dr_se2[4, ]),
                     ME = c(rep(round(me_list, 3), each = 26),
                            rep(round(me_list, 3), each = 21)),
                     Exposure = c(rep("ADCP", 104), rep("RII", 84)))

latdat$Risk[latdat$Risk < 0] <- 0
latdat$Risk_low[latdat$Risk_low < 0] <- 0

# New facet label names for ME variable
me.labs <- c(expression(paste("No ME: ", sigma^2, "=0")),
             expression(paste("Low ME: ", sigma^2, "=0.1")),
             expression(paste("Moderate ME: ", sigma^2, "=0.2")),
             expression(paste("High ME: ", sigma^2, "=0.3")))
names(me.labs) <- c("0", "0.1", "0.2", "0.3")

# New facet label names for Exposure variable
exp.labs <- c("ADCP", expression(paste("R", as.character(as.roman(2)))))
names(exp.labs) <- c("ADCP", "RII")

ggplot(latdat, aes(x = vals, y = Risk)) +
  geom_line() +
  facet_grid(ME ~ Exposure, scales = "free",
             labeller = label_bquote(sigma[me]^2 == .(ME) * sigma^2)) +
  #labeller = labeller(ME = me.labs, Exposure = exp.labs)) +
  geom_ribbon(aes(ymin = Risk_low, ymax = Risk_upp), alpha = 0.3) +
  xlab("Exposure values") + ylab("HIV risk at study end") +
  ylim(c(0, 0.25)) +
  theme_bw()
