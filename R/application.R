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

# Simple imputation of missing confounder, marker2 0 value
set.seed(1234)
fulldat$BMI.y[is.na(fulldat$BMI.y)] <- median(fulldat$BMI.y, na.rm = T)
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
         race.y, BMI.y, age.y, bhvrisk.y, ADCP1_bin, CD4PFS, CD8PFS)

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
denom_mod1 <- lm(marker1 ~ race.y + BMI.y + age.y + bhvrisk.y,
                 data = analysisdat, weights = estweights)
p_denom1 <- predict(denom_mod1, type='response')
dens_denom1 <-
  dnorm(analysisdat$marker1, p_denom1,
        sd(residuals(denom_mod1)[analysisdat$marker1 != 0 |
                                   analysisdat$HIVwk28preunbl.y == 1]))
num_mod1 <- lm(marker1 ~ 1, data = analysisdat, weights = estweights)
p_num1 <- predict(num_mod1, type='response')
dens_num1 <-
  dnorm(analysisdat$marker1, p_num1,
        sd(residuals(denom_mod1)[analysisdat$marker1 != 0 |
                                   analysisdat$HIVwk28preunbl.y == 1]))
analysisdat$w1 <- rep(0, dim(analysisdat)[1])
analysisdat$w1[analysisdat$marker1 != 0 |
                 analysisdat$HIVwk28preunbl.y == 1] <-
  dens_num1[analysisdat$marker1 != 0 |
              analysisdat$HIVwk28preunbl.y == 1] /
  dens_denom1[analysisdat$marker1 != 0 |
                analysisdat$HIVwk28preunbl.y == 1]
analysisdat$w1[is.na(analysisdat$w1)] <- 0

analysisdat$sw1 <- analysisdat$w1*analysisdat$estweights

# Second marker
denom_mod2 <- lm(marker2 ~ race.y + BMI.y + age.y + bhvrisk.y,
                 data = analysisdat, weights = estweights)
p_denom2 <- predict(denom_mod2, type='response')
dens_denom2 <-
  dnorm(analysisdat$marker2, p_denom2,
        sd(residuals(denom_mod2)[analysisdat$marker2 != 0 |
                                   analysisdat$HIVwk28preunbl.y == 1]))
num_mod2 <- lm(marker2 ~ 1, data = analysisdat, weights = estweights)
p_num2 <- predict(num_mod2, type='response')
dens_num2 <-
  dnorm(analysisdat$marker2, p_num2,
        sd(residuals(denom_mod2)[analysisdat$marker2 != 0 |
                                   analysisdat$HIVwk28preunbl.y == 1]))
analysisdat$w2 <- rep(0, dim(analysisdat)[1])
analysisdat$w2[analysisdat$marker2 != 0 |
                 analysisdat$HIVwk28preunbl.y == 1] <-
  dens_num2[analysisdat$marker2 != 0 |
              analysisdat$HIVwk28preunbl.y == 1] /
  dens_denom2[analysisdat$marker2 != 0 |
                analysisdat$HIVwk28preunbl.y == 1]
analysisdat$w2[is.na(analysisdat$w2)] <- 0

analysisdat$sw2 <- analysisdat$w2*analysisdat$estweights

# Starting values
glmmod1 <-
  glm(HIVwk28preunbl.y ~ marker1 + race.y + age.y + CD8PFS + CD4PFS +
        marker1*race.y + marker1*age.y,
      data = analysisdat, family = "binomial", weights = sw1)

glmmod2 <-
  glm(HIVwk28preunbl.y ~ marker2 + race.y + age.y + CD8PFS + CD4PFS +
        marker2*race.y + marker2*age.y,
      data = analysisdat, family = "binomial", weights = sw2)

# DR estimating equations
eefun_dr1 <- function(data, val) {
  Y <- data$HIVwk28preunbl.y
  A1star <- data$marker1
  A1star[is.na(A1star)] <- 0
  L1 <- data$race.y
  L2 <- data$age.y
  L3 <- data$CD8PFS
  L3[is.na(L3)] <- 0
  L4 <- data$CD4PFS
  L4[is.na(L4)] <- 0
  sw <- data$sw1
  delta1 <- function(beta1, beta6, beta7) {
    A1star + sigma_me1*(beta1 + beta6*L1 + beta7*L2)*Y
  }
  H <- function(x) {
    1 / (1 + exp(-x))
  }
  condexp <- function(beta0, beta1, beta2, beta3, beta4, beta5, beta6,
                      beta7) {
    H(beta0 + beta2*L1 + beta3*L2 + beta4*L3 + beta5*L4 +
        (beta1 + beta6*L1 + beta7*L2)*
        delta1(beta1, beta6, beta7) -
        ((beta1 + beta6*L1 + beta7*L2)^2 * sigma_me1) / 2)
  }
  function(theta) {
    c(sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8])),
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8]))*
        delta1(theta[2], theta[7], theta[8]),
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8]))*
        L1,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8]))*
        L2,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8]))*
        L3,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8]))*
        L4,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8]))*
        delta1(theta[2], theta[7], theta[8])*L1,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8]))*
        delta1(theta[2], theta[7], theta[8])*L2,
      theta[9] -
        exp(theta[1] + theta[2]*val + theta[3]*L1 + theta[4]*L2 +
            theta[5]*L3 + theta[6]*L4 + theta[7]*val*L1 + theta[8]*val*L2) /
        (1 + exp(theta[1] + theta[2]*val + theta[3]*L1 + theta[4]*L2 +
            theta[5]*L3 + theta[6]*L4 + theta[7]*val*L1 + theta[8]*val*L2))
    )
  }
}

eefun_dr2 <- function(data, val) {
  Y <- data$HIVwk28preunbl.y
  A2star <- data$marker2
  A2star[is.na(A2star)] <- 0
  L1 <- data$race.y
  L2 <- data$age.y
  L3 <- data$CD8PFS
  L3[is.na(L3)] <- 0
  L4 <- data$CD4PFS
  L4[is.na(L4)] <- 0
  sw <- data$sw2
  delta2 <- function(beta1, beta6, beta7) {
    A2star + sigma_me2*(beta1 + beta6*L1 + beta7*L2)*Y
  }
  H <- function(x) {
    1 / (1 + exp(-x))
  }
  condexp <- function(beta0, beta1, beta2, beta3, beta4, beta5, beta6,
                      beta7) {
    H(beta0 + beta2*L1 + beta3*L2 + beta4*L3 + beta5*L4 +
        (beta1 + beta6*L1 + beta7*L2)*
        delta2(beta1, beta6, beta7) -
        ((beta1 + beta6*L1 + beta7*L2)^2 * sigma_me2) / 2)
  }
  function(theta) {
    c(sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8])),
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8]))*
        delta2(theta[2], theta[7], theta[8]),
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8]))*
        L1,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8]))*
        L2,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8]))*
        L3,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8]))*
        L4,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8]))*
        delta2(theta[2], theta[7], theta[8])*L1,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8]))*
        delta2(theta[2], theta[7], theta[8])*L2,
      theta[9] -
        exp(theta[1] + theta[2]*val + theta[3]*L1 + theta[4]*L2 +
          theta[5]*L3 + theta[6]*L4 + theta[7]*val*L1 + theta[8]*val*L2) /
        (1 + exp(theta[1] + theta[2]*val + theta[3]*L1 + theta[4]*L2 +
            theta[5]*L3 + theta[6]*L4 + theta[7]*val*L1 + theta[8]*val*L2))
    )
  }
}

# Calculate DR estimator for 4 ME and a grid of values
dr_ests1 <- array(NA, dim = c(4, 31))
dr_ests2 <- array(NA, dim = c(4, 23))
dr_se1 <- array(NA, dim = c(4, 31))
dr_se2 <- array(NA, dim = c(4, 23))

for (me in 0:3) {

  sigma_me1 <- 0.1*me
  sigma_me2 <- 0.1*me

  for (i in seq(0, 3, 0.1)) {

    guess1 <- 0.1
    if (i > 0) {
      guess1 <- coef(results_dr1)[9]
    }

    results_dr1 <-
      m_estimate(estFUN = eefun_dr1, data = analysisdat,
                 outer_args = list(i), compute_roots = TRUE,
                 root_control =
                   setup_root_control(start = c(coef(glmmod1), guess1)))
    dr_ests1[(me + 1), i*10 + 1] <- coef(results_dr1)[9]
    dr_se1[(me + 1), i*10 + 1] <- sqrt(vcov(results_dr1)[9, 9])

  }

  for (i in seq(0, 11, 0.5)) {

    guess2 <- 0.6
    if (i > 0) {
      guess2 <- coef(results_dr2)[9]
    }

    results_dr2 <-
      m_estimate(estFUN = eefun_dr2, data = analysisdat,
                 outer_args = list(i), compute_roots = TRUE,
                 root_control =
                   setup_root_control(start = c(coef(glmmod2), guess2)))
    dr_ests2[(me + 1), i*2 + 1] <- coef(results_dr2)[9]
    dr_se2[(me + 1), i*2 + 1] <- sqrt(vcov(results_dr2)[9, 9])

    if (is.na(coef(results_dr2)[9])) {
      guess2 <- min(guess2 + rnorm(1, 0, 0.1), 0.02)
      results_dr2 <-
        m_estimate(estFUN = eefun_dr2, data = analysisdat,
                   outer_args = list(i), compute_roots = TRUE,
                   root_control =
                     setup_root_control(start = c(coef(glmmod2), guess2)))
      dr_ests2[(me + 1), (i-4)*5 + 1] <- coef(results_dr2)[9]
      dr_se2[(me + 1), (i-4)*5 + 1] <- sqrt(vcov(results_dr2)[9, 9])
    }

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
            expression(paste("Low ME: ", sigma^2, "=0.1")),
            expression(paste("Moderate ME: ", sigma^2, "=0.2")),
            expression(paste("High ME: ", sigma^2, "=0.3")))

melist <- c("No ME: sigma^2 = 0", "Low ME: sigma^2 = 0.1",
            "Moderate ME: sigma^2 = 0.2", "High ME: sigma^2 = 0.3")

latdat <- data.frame(vals = c(rep(seq(0, 3, 0.1), 4),
                              rep(seq(0, 11, 0.5), 4)),
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
                     ME = c(rep((0:3)/10, each = 31),
                            rep((0:3)/10, each = 23)),
                     Exposure = c(rep("ADCP", 124), rep("RII", 92)))

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
             labeller = label_bquote(sigma^2 == .(ME))) +
             #labeller = labeller(ME = me.labs, Exposure = exp.labs)) +
  geom_ribbon(aes(ymin = Risk_low, ymax = Risk_upp), alpha = 0.3) +
  xlab("Exposure values") + ylab("HIV risk at study end") +
  ylim(c(-0.07, 0.93)) +
  theme_bw()


# Ignore below

# Estimating equations for both markers together at fixed ME
sigma_me1 <- 0.1
sigma_me2 <- 0.1

# New starting values
glmmod3 <-
  glm(HIVwk28preunbl.y ~ marker1 + race.y + age.y + CD8PFS + CD4PFS +
        marker1*race.y + marker1*age.y + marker2 + marker2*race.y +
        marker2*age.y,
      data = analysisdat, family = "binomial", weights = sw1*w2)

coefstart <- c(coef(glmmod3)[1:6], coef(glmmod3)[8:9], coef(glmmod3)[7],
               coef(glmmod3)[10:11])

# Estimator
eefun_dr3 <- function(data, val1, val2) {
  Y <- data$HIVwk28preunbl.y
  A1star <- data$marker1
  A1star[is.na(A1star)] <- 0
  A2star <- data$marker2
  A2star[is.na(A2star)] <- 0
  L1 <- data$race.y
  L2 <- data$age.y
  L3 <- data$CD8PFS
  L3[is.na(L3)] <- 0
  L4 <- data$CD4PFS
  L4[is.na(L4)] <- 0
  sw <- data$sw1*data$w2
  delta1 <- function(beta1, beta6, beta7) {
    A1star + sigma_me1*(beta1 + beta6*L1 + beta7*L2)*Y
  }
  delta2 <- function(beta8, beta9, beta10) {
    A2star + sigma_me2*(beta8 + beta9*L1 + beta10*L2)*Y
  }
  H <- function(x) {
    1 / (1 + exp(-x))
  }
  condexp <- function(beta0, beta1, beta2, beta3, beta4, beta5, beta6,
                      beta7, beta8, beta9, beta10) {
    H(beta0 + beta2*L1 + beta3*L2 + beta4*L3 + beta5*L4 +
        (beta1 + beta6*L1 + beta7*L2)*delta1(beta1, beta6, beta7) +
        (beta8 + beta9*L1 + beta10*L2)*delta2(beta8, beta9, beta10) -
        ((beta1 + beta6*L1 + beta7*L2)^2 * sigma_me1 +
         (beta8 + beta9*L1 + beta10*L2)^2 * sigma_me2) / 2)
  }
  function(theta) {
    c(sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9], theta[10],
                      theta[11])),
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9], theta[10],
                      theta[11]))*
        delta1(theta[2], theta[7], theta[8]),
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
        delta1(theta[2], theta[7], theta[8])*L1,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9], theta[10],
                      theta[11]))*
        delta1(theta[2], theta[7], theta[8])*L2,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9], theta[10],
                      theta[11]))*
        delta2(theta[9], theta[10], theta[11]),
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9], theta[10],
                      theta[11]))*
        delta2(theta[9], theta[10], theta[11])*L1,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9], theta[10],
                      theta[11]))*
        delta2(theta[9], theta[10], theta[11])*L2,
      theta[12] -
        exp(theta[1] + theta[2]*val1 + theta[3]*L1 + theta[4]*L2 +
            theta[5]*L3 + theta[6]*L4 + theta[7]*val1*L1 +
            theta[8]*val1*L2 + theta[9]*val2 +
            theta[10]*val2*L1 + theta[11]*val2*L2) /
        (1 + exp(theta[1] + theta[2]*val1 + theta[3]*L1 + theta[4]*L2 +
                 theta[5]*L3 + theta[6]*L4 + theta[7]*val1*L1 +
                 theta[8]*val1*L2 + theta[9]*val2 +
                 theta[10]*val2*L1 + theta[11]*val2*L2))
    )
  }
}

eefun_dr3 <- function(data, val1, val2) {
  Y <- data$HIVwk28preunbl.y
  A1star <- data$marker1
  A1star[is.na(A1star)] <- 0
  A2star <- data$marker2
  A2star[is.na(A2star)] <- 0
  L1 <- data$race.y
  L2 <- data$age.y
  L3 <- data$CD8PFS
  L3[is.na(L3)] <- 0
  L4 <- data$CD4PFS
  L4[is.na(L4)] <- 0
  sw <- data$sw1*data$w2
  delta1 <- function(beta1, beta6, beta7) {
    A1star + sigma_me1*(beta1 + beta6*L1 + beta7*L2)*Y
  }
  delta2 <- function(beta8) {
    A2star + sigma_me2*(beta8)*Y
  }
  H <- function(x) {
    1 / (1 + exp(-x))
  }
  condexp <- function(beta0, beta1, beta2, beta3, beta4, beta5, beta6,
                      beta7, beta8) {
    H(beta0 + beta2*L1 + beta3*L2 + beta4*L3 + beta5*L4 +
        (beta1 + beta6*L1 + beta7*L2)*delta1(beta1, beta6, beta7) +
        (beta8)*delta2(beta8) -
        ((beta1 + beta6*L1 + beta7*L2)^2 * sigma_me1 +
           (beta8)^2 * sigma_me2) / 2)
  }
  function(theta) {
    c(sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9])),
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9]))*
        delta1(theta[2], theta[7], theta[8]),
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9]))*
        L1,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9]))*
        L2,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9]))*
        L3,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9]))*
        L4,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9]))*
        delta1(theta[2], theta[7], theta[8])*L1,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9]))*
        delta1(theta[2], theta[7], theta[8])*L2,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9]))*
        delta2(theta[9]),
      theta[10] -
        exp(theta[1] + theta[2]*val1 + theta[3]*L1 + theta[4]*L2 +
              theta[5]*L3 + theta[6]*L4 + theta[7]*val1*L1 +
              theta[8]*val1*L2 + theta[9]*val2) /
        (1 + exp(theta[1] + theta[2]*val1 + theta[3]*L1 + theta[4]*L2 +
                   theta[5]*L3 + theta[6]*L4 + theta[7]*val1*L1 +
                   theta[8]*val1*L2 + theta[9]*val2))
    )
  }
}

# Get grid of values for mild ME to make 3D plot
dr_ests3 <- array(NA, dim = c(16, 21))
dr_se3 <- array(NA, dim = c(16, 21))

for (i in seq(0, 3, 0.2)) {

  for (j in seq(7, 11, 0.2)) {

    guess3 <- 0.4
    if (i > 0) {
      guess3 <- coef(results_dr3)[12]
    }

    results_dr3 <-
      m_estimate(estFUN = eefun_dr3, data = analysisdat,
                 outer_args = list(i, j), compute_roots = TRUE,
                 root_control =
                   setup_root_control(start = c(coefstart[1:9], guess3)))

    if (is.na(coef(results_dr3)[12])) {
      guess3 <- min(guess3 + rnorm(1, 0, 0.1), 0.02)
      results_dr3 <-
        m_estimate(estFUN = eefun_dr3, data = analysisdat,
                   outer_args = list(i, j), compute_roots = TRUE,
                   root_control =
                     setup_root_control(start = c(coefstart, guess3)))
    }

    dr_ests3[i*5 + 1, (i-7)*5 + 1] <- coef(results_dr3)[12]
    dr_se3[i*5 + 1, (i-7)*5 + 1] <- sqrt(vcov(results_dr3)[12, 12])

  }

}














gform_ests <- rep(NA, length(seq(0, 3, 0.1)))
cgf <- coef(results_gform1)

for (i in seq(0, 3, 0.1)) {

  gform_ests[i*10 + 1] <-
    exp(t(coef(results_gform1)[1:8]) %*% c(1, i, 1, 1, i, i)) /
    (1 + exp(t(coef(results_gform1)[1:8]) %*% c(1, i, 1, 1, i, i)))*
    mean(analysisdat$race.y)*mean(analysisdat$age.y) +
    exp(t(coef(results_gform1)[1:8]) %*% c(1, i, 1, 0, i, 0)) /
    (1 + exp(t(coef(results_gform1)[1:8]) %*% c(1, i, 1, 0, i, 0)))*
    mean(analysisdat$race.y)*(1 - mean(analysisdat$age.y)) +
    exp(t(coef(results_gform1)[1:8]) %*% c(1, i, 0, 1, 0, i)) /
    (1 + exp(t(coef(results_gform1)[1:8]) %*% c(1, i, 0, 1, 0, i)))*
    (1 - mean(analysisdat$race.y))*mean(analysisdat$age.y) +
    exp(t(coef(results_gform1)[1:8]) %*% c(1, i, 0, 0, 0, 0)) /
    (1 + exp(t(coef(results_gform1)[1:8]) %*% c(1, i, 0, 0, 0, 0)))*
    (1 - mean(analysisdat$race.y))*(1 - mean(analysisdat$age.y))

  t(coef(results_gform1)[1:8]) %*% c(1, i, 1, 1, i, i)

  cgf[1] + cgf[2]*i + cgf[3]*analysisdat

}

# IPW estimator



# Starting values
wmod1 <- glm(HIVwk28preunbl.y ~ marker1 + CD4PFS + CD8PFS, weights = sw1,
             data=analysisdat, family = "binomial")
wmod2 <- glm(HIVwk28preunbl.y ~ marker2 + CD4PFS + CD8PFS, weights = sw2,
             data=analysisdat, family = "binomial")

eefun_ipw1 <- function(data) {
  Y <- data$HIVwk28preunbl.y
  A1star <- data$marker1
  A1star[is.na(A1star)] <- 0
  L3 <- data$CD8PFS
  L3[is.na(L3)] <- 0
  L4 <- data$CD4PFS
  L4[is.na(L4)] <- 0
  sw <- data$sw1
  delta1 <- function(beta1) {
    A1star + beta1*sigma_me1*Y
  }
  H <- function(x) {
    1 / (1 + exp(-x))
  }
  condexp <- function(beta0, beta1, beta2, beta3) {
    H(beta0 + beta1*delta1(beta1) + beta2*L3 + beta3*L4 -
        (beta1^2 * sigma_me1) / 2)
  }
  function(theta) {
    c(sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4])),
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4]))*
        delta1(theta[2]),
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4]))*
        L3,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4]))*
        L4
    )
  }
}

results_ipw_csme1 <- m_estimate(estFUN = eefun_ipw1, data = analysisdat,
                                root_control =
                                  setup_root_control(start = coef(wmod)))

# DR estimator
eefun_dr1 <- function(data) {
  Y <- data$HIVwk28preunbl.y
  A1star <- data$marker1
  A1star[is.na(A1star)] <- 0
  L1 <- data$race.y
  L2 <- data$age.y
  L3 <- data$CD8PFS
  L3[is.na(L3)] <- 0
  L4 <- data$CD4PFS
  L4[is.na(L4)] <- 0
  sw <- data$sw1
  delta1 <- function(beta1, beta6, beta7) {
    A1star + sigma_me1*(beta1 + beta6*L1 + beta7*L2)*Y
  }
  H <- function(x) {
    1 / (1 + exp(-x))
  }
  condexp <- function(beta0, beta1, beta2, beta3, beta4, beta5, beta6,
                      beta7) {
    H(beta0 + beta2*L1 + beta3*L2 + beta4*L3 + beta5*L4 +
        (beta1 + beta6*L1 + beta7*L2)*
        delta1(beta1, beta6, beta7) -
        ((beta1 + beta6*L1 + beta7*L2)^2 * sigma_me1) / 2)
  }
  function(theta) {
    c(sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8])),
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8]))*
        delta1(theta[2], theta[7], theta[8]),
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8]))*
        L1,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8]))*
        L2,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8]))*
        L3,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8]))*
        L4,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8]))*
        delta1(theta[2], theta[7], theta[8])*L1,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8]))*
        delta1(theta[2], theta[7], theta[8])*L2
    )
  }
}

results_dr1 <-
  m_estimate(estFUN = eefun_dr1, data = analysisdat,
             compute_roots = TRUE,
             root_control =
               setup_root_control(start = c(coef(glmmod1))))



























# IPW LR no ME for original model
testmodipw <- glm(HIVwk28preunbl.y ~ marker1,
               data = analysisdat, weights = analysisdat$sw1,
               family = binomial)

# DR no ME for original model
testmoddr <- glm(HIVwk28preunbl.y ~ marker1*age.y + marker1*race.y,
                 data = analysisdat, weights = analysisdat$sw1,
                 family = binomial)

testmoddr <- glm(HIVwk28preunbl.y ~ marker2*age.y + marker2*race.y,
                 data = analysisdat, weights = analysisdat$sw2,
                 family = binomial)

# DR no ME for cell-adjusted model
testmoddr <- glm(HIVwk28preunbl.y ~ marker1*age.y + marker1*race.y +
                   marker1*CD8PFS + CD4PFS,
                 data = analysisdat, weights = analysisdat$sw1,
                 family = binomial)

testmoddr <- glm(HIVwk28preunbl.y ~ marker2*age.y + marker2*race.y +
                   marker2*CD8PFS + CD4PFS,
                 data = analysisdat, weights = analysisdat$sw2,
                 family = binomial)

# Application for marker 1 and 2 separately
# Need to adjust for var est in EE stack
# Four EE stacks coeesponding to four analyses


# First EE stack
# ADCP exposure not adjusting for CD4 and CD8
eefun_csme_aipw1_nocell <- function(data) {
  Y <- data$HIVwk28preunbl.y
  X1star <- data$marker1
  X1star[is.na(X1star)] <- 0
  L1 <- data$race.y
  L2 <- data$age.y
  sw <- data$sw1
  delta1 <- function(beta1, beta4, beta5) {
    X1star + sigma_me1*(beta1 + beta4*L1 + beta5*L2)*Y
  }
  H <- function(x) {
    1 / (1 + exp(-x))
  }
  condexp <- function(beta0, beta1, beta2, beta3, beta4, beta5) {
    H(beta0 + beta2*L1 + beta3*L2 +
        (beta1 + beta4*L1 + beta5*L2)*
        delta1(beta1, beta4, beta5) -
        ((beta1 + beta4*L1 + beta5*L2)^2 * sigma_me1) / 2)
  }
  function(theta) {
    c(sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6])),
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6]))*
        delta1(theta[2], theta[5], theta[6]),
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6]))*
        L1,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6]))*
        L2,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6]))*
        delta1(theta[2], theta[5], theta[6])*L1,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6]))*
        delta1(theta[2], theta[5], theta[6])*L2,
      L1 - theta[7],
      L2 - theta[8],
      theta[2] + theta[5]*theta[7] + theta[6]*theta[8] - theta[9]
    )
  }
}


# Second EE stack
# ADCP exposure adjusting for CD4 and CD8
eefun_csme_aipw1_cell <- function(data) {
  Y <- data$HIVwk28preunbl.y
  X1star <- data$marker1
  X1star[is.na(X1star)] <- 0
  L1 <- data$race.y
  L2 <- data$age.y
  L3 <- data$CD8PFS
  L3[is.na(L3)] <- 0
  L4 <- data$CD4PFS
  L4[is.na(L4)] <- 0
  sw <- data$sw1
  delta1 <- function(beta1, beta6, beta7, beta8) {
    X1star + sigma_me1*(beta1 + beta6*L1 + beta7*L2 + beta8*L3)*Y
  }
  H <- function(x) {
    1 / (1 + exp(-x))
  }
  condexp <- function(beta0, beta1, beta2, beta3, beta4, beta5, beta6,
                      beta7, beta8) {
    H(beta0 + beta2*L1 + beta3*L2 + beta4*L3 + beta5*L4 +
        (beta1 + beta6*L1 + beta7*L2 + beta8*L3)*
        delta1(beta1, beta6, beta7, beta8) -
        ((beta1 + beta6*L1 + beta7*L2 + beta8*L3)^2 * sigma_me1) / 2)
  }
  function(theta) {
    c(sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9])),
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9]))*
        delta1(theta[2], theta[7], theta[8], theta[9]),
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9]))*
        L1,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9]))*
        L2,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9]))*
        L3,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9]))*
        L4,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9]))*
        delta1(theta[2], theta[7], theta[8], theta[9])*L1,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9]))*
        delta1(theta[2], theta[7], theta[8], theta[9])*L2,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9]))*
        delta1(theta[2], theta[7], theta[8], theta[9])*L3,
      L1 - theta[10],
      L2 - theta[11],
      sw*(L3 - theta[12]),
      theta[2] + theta[7]*theta[10] + theta[8]*theta[11] + theta[9]*theta[12] - theta[13]
    )
  }
}

# Alternate second stack
eefun_csme_aipw1_cell <- function(data) {
  Y <- data$HIVwk28preunbl.y
  X1star <- data$marker1
  X1star[is.na(X1star)] <- 0
  L1 <- data$race.y
  L2 <- data$age.y
  L3 <- data$CD8PFS
  L3[is.na(L3)] <- 0
  L4 <- data$CD4PFS
  L4[is.na(L4)] <- 0
  sw <- data$sw1
  delta1 <- function(beta1, beta6, beta7) {
    X1star + sigma_me1*(beta1 + beta6*L1 + beta7*L2)*Y
  }
  H <- function(x) {
    1 / (1 + exp(-x))
  }
  condexp <- function(beta0, beta1, beta2, beta3, beta4, beta5, beta6,
                      beta7) {
    H(beta0 + beta2*L1 + beta3*L2 + beta4*L3 + beta5*L4 +
        (beta1 + beta6*L1 + beta7*L2)*
        delta1(beta1, beta6, beta7) -
        ((beta1 + beta6*L1 + beta7*L2)^2 * sigma_me1) / 2)
  }
  function(theta) {
    c(sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8])),
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8]))*
        delta1(theta[2], theta[7], theta[8]),
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8]))*
        L1,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8]))*
        L2,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8]))*
        L3,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8]))*
        L4,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8]))*
        delta1(theta[2], theta[7], theta[8])*L1,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8]))*
        delta1(theta[2], theta[7], theta[8])*L2,
      L1 - theta[9],
      L2 - theta[10],
      theta[2] + theta[7]*theta[9] + theta[8]*theta[10] - theta[11]
    )
  }
}


# Third EE stack
# R2 exposure not adjusting for CD4 and CD8
eefun_csme_aipw2_nocell <- function(data) {
  Y <- data$HIVwk28preunbl.y
  X2star <- data$marker2
  #X2star[X2star < 5 & !is.na(X2star)] <- median(X2star[X2star > 0])
  X2star[is.na(X2star)] <- 0
  L1 <- data$race.y
  L2 <- data$age.y
  #L3 <- data$BMI.y
  sw <- data$sw2
  delta2 <- function(beta1, beta4, beta5) {
    X2star + sigma_me2*(beta1 + beta4*L1 + beta5*L2)*Y
  }
  H <- function(x) {
    1 / (1 + exp(-x))
  }
  condexp <- function(beta0, beta1, beta2, beta3, beta4, beta5) {
    H(beta0 + beta2*L1 + beta3*L2 +
        (beta1 + beta4*L1 + beta5*L2)*
        delta2(beta1, beta4, beta5) -
        ((beta1 + beta4*L1 + beta5*L2)^2 * sigma_me2) / 2)
  }
  function(theta) {
    c(sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6])),
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6]))*
        delta2(theta[2], theta[5], theta[6]),
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6]))*
        L1,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6]))*
        L2,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6]))*
        delta2(theta[2], theta[5], theta[6])*L1,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6]))*
        delta2(theta[2], theta[5], theta[6])*L2,
      L1 - theta[7],
      L2 - theta[8],
      theta[2] + theta[5]*theta[7] + theta[6]*theta[8] - theta[9]
    )
  }
}


# Fourth EE stack
# R2 exposure adjusting for CD4 and CD8
eefun_csme_aipw2_cell <- function(data) {
  Y <- data$HIVwk28preunbl.y
  X2star <- data$marker2
  X2star[is.na(X2star)] <- 0
  L1 <- data$race.y
  L2 <- data$age.y
  L3 <- data$CD8PFS
  L3[is.na(L3)] <- 0
  L4 <- data$CD4PFS
  L4[is.na(L4)] <- 0
  sw <- data$sw2
  delta2 <- function(beta1, beta6, beta7, beta8) {
    X2star + sigma_me2*(beta1 + beta6*L1 + beta7*L2 + beta8*L3)*Y
  }
  H <- function(x) {
    1 / (1 + exp(-x))
  }
  condexp <- function(beta0, beta1, beta2, beta3, beta4, beta5, beta6,
                      beta7, beta8) {
    H(beta0 + beta2*L1 + beta3*L2 + beta4*L3 + beta5*L4 +
        (beta1 + beta6*L1 + beta7*L2 + beta8*L3)*
        delta2(beta1, beta6, beta7, beta8) -
        ((beta1 + beta6*L1 + beta7*L2 + beta8*L3)^2 * sigma_me2) / 2)
  }
  function(theta) {
    c(sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9])),
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9]))*
        delta2(theta[2], theta[7], theta[8], theta[9]),
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9]))*
        L1,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9]))*
        L2,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9]))*
        L3,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9]))*
        L4,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9]))*
        delta2(theta[2], theta[7], theta[8], theta[9])*L1,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9]))*
        delta2(theta[2], theta[7], theta[8], theta[9])*L2,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8], theta[9]))*
        delta2(theta[2], theta[7], theta[8], theta[9])*L3,
      L1 - theta[10],
      L2 - theta[11],
      sw*(L3 - theta[12]),
      theta[2] + theta[7]*theta[10] + theta[8]*theta[11] + theta[9]*theta[12] - theta[13]
    )
  }
}


# Alternative fourth stack
eefun_csme_aipw2_cell <- function(data) {
  Y <- data$HIVwk28preunbl.y
  X2star <- data$marker2
  X2star[is.na(X2star)] <- 0
  L1 <- data$race.y
  L2 <- data$age.y
  L3 <- data$CD8PFS
  L3[is.na(L3)] <- 0
  L4 <- data$CD4PFS
  L4[is.na(L4)] <- 0
  sw <- data$sw2
  delta2 <- function(beta1, beta6, beta7) {
    X2star + sigma_me2*(beta1 + beta6*L1 + beta7*L2)*Y
  }
  H <- function(x) {
    1 / (1 + exp(-x))
  }
  condexp <- function(beta0, beta1, beta2, beta3, beta4, beta5, beta6,
                      beta7) {
    H(beta0 + beta2*L1 + beta3*L2 + beta4*L3 + beta5*L4 +
        (beta1 + beta6*L1 + beta7*L2)*
        delta2(beta1, beta6, beta7) -
        ((beta1 + beta6*L1 + beta7*L2)^2 * sigma_me2) / 2)
  }
  function(theta) {
    c(sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8])),
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8]))*
        delta2(theta[2], theta[7], theta[8]),
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8]))*
        L1,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8]))*
        L2,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8]))*
        L3,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8]))*
        L4,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8]))*
        delta2(theta[2], theta[7], theta[8])*L1,
      sw*(Y - condexp(theta[1], theta[2], theta[3], theta[4], theta[5],
                      theta[6], theta[7], theta[8]))*
        delta2(theta[2], theta[7], theta[8])*L2,
      L1 - theta[9],
      L2 - theta[10],
      theta[2] + theta[7]*theta[9] + theta[8]*theta[10] - theta[11]
    )
  }
}



# Make analysis plot
library(plotrix)
sigma1list <- seq(0, 0.3, 0.05)
sigma2list <- seq(0, 0.3, 0.05)
#sigma1list <- seq(0, 0.3, 0.3)
#sigma2list <- seq(0, 0.3, 0.3)

marker1ests_nocell <- rep(NA, length(sigma1list))
marker1lower_nocell <- rep(NA, length(sigma1list))
marker1upper_nocell <- rep(NA, length(sigma1list))
marker1ests_cell <- rep(NA, length(sigma1list))
marker1lower_cell <- rep(NA, length(sigma1list))
marker1upper_cell <- rep(NA, length(sigma1list))
marker1ests_cell_int <- rep(NA, length(sigma1list))
marker1lower_cell_int <- rep(NA, length(sigma1list))
marker1upper_cell_int <- rep(NA, length(sigma1list))

marker2ests_nocell <- rep(NA, length(sigma2list))
marker2lower_nocell <- rep(NA, length(sigma2list))
marker2upper_nocell <- rep(NA, length(sigma2list))
marker2ests_cell <- rep(NA, length(sigma2list))
marker2lower_cell <- rep(NA, length(sigma2list))
marker2upper_cell <- rep(NA, length(sigma2list))
marker2ests_cell_int <- rep(NA, length(sigma2list))
marker2lower_cell_int <- rep(NA, length(sigma2list))
marker2upper_cell_int <- rep(NA, length(sigma2list))

for (i in 1:7) {

  sigma_me1 <- sigma1list[i]
  sigma_me2 <- sigma2list[i]

  # Use Neidich results for starting value, then switch to previous ME mods
  startval_1_nocell <- log(0.47)
  startval_1_cell <- log(0.59)
  startval_2_nocell <- log(0.48)
  startval_2_cell <- log(0.63)

  if (i > 1) {
    startval_1_nocell <- coef(results_csme_aipw1_nocell)[9]
    startval_1_cell <- coef(results_csme_aipw1_cell)[11]
    startval_2_nocell <- coef(results_csme_aipw2_nocell)[9]
    startval_2_cell <- coef(results_csme_aipw2_cell)[11]
  }

  # ADCP, no cell adjustment
  results_csme_aipw1_nocell <-
    m_estimate(estFUN = eefun_csme_aipw1_nocell, data = analysisdat,
               compute_roots = TRUE,
               root_control =
                  setup_root_control(start = c(coef(glmmod1_nocell),
                                               mean(analysisdat$race.y),
                                               mean(analysisdat$age.y),
                                               startval_1_nocell)))

  marker1ests_nocell[i] <- exp(coef(results_csme_aipw1_nocell)[9])
  se1 <- sqrt(vcov(results_csme_aipw1_nocell)[9, 9])
  marker1lower_nocell[i] <- exp(coef(results_csme_aipw1_nocell)[9] - 1.96*se1)
  marker1upper_nocell[i] <- exp(coef(results_csme_aipw1_nocell)[9] + 1.96*se1)

  print("1")
  # ADCP, CD4 and CD8 adjustment
  results_csme_aipw1_cell <-
    m_estimate(estFUN = eefun_csme_aipw1_cell, data = analysisdat,
               compute_roots = TRUE,
               root_control =
                 setup_root_control(start = c(coef(glmmod1_cell),
                                              mean(analysisdat$race.y),
                                              mean(analysisdat$age.y),
                                              startval_1_cell)))

  marker1ests_cell[i] <- exp(coef(results_csme_aipw1_cell)[11])
  se1 <- sqrt(vcov(results_csme_aipw1_cell)[11, 11])
  marker1lower_cell[i] <- exp(coef(results_csme_aipw1_cell)[11] - 1.96*se1)
  marker1upper_cell[i] <- exp(coef(results_csme_aipw1_cell)[11] + 1.96*se1)
  print("2")
  # R2, no cell adjustment
  results_csme_aipw2_nocell <-
    m_estimate(estFUN = eefun_csme_aipw2_nocell, data = analysisdat,
               compute_roots = TRUE,
               root_control =
                 setup_root_control(start = c(coef(glmmod2_nocell),
                                              mean(analysisdat$race.y),
                                              mean(analysisdat$age.y),
                                              startval_2_nocell)))

  marker2ests_nocell[i] <- exp(coef(results_csme_aipw2_nocell)[9])
  se2 <- sqrt(vcov(results_csme_aipw2_nocell)[9, 9])
  marker2lower_nocell[i] <- exp(coef(results_csme_aipw2_nocell)[9] - 1.96*se2)
  marker2upper_nocell[i] <- exp(coef(results_csme_aipw2_nocell)[9] + 1.96*se2)
  print("3")
  # R2, CD4 and CD8 adjustment
  results_csme_aipw2_cell <-
    m_estimate(estFUN = eefun_csme_aipw2_cell, data = analysisdat,
               compute_roots = TRUE,
               root_control =
                 setup_root_control(start = c(coef(glmmod2_cell),
                                              mean(analysisdat$race.y),
                                              mean(analysisdat$age.y),
                                              startval_2_cell)))

  marker2ests_cell[i] <- exp(coef(results_csme_aipw2_cell)[11])
  se2 <- sqrt(vcov(results_csme_aipw2_cell)[11, 11])
  marker2lower_cell[i] <- exp(coef(results_csme_aipw2_cell)[11] - 1.96*se2)
  marker2upper_cell[i] <- exp(coef(results_csme_aipw2_cell)[11] + 1.96*se2)

  print("4")

}

par(mfrow = c(2,2))
plotCI(seq(0, 0.5, 0.5/6), marker1ests_nocell, li = marker1lower_nocell,
       ui = marker1upper_nocell, xlab = "ME as proportion of intra-vaccinee variance",
       ylab = "Odds ratio", main = "ADCP without cellular adjustment",
       xaxt = "n")
axis(1, at = c(0, 0.25, 0.5))
plotCI(seq(0, 0.5, 0.5/6), marker1ests_cell, li = marker1lower_cell,
       ui = marker1upper_cell, xlab = "ME as proportion of intra-vaccinee variance",
       ylab = "Odds ratio", main = "ADCP with cellular adjustment",
       xaxt = "n")
axis(1, at = c(0, 0.25, 0.5))
plotCI(seq(0, 0.5, 0.5/6), marker2ests_nocell, li = marker2lower_nocell,
       ui = marker2upper_nocell, xlab = "ME as proportion of intra-vaccinee variance",
       ylab = "Odds ratio", main = "RII without cellular adjustment",
       xaxt = "n")
axis(1, at = c(0, 0.25, 0.5))
plotCI(seq(0, 0.5, 0.5/6), marker2ests_cell, li = marker2lower_cell,
       ui = marker2upper_cell, xlab = "ME as proportion of intra-vaccinee variance",
       ylab = "Odds ratio", main = "RII with cellular adjustment",
       xaxt = "n")
axis(1, at = c(0, 0.25, 0.5))
