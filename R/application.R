# Build application data set
rm(list = ls())
load("C:/Users/bblette1/Downloads/HVTN505_2019-08-08/HVTN505/data/dat.505.rda")

assays <- subset(var.505, assay %in% c("fcrR2a", "fcrR3a", "phago"))

primarydat <- read.csv("C:/Users/bblette1/Documents/primary505_for_sharing .csv")
primarydat$ptid <- primarydat$pub_id

fulldat <- merge(dat.505, primarydat, by = "ptid", all = T)
fulldat$marker1 <- log(fulldat$ADCP1)
fulldat$marker2 <- fulldat$R2aConSgp140CFI
fulldat$marker3 <- fulldat$R3agp41
fulldat$marker3 <- fulldat$R3aConSgp140CFI

set.seed(1234)

# Simple impute missing confounder, marker2 0 value
fulldat$BMI.y[is.na(fulldat$BMI.y)] <- median(fulldat$BMI.y, na.rm = T)
imputemax <- min(fulldat$marker2[fulldat$marker2 > 0], na.rm = T)
fulldat$marker2[fulldat$marker2 == 0 & !is.na(fulldat$marker2) &
                fulldat$trt.y == 1] <- runif(1, 0, imputemax)

fulldat$sampled <- !is.na(fulldat$marker1)
sampmod <- glm(sampled ~ HIVwk28preunbl.y, data = fulldat[fulldat$trt.y == 1, ])

# Create analysis dataset
# Restrict to treatment, use all three markers
# Where is behavior risk variable?
# Do we need to do anything with CD4 and CD8
# OR from 2019 paper vs DR estimator vs DR with ME
# Different specifications for OR and IPW, leave one confounder out (BMI)?
library(dplyr)
fulldat$CD4PFS <- fulldat$CD4_ANYVRCENV_PolyfunctionalityScore_score_bin
fulldat$CD8PFS <- fulldat$CD8_ANYVRCENV_PolyfunctionalityScore_score_bin
analysisdat <- fulldat %>%
  filter(trt.y == 1) %>%
  select(HIVall.y, HIVwk28preunbl.y, marker1, marker2, marker3,
         race.y, BMI.y, age.y, bhvrisk.y, ADCP1_bin, CD4PFS, CD8PFS)

# Probability selected by treatment
# Weights based on HIV at week 28, not end of study
analysisdat$estweights <-
  (!is.na(analysisdat$marker1))*(analysisdat$HIVwk28preunbl.y == 1) / (25 / sum(fulldat$HIVwk28preunbl.y[fulldat$trt.y == 1])) +
  (!is.na(analysisdat$marker1))*(analysisdat$HIVwk28preunbl.y == 0) / (125 / sum(1 - fulldat$HIVwk28preunbl.y[fulldat$trt.y == 1]))

# Different race yields same model?

testmod <- glm(HIVwk28preunbl.y ~ marker1 + age.y + BMI.y + race.y + bhvrisk.y,
               data = analysisdat, weights = analysisdat$estweights,
               family = binomial)

testmod <- glm(HIVwk28preunbl.y ~ marker1 + CD8PFS + CD4PFS +age.y + BMI.y + race.y + bhvrisk.y,
               data = analysisdat, weights = analysisdat$estweights,
               family = binomial)

testmod <- glm(HIVwk28preunbl.y ~ marker2 + age.y + BMI.y + race.y + bhvrisk.y,
               data = analysisdat, weights = analysisdat$estweights,
               family = binomial)
# Look at both weighted and unweighted

# In a pre-specified secondary analysis reported in 2019, an HIV-1 specific Env protein
# VRC B gp140 had the strongest Fcgamma2ra/3a association with HIV acquisition
# out of an array of measured proteins. These assays may be subject to ME
# nd an additive ME model seems plausible for each of the three biomarkers.
# No supp data, snes analysis etc.


# Outcome regressions for starting values
glmmod1_nocell <-
  glm(HIVwk28preunbl.y ~ marker1 + race.y + age.y +
                         marker1*race.y + marker1*age.y,
      data = analysisdat, family = "binomial", weights = estweights)

glmmod1_cell <-
  glm(HIVwk28preunbl.y ~ marker1 + race.y + age.y + CD8PFS + CD4PFS + 
                         marker1*race.y + marker1*age.y,
      data = analysisdat, family = "binomial", weights = estweights)

glmmod2_nocell <-
  glm(HIVwk28preunbl.y ~ marker2 + race.y + age.y +
        marker2*race.y + marker2*age.y,
      data = analysisdat, family = "binomial", weights = estweights)

glmmod2_cell <-
  glm(HIVwk28preunbl.y ~ marker2 + race.y + age.y + CD8PFS + CD4PFS +
        marker2*race.y + marker2*age.y,
      data = analysisdat, family = "binomial", weights = estweights)

# Estimate IPTW
denom_mod1 <- lm(marker1 ~ race.y + BMI.y + age.y + bhvrisk.y,
                 data = analysisdat)
p_denom1 <- predict(denom_mod1, type='response')
dens_denom1 <- dnorm(analysisdat$marker1, p_denom1, summary(denom_mod1)$sigma)
num_mod1 <- lm(marker1 ~ 1, data = analysisdat)
p_num1 <- predict(num_mod1, type='response')
dens_num1 <- dnorm(analysisdat$marker1, p_num1, summary(denom_mod1)$sigma)
analysisdat$w1 <- dens_num1 / dens_denom1
analysisdat$w1[is.na(analysisdat$w1)] <- 0

analysisdat$sw1 <- analysisdat$w1*analysisdat$estweights

denom_mod2 <- lm(marker2 ~ race.y + BMI.y + age.y + bhvrisk.y,
                 data = analysisdat)
p_denom2 <- predict(denom_mod2, type='response')
dens_denom2 <- dnorm(analysisdat$marker2, p_denom2, summary(denom_mod2)$sigma)
num_mod2 <- lm(marker2 ~ 1, data = analysisdat)
p_num2 <- predict(num_mod2, type='response')
dens_num2 <- dnorm(analysisdat$marker2, p_num2, summary(denom_mod2)$sigma)
analysisdat$w2 <- dens_num2 / dens_denom2
analysisdat$w2[is.na(analysisdat$w2)] <- 0

analysisdat$sw2 <- analysisdat$w2*analysisdat$estweights

# Which HIV variable to use as outcome?

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
library(geex)


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
plotCI(seq(0, 0.3, 0.05), marker1ests_nocell, li = marker1lower_nocell,
       ui = marker1upper_nocell, xlab = "ME variance", ylab = "Odds ratio",
       main = "ADCP without cellular adjustment")
plotCI(seq(0, 0.3, 0.05), marker1ests_cell, li = marker1lower_cell,
       ui = marker1upper_cell, xlab = "ME variance", ylab = "Odds ratio",
       main = "ADCP with cellular adjustment")
plotCI(seq(0, 0.3, 0.05), marker2ests_nocell, li = marker2lower_nocell,
       ui = marker2upper_nocell, xlab = "ME variance", ylab = "Odds ratio",
       main = "R2 without cellular adjustment")
plotCI(seq(0, 0.3, 0.05), marker2ests_cell, li = marker2lower_cell,
       ui = marker2upper_cell, xlab = "ME variance", ylab = "Odds ratio",
       main = "R2 with cellular adjustment")
