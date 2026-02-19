################
# NCC Cov adj
# ANALYSIS
# January 2026
################

# rm(list=ls())
# setwd("C:/Users/marta/Dropbox/C5/RESEARCH/RCode/2025-12-NCC-covadj")
# source("simdata_uneqblock.R")

################
# lm
lmmodel <- function(data, trt = 1, dataset="ACA") {
  
  cc_periods <- which(as.vector(table(data$t, data$p)[trt+1,])>0)
  
  # analysis data 
  if(dataset=="ACA"){
    d <- subset(data, t %in% c(0,trt) & p %in% cc_periods)  # ACA data
  }
  if(dataset=="LACA"){
    d <- subset(data, t %in% c(0,trt) & p %in% max(cc_periods)) # LACA data
  }
  if(dataset=="ECE"){
    d <- subset(data, p %in% cc_periods) # ECE data
  }
  if(dataset=="NCC"){
    d <- data # NCC data
    if(trt==1){
      d <- subset(data, p %in% cc_periods) # ECE data
    }
  } 
  
  # linear model
  fit <- lm(y ~ as.factor(t) + as.factor(p) + x, data = d)
  coef_name <- paste0("as.factor(t)", trt)
  lm_coef <- as.numeric(coef(fit)[coef_name])  # ctrl = 0 baseline
  
  # Variance–covariance matrix of coefficients
  vc <- vcov(fit)
  coef_var <- vc[coef_name, coef_name] # Variance of the coefficient
  coef_sd <- sqrt(coef_var) # Standard deviation (standard error)
  
  # estimand
  esd_m <- as.numeric(colMeans(d[2+as.numeric(trt)]))-mean(d$y0)
  
  return(c(est=lm_coef, var=coef_var, se=coef_sd, esd_m=esd_m))
}

################
# g-estimation
g_estimate <- function(data, trt = 1, dataset="ACA") {
  
  cc_periods <- which(as.vector(table(data$t, data$p)[trt+1,])>0)
  
  # data to predict marginal means
  if(dataset=="ACA"){
    d <- subset(data, t %in% c(0,trt) & p %in% cc_periods)  # ACA data
  }
  if(dataset=="LACA"){
    d <- subset(data, t %in% c(0,trt) & p %in% max(cc_periods)) # LACA data
  }
  if(dataset=="ECE"){
    d <- subset(data, p %in% cc_periods) # ECE data
  }
  if(dataset=="NCC"){
    d <- data # NCC data
    if(trt==1){
      d <- subset(data, p %in% cc_periods) # ECE data
    }
  }   
  
  # outcome model
  fit <- lm(y ~ as.factor(t) + as.factor(p) + x, data = d)
  
  # Model matrix for each treatment
  X1 <- model.matrix(fit, data = transform(d, t = factor(trt, levels(d$t))))
  X0 <- model.matrix(fit, data = transform(d, t = factor(0,  levels(d$t))))
  
  # Average contrast vector
  cvec <- colMeans(X1 - X0)
  
  # Point estimate
  beta_hat <- coef(fit)
  g_diffmeans <- sum(cvec * beta_hat)
  
  # Variance via delta method
  Vbeta <- vcov(fit)
  g_var <- as.numeric(t(cvec) %*% Vbeta %*% cvec)
  g_se <- sqrt(g_var)
  
  # # Create predicted outcomes for each treatment
  # d$mu_hat_0 <- predict(fit, newdata = transform(d, t = factor(0, levels = levels(d$t))))
  # d$mu_hat_1 <- predict(fit, newdata = transform(d, t = factor(trt, levels = levels(d$t))))
  # 
  # # g_diffmeans <- mean(d$mu_hat_1) - mean(d$mu_hat_0)
  # 
  # # Individual-level predicted effects
  # d$gdiff <- d$mu_hat_1 - d$mu_hat_0
  # 
  # # Mean 
  # g_diffmeans <- mean(d$gdiff)
  # 
  # # Variance and SE of the estimator
  # n <- nrow(d)
  # g_var <- var(d$gdiff)/n
  # g_se <- sqrt(g_var) 
  
  # estimand
  esd_m <- as.numeric(colMeans(d[2+as.numeric(trt)]))-mean(d$y0)
  
  return(c(est=g_diffmeans, var=g_var, se=g_se, esd_m=esd_m))
}

################
# aipw estimation
aipw_estimate <- function(data, trt = 1, period_of_interest = NULL) {
  
  # Subset to relevant treatments
  cc_periods <- which(as.vector(table(data$t, data$p)[trt+1,])>0)
  d <- subset(data, t %in% c(0,trt) & p %in% cc_periods) # by default ACA data
  
  # cc_periods <- which(as.vector(table(data$t, data$p)[trt+1,])>0)
  # d <- subset(data, p %in% cc_periods) # ECE data
  
  # Propensity scores: known from randomization
  # table(data$t, data$p)
  # pi_hat_p1 <- table(data$t, data$p)[,1]/sum(table(data$t, data$p)[,1])
  # pi_hat_p2 <- table(data$t, data$p)[,2]/sum(table(data$t, data$p)[,2])
  # pi_hat_p3 <- table(data$t, data$p)[,3]/sum(table(data$t, data$p)[,3])
  # pi_hat <- matrix(c(pi_hat_p1, pi_hat_p2, pi_hat_p3), nrow=3)
  
  pi_hat <- numeric(nrow(d))
  # create pi_hat per ind: loop over periods
  for (period_level in levels(d$p)) {
    idx <- which(d$p == period_level)
    t_counts <- table(d$t[idx])
    prop <- t_counts / sum(t_counts)  # pi(A | S=period)
    
    # assign propensity to each ind
    pi_hat[idx] <- prop[as.character(d$t[idx])]
  }
  
  d$pi_hat <- pi_hat
  
  # Outcome model 
  fit <- lm(y ~ as.factor(t) + as.factor(p) + x, data = d)
  
  # Predicted outcomes for treatment 0 and treatment of interest
  d$mu_hat_0 <- predict(fit, newdata = transform(d, t = factor(0, levels = levels(d$t))))
  d$mu_hat_a <- predict(fit, newdata = transform(d, t = factor(trt, levels = levels(d$t))))
  
  # Compute individual contributions
  psi_i <- ((d$t == trt) / d$pi_hat) * (d$y - d$mu_hat_a) -
    ((d$t == 0) / d$pi_hat) * (d$y - d$mu_hat_0) +
    (d$mu_hat_a - d$mu_hat_0)
  
  # If period_of_interest is specified, subset and stabilize
  if (!is.null(period_of_interest)) {
    I_P <- which(d$p %in% period_of_interest) 
    psi_hat <- mean(psi_i[I_P])
  } else {
    psi_hat <- mean(psi_i)
  } 
  
  # # Mean and SD
  # psi_hat <- mean(psi_i)
  # gsd <- sd(psi_i)
  # gse <- gsd / sqrt(length(psi_i))
  
  n <- length(psi_i)
  
  # Empirical variance of influence function
  var_psi <- var(psi_i)
  
  # Variance and SE of the estimator
  var_hat <- var_psi/n
  se_hat  <- sqrt(var_hat)
  
  # estimand
  esd_m <- as.numeric(colMeans(d[2+as.numeric(trt)]))-mean(d$y0)
  
  return(c(est=psi_hat, var=var_hat, se=se_hat, esd_m=esd_m))
}


# EXAMPLES
# dat <- simdata_blocked_unequal(
#   mu0 = 0,
#   mu1 = 1,
#   mu2 = 2,
#   beta=2,
#   N1 = 120,
#   N2 = 300,
#   N3 = 200,
#   alloc1 = c(2/3, 1/3),        # 2:1
#   alloc2 = c(3/6, 2/6, 1/6),   # 3:2:1
#   alloc3 = c(4/5, 1/5),
#   # alloc1 = rep(1,2)*1/2, # eq alloc
#   # alloc2 = rep(1,3)*1/3, # eq alloc
#   # alloc3 = rep(1,2)*1/2, # eq alloc
#   lambda = 0.5
# )
# 
# # linear model
# lmmodel(dat, trt=1)
# 
# # g-estimator
# gest <- g_estimate(dat, trt=1)
# gest
# 
# # aipw
# aipw_estimate(dat, trt = 1)
# aipw_estimate(dat, trt = 1, period_of_interest = c(1,2))
# aipw_estimate(dat, trt = 1, period_of_interest = c(2))
# aipw_estimate(dat, trt = 1, period_of_interest = c(1,2,3))
