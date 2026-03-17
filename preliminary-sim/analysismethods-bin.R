################
# NCC Cov adj
# ANALYSIS
# March 2026
################

# rm(list=ls())
# setwd("C:/Users/marta/Dropbox/C5/RESEARCH/RCode/2025-12-NCC-covadj")
# source("simdata_uneqblock.R")

################
# glm
glmmodel <- function(data, trt = 1, dataset="ACA") {
  
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
  fit <- glm(y ~ as.factor(t) + as.factor(p) + x,
             data = d, family = binomial(link = "logit"))
  coef_name <- paste0("as.factor(t)", trt)
  glm_coef <- as.numeric(coef(fit)[coef_name])  # ctrl = 0 baseline
  
  # Variance–covariance matrix of coefficients
  vc <- vcov(fit)
  coef_var <- vc[coef_name, coef_name] # Variance of the coefficient
  coef_sd <- sqrt(coef_var) # Standard deviation (standard error)
  
  # difference in proportions
  p_hat    <- fitted(fit)                  # predicted probabilities under observed data
  w_bar    <- mean(p_hat * (1 - p_hat))    # mean logistic weight
  est      <- glm_coef * w_bar             # RD = β * mean(p(1-p))
  est_var  <- (w_bar^2) * coef_var         # delta method: (dRD/dβ)² * Var(β)
  est_se   <- sqrt(est_var)
  
  # estimand
  esd_m <- as.numeric(colMeans(d[2+as.numeric(trt)]))-mean(d$y0)
  
  return(c(est = est, var = est_var, se = est_se, esd_m = esd_m))
  # return(c(est=glm_coef, var=coef_var, se=coef_sd, esd_m=esd_m)) #returns beta
}

################
# g-estimation
g_estimate_binary <- function(data, trt = 1, dataset="ACA") {
  
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
  fit <- glm(y ~ as.factor(t) + as.factor(p) + x, data = d,
             family = binomial(link = "logit"))
  
  d1 <- transform(d, t = factor(trt, levels = levels(as.factor(d$t))))
  d0 <- transform(d, t = factor(0,   levels = levels(as.factor(d$t))))
  
  p1 <- predict(fit, newdata = d1, type = "response")  # P(y=1 | t=trt, x)
  p0 <- predict(fit, newdata = d0, type = "response")  # P(y=1 | t=0,   x)
  
  g_diffmeans <- mean(p1 - p0)
  
  # Variance via delta method (numerical gradient)
  beta_hat <- coef(fit)
  
  g_func <- function(beta) {
    fit_tmp <- fit
    fit_tmp$coefficients <- beta
    p1_tmp <- predict(fit_tmp, newdata = d1, type = "response")
    p0_tmp <- predict(fit_tmp, newdata = d0, type = "response")
    mean(p1_tmp - p0_tmp)
  }
  
  grad   <- numDeriv::grad(g_func, beta_hat)
  Vbeta  <- vcov(fit)
  g_var  <- as.numeric(t(grad) %*% Vbeta %*% grad)
  g_se   <- sqrt(g_var)
  
  # estimand  
  esd_m <- as.numeric(colMeans(d[2+as.numeric(trt)]))-mean(d$y0) 
  
  return(c(est=g_diffmeans, var=g_var, se=g_se, esd_m=esd_m))
}

################
# aipw estimation
aipw_estimate_binary <- function(data, trt = 1, dataset="ACA", population = NULL) {
  
  # Subset to relevant treatments
  cc_periods <- which(as.vector(table(data$t, data$p)[trt+1,])>0)
  
  # data to predict marginal means
  if(dataset=="ACA"){
    d <- subset(data, t %in% c(0,trt) & p %in% cc_periods) # by default ACA data
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
  fit <- glm(y ~ as.factor(t) + as.factor(p) + x, data = d,
             family = binomial(link = "logit"))
  
  # predicted probabilities (type = "response" applies the inverse logit)
  d$mu_hat_0 <- predict(fit, newdata = transform(d, t = factor(0,   levels = levels(d$t))), type = "response")
  d$mu_hat_a <- predict(fit, newdata = transform(d, t = factor(trt, levels = levels(d$t))), type = "response")
  
  
  # If population is specified, subset and stabilize
  if (!is.null(population)) {
    # I_P <- which(d$p %in% period_of_interest) 
    # d_int <- d[I_P,]
    
    if(population=="ACA"){
      d_int <- subset(d, t %in% c(0,trt) & p %in% cc_periods) # by default ACA data
    }
    if(population=="LACA"){
      d_int <- subset(d, t %in% c(0,trt) & p %in% max(cc_periods)) # LACA data
    }
    if(population=="ECE"){
      d_int <- subset(d, p %in% cc_periods) # ECE data
    }
    d_int[which(d_int$t==0),]$pi_hat <-  sum(d_int$t==0)/dim(d_int)[1]
    d_int[which(d_int$t==trt),]$pi_hat <-  sum(d_int$t==trt)/dim(d_int)[1]
    d <- d_int
  }
  
  # Compute individual contributions
  psi_i <- ((d$t == trt) / d$pi_hat) * (d$y - d$mu_hat_a) -
    ((d$t == 0) / d$pi_hat) * (d$y - d$mu_hat_0) +
    (d$mu_hat_a - d$mu_hat_0)
  
  psi_hat <- mean(psi_i)
  
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
dat <- simdata_blocked_unequal(
  mu0 = 0,
  mu1 = 1,
  mu2 = 2,
  beta = 2,
  N1 = 120,
  N2 = 300,
  N3 = 200,
  # alloc1 = c(2/3, 1/3),        # 2:1
  # alloc2 = c(3/6, 2/6, 1/6),   # 3:2:1
  # alloc3 = c(4/5, 1/5),
  alloc1 = rep(1,2)*1/2, # eq alloc
  alloc2 = rep(1,3)*1/3, # eq alloc
  alloc3 = rep(1,2)*1/2, # eq alloc
  lambda = 0.5
)
head(dat)
dat$y <- as.numeric(dat$y>0)
dat$y0 <- as.numeric(dat$y0>0)
dat$y1 <- as.numeric(dat$y1>0)
dat$y2 <- as.numeric(dat$y2>0)
head(dat)

trt=1
data<-dat
dataset="ACA"

# linear model
glmmodel(dat, trt=1)
# 
# g-estimator
gest <- g_estimate_binary(dat, trt=1)
gest
# 
# aipw
aipw_estimate_binary(dat, trt = 1)
aipw_estimate_binary(dat, trt = 2, dataset = "NCC", population = "ACA")

