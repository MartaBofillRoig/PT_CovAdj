# aipw estimation
aipw_estimate <- function(data, trt = 1, dataset="ACA") {
  
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
  
  # Outcome model (Here we could use another population if we would want)
  fit <- lm(y ~ as.factor(t) + as.factor(p) + x, data = d)
  
  # Predicted outcomes for treatment 0 and treatment of interest
  d$mu_hat_0 <- predict(fit, newdata = transform(d, t = factor(0, levels = levels(d$t))))
  d$mu_hat_a <- predict(fit, newdata = transform(d, t = factor(trt, levels = levels(d$t))))
  
  # Compute individual contributions
  psi_i <- ((d$t == trt) / d$pi_hat) * (d$y - d$mu_hat_a) -
    ((d$t == 0) / d$pi_hat) * (d$y - d$mu_hat_0) +
    (d$mu_hat_a - d$mu_hat_0)
  
  psi_hat <- mean(psi_i)
  
  n <- length(psi_i)
  
  # Empirical variance of influence function
  var_psi <- var(psi_i)
  
  # Variance and SE of the estimator
  var_hat <- var_psi/n
  se_hat  <- sqrt(var_hat)
  
  return(c(est=psi_hat, var=var_hat, se=se_hat))
}