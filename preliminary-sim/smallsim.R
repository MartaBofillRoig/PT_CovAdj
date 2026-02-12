################
# NCC Cov adj
# SMALL SIM
# Feb 2026
################

rm(list=ls())

setwd("C:/Users/marta/Dropbox/C5/GitHub/PT_CovAdj/preliminary-sim")
# setwd("C:/Users/marta.bofill/Dropbox/C5/RESEARCH/RCode/2026-02-NCC-covadj")
# setwd("C:/Users/marta/Dropbox/C5/RESEARCH/RCode/2026-02-NCC-covadj")

################

source("simdata_uneqblock.R")
source("analysismethods.R")
# --- helper functions ---
# g_estimate(data, trt=1)
# aipw_estimate(data, trt=1)
# simdata_blocked_unequal(...) from your previous setup

################

library(kableExtra)
library(knitr)
library(dplyr)

################
# LOOP SIM 1 trial

set.seed(1225)
n <- 1000
# --- define scenarios --- 
scenarios <- list(
  # ===== BASELINE SCENARIOS =====
  list(name = "Equal alloc, no trend (h0)",
       mu0 = 0, mu1 = 0, mu2 = 0,      # treatment main effects
       beta = 1,                        # baseline covariate effect
       gamma1 = 0, gamma2 = 0,          # treatment-covariate interactions
       delta = 0,                       # period-covariate interaction
       tau1 = 0, tau2 = 0,              # treatment-period interactions
       N1 = n* 100, N2 = n* 100, N3 = n* 100,
       alloc1 = c(0.5, 0.5),        # period 1: (0,1)
       alloc2 = c(1/3, 1/3, 1/3),   # period 2: (0,1,2)
       alloc3 = c(0.5, 0.5),        # period 3: (0,2)
       lambda = 0,
       trendp = "stepwise",
       sd=1), 
  
  list(name = "Unequal alloc, no trend (h0)",
       mu0 = 0, mu1 = 0, mu2 = 0,      # treatment main effects
       beta = 1,                        # baseline covariate effect
       gamma1 = 0, gamma2 = 0,          # treatment-covariate interactions
       delta = 0,                       # period-covariate interaction
       tau1 = 0, tau2 = 0,              # treatment-period interactions
       N1 = n* 100, N2 = n* 100, N3 = n* 100,
       alloc1 = c(2, 1)/3,
       alloc2 = c(2, 1, 1)/4, 
       alloc3 = c(2, 1)/3,
       lambda = 0,
       trendp = "stepwise",
       sd=1), 
  
  list(name = "Unequal alloc, no trend (h1)",
       mu0 = 0, mu1 = 1, mu2 = 2,      # treatment main effects
       beta = 1,                        # baseline covariate effect
       gamma1 = 0, gamma2 = 0,          # treatment-covariate interactions
       delta = 0,                       # period-covariate interaction
       tau1 = 0, tau2 = 0,              # treatment-period interactions
       N1 = n* 100, N2 = n* 100, N3 = n* 100,
       alloc1 = c(2, 1)/3,
       alloc2 = c(2, 1, 1)/4, 
       alloc3 = c(2, 1)/3,
       lambda = 0,
       trendp = "stepwise",
       sd=1), 
  
  # ===== TREATMENT-PERIOD INTERACTIONS =====
  list(name = "Equal alloc, treat-period interaction (h0)",
       mu0 = 0, mu1 = 0, mu2 = 0, beta = 0, 
       tau1 = 3, tau2 = 0,   
       N1 = n* 100, N2 = n* 100, N3 = n* 100,
       alloc1 = c(1, 1)/2,
       alloc2 = c(1, 1, 1)/3,
       alloc3 = c(1, 1)/2,
       lambda = 1,
       gamma1 = 0, gamma2 = 0,
       delta = 0, trendp = "stepwise", sd=1),
  
  list(name = "Equal alloc, treat-period interaction (h1)",
       mu0 = 0, mu1 = 1, mu2 = 2, beta = 0, 
       tau1 = 3, tau2 = 0,   
       N1 = n* 100, N2 = n* 100, N3 = n* 100,
       alloc1 = c(1, 1)/2,
       alloc2 = c(1, 1, 1)/3,
       alloc3 = c(1, 1)/2,
       lambda = 1,
       gamma1 = 0, gamma2 = 0,
       delta = 0, trendp = "stepwise", sd=1) 
)

# --- contrasts to compute: 1 vs 0 and 2 vs 0 ---
contrasts <- list(c(1,0), c(2,0))

# --- function to save results in dataframe ---
add_result <- function(results, sc, trt, method, dataset, out) {
  rbind(
    results,
    data.frame(
      Scenario = sc$name,
      Contrast = trt,
      Method   = method,
      Dataset  = dataset,
      Est      = out["est"],
      Var      = out["var"],
      SE       = out["se"],
      stringsAsFactors = FALSE
    )
  )
}

# --- run simulation: simulation of 1 trial ---
results <- data.frame()

for (sc in scenarios) {
  
  dat <- simdata_blocked_unequal(
    mu0 = sc$mu0, mu1 = sc$mu1, mu2 = sc$mu2,
    tau1 = sc$tau1, tau2 = sc$tau2,   
    beta = sc$beta,
    N1 = sc$N1, N2 = sc$N2, N3 = sc$N3,
    alloc1 = sc$alloc1, alloc2 = sc$alloc2, alloc3 = sc$alloc3,
    lambda = sc$lambda,
    gamma1 = sc$gamma1, gamma2 = sc$gamma2,
    delta = sc$delta, trendp = sc$trendp, sd=sc$sd
  )
  
  for (ctr in contrasts) {
    
    trt <- ctr[1]
    
    ## --- Linear model ---
    results <- add_result(results, sc, trt, "LM", "ACA",
                          lmmodel(dat, trt = trt, dataset = "ACA"))
    
    results <- add_result(results, sc, trt, "LM", "ECE",
                          lmmodel(dat, trt = trt, dataset = "ECE"))
    
    results <- add_result(results, sc, trt, "LM", "NCC",
                          lmmodel(dat, trt = trt, dataset = "NCC"))
    
    ## --- G-estimation ---
    results <- add_result(results, sc, trt, "G", "ACA",
                          g_estimate(dat, trt = trt, dataset = "ACA"))
    
    results <- add_result(results, sc, trt, "G", "ECE",
                          g_estimate(dat, trt = trt, dataset = "ECE"))
    
    results <- add_result(results, sc, trt, "G", "NCC",
                          g_estimate(dat, trt = trt, dataset = "NCC"))
    
    ## --- AIPW ---
    results <- add_result(results, sc, trt, "AIPW", "ACA",
                          aipw_estimate(dat, trt = trt))
    
    results <- add_result(results, sc, trt, "AIPW", "LACA",
                          aipw_estimate(dat, trt = trt,
                                        period_of_interest = trt + 1))
  }
}

################
# RESULTS
results 

kable(results, format = "latex", booktabs = TRUE, digits=3) 
# kable(results,  booktabs = TRUE, digits=3) 


