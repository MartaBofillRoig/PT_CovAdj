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
  list(name = "Equal alloc, no trend", 
       mu0 = 0, mu1 = 1, mu2 = 2, beta = 1, 
       N1 = n* 100, N2 = n* 100, N3 = n* 100,
       alloc1 = c(1, 1)/2,        # P1: 50 control, 50 trt1
       alloc2 = c(1, 1, 1)/3,     # P2: 33 control, 33 trt1, 33 trt2
       alloc3 = c(1, 1)/2,        # P3: 50 control, 50 trt2
       lambda = 0),
  # Total: Control=133, Trt1=83, Trt2=83
  
  list(name = "Unequal alloc, no trend",
       mu0 = 0, mu1 = 1, mu2 = 2, beta = 1, 
       N1 = n* 150, N2 = n* 300, N3 = n* 150,
       alloc1 = c(2, 1)/3,        # P1: 100 control, 50 trt1
       alloc2 = c(2, 1, 1)/4,     # P2: 150 control, 75 trt1, 75 trt2
       alloc3 = c(2, 1)/3,        # P3: 100 control, 50 trt2
       lambda = 0),
  # Total: Control=350, Trt1=125, Trt2=125
  
  # ===== TREND SCENARIOS =====
  list(name = "Equal alloc, strong trend",
       mu0 = 0, mu1 = 1, mu2 = 2, beta = 1, 
       N1 = n* 100, N2 = n* 100, N3 = n* 100,
       alloc1 = c(1, 1)/2,
       alloc2 = c(1, 1, 1)/3,
       alloc3 = c(1, 1)/2,
       lambda = 1),
  # Total: Control=133, Trt1=83, Trt2=83
  
  list(name = "Unequal alloc, trend + large beta",
       mu0 = 0, mu1 = 1, mu2 = 2, beta = 2, 
       N1 = n* 150, N2 = n* 300, N3 = n* 150,
       alloc1 = c(2, 1)/3,
       alloc2 = c(2, 1, 1)/4,
       alloc3 = c(2, 1)/3,
       lambda = 1),
  # Total: Control=350, Trt1=125, Trt2=125
  
  # ===== TREATMENT-COVARIATE INTERACTION =====
  list(name = "Equal alloc, trend, treat-cov interaction", 
       mu0 = 0, mu1 = 1, mu2 = 2, beta = 2, gamma1 = 2, gamma2 = 0,
       N1 = n* 100, N2 = n* 100, N3 = n* 100,
       alloc1 = c(1, 1)/2,
       alloc2 = c(1, 1, 1)/3,
       alloc3 = c(1, 1)/2,
       lambda = 1), 
  # Total: Control=133, Trt1=83, Trt2=83
  
  list(name = "Unequal alloc, trend, treat-cov interaction", 
       mu0 = 0, mu1 = 1, mu2 = 2, beta = 2, gamma1 = 2, gamma2 = 1,
       N1 = n* 150, N2 = n* 300, N3 = n* 150,
       alloc1 = c(2, 1)/3,
       alloc2 = c(2, 1, 1)/4,
       alloc3 = c(2, 1)/3,
       lambda = 1),
  # Total: Control=350, Trt1=125, Trt2=125
  
  # ===== PERIOD-COVARIATE INTERACTION =====
  list(name = "Equal alloc, trend, period-cov interaction", 
       mu0 = 0, mu1 = 1, mu2 = 2, beta = 2, delta = 2,
       N1 = n* 100, N2 = n* 100, N3 = n* 100,
       alloc1 = c(1, 1)/2,
       alloc2 = c(1, 1, 1)/3,
       alloc3 = c(1, 1)/2,
       lambda = 1),
  # Total: Control=133, Trt1=83, Trt2=83
  
  list(name = "Unequal alloc, trend, period-cov interaction", 
       mu0 = 0, mu1 = 1, mu2 = 2, beta = 2, delta = 2,
       N1 = n* 150, N2 = n* 300, N3 = n* 150,
       alloc1 = c(2, 1)/3,
       alloc2 = c(2, 1, 1)/4,
       alloc3 = c(2, 1)/3,
       lambda = 1),
  # Total: Control=350, Trt1=125, Trt2=125
  
  # ===== TREATMENT-PERIOD INTERACTIONS (NEW) =====
  list(name = "Equal alloc, treat-period interaction (waning trt1)",
       mu0 = 0, mu1 = 1, mu2 = 2, beta = 1, 
       tau1 = -1, tau2 = 0,  # Trt1 loses 1 unit per period
       N1 = n* 100, N2 = n* 100, N3 = n* 100,
       alloc1 = c(1, 1)/2,
       alloc2 = c(1, 1, 1)/3,
       alloc3 = c(1, 1)/2,
       lambda = 1),
  # Total: Control=133, Trt1=83, Trt2=83
  # Effect trajectory: Trt1 goes 3→2→1, Trt2 stays at 2
  
  list(name = "Unequal alloc, treat-period interaction (delayed trt2)",
       mu0 = 0, mu1 = 1, mu2 = 2, beta = 1,
       tau1 = 0, tau2 = 1,  # Trt2 gains 1 unit per period
       N1 = n* 150, N2 = n* 300, N3 = n* 150,
       alloc1 = c(2, 1)/3,
       alloc2 = c(2, 1, 1)/4,
       alloc3 = c(2, 1)/3,
       lambda = 1),
  # Total: Control=350, Trt1=125, Trt2=125
  # Effect trajectory: Trt1 stays at 2, Trt2 goes 0→1→2
  
  list(name = "Equal alloc, opposing treat-period effects",
       mu0 = 0, mu1 = 1, mu2 = 2, beta = 1,
       tau1 = -1, tau2 = 1,  # Trt1 wanes, Trt2 improves (strong effects)
       N1 = n* 100, N2 = n* 100, N3 = n* 100,
       alloc1 = c(1, 1)/2,
       alloc2 = c(1, 1, 1)/3,
       alloc3 = c(1, 1)/2,
       lambda = 1),
  # Total: Control=133, Trt1=83, Trt2=83
  # Effect trajectory: Trt1: 3→2→1, Trt2: 0→1→2 (crossing effects!)
  
  list(name = "Unequal alloc, extreme all interactions",
       mu0 = 0, mu1 = 1, mu2 = 2, beta = 2,
       gamma1 = 2, gamma2 = -2,     # Strong opposing covariate effects
       tau1 = 1, tau2 = -1,          # Opposing period trends
       delta = 2,                    # Strong period-covariate interaction
       N1 = n* 150, N2 = n* 300, N3 = n* 150,
       alloc1 = c(2, 1)/3,
       alloc2 = c(2, 1, 1)/4,
       alloc3 = c(2, 1)/3,
       lambda = 2),  # Even stronger time trend
  # Total: Control=350, Trt1=125, Trt2=125 
  
  # ===== ADDITIONAL EXTREME SCENARIOS =====
  list(name = "Equal alloc, extreme covariate interactions only",
       mu0 = 0, mu1 = 1, mu2 = 2, beta = 3,
       gamma1 = 3, gamma2 = -3,  # Very strong opposing interactions
       N1 = n* 100, N2 = n* 100, N3 = n* 100,
       alloc1 = c(1, 1)/2,
       alloc2 = c(1, 1, 1)/3,
       alloc3 = c(1, 1)/2,
       lambda = 0),
  # Total: Control=133, Trt1=83, Trt2=83 
  
  list(name = "Unequal alloc, extreme period effects",
       mu0 = 0, mu1 = 1, mu2 = 2, beta = 1,
       tau1 = 2, tau2 = 2,     # Both treatments improve dramatically
       delta = 3,              # Period-covariate also very strong
       N1 = n* 150, N2 = n* 300, N3 = n* 150,
       alloc1 = c(2, 1)/3,
       alloc2 = c(2, 1, 1)/4,
       alloc3 = c(2, 1)/3,
       lambda = 2)
  # Total: Control=350, Trt1=125, Trt2=125 
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
    beta = sc$beta,
    N1 = sc$N1, N2 = sc$N2, N3 = sc$N3,
    alloc1 = sc$alloc1, alloc2 = sc$alloc2, alloc3 = sc$alloc3,
    lambda = sc$lambda
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


