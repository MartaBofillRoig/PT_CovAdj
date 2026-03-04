################
# NCC Cov adj
# SMALL SIM
# March 2026
################

rm(list=ls())

setwd("C:/Users/marta.bofill/Dropbox/C5/GitHub/PT_CovAdj/preliminary-sim")
# setwd("C:/Users/marta/Dropbox/C5/GitHub/PT_CovAdj/preliminary-sim")
source("simdata_uneqblock.R")
source("analysismethods.R")

################

library(kableExtra)
library(knitr)
library(dplyr)

################


################
# SIMULATION STUDY
################
set.seed(1225)

# --- define scenarios --- 
scenarios <- list(
  # ===== BASELINE SCENARIOS =====
  list(name = "Equal alloc, no trend (h0)",
       mu0 = 0, mu1 = 0, mu2 = 0,      # treatment main effects
       beta = 0,                        # baseline covariate effect
       gamma1 = 0, gamma2 = 0,          # treatment-covariate interactions
       delta = 0, subgroup = F,         # period-covariate interaction
       tau1 = 0, tau2 = 0,              # treatment-period interactions
       N1 =  100, N2 =  100, N3 =  100,
       alloc1 = c(0.5, 0.5),        # period 1: (0,1)
       alloc2 = c(1/3, 1/3, 1/3),   # period 2: (0,1,2)
       alloc3 = c(0.5, 0.5),        # period 3: (0,2)
       lambda = 0,
       trendp = "stepwise",
       sd=1), 
  
  list(name = "Unequal alloc, no trend (h0)",
       mu0 = 0, mu1 = 0, mu2 = 0,      # treatment main effects
       beta = 0,                        # baseline covariate effect
       gamma1 = 0, gamma2 = 0,          # treatment-covariate interactions
       delta = 0, subgroup = F,         # period-covariate interaction
       tau1 = 0, tau2 = 0,              # treatment-period interactions
       N1 =  100, N2 =  100, N3 =  100,
       alloc1 = c(2, 1)/3,
       alloc2 = c(2, 1, 1)/4, 
       alloc3 = c(2, 1)/3,
       lambda = 0,
       trendp = "stepwise",
       sd=1), 
  
  list(name = "Unequal alloc, no trend (h1)",
       mu0 = 0, mu1 = 1, mu2 = 2,      # treatment main effects
       beta = 0,                        # baseline covariate effect
       gamma1 = 0, gamma2 = 0,          # treatment-covariate interactions
       delta = 0, subgroup = F,         # period-covariate interaction
       tau1 = 0, tau2 = 0,              # treatment-period interactions
       N1 =  100, N2 =  100, N3 =  100,
       alloc1 = c(2, 1)/3,
       alloc2 = c(2, 1, 1)/4, 
       alloc3 = c(2, 1)/3,
       lambda = 0,
       trendp = "stepwise",
       sd=1), 
  
  list(name = "Unequal alloc, trend (h1)",
       mu0 = 0, mu1 = 1, mu2 = 2,      # treatment main effects
       beta = 0,                        # baseline covariate effect
       gamma1 = 0, gamma2 = 0,          # treatment-covariate interactions
       delta = 0, subgroup = F,         # period-covariate interaction
       tau1 = 0, tau2 = 0,              # treatment-period interactions
       N1 =  100, N2 =  100, N3 =  100,
       alloc1 = c(2, 1)/3,
       alloc2 = c(2, 1, 1)/4, 
       alloc3 = c(2, 1)/3,
       lambda = 1,
       trendp = "stepwise",
       sd=1), 
  
  list(name = "Equal alloc, trend, diff sizes N (h1)",
       mu0 = 0, mu1 = 1, mu2 = 2,      # treatment main effects
       beta = 0,                        # baseline covariate effect
       gamma1 = 0, gamma2 = 0,          # treatment-covariate interactions
       delta = 0, subgroup = F,         # period-covariate interaction
       tau1 = 0, tau2 = 0,              # treatment-period interactions
       N1 =  50, N2 =  150, N3 =  100,
       alloc1 = c(0.5, 0.5),        # period 1: (0,1)
       alloc2 = c(1/3, 1/3, 1/3),   # period 2: (0,1,2)
       alloc3 = c(0.5, 0.5),        # period 3: (0,2)
       lambda = 1,
       trendp = "stepwise",
       sd=1), 
  
  # ===== COVARIATE EFFECTS =====
  list(name = "Equal alloc, no trend, with cov effect (h0)",
       mu0 = 0, mu1 = 0, mu2 = 0,      # treatment main effects
       beta = 2,                        # baseline covariate effect
       gamma1 = 0, gamma2 = 0,          # treatment-covariate interactions
       delta = 0, subgroup = F,         # period-covariate interaction
       tau1 = 0, tau2 = 0,              # treatment-period interactions
       N1 =  100, N2 =  100, N3 =  100,
       alloc1 = c(0.5, 0.5),        # period 1: (0,1)
       alloc2 = c(1/3, 1/3, 1/3),   # period 2: (0,1,2)
       alloc3 = c(0.5, 0.5),        # period 3: (0,2)
       lambda = 0,
       trendp = "stepwise",
       sd=1), 
  
  list(name = "Unequal alloc, trend, with cov effect (h1)",
       mu0 = 0, mu1 = 1, mu2 = 2,      # treatment main effects
       beta = 2,                        # baseline covariate effect
       gamma1 = 0, gamma2 = 0,          # treatment-covariate interactions
       delta = 0, subgroup = F,         # period-covariate interaction
       tau1 = 0, tau2 = 0,              # treatment-period interactions
       N1 =  100, N2 =  100, N3 =  100,
       alloc1 = c(2, 1)/3,
       alloc2 = c(2, 1, 1)/4, 
       alloc3 = c(2, 1)/3,
       lambda = 1,
       trendp = "stepwise",
       sd=1), 
  
  list(name = "Equal alloc, trend, diff sizes N, with cov effect (h1)",
       mu0 = 0, mu1 = 1, mu2 = 2,      # treatment main effects
       beta = 2,                        # baseline covariate effect
       gamma1 = 0, gamma2 = 0,          # treatment-covariate interactions
       delta = 0, subgroup = F,         # period-covariate interaction
       tau1 = 0, tau2 = 0,              # treatment-period interactions
       N1 =  50, N2 =  150, N3 =  100,
       alloc1 = c(0.5, 0.5),        # period 1: (0,1)
       alloc2 = c(1/3, 1/3, 1/3),   # period 2: (0,1,2)
       alloc3 = c(0.5, 0.5),        # period 3: (0,2)
       lambda = 1,
       trendp = "stepwise",
       sd=1), 
  
  # ===== COVARIATE-PERIOD INTERACTIONS =====
  list(name = "Equal alloc, no trend, cov-period interaction (h0)",
       mu0 = 0, mu1 = 0, mu2 = 0,      # treatment main effects
       beta = 2,                        # baseline covariate effect
       gamma1 = 0, gamma2 = 0,          # treatment-covariate interactions
       delta = 0, subgroup = T,         # period-covariate interaction
       tau1 = 0, tau2 = 0,              # treatment-period interactions
       N1 =  100, N2 =  100, N3 =  100,
       alloc1 = c(0.5, 0.5),        # period 1: (0,1)
       alloc2 = c(1/3, 1/3, 1/3),   # period 2: (0,1,2)
       alloc3 = c(0.5, 0.5),        # period 3: (0,2)
       lambda = 0,
       trendp = "stepwise",
       sd=1), 
  
  list(name = "Equal alloc, no trend, cov-period interaction (h1)",
       mu0 = 0, mu1 = 1, mu2 = 2,      # treatment main effects
       beta = 2,                        # baseline covariate effect
       gamma1 = 0, gamma2 = 0,          # treatment-covariate interactions
       delta = 0, subgroup = T,         # period-covariate interaction
       tau1 = 0, tau2 = 0,              # treatment-period interactions
       N1 =  100, N2 =  100, N3 =  100,
       alloc1 = c(0.5, 0.5),        # period 1: (0,1)
       alloc2 = c(1/3, 1/3, 1/3),   # period 2: (0,1,2)
       alloc3 = c(0.5, 0.5),        # period 3: (0,2)
       lambda = 0,
       trendp = "stepwise",
       sd=1), 
  
  list(name = "Unequal alloc, no trend, cov-period interaction (h0)",
       mu0 = 0, mu1 = 0, mu2 = 0,      # treatment main effects
       beta = 2,                        # baseline covariate effect
       gamma1 = 0, gamma2 = 0,          # treatment-covariate interactions
       delta = 0, subgroup = T,         # period-covariate interaction
       tau1 = 0, tau2 = 0,              # treatment-period interactions
       N1 =  100, N2 =  100, N3 =  100,
       alloc1 = c(2, 1)/3,
       alloc2 = c(2, 1, 1)/4, 
       alloc3 = c(2, 1)/3,
       lambda = 0,
       trendp = "stepwise",
       sd=1), 
  
  list(name = "Unequal alloc, no trend, cov-period interaction (h1)",
       mu0 = 0, mu1 = 1, mu2 = 2,      # treatment main effects
       beta = 2,                        # baseline covariate effect
       gamma1 = 2, gamma2 = 0,          # treatment-covariate interactions
       delta = 0, subgroup = T,         # period-covariate interaction
       tau1 = 0, tau2 = 0,              # treatment-period interactions
       N1 =  100, N2 =  100, N3 =  100,
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
       N1 =  100, N2 =  100, N3 =  100,
       alloc1 = c(1, 1)/2,
       alloc2 = c(1, 1, 1)/3,
       alloc3 = c(1, 1)/2,
       lambda = 1,
       gamma1 = 0, gamma2 = 0,
       delta = 0, subgroup = F,         # period-covariate interaction
       trendp = "stepwise", sd=1),
  
  list(name = "Equal alloc, treat-period interaction (h1)",
       mu0 = 0, mu1 = 1, mu2 = 2, beta = 0, 
       tau1 = 3, tau2 = 0,   
       N1 =  100, N2 =  100, N3 =  100,
       alloc1 = c(1, 1)/2,
       alloc2 = c(1, 1, 1)/3,
       alloc3 = c(1, 1)/2,
       lambda = 1,
       gamma1 = 0, gamma2 = 0,
       delta = 0, subgroup = F,         # period-covariate interaction
       trendp = "stepwise", sd=1) 
)

# --- contrasts to compute: 1 vs 0 and 2 vs 0 ---
contrasts <- list(c(1,0), c(2,0))

# --- helper function to add results ---
add_result <- function(results, sc, trt, method, dataset, population, out) {
  rbind(
    results,
    data.frame(
      Scenario = sc$name,
      Contrast = trt,
      Method   = method,
      Dataset  = dataset,
      Pop      = population,
      Est      = out["est"],
      Var      = out["var"],
      SE       = out["se"],
      Estimandm = out["esd_m"], 
      stringsAsFactors = FALSE
    )
  )
}

# --- function to run single simulation ---
run_single_sim <- function() {
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
      delta = sc$delta, subgroup = sc$subgroup,
      trendp = sc$trendp, sd=sc$sd
    )
    
    for (ctr in contrasts) {
      
      trt <- ctr[1]
      
      ## --- Linear model ---
      results <- add_result(results, sc, trt, "LM", "ACA", "ACA",
                            lmmodel(dat, trt = trt, dataset = "ACA"))
      
      results <- add_result(results, sc, trt, "LM", "ECE", "ECE",
                            lmmodel(dat, trt = trt, dataset = "ECE"))
      
      results <- add_result(results, sc, trt, "LM", "NCC", "NCC",
                            lmmodel(dat, trt = trt, dataset = "NCC"))
      
      ## --- G-estimation ---
      results <- add_result(results, sc, trt, "G", "ACA", "ACA",
                            g_estimate(dat, trt = trt, dataset = "ACA"))
      
      results <- add_result(results, sc, trt, "G", "ECE", "ECE",
                            g_estimate(dat, trt = trt, dataset = "ECE"))
      
      results <- add_result(results, sc, trt, "G", "NCC", "NCC", 
                            g_estimate(dat, trt = trt, dataset = "NCC"))
      
      ## --- AIPW ---
      results <- add_result(results, sc, trt, "AIPW", "ACA", "ACA",
                            aipw_estimate(dat, trt = trt, dataset="ACA", population = NULL)) 
      
      results <- add_result(results, sc, trt, "AIPW", "ECE", "ACA",
                            aipw_estimate(dat, trt = trt, dataset="ECE", population = "ACA")) 
      
      results <- add_result(results, sc, trt, "AIPW", "ECE", "ECE",
                            aipw_estimate(dat, trt = trt, dataset="ECE", population = NULL)) 
      
      results <- add_result(results, sc, trt, "AIPW", "NCC", "ECE",
                            aipw_estimate(dat, trt = trt, dataset="NCC", population = "ECE")) 
      
      results <- add_result(results, sc, trt, "AIPW", "ECE", "LACA",
                            aipw_estimate(dat, trt = trt, dataset="ECE", population = "LACA"))
    }
  }
  
  return(results)
}

# --- run full simulation study: run n_sims trials ---
n_sims <- 10000

cat("Starting simulation study with", n_sims, "iterations...\n")
start_time <- Sys.time()

# Initialize list to store results from each simulation
all_results <- list()

# Progress tracking
progress_interval <- 1000

for (i in 1:n_sims) {
  
  # Progress reporting
  if (i %% progress_interval == 0) {
    elapsed_time <- Sys.time() - start_time
    cat("Completed", i, "simulations. Elapsed time:", 
        round(elapsed_time, 2), attr(elapsed_time, "units"), "\n")
  }
  
  # Run single simulation (with error handling)
  tryCatch({
    sim_result <- run_single_sim()
    all_results[[i]] <- sim_result
  }, error = function(e) {
    cat("Error in simulation", i, ":", e$message, "\n")
    all_results[[i]] <- NA
  })
}

################
# RESULTS

# --- combine all results ---
# Remove failed simulations
valid_results <- all_results[!is.na(all_results)]
cat("Successfully completed", length(valid_results), "out of", n_sims, "simulations\n")

# Combine all results into one data frame
combined_results <- do.call(rbind, valid_results)

# Add simulation number
combined_results$sim_number <- rep(1:length(valid_results), 
                                   each = nrow(valid_results[[1]])) 

# --- summarising the simulation iterations ---
combined_results <- combined_results %>%
  mutate(Scenario = factor(Scenario, levels = unique(Scenario)))

summary_stats <- combined_results %>%
  group_by(Scenario, Contrast, Method, Dataset, Pop) %>%
  summarise(
    # Sample size
    n_sims = n(),
    
    # Estimate statistics
    mean_Est = mean(Est, na.rm = TRUE),
    sd_Est = sd(Est, na.rm = TRUE),
    median_Est = median(Est, na.rm = TRUE),
    
    # Variance statistics
    mean_Var = mean(Var, na.rm = TRUE),
    sd_Var = sd(Var, na.rm = TRUE),
    median_Var = median(Var, na.rm = TRUE),
    
    # Standard Error statistics
    mean_SE = mean(SE, na.rm = TRUE),
    sd_SE = sd(SE, na.rm = TRUE),
    median_SE = median(SE, na.rm = TRUE),
    
    # Additional useful statistics
    empirical_Var = var(Est, na.rm = TRUE), # empirical standard error
    mean_Var_vs_empirical_Var = mean(Var, na.rm = TRUE) / var(Est, na.rm = TRUE),
    
    .groups = 'drop'
  ) 

# --- Table in latex --- 
# Create line separator pattern: "" for most rows, "\\addlinespace" every 16 rows
n_rows <- nrow(summary_stats)
linesep_pattern <- rep("", n_rows)
linesep_pattern[seq(16, n_rows, by = 16)] <- "\\addlinespace"

# Create the table
kable(summary_stats[,-c(5)], 
      format = "latex", 
      booktabs = TRUE, 
      digits = 3,
      linesep = linesep_pattern)

# --- display results ---
end_time <- Sys.time()
total_time <- end_time - start_time

cat("\nSimulation completed!\n")
cat("Total time:", round(total_time, 2), attr(total_time, "units"), "\n")
cat("Results summary:\n")
summary_stats

# --- save workspace ---  
save(combined_results, summary_stats, scenarios, n_sims,
     file = "simulation_results.RData") 


