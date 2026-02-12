################
# NCC Cov adj
# SMALL SIM
# Feb 2026
################

rm(list=ls())

setwd("C:/Users/marta.bofill/Dropbox/C5/RESEARCH/RCode/2026-02-NCC-covadj")
# setwd("C:/Users/marta/Dropbox/C5/RESEARCH/RCode/2026-02-NCC-covadj")
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
  list(name = "Equal alloc, no trend", 
       mu0 = 0, mu1 = 1, mu2 = 2, beta = 1, 
       N1 = 100, N2 = 100, N3 = 100,
       alloc1 = c(1, 1)/2,        # P1: 50 control, 50 trt1
       alloc2 = c(1, 1, 1)/3,     # P2: 33 control, 33 trt1, 33 trt2
       alloc3 = c(1, 1)/2,        # P3: 50 control, 50 trt2
       lambda = 0),
  # Total: Control=133, Trt1=83, Trt2=83
  
  list(name = "Unequal alloc, no trend",
       mu0 = 0, mu1 = 1, mu2 = 2, beta = 1, 
       N1 = 150, N2 = 300, N3 = 150,
       alloc1 = c(2, 1)/3,        # P1: 100 control, 50 trt1
       alloc2 = c(2, 1, 1)/4,     # P2: 150 control, 75 trt1, 75 trt2
       alloc3 = c(2, 1)/3,        # P3: 100 control, 50 trt2
       lambda = 0),
  # Total: Control=350, Trt1=125, Trt2=125
  
  # ===== TREND SCENARIOS =====
  list(name = "Equal alloc, strong trend",
       mu0 = 0, mu1 = 1, mu2 = 2, beta = 1, 
       N1 = 100, N2 = 100, N3 = 100,
       alloc1 = c(1, 1)/2,
       alloc2 = c(1, 1, 1)/3,
       alloc3 = c(1, 1)/2,
       lambda = 1),
  # Total: Control=133, Trt1=83, Trt2=83
  
  list(name = "Unequal alloc, trend + large beta",
       mu0 = 0, mu1 = 1, mu2 = 2, beta = 2, 
       N1 = 150, N2 = 300, N3 = 150,
       alloc1 = c(2, 1)/3,
       alloc2 = c(2, 1, 1)/4,
       alloc3 = c(2, 1)/3,
       lambda = 1),
  # Total: Control=350, Trt1=125, Trt2=125
  
  # ===== TREATMENT-COVARIATE INTERACTION =====
  list(name = "Equal alloc, trend, treat-cov interaction", 
       mu0 = 0, mu1 = 1, mu2 = 2, beta = 2, gamma1 = 2, gamma2 = 0,
       N1 = 100, N2 = 100, N3 = 100,
       alloc1 = c(1, 1)/2,
       alloc2 = c(1, 1, 1)/3,
       alloc3 = c(1, 1)/2,
       lambda = 1), 
  # Total: Control=133, Trt1=83, Trt2=83
  
  list(name = "Unequal alloc, trend, treat-cov interaction", 
       mu0 = 0, mu1 = 1, mu2 = 2, beta = 2, gamma1 = 2, gamma2 = 1,
       N1 = 150, N2 = 300, N3 = 150,
       alloc1 = c(2, 1)/3,
       alloc2 = c(2, 1, 1)/4,
       alloc3 = c(2, 1)/3,
       lambda = 1),
  # Total: Control=350, Trt1=125, Trt2=125
  
  # ===== PERIOD-COVARIATE INTERACTION =====
  list(name = "Equal alloc, trend, period-cov interaction", 
       mu0 = 0, mu1 = 1, mu2 = 2, beta = 2, delta = 2,
       N1 = 100, N2 = 100, N3 = 100,
       alloc1 = c(1, 1)/2,
       alloc2 = c(1, 1, 1)/3,
       alloc3 = c(1, 1)/2,
       lambda = 1),
  # Total: Control=133, Trt1=83, Trt2=83
  
  list(name = "Unequal alloc, trend, period-cov interaction", 
       mu0 = 0, mu1 = 1, mu2 = 2, beta = 2, delta = 2,
       N1 = 150, N2 = 300, N3 = 150,
       alloc1 = c(2, 1)/3,
       alloc2 = c(2, 1, 1)/4,
       alloc3 = c(2, 1)/3,
       lambda = 1),
  # Total: Control=350, Trt1=125, Trt2=125
  
  # ===== TREATMENT-PERIOD INTERACTIONS (NEW) =====
  list(name = "Equal alloc, treat-period interaction (waning trt1)",
       mu0 = 0, mu1 = 1, mu2 = 2, beta = 1, 
       tau1 = -1, tau2 = 0,  # Trt1 loses 1 unit per period
       N1 = 100, N2 = 100, N3 = 100,
       alloc1 = c(1, 1)/2,
       alloc2 = c(1, 1, 1)/3,
       alloc3 = c(1, 1)/2,
       lambda = 1),
  # Total: Control=133, Trt1=83, Trt2=83
  # Effect trajectory: Trt1 goes 3→2→1, Trt2 stays at 2
  
  list(name = "Unequal alloc, treat-period interaction (delayed trt2)",
       mu0 = 0, mu1 = 1, mu2 = 2, beta = 1,
       tau1 = 0, tau2 = 1,  # Trt2 gains 1 unit per period
       N1 = 150, N2 = 300, N3 = 150,
       alloc1 = c(2, 1)/3,
       alloc2 = c(2, 1, 1)/4,
       alloc3 = c(2, 1)/3,
       lambda = 1),
  # Total: Control=350, Trt1=125, Trt2=125
  # Effect trajectory: Trt1 stays at 2, Trt2 goes 0→1→2
  
  list(name = "Equal alloc, opposing treat-period effects",
       mu0 = 0, mu1 = 1, mu2 = 2, beta = 1,
       tau1 = -1, tau2 = 1,  # Trt1 wanes, Trt2 improves (strong effects)
       N1 = 100, N2 = 100, N3 = 100,
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
       N1 = 150, N2 = 300, N3 = 150,
       alloc1 = c(2, 1)/3,
       alloc2 = c(2, 1, 1)/4,
       alloc3 = c(2, 1)/3,
       lambda = 2),  # Even stronger time trend
  # Total: Control=350, Trt1=125, Trt2=125 
  
  # ===== ADDITIONAL EXTREME SCENARIOS =====
  list(name = "Equal alloc, extreme covariate interactions only",
       mu0 = 0, mu1 = 1, mu2 = 2, beta = 3,
       gamma1 = 3, gamma2 = -3,  # Very strong opposing interactions
       N1 = 100, N2 = 100, N3 = 100,
       alloc1 = c(1, 1)/2,
       alloc2 = c(1, 1, 1)/3,
       alloc3 = c(1, 1)/2,
       lambda = 0),
  # Total: Control=133, Trt1=83, Trt2=83 
  
  list(name = "Unequal alloc, extreme period effects",
       mu0 = 0, mu1 = 1, mu2 = 2, beta = 1,
       tau1 = 2, tau2 = 2,     # Both treatments improve dramatically
       delta = 3,              # Period-covariate also very strong
       N1 = 150, N2 = 300, N3 = 150,
       alloc1 = c(2, 1)/3,
       alloc2 = c(2, 1, 1)/4,
       alloc3 = c(2, 1)/3,
       lambda = 2)
  # Total: Control=350, Trt1=125, Trt2=125 
)

# --- contrasts to compute: 1 vs 0 and 2 vs 0 ---
contrasts <- list(c(1,0), c(2,0))

# --- helper function to add results ---
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

# --- function to run single simulation ---
run_single_sim <- function() {
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
  group_by(Scenario, Contrast, Method, Dataset) %>%
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


