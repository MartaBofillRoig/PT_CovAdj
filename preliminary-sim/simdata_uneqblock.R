################
# NCC Cov adj
# SIMULATION
# December 2025
################

# Simulation function including unequal allocation and covariate X
# As in the prev version, we keep three-period trial and block randomisation

################
# Aux function: blocked randomisation with unequal allocation
block_randomise <- function(N, trt, alloc, block_size = NULL) {
  
  # alloc = allocation proportions (must sum to 1)
  # trt   = treatment labels
  # N     = total sample size
  
  if (is.null(block_size)) {
    block_size <- sum(alloc / min(alloc))
  }
  
  counts <- alloc * block_size 
  counts <- round(counts)
  
  # one block
  block <- rep(trt, counts)
  
  # number of full blocks
  n_blocks <- N %/% block_size
  
  t <- unlist(replicate(n_blocks, sample(block), simplify = FALSE))
  
  # remainder: sample proportionally WITHOUT replacement
  remainder_size <- N - length(t)
  if (remainder_size > 0) {
    # Create a partial block with proportional allocation
    remainder_counts <- round(alloc * remainder_size)
    
    # Adjust if rounding caused sum != remainder_size
    diff <- remainder_size - sum(remainder_counts)
    if (diff != 0) {
      idx <- which.max(alloc)
      remainder_counts[idx] <- remainder_counts[idx] + diff
    }
    
    remainder_block <- rep(trt, remainder_counts)
    t <- c(t, sample(remainder_block))
  }
  
  return(t)
}

################
# Main simulation function (with blocks)
simdata_blocked_unequal <- function(
    mu0 = 0, mu1 = 0, mu2 = 0,      # treatment main effects
    beta = 1,                        # baseline covariate effect
    gamma1 = 0, gamma2 = 0,          # treatment-covariate interactions
    delta = 0,                       # period-covariate interaction
    tau1 = 0, tau2 = 0,              # treatment-period interactions
    N1 = 200, N2 = 200, N3 = 200,
    alloc1 = c(0.5, 0.5),        # period 1: (0,1)
    alloc2 = c(1/3, 1/3, 1/3),   # period 2: (0,1,2)
    alloc3 = c(0.5, 0.5),        # period 3: (0,2)
    lambda = 0,
    trendp = "stepwise",
    sd=1,
    subgroup=F
) {
  
  stopifnot(abs(sum(alloc1)-1) < 1e-5,
            abs(sum(alloc2)-1) < 1e-5,
            abs(sum(alloc3)-1) < 1e-5)
  
  ## --- blocked randomisation ---
  t1 <- block_randomise(
    N = N1,
    trt = c(0, 1),
    alloc = alloc1
  )
  
  t2 <- block_randomise(
    N = N2,
    trt = c(0, 1, 2),
    alloc = alloc2
  )
  
  t3 <- block_randomise(
    N = N3,
    trt = c(0, 2),
    alloc = alloc3
  )
  
  t <- c(t1, t2, t3)
  p <- factor(c(rep(1, N1), rep(2, N2), rep(3, N3)))
  
  ## --- time trend and baseline covariate ---
  N <- N1 + N2 + N3
  if(trendp == "linear"){
    trend <- (0:(N - 1)) / N * lambda #linear
  }
  if(trendp == "stepwise"){
    cj <- as.numeric(p)
    trend <- as.numeric(lambda)*(cj-1) #stepwise
  }

  x <- rnorm(N, mean = 0, sd = 1)
  
  if(subgroup==T){
    # Two subgroups based on a covariate
    x1 <- rnorm(N, mean = 0, sd = 1)
    x2 <- rnorm(N, mean = 2, sd = 1)
    # Proportion of subgroup per period
    p1_x1 <- 0.5
    p2_x1 <- 0.4
    p3_x1 <- 0.2
    
    x <- c(
      sample(c(sample(x1, round(p1_x1*N1)), sample(x2, N1-round(p1_x1*N1)))),
      sample(c(sample(x1, round(p2_x1*N2)), sample(x2, N2-round(p2_x1*N2)))),
      sample(c(sample(x1, round(p3_x1*N3)), sample(x2, N3-round(p3_x1*N3))))
    )
  }
  ## --- potential outcomes ---
  mu <- c(mu0, mu1, mu2)
  
  ## --- period-specific covariate effect ---
  # gamma1, gamma2 -> treat-covariate interaction
  beta_period <- beta + delta * (as.numeric(p) - 1) #period-covariate interaction
  
  ## --- treatment-period interaction effects ---
  # tau1 and tau2 represent how treatment effects change per period
  period_effect_trt1 <- tau1 * (as.numeric(p) - 1)
  period_effect_trt2 <- tau2 * (as.numeric(p) - 1)
  
  ## --- data generation: potential outcomes ---
  # Treatment 0: baseline (no treatment-period interaction)
  y0 <- mu0 + beta_period * x + trend + rnorm(N, sd = sd)
  
  # Treatment 1: includes treatment-period interaction
  y1 <- (mu1 + period_effect_trt1) + (beta_period + gamma1) * x + trend + rnorm(N, sd = sd)
  
  # Treatment 2: includes treatment-period interaction
  y2 <- (mu2 + period_effect_trt2) + (beta_period + gamma2) * x + trend + rnorm(N, sd = sd)
  
  ## --- observed outcome ---
  y <- (t == 0) * y0 + (t == 1) * y1 + (t == 2) * y2
  
  dat <- data.frame(
    y = y,
    y0 = y0,
    y1 = y1,
    y2 = y2,
    x = x,
    t = factor(t),
    p = p
  )
  
  return(dat)
}

################
# EXAMPLES
# 
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
#   lambda = 0.4
# )
# 
# table(dat$t, dat$p)
# head(dat[,1:5])
# 
# fit <- lm(y ~ t + p + x, data = dat)
# summary(fit)
# 
# mu0 = 0;
# mu1 = 1;
# mu2 = 2;
# beta=2;
# N1 = 120;
# N2 = 300;
# N3 = 200;
# alloc1 = c(2/3, 1/3);        # 2:1
# alloc2 = c(3/6, 2/6, 1/6);   # 3:2:1
# alloc3 = c(4/5, 1/5);
# lambda = 0.4

