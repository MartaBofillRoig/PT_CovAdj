### "Equal alloc, no trend"

mu0 = 0; mu1 = 1; mu2 = 2; beta = 1; 
N1 = 100000; N2 = 100000; N3 = 100000;
alloc1 = c(1, 1)/2        # P1: 50 control, 50 trt1
alloc2 = c(1, 1, 1)/3     # P2: 33 control, 33 trt1, 33 trt2
alloc3 = c(1, 1)/2        # P3: 50 control, 50 trt2
lambda = 0

set.seed(12345)
data <- simdata_blocked_unequal(
  mu0 = mu0, mu1 = mu1, mu2 = mu2,
  beta = beta,
  N1 = N1, N2 = N2, N3 = N3,
  alloc1 = alloc1, alloc2 = alloc2, alloc3 = alloc3,
  lambda = lambda
)

# Treatment 1
trt = 1
cc_periods <- which(as.vector(table(data$t, data$p)[trt+1,])>0)

# Analysis data 
data_ACA1 <- subset(data, t %in% c(0,trt) & p %in% cc_periods)  # ACA data #WRONG, no?
data_LACA1 <-subset(data, t %in% c(0,trt) & p %in% max(cc_periods)) # LACA data
data_ECE1 <- subset(data, p %in% cc_periods) # ECE data
data_NCC1 <- data # NCC data

# Estimand values  
ACA_1 = mean(data_ACA1$y1)-mean(data_ACA1$y0)
LACA_1 = mean(data_LACA1$y1)-mean(data_LACA1$y0)
ECE_1 = mean(data_ECE1$y1)-mean(data_ECE1$y0)
NCC_1 = mean(data_NCC1$y1)-mean(data_NCC1$y0)

# Treatment 2
trt = 2
cc_periods <- which(as.vector(table(data$t, data$p)[trt+1,])>0)

# Analysis data 
data_ACA2 <- subset(data, t %in% c(0,trt) & p %in% cc_periods)  # ACA data #WRONG, no?
data_LACA2 <- subset(data, p %in% max(cc_periods)) # LACA data
data_ECE2 <- subset(data, p %in% cc_periods) # ECE data
data_NCC2 <- data # NCC data

# Estimand values  
ACA_2 = mean(data_ACA2$y2)-mean(data_ACA2$y0)
LACA_2 = mean(data_LACA2$y2)-mean(data_LACA2$y0)
ECE_2 = mean(data_ECE2$y2)-mean(data_ECE2$y0)
NCC_2 = mean(data_NCC2$y2)-mean(data_NCC2$y0)

##############################################
##############################################
##############################################

### "Unequal alloc, no trend",
mu0 = 0; mu1 = 1; mu2 = 2; beta = 1 
N1 = 150000; N2 = 300000; N3 = 150000
alloc1 = c(2, 1)/3        # P1: 100 control, 50 trt1
alloc2 = c(2, 1, 1)/4     # P2: 150 control, 75 trt1, 75 trt2
alloc3 = c(2, 1)/3        # P3: 100 control, 50 trt2
lambda = 0
# Total: Control=350, Trt1=125, Trt2=125

set.seed(12345)
data <- simdata_blocked_unequal(
  mu0 = mu0, mu1 = mu1, mu2 = mu2,
  beta = beta,
  N1 = N1, N2 = N2, N3 = N3,
  alloc1 = alloc1, alloc2 = alloc2, alloc3 = alloc3,
  lambda = lambda
)

# Treatment 1
trt = 1
cc_periods <- which(as.vector(table(data$t, data$p)[trt+1,])>0)

# Analysis data 
data_ACA1 <- subset(data, t %in% c(0,trt) & p %in% cc_periods)  # ACA data #WRONG, no?
data_LACA1 <- subset(data, p %in% max(cc_periods)) # LACA data
data_ECE1 <- subset(data, p %in% cc_periods) # ECE data
data_NCC1 <- data # NCC data

# Estimand values  
ACA_1 = mean(data_ACA1$y1)-mean(data_ACA1$y0)
LACA_1 = mean(data_LACA1$y1)-mean(data_LACA1$y0)
ECE_1 = mean(data_ECE1$y1)-mean(data_ECE1$y0)
NCC_1 = mean(data_NCC1$y1)-mean(data_NCC1$y0)

# Treatment 2
trt = 2
cc_periods <- which(as.vector(table(data$t, data$p)[trt+1,])>0)

# Analysis data 
data_ACA2 <- subset(data, t %in% c(0,trt) & p %in% cc_periods)  # ACA data #WRONG, no?
data_LACA2 <- subset(data, p %in% max(cc_periods)) # LACA data
data_ECE2 <- subset(data, p %in% cc_periods) # ECE data
data_NCC2 <- data # NCC data

# Estimand values  
ACA_2 = mean(data_ACA2$y2)-mean(data_ACA2$y0)
LACA_2 = mean(data_LACA2$y2)-mean(data_LACA2$y0)
ECE_2 = mean(data_ECE2$y2)-mean(data_ECE2$y0)
NCC_2 = mean(data_NCC2$y2)-mean(data_NCC2$y0)


##############################################
##############################################
##############################################

### "Equal alloc, strong trend"
mu0 = 0; mu1 = 1; mu2 = 2; beta = 1 
N1 = 100000; N2 = 100000; N3 = 100000
alloc1 = c(1, 1)/2
alloc2 = c(1, 1, 1)/3
alloc3 = c(1, 1)/2
lambda = 1
# Total: Control=133, Trt1=83, Trt2=83

set.seed(12345)
data <- simdata_blocked_unequal(
  mu0 = mu0, mu1 = mu1, mu2 = mu2,
  beta = beta,
  N1 = N1, N2 = N2, N3 = N3,
  alloc1 = alloc1, alloc2 = alloc2, alloc3 = alloc3,
  lambda = lambda
)

# Treatment 1
trt = 1
cc_periods <- which(as.vector(table(data$t, data$p)[trt+1,])>0)

# Analysis data 
data_ACA1 <- subset(data, t %in% c(0,trt) & p %in% cc_periods)  # ACA data #WRONG, no?
data_LACA1 <- subset(data, p %in% max(cc_periods)) # LACA data
data_ECE1 <- subset(data, p %in% cc_periods) # ECE data
data_NCC1 <- data # NCC data

# Estimand values  
ACA_1 = mean(data_ACA1$y1)-mean(data_ACA1$y0)
LACA_1 = mean(data_LACA1$y1)-mean(data_LACA1$y0)
ECE_1 = mean(data_ECE1$y1)-mean(data_ECE1$y0)
NCC_1 = mean(data_NCC1$y1)-mean(data_NCC1$y0)

# Treatment 2
trt = 2
cc_periods <- which(as.vector(table(data$t, data$p)[trt+1,])>0)

# Analysis data 
data_ACA2 <- subset(data, t %in% c(0,trt) & p %in% cc_periods)  # ACA data #WRONG, no?
data_LACA2 <- subset(data, p %in% max(cc_periods)) # LACA data
data_ECE2 <- subset(data, p %in% cc_periods) # ECE data
data_NCC2 <- data # NCC data

# Estimand values  
ACA_2 = mean(data_ACA2$y2)-mean(data_ACA2$y0)
LACA_2 = mean(data_LACA2$y2)-mean(data_LACA2$y0)
ECE_2 = mean(data_ECE2$y2)-mean(data_ECE2$y0)
NCC_2 = mean(data_NCC2$y2)-mean(data_NCC2$y0)


##############################################
##############################################
##############################################

### "Unequal alloc, trend + large beta"
mu0 = 0; mu1 = 1; mu2 = 2; beta = 2
N1 = 150000; N2 = 300000; N3 = 150000
alloc1 = c(2, 1)/3
alloc2 = c(2, 1, 1)/4
alloc3 = c(2, 1)/3
lambda = 1
# Total: Control=350, Trt1=125, Trt2=125

set.seed(12345)
data <- simdata_blocked_unequal(
  mu0 = mu0, mu1 = mu1, mu2 = mu2,
  beta = beta,
  N1 = N1, N2 = N2, N3 = N3,
  alloc1 = alloc1, alloc2 = alloc2, alloc3 = alloc3,
  lambda = lambda
)

# Treatment 1
trt = 1
cc_periods <- which(as.vector(table(data$t, data$p)[trt+1,])>0)

# Analysis data 
data_ACA1 <- subset(data, t %in% c(0,trt) & p %in% cc_periods)  # ACA data #WRONG, no?
data_LACA1 <- subset(data, p %in% max(cc_periods)) # LACA data
data_ECE1 <- subset(data, p %in% cc_periods) # ECE data
data_NCC1 <- data # NCC data

# Estimand values  
ACA_1 = mean(data_ACA1$y1)-mean(data_ACA1$y0)
LACA_1 = mean(data_LACA1$y1)-mean(data_LACA1$y0)
ECE_1 = mean(data_ECE1$y1)-mean(data_ECE1$y0)
NCC_1 = mean(data_NCC1$y1)-mean(data_NCC1$y0)

# Treatment 2
trt = 2
cc_periods <- which(as.vector(table(data$t, data$p)[trt+1,])>0)

# Analysis data 
data_ACA2 <- subset(data, t %in% c(0,trt) & p %in% cc_periods)  # ACA data #WRONG, no?
data_LACA2 <- subset(data, p %in% max(cc_periods)) # LACA data
data_ECE2 <- subset(data, p %in% cc_periods) # ECE data
data_NCC2 <- data # NCC data

# Estimand values  
ACA_2 = mean(data_ACA2$y2)-mean(data_ACA2$y0)
LACA_2 = mean(data_LACA2$y2)-mean(data_LACA2$y0)
ECE_2 = mean(data_ECE2$y2)-mean(data_ECE2$y0)
NCC_2 = mean(data_NCC2$y2)-mean(data_NCC2$y0)


##############################################
##############################################
##############################################

### "Equal alloc, trend, treat-cov interaction"
mu0 = 0; mu1 = 1; mu2 = 2; beta = 2
gamma1 = 2; gamma2 = 0;
N1 = 100000; N2 = 100000; N3 = 100000
alloc1 = c(1, 1)/2
alloc2 = c(1, 1, 1)/3
alloc3 = c(1, 1)/2
lambda = 1 
# Total: Control=133, Trt1=83, Trt2=83


set.seed(12345)
data <- simdata_blocked_unequal(
  mu0 = mu0, mu1 = mu1, mu2 = mu2,
  beta = beta,
  gamma1 = gamma1, gamma2=gamma2,
  N1 = N1, N2 = N2, N3 = N3,
  alloc1 = alloc1, alloc2 = alloc2, alloc3 = alloc3,
  lambda = lambda
)

# Treatment 1
trt = 1
cc_periods <- which(as.vector(table(data$t, data$p)[trt+1,])>0)

# Analysis data 
data_ACA1 <- subset(data, t %in% c(0,trt) & p %in% cc_periods)  # ACA data #WRONG, no?
data_LACA1 <- subset(data, p %in% max(cc_periods)) # LACA data
data_ECE1 <- subset(data, p %in% cc_periods) # ECE data
data_NCC1 <- data # NCC data

# Estimand values  
ACA_1 = mean(data_ACA1$y1)-mean(data_ACA1$y0)
LACA_1 = mean(data_LACA1$y1)-mean(data_LACA1$y0)
ECE_1 = mean(data_ECE1$y1)-mean(data_ECE1$y0)
NCC_1 = mean(data_NCC1$y1)-mean(data_NCC1$y0)

# Treatment 2
trt = 2
cc_periods <- which(as.vector(table(data$t, data$p)[trt+1,])>0)

# Analysis data 
data_ACA2 <- subset(data, t %in% c(0,trt) & p %in% cc_periods)  # ACA data #WRONG, no?
data_LACA2 <- subset(data, p %in% max(cc_periods)) # LACA data
data_ECE2 <- subset(data, p %in% cc_periods) # ECE data
data_NCC2 <- data # NCC data

# Estimand values  
ACA_2 = mean(data_ACA2$y2)-mean(data_ACA2$y0)
LACA_2 = mean(data_LACA2$y2)-mean(data_LACA2$y0)
ECE_2 = mean(data_ECE2$y2)-mean(data_ECE2$y0)
NCC_2 = mean(data_NCC2$y2)-mean(data_NCC2$y0)


##############################################
##############################################
##############################################

### "Unequal alloc, trend, treat-cov interaction", 
mu0 = 0; mu1 = 1; mu2 = 2; beta = 2 
gamma1 = 2; gamma2 = 1
N1 = 150000; N2 = 300000; N3 = 150000
alloc1 = c(2, 1)/3
alloc2 = c(2, 1, 1)/4
alloc3 = c(2, 1)/3
lambda = 1
# Total: Control=350, Trt1=125, Trt2=125


set.seed(12345)
data <- simdata_blocked_unequal(
  mu0 = mu0, mu1 = mu1, mu2 = mu2,
  beta = beta,
  gamma1 = gamma1, gamma2=gamma2,
  N1 = N1, N2 = N2, N3 = N3,
  alloc1 = alloc1, alloc2 = alloc2, alloc3 = alloc3,
  lambda = lambda
)

# Treatment 1
trt = 1
cc_periods <- which(as.vector(table(data$t, data$p)[trt+1,])>0)

# Analysis data 
data_ACA1 <- subset(data, t %in% c(0,trt) & p %in% cc_periods)  # ACA data #WRONG, no?
data_LACA1 <- subset(data, p %in% max(cc_periods)) # LACA data
data_ECE1 <- subset(data, p %in% cc_periods) # ECE data
data_NCC1 <- data # NCC data

# Estimand values  
ACA_1 = mean(data_ACA1$y1)-mean(data_ACA1$y0)
LACA_1 = mean(data_LACA1$y1)-mean(data_LACA1$y0)
ECE_1 = mean(data_ECE1$y1)-mean(data_ECE1$y0)
NCC_1 = mean(data_NCC1$y1)-mean(data_NCC1$y0)

# Treatment 2
trt = 2
cc_periods <- which(as.vector(table(data$t, data$p)[trt+1,])>0)

# Analysis data 
data_ACA2 <- subset(data, t %in% c(0,trt) & p %in% cc_periods)  # ACA data #WRONG, no?
data_LACA2 <- subset(data, p %in% max(cc_periods)) # LACA data
data_ECE2 <- subset(data, p %in% cc_periods) # ECE data
data_NCC2 <- data # NCC data

# Estimand values  
ACA_2 = mean(data_ACA2$y2)-mean(data_ACA2$y0)
LACA_2 = mean(data_LACA2$y2)-mean(data_LACA2$y0)
ECE_2 = mean(data_ECE2$y2)-mean(data_ECE2$y0)
NCC_2 = mean(data_NCC2$y2)-mean(data_NCC2$y0)


##############################################
##############################################
##############################################

### "Equal alloc, trend, period-cov interaction", 
mu0 = 0; mu1 = 1; mu2 = 2; beta = 2; delta = 2
N1 = 100000; N2 = 100000; N3 = 100000
alloc1 = c(1, 1)/2
alloc2 = c(1, 1, 1)/3
alloc3 = c(1, 1)/2
lambda = 1
# Total: Control=133, Trt1=83, Trt2=83


set.seed(12345)
data <- simdata_blocked_unequal(
  mu0 = mu0, mu1 = mu1, mu2 = mu2,
  delta=delta,
  beta = beta,
  N1 = N1, N2 = N2, N3 = N3,
  alloc1 = alloc1, alloc2 = alloc2, alloc3 = alloc3,
  lambda = lambda
)

# Treatment 1
trt = 1
cc_periods <- which(as.vector(table(data$t, data$p)[trt+1,])>0)

# Analysis data 
data_ACA1 <- subset(data, t %in% c(0,trt) & p %in% cc_periods)  # ACA data #WRONG, no?
data_LACA1 <- subset(data, p %in% max(cc_periods)) # LACA data
data_ECE1 <- subset(data, p %in% cc_periods) # ECE data
data_NCC1 <- data # NCC data

# Estimand values  
ACA_1 = mean(data_ACA1$y1)-mean(data_ACA1$y0)
LACA_1 = mean(data_LACA1$y1)-mean(data_LACA1$y0)
ECE_1 = mean(data_ECE1$y1)-mean(data_ECE1$y0)
NCC_1 = mean(data_NCC1$y1)-mean(data_NCC1$y0)

# Treatment 2
trt = 2
cc_periods <- which(as.vector(table(data$t, data$p)[trt+1,])>0)

# Analysis data 
data_ACA2 <- subset(data, t %in% c(0,trt) & p %in% cc_periods)  # ACA data #WRONG, no?
data_LACA2 <- subset(data, p %in% max(cc_periods)) # LACA data
data_ECE2 <- subset(data, p %in% cc_periods) # ECE data
data_NCC2 <- data # NCC data

# Estimand values  
ACA_2 = mean(data_ACA2$y2)-mean(data_ACA2$y0)
LACA_2 = mean(data_LACA2$y2)-mean(data_LACA2$y0)
ECE_2 = mean(data_ECE2$y2)-mean(data_ECE2$y0)
NCC_2 = mean(data_NCC2$y2)-mean(data_NCC2$y0)


##############################################
##############################################
##############################################

### "Unequal alloc, trend, period-cov interaction" 
mu0 = 0; mu1 = 1; mu2 = 2; beta = 2; delta = 2
N1 = 150000; N2 = 300000; N3 = 150000
alloc1 = c(2, 1)/3
alloc2 = c(2, 1, 1)/4
alloc3 = c(2, 1)/3
lambda = 1
# Total: Control=350, Trt1=125, Trt2=125


set.seed(12345)
data <- simdata_blocked_unequal(
  mu0 = mu0, mu1 = mu1, mu2 = mu2,
  delta=delta,
  beta = beta,
  N1 = N1, N2 = N2, N3 = N3,
  alloc1 = alloc1, alloc2 = alloc2, alloc3 = alloc3,
  lambda = lambda
)

# Treatment 1
trt = 1
cc_periods <- which(as.vector(table(data$t, data$p)[trt+1,])>0)

# Analysis data 
data_ACA1 <- subset(data, t %in% c(0,trt) & p %in% cc_periods)  # ACA data #WRONG, no?
data_LACA1 <- subset(data, p %in% max(cc_periods)) # LACA data
data_ECE1 <- subset(data, p %in% cc_periods) # ECE data
data_NCC1 <- data # NCC data

# Estimand values  
ACA_1 = mean(data_ACA1$y1)-mean(data_ACA1$y0)
LACA_1 = mean(data_LACA1$y1)-mean(data_LACA1$y0)
ECE_1 = mean(data_ECE1$y1)-mean(data_ECE1$y0)
NCC_1 = mean(data_NCC1$y1)-mean(data_NCC1$y0)

# Treatment 2
trt = 2
cc_periods <- which(as.vector(table(data$t, data$p)[trt+1,])>0)

# Analysis data 
data_ACA2 <- subset(data, t %in% c(0,trt) & p %in% cc_periods)  # ACA data #WRONG, no?
data_LACA2 <- subset(data, p %in% max(cc_periods)) # LACA data
data_ECE2 <- subset(data, p %in% cc_periods) # ECE data
data_NCC2 <- data # NCC data

# Estimand values  
ACA_2 = mean(data_ACA2$y2)-mean(data_ACA2$y0)
LACA_2 = mean(data_LACA2$y2)-mean(data_LACA2$y0)
ECE_2 = mean(data_ECE2$y2)-mean(data_ECE2$y0)
NCC_2 = mean(data_NCC2$y2)-mean(data_NCC2$y0)


##############################################
##############################################
##############################################

### "Equal alloc, treat-period interaction (waning trt1)",
mu0 = 0; mu1 = 1; mu2 = 2; beta = 1 
tau1 = -1; tau2 = 0  # Trt1 loses 1 unit per period
N1 = 100000; N2 = 100000; N3 = 100000
alloc1 = c(1, 1)/2
alloc2 = c(1, 1, 1)/3
alloc3 = c(1, 1)/2
lambda = 1
# Total: Control=133, Trt1=83, Trt2=83
# Effect trajectory: Trt1 goes 3→2→1, Trt2 stays at 2


set.seed(12345)
data <- simdata_blocked_unequal(
  mu0 = mu0, mu1 = mu1, mu2 = mu2,
  tau1=tau1, tau2=tau2,
  beta = beta,
  N1 = N1, N2 = N2, N3 = N3,
  alloc1 = alloc1, alloc2 = alloc2, alloc3 = alloc3,
  lambda = lambda
)

# Treatment 1
trt = 1
cc_periods <- which(as.vector(table(data$t, data$p)[trt+1,])>0)

# Analysis data 
data_ACA1 <- subset(data, t %in% c(0,trt) & p %in% cc_periods)  # ACA data #WRONG, no?
data_LACA1 <- subset(data, p %in% max(cc_periods)) # LACA data
data_ECE1 <- subset(data, p %in% cc_periods) # ECE data
data_NCC1 <- data # NCC data

# Estimand values  
ACA_1 = mean(data_ACA1$y1)-mean(data_ACA1$y0)
LACA_1 = mean(data_LACA1$y1)-mean(data_LACA1$y0)
ECE_1 = mean(data_ECE1$y1)-mean(data_ECE1$y0)
NCC_1 = mean(data_NCC1$y1)-mean(data_NCC1$y0)

# Treatment 2
trt = 2
cc_periods <- which(as.vector(table(data$t, data$p)[trt+1,])>0)

# Analysis data 
data_ACA2 <- subset(data, t %in% c(0,trt) & p %in% cc_periods)  # ACA data #WRONG, no?
data_LACA2 <- subset(data, p %in% max(cc_periods)) # LACA data
data_ECE2 <- subset(data, p %in% cc_periods) # ECE data
data_NCC2 <- data # NCC data

# Estimand values  
ACA_2 = mean(data_ACA2$y2)-mean(data_ACA2$y0)
LACA_2 = mean(data_LACA2$y2)-mean(data_LACA2$y0)
ECE_2 = mean(data_ECE2$y2)-mean(data_ECE2$y0)
NCC_2 = mean(data_NCC2$y2)-mean(data_NCC2$y0)



##############################################
##############################################
##############################################

### "Unequal alloc, treat-period interaction (delayed trt2)",
mu0 = 0; mu1 = 1; mu2 = 2; beta = 1
tau1 = 0; tau2 = 1  # Trt2 gains 1 unit per period
N1 = 150000; N2 = 300000; N3 = 150000
alloc1 = c(2, 1)/3
alloc2 = c(2, 1, 1)/4
alloc3 = c(2, 1)/3
lambda = 1
# Total: Control=350, Trt1=125, Trt2=125
# Effect trajectory: Trt1 stays at 2, Trt2 goes 0→1→2


set.seed(12345)
data <- simdata_blocked_unequal(
  mu0 = mu0, mu1 = mu1, mu2 = mu2,
  tau1=tau1, tau2=tau2,
  beta = beta,
  N1 = N1, N2 = N2, N3 = N3,
  alloc1 = alloc1, alloc2 = alloc2, alloc3 = alloc3,
  lambda = lambda
)

# Treatment 1
trt = 1
cc_periods <- which(as.vector(table(data$t, data$p)[trt+1,])>0)

# Analysis data 
data_ACA1 <- subset(data, t %in% c(0,trt) & p %in% cc_periods)  # ACA data #WRONG, no?
data_LACA1 <- subset(data, p %in% max(cc_periods)) # LACA data
data_ECE1 <- subset(data, p %in% cc_periods) # ECE data
data_NCC1 <- data # NCC data

# Estimand values  
ACA_1 = mean(data_ACA1$y1)-mean(data_ACA1$y0)
LACA_1 = mean(data_LACA1$y1)-mean(data_LACA1$y0)
ECE_1 = mean(data_ECE1$y1)-mean(data_ECE1$y0)
NCC_1 = mean(data_NCC1$y1)-mean(data_NCC1$y0)

# Treatment 2
trt = 2
cc_periods <- which(as.vector(table(data$t, data$p)[trt+1,])>0)

# Analysis data 
data_ACA2 <- subset(data, t %in% c(0,trt) & p %in% cc_periods)  # ACA data #WRONG, no?
data_LACA2 <- subset(data, p %in% max(cc_periods)) # LACA data
data_ECE2 <- subset(data, p %in% cc_periods) # ECE data
data_NCC2 <- data # NCC data

# Estimand values  
ACA_2 = mean(data_ACA2$y2)-mean(data_ACA2$y0)
LACA_2 = mean(data_LACA2$y2)-mean(data_LACA2$y0)
ECE_2 = mean(data_ECE2$y2)-mean(data_ECE2$y0)
NCC_2 = mean(data_NCC2$y2)-mean(data_NCC2$y0)



##############################################
##############################################
##############################################

### "Equal alloc, opposing treat-period effects",
mu0 = 0; mu1 = 1; mu2 = 2; beta = 1
tau1 = -1; tau2 = 1  # Trt1 wanes, Trt2 improves (strong effects)
N1 = 100000; N2 = 100000; N3 = 100000
alloc1 = c(1, 1)/2
alloc2 = c(1, 1, 1)/3
alloc3 = c(1, 1)/2
lambda = 1
# Total: Control=133, Trt1=83, Trt2=83
# Effect trajectory: Trt1: 3→2→1, Trt2: 0→1→2 (crossing effects!)


set.seed(12345)
data <- simdata_blocked_unequal(
  mu0 = mu0, mu1 = mu1, mu2 = mu2,
  tau1=tau1, tau2=tau2,
  beta = beta,
  N1 = N1, N2 = N2, N3 = N3,
  alloc1 = alloc1, alloc2 = alloc2, alloc3 = alloc3,
  lambda = lambda
)

# Treatment 1
trt = 1
cc_periods <- which(as.vector(table(data$t, data$p)[trt+1,])>0)

# Analysis data 
data_ACA1 <- subset(data, t %in% c(0,trt) & p %in% cc_periods)  # ACA data #WRONG, no?
data_LACA1 <- subset(data, p %in% max(cc_periods)) # LACA data
data_ECE1 <- subset(data, p %in% cc_periods) # ECE data
data_NCC1 <- data # NCC data

# Estimand values  
ACA_1 = mean(data_ACA1$y1)-mean(data_ACA1$y0)
LACA_1 = mean(data_LACA1$y1)-mean(data_LACA1$y0)
ECE_1 = mean(data_ECE1$y1)-mean(data_ECE1$y0)
NCC_1 = mean(data_NCC1$y1)-mean(data_NCC1$y0)

# Treatment 2
trt = 2
cc_periods <- which(as.vector(table(data$t, data$p)[trt+1,])>0)

# Analysis data 
data_ACA2 <- subset(data, t %in% c(0,trt) & p %in% cc_periods)  # ACA data #WRONG, no?
data_LACA2 <- subset(data, p %in% max(cc_periods)) # LACA data
data_ECE2 <- subset(data, p %in% cc_periods) # ECE data
data_NCC2 <- data # NCC data

# Estimand values  
ACA_2 = mean(data_ACA2$y2)-mean(data_ACA2$y0)
LACA_2 = mean(data_LACA2$y2)-mean(data_LACA2$y0)
ECE_2 = mean(data_ECE2$y2)-mean(data_ECE2$y0)
NCC_2 = mean(data_NCC2$y2)-mean(data_NCC2$y0)


##############################################
##############################################
##############################################

### "Unequal alloc, extreme all interactions",
mu0 = 0; mu1 = 1; mu2 = 2; beta = 2
gamma1 = 2; gamma2 = -2     # Strong opposing covariate effects
tau1 = 1; tau2 = -1          # Opposing period trends
delta = 2                   # Strong period-covariate interaction
N1 = 150000; N2 = 300000; N3 = 150000
alloc1 = c(2, 1)/3
alloc2 = c(2, 1, 1)/4
alloc3 = c(2, 1)/3
lambda = 2  # Even stronger time trend
# Total: Control=350, Trt1=125, Trt2=125 

set.seed(12345)
data <- simdata_blocked_unequal(
  mu0 = mu0, mu1 = mu1, mu2 = mu2,
  tau1=tau1, tau2=tau2,
  beta = beta,
  delta=delta,
  gamma1=gamma1, gamma2=gamma2,
  N1 = N1, N2 = N2, N3 = N3,
  alloc1 = alloc1, alloc2 = alloc2, alloc3 = alloc3,
  lambda = lambda
)

# Treatment 1
trt = 1
cc_periods <- which(as.vector(table(data$t, data$p)[trt+1,])>0)

# Analysis data 
data_ACA1 <- subset(data, t %in% c(0,trt) & p %in% cc_periods)  # ACA data #WRONG, no?
data_LACA1 <- subset(data, p %in% max(cc_periods)) # LACA data
data_ECE1 <- subset(data, p %in% cc_periods) # ECE data
data_NCC1 <- data # NCC data

# Estimand values  
ACA_1 = mean(data_ACA1$y1)-mean(data_ACA1$y0)
LACA_1 = mean(data_LACA1$y1)-mean(data_LACA1$y0)
ECE_1 = mean(data_ECE1$y1)-mean(data_ECE1$y0)
NCC_1 = mean(data_NCC1$y1)-mean(data_NCC1$y0)

# Treatment 2
trt = 2
cc_periods <- which(as.vector(table(data$t, data$p)[trt+1,])>0)

# Analysis data 
data_ACA2 <- subset(data, t %in% c(0,trt) & p %in% cc_periods)  # ACA data #WRONG, no?
data_LACA2 <- subset(data, p %in% max(cc_periods)) # LACA data
data_ECE2 <- subset(data, p %in% cc_periods) # ECE data
data_NCC2 <- data # NCC data

# Estimand values  
ACA_2 = mean(data_ACA2$y2)-mean(data_ACA2$y0)
LACA_2 = mean(data_LACA2$y2)-mean(data_LACA2$y0)
ECE_2 = mean(data_ECE2$y2)-mean(data_ECE2$y0)
NCC_2 = mean(data_NCC2$y2)-mean(data_NCC2$y0)


##############################################
##############################################
##############################################

### "Equal alloc, extreme covariate interactions only"
mu0 = 0; mu1 = 1; mu2 = 2; beta = 3
gamma1 = 3; gamma2 = -3  # Very strong opposing interactions
N1 = 100000; N2 = 100000; N3 = 100000
alloc1 = c(1, 1)/2
alloc2 = c(1, 1, 1)/3
alloc3 = c(1, 1)/2
lambda = 0
# Total: Control=133, Trt1=83, Trt2=83

set.seed(12345)
data <- simdata_blocked_unequal(
  mu0 = mu0, mu1 = mu1, mu2 = mu2,
  beta = beta,
  N1 = N1, N2 = N2, N3 = N3,
  alloc1 = alloc1, alloc2 = alloc2, alloc3 = alloc3,
  lambda = lambda
)

# Treatment 1
trt = 1
cc_periods <- which(as.vector(table(data$t, data$p)[trt+1,])>0)

# Analysis data 
data_ACA1 <- subset(data, t %in% c(0,trt) & p %in% cc_periods)  # ACA data #WRONG, no?
data_LACA1 <- subset(data, p %in% max(cc_periods)) # LACA data
data_ECE1 <- subset(data, p %in% cc_periods) # ECE data
data_NCC1 <- data # NCC data

# Estimand values  
ACA_1 = mean(data_ACA1$y1)-mean(data_ACA1$y0)
LACA_1 = mean(data_LACA1$y1)-mean(data_LACA1$y0)
ECE_1 = mean(data_ECE1$y1)-mean(data_ECE1$y0)
NCC_1 = mean(data_NCC1$y1)-mean(data_NCC1$y0)

# Treatment 2
trt = 2
cc_periods <- which(as.vector(table(data$t, data$p)[trt+1,])>0)

# Analysis data 
data_ACA2 <- subset(data, t %in% c(0,trt) & p %in% cc_periods)  # ACA data #WRONG, no?
data_LACA2 <- subset(data, p %in% max(cc_periods)) # LACA data
data_ECE2 <- subset(data, p %in% cc_periods) # ECE data
data_NCC2 <- data # NCC data

# Estimand values  
ACA_2 = mean(data_ACA2$y2)-mean(data_ACA2$y0)
LACA_2 = mean(data_LACA2$y2)-mean(data_LACA2$y0)
ECE_2 = mean(data_ECE2$y2)-mean(data_ECE2$y0)
NCC_2 = mean(data_NCC2$y2)-mean(data_NCC2$y0)


##############################################
##############################################
##############################################

### "Unequal alloc, extreme period effects"
mu0 = 0; mu1 = 1; mu2 = 2; beta = 1
tau1 = 2; tau2 = 2     # Both treatments improve dramatically
delta = 3              # Period-covariate also very strong
N1 = 150000; N2 = 300000; N3 = 150000
alloc1 = c(2, 1)/3
alloc2 = c(2, 1, 1)/4
alloc3 = c(2, 1)/3
lambda = 2
# Total: Control=350, Trt1=125, Trt2=125

set.seed(12345)
data <- simdata_blocked_unequal(
  mu0 = mu0, mu1 = mu1, mu2 = mu2,
  beta = beta, tau1 = tau1, tau2 = tau2,
  N1 = N1, N2 = N2, N3 = N3,
  alloc1 = alloc1, alloc2 = alloc2, alloc3 = alloc3,
  lambda = lambda
)

# Treatment 1
trt = 1
cc_periods <- which(as.vector(table(data$t, data$p)[trt+1,])>0)

# Analysis data 
data_ACA1 <- subset(data, t %in% c(0,trt) & p %in% cc_periods)  # ACA data #WRONG, no?
data_LACA1 <- subset(data, p %in% max(cc_periods)) # LACA data
data_ECE1 <- subset(data, p %in% cc_periods) # ECE data
data_NCC1 <- data # NCC data

# Estimand values  
ACA_1 = mean(data_ACA1$y1)-mean(data_ACA1$y0)
LACA_1 = mean(data_LACA1$y1)-mean(data_LACA1$y0)
ECE_1 = mean(data_ECE1$y1)-mean(data_ECE1$y0)
NCC_1 = mean(data_NCC1$y1)-mean(data_NCC1$y0)

# Treatment 2
trt = 2
cc_periods <- which(as.vector(table(data$t, data$p)[trt+1,])>0)

# Analysis data 
data_ACA2 <- subset(data, t %in% c(0,trt) & p %in% cc_periods)  # ACA data #WRONG, no?
data_LACA2 <- subset(data, p %in% max(cc_periods)) # LACA data
data_ECE2 <- subset(data, p %in% cc_periods) # ECE data
data_NCC2 <- data # NCC data

# Estimand values  
ACA_2 = mean(data_ACA2$y2)-mean(data_ACA2$y0)
LACA_2 = mean(data_LACA2$y2)-mean(data_LACA2$y0)
ECE_2 = mean(data_ECE2$y2)-mean(data_ECE2$y0)
NCC_2 = mean(data_NCC2$y2)-mean(data_NCC2$y0)

