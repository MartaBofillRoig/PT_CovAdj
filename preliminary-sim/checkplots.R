
n = 1000;
sc <- list(name = "Equal alloc, treat-period interaction (h1)",
           mu0 = 0, mu1 = 1, mu2 = 2, beta = 0, 
           tau1 = 3, tau2 = 0,   
           N1 = n* 100, N2 = n* 100, N3 = n* 100,
           alloc1 = c(1, 1)/2,
           alloc2 = c(1, 1, 1)/3,
           alloc3 = c(1, 1)/2,
           lambda = 1,
           gamma1 = 0, gamma2 = 0,
           delta = 0, subgroup = F,         # period-covariate interaction
           trendp = "stepwise", sd=1) 

# n=1000;
# mu0 = 0; mu1 = 1; mu2 = 2; beta = 0;
# tau1 = 1; tau2 = 3;
# N1 = n* 100; N2 = n* 100; N3 = n* 100;
# alloc1 = c(1, 1)/2;
# alloc2 = c(1, 1, 1)/3;
# alloc3 = c(1, 1)/2;
# lambda = 1
# gamma1 = 0; gamma2 = 0;
# delta = 0; trendp = "stepwise";sd=1


dat <- simdata_blocked_unequal(
  mu0 = sc$mu0, mu1 = sc$mu1, mu2 = sc$mu2,
  tau1 = sc$tau1, tau2 = sc$tau2,   
  beta = sc$beta,
  N1 = sc$N1, N2 = sc$N2, N3 = sc$N3,
  alloc1 = sc$alloc1, alloc2 = sc$alloc2, alloc3 = sc$alloc3,
  lambda = sc$lambda,
  gamma1 = 0, gamma2 = 0,
  delta = 0, trendp = "stepwise", sd=1
)

# 
library(ggplot2)
library(dplyr)

# Calculate means per period and treatment
means_summary <- dat %>%
  group_by(p, t) %>%
  summarise(
    mean_y = mean(y),
    se_y = sd(y) / sqrt(n()),
    n = n(),
    .groups = 'drop'
  )
means_summary

# 2 vs 0 -> ACA
N_ACA= sc$N2*(2/3)+sc$N3
(means_summary[7,]$mean_y*(sc$N3/N_ACA) + means_summary[5,]$mean_y*(sc$N2*(2/3)/N_ACA))-(means_summary[6,]$mean_y*(sc$N3/N_ACA) + means_summary[3,]$mean_y*(sc$N2*(2/3)/N_ACA))

# 1 vs 0 -> ACA
N_ACA= sc$N2*(2/3)+sc$N1
(means_summary[4,]$mean_y*(sc$N2*(2/3)/N_ACA) + means_summary[2,]$mean_y*(sc$N1/N_ACA))-(means_summary[3,]$mean_y*(sc$N2*(2/3)/N_ACA) + means_summary[1,]$mean_y*(sc$N1/N_ACA))

# 1 vs 0 -> ECE
N_ECE= sc$N2+sc$N1
(means_summary[4,]$mean_y*(sc$N2/N_ECE) + means_summary[2,]$mean_y*(sc$N1/N_ECE))-(means_summary[3,]$mean_y*(sc$N2/N_ECE) + means_summary[1,]$mean_y*(sc$N1/N_ECE))

# 
lmmodel(dat, trt = 1, dataset = "ACA")
lmmodel(dat, trt = 1, dataset = "ECE")
lmmodel(dat, trt = 1, dataset = "NCC")

# 
lmmodel(dat, trt = 2, dataset = "ACA")
lmmodel(dat, trt = 2, dataset = "ECE")
lmmodel(dat, trt = 2, dataset = "NCC")

############################################

# Create observation index
dat <- dat %>%
  arrange(p, t) %>%
  mutate(obs_index = row_number())

# Plot Mean response by period
plot1 <- ggplot(means_summary, aes(x = as.factor(p), y = mean_y, 
                                color = as.factor(t), group = t)) +
  geom_line(linewidth = 1) +
  geom_point(linewidth = 3) +
  geom_errorbar(aes(ymin = mean_y - se_y, ymax = mean_y + se_y), 
                width = 0.2, linewidth = 0.8) +
  scale_color_brewer(palette = "Set1", name = "Treatment") +
  labs(
    title = "Mean Response by Period and Treatment",
    subtitle = "Error bars show ± 1 SE",
    x = "Period",
    y = "Mean Response (y)"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
