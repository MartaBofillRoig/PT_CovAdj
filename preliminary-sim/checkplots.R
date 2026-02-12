n = 1000;
sc <- list(name = "Equal alloc, treat-period interaction (waning trt1)",
     mu0 = 0, mu1 = 1, mu2 = 2, beta = 0, 
     tau1 = 3, tau2 = 1,   
     N1 = n* 100, N2 = n* 100, N3 = n* 100,
     alloc1 = c(1, 1)/2,
     alloc2 = c(1, 1, 1)/3,
     alloc3 = c(1, 1)/2,
     lambda = 1,
     gamma1 = 0, gamma2 = 0,
     delta = 0, trendp = "stepwise", sd=1)

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



# Plot means per period for different treatments
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


lmmodel(dat, trt = 2, dataset = "ACA")
lmmodel(dat, trt = 2, dataset = "ECE")
lmmodel(dat, trt = 2, dataset = "NCC")

# Create observation index
dat <- dat %>%
  arrange(p, t) %>%
  mutate(obs_index = row_number())



# Plot 1: Individual observations with means overlaid
p1 <- ggplot() +
  geom_point(data = dat, 
             aes(x = obs_index, y = y, color = as.factor(t)), 
             alpha = 0.3, size = 1) +
  geom_vline(data = dat %>% 
               group_by(p) %>% 
               summarise(period_start = min(obs_index), .groups = 'drop'),
             aes(xintercept = period_start), 
             linetype = "dashed", color = "gray50", alpha = 0.5) +
  scale_color_brewer(palette = "Set1", name = "Treatment") +
  labs(
    title = "Response by Observation Index and Treatment",
    x = "Observation Index",
    y = "Response (y)"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Add mean segments manually
period_means <- dat %>% 
  group_by(p, t) %>% 
  summarise(start_idx = min(obs_index),
            end_idx = max(obs_index),
            mean_y = mean(y),
            .groups = 'drop')

for(i in 1:nrow(period_means)) {
  p1 <- p1 + geom_segment(
    x = period_means$start_idx[i], 
    xend = period_means$end_idx[i],
    y = period_means$mean_y[i], 
    yend = period_means$mean_y[i],
    color = RColorBrewer::brewer.pal(max(3, length(unique(dat$t))), "Set1")[period_means$t[i] + 1],
    size = 1.2
  )
}

# Plot 2: Mean response by period
p2 <- ggplot(means_summary, aes(x = as.factor(p), y = mean_y, 
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

# Plot 3: Faceted by period
p3 <- ggplot(dat, aes(x = obs_index, y = y, color = as.factor(t))) +
  geom_point(alpha = 0.4, size = 1) +
  stat_smooth(method = "lm", se = FALSE, size = 1) +
  facet_wrap(~ p, scales = "free_x", labeller = label_both) +
  scale_color_brewer(palette = "Set1", name = "Treatment") +
  labs(
    title = "Response by Observation Index (Faceted by Period)",
    x = "Observation Index",
    y = "Response (y)"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Display plots
print(p1)
print(p2)
print(p3)