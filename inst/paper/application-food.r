# ------------------------------------------------------------
# packages 
# ------------------------------------------------------------
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("tidyr")
# install.packages("CSB")   
# install.packages(betareg)

library(ggplot2)
library(dplyr)
library(tidyr)
library(CSB)
library(betareg)

# ------------------------------------------------------------
# Food data
# ------------------------------------------------------------

data(FoodExpenditure)

foo <- FoodExpenditure$food

inc <- FoodExpenditure$income

x <- foo/inc

# ------------------------------------------------------------
# Fit competing models
# ------------------------------------------------------------

cmp <- compareCSBmodels(x)

summary(cmp)

# Comparison of fitted models for unit data
# -----------------------------------------
#
# Call:
# compareCSBmodels(x = x)
#
# Sample size: 38 
#
# Model comparison:
#  Model npar  logLik      AIC      BIC convergence
#    CSB    2 35.8413 -67.6826 -64.4074           0
#      B    2 35.3464 -66.6929 -63.4177           0
#      K    2 33.4891 -62.9782 -59.7030           0
#
# Parameter estimates:
#  Model parameter estimate      se
#      B    shape1   6.0716  1.3586
#      B    shape2  14.8221  3.3988
#      K         a   2.9546  0.3692
#      K         b  26.9654 10.8271
#    CSB     shape   6.4869  1.8570
#    CSB         q   5.8294  1.9281
#
# Best model according to AIC: CSB 
# Best model according to BIC: CSB

gofit <- gofCSB_boot(x = x, 
                     B = 1000, 
                     method = "L-BFGS-B", 
                     multistart = FALSE, 
                     seed = 2025)

summary(gofit)

# Bootstrap goodness-of-fit for the CSB distribution
# --------------------------------------------------
#
# Call:
# gofCSB_boot(x = x, B = 1000, seed = 2025, method = "L-BFGS-B", 
#     multistart = FALSE)
# 
# Sample size: 38 
# Bootstrap replicates: 1000 
#
# Parameter estimates:
#  Parameter Estimate
#      shape   6.4869
#          q   5.8294
#
# Fit statistics:
#   logLik      AIC      BIC
#  35.8413 -67.6826 -64.4074
#
# Goodness-of-fit statistics:
#              Test Statistic p_value
#  Anderson-Darling    0.6145  0.2288
#  Cramer-von Mises    0.1106  0.2128 

1 - pchisq((5.8294/1.9281)^2,1)

# 0.002499591 
# The null hypothesis H_0: q = 1 is rejected at conventional 
# significance levels

# ------------------------------------------------------------
# Fitted histogram
# ------------------------------------------------------------

dens_grid <- seq(0.01, 0.75, length.out = 600)

dens_df <- data.frame(
  x = dens_grid,
  B = dbeta(dens_grid, shape1 = 6.0716, shape2 = 14.8221),
  K = CSB:::dKum(dens_grid, a = 2.9546, b = 26.9654),
  CSB = dCSB(dens_grid, shape = 6.4869, q = 5.8294)
) %>%
  pivot_longer(
    cols = c(B, K, CSB),
    names_to = "Distribution",
    values_to = "Density"
  ) %>%
  mutate(
    Distribution = factor(
      Distribution,
      levels = c("CSB", "B", "K")
    )
  )

ggplot(data.frame(x = x), aes(x = x)) +
  geom_histogram(
    aes(y = after_stat(density)),
    breaks = seq(0, 0.7, by = 0.05),
    fill = "white",
    color = "black",
    alpha = 0.65
  ) +
  geom_line(
    data = dens_df,
    aes(x = x, y = Density,
        linetype = Distribution,
        linewidth = Distribution),
    color = "black"
  ) +
  scale_linetype_manual(
    values = c(
      CSB = "solid",
      B = "dashed",
      K = "dotdash"
    )
  ) +
  scale_linewidth_manual(
    values = c(
      CSB = 1.1,
      B = 0.9,
      K = 0.9
    )
  ) +
  labs(
    x = "Food expenditure share",
    y = "Density function",
    linetype = "Distribution",
    linewidth = "Distribution"
  ) +
  guides(
    linewidth = "none",
    linetype = guide_legend(
      override.aes = list(linewidth = c(1.1, 0.9, 0.9))
    )
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.major = element_line(color = "grey85", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    
    axis.title = element_text(size = 15, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    
    legend.position = c(0.74, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_rect(
      fill = scales::alpha("white", 0.85),
      color = "grey70",
      linewidth = 0.3
    ),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.key.width = unit(1.6, "cm"),
    legend.key.height = unit(0.45, "cm")
  )

# ------------------------------------------------------------
# cumulatives distributions
# ------------------------------------------------------------

shape_csb <- 6.4869
q_csb     <- 5.8294

grid_x <- seq(0.01, 0.75, length.out = 600)

csb_df <- data.frame(
  x = grid_x,
  CDF = pCSB(grid_x, shape = shape_csb, q = q_csb)
)

ggplot() +
  stat_ecdf(
    data = data.frame(x = x),
    aes(x = x, linetype = "Empirical"),
    geom = "step",
    linewidth = 1,
    color = "black"
  ) +
  geom_line(
    data = csb_df,
    aes(x = x, y = CDF, linetype = "CSB"),
    linewidth = 1,
    color = "black"
  ) +
  scale_linetype_manual(
    values = c(
      Empirical = "solid",
      CSB = "longdash"
    )
  ) +
  labs(
    x = "Food expenditure share",
    y = "Cumulative function",
    linetype = "Distribution"
  ) +
  coord_cartesian(xlim = c(0.01, 0.75), ylim = c(0, 1)) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.major = element_line(color = "grey85", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    axis.title = element_text(size = 15, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_rect(
      fill = scales::alpha("white", 0.85),
      color = "grey70",
      linewidth = 0.3
    ),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11)
  )


