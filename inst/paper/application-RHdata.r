# ------------------------------------------------------------
# Required packages
# ------------------------------------------------------------
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("tidyr")
# install.packages("CSB")
# install.packages("nasapower")

library(ggplot2)
library(nasapower)
library(dplyr)
library(tidyr)
library(CSB)

# ------------------------------------------------------------
# Weekly relative humidity data for Calama, Chile
# Source: NASA POWER
# ------------------------------------------------------------

get_weekly_RH_Calama <- function(
    lon = -68.93,
    lat = -22.46,
    start_date = "2015-01-01",
    end_date = "2025-12-31",
    file = "humedadsemanal.csv"
) {
  
  months_es <- c(
    "enero", "febrero", "marzo", "abril",
    "mayo", "junio", "julio", "agosto",
    "septiembre", "octubre", "noviembre", "diciembre"
  )
  
  rh_data <- get_power(
    community = "AG",
    lonlat = c(lon, lat),
    pars = "RH2M",
    dates = c(start_date, end_date),
    temporal_api = "daily"
  )
  
  data_meses <- rh_data %>%
    mutate(
      semana = case_when(
        DD <= 7  ~ 1L,
        DD <= 14 ~ 2L,
        DD <= 21 ~ 3L,
        TRUE     ~ 4L
      ),
      mes = factor(MM, levels = 1:12, labels = months_es)
    ) %>%
    group_by(YEAR, mes, semana) %>%
    summarise(
      RH_promedio = mean(RH2M, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(YEAR, mes, semana) %>%
    group_by(mes) %>%
    mutate(id = row_number()) %>%
    ungroup() %>%
    select(id, mes, RH_promedio) %>%
    pivot_wider(
      names_from = mes,
      values_from = RH_promedio
    ) %>%
    arrange(id)
  
  write.csv(data_meses, file, row.names = FALSE)
  
  data_meses
}

data_meses <- get_weekly_RH_Calama()

# ------------------------------------------------------------
# Fit competing models
# ------------------------------------------------------------

x <- data_meses$marzo / 100

cmp <- compareCSBmodels(x)

summary(cmp)

# Comparison of fitted models for unit data
# -----------------------------------------
#
# Call:
# compareCSBmodels(x = x)
#
# Sample size: 44 
#
# Model comparison:
#  Model npar  logLik      AIC      BIC convergence
#    CSB    2 44.0968 -84.1936 -80.6252           0
#      B    2 43.5103 -83.0206 -79.4522           0
#      K    2 42.7512 -81.5024 -77.9340           0
#
# Parameter estimates:
#  Model parameter estimate      se
#      B    shape1   9.2829  1.9479
#      B    shape2  17.3724  3.6916
#      K         a   4.1763  0.5089
#      K         b  53.9623 25.1703
#    CSB     shape  20.6363  9.3062
#    CSB         q   0.6991  0.9337
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
# gofCSB_boot(x = x, B = 1000, seed = 2025)
#
# Sample size: 44 
# Bootstrap replicates: 1000 
#
# Parameter estimates:
#  Parameter Estimate
#      shape  20.6363
#          q   0.6991
#
# Fit statistics:
#   logLik      AIC      BIC
#  44.0968 -84.1936 -80.6252
#
# Goodness-of-fit statistics:
#              Test Statistic p_value
#  Anderson-Darling    0.6086   0.224
#  Cramer-von Mises    0.0836   0.266

1 - pchisq((0.6991/0.9337)^2,1)

# 0.454013 
# The null hypothesis H_0: q = 1 is not rejected at conventional 
# significance levels

fit_q1 <- fitCSB_mle_fixedq(x, q = 1, multistart = TRUE)
summary(fit_q1)

# Maximum likelihood estimation for the CSB distribution
# ------------------------------------------------------
# 
# Call:
# fitCSB_mle_fixedq(x = x, q = 1, multistart = TRUE)
# 
# Parameter estimates:
#  Parameter Estimate     SE
#      shape  20.6099 9.3425
#          q   1.0000     NA
# 
# Fit statistics:
#  logLik      AIC      BIC method
#  44.045 -86.0901 -84.3059   BFGS
# 
# Convergence: 0 

# ------------------------------------------------------------
# Fitted histogram 
# ------------------------------------------------------------

dens_grid <- seq(0.01, 0.75, length.out = 600)

dens_df <- data.frame(
  x = dens_grid,
  B = dbeta(dens_grid, shape1 = 9.2829, shape2 = 17.3724),
  K = CSB:::dKum(dens_grid, a = 4.1763, b = 53.9623),
  CSB = dCSB(dens_grid, shape = 20.6363, q = 1)
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
    breaks = seq(0.15, 0.56, by = 0.0312),
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
    x = "Relative humidity",
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
# Empirical CDF vs CSB CDF
# ------------------------------------------------------------

shape_csb <- 20.6363
q_csb     <- 0.6991

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
    x = "Relative humidity",
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

