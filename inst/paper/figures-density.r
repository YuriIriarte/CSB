# ------------------------------------------------------------
# Required packages
# ------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(CSB)
library(grid)

# ------------------------------------------------------------
# Global settings
# ------------------------------------------------------------

x_values <- seq(0.001, 0.999, length.out = 1000)

colors <- c("black", "#D55E00", "#0072B2", "#CC79A7")
linetypes <- c("solid", "dashed", "dotted", "dotdash")

# ------------------------------------------------------------
# Theme (consistente con boxplots del paper)
# ------------------------------------------------------------

theme_csb <- theme_minimal(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid.major = element_line(color = "grey80", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12),
    
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.key.width = unit(1.1, "cm"),
    
    legend.position = c(0.8, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_rect(
      fill = "white",
      color = "grey60",
      linewidth = 0.4
    )
  )

# ------------------------------------------------------------
# Figure 1a. alpha = 0.5, q = {1, 5, 10, 25}
# ------------------------------------------------------------

fig1_data <- bind_rows(
  data.frame(x = x_values, density = dCSB(x_values, 0.5, 1),  label = "1"),
  data.frame(x = x_values, density = dCSB(x_values, 0.5, 5),  label = "5"),
  data.frame(x = x_values, density = dCSB(x_values, 0.5, 10), label = "10"),
  data.frame(x = x_values, density = dCSB(x_values, 0.5, 25), label = "25")
) %>%
  mutate(label = factor(label, levels = c("1", "5", "10", "25")))

ggplot(
  fig1_data,
  aes(x = x, y = density, color = label, linetype = label)
  ) +
  geom_line(linewidth = 0.9) +
  coord_cartesian(ylim = c(0, 6)) +
  scale_color_manual(
    values = colors,
    name = expression(q)
  ) +
  scale_linetype_manual(
    values = linetypes,
    name = expression(q)
  ) +
  labs(
    x = expression(t),
    y = "Density function"
  ) +
  theme_csb

# ------------------------------------------------------------
# Figure 1b. alpha = 1, q = {1, 5, 10, 25}
# ------------------------------------------------------------

fig1_data <- bind_rows(
  data.frame(x = x_values, density = dCSB(x_values, 1, 1),  label = "1"),
  data.frame(x = x_values, density = dCSB(x_values, 1, 5),  label = "5"),
  data.frame(x = x_values, density = dCSB(x_values, 1, 10), label = "10"),
  data.frame(x = x_values, density = dCSB(x_values, 1, 25), label = "25")
) %>%
  mutate(label = factor(label, levels = c("1", "5", "10", "25")))

ggplot(
  fig1_data,
  aes(x = x, y = density, color = label, linetype = label)
  ) +
  geom_line(linewidth = 0.9) +
  coord_cartesian(ylim = c(0, 2)) +
  scale_color_manual(
    values = colors,
    name = expression(q)
  ) +
  scale_linetype_manual(
    values = linetypes,
    name = expression(q)
  ) +
  labs(
    x = expression(t),
    y = "Density function"
  ) +
  theme_csb

# ------------------------------------------------------------
# Figure 1c. alpha = 5, q = {5, 5, 10, 25}
# ------------------------------------------------------------

fig1_data <- bind_rows(
  data.frame(x = x_values, density = dCSB(x_values, 5, 1),  label = "1"),
  data.frame(x = x_values, density = dCSB(x_values, 5, 5),  label = "5"),
  data.frame(x = x_values, density = dCSB(x_values, 5, 10), label = "10"),
  data.frame(x = x_values, density = dCSB(x_values, 5, 25), label = "25")
) %>%
  mutate(label = factor(label, levels = c("1", "5", "10", "25")))

ggplot(
  fig1_data,
  aes(x = x, y = density, color = label, linetype = label)
  ) +
  geom_line(linewidth = 0.9) +
  coord_cartesian(ylim = c(0, 5)) +
  scale_color_manual(
    values = colors,
    name = expression(q)
  ) +
  scale_linetype_manual(
    values = linetypes,
    name = expression(q)
  ) +
  labs(
    x = expression(t),
    y = "Density function"
  ) +
  theme_csb



