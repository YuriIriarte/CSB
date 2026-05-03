# ------------------------------------------------------------
# Required packages
# ------------------------------------------------------------

library(CSB)
library(ggplot2)
library(dplyr)
library(grid)
library(plot3D)

# ------------------------------------------------------------
# Global graphical settings
# ------------------------------------------------------------

colors5 <- c("black", "#D55E00", "#0072B2", "#CC79A7", "#009E73")
linetypes5 <- c("solid", "dashed", "dotted", "dotdash", "longdash")

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
    legend.background = element_rect(
      fill = "white",
      color = "grey60",
      linewidth = 0.4
    )
  )

# ------------------------------------------------------------
# Figure: Skewness
# q = {5, 8, 10, 15, 20}
# ------------------------------------------------------------

shape_values <- seq(3, 150, length.out = 300)
q_skew <- c(5, 8, 10, 15, 20)

data_skew <- expand.grid(
  shape = shape_values,
  q = q_skew
) %>%
  rowwise() %>%
  mutate(
    skewness = skewCSB(shape = shape, q = q)
  ) %>%
  ungroup() %>%
  mutate(
    q_label = factor(q, levels = q_skew)
  )

ggplot(
  data_skew,
  aes(x = shape, y = skewness, color = q_label, linetype = q_label)
  ) +
  geom_line(linewidth = 0.9) +
  coord_cartesian(ylim = c(0, 1.7)) +
  scale_color_manual(
    values = colors5,
    breaks = q_skew,
    name = expression(q)
  ) +
  scale_linetype_manual(
    values = linetypes5,
    breaks = q_skew,
    name = expression(q)
  ) +
  labs(
    x = expression(alpha),
    y = "Skewness"
  ) +
  theme_csb +
  theme(
    legend.position = c(0.97, 0.05),
    legend.justification = c(1, 0)
  )

# ------------------------------------------------------------
# Figure: Kurtosis
# q = {5, 10, 25, 50, 100}
# ------------------------------------------------------------

shape_values <- seq(4, 150, length.out = 300)
q_kurt <- c(5, 10, 25, 50, 100)

data_kurt <- expand.grid(
  shape = shape_values,
  q = q_kurt
) %>%
  rowwise() %>%
  mutate(
    kurtosis = kurtCSB(
      shape = shape,
      q = q,
      excess = FALSE
    )
  ) %>%
  ungroup() %>%
  mutate(
    q_label = factor(q, levels = q_kurt)
  )

ggplot(
  data_kurt,
  aes(x = shape, y = kurtosis, color = q_label, linetype = q_label)
  ) +
  geom_line(linewidth = 0.9) +
  coord_cartesian(ylim = c(2.5, 8)) +
  scale_color_manual(
    values = colors5,
    breaks = q_kurt,
    name = expression(q)
  ) +
  scale_linetype_manual(
    values = linetypes5,
    breaks = q_kurt,
    name = expression(q)
  ) +
  labs(
    x = expression(alpha),
    y = "Kurtosis"
  ) +
  theme_csb +
  theme(
    legend.position = c(0.03, 0.97),
    legend.justification = c(0, 1)
  )

# ------------------------------------------------------------
# Skewness 2D surface (ggplot version)
# ------------------------------------------------------------

q_grid <- seq(3, 50, length.out = 80)
shape_grid <- seq(3, 50, length.out = 80)

data_skew3D <- expand.grid(
  q = q_grid,
  shape = shape_grid
) %>%
  rowwise() %>%
  mutate(
    skewness = skewCSB(shape = shape, q = q)
  ) %>%
  ungroup()

ggplot(data_skew3D, aes(x = q, y = shape)) +
  geom_tile(aes(fill = skewness)) +
  geom_contour(
    aes(z = skewness),
    color = "black",
    linewidth = 0.3
  ) +
  scale_fill_gradient(
    low = "white",
    high = "grey40",
    name = "Skewness"
  ) +
  labs(
    x = "q",
    y = expression(alpha)
  ) +
  theme_csb +
  theme(
    legend.position = c(0.95, 0.25),
    legend.justification = c(1, 0.5)
  )

# ------------------------------------------------------------
# Kurtosis 2D surface (ggplot version)
# ------------------------------------------------------------

q_grid <- seq(4, 50, length.out = 80)
shape_grid <- seq(4, 50, length.out = 80)

data_kurt3D <- expand.grid(
  q = q_grid,
  shape = shape_grid
) %>%
  rowwise() %>%
  mutate(
    kurtosis = kurtCSB(shape = shape, q = q, excess = FALSE)
  ) %>%
  ungroup()

ggplot(data_kurt3D, aes(x = q, y = shape)) +
  geom_tile(aes(fill = kurtosis)) +
  geom_contour(
    aes(z = kurtosis),
    color = "black",
    linewidth = 0.3
  ) +
  scale_fill_gradient(
    low = "white",
    high = "grey40",
    name = "Kurtosis"
  ) +
  labs(
    x = "q",
    y = expression(alpha)
  ) +
  theme_csb +
  theme(
    legend.position = c(0.95, 0.25),
    legend.justification = c(1, 0.5)
  )


