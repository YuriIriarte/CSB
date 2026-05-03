# ============================================================
# Monte Carlo study for the CSB distribution
# ============================================================

# Required packages
library(CSB)
library(dplyr)
library(tidyr)
library(ggplot2)

# ------------------------------------------------------------
# One Monte Carlo replicate
# ------------------------------------------------------------

mc_one_csb <- function(n, shape_true, q_true) {
  
  x <- rCSB(n = n, shape = shape_true, q = q_true)
  
  fit <- try(
    fitCSB_mle(
      x,
      method = "L-BFGS-B",
      multistart = FALSE
    ),
    silent = TRUE
  )
  
  if (inherits(fit, "try-error") ||
      is.null(fit$par) ||
      any(!is.finite(fit$par)) ||
      fit$convergence != 0) {
    
    return(data.frame(
      shape_true = shape_true,
      q_true = q_true,
      n = n,
      shape_hat = NA_real_,
      q_hat = NA_real_,
      se_shape = NA_real_,
      se_q = NA_real_,
      cp_shape = NA_real_,
      cp_q = NA_real_,
      success = FALSE
    ))
  }
  
  shape_hat <- unname(fit$par["shape"])
  q_hat     <- unname(fit$par["q"])
  
  se_shape <- unname(fit$se["shape"])
  se_q     <- unname(fit$se["q"])
  
  cp_shape <- NA_real_
  cp_q     <- NA_real_
  
  if (is.finite(se_shape) && se_shape > 0) {
    li_shape <- shape_hat - 1.96 * se_shape
    ls_shape <- shape_hat + 1.96 * se_shape
    cp_shape <- as.numeric(shape_true >= li_shape && shape_true <= ls_shape)
  }
  
  if (is.finite(se_q) && se_q > 0) {
    li_q <- q_hat - 1.96 * se_q
    ls_q <- q_hat + 1.96 * se_q
    cp_q <- as.numeric(q_true >= li_q && q_true <= ls_q)
  }
  
  data.frame(
    shape_true = shape_true,
    q_true = q_true,
    n = n,
    shape_hat = shape_hat,
    q_hat = q_hat,
    se_shape = se_shape,
    se_q = se_q,
    cp_shape = cp_shape,
    cp_q = cp_q,
    success = TRUE
  )
}

# ------------------------------------------------------------
# Full Monte Carlo study
# ------------------------------------------------------------

mc_CSB_study <- function(R = 1000,
                         n_values = c(50, 100, 200, 300),
                         shape_values = c(2, 4, 8),
                         q_values = c(1, 2, 3),
                         seed = 123,
                         verbose = TRUE) {
  
  set.seed(seed)
  
  scenarios <- expand.grid(
    shape_true = shape_values,
    q_true = q_values,
    n = n_values
  )
  
  results <- vector("list", nrow(scenarios) * R)
  id <- 1L
  
  for (i in seq_len(nrow(scenarios))) {
    
    shape_true <- scenarios$shape_true[i]
    q_true     <- scenarios$q_true[i]
    n          <- scenarios$n[i]
    
    if (verbose) {
      message(
        "Scenario ", i, "/", nrow(scenarios),
        ": shape = ", shape_true,
        ", q = ", q_true,
        ", n = ", n
      )
    }
    
    for (r in seq_len(R)) {
      
      aux <- mc_one_csb(
        n = n,
        shape_true = shape_true,
        q_true = q_true
      )
      
      aux$replicate <- r
      results[[id]] <- aux
      id <- id + 1L
    }
  }
  
  bind_rows(results) %>%
    select(
      replicate,
      shape_true,
      q_true,
      n,
      shape_hat,
      q_hat,
      se_shape,
      se_q,
      cp_shape,
      cp_q,
      success
    )
}

# ------------------------------------------------------------
# Summary table: mean, bias, RMSE, coverage probability
# ------------------------------------------------------------

summarize_CSB_mc <- function(mc_results) {
  
  mc_results %>%
    group_by(shape_true, q_true, n) %>%
    summarise(
      mean_shape = mean(shape_hat, na.rm = TRUE),
      mean_q     = mean(q_hat, na.rm = TRUE),
      
      Bias_shape = mean(shape_hat - shape_true, na.rm = TRUE),
      Bias_q     = mean(q_hat - q_true, na.rm = TRUE),
      
      RMSE_shape = sqrt(mean((shape_hat - shape_true)^2, na.rm = TRUE)),
      RMSE_q     = sqrt(mean((q_hat - q_true)^2, na.rm = TRUE)),
      
      CP_shape   = mean(cp_shape, na.rm = TRUE),
      CP_q       = mean(cp_q, na.rm = TRUE),
      
      SuccessRate = mean(success, na.rm = TRUE),
      .groups = "drop"
    )
}

# ------------------------------------------------------------
# Boxplots of the parameter estimates
# ------------------------------------------------------------

plot_box_csb_free_y <- function(mc_results, alpha_val) {
  
  df <- mc_results %>%
    filter(success == TRUE,
           shape_true == alpha_val) %>%
    select(shape_true, q_true, n, shape_hat, q_hat) %>%
    pivot_longer(
      cols = c(shape_hat, q_hat),
      names_to = "Parameter",
      values_to = "Estimate"
    ) %>%
    mutate(
      Parameter = recode(
        Parameter,
        shape_hat = "hat(alpha)",
        q_hat     = "hat(q)"
      ),
      Panel = paste0(Parameter, "*','~q==", q_true),
      n = factor(n)
    )
  
  ggplot(df, aes(x = n, y = Estimate)) +
    geom_boxplot(
      width = 0.60,
      fill = "white",
      color = "black",
      outlier.size = 1.1,
      outlier.alpha = 0.55
    ) +
    facet_wrap(
      ~ Panel,
      scales = "free_y",
      nrow = 2,
      labeller = label_parsed
    ) +
    labs(
      x = "Sample size (n)",
      y = "Estimated values"
    ) +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "grey85", color = NA),
      strip.text = element_text(face = "bold", size = 14),
      
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      
      panel.grid.major = element_line(color = "grey82"),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(1.2, "lines")
    )
}

# ------------------------------------------------------------
# Run the Monte Carlo study
# ------------------------------------------------------------

mc_res <- mc_CSB_study(
  R = 1000,
  n_values = c(50, 100, 200, 300, 500),
  shape_values = c(1, 4, 8),
  q_values = c(1, 2, 3),
  seed = 123
)

# ------------------------------------------------------------
# Monte Carlo summaries
# ------------------------------------------------------------

mc_summary_full <- summarize_CSB_mc(mc_res)

mc_summary_full %>%
  filter(shape_true == 1) %>%
  mutate(across(where(is.numeric), ~ round(.x, 5)))

mc_summary_full %>%
  filter(shape_true == 4) %>%
  mutate(across(where(is.numeric), ~ round(.x, 5)))

mc_summary_full %>%
  filter(shape_true == 8) %>%
  mutate(across(where(is.numeric), ~ round(.x, 5)))

# ------------------------------------------------------------
# Boxplots
# ------------------------------------------------------------

plot_box_csb_free_y(mc_res, alpha_val = 1)

plot_box_csb_free_y(mc_res, alpha_val = 4)

plot_box_csb_free_y(mc_res, alpha_val = 8)
