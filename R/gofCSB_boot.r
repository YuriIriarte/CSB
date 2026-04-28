#' Bootstrap Goodness-of-Fit Tests for the CSB Distribution
#'
#' Performs bootstrap goodness-of-fit tests for the Contracted Symmetric Beta
#' (CSB) distribution.
#'
#' The function fits the CSB distribution by maximum likelihood, computes the
#' Anderson-Darling and Cramer-von Mises statistics, and obtains bootstrap
#' p-values using a parametric bootstrap procedure.
#'
#' @param x Numeric vector of observations in \eqn{(0,1)}.
#' @param B Number of bootstrap replicates.
#' @param seed Optional seed for reproducibility.
#' @param eps Small positive constant used to avoid evaluating probabilities
#' exactly at 0 or 1.
#' @param verbose Logical. If \code{TRUE}, progress information is printed.
#'
#' @return An object of class \code{"gofCSB_boot"}, which is a list containing:
#' \itemize{
#'   \item \code{par_hat}: maximum likelihood estimates of the CSB parameters.
#'   \item \code{logLik}: maximized log-likelihood.
#'   \item \code{AIC}: Akaike information criterion.
#'   \item \code{BIC}: Bayesian information criterion.
#'   \item \code{AD}: list containing the Anderson-Darling statistic and
#'   bootstrap p-value.
#'   \item \code{CvM}: list containing the Cramer-von Mises statistic and
#'   bootstrap p-value.
#'   \item \code{B}: number of bootstrap replicates.
#'   \item \code{n}: sample size.
#'   \item \code{call}: matched function call.
#' }
#'
#' @details
#' For each bootstrap replicate, a sample is generated from the fitted CSB
#' distribution, the CSB model is refitted, and the goodness-of-fit statistics
#' are recomputed. The p-values are obtained as the proportion of bootstrap
#' statistics greater than or equal to the observed statistics.
#'
#' The printed summary of the returned object can be obtained using
#' \code{summary()}.
#'
#' @examples
#' set.seed(123)
#' x <- rCSB(n = 200, shape = 2, q = 5)
#'
#' gof <- gofCSB_boot(x, B = 49, seed = 123)
#' summary(gof)
#'
#' @export
gofCSB_boot <- function(x,
                        B = 999,
                        seed = 2026,
                        eps = 1e-10,
                        verbose = TRUE) {
  
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  
  if (length(x) < 2L) {
    stop("x must have length >= 2.", call. = FALSE)
  }
  if (any(x <= 0 | x >= 1)) {
    stop("All observations must lie in (0,1).", call. = FALSE)
  }
  
  n <- length(x)
  
  AD_stat <- function(x, par) {
    x <- sort(x)
    n <- length(x)
    
    u <- pCSB(x, shape = par["shape"], q = par["q"])
    u <- pmin(pmax(u, eps), 1 - eps)
    
    -n - mean((2 * seq_len(n) - 1) *
                (log(u) + log(1 - rev(u))))
  }
  
  CvM_stat <- function(x, par) {
    x <- sort(x)
    n <- length(x)
    
    u <- pCSB(x, shape = par["shape"], q = par["q"])
    u <- pmin(pmax(u, eps), 1 - eps)
    
    1 / (12 * n) + sum((u - (2 * seq_len(n) - 1) / (2 * n))^2)
  }
  
  fit <- fitCSB_mle(x)
  par_hat <- fit$par
  
  AD_obs  <- AD_stat(x, par_hat)
  CvM_obs <- CvM_stat(x, par_hat)
  
  AD_boot  <- rep(NA_real_, B)
  CvM_boot <- rep(NA_real_, B)
  conv_boot <- rep(FALSE, B)
  
  set.seed(seed)
  
  for (b in seq_len(B)) {
    
    xb <- rCSB(
      n = n,
      shape = par_hat["shape"],
      q = par_hat["q"]
    )
    
    fit_b <- try(fitCSB_mle(xb), silent = TRUE)
    
    if (!inherits(fit_b, "try-error") &&
        isTRUE(fit_b$success) &&
        !is.null(fit_b$par) &&
        all(is.finite(fit_b$par))) {
      
      par_b <- fit_b$par
      
      AD_boot[b]  <- AD_stat(xb, par_b)
      CvM_boot[b] <- CvM_stat(xb, par_b)
      conv_boot[b] <- TRUE
    }
    
    if (verbose && b %% 100 == 0) {
      cat("Bootstrap replication:", b, "of", B, "\n")
    }
  }
  
  AD_boot_valid  <- AD_boot[conv_boot]
  CvM_boot_valid <- CvM_boot[conv_boot]
  
  B_valid <- sum(conv_boot)
  
  if (B_valid == 0L) {
    p_AD <- NA_real_
    p_CvM <- NA_real_
  } else {
    p_AD  <- (1 + sum(AD_boot_valid >= AD_obs)) / (B_valid + 1)
    p_CvM <- (1 + sum(CvM_boot_valid >= CvM_obs)) / (B_valid + 1)
  }
  
  out <- list(
    model = "CSB",
    n = n,
    B = B,
    B_valid = B_valid,
    seed = seed,
    par_hat = par_hat,
    logLik = fit$logLik,
    AIC = fit$AIC,
    BIC = fit$BIC,
    AD = list(
      statistic = AD_obs,
      p_value = p_AD
    ),
    CvM = list(
      statistic = CvM_obs,
      p_value = p_CvM
    ),
    bootstrap = data.frame(
      AD = AD_boot,
      CvM = CvM_boot,
      convergence = conv_boot
    ),
    fit = fit
  )
  
  class(out) <- c("gofCSB_boot", "list")
  return(out)
}