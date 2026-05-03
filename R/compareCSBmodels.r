#' Compare Fitted Models for Unit Data
#'
#' Fits and compares the Beta, Kumaraswamy, and Contracted Symmetric Beta
#' (CSB) distributions for data supported on the unit interval.
#'
#' The comparison is based on the maximized log-likelihood, Akaike information
#' criterion (AIC), and Bayesian information criterion (BIC). Parameter
#' estimates and approximate standard errors are also returned for each fitted
#' model.
#'
#' @param x Numeric vector of observations in \eqn{(0,1)}.
#' @param suppress_warnings Logical. If \code{TRUE}, warnings generated during
#' model fitting are suppressed. The default is \code{TRUE}, which is convenient
#' for clean output in applications and simulation studies. Use
#' \code{suppress_warnings = FALSE} to inspect possible numerical warnings.
#'
#' @return An object of class \code{"compareCSBmodels"}, which is a list containing:
#' \itemize{
#'   \item \code{comparison}: data frame with model names, number of parameters,
#'   log-likelihood, AIC, BIC, and convergence codes.
#'   \item \code{estimates}: data frame with parameter estimates and approximate
#'   standard errors.
#'   \item \code{fits}: list containing the fitted model objects.
#'   \item \code{n}: sample size.
#'   \item \code{best_AIC}: model selected by AIC.
#'   \item \code{best_BIC}: model selected by BIC.
#'   \item \code{call}: matched function call.
#' }
#'
#' @details
#' The printed summary of the returned object can be obtained using
#' \code{summary()}.
#'
#' @examples
#' set.seed(123)
#' x <- rCSB(n = 200, shape = 2, q = 5)
#'
#' cmp <- compareCSBmodels(x)
#' summary(cmp)
#'
#' # To inspect numerical warnings, use:
#' cmp_check <- compareCSBmodels(x, suppress_warnings = FALSE)
#'
#' @export
compareCSBmodels <- function(x, suppress_warnings = TRUE) {
  
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  
  if (length(x) < 2L) {
    stop("x must have length >= 2.", call. = FALSE)
  }
  
  if (any(x <= 0 | x >= 1)) {
    stop("All observations must lie in (0,1).", call. = FALSE)
  }
  
  fit_try <- function(expr) {
    if (isTRUE(suppress_warnings)) {
      suppressWarnings(try(expr, silent = TRUE))
    } else {
      try(expr, silent = TRUE)
    }
  }
  
  n <- length(x)
  fits <- list()
  rows <- list()
  id <- 1L
  
  # ----------------------------------------------------------
  # Beta
  # ----------------------------------------------------------
  
  fit_beta <- fit_try(fitBeta_mle(x))
  
  if (!inherits(fit_beta, "try-error")) {
    
    fits$beta <- fit_beta
    
    rows[[id]] <- data.frame(
      Model = "B",
      npar = fit_beta$k,
      logLik = fit_beta$logLik,
      AIC = fit_beta$AIC,
      BIC = fit_beta$BIC,
      convergence = fit_beta$convergence,
      success = fit_beta$success,
      stringsAsFactors = FALSE
    )
    
    id <- id + 1L
  }
  
  # ----------------------------------------------------------
  # Kumaraswamy
  # ----------------------------------------------------------
  
  fit_kum <- fit_try(fitKum_mle(x))
  
  if (!inherits(fit_kum, "try-error")) {
    
    fits$kumaraswamy <- fit_kum
    
    rows[[id]] <- data.frame(
      Model = "K",
      npar = fit_kum$k,
      logLik = fit_kum$logLik,
      AIC = fit_kum$AIC,
      BIC = fit_kum$BIC,
      convergence = fit_kum$convergence,
      success = fit_kum$success,
      stringsAsFactors = FALSE
    )
    
    id <- id + 1L
  }
  
  # ----------------------------------------------------------
  # CSB
  # ----------------------------------------------------------
  
  fit_csb <- fit_try(fitCSB_mle(x))
  
  if (!inherits(fit_csb, "try-error")) {
    
    fits$csb <- fit_csb
    
    rows[[id]] <- data.frame(
      Model = "CSB",
      npar = 2L,
      logLik = fit_csb$logLik,
      AIC = fit_csb$AIC,
      BIC = fit_csb$BIC,
      convergence = fit_csb$convergence,
      success = fit_csb$success,
      stringsAsFactors = FALSE
    )
    
    id <- id + 1L
  }
  
  if (length(rows) == 0L) {
    stop("No model could be fitted.", call. = FALSE)
  }
  
  comparison <- do.call(rbind, rows)
  rownames(comparison) <- NULL
  
  comparison <- comparison[order(comparison$AIC), ]
  rownames(comparison) <- NULL
  
  estimates <- list()
  
  if (!is.null(fits$beta)) {
    estimates[[length(estimates) + 1L]] <- data.frame(
      Model = "B",
      parameter = names(fits$beta$par),
      estimate = as.numeric(fits$beta$par),
      se = as.numeric(fits$beta$se),
      stringsAsFactors = FALSE
    )
  }
  
  if (!is.null(fits$kumaraswamy)) {
    estimates[[length(estimates) + 1L]] <- data.frame(
      Model = "K",
      parameter = names(fits$kumaraswamy$par),
      estimate = as.numeric(fits$kumaraswamy$par),
      se = as.numeric(fits$kumaraswamy$se),
      stringsAsFactors = FALSE
    )
  }
  
  if (!is.null(fits$csb)) {
    estimates[[length(estimates) + 1L]] <- data.frame(
      Model = "CSB",
      parameter = names(fits$csb$par),
      estimate = as.numeric(fits$csb$par),
      se = as.numeric(fits$csb$se),
      stringsAsFactors = FALSE
    )
  }
  
  estimates <- do.call(rbind, estimates)
  rownames(estimates) <- NULL
  
  out <- list(
    comparison = comparison,
    estimates = estimates,
    fits = fits,
    n = n,
    best_AIC = comparison$Model[which.min(comparison$AIC)],
    best_BIC = comparison$Model[which.min(comparison$BIC)],
    suppress_warnings = suppress_warnings,
    call = match.call()
  )
  
  class(out) <- c("compareCSBmodels", "list")
  out
}