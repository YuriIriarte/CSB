#' @export
summary.fitCSB_mle <- function(object, digits = 4, ...) {
  
  par_tab <- data.frame(
    Parameter = names(object$par),
    Estimate  = as.numeric(object$par),
    SE        = as.numeric(object$se),
    row.names = NULL
  )
  
  fit_tab <- data.frame(
    logLik = object$logLik,
    AIC    = object$AIC,
    BIC    = object$BIC,
    method = object$method,
    stringsAsFactors = FALSE
  )
  
  out <- list(
    call = object$call,
    coefficients = par_tab,
    fit_statistics = fit_tab,
    convergence = object$convergence,
    success = object$success
  )
  
  class(out) <- c("summary.fitCSB_mle", "list")
  out
}

#' @export
print.summary.fitCSB_mle <- function(x, digits = 4, ...) {
  
  cat("\n")
  cat("Maximum likelihood estimation for the CSB distribution\n")
  cat(strrep("-", 54), "\n\n", sep = "")
  
  cat("Call:\n")
  print(x$call)
  
  cat("\nParameter estimates:\n")
  coef_tab <- x$coefficients
  coef_tab$Estimate <- round(coef_tab$Estimate, digits)
  coef_tab$SE <- round(coef_tab$SE, digits)
  print(coef_tab, row.names = FALSE)
  
  cat("\nFit statistics:\n")
  fit_tab <- x$fit_statistics
  fit_tab$logLik <- round(fit_tab$logLik, digits)
  fit_tab$AIC    <- round(fit_tab$AIC, digits)
  fit_tab$BIC    <- round(fit_tab$BIC, digits)
  print(fit_tab, row.names = FALSE)
  
  cat("\nConvergence:",
      x$convergence, "\n")
  
  invisible(x)
}