#' @export
summary.compareCSBmodels <- function(object, digits = 4, ...) {
  
  out <- list(
    call = object$call,
    n = object$n,
    comparison = object$comparison,
    estimates = object$estimates,
    best_AIC = object$best_AIC,
    best_BIC = object$best_BIC
  )
  
  class(out) <- c("summary.compareCSBmodels", "list")
  out
}

#' @export
print.summary.compareCSBmodels <- function(x, digits = 4, ...) {
  
  cat("\nComparison of fitted models for unit data\n")
  cat("-----------------------------------------\n\n")
  
  if (!is.null(x$call)) {
    cat("Call:\n")
    print(x$call)
    cat("\n")
  }
  
  cat("Sample size:", x$n, "\n\n")
  
  cat("Model comparison:\n")
  
  comp <- x$comparison
  
  if ("success" %in% names(comp)) {
    comp$success <- NULL
  }
  
  comp$logLik <- round(comp$logLik, digits)
  comp$AIC    <- round(comp$AIC, digits)
  comp$BIC    <- round(comp$BIC, digits)
  
  print(comp, row.names = FALSE)
  
  cat("\nParameter estimates:\n")
  
  est <- x$estimates
  est$estimate <- round(est$estimate, digits)
  est$se       <- round(est$se, digits)
  
  print(est, row.names = FALSE)
  
  cat("\nBest model according to AIC:", x$best_AIC, "\n")
  cat("Best model according to BIC:", x$best_BIC, "\n")
  
  invisible(x)
}