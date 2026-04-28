#' @export
summary.gofCSB_boot <- function(object, digits = 4, ...) {
  
  par_tab <- data.frame(
    Parameter = names(object$par_hat),
    Estimate = as.numeric(object$par_hat),
    row.names = NULL
  )
  
  fit_tab <- data.frame(
    logLik = object$logLik,
    AIC = object$AIC,
    BIC = object$BIC,
    stringsAsFactors = FALSE
  )
  
  gof_tab <- data.frame(
    Test = c("Anderson-Darling", "Cramer-von Mises"),
    Statistic = c(object$AD$statistic, object$CvM$statistic),
    p_value = c(object$AD$p_value, object$CvM$p_value),
    stringsAsFactors = FALSE
  )
  
  out <- list(
    call = object$call,
    n = object$n,
    B = object$B,
    par_hat = par_tab,
    fit_statistics = fit_tab,
    goodness_of_fit = gof_tab
  )
  
  class(out) <- c("summary.gofCSB_boot", "list")
  out
}

#' @export
print.summary.gofCSB_boot <- function(x, digits = 4, ...) {
  
  cat("\nBootstrap goodness-of-fit for the CSB distribution\n")
  cat("--------------------------------------------------\n\n")
  
  if (!is.null(x$call)) {
    cat("Call:\n")
    print(x$call)
    cat("\n")
  }
  
  cat("Sample size:", x$n, "\n")
  cat("Bootstrap replicates:", x$B, "\n\n")
  
  cat("Parameter estimates:\n")
  par_tab <- x$par_hat
  par_tab$Estimate <- round(par_tab$Estimate, digits)
  print(par_tab, row.names = FALSE)
  
  cat("\nFit statistics:\n")
  fit_tab <- x$fit_statistics
  fit_tab$logLik <- round(fit_tab$logLik, digits)
  fit_tab$AIC <- round(fit_tab$AIC, digits)
  fit_tab$BIC <- round(fit_tab$BIC, digits)
  print(fit_tab, row.names = FALSE)
  
  cat("\nGoodness-of-fit statistics:\n")
  gof_tab <- x$goodness_of_fit
  gof_tab$Statistic <- round(gof_tab$Statistic, digits)
  gof_tab$p_value <- round(gof_tab$p_value, digits)
  print(gof_tab, row.names = FALSE)
  
  invisible(x)
}