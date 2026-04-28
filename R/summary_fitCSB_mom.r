#' @export
summary.fitCSB_mom <- function(object, digits = 4, ...) {
  
  par_tab <- data.frame(
    Parameter = names(object$par),
    Estimate = as.numeric(object$par),
    row.names = NULL
  )
  
  out <- list(
    call = object$call,
    coefficients = par_tab,
    objective = object$objective,
    method = object$method,
    convergence = object$convergence,
    success = object$success
  )
  
  class(out) <- c("summary.fitCSB_mom", "list")
  out
}

#' @export
print.summary.fitCSB_mom <- function(x, digits = 4, ...) {
  
  cat("\n")
  cat("Moment estimation for the CSB distribution\n")
  cat(strrep("-", 42), "\n\n", sep = "")
  
  cat("Call:\n")
  print(x$call)
  
  cat("\nParameter estimates:\n")
  coef_tab <- x$coefficients
  coef_tab$Estimate <- round(coef_tab$Estimate, digits)
  print(coef_tab, row.names = FALSE)
  
  cat("\nObjective function value:",
      round(x$objective, digits), "\n")
  
  cat("Optimization method:",
      x$method, "\n")
  
  cat("Convergence:",
      x$convergence, "\n")
  
  cat("Successful fit:",
      x$success, "\n")
  
  invisible(x)
}