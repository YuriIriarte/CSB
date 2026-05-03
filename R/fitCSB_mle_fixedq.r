#' Maximum Likelihood Estimation for the CSB Distribution with Fixed q
#'
#' Estimates the shape parameter of the Contracted Symmetric Beta (CSB)
#' distribution by maximum likelihood when the contraction parameter \eqn{q}
#' is assumed to be known and fixed.
#'
#' The log-likelihood is evaluated using a stable log-density implementation.
#' The optimization is performed on the logarithmic scale to ensure positivity
#' of the shape parameter. Multiple starting values can be used when
#' \code{multistart = TRUE}.
#'
#' @param x Numeric vector of observations in \eqn{(0,1)}.
#' @param q Positive numeric value specifying the fixed contraction parameter.
#' @param start Optional numeric starting value for the shape parameter.
#' @param starts Optional list of numeric starting values.
#' @param methods Character vector with optimization methods passed to
#' \code{\link[stats]{optim}}.
#' @param multistart Logical. If \code{TRUE}, a deterministic grid of starting
#' values is used.
#' @param lower_shape,upper_shape Lower and upper bounds for the shape parameter.
#' @param control List of control parameters passed to \code{\link[stats]{optim}}.
#'
#' @return An object of class \code{"fitCSB_mle_fixedq"}, which is a list containing:
#' \itemize{
#'   \item \code{par}: estimated parameters, including \code{shape} and the fixed \code{q}.
#'   \item \code{logLik}: maximized log-likelihood.
#'   \item \code{se}: approximate standard error for the shape parameter.
#'   \item \code{convergence}: convergence code returned by \code{\link[stats]{optim}}.
#'   \item \code{method}: optimization method associated with the best fit.
#'   \item \code{start}: starting value associated with the best fit.
#'   \item \code{message}: optimizer message.
#'   \item \code{success}: logical indicator of successful convergence.
#'   \item \code{AIC}: Akaike information criterion (with one estimated parameter).
#'   \item \code{BIC}: Bayesian information criterion.
#'   \item \code{q_fixed}: value of the fixed contraction parameter.
#'   \item \code{all_fits}: list with all optimization attempts.
#'   \item \code{x}: analyzed sample.
#'   \item \code{call}: matched function call.
#' }
#'
#' @details
#' This function is useful when the contraction parameter \eqn{q} is known
#' or fixed a priori. Only the shape parameter is estimated, which may improve
#' numerical stability and reduce computational cost compared to the full
#' maximum likelihood estimation.
#'
#' The optimization is carried out on the log-scale to enforce positivity
#' of the shape parameter. Standard errors are computed using the observed
#' Hessian and the delta method.
#'
#' @examples
#' set.seed(123)
#' x <- rCSB(n = 200, shape = 2, q = 1)
#'
#' fit <- fitCSB_mle_fixedq(x, q = 1)
#' summary(fit)
#'
#' @export
fitCSB_mle_fixedq <- function(x,
                              q,
                              start = NULL,
                              starts = NULL,
                              methods = c("L-BFGS-B", "BFGS", "Nelder-Mead"),
                              multistart = TRUE,
                              lower_shape = 1e-6,
                              upper_shape = 1e6,
                              control = list(),
                              suppress_warnings = TRUE) {
  
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  
  if (length(x) < 2L) {
    stop("x must have length >= 2.", call. = FALSE)
  }
  
  if (any(x <= 0 | x >= 1)) {
    stop("x must lie strictly in (0,1).", call. = FALSE)
  }
  
  if (!is.numeric(q) || length(q) != 1L || !is.finite(q) || q <= 0) {
    stop("q must be a positive finite number.", call. = FALSE)
  }
  
  start_list <- list()
  
  if (!is.null(start)) {
    start_list[[length(start_list) + 1L]] <- start
  }
  
  if (isTRUE(multistart)) {
    grid_shape <- c(0.25, 0.5, 1, 2, 5, 10, 20)
    start_list <- c(start_list, as.list(grid_shape))
  } else if (length(start_list) == 0L) {
    start_list <- list(1)
  }
  
  if (!is.null(starts)) {
    if (!is.list(starts)) {
      stop("starts must be a list of numeric starting values.", call. = FALSE)
    }
    start_list <- c(start_list, starts)
  }
  
  clean_start <- function(s) {
    s <- as.numeric(s)[1]
    
    if (!is.finite(s) || s <= 0) {
      return(NULL)
    }
    
    min(max(s, lower_shape), upper_shape)
  }
  
  start_list <- lapply(start_list, clean_start)
  start_list <- start_list[!vapply(start_list, is.null, logical(1))]
  start_list <- unique(round(unlist(start_list), 10))
  
  if (length(start_list) == 0L) {
    stop("No valid starting values available.", call. = FALSE)
  }
  
  nll_eta <- function(eta) {
    shape <- exp(eta[1])
    
    ll <- sum(log_dCSB_core(x, shape = shape, q = q))
    
    if (!is.finite(ll)) {
      return(1e12)
    }
    
    -ll
  }
  
  fit_one <- function(st, method) {
    
    eta0 <- log(st)
    
    res <- tryCatch(
      {
        if (method == "L-BFGS-B") {
          stats::optim(
            par = eta0,
            fn = nll_eta,
            method = method,
            lower = log(lower_shape),
            upper = log(upper_shape),
            hessian = TRUE,
            control = control
          )
        } else {
          stats::optim(
            par = eta0,
            fn = nll_eta,
            method = method,
            hessian = TRUE,
            control = control
          )
        }
      },
      error = function(e) e
    )
    
    if (inherits(res, "error")) {
      return(list(
        par = c(shape = NA_real_, q = q),
        logLik = -Inf,
        se = c(shape = NA_real_, q = NA_real_),
        convergence = NA_integer_,
        method = method,
        start = c(shape = st, q = q),
        message = conditionMessage(res),
        success = FALSE,
        AIC = Inf,
        BIC = Inf
      ))
    }
    
    shape_hat <- exp(res$par[1])
    par_hat <- c(shape = shape_hat, q = q)
    
    logLik <- -res$value
    
    se <- c(shape = NA_real_, q = NA_real_)
    
    H <- try(as.matrix(res$hessian), silent = TRUE)
    
    if (!inherits(H, "try-error") &&
        all(dim(H) == c(1, 1)) &&
        all(is.finite(H))) {
      
      V <- try(solve(H), silent = TRUE)
      
      if (!inherits(V, "try-error") &&
          all(is.finite(V))) {
        se_eta <- sqrt(pmax(diag(V), 0))
        se["shape"] <- shape_hat * se_eta
      }
    }
    
    list(
      par = par_hat,
      logLik = logLik,
      se = se,
      convergence = res$convergence,
      method = method,
      start = c(shape = st, q = q),
      message = res$message,
      success = is.finite(logLik) &&
        is.finite(shape_hat) &&
        shape_hat > 0 &&
        res$convergence == 0,
      AIC = -2 * logLik + 2,
      BIC = -2 * logLik + log(length(x)) * 1
    )
  }
  
  all_fits <- list()
  id <- 1L
  
  for (st in start_list) {
    for (meth in methods) {
      if (isTRUE(suppress_warnings)) {
        all_fits[[id]] <- suppressWarnings(fit_one(st, meth))
      } else {
        all_fits[[id]] <- fit_one(st, meth)
      }
      id <- id + 1L
    }
  }
  
  valid <- vapply(all_fits, function(z) is.finite(z$logLik), logical(1))
  
  if (!any(valid)) {
    stop("All optimization attempts failed.", call. = FALSE)
  }
  
  ord <- order(
    vapply(all_fits, function(z) z$logLik, numeric(1)),
    decreasing = TRUE
  )
  
  all_fits <- all_fits[ord]
  best <- all_fits[[1]]
  
  out <- list(
    par = best$par,
    logLik = best$logLik,
    se = best$se,
    convergence = best$convergence,
    method = best$method,
    start = best$start,
    message = best$message,
    success = best$success,
    AIC = best$AIC,
    BIC = best$BIC,
    q_fixed = q,
    suppress_warnings = suppress_warnings,
    all_fits = all_fits,
    x = x,
    call = match.call()
  )
  
  class(out) <- c("fitCSB_mle_fixedq", "fitCSB_mle", "list")
  out
}