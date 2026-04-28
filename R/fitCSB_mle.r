#' Maximum Likelihood Estimation for the Contracted Symmetric Beta Distribution
#'
#' Estimates the parameters of the Contracted Symmetric Beta (CSB) distribution
#' by maximum likelihood.
#'
#' The log-likelihood is evaluated using a stable log-density implementation.
#' The optimization is performed on the logarithmic scale to ensure positivity
#' of the parameters. Moment-based estimates are used as initial values whenever
#' available, and additional deterministic starting values can be used when
#' \code{multistart = TRUE}.
#'
#' @param x Numeric vector of observations in \eqn{(0,1)}.
#' @param start Optional named numeric vector with initial values for
#' \code{shape} and \code{q}.
#' @param starts Optional list of named numeric vectors with additional starting
#' values.
#' @param methods Character vector with optimization methods passed to
#' \code{\link[stats]{optim}}.
#' @param multistart Logical. If \code{TRUE}, a deterministic grid of starting
#' values is used.
#' @param lower_shape,upper_shape Lower and upper bounds for the shape parameter.
#' @param lower_q,upper_q Lower and upper bounds for the contraction parameter.
#' @param control List of control parameters passed to \code{\link[stats]{optim}}.
#'
#' @return An object of class \code{"fitCSB_mle"}, which is a list containing:
#' \itemize{
#'   \item \code{par}: maximum likelihood estimates.
#'   \item \code{logLik}: maximized log-likelihood.
#'   \item \code{se}: approximate standard errors based on the observed Hessian
#'   and the delta method.
#'   \item \code{convergence}: convergence code returned by \code{\link[stats]{optim}}.
#'   \item \code{method}: optimization method associated with the best fit.
#'   \item \code{start}: starting values associated with the best fit.
#'   \item \code{message}: optimizer message.
#'   \item \code{success}: logical indicator of successful convergence.
#'   \item \code{AIC}: Akaike information criterion.
#'   \item \code{BIC}: Bayesian information criterion.
#'   \item \code{all_fits}: list with all optimization attempts.
#'   \item \code{x}: analyzed sample.
#'   \item \code{call}: matched function call.
#' }
#'
#' @examples
#' set.seed(123)
#' x <- rCSB(n = 200, shape = 2, q = 5)
#'
#' fit_mle <- fitCSB_mle(x)
#' summary(fit_mle)
#'
#' @export
fitCSB_mle <- function(x,
                       start = NULL,
                       starts = NULL,
                       methods = c("L-BFGS-B", "BFGS", "Nelder-Mead"),
                       multistart = TRUE,
                       lower_shape = 1e-6,
                       upper_shape = 1e6,
                       lower_q = 1e-6,
                       upper_q = 1e6,
                       control = list()) {
  
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  
  if (length(x) < 2L) {
    stop("x must have length >= 2.", call. = FALSE)
  }
  if (any(x <= 0 | x >= 1)) {
    stop("x must lie strictly in (0,1).", call. = FALSE)
  }
  
  start_list <- list()
  
  mom_fit <- try(fitCSB_mom(x, multistart = multistart), silent = TRUE)
  
  if (!inherits(mom_fit, "try-error") &&
      isTRUE(mom_fit$success) &&
      !is.null(mom_fit$par)) {
    start_list[[length(start_list) + 1L]] <- mom_fit$par
  }
  
  if (!is.null(start)) {
    start_list[[length(start_list) + 1L]] <- start
  }
  
  if (isTRUE(multistart)) {
    
    grid_shape <- c(0.5, 1, 2, 5, 10)
    grid_q     <- c(0.5, 1, 2, 5, 10)
    
    grid_starts <- expand.grid(
      shape = grid_shape,
      q = grid_q
    )
    
    grid_starts <- split(grid_starts, seq_len(nrow(grid_starts)))
    
    grid_starts <- lapply(grid_starts, function(z) {
      c(shape = z$shape, q = z$q)
    })
    
    start_list <- c(start_list, grid_starts)
    
  } else if (length(start_list) == 0L) {
    start_list <- list(c(shape = 1, q = 1))
  }
  
  if (!is.null(starts)) {
    if (!is.list(starts)) {
      stop("starts must be a list of named numeric vectors.", call. = FALSE)
    }
    start_list <- c(start_list, starts)
  }
  
  clean_start <- function(s) {
    
    if (is.null(s) || is.null(names(s)) ||
        !all(c("shape", "q") %in% names(s))) {
      return(NULL)
    }
    
    s <- as.numeric(s[c("shape", "q")])
    names(s) <- c("shape", "q")
    
    if (any(!is.finite(s)) || any(s <= 0)) {
      return(NULL)
    }
    
    s["shape"] <- min(max(s["shape"], lower_shape), upper_shape)
    s["q"]     <- min(max(s["q"], lower_q), upper_q)
    
    s
  }
  
  start_list <- lapply(start_list, clean_start)
  start_list <- start_list[!vapply(start_list, is.null, logical(1))]
  
  if (length(start_list) == 0L) {
    stop("No valid starting values available.", call. = FALSE)
  }
  
  keyfun <- function(s) paste(round(log(s), 6), collapse = "_")
  keys <- vapply(start_list, keyfun, character(1))
  start_list <- start_list[!duplicated(keys)]
  
  nll_eta <- function(eta) {
    
    shape <- exp(eta[1])
    q     <- exp(eta[2])
    
    ll <- sum(log_dCSB_core(x, shape = shape, q = q))
    
    if (!is.finite(ll)) {
      return(1e12)
    }
    
    -ll
  }
  
  fit_one <- function(st, method) {
    
    eta0 <- log(c(st["shape"], st["q"]))
    
    res <- tryCatch(
      {
        if (method == "L-BFGS-B") {
          stats::optim(
            par = eta0,
            fn = nll_eta,
            method = method,
            lower = log(c(lower_shape, lower_q)),
            upper = log(c(upper_shape, upper_q)),
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
        par = c(shape = NA_real_, q = NA_real_),
        logLik = -Inf,
        se = c(shape = NA_real_, q = NA_real_),
        convergence = NA_integer_,
        method = method,
        start = st,
        message = conditionMessage(res),
        success = FALSE,
        AIC = Inf,
        BIC = Inf
      ))
    }
    
    par_hat <- exp(res$par)
    names(par_hat) <- c("shape", "q")
    
    logLik <- -res$value
    
    se <- c(shape = NA_real_, q = NA_real_)
    
    H <- try(as.matrix(res$hessian), silent = TRUE)
    
    if (!inherits(H, "try-error") &&
        all(dim(H) == c(2, 2)) &&
        all(is.finite(H))) {
      
      V <- try(solve(H), silent = TRUE)
      
      if (!inherits(V, "try-error") &&
          all(is.finite(V))) {
        
        se_eta <- sqrt(pmax(diag(V), 0))
        se <- par_hat * se_eta
        names(se) <- c("shape", "q")
      }
    }
    
    list(
      par = par_hat,
      logLik = logLik,
      se = se,
      convergence = res$convergence,
      method = method,
      start = st,
      message = res$message,
      success = is.finite(logLik) &&
        all(is.finite(par_hat)) &&
        all(par_hat > 0) &&
        res$convergence == 0,
      AIC = -2 * logLik + 4,
      BIC = -2 * logLik + log(length(x)) * 2
    )
  }
  
  all_fits <- list()
  id <- 1L
  
  for (st in start_list) {
    for (meth in methods) {
      all_fits[[id]] <- fit_one(st, meth)
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
    all_fits = all_fits,
    x = x,
    call = match.call()
  )
  
  class(out) <- c("fitCSB_mle", "list")
  return(out)
}