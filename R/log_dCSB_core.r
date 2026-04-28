log_dCSB_core <- function(x, shape, q) {
  if (length(shape) != 1L || length(q) != 1L) {
    stop("shape and q must be scalars.", call. = FALSE)
  }
  if (!is.numeric(shape) || !is.numeric(q) ||
      !is.finite(shape) || !is.finite(q) ||
      shape <= 0 || q <= 0) {
    stop("shape and q must be positive finite scalars.", call. = FALSE)
  }
  if (!is.numeric(x)) {
    stop("x must be numeric.", call. = FALSE)
  }
  
  alpha <- shape
  beta  <- shape
  
  out <- rep(-Inf, length(x))
  inside <- is.finite(x) & (x > 0) & (x < 1)
  
  if (!any(inside)) {
    return(out)
  }
  
  xx <- x[inside]
  
  log_const <- log(q) +
    lbeta(shape + q, shape) -
    lbeta(shape, shape) -
    log(expm1(q * log(2)))
  
  log_kernel <- -(q + 1) * log(xx)
  
  log_prob <- rep(-Inf, length(xx))
  
  idx1 <- xx <= 0.5
  idx2 <- !idx1
  
  if (any(idx1)) {
    z <- xx[idx1]
    z2 <- pmin(2 * z, 1)
    
    logF_2z <- stats::pbeta(
      z2,
      shape1 = shape + q,
      shape2 = shape,
      log.p = TRUE
    )
    
    logF_z <- stats::pbeta(
      z,
      shape1 = shape + q,
      shape2 = shape,
      log.p = TRUE
    )
    
    delta <- logF_z - logF_2z
    delta[is.finite(delta) & delta > 0 & delta < 1e-12] <- 0
    
    lp <- rep(-Inf, length(delta))
    ok <- is.finite(logF_2z) & is.finite(delta) & (delta <= 0)
    
    if (any(ok)) {
      lp[ok] <- logF_2z[ok] + log1mexp(delta[ok])
    }
    
    log_prob[idx1] <- lp
  }
  
  if (any(idx2)) {
    z <- xx[idx2]
    
    log_prob[idx2] <- stats::pbeta(
      z,
      shape1 = shape + q,
      shape2 = shape,
      lower.tail = FALSE,
      log.p = TRUE
    )
  }
  
  ans <- log_const + log_kernel + log_prob
  ans[!is.finite(ans)] <- -Inf
  
  out[inside] <- ans
  out
}