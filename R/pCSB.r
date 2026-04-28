#' Distribution Function of the Contracted Symmetric Beta Distribution
#'
#' Computes the cumulative distribution function of the Contracted Symmetric
#' Beta (CSB) distribution.
#'
#' The CSB distribution is defined through the stochastic representation
#' \deqn{T = X/U^{1/q},}
#' where \eqn{X \sim Beta(\alpha,\alpha)}, \eqn{U \sim Uniform(1,2^q)},
#' and \eqn{X} and \eqn{U} are independent.
#'
#' @param x Numeric vector of quantiles.
#' @param shape Positive shape parameter \eqn{\alpha}.
#' @param q Positive contraction parameter.
#' @param lower.tail Logical. If \code{TRUE}, probabilities are \eqn{P(T \le x)};
#' otherwise, \eqn{P(T > x)}.
#' @param log.p Logical. If \code{TRUE}, probabilities are returned on the log scale.
#'
#' @return A numeric vector of probabilities.
#'
#' @examples
#' pCSB(0.5, shape = 2, q = 5)
#' pCSB(c(0.25, 0.5, 0.75), shape = 2, q = 5)
#' pCSB(0.5, shape = 2, q = 5, lower.tail = FALSE)
#'
#' @export
pCSB <- function(x, shape, q, lower.tail = TRUE, log.p = FALSE) {
  
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
  
  out <- rep(NA_real_, length(x))
  
  out[x <= 0] <- 0
  out[x >= 1] <- 1
  
  inside <- is.finite(x) & x > 0 & x < 1
  
  if (any(inside)) {
    
    xx <- x[inside]
    
    den <- expm1(q * log(2))       # 2^q - 1
    twoq <- exp(q * log(2))        # 2^q
    rho <- exp(lbeta(shape + q, shape) - lbeta(shape, shape))
    
    F <- function(z) {
      stats::pbeta(z, shape1 = shape, shape2 = shape)
    }
    
    G <- function(z) {
      stats::pbeta(z, shape1 = shape + q, shape2 = shape)
    }
    
    val <- numeric(length(xx))
    
    idx1 <- xx <= 0.5
    idx2 <- !idx1
    
    # Case 1: 0 < x <= 1/2
    if (any(idx1)) {
      z <- xx[idx1]
      z2 <- 2 * z
      
      val[idx1] <- (
        (z2^q) * F(z2) -
          (z^q) * F(z) -
          rho * (G(z2) - G(z))
      ) / (den * z^q)
    }
    
    # Case 2: 1/2 < x < 1
    if (any(idx2)) {
      z <- xx[idx2]
      
      val[idx2] <- (
        twoq -
          F(z) -
          rho * z^(-q) * (1 - G(z))
      ) / den
    }
    
    val <- pmin(pmax(val, 0), 1)
    out[inside] <- val
  }
  
  if (!lower.tail) {
    out <- 1 - out
  }
  
  if (log.p) {
    out <- log(out)
  }
  
  out
}