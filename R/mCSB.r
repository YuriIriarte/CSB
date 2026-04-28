#' Raw Moments of the Contracted Symmetric Beta Distribution
#'
#' Computes the raw moment of order \eqn{r} for the Contracted Symmetric Beta
#' (CSB) distribution.
#'
#' The r-th raw moment is given by
#'
#' \deqn{
#' \mathbb{E}(T^r)=
#' \prod_{j=0}^{r-1}\frac{\alpha+j}{2\alpha+j}
#' \times
#' \left\{
#' \begin{array}{ll}
#' \frac{q(2^{q-r}-1)}{(2^q-1)(q-r)}, & q\neq r,\\
#' \frac{q\log(2)}{2^q-1}, & q=r.
#' \end{array}
#' \right.
#' }
#'
#' @param r Non-negative order or vector of orders of the moment.
#' @param shape Positive shape parameter \eqn{\alpha}.
#' @param q Positive contraction parameter.
#'
#' @return A numeric vector with the raw moments of order \eqn{r}.
#'
#' @examples
#' mCSB(1, shape = 2, q = 5)
#' mCSB(1:4, shape = 2, q = 5)
#'
#' @export
mCSB <- function(r, shape, q) {
  
  if (!is.numeric(r) || any(!is.finite(r)) || any(r < 0)) {
    stop("r must be non-negative and finite.", call. = FALSE)
  }
  if (length(shape) != 1L || length(q) != 1L) {
    stop("shape and q must be scalars.", call. = FALSE)
  }
  if (!is.numeric(shape) || !is.numeric(q) ||
      !is.finite(shape) || !is.finite(q) ||
      shape <= 0 || q <= 0) {
    stop("shape and q must be positive finite scalars.", call. = FALSE)
  }
  
  out <- numeric(length(r))
  
  for (i in seq_along(r)) {
    
    rr <- r[i]
    
    beta_ratio <- exp(
      lbeta(shape + rr, shape) - lbeta(shape, shape)
    )
    
    if (abs(rr - q) < 1e-10) {
      c_rq <- q * log(2) / expm1(q * log(2))
    } else {
      c_rq <- q * expm1((q - rr) * log(2)) /
        (expm1(q * log(2)) * (q - rr))
    }
    
    out[i] <- c_rq * beta_ratio
  }
  
  out
}


#' Mean of the Contracted Symmetric Beta Distribution
#'
#' Computes the theoretical mean of the CSB distribution.
#'
#' @inheritParams mCSB
#'
#' @return A numeric value.
#'
#' @examples
#' meanCSB(shape = 2, q = 5)
#'
#' @export
meanCSB <- function(shape, q) {
  mCSB(1, shape = shape, q = q)
}


#' Variance of the Contracted Symmetric Beta Distribution
#'
#' Computes the theoretical variance of the CSB distribution.
#'
#' @inheritParams mCSB
#'
#' @return A numeric value.
#'
#' @examples
#' varCSB(shape = 2, q = 5)
#'
#' @export
varCSB <- function(shape, q) {
  m1 <- mCSB(1, shape = shape, q = q)
  m2 <- mCSB(2, shape = shape, q = q)
  m2 - m1^2
}


#' Coefficient of Variation of the Contracted Symmetric Beta Distribution
#'
#' Computes the theoretical coefficient of variation of the CSB distribution.
#'
#' @inheritParams mCSB
#'
#' @return A numeric value.
#'
#' @examples
#' cvCSB(shape = 2, q = 5)
#'
#' @export
cvCSB <- function(shape, q) {
  mu <- meanCSB(shape = shape, q = q)
  sig2 <- varCSB(shape = shape, q = q)
  sqrt(sig2) / mu
}


#' Skewness of the Contracted Symmetric Beta Distribution
#'
#' Computes the theoretical skewness coefficient of the CSB distribution.
#'
#' @inheritParams mCSB
#'
#' @return A numeric value.
#'
#' @examples
#' skewCSB(shape = 2, q = 5)
#'
#' @export
skewCSB <- function(shape, q) {
  m1 <- mCSB(1, shape = shape, q = q)
  m2 <- mCSB(2, shape = shape, q = q)
  m3 <- mCSB(3, shape = shape, q = q)
  
  mu3 <- m3 - 3 * m1 * m2 + 2 * m1^3
  sig2 <- m2 - m1^2
  
  mu3 / sig2^(3 / 2)
}


#' Kurtosis of the Contracted Symmetric Beta Distribution
#'
#' Computes the theoretical kurtosis or excess kurtosis of the CSB distribution.
#'
#' @inheritParams mCSB
#' @param excess Logical. If \code{TRUE}, the excess kurtosis is returned;
#' otherwise, the ordinary kurtosis is returned.
#'
#' @return A numeric value.
#'
#' @examples
#' kurtCSB(shape = 2, q = 5)
#' kurtCSB(shape = 2, q = 5, excess = FALSE)
#'
#' @export
kurtCSB <- function(shape, q, excess = TRUE) {
  m1 <- mCSB(1, shape = shape, q = q)
  m2 <- mCSB(2, shape = shape, q = q)
  m3 <- mCSB(3, shape = shape, q = q)
  m4 <- mCSB(4, shape = shape, q = q)
  
  mu4 <- m4 - 4 * m1 * m3 + 6 * m1^2 * m2 - 3 * m1^4
  sig2 <- m2 - m1^2
  
  kurt <- mu4 / sig2^2
  
  if (excess) {
    kurt <- kurt - 3
  }
  
  kurt
}