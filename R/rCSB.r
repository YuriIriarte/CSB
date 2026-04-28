#' Random Generation from the Contracted Symmetric Beta Distribution
#'
#' Generates random observations from the Contracted Symmetric Beta (CSB)
#' distribution.
#'
#' The CSB distribution is defined through the stochastic representation
#' \deqn{T = X/U^{1/q},}
#' where \eqn{X \sim Beta(\alpha,\alpha)}, \eqn{U \sim Uniform(1,2^q)},
#' and \eqn{X} and \eqn{U} are independent.
#'
#' @param n Number of observations.
#' @param shape Positive shape parameter \eqn{\alpha}.
#' @param q Positive contraction parameter.
#'
#' @return A numeric vector of random observations in \eqn{(0,1)}.
#'
#' @examples
#' set.seed(123)
#' x <- rCSB(100, shape = 2, q = 5)
#' hist(x, probability = TRUE)
#' curve(dCSB(x, shape = 2, q = 5), add = TRUE)
#'
#' @export
rCSB <- function(n, shape, q) {
  
  if (length(n) != 1L || !is.numeric(n) || !is.finite(n) || n < 0) {
    stop("n must be a non-negative finite scalar.", call. = FALSE)
  }
  if (length(shape) != 1L || length(q) != 1L) {
    stop("shape and q must be scalars.", call. = FALSE)
  }
  if (!is.numeric(shape) || !is.numeric(q) ||
      !is.finite(shape) || !is.finite(q) ||
      shape <= 0 || q <= 0) {
    stop("shape and q must be positive finite scalars.", call. = FALSE)
  }
  
  n <- as.integer(n)
  
  if (n == 0L) {
    return(numeric(0))
  }
  
  X <- stats::rbeta(n, shape1 = shape, shape2 = shape)
  U <- stats::runif(n, min = 1, max = exp(q * log(2)))
  
  T <- X / (U^(1 / q))
  
  T
}