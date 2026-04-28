#' Density of the Contracted Symmetric Beta Distribution
#'
#' Computes the probability density function of the Contracted Symmetric Beta
#' (CSB) distribution.
#'
#' The CSB distribution is defined through the stochastic representation
#' \deqn{T = X/U^{1/q},}
#' where \eqn{X \sim Beta(\alpha,\alpha)}, \eqn{U \sim Uniform(1,2^q)},
#' and \eqn{X} and \eqn{U} are independent.
#'
#' @param x Numeric vector of quantiles.
#' @param shape Positive shape parameter \eqn{\alpha}.
#' @param q Positive contraction parameter.
#' @param log Logical. If \code{TRUE}, the log-density is returned.
#'
#' @return A numeric vector with the density or log-density values.
#'
#' @examples
#' dCSB(0.4, shape = 2, q = 5)
#' dCSB(c(0.2, 0.5, 0.8), shape = 2, q = 5)
#' dCSB(0.4, shape = 2, q = 5, log = TRUE)
#'
#' curve(dCSB(x, shape = 2, q = 5), from = 0, to = 1)
#'
#' @export
dCSB <- function(x, shape, q, log = FALSE) {
  lx <- log_dCSB_core(x = x, shape = shape, q = q)
  
  if (log) {
    return(lx)
  }
  
  exp(lx)
}