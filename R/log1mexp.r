log1mexp <- function(a) {
  out <- rep(NaN, length(a))

  ok <- is.finite(a) & (a <= 0)
  if (any(ok)) {
    a_ok <- a[ok]
    idx <- a_ok > -log(2)
    tmp <- numeric(length(a_ok))
    tmp[idx]  <- log(-expm1(a_ok[idx]))
    tmp[!idx] <- log1p(-exp(a_ok[!idx]))
    out[ok] <- tmp
  }

  out[is.finite(a) & abs(a) < 1e-15] <- -Inf
  out
}
