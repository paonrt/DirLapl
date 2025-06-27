################################################################################
quad_form = function(M, x) drop(crossprod(M %*% x, x))
################################################################################


################################################################################
#' @export
nmp_rng = function(n, qn, A) {
  theta0 = double(n)
  theta0[1:qn] = A
  return(rnorm(n, mean = theta0, sd = 1))
}
################################################################################


################################################################################
rdisj_std_unif = function(n, alpha, beta) {
  u = runif(n = n, min = 0, max = alpha + 1 - beta)
  u = u + (beta - alpha) * (u > alpha)
  return(u)
}
################################################################################


################################################################################
rng_rej = function(mean, sd) {
  tbs = rep(x = TRUE, times = length(mean))
  res = double(length = length(mean))
  while (any(tbs)) {
    res[tbs] = rnorm(sum(tbs), mean = mean, sd = sd)
    tbs[tbs] = abs(res[tbs]) < .Machine$double.eps
  }
  return(res)
}
################################################################################


################################################################################
rng_trunc = function(mean, sd) {
  p = 1 / (1 + exp(
    pnorm(
      .Machine$double.eps, mean = mean, sd = sd, lower.tail = FALSE,
      log.p = TRUE
    ) - pnorm(
      -.Machine$double.eps, mean = mean, sd = sd, log.p = TRUE
    )
  ))
  tbs_neg = as.logical(rbinom(n = length(mean), size = 1, prob = p))
  res = truncnorm::rtruncnorm(
    n = length(mean),
    a = (.Machine$double.eps / !tbs_neg) * sign((2 * !tbs_neg) - 1),
    b = (-.Machine$double.eps / tbs_neg) * sign((2 * tbs_neg) - 1),
    mean = mean, sd = sd
  )
  return(res)
}
################################################################################


################################################################################
logBesselK = function(x, nu) log(besselK(x, nu, expon.scaled = TRUE)) - x
################################################################################


################################################################################
log_kernel_alpha = function(abs_theta, seq_alpha, p) {
  colSums(logBesselK(
    x = sqrt(2 * abs_theta),
    nu = 1 - matrix(seq_alpha, ncol = length(seq_alpha), nrow = p, byrow = TRUE)
  )) + colSums(
    0.5 * (matrix(
      seq_alpha, ncol = length(seq_alpha), nrow = p, byrow = TRUE
    ) - 1) * log(abs_theta)
  )
}
################################################################################