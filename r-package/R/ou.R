#' Create a stationary Ornstein--Uhlenbeck model
#'
#' A Gaussian continuous-time model with diagonal drift, simulated exactly at
#' arbitrary observation times. Component i has autocorrelation
#' `exp(-gap/tau[i])`. Correlated driving noise couples the series; there are
#' no cross-variable drift terms. No burn-in is required.
#'
#' @param tau Finite positive numeric vector of persistence times.
#' @param cov Stationary process covariance, not innovation covariance. Only
#'   the lower triangle is used. Both `cov` and the driving covariance with
#'   entries `cov[i,j]*(1/tau[i] + 1/tau[j])` must be positive semidefinite.
#'   Singular covariances are allowed.
#' @param mu Optional constant mean vector, defaulting to zero.
#'
#' @return A `VarmapackOUModel` object with fields `tau`, `cov`, `mu`, `r`,
#'   and method `sim(t, nrep = 1L, rng = NULL)`. The result has dimensions
#'   `r` by `length(t)` by `nrep`, including for one-dimensional models.
#'   Times must be finite and strictly increasing, in the units of `tau`.
#'   The first observation is drawn from `N(mu, cov)`. Supply a
#'   `randompack::randompack_rng()` through `rng` to control the stream.
#'
#' @examples
#' rng <- randompack::randompack_rng()
#' rng$seed(123)
#' t <- varmapack_gamma_times(100, shape = 16, rng = rng)
#' model <- varmapack_ou_model(c(1, 2), matrix(c(1, .5, .5, 1), 2))
#' X <- model$sim(t, nrep = 3, rng = rng)
#' @export
varmapack_ou_model <- function(tau, cov, mu = NULL) {
  VarmapackOUModel$new(tau, cov, mu)
}

.vp_rng_pointer <- function(rng) {
  if (is.null(rng)) return(NULL)
  if (!inherits(rng, "RandompackRNG") || is.null(rng$ptr))
    stop("rng must be a randompack_rng() object", call. = FALSE)
  rng$ptr
}

VarmapackOUModel <- R6::R6Class(
  "VarmapackOUModel",
  public = list(
    tau = NULL, cov = NULL, mu = NULL, r = NULL,
    initialize = function(tau, cov, mu = NULL) {
      tau <- .vp_numeric_array(tau, "tau")
      if (!is.null(dim(tau)) || length(tau) < 1L || any(tau <= 0))
        stop("tau must be a positive vector", call. = FALSE)
      r <- length(tau)
      if (!is.numeric(cov) || !is.matrix(cov) ||
          nrow(cov) != r || ncol(cov) != r)
        stop("cov must be an r by r matrix", call. = FALSE)
      cov <- matrix(as.double(cov), r, r)
      cov[upper.tri(cov)] <- t(cov)[upper.tri(cov)]
      if (is.null(mu)) mu <- numeric(r)
      mu <- .vp_numeric_array(mu, "mu")
      if (!is.null(dim(mu)) || length(mu) != r)
        stop("mu must be a vector of length r", call. = FALSE)
      .Call("varmapack_ou_validate_R", as.double(tau), cov, as.double(mu),
            PACKAGE = "varmapack")
      self$tau <- as.double(tau)
      self$cov <- cov
      self$mu <- as.double(mu)
      self$r <- r
    },
    sim = function(t, nrep = 1L, rng = NULL) {
      t <- .vp_numeric_array(t, "t")
      if (!is.null(dim(t)) || length(t) < 1L || any(diff(t) <= 0))
        stop("t must be a nonempty strictly increasing vector", call. = FALSE)
      nrep <- .vp_count(nrep, "nrep", 1L)
      .Call("varmapack_ou_sim_R", self$tau, self$cov, self$mu, as.double(t),
            nrep, .vp_rng_pointer(rng), PACKAGE = "varmapack")
    }
  )
)

#' Generate observation times with Gamma spacing
#'
#' @param n Number of times, at least one.
#' @param shape Positive Gamma shape. `Inf` gives constant spacing.
#' @param mean_gap Positive expected gap, in the same units as OU persistence.
#' @param rng Optional `randompack::randompack_rng()` object.
#'
#' @return A strictly increasing numeric vector starting at zero. The `n-1`
#'   independent gaps have scale `mean_gap/shape` and coefficient of variation
#'   `1/sqrt(shape)`. Their realized mean is not normalized. Unrepresentable
#'   positive gaps or times produce an error.
#' @examples
#' t <- varmapack_gamma_times(100, shape = 4, mean_gap = 1)
#' @export
varmapack_gamma_times <- function(n, shape = 16, mean_gap = 1, rng = NULL) {
  n <- .vp_count(n, "n", 1L)
  if (!is.numeric(shape) || length(shape) != 1L || is.na(shape) || shape <= 0)
    stop("shape must be positive or Inf", call. = FALSE)
  if (!is.numeric(mean_gap) || length(mean_gap) != 1L ||
      !is.finite(mean_gap) || mean_gap <= 0)
    stop("mean_gap must be positive and finite", call. = FALSE)
  .Call("varmapack_gamma_times_R", n, as.double(shape), as.double(mean_gap),
        .vp_rng_pointer(rng), PACKAGE = "varmapack")
}
