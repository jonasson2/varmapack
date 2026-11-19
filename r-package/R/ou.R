#' Create a stationary Ornstein--Uhlenbeck model
#'
#' A Gaussian continuous-time model, simulated exactly at arbitrary
#' observation times without burn-in. With diagonal drift, component i has
#' autocorrelation `exp(-gap/tau[i])`.
#'
#' @param tau Finite positive numeric vector of persistence times. Supply
#'   exactly one of `tau` and `drift`.
#' @param cov Stationary process covariance, not innovation covariance. Only
#'   the lower triangle is used. Both `cov` and the driving covariance
#'   `K %*% cov + cov %*% t(K)` must be positive semidefinite, where
#'   `K = diag(1/tau)` for persistence-time initialization. Singular
#'   covariances are allowed.
#' @param mu Optional constant mean vector, defaulting to zero.
#' @param drift Square matrix K in the drift term `-K %*% (X - mu)`.
#'   Its eigenvalues must have strictly positive real parts.
#'
#' @details Over a gap of length d, the autoregressive coefficient matrix is
#'   the matrix exponential \eqn{\exp(-Kd)}. With constant gaps, the sampled
#'   process is a VAR(1); with unequal gaps, the matrix varies from step to step.
#'
#' @return A `VarmapackOUModel` object with fields `tau`, `drift`, `cov`, `mu`, `r`,
#'   and method `sim(t, nrep = 1L, rng = NULL)`. The result has dimensions
#'   `r` by `length(t)` by `nrep`, including for one-dimensional models.
#'   The unused parameter field (`tau` or `drift`) is `NULL`.
#'   Times must be finite and strictly increasing.
#'   The first observation is drawn from `N(mu, cov)`. Supply a
#'   `randompack::randompack_rng()` through `rng` to control the stream.
#'   `acvf(lags)` returns Cov(X(t+h), X(t)) as an r by r by length(lags) array.
#'   `specrad(gap)` returns the sampled VAR(1) spectral radius at a fixed gap.
#'
#' @examples
#' rng <- randompack::randompack_rng()
#' rng$seed(123)
#' t <- varmapack_gamma_times(100, shape = 16, rng = rng)
#' model <- varmapack_ou_model(c(1, 2), matrix(c(1, .5, .5, 1), 2))
#' X <- model$sim(t, nrep = 3, rng = rng)
#' @export
varmapack_ou_model <- function(tau = NULL, cov, mu = NULL, drift = NULL) {
  VarmapackOUModel$new(tau, cov, mu, drift)
}

.vp_rng_pointer <- function(rng) {
  if (is.null(rng)) return(NULL)
  if (!inherits(rng, "RandompackRNG") || is.null(rng$ptr))
    stop("rng must be a randompack_rng() object", call. = FALSE)
  rng$ptr
}

.vp_lags <- function(lags) {
  lags <- .vp_numeric_array(lags, "lags")
  if (!is.null(dim(lags))) stop("lags must be a vector", call. = FALSE)
  as.double(lags)
}

VarmapackOUModel <- R6::R6Class(
  "VarmapackOUModel",
  public = list(
    tau = NULL, drift = NULL, cov = NULL, mu = NULL, r = NULL,
    initialize = function(tau, cov, mu = NULL, drift = NULL) {
      if (is.null(tau) == is.null(drift))
        stop("specify exactly one of tau and drift", call. = FALSE)
      if (is.null(drift)) {
        tau <- .vp_numeric_array(tau, "tau")
        if (!is.null(dim(tau)) || length(tau) < 1L || any(tau <= 0))
          stop("tau must be a positive vector", call. = FALSE)
        tau <- as.double(tau)
        r <- length(tau)
      } else {
        drift <- .vp_numeric_array(drift, "drift")
        if (!is.matrix(drift) || nrow(drift) < 1L || nrow(drift) != ncol(drift))
          stop("drift must be a nonempty square matrix", call. = FALSE)
        r <- nrow(drift)
        drift <- matrix(as.double(drift), r, r)
      }
      if (!is.numeric(cov) || !is.matrix(cov) ||
          nrow(cov) != r || ncol(cov) != r)
        stop("cov must be an r by r matrix", call. = FALSE)
      cov <- matrix(as.double(cov), r, r)
      cov[upper.tri(cov)] <- t(cov)[upper.tri(cov)]
      if (is.null(mu)) mu <- numeric(r)
      mu <- .vp_numeric_array(mu, "mu")
      if (!is.null(dim(mu)) || length(mu) != r)
        stop("mu must be a vector of length r", call. = FALSE)
      parameter <- if (is.null(drift)) tau else drift
      .Call("varmapack_ou_validate_R", parameter, cov, as.double(mu),
            PACKAGE = "varmapack")
      self$tau <- tau
      self$drift <- drift
      self$cov <- cov
      self$mu <- as.double(mu)
      self$r <- r
    },
    specrad = function(gap) {
      gap <- .vp_numeric_array(gap, "gap")
      if (length(gap) != 1L || !is.null(dim(gap)))
        stop("gap must be a real scalar", call. = FALSE)
      parameter <- if (is.null(self$drift)) self$tau else self$drift
      .Call("varmapack_ou_specrad_R", parameter, self$cov, self$mu, as.double(gap),
            PACKAGE = "varmapack")
    },
    acvf = function(lags) {
      parameter <- if (is.null(self$drift)) self$tau else self$drift
      .Call("varmapack_ou_acvf_R", parameter, self$cov, self$mu, .vp_lags(lags),
            PACKAGE = "varmapack")
    },
    sim = function(t, nrep = 1L, rng = NULL) {
      t <- .vp_numeric_array(t, "t")
      if (!is.null(dim(t)) || length(t) < 1L || any(diff(t) <= 0))
        stop("t must be a nonempty strictly increasing vector", call. = FALSE)
      nrep <- .vp_count(nrep, "nrep", 1L)
      parameter <- if (is.null(self$drift)) self$tau else self$drift
      .Call("varmapack_ou_sim_R", parameter, self$cov, self$mu, as.double(t),
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
