#' Create a stationary Slepian-based moving-average model
#'
#' A stationary Gaussian model with unit-width Brownian noise windows.
#' It reduces to a vector MA(1) at unit spacing and can be simulated
#' exactly at arbitrary observation times without burn-in.
#'
#' @param theta Real square unit-step MA(1) coefficient matrix.
#' @param cov Stationary process covariance; only the lower triangle is used.
#'   The driving covariance `Sig` must be positive semidefinite and
#'   determined by `cov = Sig + theta %*% Sig %*% t(theta)`.
#' @param mu Optional constant mean vector, defaulting to zero.
#'
#' @return A `VarmapackSlepianModel` with fields `theta`, `cov`, `mu`, `r`,
#'   and method `sim(t, nrep = 1L, rng = NULL)`. The result has dimensions
#'   `r` by `length(t)` by `nrep`. Times must be finite and strictly
#'   increasing. Supply a `randompack::randompack_rng()` through `rng`
#'   to control the stream.
#'
#' @examples
#' rng <- randompack::randompack_rng()
#' rng$seed(123)
#' t <- varmapack_gamma_times(100, shape = 16, rng = rng)
#' model <- varmapack_slepian_model(matrix(.5), matrix(1.25))
#' X <- model$sim(t, nrep = 3, rng = rng)
#' @export
varmapack_slepian_model <- function(theta, cov, mu = NULL) {
  VarmapackSlepianModel$new(theta, cov, mu)
}

VarmapackSlepianModel <- R6::R6Class(
  "VarmapackSlepianModel",
  private = list(params = NULL),
  active = list(
    theta = function() private$params$theta,
    cov = function() private$params$cov,
    mu = function() private$params$mu,
    r = function() nrow(private$params$theta)
  ),
  public = list(
    initialize = function(theta, cov, mu = NULL) {
      theta <- .vp_numeric_array(theta, "theta")
      if (!is.matrix(theta) || nrow(theta) < 1L || nrow(theta) != ncol(theta))
        stop("theta must be a nonempty square matrix", call. = FALSE)
      r <- nrow(theta)
      if (!is.numeric(cov) || !is.matrix(cov) ||
          nrow(cov) != r || ncol(cov) != r)
        stop("cov must be an r by r matrix", call. = FALSE)
      theta <- matrix(as.double(theta), r, r)
      cov <- matrix(as.double(cov), r, r)
      cov[upper.tri(cov)] <- t(cov)[upper.tri(cov)]
      if (is.null(mu)) mu <- numeric(r)
      mu <- .vp_numeric_array(mu, "mu")
      if (!is.null(dim(mu)) || length(mu) != r)
        stop("mu must be a vector of length r", call. = FALSE)
      mu <- as.double(mu)
      factor <- .Call("varmapack_slepian_factor_R", theta, cov, mu, PACKAGE = "varmapack")
      private$params <- list(theta = theta, cov = cov, mu = mu, factor = factor)
    },
    sim = function(t, nrep = 1L, rng = NULL) {
      t <- .vp_numeric_array(t, "t")
      if (!is.null(dim(t)) || length(t) < 1L || any(t[-1L] <= t[-length(t)]))
        stop("t must be a nonempty strictly increasing vector", call. = FALSE)
      nrep <- .vp_count(nrep, "nrep", 1L)
      .Call("varmapack_slepian_sim_R", self$theta, private$params$factor, self$mu,
            as.double(t), nrep, .vp_rng_pointer(rng), PACKAGE = "varmapack")
    }
  )
)
