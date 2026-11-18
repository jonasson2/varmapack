ou_rng <- function(seed = 423L) {
  rng <- randompack::randompack_rng()
  rng$seed(seed)
  rng
}

test_that("OU paths have stationary and lagged cross-covariances", {
  tau <- c(1, 2)
  cov <- matrix(c(1, .3, .3, 2), 2)
  mu <- c(2, -1)
  model <- varmapack_ou_model(tau, cov, mu)
  t <- c(-3, -2.9, -1, 3)
  X <- model$sim(t, nrep = 40000, rng = ou_rng())
  expect_equal(dim(X), c(2L, 4L, 40000L))
  expect_equal(model$tau, tau)
  expect_equal(model$cov, cov)
  for (k in seq_along(t)) {
    x <- X[, k, ] - mu
    expect_lt(max(abs(rowMeans(x))), .03)
    for (h in seq_len(k)) {
      empirical <- x %*% t(X[, h, ] - mu) / dim(X)[3]
      expected <- exp(-(t[k] - t[h]) / tau) * cov
      expect_lt(max(abs(empirical - expected)), .05)
    }
  }
})

test_that("OU RNG, PSD and input conventions work", {
  model <- varmapack_ou_model(c(1, 1), matrix(1, 2, 2))
  t <- c(0, 1e-15, 10000)
  X <- model$sim(t, nrep = 5, rng = ou_rng())
  expect_equal(X, model$sim(t, nrep = 5, rng = ou_rng()))
  expect_equal(X[1, , ], X[2, , ], tolerance = 1e-12)
  expect_equal(dim(model$sim(0)), c(2L, 1L, 1L))
  cov <- matrix(c(1, .2, NaN, 1), 2)
  model <- varmapack_ou_model(c(1, 2), cov)
  expect_equal(model$cov, matrix(c(1, .2, .2, 1), 2))
  model <- varmapack_ou_model(c(1, 2), matrix(0, 2, 2), c(3, 4))
  expect_equal(as.vector(model$sim(t)), rep(c(3, 4), 3))
})

test_that("OU rejects invalid models and observation times", {
  expect_error(varmapack_ou_model(0, matrix(1)), "tau")
  expect_error(varmapack_ou_model(Inf, matrix(1)), "finite")
  expect_error(varmapack_ou_model(1, diag(2)), "cov")
  expect_error(varmapack_ou_model(1, matrix(1), c(1, 2)), "mu")
  expect_error(varmapack_ou_model(c(1, 2), matrix(1, 2, 2)), "incompatible")
  expect_error(varmapack_ou_model(1, matrix(-1)), "cov")
  model <- varmapack_ou_model(1, matrix(1))
  for (t in list(numeric(), c(1, 1), c(1, 0), c(0, Inf), matrix(1)))
    expect_error(model$sim(t))
  expect_error(model$sim(0, nrep = 0))
  expect_error(model$sim(0, rng = list()))
})

test_that("OU accepts a full drift matrix", {
  K <- matrix(c(1L, 0L, -1L, 1L), 2)
  model <- varmapack_ou_model(cov = diag(2), drift = K)
  expect_null(model$tau)
  expect_equal(model$drift, K)
  X <- model$sim(c(0, .5), nrep = 20000L, rng = ou_rng())
  expect_lt(max(abs(X[, 2, ] %*% t(X[, 1, ]) / 20000 -
                    exp(-.5) * matrix(c(1, 0, .5, 1), 2))), .035)
  expect_error(varmapack_ou_model(c(1, 1), diag(2), drift = K), "exactly one")
  expect_error(varmapack_ou_model(cov = diag(2), drift = matrix(1)), "cov")
})

test_that("Gamma times have the specified gap moments and RNG behavior", {
  for (shape in c(.5, 1, 4, 16)) {
    t <- varmapack_gamma_times(50001, shape, 2, rng = ou_rng())
    d <- diff(t)
    expect_identical(t[1], 0)
    expect_true(all(d > 0))
    expect_lt(abs(mean(d) - 2), 7 * sqrt(4 / shape / length(d)))
    expect_lt(abs(var(d) - 4 / shape), .08 * 4 / shape)
    expect_equal(t, varmapack_gamma_times(length(t), shape, 2, rng = ou_rng()))
  }
  expect_equal(varmapack_gamma_times(5, Inf, 2), seq(0, 8, 2))
  expect_equal(varmapack_gamma_times(1), 0)
  expect_error(varmapack_gamma_times(0))
  expect_error(varmapack_gamma_times(3, 0))
  expect_error(varmapack_gamma_times(3, NA_real_))
  expect_error(varmapack_gamma_times(3, 4, 0))
  expect_error(varmapack_gamma_times(4, Inf, 1e308), "representable")
})
