test_that("Slepian paths have the specified lagged covariances", {
  theta <- matrix(c(0, .25, .5, 0), 2)
  cov <- diag(c(1.25, 1.0625))
  mu <- c(2, -1)
  model <- varmapack_slepian_model(theta, cov, mu)
  rng <- randompack::randompack_rng()
  rng$seed(423)
  X <- model$sim(c(0, .5, 1, 2), nrep = 20000L, rng = rng)
  expect_equal(dim(X), c(2L, 4L, 20000L))
  expect_equal(model$theta, theta)
  expect_equal(model$cov, cov)
  expect_error(model$theta <- diag(2))
  restored <- unserialize(serialize(model, NULL))
  rng$seed(423)
  expect_identical(X, restored$sim(c(0, .5, 1, 2), nrep = 20000L, rng = rng))
  expected <- list(cov, (cov + theta) / 2, theta, matrix(0, 2, 2))
  for (k in seq_along(expected)) {
    empirical <- (X[, k, ] - mu) %*% t(X[, 1, ] - mu) / dim(X)[3]
    expect_lt(max(abs(empirical - expected[[k]])), .035)
  }
  rng$seed(423)
  expect_identical(X, model$sim(c(0, .5, 1, 2), nrep = 20000L, rng = rng))
})

test_that("Slepian rejects invalid models and observation times", {
  expect_error(varmapack_slepian_model(c(0, .5), diag(2)), "square")
  expect_error(varmapack_slepian_model(matrix(.5), diag(2)), "cov")
  expect_error(varmapack_slepian_model(matrix(c(0, 0, 1, 0), 2), diag(c(1, 2))))
  model <- varmapack_slepian_model(matrix(0), matrix(0), mu = 3)
  expect_equal(model$sim(0), array(3, c(1L, 1L, 1L)))
  expect_error(model$sim(c(0, 0)), "increasing")
  expect_error(model$sim(c(0, Inf)), "finite")
  expect_error(model$sim(0, nrep = 0), "nrep")
})
