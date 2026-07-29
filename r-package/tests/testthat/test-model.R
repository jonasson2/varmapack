test_that("VARMA simulation accepts Randompack RNG objects", {
  model <- varmapack_model(A = matrix(0.4, 1, 1), B = matrix(0.2, 1, 1),
                            Sig = matrix(1, 1, 1))
  rng1 <- randompack::randompack_rng()
  rng2 <- randompack::randompack_rng()
  rng1$seed(123)
  rng2$seed(123)
  X1 <- model$sim(20, nrep = 3, rng = rng1)
  X2 <- model$sim(20, nrep = 3, rng = rng2)
  expect_identical(X1, X2)
  expect_equal(dim(X1), c(1L, 20L, 3L))
})

test_that("simulation can return shock arrays", {
  model <- varmapack_model(Sig = diag(2))
  out <- model$sim(10, nrep = 2, return_shocks = TRUE)
  expect_named(out, c("X", "E"))
  expect_equal(dim(out$X), c(2L, 10L, 2L))
  expect_equal(dim(out$E), c(2L, 10L, 2L))
})

test_that("time-dependent VARMA means are accepted", {
  model <- varmapack_model(A = matrix(0.2, 1, 1), Sig = matrix(1, 1, 1),
                            mu = matrix(c(1, 2), 1, 2))
  expect_equal(dim(model$sim(12, nrep = 2)), c(1L, 12L, 2L))
})

test_that("model input validation rejects incompatible inputs", {
  expect_error(varmapack_model(A = matrix(1, 2, 2), Sig = diag(1)))
  expect_error(varmapack_model(Sig = matrix(c(1, 2, 0, 1), 2, 2)))
  model <- varmapack_model(Sig = matrix(1, 1, 1))
  expect_error(model$sim(10, z = matrix(0, 1, 10)))
})
