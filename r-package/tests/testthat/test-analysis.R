test_that("model analysis methods return expected AR(1) quantities", {
  model <- varmapack_model(A = matrix(0.5, 1, 1), Sig = matrix(1, 1, 1))
  expect_equal(model$specrad(), 0.5)
  expect_equal(model$ma_specrad(), 0)
  expect_equal(drop(model$acvf(1)), c(4/3, 2/3), tolerance = 1e-12)
  expect_equal(drop(model$psi(2)), c(1, 0.5, 0.25), tolerance = 1e-12)
  expect_equal(drop(model$irf(2)), c(1, 0.5, 0.25), tolerance = 1e-12)
})

test_that("sample autocovariance uses observations in rows", {
  X <- cbind(1:8, 2*(1:8))
  C <- varmapack_autocov(X, maxlag = 1)
  expect_equal(dim(C), c(2L, 2L, 2L))
  expect_equal(C[1, 2, 1], 2*C[1, 1, 1])
  expect_error(varmapack_autocov(X, maxlag = 8))
})

test_that("theoretical autocovariance rejects VARMAX models", {
  model <- varmapack_model(C = matrix(0.2, 1, 1), Sig = matrix(1, 1, 1))
  expect_error(model$acvf(1))
})
