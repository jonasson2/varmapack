test_that("continuous autocovariance preserves lag order and orientation", {
  ou <- varmapack_ou_model(cov = diag(2), drift = matrix(c(1, 0, -1, 1), 2))
  theta <- matrix(c(0, .25, .5, 0), 2)
  cov <- diag(2) + theta %*% t(theta)
  ma <- varmapack_slepian_model(theta, cov)
  for (pair in list(list(ou, exp(-.5) * matrix(c(1, 0, .5, 1), 2), diag(2)),
                    list(ma, (cov + theta) / 2, cov))) {
    model <- pair[[1]]
    G <- model$acvf(c(.5, -.5, 0))
    expect_equal(dim(G), c(2L, 2L, 3L))
    expect_equal(G[, , 1], pair[[2]])
    expect_equal(G[, , 2], t(pair[[2]]))
    expect_equal(G[, , 3], pair[[3]])
    expect_equal(dim(model$acvf(numeric())), c(2L, 2L, 0L))
    expect_error(model$acvf(Inf), "lags")
  }
  expect_equal(ma$acvf(2)[, , 1], matrix(0, 2, 2))
  expect_equal(ou$specrad(0.5), exp(-0.5))
  expect_error(ou$specrad(0), "gap")
})
