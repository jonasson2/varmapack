test_that("named and indexed testcases create models", {
  named <- varmapack_testcase("smallARMA1")
  indexed <- varmapack_testcase(8)
  expect_equal(named$p, indexed$p)
  expect_equal(named$q, indexed$q)
  expect_equal(named$r, indexed$r)
})

test_that("unnamed testcases support rho and Randompack RNGs", {
  rho_case <- varmapack_testcase("rho", p = 2, q = 1, r = 2, rho = 0.8)
  expect_equal(rho_case$specrad(), 0.8, tolerance = 2e-6)
  rng1 <- randompack::randompack_rng()
  rng2 <- randompack::randompack_rng()
  rng1$seed(321)
  rng2$seed(321)
  random1 <- varmapack_testcase("random", p = 2, q = 1, r = 2, rng = rng1)
  random2 <- varmapack_testcase("random", p = 2, q = 1, r = 2, rng = rng2)
  expect_identical(random1$A, random2$A)
  expect_identical(random1$B, random2$B)
})

test_that("testcase listing has all named cases", {
  cases <- varmapack_testcases()
  expect_equal(nrow(cases), 16L)
  expect_named(cases, c("index", "name", "p", "q", "r"))
  expect_true("tinyAR" %in% cases$name)
  expect_true("largeARMA" %in% cases$name)
})
