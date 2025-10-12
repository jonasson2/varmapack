test_that("named testcases have correct shapes, metadata, and all=TRUE works", {
  check_testcase_named()
})
test_that("random testcases have correct shapes, metadata, and name", {
  check_testcase_random()
})
test_that("invalid testcase argument combinations error clearly", {
  check_testcase_invalid()
})
