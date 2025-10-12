# helper_testcase.R
# Expected catalogue of named cases for `testcase()`
.named_cases_tbl <- data.frame(
  no   = 1:12,
  name = c(
    "tinyAR","tinyMA","tinyARMA",
    "smallAR","smallMA","smallARMA1","smallARMA2",
    "mediumAR","mediumARMA1","mediumARMA2","mediumMA",
    "largeAR"
  ),
  p = as.integer(c(1,0,1, 1,0,1,1, 1,3,3,0,5)),
  q = as.integer(c(0,1,1, 0,2,1,2, 0,3,3,2,0)),
  r = as.integer(c(1,1,1, 2,2,2,2, 3,3,3,3,7)),
  stringsAsFactors = FALSE
)

.check_one_case <- function(tc, exp_row) {
  stopifnot(is.list(tc), all(c("A","B","Sig","p","q","r","name") %in% names(tc)))
  # p/q/r/name
  testthat::expect_identical(tc$p, as.integer(exp_row$p))
  testthat::expect_identical(tc$q, as.integer(exp_row$q))
  testthat::expect_identical(tc$r, as.integer(exp_row$r))
  testthat::expect_identical(tc$name, exp_row$name)

  # shapes
  testthat::expect_true(is.array(tc$A))
  testthat::expect_true(is.array(tc$B))
  testthat::expect_true(is.matrix(tc$Sig))

  testthat::expect_identical(dim(tc$A), c(exp_row$r, exp_row$r, exp_row$p))
  testthat::expect_identical(dim(tc$B), c(exp_row$r, exp_row$r, exp_row$q))
  testthat::expect_identical(dim(tc$Sig), c(exp_row$r, exp_row$r))
}

# Public helper used by tests
check_testcase_named <- function() {
  # (a) Individual named calls
  for (k in seq_len(nrow(.named_cases_tbl))) {
    exp_row <- .named_cases_tbl[k, ]
    tc <- testcase(name = exp_row$name)
    .check_one_case(tc, exp_row)
  }

  # (b) `all = TRUE` returns all 12, correctly named and shaped
  all_tc <- testcase(all = TRUE)
  testthat::expect_true(is.list(all_tc))
  testthat::expect_identical(length(all_tc), nrow(.named_cases_tbl))
  testthat::expect_identical(names(all_tc), .named_cases_tbl$name)

  for (k in seq_along(all_tc)) {
    tc <- all_tc[[k]]
    exp_row <- .named_cases_tbl[k, ]
    .check_one_case(tc, exp_row)
  }
}
check_testcase_random <- function(verbose = FALSE) {
  combos <- expand.grid(p = 0:2, q = 0:2, r = 1:3)

  for (k in seq_len(nrow(combos))) {
    p <- combos$p[k]; q <- combos$q[k]; r <- combos$r[k]

    tc <- testcase(p = p, q = q, r = r)

    # metadata
    expect_identical(tc$p, as.integer(p))
    expect_identical(tc$q, as.integer(q))
    expect_identical(tc$r, as.integer(r))

    # shapes
    expect_true(is.array(tc$A))
    expect_true(is.array(tc$B))
    expect_true(is.matrix(tc$Sig))

    expect_identical(dim(tc$A), c(r, r, p))
    expect_identical(dim(tc$B), c(r, r, q))
    expect_identical(dim(tc$Sig), c(r, r))

    # name
    expected_name <- sprintf("random(p=%d,q=%d,r=%d)", p, q, r)
    expect_identical(tc$name, expected_name)

    if (verbose) {
      cat(sprintf("✓ random(%d,%d,%d) OK\n", p, q, r))
    }
  }
}

check_testcase_invalid <- function() {
  # all=TRUE should not be combined with name
  expect_error(testcase(name = "smallAR", all = TRUE))

  # all=TRUE should not be combined with p,q,r
  expect_error(testcase(p = 1, q = 1, r = 2, all = TRUE))

  # invalid name should error
  expect_error(testcase(name = "not_a_case"))

  # negative dimensions should error
  expect_error(testcase(p = -1, q = 1, r = 2))
  expect_error(testcase(p = 1, q = -1, r = 2))
  expect_error(testcase(p = 1, q = 1, r = -3))

  # non-integer dimensions should error
  expect_error(testcase(p = 1.5, q = 1, r = 2))
  expect_error(testcase(p = 1, q = 2.7, r = 2))
  expect_error(testcase(p = 1, q = 1, r = 2.2))
}
