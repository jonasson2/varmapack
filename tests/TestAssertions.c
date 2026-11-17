// Self-test the counters separately because some assertions deliberately fail.
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "TestUtil.h"

static int failures = 0;

static void expect_count(char *name, int failed) {
  if (xCheckNTotal() != 1 || xCheckNFailures() != failed) {
    fprintf(stderr, "%s: expected 1 check and %d failures, got %d and %d\n",
            name, failed, xCheckNTotal(), xCheckNFailures());
    failures++;
  }
  xCheckInit("intentional helper self-test");
}

static void set_error_state(bool set) {
  double cov = 1, corr;
  bool ok = varmapack_cov2corr(set ? 0 : &cov, 1, 0, &corr);
  if (ok == set || (varmapack_last_error() != 0) != set) {
    fprintf(stderr, "could not set the expected Varmapack error state\n");
    failures++;
  }
}

int main(void) {
  double x[] = {1, 2, 3, 4}, y[] = {1, 2, 5, 8};
  double zero[] = {0, -0.0, 0, 0}, bad[] = {0, NAN, INFINITY, -INFINITY};
  fprintf(stderr, "The assertion diagnostics below are deliberate and expected.\n");
  xCheckInit("intentional helper self-test");
  xCheck(true);
  expect_count("scalar success", 0);
  xCheck(false);
  expect_count("scalar failure", 1);
  checkArraySame(x, x, 4);
  expect_count("equal arrays", 0);
  checkArraySame(x, y, 4);
  expect_count("multiple array mismatches", 1);
  checkArrayTol(x, y, 4, 3);
  expect_count("last element mismatch", 1);
  checkArrayTol(x, y, 4, 5);
  expect_count("array within tolerance", 0);
  checkArraySame(bad, bad, 4);
  expect_count("NaN and infinity comparison", 1);
  checkArraySame(0, 0, 0);
  expect_count("empty comparison", 0);
  checkArrayZero(zero, 4);
  expect_count("zero array", 0);
  checkArrayZero(x, 4);
  expect_count("nonzero array", 1);
  checkArrayZero(bad, 4);
  expect_count("nonfinite entries are not zero", 1);
  checkArrayZero(0, 0);
  expect_count("empty zero array", 0);
  checkArrayFinite(x, 4);
  expect_count("finite array", 0);
  checkArrayFinite(bad, 4);
  expect_count("multiple nonfinite entries", 1);
  checkArrayFinite(0, 0);
  expect_count("empty finite array", 0);
  test_check_group empty = {0}, passing = {0}, failing = {0}, independent = {0};
  checkGroup(&empty, "empty property");
  expect_count("empty property group", 0);
  groupCheck(&passing, true);
  groupCheckMessage(&passing, true, "second element");
  checkGroup(&passing, "all elements satisfy property");
  expect_count("passing property group", 0);
  groupCheck(&failing, true);
  groupCheckMessage(&failing, false, "first bad element");
  groupCheck(&independent, true);
  groupCheckMessage(&failing, false, "second bad element");
  groupCheck(&failing, true);
  if (failing.total != 4 || failing.failed != 2 || failing.first_failed != 2 ||
      xCheckNTotal() != 0 || strcmp(failing.context, "first bad element") != 0) {
    fprintf(stderr, "property group lost counts or first failure context\n");
    failures++;
  }
  checkGroup(&failing, "all elements satisfy property");
  expect_count("multiple failures count once", 1);
  checkGroup(&independent, "independent property");
  expect_count("property groups do not share state", 0);
  test_check_group no_context = {0};
  groupCheck(&no_context, false);
  checkGroup(&no_context, "failure without extra context");
  expect_count("property failure without context", 1);
  randompack_rng *rng = seededRng(123);
  reseedRng(rng, 123);
  randompack_free(rng);
  if (xCheckNTotal() != 0 || xCheckNFailures() != 0) {
    fprintf(stderr, "RNG setup should not count as a Varmapack check\n");
    failures++;
  }
  set_error_state(false);
  checkVarmapackSuccess(true);
  expect_count("successful API outcome", 0);
  checkVarmapackSuccess(false);
  expect_count("failed API outcome without error", 1);
  checkVarmapackFailure(false);
  expect_count("failure missing error message", 1);
  checkVarmapackFailure(true);
  expect_count("success instead of failure", 1);
  checkVarmapackNaN(NAN);
  expect_count("NaN missing error message", 1);
  checkVarmapackNaN(0);
  expect_count("number instead of NaN without error", 1);
  checkVarmapackClean();
  expect_count("clean error state", 0);
  set_error_state(true);
  checkVarmapackSuccess(true);
  expect_count("success with stale error", 1);
  checkVarmapackSuccess(false);
  expect_count("failure instead of success", 1);
  checkVarmapackFailure(false);
  expect_count("expected failure", 0);
  checkVarmapackFailure(true);
  expect_count("success with stale error instead of failure", 1);
  checkVarmapackNaN(NAN);
  expect_count("expected NaN", 0);
  checkVarmapackNaN(0);
  expect_count("number instead of NaN with error", 1);
  checkVarmapackClean();
  expect_count("stale error state", 1);
  return failures != 0;
}
