#include <stdbool.h>
#include <limits.h>
#include "TestUtil.h"
#include "Tests.h"
#include "varmapack.h"

static void check_values(void) {
  double C[6], z[20];
  double expectedC[] = {0.3, -0.5, -0.3, 0, 0.2, 0.5};
  double expectedZ[] = {0.6, -0.2, 0.8, 0, -1, 0.4, 0.2, -0.6, 1, -0.4,
                        0, 0.8, -0.8, 0.2, -0.2, 0.4, -0.6, 0.6, -0.2, 0.8};
  bool ok = varmapack_testcasex(3, 1, 2, 20, C, z);
  checkVarmapackSuccess(ok);
  checkArraySame(C, expectedC, 6);
  checkArraySame(z, expectedZ, 20);
}

static void check_zero_terms(void) {
  double C;
  bool ok = varmapack_testcasex(0, 0, 2, 3, 0, 0);
  checkVarmapackSuccess(ok);
  ok = varmapack_testcasex(1, 1, 1, 0, &C, 0);
  checkVarmapackSuccess(ok);
  xCheck(C == 0.3);
}

static void check_vector_inputs(void) {
  double C[8], z[6];
  double expectedC[] = {0.3, -0.5, -0.1, 0.2, -0.3, 0, 0.4, -0.4};
  double expectedZ[] = {0.6, -0.2, -0.2, 0.8, 0.8, 0};
  bool ok = varmapack_testcasex(2, 2, 2, 3, C, z);
  checkVarmapackSuccess(ok);
  checkArraySame(C, expectedC, 8);
  checkArraySame(z, expectedZ, 6);
}

static void check_failures(void) {
  double C, z;
  bool ok = varmapack_testcasex(-1, 1, 1, 1, &C, &z);
  checkVarmapackFailure(ok);
  ok = varmapack_testcasex(1, 0, 1, 1, &C, &z);
  checkVarmapackFailure(ok);
  ok = varmapack_testcasex(1, 1, -1, 1, &C, &z);
  checkVarmapackFailure(ok);
  ok = varmapack_testcasex(1, 1, 1, 1, 0, &z);
  checkVarmapackFailure(ok);
  ok = varmapack_testcasex(1, 1, 1, 1, &C, 0);
  checkVarmapackFailure(ok);
  ok = varmapack_testcasex(2, 1, INT_MAX, 0, &C, 0);
  checkVarmapackFailure(ok);
  ok = varmapack_testcasex(0, 0, 1, 0, 0, 0);
  checkVarmapackSuccess(ok);
}

void TestTestcasex(void) {
  check_values();
  check_zero_terms();
  check_vector_inputs();
  check_failures();
}
