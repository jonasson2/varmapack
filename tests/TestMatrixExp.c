#include <math.h>
#include "TestUtil.h"
#include "Tests.h"
#include "BlasGateway.h"
#include "mb05odGateway.h"

static void check_diagonal_and_scalar(void) {
  double A[] = {-2, 0, 0, 0, .5, 0, 0, 0, 3}, E[9], expected[9] = {0};
  double saved[9], gaps[] = {0, 1e-12, .25, -1, 10};
  copy(9, A, 1, saved, 1);
  for (int k=0; k<LEN(gaps); k++) {
    mb05od_accuracy accuracy;
    for (int i=0; i<3; i++) expected[i+3*i] = exp(gaps[k]*A[i+3*i]);
    checkVarmapackSuccess(mb05odGateway(A, gaps[k], 3, E, &accuracy));
    test_check_group entries = {0};
    for (int i=0; i<9; i++)
      groupCheck(&entries, fabs(E[i]-expected[i]) <= 3e-13*fmax(1, fabs(expected[i])));
    checkGroup(&entries, "diagonal exponential, relative/absolute error");
    xCheck(accuracy.mdig >= 0 && accuracy.idig >= 0);
    xCheck(accuracy.iwarn >= 0 && accuracy.iwarn <= 3);
  }
  checkArraySame(A, saved, 9);
  double a = -2, e;
  checkVarmapackSuccess(mb05odGateway(&a, .5, 1, &e, 0));
  xCheck(fabs(e-exp(-1)) < 1e-15);
  checkVarmapackSuccess(mb05odGateway(&a, 1000, 1, &e, 0));
  xCheck(e == 0); // Underflow to zero is allowed.
  checkVarmapackSuccess(mb05odGateway(0, 1, 0, 0, 0));
}

static void check_rotation_and_jordan(void) {
  double rotation[] = {-1, 2, -2, -1}, E[4], expected[4];
  double gaps[] = {1e-12, .1, 1, 10};
  for (int k=0; k<LEN(gaps); k++) {
    double d = gaps[k], decay = exp(-d);
    expected[0] = expected[3] = decay*cos(2*d);
    expected[1] = decay*sin(2*d);
    expected[2] = -expected[1];
    checkVarmapackSuccess(mb05odGateway(rotation, d, 2, E, 0));
    checkArrayTol(E, expected, 4, 2e-13);
  }
  // A defective, nonnormal matrix must not rely on eigenvector diagonalization.
  double jordan[] = {-1, 0, 0, 2, -1, 0, 0, 3, -1}, J[9];
  double d = .7, scale = exp(-d);
  double exact[] = {scale, 0, 0, 2*d*scale, scale, 0,
                    3*d*d*scale, 3*d*scale, scale};
  checkVarmapackSuccess(mb05odGateway(jordan, d, 3, J, 0));
  checkArrayTol(J, exact, 9, 2e-13);
  // Strongly unequal off-diagonal magnitudes exercise balancing/unbalancing.
  double balanced[] = {0, 1e-8, 1e8, 0};
  double exp_balanced[] = {cosh(.3), 1e-8*sinh(.3), 1e8*sinh(.3), cosh(.3)};
  checkVarmapackSuccess(mb05odGateway(balanced, .3, 2, E, 0));
  test_check_group entries = {0};
  for (int i=0; i<4; i++)
    groupCheck(&entries, fabs(E[i]-exp_balanced[i]) <=
                         3e-13*fmax(1, fabs(exp_balanced[i])));
  checkGroup(&entries, "balanced exponential, relative/absolute error");
}

static void check_composition_and_alias(void) {
  double A[] = {-.8, .2, -.1, .4, -1.2, .3, -.2, .1, -.5};
  double E1[9], E2[9], E3[9], product[9], inverse[9], identity[9];
  checkVarmapackSuccess(mb05odGateway(A, .3, 3, E1, 0));
  checkVarmapackSuccess(mb05odGateway(A, .7, 3, E2, 0));
  checkVarmapackSuccess(mb05odGateway(A, 1, 3, E3, 0));
  gemm("N", "N", 3, 3, 3, 1, E1, 3, E2, 3, 0, product, 3);
  checkArrayTol(product, E3, 9, 2e-13);
  checkVarmapackSuccess(mb05odGateway(A, -1, 3, inverse, 0));
  gemm("N", "N", 3, 3, 3, 1, E3, 3, inverse, 3, 0, product, 3);
  laset("All", 3, 3, 0, 1, identity, 3);
  checkArrayTol(product, identity, 9, 2e-13);
  checkVarmapackSuccess(mb05odGateway(A, 1, 3, A, 0));
  checkArrayTol(A, E3, 9, 2e-13);
  double zero[9] = {0};
  checkVarmapackSuccess(mb05odGateway(zero, 4, 3, E1, 0));
  checkArrayTol(E1, identity, 9, 1e-15);
}

static void check_published_example(void) {
  // MB05OD documentation example; reference values are printed to four decimals.
  double A[] = {2, 0, 1, 1, 3, 0, 1, 2, 4}, E[9];
  double expected[] = {22.5984, 24.4047, 29.4097, 17.2073, 27.6033, 12.2024,
                       53.8144, 83.2241, 81.4177};
  mb05od_accuracy accuracy;
  checkVarmapackSuccess(mb05odGateway(A, 1, 3, E, &accuracy));
  checkArrayTol(E, expected, 9, 5.1e-5);
  xCheck(accuracy.mdig > 0 && accuracy.idig > 0);
}

static void check_errors(void) {
  double A[] = {-1, 0, 0, -2}, E[] = {7, 8, 9, 10}, saved[] = {7, 8, 9, 10};
  mb05od_accuracy accuracy = {99, 99, 99};
  checkVarmapackFailure(mb05odGateway(A, NAN, 2, E, &accuracy));
  xCheck(accuracy.mdig == 0 && accuracy.idig == 0 && accuracy.iwarn == 0);
  checkVarmapackFailure(mb05odGateway(A, INFINITY, 2, E, 0));
  checkVarmapackFailure(mb05odGateway(A, 1, -1, E, 0));
  checkVarmapackFailure(mb05odGateway(A, 1, INT_MAX, E, 0));
  checkVarmapackFailure(mb05odGateway(A, 1, 32768, E, 0));
  checkVarmapackFailure(mb05odGateway(0, 1, 2, E, 0));
  checkVarmapackFailure(mb05odGateway(A, 1, 2, 0, 0));
  A[1] = NAN;
  checkVarmapackFailure(mb05odGateway(A, 1, 2, E, 0));
  A[1] = INFINITY;
  checkVarmapackFailure(mb05odGateway(A, 1, 2, E, 0));
  A[1] = 1e308;
  checkVarmapackFailure(mb05odGateway(A, 2, 2, E, 0));
  checkArraySame(E, saved, 4);
  double a = 1000, e = 42;
  checkVarmapackFailure(mb05odGateway(&a, 1, 1, &e, 0));
  xCheck(e == 42); // Scalar MB05OD returns INFO=0 even if EXP overflows.
  checkVarmapackFailure(mb05odGateway(&a, 1, 1, &a, 0));
  xCheck(a == 1000);
  double huge[] = {1000, 0, 0, 1001};
  checkVarmapackFailure(mb05odGateway(huge, 1, 2, E, 0));
  checkArraySame(E, saved, 4);
  A[1] = 0;
  checkVarmapackSuccess(mb05odGateway(A, 1, 2, E, 0)); // Clears stale errors.
}

void TestMatrixExp(void) {
  check_diagonal_and_scalar();
  check_rotation_and_jordan();
  check_composition_and_alias();
  check_published_example();
  check_errors();
}
