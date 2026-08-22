#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "ExtraUtil.h"
#include "xCheck.h"
#include "Tests.h"
#include "varmapack.h"
#define DEBUG
#include "printX.h"
#include "debugprint.h"
static int contains(const char *hay, const char *needle) {
  return hay && needle && strstr(hay, needle) != 0;
}

static void check_symmetric_positive_diag(double A[], int n) {
  for (int j=0; j<n; j++) {
    xCheck(A[j + j*n] > 0);
    for (int i=0; i<n; i++) {
      xCheck(fabs(A[i + j*n] - A[j + i*n]) < 1e-14);
    }
  }
}

static void check_max(void) {
  int pmax = 99, qmax = -2, rmax = 7, icase_max;
  bool ok;
  ok = varmapack_testcase("max", &icase_max, 0, &pmax, &qmax, &rmax, 0, 0, 0, 0);
  checkVarmapackSuccess(ok);
  xCheck(icase_max == 16);
  xCheck(pmax == 5);
  xCheck(qmax == 3);
  xCheck(rmax == 7);
}

static int get_max(void) {
  int pmax = 0, qmax = 0, rmax = 0, icase_max;
  bool ok = varmapack_testcase("max", &icase_max, 0, &pmax, &qmax, &rmax,
                               0, 0, 0, 0);
  checkVarmapackSuccess(ok);
  return icase_max;
}

static void check_failures(void) {
  char name[VARMAPACK_TESTCASE_NAME_LEN] = "";
  int p = 0, q = 0, r = 1, icase = 1;
  int mone = -1, zero = 0;
  double A[1] = {0.0}, B[1] = {0.0}, Sig[1] = {0.0};
  int maxcase = get_max();
  int icase_m = -2, icase_p = maxcase + 1;
  bool ok;
  // (1) A = B = Sig = 0, name = "", icase = 1 should be ok
  // but the following should fail:
  // (2) p, q or r null pointers
  // (3) icase null
  // (4) name and icase null
  // (5) some A, B, Sig null others not
  // (6) A, B, Sig != 0, invalid dimensions
  // (7) name "max": A, B, Sig not all null
  // (8) unknown name
  // (9) name empty: icase < -1
  // (10) name empty: icase > max-case
  // (11) name "rho" with A = B = Sig = 0
  ok = varmapack_testcase(name, &icase, 0, &p, &q, &r, 0, 0, 0, 0);
  checkVarmapackSuccess(ok); // (1)
  name[0] = 0;
  ok = varmapack_testcase(name, &icase, 0, 0, &q, &r, A, B, 0, 0);
  checkVarmapackFailure(ok); // (2)
  ok = varmapack_testcase(name, 0, 0, &p, &q, &r, 0, 0, 0, 0);
  checkVarmapackFailure(ok); // (3)
  ok = varmapack_testcase(0, 0, 0, &p, &q, &r, 0, 0, 0, 0);
  checkVarmapackFailure(ok); // (4)
  ok = varmapack_testcase(name, &icase, 0, &p, &q, &r, 0, B, 0, 0);
  checkVarmapackFailure(ok); // (5)
  ok = varmapack_testcase("", &zero, 0, &p, &q, &zero, A, B, Sig, 0);
  checkVarmapackFailure(ok); // (6)
  ok = varmapack_testcase("", &zero, 0, &mone, &q, &r, A, B, Sig, 0);
  checkVarmapackFailure(ok); // (6)
  ok = varmapack_testcase("max", &icase, 0, &p, &q, &r, A, B, Sig, 0);
  checkVarmapackFailure(ok); // (7)
  ok = varmapack_testcase("blah", &icase, 0, &p, &q, &r, 0, 0, 0, 0);
  checkVarmapackFailure(ok); // (8)
  ok = varmapack_testcase(name, &icase_m, 0, &p, &q, &r, 0, 0, 0, 0);
  checkVarmapackFailure(ok); // (9)
  ok = varmapack_testcase(name, &icase_p, 0, &p, &q, &r, 0, 0, 0, 0);
  checkVarmapackFailure(ok); // (10)
  char rho_name[] = "rho";
  int rho_index = 1;
  ok = varmapack_testcase(rho_name, &rho_index, 0, &p, &q, &r, 0, 0, 0, 0);
  checkVarmapackFailure(ok); // (11)
  xCheck(strcmp(rho_name, "rho") == 0);
  ok = varmapack_testcase("max", &icase, 0, &p, &q, &r, 0, 0, 0, 0);
  checkVarmapackSuccess(ok);
}

static void check_single_case(int k) {
  char name[VARMAPACK_TESTCASE_NAME_LEN] = "";
  bool ok;
  int p = 0, q = 0, r = 0, icase = k;
  ok = varmapack_testcase(name, &icase, 0, &p, &q, &r, 0, 0, 0, 0);
  checkVarmapackSuccess(ok);
  xCheck(r > 0);
  if (contains(name, "ARMA"))    xCheck(p > 0 && q > 0);
  else if (contains(name, "MA")) xCheck(p == 0 && q > 0);
  else if (contains(name, "AR")) xCheck(p > 0 && q == 0);
  else xCheck(0);  // illegal name

  // re-query by name and confirm consistency
  int p2 = 0, q2 = 0, r2 = 0, icase2 = 0;
  ok = varmapack_testcase(name, &icase2, 0, &p2, &q2, &r2, 0, 0, 0, 0);
  checkVarmapackSuccess(ok);
  xCheck(p2 == p && q2 == q && r2 == r && icase2 == icase);
}

static void check_dims(void) {
  int icase_max = get_max();
  for (int k = 1; k <= icase_max; ++k)
    check_single_case(k);
}

static void check_construct_all_named(void) {
  int icase_max = get_max();
  for (int k=1; k<=icase_max; k++) {
    int p = 0, q = 0, r = 0, icase = k;
    char name[VARMAPACK_TESTCASE_NAME_LEN] = "";
    bool ok = varmapack_testcase(name, &icase, 0, &p, &q, &r, 0, 0, 0, 0);
    checkVarmapackSuccess(ok);
    int nA = r*r*(p > 0 ? p : 1);
    int nB = r*r*(q > 0 ? q : 1);
    int nS = r*r;
    double *A = 0, *B = 0, *Sig = 0;
    xCheck(ALLOC(A, nA));
    xCheck(ALLOC(B, nB));
    xCheck(ALLOC(Sig, nS));
    ok = varmapack_testcase(name, &icase, 0, &p, &q, &r, A, B, Sig, 0);
    checkVarmapackSuccess(ok);
    if (p > 0) checkArrayFinite(A, r*r*p);
    if (q > 0) checkArrayFinite(B, r*r*q);
    checkArrayFinite(Sig, r*r);
    check_symmetric_positive_diag(Sig, r);
    name[0] = 0;
    ok = varmapack_testcase(name, &icase, 0, &p, &q, &r, A, B, Sig, 0);
    checkVarmapackSuccess(ok);
    xCheck(strlen(name) > 0);
    FREE(Sig);
    FREE(B);
    FREE(A);
  }
}

static void check_deterministic_unnamed(void) {
  int p = 2, q = 1, r = 3, icase = -1;
  double A[18], B[9], Sig[9];
  bool ok = varmapack_testcase("", &icase, 0, &p, &q, &r, A, B, Sig, 0);
  checkVarmapackSuccess(ok);
  xCheck(p == 2 && q == 1 && r == 3);
  for (int i=0; i<18; i++) {
    xCheck(fabs(A[i] - 0.5/(2*3)) < 1e-15);
  }
  for (int i=0; i<9; i++) {
    xCheck(fabs(B[i] - 1.0/3) < 1e-15);
  }
  for (int j=0; j<3; j++) {
    for (int i=0; i<3; i++) {
      double expected = 1.0/(i + j + 1);
      if (i == j) expected += 0.2;
      xCheck(fabs(Sig[i + j*3] - expected) < 1e-15);
    }
  }
}

static void check_random_unnamed(void) {
  int p = 2, q = 1, r = 2, icase = 0;
  double A1[8], B1[4], Sig1[4], A2[8], B2[4], Sig2[4];
  randompack_rng *rng = randompack_create(0);
  bool ok = varmapack_testcase("", &icase, 0, &p, &q, &r, A1, B1, Sig1, 0);
  checkVarmapackFailure(ok);
  xCheck(rng != 0);
  xCheck(randompack_seed(123, 0, 0, rng));
  ok = varmapack_testcase("", &icase, 0, &p, &q, &r, A1, B1, Sig1, rng);
  checkVarmapackSuccess(ok);
  double rho = varmapack_specrad(A1, r, p);
  checkVarmapackClean();
  xCheck(rho < 1);
  xCheck(randompack_seed(123, 0, 0, rng));
  p = 2;
  q = 1;
  r = 2;
  icase = 0;
  ok = varmapack_testcase("", &icase, 0, &p, &q, &r, A2, B2, Sig2, rng);
  checkVarmapackSuccess(ok);
  checkArraySame(A1, A2, 8);
  checkArraySame(B1, B2, 4);
  checkArraySame(Sig1, Sig2, 4);
  check_symmetric_positive_diag(Sig1, 2);
  randompack_free(rng);
}

static void check_rho_case(void) {
  int p = 3, q = 3, r = 5, icase = 0;
  double A[75], B[75], Sig[25];
  double targets[] = {0.5, 0.95, 0.995, 2};
  for (int k=0; k<4; k++) {
    double target = targets[k];
    bool ok = varmapack_testcase("rho", &icase, target, &p, &q, &r,
                                               A, B, Sig, 0);
    double rho;
    checkVarmapackSuccess(ok);
    rho = varmapack_specrad(A, r, p);
    checkVarmapackClean();
    xCheck(fabs(rho - target) < 1e-6);
  }
  p = 3;
  q = 1;
  r = 2;
  icase = 0;
  bool ok = varmapack_testcase("rho", &icase, 0, &p, &q, &r, A, B, Sig, 0);
  checkVarmapackSuccess(ok);
  for (int i=0; i<r*r*p; i++) xCheck(A[i] == 0);
  p = 0;
  q = 1;
  r = 2;
  icase = 0;
  ok = varmapack_testcase("rho", &icase, 0.5, &p, &q, &r, A, B, Sig, 0);
  checkVarmapackFailure(ok);
  p = 3;
  q = 1;
  r = 2;
  icase = 0;
  ok = varmapack_testcase("rho", &icase, INFINITY, &p, &q, &r, A, B, Sig, 0);
  checkVarmapackFailure(ok);
}

void check_inquiry(void) {
  check_max();
  check_failures();
  check_dims();
}

static void check_construction(void) {
  check_construct_all_named();
  check_deterministic_unnamed();
  check_random_unnamed();
  check_rho_case();
}

void TestTestcase(void) {
  check_inquiry();
  check_construction();
}
