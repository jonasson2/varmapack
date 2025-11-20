#include <stdio.h>
#include "Tests.h"
#include "allocate.h"
#include "xAssert.h"
#include "varmapack.h"
#include "xCheck.h"
#include "ExtraUtil.h"
#include "VarmaPackUtil.h"
#include "VarmaUtilities.h"
#include "BlasGateway.h"

static void TestFindPsi(void) {
  int p, q, r, n, icase, h;
  char name[16] = "";
  bool ok;
  double *A, *B, *Sig, *Psi;

  // Query dimensions for testcase 9
  icase = 9;
  ok = varmapack_testcase(0, 0, 0, name, &p, &q, &r, &icase, 0, stdout);
  xCheck(ok);
  xCheck(r == 3 && p == 3 && q == 3);
  h = imax(p, q);

  // Allocate and load testcase data
  allocate(A, r*r*p);
  allocate(B, r*r*q);
  allocate(Sig, r*r);
  allocate(Psi, r*h*r*h);

  ok = varmapack_testcase(A, B, Sig, name, &p, &q, &r, &icase, 0, stdout);
  xCheck(ok);

  vpack_FindPsi(A, B, Psi, p, q, r);

  // Expected Psi from Matlab find_Psi(A,B) for testcase 9 (3x9x3 flattened)
  double Psi_exp[] = {
    1.0000, 0.0000, 0.0000, 0.2100, 0.2400, 0.2600, 0.3485, 0.3756, 0.4027,
    0.0000, 1.0000, 0.0000, 0.1400, 0.1600, 0.1800, 0.3260, 0.3508, 0.3756,
    0.0000, 0.0000, 1.0000, 0.0600, 0.0800, 0.1000, 0.3020, 0.3244, 0.3468,
    0.0000, 0.0000, 0.0000, 1.0000, 0.0000, 0.0000, 0.2100, 0.2400, 0.2600,
    0.0000, 0.0000, 0.0000, 0.0000, 1.0000, 0.0000, 0.1400, 0.1600, 0.1800,
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 1.0000, 0.0600, 0.0800, 0.1000,
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 1.0000, 0.0000, 0.0000,
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 1.0000, 0.0000,
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 1.0000
  };

  n = r*r*h;
  xCheck(almostEqual(Psi, Psi_exp, n));

  freem(A); freem(B); freem(Sig); freem(Psi);
}

static void TestFindPsiHat(void) {
  int p, q, r, h, icase;
  char name[16] = "";
  double *A, *B, *Sig, *Psi, *Psi_hat;
  
  icase = 9;
  varmapack_testcase(0, 0, 0, name, &p, &q, &r, &icase, 0, stdout);
  h = imax(p, q);
  
  allocate(A, r*r*p);
  allocate(B, r*r*q);
  allocate(Sig, r*r);
  allocate(Psi, r*h*r*h);
  allocate(Psi_hat, r*h*r*h);

  varmapack_testcase(A, B, Sig, name, &p, &q, &r, &icase, 0, stdout);
  vpack_FindPsi(A, B, Psi, p, q, r);
  double Sig_mod[] = {
    1.0, 1.0, 1.0,
    1.0, 2.0, 1.0,
    1.0, 1.0, 2.0
  };
  vpack_FindPsiHat(Psi, Psi_hat, Sig_mod, r, h);

  // Expected Psi_hat from Matlab find_Psi_hat for testcase 9 with modified Sig:
  double Psi_hat_exp[] = {
    1.0000, 1.0000, 1.0000, 0.4100, 0.4800, 0.5400, 0.9765, 1.0508, 1.1251,
    0.0000, 1.0000, 0.0000, 0.1400, 0.1600, 0.1800, 0.3260, 0.3508, 0.3756,
    0.0000, 0.0000, 1.0000, 0.0600, 0.0800, 0.1000, 0.3020, 0.3244, 0.3468,
    0.0000, 0.0000, 0.0000, 1.0000, 1.0000, 1.0000, 0.4100, 0.4800, 0.5400,
    0.0000, 0.0000, 0.0000, 0.0000, 1.0000, 0.0000, 0.1400, 0.1600, 0.1800,
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 1.0000, 0.0600, 0.0800, 0.1000,
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 1.0000, 1.0000, 1.0000,
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 1.0000, 0.0000,
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 1.0000
  };
  xCheck(almostEqual(Psi_hat, Psi_hat_exp, r*h*r*h));

  freem(A); freem(B); freem(Sig); freem(Psi); freem(Psi_hat);
}

static void TestPsiTinyAR(void) {
  // Matlab comparison
  double A[1], B[1], Sig[1], Psi[1], PsiHat[1];
  int p=1, q=0, r=1, icase, h=imax(p, q);
  bool ok;
  ok = varmapack_testcase(A, B, Sig, "tinyAR", &p, &q, &r, &icase, 0, stdout);
  xCheck(ok);
  vpack_FindPsi(A, B, Psi, p, q, r);
  vpack_FindPsiHat(Psi, PsiHat, Sig, r, h);
  xCheck(Psi[0] == 1);
  xCheck(almostSame(PsiHat[0], 0.894427190999916));
}

static void TestPsiSimple(void) {
  // Another Matlab comparison
  int p, q, r, icase, h;
  char name[16] = "";
  bool ok;
  double *A, *B, *Sig, *Psi, *Psi_hat;
  double A2[] = {0.1, 0.3, 0.2, 0.4};

  icase = 7;
  ok = varmapack_testcase(0, 0, 0, name, &p, &q, &r, &icase, 0, stdout);
  xCheck(ok);
  xCheck(r == 2);
  h = imax(p, q);

  allocate(A, r*r*(p+1));
  allocate(B, r*r*q);
  allocate(Sig, r*r);
  allocate(Psi, r*h*r*h);
  allocate(Psi_hat, r*h*r*h);

  ok = varmapack_testcase(A, B, Sig, name, &p, &q, &r, &icase, 0, stdout);
  xCheck(ok);
  xCheck(r == 2);

  copy(4, A2, 1, A + 4, 1);  // A := [A A2]
  p = 2;
  Sig[0] = 1.0;              // So LSig is simple

  vpack_FindPsi(A, B, Psi, p, q, r);
  vpack_FindPsiHat(Psi, Psi_hat, Sig, r, h);
  
  double Psi_exp[] = {
    1,  0, .6, .4,
    0,  1, .4, .4,
    0,  0, 1,  0,
    0,  0, 0,  1};
  double Psi_hat_exp[] = {
    1, 1, 1,  .8,
    0, 1, .4, .4,
    0, 0, 1,  1,
    0, 0, 0,  1};

  xCheck(almostEqual(Psi, Psi_exp, r*h*r*h));
  xCheck(almostEqual(Psi_hat, Psi_hat_exp, r*h*r*h));
  freem(B); freem(Sig); freem(A); freem(Psi); freem(Psi_hat);
}

// -----------------------------------------------------------------------------
// Public entry for this test module
// -----------------------------------------------------------------------------
void TestPsi(void) {
  TestPsiTinyAR();
  TestPsiSimple();
  TestFindPsi();
  TestFindPsiHat();
}
