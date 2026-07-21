#include <stdio.h>
#include "Tests.h"
#include "error.h"
#include "varmapack.h"
#include "xCheck.h"
#include "ExtraUtil.h"
#include "VarmaPackUtil.h"
#include "VarmaUtilities.h"
#include "BlasGateway.h"

static void TestFindPsi(void) {
  int p, q, r, n, icase, h;
  char name[VARMAPACK_TESTCASE_NAME_LEN] = "";
  varmapack_error error;
  double *A, *B, *Sig, *Psi;

  // Query dimensions for mediumARMA1
  icase = 12;
  error = varmapack_testcase(name, &icase, 0, &p, &q, &r, 0, 0, 0, 0);
  xCheck(!error);
  xCheck(r == 3 && p == 3 && q == 3);
  h = imax(p, q);

  // Allocate and load testcase data
  xCheck(ALLOC(A, r*r*p));
  xCheck(ALLOC(B, r*r*q));
  xCheck(ALLOC(Sig, r*r));
  xCheck(ALLOC(Psi, r*h*r*h));

  error = varmapack_testcase(name, &icase, 0, &p, &q, &r, A, B, Sig, 0);
  xCheck(!error);

  FindPsi(A, B, Psi, p, q, r);

  // Expected Psi from Matlab find_Psi(A,B) for mediumARMA1 (3x9x3 flattened)
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

  FREE(A); FREE(B); FREE(Sig); FREE(Psi);
}

static void TestFindPsiHat(void) {
  int p, q, r, h, icase;
  char name[VARMAPACK_TESTCASE_NAME_LEN] = "";
  double *A, *B, *Sig, *Psi, *Psi_hat;
  varmapack_error error;
  
  icase = 12;
  error = varmapack_testcase(name, &icase, 0, &p, &q, &r, 0, 0, 0, 0);
  xCheck(!error);
  h = imax(p, q);
  
  xCheck(ALLOC(A, r*r*p));
  xCheck(ALLOC(B, r*r*q));
  xCheck(ALLOC(Sig, r*r));
  xCheck(ALLOC(Psi, r*h*r*h));
  xCheck(ALLOC(Psi_hat, r*h*r*h));

  error = varmapack_testcase(name, &icase, 0, &p, &q, &r, A, B, Sig, 0);
  xCheck(!error);
  FindPsi(A, B, Psi, p, q, r);
  double Sig_mod[] = {
    1.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 2.0
  };
  FindPsiHat(Psi, Psi_hat, Sig_mod, r, h);

  // Expected Psi_hat from Matlab find_Psi_hat for mediumARMA1 with modified Sig:
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

  FREE(A); FREE(B); FREE(Sig); FREE(Psi); FREE(Psi_hat);
}

static void TestPsiTinyAR(void) {
  // Matlab comparison
  double A[1], B[1], Sig[1], Psi[1], PsiHat[1];
  int p=1, q=0, r=1, icase, h=imax(p, q);
  varmapack_error error;
  error = varmapack_testcase("tinyAR", &icase, 0, &p, &q, &r, A, B, Sig, 0);
  xCheck(!error);
  FindPsi(A, B, Psi, p, q, r);
  FindPsiHat(Psi, PsiHat, Sig, r, h);
  xCheck(Psi[0] == 1);
  xCheck(almostSame(PsiHat[0], 0.894427190999916));
}

static void TestPsiSimple(void) {
  // Another Matlab comparison
  int p, q, r, icase, h;
  char name[VARMAPACK_TESTCASE_NAME_LEN] = "";
  varmapack_error error;
  double *A, *B, *Sig, *Psi, *Psi_hat;
  double A2[] = {0.1, 0.3, 0.2, 0.4};

  icase = 9;
  error = varmapack_testcase(name, &icase, 0, &p, &q, &r, 0, 0, 0, 0);
  xCheck(!error);
  xCheck(r == 2);
  h = imax(p, q);

  xCheck(ALLOC(A, r*r*(p+1)));
  xCheck(ALLOC(B, r*r*q));
  xCheck(ALLOC(Sig, r*r));
  xCheck(ALLOC(Psi, r*h*r*h));
  xCheck(ALLOC(Psi_hat, r*h*r*h));

  error = varmapack_testcase(name, &icase, 0, &p, &q, &r, A, B, Sig, 0);
  xCheck(!error);
  xCheck(r == 2);

  copy(4, A2, 1, A + 4, 1);  // A := [A A2]
  p = 2;
  Sig[0] = 1.0;              // So LSig is simple

  FindPsi(A, B, Psi, p, q, r);
  FindPsiHat(Psi, Psi_hat, Sig, r, h);
  
  double Psi_exp[] = {
    1,  0, .6, .4, 0,  1, .4, .4, 0,  0, 1,  0, 0,  0, 0,  1};
  double Psi_hat_exp[] = {
    1, 1, 1,  .8, 0, 1, .4, .4, 0, 0, 1,  1, 0, 0, 0,  1};

  xCheck(almostEqual(Psi, Psi_exp, r*h*r*h));
  xCheck(almostEqual(Psi_hat, Psi_hat_exp, r*h*r*h));
  FREE(B); FREE(Sig); FREE(A); FREE(Psi); FREE(Psi_hat);
}

static void TestPublicPsiScalar(void) {
  double A[] = {0.5};
  double B[] = {0.2};
  double Psi[4];
  double expected[] = {1, 0.7, 0.35, 0.175};
  varmapack_error error;
  error = varmapack_psi(A, B, 1, 1, 1, 3, Psi);
  xCheck(!error);
  xCheck(almostEqual(Psi, expected, 4));
}

static void TestPublicPsiMatrix(void) {
  double A[] = {0.1, 0.3, 0.2, 0.4};
  double B[] = {0.5, 0.7, 0.6, 0.8};
  double Psi[12];
  double expected[] = {
    1, 0, 0, 1, 0.6, 1.0, 0.8, 1.2, 0.26, 0.58, 0.32, 0.72
  };
  varmapack_error error;
  error = varmapack_psi(A, B, 1, 1, 2, 2, Psi);
  xCheck(!error);
  xCheck(almostEqual(Psi, expected, 12));
}

static void TestPublicPsiErrors(void) {
  double A[] = {0.5};
  double B[] = {0.2};
  double Psi[2];
  varmapack_error error;
  error = varmapack_psi(0, B, 1, 1, 1, 1, Psi);
  xCheck(error == VARMAPACK_INVALID_ARGUMENT);
  error = varmapack_psi(A, 0, 1, 1, 1, 1, Psi);
  xCheck(error == VARMAPACK_INVALID_ARGUMENT);
  error = varmapack_psi(A, B, -1, 1, 1, 1, Psi);
  xCheck(error == VARMAPACK_INVALID_ARGUMENT);
  error = varmapack_psi(A, B, 1, 1, 0, 1, Psi);
  xCheck(error == VARMAPACK_INVALID_ARGUMENT);
  error = varmapack_psi(A, B, 1, 1, 1, -1, Psi);
  xCheck(error == VARMAPACK_INVALID_ARGUMENT);
  error = varmapack_psi(A, B, 1, 1, 1, 1, 0);
  xCheck(error == VARMAPACK_INVALID_ARGUMENT);
}

static void TestPublicIrfMatrix(void) {
  double A[] = {0.1, 0.3, 0.2, 0.4};
  double B[] = {0.5, 0.7, 0.6, 0.8};
  double Sig[] = {4, 0, 0, 9};
  double Theta[12];
  double expected[] = {
    2, 0, 0, 3, 1.2, 2.0, 2.4, 3.6, 0.52, 1.16, 0.96, 2.16
  };
  varmapack_error error;
  error = varmapack_irf(A, B, Sig, 1, 1, 2, 2, Theta);
  xCheck(!error);
  xCheck(almostEqual(Theta, expected, 12));
}

static void TestPublicIrfSingularPSD(void) {
  double A[] = {0.1, 0.3, 0.2, 0.4};
  double B[] = {0.5, 0.7, 0.6, 0.8};
  double Sig1[] = {1, 2, 2, 4};
  double Sig2[] = {1, 1, 1, 1};
  double Theta[12];
  double Recon[4];
  varmapack_error error;
  error = varmapack_irf(A, B, Sig1, 1, 1, 2, 2, Theta);
  xCheck(!error);
  syrk("Low", "NoT", 2, 2, 1, Theta, 2, 0, Recon, 2);
  copylowertoupper(2, Recon, 2);
  xCheck(almostEqual(Recon, Sig1, 4));
  error = varmapack_irf(A, B, Sig2, 1, 1, 2, 2, Theta);
  xCheck(!error);
  syrk("Low", "NoT", 2, 2, 1, Theta, 2, 0, Recon, 2);
  copylowertoupper(2, Recon, 2);
  xCheck(almostEqual(Recon, Sig2, 4));
}

static void TestPublicIrfErrors(void) {
  double A[] = {0.5};
  double B[] = {0.2};
  double Sig[] = {-1};
  double Theta[2];
  varmapack_error error;
  error = varmapack_irf(A, B, 0, 1, 1, 1, 1, Theta);
  xCheck(error == VARMAPACK_INVALID_ARGUMENT);
  error = varmapack_irf(A, B, Sig, 1, 1, 1, 1, 0);
  xCheck(error == VARMAPACK_INVALID_ARGUMENT);
  error = varmapack_irf(A, B, Sig, 1, 1, 1, 1, Theta);
  xCheck(error == VARMAPACK_NOT_POSITIVE_SEMIDEFINITE);
}

// -----------------------------------------------------------------------------
// Public entry for this test module
// -----------------------------------------------------------------------------
void TestPsi(void) {
  TestPublicPsiScalar();
  TestPublicPsiMatrix();
  TestPublicPsiErrors();
  TestPublicIrfMatrix();
  TestPublicIrfSingularPSD();
  TestPublicIrfErrors();
  TestPsiTinyAR();
  TestPsiSimple();
  TestFindPsi();
  TestFindPsiHat();
}
