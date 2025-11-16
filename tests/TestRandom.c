// TestRandom.c — unit tests for RandomNumbers using ExtraUtil helpers

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "xCheck.h"

#include "RandomNumbers.h"
#include "VarmaUtilities.h"
#include "BlasGateway.h"
#include "ExtraUtil.h"    // for mean, var, relabsdiff, almostSame, almostEqual
#include "allocate.h"
#include "printX.h"

extern int TESTVERBOSITY;
static void msg(char *message) {
  if (TESTVERBOSITY >= 2) {printMsg(""); printMsgUpper(message);}
}

// Test uniform generator
static void test_uniform_basic(void) {
  const int N = 1e6;
  double meantol = 7*1/sqrt(12*N);  // 7 sigma where sigma = 1/sqrt(12·N)
  double vartol = 7*1/sqrt(180*N);  // 7 std of variance of U(0,1)
  double exactmu = 0.5, exactvar = 1.0/12.0;
  double *x = NULL;
  allocate(x, N);

  RandRng *rng = RandCreate();
  RandSetDefault(rng);
  Rand(x, N, rng);

  // Check range
  double xmin = x[0], xmax = x[0];
  for (int i = 1; i < N; i++) {
    if (x[i] < xmin) xmin = x[i];
    if (x[i] > xmax) xmax = x[i];
  }
  xCheck(xmin >= 0.0 && xmax < 1.0);
  double mu = mean(x, N);
  double va = var(x, N, mu);
  xCheck(fabs(mu - exactmu) < meantol);
  xCheck(fabs(va - exactvar) < vartol);

  RandFree(rng);
  freem(x);
}

// Test normal generator
static void test_normal_basic(void) {
  const int N = 1e6;
  double meantol = 7*1/sqrt(N);
  double vartol = 7*sqrt(2.0/(N-1)); // 7 std of variance of U(0,1)
  double exactmu = 0.0, exactvar = 1.0;
  double *x = NULL;
  allocate(x, N);

  RandRng *rng = RandCreate();
  RandN(x, N, rng);
  double mu = mean(x, N);
  double va = var(x, N, mu);
  printD("va", va);
  printD("exactvar", exactvar);
  printD("vartol", vartol);
  xCheck(fabs(mu - exactmu) < meantol);
  xCheck(fabs(va - exactvar) < vartol);

  RandFree(rng);
  freem(x);
}

// Test 

// Test determinism under same seed and RNG type
static void test_determinism_default_seed(void) {
  const int N = 1000;
  double *a = NULL, *b = NULL;
  allocate(a, N);
  allocate(b, N);

  RandRng *r1 = RandCreate();
  RandRng *r2 = RandCreate();
  RandSetDefault(r1);
  RandSetDefault(r2);
  RandSeed(123, r1);
  RandSeed(123, r2);

  RandN(a, N, r1);
  RandN(b, N, r2);

  xCheck(almostEqual(a, b, N));

  RandFree(r1);
  RandFree(r2);
  freem(a);
  freem(b);
}

// Test switching RNG type affects output
static void test_pm_vs_default_selection(void) {
  const int N = 256;
  double *a = NULL, *b = NULL;
  allocate(a, N);
  allocate(b, N);

  RandRng *r = RandCreate();

  RandSetDefault(r);
  RandSeed(123, r);
  Rand(a, N, r);

  RandSetPM(r);
  RandSeed(123, r);
  Rand(b, N, r);

  xCheck(!almostEqual(a, b, N));

  RandFree(r);
  freem(a);
  freem(b);
}



// Test randomize changes stream
static void test_randomize_changes_stream(void) {
  const int N = 128;
  double *a, *b, *c;
  allocate(a, N);
  allocate(b, N);
  allocate(c, N);

  RandRng *r = RandCreate();
  RandSetDefault(r);

  RandSeed(42, r);
  Rand(a, N, r);
  RandRandomize(r);
  Rand(b, N, r);
  RandThreadRandomize(42, r);
  Rand(c, N, r);

  xCheck(!almostEqual(a, b, N));
  xCheck(!almostEqual(a, c, N));
  xCheck(!almostEqual(b, c, N));

  RandFree(r);
  freem(a);
  freem(b);
  freem(c);
}

static void test_range(int m, int n, double A[], int k, double Y[]) {
  // True iff every column of the k × m Y is in the range of the m × n A
  double *B, *X, *R;
  int info;
  print3I("m, n, k", m, n, k);
  printM("A", A, m, n);
  printM("Y", Y, k, m);
  allocate(B, n*n);
  allocate(X, n*k);
  allocate(R, m*k);
  syrk("Lower", "T", n, m, 1.0, A, m, 1.0, B, n);       // B := A'A
  gemm("T", "T", n, k, m, 1.0, A, m, Y, k, 0.0, X, n);  // X := A'Y'
  posv("Lower", n, k, B, n, X, n, &info);               // Solve BX = A'Y
  xCheck(info == 0);
  gemm("N", "N", m, k, n, 1.0, A, m, X, n, 0.0, R, m); // R := AX
  printM("R", R, m, k);
  printM("Y", Y, k, m);
  subtracttranspose(m, k, Y, k, R, m);
  xCheck(almostZero(R, m*k));
  freem(R); freem(X); freem(B);
}

static void stdevs(double Sig[], double stds[], int n, int N) {
  for (int k=0, kk=0; k<n; k++, kk+=(n+1))
    stds[k] = sqrt(Sig[kk]/N);
}

static bool ok7sig(double x[], double sig[], int n) {
  // return true iff |x[k]| ≤ 7·sig for all k
  for (int k=0; k<n; k++) if (fabs(x[k]) > 7*sig[k]) return false;
  return true;
}

static void check_cov(char *transp, int N, double Sig[], double X[], RandRng *rng) {
  double C[16], Sms[4], Smsstd[4], Cmc[4], Cmcstd[4];
  int k, kk, ii, ij, jj;
  RandNM(transp, 0, Sig, 4, N, X, 0, rng);
  if (transp[0] == 'N')
    cov("N", N, 4, X, C);
  else
    cov("T", 4, N, X, C);
  for (k=0; k<4; k++) {
    kk = k*5;
    Sms[k] = fabs(C[kk] - Sig[kk]);
    Smsstd[k] = Sig[kk]*sqrt(2.0/(N-1));
  }
  xCheck(ok7sig(Sms, Smsstd, 4));
  k = 0;
  for (int i=0; i<4; i++) {
    for (int j=0; j<i; j++) {
      ii = i + 4*i;
      ij = i + 4*j;
      jj = j + 4*j;
      Cmc[k] = fabs(C[ij] - Sig[ij]);
      Cmcstd[k] = sqrt((Sig[ii]*Sig[jj] + Sig[ij]*Sig[ij])/(N-1));
      k++;
    }
  }
  printI("k", k);
  printV("Cmc", Cmc, 6);
  printV("Cmcstd", Cmcstd, 6);
  xCheck(ok7sig(Cmc, Cmcstd, 6));
}

static void test_randnm(double Sig[], int rank) {
  // Prepare:
  if (rank == 4) msg("Positive definite sigma");
  if (rank < 4)  msg("Positive semidefinite sigma");
  const int N = 10000, N1 = 10, N2 = 5;
  double mu[4] = {5, 10, 15, 20};
  double X[N*4], X1[N1*4], X2[N1*4], X3[N1*4], LSig[16], S[16];
  double means[4], meanstd_N[4];
  printM("Sig", Sig, 4, 4);
  RandRng *rng = RandCreate();
  
  stdevs(Sig, meanstd_N, 4, N); // Stddev of means:
  
  msg("Check that singular Sig gives singular LSig");
  RandSeed(9, rng);
  RandNM("NoT", mu, Sig, 4, N1, X1, LSig, rng);
  xCheck(rank == 4 ? LSig[15] > 0 : LSig[15] == 0);

  msg("Check LSig·LSig' = Sig:");
  laset("Upper", 4, 4, 0.0, 0.0, S, 4);
  syrk("Lower", "Not", 4, 4, 1.0, LSig, 4, 0.0, S, 4.0);
  printM("S", S, 4, 4);
  xCheck(almostEqual(Sig, S, 16));
  msg("LSig 0:");

  msg("Reuse LSig (Sig=0):");
  RandSeed(9, rng);
  RandNM("NoT", mu, 0, 4, N1, X2, LSig, rng);
  xCheck(almostEqual(X1, X2, N1*4));

  msg("Check setting seed to same value");
  RandSeed(9, rng);
  RandNM("NoT", mu, Sig, 4, N1, X3, 0, rng);
  xCheck(almostEqual(X2, X3, N1*4));
  
  
  msg("-Check that X is in the range of LSig:");
  RandNM("NoT", 0, Sig, 4, N2, X, LSig, rng);
  test_range(4, rank, LSig, N2, X);
  
  msg("Check correct means with specified mu:");
  RandSeed(9, rng);
  RandNM("NoT", mu, 0, 4, N, X, LSig, rng);
  meanmat("N", N, 4, X, N, means);
  //printM("(a) X", X, N2, 4);
  printV("means", means, 4);
  axpy(4, -1.0, mu, 1, means, 1);
  printV("mu", mu, 4);
  printV("meanstd", meanstd_N, 4);
  xCheck(ok7sig(means, meanstd_N, 4));
  RandSeed(9, rng);
  RandNM("T", mu, 0, 4, N, X, LSig, rng); // also for X^T
  //printM("(b) X", X, 4, N2);
  meanmat("T", 4, N, X, 4, means);
  printV("means", means, 4);
  axpy(4, -1.0, mu, 1, means, 1);
  xCheck(ok7sig(means, meanstd_N, 4));
  
  msg("– and with mu=0:");
  RandNM("NoT", 0, Sig, 4, N, X, LSig, rng);
  meanmat("N", N, 4, X, N, means);
  printV("means", means, 4);
  xCheck(ok7sig(means, meanstd_N, 4));
  
  msg("-Sample covariance");
  check_cov("NoT", N, Sig, X, rng);
  check_cov("T", N, Sig, X, rng);
  RandFree(rng);
}

static void test_multivariate_normal(void) {
  double Sig[16],
    Lnonsing[16] = {
      1, 2, 3, 2,
      0, 2, 1, 1,
      0, 0, 2, 0,
      0, 0, 0, 1},
    Lsing[16] = {
      1, 2, 3, 2,
      0, 2, 1, 1};
  //scal(16, 0.1, Lnonsing, 1);
  xCheckAddMsg("Test positive definite Sigma");
  laset("Upper", 4, 4, 0.0, 0.0, Sig, 4);
  syrk("Lower", "NoT", 4, 4, 1.0, Lnonsing, 4, 0.0, Sig, 4);
  test_randnm(Sig, 4);
  //
  xCheckAddMsg("Test positive semidefinite Sigma");
  laset("Upper", 4, 4, 0.0, 0.0, Sig, 4);
  syrk("Lower", "NoT", 4, 4, 1.0, Lsing, 4, 0.0, Sig, 4);
  //test_randnm(Sig, 2);  
}

// static void run_test(const char *name, void (*fn)(void)) {
//   xCheckAddMsg(name);
//   fn();
// }

#define RUN_TEST(x)                             \
 do {                                           \
   xCheckAddMsg(#x);                            \
   test_##x();                                  \
 } while(0)

void TestRandomNumbers(void) {
  RUN_TEST(multivariate_normal);
  RUN_TEST(uniform_basic);
  RUN_TEST(normal_basic);
  RUN_TEST(determinism_default_seed);
  RUN_TEST(pm_vs_default_selection);
  RUN_TEST(randomize_changes_stream);
}
