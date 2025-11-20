// TestRandom.c — unit tests for RandomNumbers using ExtraUtil helpers

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "xCheck.h"

#include "randompack.h"
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

  randompack_rng *rng = randompack_create("Xorshift", 0);
  randompack_u01(x, N, rng);

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

  randompack_free(rng);
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

  randompack_rng *rng = randompack_create("Xorshift", 0);
  randompack_norm(x, N, rng);
  double mu = mean(x, N);
  double va = var(x, N, mu);
  printD("va", va);
  printD("exactvar", exactvar);
  printD("vartol", vartol);
  xCheck(fabs(mu - exactmu) < meantol);
  xCheck(fabs(va - exactvar) < vartol);

  randompack_free(rng);
  freem(x);
}

// Test 

// Test determinism under same seed and RNG type
static void test_determinism_default_seed(void) {
  const int N = 1000;
  double *a = NULL, *b = NULL;
  allocate(a, N);
  allocate(b, N);

  randompack_rng *r1 = randompack_create("Xorshift", 123);
  randompack_rng *r2 = randompack_create("Xorshift", 123);

  randompack_norm(a, N, r1);
  randompack_norm(b, N, r2);

  xCheck(almostEqual(a, b, N));

  randompack_free(r1);
  randompack_free(r2);
  freem(a);
  freem(b);
}

// Test switching RNG type affects output
static void test_pm_vs_default_selection(void) {
  const int N = 256;
  double *a = NULL, *b = NULL;
  allocate(a, N);
  allocate(b, N);

  randompack_rng *r1 = randompack_create("Xorshift", 123);
  randompack_u01(a, N, r1);

  randompack_rng *r2 = randompack_create("PM", 123);
  randompack_u01(b, N, r2);

  xCheck(!almostEqual(a, b, N));

  randompack_free(r1);
  randompack_free(r2);
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

  randompack_rng *r1 = randompack_create("Xorshift", 42);
  randompack_u01(a, N, r1);

  randompack_rng *r2 = randompack_create("Xorshift", 0);
  randompack_u01(b, N, r2);

  randompack_rng *r3 = randompack_create("Xorshift", -42);
  randompack_u01(c, N, r3);

  xCheck(!almostEqual(a, b, N));
  xCheck(!almostEqual(a, c, N));
  xCheck(!almostEqual(b, c, N));

  randompack_free(r1);
  randompack_free(r2);
  randompack_free(r3);
  freem(a);
  freem(b);
  freem(c);
}

static void test_int_api(void) {
  const int N = 128;
  int a[N], b[N];
  randompack_rng *r1 = randompack_create("Xorshift", 99);
  randompack_rng *r2 = randompack_create("Xorshift", 99);
  randompack_int(a, N, -3, 8, r1);
  randompack_int(b, N, 0, 11, r2);
  for (int i = 0; i < N; i++) {
    xCheck(a[i] >= -3 && a[i] <= 8);
    xCheck(b[i] >= 0 && b[i] <= 11);
    xCheck(a[i] + 3 == b[i]);
  }
  randompack_free(r1);
  randompack_free(r2);
}

static void test_perm_api(void) {
  const int N = 32;
  int perm[N];
  int seen[N];
  for (int i = 0; i < N; i++) seen[i] = 0;
  randompack_rng *rng = randompack_create("Xorshift", 77);
  randompack_perm(perm, N, rng);
  for (int i = 0; i < N; i++) {
    int v = perm[i];
    xCheck(v >= 0 && v < N);
    xCheck(!seen[v]);
    seen[v] = 1;
  }
  randompack_free(rng);
}

static void test_sample_api(void) {
  const int N = 50;
  const int K = 10;
  int sample1[K], sample2[K];
  int used[N];
  randompack_rng *r1 = randompack_create("Xorshift", 11);
  randompack_rng *r2 = randompack_create("Xorshift", 11);
  randompack_sample(sample1, N, K, r1);
  randompack_sample(sample2, N, K, r2);
  for (int i = 0; i < N; i++) used[i] = 0;
  for (int i = 0; i < K; i++) {
    int v = sample1[i];
    xCheck(v >= 0 && v < N);
    xCheck(!used[v]);
    used[v] = 1;
    xCheck(sample1[i] == sample2[i]);
  }
  randompack_rng *r3 = randompack_create("Xorshift", 5);
  randompack_sample(0, N, 0, r3); // ensure k=0 is a no-op
  randompack_free(r1);
  randompack_free(r2);
  randompack_free(r3);
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

static void check_cov(char *transp, int N, double Sig[], double X[], randompack_rng *rng) {
  double C[16], Sms[4], Smsstd[4], Cmc[4], Cmcstd[4];
  int k, kk, ii, ij, jj;
  randompack_mvn(transp, 0, Sig, 4, N, X, 0, rng);
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
  double *X, *X1, *X2, *X3, LSig[16], S[16];
  double means[4], meanstd_N[4];
  allocate(X, N*4);
  allocate(X1, N1*4);
  allocate(X2, N1*4);
  allocate(X3, N1*4);
  printM("Sig", Sig, 4, 4);
  randompack_rng *rng;

  stdevs(Sig, meanstd_N, 4, N); // Stddev of means:

  msg("Check that singular Sig gives singular LSig");
  rng = randompack_create("Xorshift", 9);
  randompack_mvn("NoT", mu, Sig, 4, N1, X1, LSig, rng);
  xCheck(rank == 4 ? LSig[15] > 0 : LSig[15] == 0);
  randompack_free(rng);

  msg("Check LSig·LSig' = Sig:");
  laset("Upper", 4, 4, 0.0, 0.0, S, 4);
  syrk("Lower", "Not", 4, 4, 1.0, LSig, 4, 0.0, S, 4.0);
  printM("S", S, 4, 4);
  xCheck(almostEqual(Sig, S, 16));
  msg("LSig 0:");

  msg("Reuse LSig (Sig=0):");
  rng = randompack_create("Xorshift", 9);
  randompack_mvn("NoT", mu, 0, 4, N1, X2, LSig, rng);
  xCheck(almostEqual(X1, X2, N1*4));
  randompack_free(rng);

  msg("Check setting seed to same value");
  rng = randompack_create("Xorshift", 9);
  randompack_mvn("NoT", mu, Sig, 4, N1, X3, 0, rng);
  xCheck(almostEqual(X2, X3, N1*4));
  randompack_free(rng);

  msg("-Check that X is in the range of LSig:");
  rng = randompack_create("Xorshift", 0);
  randompack_mvn("NoT", 0, Sig, 4, N2, X, LSig, rng);
  test_range(4, rank, LSig, N2, X);
  randompack_free(rng);

  msg("Check correct means with specified mu:");
  rng = randompack_create("Xorshift", 9);
  randompack_mvn("NoT", mu, 0, 4, N, X, LSig, rng);
  meanmat("N", N, 4, X, N, means);
  //printM("(a) X", X, N2, 4);
  printV("means", means, 4);
  axpy(4, -1.0, mu, 1, means, 1);
  printV("mu", mu, 4);
  printV("meanstd", meanstd_N, 4);
  xCheck(ok7sig(means, meanstd_N, 4));
  randompack_free(rng);

  rng = randompack_create("Xorshift", 9);
  randompack_mvn("T", mu, 0, 4, N, X, LSig, rng); // also for X^T
  //printM("(b) X", X, 4, N2);
  meanmat("T", 4, N, X, 4, means);
  printV("means", means, 4);
  axpy(4, -1.0, mu, 1, means, 1);
  xCheck(ok7sig(means, meanstd_N, 4));
  randompack_free(rng);

  msg("– and with mu=0:");
  rng = randompack_create("Xorshift", 0);
  randompack_mvn("NoT", 0, Sig, 4, N, X, LSig, rng);
  meanmat("N", N, 4, X, N, means);
  printV("means", means, 4);
  xCheck(ok7sig(means, meanstd_N, 4));

  msg("-Sample covariance");
  check_cov("NoT", N, Sig, X, rng);
  check_cov("T", N, Sig, X, rng);
  randompack_free(rng);
  freem(X);
  freem(X1);
  freem(X2);
  freem(X3);
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
      0, 2, 1, 1,
      0, 0, 0, 0,
      0, 0, 0, 0};
  //scal(16, 0.1, Lnonsing, 1);
  xCheckAddMsg("Test positive definite Sigma");
  laset("Upper", 4, 4, 0.0, 0.0, Sig, 4);
  syrk("Lower", "NoT", 4, 4, 1.0, Lnonsing, 4, 0.0, Sig, 4);
  test_randnm(Sig, 4);

  xCheckAddMsg("Test positive semidefinite Sigma");
  laset("Upper", 4, 4, 0.0, 0.0, Sig, 4);
  syrk("Lower", "NoT", 4, 4, 1.0, Lsing, 4, 0.0, Sig, 4);
  test_randnm(Sig, 2);  
}

// Test type name aliases produce identical output
static void test_type_name_aliases(void) {
  const int N = 100;
  double *a = NULL, *b = NULL, *c = NULL, *d = NULL, *e = NULL;
  allocate(a, N);
  allocate(b, N);
  allocate(c, N);
  allocate(d, N);
  allocate(e, N);

  randompack_rng *r1 = randompack_create("Xorshift128+", 123);
  randompack_u01(a, N, r1);

  randompack_rng *r2 = randompack_create("Xorshift", 123);
  randompack_u01(b, N, r2);

  randompack_rng *r3 = randompack_create("X+", 123);
  randompack_u01(c, N, r3);

  xCheck(almostEqual(a, b, N));
  xCheck(almostEqual(a, c, N));

  randompack_free(r1);
  randompack_free(r2);
  randompack_free(r3);

  randompack_rng *r4 = randompack_create("Park-Miller", 456);
  randompack_u01(d, N, r4);

  randompack_rng *r5 = randompack_create("PM", 456);
  randompack_u01(e, N, r5);

  xCheck(almostEqual(d, e, N));

  randompack_free(r4);
  randompack_free(r5);

  randompack_rng *r6 = randompack_create("R", 789);
  randompack_u01(a, N, r6);

  randompack_rng *r7 = randompack_create("R-default", 789);
  randompack_u01(b, N, r7);

  xCheck(almostEqual(a, b, N));

  randompack_free(r6);
  randompack_free(r7);

  freem(a);
  freem(b);
  freem(c);
  freem(d);
  freem(e);
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
  RUN_TEST(int_api);
  RUN_TEST(perm_api);
  RUN_TEST(sample_api);
  RUN_TEST(type_name_aliases);
}
