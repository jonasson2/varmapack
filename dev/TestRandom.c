// TestRandom.c — unit tests for RandomNumbers using ExtraUtil helpers

#include <stdio.h>
#include <stdlib.h>
#include "xCheck.h"

#include "RandomNumbers.h"
#include "ExtraUtil.h"    // for mean, var, relabsdiff, almostSame, almostEqual
#include "allocate.h"
#include "printX.h"

// Test uniform generator
static void test_uniform_basic(void) {
  const int N = 1e6;
  double meantol = 0.002; // 7 sigma where sigma = 1/12/sqrt(N)
  double vartol = 0.0005; // 7 std of variance of U(0,1)
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
  free(x);
}

// Test normal generator
static void test_normal_basic(void) {
  const int N = 1e6;
  double meantol = 0.007; // 7 sigma where sigma = 1/12/sqrt(N)
  double vartol = 0.01; // 7 std of variance of U(0,1)
  double exactmu = 0.0, exactvar = 1.0;
  double *x = NULL;
  allocate(x, N);

  RandRng *rng = RandCreate();
  RandN(x, N, rng);
  double mu = mean(x, N);
  double va = var(x, N, mu);
  xCheck(fabs(mu - exactmu) < meantol);
  xCheck(fabs(va - exactvar) < vartol);

  RandFree(rng);
  free(x);
}

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
  free(a);
  free(b);
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
  free(a);
  free(b);
}

// Test randomize changes stream
static void test_randomize_changes_stream(void) {
  const int N = 128;
  double *a = NULL, *b = NULL;
  allocate(a, N);
  allocate(b, N);

  RandRng *r = RandCreate();
  RandSetDefault(r);

  RandSeed(42, r);
  Rand(a, N, r);
  RandRandomize(r);
  Rand(b, N, r);

  xCheck(!almostEqual(a, b, N));

  RandFree(r);
  free(a);
  free(b);
}

// Test RandNM for mean and covariance
static void test_randnm_mean_cov(void) {
  const int N = 50000;
  const int d = 3;
  double *X = NULL, *Sig = NULL, *L = NULL, *mu = NULL, *C = NULL;
  allocate(X, (size_t)N * d);
  allocate(Sig, d * d);
  allocate(L, d * d);
  allocate(mu, d);
  allocate(C, d * d);

  mu[0] = 1.0; mu[1] = -2.0; mu[2] = 0.5;
  Sig[0 + 0*d] = 1.0;
  Sig[1 + 1*d] = 2.0;
  Sig[2 + 2*d] = 0.5;

  RandRng *rng = RandCreate();
  RandSetDefault(rng);
  RandSeed(9, rng);

  double del = 0.0;
  int ok = 0;
  RandNM(mu, Sig, N, d, X, L, &del, rng, &ok);
  xCheck(ok == 1);

  // Means
  for (int j = 0; j < d; ++j) {
    double muj = mean(X + j*(size_t)N, N);
    xCheck(fabs(muj - mu[j]) < 5e-2);
  }

  // Covariance
  // Compute sample covariance
  cov("N", N, d, X, C);
  xCheck(fabs(C[0 + 0*d] - 1.0) < 1e-1);
  xCheck(fabs(C[1 + 1*d] - 2.0) < 1e-1);
  xCheck(fabs(C[2 + 2*d] - 0.5) < 1e-1);
  xCheck(fabs(C[0 + 1*d]) < 1e-1);
  xCheck(fabs(C[0 + 2*d]) < 1e-1);
  xCheck(fabs(C[1 + 2*d]) < 1e-1);
    
  RandFree(rng);
  free(X); free(Sig); free(L); free(mu); free(C);
}

void TestRandomNumbers() {
  test_uniform_basic();
  test_normal_basic();
  test_determinism_default_seed();
  test_pm_vs_default_selection();
  test_randomize_changes_stream();
  test_randnm_mean_cov();
}
