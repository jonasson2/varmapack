#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "error.h"
#include "TestUtil.h"
#include "Tests.h"
#include "varmapack.h"
#include "varmapack_config.h"

int TEST_NREP = TEST_NREP_DEFAULT;
#define NSERIES 11
#define PIT_BINS 50
#define MAX_Z 8
#define PIT_ALPHA 1e-12
// Upper 1e-12 quantile of chi-square with PIT_BINS - 1 degrees of freedom.
#define PIT_CHI2_LIMIT 153.7

static double series_value(double X[], int r, int t, int replicate, int i) {
  return X[replicate*r*NSERIES + t*r + i];
}

static double sample_mean(double X[], int r, int t, int i) {
  double sum = 0;
  for (int m=0; m<TEST_NREP; m++) sum += series_value(X, r, t, m, i);
  return sum/TEST_NREP;
}

static double sample_covariance(double X[], int r, int t1, int i, int t2,
                                int j, double mean1, double mean2) {
  double sum = 0;
  for (int m=0; m<TEST_NREP; m++) {
    double x1 = series_value(X, r, t1, m, i) - mean1;
    double x2 = series_value(X, r, t2, m, j) - mean2;
    sum += x1*x2;
  }
  return sum/(TEST_NREP - 1);
}

static void check_standard_error(double estimate, double target, double se,
                                 char *quantity, int t, int i, int j) {
  double difference = fabs(estimate - target);
  double z = se > 0 ? difference/se : (difference == 0 ? 0 : INFINITY);
  char message[160];
  snprintf(message, sizeof(message),
           "%s t=%d i=%d j=%d: estimate %.8g, target %.8g, z=%.3f",
           quantity, t, i, j, estimate, target, z);
  xCheckMessage(z <= MAX_Z, message);
}

static double normal_cdf(double x) {
  return 0.5*erfc(-x/sqrt(2));
}

static int compare_uniforms(const void *a, const void *b) {
  const double *x = a, *y = b;
  return (*x > *y) - (*x < *y);
}

static void check_uniform_ks(double U[], char *label) {
  double D = 0;
  double sample_size = TEST_NREP;
  // DKW bound: P(D > limit) <= PIT_ALPHA for independent uniform samples.
  double limit = sqrt(log(2/PIT_ALPHA)/(2*sample_size));
  char message[256];
  qsort(U, TEST_NREP, sizeof(*U), compare_uniforms);
  for (int m=0; m<TEST_NREP; m++) {
    D = fmax(D, (m + 1)/sample_size - U[m]);
    D = fmax(D, U[m] - m/sample_size);
  }
  snprintf(message, sizeof(message),
           "%s: KS D=%.6g, limit %.6g (alpha %.0e)", label, D, limit, PIT_ALPHA);
  xCheckMessage(D <= limit, message);
}

static void check_normal_distribution(double X[], double mean, double variance,
                                      int r, int t, int i) {
  int counts[PIT_BINS] = {0};
  double standard_deviation = sqrt(variance);
  double expected = TEST_NREP/PIT_BINS;
  double chi2 = 0;
  double *U = 0;
  char message[128];
  ASSERT(ALLOC(U, TEST_NREP), "allocation failed");
  for (int m=0; m<TEST_NREP; m++) {
    double z = (series_value(X, r, t, m, i) - mean)/standard_deviation;
    double u = normal_cdf(z);
    if (!isfinite(z) || !isfinite(u)) {
      snprintf(message, sizeof(message),
               "normal PIT t=%d i=%d replicate=%d: nonfinite value", t, i, m);
      xCheckMessage(false, message);
      FREE(U);
      return;
    }
    U[m] = u;
    int bin = (int)(PIT_BINS*u);
    if (bin < 0) bin = 0;
    if (bin >= PIT_BINS) bin = PIT_BINS - 1;
    counts[bin]++;
  }
  for (int bin=0; bin<PIT_BINS; bin++) {
    double difference = counts[bin] - expected;
    chi2 += difference*difference/expected;
  }
  snprintf(message, sizeof(message),
           "normal PIT t=%d i=%d: chi-square %.3f, limit %.1f",
           t, i, chi2, PIT_CHI2_LIMIT);
  xCheckMessage(chi2 <= PIT_CHI2_LIMIT, message);
  snprintf(message, sizeof(message), "normal PIT t=%d i=%d", t, i);
  check_uniform_ks(U, message);
  FREE(U);
}

static double target_covariance(double Gamma[], int r, int t, int i, int s,
                                 int j) {
  // Gamma_k = Cov(X_t, X_{t-k}); reverse times by transposing the matrix.
  if (t >= s) return Gamma[(t - s)*r*r + i + j*r];
  return Gamma[(s - t)*r*r + j + i*r];
}

static void check_sum(double X[], double Gamma[], double mu[], int r,
                       int times[], int components[], int count, double U[]) {
  double mean = 0, variance = 0;
  char label[160], message[256];
  for (int a=0; a<count; a++) {
    mean += mu[components[a]];
    for (int b=0; b<count; b++) {
      variance += target_covariance(Gamma, r, times[a], components[a],
                                     times[b], components[b]);
    }
  }
  snprintf(label, sizeof(label), "pair sum (%d,%d)+(%d,%d)",
           times[0], components[0], times[1], components[1]);
  if (count == 3) {
    snprintf(label, sizeof(label), "triple sum (%d,%d)+(%d,%d)+(%d,%d)",
             times[0], components[0], times[1], components[1],
             times[2], components[2]);
  }
  if (!(variance > 0) || !isfinite(variance)) {
    snprintf(message, sizeof(message), "%s: invalid target variance %.8g",
             label, variance);
    xCheckMessage(false, message);
    return;
  }
  double sd = sqrt(variance);
  for (int m=0; m<TEST_NREP; m++) {
    double sum = 0;
    for (int a=0; a<count; a++) {
      sum += series_value(X, r, times[a], m, components[a]);
    }
    double z = (sum - mean)/sd;
    if (!isfinite(z)) {
      snprintf(message, sizeof(message), "%s: nonfinite value in replicate %d",
               label, m);
      xCheckMessage(false, message);
      return;
    }
    U[m] = normal_cdf(z);
  }
  check_uniform_ks(U, label);
}

static int select_component(int r, randompack_rng *rng) {
  int i;
  bool ok = randompack_int(&i, 1, 0, r - 1, rng);
  ASSERT(ok, "component selection failed");
  return i;
}

static void check_pairs_and_triples(double X[], double Gamma[], double mu[],
                                    int r, int seed) {
  // Separate RNG keeps selections identical in default and comprehensive runs.
  randompack_rng *rng = seededRng(seed + 1000);
  int selected[3][2];
  int npairs = r == 1 ? 0 : (r == 2 ? 1 : 3);
  int lags[] = {1, 2, 10};
  int triple_times[][3] = {{0, 1, 2}, {0, 1, 10}, {0, 2, 10}};
  double *U = 0;
  ASSERT(ALLOC(U, TEST_NREP), "allocation failed");
  for (int k=0; k<npairs; k++) {
    int i, j;
    bool duplicate;
    do {
      i = select_component(r, rng);
      j = select_component(r, rng);
      duplicate = i >= j;
      for (int a=0; a<k; a++) {
        if (selected[a][0] == i && selected[a][1] == j) duplicate = true;
      }
    } while (duplicate);
    selected[k][0] = i;
    selected[k][1] = j;
    int times[] = {0, 0};
    check_sum(X, Gamma, mu, r, times, selected[k], 2, U);
  }
  for (int k=0; k<LEN(lags); k++) {
    int times[] = {0, lags[k]};
    int components[] = {select_component(r, rng), select_component(r, rng)};
    check_sum(X, Gamma, mu, r, times, components, 2, U);
  }
  for (int k=0; k<LEN(triple_times); k++) {
    int components[3];
    for (int a=0; a<3; a++) components[a] = select_component(r, rng);
    check_sum(X, Gamma, mu, r, triple_times[k], components, 3, U);
  }
  FREE(U);
  randompack_free(rng);
}

static void check_model_distribution(test_model *model, int seed) {
  int times[] = {0, 1, 2, 10};
  int lags[] = {1, 2, 10};
  int r = model->r;
  int r2 = r*r;
  double *Gamma = 0, *X = 0, *mu = 0, *means = 0;
  ASSERT(ALLOC(Gamma, r2*NSERIES), "allocation failed");
  ASSERT(ALLOC(X, r*NSERIES*TEST_NREP), "allocation failed");
  ASSERT(ALLOC(mu, r), "allocation failed");
  ASSERT(ALLOC(means, LEN(times)*r), "allocation failed");
  for (int i=0; i<r; i++) mu[i] = i + 1;
  bool ok = varmapack_acvf(model->A, model->B, model->Sig, model->p,
                           model->q, r, Gamma, NSERIES - 1);
  checkVarmapackSuccess(ok);
  randompack_rng *rng = seededRng(seed);
  ok = varmapack_sim(model->A, model->B, model->Sig, mu, 1, model->p,
                     model->q, r, NSERIES, TEST_NREP, 0, 0, 1, X, 0, rng);
  checkVarmapackSuccess(ok);
  for (int it=0; it<LEN(times); it++) {
    int t = times[it];
    for (int i=0; i<r; i++) {
      double variance = Gamma[i + i*r];
      double se = sqrt(variance/TEST_NREP);
      double estimate = sample_mean(X, r, t, i);
      means[it*r + i] = estimate;
      check_standard_error(estimate, mu[i], se, "mean", t, i, i);
      check_normal_distribution(X, mu[i], variance, r, t, i);
    }
    for (int j=0; j<r; j++) {
      for (int i=0; i<=j; i++) {
        double target = Gamma[i + j*r];
        double variance_i = Gamma[i + i*r];
        double variance_j = Gamma[j + j*r];
        double se = sqrt((variance_i*variance_j + target*target)/(TEST_NREP - 1));
        double estimate = sample_covariance(X, r, t, i, t, j,
                                            means[it*r + i], means[it*r + j]);
        check_standard_error(estimate, target, se, "covariance", t, i, j);
      }
    }
  }
  for (int ilag=0; ilag<LEN(lags); ilag++) {
    int lag = lags[ilag];
    int it = ilag + 1;
    for (int j=0; j<r; j++) {
      for (int i=0; i<r; i++) {
        double target = Gamma[lag*r2 + i + j*r];
        double variance_i = Gamma[i + i*r];
        double variance_j = Gamma[j + j*r];
        double se = sqrt((variance_i*variance_j + target*target)/(TEST_NREP - 1));
        double estimate = sample_covariance(X, r, lag, i, 0, j,
                                            means[it*r + i], means[j]);
        check_standard_error(estimate, target, se, "lagged covariance", lag,
                             i, j);
      }
    }
  }
  check_pairs_and_triples(X, Gamma, mu, r, seed);
  randompack_free(rng);
  FREE(means);
  FREE(mu);
  FREE(X);
  FREE(Gamma);
}

static void load_rho_model(test_model *model, double rho, int p, int q, int r) {
  char request[VARMAPACK_TESTCASE_NAME_LEN] = "rho";
  int index = 0;
  bool ok;
  memset(model, 0, sizeof(*model));
  model->p = p;
  model->q = q;
  model->r = r;
  ASSERT(ALLOC(model->A, r*r*(p > 0 ? p : 1)), "allocation failed");
  ASSERT(ALLOC(model->B, r*r*(q > 0 ? q : 1)), "allocation failed");
  ASSERT(ALLOC(model->Sig, r*r), "allocation failed");
  ok = varmapack_testcase(request, &index, rho, &model->p, &model->q,
                          &model->r, model->A, model->B, model->Sig, 0);
  ASSERT(ok, "rho testcase construction failed: %s", varmapack_last_error());
  snprintf(model->name, sizeof(model->name), "rho%.3f-p%d-q%d-r%d",
           rho, p, q, r);
}

static void load_near_invertible_model(test_model *model) {
  memset(model, 0, sizeof(*model));
  model->p = 1;
  model->q = 1;
  model->r = 2;
  ASSERT(ALLOC(model->A, 4), "allocation failed");
  ASSERT(ALLOC(model->B, 4), "allocation failed");
  ASSERT(ALLOC(model->Sig, 4), "allocation failed");
  model->A[0] = 0.2;
  model->A[3] = 0.1;
  model->B[0] = 0.999;
  model->B[3] = 0.999;
  model->Sig[0] = 1;
  model->Sig[3] = 1;
  STRSET(model->name, "near-invertible-r2");
  xCheck(fabs(varmapack_specrad(model->B, model->r, model->q) - 0.999) < 1e-12);
  checkVarmapackClean();
}

static void load_ill_conditioned_sigma_model(test_model *model) {
  memset(model, 0, sizeof(*model));
  model->p = 1;
  model->q = 1;
  model->r = 2;
  ASSERT(ALLOC(model->A, 4), "allocation failed");
  ASSERT(ALLOC(model->B, 4), "allocation failed");
  ASSERT(ALLOC(model->Sig, 4), "allocation failed");
  model->A[0] = 0.2;
  model->A[3] = 0.1;
  model->B[0] = 0.3;
  model->B[3] = 0.2;
  model->Sig[0] = 1;
  model->Sig[3] = 1e-10;
  STRSET(model->name, "ill-conditioned-sigma-r2");
}

void TestDistribution(void) {
  char request[VARMAPACK_TESTCASE_NAME_LEN] = "max";
  double rhos[] = {0.2, 0.8, 0.98, 0.999};
  int ncases = 0;
  bool ok = varmapack_testcase(request, &ncases, 0, 0, 0, 0, 0, 0, 0, 0);
  checkVarmapackSuccess(ok);
  for (int icase=1; icase<=ncases; icase++) {
    test_model model;
    loadIndexedTestModel(&model, icase);
    xCheckAddMsg(model.name);
    check_model_distribution(&model, 100 + icase);
    freeTestModel(&model);
  }
  for (int i=0; i<LEN(rhos); i++) {
    test_model model;
    load_rho_model(&model, rhos[i], 2, 2, 4);
    xCheckAddMsg(model.name);
    check_model_distribution(&model, 200 + i);
    freeTestModel(&model);
  }
  test_model model;
  load_rho_model(&model, 0.95, 3, 3, 20);
  xCheckAddMsg(model.name);
  check_model_distribution(&model, 300);
  freeTestModel(&model);
  load_near_invertible_model(&model);
  xCheckAddMsg(model.name);
  check_model_distribution(&model, 400);
  freeTestModel(&model);
  load_ill_conditioned_sigma_model(&model);
  xCheckAddMsg(model.name);
  check_model_distribution(&model, 500);
  freeTestModel(&model);
}
