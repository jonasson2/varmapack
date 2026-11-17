#include <math.h>
#include "TestUtil.h"
#include "Tests.h"
#include "error.h"

static void check_ou_statistics(void) {
  double tau[] = {1, 2, 3}, cov[] = {1, .2, -.1, .2, 2, .3, -.1, .3, 1};
  double mu[] = {2, -1, .5}, t[] = {-3, -2.9, -1, 3};
  int r = 3, n = 4, M = 100000;
  double *X = 0;
  ASSERT(ALLOC(X, r*n*M), "allocation failed");
  randompack_rng *rng = seededRng(132);
  varmapack_ou_model *model = varmapack_ou_create(tau, cov, mu, r);
  xCheck(model != 0);
  if (model) {
    bool ok = varmapack_ou_sim(model, t, n, M, X, rng);
    checkVarmapackSuccess(ok);
    test_check_group means = {0}, covariances = {0};
    if (ok) {
      for (int k=0; k<n; k++) {
        for (int i=0; i<r; i++) {
          double sum = 0;
          for (int m=0; m<M; m++) sum += X[m*r*n+k*r+i];
          groupCheck(&means, fabs(sum/M-mu[i]) < .025);
          for (int h=0; h<=k; h++) {
            for (int j=0; j<r; j++) {
              double v = 0;
              for (int m=0; m<M; m++) {
                v += (X[m*r*n+k*r+i]-mu[i])*(X[m*r*n+h*r+j]-mu[j]);
              }
              double expected = cov[i+j*r]*exp(-(t[k]-t[h])/tau[i]);
              groupCheck(&covariances, fabs(v/M-expected) < .035);
            }
          }
        }
      }
      checkGroup(&means, "stationary means at all observation times");
      checkGroup(&covariances, "stationary and lagged cross-covariances");
    }
  }
  varmapack_ou_free(model);
  randompack_free(rng);
  FREE(X);
}

static void check_ou_errors_and_psd(void) {
  double tau[] = {1, 2}, cov[] = {1, 1, 1, 1}, t[] = {0, 1}, X[4];
  randompack_rng *rng = seededRng(72);
  xCheck(!varmapack_ou_create(tau, cov, 0, 2)); // Valid covariance, invalid OU.
  tau[1] = 1;
  varmapack_ou_model *model = varmapack_ou_create(tau, cov, 0, 2);
  xCheck(model != 0);
  if (model) {
    checkVarmapackSuccess(varmapack_ou_sim(model, t, 2, 1, X, rng));
    xCheck(fabs(X[0]-X[1]) < 1e-12 && fabs(X[2]-X[3]) < 1e-12);
    t[1] = 0;
    checkVarmapackFailure(varmapack_ou_sim(model, t, 2, 1, X, rng));
    t[1] = NAN;
    checkVarmapackFailure(varmapack_ou_sim(model, t, 2, 1, X, rng));
    t[1] = 1;
    checkVarmapackFailure(varmapack_ou_sim(model, t, 2, 0, X, rng));
    checkVarmapackFailure(varmapack_ou_sim(model, t, 2, 1, X, 0));
    checkVarmapackSuccess(varmapack_ou_sim(model, t, 1, 1, X, rng));
  }
  varmapack_ou_free(model);
  varmapack_ou_free(0);
  cov[1] = 0; cov[2] = NAN;
  model = varmapack_ou_create(tau, cov, 0, 2);
  xCheck(model != 0); // Upper triangle is ignored.
  varmapack_ou_free(model);
  tau[0] = 0;
  xCheck(!varmapack_ou_create(tau, cov, 0, 2));
  tau[0] = INFINITY;
  xCheck(!varmapack_ou_create(tau, cov, 0, 2));
  xCheck(!varmapack_ou_create(tau, cov, 0, INT_MAX));
  xCheck(!varmapack_ou_create(0, cov, 0, 2));
  tau[0] = 1;
  cov[0] = cov[3] = 0;
  model = varmapack_ou_create(tau, cov, 0, 2);
  xCheck(model != 0);
  if (model) {
    checkVarmapackSuccess(varmapack_ou_sim(model, t, 2, 1, X, rng));
    checkArrayZero(X, 4);
  }
  varmapack_ou_free(model);
  randompack_free(rng);
}

static void check_ou_transitions(void) {
  double tau = 2, cov = 4, mu = 3, times[] = {0, .25, .5, 100};
  double X[4], normal[4], one = 1;
  randompack_rng *rng = seededRng(31);
  varmapack_ou_model *model = varmapack_ou_create(&tau, &cov, &mu, 1);
  xCheck(model != 0);
  if (model) {
    // Native models own their parameters, not pointers to caller storage.
    tau = cov = mu = 99;
    checkVarmapackSuccess(varmapack_ou_sim(model, times, 4, 1, X, rng));
    reseedRng(rng, 31);
    for (int i=0; i<4; i++)
      ASSERT(randompack_mvn("T", 0, &one, 1, 1, normal+i, 1, 0, rng), "RNG");
    double expected[4] = {3+2*normal[0]};
    for (int i=1; i<4; i++) {
      double gap = times[i]-times[i-1];
      expected[i] = 3+exp(-gap/2)*(expected[i-1]-3)+
        2*sqrt(-expm1(-gap))*normal[i];
    }
    checkArrayTol(X, expected, 4, 1e-13);
  }
  varmapack_ou_free(model);
  double taus[] = {1e-200, 1e200}, C[] = {1e-200, 0, NAN, 1e200};
  model = varmapack_ou_create(taus, C, 0, 2);
  xCheck(model != 0);
  if (model) {
    double t[] = {0, 1e-200, 1e200}, Y[6];
    checkVarmapackSuccess(varmapack_ou_sim(model, t, 3, 1, Y, rng));
    checkArrayFinite(Y, 6);
  }
  varmapack_ou_free(model);
  randompack_free(rng);
}

void TestOU(void) {
  check_ou_statistics();
  check_ou_errors_and_psd();
  check_ou_transitions();
}

void TestGammaTimes(void) {
  int n = 100001;
  double *t = 0;
  ASSERT(ALLOC(t, n), "allocation failed");
  randompack_rng *rng = seededRng(14);
  double shapes[] = {.5, 1, 4, 16, 100};
  for (int k=0; k<LEN(shapes); k++) {
    checkVarmapackSuccess(varmapack_gamma_times(n, shapes[k], 2, t, rng));
    double sum = 0, squares = 0;
    test_check_group increasing = {0};
    for (int i=1; i<n; i++) {
      double d = t[i]-t[i-1];
      sum += d;
      squares += d*d;
      groupCheck(&increasing, d > 0 && isfinite(t[i]));
    }
    checkGroup(&increasing, "positive Gamma gaps");
    double mean = sum/(n-1), variance = squares/(n-1)-mean*mean;
    xCheck(fabs(mean-2) < 7*sqrt(4/shapes[k]/(n-1)));
    xCheck(fabs(variance-4/shapes[k]) < .07*4/shapes[k]);
  }
  checkVarmapackSuccess(varmapack_gamma_times(4, INFINITY, 2, t, 0));
  double expected[] = {0, 2, 4, 6};
  checkArraySame(t, expected, 4);
  checkVarmapackSuccess(varmapack_gamma_times(1, 4, 2, t, 0));
  xCheck(t[0] == 0);
  checkVarmapackFailure(varmapack_gamma_times(0, 4, 2, t, rng));
  checkVarmapackFailure(varmapack_gamma_times(n, 0, 2, t, rng));
  checkVarmapackFailure(varmapack_gamma_times(n, NAN, 2, t, rng));
  checkVarmapackFailure(varmapack_gamma_times(n, 4, 0, t, rng));
  checkVarmapackFailure(varmapack_gamma_times(n, 4, 2, t, 0));
  checkVarmapackFailure(varmapack_gamma_times(4, INFINITY, 1e308, t, 0));
  randompack_free(rng);
  FREE(t);
}
