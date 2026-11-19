#include <math.h>
#include "TestUtil.h"
#include "Tests.h"
#include "error.h"

static void check_unit_spacing(void) {
  double theta = .5, cov = 5, mu = 2, t[] = {3, 4, 5, 6}, X[4], z[5], one = 1;
  varmapack_slepian_model *model = varmapack_slepian_create(&theta, &cov, &mu, 1);
  randompack_rng *rng = seededRng(242);
  xCheck(model != 0);
  if (model) {
    theta = cov = mu = 99;
    checkVarmapackSuccess(varmapack_slepian_sim(model, t, 4, 1, X, rng));
    reseedRng(rng, 242);
    ASSERT(randompack_mvn("T", 0, &one, 1, 5, z, 1, 0, rng), "RNG");
    double expected[4];
    for (int i=0; i<4; i++) expected[i] = 2+2*z[i+1]+z[i];
    checkArrayTol(X, expected, 4, 1e-13);
    t[0] = -1e100; t[1] = 1e100; t[2] = 2e100; t[3] = 3e100;
    reseedRng(rng, 242);
    checkVarmapackSuccess(varmapack_slepian_sim(model, t, 4, 1, X, rng));
    reseedRng(rng, 242);
    double independent[8];
    ASSERT(randompack_mvn("T", 0, &one, 1, 8, independent, 1, 0, rng), "RNG");
    for (int i=0; i<4; i++) expected[i] = 2+2*independent[2*i+1]+independent[2*i];
    checkArrayTol(X, expected, 4, 1e-13);
  }
  randompack_free(rng);
  varmapack_slepian_free(model);
}

static double overlap(double h) {
  return fmax(0, 1-fabs(h));
}

static void check_fractional_spacing(void) {
  int r = 2, n = 6, M = 70000;
  double theta[] = {0, .25, .5, 0}, mu[] = {2, -1};
  // Sig = I; the upper triangle of cov is ignored.
  double cov[] = {1.25, 0, NAN, 1.0625}, t[] = {0, .25, .5, 1, 1.75, 2};
  varmapack_slepian_model *model = varmapack_slepian_create(theta, cov, mu, r);
  randompack_rng *rng = seededRng(523);
  double *X = 0;
  ASSERT(ALLOC(X, r*n*M), "allocation failed");
  xCheck(model != 0);
  if (model) {
    bool ok = varmapack_slepian_sim(model, t, n, M, X, rng);
    checkVarmapackSuccess(ok);
    test_check_group means = {0}, covariances = {0};
    if (ok) {
      for (int k=0; k<n; k++) {
        for (int i=0; i<r; i++) {
          double mean = 0;
          for (int m=0; m<M; m++) mean += X[(m*n+k)*r+i];
          groupCheck(&means, fabs(mean/M-mu[i]) < .025);
          for (int l=0; l<=k; l++) {
            for (int j=0; j<r; j++) {
              double empirical = 0, h = t[k]-t[l];
              for (int m=0; m<M; m++)
                empirical += (X[(m*n+k)*r+i]-mu[i])*(X[(m*n+l)*r+j]-mu[j]);
              double expected = (i == j ? overlap(h)*cov[i+i*r] : 0)+
                overlap(h+1)*theta[j+i*r]+overlap(h-1)*theta[i+j*r];
              groupCheck(&covariances, fabs(empirical/M-expected) < .035);
            }
          }
        }
      }
      checkGroup(&means, "stationary means on uneven times");
      checkGroup(&covariances, "window overlap determines lagged covariance");
    }
  }
  FREE(X);
  randompack_free(rng);
  varmapack_slepian_free(model);
}

static void check_slepian_errors(void) {
  double theta[] = {0, 0, 0, 0}, rank_one[] = {1, 1, 1, 1};
  double incompatible[] = {0, 0, 1, 0}, target[] = {1, 0, 0, 2};
  double singular[] = {1, 0, 0, -1}, t[] = {0, .5}, X[4];
  randompack_rng *rng = seededRng(19);
  xCheck(!varmapack_slepian_create(incompatible, target, 0, 2));
  xCheck(!varmapack_slepian_create(singular, target, 0, 2));
  xCheck(!varmapack_slepian_create(0, target, 0, 2));
  varmapack_slepian_model *model = varmapack_slepian_create(theta, rank_one, 0, 2);
  xCheck(model != 0);
  if (model) {
    checkVarmapackSuccess(varmapack_slepian_sim(model, t, 2, 1, X, rng));
    xCheck(X[0] == X[1] && X[2] == X[3]);
    checkVarmapackFailure(varmapack_slepian_sim(model, t, 2, 1, X, 0));
    t[1] = t[0];
    checkVarmapackFailure(varmapack_slepian_sim(model, t, 2, 1, X, rng));
    t[1] = INFINITY;
    checkVarmapackFailure(varmapack_slepian_sim(model, t, 2, 1, X, rng));
  }
  varmapack_slepian_free(model);
  randompack_free(rng);
}

void TestSlepian(void) {
  check_unit_spacing();
  check_fractional_spacing();
  check_slepian_errors();
  double theta[] = {0, .25, .5, 0}, cov[] = {1.25, 0, NAN, 1.0625};
  double lags[] = {.5, -1.5, 2}, G[12];
  double expected[] = {.625, .125, .25, .53125, 0, .25, .125, 0, 0, 0, 0, 0};
  varmapack_slepian_model *model = varmapack_slepian_create(theta, cov, 0, 2);
  xCheck(model != 0);
  if (model) {
    checkVarmapackSuccess(varmapack_slepian_acvf(model, lags, 3, G));
    checkArrayTol(G, expected, 12, 1e-12);
    checkVarmapackSuccess(varmapack_slepian_acvf(model, 0, 0, 0));
    lags[0] = NAN;
    checkVarmapackFailure(varmapack_slepian_acvf(model, lags, 3, G));
  }
  varmapack_slepian_free(model);
}
