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
    xCheck(fabs(varmapack_ou_specrad(model, 1.5)-exp(-.5)) < 1e-14);
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

static void check_ou_matrix_transitions(void) {
  double K[] = {1, 0, -1, 1}, C[] = {1, 0, NAN, 1}, t[] = {0, .5, 1, 1001};
  double X[8], normal[8], one = 1;
  randompack_rng *rng = seededRng(81);
  varmapack_ou_model *model = varmapack_ou_create_drift(K, C, 0, 2);
  xCheck(model != 0);
  if (model) {
    xCheck(fabs(varmapack_ou_specrad(model, .5)-exp(-.5)) < 1e-14);
    checkVarmapackNaN(varmapack_ou_specrad(model, 0));
    K[2] = 99; // The native constructor copies its inputs.
    checkVarmapackSuccess(varmapack_ou_sim(model, t, 4, 1, X, rng));
    reseedRng(rng, 81);
    for (int i=0; i<8; i+=2)
      ASSERT(randompack_mvn("T", 0, &one, 1, 2, normal+i, 1, 0, rng), "RNG");
    double expected[8] = {normal[0], normal[1]};
    for (int k=1; k<4; k++) {
      double dt = t[k]-t[k-1], a = exp(-dt), a2 = a*a;
      // A = a*[1 dt; 0 1]; Q = I - A*A'.
      double l11 = sqrt(1-a2*(1+dt*dt)), l21 = -a2*dt/l11;
      double l22 = sqrt(1-a2-l21*l21);
      int i = 2*k;
      expected[i] = a*(expected[i-2]+dt*expected[i-1])+l11*normal[i];
      expected[i+1] = a*expected[i-1]+l21*normal[i]+l22*normal[i+1];
    }
    checkArrayTol(X, expected, 8, 1e-12);
    t[1] = 1e-18;
    reseedRng(rng, 81);
    checkVarmapackSuccess(varmapack_ou_sim(model, t, 2, 1, X, rng));
    // The increment covariance divided by the tiny gap must approach D.
    double increment[] = {(X[2]-X[0])/1e-9, (X[3]-X[1])/1e-9};
    double limit[] = {sqrt(2)*normal[2], -normal[2]/sqrt(2)+sqrt(1.5)*normal[3]};
    checkArrayTol(increment, limit, 2, 1e-6);
  }
  varmapack_ou_free(model);
  double unstable[] = {0, 1, -1, 0}, incompatible[] = {1, 0, 3, 1};
  double stiff[] = {1, 0, 0, 1e15}, tau[] = {1, 1e-15}, reference[8];
  t[1] = .5;
  model = varmapack_ou_create_drift(stiff, C, 0, 2);
  varmapack_ou_model *diag = varmapack_ou_create(tau, C, 0, 2);
  reseedRng(rng, 42);
  checkVarmapackSuccess(varmapack_ou_sim(diag, t, 4, 1, reference, rng));
  reseedRng(rng, 42);
  checkVarmapackSuccess(varmapack_ou_sim(model, t, 4, 1, X, rng));
  checkArrayTol(X, reference, 8, 1e-13);
  varmapack_ou_free(model); varmapack_ou_free(diag);
  stiff[2] = 1;
  model = varmapack_ou_create_drift(stiff, C, 0, 2);
  xCheck(model != 0);
  if (model) {
    reseedRng(rng, 42);
    // The unit off-diagonal changes this transition by only O(1e-15).
    if (varmapack_ou_sim(model, t, 2, 1, X, rng))
      checkArrayTol(X, reference, 4, 1e-9);
    else xCheck(strstr(varmapack_last_error(), "accuracy") != 0);
  }
  varmapack_ou_free(model);
  xCheck(!varmapack_ou_create_drift(unstable, C, 0, 2));
  xCheck(!varmapack_ou_create_drift(incompatible, C, 0, 2));
  double supported[] = {1, 0, 1, 2}, rank_one[] = {1, 1, NAN, 1};
  model = varmapack_ou_create_drift(supported, rank_one, 0, 2);
  xCheck(model != 0);
  if (model) {
    checkVarmapackSuccess(varmapack_ou_sim(model, t, 4, 1, X, rng));
    test_check_group support = {0};
    for (int i=0; i<8; i+=2) groupCheck(&support, fabs(X[i]-X[i+1]) < 1e-7);
    checkGroup(&support, "matrix OU preserves rank-one stationary support");
  }
  varmapack_ou_free(model);
  randompack_free(rng);
}

void TestOU(void) {
  check_ou_statistics();
  check_ou_errors_and_psd();
  check_ou_transitions();
  check_ou_matrix_transitions();
  double K[] = {1, 0, -1, 1}, cov[] = {1, 0, NAN, 1};
  double lags[] = {.5, -.5, 0}, G[12], a = exp(-.5);
  double expected[] = {a, 0, .5*a, a, a, .5*a, 0, a, 1, 0, 0, 1};
  varmapack_ou_model *model = varmapack_ou_create_drift(K, cov, 0, 2);
  xCheck(model != 0);
  if (model) {
    checkVarmapackSuccess(varmapack_ou_acvf(model, lags, 3, G));
    checkArrayTol(G, expected, 12, 1e-12);
    checkVarmapackSuccess(varmapack_ou_acvf(model, 0, 0, 0));
    lags[0] = INFINITY;
    checkVarmapackFailure(varmapack_ou_acvf(model, lags, 3, G));
  }
  varmapack_ou_free(model);
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
