// Exact stationary OU transitions for diagonal drift. See Sarkka & Solin,
// Applied Stochastic Differential Equations (2019), sections 6.2 and 6.5.
#include <math.h>
#include "varmapack.h"
#include "VarmaPackUtil.h"
#include "BlasGateway.h"
#include "error.h"

struct varmapack_ou_model {
  int r;
  double *tau, *mu, *sd, *corr, *factor;
};

void varmapack_ou_free(varmapack_ou_model *model) {
  if (!model) return;
  FREE(model->tau); FREE(model->mu); FREE(model->sd);
  FREE(model->corr); FREE(model->factor);
  free(model);
}

varmapack_ou_model *varmapack_ou_create(double tau[], double cov[], double mu[],
                                       int r) {
  varmapack_ou_model *model = 0;
  double *D = 0;
  int rr;
  clear_error();
  if (!tau || !cov || r <= 0 || !intProduct(r, r, &rr) || r > INT_MAX/2) {
    fail_error("invalid OU dimensions or missing input");
    return 0;
  }
  if (!ALLOC(model, 1)) goto alloc_fail;
  model->r = r;
  if (!ALLOC(model->tau, r) || !ALLOC(model->mu, r) || !ALLOC(model->sd, r) ||
      !ALLOC(model->corr, rr) || !ALLOC(model->factor, rr) || !ALLOC(D, rr))
    goto alloc_fail;
  for (int i=0; i<r; i++) {
    if (!isfinite(tau[i]) || tau[i] <= 0 || (mu && !isfinite(mu[i])) ||
        !isfinite(cov[i+i*r]) || cov[i+i*r] < 0) {
      fail_error("OU tau must be positive and finite; mu/cov must be finite");
      goto fail;
    }
    model->tau[i] = tau[i];
    model->mu[i] = mu ? mu[i] : 0;
    model->sd[i] = sqrt(cov[i+i*r]);
  }
  for (int j=0; j<r; j++) {
    for (int i=j; i<r; i++) {
      double c = cov[i+j*r], si = model->sd[i], sj = model->sd[j];
      if (!isfinite(c) || ((si == 0 || sj == 0) && c != 0)) {
        fail_error("OU covariance is not positive semidefinite");
        goto fail;
      }
      double v = si == 0 || sj == 0 ? 0 : (c/si)/sj;
      model->corr[i+j*r] = model->corr[j+i*r] = v;
      // Correlation of the continuous-time driving noise, scaled to unit diagonal.
      double h = fabs(log(tau[i])-log(tau[j]))/2;
      double d = v == 0 ? 0 : copysign(
        exp(log(fabs(v))+h+log1p(exp(-2*h))-log(2)), v);
      D[i+j*r] = D[j+i*r] = d;
      if (!isfinite(v) || !isfinite(d)) {
        fail_error("OU covariance and persistence times are incompatible");
        goto fail;
      }
    }
  }
  if (!psdFactor(D, r, model->factor, 0)) {
    fail_error("OU covariance and persistence times are incompatible (driving PSD)");
    goto fail;
  }
  if (!psdFactor(model->corr, r, model->factor, 0)) goto fail;
  FREE(D);
  return model;
alloc_fail:
  fail_error("allocation failed");
fail:
  FREE(D);
  varmapack_ou_free(model);
  return 0;
}

bool varmapack_ou_sim(varmapack_ou_model *model, double t[], int n, int M,
                       double X[], randompack_rng *rng) {
  double *Q = 0, *F = 0, *a = 0, *z = 0;
  size_t count;
  int rn;
  clear_error();
  if (!model || !t || !X || !rng || n <= 0 || M <= 0)
    return fail_error("invalid OU simulation argument");
  int r = model->r, rr = r*r;
  if (!intProduct(r, n, &rn) || !sizeProduct(rn, M, &count) ||
      count > SIZE_MAX/sizeof(double)) return fail_error("problem size too large");
  for (int k=0; k<n; k++) {
    if (!isfinite(t[k]) || (k > 0 &&
        (!(t[k] > t[k-1]) || !isfinite(t[k]-t[k-1]))))
      return fail_error("t must be finite and strictly increasing with finite gaps");
  }
  if (!ALLOC(Q, rr) || !ALLOC(F, rr) || !ALLOC(a, r) || !ALLOC(z, (size_t)r*M))
    goto alloc_fail;
  double one = 1, previous = -1;
  for (int k=0; k<n; k++) {
    double dt = k == 0 ? 0 : t[k]-t[k-1];
    if (k > 0 && dt != previous) {
      for (int i=0; i<r; i++) a[i] = exp(-dt/model->tau[i]);
      for (int j=0; j<r; j++) {
        for (int i=0; i<r; i++) {
          double v = -expm1(-dt/model->tau[i]-dt/model->tau[j]);
          Q[i+j*r] = model->corr[i+j*r]*v;
        }
      }
      if (!psdFactor(Q, r, F, 0)) goto fail;
      previous = dt;
    }
    double *factor = k == 0 ? model->factor : F;
    if (!randompack_mvn("T", 0, &one, 1, (size_t)r*M, z, 1, 0, rng)) {
      fail_error(randompack_last_error(rng));
      goto fail;
    }
    for (int m=0; m<M; m++) {
      double *x = X+(size_t)m*rn+k*r;
      for (int i=0; i<r; i++) {
        double e = 0;
        for (int j=0; j<r; j++) e += factor[i+j*r]*z[(size_t)m*r+j];
        double old = k == 0 ? 0 : a[i]*(x[i-r]-model->mu[i]);
        x[i] = model->mu[i]+old+model->sd[i]*e;
        if (!isfinite(x[i])) {
          fail_error("OU simulated value is not representable");
          goto fail;
        }
      }
    }
  }
  FREE(z); FREE(a); FREE(F); FREE(Q);
  return true;
alloc_fail:
  fail_error("allocation failed");
fail:
  FREE(z); FREE(a); FREE(F); FREE(Q);
  return false;
}
