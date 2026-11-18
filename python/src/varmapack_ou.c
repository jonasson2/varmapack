// Exact stationary OU transitions. See Sarkka & Solin,
// Applied Stochastic Differential Equations (2019), sections 6.2 and 6.5.
#include <math.h>
#include "varmapack.h"
#include "VarmaPackUtil.h"
#include "BlasGateway.h"
#include "mb05odGateway.h"
#include "error.h"

struct varmapack_ou_model {
  int r;
  double *tau, *mu, *sd, *corr, *factor;
  double *drift, *driving, *cov;
  double drift_scale, driving_scale;
};

void varmapack_ou_free(varmapack_ou_model *model) {
  if (!model) return;
  FREE(model->tau); FREE(model->mu); FREE(model->sd);
  FREE(model->corr); FREE(model->factor);
  FREE(model->drift); FREE(model->driving); FREE(model->cov);
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

static bool ou_stable_drift(double K[], int r, double scratch[]) {
  double *wr = 0, *wi = 0, *work = 0, query, dummy[1];
  int info, lwork;
  bool ok = false;
  if (!ALLOC(wr, r) || !ALLOC(wi, r)) goto alloc_fail;
  copy(r*r, K, 1, scratch, 1);
  geev("N", "N", r, scratch, r, wr, wi, dummy, 1, dummy, 1, &query, -1, &info);
  if (info != 0 || !isfinite(query) || query < 1 || query > INT_MAX) goto eig_fail;
  lwork = ceil(query);
  if (!ALLOC(work, lwork)) goto alloc_fail;
  geev("N", "N", r, scratch, r, wr, wi, dummy, 1, dummy, 1, work, lwork, &info);
  if (info != 0) goto eig_fail;
  for (int i=0; i<r; i++) {
    if (!isfinite(wr[i]) || !isfinite(wi[i]) || wr[i] <= 0) {
      fail_error("OU drift matrix must have eigenvalues with positive real parts");
      goto done;
    }
  }
  ok = true;
  goto done;
alloc_fail:
  fail_error("allocation failed");
  goto done;
eig_fail:
  fail_error("OU drift eigenvalue computation failed");
done:
  FREE(work); FREE(wi); FREE(wr);
  return ok;
}

varmapack_ou_model *varmapack_ou_create_drift(double K[], double cov[], double mu[],
                                             int r) {
  varmapack_ou_model *model = 0;
  double *C = 0, *KC = 0;
  int rr, block_size;
  clear_error();
  if (!K || !cov || r <= 0 || !intProduct(r, r, &rr) || r > INT_MAX/2 ||
      !intProduct(2*r, 2*r, &block_size) ||
      2*(size_t)block_size + 20*(size_t)r + 9 > INT_MAX) {
    fail_error("invalid OU dimensions or missing input");
    return 0;
  }
  bool diagonal = true;
  for (int j=0; j<r; j++) {
    if (!(K[j+j*r] > 0) || !isfinite(K[j+j*r]) || !isfinite(1/K[j+j*r]))
      diagonal = false;
    for (int i=0; i<r; i++) if (i != j && K[i+j*r] != 0) diagonal = false;
  }
  if (diagonal) {
    double *tau = 0;
    if (!ALLOC(tau, r)) goto alloc_fail;
    for (int i=0; i<r; i++) tau[i] = 1/K[i+i*r];
    model = varmapack_ou_create(tau, cov, mu, r);
    FREE(tau);
    return model;
  }
  if (!ALLOC(model, 1)) goto alloc_fail;
  model->r = r;
  if (!ALLOC(model->drift, rr) || !ALLOC(model->driving, rr) ||
      !ALLOC(model->factor, rr) || !ALLOC(model->mu, r) ||
      !ALLOC(C, rr) || !ALLOC(KC, rr)) goto alloc_fail;
  for (int i=0; i<rr; i++) {
    if (!isfinite(K[i])) {
      fail_error("OU drift matrix must be finite");
      goto fail;
    }
    model->drift_scale = fmax(model->drift_scale, fabs(K[i]));
    model->drift[i] = K[i];
  }
  for (int j=0; j<r; j++) {
    if (mu && !isfinite(mu[j])) {
      fail_error("OU mu must be finite");
      goto fail;
    }
    model->mu[j] = mu ? mu[j] : 0;
    for (int i=j; i<r; i++) {
      double v = cov[i+j*r];
      if (!isfinite(v) || (i == j && v < 0)) {
        fail_error("OU covariance must be finite and positive semidefinite");
        goto fail;
      }
      C[i+j*r] = C[j+i*r] = v;
    }
  }
  if (!ou_stable_drift(model->drift, r, KC)) goto fail;
  if (!psdFactor(C, r, model->factor, 0)) goto fail;
  gemm("N", "N", r, r, r, 1, model->drift, r, C, r, 0, KC, r);
  for (int j=0; j<r; j++) {
    for (int i=j; i<r; i++) {
      double v = KC[i+j*r]+KC[j+i*r];
      if (!isfinite(v)) {
        fail_error("OU driving covariance is not representable");
        goto fail;
      }
      model->driving[i+j*r] = model->driving[j+i*r] = v;
      model->driving_scale = fmax(model->driving_scale, fabs(v));
    }
  }
  if (model->driving_scale > 0) {
    for (int i=0; i<rr; i++) model->driving[i] /= model->driving_scale;
  }
  if (!psdFactor(model->driving, r, KC, 0)) {
    fail_error("OU covariance and drift matrix are incompatible (driving PSD)");
    goto fail;
  }
  model->cov = C;
  FREE(KC);
  return model;
alloc_fail:
  fail_error("allocation failed");
fail:
  FREE(KC); FREE(C);
  varmapack_ou_free(model);
  return 0;
}

static bool ou_matrix_transition(varmapack_ou_model *model, double dt,
                                 double A[], double Q[], double B[], double work[]) {
  int r = model->r, rr = r*r, n = 2*r, steps = 0;
  double h = dt;
  // Start with a short Van Loan block exponential, then double the transition.
  // This avoids both subtracting nearly equal covariances and exponentiating an
  // unstable block over a long gap. Keep Q in units of h * driving_scale.
  while (model->drift_scale*h > .5/r) {
    h *= .5;
    steps++;
  }
  memset(B, 0, n*n*sizeof(*B));
  for (int j=0; j<r; j++) {
    for (int i=0; i<r; i++) {
      B[i+j*n] = -h*model->drift[i+j*r];
      B[r+i+(r+j)*n] = h*model->drift[j+i*r];
      B[i+(r+j)*n] = model->driving[i+j*r]/2;
    }
  }
  mb05od_accuracy accuracy;
  if (!mb05odExponential(B, 1, n, B, false, &accuracy)) return false;
  if (accuracy.iwarn == 1 || accuracy.iwarn == 2)
    return fail_error("OU matrix exponential may be inaccurate");
  lacpy("All", r, r, B, n, A, r);
  gemm("N", "T", r, r, r, 2, B+r*n, n, A, r, 0, Q, r);
  for (int k=0; k<steps; k++) {
    // Average rather than sum so Q stays normalized by the current gap.
    gemm("N", "N", r, r, r, 1, A, r, Q, r, 0, work, r);
    gemm("N", "T", r, r, r, .5, work, r, A, r, .5, Q, r);
    gemm("N", "N", r, r, r, 1, A, r, A, r, 0, work, r);
    copy(rr, work, 1, A, 1);
  }
  int gap_exp, driving_exp;
  double scale = frexp(dt, &gap_exp)*frexp(model->driving_scale, &driving_exp);
  for (int j=0; j<r; j++) {
    for (int i=j; i<r; i++) {
      int q_exp;
      double q = frexp(Q[i+j*r]/2+Q[j+i*r]/2, &q_exp);
      double v = scalbn(q*scale, gap_exp+driving_exp+q_exp);
      if (!isfinite(v)) return fail_error("OU innovation covariance is not representable");
      Q[i+j*r] = Q[j+i*r] = v;
    }
  }
  for (int i=0; i<rr; i++) {
    if (!isfinite(A[i])) return fail_error("OU transition matrix is not representable");
  }
  // Squaring a short transition can erase slow modes. Check stationarity
  // in marginal SD units, independently of the Van Loan calculation.
  gemm("N", "N", r, r, r, 1, A, r, model->cov, r, 0, work, r);
  gemm("N", "T", r, r, r, 1, work, r, A, r, 0, B, r);
  for (int j=0; j<r; j++) {
    for (int i=j; i<r; i++) {
      double residual = (B[i+j*r]-model->cov[i+j*r])+Q[i+j*r];
      double tol = (1e-10*r*sqrt(model->cov[i+i*r]))*sqrt(model->cov[j+j*r]);
      if (!isfinite(residual) || fabs(residual) > tol)
        return fail_error("OU transition failed stationarity accuracy check");
    }
  }
  return true;
}

bool varmapack_ou_sim(varmapack_ou_model *model, double t[], int n, int M,
                       double X[], randompack_rng *rng) {
  double *Q = 0, *F = 0, *a = 0, *z = 0, *block = 0, *work = 0;
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
  bool matrix = model->drift != 0;
  if (!ALLOC(Q, rr) || !ALLOC(F, rr) || !ALLOC(a, matrix ? rr : r) ||
      !ALLOC(z, (size_t)r*M) ||
      (matrix && (!ALLOC(block, 4*(size_t)rr) || !ALLOC(work, rr))))
    goto alloc_fail;
  double one = 1, previous = -1;
  for (int k=0; k<n; k++) {
    double dt = k == 0 ? 0 : t[k]-t[k-1];
    if (k > 0 && dt != previous) {
      if (matrix) {
        if (!ou_matrix_transition(model, dt, a, Q, block, work)) goto fail;
      }
      else {
        for (int i=0; i<r; i++) a[i] = exp(-dt/model->tau[i]);
        for (int j=0; j<r; j++) {
          for (int i=0; i<r; i++) {
            double v = -expm1(-dt/model->tau[i]-dt/model->tau[j]);
            Q[i+j*r] = model->corr[i+j*r]*v;
          }
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
        double old = 0;
        if (k > 0) {
          if (matrix) {
            for (int j=0; j<r; j++) old += a[i+j*r]*(x[j-r]-model->mu[j]);
          }
          else old = a[i]*(x[i-r]-model->mu[i]);
        }
        x[i] = model->mu[i]+old+(matrix ? e : model->sd[i]*e);
        if (!isfinite(x[i])) {
          fail_error("OU simulated value is not representable");
          goto fail;
        }
      }
    }
  }
  FREE(z); FREE(a); FREE(F); FREE(Q);
  FREE(work); FREE(block);
  return true;
alloc_fail:
  fail_error("allocation failed");
fail:
  FREE(z); FREE(a); FREE(F); FREE(Q);
  FREE(work); FREE(block);
  return false;
}
