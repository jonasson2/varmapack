#include <math.h>
#include "slepian_internal.h"
#include "BlasGateway.h"
#include "VarmaPackUtil.h"
#include "error.h"

void varmapack_slepian_free(varmapack_slepian_model *model) {
  if (!model) return;
  FREE(model->theta);
  FREE(model->mu);
  FREE(model->factor);
  FREE(model);
}

varmapack_slepian_model *varmapack_slepian_create(double theta[], double cov[],
                                                 double mu[], int r) {
  varmapack_slepian_model *model = 0;
  double *Sig = 0;
  int rr;
  clear_error();
  if (!theta || !cov || r < 1 || !intProduct(r, r, &rr)) {
    fail_error("invalid Slepian dimensions or missing input");
    return 0;
  }
  if (!ALLOC(model, 1) || !ALLOC(model->theta, rr) || !ALLOC(model->mu, r) ||
      !ALLOC(model->factor, rr) || !ALLOC(Sig, rr)) goto alloc_fail;
  model->r = r;
  for (int i=0; i<r; i++) {
    if (mu && !isfinite(mu[i])) {
      fail_error("Slepian mu must be finite");
      goto fail;
    }
    model->mu[i] = mu ? mu[i] : 0;
  }
  if (!varmapack_cov2sig(0, theta, cov, 0, 1, r, Sig)) goto fail;
  // Factor in marginal units so rank truncation cannot discard a small series.
  double *sd = model->theta; // Scratch until theta is copied below.
  for (int i=0; i<r; i++) sd[i] = sqrt(Sig[i+i*r]);
  for (int j=0; j<r; j++) {
    for (int i=j; i<r; i++) {
      double v = Sig[i+j*r];
      if (sd[i] && sd[j]) v = (v/sd[i])/sd[j];
      Sig[i+j*r] = Sig[j+i*r] = v;
    }
  }
  if (!psdFactor(Sig, r, model->factor, 0)) goto fail;
  for (int i=0; i<r; i++) scal(r, sd[i], model->factor+i, r);
  copy(rr, theta, 1, model->theta, 1);
  FREE(Sig);
  return model;
alloc_fail:
  fail_error("allocation failed");
fail:
  FREE(Sig);
  varmapack_slepian_free(model);
  return 0;
}

typedef struct {
  double hi, lo;
} slepian_endpoint;

static slepian_endpoint endpoint(double t, double shift) {
  // TwoSum preserves shifts even when the absolute time is large.
  double hi = t+shift, z = hi-t;
  slepian_endpoint value = {hi, (t-(hi-z))+(shift-z)};
  return value;
}

static bool earlier(slepian_endpoint a, slepian_endpoint b) {
  return a.hi < b.hi || (a.hi == b.hi && a.lo < b.lo);
}

static int merge_windows(double t[], int n, int index[], double scale[]) {
  int pos[3] = {0}, active = 0, count = 0;
  slepian_endpoint next[3] = {endpoint(t[0], -2), endpoint(t[0], -1),
                            endpoint(t[0], 0)};
  slepian_endpoint prev = next[0];
  // Record endpoints as offsets into the increments; omit gaps outside all windows.
  for (int k=0; k<3*n; k++) {
    int j = -1;
    for (int h=0; h<3; h++)
      if (pos[h] < n && (j < 0 || earlier(next[h], next[j]))) j = h;
    slepian_endpoint here = next[j];
    if (active && earlier(prev, here)) {
      double gap = (here.hi-prev.hi)+(here.lo-prev.lo);
      scale[count++] = sqrt(gap);
    }
    index[j*n+pos[j]] = count;
    if (j == 0) active++;
    if (j == 2) active--;
    prev = here;
    pos[j]++;
    if (pos[j] < n) next[j] = endpoint(t[pos[j]], j-2);
  }
  return count;
}

static void slide_window(int left, int right, int *first, int *last,
                          double increments[], int r, double sum[]) {
  // Accumulate only one window, not the entire Brownian path.
  if (left >= *last) {
    for (int i=0; i<r; i++) sum[i] = 0;
    *last = left;
  }
  else {
    for (int k=*first; k<left; k++)
      for (int i=0; i<r; i++) sum[i] -= increments[(size_t)k*r+i];
  }
  for (int k=*last; k<right; k++)
    for (int i=0; i<r; i++) sum[i] += increments[(size_t)k*r+i];
  *first = left;
  *last = right;
}

bool varmapack_slepian_sim(varmapack_slepian_model *model, double t[], int n, int M,
                            double X[], randompack_rng *rng) {
  double *scale = 0, *z = 0, *increments = 0, *sums = 0;
  int *index = 0;
  int rn, nz;
  size_t count;
  bool ok = false;
  clear_error();
  if (!model || !t || !X || !rng || n < 1 || M < 1)
    return fail_error("invalid Slepian simulation argument");
  int r = model->r;
  if (n > INT_MAX/3 || !intProduct(r, n, &rn) || !sizeProduct(rn, M, &count) ||
      count > SIZE_MAX/sizeof(double)) return fail_error("problem size too large");
  for (int i=0; i<n; i++) {
    if (!isfinite(t[i]) || (i > 0 && !(t[i] > t[i-1])))
      return fail_error("t must be finite and strictly increasing");
  }
  if (!ALLOC(index, 3*n) || !ALLOC(scale, 3*n-1)) goto alloc_fail;
  int ne = merge_windows(t, n, index, scale);
  if (!intProduct(r, ne, &nz)) {
    fail_error("problem size too large");
    goto done;
  }
  if (!ALLOC(z, nz) || !ALLOC(increments, nz) || !ALLOC(sums, 2*r)) goto alloc_fail;
  double one = 1;
  for (int m=0; m<M; m++) {
    if (!randompack_mvn("T", 0, &one, 1, nz, z, 1, 0, rng)) {
      fail_error(randompack_last_error(rng));
      goto done;
    }
    gemm("N", "N", r, ne, r, 1, model->factor, r, z, r, 0, increments, r);
    for (int k=0; k<ne; k++)
      for (int i=0; i<r; i++) increments[(size_t)k*r+i] *= scale[k];
    int a = 0, b = 0, c = 0, d = 0;
    for (int k=0; k<n; k++) {
      int left = index[k], mid = index[n+k], right = index[2*n+k];
      slide_window(left, mid, &a, &b, increments, r, sums);
      slide_window(mid, right, &c, &d, increments, r, sums+r);
      double *x = X+(size_t)m*rn+k*r;
      for (int i=0; i<r; i++) {
        double value = model->mu[i]+sums[r+i];
        for (int j=0; j<r; j++) value += model->theta[i+j*r]*sums[j];
        if (!isfinite(value)) {
          fail_error("Slepian simulated value is not representable");
          goto done;
        }
        x[i] = value;
      }
    }
  }
  ok = true;
  goto done;
alloc_fail:
  fail_error("allocation failed");
done:
  FREE(sums); FREE(increments); FREE(z); FREE(scale); FREE(index);
  return ok;
}
