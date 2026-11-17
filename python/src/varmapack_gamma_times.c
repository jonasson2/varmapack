#include <math.h>
#include "varmapack.h"
#include "error.h"

// Marsaglia & Tsang (2000), doi:10.1145/358407.358414. Use the uniform
// gateway shared by C, Python and R, without requiring a newer Randompack ABI.
static bool gamma_gap(double shape, double mean, double *gap, randompack_rng *rng) {
  double d = (shape < 1 ? shape+1 : shape)-1.0/3;
  double c = (1.0/3)/sqrt(d);
  for (int k=0; k<100000; k++) {
    double u[4];
    if (!randompack_u01(u, 4, rng)) return fail_error(randompack_last_error(rng));
    double z = sqrt(-2*log1p(-u[0]))*cos(6.2831853071795864769*u[1]);
    double w = c*z;
    if (w <= -1) continue;
    double v = (1+w)*(1+w)*(1+w);
    if (u[2] > 1-0.0331*z*z*z*z &&
        log(u[2]) > z*z/2+d*(3*log1p(w)-(v-1))) continue;
    // Log scale for shape < 1 avoids an overflowing d/shape intermediate.
    *gap = shape < 1 ? exp(log(mean)+log(d)-log(shape)+3*log1p(w)+
                            log1p(-u[3])/shape) : mean*(v*(d/shape));
    if (!isfinite(*gap) || *gap <= 0)
      return fail_error("Gamma gap is not representable as a positive finite double");
    return true;
  }
  return fail_error("Gamma rejection sampler did not converge");
}

bool varmapack_gamma_times(int n, double shape, double mean_gap, double t[],
                            randompack_rng *rng) {
  clear_error();
  if (n <= 0 || !t || !(shape > 0) || !isfinite(mean_gap) || mean_gap <= 0 ||
      (n > 1 && isfinite(shape) && !rng))
    return fail_error("invalid Gamma times argument");
  t[0] = 0;
  for (int i=1; i<n; i++) {
    double gap = mean_gap;
    if (isfinite(shape) && !gamma_gap(shape, mean_gap, &gap, rng)) return false;
    t[i] = t[i-1]+gap;
    if (!isfinite(t[i]) || !(t[i] > t[i-1]))
      return fail_error("Gamma observation time is not representable");
  }
  return true;
}
