// RandomNumbers — random number generation utilities for VARMASIM
//
// See RandomNumbers.h for further information, including parameter descriptions

#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <stdbool.h>
#include "allocate.h"
#include "RandomNumbers.h"
#include "BlasGateway.h"
#include "xAssert.h"
#include "random.h"

#ifdef USING_R
#include <R.h>
#include <Rmath.h>
#endif

// TODO  Matlab mex

static uint64_t splitmix64(uint64_t *x) {
    uint64_t z = (*x += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

RandRng *RandCreate(void) {
    return rand_create();
}

void RandFree(RandRng *rng) {
    rand_free(rng);
}

void RandSetPM(RandRng *rng) {
    rand_settype(PARKMILLER, rng);
}

void RandSetDefault(RandRng *rng) {
    rand_settype(DEFAULT_RNG, rng);
}

void RandSeed(int seed, RandRng *rng) {
    uint32_t pm = (uint32_t)seed;
    uint64_t x = (uint64_t)(uint32_t)seed;
    uint64_t s0 = splitmix64(&x);
    uint64_t s1 = splitmix64(&x);
    if (s0 == 0 && s1 == 0) s1 = 1;
    rand_setPMseed(pm, rng);
    rand_setstate(s0, s1, rng);
}

void RandRandomize(RandRng *rng) {
  rand_randomize(0, rng);
}

void RandomThreadRandomize(uint64_t thread_id, RandRng *rng) {
  rand_randomize(thread_id, rng);
}

void Rand (double x[], int n, RandRng *rng) {
#ifdef USING_R
  if (rng->type == PARKMILLER)
    rand_dble(x, n, rng);
  else {
    GetRNGstate();
    for (int i = 0; i < n; i++) {
      x[i] = unif_rand(); // R-s built-in generator
    }
    PutRNGstate();
  }
#else
  rand_dble(x, n, rng);
#endif  
}

void RandN (double x[], int n, RandRng *rng) {
#ifdef USING_R
  if (rng->type == PARKMILLER)
    rand_normal(x, n, rng);
  else {
    GetRNGstate();
    for (int i = 0; i < n; i++) {
      x[i] = unif_rand(); // R-s built-in generator
    }
    PutRNGstate();
  }
#else
  rand_normal(x, n, rng);
#endif  
}

void RandNM (double mu[], double Sig[], int n, int d, double X[], double L[],
	     double *del, RandRng *rng, int *ok) {
  int imx, info, i;
  double *LSig, del0, macheps, fourthirds, onethird;
  double maxeps[] = {1.e-5, 2.e-14, 2.e-16, 1.e-17, 2.e-32};
  double mineps[] = {1.e-8, 2.e-17, 2.e-19, 1.e-20, 2.e-35};
  // Determine machine epsilon with the method of [3]:
  fourthirds = 4.0/3.0;
  onethird = fourthirds - 1.0;
  macheps = fabs(onethird + onethird + onethird - 1.0);
  // Don't make macheps too big or too small
  switch(sizeof(double)) {
    case 4: i = 0; break;
    case 8: i = 1; break;
    case 10: i = 2; break;
    case 12: i = 3; break;
    case 16: i = 4; break;
    default: i = 1; 
  }
  if (macheps < mineps[i]) macheps = mineps[i];
  if (macheps > maxeps[i]) macheps = maxeps[i];
  *ok = 1;
  if (n == 0 || d == 0) return;
  if (L==0) {
    if (Sig==0) { *ok = 0; return; }
    allocate(LSig, d*d);
  }
  else 
    LSig = L;
  if (Sig) {
    copy(d*d, Sig, 1, LSig, 1);
    potrf("Low", d, LSig, d, &info); // Cholesky factorize Sig
    if (info != 0) { // Sig may be positive semidefinite; try adding del*I to it
      imx = iamax(d, Sig, d+1); // (=blas idamax) points to max diagonal element
      del0 = macheps/10.0 + Sig[imx]*(3.0 + d/50.0)*macheps;
      *del = del0;
      while (info !=0) {
        copy(d*d, Sig, 1, LSig, 1);
        for (i=0; i<d; i++) LSig[i*(d+1)] += *del; // add del*I to Sig
        *del *= 1.5;
        potrf("Low", d, LSig, d, &info);
        if (*del > 2000*del0) { // Max del tried is
          if (L==0) freem(LSig);   // max(Sig_ii)*(6 + d/25)*1000*(machine epsilon)
          *ok = 0;
          return;
        }
      }
    }
  }  
  RandN(X, n*d, rng);
  trmm("Right", "Lower", "T", "NotUdia", n, d, 1.0, LSig, d, X, n);
  if (mu) for (i=0; i<n; i++) axpy(d, 1.0, mu, 1, X+i, n); // add mu to each row
  if (L==0) freem(LSig);
}
