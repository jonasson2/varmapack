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
#include "ExtraUtil.h"

#ifdef USING_R
#include <R.h>
#include <Rmath.h>
#endif

// TODO  Matlab mex

static int LastNonzeroColumn(int m, int n, double A[]) {
  // If they are all zero, return 0.
  int k = n-1;
  while (k >= 0) {
    int imx = iamax(m, A+k*m, 1);
    if (A[imx + k*m] > 0) return k + 1;
    k--;
  }
  return k + 2;
}

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

void RandThreadRandomize(uint64_t thread_id, RandRng *rng) {
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

// static double get_macheps(void) {
//   int i;
//   double macheps, fourthirds, onethird;
//   double maxeps[] = {1.e-5, 2.e-14, 2.e-16, 1.e-17, 2.e-32};
//   double mineps[] = {1.e-8, 2.e-17, 2.e-19, 1.e-20, 2.e-35};
//   // Determine machine epsilon with the method of [3]:
//   fourthirds = 4.0/3.0;
//   onethird = fourthirds - 1.0;
//   macheps = fabs(onethird + onethird + onethird - 1.0);
//   // Don't make macheps too big or too small
//   switch(sizeof(double)) {
//     case 4: i = 0; break;
//     case 8: i = 1; break;
//     case 10: i = 2; break;
//     case 12: i = 3; break;
//     case 16: i = 4; break;
//     default: i = 1; 
//   }
//   if (macheps < mineps[i]) macheps = mineps[i];
//   if (macheps > maxeps[i]) macheps = maxeps[i];
//   return macheps;
// }

void RandNM (char *transp, double mu[], double Sig[], int d, int n, double X[], double
	     L[], RandRng *rng) {
  // If Sig, fill L if it is supplied, otherwise use L
  // Sig is d × d, X is either n × d (if transp="N...") or d × n (if transp="T...")
  int i, info, rank;
  bool Lalloc = false;
  bool TRAN = (transp[0] == 'T');
  //double macheps = get_macheps();
  if (n == 0 || d == 0) return;
  xAssert(L != 0 || Sig != 0);
  if (Sig) {
    double *work, *A;
    int *piv;
    if (!L) {
      Lalloc = true;
      allocate(L, d*d);
    }
    else
      laset("Upper", d, d, 0.0, 0.0, L, d);
    allocate(work, 2*d);
    allocate(piv, d);
    // printM("before lacpy, L", L, d, d);
    lacpy("Lower", d, d, Sig, d, L, d);
    potrf("Low", d, L, d, &info);
    if (info > 0) {
      lacpy("Lower", d, d, Sig, d, L, d);      
      pstrf("Low", d, L, d, piv, &rank, 1e-14, work, &info); // PSD Chol fact.
      allocate(A, rank*d);
      copy(rank*d, L, 1, A, 1);
      // Set L to [A(piv, 1:rank) 0]
      for (int i=0; i<d; i++) copy(rank, A + piv[i], d, L + i, d);
      laset("All", d, d-rank, 0.0, 0.0, L + d*rank, d);
      free(A);
    }
    freem(piv); freem(work);
  }
  rank = LastNonzeroColumn(d, d, L);
  if (rank == d) {
    RandN(X, n*d, rng);
    if (TRAN) {
      trmm("Left", "Lower", "NoT", "NotUdia", d, n, 1.0, L, d, X, d);
    }
    else
      trmm("Right", "Lower", "T", "NotUdia", n, d, 1.0, L, d, X, n);
  }
  else { // Singular Sig, use workspace for the N(0,1) randoms
    double *Wrk;
    allocate(Wrk, n*rank);
    RandN(Wrk, n*rank, rng);
    if (TRAN)
      gemm("NoT", "NoT", d, n, rank, 1.0, L, d, Wrk, n, 0.0, X, d);
    else
      gemm("NoT", "T", n, d, rank, 1.0, Wrk, n, L, d, 0.0, X, n);
    free(Wrk);
  }
  if (mu) {
    if (TRAN) for (i=0; i<n; i++) axpy(d, 1.0, mu, 1, X+i*d, 1); // add mu to each col
    else      for (i=0; i<n; i++) axpy(d, 1.0, mu, 1, X+i, n);   // add mu to each row
    if (Lalloc) freem(L);
  }
}
