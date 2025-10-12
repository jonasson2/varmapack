#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "allocate.h"
#include "RandomNumbers.h"
#include "BlasGateway.h"
#include "xAssert.h"
//#include "allocate.h"
#include "random.h"
#ifdef USING_R
#include <R.h>
#include <Rmath.h>
#endif

// TODO  Matlab mex

void Rand ( // Generate uniform random numbers.
  double x[],    // out  n-vector, returns n uniform [0,1) random numbers
  int n,         // in   dimension of x
  rand_rng *rng) // in   random number generator
{
#ifdef USING_R
  if (rng->type == PARKMILLER)
    rand_dble(x, n, rng);
  else {
    GetRNGstate();
    for (int i = 0; i < n; i++) {
      Rprintf("i=",i);
      x[i] = unif_rand(); // R-s built-in generator
    }
    PutRNGstate();
  }
#else
  rand_dble(x, n, rng);
#endif  
}

void RandN ( // Generate normal random numbers.
  double x[],    // out  n-vector, returns n uniform [0,1) random numbers
  int n,         // in   dimension of x
  rand_rng *rng) // in   random number generator
{
#ifdef USING_R
  if (rng->type == PARKMILLER)
    rand_norm(x, n, rng);
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

void RandNM ( // Generate multivariate normal random vectors
  double mu[],   // in     d-vec., mean of generated vectors, null for zero-mean
  double Sig[],  // in     d × d, covariance of generated vectors, null to use L
  int n,         // in     Number of replicates
  int d,         // in     dimension of generated vectors
  double X[],    // out    n × d, generated vectors
  double L[],    // in/out d × d, lower Cholesky factor of Sig (or null to omit)
  double *del,   // out    multiple of I added to Sig to make it pos.def
  rand_rng *rng, // in     random number generator
  int *ok)       // out    0 if Sig could not be made postive definite, else 1
  
// NOTE 1: Sig and X are stored columnwise in Fortran fashion.
// NOTE 2: There are situations when Sig is indefinite but close to being
//         positive definite, for example due to rounding errors. An attempt
//         is made to remedy this by adding a small multiple, del, of the
//         identity matrix, I, to Sig. It was found empirically that the
//         distribution of X does not depend critically on del, and the
//         formula for del0 below was determined to be approximately the
//         minimum needed. The while loop further esures that.
// NOTE 3: In case the calling program needs the Cholesky factor of the
//         (possibly modified) Sig, this may be returned in L by letting it be
//         an d × d array instead of 0 in the call.
// NOTE 4: When RandNM is to be called multiple times for the same covariance
//         it is possible to save execution time by reusing the Cholesky
//         factorization of Sig on all calls but the first. Specify Sig and
//         return L on the first call, and let Sig be null and specify L on
//         subsequent calls. If both Sig and L are null, the function exits
//         with ok = 0.
  
{
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
