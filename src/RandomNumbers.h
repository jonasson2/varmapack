// RandomNumbers — random number generation utilities for VARMASIM
//
// This header provides user-facing functions for generating random values
// drawn from uniform, normal, and multivariate normal distributions, along
// with convenience wrappers for creating, seeding, and randomizing RNGs.
//
// DEPENDENCIES:
//   It builds on the low-level RNG framework defined in random.h, which is
//   included below so that users need not include it themselves (but can
//   for extra features, e.g. complete control over rng state).
//
// INTEGRATION WITH R:
//   When compiled with -DUSING_R, the default RNG redirects calls to R's
//   built-in random number generator for reproducibility inside R packages.
//   Otherwise, the standalone implementation uses Xorshift128+ as the
//   DEFAULT_RNG and Park–Miller for cross-platform comparison.
//
// NOTES:
//   Implementation details and additional references are in RandomNumbers.c.

#ifndef RANDOMNUMBERS_H
#define RANDOMNUMBERS_H

#include <stdbool.h>
#include <stdint.h>

typedef struct rand_rng randompack_rng;

randompack_rng *randompack_create( // Create RNG with specified type and seed
  char *type,  // in   Park-Miller/PM, Xorshift128+/Xorshift/X+, R/R-default
  int seed     // in   0 to randomize, >0 to seed, <0 for thread randomize
);

void randompack_free(   // Free an RNG created with randompack_create
  randompack_rng *rng   // in   Random number generator
);

void randompack_u01(    // Generate uniform random numbers in [0,1)
  double x[],           // out  n-vector: uniform random numbers in [0,1)
  int n,                // in   Number of variates
  randompack_rng *rng   // in   Random number generator
);

void randompack_norm(   // Generate standard normal random numbers N(0,1)
  double x[],           // out  n-vector: standard normal random numbers
  int n,                // in   Number of variates
  randompack_rng *rng   // in   Random number generator
);

void randompack_mvn(    // Generate multivariate normal random vectors N(mu,Sig)
  char *transp,         // in     "N" to get n×d X, "T" to get d×n X
  double mu[],          // in     d-vector: mean (NULL → zero-mean)
  double Sig[],         // in     d×d covariance matrix (NULL → use L as-is)
  int d,                // in     Dimension of each vector
  int n,                // in     Number of replicates
  double X[],           // out    n×d or d×n matrix of generated vectors
  double L[],           // in/out d×d lower Cholesky factor of Sig (or NULL)
  randompack_rng *rng   // in     Random number generator
);
// NOTE 1: Sig, X and L are stored columnwise in Fortran fashion.
// NOTE 2: There are situations when Sig is indefinite but close to being positive
//         definite, for example due to rounding errors. To remedy this a pivoted
//         Cholesky factorization of Sig is employed, which may return an L matrix
//         which is not lower and has some trailing columns 0, but does satisfy
//         L·L' = Sig, and can thus be used to generate vectors with the right
//         distribution.
// NOTE 3: In case the calling program needs the Cholesky factor of Sig this may
//         be returned in L by letting it be a d × d array instead of 0 in the
//         call.
// NOTE 4: When randompack_mvn is to be called multiple times for the same
//         covariance it is possible to save execution time by reusing the
//         Cholesky factorization of Sig on all calls but the first. Specify Sig
//         and return L on the first call, and let Sig be null and specify L on
//         subsequent calls. If both Sig and L are null, the function exits with
//         ok = 0.

#endif /* RANDOMNUMBERS_H */
