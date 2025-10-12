#include <stdbool.h>
#include "random.h"
#ifndef RANDOMNUMBERS
#define RANDOMNUMBERS

// Before using the routines Rand, RandN, and RandNM, declared below, a random
// number generator should be declared and created, and afterwards it should be
// freed. Also it is possible to select a random number generator, set its state
// and possibly randomize the generator with entropy from the system.
//
// When used from R, the default rng is the R built-in, but elsewhere it is
// xorshift128+. Both from R and elsewhere one can also select Park-Miller,
// to allow direct comparison with results obtained from Matlab. The state of
// either default generator is independent of the Park-Miller seed; they carry
// on from where they left off when Park-Miller was selected. The R generator
// should be randomized/seeded from R.
//
// Typedefs and functions for these operations, declared in random.h, are:
//
//   PARKMILLER                  enum constant for rng type Park-Miller
//   DEFAULT_RNG                 enum constant for the default rng
//   rand_rng                    typedef for a random number generator
//   rand_create()               returns an initialized default rng
//   rand_free(rng)              stop using rng and free its memory
//   rand_randomize(rng)         use system entropy to set the rng state
//   rand_settype(type, rng)     set type (PARKMILLER or DEFAULT_RNG)
//   rand_seedPM(s, rng)         seed the Park-Miller generator
//   rand_setstate(s0, s1, rng)  seed the xorshift128+ generator

void Rand( // Generate uniform random numbers.
  double x[],     // out  n-vector, returns n uniform random numbers in [0,1)
  int n,          // in   dimension of x
  rand_rng *rng); // in  random number generator

void RandN( // Generate standard normal random numbers
  double x[],     // out  n-vector, returns n normal random numbers from N(0,1)
  int n,          // in   dimension of x
  rand_rng *rng); // in   random number generator

void RandNM ( // Generate multivariate normal random vectors
  double mu[],  // in     d-vec., mean of generated vectors, null for zero-mean
  double Sig[], // in     d×d, covariance of generated vectors, null to use L
  int n,        // in     Number of replicates
  int d,        // in     dimension of generated vectors
  double X[],   // out    n×d, generated vectors
  double L[],   // in/out d×d, lower Cholesky factor of Sig (or null to omit)
  double *del,  // out    multiple of I that was added to Sig to make it pos.def
  rand_rng *rng,// in     random number generator
  int *ok);     // out    0 if Sig could not be made postive definite, else 1
                // NOTE: Sig and X are stored columnwise in Fortran fashion.
#endif
