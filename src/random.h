// RANDOM NUMBERS
//
// NOTES
// * The default random number generator (rng) is Xorshift128+ from [2]. It
//   passes many quality tests, is very fast, and requires only two 64 byte
//   words of state. It is also possible to use a different default rng, e.g.
//   when calling from Matlab or R, as demonstrated by "varmac".
//
// * Park-Miller [2] is a simple rng that allows the generation of exactly the
//   same sequence of random numbers as produced by an accompanying Matlab
//   function rando.m (useful for testing and debugging). This applies to
//   rand_int, rand_d, rand_double, rand_normal, rand_perm and rand_sample. The
//   functions rand_uint32 and rand_uint64 are unaffected. The Park-Miller RNG
//   has (much) lower quality, and random ints fail for bound greater than 2^31
//   - 1, and start to become biased for a lower bound, because they are based
//   on a simple modulus.
//
// * To set the seed of the P-M generator, use rand_setseed(s, 0) with s <
//   2^31-1
//
// * The program suite is thread-safe when seeded or randomized differently in
//   each thread.
//
// REFERENCES
// [1] S.K. Park and K.W. Miller (1988). "Random number generators: Good ones
//     are hard to find". Communications of the ACM 31 (10), 1192-1201.
// [2] Vigna, Sebastiano. "Further scramblings of Marsaglia’s xorshift
//     generators." Journal of Computational and Applied Mathematics 315 (2017):
//     175-181.

#ifndef RANDOM_H
#define RANDOM_H t

#include <inttypes.h>
#include <stdbool.h>

// TYPEDEFS
typedef enum {
  PARKMILLER,
  DEFAULT_RNG
} rng_type;

typedef struct {
    rng_type type;
    uint64_t state[2]; // State for the default Xorshift128+
    uint32_t PMseed;   // Seed for Park-Miller
} rand_rng;

// CREATE AND FREE AN RNG
rand_rng* rand_create(void);
void rand_free(rand_rng *rng);

// INLINE FUNCTIONS (FOR SPEED)
// _int gives integers in {0,...,bound}; _i gives one integer
// _dble gives doubles in [0, 1); _d gives one double
static inline void rand_int(int bound, int x[], int n, rand_rng *rng);
static inline void rand_dble(double x[], int n, rand_rng *rng);
static inline int rand_i(int, rand_rng *rng);
static inline double rand_d(rand_rng *rng);

// CREATE RANDOM NUMBERS
// _normal provides standard normal randoms
// _perm is a random permutation of 0...n-1
// _sample gives k values from {0...n-1} without replacement
// _uint64 and _uint32 are also in {0,...,bound}
void rand_normal(double x[], int n, rand_rng *rng);
void rand_perm(int p[], int n, rand_rng *rng);
void rand_sample(int p[], int n, int k, rand_rng *rng);
void rand_uint64(uint64_t bnd, uint64_t x[], int n, rand_rng *rng);
void rand_uint32(uint32_t bnd, uint32_t x[], int n, rand_rng *rng);

// RANDOMIZE FUNCTIONS (BOTH STATE AND SEED)
void rand_randomize(rand_rng* rng);
void rand_thread_randomize(uint64_t thread_id, rand_rng*rng);

// SET AND GET STATE AND TYPE
void rand_setstate(uint64_t state0, uint64_t state1, rand_rng *rng);
void rand_getstate(uint64_t state[2], rand_rng *rng);
void rand_setPMseed(uint32_t seed, rand_rng* rng);
uint32_t rand_getPMseed(rand_rng *rng);
void rand_settype(rng_type type, rand_rng *rng); // randomizes
rng_type rand_gettype(rand_rng *rng);

// IMPLEMENTATION
#include "random_inlined.h"

#endif
