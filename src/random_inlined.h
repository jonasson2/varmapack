// DEFINITIONS OF STATIC INLINE FUNCTION DECLARED IN RANDOM.H
#ifndef RANDOM_INLINED_H
#define RANDOM_INLINED_H

#include <stddef.h>
#include <inttypes.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>

static const uint32_t mersenne8 = 2147483647;

static inline uint32_t PM_rand_bits(rand_rng *rng) {
  const uint32_t a = 16807, m = mersenne8, q = m/a, r = m%a;
  uint32_t s1, s2;
  s1 = a*(rng->PMseed % q);
  s2 = r*(rng->PMseed/q);
  if (s1 > s2) rng->PMseed = s1 - s2; else rng->PMseed = (m - s2) + s1;
  return rng->PMseed;
}

#include <stdint.h>
#include <assert.h>

static inline void PM_rand_int(int bound, int x[], int n, rand_rng *rng) {
  assert(bound > 0);
//#if UINT_MAX == UINT32_MAX // int is 32 bits
  uint32_t ubound = (uint32_t)bound;
  for (int i = 0; i < n; i++) {
    uint32_t r = PM_rand_bits(rng);      // returns uint32_t
    x[i] = (int)(r % ubound);            // result is in [0, bound-1]
  }
//#else
//#error "Unsupported int size with Park-Miller"
//#endif
}

static inline void rand_int(int bound, int x[], int n, rand_rng *rng) {
  // TODO fix this (make it work for any int size, and check also for ILP64
  // by using uint64_t instead of int in the testing
  assert(bound >=0);  
#if UINT_MAX == UINT32_MAX // int is 32 bits
  if (rng->type == PARKMILLER)
    PM_rand_int(bound, x, n, rng);
  else
    rand_uint32((uint32_t) bound, (uint32_t*) x, n, rng);
#elif UINT_MAX == UINT64_MAX
  assert(rng->type != PARKMILLER); // Unsupported int size with Park-Miller
  rand_uint64((uint64_t) bound, (uint64_t*) x, n, rng);
  // Needs an ILP64 compiler to test, e.g. Cray
#else
#error "Unsupported int size"
#endif
}
// XORSHIRO128+
//
// static inline uint64_t rotl(const uint64_t x, int k) {
// 	return (x << k) | (x >> (64 - k));
// }
// 
// static uint64_t s[2];
//
// uint64_t next(void) {
// 	const uint64_t s0 = s[0];
// 	uint64_t s1 = s[1];
// 	const uint64_t result = s0 + s1;
// 	s1 ^= s0;
// 	s[0] = rotl(s0, 24) ^ s1 ^ (s1 << 16); // a, b
// 	s[1] = rotl(s1, 37); // c
// 	return result;
// }

static inline uint64_t rand_64bits(rand_rng *rng) {
  // From ref. [2] (ree random.h)
  uint64_t s1 = rng->state[0];
  uint64_t s0 = rng->state[1];
  rng->state[0] = s0;
  s1 ^= s1 << 23; // a
  uint64_t result = s0 + s1;
  rng->state[1] = s1^s0^(s1 >> 18)^(s0 >> 5); // b, c
  return result;
}

// Older version, from Wikipedia ~2019
// static inline uint64_t rand_64bits(rng) {
//   uint64_t s1 = rng->state[0];
//   uint64_t s0 = rng->state[1];
//   rng->state[0] = s0;
//   s1 ^= s1 << 23; // a
//   return (rng->state[1] = (s1^s0^(s1 >> 17 )^(s0 >> 26))) + s0; // b, c
// }

static inline uint32_t rand_32bits(rand_rng *rng) {
  return (uint32_t) rand_64bits(rng);
}

static inline void rand_dble(double x[], int n, rand_rng *rng) {
  int i;
  // Random uniform double
  if (rng->type == PARKMILLER) {
    const double maxrand = 1.0 + (double)mersenne8;
    for (i = 0; i < n; i++) {
      uint32_t u = PM_rand_bits(rng);
      double r = (double)u;
      x[i] = r / maxrand;
    }
  }
  else { // Xorshift128+
    for (i = 0; i < n; i++) {
      uint64_t v = rand_64bits(rng);
      x[i] = (double)(v >> 11) * 0x1.0p-53;
    }
  }
}

static inline double rand_d(rand_rng *rng) {
  double r;
  rand_dble(&r, 1, rng);
  return r;
}

static inline int rand_i(int bound, rand_rng *rng) {
  int i;
  rand_int(bound, &i, 1, rng);
  return i;
}

#endif

