#include <time.h>
#include <stdlib.h>
#include <stdint.h>
#include "random.h"
#include "debugprint.h"

rand_rng *rand_create(void) {
  rand_rng* rng = malloc(sizeof(rand_rng));
  if (!rng) return NULL;
  rng->type = DEFAULT_RNG;
  rand_randomize(0, rng);
  return rng;
}

void rand_free(rand_rng *rng) {
  free(rng);
}

void rand_normal(double x[], int n, rand_rng *rng) {
  // Normal random numbers from N(0,1); uses polar method.
  double u=0.0, v=0.0, uv[2], s=0.0, R;
  uint64_t k;
  if (n<=0) return;
  if (n<=2) {
    while (1) {
      rand_dble(uv, 2, rng);
      u = uv[0]*2.0 - 1.0;
      v = uv[1]*2.0 - 1.0;
      s = u*u + v*v;
      if (s > 0.0 && s < 1.0) break;
    }
    R = sqrt(-2.0*log(s)/s);
    x[0] = u*R;
    if (n==2) x[1] = v*R;
  }
  else {
    k = 0;
    while (k < (uint64_t) n) {
      rand_dble(uv, 2, rng);
      u = uv[0]*2.0 - 1.0;
      v = uv[1]*2.0 - 1.0;
      s = u*u + v*v;
      if (s > 0.0 && s < 1.0) {
        R = sqrt(-2.0*log(s)/s);
        x[k++] = u*R;
        if (k < (uint64_t) n) x[k++] = v*R;
      }
    }
  }
}

void rand_perm(int p[], int n, rand_rng *rng) {
  // Fisher-Yates algorithm, stores a random permutation of 0..n-1 in p
  int i, j;
  for (i=0; i<n; i++) {
    rand_int(i+1, &j, 1, rng);
    if (j != i) p[i] = p[j];
    p[j] = i;
  }
}

void rand_sample(int p[], int n, int k, rand_rng *rng) {
  // Select k values from 0..n-1, c.f. reservoir sampling on Wikipedia
  int i, r;
  p[0] = 0;
  for (i=0; i<k; i++) {
    rand_int(i+1, &r, 1, rng);
    if (r != i) p[i] = p[r];
    p[r] = i;
  }
  for (; i<n; i++) {
    rand_int(i+1, &r, 1, rng);
    if (r < k) p[r] = i;
  }
}

void rand_uint64(uint64_t bound, uint64_t x[], int n, rand_rng *rng) {
  // Sets x to n random integers in {0, 1,..., bound-1}
  int i = 0, j;
  uint64_t R;
  if (bound <= UINT8_MAX) {
    uint8_t r, b = (uint8_t)bound, thresh = (1<<8) % b;
    while (i < n) {
      R = rand_64bits(rng);
      for (j = 0; j<8; j++) {
        r = (uint8_t)R;
        if (r >= thresh) x[i++] = (uint64_t)(r % b);
        if (i>=n) break;
        R >>= 8;
      }
    }
  }
  else if (bound <= (uint64_t)UINT16_MAX) {
    uint16_t r, b = (uint16_t)bound, thresh = (1<<16) % b;
    while (i < n) {
      R = rand_64bits(rng);
      for (j = 0; j<4; j++) {
        r = (uint16_t)R;
        if (r >= thresh) x[i++] = (uint64_t)(r % b);
        if (i>=n) break;
        R >>= 16;
      }
    }
  }
  else if (bound <= (uint64_t)UINT32_MAX) {
    uint32_t r, b = (uint32_t)bound, M = 0-b, thresh = M % b;
    while (i < n) {
      R = rand_64bits(rng);
      r = (uint32_t)R;
      if (r >= thresh) x[i++] = (uint64_t)(r % b);
      if (i>=n) break;
      R >>= 32;
      r = (uint32_t)R;
      if (r >= thresh) x[i++] = (uint64_t)(r % b);
    }
  }
  else {
    uint64_t M = 0-bound, thresh = M % bound;
    while (i < n) {
      R = rand_64bits(rng);
      if (R >= thresh) x[i++] = R % bound;
    }
  }
}

void rand_uint32(uint32_t bound, uint32_t x[], int n, rand_rng *rng) {
  // Sets x to n random integers in {0, 1,..., bound-1}
  int i = 0, j;
  uint64_t R;
  if (bound <= (uint32_t)UINT8_MAX) {
    uint8_t r, b = (uint8_t)bound, thresh = (1<<8) % b;
    while (i < n) {
      R = rand_64bits(rng);
      for (j = 0; j<8; j++) {
        r = (uint8_t)R;
        if (r >= thresh) x[i++] = (uint32_t)(r % b);
        if (i>=n) break;
        R >>= 8;
      }
    }
  }
  else if (bound <= (uint32_t)UINT16_MAX) {
    uint16_t r, b = (uint16_t)bound, thresh = (1<<16) % b;
    while (i < n) {
      R = rand_64bits(rng);
      for (j = 0; j<4; j++) {
        r = (uint16_t)R;
        if (r >= thresh) x[i++] = (uint32_t)(r % b);
        if (i>=n) break;
        R >>= 16;
      }
    }
  }
  else {
    uint32_t r, b = (uint32_t)bound, M = 0-b, thresh = M % b;
    while (i < n) {
      R = rand_64bits(rng);
      r = (uint32_t)R;
      if (r >= thresh) x[i++] = (uint32_t)(r % b);
      if (i>=n) break;
      R >>= 32;
      r = (uint32_t)R;
      if (r >= thresh) x[i++] = (uint32_t)(r % b);
    }
  }
}

// SET AND GET STATE AND TYPE FUNCTIONS

void rand_setstate(uint64_t state0, uint64_t state1, rand_rng *rng) {
  if (state0==0 && state1==0)
    state0 = state1 = 0xAAAAAAAA << 16;  // don't seed with 0.
  rng->state[0] = state0;
  rng->state[1] = state1;
  //rand_d(rng); // spin_up
}

void rand_getstate(uint64_t state[2], rand_rng *rng) {
  if (rng->type == PARKMILLER) {
    state[0] = rng->PMseed;
    state[1] = 0;
  }
  else {
    state[0] = state[0];
    state[1] = state[1];
  }
}

void rand_setPMseed(uint32_t seed, rand_rng *rng) {
  rng->PMseed = seed;
  PM_rand_bits(rng); // spinup one number
}

uint32_t rand_getPMseed(rand_rng *rng) {
  return rng->PMseed;
}

void rand_settype(rng_type type, rand_rng *rng) {
  rng->type = type;
}

// RANDOMIZE FUNCTIONS

static uint64_t splitmix64(uint64_t *x) {
  uint64_t z = (*x += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}

static inline uint64_t os_pid(void) {
#ifdef __unix__
  return (uint64_t)getpid();
#elif defined(_WIN32)
  return (uint64_t)GetCurrentProcessId();
#else
  return 0;
#endif
}

static bool get_system_entropy(uint64_t *s0, uint64_t *s1) {
// Try to fill *s0, *s1 from OS entropy. Return true on success, false otherwise.
#ifdef __unix__
  int fd = open("/dev/urandom", O_RDONLY);
  if (fd < 0) {
    return false;
  }
  ssize_t r0 = read(fd, s0, sizeof(*s0));
  ssize_t r1 = read(fd, s1, sizeof(*s1));
  close(fd);
  if (r0 == sizeof(*s0) && r1 == sizeof(*s1)) {
    if (*s0 == 0 && *s1 == 0) {
      *s1 = 1;
    }
    return true;
  }
  return false;
#elif defined(_WIN32)
  if (BCRYPT_SUCCESS(BCryptGenRandom(NULL, (PUCHAR)s0, sizeof(*s0),
				     BCRYPT_USE_SYSTEM_PREFERRED_RNG))
      && BCRYPT_SUCCESS(BCryptGenRandom(NULL, (PUCHAR)s1, sizeof(*s1),
					BCRYPT_USE_SYSTEM_PREFERRED_RNG))) {
    if (*s0 == 0 && *s1 == 0) *s1 = 1;
    return true;
  }
  else return false;
#else
  (void) s0;
  (void) s1;
  return false;
#endif
}

void rand_randomize(uint64_t thread_id, rand_rng *rng) {
  uint64_t s0 = 0, s1 = 0;

  if (!get_system_entropy(&s0, &s1)) {
    // fallback mixing path
    uint64_t x = os_pid();
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    uint64_t t = (uint64_t)ts.tv_sec * 1000000000ULL + (uint64_t)ts.tv_nsec;

    x ^= t;
    x ^= thread_id * 0x9E3779B97F4A7C15ULL;

    s0 = splitmix64(&x);
    s1 = splitmix64(&x);
    if (s0 == 0 && s1 == 0) {
      s1 = 1;
    }
  }

  // Always set both state words and the PM seed
  rng->state[0] = s0;
  rng->state[1] = s1;
  rng->PMseed = (uint32_t)(s0 & 0xFFFFFFFFUL);
}
