#include <time.h>
#include <stdlib.h>
#include <stdint.h>
#include "random.h"
#include "debugprint.h"

static uint64_t generate_seed(uint64_t); // forward

rand_rng *rand_create(void) {
  rand_rng* rng = malloc(sizeof(rand_rng));
  if (!rng) return NULL;
  rng->type = DEFAULT_RNG;
  rand_randomize(rng);
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
// RANDOMIZE FUNCTIONS
//-------------------------------------------------------------------
void rand_randomize(rand_rng *rng) {
  if (rng->type == PARKMILLER)
    rng->PMseed = (uint32_t)generate_seed(0);
  else {
    rng->state[0] = generate_seed(0);
    rng->state[1] = generate_seed(0);
  }   
}

void rand_thread_randomize(uint64_t thread_id, rand_rng *rng) {
  if (rng->type == PARKMILLER)
    rng->PMseed = (uint32_t)generate_seed(thread_id);
  else {
    rng->state[0] = generate_seed(thread_id);
    rng->state[1] = generate_seed(thread_id);
  }   
}
// SET AND GET STATE AND TYPE FUNCTIONS
//---------------------------------------------------------------------
void rand_setstate(uint64_t state0, uint64_t state1, rand_rng *rng) {
  if (state0==0 && state1==0)
    state0 = state1 = 0xAAAAAAAA << 16;  // don't seed with 0.
  rng->state[0] = state0;
  rng->state[1] = state1;
  rand_d(rng); // spin_up
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
}

uint32_t rand_getPMseed(rand_rng *rng) {
  return rng->PMseed;
}

void rand_settype(rng_type type, rand_rng *rng) {
  rng->type = type;
  rand_randomize(rng);
}

static uint64_t generate_seed(uint64_t thread_id) {
  // Generate seed from /dev/urandom if available, othewise from
  // an xor combination of clock, pid, and thread-id
  uint64_t seed = 0;
#ifdef __unix__
  if (int fd = open("/dev/urandom", O_RDONLY) >= 0) {
    if (read(fd, &seed, sizeof(seed)) == sizeof(seed)) {
      close(fd);
      return seed;
    }
    close(fd);
  }
  seed = (uint64_t)getpid();
#elif defined(_WIN32)
  seed = (uint64_t)GetCurrentProcessId();
#endif
  struct timespec ts;
  timespec_get(&ts, TIME_UTC);
  seed ^= (uint64_t)ts.tv_sec * 1000000000 + ts.tv_nsec;
  seed ^= thread_id;
  return seed;
}
