// RandomNumbers — random number generation utilities for VARMASIM
//
// See RandomNumbers.h for further information, including parameter descriptions

#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include "randompack.h"
#include "BlasGateway.h"
#include "error.h"
//#include "ExtraUtil.h"

// These typedef-s from random.h can't be forward declared for C-technical reasons
typedef enum {   // Random number generator type, can't be forward declared
  PARKMILLER,     // Park–Miller LCG
  DEFAULT         // Xorshift128+ or R’s builtin
} rng_type;

typedef struct rand_rng {      // State holder for the active RNG
  rng_type type;   // RNG type
  uint64_t       state[2];
  uint32_t       PMseed;
  double         spare_norm; // cached spare normal or NAN if none
} rand_rng;

// Forward declarations to rand functions formerly in random.h, but now below.
static rand_rng *rand_create(void);
static void rand_free(rand_rng *rng);
static void rand_normal(double x[], int n, rand_rng *rng);
static void rand_perm(int p[], int n, rand_rng *rng);
static void rand_sample(int p[], int n, int k, rand_rng *rng);
static void rand_uint64(uint64_t bound, uint64_t x[], int n, rand_rng *rng);
static void rand_uint32(uint32_t bound, uint32_t x[], int n, rand_rng *rng);
static void rand_setstate(uint64_t state0, uint64_t state1, rand_rng *rng);
static void rand_getstate(uint64_t state[2], rand_rng *rng);
static void rand_setPMseed(uint32_t seed, rand_rng *rng);
static uint32_t rand_getPMseed(rand_rng *rng);
static void rand_settype(rng_type type, rand_rng *rng);
static uint64_t rand_splitmix64(uint64_t *x);
static inline uint64_t rand_os_pid(void);
static bool rand_get_system_entropy(uint64_t *s0, uint64_t *s1);
static void rand_randomize(uint64_t thread_id, rand_rng *rng);
static void rand_keepalive(void);

// Include file with inlined rand_ functions:
#include "random_inlined.h"
    
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

randompack_rng *randompack_create(const char *type, long long seed) {
  randompack_rng *rng = rand_create();
  if (!rng) return 0;

  // Parse type string
  rng_type rtype = DEFAULT;
  if (type) {
    if (strcmp(type, "Park-Miller")==0 || strcmp(type, "PM")==0) {
      rtype = PARKMILLER;
    }
    else if (strcmp(type, "Xorshift128+")==0 || strcmp(type, "Xorshift")==0 ||
             strcmp(type, "X+")==0) {
      rtype = DEFAULT;
    }
    else if (strcmp(type, "R")==0 || strcmp(type, "R-default")==0) {
      rtype = DEFAULT;
    }
  }
  rand_settype(rtype, rng);

  // Seed or randomize
  if (seed == 0) {
    rand_randomize(0, rng);
  }
  else if (seed < 0) {
    rand_randomize(-(uint64_t)(seed), rng); // note safe cast-before minus
  }
  else {
    uint32_t pm = (uint32_t)seed;
    uint64_t x = (uint64_t)seed;
    uint64_t s0 = rand_splitmix64(&x);
    uint64_t s1 = rand_splitmix64(&x);
    if (s0 == 0 && s1 == 0) s1 = 1;
    rand_setPMseed(pm, rng);
    rand_setstate(s0, s1, rng);
  }

  (void)rand_keepalive; // To suppress warnings of unused rand_ functions
  return rng;
}

void randompack_free(randompack_rng *rng) {
  rand_free(rng);
}
unsigned int randompack_getPMseed(randompack_rng *rng) {
  if (rng == 0) return 0;
  return (unsigned int)rng->PMseed;
}

bool randompack_u01(double x[], int len, randompack_rng *rng) {
  if (!len && !x) return true;
  ASSERT(x && rng && len >= 0, "illegal arguments to randompack_u01");
#ifdef USING_R
  if (rng->type == PARKMILLER)
    rand_dble(x, len, rng);
  else {
    GetRNGstate();
    for (int i = 0; i < len; i++) {
      x[i] = unif_rand(); // R-s built-in generator
    }
    PutRNGstate();
  }
#else
  rand_dble(x, len, rng);
#endif
  return true;
}

bool randompack_int(int x[], int len, int m, int n, randompack_rng *rng) {
  //  if (!rng || len < 0 || m > n) return;
  if (!len && !x) return true;
  int64_t span = (int64_t)n - (int64_t)m + 1;
  ASSERT(x && rng, "illegal arguments to randompack_int");
  ASSERT(len >= 0, "randompack_int: len must be ≥ 0");
  ASSERT(m <= n, "randompack_int: m must not be bigger than n");
  ASSERT(span <= INT32_MAX - 1, "maximum span for random integers is 2^31 - 1");
  rand_int((int)span, x, len, rng);
  for (int i = 0; i < len; i++)
    x[i] += m;
  return true;
}

bool randompack_perm(int x[], int len, randompack_rng *rng) {
  if (!len && !x) return true;
  ASSERT(x && rng && len >= 0, "illegal argument(s) to randompack_perm");
  rand_perm(x, len, rng);
  return true;
}

bool randompack_sample(int x[], int len, int k, randompack_rng *rng) {
  if ((!len && !x) || k == 0) return true;
  ASSERT(rng, "randompack_sample: rng must not be NULL");
  ASSERT(len >= 0, "randompack_sample: len is negative");
  ASSERT(0 <= k && k <= len, "randompack_sample: k must be in [0, len]");
  ASSERT(x, "randompack_sample: x must not be NULL when k, len > 0");
  rand_sample(x, len, k, rng);
  return true;
}

bool randompack_norm(double x[], int len, randompack_rng *rng) {
  if (!len && !x) return true;
  ASSERT(x && rng && len>=0, "illegal arguments to randompack_norm");
#ifdef USING_R
  if (rng->type == PARKMILLER)
    rand_normal(x, len, rng);
  else {
    GetRNGstate();
    for (int i = 0; i < len; i++) {
      x[i] = unif_rand(); // R-s built-in generator
    }
    PutRNGstate();
  }
#else
  rand_normal(x, len, rng);
#endif  
  return true;
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

bool randompack_mvn(char *transp, double mu[], double Sig[], int d, int n, double X[],
                    int ldX, double L[], randompack_rng *rng) {
  // If Sig, fill L if it is supplied, otherwise use L Sig is d × d, X is either n × d (if
  // transp="N...") or d × n (if transp="T...") if (!n && !X) return true. X having two
  // possible shape is modelled after dsyrk.
  ASSERT(Sig || L, "randompack_mvn: either Sig or L must be specified");
  ASSERT(X && rng && d > 0 && n >= 0, "illegal arguments to randompack_mvn");
  ASSERT(transp && strchr("NT", transp[0]), "transp must begin with \"N\" or \"T\"");
  int i, info, rank;
  bool Lalloc = false;
  bool TRAN = (transp[0] == 'T');
  //double macheps = get_macheps();
  if (Sig) {
    double *work, *A;
    int *piv;
    if (!L && X == 0) return true;
    if (!L) {
      Lalloc = true;
      ALLOC(L, d*d);
    }
    else
      laset("Upper", d, d, 0.0, 0.0, L, d);
    ALLOC(work, 2*d);
    ALLOC(piv, d);
    // printM("before lacpy, L", L, d, d);
    lacpy("Lower", d, d, Sig, d, L, d);
    potrf("Low", d, L, d, &info);
    printM("in randompack_mvn, L", L, d, d);
    printI("info", info);
    if (info > 0) {
      lacpy("Lower", d, d, Sig, d, L, d);      
      pstrf("Low", d, L, d, piv, &rank, 1e-14, work, &info); // PSD Chol fact.
      ALLOC(A, rank*d);
      copy(rank*d, L, 1, A, 1);
      // Set L to [A(piv, 1:rank) 0]
      for (int i=0; i<d; i++) copy(rank, A + piv[i], d, L + i, d);
      laset("All", d, d-rank, 0.0, 0.0, L + d*rank, d);
      FREE(A);
    }
    FREE(piv); FREE(work);
  }
  rank = LastNonzeroColumn(d, d, L);
  if (n == 0) return true;
  if (rank == d) {
    printM("L", L, d, d);
    if (TRAN) {
      for (int i=0; i<n; i++) {
	randompack_norm(X + i*ldX, d, rng);
      }
      printMP("X before transform", X, d, n, X, ldX);
      trmm("Left", "Lower", "NoT", "NotUdia", d, n, 1.0, L, d, X, ldX);
      printMP("(T) X", X, d, n, X, ldX);
    }
    else {
      for (int i=0; i<d; i++) randompack_norm(X + i*ldX, n, rng);
      printMP("X before transform", X, n, d, X, ldX);
      trmm("Right", "Lower", "T", "NotUdia", n, d, 1.0, L, d, X, ldX);
      printMP("(NT) X", X, n, d, X, ldX);
    }
  }
  else { // Singular Sig, use workspace for the N(0,1) randoms
    double *Wrk;
    ALLOC(Wrk, n*rank);
    randompack_norm(Wrk, n*rank, rng);
    if (TRAN)
      gemm("NoT", "NoT", d, n, rank, 1.0, L, d, Wrk, rank, 0.0, X, ldX);
    else
      gemm("NoT", "T", n, d, rank, 1.0, Wrk, n, L, d, 0.0, X, ldX);
    FREE(Wrk);
  }
  if (mu) {
    if (TRAN) for (i=0; i<n; i++) axpy(d, 1.0, mu, 1, X+i*d, 1); // add mu to each col
    else      for (i=0; i<n; i++) axpy(d, 1.0, mu, 1, X+i, n);   // add mu to each row
    if (Lalloc) FREE(L);
  }
  return true;
}

// ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
// The contents of random.h follow, with all functions declared static (to make the
// symbols private)

//#define static static

// These can't be forward declared for C-technical reasons
// typedef enum { // Random number genator type
//   PARKMILLER,    // Park-Miller, good for comparison across platforms
//   DEFAULT_RNG    // Fast Xorshift128+ in standalone builds, or R's builtin
// } rng_type;

// typedef struct rand_rng{ // Opaque-ish state holder for the active RNG
//   rng_type type;     // RNG type
//   uint64_t state[2]; // State for Xorshift128+
//   uint32_t PMseed;   // Seed for Park-Miller
// } rand_rng;

static rand_rng *rand_create(void) {
  rand_rng* rng = malloc(sizeof(rand_rng));
  if (!rng) return NULL;
  rng->type = DEFAULT;
  rand_randomize(0, rng);
  rng->spare_norm = NAN;
  return rng;
}

static void rand_free(rand_rng *rng) {
  FREE(rng);
}

static void rand_normal(double x[], int n, rand_rng *rng) {
  // Normal random numbers from N(0,1); uses polar method.
  double u=0.0, v=0.0, uv[2], s=0.0, R;
  if (n<=0) return;
  int k = 0;
  if (rng->spare_norm == rng->spare_norm) {
    x[k++] = rng->spare_norm;
    rng->spare_norm = NAN;
    if (k >= n) return;
  }
  while (k + 1 < n) {
    do {
      rand_dble(uv, 2, rng);
      u = uv[0]*2.0 - 1.0;
      v = uv[1]*2.0 - 1.0;
      s = u*u + v*v;
    } while (s <= 0.0 || s >= 1.0);
    R = sqrt(-2.0*log(s)/s);
    x[k++] = u*R;
    x[k++] = v*R;
  }
  if (k < n) {
    do {
      rand_dble(uv, 2, rng);
      u = uv[0]*2.0 - 1.0;
      v = uv[1]*2.0 - 1.0;
      s = u*u + v*v;
    } while (s <= 0.0 || s >= 1.0);
    R = sqrt(-2.0*log(s)/s);
    x[k] = u*R;
    rng->spare_norm = v*R;
  }
}

static void rand_perm(int p[], int n, rand_rng *rng) {
  // Fisher-Yates algorithm, stores a random permutation of 0..n-1 in p
  int i, j;
  for (i=0; i<n; i++) {
    rand_int(i+1, &j, 1, rng);
    if (j != i) p[i] = p[j];
    p[j] = i;
  }
}

static void rand_sample(int p[], int n, int k, rand_rng *rng) {
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

static void rand_uint64(uint64_t bound, uint64_t x[], int n, rand_rng *rng) {
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

static void rand_uint32(uint32_t bound, uint32_t x[], int n, rand_rng *rng) {
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

static void rand_setstate(uint64_t state0, uint64_t state1, rand_rng *rng) {
  if (state0==0 && state1==0)
    state0 = state1 = 0xAAAAAAAA << 16;  // don't seed with 0.
  rng->state[0] = state0;
  rng->state[1] = state1;
  //rand_d(rng); // spin_up
}

static void rand_getstate(uint64_t state[2], rand_rng *rng) {
  if (rng->type == PARKMILLER) {
    state[0] = rng->PMseed;
    state[1] = 0;
  }
  else {
    state[0] = state[0];
    state[1] = state[1];
  }
}

static void rand_setPMseed(uint32_t seed, rand_rng *rng) {
  rng->PMseed = seed;
  PM_rand_bits(rng); // spinup one number
}

static uint32_t rand_getPMseed(rand_rng *rng) {
  return rng->PMseed;
}

static void rand_settype(rng_type type, rand_rng *rng) {
  rng->type = type;
}

// RANDOMIZE FUNCTIONS

static uint64_t rand_splitmix64(uint64_t *x) {
  uint64_t z = (*x += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}

static inline uint64_t rand_os_pid(void) {
#ifdef __unix__
  return (uint64_t)getpid();
#elif defined(_WIN32)
  return (uint64_t)GetCurrentProcessId();
#else
  return 0;
#endif
}

static bool rand_get_system_entropy(uint64_t *s0, uint64_t *s1) {
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

static void rand_randomize(uint64_t thread_id, rand_rng *rng) {
  uint64_t s0 = 0, s1 = 0;

  if (!rand_get_system_entropy(&s0, &s1)) {
    // fallback mixing path
    uint64_t x = rand_os_pid();
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    uint64_t t = (uint64_t)ts.tv_sec * 1000000000ULL + (uint64_t)ts.tv_nsec;

    x ^= t;
    x ^= thread_id * 0x9E3779B97F4A7C15ULL;

    s0 = rand_splitmix64(&x);
    s1 = rand_splitmix64(&x);
    if (s0 == 0 && s1 == 0) {
      s0 = 1; // Assures both state and PMseed not 0
    }
  }

  // Always set both state words and the PM seed
  rng->state[0] = s0;
  rng->state[1] = s1;
  rng->PMseed = (uint32_t)(s0 & 0xFFFFFFFFUL);
}

static void rand_keepalive(void) {
  (void) rand_perm;
  (void) rand_sample;
  (void) rand_uint64;
  (void) rand_getstate;
  (void) rand_getPMseed;
  // add any other rand_* helpers you want to keep
}
