// TimeSimulate.c: time varmapack_sim on named and representative unnamed models.

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "getopt.h"
#include "varmapack.h"

typedef struct {
  char *label;
  bool rho_case;
  int p;
  int q;
  int r;
  double rho;
} bench_case;

typedef struct {
  double t;
  double w;
  int n;
  int M;
} options;

static uint64_t clock_nsec(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return 1000000000ull*(uint64_t)ts.tv_sec + (uint64_t)ts.tv_nsec;
}

static void consume_double(double x) {
  static volatile double sink;
  sink += x;
}

static void die(char *msg) {
  fprintf(stderr, "%s\n", msg);
  exit(1);
}

static void die_last_error(char *where) {
  fprintf(stderr, "%s failed: %s\n", where, varmapack_last_error());
  exit(1);
}

static void default_options(options *opts) {
  *opts = (options) {.t = .1, .w = .1, .n = 100, .M = 1000};
}

static void print_help(void) {
  printf("TimeSimulate -- time direct C-library simulation on built-in testcases\n");
  printf("Usage: TimeSimulate [options]\n\n");
  printf("This permits comparison with the corresponding MATLAB and Python interface\n");
  printf("benchmarks.\n");
  printf("Times all named testcases and an unnamed p=q=3, r=10, rho=.98 model.\n\n");
  printf("Options:\n");
  printf("  -h          show this help\n");
  printf("  -t seconds  timing target per testcase (default 0.1)\n");
  printf("  -w seconds  CPU warmup time before timing (default 0.1)\n");
  printf("  -n length   length of each series (default 100)\n");
  printf("  -M count    replicates per call (default 1000)\n");
}

static bool get_options(int argc, char **argv, options *opts, bool *help) {
  int opt;
  default_options(opts);
  *help = false;
  opterr = 0;
  optind = 1;
  while ((opt = getopt(argc, argv, "ht:w:n:M:")) != -1) {
    switch (opt) {
      case 'h':
        *help = true;
        return true;
      case 't':
        opts->t = atof(optarg);
        if (opts->t <= 0) return false;
        break;
      case 'w':
        opts->w = atof(optarg);
        if (opts->w < 0) return false;
        break;
      case 'n':
        opts->n = atoi(optarg);
        if (opts->n <= 0) return false;
        break;
      case 'M':
        opts->M = atoi(optarg);
        if (opts->M <= 0) return false;
        break;
      default:
        return false;
    }
  }
  return optind == argc;
}

static void warm_cpu(double seconds) {
  double x = 1.000001;
  uint64_t start, deadline;
  if (seconds <= 0) return;
  start = clock_nsec();
  deadline = start + (uint64_t)(seconds*1e9);
  while (clock_nsec() < deadline) {
    for (int i=0; i<1000; i++) {
      x += log(x);
      if (x > 2) x = 1.000001;
    }
  }
  consume_double(x);
}

static void allocate_case(int p, int q, int r, double **A, double **B, double **Sig) {
  *A = malloc(sizeof(**A)*r*r*(p > 0 ? p : 1));
  *B = malloc(sizeof(**B)*r*r*(q > 0 ? q : 1));
  *Sig = malloc(sizeof(**Sig)*r*r);
  if (!*A || !*B || !*Sig) die("allocation failed");
}

static void construct_case(bench_case c, randompack_rng *rng, double **A, double **B,
                           double **Sig, int *p, int *q, int *r) {
  char name[VARMAPACK_TESTCASE_NAME_LEN];
  int index = c.rho_case ? 0 : 1;
  *p = c.p;
  *q = c.q;
  *r = c.r;
  snprintf(name, sizeof(name), "%s", c.rho_case ? "rho" : c.label);
  if (!c.rho_case && !varmapack_testcase(name, &index, 0, p, q, r, 0, 0, 0, rng)) {
    die_last_error("varmapack_testcase");
  }
  allocate_case(*p, *q, *r, A, B, Sig);
  if (!varmapack_testcase(name, &index, c.rho, p, q, r, *A, *B, *Sig, rng)) {
    die_last_error("varmapack_testcase");
  }
}

static double time_simulation(double A[], double B[], double Sig[], int p, int q, int r,
                              int n, int M, double target, randompack_rng *rng) {
  double *X;
  uint64_t start, t;
  int reps = 0;
  X = malloc(sizeof(*X)*r*n*M);
  if (X == 0) die("allocation failed");
  if (!randompack_seed(12345, 0, 0, rng)) die("randompack_seed failed");
  start = clock_nsec();
  t = start;
  while ((t - start)*1e-9 < target) {
    if (!varmapack_sim(A, B, Sig, 0, 0, p, q, r, n, M, 0, 0, 1, X, 0, rng)) {
      die_last_error("varmapack_sim");
    }
    consume_double(X[r*n*M - 1]);
    reps++;
    t = clock_nsec();
  }
  free(X);
  return (double)(t - start)/((double)reps*r*n*M);
}

static void time_case(bench_case c, options *opts, randompack_rng *rng) {
  double *A, *B, *Sig;
  double rho, total;
  int p, q, r;
  construct_case(c, rng, &A, &B, &Sig, &p, &q, &r);
  rho = p == 0 ? 0 : varmapack_specrad(A, r, p);
  if (isnan(rho)) die_last_error("varmapack_specrad");
  total = time_simulation(A, B, Sig, p, q, r, opts->n, opts->M, opts->t, rng);
  printf("%-12s %2d %2d %2d %5.2f %10.1f\n", c.label, p, q, r, rho, total);
  free(Sig);
  free(B);
  free(A);
}

int main(int argc, char **argv) {
  bench_case cases[] = {
    {"tinyAR", false, 0, 0, 0, 0}, {"tinyMA", false, 0, 0, 0, 0},
    {"tinyARMA", false, 0, 0, 0, 0}, {"smallAR1", false, 0, 0, 0, 0},
    {"smallAR2", false, 0, 0, 0, 0}, {"smallMA1", false, 0, 0, 0, 0},
    {"smallMA2", false, 0, 0, 0, 0}, {"smallARMA1", false, 0, 0, 0, 0},
    {"smallARMA2", false, 0, 0, 0, 0}, {"mediumAR", false, 0, 0, 0, 0},
    {"mediumMA1", false, 0, 0, 0, 0}, {"mediumARMA1", false, 0, 0, 0, 0},
    {"mediumARMA2", false, 0, 0, 0, 0}, {"mediumMA2", false, 0, 0, 0, 0},
    {"largeAR", false, 0, 0, 0, 0},
    {"Unnamed", true, 3, 3, 10, .98}
  };
  options opts;
  randompack_rng *rng;
  bool help;
  if (!get_options(argc, argv, &opts, &help) || help) {
    print_help();
    return help ? 0 : 1;
  }
  rng = randompack_create(0);
  if (rng == 0) die("randompack_create failed");
  printf("Varmapack C benchmark\n");
  printf("Benchmark unit:          ns/value\n");
  printf("Length per series:       %d\n", opts.n);
  printf("Replicates per call:     %d\n", opts.M);
  printf("Benchmark time per case: %.1f s\n\n", opts.t);
  warm_cpu(opts.w);
  printf("%-12s %2s %2s %2s %5s %10s\n", "Testcase", "p", "q", "r", "rho",
         "Varmapack");
  for (int i=0; i<(int)(sizeof(cases)/sizeof(*cases)); i++) {
    time_case(cases[i], &opts, rng);
  }
  randompack_free(rng);
  return 0;
}
