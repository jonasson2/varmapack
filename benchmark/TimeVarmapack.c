// TimeVarmapack.c: time current C varmapack_sim implementation.

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "getopt.h"
#include "varmapack.h"

enum { MAX_VALUES = 32, MAX_CASES = 128 };

typedef struct {
  int p;
  int q;
  int r;
  int n;
  int M;
} bench_case;

typedef struct {
  double t;
  double w;
  int d;
  bench_case ref;
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

static void default_options(options *opts) {
  *opts = (options) {
    .t = 0.2, .w = 0.1, .d = 2, .ref = {3, 3, 5, 100, 100}
  };
}

static void print_help(void) {
  printf("TimeVarmapack -- time current C varmapack_sim (ns/value)\n");
  printf("Usage: TimeVarmapack [options]\n\n");
  printf("Options:\n");
  printf("  -h          show this help\n");
  printf("  -t seconds  timing target per case (default 0.2)\n");
  printf("  -w seconds  CPU warmup time before timing (default 0.1)\n");
  printf("  -d digits   printed digits (default 2)\n");
  printf("  -p p        reference p (default 3)\n");
  printf("  -q q        reference q (default 3)\n");
  printf("  -r r        reference r (default 5)\n");
  printf("  -n n        reference n (default 100)\n");
  printf("  -M M        reference M (default 100)\n");
}

static bool get_options(int argc, char **argv, options *opts, bool *help) {
  int opt;
  default_options(opts);
  *help = false;
  opterr = 0;
  optind = 1;
  while ((opt = getopt(argc, argv, "ht:w:d:p:q:r:n:M:")) != -1) {
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
      case 'd':
        opts->d = atoi(optarg);
        if (opts->d < 0) return false;
        break;
      case 'p':
        opts->ref.p = atoi(optarg);
        if (opts->ref.p < 0) return false;
        break;
      case 'q':
        opts->ref.q = atoi(optarg);
        if (opts->ref.q < 0) return false;
        break;
      case 'r':
        opts->ref.r = atoi(optarg);
        if (opts->ref.r <= 0) return false;
        break;
      case 'n':
        opts->ref.n = atoi(optarg);
        if (opts->ref.n <= 0) return false;
        break;
      case 'M':
        opts->ref.M = atoi(optarg);
        if (opts->ref.M <= 0) return false;
        break;
      default:
        return false;
    }
  }
  return optind == argc;
}

static bool same_case(bench_case a, bench_case b) {
  return a.p == b.p && a.q == b.q && a.r == b.r && a.n == b.n && a.M == b.M;
}

static void add_case(bench_case cases[], int *ncases, bench_case row) {
  for (int i = 0; i < *ncases; i++) {
    if (same_case(cases[i], row)) return;
  }
  if (*ncases >= MAX_CASES) die("too many benchmark cases");
  cases[*ncases] = row;
  (*ncases)++;
}

static bool case_less(bench_case a, bench_case b) {
  if (a.r != b.r) return a.r < b.r;
  if (a.n != b.n) return a.n < b.n;
  if (a.M != b.M) return a.M < b.M;
  if (a.p != b.p) return a.p < b.p;
  return a.q < b.q;
}

static void sort_cases(bench_case cases[], int ncases) {
  for (int i = 1; i < ncases; i++) {
    bench_case row = cases[i];
    int j = i;
    while (j > 0 && case_less(row, cases[j-1])) {
      cases[j] = cases[j-1];
      j--;
    }
    cases[j] = row;
  }
}

static int make_cases(options *opts, bench_case cases[]) {
  int p_axis[] = {1, 3, 5};
  int q_axis[] = {1, 3, 5};
  int r_axis[] = {2, 5, 12, 32};
  int n_axis[] = {5, 10, 100, 10000};
  int M_axis[] = {1, 10, 100, 10000};
  bench_case base = opts->ref;
  int ncases = 0;
  add_case(cases, &ncases, base);
  for (int i = 0; i < 3; i++) {
    bench_case row = base;
    row.p = p_axis[i];
    add_case(cases, &ncases, row);
  }
  for (int i = 0; i < 3; i++) {
    bench_case row = base;
    row.q = q_axis[i];
    add_case(cases, &ncases, row);
  }
  for (int i = 0; i < 4; i++) {
    bench_case row = base;
    row.r = r_axis[i];
    add_case(cases, &ncases, row);
  }
  for (int i = 0; i < 4; i++) {
    bench_case row = base;
    row.n = n_axis[i];
    add_case(cases, &ncases, row);
  }
  for (int i = 0; i < 4; i++) {
    bench_case row = base;
    row.M = M_axis[i];
    add_case(cases, &ncases, row);
  }
  sort_cases(cases, ncases);
  return ncases;
}

static void warm_cpu(double seconds) {
  if (seconds <= 0) return;
  double x = 1.000001;
  uint64_t start = clock_nsec();
  uint64_t deadline = start + (uint64_t)(seconds*1e9);
  while (clock_nsec() < deadline) {
    for (int i = 0; i < 1000; i++) {
      x += log(x);
      if (x > 2) x = 1.000001;
    }
  }
  consume_double(x);
}

static void alloc_problem(int p, int q, int r, int n, int M, double **A,
                          double **B, double **Sig, double **X) {
  int nA = r*r*(p > 0 ? p : 1);
  int nB = r*r*(q > 0 ? q : 1);
  *A = malloc(sizeof(double)*nA);
  *B = malloc(sizeof(double)*nB);
  *Sig = malloc(sizeof(double)*r*r);
  *X = malloc(sizeof(double)*r*n*M);
  if (!*A || !*B || !*Sig || !*X) die("allocation failed");
}

static double time_case(bench_case c, double target, randompack_rng *rng) {
  int p = c.p;
  int q = c.q;
  int r = c.r;
  int n = c.n == 0 ? (p > q ? p : q) : c.n;
  int M = c.M;
  int icase = 0;
  char name[12] = "rho";
  double rho = 0.95;
  double *A, *B, *Sig, *X;
  varmapack_error error;
  int reps = 0;
  uint64_t start, t;
  alloc_problem(p, q, r, n, M, &A, &B, &Sig, &X);
  error = varmapack_testcase(name, &icase, rho, &p, &q, &r, A, B, Sig, rng);
  if (error) {
    fprintf(stderr, "varmapack_testcase failed: %s\n", varmapack_strerror(error));
    exit(1);
  }
  if (varmapack_specrad(A, r, p) >= 1) die("rho must be below 1");
  if (!randompack_seed(12345, 0, 0, rng)) die("randompack_seed failed");
  start = clock_nsec();
  t = start;
  while ((t - start)*1e-9 < target) {
    error = varmapack_sim(A, B, Sig, 0, 0, p, q, r, n, M, 0, 0, 1, X, 0, rng);
    if (error) {
      fprintf(stderr, "varmapack_sim failed: %s\n", varmapack_strerror(error));
      exit(1);
    }
    consume_double(X[r*n*M - 1]);
    reps++;
    t = clock_nsec();
  }
  free(X);
  free(Sig);
  free(B);
  free(A);
  return (double)(t - start)/((double)reps*r*n*M);
}

static double time_setup_case(bench_case c, double target, randompack_rng *rng) {
  bench_case setup = c;
  setup.n = c.p > c.q ? c.p : c.q;
  setup.M = 1;
  return time_case(setup, target, rng);
}

int main(int argc, char **argv) {
  options opts;
  bool help;
  bench_case cases[MAX_CASES];
  int ncases;
  randompack_rng *rng;
  if (!get_options(argc, argv, &opts, &help) || help) {
    print_help();
    return help ? 0 : 1;
  }
  rng = randompack_create(0);
  if (!rng) die("randompack_create failed");
  printf("Warmup time:     %.2f s\n", opts.w);
  printf("Bench time:      %.3g s per case\n", opts.t);
  printf("Varmapack time:  ns/value\n");
  printf("\n");
  warm_cpu(opts.w);
  ncases = make_cases(&opts, cases);
  printf("%3s %5s %3s %3s %5s %10s %10s\n", "r", "n", "p", "q", "M", "Setup", "C");
  for (int i = 0; i < ncases; i++) {
    bench_case c = cases[i];
    int n = c.n == 0 ? (c.p > c.q ? c.p : c.q) : c.n;
    int h = c.p > c.q ? c.p : c.q;
    double init = time_setup_case(c, opts.t, rng)*h/(n*c.M);
    double ns = time_case(c, opts.t, rng);
    printf("%3d %5d %3d %3d %5d %10.2f %10.1f%s\n", c.r, n, c.p, c.q, c.M, init, ns,
           same_case(c, opts.ref) ? " reference" : "");
  }
  randompack_free(rng);
  return 0;
}
