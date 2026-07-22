// TimeBreakEven.c: time VYW vs Lyapunov covariance solvers only.

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "getopt.h"
#include "Lyapunov.h"
#include "varmapack.h"
#include "VYW.h"

typedef struct {
  double t;
  double w;
  int r0;
  int r1;
  int p0;
  int p1;
  int q0;
  int q1;
} options;

typedef struct {
  double vyw_ns;
  double lyap_ns;
} timing;

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
    .t = 0.01, .w = 0.1, .r0 = 7, .r1 = 30, .p0 = 1, .p1 = 7, .q0 = 0, .q1 = 7
  };
}

static void print_help(void) {
  printf("TimeBreakEven -- time VYW and Lyapunov solvers only\n");
  printf("Usage: TimeBreakEven [options]\n\n");
  printf("Options:\n");
  printf("  -h          show this help\n");
  printf("  -t seconds  timing target per case (default 0.01)\n");
  printf("  -w seconds  CPU warmup time before timing (default 0.1)\n");
}

static bool get_options(int argc, char **argv, options *opts, bool *help) {
  int opt;
  default_options(opts);
  *help = false;
  opterr = 0;
  optind = 1;
  while ((opt = getopt(argc, argv, "ht:w:")) != -1) {
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
      default:
        return false;
    }
  }
  return optind == argc;
}

static void warm_cpu(double seconds) {
  double x = 1.000001;
  if (seconds <= 0) return;
  uint64_t deadline = clock_nsec() + (uint64_t)(seconds*1e9);
  while (clock_nsec() < deadline) {
    for (int i=0; i<1000; i++) {
      x += log(x);
      if (x > 2) x = 1.000001;
    }
  }
  consume_double(x);
}

static void make_problem(int p, int q, int r, double **A, double **B,
                         double **Sig, randompack_rng *rng) {
  int icase = 0;
  char name[12] = "rho";
  double rho = p > 0 ? 0.95 : 0;
  varmapack_error error;
  *A = malloc(sizeof(double)*r*r*(p > 0 ? p : 1));
  *B = malloc(sizeof(double)*r*r*(q > 0 ? q : 1));
  *Sig = malloc(sizeof(double)*r*r);
  if (!*A || !*B || !*Sig) die("allocation failed");
  error = varmapack_testcase(name, &icase, rho, &p, &q, &r, *A, *B, *Sig, rng);
  if (error) die("varmapack_testcase failed");
}

static timing time_case(int p, int q, int r, double target, randompack_rng *rng) {
  int rr = r*r;
  int reps = 0;
  uint64_t start, now, t0, t1, t2;
  uint64_t vyw_ns = 0;
  uint64_t lyap_ns = 0;
  double *A = 0;
  double *B = 0;
  double *Sig = 0;
  double *S = 0;
  double *C = 0;
  double *G = 0;
  make_problem(p, q, r, &A, &B, &Sig, rng);
  S = malloc(sizeof(double)*rr*(p+1));
  C = malloc(sizeof(double)*rr*(q+1));
  G = malloc(sizeof(double)*rr*(q+1));
  if (!S || !C || !G) die("allocation failed");
  start = clock_nsec();
  now = start;
  while ((now - start)*1e-9 < target) {
    t0 = clock_nsec();
    if (!VYWFactorizeSolve(A, B, Sig, p, q, r, S, C, G)) die("VYW failed");
    t1 = clock_nsec();
    if (!LyapunovFactorizeSolve(A, B, Sig, p, q, r, S, C, G)) {
      die("Lyapunov failed");
    }
    t2 = clock_nsec();
    vyw_ns += t1 - t0;
    lyap_ns += t2 - t1;
    consume_double(S[0]);
    reps++;
    now = clock_nsec();
  }
  free(G);
  free(C);
  free(S);
  free(Sig);
  free(B);
  free(A);
  return (timing){(double)vyw_ns/(double)reps, (double)lyap_ns/(double)reps};
}

int main(int argc, char **argv) {
  options opts;
  bool help;
  randompack_rng *rng;
  if (!get_options(argc, argv, &opts, &help) || help) {
    print_help();
    return help ? 0 : 1;
  }
  rng = randompack_create(0);
  if (!rng) die("randompack_create failed");
  printf("Warmup time: %.2f s\n", opts.w);
  printf("Bench time:  %.3g s per case\n\n", opts.t);
  warm_cpu(opts.w);
  printf(" r  p  q  nVYW  nLyap     VYW_us    Lyap_us  Lyap/VYW\n");
  for (int p=opts.p0; p<=opts.p1; p++) {
    int r1p = 32 - 2*p;
    if (r1p > opts.r1) r1p = opts.r1;
    for (int r=opts.r0; r<=r1p; r++) {
      for (int q=opts.q0; q<=opts.q1; q++) {
        int nvyw = p == 0 ? 0 : r*r*p - r*(r-1)/2;
        int nx = p > 0 ? p + 1 : 1;
        int nlyap = r*(nx + q + 1);
        timing tm = time_case(p, q, r, opts.t, rng);
        printf("%2d %2d %2d %5d %6d %10.1f %10.1f %9.3f\n",
               r, p, q, nvyw, nlyap, tm.vyw_ns/1000, tm.lyap_ns/1000,
               tm.lyap_ns/tm.vyw_ns);
      }
    }
  }
  randompack_free(rng);
  return 0;
}
