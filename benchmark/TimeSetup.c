// TimeSetup.c: time the current VYW plus SBuild setup.

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "getopt.h"
#include "BlasGateway.h"
#include "Lyapunov.h"
#include "VarmaPackUtil.h"
#include "VarmaUtilities.h"
#include "varmapack.h"
#include "VYW.h"

enum { MAX_VALUES = 32, MAX_CASES = 128 };

typedef struct {
  int n;
  int v[MAX_VALUES];
} int_list;

typedef struct {
  int p;
  int q;
  int r;
} setup_case;

typedef struct {
  double t;
  double w;
  int d;
  int_list p;
  int_list r;
} options;

typedef struct {
  double vyw_ns;
  double lyap_ns;
  double other_vyw_ns;
  double other_lyap_ns;
} time_result;

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

static bool parse_int_list(char *s, int_list *xs) {
  char *p = s;
  char *end;
  xs->n = 0;
  while (*p) {
    if (xs->n >= MAX_VALUES) return false;
    xs->v[xs->n] = (int)strtol(p, &end, 10);
    if (end == p) return false;
    xs->n++;
    if (*end == 0) return true;
    if (*end != ',') return false;
    p = end + 1;
  }
  return xs->n > 0;
}

static void default_options(options *opts) {
  *opts = (options) {
    .t = 0.2,
    .w = 0.1,
    .d = 2,
    .p = {3, {1, 3, 5}},
    .r = {4, {2, 5, 12, 32}}
  };
}

static void print_help(void) {
  printf("TimeSetup -- time VYWFactorizeSolve + SBuild (us/setup)\n");
  printf("Usage: TimeSetup [options]\n\n");
  printf("Options:\n");
  printf("  -h          show this help\n");
  printf("  -t seconds  timing target per case (default 0.2)\n");
  printf("  -w seconds  CPU warmup time before timing (default 0.1)\n");
  printf("  -d digits   printed digits (default 2)\n");
  printf("  -p list     p values (default 1,3,5)\n");
}

static bool get_options(int argc, char **argv, options *opts, bool *help) {
  int opt;
  default_options(opts);
  *help = false;
  opterr = 0;
  optind = 1;
  while ((opt = getopt(argc, argv, "ht:w:d:p:")) != -1) {
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
        if (!parse_int_list(optarg, &opts->p)) return false;
        break;
      default:
        return false;
    }
  }
  return optind == argc;
}

static bool same_case(setup_case a, setup_case b) {
  return a.p == b.p && a.q == b.q && a.r == b.r;
}

static void add_case(setup_case cases[], int *ncases, setup_case row) {
  for (int i = 0; i < *ncases; i++) {
    if (same_case(cases[i], row)) return;
  }
  if (*ncases >= MAX_CASES) die("too many benchmark cases");
  cases[*ncases] = row;
  (*ncases)++;
}

static bool case_less(setup_case a, setup_case b) {
  if (a.r != b.r) return a.r < b.r;
  if (a.p != b.p) return a.p < b.p;
  return a.q < b.q;
}

static void sort_cases(setup_case cases[], int ncases) {
  for (int i = 1; i < ncases; i++) {
    setup_case row = cases[i];
    int j = i;
    while (j > 0 && case_less(row, cases[j-1])) {
      cases[j] = cases[j-1];
      j--;
    }
    cases[j] = row;
  }
}

static int make_cases(options *opts, setup_case cases[]) {
  int ncases = 0;
  for (int i = 0; i < opts->r.n; i++) {
    for (int j = 0; j < opts->p.n; j++) {
      if (opts->r.v[i] == 5) {
        add_case(cases, &ncases, (setup_case){opts->p.v[j], 0, opts->r.v[i]});
        add_case(cases, &ncases, (setup_case){opts->p.v[j], 3, opts->r.v[i]});
        add_case(cases, &ncases, (setup_case){opts->p.v[j], 7, opts->r.v[i]});
      }
      else {
        add_case(cases, &ncases, (setup_case){opts->p.v[j], 3, opts->r.v[i]});
      }
    }
  }
  sort_cases(cases, ncases);
  return ncases;
}

static void warm_cpu(double seconds) {
  double x = 1.000001;
  if (seconds <= 0) return;
  uint64_t deadline = clock_nsec() + (uint64_t)(seconds*1e9);
  while (clock_nsec() < deadline) {
    for (int i = 0; i < 1000; i++) {
      x += log(x);
      if (x > 2) x = 1.000001;
    }
  }
  consume_double(x);
}

static void make_problem(setup_case c, int *p, int *q, int *r, int *h,
                         double **A, double **B, double **Sig,
                         randompack_rng *rng) {
  int icase = 0;
  char name[12] = "rho";
  double rho = 0.95;
  varmapack_error error;
  *p = c.p;
  *q = c.q;
  *r = c.r;
  *h = *p > *q ? *p : *q;
  *A = malloc(sizeof(double)*(*r)*(*r)*(*p > 0 ? *p : 1));
  *B = malloc(sizeof(double)*(*r)*(*r)*(*q > 0 ? *q : 1));
  *Sig = malloc(sizeof(double)*(*r)*(*r));
  if (!*A || !*B || !*Sig) die("allocation failed");
  error = varmapack_testcase(*A, *B, *Sig, name, p, q, r, &icase, rho, rng);
  if (error) {
    fprintf(stderr, "varmapack_testcase failed: %s\n", varmapack_strerror(error));
    exit(1);
  }
}

static bool setup_from_scg(double A[], double B[], double Sig[], int p, int q,
                           int r, int h, double S[], double G[],
                           randompack_rng *rng, double X[]) {
  int rh = r*h;
  double *E = 0;
  double *Psi = 0;
  double *PsiHat = 0;
  double *R = 0;
  double *SS = 0;
  double *Wrk = 0;
  bool ok = false;
  if (!ALLOC(E, rh)) goto done;
  if (!ALLOC(Psi, rh*rh)) goto done;
  if (!ALLOC(PsiHat, rh*rh)) goto done;
  if (!ALLOC(R, rh*rh)) goto done;
  if (!ALLOC(SS, rh*rh)) goto done;
  if (!ALLOC(Wrk, rh)) goto done;
  if (!randompack_mvn("T", 0, Sig, r, h, E, r, 0, rng)) goto done;
  if (!SBuild("Low", S, A, G, p, q, r, h, SS)) goto done;
  FindPsi(A, B, Psi, p, q, r);
  FindPsiHat(Psi, PsiHat, Sig, r, h);
  lacpy("Low", rh, rh, SS, rh, R, rh);
  syrk("Low", "NoT", rh, rh, -1, PsiHat, rh, 1, R, rh);
  if (!randompack_mvn("T", 0, R, rh, 1, Wrk, rh, 0, rng)) goto done;
  lacpy("All", rh, 1, Wrk, rh, X, rh);
  gemm("NoT", "NoT", rh, 1, rh, 1, Psi, rh, E, rh, 1, X, rh);
  ok = true;
done:
  FREE(Wrk);
  FREE(SS);
  FREE(R);
  FREE(PsiHat);
  FREE(Psi);
  FREE(E);
  return ok;
}

static time_result time_case(setup_case c, double target, randompack_rng *rng) {
  int p = c.p;
  int q = c.q;
  int r = c.r;
  int h;
  double *A, *B, *Sig, *S, *C, *G, *X;
  int reps = 0;
  uint64_t start, t, t0, t1, t2, t3, t4;
  uint64_t vyw_ns = 0;
  uint64_t lyap_ns = 0;
  uint64_t other_vyw_ns = 0;
  uint64_t other_lyap_ns = 0;
  make_problem(c, &p, &q, &r, &h, &A, &B, &Sig, rng);
  S = malloc(sizeof(double)*r*r*(p+1));
  C = malloc(sizeof(double)*r*r*(q+1));
  G = malloc(sizeof(double)*r*r*(q+1));
  X = malloc(sizeof(double)*r*h);
  if (!S || !C || !G || !X) die("allocation failed");
  if (!randompack_seed(12345, 0, 0, rng)) die("randompack_seed failed");
  start = clock_nsec();
  t = start;
  while ((t - start)*1e-9 < target) {
    t0 = clock_nsec();
    if (!VYWFactorizeSolve(A, B, Sig, p, q, r, S, C, G)) {
      die("VYWFactorizeSolve failed");
    }
    t1 = clock_nsec();
    if (!setup_from_scg(A, B, Sig, p, q, r, h, S, G, rng, X)) {
      die("VYW post-solver setup failed");
    }
    t2 = clock_nsec();
    if (!LyapunovFactorizeSolve(A, B, Sig, p, q, r, S, C, G)) {
      die("LyapunovFactorizeSolve failed");
    }
    t3 = clock_nsec();
    if (!setup_from_scg(A, B, Sig, p, q, r, h, S, G, rng, X)) {
      die("Lyapunov post-solver setup failed");
    }
    t4 = clock_nsec();
    vyw_ns += t1 - t0;
    other_vyw_ns += t2 - t1;
    lyap_ns += t3 - t2;
    other_lyap_ns += t4 - t3;
    consume_double(X[r*h - 1]);
    reps++;
    t = clock_nsec();
  }
  time_result result = {
    .vyw_ns = (double)vyw_ns/(double)reps,
    .lyap_ns = (double)lyap_ns/(double)reps,
    .other_vyw_ns = (double)other_vyw_ns/(double)reps,
    .other_lyap_ns = (double)other_lyap_ns/(double)reps
  };
  free(X);
  free(G);
  free(C);
  free(S);
  free(Sig);
  free(B);
  free(A);
  return result;
}

int main(int argc, char **argv) {
  options opts;
  bool help;
  setup_case cases[MAX_CASES];
  int ncases;
  randompack_rng *rng;
  if (!get_options(argc, argv, &opts, &help) || help) {
    print_help();
    return help ? 0 : 1;
  }
  rng = randompack_create(0);
  if (!rng) die("randompack_create failed");
  printf("Warmup time %.2f s\n", opts.w);
  printf("Bench time %.3g s per case\n", opts.t);
  printf("There are r·h values set up for each case\n");
  printf("The table shows setup time in µs except r = 32 which shows ms\n");
  printf("\n");
  warm_cpu(opts.w);
  ncases = make_cases(&opts, cases);
  printf("–– Dimensions ––    ––––––– Time per case ––––––––"
         "    ––– Time per value –––\n");
  printf("   r   p   q   h     VYW   SLICOT    Other   VTotal"
         "   STotal     VYW   SLICOT\n");
  for (int i = 0; i < ncases; i++) {
    setup_case c = cases[i];
    int h = c.p > c.q ? c.p : c.q;
    time_result t = time_case(c, opts.t, rng);
    double other_ns = 0.5*(t.other_vyw_ns + t.other_lyap_ns);
    double vtotal_ns = t.vyw_ns + t.other_vyw_ns;
    double stotal_ns = t.lyap_ns + t.other_lyap_ns;
    double sc = c.r < 20 ? 1000 : 1000000;
    double scale = c.r*h*sc;
    printf("%4d %3d %3d %3d %7.1f %8.1f %8.1f %8.1f %8.1f %7.1f %8.1f\n",
           c.r, c.p, c.q, h, t.vyw_ns/sc, t.lyap_ns/sc, other_ns/sc,
           vtotal_ns/sc, stotal_ns/sc, t.vyw_ns/scale, t.lyap_ns/scale);
  }
  randompack_free(rng);
  return 0;
}
