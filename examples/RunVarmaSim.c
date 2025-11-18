#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>

#include "allocate.h"
#include "RandomNumbers.h"
#include "printX.h"
#include "varmapack.h"
#include "varmapack.h"
#include "getopt.h"
//#include "ExtraUtil.h"
#include "VarmaUtilities.h"
#include "varmapack.h"
#include "xAssert.h"

#define TCN 20  // Max chars in testcase argument (not counting \0) 

static void print_testcase_table(void) {
  int p, q, r, icase, ncases;
  bool ok;
  varmapack_testcase(0, 0, 0, "max", &p, &q, &r, &ncases, 0, 0); // how many?
  char name[32];
  printf("No. Name          p  q  r\n");
  for(icase = 1; icase <= ncases; icase++) {
    name[0] = 0;
    ok = varmapack_testcase(0, 0, 0, name, &p, &q, &r, &icase, 0, stdout);
    xAssert(ok);
    printf("%2d  %-12s %2d %2d %2d\n", icase, name, p, q, r);
  }
}

static void print_help(void) {
  printf("Runvarmapack_sim — simulate spin-up-freem VAR/VARMA time series\n");
  printf("Usage: Runvarmapack_sim testcase [options]\n\n");
  printf("  testcase may be a named testcase (e.g. smallMA), a number (1–12),\n");
  printf("  or dimensions p,q,r.\n\n");

  printf("Options:\n");
  printf("  -h         Show this help message\n");
  printf("  -n number  Number of terms to generate (default 20)\n");
  printf("  -P         Use Park–Miller RNG (default: Xorshift128+)\n");
  printf("  -e seed    RNG seed (default randomized, but 42 for Park-Miller)\n");
  printf("  -p         Print generated series (otherwise only a summary)\n\n");
  printf("Default RNG: Xorshift128+.\n");
  printf("To randomize with -P set seed to -1\n");
  printf("Summary fields: length of series, n\n");
  printf("                one-step correlation matrix Corr(x(t), x(t-1))\n");
  printf("                run-time (N/A in this build).\n\n");
  printf("Available named testcases:\n");
  print_testcase_table();
}

#define FAIL(...) do { fprintf(stderr, __VA_ARGS__); \
                       fputc('\n', stderr); return false; } while (0)

static bool all_digits(const char *s) {
  if (*s == 0) return false;
  for (const char *p = s; *p; ++p)
    if (!isdigit((unsigned char)*p)) return false;
  return true;
}

static bool parse_dims_token(const char *s, int *p, int *q, int *r) {
  // returns true if ok, else false; expects "p,q,r" p, q ≥ 0 and r ≥ 1
  int result;
  result = sscanf(s, "%d,%d,%d", p, q, r);
  if (result != 3) FAIL("Wrong number of testcase dimensions");
  if (*p < 0 || *q < 0 || *r < 1) FAIL("Illegal testcase dimensions");
  return true;
}

static bool get_options(int argc, char **argv, char *testcase, bool *print, bool
                        *ParkMiller, bool *help, int *n, int *seed) {
  opterr = 0;
  optind = 1;
  int opt;
  *print = *ParkMiller = *help = false;
  *seed = -1; // to indicate no seed
  *n = 20;

  while ((opt = getopt(argc, argv, "hn:Pe:p")) != -1) {
    switch (opt) {
      case 'h':
        *help = true;
        return true;
      case 'n':
        *n = atoi(optarg);
        if (*n <= 0)
          FAIL("Number of terms must be positive");
        break;
      case 'P':
        *ParkMiller = true;
        break;
      case 'e':
        *seed = atoi(optarg);
        break;
      case 'p':
        *print = true;
        break;
      case '?':
      default:
        FAIL("Illegal option: %c", optopt);
    }
  }
  if (*ParkMiller && *seed==-1) *seed = 42;
  // Copy testcase
  if (optind >= argc)
    FAIL("Missing testcase");
  if (optind + 1 < argc)
    FAIL("Too many arguments");
  int tclen = snprintf(testcase, TCN + 1, "%s", argv[optind]);
  if (tclen > TCN) FAIL("varmapack_testcase argument, too many characters");
  return true;
}

static bool testcase_dims(char *s, int *p, int *q, int *r, int *icase) {
  // Determine testcase type and dimensions.
  // Returns true on success, false on error.
  int ncase, P, Q, R;
  if (!s || !p || !q || !r || !icase) FAIL("Internal error");
  if (strchr(s, ',')) { // Comma-separated dimensions "p,q,r" */
    if (!parse_dims_token(s, p, q, r)) return false;    
    *icase = 0;
  } 
  else {
    bool ok;
    char name[16] = "";
    if (!strcmp(s, "")) FAIL("Empty argument");  
    char MAX[TCN + 1] = "max";
    if (!varmapack_testcase(0, 0, 0, MAX, &P, &Q, &R, &ncase, 0, 0)) FAIL("Internal error");
    if (all_digits(s)) { // Pure digits (indexed testcase) → ask varmapack_testcase() for dimensions
      *icase = atoi(s);
      if (*icase < 1 || *icase > ncase) FAIL("Illegal testcase index");
      ok = varmapack_testcase(0, 0, 0, name, p, q, r, icase, 0, stderr);
    }
    else { // Named testcase → ask varmapack_testcase() for dimensions
      ok = varmapack_testcase(0, 0, 0, s, p, q, r, icase, 0, stderr);
    }
    if (!ok) return false;
  }
  return true;
}

int main(int argc, char **argv) {
  bool ok, vsok;
  int n = 0, p, q, r, icase, seed;
  char testcase[TCN + 1];
  double *A, *B, *Sig, *X;
  bool print, ParkMiller, help;
  ok = get_options(argc, argv, testcase, &print, &ParkMiller, &help, &n, &seed);
  if (help || !ok) {
    print_help();
    return ok ? 0 : 1;
  }
  if (!testcase_dims(testcase, &p, &q, &r, &icase)) return 1;
  allocate(A, p*r*r);
  allocate(B, q*r*r);
  allocate(Sig, r*r);
  allocate(X, r*n);
  const char *rngtype = ParkMiller ? "PM" : "Xorshift";
  int rngseed = (seed < 0) ? 0 : seed;
  randompack_rng *rng = randompack_create(rngtype, rngseed);
  char name[16] = "";
  if (!varmapack_testcase(A, B, Sig, name, &p, &q, &r, &icase, rng, stderr)) return 1;
  //
  varmapack_sim(A, B, Sig, 0, p, q, r, n, 1, 0, 0, rng, X, 0, &vsok);
  //
  printM("X", X, r, n);
  double *mu, *Gamma;
  allocate(mu, r);
  allocate(Gamma, r*r);
  meanmat("T", r, n, X, r,  mu);
  varmapack_acvf(A, B, Sig, p, q, r, Gamma, 1);
  if (print) {
    print4I("p, q, r, n", p, q, r, n);
    printMsg("Model definition matrices:");
    printM("A",  A,  r, r * p);
    printM("B",  B,  r, r * q);
    printM("Sig", Sig, r, r);
    printM("Gamma", Gamma, r, 2*r);
    printM("X", X, r, n);
    printMT("X", X, r, n);
  }
  else
    printf("This will later print a summary\n");
  freem(A);
  freem(B);
  freem(Sig);
  freem(X);
  freem(mu);
  randompack_free(rng);
  return 0;
}
