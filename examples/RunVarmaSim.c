#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>

#include "error.h"
#include "randompack.h"
#include "printX.h"
#include "varmapack.h"
#include "getopt.h"
//#include "ExtraUtil.h"
#include "VarmaUtilities.h"

#define TCN 20  // Max chars in testcase argument (not counting \0) 

static void print_testcase_table(void) {
  int p, q, r, icase, ncases;
  varmapack_error error;
  error = varmapack_testcase(0, 0, 0, "max", &p, &q, &r, &ncases, 0, 0);
  xAssert(error == VARMAPACK_OK);
  char name[32];
  printf("No. Name          p  q  r\n");
  for(icase = 1; icase <= ncases; icase++) {
    name[0] = 0;
    error = varmapack_testcase(0, 0, 0, name, &p, &q, &r, &icase, 0, 0);
    xAssert(error == VARMAPACK_OK);
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
  printf("  -e seed    RNG seed (default randomized)\n");
  printf("  -p         Print generated series (otherwise only a summary)\n\n");
  printf("Default RNG: Randompack default engine.\n");
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

static bool get_options(int argc, char **argv, char *testcase, bool *print,
                        bool *help, int *n, int *seed) {
  opterr = 0;
  optind = 1;
  int opt;
  *print = *help = false;
  *seed = -1; // to indicate no seed
  *n = 20;

  while ((opt = getopt(argc, argv, "hn:e:p")) != -1) {
    switch (opt) {
      case 'h':
        *help = true;
        return true;
      case 'n':
        *n = atoi(optarg);
        if (*n <= 0)
          FAIL("Number of terms must be positive");
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
    varmapack_error error;
    char name[16] = "";
    if (!strcmp(s, "")) FAIL("Empty argument");  
    char MAX[TCN + 1] = "max";
    error = varmapack_testcase(0, 0, 0, MAX, &P, &Q, &R, &ncase, 0, 0);
    if (error != VARMAPACK_OK) FAIL("Internal error");
    if (all_digits(s)) { // Pure digits (indexed testcase) → ask varmapack_testcase() for dimensions
      *icase = atoi(s);
      if (*icase < 1 || *icase > ncase) FAIL("Illegal testcase index");
      error = varmapack_testcase(0, 0, 0, name, p, q, r, icase, 0, 0);
    }
    else { // Named testcase → ask varmapack_testcase() for dimensions
      error = varmapack_testcase(0, 0, 0, s, p, q, r, icase, 0, 0);
    }
    if (error != VARMAPACK_OK) return false;
  }
  return true;
}

int main(int argc, char **argv) {
  bool ok;
  varmapack_error error;
  int n = 0, p, q, r, icase, seed;
  char testcase[TCN + 1];
  double *A = 0, *B = 0, *Sig = 0, *X = 0, *mu = 0, *Gamma = 0;
  randompack_rng *rng = 0;
  bool print, help;
  ok = get_options(argc, argv, testcase, &print, &help, &n, &seed);
  if (help || !ok) {
    print_help();
    return ok ? 0 : 1;
  }
  if (!testcase_dims(testcase, &p, &q, &r, &icase)) return 1;
  if (!ALLOC(A, (p > 0 ? p : 1)*r*r)) goto fail;
  if (!ALLOC(B, (q > 0 ? q : 1)*r*r)) goto fail;
  if (!ALLOC(Sig, r*r)) goto fail;
  if (!ALLOC(X, r*n)) goto fail;
  rng = randompack_create(0);
  if (!rng) goto fail;
  if (seed < 0)
    randompack_randomize(rng);
  else
    randompack_seed(seed, 0, 0, rng);
  char name[16] = "";
  error = varmapack_testcase(A, B, Sig, name, &p, &q, &r, &icase, 0, rng);
  if (error != VARMAPACK_OK) goto fail;
  //
  error = varmapack_sim(A, B, Sig, 0, 0, p, q, r, n, 1, 0, 0, 1, X, 0, rng);
  if (error != VARMAPACK_OK) goto fail;
  //
  printM("X", X, r, n);
  if (!ALLOC(mu, r)) goto fail;
  if (!ALLOC(Gamma, 2*r*r)) goto fail;
  meanmat("T", r, n, X, r,  mu);
  error = varmapack_acvf(A, B, Sig, p, q, r, Gamma, 1);
  if (error != VARMAPACK_OK) goto fail;
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
  FREE(A);
  FREE(B);
  FREE(Sig);
  FREE(X);
  FREE(mu);
  FREE(Gamma);
  randompack_free(rng);
  return 0;
fail:
  FREE(A);
  FREE(B);
  FREE(Sig);
  FREE(X);
  FREE(mu);
  FREE(Gamma);
  randompack_free(rng);
  return 1;
}
