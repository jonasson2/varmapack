#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>

#include "allocate.h"
#include "RandomNumbers.h"
#include "printX.h"
#include "Testcase.h"
#include "VarmaSim.h"
#include "getopt.h"
#include "ExtraUtil.h"
#include "VarmaUtilities.h"
#include "ACVF.h"

#define TCN 20  // Max chars in testcase argument (not counting \0) 

static void print_testcase_table(void) {
  printf("No. Name          p  q  r\n");
  printf(" 1  tinyAR        1  0  1\n");
  printf(" 2  tinyMA        0  1  1\n");
  printf(" 3  tinyARMA      1  1  1\n");
  printf(" 4  smallAR       1  0  2\n");
  printf(" 5  smallMA       0  2  2\n");
  printf(" 6  smallARMA1    1  1  2\n");
  printf(" 7  smallARMA2    1  2  2\n");
  printf(" 8  mediumAR      1  0  3\n");
  printf(" 9  mediumARMA1   3  3  3\n");
  printf("10  mediumARMA2   3  3  3\n");
  printf("11  mediumMA      0  2  3\n");
  printf("12  largeAR       5  0  7\n");
}

static void print_help(void) {
  printf("RunVarmaSim — simulate spin-up-free VAR/VARMA time series\n");
  printf("Usage: RunVarmaSim testcase [options]\n\n");
  printf("  testcase may be a named testcase (e.g. smallMA), a number (1–12),\n");
  printf("  or dimensions p,q,r.\n\n");

  printf("Options:\n");
  printf("  -h         Show this help message\n");
  printf("  -n number  Number of terms to generate (default 1000)\n");
  printf("  -P         Use Park–Miller random number generator\n");
  printf("             (default: Xorshift128+)\n");
  printf("  -e seed    RNG seed (default randomized, but 42 for Park-Miller)\n");
  printf("  -p         Print generated series (otherwise only a summary)\n\n");

  printf("Default RNG: Xorshift128+.\n");
  printf("To randomize with -P set seed to -1");
  printf("Summary fields: n, one-step correlation matrix Corr[x(t), x(t-1)],\n");
  printf("and run-time (N/A in this build).\n\n");

  printf("Available named testcases:\n");
  print_testcase_table();
}

/* --- helpers to interpret the testcase token --- */

#define FAIL(...) do { fprintf(stderr, __VA_ARGS__); \
                       fputc('\n', stderr); return false; } while (0)
static bool all_digits(const char *s) {
  if (*s == 0) return false;
  for (const char *p = s; *p; ++p)
    if (!isdigit((unsigned char)*p)) return false;
  return true;
}

static bool parse_dims_token(const char *s, int *p, int *q, int *r) {
  // returns true if ok, else false; expects "p,q,r" with 0 <= p,q <= 12, 1 <= r <= 7
  int result = sscanf(s, "%d,%d,%d", p, q, r);
  if (result != 3) FAIL("Wrong number of testcase dimensions");
  if (*p < 0 || *q < 0 || *r < 1) FAIL("Illegal testcase dimensions");
  if (*p > 12 || *q > 12 || *r > 7) FAIL("Illegal testcase dimensions");
  return true;
}

// Parse options; returns true on success, false on error
static bool get_options(int argc, char **argv,
                        char *testcase,
                        bool *print,
                        bool *ParkMiller,
                        bool *help,
                        int *n,
                        int *seed) {
  
  opterr = 0;
  optind = 1;
  int opt;
  *print = *ParkMiller = *help = false;
  *seed = -1; // to indicate no seed
  *n = 1000;

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
  if (tclen > TCN) FAIL("Testcase argument, too many characters");
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
    if (!strcmp(s, "")) FAIL("Empty argument");  
    char MAX[TCN + 1] = "max";
    if (!Testcase(0, 0, 0, MAX, &P, &Q, &R, &ncase, 0, 0)) FAIL("Internal error");
    if (all_digits(s)) { // Pure digits (indexed testcase) → ask Testcase() for dimensions
      *icase = atoi(s);
      if (*icase < 1 || *icase > ncase) FAIL("Illegal testcase index");
      ok = Testcase(0, 0, 0, 0, p, q, r, icase, 0, stderr);
    }
    else { // Named testcase → ask Testcase() for dimensions
      ok = Testcase(0, 0, 0, s, p, q, r, icase, 0, stderr);
    }
    if (!ok) return false;
  }
  return true;
}

int main(int argc, char **argv) {
  bool ok;
  int vsok;
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
  RandRng *rng = RandCreate();  
  if (ParkMiller) RandSetPM(rng);
  if (seed >= 0) RandSeed(seed, rng);
  
  if (!Testcase(A, B, Sig, 0, &p, &q, &r, &icase, rng, stderr)) return 1;

  VarmaSim(A, B, Sig, 0, p, q, r, n, 1, 0, rng, X, 0, &vsok);

  double *mu, *Gamma;
  allocate(mu, r);
  allocate(Gamma, r*r);
  meanmat("T", r, n, X, r,  mu);
  printM("Mean", mu, 1, r);
  ACVF(A, B, Sig, p, q, r, Gamma, 1);
  if (print) {
    printMsg("Model definition matrices:");
    printM("A",  A,  r, r * p);
    printM("B",  B,  r, r * q);
    printM("Sig", Sig, r, r);
    printM("Gamma", Gamma, r, 2*r);
    // printM("X", X, r, n);
  }
  else
    printf("This will later print a summary\n");
  free(A);
  free(B);
  free(Sig);
  free(X);
  free(mu);
  RandFree(rng);
  return 0;
}
