#define _GNU_SOURCE
#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* Your project headers / APIs */
#include "BlasGateway.h"   /* printM prototype lives here */
#include "RandomNumbers.h" /* rand_rng (opaque); we pass NULL for now */
#include "allocate.h"      /* allocate(ptr, count) */

/* Opaque RNG forward (if not included via RandomNumbers.h) */
typedef struct rand_rng rand_rng;

/* Library prototypes (as per your docs) */
extern void Testcase(double* A, double* B, double* Sig, char* name, int* p, int* q, int* r,
                     int icase, rand_rng* rng);

/* printM: print matrix in column-major order (nr rows, nc cols) */
extern void printM(char* name, double A[], int nr, int nc);

static void print_testcase_table(void)
{
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

static void print_help(void)
{
  printf("RunVarmaSim — simulate spin-up-free VAR/VARMA time series\n");
  printf("Usage: RunVarmaSim testcase [options]\n\n");
  printf("  testcase may be a named testcase (e.g. smallMA), a number (1–12),\n");
  printf("  or dimensions p,q,r.\n\n");

  printf("Options:\n");
  printf("  -h                 Show this help message\n");
  printf("  -n number          Number of terms to generate (default 1000)\n");
  printf("  -P                 Use Park–Miller random number generator\n");
  printf("                     (default: Xorshift128+)\n");
  printf("  -e seed            RNG seed (default: 42)\n");
  printf("  -p                 Print generated series (otherwise only a summary)\n\n");

  printf("Default RNG: Xorshift128+.\n");
  printf("Summary fields: n, one-step correlation matrix Corr[x(t), x(t-1)],\n");
  printf("and run-time (N/A in this build).\n\n");

  printf("Available named testcases:\n");
  print_testcase_table();
}

/* --- helpers to interpret the testcase token --- */

static int all_digits(const char* s)
{
  if (*s == 0)
    return 0;
  for (const char* p = s; *p; ++p)
    if (!isdigit((unsigned char)*p))
      return 0;
  return 1;
}

static int parse_dims_token(const char* s, int* p, int* q, int* r)
{
  /* returns 1 if ok, else 0; expects "p,q,r" with p,q>=0, r>=1 */
  int P, Q, R;
  if (sscanf(s, "%d,%d,%d", &P, &Q, &R) == 3 && P >= 0 && Q >= 0 && R >= 1) {
    *p = P;
    *q = Q;
    *r = R;
    return 1;
  }
  return 0;
}

static int name_to_icase(const char* s)
{
  static const char* names[] = {"tinyAR",      "tinyMA",      "tinyARMA",   "smallAR",
                                "smallMA",     "smallARMA1",  "smallARMA2", "mediumAR",
                                "mediumARMA1", "mediumARMA2", "mediumMA",   "largeAR"};
  for (int i = 0; i < 12; ++i)
    if (strcmp(s, names[i]) == 0)
      return i + 1;
  return 0;
}

/* Parse options; returns 0 on success, nonzero on error. Exits on -h. */
static int get_options(int argc, char** argv, char* out_testcase, size_t out_tsize, int* out_n)
{
  if (argc < 2) {
    fprintf(stderr, "Error: missing testcase argument.\n");
    print_help();
    return 1;
  }

  strncpy(out_testcase, argv[1], out_tsize - 1);
  out_testcase[out_tsize - 1] = 0;

  int n = 1000; /* default number of terms */
  int opt;
  opterr = 0;
  optind = 2; /* skip the testcase argument */

  while ((opt = getopt(argc, argv, "hn:Pe:p")) != -1) {
    switch (opt) {
    case 'h':
      print_help();
      exit(0);
    case 'n': {
      int tmp = atoi(optarg);
      if (tmp <= 0) {
        fprintf(stderr, "Error: number of terms must be positive.\n");
        return 1;
      }
      n = tmp;
      break;
    }
    case 'P':
    case 'e':
    case 'p':
      /* reserved for later */
      break;
    case '?':
    default:
      print_help();
      return 1;
    }
  }

  if (optind < argc) {
    fprintf(stderr, "Unexpected argument: %s\n", argv[optind]);
    return 1;
  }

  *out_n = n;
  return 0;
}

int main(int argc, char** argv)
{
  int n = 0;
  char tc_token[128] = "";

  if (get_options(argc, argv, tc_token, sizeof(tc_token), &n) != 0)
    return 1;

  /* Decide icase and/or dimensions (use asserts for validity). */
  int p = 0, q = 0, r = 0;
  int icase = -1;
  char testname[32] = ""; /* Testcase() expects >= 12; 32 is fine */

  if (strchr(tc_token, ',')) {
    /* Dimensions p,q,r provided: random model (icase=0). */
    int ok = parse_dims_token(tc_token, &p, &q, &r);
    assert(ok && "Malformed dimensions: use p,q,r with p,q>=0 and r>=1");
    icase = 0;
  }
  else if (all_digits(tc_token)) {
    /* Numbered named testcase. */
    icase = atoi(tc_token);
    assert(icase >= 1 && icase <= 12 && "Testcase number must be 1..12");
    /* Enquiry call to get p,q,r (A=B=Sig=NULL). */
    Testcase(0, 0, 0, testname, &p, &q, &r, icase, 0);
  }
  else {
    /* Named testcase by string. */
    icase = name_to_icase(tc_token);
    assert(icase >= 1 && icase <= 12 && "Unknown testcase name");
    Testcase(0, 0, 0, testname, &p, &q, &r, icase, 0);
  }

  printf("Running testcase: %s\n", (icase >= 1 ? testname : tc_token));
  printf("Dims: p=%d q=%d r=%d, n=%d\n", p, q, r, n);

  /* Allocate arrays (your allocate() handles the memory). */
  double *A = 0, *B = 0, *Sig = 0, *mu = 0;
  allocate(A, (size_t)r * (size_t)r * (size_t)p);
  allocate(B, (size_t)r * (size_t)r * (size_t)q);
  allocate(Sig, (size_t)r * (size_t)r);
  allocate(mu, (size_t)r);

  /* Optional: initialize mu (not used further here). */
  for (int i = 0; i < r; ++i)
    mu[i] = (double)(i + 1);

  /* Fill parameter matrices via Testcase. */
  Testcase(A, B, Sig, testname, &p, &q, &r, icase, 0);

  /* Print A (r x r*p), B (r x r*q), and Sig (r x r). */
  printM("A", A, r, r * p);
  printM("B", B, r, r * q);
  printM("Sig", Sig, r, r);

  /* Cleanup: allocate() pairs with free(). */
  free(A);
  free(B);
  free(Sig);
  free(mu);
  return 0;
}
