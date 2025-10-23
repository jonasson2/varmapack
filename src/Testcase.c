#include <string.h>
#include <stdbool.h>
#include "BlasGateway.h"
#include "allocate.h"
#include "VYW.h"
#include "VarmaMisc.h"
#include "VarmaUtilities.h"
#include "RandomNumbers.h"
#include "Testcase.h"
#include "printX.h"
#define DEBUG
#include "debugprint.h"

static void flipud(int m, int n, double src[], double dst[]);
static void hilb(double A[], int m);
static void simplerand(double A[], int n, int seed);
static int find_named_case(const char *namev[], const char *name, int Ncase);
static bool error(FILE *errfp, const char *msg);

bool Testcase (  // Create a testcase for VARMA likelihood calculation
  double A[],    // out     r×r×p, autoregressive parameter matrices (or null)
  double B[],    // out     r×r×q, moving average parameter matrices (or null)
  double Sig[],  // out     r×r, covariance of the shock terms eps(t) (or null)
  char *name,    // out     string w. length >= 12, name of testcase, or null for unnamed
  int *p_p,      // in/out  Number of autoregressive terms
  int *q_p,      // in/out  Number of moving avg. terms
  int *r_p,      // in/out  dimension of each x(t)
  int *icase,    // in/out  index of named testcase to create or 0 or -1 to use p,q,r
  RandRng *rng,  // in      random number generator
  FILE *fp)      // in      stream to print errors (or null)
{
  // TESTCASE CREATION
  // The following kinds of testcases, suitable for testing or timing various components
  // of the Varmasim package (or the parent Varma package) may be constructed. A, B and
  // Sig return the parameter matrices of the created case and p, q and r its dimensions.
  //
  // 1. Named or indexed. When name is specified and icase is null the testcase with the
  //    specified name is returned. The case index is put in icase. When name is "" and
  //    icase > 0 the correspending named testcase, including its name, is returned.
  //
  // 3. Random unnamed. When icase == 0, p, q and r should be specified on entry, and then
  //    a random model of the specified dimensions is returned. A and B are set to random
  //    numbers and Sig is as in case 2. The matrices must be large enough.
  //
  // 4. Deterministic unnamed. When icase = -1, p, q, and r should be specified on entry,
  //    and then A returns a matrix with all elements set to 0.5/r/p, B has all elements
  //    set to 1/r/q and Sig returns Hilbert(r) * 0.2*I [with all elements of A equal to
  //    1/r/p the model would be on the stationary / non-stationary boundary]. This model
  //    is suitable for consistent timing of likelihood calculation.
  //
  // INQUIRY
  // 1. Dimensions. To inquire about the dimensions of a testcase, A, B and Sig may be
  //    null and either name or icase may be specified. The dimensions of the
  //    corresponding case are put in p, q and r, and either name or icase are assigned
  //    to. Testcase(0, 0, 0, "smallAR", &p, &q, &r, &icase, 0) will thus set p, q, r,
  //    icase to 1, 0, 2, 4.
  //
  // 2. Count + max dimensions. Testcase(0, 0, 0, "max", &p, &q, &r, &icase, 0) sets p, q,
  //    r to the maximum dimensions over all testcases and icase to the number of named
  //    testcases.

  const char *namev[] = {// nr   p  q  r
    "tinyAR",      // 1    1  0  1
    "tinyMA",      // 2    0  1  1
    "tinyARMA",    // 3    1  1  1
    "smallAR",     // 4    1  0  2
    "smallMA",     // 5    0  2  2
    "smallARMA1",  // 6    1  1  2
    "smallARMA2",  // 7    1  2  2
    "mediumAR",    // 8    1  0  3
    "mediumARMA1", // 9    3  3  3
    "mediumARMA2", // 10   3  3  3
    "mediumMA",    // 11   0  2  3
    "largeAR"      // 12   5  0  7
  };
  int    pv[]      = { 1,  0,  1,  1,  0,  1,  1,  1,  3,   3,   0,   5};
  int    qv[]      = { 0,  1,  1,  0,  2,  1,  2,  0,  3,   3,   2,   0};
  int    rv[]      = { 1,  1,  1,  2,  2,  2,  2,  3,  3,   3,   3,   7};
  //  Following would be used except that cl complains it is nonstandard
  //                   0   1   2   3   4   5   6   7   8    9   10   11
  //double *Av[]   = {A1,  0, A3, A4,  0, A6, A7, A8, A9, A10,   0, A12};
  //double *Bv[]   = { 0, B2, B3,  0, B5, B6, B7,  0, B9, B10, B11,   0};
  //double *Sigv[] = {S1, S1, S1, S2, S3, S3, S3, S4, S4,  S4,  S4,  S6};

  int Ncase = sizeof(namev)/sizeof(char*);  
  
  // SANITY CHECKS:
  if (!p_p || !q_p || !r_p) return error(fp,"Testcase: p, q and r must not be null ptrs");
  if (icase == 0) return error(fp, "icase must not be a null pointer");
  int nnull = !A + !B + !Sig;
  if (nnull != 0 && nnull != 3)
    return error(fp, "If any of A, B, Sig is null, they must all be");
  int p = *p_p, q = *q_p, r = *r_p;
  if (nnull == 0 && (p < 0 || q < 0 || r <= 0)) {
    char msg[40];
    snprintf(msg, sizeof(msg), "invalid dimensions p=%d q=%d r=%d", *p_p, *q_p, *r_p);
    return error(fp, msg);
  }
  if (name && !strcmp(name, "max")) {
    if (nnull != 3) return error(fp, "When name is max A, B, Sig must all be null");
  }
  else if (name && strlen(name) > 0) {
    *icase = find_named_case(namev, name, Ncase);
    if (*icase == 0) return error(fp, "Unknown testcase name");
  }
  else {
    if (*icase < -1 || *icase > Ncase) return error(fp, "icase out of range");
    if (*icase <= 0 && nnull == 3)
      return error(fp, "A, B, Sig must not be null if icase <= 0");
  }

  // MAX INQUIRY
  if (name && !strcmp(name, "max")) {
    *icase = Ncase;
    p = q = r = 0;
    for (int k=0; k<Ncase; k++) {
      p = max(p, pv[k]);
      q = max(q, qv[k]);
      r = max(r, rv[k]);
    }
    *p_p = p;
    *q_p = q;
    *r_p = r;
    return true;
  }      
  // FILL CASE NAME OR NUMBER
  if (name && strlen(name) == 0)
    snprintf(name, 12, "%s", namev[*icase - 1]);
  else if (name)
    *icase = find_named_case(namev, name, Ncase);

  // CONSTRUCT TESTCASE
  double
    A1[] = {0.4},
    A3[] = {0.4},
    A4[] = {0.1, 0.1,   0.1, 0.1},
    A6[] = {0.3, 0.1,   0.4, 0.2},
    A7[] = {0.2, 0.3,   0.2, 0.3},
    A8[] = {
      0.35, 0.15, 0.15,//  0.11, 0.14, 0.07,
      0.25, 0.15, 0.05,//  0.12, 0.15, 0.08,
      0.15, 0.05, 0.01,//  0.13, 0.16, 0.09
    },
    A9[] = {
      0.15, 0.10, 0.05,  0.11, 0.14, 0.17,  0.01, 0.04, 0.06,
      0.16, 0.11, 0.06,  0.12, 0.15, 0.18,  0.02, 0.05, 0.08,
      0.17, 0.12, 0.07,  0.13, 0.16, 0.19,  0.03, 0.06, 0.09
    },
    A10[9*3],
    A12[7*7*5],
    B2[] = {0.5},
    B3[] = {0.4},
    B5[] = {0.3, 0.1, 0.1, 0.1,   0.3, 0.1, 0.1, 0.1},
    B6[] = {0.2, 0.3,   0.2, 0.3},
    B7[] = {0.4, 0.1, 0.2, 0.1,   0.2, 0.1, 0.3, 0.1},
    B9[9*3],
    B10[9*3],
    B11[] =  {
      0.35, 0.25, 0.15,  0.11, 0.14, 0.17,
      0.25, 0.15, 0.05,  0.12, 0.15, 0.18,
      0.15, 0.05, 0.01,  0.13, 0.16, 0.19
    },
    S1[] = {0.8},
    S2[] = {2, 1,   1, 3},
    S3[] = {2, 1,   1, 2},
    S4[] = {2.0, 0.5, 0.0,   0.5, 2.0, 0.5,   0.0, 0.5, 1.0},
    //S5[] = {1.0, 0.0, 0.0,   0.0, 1.0, 0.0,   0.0, 0.0, 1.0},
    S6[7*7];
  double *Av[]   = {0,0,0,0,0,0,0,0,0,0,0,0};
  double *Bv[]   = {0,0,0,0,0,0,0,0,0,0,0,0};
  double *Sigv[]   = {0,0,0,0,0,0,0,0,0,0,0,0};

  int nVYW, i, j, *piv, info, ok;
  double *VYWfactors;
  // Finish defining the starting matrices:
  Av[0]=A1; Av[2]=A3; Av[3]=A4; Av[5]=A6; Av[6]=A7; Av[7]=A8; Av[8]=A9; Av[9]=A10; Av[11]=A12;
  Bv[1]=B2; Bv[2]=B3; Bv[4]=B5; Bv[5]=B6; Bv[6]=B7; Bv[8]=B9; Bv[9]=B10; Bv[10]=B11;
  Sigv[0]=Sigv[1]=Sigv[2]=S1; Sigv[3]=S2; Sigv[4]=Sigv[5]=Sigv[6]=S3;
  Sigv[7]=Sigv[8]=Sigv[9]=Sigv[10]=S4; Sigv[11]=S6;
  flipud(9, 3, A9, B9);
  flipud(9, 3, A9, A10);
  copy(9*3, A9, 1, B10, 1);
  simplerand(A12, 7*7*5, 366); // 366 is seed
  scal(7*7*5, 1.8/5/7, A12, 1);
  hilb(S6, 7);
  if (*icase == -1 || *icase == 0) {
    if (Sig) {
      hilb(Sig, r);
      for (i = 0; i < r; i++) Sig[i + i*r] += 0.2; // add 0.2 to diagonal
    }
  }
  if (*icase == -1) {
    if (A && p>0) for (i=0; i<r*r*p; i++) A[i] = 0.5/(p*r);
    if (B && q>0) for (i=0; i<r*r*q; i++) B[i] = 1.0/(q*r);     
  }
  else if (*icase == 0) {
    if (A && p>0) {
      Rand(A, r*r*p, rng);
      scal(r*r*p, 0.5/(p*r), A, 1); 
      nVYW = r*r*p - r*(r-1)/2;
      allocate(VYWfactors, nVYW*nVYW);
      allocate(piv, nVYW);
      j = 0;
      while (true) {
        VYWFactorize(A, VYWfactors, piv, p, r, &info);
        ok = info == 0 && IsStationary(A, Sig, VYWfactors, piv, p, r);
        if (ok) break;
        scal(r*p, 0.5, A, 1);
        j++;
        xAssert(j < 10);
      }
      freem(piv); freem(VYWfactors);
    }
    if (B && q>0) {
      Rand(B, r*r*q, rng);
      scal(r*r*q, 1.0/(q*r), B, 1);
    }
  }
  else if (1 <= *icase && *icase <= Ncase) {
    p = pv[*icase-1];
    q = qv[*icase-1];
    r = rv[*icase-1];
    for (i=0; i<7; i++) S6[i + i*r] += 0.2; // add 0.2 to diagonal
    if (A && Av[*icase-1]) copytranspose(p*r, r, Av[*icase-1], p*r, A, r);
    if (B && Bv[*icase-1]) copytranspose(q*r, r, Bv[*icase-1], q*r, B, r);
    if (Sig && Sigv[*icase-1]) copy(r*r, Sigv[*icase-1], 1, Sig, 1);
    if (name) strcpy(name, namev[*icase-1]); // name <---namev
  }
  else return error(fp, "Unexpected error");
    //     case {"pivotfailure"} % create almost singular vyw equations
    //       A = {
    //         -0.6250  0.2400    0.2500    0.1250    0.6250    0.3750
    //         0             0    0.1250    0.5000    0.2500    0.1250
    //         0.2500        0    0.1430         0    0.2500    0.2500};
    //       A(1,1) = -0.62999813363195223; %% as singular as can be made
    //       B={};
    //       Sig = eye(3);
  *p_p = p;
  *q_p = q;
  *r_p = r;
  return true;
}

static void flipud(int m, int n, double*src, double*dst) {
  int i, j;
  for (j=0; j<n; j++)
    for (i=0; i<m; i++) 
      dst[i + j*m] = src[m-1-i + j*m];
}

static void hilb(double A[], int n) {
  // Return n by n Hilbert matrix
  int i, j;
  for (j=0; j<n; j++)
    for (i=0; i<n; i++) 
      A[i + j*n] = 1.0/(i + j + 1);
}

static int find_named_case(const char *namev[], const char *name, int Ncase) {
  for (int i = 0; i < Ncase; i++) {
    if (strcmp(namev[i], name) == 0) return i + 1;
  }
  return 0; // not found
}

static void simplerand(double A[], int n, int seed) {
// (procedure from Posix.1, 2001)
  int i, next = seed, irand;
  for (i=0; i<n; i++) {
    next = next * 1103515245 + 12345;
    irand = (next/65536) % 32768;
    A[i] = irand/32768.0;
  }
}

static bool error(FILE *errfp, const char *msg) {
  if (errfp) {
    fprintf(errfp, "Testcase error: %s\n", msg);
    fflush(errfp);
  }
  return false;
}
