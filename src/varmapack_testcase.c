#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "BlasGateway.h"
#include "error.h"
#include "VarmaUtilities.h"
#include "randompack.h"
#include "varmapack.h"
#include "printX.h"
#define DEBUG
#include "debugprint.h"
#include "VYW.h"

static void flipud(int m, int n, double src[], double dst[]);
static void hilb(double A[], int m);
static int find_named_case(const char *namev[], const char *name, int Ncase);
static void construct_rho_A(double A[], int p, int r);
static varmapack_error scale_to_rho(double A[], int p, int r, double rho);
static double scaled_specrad(double A[], double scale, int p, int r);

varmapack_error varmapack_testcase( // Create a testcase for VARMA likelihood calculation
  double A[],    // out     r×r×p, autoregressive parameter matrices (or null)
  double B[],    // out     r×r×q, moving average parameter matrices (or null)
  double Sig[],  // out     r×r, covariance of the shock terms eps(t) (or null)
  char *name,    // out     string w. length >= 12, name of testcase, or "" for unnamed
  int *pp,       // in/out  Number of autoregressive terms (or null)
  int *qp,       // in/out  Number of moving avg. terms (or null)
  int *rp,       // in/out  dimension of each x(t) (or null)
  int *icase,    // in/out  index of named testcase to create or 0 or -1 to use p,q,r
  double rho,    // in      target spectral radius when name is "rho"
  randompack_rng *rng)  // in      random number generator
{
  // TESTCASE CREATION
  // The following kinds of testcases, suitable for testing or timing various components
  // of the Varmasim package (or the parent Varma package) may be constructed. A, B and
  // Sig return the parameter matrices of the created case and p, q and r its dimensions.
  //
  // 1. Named or indexed. When name is specified the testcase with the specified name is
  //    returned. If index is non-null, the case index is put there. When name is "" and
  //    icase > 0 the correspending named testcase, including its name, is returned. When
  //    icase = -1 or 0, name is ignored (and can be null).
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
  //    corresponding case are put in p, q and r. If name is specified icase is assigned
  //    to, and if icase is specified and name is not null, it is assigned to. varmapack_testcase(0,
  //    0, 0, "smallAR1", &p, &q, &r, &icase, 0) will thus set p, q, r, icase to 1, 0, 2,
  //    4.
  //
  // 2. Count + max dimensions. varmapack_testcase(0, 0, 0, "max", &p, &q, &r, &icase, 0) sets p, q,
  //    r to the maximum dimensions over all testcases and icase to the number of named
  //    testcases.
  //
  // SUMMARY: When icase is 0 or -1 a model with dimensions p, q, r is created; when -1
  // deterministic, and when 0 random. When 1 <= icase <= 15 and A, B and Sig are null, p,
  // q and r return dimensions of one of 15 predefined testcases, and when they are
  // non-null they return the data for these cases. See more detailes in varmapack.h.

#define p12 4
#define r12 5
#define seed12 42
#define c12 0.05
#define n12 (r12*r12*p12)
  
  const char *namev[] = {// nr   p  q  r
    "tinyAR",      // 1    1  0  1
    "tinyMA",      // 2    0  1  1
    "tinyARMA",    // 3    1  1  1
    "smallAR1",    // 4    1  0  2
    "smallAR2",    // 5    2  0  2
    "smallMA1",    // 6    0  1  2
    "smallMA2",    // 7    0  2  2
    "smallARMA1",  // 8    1  1  2
    "smallARMA2",  // 9    1  2  2
    "mediumAR",    // 10   1  0  3
    "mediumMA1",   // 11   0  1  3
    "mediumARMA1", // 12   3  3  3
    "mediumARMA2", // 13   3  3  3
    "mediumMA2",   // 14   0  2  3
    "largeAR"      // 15   4  0  5
  };
  int    pv[]      = { 1,  0,  1,  1,  2,  0,  0,  1,  1,  1,  0,  3,  3,  0, p12};
  int    qv[]      = { 0,  1,  1,  0,  0,  1,  2,  1,  2,  0,  1,  3,  3,  2, 0};
  int    rv[]      = { 1,  1,  1,  2,  2,  2,  2,  2,  2,  3,  3,  3,  3,  3, r12};
  //  Following would be used except that cl complains it is nonstandard
  //                   0   1   2   3    4   5   6   7   8   9  10   11   12   13   14
  //double *Av[]   = {A1,  0, A3, A4, A5,  0,  0, A8, A9,A10,  0, A12, A13,   0, A15};
  //double *Bv[]   = { 0, B2, B3,  0,   0, B6, B7, B8, B9,  0,B11, B12, B13, B14,   0};
  //double *Sigv[] = {S1, S1, S1, S2,  S2, S3, S3, S3, S3, S4, S4,  S4,  S4,  S4,  S6};
  int p, q, r, Ncase = sizeof(namev)/sizeof(char*);
  // SANITY CHECKS:
  if (icase == 0) return VARMAPACK_INVALID_ARGUMENT;
  int nnull = !A + !B + !Sig;
  if (nnull != 0 && nnull !=3)
    return VARMAPACK_INVALID_ARGUMENT;
  bool MAX = name && !strcmp(name, "max");
  bool RHO = name && !strcmp(name, "rho");
  bool INQUIRY = nnull == 3;
  bool NAMED = name && strlen(name) > 0 && !RHO;
  if (!INQUIRY) {
    if (NAMED) {
      *icase = find_named_case(namev, name, Ncase);
      if (*icase == 0) return VARMAPACK_UNKNOWN_TESTCASE;
    }
    if (*icase <= 0 || RHO) {
      if (!pp || !qp || !rp) return VARMAPACK_INVALID_ARGUMENT;
      if (*pp < 0 || *qp < 0 || *rp <= 0) {
        return VARMAPACK_INVALID_ARGUMENT;
      }
      if (*icase == 0 && !RHO && rng == 0) return VARMAPACK_INVALID_ARGUMENT;
    }
    if (!RHO && (*icase < -1 || *icase > Ncase)) return VARMAPACK_UNKNOWN_TESTCASE;
  }
  else if (!MAX && NAMED) {
    *icase = find_named_case(namev, name, Ncase);
    if (*icase == 0) return VARMAPACK_UNKNOWN_TESTCASE;
  }
  else if (INQUIRY && !NAMED) {
    if (*icase == 0 || *icase == -1) {
      return VARMAPACK_INVALID_ARGUMENT;
    }
    if (*icase < -1 || *icase > Ncase) {
      return VARMAPACK_UNKNOWN_TESTCASE;
    }
  }
  // MAX INQUIRY
  if (MAX) {
    *icase = Ncase;
    for (int k=0; k<Ncase; k++) {
      if (pp) *pp = imax(*pp, pv[k]);
      if (qp) *qp = imax(*qp, qv[k]);
      if (rp) *rp = imax(*rp, rv[k]);
    }
    return VARMAPACK_OK;
  }
  // CASE INQUIRY
  else if (INQUIRY) {
    if (NAMED) {
      *icase = find_named_case(namev, name, Ncase);
    }
    else if (name) {
      snprintf(name, 12, "%s", namev[*icase - 1]);
    }
    int k = *icase - 1;
    if (pp) *pp = pv[k];
    if (qp) *qp = qv[k];
    if (rp) *rp = rv[k];
    return VARMAPACK_OK;
  }
  // CONSTRUCT TESTCASE
  double
    A1[] = {0.5}, A3[] = {0.4}, A4[] = {0.1, 0.1,   0.1, 0.100000001},
    A5[] = {0.60, 0.05, 0.30, 0.02,   0.04, 0.60, 0.02, 0.30},
    A8[] = {0.3, 0.1,   0.4, 0.2}, A9[] = {0.2, 0.3,   0.2, 0.3}, A10[] = {
      0.35, 0.15, 0.15,//  0.11, 0.14, 0.07,
      0.25, 0.15, 0.05,//  0.12, 0.15, 0.08,
      0.15, 0.05, 0.01,//  0.13, 0.16, 0.09
    }, A12[9*3], A13[9*3], A33[] = {
      0.15, 0.10, 0.05,  0.11, 0.14, 0.17,  0.01, 0.04, 0.06,
      0.16, 0.11, 0.06,  0.12, 0.15, 0.18,  0.02, 0.05, 0.08,
      0.17, 0.12, 0.07,  0.13, 0.16, 0.19,  0.03, 0.06, 0.09
    }, A15[n12], B2[] = {0.5}, B3[] = {0.4},
    B6[] = {0.3, 0.1,   0.1, 0.2}, B7[] = {0.3, 0.1, 0.1, 0.1,   0.3, 0.1, 0.1, 0.1},
    B8[] = {0.2, 0.3,   0.2, 0.3}, B9[] = {0.4, 0.1, 0.2, 0.1,   0.2, 0.1, 0.3, 0.1},
    B11[] = {
      0.35, 0.25, 0.15, 0.25, 0.15, 0.05, 0.15, 0.05, 0.01
    }, B12[9*3], B13[9*3], B14[] =  {
      0.35, 0.25, 0.15,  0.11, 0.14, 0.17, 0.25, 0.15, 0.05,  0.12, 0.15, 0.18,
      0.15, 0.05, 0.01,  0.13, 0.16, 0.19
    }, S1[] = {0.8}, S2[] = {2, 1,   1, 3}, S3[] = {2, 1,   1, 2},
    S4[] = {2.0, 0.5, 0.0,   0.5, 2.0, 0.5,   0.0, 0.5, 1.0},
    //S5[] = {1.0, 0.0, 0.0,   0.0, 1.0, 0.0,   0.0, 0.0, 1.0},
    S6[r12*r12];
  copy(9*3, A33, 1, A12, 1);
  flipud(9, 3, A33, B12);
  flipud(9, 3, A33, A13);
  for (int j=0; j<3; j++)
    for (int i=3; i<9; i++) A13[i + j*9] += 0.01;
  copy(9*3, A33, 1, B13, 1);

  double *Av[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double *Bv[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double *Sigv[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  int i, j;
  bool ok;
  // Finish defining the starting matrices:
  Av[0]=A1; Av[2]=A3; Av[3]=A4; Av[4]=A5; Av[7]=A8; Av[8]=A9;
  Av[9]=A10; Av[11]=A12; Av[12]=A13; Av[14]=A15;
  Bv[1]=B2; Bv[2]=B3; Bv[5]=B6; Bv[6]=B7; Bv[7]=B8; Bv[8]=B9;
  Bv[10]=B11; Bv[11]=B12; Bv[12]=B13; Bv[13]=B14;
  Sigv[0]=Sigv[1]=Sigv[2]=S1; Sigv[3]=Sigv[4]=S2;
  Sigv[5]=Sigv[6]=Sigv[7]=Sigv[8]=S3;
  Sigv[9]=Sigv[10]=Sigv[11]=Sigv[12]=Sigv[13]=S4; Sigv[14]=S6;
  //
  randompack_rng *rng12 = randompack_create(0);
  if (rng12 == 0) return VARMAPACK_ALLOCATION;
  if (!randompack_seed(seed12, 0, 0, rng12)) {
    randompack_free(rng12);
    return VARMAPACK_INTERNAL;
  }
  if (!randompack_u01(A15, n12, rng12)) {
    randompack_free(rng12);
    return VARMAPACK_INTERNAL;
  }
  randompack_free(rng12);
  scal(n12, c12, A15, 1);
  hilb(S6, r12);
  if (RHO || *icase <= 0) {
    p = *pp;
    q = *qp;
    r = *rp;
  }
  else {
    p = pv[*icase-1];
    q = qv[*icase-1];
    r = rv[*icase-1];
  }
  if (RHO || *icase == -1 || *icase == 0) {
    if (Sig) {
      hilb(Sig, r);
      for (i = 0; i < r; i++) Sig[i + i*r] += 0.2; // add 0.2 to diagonal
    }
  }
  if (RHO) {
    if (A && p > 0) {
      construct_rho_A(A, p, r);
      varmapack_error error = scale_to_rho(A, p, r, rho);
      if (error) return error;
    }
    if (A && p == 0 && rho != 0) return VARMAPACK_INVALID_ARGUMENT;
    if (B && q>0) for (i=0; i<r*r*q; i++) B[i] = 1.0/(q*r);
  }
  else if (*icase == -1) {
    if (A && p>0) for (i=0; i<r*r*p; i++) A[i] = 0.5/(p*r);
    if (B && q>0) for (i=0; i<r*r*q; i++) B[i] = 1.0/(q*r);     
  }
  else if (*icase == 0) {
    if (A && p>0) {
      if (!randompack_u01(A, r*r*p, rng)) return VARMAPACK_INTERNAL;
      scal(r*r*p, 0.5/(p*r), A, 1); 
      double *tmpS = 0;
      if (!ALLOC(tmpS, r*r*(p+1))) return VARMAPACK_ALLOCATION;
      double *tmpC = 0;
      if (!ALLOC(tmpC, r*r*(q+1))) {
        FREE(tmpS);
        return VARMAPACK_ALLOCATION;
      }
      double *tmpG = 0;
      if (!ALLOC(tmpG, r*r*(q+1))) {
        FREE(tmpC);
        FREE(tmpS);
        return VARMAPACK_ALLOCATION;
      }
      double *tmpB = 0;
      if (q > 0) {
        if (!ALLOC(tmpB, r*r*q)) {
          FREE(tmpG);
          FREE(tmpC);
          FREE(tmpS);
          return VARMAPACK_ALLOCATION;
        }
        setzero(r*r*q, tmpB);
      }
      j = 0;
      while (true) {
        ok = VYWFactorizeSolve(A, q > 0 ? tmpB : 0, Sig, p, q, r, tmpS, tmpC, tmpG);
        if (ok) break;
        scal(r*p, 0.5, A, 1);
        j++;
        if (j >= 10) {
          if (q > 0) FREE(tmpB);
          FREE(tmpG);
          FREE(tmpC);
          FREE(tmpS);
          return VARMAPACK_INTERNAL;
        }
      }
      if (q > 0) FREE(tmpB);
      FREE(tmpG);
      FREE(tmpC);
      FREE(tmpS);
    }
    if (B && q>0) {
      if (!randompack_u01(B, r*r*q, rng)) return VARMAPACK_INTERNAL;
      scal(r*r*q, 1.0/(q*r), B, 1);
    }
  }
  else if (1 <= *icase && *icase <= Ncase) {
    if (*icase == 15) {
      for (i=0; i<r12; i++) S6[i + i*r] += 1; // add I to diagonal
      lacpy("All", r, p*r, Av[*icase - 1], r, A, r);
    }
    else {
      if (A && Av[*icase-1]) copytranspose(p*r, r, Av[*icase-1], p*r, A, r);
      if (B && Bv[*icase-1]) copytranspose(q*r, r, Bv[*icase-1], q*r, B, r);
    }
    if (Sig && Sigv[*icase-1]) copy(r*r, Sigv[*icase-1], 1, Sig, 1);
    if (name && !NAMED) strcpy(name, namev[*icase-1]); // name <---namev
  }
  else return VARMAPACK_INTERNAL;
    //     case {"pivotfailure"} % create almost singular vyw equations
    //       A = {
    //         -0.6250  0.2400    0.2500    0.1250    0.6250    0.3750
    //         0             0    0.1250    0.5000    0.2500    0.1250
    //         0.2500        0    0.1430         0    0.2500    0.2500};
    //       A(1,1) = -0.62999813363195223; %% as singular as can be made
    //       B={};
    //       Sig = eye(3);
  if (pp) *pp = p;
  if (qp) *qp = q;
  if (rp) *rp = r;
  return VARMAPACK_OK;
}

static void flipud(int m, int n, double *src, double *dst) {
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

static void construct_rho_A(double A[], int p, int r) {
  double scale = 2/(r*pow(p, 1.0/3.0));
  for (int j=0; j<p*r; j++) {
    int jj = j + 1;
    for (int i=0; i<r; i++) {
      int ii = i + 1;
      int mn = ii < jj ? ii : jj;
      int mx = ii > jj ? ii : jj;
      A[i + j*r] = scale*mn/mx;
    }
  }
}

static varmapack_error scale_to_rho(double A[], int p, int r, double rho) {
  double lo = 0, hi = 1, mid = 0, midrho = 0;
  if (rho < 0 || isnan(rho)) return VARMAPACK_INVALID_ARGUMENT;
  if (p == 0) {
    if (rho != 0) return VARMAPACK_INVALID_ARGUMENT;
    return VARMAPACK_OK;
  }
  if (rho == 0) {
    setzero(r*r*p, A);
    return VARMAPACK_OK;
  }
  for (int iter=0; scaled_specrad(A, hi, p, r) <= rho; iter++) {
    hi = 2*hi;
    if (iter >= 60) return VARMAPACK_INTERNAL;
  }
  for (int iter=0; iter<100; iter++) {
    mid = (lo + hi)/2;
    midrho = scaled_specrad(A, mid, p, r);
    if (isnan(midrho)) return VARMAPACK_INTERNAL;
    if (fabs(midrho - rho) < 1e-4) {
      scal(r*r*p, mid, A, 1);
      return VARMAPACK_OK;
    }
    if (midrho < rho) lo = mid;
    else hi = mid;
  }
  scal(r*r*p, (lo + hi)/2, A, 1);
  return VARMAPACK_OK;
}

static double scaled_specrad(double A[], double scale, int p, int r) {
  double *As = 0, rho;
  if (!ALLOC(As, r*r*p)) return NAN;
  copy(r*r*p, A, 1, As, 1);
  scal(r*r*p, scale, As, 1);
  rho = varmapack_specrad(As, r, p);
  FREE(As);
  return rho;
}
