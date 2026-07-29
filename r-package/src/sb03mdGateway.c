#include "sb03mdGateway.h"

#include <math.h>
#include "BlasGateway.h"
#include "VarmaUtilities.h"

typedef long int sb03md_fstrlen;

void sb03md_( char *dico, char *job, char *fact,
  char *trana, int *n, double a[], int *lda, double u[], int *ldu, double c[], int *ldc,
  double *scale, double *sep, double *ferr, double wr[],
  double wi[], int iwork[], double dwork[], int *ldwork,
  int *info, sb03md_fstrlen lendico, sb03md_fstrlen lenjob, sb03md_fstrlen lenfact,
  sb03md_fstrlen lentrana);

static void call_sb03md(char dico, char job, char fact, char trana, int n,
                        double A[], int lda, double U[], int ldu, double C[],
                        int ldc, double *scale, double *sep, double *ferr,
                        double wr[], double wi[], int iwork[], double dwork[],
                        int ldwork, int *info) {
  sb03md_(&dico, &job, &fact, &trana, &n, A, &lda, U, &ldu, C, &ldc, scale,
          sep, ferr, wr, wi, iwork, dwork, &ldwork, info, 1, 1, 1, 1);
}

HIDDEN bool sb03mdGateway(double F[], double Q[], int n, double P[]) {
  double scale = 1;
  double sep = 0;
  double ferr = 0;
  double query;
  double *A = 0;
  double *U = 0;
  double *wr = 0;
  double *wi = 0;
  double *work = 0;
  int *iwork = 0;
  int lwork = -1;
  int info = 0;
  bool ok = false;
  if (!ALLOC(A, n*n)) goto done;
  if (!ALLOC(U, n*n)) goto done;
  if (!ALLOC(wr, n)) goto done;
  if (!ALLOC(wi, n)) goto done;
  if (!ALLOC(iwork, n*n > 1 ? n*n : 1)) goto done;
  lacpy("All", n, n, F, n, A, n);
  lacpy("All", n, n, Q, n, P, n);
  scal(n*n, -1, P, 1);
  call_sb03md('D', 'X', 'N', 'T', n, A, n, U, n, P, n, &scale, &sep, &ferr,
              wr, wi, iwork, &query, lwork, &info);
  if (info != 0) goto done;
  lwork = (int)query;
  if (lwork < 1) lwork = 1;
  if (!ALLOC(work, lwork)) goto done;
  lacpy("All", n, n, F, n, A, n);
  lacpy("All", n, n, Q, n, P, n);
  scal(n*n, -1, P, 1);
  call_sb03md('D', 'X', 'N', 'T', n, A, n, U, n, P, n, &scale, &sep, &ferr,
              wr, wi, iwork, work, lwork, &info);
  if (info != 0 || scale == 0 || !isfinite(scale)) goto done;
  if (scale != 1) scal(n*n, 1/scale, P, 1);
  ok = true;
done:
  FREE(work);
  FREE(iwork);
  FREE(wi);
  FREE(wr);
  FREE(U);
  FREE(A);
  return ok;
}
