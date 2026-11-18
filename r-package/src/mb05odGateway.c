#include "mb05odGateway.h"

#include <math.h>
#include "BlasGateway.h"
#include "error.h"

#if defined(USING_R)
typedef fstrlen mb05od_fstrlen;
#else
typedef long int mb05od_fstrlen;
#endif

#ifdef VARMAPACK_FORTRAN_UPPERCASE
#define mb05od_ MB05OD
#endif

void mb05od_(char *balanc, int *n, int *ndiag, double *delta, double A[],
              int *lda, int *mdig, int *idig, int iwork[], double dwork[],
              int *ldwork, int *iwarn, int *info, mb05od_fstrlen lenbalanc);

HIDDEN bool mb05odGateway(double A[], double delta, int n, double E[],
                          mb05od_accuracy *accuracy) {
  return mb05odExponential(A, delta, n, E, true, accuracy);
}

HIDDEN bool mb05odExponential(double A[], double delta, int n, double E[],
                              bool balance, mb05od_accuracy *accuracy) {
  double *result = 0, *work = 0;
  int *iwork = 0;
  int nn, ndiag = 9, info = 0, lda = n > 0 ? n : 1;
  char balanc = balance ? 'S' : 'N', message[128];
  mb05od_accuracy status = {0};
  bool ok = false;
  clear_error();
  if (accuracy) *accuracy = status;
  if (!isfinite(delta) || !intProduct(n, n, &nn) || (n > 0 && (!A || !E)))
    return fail_error("invalid matrix exponential dimensions or input");
  size_t count = n > 1 ? 2*(size_t)nn + 10*(size_t)n + 9 : 1;
  if (count > INT_MAX) return fail_error("matrix exponential workspace too large");
  for (int i=0; i<nn; i++) {
    if (!isfinite(A[i]) || !isfinite(A[i]*delta))
      return fail_error("matrix exponential input and delta*A must be finite");
  }
  int lwork = count;
  if (!ALLOC(result, nn > 0 ? nn : 1) || !ALLOC(iwork, n > 0 ? n : 1) ||
      !ALLOC(work, lwork)) {
    fail_error("allocation failed");
    goto done;
  }
  if (n > 0) copy(nn, A, 1, result, 1);
  mb05od_(&balanc, &n, &ndiag, &delta, result, &lda, &status.mdig,
           &status.idig, iwork, work, &lwork, &status.iwarn, &info, 1);
  if (accuracy) *accuracy = status;
  if (info != 0) {
    STRSETF(message, "SLICOT MB05OD failed (INFO=%d)", info);
    fail_error(message);
    goto done;
  }
  for (int i=0; i<nn; i++) {
    if (!isfinite(result[i])) {
      fail_error("matrix exponential output is not finite");
      goto done;
    }
  }
  if (n > 0) copy(nn, result, 1, E, 1);
  ok = true;
done:
  FREE(work);
  FREE(iwork);
  FREE(result);
  return ok;
}
