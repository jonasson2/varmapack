// Make test data for VARMA simulation (or likelihood calculation)
#ifndef TESTCASE_H
#define TESTCASE_H
#ifdef __cplusplus
extern "C" {
#endif

#include "RandomNumbers.h"

bool Testcase( // Create a testcase for VARMA likelihood calculation
  double A[],    // out     r×r×p, autoregressive parameter matrices (or null)
  double B[],    // out     r×r×q, moving average parameter matrices (or null)
  double Sig[],  // out     r×r, covariance of the shock terms eps(t) (or null)
  char *name,    // out     string w. length >= 12, name of testcase, or ""
  int *pp,       // in/out  Number of autoregressive terms
  int *qp,       // in/out  Number of moving avg. terms
  int *rp,       // in/out  dimension of each x(t)
  int *icase,    // in/out  index of named testcase to create or 0 or -1 to use p,q,r
  RandRng *rng,  // in      random number generator
  FILE *fp);     // in      stream to print errors (or null)

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
  //    to, and if icase is specified and name is not null, it is assigned to. Testcase(0,
  //    0, 0, "smallAR", &p, &q, &r, &icase, 0) will thus set p, q, r, icase to 1, 0, 2,
  //    4.
  //
  // 2. Count + max dimensions. Testcase(0, 0, 0, "max", &p, &q, &r, &icase, 0) sets p, q,
  //    r to the maximum dimensions over all testcases and icase to the number of named
  //    testcases.

#ifdef __cplusplus
}
#endif
#endif
