// Shared assertions, comparisons, RNG setup, and testcase fixtures for C tests.

#ifndef TESTUTIL_H
#define TESTUTIL_H

#include <stdbool.h>
#include "randompack.h"
#include "varmapack.h"

#define xCheck(e) ((e) ? (void)xCheckOK() : \
  xCheckFunc(#e, __FILE__, __LINE__, __func__, 0))

#define xCheckMessage(e, msg) ((e) ? (void)xCheckOK() : \
  xCheckFunc(#e, __FILE__, __LINE__, __func__, (msg)))

#define checkArrayTol(x, y, n, tol) \
  checkArrayTolFunc((x), (y), (n), (tol), __FILE__, __LINE__, __func__)

#define checkArraySame(x, y, n) \
  checkArrayTolFunc((x), (y), (n), 1e-15, __FILE__, __LINE__, __func__)

#define checkArrayZero(x, n) \
  checkArrayZeroFunc((x), (n), __FILE__, __LINE__, __func__)

#define checkArrayFinite(x, n) \
  checkArrayFiniteFunc((x), (n), __FILE__, __LINE__, __func__)

#define checkVarmapackClean() \
  checkVarmapackCleanFunc(__FILE__, __LINE__, __func__)

#define checkVarmapackSuccess(ok) \
  checkVarmapackSuccessFunc((ok), __FILE__, __LINE__, __func__)

#define checkVarmapackFailure(ok) \
  checkVarmapackFailureFunc((ok), __FILE__, __LINE__, __func__)

#define checkVarmapackNaN(value) \
  checkVarmapackNaNFunc((value), __FILE__, __LINE__, __func__)

#define seededRng(seed) seededRngFunc((seed), __FILE__, __LINE__, __func__)

#define reseedRng(rng, seed) \
  reseedRngFunc((rng), (seed), __FILE__, __LINE__, __func__)

typedef struct {
  char name[VARMAPACK_TESTCASE_NAME_LEN];
  int index;
  int p;
  int q;
  int r;
  double *A;
  double *B;
  double *Sig;
} test_model;

void xCheckFunc(const char *message, const char *file, int line,
                const char *func, const char *context);
void xCheckOK(void);
void xCheckInit(const char *message);
void xCheckAddMsg(const char *message);
int xCheckNFailures(void);
int xCheckNTotal(void);

bool almostSame(double a, double b);
bool almostEqual(double a[], double b[], int n);
void checkArrayTolFunc(double x[], double y[], int n, double tol,
                       const char *file, int line, const char *func);
void checkArrayZeroFunc(double x[], int n, const char *file, int line,
                        const char *func);
void checkArrayFiniteFunc(double x[], int n, const char *file, int line,
                          const char *func);
void checkVarmapackCleanFunc(const char *file, int line, const char *func);
void checkVarmapackSuccessFunc(bool ok, const char *file, int line,
                               const char *func);
void checkVarmapackFailureFunc(bool ok, const char *file, int line,
                               const char *func);
void checkVarmapackNaNFunc(double value, const char *file, int line,
                           const char *func);
randompack_rng *seededRngFunc(int seed, const char *file, int line,
                              const char *func);
void reseedRngFunc(randompack_rng *rng, int seed, const char *file,
                   int line, const char *func);
void loadNamedTestModel(test_model *model, const char *name);
void loadIndexedTestModel(test_model *model, int index);
void freeTestModel(test_model *model);

#endif
