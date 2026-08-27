#include <math.h>
#include <stdio.h>
#include <string.h>
#include "TestUtil.h"
#include "VarmaUtilities.h"
#include "varmapack_config.h"

static int NTOTAL = 0;
static int NFAIL = 0;
static char commonMessage[128] = "";
static char additionalMessage[128] = "";

static const char *baseName(const char *file) {
  const char *slash = strrchr(file, '/');
  const char *backslash = strrchr(file, '\\');
  if (slash && (!backslash || slash > backslash)) return slash + 1;
  if (backslash) return backslash + 1;
  return file;
}

static void recordCheck(bool ok, const char *message, const char *file, int line,
                        const char *func, const char *context) {
  if (ok) {
    xCheckOK();
    return;
  }
  xCheckFunc(message, file, line, func, context);
}

void xCheckFunc(const char *message, const char *file, int line,
                const char *func, const char *context) {
  fprintf(stderr, "%s:%d (%s): %s test failed: %s is false", baseName(file), line,
          func, commonMessage, message);
  if (additionalMessage[0]) fprintf(stderr, " (%s)", additionalMessage);
  if (context && context[0]) fprintf(stderr, " [%s]", context);
  fprintf(stderr, "\n");
  fflush(stderr);
  NTOTAL++;
  NFAIL++;
}

void xCheckOK(void) {
  NTOTAL++;
}

void xCheckInit(const char *message) {
  snprintf(commonMessage, sizeof(commonMessage), "%s", message ? message : "");
  additionalMessage[0] = 0;
  NTOTAL = 0;
  NFAIL = 0;
}

void xCheckAddMsg(const char *message) {
  snprintf(additionalMessage, sizeof(additionalMessage), "%s",
           message ? message : "");
}

int xCheckNFailures(void) {
  return NFAIL;
}

int xCheckNTotal(void) {
  return NTOTAL;
}

bool almostSame(double a, double b) {
  return relabsdiff(&a, &b, 1) < 5.0e-14;
}

bool almostEqual(double a[], double b[], int n) {
  return relabsdiff(a, b, n) < 5.0e-14;
}

void checkArrayTolFunc(double x[], double y[], int n, double tol,
                       const char *file, int line, const char *func) {
  for (int i=0; i<n; i++) {
    char context[128];
    bool ok = fabs(x[i] - y[i]) < tol;
    snprintf(context, sizeof(context), "index %d: %.17g versus %.17g, tol %.3g",
             i, x[i], y[i], tol);
    recordCheck(ok, "array values differ", file, line, func, context);
  }
}

void checkArrayZeroFunc(double x[], int n, const char *file, int line,
                        const char *func) {
  for (int i=0; i<n; i++) {
    char context[80];
    snprintf(context, sizeof(context), "index %d: %.17g", i, x[i]);
    recordCheck(x[i] == 0, "array value is not zero", file, line, func, context);
  }
}

void checkArrayFiniteFunc(double x[], int n, const char *file, int line,
                          const char *func) {
  for (int i=0; i<n; i++) {
    char context[80];
    snprintf(context, sizeof(context), "index %d: %.17g", i, x[i]);
    recordCheck(isfinite(x[i]), "array value is not finite", file, line, func,
                context);
  }
}

void checkVarmapackCleanFunc(const char *file, int line, const char *func) {
  char *message = varmapack_last_error();
  recordCheck(message == 0 || message[0] == 0, "Varmapack error is not empty", file,
              line, func, message);
}

void checkVarmapackSuccessFunc(bool ok, const char *file, int line,
                               const char *func) {
  char *message = varmapack_last_error();
  recordCheck(ok, "Varmapack call failed", file, line, func, message);
  recordCheck(message == 0 || message[0] == 0, "Varmapack error is not empty", file,
              line, func, message);
}

void checkVarmapackFailureFunc(bool ok, const char *file, int line,
                               const char *func) {
  char *message = varmapack_last_error();
  recordCheck(!ok, "Varmapack call unexpectedly succeeded", file, line, func,
              message);
  recordCheck(message != 0 && message[0] != 0, "Varmapack error is empty", file,
              line, func, 0);
}

void checkVarmapackNaNFunc(double value, const char *file, int line,
                           const char *func) {
  char *message = varmapack_last_error();
  recordCheck(isnan(value), "Varmapack result is not NaN", file, line, func, message);
  recordCheck(message != 0 && message[0] != 0, "Varmapack error is empty", file,
              line, func, 0);
}

randompack_rng *seededRngFunc(int seed, const char *file, int line,
                              const char *func) {
  randompack_rng *rng = randompack_create(0);
  recordCheck(rng != 0, "Randompack generator creation failed", file, line, func, 0);
  if (rng) reseedRngFunc(rng, seed, file, line, func);
  return rng;
}

void reseedRngFunc(randompack_rng *rng, int seed, const char *file,
                   int line, const char *func) {
  bool ok = randompack_seed(seed, 0, 0, rng);
  recordCheck(ok, "Randompack seeding failed", file, line, func,
              randompack_last_error(rng));
}

static void loadTestModel(test_model *model, const char *name, int index) {
  bool ok;
  memset(model, 0, sizeof(*model));
  model->index = index;
  if (name) snprintf(model->name, sizeof(model->name), "%s", name);
  ok = varmapack_testcase(model->name, &model->index, 0, &model->p, &model->q,
                          &model->r, 0, 0, 0, 0);
  ASSERT(ok, "testcase inquiry failed: %s", varmapack_last_error());
  ASSERT(ALLOC(model->A, model->r*model->r*(model->p > 0 ? model->p : 1)),
         "allocation failed");
  ASSERT(ALLOC(model->B, model->r*model->r*(model->q > 0 ? model->q : 1)),
         "allocation failed");
  ASSERT(ALLOC(model->Sig, model->r*model->r), "allocation failed");
  ok = varmapack_testcase(model->name, &model->index, 0, &model->p, &model->q,
                          &model->r, model->A, model->B, model->Sig, 0);
  ASSERT(ok, "testcase construction failed: %s", varmapack_last_error());
}

void loadNamedTestModel(test_model *model, const char *name) {
  loadTestModel(model, name, 0);
}

void loadIndexedTestModel(test_model *model, int index) {
  loadTestModel(model, 0, index);
}

void freeTestModel(test_model *model) {
  FREE(model->Sig);
  FREE(model->B);
  FREE(model->A);
  memset(model, 0, sizeof(*model));
}
