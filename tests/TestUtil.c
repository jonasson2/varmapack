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

void groupCheckFunc(test_check_group *group, bool ok, const char *condition,
                    const char *file, int line, const char *func, const char *context) {
  group->total++;
  if (ok) return;
  group->failed++;
  if (group->failed != 1) return;
  group->first_failed = group->total;
  group->file = file;
  group->line = line;
  group->func = func;
  STRSET(group->condition, condition);
  STRSET(group->context, context);
}

void checkGroupFunc(test_check_group *group, const char *description,
                    const char *file, int line, const char *func) {
  char context[448];
  if (group->failed == 0) {
    recordCheck(true, description, file, line, func, 0);
    return;
  }
  snprintf(context, sizeof(context), "%d of %d checks failed; first at check %d: %s%s%s",
           group->failed, group->total, group->first_failed, group->condition,
           group->context[0] ? "; " : "", group->context);
  recordCheck(false, description, group->file, group->line, group->func, context);
}

bool almostSame(double a, double b) {
  return relabsdiff(&a, &b, 1) < 5.0e-14;
}

bool almostEqual(double a[], double b[], int n) {
  return relabsdiff(a, b, n) < 5.0e-14;
}

void checkArrayTolFunc(double x[], double y[], int n, double tol,
                       const char *file, int line, const char *func) {
  int first = -1, mismatches = 0;
  char context[160];
  for (int i=0; i<n; i++) {
    if (!(fabs(x[i] - y[i]) < tol)) {
      if (first < 0) first = i;
      mismatches++;
    }
  }
  if (mismatches > 0) {
    snprintf(context, sizeof(context),
             "%d of %d entries differ; first index %d: %.17g versus %.17g, tol %.3g",
             mismatches, n, first, x[first], y[first], tol);
  }
  recordCheck(mismatches == 0, "array values differ", file, line, func,
              mismatches ? context : 0);
}

void checkArrayZeroFunc(double x[], int n, const char *file, int line,
                        const char *func) {
  int first = -1, nonzero = 0;
  char context[128];
  for (int i=0; i<n; i++) {
    if (x[i] != 0) {
      if (first < 0) first = i;
      nonzero++;
    }
  }
  if (nonzero > 0) {
    snprintf(context, sizeof(context),
             "%d of %d entries are not zero; first index %d: %.17g",
             nonzero, n, first, x[first]);
  }
  recordCheck(nonzero == 0, "array is not zero", file, line, func,
              nonzero ? context : 0);
}

void checkArrayFiniteFunc(double x[], int n, const char *file, int line,
                          const char *func) {
  int first = -1, nonfinite = 0;
  char context[128];
  for (int i=0; i<n; i++) {
    if (!isfinite(x[i])) {
      if (first < 0) first = i;
      nonfinite++;
    }
  }
  if (nonfinite > 0) {
    snprintf(context, sizeof(context),
             "%d of %d entries are not finite; first index %d: %.17g",
             nonfinite, n, first, x[first]);
  }
  recordCheck(nonfinite == 0, "array is not finite", file, line, func,
              nonfinite ? context : 0);
}

void checkVarmapackCleanFunc(const char *file, int line, const char *func) {
  char *message = varmapack_last_error();
  recordCheck(message == 0 || message[0] == 0, "Varmapack error is not empty", file,
              line, func, message);
}

void checkVarmapackSuccessFunc(bool ok, const char *file, int line,
                               const char *func) {
  char *message = varmapack_last_error();
  recordCheck(ok && (message == 0 || message[0] == 0),
              ok ? "Varmapack error is not empty" : "Varmapack call failed",
              file, line, func, message);
}

void checkVarmapackFailureFunc(bool ok, const char *file, int line,
                               const char *func) {
  char *message = varmapack_last_error();
  recordCheck(!ok && message != 0 && message[0] != 0,
              ok ? "Varmapack call unexpectedly succeeded" : "Varmapack error is empty",
              file, line, func, message);
}

void checkVarmapackNaNFunc(double value, const char *file, int line,
                           const char *func) {
  char *message = varmapack_last_error();
  recordCheck(isnan(value) && message != 0 && message[0] != 0,
              isnan(value) ? "Varmapack error is empty" : "Varmapack result is not NaN",
              file, line, func, message);
}

randompack_rng *seededRngFunc(int seed, const char *file, int line,
                              const char *func) {
  randompack_rng *rng = randompack_create(0);
  if (!rng) {
    xCheckFunc("Randompack generator creation failed", file, line, func, 0);
    exit(EXIT_FAILURE);
  }
  reseedRngFunc(rng, seed, file, line, func);
  return rng;
}

void reseedRngFunc(randompack_rng *rng, int seed, const char *file,
                   int line, const char *func) {
  bool ok = randompack_seed(seed, 0, 0, rng);
  if (!ok) {
    xCheckFunc("Randompack seeding failed", file, line, func,
                randompack_last_error(rng));
    exit(EXIT_FAILURE);
  }
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
