#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include "error.h"
#include "VarmaUtilities.h"
#include "FromMatlab.h"

static bool FindVar(FILE *f, char *varname, int *m, int *n, int *line) {
  char tag[128], name[128];
  *line = 0;
  rewind(f);
  while (fscanf(f, " %127s", tag) == 1) {
    if (tag[0] == '#') {
      fscanf(f, "%*[^\n]");
      continue;
    }
    ASSERT(!strcmp(tag, "var"), "expected var in compare file");
    int ndim;
    int count = 1;
    int dim0 = 1;
    (*line)++;
    ASSERT(fscanf(f, " %127s %d", name, &ndim) == 2,
           "illegal variable header in compare file");
    ASSERT(ndim >= 0 && ndim <= 8, "illegal rank for %s", name);
    for (int i=0; i<ndim; i++) {
      int dim;
      ASSERT(fscanf(f, " %d", &dim) == 1, "illegal dimension for %s", name);
      ASSERT(dim >= 0, "negative dimension for %s", name);
      if (i == 0) dim0 = dim;
      count *= dim;
    }
    if (!strcmp(name, varname)) {
      *m = ndim == 0 ? 1 : dim0;
      *n = ndim <= 1 ? 1 : (dim0 == 0 ? 0 : count/dim0);
      return true;
    }
    for (int i=0; i<count; i++) {
      ASSERT(fscanf(f, " %*s") == 0, "illegal value for %s", name);
    }
  }
  return false;
}

static bool ReadVec(FILE *f, char *fmt, double x[], int n, int line) {
  for (int i = 0; i < n; i++) {
    int k = fscanf(f, fmt, &x[i]);
    ASSERT(k == 1, "illegal value when reading line %d from compare file", line);
  }
  return true;
} 

static bool ReadIntVec(FILE *f, int x[], int n, int line) {
  for (int i = 0; i < n; i++) {
    int k = fscanf(f, "%d", &x[i]);
    ASSERT(k == 1, "illegal value when reading line %d from compare file", line);
  }
  return true;
}

bool DoubleFromMatlab(FILE *f, char *var, double *val) {
  int m, n, line;
  bool found = FindVar(f, var, &m, &n, &line);
  if (!found) return false;
  if (m != 1 || n != 1) return false;
  if (!ReadVec(f, "%lf", val, 1, line)) return false;
  return true;
}

bool IntFromMatlab(FILE *f, char *var, int *val) {
  int m, n, line;
  bool found = FindVar(f, var, &m, &n, &line);
  if (!found) return false;
  if (m != 1 || n != 1) return false;
  if (!ReadIntVec(f, val, 1, line)) return false;
  return true;
}

bool IntVecFromMatlab(FILE *f, char *vector, int *x, int n) {
  int M, N, line;
  bool found = FindVar(f, vector, &M, &N, &line);
  if (!found) return false;
  bool ok = (M == n && N == 1) || (N == n && M == 1);
  if (!ok) return false;
  if (!ReadIntVec(f, x, n, line)) return false;
  return true;
}

bool MatrixFromMatlab(FILE *f, char *matrix, double A[], int m, int n) {
  int M, N, line;
  bool found = FindVar(f, matrix, &M, &N, &line);
  if (!found) return false;
  if (m*n != M*N) return false;
  if (!ReadVec(f, "%lf", A, m*n, line)) return false;
  return true;
}

bool VectorFromMatlab(FILE *f, char *vector, double *x, int n) {
  return MatrixFromMatlab(f, vector, x, n, 1);
}

bool CompareWithMatlab(FILE *f, char *matrix, double *A, int m, int n, double *diff) {
  double *B = 0;
  if (m*n == 0) {
    *diff = 0;
    return MatrixFromMatlab(f, matrix, A, m, n);
  }
  ASSERT(ALLOC(B, m*n), "allocation failed");
  bool ok = MatrixFromMatlab(f, matrix, B, m, n);
  if (!ok) {
    FREE(B);
    return false;
  }
  *diff = relabsdiff(A, B, m*n);
  FREE(B);
  return true;
}
