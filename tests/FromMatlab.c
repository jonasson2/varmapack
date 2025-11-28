#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include "error.h"
#include "VarmaUtilities.h"
#include "FromMatlab.h"

static void SkipLine(FILE *f) {
  fscanf(f, "%*[^\n]%*c");
}

static void ReadLine(FILE *f, char *buf, int n) {
  char *p = fgets(buf, n, f);
  if (!p) buf[0] = 0;
  int len = strlen(buf);
  if (!len || buf[len-1] != '\n') SkipLine(f);
}

static bool FindVar(FILE *f, char *varname, int *m, int *n, int *line) {
  char name[128], buf[128];
  int got;
  rewind(f);
  do {
    (*line)++;
    ReadLine(f, buf, 128);
    if (buf[0] == 0) break;
    got = sscanf(buf, " %127[^,] , %d , %d", name, m, n);
    ASSERT(got >= 3, "fewer than 3 items on line %d", *line);
    if (!strcmp(name, varname)) return true;
  } while(true);
  ASSERT(false, "%s not found in compare file", varname);
}

static bool ReadVec(FILE *f, char *fmt, double x[], int n, int line) {
  for (int i = 0; i < n; i++) {
    int k = fscanf(f, fmt, &x[i]);
    ASSERT(k == 1, "illegal value when reading line %d from compare file", line);
  }
  return true;
} 

bool DoubleFromMatlab(FILE *f, char *var, double *val) {
  int m, n, line;
  bool found = FindVar(f, var, &m, &n, &line);
  if (!found) return false;
  ASSERT(m == 1 && n == 1, "expected 1,1 on line %d in compare file", line);
  if (!ReadVec(f, "%lf,", val, 1, line)) return false;
  return true;
}

bool IntFromMatlab(FILE *f, char *var, int *val) {
  int m, n, line;
  bool found = FindVar(f, var, &m, &n, &line);
  if (!found) return false;
  ASSERT(m == 1 && n == 1, "expected 1,1 on line %d in compare file", line);
  if (!ReadVec(f, "%d,", val, 1, line)) return false;
  return true;
}

bool IntVecFromMatlab(FILE *f, char *vector, int *x, int n) {
  int M, N, line;
  bool found = FindVar(f, vector, &M, &N, &line);
  if (!found) return false;
  bool ok = (M == n && N == 1) || (N == n && M == 1);
  ASSERT(ok, "expected 1,%d or %d,1 on line %d in compare file", n, n, line);
  if (!ReadVec(f, "%d,", x, n, line)) return false;
  return true;
}

bool MatrixFromMatlab(FILE *f, char *matrix, double A[], int m, int n) {
  int M, N, line;
  bool found = FindVar(f, matrix, &M, &N, &line);
  if (!found) return false;
  ASSERT(m == M && n == N, "expected %d,%d on line %d in compare file", m, n, line);
  if (!ReadVec(f, "%lf,", A, m*n, line)) return false;
  return true;
}

bool VectorFromMatlab(FILE *f, char *vector, double *x, int n) {
  return MatrixFromMatlab(f, vector, x, n, 1);
}

bool CompareWithMatlab(FILE *f, char *matrix, double *A, int m, int n, double *diff) {
  double *B;
  ALLOC(B, m*n);
  bool ok = MatrixFromMatlab(f, matrix, B, m, n);
  FREE(B);
  if (!ok) return false;
  *diff = relabsdiff(A, B, m*n);
  return true;
}
