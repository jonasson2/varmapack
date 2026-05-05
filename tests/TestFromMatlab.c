#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include "xCheck.h"
#include "Tests.h"
#include "error.h"
#include "FromMatlab.h"

static void write_matlab_file(const char *path) {
  FILE *f = fopen(path, "w");
  xCheck(f != 0);
  fprintf(f, "# FromMatlab test fixture\n");
  fprintf(f, "var skip 0\n7\n");
  fprintf(f, "var B 2 2 3\n1 2 3 4 5 6\n");
  fclose(f);
}

void TestFromMatlab(void) {
  const char fname[] = "frommatlab_test.tmp";
  double A[] = {1, 2, 3, 4, 5, 6};
  write_matlab_file(fname);
  FILE *f = fopen(fname, "r");
  xCheck(f != 0);
  double s = 0.0;
  bool ok = DoubleFromMatlab(f, "skip", &s);
  xCheck(ok && fabs(s - 7.0) < 1e-15);
  double B[6] = {0};
  ok = MatrixFromMatlab(f, "B", B, 2, 3);
  xCheck(ok);
  for (int i=0; i<6; i++) xCheck(fabs(B[i] - A[i]) < 1e-15);
  double diff = 0.0;
  ok = CompareWithMatlab(f, "B", A, 2, 3, &diff);
  xCheck(ok && diff < 1e-15);
  varmapack_set_errstream(0);
  ok = CompareWithMatlab(f, "ZZ", A, 2, 3, &diff);
  xCheck(!ok);
  ok = CompareWithMatlab(f, "B", A, 2, 2, &diff);
  xCheck(!ok);
  varmapack_set_errstream(stderr);
  fclose(f);
  remove(fname);
}
