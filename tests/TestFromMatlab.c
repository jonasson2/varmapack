#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include "xCheck.h"
#include "Tests.h"
#include "FromMatlab.h"

static void write_matlab_file(const char *path) {
  FILE *f = fopen(path, "w");
  xCheck(f != 0);
  fprintf(f, "skip,1,1\n7,\n");
  fprintf(f, "B,2,3\n1,2,3,4,5,6,\n");
  fclose(f);
}

void TestFromMatlab(void) {
  const char fname[] = "frommatlab_test.txt";
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
  ok = CompareWithMatlab(f, "ZZ", A, 2, 3, &diff);
  xCheck(!ok);
  ok = CompareWithMatlab(f, "B", A, 2, 2, &diff);
  xCheck(!ok);
  fclose(f);
  remove(fname);
}
