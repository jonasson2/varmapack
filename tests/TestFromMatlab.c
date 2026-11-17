#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include "TestUtil.h"
#include "Tests.h"
#include "FromMatlab.h"
#include "error.h"

static void write_matlab_file(const char *path) {
  FILE *f = fopen(path, "w");
  ASSERT(f != 0, "could not create test fixture");
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
  ASSERT(f != 0, "could not open test fixture");
  double s = 0.0;
  bool ok = DoubleFromMatlab(f, "skip", &s);
  xCheck(ok && fabs(s - 7.0) < 1e-15);
  double B[6] = {0};
  ok = MatrixFromMatlab(f, "B", B, 2, 3);
  xCheck(ok);
  checkArrayTol(B, A, 6, 1e-15);
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
