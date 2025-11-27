#include <math.h>
#include <stdio.h>
#include "xCheck.h"
#include "Tests.h"

double ScalarFromMatlab(char *file, char *var);
double CompareWithMatlab(char *filename, char *Bname, double *A, int m, int n,
			 double tol);

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
  double s = ScalarFromMatlab((char *)fname, "skip");
  xCheck(fabs(s - 7.0) < 1e-15);
  double diff = CompareWithMatlab((char *)fname, "B", A, 2, 3, 1e-12);
  xCheck(diff < 1e-15);
  double missing = CompareWithMatlab((char *)fname, "ZZ", A, 2, 3, 1e-12);
  xCheck(isnan(missing));
  double bad_dims = CompareWithMatlab((char *)fname, "B", A, 2, 2, 1e-12);
  xCheck(isnan(bad_dims));
  remove(fname);
}
