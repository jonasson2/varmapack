#include <stdio.h>
#include <stdlib.h>
#include "error.h"
#include "printX.h"
#include "DebugUtil.h"
#include "VarmaUtilities.h"

char matmatfolder[256] = "";

bool TestMatlabMatrix(char *filename, double *A, int m, int n) {
  double *B;
  char filepath[512], buf[256], name[256];
  snprintf(filepath, sizeof(filepath), "%s/%s", matmatfolder, filename);
  FILE *fp = fopen(filepath, "r");
  ASSERT(fp, "TestMatlabMatrix error: cannot open file %s", filename);
  //
  int file_m, file_n;
  int nread = fscanf(fp, "%255[^,],%d,%d\n", name, &file_m, &file_n);
  ASSERT(nread == 3, "TestMatlabMatrix error: invalid file format in %s", filename);
  ASSERT(file_m == m && file_n == n, "%s: dimension mismatch, file (%d,%d) ≠ C (%d,%d)",
	 basename_only(filepath), file_m, file_n, m, n);
  ALLOC(B, m*n);
  for (int i = 0; i < m*n; i++) {
    nread = fscanf(fp, "%lf,", &B[i]);
    ASSERT(nread == 1, "Error: failed to read value %d from '%s'\n", i, filename);
  }
  double rd = relabsdiff(A, B, m * n);
  snprintf(buf, sizeof(buf), "Matrix %s comparison, rel.abs.diff = %.3e", name, rd);
  printMsg(buf);
  FREE(B);
  fclose(fp);
  return true;
}
