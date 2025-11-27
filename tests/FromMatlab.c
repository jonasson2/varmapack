#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include "error.h"
#include "VarmaUtilities.h"

#define FASSERT(f, B, cond, ...)				\
 do {							\
   if (!(cond)) {					\
     error_report(__FILE__, __LINE__, __VA_ARGS__);	\
     if (f) fclose(f);					\
     if (B) FREE(B);					\
     varmapack_set_errprefix("");			\
     return NAN;					\
   }							\
 } while (0)

static void SkipLine(FILE *f) { fscanf(f, "%*[^\n]%*c"); }

static void ReadLine(FILE *f, char *buf, int n) {
  char *p = fgets(buf, n, f);
  if (!p) buf[0] = 0;
  int len = strlen(buf);
  if (!len || buf[len-1] != '\n') SkipLine(f);
}

double ScalarFromMatlab(char *file, char *var) {
  double val = NAN;
  int line, got;
  varmapack_set_errprefix("ScalarFromMatlab");
  FILE *f = fopen(file, "r");
  FASSERT(f, 0, f != 0, "cannot open file");
  char name[128], buf[128];
  int m, n;
  line = 0;
  do {
    line++;
    ReadLine(f, buf, 128);
    if (!buf[0]) break;
    got = sscanf(buf, " %127[^,] , %d , %d", name, &m, &n);
    FASSERT(f, 0, got >= 3, "fewer than 3 items on line %d in file %s",
	    line, file);
    if (strcmp(name, var) == 0) {
      FASSERT(f, 0, m==1 && n==1, "expected m = n = 1 on line %d in file %s",
	      line, file);
      line++;
      ReadLine(f, buf, 128);
      int k = sscanf(buf, "%lf", &val);
      FASSERT(f, 0, k == 1, "illegal value on line %d in file %s", line, file);
      return val;
    }
    else {
      SkipLine(f);
      line++;
    }
  } while (true);
  FASSERT(f, 0, false, "variable %s not found in file %s", var, file);
}

double CompareWithMatlab(char *filename, char *Bname, double *A, int m, int n,
			 double tol) {
  char name[128], buf[128];
  int file_m, file_n, line = 0, got;
  double *B = 0;
  varmapack_set_errprefix("CompareWithMatlab");
  FILE *f = fopen(filename, "r");
  FASSERT(f, 0, f != 0, "cannot open file %s", filename);
  do {
    line++;
    ReadLine(f, buf, 128);
    if (!buf[0]) break;
    got = sscanf(buf, " %127[^,] , %d , %d", name, &file_m, &file_n);
    FASSERT(f, 0, got >= 3, "fewer than 3 items on line %d in file %s",
	    line, filename);
    if (strcmp(name, Bname) == 0) {
      FASSERT(f, 0, file_m == m && file_n == n,
	      "expected (%d,%d) for %s on line %d in file %s",
	      m, n, Bname, line, filename);
      B = calloc((size_t)(m*n), sizeof(*B));
      FASSERT(f, B, B || m*n == 0, "allocation failed for B (len=%d)", m*n);
      for (int i=0; i<m*n; i++) {
        int k = fscanf(f, "%lf,", &B[i]);
        FASSERT(f, B, k == 1, "illegal value %d on line %d in file %s",
		i + 1, line + 1, filename);
      }
      fclose(f);
      double diff = relabsdiff(A, B, m*n);
      FREE(B);
      varmapack_set_errprefix("");
      (void)tol;
      return diff;
    }
    SkipLine(f);
    line++;
  } while (true);
  FASSERT(f, B, false, "matrix %s not found in file %s", Bname, filename);
}
