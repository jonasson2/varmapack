#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "error.h"

#define FASSERT(f, cond, ...)				\
 do {							\
   if (!(cond)) {					\
     error_report(__FILE__, __LINE__, __VA_ARGS__);	\
     if (f) fclose(f);					\
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
  FASSERT(f, f != 0, "cannot open file");
  char name[128], buf[128];
  int m, n;
  line = 0;
  do {
    line++;
    ReadLine(f, buf, 128);
    if (!buf[0]) break;
    got = sscanf(buf, " %127[^,] , %d , %d", name, &m, &n);
    FASSERT(f, got >= 3, "fewer than 3 items on line %d in file %s", line, file);
    if (strcmp(name, var) == 0) {
      FASSERT(f, m==1 && n==1, "expected m = n = 1 on line %d in file %s", line, file);
      line++;
      ReadLine(f, buf, 128);
      int k = sscanf(buf, "%lf", &val);
      FASSERT(f, k == 1, "illegal value on line %d in file %s", line, file);
      return val;
    }
    else {
      SkipLine(f);
      line++;
    }
  } while (true);
  FASSERT(f, false, "variable %s not found in file %s", var, file);
}
