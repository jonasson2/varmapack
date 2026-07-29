#include "error.h"
#include "varmapack.h"

bool varmapack_testcasex(int s, int r, int n, double C[], double z[]) {
  double pattern[] = {3, -1, 4, 0, -5, 2, 1, -3, 5,
                      -2, 0, 4, -4, 1, -1, 2, -3};
  clear_error();
  if (s < 0 || r <= 0 || n < 0) return fail_error("invalid dimension(s)");
  if (s > 0 && r > INT_MAX/s) return fail_error("invalid dimension(s)");
  if ((s > 0 && C == 0) || (n > 0 && z == 0)) return fail_error("invalid argument");
  for (int j=0; j<s; j++) {
    for (int i=0; i<r; i++) {
      double value = (3*(i + 1) + 5*(j + 1)) % 11 - 5;
      C[i + j*r] = value/10;
    }
  }
  for (int i=0; i<n; i++) z[i] = pattern[i % LEN(pattern)]/5;
  return true;
}
