#include "error.h"
#include "varmapack.h"

bool varmapack_testcasex(int s, int d, int r, int n, double C[], double z[]) {
  double pattern[] = {3, -1, 4, 0, -5, 2, 1, -3, 5,
                      -2, 0, 4, -4, 1, -1, 2, -3};
  clear_error();
  if (s < 0 || d < 0 || r <= 0 || n < 0 || (s == 0 && d != 0) ||
      (s > 0 && d == 0)) return fail_error("invalid dimension(s)");
  if (s > 0 && (r > INT_MAX/d || r*d > INT_MAX/s || (n > 0 && d > INT_MAX/n))) {
    return fail_error("invalid dimension(s)");
  }
  if ((s > 0 && C == 0) || (d > 0 && n > 0 && z == 0)) {
    return fail_error("invalid argument");
  }
  for (int k=0; k<s; k++) {
    for (int j=0; j<d; j++) {
      for (int i=0; i<r; i++) {
        double value = (3*(i + 1) + 5*(k + 1) + 7*j) % 11 - 5;
        C[i + j*r + k*r*d] = value/10;
      }
    }
  }
  for (int j=0; j<n; j++) {
    for (int i=0; i<d; i++) z[i + j*d] = pattern[(i + j) % LEN(pattern)]/5;
  }
  return true;
}
