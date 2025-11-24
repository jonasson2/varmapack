#include <stdbool.h>
#include "randompack.h"
#include <stdio.h>
#include "printX.h"
#include "BlasGateway.h"

int main(void) {
  // Test dsyev_ interface with a simple 2x2 symmetric matrix
  // Matrix: [2 1]
  //         [1 2]
  // Expected eigenvalues: 1 and 3
  int n = 2;
  double a[4] = {2.0, 1.0, 1.0, 2.0};
  double w[2];
  double work[20];
  int lwork = 20;
  int info;

  printf("Testing syev interface\n");
  printf("Matrix (column-major):\n");
  printf("  %6.3f %6.3f\n", a[0], a[2]);
  printf("  %6.3f %6.3f\n", a[1], a[3]);

  syev("N", "U", n, a, n, w, work, lwork, &info);

  if (info == 0) {
    printf("Eigenvalues: %6.3f %6.3f\n", w[0], w[1]);
  }
  else {
    printf("Error: info = %d\n", info);
  }

  return 0;
}
