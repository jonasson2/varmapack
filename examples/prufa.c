#include <stdio.h>
#include "BlasGateway.h"
#include "printX.h"
#include "RandomNumbers.h"

int main(void) {
  double A[6] = {1, 2, 3, 4, 5, 6};
  printM("A", A, 3, 2);
  printM("vec(A)", A, 6, 1);
  printMT("A'", A, 3, 2);
  printMT("vec(A)'", A, 6, 1);
  return 0;
}
