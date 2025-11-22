#include <stdbool.h>
#include "randompack.h"
#include <stdio.h>

int main(void) {
  double x;
  randompack_rng *rng = randompack_create("Park-Miller", 42);
  randompack_u01(&x, 1, rng);
  printf("u01 random number: %.4f\n", x);
}
