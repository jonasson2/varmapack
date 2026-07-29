#include <stdio.h>
#include "varmapack.h"

int main(void) {
  int status = 0;
  int r = 2, n = 200, M = 10;
  double A[] = {0.6, 0, 0.1, 0.4};
  double Sig[] = {2, 0, 0, 1};
  double X[2*200*10];
  randompack_rng *rng = randompack_create(0);
  if (!rng) {
    fprintf(stderr, "could not create RNG\n");
    return 1;
  }
  if (!randompack_seed(123, 0, 0, rng)) {
    fprintf(stderr, "Randompack: %s\n", randompack_last_error(rng));
    status = 1;
  }
  else if (!varmapack_sim(A, 0, Sig, 0, 0, 1, 0, r, n, M, 0, 0, 1, X, 0,
                           rng)) {
    fprintf(stderr, "Varmapack: %s\n", varmapack_last_error());
    status = 1;
  }
  randompack_free(rng);
  return status;
}
