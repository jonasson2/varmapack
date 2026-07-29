#include <stdio.h>
#include <stdlib.h>
#include "varmapack.h"

int main(void) {
  int p = 1, q = 0, r = 2, n = 200, M = 10;
  double Avar[] = {0.6, 0, 0.1, 0.4};
  double Sigvar[] = {2, 0, 0, 1};
  double X[2*200*10];
  randompack_rng *rng = randompack_create(0);
  randompack_seed(123, 0, 0, rng);
  // Simulate a VAR(1) model.
  varmapack_sim(Avar, 0, Sigvar, 0, 0, p, q, r, n, M, 0, 0, 1, X, 0, rng);
  for (int t=0; t<5; t++) printf("%8.4f %8.4f\n", X[r*t], X[1 + r*t]);
  // Build the smallARMA1 testcase after first inquiring about its dimensions.
  char name[VARMAPACK_TESTCASE_NAME_LEN] = "smallARMA1";
  int index = 0;
  varmapack_testcase(name, &index, 0, &p, &q, &r, 0, 0, 0, 0);
  double *A = malloc((p ? p : 1)*r*r*sizeof(*A));
  double *B = malloc((q ? q : 1)*r*r*sizeof(*B));
  double *Sig = malloc(r*r*sizeof(*Sig));
  varmapack_testcase(name, &index, 0, &p, &q, &r, A, B, Sig, 0);
  // Simulate the testcase with a time-dependent mean and supplied start.
  double mu[] = {0, 0, 1, 1};
  double X0[] = {0, 0, 0.1, -0.1};
  double E[2*200*10];
  varmapack_sim(A, B, Sig, mu, 2, p, q, r, n, M, X0, 2, 1, X, E, rng);
  // Simulate and analyze a zero-mean realization of the testcase.
  varmapack_sim(A, B, Sig, 0, 0, p, q, r, n, M, 0, 0, 1, X, 0, rng);
  enum { maxlag = 3 };
  double Gamma[2*2*(maxlag+1)], Psi[2*2*(maxlag+1)];
  double Theta[2*2*(maxlag+1)], GammaHat[2*2*(maxlag+1)];
  double rho = varmapack_specrad(A, r, p);
  double rhoMA = varmapack_ma_specrad(B, r, q);
  varmapack_acvf(A, B, Sig, p, q, r, Gamma, maxlag);
  varmapack_autocov("N", "ML", r, n, X, maxlag, GammaHat);
  varmapack_psi(A, B, p, q, r, maxlag, Psi);
  varmapack_irf(A, B, Sig, p, q, r, maxlag, Theta);
  printf("AR spectral radius: %.4f\n", rho);
  printf("MA spectral radius: %.4f\n", rhoMA);
  // Add deterministic exogenous testcase data and simulate a VARMAX model.
  int s = 2, d = 2;
  double C[2*2*2];
  double z[2*200];
  varmapack_testcasex(s, d, r, n, C, z);
  varmapack_simx(A, B, C, Sig, z, 1, p, q, s, d, r, n, M, X0, 2, 1, X, E, rng);
  randompack_free(rng);
  free(Sig);
  free(B);
  free(A);
  return 0;
}
