#include <stdio.h>
#include <stdlib.h>
#include <varmapack.h>

int main(void) {
  int p = 1, q = 0, r = 2, n = 200, M = 10;
  double Avar[] = {0.6, 0, 0.1, 0.4}; // r by r by p
  double Sigvar[] = {2, 0, 0, 1};     // r by r
  double X[2*200*10];                 // r by n by M
  randompack_rng *rng = randompack_create(0); // select default Randompack engine
  randompack_seed(123, 0, 0, rng); // optional, but gives reproducible results
  varmapack_sim(Avar, 0, Sigvar, 0, 0, p, q, r, n, M, 0, 0, 1, X, 0, rng);
  for (int t=0; t<5; t++) printf("%8.4f %8.4f\n", X[r*t], X[1 + r*t]);
  char name[VARMAPACK_TESTCASE_NAME_LEN] = "smallARMA1";
  int index = 0;
  varmapack_testcase(name, &index, 0, &p, &q, &r, 0, 0, 0, 0);
  double *A = malloc((p ? p : 1)*r*r*sizeof(*A));
  double *B = malloc((q ? q : 1)*r*r*sizeof(*B));
  double *Sig = malloc(r*r*sizeof(*Sig));
  varmapack_testcase(name, &index, 0, &p, &q, &r, A, B, Sig, 0);
  double mu[] = {0, 0, 1, 1};        // r by nmu (time dependent means)
  double X0[] = {0, 0, 0.1, -0.1};   // r by nX0 (common start for all replicates)
  double E[2*200*10];                // r by n by M (returns shocks)
  varmapack_sim(A, B, Sig, mu, 2, p, q, r, n, M, X0, 2, 1, X, E, rng);
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
  for (int k=0; k<=maxlag; k++) printf("Gamma[%d]: %.4f %.4f %.4f %.4f\n",
    k, Gamma[4*k], Gamma[4*k+1], Gamma[4*k+2], Gamma[4*k+3]);
  for (int k=0; k<=maxlag; k++) printf("GammaHat[%d]: %.4f %.4f %.4f %.4f\n",
    k, GammaHat[4*k], GammaHat[4*k+1], GammaHat[4*k+2], GammaHat[4*k+3]);
  for (int k=0; k<=maxlag; k++) printf("Psi[%d]: %.4f %.4f %.4f %.4f\n",
    k, Psi[4*k], Psi[4*k+1], Psi[4*k+2], Psi[4*k+3]);
  for (int k=0; k<=maxlag; k++) printf("Theta[%d]: %.4f %.4f %.4f %.4f\n",
    k, Theta[4*k], Theta[4*k+1], Theta[4*k+2], Theta[4*k+3]);
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
