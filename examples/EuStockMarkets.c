#include <math.h>
#include <stdio.h>
#include <varmapack.h>

int main(void) {
  int p = 1, q = 0, r = 4, n = 260, M = 10;
  double A[] = {
    0.0046240972396264, -0.00930516507225788,
    -0.0265233470280952, -0.0102959032743018,
    -0.0957618300123132, -0.00717196869314693,
    -0.11365836007501, -0.0892451180978641,
    0.0399411319128781, 0.0378105515568093,
    0.0637551035045599, -0.0031969313874636,
    0.0485658203952483, 0.0682577473923973,
    0.0915506339486754, 0.164089912507651
  };
  double Sig[] = {
    0.000106041685893524, 0.0000670486580566326,
    0.0000831286412891527, 0.0000520859447746431,
    0.0000670486580566325, 0.0000853006845924284,
    0.0000627155169059886, 0.0000427156427225213,
    0.0000831286412891527, 0.0000627155169059886,
    0.000121206588007163, 0.0000563191801322241,
    0.0000520859447746431, 0.0000427156427225213,
    0.0000563191801322241, 0.0000624949277411293
  };
  double mu[] = {
    0.000652041747691327, 0.000817899655305225,
    0.000437053986900166, 0.000431985076649575
  };
  double last[] = {5473.72, 7676.3, 3995, 5455};
  double X[4*260*10];
  char *names[] = {"DAX", "SMI", "CAC", "FTSE"};
  int status = 0;
  randompack_rng *rng = randompack_create(0);
  if (!rng) {
    fprintf(stderr, "could not create RNG\n");
    return 1;
  }
  if (!randompack_seed(2026, 0, 0, rng)) {
    fprintf(stderr, "Randompack: %s\n", randompack_last_error(rng));
    status = 1;
  }
  else if (!varmapack_sim(A, 0, Sig, mu, 1, p, q, r, n, M, 0, 0, 1,
                           X, 0, rng)) {
    fprintf(stderr, "Varmapack: %s\n", varmapack_last_error());
    status = 1;
  }
  else {
    printf("AR spectral radius: %.6f\n", varmapack_specrad(A, r, p));
    printf("First simulated path, final synthetic index levels:\n");
    for (int i=0; i<r; i++) {
      double totalReturn = 0;
      for (int t=0; t<n; t++) totalReturn += X[i + r*t];
      printf("  %-4s %9.2f\n", names[i], last[i]*exp(totalReturn));
    }
  }
  randompack_free(rng);
  return status;
}
