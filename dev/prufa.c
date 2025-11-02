#include "printX.h"
#include "RandomNumbers.h"
int main() {
  double x[2];
  RandRng *rng = RandCreate();
  printMsg("a");
  RandSetPM(rng);
  printMsg("b");
  RandSeed(42, rng);
  printMsg("c");
  Rand(x, 2, rng);
  printV("x", x, 2);
}
