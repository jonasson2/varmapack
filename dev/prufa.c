#include "printX.h"
int main() {
#include "matlabcompare.txt"
  for (int k=0; k<ncase; k++) {
    print3I("p, q, r", p[k], q[k], r[k]);
    printM("A", A[k], r[k], r[k]*p[k]);
    printM("B", B[k], r[k], r[k]*q[k]);
    printM("Sig", Sig[k], r[k], r[k]);
    printM("X", X[k], r[k], n[k]);
    
  }
}
