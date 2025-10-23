// mds.h: Declarations related to mds (matrix-derivative-sparsity) structure
#ifndef MDS_H
#define MDS_H
#include "VarmaUtilities.h"
#include "xAssert.h"

struct mds { // A structure of a matrix, its derivatives, and sparsity codes
  double *mat;  // An r×ncol matrix
  double *der;  // r×r×(r·r·nparmat) array with derivatives of mat
  char *spcod;  // sparsity codes (see below)
  int r;        // dimension a time series observations (mat must have r rows)
  int ncol;     // column count (mat need not have r columns)
  int nparmat;  // total number of parameter matrices, p + q + 1
};

// NOTE: spcod[k]  Meaning that for each (l,c) the following holds:
//         'z'     der(,,l,c,k) is the r×ncol zero matrix
//         'e'     der(l,c,l,c,k) is the only nonzero in der(,,l,c,k)
//         'r'     The nonzeros of der(,,l,c,k) are all in its l-th row
//         'c'     The nonzeros of der(,,l,c,k) are all in its j-th column
//         'x'     The nonzeros of der(,,l,c,k) are all in its l-th column
//         'f'     der(,,l,c,k) is a full matrix for each (l,c)
//
// For sparsitycodes 'e','r','c' and 'x' ncol must be equal to r

void MdsAddProdXY(char *trnsp, struct mds F, double X[], double Y[], int nY, 
                   int kX, int kY);
void MdsAddProdXG(char *trnsp, struct mds F, double X[], struct mds G, int kX);
void MdsSetZero(struct mds *F, int r, int ncol, int nparmat);
void MdsSetParmat(struct mds *F, double X[], int kX, int r, int nparmat);
void MdsDestroy(struct mds F);
void MdsCopy(struct mds F, struct mds *G);
void der2array ( // Convert matrices-derivatives in MDS format to double arrays
  struct mds FS[], // in    n-vector of MDS structures to be converted
  double F[],      // out   r×r×n, returns matrix parts of FS (unless F=0)
  double Fd[],     // out   r×r×nPar×n, returns derivative parts of FS
  int n,           // in    number of components in FS
  int nPar);       // in    number of parameters
#endif

