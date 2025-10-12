#ifndef EXAMPLEUTIL_H
#define EXAMPLEUTIL_H

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "RandomNumbers.h"

#define MISS_0     0
#define MISS_1ST   -1
#define MISS_LAST  -2
#define MISS_1ST4  -3
#define MISS_LAST4 -4
#define MISS_5A    -5
#define MISS_5B    -6
#define MISS_25A   -7

static inline int GetMissPat(const char *code) {
    const int vals[] = {0, -1, -2, -3, -4, -5, -6, -7};
    const char *keys[] = {"0", "1st", "last", "1st4", "last4", "5a","5b","25a"};
    for (int i=0; i<8; i++) if (strcmp(code, keys[i]) == 0) return vals[i];
    return atoi(code);
}

int MakeMissing ( // Make some observations missing, return missing value count
  int n,         // in   Number of observation times
  int r,         // in   Dimension of vector observed at each time
  int mpatt,     // in   Code (pattern) controlling which values to drop
  rand_rng *rng, // in   Random number generator
  int miss[]);    // out  r by n array; 1 where value is dropped and 0 elsewhere

double loglik ( // Evaluate minus log-likelihood and gradient for VARMA model
               // when model parameters are combined in a single vector.
  int     ntheta,  // in  parameter count
  double  theta[], // in  ntheta, (info on) model parameters: A, B, Sig, mu
  double  grad[],  // out nPar, gradient of log-likelihood (or null)
  int     ipar[],  // in  7 + r*n, additional integer information
  double  par[]);  // in  r*n + mJ*nJ, additional double information (X and J)
void theta2mat(int ntheta, double theta[], int ipar[], double par[], double A[],
               double B[], double Sig[], double mu[]);
void Xdim2par(double X[], int miss[], int p, int q, int r, int n, int M,
              int ipar[], double par[]);
void mat2theta(double A[], double B[], double Sig[], double mu[], int p, int q,
               int r, double theta[]);
void Jdim2par(double J[], int r, int n, int mJ, int nJ, int ipar[], double par[]);
void vech(int n, double A[], double avec[]); // Lower triangle as a vector

#endif
