#ifndef FROMMATLAB_H
#define FROMMATLAB_H

#include <stdio.h>
#include <stdbool.h>

bool IntFromMatlab(FILE *f, char *var, int *val);
bool DoubleFromMatlab(FILE *f, char *var, double *val);
bool IntVecFromMatlab(FILE *f, char *vector, int *x, int n);
bool VectorFromMatlab(FILE *f, char *vector, double *x, int n);
bool MatrixFromMatlab(FILE *f, char *matrix, double A[], int m, int n);
bool CompareWithMatlab(FILE *f, char *matrix, double *A, int m, int n, double *diff);

#endif
