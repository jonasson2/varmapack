#ifndef TESTMATLABMATRIX_H
#define TESTMATLABMATRIX_H

#include <stdbool.h>

extern char matmatfolder[];
bool TestMatlabMatrix(char *filename, double *A, // Compares A with the matrix in the
                      int m, int n);             // specified file.

#endif
