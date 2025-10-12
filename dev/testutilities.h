// Utilities used by tests of VarmaC
#ifndef TESTUTILITIES_H
#define TESTUTILITIES_H

typedef double fun(int n, double x[], double g[], int ipar[], double par[]);

double relabsdiff(double a[], double b[], int n);
int almostSame(double a, double b);
int almostEqual(double a[], double b[], int n);
int almostAllSame(double a[], int n);
void mean3(double X[], int m, int n, int k, int idim, double mu[]);
double mean(int n, double x[]);
void cov(char *transp, int n, int m, double X[], double C[]);
void lehmer(int n, double A[]);
void minij(int n, double A[]);
double difftest(fun f, int n, double x[], int ipar[], double par[]);
void checktestutilities(void);
void checkloglik(void);

#endif
