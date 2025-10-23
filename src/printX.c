#include "printX.h"
#include <math.h>
#include <string.h>

#ifdef MEX
#include "mex.h"
#define printf mexPrintf
extern int mexPrintf(const char *fmt, ...);
#define fflush(x)

#elif defined(USING_R)
#include <R.h>  // For Rprintf
#define printf Rprintf
#define fflush(x)  // No-op, as Rprintf automatically flushes

#else
#include <stdio.h>
#endif

static int noprint = 0;
static int width = 7;
static int ndec = 3;
static char fmtc = 'f';

void printOff(void)                    { noprint = 1; } // useful for debugging
void printOn(void)                     { noprint = 0; } // useful for debugging
int printIsOff(void)                   { return noprint; }
void printSetNdec(int n)               { ndec = n; }
void printSetWidth(int w)              { width = w; }
void printSetFmtChar(char f)           { fmtc = f; }
void printSetFmt(int w, int n, char f) { width = w; ndec = n; fmtc = f; }

void printNewl(void) { // print new line
  if (noprint) return;
  printf("\n");
  fflush(0);
}

void printMsg(char *s) { // print message (or string without name)
  if (noprint) return;
  printf("%s\n", s);
  fflush(0);
}

void printMsgUnderl(char *s) { // print underlined message
  int i, nm = strlen(s);
  printMsg(s);
  for (i=0; i<nm; i++) printf("%c", '-');
  printMsg("");
  fflush(0);
}

void printI(char *name, int i) { // print named int
  if (noprint) return;
  printf("%s = %d\n", name, i);
  fflush(0);
}

void printP(char *name, void *p) {
  if (noprint) return;
  printf("%s = %p\n",name, p);
  fflush(0);
}

void printD(char *name, double d) { // print named double
  char fmt[] = "%s = %1.*_\n", *uscore;
  if (noprint) return;
  uscore = strchr(fmt, '_'); *uscore = fmtc;
  printf(fmt,name,ndec,d);
  fflush(0);
}

void printV(char *name, double x[], int n) { // print named vector
  char fmt[] = "%1.*_ ", *uscore;
  int j;
  if (noprint) return;
  uscore = strchr(fmt, '_'); *uscore = fmtc;
  printf("%s = ",name);
  if (x==0)
    printf("NULL");
  else
    for (j=0; j<n; j++) printf(fmt, ndec, x[j]);
  printf("\n");
  fflush(0);
}

void printIV(char *name, int iv[], int n) { // print named integer vector
  int j;
  if (noprint) return;
  printf("%s = ",name);
  for (j=0; j<n; j++) printf("%d ", iv[j]);
  printf("\n");
  fflush(0);
}

void printS(char *name, char *s) { // print named string
  if (noprint) return;
  printf("%s = %s\n",name,s);
  fflush(0);
}

// Stop cl from griping about size_t to int conversion
static int len(char *s) { return (int)strlen(s); }

static void printMat(char *name, double d[], int ldd, int nr, int nc) {
  // local function to print matrix, called by printM and printMP
  char fmt[] = "%*.*_ ", *uscore;
  int i, j, is, ns, idx;
  if (noprint) return;
  uscore = strchr(fmt, '_'); *uscore = fmtc;
  ns = len(name);
  printf("%s = ",name);
  if (d==0) 
    printf("NULL\n");
  else {
    for (i=0; i<nr; i++) {
      if (i>0) for (is=0; is<ns+3; is++) printf(" ");
      for (j=0; j<nc; j++) {
        idx = j*ldd + i;
        if (fmtc == 'f') { // when fmt is f print < 1e-14 as 0
          if (fabs(d[idx]) >= 1e-14 && fmtc == 'f')
            printf("%*.*f ", width, ndec, d[idx]);
          else
            printf("%*.0f ", width, d[idx]);
        }
        else {
          printf(fmt, width, ndec, d[idx]);
        }
      }
      printf("\n");
    }
    if (nr==0) printf("\n");
  }
  fflush(0);
}

void printM(char *name, double A[], int nr, int nc) { // print named matrix
  printMat(name, A, nr, nr, nc);
}

void printMP(char *name, double *ap, int nr, int nc, double A[], int ldA) {
  // print submatrix
  char s[100];
  int r = (int)((ap-A) % ldA);
  int c = (int)((ap-A) / ldA);
  sprintf(s, "%s[%d:%d,%d:%d]", name, r, r+nr-1, c, c+nc-1);
  printMat(s, ap, ldA, nr, nc);
}

// Following functions are useful for debugging:
void print2I(char *s, int i1, int i2) { // print 2 integers
  if (noprint) return;
  printf("%s = %d, %d\n", s, i1, i2);
  fflush(0);
}

void print3I(char *s, int i1, int i2, int i3) { // print 3 integers
  if (noprint) return;
  printf("%s = %d, %d, %d\n", s, i1, i2, i3);
  fflush(0);
}

void print4I(char *s, int i1, int i2, int i3, int i4) { // print 4 integers
  if (noprint) return;
  printf("%s = %d, %d, %d, %d\n", s, i1, i2, i3, i4);
  fflush(0);
}

void print5I(char *s, int i1, int i2, int i3, int i4, int i5) { // print 5 ints
  if (noprint) return;
  printf("%s = %d, %d, %d, %d, %d\n", s, i1, i2, i3, i4, i5);
  fflush(0);
}

void print6I(char *s, int i1, int i2, int i3, int i4, int i5, int i6) { // 6 ints
  if (noprint) return;
  printf("%s = %d, %d, %d, %d, %d, %d\n", s, i1, i2, i3, i4, i5, i6);
  fflush(0);
}
