#ifndef PRINTX_H
#define PRINTX_H

static inline void printOff(void);
static inline void printOn(void);
static inline int  printIsOff(void);
static inline void printSetStatus(int noprint);
static inline void printSetNdec(int n);
static inline void printSetWidth(int w);
static inline void printSetFmtChar(char f);
static inline void printSetFmt(int w, int n, char f);
static inline int printGetNdec(void);

static inline void printNewl(void);                                  // print new line
static inline void printMsg(const char *s);                          // print message
static inline void printMsgUpper(const char *s);                     // print uppercased message
static inline void printMsgUnderl(const char *s);                    // print underlined message

static inline void printI(const char *name, int i);                  // print named int
static inline void printP(const char *name, const void *p);          // print named pointer
static inline void printD(const char *name, double d);               // print named double
static inline void printV(const char *name, const double x[], int n);// print named vector
static inline void printIV(const char *name, const int iv[], int n); // print named int vector
static inline void printS(const char *name, const char *s);          // print named string

static inline void printM(const char *name, const double A[], int nr, int nc);  // print matrix
static inline void printMP(const char *name, const double *ap, int nr, int nc,  // print submat
             const double A[], int ldA);                          // (ap points into A)
static inline void printMT(const char *name, const double A[], int nr, int nc); // print transpose

static inline void print2I(const char *s, int i1, int i2);          // print 2 integers
static inline void print3I(const char *s, int i1, int i2, int i3);  // print 3 integers
static inline void print4I(const char *s, int i1, int i2, int i3, int i4); // print 4 integers
static inline void print5I(const char *s, int i1, int i2, int i3, int i4, int i5); // 5 ints
static inline void print6I(const char *s, int i1, int i2, int i3, int i4, int i5,
                           int i6); // 6 ints

#include <math.h>
#include <string.h>
#include <ctype.h>

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

static inline void printOff(void)                    { noprint = 1; }
static inline void printOn(void)                     { noprint = 0; }
static inline void printSetStatus(int status)        { noprint = status; }
static inline int printIsOff(void)                   { return noprint; }
static inline void printSetNdec(int n)               { ndec = n; }
static inline void printSetWidth(int w)              { width = w; }
static inline void printSetFmtChar(char f)           { fmtc = f; }
static inline void printSetFmt(int w, int n, char f) { width = w; ndec = n; fmtc = f; }
static inline int printGetNdec(void)                 { return ndec; }

static inline void printNewl(void) { // print new line
  if (noprint) return;
  printf("\n");
  fflush(0);
}

static inline void printMsg(const char *s) { // print message (or string without name)
  if (noprint) return;
  printf("%s\n", s);
  fflush(0);
}

static inline void printMsgUpper(const char *s) {
  while (*s) putchar(toupper((unsigned char)*s++));
  putchar('\n');
  fflush(0);
}

static inline void printMsgUnderl(const char *s) { // print underlined message
  int i, nm = strlen(s);
  putchar('\n');
  printMsg(s);
  for (i=0; i<nm; i++) printf("%c", '-');
  putchar('\n');
  fflush(0);
}

static inline void printI(const char *name, int i) { // print named int
  if (noprint) return;
  printf("%s = %d\n", name, i);
  fflush(0);
}

static inline void printP(const char *name, const void *p) {
  if (noprint) return;
  printf("%s = %p\n",name, p);
  fflush(0);
}

static inline void printD(const char *name, double d) { // print named double
  char fmt[] = "%s = %1.*_\n", *uscore;
  if (noprint) return;
  uscore = strchr(fmt, '_'); *uscore = fmtc;
  printf(fmt,name,ndec,d);
  fflush(0);
}

static inline void printV(const char *name, const double x[], int n) { // print named vector
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

static inline void printIV(const char *name, const int iv[], int n) { // print named integer vector
  int j;
  if (noprint) return;
  printf("%s = ",name);
  for (j=0; j<n; j++) printf("%d ", iv[j]);
  printf("\n");
  fflush(0);
}

static inline void printS(const char *name, const char *s) { // print named string
  if (noprint) return;
  printf("%s = %s\n",name,s);
  fflush(0);
}

static inline void printMat(char *transp, const char *name, const double d[], int ldd,
		     int nr, int nc) {
  // local function to print matrix, called by printM and printMP
  //printf("In printMat, d[1]=%.3f\n", d[1]);
  char fmt[] = "%*.*_ ", *uscore;
  int i, j, is, ns, idx;
  if (noprint) return;
  uscore = strchr(fmt, '_'); *uscore = fmtc;
  ns = strlen(name);
  printf("%s = ",name);
  if (d==0) 
    printf("NULL\n");
  else {
    for (i=0; i<nr; i++) {
      if (i>0) for (is=0; is<ns+3; is++) printf(" ");
      for (j=0; j<nc; j++) {
        idx = *transp == 'N' ? j*ldd + i : i*ldd + j;
        //printf("In printMat, idx=%d\n", idx);
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

static inline void printM(const char *name, const double A[], int nr, int nc) { // print named matrix
  //printf("In printM, A[1]=%.3f\n", A[1]);
  printMat("N", name, A, nr, nr, nc);
}

static inline void printMT(const char *name, const double A[], int nr, int nc) { // print transpose
  //printf("In printMT, A[1]=%.3f\n", A[1]);
  printMat("T", name, A, nr, nc, nr);
}

static inline void printMP(const char *name, const double *ap, int nr, int nc,
             const double A[], int ldA) {
  // print submatrix
  char s[100];
  int r = (int)((ap-A) % ldA);
  int c = (int)((ap-A) / ldA);
  sprintf(s, "%s[%d:%d,%d:%d]", name, r, r+nr-1, c, c+nc-1);
  printMat("N", s, ap, ldA, nr, nc);
}

// Following functions are useful for debugging:
static inline void print2I(const char *s, int i1, int i2) { // print 2 integers
  if (noprint) return;
  printf("%s = %d, %d\n", s, i1, i2);
  fflush(0);
}

static inline void print3I(const char *s, int i1, int i2, int i3) { // print 3 integers
  if (noprint) return;
  printf("%s = %d, %d, %d\n", s, i1, i2, i3);
  fflush(0);
}

static inline void print4I(const char *s, int i1, int i2, int i3, int i4) { // print 4 integers
  if (noprint) return;
  printf("%s = %d, %d, %d, %d\n", s, i1, i2, i3, i4);
  fflush(0);
}

static inline void print5I(const char *s, int i1, int i2, int i3, int i4, int i5) { // print 5 ints
  if (noprint) return;
  printf("%s = %d, %d, %d, %d, %d\n", s, i1, i2, i3, i4, i5);
  fflush(0);
}

static inline void print6I(const char *s, int i1, int i2, int i3, int i4, int i5, int i6) { // 6 ints
  if (noprint) return;
  printf("%s = %d, %d, %d, %d, %d, %d\n", s, i1, i2, i3, i4, i5, i6);
  fflush(0);
}


#endif /* PRINTX_H */
