#ifndef PRINTX_H
#define PRINTX_H

void printOff(void);
void printOn(void);
int  printIsOff(void);
void printSetStatus(int noprint);
void printSetNdec(int n);
void printSetWidth(int w);
void printSetFmtChar(char f);
void printSetFmt(int w, int n, char f);
int printGetNdec(void);

void printNewl(void);                                  // print new line
void printMsg(const char *s);                          // print message
void printMsgUpper(const char *s);                     // print uppercased message
void printMsgUnderl(const char *s);                    // print underlined message

void printI(const char *name, int i);                  // print named int
void printP(const char *name, const void *p);          // print named pointer
void printD(const char *name, double d);               // print named double
void printV(const char *name, const double x[], int n);// print named vector
void printIV(const char *name, const int iv[], int n); // print named int vector
void printS(const char *name, const char *s);          // print named string

void printM(const char *name, const double A[], int nr, int nc);  // print matrix
void printMP(const char *name, const double *ap, int nr, int nc,  // print submat
             const double A[], int ldA);                          // (ap points into A)
void printMT(const char *name, const double A[], int nr, int nc); // print transpose

void print2I(const char *s, int i1, int i2);          // print 2 integers
void print3I(const char *s, int i1, int i2, int i3);  // print 3 integers
void print4I(const char *s, int i1, int i2, int i3, int i4); // print 4 integers
void print5I(const char *s, int i1, int i2, int i3, int i4, int i5); // 5 ints
void print6I(const char *s, int i1, int i2, int i3, int i4, int i5,
                           int i6); // 6 ints

#endif /* PRINTX_H */
