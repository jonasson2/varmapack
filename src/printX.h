#ifndef PRINTX_H
#define PRINTX_H
// PRINTX.H SIMPLE PRINTING PACKAGE FOR C-FROM-MATLAB
// (primarily intended for debugging numerical methods)

void printOff(void);                     // turn off printing
void printOn(void);                      // turn on printing
int printIsOff(void);                    // return 1 (true) if printing is off
void printSetNdec(int n);                // set number of fraction decimals
void printSetWidth(int w);               // set field width for matrix printing
void printSetFmtChar(char f);            // set to 'f' or 'e' format descriptor
void printSetFmt(int w, int n, char f);  // set width, ndec and fmt-char
void printNewl(void);                    // Print newline
void printMsg(char *s);                                      // print message
void printMsgUnderl(char *s);                                // print underlined msg
void printI(char *name, int i);                              // print int
void printP(char *name, void *p);                            // print pointer
void printD(char *name, double d);                           // print double
void printV(char *name, double x[], int n);                  // print vector
void printS(char *name, char s[]);                           // print string
void printM(char *name, double A[], int nr, int nc);         // print matrix
void printMP(char *name,double *ap,int nr, int nc, double A[], int ldA); //print
                                                             // matrix partition
void printIV(char *name, int iv[], int n);                   // print intvector
void print2I(char *s, int i1, int i2);                        // print 2 ints
void print3I(char *s, int i1, int i2, int i3);                // print 3 ints
void print4I(char *s, int i1, int i2, int i3, int i4);        // print 4 ints
void print5I(char *s, int i1, int i2, int i3, int i4, int i5);// print 5 ints
void print6I(char *s, int i1, int i2, int i3, int i4, int i5, int i6); // 6 ints

// NOTES:
// default double and vector format is "%.3f"
// default matrix format is "%7.3f"
// incx = increment (used to print matrix rows assuming column major order)
// ldA = leading dimension of A
// nr & nc = number of rows & columns
// ap is a pointer to upper left corner of nr×nc partition (i.e. submatrix) in A

#endif
