// Assert that works from within Matlab and (maybe) R. 
#ifndef XASSERT_H
#define XASSERT_H

#include <string.h>
#include <stdarg.h>

#define xAssert(e)  ((e) ? (void)0 : xPrintAssertion(#e, __FILE__, __LINE__))
// xAssert(expression) evaluates expression, and if it is 0 (false), a message
// containing the failed expression (literally) is printed along with the source
// file name and line number where the assertion was made, and the execution of
// the program is stopped. For example, xAssert(2+2==3) on line 12 in file tst.c
// will cause program termination after printing (something like): "Assertion
// 2+2==3 failed on line 12 in tst.c".

#ifdef MEX // FROM MATLAB MEX FILE
#include "mex.h"
//extern void mexErrMsgTxt(char *msg);
static inline void xErrorExit(char *message) {
  // Print message to console window and exit to Matlab prompt
  mexErrMsgTxt(message);
}

#elif defined(USING_R) // FROM R
#include <R.h>  // Provides Rf_error
static inline void xErrorExit(char *message) {
  // Raise an error in R and return control to the R interpreter
  Rf_error("%s", message);
}

#else // NOT FROM MATLAB OR R
#include <stdio.h>
#include <stdlib.h>
static inline void xErrorExit(char *message) {
  // Print message to stderr and exit with a status of 1
  fprintf(stderr, "%s\n", message);
  exit(1);
}
#endif

static inline void xPrintAssertion(char* assertion, char* file, int line) {
  // The macro xAssert is replaced by a conditional call to this function, which (if
  // called) prints (literally) the expression which caused the assertion failure and
  // stops program execution.
  char msg[200]; // should be large enough
  char m1[] = "Assertion failed:";
  char m2[] = ", at line";
  char m3[] = "of file";
  int nmax, n;
  nmax = 200 - (strlen(file) + strlen(m1) + strlen(m2) + 25 + strlen(m3));
    // 25 >= largest int length + spaces and '...' in msg
  n = strlen(assertion);
  if (n <= nmax)
    sprintf(msg,"%s %s%s %d %s %s", m1, assertion, m2, line, m3, file);
  else {
    sprintf(msg, "%s %.*s...%s %d %s %s", m1, nmax, assertion, m2, line, m3, file);
  }
  xErrorExit(msg);
}

static inline void xAssertMessage(int expression, char *message) {
  // xAssertMessage(expression, "message") checks that expression is non-zero (true),
  // and if not the supplied message is printed and program execution is stopped.
  // For example, xAssertMessage(2+2==3, "No, 2+2 is not 3.") will terminate the 
  // program after printing "No, 2+2 is not 3.".

  // The reason for not using C's built in assert macro is to enable use from
  // environments where special functions must be called to terminate execution,
  // such as when called from Matlab mex files. Also, the xAssert macro can be
  // replaced with an empty expansion, if the program is error-free.
  if (!expression) xErrorExit(message);
}

#endif
