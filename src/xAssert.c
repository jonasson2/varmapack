#include <string.h>
#include <stdarg.h>
#include "xAssert.h"

// For more detailed description of the functions (and description of the
// xAssert macro) see comments in xAssert.h
#ifdef MEX // FROM MATLAB MEX FILE
#include "mex.h"
extern void mexErrMsgTxt(char *msg);
void xErrorExit(char *message) {
  // Print message to console window and exit to Matlab prompt
  mexErrMsgTxt(message);
}

#elif defined(USING_R) // FROM R
#include <R.h>  // Provides Rf_error
void xErrorExit(char *message) {
  // Raise an error in R and return control to the R interpreter
  Rf_error("%s", message);
}

#else // NOT FROM MATLAB OR R
#include <stdio.h>
#include <stdlib.h>
void xErrorExit(char *message) {
  // Print message to stderr and exit with a status of 1
  fprintf(stderr, "%s\n", message);
  exit(1);
}
#endif

// Stop cl from griping about size_t to int conversion
static int len(char *s) { return (int)strlen(s); }

void xPrintAssertion(char* assertion, char* file, int line) {
  // Print (literally) the expression which caused assertion failure and stop
  // program execution. This function is called by the xAssert macro defined in
  // xAssert.h
  char msg[200]; // should be large enough
  char m1[] = "Assertion failed:";
  char m2[] = ", at line";
  char m3[] = "of file";
  int nmax, n;
  nmax = 200 - (len(file) + len(m1) + len(m2) + 25 + len(m3));
    // 25 >= largest int length + spaces and '...' in msg
  n = len(assertion);
  if (n <= nmax)
    sprintf(msg,"%s %s%s %d %s %s", m1, assertion, m2, line, m3, file);
  else {
    sprintf(msg, "%s %.*s...%s %d %s %s", m1, nmax, assertion, m2, line, m3, file);
  }
  xErrorExit(msg);
}

void xAssertMessage(int expression, char *message) {
  // Assert that expression is true; if not print message and stop program.
  if (!expression) xErrorExit(message);
}

void option_error(const char *fmt, ...) {
    char msg_opt[256];
    va_list args;
    va_start(args, fmt);
    vsnprintf(msg_opt, sizeof(msg_opt), fmt, args);
    va_end(args);
    xErrorExit(msg_opt);
}
