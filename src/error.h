// Assert helpers usable from Matlab, R, or standalone C.
#ifndef ASSERT_H
#define ASSERT_H

#include <string.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#define xAssert(e)  ((e) ? (void)0 : xPrintAssertion(#e, __FILE__, __LINE__))

#define ASSERT(cond, ...)				\
 do {							\
   if (!(cond)) {					\
     error_report(__FILE__, __LINE__, __VA_ARGS__);	\
     return false;					\
   }							\
 } while (0)

#define ALLOC(ptr, count)					\
 do {								\
   if ((count) < 0) {						\
     error_report(__FILE__, __LINE__,				\
		  "negative allocation request for %s", #ptr);	\
     return false;						\
   }								\
   (ptr) = calloc((size_t)(count), sizeof(*(ptr)));		\
   if ((count) > 0 && !(ptr)) {					\
     error_report(__FILE__, __LINE__,				\
		  "allocation failed for %s (len=%zu)",		\
		  #ptr, (size_t)(count));			\
     return false;						\
   }								\
 } while (0)

#define FREE(p) free(p)

void varmapack_set_errstream(FILE *stream);
void error_report(const char *file, int line, const char *fmt, ...);

#ifdef MEX
#include "mex.h"
static inline void xErrorExit(char *message) { mexErrMsgTxt(message); }
#elif defined(USING_R)
#include <R.h>
static inline void xErrorExit(char *message) { Rf_error("%s", message); }
#else
static inline void xErrorExit(char *message) {
  fprintf(stderr, "%s\n", message);
  exit(1);
}
#endif

static inline void xPrintAssertion(char *assertion, char *file, int line) {
  char msg[200];
  char m1[] = "Assertion failed:";
  char m2[] = ", at line";
  char m3[] = "of file";
  int nmax = 200 - (strlen(file) + strlen(m1) + strlen(m2) + 25 + strlen(m3));
  int n = (int)strlen(assertion);
  if (n <= nmax) {
    sprintf(msg, "%s %s%s %d %s %s", m1, assertion, m2, line, m3, file);
  }
  else {
    sprintf(msg, "%s %.*s...%s %d %s %s", m1, nmax, assertion, m2, line, m3, file);
  }
  xErrorExit(msg);
}

static inline void xAssertMessage(int expression, char *message) {
  if (!expression) xErrorExit(message);
}

#endif
