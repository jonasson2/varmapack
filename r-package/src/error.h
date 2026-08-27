// Assert helpers usable from Matlab, R, or standalone C.
#ifndef ASSERT_H
#define ASSERT_H

#include "varmapack_config.h"

HIDDEN void clear_error(void);
HIDDEN void set_error(char *message);
HIDDEN bool fail_error(char *message);

#define xAssert(e)  ((e) ? (void)0 : xPrintAssertion(#e, __FILE__, __LINE__))

#define ASSERT(cond, ...)				\
 do {							\
   if (!(cond)) {					\
     fprintf(stderr, "error in %s:%d: ", basename_only(__FILE__), __LINE__); \
     fprintf(stderr, __VA_ARGS__);			\
     fprintf(stderr, "\n");				\
     abort();						\
   }							\
 } while (0)

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
  char *ellipsis = "";
  int fileLength = 75;
  int assertionLength;
  if (strlen(file) < (size_t)fileLength) fileLength = (int)strlen(file);
  assertionLength = 148 - fileLength;
  if (strlen(assertion) < (size_t)assertionLength) {
    assertionLength = (int)strlen(assertion);
  }
  else {
    ellipsis = "...";
  }
  snprintf(msg, sizeof(msg), "%s %.*s%s%s %d %s %.*s", m1, assertionLength,
           assertion, ellipsis, m2, line, m3, fileLength, file);
  xErrorExit(msg);
}

static inline const char *basename_only(const char *path) {
  if (!path) return "";
  const char *slash = strrchr(path, '/');
#ifdef _WIN32
  const char *bslash = strrchr(path, '\\');
  if (!slash || (bslash && bslash > slash)) slash = bslash;
#endif
  return slash ? slash + 1 : path;
}

#endif
