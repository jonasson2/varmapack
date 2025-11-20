#include <stdio.h>
#include <stdarg.h>
#include <stdbool.h>
#include <string.h>
#include "error.h"

static FILE *error_stream = 0;
static bool error_stream_explicit = false;

void varmapack_set_errstream(FILE *stream) {
  error_stream = stream;
  error_stream_explicit = true;
}

static void error_init(void) {
  if (!error_stream_explicit && error_stream == 0) {
    error_stream = stderr;
  }
}

static const char *basename_only(const char *path) {
  if (!path) return "";
  const char *slash = strrchr(path, '/');
#ifdef _WIN32
  const char *bslash = strrchr(path, '\\');
  if (!slash || (bslash && bslash > slash)) slash = bslash;
#endif
  return slash ? slash + 1 : path;
}

void error_report(const char *file, int line, const char *fmt, ...) {
  error_init();
  if (error_stream == 0) return;
  va_list ap;
  va_start(ap, fmt);
  fprintf(error_stream, "error in %s:%d: ", basename_only(file), line);
  vfprintf(error_stream, fmt, ap);
  fprintf(error_stream, "\n");
  va_end(ap);
}
