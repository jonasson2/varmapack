#include "error.h"
#include "varmapack.h"

static THREADLOCAL char last_error[256];

HIDDEN void clear_error(void) {
  last_error[0] = 0;
}

HIDDEN void set_error(char *message) {
  STRSET(last_error, message ? message : "internal error");
}

HIDDEN bool fail_error(char *message) {
  set_error(message);
  return false;
}

char *varmapack_last_error(void) {
  return last_error[0] ? last_error : 0;
}
