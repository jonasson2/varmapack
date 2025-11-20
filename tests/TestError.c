#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include "error.h"
#include "xCheck.h"

static bool fail_helper(FILE *stream) {
  varmapack_set_errstream(stream);
  XASSERT(false, "forced failure %d", 42);
  return true;
}

static bool pass_helper(FILE *stream) {
  varmapack_set_errstream(stream);
  XASSERT(true, "should not trigger");
  return true;
}

static bool alloc_helper(int len) {
  double *buf = 0;
  ALLOC(buf, len);
  free(buf);
  return true;
}

void TestError(void) {
  FILE *tmp = tmpfile();
  xCheck(tmp != 0);
  bool ok = fail_helper(tmp);
  xCheck(!ok);
  fflush(tmp);
  rewind(tmp);
  char buf[128] = {0};
  size_t n = fread(buf, 1, sizeof(buf)-1, tmp);
  xCheck(n > 0);
  xCheck(strstr(buf, "forced failure 42") != 0);
  fclose(tmp);

  tmp = tmpfile();
  xCheck(tmp != 0);
  ok = pass_helper(tmp);
  xCheck(ok);
  fflush(tmp);
  rewind(tmp);
  n = fread(buf, 1, sizeof(buf)-1, tmp);
  xCheck(n == 0);
  fclose(tmp);

  varmapack_set_errstream(0);
  ok = fail_helper(0);
  xCheck(!ok);

  FILE *tmp2 = tmpfile();
  xCheck(tmp2 != 0);
  varmapack_set_errstream(tmp2);
  ok = alloc_helper(0);
  xCheck(ok);
  ok = alloc_helper(16);
  xCheck(ok);
  fclose(tmp2);
}
