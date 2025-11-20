#include <stdio.h>
#include <stdbool.h>
#include <string.h>

#include "xAssert.h"
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

void TestXAssert(void) {
  FILE *tmp = tmpfile();
  xCheck(tmp != 0);
  bool ok = fail_helper(tmp);
  xCheck(!ok);
  fflush(tmp);
  rewind(tmp);
  char buf[128] = {0};
  size_t n = fread(buf, 1, sizeof(buf)-1, tmp);
  xCheck(n > 0);
  xCheck(strstr(buf, "error in TestXAssert.c") != 0);
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
}
