#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include "error.h"

static bool assert_demo(void) {
  ASSERT(2 + 2 == 5, "math is hard: %d", 42);
  return true;
}

static bool alloc_demo(void) {
  double *x;
  ALLOC(x, INT64_MAX);
  ASSERT(2 + 2 == 5, "math is hard: %d", 42);
  return true;
}

int main(void) {
  bool ok;
  ok = assert_demo();
  printf("assert_demo -> %s\n", ok ? "true" : "false");
  ok = alloc_demo();
  printf("alloc_demo -> %s\n", ok ? "true" : "false");
  return 0;
}
