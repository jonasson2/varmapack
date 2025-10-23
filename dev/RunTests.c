#include <getopt.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include "printX.h"
#include "Tests.h"
#include "xAssert.h"
#include "xCheck.h"

static int VERBOSE_MODE = 0, NTOTAL = 0, NFAIL = 0;
static const char *headr_fmt = "%-15s %8s %8s\n";
static const char *table_fmt = "%-15s %8d %8d\n";

static void print_help(void) {
  puts("Usage: RunTests [options]\n"
       "Options:\n"
       "  -h  Show this help message\n"
       "  -v  Verbose tests\n"
       );
}

static void vprint(const char *fmt, ...) {
  if (!VERBOSE_MODE) return;
  va_list ap;
  va_start(ap, fmt);
  vprintf(fmt, ap);
  va_end(ap);
  fflush(stdout);
}

static void run_test(const char *name, void (*fn)(void)) {
  int ntotal, nfail;
  xCheckInit(name);
  fn();
  ntotal = xCheckNTotal();
  nfail  = xCheckNFailures();
  NTOTAL += ntotal;
  NFAIL  += nfail;
  vprint(table_fmt, name, ntotal - nfail, nfail);
}

int main(int argc, char **argv) {
  char optstring[10] = ":vh", c;
  while ((c = getopt(argc, argv, optstring)) != -1) {
    switch (c) {
      case 'h': print_help(); return 0;
      case 'v': VERBOSE_MODE = 1; break;
      case '?': option_error("Unknown option -%c", optopt);
    }
  }
  if (optind < argc) option_error("Unexpected argument %s", argv[optind]);
  vprint(headr_fmt, "TEST OF", "PASSED", "FAILED");
  run_test("ExtraUtil",      TestExtraUtil);
  run_test("RandomNumbers",  TestRandomNumbers);
  run_test("Testcase",       TestTestcase);
  vprint(table_fmt, "TOTAL", NTOTAL - NFAIL, NFAIL);
  return (NFAIL > 0);
}

#include <stdio.h>
#include <stdarg.h>
#include <string.h>
