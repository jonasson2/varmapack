#include "Tests.h"
#include "xAssert.h"
#include <getopt.h>
#include <stdio.h>
#include <string.h>
#include "printX.h"

void print_help(void) {
  puts("Usage: RunVarmaSim [options]\n"
       "Options:\n"
       "  -h, --help    Show this help message\n"
       "  -f            Run full tests\n"
       "  -p            Print name of each testcase");
}

int main(int argc, char **argv) {
  int fail;
  char optstring[10] = ":fph", code[3] = "", c;
  while ((c = getopt(argc, argv, optstring)) != -1) {
    switch (c) {
      case 'h': print_help(); return 0;
      case 'p': strcat(code, "p"); break;
      case 'f': strcat(code, "f"); break;
      case '?': option_error("Unknown option -%c", optopt);
    }
  }
  if (optind < argc) option_error("Unexpected argument %s", argv[optind]);
  printf("Testing VarmaSim...\n");
  fflush(0);
  printS("code", code);
  fail = TestVarmaSim(code);
  if (fail) {
    printf("\n  failed\n");
    return 1;
  }
  else {
    printf("\n  OK\n");
    return 0;
  }
}         
