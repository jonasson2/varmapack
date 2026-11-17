#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include "error.h"
#include "getopt.h"
#include "printX.h"
#include "Tests.h"
#include "TestUtil.h"

static int NTOTAL = 0, NFAIL = 0;
static char headr_fmt[] = "%-24s %-39s %8s %8s\n";
static char table_fmt[] = "%-24s %-39s %8d %8d\n";
static bool verbose = false;

static void print_help(void) {
  puts("Usage: RunTests [options]\n"
       "Options:\n"
       "  -h    Show this help message\n"
       "  -C    Comprehensive statistics: 200000 replicates (default 10000)\n"
       "  -v    Show test groups, coverage, and check counts\n");
}

static void vprint(const char *fmt, ...) {
  if (!verbose) return;
  va_list ap;
  va_start(ap, fmt);
  vprintf(fmt, ap);
  va_end(ap);
  fflush(stdout);
}

static void option_error(const char *fmt, ...) {
  // Prints option error message formatted by fmt and exits
  char msg_opt[256];
  va_list args;
  va_start(args, fmt);
  vsnprintf(msg_opt, sizeof(msg_opt), fmt, args);
  va_end(args);
  xErrorExit(msg_opt);
}

static void run_test(const char *name, const char *coverage,
                     void (*fn)(void)) {
  int ntotal, nfail;
  xCheckInit(name);
  fn();
  ntotal = xCheckNTotal();
  nfail  = xCheckNFailures();
  NTOTAL += ntotal;
  NFAIL  += nfail;
  vprint(table_fmt, name, coverage, ntotal - nfail, nfail);
}

int main(int argc, char **argv) {
  char optstring[] = ":vCh";
  int c;
  while ((c = getopt(argc, argv, optstring)) != -1) {
    switch (c) {
      case 'h': print_help(); return 0;
      case 'v': verbose = true; break;
      case 'C': TEST_NREP = 20*TEST_NREP_DEFAULT; break;
      case '?': option_error("Unknown option -%c", optopt);
    }
  }
  if (optind < argc) option_error("Unexpected argument %s", argv[optind]);
  printOff();
  vprint("\n");
  vprint("Statistical tests: %d independent replicates per model\n", TEST_NREP);
  vprint(headr_fmt, "TEST GROUP", "COVERAGE", "PASSED", "FAILED");
  run_test("FromMatlab helpers", "Matlab fixture parsing", TestFromMatlab);
  run_test("FindC/FindG", "Startup matrices and recursions", TestFindCFindG);
  run_test("varmapack_testcase", "Named VARMA model construction", TestTestcase);
  run_test("varmapack_testcasex", "VARMAX test-data construction", TestTestcasex);
  run_test("spectral radius", "AR and MA spectral radii", Testvarmapack_specrad);
  run_test("varmapack_acvf", "Theoretical autocovariances", TestAcvf);
  run_test("varmapack_autocov", "Sample autocovariances", TestAutocov);
  run_test("autocov edge cases", "Autocovariance error handling",
           TestAutocovEdgeCases);
  run_test("varmapack_cov2corr", "Covariance-to-correlation conversion",
           TestCov2corr);
  run_test("varmapack_cov2sig", "Inverse VYW and covariance validation",
           TestCov2sig);
  run_test("varmapack_sim", "VARMA simulation and errors", TestSimEdgeCases);
  run_test("OU simulation", "Stationary and lagged covariances", TestOU);
  run_test("Gamma times", "Gap moments and validation", TestGammaTimes);
  run_test("varmapack_simx", "VARMAX simulation and errors", TestSimxEdgeCases);
  run_test("PSD cond covariance", "Rank-deficient startup covariance",
           TestPsdCondCov);
  run_test("simulation statistics", "Means, covariances, and normality",
           TestDistribution);
  run_test("impulse responses", "Psi coefficients and orthogonalized IRFs",
           TestPsi);
  run_test("Lyapunov solver", "Internal equation solver", TestLyapunov);
  vprint(table_fmt, "TOTAL", "", NTOTAL - NFAIL, NFAIL);
  return (NFAIL > 0);
}
