#include <string.h>
#include <stdbool.h>
#include "xCheck.h"
#include "Tests.h"
#include "Testcase.h"
#define DEBUG
#include "printX.h"
#include "debugprint.h"
/*-------------------------------------------------------------*/

static FILE *err = 0;

static int contains(const char *hay, const char *needle) {
  return hay && needle && strstr(hay, needle) != 0;
}

static void check_max(void) {
  int pmax = 0, qmax = 0, rmax = 0, icase_max;
  bool ok;
  ok = Testcase(0, 0, 0, "max", &pmax, &qmax, &rmax, &icase_max, 0, err);
  xCheck(ok);
  // Realistic values:
  xCheck(icase_max >= 10);
  xCheck(pmax >= 1 && pmax <= 10);
  xCheck(qmax >= 1 && qmax <= 10);
  xCheck(rmax >= 1 && rmax <= 10);
}

static int get_max(void) {
  int pmax = 0, qmax = 0, rmax = 0, icase_max;
  Testcase(0, 0, 0, "max", &pmax, &qmax, &rmax, &icase_max, 0, err);
  return icase_max;
}

static void check_failures(void) {
  char name[16] = "";
  int p = 0, q = 0, r = 1, icase = 1;
  int mone = -1, zero = 0;
  double A[1] = {0.0}, B[1] = {0.0}, Sig[1] = {0.0};
  int maxcase = get_max();
  int icase_m = -2, icase_p = maxcase + 1;
  bool ok;
  // (1) A = B = Sig = 0, name = "", icase = 1 should be ok
  // but the following should fail:
  // (2) p, q or r null pointers
  // (3) icase null
  // (4) name and icase null
  // (5) some A, B, Sig null others not
  // (6) A, B, Sig != 0, invalid dimensions
  // (7) name "max": A, B, Sig not all null
  // (8) unknown name
  // (9) name empty: icase < -1
  // (10) name empty: icase > max-case
  ok = Testcase(0, 0, 0, name, &p, &q, &r, &icase, 0, err);    xCheck(ok);  // (1)
  strcpy(name, "");                                                
  ok = Testcase(A, B, 0, name, 0, &q, &r, &icase, 0, err);     xCheck(!ok); // (2)
  ok = Testcase(0, 0, 0, name, &p, &q, &r, 0, 0, err);         xCheck(!ok); // (3)
  ok = Testcase(0, 0, 0, 0, &p, &q, &r, 0, 0, err);            xCheck(!ok); // (4)
  ok = Testcase(0, B, 0, name, &p, &q, &r, &icase, 0, err);    xCheck(!ok); // (5)
  ok = Testcase(A, B, Sig, "", &p, &q, &zero, &zero, 0, err);  xCheck(!ok); // (6)
  ok = Testcase(A, B, Sig, "", &mone, &q, &r, &zero, 0, err);  xCheck(!ok); // (6)
  ok = Testcase(A, B, Sig, "max", &p, &q, &r, &icase, 0, err); xCheck(!ok); // (7)
  ok = Testcase(0, 0, 0, "blah", &p, &q, &r, &icase, 0, err);  xCheck(!ok); // (8)
  ok = Testcase(0, 0, 0, name, &p, &q, &r, &icase_m, 0, err);  xCheck(!ok); // (9)
  ok = Testcase(0, 0, 0, name, &p, &q, &r, &icase_p, 0, err);  xCheck(!ok); // (10)
}

static void check_single_case(int k) {
  char name[64] = "";
  bool ok;
  int p = 0, q = 0, r = 0, icase = k;
  ok = Testcase(0, 0, 0, name, &p, &q, &r, &icase, 0, err);
  xCheck(ok);
  xCheck(r > 0);
  if (contains(name, "ARMA"))    xCheck(p > 0 && q > 0);
  else if (contains(name, "MA")) xCheck(p == 0 && q > 0);
  else if (contains(name, "AR")) xCheck(p > 0 && q == 0);
  else xCheck(0);  // illegal name

  // re-query by name and confirm consistency
  int p2 = 0, q2 = 0, r2 = 0, icase2 = 0;
  ok = Testcase(0, 0, 0, name, &p2, &q2, &r2, &icase2, 0, err); xCheck(ok);
  xCheck(p2 == p && q2 == q && r2 == r && icase2 == icase);
}

static void check_dims(void) {
  int icase_max = get_max();
  for (int k = 1; k <= icase_max; ++k)
    check_single_case(k);
}

void check_inquiry(void) {
  check_max();
  check_failures();
  check_dims();
}

void TestTestcase(void) {
  //err = stderr;
  err = 0;
  check_inquiry();
}
