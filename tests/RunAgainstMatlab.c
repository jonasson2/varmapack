#include <stdbool.h>
#include "printX.h"
#include "Tests.h"
#include "xCheck.h"

int main(void) {
  bool ok;
  printOff();
  xCheckInit("AgainstMatlab");
  ok = TestAgainstMatlab();
  xCheck(ok);
  return xCheckNFailures() > 0;
}
