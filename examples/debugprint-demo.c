#define DEBUG
#include "debugprint.h"
// Not a very comprehensive test.
// Try compiling with and without -DNODEBUG.
int main(void) {
  int a = 64, b=-55, c=888, d=9999;
  size_t s = 12345678901234567890LLU;
  double x[] = {1,1.333, 3.14,2.45, 1.22,3.41, 8.22};
  double y[] = {100,11.333,-33.14,2.45,1.22,3.41,18.22};
  print0();
  printh("This is a string");
  print(d,a);
  print2(d,a,b);
  print3(u,a,b,c);
  print4(d,a,b,c,d);
  print5(d,a,b,c,d,a);
  print(zu,s);
  printv(.4f,x,6);
  printnl();
  printh("printm:");
  printm(.2f,x,2,3,2);
  printh("printmR:");
  printmR(.2f,x,2,3,3);
  printh("printmL:");
  printmL(4.2f,x,2,2,2);
  printh("printmU:");
  printmU(4.2f,x,2,2,2);
  printnl();
  printmR(6.2f,y,2,2,2);
  debugprint("x[0:2] = %4.2f %4.2f\n",x[0],x[1]);
  return 0;
}
