// Assert that works from within Matlab and (maybe) R. 
#ifndef XASSERT_H
#define XASSERT_H

#define xAssert(e)  ((e) ? (void)0 : xPrintAssertion(#e, __FILE__, __LINE__))
// xAssert(expression) evaluates expression, and if it is 0 (false), a message
// containing the failed expression (literally) is printed along with the source
// file name and line number where the assertion was made, and the execution of
// the program is stopped. For example, xAssert(2+2==3) on line 12 in file tst.c
// will cause program termination after printing (something like): "Assertion
// 2+2==3 failed on line 12 in tst.c".

void xErrorExit(char *message);
// xErrorExit("message") prints message to console window (or stderr) and stops
// program execution.
  
void xPrintAssertion(char* assertion, char* file, int line);
// The macro xAssert is replaced by a conditional call to this function, which
// (if called) prints a failed assertion (literally) and terminates the program.

void xAssertMessage(int expression, char *message);
// xAssertMessage(expression, "message") checks that expression is non-zero (true),
// and if not the supplied message is printed and program execution is stopped.
// For example, xAssertMessage(2+2==3, "No, 2+2 is not 3.") will terminate the 
// program after printing "No, 2+2 is not 3.".

// The reason for not using C's built in assert macro is to enable use from
// environments where special functions must be called to terminate execution,
// such as when called from Matlab mex files. Also, the xAssert macro can be
// replaced with an empty expansion, if the program is error-free.

void option_error(const char *fmt, ...);
// Prints option error message formatted by fmt and exits

#endif
