#ifndef XCHECK_H
#define XCHECK_H

#define xCheck(ok, msg)  ((ok) ? (void)0 : xCheckFunc(msg, __FILE__, __LINE__))
  // xCheck(condition, message) returns silently if condition is true, but if it
  // is false, then message is printed on stderr along with the source file name
  // and line number where the check was made and the common message set with 
  // xCheckInit. In addition, when condition is false the satic variable FAIL
  // (defined in xCheck.c) is set to true.

void xCheckFunc(char *message, char* file, int line);
  // The macro CHECK is replaced by a call to this function.
void xCheckInit(char *msg);
  // Initialize static variables commonmsg and FAIL of xCheck.c to msg and false
void xCheckAddMsg(char *amsg);
  // set additional message to amsg
int xCheckFailed(void);
  // Return true if FAIL is true.

#endif
