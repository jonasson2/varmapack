#ifndef XCHECK_H
#define XCHNECK_H

#include <stdio.h>
#include <string.h>
#include "xCheck.h"

static int FAIL = 0;
static char commonmsg[60] = "";
static char addmsg[60] = "";

void xCheckFunc(char *message, char* file, int line) {
  // Print filename, line number, and upto three messages.
  fprintf(stderr, "Test failed on line %d in file %s\n", line, file);
  if (commonmsg[0] != 0) fprintf(stderr, "%s\n", commonmsg);
  if (addmsg[0] != 0) fprintf(stderr, "%s\n", addmsg);
  fprintf(stderr, "%s\n", message);
  FAIL = 1;
}

void xCheckInit(char *msg) { // set common message to msg and FAIL to false
  strcpy(commonmsg, msg);
  addmsg[0] = 0;
  FAIL = 0;
}

void xCheckAddMsg(char *amsg) { // set additional message to amsg
  strcpy(addmsg, amsg);
}

int xCheckFailed(void) {
  return FAIL;
}

#endif
