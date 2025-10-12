#ifndef ALLOCATE_H
#define ALLOCATE_H
#include <stdlib.h>
#include "xAssert.h"
// Following are declared in stdlib.h and need (should) not be redeclared here
#define freem(p) free(p)
#define allocate(p, len) {p=calloc(len,sizeof(*(p)));xAssert((p)!=0);}
// Following are to be ignored:
#define checkStatus()
#define checkLeak()
#endif
