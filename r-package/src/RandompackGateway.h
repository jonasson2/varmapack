#ifndef RANDOMPACKGATEWAY_H
#define RANDOMPACKGATEWAY_H

#ifdef USING_R
#include "RandompackRGateway.h"
#elif defined(USING_PYTHON)
#include "RandompackPythonGateway.h"
#else
#include "randompack.h"
#endif

#endif // RANDOMPACKGATEWAY_H
