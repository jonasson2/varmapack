#ifndef RANDOMPACKPYTHONGATEWAY_H
#define RANDOMPACKPYTHONGATEWAY_H

#include "varmapack_config.h"
#include <randompack_python_api.h>

HIDDEN int randompack_PY_initialize(void);
HIDDEN randompack_rng *randompack_PY_create(const char *engine);
HIDDEN void randompack_PY_free(randompack_rng *rng);
HIDDEN char *randompack_PY_last_error(randompack_rng *rng);
HIDDEN bool randompack_PY_mvn(char *transp, double mu[], double Sig[], int d,
                              size_t len, double X[], int ldx, double L[],
                              randompack_rng *rng);
HIDDEN bool randompack_PY_seed(int seed, uint32_t *spawn_key, int n_key,
                               randompack_rng *rng);
HIDDEN bool randompack_PY_u01(double x[], size_t len, randompack_rng *rng);

#define randompack_create randompack_PY_create
#define randompack_free randompack_PY_free
#define randompack_last_error randompack_PY_last_error
#define randompack_mvn randompack_PY_mvn
#define randompack_seed randompack_PY_seed
#define randompack_u01 randompack_PY_u01

#endif // RANDOMPACKPYTHONGATEWAY_H
