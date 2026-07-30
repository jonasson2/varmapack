#include <Python.h>
#include "RandompackPythonGateway.h"

static randompack_python_api *api = 0;

int randompack_PY_initialize(void) {
  uint32_t found_version;
  if (api) return 0;
  api = PyCapsule_Import(RANDOMPACK_PYTHON_API_CAPSULE, 0);
  if (!api) {
    PyErr_Clear();
    PyErr_Format(PyExc_ImportError,
                 "Varmapack requires Randompack Python C API version %d",
                 RANDOMPACK_PYTHON_API_VERSION);
    return -1;
  }
  if (api->abi_version != RANDOMPACK_PYTHON_API_VERSION) {
    found_version = api->abi_version;
    api = 0;
    PyErr_Format(PyExc_ImportError,
                 "Randompack Python C API version %u is incompatible; expected %d",
                 found_version, RANDOMPACK_PYTHON_API_VERSION);
    return -1;
  }
  if (api->struct_size < sizeof(randompack_python_api)) {
    api = 0;
    PyErr_SetString(PyExc_ImportError,
                    "Randompack Python C API function table is incomplete");
    return -1;
  }
  return 0;
}

randompack_rng *randompack_PY_create(const char *engine) {
  return api->create(engine);
}

void randompack_PY_free(randompack_rng *rng) {
  api->free_rng(rng);
}

char *randompack_PY_last_error(randompack_rng *rng) {
  return api->last_error(rng);
}

bool randompack_PY_mvn(char *transp, double mu[], double Sig[], int d,
                       size_t len, double X[], int ldx, double L[],
                       randompack_rng *rng) {
  return api->mvn(transp, mu, Sig, d, len, X, ldx, L, rng);
}

bool randompack_PY_seed(int seed, uint32_t *spawn_key, int n_key,
                        randompack_rng *rng) {
  return api->seed(seed, spawn_key, n_key, rng);
}

bool randompack_PY_u01(double x[], size_t len, randompack_rng *rng) {
  return api->u01(x, len, rng);
}
