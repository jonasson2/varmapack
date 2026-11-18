#pragma once
#include "varmapack.h"

// Source-tree bindings may borrow these arrays for a sim call. Only a model
// returned by varmapack_slepian_create owns its arrays and may be freed.
struct varmapack_slepian_model {
  int r;
  double *theta, *mu, *factor;
};
