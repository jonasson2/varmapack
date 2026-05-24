#include "varmapack.h"

const char *varmapack_strerror(varmapack_error error) {
  switch (error) {
    case VARMAPACK_OK:
      return "ok";
    case VARMAPACK_INVALID_ARGUMENT:
      return "invalid argument";
    case VARMAPACK_DIMENSION:
      return "invalid dimension(s)";
    case VARMAPACK_ALLOCATION:
      return "allocation failed";
    case VARMAPACK_UNKNOWN_TESTCASE:
      return "unknown testcase";
    case VARMAPACK_NONSTATIONARY:
      return "nonstationary model";
    case VARMAPACK_NONSTATIONARY_MA:
      return "nonstationary model with MA term(s) and specified X0";
    case VARMAPACK_SINGULAR:
      return "singular matrix";
    case VARMAPACK_NOT_POSITIVE_SEMIDEFINITE:
      return "matrix is not positive semidefinite";
    case VARMAPACK_INTERNAL:
      return "internal error";
    default:
      return "unknown varmapack error";
  }
}
