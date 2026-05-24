// varmapack_ma_specrad  Spectral radius of companion matrix of a VMA process
#include "VarmaPackUtil.h"

double varmapack_ma_specrad(double B[], int r, int q) {
  return CompanionSpecrad(B, r, q, -1);
}
