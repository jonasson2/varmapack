src_files = files(
  'RandomNumbers.c',
  'VarmaPackUtil.c',
  'error.c',
  'varmapack_VYW.c',
  'varmapack_sim.c',
  'varmapack_covar.c',
  'varmapack_specrad.c',
  'varmapack_testcase.c',
  'lapack_dpstf2.f',
  'lapack_dpstrf.f',
  'printX.c',
  'DebugUtil.c'
)

libvarmapack = library(
  'varmapack',
  src_files,
  include_directories : inc,
  dependencies : [blas_dep],
)

varmapack_dep = declare_dependency(
  link_with : libvarmapack,
  include_directories : inc,
  dependencies : [blas_dep],
)
