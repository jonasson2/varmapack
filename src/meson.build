src_files = files(
  'RandomNumbers.c',
  'VarmaPackUtil.c',
  'VYWFactorize.c',
  'VYWSolve.c',
  'varmapack_sim.c',
  'varmapack_acvf.c',
  'varmapack_spec_rad.c',
  'varmapack_testcase.c',
  'lapack_dpstf2.f',
  'lapack_dpstrf.f',
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
