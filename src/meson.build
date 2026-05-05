src_files = files(
  'VarmaPackUtil.c',
  'error.c',
  'varmapack_VYW.c',
  'varmapack_sim.c',
  'varmapack_acvf.c',
  'varmapack_autocov.c',
  'varmapack_specrad.c',
  'varmapack_testcase.c',
  'rp_dpstrf.c',
  'printX.c',
  'DebugUtil.c',
)

libvarmapack = library(
  'varmapack',
  src_files,
  include_directories : inc,
  dependencies : [blas_dep, randompack_dep],
  c_args : local_dpstrf_args,
)

varmapack_dep = declare_dependency(
  link_with : libvarmapack,
  include_directories : inc,
  dependencies : [blas_dep, randompack_dep],
)
