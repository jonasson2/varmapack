src_files = files( 'FindW.c', 'VarmaPackUtil.c', 'error.c', 'exog.c', 'sb03mdGateway.c',
  'Lyapunov.c', 'VYW.c', 'varmapack_sim.c', 'varmapack_acvf.c',
  'varmapack_autocov.c', 'varmapack_psi.c', 'varmapack_specrad.c', 'varmapack_testcase.c',
  'rp_dpstrf.c', )

sb03md_fortran_files = files( 'sb03md-complete.f', )

lib_kwargs = {
  'include_directories' : inc, 'dependencies' : [blas_dep, randompack_dep],
  'c_args' : local_dpstrf_args, }
if host_system == 'windows'
  lib_kwargs += {'vs_module_defs' : 'varmapack.def'}
endif
libvarmapack = library('varmapack', src_files + sb03md_fortran_files, kwargs : lib_kwargs)

varmapack_dep = declare_dependency( link_with : libvarmapack,
  include_directories : inc, dependencies : [blas_dep, randompack_dep], )
