time_util = files('TimeUtil.c')

time_scalability = executable( 'TimeScalability',
  ['TimeScalability.c', '../examples/getopt.c'] + time_util,
  include_directories : [inc, include_directories('../examples')],
  dependencies : [varmapack_dep], )

time_simulate = executable( 'TimeSimulate',
  ['TimeSimulate.c', '../examples/getopt.c'] + time_util,
  include_directories : [inc, include_directories('../examples')],
  dependencies : [varmapack_dep], )

time_setup = executable( 'TimeSetup',
  ['TimeSetup.c', '../examples/getopt.c'] + time_util + src_files +
    sb03md_fortran_files,
  include_directories : [inc, include_directories('../examples')],
  dependencies : [blas_dep, randompack_dep],
  c_args : local_dpstrf_args + ['-DVARMAPACK_BUILD'], )

time_break_even = executable( 'TimeBreakEven',
  ['TimeBreakEven.c', '../examples/getopt.c'] + time_util + src_files +
    sb03md_fortran_files,
  include_directories : [inc, include_directories('../examples')],
  dependencies : [blas_dep, randompack_dep],
  c_args : local_dpstrf_args + ['-DVARMAPACK_BUILD'], )

# No `test()` call here; this is a benchmark.
