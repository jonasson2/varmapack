example_util = files( 'getopt.c', )

run_varmasim = executable( 'RunVarmaSim',
  ['RunVarmaSim.c', '../src/printX.c'] + example_util, dependencies : [varmapack_dep],
  c_args : debug_print_args, )

debugprint_demo = executable( 'debugprint-demo',
  ['debugprint-demo.c'] + example_util, dependencies : [varmapack_dep], )

# No `test()` calls here; these are just demos.
