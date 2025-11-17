example_util = files(
  'getopt.c',
)

run_varmasim = executable(
  'RunVarmaSim',
  ['RunVarmaSim.c'] + example_util,
  dependencies : [varmapack_dep],
)

prufa = executable(
  'prufa',
  ['prufa.c'] + example_util,
  dependencies : [varmapack_dep],
)

debugprint_demo = executable(
  'debugprint-demo',
  ['debugprint-demo.c'] + example_util,
  dependencies : [varmapack_dep],
)

# No `test()` calls here; these are just demos.
