example_util = files(
  'getopt.c',
)

configure_file(
  input : 'prufa_data.txt',
  output : 'prufa_data.txt',
  copy : true,
)

run_varmasim = executable(
  'RunVarmaSim',
  ['RunVarmaSim.c'] + example_util,
  dependencies : [varmapack_dep],
)

run_cryptorandom = executable(
  'RunCryptoRandom',
  ['RunCryptoRandom.c'] + example_util,
  dependencies : [varmapack_dep],
)

prufa = executable(
  'prufa',
  ['prufa.c', '../tests/FromMatlab.c'] + example_util,
  dependencies : [varmapack_dep],
)
alias_target('prufa', prufa)

debugprint_demo = executable(
  'debugprint-demo',
  ['debugprint-demo.c'] + example_util,
  dependencies : [varmapack_dep],
)

xassert_demo = executable(
  'error-demo',
  ['error_demo.c'],
  dependencies : [varmapack_dep],
)

# No `test()` calls here; these are just demos.
