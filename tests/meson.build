fs = import('fs')

matlab_compare = 'matlabcompare.txt'
if fs.exists(matlab_compare)
  configure_file( input : matlab_compare, output : 'matlabcompare.txt', copy : true, )
else
  message('tests/matlabcompare.txt not found; skipping copy')
endif

matlab_compare_rolling = 'matlabcompare-rolling.txt'
if fs.exists(matlab_compare_rolling)
  configure_file( input : matlab_compare_rolling,
    output : 'matlabcompare-rolling.txt', copy : true, )
else
  message('tests/matlabcompare-rolling.txt not found; skipping copy')
endif

test_sources = files( 'RunTests.c', 'FromMatlab.c', 'TestExtraUtil.c',
  'TestAcvf.c', 'TestAutocov.c', 'TestAutocovEdgeCases.c', 'TestSimEdgeCases.c',
  'TestTestcase.c', 'TestVarSpecRad.c', 'TestAgainstMatlab.c', 'TestFromMatlab.c',
  'TestFindCG.c', 'TestPsdCondCov.c', 'TestPsi.c', )

extra_util = files( 'ExtraUtil.c', 'xCheck.c', )

run_tests = executable( 'RunTests',
  test_sources + extra_util + files( '../examples/getopt.c',
    '../src/FindW.c', '../src/printX.c', '../src/VarmaPackUtil.c', '../src/VYW.c', ),
  include_directories : include_directories('../examples'),
  dependencies : [varmapack_dep], c_args : debug_print_args, )

test('RunTests', run_tests, args : ['-v'], workdir : meson.current_build_dir())

test_lyapunov = executable( 'TestLyapunov',
  files('TestLyapunov.c') + extra_util + src_files + sb03md_fortran_files,
  include_directories : inc, dependencies : [blas_dep, randompack_dep],
  c_args : local_dpstrf_args + ['-DVARMAPACK_BUILD'], )

test('TestLyapunov', test_lyapunov)
