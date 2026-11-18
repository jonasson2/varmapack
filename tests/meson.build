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

test_sources = files( 'RunTests.c', 'FromMatlab.c', 'TestOU.c', 'TestSlepian.c',
  'TestAcvf.c', 'TestAutocov.c', 'TestAutocovEdgeCases.c', 'TestSimEdgeCases.c',
  'TestCov2corr.c', 'TestCov2sig.c', 'TestSimxEdgeCases.c', 'TestTestcase.c',
  'TestVarSpecRad.c',
  'TestFromMatlab.c', 'TestFindCFindG.c', 'TestPsdCondCov.c', 'TestPsi.c',
  'TestTestcasex.c', 'TestLyapunov.c', 'TestMatrixExp.c', 'TestDistribution.c', )

test_util = files('TestUtil.c')

strict_fortran_checks = []
if fc.get_id() == 'gcc'
  strict_fortran_checks += custom_target('mb05od-strict-fortran',
    input : '../src/mb05od-complete.F',
    output : 'mb05od-strict.o',
    command : fc.cmd_array() + ['-c', '-std=f2008', '-pedantic-errors', '-Werror',
                                '@INPUT@', '-o', '@OUTPUT@'], )
endif

test_assertions = executable('TestAssertions', files('TestAssertions.c') + test_util,
  dependencies : varmapack_dep)
test('TestAssertions', test_assertions)

run_tests = executable( 'RunTests',
  test_sources + test_util + files( '../examples/getopt.c',
    '../src/printX.c', ) + src_files + sb03md_fortran_files,
  include_directories : [inc, include_directories('../examples')],
  dependencies : [math_dep, blas_dep, randompack_dep],
  c_args : debug_print_args + local_dpstrf_args + ['-DVARMAPACK_BUILD'], )

test('RunTests', run_tests, args : ['-v'], workdir : meson.current_build_dir(),
  depends : strict_fortran_checks)

against_matlab = executable( 'AgainstMatlab',
  files('RunAgainstMatlab.c', 'TestAgainstMatlab.c', 'FromMatlab.c') + test_util +
  files('../src/printX.c') + src_files + sb03md_fortran_files,
  include_directories : [inc, include_directories('../examples')],
  dependencies : [math_dep, blas_dep, randompack_dep],
  c_args : debug_print_args + local_dpstrf_args + ['-DVARMAPACK_BUILD'], )

test('AgainstMatlab', against_matlab, workdir : meson.current_build_dir())
