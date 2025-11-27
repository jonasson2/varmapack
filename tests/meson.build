fs = import('fs')

matlab_compare = '../matlab/matlabcompare.txt'
if fs.exists(matlab_compare)
  configure_file(
    input : matlab_compare,
    output : 'matlabcompare.txt',
    copy : true,
  )
else
  message('matlabcompare.txt not found; skipping copy')
endif

test_sources = files(
  'RunTests.c',
  'FromMatlab.c',
  'TestRandom.c',
  'TestRandomMvn.c',
  'TestError.c',
  'TestExtraUtil.c',
  'TestCovar.c',
  'TestTestcase.c',
  'TestVarSpecRad.c',
  'TestAgainstMatlab.c',
  'TestFromMatlab.c',
  'TestFindCG.c',
  'TestPsi.c',
)

extra_util = files(
  'ExtraUtil.c',
  'xCheck.c',
)

run_tests = executable(
  'RunTests',
  test_sources + extra_util,
  dependencies : [varmapack_dep],
)

test('RunTests', run_tests, args : ['-v'])
