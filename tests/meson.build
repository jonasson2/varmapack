test_sources = files(
  'RunTests.c',
  'TestRandom.c',
  'TestExtraUtil.c',
  'TestTestcase.c',
  'TestVarSpecRad.c',
  'TestAgainstMatlab.c',
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
