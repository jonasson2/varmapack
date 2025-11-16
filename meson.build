project('varmapack', 'c', 'fortran',
  version : '1.0',
  default_options : [
    'c_std=c11',
    'warning_level=3',
    'buildtype=debug',
  ]
)

blas_option = get_option('blas')
add_project_arguments('-Wno-do-subscript', language : 'fortran')

# ---- Core sources from ../src ----
src_names = [
  'BlasGateway.c', 'RandomNumbers.c', 'VYWFactorize.c', 'VarmaSim.c',
  'VarmaUtilities.c', 'CCBuild.c', 'SBuild.c', 'VYWSolve.c', 'SExtend.c',
  'FindCGW.c', 'xAssert.c', 'printX.c', 'random.c', 'IsStationary.c',
  'OmegaBuild.c', 'Testcase.c', 'VarSpecRad.c', 'VarmaPackUtil.c',
  'lapack_dpstrf.f', 'lapack_dpstf2.f'
]

src_dir = '../src'
src_files = []
foreach f : src_names
  src_files += files(join_paths(src_dir, f))
endforeach

# ---- Test files and utilities ----
test_files = files('RunTests.c', 'TestRandom.c', 'TestExtraUtil.c', 'TestTestcase.c',
                   'TestVarSpecRad.c', 'TestAgainstMatlab.c', 'TestFindCG.c',
                   'TestPsi.c')
extra_util = ['ExtraUtil.c', 'xCheck.c', 'getopt.c', 'ACVF.c']

# ---- Include directories ----
inc = include_directories('.', '../src')

# ---- BLAS selection ----
if blas_option == 'accelerate'
  blas_lib = declare_dependency(link_args: ['-framework', 'Accelerate'])
elif blas_option == 'openblas'
  openblas_dep = dependency('openblas', required : false)
  if openblas_dep.found()
    blas_lib = openblas_dep
  else
    blas_lib = declare_dependency(
      include_directories: include_directories('/opt/homebrew/opt/openblas/include'),
      link_args: [
        '-L/opt/homebrew/opt/openblas/lib',
        '-Wl,-rpath,/opt/homebrew/opt/openblas/lib',
        '-lopenblas'
      ]
    )
  endif
elif blas_option == 'mkl'
  blas_lib = declare_dependency(
    link_args: [
      '-L/opt/intel/mkl/lib', '-lmkl_rt',
      '-Wl,-rpath,/opt/intel/mkl/lib'
    ],
    include_directories: include_directories('/opt/intel/mkl/include')
  )
elif blas_option == 'clapack'
  subdir('clapack')
  blas_lib = clapack_dep
elif blas_option == 'none'
  blas_lib = declare_dependency()  # no BLAS, for symbol inspection
else
  error('Unknown BLAS option: ' + blas_option + '. Supported: accelerate, openblas, mkl, clapack.')
endif

# ---- Test and demo executables ----
run_tests = executable(
  'RunTests',
  test_files + extra_util + src_files,
  include_directories: inc,
  dependencies: [blas_lib],
)

run_varmasim = executable(
  'RunVarmaSim',
  ['RunVarmaSim.c'] + extra_util + src_files,
  include_directories: inc,
  dependencies: [blas_lib],
)

prufa = executable(
  'prufa',
  ['prufa.c'] + extra_util + src_files,
  include_directories: inc,
  dependencies: [blas_lib],
)

# ---- Register tests ----
test('RunTests', run_tests)
# (RunVarmaSim is a demo, so we don’t register it as a Meson test.)
