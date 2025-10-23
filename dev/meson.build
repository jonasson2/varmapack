project('varmasim', 'c',
  version : '1.0',
  default_options : ['warning_level=2', 'c_std=c11']
)

blas_option = get_option('blas')

# ---- Core sources from ../src ----
src_names = [
  'BlasGateway.c', 'RandomNumbers.c', 'VYWFactorize.c', 'VarmaSim.c',
  'VarmaUtilities.c', 'CCBuild.c', 'SBuild.c', 'VYWSolve.c',
  'SExtend.c', 'FindCGW.c', 'xAssert.c', 'printX.c', 
  'random.c', 'IsStationary.c', 'OmegaBuild.c', 'Testcase.c',
]

src_dir = '../src'
src_files = []
foreach f : src_names
  src_files += files(join_paths(src_dir, f))
endforeach

# ---- Test files and utilities ----
test_files = files('RunTests.c', 'TestRandom.c', 'TestExtraUtil.c', 'TestTestcase.c')
extra_util = ['ExtraUtil.c', 'xCheck.c', 'getopt.c']

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
else
  error('Unknown BLAS option: ' + blas_option + '. Supported: accelerate, openblas, mkl, clapack.')
endif

# ---- Build flags ----
buildtype = get_option('buildtype')
c_args_local = []
if buildtype == 'release'
  c_args_local += ['-funroll-loops', '-fstrict-aliasing', '-march=native', '-flto']
endif

# ---- Test and demo executables ----
run_tests = executable(
  'RunTests',
  test_files + extra_util + src_files,
  include_directories: inc,
  dependencies: [blas_lib],
  c_args: c_args_local
)

run_varmasim = executable(
  'RunVarmaSim',
  ['RunVarmaSim.c'] + extra_util + src_files,
  include_directories: inc,
  dependencies: [blas_lib],
  c_args: c_args_local
)

# ---- Register tests ----
test('RunTests', run_tests)
# (RunVarmaSim is a demo, so we don’t register it as a Meson test.)
