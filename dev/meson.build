project('varmasim', 'c',
  version : '1.0',
  default_options : ['warning_level=2']
)

blas_option = get_option('blas')
# Source files for the main executable
src_names = [
  'BlasGateway.c',    'RandomNumbers.c', 'VYWFactorize.c', 'VarmaSim.c',
  'VarmaUtilities.c', 'CCBuild.c',       'SBuild.c',       'VYWSolve.c',
  'SExtend.c',        'FindCGW.c',       'xAssert.c',      'printX.c',
  'random.c',         'IsStationary.c',  'OmegaBuild.c',   'Testcase.c',
]
src_dir = '../src'
src_names = [
  'BlasGateway.c', 'RandomNumbers.c', 'VYWFactorize.c', 'VarmaSim.c',
  'VarmaUtilities.c', 'CCBuild.c', 'SBuild.c', 'VYWSolve.c',
  'SExtend.c', 'FindCGW.c', 'xAssert.c', 'printX.c',
  'random.c', 'IsStationary.c', 'OmegaBuild.c', 'Testcase.c',
]
src_files = []
foreach f : src_names
  src_files += files(join_paths(src_dir, f))
endforeach

# Add files specific to the test executable
test_src_files = [
  'testutilities.c',   # Testing utility functions
  'TestVarmaSim.c',    # Test VarmaSim
  'RunVarmaSim.c', # Main program that calls TestVarmaSim
  'xCheck.c',
  'SCbuild.c',
]

inc_dir = '../src'

# Define executable names
#executable_name = 'VarmaSim'
test_executable_name = 'RunVarmaSim'

# Get the current buildtype
buildtype = get_option('buildtype')

# Base flags
cpp_args = []

# Add additional flags for release builds
if buildtype == 'release'
  cpp_args += ['-funroll-loops', '-fstrict-aliasing', '-march=native', '-flto']
endif

if blas_option == 'accelerate'
  blas_lib = declare_dependency(link_args: ['-framework', 'Accelerate'])
elif blas_option == 'openblas'
  blas_lib = declare_dependency(
    include_directories: include_directories('/opt/homebrew/opt/openblas/include'),
    link_args: ['-L/opt/homebrew/opt/openblas/lib', '-lopenblas']
  )
elif blas_option == 'mkl'
  blas_lib = declare_dependency(
    link_args: ['-lmkl_rt', '-L/opt/intel/mkl/lib', '-I/opt/intel/mkl/include']
  )
elif blas_option == 'clapack'
  subdir('clapack')  # Process clapack subdirectory
  blas_lib = clapack_dep  # Assign clapack_dep to blas_lib
else
  error('Unknown BLAS option: ' + blas_option + '. ' +
        'Supported options: accelerate, openblas, mkl, clapack.')
endif

# # Build the main executable
# executable(
#   executable_name,
#   src_files,
#   include_directories: include_directories('.'),
#   dependencies: [blas_lib],  # Add the BLAS dependency
#   cpp_args: cpp_args,
#   install: true
# )

test_executable = executable(
  test_executable_name,
  test_src_files + src_files, # Includes shared and test-specific files
  include_directories: include_directories('.', '../src'),
  dependencies: [blas_lib],  # Add the BLAS dependency
  cpp_args: cpp_args
)

test('RunVarmaSim', test_executable)
