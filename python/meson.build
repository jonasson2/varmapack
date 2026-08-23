project('varmapack', 'c', 'fortran',
  version: '0.1.0',
  default_options: ['c_std=c11']
)

host_system = host_machine.system()
if host_system != 'windows'
  add_project_arguments('-D_POSIX_C_SOURCE=200809L', language: 'c')
  if host_system == 'darwin'
    add_project_arguments('-D_DARWIN_C_SOURCE', language: 'c')
  endif
endif

py = import('python').find_installation(pure: false)
cc = meson.get_compiler('c')
py_dep = py.dependency()

numpy_inc = run_command(py, ['-c', 'import numpy; print(numpy.get_include())'],
  check: true).stdout().strip()
randompack_include_command = [
  '-c',
  'from pathlib import Path; import randompack; ' +
  'print(Path(randompack.__file__).parent / "include")',
]
randompack_inc = run_command(
  py,
  randompack_include_command,
  check: true,
).stdout().strip()
inc = include_directories('src', 'varmapack')
numpy_dep = declare_dependency(compile_args: ['-I' + numpy_inc])
randompack_api_dep = declare_dependency(compile_args: ['-I' + randompack_inc])

cython = find_program('cython', required: true)

blas_order_option = get_option('blas_order')
if blas_order_option == 'auto'
  if host_system == 'windows'
    blas_search_order = ['scipy-openblas', 'mkl', 'openblas']
  elif host_system == 'darwin'
    blas_search_order = ['accelerate', 'openblas', 'mkl']
  else
    blas_search_order = ['scipy-openblas', 'openblas', 'mkl']
  endif
else
  blas_search_order = blas_order_option.split(',')
endif

blas_dep = disabler()
blas_found = 'none'
scipy_openblas_args = []
scipy_openblas_fortran_args = []
foreach blas_name : blas_search_order
  if blas_name == 'scipy-openblas'
    scipy_openblas_result = run_command(
      py,
      [
        '-c',
        'from pathlib import Path; import scipy_openblas32 as so; ' +
        'print(Path(so.get_lib_dir()) / so.get_library(True))',
      ],
      check: false,
    )
    if scipy_openblas_result.returncode() == 0
      scipy_openblas_library = scipy_openblas_result.stdout().strip()
      if scipy_openblas_library != ''
        blas_dep = declare_dependency(link_args: [scipy_openblas_library])
        blas_found = 'SciPy OpenBLAS'
        scipy_openblas_args = ['-DSCIPY_OPENBLAS']
        scipy_openblas_fortran_args = [
          '-DSCIPY_OPENBLAS',
          '-ffixed-line-length-none',
        ]
      endif
    endif
  elif blas_name == 'mkl'
    foreach mkl_variant : ['mkl-dynamic-lp64-iomp', 'mkl-static-lp64-iomp', 'mkl']
      blas_dep = dependency(mkl_variant, required: false)
      if blas_dep.found()
        break
      endif
    endforeach
    if not blas_dep.found()
      blas_dep = cc.find_library('mkl_rt', required: false)
    endif
    if blas_dep.found()
      blas_found = 'MKL'
    endif
  elif blas_name == 'accelerate'
    if host_system == 'darwin'
      blas_dep = dependency('appleframeworks', modules: 'Accelerate', required: false)
      if blas_dep.found()
        blas_found = 'Accelerate'
      endif
    endif
  elif blas_name == 'openblas'
    blas_dep = dependency('openblas', required: false)
    if not blas_dep.found()
      blas_dep = cc.find_library('openblas', required: false)
    endif
    if blas_dep.found()
      blas_found = 'OpenBLAS'
    endif
  else
    error('Unknown BLAS option: ' + blas_name)
  endif
  if blas_dep.found()
    break
  endif
endforeach
if not blas_dep.found()
  error('No optimized BLAS found, reference implementation not yet available')
endif
message('Using BLAS library: ' + blas_found)

using_accelerate = (host_system == 'darwin') and (blas_found == 'Accelerate')
local_dpstrf_args = using_accelerate ? ['-DLOCAL_DPSTRF'] : []

core_sources = files(
  'src/FindW.c',
  'src/VarmaPackUtil.c',
  'src/error.c',
  'src/exog.c',
  'src/sb03mdGateway.c',
  'src/Lyapunov.c',
  'src/VYW.c',
  'src/varmapack_sim.c',
  'src/varmapack_simx.c',
  'src/varmapack_acvf.c',
  'src/varmapack_autocov.c',
  'src/varmapack_cov2corr.c',
  'src/varmapack_psi.c',
  'src/varmapack_specrad.c',
  'src/varmapack_ma_specrad.c',
  'src/varmapack_testcase.c',
  'src/rp_dpstrf.c',
  'src/sb03md-complete.F',
)
python_gateway_sources = files('varmapack/RandompackPythonGateway.c')

subdir('varmapack')
py.install_sources('varmapack/__init__.py', subdir: 'varmapack')

foreach t : [
  'test_basic.py',
  'test_data.py',
  'test_errors.py',
  'test_interop.py',
  'test_moments.py',
  'test_testcases.py',
  'test_varmax.py',
]
  test('python-' + t.split('.')[0], py,
    args: ['-m', 'pytest', '-q', files('tests' / t)],
    env: {'PYTHONPATH': meson.current_build_dir()}, )
endforeach
