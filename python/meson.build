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
inc = include_directories('src')
numpy_dep = declare_dependency(compile_args: ['-I' + numpy_inc])

cython = find_program('cython', required: true)
randompack_dep = dependency('randompack', required : true)

blas_order_option = get_option('blas_order')
if blas_order_option == 'auto'
  if host_system == 'windows'
    blas_search_order = ['mkl', 'openblas']
  elif host_system == 'darwin'
    blas_search_order = ['accelerate', 'openblas', 'mkl']
  else
    blas_search_order = ['openblas', 'mkl']
  endif
else
  blas_search_order = blas_order_option.split(',')
endif

blas_dep = disabler()
blas_found = 'none'
foreach blas_name : blas_search_order
  if blas_name == 'mkl'
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
    if not blas_dep.found()
      blas_dep = dependency('scipy-openblas', method: 'pkg-config', required: false)
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
  'src/varmapack_psi.c',
  'src/varmapack_specrad.c',
  'src/varmapack_ma_specrad.c',
  'src/varmapack_testcase.c',
  'src/rp_dpstrf.c',
  'src/sb03md-complete.f',
)

subdir('varmapack')
py.install_sources('varmapack/__init__.py', subdir: 'varmapack')

foreach t : [
  'test_basic.py',
  'test_errors.py',
  'test_moments.py',
  'test_varmax.py',
]
  test('python-' + t.split('.')[0], py, args: [files('tests' / t)],
    env: {'PYTHONPATH': meson.current_build_dir()}, )
endforeach
