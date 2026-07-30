configure_file(input: '__init__.py', output: '__init__.py', copy: true)

wrap_sources = [custom_target( 'cython_varmapack_core',
  input: '_core.pyx', output: 'varmapack__core.c',
  command: [cython, '-3', '@INPUT@', '-o', '@OUTPUT@'], )]

py.extension_module(
  '_core',
  sources: wrap_sources + core_sources + python_gateway_sources,
  include_directories: inc,
  dependencies: [py_dep, numpy_dep, randompack_api_dep, blas_dep],
  c_args: local_dpstrf_args + scipy_openblas_args + ['-DUSING_PYTHON'],
  fortran_args: scipy_openblas_fortran_args,
  install: true,
  subdir: 'varmapack',
)
