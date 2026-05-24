configure_file(input: '__init__.py', output: '__init__.py', copy: true)

wrap_sources = [custom_target( 'cython_varmapack_core',
  input: '_core.pyx', output: 'varmapack__core.c',
  command: [cython, '-3', '@INPUT@', '-o', '@OUTPUT@'], )]

py.extension_module( '_core', sources: wrap_sources + core_sources, include_directories: inc,
  dependencies: [py_dep, blas_dep, randompack_dep], c_args: local_dpstrf_args,
  install: true, subdir: 'varmapack', )
