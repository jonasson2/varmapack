configure_file(input: '__init__.py', output: '__init__.py', copy: true)

wrap_sources = [custom_target( 'cython_varmapack_core',
  input: '_core.pyx', output: 'varmapack__core.c',
  command: [cython, '-3', '@INPUT@', '-o', '@OUTPUT@'], )]

py.extension_module( '_core', sources: wrap_sources, include_directories: [inc, py_inc],
  dependencies: [py_dep, varmapack_dep], install: true, )
