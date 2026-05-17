py = import('python').find_installation(pure: false)
py_dep = py.dependency()

numpy_inc = run_command(py, ['-c', 'import numpy; print(numpy.get_include())'],
  check: true).stdout().strip()
py_inc = include_directories(numpy_inc)

cython = find_program('cython', required: true)
subdir('varmapack')
py.install_sources('varmapack/__init__.py', subdir: 'varmapack')

test('python-sim-smoke', py, args: [files('tests/test_sim_smoke.py')],
  env: {'PYTHONPATH': meson.current_build_dir()}, )
