project('varmapack', 'c', 'fortran',
  version : '1.0',
  default_options : [
    'c_std=c11',
    'warning_level=3',
    'buildtype=debug',
  ]
)

add_project_arguments('-Wno-do-subscript', language : 'fortran')

blas_option = get_option('blas')

# ---- BLAS selection ----
if blas_option == 'accelerate'
  blas_dep = declare_dependency(link_args: ['-framework', 'Accelerate'])

elif blas_option == 'openblas'
  openblas_dep = dependency('openblas', required : false)
  if openblas_dep.found()
    blas_dep = openblas_dep
  else
    blas_dep = declare_dependency(
      include_directories: include_directories('/opt/homebrew/opt/openblas/include'),
      link_args: [
        '-L/opt/homebrew/opt/openblas/lib',
        '-Wl,-rpath,/opt/homebrew/opt/openblas/lib',
        '-lopenblas'
      ]
    )
  endif

elif blas_option == 'mkl'
  blas_dep = declare_dependency(
    link_args: [
      '-L/opt/intel/mkl/lib', '-lmkl_rt',
      '-Wl,-rpath,/opt/intel/mkl/lib'
    ],
    include_directories: include_directories('/opt/intel/mkl/include')
  )

elif blas_option == 'clapack'
  subdir('clapack')
  blas_dep = clapack_dep

elif blas_option == 'none'
  blas_dep = declare_dependency()  # no BLAS, for symbol inspection

else
  error('Unknown BLAS option: ' + blas_option +
        '. Supported: accelerate, openblas, mkl, clapack, none.')
endif

# Your own headers live in src/ for now
inc = include_directories('src')

subdir('src')
subdir('tests')
subdir('examples')
