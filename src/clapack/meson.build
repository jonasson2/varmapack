clapack_sources = files(
  'daxpy.c',  'dcopy.c',  'ddot.c',   'dgemm.c',  'dgemv.c',  'dger.c',
  'dgetf2.c', 'dgetrf.c', 'dgetrs.c', 'dlamch.c', 'dlaswp.c', 'dnrm2.c',
  'dpotf2.c', 'dpotrf.c', 'dscal.c',  'dsymm.c',  'dsymv.c',  'dsyr.c',
  'dsyr2k.c', 'dsyrk.c',  'dswap.c',  'dtrmm.c',  'dtrmv.c',  'dtrsm.c',
  'dtrsv.c',  'ieeeck.c', 'ilaenv.c', 'iparmq.c', 'lsame.c',  'pow_di.c',
  's_cmp.c',  's_copy.c', 'xerbla.c', 'i_nint.c', 'idamax.c'
)

c_args = [
  '-Wno-logical-op-parentheses',
  '-Wno-unused-variable',
  '-Wno-unused-parameter',
  '-Wno-unused-but-set-variable',
  '-Wno-sometimes-uninitialized'
]

lib_clapack = static_library(
  'clapack',
  sources:             clapack_sources,
  c_args:              c_args
)

clapack_dep = declare_dependency(
  link_with: lib_clapack,
)
