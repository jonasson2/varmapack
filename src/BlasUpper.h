// Include file to change blas names to uppercase if required
#ifndef BLASUPPER_H
#define BLASUPPER_H

#ifdef UPPER
#define daxpy_ DAXPY
#define dcopy_ DCOPY
#define ddot_ DDOT
#define dgemm_ DGEMM
#define dgemv_ DGEMV
#define dger_ DGER
#define dgetrf_ DGETRF
#define dgetrs_ DGETRS
#define dnrm2_ DNRM2
#define idamax_ IDAMAX
#define dpotrf_ DPOTRF
#define dscal_ DSCAL
#define dspr_ DSPR
#define dspr2_ DSPR2
#define dsymm_ DSYMM
#define dsymv_ DSYMV
#define dsyr_ DSYR
#define dsyr2k_ DSYR2K
#define dsyrk_ DSYRK
#define dtrmm_ DTRMM
#define dtrmv_ DTRMV
#define dtrsm_ DTRSM
#define dtrsv_ DTRSV
#endif

#endif
