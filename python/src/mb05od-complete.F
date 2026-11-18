C     Bundled SLICOT MB05OD, MB04MD and MB05OY, unchanged upstream.
C     https://github.com/SLICOT/SLICOT-Reference
C     Commit 5eeb845d84f5a60f6069754b7819dad468951933
C     Only this provenance, license and BLAS symbol preamble are added.
C
C     Copyright (c) 1996-2026, The SLICOT Team
C     All rights reserved.
C
C     Redistribution and use in source and binary forms, with or without
C     modification, are permitted provided that the following conditions are met:
C
C     1. Redistributions of source code must retain the above copyright notice, this
C        list of conditions and the following disclaimer.
C
C     2. Redistributions in binary form must reproduce the above copyright notice,
C        this list of conditions and the following disclaimer in the documentation
C        and/or other materials provided with the distribution.
C
C     3. Neither the name of the copyright holder nor the names of its
C        contributors may be used to endorse or promote products derived from
C        this software without specific prior written permission.
C
C     THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
C     AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
C     IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
C     DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
C     FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
C     DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
C     SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
C     CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
C     OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C     OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#if defined(SCIPY_OPENBLAS)
#define DASUM SCIPY_DASUM
#define DCOPY SCIPY_DCOPY
#define DGEMM SCIPY_DGEMM
#define DGEMV SCIPY_DGEMV
#define DGETRF SCIPY_DGETRF
#define DGETRS SCIPY_DGETRS
#define DLACPY SCIPY_DLACPY
#define DLAMCH SCIPY_DLAMCH
#define DLANGE SCIPY_DLANGE
#define DLASCL SCIPY_DLASCL
#define DLASET SCIPY_DLASET
#define DNRM2 SCIPY_DNRM2
#define DSCAL SCIPY_DSCAL
#define DSWAP SCIPY_DSWAP
#define IDAMAX SCIPY_IDAMAX
#define LSAME SCIPY_LSAME
#define XERBLA SCIPY_XERBLA
#endif
C
C SPDX-License-Identifier: BSD-3-Clause
C
      SUBROUTINE MB05OD( BALANC, N, NDIAG, DELTA, A, LDA, MDIG, IDIG,
     $                   IWORK, DWORK, LDWORK, IWARN, INFO )
C
C     PURPOSE
C
C     To compute exp(A*delta) where A is a real N-by-N matrix and delta
C     is a scalar value. The routine also returns the minimal number of
C     accurate digits in the 1-norm of exp(A*delta) and the number of
C     accurate digits in the 1-norm of exp(A*delta) at 95% confidence
C     level.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     BALANC  (input) CHARACTER*1
C             Specifies whether or not a balancing transformation (done
C             by SLICOT Library routine MB04MD) is required, as follows:
C             = 'N', do not use balancing;
C             = 'S', use balancing (scaling).
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix A.  N >= 0.
C
C     NDIAG   (input) INTEGER
C             The specified order of the diagonal Pade approximant.
C             In the absence of further information NDIAG should
C             be set to 9.  NDIAG should not exceed 15.  NDIAG >= 1.
C
C     DELTA   (input) DOUBLE PRECISION
C             The scalar value delta of the problem.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On input, the leading N-by-N part of this array must
C             contain the matrix A of the problem. (This is not needed
C             if DELTA = 0.)
C             On exit, if INFO = 0, the leading N-by-N part of this
C             array contains the solution matrix exp(A*delta).
C
C     LDA     (input) INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     MDIG    (output) INTEGER
C             The minimal number of accurate digits in the 1-norm of
C             exp(A*delta).
C
C     IDIG    (output) INTEGER
C             The number of accurate digits in the 1-norm of
C             exp(A*delta) at 95% confidence level.
C
C     Workspace
C
C     IWORK   (input/output) INTEGER array, dimension (N)
C
C     DWORK   (input/output) DOUBLE PRECISION array, dimension (LDWORK)
C
C     LDWORK  (input) INTEGER
C             The length of the array DWORK.
C             LDWORK >= N*(2*N+NDIAG+1)+NDIAG, if N >  1.
C             LDWORK >= 1,                     if N <= 1.
C
C     Warning Indicator
C
C     IWARN   (output) INTEGER
C             = 0:  no warning;
C             = 1:  if MDIG = 0 and IDIG > 0, warning for possible
C                   inaccuracy (the exponential has been computed);
C             = 2:  if MDIG = 0 and IDIG = 0, warning for severe
C                   inaccuracy (the exponential has been computed);
C             = 3:  if balancing has been requested, but it failed to
C                   reduce the matrix norm and was not actually used.
C
C     Error Indicator
C
C     INFO    (output) INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  if the norm of matrix A*delta (after a possible
C                   balancing) is too large to obtain an accurate
C                   result;
C             = 2:  if the coefficient matrix (the denominator of the
C                   Pade approximant) is exactly singular; try a
C                   different value of NDIAG;
C             = 3:  if the solution exponential would overflow, possibly
C                   due to a too large value DELTA; the calculations
C                   stopped prematurely. This error is not likely to
C                   appear.
C
C     METHOD
C
C     The exponential of the matrix A is evaluated from a diagonal Pade
C     approximant. This routine is a modification of the subroutine
C     PADE, described in reference [1]. The routine implements an
C     algorithm which exploits the identity
C
C         (exp[(2**-m)*A]) ** (2**m) = exp(A),
C
C     where m is an integer determined by the algorithm, to improve the
C     accuracy for matrices with large norms.
C
C     REFERENCES
C
C     [1] Ward, R.C.
C         Numerical computation of the matrix exponential with accuracy
C         estimate.
C         SIAM J. Numer. Anal., 14, pp. 600-610, 1977.
C
C     NUMERICAL ASPECTS
C                               3
C     The algorithm requires 0(N ) operations.
C
C     CONTRIBUTORS
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Apr. 1997.
C     Supersedes Release 2.0 routine MB05CD by T.W.C. Williams, Kingston
C     Polytechnic, March 1982.
C
C     REVISIONS
C
C     June 14, 1997, April 25, 2003, December 12, 2004.
C
C     KEYWORDS
C
C     Continuous-time system, matrix algebra, matrix exponential,
C     matrix operations, Pade approximation.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, HALF, ONE, TWO, FOUR, EIGHT, TEN, TWELVE,
     $                  NINTEN, TWO4, FOUR7, TWOHND
      PARAMETER         ( ZERO = 0.0D0,    HALF = 0.5D0,   ONE = 1.0D0,
     $                    TWO  = 2.0D0,    FOUR = 4.0D0, EIGHT = 8.0D0,
     $                    TEN  = 10.0D0, TWELVE = 12.0D0,
     $                    NINTEN = 19.0D0, TWO4 = 24.0D0,
     $                    FOUR7  = 47.0D0, TWOHND = 200.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         BALANC
      INTEGER           IDIG, INFO, IWARN, LDA, LDWORK, MDIG, N,
     $                  NDIAG
      DOUBLE PRECISION  DELTA
C     .. Array Arguments ..
      INTEGER           IWORK(*)
      DOUBLE PRECISION  A(LDA,*), DWORK(*)
C     .. Local Scalars ..
      LOGICAL           LBALS
      CHARACTER         ACTBAL
      INTEGER           BASE, I, IFAIL, IJ, IK, IM1, J, JWORA1, JWORA2,
     $                  JWORA3, JWORV1, JWORV2, K, M, MPOWER, NDAGM1,
     $                  NDAGM2, NDEC, NDECM1
      DOUBLE PRECISION  ANORM, AVGEV, BD, BIG, EABS, EAVGEV, EMNORM,
     $                  EPS, FACTOR, FN, GN, MAXRED, OVRTH2, OVRTHR, P,
     $                  RERL, RERR, S, SD2, SIZE, SMALL, SS, SUM2D,
     $                  TEMP, TMP1, TR, U, UNDERF, VAR, VAREPS, XN
C     .. External Functions ..
      LOGICAL           LSAME
      DOUBLE PRECISION  DASUM, DLAMCH, DLANGE, DNRM2
      EXTERNAL          DASUM, DLAMCH, DLANGE, DNRM2, LSAME
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEMM, DGEMV, DGETRF, DGETRS, DLACPY,
     $                  DLASCL, DLASET, DSCAL, MB04MD, MB05OY, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, EXP, INT, LOG, LOG10, MAX, MIN, MOD, SQRT
C     .. Executable Statements ..
C
      IWARN = 0
      INFO  = 0
      LBALS = LSAME( BALANC, 'S' )
C
C     Test the input scalar arguments.
C
      IF( .NOT.( LSAME( BALANC, 'N' ) .OR. LBALS ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NDIAG.LT.1 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDWORK.LT.1 .OR.
     $       ( LDWORK.LT.N*( 2*N + NDIAG + 1 ) + NDIAG .AND. N.GT.1 )
     $       ) THEN
         INFO = -11
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'MB05OD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      EPS  = DLAMCH( 'Epsilon' )
      NDEC = INT( LOG10( ONE/EPS ) + ONE )
C
      IF ( N.EQ.0 ) THEN
         MDIG = NDEC
         IDIG = NDEC
         RETURN
      END IF
C
C     Set some machine parameters.
C
      BASE = DLAMCH( 'Base' )
      NDECM1 = NDEC - 1
      UNDERF = DLAMCH( 'Underflow' )
      OVRTHR = DLAMCH( 'Overflow' )
      OVRTH2 = SQRT( OVRTHR )
C
      IF ( DELTA.EQ.ZERO ) THEN
C
C        The DELTA = 0 case.
C
         CALL DLASET( 'Full', N, N, ZERO, ONE, A, LDA )
         MDIG = NDECM1
         IDIG = NDECM1
         RETURN
      END IF
C
      IF ( N.EQ.1 ) THEN
C
C        The 1-by-1 case.
C
         A(1,1) = EXP( A(1,1)*DELTA )
         MDIG = NDECM1
         IDIG = NDECM1
         RETURN
      END IF
C
C     Set pointers for the workspace.
C
      JWORA1 = 1
      JWORA2 = JWORA1 + N*N
      JWORA3 = JWORA2 + N*NDIAG
      JWORV1 = JWORA3 + N*N
      JWORV2 = JWORV1 + N
C
C     Compute Pade coefficients in DWORK(JWORV2:JWORV2+NDIAG-1).
C
      DWORK(JWORV2) = HALF
C
      DO 20 I = 2, NDIAG
         IM1 = I - 1
         DWORK(JWORV2+IM1) = DWORK(JWORV2+I-2)*DBLE( NDIAG - IM1 )/
     $                                   DBLE( I*( 2*NDIAG - IM1 ) )
   20 CONTINUE
C
      VAREPS = EPS**2*( ( DBLE( BASE )**2 - ONE )/
     $        ( TWO4*LOG( DBLE( BASE ) ) ) )
      XN = DBLE( N )
      TR = ZERO
C
C     Apply a translation with the mean of the eigenvalues of A*DELTA.
C
      DO 40 I = 1, N
         CALL DSCAL( N, DELTA, A(1,I), 1 )
         TR = TR + A(I,I)
   40 CONTINUE
C
      AVGEV = TR/XN
      IF ( AVGEV.GT.LOG( OVRTHR ) .OR. AVGEV.LT.LOG( UNDERF ) )
     $   AVGEV = ZERO
      IF ( AVGEV.NE.ZERO ) THEN
         ANORM = DLANGE( '1-norm', N, N, A, LDA, DWORK(JWORA1) )
C
         DO 60 I = 1, N
            A(I,I) = A(I,I) - AVGEV
   60    CONTINUE
C
         TEMP = DLANGE( '1-norm', N, N, A, LDA, DWORK(JWORA1) )
         IF ( TEMP.GT.HALF*ANORM ) THEN
C
            DO 80 I = 1, N
               A(I,I) = A(I,I) + AVGEV
   80       CONTINUE
C
            AVGEV = ZERO
         END IF
      END IF
      ACTBAL = BALANC
      IF ( LBALS ) THEN
C
C        Balancing (scaling) has been requested.  First, save A.
C
         CALL DLACPY( 'Full', N, N, A, LDA, DWORK(JWORA1), N )
         MAXRED = TWOHND
         CALL MB04MD( N, MAXRED, A, LDA, DWORK(JWORV1), INFO )
         IF ( MAXRED.LT.ONE ) THEN
C
C           Recover the matrix and reset DWORK(JWORV1,...,JWORV1+N-1)
C           to 1, as no reduction of the norm occured (unlikely event).
C
            CALL DLACPY( 'Full', N, N, DWORK(JWORA1), N, A, LDA )
            ACTBAL = 'N'
            DWORK(JWORV1) = ONE
            CALL DCOPY( N-1, DWORK(JWORV1), 0, DWORK(JWORV1+1), 1 )
            IWARN = 3
         END IF
      END IF
C
C     Scale the matrix by 2**(-M), where M is the minimum integer
C     so that the resulted matrix has the 1-norm less than 0.5.
C
      ANORM = DLANGE( '1-norm', N, N, A, LDA, DWORK(JWORA1) )
      M = 0
      IF ( ANORM.GE.HALF ) THEN
         MPOWER = INT( LOG( OVRTHR )/LOG( TWO ) )
         M = INT( LOG( ANORM )/LOG( TWO ) ) + 1
         IF ( M.GT.MPOWER ) THEN
C
C           Error return: The norm of A*DELTA is too large.
C
            INFO = 1
            RETURN
         END IF
         FACTOR = TWO**M
         IF ( M+1.LT.MPOWER ) THEN
            M = M + 1
            FACTOR = FACTOR*TWO
         END IF
C
         DO 120 I = 1, N
            CALL DSCAL( N, ONE/FACTOR, A(1,I), 1 )
  120    CONTINUE
C
      END IF
      NDAGM1 = NDIAG - 1
      NDAGM2 = NDAGM1 - 1
      IJ = 0
C
C     Compute the factors of the diagonal Pade approximant.
C     The loop 200 takes the accuracy requirements into account:
C     Pade coefficients decrease with K, so the calculations should
C     be performed in backward order, one column at a time.
C     (A BLAS 3 implementation in forward order, using DGEMM, could
C     possibly be less accurate.)
C
      DO 200 J = 1, N
         CALL DGEMV( 'No transpose', N, N, ONE, A, LDA, A(1,J), 1, ZERO,
     $               DWORK(JWORA2), 1 )
         IK = 0
C
         DO 140 K = 1, NDAGM2
            CALL DGEMV( 'No transpose', N, N, ONE, A, LDA,
     $                  DWORK(JWORA2+IK), 1, ZERO, DWORK(JWORA2+IK+N),
     $                  1 )
            IK = IK + N
  140    CONTINUE
C
         DO 180 I = 1, N
            S = ZERO
            U = ZERO
            IK = NDAGM2*N + I - 1
C
            DO 160 K = NDAGM1, 1, -1
               P = DWORK(JWORV2+K)*DWORK(JWORA2+IK)
               IK = IK - N
               S = S + P
               IF ( MOD( K+1, 2 ).EQ.0 ) THEN
                  U = U + P
               ELSE
                  U = U - P
               END IF
  160       CONTINUE
C
            P = DWORK(JWORV2)*A(I,J)
            S = S + P
            U = U - P
            IF ( I.EQ.J ) THEN
               S = S + ONE
               U = U + ONE
            END IF
            DWORK(JWORA3+IJ) = S
            DWORK(JWORA1+IJ) = U
            IJ = IJ + 1
  180    CONTINUE
C
  200 CONTINUE
C
C     Compute the exponential of the scaled matrix, using diagonal Pade
C     approximants.  As, in theory [1], the denominator of the Pade
C     approximant should be very well conditioned, no condition estimate
C     is computed.
C
      CALL DGETRF( N, N, DWORK(JWORA1), N, IWORK, IFAIL )
      IF ( IFAIL.GT.0 ) THEN
C
C        Error return: The matrix is exactly singular.
C
         INFO = 2
         RETURN
      END IF
C
      CALL DLACPY( 'Full', N, N, DWORK(JWORA3), N, A, LDA )
      CALL DGETRS( 'No transpose', N, N, DWORK(JWORA1), N, IWORK, A,
     $             LDA, IFAIL )
C
C     Prepare for the calculation of the accuracy estimates.
C     Note that ANORM here is in the range [1, e].
C
      ANORM = DLANGE( '1-norm', N, N, A, LDA, DWORK(JWORA1) )
      IF ( ANORM.GE.ONE ) THEN
         EABS = ( NINTEN*XN + FOUR7 )*( EPS*ANORM )
      ELSE
         EABS = ( ( NINTEN*XN + FOUR7 )*EPS )*ANORM
      END IF
      IF ( M.NE.0 ) THEN
         VAR = XN*VAREPS
         FN = ( FOUR*XN )/( ( XN + TWO )*( XN + ONE ) )
         GN = ( ( TWO*XN + TEN )*XN - FOUR )/( ( ( XN + TWO )**2 )
     $                                        *( ( XN + ONE )**2 ) )
C
C        Square-up the computed exponential matrix M times, with caution
C        for avoiding overflows.
C
         DO 220 K = 1, M
            IF ( ANORM.GT.OVRTH2 ) THEN
C
C              The solution could overflow.
C
               CALL DGEMM( 'No transpose', 'No transpose', N, N, N,
     $                     ONE/ANORM, A, LDA, A, LDA, ZERO,
     $                     DWORK(JWORA1), N )
               S = DLANGE( '1-norm', N, N, DWORK(JWORA1), N,
     $                     DWORK(JWORA1) )
               IF ( ANORM.LE.OVRTHR/S ) THEN
                  CALL DLASCL( 'General', N, N, ONE, ANORM, N, N,
     $                         DWORK(JWORA1), N, INFO )
                  TEMP = OVRTHR
               ELSE
C
C                 Error return: The solution would overflow.
C                 This will not happen on most machines, due to the
C                 selection of M.
C
                  INFO = 3
                  RETURN
               END IF
            ELSE
               CALL DGEMM( 'No transpose', 'No transpose', N, N, N, ONE,
     $                     A, LDA, A, LDA, ZERO, DWORK(JWORA1), N )
               TEMP = ANORM**2
            END IF
            IF ( EABS.LT.ONE ) THEN
               EABS = ( TWO*ANORM + EABS )*EABS  + XN*( EPS*TEMP )
            ELSE IF ( EABS.LT.SQRT( ONE - XN*EPS + OVRTHR/TEMP )*ANORM -
     $                        ANORM ) THEN
               EABS = XN*( EPS*TEMP ) + TWO*( ANORM*EABS ) + EABS**2
            ELSE
               EABS = OVRTHR
            END IF
C
            TMP1 = FN*VAR + GN*( TEMP*VAREPS )
            IF ( TMP1.GT.OVRTHR/TEMP ) THEN
               VAR = OVRTHR
            ELSE
               VAR = TMP1*TEMP
            END IF
C
            CALL DLACPY( 'Full', N, N, DWORK(JWORA1), N, A, LDA )
            ANORM = DLANGE( '1-norm', N, N, A, LDA, DWORK(JWORA1) )
  220    CONTINUE
C
      ELSE
         VAR = ( TWELVE*XN )*VAREPS
      END IF
C
C     Apply back transformations, if balancing was effectively used.
C
      CALL MB05OY( ACTBAL, N, 1, N, A, LDA, DWORK(JWORV1), INFO )
      EAVGEV = EXP( AVGEV )
      EMNORM = DLANGE( '1-norm', N, N, A, LDA, DWORK(JWORA1) )
C
C     Compute auxiliary quantities needed for the accuracy estimates.
C
      BIG = ONE
      SMALL = ONE
      IF ( LBALS ) THEN
C
C        Compute norms of the diagonal scaling matrix and its inverse.
C
         DO 240 I = 1, N
            U = DWORK(JWORV1+I-1)
            IF (   BIG.LT.U )   BIG = U
            IF ( SMALL.GT.U ) SMALL = U
  240    CONTINUE
C
         SUM2D = DNRM2( N, DWORK(JWORV1), 1 )
      ELSE
         SUM2D = SQRT( XN )
      END IF
C
C     Update the exponential for the initial translation, and update the
C     auxiliary quantities needed for the accuracy estimates.
C
      SD2 = SQRT( EIGHT*XN*VAREPS )*ANORM
      BD  = SQRT( VAR )
      SS  = MAX( BD, SD2 )
      BD  = MIN( BD, SD2 )
      SD2 = SS*SQRT( ONE + ( BD/SS )**2 )
      IF ( SD2.LE.ONE ) THEN
         SD2 = ( TWO/XN )*SUM2D*SD2
      ELSE IF ( SUM2D/XN.LT.OVRTHR/TWO/SD2 ) THEN
         SD2 = ( TWO/XN )*SUM2D*SD2
      ELSE
         SD2 = OVRTHR
      END IF
      IF ( LBALS ) THEN
         SIZE = ZERO
      ELSE
         IF ( SD2.LT.OVRTHR - EMNORM ) THEN
            SIZE = EMNORM + SD2
         ELSE
            SIZE = OVRTHR
         END IF
      END IF
C
      DO 260 J = 1, N
         SS = DASUM( N, A(1,J), 1 )
         CALL DSCAL( N, EAVGEV, A(1,J), 1 )
         IF ( LBALS ) THEN
            BD = DWORK(JWORV1+J-1)
            SIZE = MAX( SIZE, SS + SD2/BD )
         END IF
  260 CONTINUE
C
C     Set the accuracy estimates and warning errors, if any.
C
      RERR = LOG10( BIG ) + LOG10( EABS ) - LOG10( SMALL ) -
     $       LOG10( EMNORM ) - LOG10( EPS )
      IF ( SIZE.GT.EMNORM ) THEN
         RERL = LOG10( ( SIZE/EMNORM - ONE )/EPS )
      ELSE
         RERL = ZERO
      END IF
      MDIG = MIN( NDEC - INT( RERR + HALF ), NDECM1 )
      IDIG = MIN( NDEC - INT( RERL + HALF ), NDECM1 )
C
      IF ( MDIG.LE.0 ) THEN
         MDIG = 0
         IWARN = 1
      END IF
      IF ( IDIG.LE.0 ) THEN
         IDIG = 0
         IWARN = 2
      END IF
C
      RETURN
C *** Last line of MB05OD ***
      END
C
C SPDX-License-Identifier: BSD-3-Clause
C
      SUBROUTINE MB04MD( N, MAXRED, A, LDA, SCALE, INFO )
C
C     PURPOSE
C
C     To reduce the 1-norm of a general real matrix A by balancing.
C     This involves diagonal similarity transformations applied
C     iteratively to A to make the rows and columns as close in norm as
C     possible.
C
C     This routine can be used instead LAPACK Library routine DGEBAL,
C     when no reduction of the 1-norm of the matrix is possible with
C     DGEBAL, as for upper triangular matrices. LAPACK Library routine
C     DGEBAK, with parameters ILO = 1, IHI = N, and JOB = 'S', should
C     be used to apply the backward transformation.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix A.  N >= 0.
C
C     MAXRED  (input/output) DOUBLE PRECISION
C             On entry, the maximum allowed reduction in the 1-norm of
C             A (in an iteration) if zero rows or columns are
C             encountered.
C             If MAXRED > 0.0, MAXRED must be larger than one (to enable
C             the norm reduction).
C             If MAXRED <= 0.0, then the value 10.0 for MAXRED is
C             used.
C             On exit, if the 1-norm of the given matrix A is non-zero,
C             the ratio between the 1-norm of the given matrix and the
C             1-norm of the balanced matrix. Usually, this ratio will be
C             larger than one, but it can sometimes be one, or even less
C             than one (for instance, for some companion matrices).
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the input matrix A.
C             On exit, the leading N-by-N part of this array contains
C             the balanced matrix.
C
C     LDA     (input) INTEGER
C             The leading dimension of the array A.  LDA >= max(1,N).
C
C     SCALE   (output) DOUBLE PRECISION array, dimension (N)
C             The scaling factors applied to A.  If D(j) is the scaling
C             factor applied to row and column j, then SCALE(j) = D(j),
C             for j = 1,...,N.
C
C     Error Indicator
C
C     INFO    (output) INTEGER
C             = 0:  successful exit.
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value.
C
C     METHOD
C
C     Balancing consists of applying a diagonal similarity
C     transformation inv(D) * A * D to make the 1-norms of each row
C     of A and its corresponding column nearly equal.
C
C     Information about the diagonal matrix D is returned in the vector
C     SCALE.
C
C     REFERENCES
C
C     [1] Anderson, E., Bai, Z., Bischof, C., Demmel, J., Dongarra, J.,
C         Du Croz, J., Greenbaum, A., Hammarling, S., McKenney, A.,
C         Ostrouchov, S., and Sorensen, D.
C         LAPACK Users' Guide: Second Edition.
C         SIAM, Philadelphia, 1995.
C
C     NUMERICAL ASPECTS
C
C     None.
C
C     CONTRIBUTORS
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, June 1997.
C     Supersedes Release 2.0 routine MB04AD by T.W.C. Williams,
C     Kingston Polytechnic, United Kingdom, October 1984.
C     This subroutine is based on LAPACK routine DGEBAL, and routine
C     BALABC (A. Varga, German Aerospace Research Establishment, DLR).
C
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Balancing, eigenvalue, matrix algebra, matrix operations,
C     similarity transformation.
C
C  *********************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   SCLFAC
      PARAMETER          ( SCLFAC = 1.0D+1 )
      DOUBLE PRECISION   FACTOR, MAXR
      PARAMETER          ( FACTOR = 0.95D+0, MAXR = 10.0D+0 )
C     ..
C     .. Scalar Arguments ..
      INTEGER            INFO, LDA, N
      DOUBLE PRECISION   MAXRED
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), SCALE( * )
C     ..
C     .. Local Scalars ..
      LOGICAL            NOCONV
      INTEGER            I, ICA, IRA, J
      DOUBLE PRECISION   ANORM, C, CA, F, G, MAXNRM, R, RA, S, SFMAX1,
     $                   SFMAX2, SFMIN1, SFMIN2, SRED
C     ..
C     .. External Functions ..
      INTEGER            IDAMAX
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           DLAMCH, DLANGE, IDAMAX
C     ..
C     .. External Subroutines ..
      EXTERNAL           DSCAL, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
C     ..
C     .. Executable Statements ..
C
C     Test the scalar input arguments.
C
      INFO  = 0
C
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( MAXRED.GT.ZERO .AND. MAXRED.LT.ONE ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB04MD', -INFO )
         RETURN
      END IF
C
      IF( N.EQ.0 )
     $   RETURN
C
      DO 10 I = 1, N
         SCALE( I ) = ONE
   10 CONTINUE
C
C     Compute the 1-norm of matrix A and exit if it is zero.
C
      ANORM = DLANGE( '1-norm', N, N, A, LDA, SCALE )
      IF( ANORM.EQ.ZERO )
     $   RETURN
C
C     Set some machine parameters and the maximum reduction in the
C     1-norm of A if zero rows or columns are encountered.
C
      SFMIN1 = DLAMCH( 'S' ) / DLAMCH( 'P' )
      SFMAX1 = ONE / SFMIN1
      SFMIN2 = SFMIN1*SCLFAC
      SFMAX2 = ONE / SFMIN2
C
      SRED = MAXRED
      IF( SRED.LE.ZERO ) SRED = MAXR
C
      MAXNRM = MAX( ANORM/SRED, SFMIN1 )
C
C     Balance the matrix.
C
C     Iterative loop for norm reduction.
C
   20 CONTINUE
      NOCONV = .FALSE.
C
      DO 80 I = 1, N
         C = ZERO
         R = ZERO
C
         DO 30 J = 1, N
            IF( J.EQ.I )
     $         GO TO 30
            C = C + ABS( A( J, I ) )
            R = R + ABS( A( I, J ) )
   30    CONTINUE
         ICA = IDAMAX( N, A( 1, I ), 1 )
         CA = ABS( A( ICA, I ) )
         IRA = IDAMAX( N, A( I, 1 ), LDA )
         RA = ABS( A( I, IRA ) )
C
C        Special case of zero C and/or R.
C
         IF( C.EQ.ZERO .AND. R.EQ.ZERO )
     $      GO TO 80
         IF( C.EQ.ZERO ) THEN
            IF( R.LE.MAXNRM)
     $         GO TO 80
            C = MAXNRM
         END IF
         IF( R.EQ.ZERO ) THEN
            IF( C.LE.MAXNRM )
     $         GO TO 80
            R = MAXNRM
         END IF
C
C        Guard against zero C or R due to underflow.
C
         G = R / SCLFAC
         F = ONE
         S = C + R
   40    CONTINUE
         IF( C.GE.G .OR. MAX( F, C, CA ).GE.SFMAX2 .OR.
     $       MIN( R, G, RA ).LE.SFMIN2 )GO TO 50
         F = F*SCLFAC
         C = C*SCLFAC
         CA = CA*SCLFAC
         R = R / SCLFAC
         G = G / SCLFAC
         RA = RA / SCLFAC
         GO TO 40
C
   50    CONTINUE
         G = C / SCLFAC
   60    CONTINUE
         IF( G.LT.R .OR. MAX( R, RA ).GE.SFMAX2 .OR.
     $       MIN( F, C, G, CA ).LE.SFMIN2 )GO TO 70
         F = F / SCLFAC
         C = C / SCLFAC
         G = G / SCLFAC
         CA = CA / SCLFAC
         R = R*SCLFAC
         RA = RA*SCLFAC
         GO TO 60
C
C        Now balance.
C
   70    CONTINUE
         IF( ( C+R ).GE.FACTOR*S )
     $      GO TO 80
         IF( F.LT.ONE .AND. SCALE( I ).LT.ONE ) THEN
            IF( F*SCALE( I ).LE.SFMIN1 )
     $         GO TO 80
         END IF
         IF( F.GT.ONE .AND. SCALE( I ).GT.ONE ) THEN
            IF( SCALE( I ).GE.SFMAX1 / F )
     $         GO TO 80
         END IF
         G = ONE / F
         SCALE( I ) = SCALE( I )*F
         NOCONV = .TRUE.
C
         CALL DSCAL( N, G, A( I, 1 ), LDA )
         CALL DSCAL( N, F, A( 1, I ), 1 )
C
   80 CONTINUE
C
      IF( NOCONV )
     $   GO TO 20
C
C     Set the norm reduction parameter.
C
      MAXRED = ANORM/DLANGE( '1-norm', N, N, A, LDA, SCALE )
C
      RETURN
C *** End of MB04MD ***
      END
C
C SPDX-License-Identifier: BSD-3-Clause
C
      SUBROUTINE MB05OY( JOB, N, LOW, IGH, A, LDA, SCALE, INFO )
C
C     PURPOSE
C
C     To restore a matrix after it has been transformed by applying
C     balancing transformations (permutations and scalings), as
C     determined by LAPACK Library routine DGEBAL.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOB     (input) CHARACTER*1
C             Specifies the type of backward transformation required,
C             as follows:
C             = 'N', do nothing, return immediately;
C             = 'P', do backward transformation for permutation only;
C             = 'S', do backward transformation for scaling only;
C             = 'B', do backward transformations for both permutation
C                    and scaling.
C             JOB must be the same as the argument JOB supplied
C             to DGEBAL.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix A.  N >= 0.
C
C     LOW     (input) INTEGER
C             See IGH
C
C     IGH     (input) INTEGER
C             The integers LOW and IGH determined by DGEBAL.
C             1 <= LOW <= IGH <= N, if N > 0; LOW=1 and IGH=0, if N=0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the matrix to be back-transformed.
C             On exit, the leading N-by-N part of this array contains
C             the transformed matrix.
C
C     LDA     (input) INTEGER
C             The leading dimension of the array A.  LDA >= max(1,N).
C
C     SCALE   (input) DOUBLE PRECISION array, dimension (N)
C             Details of the permutation and scaling factors, as
C             returned by DGEBAL.
C
C     Error Indicator
C
C     INFO    (output) INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value.
C
C     METHOD
C
C     Let P be a permutation matrix, and D a diagonal matrix of scaling
C     factors, both of order N. The routine computes
C                     -1
C        A <-- P D A D  P'.
C
C     where the permutation and scaling factors are encoded in the
C     array SCALE.
C
C     REFERENCES
C
C     None.
C
C     NUMERICAL ASPECTS
C                               2
C     The algorithm requires O(N ) operations.
C
C     CONTRIBUTORS
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Apr. 1997.
C     Supersedes Release 2.0 routine MB05CY.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Elementary matrix operations, matrix algebra, matrix operations.
C
C     ******************************************************************
C
      DOUBLE PRECISION  ONE
      PARAMETER         ( ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         JOB
      INTEGER           IGH, INFO, LDA, LOW, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), SCALE(*)
C     .. Local Scalars ..
      INTEGER           I, II, J, K
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          DSCAL, DSWAP, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Executable Statements ..
C
C     Test the input scalar arguments.
C
      INFO = 0
      IF( .NOT.LSAME( JOB, 'N' ) .AND. .NOT.LSAME( JOB, 'P' ) .AND.
     $    .NOT.LSAME( JOB, 'S' ) .AND. .NOT.LSAME( JOB, 'B' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 )THEN
         INFO = -2
      ELSE IF( LOW.LT.1 .OR. LOW.GT.MAX( 1, N ) ) THEN
         INFO = -3
      ELSE IF( IGH.LT.MIN( LOW, N ) .OR. IGH.GT.N ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = -6
      END IF
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB05OY', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 .OR. LSAME( JOB, 'N' ) )
     $   RETURN
C
      IF ( .NOT.LSAME( JOB, 'P' ) .AND. IGH.NE.LOW ) THEN
C
         DO 20 I = LOW, IGH
            CALL DSCAL( N, SCALE(I), A(I,1), LDA )
   20    CONTINUE
C
         DO 40 J = LOW, IGH
            CALL DSCAL( N, ONE/SCALE(J), A(1,J), 1 )
   40    CONTINUE
C
      END IF
C
      IF( .NOT.LSAME( JOB, 'S' ) ) THEN
C
         DO 60 II = 1, N
            I = II
            IF ( I.LT.LOW .OR. I.GT.IGH ) THEN
               IF ( I.LT.LOW ) I = LOW - II
               K = SCALE(I)
               IF ( K.NE.I ) THEN
                  CALL DSWAP( N, A(I,1), LDA, A(K,1), LDA )
                  CALL DSWAP( N, A(1,I), 1, A(1,K), 1 )
               END IF
            END IF
   60    CONTINUE
C
      END IF
C
      RETURN
C *** Last line of MB05OY ***
      END
