      SUBROUTINE SB03MD( DICO, JOB, FACT, TRANA, N, A, LDA, U, LDU, C,
     $                   LDC, SCALE, SEP, FERR, WR, WI, IWORK, DWORK,
     $                   LDWORK, INFO )
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
      CHARACTER         DICO, FACT, JOB, TRANA
      INTEGER           INFO, LDA, LDC, LDU, LDWORK, N
      DOUBLE PRECISION  FERR, SCALE, SEP
      INTEGER           IWORK( * )
      DOUBLE PRECISION  A( LDA, * ), C( LDC, * ), DWORK( * ),
     $                  U( LDU, * ), WI( * ), WR( * )
      LOGICAL           CONT, LQUERY, NOFACT, NOTA, WANTBH, WANTSP,
     $                  WANTX
      CHARACTER         NOTRA, NTRNST, TRANST, UPLO
      INTEGER           I, IERR, KASE, LWA, MINWRK, NN, NN2, SDIM
      DOUBLE PRECISION  EPS, EST, SCALEF
      LOGICAL           BWORK( 1 )
      INTEGER           ISAVE( 3 )
      LOGICAL           LSAME, SELECT
      DOUBLE PRECISION  DLAMCH, DLANHS
      EXTERNAL          DLAMCH, DLANHS, LSAME, SELECT
      EXTERNAL          DCOPY, DGEES, DLACN2, MB01RD, SB03MX, SB03MY,
     $                  XERBLA
      INTRINSIC         DBLE, INT, MAX
      CONT   = LSAME( DICO,  'C' )
      WANTX  = LSAME( JOB,   'X' )
      WANTSP = LSAME( JOB,   'S' )
      WANTBH = LSAME( JOB,   'B' )
      NOFACT = LSAME( FACT,  'N' )
      NOTA   = LSAME( TRANA, 'N' )
      LQUERY = LDWORK.EQ.-1
      NN  = N*N
      NN2 = 2*NN
      INFO = 0
      IF( .NOT.CONT .AND. .NOT.LSAME( DICO, 'D' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.WANTBH .AND. .NOT.WANTSP .AND. .NOT.WANTX ) THEN
         INFO = -2
      ELSE IF( .NOT.NOFACT .AND. .NOT.LSAME( FACT, 'F' ) ) THEN
         INFO = -3
      ELSE IF( .NOT.NOTA .AND. .NOT.LSAME( TRANA, 'T' ) .AND.
     $                         .NOT.LSAME( TRANA, 'C' ) ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDU.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( WANTSP .AND. LDC.LT.1 .OR.
     $    .NOT.WANTSP .AND. LDC.LT.MAX( 1, N ) ) THEN
         INFO = -11
      ELSE
         IF ( WANTX ) THEN
            IF ( NOFACT ) THEN
               MINWRK = MAX( NN, 3*N )
            ELSE IF ( CONT ) THEN
               MINWRK = NN
            ELSE
               MINWRK = MAX( NN, 2*N )
            END IF
         ELSE
            IF ( CONT ) THEN
               IF ( NOFACT ) THEN
                  MINWRK = MAX( NN2, 3*N )
               ELSE
                  MINWRK = NN2
               END IF
            ELSE
               MINWRK = NN2 + 2*N
            END IF
         END IF
         MINWRK = MAX( 1, MINWRK )
         IF( LQUERY ) THEN
            IF( NOFACT ) THEN
               CALL DGEES( 'Vectors', 'Not ordered', SELECT, N, A, LDA,
     $                     SDIM, WR, WI, U, LDU, DWORK, -1, BWORK,
     $                     INFO )
               LWA = MAX( MINWRK, INT( DWORK( 1 ) ) )
            ELSE
               LWA = MINWRK
            END IF
         ELSE IF( LDWORK.LT.MINWRK ) THEN
            INFO = -19
         END IF
      END IF
      IF ( INFO.NE.0 ) THEN
         CALL XERBLA( 'SB03MD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         DWORK(1) = LWA
         RETURN
      END IF
      IF( N.EQ.0 ) THEN
         SCALE = ONE
         IF( .NOT.WANTX )
     $      SEP   = ZERO
         IF( WANTBH )
     $      FERR  = ZERO
         DWORK(1) = ONE
         RETURN
      END IF
      LWA = 0
      IF( NOFACT ) THEN
         CALL DGEES( 'Vectors', 'Not ordered', SELECT, N, A, LDA, SDIM,
     $               WR, WI, U, LDU, DWORK, LDWORK, BWORK, INFO )
         IF( INFO.GT.0 )
     $      RETURN
         LWA = INT( DWORK( 1 ) )
      END IF
      IF( .NOT.WANTSP ) THEN
         NTRNST = 'N'
         TRANST = 'T'
         UPLO   = 'U'
         CALL MB01RD( UPLO, TRANST, N, N, ZERO, ONE, C, LDC, U, LDU, C,
     $                LDC, DWORK, LDWORK, INFO )
         DO 10 I = 2, N
            CALL DCOPY( I-1, C(1,I), 1, C(I,1), LDC )
   10    CONTINUE
         LWA = MAX( LWA, NN )
         IF ( CONT ) THEN
            CALL SB03MY( TRANA, N, A, LDA, C, LDC, SCALE, INFO )
         ELSE
            CALL SB03MX( TRANA, N, A, LDA, C, LDC, SCALE, DWORK, INFO )
         END IF
         IF( INFO.GT.0 )
     $      INFO = N + 1
         CALL MB01RD( UPLO, NTRNST, N, N, ZERO, ONE, C, LDC, U, LDU, C,
     $                LDC, DWORK, LDWORK, IERR )
         DO 20 I = 2, N
            CALL DCOPY( I-1, C(1,I), 1, C(I,1), LDC )
   20    CONTINUE
      END IF
      IF( .NOT.WANTX ) THEN
         IF( NOTA ) THEN
            NOTRA = 'T'
         ELSE
            NOTRA = 'N'
         END IF
         EST = ZERO
         KASE = 0
   30    CONTINUE
         CALL DLACN2( NN, DWORK(NN+1), DWORK, IWORK, EST, KASE, ISAVE )
         IF( KASE.NE.0 ) THEN
            IF( KASE.EQ.1 ) THEN
               IF( CONT ) THEN
                  CALL SB03MY( TRANA, N, A, LDA, DWORK, N, SCALEF,
     $                         IERR )
               ELSE
                  CALL SB03MX( TRANA, N, A, LDA, DWORK, N, SCALEF,
     $                         DWORK(NN2+1), IERR )
               END IF
            ELSE
               IF( CONT ) THEN
                  CALL SB03MY( NOTRA, N, A, LDA, DWORK, N, SCALEF,
     $                         IERR )
               ELSE
                  CALL SB03MX( NOTRA, N, A, LDA, DWORK, N, SCALEF,
     $                         DWORK(NN2+1), IERR )
               END IF
            END IF
            GO TO 30
         END IF
         SEP = SCALEF / EST
         IF( WANTBH ) THEN
            EPS = DLAMCH( 'P' )
            IF ( CONT ) THEN
               FERR = EPS*DLANHS( 'Frobenius', N, A, LDA, DWORK )/SEP
            ELSE
               FERR = EPS*DLANHS( 'Frobenius', N, A, LDA, DWORK )**2/SEP
            END IF
         END IF
      END IF
      DWORK( 1 ) = DBLE( MAX( LWA, MINWRK ) )
      RETURN
      END
      SUBROUTINE MB01RD( UPLO, TRANS, M, N, ALPHA, BETA, R, LDR, A, LDA,
     $                   X, LDX, DWORK, LDWORK, INFO )
      DOUBLE PRECISION  ZERO, ONE, HALF
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0, HALF = 0.5D0 )
      CHARACTER         TRANS, UPLO
      INTEGER           INFO, LDA, LDR, LDWORK, LDX, M, N
      DOUBLE PRECISION  ALPHA, BETA
      DOUBLE PRECISION  A(LDA,*), DWORK(*), R(LDR,*), X(LDX,*)
      CHARACTER*12      NTRAN
      LOGICAL           LTRANS, LUPLO
      INTEGER           J, JWORK, LDW, NROWA
      LOGICAL           LSAME
      EXTERNAL          LSAME
      EXTERNAL          DAXPY, DCOPY, DGEMM, DLACPY, DLASCL, DLASET,
     $                  DSCAL, DTRMM, XERBLA
      INTRINSIC         MAX
      INFO = 0
      LUPLO  = LSAME( UPLO,  'U' )
      LTRANS = LSAME( TRANS, 'T' ) .OR. LSAME( TRANS, 'C' )
      IF ( LTRANS ) THEN
         NROWA = N
         NTRAN = 'No transpose'
      ELSE
         NROWA = M
         NTRAN = 'Transpose'
      END IF
      LDW = MAX( 1, M )
      IF(      ( .NOT.LUPLO  ).AND.( .NOT.LSAME( UPLO,  'L' ) ) )THEN
         INFO = -1
      ELSE IF( ( .NOT.LTRANS ).AND.( .NOT.LSAME( TRANS, 'N' ) ) )THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDR.LT.LDW ) THEN
         INFO = -8
      ELSE IF( LDA.LT.MAX( 1, NROWA ) ) THEN
         INFO = -10
      ELSE IF( LDX.LT.MAX( 1, N ) ) THEN
         INFO = -12
      ELSE IF( ( BETA.NE.ZERO .AND. LDWORK.LT.MAX( 1, M*N ) )
     $     .OR.( BETA.EQ.ZERO .AND. LDWORK.LT.1 ) ) THEN
         INFO = -14
      END IF
      IF ( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB01RD', -INFO )
         RETURN
      END IF
      CALL DSCAL( N, HALF, X, LDX+1 )
      IF ( M.EQ.0 )
     $   RETURN
      IF ( BETA.EQ.ZERO .OR. N.EQ.0 ) THEN
         IF ( ALPHA.EQ.ZERO ) THEN
            CALL DLASET( UPLO, M, M, ZERO, ZERO, R, LDR )
         ELSE
            IF ( ALPHA.NE.ONE )
     $         CALL DLASCL( UPLO, 0, 0, ONE, ALPHA, M, M, R, LDR, INFO )
         END IF
         RETURN
      END IF
      IF( LTRANS ) THEN
         JWORK = 1
         DO 10 J = 1, N
            CALL DCOPY( M, A(J,1), LDA, DWORK(JWORK), 1 )
            JWORK = JWORK + LDW
 10      CONTINUE
      ELSE
         CALL DLACPY( 'Full', M, N, A, LDA, DWORK, LDW )
      END IF
      CALL DTRMM( 'Right', UPLO, 'No transpose', 'Non-unit', M, N, BETA,
     $            X, LDX, DWORK, LDW )
      IF ( ALPHA.NE.ZERO ) THEN
         IF ( M.GT.1 ) THEN
            IF ( LUPLO ) THEN
               CALL DLASET( 'Lower', M-1, M-1, ZERO, ZERO, R(2,1), LDR )
            ELSE
               CALL DLASET( 'Upper', M-1, M-1, ZERO, ZERO, R(1,2), LDR )
            END IF
         END IF
         CALL DSCAL( M, HALF, R, LDR+1 )
      END IF
      CALL DGEMM( 'No transpose', NTRAN, M, M, N, ONE, DWORK, LDW, A,
     $            LDA, ALPHA, R, LDR )
      IF( LUPLO ) THEN
         DO 20 J = 1, M
            CALL DAXPY( J, ONE, R(J,1), LDR, R(1,J), 1 )
   20    CONTINUE
      ELSE
         DO 30 J = 1, M
            CALL DAXPY( J, ONE, R(1,J), 1, R(J,1), LDR )
 30      CONTINUE
      END IF
      RETURN
      END
      SUBROUTINE SB03MX( TRANA, N, A, LDA, C, LDC, SCALE, DWORK, INFO )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      CHARACTER          TRANA
      INTEGER            INFO, LDA, LDC, N
      DOUBLE PRECISION   SCALE
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), DWORK( * )
      LOGICAL            NOTRNA, LUPPER
      INTEGER            IERR, J, K, K1, K2, KNEXT, L, L1, L2, LNEXT,
     $                   MINK1N, MINK2N, MINL1N, MINL2N, NP1
      DOUBLE PRECISION   A11, BIGNUM, DA11, DB, EPS, P11, P12, P21, P22,
     $                   SCALOC, SMIN, SMLNUM, XNORM
      DOUBLE PRECISION   VEC( 2, 2 ), X( 2, 2 )
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT, DLAMCH, DLANHS
      EXTERNAL           DDOT, DLAMCH, DLANHS, LSAME
      EXTERNAL           DLABAD, DLALN2, DSCAL, DSYMV, SB03MV, SB04PX,
     $                   XERBLA
      INTRINSIC          ABS, DBLE, MAX, MIN
      NOTRNA = LSAME( TRANA, 'N' )
      LUPPER = .TRUE.
      INFO = 0
      IF( .NOT.NOTRNA .AND. .NOT.LSAME( TRANA, 'T' ) .AND.
     $                      .NOT.LSAME( TRANA, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LDC.LT.MAX( 1, N ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SB03MX', -INFO )
         RETURN
      END IF
      SCALE = ONE
      IF( N.EQ.0 )
     $   RETURN
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
      SMLNUM = SMLNUM*DBLE( N*N ) / EPS
      BIGNUM = ONE / SMLNUM
      SMIN = MAX( SMLNUM, EPS*DLANHS( 'Max', N, A, LDA, DWORK ) )
      NP1  = N + 1
      IF( NOTRNA ) THEN
         LNEXT = 1
         DO 60 L = 1, N
            IF( L.LT.LNEXT )
     $         GO TO 60
            L1 = L
            L2 = L
            IF( L.LT.N ) THEN
               IF( A( L+1, L ).NE.ZERO )
     $            L2 = L2 + 1
               LNEXT = L2 + 1
            END IF
            DWORK( L1 )   = ZERO
            DWORK( N+L1 ) = ZERO
            CALL DSYMV( 'Lower', L1-1, ONE, C, LDC, A( 1, L1 ), 1, ZERO,
     $                  DWORK, 1 )
            CALL DSYMV( 'Lower', L1-1, ONE, C, LDC, A( 1, L2 ), 1, ZERO,
     $                  DWORK( NP1 ), 1 )
            KNEXT = L
            DO 50 K = L, N
               IF( K.LT.KNEXT )
     $            GO TO 50
               K1 = K
               K2 = K
               IF( K.LT.N ) THEN
                  IF( A( K+1, K ).NE.ZERO )
     $               K2 = K2 + 1
                  KNEXT = K2 + 1
               END IF
               IF( L1.EQ.L2 .AND. K1.EQ.K2 ) THEN
                  DWORK( K1 ) = DDOT( L1-1, C( K1, 1 ), LDC, A( 1, L1 ),
     $                                1 )
                  VEC( 1, 1 ) = C( K1, L1 ) -
     $               ( DDOT( K1, A( 1, K1 ), 1, DWORK, 1 ) + A( L1, L1 )
     $                *DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 ) )
                  SCALOC = ONE
                  A11 = A( K1, K1 )*A( L1, L1 ) - ONE
                  DA11 = ABS( A11 )
                  IF( DA11.LE.SMIN ) THEN
                     A11 = SMIN
                     DA11 = SMIN
                     INFO = 1
                  END IF
                  DB = ABS( VEC( 1, 1 ) )
                  IF( DA11.LT.ONE .AND. DB.GT.ONE ) THEN
                     IF( DB.GT.BIGNUM*DA11 )
     $                  SCALOC = ONE / DB
                  END IF
                  X( 1, 1 ) = ( VEC( 1, 1 )*SCALOC ) / A11
                  IF( SCALOC.NE.ONE ) THEN
                     DO 10 J = 1, N
                        CALL DSCAL( N, SCALOC, C( 1, J ), 1 )
   10                CONTINUE
                     CALL DSCAL( N, SCALOC, DWORK, 1 )
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  IF( K1.NE.L1 ) THEN
                     C( L1, K1 ) = X( 1, 1 )
                  END IF
               ELSE IF( L1.EQ.L2 .AND. K1.NE.K2 ) THEN
                  DWORK( K1 ) = DDOT( L1-1, C( K1, 1 ), LDC, A( 1, L1 ),
     $                                1 )
                  DWORK( K2 ) = DDOT( L1-1, C( K2, 1 ), LDC, A( 1, L1 ),
     $                                1 )
                  VEC( 1, 1 ) = C( K1, L1 ) -
     $               ( DDOT( K2, A( 1, K1 ), 1, DWORK, 1 ) + A( L1, L1 )
     $                *DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 ) )
                  VEC( 2, 1 ) = C( K2, L1 ) -
     $               ( DDOT( K2, A( 1, K2 ), 1, DWORK, 1 ) + A( L1, L1 )
     $                *DDOT( K1-1, A( 1, K2 ), 1, C( 1, L1 ), 1 ) )
                  CALL DLALN2( .TRUE., 2, 1, SMIN, A( L1, L1 ),
     $                         A( K1, K1 ), LDA, ONE, ONE, VEC, 2, ONE,
     $                         ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 1
                  IF( SCALOC.NE.ONE ) THEN
                     DO 20 J = 1, N
                        CALL DSCAL( N, SCALOC, C( 1, J ), 1 )
   20                CONTINUE
                     CALL DSCAL( N, SCALOC, DWORK, 1 )
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K2, L1 ) = X( 2, 1 )
                  C( L1, K1 ) = X( 1, 1 )
                  C( L1, K2 ) = X( 2, 1 )
               ELSE IF( L1.NE.L2 .AND. K1.EQ.K2 ) THEN
                  DWORK( K1 ) = DDOT( L1-1, C( K1, 1 ), LDC, A( 1, L1 ),
     $                                1 )
                  DWORK( N+K1 ) = DDOT( L1-1, C( K1, 1 ), LDC,
     $                                  A( 1, L2 ), 1 )
                  P11 = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  P12 = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L2 ), 1 )
                  VEC( 1, 1 ) = C( K1, L1 ) -
     $               ( DDOT( K1, A( 1, K1 ), 1, DWORK, 1 ) +
     $                 P11*A( L1, L1 ) + P12*A( L2, L1 ) )
                  VEC( 2, 1 ) = C( K1, L2 ) -
     $               ( DDOT( K1, A( 1, K1 ), 1, DWORK( NP1 ), 1 ) +
     $                 P11*A( L1, L2 ) + P12*A( L2, L2 ) )
                  CALL DLALN2( .TRUE., 2, 1, SMIN, A( K1, K1 ),
     $                         A( L1, L1 ), LDA, ONE, ONE, VEC, 2, ONE,
     $                         ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 1
                  IF( SCALOC.NE.ONE ) THEN
                     DO 30 J = 1, N
                        CALL DSCAL( N, SCALOC, C( 1, J ), 1 )
   30                CONTINUE
                     CALL DSCAL( N, SCALOC, DWORK, 1 )
                     CALL DSCAL( N, SCALOC, DWORK( NP1 ), 1 )
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 2, 1 )
                  C( L1, K1 ) = X( 1, 1 )
                  C( L2, K1 ) = X( 2, 1 )
               ELSE IF( L1.NE.L2 .AND. K1.NE.K2 ) THEN
                  DWORK( K1 ) = DDOT( L1-1, C( K1, 1 ), LDC, A( 1, L1 ),
     $                                1 )
                  DWORK( K2 ) = DDOT( L1-1, C( K2, 1 ), LDC, A( 1, L1 ),
     $                                1 )
                  DWORK( N+K1 ) = DDOT( L1-1, C( K1, 1 ), LDC,
     $                                  A( 1, L2 ), 1 )
                  DWORK( N+K2 ) = DDOT( L1-1, C( K2, 1 ), LDC,
     $                                  A( 1, L2 ), 1 )
                  P11 = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  P12 = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L2 ), 1 )
                  P21 = DDOT( K1-1, A( 1, K2 ), 1, C( 1, L1 ), 1 )
                  P22 = DDOT( K1-1, A( 1, K2 ), 1, C( 1, L2 ), 1 )
                  VEC( 1, 1 ) = C( K1, L1 ) -
     $               ( DDOT( K2, A( 1, K1 ), 1, DWORK, 1 ) +
     $                 P11*A( L1, L1 ) + P12*A( L2, L1 ) )
                  VEC( 1, 2 ) = C( K1, L2 ) -
     $               ( DDOT( K2, A( 1, K1 ), 1, DWORK( NP1 ), 1 ) +
     $                 P11*A( L1, L2 ) + P12*A( L2, L2 ) )
                  VEC( 2, 1 ) = C( K2, L1 ) -
     $               ( DDOT( K2, A( 1, K2 ), 1, DWORK, 1 ) +
     $                 P21*A( L1, L1 ) + P22*A( L2, L1 ) )
                  VEC( 2, 2 ) = C( K2, L2 ) -
     $               ( DDOT( K2, A( 1, K2 ), 1, DWORK( NP1 ), 1 ) +
     $                 P21*A( L1, L2 ) + P22*A( L2, L2 ) )
                  IF( K1.EQ.L1 ) THEN
                     CALL SB03MV( .FALSE., LUPPER, A( K1, K1 ), LDA,
     $                            VEC, 2, SCALOC, X, 2, XNORM, IERR )
                     IF( LUPPER ) THEN
                        X( 2, 1 ) = X( 1, 2 )
                     ELSE
                        X( 1, 2 ) = X( 2, 1 )
                     END IF
                  ELSE
                     CALL SB04PX( .TRUE., .FALSE., -1, 2, 2,
     $                            A( K1, K1 ), LDA, A( L1, L1 ), LDA,
     $                            VEC, 2, SCALOC, X, 2, XNORM, IERR )
                  END IF
                  IF( IERR.NE.0 )
     $               INFO = 1
                  IF( SCALOC.NE.ONE ) THEN
                     DO 40 J = 1, N
                        CALL DSCAL( N, SCALOC, C( 1, J ), 1 )
   40                CONTINUE
                     CALL DSCAL( N, SCALOC, DWORK, 1 )
                     CALL DSCAL( N, SCALOC, DWORK( NP1 ), 1 )
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 1, 2 )
                  C( K2, L1 ) = X( 2, 1 )
                  C( K2, L2 ) = X( 2, 2 )
                  IF( K1.NE.L1 ) THEN
                     C( L1, K1 ) = X( 1, 1 )
                     C( L2, K1 ) = X( 1, 2 )
                     C( L1, K2 ) = X( 2, 1 )
                     C( L2, K2 ) = X( 2, 2 )
                  END IF
               END IF
   50       CONTINUE
   60    CONTINUE
      ELSE
         LNEXT = N
         DO 120 L = N, 1, -1
            IF( L.GT.LNEXT )
     $         GO TO 120
            L1 = L
            L2 = L
            IF( L.GT.1 ) THEN
               IF( A( L, L-1 ).NE.ZERO ) THEN
                  L1 = L1 - 1
                  DWORK( L1 ) = ZERO
                  DWORK( N+L1 ) = ZERO
               END IF
               LNEXT = L1 - 1
            END IF
            MINL1N = MIN( L1+1, N )
            MINL2N = MIN( L2+1, N )
            IF( L2.LT.N ) THEN
               CALL DSYMV( 'Upper', N-L2, ONE, C( L2+1, L2+1 ), LDC,
     $                     A( L1, L2+1 ), LDA, ZERO, DWORK( L2+1 ), 1 )
               CALL DSYMV( 'Upper', N-L2, ONE, C( L2+1, L2+1 ), LDC,
     $                     A( L2, L2+1 ), LDA, ZERO, DWORK( NP1+L2 ), 1)
            END IF
            KNEXT = L
            DO 110 K = L, 1, -1
               IF( K.GT.KNEXT )
     $            GO TO 110
               K1 = K
               K2 = K
               IF( K.GT.1 ) THEN
                  IF( A( K, K-1 ).NE.ZERO )
     $               K1 = K1 - 1
                  KNEXT = K1 - 1
               END IF
               MINK1N = MIN( K1+1, N )
               MINK2N = MIN( K2+1, N )
               IF( L1.EQ.L2 .AND. K1.EQ.K2 ) THEN
                  DWORK( K1 ) = DDOT( N-L1, C( K1, MINL1N ), LDC,
     $                                A( L1, MINL1N ), LDA )
                  VEC( 1, 1 ) = C( K1, L1 ) -
     $               ( DDOT( N-K1+1, A( K1, K1 ), LDA, DWORK( K1 ), 1 )
     $               + DDOT( N-K1, A( K1, MINK1N ), LDA,
     $                       C( MINK1N, L1 ), 1 )*A( L1, L1 ) )
                  SCALOC = ONE
                  A11 = A( K1, K1 )*A( L1, L1 ) - ONE
                  DA11 = ABS( A11 )
                  IF( DA11.LE.SMIN ) THEN
                     A11 = SMIN
                     DA11 = SMIN
                     INFO = 1
                  END IF
                  DB = ABS( VEC( 1, 1 ) )
                  IF( DA11.LT.ONE .AND. DB.GT.ONE ) THEN
                     IF( DB.GT.BIGNUM*DA11 )
     $                  SCALOC = ONE / DB
                  END IF
                  X( 1, 1 ) = ( VEC( 1, 1 )*SCALOC ) / A11
                  IF( SCALOC.NE.ONE ) THEN
                     DO 70 J = 1, N
                        CALL DSCAL( N, SCALOC, C( 1, J ), 1 )
   70                CONTINUE
                     CALL DSCAL( N, SCALOC, DWORK, 1 )
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  IF( K1.NE.L1 ) THEN
                     C( L1, K1 ) = X( 1, 1 )
                  END IF
               ELSE IF( L1.EQ.L2 .AND. K1.NE.K2 ) THEN
                  DWORK( K1 ) = DDOT( N-L1, C( K1, MINL1N ), LDC,
     $                                A( L1, MINL1N ), LDA )
                  DWORK( K2 ) = DDOT( N-L1, C( K2, MINL1N ), LDC,
     $                                A( L1, MINL1N ), LDA )
                  VEC( 1, 1 ) = C( K1, L1 ) -
     $               ( DDOT( NP1-K1, A( K1, K1 ), LDA, DWORK( K1 ), 1 )
     $               + DDOT( N-K2, A( K1, MINK2N ), LDA,
     $                       C( MINK2N, L1 ), 1 )*A( L1, L1 ) )
                  VEC( 2, 1 ) = C( K2, L1 ) -
     $               ( DDOT( NP1-K1, A( K2, K1 ), LDA, DWORK( K1 ), 1 )
     $               + DDOT( N-K2, A( K2, MINK2N ), LDA,
     $                       C( MINK2N, L1 ), 1 )*A( L1, L1 ) )
                  CALL DLALN2( .FALSE., 2, 1, SMIN, A( L1, L1 ),
     $                         A( K1, K1 ), LDA, ONE, ONE, VEC, 2, ONE,
     $                         ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 1
                  IF( SCALOC.NE.ONE ) THEN
                     DO 80 J = 1, N
                        CALL DSCAL( N, SCALOC, C( 1, J ), 1 )
   80                CONTINUE
                     CALL DSCAL( N, SCALOC, DWORK, 1 )
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K2, L1 ) = X( 2, 1 )
                  C( L1, K1 ) = X( 1, 1 )
                  C( L1, K2 ) = X( 2, 1 )
               ELSE IF( L1.NE.L2 .AND. K1.EQ.K2 ) THEN
                  DWORK( K1 ) = DDOT( N-L2, C( K1, MINL2N ), LDC,
     $                                A( L1, MINL2N ), LDA )
                  DWORK( N+K1 ) = DDOT( N-L2, C( K1, MINL2N ), LDC,
     $                                  A( L2, MINL2N ), LDA )
                  P11 = DDOT( N-K1, A( K1, MINK1N ), LDA,
     $                        C( MINK1N, L1 ), 1 )
                  P12 = DDOT( N-K1, A( K1, MINK1N ), LDA,
     $                        C( MINK1N, L2 ), 1 )
                  VEC( 1, 1 ) = C( K1, L1 ) -
     $               ( DDOT( NP1-K1, A( K1, K1 ), LDA, DWORK( K1 ), 1 )
     $               + P11*A( L1, L1 ) + P12*A( L1, L2 ) )
                  VEC( 2, 1 ) = C( K1, L2 ) -
     $               ( DDOT( NP1-K1, A( K1, K1 ), LDA, DWORK( N+K1 ), 1)
     $               + P11*A( L2, L1 ) + P12*A( L2, L2 ) )
                  CALL DLALN2( .FALSE., 2, 1, SMIN, A( K1, K1 ),
     $                         A( L1, L1 ), LDA, ONE, ONE, VEC, 2, ONE,
     $                         ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 1
                  IF( SCALOC.NE.ONE ) THEN
                     DO 90 J = 1, N
                        CALL DSCAL( N, SCALOC, C( 1, J ), 1 )
   90                CONTINUE
                     CALL DSCAL( N, SCALOC, DWORK, 1 )
                     CALL DSCAL( N, SCALOC, DWORK( NP1 ), 1 )
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 2, 1 )
                  C( L1, K1 ) = X( 1, 1 )
                  C( L2, K1 ) = X( 2, 1 )
               ELSE IF( L1.NE.L2 .AND. K1.NE.K2 ) THEN
                  DWORK( K1 ) = DDOT( N-L2, C( K1, MINL2N ), LDC,
     $                                A( L1, MINL2N ), LDA )
                  DWORK( K2 ) = DDOT( N-L2, C( K2, MINL2N ), LDC,
     $                                A( L1, MINL2N ), LDA )
                  DWORK( N+K1 ) = DDOT( N-L2, C( K1, MINL2N ), LDC,
     $                                  A( L2, MINL2N ), LDA )
                  DWORK( N+K2 ) = DDOT( N-L2, C( K2, MINL2N ), LDC,
     $                                  A( L2, MINL2N ), LDA )
                  P11 = DDOT( N-K2, A( K1, MINK2N ), LDA,
     $                        C( MINK2N, L1 ), 1 )
                  P12 = DDOT( N-K2, A( K1, MINK2N ), LDA,
     $                        C( MINK2N, L2 ), 1 )
                  P21 = DDOT( N-K2, A( K2, MINK2N ), LDA,
     $                        C( MINK2N, L1 ), 1 )
                  P22 = DDOT( N-K2, A( K2, MINK2N ), LDA,
     $                        C( MINK2N, L2 ), 1 )
                  VEC( 1, 1 ) = C( K1, L1 ) -
     $               ( DDOT( NP1-K1, A( K1, K1 ), LDA, DWORK( K1 ), 1 )
     $               + P11*A( L1, L1 ) + P12*A( L1, L2 ) )
                  VEC( 1, 2 ) = C( K1, L2 ) -
     $               ( DDOT( NP1-K1, A( K1, K1 ), LDA, DWORK( N+K1 ),
     $                       1) + P11*A( L2, L1 ) + P12*A( L2, L2 ) )
                  VEC( 2, 1 ) = C( K2, L1 ) -
     $               ( DDOT( NP1-K1, A( K2, K1 ), LDA, DWORK( K1 ),
     $                       1) + P21*A( L1, L1 ) + P22*A( L1, L2 ) )
                  VEC( 2, 2 ) = C( K2, L2 ) -
     $               ( DDOT( NP1-K1, A( K2, K1 ), LDA, DWORK( N+K1 ), 1)
     $               + P21*A( L2, L1 ) + P22*A( L2, L2 ) )
                  IF( K1.EQ.L1 ) THEN
                     CALL SB03MV( .TRUE., LUPPER, A( K1, K1 ), LDA, VEC,
     $                            2, SCALOC, X, 2, XNORM, IERR )
                     IF( LUPPER ) THEN
                        X( 2, 1 ) = X( 1, 2 )
                     ELSE
                        X( 1, 2 ) = X( 2, 1 )
                     END IF
                  ELSE
                     CALL SB04PX( .FALSE., .TRUE., -1, 2, 2,
     $                            A( K1, K1 ), LDA, A( L1, L1 ), LDA,
     $                            VEC, 2, SCALOC, X, 2, XNORM, IERR )
                  END IF
                  IF( IERR.NE.0 )
     $               INFO = 1
                  IF( SCALOC.NE.ONE ) THEN
                     DO 100 J = 1, N
                        CALL DSCAL( N, SCALOC, C( 1, J ), 1 )
  100                CONTINUE
                     CALL DSCAL( N, SCALOC, DWORK, 1 )
                     CALL DSCAL( N, SCALOC, DWORK( NP1 ), 1 )
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 1, 2 )
                  C( K2, L1 ) = X( 2, 1 )
                  C( K2, L2 ) = X( 2, 2 )
                  IF( K1.NE.L1 ) THEN
                     C( L1, K1 ) = X( 1, 1 )
                     C( L2, K1 ) = X( 1, 2 )
                     C( L1, K2 ) = X( 2, 1 )
                     C( L2, K2 ) = X( 2, 2 )
                  END IF
               END IF
  110       CONTINUE
  120    CONTINUE
      END IF
      RETURN
      END
      SUBROUTINE SB03MY( TRANA, N, A, LDA, C, LDC, SCALE, INFO )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      CHARACTER          TRANA
      INTEGER            INFO, LDA, LDC, N
      DOUBLE PRECISION   SCALE
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * )
      LOGICAL            NOTRNA, LUPPER
      INTEGER            IERR, J, K, K1, K2, KNEXT, L, L1, L2, LNEXT,
     $                   MINK1N, MINK2N, MINL1N, MINL2N
      DOUBLE PRECISION   A11, BIGNUM, DA11, DB, EPS, SCALOC, SMIN,
     $                   SMLNUM, XNORM
      DOUBLE PRECISION   DUM( 1 ), VEC( 2, 2 ), X( 2, 2 )
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT, DLAMCH, DLANHS
      EXTERNAL           DDOT, DLAMCH, DLANHS, LSAME
      EXTERNAL           DLABAD, DLALN2, DLASY2, DSCAL, SB03MW, XERBLA
      INTRINSIC          ABS, DBLE, MAX, MIN
      NOTRNA = LSAME( TRANA, 'N' )
      LUPPER = .TRUE.
      INFO = 0
      IF( .NOT.NOTRNA .AND. .NOT.LSAME( TRANA, 'T' ) .AND.
     $                      .NOT.LSAME( TRANA, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LDC.LT.MAX( 1, N ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SB03MY', -INFO )
         RETURN
      END IF
      SCALE = ONE
      IF( N.EQ.0 )
     $   RETURN
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
      SMLNUM = SMLNUM*DBLE( N*N ) / EPS
      BIGNUM = ONE / SMLNUM
      SMIN = MAX( SMLNUM, EPS*DLANHS( 'Max', N, A, LDA, DUM ) )
      IF( NOTRNA ) THEN
         LNEXT = 1
         DO 60 L = 1, N
            IF( L.LT.LNEXT )
     $         GO TO 60
            L1 = L
            L2 = L
            IF( L.LT.N ) THEN
               IF( A( L+1, L ).NE.ZERO )
     $            L2 = L2 + 1
               LNEXT = L2 + 1
            END IF
            KNEXT = L
            DO 50 K = L, N
               IF( K.LT.KNEXT )
     $            GO TO 50
               K1 = K
               K2 = K
               IF( K.LT.N ) THEN
                  IF( A( K+1, K ).NE.ZERO )
     $               K2 = K2 + 1
                  KNEXT = K2 + 1
               END IF
               IF( L1.EQ.L2 .AND. K1.EQ.K2 ) THEN
                  VEC( 1, 1 ) = C( K1, L1 ) -
     $               ( DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 ) +
     $                 DDOT( L1-1, C( K1, 1 ), LDC, A( 1, L1 ), 1 ) )
                  SCALOC = ONE
                  A11 = A( K1, K1 ) + A( L1, L1 )
                  DA11 = ABS( A11 )
                  IF( DA11.LE.SMIN ) THEN
                     A11 = SMIN
                     DA11 = SMIN
                     INFO = 1
                  END IF
                  DB = ABS( VEC( 1, 1 ) )
                  IF( DA11.LT.ONE .AND. DB.GT.ONE ) THEN
                     IF( DB.GT.BIGNUM*DA11 )
     $                  SCALOC = ONE / DB
                  END IF
                  X( 1, 1 ) = ( VEC( 1, 1 )*SCALOC ) / A11
                  IF( SCALOC.NE.ONE ) THEN
                     DO 10 J = 1, N
                        CALL DSCAL( N, SCALOC, C( 1, J ), 1 )
   10                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  IF( K1.NE.L1 ) THEN
                     C( L1, K1 ) = X( 1, 1 )
                  END IF
               ELSE IF( L1.EQ.L2 .AND. K1.NE.K2 ) THEN
                  VEC( 1, 1 ) = C( K1, L1 ) -
     $               ( DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 ) +
     $                 DDOT( L1-1, C( K1, 1 ), LDC, A( 1, L1 ), 1 ) )
                  VEC( 2, 1 ) = C( K2, L1 ) -
     $               ( DDOT( K1-1, A( 1, K2 ), 1, C( 1, L1 ), 1 ) +
     $                 DDOT( L1-1, C( K2, 1 ), LDC, A( 1, L1 ), 1 ) )
                  CALL DLALN2( .TRUE., 2, 1, SMIN, ONE, A( K1, K1 ),
     $                         LDA, ONE, ONE, VEC, 2, -A( L1, L1 ),
     $                         ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 1
                  IF( SCALOC.NE.ONE ) THEN
                     DO 20 J = 1, N
                        CALL DSCAL( N, SCALOC, C( 1, J ), 1 )
   20                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K2, L1 ) = X( 2, 1 )
                  C( L1, K1 ) = X( 1, 1 )
                  C( L1, K2 ) = X( 2, 1 )
               ELSE IF( L1.NE.L2 .AND. K1.EQ.K2 ) THEN
                  VEC( 1, 1 ) = C( K1, L1 ) -
     $               ( DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 ) +
     $                 DDOT( L1-1, C( K1, 1 ), LDC, A( 1, L1 ), 1 ) )
                  VEC( 2, 1 ) = C( K1, L2 ) -
     $               ( DDOT( K1-1, A( 1, K1 ), 1, C( 1, L2 ), 1 ) +
     $                 DDOT( L1-1, C( K1, 1 ), LDC, A( 1, L2 ), 1 ) )
                  CALL DLALN2( .TRUE., 2, 1, SMIN, ONE, A( L1, L1 ),
     $                         LDA, ONE, ONE, VEC, 2, -A( K1, K1 ),
     $                         ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 1
                  IF( SCALOC.NE.ONE ) THEN
                     DO 30 J = 1, N
                        CALL DSCAL( N, SCALOC, C( 1, J ), 1 )
   30                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 2, 1 )
                  C( L1, K1 ) = X( 1, 1 )
                  C( L2, K1 ) = X( 2, 1 )
               ELSE IF( L1.NE.L2 .AND. K1.NE.K2 ) THEN
                  VEC( 1, 1 ) = C( K1, L1 ) -
     $               ( DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 ) +
     $                 DDOT( L1-1, C( K1, 1 ), LDC, A( 1, L1 ), 1 ) )
                  VEC( 1, 2 ) = C( K1, L2 ) -
     $               ( DDOT( K1-1, A( 1, K1 ), 1, C( 1, L2 ), 1 ) +
     $                 DDOT( L1-1, C( K1, 1 ), LDC, A( 1, L2 ), 1 ) )
                  VEC( 2, 1 ) = C( K2, L1 ) -
     $               ( DDOT( K1-1, A( 1, K2 ), 1, C( 1, L1 ), 1 ) +
     $                 DDOT( L1-1, C( K2, 1 ), LDC, A( 1, L1 ), 1 ) )
                  VEC( 2, 2 ) = C( K2, L2 ) -
     $               ( DDOT( K1-1, A( 1, K2 ), 1, C( 1, L2 ), 1 ) +
     $                 DDOT( L1-1, C( K2, 1 ), LDC, A( 1, L2 ), 1 ) )
                  IF( K1.EQ.L1 ) THEN
                     CALL SB03MW( .FALSE., LUPPER, A( K1, K1 ), LDA,
     $                            VEC, 2, SCALOC, X, 2, XNORM, IERR )
                     IF( LUPPER ) THEN
                        X( 2, 1 ) = X( 1, 2 )
                     ELSE
                        X( 1, 2 ) = X( 2, 1 )
                     END IF
                  ELSE
                     CALL DLASY2( .TRUE., .FALSE., 1, 2, 2, A( K1, K1 ),
     $                            LDA, A( L1, L1 ), LDA, VEC, 2, SCALOC,
     $                            X, 2, XNORM, IERR )
                  END IF
                  IF( IERR.NE.0 )
     $               INFO = 1
                  IF( SCALOC.NE.ONE ) THEN
                     DO 40 J = 1, N
                        CALL DSCAL( N, SCALOC, C( 1, J ), 1 )
   40                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 1, 2 )
                  C( K2, L1 ) = X( 2, 1 )
                  C( K2, L2 ) = X( 2, 2 )
                  IF( K1.NE.L1 ) THEN
                     C( L1, K1 ) = X( 1, 1 )
                     C( L2, K1 ) = X( 1, 2 )
                     C( L1, K2 ) = X( 2, 1 )
                     C( L2, K2 ) = X( 2, 2 )
                  END IF
               END IF
   50       CONTINUE
   60    CONTINUE
      ELSE
         LNEXT = N
         DO 120 L = N, 1, -1
            IF( L.GT.LNEXT )
     $         GO TO 120
            L1 = L
            L2 = L
            IF( L.GT.1 ) THEN
               IF( A( L, L-1 ).NE.ZERO )
     $            L1 = L1 - 1
               LNEXT = L1 - 1
            END IF
            MINL1N = MIN( L1+1, N )
            MINL2N = MIN( L2+1, N )
            KNEXT = L
            DO 110 K = L, 1, -1
               IF( K.GT.KNEXT )
     $            GO TO 110
               K1 = K
               K2 = K
               IF( K.GT.1 ) THEN
                  IF( A( K, K-1 ).NE.ZERO )
     $               K1 = K1 - 1
                  KNEXT = K1 - 1
               END IF
               MINK1N = MIN( K1+1, N )
               MINK2N = MIN( K2+1, N )
               IF( L1.EQ.L2 .AND. K1.EQ.K2 ) THEN
                  VEC( 1, 1 ) = C( K1, L1 ) -
     $               ( DDOT( N-K1, A( K1, MINK1N ), LDA,
     $                       C( MINK1N, L1 ), 1 ) +
     $                 DDOT( N-L1, C( K1, MINL1N ), LDC,
     $                       A( L1, MINL1N ), LDA ) )
                  SCALOC = ONE
                  A11 = A( K1, K1 ) + A( L1, L1 )
                  DA11 = ABS( A11 )
                  IF( DA11.LE.SMIN ) THEN
                     A11 = SMIN
                     DA11 = SMIN
                     INFO = 1
                  END IF
                  DB = ABS( VEC( 1, 1 ) )
                  IF( DA11.LT.ONE .AND. DB.GT.ONE ) THEN
                     IF( DB.GT.BIGNUM*DA11 )
     $                  SCALOC = ONE / DB
                  END IF
                  X( 1, 1 ) = ( VEC( 1, 1 )*SCALOC ) / A11
                  IF( SCALOC.NE.ONE ) THEN
                     DO 70 J = 1, N
                        CALL DSCAL( N, SCALOC, C( 1, J ), 1 )
   70                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  IF( K1.NE.L1 ) THEN
                     C( L1, K1 ) = X( 1, 1 )
                  END IF
               ELSE IF( L1.EQ.L2 .AND. K1.NE.K2 ) THEN
                  VEC( 1, 1 ) = C( K1, L1 ) -
     $               ( DDOT( N-K2, A( K1, MINK2N ), LDA,
     $                       C( MINK2N, L1 ), 1 ) +
     $                 DDOT( N-L2, C( K1, MINL2N ), LDC,
     $                       A( L1, MINL2N ), LDA ) )
                  VEC( 2, 1 ) = C( K2, L1 ) -
     $               ( DDOT( N-K2, A( K2, MINK2N ), LDA,
     $                     C( MINK2N, L1 ), 1 ) +
     $                 DDOT( N-L2, C( K2, MINL2N ), LDC,
     $                     A( L1, MINL2N ), LDA ) )
                  CALL DLALN2( .FALSE., 2, 1, SMIN, ONE, A( K1, K1 ),
     $                         LDA, ONE, ONE, VEC, 2, -A( L1, L1 ),
     $                         ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 1
                  IF( SCALOC.NE.ONE ) THEN
                     DO 80 J = 1, N
                        CALL DSCAL( N, SCALOC, C( 1, J ), 1 )
   80                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K2, L1 ) = X( 2, 1 )
                  C( L1, K1 ) = X( 1, 1 )
                  C( L1, K2 ) = X( 2, 1 )
               ELSE IF( L1.NE.L2 .AND. K1.EQ.K2 ) THEN
                  VEC( 1, 1 ) = C( K1, L1 ) -
     $               ( DDOT( N-K1, A( K1, MINK1N ), LDA,
     $                       C( MINK1N, L1 ), 1 ) +
     $                 DDOT( N-L2, C( K1, MINL2N ), LDC,
     $                       A( L1, MINL2N ), LDA ) )
                  VEC( 2, 1 ) = C( K1, L2 ) -
     $               ( DDOT( N-K1, A( K1, MINK1N ), LDA,
     $                       C( MINK1N, L2 ), 1 ) +
     $                 DDOT( N-L2, C( K1, MINL2N ), LDC,
     $                       A( L2, MINL2N ), LDA ) )
                  CALL DLALN2( .FALSE., 2, 1, SMIN, ONE, A( L1, L1 ),
     $                         LDA, ONE, ONE, VEC, 2, -A( K1, K1 ),
     $                         ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 1
                  IF( SCALOC.NE.ONE ) THEN
                     DO 90 J = 1, N
                        CALL DSCAL( N, SCALOC, C( 1, J ), 1 )
   90                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 2, 1 )
                  C( L1, K1 ) = X( 1, 1 )
                  C( L2, K1 ) = X( 2, 1 )
               ELSE IF( L1.NE.L2 .AND. K1.NE.K2 ) THEN
                  VEC( 1, 1 ) = C( K1, L1 ) -
     $               ( DDOT( N-K2, A( K1, MINK2N ), LDA,
     $                       C( MINK2N, L1 ), 1 ) +
     $                 DDOT( N-L2, C( K1, MINL2N ), LDC,
     $                       A( L1, MINL2N ), LDA ) )
                  VEC( 1, 2 ) = C( K1, L2 ) -
     $               ( DDOT( N-K2, A( K1, MINK2N ), LDA,
     $                       C( MINK2N, L2 ), 1 ) +
     $                 DDOT( N-L2, C( K1, MINL2N ), LDC,
     $                       A( L2, MINL2N ), LDA ) )
                  VEC( 2, 1 ) = C( K2, L1 ) -
     $               ( DDOT( N-K2, A( K2, MINK2N ), LDA,
     $                       C( MINK2N, L1 ), 1 ) +
     $                 DDOT( N-L2, C( K2, MINL2N ), LDC,
     $                       A( L1, MINL2N ), LDA ) )
                  VEC( 2, 2 ) = C( K2, L2 ) -
     $               ( DDOT( N-K2, A( K2, MINK2N ), LDA,
     $                       C( MINK2N, L2 ), 1 ) +
     $                 DDOT( N-L2, C( K2, MINL2N ), LDC,
     $                       A( L2, MINL2N ), LDA ) )
                  IF( K1.EQ.L1 ) THEN
                     CALL SB03MW( .TRUE., LUPPER, A( K1, K1 ), LDA, VEC,
     $                            2, SCALOC, X, 2, XNORM, IERR )
                     IF( LUPPER ) THEN
                        X( 2, 1 ) = X( 1, 2 )
                     ELSE
                        X( 1, 2 ) = X( 2, 1 )
                     END IF
                  ELSE
                     CALL DLASY2( .FALSE., .TRUE., 1, 2, 2, A( K1, K1 ),
     $                            LDA, A( L1, L1 ), LDA, VEC, 2, SCALOC,
     $                            X, 2, XNORM, IERR )
                  END IF
                  IF( IERR.NE.0 )
     $               INFO = 1
                  IF( SCALOC.NE.ONE ) THEN
                     DO 100 J = 1, N
                        CALL DSCAL( N, SCALOC, C( 1, J ), 1 )
  100                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 1, 2 )
                  C( K2, L1 ) = X( 2, 1 )
                  C( K2, L2 ) = X( 2, 2 )
                  IF( K1.NE.L1 ) THEN
                     C( L1, K1 ) = X( 1, 1 )
                     C( L2, K1 ) = X( 1, 2 )
                     C( L1, K2 ) = X( 2, 1 )
                     C( L2, K2 ) = X( 2, 2 )
                  END IF
               END IF
  110       CONTINUE
  120    CONTINUE
      END IF
      RETURN
      END
      SUBROUTINE SB03MV( LTRAN, LUPPER, T, LDT, B, LDB, SCALE, X, LDX,
     $                   XNORM, INFO )
      DOUBLE PRECISION   ZERO, ONE, FOUR
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, FOUR = 4.0D+0 )
      LOGICAL            LTRAN, LUPPER
      INTEGER            INFO, LDB, LDT, LDX
      DOUBLE PRECISION   SCALE, XNORM
      DOUBLE PRECISION   B( LDB, * ), T( LDT, * ), X( LDX, * )
      INTEGER            I, IP, IPSV, J, JP, JPSV, K
      DOUBLE PRECISION   EPS, SMIN, SMLNUM, TEMP, XMAX
      INTEGER            JPIV( 3 )
      DOUBLE PRECISION   BTMP( 3 ), T9( 3, 3 ), TMP( 3 )
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
      EXTERNAL           DSWAP
      INTRINSIC          ABS, MAX
      INFO = 0
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' ) / EPS
      SMIN = MAX( ABS( T( 1, 1 ) ), ABS( T( 1, 2 ) ),
     $            ABS( T( 2, 1 ) ), ABS( T( 2, 2 ) ) )
      SMIN = MAX( EPS*SMIN, SMLNUM )
      T9( 1, 1 ) = T( 1, 1 )*T( 1, 1 ) - ONE
      T9( 2, 2 ) = T( 1, 1 )*T( 2, 2 ) + T( 1, 2 )*T( 2, 1 ) - ONE
      T9( 3, 3 ) = T( 2, 2 )*T( 2, 2 ) - ONE
      IF( LTRAN ) THEN
         T9( 1, 2 ) = T( 1, 1 )*T( 1, 2 ) + T( 1, 1 )*T( 1, 2 )
         T9( 1, 3 ) = T( 1, 2 )*T( 1, 2 )
         T9( 2, 1 ) = T( 1, 1 )*T( 2, 1 )
         T9( 2, 3 ) = T( 1, 2 )*T( 2, 2 )
         T9( 3, 1 ) = T( 2, 1 )*T( 2, 1 )
         T9( 3, 2 ) = T( 2, 1 )*T( 2, 2 ) + T( 2, 1 )*T( 2, 2 )
      ELSE
         T9( 1, 2 ) = T( 1, 1 )*T( 2, 1 ) + T( 1, 1 )*T( 2, 1 )
         T9( 1, 3 ) = T( 2, 1 )*T( 2, 1 )
         T9( 2, 1 ) = T( 1, 1 )*T( 1, 2 )
         T9( 2, 3 ) = T( 2, 1 )*T( 2, 2 )
         T9( 3, 1 ) = T( 1, 2 )*T( 1, 2 )
         T9( 3, 2 ) = T( 1, 2 )*T( 2, 2 ) + T( 1, 2 )*T( 2, 2 )
      END IF
      BTMP( 1 ) = B( 1, 1 )
      IF ( LUPPER ) THEN
         BTMP( 2 ) = B( 1, 2 )
      ELSE
         BTMP( 2 ) = B( 2, 1 )
      END IF
      BTMP( 3 ) = B( 2, 2 )
      DO 50 I = 1, 2
         XMAX = ZERO
         DO 20 IP = I, 3
            DO 10 JP = I, 3
               IF( ABS( T9( IP, JP ) ).GE.XMAX ) THEN
                  XMAX = ABS( T9( IP, JP ) )
                  IPSV = IP
                  JPSV = JP
               END IF
   10       CONTINUE
   20    CONTINUE
         IF( IPSV.NE.I ) THEN
            CALL DSWAP( 3, T9( IPSV, 1 ), 3, T9( I, 1 ), 3 )
            TEMP = BTMP( I )
            BTMP( I ) = BTMP( IPSV )
            BTMP( IPSV ) = TEMP
         END IF
         IF( JPSV.NE.I )
     $      CALL DSWAP( 3, T9( 1, JPSV ), 1, T9( 1, I ), 1 )
         JPIV( I ) = JPSV
         IF( ABS( T9( I, I ) ).LT.SMIN ) THEN
            INFO = 1
            T9( I, I ) = SMIN
         END IF
         DO 40 J = I + 1, 3
            T9( J, I ) = T9( J, I ) / T9( I, I )
            BTMP( J ) = BTMP( J ) - T9( J, I )*BTMP( I )
            DO 30 K = I + 1, 3
               T9( J, K ) = T9( J, K ) - T9( J, I )*T9( I, K )
   30       CONTINUE
   40    CONTINUE
   50 CONTINUE
      IF( ABS( T9( 3, 3 ) ).LT.SMIN )
     $   T9( 3, 3 ) = SMIN
      SCALE = ONE
      IF( ( FOUR*SMLNUM )*ABS( BTMP( 1 ) ).GT.ABS( T9( 1, 1 ) ) .OR.
     $    ( FOUR*SMLNUM )*ABS( BTMP( 2 ) ).GT.ABS( T9( 2, 2 ) ) .OR.
     $    ( FOUR*SMLNUM )*ABS( BTMP( 3 ) ).GT.ABS( T9( 3, 3 ) ) ) THEN
         SCALE = ( ONE / FOUR ) / MAX( ABS( BTMP( 1 ) ),
     $               ABS( BTMP( 2 ) ), ABS( BTMP( 3 ) ) )
         BTMP( 1 ) = BTMP( 1 )*SCALE
         BTMP( 2 ) = BTMP( 2 )*SCALE
         BTMP( 3 ) = BTMP( 3 )*SCALE
      END IF
      DO 70 I = 1, 3
         K = 4 - I
         TEMP = ONE / T9( K, K )
         TMP( K ) = BTMP( K )*TEMP
         DO 60 J = K + 1, 3
            TMP( K ) = TMP( K ) - ( TEMP*T9( K, J ) )*TMP( J )
  60     CONTINUE
  70  CONTINUE
      DO 80 I = 1, 2
         IF( JPIV( 3-I ).NE.3-I ) THEN
            TEMP = TMP( 3-I )
            TMP( 3-I ) = TMP( JPIV( 3-I ) )
            TMP( JPIV( 3-I ) ) = TEMP
         END IF
  80  CONTINUE
      X( 1, 1 ) = TMP( 1 )
      IF ( LUPPER ) THEN
         X( 1, 2 ) = TMP( 2 )
      ELSE
         X( 2, 1 ) = TMP( 2 )
      END IF
      X( 2, 2 ) = TMP( 3 )
      XNORM = MAX( ABS( TMP( 1 ) ) + ABS( TMP( 2 ) ),
     $             ABS( TMP( 2 ) ) + ABS( TMP( 3 ) ) )
      RETURN
      END
      SUBROUTINE SB03MW( LTRAN, LUPPER, T, LDT, B, LDB, SCALE, X, LDX,
     $                   XNORM, INFO )
      DOUBLE PRECISION   ZERO, ONE, TWO, FOUR
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0,
     $                     FOUR = 4.0D+0 )
      LOGICAL            LTRAN, LUPPER
      INTEGER            INFO, LDB, LDT, LDX
      DOUBLE PRECISION   SCALE, XNORM
      DOUBLE PRECISION   B( LDB, * ), T( LDT, * ), X( LDX, * )
      INTEGER            I, IP, IPSV, J, JP, JPSV, K
      DOUBLE PRECISION   EPS, SMIN, SMLNUM, TEMP, XMAX
      INTEGER            JPIV( 3 )
      DOUBLE PRECISION   BTMP( 3 ), T9( 3, 3 ), TMP( 3 )
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
      EXTERNAL           DSWAP
      INTRINSIC          ABS, MAX
      INFO = 0
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' ) / EPS
      SMIN = MAX( MAX( ABS( T( 1, 1 ) ), ABS( T( 1, 2 ) ),
     $                 ABS( T( 2, 1 ) ), ABS( T( 2, 2 ) ) )*EPS,
     $            SMLNUM )
      T9( 1, 3 ) = ZERO
      T9( 3, 1 ) = ZERO
      T9( 1, 1 ) = T( 1, 1 )
      T9( 2, 2 ) = T( 1, 1 ) + T( 2, 2 )
      T9( 3, 3 ) = T( 2, 2 )
      IF( LTRAN ) THEN
         T9( 1, 2 ) = T( 1, 2 )
         T9( 2, 1 ) = T( 2, 1 )
         T9( 2, 3 ) = T( 1, 2 )
         T9( 3, 2 ) = T( 2, 1 )
      ELSE
         T9( 1, 2 ) = T( 2, 1 )
         T9( 2, 1 ) = T( 1, 2 )
         T9( 2, 3 ) = T( 2, 1 )
         T9( 3, 2 ) = T( 1, 2 )
      END IF
      BTMP( 1 ) = B( 1, 1 )/TWO
      IF ( LUPPER ) THEN
         BTMP( 2 ) = B( 1, 2 )
      ELSE
         BTMP( 2 ) = B( 2, 1 )
      END IF
      BTMP( 3 ) = B( 2, 2 )/TWO
      DO 50 I = 1, 2
         XMAX = ZERO
         DO 20 IP = I, 3
            DO 10 JP = I, 3
               IF( ABS( T9( IP, JP ) ).GE.XMAX ) THEN
                  XMAX = ABS( T9( IP, JP ) )
                  IPSV = IP
                  JPSV = JP
               END IF
   10       CONTINUE
   20    CONTINUE
         IF( IPSV.NE.I ) THEN
            CALL DSWAP( 3, T9( IPSV, 1 ), 3, T9( I, 1 ), 3 )
            TEMP = BTMP( I )
            BTMP( I ) = BTMP( IPSV )
            BTMP( IPSV ) = TEMP
         END IF
         IF( JPSV.NE.I )
     $      CALL DSWAP( 3, T9( 1, JPSV ), 1, T9( 1, I ), 1 )
         JPIV( I ) = JPSV
         IF( ABS( T9( I, I ) ).LT.SMIN ) THEN
            INFO = 1
            T9( I, I ) = SMIN
         END IF
         DO 40 J = I + 1, 3
            T9( J, I ) = T9( J, I ) / T9( I, I )
            BTMP( J ) = BTMP( J ) - T9( J, I )*BTMP( I )
            DO 30 K = I + 1, 3
               T9( J, K ) = T9( J, K ) - T9( J, I )*T9( I, K )
   30       CONTINUE
   40    CONTINUE
   50 CONTINUE
      IF( ABS( T9( 3, 3 ) ).LT.SMIN )
     $   T9( 3, 3 ) = SMIN
      SCALE = ONE
      IF( ( FOUR*SMLNUM )*ABS( BTMP( 1 ) ).GT.ABS( T9( 1, 1 ) ) .OR.
     $    ( FOUR*SMLNUM )*ABS( BTMP( 2 ) ).GT.ABS( T9( 2, 2 ) ) .OR.
     $    ( FOUR*SMLNUM )*ABS( BTMP( 3 ) ).GT.ABS( T9( 3, 3 ) ) ) THEN
         SCALE = ( ONE / FOUR ) / MAX( ABS( BTMP( 1 ) ),
     $               ABS( BTMP( 2 ) ), ABS( BTMP( 3 ) ) )
         BTMP( 1 ) = BTMP( 1 )*SCALE
         BTMP( 2 ) = BTMP( 2 )*SCALE
         BTMP( 3 ) = BTMP( 3 )*SCALE
      END IF
      DO 70 I = 1, 3
         K = 4 - I
         TEMP = ONE / T9( K, K )
         TMP( K ) = BTMP( K )*TEMP
         DO 60 J = K + 1, 3
            TMP( K ) = TMP( K ) - ( TEMP*T9( K, J ) )*TMP( J )
  60     CONTINUE
  70  CONTINUE
      DO 80 I = 1, 2
         IF( JPIV( 3-I ).NE.3-I ) THEN
            TEMP = TMP( 3-I )
            TMP( 3-I ) = TMP( JPIV( 3-I ) )
            TMP( JPIV( 3-I ) ) = TEMP
         END IF
  80  CONTINUE
      X( 1, 1 ) = TMP( 1 )
      IF ( LUPPER ) THEN
         X( 1, 2 ) = TMP( 2 )
      ELSE
         X( 2, 1 ) = TMP( 2 )
      END IF
      X( 2, 2 ) = TMP( 3 )
      XNORM = MAX( ABS( TMP( 1 ) ) + ABS( TMP( 2 ) ),
     $             ABS( TMP( 2 ) ) + ABS( TMP( 3 ) ) )
      RETURN
      END
      SUBROUTINE SB04PX( LTRANL, LTRANR, ISGN, N1, N2, TL, LDTL, TR,
     $                   LDTR, B, LDB, SCALE, X, LDX, XNORM, INFO )
      DOUBLE PRECISION   ZERO, ONE, TWO, HALF, EIGHT
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0,
     $                     TWO = 2.0D+0, HALF = 0.5D+0, EIGHT = 8.0D+0 )
      LOGICAL            LTRANL, LTRANR
      INTEGER            INFO, ISGN, LDB, LDTL, LDTR, LDX, N1, N2
      DOUBLE PRECISION   SCALE, XNORM
      DOUBLE PRECISION   B( LDB, * ), TL( LDTL, * ), TR( LDTR, * ),
     $                   X( LDX, * )
      LOGICAL            BSWAP, XSWAP
      INTEGER            I, IP, IPIV, IPSV, J, JP, JPSV, K
      DOUBLE PRECISION   BET, EPS, GAM, L21, SGN, SMIN, SMLNUM, TAU1,
     $                   TEMP, U11, U12, U22, XMAX
      LOGICAL            BSWPIV( 4 ), XSWPIV( 4 )
      INTEGER            JPIV( 4 ), LOCL21( 4 ), LOCU12( 4 ),
     $                   LOCU22( 4 )
      DOUBLE PRECISION   BTMP( 4 ), T16( 4, 4 ), TMP( 4 ), X2( 2 )
      INTEGER            IDAMAX
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH, IDAMAX
      EXTERNAL           DSWAP
      INTRINSIC          ABS, MAX
      DATA               LOCU12 / 3, 4, 1, 2 / , LOCL21 / 2, 1, 4, 3 / ,
     $                   LOCU22 / 4, 3, 2, 1 /
      DATA               XSWPIV / .FALSE., .FALSE., .TRUE., .TRUE. /
      DATA               BSWPIV / .FALSE., .TRUE., .FALSE., .TRUE. /
      INFO  = 0
      SCALE = ONE
      IF( N1.EQ.0 .OR. N2.EQ.0 ) THEN
         XNORM = ZERO
         RETURN
      END IF
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' ) / EPS
      SGN = ISGN
      K = N1 + N1 + N2 - 2
      GO TO ( 10, 20, 30, 50 )K
   10 CONTINUE
      TAU1 = TL( 1, 1 )*TR( 1, 1 ) + SGN
      BET  = ABS( TAU1 )
      IF( BET.LE.SMLNUM ) THEN
         TAU1 = SMLNUM
         BET  = SMLNUM
         INFO = 1
      END IF
      GAM = ABS( B( 1, 1 ) )
      IF( SMLNUM*GAM.GT.BET )
     $   SCALE = ONE / GAM
      X( 1, 1 ) = ( B( 1, 1 )*SCALE ) / TAU1
      XNORM = ABS( X( 1, 1 ) )
      RETURN
   20 CONTINUE
      SMIN = MAX( MAX( ABS( TR( 1, 1 ) ), ABS( TR( 1, 2 ) ),
     $                 ABS( TR( 2, 1 ) ), ABS( TR( 2, 2 ) ) )
     $                *ABS( TL( 1, 1 ) )*EPS,
     $            SMLNUM )
      TMP( 1 ) = TL( 1, 1 )*TR( 1, 1 ) + SGN
      TMP( 4 ) = TL( 1, 1 )*TR( 2, 2 ) + SGN
      IF( LTRANR ) THEN
         TMP( 2 ) = TL( 1, 1 )*TR( 2, 1 )
         TMP( 3 ) = TL( 1, 1 )*TR( 1, 2 )
      ELSE
         TMP( 2 ) = TL( 1, 1 )*TR( 1, 2 )
         TMP( 3 ) = TL( 1, 1 )*TR( 2, 1 )
      END IF
      BTMP( 1 ) = B( 1, 1 )
      BTMP( 2 ) = B( 1, 2 )
      GO TO 40
   30 CONTINUE
      SMIN = MAX( MAX( ABS( TL( 1, 1 ) ), ABS( TL( 1, 2 ) ),
     $                 ABS( TL( 2, 1 ) ), ABS( TL( 2, 2 ) ) )
     $                *ABS( TR( 1, 1 ) )*EPS,
     $            SMLNUM )
      TMP( 1 ) = TL( 1, 1 )*TR( 1, 1 ) + SGN
      TMP( 4 ) = TL( 2, 2 )*TR( 1, 1 ) + SGN
      IF( LTRANL ) THEN
         TMP( 2 ) = TL( 1, 2 )*TR( 1, 1 )
         TMP( 3 ) = TL( 2, 1 )*TR( 1, 1 )
      ELSE
         TMP( 2 ) = TL( 2, 1 )*TR( 1, 1 )
         TMP( 3 ) = TL( 1, 2 )*TR( 1, 1 )
      END IF
      BTMP( 1 ) = B( 1, 1 )
      BTMP( 2 ) = B( 2, 1 )
   40 CONTINUE
      IPIV = IDAMAX( 4, TMP, 1 )
      U11  = TMP( IPIV )
      IF( ABS( U11 ).LE.SMIN ) THEN
         INFO = 1
         U11  = SMIN
      END IF
      U12 = TMP( LOCU12( IPIV ) )
      L21 = TMP( LOCL21( IPIV ) ) / U11
      U22 = TMP( LOCU22( IPIV ) ) - U12*L21
      XSWAP = XSWPIV( IPIV )
      BSWAP = BSWPIV( IPIV )
      IF( ABS( U22 ).LE.SMIN ) THEN
         INFO = 1
         U22  = SMIN
      END IF
      IF( BSWAP ) THEN
         TEMP = BTMP( 2 )
         BTMP( 2 ) = BTMP( 1 ) - L21*TEMP
         BTMP( 1 ) = TEMP
      ELSE
         BTMP( 2 ) = BTMP( 2 ) - L21*BTMP( 1 )
      END IF
      IF( ( TWO*SMLNUM )*ABS( BTMP( 2 ) ).GT.ABS( U22 ) .OR.
     $    ( TWO*SMLNUM )*ABS( BTMP( 1 ) ).GT.ABS( U11 ) ) THEN
         SCALE = HALF / MAX( ABS( BTMP( 1 ) ), ABS( BTMP( 2 ) ) )
         BTMP( 1 ) = BTMP( 1 )*SCALE
         BTMP( 2 ) = BTMP( 2 )*SCALE
      END IF
      X2( 2 ) = BTMP( 2 ) / U22
      X2( 1 ) = BTMP( 1 ) / U11 - ( U12 / U11 )*X2( 2 )
      IF( XSWAP ) THEN
         TEMP = X2( 2 )
         X2( 2 ) = X2( 1 )
         X2( 1 ) = TEMP
      END IF
      X( 1, 1 ) = X2( 1 )
      IF( N1.EQ.1 ) THEN
         X( 1, 2 ) = X2( 2 )
         XNORM = ABS( X2( 1 ) ) + ABS( X2( 2 ) )
      ELSE
         X( 2, 1 ) = X2( 2 )
         XNORM = MAX( ABS( X2( 1 ) ), ABS( X2( 2 ) ) )
      END IF
      RETURN
   50 CONTINUE
      SMIN = MAX( ABS( TR( 1, 1 ) ), ABS( TR( 1, 2 ) ),
     $            ABS( TR( 2, 1 ) ), ABS( TR( 2, 2 ) ) )
      SMIN = MAX( ABS( TL( 1, 1 ) ), ABS( TL( 1, 2 ) ),
     $            ABS( TL( 2, 1 ) ), ABS( TL( 2, 2 ) ) )*SMIN
      SMIN = MAX( EPS*SMIN, SMLNUM )
      T16( 1, 1 ) = TL( 1, 1 )*TR( 1, 1 ) + SGN
      T16( 2, 2 ) = TL( 2, 2 )*TR( 1, 1 ) + SGN
      T16( 3, 3 ) = TL( 1, 1 )*TR( 2, 2 ) + SGN
      T16( 4, 4 ) = TL( 2, 2 )*TR( 2, 2 ) + SGN
      IF( LTRANL ) THEN
         T16( 1, 2 ) = TL( 2, 1 )*TR( 1, 1 )
         T16( 2, 1 ) = TL( 1, 2 )*TR( 1, 1 )
         T16( 3, 4 ) = TL( 2, 1 )*TR( 2, 2 )
         T16( 4, 3 ) = TL( 1, 2 )*TR( 2, 2 )
      ELSE
         T16( 1, 2 ) = TL( 1, 2 )*TR( 1, 1 )
         T16( 2, 1 ) = TL( 2, 1 )*TR( 1, 1 )
         T16( 3, 4 ) = TL( 1, 2 )*TR( 2, 2 )
         T16( 4, 3 ) = TL( 2, 1 )*TR( 2, 2 )
      END IF
      IF( LTRANR ) THEN
         T16( 1, 3 ) = TL( 1, 1 )*TR( 1, 2 )
         T16( 2, 4 ) = TL( 2, 2 )*TR( 1, 2 )
         T16( 3, 1 ) = TL( 1, 1 )*TR( 2, 1 )
         T16( 4, 2 ) = TL( 2, 2 )*TR( 2, 1 )
      ELSE
         T16( 1, 3 ) = TL( 1, 1 )*TR( 2, 1 )
         T16( 2, 4 ) = TL( 2, 2 )*TR( 2, 1 )
         T16( 3, 1 ) = TL( 1, 1 )*TR( 1, 2 )
         T16( 4, 2 ) = TL( 2, 2 )*TR( 1, 2 )
      END IF
      IF( LTRANL .AND. LTRANR ) THEN
         T16( 1, 4 ) = TL( 2, 1 )*TR( 1, 2 )
         T16( 2, 3 ) = TL( 1, 2 )*TR( 1, 2 )
         T16( 3, 2 ) = TL( 2, 1 )*TR( 2, 1 )
         T16( 4, 1 ) = TL( 1, 2 )*TR( 2, 1 )
      ELSE IF( LTRANL .AND. .NOT.LTRANR ) THEN
         T16( 1, 4 ) = TL( 2, 1 )*TR( 2, 1 )
         T16( 2, 3 ) = TL( 1, 2 )*TR( 2, 1 )
         T16( 3, 2 ) = TL( 2, 1 )*TR( 1, 2 )
         T16( 4, 1 ) = TL( 1, 2 )*TR( 1, 2 )
      ELSE IF( .NOT.LTRANL .AND. LTRANR ) THEN
          T16( 1, 4 ) = TL( 1, 2 )*TR( 1, 2 )
          T16( 2, 3 ) = TL( 2, 1 )*TR( 1, 2 )
          T16( 3, 2 ) = TL( 1, 2 )*TR( 2, 1 )
          T16( 4, 1 ) = TL( 2, 1 )*TR( 2, 1 )
      ELSE
          T16( 1, 4 ) = TL( 1, 2 )*TR( 2, 1 )
          T16( 2, 3 ) = TL( 2, 1 )*TR( 2, 1 )
          T16( 3, 2 ) = TL( 1, 2 )*TR( 1, 2 )
          T16( 4, 1 ) = TL( 2, 1 )*TR( 1, 2 )
      END IF
      BTMP( 1 ) = B( 1, 1 )
      BTMP( 2 ) = B( 2, 1 )
      BTMP( 3 ) = B( 1, 2 )
      BTMP( 4 ) = B( 2, 2 )
      DO 100 I = 1, 3
         XMAX = ZERO
         DO 70 IP = I, 4
            DO 60 JP = I, 4
               IF( ABS( T16( IP, JP ) ).GE.XMAX ) THEN
                  XMAX = ABS( T16( IP, JP ) )
                  IPSV = IP
                  JPSV = JP
               END IF
   60       CONTINUE
   70    CONTINUE
         IF( IPSV.NE.I ) THEN
            CALL DSWAP( 4, T16( IPSV, 1 ), 4, T16( I, 1 ), 4 )
            TEMP = BTMP( I )
            BTMP( I ) = BTMP( IPSV )
            BTMP( IPSV ) = TEMP
         END IF
         IF( JPSV.NE.I )
     $      CALL DSWAP( 4, T16( 1, JPSV ), 1, T16( 1, I ), 1 )
         JPIV( I ) = JPSV
         IF( ABS( T16( I, I ) ).LT.SMIN ) THEN
            INFO = 1
            T16( I, I ) = SMIN
         END IF
         DO 90 J = I + 1, 4
            T16( J, I ) = T16( J, I ) / T16( I, I )
            BTMP( J ) = BTMP( J ) - T16( J, I )*BTMP( I )
            DO 80 K = I + 1, 4
               T16( J, K ) = T16( J, K ) - T16( J, I )*T16( I, K )
   80       CONTINUE
   90    CONTINUE
  100 CONTINUE
      IF( ABS( T16( 4, 4 ) ).LT.SMIN )
     $   T16( 4, 4 ) = SMIN
      IF( ( EIGHT*SMLNUM )*ABS( BTMP( 1 ) ).GT.ABS( T16( 1, 1 ) ) .OR.
     $    ( EIGHT*SMLNUM )*ABS( BTMP( 2 ) ).GT.ABS( T16( 2, 2 ) ) .OR.
     $    ( EIGHT*SMLNUM )*ABS( BTMP( 3 ) ).GT.ABS( T16( 3, 3 ) ) .OR.
     $    ( EIGHT*SMLNUM )*ABS( BTMP( 4 ) ).GT.ABS( T16( 4, 4 ) ) ) THEN
         SCALE = ( ONE / EIGHT ) / MAX( ABS( BTMP( 1 ) ),
     $                ABS( BTMP( 2 ) ), ABS( BTMP( 3 ) ),
     $                ABS( BTMP( 4 ) ) )
         BTMP( 1 ) = BTMP( 1 )*SCALE
         BTMP( 2 ) = BTMP( 2 )*SCALE
         BTMP( 3 ) = BTMP( 3 )*SCALE
         BTMP( 4 ) = BTMP( 4 )*SCALE
      END IF
      DO 120 I = 1, 4
         K = 5 - I
         TEMP = ONE / T16( K, K )
         TMP( K ) = BTMP( K )*TEMP
         DO 110 J = K + 1, 4
            TMP( K ) = TMP( K ) - ( TEMP*T16( K, J ) )*TMP( J )
  110    CONTINUE
  120 CONTINUE
      DO 130 I = 1, 3
         IF( JPIV( 4-I ).NE.4-I ) THEN
            TEMP = TMP( 4-I )
            TMP( 4-I ) = TMP( JPIV( 4-I ) )
            TMP( JPIV( 4-I ) ) = TEMP
         END IF
  130 CONTINUE
      X( 1, 1 ) = TMP( 1 )
      X( 2, 1 ) = TMP( 2 )
      X( 1, 2 ) = TMP( 3 )
      X( 2, 2 ) = TMP( 4 )
      XNORM = MAX( ABS( TMP( 1 ) ) + ABS( TMP( 3 ) ),
     $             ABS( TMP( 2 ) ) + ABS( TMP( 4 ) ) )
      RETURN
      END
      LOGICAL FUNCTION  SELECT( PAR1, PAR2 )
      DOUBLE PRECISION  PAR1, PAR2
      SELECT = .TRUE.
      RETURN
      END
