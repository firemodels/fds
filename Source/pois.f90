MODULE POIS

! CODE CONVERTED USING TO_F90 BY ALAN MILLER
! DATE: 2006-01-18  TIME: 13:48:46

! POISSON SOLVER ROUTINES

USE PRECISION_PARAMETERS
IMPLICIT NONE
PRIVATE
REAL(EB) SCALE
INTEGER :: KAPPA,NMAX,IKPWR
LOGICAL :: TPOSE,NOCOPY,OUTARY

PUBLIC H3CZIS,H3CZSS,H2CZSS,H2CYSS,H3CSSS,H2CZIS,H3CSIS,H2CYIS

CONTAINS

SUBROUTINE H3CZIS(XS,XF,L,LBDCND,YS,YF,M,MBDCND,ZS,ZF,N,NBDCND,  &
    H,ELMBDA,LDIMF,MDIMF,IERROR,SAVE)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

REAL(EB)::  XS, XF, YS, YF, ZS, ZF, ELMBDA, DLZSQR
INTEGER:: L, LBDCND, M, MBDCND, N, NBDCND, LDIMF, MDIMF, IERROR, LPEROD
REAL(EB) SAVE(-3:*),H(0:*)

REAL(EB) :: DX, DY, DZ, DLXSQR, DLYSQR, HM, HP
INTEGER :: LP, MP, NP, IA, IB, IC, ID, IS, I, IR

!                               CHECK FOR INVALID INPUT

IERROR = 0

IF (XS>XF) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 1._EB
END IF

IF (L<3) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 2._EB
END IF

IF (LBDCND<0 .OR. LBDCND>4) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 3._EB
END IF

IF (YS>YF) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 4._EB
END IF

IF (M<3) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 5._EB
END IF

IF (MBDCND<0 .OR. MBDCND>4) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 6._EB
END IF

IF (ZS>ZF) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 7._EB
END IF

IF (N<3) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 8._EB
END IF

IF (NBDCND<0 .OR. NBDCND>4) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 9._EB
END IF

IF (LDIMF<L) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 10._EB
END IF

IF (MDIMF<M) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 11._EB
END IF

IF (IERROR/=0) THEN
  RETURN
ELSE
  SAVE(1) = 0._EB
END IF

!                               DEFINE GRID PARAMETERS

DX = (XF-XS)/L
DY = (YF-YS)/M
DZ = (ZF-ZS)/N
DLXSQR = 1._EB/DX**2
DLYSQR = 1._EB/DY**2
DLZSQR = 1._EB/DZ**2
LP = LBDCND + 1
MP = MBDCND + 1
NP = NBDCND + 1

!                               ALLOCATE SAVE ARRAY

IA = 12
IB = IA + L
IC = IB + L
ID = IC + L
IS = ID + L

!                               DEFINE THE A,B,C COEFFICIENTS
!                               IN SAVE ARRAY

DO  I = 0,L - 1
  HM = .5_EB*(H(I)+H(I+1))
  HP = .5_EB*(H(I+1)+H(I+2))
  SAVE(IA+I) = DLXSQR/(H(I+1)*HM)
  SAVE(IC+I) = DLXSQR/(H(I+1)*HP)
  SAVE(IB+I) = - (SAVE(IA+I)+SAVE(IC+I)) + ELMBDA
  SAVE(ID+I) = DLYSQR
END DO

SELECT CASE(LP)
CASE(2)
SAVE(IB) = SAVE(IB) - SAVE(IA)
SAVE(IC-1) = SAVE(IC-1) - SAVE(ID-1)
CASE(3)
SAVE(IB) = SAVE(IB) - SAVE(IA)
SAVE(IC-1) = SAVE(IC-1) + SAVE(ID-1)
CASE(4)
SAVE(IB) = SAVE(IB) + SAVE(IA)
SAVE(IC-1) = SAVE(IC-1) + SAVE(ID-1)
CASE(5)
SAVE(IB) = SAVE(IB) + SAVE(IA)
SAVE(IC-1) = SAVE(IC-1) - SAVE(ID-1)
END SELECT


!                               DETERMINE WHETHER OR NOT BOUNDARY
!                               CONDITION IS PERIODIC IN X

IF (LBDCND==0) THEN
  LPEROD = 0
ELSE
  LPEROD = 1
END IF

!                               INITIALIZE SOLVER ROUTINE S3CFIS.
CALL S3CFIS(LPEROD,L,MBDCND,M,NBDCND,N,DLZSQR,SAVE(IA),SAVE(IB),  &
    SAVE(IC),SAVE(ID),LDIMF,MDIMF,IR,SAVE(IS))

!                               TEST ERROR FLAG FROM S3CFIS FOR
!                               INTERNAL ERROR

IF (IR/=0) THEN
  SAVE(1) = 99._EB
  IERROR = 1

  RETURN
END IF

!                               SAVE PARAMETERS FOR H3CZSS IN SAVE ARRAY

SAVE(2) = DX
SAVE(3) = REAL(L,EB)
SAVE(4) = REAL(LP,EB)
SAVE(5) = DY
SAVE(6) = REAL(M,EB)
SAVE(7) = REAL(MP,EB)
SAVE(8) = DZ
SAVE(9) = REAL(N,EB)
SAVE(10) = REAL(NP,EB)
SAVE(11) = ELMBDA

RETURN
END SUBROUTINE H3CZIS

SUBROUTINE H3CZSS(BDXS,BDXF,BDYS,BDYF,BDZS,BDZF,LDIMF,MDIMF,F,PERTRB,SAVE,W,H)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+


USE GLOBAL_CONSTANTS, ONLY: PRES_METHOD
INTEGER:: LDIMF, MDIMF
INTEGER :: L, LP, M, MP, N, NP, IA, IB, IC, ID, IS, K, J, I, ISING
REAL(EB):: BDXS(MDIMF,*), BDXF(MDIMF,*), BDYS(LDIMF,*), BDYF(LDIMF,*), BDZS(LDIMF,*), BDZF(LDIMF,*), &
           F(LDIMF,MDIMF,*),SAVE(-3:*),W(*),H(0:*), PERTRB
REAL(EB) :: DX, DY, DZ, ELMBDA, DLYRCP, TWDYSQ, DLZRCP, TWDZSQ, PERT, S1, S3, PRTSAV


!                               CHECK VALUE OF IERROR (=SAVE(1)).
!                               IF NON-ZERO, RETURN.

IF (ABS(SAVE(1))>=TWO_EPSILON_EB) RETURN

!                               GET PARAMETERS FOR H3CZSS FROM SAVE
!                               ARRAY WHERE THEY WERE STORED IN
!                               INITIALIZATION SUBROUTINE H3CZIS.

DX = SAVE(2)
L  = NINT(SAVE(3))
LP = NINT(SAVE(4))
DY = SAVE(5)
M  = NINT(SAVE(6))
MP = NINT(SAVE(7))
DZ = SAVE(8)
N  = NINT(SAVE(9))
NP = NINT(SAVE(10))
ELMBDA = SAVE(11)

DLYRCP = 1._EB/DY
TWDYSQ = 2._EB/DY**2
DLZRCP = 1._EB/DZ
TWDZSQ = 2._EB/DZ**2

!                               ALLOCATE SAVE ARRAY

IA = 12
IB = IA + L
IC = IB + L
ID = IC + L
IS = 12 + 4*L


IF (PRES_METHOD /= 'SCARC' .AND. PRES_METHOD /= 'USCARC') THEN

!                               ENTER BOUNDARY DATA FOR X-BOUNDARIES

   IF (LP==2 .OR. LP==3) THEN
     DO K = 1,N
       DO J = 1,M
         F(1,J,K) = F(1,J,K) - 2._EB*BDXS(J,K)*SAVE(IA)
       END DO
     END DO
   END IF

   IF (LP==4 .OR. LP==5) THEN
     DO K = 1,N
       DO J = 1,M
         F(1,J,K) = F(1,J,K) + SAVE(IA)*DX*BDXS(J,K)
       END DO
     END DO
   END IF

   IF (LP==2 .OR. LP==5) THEN
     DO K = 1,N
       DO J = 1,M
         F(L,J,K) = F(L,J,K) - 2._EB*BDXF(J,K)*SAVE(ID-1)
       END DO
     END DO
   END IF

   IF (LP==3 .OR. LP==4) THEN
     DO K = 1,N
       DO J = 1,M
         F(L,J,K) = F(L,J,K) - SAVE(ID-1)*DX*BDXF(J,K)
       END DO
     END DO
   END IF

   !                               ENTER BOUNDARY DATA FOR Y-BOUNDARIES

   IF (MP==2 .OR. MP==3) THEN
     DO K = 1,N
       DO I = 1,L
         F(I,1,K) = F(I,1,K) - BDYS(I,K)*TWDYSQ
       END DO
     END DO
   END IF

   IF (MP==4 .OR. MP==5) THEN
     DO K = 1,N
       DO I = 1,L
         F(I,1,K) = F(I,1,K) + BDYS(I,K)*DLYRCP
       END DO
     END DO
   END IF

   IF (MP==2 .OR. MP==5) THEN
     DO K = 1,N
       DO I = 1,L
         F(I,M,K) = F(I,M,K) - BDYF(I,K)*TWDYSQ
       END DO
     END DO
   END IF

   IF (MP==3 .OR. MP==4) THEN
     DO K = 1,N
       DO I = 1,L
         F(I,M,K) = F(I,M,K) - BDYF(I,K)*DLYRCP
       END DO
     END DO
   END IF

   !                               ENTER BOUNDARY DATA FOR Z-BOUNDARIES

   IF (NP==2 .OR. NP==3) THEN
     DO J = 1,M
       DO I = 1,L
         F(I,J,1) = F(I,J,1) - BDZS(I,J)*TWDZSQ
       END DO
     END DO
   END IF

   IF (NP==4 .OR. NP==5) THEN
     DO J = 1,M
       DO I = 1,L
         F(I,J,1) = F(I,J,1) + BDZS(I,J)*DLZRCP
       END DO
     END DO
   END IF

   IF (NP==2 .OR. NP==5) THEN
     DO J = 1,M
       DO I = 1,L
         F(I,J,N) = F(I,J,N) - BDZF(I,J)*TWDZSQ
       END DO
     END DO
   END IF

   IF (NP==3 .OR. NP==4) THEN
     DO J = 1,M
       DO I = 1,L
         F(I,J,N) = F(I,J,N) - BDZF(I,J)*DLZRCP
       END DO
     END DO
   END IF

END IF

PERTRB = 0._EB
PERT   = 0._EB
ISING = 0

!                               FOR SINGULAR PROBLEMS ADJUST DATA TO
!                               INSURE A SOLUTION WILL EXIST.  GO THRU
!                               THIS CODE TWICE: ISING=1 FOR CALCULATING
!                               PERTRB; ISING=2 FOR NORMALIZING SOLUTION
!                               AFTER IT IS COMPUTED.

SELECT CASE(LP)
CASE(1)   ; GO TO 630
CASE(2:3) ; GO TO 750
CASE(4)   ; GO TO 630
CASE(5)   ; GO TO 750
END SELECT

630 CONTINUE
SELECT CASE(MP)
CASE(1)   ; GO TO 640
CASE(2:3) ; GO TO 750
CASE(4)   ; GO TO 640
CASE(5)   ; GO TO 750
END SELECT

640 CONTINUE
SELECT CASE(NP)
CASE(1)   ; GO TO 650
CASE(2:3) ; GO TO 750
CASE(4)   ; GO TO 650
CASE(5)   ; GO TO 750
END SELECT

650 CONTINUE
IF (ABS(ELMBDA)>=TWO_EPSILON_EB) GO TO 750
ISING = 1
660 CONTINUE
PERT = 0._EB
DO  I = 1,L
  W(I) = 0._EB
END DO
DO  K = 1,N
  DO  J = 1,M
    DO  I = 1,L
      W(I) = W(I) + F(I,J,K)
    END DO
  END DO
END DO
S1 = 0._EB
S3 = 0._EB
DO  I = 1,L
  S3 = S3 + H(I)
  S1 = S1 + H(I)*W(I)
END DO

S3   = S3*M*N
PERT = S1/S3


!                               ADJUST F ARRAY BY PERT

DO  K = 1,N
  DO  J = 1,M
    DO  I = 1,L
      F(I,J,K) = F(I,J,K) - PERT
    END DO
  END DO
END DO
750 CONTINUE

!                               IF NORMALIZING SOLUTION, RESTORE PERTRB
!                               AND JUMP TO END

IF (ISING==2) THEN
  PERTRB = PRTSAV

  GO TO 800

END IF

PRTSAV = PERT

!                               SOLVE THE EQUATION

CALL S3CFSS(LDIMF,MDIMF,F,SAVE(IS),W)

!                               IF A SINGULAR PROBLEM,
!                               RE-NORMALIZE SOLUTION (ISING=2)

IF (ISING==1) THEN
  ISING = 2

  GO TO 660

END IF

800 CONTINUE
RETURN
END SUBROUTINE H3CZSS


SUBROUTINE S3CFIS(LPEROD,L,MPEROD,M,NPEROD,N,SCAL,A,B,C,D,LDIMF,  &
    MDIMF,IERROR,SAVE)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+



INTEGER:: LPEROD, L, MPEROD, M, NPEROD, N
REAL(EB):: SCAL
REAL(EB):: A(L), B(L), C(L), D(L), SAVE(-3:*)
INTEGER:: LDIMF, MDIMF, IERROR
INTEGER :: I, LDIMFC, IGRID
!                               CHECK FOR INVALID INPUT

IERROR = 0

IF (LPEROD/=0 .AND. LPEROD/=1) THEN
  IERROR = 1
  SAVE(1) = 1._EB
END IF

IF (L<3) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 2._EB
END IF

IF (MPEROD<0 .AND. MPEROD>4) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 3._EB
END IF

IF (M<3) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 4._EB
END IF

IF (NPEROD<0 .AND. NPEROD>4) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 5._EB
END IF

IF (N<3) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 6._EB
END IF

IF (LDIMF<L) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 7._EB
END IF

IF (MDIMF<M) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 8._EB
END IF

IF (LPEROD==0) THEN
  DO  I = 1,L

    IF (ABS(A(I)-A(1))>=TWO_EPSILON_EB) GO TO 110
    IF (ABS(B(I)-B(1))>=TWO_EPSILON_EB) GO TO 110
    IF (ABS(C(I)-A(1))>=TWO_EPSILON_EB) GO TO 110
    IF (ABS(D(I)-D(1))>=TWO_EPSILON_EB) GO TO 110

  END DO

  GO TO 120

  110     CONTINUE
  IERROR = IERROR + 1
  SAVE(IERROR) = 9._EB
END IF

120 CONTINUE

IF (IERROR/=0) THEN
  RETURN
ELSE
  SAVE(1) = IERROR
END IF

LDIMFC=L
IF (LDIMF>L .AND. MOD(L,2)==0) LDIMFC=L+1

IGRID = 2

CALL FSH00S(IGRID,LPEROD,L,MPEROD,M,NPEROD,N,LDIMFC, SCAL,A,B,C,D,SAVE)

RETURN
END SUBROUTINE S3CFIS


SUBROUTINE S3CFSS(LDIMF,MDIMF,F,SAVE,W)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

INTEGER:: LDIMF, MDIMF
REAL(EB)   F(LDIMF,MDIMF,*), SAVE(-3:*), W(*)

!                               CHECK VALUE OF IERROR (=SAVE(1)).
!                               IF NON-ZERO, RETURN.

IF (ABS(SAVE(1))>=TWO_EPSILON_EB) RETURN

CALL FSH02S(LDIMF,MDIMF,F,SAVE,W)

RETURN
END SUBROUTINE S3CFSS


SUBROUTINE FSH00S(IGRID,LPEROD,L,MPEROD,M,NPEROD,N,LDIMFC,C2, A,B,C,D,SAVE)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+
INTEGER:: IGRID, LPEROD, L,  MPEROD,  M, NPEROD, N, LDIMFC
REAL(EB)   a(l),b(l),c(l),d(l),save(-3:*), C2
INTEGER :: IA, IB, IC, ICFY, ICFZ, IFCTRD, IWSY, IWSZ, IZRT, I, LP, MP, NP
!                               THIS SUBROUTINE INITIALIZES FFT SOLVER

!                               ALLOCATE SAVE ARRAY

IA = 12
IB = IA + L
IC = IB + L
ICFY = IC + L

IF (IGRID==1) THEN
  ICFZ = ICFY + 2*M
  IFCTRD = ICFZ + 2*N
ELSE
  ICFZ = ICFY + 4*M
  IFCTRD = ICFZ + 4*N
END IF

IWSY = IFCTRD + LDIMFC*M*N
IWSZ = IWSY + M + 16
IZRT = IWSZ + N + 16

!                               COPY COEFFICIENT ARRAYS A,B, AND C INTO
!                               SAVE ARRAY.  A COPY OF B IS MADE BECAUSE
!                               IN THE NEXT LEVEL ROUTINE, BOUNDARY
!                               ELEMENTS OF B MAY BE CHANGED.

DO  I = 0,L - 1
  SAVE(IA+I) = A(I+1)
  SAVE(IB+I) = B(I+1)
  SAVE(IC+I) = C(I+1)
END DO
LP = LPEROD + 1
MP = MPEROD + 1
NP = NPEROD + 1

!                               CALL LOWER LEVEL INITIALIZATION ROUTINE

CALL FSH01S(IGRID,L,LP,M,MP,D,N,NP,LDIMFC,C2,SAVE(IA),SAVE(IB),  &
    SAVE(IC),SAVE(ICFY),SAVE(ICFZ),SAVE(IFCTRD), SAVE(IWSY),SAVE(IWSZ),SAVE(IZRT))

!                               SAVE PARAMETERS FOR SUBROUTINE SOLVER

SAVE(2) = REAL(L,EB)
SAVE(3) = REAL(LP,EB)
SAVE(4) = REAL(M,EB)
SAVE(5) = REAL(MP,EB)
SAVE(6) = REAL(N,EB)
SAVE(7) = REAL(NP,EB)
SAVE(8) = REAL(ICFZ,EB)
SAVE(9) = REAL(IWSZ,EB)
SAVE(10) = REAL(IZRT,EB)
SAVE(11) = REAL(IGRID,EB)

RETURN
END SUBROUTINE FSH00S


SUBROUTINE FSH01S(IGRID,L,LP,M,MP,D,N,NP,LDIMFC,C2,A,B,C,CFY,CFZ,  &
    FCTRD,WSAVEY,WSAVEZ,ZRT)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

INTEGER:: IGRID, L, LP, M, MP, N, NP, LDIMFC, K, I
REAL(EB) a(l),b(l),c(l),d(l),cfy(4*m),cfz(4*n),fctrd(ldimfc,n,m),wsavey(m+16),wsavez(n+16),zrt(n)
REAL(EB):: C2
REAL(EB) :: EPS, DEL, DEN, BMAX
INTEGER :: LH, LODD, MRDEL, J, NRDEL


!                               INITIALIZATION ROUTINE FOR FFT SOLVERS

EPS = FSH20S()

IF (L>1 .AND. LP==1) THEN
  LH = (L+1)/2
  LODD = 1

  IF (2*LH==L) LODD = 2
  C(LH-1) = 0._EB
  A(LH) = 0._EB
  C(LH) = 2._EB*C(LH)

  IF (LODD==1) THEN
    B(LH-1) = B(LH-1) - A(LH-1)
    B(L) = B(L) + A(L)
  END IF

  IF (LODD==2) A(L) = C(LH)
END IF

!                               COMPUTE TRANSFORM ROOTS FOR J-DIRECTION

IF (M==1) THEN
  CFY(1) = 0._EB
ELSE

  IF (IGRID==1) THEN
    MRDEL = ((MP-1)* (MP-3)* (MP-5))/3
    DEL = PI/ (2._EB* (M+MRDEL))
  ELSE
    DEL = PI/ (2*M)
  END IF

  IF (MP==1) THEN
    CFY(1) = 0._EB
    CFY(M) = -4._EB
    DO J = 2,M - 1,2
      CFY(J) = -4._EB*SIN(J*DEL)**2
      CFY(J+1) = CFY(J)
    END DO
  END IF

  IF (MP==2) THEN
    DO J = 1,M
      CFY(J) = -4._EB*SIN(J*DEL)**2
    END DO
  END IF

  IF (MP==3 .OR. MP==5) THEN
    DO J = 1,M
      CFY(J) = -4._EB*SIN((J-.5_EB)*DEL)**2
    END DO
  END IF

  IF (MP==4) THEN
    DO J = 1,M
      CFY(J) = -4._EB*SIN((J-1)*DEL)**2
    END DO
  END IF

END IF

!                               COMPUTE TRANSFORM ROOTS IN K-DIRECTION

IF (N==1) THEN
  ZRT(1) = 0._EB
ELSE

  IF (IGRID==1) THEN
    NRDEL = ((NP-1)* (NP-3)* (NP-5))/3
    DEL = PI/ (2._EB* (N+NRDEL))
  ELSE
    DEL = PI/ (2*N)
  END IF

  IF (NP==1) THEN
    ZRT(1) = 0._EB
    ZRT(N) = -4._EB*C2
    DO K = 2,N - 1,2
      ZRT(K) = -4._EB*C2*SIN(K*DEL)**2
      ZRT(K+1) = ZRT(K)
    END DO
  END IF

  IF (NP==2) THEN
    DO K = 1,N
      ZRT(K) = -4._EB*C2*SIN(K*DEL)**2
    END DO
  END IF

  IF (NP==3 .OR. NP==5) THEN
    DO K = 1,N
      ZRT(K) = -4._EB*C2*SIN((K-.5_EB)*DEL)**2
    END DO
  END IF

  IF (NP==4) THEN
    DO K = 1,N
      ZRT(K) = -4._EB*C2*SIN((K-1)*DEL)**2
    END DO
  END IF

END IF

IF (L>1) THEN

!                               FACTOR M*N TRIDIAGONAL SYSTEMS.
!                               FIRST, DO THE POSSIBLY SINGULAR
!                               CASE CORRESPONDING TO J = K = 1.

  FCTRD(1,1,1) = 1._EB/ (B(1)+D(1)*CFY(1)+ZRT(1))
  DO  I = 2,L - 1
    FCTRD(I,1,1) = 1._EB/ (B(I)+D(I)*CFY(1)+ZRT(1)- A(I)*C(I-1)*FCTRD(I-1,1,1))
  END DO

!                               IF TRIDIAGONAL SYSTEM
!                               (...,A(I),B(I),C(I),...) IS SINGULAR
!                               THEN FCTRD(1,1,L) IS 1._EB/0.  IF
!                               DENOMINATOR IS WITHIN ROUND-OFF OF 0,
!                               SET FCTRD(1,1,L) ARBITRARILY

  DEN = B(L) + D(L)*CFY(1) + ZRT(1) - A(L)*C(L-1)*FCTRD(L-1,1,1)
  BMAX = ABS(B(1))
  DO  I = 2,L
    BMAX = MAX(BMAX,ABS(B(I)))
  END DO

  IF (ABS(DEN/BMAX)<=10._EB *EPS) DEN = BMAX
  FCTRD(L,1,1) = 1._EB/DEN

!                               FACTOR CASES J=1, K=2,...,N.

  DO  K = 2,N
    FCTRD(1,K,1) = 1._EB/ (B(1)+D(1)*CFY(1)+ZRT(K))
  END DO
  DO  I = 2,L
    DO  K = 2,N
      FCTRD(I,K,1) = 1._EB/ (B(I)+D(I)*CFY(1)+ZRT(K)- A(I)*C(I-1)*FCTRD(I-1,K,1))
    END DO
  END DO

!                               FACTOR CASES K=1, J=2,...,M.

  DO  J = 2,M
    FCTRD(1,1,J) = 1._EB/ (B(1)+D(1)*CFY(J)+ZRT(1))
  END DO
  DO  I = 2,L
    DO  J = 2,M
      FCTRD(I,1,J) = 1._EB/ (B(I)+D(I)*CFY(J)+ZRT(1)- A(I)*C(I-1)*FCTRD(I-1,1,J))
    END DO
  END DO

!                               FACTOR REMAINING CASES.

  DO  K = 2,N
    DO  J = 2,M
      FCTRD(1,K,J) = 1._EB/ (B(1)+D(1)*CFY(J)+ZRT(K))
    END DO
  END DO
  DO  K = 2,N
    DO  I = 2,L
      DO  J = 2,M
        FCTRD(I,K,J) = 1._EB/ (B(I)+D(I)*CFY(J)+ZRT(K)-  &
            A(I)*C(I-1)*FCTRD(I-1,K,J))
      END DO
    END DO
  END DO
END IF

!                               INITIALIZE FFT TRANSFORMS AND
!                               PRE-PROCESSING COEFFICIENTS IN J

IF (M/=1) THEN

  SELECT CASE(IGRID)
  CASE(1)   ; GO TO 440
  CASE(2)   ; GO TO 450
END SELECT

440     CONTINUE
SELECT CASE(MP)
CASE(1)   ; GO TO 460
CASE(2)   ; GO TO 470
CASE(3)   ; GO TO 480
CASE(4)   ; GO TO 500
CASE(5)   ; GO TO 510
END SELECT

450     CONTINUE
SELECT CASE(MP)
CASE(1)   ; GO TO 460
CASE(2)   ; GO TO 480
CASE(3)   ; GO TO 490
CASE(4)   ; GO TO 510
CASE(5)   ; GO TO 520
END SELECT

460     CONTINUE
CALL VSRFTI(M,WSAVEY)

GO TO 530

470     CONTINUE
CALL VSINTI(M,CFY,WSAVEY)

GO TO 530

480     CONTINUE
!EB          CALL VSSINI(M,CFY(1),CFY(M+1),WSAVEY)
CALL VSCOSI(M,CFY(1),CFY(M+1),WSAVEY)

GO TO 530

490     CONTINUE
!EB          CALL VSSNQI(M,CFY(1),CFY(M+1),CFY(2*M+1),CFY(3*M+1),WSAVEY)
CALL VSCSQI(M,CFY(1),CFY(M+1),CFY(2*M+1),CFY(3*M+1),WSAVEY)

GO TO 530

500     CONTINUE
CALL VCOSTI(M,CFY,WSAVEY)

GO TO 530

510     CONTINUE
CALL VSCOSI(M,CFY(1),CFY(M+1),WSAVEY)

GO TO 530

520     CONTINUE
CALL VSCSQI(M,CFY(1),CFY(M+1),CFY(2*M+1),CFY(3*M+1),WSAVEY)
530     CONTINUE
END IF

!                               INITIALIZE FFT TRANSFORMS AND
!                               PRE-PROCESSING COEFFICIENTS IN K

IF (N/=1) THEN

  SELECT CASE(IGRID)
  CASE(1)   ; GO TO 540
  CASE(2)   ; GO TO 550
END SELECT

540     CONTINUE
SELECT CASE(NP)
CASE(1)   ; GO TO 560
CASE(2)   ; GO TO 570
CASE(3)   ; GO TO 580
CASE(4)   ; GO TO 600
CASE(5)   ; GO TO 610
END SELECT

550     CONTINUE
SELECT CASE(NP)
CASE(1)   ; GO TO 560
CASE(2)   ; GO TO 580
CASE(3)   ; GO TO 590
CASE(4)   ; GO TO 610
CASE(5)   ; GO TO 620
END SELECT

560     CONTINUE
CALL VSRFTI(N,WSAVEZ)

GO TO 630

570     CONTINUE
CALL VSINTI(N,CFZ,WSAVEZ)

GO TO 630

580     CONTINUE
!EB          CALL VSSINI(N,CFZ(1),CFZ(N+1),WSAVEZ)
CALL VSCOSI(N,CFZ(1),CFZ(N+1),WSAVEZ)

GO TO 630

590     CONTINUE
!EB          CALL VSSNQI(N,CFZ(1),CFZ(N+1),CFZ(2*N+1),CFZ(3*N+1),WSAVEZ)
CALL VSCSQI(N,CFZ(1),CFZ(N+1),CFZ(2*N+1),CFZ(3*N+1),WSAVEZ)

GO TO 630

600     CONTINUE
CALL VCOSTI(N,CFZ,WSAVEZ)

GO TO 630

610     CONTINUE
CALL VSCOSI(N,CFZ(1),CFZ(N+1),WSAVEZ)

GO TO 630

620     CONTINUE
CALL VSCSQI(N,CFZ(1),CFZ(N+1),CFZ(2*N+1),CFZ(3*N+1),WSAVEZ)
630     CONTINUE
END IF

RETURN
END SUBROUTINE FSH01S


SUBROUTINE FSH02S(LDIMF,MDIMF,F,SAVE,W)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+



INTEGER:: LDIMF, MDIMF
REAL(EB)   F(LDIMF,MDIMF,*), W(*), SAVE(-3:*)
INTEGER :: L, LP, M, MP, N, NP, IGRID, IA, IC, ICFY, ICFZ, IFCTRD, LDIMFC, IWSY, IWSZ, LDIMFT

!                               RETRIEVE CONSTANTS FROM SAVE ARRAY

L  = NINT(SAVE(2))
LP = NINT(SAVE(3))
M  = NINT(SAVE(4))
MP = NINT(SAVE(5))
N  = NINT(SAVE(6))
NP = NINT(SAVE(7))
IGRID = NINT(SAVE(11))

!                               ALLOCATION OF SAVE ARRAY
IA = 12
IC = IA + 2*L
ICFY = IC + L

IF (IGRID==1) THEN
  ICFZ = ICFY + 2*M
  IFCTRD = ICFZ + 2*N
ELSE
  ICFZ = ICFY + 4*M
  IFCTRD = ICFZ + 4*N
END IF
LDIMFC=L
IF (LDIMF>L .AND. MOD(L,2)==0) LDIMFC=L+1
IWSY = IFCTRD + LDIMFC*M*N
IWSZ = IWSY + M + 16
LDIMFT=L

IF (LDIMF==L .AND. MDIMF==M) THEN

!                               NO HOLES IN DATA ARRAY, SO CALL SOLVER

  CALL FSH03S(IGRID,L,LP,M,MP,N,NP,LDIMFC,LDIMFT,F,SAVE(ICFY), SAVE(ICFZ),  &
      W,SAVE(IA),SAVE(IC),SAVE(IFCTRD),SAVE(IWSY), SAVE(IWSZ))
ELSE
  IF (LDIMF>L .AND. MOD(L,2)==0) LDIMFT=L+1

!                               PACK DATA ARRAY, CALL SOLVER,
!                               AND THEN UNPACK SOLUTION ARRAY

  CALL FSH04S(L,M,N,LDIMF,MDIMF,LDIMFT,F,W)
  CALL FSH03S(IGRID,L,LP,M,MP,N,NP,LDIMFC,LDIMFT,W,SAVE(ICFY),  &
      SAVE(ICFZ),F, SAVE(IA),SAVE(IC),SAVE(IFCTRD),SAVE(IWSY),  &
      SAVE(IWSZ))
  CALL FSH05S(L,M,N,LDIMF,MDIMF,LDIMFT,F,W)
END IF

RETURN
END SUBROUTINE FSH02S


SUBROUTINE FSH03S(IGRID,L,LP,M,MP,N,NP,LDIMFC,LDIMFT,F,CFY,CFZ,  &
    FT,A,C,FCTRD,WSAVEY,WSAVEZ)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+


INTEGER                   :: IGRID
INTEGER                       :: L
INTEGER                   :: LP
INTEGER                       :: M
INTEGER                   :: MP
INTEGER                       :: N
INTEGER                   :: NP
INTEGER                   :: LDIMFC
INTEGER                   :: LDIMFT
REAL(EB)  CFY(4*M)
REAL(EB)  CFZ(4*N)
REAL(EB)  FT(LDIMFT,M,N)
REAL(EB)   FCTRD(LDIMFC,M,N)
REAL(EB)   WSAVEY(M+16)
REAL(EB)   WSAVEZ(N+16)
REAL(EB)   A(L),C(L),F(LDIMFT,M,N)
INTEGER :: K, J, IFWRD, I


LOGICAL :: DATARY,DATASW
!                               ZERO OUT BOTTOM PLANE OF ARRAY FT

DO  K=1,N
  DO  J=1,M
    FT(LDIMFT,J,K)=0._EB
  END DO
END DO

NOCOPY=.TRUE.
DATARY=.TRUE.
SCALE=1._EB
IFWRD = 1
100 CONTINUE

IF (N/=1) THEN
  TPOSE=.FALSE.
  IF (IFWRD==2) TPOSE=.TRUE.

!                               TRANSFORM IN Z
  IF (DATARY) THEN
    CALL FSH26S(IGRID,IFWRD,NP,L,N,M,LDIMFT,F,FT,CFZ,WSAVEZ)
    DATARY=OUTARY

  ELSE
    CALL FSH26S(IGRID,IFWRD,NP,L,N,M,LDIMFT,FT,F,CFZ,WSAVEZ)
    DATARY=.NOT.OUTARY
  END IF

END IF
SELECT CASE(IFWRD)
CASE(1)   ; GO TO 490
CASE(2)   ; GO TO 510
END SELECT

490 CONTINUE
IF (M/=1) THEN
  TPOSE=.TRUE.
  IF (IFWRD==2) TPOSE=.FALSE.

!                               TRANSFORM Y


  IF (DATARY) THEN
    CALL FSH26S(IGRID,IFWRD,MP,L,M,N,LDIMFT,F,FT,CFY,WSAVEY)
    DATARY=OUTARY
  ELSE
    CALL FSH26S(IGRID,IFWRD,MP,L,M,N,LDIMFT,FT,F,CFY,WSAVEY)
    DATARY=.NOT.OUTARY
  END IF

END IF

SELECT CASE(IFWRD)
CASE(1)   ; GO TO 285
CASE(2)   ; GO TO 100
END SELECT
285 CONTINUE

IF (L>1) THEN

!                               SOLVE TRIDIAGONAL SYSTEMS IN X THAT WERE
!                               PREVIOUSLY FACTORED IN FSH01S

!                               CALL VECTORIZED TRIDIAGONAL SOLVER

  DATASW=.FALSE.
  IF (NP==1) DATASW=.NOT.DATASW
  IF (MP==1) DATASW=.NOT.DATASW
  IF (DATARY) THEN
    IF (DATASW) THEN
      CALL FSH06S(L,LP,M*N,LDIMFC,LDIMFT,SCALE,A,C,F,FT,FCTRD)
      DATARY=.FALSE.
    ELSE
      CALL FSH06S(L,LP,M*N,LDIMFC,LDIMFT,SCALE,A,C,F,F,FCTRD)
    END IF
  ELSE
    IF (DATASW) THEN
      CALL FSH06S(L,LP,M*N,LDIMFC,LDIMFT,SCALE,A,C,FT,F,FCTRD)
      DATARY=.TRUE.
    ELSE
      CALL FSH06S(L,LP,M*N,LDIMFC,LDIMFT,SCALE,A,C,FT,FT,FCTRD)
    END IF
  END IF

END IF

IFWRD = 2

GO TO 490

510 CONTINUE

IF (.NOT.DATARY) THEN
  DO  K=1,N
    DO  J=1,M
      DO  I=1,L
        F(I,J,K)=FT(I,J,K)
      END DO
    END DO
  END DO
END IF
RETURN
END SUBROUTINE FSH03S


SUBROUTINE FSH04S(L,M,N,LDIMF,MDIMF,LDIMG,F,G)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+



INTEGER                       :: L
INTEGER                       :: M
INTEGER                       :: N
INTEGER                   :: LDIMF
INTEGER                   :: MDIMF
INTEGER                   :: LDIMG
REAL(EB)   F(LDIMF,MDIMF, N), G(LDIMG,M,N)
INTEGER :: K, J, I

!                              THIS SUBROUTINE PACKS THE SUB-ARRAY
!                              F(I,J,K), I=1,...,L, J=1,...,M, K=1,...,N
!                              INTO THE ARRAY G.

DO  K = 1,N
  DO  J = 1,M
    DO  I = 1,L
      G(I,J,K) = F(I,J,K)
    END DO
  END DO
END DO
IF (LDIMG>L) THEN
  DO  K=1,N
    DO  J=1,M
      G(LDIMG,J,K)=0._EB
    END DO
  END DO
END IF
RETURN
END SUBROUTINE FSH04S


SUBROUTINE FSH05S(L,M,N,LDIMF,MDIMF,LDIMG,F,G)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+



INTEGER                       :: L
INTEGER                       :: M
INTEGER                       :: N
INTEGER                   :: LDIMF
INTEGER                   :: MDIMF
INTEGER                   :: LDIMG
REAL(EB)   G(LDIMG,M,N)
REAL(EB)   F(LDIMF,MDIMF,N)
INTEGER :: K, J, I

!                               THIS SUBROUTINE EXPANDS THE ARRAY G OF
!                               DIMENSION L X M X N INTO THE ARRAY F OF
!                               DIMENSION LDIMF X MDIMF X N.

DO  K = 1,N
  DO  J = 1,M
    DO  I = 1,L
      F(I,J,K) = G(I,J,K)
    END DO
  END DO
END DO
RETURN
END SUBROUTINE FSH05S


SUBROUTINE FSH06S(L,LP,M,LDIMFC,LDIMFT,SCALE,A,C,F,FT,FCTRD)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+



INTEGER                       :: L
INTEGER                   :: LP
INTEGER                       :: M
INTEGER                   :: LDIMFC
INTEGER                   :: LDIMFT
REAL(EB) :: SCALE
REAL(EB)   FT(LDIMFT,M)
REAL(EB)  FCTRD(LDIMFC,M)
REAL(EB)   A(L),C(L),F(LDIMFT,M)
INTEGER :: LH, LQ, I, J

!                               THIS SUBROUTINE SOLVES M TRIDIAGONAL
!                               SYSTEMS OF ORDER L THAT HAVE
!                               FACTORIZATION STORED IN ARRAY FCTRD AND
!                               RIGHT SIDES IN FT
SCALE=SCALE**2
IF (LP==1) THEN
  LH=(L+1)/2
  LQ=(LH-1)/2
  DO  I=1,LQ
    DO  J=1,M
      F(I,J)=F(I,J)-F(2*LH-I,J)
      F(LH-I,J)=F(LH-I,J)-F(LH+I,J)
      F(LH+I,J)=F(LH-I,J)+2._EB*F(LH+I,J)
      F(2*LH-I,J)=F(I,J)+2._EB*F(2*LH-I,J)
    END DO
    CALL SSWAP(M,F(I,1),LDIMFT,F(LH-I,1),LDIMFT)
  END DO
  DO  J=1,M
    F(LH,J)=2._EB*F(LH,J)
  END DO
IF (MOD(L,2)==0)  THEN
  DO  J=1,M
    F(L,J)=2._EB*F(L,J)
  END DO
END IF
IF (MOD(LH,2)==0) THEN
  DO  J=1,M
    F(  LH/2,J)=F(LH/2,J)-F(3*LH/2,J)
    F(3*LH/2,J)=F(LH/2,J)+2._EB*F(3*LH/2,J)
  END DO
END IF
SCALE=0.5_EB*SCALE
END IF
!                               FORWARD SUBSTITUTION

DO  J = 1,M
  FT(1,J) = SCALE*F(1,J)*FCTRD(1,J)
  DO  I = 2,L
    FT(I,J) = (SCALE*F(I,J)-A(I)*FT(I-1,J))*FCTRD(I,J)
  END DO
END DO

!                               BACKWARD SUBSTITUTION
DO  J = 1,M
  DO  I = L - 1,1,-1
    FT(I,J) = FT(I,J) - C(I)*FCTRD(I,J)*FT(I+1,J)
  END DO
END DO

IF (LP==1) THEN

  DO  I=1,LQ
    DO  J=1,M
      FT(I,J)=FT(LH+I,J)+FT(I,J)
      FT(LH-I,J)=FT(2*LH-I,J)+FT(LH-I,J)
      FT(2*LH-I,J)=2._EB*FT(2*LH-I,J)-FT(LH-I,J)
      FT(  LH+I,J)=2._EB*FT(LH+I,J)-FT(I,J)
    END DO
    CALL SSWAP(M,FT(I,1),LDIMFT,FT(LH-I,1),LDIMFT)
  END DO
  IF (MOD(LH,2)==0) THEN
    DO  J=1,M
      FT(  LH/2,J)=FT(3*LH/2,J)+FT(LH/2,J)
      FT(3*LH/2,J)=2._EB*FT(3*LH/2,J)-FT(LH/2,J)
    END DO
  END IF

END IF

RETURN
END SUBROUTINE FSH06S



FUNCTION FSH20S()

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!                               THIS FUNCTION FURNISHES THE MACHINE
!                               PREC EPS--THE SMALLEST POSITIVE
!                               MACHINE NUMBER SATISFYING

!                               FL(1._EB+EPS) > 1.

!     FSH20S = 2.4E-7_EB
REAL(EB) FSH20S
FSH20S = SPACING (1._EB)!4.E-15_EB

RETURN
END FUNCTION FSH20S


SUBROUTINE FSH26S(IGRID,IFWRD,MP,L,M,N,LDIMFT,F,FT,CFY,WSAVEY)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+


INTEGER                   :: IGRID
INTEGER                   :: IFWRD
INTEGER                   :: MP
INTEGER                   :: L
INTEGER                   :: M
INTEGER                   :: N
INTEGER                   :: LDIMFT
REAL(EB) FT(LDIMFT,M,N)
REAL(EB) CFY(4*M)
REAL(EB) WSAVEY(M+16)
REAL(EB) F(LDIMFT,M,N)

SELECT CASE(IGRID)
CASE(1)   ; GO TO 110
CASE(2)   ; GO TO 120
END SELECT

110     CONTINUE
SELECT CASE(MP)
CASE(1)   ; GO TO 130
CASE(2)   ; GO TO 160
CASE(3)   ; GO TO 170
CASE(4)   ; GO TO 220
CASE(5)   ; GO TO 230
END SELECT
120     CONTINUE
SELECT CASE(MP)
CASE(1)   ; GO TO 130
CASE(2)   ; GO TO 180
CASE(3)   ; GO TO 210
CASE(4)   ; GO TO 240
CASE(5)   ; GO TO 270
END SELECT
130     CONTINUE
SELECT CASE(IFWRD)
CASE(1)   ; GO TO 140
CASE(2)   ; GO TO 150
END SELECT
140     CONTINUE

CALL VSRFTF(F,L,M,N,LDIMFT,FT,WSAVEY)
GO TO 280
150     CONTINUE
CALL VSRFTB(F,L,M,N,LDIMFT,FT,WSAVEY)
GO TO 280
160     CONTINUE
CALL VSINT(F,L,M,N,LDIMFT,FT,CFY,WSAVEY)
GO TO 280
170     CONTINUE
SELECT CASE(IFWRD)
CASE(1)   ; GO TO 200
CASE(2)   ; GO TO 190
END SELECT
180     CONTINUE
SELECT CASE(IFWRD)
CASE(1)   ; GO TO 190
CASE(2)   ; GO TO 200
END SELECT
190     CONTINUE
CALL VSSINF(F,L,M,N,LDIMFT,FT,CFY(1),CFY(M+1),WSAVEY)
GO TO 280
200     CONTINUE
CALL VSSINB(F,L,M,N,LDIMFT,FT,CFY(1),CFY(M+1),WSAVEY)
GO TO 280
210     CONTINUE
CALL VSSINQ(F,L,M,N,LDIMFT,FT,CFY(1),CFY(M+1), CFY(2*M+1),CFY(3*M+1),WSAVEY)
GO TO 280
220     CONTINUE
CALL VCOST(F,L,M,N,LDIMFT,FT,CFY,WSAVEY)
GO TO 280
230     CONTINUE
SELECT CASE(IFWRD)
CASE(1)   ; GO TO 260
CASE(2)   ; GO TO 250
END SELECT
240     CONTINUE
SELECT CASE(IFWRD)
CASE(1)   ; GO TO 250
CASE(2)   ; GO TO 260
END SELECT
250     CONTINUE
CALL VSCOSF(F,L,M,N,LDIMFT,FT,CFY(1),CFY(M+1),WSAVEY)
GO TO 280
260     CONTINUE
CALL VSCOSB(F,L,M,N,LDIMFT,FT,CFY(1),CFY(M+1),WSAVEY)
GO TO 280
270     CONTINUE
CALL VSCOSQ(F,L,M,N,LDIMFT,FT,CFY(1),CFY(M+1), CFY(2*M+1),CFY(3*M+1),WSAVEY)
280     CONTINUE

RETURN
END SUBROUTINE FSH26S


SUBROUTINE VCOST(X,L,M,N,LDIMX,XT,C,WSAVE)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

INTEGER                   :: L
INTEGER                       :: M
INTEGER                       :: N
INTEGER                       :: LDIMX
REAL(EB) XT(LDIMX*N,M)
REAL(EB) C(M)
REAL(EB) WSAVE(M+15)
REAL(EB) X(LDIMX*N,M)
INTEGER :: MM1, MS2, I, J, JC

!     TPOSE = .TRUE. IF TRANSFORMING SECOND INDEX (SO TRANSPOSE)
!           = .FALSE. IF TRANSFORMING THIRD INDEX

MM1 = M-1
MS2 = M/2
SCALE=SCALE*SQRT(0.5_EB)
IF (TPOSE) THEN
  CALL VCOST1(L,M,N,LDIMX,MS2,C,XT,X)
ELSE
  DO  I=1,LDIMX*N
    XT(I,M)=  X(I,1)-X(I,M)
    XT(I,1) = X(I,1)+X(I,M)
  END DO
DO  J=2,MS2
  JC = M+1-J
  DO  I=1,LDIMX*N
    XT(I,M) = XT(I,M)+C(JC)*(X(I,J)-X(I,JC))
    XT(I,J) =  (X(I,J)+X(I,JC))-(C(J)*(X(I,J)-X(I,JC)))
    XT(I,JC) = (X(I,J)+X(I,JC))+(C(J)*(X(I,J)-X(I,JC)))
  END DO
END DO
IF (MOD(M,2) == 0) GO TO 404
DO  I=1,LDIMX*N
  XT(I,MS2+1) = 2._EB*X(I,MS2+1)
END DO
404    CONTINUE
END IF
IF (M>3) THEN
  CALL VRFFTF (LDIMX*N,MM1,XT,LDIMX*N,X,WSAVE)
  IF (OUTARY) THEN
    CALL VCOSTA(M,N*LDIMX,MM1,XT(1,M),X,XT)
  ELSE
    CALL VCOSTA(M,N*LDIMX,MM1,XT(1,M),XT,X)
  END IF
ELSE IF (M==2) THEN
  OUTARY=.FALSE.
ELSE
  DO  I=1,LDIMX*N
    X(I,2)=XT(I,3)
    X(I,1)=XT(I,1)+XT(I,2)
    X(I,3)=XT(I,1)-XT(I,2)
  END DO
OUTARY=.TRUE.
SCALE=SCALE*SQRT(0.5_EB)
END IF
RETURN
END SUBROUTINE VCOST


SUBROUTINE VCOST1(L,M,N,LDIMX,MS2,C,XT,X)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989



INTEGER                       :: L
INTEGER                       :: M
INTEGER                       :: N
INTEGER                   :: LDIMX
INTEGER                   :: MS2
REAL(EB)   C(*)
REAL(EB)  XT(LDIMX,N,M)
REAL(EB)   X(LDIMX,M,N)
INTEGER :: K, I, J, JC

DO  K=1,N
  DO  I=1,L
    XT(I,K,M)=X(I,1,K)-X(I,M,K)
    XT(I,K,1) = X(I,1,K)+X(I,M,K)
  END DO
END DO
DO  J=2,MS2
  JC = M+1-J
  DO  K=1,N
    DO  I=1,L
      XT(I,K,M) = XT(I,K,M)+C(JC)*(X(I,J,K)-X(I,JC,K))
      XT(I,K,J) =  (X(I,J,K)+X(I,JC,K))-(C(J)*(X(I,J,K)-X(I,JC,K)))
      XT(I,K,JC) = (X(I,J,K)+X(I,JC,K))+(C(J)*(X(I,J,K)-X(I,JC,K)))
    END DO
  END DO
END DO
IF (MOD(M,2) == 0) GO TO 404
DO  K=1,N
  DO  I=1,L
    XT(I,K,MS2+1) = 2._EB*X(I,MS2+1,K)
  END DO
END DO
404 CONTINUE
RETURN
END SUBROUTINE VCOST1


SUBROUTINE VCOSTA(M,LDIMX,MM1,PL,X,XT)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989



INTEGER                       :: M
INTEGER                       :: LDIMX
INTEGER                   :: MM1
REAL(EB)   X(LDIMX,M), PL(LDIMX),XT(LDIMX,M)
INTEGER :: I, J

DO  I=1,LDIMX
  X(I,1) = XT(I,1)
  X(I,2) = PL(I)
END DO
DO  J=4,M,2
  DO  I=1,LDIMX
    X(I,J) = X(I,J-2)-XT(I,J-1)
    X(I,J-1) = XT(I,J-2)
  END DO
END DO
IF (MOD(M,2) == 0) GO TO 409
DO  I=1,LDIMX
  X(I,M) = XT(I,MM1)
END DO
409 CONTINUE
RETURN
END SUBROUTINE VCOSTA


SUBROUTINE VCOSTI(N,C,WSAVE)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989


!     *******   NOTE   ******

!     SPECIAL DOCUMENTATION FOR USE OF INTERNAL PARAMETER TPOSE.  IF
!     TPOSE = .TRUE., THEN ALL PREPROCESSING ROUTINES TRANSPOSE SECOND
!                     AND THIRD INDICES IN PREPARATION FOR TRANSFORM OF
!                     ORIGINAL SECOND INDEX.  IF
!     TPOSE = .FALSE., THEN ALL PREPROCESSING ROUTINES PREPARE ONLY FOR
!                      TRANSFORM OF THIRD INDEX.  TO USE IN THIS CASE,
!                      CALL ROUTINE WITH M AND N INTERCHANGED.



INTEGER                       :: N
REAL(EB)    C(N), WSAVE(N+15)
INTEGER :: NP1, NS2, K, KC
REAL(EB) :: DT


!                               INITIALIZE NOCOPY AND TPOSE TO DEFAULT
!                               VALUES

NOCOPY = .FALSE.
TPOSE = .TRUE.

IF (N <= 3) RETURN
NP1 = N+1
NS2 = N/2
DT = PI/(N-1)

DO  K=2,NS2
  KC = NP1-K
  C(K) = 2._EB*SIN((K-1)*DT)
  C(KC) = 2._EB*COS((K-1)*DT)
END DO
CALL VRFFTI (N-1,WSAVE)
RETURN
END SUBROUTINE VCOSTI


SUBROUTINE VRFFTF (M,N,R,MDIMR,RT,WSAVE)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989



INTEGER                   :: M
INTEGER                   :: N
INTEGER                   :: MDIMR
REAL(EB)   RT(M,N)
REAL(EB)   WSAVE(N+15)
REAL(EB)         R(MDIMR,N)
IF (N == 1) RETURN
CALL VRFTF1 (M,N,R,MDIMR,RT,WSAVE(1),WSAVE(N+1))
RETURN
END SUBROUTINE VRFFTF


SUBROUTINE VRFFTI (N,WSAVE)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989



INTEGER                   :: N
REAL(EB)         WSAVE(N+15)

!                               INITIALIZE NOCOPY AND TPOSE TO DEFAULT
!                               VALUES

NOCOPY = .FALSE.
TPOSE = .TRUE.

IF (N <= 1) RETURN
CALL VRFTI1 (N,WSAVE(1),WSAVE(N+1))
RETURN
END SUBROUTINE VRFFTI


SUBROUTINE VRFTF1 (M,N,C,MDIMC,CH,WA,FAC)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989



INTEGER                       :: M
INTEGER                       :: N
INTEGER                   :: MDIMC
REAL(EB)    CH(M,N)
REAL(EB)   WA(N)
REAL(EB)  FAC(15)
REAL(EB)         C(MDIMC,N)
INTEGER :: NF, NA, L2, IW, K1, KH, IP, L1, IDO, IDL1, IX2, IX3, IX4, J, I

NF = INT(FAC(2))
NA = 1
L2 = N
IW = N
DO  K1=1,NF
  KH = NF-K1
  IP = INT(FAC(KH+3))
  L1 = L2/IP
  IDO = N/L2
  IDL1 = IDO*L1
  IW = IW-(IP-1)*IDO
  NA = 1-NA
  IF (IP /= 4) GO TO 102
  IX2 = IW+IDO
  IX3 = IX2+IDO
  IF (NA /= 0) GO TO 101
  CALL VRADF4 (M,IDO,L1,C,MDIMC,CH,M,WA(IW),WA(IX2),WA(IX3))
  GO TO 110
  101    CALL VRADF4 (M,IDO,L1,CH,M,C,MDIMC,WA(IW),WA(IX2),WA(IX3))
  GO TO 110
  102    IF (IP /= 2) GO TO 104
  IF (NA /= 0) GO TO 103
  CALL VRADF2 (M,IDO,L1,C,MDIMC,CH,M,WA(IW))
  GO TO 110
  103    CALL VRADF2 (M,IDO,L1,CH,M,C,MDIMC,WA(IW))
  GO TO 110
  104    IF (IP /= 3) GO TO 106
  IX2 = IW+IDO
  IF (NA /= 0) GO TO 105
  CALL VRADF3 (M,IDO,L1,C,MDIMC,CH,M,WA(IW),WA(IX2))
  GO TO 110
  105    CALL VRADF3 (M,IDO,L1,CH,M,C,MDIMC,WA(IW),WA(IX2))
  GO TO 110
  106    IF (IP /= 5) GO TO 108
  IX2 = IW+IDO
  IX3 = IX2+IDO
  IX4 = IX3+IDO
  IF (NA /= 0) GO TO 107
  CALL VRADF5(M,IDO,L1,C,MDIMC,CH,M,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  GO TO 110
  107 CALL VRADF5 (M,IDO,L1,CH,M,C,MDIMC,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  GO TO 110
  108    IF (IDO == 1) NA = 1-NA
  IF (NA /= 0) GO TO 109
  CALL VRADFG (M,IDO,IP,L1,IDL1,C,C,C,MDIMC,CH,CH,M,WA(IW))
  NA = 1
  GO TO 110
  109    CALL VRADFG (M,IDO,IP,L1,IDL1,CH,CH,CH,M,C,C,MDIMC,WA(IW))
  NA = 0
  110    L2 = L1
END DO
OUTARY=.TRUE.
IF (NOCOPY) THEN
  SCALE=SCALE*SQRT(1.0_EB/REAL(N,EB))
  IF (NA==0) THEN
    OUTARY=.FALSE.
  END IF
ELSE
  SCALE=SQRT(1.0_EB/REAL(N,EB))
  IF (NA == 1) GO TO 113
  DO  J=1,N
    DO  I=1,M
      C(I,J) = SCALE*CH(I,J)
    END DO
  END DO
  RETURN
  113    DO  J=1,N
    DO  I=1,M
      C(I,J)=SCALE*C(I,J)
    END DO
  END DO
END IF
RETURN
END SUBROUTINE VRFTF1


SUBROUTINE VRFTI1 (N,WA,FAC)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989


INTEGER :: N
REAL(EB) :: FAC(15), WA(N)
INTEGER :: NL, NF, J, NTRY, NQ, NR, I, IB, IS, NFM1, L1, K1, IP, LD, L2, IDO, IPM, II, NTRYH(4)
REAL(EB) :: ARGH, TPI, ARGLD, FI, ARG
DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/4,2,3,5/

NL = N
NF = 0
J = 0
101 J = J+1
IF (J-4<0) GO TO 102
IF (J-4==0) GO TO 102
IF (J-4>0) GO TO 103
102 NTRY = NTRYH(J)
GO TO 104
103 NTRY = NTRY+2
104 NQ = NL/NTRY
NR = NL-NTRY*NQ
IF (NR<0) GO TO 101
IF (NR==0) GO TO 105
IF (NR>0) GO TO 101
105 NF = NF+1
FAC(NF+2) = NTRY
NL = NQ
IF (NTRY /= 2) GO TO 107
IF (NF == 1) GO TO 107
DO  I=2,NF
  IB = NF-I+2
  FAC(IB+2) = FAC(IB+1)
END DO
FAC(3) = 2
107 IF (NL /= 1) GO TO 104
FAC(1) = N
FAC(2) = NF
TPI = 2._EB*PI
ARGH = TPI/REAL(N,EB)
IS = 0
NFM1 = NF-1
L1 = 1
IF (NFM1 == 0) RETURN
DO  K1=1,NFM1
  IP = INT(FAC(K1+2))
  LD = 0
  L2 = L1*IP
  IDO = N/L2
  IPM = IP-1
  DO  J=1,IPM
    LD = LD+L1
    I = IS
    ARGLD = REAL(LD,EB)*ARGH
    FI = 0._EB
    DO  II=3,IDO,2
      I = I+2
      FI = FI+1._EB
      ARG = FI*ARGLD
      WA(I-1) = COS(ARG)
      WA(I) = SIN(ARG)
    END DO
    IS = IS+IDO
  END DO
  L1 = L2
END DO
RETURN
END SUBROUTINE VRFTI1


SUBROUTINE VSCOSB(F,L,M,N,LDIMF,FT,C1,C2,WORK)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

INTEGER                   :: L
INTEGER                       :: M
INTEGER                       :: N
INTEGER                       :: LDIMF
REAL(EB)   FT(LDIMF*N,M)
REAL(EB)  C1(M)
REAL(EB)    C2(M)
REAL(EB)   WORK(M+15)
REAL(EB)   F(LDIMF*N,M)
INTEGER :: I, J

!     TPOSE = .TRUE. IF TRANSFORMING SECOND INDEX (SO TRANSPOSE)
!           = .FALSE. IF TRANSFORMING THIRD INDEX


!     PREPROCESSING

IF (TPOSE) THEN
  CALL VSCSB1(L,M,N,LDIMF,F,FT,C1,C2)
ELSE
  DO  I=1,LDIMF*N
    FT(I,1)=.5_EB*F(I,1)
  END DO
DO  J=2,M
  DO  I=1,LDIMF*N
    FT(I,J) =(C1(J)*F(I,J)+C2(J)*F(I,M-J+2))
  END DO
END DO
END IF

!     REAL(EB),PERIODIC ANALYSIS

CALL VRFFTF(LDIMF*N,M,FT,LDIMF*N,F,WORK)

!     POSTPROCESSING

SCALE=SQRT(2.0_EB)*SCALE
IF (OUTARY) THEN
  CALL VSCSBA(M,N*LDIMF,FT,F)
ELSE
  CALL VSCSBA(M,N*LDIMF,F,FT)
END IF
RETURN
END SUBROUTINE VSCOSB


SUBROUTINE VSCOSF(F,L,M,N,LDIMF,FT,C1,C2,WORK)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989


INTEGER                   :: L
INTEGER                       :: M
INTEGER                       :: N
INTEGER                       :: LDIMF
REAL(EB)    FT(LDIMF*N,M)
REAL(EB)  C1(M)
REAL(EB)   C2(M)
REAL(EB)   WORK(M+15)
REAL(EB)   F(LDIMF*N,M)
INTEGER :: I, J

!     TPOSE = .TRUE. IF TRANSFORMING SECOND INDEX (SO TRANSPOSE)
!           = .FALSE. IF TRANSFORMING THIRD INDEX


!     PREPROCESSING

SCALE=SQRT(2.0_EB)*SCALE
IF (TPOSE) THEN
  CALL VSCSF1(L,M,N,LDIMF,F,FT)
ELSE
  DO  I=1,LDIMF*N
    FT(I,1)=F(I,1)
  END DO
CONTINUE
IF (MOD(M,2)==0) THEN
  DO  I=1,LDIMF*N
    FT(I,M)=F(I,M)
  END DO
END IF
DO  J=2,M-1,2
  DO  I=1,LDIMF*N
    FT(I,J)   =  0.5_EB*(F(I,J)+F(I,J+1))
    FT(I,J+1) = -0.5_EB*(F(I,J)-F(I,J+1))
  END DO
END DO
END IF

!     REAL(EB),PERIODIC SYNTHESIS

CALL VRFFTB(LDIMF*N,M,FT,LDIMF*N,F,WORK)

!     POSTPROCESSING

IF (OUTARY) THEN
  CALL VSCSFA(M,N*LDIMF,F,FT,C1,C2)
ELSE
  CALL VSCSFA(M,N*LDIMF,FT,F,C1,C2)
END IF
RETURN
END SUBROUTINE VSCOSF


SUBROUTINE VSCOSI(N,C1,C2,WSAVE)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

!      ENTRY VSSINI(N,C1,C2,WSAVE)


INTEGER                       :: N
REAL(EB)    WSAVE(N+15)
REAL(EB)   C1(N),C2(N)
INTEGER :: I
REAL(EB) :: DX, C, S

DX=PI/(2*N)

!     GENERATE A(I)+-B(I)

DO  I=1,N
  C=COS((I-1)*DX)
  S=SIN((I-1)*DX)
  C1(I)=.5_EB*(S+C)
  C2(I)=.5_EB*(S-C)
END DO

!     INITIALIZE VRFFTPK ROUTINE

CALL VRFFTI(N,WSAVE)
RETURN
END SUBROUTINE VSCOSI


SUBROUTINE VSCOSQ(F,L,M,N,LDIMF,FT,C1,C2,C3,C4,WORK)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989


INTEGER                   :: L
INTEGER                       :: M
INTEGER                       :: N
INTEGER                       :: LDIMF
REAL(EB)  FT(LDIMF*N,M)
REAL(EB)  C1(M)
REAL(EB)  C2(M)
REAL(EB)  C3(M)
REAL(EB)  C4(M)
REAL(EB)  WORK(M+15)
REAL(EB)   F(LDIMF*N,M)
INTEGER :: I, J, JBY2

!     TPOSE = .TRUE. IF TRANSFORMING SECOND INDEX (SO TRANSPOSE)
!           = .FALSE. IF TRANSFORMING THIRD INDEX


!     PREPROCESSING

IF (TPOSE) THEN
   CALL VSCSQ1(L,M,N,LDIMF,F,FT,C1,C2)
ELSE
   DO  I=1,LDIMF*N
      FT(I,1)=F(I,1)
   END DO
   IF (MOD(M,2)==0) THEN
      DO  I=1,LDIMF*N
         FT(I,M)=-F(I,M)
      END DO
   END IF
      DO  J=2,M-1,2
         JBY2=J/2
         DO  I=1,LDIMF*N
            FT(I,J)   =  F(I,J+1)*C1(JBY2)+F(I,J)*C2(JBY2)
            FT(I,J+1) = -F(I,J+1)*C2(JBY2)+F(I,J)*C1(JBY2)
         END DO
      END DO
END IF

!     REAL(EB),PERIODIC SYNTHESIS

CALL VRFFTB(LDIMF*N,M,FT,LDIMF*N,F,WORK)

!     POSTPROCESSING

IF (OUTARY) THEN
  DO  J=1,M
    DO  I=1,LDIMF*N
      F(I,J)=C4(J)*FT(I,J)+C3(J)*FT(I,M+1-J)
    END DO
  END DO
ELSE
   DO  J=1,M
      DO  I=1,LDIMF*N
         FT(I,J)=C4(J)*F(I,J)+C3(J)*F(I,M+1-J)
      END DO
   END DO
END IF
RETURN
END SUBROUTINE VSCOSQ


SUBROUTINE VSCSB1(L,M,N,LDIMF,F,FT,C1,C2)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989



INTEGER                       :: L
INTEGER                       :: M
INTEGER                       :: N
INTEGER                   :: LDIMF
REAL(EB)  FT(LDIMF,N,M)
REAL(EB)  C1(M)
REAL(EB)  C2(M)
REAL(EB)   F(LDIMF,M,N)
INTEGER :: K, I, J

DO  K=1,N
  DO  I=1,L
    FT(I,K,1)=.5_EB*F(I,1,K)
  END DO
END DO
DO  K=1,N
  DO  J=2,M
    DO  I=1,L
      FT(I,K,J) =(C1(J)*F(I,J,K)+C2(J)*F(I,M-J+2,K))
    END DO
  END DO
END DO
RETURN
END SUBROUTINE VSCSB1


SUBROUTINE VSCSBA(M,LDIMF,FT,F)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989



INTEGER                       :: M
INTEGER                       :: LDIMF
REAL(EB)   F(LDIMF,M), FT(LDIMF,M)
INTEGER :: I, J

DO  I=1,LDIMF
  F(I,1) = FT(I,1)
END DO
IF (MOD(M,2)==0) THEN
  DO  I=1,LDIMF
    F(I,M) = FT(I,M)
  END DO
END IF
DO  J=2,M-1,2
  DO  I=1,LDIMF
    F(I,J)=FT(I,J)-FT(I,J+1)
    F(I,J+1)=FT(I,J)+FT(I,J+1)
  END DO
END DO
RETURN
END SUBROUTINE VSCSBA


SUBROUTINE VSCSF1(L,M,N,LDIMF,F,FT)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989



INTEGER                       :: L
INTEGER                       :: M
INTEGER                       :: N
INTEGER                   :: LDIMF
REAL(EB)   FT(LDIMF,N,M)
REAL(EB)   F(LDIMF,M,N)
INTEGER :: K, I, J
DO  K=1,N
  DO  I=1,L
    FT(I,K,1)=F(I,1,K)
  END DO
END DO
IF (MOD(M,2)==0) THEN
  DO  K=1,N
    DO  I=1,L
      FT(I,K,M)=F(I,M,K)
    END DO
  END DO
END IF
DO  K=1,N
  DO  J=2,M-1,2
    DO  I=1,L
      FT(I,K,J)   =  0.5_EB*(F(I,J,K)+F(I,J+1,K))
      FT(I,K,J+1) = -0.5_EB*(F(I,J,K)-F(I,J+1,K))
    END DO
  END DO
END DO
RETURN
END SUBROUTINE VSCSF1


SUBROUTINE VSCSFA(M,LDIMF,F,FT,C1,C2)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989



INTEGER                       :: M
INTEGER                       :: LDIMF
REAL(EB)    FT(LDIMF,M)
REAL(EB)    C1(M)
REAL(EB)   F(LDIMF,M),  C2(M)
INTEGER :: J, I

DO  J=2,M
  DO  I=1,LDIMF
    F(I,J)=C1(J)*FT(I,J)-C2(J)*FT(I,M+2-J)
  END DO
END DO
DO  I=1,LDIMF
  F(I,1)=FT(I,1)
END DO
RETURN
END SUBROUTINE VSCSFA


SUBROUTINE VSCSQ1(L,M,N,LDIMF,F,FT,C1,C2)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989



INTEGER                       :: L
INTEGER                       :: M
INTEGER                       :: N
INTEGER                   :: LDIMF
REAL(EB)  FT(LDIMF,N,M)
REAL(EB)  C1(M)
REAL(EB)    C2(M)
REAL(EB)   F(LDIMF,M,N)
INTEGER :: K, I, J, JBY2

DO  K=1,N
  DO  I=1,L
    FT(I,K,1)=F(I,1,K)
  END DO
END DO
IF (MOD(M,2)==0) THEN
  DO  K=1,N
    DO  I=1,L
      FT(I,K,M)=-F(I,M,K)
    END DO
  END DO
END IF
DO  J=2,M-1,2
  JBY2=J/2
  DO  K=1,N
    DO  I=1,L
      FT(I,K,J)   =  F(I,J+1,K)*C1(JBY2)+F(I,J,K)*C2(JBY2)
      FT(I,K,J+1) = -F(I,J+1,K)*C2(JBY2)+F(I,J,K)*C1(JBY2)
    END DO
  END DO
END DO
RETURN
END SUBROUTINE VSCSQ1


SUBROUTINE VSCSQI(N,C1,C2,C3,C4,WSAVE)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

!      ENTRY VSSNQI(N,C1,C2,C3,C4,WSAVE)


INTEGER                       :: N
REAL(EB)  C3(N)
REAL(EB)   C4(N)
REAL(EB)   WSAVE(N+15)
REAL(EB)   C1(N),C2(N)
INTEGER :: I
REAL(EB) :: DX,C,S

DX=PI/N
SCALE=SQRT(.5_EB)

!     GENERATE A(I)+-B(I)

DO  I=1,(N-1)/2
  C=COS(I*DX)
  S=SIN(I*DX)
  C1(I)=.5_EB*(S+C)
  C2(I)=.5_EB*(C-S)
END DO

DX=PI/(2*N)
DO  I=1,N
  C=COS((I-.5_EB)*DX)
  S=SIN((I-.5_EB)*DX)
  C3(I)=SCALE*(C+S)
  C4(I)=SCALE*(C-S)
END DO

!     INITIALIZE VRFFTPK ROUTINE

CALL VRFFTI(N,WSAVE)
RETURN
END SUBROUTINE VSCSQI


SUBROUTINE VSINT(X,L,M,N,LDIMX,XT,C,WSAVE)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989


INTEGER                   :: L
INTEGER                   :: M
INTEGER                       :: N
INTEGER                       :: LDIMX
REAL(EB)  XT(LDIMX*N,M+1)
REAL(EB)  C(M/2)
REAL(EB)  WSAVE(M+16)
REAL(EB)  X(LDIMX*N,M)
INTEGER MODM,MP1,MS2,I,J,JC
REAL(EB) :: SQRT2I

!     TPOSE = .TRUE. IF TRANSFORMING SECOND INDEX (SO TRANSPOSE)
!           = .FALSE. IF TRANSFORMING THIRD INDEX

SQRT2I=SQRT(.5_EB)
MODM = MOD(M,2)
MP1 = M+1
MS2 = M/2
DO  I=1,LDIMX*N
  XT(I,1) = 0._EB
  XT(I,M+1)=0._EB
END DO
!                         ZERO OUT LAST PLANE BECAUSE FOR SINE
!                         IT DOESN'T GET DONE IN FSH02
IF (TPOSE) THEN
  CALL VSINT1(L,M,N,LDIMX,MODM,MS2,C,XT,X)
ELSE
  DO  J=1,MS2
    JC=M+1-J
    DO  I=1,LDIMX*N
      XT(I,J+1) = (X(I,J)-X(I,JC))+C(J)*(X(I,J)+X(I,JC))
      XT(I,JC+1) = C(J)*(X(I,J)+X(I,JC))-(X(I,J)-X(I,JC))
    END DO
  END DO
  IF (MODM == 0) GO TO 405
  DO  I=1,LDIMX*N
    XT(I,MS2+2) = 4._EB*X(I,MS2+1)
  END DO
405    CONTINUE
END IF
CALL VRFFTF (LDIMX*N,MP1,XT,LDIMX*N,X,WSAVE)
SCALE=SCALE*SQRT2I
IF (OUTARY) THEN
  CALL VSINTA(M,LDIMX*N,MODM,X,XT)
ELSE
  CALL VSINTA(M,LDIMX*N,MODM,XT,X)
END IF
RETURN
END SUBROUTINE VSINT


SUBROUTINE VSINT1(L,M,N,LDIMX,MODM,MS2,C,XT,X)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989



INTEGER                       :: L
INTEGER                       :: M
INTEGER                       :: N
INTEGER                   :: LDIMX
INTEGER                   :: MODM
INTEGER                   :: MS2
REAL(EB)   C(*)
REAL(EB)   X(LDIMX,M,N)
REAL(EB)   XT(LDIMX,N,M+1)
INTEGER J, JC, K, I

DO  J=1,MS2
  JC = M+1-J
  DO  K=1,N
    DO  I=1,L
      XT(I,K,J+1) = (X(I,J,K)-X(I,JC,K))+C(J)*(X(I,J,K)+X(I,JC,K))
      XT(I,K,JC+1) = C(J)*(X(I,J,K)+X(I,JC,K))-(X(I,J,K)-X(I,JC,K))
    END DO
  END DO
END DO
IF (MODM == 0) GO TO 405
DO  K=1,N
  DO  I=1,L
    XT(I,K,MS2+2) = 4._EB*X(I,MS2+1,K)
  END DO
END DO
405 CONTINUE
RETURN
END SUBROUTINE VSINT1


SUBROUTINE VSINTA(M,LDIMX,MODM,X,XT)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989



INTEGER                       :: M
INTEGER                       :: LDIMX
INTEGER                   :: MODM
REAL(EB)   X(LDIMX,M),  XT(LDIMX,M+1)
INTEGER :: I,J

DO  I=1,LDIMX
  X(I,1) = 0.5_EB*XT(I,1)
END DO
DO  J=2,M-1,2
  DO  I=1,LDIMX
    X(I,J) = -XT(I,J+1)
    X(I,J+1) = X(I,J-1)+XT(I,J)
  END DO
END DO
IF (MODM /= 0) GO TO 900
DO  I=1,LDIMX
  X(I,M) = -XT(I,M+1)
END DO
900 CONTINUE
RETURN
END SUBROUTINE VSINTA


SUBROUTINE VSINTI(N,C,WSAVE)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

!     *******   NOTE   ******

!     SPECIAL DOCUMENTATION FOR USE OF INTERNAL PARAMETER TPOSE.  IF
!     TPOSE = .TRUE., THEN ALL PREPROCESSING ROUTINES TRANSPOSE SECOND
!                     AND THIRD INDICES IN PREPARATION FOR TRANSFORM OF
!                     ORIGINAL SECOND INDEX.  IF
!     TPOSE = .FALSE., THEN ALL PREPROCESSING ROUTINES PREPARE ONLY FOR
!                      TRANSFORM OF THIRD INDEX.  TO USE IN THIS CASE,
!                      CALL ROUTINE WITH M AND N INTERCHANGED.



INTEGER                       :: N
REAL(EB)   C(N/2),WSAVE (N+16)
INTEGER :: NP1,NS2,K
REAL(EB) :: DT

!                               INITIALIZE NOCOPY AND TPOSE TO DEFAULT
!                               VALUES

NOCOPY = .FALSE.
TPOSE = .TRUE.

IF (N <= 1) RETURN
NP1 = N+1
NS2 = N/2
DT = PI/REAL(NP1,EB)
DO  K=1,NS2
  C(K) = 2._EB*SIN(K*DT)
END DO
CALL VRFFTI (NP1,WSAVE)
RETURN
END SUBROUTINE VSINTI


SUBROUTINE VSRFTB(F,L,M,N,LDIMF,FT,WSAVE)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989


!     *******   NOTE   ******

!     SPECIAL DOCUMENTATION FOR USE OF INTERNAL PARAMETER TPOSE.  IF
!     TPOSE = .TRUE., THEN ALL PREPROCESSING ROUTINES TRANSPOSE SECOND
!                     AND THIRD INDICES IN PREPARATION FOR TRANSFORM OF
!                     ORIGINAL SECOND INDEX.  IF
!     TPOSE = .FALSE., THEN ALL PREPROCESSING ROUTINES PREPARE ONLY FOR
!                      TRANSFORM OF THIRD INDEX.  TO USE IN THIS CASE,
!                      CALL ROUTINE WITH M AND N INTERCHANGED.


!     VSFFTPK, VERSION 2, JUNE 1988


INTEGER                       :: L
INTEGER                       :: M
INTEGER                       :: N
INTEGER                   :: LDIMF
REAL(EB)   FT(LDIMF,N,M)
REAL(EB)   WSAVE(M+15)
REAL(EB)   F(LDIMF,M,N)
INTEGER :: I,K,J

!     TPOSE = .TRUE. IF TRANSFORMING SECOND INDEX (SO TRANSPOSE)
!           = .FALSE. IF TRANSFORMING THIRD INDEX

IF (TPOSE) THEN

!        RE-ORDER INPUT

  DO  K=1,N
    DO  J=1,M
      DO  I=1,L
        FT(I,K,J)=F(I,J,K)
      END DO
    END DO
  END DO

!        REAL(EB), PERIODIC TRANSFORM

  CALL VRFFTB(LDIMF*N,M,FT,LDIMF*N,F,WSAVE)
  OUTARY=.NOT.OUTARY
  IF (.NOT.NOCOPY) THEN
    CALL VSRTB1(M,N*LDIMF,F,FT)
  END IF
ELSE
  CALL VRFFTB(LDIMF*N,M,F,LDIMF*N,FT,WSAVE)
END IF
RETURN
END SUBROUTINE VSRFTB


SUBROUTINE VSRFTF(F,L,M,N,LDIMF,FT,WSAVE)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

INTEGER                   :: L
INTEGER                   :: M
INTEGER                   :: N
INTEGER                   :: LDIMF
REAL(EB)  FT(LDIMF,N,M)
REAL(EB)   WSAVE(M+15)
REAL(EB)   F(LDIMF,N,M)

!     TPOSE = .TRUE. IF TRANSFORMING SECOND INDEX (SO TRANSPOSE)
!           = .FALSE. IF TRANSFORMING THIRD INDEX

IF (TPOSE) THEN

!        RE-ORDER INPUT

  CALL VSRTF1(L,M,N,LDIMF,F,FT)

!        REAL(EB), PERIODIC TRANSFORM

  CALL VRFFTF(LDIMF*N,M,FT,LDIMF*N,F,WSAVE)
  OUTARY=.NOT.OUTARY
  IF (.NOT.NOCOPY) THEN
    CALL VSRTB1(M,N*LDIMF,F,FT)
  END IF
ELSE
  CALL VRFFTF(LDIMF*N,M,F,LDIMF*N,FT,WSAVE)
END IF
RETURN
END SUBROUTINE VSRFTF


SUBROUTINE VSRFTI(N,WSAVE)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989


!     *******   NOTE   ******

!     SPECIAL DOCUMENTATION FOR USE OF INTERNAL PARAMETER TPOSE.  IF
!     TPOSE = .TRUE., THEN ALL PREPROCESSING ROUTINES TRANSPOSE SECOND
!                     AND THIRD INDICES IN PREPARATION FOR TRANSFORM OF
!                     ORIGINAL SECOND INDEX.  IF
!     TPOSE = .FALSE., THEN ALL PREPROCESSING ROUTINES PREPARE ONLY FOR
!                      TRANSFORM OF THIRD INDEX.  TO USE IN THIS CASE,
!                      CALL ROUTINE WITH M AND N INTERCHANGED.



INTEGER                   :: N
REAL(EB)   WSAVE(N+15)

!     INITIALIZE VRFFTPK ROUTINE

CALL VRFFTI(N,WSAVE)
RETURN
END SUBROUTINE VSRFTI


SUBROUTINE VSRTB1(M,LDIMF,F,FT)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989



INTEGER                       :: M
INTEGER                       :: LDIMF
REAL(EB)   F(LDIMF,M),  FT(LDIMF,M)
INTEGER :: I,J

DO  J=1,M
  DO  I=1,LDIMF
    F(I,J)=FT(I,J)
  END DO
END DO
RETURN
END SUBROUTINE VSRTB1


SUBROUTINE VSRTF1(L,M,N,LDIMF,F,FT)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989



INTEGER                       :: L
INTEGER                       :: M
INTEGER                       :: N
INTEGER                   :: LDIMF
REAL(EB)   FT(LDIMF,N,M)
REAL(EB)   F(LDIMF,M,N)
INTEGER :: I,K,J

DO  K=1,N
  DO  J=1,M
    DO  I=1,L
      FT(I,K,J)=F(I,J,K)
    END DO
  END DO
END DO
RETURN
END SUBROUTINE VSRTF1


SUBROUTINE VSSINB(F,L,M,N,LDIMF,FT,C1,C2,WORK)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

INTEGER                   :: L
INTEGER                       :: M
INTEGER                       :: N
INTEGER                       :: LDIMF
REAL(EB) FT(LDIMF*N,M)
REAL(EB) C1(M)
REAL(EB)  C2(M)
REAL(EB)  WORK(M+15)
REAL(EB)   F(LDIMF*N,M)
INTEGER :: I,J

!     TPOSE = .TRUE. IF TRANSFORMING SECOND INDEX (SO TRANSPOSE)
!           = .FALSE. IF TRANSFORMING THIRD INDEX


!     PREPROCESSING

IF (TPOSE) THEN
  CALL VSSNB1(L,M,N,LDIMF,F,FT,C1,C2)
ELSE
  DO  J=2,M
    DO  I=1,LDIMF*N
      FT(I,J)=C1(J)*F(I,J-1)-C2(J)*F(I,M-J+1)
    END DO
  END DO
DO  I=1,LDIMF*N
  FT(I,1) = .5_EB*F(I,M)
END DO
END IF

!     REAL(EB),PERIODIC ANALYSIS

CALL VRFFTF(LDIMF*N,M,FT,LDIMF*N,F,WORK)

!     POSTPROCESSING

SCALE=SQRT(2.0_EB)*SCALE
IF (OUTARY) THEN
  CALL VSSNBA(M,N*LDIMF,F,FT)
ELSE
  CALL VSSNBA(M,N*LDIMF,FT,F)
END IF
RETURN
END SUBROUTINE VSSINB


SUBROUTINE VSSINF(F,L,M,N,LDIMF,FT,C1,C2,WORK)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989



INTEGER                   :: L
INTEGER                       :: M
INTEGER                       :: N
INTEGER                       :: LDIMF
REAL(EB)  FT(LDIMF*N,M)
REAL(EB) C1(M)
REAL(EB) C2(M)
REAL(EB)  WORK(M+15)
REAL(EB)   F(LDIMF*N,M)
INTEGER :: I,J

!     TPOSE = .TRUE. IF TRANSFORMING SECOND INDEX (SO TRANSPOSE)
!           = .FALSE. IF TRANSFORMING THIRD INDEX


!     PREPROCESSING

SCALE=SQRT(2.0_EB)*SCALE
IF (TPOSE) THEN
  CALL VSSNF1(L,M,N,LDIMF,F,FT)
ELSE
  DO  I=1,LDIMF*N
    FT(I,1)=F(I,1)
  END DO
IF (MOD(M,2)==0) THEN
  DO  I=1,LDIMF*N
    FT(I,M)=-F(I,M)
  END DO
END IF
DO  J=2,M-1,2
  DO  I=1,LDIMF*N
    FT(I,J)   = 0.5_EB*(F(I,J+1)-F(I,J))
    FT(I,J+1) =-0.5_EB*(F(I,J+1)+F(I,J))
  END DO
END DO
END IF

!     REAL(EB),PERIODIC SYNTHESIS

CALL VRFFTB(LDIMF*N,M,FT,LDIMF*N,F,WORK)

!     POSTPROCESSING

IF (OUTARY) THEN
  CALL VSSNFA(M,N*LDIMF,F,FT,C1,C2)
ELSE
  CALL VSSNFA(M,N*LDIMF,FT,F,C1,C2)
END IF
RETURN
END SUBROUTINE VSSINF


SUBROUTINE VSSINQ(F,L,M,N,LDIMF,FT,C1,C2,C3,C4,WORK)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

INTEGER                   :: L
INTEGER                       :: M
INTEGER                       :: N
INTEGER                       :: LDIMF
REAL(EB) FT(LDIMF*N,M)
REAL(EB)  C1(M)
REAL(EB)  C2(M)
REAL(EB)  C3(M)
REAL(EB)  C4(M)
REAL(EB) WORK(M+15)
REAL(EB)   F(LDIMF*N,M)
INTEGER :: I,J,JBY2

!     TPOSE = .TRUE. IF TRANSFORMING SECOND INDEX (SO TRANSPOSE)
!           = .FALSE. IF TRANSFORMING THIRD INDEX


!     PREPROCESSING

IF (TPOSE) THEN
  CALL VSSNQ1(L,M,N,LDIMF,F,FT,C1,C2)
ELSE
  DO  I=1,LDIMF*N
    FT(I,1)=F(I,1)
  END DO
IF (MOD(M,2)==0) THEN
  DO  I=1,LDIMF*N
    FT(I,M)=F(I,M)
  END DO
END IF
DO  J=2,M-1,2
  JBY2=J/2
  DO  I=1,LDIMF*N
    FT(I,J)   =   F(I,J+1)*C1(JBY2)-F(I,J)*C2(JBY2)
    FT(I,J+1) = F(I,J+1)*(-C2(JBY2))-F(I,J)*C1(JBY2)
  END DO
END DO
END IF

!     REAL(EB),PERIODIC SYNTHESIS

CALL VRFFTB(LDIMF*N,M,FT,LDIMF*N,F,WORK)

!     POSTPROCESSING

IF (OUTARY) THEN
  DO  J=1,M
    DO  I=1,LDIMF*N
      F(I,J)=C3(J)*FT(I,J)-C4(J)*FT(I,M+1-J)
    END DO
  END DO
ELSE
DO  J=1,M
  DO  I=1,LDIMF*N
    FT(I,J)=C3(J)*F(I,J)-C4(J)*F(I,M+1-J)
  END DO
END DO
END IF
RETURN
END SUBROUTINE VSSINQ


SUBROUTINE VSSNB1(L,M,N,LDIMF,F,FT,C1,C2)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989



INTEGER                       :: L
INTEGER                       :: M
INTEGER                       :: N
INTEGER                   :: LDIMF
REAL(EB)  FT(LDIMF,N,M)
REAL(EB)  C1(M)
REAL(EB) C2(M)
REAL(EB)   F(LDIMF,M,N)
INTEGER :: I,J,K

DO  J=2,M
  DO  K=1,N
    DO  I=1,L
      FT(I,K,J)=C1(J)*F(I,J-1,K)-C2(J)*F(I,M-J+1,K)
    END DO
  END DO
END DO
DO  K=1,N
  DO  I=1,L
    FT(I,K,1) = .5_EB*F(I,M,K)
  END DO
END DO
RETURN
END SUBROUTINE VSSNB1


SUBROUTINE VSSNBA(M,LDIMF,F,FT)
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989



INTEGER                       :: M
INTEGER                       :: LDIMF
INTEGER :: I,J

REAL(EB)   F(LDIMF,M),  FT(LDIMF,M)

DO  I=1,LDIMF
  F(I,1) = FT(I,1)
END DO
IF (MOD(M,2)==0) THEN
  DO  I=1,LDIMF
    F(I,M) = -FT(I,M)
  END DO
END IF
DO  J=2,M-1,2
  DO  I=1,LDIMF
    F(I,J)    = -(FT(I,J+1)+FT(I,J))
    F(I,J+1)  =   FT(I,J)-FT(I,J+1)
  END DO
END DO

RETURN
END SUBROUTINE VSSNBA


SUBROUTINE VSSNF1(L,M,N,LDIMF,F,FT)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989



INTEGER                       :: L
INTEGER                       :: M
INTEGER                       :: N
INTEGER                   :: LDIMF
REAL(EB)   FT(LDIMF,N,M)
REAL(EB)   F(LDIMF,M,N)
INTEGER :: I,K,J

DO  K=1,N
  DO  I=1,L
    FT(I,K,1)=F(I,1,K)
  END DO
END DO
IF (MOD(M,2)==0) THEN
  DO  K=1,N
    DO  I=1,L
      FT(I,K,M)=-F(I,M,K)
    END DO
  END DO
END IF
DO  J=2,M-1,2
  DO  K=1,N
    DO  I=1,L
      FT(I,K,J)   = 0.5_EB*(F(I,J+1,K)-F(I,J,K))
      FT(I,K,J+1) =-0.5_EB*(F(I,J+1,K)+F(I,J,K))
    END DO
  END DO
END DO
RETURN
END SUBROUTINE VSSNF1


SUBROUTINE VSSNFA(M,LDIMF,F,FT,C1,C2)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989



INTEGER                       :: M
INTEGER                       :: LDIMF
REAL(EB)   FT(LDIMF,M)
REAL(EB)    C1(M)
REAL(EB)   F(LDIMF,M), C2(M)
INTEGER :: I, J

DO  J=2,M
  DO  I=1,LDIMF
    F(I,J-1)=C1(J)*FT(I,J)+C2(J)*FT(I,M+2-J)
  END DO
END DO
DO  I=1,LDIMF
  F(I,M)=FT(I,1)
END DO
RETURN
END SUBROUTINE VSSNFA


SUBROUTINE VSSNQ1(L,M,N,LDIMF,F,FT,C1,C2)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989



INTEGER                       :: L
INTEGER                       :: M
INTEGER                       :: N
INTEGER                   :: LDIMF
REAL(EB) FT(LDIMF,N,M)
REAL(EB)  C1(M)
REAL(EB)  C2(M)
REAL(EB)   F(LDIMF,M,N)
INTEGER :: I,K,J,JBY2

DO  K=1,N
  DO  I=1,L
    FT(I,K,1)=F(I,1,K)
  END DO
END DO
IF (MOD(M,2)==0) THEN
  DO  K=1,N
    DO  I=1,L
      FT(I,K,M)=F(I,M,K)
    END DO
  END DO
END IF
DO  J=2,M-1,2
  JBY2=J/2
  DO  K=1,N
    DO  I=1,L
      FT(I,K,J)   =   F(I,J+1,K)*C1(JBY2)-F(I,J,K)*C2(JBY2)
      FT(I,K,J+1) = -(F(I,J+1,K)*C2(JBY2)+F(I,J,K)*C1(JBY2))
    END DO
  END DO
END DO
RETURN
END SUBROUTINE VSSNQ1


SUBROUTINE VRADF2 (MP,IDO,L1,CC,MDIMC,CH,MDIMCH,WA1)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989

INTEGER :: MDIMCH
INTEGER  :: MP, IDO, L1
INTEGER  :: MDIMC
REAL(EB) CH(MDIMCH,IDO,2,L1), CC(MDIMC,IDO,L1,2),  WA1(IDO)
INTEGER :: I,K,M,IDP2,IC

DO  K=1,L1
  DO  M=1,MP
    CH(M,1,1,K) = CC(M,1,K,1)+CC(M,1,K,2)
    CH(M,IDO,2,K) = CC(M,1,K,1)-CC(M,1,K,2)
  END DO
END DO
IF (IDO-2<0) GO TO 107
IF (IDO-2==0) GO TO 105
IF (IDO-2>0) GO TO 102
102 IDP2 = IDO+2
DO  K=1,L1
  DO  I=3,IDO,2
    IC = IDP2-I
    DO  M=1,MP
      CH(M,I,1,K) = CC(M,I,K,1)+(WA1(I-2)*CC(M,I,K,2)- WA1(I-1)*CC(M,I-1,K,2))
      CH(M,IC,2,K) = (WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*  &
          CC(M,I-1,K,2))-CC(M,I,K,1)
      CH(M,I-1,1,K) = CC(M,I-1,K,1)+(WA1(I-2)*CC(M,I-1,K,2)+  &
          WA1(I-1)*CC(M,I,K,2))
      CH(M,IC-1,2,K) = CC(M,I-1,K,1)-(WA1(I-2)*CC(M,I-1,K,2)+  &
          WA1(I-1)*CC(M,I,K,2))
    END DO
  END DO
END DO
IF (MOD(IDO,2) == 1) RETURN
105 DO  K=1,L1
  DO  M=1,MP
    CH(M,1,2,K) = -CC(M,IDO,K,2)
    CH(M,IDO,1,K) = CC(M,IDO,K,1)
  END DO
END DO
107 RETURN
END SUBROUTINE VRADF2


SUBROUTINE VRADF3 (MP,IDO,L1,CC,MDIMC,CH,MDIMCH,WA1,WA2)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989


INTEGER                       :: MP, IDO, L1
INTEGER                   :: MDIMC
INTEGER                   :: MDIMCH
REAL(EB)    CH(MDIMCh,IDO,3,L1)  ,CC(MDIMC,IDO,L1,3)
REAL(EB)  WA1(IDO)
REAL(EB)   WA2(IDO)
INTEGER :: I,M,K,IDP2,IC
REAL(EB) :: ARG,TAUR,TAUI

ARG=2._EB*PI/3._EB
TAUR=COS(ARG)
TAUI=SIN(ARG)
DO  M=1,MP
   DO  K=1,L1
       CH(M,1,1,K) = CC(M,1,K,1)+(CC(M,1,K,2)+CC(M,1,K,3))
       CH(M,1,3,K) = TAUI*(CC(M,1,K,3)-CC(M,1,K,2))
       CH(M,IDO,2,K) = CC(M,1,K,1)+TAUR* (CC(M,1,K,2)+CC(M,1,K,3))
   END DO
END DO
IF (IDO == 1) RETURN
IDP2 = IDO+2
DO  K=1,L1
   DO  I=3,IDO,2
      IC = IDP2-I
      DO  M=1,MP
         CH(M,I-1,1,K) = CC(M,I-1,K,1)+((WA1(I-2)*CC(M,I-1,K,2)+  &
               WA1(I-1)*CC(M,I,K,2))+(WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)* CC(M,I,K,3)))
         CH(M,I,1,K) = CC(M,I,K,1)+((WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*  &
               CC(M,I-1,K,2))+(WA2(I-2)*CC(M,I,K,3)-WA2(I-1)* CC(M,I-1,K,3)))
         CH(M,I-1,3,K) = (CC(M,I-1,K,1)+TAUR*((WA1(I-2)*  &
               CC(M,I-1,K,2)+WA1(I-1)*CC(M,I,K,2))+(WA2(I-2)*  &
               CC(M,I-1,K,3)+WA2(I-1)*CC(M,I,K,3))))+(TAUI*((WA1(I-2)*  &
               CC(M,I,K,2)-WA1(I-1)*CC(M,I-1,K,2))-(WA2(I-2)*  &
               CC(M,I,K,3)-WA2(I-1)*CC(M,I-1,K,3))))
         CH(M,IC-1,2,K) = (CC(M,I-1,K,1)+TAUR*((WA1(I-2)*  &
               CC(M,I-1,K,2)+WA1(I-1)*CC(M,I,K,2))+(WA2(I-2)*  &
               CC(M,I-1,K,3)+WA2(I-1)*CC(M,I,K,3))))-(TAUI*((WA1(I-2)*  &
               CC(M,I,K,2)-WA1(I-1)*CC(M,I-1,K,2))-(WA2(I-2)*  &
               CC(M,I,K,3)-WA2(I-1)*CC(M,I-1,K,3))))
         CH(M,I,3,K) = (CC(M,I,K,1)+TAUR*((WA1(I-2)*CC(M,I,K,2)-  &
               WA1(I-1)*CC(M,I-1,K,2))+(WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*  &
               CC(M,I-1,K,3))))+(TAUI*((WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*  &
               CC(M,I,K,3))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)* CC(M,I,K,2))))
         CH(M,IC,2,K) = (TAUI*((WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*  &
               CC(M,I,K,3))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*  &
               CC(M,I,K,2))))-(CC(M,I,K,1)+TAUR*((WA1(I-2)*CC(M,I,K,2)-  &
               WA1(I-1)*CC(M,I-1,K,2))+(WA2(I-2)*CC(M,I,K,3)-WA2(I-1)* CC(M,I-1,K,3))))
      END DO
   END DO
END DO
RETURN
END SUBROUTINE VRADF3


SUBROUTINE VRADF4 (MP,IDO,L1,CC,MDIMC,CH,MDIMCH,WA1,WA2,WA3)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989


INTEGER                       :: MP, IDO, L1
INTEGER                   :: MDIMC
INTEGER                   :: MDIMCH
REAL(EB)  CC(MDIMC,IDO,L1,4)   ,CH(MDIMCH,IDO,4,L1)
REAL(EB)    WA1(IDO)
REAL(EB)   WA2(IDO)
REAL(EB)   WA3(IDO)
INTEGER :: I,M,K,IDP2,IC
REAL(EB) :: HSQT2

HSQT2=SQRT(2.0_EB)*0.5_EB
DO  M=1,MP
   DO  K=1,L1
      CH(M,1,1,K) = (CC(M,1,K,2)+CC(M,1,K,4)) +(CC(M,1,K,1)+CC(M,1,K,3))
      CH(M,IDO,4,K) = (CC(M,1,K,1)+CC(M,1,K,3)) -(CC(M,1,K,2)+CC(M,1,K,4))
      CH(M,IDO,2,K) = CC(M,1,K,1)-CC(M,1,K,3)
      CH(M,1,3,K) = CC(M,1,K,4)-CC(M,1,K,2)
   END DO
END DO

IF (IDO < 2) RETURN
IF (IDO > 2) THEN
   IDP2 = IDO+2
   DO  M=1,MP
      DO  I=3,IDO,2
         IC = IDP2-I
         DO  K=1,L1
            CH(M,I-1,1,K) = ((WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*  &
                CC(M,I,K,2))+(WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*  &
                CC(M,I,K,4)))+(CC(M,I-1,K,1)+(WA2(I-2)*CC(M,I-1,K,3)+  &
                WA2(I-1)*CC(M,I,K,3)))
            CH(M,IC-1,4,K) = (CC(M,I-1,K,1)+(WA2(I-2)*CC(M,I-1,K,3)+  &
                WA2(I-1)*CC(M,I,K,3)))-((WA1(I-2)*CC(M,I-1,K,2)+  &
                WA1(I-1)*CC(M,I,K,2))+(WA3(I-2)*CC(M,I-1,K,4)+ WA3(I-1)*CC(M,I,K,4)))
            CH(M,I,1,K) = ((WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*  &
                CC(M,I-1,K,2))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*  &
                CC(M,I-1,K,4)))+(CC(M,I,K,1)+(WA2(I-2)*CC(M,I,K,3)-  &
                WA2(I-1)*CC(M,I-1,K,3)))
            CH(M,IC,4,K) = ((WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*  &
                CC(M,I-1,K,2))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*  &
                CC(M,I-1,K,4)))-(CC(M,I,K,1)+(WA2(I-2)*CC(M,I,K,3)-  &
                WA2(I-1)*CC(M,I-1,K,3)))
            CH(M,I-1,3,K) = ((WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*  &
                CC(M,I-1,K,2))-(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*  &
                CC(M,I-1,K,4)))+(CC(M,I-1,K,1)-(WA2(I-2)*CC(M,I-1,K,3)+  &
                WA2(I-1)*CC(M,I,K,3)))
            CH(M,IC-1,2,K) = (CC(M,I-1,K,1)-(WA2(I-2)*CC(M,I-1,K,3)+  &
                WA2(I-1)*CC(M,I,K,3)))-((WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*  &
                CC(M,I-1,K,2))-(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)* CC(M,I-1,K,4)))
            CH(M,I,3,K) = ((WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*  &
                CC(M,I,K,4))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*  &
                CC(M,I,K,2)))+(CC(M,I,K,1)-(WA2(I-2)*CC(M,I,K,3)-  &
                WA2(I-1)*CC(M,I-1,K,3)))
            CH(M,IC,2,K) = ((WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*  &
                CC(M,I,K,4))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*  &
                CC(M,I,K,2)))-(CC(M,I,K,1)-(WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*  &
                CC(M,I-1,K,3)))
         END DO
      END DO
   END DO
   IF (MOD(IDO,2) == 1) RETURN
ENDIF

DO  M=1,MP
   DO  K=1,L1
      CH(M,IDO,1,K) = (HSQT2*(CC(M,IDO,K,2)-CC(M,IDO,K,4)))+ CC(M,IDO,K,1)
      CH(M,IDO,3,K) = CC(M,IDO,K,1)-(HSQT2*(CC(M,IDO,K,2)- CC(M,IDO,K,4)))
      CH(M,1,2,K) = (-HSQT2*(CC(M,IDO,K,2)+CC(M,IDO,K,4)))- CC(M,IDO,K,3)
      CH(M,1,4,K) = (-HSQT2*(CC(M,IDO,K,2)+CC(M,IDO,K,4)))+ CC(M,IDO,K,3)
   END DO
END DO

END SUBROUTINE VRADF4


SUBROUTINE VRADF5 (MP,IDO,L1,CC,MDIMC,CH,MDIMCH,WA1,WA2,WA3,WA4)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989


INTEGER                       :: MP
INTEGER                      :: IDO,L1
INTEGER                   :: MDIMC
INTEGER                   :: MDIMCH
REAL(EB)   CC(MDIMC,IDO,L1,5)    ,CH(MDIMCH,IDO,5,L1)
REAL(EB)   WA1(IDO)
REAL(EB)   WA2(IDO)
REAL(EB)   WA3(IDO)
REAL(EB)    WA4(IDO)
INTEGER :: K,M,IDP2,I,IC
REAL(EB) :: ARG,TR11,TI11,TR12,TI12

ARG=2._EB*PI/5._EB
TR11=COS(ARG)
TI11=SIN(ARG)
TR12=COS(2._EB*ARG)
TI12=SIN(2._EB*ARG)
DO  K=1,L1
  DO  M=1,MP
    CH(M,1,1,K) = CC(M,1,K,1)+(CC(M,1,K,5)+CC(M,1,K,2))+  &
        (CC(M,1,K,4)+CC(M,1,K,3))
    CH(M,IDO,2,K) = CC(M,1,K,1)+TR11*(CC(M,1,K,5)+CC(M,1,K,2))+  &
        TR12*(CC(M,1,K,4)+CC(M,1,K,3))
    CH(M,1,3,K) = TI11*(CC(M,1,K,5)-CC(M,1,K,2))+TI12*  &
        (CC(M,1,K,4)-CC(M,1,K,3))
    CH(M,IDO,4,K) = CC(M,1,K,1)+TR12*(CC(M,1,K,5)+CC(M,1,K,2))+  &
        TR11*(CC(M,1,K,4)+CC(M,1,K,3))
    CH(M,1,5,K) = TI12*(CC(M,1,K,5)-CC(M,1,K,2))-TI11*  &
        (CC(M,1,K,4)-CC(M,1,K,3))
  END DO
END DO
IF (IDO == 1) RETURN
IDP2 = IDO+2
DO  K=1,L1
  DO  I=3,IDO,2
    IC = IDP2-I
    DO  M=1,MP
      CH(M,I-1,1,K) = CC(M,I-1,K,1)+((WA1(I-2)*CC(M,I-1,K,2)+  &
          WA1(I-1)*CC(M,I,K,2))+(WA4(I-2)*CC(M,I-1,K,5)+WA4(I-1)*  &
          CC(M,I,K,5)))+((WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*  &
          CC(M,I,K,3))+(WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*CC(M,I,K,4)))
      CH(M,I,1,K) = CC(M,I,K,1)+((WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*  &
          CC(M,I-1,K,2))+(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*  &
          CC(M,I-1,K,5)))+((WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*  &
          CC(M,I-1,K,3))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)* CC(M,I-1,K,4)))
      CH(M,I-1,3,K) = CC(M,I-1,K,1)+TR11*  &
          ( WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*CC(M,I,K,2)  &
          +WA4(I-2)*CC(M,I-1,K,5)+WA4(I-1)*CC(M,I,K,5))+TR12*  &
          ( WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*CC(M,I,K,3)  &
          +WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*CC(M,I,K,4))+TI11*  &
          ( WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*CC(M,I-1,K,2)  &
          -(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*CC(M,I-1,K,5)))+TI12*  &
          ( WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*CC(M,I-1,K,3)  &
          -(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*CC(M,I-1,K,4)))
      CH(M,IC-1,2,K) = CC(M,I-1,K,1)+TR11*  &
          ( WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*CC(M,I,K,2)  &
          +WA4(I-2)*CC(M,I-1,K,5)+WA4(I-1)*CC(M,I,K,5))+TR12*  &
          ( WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*CC(M,I,K,3)  &
          +WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*CC(M,I,K,4))-(TI11*  &
          ( WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*CC(M,I-1,K,2)  &
          -(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*CC(M,I-1,K,5)))+TI12*  &
          ( WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*CC(M,I-1,K,3)  &
          -(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*CC(M,I-1,K,4))))
      CH(M,I,3,K) = (CC(M,I,K,1)+TR11*((WA1(I-2)*CC(M,I,K,2)-  &
          WA1(I-1)*CC(M,I-1,K,2))+(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*  &
          CC(M,I-1,K,5)))+TR12*((WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*  &
          CC(M,I-1,K,3))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*  &
          CC(M,I-1,K,4))))+(TI11*((WA4(I-2)*CC(M,I-1,K,5)+  &
          WA4(I-1)*CC(M,I,K,5))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*  &
          CC(M,I,K,2)))+TI12*((WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*  &
          CC(M,I,K,4))-(WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)* CC(M,I,K,3))))
      CH(M,IC,2,K) = (TI11*((WA4(I-2)*CC(M,I-1,K,5)+WA4(I-1)*  &
          CC(M,I,K,5))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*  &
          CC(M,I,K,2)))+TI12*((WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*  &
          CC(M,I,K,4))-(WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*  &
          CC(M,I,K,3))))-(CC(M,I,K,1)+TR11*((WA1(I-2)*CC(M,I,K,2)-  &
          WA1(I-1)*CC(M,I-1,K,2))+(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*  &
          CC(M,I-1,K,5)))+TR12*((WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*  &
          CC(M,I-1,K,3))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)* CC(M,I-1,K,4))))
      CH(M,I-1,5,K) = (CC(M,I-1,K,1)+TR12*((WA1(I-2)*  &
          CC(M,I-1,K,2)+WA1(I-1)*CC(M,I,K,2))+(WA4(I-2)*  &
          CC(M,I-1,K,5)+WA4(I-1)*CC(M,I,K,5)))+TR11*((WA2(I-2)*  &
          CC(M,I-1,K,3)+WA2(I-1)*CC(M,I,K,3))+(WA3(I-2)*  &
          CC(M,I-1,K,4)+WA3(I-1)*CC(M,I,K,4))))+(TI12*((WA1(I-2)*  &
          CC(M,I,K,2)-WA1(I-1)*CC(M,I-1,K,2))-(WA4(I-2)*CC(M,I,K,5)-  &
          WA4(I-1)*CC(M,I-1,K,5)))-TI11*((WA2(I-2)*CC(M,I,K,3)-  &
          WA2(I-1)*CC(M,I-1,K,3))-(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)* CC(M,I-1,K,4))))
      CH(M,IC-1,4,K) = (CC(M,I-1,K,1)+TR12*((WA1(I-2)*  &
          CC(M,I-1,K,2)+WA1(I-1)*CC(M,I,K,2))+(WA4(I-2)*  &
          CC(M,I-1,K,5)+WA4(I-1)*CC(M,I,K,5)))+TR11*((WA2(I-2)*  &
          CC(M,I-1,K,3)+WA2(I-1)*CC(M,I,K,3))+(WA3(I-2)*  &
          CC(M,I-1,K,4)+WA3(I-1)*CC(M,I,K,4))))-(TI12*((WA1(I-2)*  &
          CC(M,I,K,2)-WA1(I-1)*CC(M,I-1,K,2))-(WA4(I-2)*CC(M,I,K,5)-  &
          WA4(I-1)*CC(M,I-1,K,5)))-TI11*((WA2(I-2)*CC(M,I,K,3)-  &
          WA2(I-1)*CC(M,I-1,K,3))-(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)* CC(M,I-1,K,4))))
      CH(M,I,5,K) = (CC(M,I,K,1)+TR12*((WA1(I-2)*CC(M,I,K,2)-  &
          WA1(I-1)*CC(M,I-1,K,2))+(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*  &
          CC(M,I-1,K,5)))+TR11*((WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*  &
          CC(M,I-1,K,3))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*  &
          CC(M,I-1,K,4))))+(TI12*((WA4(I-2)*CC(M,I-1,K,5)+  &
          WA4(I-1)*CC(M,I,K,5))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*  &
          CC(M,I,K,2)))-TI11*((WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*  &
          CC(M,I,K,4))-(WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)* CC(M,I,K,3))))
      CH(M,IC,4,K) = (TI12*((WA4(I-2)*CC(M,I-1,K,5)+WA4(I-1)*  &
          CC(M,I,K,5))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*  &
          CC(M,I,K,2)))-TI11*((WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*  &
          CC(M,I,K,4))-(WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*  &
          CC(M,I,K,3))))-(CC(M,I,K,1)+TR12*((WA1(I-2)*CC(M,I,K,2)-  &
          WA1(I-1)*CC(M,I-1,K,2))+(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*  &
          CC(M,I-1,K,5)))+TR11*((WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*  &
          CC(M,I-1,K,3))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)* CC(M,I-1,K,4))))
    END DO
  END DO
END DO
RETURN
END SUBROUTINE VRADF5


SUBROUTINE VRADFG(MP,IDO,IP,L1,IDL1,CC,C1,C2,MDIMC, CH,CH2,MDIMCH,WA)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989


INTEGER                       :: MP, IDO
INTEGER                       :: IP, L1
INTEGER                       :: IDL1
INTEGER                   :: MDIMC
INTEGER                   :: MDIMCH
REAL(EB)   CC(MDIMC,IDO,IP,L1)
REAL(EB)   C1(MDIMC,IDO,L1,IP)
REAL(EB)    C2(MDIMC,IDL1,IP)
REAL(EB)   CH2(MDIMCH,IDL1,IP)
REAL(EB)    WA(IDO)
REAL(EB)    CH(MDIMCH,IDO,L1,IP)
INTEGER :: IPPH, IPP2, IDP2, NBD, IK, M, J, K, IS, IDIJ, I, JC, L, LC, J2, IC
REAL(EB) :: TPI, ARG, DCP, DSP, AR1, AI1, AR1H, DC2, DS2, AR2, AI2, AR2H


TPI=2._EB*PI
ARG = TPI/REAL(IP,EB)
DCP = COS(ARG)
DSP = SIN(ARG)
IPPH = (IP+1)/2
IPP2 = IP+2
IDP2 = IDO+2
NBD = (IDO-1)/2
IF (IDO == 1) GO TO 119
DO  IK=1,IDL1
  DO  M=1,MP
    CH2(M,IK,1) = C2(M,IK,1)
  END DO
END DO
DO  J=2,IP
  DO  K=1,L1
    DO  M=1,MP
      CH(M,1,K,J) = C1(M,1,K,J)
    END DO
  END DO
END DO
IF (NBD > L1) GO TO 107
IS = -IDO
DO  J=2,IP
  IS = IS+IDO
  IDIJ = IS
  DO  I=3,IDO,2
    IDIJ = IDIJ+2
    DO  K=1,L1
      DO  M=1,MP
        CH(M,I-1,K,J) = WA(IDIJ-1)*C1(M,I-1,K,J)+WA(IDIJ) *C1(M,I,K,J)
        CH(M,I,K,J) = WA(IDIJ-1)*C1(M,I,K,J)-WA(IDIJ) *C1(M,I-1,K,J)
      END DO
    END DO
  END DO
END DO
GO TO 111
107 IS = -IDO
DO  J=2,IP
  IS = IS+IDO
  DO  K=1,L1
    IDIJ = IS
    DO  I=3,IDO,2
      IDIJ = IDIJ+2
      DO  M=1,MP
        CH(M,I-1,K,J) = WA(IDIJ-1)*C1(M,I-1,K,J)+WA(IDIJ) *C1(M,I,K,J)
        CH(M,I,K,J) = WA(IDIJ-1)*C1(M,I,K,J)-WA(IDIJ) *C1(M,I-1,K,J)
      END DO
    END DO
  END DO
END DO
111 IF (NBD < L1) GO TO 115
DO  J=2,IPPH
  JC = IPP2-J
  DO  K=1,L1
    DO  I=3,IDO,2
      DO  M=1,MP
        C1(M,I-1,K,J) = CH(M,I-1,K,J)+CH(M,I-1,K,JC)
        C1(M,I-1,K,JC) = CH(M,I,K,J)-CH(M,I,K,JC)
        C1(M,I,K,J) = CH(M,I,K,J)+CH(M,I,K,JC)
        C1(M,I,K,JC) = CH(M,I-1,K,JC)-CH(M,I-1,K,J)
      END DO
    END DO
  END DO
END DO
GO TO 121
115 CONTINUE
DO  J=2,IPPH
  JC = IPP2-J
  DO  I=3,IDO,2
    DO  K=1,L1
      DO  M=1,MP
        C1(M,I-1,K,J) = CH(M,I-1,K,J)+CH(M,I-1,K,JC)
        C1(M,I-1,K,JC) = CH(M,I,K,J)-CH(M,I,K,JC)
        C1(M,I,K,J) = CH(M,I,K,J)+CH(M,I,K,JC)
        C1(M,I,K,JC) = CH(M,I-1,K,JC)-CH(M,I-1,K,J)
      END DO
    END DO
  END DO
END DO
GO TO 121
119 CONTINUE
DO  IK=1,IDL1
  DO  M=1,MP
    C2(M,IK,1) = CH2(M,IK,1)
  END DO
END DO
121 CONTINUE
DO  J=2,IPPH
  JC = IPP2-J
  DO  K=1,L1
    DO  M=1,MP
      C1(M,1,K,J) = CH(M,1,K,J)+CH(M,1,K,JC)
      C1(M,1,K,JC) = CH(M,1,K,JC)-CH(M,1,K,J)
    END DO
  END DO
END DO

AR1 = 1._EB
AI1 = 0._EB
DO  L=2,IPPH
  LC = IPP2-L
  AR1H = DCP*AR1-DSP*AI1
  AI1 = DCP*AI1+DSP*AR1
  AR1 = AR1H
  DO  IK=1,IDL1
    DO  M=1,MP
      CH2(M,IK,L) = C2(M,IK,1)+AR1*C2(M,IK,2)
      CH2(M,IK,LC) = AI1*C2(M,IK,IP)
    END DO
  END DO
  DC2 = AR1
  DS2 = AI1
  AR2 = AR1
  AI2 = AI1
  DO  J=3,IPPH
    JC = IPP2-J
    AR2H = DC2*AR2-DS2*AI2
    AI2 = DC2*AI2+DS2*AR2
    AR2 = AR2H
    DO  IK=1,IDL1
      DO  M=1,MP
        CH2(M,IK,L) = CH2(M,IK,L)+AR2*C2(M,IK,J)
        CH2(M,IK,LC) = CH2(M,IK,LC)+AI2*C2(M,IK,JC)
      END DO
    END DO
  END DO
END DO
DO  M=1,MP
  DO  IK=1,IDL1
    DO  J=2,IPPH
      CH2(M,IK,1) = CH2(M,IK,1)+C2(M,IK,J)
    END DO
  END DO
END DO

IF (IDO < L1) GO TO 132
DO  K=1,L1
  DO  I=1,IDO
    DO  M=1,MP
      CC(M,I,1,K) = CH(M,I,K,1)
    END DO
  END DO
END DO
GO TO 135
132 CONTINUE
DO  I=1,IDO
  DO  K=1,L1
    DO  M=1,MP
      CC(M,I,1,K) = CH(M,I,K,1)
    END DO
  END DO
END DO
135 CONTINUE
DO  J=2,IPPH
  JC = IPP2-J
  J2 = J+J
  DO  K=1,L1
    DO  M=1,MP
      CC(M,IDO,J2-2,K) = CH(M,1,K,J)
      CC(M,1,J2-1,K) = CH(M,1,K,JC)
    END DO
  END DO
END DO
IF (IDO == 1) RETURN
IF (NBD < L1) GO TO 141
DO  J=2,IPPH
  JC = IPP2-J
  J2 = J+J
  DO  K=1,L1
    DO  I=3,IDO,2
      IC = IDP2-I
      DO  M=1,MP
        CC(M,I-1,J2-1,K) = CH(M,I-1,K,J)+CH(M,I-1,K,JC)
        CC(M,IC-1,J2-2,K) = CH(M,I-1,K,J)-CH(M,I-1,K,JC)
        CC(M,I,J2-1,K) = CH(M,I,K,J)+CH(M,I,K,JC)
        CC(M,IC,J2-2,K) = CH(M,I,K,JC)-CH(M,I,K,J)
      END DO
    END DO
  END DO
END DO
RETURN
141 CONTINUE
DO  J=2,IPPH
  JC = IPP2-J
  J2 = J+J
  DO  I=3,IDO,2
    IC = IDP2-I
    DO  K=1,L1
      DO  M=1,MP
        CC(M,I-1,J2-1,K) = CH(M,I-1,K,J)+CH(M,I-1,K,JC)
        CC(M,IC-1,J2-2,K) = CH(M,I-1,K,J)-CH(M,I-1,K,JC)
        CC(M,I,J2-1,K) = CH(M,I,K,J)+CH(M,I,K,JC)
        CC(M,IC,J2-2,K) = CH(M,I,K,JC)-CH(M,I,K,J)
      END DO
    END DO
  END DO
END DO
RETURN
END SUBROUTINE VRADFG


SUBROUTINE VRFFTB(M,N,R,MDIMR,RT,WSAVE)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989



INTEGER                   :: M
INTEGER                   :: N
INTEGER                   :: MDIMR
REAL(EB)       RT(M,N)
REAL(EB)       WSAVE(N+15)
REAL(EB)       R(MDIMR,N)
IF (N == 1) RETURN
CALL VRFTB1 (M,N,R,MDIMR,RT,WSAVE(1),WSAVE(N+1))
RETURN
END SUBROUTINE VRFFTB


SUBROUTINE VRFTB1 (M,N,C,MDIMC,CH,WA,FAC)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989



INTEGER                       :: M
INTEGER                       :: N
INTEGER   :: MDIMC
REAL(EB)  C(MDIMC,N)
REAL(EB)    WA(N)
REAL(EB)    FAC(15)
REAL(EB)    CH(M,N)
INTEGER :: I, NF, NA, L1, IW, K1, IP, L2, IDO, IDL1, IX2, IX3, IX4, J

NF = INT(FAC(2))
NA = 0
L1 = 1
IW = 1
DO  K1=1,NF
  IP = INT(FAC(K1+2))
  L2 = IP*L1
  IDO = N/L2
  IDL1 = IDO*L1
  IF (IP /= 4) GO TO 103
  IX2 = IW+IDO
  IX3 = IX2+IDO
  IF (NA /= 0) GO TO 101
  CALL VRADB4 (M,IDO,L1,C,MDIMC,CH,M,WA(IW),WA(IX2),WA(IX3))
  GO TO 102
  101    CALL VRADB4 (M,IDO,L1,CH,M,C,MDIMC,WA(IW),WA(IX2),WA(IX3))
  102    NA = 1-NA
  GO TO 115
  103    IF (IP /= 2) GO TO 106
  IF (NA /= 0) GO TO 104
  CALL VRADB2 (M,IDO,L1,C,MDIMC,CH,M,WA(IW))
  GO TO 105
  104    CALL VRADB2 (M,IDO,L1,CH,M,C,MDIMC,WA(IW))
  105    NA = 1-NA
  GO TO 115
  106    IF (IP /= 3) GO TO 109
  IX2 = IW+IDO
  IF (NA /= 0) GO TO 107
  CALL VRADB3 (M,IDO,L1,C,MDIMC,CH,M,WA(IW),WA(IX2))
  GO TO 108
  107    CALL VRADB3 (M,IDO,L1,CH,M,C,MDIMC,WA(IW),WA(IX2))
  108    NA = 1-NA
  GO TO 115
  109    IF (IP /= 5) GO TO 112
  IX2 = IW+IDO
  IX3 = IX2+IDO
  IX4 = IX3+IDO
  IF (NA /= 0) GO TO 110
  CALL VRADB5 (M,IDO,L1,C,MDIMC,CH,M,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  GO TO 111
  110 CALL VRADB5 (M,IDO,L1,CH,M,C,MDIMC,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  111    NA = 1-NA
  GO TO 115
  112    IF (NA /= 0) GO TO 113
  CALL VRADBG (M,IDO,IP,L1,IDL1,C,C,C,MDIMC,CH,CH,M,WA(IW))
  GO TO 114
  113    CALL VRADBG (M,IDO,IP,L1,IDL1,CH,CH,CH,M,C,C,MDIMC,WA(IW))
  114    IF (IDO == 1) NA = 1-NA
  115    L1 = L2
  IW = IW+(IP-1)*IDO
END DO
OUTARY=.TRUE.
IF (NOCOPY) THEN
  SCALE=SCALE*SQRT(1.0_EB/REAL(N,EB))
  IF (NA==1) THEN
    OUTARY=.FALSE.
  END IF
ELSE
  SCALE=SQRT(1.0_EB/REAL(N,EB))
  IF (NA == 0) GO TO 118
  DO  J=1,N
    DO  I=1,M
      C(I,J) = SCALE*CH(I,J)
    END DO
  END DO
  RETURN
  118    DO  J=1,N
    DO  I=1,M
      C(I,J)=SCALE*C(I,J)
    END DO
  END DO
END IF
RETURN
END SUBROUTINE VRFTB1


SUBROUTINE VRADB2 (MP,IDO,L1,CC,MDIMC,CH,MDIMCH,WA1)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989


INTEGER                       :: MP
INTEGER                   :: IDO
INTEGER                       :: L1
INTEGER                   :: MDIMC
INTEGER                   :: MDIMCH
REAL(EB)   CH(MDIMCH,IDO,L1,2)
REAL(EB) WA1(IDO)
REAL(EB)    CC(MDIMC,IDO,2,L1)
INTEGER :: I, K, M, IDP2, IC

DO  K=1,L1
  DO  M=1,MP
    CH(M,1,K,1) = CC(M,1,1,K)+CC(M,IDO,2,K)
    CH(M,1,K,2) = CC(M,1,1,K)-CC(M,IDO,2,K)
  END DO
END DO
IF (IDO-2<0) GO TO 107
IF (IDO-2==0) GO TO 105
IF (IDO-2>0) GO TO 102
102 IDP2 = IDO+2
DO  K=1,L1
  DO  I=3,IDO,2
    IC = IDP2-I
    DO  M=1,MP
      CH(M,I-1,K,1) = CC(M,I-1,1,K)+CC(M,IC-1,2,K)
      CH(M,I,K,1) = CC(M,I,1,K)-CC(M,IC,2,K)
      CH(M,I-1,K,2) = WA1(I-2)*(CC(M,I-1,1,K)-CC(M,IC-1,2,K))  &
          -WA1(I-1)*(CC(M,I,1,K)+CC(M,IC,2,K))
      CH(M,I,K,2) = WA1(I-2)*(CC(M,I,1,K)+CC(M,IC,2,K))+WA1(I-1)  &
          *(CC(M,I-1,1,K)-CC(M,IC-1,2,K))
    END DO
  END DO
END DO
IF (MOD(IDO,2) == 1) RETURN
105 DO  K=1,L1
  DO  M=1,MP
    CH(M,IDO,K,1) = CC(M,IDO,1,K)+CC(M,IDO,1,K)
    CH(M,IDO,K,2) = -(CC(M,1,2,K)+CC(M,1,2,K))
  END DO
END DO
107 RETURN
END SUBROUTINE VRADB2


SUBROUTINE VRADB3 (MP,IDO,L1,CC,MDIMC,CH,MDIMCH,WA1,WA2)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989


INTEGER                       :: MP
INTEGER                       :: IDO
INTEGER                       :: L1
INTEGER                   :: MDIMC
INTEGER                   :: MDIMCH
REAL(EB)   CH(MDIMCH,IDO,L1,3)
REAL(EB)   WA1(IDO)
REAL(EB)    WA2(IDO)
REAL(EB)    CC(MDIMC,IDO,3,L1)
INTEGER :: I, M, K, IDP2, IC
REAL(EB) :: ARG, TAUR, TAUI

ARG=2._EB*PI/3._EB
TAUR=COS(ARG)
TAUI=SIN(ARG)
DO  M=1,MP
   DO  K=1,L1
      CH(M,1,K,1) = CC(M,1,1,K)+2._EB*CC(M,IDO,2,K)
      CH(M,1,K,2) = CC(M,1,1,K)+(2._EB*TAUR)*CC(M,IDO,2,K) -(2._EB*TAUI)*CC(M,1,3,K)
      CH(M,1,K,3) = CC(M,1,1,K)+(2._EB*TAUR)*CC(M,IDO,2,K) +2._EB*TAUI*CC(M,1,3,K)
   END DO
END DO
IF (IDO == 1) RETURN
IDP2 = IDO+2
DO  M=1,MP
   DO  I=3,IDO,2
      IC = IDP2-I
      DO  K=1,L1
         CH(M,I-1,K,1) = CC(M,I-1,1,K)+(CC(M,I-1,3,K)+CC(M,IC-1,2,K))
         CH(M,I,K,1) = CC(M,I,1,K)+(CC(M,I,3,K)-CC(M,IC,2,K))
         CH(M,I-1,K,2) = WA1(I-2)*  &
               ((CC(M,I-1,1,K)+TAUR*(CC(M,I-1,3,K)+CC(M,IC-1,2,K)))-  &
               (TAUI*(CC(M,I,3,K)+CC(M,IC,2,K)))) -WA1(I-1)*  &
               ((CC(M,I,1,K)+TAUR*(CC(M,I,3,K)-CC(M,IC,2,K)))+  &
               (TAUI*(CC(M,I-1,3,K)-CC(M,IC-1,2,K))))
         CH(M,I,K,2) = WA1(I-2)*  &
               ((CC(M,I,1,K)+TAUR*(CC(M,I,3,K)-CC(M,IC,2,K)))+  &
               (TAUI*(CC(M,I-1,3,K)-CC(M,IC-1,2,K)))) +WA1(I-1)*  &
               ((CC(M,I-1,1,K)+TAUR*(CC(M,I-1,3,K)+CC(M,IC-1,2,K)))-  &
               (TAUI*(CC(M,I,3,K)+CC(M,IC,2,K))))
         CH(M,I-1,K,3) = WA2(I-2)*  &
               ((CC(M,I-1,1,K)+TAUR*(CC(M,I-1,3,K)+CC(M,IC-1,2,K)))+  &
               (TAUI*(CC(M,I,3,K)+CC(M,IC,2,K)))) -WA2(I-1)*  &
               ((CC(M,I,1,K)+TAUR*(CC(M,I,3,K)-CC(M,IC,2,K)))-  &
               (TAUI*(CC(M,I-1,3,K)-CC(M,IC-1,2,K))))
         CH(M,I,K,3) = WA2(I-2)*  &
               ((CC(M,I,1,K)+TAUR*(CC(M,I,3,K)-CC(M,IC,2,K)))-  &
               (TAUI*(CC(M,I-1,3,K)-CC(M,IC-1,2,K)))) +WA2(I-1)*  &
               ((CC(M,I-1,1,K)+TAUR*(CC(M,I-1,3,K)+CC(M,IC-1,2,K)))+  &
               (TAUI*(CC(M,I,3,K)+CC(M,IC,2,K))))
      END DO
   END DO
END DO
RETURN
END SUBROUTINE VRADB3


SUBROUTINE VRADB4 (MP,IDO,L1,CC,MDIMC,CH,MDIMCH,WA1,WA2,WA3)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989


INTEGER                       :: MP
INTEGER                   :: IDO
INTEGER                       :: L1
INTEGER                   :: MDIMC
INTEGER                   :: MDIMCH
REAL(EB)     CH(MDIMCH,IDO,L1,4)
REAL(EB)   WA1(IDO)
REAL(EB)    WA2(IDO)
REAL(EB)   WA3(IDO)
REAL(EB)    CC(MDIMC,IDO,4,L1)
INTEGER :: I, M, K, IDP2, IC
REAL(EB) :: SQRT2

SQRT2=SQRT(2.0_EB)
DO  M=1,MP
   DO  K=1,L1
      CH(M,1,K,3) = (CC(M,1,1,K)+CC(M,IDO,4,K)) -(CC(M,IDO,2,K)+CC(M,IDO,2,K))
      CH(M,1,K,1) = (CC(M,1,1,K)+CC(M,IDO,4,K)) +(CC(M,IDO,2,K)+CC(M,IDO,2,K))
      CH(M,1,K,4) = (CC(M,1,1,K)-CC(M,IDO,4,K)) +(CC(M,1,3,K)+CC(M,1,3,K))
      CH(M,1,K,2) = (CC(M,1,1,K)-CC(M,IDO,4,K)) -(CC(M,1,3,K)+CC(M,1,3,K))
   END DO
END DO
IF (IDO<2) RETURN
IF (IDO>2) THEN
   IDP2 = IDO+2
   DO  M=1,MP
      DO  I=3,IDO,2
         IC = IDP2-I
         DO  K=1,L1
            CH(M,I-1,K,1) = (CC(M,I-1,1,K)+CC(M,IC-1,4,K))  &
                +(CC(M,I-1,3,K)+CC(M,IC-1,2,K))
            CH(M,I,K,1) = (CC(M,I,1,K)-CC(M,IC,4,K)) +(CC(M,I,3,K)-CC(M,IC,2,K))
            CH(M,I-1,K,2)=WA1(I-2)*((CC(M,I-1,1,K)-CC(M,IC-1,4,K))  &
                -(CC(M,I,3,K)+CC(M,IC,2,K)))-WA1(I-1)  &
                *((CC(M,I,1,K)+CC(M,IC,4,K))+(CC(M,I-1,3,K)-CC(M,IC-1,2,K)))
            CH(M,I,K,2)=WA1(I-2)*((CC(M,I,1,K)+CC(M,IC,4,K))  &
                +(CC(M,I-1,3,K)-CC(M,IC-1,2,K)))+WA1(I-1)  &
                *((CC(M,I-1,1,K)-CC(M,IC-1,4,K))-(CC(M,I,3,K)+CC(M,IC,2,K)))
            CH(M,I-1,K,3)=WA2(I-2)*((CC(M,I-1,1,K)+CC(M,IC-1,4,K))  &
                -(CC(M,I-1,3,K)+CC(M,IC-1,2,K)))-WA2(I-1)  &
                *((CC(M,I,1,K)-CC(M,IC,4,K))-(CC(M,I,3,K)-CC(M,IC,2,K)))
            CH(M,I,K,3)=WA2(I-2)*((CC(M,I,1,K)-CC(M,IC,4,K))  &
                -(CC(M,I,3,K)-CC(M,IC,2,K)))+WA2(I-1)  &
                *((CC(M,I-1,1,K)+CC(M,IC-1,4,K))-(CC(M,I-1,3,K) +CC(M,IC-1,2,K)))
            CH(M,I-1,K,4)=WA3(I-2)*((CC(M,I-1,1,K)-CC(M,IC-1,4,K))  &
                +(CC(M,I,3,K)+CC(M,IC,2,K)))-WA3(I-1)  &
                *((CC(M,I,1,K)+CC(M,IC,4,K))-(CC(M,I-1,3,K)-CC(M,IC-1,2,K)))
            CH(M,I,K,4)=WA3(I-2)*((CC(M,I,1,K)+CC(M,IC,4,K))  &
                -(CC(M,I-1,3,K)-CC(M,IC-1,2,K)))+WA3(I-1)  &
                *((CC(M,I-1,1,K)-CC(M,IC-1,4,K))+(CC(M,I,3,K)+CC(M,IC,2,K)))
         END DO
      END DO
   END DO
   IF (MOD(IDO,2) == 1) RETURN
ENDIF

DO  M=1,MP
   DO  K=1,L1
      CH(M,IDO,K,1) = (CC(M,IDO,1,K)+CC(M,IDO,3,K))  &
         +(CC(M,IDO,1,K)+CC(M,IDO,3,K))
      CH(M,IDO,K,2) = SQRT2*((CC(M,IDO,1,K)-CC(M,IDO,3,K))  &
         -(CC(M,1,2,K)+CC(M,1,4,K)))
      CH(M,IDO,K,3) = (CC(M,1,4,K)-CC(M,1,2,K)) +(CC(M,1,4,K)-CC(M,1,2,K))
      CH(M,IDO,K,4) = -SQRT2*((CC(M,IDO,1,K)-CC(M,IDO,3,K))  &
         +(CC(M,1,2,K)+CC(M,1,4,K)))
   END DO
END DO

END SUBROUTINE VRADB4


SUBROUTINE VRADB5 (MP,IDO,L1,CC,MDIMC,CH,MDIMCH,WA1,WA2,WA3,WA4)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989


INTEGER                       :: MP
INTEGER                       :: IDO
INTEGER                       :: L1
INTEGER                   :: MDIMC
INTEGER                   :: MDIMCH
REAL(EB)   CH(MDIMCH,IDO,L1,5)
REAL(EB)   WA1(IDO)
REAL(EB)  WA2(IDO)
REAL(EB)   WA3(IDO)
REAL(EB)   WA4(IDO)
REAL(EB)    CC(MDIMC,IDO,5,L1)
INTEGER :: I, K, M, IDP2, IC
REAL(EB) :: ARG, TR11, TI11, TR12, TI12

ARG=2._EB*PI/5._EB
TR11=COS(ARG)
TI11=SIN(ARG)
TR12=COS(2._EB*ARG)
TI12=SIN(2._EB*ARG)
DO  K=1,L1
  DO  M=1,MP
    CH(M,1,K,1) = CC(M,1,1,K)+2._EB*CC(M,IDO,2,K)+2._EB*CC(M,IDO,4,K)
    CH(M,1,K,2) = (CC(M,1,1,K)+TR11*2._EB*CC(M,IDO,2,K)  &
        +TR12*2._EB*CC(M,IDO,4,K))-(TI11*2._EB*CC(M,1,3,K) +TI12*2._EB*CC(M,1,5,K))
    CH(M,1,K,3) = (CC(M,1,1,K)+TR12*2._EB*CC(M,IDO,2,K)  &
        +TR11*2._EB*CC(M,IDO,4,K))-(TI12*2._EB*CC(M,1,3,K) -TI11*2._EB*CC(M,1,5,K))
    CH(M,1,K,4) = (CC(M,1,1,K)+TR12*2._EB*CC(M,IDO,2,K)  &
        +TR11*2._EB*CC(M,IDO,4,K))+(TI12*2._EB*CC(M,1,3,K) -TI11*2._EB*CC(M,1,5,K))
    CH(M,1,K,5) = (CC(M,1,1,K)+TR11*2._EB*CC(M,IDO,2,K)  &
        +TR12*2._EB*CC(M,IDO,4,K))+(TI11*2._EB*CC(M,1,3,K) +TI12*2._EB*CC(M,1,5,K))
  END DO
END DO
IF (IDO == 1) RETURN
IDP2 = IDO+2
DO  K=1,L1
  DO  I=3,IDO,2
    IC = IDP2-I
    DO  M=1,MP
      CH(M,I-1,K,1) = CC(M,I-1,1,K)+(CC(M,I-1,3,K)+CC(M,IC-1,2,K))  &
          +(CC(M,I-1,5,K)+CC(M,IC-1,4,K))
      CH(M,I,K,1) = CC(M,I,1,K)+(CC(M,I,3,K)-CC(M,IC,2,K))  &
          +(CC(M,I,5,K)-CC(M,IC,4,K))
      CH(M,I-1,K,2) = WA1(I-2)*((CC(M,I-1,1,K)+TR11*  &
          (CC(M,I-1,3,K)+CC(M,IC-1,2,K))+TR12  &
          *(CC(M,I-1,5,K)+CC(M,IC-1,4,K)))-(TI11*(CC(M,I,3,K)  &
          +CC(M,IC,2,K))+TI12*(CC(M,I,5,K)+CC(M,IC,4,K))))  &
          -WA1(I-1)*((CC(M,I,1,K)+TR11*(CC(M,I,3,K)-CC(M,IC,2,K))  &
          +TR12*(CC(M,I,5,K)-CC(M,IC,4,K)))+(TI11*(CC(M,I-1,3,K)  &
          -CC(M,IC-1,2,K))+TI12*(CC(M,I-1,5,K)-CC(M,IC-1,4,K))))
      CH(M,I,K,2) = WA1(I-2)*((CC(M,I,1,K)+TR11*(CC(M,I,3,K)  &
          -CC(M,IC,2,K))+TR12*(CC(M,I,5,K)-CC(M,IC,4,K)))  &
          +(TI11*(CC(M,I-1,3,K)-CC(M,IC-1,2,K))+TI12  &
          *(CC(M,I-1,5,K)-CC(M,IC-1,4,K))))+WA1(I-1)  &
          *((CC(M,I-1,1,K)+TR11*(CC(M,I-1,3,K)  &
          +CC(M,IC-1,2,K))+TR12*(CC(M,I-1,5,K)+CC(M,IC-1,4,K)))  &
          -(TI11*(CC(M,I,3,K)+CC(M,IC,2,K))+TI12 *(CC(M,I,5,K)+CC(M,IC,4,K))))
      CH(M,I-1,K,3) = WA2(I-2)  &
          *((CC(M,I-1,1,K)+TR12*(CC(M,I-1,3,K)+CC(M,IC-1,2,K))  &
          +TR11*(CC(M,I-1,5,K)+CC(M,IC-1,4,K)))-(TI12*(CC(M,I,3,K)  &
          +CC(M,IC,2,K))-TI11*(CC(M,I,5,K)+CC(M,IC,4,K)))) -WA2(I-1)  &
          *((CC(M,I,1,K)+TR12*(CC(M,I,3,K)-  &
          CC(M,IC,2,K))+TR11*(CC(M,I,5,K)-CC(M,IC,4,K)))  &
          +(TI12*(CC(M,I-1,3,K)-CC(M,IC-1,2,K))-TI11  &
          *(CC(M,I-1,5,K)-CC(M,IC-1,4,K))))
      CH(M,I,K,3) = WA2(I-2) *((CC(M,I,1,K)+TR12*(CC(M,I,3,K)-  &
          CC(M,IC,2,K))+TR11*(CC(M,I,5,K)-CC(M,IC,4,K)))  &
          +(TI12*(CC(M,I-1,3,K)-CC(M,IC-1,2,K))-TI11  &
          *(CC(M,I-1,5,K)-CC(M,IC-1,4,K)))) +WA2(I-1)  &
          *((CC(M,I-1,1,K)+TR12*(CC(M,I-1,3,K)+CC(M,IC-1,2,K))  &
          +TR11*(CC(M,I-1,5,K)+CC(M,IC-1,4,K)))-(TI12*(CC(M,I,3,K)  &
          +CC(M,IC,2,K))-TI11*(CC(M,I,5,K)+CC(M,IC,4,K))))
      CH(M,I-1,K,4) = WA3(I-2)  &
          *((CC(M,I-1,1,K)+TR12*(CC(M,I-1,3,K)+CC(M,IC-1,2,K))  &
          +TR11*(CC(M,I-1,5,K)+CC(M,IC-1,4,K)))+(TI12*(CC(M,I,3,K)  &
          +CC(M,IC,2,K))-TI11*(CC(M,I,5,K)+CC(M,IC,4,K)))) -WA3(I-1)  &
          *((CC(M,I,1,K)+TR12*(CC(M,I,3,K)-  &
          CC(M,IC,2,K))+TR11*(CC(M,I,5,K)-CC(M,IC,4,K)))  &
          -(TI12*(CC(M,I-1,3,K)-CC(M,IC-1,2,K))-TI11  &
          *(CC(M,I-1,5,K)-CC(M,IC-1,4,K))))
      CH(M,I,K,4) = WA3(I-2) *((CC(M,I,1,K)+TR12*(CC(M,I,3,K)-  &
          CC(M,IC,2,K))+TR11*(CC(M,I,5,K)-CC(M,IC,4,K)))  &
          -(TI12*(CC(M,I-1,3,K)-CC(M,IC-1,2,K))-TI11  &
          *(CC(M,I-1,5,K)-CC(M,IC-1,4,K)))) +WA3(I-1)  &
          *((CC(M,I-1,1,K)+TR12*(CC(M,I-1,3,K)+CC(M,IC-1,2,K))  &
          +TR11*(CC(M,I-1,5,K)+CC(M,IC-1,4,K)))+(TI12*(CC(M,I,3,K)  &
          +CC(M,IC,2,K))-TI11*(CC(M,I,5,K)+CC(M,IC,4,K))))
      CH(M,I-1,K,5) = WA4(I-2)  &
          *((CC(M,I-1,1,K)+TR11*(CC(M,I-1,3,K)+CC(M,IC-1,2,K))  &
          +TR12*(CC(M,I-1,5,K)+CC(M,IC-1,4,K)))+(TI11*(CC(M,I,3,K)  &
          +CC(M,IC,2,K))+TI12*(CC(M,I,5,K)+CC(M,IC,4,K)))) -WA4(I-1)  &
          *((CC(M,I,1,K)+TR11*(CC(M,I,3,K)-CC(M,IC,2,K))  &
          +TR12*(CC(M,I,5,K)-CC(M,IC,4,K)))-(TI11*(CC(M,I-1,3,K)  &
          -CC(M,IC-1,2,K))+TI12*(CC(M,I-1,5,K)-CC(M,IC-1,4,K))))
      CH(M,I,K,5) = WA4(I-2) *((CC(M,I,1,K)+TR11*(CC(M,I,3,K)-CC(M,IC,2,K))  &
          +TR12*(CC(M,I,5,K)-CC(M,IC,4,K)))-(TI11*(CC(M,I-1,3,K)  &
          -CC(M,IC-1,2,K))+TI12*(CC(M,I-1,5,K)-CC(M,IC-1,4,K)))) +WA4(I-1)  &
          *((CC(M,I-1,1,K)+TR11*(CC(M,I-1,3,K)+CC(M,IC-1,2,K))  &
          +TR12*(CC(M,I-1,5,K)+CC(M,IC-1,4,K)))+(TI11*(CC(M,I,3,K)  &
          +CC(M,IC,2,K))+TI12*(CC(M,I,5,K)+CC(M,IC,4,K))))
    END DO
  END DO
END DO
RETURN
END SUBROUTINE VRADB5


SUBROUTINE VRADBG (MP,IDO,IP,L1,IDL1,CC,C1,C2,MDIMC, CH,CH2,MDIMCH,WA)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!     PACKAGE VFFTPAK, VERSION 1, JUNE 1989
INTEGER                       :: MP
INTEGER      :: IDO, L1
INTEGER                       :: IP
INTEGER                       :: IDL1
INTEGER                   :: MDIMC
INTEGER                   :: MDIMCH
REAL(EB)      CH(MDIMCH,IDO,L1,IP)    ,CC(MDIMC,IDO,IP,L1) ,C1(MDIMC,IDO,L1,IP)     ,C2(MDIMC,IDL1,IP), &
               CH2(MDIMCH,IDL1,IP)       ,WA(IDO)
INTEGER :: I, IDP2, NBD, IPP2, IPPH, K, M, J, JC, J2, IC, L, LC, IK, IS, IDIJ
REAL(EB) :: TPI, ARG, DCP, DSP, AR1, AI1, AR1H, DC2, DS2, AR2, AI2, AR2H


TPI=2._EB*PI
ARG = TPI/REAL(IP,EB)
DCP = COS(ARG)
DSP = SIN(ARG)
IDP2 = IDO+2
NBD = (IDO-1)/2
IPP2 = IP+2
IPPH = (IP+1)/2
IF (IDO < L1) GO TO 103
DO  K=1,L1
  DO  I=1,IDO
    DO  M=1,MP
      CH(M,I,K,1) = CC(M,I,1,K)
    END DO
  END DO
END DO
GO TO 106
103 CONTINUE
DO  I=1,IDO
  DO  K=1,L1
    DO  M=1,MP
      CH(M,I,K,1) = CC(M,I,1,K)
    END DO
  END DO
END DO
106 CONTINUE
DO  J=2,IPPH
  JC = IPP2-J
  J2 = J+J
  DO  K=1,L1
    DO  M=1,MP
      CH(M,1,K,J) = CC(M,IDO,J2-2,K)+CC(M,IDO,J2-2,K)
      CH(M,1,K,JC) = CC(M,1,J2-1,K)+CC(M,1,J2-1,K)
    END DO
  END DO
END DO
IF (IDO == 1) GO TO 116
IF (NBD < L1) GO TO 112
DO  J=2,IPPH
  JC = IPP2-J
  DO  K=1,L1
    DO  I=3,IDO,2
      IC = IDP2-I
      DO  M=1,MP
        CH(M,I-1,K,J) = CC(M,I-1,2*J-1,K)+CC(M,IC-1,2*J-2,K)
        CH(M,I-1,K,JC) = CC(M,I-1,2*J-1,K)-CC(M,IC-1,2*J-2,K)
        CH(M,I,K,J) = CC(M,I,2*J-1,K)-CC(M,IC,2*J-2,K)
        CH(M,I,K,JC) = CC(M,I,2*J-1,K)+CC(M,IC,2*J-2,K)
      END DO
    END DO
  END DO
END DO
GO TO 116
112 CONTINUE
DO  J=2,IPPH
  JC = IPP2-J
  DO  I=3,IDO,2
    IC = IDP2-I
    DO  K=1,L1
      DO  M=1,MP
        CH(M,I-1,K,J) = CC(M,I-1,2*J-1,K)+CC(M,IC-1,2*J-2,K)
        CH(M,I-1,K,JC) = CC(M,I-1,2*J-1,K)-CC(M,IC-1,2*J-2,K)
        CH(M,I,K,J) = CC(M,I,2*J-1,K)-CC(M,IC,2*J-2,K)
        CH(M,I,K,JC) = CC(M,I,2*J-1,K)+CC(M,IC,2*J-2,K)
      END DO
    END DO
  END DO
END DO
116 AR1 = 1._EB
AI1 = 0._EB
DO  L=2,IPPH
  LC = IPP2-L
  AR1H = DCP*AR1-DSP*AI1
  AI1 = DCP*AI1+DSP*AR1
  AR1 = AR1H
  DO  IK=1,IDL1
    DO  M=1,MP
      C2(M,IK,L) = CH2(M,IK,1)+AR1*CH2(M,IK,2)
      C2(M,IK,LC) = AI1*CH2(M,IK,IP)
    END DO
  END DO
  DC2 = AR1
  DS2 = AI1
  AR2 = AR1
  AI2 = AI1
  DO  J=3,IPPH
    JC = IPP2-J
    AR2H = DC2*AR2-DS2*AI2
    AI2 = DC2*AI2+DS2*AR2
    AR2 = AR2H
    DO  IK=1,IDL1
      DO  M=1,MP
        C2(M,IK,L) = C2(M,IK,L)+AR2*CH2(M,IK,J)
        C2(M,IK,LC) = C2(M,IK,LC)+AI2*CH2(M,IK,JC)
      END DO
    END DO
  END DO
END DO
DO  M=1,MP
  DO  IK=1,IDL1
    DO  J=2,IPPH
      CH2(M,IK,1) = CH2(M,IK,1)+CH2(M,IK,J)
    END DO
  END DO
END DO
DO  J=2,IPPH
  JC = IPP2-J
  DO  K=1,L1
    DO  M=1,MP
      CH(M,1,K,J) = C1(M,1,K,J)-C1(M,1,K,JC)
      CH(M,1,K,JC) = C1(M,1,K,J)+C1(M,1,K,JC)
    END DO
  END DO
END DO
IF (IDO == 1) GO TO 132
IF (NBD < L1) GO TO 128
DO  J=2,IPPH
  JC = IPP2-J
  DO  K=1,L1
    DO  I=3,IDO,2
      DO  M=1,MP
        CH(M,I-1,K,J) = C1(M,I-1,K,J)-C1(M,I,K,JC)
        CH(M,I-1,K,JC) = C1(M,I-1,K,J)+C1(M,I,K,JC)
        CH(M,I,K,J) = C1(M,I,K,J)+C1(M,I-1,K,JC)
        CH(M,I,K,JC) = C1(M,I,K,J)-C1(M,I-1,K,JC)
      END DO
    END DO
  END DO
END DO
GO TO 132
128 CONTINUE
DO  J=2,IPPH
  JC = IPP2-J
  DO  I=3,IDO,2
    DO  K=1,L1
      DO  M=1,MP
        CH(M,I-1,K,J) = C1(M,I-1,K,J)-C1(M,I,K,JC)
        CH(M,I-1,K,JC) = C1(M,I-1,K,J)+C1(M,I,K,JC)
        CH(M,I,K,J) = C1(M,I,K,J)+C1(M,I-1,K,JC)
        CH(M,I,K,JC) = C1(M,I,K,J)-C1(M,I-1,K,JC)
      END DO
    END DO
  END DO
END DO
132 CONTINUE
IF (IDO == 1) RETURN
DO  IK=1,IDL1
  DO  M=1,MP
    C2(M,IK,1) = CH2(M,IK,1)
  END DO
END DO
DO  J=2,IP
  DO  K=1,L1
    DO  M=1,MP
      C1(M,1,K,J) = CH(M,1,K,J)
    END DO
  END DO
END DO
IF (NBD > L1) GO TO 139
DO  J=2,IP
  IS = (J-2)*IDO
  IDIJ = IS
  DO  I=3,IDO,2
    IDIJ = IDIJ+2
    DO  K=1,L1
      DO  M=1,MP
        C1(M,I-1,K,J) = WA(IDIJ-1)*CH(M,I-1,K,J)-WA(IDIJ)* CH(M,I,K,J)
        C1(M,I,K,J) = WA(IDIJ-1)*CH(M,I,K,J)+WA(IDIJ)* CH(M,I-1,K,J)
      END DO
    END DO
  END DO
END DO
GO TO 143
139 CONTINUE
DO  J=2,IP
  IS = (J-2)*IDO
  DO  K=1,L1
    IDIJ = IS
    DO  I=3,IDO,2
      IDIJ = IDIJ+2
      DO  M=1,MP
        C1(M,I-1,K,J) = WA(IDIJ-1)*CH(M,I-1,K,J)-WA(IDIJ)* CH(M,I,K,J)
        C1(M,I,K,J) = WA(IDIJ-1)*CH(M,I,K,J)+WA(IDIJ)* CH(M,I-1,K,J)
      END DO
    END DO
  END DO
END DO
143 RETURN
END SUBROUTINE VRADBG


SUBROUTINE SSWAP(N,SX,INCX,SY,INCY)
!***BEGIN PROLOGUE  SSWAP
!***DATE WRITTEN   791001   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  D1A5
!***KEYWORDS  BLAS,INTERCHANGE,LINEAR ALGEBRA,VECTOR
!***AUTHOR  LAWSON, C. L., (JPL)
!           HANSON, R. J., (SNLA)
!           KINCAID, D. R., (U. OF TEXAS)
!           KROGH, F. T., (JPL)
!***PURPOSE  INTERCHANGE S.P VECTORS
!***DESCRIPTION

!                B L A S  SUBPROGRAM
!    DESCRIPTION OF PARAMETERS

!     --INPUT--
!        N  NUMBER OF ELEMENTS IN INPUT VECTOR(S)
!       SX  SINGLE PRECISION VECTOR WITH N ELEMENTS
!     INCX  STORAGE SPACING BETWEEN ELEMENTS OF SX
!       SY  SINGLE PRECISION VECTOR WITH N ELEMENTS
!     INCY  STORAGE SPACING BETWEEN ELEMENTS OF SY

!     --OUTPUT--
!       SX  INPUT VECTOR SY (UNCHANGED IF N <= 0)
!       SY  INPUT VECTOR SX (UNCHANGED IF N <= 0)

!     INTERCHANGE SINGLE PRECISION SX AND SINGLE PRECISION SY.
!     FOR I = 0 TO N-1, INTERCHANGE  SX(LX+I*INCX) AND SY(LY+I*INCY),
!     WHERE LX = 1 IF INCX >= 0, ELSE LX = (-INCX)*N, AND LY IS
!     DEFINED IN A SIMILAR WAY USING INCY.
!***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
!                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
!                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
!                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  SSWAP



INTEGER                       :: N
INTEGER                       :: INCX
INTEGER                       :: INCY
REAL(EB) SX(1),SY(1),STEMP1,STEMP2,STEMP3
INTEGER :: I, IX, IY, M,  MP1, NS

!***FIRST EXECUTABLE STATEMENT  SSWAP
IF(N<=0)RETURN
IF(INCX==INCY) THEN
  IF(INCX-1<0) GO TO 5
  IF(INCX-1==0) GO TO 20
  IF(INCX-1>0) GO TO 60
END IF
5 CONTINUE

!       CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.

IX = 1
IY = 1
IF(INCX<0)IX = (-N+1)*INCX + 1
IF(INCY<0)IY = (-N+1)*INCY + 1
DO  I = 1,N
  STEMP1 = SX(IX)
  SX(IX) = SY(IY)
  SY(IY) = STEMP1
  IX = IX + INCX
  IY = IY + INCY
END DO
RETURN

!       CODE FOR BOTH INCREMENTS EQUAL TO 1


!       CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 3.

20 M = MOD(N,3)
IF( M == 0 ) GO TO 40
DO  I = 1,M
  STEMP1 = SX(I)
  SX(I) = SY(I)
  SY(I) = STEMP1
END DO
IF( N < 3 ) RETURN
40 MP1 = M + 1
DO  I = MP1,N,3
  STEMP1 = SX(I)
  STEMP2 = SX(I+1)
  STEMP3 = SX(I+2)
  SX(I) = SY(I)
  SX(I+1) = SY(I+1)
  SX(I+2) = SY(I+2)
  SY(I) = STEMP1
  SY(I+1) = STEMP2
  SY(I+2) = STEMP3
END DO
RETURN
60 CONTINUE

!     CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.

NS = N*INCX
DO  I=1,NS,INCX
  STEMP1 = SX(I)
  SX(I) = SY(I)
  SY(I) = STEMP1
END DO
RETURN
END SUBROUTINE SSWAP


SUBROUTINE H3CSIS(RS,RF,L,LBDCND,TS,TF,M,MBDCND,PS,PF,N,NBDCND,  &
    ELMBDA,LDIMF,MDIMF,IERROR,SAVE,W,HX,HY)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+


REAL(EB)   SAVE(-3:*),W(*),HX(0:*),HY(0:*)
REAL(EB)  RS, RF, TS, TF, PS, PF, ELMBDA
INTEGER                       :: L
INTEGER                       :: LBDCND
INTEGER                       :: M
INTEGER                       :: MBDCND
INTEGER                       :: N
INTEGER                       :: NBDCND
INTEGER                   :: LDIMF
INTEGER                   :: MDIMF
INTEGER                      :: IERROR
INTEGER :: I, IAL, IBL, ICL, IDL, ISL, IAM, IBM, ICM, IDM, ISM, ISVPS, J, IERR1
REAL(EB) :: DR, DRBY2, DRSQR, DT, DTBY2, DTSQR, DP, DPSQR, HXM, HXP, HYM, HYP, SUM, S3

!                               CHECK FOR INVALID INPUT

IERROR = 0

!MCG  IF (RS<0.0) THEN
!MCG      IERROR = IERROR + 1
!MCG      SAVE(IERROR) = 1._EB
!MCG  END IF

IF (RF<=RS) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 2._EB
END IF

IF (L<3) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 3._EB
END IF

IF (LBDCND<1 .OR. LBDCND>6) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 4._EB
END IF

!MCG  IF (RS==0.0 .AND. LBDCND<=4) THEN
!MCG      IERROR = IERROR + 1
!MCG      SAVE(IERROR) = 5._EB
!MCG  END IF

IF (ABS(RS)>=TWO_EPSILON_EB .AND. LBDCND>=5) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 6._EB
END IF

IF (TF<=TS) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 7._EB
END IF

IF (M<3) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 8._EB
END IF

IF (MBDCND<1 .OR. MBDCND>9) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 9._EB
END IF

!MCG  IF (TS==0. .AND. (MBDCND==3.OR.MBDCND==4.OR.
!MCG *    MBDCND==8)) THEN
!MCG      IERROR = IERROR + 1
!MCG      SAVE(IERROR) = 10.
!MCG  END IF

IF (ABS(TF-PI)<=TWO_EPSILON_EB .AND. (MBDCND==2.OR.MBDCND==3.OR. MBDCND==6)) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 11._EB
END IF

IF (PF<=PS) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 12._EB
END IF

IF (N<3) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 13._EB
END IF

IF (NBDCND<0 .OR. NBDCND>4) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 14._EB
END IF

IF (LDIMF<L) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 15._EB
END IF

IF (MDIMF<M) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 16._EB
END IF

IF (MBDCND>4 .AND. (NBDCND==1.OR.NBDCND==2.OR. NBDCND==4)) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 17._EB
END IF

IF (LBDCND>4 .AND. .NOT. (MBDCND==3.OR.MBDCND==6.OR.  &
      MBDCND==8.OR.MBDCND==9)) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 18._EB
END IF

IF (LBDCND>4 .AND. .NOT. (NBDCND==0.OR.NBDCND==3)) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 19._EB
END IF

IF (IERROR/=0) THEN
  RETURN
ELSE
  SAVE(1) = IERROR
END IF

!                               DEFINE GRID PARAMETERS

DR = (RF-RS)/FLOAT(L)
DRBY2 = DR/2._EB
DRSQR = 1._EB/ (DR**2)
DT = (TF-TS)/FLOAT(M)
DTBY2 = DT/2._EB
DTSQR = 1._EB/ (DT**2)
DP = (PF-PS)/FLOAT(N)
DPSQR = 1._EB/ (DP**2)

!                               ALLOCATE SAVE ARRAY

IAL = 13
IBL = IAL + L
ICL = IBL + L
IDL = ICL + L
ISL = IDL + L
IAM = ISL + L
IBM = IAM + M
ICM = IBM + M
IDM = ICM + M
ISM = IDM + M
ISVPS = ISM + M

!                               DEFINE AL,BL,CL,DL,SL ARRAYS IN
!                               ARRAY SAVE.  SL IS SYMMETRIZER FOR
!                               R-OPERATOR.
DO  I = 1,L
  HXM = .5_EB*(HX(I-1)+HX(I))
  HXP = .5_EB*(HX(I)+HX(I+1))
!MCG      RI = RS + (I-.5)*DR
!MCG      SAVE(IDL+I) = 1._EB/RI**2
  SAVE(IDL+I) = 1._EB
!MCG      SAVE(IAL+I) = DRSQR* (RI-DRBY2)**2
!MCG      SAVE(ICL+I) = DRSQR* (RI+DRBY2)**2
  SAVE(IAL+I) = DRSQR/(HX(I)*HXM)
  SAVE(ICL+I) = DRSQR/(HX(I)*HXP)
  SAVE(IBL+I) = - (SAVE(IAL+I)+SAVE(ICL+I)) + ELMBDA/SAVE(IDL+I)
!MCG      SAVE(ISL+I) = RI**2
  SAVE(ISL+I) = 1._EB
END DO

!                               DEFINE BOUNDARY COEFFICIENTS

SELECT CASE(LBDCND)
CASE(1:2) ; GO TO 110
CASE(3:6) ; GO TO 120
END SELECT

110 CONTINUE
SAVE(IBL+1) = SAVE(IBL+1) - SAVE(IAL+1)

GO TO 130

120 CONTINUE
SAVE(IBL+1) = SAVE(IBL+1) + SAVE(IAL+1)
130 CONTINUE

SELECT CASE(LBDCND)
CASE(1)   ; GO TO 140
CASE(2:3) ; GO TO 150
CASE(4:5) ; GO TO 140
CASE(6)   ; GO TO 150
END SELECT

140 CONTINUE
SAVE(IBL+L) = SAVE(IBL+L) - SAVE(ICL+L)

GO TO 160

150 CONTINUE
SAVE(IBL+L) = SAVE(IBL+L) + SAVE(ICL+L)
160 CONTINUE

!                               DEFINE ARRAYS AM,BM,CM,DM,SM.  SM IS
!                               SYMMETRIZER FOR THETA-OPERATOR.
DO  J = 1,M
  HYM = .5_EB*(HY(J-1)+HY(J))
  HYP = .5_EB*(HY(J)+HY(J+1))
!MCG      TJ = TS + (J-.5)*DT
!MCG      SAVE(ISM+J) = SIN(TJ)
  SAVE(ISM+J) = 1._EB
!MCG      SAVE(IDM+J) = DPSQR/SIN(TJ)**2
  SAVE(IDM+J) = DPSQR
!MCG      SAVE(IAM+J) = DTSQR*SIN(TJ-DTBY2)/SIN(TJ)
!MCG      SAVE(ICM+J) = DTSQR*SIN(TJ+DTBY2)/SIN(TJ)
  SAVE(IAM+J) = DTSQR/(HY(J)*HYM)
  SAVE(ICM+J) = DTSQR/(HY(J)*HYP)
  SAVE(IBM+J) = - (SAVE(IAM+J)+SAVE(ICM+J))
END DO

!                               DEFINE BOUNDARY COEFFICIENTS

SELECT CASE(MBDCND)
CASE(1:2) ; GO TO 180
CASE(3:6) ; GO TO 190
CASE(7)   ; GO TO 180
CASE(8:9) ; GO TO 190
END SELECT

180 CONTINUE
SAVE(IBM+1) = SAVE(IBM+1) - SAVE(IAM+1)

GO TO 200

190 CONTINUE
SAVE(IBM+1) = SAVE(IBM+1) + SAVE(IAM+1)
200 CONTINUE

SELECT CASE(MBDCND)
CASE(1)   ; GO TO 210
CASE(2:3) ; GO TO 220
CASE(4:5) ; GO TO 210
CASE(6:9) ; GO TO 220
END SELECT

210 CONTINUE
SAVE(IBM+M) = SAVE(IBM+M) - SAVE(ICM+M)

GO TO 230

220 CONTINUE
SAVE(IBM+M) = SAVE(IBM+M) + SAVE(ICM+M)
230 CONTINUE

!                               INITIALIZE SOLVER ROUTINE S3CCIS

CALL S3CCIS(L,SAVE(IAL+1),SAVE(IBL+1),SAVE(ICL+1), M,SAVE(IAM+1:IAM+M),  &
    SAVE(IBM+1),SAVE(ICM+1),SAVE(IDM+1),NBDCND,N,LDIMF,  &
    MDIMF,IERR1,SAVE(ISVPS+1),W)

!                               TEST ERROR FLAG FROM S3CCIS FOR
!                               INTERNAL ERROR

IF (IERR1/=0) THEN
  SAVE(1) = 99._EB
  IERROR = 1
  RETURN
END IF

!                               SCALE RADIAL COEFFICIENTS
DO  I = 1,L
  SAVE(IAL+I) = SAVE(IAL+I)*SAVE(IDL+I)
  SAVE(ICL+I) = SAVE(ICL+I)*SAVE(IDL+I)
  SAVE(IBL+I) = SAVE(IBL+I)*SAVE(IDL+I)
END DO

!                               COMPUTE SCALING FOR SINGULAR PROBLEMS

SUM = 0._EB
DO  J = 1,M
!MCG      SUM = SUM + SAVE(ISM+J)
  SUM = SUM + SAVE(ISM+J)*HY(J)
END DO

S3 = 0._EB
DO  I = 1,L
!MCG      S3 = S3 + SAVE(ISL+I)
  S3 = S3 + SAVE(ISL+I)*HX(I)
END DO
SAVE(11) = SUM*S3*N

!                               RESTORE ARRAY DM

DO  J = 1,M
  SAVE(IDM+J) = SAVE(IDM+J)/DPSQR
END DO

!                               SAVE PARAMETERS FOR HS3SPH IN SAVE ARRAY

SAVE(2) = DR
SAVE(3) = REAL(L,EB)
SAVE(4) = REAL(LBDCND,EB)
SAVE(5) = DT
SAVE(6) = REAL(M,EB)
SAVE(7) = REAL(MBDCND,EB)
SAVE(8) = DP
SAVE(9) = REAL(N,EB)
SAVE(10) = REAL(NBDCND,EB)
SAVE(12) = ELMBDA

SAVE(-1) = REAL(KAPPA,EB)
SAVE(-2) = REAL(NMAX,EB)
SAVE(-3) = REAL(IKPWR,EB)

RETURN
END SUBROUTINE H3CSIS


SUBROUTINE H3CSSS(BDRS,BDRF,BDTS,BDTF,BDPS,BDPF,LDIMF,MDIMF,F,  &
    PERTRB,SAVE,W,HX,HY)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

INTEGER                   :: LDIMF
INTEGER                   :: MDIMF
REAL(EB) BDRS(MDIMF,*), BDRF(MDIMF,*), BDTS(LDIMF,*),BDTF(LDIMF,*),BDPS(LDIMF,*),BDPF(LDIMF,*)
REAL(EB) F(LDIMF,MDIMF,*)
REAL(EB) PERTRB
REAL(EB) SAVE(-3:*)
REAL(EB) W(*)
REAL(EB) HX(0:*)
REAL(EB) HY(0:*)
INTEGER :: I, L, LBDCND, M, MBDCND, N, NBDCND, IAL, IBL, ICL, IDL, ISL, IAM, IBM, ICM, IDM, ISM, ISVPS, K, J, ISING
REAL(EB) :: DR, DT, DP, DPR, DPSQR, ELMBDA, SUM, PERT, PRTSAV, SCAL

!                               CHECK VALUE OF IERROR (=SAVE(1)).
!                               IF NON-ZERO, RETURN.

IF (ABS(SAVE(1))>=TWO_EPSILON_EB) RETURN

!                               GET PARAMETERS FOR H3CSSS FROM SAVE
!                               ARRAY WHERE THEY WERE STORED IN
!                               INITIALIZATION SUBROUTINE H3CSIS.

KAPPA = NINT(SAVE(-1))   ! EXTRA VARIABLES ADDED TO SAVE ARRAY
NMAX  = NINT(SAVE(-2))
IKPWR = NINT(SAVE(-3))

DR     = SAVE(2)
L      = NINT(SAVE(3))
LBDCND = NINT(SAVE(4))
DT     = SAVE(5)
M      = NINT(SAVE(6))
MBDCND = NINT(SAVE(7))
DP     = SAVE(8)
DPR    = 1._EB/DP
DPSQR  = 2._EB/ (DP**2)
N      = NINT(SAVE(9))
NBDCND = NINT(SAVE(10))
ELMBDA = SAVE(12)

!                               ALLOCATE SAVE ARRAY

IAL = 13
IBL = IAL + L
ICL = IBL + L
IDL = ICL + L
ISL = IDL + L
IAM = ISL + L
IBM = IAM + M
ICM = IBM + M
IDM = ICM + M
ISM = IDM + M
ISVPS = ISM + M

!                               ENTER BOUNDARY DATA FOR R-BOUNDARIES

SELECT CASE(LBDCND)
CASE(1:2) ; GO TO 180
CASE(3:4) ; GO TO 210
CASE(5:6) ; GO TO 240
END SELECT

180 CONTINUE
DO  K = 1,N
  DO  J = 1,M
    F(1,J,K) = F(1,J,K) - 2._EB*SAVE(IAL+1)*BDRS(J,K)
  END DO
END DO

GO TO 240

210 CONTINUE
DO  K = 1,N
  DO  J = 1,M
    F(1,J,K) = F(1,J,K) + DR*SAVE(IAL+1)*BDRS(J,K)
  END DO
END DO

240 CONTINUE
SELECT CASE(LBDCND)
CASE(1)   ; GO TO 250
CASE(2:3) ; GO TO 280
CASE(4:5) ; GO TO 250
CASE(6)   ; GO TO 280
END SELECT

250 CONTINUE
DO  K = 1,N
  DO  J = 1,M
    F(L,J,K) = F(L,J,K) - 2._EB*SAVE(ICL+L)*BDRF(J,K)
  END DO
END DO

GO TO 310

280 CONTINUE
DO  K = 1,N
  DO  J = 1,M
    F(L,J,K) = F(L,J,K) - DR*SAVE(ICL+L)*BDRF(J,K)
  END DO
END DO
310 CONTINUE

!                               ENTER BOUNDARY DATA FOR THETA-BOUNDARIES

SELECT CASE(MBDCND)
CASE(1:2) ; GO TO 320
CASE(3:4) ; GO TO 350
CASE(5:6) ; GO TO 380
CASE(7)   ; GO TO 320
CASE(8)   ; GO TO 350
CASE(9)   ; GO TO 380
END SELECT

320 CONTINUE
DO  K = 1,N
  DO  I = 1,L
    F(I,1,K) = F(I,1,K) - 2._EB*SAVE(IAM+1)*SAVE(IDL+I)*BDTS(I,K)
  END DO
END DO

GO TO 380

350 CONTINUE
SCAL = DT*SAVE(IAM+1)
DO  K = 1,N
  DO  I = 1,L
    F(I,1,K) = F(I,1,K) + SCAL*SAVE(IDL+I)*BDTS(I,K)
  END DO
END DO

380 CONTINUE
SELECT CASE(MBDCND)
CASE(1)   ; GO TO 390
CASE(2:3) ; GO TO 420
CASE(4:5) ; GO TO 390
CASE(6)   ; GO TO 420
CASE(7:9) ; GO TO 450
END SELECT

390 CONTINUE
DO  K = 1,N
  DO  I = 1,L
    F(I,M,K) = F(I,M,K) - 2._EB*SAVE(ICM+M)*SAVE(IDL+I)*BDTF(I,K)
  END DO
END DO

GO TO 450

420 CONTINUE
SCAL = DT*SAVE(ICM+M)
DO  K = 1,N
  DO  I = 1,L
    F(I,M,K) = F(I,M,K) - SCAL*SAVE(IDL+I)*BDTF(I,K)
  END DO
END DO
450 CONTINUE

!                               ENTER BOUNDARY DATA FOR PHI-BOUNDARIES

SELECT CASE(NBDCND+1)
CASE(1)   ; GO TO 590
CASE(2:3) ; GO TO 460
CASE(4:5) ; GO TO 490
END SELECT

460 CONTINUE
DO  J = 1,M
  SCAL = DPSQR*SAVE(IDM+J)
  DO  I = 1,L
    F(I,J,1) = F(I,J,1) - SCAL*SAVE(IDL+I)*BDPS(I,J)
  END DO
END DO

GO TO 520

490 CONTINUE
DO  J = 1,M
  SCAL = DPR*SAVE(IDM+J)
  DO  I = 1,L
    F(I,J,1) = F(I,J,1) + SCAL*SAVE(IDL+I)*BDPS(I,J)
  END DO
END DO

520 CONTINUE
SELECT CASE(NBDCND+1)
CASE(1)   ; GO TO 590
CASE(2)   ; GO TO 530
CASE(3:4) ; GO TO 560
CASE(5)   ; GO TO 530
END SELECT

530 CONTINUE
DO  J = 1,M
  SCAL = DPSQR*SAVE(IDM+J)
  DO  I = 1,L
    F(I,J,N) = F(I,J,N) - SCAL*SAVE(IDL+I)*BDPF(I,J)
  END DO
END DO

GO TO 590

560 CONTINUE
DO  J = 1,M
  SCAL = DPR*SAVE(IDM+J)
  DO  I = 1,L
    F(I,J,N) = F(I,J,N) - SCAL*SAVE(IDL+I)*BDPF(I,J)
  END DO
END DO
590 CONTINUE

PERTRB = 0._EB
ISING = 0

!                               FOR SINGULAR PROBLEMS ADJUST DATA TO
!                               INSURE A SOLUTION WILL EXIST.  GO THRU
!                               THIS CODE TWICE: ISING=1 FOR CALCULATING
!                               PERTRB; ISING=2 FOR NORMALIZING SOLUTION
!                               AFTER IT IS COMPUTED.

SELECT CASE(LBDCND)
CASE(1:2) ; GO TO 770
CASE(3)   ; GO TO 620
CASE(4:5) ; GO TO 770
CASE(6)   ; GO TO 620
END SELECT

620 CONTINUE
SELECT CASE(MBDCND)
CASE(1:2) ; GO TO 770
CASE(3)   ; GO TO 630
CASE(4:5) ; GO TO 770
CASE(6)   ; GO TO 630
CASE(7)   ; GO TO 770
CASE(8:9) ; GO TO 630
END SELECT

630 CONTINUE
SELECT CASE(NBDCND+1)
CASE(1)   ; GO TO 640
CASE(2:3) ; GO TO 770
CASE(4)   ; GO TO 640
CASE(5)   ; GO TO 770
END SELECT

640 CONTINUE
IF (ABS(ELMBDA)<=TWO_EPSILON_EB) THEN
  GO TO 650
ELSE
  GO TO 770
END IF
650 CONTINUE
ISING = 1
660 CONTINUE
DO  I = 1,L
  W(I) = 0._EB
END DO
DO  I = 1,L
  DO  J = 1,M
    DO  K = 1,N
!MCG              W(I) = W(I) + SAVE(ISM+J)*F(I,J,K)
      W(I) = W(I) + HY(J)*SAVE(ISM+J)*F(I,J,K)
    END DO
  END DO
END DO

SUM = 0._EB
DO  I = 1,L
!MCG      SUM = SUM + W(I)*SAVE(ISL+I)
  SUM = SUM + W(I)*SAVE(ISL+I)*HX(I)
END DO

PERT = SUM/SAVE(11)
!                               ADJUST F ARRAY BY PERT

DO  K = 1,N
  DO  J = 1,M
    DO  I = 1,L
      F(I,J,K) = F(I,J,K) - PERT
    END DO
  END DO
END DO

!                               IF NORMALIZING SOLUTION, RESTORE PERTRB
!                               AND JUMP TO END

IF (ISING==2) THEN
  PERTRB = PRTSAV

  GO TO 850

END IF

PRTSAV = PERT

770 CONTINUE

!                               SCALE RIGHT SIDE OF EQUATION BEFORE
!                               CALL TO S3CCSS

DO  K = 1,N
  DO  J = 1,M
    DO  I = 1,L
      F(I,J,K) = F(I,J,K)*SAVE(ISL+I)
    END DO
  END DO
END DO

!                               SOLVE SYSTEM USING S3CCSS

CALL S3CCSS(LDIMF,MDIMF,F,SAVE(ISVPS+1),W)

!                               IF A SINGULAR PROBLEM,
!                               RE-NORMALIZE SOLUTION (ISING=2)

IF (ISING==1) THEN
  ISING = 2

  GO TO 660

END IF

850 CONTINUE

RETURN
END SUBROUTINE H3CSSS


SUBROUTINE S3CCIS(L,AL,BL,CL,M,AM,BM,CM,DM,NPEROD,N,LDIMF,MDIMF,  &
    IERROR,SAVE,W)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+



INTEGER                       :: M
INTEGER                       :: NPEROD
INTEGER                   :: N
INTEGER                   :: LDIMF
INTEGER                   :: MDIMF
INTEGER                      :: IERROR
INTEGER :: L, J, IGRID
REAL(EB) AL(L),BL(L),CL(L),AM(M),BM(M),CM(M),DM(M),SAVE(-3:*),W(*)

!                               CHECK FOR INVALID INPUT

IERROR = 0

IF (L<2) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 1._EB
END IF

IF (M<3) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 2._EB
END IF

IF (LDIMF<L) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 3._EB
END IF

DO  J = 2,M

  IF (AM(J)*CM(J-1)<0.0_EB) THEN
    IERROR = IERROR + 1
    SAVE(IERROR) = 5._EB

    EXIT

  END IF

END DO

IF (MDIMF<M) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 6._EB
END IF

IF (NPEROD<0 .AND. NPEROD>4) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 7._EB
END IF

IF (N<3) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 8._EB
END IF

IF (IERROR/=0) THEN
  RETURN
ELSE
  SAVE(1) = IERROR
END IF

IGRID = 2

CALL FSH21S(IGRID,L,AL,BL,CL,M,AM,BM,CM,DM,NPEROD,N,IERROR,SAVE,W)

RETURN
END SUBROUTINE S3CCIS


SUBROUTINE S3CCSS(LDIMF,MDIMF,F,SAVE,W)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+



INTEGER                   :: LDIMF
INTEGER                   :: MDIMF
REAL(EB)   SAVE(-3:*)
REAL(EB)   W(*)
REAL(EB)   F(LDIMF,MDIMF,*)


!                               CHECK VALUE OF IERROR (=SAVE(1)).
!                               IF NON-ZERO, RETURN.

IF (ABS(SAVE(1))>=TWO_EPSILON_EB) RETURN

CALL FSH15S(LDIMF,MDIMF,F,SAVE,W)

RETURN
END SUBROUTINE S3CCSS


SUBROUTINE FSH15S(LDIMY,MDIMY,Y,SAVE,W)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+



INTEGER                   :: LDIMY
INTEGER                       :: MDIMY
REAL(EB)    SAVE(-3:*)
REAL(EB)    W(*)
REAL(EB)   Y(*)
INTEGER :: L, M, N, NP, IGRID, IAL, IBL, ICL, IAM, ICM, ICFZ, IWSZ, IB, ICF, LENY

!                               SOLVER BASED ON CYCLIC REDUCTION AND
!                               FAST FOURIER TRANSFORMS

L  = NINT(SAVE(2))
M  = NINT(SAVE(4))
N  = NINT(SAVE(5))
NP = NINT(SAVE(6))
IGRID = NINT(SAVE(7))

!                               ALLOCATE SAVE ARRAY

IAL = 8
IBL = IAL + L
ICL = IBL + L
IAM = ICL + L
ICM = IAM + M
ICFZ = ICM + M

IF (IGRID==1) THEN
  IWSZ = ICFZ + 2*N
ELSE
  IWSZ = ICFZ + 4*N
END IF

IB = IWSZ + N + 16
ICF = IB + ((KAPPA-2)*IKPWR+KAPPA+6)*N

IF (LDIMY==L .AND. MDIMY==M) THEN

!                               DATA ARRAY HAS NO HOLES, SO CALL SOLVER

  CALL FSH16S(IGRID,L,M,N,NP,SAVE(IAL),SAVE(IBL),SAVE(ICL),  &
      SAVE(IAM),SAVE(ICM),SAVE(ICFZ),SAVE(IWSZ),  &
      SAVE(IB),SAVE(ICF),Y,W,W(1+M),W(1+M+L*M))

ELSE

!                               PACK DATA ARRAY, CALL SOLVER, AND UNPACK

  CALL FSH04S(L,M,N,LDIMY,MDIMY,L,Y,W)

  IF (N>1) THEN
    LENY = L*M* (N+1)
  ELSE
    LENY = L*M
  END IF

  CALL FSH16S(IGRID,L,M,N,NP,SAVE(IAL),SAVE(IBL),SAVE(ICL),  &
      SAVE(IAM),SAVE(ICM),SAVE(ICFZ),SAVE(IWSZ),  &
      SAVE(IB),SAVE(ICF),W,W(1+LENY),W(1+LENY+M),Y)

  CALL FSH05S(L,M,N,LDIMY,MDIMY,L,Y,W)

END IF

RETURN
END SUBROUTINE FSH15S


SUBROUTINE FSH16S(IGRID,L,M,N,NP,AL,BL,CL,AM,CM,CFZ,WSAVEZ,B, COEF,F,W1,W2,FT)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+



INTEGER                   :: IGRID
INTEGER                       :: L
INTEGER                       :: M
INTEGER                       :: N
INTEGER                   :: NP
REAL(EB)   AL(L),BL(L), CL(L)
REAL(EB)    AM(M)
REAL(EB)    CM(M)
REAL(EB)   CFZ(4*N)
REAL(EB)    WSAVEZ(N+16)
REAL(EB)    B(*)
REAL(EB)   COEF(*)
REAL(EB)    F(L,M,N)
REAL(EB)    W1(M)
REAL(EB)    W2(L*M)
REAL(EB)   FT(L,M,N)
INTEGER :: K,J,I
INTEGER :: IFWRD, LDIMFT, IB, ICF

LOGICAL :: DATARY

!                               BEGIN SOLUTION

NOCOPY=.TRUE.
DATARY=.TRUE.
SCALE=1._EB
IFWRD = 1
TPOSE=.FALSE.
LDIMFT=L
100 CONTINUE

IF (N/=1) THEN

!                               TRANSFORM IN Z

  IF (DATARY) THEN
    CALL FSH26S(IGRID,IFWRD,NP,L,N,M,LDIMFT,F,FT,CFZ,WSAVEZ)
    DATARY=OUTARY
  ELSE
    CALL FSH26S(IGRID,IFWRD,NP,L,N,M,LDIMFT,FT,F,CFZ,WSAVEZ)
    DATARY=.NOT.OUTARY
  END IF

END IF

IF (IFWRD/=1) GO TO 900
IB = 1
ICF = 1
DO  K = 1,N
  IF (DATARY) THEN
    CALL FSH17S(AM,CM,L,AL,BL,CL,L,F(1,1,K), B(IB),COEF(ICF),W1,W2,FT)
  ELSE
    CALL FSH17S(AM,CM,L,AL,BL,CL,L,FT(1,1,K), B(IB),COEF(ICF),W1,W2,F)
  END IF
  IB = IB + (KAPPA-2)*IKPWR + KAPPA + 6
  ICF = INT(ICF + (2*KAPPA-4.5)*IKPWR + KAPPA + 8)
END DO
IFWRD = 2

GO TO 100

900 CONTINUE

IF (DATARY) THEN
  DO K=1,N
    DO J=1,M
      DO I=1,L
        F(I,J,K)=SCALE*F(I,J,K)
      END DO
    END DO
  END DO
ELSE
  DO K=1,N
    DO J=1,M
      DO I=1,L
        F(I,J,K)=SCALE*FT(I,J,K)
      END DO
    END DO
  END DO
END IF

RETURN
END SUBROUTINE FSH16S


SUBROUTINE FSH17S(AN,CN,M,AM,BM,CM,IDIMY,Y,B,COEF,RT,WS,D)

!                               SUBROUTINE FSH17S SOLVES THE LINEAR
!                               SYSTEM BY VECTORIZED CYCLIC REDUCTION


INTEGER                       :: M
REAL(EB) AM(M)
REAL(EB) BM(M)
REAL(EB) CM(M)
INTEGER                   :: IDIMY
REAL(EB)   Y(IDIMY,NMAX)
REAL(EB)    B(*)
REAL(EB)   COEF(*)
REAL(EB)   RT(NMAX)
REAL(EB)    WS(NMAX,M)
REAL(EB)    D(NMAX,M)
REAL(EB)   AN(NMAX),CN(NMAX)
INTEGER :: I, KDO, IF, IC, IR, IRM1, I2, I3, I4, IS, IM2, NM2, IP, J, IPI2, IMI2, IP2, NP2, I1, IMI1, IPI1, IM1, NM1
INTEGER :: IP1, NP1, IZ, NZ

!                               LET KAPPA = LOG2(NMAX) + 1, THEN
! LENGTH OF B ARRAY = (KAPPA-2)*2**(KAPPA+1) + KAPPA + 5
! LENGTH OF COEF ARRAY = (2*KAPPA-4.5)*2**(KAPPA+1) + KAPPA + 8


!                               BEGIN REDUCTION PHASE

KDO = KAPPA - 1
IF = 2**KAPPA
IC = 0
DO  IR = 0,KAPPA - 2
  IRM1 = IR - 1
  I2 = 2**IR
  I3 = I2 + I2/2
  I4 = I2 + I2
  IS = 0
  DO  I = I2,NMAX,I4
    CALL FSH10S(I,IR,IM2,NM2)
    DO  IP = 1,NM2
      RT(IS+IP) = B(IM2+IP-1)
    END DO
    DO  IP = 1,NM2
      DO  J = 1,M
        WS(IS+IP,J) = Y(J,I)
      END DO
    END DO
    IS = IS + NM2
  END DO
  CALL FSH18S(IS,M,AM,BM,CM,RT,WS,D)
  IS = 0
  DO  I = I4,NMAX,I4
    IPI2 = I + I2
    IMI2 = I - I2
    CALL FSH10S(IMI2,IR,IM2,NM2)
    DO  IP = 1,NM2
      DO  J = 1,M
        Y(J,I) = Y(J,I) + COEF(IC+IP)*WS(IS+IP,J)
      END DO
    END DO
    IC = IC + NM2

    IF (IPI2>NMAX) CYCLE
    IS = IS + NM2
    CALL FSH10S(IPI2,IR,IP2,NP2)
    DO  IP = 1,NP2
      DO  J = 1,M
        Y(J,I) = Y(J,I) + COEF(IC+IP)*WS(IS+IP,J)
      END DO
    END DO
    IC = IC + NP2
  END DO
END DO

!                               BEGIN BACK SUBSTITUTION PHASE

DO  IR = KDO,0,-1
  IRM1 = IR - 1
  I2 = 2**IR
  I1 = I2/2
  I4 = I2 + I2

  IF (IR==KDO) GO TO 430

  IF (IR/=0) GO TO 280
  DO  I = 1,NMAX,2

    IF (I/=1) GO TO 230
    DO  J = 1,M
      Y(J,I) = Y(J,I) - CN(I)*Y(J,I+1)
    END DO

    CYCLE

    230         CONTINUE
    IF (I/=NMAX) GO TO 250
    DO  J = 1,M
      Y(J,I) = Y(J,I) - AN(I)*Y(J,I-1)
    END DO

    CYCLE

    250         CONTINUE
    DO  J = 1,M
      Y(J,I) = Y(J,I) - AN(I)*Y(J,I-1) - CN(I)*Y(J,I+1)
    END DO
  END DO

  GO TO 430

  280     CONTINUE
  IS = 0
  DO  I = I2,NMAX,I4
    IMI1 = I - I1
    IPI1 = I + I1
    IMI2 = I - I2
    IPI2 = I + I2

    IF (I==I2) GO TO 320
    CALL FSH10S(IMI1,IRM1,IM1,NM1)
    DO  IP = 1,NM1
      RT(IS+IP) = B(IM1+IP-1)
    END DO
    DO  IP = 1,NM1
      DO  J = 1,M
        WS(IS+IP,J) = Y(J,IMI2)
      END DO
    END DO
    IS = IS + NM1

    320         CONTINUE
    IF (IPI2>NMAX) CYCLE
    CALL FSH10S(IPI1,IRM1,IP1,NP1)
    DO  IP = 1,NP1
      RT(IS+IP) = B(IP1+IP-1)
    END DO
    DO  IP = 1,NP1
      DO  J = 1,M
        WS(IS+IP,J) = Y(J,IPI2)
      END DO
    END DO
    IS = IS + NP1
  END DO
  CALL FSH18S(IS,M,AM,BM,CM,RT,WS,D)
  IS = 0
  DO  I = I2,NMAX,I4
    IMI1 = I - I1
    IPI1 = I + I1

    IF (I==I2) GO TO 390
    CALL FSH10S(IMI1,IRM1,IM1,NM1)
    DO  IP = 1,NM1
      DO  J = 1,M
        Y(J,I) = Y(J,I) + COEF(IC+IP)*WS(IS+IP,J)
      END DO
    END DO
    IC = IC + NM1
    IS = IS + NM1

    390         CONTINUE
    IF (I+I2>NMAX) CYCLE
    CALL FSH10S(IPI1,IRM1,IP1,NP1)
    DO  IP = 1,NP1
      DO  J = 1,M
        Y(J,I) = Y(J,I) + COEF(IC+IP)*WS(IS+IP,J)
      END DO
    END DO
    IC = IC + NP1
    IS = IS + NP1
  END DO
  430     CONTINUE
  IS = 0
  DO  I = I2,NMAX,I4
    CALL FSH10S(I,IR,IZ,NZ)
    DO  IP = 1,NZ
      RT(IS+IP) = B(IZ+IP-1)
    END DO
    DO  IP = 1,NZ
      DO  J = 1,M
        WS(IS+IP,J) = Y(J,I)
      END DO
    END DO
    IS = IS + NZ
  END DO
  CALL FSH18S(IS,M,AM,BM,CM,RT,WS,D)
  IS = 0
  DO  I = I2,NMAX,I4
    CALL FSH10S(I,IR,IZ,NZ)
    DO  J = 1,M
      Y(J,I) = COEF(IC+1)*WS(IS+1,J)
    END DO
    DO  IP = 2,NZ
      DO  J = 1,M
        Y(J,I) = Y(J,I) + COEF(IC+IP)*WS(IS+IP,J)
      END DO
    END DO
    IC = IC + NZ
    IS = IS + NZ
  END DO
END DO
RETURN
END SUBROUTINE FSH17S


SUBROUTINE FSH18S(IDEG,M,A,B,C,XL,Y,D)



INTEGER                       :: IDEG
INTEGER                       :: M
REAL(EB)   XL(IDEG)
REAL(EB)   Y(NMAX,M)
REAL(EB)   D(NMAX,M)
REAL(EB)   A(M),B(M),C(M)
INTEGER :: I, MM1, K

!                              FSH18S IS A VECTORIZED TRIDIAGONAL SOLVER
!                              THAT ALSO DOES FACTORIZATION

MM1 = M - 1
DO  K = 1,IDEG
  D(K,1) = C(1)/ (B(1)-XL(K))
  Y(K,1) = Y(K,1)/ (B(1)-XL(K))
END DO
DO  I = 2,MM1
  DO  K = 1,IDEG
    D(K,I) = C(I)/ (B(I)-XL(K)-A(I)*D(K,I-1))
    Y(K,I) = (Y(K,I)-A(I)*Y(K,I-1))/ (B(I)-XL(K)-A(I)*D(K,I-1))
  END DO
END DO
DO  K = 1,IDEG
  D(K,M) = Y(K,M) - A(M)*Y(K,MM1)
  Y(K,M) = B(M) - XL(K) - A(M)*D(K,MM1)
END DO
DO  K = 1,IDEG

!                               Y(K,M) = CVMGZ(0.,D(K,M)/Y(K,M),Y(K,M))
!                               ON A CRAY-1

  IF (ABS(Y(K,M))>=TWO_EPSILON_EB) Y(K,M) = D(K,M)/Y(K,M)
END DO
DO  I = M - 1,1,-1
  DO  K = 1,IDEG
    Y(K,I) = Y(K,I) - D(K,I)*Y(K,I+1)
  END DO
END DO
RETURN
END SUBROUTINE FSH18S


SUBROUTINE FSH21S(IGRID,L,AL,BL,CL,M,AM,BM,CM,DM,NPEROD,N,IERROR,SAVE,W)
INTEGER :: IGRID, L , M, N, IERROR, NPEROD
REAL(EB) AL(L),BL(L), CL(L)
REAL(EB)  AM(M)
REAL(EB)  BM(M)
REAL(EB)  CM(M)
REAL(EB)  DM(M)
REAL(EB) SAVE(-3:*)
REAL(EB)  W(*)
INTEGER :: I, ML, NP, IAL, IBL, ICL, IAM, ICM, ICFZ, IWSZ, IB, ICF, NRDEL, K, J, IERR1
REAL(EB) :: DEL


! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+




!                               THIS ROUTINE INITIALIZES SOLVER BASED ON
!                               CYCLIC REDUCTION AND FAST FOURIER TRANS.

!                               COMPUTE SMALLEST INTEGER KAPPA SUCH THAT
!                                         2**KAPPA >= M+1

KAPPA = 2
ML = 4
110 CONTINUE

IF (ML<M+1) THEN
  KAPPA = KAPPA + 1
  ML = 2*ML

  GO TO 110

END IF

IKPWR = 2*ML
ML = ML - 1

NMAX = M
NP = NPEROD + 1

!                               ALLOCATE SAVE ARRAY

IAL = 8
IBL = IAL + L
ICL = IBL + L
IAM = ICL + L
ICM = IAM + M
ICFZ = ICM + M

IF (IGRID==1) THEN
  IWSZ = ICFZ + 2*N
ELSE
  IWSZ = ICFZ + 4*N
END IF

IB = IWSZ + N + 16
ICF = IB + ((KAPPA-2)*IKPWR+KAPPA+6)*N

!                               COMPUTE TRANSFORM ROOTS IN K-DIRECTION
!                               AND STORE IN W(1),...,W(N)

IF (N==1) THEN

  W(1) = 0._EB

ELSE

  IF (IGRID==1) THEN
    NRDEL = ((NP-1)* (NP-3)* (NP-5))/3
    DEL = PI/ (2._EB* (N+NRDEL))
  ELSE
    DEL = PI/ (2*N)
  END IF

  IF (NP==1) THEN
    W(1) = 0._EB
    W(N) = -4._EB
    DO K = 2,N - 1,2
      W(K) = -4._EB*SIN(K*DEL)**2
      W(K+1) = W(K)
    END DO
  END IF

  IF (NP==2) THEN
    DO K = 1,N
      W(K) = -4._EB*SIN(K*DEL)**2
    END DO
  END IF

  IF (NP==3 .OR. NP==5) THEN
    DO K = 1,N
      W(K) = -4._EB*SIN((K-.5_EB)*DEL)**2
    END DO
  END IF

  IF (NP==4) THEN
    DO K = 1,N
      W(K) = -4._EB*SIN((K-1)*DEL)**2
    END DO
  END IF

END IF

!                               INITIALIZE FFT TRANSFORMS AND
!                               PRE-PROCESSING COEFFICIENTS IN K

IF (N/=1) THEN

  SELECT CASE(IGRID)
  CASE(1)   ; GO TO 210
  CASE(2)   ; GO TO 220
END SELECT

210     CONTINUE
SELECT CASE(NP)
CASE(1)   ; GO TO 230
CASE(2)   ; GO TO 240
CASE(3)   ; GO TO 250
CASE(4)   ; GO TO 270
CASE(5)   ; GO TO 280
END SELECT

220     CONTINUE
SELECT CASE(NP)
CASE(1)   ; GO TO 230
CASE(2)   ; GO TO 250
CASE(3)   ; GO TO 260
CASE(4)   ; GO TO 280
CASE(5)   ; GO TO 290
END SELECT

230     CONTINUE
CALL VSRFTI(N,SAVE(IWSZ))

GO TO 300

240     CONTINUE
CALL VSINTI(N,SAVE(ICFZ),SAVE(IWSZ))

GO TO 300

250     CONTINUE
CALL VSCOSI(N,SAVE(ICFZ),SAVE(ICFZ+N),SAVE(IWSZ))

GO TO 300

260     CONTINUE
!EB          CALL VSSNQI(N,SAVE(ICFZ),SAVE(ICFZ+N),SAVE(ICFZ+2*N),
!EB     *                SAVE(ICFZ+3*N),SAVE(IWSZ))
CALL VSCSQI(N,SAVE(ICFZ),SAVE(ICFZ+N),SAVE(ICFZ+2*N),  &
    SAVE(ICFZ+3*N),SAVE(IWSZ))

GO TO 300

270     CONTINUE
CALL VCOSTI(N,SAVE(ICFZ),SAVE(IWSZ))

GO TO 300

280     CONTINUE
CALL VSCOSI(N,SAVE(ICFZ),SAVE(ICFZ+N),SAVE(IWSZ))

GO TO 300

290     CONTINUE
CALL VSCSQI(N,SAVE(ICFZ),SAVE(ICFZ+N),SAVE(ICFZ+2*N),  &
    SAVE(ICFZ+3*N),SAVE(IWSZ))
300     CONTINUE
END IF

DO  K = 1,N

!                               COMPUTE NEW BM ARRAY AND STORE IN
!                               W(N+1),...,W(N+M)
  DO  J = 1,M
    W(N+J) = BM(J) + W(K)*DM(J)
  END DO

!                               SUBROUTINE FSH07S COMPUTES THE ROOTS OF
!                               THE B POLYNOMIALS

  CALL FSH07S(IERR1,AM,W(N+1),CM,SAVE(IB),W(M+N+1), W(3*M+N+1))


  IF (IERR1/=0) THEN
    IERROR = 1
    SAVE(1) = IERR1
    RETURN
  END IF

!                               FSH08S COMPUTES COEFFICIENTS OF PARTIAL
!                               FRACTION EXPANSIONS

  CALL FSH08S(AM,CM,SAVE(IB),SAVE(ICF),W(M+N+1))

  IB = IB + (KAPPA-2)*IKPWR + KAPPA + 6
  ICF = INT(ICF + (2*KAPPA-4.5)*IKPWR + KAPPA + 8)
END DO

!                               SAVE QUANTITIES FOR USE IN SOLVERS

SAVE(2) = REAL(L,EB)
SAVE(4) = REAL(M,EB)
SAVE(5) = REAL(N,EB)
SAVE(6) = REAL(NP,EB)
SAVE(7) = REAL(IGRID,EB)

DO  I = 1,L
  SAVE(IAL+I-1) = AL(I)
  SAVE(IBL+I-1) = BL(I)
  SAVE(ICL+I-1) = CL(I)
END DO

DO  J = 1,M
  SAVE(IAM+J-1) = AM(J)
  SAVE(ICM+J-1) = CM(J)
END DO

RETURN
END SUBROUTINE FSH21S


SUBROUTINE FSH07S(IERROR,AN,BN,CN,B,AH,BH)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!                               FSH07S COMPUTES THE ROOTS OF THE B
!                               POLYNOMIALS USING SUBROUTINE FSH14S
!                               WHICH IS A MODIFICATION THE EISPACK
!                               SUBROUTINE TQLRAT.  IERROR IS SET TO 4
!                               IF EITHER FSH14S FAILS OR IF
!                               A(J+1)*C(J) < 0 FOR SOME J.
!                               AH,BH ARE TEMPORARY WORK ARRAYS.



INTEGER                      :: IERROR, J
REAL(EB) BN(*)
REAL(EB) CN(*)
REAL(EB) B(*)
REAL(EB) AH(*)
REAL(EB) BH(*)
REAL(EB) AN(*)
INTEGER :: IF, KDO, L, IR, I2, I4, IPL, IFD, I, IB, NB, JS, JF, LS, LH
REAL(EB) :: ARG

IERROR = 0
DO  J = 2,NMAX
  ARG = AN(J)*CN(J-1)
  B(J) = SIGN(SQRT(ARG),AN(J))
END DO
IF = 2**KAPPA
KDO = KAPPA - 1
LOOP160:  DO  L = 1,KDO
  IR = L - 1
  I2 = 2**IR
  I4 = I2 + I2
  IPL = I4 - 1
  IFD = IF - I4
  DO  I = I4,IFD,I4
    CALL FSH10S(I,L,IB,NB)
    IF (NB<=0) CYCLE LOOP160
    IF (NB>0) GO TO 110
    110         CONTINUE
    JS = I - IPL
    JF = JS + NB - 1
    LS = 0
    DO  J = JS,JF
      LS = LS + 1
      BH(LS) = BN(J)
      AH(LS) = B(J)
    END DO
    CALL FSH14S(NB,BH,AH,IERROR)
    IF (IERROR==0) THEN
      GO TO 130
    ELSE
      GO TO 190
    END IF
    130         CONTINUE
    LH = IB - 1
    DO  J = 1,NB
      LH = LH + 1
      B(LH) = -BH(J)
    END DO
  END DO
END DO LOOP160
DO  J = 1,NMAX
  B(J) = -BN(J)
END DO
RETURN
190 CONTINUE
IERROR = 4
RETURN
END SUBROUTINE FSH07S


SUBROUTINE FSH08S(AN,CN,B,COEF,T)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!                               SUBROUTINE FSH08S GENERATES COEFFICIENTS
!                               FOR PARTIAL FRACTION EXPANSION OF MATRIX
!                               POLYNOMIAL PRODUCTS


REAL(EB) B(*)
REAL(EB) COEF(*)
REAL(EB) T(NMAX,4)
REAL(EB) AN(NMAX),CN(NMAX), DUM(0:0)
INTEGER :: I, IF, KDO, IS, IR, IRM1, I2, I3, I4, IMI2, IMI3, IDXA, NA, IM2, NM2, IM3
INTEGER :: NM3, IPI2, IPI3, IDXC, NC, IP2, NP2, IP3, NP3, I1, IMI1, IM1, NM1, IPI1, IP1, NP1, IZ, NZ

IF = 2**KAPPA
KDO = KAPPA - 1
IS = 1
DO  IR = 0,KAPPA - 2
  IRM1 = IR - 1
  I2 = 2**IR
  I3 = I2 + I2/2
  I4 = I2 + I2
  DO  I = I4,NMAX,I4
    IMI2 = I - I2
    IMI3 = I - I3
    CALL FSH09S(I,IR,IDXA,NA)
    CALL FSH10S(IMI2,IR,IM2,NM2)
    CALL FSH10S(IMI3,IRM1,IM3,NM3)
    IF (IM3<1) IM3=1
    CALL FSH12S(IR,NA,AN(IDXA),NM2,B(IM2),NM3,B(IM3),0,DUM,  &
        COEF(IS),T,T(1,2),T(1,3),T(1,4))
    IS = IS + NM2
    IPI2 = I + I2

    IF (IPI2>NMAX) CYCLE
    IPI3 = I + I3
    CALL FSH11S(I,IR,IDXC,NC)
    CALL FSH10S(IPI2,IR,IP2,NP2)
    CALL FSH10S(IPI3,IRM1,IP3,NP3)
    IF (IP3<=0) IP3 = 1   ! PREVENT OUT OF BOUNDS
    CALL FSH12S(IR,NC,CN(IDXC),NP2,B(IP2),NP3,B(IP3),0,DUM,  &
        COEF(IS),T,T(1,2),T(1,3),T(1,4))
    IS = IS + NP2
  END DO
END DO

!                               BEGIN BACK SUBSTITUTION PHASE

DO  IR = KDO,0,-1
  IRM1 = IR - 1
  I2 = 2**IR
  I1 = I2/2
  I4 = I2 + I2

  IF (IR==KDO) GO TO 160

  IF (IR==0) GO TO 160
  DO  I = I2,NMAX,I4

    IF (I==I2) GO TO 140
    IMI1 = I - I1
    CALL FSH09S(I,IR,IDXA,NA)
    CALL FSH10S(IMI1,IRM1,IM1,NM1)
    CALL FSH12S(IRM1,NA,AN(IDXA),NM1,B(IM1),0,DUM,0,DUM,  &
        COEF(IS),T,T(1,2),T(1,3),T(1,4))
    IS = IS + NM1

    140         CONTINUE
    IF (I+I2>NMAX) CYCLE
    IPI1 = I + I1
    CALL FSH11S(I,IR,IDXC,NC)
    CALL FSH10S(IPI1,IRM1,IP1,NP1)
    CALL FSH12S(IRM1,NC,CN(IDXC),NP1,B(IP1),0,DUM,0,DUM,  &
        COEF(IS),T,T(1,2),T(1,3),T(1,4))
    IS = IS + NP1
  END DO
  160     CONTINUE
  DO  I = I2,NMAX,I4
    IMI1 = I - I1
    IPI1 = I + I1
    CALL FSH10S(I,IR,IZ,NZ)
    CALL FSH10S(IMI1,IRM1,IM1,NM1)
    CALL FSH10S(IPI1,IRM1,IP1,NP1)
    IF (IM1<1) IM1=1
    IF (IP1<1) IP1=1
    CALL FSH12S(IR,0,DUM,NZ,B(IZ),NM1,B(IM1),NP1,B(IP1),  &
        COEF(IS),T,T(1,2),T(1,3),T(1,4))
    IS = IS + NZ
  END DO
END DO
RETURN
END SUBROUTINE FSH08S


SUBROUTINE FSH09S(I,IR,IDXA,NA)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+
INTEGER :: I, IR, IDXA, NA

NA = 2**IR
IDXA = I - NA + 1

IF (I-NMAX<=0) GO TO 110
IF (I-NMAX>0) GO TO 100
100 CONTINUE
NA = 0
110 CONTINUE
RETURN
END SUBROUTINE FSH09S


SUBROUTINE FSH10S(I,IR,IDX,IDP)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+


!                               B(IDX) IS THE LOCATION OF THE FIRST ROOT
!                               OF THE B(I,IR) POLYNOMIAL

INTEGER :: I, IR, IDX, IDP
INTEGER :: IZH, ID, IPL

IDP = 0
IDX = 0

IF (IR<0) GO TO 160
IF (IR==0) GO TO 100
IF (IR>0) GO TO 120

100 CONTINUE
IF (I-NMAX<=0) GO TO 110
IF (I-NMAX>0) GO TO 160
110 CONTINUE
IDX = I
IDP = 1
RETURN
120 CONTINUE
IZH = 2**IR
ID = I - IZH - IZH
IDX = ID + ID + (IR-1)*IKPWR + IR + (IKPWR-I)/IZH + 4
IPL = IZH - 1
IDP = IZH + IZH - 1

IF (I-IPL-NMAX<=0) GO TO 140
IF (I-IPL-NMAX>0) GO TO 130
130 CONTINUE
IDP = 0
RETURN

140 CONTINUE
IF (I+IPL-NMAX<=0) GO TO 160
IF (I+IPL-NMAX>0) GO TO 150
150 CONTINUE
IDP = NMAX + IPL - I + 1
160 CONTINUE
RETURN
END SUBROUTINE FSH10S


SUBROUTINE FSH11S(I,IR,IDXC,NC)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

INTEGER :: I, IR, IDXC, NC

NC = 2**IR
IDXC = I

IF (IDXC+NC-1-NMAX<=0) GO TO 110
IF (IDXC+NC-1-NMAX>0) GO TO 100
100 CONTINUE
NC = 0
110 CONTINUE
RETURN
END SUBROUTINE FSH11S


SUBROUTINE FSH12S(IR,NA,A,NC,C,NB1,B1,NB2,B2,COEF,W1,W2,W3,W4)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

!                               SUBROUTINE FSH12S ACTUALLY GENERATES THE
!                               COEFFICIENTS


INTEGER                   :: IR
INTEGER                       :: NA
INTEGER                       :: NC
INTEGER                       :: NB1
INTEGER                       :: NB2
REAL(EB) B2(NB2)
REAL(EB) COEF(NC)
REAL(EB) W1(*)
REAL(EB) W2(*)
REAL(EB) W3(*)
REAL(EB) W4(*)
REAL(EB) A(NA),C(NC),B1(NB1)
INTEGER :: NB, IP, M, I, K
REAL(EB) :: SIGN
!                               MAXIMUM LENGTH OF W1, W2, W3, W4 IS NM.


NB = NB1 + NB2
DO  IP = 1,NB1
  W1(IP) = B1(IP)
END DO
DO  IP = 1,NB2
  W1(NB1+IP) = B2(IP)
END DO

!                               MERGE TWO STRINGS

CALL FSH13S(W1,0,NB1,NB1,NB2,NB)
DO  IP = 1,NB
  W4(IP) = W1(NB+IP)
END DO
M = MAX0(NA,NC,NB)
DO  I = 1,NA
  W1(I) = A(I)
END DO
DO  I = NA + 1,M
  W1(I) = 1._EB
END DO
DO  I = NB + 1,M
  W2(I) = 1._EB
END DO
DO  I = NC + 1,M
  W3(I) = 1._EB
END DO
SIGN = 1._EB

IF (IR==0 .AND. NA==1) SIGN = -1._EB
DO  K = 1,NC
  DO  I = 1,NB
    W2(I) = C(K) - W4(I)
  END DO
  DO  I = 1,NC
    W3(I) = C(K) - C(I)
  END DO
  W3(K) = 1._EB
  COEF(K) = SIGN
  DO  I = 1,M
    COEF(K) = COEF(K)*W1(I)*W2(I)/W3(I)
  END DO
END DO
RETURN
END SUBROUTINE FSH12S


SUBROUTINE FSH13S(TCOS,I1,M1,I2,M2,I3)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

INTEGER                   :: I1
INTEGER                       :: M1
INTEGER                   :: I2
INTEGER                       :: M2
INTEGER                       :: I3
REAL(EB)   TCOS(*)
INTEGER :: J1, J2, J, K
REAL(EB) :: X, Y

!                               THIS SUBROUTINE MERGES TWO ASCENDING
!                               STRINGS OF NUMBERS IN THE ARRAY TCOS.
!                               THE FIRST STRING IS OF LENGTH M1 AND
!                               STARTS AT TCOS(I1+1).  THE SECOND
!                               STRING IS OF LENGTH M2 AND STARTS AT
!                               TCOS(I2+1).  THE MERGED STRING GOES INTO
!                               TCOS(I3+1).  IT IS ALMOST A STRAIGHT
!                               COPY FROM THE ORIGINAL FISHPAK VERSION.

J1 = 1
J2 = 1
J = I3

IF (M1+M2==0) GO TO 180

IF (M1==0) GO TO 160

IF (M2==0) GO TO 130
100 CONTINUE
J = J + 1
X = TCOS(I1+J1)
Y = TCOS(I2+J2)

IF (X<=Y) THEN
  TCOS(J) = X
  J1 = J1 + 1

  IF (J1>M1) GO TO 150

  GO TO 100

ELSE
  TCOS(J) = Y
  J2 = J2 + 1
END IF

IF (J2<=M2) GO TO 100

IF (J1>M1) GO TO 180
130 CONTINUE
K = J - J1 + 1

DO  J = J1,M1
  TCOS(K+J) = TCOS(J+I1)
END DO

GO TO 180

150 CONTINUE

IF (J2>M2) GO TO 180
160 CONTINUE
K = J - J2 + 1

DO  J = J2,M2
  TCOS(K+J) = TCOS(J+I2)
END DO
180 CONTINUE
RETURN
END SUBROUTINE FSH13S


SUBROUTINE FSH14S(N,D,E2,IERR)



INTEGER                       :: N
INTEGER                      :: IERR
REAL(EB)   D(N),E2(N)
INTEGER :: I, L, J, M, L1, MML, II, NHALF, NTOP
REAL(EB) :: EPS, F, B, H, C, S, G, P, R, DHOLD

!                               THIS SUBROUTINE IS A MODIFICATION OF THE
!                               EISPACK SUBROUTINE TQLRAT ALGORITHM 464,
!                               COMM. ACM 16, 689(1973) BY REINSCH.

! THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC
! TRIDIAGONAL MATRIX BY THE RATIONAL QL METHOD.

!                               ON INPUT-

!                               N IS THE ORDER OF THE MATRIX,

! D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,

! E2 CONTAINS THE                SUBDIAGONAL ELEMENTS OF THE
! INPUT MATRIX IN ITS LAST N-1 POSITIONS.  E2(1) IS ARBITRARY.

!                               ON OUTPUT-

! D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
! ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND
! ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE
!                               THE SMALLEST EIGENVALUES,

!                               E2 HAS BEEN DESTROYED,

!                               IERR IS SET TO
!                               ZERO       FOR NORMAL RETURN,
! J          IF THE J-TH EIGENVALUE HAS NOT BEEN
!                               DETERMINED AFTER 30 ITERATIONS.

! ********** EPS IS A MACHINE DEPENDENT PARAMETER SPECIFYING
! THE RELATIVE PREC OF FLOATING POINT ARITHMETIC.

IERR = 0

IF (N==1) GO TO 240

EPS = FSH20S()
!EPS = 3.55E-15

DO  I = 2,N
  E2(I-1) = E2(I)*E2(I)
END DO

F = 0.0_EB
B = 0.0_EB
E2(N) = 0.0_EB

DO  L = 1,N
  J = 0
  H = EPS* (ABS(D(L))+SQRT(E2(L)))

  IF (B>H) GO TO 110
  B = H
  C = B*B

!            LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT

  110     CONTINUE
  DO  M = L,N

    IF (E2(M)<=C) GO TO 130

!            E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
!            THROUGH THE BOTTOM OF THE LOOP

  END DO

  130     CONTINUE
  IF (M==L) GO TO 170

  140     CONTINUE
  IF (J==30) GO TO 230
  J = J + 1

!                               FORM SHIFT

  L1 = L + 1
  S = SQRT(E2(L))
  G = D(L)
  P = (D(L1)-G)/ (2.0_EB*S)
  R = SQRT(P*P+1.0_EB)
  D(L) = S/ (P+SIGN(R,P))
  H = G - D(L)

  DO  I = L1,N
    D(I) = D(I) - H
  END DO

  F = F + H

!                               RATIONAL QL TRANSFORMATION

  G = D(M)

  IF (ABS(G)<=TWO_EPSILON_EB) G = B
  H = G
  S = 0.0_EB
  MML = M - L

!                               FOR I=M-1 STEP -1 UNTIL L DO --

  DO  II = 1,MML
    I = M - II
    P = G*H
    R = P + E2(I)
    E2(I+1) = S*R
    S = E2(I)/R
    D(I+1) = H + S* (H+D(I))
    G = D(I) - E2(I)/G

    IF (ABS(G)<=TWO_EPSILON_EB) G = B
    H = G*P/R
  END DO

  E2(L) = S*G
  D(L) = H

!                               GUARD AGAINST UNDERFLOWED H

  IF (ABS(H)<=TWO_EPSILON_EB) GO TO 170

  IF (ABS(E2(L))<=ABS(C/H)) GO TO 170
  E2(L) = H*E2(L)

  IF (ABS(E2(L))>=TWO_EPSILON_EB) GO TO 140
  170     CONTINUE
  P = D(L) + F

!                               ORDER EIGENVALUES

  IF (L==1) GO TO 190

!                               FOR I=L STEP -1 UNTIL 2 DO --

  DO  II = 2,L
    I = L + 2 - II

    IF (P>=D(I-1)) GO TO 200
    D(I) = D(I-1)
  END DO

  190     CONTINUE
  I = 1
  200     CONTINUE
  D(I) = P
END DO

IF (ABS(D(N))>=ABS(D(1))) GO TO 240
NHALF = N/2
DO  I = 1,NHALF
  NTOP = N - I
  DHOLD = D(I)
  D(I) = D(NTOP+1)
  D(NTOP+1) = DHOLD
END DO

GO TO 240

!                               SET ERROR -- NO CONVERGENCE TO AN
!                               EIGENVALUE AFTER 30 ITERATIONS

230 CONTINUE
IERR = L
240 CONTINUE
RETURN
END SUBROUTINE FSH14S


SUBROUTINE H2CZIS(XS,XF,L,LBDCND,YS,YF,M,MBDCND,H, ELMBDA,LDIMF,IERROR,SAVE)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+


REAL(EB)  XS, XF, YS, YF
INTEGER                       :: L
INTEGER                       :: LBDCND
INTEGER                       :: M
INTEGER                       :: MBDCND
REAL(EB) ELMBDA
INTEGER                   :: LDIMF
INTEGER                      :: IERROR
REAL(EB) H(0:*)
REAL(EB) SAVE(-3:*)
INTEGER :: I
INTEGER :: LP, MP, IA, IB, IC, ID, IS, LPEROD, IR
REAL(EB) :: DX, DY, DLXSQR, DLYSQR, HM, HP



!                               CHECK FOR INVALID INPUT

IERROR = 0

IF (XS>XF) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 1._EB
END IF

IF (L<3) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 2._EB
END IF

IF (LBDCND<0 .OR. LBDCND>4) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 3._EB
END IF

IF (YS>YF) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 4._EB
END IF

IF (M<3) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 5._EB
END IF

IF (MBDCND<0 .OR. MBDCND>4) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 6._EB
END IF

IF (LDIMF<L) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 7._EB
END IF

IF (IERROR/=0) THEN
  RETURN
ELSE
  SAVE(1) = 0._EB
END IF

!                               DEFINE GRID PARAMETERS

DX = (XF-XS)/L
DY = (YF-YS)/M
DLXSQR = 1._EB/DX**2
DLYSQR = 1._EB/DY**2
LP = LBDCND + 1
MP = MBDCND + 1

!                               ALLOCATE SAVE ARRAY

IA = 9
IB = IA + L
IC = IB + L
ID = IC + L
IS = ID + L

!                               DEFINE THE A,B,C,D COEFFICIENTS
!                               IN SAVE ARRAY
DO  I = 0,L - 1
  HM = .5_EB*(H(I)+H(I+1))
  HP = .5_EB*(H(I+1)+H(I+2))
  SAVE(IA+I) = DLXSQR/(H(I+1)*HM)
  SAVE(IC+I) = DLXSQR/(H(I+1)*HP)
  SAVE(IB+I) = - (SAVE(IA+I)+SAVE(IC+I)) + ELMBDA
  SAVE(ID+I) = DLYSQR
END DO

SELECT CASE(LP)
CASE(2)
SAVE(IB) = SAVE(IB) - SAVE(IA)
SAVE(IC-1) = SAVE(IC-1) - SAVE(ID-1)
CASE(3)
SAVE(IB) = SAVE(IB) - SAVE(IA)
SAVE(IC-1) = SAVE(IC-1) + SAVE(ID-1)
CASE(4)
SAVE(IB) = SAVE(IB) + SAVE(IA)
SAVE(IC-1) = SAVE(IC-1) + SAVE(ID-1)
CASE(5)
SAVE(IB) = SAVE(IB) + SAVE(IA)
SAVE(IC-1) = SAVE(IC-1) - SAVE(ID-1)
END SELECT

!                               DETERMINE WHETHER OR NOT BOUNDARY
!                               CONDITION IS PERIODIC IN X

IF (LBDCND==0) THEN
  LPEROD = 0
ELSE
  LPEROD = 1
END IF

!                               INITIALIZE SOLVER ROUTINE S2CFIS.

CALL S2CFIS(LPEROD,L,MBDCND,M,SAVE(IA),SAVE(IB),  &
    SAVE(IC),SAVE(ID),LDIMF,IR,SAVE(IS))

!                               TEST ERROR FLAG FROM S2CFIS FOR
!                               INTERNAL ERROR

IF (IR/=0) THEN
  SAVE(1) = 99._EB
  IERROR = 1

  RETURN
END IF

!                               SAVE PARAMETERS FOR H2CCSS IN SAVE ARRAY

SAVE(2) = DX
SAVE(3) = REAL(L,EB)
SAVE(4) = REAL(LP,EB)
SAVE(5) = DY
SAVE(6) = REAL(M,EB)
SAVE(7) = REAL(MP,EB)
SAVE(8) = ELMBDA

RETURN
END SUBROUTINE H2CZIS


SUBROUTINE H2CZSS(BDXS,BDXF,BDYS,BDYF,LDIMF,F,PERTRB,SAVE,W,H)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

USE GLOBAL_CONSTANTS, ONLY: PRES_METHOD
REAL(EB)  BDXS(*)
REAL(EB)  BDXF(*)
REAL(EB)  BDYS(*)
REAL(EB) BDYF(*)
INTEGER  :: LDIMF
REAL(EB)  F(LDIMF,*)
REAL(EB) PERTRB
REAL(EB)  SAVE(-3:*)
REAL(EB)  W(*)
REAL(EB) H(*)
INTEGER :: I, L, LP, M, MP, IA, IB, IC, ID, IS, J, ISING
REAL(EB) :: DX, DY, ELMBDA, DLYRCP, TWDYSQ, PERT, S1, S3, PRTSAV

!                               CHECK VALUE OF IERROR (=SAVE(1)).
!                               IF NON-ZERO, RETURN.

IF (ABS(SAVE(1))>=TWO_EPSILON_EB) RETURN

!                               GET PARAMETERS FOR H2CCSS FROM SAVE
!                               ARRAY WHERE THEY WERE STORED IN
!                               INITIALIZATION SUBROUTINE H2CCIS.

DX = SAVE(2)
L  = NINT(SAVE(3))
LP = NINT(SAVE(4))
DY = SAVE(5)
M  = NINT(SAVE(6))
MP = NINT(SAVE(7))
ELMBDA = SAVE(8)

DLYRCP = 1._EB/DY
TWDYSQ = 2._EB/DY**2

!                               ALLOCATE SAVE ARRAY

IA = 9
IB = IA + L
IC = IB + L
ID = IC + L
IS = 9 + 4*L


IF (PRES_METHOD /= 'SCARC' .AND. PRES_METHOD /= 'USCARC') THEN

!                               ENTER BOUNDARY DATA FOR X-BOUNDARIES

   IF (LP==2 .OR. LP==3) THEN
     DO J = 1,M
       F(1,J) = F(1,J) - 2._EB*BDXS(J)*SAVE(IA)
     END DO
   END IF

   IF (LP==4 .OR. LP==5) THEN
     DO J = 1,M
       F(1,J) = F(1,J) + SAVE(IA)*DX*BDXS(J)
     END DO
   END IF

   IF (LP==2 .OR. LP==5) THEN
     DO J = 1,M
       F(L,J) = F(L,J) - 2._EB*BDXF(J)*SAVE(ID-1)
     END DO
   END IF

   IF (LP==3 .OR. LP==4) THEN
     DO J = 1,M
       F(L,J) = F(L,J) - SAVE(ID-1)*DX*BDXF(J)
     END DO
   END IF

   !                               ENTER BOUNDARY DATA FOR Y-BOUNDARIES

   IF (MP==2 .OR. MP==3) THEN
     DO I = 1,L
       F(I,1) = F(I,1) - BDYS(I)*TWDYSQ
     END DO
   END IF

   IF (MP==4 .OR. MP==5) THEN
     DO I = 1,L
       F(I,1) = F(I,1) + BDYS(I)*DLYRCP
     END DO
   END IF

   IF (MP==2 .OR. MP==5) THEN
     DO I = 1,L
       F(I,M) = F(I,M) - BDYF(I)*TWDYSQ
     END DO
   END IF

   IF (MP==3 .OR. MP==4) THEN
     DO I = 1,L
       F(I,M) = F(I,M) - BDYF(I)*DLYRCP
     END DO
   END IF

END IF

PERT=0._EB
PERTRB = 0._EB
ISING = 0

!                               FOR SINGULAR PROBLEMS ADJUST DATA TO
!                               INSURE A SOLUTION WILL EXIST.  GO THRU
!                               THIS CODE TWICE: ISING=1 FOR CALCULATING
!                               PERTRB; ISING=2 FOR NORMALIZING SOLUTION
!                               AFTER IT IS COMPUTED.

SELECT CASE(LP)
CASE(1)   ; GO TO 630
CASE(2:3) ; GO TO 750
CASE(4)   ; GO TO 630
CASE(5)   ; GO TO 750
END SELECT

630 CONTINUE
SELECT CASE(MP)
CASE(1)   ; GO TO 640
CASE(2:3) ; GO TO 750
CASE(4)   ; GO TO 640
CASE(5)   ; GO TO 750
END SELECT

640 CONTINUE

IF (ABS(ELMBDA)>=TWO_EPSILON_EB) GO TO 750
ISING = 1
660 CONTINUE
PERT = 0._EB
DO  I = 1,L
  W(I) = 0._EB
END DO
DO  J = 1,M
  DO  I = 1,L
    W(I) = W(I) + F(I,J)
  END DO
END DO
S1 = 0._EB
S3 = 0._EB
DO  I = 1,L
  S3 = S3 + H(I+1)
  S1 = S1 + H(I+1)*W(I)
END DO

S3   = S3*M
PERT = S1/S3

!                               ADJUST F ARRAY BY PERT

DO  J = 1,M
  DO  I = 1,L
    F(I,J) = F(I,J) - PERT
  END DO
END DO
750 CONTINUE

!                               IF NORMALIZING SOLUTION, RESTORE PERTRB
!                               AND JUMP TO END

IF (ISING==2) THEN
  PERTRB = PRTSAV

  GO TO 800

END IF

PRTSAV = PERT

!                               SOLVE THE EQUATION

CALL S2CFSS(LDIMF,F,SAVE(IS),W)

!                               IF A SINGULAR PROBLEM,
!                               RE-NORMALIZE SOLUTION (ISING=2)

IF (ISING==1) THEN
  ISING = 2

  GO TO 660

END IF

800 CONTINUE
RETURN
END SUBROUTINE H2CZSS


SUBROUTINE S2CFIS(LPEROD,L,MPEROD,M,A,B,C,D,LDIMF,IERROR,SAVE)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+

INTEGER                       :: LPEROD
INTEGER                       :: L
INTEGER                       :: MPEROD
INTEGER                   :: M
INTEGER                   :: LDIMF, LDIMFC
INTEGER                      :: IERROR
REAL(EB) SAVE(-3:*),SCAL
REAL(EB) A(L),B(L),C(L), D(L)
INTEGER :: I, IGRID, N, NPEROD


!                               CHECK FOR INVALID INPUT

IERROR = 0

IF (LPEROD/=0 .AND. LPEROD/=1) THEN
  IERROR = 1
  SAVE(1) = 1._EB
END IF

IF (L<3) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 2._EB
END IF

IF (MPEROD<0 .AND. MPEROD>4) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 3._EB
END IF

IF (M<3) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 4._EB
END IF

IF (LDIMF<L) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 5._EB
END IF

IF (LPEROD==0) THEN
  DO  I = 1,L

    IF (ABS(A(I)-A(1))>=TWO_EPSILON_EB) GO TO 110
    IF (ABS(B(I)-B(1))>=TWO_EPSILON_EB) GO TO 110
    IF (ABS(C(I)-A(1))>=TWO_EPSILON_EB) GO TO 110
    IF (ABS(D(I)-D(1))>=TWO_EPSILON_EB) GO TO 110

  END DO

  GO TO 120

  110     CONTINUE
  IERROR = IERROR + 1
  SAVE(IERROR) = 6._EB
END IF

120 CONTINUE

IF (IERROR/=0) THEN
  RETURN
ELSE
  SAVE(1) = IERROR
END IF

LDIMFC=L
IF (LDIMF>L .AND. MOD(L,2)==0) LDIMFC=L+1
IGRID = 2
N = 1
NPEROD = 0
SCAL = 0._EB

CALL FSH00S(IGRID,LPEROD,L,MPEROD,M,NPEROD,N,LDIMFC, SCAL,A,B,C,D,SAVE)

RETURN
END SUBROUTINE S2CFIS


SUBROUTINE S2CFSS(LDIMF,F,SAVE,W)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+



INTEGER :: LDIMF
REAL(EB) :: SAVE(-3:*)
REAL(EB) :: W(*)
REAL(EB) :: F(LDIMF,*)
INTEGER :: M


!                               CHECK VALUE OF IERROR (=SAVE(1)).
!                               IF NON-ZERO, RETURN.

IF (ABS(SAVE(1))>=TWO_EPSILON_EB) RETURN

!                               DEBUG PRINTOUTS UNDER CONTROL OF
!                               PARAMETER LVLPRN

M = NINT(SAVE(4))

CALL FSH02S(LDIMF,M,F,SAVE,W)

RETURN
END SUBROUTINE S2CFSS


SUBROUTINE H2CYIS(RS,RF,L,LBDCND,ZS,ZF,N,NBDCND, ELMBDA,XMU,LDIMF,IERROR,SAVE)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+



REAL(EB)RS
REAL(EB) RF
INTEGER                       :: L
INTEGER                       :: LBDCND
REAL(EB) ZS
REAL(EB) ZF
INTEGER                       :: N
INTEGER                       :: NBDCND
REAL(EB) ELMBDA
REAL(EB) XMU
INTEGER                   :: LDIMF
INTEGER                      :: IERROR
REAL(EB) SAVE(-3:*)
INTEGER :: NP, IA, IB, IC, ID, IS, ISVPS, I, IR
REAL(EB) :: DZ, DZSQR, DR, DRBY2, DRSQR, RI



!                               CHECK FOR INVALID INPUT

IERROR = 0

IF (RS<0.0_EB) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 1._EB
END IF

IF (RF<=RS) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 2._EB
END IF

IF (L<3) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 3._EB
END IF

IF (LBDCND<1 .OR. LBDCND>6) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 4._EB
END IF

IF (ABS(RS)<=TWO_EPSILON_EB .AND. LBDCND<=4) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 5._EB
END IF

IF (ABS(RS)<=TWO_EPSILON_EB .AND. ABS(XMU)>=TWO_EPSILON_EB) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 6._EB
END IF

IF (ZF<=ZS) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 7._EB
END IF

IF (N<3) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 8._EB
END IF

IF (NBDCND<0 .OR. NBDCND>4) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 9._EB
END IF

IF (LDIMF<L) THEN
  IERROR = IERROR + 1
  SAVE(IERROR) = 10._EB
END IF

IF (IERROR/=0) THEN
  RETURN
ELSE
  SAVE(1) = IERROR
END IF

!                               DEFINE GRID PARAMETERS

DZ = (ZF-ZS)/FLOAT(N)
DZSQR = 1._EB/ (DZ**2)
NP = NBDCND + 1

DR = (RF-RS)/FLOAT(L)
DRBY2 = DR/2._EB
DRSQR = 1._EB/ (DR**2)

!                               ALLOCATE SAVE ARRAY

IA = 13
IB = IA + L
IC = IB + L
ID = IC + L
IS = ID + L
ISVPS = IS + L

!                               DEFINE A,B,C,D,S ARRAYS IN ARRAY SAVE.
!                               S ARRAY IS SYMMETRIZER FOR SING. CASES

DO  I = 1,L
  RI = RS + (I-.5_EB)*DR
  SAVE(ID+I) = DZSQR
  SAVE(IA+I) = DRSQR* (RI-DRBY2)/RI
  SAVE(IC+I) = DRSQR* (RI+DRBY2)/RI
  SAVE(IB+I) = - (SAVE(IA+I)+SAVE(IC+I)) + ELMBDA + XMU/RI**2
  SAVE(IS+I) = RI
END DO

!                               DEFINE BOUNDARY COEFFICIENTS

SELECT CASE(LBDCND)
CASE(1:2) ; GO TO 110
CASE(3:6) ; GO TO 120
END SELECT

110 CONTINUE
SAVE(IB+1) = SAVE(IB+1) - SAVE(IA+1)

GO TO 130

120 CONTINUE
SAVE(IB+1) = SAVE(IB+1) + SAVE(IA+1)
130 CONTINUE

SELECT CASE(LBDCND)
CASE(1)   ; GO TO 140
CASE(2:3) ; GO TO 150
CASE(4:5) ; GO TO 140
CASE(6)   ; GO TO 150
END SELECT

140 CONTINUE
SAVE(IB+L) = SAVE(IB+L) - SAVE(IC+L)

GO TO 160

150 CONTINUE
SAVE(IB+L) = SAVE(IB+L) + SAVE(IC+L)
160 CONTINUE

!                               INITIALIZE SOLVER ROUTINE S2CFIS

CALL S2CFIS(1,L,NBDCND,N,SAVE(IA+1),SAVE(IB+1),  &
    SAVE(IC+1),SAVE(ID+1),LDIMF,IR,SAVE(ISVPS+1))

!                               TEST ERROR FLAG FROM S2CFIS FOR
!                               INTERNAL ERROR

IF (IR/=0) THEN
  SAVE(1) = 99._EB
  IERROR = 1

  RETURN
END IF

!                               SAVE PARAMETERS FOR H3CYSS IN SAVE ARRAY

SAVE(2) = DR
SAVE(3) = REAL(L,EB)
SAVE(4) = REAL(LBDCND,EB)
SAVE(6) = DZ
SAVE(7) = REAL(N,EB)
SAVE(8) = REAL(NP,EB)
SAVE(9) = DZSQR
SAVE(10) = 1._EB/DZ
SAVE(11) = ELMBDA
SAVE(12) = XMU

RETURN
END SUBROUTINE H2CYIS


SUBROUTINE H2CYSS(BDRS,BDRF,BDZS,BDZF,LDIMF,F,PERTRB,SAVE,W)

! +--------------------------------------------------------------------+
! |                                                                    |
! |                       COPYRIGHT (C) 1989 BY                        |
! |                          ROLAND A. SWEET                           |
! |                        ALL RIGHTS RESERVED                         |
! |                                                                    |
! +--------------------------------------------------------------------+



REAL(EB)BDRS(*)
REAL(EB) BDRF(*)
REAL(EB) BDZS(*)
REAL(EB)BDZF(*)
INTEGER                   :: LDIMF
REAL(EB)F(LDIMF,*)
REAL(EB)PERTRB
REAL(EB) SAVE(-3:*)
REAL(EB) W(*)
INTEGER :: I, L, LBDCND, N, NP, IA, IB, IC, ID, IS, ISVPS, K, ISING
REAL(EB) :: DR, DZSQR, DZR, ELMBDA, XMU, S3, S1, PERT, PRTSAV



!                               CHECK VALUE OF IERROR (=SAVE(1)).
!                               IF NON-ZERO, RETURN.

IF (ABS(SAVE(1))>=TWO_EPSILON_EB) RETURN

!                               GET PARAMETERS FOR H2CYSS FROM SAVE
!                               ARRAYWHERE THEY WERE STORED IN
!                               INITIALIZATION SUBROUTINE H2CYIS.

DR     = SAVE(2)
L      = NINT(SAVE(3))
LBDCND = NINT(SAVE(4))
N      = NINT(SAVE(7))
NP     = NINT(SAVE(8))
DZSQR  = SAVE(9)
DZR    = SAVE(10)
ELMBDA = SAVE(11)
XMU    = SAVE(12)

!                               ALLOCATE SAVE ARRAY

IA = 13
IB = IA + L
IC = IB + L
ID = IC + L
IS = ID + L
ISVPS = IS + L

!                               ENTER BOUNDARY DATA FOR R-BOUNDARIES.

IF (LBDCND==1 .OR. LBDCND==2) THEN
  DO K = 1,N
    F(1,K) = F(1,K) - 2._EB*SAVE(IA+1)*BDRS(K)
  END DO
END IF

IF (LBDCND==3 .OR. LBDCND==4) THEN
  DO K = 1,N
    F(1,K) = F(1,K) + DR*SAVE(IA+1)*BDRS(K)
  END DO
END IF

IF (LBDCND==1 .OR. LBDCND==4 .OR. LBDCND==5) THEN
  DO K = 1,N
    F(L,K) = F(L,K) - 2._EB*SAVE(IC+L)*BDRF(K)
  END DO
END IF

IF (LBDCND==2 .OR. LBDCND==3 .OR. LBDCND==6) THEN
  DO K = 1,N
    F(L,K) = F(L,K) - DR*SAVE(IC+L)*BDRF(K)
  END DO
END IF

!                               ENTER BOUNDARY DATA FOR Z-BOUNDARIES.

IF (NP==2 .OR. NP==3) THEN
  DO I = 1,L
    F(I,1) = F(I,1) - 2._EB*DZSQR*BDZS(I)
  END DO
END IF

IF (NP==4 .OR. NP==5) THEN
  DO I = 1,L
    F(I,1) = F(I,1) + DZR*BDZS(I)
  END DO
END IF

IF (NP==2 .OR. NP==5) THEN
  DO I = 1,L
    F(I,N) = F(I,N) - 2._EB*DZSQR*BDZF(I)
  END DO
END IF

IF (NP==3 .OR. NP==4) THEN
  DO I = 1,L
    F(I,N) = F(I,N) - DZR*BDZF(I)
  END DO
END IF

PERTRB = 0._EB
ISING = 0

!                               FOR SINGULAR PROBLEMS ADJUST DATA TO
!                               INSURE A SOLUTION WILL EXIST.  GO THRU
!                               THIS CODE TWICE: ISING=1 FOR CALCULATING
!                               PERTRB; ISING=2 FOR NORMALIZING SOLUTION
!                               AFTER IT IS COMPUTED.

SELECT CASE(LBDCND)
CASE(1:2) ; GO TO 770
CASE(3)   ; GO TO 620
CASE(4:5) ; GO TO 770
CASE(6)   ; GO TO 620
END SELECT

620 CONTINUE
SELECT CASE(NP)
CASE(1)   ; GO TO 640
CASE(2:3) ; GO TO 770
CASE(4)   ; GO TO 640
CASE(5)   ; GO TO 770
END SELECT

640 CONTINUE
IF (ABS(ELMBDA)<=TWO_EPSILON_EB) THEN
  GO TO 650
ELSE
  GO TO 770
END IF

650 CONTINUE
IF (ABS(XMU)<=TWO_EPSILON_EB) THEN
  GO TO 665
ELSE
  GO TO 770
END IF
665 CONTINUE
ISING = 1
660 CONTINUE
W(1:L) = 0._EB
DO K = 1,N
  DO I = 1,L
    W(I) = W(I) + F(I,K)
  END DO
END DO

S3 = 0._EB
S1 = 0._EB
DO I = 1,L
  S3 = S3 + SAVE(IS+I)
  S1 = S1 + W(I)*SAVE(IS+I)
END DO
S3 = S3*N
PERT = S1/S3

DO K = 1,N
  DO I = 1,L
    F(I,K) = F(I,K) - PERT
  END DO
END DO

!                               IF NORMALIZING SOLUTION, RESTORE PERTRB
!                               AND JUMP TO END

IF (ISING==2) THEN
  PERTRB = PRTSAV
  GO TO 800
END IF

PRTSAV = PERT

770 CONTINUE

!                               SOLVE SYSTEM USING S2CFSS

CALL S2CFSS(LDIMF,F,SAVE(ISVPS+1),W)

!                               IF A SINGULAR PROBLEM,
!                               RE-NORMALIZE SOLUTION (ISING=2)

IF (ISING==1) THEN
  ISING = 2

  GO TO 660

END IF

800 CONTINUE
RETURN
END SUBROUTINE H2CYSS

END MODULE POIS
