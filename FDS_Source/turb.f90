! This module is useful for verification tests and development of
! turbulence models.

MODULE TURBULENCE

USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_POINTERS
USE MESH_VARIABLES
USE COMP_FUNCTIONS

IMPLICIT NONE

CHARACTER(255), PARAMETER :: turbid='$Id: turb.f90 2900 2008-12-18 18:26:05Z drjfloyd $'
CHARACTER(255), PARAMETER :: turbrev='$Revision: 2900 $'
CHARACTER(255), PARAMETER :: turbdate='$Date: 2008-12-18 13:26:05 -0500 (Thu, 18 Dec 2008) $'

PRIVATE
PUBLIC :: NS_ANALYTICAL_SOLUTION, INIT_TURB_ARRAYS, VARDEN_DYNSMAG, &
          GET_REV_turb, MEASURE_TURBULENCE_RESOLUTION, WERNER_WENGLE_WALL_MODEL, COMPRESSION_WAVE, &
          SURFACE_HEAT_FLUX_MODEL, SYNTHETIC_TURBULENCE, SYNTHETIC_EDDY_SETUP, TEST_FILTER
 
CONTAINS


SUBROUTINE INIT_TURB_ARRAYS(NM)
USE MEMORY_FUNCTIONS, ONLY: ChkMemErr
IMPLICIT NONE
INTEGER, INTENT(IN) :: NM
INTEGER :: IZERO
TYPE (MESH_TYPE), POINTER :: M

CALL POINT_TO_MESH(NM)
M => MESHES(NM)

IF (PERIODIC_TEST==2 .OR. DYNSMAG) THEN
   ! real work arrays
   ALLOCATE(M%TURB_WORK1(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
   CALL ChkMemErr('INIT_TURB_ARRAYS','TURB_WORK1',IZERO)
   M%TURB_WORK1 = 0._EB
   ALLOCATE(M%TURB_WORK2(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
   CALL ChkMemErr('INIT_TURB_ARRAYS','TURB_WORK2',IZERO)
   M%TURB_WORK2 = 0._EB
   ALLOCATE(M%TURB_WORK3(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
   CALL ChkMemErr('INIT_TURB_ARRAYS','TURB_WORK3',IZERO)
   M%TURB_WORK3 = 0._EB
   ALLOCATE(M%TURB_WORK4(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
   CALL ChkMemErr('INIT_TURB_ARRAYS','TURB_WORK4',IZERO)
   M%TURB_WORK4 = 0._EB
ENDIF

IF (DYNSMAG) THEN
   ! real work arrays
   ALLOCATE(M%TURB_WORK5(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
   CALL ChkMemErr('INIT_TURB_ARRAYS','TURB_WORK5',IZERO)
   M%TURB_WORK5 = 0._EB
   ALLOCATE(M%TURB_WORK6(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
   CALL ChkMemErr('INIT_TURB_ARRAYS','TURB_WORK6',IZERO)
   M%TURB_WORK6 = 0._EB
   ALLOCATE(M%TURB_WORK7(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
   CALL ChkMemErr('INIT_TURB_ARRAYS','TURB_WORK7',IZERO)
   M%TURB_WORK7 = 0._EB
   ALLOCATE(M%TURB_WORK8(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
   CALL ChkMemErr('INIT_TURB_ARRAYS','TURB_WORK8',IZERO)
   M%TURB_WORK8 = 0._EB
   
   ALLOCATE(M%TURB_WORK9(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
   CALL ChkMemErr('INIT_TURB_ARRAYS','TURB_WORK9',IZERO)
   M%TURB_WORK9 = 0._EB
   ALLOCATE(M%TURB_WORK10(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
   CALL ChkMemErr('INIT_TURB_ARRAYS','TURB_WORK10',IZERO)
   M%TURB_WORK10 = 0._EB
ENDIF
   
! 1D working arrays
IF (DYNSMAG .OR. CHECK_KINETIC_ENERGY) THEN
   ALLOCATE(M%TURB_WORK11(MAX(IBAR,JBAR,KBAR)),STAT=IZERO)
   CALL ChkMemErr('INIT_TURB_ARRAYS','TURB_WORK11',IZERO)
   M%TURB_WORK11 = 0._EB
   ALLOCATE(M%TURB_WORK12(MAX(IBAR,JBAR,KBAR)),STAT=IZERO)
   CALL ChkMemErr('INIT_TURB_ARRAYS','TURB_WORK12',IZERO)
   M%TURB_WORK12 = 0._EB
ENDIF

END SUBROUTINE INIT_TURB_ARRAYS


SUBROUTINE NS_ANALYTICAL_SOLUTION(NM)
IMPLICIT NONE
! Initialize flow variables with an analytical solution of the governing equations

INTEGER, INTENT(IN) :: NM
INTEGER :: I,J,K
REAL(EB) :: UU,WW

CALL POINT_TO_MESH(NM)

DO K=1,KBAR
   DO J=1,JBAR
      DO I=0,IBAR
         U(I,J,K) = 1._EB - 2._EB*COS(X(I))*SIN(ZC(K))
      ENDDO
   ENDDO
ENDDO
DO K=0,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         W(I,J,K) = 1._EB + 2._EB*SIN(XC(I))*COS(Z(K))
      ENDDO
   ENDDO
ENDDO
DO K=0,KBP1
   DO J=0,JBP1
      DO I=0,IBP1
         UU = 1._EB - 2._EB*COS(XC(I))*SIN(ZC(K))
         WW = 1._EB + 2._EB*SIN(XC(I))*COS(ZC(K))
         H(I,J,K) = -( COS(2._EB*XC(I)) + COS(2._EB*ZC(K)) ) + 0.5_EB*(UU**2+WW**2)
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE NS_ANALYTICAL_SOLUTION


SUBROUTINE COMPRESSION_WAVE(NM,T,ITEST)
IMPLICIT NONE

INTEGER, INTENT(IN) :: NM,ITEST
REAL(EB), INTENT(IN) :: T
INTEGER :: I,J,K

CALL POINT_TO_MESH(NM)

SELECT CASE(ITEST)
   CASE(3) ! stationary compression wave
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=0,IBAR
               U(I,J,K)  = 2._EB + SIN(X(I))
               US(I,J,K) = U(I,J,K)
            ENDDO
         ENDDO
      ENDDO
      DO K=1,KBAR
         DO J=0,JBAR
            DO I=1,IBAR
               V(I,J,K)  = 0._EB
               VS(I,J,K) = 0._EB
            ENDDO
         ENDDO
      ENDDO
      DO K=0,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               W(I,J,K)  = 3._EB + SIN(Z(K))
               WS(I,J,K) = W(I,J,K)
            ENDDO
         ENDDO
      ENDDO
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               D(I,J,K) = (U(I,J,K)-U(I-1,J,K))*RDX(I) + (W(I,J,K)-W(I,J,K-1))*RDZ(K)
               DS(I,J,K) = D(I,J,K)
            ENDDO
         ENDDO
      ENDDO
   CASE(4) ! pulsating dilation
      PREDICTOR_IF: IF (PREDICTOR) THEN
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=0,IBAR
                  US(I,J,K) = SIN(X(I))*COS(T)
               ENDDO
            ENDDO
         ENDDO
         DO K=1,KBAR
            DO J=0,JBAR
               DO I=1,IBAR
                  VS(I,J,K) = 0._EB
               ENDDO
            ENDDO
         ENDDO
         DO K=0,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  WS(I,J,K) = SIN(Z(K))*COS(T)
               ENDDO
            ENDDO
         ENDDO
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  DS(I,J,K) = (US(I,J,K)-US(I-1,J,K))*RDX(I) + (WS(I,J,K)-WS(I,J,K-1))*RDZ(K)
               ENDDO
            ENDDO
         ENDDO
      ELSE PREDICTOR_IF
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=0,IBAR
                  U(I,J,K) = SIN(X(I))*COS(T)
               ENDDO
            ENDDO
         ENDDO
         DO K=1,KBAR
            DO J=0,JBAR
               DO I=1,IBAR
                  V(I,J,K) = 0._EB
               ENDDO
            ENDDO
         ENDDO
         DO K=0,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  W(I,J,K) = SIN(Z(K))*COS(T)
               ENDDO
            ENDDO
         ENDDO
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  D(I,J,K) = (U(I,J,K)-U(I-1,J,K))*RDX(I) + (W(I,J,K)-W(I,J,K-1))*RDZ(K)
               ENDDO
            ENDDO
         ENDDO
      ENDIF PREDICTOR_IF
END SELECT

END SUBROUTINE COMPRESSION_WAVE


SUBROUTINE VARDEN_DYNSMAG(NM)
IMPLICIT NONE

!--------------------------------------------------------------
!     for all tensors, the indices are defined as follows...
!
!     |  11    12    13  |
!     |                  |
!     |  21    22    23  |
!     |                  |
!     |  31    32    33  |
!
!     I"name"I is the 'magnitude' of "name".
!--------------------------------------------------------------

INTEGER, INTENT(IN) :: NM

! Velocities relative to the p-cell center
REAL(EB) :: U_E,U_W,U_N,U_S,U_T,U_B,U_P2
REAL(EB) :: V_E,V_W,V_N,V_S,V_T,V_B,V_P2
REAL(EB) :: W_E,W_W,W_N,W_S,W_T,W_B,W_P2
REAL(EB) :: DELTA,SKK,TEMP_TERM
INTEGER :: I,J,K,N_LO(3),N_HI(3),ARRAY_LO(3),ARRAY_HI(3)

REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW,UP,VP,WP,RHOP,RHOPHAT
REAL(EB), POINTER, DIMENSION(:,:,:) :: S11,S22,S33,S12,S13,S23,SS
REAL(EB), POINTER, DIMENSION(:,:,:) :: SHAT11,SHAT22,SHAT33,SHAT12,SHAT13,SHAT23,SSHAT
REAL(EB), POINTER, DIMENSION(:,:,:) :: BETA11,BETA22,BETA33,BETA12,BETA13,BETA23
REAL(EB), POINTER, DIMENSION(:,:,:) :: BETAHAT11,BETAHAT22,BETAHAT33,BETAHAT12,BETAHAT13,BETAHAT23
REAL(EB), POINTER, DIMENSION(:,:,:) :: M11,M22,M33,M12,M13,M23,MM,MMHAT
REAL(EB), POINTER, DIMENSION(:,:,:) :: L11,L22,L33,L12,L13,L23,ML,MLHAT

REAL(EB), PARAMETER :: ALPHA = 6.0_EB ! See Lund, 1997 CTR briefs.

! References:
!
! M. Germano, U. Piomelli, P. Moin, and W. Cabot.  A dynamic subgrid-scale eddy viscosity model.
! Phys. Fluids A, 3(7):1760-1765, 1991.
!
! M. Pino Martin, U. Piomelli, and G. Candler. Subgrid-scale models for compressible large-eddy
! simulation. Theoret. Comput. Fluid Dynamics, 13:361-376, 2000.
!
! P. Moin, K. Squires, W. Cabot, and S. Lee.  A dynamic subgrid-scale model for compressible
! turbulence and scalar transport. Phys. Fluids A, 3(11):2746-2757, 1991.
!
! T. S. Lund. On the use of discrete filters for large eddy simulation.  Center for Turbulence
! Research Annual Research Briefs, 1997.
!
! R. McDermott. Variable density formulation of the dynamic Smagorinsky model.
! http://randy.mcdermott.googlepages.com/dynsmag_comp.pdf

! *****************************************************************************
! CAUTION WHEN MODIFYING: The order in which the tensor components are computed
! is important because we overwrite pointers several times to conserve memory.
! *****************************************************************************

CALL POINT_TO_MESH(NM)

N_LO = 1
N_HI = (/IBAR,JBAR,KBAR/)

ARRAY_LO = 0
ARRAY_HI = (/IBP1,JBP1,KBP1/)

RHOP=>WORK8 
IF (PREDICTOR) THEN
   UU=>U
   VV=>V
   WW=>W
   RHOP=RHO(0:IBP1,0:JBP1,0:KBP1)
ELSE
   UU=>US
   VV=>VS
   WW=>WS
   RHOP=RHOS(0:IBP1,0:JBP1,0:KBP1)
ENDIF

UP => TURB_WORK1
VP => TURB_WORK2
WP => TURB_WORK3

S11 => WORK1
S22 => WORK2
S33 => WORK3
S12 => WORK4
S13 => WORK5
S23 => WORK6
SS  => WORK7


!$OMP PARALLEL DO COLLAPSE(3) PRIVATE(K,J,I,U_E,U_W,U_N,U_S,U_T,U_B,V_E,V_W,V_N,V_S,V_T,V_B,W_E,W_W,W_N,W_S,W_T,W_B,SKK) &
!$OMP PRIVATE(U_P2,V_P2,W_P2)
DO K = N_LO(3),N_HI(3)
   DO J = N_LO(2),N_HI(2)
      DO I = N_LO(1),N_HI(1)

         U_E = UU(I,J,K)
         U_W = UU(I-1,J,K)
         UP(I,J,K) = 0.5_EB*(U_E + U_W)
         U_P2 = 0.5_EB*UP(I,J,K)
         U_N = U_P2+0.25_EB*(UU(I,J+1,K) + UU(I-1,J+1,K) )
         U_S = U_P2+0.25_EB*(UU(I,J-1,K) + UU(I-1,J-1,K) )
         U_T = U_P2+0.25_EB*(UU(I,J,K+1) + UU(I-1,J,K+1) )
         U_B = U_P2+0.25_EB*(UU(I,J,K-1) + UU(I-1,J,K-1) )

         V_N = VV(I,J,K)
         V_S = VV(I,J-1,K)
         VP(I,J,K) = 0.5_EB*(V_N + V_S)
         V_P2 = 0.5_EB*VP(I,J,K)         
         V_E = V_P2+0.25_EB*(VV(I+1,J,K) + VV(I+1,J-1,K) )
         V_W = V_P2+0.25_EB*(VV(I-1,J,K) + VV(I-1,J-1,K) )
         V_T = V_P2+0.25_EB*(VV(I,J,K+1) + VV(I,J-1,K+1) )
         V_B = V_P2+0.25_EB*(VV(I,J,K-1) + VV(I,J-1,K-1) )

         W_T = WW(I,J,K)
         W_B = WW(I,J,K-1)
         WP(I,J,K) = 0.5_EB*(W_T + W_B)
         W_P2 = 0.5_EB*WP(I,J,K)
         W_E = W_P2+0.25_EB*(WW(I+1,J,K) + WW(I+1,J,K-1) ) 
         W_W = W_P2+0.25_EB*(WW(I-1,J,K) + WW(I-1,J,K-1) )
         W_N = W_P2+0.25_EB*(WW(I,J+1,K) + WW(I,J+1,K-1) )
         W_S = W_P2+0.25_EB*(WW(I,J-1,K) + WW(I,J-1,K-1) )

         ! calculate the grid strain rate tensor
         
         S11(I,J,K) = (U_E - U_W)/DX(I)
         S22(I,J,K) = (V_N - V_S)/DY(J)
         S33(I,J,K) = (W_T - W_B)/DZ(K)
         SKK = S11(I,J,K) + S22(I,J,K) + S33(I,J,K)
         S11(I,J,K) = S11(I,J,K) - ONTH*SKK
         S22(I,J,K) = S22(I,J,K) - ONTH*SKK
         S33(I,J,K) = S33(I,J,K) - ONTH*SKK
         S12(I,J,K) = 0.5_EB*( (U_N - U_S)/DY(J) + (V_E - V_W)/DX(I) )
         S13(I,J,K) = 0.5_EB*( (U_T - U_B)/DZ(K) + (W_E - W_W)/DX(I) )
         S23(I,J,K) = 0.5_EB*( (V_T - V_B)/DZ(K) + (W_N - W_S)/DY(J) )
         
         ! calculate magnitude of the grid strain rate

         SS(I,J,K) = SQRT(2._EB*(S11(I,J,K)*S11(I,J,K) + &
                                 S22(I,J,K)*S22(I,J,K) + &
                                 S33(I,J,K)*S33(I,J,K) + &
                          2._EB*(S12(I,J,K)*S12(I,J,K) + &
                                 S13(I,J,K)*S13(I,J,K) + &
                                 S23(I,J,K)*S23(I,J,K)) ) )

      ENDDO
   ENDDO
ENDDO
!$OMP END PARALLEL DO

! test filter the strain rate

SHAT11 => TURB_WORK4
SHAT22 => TURB_WORK5
SHAT33 => TURB_WORK6
SHAT12 => TURB_WORK7
SHAT13 => TURB_WORK8
SHAT23 => TURB_WORK9
SSHAT  => TURB_WORK10

CALL TEST_FILTER(SHAT11,S11,N_LO,N_HI,ARRAY_LO,ARRAY_HI)
CALL TEST_FILTER(SHAT22,S22,N_LO,N_HI,ARRAY_LO,ARRAY_HI)
CALL TEST_FILTER(SHAT33,S33,N_LO,N_HI,ARRAY_LO,ARRAY_HI)
CALL TEST_FILTER(SHAT12,S12,N_LO,N_HI,ARRAY_LO,ARRAY_HI)
CALL TEST_FILTER(SHAT13,S13,N_LO,N_HI,ARRAY_LO,ARRAY_HI)
CALL TEST_FILTER(SHAT23,S23,N_LO,N_HI,ARRAY_LO,ARRAY_HI)

! calculate magnitude of test filtered strain rate

!$OMP PARALLEL DO COLLAPSE(3) PRIVATE(K,J,I)
DO K = N_LO(3),N_HI(3)
   DO J = N_LO(2),N_HI(2)
      DO I = N_LO(1),N_HI(1)
      
         SSHAT(I,J,K) = SQRT(2._EB*(SHAT11(I,J,K)*SHAT11(I,J,K) + &
                                    SHAT22(I,J,K)*SHAT22(I,J,K) + &
                                    SHAT33(I,J,K)*SHAT33(I,J,K) + &
                             2._EB*(SHAT12(I,J,K)*SHAT12(I,J,K) + &
                                    SHAT13(I,J,K)*SHAT13(I,J,K) + &
                                    SHAT23(I,J,K)*SHAT23(I,J,K)) ) )

      ENDDO
   ENDDO
ENDDO
!$OMP END PARALLEL DO

! calculate the grid filtered stress tensor, beta

BETA11 => WORK1
BETA22 => WORK2
BETA33 => WORK3
BETA12 => WORK4
BETA13 => WORK5
BETA23 => WORK6

!$OMP PARALLEL WORKSHARE
BETA11 = RHOP*SS*S11
BETA22 = RHOP*SS*S22
BETA33 = RHOP*SS*S33
BETA12 = RHOP*SS*S12
BETA13 = RHOP*SS*S13
BETA23 = RHOP*SS*S23
!$OMP END PARALLEL WORKSHARE

! test filter the grid filtered stress tensor

BETAHAT11 => WORK1
BETAHAT22 => WORK2
BETAHAT33 => WORK3
BETAHAT12 => WORK4
BETAHAT13 => WORK5
BETAHAT23 => WORK6

CALL TEST_FILTER(BETAHAT11,BETA11,N_LO,N_HI,ARRAY_LO,ARRAY_HI)
CALL TEST_FILTER(BETAHAT22,BETA22,N_LO,N_HI,ARRAY_LO,ARRAY_HI)
CALL TEST_FILTER(BETAHAT33,BETA33,N_LO,N_HI,ARRAY_LO,ARRAY_HI)
CALL TEST_FILTER(BETAHAT12,BETA12,N_LO,N_HI,ARRAY_LO,ARRAY_HI)
CALL TEST_FILTER(BETAHAT13,BETA13,N_LO,N_HI,ARRAY_LO,ARRAY_HI)
CALL TEST_FILTER(BETAHAT23,BETA23,N_LO,N_HI,ARRAY_LO,ARRAY_HI)

! test filter the density

RHOPHAT => WORK7
CALL TEST_FILTER(RHOPHAT,RHOP,N_LO,N_HI,ARRAY_LO,ARRAY_HI)
! calculate the Mij tensor

M11 => WORK1
M22 => WORK2
M33 => WORK3
M12 => WORK4
M13 => WORK5
M23 => WORK6

!$OMP PARALLEL DO COLLAPSE(3) PRIVATE(K,J,I,TEMP_TERM)
DO K = N_LO(3),N_HI(3)
   DO J = N_LO(2),N_HI(2)
      DO I = N_LO(1),N_HI(1)
         TEMP_TERM = ALPHA*RHOPHAT(I,J,K)*SSHAT(I,J,K)
         M11(I,J,K) = 2._EB*(BETAHAT11(I,J,K) - TEMP_TERM*SHAT11(I,J,K))
         M22(I,J,K) = 2._EB*(BETAHAT22(I,J,K) - TEMP_TERM*SHAT22(I,J,K))
         M33(I,J,K) = 2._EB*(BETAHAT33(I,J,K) - TEMP_TERM*SHAT33(I,J,K))
         M12(I,J,K) = 2._EB*(BETAHAT12(I,J,K) - TEMP_TERM*SHAT12(I,J,K))
         M13(I,J,K) = 2._EB*(BETAHAT13(I,J,K) - TEMP_TERM*SHAT13(I,J,K))
         M23(I,J,K) = 2._EB*(BETAHAT23(I,J,K) - TEMP_TERM*SHAT23(I,J,K))
         
      ENDDO
   ENDDO
ENDDO
!$OMP END PARALLEL DO

! calculate the Leonard term, Lij

L11 => TURB_WORK4
L22 => TURB_WORK5
L33 => TURB_WORK6
L12 => TURB_WORK7
L13 => TURB_WORK8
L23 => TURB_WORK9

CALL CALC_VARDEN_LEONARD_TERM(NM)

! calculate Mij*Lij & Mij*Mij

MM    => TURB_WORK1
MMHAT => TURB_WORK1

ML    => TURB_WORK2
MLHAT => TURB_WORK2

!$OMP PARALLEL DO COLLAPSE(3) PRIVATE(K,J,I)
DO K = N_LO(3),N_HI(3)
   DO J = N_LO(2),N_HI(2)
      DO I = N_LO(1),N_HI(1)
      
         ML(I,J,K) = M11(I,J,K)*L11(I,J,K) + M22(I,J,K)*L22(I,J,K) + M33(I,J,K)*L33(I,J,K) + &
              2._EB*(M12(I,J,K)*L12(I,J,K) + M13(I,J,K)*L13(I,J,K) + M23(I,J,K)*L23(I,J,K))
       
         MM(I,J,K) = M11(I,J,K)*M11(I,J,K) + M22(I,J,K)*M22(I,J,K) + M33(I,J,K)*M33(I,J,K) + &
              2._EB*(M12(I,J,K)*M12(I,J,K) + M13(I,J,K)*M13(I,J,K) + M23(I,J,K)*M23(I,J,K))
              
      ENDDO
   ENDDO
ENDDO
!$OMP END PARALLEL DO

! do some smoothing

CALL TEST_FILTER(MLHAT,ML,N_LO,N_HI,ARRAY_LO,ARRAY_HI)
CALL TEST_FILTER(MMHAT,MM,N_LO,N_HI,ARRAY_LO,ARRAY_HI)

!$OMP PARALLEL DO COLLAPSE(3) PRIVATE(K,J,I,DELTA)
DO K = N_LO(3),N_HI(3)
   DO J = N_LO(2),N_HI(2)
      DO I = N_LO(1),N_HI(1)

         ! calculate the local Smagorinsky coefficient

         ! perform "clipping" in case MLij is negative...
         IF (MLHAT(I,J,K) < 0._EB) MLHAT(I,J,K) = 0._EB

         ! calculate the effective viscosity

         ! handle the case where we divide by zero, note MMHAT is positive semi-definite
         IF (MMHAT(I,J,K) == 0._EB) THEN
            C_DYNSMAG(I,J,K) = 0._EB
         ELSE
            ! filter width
            IF (TWO_D) THEN
               DELTA = (DX(I)*DZ(K))**0.5_EB
            ELSE
               DELTA = (DX(I)*DY(J)*DZ(K))**ONTH
            ENDIF
            IF (USE_MAX_FILTER_WIDTH) DELTA=MAX(DX(I),DY(J),DZ(K))
            C_DYNSMAG(I,J,K) = SQRT(MLHAT(I,J,K)/MMHAT(I,J,K))/DELTA
         ENDIF
         
         ! clip max value of CS, note that CS*DELTA is the "mixing length", so DELTA
         ! is a reasonable upper bound
         C_DYNSMAG(I,J,K) = MIN(C_DYNSMAG(I,J,K),1._EB)
         
      END DO
   END DO
END DO
!$OMP END PARALLEL DO

END SUBROUTINE VARDEN_DYNSMAG


SUBROUTINE CALC_VARDEN_LEONARD_TERM(NM)
IMPLICIT NONE

INTEGER, INTENT(IN) :: NM

REAL(EB), POINTER, DIMENSION(:,:,:) :: L11,L22,L33,L12,L13,L23
REAL(EB), POINTER, DIMENSION(:,:,:) :: RHOP,RHOPHAT
REAL(EB), POINTER, DIMENSION(:,:,:) :: UP,VP,WP
REAL(EB), POINTER, DIMENSION(:,:,:) :: RUU,RVV,RWW,RUV,RUW,RVW
REAL(EB), POINTER, DIMENSION(:,:,:) :: RU,RV,RW
REAL(EB), POINTER, DIMENSION(:,:,:) :: RUU_HAT,RVV_HAT,RWW_HAT,RUV_HAT,RUW_HAT,RVW_HAT
REAL(EB), POINTER, DIMENSION(:,:,:) :: RU_HAT,RV_HAT,RW_HAT
REAL(EB) :: INV_RHOPHAT
TYPE(MESH_TYPE), POINTER :: M
INTEGER :: I,J,K,N_LO(3),N_HI(3),ARRAY_LO(3),ARRAY_HI(3)

! *****************************************************************************
! CAUTION WHEN MODIFYING: The order in which the tensor components are computed
! is important because we overwrite pointers several times to conserve memory.
! *****************************************************************************

M => MESHES(NM)

IF (PREDICTOR) THEN
   RHOP=>M%RHO
ELSE
   RHOP=>M%RHOS
ENDIF
RHOPHAT => M%WORK7

N_LO = 1
N_HI = (/M%IBAR,M%JBAR,M%KBAR/)

ARRAY_LO = 0
ARRAY_HI = (/M%IBP1,M%JBP1,M%KBP1/)

! Compute rho*UiUj

UP => M%TURB_WORK1 ! will be overwritten by RU
VP => M%TURB_WORK2
WP => M%TURB_WORK3

RUU => M%TURB_WORK4 ! will be overwritten by RUU_HAT
RVV => M%TURB_WORK5
RWW => M%TURB_WORK6
RUV => M%TURB_WORK7
RUW => M%TURB_WORK8
RVW => M%TURB_WORK9

!$OMP PARALLEL DO COLLAPSE(3) PRIVATE(K,J,I)
DO K=N_LO(3),N_HI(3)
   DO J=N_LO(2),N_HI(2)
      DO I=N_LO(1),N_HI(1)

         RUU(I,J,K) = RHOP(I,J,K)*UP(I,J,K)*UP(I,J,K)
         RVV(I,J,K) = RHOP(I,J,K)*VP(I,J,K)*VP(I,J,K)
         RWW(I,J,K) = RHOP(I,J,K)*WP(I,J,K)*WP(I,J,K)
         RUV(I,J,K) = RHOP(I,J,K)*UP(I,J,K)*VP(I,J,K)
         RUW(I,J,K) = RHOP(I,J,K)*UP(I,J,K)*WP(I,J,K)
         RVW(I,J,K) = RHOP(I,J,K)*VP(I,J,K)*WP(I,J,K)

      ENDDO
   ENDDO
ENDDO
!$OMP END PARALLEL DO

! Test filter rho*UiUj

RUU_HAT => M%TURB_WORK4 ! will be overwritten by Lij
RVV_HAT => M%TURB_WORK5
RWW_HAT => M%TURB_WORK6
RUV_HAT => M%TURB_WORK7
RUW_HAT => M%TURB_WORK8
RVW_HAT => M%TURB_WORK9

CALL TEST_FILTER(RUU_HAT,RUU,N_LO,N_HI,ARRAY_LO,ARRAY_HI)
CALL TEST_FILTER(RVV_HAT,RVV,N_LO,N_HI,ARRAY_LO,ARRAY_HI)
CALL TEST_FILTER(RWW_HAT,RWW,N_LO,N_HI,ARRAY_LO,ARRAY_HI)
CALL TEST_FILTER(RUV_HAT,RUV,N_LO,N_HI,ARRAY_LO,ARRAY_HI)
CALL TEST_FILTER(RUW_HAT,RUW,N_LO,N_HI,ARRAY_LO,ARRAY_HI)
CALL TEST_FILTER(RVW_HAT,RVW,N_LO,N_HI,ARRAY_LO,ARRAY_HI)

! Compute rho*Ui

RU => M%TURB_WORK1 ! will be overwritten by RU_HAT
RV => M%TURB_WORK2
RW => M%TURB_WORK3

!$OMP PARALLEL WORKSHARE
RU = RHOP(0:IBP1,0:JBP1,0:KBP1)*UP
RV = RHOP(0:IBP1,0:JBP1,0:KBP1)*VP
RW = RHOP(0:IBP1,0:JBP1,0:KBP1)*WP
!$OMP END PARALLEL WORKSHARE

! Test filter rho*Ui

RU_HAT => M%TURB_WORK1
RV_HAT => M%TURB_WORK2
RW_HAT => M%TURB_WORK3

CALL TEST_FILTER(RU_HAT,RU,N_LO,N_HI,ARRAY_LO,ARRAY_HI)
CALL TEST_FILTER(RV_HAT,RV,N_LO,N_HI,ARRAY_LO,ARRAY_HI)
CALL TEST_FILTER(RW_HAT,RW,N_LO,N_HI,ARRAY_LO,ARRAY_HI)

! Compute variable density Leonard stress

L11 => M%TURB_WORK4
L22 => M%TURB_WORK5
L33 => M%TURB_WORK6
L12 => M%TURB_WORK7
L13 => M%TURB_WORK8
L23 => M%TURB_WORK9

!$OMP PARALLEL DO COLLAPSE(3) PRIVATE(K,J,I,INV_RHOPHAT)
DO K=N_LO(3),N_HI(3)
   DO J=N_LO(2),N_HI(2)
      DO I=N_LO(1),N_HI(1)
         IF (RHOPHAT(I,J,K)>0._EB) THEN
            INV_RHOPHAT = 1._EB/RHOPHAT(I,J,K)
            L11(I,J,K) = RUU_HAT(I,J,K) - RU_HAT(I,J,K)*RU_HAT(I,J,K)*INV_RHOPHAT 
            L22(I,J,K) = RVV_HAT(I,J,K) - RV_HAT(I,J,K)*RV_HAT(I,J,K)*INV_RHOPHAT 
            L33(I,J,K) = RWW_HAT(I,J,K) - RW_HAT(I,J,K)*RW_HAT(I,J,K)*INV_RHOPHAT 
            L12(I,J,K) = RUV_HAT(I,J,K) - RU_HAT(I,J,K)*RV_HAT(I,J,K)*INV_RHOPHAT 
            L13(I,J,K) = RUW_HAT(I,J,K) - RU_HAT(I,J,K)*RW_HAT(I,J,K)*INV_RHOPHAT 
            L23(I,J,K) = RVW_HAT(I,J,K) - RV_HAT(I,J,K)*RW_HAT(I,J,K)*INV_RHOPHAT 
         ELSE
            L11(I,J,K) = 0._EB
            L22(I,J,K) = 0._EB
            L33(I,J,K) = 0._EB
            L12(I,J,K) = 0._EB
            L13(I,J,K) = 0._EB
            L23(I,J,K) = 0._EB
         ENDIF
      ENDDO
   ENDDO
ENDDO
!$OMP END PARALLEL DO

END SUBROUTINE CALC_VARDEN_LEONARD_TERM


SUBROUTINE TEST_FILTER(PHIBAR,PHI,N_LO,N_HI,ARRAY_LO,ARRAY_HI)
IMPLICIT NONE

INTEGER, INTENT(IN) :: N_LO(3),N_HI(3),ARRAY_LO(3),ARRAY_HI(3)
REAL(EB), INTENT(IN) :: PHI(ARRAY_LO(1):ARRAY_HI(1),ARRAY_LO(2):ARRAY_HI(2),ARRAY_LO(3):ARRAY_HI(3))
REAL(EB), INTENT(OUT) :: PHIBAR(ARRAY_LO(1):ARRAY_HI(1),ARRAY_LO(2):ARRAY_HI(2),ARRAY_LO(3):ARRAY_HI(3))
REAL(EB), POINTER, DIMENSION(:) :: PHI1,PHI2
INTEGER I,J,K

PHI1 => TURB_WORK11
PHI2 => TURB_WORK12

PHIBAR = PHI

! filter in x:
DO K = N_LO(3),N_HI(3)
   DO J = N_LO(2),N_HI(2)
      PHI1(N_LO(1):N_HI(1)) = PHIBAR(N_LO(1):N_HI(1),J,K)
      CALL TOPHAT_FILTER_1D(PHI2(N_LO(1):N_HI(1)),PHI1(N_LO(1):N_HI(1)),N_LO(1),N_HI(1))
      PHIBAR(N_LO(1):N_HI(1),J,K) = PHI2(N_LO(1):N_HI(1))
   ENDDO
ENDDO

IF (.NOT.TWO_D) THEN
   ! filter in y:
   DO K = N_LO(3),N_HI(3)
      DO I = N_LO(1),N_HI(1)
         PHI1(N_LO(2):N_HI(2)) = PHIBAR(I,N_LO(2):N_HI(2),K)
         CALL TOPHAT_FILTER_1D(PHI2(N_LO(2):N_HI(2)),PHI1(N_LO(2):N_HI(2)),N_LO(2),N_HI(2))
         PHIBAR(I,N_LO(2):N_HI(2),K) = PHI2(N_LO(2):N_HI(2))
      ENDDO
   ENDDO
ENDIF

! filter in z:
DO J = N_LO(2),N_HI(2)
   DO I = N_LO(1),N_HI(1)
      PHI1(N_LO(3):N_HI(3)) = PHIBAR(I,J,N_LO(3):N_HI(3))
      CALL TOPHAT_FILTER_1D(PHI2(N_LO(3):N_HI(3)),PHI1(N_LO(3):N_HI(3)),N_LO(3),N_HI(3))
      PHIBAR(I,J,N_LO(3):N_HI(3)) = PHI2(N_LO(3):N_HI(3))
   ENDDO
ENDDO

END SUBROUTINE TEST_FILTER


SUBROUTINE TOPHAT_FILTER_1D(UBAR,U,N_LO,N_HI)
IMPLICIT NONE

INTEGER, INTENT(IN) :: N_LO,N_HI
REAL(EB), INTENT(IN) :: U(N_LO:N_HI)
REAL(EB), INTENT(OUT) :: UBAR(N_LO:N_HI)
INTEGER :: J
!REAL(EB), POINTER, DIMENSION(:) :: UU
REAL(EB),PARAMETER:: W(-1:1) = (/0.25_EB,0.5_EB,0.25_EB/)   ! trapezoid rule
!REAL(EB),PARAMETER::W(-1:1) = (/ONSI,TWTH,ONSI/)           ! Simpson's rule


!UU => WORK
!UU(N_LO:N_HI) = U

! Filter the u field to obtain ubar
DO J=N_LO+1,N_HI-1
   UBAR(J) = DOT_PRODUCT(W(-1:1),U(J-1:J+1))
ENDDO
! set boundary values (not ideal, but fast and simple)
UBAR(N_LO) = UBAR(N_LO+1)
UBAR(N_HI) = UBAR(N_HI-1)

END SUBROUTINE TOPHAT_FILTER_1D


SUBROUTINE MEASURE_TURBULENCE_RESOLUTION(NM)
IMPLICIT NONE

INTEGER, INTENT(IN) :: NM
REAL(EB) :: KSGS
INTEGER :: I,J,K,N_LO(3),N_HI(3),ARRAY_LO(3),ARRAY_HI(3)
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW,UP,VP,WP,UP_HAT,VP_HAT,WP_HAT

CALL POINT_TO_MESH(NM)

N_LO = 1
N_HI = (/IBAR,JBAR,KBAR/)

ARRAY_LO = 0
ARRAY_HI = (/IBP1,JBP1,KBP1/)

IF (PREDICTOR) THEN
   UU=>US
   VV=>VS
   WW=>WS
ELSE
   UU=>U
   VV=>V
   WW=>W
ENDIF

! Velocities relative to the p-cell center
UP => WORK1
VP => WORK2
WP => WORK3

DO K = N_LO(3),N_HI(3)
   DO J = N_LO(2),N_HI(2)
      DO I = N_LO(1),N_HI(1)
         UP(I,J,K) = 0.5_EB*(UU(I,J,K) + UU(I-1,J,K))
         VP(I,J,K) = 0.5_EB*(VV(I,J,K) + VV(I,J-1,K))
         WP(I,J,K) = 0.5_EB*(WW(I,J,K) + WW(I,J,K-1))
      ENDDO
   ENDDO
ENDDO

UP_HAT => WORK4
VP_HAT => WORK5
WP_HAT => WORK6

CALL TEST_FILTER(UP_HAT,UP,N_LO,N_HI,ARRAY_LO,ARRAY_HI)
CALL TEST_FILTER(VP_HAT,VP,N_LO,N_HI,ARRAY_LO,ARRAY_HI)
CALL TEST_FILTER(WP_HAT,WP,N_LO,N_HI,ARRAY_LO,ARRAY_HI)

DO K = N_LO(3),N_HI(3)
   DO J = N_LO(2),N_HI(2)
      DO I = N_LO(1),N_HI(1)

         IF (KRES(I,J,K)>TKE_TOLERANCE) THEN
            KSGS = 0.5_EB*( (UP(I,J,K)-UP_HAT(I,J,K))**2 + (VP(I,J,K)-VP_HAT(I,J,K))**2 + (WP(I,J,K)-WP_HAT(I,J,K))**2 )
            MTR(I,J,K) = KSGS/(KRES(I,J,K)+KSGS)
         ELSE
            MTR(I,J,K) = 0._EB
         ENDIF
         
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE MEASURE_TURBULENCE_RESOLUTION


SUBROUTINE WERNER_WENGLE_WALL_MODEL(SF,U_TAU,U1,NU,DZ,ROUGHNESS)
IMPLICIT NONE

REAL(EB), INTENT(OUT) :: SF
REAL(EB), INTENT(IN) :: U1,NU,DZ,ROUGHNESS

REAL(EB), PARAMETER :: A=8.3_EB,B=1._EB/7._EB
REAL(EB), PARAMETER :: Z_PLUS_TURBULENT = 11.81_EB
REAL(EB), PARAMETER :: ALPHA=7.202125273562269_EB !! ALPHA=(1._EB-B)/2._EB*A**((1._EB+B)/(1._EB-B))
REAL(EB), PARAMETER :: BETA=1._EB+B
REAL(EB), PARAMETER :: ETA=(1._EB+B)/A
REAL(EB), PARAMETER :: GAMMA=2._EB/(1._EB+B)
REAL(EB), PARAMETER :: RKAPPA=2.44_EB ! 1./von Karman constant
REAL(EB), PARAMETER :: BTILDE=7.44_EB ! see Pope p. 297 (constant has been modified)

REAL(EB) :: U_TAU,TAU_W,NU_OVER_DZ,Z_PLUS,TAU_ROUGH

! References (for smooth walls):
!
! Werner, H., Wengle, H. (1991) Large-eddy simulation of turbulent flow over
! and around a cube in a plate channel. 8th Symposium on Turbulent Shear
! Flows, Munich, Germany.
!
! Pierre Sagaut. Large Eddy Simulation for Incompressible Flows: An Introduction.
! Springer, 2001.
!
! Temmerman, L., Leschziner, M.A., Mellen, C.P., and Frohlich, J. (2003)
! Investigation of wall-function approximations and subgrid-scale models in
! Large Eddy Simulation of separated flow in a channel with streamwise
! periodic constrictions. International Journal of Heat and Fluid Flow,
! Vol. 24, No. 2, pp. 157-180.
!
! Breuer, M., Kniazev, B., and Abel, M. (2007) Development of wall models
! for LES of separated flows using statistical evaluations. Computers and
! Fluids, Vol. 36, pp. 817-837.
!
! McDermott, R. (2009) FDS Wall Flows, Part I: Straight Channels, NIST Technical Note.
!
! References (for rough surfaces):
!
! S. B. Pope (2000) Turbulent Flows, Cambridge.
!
! Moeng, C.-H. (1984) A Large-Eddy Simulation Model for the Study of Planetary
! Boundary-Layer Turbulence. Journal of the Atmospheric Sciences, Vol. 41, No. 13,
! pp. 2052-2062.
!
! Stoll, R., Porte-Agel, F. (2008) Large-Eddy Simulation of the Stable Atmospheric
! Boundary Layer using Dynamic Models with Different Averaging Schemes. Boundary-Layer
! Meteorology, 126:1-28.
!
! Comments:
!
! The slip factor (SF) is based on the following approximation to the wall stress
! (note that u0 is the ghost cell value of the streamwise velocity component and
! z is the wall-normal direction):
! tau_w = mu*(u1-u0)/dz = mu*(u1-SF*u1)/dz = mu*u1/dz*(1-SF)
! note that tau_w/rho = nu*u1/dz*(1-SF)

TAU_ROUGH = 0._EB
IF (ROUGHNESS>0._EB) THEN
   ! Pope (2000)
   TAU_ROUGH = ( U1/(RKAPPA*LOG(0.5_EB*DZ/ROUGHNESS)+BTILDE) )**2 ! actually tau_w/rho
ENDIF
! Werner-Wengle
NU_OVER_DZ = NU/DZ
TAU_W = (ALPHA*(NU_OVER_DZ)**BETA + ETA*(NU_OVER_DZ)**B*ABS(U1))**GAMMA ! actually tau_w/rho
TAU_W = MAX(TAU_W,TAU_ROUGH)
U_TAU = SQRT(TAU_W)
Z_PLUS = DZ/(NU/(U_TAU+1.E-10_EB))
IF (Z_PLUS>Z_PLUS_TURBULENT) THEN
   SF = 1._EB-TAU_W/(NU/DZ*ABS(U1)) ! log layer
ELSE
   SF = -1._EB ! viscous sublayer
ENDIF

!! check values...
!IF (Z_PLUS>Z_PLUS_TURBULENT) THEN
!   print *,'A = ',A
!   print *,'B = ',B
!   print *,'ALPHA = ',ALPHA
!   print *,'BETA = ',BETA
!   print *,'ETA = ',ETA
!   print *,'GAMMA = ',GAMMA
!   print *,'U1 = ',U1
!   print *,'NU/DZ = ',NU_OVER_DZ
!   print *,'TAU_W/RHO = ',TAU_W
!   print *,'Z_PLUS = ',Z_PLUS
!   print *,'SF = ',SF
!   print *
!ENDIF

END SUBROUTINE WERNER_WENGLE_WALL_MODEL


SUBROUTINE SURFACE_HEAT_FLUX_MODEL(H,U_TAU,DZ,ROUGHNESS,IOR,RHO,CP)

REAL(EB), INTENT(OUT) :: H ! heat transfer coefficient
REAL(EB), INTENT(IN) :: U_TAU,DZ,ROUGHNESS,RHO,CP
INTEGER, INTENT(IN) :: IOR
REAL(EB), PARAMETER :: KAPPA=0.41_EB ! von Karman constant
REAL(EB) :: PSI,MOL,Z0

! References:
!
! Stoll, R., Porte-Agel, F. (2008) Large-Eddy Simulation of the Stable Atmospheric
! Boundary Layer using Dynamic Models with Different Averaging Schemes. Boundary-Layer
! Meteorology, 126:1-28.

PSI = 0._EB
MOL = 0._EB
Z0 = MAX(ROUGHNESS,1.E-6_EB)

! atmospheric stability correction (use later)
IF (IOR==3) THEN
   MOL = 0._EB !! -U_TAU**3*THETA/(KAPPA*GRAV*HEAT_FLUX)
   PSI = 0._EB !! -7.8_EB*0.5*DZ/MOL
ENDIF

H = RHO*CP*U_TAU*KAPPA/(LOG(0.5_EB*DZ/Z0)-PSI)

END SUBROUTINE SURFACE_HEAT_FLUX_MODEL


SUBROUTINE SYNTHETIC_EDDY_SETUP(NM)
IMPLICIT NONE

INTEGER, INTENT(IN) :: NM
TYPE(VENTS_TYPE), POINTER :: VT=>NULL()
INTEGER :: NE,NV,IERROR
REAL(EB), POINTER, DIMENSION(:,:) :: A_IJ=>NULL(),R_IJ=>NULL()

VENT_LOOP: DO NV=1,MESHES(NM)%N_VENT
   VT => MESHES(NM)%VENTS(NV)
   IF (VT%N_EDDY==0) CYCLE VENT_LOOP
   
   VT%X_EDDY_MIN = VT%X1-VT%L_EDDY
   VT%X_EDDY_MAX = VT%X2+VT%L_EDDY
   VT%Y_EDDY_MIN = VT%Y1-VT%L_EDDY
   VT%Y_EDDY_MAX = VT%Y2+VT%L_EDDY
   VT%Z_EDDY_MIN = VT%Z1-VT%L_EDDY
   VT%Z_EDDY_MAX = VT%Z2+VT%L_EDDY
   VT%EDDY_BOX_VOLUME = (VT%X_EDDY_MAX-VT%X_EDDY_MIN)*(VT%Y_EDDY_MAX-VT%Y_EDDY_MIN)*(VT%Z_EDDY_MAX-VT%Z_EDDY_MIN)
   
   EDDY_LOOP: DO NE=1,VT%N_EDDY
      CALL EDDY_POSITION(NE,NV,NM,IERROR)
      CALL EDDY_AMPLITUDE(NE,NV,NM)
   ENDDO EDDY_LOOP
   
   ! Cholesky decomposition of Reynolds stress tensor
   A_IJ => VT%A_IJ
   R_IJ => VT%R_IJ
   A_IJ = 0._EB
   A_IJ(1,1) = SQRT(R_IJ(1,1))
   A_IJ(2,1) = R_IJ(2,1)/A_IJ(1,1)
   A_IJ(2,2) = SQRT(R_IJ(2,2)-A_IJ(2,1)**2)
   A_IJ(3,1) = R_IJ(3,1)/A_IJ(1,1)
   A_IJ(3,2) = (R_IJ(3,2)-A_IJ(2,1)*A_IJ(3,1))/A_IJ(2,2)
   A_IJ(3,3) = SQRT(R_IJ(3,3)-A_IJ(3,1)**2-A_IJ(3,2)**2)

ENDDO VENT_LOOP

END SUBROUTINE SYNTHETIC_EDDY_SETUP


SUBROUTINE SYNTHETIC_TURBULENCE(DT,NM)
IMPLICIT NONE

REAL(EB), INTENT(IN) :: DT
INTEGER, INTENT(IN) :: NM
INTEGER :: NE,NV,II,JJ,KK,IERROR
TYPE(VENTS_TYPE), POINTER :: VT=>NULL()
TYPE(SURFACE_TYPE), POINTER :: SF=>NULL()
REAL(EB) :: XX,YY,ZZ,SHAPE_FACTOR,VOLUME_WEIGHTING_FACTOR

! Reference:
!
! Nicolas Jarrin. Synthetic Inflow Boundary Conditions for the Numerical Simulation of Turbulence. PhD Thesis,
! The University of Manchester, 2008.
!
! See Chapter 4: The Synthetic Eddy Method

VENT_LOOP: DO NV=1,MESHES(NM)%N_VENT
   VT => MESHES(NM)%VENTS(NV)
   IF (VT%N_EDDY==0) CYCLE VENT_LOOP
   
   VT%U_EDDY = 0._EB
   VT%V_EDDY = 0._EB
   VT%W_EDDY = 0._EB
   SF => SURFACE(VT%IBC)
   
   IOR_SELECT: SELECT CASE(ABS(VT%IOR))
      CASE(1)
         EDDY_LOOP_1: DO NE=1,VT%N_EDDY ! loop over eddies
            VT%X_EDDY(NE) = VT%X_EDDY(NE)-DT*SF%VEL*SIGN(1._EB,REAL(VT%IOR,EB))
            VT%Y_EDDY(NE) = VT%Y_EDDY(NE)+DT*SF%VEL_T(1)
            VT%Z_EDDY(NE) = VT%Z_EDDY(NE)+DT*SF%VEL_T(2)   
            IERROR=0;      CALL EDDY_POSITION(NE,NV,NM,IERROR)
            IF (IERROR==1) CALL EDDY_AMPLITUDE(NE,NV,NM)
            DO KK=VT%K1+1,VT%K2 ! this block can be made more efficient
               DO JJ=VT%J1+1,VT%J2
                  XX = (        VT%X1     - VT%X_EDDY(NE))/VT%L_EDDY
                  YY = (MESHES(NM)%YC(JJ) - VT%Y_EDDY(NE))/VT%L_EDDY
                  ZZ = (MESHES(NM)%ZC(KK) - VT%Z_EDDY(NE))/VT%L_EDDY
                  SHAPE_FACTOR = SHAPE_FUNCTION(XX,1)*SHAPE_FUNCTION(YY,1)*SHAPE_FUNCTION(ZZ,1)
                  VT%U_EDDY(JJ,KK) = VT%U_EDDY(JJ,KK) + VT%CU_EDDY(NE)*SHAPE_FACTOR
                  VT%V_EDDY(JJ,KK) = VT%V_EDDY(JJ,KK) + VT%CV_EDDY(NE)*SHAPE_FACTOR
                  VT%W_EDDY(JJ,KK) = VT%W_EDDY(JJ,KK) + VT%CW_EDDY(NE)*SHAPE_FACTOR
               ENDDO
            ENDDO
         ENDDO EDDY_LOOP_1
      CASE(2)
         EDDY_LOOP_2: DO NE=1,VT%N_EDDY
            VT%X_EDDY(NE) = VT%X_EDDY(NE)+DT*SF%VEL_T(2) 
            VT%Y_EDDY(NE) = VT%Y_EDDY(NE)-DT*SF%VEL*SIGN(1._EB,REAL(VT%IOR,EB))
            VT%Z_EDDY(NE) = VT%Z_EDDY(NE)+DT*SF%VEL_T(1)  
            IERROR=0;      CALL EDDY_POSITION(NE,NV,NM,IERROR)
            IF (IERROR==1) CALL EDDY_AMPLITUDE(NE,NV,NM)
            DO KK=VT%K1+1,VT%K2
               DO II=VT%I1+1,VT%I2
                  XX = (MESHES(NM)%XC(II) - VT%X_EDDY(NE))/VT%L_EDDY
                  YY = (        VT%Y1     - VT%Y_EDDY(NE))/VT%L_EDDY
                  ZZ = (MESHES(NM)%ZC(KK) - VT%Z_EDDY(NE))/VT%L_EDDY
                  SHAPE_FACTOR = SHAPE_FUNCTION(XX,1)*SHAPE_FUNCTION(YY,1)*SHAPE_FUNCTION(ZZ,1)
                  VT%U_EDDY(II,KK) = VT%U_EDDY(II,KK) + VT%CU_EDDY(NE)*SHAPE_FACTOR
                  VT%V_EDDY(II,KK) = VT%V_EDDY(II,KK) + VT%CV_EDDY(NE)*SHAPE_FACTOR
                  VT%W_EDDY(II,KK) = VT%W_EDDY(II,KK) + VT%CW_EDDY(NE)*SHAPE_FACTOR 
               ENDDO
            ENDDO  
         ENDDO EDDY_LOOP_2
      CASE(3)
         EDDY_LOOP_3: DO NE=1,VT%N_EDDY
            VT%X_EDDY(NE) = VT%X_EDDY(NE)+DT*SF%VEL_T(1)
            VT%Y_EDDY(NE) = VT%Y_EDDY(NE)+DT*SF%VEL_T(2)
            VT%Z_EDDY(NE) = VT%Z_EDDY(NE)-DT*SF%VEL*SIGN(1._EB,REAL(VT%IOR,EB))
            IERROR=0;      CALL EDDY_POSITION(NE,NV,NM,IERROR)
            IF (IERROR==1) CALL EDDY_AMPLITUDE(NE,NV,NM)
            DO JJ=VT%J1+1,VT%J2
               DO II=VT%I1+1,VT%I2
                  XX = (MESHES(NM)%XC(II) - VT%X_EDDY(NE))/VT%L_EDDY
                  YY = (MESHES(NM)%YC(JJ) - VT%Y_EDDY(NE))/VT%L_EDDY
                  ZZ = (        VT%Z1     - VT%Z_EDDY(NE))/VT%L_EDDY
                  SHAPE_FACTOR = SHAPE_FUNCTION(XX,1)*SHAPE_FUNCTION(YY,1)*SHAPE_FUNCTION(ZZ,1)
                  VT%U_EDDY(II,JJ) = VT%U_EDDY(II,JJ) + VT%CU_EDDY(NE)*SHAPE_FACTOR
                  VT%V_EDDY(II,JJ) = VT%V_EDDY(II,JJ) + VT%CV_EDDY(NE)*SHAPE_FACTOR
                  VT%W_EDDY(II,JJ) = VT%W_EDDY(II,JJ) + VT%CW_EDDY(NE)*SHAPE_FACTOR 
               ENDDO
            ENDDO  
         ENDDO EDDY_LOOP_3
   END SELECT IOR_SELECT
   
   VOLUME_WEIGHTING_FACTOR = MIN(1._EB,SQRT(VT%EDDY_BOX_VOLUME/REAL(VT%N_EDDY,EB)/VT%L_EDDY**3))
   ! note: L_EDDY included in SQRT based on Jung-il Choi write up.
   VT%U_EDDY = VT%U_EDDY*VOLUME_WEIGHTING_FACTOR
   VT%V_EDDY = VT%V_EDDY*VOLUME_WEIGHTING_FACTOR
   VT%W_EDDY = VT%W_EDDY*VOLUME_WEIGHTING_FACTOR

ENDDO VENT_LOOP

END SUBROUTINE SYNTHETIC_TURBULENCE


SUBROUTINE EDDY_POSITION(NE,NV,NM,IERROR)
IMPLICIT NONE

INTEGER, INTENT(IN) :: NE,NV,NM
INTEGER, INTENT(INOUT) :: IERROR
REAL(EB) :: RN
TYPE(VENTS_TYPE), POINTER :: VT=>NULL()

VT => MESHES(NM)%VENTS(NV)

IF (IERROR==0) THEN
   ! check to see if eddy is outside box
   IF (VT%X_EDDY(NE)<VT%X_EDDY_MIN .OR. VT%X_EDDY(NE)>VT%X_EDDY_MAX .OR. &
       VT%Y_EDDY(NE)<VT%Y_EDDY_MIN .OR. VT%Y_EDDY(NE)>VT%Y_EDDY_MAX .OR. &
       VT%Z_EDDY(NE)<VT%Z_EDDY_MIN .OR. VT%Z_EDDY(NE)>VT%Z_EDDY_MAX)       THEN
       
       IERROR=1 ! generate new positions and amplitudes (see EDDY_AMPLITUDE)
    ENDIF
ENDIF

IF (IERROR==1) THEN
    CALL RANDOM_NUMBER(RN); VT%X_EDDY(NE) = VT%X_EDDY_MIN + RN*(VT%X_EDDY_MAX-VT%X_EDDY_MIN)
    CALL RANDOM_NUMBER(RN); VT%Y_EDDY(NE) = VT%Y_EDDY_MIN + RN*(VT%Y_EDDY_MAX-VT%Y_EDDY_MIN)
    CALL RANDOM_NUMBER(RN); VT%Z_EDDY(NE) = VT%Z_EDDY_MIN + RN*(VT%Z_EDDY_MAX-VT%Z_EDDY_MIN)
ENDIF

END SUBROUTINE EDDY_POSITION


SUBROUTINE EDDY_AMPLITUDE(NE,NV,NM)
IMPLICIT NONE

INTEGER, INTENT(IN) :: NE,NV,NM
REAL(EB) :: RN,EPS_EDDY(3)
TYPE(VENTS_TYPE), POINTER :: VT=>NULL()
INTEGER :: J
    
EPS_EDDY=-1._EB
CALL RANDOM_NUMBER(RN); IF (RN>0.5_EB) EPS_EDDY(1)=1._EB
CALL RANDOM_NUMBER(RN); IF (RN>0.5_EB) EPS_EDDY(2)=1._EB
CALL RANDOM_NUMBER(RN); IF (RN>0.5_EB) EPS_EDDY(3)=1._EB

VT => MESHES(NM)%VENTS(NV)
VT%CU_EDDY(NE)=0._EB
VT%CV_EDDY(NE)=0._EB
VT%CW_EDDY(NE)=0._EB
! A_IJ is the Cholesky decomposition of R_IJ, see SYNTHETIC_EDDY_SETUP
DO J=1,3
   VT%CU_EDDY(NE)=VT%CU_EDDY(NE)+VT%A_IJ(1,J)*EPS_EDDY(J)
   VT%CV_EDDY(NE)=VT%CV_EDDY(NE)+VT%A_IJ(2,J)*EPS_EDDY(J)
   VT%CW_EDDY(NE)=VT%CW_EDDY(NE)+VT%A_IJ(3,J)*EPS_EDDY(J)
ENDDO

END SUBROUTINE EDDY_AMPLITUDE


REAL(EB) FUNCTION SHAPE_FUNCTION(X,CODE)
IMPLICIT NONE

REAL(EB), INTENT(IN) :: X
INTEGER, INTENT(IN) :: CODE

SHAPE_FUNCTION = 0._EB
SELECT CASE(CODE)
   CASE(1) ! tent function, Jarrin Eq. (4.59)
      IF (ABS(X)<1._EB) SHAPE_FUNCTION = SQRT(1.5_EB)*(1._EB-ABS(X))
   CASE(2) ! top hat function
      IF (ABS(X)<1._EB) SHAPE_FUNCTION = 0.707106781186547_EB ! 1/sqrt(2)
   !CASE(3) ! truncated Gaussian
   !   IF (ABS(X)<1._EB) SHAPE_FUNCTION = C*EXP(-4.5_EB*X**2)
END SELECT

END FUNCTION SHAPE_FUNCTION


SUBROUTINE GET_REV_turb(MODULE_REV,MODULE_DATE)
INTEGER,INTENT(INOUT) :: MODULE_REV
CHARACTER(255),INTENT(INOUT) :: MODULE_DATE

WRITE(MODULE_DATE,'(A)') turbrev(INDEX(turbrev,':')+1:LEN_TRIM(turbrev)-2)
READ (MODULE_DATE,'(I5)') MODULE_REV
WRITE(MODULE_DATE,'(A)') turbdate

END SUBROUTINE GET_REV_turb


END MODULE TURBULENCE

! NOTE: The embedded mesh module has been commented temporarily to allow cleanup for FDS v6 release.
! Development will continue at a later date.  But this module requires significant architectural
! changes for proper implementation and this would derail an already behind schedule v6 release. -RJM

!! This module is an experimental implementation of my embedded mesh method (EMB),
!! a prelude to adaptive mesh refinement.
!
!MODULE EMBEDDED_MESH_METHOD
!
!USE PRECISION_PARAMETERS
!USE GLOBAL_CONSTANTS
!USE MESH_VARIABLES
!USE MESH_POINTERS
!
!IMPLICIT NONE
!
!PRIVATE
!PUBLIC SCALARF_EMB,VELOCITY_EMB,RESTRICT_MASS_EMB,RESTRICT_DIV_EMB,PROJECT_VELOCITY, &
!       SORT_MESH_LEVEL,MATCH_VELOCITY_EMB,SCALAR_GHOST_EMB
! 
!CONTAINS
!
!
!SUBROUTINE SCALARF_EMB(NM1,NM2,IERROR)
!IMPLICIT NONE
!
!INTEGER, INTENT(IN) :: NM1,NM2
!
!TYPE(MESH_TYPE), POINTER :: M1,M2
!INTEGER :: N,I,J,K,I_LO,I_HI,J_LO,J_HI,K_LO,K_HI,II_0,JJ_0,KK_0,II,JJ,KK, &
!           NRX,NRY,NRZ,N2X,N2Y,N2Z,II_LO,JJ_LO,KK_LO,INDEX_LIST(12),IERROR
!REAL(EB) :: VOLUME_LIST(3)
!REAL(EB), POINTER, DIMENSION(:,:,:,:) :: FX1,FY1,FZ1,FX2,FY2,FZ2
!
!!   Comments:
!!
!!   Assumes uniform grid in each direction and that M2 lies within M1.
!!
!!   -------------------------------
!!   |         |         |         |
!!   |         |         |         |<---MESHES(M1)
!!   |         |         |         |
!!   |         |         |         |
!!   -------------------------------
!!   |         |-|-|-|-|-|         |
!!   |         |-|-|-|-|-|<-------------MESHES(M2)
!!   |         |-|-|-|-|-|         |
!!   |         |-|-|-|-|-|         |
!!   -------------------------------
!!   |         |         |         |
!!   |         |         |         |
!!   |         |         |         |
!!   |         |         |         |
!!   -------------------------------
!
!CALL LOCATE_MESH(INDEX_LIST,VOLUME_LIST,NM1,NM2,IERROR)
!SELECT CASE (IERROR)
!   CASE(0)
!      I_LO = INDEX_LIST(1)
!      I_HI = INDEX_LIST(2)
!      J_LO = INDEX_LIST(3)
!      J_HI = INDEX_LIST(4)
!      K_LO = INDEX_LIST(5)
!      K_HI = INDEX_LIST(6)
!      II_LO = INDEX_LIST(7)
!      JJ_LO = INDEX_LIST(8)
!      KK_LO = INDEX_LIST(9)
!      NRX = INDEX_LIST(10)
!      NRY = INDEX_LIST(11)
!      NRZ = INDEX_LIST(12)
!   CASE(1)
!      RETURN
!END SELECT
!
!M1=>MESHES(NM1) ! coarse mesh
!M2=>MESHES(NM2) ! fine mesh
!
!N2X = NRY*NRZ
!N2Y = NRX*NRZ
!N2Z = NRX*NRY
!
!! Fluxes
!
!FX1=>M1%SCALAR_SAVE1
!FY1=>M1%SCALAR_SAVE2
!FZ1=>M1%SCALAR_SAVE3
!
!FX2=>M2%SCALAR_SAVE1
!FY2=>M2%SCALAR_SAVE2
!FZ2=>M2%SCALAR_SAVE3
!
!! Restrict fine mesh to coarse mesh for embedded cells
!
!SPECIES_LOOP: DO N=0,N_SPECIES
!
!   ! x-direction fluxes
!
!   DO K = K_LO,K_HI
!      KK_0 = KK_LO + (K-K_LO)*NRZ
!      DO J = J_LO,J_HI
!         JJ_0 = JJ_LO + (J-J_LO)*NRY
!         DO I = I_LO-1,I_HI !! note: this includes fine mesh boundary
!            II_0 = II_LO + (I-I_LO+1)*NRX !!
!                  
!            FX1(I,J,K,N) = 0._EB
!            DO KK = KK_0+1,KK_0+NRZ
!               DO JJ = JJ_0+1,JJ_0+NRY
!                  FX1(I,J,K,N) = FX1(I,J,K,N) + FX2(II_0,JJ,KK,N)
!               ENDDO
!            ENDDO
!            FX1(I,J,K,N) = FX1(I,J,K,N)/N2X
!         
!         ENDDO
!      ENDDO
!   ENDDO
!   
!   ! y-direction fluxes
!
!   DO K = K_LO,K_HI
!      KK_0 = KK_LO + (K-K_LO)*NRZ
!      DO J = J_LO-1,J_HI !!
!         JJ_0 = JJ_LO + (J-J_LO+1)*NRY !!
!         DO I = I_LO,I_HI 
!            II_0 = II_LO + (I-I_LO)*NRX
!                  
!            FY1(I,J,K,N) = 0._EB
!            DO KK = KK_0+1,KK_0+NRZ
!               DO II = II_0+1,II_0+NRX
!                  FY1(I,J,K,N) = FY1(I,J,K,N) + FY2(II,JJ_0,KK,N)
!               ENDDO
!            ENDDO
!            FY1(I,J,K,N) = FY1(I,J,K,N)/N2Y
!         
!         ENDDO
!      ENDDO
!   ENDDO
!   
!   ! z-direction fluxes
!
!   DO K = K_LO-1,K_HI !!
!      KK_0 = KK_LO + (K-K_LO+1)*NRZ !!
!      DO J = J_LO,J_HI
!         JJ_0 = JJ_LO + (J-J_LO)*NRY
!         DO I = I_LO,I_HI 
!            II_0 = II_LO + (I-I_LO)*NRX
!                  
!            FZ1(I,J,K,N) = 0._EB
!            DO JJ = JJ_0+1,JJ_0+NRY
!               DO II = II_0+1,II_0+NRX
!                  FZ1(I,J,K,N) = FZ1(I,J,K,N) + FZ2(II,JJ,KK_0,N)
!               ENDDO
!            ENDDO
!            FZ1(I,J,K,N) = FZ1(I,J,K,N)/N2Z
!         
!         ENDDO
!      ENDDO
!   ENDDO
!   
!ENDDO SPECIES_LOOP
!
!END SUBROUTINE SCALARF_EMB
!
!
!SUBROUTINE VELOCITY_EMB(NM1,NM2,IERROR)
!IMPLICIT NONE
!
!INTEGER, INTENT(IN) :: NM1,NM2
!
!TYPE(MESH_TYPE), POINTER :: M1,M2
!INTEGER :: I,J,K,I_LO,I_HI,J_LO,J_HI,K_LO,K_HI,II_0,JJ_0,KK_0,II,JJ,KK, &
!           NRX,NRY,NRZ,N2X,N2Y,N2Z,II_LO,JJ_LO,KK_LO,INDEX_LIST(12),IERROR
!REAL(EB) :: VOLUME_LIST(3)
!REAL(EB), POINTER, DIMENSION(:,:,:) :: UU1,VV1,WW1,UU2,VV2,WW2
!
!CALL LOCATE_MESH(INDEX_LIST,VOLUME_LIST,NM1,NM2,IERROR)
!SELECT CASE (IERROR)
!   CASE(0)
!      I_LO = INDEX_LIST(1)
!      I_HI = INDEX_LIST(2)
!      J_LO = INDEX_LIST(3)
!      J_HI = INDEX_LIST(4)
!      K_LO = INDEX_LIST(5)
!      K_HI = INDEX_LIST(6)
!      II_LO = INDEX_LIST(7)
!      JJ_LO = INDEX_LIST(8)
!      KK_LO = INDEX_LIST(9)
!      NRX = INDEX_LIST(10)
!      NRY = INDEX_LIST(11)
!      NRZ = INDEX_LIST(12)
!   CASE(1)
!      RETURN
!END SELECT
!
!M1=>MESHES(NM1) ! coarse mesh
!M2=>MESHES(NM2) ! fine mesh
!
!N2X = NRY*NRZ
!N2Y = NRX*NRZ
!N2Z = NRX*NRY
!
!IF (PREDICTOR) THEN
!   UU1=>M1%U
!   VV1=>M1%V
!   WW1=>M1%W
!   UU2=>M2%U
!   VV2=>M2%V
!   WW2=>M2%W
!ELSEIF (CORRECTOR) THEN
!   UU1=>M1%US
!   VV1=>M1%VS
!   WW1=>M1%WS
!   UU2=>M2%US
!   VV2=>M2%VS
!   WW2=>M2%WS
!ENDIF
!
!! Restrict fine mesh to coarse mesh for embedded cells
!
!! U-VELOCITY
!
!DO K = K_LO,K_HI
!   KK_0 = KK_LO + (K-K_LO)*NRZ
!   DO J = J_LO,J_HI
!      JJ_0 = JJ_LO + (J-J_LO)*NRY
!      DO I = I_LO,I_HI-1 ! excludes boundary values
!         II_0 = II_LO + (I-I_LO+1)*NRX
!                  
!         UU1(I,J,K) = 0._EB
!         DO KK = KK_0+1,KK_0+NRZ
!            DO JJ = JJ_0+1,JJ_0+NRY
!               UU1(I,J,K) = UU1(I,J,K) + UU2(II_0,JJ,KK)
!            ENDDO
!         ENDDO
!         UU1(I,J,K) = UU1(I,J,K)/N2X
!         
!      ENDDO
!   ENDDO
!ENDDO
!   
!! V-VELOCITY
!
!DO K = K_LO,K_HI
!   KK_0 = KK_LO + (K-K_LO)*NRZ
!   DO J = J_LO,J_HI-1 ! excludes boundary values
!      JJ_0 = JJ_LO + (J-J_LO+1)*NRY
!      DO I = I_LO,I_HI 
!         II_0 = II_LO + (I-I_LO)*NRX
!                  
!         VV1(I,J,K) = 0._EB
!         DO KK = KK_0+1,KK_0+NRZ
!            DO II = II_0+1,II_0+NRX
!               VV1(I,J,K) = VV1(I,J,K) + VV2(II,JJ_0,KK)
!            ENDDO
!         ENDDO
!         VV1(I,J,K) = VV1(I,J,K)/N2Y
!         
!      ENDDO
!   ENDDO
!ENDDO
!   
!! W-VELOCITY
!
!DO K = K_LO,K_HI-1 ! excludes boundary values
!   KK_0 = KK_LO + (K-K_LO+1)*NRZ
!   DO J = J_LO,J_HI
!      JJ_0 = JJ_LO + (J-J_LO)*NRY
!      DO I = I_LO,I_HI 
!         II_0 = II_LO + (I-I_LO)*NRX
!                  
!         WW1(I,J,K) = 0._EB
!         DO JJ = JJ_0+1,JJ_0+NRY
!            DO II = II_0+1,II_0+NRX
!               WW1(I,J,K) = WW1(I,J,K) + WW2(II,JJ,KK_0)
!            ENDDO
!         ENDDO
!         WW1(I,J,K) = WW1(I,J,K)/N2Z
!         
!      ENDDO
!   ENDDO
!ENDDO
!
!END SUBROUTINE VELOCITY_EMB
!
!
!SUBROUTINE RESTRICT_MASS_EMB(NM1,NM2,IERROR)
!IMPLICIT NONE
!
!INTEGER, INTENT(IN) :: NM1,NM2
!
!TYPE(MESH_TYPE), POINTER :: M1,M2
!INTEGER :: N,I,J,K,I_LO,I_HI,J_LO,J_HI,K_LO,K_HI,II_0,JJ_0,KK_0,II,JJ,KK, &
!           NRX,NRY,NRZ, II_LO,JJ_LO,KK_LO,INDEX_LIST(12),IERROR
!REAL(EB) :: DV1,DV2,DVRAT,VOLUME_LIST(3)
!REAL(EB), POINTER, DIMENSION(:,:,:) :: RHO1,RHO2
!REAL(EB), POINTER, DIMENSION(:,:,:,:) :: YY1,YY2
!
!CALL LOCATE_MESH(INDEX_LIST,VOLUME_LIST,NM1,NM2,IERROR)
!SELECT CASE (IERROR)
!   CASE(0)
!      I_LO = INDEX_LIST(1)
!      I_HI = INDEX_LIST(2)
!      J_LO = INDEX_LIST(3)
!      J_HI = INDEX_LIST(4)
!      K_LO = INDEX_LIST(5)
!      K_HI = INDEX_LIST(6)
!      II_LO = INDEX_LIST(7)
!      JJ_LO = INDEX_LIST(8)
!      KK_LO = INDEX_LIST(9)
!      NRX = INDEX_LIST(10)
!      NRY = INDEX_LIST(11)
!      NRZ = INDEX_LIST(12)
!      DV1 = VOLUME_LIST(1)
!      DV2 = VOLUME_LIST(2)
!      DVRAT = VOLUME_LIST(3)
!   CASE(1)
!      RETURN
!END SELECT
!
!M1=>MESHES(NM1) ! coarse mesh
!M2=>MESHES(NM2) ! fine mesh
!
!IF (PREDICTOR) THEN
!   RHO1 => M1%RHOS
!   RHO2 => M2%RHOS
!   IF (N_SPECIES>0) YY1  => M1%YYS
!   IF (N_SPECIES>0) YY2  => M2%YYS
!ELSEIF (CORRECTOR) THEN
!   RHO1 => M1%RHO
!   RHO2 => M2%RHO
!   IF (N_SPECIES>0) YY1  => M1%YY
!   IF (N_SPECIES>0) YY2  => M2%YY
!ENDIF
!
!DO K = K_LO,K_HI
!   KK_0 = KK_LO + (K-K_LO)*NRZ
!   DO J = J_LO,J_HI
!      JJ_0 = JJ_LO + (J-J_LO)*NRY
!      DO I = I_LO,I_HI
!         II_0 = II_LO + (I-I_LO)*NRX
!            
!         RHO1(I,J,K) = 0._EB
!         
!         DO KK = KK_0+1,KK_0+NRZ
!            DO JJ = JJ_0+1,JJ_0+NRY
!               DO II = II_0+1,II_0+NRX
!                 
!                  RHO1(I,J,K) = RHO1(I,J,K) + RHO2(II,JJ,KK)*DVRAT
!                     
!               ENDDO
!            ENDDO
!         ENDDO
!      
!      ENDDO
!   ENDDO
!ENDDO
!
!IF (N_SPECIES>0) THEN
!
!   SPECIES_LOOP: DO N=1,N_SPECIES
!   
!      DO K = K_LO,K_HI
!         KK_0 = KK_LO + (K-K_LO)*NRZ
!         DO J = J_LO,J_HI
!            JJ_0 = JJ_LO + (J-J_LO)*NRY
!            DO I = I_LO,I_HI
!               II_0 = II_LO + (I-I_LO)*NRX
!            
!               YY1(I,J,K,N) = 0._EB
!         
!               DO KK = KK_0+1,KK_0+NRZ
!                  DO JJ = JJ_0+1,JJ_0+NRY
!                     DO II = II_0+1,II_0+NRX
!                 
!                        YY1(I,J,K,N) = YY1(I,J,K,N) + RHO2(II,JJ,KK)*YY2(II,JJ,KK,N)*DV2
!                     
!                     ENDDO
!                  ENDDO
!               ENDDO
!               
!               YY1(I,J,K,N) = YY1(I,J,K,N)/(RHO1(I,J,K)*DV1)
!      
!            ENDDO
!         ENDDO
!      ENDDO
!
!   ENDDO SPECIES_LOOP
!
!ENDIF
!
!END SUBROUTINE RESTRICT_MASS_EMB
!
!
!SUBROUTINE RESTRICT_DIV_EMB(NM1,NM2,IERROR)
!IMPLICIT NONE
!
!INTEGER, INTENT(IN) :: NM1,NM2
!
!TYPE(MESH_TYPE), POINTER :: M1,M2
!INTEGER :: I,J,K,I_LO,I_HI,J_LO,J_HI,K_LO,K_HI,II_0,JJ_0,KK_0,II,JJ,KK, &
!           NRX,NRY,NRZ,II_LO,JJ_LO,KK_LO,INDEX_LIST(12),IERROR
!REAL(EB) :: DVRAT,VOLUME_LIST(3)
!REAL(EB), POINTER, DIMENSION(:,:,:) :: DP1,DP2
!
!CALL LOCATE_MESH(INDEX_LIST,VOLUME_LIST,NM1,NM2,IERROR)
!SELECT CASE (IERROR)
!   CASE(0)
!      I_LO = INDEX_LIST(1)
!      I_HI = INDEX_LIST(2)
!      J_LO = INDEX_LIST(3)
!      J_HI = INDEX_LIST(4)
!      K_LO = INDEX_LIST(5)
!      K_HI = INDEX_LIST(6)
!      II_LO = INDEX_LIST(7)
!      JJ_LO = INDEX_LIST(8)
!      KK_LO = INDEX_LIST(9)
!      NRX = INDEX_LIST(10)
!      NRY = INDEX_LIST(11)
!      NRZ = INDEX_LIST(12)
!      DVRAT = VOLUME_LIST(3)
!   CASE(1)
!      RETURN
!END SELECT
!
!M1=>MESHES(NM1) ! coarse mesh
!M2=>MESHES(NM2) ! fine mesh
!
!IF (PREDICTOR) THEN
!   DP1 => M1%DS
!   DP2 => M2%DS
!ELSEIF (CORRECTOR) THEN
!   DP1 => M1%DDDT
!   DP2 => M2%DDDT
!ENDIF
!
!! Restrict divergence
!   
!DO K = K_LO,K_HI
!   KK_0 = KK_LO + (K-K_LO)*NRZ
!   DO J = J_LO,J_HI
!      JJ_0 = JJ_LO + (J-J_LO)*NRY
!      DO I = I_LO,I_HI
!         II_0 = II_LO + (I-I_LO)*NRX
!            
!         DP1(I,J,K) = 0._EB
!         
!         DO KK = KK_0+1,KK_0+NRZ
!            DO JJ = JJ_0+1,JJ_0+NRY
!               DO II = II_0+1,II_0+NRX
!                 
!                  DP1(I,J,K) = DP1(I,J,K) + DP2(II,JJ,KK)*DVRAT
!                     
!               ENDDO
!            ENDDO
!         ENDDO
!      
!      ENDDO
!   ENDDO
!ENDDO
!
!END SUBROUTINE RESTRICT_DIV_EMB
!
!
!SUBROUTINE PROJECT_VELOCITY(NM)
!USE POIS, ONLY: H3CZSS,H2CZSS
!IMPLICIT NONE
!
!INTEGER, INTENT(IN) :: NM
!INTEGER :: I,J,K
!REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW,DP,PP,PRHS_SAVE
!REAL(EB) :: DIV,LHSS,RES,POIS_ERR
!
!CALL POINT_TO_MESH(NM)
!
!IF (PREDICTOR) THEN
!   ! note: PROJECT_VELOCITY is called AFTER the predictor update of velocity
!   UU=>US
!   VV=>VS
!   WW=>WS
!   DP=>D
!ELSEIF (CORRECTOR) THEN
!   UU=>U
!   VV=>V
!   WW=>W
!   DP=>DS
!ENDIF
!PP=>WORK1
!PRHS_SAVE=>WORK2
!
!! build source
!
!DO K=1,KBAR
!   DO J=1,JBAR
!      DO I=1,IBAR
!         DIV = (UU(I,J,K)-UU(I-1,J,K))*RDX(I) + (VV(I,J,K)-VV(I,J-1,K))*RDY(J) + (WW(I,J,K)-WW(I,J,K-1))*RDZ(K)
!         PRHS(I,J,K) = DIV-DP(I,J,K)
!      ENDDO
!   ENDDO
!ENDDO
!
!! solve Poisson equation
!
!BXS = 0._EB
!BXF = 0._EB
!BYS = 0._EB
!BYF = 0._EB
!BZS = 0._EB
!BZF = 0._EB
!
!PRHS_SAVE(1:IBAR,1:JBAR,1:KBAR) = PRHS
!IF (.NOT.TWO_D) CALL H3CZSS(BXS,BXF,BYS,BYF,BZS,BZF,ITRN,JTRN,PRHS,POIS_PTB,SAVE2,WORK,HX)
!IF (TWO_D .AND. .NOT. CYLINDRICAL) CALL H2CZSS(BXS,BXF,BZS,BZF,ITRN,PRHS,POIS_PTB,SAVE2,WORK,HX)
!PP(1:IBAR,1:JBAR,1:KBAR) = PRHS
!
!! Apply boundary conditions to PP
! 
!DO K=1,KBAR
!   DO J=1,JBAR
!      PP(0,J,K)    = -PP(1,J,K) ! use minus if Dirichlet, plus if Neumann, see init of SAVE2
!      PP(IBP1,J,K) = -PP(IBAR,J,K)
!   ENDDO
!ENDDO
! 
!DO K=1,KBAR
!   DO I=1,IBAR
!      PP(I,0,K)    = -PP(I,1,K)
!      PP(I,JBP1,K) = -PP(I,JBAR,K)
!   ENDDO
!ENDDO
!
!DO J=1,JBAR
!   DO I=1,IBAR
!      PP(I,J,0)    = -PP(I,J,1)
!      PP(I,J,KBP1) = -PP(I,J,KBAR)
!   ENDDO
!ENDDO
!
!! ************************* Check the Solution *************************
! 
!IF (.FALSE.) THEN
!
!   POIS_ERR = 0._EB
!   DO K=1,KBAR
!      DO J=1,JBAR
!         DO I=1,IBAR
!            LHSS = ((PP(I+1,J,K)-PP(I,J,K))*RDXN(I) - (PP(I,J,K)-PP(I-1,J,K))*RDXN(I-1) )*RDX(I) &
!                 + ((PP(I,J+1,K)-PP(I,J,K))*RDYN(J) - (PP(I,J,K)-PP(I,J-1,K))*RDYN(J-1) )*RDY(J) &
!                 + ((PP(I,J,K+1)-PP(I,J,K))*RDZN(K) - (PP(I,J,K)-PP(I,J,K-1))*RDZN(K-1) )*RDZ(K)
!            RES = ABS(PRHS_SAVE(I,J,K)-LHSS)
!            POIS_ERR = MAX(RES,POIS_ERR)
!         ENDDO
!      ENDDO
!   ENDDO
!   WRITE(0,*) 'POIS ERROR:',pois_ptb,pois_err
!
!ENDIF
!
!! correct velocities
!
!DO K=1,KBAR
!   DO J=1,JBAR
!      DO I=0,IBAR
!         UU(I,J,K) = UU(I,J,K) - RDXN(I)*(PP(I+1,J,K)-PP(I,J,K))
!      ENDDO
!   ENDDO
!ENDDO
!
!DO K=1,KBAR
!   DO J=0,JBAR
!      DO I=1,IBAR
!         VV(I,J,K) = VV(I,J,K) - RDYN(J)*(PP(I,J+1,K)-PP(I,J,K))
!      ENDDO
!   ENDDO
!ENDDO
!
!DO K=0,KBAR
!   DO J=1,JBAR
!      DO I=1,IBAR
!         WW(I,J,K) = WW(I,J,K) - RDZN(K)*(PP(I,J,K+1)-PP(I,J,K))
!      ENDDO
!   ENDDO
!ENDDO
!
!! check divergence
!
!IF (.FALSE.) THEN
!   POIS_ERR = 0._EB
!   DO K=1,KBAR
!      DO J=1,JBAR
!         DO I=1,IBAR
!            DIV = (UU(I,J,K)-UU(I-1,J,K))*RDX(I) + (VV(I,J,K)-VV(I,J-1,K))*RDY(J) + (WW(I,J,K)-WW(I,J,K-1))*RDZ(K)
!            RES = ABS(DIV-DP(I,J,K))
!            POIS_ERR = MAX(RES,POIS_ERR)
!         ENDDO
!      ENDDO
!   ENDDO
!   WRITE(0,*) NM,MAXVAL(ABS(DP)),POIS_ERR
!
!ENDIF
!
!END SUBROUTINE PROJECT_VELOCITY
!
!
!SUBROUTINE SORT_MESH_LEVEL
!IMPLICIT NONE
!
!INTEGER :: IRANK,NM,ML,MLMIN,MLMAX
!
!MESH_LIST_EMB = 0
!
!MLMAX = MAXVAL(MESHES(1:NMESHES)%MESH_LEVEL)
!MLMIN = MINVAL(MESHES(1:NMESHES)%MESH_LEVEL)
!
!IRANK=0
!
!DO ML=MLMAX,MLMIN,-1
!   DO NM=1,NMESHES
!   
!      IF (MESHES(NM)%MESH_LEVEL==ML) THEN
!         IRANK=IRANK+1
!         MESH_LIST_EMB(IRANK)=NM
!      ENDIF
!
!   ENDDO
!ENDDO
!
!!PRINT *,MLMIN,MLMAX
!!DO IRANK=1,NMESHES
!!   PRINT *,MESH_LIST_EMB(IRANK)
!!ENDDO
!!STOP
!
!END SUBROUTINE SORT_MESH_LEVEL
!
!
!SUBROUTINE MATCH_VELOCITY_EMB(NM1,NM2,IERROR)
!IMPLICIT NONE
!
!INTEGER, INTENT(IN) :: NM1,NM2
!
!TYPE(MESH_TYPE), POINTER :: M1,M2
!INTEGER :: I,J,K,I_LO,I_HI,J_LO,J_HI,K_LO,K_HI,II_0,JJ_0,KK_0,II,JJ,KK, &
!           NRX,NRY,NRZ,II_LO,JJ_LO,KK_LO,INDEX_LIST(12),IERROR,IW,IOR
!REAL(EB) :: VOLUME_LIST(3)
!REAL(EB), POINTER, DIMENSION(:,:,:) :: UU1,VV1,WW1,UU2,VV2,WW2
!REAL(EB), PARAMETER :: RF=0.5_EB,OMRF=0.5_EB
!
!CALL LOCATE_MESH(INDEX_LIST,VOLUME_LIST,NM1,NM2,IERROR)
!SELECT CASE (IERROR)
!   CASE(0)
!      I_LO = INDEX_LIST(1)
!      I_HI = INDEX_LIST(2)
!      J_LO = INDEX_LIST(3)
!      J_HI = INDEX_LIST(4)
!      K_LO = INDEX_LIST(5)
!      K_HI = INDEX_LIST(6)
!      II_LO = INDEX_LIST(7)
!      JJ_LO = INDEX_LIST(8)
!      KK_LO = INDEX_LIST(9)
!      NRX = INDEX_LIST(10)
!      NRY = INDEX_LIST(11)
!      NRZ = INDEX_LIST(12)
!   CASE(1)
!      RETURN
!END SELECT
!
!M1=>MESHES(NM1) ! coarse mesh
!M2=>MESHES(NM2) ! fine mesh
!
!IF (PREDICTOR) THEN
!   UU1=>M1%US
!   VV1=>M1%VS
!   WW1=>M1%WS
!   UU2=>M2%US
!   VV2=>M2%VS
!   WW2=>M2%WS
!ELSEIF (CORRECTOR) THEN
!   UU1=>M1%U
!   VV1=>M1%V
!   WW1=>M1%W
!   UU2=>M2%U
!   VV2=>M2%V
!   WW2=>M2%W
!ENDIF
!
!! Set fine mesh boundary value to corresponding coarse mesh value
!
!! U-VELOCITY
!
!DO K = K_LO,K_HI
!   KK_0 = KK_LO + (K-K_LO)*NRZ
!   DO J = J_LO,J_HI
!      JJ_0 = JJ_LO + (J-J_LO)*NRY
!
!      ! east face
!      I = I_HI
!      II_0 = II_LO + (I-I_LO+1)*NRX
!      IF (II_0==M2%IBAR) THEN
!         DO KK = KK_0+1,KK_0+NRZ
!            DO JJ = JJ_0+1,JJ_0+NRY
!               UU2(II_0,JJ,KK) = UU1(I,J,K)
!            ENDDO
!         ENDDO
!      ENDIF
!         
!      ! west face
!      I = I_LO-1
!      II_0 = II_LO + (I-I_LO+1)*NRX
!      IF (II_0==0) THEN
!         DO KK = KK_0+1,KK_0+NRZ
!            DO JJ = JJ_0+1,JJ_0+NRY
!               UU2(II_0,JJ,KK) = UU1(I,J,K)
!            ENDDO
!         ENDDO
!      ENDIF
!         
!   ENDDO
!ENDDO
!   
!! V-VELOCITY
!
!DO K = K_LO,K_HI
!   KK_0 = KK_LO + (K-K_LO)*NRZ
!   DO I = I_LO,I_HI
!      II_0 = II_LO + (I-I_LO)*NRX
!
!      ! north face
!      J = J_HI
!      JJ_0 = JJ_LO + (J-J_LO+1)*NRY
!      IF (JJ_0==M2%JBAR) THEN
!         DO KK = KK_0+1,KK_0+NRZ
!            DO II = II_0+1,II_0+NRX
!               VV2(II,JJ_0,KK) = VV1(I,J,K)
!            ENDDO
!         ENDDO
!      ENDIF
!         
!      ! south face
!      J = J_LO-1
!      JJ_0 = JJ_LO + (J-J_LO+1)*NRY
!      IF (JJ_0==0) THEN
!         DO KK = KK_0+1,KK_0+NRZ
!            DO II = II_0+1,II_0+NRX
!               VV2(II,JJ_0,KK) = VV1(I,J,K)
!            ENDDO
!         ENDDO
!      ENDIF
!         
!   ENDDO
!ENDDO
!   
!! W-VELOCITY
!
!DO J = J_LO,J_HI
!   JJ_0 = JJ_LO + (J-J_LO)*NRY
!   DO I = I_LO,I_HI
!      II_0 = II_LO + (I-I_LO)*NRX
!
!      ! top face
!      K = K_HI
!      KK_0 = KK_LO + (K-K_LO+1)*NRZ
!      IF (KK_0==M2%KBAR) THEN
!         DO JJ = JJ_0+1,JJ_0+NRY
!            DO II = II_0+1,II_0+NRX
!               WW2(II,JJ,KK_0) = WW1(I,J,K)
!            ENDDO
!         ENDDO
!      ENDIF
!         
!      ! bottom face
!      K = K_LO-1
!      KK_0 = KK_LO + (K-K_LO+1)*NRZ
!      IF (KK_0==0) THEN
!         DO JJ = JJ_0+1,JJ_0+NRY
!            DO II = II_0+1,II_0+NRX
!               WW2(II,JJ,KK_0) = WW1(I,J,K)
!            ENDDO
!         ENDDO
!      ENDIF
!         
!   ENDDO
!ENDDO
!
!! fine mesh boundary loop
!
!FINE_MESH_WALL_LOOP: DO IW=1,M2%NEWC
!   II  = M2%IJKW(1,IW)
!   JJ  = M2%IJKW(2,IW)
!   KK  = M2%IJKW(3,IW)
!   IOR = M2%IJKW(4,IW)
!   SELECT CASE (IOR)
!      CASE(1)
!         M2%UVW_SAVE(IW)=UU2(0,JJ,KK)
!      CASE(-1)
!         M2%UVW_SAVE(IW)=UU2(M2%IBAR,JJ,KK)
!      CASE(2)
!         M2%UVW_SAVE(IW)=UU2(II,0,KK)
!      CASE(-2)
!         M2%UVW_SAVE(IW)=UU2(II,M2%JBAR,KK)
!      CASE(3)
!         M2%UVW_SAVE(IW)=UU2(II,JJ,0)
!      CASE(-3)
!         M2%UVW_SAVE(IW)=UU2(II,JJ,M2%KBAR)
!   END SELECT
!ENDDO FINE_MESH_WALL_LOOP
!
!END SUBROUTINE MATCH_VELOCITY_EMB
!
!
!SUBROUTINE SCALAR_GHOST_EMB(NM1,NM2,IERROR)
!IMPLICIT NONE
!
!INTEGER, INTENT(IN) :: NM1,NM2
!
!TYPE(MESH_TYPE), POINTER :: M1,M2
!INTEGER :: N,I,J,K,I_LO,I_HI,J_LO,J_HI,K_LO,K_HI,II_0,JJ_0,KK_0,II,JJ,KK, &
!           NRX,NRY,NRZ,II_LO,JJ_LO,KK_LO,INDEX_LIST(12),IERROR,IW
!REAL(EB) :: VOLUME_LIST(3)
!REAL(EB), POINTER, DIMENSION(:,:,:) :: RHOP1,RHOP2,TMP1,TMP2,HH1,HH2
!REAL(EB), POINTER, DIMENSION(:,:,:,:) :: YYP1,YYP2
!
!CALL LOCATE_MESH(INDEX_LIST,VOLUME_LIST,NM1,NM2,IERROR)
!SELECT CASE (IERROR)
!   CASE(0)
!      I_LO = INDEX_LIST(1)
!      I_HI = INDEX_LIST(2)
!      J_LO = INDEX_LIST(3)
!      J_HI = INDEX_LIST(4)
!      K_LO = INDEX_LIST(5)
!      K_HI = INDEX_LIST(6)
!      II_LO = INDEX_LIST(7)
!      JJ_LO = INDEX_LIST(8)
!      KK_LO = INDEX_LIST(9)
!      NRX = INDEX_LIST(10)
!      NRY = INDEX_LIST(11)
!      NRZ = INDEX_LIST(12)
!   CASE(1)
!      RETURN
!END SELECT
!
!M1=>MESHES(NM1) ! coarse mesh
!M2=>MESHES(NM2) ! fine mesh
!
!IF (PREDICTOR) THEN
!   RHOP1 => M1%RHOS
!   YYP1  => M1%YYS
!   HH1   => M1%H
!   
!   RHOP2 => M2%RHOS
!   YYP2  => M2%YYS
!   HH2   => M2%H
!ELSEIF (CORRECTOR) THEN
!   RHOP1 => M1%RHO
!   YYP1  => M1%YY
!   HH1   => M1%HS
!   
!   RHOP2 => M2%RHO
!   YYP2  => M2%YY
!   HH2   => M2%HS
!ENDIF
!TMP1 => M1%TMP
!TMP2 => M2%TMP
!
!
!! Set fine mesh boundary value to corresponding coarse mesh value
!
!SPECIES_LOOP: DO N=1,N_SPECIES
!
!   DO K = K_LO,K_HI
!      KK_0 = KK_LO + (K-K_LO)*NRZ
!      DO J = J_LO,J_HI
!         JJ_0 = JJ_LO + (J-J_LO)*NRY
!
!         ! east face
!         I = I_HI+1
!         II_0 = II_LO + (I-I_LO)*NRX + 1
!         IF (II_0==M2%IBP1 .AND. I_HI/=M1%IBAR) THEN
!         ! if I_HI==M1%IBAR then this might be an external boundary and the ghost cell value
!         ! is handled in WALL_BC; similar conditions apply below
!            DO KK = KK_0+1,KK_0+NRZ
!               DO JJ = JJ_0+1,JJ_0+NRY
!                  RHOP2(II_0,JJ,KK) = RHOP1(I,J,K)
!                  TMP2(II_0,JJ,KK) = TMP1(I,J,K)
!                  HH2(II_0,JJ,KK) = HH1(I,J,K)
!                  IF (N_SPECIES>0) YYP2(II_0,JJ,KK,N) = YYP1(I,J,K,N)
!               ENDDO
!            ENDDO
!         ENDIF
!         
!         ! west face
!         I = I_LO-1
!         II_0 = II_LO + (I-I_LO+1)*NRX
!         IF (II_0==0 .AND. I_LO/=1) THEN
!            DO KK = KK_0+1,KK_0+NRZ
!               DO JJ = JJ_0+1,JJ_0+NRY
!                  RHOP2(II_0,JJ,KK) = RHOP1(I,J,K)
!                  TMP2(II_0,JJ,KK) = TMP1(I,J,K)
!                  HH2(II_0,JJ,KK) = HH1(I,J,K)
!                  IF (N_SPECIES>0) YYP2(II_0,JJ,KK,N) = YYP1(I,J,K,N)
!               ENDDO
!            ENDDO
!         ENDIF
!         
!      ENDDO
!   ENDDO
!
!   DO K = K_LO,K_HI
!      KK_0 = KK_LO + (K-K_LO)*NRZ
!      DO I = I_LO,I_HI
!         II_0 = II_LO + (I-I_LO)*NRX
!
!         ! north face
!         J = J_HI+1
!         JJ_0 = JJ_LO + (J-J_LO)*NRY + 1
!         IF (JJ_0==M2%JBP1 .AND. J_HI/=M1%JBAR) THEN
!            DO KK = KK_0+1,KK_0+NRZ
!               DO II = II_0+1,II_0+NRX
!                  RHOP2(II,JJ_0,KK) = RHOP1(I,J,K)
!                  TMP2(II,JJ_0,KK) = TMP1(I,J,K)
!                  HH2(II,JJ_0,KK) = HH1(I,J,K)
!                  IF (N_SPECIES>0) YYP2(II,JJ_0,KK,N) = YYP1(I,J,K,N)
!               ENDDO
!            ENDDO
!         ENDIF
!         
!         ! south face
!         J = J_LO-1
!         JJ_0 = JJ_LO + (J-J_LO+1)*NRY
!         IF (JJ_0==0 .AND. J_LO/=1) THEN
!            DO KK = KK_0+1,KK_0+NRZ
!               DO II = II_0+1,II_0+NRX
!                  RHOP2(II,JJ_0,KK) = RHOP1(I,J,K)
!                  TMP2(II,JJ_0,KK) = TMP1(I,J,K)
!                  HH2(II,JJ_0,KK) = HH1(I,J,K)
!                  IF (N_SPECIES>0) YYP2(II,JJ_0,KK,N) = YYP1(I,J,K,N)
!               ENDDO
!            ENDDO
!         ENDIF
!         
!      ENDDO
!   ENDDO
!
!   DO J = J_LO,J_HI
!      JJ_0 = JJ_LO + (J-J_LO)*NRY
!      DO I = I_LO,I_HI
!         II_0 = II_LO + (I-I_LO)*NRX
!
!         ! top face
!         K = K_HI+1
!         KK_0 = KK_LO + (K-K_LO)*NRZ + 1
!         IF (KK_0==M2%KBP1  .AND. K_HI/=M1%KBAR) THEN
!            DO JJ = JJ_0+1,JJ_0+NRY
!               DO II = II_0+1,II_0+NRX
!                  RHOP2(II,JJ,KK_0) = RHOP1(I,J,K)
!                  TMP2(II,JJ,KK_0) = TMP1(I,J,K)
!                  HH2(II,JJ,KK_0) = HH1(I,J,K)
!                  IF (N_SPECIES>0) YYP2(II,JJ,KK_0,N) = YYP1(I,J,K,N)
!               ENDDO
!            ENDDO
!         ENDIF
!         
!         ! bottom face
!         K = K_LO-1
!         KK_0 = KK_LO + (K-K_LO+1)*NRZ
!         IF (KK_0==0 .AND. K_LO/=1) THEN
!            DO JJ = JJ_0+1,JJ_0+NRY
!               DO II = II_0+1,II_0+NRX
!                  RHOP2(II,JJ,KK_0) = RHOP1(I,J,K)
!                  TMP2(II,JJ,KK_0) = TMP1(I,J,K)
!                  HH2(II,JJ,KK_0) = HH1(I,J,K)
!                  IF (N_SPECIES>0) YYP2(II,JJ,KK_0,N) = YYP1(I,J,K,N)
!               ENDDO
!            ENDDO
!         ENDIF
!         
!      ENDDO
!   ENDDO
!
!   WALL_LOOP: DO IW=1,M2%NWC
!      IF (M2%BOUNDARY_TYPE(IW)/=INTERPOLATED_BOUNDARY) CYCLE WALL_LOOP
!      II = M2%IJKW(1,IW)
!      JJ = M2%IJKW(2,IW)
!      KK = M2%IJKW(3,IW)
!      M2%RHO_F(IW)  = RHOP2(II,JJ,KK) 
!      M2%YY_F(IW,N) = YYP2(II,JJ,KK,N)
!   ENDDO WALL_LOOP
!
!ENDDO SPECIES_LOOP
!
!END SUBROUTINE SCALAR_GHOST_EMB
!
!
!SUBROUTINE LOCATE_MESH(INDEX_LIST,VOLUME_LIST,NM1,NM2,IERROR)
!IMPLICIT NONE
!
!INTEGER, INTENT(IN) :: NM1,NM2
!
!INTEGER, INTENT(OUT) :: IERROR,INDEX_LIST(12)
!REAL(EB), INTENT(OUT) :: VOLUME_LIST(3)
!TYPE (MESH_TYPE), POINTER :: M1,M2
!
!INTEGER :: I_LO,I_HI,J_LO,J_HI,K_LO,K_HI,II_LO,JJ_LO,KK_LO,NRX,NRY,NRZ
!REAL(EB) :: DV1,DV2,DVRAT
!
!IERROR=0
!INDEX_LIST=0
!VOLUME_LIST=0._EB
!
!M1=>MESHES(NM1) ! coarse mesh
!M2=>MESHES(NM2) ! fine mesh
!
!! Locate fine mesh within coarse mesh
!
!I_LO = MAX(1,       NINT((M2%XS-M1%XS)/M1%DX(1))+1 )
!I_HI = MIN(M1%IBAR, NINT((M2%XF-M1%XS)/M1%DX(1))   )
!IF (I_LO>M1%IBAR .OR. I_HI<1) THEN ! meshes do not overlap
!   IERROR=1
!   RETURN
!ENDIF
!
!J_LO = MAX(1,       NINT((M2%YS-M1%YS)/M1%DY(1))+1 )
!J_HI = MIN(M1%JBAR, NINT((M2%YF-M1%YS)/M1%DY(1))   )
!IF (J_LO>M1%JBAR .OR. J_HI<1) THEN ! meshes do not overlap
!   IERROR=1
!   RETURN
!ENDIF
!
!K_LO = MAX(1,       NINT((M2%ZS-M1%ZS)/M1%DZ(1))+1 )
!K_HI = MIN(M1%KBAR, NINT((M2%ZF-M1%ZS)/M1%DZ(1))   )
!IF (K_LO>M1%KBAR .OR. K_HI<1) THEN ! meshes do not overlap
!   IERROR=1
!   RETURN
!ENDIF
!
!! Find fine mesh off-set
!
!II_LO = MAX(0, NINT((M1%XS-M2%XS)/M2%DX(1)) )
!JJ_LO = MAX(0, NINT((M1%YS-M2%YS)/M2%DY(1)) )
!KK_LO = MAX(0, NINT((M1%ZS-M2%ZS)/M2%DZ(1)) )
!
!! Compute refinment ratio in each direction
!
!NRX = NINT(M1%DX(1)/M2%DX(1))
!NRY = NINT(M1%DY(1)/M2%DY(1))
!NRZ = NINT(M1%DZ(1)/M2%DZ(1))
!
!! Cell volumes
!
!DV1 = M1%DX(1)*M1%DY(1)*M1%DZ(1)
!DV2 = M2%DX(1)*M2%DY(1)*M2%DZ(1)
!DVRAT = DV2/DV1
!
!INDEX_LIST = (/I_LO,I_HI,J_LO,J_HI,K_LO,K_HI,II_LO,JJ_LO,KK_LO,NRX,NRY,NRZ/)
!VOLUME_LIST = (/DV1,DV2,DVRAT/)
!
!END SUBROUTINE LOCATE_MESH
!
!END MODULE EMBEDDED_MESH_METHOD


MODULE COMPLEX_GEOMETRY

USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_VARIABLES
USE MESH_POINTERS
USE TURBULENCE

IMPLICIT NONE

PRIVATE
PUBLIC :: INIT_IBM,VELTAN2D,VELTAN3D,TRILINEAR,GETX,GETU,GETGRAD
 
CONTAINS

SUBROUTINE INIT_IBM(T,NM)
USE MEMORY_FUNCTIONS, ONLY: ChkMemErr
IMPLICIT NONE
INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: T
INTEGER :: IZERO,I,J,K,NG
TYPE (MESH_TYPE), POINTER :: M
TYPE (GEOMETRY_TYPE), POINTER :: G
REAL(EB) :: DELTA,D2,R2,RP2,XU(3),PP(3),DP,TIME,TOL=1.E-9_EB
REAL(EB) :: X_MIN,Y_MIN,Z_MIN,X_MAX,Y_MAX,Z_MAX

TIME = T
M => MESHES(NM)

! geometry loop

GEOM_LOOP: DO NG=1,N_GEOM

   G => GEOMETRY(NG)

   IF (ICYC>1 .AND. (.NOT. G%TRANSLATE) .AND. (.NOT. G%ROTATE)) CYCLE GEOM_LOOP
   
   ! acceleration (not implemented yet)
   
   G%U = G%U0
   G%V = G%V0
   G%W = G%W0
   
   ! translation (linear for now)
   
   G%X = G%X0 + G%U*TIME
   G%Y = G%Y0 + G%V*TIME
   G%Z = G%Z0 + G%W*TIME
   
   IBM_UVWMAX = MAXVAL((/ABS(G%U),ABS(G%V),ABS(G%W)/))*RDX(1)
      
   IF (TWO_D) THEN
      !DELTA = SQRT(M%DX(1)*M%DZ(1))
      DELTA = MIN(M%DX(1),M%DZ(1))
   ELSE
      !DELTA = (M%DX(1)*M%DY(1)*M%DZ(1))**ONTH
      DELTA = MIN(M%DX(1),M%DY(1),M%DZ(1))
   ENDIF
   D2 = DELTA**2

   ! find bounding box

   SELECT_SHAPE: SELECT CASE(G%ISHAPE)
      CASE(IBOX) SELECT_SHAPE
         G%X1 = G%X1 + G%U*DT
         G%X2 = G%X2 + G%U*DT
         G%Y1 = G%Y1 + G%V*DT
         G%Y2 = G%Y2 + G%V*DT
         G%Z1 = G%Z1 + G%W*DT
         G%Z2 = G%Z2 + G%W*DT
         X_MIN = G%X1
         X_MAX = G%X2
         Y_MIN = G%Y1
         Y_MAX = G%Y2
         Z_MIN = G%Z1
         Z_MAX = G%Z2
         G%HL(1) = 0.5_EB*(X_MAX-X_MIN) + TOL
         G%HL(2) = 0.5_EB*(X_MAX-X_MIN) + TOL
         G%HL(3) = 0.5_EB*(X_MAX-X_MIN) + TOL
      CASE(ISPHERE) SELECT_SHAPE
         X_MIN = G%X-G%RADIUS
         Y_MIN = G%Y-G%RADIUS
         Z_MIN = G%Z-G%RADIUS
         X_MAX = G%X+G%RADIUS
         Y_MAX = G%Y+G%RADIUS
         Z_MAX = G%Z+G%RADIUS
         R2 = G%RADIUS**2
         IBM_UVWMAX = IBM_UVWMAX + G%RADIUS*MAXVAL((/ABS(G%OMEGA_X),ABS(G%OMEGA_Y),ABS(G%OMEGA_Z)/))*RDX(1)
      CASE(ICYLINDER) SELECT_SHAPE ! assume aligned with y axis
         G%HL(1) = ABS(G%XOR-G%X)
         G%HL(2) = ABS(G%YOR-G%Y)
         G%HL(3) = ABS(G%ZOR-G%Z)
         X_MIN = G%X-G%RADIUS
         Y_MIN = G%Y-G%HL(2)
         Z_MIN = G%Z-G%RADIUS
         X_MAX = G%X+G%RADIUS
         Y_MAX = G%Y+G%HL(2)
         Z_MAX = G%Z+G%RADIUS
         R2  = G%RADIUS**2
      CASE(IPLANE)
         X_MIN = M%XS
         Y_MIN = M%YS
         Z_MIN = M%ZS
         X_MAX = M%XF
         Y_MAX = M%YF
         Z_MAX = M%ZF
         PP   = (/G%X,G%Y,G%Z/)
         G%NN = (/G%XOR,G%YOR,G%ZOR/) - PP          ! normal vector to plane
         G%NN = G%NN/SQRT(DOT_PRODUCT(G%NN,G%NN))   ! unit normal
   END SELECT SELECT_SHAPE

   G%MIN_I(NM) = M%IBAR
   G%MIN_J(NM) = M%JBAR
   G%MIN_K(NM) = M%KBAR

   IF (X_MIN>=M%XS .AND. X_MIN<=M%XF) G%MIN_I(NM) = MAX(0,FLOOR((X_MIN-M%XS)/M%DX(1))-1)
   IF (Y_MIN>=M%YS .AND. Y_MIN<=M%YF) G%MIN_J(NM) = MAX(0,FLOOR((Y_MIN-M%YS)/M%DY(1))-1)
   IF (Z_MIN>=M%ZS .AND. Z_MIN<=M%ZF) G%MIN_K(NM) = MAX(0,FLOOR((Z_MIN-M%ZS)/M%DZ(1))-1)
   
   G%MAX_I(NM) = 0
   G%MAX_J(NM) = 0
   G%MAX_K(NM) = 0

   IF (X_MAX>=M%XS .AND. X_MAX<=M%XF) G%MAX_I(NM) = MIN(M%IBAR,CEILING((X_MAX-M%XS)/M%DX(1))+1)
   IF (Y_MAX>=M%YS .AND. Y_MAX<=M%YF) G%MAX_J(NM) = MIN(M%JBAR,CEILING((Y_MAX-M%YS)/M%DY(1))+1)
   IF (Z_MAX>=M%ZS .AND. Z_MAX<=M%ZF) G%MAX_K(NM) = MIN(M%KBAR,CEILING((Z_MAX-M%ZS)/M%DZ(1))+1)
   
   IF (TWO_D) THEN
      G%MIN_J(NM)=1
      G%MAX_J(NM)=1
   ENDIF
   
   IF ( G%MAX_I(NM)<G%MIN_I(NM) .OR. &
        G%MAX_J(NM)<G%MIN_J(NM) .OR. &
        G%MAX_K(NM)<G%MIN_K(NM) ) CYCLE GEOM_LOOP
   
   IF (G%IBM_ALLOCATED) DEALLOCATE(G%U_MASK,G%V_MASK,G%W_MASK)
   G%IBM_ALLOCATED=.FALSE.
   
   ALLOCATE(G%U_MASK(G%MIN_I(NM):G%MAX_I(NM),G%MIN_J(NM):G%MAX_J(NM),G%MIN_K(NM):G%MAX_K(NM)),STAT=IZERO)
   CALL ChkMemErr('INIT_IBM','U_MASK',IZERO)
   ALLOCATE(G%V_MASK(G%MIN_I(NM):G%MAX_I(NM),G%MIN_J(NM):G%MAX_J(NM),G%MIN_K(NM):G%MAX_K(NM)),STAT=IZERO)
   CALL ChkMemErr('INIT_IBM','V_MASK',IZERO)
   ALLOCATE(G%W_MASK(G%MIN_I(NM):G%MAX_I(NM),G%MIN_J(NM):G%MAX_J(NM),G%MIN_K(NM):G%MAX_K(NM)),STAT=IZERO)
   CALL ChkMemErr('INIT_IBM','W_MASK',IZERO)
   IF (IZERO==0) G%IBM_ALLOCATED=.TRUE.
   
   G%U_MASK = 1 ! default to gas phase
   G%V_MASK = 1
   G%W_MASK = 1
   
   ! mask cells

   DO K=G%MIN_K(NM),G%MAX_K(NM)
      DO J=G%MIN_J(NM),G%MAX_J(NM)
         DO I=G%MIN_I(NM),G%MAX_I(NM)
         
            MASK_SHAPE: SELECT CASE(G%ISHAPE)
            
               CASE(IBOX) MASK_SHAPE
                  
                  ! see if the point is inside geometry
                  IF (ABS( M%X(I)-G%X)<G%HL(1) .AND. &
                      ABS(M%YC(J)-G%Y)<G%HL(2) .AND. &
                      ABS(M%ZC(K)-G%Z)<G%HL(3)) G%U_MASK(I,J,K) = -1
                  
                  IF (ABS(M%XC(I)-G%X)<G%HL(1) .AND. &
                      ABS( M%Y(J)-G%Y)<G%HL(2) .AND. &
                      ABS(M%ZC(K)-G%Z)<G%HL(3)) G%V_MASK(I,J,K) = -1
                  
                  IF (ABS(M%XC(I)-G%X)<G%HL(1) .AND. &
                      ABS(M%YC(J)-G%Y)<G%HL(2) .AND. &
                      ABS( M%Z(K)-G%Z)<G%HL(3)) G%W_MASK(I,J,K) = -1
                  
                  ! see if the point is in surface layer
                  IF (X_MAX<M%X(I) .AND. M%X(I)<X_MAX+DELTA) G%U_MASK(I,J,K) = 0
                  IF (Y_MAX<M%Y(J) .AND. M%Y(J)<Y_MAX+DELTA) G%V_MASK(I,J,K) = 0
                  IF (Z_MAX<M%Z(K) .AND. M%Z(K)<Z_MAX+DELTA) G%W_MASK(I,J,K) = 0
                  
                  IF (X_MIN-DELTA<M%X(I) .AND. M%X(I)<X_MIN) G%U_MASK(I,J,K) = 0
                  IF (Y_MIN-DELTA<M%Y(J) .AND. M%Y(J)<Y_MIN) G%V_MASK(I,J,K) = 0
                  IF (Z_MIN-DELTA<M%Z(K) .AND. M%Z(K)<Z_MIN) G%W_MASK(I,J,K) = 0
                  
               CASE(ISPHERE) MASK_SHAPE
               
                  RP2 = (M%X(I)-G%X)**2+(M%YC(J)-G%Y)**2+(M%ZC(K)-G%Z)**2
                  IF (RP2-R2 < D2 ) G%U_MASK(I,J,K) = 0
                  IF (RP2-R2 < TOL) G%U_MASK(I,J,K) = -1
                  
                  RP2 = (M%XC(I)-G%X)**2+(M%Y(J)-G%Y)**2+(M%ZC(K)-G%Z)**2
                  IF (RP2-R2 < D2 ) G%V_MASK(I,J,K) = 0
                  IF (RP2-R2 < TOL) G%V_MASK(I,J,K) = -1
                  
                  RP2 = (M%XC(I)-G%X)**2+(M%YC(J)-G%Y)**2+(M%Z(K)-G%Z)**2
                  IF (RP2-R2 < D2 ) G%W_MASK(I,J,K) = 0
                  IF (RP2-R2 < TOL) G%W_MASK(I,J,K) = -1
                  
               CASE(ICYLINDER) MASK_SHAPE ! align with y axis for now
               
                  RP2 = (M%X(I)-G%X)**2+(M%ZC(K)-G%Z)**2
                  IF (RP2-R2 < D2 ) G%U_MASK(I,J,K) = 0
                  IF (RP2-R2 < TOL) G%U_MASK(I,J,K) = -1
                  
                  RP2 = (M%XC(I)-G%X)**2+(M%ZC(K)-G%Z)**2
                  IF (RP2-R2 < D2 ) G%V_MASK(I,J,K) = 0
                  IF (RP2-R2 < TOL) G%V_MASK(I,J,K) = -1
                  
                  RP2 = (M%XC(I)-G%X)**2+(M%Z(K)-G%Z)**2
                  IF (RP2-R2 < D2 ) G%W_MASK(I,J,K) = 0
                  IF (RP2-R2 < TOL) G%W_MASK(I,J,K) = -1
                  
               CASE(IPLANE) MASK_SHAPE
               
                  ! see Section 10.3 Schneider and Eberly
                  
                  XU = (/M%X(I),M%YC(J),M%ZC(K)/)        ! point of interest
                  IF (G%TWO_SIDED) THEN
                     DP = ABS(DOT_PRODUCT(G%NN,XU-PP))   ! distance to plane
                  ELSE
                     DP = DOT_PRODUCT(G%NN,XU-PP)        ! signed distance to plane
                  ENDIF
                  IF (DP<DELTA) G%U_MASK(I,J,K) = 0
                  IF (DP<TOL)   G%U_MASK(I,J,K) = -1
                  
                  XU = (/M%XC(I),M%Y(J),M%ZC(K)/)
                  IF (G%TWO_SIDED) THEN
                     DP = ABS(DOT_PRODUCT(G%NN,XU-PP))
                  ELSE
                     DP = DOT_PRODUCT(G%NN,XU-PP)
                  ENDIF
                  IF (DP<DELTA) G%V_MASK(I,J,K) = 0
                  IF (DP<TOL)   G%V_MASK(I,J,K) = -1
                  
                  XU = (/M%XC(I),M%YC(J),M%Z(K)/)
                  IF (G%TWO_SIDED) THEN
                     DP = ABS(DOT_PRODUCT(G%NN,XU-PP))
                  ELSE
                     DP = DOT_PRODUCT(G%NN,XU-PP)
                  ENDIF
                  IF (DP<DELTA) G%W_MASK(I,J,K) = 0
                  IF (DP<TOL)   G%W_MASK(I,J,K) = -1
                  
            END SELECT MASK_SHAPE
      
         ENDDO
      ENDDO
   ENDDO

ENDDO GEOM_LOOP

END SUBROUTINE INIT_IBM


REAL(EB) FUNCTION TRILINEAR(UU,DXI,LL)
IMPLICIT NONE

REAL(EB), INTENT(IN) :: UU(0:1,0:1,0:1),DXI(3),LL(3)

! Comments:
!
! see http://local.wasp.uwa.edu.au/~pbourke/miscellaneous/interpolation/index.html
! with appropriate scaling. LL is length of side.
!
!                       UU(1,1,1)
!        z /----------/
!        ^/          / |
!        ------------  |    Particle position
!        |          |  |
!  LL(3) |   o<-----|------- DXI = [DXI(1),DXI(2),DXI(3)]
!        |          | /        
!        |          |/      Particle property at XX = TRILINEAR
!        ------------> x
!        ^
!        |
!   X0 = [0,0,0]
!
!    UU(0,0,0)
!
!===========================================================

TRILINEAR = UU(0,0,0)*(LL(1)-DXI(1))*(LL(2)-DXI(2))*(LL(3)-DXI(3)) +    &
            UU(1,0,0)*DXI(1)*(LL(2)-DXI(2))*(LL(3)-DXI(3)) +            &
            UU(0,1,0)*(LL(1)-DXI(1))*DXI(2)*(LL(3)-DXI(3)) +            &
            UU(0,0,1)*(LL(1)-DXI(1))*(LL(2)-DXI(2))*DXI(3) +            &
            UU(1,0,1)*DXI(1)*(LL(2)-DXI(2))*DXI(3) +                    &
            UU(0,1,1)*(LL(1)-DXI(1))*DXI(2)*DXI(3) +                    &
            UU(1,1,0)*DXI(1)*DXI(2)*(LL(3)-DXI(3)) +                    &
            UU(1,1,1)*DXI(1)*DXI(2)*DXI(3)

TRILINEAR = TRILINEAR/(LL(1)*LL(2)*LL(3))

END FUNCTION TRILINEAR


SUBROUTINE GETX(XI,XU,NG)
IMPLICIT NONE

REAL(EB), INTENT(OUT) :: XI(3)
REAL(EB), INTENT(IN) :: XU(3)
INTEGER, INTENT(IN) :: NG
TYPE(GEOMETRY_TYPE), POINTER :: G
REAL(EB) :: PP(3),RU(3),RUMAG,DR,DP

G => GEOMETRY(NG)
SELECT CASE(G%ISHAPE)
   CASE(IBOX)
      XI = XU
      IF (XU(1)<G%X1) XI(1) = XU(1) + (XU(1)-G%X1)
      IF (XU(1)>G%X2) XI(1) = XU(1) + (XU(1)-G%X2)
      IF (XU(2)<G%Y1) XI(2) = XU(2) + (XU(2)-G%Y1)
      IF (XU(2)>G%Y2) XI(2) = XU(2) + (XU(2)-G%Y2)
      IF (XU(3)<G%Z1) XI(3) = XU(3) + (XU(3)-G%Z1)
      IF (XU(3)>G%Z2) XI(3) = XU(3) + (XU(3)-G%Z2)
   CASE(ISPHERE)
      RU     = XU-(/G%X,G%Y,G%Z/)
      RUMAG  = SQRT(DOT_PRODUCT(RU,RU))
      DR     = RUMAG-G%RADIUS
      XI     = XU + DR*RU/RUMAG
   CASE(ICYLINDER)
      RU     = (/XU(1),0._EB,XU(3)/)-(/G%X,0._EB,G%Z/)
      RUMAG  = SQRT(DOT_PRODUCT(RU,RU))
      DR     = RUMAG-G%RADIUS
      XI     = XU + DR*RU/RUMAG
   CASE(IPLANE)
      PP = (/G%X,G%Y,G%Z/)           ! point in the plane
      DP = DOT_PRODUCT(G%NN,XU-PP)   ! signed distance to plane
      XI = XU + DP*G%NN
END SELECT

END SUBROUTINE GETX


SUBROUTINE GETU(U_DATA,DXI,XI,XU,INDU,I_VEL,NM)
IMPLICIT NONE

REAL(EB), INTENT(OUT) :: U_DATA(0:1,0:1,0:1),DXI(3)
REAL(EB), INTENT(IN) :: XI(3),XU(3)
INTEGER, INTENT(IN) :: INDU(3),I_VEL,NM
TYPE(MESH_TYPE), POINTER :: M
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW
INTEGER :: II,JJ,KK

M=>MESHES(NM)
IF (PREDICTOR) THEN
   UU => M%U
   VV => M%V
   WW => M%W
ELSE
   UU => M%US
   VV => M%VS
   WW => M%WS
ENDIF

! first assume XI >= XU
II = INDU(1)
JJ = INDU(2)
KK = INDU(3)
DXI = XI-XU

SELECT CASE(I_VEL)
   CASE(1)
      IF (XI(1)<XU(1)) THEN
         II=MAX(0,II-1)
         DXI(1)=DXI(1)+M%DX(II)
      ENDIF
      IF (XI(2)<XU(2)) THEN
         JJ=MAX(0,JJ-1)
         DXI(2)=DXI(2)+M%DYN(JJ)
      ENDIF
      IF (XI(3)<XU(3)) THEN
         KK=MAX(0,KK-1)
         DXI(3)=DXI(3)+M%DZN(KK)
      ENDIF
      U_DATA(0,0,0) = UU(II,JJ,KK)
      U_DATA(1,0,0) = UU(II+1,JJ,KK)
      U_DATA(0,1,0) = UU(II,JJ+1,KK)
      U_DATA(0,0,1) = UU(II,JJ,KK+1)
      U_DATA(1,0,1) = UU(II+1,JJ,KK+1)
      U_DATA(0,1,1) = UU(II,JJ+1,KK+1)
      U_DATA(1,1,0) = UU(II+1,JJ+1,KK)
      U_DATA(1,1,1) = UU(II+1,JJ+1,KK+1)
   CASE(2)
      IF (XI(1)<XU(1)) THEN
         II=MAX(0,II-1)
         DXI(1)=DXI(1)+M%DXN(II)
      ENDIF
      IF (XI(2)<XU(2)) THEN
         JJ=MAX(0,JJ-1)
         DXI(2)=DXI(2)+M%DY(JJ)
      ENDIF
      IF (XI(3)<XU(3)) THEN
         KK=MAX(0,KK-1)
         DXI(3)=DXI(3)+M%DZN(KK)
      ENDIF
      U_DATA(0,0,0) = VV(II,JJ,KK)
      U_DATA(1,0,0) = VV(II+1,JJ,KK)
      U_DATA(0,1,0) = VV(II,JJ+1,KK)
      U_DATA(0,0,1) = VV(II,JJ,KK+1)
      U_DATA(1,0,1) = VV(II+1,JJ,KK+1)
      U_DATA(0,1,1) = VV(II,JJ+1,KK+1)
      U_DATA(1,1,0) = VV(II+1,JJ+1,KK)
      U_DATA(1,1,1) = VV(II+1,JJ+1,KK+1)
   CASE(3)
      IF (XI(1)<XU(1)) THEN
         II=MAX(0,II-1)
         DXI(1)=DXI(1)+M%DXN(II)
      ENDIF
      IF (XI(2)<XU(2)) THEN
         JJ=MAX(0,JJ-1)
         DXI(2)=DXI(2)+M%DYN(JJ)
      ENDIF
      IF (XI(3)<XU(3)) THEN
         KK=MAX(0,KK-1)
         DXI(3)=DXI(3)+M%DZ(KK)
      ENDIF
      U_DATA(0,0,0) = WW(II,JJ,KK)
      U_DATA(1,0,0) = WW(II+1,JJ,KK)
      U_DATA(0,1,0) = WW(II,JJ+1,KK)
      U_DATA(0,0,1) = WW(II,JJ,KK+1)
      U_DATA(1,0,1) = WW(II+1,JJ,KK+1)
      U_DATA(0,1,1) = WW(II,JJ+1,KK+1)
      U_DATA(1,1,0) = WW(II+1,JJ+1,KK)
      U_DATA(1,1,1) = WW(II+1,JJ+1,KK+1)
   CASE(4) ! viscosity
      IF (XI(1)<XU(1)) THEN
         II=MAX(0,II-1)
         DXI(1)=DXI(1)+M%DXN(II)
      ENDIF
      IF (XI(2)<XU(2)) THEN
         JJ=MAX(0,JJ-1)
         DXI(2)=DXI(2)+M%DYN(JJ)
      ENDIF
      IF (XI(3)<XU(3)) THEN
         KK=MAX(0,KK-1)
         DXI(3)=DXI(3)+M%DZN(KK)
      ENDIF
      U_DATA(0,0,0) = M%MU(II,JJ,KK)
      U_DATA(1,0,0) = M%MU(II+1,JJ,KK)
      U_DATA(0,1,0) = M%MU(II,JJ+1,KK)
      U_DATA(0,0,1) = M%MU(II,JJ,KK+1)
      U_DATA(1,0,1) = M%MU(II+1,JJ,KK+1)
      U_DATA(0,1,1) = M%MU(II,JJ+1,KK+1)
      U_DATA(1,1,0) = M%MU(II+1,JJ+1,KK)
      U_DATA(1,1,1) = M%MU(II+1,JJ+1,KK+1)
END SELECT

END SUBROUTINE GETU


SUBROUTINE GETGRAD(G_DATA,DXI,XI,XU,INDU,COMP_I,COMP_J,NM)
IMPLICIT NONE

REAL(EB), INTENT(OUT) :: G_DATA(0:1,0:1,0:1),DXI(3)
REAL(EB), INTENT(IN) :: XI(3),XU(3)
INTEGER, INTENT(IN) :: INDU(3),COMP_I,COMP_J,NM
TYPE(MESH_TYPE), POINTER :: M
REAL(EB), POINTER, DIMENSION(:,:,:) :: DUDX
INTEGER :: II,JJ,KK

M=>MESHES(NM)

IF (COMP_I==1 .AND. COMP_J==1) DUDX => M%WORK5
IF (COMP_I==1 .AND. COMP_J==2) DUDX => M%IBM_SAVE1
IF (COMP_I==1 .AND. COMP_J==3) DUDX => M%IBM_SAVE2
IF (COMP_I==2 .AND. COMP_J==1) DUDX => M%IBM_SAVE3
IF (COMP_I==2 .AND. COMP_J==2) DUDX => M%WORK6
IF (COMP_I==2 .AND. COMP_J==3) DUDX => M%IBM_SAVE4
IF (COMP_I==3 .AND. COMP_J==1) DUDX => M%IBM_SAVE5
IF (COMP_I==3 .AND. COMP_J==2) DUDX => M%IBM_SAVE6
IF (COMP_I==3 .AND. COMP_J==3) DUDX => M%WORK7

! first assume XI >= XU
II = INDU(1)
JJ = INDU(2)
KK = INDU(3)
DXI = XI-XU

IF (XI(1)<XU(1)) THEN
   II=MAX(0,II-1)
   DXI(1)=DXI(1)+M%DX(II)
ENDIF
IF (XI(2)<XU(2)) THEN
   JJ=MAX(0,JJ-1)
   DXI(2)=DXI(2)+M%DY(JJ)
ENDIF
IF (XI(3)<XU(3)) THEN
   KK=MAX(0,KK-1)
   DXI(3)=DXI(3)+M%DZ(KK)
ENDIF
G_DATA(0,0,0) = DUDX(II,JJ,KK)
G_DATA(1,0,0) = DUDX(II+1,JJ,KK)
G_DATA(0,1,0) = DUDX(II,JJ+1,KK)
G_DATA(0,0,1) = DUDX(II,JJ,KK+1)
G_DATA(1,0,1) = DUDX(II+1,JJ,KK+1)
G_DATA(0,1,1) = DUDX(II,JJ+1,KK+1)
G_DATA(1,1,0) = DUDX(II+1,JJ+1,KK)
G_DATA(1,1,1) = DUDX(II+1,JJ+1,KK+1)

END SUBROUTINE GETGRAD


REAL(EB) FUNCTION VELTAN2D(U_VELO,U_SURF,NN,DN,DIVU,GRADU,GRADP,TAU_IJ,DT,RRHO,MU,I_VEL)
IMPLICIT NONE

REAL(EB), INTENT(IN) :: U_VELO(2),U_SURF(2),NN(2),DN,DIVU,GRADU(2,2),GRADP(2),TAU_IJ(2,2),DT,RRHO,MU
INTEGER, INTENT(IN) :: I_VEL
REAL(EB) :: C(2,2),SS(2),SLIP_COEF,ETA,AA,BB,U_STRM_0,DUMMY, &
            U_STRM,U_NORM,U_STRM_WALL,U_NORM_WALL,DPDS,DUSDS,DUSDN,TSN
INTEGER :: SUBIT

! Cartesian grid coordinate system orthonormal basis vectors
REAL(EB), DIMENSION(2), PARAMETER :: XX=(/1._EB, 0._EB/),YY=(/0._EB, 1._EB/)


! streamwise unit vector
SS = (/NN(2),-NN(1)/)

! directional cosines (see Pope, Eq. A.11)
C(1,1) = DOT_PRODUCT(XX,SS)
C(1,2) = DOT_PRODUCT(XX,NN)
C(2,1) = DOT_PRODUCT(YY,SS)
C(2,2) = DOT_PRODUCT(YY,NN)

! transform velocity (see Pope, Eq. A.17)
U_STRM = C(1,1)*U_VELO(1) + C(2,1)*U_VELO(2)
U_NORM = C(1,2)*U_VELO(1) + C(2,2)*U_VELO(2)

! transform wall velocity
U_STRM_WALL = C(1,1)*U_SURF(1) + C(2,1)*U_SURF(2)
U_NORM_WALL = C(1,2)*U_SURF(1) + C(2,2)*U_SURF(2)

! transform pressure gradient
DPDS = C(1,1)*GRADP(1) + C(2,1)*GRADP(2)

! transform velocity gradient tensor (Pope A.23)
DUSDS = C(1,1)*C(1,1)*GRADU(1,1) + C(1,1)*C(2,1)*GRADU(1,2) &
      + C(2,1)*C(1,1)*GRADU(2,1) + C(2,1)*C(2,1)*GRADU(2,2)
      
DUSDN = C(1,1)*C(1,2)*GRADU(1,1) + C(1,1)*C(2,2)*GRADU(1,2) &
      + C(2,1)*C(1,2)*GRADU(2,1) + C(2,1)*C(2,2)*GRADU(2,2)

! transform stress tensor
TSN = C(1,1)*C(1,2)*TAU_IJ(1,1) + C(1,1)*C(2,2)*TAU_IJ(1,2) &
    + C(2,1)*C(1,2)*TAU_IJ(2,1) + C(2,1)*C(2,2)*TAU_IJ(2,2)
    
! update boundary layer equations

! update wall-normal velocity
U_NORM = U_NORM_WALL + DN*(DIVU-0.5_EB*DUSDS)

! ODE solution
IF (DNS) THEN
   ETA = U_NORM + RRHO*MU/DN
   AA  = -(0.5_EB*DUSDS + TWTH*ETA/DN)
   BB  = (TWTH*U_STRM_WALL/DN + ONSI*DUSDN)*ETA - (U_NORM*0.5_EB*DUSDN + RRHO*( DPDS + TSN/(2._EB*DN) ))
   !AA  = -0.5_EB*(DUSDS + ETA/DN)
   !BB  = 0.5_EB*US_WALL/DN*ETA - (UN*0.5_EB*DUSDN + RRHO*( DPDS + TSN/(2._EB*DN) ))
   U_STRM = ((AA*U_STRM + BB)*EXP(AA*DT) - BB)/AA
ELSE
   U_STRM_0 = U_STRM
   DO SUBIT=1,1
      CALL WERNER_WENGLE_WALL_MODEL(SLIP_COEF,DUMMY,U_STRM-U_STRM_WALL,MU*RRHO,DN,0._EB)
      !IF (SLIP_COEF< -1._EB .OR. SLIP_COEF>-1._EB) THEN
      !   PRINT *,SUBIT,'WARNING: SLIP_COEF=',SLIP_COEF
      !ENDIF
      ETA = RRHO*(1-SLIP_COEF)*MU/(2._EB*DN**2)
      AA  = -(0.5_EB*DUSDS + TWTH*U_NORM/DN + ETA)
      BB  = ETA*U_STRM_WALL - (U_NORM*ONTH*DUSDN + RRHO*( DPDS + TSN/(2._EB*DN) ))
      U_STRM = ((AA*U_STRM_0 + BB)*EXP(AA*DT) - BB)/AA
   ENDDO
ENDIF

! transform velocity back to Cartesian component I_VEL
VELTAN2D = C(I_VEL,1)*U_STRM + C(I_VEL,2)*U_NORM

END FUNCTION VELTAN2D


REAL(EB) FUNCTION VELTAN3D(U_VELO,U_SURF,NN,DN,DIVU,GRADU,GRADP,TAU_IJ,DT,RRHO,MU,I_VEL,ROUGHNESS)
USE MATH_FUNCTIONS, ONLY: CROSS_PRODUCT, NORM2
IMPLICIT NONE

REAL(EB), INTENT(IN) :: U_VELO(3),U_SURF(3),NN(3),DN,DIVU,GRADU(3,3),GRADP(3),TAU_IJ(3,3),DT,RRHO,MU,ROUGHNESS
INTEGER, INTENT(IN) :: I_VEL
REAL(EB) :: C(3,3),SS(3),PP(3),SLIP_COEF,ETA,AA,BB,U_STRM_0,DUMMY,U_RELA(3), &
            U_STRM,U_ORTH,U_NORM,DPDS,DUSDS,DUSDN,TSN
INTEGER :: SUBIT,I,J

! Cartesian grid coordinate system orthonormal basis vectors
REAL(EB), DIMENSION(3), PARAMETER :: E1=(/1._EB,0._EB,0._EB/),E2=(/0._EB,1._EB,0._EB/),E3=(/0._EB,0._EB,1._EB/)

U_RELA = U_VELO-U_SURF
IF (DOT_PRODUCT(NN,U_RELA)<1.E-6_EB) THEN
   VELTAN3D = U_VELO(I_VEL)
   RETURN
ENDIF

! find a vector PP in the tangent plane of the surface and orthogonal to U_VELO-U_SURF
CALL CROSS_PRODUCT(PP,NN,U_RELA)
PP = PP/NORM2(PP) ! normalize to unit vector

! define the streamwise unit vector SS
CALL CROSS_PRODUCT(SS,PP,NN)

! directional cosines (see Pope, Eq. A.11)
C(1,1) = DOT_PRODUCT(E1,SS)
C(1,2) = DOT_PRODUCT(E1,PP)
C(1,3) = DOT_PRODUCT(E1,NN)
C(2,1) = DOT_PRODUCT(E2,SS)
C(2,2) = DOT_PRODUCT(E2,PP)
C(2,3) = DOT_PRODUCT(E2,NN)
C(3,1) = DOT_PRODUCT(E3,SS)
C(3,2) = DOT_PRODUCT(E3,PP)
C(3,3) = DOT_PRODUCT(E3,NN)

! transform velocity (see Pope, Eq. A.17)
U_STRM = C(1,1)*U_RELA(1) + C(2,1)*U_RELA(2) + C(3,1)*U_RELA(3)
U_ORTH = C(1,2)*U_RELA(1) + C(2,2)*U_RELA(2) + C(3,2)*U_RELA(3)
U_NORM = C(1,3)*U_RELA(1) + C(2,3)*U_RELA(2) + C(3,3)*U_RELA(3)

!! check UP, should be zero
!print *, U_ORTH

! transform pressure gradient
DPDS = C(1,1)*GRADP(1) + C(2,1)*GRADP(2) + C(3,1)*GRADP(3)

! transform tensors (Pope A.23)
DUSDS = 0._EB
DUSDN = 0._EB
TSN = 0._EB
DO I=1,3
   DO J=1,3
      DUSDS = DUSDS + C(I,1)*C(J,1)*GRADU(I,J)
      DUSDN = DUSDN + C(I,1)*C(J,3)*GRADU(I,J)
      TSN = TSN + C(I,1)*C(J,3)*TAU_IJ(I,J)
   ENDDO
ENDDO
    
! update boundary layer equations

! update wall-normal velocity
U_NORM = DN*(DIVU-0.5_EB*DUSDS)

! ODE solution
IF (DNS) THEN
   ETA = U_NORM + RRHO*MU/DN
   AA  = -(0.5_EB*DUSDS + TWTH*ETA/DN)
   BB  = ONSI*DUSDN*ETA - (U_NORM*0.5_EB*DUSDN + RRHO*( DPDS + TSN/(2._EB*DN) ))
   !AA  = -0.5_EB*(DUSDS + ETA/DN)
   !BB  = 0.5_EB*US_WALL/DN*ETA - (UN*0.5_EB*DUSDN + RRHO*( DPDS + TSN/(2._EB*DN) ))
   U_STRM = ((AA*U_STRM + BB)*EXP(AA*DT) - BB)/AA
ELSE
   U_STRM_0 = U_STRM
   DO SUBIT=1,1
      CALL WERNER_WENGLE_WALL_MODEL(SLIP_COEF,DUMMY,U_STRM,MU*RRHO,DN,ROUGHNESS)
      !IF (SLIP_COEF< -1._EB .OR. SLIP_COEF>-1._EB) THEN
      !   PRINT *,SUBIT,'WARNING: SLIP_COEF=',SLIP_COEF
      !ENDIF
      ETA = RRHO*(1-SLIP_COEF)*MU/(2._EB*DN**2)
      AA  = -(0.5_EB*DUSDS + TWTH*U_NORM/DN + ETA)
      BB  = -(U_NORM*ONTH*DUSDN + RRHO*( DPDS + TSN/(2._EB*DN) ))
      U_STRM = ((AA*U_STRM_0 + BB)*EXP(AA*DT) - BB)/AA
   ENDDO
ENDIF

! transform velocity back to Cartesian component I_VEL
VELTAN3D = C(I_VEL,1)*U_STRM + C(I_VEL,3)*U_NORM + U_SURF(I_VEL)

END FUNCTION VELTAN3D



END MODULE COMPLEX_GEOMETRY

