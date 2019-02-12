! This module is useful for verification tests and development of
! turbulence models.

MODULE TURBULENCE

USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_POINTERS
USE MESH_VARIABLES
USE COMP_FUNCTIONS

IMPLICIT NONE
PRIVATE

PUBLIC :: INIT_TURB_ARRAYS, VARDEN_DYNSMAG, WANNIER_FLOW, &
          WALL_MODEL, COMPRESSION_WAVE, VELTAN2D,VELTAN3D, &
          SYNTHETIC_TURBULENCE, SYNTHETIC_EDDY_SETUP, TEST_FILTER, EX2G3D, TENSOR_DIFFUSIVITY_MODEL, &
          TWOD_VORTEX_CERFACS, TWOD_VORTEX_UMD, TWOD_SOBOROT_UMD, &
          LOGLAW_HEAT_FLUX_MODEL, ABL_HEAT_FLUX_MODEL, RNG_EDDY_VISCOSITY, &
          NS_ANALYTICAL_SOLUTION, NS_U_EXACT, NS_V_EXACT, NS_H_EXACT, SANDIA_DAT, SPECTRAL_OUTPUT, SANDIA_OUT, &
          FILL_EDGES, NATURAL_CONVECTION_MODEL, FORCED_CONVECTION_MODEL, RAYLEIGH_HEAT_FLUX_MODEL, YUAN_HEAT_FLUX_MODEL, &
          WALE_VISCOSITY, TAU_WALL_IJ

CONTAINS


SUBROUTINE INIT_TURB_ARRAYS(NM)

USE MEMORY_FUNCTIONS, ONLY: ChkMemErr
INTEGER, INTENT(IN) :: NM
INTEGER :: IZERO
TYPE (MESH_TYPE), POINTER :: M

IF (EVACUATION_ONLY(NM)) RETURN
CALL POINT_TO_MESH(NM)
M => MESHES(NM)

!IF (PERIODIC_TEST==2 .OR. TURB_MODEL==DYNSMAG) THEN
!   ! real work arrays
!   ALLOCATE(M%TURB_WORK1(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
!   CALL ChkMemErr('INIT_TURB_ARRAYS','TURB_WORK1',IZERO)
!   M%TURB_WORK1 = 0._EB
!   ALLOCATE(M%TURB_WORK2(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
!   CALL ChkMemErr('INIT_TURB_ARRAYS','TURB_WORK2',IZERO)
!   M%TURB_WORK2 = 0._EB
!   ALLOCATE(M%TURB_WORK3(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
!   CALL ChkMemErr('INIT_TURB_ARRAYS','TURB_WORK3',IZERO)
!   M%TURB_WORK3 = 0._EB
!   ALLOCATE(M%TURB_WORK4(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
!   CALL ChkMemErr('INIT_TURB_ARRAYS','TURB_WORK4',IZERO)
!   M%TURB_WORK4 = 0._EB
!ENDIF

IF (TURB_MODEL==DYNSMAG) THEN
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

! Decaying isotropic turbulence (high resolution)

IF (PERIODIC_TEST==9) CALL INIT_SPECTRAL_DATA(NM)

END SUBROUTINE INIT_TURB_ARRAYS


SUBROUTINE NS_ANALYTICAL_SOLUTION(NM,T,RK_STAGE)

! Initialize flow variables with an analytical solution of the governing equations

INTEGER, INTENT(IN) :: NM,RK_STAGE
REAL(EB), INTENT(IN) :: T
INTEGER :: I,J,K
REAL(EB) :: UP,WP
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU=>NULL(),WW=>NULL(),RHOP=>NULL()
REAL(EB), PARAMETER :: AA=2._EB

IF (EVACUATION_ONLY(NM)) RETURN
CALL POINT_TO_MESH(NM)

SELECT CASE(RK_STAGE)
   CASE(1)
      UU=>US
      WW=>WS
      RHOP=>RHOS
   CASE(2)
      UU=>U
      WW=>W
      RHOP=>RHO
END SELECT

DO K=1,KBAR
   DO J=1,JBAR
      DO I=0,IBAR
         UU(I,J,K) = NS_U_EXACT(X(I),ZC(K),T,MU(I,J,K),RHOP(I,J,K),AA)
      ENDDO
   ENDDO
ENDDO
DO K=0,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         WW(I,J,K) = NS_V_EXACT(XC(I),Z(K),T,MU(I,J,K),RHOP(I,J,K),AA)
      ENDDO
   ENDDO
ENDDO
DO K=0,KBP1
   DO J=0,JBP1
      DO I=0,IBP1
         H(I,J,K)  = NS_H_EXACT(XC(I),ZC(K),T,MU(I,J,K),RHOP(I,J,K),AA)
         HS(I,J,K) = NS_H_EXACT(XC(I),ZC(K),T,MU(I,J,K),RHOP(I,J,K),AA)
      ENDDO
   ENDDO
ENDDO

! KRES loop

DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         UP = 0.5_EB*(UU(I,J,K) + UU(I-1,J,K))
         WP = 0.5_EB*(WW(I,J,K) + WW(I,J,K-1))
         KRES(I,J,K) = 0.5_EB*(UP**2 + WP**2)
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE NS_ANALYTICAL_SOLUTION


REAL(EB) FUNCTION NS_U_EXACT(XX,YY,TT,MUP,RHOP,AA)
REAL(EB), INTENT(IN) :: XX,YY,TT,MUP,RHOP,AA
NS_U_EXACT = 1._EB - AA*COS(XX - TT)*SIN(YY - TT)*EXP(-2._EB * MUP / RHOP * TT)
END FUNCTION NS_U_EXACT


REAL(EB) FUNCTION NS_V_EXACT(XX,YY,TT,MUP,RHOP,AA)
REAL(EB), INTENT(IN) :: XX,YY,TT,MUP,RHOP,AA
NS_V_EXACT = 1._EB + AA*SIN(XX - TT)*COS(YY - TT)*EXP(-2._EB * MUP / RHOP * TT)
END FUNCTION NS_V_EXACT


REAL(EB) FUNCTION NS_H_EXACT(XX,YY,TT,MUP,RHOP,AA)
REAL(EB), INTENT(IN) :: XX,YY,TT,MUP,RHOP,AA
REAL(EB) :: UP,VP
UP = NS_U_EXACT(XX,YY,TT,MUP,RHOP,AA)
VP = NS_V_EXACT(XX,YY,TT,MUP,RHOP,AA)
NS_H_EXACT = -AA**2/4._EB*( COS(2._EB*(XX-TT)) + COS(2._EB*(YY-TT)) )*EXP(-4._EB * MUP / RHOP * TT) + 0.5_EB*(UP**2 + VP**2)
END FUNCTION NS_H_EXACT


REAL(EB) FUNCTION WANNIER_FLOW(XX,YY,IVEL)

! References:
!
! G.H. Wannier. A contribution to the hydrodynamics of lubrication. Quart. Appl.
! Math. 8(1) (1950).
!
! S.J. Sherwin, G.E. Karniadakis. Comput. Methods Appl. Mech. Engrg. 123
! (1995) 189-229. (note: errors in this write up)
!
! T. Ye, R. Mittal, H.S. Udaykumar, and W. Shyy. An Accurate
! Cartesian Grid Method for Viscous Incompressible Flows with Complex Immersed
! Boundaries. J. Comp. Phys. 156:209-240 (1999). Appendix 2

REAL(EB), INTENT(IN) :: XX,YY ! position
INTEGER, INTENT(IN) :: IVEL ! velocity component
REAL(EB), PARAMETER :: RR=0.5_EB,HH=1.75_EB,UU=1._EB,X0=0._EB,Y0=HH
REAL(EB) :: AA,BB,CC,FF,GG,SS,OMEGA,K1,K2

IF ( (XX-X0)**2 + (YY-Y0)**2 < RR**2 ) THEN
   WANNIER_FLOW = 0._EB
   RETURN
ENDIF
OMEGA = -GEOMETRY(1)%OMEGA
SS = SQRT(HH**2-RR**2)
GG = (HH+SS)/(HH-SS)
K1 = XX**2+(SS+YY)**2
K2 = XX**2+(SS-YY)**2

FF = UU/LOG(GG)
AA = -HH*FF - 0.5_EB*HH*OMEGA/SS
BB = 2._EB*(HH+SS)*FF + (HH+SS)*OMEGA/SS
CC = 2._EB*(HH-SS)*FF + (HH-SS)*OMEGA/SS

SELECT CASE(IVEL)
   CASE(1)
      WANNIER_FLOW = UU - 2._EB*(AA+FF*YY)/K1*((SS+YY)+K1/K2*(SS-YY)) - FF*LOG(K1/K2) &
                   - BB/K1*(SS+2._EB*YY-2*YY*(SS+YY)**2._EB/K1) &
                   - CC/K2*(SS-2._EB*YY+2*YY*(SS-YY)**2._EB/K2)
   CASE(2)
      WANNIER_FLOW = 2._EB*XX/(K1*K2)*(AA+FF*YY)*(K2-K1) - 2._EB*BB*XX*YY*(SS+YY)/K1**2 - 2._EB*CC*XX*YY*(SS-YY)/K2**2
END SELECT

END FUNCTION WANNIER_FLOW


SUBROUTINE COMPRESSION_WAVE(NM,T,ITEST)

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


SUBROUTINE TWOD_VORTEX_CERFACS(NM)
!-------------------------------------------------------------------------------
! Ragini Acharya, Nov-Dec 2011, United Technologies Research Center
! Laterally moving 2d-vortex
! http://elearning.cerfacs.fr/pdfs/numerical/TestCaseVortex2D.pdf
!-------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER, INTENT(IN) :: NM
INTEGER :: I,J,K
REAL(EB) :: UINF,PINF,RHOINF,MA,RC,GAMMA,BIGGAMMA,EXP_VAL

CALL POINT_TO_MESH(NM)

UINF     = 35.0_EB         ! REFERENCE VELOCITY, M/S
PINF     = P_INF !101300._EB     ! AMBIENT PRESSURE, PA
RHOINF   = RHOA  !1.1717_EB      ! AMBIENT DENSITY, KG/M3
GAMMA    = 1.4_EB         ! RATIO OF SPECIFIC HEAT
RC       = 0.01556_EB     ! VORTEX RADIUS, M
MA       = UINF/SQRT(GAMMA*PINF/RHOINF) ! REFERENCE MACH NUMBER
EXP_VAL  = EXP(1.0_EB)
BIGGAMMA = 0.04*UINF*RC*SQRT(EXP_VAL)

DO K=0,KBAR
   DO J=0,JBAR
      DO I=0,IBAR
         U(I,J,K) = UINF-BIGGAMMA*(ZC(K)/(RC*RC))*EXP(-(X(I)*X(I)+ZC(K)*ZC(K))/(2._EB*RC*RC))
      ENDDO
   ENDDO
ENDDO

V=0._EB

DO K=0,KBAR
   DO J=0,JBAR
      DO I=0,IBAR
         W(I,J,K) = BIGGAMMA*(XC(I)/(RC*RC))*EXP(-(XC(I)*XC(I)+Z(K)*Z(K))/(2._EB*RC*RC))
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE TWOD_VORTEX_CERFACS


SUBROUTINE TWOD_VORTEX_UMD(NM)
!-------------------------------------------------------------------------------
! James White, University of Maryland
!-------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER, INTENT(IN) :: NM
INTEGER :: I,J,K
REAL(EB), PARAMETER :: RC    = 0.005_EB  ! Core radius of vortex (m)
REAL(EB), PARAMETER :: UMAX  = 0.5_EB    ! Vortex maximum velocity (m/s)
REAL(EB), PARAMETER :: XCLOC = 0._EB     ! Center of vortex, x (m)
REAL(EB), PARAMETER :: ZCLOC = 0._EB     ! Center of vortex, z (m)
REAL(EB), PARAMETER :: UCOF  = 0.1_EB    ! Coflow velocity (m/s)

CALL POINT_TO_MESH(NM)

DO K=0,KBAR
   DO J=0,JBAR
      DO I=0,IBAR
         U(I,J,K) = UMAX*((ZC(K)-ZCLOC)/RC)*EXP(0.5_EB- &
                    ((((X(I)-XCLOC)**2)+((ZC(K)-ZCLOC)**2))/(2._EB*RC**2)))+UCOF
      ENDDO
   ENDDO
ENDDO

V=0._EB

DO K=0,KBAR
   DO J=0,JBAR
      DO I=0,IBAR
         W(I,J,K) = -UMAX*((XC(I)-XCLOC)/RC)*EXP(0.5_EB- &
                    ((((XC(I)-XCLOC)**2)+((Z(K)-ZCLOC)**2))/(2._EB*RC**2)))
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE TWOD_VORTEX_UMD


SUBROUTINE TWOD_SOBOROT_UMD(NM)
!-------------------------------------------------------------------------------
! Salman Verma, University of Maryland
!
! Solid body rotation velocity field, u = z and w = -x.
!
! Used for PERIODIC_TEST==12,13
!-------------------------------------------------------------------------------
IMPLICIT NONE
INTEGER, INTENT(IN) :: NM
INTEGER :: I,J,K,IOR,II,JJ,KK,IW,IC
REAL(EB), PARAMETER :: USCAL = 1._EB     ! scale velocity (m/s)
REAL(EB), PARAMETER :: WSCAL = 1._EB     ! scale velocity (m/s)
REAL(EB), PARAMETER :: XCLOC = 0._EB     ! Center of vortex, x (m)
REAL(EB), PARAMETER :: ZCLOC = 0._EB     ! Center of vortex, z (m)

CALL POINT_TO_MESH(NM)

DO K=0,KBP1
   DO J=0,JBP1
      DO I=0,IBP1
         U(I,J,K) = USCAL*(ZC(K)-ZCLOC)
         US(I,J,K) = U(I,J,K)
      ENDDO
   ENDDO
ENDDO

V=0._EB
VS=0._EB

DO K=-0,KBP1
   DO J=0,JBP1
      DO I=0,IBP1
         W(I,J,K) = -WSCAL*(XC(I)-XCLOC)
         WS(I,J,K) = W(I,J,K)
      ENDDO
   ENDDO
ENDDO

! fill ghost values for smokeview

DO K=0,KBP1
   DO J=0,JBP1
      DO I=0,IBAR
         IC = CELL_INDEX(I,J,K)
         IF (IC==0) CYCLE
         V_EDGE_X(IC) = 0.5_EB*(V(I,J,K)+V(I+1,J,K))
         W_EDGE_X(IC) = 0.5_EB*(W(I,J,K)+W(I+1,J,K))
      ENDDO
   ENDDO
ENDDO

DO K=0,KBP1
   DO J=0,JBAR
      DO I=0,IBP1
         IC = CELL_INDEX(I,J,K)
         IF (IC==0) CYCLE
         U_EDGE_Y(IC) = 0.5_EB*(U(I,J,K)+U(I,J+1,K))
         W_EDGE_Y(IC) = 0.5_EB*(W(I,J,K)+W(I,J+1,K))
      ENDDO
   ENDDO
ENDDO

DO K=0,KBAR
   DO J=0,JBP1
      DO I=0,IBP1
         IC = CELL_INDEX(I,J,K)
         IF (IC==0) CYCLE
         U_EDGE_Z(IC) = 0.5_EB*(U(I,J,K)+U(I,J,K+1))
         V_EDGE_Z(IC) = 0.5_EB*(V(I,J,K)+V(I,J,K+1))
      ENDDO
   ENDDO
ENDDO

! Set normal velocity on external and internal boundaries (follows divg)

DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
   IOR = WALL(IW)%ONE_D%IOR
   II  = WALL(IW)%ONE_D%II
   JJ  = WALL(IW)%ONE_D%JJ
   KK  = WALL(IW)%ONE_D%KK
   SELECT CASE(IOR)
      CASE( 1)
         WALL(IW)%ONE_D%U_NORMAL_S = -U(II,JJ,KK)
      CASE(-1)
         WALL(IW)%ONE_D%U_NORMAL_S =  U(II-1,JJ,KK)
      CASE( 2)
         WALL(IW)%ONE_D%U_NORMAL_S = -V(II,JJ,KK)
      CASE(-2)
         WALL(IW)%ONE_D%U_NORMAL_S =  V(II,JJ-1,KK)
      CASE( 3)
         WALL(IW)%ONE_D%U_NORMAL_S = -W(II,JJ,KK)
      CASE(-3)
         WALL(IW)%ONE_D%U_NORMAL_S =  W(II,JJ,KK-1)
   END SELECT
   WALL(IW)%ONE_D%U_NORMAL = WALL(IW)%ONE_D%U_NORMAL_S
ENDDO

END SUBROUTINE TWOD_SOBOROT_UMD


SUBROUTINE VARDEN_DYNSMAG(NM)

INTEGER, INTENT(IN) :: NM

REAL(EB) :: TEMP_TERM,DUDY,DUDZ,DVDX,DVDZ,DWDX,DWDY,ONTHDIV
INTEGER :: I,J,K

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

IF (EVACUATION_ONLY(NM)) RETURN
CALL POINT_TO_MESH(NM)

IF (PREDICTOR) THEN
   UU=>U
   VV=>V
   WW=>W
   RHOP=>RHO
ELSE
   UU=>US
   VV=>VS
   WW=>WS
   RHOP=>RHOS
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

DO K = 1,KBAR
   DO J = 1,JBAR
      DO I = 1,IBAR

         UP(I,J,K) = 0.5_EB*(UU(I,J,K) + UU(I-1,J,K))
         VP(I,J,K) = 0.5_EB*(VV(I,J,K) + VV(I,J-1,K))
         WP(I,J,K) = 0.5_EB*(WW(I,J,K) + WW(I,J,K-1))

         S11(I,J,K) = RDX(I)*(UU(I,J,K)-UU(I-1,J,K))
         S22(I,J,K) = RDY(J)*(VV(I,J,K)-VV(I,J-1,K))
         S33(I,J,K) = RDZ(K)*(WW(I,J,K)-WW(I,J,K-1))

         ONTHDIV = ONTH*(S11(I,J,K)+S22(I,J,K)+S33(I,J,K))
         S11(I,J,K) = S11(I,J,K)-ONTHDIV
         S22(I,J,K) = S22(I,J,K)-ONTHDIV
         S33(I,J,K) = S33(I,J,K)-ONTHDIV

         DUDY = 0.25_EB*RDY(J)*(UU(I,J+1,K)-UU(I,J-1,K)+UU(I-1,J+1,K)-UU(I-1,J-1,K))
         DUDZ = 0.25_EB*RDZ(K)*(UU(I,J,K+1)-UU(I,J,K-1)+UU(I-1,J,K+1)-UU(I-1,J,K-1))
         DVDX = 0.25_EB*RDX(I)*(VV(I+1,J,K)-VV(I-1,J,K)+VV(I+1,J-1,K)-VV(I-1,J-1,K))
         DVDZ = 0.25_EB*RDZ(K)*(VV(I,J,K+1)-VV(I,J,K-1)+VV(I,J-1,K+1)-VV(I,J-1,K-1))
         DWDX = 0.25_EB*RDX(I)*(WW(I+1,J,K)-WW(I-1,J,K)+WW(I+1,J,K-1)-WW(I-1,J,K-1))
         DWDY = 0.25_EB*RDY(J)*(WW(I,J+1,K)-WW(I,J-1,K)+WW(I,J+1,K-1)-WW(I,J-1,K-1))
         S12(I,J,K) = 0.5_EB*(DUDY+DVDX)
         S13(I,J,K) = 0.5_EB*(DUDZ+DWDX)
         S23(I,J,K) = 0.5_EB*(DVDZ+DWDY)

         ! calculate magnitude of the grid strain rate

         SS(I,J,K) = RHOP(I,J,K)*STRAIN_RATE(I,J,K)

      ENDDO
   ENDDO
ENDDO

! second-order extrapolation to ghost cells

CALL EX2G3D(UP,-1.E10_EB,1.E10_EB)
CALL EX2G3D(VP,-1.E10_EB,1.E10_EB)
CALL EX2G3D(WP,-1.E10_EB,1.E10_EB)

CALL EX2G3D(S11,-1.E10_EB,1.E10_EB)
CALL EX2G3D(S22,-1.E10_EB,1.E10_EB)
CALL EX2G3D(S33,-1.E10_EB,1.E10_EB)
CALL EX2G3D(S12,-1.E10_EB,1.E10_EB)
CALL EX2G3D(S13,-1.E10_EB,1.E10_EB)
CALL EX2G3D(S23,-1.E10_EB,1.E10_EB)
CALL EX2G3D(SS,0._EB,1.E10_EB)

! test filter the strain rate

SHAT11 => TURB_WORK4
SHAT22 => TURB_WORK5
SHAT33 => TURB_WORK6
SHAT12 => TURB_WORK7
SHAT13 => TURB_WORK8
SHAT23 => TURB_WORK9
SSHAT  => TURB_WORK10

CALL TEST_FILTER(SHAT11,S11)
CALL TEST_FILTER(SHAT22,S22)
CALL TEST_FILTER(SHAT33,S33)
CALL TEST_FILTER(SHAT12,S12)
CALL TEST_FILTER(SHAT13,S13)
CALL TEST_FILTER(SHAT23,S23)


! calculate magnitude of test filtered strain rate

DO K = 1,KBAR
   DO J = 1,JBAR
      DO I = 1,IBAR
         SSHAT(I,J,K) = SQRT(2._EB*(SHAT11(I,J,K)*SHAT11(I,J,K) + &
                                    SHAT22(I,J,K)*SHAT22(I,J,K) + &
                                    SHAT33(I,J,K)*SHAT33(I,J,K) + &
                             2._EB*(SHAT12(I,J,K)*SHAT12(I,J,K) + &
                                    SHAT13(I,J,K)*SHAT13(I,J,K) + &
                                    SHAT23(I,J,K)*SHAT23(I,J,K)) ) )
      ENDDO
   ENDDO
ENDDO

! calculate the grid filtered stress tensor, beta

BETA11 => WORK1
BETA22 => WORK2
BETA33 => WORK3
BETA12 => WORK4
BETA13 => WORK5
BETA23 => WORK6

BETA11 = SS*S11
BETA22 = SS*S22
BETA33 = SS*S33
BETA12 = SS*S12
BETA13 = SS*S13
BETA23 = SS*S23

! ghost values for beta_ij should be filled already

! test filter the grid filtered stress tensor

BETAHAT11 => WORK1
BETAHAT22 => WORK2
BETAHAT33 => WORK3
BETAHAT12 => WORK4
BETAHAT13 => WORK5
BETAHAT23 => WORK6

WORK9=BETA11; CALL TEST_FILTER(BETAHAT11,WORK9)
WORK9=BETA22; CALL TEST_FILTER(BETAHAT22,WORK9)
WORK9=BETA33; CALL TEST_FILTER(BETAHAT33,WORK9)
WORK9=BETA12; CALL TEST_FILTER(BETAHAT12,WORK9)
WORK9=BETA13; CALL TEST_FILTER(BETAHAT13,WORK9)
WORK9=BETA23; CALL TEST_FILTER(BETAHAT23,WORK9)

! test filter the density

RHOPHAT => WORK7
CALL TEST_FILTER(RHOPHAT,RHOP)

! calculate the Mij tensor

M11 => WORK1
M22 => WORK2
M33 => WORK3
M12 => WORK4
M13 => WORK5
M23 => WORK6

DO K = 1,KBAR
   DO J = 1,JBAR
      DO I = 1,IBAR
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

! calculate the Leonard term, Lij

L11 => TURB_WORK4
L22 => TURB_WORK5
L33 => TURB_WORK6
L12 => TURB_WORK7
L13 => TURB_WORK8
L23 => TURB_WORK9

CALL CALC_VARDEN_LEONARD_TERM

! calculate Mij*Lij & Mij*Mij

MM    => TURB_WORK1
MMHAT => TURB_WORK1

ML    => TURB_WORK2
MLHAT => TURB_WORK2

DO K = 1,KBAR
   DO J = 1,JBAR
      DO I = 1,IBAR

         ML(I,J,K) = M11(I,J,K)*L11(I,J,K) + M22(I,J,K)*L22(I,J,K) + M33(I,J,K)*L33(I,J,K) + &
              2._EB*(M12(I,J,K)*L12(I,J,K) + M13(I,J,K)*L13(I,J,K) + M23(I,J,K)*L23(I,J,K))

         MM(I,J,K) = M11(I,J,K)*M11(I,J,K) + M22(I,J,K)*M22(I,J,K) + M33(I,J,K)*M33(I,J,K) + &
              2._EB*(M12(I,J,K)*M12(I,J,K) + M13(I,J,K)*M13(I,J,K) + M23(I,J,K)*M23(I,J,K))

      ENDDO
   ENDDO
ENDDO

! extrapolate to ghost

CALL EX2G3D(ML,0._EB,1.E10_EB)
CALL EX2G3D(MM,0._EB,1.E10_EB)

! do some smoothing

WORK9=ML; CALL TEST_FILTER(MLHAT,WORK9)
WORK9=MM; CALL TEST_FILTER(MMHAT,WORK9)

DO K = 1,KBAR
   DO J = 1,JBAR
      DO I = 1,IBAR

         ! calculate the local Smagorinsky coefficient

         ! perform "clipping" in case MLij is negative
         IF (MLHAT(I,J,K) < 0._EB) MLHAT(I,J,K) = 0._EB

         ! calculate the effective viscosity

         ! handle the case where we divide by zero, note MMHAT is positive semi-definite
         IF (MMHAT(I,J,K)<=TWO_EPSILON_EB) THEN
            CSD2(I,J,K) = 0._EB
         ELSE
            CSD2(I,J,K) = MLHAT(I,J,K)/MMHAT(I,J,K) ! (Cs*Delta)**2
         ENDIF

      END DO
   END DO
END DO

END SUBROUTINE VARDEN_DYNSMAG


SUBROUTINE CALC_VARDEN_LEONARD_TERM

REAL(EB), POINTER, DIMENSION(:,:,:) :: L11,L22,L33,L12,L13,L23
REAL(EB), POINTER, DIMENSION(:,:,:) :: RHOP,RHOPHAT
REAL(EB), POINTER, DIMENSION(:,:,:) :: UP,VP,WP
REAL(EB), POINTER, DIMENSION(:,:,:) :: RUU,RVV,RWW,RUV,RUW,RVW
REAL(EB), POINTER, DIMENSION(:,:,:) :: RU,RV,RW
REAL(EB), POINTER, DIMENSION(:,:,:) :: RUU_HAT,RVV_HAT,RWW_HAT,RUV_HAT,RUW_HAT,RVW_HAT
REAL(EB), POINTER, DIMENSION(:,:,:) :: RU_HAT,RV_HAT,RW_HAT
REAL(EB) :: INV_RHOPHAT
INTEGER :: I,J,K

! *****************************************************************************
! CAUTION WHEN MODIFYING: The order in which the tensor components are computed
! is important because we overwrite pointers several times to conserve memory.
! *****************************************************************************

IF (PREDICTOR) THEN
   RHOP=>RHO
ELSE
   RHOP=>RHOS
ENDIF
RHOPHAT => WORK7

! Compute rho*UiUj

UP => TURB_WORK1 ! will be overwritten by RU
VP => TURB_WORK2
WP => TURB_WORK3

RUU => TURB_WORK4 ! will be overwritten by RUU_HAT
RVV => TURB_WORK5
RWW => TURB_WORK6
RUV => TURB_WORK7
RUW => TURB_WORK8
RVW => TURB_WORK9

RUU = RHOP*UP*UP
RVV = RHOP*VP*VP
RWW = RHOP*WP*WP
RUV = RHOP*UP*VP
RUW = RHOP*UP*WP
RVW = RHOP*VP*WP

! extrapolate to ghost cells

CALL EX2G3D(RUU,0.E10_EB,1.E10_EB)
CALL EX2G3D(RVV,0.E10_EB,1.E10_EB)
CALL EX2G3D(RWW,0.E10_EB,1.E10_EB)
CALL EX2G3D(RUV,-1.E10_EB,1.E10_EB)
CALL EX2G3D(RUW,-1.E10_EB,1.E10_EB)
CALL EX2G3D(RVW,-1.E10_EB,1.E10_EB)

! Test filter rho*UiUj

RUU_HAT => TURB_WORK4 ! will be overwritten by Lij
RVV_HAT => TURB_WORK5
RWW_HAT => TURB_WORK6
RUV_HAT => TURB_WORK7
RUW_HAT => TURB_WORK8
RVW_HAT => TURB_WORK9

WORK9=RUU; CALL TEST_FILTER(RUU_HAT,WORK9)
WORK9=RVV; CALL TEST_FILTER(RVV_HAT,WORK9)
WORK9=RWW; CALL TEST_FILTER(RWW_HAT,WORK9)
WORK9=RUV; CALL TEST_FILTER(RUV_HAT,WORK9)
WORK9=RUW; CALL TEST_FILTER(RUW_HAT,WORK9)
WORK9=RVW; CALL TEST_FILTER(RVW_HAT,WORK9)

! Compute rho*Ui

RU => TURB_WORK1 ! will be overwritten by RU_HAT
RV => TURB_WORK2
RW => TURB_WORK3

RU = RHOP*UP
RV = RHOP*VP
RW = RHOP*WP

! extrapolate to ghost cells

CALL EX2G3D(RU,-1.E10_EB,1.E10_EB)
CALL EX2G3D(RV,-1.E10_EB,1.E10_EB)
CALL EX2G3D(RW,-1.E10_EB,1.E10_EB)

! Test filter rho*Ui

RU_HAT => TURB_WORK1
RV_HAT => TURB_WORK2
RW_HAT => TURB_WORK3

WORK9=RU; CALL TEST_FILTER(RU_HAT,WORK9)
WORK9=RV; CALL TEST_FILTER(RV_HAT,WORK9)
WORK9=RW; CALL TEST_FILTER(RW_HAT,WORK9)

! Compute variable density Leonard stress

L11 => TURB_WORK4
L22 => TURB_WORK5
L33 => TURB_WORK6
L12 => TURB_WORK7
L13 => TURB_WORK8
L23 => TURB_WORK9

DO K = 1,KBAR
   DO J = 1,JBAR
      DO I = 1,IBAR
         INV_RHOPHAT = 1._EB/RHOPHAT(I,J,K)
         L11(I,J,K) = RUU_HAT(I,J,K) - RU_HAT(I,J,K)*RU_HAT(I,J,K)*INV_RHOPHAT
         L22(I,J,K) = RVV_HAT(I,J,K) - RV_HAT(I,J,K)*RV_HAT(I,J,K)*INV_RHOPHAT
         L33(I,J,K) = RWW_HAT(I,J,K) - RW_HAT(I,J,K)*RW_HAT(I,J,K)*INV_RHOPHAT
         L12(I,J,K) = RUV_HAT(I,J,K) - RU_HAT(I,J,K)*RV_HAT(I,J,K)*INV_RHOPHAT
         L13(I,J,K) = RUW_HAT(I,J,K) - RU_HAT(I,J,K)*RW_HAT(I,J,K)*INV_RHOPHAT
         L23(I,J,K) = RVW_HAT(I,J,K) - RV_HAT(I,J,K)*RW_HAT(I,J,K)*INV_RHOPHAT
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE CALC_VARDEN_LEONARD_TERM


SUBROUTINE EX2G3D(A,A_MIN,A_MAX)

! Second order extrapolation of 3D array to ghost cells

REAL(EB), INTENT(IN) :: A_MIN,A_MAX
REAL(EB), INTENT(INOUT) :: A(0:IBP1,0:JBP1,0:KBP1)

A(0,:,:) = MIN(A_MAX,MAX(A_MIN,2._EB*A(1,:,:)-A(2,:,:)))
A(:,0,:) = MIN(A_MAX,MAX(A_MIN,2._EB*A(:,1,:)-A(:,2,:)))
A(:,:,0) = MIN(A_MAX,MAX(A_MIN,2._EB*A(:,:,1)-A(:,:,2)))

A(IBP1,:,:) = MIN(A_MAX,MAX(A_MIN,2._EB*A(IBAR,:,:)-A(IBM1,:,:)))
A(:,JBP1,:) = MIN(A_MAX,MAX(A_MIN,2._EB*A(:,JBAR,:)-A(:,JBM1,:)))
A(:,:,KBP1) = MIN(A_MAX,MAX(A_MIN,2._EB*A(:,:,KBAR)-A(:,:,KBM1)))

END SUBROUTINE EX2G3D


SUBROUTINE FILL_EDGES(A)

! Extrapolate A to egdes and corners of 3D arrays

REAL(EB), INTENT(INOUT) :: A(0:IBP1,0:JBP1,0:KBP1)
INTEGER :: I,J,K

!    ---------------------
!    |         |         |
!    |    A    |    C    |
!    |         |         |
!    |---------x----------
!    |         |         |
!    |    O    |    B    |
!    |         |         |
!    ---------------------
!
!  For the edge above (x pointing out of the screen) O is the unknown.
!  It is obtained from a 2nd-order extrapolation.  First, A and B are
!  averaged to point x, then O = 2*(.5*(A+B)) - C = (A+B) - C.
!
!  Note that for edges of an interpolated mesh, the values A and B
!  are obtained from a mesh exchange and populated as wall ghost cell
!  values in U_GHOST(IW), for example.

! x edges

J=0; K=0
DO I=1,IBAR
   A(I,J,K) = ( A(I,J+1,K) + A(I,J,K+1) ) - A(I,J+1,K+1)
ENDDO

J=0; K=KBP1
DO I=1,IBAR
   A(I,J,K) = ( A(I,J+1,K) + A(I,J,K-1) ) - A(I,J+1,K-1)
ENDDO

J=JBP1; K=0
DO I=1,IBAR
   A(I,J,K) = ( A(I,J-1,K) + A(I,J,K+1) ) - A(I,J-1,K+1)
ENDDO

J=JBP1; K=KBP1
DO I=1,IBAR
   A(I,J,K) = ( A(I,J-1,K) + A(I,J,K-1) ) - A(I,J-1,K-1)
ENDDO

! y edges

I=0; K=0
DO J=1,JBAR
   A(I,J,K) = ( A(I+1,J,K) + A(I,J,K+1) ) - A(I+1,J,K+1)
ENDDO

I=0; K=KBP1
DO J=1,JBAR
   A(I,J,K) = ( A(I+1,J,K) + A(I,J,K-1) ) - A(I+1,J,K-1)
ENDDO

I=IBP1; K=0
DO J=1,JBAR
   A(I,J,K) = ( A(I-1,J,K) + A(I,J,K+1) ) - A(I-1,J,K+1)
ENDDO

I=IBP1; K=KBP1
DO J=1,JBAR
   A(I,J,K) = ( A(I-1,J,K) + A(I,J,K-1) ) - A(I-1,J,K-1)
ENDDO

! z edges

I=0; J=0
DO K=1,KBAR
   A(I,J,K) = ( A(I+1,J,K) + A(I,J+1,K) ) - A(I+1,J+1,K)
ENDDO

I=0; J=JBP1
DO K=1,KBAR
   A(I,J,K) = ( A(I+1,J,K) + A(I,J-1,K) ) - A(I+1,J-1,K)
ENDDO

I=IBP1; J=0
DO K=1,KBAR
   A(I,J,K) = ( A(I-1,J,K) + A(I,J+1,K) ) - A(I-1,J+1,K)
ENDDO

I=IBP1; J=JBP1
DO K=1,KBAR
   A(I,J,K) = ( A(I-1,J,K) + A(I,J-1,K) ) - A(I-1,J-1,K)
ENDDO

! Corners

A(0,0,0) = 2._EB*A(1,1,1) - A(2,2,2)
A(IBP1,0,0) = 2._EB*A(IBP1-1,1,1) - A(IBP1-2,2,2)
A(0,JBP1,0) = 2._EB*A(1,JBP1-1,1) - A(2,JBP1-2,2)
A(0,0,KBP1) = 2._EB*A(1,1,KBP1-1) - A(2,2,KBP1-2)
A(IBP1,JBP1,0) = 2._EB*A(IBP1-1,JBP1-1,1) - A(IBP1-2,JBP1-2,2)
A(IBP1,0,KBP1) = 2._EB*A(IBP1-1,1,KBP1-1) - A(IBP1-2,2,KBP1-2)
A(0,JBP1,KBP1) = 2._EB*A(1,JBP1-1,KBP1-1) - A(2,JBP1-2,KBP1-2)
A(IBP1,JBP1,KBP1) = 2._EB*A(IBP1-1,JBP1-1,KBP1-1) - A(IBP1-2,JBP1-2,KBP1-2)

END SUBROUTINE FILL_EDGES


SUBROUTINE TEST_FILTER(HAT,ORIG)

! Tophat filter

REAL(EB), INTENT(IN) :: ORIG(0:IBP1,0:JBP1,0:KBP1)
REAL(EB), INTENT(OUT) :: HAT(0:IBP1,0:JBP1,0:KBP1)
INTEGER :: I, J, K, L, M, N
REAL(EB), PARAMETER :: K1DM(3)=(/1.0_EB,1.0_EB,1.0_EB/)
REAL(EB), PARAMETER :: K3DM(-1:1,-1:1,-1:1)=RESHAPE((/(((K1DM(I)*K1DM(J)*K1DM(K)/27._EB,I=1,3),J=1,3),K=1,3)/),(/3,3,3/))
REAL(EB), PARAMETER :: K1DT(3)=(/1.0_EB,2.0_EB,1.0_EB/)
REAL(EB), PARAMETER :: K3DT(-1:1,-1:1,-1:1)=RESHAPE((/(((K1DT(I)*K1DT(J)*K1DT(K)/64._EB,I=1,3),J=1,3),K=1,3)/),(/3,3,3/))
REAL(EB), PARAMETER :: K1DS(3)=(/1.0_EB,4.0_EB,1.0_EB/)
REAL(EB), PARAMETER :: K3DS(-1:1,-1:1,-1:1)=RESHAPE((/(((K1DS(I)*K1DS(J)*K1DS(K)/216._EB,I=1,3),J=1,3),K=1,3)/),(/3,3,3/))

! Traverse bulk of mesh

QUADRATURE_SELECT: SELECT CASE(TEST_FILTER_QUADRATURE)

   CASE(TRAPAZOID_QUADRATURE) ! default

      !$OMP PARALLEL
      !$OMP DO SCHEDULE(static)
      DO K = 1,KBP1-1
         DO J = 1,JBP1-1
            DO I = 1,IBP1-1

               ! Apply 3x3x3 Kernel; this is faster than elementwise array multiplication.
               HAT(I,J,K) = 0._EB
               DO N = -1,1
                  DO M = -1,1
                     DO L = -1,1
                        HAT(I,J,K) = HAT(I,J,K) + ORIG(I+L,J+M,K+N) * K3DT(L,M,N)
                     ENDDO
                  ENDDO
               ENDDO

            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO
      !$OMP END PARALLEL

   CASE(SIMPSON_QUADRATURE)

      !$OMP PARALLEL
      !$OMP DO SCHEDULE(static)
      DO K = 1,KBP1-1
         DO J = 1,JBP1-1
            DO I = 1,IBP1-1

               ! Apply 3x3x3 Kernel; this is faster than elementwise array multiplication.
               HAT(I,J,K) = 0._EB
               DO N = -1,1
                  DO M = -1,1
                     DO L = -1,1
                        HAT(I,J,K) = HAT(I,J,K) + ORIG(I+L,J+M,K+N) * K3DS(L,M,N)
                     ENDDO
                  ENDDO
               ENDDO

            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO
      !$OMP END PARALLEL

   CASE(MIDPOINT_QUADRATURE)

      !$OMP PARALLEL
      !$OMP DO SCHEDULE(static)
      DO K = 1,KBP1-1
         DO J = 1,JBP1-1
            DO I = 1,IBP1-1

               ! Apply 3x3x3 Kernel; this is faster than elementwise array multiplication.
               HAT(I,J,K) = 0._EB
               DO N = -1,1
                  DO M = -1,1
                     DO L = -1,1
                        HAT(I,J,K) = HAT(I,J,K) + ORIG(I+L,J+M,K+N) * K3DM(L,M,N)
                     ENDDO
                  ENDDO
               ENDDO

            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO
      !$OMP END PARALLEL

END SELECT QUADRATURE_SELECT

! Traverse shell of mesh rather crudely.
! Edges and corners are calculated several times.

!$OMP PARALLEL
!$OMP DO SCHEDULE(static)
DO K = 0,KBP1
   DO J = 0,JBP1
      HAT(0,J,K) = 2._EB * HAT(0+1,J,K) - HAT(0+2,J,K)
      HAT(IBP1,J,K) = 2._EB * HAT(IBP1-1,J,K) - HAT(IBP1-2,J,K)
   END DO
END DO
!$OMP END DO

!$OMP DO SCHEDULE(static)
DO K = 0,KBP1
   DO I = 0,IBP1
      HAT(I,0,K) = 2._EB * HAT(I,0+1,K) - HAT(I,0+2,K)
      HAT(I,JBP1,K) = 2._EB * HAT(I,JBP1-1,K) - HAT(I,JBP1-2,K)
   END DO
END DO
!$OMP END DO

!$OMP DO SCHEDULE(static)
DO J = 0,JBP1
   DO I = 0,IBP1
      HAT(I,J,0) = 2._EB * HAT(I,J,0+1) - HAT(I,J,0+2)
      HAT(I,J,KBP1) = 2._EB * HAT(I,J,KBP1-1) - HAT(I,J,KBP1-2)
   END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE TEST_FILTER


SUBROUTINE RNG_EDDY_VISCOSITY(MU_EFF,MU_DNS_RNG,RHO,STRAIN_RATE,DELTA)

! A. Yakhot, S. A. Orszag, V. Yakhot, and M. Israeli. Renormalization Group Formulation of Large-Eddy Simulation.
! Journal of Scientific Computing, 1(1):1-51, 1989.

REAL(EB), INTENT(INOUT) :: MU_EFF
REAL(EB), INTENT(IN) :: MU_DNS_RNG,STRAIN_RATE,DELTA,RHO
REAL(EB) :: MU_S,MU2,MU3,MU_EFF_OLD
INTEGER :: ITER
INTEGER, PARAMETER :: MAX_IT=20
REAL(EB), PARAMETER :: TOL=1.E-3_EB

MU_S = RHO*(C_RNG*DELTA)**2*STRAIN_RATE
MU_EFF = MU_DNS_RNG + MU_S
MU2 = MU_S**2
MU3 = MU_DNS_RNG**3
IF (MU3>TWO_EPSILON_EB) THEN
   DO ITER=1,MAX_IT
      MU_EFF_OLD = MU_EFF
      MU_EFF = MU_DNS_RNG*( 1._EB + MAX(0._EB, MU2*MU_EFF/MU3 - C_RNG_CUTOFF) )**ONTH
      !print *, ITER, MU_EFF
      IF (ABS((MU_EFF-MU_EFF_OLD)/MU_EFF_OLD)<TOL) EXIT
   ENDDO
ENDIF

END SUBROUTINE RNG_EDDY_VISCOSITY


SUBROUTINE WALE_VISCOSITY(NU_T,G_IJ,DELTA)

! Wall Adapting Local Eddy-viscosity (WALE)

! F. Nicoud, F. Ducros. Subgrid-scale stress modelling based on the square of the velocity gradient tensor.
! Flow, Turbulence, and Combustion, Vol. 62, pp. 183-200, 1999.

REAL(EB), INTENT(OUT) :: NU_T
REAL(EB), INTENT(IN) :: G_IJ(3,3),DELTA
REAL(EB) :: S_IJ(3,3),O_IJ(3,3),S2,O2,IV_SO,SD2,DENOM
INTEGER :: I,J,K,L

! compute strain and rotation tensors

DO J=1,3
   DO I=1,3
      S_IJ(I,J) = 0.5_EB * ( G_IJ(I,J) + G_IJ(J,I) )
      O_IJ(I,J) = 0.5_EB * ( G_IJ(I,J) - G_IJ(J,I) )
   ENDDO
ENDDO

! contraction of strain and rotation tensors

S2 = 0._EB
O2 = 0._EB
DO J=1,3
   DO I=1,3
      S2 = S2 + S_IJ(I,J)*S_IJ(I,J)
      O2 = O2 + O_IJ(I,J)*O_IJ(I,J)
   ENDDO
ENDDO

! fourth order contraction

IV_SO = 0._EB
DO L=1,3
   DO K=1,3
      DO J=1,3
         DO I=1,3
            IV_SO = IV_SO + S_IJ(I,K)*S_IJ(K,J)*O_IJ(J,L)*O_IJ(L,I)
         ENDDO
      ENDDO
   ENDDO
ENDDO

! using Caley-Hamilton theorem

SD2 = ONSI*(S2*S2 + O2*O2) + TWTH*S2*O2 + 2._EB*IV_SO
IF (SD2 < 0.0) SD2 = 0._EB

DENOM = S2**2.5_EB + SD2**1.25_EB
IF (DENOM>TWO_EPSILON_EB) THEN
   NU_T = (C_WALE*DELTA)**2 * SD2**1.5_EB / DENOM
ELSE
   NU_T = 0._EB
ENDIF

END SUBROUTINE WALE_VISCOSITY


SUBROUTINE WALL_MODEL(SLIP_FACTOR,U_TAU,Y_PLUS,U,NU,DY,S,DY2,U_DY2)

REAL(EB), INTENT(OUT) :: SLIP_FACTOR,U_TAU,Y_PLUS
REAL(EB), INTENT(IN) :: U,NU,DY,S ! S is the roughness length scale (Pope's notation)
REAL(EB), OPTIONAL, INTENT(IN)  :: DY2
REAL(EB), OPTIONAL, INTENT(OUT) :: U_DY2

REAL(EB), PARAMETER :: RKAPPA=1._EB/0.41_EB ! 1/von Karman constant
REAL(EB), PARAMETER :: B=5.2_EB,BTILDE_ROUGH=8.5_EB,BTILDE_MAX=9.5_EB ! see Pope (2000) pp. 294,297,298
REAL(EB), PARAMETER :: S0=1._EB,S1=5.83_EB,S2=30._EB ! approx piece-wise function for Fig. 7.24, Pope (2000) p. 297
REAL(EB), PARAMETER :: Y1=5._EB,Y2=30._EB
REAL(EB), PARAMETER :: U1=5._EB,U2=RKAPPA*LOG(Y2)+B
REAL(EB), PARAMETER :: EPS=1.E-10_EB

REAL(EB) :: Y_CELL_CENTER,TAU_W,BTILDE,DELTA_NU,S_PLUS,DUDY,Y_CELL_CENTER2,Y_PLUS2
INTEGER :: ITER

! References:
!
! S. B. Pope (2000) Turbulent Flows, Cambridge.

! Step 1: compute laminar (DNS) stress, and initial guess for LES stress

Y_CELL_CENTER = 0.5_EB*DY
DUDY = ABS(U)/Y_CELL_CENTER
TAU_W = NU*DUDY                         ! actually tau_w/rho
U_TAU = SQRT(ABS(TAU_W))                ! friction velocity
DELTA_NU = NU/(U_TAU+EPS)               ! viscous length scale
Y_PLUS = Y_CELL_CENTER/(DELTA_NU+EPS)
SLIP_FACTOR = -1._EB

! Step 2: compute turbulent (LES) stress

LES_IF: IF (SIM_MODE/=DNS_MODE) THEN

   ! NOTE: 2 iterations converges TAU_W to roughly 5 % residual error
   !       3 iterations converges TAU_W to roughly 1 % residual error

   DO ITER=1,3

      S_PLUS = S/(DELTA_NU+EPS) ! roughness in viscous units

      IF (S_PLUS < S0) THEN
         ! smooth wall
         Y_PLUS = Y_CELL_CENTER/(DELTA_NU+EPS)
         IF (Y_PLUS < Y_WERNER_WENGLE) THEN
            ! viscous sublayer
            TAU_W = ( U/Y_PLUS )**2
            U_TAU = SQRT(TAU_W)
            DUDY = ABS(U)/Y_CELL_CENTER
         ! ELSEIF (Y_PLUS < Y2) THEN
         !   ! buffer layer
         !   TAU_W = ( U/U_PLUS_BUFFER_SEMILOG(Y_PLUS) )**2
         !   U_TAU = SQRT(TAU_W)
         !   DUDY = 0.5_EB*(ABS(U)/Y_CELL_CENTER + U_TAU*RKAPPA/Y_CELL_CENTER)
         ELSE
            ! log layer
            TAU_W = ( U/(RKAPPA*LOG(Y_PLUS)+B) )**2
            U_TAU = SQRT(TAU_W)
            DUDY = U_TAU*RKAPPA/Y_CELL_CENTER
         ENDIF
      ELSE
         ! rough wall
         IF (S_PLUS < S1) THEN
            BTILDE = B + RKAPPA*LOG(S_PLUS) ! Pope (2000) p. 297, Eq. (7.122)
         ELSE IF (S_PLUS < S2) THEN
            BTILDE = BTILDE_MAX ! approximation from Fig. 7.24, Pope (2000) p. 297
         ELSE
            BTILDE = BTILDE_ROUGH ! fully rough
         ENDIF
         Y_PLUS = Y_CELL_CENTER/S
         TAU_W = ( U/(RKAPPA*LOG(Y_PLUS)+BTILDE) )**2  ! Pope (2000) p. 297, Eq. (7.121)
         U_TAU = SQRT(TAU_W)
         DUDY = U_TAU*RKAPPA/Y_CELL_CENTER
      ENDIF

      DELTA_NU = NU/(U_TAU+EPS)

   ENDDO

   ! NOTE: SLIP_FACTOR is no longer used to compute the wall stress, see VELOCITY_BC.
   ! The stress is taken directly from U_TAU. SLIP_FACTOR is, however, still used to
   ! compute the velocity gradient at the wall that feeds into the wall vorticity.
   ! Since the gradients implied by the wall function can be large and lead to instabilities,
   ! we bound the wall slip between no slip and free slip.

   ! The slip factor (SF) is based on the following approximation to the wall gradient
   ! (note that u0 is the ghost cell value of the streamwise velocity component and
   ! y is the wall-normal direction):
   ! dudy = (u-u0)/dy = (u-SF*u)/dy = u/dy*(1-SF) => SF = 1 - dudy*dy/u
   ! In this routine, dudy is sampled from the wall model at the location y_cell_center.

   SLIP_FACTOR = MAX(-1._EB,MIN(1._EB,1._EB-DUDY*DY/(ABS(U)+EPS))) ! -1.0 <= SLIP_FACTOR <= 1.0

ENDIF LES_IF

DY2_IF : IF (PRESENT(DY2)) THEN
   IF (SIM_MODE==DNS_MODE) THEN
      Y_CELL_CENTER2 = 0.5_EB*DY2
      U_DY2 = Y_CELL_CENTER2/Y_CELL_CENTER * U ! Linear Variation of velocities is assumed.
   ELSE
      Y_CELL_CENTER2 = 0.5_EB*DY2
      S_PLUS = S/(DELTA_NU+EPS) ! roughness in viscous units
      IF (S_PLUS < S0) THEN
         ! smooth wall
         Y_PLUS2 = Y_CELL_CENTER2/(DELTA_NU+EPS)
         IF (Y_PLUS2 < Y_WERNER_WENGLE) THEN
            ! viscous sublayer
            U_DY2 = Y_CELL_CENTER2/Y_CELL_CENTER * U ! Linear Variation of velocities is assumed.
         ELSE
            ! log layer
            U_DY2   = U_TAU*(RKAPPA*LOG(Y_PLUS2)+B) ! U_TAU*U_PLUS2
         ENDIF
      ELSE
         ! rough wall
         IF (S_PLUS < S1) THEN
            BTILDE = B + RKAPPA*LOG(S_PLUS) ! Pope (2000) p. 297, Eq. (7.122)
         ELSE IF (S_PLUS < S2) THEN
            BTILDE = BTILDE_MAX ! approximation from Fig. 7.24, Pope (2000) p. 297
         ELSE
            BTILDE = BTILDE_ROUGH ! fully rough
         ENDIF
         Y_PLUS2 = Y_CELL_CENTER2/S
         U_DY2   = U_TAU*(RKAPPA*LOG(Y_PLUS2)+BTILDE) ! U_TAU*U_PLUS2 Pope (2000) p. 297, Eq. (7.121)
      ENDIF
   ENDIF
ENDIF DY2_IF

CONTAINS

REAL(EB) FUNCTION U_PLUS_BUFFER_SEMILOG(YP)

REAL(EB), INTENT(IN) :: YP
REAL(EB), PARAMETER :: RKAPPA_BUFFER=(U2-U1)/(LOG(Y2)-LOG(Y1))
REAL(EB), PARAMETER :: B_BUFFER=U1-RKAPPA_BUFFER*LOG(Y1)

! semi-log fit connecting U1=Y1=5 to U2=RKAPPA*LOG(Y2)+B at Y2=30

U_PLUS_BUFFER_SEMILOG = RKAPPA_BUFFER*LOG(YP)+B_BUFFER

END FUNCTION U_PLUS_BUFFER_SEMILOG

REAL(EB) FUNCTION U_PLUS_BUFFER_POLY4(YP)

REAL(EB), INTENT(IN) :: YP
REAL(EB) :: DYP
REAL(EB), PARAMETER :: DYPLUS=25._EB
REAL(EB), PARAMETER :: B1 = (U2-Y2)/DYPLUS**2
REAL(EB), PARAMETER :: B2 = (RKAPPA/Y2-1._EB)/DYPLUS
REAL(EB), PARAMETER :: B3 = (-RKAPPA/Y2**2)*0.5_EB
REAL(EB), PARAMETER :: C3 = 6._EB*B1-3._EB*B2+B3
REAL(EB), PARAMETER :: C2 = (4._EB*B1-B2 - 2._EB*C3)/DYPLUS
REAL(EB), PARAMETER :: C1 = (B1-C3-DYPLUS*C2)/DYPLUS**2

! Jung-il Choi, Yonsei University
! 4th-order polynomial fit connecting U1=Y1=5 to U2=RKAPPA*LOG(Y2)+B at Y2=30

DYP = YP-Y1
U_PLUS_BUFFER_POLY4 = C1*DYP**4 + C2*DYP**3 + C3*DYP**2 + YP

END FUNCTION U_PLUS_BUFFER_POLY4

END SUBROUTINE WALL_MODEL


SUBROUTINE NATURAL_CONVECTION_MODEL(H_NATURAL,DELTA_TMP,C_VERTICAL,C_HORIZONTAL,SURF_GEOMETRY_INDEX,IOR,K_G,DN)

REAL(EB), INTENT(OUT) :: H_NATURAL
REAL(EB), INTENT(IN) :: DELTA_TMP,C_VERTICAL,C_HORIZONTAL,K_G,DN
INTEGER, INTENT(IN) :: SURF_GEOMETRY_INDEX,IOR

! Calculate the HTC for natural/free convection (Holman, 1990, Table 7-2)

SELECT CASE(SURF_GEOMETRY_INDEX)
   CASE (SURF_CARTESIAN)
      SELECT CASE(ABS(IOR))
         CASE(0:2)
            H_NATURAL = C_VERTICAL*ABS(DELTA_TMP)**ONTH
         CASE(3)
            H_NATURAL = C_HORIZONTAL*ABS(DELTA_TMP)**ONTH
      END SELECT
      H_NATURAL = MAX(H_NATURAL,2._EB*K_G/DN)

   CASE (SURF_CYLINDRICAL)
      H_NATURAL = C_VERTICAL*ABS(DELTA_TMP)**ONTH

   CASE (SURF_SPHERICAL) ! It is assumed that the forced HTC represents natural convection as well
      H_NATURAL = 0._EB
END SELECT

END SUBROUTINE NATURAL_CONVECTION_MODEL


SUBROUTINE FORCED_CONVECTION_MODEL(H_FORCED,RE,K_G,CONV_LENGTH,SURF_GEOMETRY_INDEX)

REAL(EB), INTENT(OUT) :: H_FORCED
REAL(EB), INTENT(IN) :: RE,K_G,CONV_LENGTH
INTEGER, INTENT(IN) :: SURF_GEOMETRY_INDEX
REAL(EB) :: NUSSELT

SELECT CASE(SURF_GEOMETRY_INDEX)
   CASE (SURF_CARTESIAN)
      ! Incropera and DeWitt, 3rd, 1990, Eq. 7.44
      NUSSELT = 0.037_EB*RE**0.8_EB*PR_ONTH
   CASE (SURF_CYLINDRICAL)
      ! Incropera and DeWitt, 3rd, 1990, Eq. 7.55
      IF (RE >= 40._EB) THEN
         NUSSELT = 0.683_EB*RE**0.466_EB*PR_ONTH
      ELSEIF (RE >= 4._EB) THEN
         NUSSELT = 0.911_EB*RE**0.385_EB*PR_ONTH
      ELSE
         NUSSELT = 0.989_EB*RE**0.330_EB*PR_ONTH
      ENDIF

   CASE (SURF_SPHERICAL)
      ! Incropera and DeWitt, 3rd, 1990, Eq. 7.59
      NUSSELT = 2._EB + 0.6_EB*SQRT(RE)*PR_ONTH
END SELECT
H_FORCED = MAX(1._EB,NUSSELT)*K_G/CONV_LENGTH

END SUBROUTINE FORCED_CONVECTION_MODEL


SUBROUTINE RAYLEIGH_HEAT_FLUX_MODEL(H,Z_STAR,DZ,TMP_W,TMP_G,K_G,RHO_G,CP_G,MU_G)

!!!!! EXPERIMENTAL !!!!!

! Rayleigh number scaling in nondimensional thermal wall units
!
! The formulation is based on the discussion of natural convection systems in
! J.P. Holman, Heat Transfer, 7th Ed., McGraw-Hill, 1990, p. 346.

REAL(EB), INTENT(OUT) :: H,Z_STAR
REAL(EB), INTENT(IN) :: DZ,TMP_W,TMP_G,K_G,RHO_G,CP_G,MU_G
REAL(EB) :: NUSSELT,Q,ZC,NU_G,DS,ALPHA,THETA,Q_OLD,ERROR
INTEGER :: ITER
INTEGER, PARAMETER :: MAX_ITER=10
! C_L = Z_L**(-0.8_EB)
! C_T = C_L*Z_T**(-0.2_EB)
REAL(EB), PARAMETER :: Z_L = 3.2_EB, Z_T=17._EB
REAL(EB), PARAMETER :: C_L = 3.2_EB**(-0.8_EB), C_T = 0.394_EB*17._EB**(-0.2_EB)

IF (ABS(TMP_W-TMP_G)<TWO_EPSILON_EB) THEN
   H = 0._EB
   Z_STAR = 0._EB
   RETURN
ENDIF

ZC = 0.5_EB*DZ
NU_G = MU_G/RHO_G
ALPHA = K_G/(RHO_G*CP_G)
THETA = TMP_W*K_G*ALPHA*NU_G/GRAV

! Step 1: assume a heat transfer coefficient

H = K_G/ZC ! initial guess
Q = H*ABS(TMP_W-TMP_G)

RAYLEIGH_LOOP: DO ITER=1,MAX_ITER

   ! Step 2: compute new thermal diffusive length scale, delta*, from modified Grashof number * Pr

   DS = (THETA/Q)**0.25_EB

   ! Step 3: compute new z* (thermal)

   Z_STAR = ZC/DS ! Ra* = (z*)**4

   ! Step 4: based on z*, choose Ra scaling law

   IF (Z_STAR<=Z_L) THEN
      NUSSELT = 1._EB
   ELSEIF (Z_STAR>Z_L .AND. Z_STAR<=Z_T) THEN
      NUSSELT = C_L * Z_STAR**0.8_EB
   ELSE
      NUSSELT = C_T * Z_STAR
   ENDIF

   ! Step 5: update heat transfer coefficient

   H = NUSSELT*K_G/ZC
   Q_OLD = Q
   Q = H*ABS(TMP_W-TMP_G)

   ERROR = ABS(Q-Q_OLD)/MAX(Q_OLD,TWO_EPSILON_EB)

   IF (ERROR<0.001_EB) EXIT RAYLEIGH_LOOP

ENDDO RAYLEIGH_LOOP

END SUBROUTINE RAYLEIGH_HEAT_FLUX_MODEL


SUBROUTINE YUAN_HEAT_FLUX_MODEL(H,Y_STAR,DY,TMP_W,TMP_G,K_G,RHO_G,CP_G)

!!!!! EXPERIMENTAL !!!!!

! This model is very similar to RAYLEIGH_HEAT_FLUX_MODEL
!
! X. Yuan, A. Moser, P. Suter. Wall functions for numerical simulation of turbulent
! natural convection along vertical plates. Int. J. Heat Mass Transfer, Vol. 36,
! No. 18 pp. 4477-4485, 1993.

REAL(EB), INTENT(OUT) :: H,Y_STAR
REAL(EB), INTENT(IN) :: DY,TMP_W,TMP_G,K_G,RHO_G,CP_G
REAL(EB) :: T_STAR,Q,YC,TQ,UQ,ALPHA,THETA,GAMMA,Q_OLD,ERROR
INTEGER :: ITER
INTEGER, PARAMETER :: MAX_ITER=10

IF (ABS(TMP_W-TMP_G)<TWO_EPSILON_EB) THEN
   H = 0._EB
   Y_STAR = 0._EB
   RETURN
ENDIF

YC = 0.5_EB*DY
ALPHA = K_G/(RHO_G*CP_G)
THETA = GRAV/TMP_W*ALPHA*(RHO_G*CP_G)**3
GAMMA = GRAV/TMP_W*ALPHA/(RHO_G*CP_G)

H = K_G/YC ! initial guess
Q = H*ABS(TMP_W-TMP_G)

YUAN_LOOP: DO ITER=1,MAX_ITER

   UQ = ( Q / GAMMA )**0.25_EB

   Y_STAR = YC*UQ/ALPHA ! Yuan et al. Eq. (15)

   ! Yuan et al. Eqs. (22)-(24)
   IF (Y_STAR<=1._EB) THEN
      T_STAR = Y_STAR
   ELSEIF (Y_STAR>1._EB .AND. Y_STAR<=100._EB) THEN
      T_STAR = 1._EB + 1.36_EB*LOG(Y_STAR) - 0.135_EB*LOG(Y_STAR)**2
   ELSE
      T_STAR = 4.4_EB ! curvefit for AIR
   ENDIF

   TQ = ABS(TMP_W-TMP_G)/MAX(T_STAR,TWO_EPSILON_EB)
   Q_OLD = Q
   Q = (THETA*TQ**4)**ONTH

   ERROR = ABS(Q-Q_OLD)/MAX(Q_OLD,TWO_EPSILON_EB)

   IF (ERROR<0.001_EB) EXIT YUAN_LOOP

ENDDO YUAN_LOOP

H = MAX( K_G/YC, Q/ABS(TMP_W-TMP_G) )

END SUBROUTINE YUAN_HEAT_FLUX_MODEL


SUBROUTINE LOGLAW_HEAT_FLUX_MODEL(H,YPLUS,U_TAU,K_G,RHO_G,CP_G,MU_G)

! Kiyoung Moon, Yonsei University
! Ezgi Oztekin, Technology and Manangement International

REAL(EB), INTENT(OUT) :: H
REAL(EB), INTENT(IN) :: YPLUS,U_TAU,K_G,RHO_G,CP_G,MU_G
REAL(EB) :: PR_M,TPLUS,B_T
REAL(EB), PARAMETER :: RKAPPA=1._EB/0.41_EB

PR_M = CP_G*MU_G/K_G

IF (YPLUS < Y_WERNER_WENGLE) THEN
   TPLUS = PR_M*YPLUS
ELSE
   B_T = (3.85_EB*PR_M**ONTH-1.3_EB)**2 + 2.12_EB*LOG(PR_M) ! Kader, 1981
   TPLUS = PR*RKAPPA*LOG(YPLUS)+B_T
ENDIF

H = RHO_G*U_TAU*CP_G/TPLUS

END SUBROUTINE LOGLAW_HEAT_FLUX_MODEL


SUBROUTINE ABL_HEAT_FLUX_MODEL(H,U_TAU,DZ,Z0,TMP_G,TMP_S,RHO,CP)

REAL(EB), INTENT(OUT) :: H ! heat transfer coefficient
REAL(EB), INTENT(IN) :: U_TAU,DZ,Z0,TMP_G,TMP_S,RHO,CP
REAL(EB), PARAMETER :: KAPPA=0.41_EB ! von Karman constant
REAL(EB) :: PSI_H,L,ZP,Q3S,A,B,C,DET

! References:
!
! Stoll, R., Porte-Agel, F. (2008) Large-Eddy Simulation of the Stable Atmospheric
! Boundary Layer using Dynamic Models with Different Averaging Schemes. Boundary-Layer
! Meteorology, 126:1-28.

ZP = 0.5_EB*DZ/Z0
PSI_H = 0._EB

IF (TMP_S<TMP_G) THEN ! stability correction needed

   ! Stoll and Porte-Agel Eq. (27) may be rewritten in quadratic form as follows:
   ! A*Q3S**2 + B*Q3S + C = 0

   A = 7.8_EB*KAPPA*GRAV
   B = -U_TAU**3*TMP_G*LOG(ZP)
   C = (TMP_S-TMP_G)*U_TAU**4*KAPPA*GRAV   ! negative, so DET should be positive
   DET = B**2-4._EB*A*C
   IF (DET>TWO_EPSILON_EB) THEN
      Q3S = (-B-SQRT(DET))/(2._EB*A)       ! take negative root, else Q3S>0, which we know is not correct here
      L = -U_TAU**3*TMP_G/(KAPPA*GRAV*Q3S) ! Eq. (28), Obukhov length
      PSI_H = -3.9_EB*DZ/MAX(L,MICRON)     ! Eq. (30), 3.9 = 7.8*0.5
   ELSE
      Q3S = 0._EB
      L = 0._EB
   ENDIF

ENDIF

H = RHO*CP*U_TAU*KAPPA/(LOG(ZP)-PSI_H)

END SUBROUTINE ABL_HEAT_FLUX_MODEL


REAL(EB) FUNCTION VELTAN2D(U_VELO,U_SURF,NN,DN,DIVU,GRADU,GRADP,TAU_IJ,DT,RRHO,MU,I_VEL)

REAL(EB), INTENT(IN) :: U_VELO(2),U_SURF(2),NN(2),DN,DIVU,GRADU(2,2),GRADP(2),TAU_IJ(2,2),DT,RRHO,MU
INTEGER, INTENT(IN) :: I_VEL
REAL(EB) :: C(2,2),SS(2),SLIP_COEF,ETA,AA,BB,U_STRM_0,DUMMY(3),&
            U_STRM,U_NORM,U_STRM_WALL,U_NORM_WALL,DPDS,DUSDS,DUSDN,TSN,RDN
INTEGER :: SUBIT

! Cartesian grid coordinate system orthonormal basis vectors
REAL(EB), DIMENSION(2), PARAMETER :: XX=(/1._EB, 0._EB/),YY=(/0._EB, 1._EB/)


! streamwise unit vector
SS = (/NN(2),-NN(1)/)

! direction cosines (see Pope, Eq. A.11)
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

RDN = 1._EB/DN

! ODE solution
IF (SIM_MODE==DNS_MODE) THEN
   ETA = U_NORM + RRHO*MU*RDN
   AA  = -(0.5_EB*DUSDS + TWTH*ETA*RDN)
   BB  = (TWTH*U_STRM_WALL*RDN + ONSI*DUSDN)*ETA - (U_NORM*0.5_EB*DUSDN + RRHO*( DPDS + TSN*0.5_EB*RDN ))
   !AA  = -0.5_EB*(DUSDS + ETA/DN)
   !BB  = 0.5_EB*US_WALL/DN*ETA - (UN*0.5_EB*DUSDN + RRHO*( DPDS + TSN/(2._EB*DN) ))
   U_STRM = ((AA*U_STRM + BB)*EXP(AA*DT) - BB)/AA
ELSE
   U_STRM_0 = U_STRM
   DO SUBIT=1,1
      CALL WALL_MODEL(SLIP_COEF,DUMMY(1),DUMMY(2),U_STRM-U_STRM_WALL,MU*RRHO,DN,0._EB)
      !IF (SLIP_COEF< -1._EB .OR. SLIP_COEF>-1._EB) THEN
      !   PRINT *,SUBIT,'WARNING: SLIP_COEF=',SLIP_COEF
      !ENDIF
      ETA = RRHO*(1-SLIP_COEF)*MU*0.5_EB*RDN**2
      AA  = -(0.5_EB*DUSDS + TWTH*U_NORM*RDN + ETA)
      BB  = ETA*U_STRM_WALL - (U_NORM*ONTH*DUSDN + RRHO*( DPDS + TSN*0.5_EB*RDN ))
      U_STRM = ((AA*U_STRM_0 + BB)*EXP(AA*DT) - BB)/AA
   ENDDO
ENDIF

! transform velocity back to Cartesian component I_VEL
VELTAN2D = C(I_VEL,1)*U_STRM + C(I_VEL,2)*U_NORM

END FUNCTION VELTAN2D


REAL(EB) FUNCTION VELTAN3D(U_VELO,U_SURF,NN,DN,DIVU,GRADU,GRADP,TAU_IJ,DT,RRHO,MU,I_VEL,ROUGHNESS,U_INT)
USE MATH_FUNCTIONS, ONLY: CROSS_PRODUCT

REAL(EB), INTENT(IN) :: U_VELO(3),U_SURF(3),NN(3),DN,DIVU,GRADU(3,3),GRADP(3),TAU_IJ(3,3),DT,RRHO,MU,ROUGHNESS,U_INT
INTEGER, INTENT(IN) :: I_VEL
REAL(EB) :: C(3,3),SS(3),PP(3),SLIP_COEF,ETA,AA,BB,U_STRM_0,U_RELA(3),DUMMY(3),&
            U_STRM,U_ORTH,U_NORM,DPDS,DUSDS,DUSDN,TSN,DUPDP,DUNDN,RDN
INTEGER :: SUBIT,I,J

! Cartesian grid coordinate system orthonormal basis vectors
REAL(EB), DIMENSION(3), PARAMETER :: E1=(/1._EB,0._EB,0._EB/),E2=(/0._EB,1._EB,0._EB/),E3=(/0._EB,0._EB,1._EB/)

! find a vector PP in the tangent plane of the surface and orthogonal to U_VELO-U_SURF
U_RELA = U_VELO-U_SURF
CALL CROSS_PRODUCT(PP,NN,U_RELA) ! PP = NN x U_RELA
IF (ABS(NORM2(PP))<=TWO_EPSILON_EB) THEN
   ! tangent vector is completely arbitrary, just perpendicular to NN
   IF (ABS(NN(1))>=TWO_EPSILON_EB .OR.  ABS(NN(2))>=TWO_EPSILON_EB) PP = (/NN(2),-NN(1),0._EB/)
   IF (ABS(NN(1))<=TWO_EPSILON_EB .AND. ABS(NN(2))<=TWO_EPSILON_EB) PP = (/NN(3),0._EB,-NN(1)/)
ENDIF
PP = PP/NORM2(PP) ! normalize to unit vector
CALL CROSS_PRODUCT(SS,PP,NN) ! define the streamwise unit vector SS

!! check unit normal vectors
!print *,DOT_PRODUCT(SS,SS) ! should be 1
!print *,DOT_PRODUCT(SS,PP) ! should be 0
!print *,DOT_PRODUCT(SS,NN) ! should be 0
!print *,DOT_PRODUCT(PP,PP) ! should be 1
!print *,DOT_PRODUCT(PP,NN) ! should be 0
!print *,DOT_PRODUCT(NN,NN) ! should be 1
!print *                    ! blank line

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
DUPDP = 0._EB
DUNDN = 0._EB
DUSDN = 0._EB
TSN = 0._EB
DO I=1,3
   DO J=1,3
      DUSDS = DUSDS + C(I,1)*C(J,1)*GRADU(I,J)
      DUPDP = DUPDP + C(I,2)*C(J,2)*GRADU(I,J)
      DUNDN = DUNDN + C(I,3)*C(J,3)*GRADU(I,J)

      DUSDN = DUSDN + C(I,1)*C(J,3)*GRADU(I,J)
      TSN   = TSN   + C(I,1)*C(J,3)*TAU_IJ(I,J)
   ENDDO
ENDDO

!! check trace of transformed tensor
!print *
!print *,DIVU,DUSDS+DUPDP+DUNDN
!print *,GRADU

! update boundary layer equations

! update wall-normal velocity
U_NORM = DN*(DIVU-0.5_EB*DUSDS)
RDN = 1._EB/DN
! ODE solution
IF (SIM_MODE==DNS_MODE) THEN
   ETA = U_NORM + RRHO*MU*RDN
   AA  = -(0.5_EB*DUSDS + TWTH*ETA*RDN)
   BB  = ONSI*DUSDN*ETA - (U_NORM*0.5_EB*DUSDN + RRHO*( DPDS + TSN*0.5_EB*RDN ))
   IF (ABS(AA)>=TWO_EPSILON_EB) THEN
      U_STRM = ((AA*U_STRM + BB)*EXP(AA*DT) - BB)/AA
   ELSE
      VELTAN3D = U_INT
      RETURN
   ENDIF
ELSE
   U_STRM_0 = U_STRM
   DO SUBIT=1,1
      CALL WALL_MODEL(SLIP_COEF,DUMMY(1),DUMMY(2),U_STRM,MU*RRHO,DN,ROUGHNESS)
      !IF (SLIP_COEF<-100._EB .OR. SLIP_COEF>100._EB) THEN
      !   PRINT *,SUBIT,'WARNING: SLIP_COEF=',SLIP_COEF
      !ENDIF
      ETA = RRHO*(1-SLIP_COEF)*MU*0.5_EB*RDN**2
      AA  = -(0.5_EB*DUSDS + TWTH*U_NORM*RDN + ETA)
      BB  = -(U_NORM*ONTH*DUSDN + RRHO*( DPDS + TSN*0.5_EB*RDN))
      !print *,MU*RRHO*DT/(DN**2)
      IF (ABS(AA)>=TWO_EPSILON_EB) THEN
         U_STRM = ((AA*U_STRM_0 + BB)*EXP(AA*DT) - BB)/AA
      ELSE
         VELTAN3D = U_INT
         RETURN
      ENDIF
   ENDDO
ENDIF

! transform velocity back to Cartesian component I_VEL
VELTAN3D = C(I_VEL,1)*U_STRM + C(I_VEL,3)*U_NORM + U_SURF(I_VEL)

END FUNCTION VELTAN3D


SUBROUTINE TAU_WALL_IJ(TAU_IJ,SS,U_VELO,U_SURF,NN,DN,DIVU,MU,RHO,ROUGHNESS)
USE MATH_FUNCTIONS, ONLY: CROSS_PRODUCT

REAL(EB), INTENT(OUT) :: TAU_IJ(3,3),SS(3)
REAL(EB), INTENT(IN) :: U_VELO(3),U_SURF(3),NN(3),DN,DIVU,MU,RHO,ROUGHNESS
REAL(EB) :: C(3,3),TT(3),U_RELA(3),SLIP_COEF,Y_PLUS,U_STRM,U_ORTH,U_NORM,TAUBAR_IJ(3,3),U_TAU
INTEGER :: K,L,M,N

! Cartesian grid coordinate system orthonormal basis vectors
REAL(EB), DIMENSION(3), PARAMETER :: E1=(/1._EB,0._EB,0._EB/),E2=(/0._EB,1._EB,0._EB/),E3=(/0._EB,0._EB,1._EB/)

! find a vector TT in the tangent plane of the surface and orthogonal to U_VELO-U_SURF
U_RELA = U_VELO-U_SURF
CALL CROSS_PRODUCT(TT,NN,U_RELA) ! TT = NN x U_RELA
IF (ABS(NORM2(TT))<=TWO_EPSILON_EB) THEN
   ! tangent vector is completely arbitrary, just perpendicular to NN
   IF (ABS(NN(1))>=TWO_EPSILON_EB .OR.  ABS(NN(2))>=TWO_EPSILON_EB) TT = (/NN(2),-NN(1),0._EB/)
   IF (ABS(NN(1))<=TWO_EPSILON_EB .AND. ABS(NN(2))<=TWO_EPSILON_EB) TT = (/NN(3),0._EB,-NN(1)/)
ENDIF
TT = TT/NORM2(TT) ! normalize to unit vector
CALL CROSS_PRODUCT(SS,TT,NN) ! define the streamwise unit vector SS

!! check unit normal vectors
!print *,DOT_PRODUCT(SS,SS) ! should be 1
!print *,DOT_PRODUCT(SS,TT) ! should be 0
!print *,DOT_PRODUCT(SS,NN) ! should be 0
!print *,DOT_PRODUCT(TT,TT) ! should be 1
!print *,DOT_PRODUCT(TT,NN) ! should be 0
!print *,DOT_PRODUCT(NN,NN) ! should be 1
!print *                    ! blank line

! directional cosines (see Pope, Eq. A.11)
C(1,1) = DOT_PRODUCT(E1,SS)
C(1,2) = DOT_PRODUCT(E1,TT)
C(1,3) = DOT_PRODUCT(E1,NN)
C(2,1) = DOT_PRODUCT(E2,SS)
C(2,2) = DOT_PRODUCT(E2,TT)
C(2,3) = DOT_PRODUCT(E2,NN)
C(3,1) = DOT_PRODUCT(E3,SS)
C(3,2) = DOT_PRODUCT(E3,TT)
C(3,3) = DOT_PRODUCT(E3,NN)

! transform velocity (see Pope, Eq. A.17)
U_STRM = C(1,1)*U_RELA(1) + C(2,1)*U_RELA(2) + C(3,1)*U_RELA(3)
U_ORTH = C(1,2)*U_RELA(1) + C(2,2)*U_RELA(2) + C(3,2)*U_RELA(3)
U_NORM = C(1,3)*U_RELA(1) + C(2,3)*U_RELA(2) + C(3,3)*U_RELA(3)

! ! check U_ORTH, should be zero
! print *, 'U_STRM: ',U_STRM
! print *, 'U_ORTH: ',U_ORTH
! print *, 'U_NORM: ',U_NORM

! in the streamwise coordinate system, the stress tensor simplifies to the symmetric tensor
! T = [0      0 T(1,3)]
!     [0      0      0]
!     [T(3,1) 0 T(3,3)]

TAUBAR_IJ      = 0._EB
TAUBAR_IJ(3,3) = -2._EB*MU*(U_NORM*2._EB/DN - ONTH*DIVU)

IF (SIM_MODE==DNS_MODE) THEN
   TAUBAR_IJ(1,3) = MU*U_STRM*2._EB/DN
ELSE
   CALL WALL_MODEL(SLIP_COEF,U_TAU,Y_PLUS,U_STRM,MU/RHO,DN,ROUGHNESS)
   TAUBAR_IJ(1,3) = RHO*U_TAU**2
ENDIF
TAUBAR_IJ(3,1) = TAUBAR_IJ(1,3)

! transform tensors (Pope A.23)
TAU_IJ = 0._EB
DO M=1,3
   DO N=1,3
      ! inner summation for component m,n ---------------------
      DO L=1,3
         DO K=1,3
            TAU_IJ(M,N) = TAU_IJ(M,N) + C(M,K)*C(N,L)*TAUBAR_IJ(K,L)
         ENDDO
      ENDDO
      !--------------------------------------------------------
   ENDDO
ENDDO

RETURN
END SUBROUTINE TAU_WALL_IJ


SUBROUTINE SYNTHETIC_EDDY_SETUP(NM)

INTEGER, INTENT(IN) :: NM
TYPE(VENTS_TYPE), POINTER :: VT=>NULL()
INTEGER :: NE,NV,IERROR
REAL(EB), POINTER, DIMENSION(:,:) :: A_IJ=>NULL(),R_IJ=>NULL()

IF (EVACUATION_ONLY(NM)) RETURN
VENT_LOOP: DO NV=1,MESHES(NM)%N_VENT
   VT => MESHES(NM)%VENTS(NV)
   IF (VT%N_EDDY==0) CYCLE VENT_LOOP

   SELECT CASE(ABS(VT%IOR))
      CASE(1)
         VT%X_EDDY_MIN = VT%X1-MAXVAL(VT%SIGMA_IJ(:,1))
         VT%X_EDDY_MAX = VT%X2+MAXVAL(VT%SIGMA_IJ(:,1))
         VT%Y_EDDY_MIN = VT%Y1
         VT%Y_EDDY_MAX = VT%Y2
         VT%Z_EDDY_MIN = VT%Z1
         VT%Z_EDDY_MAX = VT%Z2
      CASE(2)
         VT%X_EDDY_MIN = VT%X1
         VT%X_EDDY_MAX = VT%X2
         VT%Y_EDDY_MIN = VT%Y1-MAXVAL(VT%SIGMA_IJ(:,2))
         VT%Y_EDDY_MAX = VT%Y2+MAXVAL(VT%SIGMA_IJ(:,2))
         VT%Z_EDDY_MIN = VT%Z1
         VT%Z_EDDY_MAX = VT%Z2
      CASE(3)
         VT%X_EDDY_MIN = VT%X1
         VT%X_EDDY_MAX = VT%X2
         VT%Y_EDDY_MIN = VT%Y1
         VT%Y_EDDY_MAX = VT%Y2
         VT%Z_EDDY_MIN = VT%Z1-MAXVAL(VT%SIGMA_IJ(:,3))
         VT%Z_EDDY_MAX = VT%Z2+MAXVAL(VT%SIGMA_IJ(:,3))
   END SELECT

   VT%EDDY_BOX_VOLUME = (VT%X_EDDY_MAX-VT%X_EDDY_MIN)*(VT%Y_EDDY_MAX-VT%Y_EDDY_MIN)*(VT%Z_EDDY_MAX-VT%Z_EDDY_MIN)

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

   EDDY_LOOP: DO NE=1,VT%N_EDDY
      IERROR=1; CALL EDDY_POSITION(NE,NV,NM,IERROR)
      CALL EDDY_AMPLITUDE(NE,NV,NM)
   ENDDO EDDY_LOOP

ENDDO VENT_LOOP

END SUBROUTINE SYNTHETIC_EDDY_SETUP


SUBROUTINE SYNTHETIC_TURBULENCE(DT,T,NM)

USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
USE TRAN, ONLY: GET_IJK

REAL(EB), INTENT(IN) :: DT,T
INTEGER, INTENT(IN) :: NM
INTEGER :: NE,NV,II,JJ,KK,IERROR
TYPE(VENTS_TYPE), POINTER :: VT=>NULL()
TYPE(SURFACE_TYPE), POINTER :: SF=>NULL()
REAL(EB) :: XX,YY,ZZ,SHAPE_FACTOR,VOLUME_WEIGHTING_FACTOR(3),EDDY_VOLUME(3),PROFILE_FACTOR,RAMP_T,TSI
REAL(EB), PARAMETER :: EPSDX=1.E-10_EB
INTEGER, PARAMETER :: SHAPE_CODE=1 ! 1=tent, 2=tophat

! Reference:
!
! Nicolas Jarrin. Synthetic Inflow Boundary Conditions for the Numerical Simulation of Turbulence. PhD Thesis,
! The University of Manchester, 2008.
!
! See Chapter 4: The Synthetic Eddy Method

IF (EVACUATION_ONLY(NM)) RETURN
CALL POINT_TO_MESH(NM)

VENT_LOOP: DO NV=1,N_VENT
   VT => VENTS(NV)
   IF (VT%N_EDDY==0) CYCLE VENT_LOOP

   VT%U_EDDY = 0._EB
   VT%V_EDDY = 0._EB
   VT%W_EDDY = 0._EB
   SF => SURFACE(VT%SURF_INDEX)

   IF (ABS(SF%T_IGN-T_BEGIN)<=SPACING(SF%T_IGN) .AND. SF%RAMP_INDEX(TIME_VELO)>=1) THEN
      TSI = T
   ELSE
      TSI=T-SF%T_IGN
   ENDIF
   PROFILE_FACTOR = 1._EB
   RAMP_T = EVALUATE_RAMP(TSI,SF%TAU(TIME_VELO),SF%RAMP_INDEX(TIME_VELO))

   IOR_SELECT: SELECT CASE(ABS(VT%IOR))
      CASE(1)
         EDDY_LOOP_1: DO NE=1,VT%N_EDDY ! loop over eddies
            PROFILE_FACTOR=1._EB
            VT%X_EDDY(NE) = VT%X_EDDY(NE) - DT*SF%VEL*SIGN(1._EB,REAL(VT%IOR,EB))*PROFILE_FACTOR*RAMP_T
            VT%Y_EDDY(NE) = VT%Y_EDDY(NE) + DT*SF%VEL_T(1)*PROFILE_FACTOR*RAMP_T
            VT%Z_EDDY(NE) = VT%Z_EDDY(NE) + DT*SF%VEL_T(2)*PROFILE_FACTOR*RAMP_T
            IERROR=0;      CALL EDDY_POSITION(NE,NV,NM,IERROR)
            IF (IERROR==1) CALL EDDY_AMPLITUDE(NE,NV,NM)
            DO KK=VT%K1+1,VT%K2 ! this block can be made more efficient
               DO JJ=VT%J1+1,VT%J2
                  XX = (VT%X1  - VT%X_EDDY(NE))/VT%SIGMA_IJ(1,1)
                  YY = (YC(JJ) - VT%Y_EDDY(NE))/VT%SIGMA_IJ(1,2)
                  ZZ = (ZC(KK) - VT%Z_EDDY(NE))/VT%SIGMA_IJ(1,3)
                  SHAPE_FACTOR = SHAPE_FUNCTION(XX,SHAPE_CODE)*SHAPE_FUNCTION(YY,SHAPE_CODE)*SHAPE_FUNCTION(ZZ,SHAPE_CODE)
                  VT%U_EDDY(JJ,KK) = VT%U_EDDY(JJ,KK) + VT%CU_EDDY(NE)*SHAPE_FACTOR

                  XX = (VT%X1  - VT%X_EDDY(NE))/VT%SIGMA_IJ(2,1)
                  YY = (YC(JJ) - VT%Y_EDDY(NE))/VT%SIGMA_IJ(2,2)
                  ZZ = (ZC(KK) - VT%Z_EDDY(NE))/VT%SIGMA_IJ(2,3)
                  SHAPE_FACTOR = SHAPE_FUNCTION(XX,SHAPE_CODE)*SHAPE_FUNCTION(YY,SHAPE_CODE)*SHAPE_FUNCTION(ZZ,SHAPE_CODE)
                  VT%V_EDDY(JJ,KK) = VT%V_EDDY(JJ,KK) + VT%CV_EDDY(NE)*SHAPE_FACTOR

                  XX = (VT%X1  - VT%X_EDDY(NE))/VT%SIGMA_IJ(3,1)
                  YY = (YC(JJ) - VT%Y_EDDY(NE))/VT%SIGMA_IJ(3,2)
                  ZZ = (ZC(KK) - VT%Z_EDDY(NE))/VT%SIGMA_IJ(3,3)
                  SHAPE_FACTOR = SHAPE_FUNCTION(XX,SHAPE_CODE)*SHAPE_FUNCTION(YY,SHAPE_CODE)*SHAPE_FUNCTION(ZZ,SHAPE_CODE)
                  VT%W_EDDY(JJ,KK) = VT%W_EDDY(JJ,KK) + VT%CW_EDDY(NE)*SHAPE_FACTOR
               ENDDO
            ENDDO
         ENDDO EDDY_LOOP_1
      CASE(2)
         EDDY_LOOP_2: DO NE=1,VT%N_EDDY
            PROFILE_FACTOR=1._EB
            VT%X_EDDY(NE) = VT%X_EDDY(NE) + DT*SF%VEL_T(1)*PROFILE_FACTOR*RAMP_T
            VT%Y_EDDY(NE) = VT%Y_EDDY(NE) - DT*SF%VEL*SIGN(1._EB,REAL(VT%IOR,EB))*PROFILE_FACTOR*RAMP_T
            VT%Z_EDDY(NE) = VT%Z_EDDY(NE) + DT*SF%VEL_T(2)*PROFILE_FACTOR*RAMP_T
            IERROR=0;      CALL EDDY_POSITION(NE,NV,NM,IERROR)
            IF (IERROR==1) CALL EDDY_AMPLITUDE(NE,NV,NM)
            DO KK=VT%K1+1,VT%K2
               DO II=VT%I1+1,VT%I2
                  XX = (XC(II) - VT%X_EDDY(NE))/VT%SIGMA_IJ(1,1)
                  YY = (VT%Y1  - VT%Y_EDDY(NE))/VT%SIGMA_IJ(1,2)
                  ZZ = (ZC(KK) - VT%Z_EDDY(NE))/VT%SIGMA_IJ(1,3)
                  SHAPE_FACTOR = SHAPE_FUNCTION(XX,SHAPE_CODE)*SHAPE_FUNCTION(YY,SHAPE_CODE)*SHAPE_FUNCTION(ZZ,SHAPE_CODE)
                  VT%U_EDDY(II,KK) = VT%U_EDDY(II,KK) + VT%CU_EDDY(NE)*SHAPE_FACTOR

                  XX = (XC(II) - VT%X_EDDY(NE))/VT%SIGMA_IJ(2,1)
                  YY = (VT%Y1  - VT%Y_EDDY(NE))/VT%SIGMA_IJ(2,2)
                  ZZ = (ZC(KK) - VT%Z_EDDY(NE))/VT%SIGMA_IJ(2,3)
                  SHAPE_FACTOR = SHAPE_FUNCTION(XX,SHAPE_CODE)*SHAPE_FUNCTION(YY,SHAPE_CODE)*SHAPE_FUNCTION(ZZ,SHAPE_CODE)
                  VT%V_EDDY(II,KK) = VT%V_EDDY(II,KK) + VT%CV_EDDY(NE)*SHAPE_FACTOR

                  XX = (XC(II) - VT%X_EDDY(NE))/VT%SIGMA_IJ(3,1)
                  YY = (VT%Y1  - VT%Y_EDDY(NE))/VT%SIGMA_IJ(3,2)
                  ZZ = (ZC(KK) - VT%Z_EDDY(NE))/VT%SIGMA_IJ(3,3)
                  SHAPE_FACTOR = SHAPE_FUNCTION(XX,SHAPE_CODE)*SHAPE_FUNCTION(YY,SHAPE_CODE)*SHAPE_FUNCTION(ZZ,SHAPE_CODE)
                  VT%W_EDDY(II,KK) = VT%W_EDDY(II,KK) + VT%CW_EDDY(NE)*SHAPE_FACTOR
               ENDDO
            ENDDO
         ENDDO EDDY_LOOP_2
      CASE(3)
         EDDY_LOOP_3: DO NE=1,VT%N_EDDY
            PROFILE_FACTOR=1._EB
            VT%X_EDDY(NE) = VT%X_EDDY(NE) + DT*SF%VEL_T(1)*PROFILE_FACTOR*RAMP_T
            VT%Y_EDDY(NE) = VT%Y_EDDY(NE) + DT*SF%VEL_T(2)*PROFILE_FACTOR*RAMP_T
            VT%Z_EDDY(NE) = VT%Z_EDDY(NE) - DT*SF%VEL*SIGN(1._EB,REAL(VT%IOR,EB))*PROFILE_FACTOR*RAMP_T
            IERROR=0;      CALL EDDY_POSITION(NE,NV,NM,IERROR)
            IF (IERROR==1) CALL EDDY_AMPLITUDE(NE,NV,NM)
            DO JJ=VT%J1+1,VT%J2
               DO II=VT%I1+1,VT%I2
                  XX = (XC(II) - VT%X_EDDY(NE))/VT%SIGMA_IJ(1,1)
                  YY = (YC(JJ) - VT%Y_EDDY(NE))/VT%SIGMA_IJ(1,2)
                  ZZ = (VT%Z1  - VT%Z_EDDY(NE))/VT%SIGMA_IJ(1,3)
                  SHAPE_FACTOR = SHAPE_FUNCTION(XX,SHAPE_CODE)*SHAPE_FUNCTION(YY,SHAPE_CODE)*SHAPE_FUNCTION(ZZ,SHAPE_CODE)
                  VT%U_EDDY(II,JJ) = VT%U_EDDY(II,JJ) + VT%CU_EDDY(NE)*SHAPE_FACTOR

                  XX = (XC(II) - VT%X_EDDY(NE))/VT%SIGMA_IJ(2,1)
                  YY = (YC(JJ) - VT%Y_EDDY(NE))/VT%SIGMA_IJ(2,2)
                  ZZ = (VT%Z1  - VT%Z_EDDY(NE))/VT%SIGMA_IJ(2,3)
                  SHAPE_FACTOR = SHAPE_FUNCTION(XX,SHAPE_CODE)*SHAPE_FUNCTION(YY,SHAPE_CODE)*SHAPE_FUNCTION(ZZ,SHAPE_CODE)
                  VT%V_EDDY(II,JJ) = VT%V_EDDY(II,JJ) + VT%CV_EDDY(NE)*SHAPE_FACTOR

                  XX = (XC(II) - VT%X_EDDY(NE))/VT%SIGMA_IJ(3,1)
                  YY = (YC(JJ) - VT%Y_EDDY(NE))/VT%SIGMA_IJ(3,2)
                  ZZ = (VT%Z1  - VT%Z_EDDY(NE))/VT%SIGMA_IJ(3,3)
                  SHAPE_FACTOR = SHAPE_FUNCTION(XX,SHAPE_CODE)*SHAPE_FUNCTION(YY,SHAPE_CODE)*SHAPE_FUNCTION(ZZ,SHAPE_CODE)
                  VT%W_EDDY(II,JJ) = VT%W_EDDY(II,JJ) + VT%CW_EDDY(NE)*SHAPE_FACTOR
               ENDDO
            ENDDO
         ENDDO EDDY_LOOP_3
   END SELECT IOR_SELECT

   EDDY_VOLUME(1) = VT%SIGMA_IJ(1,1)*VT%SIGMA_IJ(1,2)*VT%SIGMA_IJ(1,3)
   EDDY_VOLUME(2) = VT%SIGMA_IJ(2,1)*VT%SIGMA_IJ(2,2)*VT%SIGMA_IJ(2,3)
   EDDY_VOLUME(3) = VT%SIGMA_IJ(3,1)*VT%SIGMA_IJ(3,2)*VT%SIGMA_IJ(3,3)

   VOLUME_WEIGHTING_FACTOR(1) = MIN(1._EB,SQRT(VT%EDDY_BOX_VOLUME/REAL(VT%N_EDDY,EB)/EDDY_VOLUME(1)))
   VOLUME_WEIGHTING_FACTOR(2) = MIN(1._EB,SQRT(VT%EDDY_BOX_VOLUME/REAL(VT%N_EDDY,EB)/EDDY_VOLUME(2)))
   VOLUME_WEIGHTING_FACTOR(3) = MIN(1._EB,SQRT(VT%EDDY_BOX_VOLUME/REAL(VT%N_EDDY,EB)/EDDY_VOLUME(3)))

   ! note: EDDY_VOLUME included in SQRT based on Jung-il Choi write up.

   VT%U_EDDY = VT%U_EDDY*VOLUME_WEIGHTING_FACTOR(1)
   VT%V_EDDY = VT%V_EDDY*VOLUME_WEIGHTING_FACTOR(2)
   VT%W_EDDY = VT%W_EDDY*VOLUME_WEIGHTING_FACTOR(3)

   ! subtract mean from normal components so that fluctuations do not affect global volume flow

   SELECT CASE (ABS(VT%IOR))
      CASE(1)
         VT%U_EDDY = VT%U_EDDY - SUM(VT%U_EDDY)/SIZE(VT%U_EDDY)
      CASE(2)
         VT%V_EDDY = VT%V_EDDY - SUM(VT%V_EDDY)/SIZE(VT%V_EDDY)
      CASE(3)
         VT%W_EDDY = VT%W_EDDY - SUM(VT%W_EDDY)/SIZE(VT%W_EDDY)
   END SELECT


ENDDO VENT_LOOP

END SUBROUTINE SYNTHETIC_TURBULENCE


SUBROUTINE EDDY_POSITION(NE,NV,NM,IERROR)

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

INTEGER, INTENT(IN) :: NE,NV,NM
REAL(EB) :: EPS_EDDY(3)
REAL(EB) :: RN
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


SUBROUTINE TENSOR_DIFFUSIVITY_MODEL(NM,N)

USE PHYSICAL_FUNCTIONS, ONLY: LES_FILTER_WIDTH_FUNCTION

INTEGER, INTENT(IN) :: NM,N
INTEGER :: I,J,K
REAL(EB) :: DELTA,DRHOZDX,DRHOZDY,DRHOZDZ,DUDX,DUDY,DUDZ,DVDX,DVDY,DVDZ,DWDX,DWDY,DWDZ,DYN1,DZN1,OO12=1._EB/12._EB
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZP=>NULL()
REAL(EB), POINTER, DIMENSION(:,:,:) :: RHOP=>NULL(),RHO_D_DZDX=>NULL(),RHO_D_DZDY=>NULL(),RHO_D_DZDZ=>NULL(),&
                                       UU=>NULL(),VV=>NULL(),WW=>NULL()

IF (EVACUATION_ONLY(NM)) RETURN
CALL POINT_TO_MESH(NM)

! SGS scalar flux
! CAUTION: The flux arrays must point to the same work arrays used in DIVERGENCE_PART_1
RHO_D_DZDX=>WORK1
RHO_D_DZDY=>WORK2
RHO_D_DZDZ=>WORK3

IF (PREDICTOR) THEN
   UU=>U
   VV=>V
   WW=>W
   RHOP=>RHOS
   ZZP=>ZZS
ELSE
   UU=>US
   VV=>VS
   WW=>WS
   RHOP=>RHO
   ZZP=>ZZ
ENDIF

DO K=1,KBAR
   DZN1 = 1._EB/(DZN(K-1)+DZN(K))
   DO J=1,JBAR
      DYN1 = 1._EB/(DYN(J-1)+DYN(J))
      DO I=0,IBAR

         DELTA = LES_FILTER_WIDTH_FUNCTION(DXN(I),DY(J),DZ(K))

         DUDX = (UU(I+1,J,K)-UU(I-1,J,K))/(DX(I)+DX(I+1))
         DUDY = (UU(I,J+1,K)-UU(I,J-1,K))*DYN1
         DUDZ = (UU(I,J,K+1)-UU(I,J,K-1))*DZN1

         DRHOZDX = RDXN(I)*(RHOP(I+1,J,K)*ZZP(I+1,J,K,N)-RHOP(I,J,K)*ZZP(I,J,K,N))

         DRHOZDY = 0.25_EB*RDY(J)*( RHOP(I,J+1,K)*ZZP(I,J+1,K,N) + RHOP(I+1,J+1,K)*ZZP(I+1,J+1,K,N) &
                                  - RHOP(I,J-1,K)*ZZP(I,J-1,K,N) - RHOP(I+1,J-1,K)*ZZP(I+1,J-1,K,N) )

         DRHOZDZ = 0.25_EB*RDZ(K)*( RHOP(I,J,K+1)*ZZP(I,J,K+1,N) + RHOP(I+1,J,K+1)*ZZP(I+1,J,K+1,N) &
                                  - RHOP(I,J,K-1)*ZZP(I,J,K-1,N) - RHOP(I+1,J,K-1)*ZZP(I+1,J,K-1,N) )

         RHO_D_DZDX(I,J,K) = RHO_D_DZDX(I,J,K) - DELTA**2*OO12*(DUDX*DRHOZDX + DUDY*DRHOZDY + DUDZ*DRHOZDZ)

      ENDDO
   ENDDO
ENDDO

DO K=1,KBAR
   DZN1 = 1._EB/(DZN(K-1)+DZN(K))
   DO J=0,JBAR
      DYN1 = 1._EB/(DY(J)+DY(J+1))
      DO I=1,IBAR

         DELTA = LES_FILTER_WIDTH_FUNCTION(DX(I),DYN(J),DZ(K))

         DVDX = (VV(I+1,J,K)-VV(I-1,J,K))/(DXN(I-1)+DXN(I))
         DVDY = (VV(I,J+1,K)-VV(I,J-1,K))*DYN1
         DVDZ = (VV(I,J,K+1)-VV(I,J,K-1))*DZN1

         DRHOZDX = 0.25_EB*RDX(I)*( RHOP(I+1,J,K)*ZZP(I+1,J,K,N) + RHOP(I+1,J+1,K)*ZZP(I+1,J+1,K,N) &
                                  - RHOP(I-1,J,K)*ZZP(I-1,J,K,N) - RHOP(I-1,J+1,K)*ZZP(I-1,J+1,K,N) )

         DRHOZDY = RDYN(J)*(RHOP(I,J+1,K)*ZZP(I,J+1,K,N)-RHOP(I,J,K)*ZZP(I,J,K,N))

         DRHOZDZ = 0.25_EB*RDZ(K)*( RHOP(I,J,K+1)*ZZP(I,J,K+1,N) + RHOP(I,J+1,K+1)*ZZP(I,J+1,K+1,N) &
                                  - RHOP(I,J,K-1)*ZZP(I,J,K-1,N) - RHOP(I,J+1,K-1)*ZZP(I,J+1,K-1,N) )

         RHO_D_DZDY(I,J,K) = RHO_D_DZDY(I,J,K) - DELTA**2*OO12*(DVDX*DRHOZDX + DVDY*DRHOZDY + DVDZ*DRHOZDZ)

      ENDDO
   ENDDO
ENDDO

DO K=0,KBAR
   DZN1 = 1._EB/(DZ(K)+DZ(K+1))
   DO J=1,JBAR
      DYN1 = 1._EB/(DYN(J-1)+DYN(J))
      DO I=1,IBAR

         DELTA = LES_FILTER_WIDTH_FUNCTION(DX(I),DY(J),DZN(K))

         DWDX = (WW(I+1,J,K)-WW(I-1,J,K))/(DXN(I-1)+DXN(I))
         DWDY = (WW(I,J+1,K)-WW(I,J-1,K))*DYN1
         DWDZ = (WW(I,J,K+1)-WW(I,J,K-1))*DZN1

         DRHOZDX = 0.25_EB*RDX(I)*( RHOP(I+1,J,K)*ZZP(I+1,J,K,N) + RHOP(I+1,J,K+1)*ZZP(I+1,J,K+1,N) &
                                  - RHOP(I-1,J,K)*ZZP(I-1,J,K,N) - RHOP(I-1,J,K+1)*ZZP(I-1,J,K+1,N) )

         DRHOZDY = 0.25_EB*RDY(J)*( RHOP(I,J+1,K)*ZZP(I,J+1,K,N) + RHOP(I,J+1,K+1)*ZZP(I,J+1,K+1,N) &
                                  - RHOP(I,J-1,K)*ZZP(I,J-1,K,N) - RHOP(I,J-1,K+1)*ZZP(I,J-1,K+1,N) )

         DRHOZDZ = RDZN(K)*(RHOP(I,J,K+1)*ZZP(I,J,K+1,N)-RHOP(I,J,K)*ZZP(I,J,K,N))

         RHO_D_DZDZ(I,J,K) = RHO_D_DZDZ(I,J,K) - DELTA**2*OO12*(DWDX*DRHOZDX + DWDY*DRHOZDY + DWDZ*DRHOZDZ)

      ENDDO
   ENDDO
ENDDO

END SUBROUTINE TENSOR_DIFFUSIVITY_MODEL


SUBROUTINE SANDIA_DAT(NM,FN_ISO)
IMPLICIT NONE

! This routine reads the file 'iso_ini.dat', which is generated by turb_init.
! This exe generates a random velocity field with a spectrum that matches the
! Comte-Bellot/Corrsin 1971 experimental data.

REAL(EB) :: MEANU,MEANV,MEANW,DUMMY(3)
INTEGER :: I,J,K,II,JJ,KK,FILE_NUM,IM,IW,IOR,IIO,JJO,KKO,NX,NX_BLOCK,I_MESH,J_MESH,K_MESH,MX
INTEGER, INTENT(IN) :: NM
CHARACTER(80), INTENT(IN) :: FN_ISO
TYPE (MESH_TYPE), POINTER :: M=>NULL()
TYPE (MESH_TYPE), POINTER :: M2=>NULL()
TYPE (WALL_TYPE), POINTER :: WC=>NULL()
TYPE (EXTERNAL_WALL_TYPE), POINTER :: EWC=>NULL()

IF (NM>1) RETURN

M=>MESHES(NM)
NX = NINT((XF_MAX-XS_MIN)/M%DX(1))
NX_BLOCK = M%IBAR
MX = NX/NX_BLOCK

! initialize velocity arrays
DO IM = 1,NMESHES
   M=>MESHES(IM)
   M%U=0._EB
   M%V=0._EB
   M%W=0._EB
ENDDO

! zero out mean values
MEANU=0._EB
MEANV=0._EB
MEANW=0._EB

FILE_NUM = GET_FILE_NUMBER()

OPEN (UNIT=FILE_NUM,FILE=FN_ISO,FORM='formatted',STATUS='old')

READ (FILE_NUM,*) II, JJ, KK  ! reads number of points in each direction

IF (II/=NX) CALL SHUTDOWN('ERROR: wrong iso_ini.dat file')

! read physical dimensions

READ (FILE_NUM,*) DUMMY(1), DUMMY(2), DUMMY(3)
READ (FILE_NUM,*) DUMMY(1), DUMMY(2), DUMMY(3)

DO KK=1,NX
   K_MESH = CEILING(REAL(KK,EB)/REAL(NX_BLOCK,EB)) ! "k" position of mesh block in global domain
   K = KK-(K_MESH-1)*NX_BLOCK ! local k in mesh block
   DO JJ=1,NX
      J_MESH = CEILING(REAL(JJ,EB)/REAL(NX_BLOCK,EB))
      J = JJ-(J_MESH-1)*NX_BLOCK
      DO II=1,NX
         I_MESH = CEILING(REAL(II,EB)/REAL(NX_BLOCK,EB))
         I = II-(I_MESH-1)*NX_BLOCK

         ! lexicographic ordering of mesh blocks
         IM = (K_MESH-1)*MX*MX + (J_MESH-1)*MX + I_MESH
         M=>MESHES(IM)

         READ (FILE_NUM,*) M%U(I,J,K), M%V(I,J,K), M%W(I,J,K)

         MEANU = MEANU + M%U(I,J,K)
         MEANV = MEANV + M%V(I,J,K)
         MEANW = MEANW + M%W(I,J,K)

      ENDDO
   ENDDO
ENDDO

! subtract mean

MEANU = MEANU/NX**3
MEANV = MEANV/NX**3
MEANW = MEANW/NX**3
DO IM = 1,NMESHES
   M=>MESHES(IM)
   M%U=M%U-MEANU
   M%V=M%V-MEANV
   M%W=M%W-MEANW
ENDDO


! handle boundaries and I=0,J=0,K=0 cells for each mesh

MESH_LOOP_2: DO IM = 1,NMESHES

   M=>MESHES(IM)

   WALL_LOOP: DO IW=1,M%N_EXTERNAL_WALL_CELLS

      WC=>M%WALL(IW)
      EWC=>M%EXTERNAL_WALL(IW)
      IOR = WC%ONE_D%IOR
      II = WC%ONE_D%II
      JJ = WC%ONE_D%JJ
      KK = WC%ONE_D%KK

      ! assumes no refinement
      IIO = EWC%IIO_MIN
      JJO = EWC%JJO_MIN
      KKO = EWC%KKO_MIN

      M2=>MESHES(EXTERNAL_WALL(IW)%NOM)

      ! for I=0, IOR=+1, etc.
      SELECT CASE(IOR)
         CASE(1)
            M%U(II,JJ,KK) = M2%U(IIO,JJO,KKO)
         CASE(2)
            M%V(II,JJ,KK) = M2%V(IIO,JJO,KKO)
         CASE(3)
            M%W(II,JJ,KK) = M2%W(IIO,JJO,KKO)
      END SELECT

   ENDDO WALL_LOOP

ENDDO MESH_LOOP_2

END SUBROUTINE SANDIA_DAT


SUBROUTINE init_spectral_data(NM)
USE MEMORY_FUNCTIONS, ONLY: ChkMemErr
IMPLICIT NONE
INTEGER, INTENT(IN) :: NM
INTEGER :: IZERO, NX
TYPE (MESH_TYPE), POINTER :: M

IF (NM>1) RETURN

M => MESHES(NM)
NX = NINT((XF_MAX-XS_MIN)/M%DX(1))

! real work arrays
ALLOCATE(M%PWORK1(NX,NX,NX),STAT=IZERO)
CALL ChkMemErr('init_spectral_data','PWORK1',IZERO)
M%PWORK1 = 0._EB
ALLOCATE(M%PWORK2(NX,NX,NX),STAT=IZERO)
CALL ChkMemErr('init_spectral_data','PWORK2',IZERO)
M%PWORK2 = 0._EB
ALLOCATE(M%PWORK3(NX,NX,NX),STAT=IZERO)
CALL ChkMemErr('init_spectral_data','PWORK3',IZERO)
M%PWORK3 = 0._EB

! complex work arrays
ALLOCATE(M%PWORK5(NX,NX,NX),STAT=IZERO)
CALL ChkMemErr('init_spectral_data','PWORK5',IZERO)
M%PWORK5 = 0._EB
ALLOCATE(M%PWORK6(NX,NX,NX),STAT=IZERO)
CALL ChkMemErr('init_spectral_data','PWORK6',IZERO)
M%PWORK6 = 0._EB
ALLOCATE(M%PWORK7(NX,NX,NX),STAT=IZERO)
CALL ChkMemErr('init_spectral_data','PWORK7',IZERO)
M%PWORK7 = 0._EB
ALLOCATE(M%PWORK8(NX,NX,NX),STAT=IZERO)
CALL ChkMemErr('init_spectral_data','PWORK8',IZERO)
M%PWORK8 = 0._EB

END SUBROUTINE init_spectral_data


SUBROUTINE spectral_output(NM,FN_SPEC)
IMPLICIT NONE
INTEGER, INTENT(IN) :: NM
CHARACTER(80), INTENT(IN) :: FN_SPEC
INTEGER :: NN(3),I,J,K,IM,II,JJ,KK
REAL(EB),     POINTER, DIMENSION(:,:,:) :: UU=>NULL(),VV=>NULL(),WW=>NULL()
COMPLEX(DPC), POINTER, DIMENSION(:,:,:) :: UUHT=>NULL(),VVHT=>NULL(),WWHT=>NULL(),KKHT=>NULL()
TYPE (MESH_TYPE), POINTER :: MM,M

IF (NM>1) RETURN

MM => MESHES(NM)
NN = NINT((XF_MAX-XS_MIN)/MM%DX(1)) ! obviously, assumes uniform grid spacing

UU => MM%PWORK1
VV => MM%PWORK2
WW => MM%PWORK3
UUHT => MM%PWORK5
VVHT => MM%PWORK6
WWHT => MM%PWORK7
KKHT => MM%PWORK8

MESH_LOOP: DO IM = 1,NMESHES

   M=>MESHES(IM)

   DO K=1,M%KBAR
      KK = NINT(M%ZS/M%DZ(1)) + K
      DO J=1,M%JBAR
         JJ = NINT(M%YS/M%DY(1)) + J
         DO I=1,M%IBAR
            II = NINT(M%XS/M%DX(1)) + I

            UU(II,JJ,KK) = M%U(I,J,K)
            VV(II,JJ,KK) = M%V(I,J,K)
            WW(II,JJ,KK) = M%W(I,J,K)

         ENDDO
      ENDDO
   ENDDO

ENDDO MESH_LOOP

! take fourier transform of velocities in 3d...

call fft3d_f90(UU, UUHT, NN)
call fft3d_f90(VV, VVHT, NN)
call fft3d_f90(WW, WWHT, NN)

! calc the spectral kinetic energy

call complex_tke_f90(KKHT, UUHT, VVHT, WWHT, NN(1))

! total up the spectral energy for each mode and integrate over
! the resolved modes...

call spectrum_f90(KKHT, NN(1), XF_MAX-XS_MIN, FN_SPEC)

END SUBROUTINE spectral_output


SUBROUTINE sandia_out(NM)
IMPLICIT NONE

INTEGER, INTENT(IN) :: NM
INTEGER :: NN(3),I,J,K,IM,II,JJ,KK,file_num
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU=>NULL(),VV=>NULL(),WW=>NULL()
TYPE (MESH_TYPE), POINTER :: MM,M

! This subroutine writes out the velocity, relative pressure, and
! turbulent kinetic energy to the file 'ini_salsa.dat'.  This is
! then used as input to the 'turb_init' program.

IF (NM>1) RETURN

MM => MESHES(NM)
NN = NINT((XF_MAX-XS_MIN)/MM%DX(1)) ! obviously, assumes uniform grid spacing

UU => MM%PWORK1
VV => MM%PWORK2
WW => MM%PWORK3

MESH_LOOP: DO IM = 1,NMESHES

   M=>MESHES(IM)

   DO K=1,M%KBAR
      KK = NINT(M%ZS/M%DZ(1)) + K
      DO J=1,M%JBAR
         JJ = NINT(M%YS/M%DY(1)) + J
         DO I=1,M%IBAR
            II = NINT(M%XS/M%DX(1)) + I

            UU(II,JJ,KK) = M%U(I,J,K)
            VV(II,JJ,KK) = M%V(I,J,K)
            WW(II,JJ,KK) = M%W(I,J,K)

         ENDDO
      ENDDO
   ENDDO

ENDDO MESH_LOOP

file_num = GET_FILE_NUMBER()
OPEN (unit=file_num, file='ini_salsa.dat', form='formatted', status='unknown', position='rewind')

WRITE (file_num,997) nn(1),nn(2),nn(3)
WRITE (file_num,998) 0._EB,0._EB,0._EB
WRITE (file_num,998) 2._EB*pi,2._EB*pi,2._EB*pi

DO k = 1,nn(3)
   DO j = 1,nn(2)
      DO i = 1,nn(1)

         WRITE (file_num,999) UU(I,J,K),VV(I,J,K),WW(I,J,K),0._EB,0._EB

      END DO
   END DO
END DO

CLOSE (unit=file_num)

997   FORMAT(3(i6,8x))
998   FORMAT(3(f12.6,2x))
999   FORMAT(5(f12.6,2x))

END SUBROUTINE sandia_out


SUBROUTINE complex_tke_f90(tkeht, upht, vpht, wpht, n)
IMPLICIT NONE
INTEGER, INTENT(IN) :: n
COMPLEX(DPC), INTENT(OUT) :: tkeht(n,n,n)
COMPLEX(DPC), INTENT(IN) :: upht(n,n,n),vpht(n,n,n),wpht(n,n,n)
INTEGER i,j,k

do k = 1,n
   do j = 1,n
      do i = 1,n
         tkeht(i,j,k) = 0.5_EB*( upht(i,j,k)*conjg(upht(i,j,k)) + &
                                 vpht(i,j,k)*conjg(vpht(i,j,k)) + &
                                 wpht(i,j,k)*conjg(wpht(i,j,k)) )
      end do
   end do
end do

END SUBROUTINE complex_tke_f90


SUBROUTINE fft3d_f90(v, vht, nn)
IMPLICIT NONE
INTEGER, INTENT(IN) :: nn(3)
REAL(EB), INTENT(IN) :: v(nn(1),nn(2),nn(3))
COMPLEX(DPC), INTENT(INOUT) :: vht(nn(1),nn(2),nn(3))

! This routine performs an FFT on the real array v and places the
! result in the complex array vht.

INTEGER :: i,j,k
REAL(EB) :: z2(2*nn(1))
!REAL(EB) :: tke

! convert v to spectral space

!tke = 0._EB

do k = 1,nn(3)
   do j = 1,nn(2)
      do i = 1,nn(1)
         vht(i,j,k) = cmplx(v(i,j,k),0.0_EB,kind=DPC)
         !tke = tke + 0.5_EB*v(i,j,k)**2
      end do
   end do
end do

!tke = tke/real(nn(1)*nn(2)*nn(3),EB)

!print*, ' fft3d internal check-'
!print*, ' tkeave (physical) =', tke

!do k=1,nn(3)
!   do j=1,nn(2)
!      do i=1,nn(1)
!         vht(i,j,k) = cmplx(float(i)+k*100,float(j)+.5+k*100)
!      enddo
!   enddo
!enddo

call fourier_f2003(vht,nn,3,-1,0,z2)
!call four3(vht,-1)

do k = 1,nn(3)
   do j = 1,nn(2)
      do i = 1,nn(1)
         vht(i,j,k) = vht(i,j,k)/nn(1)/nn(2)/nn(3)
      end do
   end do
end do

!tke = 0._EB
!do k = 1,nn(3)
!   do j = 1,nn(2)
!      do i = 1,nn(1)
!         tke = tke + 0.5_EB*real(vht(i,j,k)*conjg(vht(i,j,k)),kind=EB)
!      end do
!   end do
!end do

!print*, ' tkeave (spectral) =', tke

END SUBROUTINE fft3d_f90


SUBROUTINE spectrum_f90(vht, n, Lm, filename)
IMPLICIT NONE
INTEGER, INTENT(IN) :: n
CHARACTER(80), INTENT(IN) :: filename
COMPLEX(DPC), INTENT(IN) :: vht(n,n,n)
REAL(EB), INTENT(IN) :: Lm

! This routine is copied from SNL and is intended to compute the
! kinetic energy spectrum in wavenumber space.

INTEGER :: kmax, kx, ky, kz, k, ksum, file_num
INTEGER :: num(0:n)

REAL(EB) :: rk, temp, vsum, rkx, rky, rkz, etot
REAL(EB) :: vt(0:n)

! for dimensional wavenumbers
REAL(EB) :: wn(0:n)
REAL(EB) :: L, k0

L = Lm

k0 = 2*PI/L
kmax = n/2

wn(0) = 0._EB
do k = 1,n
   wn(k) = k0*k
end do

do k = 0,n
   vt(k) = 0._EB
   num(k) = 0
end do

etot = 0._EB

do kx = 1,n

   rkx = real(kx-1,kind=EB)
   if (rkx .gt. kmax) then
      rkx = n - rkx
   end if

   do ky = 1,n

      rky = real(ky-1,kind=EB)
      if (rky .gt. kmax) then
         rky = n - rky
      end if

      do kz = 1,n

         rkz = real(kz-1,kind=EB)
         if (rkz .gt. kmax) then
            rkz = n - rkz
         end if

         rk     = sqrt(rkx*rkx + rky*rky + rkz*rkz)
         k      = nint(rk)

         num(k) = num(k) + 1
         temp   = real(vht(kx,ky,kz)*conjg(vht(kx,ky,kz)),kind=EB)
         etot   = etot + sqrt(temp)
         vt(k)  = vt(k) + sqrt(temp)*(L/(2*PI))

      end do
   end do
end do

write(6,*) ' '
write(6,*) ' Spectrum Internal Check-'
write(6,*) ' Total Energy (-) in 3D field = ', etot
write(6,*) ' k(-), Num(k), k(1/cm), E(cm3/s2) '
ksum = 0
vsum = 0._EB
do k = 0,n
   write(6,*) k, num(k), wn(k), vt(k)
   ksum = ksum + num(k)
   vsum = vsum + vt(k)
end do

write(6,*) ' ksum: ', ksum
write(6,*) ' Total Energy (-) in spectrum: ', vsum
write(6,*) ' '

! write the spectral data to a file
! units are (1/cm) for wavenumber and (cm3/s2) for energy
! this matches the Comte-Bellot/Corrsin units

file_num = GET_FILE_NUMBER()
open (unit=file_num, file=filename, status='unknown', form='formatted')

do k = 0,n
   write (file_num,*) wn(k),',',vt(k)
end do

close (unit = file_num)

END SUBROUTINE spectrum_f90


SUBROUTINE fourier_f2003(data3,nn,ndim,isign,iform,work)
IMPLICIT NONE
! Converted to F90 10/8/2008 by RJM
!
! Converted to F2003 11/18/2014 by RJM
!
!     the cooley-tukey fast fourier transform in usasi basic fortran
!     transform(j1,j2,,,,) = sum(data(i1,i2,,,,)*w1**((i2-1)*(j2-1))
!                                 *w2**((i2-1)*(j2-1))*,,,),
!     where i1 and j1 run from 1 to nn(1) and w1=exp(isign*2*pi=
!     sqrt(-1)/nn(1)), etc.  there is no limit on the dimensionality
!     (number of subscripts) of the data array.  if an inverse
!     transform (isign=+1) is performed upon an array of transformed
!     (isign=-1) data, the original data will reappear.
!     multiplied by nn(1)*nn(2)*,,,  the array of input data must be
!     in complex format.  however, if all imaginary parts are zero (i.e.
!     the data are disguised real) running time is cut up to forty per-
!     cent.  (for fastest transform of real data, nn(1) should be even.)
!     the transform values are always complex and are returned in the
!     original array of data, replacing the input data.  the length
!     of each dimension of the data array may be any integer.  the
!     program runs faster on composite integers than on primes, and is
!     particularly fast on numbers rich in factors of two.
!
!     timing is in fact given by the following formula.  let ntot be the
!     total number of points (real or complex) in the data array, that
!     is, ntot=nn(1)*nn(2)*...  decompose ntot into its prime factors,
!     such as 2**k2 * 3**k3 * 5**k5 * ...  let sum2 be the sum of all
!     the factors of two in ntot, that is, sum2 = 2*k2.  let sumf be
!     the sum of all other factors of ntot, that is, sumf = 3*k3*5*k5*..
!     the time taken by a multidimensional transform on these ntot data
!     is t = t0 + ntot*(t1+t2*sum2+t3*sumf).  on the cdc 3300 (floating
!     point add time = six microseconds), t = 3000 + ntot*(600+40*sum2+
!     175*sumf) microseconds on complex data.
!
!     implementation of the definition by summation will run in a time
!     proportional to ntot*(nn(1)+nn(2)+...).  for highly composite ntot
!     the savings offered by this program can be dramatic.  a one-dimen-
!     sional array 4000 in length will be transformed in 4000*(600+
!     40*(2+2+2+2+2)+175*(5+5+5)) = 14.5 seconds versus about 4000*
!     4000*175 = 2800 seconds for the straightforward technique.
!
!     the fast fourier transform places three restrictions upon the
!     data.
!     1.  the number of input data and the number of transform values
!     must be the same.
!     2.  both the input data and the transform values must represent
!     equispaced points in their respective domains of time and
!     frequency.  calling these spacings deltat and deltaf, it must be
!     true that deltaf=2*pi/(nn(i)*deltat).  of course, deltat need not
!     be the same for every dimension.
!     3.  conceptually at least, the input data and the transform output
!     represent single cycles of periodic functions.
!
!     the calling sequence is--
!     call fourt(data,nn,ndim,isign,iform,work)
!
!     data is the array used to hold the real and imaginary parts
!     of the data on input and the transform values on output.  it
!     is a multidimensional floating point array, with the real and
!     imaginary parts of a datum stored immediately adjacent in storage
!     (such as fortran iv places them).  normal fortran ordering is
!     expected, the first subscript changing fastest.  the dimensions
!     are given in the integer array nn, of length ndim.  isign is -1
!     to indicate a forward transform (exponential sign is -) and +1
!     for an inverse transform (sign is +).  iform is +1 if the data are
!     complex, 0 if the data are real.  if it is 0, the imaginary
!     parts of the data must be set to zero.  as explained above, the
!     transform values are always complex and are stored in array data.
!     work is an array used for working storage.  it is floating point
!     real, one dimensional of length equal to twice the largest array
!     dimension nn(i) that is not a power of two.  if all nn(i) are
!     powers of two, it is not needed and may be replaced by zero in the
!     calling sequence.  thus, for a one-dimensional array, nn(1) odd,
!     work occupies as many storage locations as data.  if supplied,
!     work must not be the same array as data.  all subscripts of all
!     arrays begin at one.
!
!     example 1.  three-dimensional forward fourier transform of a
!     complex array dimensioned 32 by 25 by 13 in fortran iv.
!     dimension data(32,25,13),work(50),nn(3)
!     complex data
!     data nn/32,25,13/
!     do 1 i=1,32
!     do 1 j=1,25
!     do 1 k=1,13
!  1  data(i,j,k)=complex value
!     call fourt(data,nn,3,-1,1,work)
!
!     example 2.  one-dimensional forward transform of a real array of
!     length 64 in fortran ii,
!     dimension data(2,64)
!     do 2 i=1,64
!     data(1,i)=real part
!  2  data(2,i)=0.
!     call fourt(data,64,1,-1,0,0)
!
!     there are no error messages or error halts in this program.  the
!     program returns immediately if ndim or any nn(i) is less than one.
!
!     program by norman brenner from the basic program by charles
!     rader,  june 1967.  the idea for the digit reversal was
!     suggested by ralph alter.
!
!     this is the fastest and most versatile version of the fft known
!     to the author.  a program called four2 is available that also
!     performs the fast fourier transform and is written in usasi basic
!     fortran.  it is about one third as long and restricts the
!     dimensions of the input array (which must be complex) to be powers
!     of two.  another program, called four1, is one tenth as long and
!     runs two thirds as fast on a one-dimensional complex array whose
!     length is a power of two.
!
!     reference--
!     ieee audio transactions (june 1967), special issue on the fft.

      INTEGER, INTENT(IN) :: nn(3),ndim,isign,iform
      COMPLEX(DPC), INTENT(INOUT) :: data3(nn(1),nn(2),nn(3))
      REAL(EB), INTENT (INOUT) :: work(2*nn(1))

      INTEGER :: ifact(32),ntot,idim,np1,n,np2,m,ntwo,iif,idiv,iquot,irem,inon2,     &
                 icase,ifmin,i1rng,i,j,k,np2hf,i2,i1max,i1,i3,j3,nwork,ifp2,ifp1,    &
                 i2max,np1tw,ipar,k1,k2,mmax,lmax,l,kmin,kdif,kstep,k3,k4,np1hf,     &
                 j1min,j1,j2min,j2max,j2,j2rng,j1max,j3max,jmin,jmax,iconj,nhalf,    &
                 imin,imax,nprev=0,np0=0

      REAL(EB) :: data(2*nn(1)*nn(2)*nn(3)),tempr,tempi,u1r,u1i,u2r,u2i,u3r,         &
                  u3i,u4r,u4i,t2r,t2i,t3r,t3i,t4r,t4i,sumr,sumi,oldsr,oldsi,         &
                  difr,difi,theta,wr,wi,w2r,w2i,w3r,w3i,wstpr,wstpi,twowr

      REAL(EB), PARAMETER :: twopi=6.2831853071796, rthlf=0.70710678118655

      ! reshape data3 to 1D array
      data=0._EB
      n=1
      do k=1,nn(3)
        do j=1,nn(2)
          do i=1,nn(1)
            data(n)=real(data3(i,j,k),EB)
            data(n+1)=aimag(data3(i,j,k))
            n=n+2
          enddo
        enddo
      enddo

      if(ndim-1.lt.0) goto 920
      ntot=2
      do idim=1,ndim
      if(nn(idim).le.0) goto 920
      ntot=ntot*nn(idim)
      enddo
!
!     main loop for each dimension
!
      np1=2
      do idim=1,ndim
      n=nn(idim)
      np2=np1*n
      if(n-1.lt.0) then
         goto 920
      elseif(n-1.eq.0) then
         goto 900
      endif
!
!     is n a power of two and if not, what are its factors
!
      m=n
      ntwo=np1
      iif=1
      idiv=2
10    iquot=m/idiv
      irem=m-idiv*iquot
      if(iquot-idiv.lt.0) goto 50
      if(irem.ne.0) goto 20
      ntwo=ntwo+ntwo
      ifact(iif)=idiv
      iif=iif+1
      m=iquot
      go to 10
20    idiv=3
      inon2=iif
30    iquot=m/idiv
      irem=m-idiv*iquot
      if(iquot-idiv.lt.0) goto 60
      if(irem.ne.0) goto 40
      ifact(iif)=idiv
      iif=iif+1
      m=iquot
      go to 30
40    idiv=idiv+2
      go to 30
50    inon2=iif
      if(irem.ne.0) goto 60
      ntwo=ntwo+ntwo
      go to 70
60    ifact(iif)=m
!
!     separate four cases--
!        1. complex transform or real transform for the 4th, 9th,etc.
!           dimensions.
!        2. real transform for the 2nd or 3rd dimension.  method--
!           transform half the data, supplying the other half by con-
!           jugate symmetry.
!        3. real transform for the 1st dimension, n odd.  method--
!           set the imaginary parts to zero.
!        4. real transform for the 1st dimension, n even.  method--
!           transform a complex array of length n/2 whose real parts
!           are the even numbered real values and whose imaginary parts
!           are the odd numbered real values.  separate and supply
!           the second half by conjugate symmetry.
!
70    icase=1
      ifmin=1
      i1rng=np1
      if(idim-4.ge.0) goto 100
      if(iform.gt.0) goto 100
      icase=2
      i1rng=np0*(1+nprev/2)
      if(idim-1.gt.0) goto 100
      icase=3
      i1rng=np1
      if(ntwo-np1.le.0) goto 100
      icase=4
      ifmin=2
      ntwo=ntwo/2
      n=n/2
      np2=np2/2
      ntot=ntot/2
      i=1
      do j=1,ntot
      data(j)=data(i)
      i=i+2
      enddo
!
!     shuffle data by bit reversal, since n=2**k.  as the shuffling
!     can be done by simple interchange, no working array is needed
!
100   if(ntwo-np2.lt.0) goto 200
      np2hf=np2/2
      j=1
      do i2=1,np2,np1
      if(j-i2.ge.0) goto 130
      i1max=i2+np1-2
      do i1=i2,i1max,2
         do i3=i1,ntot,np2
            j3=j+i3-i2
            tempr=data(i3)
            tempi=data(i3+1)
            data(i3)=data(j3)
            data(i3+1)=data(j3+1)
            data(j3)=tempr
            data(j3+1)=tempi
         enddo
      enddo
130   m=np2hf
140   if(j-m.le.0) goto 150
      j=j-m
      m=m/2
      if(m-np1.ge.0) goto 140
150   j=j+m
      enddo
      goto 300
!
!     shuffle data by digit reversal for general n
!
200   nwork=2*n
      do i1=1,np1,2
      do i3=i1,ntot,np2
      j=i3
      do 260 i=1,nwork,2
      if(icase-3.ne.0) goto 220
      work(i)=data(j)
      work(i+1)=data(j+1)
      go to 230
220   work(i)=data(j)
      work(i+1)=0._EB
230   ifp2=np2
      iif=ifmin
240   ifp1=ifp2/ifact(iif)
      j=j+ifp1
      if(j-i3-ifp2.lt.0) goto 260
      j=j-ifp2
      ifp2=ifp1
      iif=iif+1
      if(ifp2-np1.gt.0) goto 240
260   continue
      i2max=i3+np2-np1
      i=1
      do i2=i3,i2max,np1
      data(i2)=work(i)
      data(i2+1)=work(i+1)
      i=i+2
      enddo
      enddo
      enddo
!
!     main loop for factors of two.  perform fourier transforms of
!     length four, with one of length two if needed.  the twiddle factor
!     w=exp(isign*2*pi*sqrt(-1)*m/(4*mmax)).  check for w=isign*sqrt(-1)
!     and repeat for w=w*(1+isign*sqrt(-1))/sqrt(2).
!
300   if(ntwo-np1.le.0) goto 600
      np1tw=np1+np1
      ipar=ntwo/np1
310   if(ipar-2.lt.0) then
         goto 350
      elseif(ipar-2.eq.0) then
         goto 330
      endif
      ipar=ipar/4
      go to 310
330   do i1=1,i1rng,2
      do k1=i1,ntot,np1tw
      k2=k1+np1
      tempr=data(k2)
      tempi=data(k2+1)
      data(k2)=data(k1)-tempr
      data(k2+1)=data(k1+1)-tempi
      data(k1)=data(k1)+tempr
      data(k1+1)=data(k1+1)+tempi
      enddo
      enddo
350   mmax=np1
360   if(mmax-ntwo/2.ge.0) goto 600
      lmax=max0(np1tw,mmax/2)
      do 570 l=np1,lmax,np1tw
      m=l
      if(mmax-np1.le.0) goto 420
      theta=-twopi*REAL(l,EB)/REAL(4*mmax,EB)
      if(isign.lt.0) goto 400
      theta=-theta
400   wr=cos(theta)
      wi=sin(theta)
410   w2r=wr*wr-wi*wi
      w2i=2._EB*wr*wi
      w3r=w2r*wr-w2i*wi
      w3i=w2r*wi+w2i*wr
420   do 530 i1=1,i1rng,2
      kmin=i1+ipar*m
      if(mmax-np1.gt.0) goto 440
      kmin=i1
440   kdif=ipar*mmax
450   kstep=4*kdif
      if(kstep-ntwo.gt.0) goto 530
      do k1=kmin,ntot,kstep
      k2=k1+kdif
      k3=k2+kdif
      k4=k3+kdif
      if(mmax-np1.gt.0) goto 480
      u1r=data(k1)+data(k2)
      u1i=data(k1+1)+data(k2+1)
      u2r=data(k3)+data(k4)
      u2i=data(k3+1)+data(k4+1)
      u3r=data(k1)-data(k2)
      u3i=data(k1+1)-data(k2+1)
      if(isign.ge.0) goto 472
      u4r=data(k3+1)-data(k4+1)
      u4i=data(k4)-data(k3)
      go to 510
472   u4r=data(k4+1)-data(k3+1)
      u4i=data(k3)-data(k4)
      go to 510
480   t2r=w2r*data(k2)-w2i*data(k2+1)
      t2i=w2r*data(k2+1)+w2i*data(k2)
      t3r=wr*data(k3)-wi*data(k3+1)
      t3i=wr*data(k3+1)+wi*data(k3)
      t4r=w3r*data(k4)-w3i*data(k4+1)
      t4i=w3r*data(k4+1)+w3i*data(k4)
      u1r=data(k1)+t2r
      u1i=data(k1+1)+t2i
      u2r=t3r+t4r
      u2i=t3i+t4i
      u3r=data(k1)-t2r
      u3i=data(k1+1)-t2i
      if(isign.ge.0) goto 500
      u4r=t3i-t4i
      u4i=t4r-t3r
      go to 510
500   u4r=t4i-t3i
      u4i=t3r-t4r
510   data(k1)=u1r+u2r
      data(k1+1)=u1i+u2i
      data(k2)=u3r+u4r
      data(k2+1)=u3i+u4i
      data(k3)=u1r-u2r
      data(k3+1)=u1i-u2i
      data(k4)=u3r-u4r
      data(k4+1)=u3i-u4i
      enddo
      kdif=kstep
      kmin=4*(kmin-i1)+i1
      go to 450
530   continue
      m=m+lmax
      if(m-mmax.gt.0) goto 570
      if(isign.ge.0) goto 560
      tempr=wr
      wr=(wr+wi)*rthlf
      wi=(wi-tempr)*rthlf
      go to 410
560   tempr=wr
      wr=(wr-wi)*rthlf
      wi=(tempr+wi)*rthlf
      go to 410
570   continue
      ipar=3-ipar
      mmax=mmax+mmax
      go to 360
!
!     main loop for factors not equal to two.  apply the twiddle factor
!     w=exp(isign*2*pi*sqrt(-1)*(j1-1)*(j2-j1)/(ifp1+ifp2)), then
!     perform a fourier transform of length ifact(iif), making use of
!     conjugate symmetries.
!
600   if(ntwo-np2.ge.0) goto 700
      ifp1=ntwo
      iif=inon2
      np1hf=np1/2
610   ifp2=ifact(iif)*ifp1
      j1min=np1+1
      if(j1min-ifp1.gt.0) goto 640
      do j1=j1min,ifp1,np1
      theta=-twopi*REAL(j1-1,EB)/REAL(ifp2,EB)
      if(isign.lt.0) goto 625
      theta=-theta
625   wstpr=cos(theta)
      wstpi=sin(theta)
      wr=wstpr
      wi=wstpi
      j2min=j1+ifp1
      j2max=j1+ifp2-ifp1
      do j2=j2min,j2max,ifp1
      i1max=j2+i1rng-2
      do i1=j2,i1max,2
      do j3=i1,ntot,ifp2
      tempr=data(j3)
      data(j3)=data(j3)*wr-data(j3+1)*wi
      data(j3+1)=tempr*wi+data(j3+1)*wr
      enddo
      enddo
      tempr=wr
      wr=wr*wstpr-wi*wstpi
      wi=tempr*wstpi+wi*wstpr
      enddo
      enddo
640   theta=-twopi/REAL(ifact(iif),EB)
      if(isign.lt.0) goto 650
      theta=-theta
650   wstpr=cos(theta)
      wstpi=sin(theta)
      j2rng=ifp1*(1+ifact(iif)/2)
      do i1=1,i1rng,2
      do i3=i1,ntot,np2
      j2max=i3+j2rng-ifp1
      do j2=i3,j2max,ifp1
      j1max=j2+ifp1-np1
      do j1=j2,j1max,np1
      j3max=j1+np2-ifp2
      do j3=j1,j3max,ifp2
      jmin=j3-j2+i3
      jmax=jmin+ifp2-ifp1
      i=1+(j3-i3)/np1hf
      if(j2-i3.gt.0) goto 665
      sumr=0._EB
      sumi=0._EB
      do j=jmin,jmax,ifp1
      sumr=sumr+data(j)
      sumi=sumi+data(j+1)
      enddo
      work(i)=sumr
      work(i+1)=sumi
      go to 680
665   iconj=1+(ifp2-2*j2+i3+j3)/np1hf
      j=jmax
      sumr=data(j)
      sumi=data(j+1)
      oldsr=0._EB
      oldsi=0._EB
      j=j-ifp1
670   tempr=sumr
      tempi=sumi
      sumr=twowr*sumr-oldsr+data(j)
      sumi=twowr*sumi-oldsi+data(j+1)
      oldsr=tempr
      oldsi=tempi
      j=j-ifp1
      if(j-jmin.gt.0) goto 670
      tempr=wr*sumr-oldsr+data(j)
      tempi=wi*sumi
      work(i)=tempr-tempi
      work(iconj)=tempr+tempi
      tempr=wr*sumi-oldsi+data(j+1)
      tempi=wi*sumr
      work(i+1)=tempr+tempi
      work(iconj+1)=tempr-tempi
680   continue
      enddo
      enddo
      if(j2-i3.gt.0) goto 686
      wr=wstpr
      wi=wstpi
      go to 690
686   tempr=wr
      wr=wr*wstpr-wi*wstpi
      wi=tempr*wstpi+wi*wstpr
690   twowr=wr+wr
      enddo
      i=1
      i2max=i3+np2-np1
      do i2=i3,i2max,np1
      data(i2)=work(i)
      data(i2+1)=work(i+1)
      i=i+2
      enddo
      enddo
      enddo
      iif=iif+1
      ifp1=ifp2
      if(ifp1-np2.lt.0) goto 610
!
!     complete a real transform in the 1st dimension, n even, by con-
!     jugate symmetries.
!
700   select case(icase)
         case(1,3)
            goto 900
         case(2)
            goto 800
         case(4)
            goto 701
      end select
701   nhalf=n
      n=n+n
      theta=-twopi/REAL(n,EB)
      if(isign.lt.0) goto 703
      theta=-theta
703   wstpr=cos(theta)
      wstpi=sin(theta)
      wr=wstpr
      wi=wstpi
      imin=3
      jmin=2*nhalf-1
      go to 725
710   j=jmin
      do i=imin,ntot,np2
      sumr=(data(i)+data(j))/2._EB
      sumi=(data(i+1)+data(j+1))/2._EB
      difr=(data(i)-data(j))/2._EB
      difi=(data(i+1)-data(j+1))/2._EB
      tempr=wr*sumi+wi*difr
      tempi=wi*sumi-wr*difr
      data(i)=sumr+tempr
      data(i+1)=difi+tempi
      data(j)=sumr-tempr
      data(j+1)=-difi+tempi
      j=j+np2
      enddo
      imin=imin+2
      jmin=jmin-2
      tempr=wr
      wr=wr*wstpr-wi*wstpi
      wi=tempr*wstpi+wi*wstpr
725   if(imin-jmin.lt.0) then
         goto 710
      elseif(imin-jmin.gt.0) then
         goto 740
      endif
      if(isign.ge.0) goto 740
      do i=imin,ntot,np2
      data(i+1)=-data(i+1)
      enddo
740   np2=np2+np2
      ntot=ntot+ntot
      j=ntot+1
      imax=ntot/2+1
745   imin=imax-2*nhalf
      i=imin
      go to 755
750   data(j)=data(i)
      data(j+1)=-data(i+1)
755   i=i+2
      j=j-2
      if(i-imax.lt.0) goto 750
      data(j)=data(imin)-data(imin+1)
      data(j+1)=0._EB
      if(i-j.lt.0) then
         goto 770
      else
         goto 780
      endif
765   data(j)=data(i)
      data(j+1)=data(i+1)
770   i=i-2
      j=j-2
      if(i-imin.gt.0) goto 765
      data(j)=data(imin)+data(imin+1)
      data(j+1)=0._EB
      imax=imin
      go to 745
780   data(1)=data(1)+data(2)
      data(2)=0._EB
      go to 900
!
!     complete a real transform for the 2nd or 3rd dimension by
!     conjugate symmetries.
!
800   if(i1rng-np1.ge.0) goto 900
      do i3=1,ntot,np2
      i2max=i3+np2-np1
      do i2=i3,i2max,np1
      imin=i2+i1rng
      imax=i2+np1-2
      jmax=2*i3+np1-imin
      if(i2-i3.le.0) goto 820
      jmax=jmax+np2
820   if(idim-2.le.0) goto 850
      j=jmax+np0
      do i=imin,imax,2
      data(i)=data(j)
      data(i+1)=-data(j+1)
      j=j-2
      enddo
850   j=jmax
      do i=imin,imax,np0
      data(i)=data(j)
      data(i+1)=-data(j+1)
      j=j-np0
      enddo
      enddo
      enddo
!
!     end of loop on each dimension
!
900   np0=np1
      np1=np2
      nprev=n
      enddo

      ! reshape data back to 3D complex array

      !! for debug purposes (move to 920)
      !print *,size(data)
      !do i=1,size(data)
      !   print *,data(i)
      !enddo
      !stop

920   n=1
      do k=1,nn(3)
        do j=1,nn(2)
          do i=1,nn(1)
            data3(i,j,k)=cmplx(data(n),data(n+1),kind=DPC)
            n=n+2
          enddo
        enddo
      enddo
      return

END SUBROUTINE fourier_f2003

END MODULE TURBULENCE



MODULE MANUFACTURED_SOLUTIONS

USE PRECISION_PARAMETERS

IMPLICIT NONE

! L. Shunn, F. Ham, P. Moin. Verification of variable-density flow solvers using manufactured solutions.
!    Journal of Computational Physics 231 (2012) 3801-3827. (below is the code for Shunn Problem 3)

REAL(EB), PRIVATE :: K,W,D,MU,R0,R1,UF,VF
REAL(EB), PUBLIC :: RHO_0_MMS,RHO_1_MMS
REAL(EB), PUBLIC, PARAMETER :: UF_MMS=0.5_EB,WF_MMS=0.5_EB,VISC_MMS=0.001_EB,DIFF_MMS=0.001_EB,WAVE_NUM_MMS=2._EB,FREQ_MMS=2._EB

PRIVATE
PUBLIC :: SHUNN_MMS_3,VD2D_MMS_Z,VD2D_MMS_RHO,VD2D_MMS_RHO_OF_Z,VD2D_MMS_U,VD2D_MMS_V,VD2D_MMS_P, &
          VD2D_MMS_Z_SRC,VD2D_MMS_U_SRC,VD2D_MMS_V_SRC,VD2D_MMS_Z_OF_RHO,VD2D_MMS_DIV,SAAD_MMS_1,VD2D_MMS_H, &
          VD2D_MMS_U_SRC_3,VD2D_MMS_V_SRC_3,VD2D_MMS_H_3,VD2D_MMS_P_3

CONTAINS

SUBROUTINE SAAD_MMS_1(NM)
USE MESH_POINTERS
USE MESH_VARIABLES
USE GLOBAL_CONSTANTS, ONLY: P_INF,TMPA,N_TRACKED_SPECIES
USE PHYSICAL_FUNCTIONS, ONLY: GET_SPECIFIC_GAS_CONSTANT
IMPLICIT NONE
INTEGER, INTENT(IN) :: NM
INTEGER :: I,J,K
REAL(EB) :: ZZ_GET(1:N_TRACKED_SPECIES)
REAL(EB), PARAMETER :: L=2._EB, X0=0._EB, SIGMA=0.1_EB
REAL(EB), PARAMETER :: GAS_CONSTANT=8314.472_EB

! set up problem parameters

RHO_0_MMS = P_INF*SPECIES_MIXTURE(1)%MW/(GAS_CONSTANT*TMPA)
RHO_1_MMS = P_INF*SPECIES_MIXTURE(2)%MW/(GAS_CONSTANT*TMPA)

CALL POINT_TO_MESH(NM)

DO K=0,KBP1
   DO J=0,JBP1
      DO I=0,IBP1
         !ZZ(I,J,K,2) = EXP(-.5_EB*(XC(I)-X0)**2/SIGMA**2)
         ZZ(I,J,K,2) = 0.5_EB*(1._EB+SIN(2._EB*PI*XC(I)/L))
         ZZ(I,J,K,1) = 1._EB - ZZ(I,J,K,2)
         ZZ_GET(1:2) = ZZ(I,J,K,1:2)
         CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,RSUM(I,J,K))
         RHO(I,J,K) = P_INF/(TMPA*RSUM(I,J,K))
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE SAAD_MMS_1

SUBROUTINE SHUNN_MMS_3(DT,NM)
USE MESH_POINTERS
USE MESH_VARIABLES
USE GLOBAL_CONSTANTS, ONLY: P_INF,TMPA,N_TRACKED_SPECIES
USE PHYSICAL_FUNCTIONS, ONLY: GET_SPECIFIC_GAS_CONSTANT
IMPLICIT NONE
INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: DT
INTEGER :: I,J,K
REAL(EB) :: ZZ_GET(1:N_TRACKED_SPECIES)
REAL(EB), PARAMETER :: GAS_CONSTANT=8314.472_EB

! set up problem parameters

RHO_0_MMS = P_INF*SPECIES_MIXTURE(1)%MW/(GAS_CONSTANT*TMPA)
RHO_1_MMS = P_INF*SPECIES_MIXTURE(2)%MW/(GAS_CONSTANT*TMPA)

! initialize private variables

CALL VD2D_MMS_INIT

! initialize flow fields

CALL POINT_TO_MESH(NM)
U = UF_MMS; US = UF_MMS
V = 0._EB
W = WF_MMS; WS = WF_MMS
DO K=0,KBP1
   DO J=0,JBP1
      DO I=0,IBP1
         ZZ(I,J,K,2) = VD2D_MMS_Z(XC(I),ZC(K),0._EB)
         ZZ(I,J,K,1) = 1._EB - ZZ(I,J,K,2)
         ZZ_GET(1:2) = ZZ(I,J,K,1:2)
         CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,RSUM(I,J,K))
         RHO(I,J,K) = P_INF/(TMPA*RSUM(I,J,K)) !; print *, RHO(I,J,K), VD2D_MMS_RHO(XC(I),ZC(K),0._EB)
         H(I,J,K)  = VD2D_MMS_H_3(XC(I),ZC(K),0._EB)
         HS(I,J,K) = VD2D_MMS_H_3(XC(I),ZC(K),DT)
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE SHUNN_MMS_3

SUBROUTINE VD2D_MMS_INIT
IMPLICIT NONE
K = WAVE_NUM_MMS
W = FREQ_MMS
R0 = RHO_0_MMS
R1 = RHO_1_MMS
D = DIFF_MMS
MU = VISC_MMS
UF = UF_MMS
VF = WF_MMS
END SUBROUTINE VD2D_MMS_INIT

FUNCTION VD2D_MMS_RHO_OF_Z(Z)
IMPLICIT NONE
REAL(EB),INTENT (IN) :: Z
REAL(EB) :: VD2D_MMS_RHO_OF_Z
VD2D_MMS_RHO_OF_Z = 1.0_EB/(Z/R1+(1.0_EB-Z)/R0)
END FUNCTION VD2D_MMS_RHO_OF_Z

REAL(EB) FUNCTION VD2D_MMS_Z_OF_RHO(RHO)
IMPLICIT NONE
REAL(EB),INTENT (IN) :: RHO
VD2D_MMS_Z_OF_RHO = (R0-RHO)/(R0-R1)*R1/RHO
END FUNCTION VD2D_MMS_Z_OF_RHO

FUNCTION VD2D_MMS_Z(X,Y,T)
IMPLICIT NONE
REAL(EB),INTENT (IN) :: X,Y,T
REAL(EB) :: VD2D_MMS_Z
VD2D_MMS_Z = &
   (1.0_EB + SIN(K * PI * X) * SIN(K * PI * Y) *COS (W * PI * T))/ &
   (1.0_EB + R0/R1+(1.0_EB-R0/R1) &
   * SIN(K * PI * X) * SIN(K * PI * Y) * COS(W * PI * T))
END FUNCTION VD2D_MMS_Z

FUNCTION VD2D_MMS_RHO(X,Y,T)
IMPLICIT NONE
REAL(EB),INTENT (IN) :: X,Y,T
REAL(EB) :: VD2D_MMS_RHO
VD2D_MMS_RHO = VD2D_MMS_RHO_OF_Z(VD2D_MMS_Z(X,Y,T))
END FUNCTION VD2D_MMS_RHO

FUNCTION VD2D_MMS_U(X,Y,T)
IMPLICIT NONE
REAL(EB),INTENT (IN) :: X,Y,T
REAL(EB) :: VD2D_MMS_U
! VD2D_MMS_U=-W/K/4.0_EB * COS(K * PI * X) * SIN(K * PI * Y) * SIN(W * PI * T) &
!    * (R1-R0)/VD2D_MMS_RHO(X,Y,T)
VD2D_MMS_U = UF_MMS + (R1-R0)/VD2D_MMS_RHO(X,Y,T)*(-W/(4._EB*K))*COS(PI*K*X)*SIN(PI*K*Y)*SIN(PI*W*T)
END FUNCTION VD2D_MMS_U

FUNCTION VD2D_MMS_V(X,Y,T)
IMPLICIT NONE
REAL(EB),INTENT (IN) :: X,Y,T
REAL(EB) :: VD2D_MMS_V
! VD2D_MMS_V = VD2D_MMS_U(Y,X,T)
! !! NOTE: V OBTAINED FROM U BY TRANSPOSING X <--> Y
VD2D_MMS_V = WF_MMS + (R1-R0)/VD2D_MMS_RHO(X,Y,T)*(-W/(4._EB*K))*SIN(PI*K*X)*COS(PI*K*Y)*SIN(PI*W*T)
END FUNCTION VD2D_MMS_V

FUNCTION VD2D_MMS_P(X,Y,T)
IMPLICIT NONE
REAL(EB),INTENT (IN) :: X,Y,T
REAL(EB) :: VD2D_MMS_P
VD2D_MMS_P = &
   VD2D_MMS_RHO(X,Y,T) * VD2D_MMS_U(X,Y,T) * VD2D_MMS_V(X,Y,T) / 2.0_EB
   !! NOTE: PRESSURE CAN VARY BY A CONSTANT AND STILL SATISFY THE MMS
END FUNCTION VD2D_MMS_P

FUNCTION VD2D_MMS_P_3(X,Y,T)
IMPLICIT NONE
REAL(EB),INTENT (IN) :: X,Y,T
REAL(EB) :: VD2D_MMS_P_3
VD2D_MMS_P_3 = &
   VD2D_MMS_RHO(X,Y,T) * ( VD2D_MMS_H_3(X,Y,T) - 0.5_EB*(VD2D_MMS_U(X,Y,T)**2 + VD2D_MMS_V(X,Y,T)**2) )
END FUNCTION VD2D_MMS_P_3

FUNCTION VD2D_MMS_H(X,Y,T)
IMPLICIT NONE
REAL(EB),INTENT (IN) :: X,Y,T
REAL(EB) :: VD2D_MMS_H
VD2D_MMS_H = &
   VD2D_MMS_P(X,Y,T)/VD2D_MMS_RHO(X,Y,T) + 0.5_EB*(VD2D_MMS_U(X,Y,T)**2 + VD2D_MMS_V(X,Y,T)**2)
END FUNCTION VD2D_MMS_H

FUNCTION VD2D_MMS_H_3(X,Y,T)
IMPLICIT NONE
REAL(EB),INTENT (IN) :: X,Y,T
REAL(EB) :: VD2D_MMS_H_3
VD2D_MMS_H_3 = 0.5_EB * (VD2D_MMS_U(X,Y,T)-UF_MMS) * (VD2D_MMS_V(X,Y,T)-WF_MMS)
END FUNCTION VD2D_MMS_H_3

FUNCTION VD2D_MMS_Z_SRC(X,Y,T)
IMPLICIT NONE
REAL(EB),INTENT (IN) :: X,Y,T
REAL(EB) :: VD2D_MMS_Z_SRC
REAL(EB) :: S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S16
S0 = COS(W * PI * T)
S1 = SIN(W * PI * T)
S2 = COS(K * PI * X)
S3 = SIN(K * PI * X)
S4 = COS(K * PI * Y)
S5 = SIN(K * PI * Y)
S6 = R1-R0
S7 = 1.0_EB + S3 * S5 * S0
S8 = (0.5_EB + 0.5_EB * S3 * S5 * S0) * S6 + R0 ! RHO
S9 = R1 * S3 * S5 * S0 + R1-S3 * S5 * S0 * R0 + R0
S10 = R1 * S2 * K * PI * S5 * S0-S2 * K * PI * S5 * S0 * R0
S11 = R1 * S3 * S4 * K * PI * S0-S3 * S4 * K * PI * S0 * R0
S12 = R0 * S3 * S5 * S1 * W * PI-S3 * S5 * S1 * W * PI * R1
S13 = S3 * K ** 2 * PI ** 2 * S5 * S0 * R0-R1 * S3 * K ** 2 * PI ** 2 * S5 * S0
S14 = 1.0_EB/4.0_EB
S15 = 2.0_EB * D
S16 = S15 * PI * PI
VD2D_MMS_Z_SRC=&
   -R1 * (S8 * S3 * S5 * S1 * W * PI * S9 ** 2 * K + S8 * S7 * S12 * S9 * K + S14 * S2 ** 2 * PI * S5 ** 2 &
   * S0 * S1 * W * S6 * S9 ** 2 * K-S14 * S7 * S5 * S1 * W * S6 * S2 * S10 * S9 + S14 * S3 ** 2 * S4 ** 2 * PI &
   * S0 * S1 * W * S6 * S9 ** 2 * K-S14 * S7 * S3 * S1 * W * S6 * S4 * S11 * S9-S16 * K ** 3 * S3 * S5 * S0 &
   * S9 ** 2-S15 * K ** 2 * S2 * PI * S5 * S0 * S10 * S9 + S15 * K * S7 * S10 ** 2-S15 * K * S7 * S13 * S9 &
   -S15 * K ** 2 * S3 * S4 * PI * S0 * S11 * S9 + S15 * K * S7 * S11 ** 2)/S9 ** 3/K
END FUNCTION VD2D_MMS_Z_SRC

FUNCTION VD2D_MMS_U_SRC(X,Y,T)
IMPLICIT NONE
REAL(EB),INTENT (IN) :: X,Y,T
REAL(EB) :: VD2D_MMS_U_SRC
REAL(EB) :: S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15
S0 = COS(W * PI * T)
S1 = SIN(W * PI * T)
S2 = COS(K * PI * X)
S3 = SIN(K * PI * X)
S4 = COS(K * PI * Y)
S5 = SIN(K * PI * Y)
S6 = R1-R0
S7=(0.5_EB + 0.5_EB * S3 * S5 * S0) * S6 + R0 ! RHO
S8 = 1.0_EB/64.0_EB
S9 = S8 * 2.0_EB
S10 = S9 * 2.0_EB
S11 = S10 * 3.0_EB
S12 = S10 * 4.0_EB
S13 = PI * PI/6.0_EB
S14 = S13 * 2.0_EB
S15 = S14 * 2.0_EB
VD2D_MMS_U_SRC=&
   W * S6 * (-S12 * S5 * S0 * W * PI * S2 * S7 ** 3-S11 * S5 ** 2 * S1 ** 2 * W * S6 * S2 * S3 * PI * S7 ** 2 &
   -S9 * S5 ** 3 * S1 ** 2 * W * S6 ** 2 * S2 ** 3 * PI * S0 * S7 + S10 * S4 ** 2 * PI * S1 ** 2 * W * S6 * S2 * S3 &
   * S7 ** 2-S9 * S5 * S1 ** 2 * W * S6 ** 2 * S2 * S3 ** 2 * S4 ** 2 * PI * S0 * S7-S15 * MU * S1 * S2 * K ** 2 &
   * S5 * S7 ** 2 + S15 * MU * S1 * S2 * K ** 2 * S5 ** 2 * S6 * S3 * S0 * S7 + S13 * MU * S1 * S2 ** 3 * K ** 2 &
   * S5 ** 3 * S6 ** 2 * S0 ** 2-S14 * MU * S1 * S2 * K ** 2 * S3 * S6 * S4 ** 2 * S0 * S7 + S13 * MU * S1 * S2 &
   * K ** 2 * S3 ** 2 * S6 ** 2 * S4 ** 2 * S0 ** 2 * S5-S9 * S5 * S1 ** 2 * W * S6 * S3 ** 2 * PI * S4 * S7 ** 2 &
   + S9 * S5 * S1 ** 2 * W * S6 * S2 ** 2 * PI * S4 * S7 ** 2-S8 * S5 ** 2 * S1 ** 2 * W * S6 ** 2 * S2 ** 2 * S3 &
   * S4 * PI * S0 * S7)/K/S7 ** 3
END FUNCTION VD2D_MMS_U_SRC

FUNCTION VD2D_MMS_V_SRC(X,Y,T)
IMPLICIT NONE
REAL(EB),INTENT (IN) :: X,Y,T
REAL(EB) :: VD2D_MMS_V_SRC
VD2D_MMS_V_SRC = VD2D_MMS_U_SRC(Y,X,T)
!! NOTE: V OBTAINED FROM U BY TRANSPOSING X <--> Y
END FUNCTION VD2D_MMS_V_SRC

REAL(EB) FUNCTION VD2D_MMS_DIV(X,Y,T)
IMPLICIT NONE
REAL(EB),INTENT (IN) :: X,Y,T ! X=>XHAT, Y=>YHAT
REAL(EB) :: S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,DUDX,DVDY,DZDX,DZDY
S0 = COS(W * PI * T)
S1 = SIN(W * PI * T)
S2 = COS(K * PI * X)
S3 = SIN(K * PI * X)
S4 = COS(K * PI * Y)
S5 = SIN(K * PI * Y)
S6 = R1-R0
S7 = VD2D_MMS_RHO(X,Y,T)
S8 = -W/(4._EB*K)

S9 = (1._EB + R0/R1)
S10 = S3 * S5 * S0
S11 = 1._EB + S10
S12 = S9 + S9 * S10

DZDX = S11 * (-S9 * S2 * S5 * S0 * PI * K) / S12**2 + (S2 * S5 * S0 * PI * K) / S12
DZDY = S11 * (-S9 * S3 * S4 * S0 * PI * K) / S12**2 + (S3 * S4 * S0 * PI * K) / S12

S13 = (1._EB/R1 - 1._EB/R0)*DZDX
S14 = (1._EB/R1 - 1._EB/R0)*DZDY

DUDX = S6 * S8 * S5 * S1 * (S2 * S13 - PI * K * S3 / S7)
DVDY = S6 * S8 * S3 * S1 * (S4 * S14 - PI * K * S5 / S7)

VD2D_MMS_DIV = DUDX + DVDY

END FUNCTION VD2D_MMS_DIV

FUNCTION VD2D_MMS_U_SRC_3(X,Y,T)
! This source term is obtained by running the symbolic python script
! /Utilities/Python/shunn3_stokes_mms_sym.py
IMPLICIT NONE
REAL(EB),INTENT (IN) :: X,Y,T
REAL(EB) :: VD2D_MMS_U_SRC_3, UF,VF
UF=UF_MMS
VF=WF_MMS

VD2D_MMS_U_SRC_3 = -pi*uf*w*(-r0 + r1)*((sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + &
1)/(r1*(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)) + &
(-(sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)/(r0/r1 + (-r0/r1 + &
1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1) + 1)/r0)*sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*sin(pi*t*w)/4 + pi*vf*w*(-r0 + r1)*((sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf &
+ y))*cos(pi*t*w) + 1)/(r1*(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + &
y))*cos(pi*t*w) + 1)) + (-(sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)/(r0/r1 + &
(-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1) + &
1)/r0)*sin(pi*t*w)*cos(pi*k*(-t*uf + x))*cos(pi*k*(-t*vf + y))/4 - (vf - w*(-r0 + &
r1)*((sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)/(r1*(r0/r1 + (-r0/r1 + &
1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)) + (-(sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)/(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1) + 1)/r0)*sin(pi*k*(-t*uf + &
x))*sin(pi*t*w)*cos(pi*k*(-t*vf + y))/(4*k))*(w*(-r0 + r1)*(-pi*k*(-r0/r1 + 1)*(sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)*sin(pi*k*(-t*uf + x))*cos(pi*k*(-t*vf + &
y))*cos(pi*t*w)/(r1*(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + &
1)**2) + pi*k*sin(pi*k*(-t*uf + x))*cos(pi*k*(-t*vf + y))*cos(pi*t*w)/(r1*(r0/r1 + (-r0/r1 + &
1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)) + (pi*k*(-r0/r1 + &
1)*(sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)*sin(pi*k*(-t*uf + &
x))*cos(pi*k*(-t*vf + y))*cos(pi*t*w)/(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + &
y))*cos(pi*t*w) + 1)**2 - pi*k*sin(pi*k*(-t*uf + x))*cos(pi*k*(-t*vf + y))*cos(pi*t*w)/(r0/r1 + &
(-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1))/r0)*sin(pi*k*(-t*vf + &
y))*sin(pi*t*w)*cos(pi*k*(-t*uf + x))/(4*k) - w*(-r0 + r1)*(-pi*k*(-r0/r1 + 1)*(sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)*sin(pi*k*(-t*vf + y))*cos(pi*k*(-t*uf + &
x))*cos(pi*t*w)/(r1*(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + &
1)**2) + pi*k*sin(pi*k*(-t*vf + y))*cos(pi*k*(-t*uf + x))*cos(pi*t*w)/(r1*(r0/r1 + (-r0/r1 + &
1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)) + (pi*k*(-r0/r1 + &
1)*(sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)*sin(pi*k*(-t*vf + &
y))*cos(pi*k*(-t*uf + x))*cos(pi*t*w)/(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + &
y))*cos(pi*t*w) + 1)**2 - pi*k*sin(pi*k*(-t*vf + y))*cos(pi*k*(-t*uf + x))*cos(pi*t*w)/(r0/r1 + &
(-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1))/r0)*sin(pi*k*(-t*uf + &
x))*sin(pi*t*w)*cos(pi*k*(-t*vf + y))/(4*k)) - (2*mu*(pi**2*k*w*(-r0 + r1)*((sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)/(r1*(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)) + (-(sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + &
y))*cos(pi*t*w) + 1)/(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + &
1) + 1)/r0)*sin(pi*k*(-t*vf + y))*sin(pi*t*w)*cos(pi*k*(-t*uf + x))/12 + pi*w*(-r0 + &
r1)*(-pi*k*(-r0/r1 + 1)*(sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + &
1)*sin(pi*k*(-t*uf + x))*cos(pi*k*(-t*vf + y))*cos(pi*t*w)/(r1*(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf &
+ x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)**2) + pi*k*sin(pi*k*(-t*uf + x))*cos(pi*k*(-t*vf + &
y))*cos(pi*t*w)/(r1*(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + &
1)) + (pi*k*(-r0/r1 + 1)*(sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + &
1)*sin(pi*k*(-t*uf + x))*cos(pi*k*(-t*vf + y))*cos(pi*t*w)/(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)**2 - pi*k*sin(pi*k*(-t*uf + x))*cos(pi*k*(-t*vf + &
y))*cos(pi*t*w)/(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + &
1))/r0)*sin(pi*t*w)*cos(pi*k*(-t*uf + x))*cos(pi*k*(-t*vf + y))/12 + pi*w*(-r0 + r1)*(-pi*k*(-r0/r1 &
+ 1)*(sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)*sin(pi*k*(-t*vf + &
y))*cos(pi*k*(-t*uf + x))*cos(pi*t*w)/(r1*(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)**2) + pi*k*sin(pi*k*(-t*vf + y))*cos(pi*k*(-t*uf + &
x))*cos(pi*t*w)/(r1*(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + &
1)) + (pi*k*(-r0/r1 + 1)*(sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + &
1)*sin(pi*k*(-t*vf + y))*cos(pi*k*(-t*uf + x))*cos(pi*t*w)/(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)**2 - pi*k*sin(pi*k*(-t*vf + y))*cos(pi*k*(-t*uf + &
x))*cos(pi*t*w)/(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + &
1))/r0)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*sin(pi*t*w)/4 - w*(-r0 + &
r1)*(2*pi**2*k**2*(-r0/r1 + 1)**2*(sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + &
1)*sin(pi*k*(-t*vf + y))**2*cos(pi*k*(-t*uf + x))**2*cos(pi*t*w)**2/(r1*(r0/r1 + (-r0/r1 + &
1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)**3) + pi**2*k**2*(-r0/r1 + &
1)*(sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)*sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w)/(r1*(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)**2) - 2*pi**2*k**2*(-r0/r1 + 1)*sin(pi*k*(-t*vf + &
y))**2*cos(pi*k*(-t*uf + x))**2*cos(pi*t*w)**2/(r1*(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)**2) - pi**2*k**2*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + &
y))*cos(pi*t*w)/(r1*(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + &
1)) + (-2*pi**2*k**2*(-r0/r1 + 1)**2*(sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + &
1)*sin(pi*k*(-t*vf + y))**2*cos(pi*k*(-t*uf + x))**2*cos(pi*t*w)**2/(r0/r1 + (-r0/r1 + &
1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)**3 - pi**2*k**2*(-r0/r1 + &
1)*(sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)*sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w)/(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + &
y))*cos(pi*t*w) + 1)**2 + 2*pi**2*k**2*(-r0/r1 + 1)*sin(pi*k*(-t*vf + y))**2*cos(pi*k*(-t*uf + &
x))**2*cos(pi*t*w)**2/(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) &
+ 1)**2 + pi**2*k**2*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w)/(r0/r1 + (-r0/r1 + &
1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1))/r0)*sin(pi*k*(-t*vf + &
y))*sin(pi*t*w)*cos(pi*k*(-t*uf + x))/(6*k) + w*(-r0 + r1)*(2*pi**2*k**2*(-r0/r1 + &
1)**2*(sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)*sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*k*(-t*uf + x))*cos(pi*k*(-t*vf + y))*cos(pi*t*w)**2/(r1*(r0/r1 + &
(-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)**3) - pi**2*k**2*(-r0/r1 + &
1)*(sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)*cos(pi*k*(-t*uf + &
x))*cos(pi*k*(-t*vf + y))*cos(pi*t*w)/(r1*(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)**2) - 2*pi**2*k**2*(-r0/r1 + 1)*sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*k*(-t*uf + x))*cos(pi*k*(-t*vf + y))*cos(pi*t*w)**2/(r1*(r0/r1 + &
(-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)**2) + &
pi**2*k**2*cos(pi*k*(-t*uf + x))*cos(pi*k*(-t*vf + y))*cos(pi*t*w)/(r1*(r0/r1 + (-r0/r1 + &
1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)) + (-2*pi**2*k**2*(-r0/r1 + &
1)**2*(sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)*sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*k*(-t*uf + x))*cos(pi*k*(-t*vf + y))*cos(pi*t*w)**2/(r0/r1 + &
(-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)**3 + pi**2*k**2*(-r0/r1 + &
1)*(sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)*cos(pi*k*(-t*uf + &
x))*cos(pi*k*(-t*vf + y))*cos(pi*t*w)/(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + &
y))*cos(pi*t*w) + 1)**2 + 2*pi**2*k**2*(-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + &
y))*cos(pi*k*(-t*uf + x))*cos(pi*k*(-t*vf + y))*cos(pi*t*w)**2/(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf &
+ x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)**2 - pi**2*k**2*cos(pi*k*(-t*uf + x))*cos(pi*k*(-t*vf + &
y))*cos(pi*t*w)/(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + &
1))/r0)*sin(pi*k*(-t*uf + x))*sin(pi*t*w)*cos(pi*k*(-t*vf + y))/(12*k)) + mu*(pi**2*k*w*(-r0 + &
r1)*((sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)/(r1*(r0/r1 + (-r0/r1 + &
1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)) + (-(sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)/(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1) + 1)/r0)*sin(pi*k*(-t*vf + &
y))*sin(pi*t*w)*cos(pi*k*(-t*uf + x))/2 - 3*pi*w*(-r0 + r1)*(-pi*k*(-r0/r1 + 1)*(sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)*sin(pi*k*(-t*uf + x))*cos(pi*k*(-t*vf + &
y))*cos(pi*t*w)/(r1*(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + &
1)**2) + pi*k*sin(pi*k*(-t*uf + x))*cos(pi*k*(-t*vf + y))*cos(pi*t*w)/(r1*(r0/r1 + (-r0/r1 + &
1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)) + (pi*k*(-r0/r1 + &
1)*(sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)*sin(pi*k*(-t*uf + &
x))*cos(pi*k*(-t*vf + y))*cos(pi*t*w)/(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + &
y))*cos(pi*t*w) + 1)**2 - pi*k*sin(pi*k*(-t*uf + x))*cos(pi*k*(-t*vf + y))*cos(pi*t*w)/(r0/r1 + &
(-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + &
1))/r0)*sin(pi*t*w)*cos(pi*k*(-t*uf + x))*cos(pi*k*(-t*vf + y))/4 + pi*w*(-r0 + r1)*(-pi*k*(-r0/r1 + &
1)*(sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)*sin(pi*k*(-t*vf + &
y))*cos(pi*k*(-t*uf + x))*cos(pi*t*w)/(r1*(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)**2) + pi*k*sin(pi*k*(-t*vf + y))*cos(pi*k*(-t*uf + &
x))*cos(pi*t*w)/(r1*(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + &
1)) + (pi*k*(-r0/r1 + 1)*(sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + &
1)*sin(pi*k*(-t*vf + y))*cos(pi*k*(-t*uf + x))*cos(pi*t*w)/(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)**2 - pi*k*sin(pi*k*(-t*vf + y))*cos(pi*k*(-t*uf + &
x))*cos(pi*t*w)/(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + &
1))/r0)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*sin(pi*t*w)/4 - w*(-r0 + &
r1)*(2*pi**2*k**2*(-r0/r1 + 1)**2*(sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + &
1)*sin(pi*k*(-t*uf + x))**2*cos(pi*k*(-t*vf + y))**2*cos(pi*t*w)**2/(r1*(r0/r1 + (-r0/r1 + &
1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)**3) + pi**2*k**2*(-r0/r1 + &
1)*(sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)*sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w)/(r1*(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)**2) - 2*pi**2*k**2*(-r0/r1 + 1)*sin(pi*k*(-t*uf + &
x))**2*cos(pi*k*(-t*vf + y))**2*cos(pi*t*w)**2/(r1*(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)**2) - pi**2*k**2*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + &
y))*cos(pi*t*w)/(r1*(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + &
1)) + (-2*pi**2*k**2*(-r0/r1 + 1)**2*(sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + &
1)*sin(pi*k*(-t*uf + x))**2*cos(pi*k*(-t*vf + y))**2*cos(pi*t*w)**2/(r0/r1 + (-r0/r1 + &
1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)**3 - pi**2*k**2*(-r0/r1 + &
1)*(sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)*sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w)/(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + &
y))*cos(pi*t*w) + 1)**2 + 2*pi**2*k**2*(-r0/r1 + 1)*sin(pi*k*(-t*uf + x))**2*cos(pi*k*(-t*vf + &
y))**2*cos(pi*t*w)**2/(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) &
+ 1)**2 + pi**2*k**2*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w)/(r0/r1 + (-r0/r1 + &
1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1))/r0)*sin(pi*k*(-t*vf + &
y))*sin(pi*t*w)*cos(pi*k*(-t*uf + x))/(4*k) - w*(-r0 + r1)*(2*pi**2*k**2*(-r0/r1 + &
1)**2*(sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)*sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*k*(-t*uf + x))*cos(pi*k*(-t*vf + y))*cos(pi*t*w)**2/(r1*(r0/r1 + &
(-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)**3) - pi**2*k**2*(-r0/r1 + &
1)*(sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)*cos(pi*k*(-t*uf + &
x))*cos(pi*k*(-t*vf + y))*cos(pi*t*w)/(r1*(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)**2) - 2*pi**2*k**2*(-r0/r1 + 1)*sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*k*(-t*uf + x))*cos(pi*k*(-t*vf + y))*cos(pi*t*w)**2/(r1*(r0/r1 + &
(-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)**2) + &
pi**2*k**2*cos(pi*k*(-t*uf + x))*cos(pi*k*(-t*vf + y))*cos(pi*t*w)/(r1*(r0/r1 + (-r0/r1 + &
1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)) + (-2*pi**2*k**2*(-r0/r1 + &
1)**2*(sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)*sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*k*(-t*uf + x))*cos(pi*k*(-t*vf + y))*cos(pi*t*w)**2/(r0/r1 + &
(-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)**3 + pi**2*k**2*(-r0/r1 + &
1)*(sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)*cos(pi*k*(-t*uf + &
x))*cos(pi*k*(-t*vf + y))*cos(pi*t*w)/(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + &
y))*cos(pi*t*w) + 1)**2 + 2*pi**2*k**2*(-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + &
y))*cos(pi*k*(-t*uf + x))*cos(pi*k*(-t*vf + y))*cos(pi*t*w)**2/(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf &
+ x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)**2 - pi**2*k**2*cos(pi*k*(-t*uf + x))*cos(pi*k*(-t*vf + &
y))*cos(pi*t*w)/(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + &
1))/r0)*sin(pi*k*(-t*uf + x))*sin(pi*t*w)*cos(pi*k*(-t*vf + y))/(4*k)))*((sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)/(r1*(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)) + (-(sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + &
y))*cos(pi*t*w) + 1)/(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + &
1) + 1)/r0) - (-pi*k*(-r0/r1 + 1)*(sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + &
1)*sin(pi*k*(-t*vf + y))*cos(pi*k*(-t*uf + x))*cos(pi*t*w)/(r1*(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf &
+ x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)**2) + pi*k*sin(pi*k*(-t*vf + y))*cos(pi*k*(-t*uf + &
x))*cos(pi*t*w)/(r1*(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + &
1)) + (pi*k*(-r0/r1 + 1)*(sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + &
1)*sin(pi*k*(-t*vf + y))*cos(pi*k*(-t*uf + x))*cos(pi*t*w)/(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)**2 - pi*k*sin(pi*k*(-t*vf + y))*cos(pi*k*(-t*uf + &
x))*cos(pi*t*w)/(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + &
1))/r0)*(-(uf - w*(-r0 + r1)*((sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + &
1)/(r1*(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)) + &
(-(sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)/(r0/r1 + (-r0/r1 + &
1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1) + 1)/r0)*sin(pi*k*(-t*vf + &
y))*sin(pi*t*w)*cos(pi*k*(-t*uf + x))/(4*k))**2/2 - (vf - w*(-r0 + r1)*((sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)/(r1*(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)) + (-(sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + &
y))*cos(pi*t*w) + 1)/(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + &
1) + 1)/r0)*sin(pi*k*(-t*uf + x))*sin(pi*t*w)*cos(pi*k*(-t*vf + y))/(4*k))**2/2 + w**2*(-r0 + &
r1)**2*((sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)/(r1*(r0/r1 + (-r0/r1 + &
1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)) + (-(sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)/(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1) + 1)/r0)**2*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + &
y))*sin(pi*t*w)**2*cos(pi*k*(-t*uf + x))*cos(pi*k*(-t*vf + y))/(32*k**2))/((sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)/(r1*(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)) + (-(sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + &
y))*cos(pi*t*w) + 1)/(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + &
1) + 1)/r0) - pi*w**2*(-r0 + r1)**2*((sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + &
1)/(r1*(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)) + &
(-(sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)/(r0/r1 + (-r0/r1 + &
1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1) + 1)/r0)**2*sin(pi*k*(-t*uf + &
x))**2*sin(pi*k*(-t*vf + y))*sin(pi*t*w)**2*cos(pi*k*(-t*vf + y))/(32*k) + pi*w**2*(-r0 + &
r1)**2*((sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)/(r1*(r0/r1 + (-r0/r1 + &
1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)) + (-(sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)/(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1) + 1)/r0)**2*sin(pi*k*(-t*vf + &
y))*sin(pi*t*w)**2*cos(pi*k*(-t*uf + x))**2*cos(pi*k*(-t*vf + y))/(32*k) - pi*w**2*(-r0 + &
r1)*((sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)/(r1*(r0/r1 + (-r0/r1 + &
1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)) + (-(sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)/(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1) + 1)/r0)*sin(pi*k*(-t*vf + y))*cos(pi*k*(-t*uf + &
x))*cos(pi*t*w)/(4*k) - w*(-r0 + r1)*((sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + &
1)*(pi*k*uf*(-r0/r1 + 1)*sin(pi*k*(-t*vf + y))*cos(pi*k*(-t*uf + x))*cos(pi*t*w) + pi*k*vf*(-r0/r1 + &
1)*sin(pi*k*(-t*uf + x))*cos(pi*k*(-t*vf + y))*cos(pi*t*w) + pi*w*(-r0/r1 + 1)*sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*sin(pi*t*w))/(r1*(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)**2) + (-pi*k*uf*sin(pi*k*(-t*vf + y))*cos(pi*k*(-t*uf + &
x))*cos(pi*t*w) - pi*k*vf*sin(pi*k*(-t*uf + x))*cos(pi*k*(-t*vf + y))*cos(pi*t*w) - &
pi*w*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*sin(pi*t*w))/(r1*(r0/r1 + (-r0/r1 + &
1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)) + (-(sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)*(pi*k*uf*(-r0/r1 + 1)*sin(pi*k*(-t*vf + &
y))*cos(pi*k*(-t*uf + x))*cos(pi*t*w) + pi*k*vf*(-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*cos(pi*k*(-t*vf + &
y))*cos(pi*t*w) + pi*w*(-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*sin(pi*t*w))/(r0/r1 &
+ (-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)**2 - &
(-pi*k*uf*sin(pi*k*(-t*vf + y))*cos(pi*k*(-t*uf + x))*cos(pi*t*w) - pi*k*vf*sin(pi*k*(-t*uf + &
x))*cos(pi*k*(-t*vf + y))*cos(pi*t*w) - pi*w*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + &
y))*sin(pi*t*w))/(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + &
1))/r0)*sin(pi*k*(-t*vf + y))*sin(pi*t*w)*cos(pi*k*(-t*uf + x))/(4*k) + w**2*(-r0 + &
r1)**2*((sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)/(r1*(r0/r1 + (-r0/r1 + &
1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)) + (-(sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)/(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1) + 1)/r0)*(-2*pi*k*(-r0/r1 + 1)*(sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)*sin(pi*k*(-t*vf + y))*cos(pi*k*(-t*uf + &
x))*cos(pi*t*w)/(r1*(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + &
1)**2) + 2*pi*k*sin(pi*k*(-t*vf + y))*cos(pi*k*(-t*uf + x))*cos(pi*t*w)/(r1*(r0/r1 + (-r0/r1 + &
1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)) + 2*(pi*k*(-r0/r1 + &
1)*(sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1)*sin(pi*k*(-t*vf + &
y))*cos(pi*k*(-t*uf + x))*cos(pi*t*w)/(r0/r1 + (-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + &
y))*cos(pi*t*w) + 1)**2 - pi*k*sin(pi*k*(-t*vf + y))*cos(pi*k*(-t*uf + x))*cos(pi*t*w)/(r0/r1 + &
(-r0/r1 + 1)*sin(pi*k*(-t*uf + x))*sin(pi*k*(-t*vf + y))*cos(pi*t*w) + 1))/r0)*sin(pi*k*(-t*uf + &
x))*sin(pi*k*(-t*vf + y))*sin(pi*t*w)**2*cos(pi*k*(-t*uf + x))*cos(pi*k*(-t*vf + y))/(32*k**2)

END FUNCTION VD2D_MMS_U_SRC_3

FUNCTION VD2D_MMS_V_SRC_3(X,Y,T)
IMPLICIT NONE
REAL(EB),INTENT (IN) :: X,Y,T
REAL(EB) :: VD2D_MMS_V_SRC_3

VD2D_MMS_V_SRC_3 = VD2D_MMS_U_SRC_3(Y,X,T)

END FUNCTION VD2D_MMS_V_SRC_3

END MODULE MANUFACTURED_SOLUTIONS



MODULE ONE_DIMENSIONAL_TURBULENCE

USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS, ONLY: LU_ERR,CHID
USE COMP_FUNCTIONS, ONLY: SHUTDOWN

IMPLICIT NONE

PRIVATE
PUBLIC :: GODT,SHEARLAYERODT

CONTAINS

!=======================================================================!
!                     This module is based on ...                       !
!                                                                       !
!                         lesodttoolsvalpha.f                           !
!                                by                                     !
!                           R. McDermott                                !
!                          Version Alpha 13                             !
!                    Last modified: 4 Jan 2006                          !
!                                                                       !
! Funding for this work provided through the Department of Energy Comp- !
! utational Science Graduate Fellowship (DE-FG02-97ER25308).            !
!                                                                       !
!=======================================================================!
!                                                                       !
! Conversion to Fortran 2008 and merged to FDS 18 May 2015. ~RJM        !
!                                                                       !
!=======================================================================!
! References:                                                           !
!                                                                       !
! [1] A.R. Kerstein, W.T. Ashurst, S. Wunsch, and V. Nilsen. One-dimen- !
! sional turbulence: vector formulation and application to free shear   !
! flows. J. Fluid Mech., 447:85-109, 2001.                              !
!                                                                       !
! [2] V. Nilsen. ODTN code description. Sandia National Laboratories    !
! internal Report.                                                      !
!                                                                       !
! [3] S. Wunsch and A. Kerstein. A model for layer formation in stably  !
! stratified turbulence. Phys. Fluids, 13(3):702, 2001.                 !
!                                                                       !
! [4] W.T. Ashurst and A.R. Kerstein. One-dimensional turbulence:       !
! Variable density formulation and application to mixing layers. Phys.  !
! Fluids, 17(2):025017, 2005.                                           !
!                                                                       !
! [5] A. Kerstein and S. Wunsch. Simulation of a stably stratified      !
! boundary layer using One-Dimensional Turbulence. Boundary Layer       !
! Meteorology, in press (2005).                                         !
!                                                                       !
! [6] A.R. Kerstein. One-dimensional turbulence: model formulation and  !
! application to homogeneous turbulence, shear flows, and buoyant       !
! stratified flows. J. Fluid Mech. 392:277-334, 1999.                   !
!=======================================================================!

SUBROUTINE GODT(U,V,W,RHO,PHI,MU,DTS,PA,LEDDY,Y0,L, &
                J0,KE,N_ODT,N_SCALAR,H,N_LMAX,G,GC, &
                BC,DEBUG,LE_SUP,PDFTAG,RET,C_EDDY,Z_EDDY,PAMAX,DTS_MIN)
IMPLICIT NONE

! INTENT IN:
INTEGER :: N_ODT,N_SCALAR,N_LMAX,BC,DEBUG,LE_SUP,PDFTAG
REAL(EB) :: DTS,H,G,GC
REAL(EB) :: C_EDDY,Z_EDDY,RET,PAMAX,DTS_MIN
! INTENT INOUT:
REAL(EB) :: U(N_ODT)
REAL(EB) :: V(N_ODT)
REAL(EB) :: W(N_ODT)
REAL(EB) :: RHO(N_ODT)
REAL(EB) :: PHI(N_SCALAR,N_ODT)
REAL(EB) :: MU(N_ODT)
! INTENT OUT:
INTEGER J0,KE,LEDDY
REAL(EB) :: PA,Y0,L

! LOCAL VARIABLES:
REAL(EB) :: D,L_MIN,L_MAX,L_MP
REAL(EB) :: U_EDDY(N_LMAX),V_EDDY(N_LMAX),W_EDDY(N_LMAX)
REAL(EB) :: RHO_EDDY(N_LMAX),MU_EDDY(N_LMAX)
REAL(EB) :: PHI_EDDY(N_SCALAR,N_LMAX)
REAL(EB) :: U_RHOK,V_RHOK,W_RHOK,RHO_K,RHO_KK
REAL(EB) :: P_U,P_V,P_W,S
REAL(EB) :: B_U,B_V,B_W,A
REAL(EB) :: C_U,C_V,C_W
INTEGER F(N_LMAX),K(N_LMAX)
INTEGER JL
INTEGER I,J,M,JMN_ODT
REAL(EB) :: E0,EF

!   N_ODT    = NUMBER OF ODT PODINTS ON THE ODT LINE
!   N_LMAX   = NUMBER OF ODT PODINTS IN MAX EDDY (> 18 DUE TO L_MP)
!   N_SCALAR = NUMBER OF SCALARS, E.G. TEMPERATURE, MIXTURE FRACTION
!   BC       = BOUNDARY CONDITION: 1=PERIODIC, 2=BOUNDED
!   RNTAG    = WHICH RANDOM NUMBER GENERATOR ?
!      0 = F77 STD "RAND", 1 = TOM DREEBEN "RNG"
!   DEBUG    = DEBUG LEVEL: 0 = NO CHECKS, 1 = CHECK ENERGY CONSERVATION
!   LE_SUP   = LARGE-EDDY SUPPRESSION: 0 = NONE, 1 = SUBINTERVAL, 2 = MEDIAN
!   PDFTAG   = WHICH SAMPLE LENGTH SCALE PDF? 1 = VEBJORN NILSEN, 2 = LEM (H. SHIHN)
!   RET      = TURBULENCE REYNOLDS NUMBER PARAMETER FOR LEM PDF (H. SHIHN)
!   DTS      = SAMPLING PERIOD (SECONDS)
!   H        = ODT GRID SPACING (METERS)
!   G        = ACCELERATION DUE TO GRAVITY (E.G. -9.8 M/S^2)
!   GC       = ABS(G) FOR BAROCLINIC TERM (I.E. HORIZONTAL LINES)
!   NOTE: G AND GC SHOULD SATISFY THE ORTHOGONALITY REQUIREMENT
!      THAT GC^2 + G^2 = (9.8)^2, FOR EXAMPLE, DEPENDING ON YOUR UNITS.
!   MU       = DYNAMIC VISCOSITY (N*S/M^2)
!   U,V,W    = VELOCITY (M/S^2)
!   RHO      = DENSITY (KG/M^3)
!   PHI      = ARBITRARY SCALARS NOT SUBJECT TO KERNEL TRANSFORM
!
!   D        = DOMAIN LENGTH (METERS)
!   L_MIN    = MINIMUM EDDY SIZE, = 6*H DUE TO DISCRETE TRIPLET MAP
!   L_MAX    = MAXIMUM EDDY SIZE, USER SPECIFIED
!   L_MP     = "MOST PROBABLE" EDDY SIZE IN SAMPLE DISTRIBUTION
!
!   UI_EDDY  = CANDIDATE EDDY VELOCITIES
!   RHO_EDDY = CANDIDATE EDDY DENSITY
!   MU_EDDY  = CANDIDATE EDDY VISCOSITY
!   PHI_EDDY = USER SPECIFIED SCALARS
!
!   C_EDDY   = EDDY RATE CONSTANT, ROUGHLY 1.4*DSQRT(54)/3 = 3.43
!      NOTE: C = 3.68 IS EMPIRICAL VALUE FROM SHEAR LAYER IN [1].
!   Z_EDDY   = VISCOUS CUTOFF, ROUGHLY 3.0
!   ALPHA    = FRACTION OF MAX ENERGY EXCHANGE IN PRESSURE-SCRAMBLING
!      MODEL, ONLY USED IF FLAVOR = 0 (RECOMMEND FLAVOR = 1)
!   PA       = ACCEPTANCE PROBABILITY, ADJUST SAMPLING PERIOD (DTS)
!      TO KEEP PA \APPROX 0.1
!   PAMAX    = MAXIMUM ALLOWABLE ACCEPTANCE PROBABILITY, SCALES DTS WITHIN GODT
!   Y0       = EDDY STARTING LOCATION
!   L        = SAMPLED CONTINUOUS EDDY SIZE
!   UI_RHOK  = EDDY VELOCITY SCALE
!   C_I      = KERNEL AMPLITUDE
!   F(J)     = DISCRETE TRIPLET MAP PERMUTATION VECTOR
!   K(J)     = DISCRETE KERNEL PERMUTATION VECTOR
!   FLAVOR   = 1 SPEEDS CODE FOR ALPHA = 2/3 CASE (RECOMMENDED)
!   J0       = INDEX OF CANDIDATE EDDY STARTING LOCATION (~Y0)
!   JL       = INDEX OF CANDIDATE EDDY ENDING LOCATION (~Y0+L)
!   KE       = NUMBER OF PODINTS IN CANDIDATE EDDY DISCRETE EDDY LENGTH, LD = KE-1
!   I,J,M,JMN_ODT   = COUNTING INDICES
!   LEDDY    = IF LEDDY.EQ.1 IMPLEMENT EDDY
!   E_INITIAL, E_FINAL = USED FOR ENERGY CONSERVATION CHECK
!
!-----------------------------------------------------------------------
!   C_EDDY = 3.68 FOR JFM 2001 SHEAR LAYER SIMULATION
!   PARAMETER (C_EDDY = 3.68D0)
!   C_EDDY = 5.60 FOR PHYS. FLUIDS 2005 MIXING LAYER
!   PARAMETER (C_EDDY = 5.6D0)
!   PARAMETER (Z_EDDY = 0.14D0)
INTEGER, PARAMETER :: FLAVOR = 1
REAL(EB), PARAMETER :: ALPHA = 2._EB/3._EB

!   CHECK TO MAKE SURE L_MAX IS GREATER THAN MOST PROBABLE EDDY SIZE
IF (N_LMAX < 18) THEN
   WRITE(LU_ERR,*) 'ERROR: L_MAX < L_MP'
   STOP
ENDIF
!   CHECK TO MAKE SURE PERIODIC BCS ARE NOT USED WITH GRAVITY ON
IF (BC==1 .AND. ABS(G)>1.E-10_EB) THEN
   WRITE(LU_ERR,*) 'ERROR: CANNOT USE GRAVITY WITH PERIODIC BC'
   STOP
ENDIF

D = N_ODT*H
L_MIN = 6._EB*H
L_MAX = N_LMAX*H
L_MP = 3._EB*L_MIN

!       SELECT A RANDOM EDDY LOCATION (Y0) AND LENGTH (L)
!       J0 = INDEX LOCATION OF PHYSICAL POSITION, Y0
!       JL = INDEX LOCATION OF PHYSICAL POSITION, Y0 + L
!       KE = NUMBER OF DISCRETE PODINTS IN THE EDDY, MUST BE MULTIPLE OF 3!!
!          |<------------- L ------------->|
!          4   5   6   7   8   9   10  11  12
!        | O | O | O | O | O | O | O | O | O |   FOR THIS EXAMPLE, KE = 9
!          ^                               ^
!          J0                              JL    JL = J0 + KE - 1

CALL SELECTEDDY(J0,KE,Y0,L,0._EB,D,L_MIN,L_MAX,L_MP,H,1,BC,PDFTAG,RET,N_LMAX)

JL = J0+KE-1
I = 0
JL_IF: IF (JL < N_ODT) THEN
   DO J = J0,JL
      I = I+1
      U_EDDY(I) = U(J)
      V_EDDY(I) = V(J)
      W_EDDY(I) = W(J)
      RHO_EDDY(I) = RHO(J)
      MU_EDDY(I) = MU(J)
      DO M = 1,N_SCALAR
         PHI_EDDY(M,I) = PHI(M,J)
      ENDDO
   ENDDO
ELSE JL_IF
   ! should only possibly use this if bc == 1 (periodic)
   IF (BC/=1) THEN
      WRITE(LU_ERR,*) 'WARNING: EDDY EXTENDS OUTSIDE BOUNDARY'
      WRITE(LU_ERR,*) 'J0 = ',J0, ' KE = ',KE
      WRITE(LU_ERR,*) 'JL = ',JL, ' N_ODT = ',N_ODT
      IF (DEBUG==1) THEN
         STOP
      ELSE
         J0 = J0-1
         JL = JL-1
      ENDIF
   ENDIF
   DO J = J0,N_ODT
      I = I+1
      U_EDDY(I) = U(J)
      V_EDDY(I) = V(J)
      W_EDDY(I) = W(J)
      RHO_EDDY(I) = RHO(J)
      MU_EDDY(I) = MU(J)
      DO M = 1,N_SCALAR
         PHI_EDDY(M,I) = PHI(M,J)
      ENDDO
   ENDDO
   DO J = N_ODT+1,JL
      JMN_ODT = J-N_ODT
      I = I+1
      U_EDDY(I) = U(JMN_ODT)
      V_EDDY(I) = V(JMN_ODT)
      W_EDDY(I) = W(JMN_ODT)
      RHO_EDDY(I) = RHO(JMN_ODT)
      MU_EDDY(I) = MU(JMN_ODT)
      DO M = 1,N_SCALAR
         PHI_EDDY(M,I) = PHI(M,JMN_ODT)
      ENDDO
   ENDDO
ENDIF JL_IF

! DECIDE IF EDDY SHOULD BE IMPLEMENTED
CALL EVALUATEEDDY(LEDDY,U_RHOK,V_RHOK,W_RHOK, &
   RHO_K,RHO_KK, F,K, PA, KE,H, DTS, &
   U_EDDY,V_EDDY,W_EDDY,RHO_EDDY,MU_EDDY, &
   L_MIN,L_MAX,L_MP,L,0.D0,D,C_EDDY,Z_EDDY,G,GC,BC,PDFTAG,RET, &
   P_U,P_V,P_W,S,A,PAMAX,DTS_MIN)

! LARGE-EDDY SUPPRESSION
IF (LEDDY==1 .AND. KE>18) THEN
   SELECT CASE(LE_SUP)
      CASE(1); CALL SUBINTERVALEVALUATION(LEDDY,PA,U_EDDY,V_EDDY,W_EDDY,RHO_EDDY,MU_EDDY, KE,H, Z_EDDY,G,GC)
      !CASE(2); CALL MEDIANMETHOD(LEDDY,PA,U_EDDY,V_EDDY,W_EDDY,RHO_EDDY,MU_EDDY, KE,H, Z_EDDY)
   END SELECT
ENDIF

! IF LEDDY IS TRUE, IMPLEMENT THE EDDY AND APPLY KERNEL TRANSFORMATION
LEDDY_IF: IF (LEDDY==1) THEN

   IF (FLAVOR==0 .AND. ABS(ALPHA)<TWO_EPSILON_EB) THEN
      C_U=0._EB
      C_V=0._EB
      C_W=0._EB
   ELSE
      CALL KERNELAMP(C_U,C_V,C_W, &
         RHO_K, &
         FLAVOR,ALPHA,KE,H,G, &
         P_U,P_V,P_W,S)
   ENDIF

   B_U = -A*C_U
   B_V = -A*C_V
   B_W = -A*C_W

   IF (DEBUG==1) THEN
      ! CHECK FOR ENERGY CONSERVATION
      ! TOTAL INITIAL ENERGY
      E0 = 0._EB
      DO J = J0,JL
         I = J
         IF (BC==1 .AND. I>N_ODT) I = I-N_ODT
         E0 = E0 + 0.5*RHO(I)*(U(I)**2+V(I)**2+W(I)**2) - RHO(I)*G*(J*H)
      ENDDO
   ENDIF

   ! APPLY STIRRING EVENT
   DO J = 1,KE
      I = J0+J-1
      IF (BC==1 .AND. I>N_ODT) I = I-N_ODT
      U(I) = U_EDDY(F(J)) + B_U*ABS(K(J)) + C_U*K(J)
      V(I) = V_EDDY(F(J)) + B_V*ABS(K(J)) + C_V*K(J)
      W(I) = W_EDDY(F(J)) + B_W*ABS(K(J)) + C_W*K(J)
      RHO(I) = RHO_EDDY(F(J))
      DO M = 1,N_SCALAR
         PHI(M,I) = PHI_EDDY(M,F(J))
      ENDDO
      MU(I) = MU_EDDY(F(J))
   ENDDO

   IF (DEBUG==1 .AND. KE>6) THEN
      ! TOTAL FINAL ENERGY
      EF = 0._EB
      DO J = J0,JL
         I = J
         IF (BC==1 .AND. I>N_ODT) I = I-N_ODT
         EF = EF + 0.5*RHO(I)*(U(I)**2+V(I)**2+W(I)**2) - RHO(I)*G*(J*H)
      ENDDO
      IF ( ABS((EF - E0)/E0) > 1.E-10_EB ) THEN
         WRITE(LU_ERR,*) 'ENERGY CONSERVATION VIOLATED!'
         WRITE(LU_ERR,*) 'INITIAL ENERGY = ',E0
         WRITE(LU_ERR,*) 'FINAL ENERGY = ',EF
         WRITE(LU_ERR,*) C_U,C_V,C_W,G,J0,KE,JL
         STOP
      ENDIF
   ENDIF

ENDIF LEDDY_IF

END SUBROUTINE GODT


SUBROUTINE SELECTEDDY(J0,KE,Y0,L,R0,RD,L_MIN,L_MAX,L_MP,H,J1,BC,PDFTAG,RET,N_LMAX)
IMPLICIT NONE

! INTENT OUT:
INTEGER :: J0,KE
REAL(EB) :: Y0,L
! INTENT IN:
REAL(EB) :: R0,RD,L_MIN,L_MAX,L_MP,H,RET
INTEGER :: J1,BC,PDFTAG,N_LMAX

! OUTPUTS:
! J0 = FIRST GRID INDEX OF CANDIDATE EDDY
! KE = NUMBER OF PODINTS IN CANDIDATE EDDY, MUST BE MULTIPLE OF 3
! Y0 = STARTING LOCATION OF EDDY SAMPLED FROM CONTINUOUS DISTRIBUTION
! L = PHYSICAL SIZE OF EDDY SELECTED FROM DISTRIBUTION; THIS IS
!     USED LATER IN "EVALUATE_EDDY" TO CORRECT FOR THE DISCRETE IMPLEMENTATION
! INPUTS:
! R0 = PHYSICAL LOCATION OF LOWER END OF THE DOMAIN (E.G. USUALLY R0 = 0)
! RD = PHYSICAL LOCATION OF UPPER END OF DOMAIN
! L_MIN = PHYSICAL MINIMUM EDDY SIZE
! L_MAX = PHYSICAL MAXIMUM EDDY SIZE
! L_MP = MOST PROBABLE EDDY SIZE
! H  = PHYSICAL GRID SPACING
! J1 = FIRST GRID INDEX ON ODT DOMAIN
! BC: 1 = PERIODIC, 2 = BOUNDED
! RNTAG: 0 = F77 STD RN = RAND, 1 = TOM DREEBEN RNG(RN)
! PDFTAG: 1 = USE V. NILSEN PDF, 2 = USE LEM PDF (H. SHIHN)
! RET = TURBULENCE REYNOLDS NUMBER PARAMETER FOR LEM PDF
!
! LOCAL VARIABLES:
REAL(EB) :: RN
! GAMMA = -8/5, SEE (13) IN [2]
! RN = PSUEDO-RANDOM NUMBER IN [0,1]
! FIVETHIRDS,THREEFIFTHS = HOPEFULLY SELF EXPLANATORY
REAL(EB) :: P,REFAC,LK,NK ! HARMAN'S DISTRIBUTION CONSTANT-LEM

!----------------------------------------------------------------------------
!  J0 = INDEX LOCATION OF PHYSICAL POSITION, Y0
!  JL = INDEX LOCATION OF PHYSICAL POSITION, Y0 + L
!  KE = NUMBER OF DISCRETE PODINTS IN THE EDDY, MUST BE MULTIPLE OF 3!!
!  |<--------------- L --------------->|
!  | 4   5   6   7   8   9   10  11  12|
!  | O | O | O | O | O | O | O | O | O |   FOR THIS EXAMPLE, KE = 9
!  | ^                               ^ |
!    J0                              JL    JL = J0 + KE - 1

REAL(EB), PARAMETER :: GAMMA = -1.6_EB
REAL(EB), PARAMETER :: FIVETHIRDS = 5._EB/3._EB, THREEFIFTHS = 3._EB/5._EB

!  SELECT THE LENGTH SCALE, EQN (14) IN [2]
CALL RANDOM_NUMBER(RN)
SELECT CASE(PDFTAG)
   CASE(0)
      L = RN*(1./(L_MAX-L_MIN))
   CASE(1)
      L = L_MP/(((1./GAMMA)* &
          LOG(EXP(GAMMA*((L_MP/L_MIN)**(FIVETHIRDS)))+ &
          (EXP(GAMMA*((L_MP/L_MAX)**(FIVETHIRDS)))- &
          EXP(GAMMA*((L_MP/L_MIN)**(FIVETHIRDS))))*RN))**(THREEFIFTHS))
   CASE(2)
      P = 1.334_EB
      ! RET = 80._EB
      REFAC = RET**((3._EB/P)-1._EB)
      LK = 0.4_EB/(RD-R0)
      NK = (LK)/(RET**(1._EB/P))
      L = (LK)*((RN*(1._EB-(RET**((3._EB/P)-1._EB)))+(RET**((3._EB/P)-1._EB)))**(1._EB/(P-3._EB)))
END SELECT

!  CHOOSE HOW MANY CELLS THE EDDY SHOULD HAVE: SHRINK OR EXPAND IT TO AN
!  INTEGER DIVISIBLE BY THREE
KE=3*NINT((L+H)/(3*H))
KE=MIN(N_LMAX,KE)
KE=KE-MOD(KE,3)
KE=MAX(6,KE)

!  SELECT POSITION FROM A UNIFORM DISTRIBUTION
CALL RANDOM_NUMBER(RN)

IF (BC==1) THEN
!  PERIODIC
   Y0 = R0 + RN*(RD-R0)
ELSE IF (BC==2) THEN
!  BOUNDED
   Y0 = R0 + RN*(RD-R0-KE*H)
ENDIF

!   THE SUBTRACTION OF 1.E-10 HANDLES THE CASE WHERE RN = 1.D0, WHICH
!   WOULD PUT J0+KE-1 = N_ODT + 1, AND (BELIEVE IT OR NOT) DOES OCCUR.
J0 = J1 + INT((Y0-R0-1.E-10_EB)/H)

END SUBROUTINE SELECTEDDY


SUBROUTINE EVALUATEEDDY(LEDDY,U_RHOK,V_RHOK,W_RHOK, &
                        RHO_K,RHO_KK, F,K, PA, KE,H, DTS, &
                        U,V,W,RHO,MU, &
                        L_MIN,L_MAX,L_MP,L,R0,RD, &
                        C_EDDY,Z_EDDY,G,GC,BC,PDFTAG,RET, &
                        P_U,P_V,P_W,S,A,PAMAX,DTS_MIN)
IMPLICIT NONE

! INTENT IN:
INTEGER :: KE,BC,PDFTAG
REAL(EB) :: U(KE),V(KE),W(KE)
REAL(EB) :: RHO(KE),MU(KE),RET
REAL(EB) :: H,L_MIN,L_MAX,L_MP,L,R0,RD
REAL(EB) :: C_EDDY,Z_EDDY,G,GC,PAMAX,DTS_MIN

! INTENT OUT:
REAL(EB) :: U_RHOK,V_RHOK,W_RHOK
REAL(EB) :: U_RHOJ,V_RHOJ,W_RHOJ
REAL(EB) :: RHO_J,RHO_K,RHO_KK,RHO_JK,PA
REAL(EB) :: P_U,P_V,P_W,S,A
INTEGER :: F(KE),K(KE),LEDDY

! INTENT INOUT:
REAL(EB) :: DTS

! OUTPUTS:
! LEDDY = INTEGER, IF LEDDY=1 EDDY IS TO BE IMPLEMENTED
! U_RHOK = SCALAR, MEASURE OF EDDY U-MOMENTUM CHANGE (KG/M^2*S)
! V_RHOK = SCALAR, MEASURE OF EDDY V-MOMENTUM CHANGE (KG/M^S*S)
! W_RHOK = SCALAR, MEASURE OF EDDY W-MOMENTUM CHANGE (KG/M^2*S)
! RHO_K = SCALAR, MEASURE OF DENSITY VARIATION (KG/M^3)
! RHO_KK = SCALAR, MEASURE OF STD DEV OF DENSITY VARIATION (KG/M^3)
! MU_AVE = AVERAGE DYNAMIC VISCOSITY
! F(KE) = TRIPLET MAP PERMUATION VECTOR
! K(KE) = KERNEL PERMUTATION VECTOR
! PA = ACCEPTANCE PROBABILITY
! INPUTS:
! U(KE) = X VELOCITY (M/S)
! V(KE) = Y VELOCITY (M/S)
! W(KE) = Z VELOCITY (M/S)
! RHO(KE) = DENSITY (KG/M^3)
! MU(KE) = DYNAMIC VISCOSITY (N*S/M^2)
! KE = NUMBER OF PODINTS IN CANDIDATE EDDY, SHOULD BE MULTIPLE OF 3
! H = GRID SPACING (M)
! L_MIN = MIMIMUM EDDY SIZE (M)
! L_MAX = MAXIMUM EDDY SIZE (M)
! L = CONTINUOUS EDDY SIZE (M)
! R0 = DOMAIN LOWER EXTENT (M)
! RD = DOMAIN UPPER EXTENT (M)
! CE = EDDY RATE CONSTANT \APPROX CE = 1.4/3*SQRT(54)
! Z = VISCOUS CUTOFF \APPROX Z = 3.0
! G = GRAVITATIONAL ACCELERATION (E.G. -9.8 M/S^2)
! GC = ABS(G) FOR BAROCLINIC, SHOULD SATISFY GC^2 + G2 = (9.8)^2
! BC = TAG FOR BOUNDARY CONDITION: 1 = PERIODIC, 2 = WALLS
! RNTAG = TAG FOR RANDOM NUMBER GENERATOR: 0 = RAND, 1 = RNG
! PDFTAG = TAG FOR LENGTH SCALE PDF: 1 = V. NILSEN, 2 = LEM (H. SHIHN)

! LOCAL VARIABLES:
REAL(EB) :: LNORM,FL,FY0
INTEGER M,I
REAL(EB) :: RE_U,RE_V,RE_W,GR,MU_BAR
REAL(EB) :: LD,DISCOR,LAMBDA,ENERGY
REAL(EB) :: RN
REAL(EB) :: LL,LLL,J(KE)
! FOR HARMAN'S LEM DISTRIBUTION
REAL(EB) :: P,LK,NK
! FOR HARMAN'S BAROCLINIC TERM
REAL(EB) :: RI_RE2,RHO_DIFF,RHO_LEFT,RHO_RIGHT

! GAMMA = -8/5
! FIVETHIRDS = 5/3
! LNORM = PRE-FACTOR IN LENGTH-SCALE PROBABILITY DENSITY
! FL = LENGTH-SCALE PROBABILITY DENSITY (CONTINUOUS DISTRIBUTION)
! FY0 = LOCATION PROBABILITY DENSITY (CONTINUOUS DISTRIBUTION)
! RN = RANDOM NUMBER ON DINTERVAL [0,1]
! M = KE/3
! J = COUNTING INDEX
! RE_I = COMPONENT EDDY REYNOLDS NUMBER
! GR = EDDY GRASHOF NUMBER (NON-DIMENSIONAL POTENTIAL ENERGY)
! LD = DISCRETE EDDY PHYSICAL LENGTH (METERS)
! DISCOR = CORRECTION BECAUSE THE MEAN SQUARE DISPLACEMENT OF A
!          DISCRETE EDDY IS DIFFERENT FROM A CONTINUOUS EDDY
! LAMBDA = EDDY EVENT RATE DENSITY [=] EVENTS/LOCATION/LENGTH/TIME
! ENERGY = ARGUMENT OF SQUARE-ROOT TERM IN RATE DENSITY
!         A NEGATIVE ENERGY INDICATES THE EDDY IS NOT ENERGETIC
!        ENOUGH TO TURN OVER AND IS HENCE NOT IMPLEMENTED
! RN = RANDOM NUMBER ON INTERVAL [0,1]
! ======================================================================
REAL(EB), PARAMETER :: GAMMA = -1.6_EB
REAL(EB), PARAMETER :: FIVETHIRDS = 5._EB/3._EB
REAL(EB), PARAMETER :: B = 2.2_EB ! HARMAN'S EMPIRICALLY DETERMINED CONSTANT B=0.05 CONVERTS TO 2.2 FOR VAR DEN
REAL(EB), PARAMETER :: EPS_LOC=1.E-10_EB

LEDDY = 1

! EVALUATE LOCATION PROBABILITY DENSITY AT L
SELECT CASE(PDFTAG)
   CASE(0)
      FL = 1._EB/(L_MAX-L_MIN) ! uniform distribution
   CASE(1)
      ! LNORM can be precomputed...
      LNORM =(-(FIVETHIRDS)*GAMMA*(L_MP**(FIVETHIRDS)))/ &
             (EXP(GAMMA*((L_MP/L_MAX)**(FIVETHIRDS)))- &
              EXP(GAMMA*((L_MP/L_MIN)**(FIVETHIRDS))))
      FL = LNORM*(L**(-8._EB/3._EB))*EXP(GAMMA*((L_MP/L)**(FIVETHIRDS)))
   CASE(2)
      !RET = 80.D0
      P = 1.334D0
      LK = 0.4D0/(RD-R0)
      NK = (LK)/(RET**(1.D0/P))
      FL =  ((3-P)/(LK*(((LK/NK)**(3-P))-1))) ! LEM PDF
END SELECT
! EVALUATE LOCATION PROBABILITY DENSITY AT Y0 (UNIFORM)
IF (BC==1) THEN
   FY0 = 1._EB/(RD-R0)
ELSE IF (BC==2) THEN
   FY0 = 1._EB/(RD-L-R0)
ENDIF

! IF KE < 6 USE A DIFFERENT PERMUTATION (SEE SUBINTERVAL EVAL)
IF (KE==4) THEN
   F(1) = 1
   F(2) = 4
   F(3) = 1
   F(4) = 4
ELSEIF (KE==5) THEN
   F(1) = 1
   F(2) = 4
   F(3) = 3
   F(4) = 2
   F(5) = 5
ELSE
   IF (MOD(KE,3)/=0) THEN
      CALL SHUTDOWN('ERROR: eddy not divisible by 3!') ; RETURN
   ENDIF
   M = KE/3
   DO I = 1,M
      F(I) = 3*I - 2
   ENDDO
   DO I = M+1,2*M
      F(I) = 2*(KE-1) - 3*I + 4
   ENDDO
   DO I = 2*M+1,KE
      F(I) = 3*I - 2 - 2*(KE-1)
   ENDDO
ENDIF

! INTEGRAL EDDY PARAMETERS
U_RHOK = 0._EB
V_RHOK = 0._EB
W_RHOK = 0._EB
U_RHOJ = 0._EB
V_RHOJ = 0._EB
W_RHOJ = 0._EB
RHO_J  = 0._EB
RHO_K  = 0._EB
RHO_KK = 0._EB
RHO_JK = 0._EB

DO I = 1,KE
   K(I) = I - F(I)
   J(I) = ABS(K(I))
   U_RHOK = U_RHOK + RHO(F(I))*U(F(I))*K(I)
   V_RHOK = V_RHOK + RHO(F(I))*V(F(I))*K(I)
   W_RHOK = W_RHOK + RHO(F(I))*W(F(I))*K(I)
   U_RHOJ = U_RHOJ + RHO(F(I))*U(F(I))*J(I)
   V_RHOJ = V_RHOJ + RHO(F(I))*V(F(I))*J(I)
   W_RHOJ = W_RHOJ + RHO(F(I))*W(F(I))*J(I)
   RHO_J  = RHO_J  + RHO(F(I))*J(I)
   RHO_K  = RHO_K  + RHO(F(I))*K(I)
   RHO_KK = RHO_KK + RHO(F(I))*K(I)*K(I)
   RHO_JK = RHO_JK + RHO(F(I))*J(I)*K(I)
ENDDO

LL = (REAL(KE,EB)-1._EB)**2
LLL = (REAL(KE,EB)-1._EB)**3
U_RHOK = U_RHOK/LL
V_RHOK = V_RHOK/LL
W_RHOK = W_RHOK/LL
U_RHOJ = U_RHOJ/LL
V_RHOJ = V_RHOJ/LL
W_RHOJ = W_RHOJ/LL
RHO_J = RHO_J/LL
RHO_K = RHO_K/LL
RHO_KK = RHO_KK/LLL
RHO_JK = RHO_JK/LLL

! ASHURST VARIABLE DENSITY PARAMETERS
LD = (REAL(KE,EB)-1._EB)*H
A = RHO_K/RHO_J
P_U = (U_RHOK - A*U_RHOJ)*LD**2
P_V = (V_RHOK - A*V_RHOJ)*LD**2
P_W = (W_RHOK - A*W_RHOJ)*LD**2
S = (0.5_EB*(A**2+1._EB)*RHO_KK - A*RHO_JK)*LD**3

! AVERAGE DYNAMIC VISCOSITY
! SEE E.G., ECHEKKI ET AL. COMBUSTION AND FLAME (2001)
MU_BAR = 0.
DO I = 1,KE
   MU_BAR = MU_BAR+1._EB/MU(I)
ENDDO
MU_BAR = MU_BAR/KE
MU_BAR = 1._EB/MU_BAR

! EDDY REYNOLDS NUMBERS
! NOTE THAT MY DISCRETE EDDY LENGTH IS DIFFERENT FROM VEBJORN'S
IF (MU_BAR<EPS_LOC) THEN
   MU_BAR = 1._EB
   Z_EDDY = 0._EB
ENDIF
RE_U = U_RHOK*LD/MU_BAR
RE_V = V_RHOK*LD/MU_BAR
RE_W = W_RHOK*LD/MU_BAR

! EDDY GRASHOF NUMBER
GR = (2._EB*LD**3*RHO_KK/MU_BAR**2)*G*RHO_K
print *, GR

! HARMANJEET SHIHN (SUNY, BUFFALO) BAROCLINIC TORQUE TERM
IF (GC>EPS_LOC) THEN

   ! THE DELTA RHO TERM IN THE RICHARDSON NUMBER IS COMPUTED USING
   ! EQUATION (2.4) FROM REFERCNE [6], EXCEPT HARMAN LEAVES OUT THE
   ! FACTOR OF TWO IN THE FINAL DIFFERENCE.
   RHO_LEFT = 0._EB
   RHO_RIGHT = 0._EB
   DO I = 1,KE/2
      RHO_LEFT = RHO_LEFT + RHO(I)
   ENDDO
   RHO_LEFT = 2._EB*RHO_LEFT/KE
   DO I = KE/2+1,KE
      RHO_RIGHT = RHO_RIGHT + RHO(I)
   ENDDO
   RHO_RIGHT = 2._EB*RHO_RIGHT/KE
   RHO_DIFF = ABS(RHO_RIGHT-RHO_LEFT)
   RI_RE2 = 2._EB*LD**3*RHO_KK/MU_BAR**2*GC*RHO_DIFF
ELSE
   RI_RE2 = 0._EB
ENDIF

! EVENT RATE DENSITY
! ENERGY = RE_U**2+RE_V**2+RE_W**2 + GR + B*Ri*Re**2 - Z**2
ENERGY = RHO_KK*LD/(2._EB*S*MU_BAR**2)*(P_U**2 + P_V**2 + P_W**2) + GR + B*RI_RE2 - Z_EDDY**2

IF (ENERGY>TWO_EPSILON_EB) THEN
   LAMBDA = C_EDDY*MU_BAR/(RHO_KK*LD**2)*SQRT(ENERGY)
ELSE
   PA = 0._EB
   LEDDY = 0
   RETURN
ENDIF

! CONTINUOUS--DISCRETE CORRECTION FACTOR NEEDED TO MAKE THE
! CONTINUOUS RANDOM-WALK DIFFUSIVITY EQUAL TO THE DISCRETE
! RANDOM-WALK DIFFUSIVITY
DISCOR = (L/LD)**2 * REAL((KE-1)**3,EB)/REAL(KE**3 - 3*KE**2,EB)
LAMBDA = DISCOR*LAMBDA

! ACCEPTANCE PROBABILITY
PA=LAMBDA*DTS/(FY0*FL)

! CHECK THAT PA IS POSITIVE
IF (PA<0._EB) THEN
   CALL SHUTDOWN('ERROR: acceptance probability < 0')
   RETURN
ENDIF

! SAMPLING PERIOD CONTROLLER:
IF (PA > PAMAX) THEN
   DTS = MAX(DTS*PAMAX/PA,DTS_MIN)
   PA = PAMAX
   LEDDY = 0
   RETURN
END IF

! BERNOULLI TRIAL: REJECT CANDIDATE EDDY IF RANDOM PROBABILITY GREATER THAN ACCEPTANCE PROBABILITY
CALL RANDOM_NUMBER(RN)
IF (RN > PA) THEN
   LEDDY = 0
ENDIF

END SUBROUTINE EVALUATEEDDY


SUBROUTINE SUBINTERVALEVALUATION(LEDDY,PA,U,V,W,RHO,MU, KE,H, Z,G,GC)
IMPLICIT NONE

! INTENT IN:
INTEGER :: KE
REAL(EB) :: U(KE),V(KE),W(KE),RHO(KE),MU(KE)
REAL(EB) :: H,Z,G,GC

! INTENT INOUT:
INTEGER :: LEDDY
REAL(EB) :: PA

! OUTPUTS:
! LEDDY = LOGICAL, IF TRUE EDDY IS TO BE IMPLEMENTED
! INPUTS:
! U(KE) = X VELOCITY (M/S)
! V(KE) = Y VELOCITY (M/S)
! W(KE) = Z VELOCITY (M/S)
! RHO(KE) = DENSITY (KG/M^3)
! MU(KE) = DYNAMIC VISCOSITY (N*S/M^2)
! KE = NUMBER OF POINTS IN CANDIDATE EDDY, SHOULD BE MULTIPLE OF 3
! H = GRID SPACING (M)
! Z = VISCOUS CUTOFF \APPROX Z = 3.0
! G = GRAVITATIONAL ACCELERATION (E.G. -9.8 M/S^2)
! GC = ABS(G) FOR BAROCLINIC TERM (I.E. HORIZNONTAL LINES)
! NOTE: G AND GC SHOULD SATISFY THE ORTHOGONALITY REQUIREMENT
! THAT GC^2 + G^2 = (9.8)^2, FOR EXAMPLE, DEPENDING ON YOUR UNITS.

! LOCAL VARIABLES:
INTEGER :: M,MM,I,J,J_SUB
REAL(EB) :: RE_U,RE_V,RE_W,GR,MU0
REAL(EB) :: LD,LS,ENERGY
REAL(EB) :: U_SUB(KE/3)
REAL(EB) :: V_SUB(KE/3)
REAL(EB) :: W_SUB(KE/3)
REAL(EB) :: RHO_SUB(KE/3)
INTEGER F(KE/3),K(KE/3),JJ(KE/3),FFF(KE)
REAL(EB) :: U_RHOK,V_RHOK,W_RHOK,RHO_K,RHO_KK
REAL(EB) :: U_RHOJ,V_RHOJ,W_RHOJ,RHO_J,RHO_JK
REAL(EB) :: LL,LLL,A,S,P_U,P_V,P_W
! FOR HARMAN'S BAROCLINIC TERM
REAL(EB) :: RI_RE2,RHO_DIFF,RHO_LEFT,RHO_RIGHT
REAL(EB), PARAMETER :: B=2.2_EB
REAL(EB), PARAMETER :: EPS_LOC=1.E-10_EB
! M = KE/3
! J = COUNTING INDEX
! MU_AVE = AVERAGE DYNAMIC VISCOSITY
! RE_I = COMPONENT EDDY REYNOLDS NUMBER
! GR = EDDY GRASHOF NUMBER (NON-DIMENSIONAL POTENTIAL ENERGY)
! LD = DISCRETE EDDY PHYSICAL LENGTH (METERS)
! LS = DISCRETE LENGTH OF EDDY SUB INTERVAL ~ L/3
! ENERGY = ARGUMENT OF SQUARE-ROOT TERM IN RATE DENSITY
! A NEGATIVE ENERGY INDICATES THE EDDY IS NOT ENERGETIC
! ENOUGH TO TURN OVER AND IS HENCE NOT IMPLEMENTED

M = KE/3
MM = M/3

! CREATE SUB-INTERVAL TRIPLET MAP FOR ARBITRARY NUMBER
! OF POINTS; NOTE: THIS MAPPING IS NOT CONSERVATIVE,
! BUT THIS IS NOT THE MAPPING USED IN THE EDDY IMPLEMENTATION,
! AND SINCE THE INVERVAL, 'M', IS NOT GENERALLY DIVISIBLE
! BY 3, THE STANDARD DISCRETE MAPPING DOES NOT WORK.  THE
! WAY TO THINK ABOUT THIS MAPPING IS TO WRITE OUT THE
! NUMBERS IN ASCENDING ORDER, THEN DESCENDING ORDER, THEN
! ASCENDING ORDER AGAIN, WITHOUT WRITING A NUMBER TWICE IN
! A ROW.  THEN TAKE EVERY THIRD NUMBER AS THE MAPPING. AS
! AN EXAMPLE, SUPPOSE KE = 21, THEN M = 7, SO WE WRITE
!
!      1  2  3  4  5  6  7  6  5  4  3  2  1  2  3  4  5  6  7
!      ~        ~        ~        ~        ~        ~        ~
!
! HENCE, FOR        J  =  1  2  3  4  5  6  7
!                 F(J) =  1  4  7  4  1  4  7
!                 K(J) =  0 -2 -4  0  4  2  0
!                JJ(J) =  0  2  4  0  4  2  0
!
! NOTICE THAT SOME OF THE INDICES HAVE BEEN LOST AND SOME
! HAVE BEEN REPEATED, SO THIS   IS NOT CONSERVATIVE, BUT WILL
! DO FOR SUB-INTERVAL EVALUATION.

DO I = 1,M
   FFF(I) = I
ENDDO
I = M
DO J = M+1,2*M-1
   I = I-1
   FFF(J) = I
ENDDO
I = 1
DO J = 2*M,KE-2
   I = I+1
   FFF(J) = I
ENDDO
I = 0
DO J = 1,KE-2,3
   I = I+1
   F(I) = FFF(J)
   K(I) = I-F(I)
   JJ(I) = ABS(K(I))
ENDDO

! EVALUATE SUB-INTERVAL I
SUBINTERVAL_LOOP: DO I = 1,3

   DO J = 1,M
      J_SUB = (I-1)*M+J
      U_SUB(J) = U(J_SUB)
      V_SUB(J) = V(J_SUB)
      W_SUB(J) = W(J_SUB)
      RHO_SUB(J) = RHO(J_SUB)
   ENDDO

   ! INTEGRAL EDDY PARAMETERS
   U_RHOK = 0._EB
   V_RHOK = 0._EB
   W_RHOK = 0._EB
   U_RHOJ = 0._EB
   V_RHOJ = 0._EB
   W_RHOJ = 0._EB
   RHO_K  = 0._EB
   RHO_J  = 0._EB
   RHO_KK = 0._EB
   RHO_JK = 0._EB

   ! NOTE: JJ(J) IS SUBSTITUTED FOR J(J) BECAUSE WE ARE
   ! RUNNING OUT OF SYMBOLS.
   DO J = 1,MM
      U_RHOK = U_RHOK + RHO(F(J))*U_SUB(F(J))*K(J)
      V_RHOK = V_RHOK + RHO(F(J))*V_SUB(F(J))*K(J)
      W_RHOK = W_RHOK + RHO(F(J))*W_SUB(F(J))*K(J)
      U_RHOJ = U_RHOJ + RHO(F(J))*U_SUB(F(J))*JJ(J)
      V_RHOJ = V_RHOJ + RHO(F(J))*V_SUB(F(J))*JJ(J)
      W_RHOJ = W_RHOJ + RHO(F(J))*W_SUB(F(J))*JJ(J)
      RHO_K  = RHO_K  + RHO_SUB(F(J))*K(J)
      RHO_J  = RHO_J  + RHO_SUB(F(J))*JJ(J)
      RHO_KK = RHO_KK + RHO_SUB(F(J))*K(J)*K(J)
      RHO_JK = RHO_JK + RHO_SUB(F(J))*JJ(J)*K(J)
   ENDDO
   DO J = MM+1,2*MM
      U_RHOK = U_RHOK + RHO(F(J))*U_SUB(F(J))*K(J)
      V_RHOK = V_RHOK + RHO(F(J))*V_SUB(F(J))*K(J)
      W_RHOK = W_RHOK + RHO(F(J))*W_SUB(F(J))*K(J)
      U_RHOJ = U_RHOJ + RHO(F(J))*U_SUB(F(J))*JJ(J)
      V_RHOJ = V_RHOJ + RHO(F(J))*V_SUB(F(J))*JJ(J)
      W_RHOJ = W_RHOJ + RHO(F(J))*W_SUB(F(J))*JJ(J)
      RHO_K  = RHO_K  + RHO_SUB(F(J))*K(J)
      RHO_J  = RHO_J  + RHO_SUB(F(J))*JJ(J)
      RHO_KK = RHO_KK + RHO_SUB(F(J))*K(J)*K(J)
      RHO_JK = RHO_JK + RHO_SUB(F(J))*JJ(J)*K(J)
   ENDDO
   DO J = 2*MM+1,M
      U_RHOK = U_RHOK + RHO(F(J))*U_SUB(F(J))*K(J)
      V_RHOK = V_RHOK + RHO(F(J))*V_SUB(F(J))*K(J)
      W_RHOK = W_RHOK + RHO(F(J))*W_SUB(F(J))*K(J)
      U_RHOJ = U_RHOJ + RHO(F(J))*U_SUB(F(J))*JJ(J)
      V_RHOJ = V_RHOJ + RHO(F(J))*V_SUB(F(J))*JJ(J)
      W_RHOJ = W_RHOJ + RHO(F(J))*W_SUB(F(J))*JJ(J)
      RHO_K  = RHO_K  + RHO_SUB(F(J))*K(J)
      RHO_J  = RHO_J  + RHO_SUB(F(J))*JJ(J)
      RHO_KK = RHO_KK + RHO_SUB(F(J))*K(J)*K(J)
      RHO_JK = RHO_JK + RHO_SUB(F(J))*JJ(J)*K(J)
   ENDDO
   LL = (REAL(M,EB)-1._EB)**2
   LLL = (REAL(M,EB)-1._EB)**3
   U_RHOK = U_RHOK/LL
   V_RHOK = V_RHOK/LL
   W_RHOK = W_RHOK/LL
   U_RHOJ = U_RHOJ/LL
   V_RHOJ = V_RHOJ/LL
   W_RHOJ = W_RHOJ/LL
   RHO_K  = RHO_K/LL
   RHO_J  = RHO_J/LL
   RHO_KK = RHO_KK/LLL
   RHO_JK = RHO_JK/LLL

   ! AVERAGE DYNAMIC VISCOSITY
   MU0 = 0._EB
   DO J = (I-1)*M+1,I*M
      MU0 = MU0+MU(I)
   ENDDO
   MU0 = MU0/M

   ! EDDY REYNOLDS NUMBERS
   ! NOTE THAT MY DISCRETE EDDY LENGTH IS DIFFERENT FROM VEBJORN'S
   LD = (REAL(KE,EB)-1._EB)*H
   LS = (M-1)*H
   IF (MU0<EPS_LOC) THEN
      MU0 = 1._EB
      Z = 0._EB
   ENDIF
   RE_U = U_RHOK*LD/MU0
   RE_V = V_RHOK*LD/MU0
   RE_W = W_RHOK*LD/MU0

   ! ASHURST VARIABLE DENSITY PARAMETERS
   A = RHO_K/RHO_J
   P_U = (U_RHOK - A*U_RHOJ)*LD**2
   P_V = (V_RHOK - A*V_RHOJ)*LD**2
   P_W = (W_RHOK - A*W_RHOJ)*LD**2
   S = (0.5_EB*(A**2+1._EB)*RHO_KK - A*RHO_JK)*LD**3

   ! SUB-EDDY GRASHOF NUMBER
   ! NOTE THAT THE LENGTH SCALE THAT MULTIPLIES THE VISCOUS TERM
   ! IS LEFT AS LD, THE TOTAL EDDY LENGTH; LS IS THE SUB-EDDY
   ! INTERVAL AND SO RHO_K*LS IS THE MASS IN THE SUB-INTERVAL
   IF (G>EPS_LOC) THEN
      GR = (2._EB*LD**2*RHO_KK/MU0**2)*G*RHO_K*LS
   ELSE
      GR = 0._EB
   ENDIF

   ! HARMANJEET SHIHN (SUNY, BUFFALO) BAROCLINIC TORQUE TERM
   ! SEE COMMENTS IN 'EVALUATEEDDY'
   IF (GC>EPS_LOC) THEN
      RHO_LEFT = 0._EB
      RHO_RIGHT = 0._EB
      DO J = 1,M/2
         J_SUB = (I-1)*M+J
         RHO_LEFT = RHO_LEFT + RHO_SUB(J_SUB)
      ENDDO
      RHO_LEFT = 2._EB*RHO_LEFT/M
      DO J = M/2+1,M
         J_SUB = (I-1)*M+J
         RHO_RIGHT = RHO_RIGHT + RHO_SUB(J_SUB)
      ENDDO
      RHO_RIGHT = 2._EB*RHO_RIGHT/M
      RHO_DIFF = ABS(RHO_RIGHT-RHO_LEFT)
      RI_RE2 = 2._EB*LD**2*RHO_KK/MU0**2*GC*RHO_DIFF*LS
   ELSE
      RI_RE2 = 0._EB
   ENDIF

   ! ARGUMENT OF SQRT TERM...
   ! ENERGY = RE_U**2+RE_V**2+RE_W**2 + GR + B*Ri*Re**2 - Z**2
   ENERGY = RHO_KK*LD/(2*S*MU0**2)*(P_U**2 + P_V**2 + P_W**2) + GR + B*RI_RE2 - Z**2
   IF (ENERGY<TWO_EPSILON_EB) THEN
      LEDDY = 0
      PA = 0._EB
      RETURN
   ENDIF

ENDDO SUBINTERVAL_LOOP

END SUBROUTINE SUBINTERVALEVALUATION


SUBROUTINE KERNELAMP(C_U,C_V,C_W, &
                     RHO_K, &
                     FLAVOR,ALPHA,KE,H,G, &
                     P_U,P_V,P_W,S)
IMPLICIT NONE

! INTENT OUT:
REAL(EB) :: C_U,C_V,C_W

! INTENT IN:
REAL(EB) :: RHO_K
INTEGER :: FLAVOR,KE
REAL(EB) :: ALPHA,H,G
REAL(EB) :: P_U,P_V,P_W,S

! OUTPUTS:
! C_U = U VELOCITY KERNEL AMPLITUDE (M/S)
! C_V = V VELOCITY KERNEL AMPLITUDE (M/S)
! C_W = W VELOCITY KERNEL AMPLITUDE (M/S)
! INPUTS:
! U_RHOK = MEASURE OF U MASS FLUX (KG/M^S*S)
! V_RHOK = MEASURE OF V MASS FLUX (KG/M^2*S)
! W_RHOK = MEASURE OF W MASS FLUX (KG/M^2*S)
! RHO_K = MEASURE OF MASS DISPLACEMENT (KG/M^3)
! RHO_KK = STD DEV OF MASS DISPLACEMENT (KG/M^3)
! FLAVOR = INTEGER: 0 = SPECIFY ALPHA, 1 = USE ALPHA = 2/3 SPECIAL CASE
!        WHICH IS USED MOST OFTEN AND WE CAN MAKE CODING MORE EFFICIENT
! ALPHA = FRACTION OF MAX ENERGY EXCHANGE IN PRESSURE-SCRAMBLING MODEL
! KE = NUMBER OF PODINTS IN DISCRETE EDDY
! H = ODT GRID SPACING
! G = GRAVITATIONAL ACCELERATION (E.G. -9.8 M/S^2)
! RNTAG = RANDOM NUMBER GENERATOR: 0 = RAND, 1 = RNG
!
! LOCAL VARIABLES:
REAL(EB) :: ARG,SGN,COEF_C,LD
REAL(EB) :: PU2,PV2,PW2,SUMP2
REAL(EB) :: BETA_U,BETA_V,BETA_W
REAL(EB) :: GRAVITY_TERM
REAL(EB) :: RN
REAL(EB), PARAMETER :: EPS_LOC=1.E-10_EB

! ARG = USED TO SPEED EXECUTION
! SGN = HOLDS +1. OR -1.
! COEF_C = USED TO SPEED EXECUTION
! LD = DISCRETE PHYSICAL EDDY LENGTH
! RN = RANDOM NUMBER ON DINTERVAL [0,1]
! UIK2 = EDDY KINETIC ENERGY AFTER PRESSURE SCRAMBLING
! BETA_I = COEFFICIENT OF GRAVITY TERM, MUST SUM TO 1
! GRAVITY_TERM = USED TO SPEED CODING IF G = 0
!
! ------------------------------------------------------------------------
! NOTE: THIS COEFFICIENT DOES NOT CONTAIN THE LENGTH PARAMETER, H, AND
! SO WHEN WE ADD C_I*K(Y) TO THE VELOCITY U_I WE DO NOT USE H EITHER.
! COEF_C = 1./((KE-1)*RHO_KK)

! NOTE: IN THE VARIABLE DENSITY FORMULATION
!       S = ((KE-1)*H)^3*(0.5*(A^2+1)*RHO_KK - A*RHO_JK)
!       THE "H" IN THE NUMERATOR OF THE COEFFICIENT
!       IN NEEDED SO THAT WE DO NOT NEED TO MULTIPLY
!       BY THE KERNEL BY H IN THE STIRRING EVENT.
COEF_C = H/(2._EB*S)

LD = (REAL(KE,EB)-1._EB)*H
GRAVITY_TERM = 0._EB

IF (FLAVOR==1) THEN

   IF (ABS(G) > EPS_LOC) THEN
      GRAVITY_TERM = 4._EB*S*LD**2*G*RHO_K
   ENDIF
   ARG = (P_U**2 + P_V**2 + P_W**2 + GRAVITY_TERM)/3._EB

   IF (ABS(P_U) < EPS_LOC) THEN
      CALL RANDOM_NUMBER(RN)
      IF (RN > 0.5_EB) THEN
         SGN = 1._EB
      ELSE
         SGN = -1._EB
      ENDIF
   ELSE
      SGN = SIGN(1._EB,P_U)
   ENDIF
   C_U = COEF_C*( -P_U + SGN*SQRT(ARG) )

   IF (ABS(P_V) < EPS_LOC) THEN
      CALL RANDOM_NUMBER(RN)
      IF (RN > 0.5_EB) THEN
         SGN = 1._EB
      ELSE
         SGN = -1._EB
      ENDIF
   ELSE
      SGN = SIGN(1._EB,P_V)
   ENDIF
   C_V = COEF_C*( -P_V + SGN*SQRT(ARG) )

   IF (ABS(P_W) < EPS_LOC) THEN
      CALL RANDOM_NUMBER(RN)
      IF (RN > 0.5_EB) THEN
         SGN = 1._EB
      ELSE
         SGN = -1._EB
      ENDIF
   ELSE
      SGN = SIGN(1._EB,P_W)
   ENDIF
   C_W = COEF_C*( -P_W + SGN*SQRT(ARG) )

ELSE

   PU2=P_U**2 + ALPHA*(-P_U**2 + 0.5*P_V**2 + 0.5*P_W**2)
   PV2=P_V**2 + ALPHA*( 0.5*P_U**2 - P_V**2 + 0.5*P_W**2)
   PW2=P_W**2 + ALPHA*( 0.5*P_U**2 + 0.5*P_V**2 - P_W**2)
   SUMP2 = PU2 + PV2 + PW2
   BETA_U = PU2/SUMP2
   BETA_V = PV2/SUMP2
   BETA_W = PW2/SUMP2
   ! SUM_BETA = BETA_U + BETA_V + BETA_W; ! CHECK, MUST EQUAL 1

   IF (ABS(P_U) < EPS_LOC) THEN
      CALL RANDOM_NUMBER(RN)
      IF (RN > 0.5_EB) THEN
         SGN = 1._EB
      ELSE
         SGN = -1._EB
      ENDIF
   ELSE
      SGN = SIGN(1._EB,P_U)
   ENDIF
   IF (ABS(G) > EPS_LOC) THEN
      GRAVITY_TERM = 4._EB*BETA_U*S*LD**2*G*RHO_K
   ENDIF
   C_U = COEF_C*( -P_U + SGN*SQRT( PU2 + GRAVITY_TERM ) )

   IF (ABS(P_V) < EPS_LOC) THEN
      CALL RANDOM_NUMBER(RN)
      IF (RN > 0.5) THEN
         SGN = 1._EB
      ELSE
         SGN = -1._EB
      ENDIF
   ELSE
      SGN = SIGN(1._EB,P_V)
   ENDIF
   IF (ABS(G) > EPS_LOC) THEN
      GRAVITY_TERM = 4._EB*BETA_V*S*LD**2*G*RHO_K
   ENDIF
   C_V = COEF_C*( -P_V + SGN*SQRT( PV2 + GRAVITY_TERM ) )

   IF (ABS(P_W) < EPS_LOC) THEN
      CALL RANDOM_NUMBER(RN)
      IF (RN > 0.5) THEN
         SGN = 1._EB
      ELSE
         SGN = -1._EB
      ENDIF
   ELSE
      SGN = SIGN(1._EB,P_W)
   ENDIF
   IF (ABS(G) > EPS_LOC) THEN
      GRAVITY_TERM = 4._EB*BETA_W*S*LD**2*G*RHO_K
   ENDIF
   C_W = COEF_C*( -P_W + SGN*SQRT( PW2 + GRAVITY_TERM ) )

ENDIF

END SUBROUTINE KERNELAMP


SUBROUTINE SHEARLAYERODT
USE COMP_FUNCTIONS, ONLY: GET_FILE_NUMBER
IMPLICIT NONE

! THIS ROUTINE IS INTENDED AS A VERIFICATION TEST FOR THE ODT CODE.
! THE TEMPORALLY DEVELOPING SHEAR LAYER PROBLEM WAS RUN BY [1] WITH
! C = 3.68 AND Z = 0.14. [1] THEN COMPARES THE EVOLUTION OF THE
! MOMENTUM DEFICIT, \DELTA_M, WITH THE DNS RESULTS OF MOSER AND
! ROGERS (JFM, 1998), WHO OBTAIN D(\DELTA_M)/DT = 0.014*(\DELTA U)
! FOR THE SELF-SIMILAR REGION.  IN [1], THE ODT FIELD IS INITIALIZED
! WITH A PIECEWISE-LINEAR PROFILE SUCH THAT THE RE BASED ON \DELTA_M
! WAS 427.  THAT IS, RE_D = (\DELTA U)*(\DELTA_M)/NU = 427.
!
!                            Y ^
!                              :---->| U(Y) = (\DELTA U)/2
!                              :     |
!                              :     |
!                              :     |
!                              :     |
!                              :     | _ _ _ _ _ _ _
!                              :    /           ^
!   U(Y) = ((\DELTA U)/DY)*Y   :  /             |
!                              :/____> Y = 0    |
!                             /:                DY
!                           /  :                |
!                         /    :       _ _ _ _ _V_ _
!                        |     :
!                        |     :
!                        |     :
!                        |     :
!                        |     :
!   U(Y) = -(\DELTA U)/2 |<----:
!                              :
!
!
!
! HERE WE PICK (\DELTA U) = 2, AND DY = 2.  THE MOMENTUM DEFICIT FOR
! THE SHEAR LAYER IS GIVEN BY
!
! \DELTA_M \EQUIV \INT_{-\INFTY}^{+\INFTY} ( 1/4 - (U(Y)/(\DELTA U))^2 ) DY
!
! HENCE, WE OBTAIN FOR RE_D = 427, NU = 2/(3*427).

INTEGER, PARAMETER :: N_PTS = 100000
INTEGER, PARAMETER :: N_SCALAR = 1

INTEGER :: J,J1,JD,BC,JD1,JD2,K
REAL(EB) :: NU,SC,RE,DU,DY,D,H,H2,L_MAX,UL,UH,ETA
REAL(EB) :: R0,RD
REAL(EB) :: R(N_PTS),RHO(N_PTS),PHI(N_SCALAR,N_PTS)
REAL(EB) :: U(N_PTS),V(N_PTS),W(N_PTS),MU(N_PTS)
REAL(EB) :: U_N(N_PTS),V_N(N_PTS),W_N(N_PTS),PHI_N(N_PTS)
REAL(EB) :: MW0,MW(N_SCALAR+1),PRT

REAL(EB) :: DT,DTS,T,TNP1,TOTAL_TIME,DM
REAL(EB) :: PA,PA_SUM,PA_MEAN,PA_MAX,PA_MAXIMP,PA_CUM
INTEGER :: N_TRIALS,N_EDDIES,N_LMAX,N_SAMPLES
INTEGER :: LEDDY
REAL(EB) :: Y0,L
INTEGER :: J0,KE
! TRIDIAGONAL SOLVE
REAL(EB) :: A1(N_PTS),A0(N_PTS)
REAL(EB) :: B_U(N_PTS),B_V(N_PTS),B_W(N_PTS)
REAL(EB) :: UBC0,VBC0,WBC0,UBCN,VBCN,WBCN
! PROFILE OUTPUT
CHARACTER(6) :: EXT
LOGICAL :: LIMPLICIT

INTEGER :: LU_MOM,LU_EDD,LU_PRO
CHARACTER(80) :: FN_MOM,FN_EDD

! PROBLEM PARAMETERS
! BC = 2
! RE = 427._EB
! DU = 2._EB
! DY = 2._EB
! NU = DU*DY/RE
! D = 3000._EB*(20._EB*NU)/DU
! J1 = 1
! JD = N_PTS
! H = D/(JD-J1+1)
! L_MAX = D
! N_LMAX = INT(D/H)
! LIMPLICIT = .TRUE.

! ASHURST PARAMETERS
! THE SHEAR LAYER SIMULATION IN THE 2001 JFM PAPER WAS RUN BY BILL ASHURST
! AND HE HAS INFORMED ME THAT THESE ARE THE PARAMETERS HE USED IN RUNNING
! THE SIMULATION.  THE NUMBER OF ODT POINTS WAS 16383, VISCOSITY = 0.05,
! GRID SPACING = UNITY.  ALSO, THE DIFFUSION EQUATION WAS SOLVED WITH
! FORWARD EULER (EXPLICIT) WITH CFL = 0.2.  THE REST OF THE DETAILS ARE
! GIVEN IN THE 2001 JFM PAPER, BUT NOTE THE RATE CONSTANT SHOULD BE 3.68.
! BC = 2
! RE = 427._EB
! DU = 2._EB
! DY = 64.05_EB
! NU = 0.05_EB
! SC = 1._EB
! H = 1._EB
! D = N_PTS*H
! L_MAX = D
! N_LMAX = INT(D/H)
! J1 = 1
! JD = N_PTS
! LIMPLICIT = .FALSE.

! ROGERS AND MOSER (1994)
BC = 2
UL = -0.5_EB
UH =  0.5_EB
DU = UH-UL ! 1
RE = 800._EB
DY = 100._EB
NU = RE/(DU*DY)
SC = 1._EB
ETA = DY*RE**(-0.75_EB) ! KOLMOGOROV LENGTH SCALE
H = ETA/6._EB
D = N_PTS*H
L_MAX = D
N_LMAX = INT(D/H)
J1 = 1
JD = N_PTS
LIMPLICIT = .TRUE.

!WRITE(LU_ERR,*) 'KOLMOGOROV = ', ETA

H2 = H**2

! INITIALIZE DOMAIN
R(J1) = -D/2._EB + H/2._EB
DO J = J1+1,JD
   R(J) = R(J-1) + H
ENDDO
R0 = R(J1)-H/2._EB
RD = R(JD)+H/2._EB

! INITIALIZE VELOCITY FIELD (SEE SKETCH ABOVE)
JD1 = JD/2-INT((DY/2)/H)
JD2 = JD/2+INT((DY/2)/H)
DO J = J1,JD1
   U(J) = UL
ENDDO
DO J = JD1+1,JD2
   U(J) = DU/DY*R(J)
ENDDO
DO J = JD2+1,JD
   U(J) = UH
ENDDO

DO J = J1,JD
   V(J) = 0._EB
   W(J) = 0._EB
ENDDO
! LET PHI(1) BE MASS FRACTION OF AIR
! BOTTOM HALF OF DOMAIN IS INITIALLY AIR
DO J = J1,(JD-J1+1)/2
   PHI(1,J) = 1._EB
ENDDO
! TOP HALF OF DOMAIN IS INITIALLY HE (STABLE STRATIFICATION)
DO J = (JD-J1+1)/2+1,JD
   PHI(1,J) = 0._EB
ENDDO

! COMPUTE INITIAL DENSITY
! FOR FUN, LET COMP 1 BE AIR, COMP 2 BE HELIUM
MW(1) = 24.4530_EB
MW(2) = 24.4530_EB
! P = 1. ATM
! R = 0.082057 KG*ATM/(KMOL*K)
! T = 298. K
PRT = 1._EB/(0.082057_EB*298._EB)
DO J = J1,JD
   MW0 = 1._EB/( PHI(1,J)/MW(1) + (1._EB-PHI(1,J))/MW(2) )
   RHO(J) = MW0*PRT
   MU(J) = NU*RHO(J)
ENDDO

! SIMULATION TIME
TOTAL_TIME = 200._EB*DY/DU
T = 0._EB
K = 0

! DIFFUSION TIME SCALE
DT = 0.2*H2/NU                            ! EXPLICIT SOLVER
IF (LIMPLICIT) DT = TOTAL_TIME/1000._EB    ! IMPLICIT SOLVER
! INITIAL SAMPLING TIME PERIOD (ADJUSTED LATER)
DTS = 0.01_EB*DT

IF (LIMPLICIT) THEN
   ! SET COEFFICIENTS
   A1(1) = 0.5_EB
   A0(1) = -(1.5+H2/(4*DT*NU))
   DO J = 2,N_PTS-1
      A1(J) = 1._EB
      A0(J) = -(2._EB+H2/(DT*NU))
   ENDDO
   A1(N_PTS) = 0.5_EB
   A0(N_PTS) = -(1.5_EB+H2/(4._EB*DT*NU))
   ! SET DIRICHLET BCS
   UBC0 = U(1)
   VBC0 = V(1)
   WBC0 = W(1)
   UBCN = U(N_PTS)
   VBCN = V(N_PTS)
   WBCN = W(N_PTS)
ENDIF

PA_SUM = 0._EB
PA_MEAN = 0._EB
PA_MAX = 1._EB
PA_MAXIMP = 0._EB
PA_CUM = 0._EB
N_TRIALS = 0
N_EDDIES = 0
N_SAMPLES = 0

LU_MOM = GET_FILE_NUMBER()
LU_EDD = GET_FILE_NUMBER()
FN_MOM = TRIM(CHID)//'_mom.dat'
FN_EDD = TRIM(CHID)//'_edd.dat'
OPEN(LU_MOM,FILE=FN_MOM,STATUS='UNKNOWN',FORM='FORMATTED')
OPEN(LU_EDD,FILE=FN_EDD,STATUS='UNKNOWN',FORM='FORMATTED')

TIME_LOOP: DO

   TNP1 = T + DT
   K = K + 1
   IF (TNP1 > TOTAL_TIME) THEN
      EXIT TIME_LOOP
   ENDIF

   EDDY_SAMPLING_LOOP: DO

      LEDDY = 0
      CALL GODT(U,V,W,RHO,PHI,MU,DTS,PA,LEDDY,Y0,L, &
                J0,KE,N_PTS,N_SCALAR,H,N_LMAX,0._EB,0._EB, &
                BC,1,1,0,0._EB,3.68_EB,0._EB,PA_MAX,1.E-6_EB)

          !GODT(U,V,W,RHO,PHI,MU,DTS,PA,LEDDY,Y0,L, &
          !     J0,KE,N_ODT,N_SCALAR,H,N_LMAX,G,GC, &
          !     BC,DEBUG,LE_SUP,PDFTAG,RET,C_EDDY,Z_EDDY,PAMAX,DTS_MIN)

      N_TRIALS = N_TRIALS + 1
      N_SAMPLES = N_SAMPLES + 1
      PA_SUM = PA_SUM + PA
      PA_CUM = PA_CUM + PA
      IF (PA > PA_MAXIMP) PA_MAXIMP = PA

      IF (N_TRIALS>10000 .AND. PA_SUM>TWO_EPSILON_EB) THEN
         DTS = 0.1_EB/(PA_SUM/REAL(N_TRIALS,EB))*DTS
         N_TRIALS = 0
         PA_SUM = 0._EB
      ENDIF

      IF (LEDDY==1) THEN
         N_EDDIES = N_EDDIES + 1
         WRITE(LU_EDD,*) T,Y0-D/2,L
      ENDIF

      T = T + DTS
      IF (T > TNP1) THEN
         T = TNP1
         EXIT EDDY_SAMPLING_LOOP
      ENDIF

   ENDDO EDDY_SAMPLING_LOOP

   ! INTEGRATE DIFFUSION EQUATION
   IMPLICIT_IF: IF (LIMPLICIT) THEN
      ! IMPLICIT BACKWARD EULER, ONLY WITH DIRICHLET B.C.S
      IF (BC/=2) CALL SHUTDOWN('ERROR: must use Dirichlet bc with LIMPLICIT ODT')
      ! SET RIGHT-HAND-SIDE
      B_U(1) = -(U(1)*H2/(4._EB*DT*NU) + UBC0)
      B_V(1) = -(V(1)*H2/(4._EB*DT*NU) + VBC0)
      B_W(1) = -(W(1)*H2/(4._EB*DT*NU) + WBC0)
      DO J = 2,N_PTS-1
         B_U(J) = -U(J)*H2/(DT*NU)
         B_V(J) = -V(J)*H2/(DT*NU)
         B_W(J) = -W(J)*H2/(DT*NU)
      ENDDO
      B_U(N_PTS) = -(U(N_PTS)*H2/(4._EB*DT*NU) + UBCN)
      B_V(N_PTS) = -(V(N_PTS)*H2/(4._EB*DT*NU) + VBCN)
      B_W(N_PTS) = -(W(N_PTS)*H2/(4._EB*DT*NU) + WBCN)
      CALL TRIDIAGODT(A1,A0,A1,B_U,U,N_PTS)
      CALL TRIDIAGODT(A1,A0,A1,B_V,V,N_PTS)
      CALL TRIDIAGODT(A1,A0,A1,B_W,W,N_PTS)
   ELSE IMPLICIT_IF
      ! EXPLICIT FORWARD EULER
      ! NOTICE THAT THE FOLLOWING ASSUMES YOUR DOMAIN IS LARGE ENOUGH
      ! THAT THE PROFILE IS FLAT AT THE DOMAIN EXTENTS.
      DO J = J1,JD
         U_N(J) = U(J)
         V_N(J) = V(J)
         W_N(J) = W(J)
         PHI_N(J) = PHI(1,J)
      ENDDO
      DO J = J1+1,JD-1
         U(J)=U_N(J) + DT*NU*(U_N(J+1)-2._EB*U_N(J)+U_N(J-1))/H2
         V(J)=V_N(J) + DT*NU*(V_N(J+1)-2._EB*V_N(J)+V_N(J-1))/H2
         W(J)=W_N(J) + DT*NU*(W_N(J+1)-2._EB*W_N(J)+W_N(J-1))/H2
         PHI(1,J) = PHI_N(J) + DT*NU/SC*(PHI_N(J+1)-2._EB*PHI_N(J)+PHI_N(J-1))/H2
      ENDDO
      DO J = J1,JD
         MW0 = 1._EB/( PHI(1,J)/MW(1) + (1._EB-PHI(1,J))/MW(2) )
         RHO(J) = MW0*PRT
      ENDDO
   ENDIF IMPLICIT_IF

   IF (MOD(K,1)==0) THEN

      ! COMPUTE MOMENTUM DEFICIT
      DM = 0._EB
      DO J = J1,JD
         DM = DM + (0.25_EB - (U(J)/DU)**2)*H
      ENDDO
      WRITE(LU_MOM,*) T,DM

      IF (N_TRIALS>0) THEN
         PA_MEAN = PA_SUM/N_TRIALS
      ENDIF

      WRITE(LU_ERR,201) T,DTS,PA_MEAN,PA_MAXIMP,N_EDDIES

   ENDIF
   IF (.TRUE. .AND. MOD(K,1)==0) THEN
      ! WRITE VELOCITY PROFILES
      LU_PRO=GET_FILE_NUMBER()
      WRITE(EXT,202) K
      OPEN(LU_PRO,FILE=TRIM(CHID)//'_profiles_'//TRIM(EXT)//'.dat',STATUS='UNKNOWN',FORM='FORMATTED')
      DO J = J1,JD
         WRITE(LU_PRO,203) R(J),U(J),V(J),RHO(J)
      ENDDO
      CLOSE(LU_PRO)
   ENDIF

ENDDO TIME_LOOP
CLOSE(LU_MOM)
CLOSE(LU_EDD)

WRITE(LU_ERR,*) 'SHEAR_LAYER_ODT FINISHED!!'
WRITE(LU_ERR,*) 'PA_CUM = ',PA_CUM/N_SAMPLES
STOP

201 FORMAT (F12.6,F12.6,F12.6,F12.6,I10)
202 FORMAT (I5.5)
203 FORMAT (E14.6,F12.6,F12.6,F12.6)

END SUBROUTINE SHEARLAYERODT


SUBROUTINE TRIDIAGODT(A,B,C,R,U,N)
IMPLICIT NONE
INTEGER, INTENT(IN) :: N
REAL(EB) :: A(N),B(N),C(N),R(N),U(N)
INTEGER, PARAMETER :: NMAX = 1000000
INTEGER J
REAL(EB) :: BB,AA(NMAX)

BB = B(1)
U(1) = R(1)/BB
! forward sub
DO J = 2,N
   AA(J) = C(J-1)/BB
   BB = B(J)-A(J)*AA(J)
   IF (ABS(BB)<TWO_EPSILON_EB) THEN
      CALL SHUTDOWN('ERROR: tridiagodt failed') ; RETURN
   ENDIF
   U(J) = (R(J)-A(J)*U(J-1))/BB
ENDDO
! backward sub
DO J = N-1,1,-1
   U(J) = U(J)-AA(J+1)*U(J+1)
ENDDO

END SUBROUTINE TRIDIAGODT


END MODULE ONE_DIMENSIONAL_TURBULENCE
