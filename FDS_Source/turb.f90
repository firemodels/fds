! This module is useful for verification tests and development of
! turbulence models.

MODULE TURBULENCE

USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_POINTERS
USE MESH_VARIABLES
USE COMP_FUNCTIONS

IMPLICIT NONE

CHARACTER(255), PARAMETER :: turbid='$Id$'
CHARACTER(255), PARAMETER :: turbrev='$Revision$'
CHARACTER(255), PARAMETER :: turbdate='$Date$'

PRIVATE
PUBLIC :: INIT_TURB_ARRAYS, VARDEN_DYNSMAG, WANNIER_FLOW, &
          GET_REV_turb, WALL_MODEL, COMPRESSION_WAVE, VELTAN2D,VELTAN3D, &
          SYNTHETIC_TURBULENCE, SYNTHETIC_EDDY_SETUP, TEST_FILTER, EX2G3D, TENSOR_DIFFUSIVITY_MODEL, &
          TWOD_VORTEX_CERFACS, HEAT_FLUX_MODEL, ABL_HEAT_FLUX_MODEL, RNG_EDDY_VISCOSITY, &
          NS_ANALYTICAL_SOLUTION, NS_U_EXACT, NS_V_EXACT, NS_H_EXACT, SANDIA_DAT, SPECTRAL_OUTPUT, SANDIA_OUT
 
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


SUBROUTINE WALL_MODEL(SLIP_FACTOR,U_TAU,Y_PLUS,U,NU,DY,S)

REAL(EB), INTENT(OUT) :: SLIP_FACTOR,U_TAU,Y_PLUS
REAL(EB), INTENT(IN) :: U,NU,DY,S ! S is the roughness length scale (Pope's notation)

REAL(EB), PARAMETER :: RKAPPA=1._EB/0.41_EB ! 1/von Karman constant
REAL(EB), PARAMETER :: B=5.2_EB,BTILDE_ROUGH=8.5_EB,BTILDE_MAX=9.5_EB ! see Pope (2000) pp. 294,297,298
REAL(EB), PARAMETER :: S0=1._EB,S1=5.83_EB,S2=30._EB ! approx piece-wise function for Fig. 7.24, Pope (2000) p. 297
REAL(EB), PARAMETER :: Y1=5._EB,Y2=30._EB
REAL(EB), PARAMETER :: U1=5._EB,U2=RKAPPA*LOG(Y2)+B
REAL(EB), PARAMETER :: EPS=1.E-10_EB

REAL(EB) :: Y_CELL_CENTER,TAU_W,BTILDE,DELTA_NU,S_PLUS,DUDY
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

LES_IF: IF (LES) THEN

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
            U_TAU = SQRT(ABS(TAU_W))
            DUDY = ABS(U)/Y_CELL_CENTER
         !ELSE IF (Y_PLUS < Y2) THEN
         !   ! buffer layer
         !   TAU_W = ( U/U_PLUS_BUFFER_SEMILOG(Y_PLUS) )**2
         !   U_TAU = SQRT(TAU_W)
         !   DUDY = 0.5_EB*(ABS(U)/Y_CELL_CENTER + U_TAU*RKAPPA/Y_CELL_CENTER)
         ELSE
            ! log layer
            TAU_W = ( U/(RKAPPA*LOG(Y_PLUS)+B) )**2
            U_TAU = SQRT(ABS(TAU_W))
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
         U_TAU = SQRT(ABS(TAU_W))
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


SUBROUTINE HEAT_FLUX_MODEL(H,YPLUS,U_TAU,K,RHO,CP,MU)

! Kiyoung Moon, Yonsei University
! Ezgi Oztekin, Technology and Manangement International

REAL(EB), INTENT(OUT) :: H
REAL(EB), INTENT(IN) :: YPLUS,U_TAU,K,RHO,CP,MU
REAL(EB) :: PR_M,TPLUS,B_T !,T1,T2,CA,CB
REAL(EB), PARAMETER :: RKAPPA=1._EB/0.41_EB
REAL(EB), PARAMETER :: Y1=5._EB,Y2=30._EB
REAL(EB), PARAMETER :: LOG_Y1=LOG(Y1),LOG_Y2=LOG(Y2),DLOGY=LOG(Y2/Y1)

PR_M = CP*MU/K

IF (YPLUS < Y_WERNER_WENGLE) THEN
   ! viscous sublayer
   TPLUS = PR_M*YPLUS
ELSE
   B_T = (3.85_EB*PR_M**ONTH-1.3_EB)**2 + 2.12_EB*LOG(PR_M) ! Kader, 1981
   !IF (YPLUS < Y2) THEN
   !   ! buffer layer
   !   T2 = PR*RKAPPA*LOG_Y2 + B_T
   !   T1 = PR_M*Y1
   !   T2 = MAX(T1,T2)
   !   CA = (T2-T1)/DLOGY
   !   CB = T1-CA*LOG_Y1
   !   TPLUS = CA*LOG(YPLUS)+CB
   !ELSE
      ! log layer
      TPLUS = PR*RKAPPA*LOG(YPLUS)+B_T
   !ENDIF
ENDIF

H = RHO*U_TAU*CP/TPLUS

END SUBROUTINE HEAT_FLUX_MODEL


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
REAL(EB) :: C(2,2),SS(2),SLIP_COEF,ETA,AA,BB,U_STRM_0,DUMMY, &
            U_STRM,U_NORM,U_STRM_WALL,U_NORM_WALL,DPDS,DUSDS,DUSDN,TSN,RDN
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

RDN = 1._EB/DN

! ODE solution
IF (DNS) THEN
   ETA = U_NORM + RRHO*MU*RDN
   AA  = -(0.5_EB*DUSDS + TWTH*ETA*RDN)
   BB  = (TWTH*U_STRM_WALL*RDN + ONSI*DUSDN)*ETA - (U_NORM*0.5_EB*DUSDN + RRHO*( DPDS + TSN*0.5_EB*RDN ))
   !AA  = -0.5_EB*(DUSDS + ETA/DN)
   !BB  = 0.5_EB*US_WALL/DN*ETA - (UN*0.5_EB*DUSDN + RRHO*( DPDS + TSN/(2._EB*DN) ))
   U_STRM = ((AA*U_STRM + BB)*EXP(AA*DT) - BB)/AA
ELSE
   U_STRM_0 = U_STRM
   DO SUBIT=1,1
      CALL WALL_MODEL(SLIP_COEF,DUMMY,DUMMY,U_STRM-U_STRM_WALL,MU*RRHO,DN,0._EB)
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
USE MATH_FUNCTIONS, ONLY: CROSS_PRODUCT, NORM2

REAL(EB), INTENT(IN) :: U_VELO(3),U_SURF(3),NN(3),DN,DIVU,GRADU(3,3),GRADP(3),TAU_IJ(3,3),DT,RRHO,MU,ROUGHNESS,U_INT
INTEGER, INTENT(IN) :: I_VEL
REAL(EB) :: C(3,3),SS(3),PP(3),SLIP_COEF,ETA,AA,BB,U_STRM_0,DUMMY,U_RELA(3), &
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
IF (DNS) THEN
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
      CALL WALL_MODEL(SLIP_COEF,DUMMY,DUMMY,U_STRM,MU*RRHO,DN,ROUGHNESS)
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
   
   EDDY_LOOP: DO NE=1,VT%N_EDDY
      IERROR=1; CALL EDDY_POSITION(NE,NV,NM,IERROR)
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

CALL POINT_TO_MESH(NM)

IF (EVACUATION_ONLY(NM)) RETURN
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
            VT%X_EDDY(NE) = VT%X_EDDY(NE) + DT*SF%VEL_T(2)*PROFILE_FACTOR*RAMP_T
            VT%Y_EDDY(NE) = VT%Y_EDDY(NE) - DT*SF%VEL*SIGN(1._EB,REAL(VT%IOR,EB))*PROFILE_FACTOR*RAMP_T
            VT%Z_EDDY(NE) = VT%Z_EDDY(NE) + DT*SF%VEL_T(1)*PROFILE_FACTOR*RAMP_T
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
REAL :: RN
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
    CALL RANDOM_NUMBER(RN); VT%X_EDDY(NE) = VT%X_EDDY_MIN + REAL(RN,EB)*(VT%X_EDDY_MAX-VT%X_EDDY_MIN)
    CALL RANDOM_NUMBER(RN); VT%Y_EDDY(NE) = VT%Y_EDDY_MIN + REAL(RN,EB)*(VT%Y_EDDY_MAX-VT%Y_EDDY_MIN)
    CALL RANDOM_NUMBER(RN); VT%Z_EDDY(NE) = VT%Z_EDDY_MIN + REAL(RN,EB)*(VT%Z_EDDY_MAX-VT%Z_EDDY_MIN)
ENDIF

END SUBROUTINE EDDY_POSITION


SUBROUTINE EDDY_AMPLITUDE(NE,NV,NM)

INTEGER, INTENT(IN) :: NE,NV,NM
REAL(EB) :: EPS_EDDY(3)
REAL     :: RN
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

REAL(EB) :: DUMMY,MEANU,MEANV,MEANW
INTEGER :: I,J,K,II,JJ,KK,FILE_NUM,IM,IW,IOR,IIO,JJO,KKO,NX,NX_BLOCK,I_MESH,J_MESH,K_MESH,MX
INTEGER, INTENT(IN) :: NM
CHARACTER(80), INTENT(IN) :: FN_ISO
TYPE (MESH_TYPE), POINTER :: M=>NULL()
TYPE (MESH_TYPE), POINTER :: M2=>NULL()
TYPE (WALL_TYPE), POINTER :: WC=>NULL()

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

READ (FILE_NUM,*) DUMMY, DUMMY, DUMMY
READ (FILE_NUM,*) DUMMY, DUMMY, DUMMY

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
      IOR = WC%ONE_D%IOR
      II = WC%ONE_D%II
      JJ = WC%ONE_D%JJ
      KK = WC%ONE_D%KK

      ! assumes no refinement
      IIO = WC%NOM_IB(1)
      JJO = WC%NOM_IB(2)
      KKO = WC%NOM_IB(3)

      M2=>MESHES(WC%NOM)

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
1     ntot=2
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
5     m=n
      ntwo=np1
      iif=1
      idiv=2
10    iquot=m/idiv
      irem=m-idiv*iquot
      if(iquot-idiv.lt.0) goto 50
11    if(irem.ne.0) goto 20
12    ntwo=ntwo+ntwo
      ifact(iif)=idiv
      iif=iif+1
      m=iquot
      go to 10
20    idiv=3
      inon2=iif
30    iquot=m/idiv
      irem=m-idiv*iquot
      if(iquot-idiv.lt.0) goto 60
31    if(irem.ne.0) goto 40
32    ifact(iif)=idiv
      iif=iif+1
      m=iquot
      go to 30
40    idiv=idiv+2
      go to 30
50    inon2=iif
      if(irem.ne.0) goto 60
51    ntwo=ntwo+ntwo
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
71    if(iform.gt.0) goto 100
72    icase=2
      i1rng=np0*(1+nprev/2)
      if(idim-1.gt.0) goto 100
73    icase=3
      i1rng=np1
      if(ntwo-np1.le.0) goto 100
74    icase=4
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
110   np2hf=np2/2
      j=1
      do i2=1,np2,np1
      if(j-i2.ge.0) goto 130
120   i1max=i2+np1-2
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
145   j=j-m
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
210   work(i)=data(j)
      work(i+1)=data(j+1)
      go to 230
220   work(i)=data(j)
      work(i+1)=0._EB
230   ifp2=np2
      iif=ifmin
240   ifp1=ifp2/ifact(iif)
      j=j+ifp1
      if(j-i3-ifp2.lt.0) goto 260
250   j=j-ifp2
      ifp2=ifp1
      iif=iif+1
      if(ifp2-np1.gt.0) goto 240
260   continue
      i2max=i3+np2-np1
      i=1
      do i2=i3,i2max,np1
      data(i2)=work(i)
      data(i2+1)=work(i+1)
270   i=i+2
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
305   np1tw=np1+np1
      ipar=ntwo/np1
310   if(ipar-2.lt.0) then
         goto 350
      elseif(ipar-2.eq.0) then
         goto 330
      endif
320   ipar=ipar/4
      go to 310
330   do i1=1,i1rng,2
      do k1=i1,ntot,np1tw
      k2=k1+np1
      tempr=data(k2)
      tempi=data(k2+1)
      data(k2)=data(k1)-tempr
      data(k2+1)=data(k1+1)-tempi
      data(k1)=data(k1)+tempr
340   data(k1+1)=data(k1+1)+tempi
      enddo
      enddo
350   mmax=np1
360   if(mmax-ntwo/2.ge.0) goto 600
370   lmax=max0(np1tw,mmax/2)
      do 570 l=np1,lmax,np1tw
      m=l
      if(mmax-np1.le.0) goto 420
380   theta=-twopi*REAL(l,EB)/REAL(4*mmax,EB)
      if(isign.lt.0) goto 400
390   theta=-theta
400   wr=cos(theta)
      wi=sin(theta)
410   w2r=wr*wr-wi*wi
      w2i=2._EB*wr*wi
      w3r=w2r*wr-w2i*wi
      w3i=w2r*wi+w2i*wr
420   do 530 i1=1,i1rng,2
      kmin=i1+ipar*m
      if(mmax-np1.gt.0) goto 440
430   kmin=i1
440   kdif=ipar*mmax
450   kstep=4*kdif
      if(kstep-ntwo.gt.0) goto 530
460   do k1=kmin,ntot,kstep
      k2=k1+kdif
      k3=k2+kdif
      k4=k3+kdif
      if(mmax-np1.gt.0) goto 480
470   u1r=data(k1)+data(k2)
      u1i=data(k1+1)+data(k2+1)
      u2r=data(k3)+data(k4)
      u2i=data(k3+1)+data(k4+1)
      u3r=data(k1)-data(k2)
      u3i=data(k1+1)-data(k2+1)
      if(isign.ge.0) goto 472
471   u4r=data(k3+1)-data(k4+1)
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
490   u4r=t3i-t4i
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
520   data(k4+1)=u3i-u4i
      enddo
      kdif=kstep
      kmin=4*(kmin-i1)+i1
      go to 450
530   continue
      m=m+lmax
      if(m-mmax.gt.0) goto 570
540   if(isign.ge.0) goto 560
550   tempr=wr
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
605   ifp1=ntwo
      iif=inon2
      np1hf=np1/2
610   ifp2=ifact(iif)*ifp1
      j1min=np1+1
      if(j1min-ifp1.gt.0) goto 640
615   do j1=j1min,ifp1,np1
      theta=-twopi*REAL(j1-1,EB)/REAL(ifp2,EB)
      if(isign.lt.0) goto 625
620   theta=-theta
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
630   data(j3+1)=tempr*wi+data(j3+1)*wr
      enddo
      enddo
      tempr=wr
      wr=wr*wstpr-wi*wstpi
635   wi=tempr*wstpi+wi*wstpr
      enddo
      enddo
640   theta=-twopi/REAL(ifact(iif),EB)
      if(isign.lt.0) goto 650
645   theta=-theta
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
655   sumr=0._EB
      sumi=0._EB
      do j=jmin,jmax,ifp1
659   sumr=sumr+data(j)
660   sumi=sumi+data(j+1)
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
675   tempr=wr*sumr-oldsr+data(j)
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
685   wr=wstpr
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
695   i=i+2
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
702   theta=-theta
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
720   j=j+np2
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
730   if(isign.ge.0) goto 740
731   do i=imin,ntot,np2
735   data(i+1)=-data(i+1)
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
760   data(j)=data(imin)-data(imin+1)
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
775   data(j)=data(imin)+data(imin+1)
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
805   do i3=1,ntot,np2
      i2max=i3+np2-np1
      do i2=i3,i2max,np1
      imin=i2+i1rng
      imax=i2+np1-2
      jmax=2*i3+np1-imin
      if(i2-i3.le.0) goto 820
810   jmax=jmax+np2
820   if(idim-2.le.0) goto 850
830   j=jmax+np0
      do i=imin,imax,2
      data(i)=data(j)
      data(i+1)=-data(j+1)
840   j=j-2
      enddo
850   j=jmax
      do i=imin,imax,np0
      data(i)=data(j)
      data(i+1)=-data(j+1)
860   j=j-np0
      enddo
      enddo
      enddo
!
!     end of loop on each dimension
!
900   np0=np1
      np1=np2
910   nprev=n
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


SUBROUTINE GET_REV_turb(MODULE_REV,MODULE_DATE)
INTEGER,INTENT(INOUT) :: MODULE_REV
CHARACTER(255),INTENT(INOUT) :: MODULE_DATE

WRITE(MODULE_DATE,'(A)') turbrev(INDEX(turbrev,':')+2:LEN_TRIM(turbrev)-2)
READ (MODULE_DATE,'(I5)') MODULE_REV
WRITE(MODULE_DATE,'(A)') turbdate

END SUBROUTINE GET_REV_turb


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
PUBLIC :: SHUNN_MMS_3,VD2D_MMS_Z,VD2D_MMS_RHO,VD2D_MMS_RHO_OF_Z,VD2D_MMS_U,VD2D_MMS_V, &
          VD2D_MMS_Z_SRC,VD2D_MMS_U_SRC,VD2D_MMS_V_SRC,VD2D_MMS_Z_OF_RHO,VD2D_MMS_DIV

CONTAINS

SUBROUTINE SHUNN_MMS_3(NM)
USE MESH_POINTERS
USE MESH_VARIABLES
USE GLOBAL_CONSTANTS, ONLY: P_INF,TMPA,N_TRACKED_SPECIES
USE PHYSICAL_FUNCTIONS, ONLY: GET_SPECIFIC_GAS_CONSTANT
IMPLICIT NONE
INTEGER, INTENT(IN) :: NM
INTEGER :: I,J,K
REAL(EB) :: ZZ_GET(1:N_TRACKED_SPECIES)
REAL(EB), PARAMETER :: GAS_CONSTANT=8314.472_EB

! set up problem parameters

RHO_0_MMS = P_INF*SPECIES_MIXTURE(2)%MW/(GAS_CONSTANT*TMPA)
RHO_1_MMS = P_INF*SPECIES_MIXTURE(1)%MW/(GAS_CONSTANT*TMPA)

! initialize private variables

CALL VD2D_MMS_INIT

! initialize flow fields

CALL POINT_TO_MESH(NM)
U = UF_MMS
V = 0._EB
W = WF_MMS
DO K=0,KBP1
   DO J=0,JBP1
      DO I=0,IBP1
         ZZ(I,J,K,1) = VD2D_MMS_Z(XC(I),ZC(K),0._EB)
         ZZ_GET(1) = ZZ(I,J,K,1)
         ZZ(I,J,K,2) = 1._EB - ZZ(I,J,K,1)
         ZZ_GET(2) = ZZ(I,J,K,2)
         CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,RSUM(I,J,K))
         RHO(I,J,K) = P_INF/(TMPA*RSUM(I,J,K))
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
VD2D_MMS_U=-W/K/4.0_EB * COS(K * PI * X) * SIN(K * PI * Y) * SIN(W * PI * T) &
   * (R1-R0)/VD2D_MMS_RHO(X,Y,T)
END FUNCTION VD2D_MMS_U

FUNCTION VD2D_MMS_V(X,Y,T)
IMPLICIT NONE
REAL(EB),INTENT (IN) :: X,Y,T
REAL(EB) :: VD2D_MMS_V
VD2D_MMS_V = VD2D_MMS_U(Y,X,T)
!! NOTE: V OBTAINED FROM U BY TRANSPOSING X <--> Y
END FUNCTION VD2D_MMS_V

FUNCTION VD2D_MMS_P(X,Y,T)
IMPLICIT NONE
REAL(EB),INTENT (IN) :: X,Y,T
REAL(EB) :: VD2D_MMS_P
VD2D_MMS_P = &
   VD2D_MMS_RHO(X,Y,T) * VD2D_MMS_U(X,Y,T) * VD2D_MMS_V(X,Y,T)/2.0_EB
   !! NOTE: PRESSURE CAN VARY BY A CONSTANT AND STILL SATISFY THE MMS
END FUNCTION VD2D_MMS_P

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

END MODULE MANUFACTURED_SOLUTIONS


