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
PUBLIC :: NS_ANALYTICAL_SOLUTION, INIT_TURB_ARRAYS, VARDEN_DYNSMAG, WANNIER_FLOW, &
          GET_REV_turb, WALL_MODEL, COMPRESSION_WAVE, VELTAN2D,VELTAN3D, &
          SYNTHETIC_TURBULENCE, SYNTHETIC_EDDY_SETUP, TEST_FILTER, EX2G3D, TENSOR_DIFFUSIVITY_MODEL, &
          TWOD_VORTEX_CERFACS, HEAT_FLUX_MODEL, ABL_HEAT_FLUX_MODEL
 
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
   
! 1D working arrays
IF (LES) THEN
   ALLOCATE(M%TURB_WORK11(0:MAX(IBP1,JBP1,KBP1)),STAT=IZERO)
   CALL ChkMemErr('INIT_TURB_ARRAYS','TURB_WORK11',IZERO)
   M%TURB_WORK11 = 0._EB
   ALLOCATE(M%TURB_WORK12(0:MAX(IBP1,JBP1,KBP1)),STAT=IZERO)
   CALL ChkMemErr('INIT_TURB_ARRAYS','TURB_WORK12',IZERO)
   M%TURB_WORK12 = 0._EB
ENDIF

END SUBROUTINE INIT_TURB_ARRAYS


SUBROUTINE NS_ANALYTICAL_SOLUTION(NM)

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

CALL TEST_FILTER(SHAT11,S11,-1.E10_EB,1.E10_EB)
CALL TEST_FILTER(SHAT22,S22,-1.E10_EB,1.E10_EB)
CALL TEST_FILTER(SHAT33,S33,-1.E10_EB,1.E10_EB)
CALL TEST_FILTER(SHAT12,S12,-1.E10_EB,1.E10_EB)
CALL TEST_FILTER(SHAT13,S13,-1.E10_EB,1.E10_EB)
CALL TEST_FILTER(SHAT23,S23,-1.E10_EB,1.E10_EB)


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

CALL TEST_FILTER(BETAHAT11,BETA11,-1.E10_EB,1.E10_EB)
CALL TEST_FILTER(BETAHAT22,BETA22,-1.E10_EB,1.E10_EB)
CALL TEST_FILTER(BETAHAT33,BETA33,-1.E10_EB,1.E10_EB)
CALL TEST_FILTER(BETAHAT12,BETA12,-1.E10_EB,1.E10_EB)
CALL TEST_FILTER(BETAHAT13,BETA13,-1.E10_EB,1.E10_EB)
CALL TEST_FILTER(BETAHAT23,BETA23,-1.E10_EB,1.E10_EB)

! test filter the density

RHOPHAT => WORK7
CALL TEST_FILTER(RHOPHAT,RHOP,RHOMIN,RHOMAX)

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

CALL TEST_FILTER(MLHAT,ML,0._EB,1.E10_EB)
CALL TEST_FILTER(MMHAT,MM,0._EB,1.E10_EB)

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

CALL TEST_FILTER(RUU_HAT,RUU,0.E10_EB,1.E10_EB)
CALL TEST_FILTER(RVV_HAT,RVV,0.E10_EB,1.E10_EB)
CALL TEST_FILTER(RWW_HAT,RWW,0.E10_EB,1.E10_EB)
CALL TEST_FILTER(RUV_HAT,RUV,-1.E10_EB,1.E10_EB)
CALL TEST_FILTER(RUW_HAT,RUW,-1.E10_EB,1.E10_EB)
CALL TEST_FILTER(RVW_HAT,RVW,-1.E10_EB,1.E10_EB)

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

CALL TEST_FILTER(RU_HAT,RU,-1.E10_EB,1.E10_EB)
CALL TEST_FILTER(RV_HAT,RV,-1.E10_EB,1.E10_EB)
CALL TEST_FILTER(RW_HAT,RW,-1.E10_EB,1.E10_EB)

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


SUBROUTINE TEST_FILTER(PHIBAR,PHI,PHI_MIN,PHI_MAX)

REAL(EB), INTENT(IN) :: PHI_MIN,PHI_MAX
REAL(EB), INTENT(IN) :: PHI(0:IBP1,0:JBP1,0:KBP1)
REAL(EB), INTENT(OUT) :: PHIBAR(0:IBP1,0:JBP1,0:KBP1)
REAL(EB), POINTER, DIMENSION(:) :: PHI1,PHI2
INTEGER I,J,K

PHI1 => TURB_WORK11
PHI2 => TURB_WORK12

! filter in x:
DO K = 0,KBP1
   DO J = 0,JBP1
      PHI1(0:IBP1) = PHI(0:IBP1,J,K)
      CALL TOPHAT_FILTER_1D(PHI2(0:IBP1),PHI1(0:IBP1),0,IBP1,PHI_MIN,PHI_MAX)
      PHIBAR(0:IBP1,J,K) = PHI2(0:IBP1)
   ENDDO
ENDDO

IF (.NOT.TWO_D) THEN
   ! filter in y:
   DO K = 0,KBP1
      DO I = 0,IBP1
         PHI1(0:JBP1) = PHIBAR(I,0:JBP1,K)
         CALL TOPHAT_FILTER_1D(PHI2(0:JBP1),PHI1(0:JBP1),0,JBP1,PHI_MIN,PHI_MAX)
         PHIBAR(I,0:JBP1,K) = PHI2(0:JBP1)
      ENDDO
   ENDDO
ENDIF

! filter in z:
DO J = 0,JBP1
   DO I = 0,IBP1
      PHI1(0:KBP1) = PHIBAR(I,J,0:KBP1)
      CALL TOPHAT_FILTER_1D(PHI2(0:KBP1),PHI1(0:KBP1),0,KBP1,PHI_MIN,PHI_MAX)
      PHIBAR(I,J,0:KBP1) = PHI2(0:KBP1)
   ENDDO
ENDDO

END SUBROUTINE TEST_FILTER


SUBROUTINE TOPHAT_FILTER_1D(UBAR,U,N_LO,N_HI,U_MIN,U_MAX)

INTEGER, INTENT(IN) :: N_LO,N_HI
REAL(EB), INTENT(IN) :: U(N_LO:N_HI),U_MIN,U_MAX
REAL(EB), INTENT(OUT) :: UBAR(N_LO:N_HI)
INTEGER :: J
!REAL(EB),PARAMETER:: W(-1:1) = (/0.25_EB,0.5_EB,0.25_EB/)   ! trapezoid rule
!REAL(EB),PARAMETER:: W(-1:1) = (/ONSI,TWTH,ONSI/)           ! Simpson's rule

! Filter the u field to obtain ubar
DO J=N_LO+1,N_HI-1
   !UBAR(J) = DOT_PRODUCT(W(-1:1),U(J-1:J+1))
   UBAR(J) = 0.5_EB*U(J) + 0.25_EB*(U(J-1)+U(J+1))
ENDDO

! set boundary values (not ideal, but fast and simple)
!UBAR(N_LO) = UBAR(N_LO+1)
!UBAR(N_HI) = UBAR(N_HI-1)
UBAR(N_LO) = MIN(U_MAX,MAX(U_MIN,2._EB*UBAR(N_LO+1)-UBAR(N_LO+2)))
UBAR(N_HI) = MIN(U_MAX,MAX(U_MIN,2._EB*UBAR(N_HI-1)-UBAR(N_HI-2)))

END SUBROUTINE TOPHAT_FILTER_1D


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
!
! Comments:
!
! The slip factor (SF) is based on the following approximation to the wall gradient
! (note that u0 is the ghost cell value of the streamwise velocity component and
! y is the wall-normal direction):
! dudy = (u-u0)/dy = (u-SF*u)/dy = u/dy*(1-SF) => SF = 1 - dudy*dy/u
! In this routine, dudy is sampled from the wall model at the location y_cell_center.

! New scheme

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
USE TRAN, ONLY : GET_IJK

REAL(EB), INTENT(IN) :: DT,T
INTEGER, INTENT(IN) :: NM
INTEGER :: NE,NV,II,JJ,KK,IERROR,IC,IW
TYPE(VENTS_TYPE), POINTER :: VT=>NULL()
TYPE(SURFACE_TYPE), POINTER :: SF=>NULL()
REAL(EB) :: XX,YY,ZZ,SHAPE_FACTOR,VOLUME_WEIGHTING_FACTOR(3),EDDY_VOLUME(3),PROFILE_FACTOR,RAMP_T,TSI,XI,YJ,ZK
REAL(EB), PARAMETER :: EPSDX=1.E-10_EB

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
            ! This IF block is needed because eddies may be centered outside the vent.
            IF (VT%Y_EDDY(NE)>VT%Y1 .AND. VT%Y_EDDY(NE)<VT%Y2 .AND. &
                VT%Z_EDDY(NE)>VT%Z1 .AND. VT%Z_EDDY(NE)<VT%Z2 .AND. &
                VT%Y_EDDY(NE)>YS .AND. VT%Y_EDDY(NE)<YF .AND. &
                VT%Z_EDDY(NE)>ZS .AND. VT%Z_EDDY(NE)<ZF .AND. &
                ABS(SF%VEL)>TWO_EPSILON_EB) THEN
               SELECT CASE (VT%IOR)
                  CASE( 1); CALL GET_IJK(VT%X1+EPSDX,VT%Y_EDDY(NE),VT%Z_EDDY(NE),NM,XI,YJ,ZK,II,JJ,KK)
                  CASE(-1); CALL GET_IJK(VT%X1-EPSDX,VT%Y_EDDY(NE),VT%Z_EDDY(NE),NM,XI,YJ,ZK,II,JJ,KK)
               END SELECT
               IC = CELL_INDEX(II,JJ,KK)
               IW = WALL_INDEX(IC,-VT%IOR)
               PROFILE_FACTOR = ABS(WALL(IW)%UW0/SF%VEL)
            ELSE
               PROFILE_FACTOR = 0._EB
            ENDIF
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
                  SHAPE_FACTOR = SHAPE_FUNCTION(XX,1)*SHAPE_FUNCTION(YY,1)*SHAPE_FUNCTION(ZZ,1)
                  VT%U_EDDY(JJ,KK) = VT%U_EDDY(JJ,KK) + VT%CU_EDDY(NE)*SHAPE_FACTOR
                  
                  XX = (VT%X1  - VT%X_EDDY(NE))/VT%SIGMA_IJ(2,1)
                  YY = (YC(JJ) - VT%Y_EDDY(NE))/VT%SIGMA_IJ(2,2)
                  ZZ = (ZC(KK) - VT%Z_EDDY(NE))/VT%SIGMA_IJ(2,3)
                  SHAPE_FACTOR = SHAPE_FUNCTION(XX,1)*SHAPE_FUNCTION(YY,1)*SHAPE_FUNCTION(ZZ,1)
                  VT%V_EDDY(JJ,KK) = VT%V_EDDY(JJ,KK) + VT%CV_EDDY(NE)*SHAPE_FACTOR
                  
                  XX = (VT%X1  - VT%X_EDDY(NE))/VT%SIGMA_IJ(3,1)
                  YY = (YC(JJ) - VT%Y_EDDY(NE))/VT%SIGMA_IJ(3,2)
                  ZZ = (ZC(KK) - VT%Z_EDDY(NE))/VT%SIGMA_IJ(3,3)
                  SHAPE_FACTOR = SHAPE_FUNCTION(XX,1)*SHAPE_FUNCTION(YY,1)*SHAPE_FUNCTION(ZZ,1)
                  VT%W_EDDY(JJ,KK) = VT%W_EDDY(JJ,KK) + VT%CW_EDDY(NE)*SHAPE_FACTOR
               ENDDO
            ENDDO
         ENDDO EDDY_LOOP_1
      CASE(2)
         EDDY_LOOP_2: DO NE=1,VT%N_EDDY
            IF (VT%X_EDDY(NE)>VT%X1 .AND. VT%X_EDDY(NE)<VT%X2 .AND. &
                VT%Z_EDDY(NE)>VT%Z1 .AND. VT%Z_EDDY(NE)<VT%Z2 .AND. &
                VT%X_EDDY(NE)>XS .AND. VT%X_EDDY(NE)<XF .AND. &
                VT%Z_EDDY(NE)>ZS .AND. VT%Z_EDDY(NE)<ZF .AND. &
                ABS(SF%VEL)>TWO_EPSILON_EB) THEN
               SELECT CASE (VT%IOR)
                  CASE( 2); CALL GET_IJK(VT%X_EDDY(NE),VT%Y1+EPSDX,VT%Z_EDDY(NE),NM,XI,YJ,ZK,II,JJ,KK)
                  CASE(-2); CALL GET_IJK(VT%X_EDDY(NE),VT%Y1-EPSDX,VT%Z_EDDY(NE),NM,XI,YJ,ZK,II,JJ,KK)
               END SELECT
               IC = CELL_INDEX(II,JJ,KK)
               IW = WALL_INDEX(IC,-VT%IOR)
               PROFILE_FACTOR = ABS(WALL(IW)%UW0/SF%VEL)
            ELSE
               PROFILE_FACTOR = 0._EB
            ENDIF
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
                  SHAPE_FACTOR = SHAPE_FUNCTION(XX,1)*SHAPE_FUNCTION(YY,1)*SHAPE_FUNCTION(ZZ,1)
                  VT%U_EDDY(II,KK) = VT%U_EDDY(II,KK) + VT%CU_EDDY(NE)*SHAPE_FACTOR
                  
                  XX = (XC(II) - VT%X_EDDY(NE))/VT%SIGMA_IJ(2,1)
                  YY = (VT%Y1  - VT%Y_EDDY(NE))/VT%SIGMA_IJ(2,2)
                  ZZ = (ZC(KK) - VT%Z_EDDY(NE))/VT%SIGMA_IJ(2,3)
                  SHAPE_FACTOR = SHAPE_FUNCTION(XX,1)*SHAPE_FUNCTION(YY,1)*SHAPE_FUNCTION(ZZ,1)
                  VT%V_EDDY(II,KK) = VT%V_EDDY(II,KK) + VT%CV_EDDY(NE)*SHAPE_FACTOR
                  
                  XX = (XC(II) - VT%X_EDDY(NE))/VT%SIGMA_IJ(3,1)
                  YY = (VT%Y1  - VT%Y_EDDY(NE))/VT%SIGMA_IJ(3,2)
                  ZZ = (ZC(KK) - VT%Z_EDDY(NE))/VT%SIGMA_IJ(3,3)
                  SHAPE_FACTOR = SHAPE_FUNCTION(XX,1)*SHAPE_FUNCTION(YY,1)*SHAPE_FUNCTION(ZZ,1)
                  VT%W_EDDY(II,KK) = VT%W_EDDY(II,KK) + VT%CW_EDDY(NE)*SHAPE_FACTOR 
               ENDDO
            ENDDO  
         ENDDO EDDY_LOOP_2
      CASE(3)
         EDDY_LOOP_3: DO NE=1,VT%N_EDDY
            IF (VT%X_EDDY(NE)>VT%X1 .AND. VT%X_EDDY(NE)<VT%X2 .AND. &
                VT%Y_EDDY(NE)>VT%Y1 .AND. VT%Y_EDDY(NE)<VT%Y2 .AND. &
                VT%X_EDDY(NE)>XS .AND. VT%X_EDDY(NE)<XF .AND. &
                VT%Y_EDDY(NE)>YS .AND. VT%Y_EDDY(NE)<YF .AND. &
                ABS(SF%VEL)>TWO_EPSILON_EB) THEN
               SELECT CASE (VT%IOR)
                  CASE( 3); CALL GET_IJK(VT%X_EDDY(NE),VT%Y_EDDY(NE),VT%Z1+EPSDX,NM,XI,YJ,ZK,II,JJ,KK)
                  CASE(-3); CALL GET_IJK(VT%X_EDDY(NE),VT%Y_EDDY(NE),VT%Z1-EPSDX,NM,XI,YJ,ZK,II,JJ,KK)
               END SELECT
               IC = CELL_INDEX(II,JJ,KK)
               IW = WALL_INDEX(IC,-VT%IOR)
               PROFILE_FACTOR = ABS(WALL(IW)%UW0/SF%VEL)
            ELSE
               PROFILE_FACTOR = 0._EB
            ENDIF
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
                  SHAPE_FACTOR = SHAPE_FUNCTION(XX,1)*SHAPE_FUNCTION(YY,1)*SHAPE_FUNCTION(ZZ,1)
                  VT%U_EDDY(II,JJ) = VT%U_EDDY(II,JJ) + VT%CU_EDDY(NE)*SHAPE_FACTOR
                  
                  XX = (XC(II) - VT%X_EDDY(NE))/VT%SIGMA_IJ(2,1)
                  YY = (YC(JJ) - VT%Y_EDDY(NE))/VT%SIGMA_IJ(2,2)
                  ZZ = (VT%Z1  - VT%Z_EDDY(NE))/VT%SIGMA_IJ(2,3)
                  SHAPE_FACTOR = SHAPE_FUNCTION(XX,1)*SHAPE_FUNCTION(YY,1)*SHAPE_FUNCTION(ZZ,1)
                  VT%V_EDDY(II,JJ) = VT%V_EDDY(II,JJ) + VT%CV_EDDY(NE)*SHAPE_FACTOR
                  
                  XX = (XC(II) - VT%X_EDDY(NE))/VT%SIGMA_IJ(3,1)
                  YY = (YC(JJ) - VT%Y_EDDY(NE))/VT%SIGMA_IJ(3,2)
                  ZZ = (VT%Z1  - VT%Z_EDDY(NE))/VT%SIGMA_IJ(3,3)
                  SHAPE_FACTOR = SHAPE_FUNCTION(XX,1)*SHAPE_FUNCTION(YY,1)*SHAPE_FUNCTION(ZZ,1)
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
CALL RANDOM_NUMBER(RN); IF (RN>0.5_8) EPS_EDDY(1)=1._EB
CALL RANDOM_NUMBER(RN); IF (RN>0.5_8) EPS_EDDY(2)=1._EB
CALL RANDOM_NUMBER(RN); IF (RN>0.5_8) EPS_EDDY(3)=1._EB

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

REAL(EB), PRIVATE :: K,W,RHO0,RHO1,D,MU
REAL(EB), PUBLIC :: RHO_0_MMS,RHO_1_MMS
REAL(EB), PUBLIC, PARAMETER :: UF_MMS=0.5_EB,WF_MMS=0.5_EB,VISC_MMS=0.001_EB,DIFF_MMS=0.001_EB,WAVE_NUM_MMS=2._EB,FREQ_MMS=2._EB

PRIVATE
PUBLIC :: SHUNN_MMS_3,VD2D_MMS_INIT,VD2D_MMS_Z,VD2D_MMS_RHO,VD2D_MMS_RHO_OF_Z,VD2D_MMS_U,VD2D_MMS_V, &
          VD2D_MMS_Z_SRC,VD2D_MMS_U_SRC,VD2D_MMS_V_SRC,VD2D_MMS_Z_OF_RHO,VD2D_MMS_DIV,VD2D_MMS_RHO_SRC

CONTAINS

SUBROUTINE SHUNN_MMS_3(NM)
USE MESH_POINTERS
USE MESH_VARIABLES
USE GLOBAL_CONSTANTS
USE PHYSICAL_FUNCTIONS, ONLY: GET_SPECIFIC_GAS_CONSTANT
IMPLICIT NONE
INTEGER, INTENT(IN) :: NM
INTEGER :: I,J,K
REAL(EB) :: ZZ_GET(0:N_TRACKED_SPECIES)
REAL(EB), PARAMETER :: WAVE_NUM=2._EB,FREQ=2._EB

! set up problem parameters

RHO_0_MMS = P_INF*SPECIES_MIXTURE(0)%MW/(R0*TMPA)
RHO_1_MMS = P_INF*SPECIES_MIXTURE(1)%MW/(R0*TMPA)
ZZ_GET(1) = 1._EB
CALL VD2D_MMS_INIT(WAVE_NUM,FREQ,RHO_0_MMS,RHO_1_MMS,VISC_MMS,DIFF_MMS)

! initialize flow fields

CALL POINT_TO_MESH(NM)
U = 0._EB
V = 0._EB
W = 0._EB
US = 0._EB
VS = 0._EB
WS = 0._EB
DO K=0,KBP1
   DO J=0,JBP1
      DO I=0,IBP1
         ZZ(I,J,K,1) = VD2D_MMS_Z(XC(I),ZC(K),0._EB)
         ZZ_GET(1) = ZZ(I,J,K,1)
         CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,RSUM(I,J,K))
         RHO(I,J,K) = P_INF/(TMPA*RSUM(I,J,K))
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE SHUNN_MMS_3

SUBROUTINE VD2D_MMS_INIT(WAVE_NUM,FREQ,RHO_0,RHO_1,VISC,DIFF)
IMPLICIT NONE
REAL(EB),INTENT (IN) :: WAVE_NUM,FREQ,RHO_0,RHO_1,VISC,DIFF
K = WAVE_NUM
W = FREQ
RHO0 = RHO_0
RHO1 = RHO_1
MU = VISC
D = DIFF
END SUBROUTINE VD2D_MMS_INIT

FUNCTION VD2D_MMS_RHO_OF_Z(Z)
IMPLICIT NONE
REAL(EB),INTENT (IN) :: Z
REAL(EB) :: VD2D_MMS_RHO_OF_Z
VD2D_MMS_RHO_OF_Z = 1.0_EB/(Z/RHO1+(1.0_EB-Z)/RHO0)
END FUNCTION VD2D_MMS_RHO_OF_Z

REAL(EB) FUNCTION VD2D_MMS_Z_OF_RHO(RHO)
IMPLICIT NONE
REAL(EB),INTENT (IN) :: RHO
VD2D_MMS_Z_OF_RHO = (RHO0-RHO)/(RHO0-RHO1)*RHO1/RHO
END FUNCTION VD2D_MMS_Z_OF_RHO

FUNCTION VD2D_MMS_Z(X,Y,T)
IMPLICIT NONE
REAL(EB),INTENT (IN) :: X,Y,T
REAL(EB) :: VD2D_MMS_Z
VD2D_MMS_Z = &
   (1.0_EB + SIN(K * PI * X) * SIN(K * PI * Y) *COS (W * PI * T))/ &
   (1.0_EB + RHO0/RHO1+(1.0_EB-RHO0/RHO1) &
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
   * (RHO1-RHO0)/VD2D_MMS_RHO(X,Y,T)
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
S6 = RHO1-RHO0
S7 = 1.0_EB + S3 * S5 * S0
S8 = (0.5_EB + 0.5_EB * S3 * S5 * S0) * S6 + RHO0! RHO
S9 = RHO1 * S3 * S5 * S0 + RHO1-S3 * S5 * S0 * RHO0 + RHO0
S10 = RHO1 * S2 * K * PI * S5 * S0-S2 * K * PI * S5 * S0 * RHO0
S11 = RHO1 * S3 * S4 * K * PI * S0-S3 * S4 * K * PI * S0 * RHO0
S12 = RHO0 * S3 * S5 * S1 * W * PI-S3 * S5 * S1 * W * PI * RHO1
S13 = S3 * K ** 2 * PI ** 2 * S5 * S0 * RHO0-RHO1 * S3 * K ** 2 * PI ** 2 * S5 * S0
S14 = 1.0_EB/4.0_EB
S15 = 2.0_EB * D
S16 = S15 * PI * PI
VD2D_MMS_Z_SRC=&
   -RHO1 * (S8 * S3 * S5 * S1 * W * PI * S9 ** 2 * K + S8 * S7 * S12 * S9 * K + S14 * S2 ** 2 * PI * S5 ** 2 &
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
S6 = RHO1-RHO0
S7=(0.5_EB + 0.5_EB * S3 * S5 * S0) * S6 + RHO0 ! RHO
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
S6 = RHO1-RHO0
S7 = VD2D_MMS_RHO(X,Y,T)
S8 = -W/(4._EB*K)

S9 = (1._EB + RHO0/RHO1)
S10 = S3 * S5 * S0
S11 = 1._EB + S10
S12 = S9 + S9 * S10

DZDX = S11 * (-S9 * S2 * S5 * S0 * PI * K) / S12**2 + (S2 * S5 * S0 * PI * K) / S12
DZDY = S11 * (-S9 * S3 * S4 * S0 * PI * K) / S12**2 + (S3 * S4 * S0 * PI * K) / S12

S13 = (1._EB/RHO1 - 1._EB/RHO0)*DZDX
S14 = (1._EB/RHO1 - 1._EB/RHO0)*DZDY

DUDX = S6 * S8 * S5 * S1 * (S2 * S13 - PI * K * S3 / S7)
DVDY = S6 * S8 * S3 * S1 * (S4 * S14 - PI * K * S5 / S7)

VD2D_MMS_DIV = DUDX + DVDY

END FUNCTION VD2D_MMS_DIV

REAL(EB) FUNCTION VD2D_MMS_RHO_SRC(X,Y,T)
IMPLICIT NONE
REAL(EB),INTENT (IN) :: X,Y,T
REAL(EB) :: UF,VF,R0,R1,K,W

UF = UF_MMS
VF = WF_MMS
R0 = RHO_0_MMS
R1 = RHO_1_MMS
K = WAVE_NUM_MMS
W = FREQ_MMS

VD2D_MMS_RHO_SRC = &
   PI*K*(-R0*UF*SIN(PI*K*(T*VF - Y))*COS(PI*K*(T*UF - X)) - R0*VF*SIN(PI*K*(T*UF - X))*COS(PI*K*(T*VF - Y)) + &
   R1*UF*SIN(PI*K*(T*VF - Y))*COS(PI*K*(T*UF - X)) + R1*VF*SIN(PI*K*(T*UF - X))*COS(PI*K*(T*VF - Y)))*COS(PI*T*W)/2._EB

END FUNCTION VD2D_MMS_RHO_SRC

END MODULE MANUFACTURED_SOLUTIONS



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
!IMPLICIT NONE
!
!
!PRIVATE
!PUBLIC SCALARF_EMB,VELOCITY_EMB,RESTRICT_MASS_EMB,RESTRICT_DIV_EMB,PROJECT_VELOCITY, &
!       SORT_MESH_LEVEL,MATCH_VELOCITY_EMB,SCALAR_GHOST_EMB
! 
!CONTAINS
!
!
!SUBROUTINE SCALARF_EMB(NM1,NM2,IERROR)
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
!SPECIES_LOOP: DO N=0,N_TRACKED_SPECIES
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
!
!INTEGER, INTENT(IN) :: NM1,NM2
!
!TYPE(MESH_TYPE), POINTER :: M1,M2
!INTEGER :: N,I,J,K,I_LO,I_HI,J_LO,J_HI,K_LO,K_HI,II_0,JJ_0,KK_0,II,JJ,KK, &
!           NRX,NRY,NRZ, II_LO,JJ_LO,KK_LO,INDEX_LIST(12),IERROR
!REAL(EB) :: DV1,DV2,DVRAT,VOLUME_LIST(3)
!REAL(EB), POINTER, DIMENSION(:,:,:) :: RHO1,RHO2
!REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZ1,ZZ2
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
!   IF (N_TRACKED_SPECIES>0) ZZ1  => M1%ZZS
!   IF (N_TRACKED_SPECIES>0) ZZ2  => M2%ZZS
!ELSEIF (CORRECTOR) THEN
!   RHO1 => M1%RHO
!   RHO2 => M2%RHO
!   IF (N_TRACKED_SPECIES>0) ZZ1  => M1%ZZ
!   IF (N_TRACKED_SPECIES>0) ZZ2  => M2%ZZ
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
!IF (N_TRACKED_SPECIES>0) THEN
!
!   SPECIES_LOOP: DO N=1,N_TRACKED_SPECIES
!   
!      DO K = K_LO,K_HI
!         KK_0 = KK_LO + (K-K_LO)*NRZ
!         DO J = J_LO,J_HI
!            JJ_0 = JJ_LO + (J-J_LO)*NRY
!            DO I = I_LO,I_HI
!               II_0 = II_LO + (I-I_LO)*NRX
!            
!               ZZ1(I,J,K,N) = 0._EB
!         
!               DO KK = KK_0+1,KK_0+NRZ
!                  DO JJ = JJ_0+1,JJ_0+NRY
!                     DO II = II_0+1,II_0+NRX
!                 
!                        ZZ1(I,J,K,N) = ZZ1(I,J,K,N) + RHO2(II,JJ,KK)*ZZ2(II,JJ,KK,N)*DV2
!                     
!                     ENDDO
!                  ENDDO
!               ENDDO
!               
!               ZZ1(I,J,K,N) = ZZ1(I,J,K,N)/(RHO1(I,J,K)*DV1)
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
!FINE_MESH_WALL_LOOP: DO IW=1,M2%N_EXTERNAL_WALL_CELLS
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
!
!INTEGER, INTENT(IN) :: NM1,NM2
!
!TYPE(MESH_TYPE), POINTER :: M1,M2
!INTEGER :: N,I,J,K,I_LO,I_HI,J_LO,J_HI,K_LO,K_HI,II_0,JJ_0,KK_0,II,JJ,KK, &
!           NRX,NRY,NRZ,II_LO,JJ_LO,KK_LO,INDEX_LIST(12),IERROR,IW
!REAL(EB) :: VOLUME_LIST(3)
!REAL(EB), POINTER, DIMENSION(:,:,:) :: RHOP1,RHOP2,TMP1,TMP2,HH1,HH2
!REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZP1,ZZP2
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
!   ZZP1  => M1%ZZS
!   HH1   => M1%H
!   
!   RHOP2 => M2%RHOS
!   ZZP2  => M2%ZZS
!   HH2   => M2%H
!ELSEIF (CORRECTOR) THEN
!   RHOP1 => M1%RHO
!   ZZP1  => M1%ZZ
!   HH1   => M1%HS
!   
!   RHOP2 => M2%RHO
!   ZZP2  => M2%ZZ
!   HH2   => M2%HS
!ENDIF
!TMP1 => M1%TMP
!TMP2 => M2%TMP
!
!
!! Set fine mesh boundary value to corresponding coarse mesh value
!
!SPECIES_LOOP: DO N=1,N_TRACKED_SPECIES
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
!                  IF (N_TRACKED_SPECIES>0) ZZP2(II_0,JJ,KK,N) = ZZP1(I,J,K,N)
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
!                  IF (N_TRACKED_SPECIES>0) ZZP2(II_0,JJ,KK,N) = ZZP1(I,J,K,N)
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
!                  IF (N_TRACKED_SPECIES>0) ZZP2(II,JJ_0,KK,N) = ZZP1(I,J,K,N)
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
!                  IF (N_TRACKED_SPECIES>0) ZZP2(II,JJ_0,KK,N) = ZZP1(I,J,K,N)
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
!                  IF (N_TRACKED_SPECIES>0) ZZP2(II,JJ,KK_0,N) = ZZP1(I,J,K,N)
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
!                  IF (N_TRACKED_SPECIES>0) ZZP2(II,JJ,KK_0,N) = ZZP1(I,J,K,N)
!               ENDDO
!            ENDDO
!         ENDIF
!         
!      ENDDO
!   ENDDO
!
!   WALL_LOOP: DO IW=1,M2%N_EXTERNAL_WALL_CELLS+M2%N_INTERNAL_WALL_CELLS
!      IF (M2%BOUNDARY_TYPE(IW)/=INTERPOLATED_BOUNDARY) CYCLE WALL_LOOP
!      II = M2%IJKW(1,IW)
!      JJ = M2%IJKW(2,IW)
!      KK = M2%IJKW(3,IW)
!      M2%RHO_F(IW)  = RHOP2(II,JJ,KK) 
!      M2%ZZ_F(IW,N) = ZZP2(II,JJ,KK,N)
!   ENDDO WALL_LOOP
!
!ENDDO SPECIES_LOOP
!
!END SUBROUTINE SCALAR_GHOST_EMB
!
!
!SUBROUTINE LOCATE_MESH(INDEX_LIST,VOLUME_LIST,NM1,NM2,IERROR)
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

