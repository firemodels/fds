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
PUBLIC :: NS_ANALYTICAL_SOLUTION, INIT_TURB_ARRAYS, VARDEN_DYNSMAG, &
          GET_REV_turb, WERNER_WENGLE_WALL_MODEL, COMPRESSION_WAVE, VELTAN2D,VELTAN3D,STRATIFIED_MIXING_LAYER, &
          SURFACE_HEAT_FLUX_MODEL, SYNTHETIC_TURBULENCE, SYNTHETIC_EDDY_SETUP, TEST_FILTER, EX2G3D, TENSOR_DIFFUSIVITY_MODEL, &
          TWOD_VORTEX_CERFACS
 
CONTAINS


SUBROUTINE INIT_TURB_ARRAYS(NM)

USE MEMORY_FUNCTIONS, ONLY: ChkMemErr
INTEGER, INTENT(IN) :: NM
INTEGER :: IZERO
TYPE (MESH_TYPE), POINTER :: M

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


SUBROUTINE STRATIFIED_MIXING_LAYER(NM)

USE PHYSICAL_FUNCTIONS, ONLY: GET_SPECIFIC_GAS_CONSTANT
INTEGER, INTENT(IN) :: NM
INTEGER :: I,J,K
REAL(EB) :: ZZ_GET(0:N_TRACKED_SPECIES)

CALL POINT_TO_MESH(NM)

DO K=0,KBP1
   U(:,:,K)   =               0.5_EB*MIXING_LAYER_U0*    TANH(2._EB*ZC(K)/MIXING_LAYER_H0)
   RHO(:,:,K) = RHOA*(1._EB - 0.5_EB*MIXING_LAYER_THETA0*TANH(2._EB*ZC(K)/MIXING_LAYER_H0))
ENDDO

! initialize temperature
DO K=0,KBP1
   DO J=0,JBP1
      DO I=0,IBP1
         IF (N_TRACKED_SPECIES>0) ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(I,J,K,1:N_TRACKED_SPECIES)
         CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,RSUM(I,J,K)) 
         TMP(I,J,K) = PBAR(K,PRESSURE_ZONE(I,J,K))/(RSUM(I,J,K)*RHO(I,J,K))
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE STRATIFIED_MIXING_LAYER


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
         IF (MMHAT(I,J,K)<=ZERO_P) THEN
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


SUBROUTINE WERNER_WENGLE_WALL_MODEL(SF,U_TAU,U1,NU,DZ,ROUGHNESS)

REAL(EB), INTENT(OUT) :: SF,U_TAU
REAL(EB), INTENT(IN) :: U1,NU,DZ,ROUGHNESS

REAL(EB), PARAMETER :: A=8.3_EB,B=1._EB/7._EB
REAL(EB), PARAMETER :: Z_PLUS_TURBULENT = 11.81_EB
REAL(EB), PARAMETER :: ALPHA=7.202125273562269_EB !! ALPHA=(1._EB-B)/2._EB*A**((1._EB+B)/(1._EB-B))
REAL(EB), PARAMETER :: BETA=1._EB+B
REAL(EB), PARAMETER :: ETA=(1._EB+B)/A
REAL(EB), PARAMETER :: GAMMA=2._EB/(1._EB+B)
REAL(EB), PARAMETER :: RKAPPA=2.44_EB  ! 1./von Karman constant
REAL(EB), PARAMETER :: B_LOGLAW=5.2_EB, B2=8.5_EB, BTILDE_MAX = 9.5_EB ! see Pope (2000) pp. 294,297,298
REAL(EB), PARAMETER :: R_PLUS_SMOOTH=5.83_EB, R_PLUS_ROUGH=30._EB ! approx piece-wise function for Fig. 7.24, Pope (2000) p. 297
REAL(EB), PARAMETER :: EPS=1.E-10_EB

REAL(EB) :: TAU_W,NUODZ,Z_PLUS,TAU_ROUGH,BTILDE,RD_NU,R_PLUS,TAU_SMOOTH
INTEGER :: ITER

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

! Werner-Wengle
NUODZ = NU/DZ
TAU_W = (ALPHA*(NUODZ**BETA) + ETA*(NUODZ**B)*ABS(U1))**GAMMA ! actually tau_w/rho
U_TAU = SQRT(TAU_W)
RD_NU  = U_TAU/(NU+EPS) ! viscous length scale

! Pope (2000)
IF (ROUGHNESS>0._EB) THEN
   TAU_SMOOTH=TAU_W
   DO ITER=1,2
      R_PLUS = ROUGHNESS*RD_NU ! roughness in viscous units
      IF (R_PLUS < R_PLUS_SMOOTH) THEN
         BTILDE = B_LOGLAW + RKAPPA*LOG(R_PLUS) ! Pope (2000) p. 297, Eq. (7.122)
      ELSEIF (R_PLUS < R_PLUS_ROUGH) THEN
         BTILDE = BTILDE_MAX ! approximation from Fig. 7.24, Pope (2000) p. 297
      ELSE
         BTILDE = B2 ! fully rough
      ENDIF
      TAU_ROUGH = ( U1/(RKAPPA*LOG(0.5_EB*DZ/ROUGHNESS)+BTILDE) )**2 ! actually tau_w/rho
      TAU_W = MAX(TAU_SMOOTH,TAU_ROUGH)
      U_TAU = SQRT(TAU_W)
      RD_NU  = U_TAU/(NU+EPS)
   ENDDO
ENDIF

Z_PLUS = DZ*RD_NU
IF (Z_PLUS>Z_PLUS_TURBULENT) THEN
   SF = 1._EB-TAU_W/(NUODZ*ABS(U1)+EPS) ! log layer
ELSE
   SF = -1._EB ! viscous sublayer
   TAU_W = 2._EB*NUODZ*ABS(U1)
   U_TAU = SQRT(TAU_W)
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
      CALL WERNER_WENGLE_WALL_MODEL(SLIP_COEF,DUMMY,U_STRM-U_STRM_WALL,MU*RRHO,DN,0._EB)
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
IF (ABS(NORM2(PP))<=ZERO_P) THEN
   ! tangent vector is completely arbitrary, just perpendicular to NN
   IF (ABS(NN(1))>=ZERO_P .OR.  ABS(NN(2))>=ZERO_P) PP = (/NN(2),-NN(1),0._EB/)
   IF (ABS(NN(1))<=ZERO_P .AND. ABS(NN(2))<=ZERO_P) PP = (/NN(3),0._EB,-NN(1)/)
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
   IF (ABS(AA)>=ZERO_P) THEN
      U_STRM = ((AA*U_STRM + BB)*EXP(AA*DT) - BB)/AA
   ELSE
      VELTAN3D = U_INT
      RETURN
   ENDIF
ELSE
   U_STRM_0 = U_STRM
   DO SUBIT=1,1
      CALL WERNER_WENGLE_WALL_MODEL(SLIP_COEF,DUMMY,U_STRM,MU*RRHO,DN,ROUGHNESS)
      !IF (SLIP_COEF<-100._EB .OR. SLIP_COEF>100._EB) THEN
      !   PRINT *,SUBIT,'WARNING: SLIP_COEF=',SLIP_COEF
      !ENDIF
      ETA = RRHO*(1-SLIP_COEF)*MU*0.5_EB*RDN**2
      AA  = -(0.5_EB*DUSDS + TWTH*U_NORM*RDN + ETA)
      BB  = -(U_NORM*ONTH*DUSDN + RRHO*( DPDS + TSN*0.5_EB*RDN))
      !print *,MU*RRHO*DT/(DN**2)
      IF (ABS(AA)>=ZERO_P) THEN
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

REAL(EB), INTENT(IN) :: DT,T
INTEGER, INTENT(IN) :: NM
INTEGER :: NE,NV,II,JJ,KK,IERROR
TYPE(VENTS_TYPE), POINTER :: VT=>NULL()
TYPE(SURFACE_TYPE), POINTER :: SF=>NULL()
REAL(EB) :: XX,YY,ZZ,SHAPE_FACTOR,VOLUME_WEIGHTING_FACTOR(3),EDDY_VOLUME(3),PROFILE_FACTOR,RAMP_T,TSI

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
            IF (SF%PROFILE==ATMOSPHERIC) PROFILE_FACTOR = (MAX(0._EB,VT%Z_EDDY(NE)-GROUND_LEVEL)/SF%Z0)**SF%PLE
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
            IF (SF%PROFILE==ATMOSPHERIC) PROFILE_FACTOR = (MAX(0._EB,VT%Z_EDDY(NE)-GROUND_LEVEL)/SF%Z0)**SF%PLE
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
            IF (SF%PROFILE==ATMOSPHERIC) PROFILE_FACTOR = (MAX(0._EB,VT%Z_EDDY(NE)-GROUND_LEVEL)/SF%Z0)**SF%PLE
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

USE PHYSICAL_FUNCTIONS, ONLY: LES_FILTER_WIDTH

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

         DELTA = LES_FILTER_WIDTH(DXN(I),DY(J),DZ(K))
               
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

         DELTA = LES_FILTER_WIDTH(DX(I),DYN(J),DZ(K))
               
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

         DELTA = LES_FILTER_WIDTH(DX(I),DY(J),DZN(K))
               
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

