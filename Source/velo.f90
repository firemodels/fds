MODULE VELO

! Module computes the velocity flux terms, baroclinic torque correction terms, and performs the CFL Check

USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_POINTERS
USE COMP_FUNCTIONS, ONLY: CURRENT_TIME

IMPLICIT NONE
PRIVATE

PUBLIC COMPUTE_VELOCITY_FLUX,VELOCITY_PREDICTOR,VELOCITY_CORRECTOR,NO_FLUX,BAROCLINIC_CORRECTION, &
       MATCH_VELOCITY,MATCH_VELOCITY_FLUX,VELOCITY_BC,COMPUTE_VISCOSITY,VISCOSITY_BC
PRIVATE VELOCITY_FLUX,VELOCITY_FLUX_CYLINDRICAL


CONTAINS


SUBROUTINE COMPUTE_VELOCITY_FLUX(T,DT,NM,FUNCTION_CODE)

REAL(EB), INTENT(IN) :: T,DT
REAL(EB) :: TNOW
INTEGER, INTENT(IN) :: NM,FUNCTION_CODE

IF (SOLID_PHASE_ONLY .OR. FREEZE_VELOCITY) RETURN

TNOW = CURRENT_TIME()

SELECT CASE(FUNCTION_CODE)
   CASE(1)
      CALL COMPUTE_VISCOSITY(T,NM)
   CASE(2)
      MESHES(NM)%BAROCLINIC_TERMS_ATTACHED = .FALSE.
      CALL VISCOSITY_BC(NM)
      IF (.NOT.CYLINDRICAL) CALL VELOCITY_FLUX(T,DT,NM)
      IF (     CYLINDRICAL) CALL VELOCITY_FLUX_CYLINDRICAL(T,NM)
END SELECT

T_USED(4) = T_USED(4) + CURRENT_TIME() - TNOW
END SUBROUTINE COMPUTE_VELOCITY_FLUX


SUBROUTINE COMPUTE_VISCOSITY(T,NM)

USE PHYSICAL_FUNCTIONS, ONLY: GET_VISCOSITY,LES_FILTER_WIDTH_FUNCTION,GET_POTENTIAL_TEMPERATURE,TURBULENT_VISCOSITY_INTERP1D
USE TURBULENCE, ONLY: VARDEN_DYNSMAG,TEST_FILTER,FILL_EDGES,WALL_MODEL,RNG_EDDY_VISCOSITY,WALE_VISCOSITY
USE MATH_FUNCTIONS, ONLY:EVALUATE_RAMP
USE COMPLEX_GEOMETRY, ONLY : CCREGION_COMPUTE_VISCOSITY
REAL(EB), INTENT(IN) :: T
INTEGER, INTENT(IN) :: NM
REAL(EB) :: ZZ_GET(1:N_TRACKED_SPECIES),NU_EDDY,DELTA,KSGS,U2,V2,W2,AA,A_IJ(3,3),BB,B_IJ(3,3),&
            DUDX,DUDY,DUDZ,DVDX,DVDY,DVDZ,DWDX,DWDY,DWDZ,MU_EFF,SLIP_COEF,VEL_GAS,VEL_T,RAMP_T,TSI,&
            VDF,LS,THETA_0,THETA_1,THETA_2,DTDZBAR,WGT,NU_2,NU_3,DX_1,DX_2,DX_3
REAL(EB), PARAMETER :: RAPLUS=1._EB/26._EB, C_LS=0.76_EB
INTEGER :: I,J,K,IIG,JJG,KKG,II,JJ,KK,IW,IOR
REAL(EB), POINTER, DIMENSION(:,:,:) :: RHOP=>NULL(),UP=>NULL(),VP=>NULL(),WP=>NULL(), &
                                       UP_HAT=>NULL(),VP_HAT=>NULL(),WP_HAT=>NULL(), &
                                       UU=>NULL(),VV=>NULL(),WW=>NULL(),DTDZ=>NULL()
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZP=>NULL()
INTEGER, POINTER, DIMENSION(:,:,:) :: CELL_COUNTER=>NULL()
TYPE(WALL_TYPE), POINTER :: WC=>NULL()
TYPE(SURFACE_TYPE), POINTER :: SF=>NULL()

IF (EVACUATION_ONLY(NM)) RETURN ! No need to update viscosity, use initial one

CALL POINT_TO_MESH(NM)

IF (PREDICTOR) THEN
   RHOP => RHO
   UU   => U
   VV   => V
   WW   => W
   ZZP  => ZZ
ELSE
   RHOP => RHOS
   UU   => US
   VV   => VS
   WW   => WS
   ZZP  => ZZS
ENDIF

! Compute viscosity for DNS using primitive species

IF (SIM_MODE==SVLES_MODE) THEN

   MU_DNS = MU_AIR_0

ELSE

   !$OMP PARALLEL DO FIRSTPRIVATE(ZZ_GET) SCHEDULE(guided)
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
            ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I,J,K,1:N_TRACKED_SPECIES)
            CALL GET_VISCOSITY(ZZ_GET,MU_DNS(I,J,K),TMP(I,J,K))
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO

ENDIF

CALL COMPUTE_STRAIN_RATE(NM)

SELECT_TURB: SELECT CASE (TURB_MODEL)

   CASE (NO_TURB_MODEL)

      MU = MU_DNS

   CASE (CONSMAG,DYNSMAG) SELECT_TURB ! Smagorinsky (1963) eddy viscosity

      IF (PREDICTOR .AND. TURB_MODEL==DYNSMAG) CALL VARDEN_DYNSMAG(NM) ! dynamic procedure, Moin et al. (1991)

      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
               MU(I,J,K) = MU_DNS(I,J,K) + RHOP(I,J,K)*CSD2(I,J,K)*STRAIN_RATE(I,J,K)
            ENDDO
         ENDDO
      ENDDO

   CASE (DEARDORFF) SELECT_TURB ! Deardorff (1980) eddy viscosity model (current default)

      ! Velocities relative to the p-cell center

      UP => WORK1
      VP => WORK2
      WP => WORK3
      UP=0._EB
      VP=0._EB
      WP=0._EB

      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
               UP(I,J,K) = 0.5_EB*(UU(I,J,K) + UU(I-1,J,K))
               VP(I,J,K) = 0.5_EB*(VV(I,J,K) + VV(I,J-1,K))
               WP(I,J,K) = 0.5_EB*(WW(I,J,K) + WW(I,J,K-1))
            ENDDO
         ENDDO
      ENDDO

      ! fill mesh boundary ghost cells

      DO IW=1,N_EXTERNAL_WALL_CELLS
         WC=>WALL(IW)
         SELECT CASE(WC%BOUNDARY_TYPE)
            CASE(INTERPOLATED_BOUNDARY)
               II = WC%ONE_D%II
               JJ = WC%ONE_D%JJ
               KK = WC%ONE_D%KK
               UP(II,JJ,KK) = U_GHOST(IW)
               VP(II,JJ,KK) = V_GHOST(IW)
               WP(II,JJ,KK) = W_GHOST(IW)
            CASE(OPEN_BOUNDARY,MIRROR_BOUNDARY)
               II = WC%ONE_D%II
               JJ = WC%ONE_D%JJ
               KK = WC%ONE_D%KK
               IIG = WC%ONE_D%IIG
               JJG = WC%ONE_D%JJG
               KKG = WC%ONE_D%KKG
               UP(II,JJ,KK) = UP(IIG,JJG,KKG)
               VP(II,JJ,KK) = VP(IIG,JJG,KKG)
               WP(II,JJ,KK) = WP(IIG,JJG,KKG)
         END SELECT
      ENDDO

      ! fill edge and corner ghost cells

      CALL FILL_EDGES(UP)
      CALL FILL_EDGES(VP)
      CALL FILL_EDGES(WP)

      UP_HAT => WORK4
      VP_HAT => WORK5
      WP_HAT => WORK6
      UP_HAT=0._EB
      VP_HAT=0._EB
      WP_HAT=0._EB

      CALL TEST_FILTER(UP_HAT,UP)
      CALL TEST_FILTER(VP_HAT,VP)
      CALL TEST_FILTER(WP_HAT,WP)

      POTENTIAL_TEMPERATURE_IF: IF (.NOT.POTENTIAL_TEMPERATURE_CORRECTION) THEN
         !$OMP PARALLEL DO PRIVATE(DELTA, KSGS, NU_EDDY) SCHEDULE(static)
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
                  DELTA = LES_FILTER_WIDTH_FUNCTION(DX(I),DY(J),DZ(K))
                  KSGS = 0.5_EB*( (UP(I,J,K)-UP_HAT(I,J,K))**2 + (VP(I,J,K)-VP_HAT(I,J,K))**2 + (WP(I,J,K)-WP_HAT(I,J,K))**2 )

                  NU_EDDY = C_DEARDORFF*DELTA*SQRT(KSGS)
                  MU(I,J,K) = MU_DNS(I,J,K) + RHOP(I,J,K)*NU_EDDY
               ENDDO
            ENDDO
         ENDDO
         !$OMP END PARALLEL DO
      ELSE POTENTIAL_TEMPERATURE_IF
         DTDZ => WORK7
         DO K=0,KBAR
            DO J=0,JBAR
               DO I=0,IBAR
                  THETA_1 = GET_POTENTIAL_TEMPERATURE(TMP(I,J,K),ZC(K))
                  THETA_2 = GET_POTENTIAL_TEMPERATURE(TMP(I,J,K+1),ZC(K+1))
                  DTDZ(I,J,K) = (THETA_2-THETA_1)*RDZN(K)
               ENDDO
            ENDDO
         ENDDO
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
                  DELTA = LES_FILTER_WIDTH_FUNCTION(DX(I),DY(J),DZ(K))
                  LS = DELTA
                  KSGS = 0.5_EB*( (UP(I,J,K)-UP_HAT(I,J,K))**2 + (VP(I,J,K)-VP_HAT(I,J,K))**2 + (WP(I,J,K)-WP_HAT(I,J,K))**2 )
                  DTDZBAR = 0.5_EB*(DTDZ(I,J,K)+DTDZ(I,J,K+1))
                  IF (DTDZBAR>0._EB) THEN
                     THETA_0 = GET_POTENTIAL_TEMPERATURE(TMP_0(K),ZC(K))
                     LS = C_LS*SQRT(KSGS)/SQRT(ABS(GVEC(3))/THETA_0*DTDZBAR) ! von Schoenberg Eq. (3.19)
                  ENDIF
                  NU_EDDY = C_DEARDORFF*MIN(LS,DELTA)*SQRT(KSGS)
                  MU(I,J,K) = MU_DNS(I,J,K) + RHOP(I,J,K)*NU_EDDY
                  PR_T(I,J,K) = 1._EB/(1._EB + (2._EB*MIN(LS,DELTA)/DELTA)) ! von Schoenberg Eq. (3.21)
               ENDDO
            ENDDO
         ENDDO
      ENDIF POTENTIAL_TEMPERATURE_IF

   CASE (VREMAN) SELECT_TURB ! Vreman (2004) eddy viscosity model (experimental)

      ! A. W. Vreman. An eddy-viscosity subgrid-scale model for turbulent shear flow: Algebraic theory and applications.
      ! Phys. Fluids, 16(10):3670-3681, 2004.

      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
               DUDX = RDX(I)*(UU(I,J,K)-UU(I-1,J,K))
               DVDY = RDY(J)*(VV(I,J,K)-VV(I,J-1,K))
               DWDZ = RDZ(K)*(WW(I,J,K)-WW(I,J,K-1))
               DUDY = 0.25_EB*RDY(J)*(UU(I,J+1,K)-UU(I,J-1,K)+UU(I-1,J+1,K)-UU(I-1,J-1,K))
               DUDZ = 0.25_EB*RDZ(K)*(UU(I,J,K+1)-UU(I,J,K-1)+UU(I-1,J,K+1)-UU(I-1,J,K-1))
               DVDX = 0.25_EB*RDX(I)*(VV(I+1,J,K)-VV(I-1,J,K)+VV(I+1,J-1,K)-VV(I-1,J-1,K))
               DVDZ = 0.25_EB*RDZ(K)*(VV(I,J,K+1)-VV(I,J,K-1)+VV(I,J-1,K+1)-VV(I,J-1,K-1))
               DWDX = 0.25_EB*RDX(I)*(WW(I+1,J,K)-WW(I-1,J,K)+WW(I+1,J,K-1)-WW(I-1,J,K-1))
               DWDY = 0.25_EB*RDY(J)*(WW(I,J+1,K)-WW(I,J-1,K)+WW(I,J+1,K-1)-WW(I,J-1,K-1))

               ! Vreman, Eq. (6)
               A_IJ(1,1)=DUDX; A_IJ(2,1)=DUDY; A_IJ(3,1)=DUDZ
               A_IJ(1,2)=DVDX; A_IJ(2,2)=DVDY; A_IJ(3,2)=DVDZ
               A_IJ(1,3)=DWDX; A_IJ(2,3)=DWDY; A_IJ(3,3)=DWDZ

               AA=0._EB
               DO JJ=1,3
                  DO II=1,3
                     AA = AA + A_IJ(II,JJ)*A_IJ(II,JJ)
                  ENDDO
               ENDDO

               ! Vreman, Eq. (7)
               B_IJ(1,1)=(DX(I)*A_IJ(1,1))**2 + (DY(J)*A_IJ(2,1))**2 + (DZ(K)*A_IJ(3,1))**2
               B_IJ(2,2)=(DX(I)*A_IJ(1,2))**2 + (DY(J)*A_IJ(2,2))**2 + (DZ(K)*A_IJ(3,2))**2
               B_IJ(3,3)=(DX(I)*A_IJ(1,3))**2 + (DY(J)*A_IJ(2,3))**2 + (DZ(K)*A_IJ(3,3))**2

               B_IJ(1,2)=DX(I)**2*A_IJ(1,1)*A_IJ(1,2) + DY(J)**2*A_IJ(2,1)*A_IJ(2,2) + DZ(K)**2*A_IJ(3,1)*A_IJ(3,2)
               B_IJ(1,3)=DX(I)**2*A_IJ(1,1)*A_IJ(1,3) + DY(J)**2*A_IJ(2,1)*A_IJ(2,3) + DZ(K)**2*A_IJ(3,1)*A_IJ(3,3)
               B_IJ(2,3)=DX(I)**2*A_IJ(1,2)*A_IJ(1,3) + DY(J)**2*A_IJ(2,2)*A_IJ(2,3) + DZ(K)**2*A_IJ(3,2)*A_IJ(3,3)

               BB = B_IJ(1,1)*B_IJ(2,2) - B_IJ(1,2)**2 &
                  + B_IJ(1,1)*B_IJ(3,3) - B_IJ(1,3)**2 &
                  + B_IJ(2,2)*B_IJ(3,3) - B_IJ(2,3)**2    ! Vreman, Eq. (8)

               IF (ABS(AA)>TWO_EPSILON_EB .AND. BB>TWO_EPSILON_EB) THEN
                  NU_EDDY = C_VREMAN*SQRT(BB/AA)  ! Vreman, Eq. (5)
               ELSE
                  NU_EDDY=0._EB
               ENDIF

               MU(I,J,K) = MU_DNS(I,J,K) + RHOP(I,J,K)*NU_EDDY

            ENDDO
         ENDDO
      ENDDO

   CASE (RNG) SELECT_TURB

      ! A. Yakhot, S. A. Orszag, V. Yakhot, and M. Israeli. Renormalization Group Formulation of Large-Eddy Simulation.
      ! Journal of Scientific Computing, 1(1):1-51, 1989.

      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
               DELTA = LES_FILTER_WIDTH_FUNCTION(DX(I),DY(J),DZ(K))
               CALL RNG_EDDY_VISCOSITY(MU_EFF,MU_DNS(I,J,K),RHOP(I,J,K),STRAIN_RATE(I,J,K),DELTA)
               MU(I,J,K) = MU_EFF
            ENDDO
         ENDDO
      ENDDO

   CASE (WALE) SELECT_TURB

      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
               DELTA = LES_FILTER_WIDTH_FUNCTION(DX(I),DY(J),DZ(K))
               ! compute velocity gradient tensor
               DUDX = RDX(I)*(UU(I,J,K)-UU(I-1,J,K))
               DVDY = RDY(J)*(VV(I,J,K)-VV(I,J-1,K))
               DWDZ = RDZ(K)*(WW(I,J,K)-WW(I,J,K-1))
               DUDY = 0.25_EB*RDY(J)*(UU(I,J+1,K)-UU(I,J-1,K)+UU(I-1,J+1,K)-UU(I-1,J-1,K))
               DUDZ = 0.25_EB*RDZ(K)*(UU(I,J,K+1)-UU(I,J,K-1)+UU(I-1,J,K+1)-UU(I-1,J,K-1))
               DVDX = 0.25_EB*RDX(I)*(VV(I+1,J,K)-VV(I-1,J,K)+VV(I+1,J-1,K)-VV(I-1,J-1,K))
               DVDZ = 0.25_EB*RDZ(K)*(VV(I,J,K+1)-VV(I,J,K-1)+VV(I,J-1,K+1)-VV(I,J-1,K-1))
               DWDX = 0.25_EB*RDX(I)*(WW(I+1,J,K)-WW(I-1,J,K)+WW(I+1,J,K-1)-WW(I-1,J,K-1))
               DWDY = 0.25_EB*RDY(J)*(WW(I,J+1,K)-WW(I,J-1,K)+WW(I,J+1,K-1)-WW(I,J-1,K-1))
               A_IJ(1,1)=DUDX; A_IJ(1,2)=DUDY; A_IJ(1,3)=DUDZ
               A_IJ(2,1)=DVDX; A_IJ(2,2)=DVDY; A_IJ(2,3)=DVDZ
               A_IJ(3,1)=DWDX; A_IJ(3,2)=DWDY; A_IJ(3,3)=DWDZ

               CALL WALE_VISCOSITY(NU_EDDY,A_IJ,DELTA)

               MU(I,J,K) = MU_DNS(I,J,K) + RHOP(I,J,K)*NU_EDDY
            ENDDO
         ENDDO
      ENDDO

END SELECT SELECT_TURB

! Compute resolved kinetic energy per unit mass

DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
         U2 = 0.25_EB*(UU(I-1,J,K)+UU(I,J,K))**2
         V2 = 0.25_EB*(VV(I,J-1,K)+VV(I,J,K))**2
         W2 = 0.25_EB*(WW(I,J,K-1)+WW(I,J,K))**2
         KRES(I,J,K) = 0.5_EB*(U2+V2+W2)
      ENDDO
   ENDDO
ENDDO

! Mirror viscosity into solids and exterior boundary cells

CELL_COUNTER => IWORK1 ; CELL_COUNTER = 0

WALL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS

   WC=>WALL(IW)
   IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE WALL_LOOP
   II  = WC%ONE_D%II
   JJ  = WC%ONE_D%JJ
   KK  = WC%ONE_D%KK
   IOR = WC%ONE_D%IOR
   IIG = WC%ONE_D%IIG
   JJG = WC%ONE_D%JJG
   KKG = WC%ONE_D%KKG
   SF=>SURFACE(WC%SURF_INDEX)

   SELECT CASE(WC%BOUNDARY_TYPE)

      CASE(SOLID_BOUNDARY)

         IF (ABS(SF%T_IGN-T_BEGIN)<=SPACING(SF%T_IGN) .AND. SF%RAMP_INDEX(TIME_VELO)>=1) THEN
            TSI = T
         ELSE
            TSI = T-SF%T_IGN
         ENDIF
         RAMP_T = EVALUATE_RAMP(TSI,SF%TAU(TIME_VELO),SF%RAMP_INDEX(TIME_VELO))
         VEL_T = RAMP_T*SQRT(SF%VEL_T(1)**2 + SF%VEL_T(2)**2)

         SELECT CASE(ABS(IOR))
            CASE(1)
               VEL_GAS = SQRT( 0.25_EB*( (VV(IIG,JJG,KKG)+VV(IIG,JJG-1,KKG))**2 + (WW(IIG,JJG,KKG)+WW(IIG,JJG,KKG-1))**2 ) )
            CASE(2)
               VEL_GAS = SQRT( 0.25_EB*( (UU(IIG,JJG,KKG)+UU(IIG-1,JJG,KKG))**2 + (WW(IIG,JJG,KKG)+WW(IIG,JJG,KKG-1))**2 ) )
            CASE(3)
               VEL_GAS = SQRT( 0.25_EB*( (UU(IIG,JJG,KKG)+UU(IIG-1,JJG,KKG))**2 + (VV(IIG,JJG,KKG)+VV(IIG,JJG-1,KKG))**2 ) )
         END SELECT

         CALL WALL_MODEL(SLIP_COEF,WC%ONE_D%U_TAU,WC%ONE_D%Y_PLUS,VEL_GAS-VEL_T,&
                         MU_DNS(IIG,JJG,KKG)/RHO(IIG,JJG,KKG),1._EB/WC%ONE_D%RDN,SURFACE(WC%SURF_INDEX)%ROUGHNESS)

         IF (SIM_MODE/=DNS_MODE) THEN
            DELTA = LES_FILTER_WIDTH_FUNCTION(DX(IIG),DY(JJG),DZ(KKG))
            SELECT CASE(NEAR_WALL_TURB_MODEL)
               CASE DEFAULT ! Constant Smagorinsky with Van Driest damping
                  VDF = 1._EB-EXP(-WC%ONE_D%Y_PLUS*RAPLUS)
                  NU_EDDY = (VDF*C_SMAGORINSKY*DELTA)**2*STRAIN_RATE(IIG,JJG,KKG)
               CASE(WALE)
                  ! compute velocity gradient tensor
                  DUDX = RDX(IIG)*(UU(IIG,JJG,KKG)-UU(IIG-1,JJG,KKG))
                  DVDY = RDY(JJG)*(VV(IIG,JJG,KKG)-VV(IIG,JJG-1,KKG))
                  DWDZ = RDZ(KKG)*(WW(IIG,JJG,KKG)-WW(IIG,JJG,KKG-1))
                  DUDY = 0.25_EB*RDY(JJG)*(UU(IIG,JJG+1,KKG)-UU(IIG,JJG-1,KKG)+UU(IIG-1,JJG+1,KKG)-UU(IIG-1,JJG-1,KKG))
                  DUDZ = 0.25_EB*RDZ(KKG)*(UU(IIG,JJG,KKG+1)-UU(IIG,JJG,KKG-1)+UU(IIG-1,JJG,KKG+1)-UU(IIG-1,JJG,KKG-1))
                  DVDX = 0.25_EB*RDX(IIG)*(VV(IIG+1,JJG,KKG)-VV(IIG-1,JJG,KKG)+VV(IIG+1,JJG-1,KKG)-VV(IIG-1,JJG-1,KKG))
                  DVDZ = 0.25_EB*RDZ(KKG)*(VV(IIG,JJG,KKG+1)-VV(IIG,JJG,KKG-1)+VV(IIG,JJG-1,KKG+1)-VV(IIG,JJG-1,KKG-1))
                  DWDX = 0.25_EB*RDX(IIG)*(WW(IIG+1,JJG,KKG)-WW(IIG-1,JJG,KKG)+WW(IIG+1,JJG,KKG-1)-WW(IIG-1,JJG,KKG-1))
                  DWDY = 0.25_EB*RDY(JJG)*(WW(IIG,JJG+1,KKG)-WW(IIG,JJG-1,KKG)+WW(IIG,JJG+1,KKG-1)-WW(IIG,JJG-1,KKG-1))
                  A_IJ(1,1)=DUDX; A_IJ(1,2)=DUDY; A_IJ(1,3)=DUDZ
                  A_IJ(2,1)=DVDX; A_IJ(2,2)=DVDY; A_IJ(2,3)=DVDZ
                  A_IJ(3,1)=DWDX; A_IJ(3,2)=DWDY; A_IJ(3,3)=DWDZ
                  CALL WALE_VISCOSITY(NU_EDDY,A_IJ,DELTA)
               CASE(MU_TURB_INTERP)
                  !! Experimental !!
                  ! avoids jump in viscosity model from bulk to near-wall cells
                  DX_1 = 1._EB/WC%ONE_D%RDN
                  SELECT CASE(IOR)
                     CASE( 1)
                        NU_2 = RHO(MIN(IIG+1,IBP1),JJG,KKG)*MU(MIN(IIG+1,IBP1),JJG,KKG)
                        NU_3 = RHO(MIN(IIG+2,IBP1),JJG,KKG)*MU(MIN(IIG+2,IBP1),JJG,KKG)
                        DX_2 = DX(MIN(IIG+1,IBP1))
                        DX_3 = DX(MIN(IIG+2,IBP1))
                     CASE(-1)
                        NU_2 = RHO(MAX(IIG-1,0),JJG,KKG)*MU(MAX(IIG-1,0),JJG,KKG)
                        NU_3 = RHO(MAX(IIG-2,0),JJG,KKG)*MU(MAX(IIG-2,0),JJG,KKG)
                        DX_2 = DX(MAX(IIG-1,0))
                        DX_3 = DX(MAX(IIG-2,0))
                     CASE( 2)
                        NU_2 = RHO(IIG,MIN(JJG+1,JBP1),KKG)*MU(IIG,MIN(JJG+1,JBP1),KKG)
                        NU_3 = RHO(IIG,MIN(JJG+2,JBP1),KKG)*MU(IIG,MIN(JJG+2,JBP1),KKG)
                        DX_2 = DY(MIN(JJG+1,JBP1))
                        DX_3 = DY(MIN(JJG+2,JBP1))
                     CASE(-2)
                        NU_2 = RHO(IIG,MAX(JJG-1,0),KKG)*MU(IIG,MAX(JJG-1,0),KKG)
                        NU_3 = RHO(IIG,MAX(JJG-2,0),KKG)*MU(IIG,MAX(JJG-2,0),KKG)
                        DX_2 = DY(MAX(JJG-1,0))
                        DX_3 = DY(MAX(JJG-2,0))
                     CASE( 3)
                        NU_2 = RHO(IIG,JJG,MIN(KKG+1,KBP1))*MU(IIG,JJG,MIN(KKG+1,KBP1))
                        NU_3 = RHO(IIG,JJG,MIN(KKG+2,KBP1))*MU(IIG,JJG,MIN(KKG+2,KBP1))
                        DX_2 = DZ(MIN(KKG+1,KBP1))
                        DX_3 = DZ(MIN(KKG+2,KBP1))
                     CASE(-3)
                        NU_2 = RHO(IIG,JJG,MAX(KKG-1,0))*MU(IIG,JJG,MAX(KKG-1,0))
                        NU_3 = RHO(IIG,JJG,MAX(KKG-2,0))*MU(IIG,JJG,MAX(KKG-2,0))
                        DX_2 = DZ(MAX(KKG-1,0))
                        DX_3 = DZ(MAX(KKG-2,0))
                  END SELECT
                  CALL TURBULENT_VISCOSITY_INTERP1D(NU_EDDY,NU_2,NU_3,DX_1,DX_2,DX_3)
            END SELECT
            IF (CELL_COUNTER(IIG,JJG,KKG)==0) MU(IIG,JJG,KKG) = 0._EB
            CELL_COUNTER(IIG,JJG,KKG) = CELL_COUNTER(IIG,JJG,KKG) + 1
            WGT = 1._EB/REAL(CELL_COUNTER(IIG,JJG,KKG),EB)
            MU(IIG,JJG,KKG) = (1._EB-WGT)*MU(IIG,JJG,KKG) + WGT*(MU_DNS(IIG,JJG,KKG) + RHOP(IIG,JJG,KKG)*NU_EDDY)
         ELSE
            MU(IIG,JJG,KKG) = MU_DNS(IIG,JJG,KKG)
         ENDIF

         IF (SOLID(CELL_INDEX(II,JJ,KK))) MU(II,JJ,KK) = MU(IIG,JJG,KKG)

      CASE(OPEN_BOUNDARY,MIRROR_BOUNDARY)

         MU(II,JJ,KK) = MU(IIG,JJG,KKG)
         KRES(II,JJ,KK) = KRES(IIG,JJG,KKG)

   END SELECT

ENDDO WALL_LOOP

IF(CC_IBM) CALL CCREGION_COMPUTE_VISCOSITY(0._EB,NM)

MU(   0,0:JBP1,   0) = MU(   1,0:JBP1,1)
MU(IBP1,0:JBP1,   0) = MU(IBAR,0:JBP1,1)
MU(IBP1,0:JBP1,KBP1) = MU(IBAR,0:JBP1,KBAR)
MU(   0,0:JBP1,KBP1) = MU(   1,0:JBP1,KBAR)
MU(0:IBP1,   0,   0) = MU(0:IBP1,   1,1)
MU(0:IBP1,JBP1,0)    = MU(0:IBP1,JBAR,1)
MU(0:IBP1,JBP1,KBP1) = MU(0:IBP1,JBAR,KBAR)
MU(0:IBP1,0,KBP1)    = MU(0:IBP1,   1,KBAR)
MU(0,   0,0:KBP1)    = MU(   1,   1,0:KBP1)
MU(IBP1,0,0:KBP1)    = MU(IBAR,   1,0:KBP1)
MU(IBP1,JBP1,0:KBP1) = MU(IBAR,JBAR,0:KBP1)
MU(0,JBP1,0:KBP1)    = MU(   1,JBAR,0:KBP1)

END SUBROUTINE COMPUTE_VISCOSITY


SUBROUTINE COMPUTE_STRAIN_RATE(NM)

INTEGER, INTENT(IN) :: NM
REAL(EB) :: DUDX,DUDY,DUDZ,DVDX,DVDY,DVDZ,DWDX,DWDY,DWDZ,S11,S22,S33,S12,S13,S23,ONTHDIV
INTEGER :: I,J,K,IOR,IIG,JJG,KKG,IW,SURF_INDEX
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU=>NULL(),VV=>NULL(),WW=>NULL()
TYPE(WALL_TYPE), POINTER :: WC=>NULL()

CALL POINT_TO_MESH(NM)

IF (PREDICTOR) THEN
   UU => U
   VV => V
   WW => W
ELSE
   UU => US
   VV => VS
   WW => WS
ENDIF

SELECT CASE (TURB_MODEL)
   CASE DEFAULT
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
               DUDX = RDX(I)*(UU(I,J,K)-UU(I-1,J,K))
               DVDY = RDY(J)*(VV(I,J,K)-VV(I,J-1,K))
               DWDZ = RDZ(K)*(WW(I,J,K)-WW(I,J,K-1))
               DUDY = 0.25_EB*RDY(J)*(UU(I,J+1,K)-UU(I,J-1,K)+UU(I-1,J+1,K)-UU(I-1,J-1,K))
               DUDZ = 0.25_EB*RDZ(K)*(UU(I,J,K+1)-UU(I,J,K-1)+UU(I-1,J,K+1)-UU(I-1,J,K-1))
               DVDX = 0.25_EB*RDX(I)*(VV(I+1,J,K)-VV(I-1,J,K)+VV(I+1,J-1,K)-VV(I-1,J-1,K))
               DVDZ = 0.25_EB*RDZ(K)*(VV(I,J,K+1)-VV(I,J,K-1)+VV(I,J-1,K+1)-VV(I,J-1,K-1))
               DWDX = 0.25_EB*RDX(I)*(WW(I+1,J,K)-WW(I-1,J,K)+WW(I+1,J,K-1)-WW(I-1,J,K-1))
               DWDY = 0.25_EB*RDY(J)*(WW(I,J+1,K)-WW(I,J-1,K)+WW(I,J+1,K-1)-WW(I,J-1,K-1))
               ONTHDIV = ONTH*(DUDX+DVDY+DWDZ)
               S11 = DUDX - ONTHDIV
               S22 = DVDY - ONTHDIV
               S33 = DWDZ - ONTHDIV
               S12 = 0.5_EB*(DUDY+DVDX)
               S13 = 0.5_EB*(DUDZ+DWDX)
               S23 = 0.5_EB*(DVDZ+DWDY)
               STRAIN_RATE(I,J,K) = SQRT(2._EB*(S11**2 + S22**2 + S33**2 + 2._EB*(S12**2 + S13**2 + S23**2)))
            ENDDO
         ENDDO
      ENDDO
   CASE (DEARDORFF)
      ! Here we omit the 3D loop, we only need the wall cell values of STRAIN_RATE
END SELECT

WALL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
   WC=>WALL(IW)
   IF (WC%BOUNDARY_TYPE/=SOLID_BOUNDARY) CYCLE WALL_LOOP

   SURF_INDEX = WC%SURF_INDEX
   IIG = WC%ONE_D%IIG
   JJG = WC%ONE_D%JJG
   KKG = WC%ONE_D%KKG
   IOR = WC%ONE_D%IOR

   ! Handle the case where OBST lives on an external boundary
   IF (IW>N_EXTERNAL_WALL_CELLS) THEN
      SELECT CASE(IOR)
         CASE( 1); IF (IIG>IBAR) CYCLE WALL_LOOP
         CASE(-1); IF (IIG<1)    CYCLE WALL_LOOP
         CASE( 2); IF (JJG>JBAR) CYCLE WALL_LOOP
         CASE(-2); IF (JJG<1)    CYCLE WALL_LOOP
         CASE( 3); IF (KKG>KBAR) CYCLE WALL_LOOP
         CASE(-3); IF (KKG<1)    CYCLE WALL_LOOP
      END SELECT
   ENDIF

   DUDX = RDX(IIG)*(UU(IIG,JJG,KKG)-UU(IIG-1,JJG,KKG))
   DVDY = RDY(JJG)*(VV(IIG,JJG,KKG)-VV(IIG,JJG-1,KKG))
   DWDZ = RDZ(KKG)*(WW(IIG,JJG,KKG)-WW(IIG,JJG,KKG-1))
   ONTHDIV = ONTH*(DUDX+DVDY+DWDZ)
   S11 = DUDX - ONTHDIV
   S22 = DVDY - ONTHDIV
   S33 = DWDZ - ONTHDIV

   DUDY = 0.25_EB*RDY(JJG)*(UU(IIG,JJG+1,KKG)-UU(IIG,JJG-1,KKG)+UU(IIG-1,JJG+1,KKG)-UU(IIG-1,JJG-1,KKG))
   DUDZ = 0.25_EB*RDZ(KKG)*(UU(IIG,JJG,KKG+1)-UU(IIG,JJG,KKG-1)+UU(IIG-1,JJG,KKG+1)-UU(IIG-1,JJG,KKG-1))
   DVDX = 0.25_EB*RDX(IIG)*(VV(IIG+1,JJG,KKG)-VV(IIG-1,JJG,KKG)+VV(IIG+1,JJG-1,KKG)-VV(IIG-1,JJG-1,KKG))
   DVDZ = 0.25_EB*RDZ(KKG)*(VV(IIG,JJG,KKG+1)-VV(IIG,JJG,KKG-1)+VV(IIG,JJG-1,KKG+1)-VV(IIG,JJG-1,KKG-1))
   DWDX = 0.25_EB*RDX(IIG)*(WW(IIG+1,JJG,KKG)-WW(IIG-1,JJG,KKG)+WW(IIG+1,JJG,KKG-1)-WW(IIG-1,JJG,KKG-1))
   DWDY = 0.25_EB*RDY(JJG)*(WW(IIG,JJG+1,KKG)-WW(IIG,JJG-1,KKG)+WW(IIG,JJG+1,KKG-1)-WW(IIG,JJG-1,KKG-1))

   S12 = 0.5_EB*(DUDY+DVDX)
   S13 = 0.5_EB*(DUDZ+DWDX)
   S23 = 0.5_EB*(DVDZ+DWDY)

   STRAIN_RATE(IIG,JJG,KKG) = SQRT(2._EB*(S11**2 + S22**2 + S33**2 + 2._EB*(S12**2 + S13**2 + S23**2)))
ENDDO WALL_LOOP

END SUBROUTINE COMPUTE_STRAIN_RATE


SUBROUTINE VISCOSITY_BC(NM)

! Specify ghost cell values of the viscosity array MU

INTEGER, INTENT(IN) :: NM
REAL(EB) :: MU_OTHER,DP_OTHER,KRES_OTHER
INTEGER :: II,JJ,KK,IW,IIO,JJO,KKO,NOM,N_INT_CELLS
TYPE(WALL_TYPE),POINTER :: WC=>NULL()
TYPE(EXTERNAL_WALL_TYPE),POINTER :: EWC=>NULL()

CALL POINT_TO_MESH(NM)

! Mirror viscosity into solids and exterior boundary cells

WALL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS
   WC =>WALL(IW)
   EWC=>EXTERNAL_WALL(IW)
   IF (EWC%NOM==0) CYCLE WALL_LOOP
   II  = WC%ONE_D%II
   JJ  = WC%ONE_D%JJ
   KK  = WC%ONE_D%KK
   NOM = EWC%NOM
   MU_OTHER   = 0._EB
   DP_OTHER   = 0._EB
   KRES_OTHER = 0._EB
   DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
      DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
         DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
            MU_OTHER = MU_OTHER + OMESH(NOM)%MU(IIO,JJO,KKO)
            KRES_OTHER = KRES_OTHER + OMESH(NOM)%KRES(IIO,JJO,KKO)
            IF (PREDICTOR) THEN
               DP_OTHER = DP_OTHER + OMESH(NOM)%D(IIO,JJO,KKO)
            ELSE
               DP_OTHER = DP_OTHER + OMESH(NOM)%DS(IIO,JJO,KKO)
            ENDIF
         ENDDO
      ENDDO
   ENDDO
   N_INT_CELLS = (EWC%IIO_MAX-EWC%IIO_MIN+1) * (EWC%JJO_MAX-EWC%JJO_MIN+1) * (EWC%KKO_MAX-EWC%KKO_MIN+1)
   MU_OTHER = MU_OTHER/REAL(N_INT_CELLS,EB)
   KRES_OTHER = KRES_OTHER/REAL(N_INT_CELLS,EB)
   DP_OTHER = DP_OTHER/REAL(N_INT_CELLS,EB)
   MU(II,JJ,KK) = MU_OTHER
   KRES(II,JJ,KK) = KRES_OTHER
   IF (PREDICTOR) THEN
      D(II,JJ,KK) = DP_OTHER
   ELSE
      DS(II,JJ,KK) = DP_OTHER
   ENDIF
ENDDO WALL_LOOP

END SUBROUTINE VISCOSITY_BC


SUBROUTINE VELOCITY_FLUX(T,DT,NM)

! Compute convective and diffusive terms of the momentum equations

USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
USE COMPLEX_GEOMETRY, ONLY : ROTATED_CUBE_VELOCITY_FLUX,CCIBM_INTERP_FACE_VEL
INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: T,DT
REAL(EB) :: MUX,MUY,MUZ,UP,UM,VP,VM,WP,WM,VTRM,OMXP,OMXM,OMYP,OMYM,OMZP,OMZM,TXYP,TXYM,TXZP,TXZM,TYZP,TYZM, &
            DTXYDY,DTXZDZ,DTYZDZ,DTXYDX,DTXZDX,DTYZDY, &
            DUDX,DVDY,DWDZ,DUDY,DUDZ,DVDX,DVDZ,DWDX,DWDY, &
            VOMZ,WOMY,UOMY,VOMX,UOMZ,WOMX, &
            RRHO,GX(0:IBAR_MAX),GY(0:IBAR_MAX),GZ(0:IBAR_MAX),TXXP,TXXM,TYYP,TYYM,TZZP,TZZM,DTXXDX,DTYYDY,DTZZDZ, &
            DUMMY=0._EB
INTEGER :: I,J,K,IEXP,IEXM,IEYP,IEYM,IEZP,IEZM,IC,IC1,IC2
REAL(EB), POINTER, DIMENSION(:,:,:) :: TXY=>NULL(),TXZ=>NULL(),TYZ=>NULL(),OMX=>NULL(),OMY=>NULL(),OMZ=>NULL(), &
                                       UU=>NULL(),VV=>NULL(),WW=>NULL(),RHOP=>NULL(),DP=>NULL()

CALL POINT_TO_MESH(NM)

IF (PREDICTOR) THEN
   UU => U
   VV => V
   WW => W
   DP => D
   RHOP => RHO
ELSE
   UU => US
   VV => VS
   WW => WS
   DP => DS
   RHOP => RHOS
ENDIF

TXY => WORK1
TXZ => WORK2
TYZ => WORK3
OMX => WORK4
OMY => WORK5
OMZ => WORK6

! Define velocities on gas cut-faces underlaying Cartesian faces.
IF (CC_IBM) CALL CCIBM_INTERP_FACE_VEL(DT,NM,.TRUE.)

! Compute vorticity and stress tensor components

!$OMP PARALLEL DO PRIVATE(DUDY, DVDX, DUDZ, DWDX, DVDZ, DWDY, &
!$OMP& MUX, MUY, MUZ) SCHEDULE(STATIC)
DO K=0,KBAR
   DO J=0,JBAR
      DO I=0,IBAR
         DUDY = RDYN(J)*(UU(I,J+1,K)-UU(I,J,K))
         DVDX = RDXN(I)*(VV(I+1,J,K)-VV(I,J,K))
         DUDZ = RDZN(K)*(UU(I,J,K+1)-UU(I,J,K))
         DWDX = RDXN(I)*(WW(I+1,J,K)-WW(I,J,K))
         DVDZ = RDZN(K)*(VV(I,J,K+1)-VV(I,J,K))
         DWDY = RDYN(J)*(WW(I,J+1,K)-WW(I,J,K))
         OMX(I,J,K) = DWDY - DVDZ
         OMY(I,J,K) = DUDZ - DWDX
         OMZ(I,J,K) = DVDX - DUDY
         MUX = 0.25_EB*(MU(I,J+1,K)+MU(I,J,K)+MU(I,J,K+1)+MU(I,J+1,K+1))
         MUY = 0.25_EB*(MU(I+1,J,K)+MU(I,J,K)+MU(I,J,K+1)+MU(I+1,J,K+1))
         MUZ = 0.25_EB*(MU(I+1,J,K)+MU(I,J,K)+MU(I,J+1,K)+MU(I+1,J+1,K))
         TXY(I,J,K) = MUZ*(DVDX + DUDY)
         TXZ(I,J,K) = MUY*(DUDZ + DWDX)
         TYZ(I,J,K) = MUX*(DVDZ + DWDY)
      ENDDO
   ENDDO
ENDDO
!$OMP END PARALLEL DO

! Wannier Flow (Stokes flow) test case

IF (PERIODIC_TEST==5) THEN
   OMX=0._EB
   OMY=0._EB
   OMZ=0._EB
   OME_E = -1.E6_EB
ENDIF

! Compute gravity components

IF (.NOT.SPATIAL_GRAVITY_VARIATION) THEN
   GX(0:IBAR) = EVALUATE_RAMP(T,DUMMY,I_RAMP_GX)*GVEC(1)
   GY(0:IBAR) = EVALUATE_RAMP(T,DUMMY,I_RAMP_GY)*GVEC(2)
   GZ(0:IBAR) = EVALUATE_RAMP(T,DUMMY,I_RAMP_GZ)*GVEC(3)
ELSE
   DO I=0,IBAR
      GX(I) = EVALUATE_RAMP(X(I),DUMMY,I_RAMP_GX)*GVEC(1)
      GY(I) = EVALUATE_RAMP(X(I),DUMMY,I_RAMP_GY)*GVEC(2)
      GZ(I) = EVALUATE_RAMP(X(I),DUMMY,I_RAMP_GZ)*GVEC(3)
   ENDDO
ENDIF

! Compute x-direction flux term FVX

!$OMP PARALLEL PRIVATE(WP, WM, VP, VM, UP, UM, &
!$OMP& OMXP, OMXM, OMYP, OMYM, OMZP, OMZM, &
!$OMP& TXZP, TXZM, TXYP, TXYM, TYZP, TYZM, &
!$OMP& IC, IEXP, IEXM, IEYP, IEYM, IEZP, IEZM, &
!$OMP& RRHO, DUDX, DVDY, DWDZ, VTRM)
!$OMP DO SCHEDULE(static) &
!$OMP& PRIVATE(WOMY, VOMZ, TXXP, TXXM, DTXXDX, DTXYDY, DTXZDZ)
DO K=1,KBAR
   DO J=1,JBAR
      DO I=0,IBAR
         WP    = WW(I,J,K)   + WW(I+1,J,K)
         WM    = WW(I,J,K-1) + WW(I+1,J,K-1)
         VP    = VV(I,J,K)   + VV(I+1,J,K)
         VM    = VV(I,J-1,K) + VV(I+1,J-1,K)
         OMYP  = OMY(I,J,K)
         OMYM  = OMY(I,J,K-1)
         OMZP  = OMZ(I,J,K)
         OMZM  = OMZ(I,J-1,K)
         TXZP  = TXZ(I,J,K)
         TXZM  = TXZ(I,J,K-1)
         TXYP  = TXY(I,J,K)
         TXYM  = TXY(I,J-1,K)
         IC    = CELL_INDEX(I,J,K)
         IEYP  = EDGE_INDEX(8,IC)
         IEYM  = EDGE_INDEX(6,IC)
         IEZP  = EDGE_INDEX(12,IC)
         IEZM  = EDGE_INDEX(10,IC)
         IF (OME_E(-1,IEYP)>-1.E5_EB) THEN
            OMYP = OME_E(-1,IEYP)
            TXZP = TAU_E(-1,IEYP)
         ENDIF
         IF (OME_E( 1,IEYM)>-1.E5_EB) THEN
            OMYM = OME_E( 1,IEYM)
            TXZM = TAU_E( 1,IEYM)
         ENDIF
         IF (OME_E(-2,IEZP)>-1.E5_EB) THEN
            OMZP = OME_E(-2,IEZP)
            TXYP = TAU_E(-2,IEZP)
         ENDIF
         IF (OME_E( 2,IEZM)>-1.E5_EB) THEN
            OMZM = OME_E( 2,IEZM)
            TXYM = TAU_E( 2,IEZM)
         ENDIF
         WOMY  = WP*OMYP + WM*OMYM
         VOMZ  = VP*OMZP + VM*OMZM
         RRHO  = 2._EB/(RHOP(I,J,K)+RHOP(I+1,J,K))
         DVDY  = (VV(I+1,J,K)-VV(I+1,J-1,K))*RDY(J)
         DWDZ  = (WW(I+1,J,K)-WW(I+1,J,K-1))*RDZ(K)
         TXXP  = MU(I+1,J,K)*( FOTH*DP(I+1,J,K) - 2._EB*(DVDY+DWDZ) )
         DVDY  = (VV(I,J,K)-VV(I,J-1,K))*RDY(J)
         DWDZ  = (WW(I,J,K)-WW(I,J,K-1))*RDZ(K)
         TXXM  = MU(I,J,K)  *( FOTH*DP(I,J,K)   - 2._EB*(DVDY+DWDZ) )
         DTXXDX= RDXN(I)*(TXXP-TXXM)
         DTXYDY= RDY(J) *(TXYP-TXYM)
         DTXZDZ= RDZ(K) *(TXZP-TXZM)
         VTRM  = DTXXDX + DTXYDY + DTXZDZ
         FVX(I,J,K) = 0.25_EB*(WOMY - VOMZ) - GX(I) + RRHO*(GX(I)*RHO_0(K) - VTRM)
      ENDDO
   ENDDO
ENDDO
!$OMP END DO NOWAIT

! Compute y-direction flux term FVY

!$OMP DO SCHEDULE(static) &
!$OMP& PRIVATE(WOMX, UOMZ, TYYP, TYYM, DTXYDX, DTYYDY, DTYZDZ)
DO K=1,KBAR
   DO J=0,JBAR
      DO I=1,IBAR
         UP    = UU(I,J,K)   + UU(I,J+1,K)
         UM    = UU(I-1,J,K) + UU(I-1,J+1,K)
         WP    = WW(I,J,K)   + WW(I,J+1,K)
         WM    = WW(I,J,K-1) + WW(I,J+1,K-1)
         OMXP  = OMX(I,J,K)
         OMXM  = OMX(I,J,K-1)
         OMZP  = OMZ(I,J,K)
         OMZM  = OMZ(I-1,J,K)
         TYZP  = TYZ(I,J,K)
         TYZM  = TYZ(I,J,K-1)
         TXYP  = TXY(I,J,K)
         TXYM  = TXY(I-1,J,K)
         IC    = CELL_INDEX(I,J,K)
         IEXP  = EDGE_INDEX(4,IC)
         IEXM  = EDGE_INDEX(2,IC)
         IEZP  = EDGE_INDEX(12,IC)
         IEZM  = EDGE_INDEX(11,IC)
         IF (OME_E(-2,IEXP)>-1.E5_EB) THEN
            OMXP = OME_E(-2,IEXP)
            TYZP = TAU_E(-2,IEXP)
         ENDIF
         IF (OME_E( 2,IEXM)>-1.E5_EB) THEN
            OMXM = OME_E( 2,IEXM)
            TYZM = TAU_E( 2,IEXM)
         ENDIF
         IF (OME_E(-1,IEZP)>-1.E5_EB) THEN
            OMZP = OME_E(-1,IEZP)
            TXYP = TAU_E(-1,IEZP)
         ENDIF
         IF (OME_E( 1,IEZM)>-1.E5_EB) THEN
            OMZM = OME_E( 1,IEZM)
            TXYM = TAU_E( 1,IEZM)
         ENDIF
         WOMX  = WP*OMXP + WM*OMXM
         UOMZ  = UP*OMZP + UM*OMZM
         RRHO  = 2._EB/(RHOP(I,J,K)+RHOP(I,J+1,K))
         DUDX  = (UU(I,J+1,K)-UU(I-1,J+1,K))*RDX(I)
         DWDZ  = (WW(I,J+1,K)-WW(I,J+1,K-1))*RDZ(K)
         TYYP  = MU(I,J+1,K)*( FOTH*DP(I,J+1,K) - 2._EB*(DUDX+DWDZ) )
         DUDX  = (UU(I,J,K)-UU(I-1,J,K))*RDX(I)
         DWDZ  = (WW(I,J,K)-WW(I,J,K-1))*RDZ(K)
         TYYM  = MU(I,J,K)  *( FOTH*DP(I,J,K)   - 2._EB*(DUDX+DWDZ) )
         DTXYDX= RDX(I) *(TXYP-TXYM)
         DTYYDY= RDYN(J)*(TYYP-TYYM)
         DTYZDZ= RDZ(K) *(TYZP-TYZM)
         VTRM  = DTXYDX + DTYYDY + DTYZDZ
         FVY(I,J,K) = 0.25_EB*(UOMZ - WOMX) - GY(I) + RRHO*(GY(I)*RHO_0(K) - VTRM)
      ENDDO
   ENDDO
ENDDO
!$OMP END DO NOWAIT

! Compute z-direction flux term FVZ

!$OMP DO SCHEDULE(static) &
!$OMP& PRIVATE(UOMY, VOMX, TZZP, TZZM, DTXZDX, DTYZDY, DTZZDZ)
DO K=0,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         UP    = UU(I,J,K)   + UU(I,J,K+1)
         UM    = UU(I-1,J,K) + UU(I-1,J,K+1)
         VP    = VV(I,J,K)   + VV(I,J,K+1)
         VM    = VV(I,J-1,K) + VV(I,J-1,K+1)
         OMYP  = OMY(I,J,K)
         OMYM  = OMY(I-1,J,K)
         OMXP  = OMX(I,J,K)
         OMXM  = OMX(I,J-1,K)
         TXZP  = TXZ(I,J,K)
         TXZM  = TXZ(I-1,J,K)
         TYZP  = TYZ(I,J,K)
         TYZM  = TYZ(I,J-1,K)
         IC    = CELL_INDEX(I,J,K)
         IEXP  = EDGE_INDEX(4,IC)
         IEXM  = EDGE_INDEX(3,IC)
         IEYP  = EDGE_INDEX(8,IC)
         IEYM  = EDGE_INDEX(7,IC)
         IF (OME_E(-1,IEXP)>-1.E5_EB) THEN
            OMXP = OME_E(-1,IEXP)
            TYZP = TAU_E(-1,IEXP)
         ENDIF
         IF (OME_E( 1,IEXM)>-1.E5_EB) THEN
            OMXM = OME_E( 1,IEXM)
            TYZM = TAU_E( 1,IEXM)
         ENDIF
         IF (OME_E(-2,IEYP)>-1.E5_EB) THEN
            OMYP = OME_E(-2,IEYP)
            TXZP = TAU_E(-2,IEYP)
         ENDIF
         IF (OME_E( 2,IEYM)>-1.E5_EB) THEN
            OMYM = OME_E( 2,IEYM)
            TXZM = TAU_E( 2,IEYM)
         ENDIF
         UOMY  = UP*OMYP + UM*OMYM
         VOMX  = VP*OMXP + VM*OMXM
         RRHO  = 2._EB/(RHOP(I,J,K)+RHOP(I,J,K+1))
         DUDX  = (UU(I,J,K+1)-UU(I-1,J,K+1))*RDX(I)
         DVDY  = (VV(I,J,K+1)-VV(I,J-1,K+1))*RDY(J)
         TZZP  = MU(I,J,K+1)*( FOTH*DP(I,J,K+1) - 2._EB*(DUDX+DVDY) )
         DUDX  = (UU(I,J,K)-UU(I-1,J,K))*RDX(I)
         DVDY  = (VV(I,J,K)-VV(I,J-1,K))*RDY(J)
         TZZM  = MU(I,J,K)  *( FOTH*DP(I,J,K)   - 2._EB*(DUDX+DVDY) )
         DTXZDX= RDX(I) *(TXZP-TXZM)
         DTYZDY= RDY(J) *(TYZP-TYZM)
         DTZZDZ= RDZN(K)*(TZZP-TZZM)
         VTRM  = DTXZDX + DTYZDY + DTZZDZ
         FVZ(I,J,K) = 0.25_EB*(VOMX - UOMY) - GZ(I) + RRHO*(GZ(I)*0.5_EB*(RHO_0(K)+RHO_0(K+1)) - VTRM)
      ENDDO
   ENDDO
ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

IF (EVACUATION_ONLY(NM)) THEN
   FVZ = 0._EB
   RETURN
END IF

! Restore previous substep velocities to gas cut-faces underlaying Cartesian faces.
IF (CC_IBM) CALL CCIBM_INTERP_FACE_VEL(DT,NM,.FALSE.)

! Additional force terms

IF (ANY(MEAN_FORCING))             CALL MOMENTUM_NUDGING           ! Mean forcing
IF (ANY(ABS(FVEC)>TWO_EPSILON_EB)) CALL DIRECT_FORCE               ! Direct force
IF (ANY(ABS(OVEC)>TWO_EPSILON_EB)) CALL CORIOLIS_FORCE             ! Coriolis force
IF (PATCH_VELOCITY)                CALL PATCH_VELOCITY_FLUX(DT,NM) ! Specified patch velocity
IF (PERIODIC_TEST==7)              CALL MMS_VELOCITY_FLUX(NM,T)    ! Source term in manufactured solution
IF (PERIODIC_TEST==21 .OR. PERIODIC_TEST==22 .OR. PERIODIC_TEST==23) CALL ROTATED_CUBE_VELOCITY_FLUX(NM,T)


CONTAINS

SUBROUTINE MOMENTUM_NUDGING

USE COMPLEX_GEOMETRY, ONLY : IBM_GASPHASE, IBM_FGSC

! Add a force vector to the momentum equation that moves the flow field towards the direction of the mean flow.

REAL(EB) :: UBAR,VBAR,WBAR,INTEGRAL,SUM_VOLUME,VC,UMEAN,VMEAN,WMEAN,DU_FORCING,DV_FORCING,DW_FORCING,DT_LOC
INTEGER  :: NSC,I_LO,J_LO,I_HI,J_HI

IF (ICYC==1) RETURN ! need one cycle to initialize forcing arrays

DT_LOC = MAX(DT,DT_MEAN_FORCING)
NSC    = SPONGE_CELLS

MEAN_FORCING_X: IF (MEAN_FORCING(1)) THEN
   SELECT_RAMP_U: SELECT CASE(I_RAMP_U0_Z)
      CASE(0) SELECT_RAMP_U
         PREDICTOR_IF_U: IF (PREDICTOR) THEN
            INTEGRAL = 0._EB
            SUM_VOLUME = 0._EB
            DO K=1,KBAR
               DO J=1,JBAR
                  DO I=0,IBAR
                     IC1 = CELL_INDEX(I,J,K)
                     IC2 = CELL_INDEX(I+1,J,K)
                     IF (SOLID(IC1)) CYCLE
                     IF (SOLID(IC2)) CYCLE
                     IF (CC_IBM) THEN
                        IF(FCVAR(I,J,K,IBM_FGSC,IAXIS) /= IBM_GASPHASE) CYCLE ! If face not regular gasphase type cycle.
                     ENDIF
                     IF (.NOT.MEAN_FORCING_CELL(I,J,K)  ) CYCLE
                     IF (.NOT.MEAN_FORCING_CELL(I+1,J,K)) CYCLE
                     VC = DXN(I)*DY(J)*DZ(K)
                     INTEGRAL = INTEGRAL + UU(I,J,K)*VC
                     SUM_VOLUME = SUM_VOLUME + VC
                  ENDDO
               ENDDO
            ENDDO
            IF (SUM_VOLUME>TWO_EPSILON_EB) THEN
               U_MEAN_LOC(NM) = INTEGRAL
            ELSE
               U_MEAN_LOC(NM) = 0._EB
            ENDIF
            M_VOLU_LOC(NM) = SUM_VOLUME
         ENDIF PREDICTOR_IF_U
         UBAR = U0*EVALUATE_RAMP(T,DUMMY,I_RAMP_U0_T)
         DU_FORCING = (UBAR-U_MEAN_GLOBAL)/DT_LOC
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=0,IBAR
                  IF (.NOT.MEAN_FORCING_CELL(I,J,K)  ) CYCLE
                  IF (.NOT.MEAN_FORCING_CELL(I+1,J,K)) CYCLE
                  FVX(I,J,K) = FVX(I,J,K) - DU_FORCING
               ENDDO
            ENDDO
         ENDDO
      CASE(1:) SELECT_RAMP_U
         K_LOOP_U: DO K=1,KBAR
            INTEGRAL = 0._EB
            SUM_VOLUME = 0._EB
            DO J=1,JBAR
               DO I=0,IBAR
                  IC1 = CELL_INDEX(I,J,K)
                  IC2 = CELL_INDEX(I+1,J,K)
                  IF (SOLID(IC1)) CYCLE
                  IF (SOLID(IC2)) CYCLE
                  IF (CC_IBM) THEN
                     IF(FCVAR(I,J,K,IBM_FGSC,IAXIS) /= IBM_GASPHASE) CYCLE
                  ENDIF
                  VC = DXN(I)*DY(J)*DZ(K)
                  INTEGRAL = INTEGRAL + UU(I,J,K)*VC
                  SUM_VOLUME = SUM_VOLUME + VC
               ENDDO
            ENDDO
            IF (SUM_VOLUME>TWO_EPSILON_EB) THEN
               UMEAN = INTEGRAL/SUM_VOLUME
            ELSE
               ! this can happen if all cells in a given row, k, are solid
               UMEAN = 0._EB
            ENDIF
            UBAR = U0*EVALUATE_RAMP(T,DUMMY,I_RAMP_U0_T)*EVALUATE_RAMP(ZC(K),DUMMY,I_RAMP_U0_Z)
            DU_FORCING = (UBAR-UMEAN)/DT_LOC
            ! Apply the average force term to bulk of domain, and apply more aggressive forcing at boundary
            I_LO = 0
            I_HI = IBAR
            IF (APPLY_SPONGE_LAYER(1)) THEN
               FVX(0:NSC-1,:,K)         = FVX(0:NSC-1,:,K)         - (UBAR-UU(0:NSC-1,:,K))/DT_LOC
               I_LO = NSC
            ENDIF
            IF (APPLY_SPONGE_LAYER(-1)) THEN
               FVX(IBAR-NSC+1:IBAR,:,K) = FVX(IBAR-NSC+1:IBAR,:,K) - (UBAR-UU(IBAR-NSC+1:IBAR,:,K))/DT_LOC
               I_HI = IBAR-NSC
            ENDIF
            FVX(I_LO:I_HI,:,K) = FVX(I_LO:I_HI,:,K) - DU_FORCING
         ENDDO K_LOOP_U
   END SELECT SELECT_RAMP_U
ENDIF MEAN_FORCING_X

MEAN_FORCING_Y: IF (MEAN_FORCING(2)) THEN
   SELECT_RAMP_V: SELECT CASE(I_RAMP_V0_Z)
      CASE(0) SELECT_RAMP_V
         PREDICTOR_IF_V: IF (PREDICTOR) THEN
            INTEGRAL = 0._EB
            SUM_VOLUME = 0._EB
            DO K=1,KBAR
               DO J=0,JBAR
                  DO I=1,IBAR
                     IC1 = CELL_INDEX(I,J,K)
                     IC2 = CELL_INDEX(I,J+1,K)
                     IF (SOLID(IC1)) CYCLE
                     IF (SOLID(IC2)) CYCLE
                     IF (CC_IBM) THEN
                        IF(FCVAR(I,J,K,IBM_FGSC,JAXIS) /= IBM_GASPHASE) CYCLE
                     ENDIF
                     IF (.NOT.MEAN_FORCING_CELL(I,J,K)  ) CYCLE
                     IF (.NOT.MEAN_FORCING_CELL(I,J+1,K)) CYCLE
                     VC = DX(I)*DYN(J)*DZ(K)
                     INTEGRAL = INTEGRAL + VV(I,J,K)*VC
                     SUM_VOLUME = SUM_VOLUME + VC
                  ENDDO
               ENDDO
            ENDDO
            IF (SUM_VOLUME>TWO_EPSILON_EB) THEN
               V_MEAN_LOC(NM) = INTEGRAL
            ELSE
               V_MEAN_LOC(NM) = 0._EB
            ENDIF
            M_VOLV_LOC(NM) = SUM_VOLUME
         ENDIF PREDICTOR_IF_V
         VBAR = V0*EVALUATE_RAMP(T,DUMMY,I_RAMP_V0_T)
         DV_FORCING = (VBAR-V_MEAN_GLOBAL)/DT_LOC
         DO K=1,KBAR
            DO J=0,JBAR
               DO I=1,IBAR
                  IF (.NOT.MEAN_FORCING_CELL(I,J,K)  ) CYCLE
                  IF (.NOT.MEAN_FORCING_CELL(I,J+1,K)) CYCLE
                  FVY(I,J,K) = FVY(I,J,K) - DV_FORCING
               ENDDO
            ENDDO
         ENDDO
      CASE(1:) SELECT_RAMP_V
         K_LOOP_V: DO K=1,KBAR
            INTEGRAL = 0._EB
            SUM_VOLUME = 0._EB
            DO J=0,JBAR
               DO I=1,IBAR
                  IC1 = CELL_INDEX(I,J,K)
                  IC2 = CELL_INDEX(I,J+1,K)
                  IF (SOLID(IC1)) CYCLE
                  IF (SOLID(IC2)) CYCLE
                  IF (CC_IBM) THEN
                     IF(FCVAR(I,J,K,IBM_FGSC,JAXIS) /= IBM_GASPHASE) CYCLE
                  ENDIF
                  VC = DX(I)*DYN(J)*DZ(K)
                  INTEGRAL = INTEGRAL + VV(I,J,K)*VC
                  SUM_VOLUME = SUM_VOLUME + VC
               ENDDO
            ENDDO
            IF (SUM_VOLUME>TWO_EPSILON_EB) THEN
               VMEAN = INTEGRAL/SUM_VOLUME
            ELSE
               VMEAN = 0._EB
            ENDIF
            VBAR = V0*EVALUATE_RAMP(T,DUMMY,I_RAMP_V0_T)*EVALUATE_RAMP(ZC(K),DUMMY,I_RAMP_V0_Z)
            DV_FORCING = (VBAR-VMEAN)/DT_LOC
            ! Apply the average force term to bulk of domain, and apply more aggressive forcing at boundary
            J_LO = 0
            J_HI = JBAR
            IF (APPLY_SPONGE_LAYER(2)) THEN
               FVY(:,0:NSC-1,K)         = FVY(:,0:NSC-1,K)         - (VBAR-VV(:,0:NSC-1,K))/DT_LOC
               J_LO = NSC
            ENDIF
            IF (APPLY_SPONGE_LAYER(-2)) THEN
               FVY(:,JBAR-NSC+1:JBAR,K) = FVY(:,JBAR-NSC+1:JBAR,K) - (VBAR-VV(:,JBAR-NSC+1:JBAR,K))/DT_LOC
               J_HI = JBAR-NSC
            ENDIF
            FVY(:,J_LO:J_HI,K) = FVY(:,J_LO:J_HI,K) - DV_FORCING
         ENDDO K_LOOP_V
   END SELECT SELECT_RAMP_V
ENDIF MEAN_FORCING_Y

MEAN_FORCING_Z: IF (MEAN_FORCING(3)) THEN
   SELECT_RAMP_W: SELECT CASE(I_RAMP_W0_Z)
      CASE(0) SELECT_RAMP_W
         PREDICTOR_IF_W: IF (PREDICTOR) THEN
            INTEGRAL = 0._EB
            SUM_VOLUME = 0._EB
            DO K=0,KBAR
               DO J=1,JBAR
                  DO I=1,IBAR
                     IC1 = CELL_INDEX(I,J,K)
                     IC2 = CELL_INDEX(I,J,K+1)
                     IF (SOLID(IC1)) CYCLE
                     IF (SOLID(IC2)) CYCLE
                     IF (CC_IBM) THEN
                        IF(FCVAR(I,J,K,IBM_FGSC,KAXIS) /= IBM_GASPHASE) CYCLE
                     ENDIF
                     IF (.NOT.MEAN_FORCING_CELL(I,J,K)  ) CYCLE
                     IF (.NOT.MEAN_FORCING_CELL(I,J,K+1)) CYCLE
                     VC = DX(I)*DY(J)*DZN(K)
                     INTEGRAL = INTEGRAL + WW(I,J,K)*VC
                     SUM_VOLUME = SUM_VOLUME + VC
                  ENDDO
               ENDDO
            ENDDO
            IF (SUM_VOLUME>TWO_EPSILON_EB) THEN
               W_MEAN_LOC(NM) = INTEGRAL
            ELSE
               W_MEAN_LOC(NM) = 0._EB
            ENDIF
            M_VOLW_LOC(NM) = SUM_VOLUME
         ENDIF PREDICTOR_IF_W
         WBAR = W0*EVALUATE_RAMP(T,DUMMY,I_RAMP_W0_T)
         DW_FORCING = (WBAR-W_MEAN_GLOBAL)/DT_LOC
         DO K=0,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  IF (.NOT.MEAN_FORCING_CELL(I,J,K)  ) CYCLE
                  IF (.NOT.MEAN_FORCING_CELL(I,J,K+1)) CYCLE
                  FVZ(I,J,K) = FVZ(I,J,K) - DW_FORCING
               ENDDO
            ENDDO
         ENDDO
      CASE(1:)
         K_LOOP_W: DO K=0,KBAR
            INTEGRAL = 0._EB
            SUM_VOLUME = 0._EB
            DO J=1,JBAR
               DO I=1,IBAR
                  IC1 = CELL_INDEX(I,J,K)
                  IC2 = CELL_INDEX(I,J,K+1)
                  IF (SOLID(IC1)) CYCLE
                  IF (SOLID(IC2)) CYCLE
                  IF (CC_IBM) THEN
                     IF(FCVAR(I,J,K,IBM_FGSC,KAXIS) /= IBM_GASPHASE) CYCLE
                  ENDIF
                  VC = DX(I)*DY(J)*DZN(K)
                  INTEGRAL = INTEGRAL + WW(I,J,K)*VC
                  SUM_VOLUME = SUM_VOLUME + VC
               ENDDO
            ENDDO
            IF (SUM_VOLUME>TWO_EPSILON_EB) THEN
               WMEAN = INTEGRAL/SUM_VOLUME
            ELSE
               WMEAN = 0._EB
            ENDIF
            WBAR = W0*EVALUATE_RAMP(T,DUMMY,I_RAMP_W0_T)*EVALUATE_RAMP(Z(K),DUMMY,I_RAMP_W0_Z)
            DW_FORCING = (WBAR-WMEAN)/DT_LOC
            ! Apply the average force term to bulk of domain, and apply more aggressive forcing at boundary
            IF (APPLY_SPONGE_LAYER(-3) .AND. K==KBAR) THEN
               DO J=1,JBAR
                  DO I=1,IBAR
                     FVZ(I,J,K) = FVZ(I,J,K) - (WBAR-WW(I,J,K))/DT_LOC
                  ENDDO
               ENDDO
            ELSE
               FVZ(:,:,K) = FVZ(:,:,K) - DW_FORCING
            ENDIF
         ENDDO K_LOOP_W
   END SELECT SELECT_RAMP_W
ENDIF MEAN_FORCING_Z

END SUBROUTINE MOMENTUM_NUDGING


SUBROUTINE DIRECT_FORCE()
REAL(EB) :: TIME_RAMP_FACTOR

TIME_RAMP_FACTOR = EVALUATE_RAMP(T,DUMMY,I_RAMP_FVX_T)
!$OMP PARALLEL DO PRIVATE(RRHO) SCHEDULE(STATIC)
DO K=1,KBAR
   DO J=1,JBAR
      DO I=0,IBAR
         RRHO = 2._EB/(RHOP(I,J,K)+RHOP(I+1,J,K))
         FVX(I,J,K) = FVX(I,J,K) - RRHO*FVEC(1)*TIME_RAMP_FACTOR
      ENDDO
   ENDDO
ENDDO
!$OMP END PARALLEL DO

TIME_RAMP_FACTOR = EVALUATE_RAMP(T,DUMMY,I_RAMP_FVY_T)
!$OMP PARALLEL DO PRIVATE(RRHO) SCHEDULE(STATIC)
DO K=1,KBAR
   DO J=0,JBAR
      DO I=1,IBAR
         RRHO = 2._EB/(RHOP(I,J,K)+RHOP(I,J+1,K))
         FVY(I,J,K) = FVY(I,J,K) - RRHO*FVEC(2)*TIME_RAMP_FACTOR
      ENDDO
   ENDDO
ENDDO
!$OMP END PARALLEL DO

TIME_RAMP_FACTOR = EVALUATE_RAMP(T,DUMMY,I_RAMP_FVZ_T)
!$OMP PARALLEL DO PRIVATE(RRHO) SCHEDULE(STATIC)
DO K=0,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         RRHO = 2._EB/(RHOP(I,J,K)+RHOP(I,J,K+1))
         FVZ(I,J,K) = FVZ(I,J,K) - RRHO*FVEC(3)*TIME_RAMP_FACTOR
      ENDDO
   ENDDO
ENDDO
!$OMP END PARALLEL DO

END SUBROUTINE DIRECT_FORCE


SUBROUTINE CORIOLIS_FORCE()

REAL(EB), POINTER, DIMENSION(:,:,:) :: UP=>NULL(),VP=>NULL(),WP=>NULL()
REAL(EB) :: UBAR,VBAR,WBAR
INTEGER :: II,JJ,KK,IW
TYPE(WALL_TYPE), POINTER :: WC=>NULL()

! Velocities relative to the p-cell center (same work done in Deardorff eddy viscosity)

UP => WORK7
VP => WORK8
WP => WORK9
UP=0._EB
VP=0._EB
WP=0._EB

!$OMP PARALLEL DO SCHEDULE(static)
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
         UP(I,J,K) = 0.5_EB*(UU(I,J,K) + UU(I-1,J,K))
         VP(I,J,K) = 0.5_EB*(VV(I,J,K) + VV(I,J-1,K))
         WP(I,J,K) = 0.5_EB*(WW(I,J,K) + WW(I,J,K-1))
      ENDDO
   ENDDO
ENDDO
!$OMP END PARALLEL DO

DO IW=1,N_EXTERNAL_WALL_CELLS
   WC=>WALL(IW)
   II = WC%ONE_D%II
   JJ = WC%ONE_D%JJ
   KK = WC%ONE_D%KK
   UP(II,JJ,KK) = U_GHOST(IW)
   VP(II,JJ,KK) = V_GHOST(IW)
   WP(II,JJ,KK) = W_GHOST(IW)
ENDDO

! x momentum

!$OMP PARALLEL DO PRIVATE(VBAR,WBAR) SCHEDULE(STATIC)
DO K=1,KBAR
   DO J=1,JBAR
      DO I=0,IBAR
         VBAR = 0.5_EB*(VP(I,J,K)+VP(I+1,J,K))
         WBAR = 0.5_EB*(WP(I,J,K)+WP(I+1,J,K))
         FVX(I,J,K) = FVX(I,J,K) + 2._EB*(OVEC(2)*WBAR-OVEC(3)*VBAR)
      ENDDO
   ENDDO
ENDDO
!$OMP END PARALLEL DO

! y momentum

!$OMP PARALLEL DO PRIVATE(UBAR,WBAR) SCHEDULE(STATIC)
DO K=1,KBAR
   DO J=0,JBAR
      DO I=1,IBAR
         UBAR = 0.5_EB*(UP(I,J,K)+UP(I,J+1,K))
         WBAR = 0.5_EB*(WP(I,J,K)+WP(I,J+1,K))
         FVY(I,J,K) = FVY(I,J,K) + 2._EB*(OVEC(3)*UBAR - OVEC(1)*WBAR)
      ENDDO
   ENDDO
ENDDO
!$OMP END PARALLEL DO

! z momentum

!$OMP PARALLEL DO PRIVATE(UBAR,VBAR) SCHEDULE(STATIC)
DO K=0,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         UBAR = 0.5_EB*(UP(I,J,K)+UP(I,J,K+1))
         VBAR = 0.5_EB*(VP(I,J,K)+VP(I,J,K+1))
         FVZ(I,J,K) = FVZ(I,J,K) + 2._EB*(OVEC(1)*VBAR - OVEC(2)*UBAR)
      ENDDO
   ENDDO
ENDDO
!$OMP END PARALLEL DO

END SUBROUTINE CORIOLIS_FORCE

END SUBROUTINE VELOCITY_FLUX


SUBROUTINE MMS_VELOCITY_FLUX(NM,T)

! Shunn et al., JCP (2012) prob 3

USE MANUFACTURED_SOLUTIONS, ONLY: VD2D_MMS_U_SRC_3,VD2D_MMS_V_SRC_3
INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: T
INTEGER :: I,J,K
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU=>NULL(),WW=>NULL()

CALL POINT_TO_MESH(NM)

IF (PREDICTOR) THEN
   UU=>U
   WW=>W
ELSE
   UU=>US
   WW=>WS
ENDIF

DO K=1,KBAR
   DO J=1,JBAR
      DO I=0,IBAR
         FVX(I,J,K) = FVX(I,J,K) - VD2D_MMS_U_SRC_3(X(I),ZC(K),T)
      ENDDO
   ENDDO
ENDDO

DO K=0,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         FVZ(I,J,K) = FVZ(I,J,K) - VD2D_MMS_V_SRC_3(XC(I),Z(K),T)
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE MMS_VELOCITY_FLUX


SUBROUTINE VELOCITY_FLUX_CYLINDRICAL(T,NM)

! Compute convective and diffusive terms for 2D axisymmetric

USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
REAL(EB) :: T,DMUDX
INTEGER :: I0
INTEGER, INTENT(IN) :: NM
REAL(EB) :: MUY,UP,UM,WP,WM,VTRM,DTXZDZ,DTXZDX,DUDX,DWDZ,DUDZ,DWDX,WOMY,UOMY,OMYP,OMYM,TXZP,TXZM, &
            AH,RRHO,GX,GZ,TXXP,TXXM,TZZP,TZZM,DTXXDX,DTZZDZ,DUMMY=0._EB
INTEGER :: I,J,K,IEYP,IEYM,IC
REAL(EB), POINTER, DIMENSION(:,:,:) :: TXZ=>NULL(),OMY=>NULL(),UU=>NULL(),WW=>NULL(),RHOP=>NULL(),DP=>NULL()

CALL POINT_TO_MESH(NM)

IF (PREDICTOR) THEN
   UU => U
   WW => W
   DP => D
   RHOP => RHO
ELSE
   UU => US
   WW => WS
   DP => DS
   RHOP => RHOS
ENDIF

TXZ => WORK2
OMY => WORK5

! Compute vorticity and stress tensor components

DO K=0,KBAR
   DO J=0,JBAR
      DO I=0,IBAR
         DUDZ = RDZN(K)*(UU(I,J,K+1)-UU(I,J,K))
         DWDX = RDXN(I)*(WW(I+1,J,K)-WW(I,J,K))
         OMY(I,J,K) = DUDZ - DWDX
         MUY = 0.25_EB*(MU(I+1,J,K)+MU(I,J,K)+MU(I,J,K+1)+MU(I+1,J,K+1))
         TXZ(I,J,K) = MUY*(DUDZ + DWDX)
      ENDDO
   ENDDO
ENDDO

! Compute gravity components

GX  = 0._EB
GZ  = EVALUATE_RAMP(T,DUMMY,I_RAMP_GZ)*GVEC(3)

! Compute r-direction flux term FVX

IF (ABS(XS)<=TWO_EPSILON_EB) THEN
   I0 = 1
ELSE
   I0 = 0
ENDIF

J = 1

DO K= 1,KBAR
   DO I=I0,IBAR
      WP    = WW(I,J,K)   + WW(I+1,J,K)
      WM    = WW(I,J,K-1) + WW(I+1,J,K-1)
      OMYP  = OMY(I,J,K)
      OMYM  = OMY(I,J,K-1)
      TXZP  = TXZ(I,J,K)
      TXZM  = TXZ(I,J,K-1)
      IC    = CELL_INDEX(I,J,K)
      IEYP  = EDGE_INDEX(8,IC)
      IEYM  = EDGE_INDEX(6,IC)
      IF (OME_E(-1,IEYP)>-1.E5_EB) THEN
         OMYP = OME_E(-1,IEYP)
         TXZP = TAU_E(-1,IEYP)
      ENDIF
      IF (OME_E( 1,IEYM)>-1.E5_EB) THEN
         OMYM = OME_E( 1,IEYM)
         TXZM = TAU_E( 1,IEYM)
      ENDIF
      WOMY  = WP*OMYP + WM*OMYM
      RRHO  = 2._EB/(RHOP(I,J,K)+RHOP(I+1,J,K))
      AH    = RHO_0(K)*RRHO - 1._EB
      DWDZ  = (WW(I+1,J,K)-WW(I+1,J,K-1))*RDZ(K)
      TXXP  = MU(I+1,J,K)*( FOTH*DP(I+1,J,K) - 2._EB*DWDZ )
      DWDZ  = (WW(I,J,K)-WW(I,J,K-1))*RDZ(K)
      TXXM  = MU(I,J,K)  *( FOTH*DP(I,J,K) -2._EB*DWDZ )
      DTXXDX= RDXN(I)*(TXXP-TXXM)
      DTXZDZ= RDZ(K) *(TXZP-TXZM)
      DMUDX = (MU(I+1,J,K)-MU(I,J,K))*RDXN(I)
      VTRM  = RRHO*( DTXXDX + DTXZDZ - 2._EB*UU(I,J,K)*DMUDX/R(I) )
      FVX(I,J,K) = 0.25_EB*WOMY + GX*AH - VTRM
   ENDDO
ENDDO

! Compute z-direction flux term FVZ

DO K=0,KBAR
   DO I=1,IBAR
      UP    = UU(I,J,K)   + UU(I,J,K+1)
      UM    = UU(I-1,J,K) + UU(I-1,J,K+1)
      OMYP  = OMY(I,J,K)
      OMYM  = OMY(I-1,J,K)
      TXZP  = TXZ(I,J,K)
      TXZM  = TXZ(I-1,J,K)
      IC    = CELL_INDEX(I,J,K)
      IEYP  = EDGE_INDEX(8,IC)
      IEYM  = EDGE_INDEX(7,IC)
      IF (OME_E(-2,IEYP)>-1.E5_EB) THEN
         OMYP = OME_E(-2,IEYP)
         TXZP = TAU_E(-2,IEYP)
      ENDIF
      IF (OME_E( 2,IEYM)>-1.E5_EB) THEN
         OMYM = OME_E( 2,IEYM)
         TXZM = TAU_E( 2,IEYM)
      ENDIF
      UOMY  = UP*OMYP + UM*OMYM
      RRHO  = 2._EB/(RHOP(I,J,K)+RHOP(I,J,K+1))
      AH    = 0.5_EB*(RHO_0(K)+RHO_0(K+1))*RRHO - 1._EB
      DUDX  = (R(I)*UU(I,J,K+1)-R(I-1)*UU(I-1,J,K+1))*RDX(I)*RRN(I)
      TZZP  = MU(I,J,K+1)*( FOTH*DP(I,J,K+1) - 2._EB*DUDX )
      DUDX  = (R(I)*UU(I,J,K)-R(I-1)*UU(I-1,J,K))*RDX(I)*RRN(I)
      TZZM  = MU(I,J,K)  *( FOTH*DP(I,J,K)   - 2._EB*DUDX )
      DTXZDX= RDX(I) *(R(I)*TXZP-R(I-1)*TXZM)*RRN(I)
      DTZZDZ= RDZN(K)*(     TZZP       -TZZM)
      VTRM  = RRHO*(DTXZDX + DTZZDZ)
      FVZ(I,J,K) = -0.25_EB*UOMY + GZ*AH - VTRM
   ENDDO
ENDDO

! Adjust FVX and FVZ at solid, internal obstructions for no flux

END SUBROUTINE VELOCITY_FLUX_CYLINDRICAL


SUBROUTINE NO_FLUX(DT,NM)

! Set FVX,FVY,FVZ inside and on the surface of solid obstructions to maintain no flux

USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: DT
REAL(EB), POINTER, DIMENSION(:,:,:) :: HP=>NULL(),OM_HP=>NULL()
REAL(EB) :: RFODT,H_OTHER,DUUDT,DVVDT,DWWDT,UN,TNOW,DHFCT
INTEGER  :: IC2,IC1,N,I,J,K,IW,II,JJ,KK,IOR,N_INT_CELLS,IIO,JJO,KKO,NOM
TYPE (OBSTRUCTION_TYPE), POINTER :: OB=>NULL()
TYPE (WALL_TYPE), POINTER :: WC=>NULL()
TYPE (EXTERNAL_WALL_TYPE), POINTER :: EWC=>NULL()
LOGICAL :: GLMAT_ON_WHOLE_DOMAIN

TNOW=CURRENT_TIME()
CALL POINT_TO_MESH(NM)

RFODT = RELAXATION_FACTOR/DT

IF (PREDICTOR) HP => H
IF (CORRECTOR) HP => HS

! Exchange H at interpolated boundaries

NO_SCARC_IF: IF (TRIM(PRES_METHOD) /= 'SCARC' .AND. TRIM(PRES_METHOD) /= 'USCARC') THEN

   DO IW=1,N_EXTERNAL_WALL_CELLS
      WC=>WALL(IW)
      EWC=>EXTERNAL_WALL(IW)
      NOM =EWC%NOM
      IF (NOM==0) CYCLE
      IF (PREDICTOR) THEN
         OM_HP=>OMESH(NOM)%H
      ELSE
         OM_HP=>OMESH(NOM)%HS
      ENDIF
      II = WC%ONE_D%II
      JJ = WC%ONE_D%JJ
      KK = WC%ONE_D%KK
      H_OTHER = 0._EB
      DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
         DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
            DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
               H_OTHER = H_OTHER + OM_HP(IIO,JJO,KKO)
            ENDDO
         ENDDO
      ENDDO
      N_INT_CELLS = (EWC%IIO_MAX-EWC%IIO_MIN+1) * (EWC%JJO_MAX-EWC%JJO_MIN+1) * (EWC%KKO_MAX-EWC%KKO_MIN+1)
      HP(II,JJ,KK)  = H_OTHER/REAL(N_INT_CELLS,EB)
   ENDDO

ENDIF NO_SCARC_IF

! Set FVX, FVY and FVZ to drive velocity components at solid boundaries within obstructions towards zero

OBST_LOOP: DO N=1,N_OBST

   OB=>OBSTRUCTION(N)

   DO K=OB%K1+1,OB%K2
      DO J=OB%J1+1,OB%J2
         DO I=OB%I1  ,OB%I2
            IC1 = CELL_INDEX(I,J,K)
            IC2 = CELL_INDEX(I+1,J,K)
            IF (SOLID(IC1) .AND. SOLID(IC2)) THEN
               IF (PREDICTOR) THEN
                  DUUDT = -RFODT*U(I,J,K)
               ELSE
                  DUUDT = -RFODT*(U(I,J,K)+US(I,J,K))
               ENDIF
               FVX(I,J,K) = -RDXN(I)*(HP(I+1,J,K)-HP(I,J,K)) - DUUDT
            ENDIF
         ENDDO
      ENDDO
   ENDDO

   DO K=OB%K1+1,OB%K2
      DO J=OB%J1  ,OB%J2
         DO I=OB%I1+1,OB%I2
            IC1 = CELL_INDEX(I,J,K)
            IC2 = CELL_INDEX(I,J+1,K)
            IF (SOLID(IC1) .AND. SOLID(IC2)) THEN
               IF (PREDICTOR) THEN
                  DVVDT = -RFODT*V(I,J,K)
               ELSE
                  DVVDT = -RFODT*(V(I,J,K)+VS(I,J,K))
               ENDIF
               FVY(I,J,K) = -RDYN(J)*(HP(I,J+1,K)-HP(I,J,K)) - DVVDT
            ENDIF
         ENDDO
      ENDDO
   ENDDO

   DO K=OB%K1  ,OB%K2
      DO J=OB%J1+1,OB%J2
         DO I=OB%I1+1,OB%I2
            IC1 = CELL_INDEX(I,J,K)
            IC2 = CELL_INDEX(I,J,K+1)
            IF (SOLID(IC1) .AND. SOLID(IC2)) THEN
               IF (PREDICTOR) THEN
                  DWWDT = -RFODT*W(I,J,K)
               ELSE
                  DWWDT = -RFODT*(W(I,J,K)+WS(I,J,K))
               ENDIF
               FVZ(I,J,K) = -RDZN(K)*(HP(I,J,K+1)-HP(I,J,K)) - DWWDT
            ENDIF
         ENDDO
      ENDDO
   ENDDO

ENDDO OBST_LOOP

! Set FVX, FVY and FVZ to drive the normal velocity at solid boundaries towards the specified value (U_NORMAL or U_NORMAL_S)
! Logical to define not to apply pressure gradient on external mesh boundaries for GLMAT.

GLMAT_ON_WHOLE_DOMAIN = (PRES_METHOD=='GLMAT') .AND. PRES_ON_WHOLE_DOMAIN

WALL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS

   WC => WALL(IW)

   IF (WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. WC%BOUNDARY_TYPE==OPEN_BOUNDARY) CYCLE WALL_LOOP

   IF (IW<=N_EXTERNAL_WALL_CELLS) THEN
      NOM = EXTERNAL_WALL(IW)%NOM
   ELSE
      NOM = 0
   ENDIF

   IF (IW>N_EXTERNAL_WALL_CELLS .AND. WC%BOUNDARY_TYPE==NULL_BOUNDARY .AND. NOM==0) CYCLE WALL_LOOP

   II  = WC%ONE_D%II
   JJ  = WC%ONE_D%JJ
   KK  = WC%ONE_D%KK
   IOR = WC%ONE_D%IOR

   DHFCT=1._EB
   IF ( (.NOT. PRES_ON_WHOLE_DOMAIN) .OR. (GLMAT_ON_WHOLE_DOMAIN .AND.  IW<=N_EXTERNAL_WALL_CELLS) ) DHFCT=0._EB

   IF (NOM/=0 .OR. WC%BOUNDARY_TYPE==SOLID_BOUNDARY .OR. WC%BOUNDARY_TYPE==NULL_BOUNDARY) THEN
      IF (PREDICTOR) THEN
         UN = -SIGN(1._EB,REAL(IOR,EB))*WC%ONE_D%U_NORMAL_S
      ELSE
         UN = -SIGN(1._EB,REAL(IOR,EB))*WC%ONE_D%U_NORMAL
      ENDIF
      SELECT CASE(IOR)
         CASE( 1)
            IF (PREDICTOR) THEN
               DUUDT = RFODT*(UN-U(II,JJ,KK))
            ELSE
               DUUDT = 2._EB*RFODT*(UN-0.5_EB*(U(II,JJ,KK)+US(II,JJ,KK)) )
            ENDIF
            FVX(II,JJ,KK) = -RDXN(II)*(HP(II+1,JJ,KK)-HP(II,JJ,KK))*DHFCT - DUUDT
         CASE(-1)
            IF (PREDICTOR) THEN
               DUUDT = RFODT*(UN-U(II-1,JJ,KK))
            ELSE
               DUUDT = 2._EB*RFODT*(UN-0.5_EB*(U(II-1,JJ,KK)+US(II-1,JJ,KK)) )
            ENDIF
            FVX(II-1,JJ,KK) = -RDXN(II-1)*(HP(II,JJ,KK)-HP(II-1,JJ,KK))*DHFCT - DUUDT
         CASE( 2)
            IF (PREDICTOR) THEN
               DVVDT = RFODT*(UN-V(II,JJ,KK))
            ELSE
               DVVDT = 2._EB*RFODT*(UN-0.5_EB*(V(II,JJ,KK)+VS(II,JJ,KK)) )
            ENDIF
            FVY(II,JJ,KK) = -RDYN(JJ)*(HP(II,JJ+1,KK)-HP(II,JJ,KK))*DHFCT - DVVDT
         CASE(-2)
            IF (PREDICTOR) THEN
               DVVDT = RFODT*(UN-V(II,JJ-1,KK))
            ELSE
               DVVDT = 2._EB*RFODT*(UN-0.5_EB*(V(II,JJ-1,KK)+VS(II,JJ-1,KK)) )
            ENDIF
            FVY(II,JJ-1,KK) = -RDYN(JJ-1)*(HP(II,JJ,KK)-HP(II,JJ-1,KK))*DHFCT - DVVDT
         CASE( 3)
            IF (PREDICTOR) THEN
               DWWDT = RFODT*(UN-W(II,JJ,KK))
            ELSE
               DWWDT = 2._EB*RFODT*(UN-0.5_EB*(W(II,JJ,KK)+WS(II,JJ,KK)) )
            ENDIF
            FVZ(II,JJ,KK) = -RDZN(KK)*(HP(II,JJ,KK+1)-HP(II,JJ,KK))*DHFCT - DWWDT
         CASE(-3)
            IF (PREDICTOR) THEN
               DWWDT = RFODT*(UN-W(II,JJ,KK-1))
            ELSE
               DWWDT = 2._EB*RFODT*(UN-0.5_EB*(W(II,JJ,KK-1)+WS(II,JJ,KK-1)) )
            ENDIF
            FVZ(II,JJ,KK-1) = -RDZN(KK-1)*(HP(II,JJ,KK)-HP(II,JJ,KK-1))*DHFCT - DWWDT
      END SELECT
   ENDIF

   IF (WC%BOUNDARY_TYPE==MIRROR_BOUNDARY) THEN
      SELECT CASE(IOR)
         CASE( 1)
            FVX(II  ,JJ,KK) = 0._EB
         CASE(-1)
            FVX(II-1,JJ,KK) = 0._EB
         CASE( 2)
            FVY(II  ,JJ,KK) = 0._EB
         CASE(-2)
            FVY(II,JJ-1,KK) = 0._EB
         CASE( 3)
            FVZ(II  ,JJ,KK) = 0._EB
         CASE(-3)
            FVZ(II,JJ,KK-1) = 0._EB
      END SELECT
   ENDIF

ENDDO WALL_LOOP

T_USED(4)=T_USED(4)+CURRENT_TIME()-TNOW
END SUBROUTINE NO_FLUX


SUBROUTINE VELOCITY_PREDICTOR(T,DT,DT_NEW,NM)

USE TURBULENCE, ONLY: COMPRESSION_WAVE
USE MANUFACTURED_SOLUTIONS, ONLY: UF_MMS,WF_MMS,VD2D_MMS_U,VD2D_MMS_V
USE COMPLEX_GEOMETRY, ONLY : CCIBM_VELOCITY_NO_GRADH

! Estimates the velocity components at the next time step

REAL(EB) :: TNOW,XHAT,ZHAT
INTEGER  :: I,J,K
INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: T,DT
REAL(EB) :: DT_NEW(NMESHES)

IF (SOLID_PHASE_ONLY) RETURN
IF (PERIODIC_TEST==4) THEN
   CALL COMPRESSION_WAVE(NM,T,4)
   CALL CHECK_STABILITY(DT,DT_NEW,NM)
   RETURN
ENDIF

TNOW=CURRENT_TIME()
CALL POINT_TO_MESH(NM)

FREEZE_VELOCITY_IF: IF (FREEZE_VELOCITY) THEN
   US = U
   VS = V
   WS = W
ELSE FREEZE_VELOCITY_IF

   DO K=1,KBAR
      DO J=1,JBAR
         DO I=0,IBAR
            US(I,J,K) = U(I,J,K) - DT*( FVX(I,J,K) + RDXN(I)*(H(I+1,J,K)-H(I,J,K)) )
         ENDDO
      ENDDO
   ENDDO

   DO K=1,KBAR
      DO J=0,JBAR
         DO I=1,IBAR
            VS(I,J,K) = V(I,J,K) - DT*( FVY(I,J,K) + RDYN(J)*(H(I,J+1,K)-H(I,J,K)) )
         ENDDO
      ENDDO
   ENDDO

   DO K=0,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            WS(I,J,K) = W(I,J,K) - DT*( FVZ(I,J,K) + RDZN(K)*(H(I,J,K+1)-H(I,J,K)) )
         ENDDO
      ENDDO
   ENDDO

   IF (PRES_METHOD == 'GLMAT' .OR. PRES_METHOD == 'USCARC') CALL WALL_VELOCITY_NO_GRADH(DT,.FALSE.)
   IF (CC_IBM) CALL CCIBM_VELOCITY_NO_GRADH(DT,.FALSE.)

ENDIF FREEZE_VELOCITY_IF

! Manufactured solution (debug)

IF (PERIODIC_TEST==7 .AND. .FALSE.) THEN
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=0,IBAR
            XHAT =  X(I) - UF_MMS*(T)
            ZHAT = ZC(K) - WF_MMS*(T)
            US(I,J,K) = VD2D_MMS_U(XHAT,ZHAT,T)
         ENDDO
      ENDDO
   ENDDO
   DO K=0,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            XHAT = XC(I) - UF_MMS*(T)
            ZHAT =  Z(K) - WF_MMS*(T)
            WS(I,J,K) = VD2D_MMS_V(XHAT,ZHAT,T)
         ENDDO
      ENDDO
   ENDDO
ENDIF

! No vertical velocity in Evacuation meshes

IF (EVACUATION_ONLY(NM)) WS = 0._EB

! Check the stability criteria, and if the time step is too small, send back a signal to kill the job

CALL CHECK_STABILITY(DT,DT_NEW,NM)

IF (DT_NEW(NM)<DT_INITIAL*LIMITING_DT_RATIO .AND. (T+DT_NEW(NM)<(T_END-TWO_EPSILON_EB))) STOP_STATUS = INSTABILITY_STOP

T_USED(4)=T_USED(4)+CURRENT_TIME()-TNOW
END SUBROUTINE VELOCITY_PREDICTOR


SUBROUTINE VELOCITY_CORRECTOR(T,DT,NM)

USE TURBULENCE, ONLY: COMPRESSION_WAVE
USE MANUFACTURED_SOLUTIONS, ONLY: UF_MMS,WF_MMS,VD2D_MMS_U,VD2D_MMS_V
USE COMPLEX_GEOMETRY, ONLY : CCIBM_VELOCITY_NO_GRADH

! Correct the velocity components

REAL(EB) :: TNOW,XHAT,ZHAT
INTEGER  :: I,J,K
INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: T,DT

IF (SOLID_PHASE_ONLY) RETURN
IF (PERIODIC_TEST==4) THEN
   CALL COMPRESSION_WAVE(NM,T,4)
   RETURN
ENDIF

TNOW=CURRENT_TIME()
CALL POINT_TO_MESH(NM)

FREEZE_VELOCITY_IF: IF (FREEZE_VELOCITY) THEN
   U = US
   V = VS
   W = WS
ELSE FREEZE_VELOCITY_IF

   IF (STORE_OLD_VELOCITY) THEN
      U_OLD = U
      V_OLD = V
      W_OLD = W
   ENDIF

   IF (PRES_METHOD == 'GLMAT' .OR. PRES_METHOD == 'USCARC') THEN
      CALL WALL_VELOCITY_NO_GRADH(DT,.TRUE.)                    ! Store U velocities on OBST surfaces.
      IF (CC_IBM) CALL CCIBM_VELOCITY_NO_GRADH(DT,.TRUE.)       ! Store velocities on GEOM SOLID faces.
   ENDIF

   DO K=1,KBAR
      DO J=1,JBAR
         DO I=0,IBAR
            U(I,J,K) = 0.5_EB*( U(I,J,K) + US(I,J,K) - DT*(FVX(I,J,K) + RDXN(I)*(HS(I+1,J,K)-HS(I,J,K))) )
         ENDDO
      ENDDO
   ENDDO

   DO K=1,KBAR
      DO J=0,JBAR
         DO I=1,IBAR
            V(I,J,K) = 0.5_EB*( V(I,J,K) + VS(I,J,K) - DT*(FVY(I,J,K) + RDYN(J)*(HS(I,J+1,K)-HS(I,J,K))) )
         ENDDO
      ENDDO
   ENDDO

   DO K=0,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            W(I,J,K) = 0.5_EB*( W(I,J,K) + WS(I,J,K) - DT*(FVZ(I,J,K) + RDZN(K)*(HS(I,J,K+1)-HS(I,J,K))) )
         ENDDO
      ENDDO
   ENDDO

   IF (PRES_METHOD == 'GLMAT' .OR. PRES_METHOD == 'USCARC') THEN
      CALL WALL_VELOCITY_NO_GRADH(DT,.FALSE.)
      IF (CC_IBM) CALL CCIBM_VELOCITY_NO_GRADH(DT,.FALSE.)
   ENDIF

ENDIF FREEZE_VELOCITY_IF

! Manufactured solution (debug)

IF (PERIODIC_TEST==7 .AND. .FALSE.) THEN
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=0,IBAR
            XHAT =  X(I) - UF_MMS*T
            ZHAT = ZC(K) - WF_MMS*T
            U(I,J,K) = VD2D_MMS_U(XHAT,ZHAT,T)
         ENDDO
      ENDDO
   ENDDO
   DO K=0,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            XHAT = XC(I) - UF_MMS*T
            ZHAT =  Z(K) - WF_MMS*T
            W(I,J,K) = VD2D_MMS_V(XHAT,ZHAT,T)
         ENDDO
      ENDDO
   ENDDO
ENDIF

! No vertical velocity in Evacuation meshes

IF (EVACUATION_ONLY(NM)) W = 0._EB

T_USED(4)=T_USED(4)+CURRENT_TIME()-TNOW
END SUBROUTINE VELOCITY_CORRECTOR


SUBROUTINE VELOCITY_BC(T,NM)

! Assert tangential velocity boundary conditions

USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
USE TURBULENCE, ONLY: WALL_MODEL,WANNIER_FLOW
REAL(EB), INTENT(IN) :: T
REAL(EB) :: MUA,TSI,WGT,TNOW,RAMP_T,OMW,MU_WALL,RHO_WALL,SLIP_COEF,VEL_T,UBAR,VBAR,WBAR, &
            UUP(2),UUM(2),DXX(2),MU_DUIDXJ(-2:2),DUIDXJ(-2:2),PROFILE_FACTOR,VEL_GAS,VEL_GHOST, &
            MU_DUIDXJ_USE(2),DUIDXJ_USE(2),VEL_EDDY,U_TAU,Y_PLUS,WT1,WT2,DUMMY
INTEGER :: I,J,K,NOM(2),IIO(2),JJO(2),KKO(2),IE,II,JJ,KK,IEC,IOR,IWM,IWP,ICMM,ICMP,ICPM,ICPP,IC,ICD,ICDO,IVL,I_SGN,IS, &
           VELOCITY_BC_INDEX,IIGM,JJGM,KKGM,IIGP,JJGP,KKGP,SURF_INDEXM,SURF_INDEXP,ITMP,ICD_SGN,ICDO_SGN, &
           BOUNDARY_TYPE_M,BOUNDARY_TYPE_P,IS2,IWPI,IWMI,VENT_INDEX
LOGICAL :: ALTERED_GRADIENT(-2:2),PROCESS_EDGE,SYNTHETIC_EDDY_METHOD,HVAC_TANGENTIAL,INTERPOLATED_EDGE
INTEGER, INTENT(IN) :: NM
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU=>NULL(),VV=>NULL(),WW=>NULL(),U_Y=>NULL(),U_Z=>NULL(), &
                                       V_X=>NULL(),V_Z=>NULL(),W_X=>NULL(),W_Y=>NULL(),RHOP=>NULL(),VEL_OTHER=>NULL()
TYPE (SURFACE_TYPE), POINTER :: SF=>NULL()
TYPE (OMESH_TYPE), POINTER :: OM=>NULL()
TYPE (VENTS_TYPE), POINTER :: VT
TYPE (WALL_TYPE), POINTER :: WCM,WCP

IF (SOLID_PHASE_ONLY) RETURN
IF (PERIODIC_TEST==12) RETURN
IF (PERIODIC_TEST==13) RETURN

TNOW = CURRENT_TIME()

! Assign local names to variables

CALL POINT_TO_MESH(NM)

! Point to the appropriate velocity field

IF (PREDICTOR) THEN
   UU => US
   VV => VS
   WW => WS
   RHOP => RHOS
ELSE
   UU => U
   VV => V
   WW => W
   RHOP => RHO
ENDIF

! Set the boundary velocity place holder to some large negative number

IF (CORRECTOR) THEN
   U_Y => WORK1
   U_Z => WORK2
   V_X => WORK3
   V_Z => WORK4
   W_X => WORK5
   W_Y => WORK6
   U_Y = -1.E6_EB
   U_Z = -1.E6_EB
   V_X = -1.E6_EB
   V_Z = -1.E6_EB
   W_X = -1.E6_EB
   W_Y = -1.E6_EB
   U_EDGE_Y = -1.E6_EB
   U_EDGE_Z = -1.E6_EB
   V_EDGE_X = -1.E6_EB
   V_EDGE_Z = -1.E6_EB
   W_EDGE_X = -1.E6_EB
   W_EDGE_Y = -1.E6_EB
ENDIF

! Set OME_E and TAU_E to very negative number

TAU_E = -1.E6_EB
OME_E = -1.E6_EB

! Loop over all cell edges and determine the appropriate velocity BCs

EDGE_LOOP: DO IE=1,N_EDGES

   INTERPOLATED_EDGE = .FALSE.

   ! Throw out edges that are completely surrounded by blockages or the exterior of the domain

   PROCESS_EDGE = .FALSE.
   DO IS=5,8
      IF (.NOT.EXTERIOR(IJKE(IS,IE)) .AND. .NOT.SOLID(IJKE(IS,IE))) THEN
         PROCESS_EDGE = .TRUE.
         EXIT
      ENDIF
   ENDDO
   IF (.NOT.PROCESS_EDGE) CYCLE EDGE_LOOP

   ! If the edge is to be "smoothed," set tau and omega to zero and cycle

   IF (EVACUATION_ONLY(NM)) THEN
      OME_E(:,IE) = 0._EB
      TAU_E(:,IE) = 0._EB
      CYCLE EDGE_LOOP
   ENDIF

   ! Unpack indices for the edge

   II     = IJKE( 1,IE)
   JJ     = IJKE( 2,IE)
   KK     = IJKE( 3,IE)
   IEC    = IJKE( 4,IE)
   ICMM   = IJKE( 5,IE)
   ICPM   = IJKE( 6,IE)
   ICMP   = IJKE( 7,IE)
   ICPP   = IJKE( 8,IE)
   NOM(1) = IJKE( 9,IE)
   IIO(1) = IJKE(10,IE)
   JJO(1) = IJKE(11,IE)
   KKO(1) = IJKE(12,IE)
   NOM(2) = IJKE(13,IE)
   IIO(2) = IJKE(14,IE)
   JJO(2) = IJKE(15,IE)
   KKO(2) = IJKE(16,IE)

   ! Get the velocity components at the appropriate cell faces

   COMPONENT: SELECT CASE(IEC)
      CASE(1) COMPONENT
         UUP(1)  = VV(II,JJ,KK+1)
         UUM(1)  = VV(II,JJ,KK)
         UUP(2)  = WW(II,JJ+1,KK)
         UUM(2)  = WW(II,JJ,KK)
         DXX(1)  = DY(JJ)
         DXX(2)  = DZ(KK)
         MUA      = 0.25_EB*(MU(II,JJ,KK) + MU(II,JJ+1,KK) + MU(II,JJ+1,KK+1) + MU(II,JJ,KK+1) )
      CASE(2) COMPONENT
         UUP(1)  = WW(II+1,JJ,KK)
         UUM(1)  = WW(II,JJ,KK)
         UUP(2)  = UU(II,JJ,KK+1)
         UUM(2)  = UU(II,JJ,KK)
         DXX(1)  = DZ(KK)
         DXX(2)  = DX(II)
         MUA      = 0.25_EB*(MU(II,JJ,KK) + MU(II+1,JJ,KK) + MU(II+1,JJ,KK+1) + MU(II,JJ,KK+1) )
      CASE(3) COMPONENT
         UUP(1)  = UU(II,JJ+1,KK)
         UUM(1)  = UU(II,JJ,KK)
         UUP(2)  = VV(II+1,JJ,KK)
         UUM(2)  = VV(II,JJ,KK)
         DXX(1)  = DX(II)
         DXX(2)  = DY(JJ)
         MUA      = 0.25_EB*(MU(II,JJ,KK) + MU(II+1,JJ,KK) + MU(II+1,JJ+1,KK) + MU(II,JJ+1,KK) )
   END SELECT COMPONENT

   ! Indicate that the velocity gradients in the two orthogonal directions have not been changed yet

   ALTERED_GRADIENT = .FALSE.

   ! Loop over all possible orientations of edge and reassign velocity gradients if appropriate
   
   SIGN_LOOP: DO I_SGN=-1,1,2
      ORIENTATION_LOOP: DO IS=1,3

         IF (IS==IEC) CYCLE ORIENTATION_LOOP

         ! IOR is the orientation of the wall cells adjacent to the edge

         IOR = I_SGN*IS

         ! IS2 is the other coordinate direction besides IOR.

         SELECT CASE(IEC)
            CASE(1)
               IF (IS==2) IS2 = 3
               IF (IS==3) IS2 = 2
            CASE(2)
               IF (IS==1) IS2 = 3
               IF (IS==3) IS2 = 1
            CASE(3)
               IF (IS==1) IS2 = 2
               IF (IS==2) IS2 = 1
            END SELECT

         ! Determine Index_Coordinate_Direction
         ! IEC=1, ICD=1 refers to DWDY; ICD=2 refers to DVDZ
         ! IEC=2, ICD=1 refers to DUDZ; ICD=2 refers to DWDX
         ! IEC=3, ICD=1 refers to DVDX; ICD=2 refers to DUDY

         IF (IS>IEC) ICD = IS-IEC
         IF (IS<IEC) ICD = IS-IEC+3
         ICD_SGN = I_SGN * ICD

         ! IWM and IWP are the wall cell indices of the boundary on either side of the edge.

         IF (IOR<0) THEN
            IWM  = WALL_INDEX(ICMM,-IOR)
            IWMI = WALL_INDEX(ICMM,IS2)
            IF (ICD==1) THEN
               IWP  = WALL_INDEX(ICMP,-IOR)
               IWPI = WALL_INDEX(ICMP,-IS2)
            ELSE ! ICD==2
               IWP  = WALL_INDEX(ICPM,-IOR)
               IWPI = WALL_INDEX(ICPM,-IS2)
            ENDIF
         ELSE
            IF (ICD==1) THEN
               IWM  = WALL_INDEX(ICPM,-IOR)
               IWMI = WALL_INDEX(ICPM,IS2)
            ELSE ! ICD==2
               IWM  = WALL_INDEX(ICMP,-IOR)
               IWMI = WALL_INDEX(ICMP,IS2)
            ENDIF
            IWP  = WALL_INDEX(ICPP,-IOR)
            IWPI = WALL_INDEX(ICPP,-IS2)
         ENDIF

         ! If both adjacent wall cells are undefined, cycle out of the loop.

         IF (IWM==0 .AND. IWP==0) CYCLE ORIENTATION_LOOP

         ! If there is a solid wall separating the two adjacent wall cells, cycle out of the loop.

         IF (WALL(IWMI)%BOUNDARY_TYPE==SOLID_BOUNDARY .OR. WALL(IWPI)%BOUNDARY_TYPE==SOLID_BOUNDARY) CYCLE ORIENTATION_LOOP

         ! If only one adjacent wall cell is defined, use its properties.

         IF (IWM>0) THEN
            WCM => WALL(IWM)
         ELSE
            WCM => WALL(IWP)
         ENDIF

         IF (IWP>0) THEN
            WCP => WALL(IWP)
         ELSE
            WCP => WALL(IWM)
         ENDIF

         ! If both adjacent wall cells are NULL, cycle out.

         BOUNDARY_TYPE_M = WCM%BOUNDARY_TYPE
         BOUNDARY_TYPE_P = WCP%BOUNDARY_TYPE

         IF (BOUNDARY_TYPE_M==NULL_BOUNDARY .AND. BOUNDARY_TYPE_P==NULL_BOUNDARY) CYCLE ORIENTATION_LOOP

         ! OPEN boundary conditions, both varieties, with and without a wind

         OPEN_AND_WIND_BC: IF ((IWM==0.OR.WALL(IWM)%BOUNDARY_TYPE==OPEN_BOUNDARY) .AND. &
                               (IWP==0.OR.WALL(IWP)%BOUNDARY_TYPE==OPEN_BOUNDARY)) THEN

            VENT_INDEX = MAX(WCM%VENT_INDEX,WCP%VENT_INDEX)
            VT => VENTS(VENT_INDEX)

            WIND_NO_WIND_IF: IF (.NOT.ANY(MEAN_FORCING)) THEN  ! For regular OPEN boundary, (free-slip) BCs

               SELECT CASE(IEC)
                  CASE(1)
                     IF (JJ==0    .AND. IOR== 2) WW(II,0,KK)    = WW(II,1,KK)
                     IF (JJ==JBAR .AND. IOR==-2) WW(II,JBP1,KK) = WW(II,JBAR,KK)
                     IF (KK==0    .AND. IOR== 3) VV(II,JJ,0)    = VV(II,JJ,1)
                     IF (KK==KBAR .AND. IOR==-3) VV(II,JJ,KBP1) = VV(II,JJ,KBAR)
                  CASE(2)
                     IF (II==0    .AND. IOR== 1) WW(0,JJ,KK)    = WW(1,JJ,KK)
                     IF (II==IBAR .AND. IOR==-1) WW(IBP1,JJ,KK) = WW(IBAR,JJ,KK)
                     IF (KK==0    .AND. IOR== 3) UU(II,JJ,0)    = UU(II,JJ,1)
                     IF (KK==KBAR .AND. IOR==-3) UU(II,JJ,KBP1) = UU(II,JJ,KBAR)
                  CASE(3)
                     IF (II==0    .AND. IOR== 1) VV(0,JJ,KK)    = VV(1,JJ,KK)
                     IF (II==IBAR .AND. IOR==-1) VV(IBP1,JJ,KK) = VV(IBAR,JJ,KK)
                     IF (JJ==0    .AND. IOR== 2) UU(II,0,KK)    = UU(II,1,KK)
                     IF (JJ==JBAR .AND. IOR==-2) UU(II,JBP1,KK) = UU(II,JBAR,KK)
               END SELECT

            ELSE WIND_NO_WIND_IF  ! For wind, use prescribed far-field velocity all around

               UBAR = U0*EVALUATE_RAMP(T,DUMMY,I_RAMP_U0_T)*EVALUATE_RAMP(ZC(KK),DUMMY,I_RAMP_U0_Z)
               VBAR = V0*EVALUATE_RAMP(T,DUMMY,I_RAMP_V0_T)*EVALUATE_RAMP(ZC(KK),DUMMY,I_RAMP_V0_Z)
               WBAR = W0*EVALUATE_RAMP(T,DUMMY,I_RAMP_W0_T)*EVALUATE_RAMP(ZC(KK),DUMMY,I_RAMP_W0_Z)

               SELECT CASE(IEC)
                  CASE(1)
                     IF (JJ==0    .AND. IOR== 2) WW(II,0,KK)    = WBAR
                     IF (JJ==JBAR .AND. IOR==-2) WW(II,JBP1,KK) = WBAR
                     IF (KK==0    .AND. IOR== 3) VV(II,JJ,0)    = VBAR
                     IF (KK==KBAR .AND. IOR==-3) VV(II,JJ,KBP1) = VBAR
                  CASE(2)
                     IF (II==0    .AND. IOR== 1) WW(0,JJ,KK)    = WBAR
                     IF (II==IBAR .AND. IOR==-1) WW(IBP1,JJ,KK) = WBAR
                     IF (KK==0    .AND. IOR== 3) UU(II,JJ,0)    = UBAR
                     IF (KK==KBAR .AND. IOR==-3) UU(II,JJ,KBP1) = UBAR
                  CASE(3)
                     IF (II==0    .AND. IOR== 1) VV(0,JJ,KK)    = VBAR
                     IF (II==IBAR .AND. IOR==-1) VV(IBP1,JJ,KK) = VBAR
                     IF (JJ==0    .AND. IOR== 2) UU(II,0,KK)    = UBAR
                     IF (JJ==JBAR .AND. IOR==-2) UU(II,JBP1,KK) = UBAR
               END SELECT

            ENDIF WIND_NO_WIND_IF

            IF (IWM/=0 .AND. IWP/=0) THEN
               CYCLE EDGE_LOOP  ! Do no further processing of this edge if both cell faces are OPEN
            ELSE
               CYCLE ORIENTATION_LOOP
            ENDIF

         ENDIF OPEN_AND_WIND_BC

         ! Define the appropriate gas and ghost velocity

         IF (ICD==1) THEN ! Used to pick the appropriate velocity component
            IVL=2
         ELSE !ICD==2
            IVL=1
         ENDIF

         IF (IOR<0) THEN
            VEL_GAS   = UUM(IVL)
            VEL_GHOST = UUP(IVL)
            IIGM = I_CELL(ICMM)
            JJGM = J_CELL(ICMM)
            KKGM = K_CELL(ICMM)
            IF (ICD==1) THEN
               IIGP = I_CELL(ICMP)
               JJGP = J_CELL(ICMP)
               KKGP = K_CELL(ICMP)
            ELSE ! ICD==2
               IIGP = I_CELL(ICPM)
               JJGP = J_CELL(ICPM)
               KKGP = K_CELL(ICPM)
            ENDIF
         ELSE
            VEL_GAS   = UUP(IVL)
            VEL_GHOST = UUM(IVL)
            IF (ICD==1) THEN
               IIGM = I_CELL(ICPM)
               JJGM = J_CELL(ICPM)
               KKGM = K_CELL(ICPM)
            ELSE ! ICD==2
               IIGM = I_CELL(ICMP)
               JJGM = J_CELL(ICMP)
               KKGM = K_CELL(ICMP)
            ENDIF
            IIGP = I_CELL(ICPP)
            JJGP = J_CELL(ICPP)
            KKGP = K_CELL(ICPP)
         ENDIF

         ! Decide whether or not to process edge using data interpolated from another mesh

         INTERPOLATION_IF: IF (NOM(ICD)==0 .OR. &
                   (BOUNDARY_TYPE_M==SOLID_BOUNDARY .OR. BOUNDARY_TYPE_P==SOLID_BOUNDARY) .OR. &
                   (BOUNDARY_TYPE_M/=INTERPOLATED_BOUNDARY .AND. BOUNDARY_TYPE_P/=INTERPOLATED_BOUNDARY)) THEN

            ! Determine appropriate velocity BC by assessing each adjacent wall cell. If the BCs are different on each
            ! side of the edge, choose the one with the specified velocity or velocity gradient, if there is one.
            ! If not, choose the max value of boundary condition index, simply for consistency.

            SURF_INDEXM = WCM%SURF_INDEX
            SURF_INDEXP = WCP%SURF_INDEX
            IF (SURFACE(SURF_INDEXM)%SPECIFIED_NORMAL_VELOCITY .OR. SURFACE(SURF_INDEXM)%SPECIFIED_NORMAL_GRADIENT) THEN
               SF=>SURFACE(SURF_INDEXM)
            ELSEIF (SURFACE(SURF_INDEXP)%SPECIFIED_NORMAL_VELOCITY .OR. SURFACE(SURF_INDEXP)%SPECIFIED_NORMAL_GRADIENT) THEN
               SF=>SURFACE(SURF_INDEXP)
            ELSE
               SF=>SURFACE(MAX(SURF_INDEXM,SURF_INDEXP))
            ENDIF
            VELOCITY_BC_INDEX = SF%VELOCITY_BC_INDEX
            IF (WCM%VENT_INDEX==WCP%VENT_INDEX .AND. WCP%VENT_INDEX > 0) THEN
               IF(VENTS(WCM%VENT_INDEX)%NODE_INDEX>0 .AND. WCM%ONE_D%U_NORMAL >= 0._EB) VELOCITY_BC_INDEX=FREE_SLIP_BC
            ENDIF

            ! Compute the viscosity in the two adjacent gas cells

            MUA = 0.5_EB*(MU(IIGM,JJGM,KKGM) + MU(IIGP,JJGP,KKGP))

            ! Set up synthetic eddy method (experimental)

            SYNTHETIC_EDDY_METHOD = .FALSE.
            HVAC_TANGENTIAL = .FALSE.
            IF (IWM>0 .AND. IWP>0) THEN
               IF (WCM%VENT_INDEX==WCP%VENT_INDEX) THEN
                  IF (WCM%VENT_INDEX>0) THEN
                     VT=>VENTS(WCM%VENT_INDEX)
                     IF (VT%N_EDDY>0) SYNTHETIC_EDDY_METHOD=.TRUE.
                     IF (ALL(VT%UVW > -1.E12_EB) .AND. VT%NODE_INDEX > 0) HVAC_TANGENTIAL = .TRUE.
                  ENDIF
               ENDIF
            ENDIF

            ! Determine if there is a tangential velocity component

            VEL_T_IF: IF (.NOT.SF%SPECIFIED_TANGENTIAL_VELOCITY .AND. .NOT.SYNTHETIC_EDDY_METHOD .AND. .NOT. HVAC_TANGENTIAL) THEN
               VEL_T = 0._EB
            ELSE VEL_T_IF
               VEL_EDDY = 0._EB
               SYNTHETIC_EDDY_IF: IF (SYNTHETIC_EDDY_METHOD) THEN
                  IS_SELECT: SELECT CASE(IS) ! unsigned vent orientation
                     CASE(1) ! yz plane
                        SELECT CASE(IEC) ! edge orientation
                           CASE(2)
                              IF (ICD==1) VEL_EDDY = 0.5_EB*(VT%U_EDDY(JJ,KK)+VT%U_EDDY(JJ,KK+1))
                              IF (ICD==2) VEL_EDDY = 0.5_EB*(VT%W_EDDY(JJ,KK)+VT%W_EDDY(JJ,KK+1))
                           CASE(3)
                              IF (ICD==1) VEL_EDDY = 0.5_EB*(VT%V_EDDY(JJ,KK)+VT%V_EDDY(JJ+1,KK))
                              IF (ICD==2) VEL_EDDY = 0.5_EB*(VT%U_EDDY(JJ,KK)+VT%U_EDDY(JJ+1,KK))
                        END SELECT
                     CASE(2) ! zx plane
                        SELECT CASE(IEC)
                           CASE(3)
                              IF (ICD==1) VEL_EDDY = 0.5_EB*(VT%V_EDDY(II,KK)+VT%V_EDDY(II+1,KK))
                              IF (ICD==2) VEL_EDDY = 0.5_EB*(VT%U_EDDY(II,KK)+VT%U_EDDY(II+1,KK))
                           CASE(1)
                              IF (ICD==1) VEL_EDDY = 0.5_EB*(VT%W_EDDY(II,KK)+VT%W_EDDY(II,KK+1))
                              IF (ICD==2) VEL_EDDY = 0.5_EB*(VT%V_EDDY(II,KK)+VT%V_EDDY(II,KK+1))
                        END SELECT
                     CASE(3) ! xy plane
                        SELECT CASE(IEC)
                           CASE(1)
                              IF (ICD==1) VEL_EDDY = 0.5_EB*(VT%W_EDDY(II,JJ)+VT%W_EDDY(II,JJ+1))
                              IF (ICD==2) VEL_EDDY = 0.5_EB*(VT%V_EDDY(II,JJ)+VT%V_EDDY(II,JJ+1))
                           CASE(2)
                              IF (ICD==1) VEL_EDDY = 0.5_EB*(VT%U_EDDY(II,JJ)+VT%U_EDDY(II+1,JJ))
                              IF (ICD==2) VEL_EDDY = 0.5_EB*(VT%W_EDDY(II,JJ)+VT%W_EDDY(II+1,JJ))
                        END SELECT
                  END SELECT IS_SELECT
               ENDIF SYNTHETIC_EDDY_IF
               IF (ABS(SF%T_IGN-T_BEGIN)<=SPACING(SF%T_IGN) .AND. SF%RAMP_INDEX(TIME_VELO)>=1) THEN
                  TSI = T
               ELSE
                  TSI=T-SF%T_IGN
               ENDIF
               PROFILE_FACTOR = 1._EB
               IF (HVAC_TANGENTIAL .AND. 0.5_EB*(WCM%ONE_D%U_NORMAL_S+WCP%ONE_D%U_NORMAL_S) > 0._EB) HVAC_TANGENTIAL = .FALSE.
               IF (HVAC_TANGENTIAL) THEN
                  VEL_T = 0._EB
                  IEC_SELECT: SELECT CASE(IEC) ! edge orientation
                     CASE (1)
                        IF (ICD==1) VEL_T = 0.5_EB*ABS((WCM%ONE_D%U_NORMAL_S+WCP%ONE_D%U_NORMAL_S)/VT%UVW(ABS(VT%IOR)))*VT%UVW(3)
                        IF (ICD==2) VEL_T = 0.5_EB*ABS((WCM%ONE_D%U_NORMAL_S+WCP%ONE_D%U_NORMAL_S)/VT%UVW(ABS(VT%IOR)))*VT%UVW(2)
                     CASE (2)
                        IF (ICD==1) VEL_T = 0.5_EB*ABS((WCM%ONE_D%U_NORMAL_S+WCP%ONE_D%U_NORMAL_S)/VT%UVW(ABS(VT%IOR)))*VT%UVW(1)
                        IF (ICD==2) VEL_T = 0.5_EB*ABS((WCM%ONE_D%U_NORMAL_S+WCP%ONE_D%U_NORMAL_S)/VT%UVW(ABS(VT%IOR)))*VT%UVW(3)
                     CASE (3)
                        IF (ICD==1) VEL_T = 0.5_EB*ABS((WCM%ONE_D%U_NORMAL_S+WCP%ONE_D%U_NORMAL_S)/VT%UVW(ABS(VT%IOR)))*VT%UVW(2)
                        IF (ICD==2) VEL_T = 0.5_EB*ABS((WCM%ONE_D%U_NORMAL_S+WCP%ONE_D%U_NORMAL_S)/VT%UVW(ABS(VT%IOR)))*VT%UVW(1)
                  END SELECT IEC_SELECT
               ELSE
                  IF (SF%PROFILE/=0 .AND. SF%VEL>TWO_EPSILON_EB) &
                     PROFILE_FACTOR = ABS(0.5_EB*(WCM%ONE_D%U_NORMAL_0+WCP%ONE_D%U_NORMAL_0)/SF%VEL)
                  RAMP_T = EVALUATE_RAMP(TSI,SF%TAU(TIME_VELO),SF%RAMP_INDEX(TIME_VELO))
                  IF (IEC==1 .OR. (IEC==2 .AND. ICD==2)) VEL_T = RAMP_T*(PROFILE_FACTOR*(SF%VEL_T(2) + VEL_EDDY))
                  IF (IEC==3 .OR. (IEC==2 .AND. ICD==1)) VEL_T = RAMP_T*(PROFILE_FACTOR*(SF%VEL_T(1) + VEL_EDDY))
               ENDIF
            ENDIF VEL_T_IF

            ! Choose the appropriate boundary condition to apply

            HVAC_IF: IF (HVAC_TANGENTIAL)  THEN

               VEL_GHOST = 2._EB*VEL_T - VEL_GAS
               DUIDXJ(ICD_SGN) = I_SGN*(VEL_GAS-VEL_GHOST)/DXX(ICD)
               MU_DUIDXJ(ICD_SGN) = MUA*DUIDXJ(ICD_SGN)
               ALTERED_GRADIENT(ICD_SGN) = .TRUE.

            ELSE HVAC_IF

               BOUNDARY_CONDITION: SELECT CASE(VELOCITY_BC_INDEX)

                  CASE (FREE_SLIP_BC) BOUNDARY_CONDITION

                     VEL_GHOST = VEL_GAS
                     DUIDXJ(ICD_SGN) = I_SGN*(VEL_GAS-VEL_GHOST)/DXX(ICD)
                     MU_DUIDXJ(ICD_SGN) = MUA*DUIDXJ(ICD_SGN)
                     ALTERED_GRADIENT(ICD_SGN) = .TRUE.

                  CASE (NO_SLIP_BC) BOUNDARY_CONDITION

                     WANNIER_BC: IF (PERIODIC_TEST==5) THEN
                        SELECT CASE(IOR)
                           CASE( 1)
                              VEL_T = WANNIER_FLOW(X(II),Z(KK),2)
                           CASE(-1)
                              VEL_T = WANNIER_FLOW(X(II),Z(KK),2)
                           CASE( 3)
                              VEL_T = WANNIER_FLOW(X(II),Z(KK),1)
                           CASE(-3)
                              VEL_T = WANNIER_FLOW(X(II),Z(KK),1)
                        END SELECT
                     ENDIF WANNIER_BC

                     VEL_GHOST = 2._EB*VEL_T - VEL_GAS
                     DUIDXJ(ICD_SGN) = I_SGN*(VEL_GAS-VEL_GHOST)/DXX(ICD)
                     MU_DUIDXJ(ICD_SGN) = MUA*DUIDXJ(ICD_SGN)
                     ALTERED_GRADIENT(ICD_SGN) = .TRUE.

                  CASE (WALL_MODEL_BC) BOUNDARY_CONDITION

                     ITMP = MIN(5000,NINT(0.5_EB*(TMP(IIGM,JJGM,KKGM)+TMP(IIGP,JJGP,KKGP))))
                     MU_WALL = MU_RSQMW_Z(ITMP,1)/RSQ_MW_Z(1)
                     RHO_WALL = 0.5_EB*( RHOP(IIGM,JJGM,KKGM) + RHOP(IIGP,JJGP,KKGP) )
                     CALL WALL_MODEL(SLIP_COEF,U_TAU,Y_PLUS,VEL_GAS-VEL_T,MU_WALL/RHO_WALL,DXX(ICD),SF%ROUGHNESS)
                     SELECT CASE(SLIP_CONDITION)
                        CASE(0)
                           SLIP_COEF = -1._EB
                        CASE(1)
                           SLIP_COEF = SLIP_COEF
                        CASE(2)
                           SLIP_COEF = 0.5_EB*(SLIP_COEF-1._EB)
                        CASE(3)
                           WT1 = MAX(0._EB,MIN(1._EB,(Y_PLUS-Y_WERNER_WENGLE)/(Y_PLUS+TWO_EPSILON_EB)))
                           WT2 = 1._EB-WT1
                           SLIP_COEF = WT1*SLIP_COEF-WT2
                        CASE(4)
                           IF ( ABS(0.5_EB*(WCM%ONE_D%U_NORMAL_S+WCP%ONE_D%U_NORMAL_S))>ABS(VEL_GAS-VEL_T) ) THEN
                              SLIP_COEF = -1._EB
                           ELSE
                              SLIP_COEF = 0.5_EB*(SLIP_COEF-1._EB)
                           ENDIF
                     END SELECT
                     VEL_GHOST = VEL_T + SLIP_COEF*(VEL_GAS-VEL_T)
                     DUIDXJ(ICD_SGN) = I_SGN*(VEL_GAS-VEL_GHOST)/DXX(ICD)
                     MU_DUIDXJ(ICD_SGN) = RHO_WALL*U_TAU**2 * SIGN(1._EB,I_SGN*(VEL_GAS-VEL_T))
                     ALTERED_GRADIENT(ICD_SGN) = .TRUE.

                  CASE (BOUNDARY_FUEL_MODEL_BC) BOUNDARY_CONDITION

                     RHO_WALL = 0.5_EB*( RHOP(IIGM,JJGM,KKGM) + RHOP(IIGP,JJGP,KKGP) )
                     VEL_T = SQRT(UU(IIGM,JJGM,KKGM)**2 + VV(IIGM,JJGM,KKGM)**2)
                     VEL_GHOST = VEL_GAS
                     DUIDXJ(ICD_SGN) = I_SGN*(VEL_GAS-VEL_GHOST)/DXX(ICD)
                     MU_DUIDXJ(ICD_SGN) = I_SGN*0.5_EB*RHO_WALL*SF%DRAG_COEFFICIENT*SF%SHAPE_FACTOR*SF%LAYER_THICKNESS(1)*&
                                          SF%PACKING_RATIO(1)*SF%SURFACE_VOLUME_RATIO(1)*VEL_GAS*VEL_T
                     ALTERED_GRADIENT(ICD_SGN) = .TRUE.

               END SELECT BOUNDARY_CONDITION

            ENDIF HVAC_IF

         ELSE INTERPOLATION_IF  ! Use data from another mesh

            INTERPOLATED_EDGE = .TRUE.
            OM => OMESH(ABS(NOM(ICD)))

            IF (PREDICTOR) THEN
               SELECT CASE(IEC)
                  CASE(1)
                     IF (ICD==1) THEN
                        VEL_OTHER => OM%WS
                     ELSE ! ICD=2
                        VEL_OTHER => OM%VS
                     ENDIF
                  CASE(2)
                     IF (ICD==1) THEN
                        VEL_OTHER => OM%US
                     ELSE ! ICD=2
                        VEL_OTHER => OM%WS
                     ENDIF
                  CASE(3)
                     IF (ICD==1) THEN
                        VEL_OTHER => OM%VS
                     ELSE ! ICD=2
                        VEL_OTHER => OM%US
                     ENDIF
               END SELECT
            ELSE
               SELECT CASE(IEC)
                  CASE(1)
                     IF (ICD==1) THEN
                        VEL_OTHER => OM%W
                     ELSE ! ICD=2
                        VEL_OTHER => OM%V
                     ENDIF
                  CASE(2)
                     IF (ICD==1) THEN
                        VEL_OTHER => OM%U
                     ELSE ! ICD=2
                        VEL_OTHER => OM%W
                     ENDIF
                  CASE(3)
                     IF (ICD==1) THEN
                        VEL_OTHER => OM%V
                     ELSE ! ICD=2
                        VEL_OTHER => OM%U
                     ENDIF
               END SELECT
            ENDIF

            WGT = EDGE_INTERPOLATION_FACTOR(IE,ICD)
            OMW = 1._EB-WGT

            SELECT CASE(IEC)
               CASE(1)
                  IF (ICD==1) THEN
                     VEL_GHOST = WGT*VEL_OTHER(IIO(ICD),JJO(ICD),KKO(ICD)) + OMW*VEL_OTHER(IIO(ICD),JJO(ICD),KKO(ICD)-1)
                  ELSE ! ICD=2
                     VEL_GHOST = WGT*VEL_OTHER(IIO(ICD),JJO(ICD),KKO(ICD)) + OMW*VEL_OTHER(IIO(ICD),JJO(ICD)-1,KKO(ICD))
                  ENDIF
               CASE(2)
                  IF (ICD==1) THEN
                     VEL_GHOST = WGT*VEL_OTHER(IIO(ICD),JJO(ICD),KKO(ICD)) + OMW*VEL_OTHER(IIO(ICD)-1,JJO(ICD),KKO(ICD))
                  ELSE ! ICD=2
                     VEL_GHOST = WGT*VEL_OTHER(IIO(ICD),JJO(ICD),KKO(ICD)) + OMW*VEL_OTHER(IIO(ICD),JJO(ICD),KKO(ICD)-1)
                  ENDIF
               CASE(3)
                  IF (ICD==1) THEN
                     VEL_GHOST = WGT*VEL_OTHER(IIO(ICD),JJO(ICD),KKO(ICD)) + OMW*VEL_OTHER(IIO(ICD),JJO(ICD)-1,KKO(ICD))
                  ELSE ! ICD==2
                     VEL_GHOST = WGT*VEL_OTHER(IIO(ICD),JJO(ICD),KKO(ICD)) + OMW*VEL_OTHER(IIO(ICD)-1,JJO(ICD),KKO(ICD))
                  ENDIF
            END SELECT

         ENDIF INTERPOLATION_IF

         ! Set ghost cell values at edge of computational domain

         SELECT CASE(IEC)
            CASE(1)
               IF (JJ==0    .AND. IOR== 2) WW(II,JJ,KK)   = VEL_GHOST
               IF (JJ==JBAR .AND. IOR==-2) WW(II,JJ+1,KK) = VEL_GHOST
               IF (KK==0    .AND. IOR== 3) VV(II,JJ,KK)   = VEL_GHOST
               IF (KK==KBAR .AND. IOR==-3) VV(II,JJ,KK+1) = VEL_GHOST
               IF (CORRECTOR .AND. JJ>0 .AND. JJ<JBAR .AND. KK>0 .AND. KK<KBAR) THEN
                 IF (ICD==1) THEN
                    W_Y(II,JJ,KK) = 0.5_EB*(VEL_GHOST+VEL_GAS)
                 ELSE ! ICD=2
                    V_Z(II,JJ,KK) = 0.5_EB*(VEL_GHOST+VEL_GAS)
                 ENDIF
               ENDIF
            CASE(2)
               IF (II==0    .AND. IOR== 1) WW(II,JJ,KK)   = VEL_GHOST
               IF (II==IBAR .AND. IOR==-1) WW(II+1,JJ,KK) = VEL_GHOST
               IF (KK==0    .AND. IOR== 3) UU(II,JJ,KK)   = VEL_GHOST
               IF (KK==KBAR .AND. IOR==-3) UU(II,JJ,KK+1) = VEL_GHOST
               IF (CORRECTOR .AND. II>0 .AND. II<IBAR .AND. KK>0 .AND. KK<KBAR) THEN
                 IF (ICD==1) THEN
                    U_Z(II,JJ,KK) = 0.5_EB*(VEL_GHOST+VEL_GAS)
                 ELSE ! ICD=2
                    W_X(II,JJ,KK) = 0.5_EB*(VEL_GHOST+VEL_GAS)
                 ENDIF
               ENDIF
            CASE(3)
               IF (II==0    .AND. IOR== 1) VV(II,JJ,KK)   = VEL_GHOST
               IF (II==IBAR .AND. IOR==-1) VV(II+1,JJ,KK) = VEL_GHOST
               IF (JJ==0    .AND. IOR== 2) UU(II,JJ,KK)   = VEL_GHOST
               IF (JJ==JBAR .AND. IOR==-2) UU(II,JJ+1,KK) = VEL_GHOST
               IF (CORRECTOR .AND. II>0 .AND. II<IBAR .AND. JJ>0 .AND. JJ<JBAR) THEN
                 IF (ICD==1) THEN
                    V_X(II,JJ,KK) = 0.5_EB*(VEL_GHOST+VEL_GAS)
                 ELSE ! ICD=2
                    U_Y(II,JJ,KK) = 0.5_EB*(VEL_GHOST+VEL_GAS)
                 ENDIF
               ENDIF
         END SELECT

      ENDDO ORIENTATION_LOOP
   ENDDO SIGN_LOOP

   ! Cycle out of the EDGE_LOOP if no tangential gradients have been altered.

   IF (.NOT.ANY(ALTERED_GRADIENT)) CYCLE EDGE_LOOP

   ! If the edge is on an interpolated boundary, and all cells around it are not solid, cycle

   IF (INTERPOLATED_EDGE) THEN
      PROCESS_EDGE = .FALSE.
      DO IS=5,8
         IF (SOLID(IJKE(IS,IE))) PROCESS_EDGE = .TRUE.
      ENDDO
      IF (.NOT.PROCESS_EDGE) CYCLE EDGE_LOOP
   ENDIF

   ! Loop over all 4 normal directions and compute vorticity and stress tensor components for each

   SIGN_LOOP_2: DO I_SGN=-1,1,2
      ORIENTATION_LOOP_2: DO ICD=1,2
         IF (ICD==1) THEN
            ICDO=2
         ELSE ! ICD=2
            ICDO=1
         ENDIF
         ICD_SGN = I_SGN*ICD
         IF (ALTERED_GRADIENT(ICD_SGN)) THEN
               DUIDXJ_USE(ICD) =    DUIDXJ(ICD_SGN)
            MU_DUIDXJ_USE(ICD) = MU_DUIDXJ(ICD_SGN)
         ELSEIF (ALTERED_GRADIENT(-ICD_SGN)) THEN
               DUIDXJ_USE(ICD) =    DUIDXJ(-ICD_SGN)
            MU_DUIDXJ_USE(ICD) = MU_DUIDXJ(-ICD_SGN)
         ELSE
            CYCLE ORIENTATION_LOOP_2
         ENDIF
         ICDO_SGN = I_SGN*ICDO
         IF (ALTERED_GRADIENT(ICDO_SGN)) THEN
               DUIDXJ_USE(ICDO) =    DUIDXJ(ICDO_SGN)
            MU_DUIDXJ_USE(ICDO) = MU_DUIDXJ(ICDO_SGN)
         ELSEIF (ALTERED_GRADIENT(-ICDO_SGN)) THEN
               DUIDXJ_USE(ICDO) =    DUIDXJ(-ICDO_SGN)
            MU_DUIDXJ_USE(ICDO) = MU_DUIDXJ(-ICDO_SGN)
         ELSE
               DUIDXJ_USE(ICDO) = 0._EB
            MU_DUIDXJ_USE(ICDO) = 0._EB
         ENDIF
         OME_E(ICD_SGN,IE) =    DUIDXJ_USE(1) -    DUIDXJ_USE(2)
         TAU_E(ICD_SGN,IE) = MU_DUIDXJ_USE(1) + MU_DUIDXJ_USE(2)
      ENDDO ORIENTATION_LOOP_2
   ENDDO SIGN_LOOP_2

ENDDO EDGE_LOOP

! Store cell edge velocity averages of the velocity components for use in Smokeview only

IF (CORRECTOR) THEN
   DO K=0,KBAR
      DO J=0,JBAR
         DO I=0,IBAR
            IC = CELL_INDEX(I,J,K)
            IF (IC==0) CYCLE
            IF (U_Y(I,J,K)>-1.E5_EB) U_EDGE_Y(IC) = U_Y(I,J,K)
            IF (U_Z(I,J,K)>-1.E5_EB) U_EDGE_Z(IC) = U_Z(I,J,K)
            IF (V_X(I,J,K)>-1.E5_EB) V_EDGE_X(IC) = V_X(I,J,K)
            IF (V_Z(I,J,K)>-1.E5_EB) V_EDGE_Z(IC) = V_Z(I,J,K)
            IF (W_X(I,J,K)>-1.E5_EB) W_EDGE_X(IC) = W_X(I,J,K)
            IF (W_Y(I,J,K)>-1.E5_EB) W_EDGE_Y(IC) = W_Y(I,J,K)
         ENDDO
      ENDDO
   ENDDO
ENDIF

T_USED(4)=T_USED(4)+CURRENT_TIME()-TNOW
END SUBROUTINE VELOCITY_BC


SUBROUTINE MATCH_VELOCITY(NM)

! Force normal component of velocity to match at interpolated boundaries

INTEGER  :: NOM,II,JJ,KK,IOR,IW,IIO,JJO,KKO
INTEGER, INTENT(IN) :: NM
REAL(EB) :: UU_AVG,VV_AVG,WW_AVG,TNOW,DA_OTHER,UU_OTHER,VV_OTHER,WW_OTHER,NOM_CELLS
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU=>NULL(),VV=>NULL(),WW=>NULL(),OM_UU=>NULL(),OM_VV=>NULL(),OM_WW=>NULL()
TYPE (OMESH_TYPE), POINTER :: OM=>NULL()
TYPE (MESH_TYPE), POINTER :: M2=>NULL()
TYPE (WALL_TYPE), POINTER :: WC=>NULL()
TYPE (EXTERNAL_WALL_TYPE), POINTER :: EWC=>NULL()
IF (SOLID_PHASE_ONLY) RETURN
IF (EVACUATION_ONLY(NM)) RETURN

TNOW = CURRENT_TIME()

! Assign local variable names

CALL POINT_TO_MESH(NM)

! Point to the appropriate velocity field

IF (PREDICTOR) THEN
   UU => US
   VV => VS
   WW => WS
   D_CORR = 0._EB
ELSE
   UU => U
   VV => V
   WW => W
   DS_CORR = 0._EB
ENDIF

! Loop over all external wall cells and force adjacent normal components of velocty at interpolated boundaries to match.

EXTERNAL_WALL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS

   WC=>WALL(IW)

   IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY) CYCLE EXTERNAL_WALL_LOOP

   EWC=>EXTERNAL_WALL(IW)

   II  = WC%ONE_D%II
   JJ  = WC%ONE_D%JJ
   KK  = WC%ONE_D%KK
   IOR = WC%ONE_D%IOR
   NOM = EWC%NOM
   OM => OMESH(NOM)
   M2 => MESHES(NOM)

   ! Determine the area of the interpolated cell face

   DA_OTHER = 0._EB

   SELECT CASE(ABS(IOR))
      CASE(1)
         IF (PREDICTOR) OM_UU => OM%US
         IF (CORRECTOR) OM_UU => OM%U
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  DA_OTHER = DA_OTHER + M2%DY(JJO)*M2%DZ(KKO)
               ENDDO
            ENDDO
         ENDDO
      CASE(2)
         IF (PREDICTOR) OM_VV => OM%VS
         IF (CORRECTOR) OM_VV => OM%V
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  DA_OTHER = DA_OTHER + M2%DX(IIO)*M2%DZ(KKO)
               ENDDO
            ENDDO
         ENDDO
      CASE(3)
         IF (PREDICTOR) OM_WW => OM%WS
         IF (CORRECTOR) OM_WW => OM%W
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  DA_OTHER = DA_OTHER + M2%DX(IIO)*M2%DY(JJO)
               ENDDO
            ENDDO
         ENDDO
   END SELECT

   ! Determine the normal component of velocity from the other mesh and use it for average

   SELECT CASE(IOR)

      CASE( 1)

         UU_OTHER = 0._EB
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  UU_OTHER = UU_OTHER + OM_UU(IIO,JJO,KKO)*M2%DY(JJO)*M2%DZ(KKO)/DA_OTHER
                  IF (EWC%AREA_RATIO>0.9_EB) OM_UU(IIO,JJO,KKO) = 0.5_EB*(OM_UU(IIO,JJO,KKO)+UU(0,JJ,KK))
               ENDDO
            ENDDO
         ENDDO
         UU_AVG = 0.5_EB*(UU(0,JJ,KK) + UU_OTHER)
         IF (PREDICTOR) D_CORR(IW) = DS_CORR(IW) + 0.5*(UU_AVG-UU(0,JJ,KK))*R(0)*RDX(1)*RRN(1)
         IF (CORRECTOR) DS_CORR(IW) = (UU_AVG-UU(0,JJ,KK))*R(0)*RDX(1)*RRN(1)
         UVW_SAVE(IW) = UU(0,JJ,KK)
         UU(0,JJ,KK)  = UU_AVG

      CASE(-1)

         UU_OTHER = 0._EB
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  UU_OTHER = UU_OTHER + OM_UU(IIO-1,JJO,KKO)*M2%DY(JJO)*M2%DZ(KKO)/DA_OTHER
                  IF (EWC%AREA_RATIO>0.9_EB) OM_UU(IIO-1,JJO,KKO) = 0.5_EB*(OM_UU(IIO-1,JJO,KKO)+UU(IBAR,JJ,KK))
               ENDDO
            ENDDO
         ENDDO
         UU_AVG = 0.5_EB*(UU(IBAR,JJ,KK) + UU_OTHER)
         IF (PREDICTOR) D_CORR(IW) = DS_CORR(IW) - 0.5*(UU_AVG-UU(IBAR,JJ,KK))*R(IBAR)*RDX(IBAR)*RRN(IBAR)
         IF (CORRECTOR) DS_CORR(IW) = -(UU_AVG-UU(IBAR,JJ,KK))*R(IBAR)*RDX(IBAR)*RRN(IBAR)
         UVW_SAVE(IW) = UU(IBAR,JJ,KK)
         UU(IBAR,JJ,KK) = UU_AVG

      CASE( 2)

         VV_OTHER = 0._EB
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  VV_OTHER = VV_OTHER + OM_VV(IIO,JJO,KKO)*M2%DX(IIO)*M2%DZ(KKO)/DA_OTHER
                  IF (EWC%AREA_RATIO>0.9_EB) OM_VV(IIO,JJO,KKO) = 0.5_EB*(OM_VV(IIO,JJO,KKO)+VV(II,0,KK))
               ENDDO
            ENDDO
         ENDDO
         VV_AVG = 0.5_EB*(VV(II,0,KK) + VV_OTHER)
         IF (PREDICTOR) D_CORR(IW) = DS_CORR(IW) + 0.5*(VV_AVG-VV(II,0,KK))*RDY(1)
         IF (CORRECTOR) DS_CORR(IW) = (VV_AVG-VV(II,0,KK))*RDY(1)
         UVW_SAVE(IW) = VV(II,0,KK)
         VV(II,0,KK)  = VV_AVG

      CASE(-2)

         VV_OTHER = 0._EB
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  VV_OTHER = VV_OTHER + OM_VV(IIO,JJO-1,KKO)*M2%DX(IIO)*M2%DZ(KKO)/DA_OTHER
                  IF (EWC%AREA_RATIO>0.9_EB) OM_VV(IIO,JJO-1,KKO) = 0.5_EB*(OM_VV(IIO,JJO-1,KKO)+VV(II,JBAR,KK))
               ENDDO
            ENDDO
         ENDDO
         VV_AVG = 0.5_EB*(VV(II,JBAR,KK) + VV_OTHER)
         IF (PREDICTOR) D_CORR(IW) = DS_CORR(IW) - 0.5*(VV_AVG-VV(II,JBAR,KK))*RDY(JBAR)
         IF (CORRECTOR) DS_CORR(IW) = -(VV_AVG-VV(II,JBAR,KK))*RDY(JBAR)
         UVW_SAVE(IW)   = VV(II,JBAR,KK)
         VV(II,JBAR,KK) = VV_AVG

      CASE( 3)

         WW_OTHER = 0._EB
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  WW_OTHER = WW_OTHER + OM_WW(IIO,JJO,KKO)*M2%DX(IIO)*M2%DY(JJO)/DA_OTHER
                  IF (EWC%AREA_RATIO>0.9_EB) OM_WW(IIO,JJO,KKO) = 0.5_EB*(OM_WW(IIO,JJO,KKO)+WW(II,JJ,0))
               ENDDO
            ENDDO
         ENDDO
         WW_AVG = 0.5_EB*(WW(II,JJ,0) + WW_OTHER)
         IF (PREDICTOR) D_CORR(IW) = DS_CORR(IW) + 0.5*(WW_AVG-WW(II,JJ,0))*RDZ(1)
         IF (CORRECTOR) DS_CORR(IW) = (WW_AVG-WW(II,JJ,0))*RDZ(1)
         UVW_SAVE(IW) = WW(II,JJ,0)
         WW(II,JJ,0)  = WW_AVG

      CASE(-3)

         WW_OTHER = 0._EB
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  WW_OTHER = WW_OTHER + OM_WW(IIO,JJO,KKO-1)*M2%DX(IIO)*M2%DY(JJO)/DA_OTHER
                  IF (EWC%AREA_RATIO>0.9_EB) OM_WW(IIO,JJO,KKO-1) = 0.5_EB*(OM_WW(IIO,JJO,KKO-1)+WW(II,JJ,KBAR))
               ENDDO
            ENDDO
         ENDDO
         WW_AVG = 0.5_EB*(WW(II,JJ,KBAR) + WW_OTHER)
         IF (PREDICTOR) D_CORR(IW) = DS_CORR(IW) - 0.5*(WW_AVG-WW(II,JJ,KBAR))*RDZ(KBAR)
         IF (CORRECTOR) DS_CORR(IW) = -(WW_AVG-WW(II,JJ,KBAR))*RDZ(KBAR)
         UVW_SAVE(IW)   = WW(II,JJ,KBAR)
         WW(II,JJ,KBAR) = WW_AVG

   END SELECT

   ! Save velocity components at the ghost cell midpoint

   U_GHOST(IW) = 0._EB
   V_GHOST(IW) = 0._EB
   W_GHOST(IW) = 0._EB

   IF (PREDICTOR) OM_UU => OM%US
   IF (CORRECTOR) OM_UU => OM%U
   IF (PREDICTOR) OM_VV => OM%VS
   IF (CORRECTOR) OM_VV => OM%V
   IF (PREDICTOR) OM_WW => OM%WS
   IF (CORRECTOR) OM_WW => OM%W

   DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
      DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
         DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
            U_GHOST(IW) = U_GHOST(IW) + 0.5_EB*(OM_UU(IIO,JJO,KKO)+OM_UU(IIO-1,JJO,KKO))
            V_GHOST(IW) = V_GHOST(IW) + 0.5_EB*(OM_VV(IIO,JJO,KKO)+OM_VV(IIO,JJO-1,KKO))
            W_GHOST(IW) = W_GHOST(IW) + 0.5_EB*(OM_WW(IIO,JJO,KKO)+OM_WW(IIO,JJO,KKO-1))
         ENDDO
      ENDDO
   ENDDO
   NOM_CELLS = REAL((EWC%IIO_MAX-EWC%IIO_MIN+1)*(EWC%JJO_MAX-EWC%JJO_MIN+1)*(EWC%KKO_MAX-EWC%KKO_MIN+1),EB)
   U_GHOST(IW) = U_GHOST(IW)/NOM_CELLS
   V_GHOST(IW) = V_GHOST(IW)/NOM_CELLS
   W_GHOST(IW) = W_GHOST(IW)/NOM_CELLS

ENDDO EXTERNAL_WALL_LOOP

T_USED(4)=T_USED(4)+CURRENT_TIME()-TNOW
END SUBROUTINE MATCH_VELOCITY


SUBROUTINE MATCH_VELOCITY_FLUX(NM)

! Force normal component of velocity flux to match at interpolated boundaries

INTEGER  :: NOM,II,JJ,KK,IOR,IW,IIO,JJO,KKO
INTEGER, INTENT(IN) :: NM
REAL(EB) :: TNOW,DA_OTHER,FVX_OTHER,FVY_OTHER,FVZ_OTHER
TYPE (OMESH_TYPE), POINTER :: OM
TYPE (MESH_TYPE), POINTER :: M2
TYPE (WALL_TYPE), POINTER :: WC
TYPE (EXTERNAL_WALL_TYPE), POINTER :: EWC

IF (NMESHES==1) RETURN
IF (SOLID_PHASE_ONLY) RETURN
IF (EVACUATION_ONLY(NM)) RETURN

TNOW = CURRENT_TIME()

! Assign local variable names

CALL POINT_TO_MESH(NM)

! Loop over all cell edges and determine the appropriate velocity BCs

EXTERNAL_WALL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS

   WC=>WALL(IW)
   EWC=>EXTERNAL_WALL(IW)
   IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY) CYCLE EXTERNAL_WALL_LOOP

   II  = WC%ONE_D%II
   JJ  = WC%ONE_D%JJ
   KK  = WC%ONE_D%KK
   IOR = WC%ONE_D%IOR
   NOM = EWC%NOM
   OM => OMESH(NOM)
   M2 => MESHES(NOM)

   ! Determine the area of the interpolated cell face

   DA_OTHER = 0._EB

   SELECT CASE(ABS(IOR))
      CASE(1)
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  DA_OTHER = DA_OTHER + M2%DY(JJO)*M2%DZ(KKO)
               ENDDO
            ENDDO
         ENDDO
      CASE(2)
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  DA_OTHER = DA_OTHER + M2%DX(IIO)*M2%DZ(KKO)
               ENDDO
            ENDDO
         ENDDO
      CASE(3)
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  DA_OTHER = DA_OTHER + M2%DX(IIO)*M2%DY(JJO)
               ENDDO
            ENDDO
         ENDDO
   END SELECT

   ! Determine the normal component of velocity from the other mesh and use it for average

   SELECT CASE(IOR)

      CASE( 1)

         FVX_OTHER = 0._EB
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  FVX_OTHER = FVX_OTHER + OM%FVX(IIO,JJO,KKO)*M2%DY(JJO)*M2%DZ(KKO)/DA_OTHER
               ENDDO
            ENDDO
         ENDDO
         FVX(0,JJ,KK) = 0.5_EB*(FVX(0,JJ,KK) + FVX_OTHER)

      CASE(-1)

         FVX_OTHER = 0._EB
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  FVX_OTHER = FVX_OTHER + OM%FVX(IIO-1,JJO,KKO)*M2%DY(JJO)*M2%DZ(KKO)/DA_OTHER
               ENDDO
            ENDDO
         ENDDO
         FVX(IBAR,JJ,KK) = 0.5_EB*(FVX(IBAR,JJ,KK) + FVX_OTHER)

      CASE( 2)

         FVY_OTHER = 0._EB
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  FVY_OTHER = FVY_OTHER + OM%FVY(IIO,JJO,KKO)*M2%DX(IIO)*M2%DZ(KKO)/DA_OTHER
               ENDDO
            ENDDO
         ENDDO
         FVY(II,0,KK) = 0.5_EB*(FVY(II,0,KK) + FVY_OTHER)

      CASE(-2)

         FVY_OTHER = 0._EB
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  FVY_OTHER = FVY_OTHER + OM%FVY(IIO,JJO-1,KKO)*M2%DX(IIO)*M2%DZ(KKO)/DA_OTHER
               ENDDO
            ENDDO
         ENDDO
         FVY(II,JBAR,KK) = 0.5_EB*(FVY(II,JBAR,KK) + FVY_OTHER)

      CASE( 3)

         FVZ_OTHER = 0._EB
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  FVZ_OTHER = FVZ_OTHER + OM%FVZ(IIO,JJO,KKO)*M2%DX(IIO)*M2%DY(JJO)/DA_OTHER
               ENDDO
            ENDDO
         ENDDO
         FVZ(II,JJ,0) = 0.5_EB*(FVZ(II,JJ,0) + FVZ_OTHER)

      CASE(-3)

         FVZ_OTHER = 0._EB
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  FVZ_OTHER = FVZ_OTHER + OM%FVZ(IIO,JJO,KKO-1)*M2%DX(IIO)*M2%DY(JJO)/DA_OTHER
               ENDDO
            ENDDO
         ENDDO
         FVZ(II,JJ,KBAR) = 0.5_EB*(FVZ(II,JJ,KBAR) + FVZ_OTHER)

   END SELECT

ENDDO EXTERNAL_WALL_LOOP

T_USED(4)=T_USED(4)+CURRENT_TIME()-TNOW
END SUBROUTINE MATCH_VELOCITY_FLUX


SUBROUTINE CHECK_STABILITY(DT,DT_NEW,NM)

! Checks the Courant and Von Neumann stability criteria, and if necessary, reduces the time step accordingly

INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: DT
REAL(EB) :: UODX,VODY,WODZ,UVW,UVWMAX,R_DX2,MU_MAX,MUTRM,PART_CFL,MU_TMP, UVWMAX_TMP
REAL(EB) :: DT_NEW(NMESHES)
INTEGER  :: I,J,K,IW,IIG,JJG,KKG, ICFL_TMP, JCFL_TMP, KCFL_TMP
TYPE(WALL_TYPE), POINTER :: WC=>NULL()
REAL(EB), PARAMETER :: DT_EPS = 1.E-10_EB

IF (EVACUATION_ONLY(NM)) RETURN

UVWMAX = 0._EB
UVWMAX_TMP = 0._EB
VN     = 0._EB
MUTRM  = 1.E-9_EB
R_DX2  = 1.E-9_EB

! Determine max CFL number from all grid cells
!$OMP PARALLEL PRIVATE(ICFL_TMP, JCFL_TMP, KCFL_TMP, UVWMAX_TMP, UODX, VODY, WODZ, UVW) SHARED(UVWMAX, ICFL, JCFL, KCFL)

!$OMP DO SCHEDULE(STATIC)
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
         UODX = MAXVAL(ABS(US(I-1:I,J,K)))*RDX(I)
         VODY = MAXVAL(ABS(VS(I,J-1:J,K)))*RDY(J)
         WODZ = MAXVAL(ABS(WS(I,J,K-1:K)))*RDZ(K)
         SELECT CASE (CFL_VELOCITY_NORM)
            CASE(0) ; UVW = MAX(UODX,VODY,WODZ) + ABS(DS(I,J,K))
            CASE(1) ; UVW = UODX + VODY + WODZ  + ABS(DS(I,J,K))
            CASE(2) ; UVW = SQRT(UODX**2+VODY**2+WODZ**2) + ABS(DS(I,J,K))
            CASE(3) ; UVW = MAX(UODX,VODY,WODZ)
         END SELECT
         IF (UVW>=UVWMAX_TMP) THEN
            UVWMAX_TMP = UVW
            ICFL_TMP = I
            JCFL_TMP = J
            KCFL_TMP = K
         ENDIF
      ENDDO
   ENDDO
ENDDO
!$OMP END DO
!$OMP CRITICAL
IF(UVWMAX_TMP>UVWMAX) THEN
   UVWMAX = UVWMAX_TMP
   ICFL = ICFL_TMP
   JCFL = JCFL_TMP
   KCFL = KCFL_TMP
ENDIF
!$OMP END CRITICAL
!$OMP END PARALLEL
HEAT_TRANSFER_IF: IF (CHECK_HT) THEN
   WALL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
      WC=>WALL(IW)
      IF (WC%BOUNDARY_TYPE/=SOLID_BOUNDARY) CYCLE WALL_LOOP
      IIG = WC%ONE_D%IIG
      JJG = WC%ONE_D%JJG
      KKG = WC%ONE_D%KKG
      UVW = (ABS(WC%ONE_D%Q_CON_F)/WC%ONE_D%RHO_F)**ONTH * 2._EB*WC%ONE_D%RDN
      IF (UVW>=UVWMAX) THEN
         UVWMAX = UVW
         ICFL=IIG
         JCFL=JJG
         KCFL=KKG
      ENDIF
   ENDDO WALL_LOOP
ENDIF HEAT_TRANSFER_IF

CFL = DT*UVWMAX
PART_CFL = DT*PART_UVWMAX

! Determine max Von Neumann Number for fine grid calcs

PARABOLIC_IF: IF (CHECK_VN) THEN

   MU_MAX = 0._EB
   DO K=1,KBAR
      DO J=1,JBAR
         I_LOOP: DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE I_LOOP
            MU_TMP = MAX(D_Z_MAX(I,J,K),MU(I,J,K)/RHOS(I,J,K))
            IF (MU_TMP>=MU_MAX) THEN
               MU_MAX = MU_TMP
               I_VN=I
               J_VN=J
               K_VN=K
            ENDIF
         ENDDO I_LOOP
      ENDDO
   ENDDO

   IF (TWO_D) THEN
      R_DX2 = RDX(I_VN)**2 + RDZ(K_VN)**2
   ELSE
      R_DX2 = RDX(I_VN)**2 + RDY(J_VN)**2 + RDZ(K_VN)**2
   ENDIF

   MUTRM = MU_MAX
   VN = DT*2._EB*R_DX2*MUTRM

ENDIF PARABOLIC_IF

! Adjust time step size if necessary

IF ((CFL<CFL_MAX .AND. VN<VN_MAX .AND. PART_CFL<PARTICLE_CFL_MAX) .OR. LOCK_TIME_STEP) THEN
   DT_NEW(NM) = DT
   IF (CFL<=CFL_MIN .AND. VN<VN_MIN .AND. PART_CFL<PARTICLE_CFL_MIN .AND. .NOT.LOCK_TIME_STEP) THEN
      SELECT CASE (RESTRICT_TIME_STEP)
         CASE (.TRUE.);  DT_NEW(NM) = MIN(1.1_EB*DT,DT_INITIAL)
         CASE (.FALSE.); DT_NEW(NM) =     1.1_EB*DT
      END SELECT
      CHANGE_TIME_STEP_INDEX(NM) = 1
   ENDIF
ELSE
   DT_NEW(NM) = 0.9_EB*MIN( CFL_MAX/MAX(UVWMAX,DT_EPS)               , &
                            VN_MAX/(2._EB*R_DX2*MAX(MUTRM,DT_EPS))   , &
                            PARTICLE_CFL_MAX/MAX(PART_UVWMAX,DT_EPS))
   CHANGE_TIME_STEP_INDEX(NM) = -1
ENDIF

END SUBROUTINE CHECK_STABILITY


SUBROUTINE BAROCLINIC_CORRECTION(T,NM)

! Add baroclinic term to the momentum equation

USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
REAL(EB), INTENT(IN) :: T
INTEGER, INTENT(IN) :: NM
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU=>NULL(),VV=>NULL(),WW=>NULL(),RHOP=>NULL(),HP=>NULL(),RHMK=>NULL(),RRHO=>NULL()
INTEGER  :: I,J,K,II,JJ,KK,IIG,JJG,KKG,IOR,IW
REAL(EB) :: P_EXTERNAL,TSI,TIME_RAMP_FACTOR,DUMMY,UN,TNOW
LOGICAL  :: INFLOW
TYPE(VENTS_TYPE), POINTER :: VT=>NULL()
TYPE(WALL_TYPE), POINTER :: WC=>NULL()

TNOW=CURRENT_TIME()
CALL POINT_TO_MESH(NM)

! If the baroclinic torque term has been added to the momentum equation RHS, subtract it off.

IF (BAROCLINIC_TERMS_ATTACHED) THEN
   FVX = FVX - FVX_B
   FVY = FVY - FVY_B
   FVZ = FVZ - FVZ_B
ENDIF

BAROCLINIC_TERMS_ATTACHED = .TRUE.

RHMK => WORK1 ! p=rho*(H-K)
RRHO => WORK2 ! reciprocal of rho

IF (PREDICTOR) THEN
   UU => U
   VV => V
   WW => W
   RHOP=>RHO
   IF (PRESSURE_ITERATIONS>1) THEN
      HP => H
   ELSE
      HP => HS
   ENDIF
   ! Note: this ordering of HP=HS in PREDICTOR is required to achieve 2nd order temporal convergence.
   ! We should rethink our notation and re-examine whether both H and HS are required.
ELSE
   UU => US
   VV => VS
   WW => WS
   RHOP=>RHOS
   IF (PRESSURE_ITERATIONS>1) THEN
      HP => HS
   ELSE
      HP => H
   ENDIF
ENDIF

! Compute pressure and 1/rho in each grid cell

!$OMP PARALLEL PRIVATE(WC, VT, TSI, TIME_RAMP_FACTOR, P_EXTERNAL, &
!$OMP& II, JJ, KK, IOR, IIG, JJG, KKG, UN, INFLOW)
!$OMP DO SCHEDULE(static)
DO K=0,KBP1
   DO J=0,JBP1
      DO I=0,IBP1
         RHMK(I,J,K) = RHOP(I,J,K)*(HP(I,J,K)-KRES(I,J,K))
         RRHO(I,J,K) = 1._EB/RHOP(I,J,K)
      ENDDO
   ENDDO
ENDDO
!$OMP END DO

! Set baroclinic term to zero at outflow boundaries and P_EXTERNAL at inflow boundaries

!$OMP MASTER
EXTERNAL_WALL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS
   WC=>WALL(IW)
   IF (WC%BOUNDARY_TYPE/=OPEN_BOUNDARY) CYCLE EXTERNAL_WALL_LOOP
   IF (WC%VENT_INDEX>0) THEN
      VT => VENTS(WC%VENT_INDEX)
      IF (ABS(WC%ONE_D%T_IGN-T_BEGIN)<=SPACING(WC%ONE_D%T_IGN) .AND. VT%PRESSURE_RAMP_INDEX>=1) THEN
         TSI = T
      ELSE
         TSI = T - T_BEGIN
      ENDIF
      TIME_RAMP_FACTOR = EVALUATE_RAMP(TSI,DUMMY,VT%PRESSURE_RAMP_INDEX)
      P_EXTERNAL = TIME_RAMP_FACTOR*VT%DYNAMIC_PRESSURE
   ENDIF
   II  = WC%ONE_D%II
   JJ  = WC%ONE_D%JJ
   KK  = WC%ONE_D%KK
   IOR = WC%ONE_D%IOR
   IIG = WC%ONE_D%IIG
   JJG = WC%ONE_D%JJG
   KKG = WC%ONE_D%KKG
   INFLOW = .FALSE.
   IOR_SELECT: SELECT CASE(IOR)
      CASE( 1); UN = UU(II,JJ,KK)
      CASE(-1); UN = UU(II-1,JJ,KK)
      CASE( 2); UN = VV(II,JJ,KK)
      CASE(-2); UN = VV(II,JJ-1,KK)
      CASE( 3); UN = WW(II,JJ,KK)
      CASE(-3); UN = WW(II,JJ,KK-1)
   END SELECT IOR_SELECT
   IF (UN*SIGN(1._EB,REAL(IOR,EB))>TWO_EPSILON_EB) INFLOW=.TRUE.
   IF (INFLOW) THEN
      RHMK(II,JJ,KK) = 2._EB*P_EXTERNAL - RHMK(IIG,JJG,KKG)  ! Pressure at inflow boundary is P_EXTERNAL
   ELSE
      RHMK(II,JJ,KK) = -RHMK(IIG,JJG,KKG)                    ! No baroclinic correction for outflow boundary
   ENDIF
ENDDO EXTERNAL_WALL_LOOP
!$OMP END MASTER
!$OMP BARRIER

! Compute baroclinic term in the x momentum equation, p*d/dx(1/rho)

!$OMP DO SCHEDULE(static)
DO K=1,KBAR
   DO J=1,JBAR
      DO I=0,IBAR
         FVX_B(I,J,K) = -(RHMK(I,J,K)*RHOP(I+1,J,K)+RHMK(I+1,J,K)*RHOP(I,J,K))*(RRHO(I+1,J,K)-RRHO(I,J,K))*RDXN(I)/ &
                         (RHOP(I+1,J,K)+RHOP(I,J,K))
         FVX(I,J,K) = FVX(I,J,K) + FVX_B(I,J,K)
      ENDDO
   ENDDO
ENDDO
!$OMP END DO nowait

! Compute baroclinic term in the y momentum equation, p*d/dy(1/rho)

IF (.NOT.TWO_D) THEN
!$OMP DO SCHEDULE(static)
   DO K=1,KBAR
      DO J=0,JBAR
         DO I=1,IBAR
            FVY_B(I,J,K) = -(RHMK(I,J,K)*RHOP(I,J+1,K)+RHMK(I,J+1,K)*RHOP(I,J,K))*(RRHO(I,J+1,K)-RRHO(I,J,K))*RDYN(J)/ &
                            (RHOP(I,J+1,K)+RHOP(I,J,K))
            FVY(I,J,K) = FVY(I,J,K) + FVY_B(I,J,K)
         ENDDO
      ENDDO
   ENDDO
!$OMP END DO nowait
ENDIF

! Compute baroclinic term in the z momentum equation, p*d/dz(1/rho)

!$OMP DO SCHEDULE(static)
DO K=0,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         FVZ_B(I,J,K) = -(RHMK(I,J,K)*RHOP(I,J,K+1)+RHMK(I,J,K+1)*RHOP(I,J,K))*(RRHO(I,J,K+1)-RRHO(I,J,K))*RDZN(K)/ &
                         (RHOP(I,J,K+1)+RHOP(I,J,K))
         FVZ(I,J,K) = FVZ(I,J,K) + FVZ_B(I,J,K)
      ENDDO
   ENDDO
ENDDO
!$OMP END DO nowait
!$OMP END PARALLEL

T_USED(4) = T_USED(4) + CURRENT_TIME() - TNOW
END SUBROUTINE BAROCLINIC_CORRECTION


! ----------------------------- PATCH_VELOCITY_FLUX -----------------------------------

SUBROUTINE PATCH_VELOCITY_FLUX(DT,NM)

! The user may specify a polynomial profile using the PROP and DEVC lines. This routine
! specifies the source term in the momentum equation to drive the local velocity toward
! this user-specified value, in much the same way as the immersed boundary method
! (see IBM_VELOCITY_FLUX).

USE DEVICE_VARIABLES, ONLY: DEVICE_TYPE,PROPERTY_TYPE,N_DEVC,DEVICE,PROPERTY
USE TRAN, ONLY: GINV
REAL(EB), INTENT(IN) :: DT
TYPE(DEVICE_TYPE), POINTER :: DV=>NULL()
TYPE(PROPERTY_TYPE), POINTER :: PY=>NULL()
INTEGER, INTENT(IN):: NM
INTEGER :: N,I,J,K,IC1,IC2,I1,I2,J1,J2,K1,K2
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU=>NULL(),VV=>NULL(),WW=>NULL(),HP=>NULL()
REAL(EB) :: VELP,DX0,DY0,DZ0

IF (PREDICTOR) THEN
   UU => U
   VV => V
   WW => W
   HP => H
ELSE
   UU => US
   VV => VS
   WW => WS
   HP => HS
ENDIF

DEVC_LOOP: DO N=1,N_DEVC

   DV=>DEVICE(N)
   IF (DV%QUANTITY/='VELOCITY PATCH') CYCLE DEVC_LOOP
   IF (DV%PROP_INDEX<1)               CYCLE DEVC_LOOP
   IF (.NOT.DEVICE(DV%DEVC_INDEX(1))%CURRENT_STATE) CYCLE DEVC_LOOP

   IF (DV%X1 > XF .OR. DV%X2 < XS .OR. &
       DV%Y1 > YF .OR. DV%Y2 < YS .OR. &
       DV%Z1 > ZF .OR. DV%Z2 < ZS) CYCLE DEVC_LOOP

   PY=>PROPERTY(DV%PROP_INDEX)

   I_VEL_SELECT: SELECT CASE(PY%I_VEL)

      CASE(1) I_VEL_SELECT

         I1 = MAX(0,   NINT( GINV(DV%X1-XS,1,NM)*RDXI   )-1)
         I2 = MIN(IBAR,NINT( GINV(DV%X2-XS,1,NM)*RDXI   )+1)
         J1 = MAX(0,   NINT( GINV(DV%Y1-YS,2,NM)*RDETA  )-1)
         J2 = MIN(JBAR,NINT( GINV(DV%Y2-YS,2,NM)*RDETA  )+1)
         K1 = MAX(0,   NINT( GINV(DV%Z1-ZS,3,NM)*RDZETA )-1)
         K2 = MIN(KBAR,NINT( GINV(DV%Z2-ZS,3,NM)*RDZETA )+1)

         DO K=K1,K2
            DO J=J1,J2
               DO I=I1,I2

                  IC1 = CELL_INDEX(I,J,K)
                  IC2 = CELL_INDEX(I+1,J,K)
                  IF (SOLID(IC1) .OR. SOLID(IC2)) CYCLE

                  IF ( X(I)<DV%X1 .OR.  X(I)>DV%X2) CYCLE ! Inefficient but simple
                  IF (YC(J)<DV%Y1 .OR. YC(J)>DV%Y2) CYCLE
                  IF (ZC(K)<DV%Z1 .OR. ZC(K)>DV%Z2) CYCLE

                  DX0 =  X(I)-DV%X
                  DY0 = YC(J)-DV%Y
                  DZ0 = ZC(K)-DV%Z
                  VELP = PY%P0 + DX0*PY%PX(1) + 0.5_EB*(DX0*DX0*PY%PXX(1,1)+DX0*DY0*PY%PXX(1,2)+DX0*DZ0*PY%PXX(1,3)) &
                               + DY0*PY%PX(2) + 0.5_EB*(DY0*DX0*PY%PXX(2,1)+DY0*DY0*PY%PXX(2,2)+DY0*DZ0*PY%PXX(2,3)) &
                               + DZ0*PY%PX(3) + 0.5_EB*(DZ0*DX0*PY%PXX(3,1)+DZ0*DY0*PY%PXX(3,2)+DZ0*DZ0*PY%PXX(3,3))

                  FVX(I,J,K) = -RDXN(I)*(HP(I+1,J,K)-HP(I,J,K)) - (VELP-UU(I,J,K))/DT
               ENDDO
            ENDDO
         ENDDO

      CASE(2) I_VEL_SELECT

         I1 = MAX(0,   NINT( GINV(DV%X1-XS,1,NM)*RDXI   )-1)
         I2 = MIN(IBAR,NINT( GINV(DV%X2-XS,1,NM)*RDXI   )+1)
         J1 = MAX(0,   NINT( GINV(DV%Y1-YS,2,NM)*RDETA  )-1)
         J2 = MIN(JBAR,NINT( GINV(DV%Y2-YS,2,NM)*RDETA  )+1)
         K1 = MAX(0,   NINT( GINV(DV%Z1-ZS,3,NM)*RDZETA )-1)
         K2 = MIN(KBAR,NINT( GINV(DV%Z2-ZS,3,NM)*RDZETA )+1)

         DO K=K1,K2
            DO J=J1,J2
               DO I=I1,I2

                  IC1 = CELL_INDEX(I,J,K)
                  IC2 = CELL_INDEX(I,J+1,K)

                  IF (SOLID(IC1) .OR. SOLID(IC2)) CYCLE

                  IF (XC(I)<DV%X1 .OR. XC(I)>DV%X2) CYCLE
                  IF ( Y(J)<DV%Y1 .OR.  Y(J)>DV%Y2) CYCLE
                  IF (ZC(K)<DV%Z1 .OR. ZC(K)>DV%Z2) CYCLE

                  DX0 = XC(I)-DV%X
                  DY0 =  Y(J)-DV%Y
                  DZ0 = ZC(K)-DV%Z
                  VELP = PY%P0 + DX0*PY%PX(1) + 0.5_EB*(DX0*DX0*PY%PXX(1,1)+DX0*DY0*PY%PXX(1,2)+DX0*DZ0*PY%PXX(1,3)) &
                               + DY0*PY%PX(2) + 0.5_EB*(DY0*DX0*PY%PXX(2,1)+DY0*DY0*PY%PXX(2,2)+DY0*DZ0*PY%PXX(2,3)) &
                               + DZ0*PY%PX(3) + 0.5_EB*(DZ0*DX0*PY%PXX(3,1)+DZ0*DY0*PY%PXX(3,2)+DZ0*DZ0*PY%PXX(3,3))

                  FVY(I,J,K) = -RDYN(J)*(HP(I,J+1,K)-HP(I,J,K)) - (VELP-VV(I,J,K))/DT
               ENDDO
            ENDDO
         ENDDO

      CASE(3) I_VEL_SELECT

         I1 = MAX(0,   NINT( GINV(DV%X1-XS,1,NM)*RDXI   )-1)
         I2 = MIN(IBAR,NINT( GINV(DV%X2-XS,1,NM)*RDXI   )+1)
         J1 = MAX(0,   NINT( GINV(DV%Y1-YS,2,NM)*RDETA  )-1)
         J2 = MIN(JBAR,NINT( GINV(DV%Y2-YS,2,NM)*RDETA  )+1)
         K1 = MAX(0,   NINT( GINV(DV%Z1-ZS,3,NM)*RDZETA )-1)
         K2 = MIN(KBAR,NINT( GINV(DV%Z2-ZS,3,NM)*RDZETA )+1)

         DO K=K1,K2
            DO J=J1,J2
               DO I=I1,I2

                  IC1 = CELL_INDEX(I,J,K)
                  IC2 = CELL_INDEX(I,J,K+1)
                  IF (SOLID(IC1) .OR. SOLID(IC2)) CYCLE

                  IF (XC(I)<DV%X1 .OR. XC(I)>DV%X2) CYCLE
                  IF (YC(J)<DV%Y1 .OR. YC(J)>DV%Y2) CYCLE
                  IF ( Z(K)<DV%Z1 .OR.  Z(K)>DV%Z2) CYCLE

                  DX0 = XC(I)-DV%X
                  DY0 = YC(J)-DV%Y
                  DZ0 =  Z(K)-DV%Z
                  VELP = PY%P0 + DX0*PY%PX(1) + 0.5_EB*(DX0*DX0*PY%PXX(1,1)+DX0*DY0*PY%PXX(1,2)+DX0*DZ0*PY%PXX(1,3)) &
                               + DY0*PY%PX(2) + 0.5_EB*(DY0*DX0*PY%PXX(2,1)+DY0*DY0*PY%PXX(2,2)+DY0*DZ0*PY%PXX(2,3)) &
                               + DZ0*PY%PX(3) + 0.5_EB*(DZ0*DX0*PY%PXX(3,1)+DZ0*DY0*PY%PXX(3,2)+DZ0*DZ0*PY%PXX(3,3))

                  FVZ(I,J,K) = -RDZN(K)*(HP(I,J,K)-HP(I,J,K+1)) - (VELP-WW(I,J,K))/DT
               ENDDO
            ENDDO
         ENDDO

   END SELECT I_VEL_SELECT

ENDDO DEVC_LOOP

END SUBROUTINE PATCH_VELOCITY_FLUX


! ------------------------ WALL_VELOCITY_NO_GRADH ---------------------------------

SUBROUTINE WALL_VELOCITY_NO_GRADH(DT,STORE_UN)

! This routine recomputes velocities on wall cells, such that the correct
! normal derivative of H is used on the projection. It is only used when the Poisson equation
! for the pressure is solved .NOT. PRES_ON_WHOLE_DOMAIN (i.e. using the GLMAT solver).

REAL(EB), INTENT(IN) :: DT
LOGICAL, INTENT(IN) :: STORE_UN

! Local variables:
INTEGER :: IIG,JJG,KKG,IOR,IW,N_INTERNAL_WALL_CELLS_AUX
REAL(EB) :: DHDN, VEL_N
TYPE (WALL_TYPE), POINTER :: WC
REAL(EB), SAVE, ALLOCATABLE, DIMENSION(:) :: UN_WALLS

N_INTERNAL_WALL_CELLS_AUX=0
IF (.NOT.PRES_ON_WHOLE_DOMAIN) N_INTERNAL_WALL_CELLS_AUX=N_INTERNAL_WALL_CELLS

STORE_UN_COND : IF ( STORE_UN .AND. CORRECTOR) THEN

   ! These velocities from the beginning of step are needed for the velocity fix on wall cells at the corrector
   ! phase (i.e. the loops in VELOCITY_CORRECTOR will change U,V,W to wrong reults using (HP1-HP)/DX gradients,
   ! when the pressure solver in the GLMAT solver.
   IF(ALLOCATED(UN_WALLS)) DEALLOCATE(UN_WALLS)
   ALLOCATE( UN_WALLS(1:N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS_AUX) )
   UN_WALLS(:) = 0._EB

   STORE_LOOP : DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS_AUX

      WC => WALL(IW)
      IIG   = WC%ONE_D%IIG
      JJG   = WC%ONE_D%JJG
      KKG   = WC%ONE_D%KKG
      IOR   = WC%ONE_D%IOR

      SELECT CASE(IOR)
      CASE( IAXIS)
         UN_WALLS(IW) = U(IIG-1,JJG  ,KKG  )
      CASE(-IAXIS)
         UN_WALLS(IW) = U(IIG  ,JJG  ,KKG  )
      CASE( JAXIS)
         UN_WALLS(IW) = V(IIG  ,JJG-1,KKG  )
      CASE(-JAXIS)
         UN_WALLS(IW) = V(IIG  ,JJG  ,KKG  )
      CASE( KAXIS)
         UN_WALLS(IW) = W(IIG  ,JJG  ,KKG-1)
      CASE(-KAXIS)
         UN_WALLS(IW) = W(IIG  ,JJG  ,KKG  )
      END SELECT

   ENDDO STORE_LOOP

   RETURN

ENDIF STORE_UN_COND

! Case of not storing, recompute INTERNAL_WALL_CELL velocities, taking into acct that DHDN=0._EB:
PREDICTOR_COND : IF (PREDICTOR) THEN

  ! Loop internal wall cells -> on OBST surfaces:
  WALL_CELL_LOOP_1: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS_AUX

     WC => WALL(IW)

     IF (WC%BOUNDARY_TYPE/=SOLID_BOUNDARY .AND. WC%BOUNDARY_TYPE/=NULL_BOUNDARY) CYCLE

     IIG   = WC%ONE_D%IIG
     JJG   = WC%ONE_D%JJG
     KKG   = WC%ONE_D%KKG
     IOR   = WC%ONE_D%IOR

     DHDN=0._EB ! Set the normal derivative of H to zero for solids.

     SELECT CASE(IOR)
     CASE( IAXIS)
        US(IIG-1,JJG  ,KKG  ) = (U(IIG-1,JJG  ,KKG  ) - DT*( FVX(IIG-1,JJG  ,KKG  ) + DHDN ))
     CASE(-IAXIS)
        US(IIG  ,JJG  ,KKG  ) = (U(IIG  ,JJG  ,KKG  ) - DT*( FVX(IIG  ,JJG  ,KKG  ) + DHDN ))
     CASE( JAXIS)
        VS(IIG  ,JJG-1,KKG  ) = (V(IIG  ,JJG-1,KKG  ) - DT*( FVY(IIG  ,JJG-1,KKG  ) + DHDN ))
     CASE(-JAXIS)
        VS(IIG  ,JJG  ,KKG  ) = (V(IIG  ,JJG  ,KKG  ) - DT*( FVY(IIG  ,JJG  ,KKG  ) + DHDN ))
     CASE( KAXIS)
        WS(IIG  ,JJG  ,KKG-1) = (W(IIG  ,JJG  ,KKG-1) - DT*( FVZ(IIG  ,JJG  ,KKG-1) + DHDN ))
     CASE(-KAXIS)
        WS(IIG  ,JJG  ,KKG  ) = (W(IIG  ,JJG  ,KKG  ) - DT*( FVZ(IIG  ,JJG  ,KKG  ) + DHDN ))
     END SELECT

  ENDDO WALL_CELL_LOOP_1

ELSE ! Corrector

  ! Loop internal wall cells -> on OBST surfaces:
  WALL_CELL_LOOP_2: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS_AUX

     WC => WALL(IW)

     IF (WC%BOUNDARY_TYPE/=SOLID_BOUNDARY .AND. WC%BOUNDARY_TYPE/=NULL_BOUNDARY) CYCLE

     IIG   = WC%ONE_D%IIG
     JJG   = WC%ONE_D%JJG
     KKG   = WC%ONE_D%KKG
     IOR   = WC%ONE_D%IOR

     DHDN=0._EB ! Set the normal derivative of H to zero for solids.

     VEL_N = UN_WALLS(IW)

     SELECT CASE(IOR)
     CASE( IAXIS)                                 ! | - Problem with this is it was modified in VELOCITY_CORRECTOR,
                                                  ! V   => Store the untouched U normal on internal WALLs.
         U(IIG-1,JJG  ,KKG  ) = 0.5_EB*(                      VEL_N + US(IIG-1,JJG  ,KKG  ) - &
                                        DT*( FVX(IIG-1,JJG  ,KKG  ) + DHDN ))
     CASE(-IAXIS)
         U(IIG  ,JJG  ,KKG  ) = 0.5_EB*(                      VEL_N + US(IIG  ,JJG  ,KKG  ) - &
                                        DT*( FVX(IIG  ,JJG  ,KKG  ) + DHDN ))
     CASE( JAXIS)
         V(IIG  ,JJG-1,KKG  ) = 0.5_EB*(                      VEL_N + VS(IIG  ,JJG-1,KKG  ) - &
                                        DT*( FVY(IIG  ,JJG-1,KKG  ) + DHDN ))
     CASE(-JAXIS)
         V(IIG  ,JJG  ,KKG  ) = 0.5_EB*(                      VEL_N + VS(IIG  ,JJG  ,KKG  ) - &
                                        DT*( FVY(IIG  ,JJG  ,KKG  ) + DHDN ))
     CASE( KAXIS)
         W(IIG  ,JJG  ,KKG-1) = 0.5_EB*(                      VEL_N + WS(IIG  ,JJG  ,KKG-1) - &
                                        DT*( FVZ(IIG  ,JJG  ,KKG-1) + DHDN ))
     CASE(-KAXIS)
         W(IIG  ,JJG  ,KKG  ) = 0.5_EB*(                      VEL_N + WS(IIG  ,JJG  ,KKG  ) - &
                                        DT*( FVZ(IIG  ,JJG  ,KKG  ) + DHDN ))
     END SELECT

  ENDDO WALL_CELL_LOOP_2

  DEALLOCATE(UN_WALLS)

ENDIF PREDICTOR_COND

RETURN
END SUBROUTINE WALL_VELOCITY_NO_GRADH

END MODULE VELO
