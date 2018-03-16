MODULE WALL_ROUTINES

! Compute the wall boundary conditions

USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_POINTERS

IMPLICIT NONE
PRIVATE

PUBLIC WALL_BC,TGA_ANALYSIS

TYPE(WALL_TYPE), POINTER :: WC
TYPE(CFACE_TYPE), POINTER :: CFA
TYPE(EXTERNAL_WALL_TYPE), POINTER :: EWC
TYPE(LAGRANGIAN_PARTICLE_TYPE), POINTER :: LP
TYPE(LAGRANGIAN_PARTICLE_CLASS_TYPE), POINTER :: LPC
TYPE(ONE_D_M_AND_E_XFER_TYPE), POINTER :: ONE_D
TYPE(SURFACE_TYPE), POINTER :: SF
TYPE(VENTS_TYPE), POINTER :: VT
TYPE(OMESH_TYPE), POINTER :: OM
TYPE(MESH_TYPE), POINTER :: MM
TYPE(MATERIAL_TYPE), POINTER :: ML
LOGICAL :: CALL_HT_1D

CONTAINS


SUBROUTINE WALL_BC(T,DT,NM)

! This is the main control routine for this module

USE COMP_FUNCTIONS, ONLY: CURRENT_TIME
USE SOOT_ROUTINES, ONLY: CALC_DEPOSITION,SURFACE_OXIDATION
REAL(EB) :: TNOW
REAL(EB), INTENT(IN) :: T,DT
INTEGER, INTENT(IN) :: NM

IF (EVACUATION_ONLY(NM)) RETURN

TNOW=CURRENT_TIME()

CALL POINT_TO_MESH(NM)

! Compute the temperature TMP_F at all boundary cells, including PYROLYSIS and 1-D heat transfer

CALL THERMAL_BC(T,NM)

! Compute rho*D at WALL cells

CALL DIFFUSIVITY_BC

! Special boundary routines

IF (DEPOSITION .AND. .NOT.INITIALIZATION_PHASE) CALL CALC_DEPOSITION(DT,NM)
IF (SOOT_OXIDATION) CALL SURFACE_OXIDATION(DT,NM)
IF (HVAC_SOLVE .AND. .NOT.INITIALIZATION_PHASE) CALL HVAC_BC

! Compute the species mass fractions, ZZ_F, at all boundary cells

CALL SPECIES_BC(T,DT,NM)

! Compute the density, RHO_F, at WALL cells only

CALL DENSITY_BC

! Special transport equation BCs

IF (TRANSPORT_UNMIXED_FRACTION) CALL ZETA_BC

T_USED(6)=T_USED(6)+CURRENT_TIME()-TNOW
END SUBROUTINE WALL_BC


SUBROUTINE THERMAL_BC(T,NM)

! Thermal boundary conditions for all boundaries.
! One dimensional heat transfer and pyrolysis is done in PYROLYSIS.
! Note also that gas phase values are assigned here to be used for all subsequent BCs.

USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
USE PHYSICAL_FUNCTIONS, ONLY : GET_SPECIFIC_GAS_CONSTANT,GET_SPECIFIC_HEAT,GET_VISCOSITY,&
                               GET_SOLID_CONDUCTIVITY,GET_SOLID_RHOCBAR,GET_SOLID_ABSORPTION_COEFFICIENT
USE COMPLEX_GEOMETRY, ONLY : CFACE_THERMAL_GASVARS
REAL(EB), INTENT(IN) :: T
REAL(EB) :: DT_BC,DTMP,DT_BC_HT3D
INTEGER  :: SURF_INDEX,IW,IP,ICF
INTEGER, INTENT(IN) :: NM
REAL(EB), POINTER, DIMENSION(:,:) :: PBAR_P
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU=>NULL(),VV=>NULL(),WW=>NULL(),RHOP=>NULL(),OM_RHOP=>NULL(),OM_TMP=>NULL()
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZP=>NULL()

IF (VEG_LEVEL_SET_UNCOUPLED) RETURN

IF (PREDICTOR) THEN
   UU => US
   VV => VS
   WW => WS
   RHOP => RHOS
   ZZP  => ZZS
   PBAR_P => PBAR_S
ELSE
   UU => U
   VV => V
   WW => W
   RHOP => RHO
   ZZP  => ZZ
   PBAR_P => PBAR
ENDIF

! For thermally-thick boundary conditions, set the flag to call the routine PYROLYSIS

CALL_HT_1D = .FALSE.
IF (.NOT.INITIALIZATION_PHASE .AND. CORRECTOR) THEN
   WALL_COUNTER = WALL_COUNTER + 1
   IF (WALL_COUNTER==WALL_INCREMENT) THEN
      DT_BC    = T - BC_CLOCK
      BC_CLOCK = T
      CALL_HT_1D = .TRUE.
      WALL_COUNTER = 0
   ENDIF
ENDIF

! Loop through all wall cells and apply heat transfer BCs

WALL_CELL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
   WC=>WALL(IW)
   IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY .AND. .NOT.SOLID_HT3D) CYCLE WALL_CELL_LOOP
   SURF_INDEX = WC%SURF_INDEX
   ONE_D => WC%ONE_D
   ONE_D%TMP_G   =  TMP(ONE_D%IIG,ONE_D%JJG,ONE_D%KKG)
   ONE_D%RHO_G   = RHOP(ONE_D%IIG,ONE_D%JJG,ONE_D%KKG)
   ONE_D%ZZ_G(1:N_TRACKED_SPECIES) = ZZP(ONE_D%IIG,ONE_D%JJG,ONE_D%KKG,1:N_TRACKED_SPECIES)
   ONE_D%RSUM_G  = RSUM(ONE_D%IIG,ONE_D%JJG,ONE_D%KKG)
   ONE_D%MU_G    = MU(ONE_D%IIG,ONE_D%JJG,ONE_D%KKG)
   ONE_D%U_TANG  = SQRT(2._EB*KRES(ONE_D%IIG,ONE_D%JJG,ONE_D%KKG))
   CALL CALCULATE_TMP_F(WALL_INDEX=IW)
   IF (SURFACE(SURF_INDEX)%THERMAL_BC_INDEX==THERMALLY_THICK .AND. CALL_HT_1D) &
      CALL SOLID_HEAT_TRANSFER_1D(NM,T,DT_BC,WALL_INDEX=IW)
ENDDO WALL_CELL_LOOP

! Loop through all CFACEs and apply heat transfer BCs

CFACE_LOOP: DO ICF=1,N_CFACE_CELLS
   CFA=>CFACE(ICF)
   IF (CFA%BOUNDARY_TYPE==NULL_BOUNDARY .AND. .NOT.SOLID_HT3D) CYCLE CFACE_LOOP
   SURF_INDEX = CFA%SURF_INDEX
   ONE_D => CFA%ONE_D
   ! Populate ONE_D%TMP_G, ONE_D%RHO_G, ONE_D%ZZ_G(:), ONE_D%RSUM_G, ONE_D%U_TANG
   CALL CFACE_THERMAL_GASVARS(ICF,ONE_D)
   CALL CALCULATE_TMP_F(CFACE_INDEX=ICF)
   IF (SURFACE(SURF_INDEX)%THERMAL_BC_INDEX==THERMALLY_THICK .AND. CALL_HT_1D) &
      CALL SOLID_HEAT_TRANSFER_1D(NM,T,DT_BC,CFACE_INDEX=ICF)
ENDDO CFACE_LOOP

! Loop through all particles and apply heat transfer BCs

IF (SOLID_PARTICLES) THEN
   DO IP = 1, NLP
      LP => LAGRANGIAN_PARTICLE(IP)
      LPC => LAGRANGIAN_PARTICLE_CLASS(LP%CLASS_INDEX)
      IF (LPC%SOLID_PARTICLE .OR. LPC%MASSLESS_TARGET) THEN  ! Target particles are included to get gas phase values
         SURF_INDEX = LPC%SURF_INDEX
         ONE_D => LP%ONE_D
         ONE_D%TMP_G   =  TMP(ONE_D%IIG,ONE_D%JJG,ONE_D%KKG)
         ONE_D%RHO_G   = RHOP(ONE_D%IIG,ONE_D%JJG,ONE_D%KKG)
         ONE_D%ZZ_G(1:N_TRACKED_SPECIES) = ZZP(ONE_D%IIG,ONE_D%JJG,ONE_D%KKG,1:N_TRACKED_SPECIES)
         ONE_D%RSUM_G  = RSUM(ONE_D%IIG,ONE_D%JJG,ONE_D%KKG)
         ONE_D%MU_G    = MU(ONE_D%IIG,ONE_D%JJG,ONE_D%KKG)
         ONE_D%U_TANG  = SQRT(2._EB*KRES(ONE_D%IIG,ONE_D%JJG,ONE_D%KKG))
         IF (LPC%SOLID_PARTICLE) CALL CALCULATE_TMP_F(PARTICLE_INDEX=IP)
         IF (SURFACE(SURF_INDEX)%THERMAL_BC_INDEX==THERMALLY_THICK .AND. CALL_HT_1D) &
            CALL SOLID_HEAT_TRANSFER_1D(NM,T,DT_BC,PARTICLE_INDEX=IP)
      ENDIF
   ENDDO
ENDIF

! *********************** UNDER CONSTRUCTION **************************
IF (.NOT.INITIALIZATION_PHASE .AND. SOLID_HT3D .AND. CORRECTOR) THEN
   WALL_COUNTER_HT3D = WALL_COUNTER_HT3D + 1
   IF (WALL_COUNTER_HT3D==WALL_INCREMENT_HT3D) THEN
      DT_BC_HT3D    = T - BC_CLOCK_HT3D
      BC_CLOCK_HT3D = T
      CALL SOLID_HEAT_TRANSFER_3D
      WALL_COUNTER_HT3D = 0
   ENDIF
ENDIF
! *********************************************************************

CONTAINS


SUBROUTINE CALCULATE_TMP_F(WALL_INDEX,CFACE_INDEX,PARTICLE_INDEX)

USE MASS, ONLY: SCALAR_FACE_VALUE

INTEGER, INTENT(IN), OPTIONAL :: WALL_INDEX,CFACE_INDEX,PARTICLE_INDEX
REAL(EB) :: ARO,FDERIV,QEXTRA,QNET,RAMP_FACTOR,RHO_G_2,RSUM_F,PBAR_F,TMP_OTHER,TSI,UN, &
            RHO_ZZ_F(1:N_TOTAL_SCALARS),ZZ_GET(1:N_TRACKED_SPECIES),DUMMY, &
            ZZZ(1:4),RHO_OTHER,RHO_OTHER_2,RHO_ZZ_OTHER(1:N_TOTAL_SCALARS),RHO_ZZ_OTHER_2,RHO_ZZ_G,RHO_ZZ_G_2, &
            DDO,PBAR_G,PBAR_OTHER,DENOM

LOGICAL :: INFLOW,SECOND_ORDER_INTERPOLATED_BOUNDARY,SOLID_OTHER,ATMOSPHERIC_INTERPOLATION
INTEGER :: II,JJ,KK,IIG,JJG,KKG,IOR,IIO,JJO,KKO,N,ADCOUNT,ICG,ICO
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: OM_ZZP=>NULL()

SF  => SURFACE(SURF_INDEX)

IF (PRESENT(WALL_INDEX)) THEN
   WC=>WALL(WALL_INDEX)
   ONE_D => WC%ONE_D
   II  = ONE_D%II
   JJ  = ONE_D%JJ
   KK  = ONE_D%KK
ELSEIF (PRESENT(CFACE_INDEX)) THEN
   CFA=>CFACE(CFACE_INDEX)
   ONE_D => CFA%ONE_D
   KK  = CUT_FACE(CFA%CUT_FACE_IND1)%IJK(KAXIS) ! CUT_FACE type INBOUNDARY -> KK is under-laying Cartesian cell index.
ELSEIF (PRESENT(PARTICLE_INDEX)) THEN
   LP=>LAGRANGIAN_PARTICLE(PARTICLE_INDEX)
   ONE_D => LP%ONE_D
ENDIF

IIG = ONE_D%IIG
JJG = ONE_D%JJG
KKG = ONE_D%KKG
IOR = ONE_D%IOR

! Compute surface temperature, TMP_F, and convective heat flux, QCONF, for various boundary conditions

METHOD_OF_HEAT_TRANSFER: SELECT CASE(SF%THERMAL_BC_INDEX)

   CASE (NO_CONVECTION) METHOD_OF_HEAT_TRANSFER

      ONE_D%TMP_F  = ONE_D%TMP_G

   CASE (INFLOW_OUTFLOW) METHOD_OF_HEAT_TRANSFER  ! Only for WALL cells

      ! Base inflow/outflow decision on velocity component with same predictor/corrector attribute

      INFLOW = .FALSE.
      SELECT CASE(IOR)
         CASE( 1)
            UN = UU(II,JJ,KK)
         CASE(-1)
            UN = -UU(II-1,JJ,KK)
         CASE( 2)
            UN = VV(II,JJ,KK)
         CASE(-2)
            UN = -VV(II,JJ-1,KK)
         CASE( 3)
            UN = WW(II,JJ,KK)
         CASE(-3)
            UN = -WW(II,JJ,KK-1)
      END SELECT
      IF (UN>TWO_EPSILON_EB) INFLOW = .TRUE.

      IF (INFLOW) THEN
         ONE_D%TMP_F = TMP_0(KK)
         IF (WC%VENT_INDEX>0) THEN
            VT => VENTS(WC%VENT_INDEX)
            IF (VT%TMP_EXTERIOR>0._EB) THEN
               TSI = T - T_BEGIN
               ONE_D%TMP_F = TMP_0(KK) + EVALUATE_RAMP(TSI,DUMMY,VT%TMP_EXTERIOR_RAMP_INDEX)*(VT%TMP_EXTERIOR-TMP_0(KK))
            ENDIF
         ENDIF
         ONE_D%ZZ_F(1:N_TRACKED_SPECIES)=SPECIES_MIXTURE(1:N_TRACKED_SPECIES)%ZZ0
      ELSE
         ONE_D%TMP_F = ONE_D%TMP_G
         ONE_D%ZZ_F(1:N_TRACKED_SPECIES)=ONE_D%ZZ_G(1:N_TRACKED_SPECIES)
      ENDIF

      ! Ghost cell values

      TMP(II,JJ,KK) = ONE_D%TMP_F
      ZZP(II,JJ,KK,1:N_TRACKED_SPECIES) = ONE_D%ZZ_F(1:N_TRACKED_SPECIES)
      ZZ_GET(1:N_TRACKED_SPECIES) = MAX(0._EB,ZZP(II,JJ,KK,1:N_TRACKED_SPECIES))
      CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,RSUM(II,JJ,KK))
      RHOP(II,JJ,KK) = PBAR_P(KK,ONE_D%PRESSURE_ZONE)/(RSUM(II,JJ,KK)*TMP(II,JJ,KK))

      ONE_D%QCONF = 2._EB*WC%ONE_D%K_G*(ONE_D%TMP_G-ONE_D%TMP_F)*ONE_D%RDN

   CASE (SPECIFIED_TEMPERATURE) METHOD_OF_HEAT_TRANSFER

      IF (ABS(ONE_D%T_IGN-T_BEGIN) <= SPACING(ONE_D%T_IGN) .AND. SF%RAMP_INDEX(TIME_TEMP)>=1) THEN
         TSI = T
      ELSE
         TSI = T - ONE_D%T_IGN
      ENDIF

      IF (ONE_D%UW<=0._EB) THEN
         IF (SF%TMP_FRONT>0._EB) THEN
            ONE_D%TMP_F = TMP_0(KKG) + EVALUATE_RAMP(TSI,SF%TAU(TIME_TEMP),SF%RAMP_INDEX(TIME_TEMP))*(SF%TMP_FRONT-TMP_0(KKG))
         ELSE
            ONE_D%TMP_F = TMP_0(KKG)
         ENDIF
      ELSE
         ONE_D%TMP_F = ONE_D%TMP_G ! If gas is being drawn from the domain, set the boundary temperature to the gas temperature
      ENDIF

      DTMP = ONE_D%TMP_G - ONE_D%TMP_F
      IF (PRESENT(WALL_INDEX)) THEN
         ONE_D%HEAT_TRANS_COEF = HEAT_TRANSFER_COEFFICIENT(DTMP,SF%H_FIXED,SURF_INDEX,WALL_INDEX=WALL_INDEX)
      ELSEIF (PRESENT(CFACE_INDEX)) THEN
         ONE_D%HEAT_TRANS_COEF = HEAT_TRANSFER_COEFFICIENT(DTMP,SF%H_FIXED,SURF_INDEX,CFACE_INDEX=CFACE_INDEX)
      ELSEIF (PRESENT(PARTICLE_INDEX)) THEN
         ONE_D%HEAT_TRANS_COEF = HEAT_TRANSFER_COEFFICIENT(DTMP,SF%H_FIXED,SURF_INDEX,PARTICLE_INDEX=PARTICLE_INDEX)
      ENDIF
      ONE_D%QCONF = ONE_D%HEAT_TRANS_COEF*DTMP

   CASE (NET_FLUX_BC) METHOD_OF_HEAT_TRANSFER

      IF (ABS(ONE_D%T_IGN-T_BEGIN)<= SPACING(ONE_D%T_IGN) .AND. SF%RAMP_INDEX(TIME_HEAT)>=1) THEN
         TSI = T
      ELSE
         TSI = T - ONE_D%T_IGN
      ENDIF
      TMP_OTHER = ONE_D%TMP_F
      RAMP_FACTOR = EVALUATE_RAMP(TSI,SF%TAU(TIME_HEAT),SF%RAMP_INDEX(TIME_HEAT))
      QNET = -RAMP_FACTOR*SF%NET_HEAT_FLUX*ONE_D%AREA_ADJUST
      ADCOUNT = 0
      ADLOOP: DO
         ADCOUNT = ADCOUNT + 1
         DTMP = ONE_D%TMP_G - TMP_OTHER
         IF (ABS(QNET) > 0._EB .AND. ABS(DTMP) <TWO_EPSILON_EB) DTMP=1._EB
         IF (PRESENT(WALL_INDEX)) THEN
            ONE_D%HEAT_TRANS_COEF = HEAT_TRANSFER_COEFFICIENT(DTMP,SF%H_FIXED,SURF_INDEX,WALL_INDEX=WALL_INDEX)
         ELSEIF (PRESENT(CFACE_INDEX)) THEN
            ONE_D%HEAT_TRANS_COEF = HEAT_TRANSFER_COEFFICIENT(DTMP,SF%H_FIXED,SURF_INDEX,CFACE_INDEX=CFACE_INDEX)
         ELSEIF (PRESENT(PARTICLE_INDEX)) THEN
            ONE_D%HEAT_TRANS_COEF = HEAT_TRANSFER_COEFFICIENT(DTMP,SF%H_FIXED,SURF_INDEX,PARTICLE_INDEX=PARTICLE_INDEX)
         ENDIF
         IF (RADIATION) THEN
            QEXTRA = ONE_D%HEAT_TRANS_COEF*DTMP + ONE_D%QRADIN - ONE_D%EMISSIVITY * SIGMA * TMP_OTHER ** 4 - QNET
            FDERIV = -ONE_D%HEAT_TRANS_COEF -  4._EB * ONE_D%EMISSIVITY * SIGMA * TMP_OTHER ** 3
         ELSE
            QEXTRA = ONE_D%HEAT_TRANS_COEF*DTMP - QNET
            FDERIV = -ONE_D%HEAT_TRANS_COEF
         ENDIF
         IF (ABS(FDERIV) > TWO_EPSILON_EB) TMP_OTHER = TMP_OTHER - QEXTRA / FDERIV
         IF (ABS(TMP_OTHER - ONE_D%TMP_F) / ONE_D%TMP_F < 1.E-4_EB .OR. ADCOUNT > 20) THEN
            ONE_D%TMP_F = MIN(TMPMAX,TMP_OTHER)
            EXIT ADLOOP
         ELSE
            ONE_D%TMP_F = MIN(TMPMAX,TMP_OTHER)
            CYCLE ADLOOP
         ENDIF
      ENDDO ADLOOP

      ONE_D%QCONF = ONE_D%HEAT_TRANS_COEF*DTMP

   CASE (CONVECTIVE_FLUX_BC) METHOD_OF_HEAT_TRANSFER

      IF (ABS(ONE_D%T_IGN-T_BEGIN) <= SPACING(ONE_D%T_IGN) .AND. SF%RAMP_INDEX(TIME_HEAT)>=1) THEN
         TSI = T
      ELSE
         TSI = T - ONE_D%T_IGN
      ENDIF
      RAMP_FACTOR = EVALUATE_RAMP(TSI,SF%TAU(TIME_HEAT),SF%RAMP_INDEX(TIME_HEAT))
      IF (SF%TMP_FRONT>0._EB) THEN
         ONE_D%TMP_F =  TMPA + RAMP_FACTOR*(SF%TMP_FRONT-TMPA)
      ELSE
         ONE_D%TMP_F =  TMP_0(KK)
      ENDIF
      ONE_D%QCONF = -RAMP_FACTOR*SF%CONVECTIVE_HEAT_FLUX*ONE_D%AREA_ADJUST

   CASE (INTERPOLATED_BC) METHOD_OF_HEAT_TRANSFER  ! Only for EXTERNAL_WALL_CELLs

      EWC => EXTERNAL_WALL(WALL_INDEX)
      OM => OMESH(EWC%NOM)
      IF (PREDICTOR) THEN
         OM_RHOP => OM%RHOS
         OM_ZZP => OM%ZZS
      ELSE
         OM_RHOP => OM%RHO
         OM_ZZP => OM%ZZ
      ENDIF
      IF (SOLID_HT3D) OM_TMP => OM%TMP
      MM => MESHES(EWC%NOM)

      ! Gather data from other mesh

      RHO_OTHER=0._EB
      RHO_ZZ_OTHER=0._EB
      TMP_OTHER=0._EB
      SOLID_OTHER=.FALSE.
      DDO=1._EB

      DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
         DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
            DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
               SELECT CASE(IOR)
                  CASE( 1)
                     ARO = MIN(1._EB , (MM%DY(JJO)*MM%DZ(KKO))/(DY(JJ)*DZ(KK)) )
                  CASE(-1)
                     ARO = MIN(1._EB , (MM%DY(JJO)*MM%DZ(KKO))/(DY(JJ)*DZ(KK)) )
                  CASE( 2)
                     ARO = MIN(1._EB , (MM%DX(IIO)*MM%DZ(KKO))/(DX(II)*DZ(KK)) )
                  CASE(-2)
                     ARO = MIN(1._EB , (MM%DX(IIO)*MM%DZ(KKO))/(DX(II)*DZ(KK)) )
                  CASE( 3)
                     ARO = MIN(1._EB , (MM%DX(IIO)*MM%DY(JJO))/(DX(II)*DY(JJ)) )
                     DDO = (DZ(KK)+DZ(KKG))/(MM%DZ(KKO)+DZ(KKG))
                  CASE(-3)
                     ARO = MIN(1._EB , (MM%DX(IIO)*MM%DY(JJO))/(DX(II)*DY(JJ)) )
                     DDO = (DZ(KK)+DZ(KKG))/(MM%DZ(KKO)+DZ(KKG))
               END SELECT
               RHO_OTHER = RHO_OTHER + ARO*OM_RHOP(IIO,JJO,KKO)      ! average multiple face values
               RHO_ZZ_OTHER(1:N_TOTAL_SCALARS) = RHO_ZZ_OTHER(1:N_TOTAL_SCALARS) &
                  + ARO*OM_RHOP(IIO,JJO,KKO)*OM_ZZP(IIO,JJO,KKO,1:N_TOTAL_SCALARS)
               IF (SOLID_HT3D) THEN
                  TMP_OTHER = TMP_OTHER + ARO*OM_TMP(IIO,JJO,KKO)
                  ICO = MM%CELL_INDEX(IIO,JJO,KKO)
                  IF (MM%SOLID(ICO)) SOLID_OTHER=.TRUE.
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      ! Determine if there are 4 equally sized cells spanning the interpolated boundary

      SECOND_ORDER_INTERPOLATED_BOUNDARY = .FALSE.
      IF (ABS(EWC%AREA_RATIO-1._EB)<0.01_EB) THEN
         IIO = EWC%IIO_MIN
         JJO = EWC%JJO_MIN
         KKO = EWC%KKO_MIN
         SELECT CASE(IOR)
            CASE( 1) ; ICG = CELL_INDEX(IIG+1,JJG,KKG) ; ICO = MM%CELL_INDEX(IIO-1,JJO,KKO)
            CASE(-1) ; ICG = CELL_INDEX(IIG-1,JJG,KKG) ; ICO = MM%CELL_INDEX(IIO+1,JJO,KKO)
            CASE( 2) ; ICG = CELL_INDEX(IIG,JJG+1,KKG) ; ICO = MM%CELL_INDEX(IIO,JJO-1,KKO)
            CASE(-2) ; ICG = CELL_INDEX(IIG,JJG-1,KKG) ; ICO = MM%CELL_INDEX(IIO,JJO+1,KKO)
            CASE( 3) ; ICG = CELL_INDEX(IIG,JJG,KKG+1) ; ICO = MM%CELL_INDEX(IIO,JJO,KKO-1)
            CASE(-3) ; ICG = CELL_INDEX(IIG,JJG,KKG-1) ; ICO = MM%CELL_INDEX(IIO,JJO,KKO+1)
         END SELECT
         IF (.NOT.SOLID(ICG) .AND. .NOT.MM%SOLID(ICO)) SECOND_ORDER_INTERPOLATED_BOUNDARY = .TRUE.
      ENDIF

      ! Density

      ATMOSPHERIC_INTERPOLATION = .FALSE.
      IF (USE_ATMOSPHERIC_INTERPOLATION .AND. STRATIFICATION .AND. ABS(DDO-1._EB)>0.01_EB .AND. ABS(IOR)==3) &
         ATMOSPHERIC_INTERPOLATION = .TRUE.

      IF (ATMOSPHERIC_INTERPOLATION) THEN
         ! interp or extrap RHO_OTHER for jump in vertical grid resolution, linear in temperature to match heat flux in divg
         PBAR_G = PBAR_P(KKG,ONE_D%PRESSURE_ZONE)
         PBAR_OTHER = MM%PBAR(EWC%KKO_MIN,ONE_D%PRESSURE_ZONE)
         DENOM = PBAR_G/ONE_D%RHO_G + DDO*(PBAR_OTHER/RHO_OTHER - PBAR_G/ONE_D%RHO_G)
         RHOP(II,JJ,KK) = PBAR_P(KK,ONE_D%PRESSURE_ZONE)/DENOM
      ELSE
         RHOP(II,JJ,KK) = RHO_OTHER
      ENDIF

      RHO_G_2        = ONE_D%RHO_G ! first-order
      RHO_OTHER      = RHOP(II,JJ,KK)
      RHO_OTHER_2    = RHOP(II,JJ,KK)

      SELECT CASE(IOR)
         CASE( 1)
            IF (SECOND_ORDER_INTERPOLATED_BOUNDARY) THEN
               RHO_G_2 = RHOP(IIG+1,JJG,KKG)
               RHO_OTHER_2 = OM_RHOP(IIO-1,JJO,KKO)
            ENDIF
            ZZZ(1:4) = (/RHO_OTHER_2,RHO_OTHER,ONE_D%RHO_G,RHO_G_2/)
            ONE_D%RHO_F = SCALAR_FACE_VALUE(UU(II,JJ,KK),ZZZ,FLUX_LIMITER)
         CASE(-1)
            IF (SECOND_ORDER_INTERPOLATED_BOUNDARY) THEN
               RHO_G_2 = RHOP(IIG-1,JJG,KKG)
               RHO_OTHER_2 = OM_RHOP(IIO+1,JJO,KKO)
            ENDIF
            ZZZ(1:4) = (/RHO_G_2,ONE_D%RHO_G,RHO_OTHER,RHO_OTHER_2/)
            ONE_D%RHO_F = SCALAR_FACE_VALUE(UU(II-1,JJ,KK),ZZZ,FLUX_LIMITER)
         CASE( 2)
            IF (SECOND_ORDER_INTERPOLATED_BOUNDARY) THEN
               RHO_G_2 = RHOP(IIG,JJG+1,KKG)
               RHO_OTHER_2 = OM_RHOP(IIO,JJO-1,KKO)
            ENDIF
            ZZZ(1:4) = (/RHO_OTHER_2,RHO_OTHER,ONE_D%RHO_G,RHO_G_2/)
            ONE_D%RHO_F = SCALAR_FACE_VALUE(VV(II,JJ,KK),ZZZ,FLUX_LIMITER)
         CASE(-2)
            IF (SECOND_ORDER_INTERPOLATED_BOUNDARY) THEN
               RHO_G_2 = RHOP(IIG,JJG-1,KKG)
               RHO_OTHER_2 = OM_RHOP(IIO,JJO+1,KKO)
            ENDIF
            ZZZ(1:4) = (/RHO_G_2,ONE_D%RHO_G,RHO_OTHER,RHO_OTHER_2/)
            ONE_D%RHO_F = SCALAR_FACE_VALUE(VV(II,JJ-1,KK),ZZZ,FLUX_LIMITER)
         CASE( 3)
            IF (SECOND_ORDER_INTERPOLATED_BOUNDARY) THEN
               RHO_G_2 = RHOP(IIG,JJG,KKG+1)
               RHO_OTHER_2 = OM_RHOP(IIO,JJO,KKO-1)
            ENDIF
            ZZZ(1:4) = (/RHO_OTHER_2,RHO_OTHER,ONE_D%RHO_G,RHO_G_2/)
            ONE_D%RHO_F = SCALAR_FACE_VALUE(WW(II,JJ,KK),ZZZ,FLUX_LIMITER)
         CASE(-3)
            IF (SECOND_ORDER_INTERPOLATED_BOUNDARY) THEN
               RHO_G_2 = RHOP(IIG,JJG,KKG-1)
               RHO_OTHER_2 = OM_RHOP(IIO,JJO,KKO+1)
            ENDIF
            ZZZ(1:4) = (/RHO_G_2,ONE_D%RHO_G,RHO_OTHER,RHO_OTHER_2/)
            ONE_D%RHO_F = SCALAR_FACE_VALUE(WW(II,JJ,KK-1),ZZZ,FLUX_LIMITER)
      END SELECT

      ! Species and temperature

      SINGLE_SPEC_IF: IF (N_TOTAL_SCALARS > 1) THEN
         SPECIES_LOOP: DO N=1,N_TOTAL_SCALARS

            RHO_ZZ_G = ONE_D%RHO_G*ONE_D%ZZ_G(N)
            RHO_ZZ_G_2 = RHO_ZZ_G ! first-order (default)
            RHO_ZZ_OTHER_2 = RHO_ZZ_OTHER(N)

            SELECT CASE(IOR)
               CASE( 1)
                  IF (SECOND_ORDER_INTERPOLATED_BOUNDARY) THEN
                     RHO_ZZ_G_2 = RHOP(IIG+1,JJG,KKG)*ZZP(IIG+1,JJG,KKG,N)
                     RHO_ZZ_OTHER_2 = OM_RHOP(IIO-1,JJO,KKO)*OM_ZZP(IIO-1,JJO,KKO,N)
                  ENDIF
                  ZZZ(1:4) = (/RHO_ZZ_OTHER_2,RHO_ZZ_OTHER(N),RHO_ZZ_G,RHO_ZZ_G_2/)
                  RHO_ZZ_F(N) = SCALAR_FACE_VALUE(UU(II,JJ,KK),ZZZ,FLUX_LIMITER)
                  PBAR_F = PBAR_P(KKG,ONE_D%PRESSURE_ZONE)
               CASE(-1)
                  IF (SECOND_ORDER_INTERPOLATED_BOUNDARY) THEN
                     RHO_ZZ_G_2 = RHOP(IIG-1,JJG,KKG)*ZZP(IIG-1,JJG,KKG,N)
                     RHO_ZZ_OTHER_2 = OM_RHOP(IIO+1,JJO,KKO)*OM_ZZP(IIO+1,JJO,KKO,N)
                  ENDIF
                  ZZZ(1:4) = (/RHO_ZZ_G_2,RHO_ZZ_G,RHO_ZZ_OTHER(N),RHO_ZZ_OTHER_2/)
                  RHO_ZZ_F(N) = SCALAR_FACE_VALUE(UU(II-1,JJ,KK),ZZZ,FLUX_LIMITER)
                  PBAR_F = PBAR_P(KKG,ONE_D%PRESSURE_ZONE)
               CASE( 2)
                  IF (SECOND_ORDER_INTERPOLATED_BOUNDARY) THEN
                     RHO_ZZ_G_2 = RHOP(IIG,JJG+1,KKG)*ZZP(IIG,JJG+1,KKG,N)
                     RHO_ZZ_OTHER_2 = OM_RHOP(IIO,JJO-1,KKO)*OM_ZZP(IIO,JJO-1,KKO,N)
                  ENDIF
                  ZZZ(1:4) = (/RHO_ZZ_OTHER_2,RHO_ZZ_OTHER(N),RHO_ZZ_G,RHO_ZZ_G_2/)
                  RHO_ZZ_F(N) = SCALAR_FACE_VALUE(VV(II,JJ,KK),ZZZ,FLUX_LIMITER)
                  PBAR_F = PBAR_P(KKG,ONE_D%PRESSURE_ZONE)
               CASE(-2)
                  IF (SECOND_ORDER_INTERPOLATED_BOUNDARY) THEN
                     RHO_ZZ_G_2 = RHOP(IIG,JJG-1,KKG)*ZZP(IIG,JJG-1,KKG,N)
                     RHO_ZZ_OTHER_2 = OM_RHOP(IIO,JJO+1,KKO)*OM_ZZP(IIO,JJO+1,KKO,N)
                  ENDIF
                  ZZZ(1:4) = (/RHO_ZZ_G_2,RHO_ZZ_G,RHO_ZZ_OTHER(N),RHO_ZZ_OTHER_2/)
                  RHO_ZZ_F(N) = SCALAR_FACE_VALUE(VV(II,JJ-1,KK),ZZZ,FLUX_LIMITER)
                  PBAR_F = PBAR_P(KKG,ONE_D%PRESSURE_ZONE)
               CASE( 3)
                  IF (SECOND_ORDER_INTERPOLATED_BOUNDARY) THEN
                     RHO_ZZ_G_2 = RHOP(IIG,JJG,KKG+1)*ZZP(IIG,JJG,KKG+1,N)
                     RHO_ZZ_OTHER_2 = OM_RHOP(IIO,JJO,KKO-1)*OM_ZZP(IIO,JJO,KKO-1,N)
                  ENDIF
                  ZZZ(1:4) = (/RHO_ZZ_OTHER_2,RHO_ZZ_OTHER(N),RHO_ZZ_G,RHO_ZZ_G_2/)
                  RHO_ZZ_F(N) = SCALAR_FACE_VALUE(WW(II,JJ,KK),ZZZ,FLUX_LIMITER)
                  PBAR_F = (PBAR_P(KK,ONE_D%PRESSURE_ZONE)*DZ(KKG) + PBAR_P(KKG,ONE_D%PRESSURE_ZONE)*DZ(KK)) / (DZ(KK)+DZ(KKG))
               CASE(-3)
                  IF (SECOND_ORDER_INTERPOLATED_BOUNDARY) THEN
                     RHO_ZZ_G_2 = RHOP(IIG,JJG,KKG-1)*ZZP(IIG,JJG,KKG-1,N)
                     RHO_ZZ_OTHER_2 = OM_RHOP(IIO,JJO,KKO+1)*OM_ZZP(IIO,JJO,KKO+1,N)
                  ENDIF
                  ZZZ(1:4) = (/RHO_ZZ_G_2,RHO_ZZ_G,RHO_ZZ_OTHER(N),RHO_ZZ_OTHER_2/)
                  RHO_ZZ_F(N) = SCALAR_FACE_VALUE(WW(II,JJ,KK-1),ZZZ,FLUX_LIMITER)
                  PBAR_F = (PBAR_P(KK,ONE_D%PRESSURE_ZONE)*DZ(KKG) + PBAR_P(KKG,ONE_D%PRESSURE_ZONE)*DZ(KK)) / (DZ(KK)+DZ(KKG))
            END SELECT
         ENDDO SPECIES_LOOP

         ! ghost cell value of temperature
         ZZP(II,JJ,KK,1:N_TOTAL_SCALARS) = MAX(0._EB,MIN(1._EB,RHO_ZZ_OTHER(1:N_TOTAL_SCALARS)/RHO_OTHER))
         ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(II,JJ,KK,1:N_TRACKED_SPECIES)
         CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,RSUM(II,JJ,KK))
         TMP(II,JJ,KK) = PBAR_P(KK,ONE_D%PRESSURE_ZONE)/(RSUM(II,JJ,KK)*RHOP(II,JJ,KK))

         ! face value of temperature
         IF (ATMOSPHERIC_INTERPOLATION) THEN
            ONE_D%TMP_F = (TMP(II,JJ,KK)*DZ(KKG) + TMP(IIG,JJG,KKG)*DZ(KK)) / (DZ(KK)+DZ(KKG))
            ONE_D%ZZ_F(1:N_TOTAL_SCALARS) = (ZZP(II,JJ,KK,1:N_TOTAL_SCALARS)*DZ(KKG) + ZZP(IIG,JJG,KKG,1:N_TOTAL_SCALARS)*DZ(KK)) &
                                          / (DZ(KK)+DZ(KKG))
            ZZ_GET(1:N_TRACKED_SPECIES) = ONE_D%ZZ_F(1:N_TRACKED_SPECIES)
            CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,RSUM_F)
            ONE_D%RHO_F = PBAR_F/(RSUM_F*ONE_D%TMP_F)
         ELSE
            ONE_D%ZZ_F(1:N_TOTAL_SCALARS) = MAX(0._EB,MIN(1._EB,RHO_ZZ_F(1:N_TOTAL_SCALARS)/ONE_D%RHO_F))
            ZZ_GET(1:N_TRACKED_SPECIES) = ONE_D%ZZ_F(1:N_TRACKED_SPECIES)
            CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,RSUM_F)
            ONE_D%TMP_F = PBAR_F/(RSUM_F*ONE_D%RHO_F)
         ENDIF

      ELSE SINGLE_SPEC_IF
         ONE_D%ZZ_F(1) = 1._EB
         TMP(II,JJ,KK) = PBAR_P(KK,ONE_D%PRESSURE_ZONE)/(RSUM0*RHOP(II,JJ,KK))
         SELECT CASE(IOR)
            CASE DEFAULT
               PBAR_F = PBAR_P(KKG,ONE_D%PRESSURE_ZONE)
            CASE (-3,3)
               PBAR_F = (PBAR_P(KK,ONE_D%PRESSURE_ZONE)*DZ(KKG) + PBAR_P(KKG,ONE_D%PRESSURE_ZONE)*DZ(KK)) / (DZ(KK)+DZ(KKG))
         END SELECT
         IF (ATMOSPHERIC_INTERPOLATION) THEN
            ONE_D%TMP_F = (TMP(II,JJ,KK)*DZ(KKG) + TMP(IIG,JJG,KKG)*DZ(KK)) / (DZ(KK)+DZ(KKG))
            ONE_D%RHO_F = PBAR_F/(RSUM0*ONE_D%TMP_F)
         ELSE
            ONE_D%TMP_F = PBAR_F/(RSUM0*ONE_D%RHO_F)
         ENDIF
      ENDIF SINGLE_SPEC_IF

      IF (SOLID_HT3D) THEN
         IF (SOLID_OTHER) TMP(II,JJ,KK) = TMP_OTHER
      ENDIF

      ONE_D%QCONF = 0._EB ! no convective heat transfer at interpolated boundary

END SELECT METHOD_OF_HEAT_TRANSFER

END SUBROUTINE CALCULATE_TMP_F


SUBROUTINE SOLID_HEAT_TRANSFER_3D

! Solves the 3D heat conduction equation internal to OBSTs.

REAL(EB) :: DT_SUB,T_LOC,K_S,K_S_M,K_S_P,TMP_G,TMP_F,TMP_S,RDN,HTC,TMP_OTHER,RAMP_FACTOR,&
            QNET,TSI,FDERIV,QEXTRA,K_S_MAX,VN_HT3D,R_K_S,TMP_I,TH_EST4,FO_EST3,&
            RHO_GET(N_MATL),K_GET,K_OTHER,RHOCBAR_S,VC,VSRVC_LOC,RDS,KDTDN_S
INTEGER  :: II,JJ,KK,I,J,K,IOR,IC,ICM,ICP,IIG,JJG,KKG,ADCOUNT,IIO,JJO,KKO,NOM,N_INT_CELLS,NN,ITER
LOGICAL :: CONT_MATL_PROP,IS_STABLE_DT_SUB
INTEGER, PARAMETER :: N_JACOBI_ITERATIONS=1,SURFACE_HEAT_FLUX_MODEL=1
REAL(EB), PARAMETER :: DT_SUB_MIN_HT3D=1.E-9_EB
REAL(EB), POINTER, DIMENSION(:,:,:) :: KDTDX=>NULL(),KDTDY=>NULL(),KDTDZ=>NULL(),TMP_NEW=>NULL(),KP=>NULL(),&
                                       VSRVC_X=>NULL(),VSRVC_Y=>NULL(),VSRVC_Z=>NULL(),VSRVC=>NULL()
TYPE(OBSTRUCTION_TYPE), POINTER :: OB=>NULL(),OBM=>NULL(),OBP=>NULL()
TYPE(MESH_TYPE), POINTER :: OM=>NULL()
TYPE(SURFACE_TYPE), POINTER :: MS=>NULL()

! Initialize verification tests

IF (ICYC==1) THEN
   SELECT CASE(HT3D_TEST)
      CASE(1); CALL CRANK_TEST_1(1)
      CASE(2); CALL CRANK_TEST_1(2)
      CASE(3); CALL CRANK_TEST_1(3)
   END SELECT
ENDIF

KDTDX=>WORK1; KDTDX=0._EB
KDTDY=>WORK2; KDTDY=0._EB
KDTDZ=>WORK3; KDTDZ=0._EB
TMP_NEW=>WORK4
KP=>WORK5; KP=0._EB
VSRVC_X=>WORK6; VSRVC_X=1._EB
VSRVC_Y=>WORK7; VSRVC_Y=1._EB
VSRVC_Z=>WORK8; VSRVC_Z=1._EB
VSRVC  =>WORK9; VSRVC  =1._EB

DT_SUB = DT_BC_HT3D
T_LOC = 0._EB

SUBSTEP_LOOP: DO WHILE ( ABS(T_LOC-DT_BC_HT3D)>TWO_EPSILON_EB )
   DT_SUB  = MIN(DT_SUB,DT_BC_HT3D-T_LOC)
   K_S_MAX = 0._EB
   VN_HT3D = 0._EB

   IS_STABLE_DT_SUB = .FALSE.
   TMP_UPDATE_LOOP: DO WHILE (.NOT.IS_STABLE_DT_SUB)

      TMP_NEW=TMP
      JACOBI_ITERATION_LOOP: DO ITER=1,N_JACOBI_ITERATIONS

         ! compute material thermal conductivity
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  IC = CELL_INDEX(I,J,K);              IF (.NOT.SOLID(IC)) CYCLE
                  OB => OBSTRUCTION(OBST_INDEX_C(IC)); IF (.NOT.OB%HT3D)   CYCLE
                  IF (OB%MATL_INDEX>0) THEN
                     CALL GET_SOLID_CONDUCTIVITY(KP(I,J,K),TMP_NEW(I,J,K),OPT_MATL_INDEX=OB%MATL_INDEX)
                  ELSEIF (OB%MATL_SURF_INDEX>0) THEN
                     MS => SURFACE(OB%MATL_SURF_INDEX)
                     IF (TWO_D) THEN
                        VC = DX(I)*DZ(K)
                     ELSE
                        VC = DX(I)*DY(J)*DZ(K)
                     ENDIF
                     VSRVC(I,J,K) = 0._EB
                     DO NN=1,MS%N_MATL
                        ML => MATERIAL(MS%MATL_INDEX(NN))
                        VSRVC(I,J,K) = VSRVC(I,J,K) + OB%RHO(I,J,K,NN)/ML%RHO_S
                     ENDDO
                     IF (VSRVC(I,J,K)>TWO_EPSILON_EB) THEN
                        RHO_GET(1:MS%N_MATL) = OB%RHO(I,J,K,1:MS%N_MATL) / VSRVC(I,J,K)
                     ELSE
                        RHO_GET(1:MS%N_MATL) = 0._EB
                     ENDIF
                     CALL GET_SOLID_CONDUCTIVITY(KP(I,J,K),TMP_NEW(I,J,K),OPT_SURF_INDEX=OB%MATL_SURF_INDEX,OPT_RHO_IN=RHO_GET)
                     SELECT CASE(ABS(OB%PYRO3D_IOR))
                        CASE DEFAULT
                           ! isotropic shrinking and swelling
                           VSRVC_X(I,J,K) = VSRVC(I,J,K)**ONTH
                           VSRVC_Y(I,J,K) = VSRVC_X(I,J,K)
                           VSRVC_Z(I,J,K) = VSRVC_X(I,J,K)
                        CASE(1)
                           VSRVC_X(I,J,K) = VSRVC(I,J,K)
                           VSRVC_Y(I,J,K) = 1._EB
                           VSRVC_Z(I,J,K) = 1._EB
                        CASE(2)
                           VSRVC_X(I,J,K) = 1._EB
                           VSRVC_Y(I,J,K) = VSRVC(I,J,K)
                           VSRVC_Z(I,J,K) = 1._EB
                        CASE(3)
                           VSRVC_X(I,J,K) = 1._EB
                           VSRVC_Y(I,J,K) = 1._EB
                           VSRVC_Z(I,J,K) = VSRVC(I,J,K)
                     END SELECT
                  ENDIF
               ENDDO
            ENDDO
         ENDDO

         KP_WALL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS
            WC => WALL(IW)
            IF (WC%BOUNDARY_TYPE/=NULL_BOUNDARY) CYCLE KP_WALL_LOOP

            II = WC%ONE_D%II
            JJ = WC%ONE_D%JJ
            KK = WC%ONE_D%KK

            EWC=>EXTERNAL_WALL(IW)
            NOM=EWC%NOM
            IF (NOM<1) CYCLE KP_WALL_LOOP
            OM=>MESHES(NOM)

            K_OTHER = 0._EB
            DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
               DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
                  EWC_IIO_LOOP: DO IIO=EWC%IIO_MIN,EWC%IIO_MAX

                     IC = OM%CELL_INDEX(IIO,JJO,KKO);           IF (.NOT.OM%SOLID(IC)) CYCLE EWC_IIO_LOOP
                     OB => OM%OBSTRUCTION(OM%OBST_INDEX_C(IC)); IF (.NOT.OB%HT3D) CYCLE EWC_IIO_LOOP
                     TMP_OTHER = OMESH(NOM)%TMP(IIO,JJO,KKO)

                     K_GET = 0._EB
                     IF (OB%MATL_INDEX>0) THEN
                        CALL GET_SOLID_CONDUCTIVITY(K_GET,TMP_OTHER,OPT_MATL_INDEX=OB%MATL_INDEX)
                     ELSEIF (OB%MATL_SURF_INDEX>0) THEN
                        VSRVC_LOC = 0._EB
                        DO NN=1,SURFACE(OB%MATL_SURF_INDEX)%N_MATL
                           ML => MATERIAL(SURFACE(OB%MATL_SURF_INDEX)%MATL_INDEX(NN))
                           VSRVC_LOC = VSRVC_LOC + OB%RHO(IIO,JJO,KKO,NN)/ML%RHO_S
                        ENDDO
                        RHO_GET(1:SURFACE(OB%MATL_SURF_INDEX)%N_MATL) = OB%RHO(IIO,JJO,KKO,1:SURFACE(OB%MATL_SURF_INDEX)%N_MATL)&
                                                                      / VSRVC_LOC
                        CALL GET_SOLID_CONDUCTIVITY(K_GET,TMP_OTHER,OPT_SURF_INDEX=OB%MATL_SURF_INDEX,OPT_RHO_IN=RHO_GET)
                     ENDIF
                     K_OTHER = K_OTHER + K_GET

                  ENDDO EWC_IIO_LOOP
               ENDDO
            ENDDO
            N_INT_CELLS = (EWC%IIO_MAX-EWC%IIO_MIN+1) * (EWC%JJO_MAX-EWC%JJO_MIN+1) * (EWC%KKO_MAX-EWC%KKO_MIN+1)
            KP(II,JJ,KK) = K_OTHER/REAL(N_INT_CELLS,EB)

         ENDDO KP_WALL_LOOP

         ! build heat flux vectors
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=0,IBAR
                  ICM = CELL_INDEX(I,J,K)
                  ICP = CELL_INDEX(I+1,J,K)
                  IF (.NOT.(SOLID(ICM).AND.SOLID(ICP))) CYCLE

                  OBM => OBSTRUCTION(OBST_INDEX_C(ICM))
                  OBP => OBSTRUCTION(OBST_INDEX_C(ICP))
                  ! At present OBST_INDEX_C is not defined for ghost cells.
                  ! This means that:
                  !    1. continuous material properties will be assumed at a mesh boundary
                  !    2. we assume that if either OBM%HT3D .OR. OBP%HT3D we should process the boundary
                  IF (.NOT.(OBM%HT3D.OR.OBP%HT3D)) CYCLE

                  K_S_M = KP(I,J,K)
                  K_S_P = KP(I+1,J,K)

                  IF (K_S_M<TWO_EPSILON_EB .OR. K_S_P<TWO_EPSILON_EB) THEN
                     KDTDX(I,J,K) = 0._EB
                     CYCLE
                  ENDIF

                  ! determine if we have continuous material properties
                  CONT_MATL_PROP=.TRUE.
                  IF (OBM%MATL_INDEX>0 .AND. OBP%MATL_INDEX>0 .AND. OBM%MATL_INDEX/=OBP%MATL_INDEX) THEN
                     CONT_MATL_PROP=.FALSE.
                  ELSEIF (OBM%MATL_SURF_INDEX>0 .AND. OBP%MATL_SURF_INDEX>0 .AND. OBM%MATL_SURF_INDEX/=OBP%MATL_SURF_INDEX) THEN
                     CONT_MATL_PROP=.FALSE.
                  ELSEIF (OBM%MATL_INDEX>0 .AND. OBP%MATL_SURF_INDEX>0) THEN
                     CONT_MATL_PROP=.FALSE.
                  ELSEIF (OBM%MATL_SURF_INDEX>0 .AND. OBP%MATL_INDEX>0) THEN
                     CONT_MATL_PROP=.FALSE.
                  ENDIF

                  IF (CONT_MATL_PROP) THEN
                     ! use linear average from inverse lever rule
                     K_S = ( K_S_M*DX(I+1) + K_S_P*DX(I) )/( DX(I) + DX(I+1) )
                     K_S_MAX = MAX(K_S_MAX,K_S)
                     KDTDX(I,J,K) = K_S*(TMP_NEW(I+1,J,K)-TMP_NEW(I,J,K))*2._EB/(DX(I+1)*VSRVC_X(I+1,J,K)+DX(I)*VSRVC_X(I,J,K))
                  ELSE
                     ! for discontinuous material properties maintain continuity of flux, C0 continuity of temperature
                     ! (allow C1 discontinuity of temperature due to jump in thermal properties across interface)
                     R_K_S = K_S_P/K_S_M * DX(I)/DX(I+1) * VSRVC_X(I,J,K)/VSRVC_X(I+1,J,K)
                     TMP_I = (TMP_NEW(I,J,K) + R_K_S*TMP_NEW(I+1,J,K))/(1._EB + R_K_S) ! interface temperature
                     !! KDTDX(I,J,K) = K_S_P * (TMP_NEW(I+1,J,K)-TMP_I) * 2._EB/(DX(I+1)*VSRVC_X(I+1,J,K)) !! should be identical
                     KDTDX(I,J,K) = K_S_M * (TMP_I-TMP_NEW(I,J,K)) * 2._EB/(DX(I)*VSRVC_X(I,J,K))
                     K_S_MAX = MAX(K_S_MAX,MAX(K_S_M,K_S_P))
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
         TWO_D_IF: IF (.NOT.TWO_D) THEN
            DO K=1,KBAR
               DO J=0,JBAR
                  DO I=1,IBAR
                     ICM = CELL_INDEX(I,J,K)
                     ICP = CELL_INDEX(I,J+1,K)
                     IF (.NOT.(SOLID(ICM).AND.SOLID(ICP))) CYCLE
                     OBM => OBSTRUCTION(OBST_INDEX_C(ICM))
                     OBP => OBSTRUCTION(OBST_INDEX_C(ICP))
                     IF (.NOT.(OBM%HT3D.OR.OBP%HT3D)) CYCLE

                     K_S_M = KP(I,J,K)
                     K_S_P = KP(I,J+1,K)

                     IF (K_S_M<TWO_EPSILON_EB .OR. K_S_P<TWO_EPSILON_EB) THEN
                        KDTDY(I,J,K) = 0._EB
                        CYCLE
                     ENDIF

                     CONT_MATL_PROP=.TRUE.
                     IF (OBM%MATL_INDEX>0 .AND. OBP%MATL_INDEX>0 .AND. OBM%MATL_INDEX/=OBP%MATL_INDEX) THEN
                        CONT_MATL_PROP=.FALSE.
                     ELSEIF (OBM%MATL_SURF_INDEX>0 .AND. OBP%MATL_SURF_INDEX>0 .AND. OBM%MATL_SURF_INDEX/=OBP%MATL_SURF_INDEX) THEN
                        CONT_MATL_PROP=.FALSE.
                     ELSEIF (OBM%MATL_INDEX>0 .AND. OBP%MATL_SURF_INDEX>0) THEN
                        CONT_MATL_PROP=.FALSE.
                     ELSEIF (OBM%MATL_SURF_INDEX>0 .AND. OBP%MATL_INDEX>0) THEN
                        CONT_MATL_PROP=.FALSE.
                     ENDIF

                     IF (CONT_MATL_PROP) THEN
                        K_S = ( K_S_M*DY(J+1) + K_S_P*DY(J) )/( DY(J) + DY(J+1) )
                        K_S_MAX = MAX(K_S_MAX,K_S)
                        KDTDY(I,J,K) = K_S*(TMP_NEW(I,J+1,K)-TMP_NEW(I,J,K))*2._EB/(DY(J+1)*VSRVC_Y(I,J+1,K)+DY(J)*VSRVC_Y(I,J,K))
                     ELSE
                        R_K_S = K_S_P/K_S_M * DY(J)/DY(J+1) * VSRVC_Y(I,J,K)/VSRVC_Y(I,J+1,K)
                        TMP_I = (TMP_NEW(I,J,K) + R_K_S*TMP_NEW(I,J+1,K))/(1._EB + R_K_S)
                        KDTDY(I,J,K) = K_S_M * (TMP_I-TMP_NEW(I,J,K)) * 2._EB/(DY(J)*VSRVC_Y(I,J,K))
                        K_S_MAX = MAX(K_S_MAX,MAX(K_S_M,K_S_P))
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ELSE TWO_D_IF
            KDTDY(I,J,K) = 0._EB
         ENDIF TWO_D_IF
         DO K=0,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  ICM = CELL_INDEX(I,J,K)
                  ICP = CELL_INDEX(I,J,K+1)
                  IF (.NOT.(SOLID(ICM).AND.SOLID(ICP))) CYCLE
                  OBM => OBSTRUCTION(OBST_INDEX_C(ICM))
                  OBP => OBSTRUCTION(OBST_INDEX_C(ICP))
                  IF (.NOT.(OBM%HT3D.OR.OBP%HT3D)) CYCLE

                  K_S_M = KP(I,J,K)
                  K_S_P = KP(I,J,K+1)

                  IF (K_S_M<TWO_EPSILON_EB .OR. K_S_P<TWO_EPSILON_EB) THEN
                     KDTDZ(I,J,K) = 0._EB
                     CYCLE
                  ENDIF

                  CONT_MATL_PROP=.TRUE.
                  IF (OBM%MATL_INDEX>0 .AND. OBP%MATL_INDEX>0 .AND. OBM%MATL_INDEX/=OBP%MATL_INDEX) THEN
                     CONT_MATL_PROP=.FALSE.
                  ELSEIF (OBM%MATL_SURF_INDEX>0 .AND. OBP%MATL_SURF_INDEX>0 .AND. OBM%MATL_SURF_INDEX/=OBP%MATL_SURF_INDEX) THEN
                     CONT_MATL_PROP=.FALSE.
                  ELSEIF (OBM%MATL_INDEX>0 .AND. OBP%MATL_SURF_INDEX>0) THEN
                     CONT_MATL_PROP=.FALSE.
                  ELSEIF (OBM%MATL_SURF_INDEX>0 .AND. OBP%MATL_INDEX>0) THEN
                     CONT_MATL_PROP=.FALSE.
                  ENDIF

                  IF (CONT_MATL_PROP) THEN
                     K_S = ( K_S_M*DZ(K+1) + K_S_P*DZ(K) )/( DZ(K) + DZ(K+1) )
                     K_S_MAX = MAX(K_S_MAX,K_S)
                     KDTDZ(I,J,K) = K_S*(TMP_NEW(I,J,K+1)-TMP_NEW(I,J,K))*2._EB/(DZ(K+1)*VSRVC_Z(I,J,K+1)+DZ(K)*VSRVC_Z(I,J,K))
                  ELSE
                     R_K_S = K_S_P/K_S_M * DZ(K)/DZ(K+1) * VSRVC_Z(I,J,K)/VSRVC_Z(I,J,K+1)
                     TMP_I = (TMP_NEW(I,J,K) + R_K_S*TMP_NEW(I,J,K+1))/(1._EB + R_K_S)
                     KDTDZ(I,J,K) = K_S_M * (TMP_I-TMP_NEW(I,J,K)) * 2._EB/(DZ(K)*VSRVC_Z(I,J,K))
                     K_S_MAX = MAX(K_S_MAX,MAX(K_S_M,K_S_P))
                  ENDIF
               ENDDO
            ENDDO
         ENDDO

         ! build fluxes on boundaries of INTERNAL WALL CELLS

         HT3D_WALL_LOOP: DO IW=N_EXTERNAL_WALL_CELLS+1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
            WC => WALL(IW)
            IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE HT3D_WALL_LOOP

            SURF_INDEX = WC%SURF_INDEX
            SF => SURFACE(SURF_INDEX)
            II = WC%ONE_D%II
            JJ = WC%ONE_D%JJ
            KK = WC%ONE_D%KK
            IIG = WC%ONE_D%IIG
            JJG = WC%ONE_D%JJG
            KKG = WC%ONE_D%KKG
            IOR = WC%ONE_D%IOR

            IC = CELL_INDEX(II,JJ,KK);           IF (.NOT.SOLID(IC)) CYCLE HT3D_WALL_LOOP
            OB => OBSTRUCTION(OBST_INDEX_C(IC)); IF (.NOT.OB%HT3D  ) CYCLE HT3D_WALL_LOOP

            MATL_IF: IF (OB%MATL_INDEX>0) THEN
               CALL GET_SOLID_CONDUCTIVITY(K_S,WC%ONE_D%TMP_F,OPT_MATL_INDEX=OB%MATL_INDEX)
            ELSEIF (OB%MATL_SURF_INDEX>0) THEN
               MS => SURFACE(OB%MATL_SURF_INDEX)
               IF (VSRVC(II,JJ,KK)>TWO_EPSILON_EB) THEN
                  RHO_GET(1:MS%N_MATL) = OB%RHO(II,JJ,KK,1:MS%N_MATL) / VSRVC(II,JJ,KK)
               ELSE
                  RHO_GET(1:MS%N_MATL) = 0._EB
               ENDIF
               CALL GET_SOLID_CONDUCTIVITY(K_S,WC%ONE_D%TMP_F,OPT_SURF_INDEX=OB%MATL_SURF_INDEX,OPT_RHO_IN=RHO_GET)
            ENDIF MATL_IF
            K_S_MAX = MAX(K_S_MAX,K_S)

            METHOD_OF_HEAT_TRANSFER: SELECT CASE(SF%THERMAL_BC_INDEX)

               CASE DEFAULT METHOD_OF_HEAT_TRANSFER ! includes SPECIFIED_TEMPERATURE

                  SELECT CASE(IOR)
                     CASE( 1); KDTDX(II,JJ,KK)   = K_S * 2._EB*(WC%ONE_D%TMP_F-TMP_NEW(II,JJ,KK))*RDX(II)/VSRVC_X(II,JJ,KK)
                     CASE(-1); KDTDX(II-1,JJ,KK) = K_S * 2._EB*(TMP_NEW(II,JJ,KK)-WC%ONE_D%TMP_F)*RDX(II)/VSRVC_X(II,JJ,KK)
                     CASE( 2); KDTDY(II,JJ,KK)   = K_S * 2._EB*(WC%ONE_D%TMP_F-TMP_NEW(II,JJ,KK))*RDY(JJ)/VSRVC_Y(II,JJ,KK)
                     CASE(-2); KDTDY(II,JJ-1,KK) = K_S * 2._EB*(TMP_NEW(II,JJ,KK)-WC%ONE_D%TMP_F)*RDY(JJ)/VSRVC_Y(II,JJ,KK)
                     CASE( 3); KDTDZ(II,JJ,KK)   = K_S * 2._EB*(WC%ONE_D%TMP_F-TMP_NEW(II,JJ,KK))*RDZ(KK)/VSRVC_Z(II,JJ,KK)
                     CASE(-3); KDTDZ(II,JJ,KK-1) = K_S * 2._EB*(TMP_NEW(II,JJ,KK)-WC%ONE_D%TMP_F)*RDZ(KK)/VSRVC_Z(II,JJ,KK)
                  END SELECT

               CASE (NET_FLUX_BC) METHOD_OF_HEAT_TRANSFER
                  SELECT CASE(IOR)
                     CASE( 1); KDTDX(II,JJ,KK)   = -SF%NET_HEAT_FLUX*WC%ONE_D%AREA_ADJUST
                     CASE(-1); KDTDX(II-1,JJ,KK) =  SF%NET_HEAT_FLUX*WC%ONE_D%AREA_ADJUST
                     CASE( 2); KDTDY(II,JJ,KK)   = -SF%NET_HEAT_FLUX*WC%ONE_D%AREA_ADJUST
                     CASE(-2); KDTDY(II,JJ-1,KK) =  SF%NET_HEAT_FLUX*WC%ONE_D%AREA_ADJUST
                     CASE( 3); KDTDZ(II,JJ,KK)   = -SF%NET_HEAT_FLUX*WC%ONE_D%AREA_ADJUST
                     CASE(-3); KDTDZ(II,JJ,KK-1) =  SF%NET_HEAT_FLUX*WC%ONE_D%AREA_ADJUST
                  END SELECT

                  SOLID_PHASE_ONLY_IF: IF (SOLID_PHASE_ONLY) THEN
                     SELECT CASE(IOR)
                        CASE( 1); WC%ONE_D%TMP_F = TMP_NEW(II,JJ,KK) + KDTDX(II,JJ,KK)   / (K_S*2._EB*RDX(II)/VSRVC_X(II,JJ,KK))
                        CASE(-1); WC%ONE_D%TMP_F = TMP_NEW(II,JJ,KK) - KDTDX(II-1,JJ,KK) / (K_S*2._EB*RDX(II)/VSRVC_X(II,JJ,KK))
                        CASE( 2); WC%ONE_D%TMP_F = TMP_NEW(II,JJ,KK) + KDTDY(II,JJ,KK)   / (K_S*2._EB*RDY(JJ)/VSRVC_Y(II,JJ,KK))
                        CASE(-2); WC%ONE_D%TMP_F = TMP_NEW(II,JJ,KK) - KDTDY(II,JJ-1,KK) / (K_S*2._EB*RDY(JJ)/VSRVC_Y(II,JJ,KK))
                        CASE( 3); WC%ONE_D%TMP_F = TMP_NEW(II,JJ,KK) + KDTDZ(II,JJ,KK)   / (K_S*2._EB*RDZ(KK)/VSRVC_Z(II,JJ,KK))
                        CASE(-3); WC%ONE_D%TMP_F = TMP_NEW(II,JJ,KK) - KDTDZ(II,JJ,KK-1) / (K_S*2._EB*RDZ(KK)/VSRVC_Z(II,JJ,KK))
                     END SELECT
                  ELSE
                     TMP_G = TMP_NEW(IIG,JJG,KKG)
                     TMP_F = WC%ONE_D%TMP_F
                     TMP_OTHER = TMP_F
                     DTMP = TMP_G - TMP_F
                     IF (ABS(ONE_D%T_IGN-T_BEGIN)<= SPACING(ONE_D%T_IGN) .AND. SF%RAMP_INDEX(TIME_HEAT)>=1) THEN
                        TSI = T
                     ELSE
                        TSI = T - ONE_D%T_IGN
                     ENDIF
                     RAMP_FACTOR = EVALUATE_RAMP(TSI,SF%TAU(TIME_HEAT),SF%RAMP_INDEX(TIME_HEAT))
                     QNET = -RAMP_FACTOR*SF%NET_HEAT_FLUX*ONE_D%AREA_ADJUST
                     ADCOUNT = 0
                     ADLOOP: DO
                        ADCOUNT = ADCOUNT + 1
                        DTMP = TMP_G - TMP_OTHER
                        IF (ABS(QNET) > 0._EB .AND. ABS(DTMP) <TWO_EPSILON_EB) DTMP=1._EB
                        WC%ONE_D%HEAT_TRANS_COEF = HEAT_TRANSFER_COEFFICIENT(DTMP,SF%H_FIXED,SURF_INDEX,WALL_INDEX=IW)
                        HTC = WC%ONE_D%HEAT_TRANS_COEF
                        IF (RADIATION) THEN
                           QEXTRA = HTC*DTMP + WC%ONE_D%QRADIN - WC%ONE_D%EMISSIVITY*SIGMA*TMP_OTHER**4 - QNET
                           FDERIV = -HTC - 4._EB*WC%ONE_D%EMISSIVITY*SIGMA*TMP_OTHER**3
                        ELSE
                           QEXTRA = HTC*DTMP - QNET
                           FDERIV = -HTC
                        ENDIF
                        IF (ABS(FDERIV) > TWO_EPSILON_EB) TMP_OTHER = TMP_OTHER - QEXTRA / FDERIV
                        IF (ABS(TMP_OTHER - TMP_F) / TMP_F < 1.E-4_EB .OR. ADCOUNT > 20) THEN
                           TMP_F = MIN(TMPMAX,TMP_OTHER)
                           EXIT ADLOOP
                        ELSE
                           TMP_F = MIN(TMPMAX,TMP_OTHER)
                           CYCLE ADLOOP
                        ENDIF
                     ENDDO ADLOOP
                     WC%ONE_D%TMP_F = TMP_F
                     WC%ONE_D%QCONF = HTC*DTMP
                  ENDIF SOLID_PHASE_ONLY_IF

               CASE (THERMALLY_THICK_HT3D) ! thermally thick, continuous heat flux

                  IIG = WC%ONE_D%IIG
                  JJG = WC%ONE_D%JJG
                  KKG = WC%ONE_D%KKG
                  TMP_G = TMP_NEW(IIG,JJG,KKG)
                  TMP_S = TMP_NEW(II,JJ,KK)
                  TMP_F = WC%ONE_D%TMP_F
                  RDS = 0._EB
                  TMP_F_LOOP: DO ADCOUNT=1,3
                     SELECT CASE(ABS(IOR))
                        CASE( 1); RDN = MAX( RDS, RDX(II)/VSRVC_X(II,JJ,KK) )
                        CASE( 2); RDN = MAX( RDS, RDY(JJ)/VSRVC_Y(II,JJ,KK) )
                        CASE( 3); RDN = MAX( RDS, RDZ(KK)/VSRVC_Z(II,JJ,KK) )
                     END SELECT
                     DTMP = TMP_G - TMP_F
                     WC%ONE_D%HEAT_TRANS_COEF = HEAT_TRANSFER_COEFFICIENT(DTMP,SF%H_FIXED,SURF_INDEX,WALL_INDEX=IW)
                     HTC = WC%ONE_D%HEAT_TRANS_COEF
                     IF (RADIATION) THEN
                        TH_EST4 = 3._EB*WC%ONE_D%EMISSIVITY*SIGMA*TMP_F**4
                        FO_EST3 = 4._EB*WC%ONE_D%EMISSIVITY*SIGMA*TMP_F**3
                        TMP_F = ( WC%ONE_D%QRADIN + TH_EST4 + HTC*TMP_G + 2._EB*K_S*RDN*TMP_S ) / &
                                (                   FO_EST3 + HTC       + 2._EB*K_S*RDN       )
                     ELSE
                        TMP_F = ( HTC*TMP_G + 2._EB*K_S*RDN*TMP_S ) / &
                                ( HTC       + 2._EB*K_S*RDN       )
                     ENDIF
                     IF (OB%MATL_INDEX>0) THEN
                        CALL GET_SOLID_RHOCBAR(RHOCBAR_S,TMP_S,OPT_MATL_INDEX=OB%MATL_INDEX)
                     ELSEIF (OB%MATL_SURF_INDEX>0) THEN
                        CALL GET_SOLID_RHOCBAR(RHOCBAR_S,TMP_S,OPT_SURF_INDEX=OB%MATL_SURF_INDEX,OPT_RHO_IN=RHO_GET)
                     ENDIF
                     SELECT CASE(SURFACE_HEAT_FLUX_MODEL)
                        CASE DEFAULT
                           RDS = 0._EB
                        CASE(1)
                           ! FDS Tech Guide (M.3), gives same length scale as 1D pyrolysis model
                           RDS = SQRT(RHOCBAR_S/K_S)
                        CASE(2)
                           ! experimental wall model, generally gives smaller length scale than 1D pyro model
                           KDTDN_S = ABS(K_S*2._EB*(TMP_F-TMP_S)*RDN)
                           RDS = 0.5_EB * ( KDTDN_S / (K_S/RHOCBAR_S)**3 / SUM(RHO_GET(1:MS%N_MATL)) )**ONTH
                     END SELECT
                  ENDDO TMP_F_LOOP
                  WC%ONE_D%TMP_F = TMP_F
                  WC%ONE_D%QCONF = HTC*(TMP_G-TMP_F)

                  SELECT CASE(IOR)
                     CASE( 1); KDTDX(II,JJ,KK)   = K_S * 2._EB*(TMP_F-TMP_S)*RDN
                     CASE(-1); KDTDX(II-1,JJ,KK) = K_S * 2._EB*(TMP_S-TMP_F)*RDN
                     CASE( 2); KDTDY(II,JJ,KK)   = K_S * 2._EB*(TMP_F-TMP_S)*RDN
                     CASE(-2); KDTDY(II,JJ-1,KK) = K_S * 2._EB*(TMP_S-TMP_F)*RDN
                     CASE( 3); KDTDZ(II,JJ,KK)   = K_S * 2._EB*(TMP_F-TMP_S)*RDN
                     CASE(-3); KDTDZ(II,JJ,KK-1) = K_S * 2._EB*(TMP_S-TMP_F)*RDN
                  END SELECT

            END SELECT METHOD_OF_HEAT_TRANSFER

         ENDDO HT3D_WALL_LOOP

         ! Note: for 2D cylindrical KDTDX at X=0 remains zero after initialization

         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  IC = CELL_INDEX(I,J,K)
                  IF (.NOT.SOLID(IC)) CYCLE
                  OB => OBSTRUCTION(OBST_INDEX_C(IC)); IF (.NOT.OB%HT3D) CYCLE
                  IF (OB%MATL_INDEX>0) THEN
                     CALL GET_SOLID_RHOCBAR(RHOCBAR_S,TMP_NEW(I,J,K),OPT_MATL_INDEX=OB%MATL_INDEX)
                  ELSEIF (OB%MATL_SURF_INDEX>0) THEN
                     MS => SURFACE(OB%MATL_SURF_INDEX)
                     RHO_GET(1:MS%N_MATL) = OB%RHO(I,J,K,1:MS%N_MATL)
                     CALL GET_SOLID_RHOCBAR(RHOCBAR_S,TMP_NEW(I,J,K),OPT_SURF_INDEX=OB%MATL_SURF_INDEX,OPT_RHO_IN=RHO_GET)
                  ENDIF
                  IF (TWO_D) THEN
                     VN_HT3D = MAX( VN_HT3D, 2._EB*K_S_MAX/RHOCBAR_S*( RDX(I)**2 + RDZ(K)**2 ) )
                  ELSE
                     VN_HT3D = MAX( VN_HT3D, 2._EB*K_S_MAX/RHOCBAR_S*( RDX(I)**2 + RDY(J)**2 + RDZ(K)**2 ) )
                  ENDIF

                  TMP_NEW(I,J,K) = TMP(I,J,K) + DT_SUB/RHOCBAR_S * ( (KDTDX(I,J,K)*R(I)-KDTDX(I-1,J,K)*R(I-1))*RDX(I)*RRN(I) + &
                                                                     (KDTDY(I,J,K)     -KDTDY(I,J-1,K)       )*RDY(J) + &
                                                                     (KDTDZ(I,J,K)     -KDTDZ(I,J,K-1)       )*RDZ(K) + &
                                                                     Q(I,J,K) + Q_DOT_PPP_S(I,J,K) )

                  TMP_NEW(I,J,K) = MIN(TMPMAX,MAX(TMPMIN,TMP_NEW(I,J,K)))
               ENDDO
            ENDDO
         ENDDO

      ENDDO JACOBI_ITERATION_LOOP

      ! time step adjustment

      IF (DT_SUB*VN_HT3D < VN_MAX .OR. LOCK_TIME_STEP) THEN
         IS_STABLE_DT_SUB = .TRUE.
         TMP = TMP_NEW
         IF (SOLID_PYRO3D) CALL SOLID_PYROLYSIS_3D(DT_SUB,T_LOC)
         T_LOC = T_LOC + DT_SUB
         IF (.NOT.LOCK_TIME_STEP) DT_SUB = MAX( DT_SUB, VN_MIN / MAX(VN_HT3D,TWO_EPSILON_EB) )
      ELSE
         DT_SUB = 0.5_EB*(VN_MIN+VN_MAX) / MAX(VN_HT3D,TWO_EPSILON_EB)
      ENDIF
      IF (DT_SUB < DT_SUB_MIN_HT3D .AND. (T+DT_SUB < (T_END-TWO_EPSILON_EB))) THEN
         WRITE(LU_ERR,'(A)') 'HT3D Instability: DT_SUB < 1e-9 s'
         STOP_STATUS = INSTABILITY_STOP
         RETURN
      ENDIF

   ENDDO TMP_UPDATE_LOOP

ENDDO SUBSTEP_LOOP

END SUBROUTINE SOLID_HEAT_TRANSFER_3D


SUBROUTINE SOLID_PYROLYSIS_3D(DT_SUB,T_LOC)

REAL(EB), INTENT(IN) :: DT_SUB,T_LOC
INTEGER :: N,NN,NS,I,J,K,IC,IIG,JJG,KKG,II2,JJ2,KK2,IOR,OBST_INDEX
REAL(EB) :: M_DOT_G_PPP_ADJUST(N_TRACKED_SPECIES),M_DOT_G_PPP_ACTUAL(N_TRACKED_SPECIES),M_DOT_S_PPP(MAX_MATERIALS),&
            RHO_IN(N_MATL),RHO_OUT(N_MATL),GEOM_FACTOR,TIME_FACTOR,VC,VC2,TMP_S,VSRVC_LOC,RHOCBAR,RHOCBAR2,DUMMY
LOGICAL :: OB2_FOUND
REAL(EB), PARAMETER :: SOLID_VOLUME_MERGE_THRESHOLD=0.1_EB, SOLID_VOLUME_CLIP_THRESHOLD=1.E-6_EB
TYPE(OBSTRUCTION_TYPE), POINTER :: OB=>NULL(),OB2=>NULL()
TYPE(SURFACE_TYPE), POINTER :: SF=>NULL(),MS=>NULL()
TYPE(WALL_TYPE), POINTER :: WC=>NULL()

TIME_FACTOR = DT_SUB/DT_BC_HT3D

INIT_IF: IF (T_LOC<TWO_EPSILON_EB) THEN
   OBST_LOOP_1: DO N=1,N_OBST
      OB => OBSTRUCTION(N)
      IF (.NOT.OB%PYRO3D) CYCLE OBST_LOOP_1
      ! Set mass fluxes to 0
      DO K=OB%K1+1,OB%K2
         DO J=OB%J1+1,OB%J2
            I_LOOP: DO I=OB%I1+1,OB%I2
               IC = CELL_INDEX(I,J,K)
               IF (.NOT.SOLID(IC)) CYCLE I_LOOP
               WC=>WALL(WALL_INDEX_HT3D(IC,OB%PYRO3D_IOR))
               IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE I_LOOP
               SF=>SURFACE(WC%SURF_INDEX)
               WC%ONE_D%MASSFLUX(1:N_TRACKED_SPECIES)      = 0._EB
               WC%ONE_D%MASSFLUX_SPEC(1:N_TRACKED_SPECIES) = 0._EB
               WC%ONE_D%MASSFLUX_MATL(1:SF%N_MATL)         = 0._EB
            ENDDO I_LOOP
         ENDDO
      ENDDO
   ENDDO OBST_LOOP_1
ENDIF INIT_IF

OBST_LOOP_2: DO N=1,N_OBST
   OB => OBSTRUCTION(N)
   IF (.NOT.OB%PYRO3D) CYCLE OBST_LOOP_2

   DO K=OB%K1+1,OB%K2
      DO J=OB%J1+1,OB%J2
         I_LOOP_2: DO I=OB%I1+1,OB%I2
            IC = CELL_INDEX(I,J,K)
            IF (.NOT.SOLID(IC)) CYCLE I_LOOP_2

            WC=>WALL(WALL_INDEX_HT3D(IC,OB%PYRO3D_IOR))
            IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE I_LOOP_2

            SF=>SURFACE(WC%SURF_INDEX)      ! PYROLYSIS SURFACE (ejection of pyrolyzate gas)
            MS=>SURFACE(OB%MATL_SURF_INDEX) ! MATERIAL SURFACE (supplies material properties)

            IIG = WC%ONE_D%IIG
            JJG = WC%ONE_D%JJG
            KKG = WC%ONE_D%KKG
            IOR = WC%ONE_D%IOR
            TMP_S = TMP(I,J,K)
            SELECT CASE(ABS(IOR))
               CASE(1); GEOM_FACTOR = DX(I)
               CASE(2); GEOM_FACTOR = DY(J)
               CASE(3); GEOM_FACTOR = DZ(K)
            END SELECT

            ! cell volume
            IF (TWO_D) THEN
               VC = DX(I)*DZ(K)
            ELSE
               VC = DX(I)*DY(J)*DZ(K)
            ENDIF

            ! update density
            RHO_IN(1:MS%N_MATL) = OB%RHO(I,J,K,1:MS%N_MATL)
            RHO_OUT(1:MS%N_MATL) = RHO_IN(1:MS%N_MATL)

            CALL PYROLYSIS(MS%N_MATL,MS%MATL_INDEX,OB%MATL_SURF_INDEX,IIG,JJG,KKG,TMP_S,DUMMY,&
                           RHO_OUT(1:MS%N_MATL),MS%LAYER_DENSITY(1),DUMMY,DT_SUB,&
                           M_DOT_G_PPP_ADJUST,M_DOT_G_PPP_ACTUAL,M_DOT_S_PPP,Q_DOT_PPP_S(I,J,K))

            OB%RHO(I,J,K,1:MS%N_MATL) = OB%RHO(I,J,K,1:MS%N_MATL) + RHO_OUT(1:MS%N_MATL)-RHO_IN(1:MS%N_MATL)

            ! simple model (no transport): pyrolyzed mass is ejected via wall cell index WALL_INDEX_HT3D(IC,OB%PYRO3D_IOR)

            DO NS=1,N_TRACKED_SPECIES
               WC%ONE_D%MASSFLUX(NS)      = WC%ONE_D%MASSFLUX(NS)      + M_DOT_G_PPP_ADJUST(NS)*GEOM_FACTOR*TIME_FACTOR
               WC%ONE_D%MASSFLUX_SPEC(NS) = WC%ONE_D%MASSFLUX_SPEC(NS) + M_DOT_G_PPP_ACTUAL(NS)*GEOM_FACTOR*TIME_FACTOR
            ENDDO
            DO NN=1,SF%N_MATL
               WC%ONE_D%MASSFLUX_MATL(NN) = WC%ONE_D%MASSFLUX_MATL(NN) + M_DOT_S_PPP(NN)*GEOM_FACTOR*TIME_FACTOR
            ENDDO

            ! If the fuel or water massflux is non-zero, set the ignition time

            IF (WC%ONE_D%T_IGN > T) THEN
               IF (SUM(WC%ONE_D%MASSFLUX(1:N_TRACKED_SPECIES)) > 0._EB) WC%ONE_D%T_IGN = T
            ENDIF

            CONSUMABLE_IF: IF (OB%CONSUMABLE) THEN

               ! recompute solid volume ratio, VS/VC, for cell (I,J,K)
               VSRVC_LOC = 0._EB
               DO NN=1,MS%N_MATL
                  ML => MATERIAL(MS%MATL_INDEX(NN))
                  VSRVC_LOC = VSRVC_LOC + OB%RHO(I,J,K,NN)/ML%RHO_S
               ENDDO

               ! if local cell volume becomes too small, put the mass in the adjacent cell and remove local cell
               THRESHOLD_IF: IF (VSRVC_LOC<SOLID_VOLUME_MERGE_THRESHOLD) THEN
                  OB2_FOUND = .FALSE.
                  II2 = I
                  JJ2 = J
                  KK2 = K
                  ! first, see if cell exists in the minus IOR direction
                  SELECT CASE(IOR)
                     CASE ( 1); II2=I-1
                     CASE (-1); II2=I+1
                     CASE ( 2); JJ2=J-1
                     CASE (-2); JJ2=J+1
                     CASE ( 3); KK2=K-1
                     CASE (-3); KK2=K+1
                  END SELECT
                  OBST_INDEX = OBST_INDEX_C(CELL_INDEX(II2,JJ2,KK2))
                  OB2 => OBSTRUCTION(OBST_INDEX)
                  IF (OB2%PYRO3D) OB2_FOUND = .TRUE.
                  ! next, check surrounding cells
                  KK2_IF: IF (.NOT.OB2_FOUND) THEN
                     II2 = I
                     JJ2 = J
                     KK2_LOOP: DO KK2=K-1,K+1,2
                        OBST_INDEX = OBST_INDEX_C(CELL_INDEX(II2,JJ2,KK2))
                        OB2 => OBSTRUCTION(OBST_INDEX)
                        IF (OB2%PYRO3D) THEN
                           OB2_FOUND = .TRUE.
                           EXIT KK2_LOOP
                        ENDIF
                     ENDDO KK2_LOOP
                  ENDIF KK2_IF
                  JJ2_IF: IF (.NOT.OB2_FOUND) THEN
                     II2 = I
                     KK2 = K
                     JJ2_LOOP: DO JJ2=J-1,J+1,2
                        OBST_INDEX = OBST_INDEX_C(CELL_INDEX(II2,JJ2,KK2))
                        OB2 => OBSTRUCTION(OBST_INDEX)
                        IF (OB2%PYRO3D) THEN
                           OB2_FOUND = .TRUE.
                           EXIT JJ2_LOOP
                        ENDIF
                     ENDDO JJ2_LOOP
                  ENDIF JJ2_IF
                  II2_IF: IF (.NOT.OB2_FOUND) THEN
                     JJ2 = J
                     KK2 = K
                     II2_LOOP: DO II2=I-1,I+1,2
                        OBST_INDEX = OBST_INDEX_C(CELL_INDEX(II2,JJ2,KK2))
                        OB2 => OBSTRUCTION(OBST_INDEX)
                        IF (OB2%PYRO3D) THEN
                           OB2_FOUND = .TRUE.
                           EXIT II2_LOOP
                        ENDIF
                     ENDDO II2_LOOP
                  ENDIF II2_IF

                  OB2_IF: IF (OB2_FOUND) THEN
                     ! if an accepting cell exists, transfer energy and mass
                     IF (TWO_D) THEN
                        VC2 = DX(II2)*DZ(KK2)
                     ELSE
                        VC2 = DX(II2)*DY(JJ2)*DZ(KK2)
                     ENDIF
                     ! get rho*c for each cell before merge
                     RHO_IN(1:MS%N_MATL) = OB%RHO(I,J,K,1:MS%N_MATL)
                     CALL GET_SOLID_RHOCBAR(RHOCBAR,TMP(I,J,K),OPT_SURF_INDEX=OB%MATL_SURF_INDEX,OPT_RHO_IN=RHO_IN)
                     RHO_IN(1:MS%N_MATL) = OB2%RHO(II2,JJ2,KK2,1:MS%N_MATL)
                     CALL GET_SOLID_RHOCBAR(RHOCBAR2,TMP(II2,JJ2,KK2),OPT_SURF_INDEX=OB2%MATL_SURF_INDEX,OPT_RHO_IN=RHO_IN)
                     ! transfer mass
                     OB2%RHO(II2,JJ2,KK2,1:MS%N_MATL) = OB2%RHO(II2,JJ2,KK2,1:MS%N_MATL) + OB%RHO(I,J,K,1:MS%N_MATL)*VC/VC2
                     ! compute new cell temperature
                     TMP(II2,JJ2,KK2) = (VC*RHOCBAR*TMP(I,J,K)+VC2*RHOCBAR2*TMP(II2,JJ2,KK2))/(VC*RHOCBAR+VC2*RHOCBAR2)
                     TMP(I,J,K) = TMP(IIG,JJG,KKG) ! replace solid cell tmp with nearest gas phase tmp
                     OB%RHO(I,J,K,1:MS%N_MATL) = 0._EB
                  ELSEIF (VSRVC_LOC<SOLID_VOLUME_CLIP_THRESHOLD) THEN OB2_IF
                     ! VS/VC is small, but there are no more cells to accept the mass, clip the mass
                     OB%RHO(I,J,K,1:MS%N_MATL) = 0._EB
                  ENDIF OB2_IF

               ENDIF THRESHOLD_IF
            ENDIF CONSUMABLE_IF

            OB%MASS = SUM(OB%RHO(I,J,K,1:MS%N_MATL))*VC
            IF (OB%MASS<TWO_EPSILON_EB) THEN
               OB%HT3D   = .FALSE.
               OB%PYRO3D = .FALSE.
               Q_DOT_PPP_S(I,J,K) = 0._EB
            ENDIF

         ENDDO I_LOOP_2
      ENDDO
   ENDDO
ENDDO OBST_LOOP_2

END SUBROUTINE SOLID_PYROLYSIS_3D


SUBROUTINE CRANK_TEST_1(DIM)
! Initialize solid temperature profile for simple 1D verification test
! J. Crank, The Mathematics of Diffusion, 2nd Ed., Oxford Press, 1975, Sec 2.3.
INTEGER, INTENT(IN) :: DIM ! DIM=1,2,3 for x,y,z dimensions
INTEGER :: I,J,K,IC
REAL(EB), PARAMETER :: LL=1._EB, AA=100._EB, NN=2._EB, X_0=-.5_EB

DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IC = CELL_INDEX(I,J,K)
         IF (.NOT.SOLID(IC)) CYCLE
         SELECT CASE(DIM)
            CASE(1)
               TMP(I,J,K) = TMPA + AA * SIN(NN*PI*(XC(I)-X_0)/LL) ! TMPA = 293.15 K
            CASE(2)
               TMP(I,J,K) = TMPA + AA * SIN(NN*PI*(YC(J)-X_0)/LL)
            CASE(3)
               TMP(I,J,K) = TMPA + AA * SIN(NN*PI*(ZC(K)-X_0)/LL)
         END SELECT
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE CRANK_TEST_1


END SUBROUTINE THERMAL_BC


SUBROUTINE DIFFUSIVITY_BC

! Calculate the term RHO_D_F=RHO*D at the wall

INTEGER :: IW,ICF

IF (N_TRACKED_SPECIES==1) RETURN

! Loop over all WALL cells

WALL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
   WC=>WALL(IW)
   IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY .OR. &
       WC%BOUNDARY_TYPE==OPEN_BOUNDARY .OR. &
       WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY) CYCLE WALL_LOOP
   ONE_D => WC%ONE_D
   CALL CALCULATE_RHO_D_F
ENDDO WALL_LOOP

! Loop over all cut face cells

CFACE_LOOP: DO ICF=1,N_CFACE_CELLS
   ONE_D => CFACE(ICF)%ONE_D
   CALL CALCULATE_RHO_D_F
ENDDO CFACE_LOOP

CONTAINS

SUBROUTINE CALCULATE_RHO_D_F

INTEGER :: N,ITMP,IIG,JJG,KKG

IF (LES .AND. .NOT. RESEARCH_MODE) THEN ! default
   DO N=1,N_TRACKED_SPECIES
      ONE_D%RHO_D_F(N) = ONE_D%MU_G*RSC*ONE_D%RHO_F/ONE_D%RHO_G
   ENDDO
ELSEIF (LES .AND. RESEARCH_MODE) THEN
   ITMP = MIN(4999,NINT(ONE_D%TMP_F))
   IIG = ONE_D%IIG
   JJG = ONE_D%JJG
   KKG = ONE_D%KKG
   DO N=1,N_TRACKED_SPECIES
      ONE_D%RHO_D_F(N) = ONE_D%RHO_F*( D_Z(ITMP,N) + (ONE_D%MU_G-MU_DNS(IIG,JJG,KKG))/ONE_D%RHO_G*RSC )
   ENDDO
ELSE ! DNS
   ITMP = MIN(4999,NINT(ONE_D%TMP_F))
   DO N=1,N_TRACKED_SPECIES
      ONE_D%RHO_D_F(N) = ONE_D%RHO_F*D_Z(ITMP,N)
   ENDDO
ENDIF

END SUBROUTINE CALCULATE_RHO_D_F

END SUBROUTINE DIFFUSIVITY_BC


SUBROUTINE SPECIES_BC(T,DT,NM)

! Compute the species mass fractions at the boundary, ZZ_F

USE PHYSICAL_FUNCTIONS, ONLY: GET_AVERAGE_SPECIFIC_HEAT,GET_SPECIFIC_HEAT,GET_SENSIBLE_ENTHALPY,SURFACE_DENSITY, &
                              GET_SPECIFIC_GAS_CONSTANT
USE TRAN, ONLY: GET_IJK
USE OUTPUT_DATA, ONLY: M_DOT
REAL(EB) :: RADIUS,AREA_SCALING,RVC,M_DOT_PPP_SINGLE,CP,CPBAR,MW_RATIO,H_G,DELTA_H_G,ZZ_GET(1:N_TRACKED_SPECIES),CPBAR2,DENOM
REAL(EB), INTENT(IN) :: T,DT
INTEGER, INTENT(IN) :: NM
INTEGER :: II,JJ,KK,IIG,JJG,KKG,IW,ICF,NS,IP,I_FUEL,SPECIES_BC_INDEX
REAL(EB), POINTER, DIMENSION(:,:) :: PBAR_P=>NULL()
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZP=>NULL()

IF (VEG_LEVEL_SET) RETURN

IF (PREDICTOR) THEN
   PBAR_P => PBAR_S
   ZZP => ZZS
ELSE
   PBAR_P => PBAR
   ZZP => ZZ
ENDIF

! Add evaporating gases from solid particles to the mesh using a volumetric source term

PARTICLE_LOOP: DO IP=1,NLP

   LP  => LAGRANGIAN_PARTICLE(IP)
   LPC => LAGRANGIAN_PARTICLE_CLASS(LP%CLASS_INDEX)

   IF (.NOT.LPC%SOLID_PARTICLE) CYCLE PARTICLE_LOOP

   SF => SURFACE(LPC%SURF_INDEX)
   ONE_D => LP%ONE_D
   IIG = ONE_D%IIG
   JJG = ONE_D%JJG
   KKG = ONE_D%KKG
   II  = IIG
   JJ  = JJG
   KK  = KKG

   CALL CALCULATE_ZZ_F

   ! Only do basic boundary conditions during the PREDICTOR stage of time step.

   IF (PREDICTOR) CYCLE PARTICLE_LOOP

   ! Get particle radius and surface area

   IF (SF%PYROLYSIS_MODEL==PYROLYSIS_PREDICTED) THEN
      RADIUS = SF%INNER_RADIUS + SUM(ONE_D%LAYER_THICKNESS(1:SF%N_LAYERS))
   ELSE
      RADIUS = SF%INNER_RADIUS + SF%THICKNESS
   ENDIF

   IF (ABS(RADIUS)<TWO_EPSILON_EB) CYCLE PARTICLE_LOOP

   AREA_SCALING = 1._EB
   IF (LPC%DRAG_LAW /= SCREEN_DRAG .AND. LPC%DRAG_LAW /= POROUS_DRAG) THEN
      SELECT CASE(SF%GEOMETRY)
         CASE(SURF_CARTESIAN)
            ONE_D%AREA = 2._EB*SF%LENGTH*SF%WIDTH
         CASE(SURF_CYLINDRICAL)
            ONE_D%AREA  = TWOPI*RADIUS*SF%LENGTH
            IF (SF%THERMAL_BC_INDEX == THERMALLY_THICK) AREA_SCALING = (SF%INNER_RADIUS+SF%THICKNESS)/RADIUS
         CASE(SURF_SPHERICAL)
            ONE_D%AREA  = 4._EB*PI*RADIUS**2
            IF (SF%THERMAL_BC_INDEX == THERMALLY_THICK) AREA_SCALING = ((SF%INNER_RADIUS+SF%THICKNESS)/RADIUS)**2
      END SELECT
   ENDIF

   ! In PYROLYSIS, all the mass fluxes are normalized by a virtual area based on the INITIAL radius.
   ! Here, correct the mass flux using the CURRENT radius.

   IF (CALL_HT_1D) THEN
      ONE_D%MASSFLUX(1:N_TRACKED_SPECIES)      = ONE_D%MASSFLUX(1:N_TRACKED_SPECIES)     *AREA_SCALING
      ONE_D%MASSFLUX_SPEC(1:N_TRACKED_SPECIES) = ONE_D%MASSFLUX_SPEC(1:N_TRACKED_SPECIES)*AREA_SCALING
      ONE_D%MASSFLUX_MATL(1:SF%N_MATL)         = ONE_D%MASSFLUX_MATL(1:SF%N_MATL)        *AREA_SCALING
   ENDIF

   ! Add evaporated particle species to gas phase and compute resulting contribution to the divergence

   RVC = RDX(IIG)*RRN(IIG)*RDY(JJG)*RDZ(KKG)
   ZZ_GET(1:N_TRACKED_SPECIES) = ONE_D%ZZ_G(1:N_TRACKED_SPECIES)

   CALL GET_SPECIFIC_HEAT(ZZ_GET,CP,ONE_D%TMP_G)
   H_G = CP*ONE_D%TMP_G

   DO NS=1,N_TRACKED_SPECIES
      IF (ABS(ONE_D%MASSFLUX(NS))<=TWO_EPSILON_EB) CYCLE
      MW_RATIO = SPECIES_MIXTURE(NS)%RCON/ONE_D%RSUM_G
      M_DOT_PPP_SINGLE = LP%PWT*ONE_D%MASSFLUX(NS)*ONE_D%AREA*RVC
      LP%M_DOT = ONE_D%MASSFLUX(NS)*ONE_D%AREA
      ZZ_GET=0._EB
      IF (NS>0) ZZ_GET(NS)=1._EB
      CALL GET_AVERAGE_SPECIFIC_HEAT(ZZ_GET,CPBAR,ONE_D%TMP_G)
      CALL GET_AVERAGE_SPECIFIC_HEAT(ZZ_GET,CPBAR2,LP%ONE_D%TMP_F)
      DELTA_H_G = CPBAR2*LP%ONE_D%TMP_F-CPBAR*ONE_D%TMP_G
      D_SOURCE(IIG,JJG,KKG) = D_SOURCE(IIG,JJG,KKG) + M_DOT_PPP_SINGLE*(MW_RATIO + DELTA_H_G/H_G)/ONE_D%RHO_G
      M_DOT_PPP(IIG,JJG,KKG,NS) = M_DOT_PPP(IIG,JJG,KKG,NS) + M_DOT_PPP_SINGLE
   ENDDO

   ! Calculate contribution to divergence term due to convective heat transfer from particle

   D_SOURCE(IIG,JJG,KKG) = D_SOURCE(IIG,JJG,KKG) - ONE_D%QCONF*ONE_D%AREA*RVC/(ONE_D%RHO_G*H_G) * LP%PWT

   ! Calculate the mass flux of fuel gas from particles

   I_FUEL = 0
   IF (N_REACTIONS>0) I_FUEL = REACTION(1)%FUEL_SMIX_INDEX

   IF (CORRECTOR) THEN
      IF (I_FUEL>0) &
      M_DOT(2,NM) = M_DOT(2,NM) +     ONE_D%MASSFLUX(I_FUEL)*ONE_D%AREA*LP%PWT
      M_DOT(4,NM) = M_DOT(4,NM) + SUM(ONE_D%MASSFLUX)       *ONE_D%AREA*LP%PWT
   ENDIF

   ! Calculate particle mass

   CALC_LP_MASS:IF (SF%THERMALLY_THICK) THEN
      SELECT CASE (SF%GEOMETRY)
         CASE (SURF_CARTESIAN)
            LP%MASS = 2._EB*SF%LENGTH*SF%WIDTH*SF%THICKNESS*SURFACE_DENSITY(NM,1,LAGRANGIAN_PARTICLE_INDEX=IP)
          CASE (SURF_CYLINDRICAL)
            LP%MASS = SF%LENGTH*PI*(SF%INNER_RADIUS+SF%THICKNESS)**2*SURFACE_DENSITY(NM,1,LAGRANGIAN_PARTICLE_INDEX=IP)
         CASE (SURF_SPHERICAL)
            LP%MASS = FOTHPI*(SF%INNER_RADIUS+SF%THICKNESS)**3*SURFACE_DENSITY(NM,1,LAGRANGIAN_PARTICLE_INDEX=IP)
      END SELECT
   ENDIF CALC_LP_MASS

ENDDO PARTICLE_LOOP

! Loop through the wall cells, apply mass boundary conditions

WALL_CELL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS

   WC => WALL(IW)

   IF (WC%BOUNDARY_TYPE==OPEN_BOUNDARY)         CYCLE WALL_CELL_LOOP
   IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY)         CYCLE WALL_CELL_LOOP
   IF (WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY) CYCLE WALL_CELL_LOOP

   ONE_D => WC%ONE_D
   SF => SURFACE(WC%SURF_INDEX)
   II  = ONE_D%II
   JJ  = ONE_D%JJ
   KK  = ONE_D%KK
   IIG = ONE_D%IIG
   JJG = ONE_D%JJG
   KKG = ONE_D%KKG

   CALL CALCULATE_ZZ_F(WALL_INDEX=IW)

   ! Only set species mass fraction in the ghost cell if it is not solid

   IF (IW<=N_EXTERNAL_WALL_CELLS .AND. .NOT.SOLID(CELL_INDEX(II,JJ,KK)) .AND. .NOT.SOLID(CELL_INDEX(IIG,JJG,KKG))) &
       ZZP(II,JJ,KK,1:N_TRACKED_SPECIES) = 2._EB*ONE_D%ZZ_F(1:N_TRACKED_SPECIES) - ONE_D%ZZ_G(1:N_TRACKED_SPECIES)

ENDDO WALL_CELL_LOOP

! Loop through the cut face cells, apply mass boundary conditions

CFACE_LOOP: DO ICF=1,N_CFACE_CELLS

   CFA => CFACE(ICF)
   ONE_D => CFA%ONE_D
   SF => SURFACE(CFA%SURF_INDEX)
   KK  = 1

   CALL CALCULATE_ZZ_F

ENDDO CFACE_LOOP

CONTAINS

SUBROUTINE CALCULATE_ZZ_F(WALL_INDEX)

USE HVAC_ROUTINES, ONLY : DUCT_MF
USE PHYSICAL_FUNCTIONS, ONLY: GET_SPECIFIC_GAS_CONSTANT, GET_REALIZABLE_MF
USE MATH_FUNCTIONS, ONLY : EVALUATE_RAMP, BOX_MULLER
REAL(EB) :: UN,DD,MFT,TSI,ZZ_GET(1:N_TRACKED_SPECIES),RSUM_F,MPUA_SUM,RHO_F_PREVIOUS,RN1,RN2,TWOMFT
INTEGER :: N,ITER
INTEGER, INTENT(IN), OPTIONAL :: WALL_INDEX
TYPE(OBSTRUCTION_TYPE), POINTER :: OB=>NULL()

! Special cases for N_TRACKED_SPECIES==1

IF (N_TRACKED_SPECIES==1) THEN

   IF (PRESENT(WALL_INDEX)) THEN
      IF ( WC%NODE_INDEX < 0 .AND. .NOT.SF%SPECIES_BC_INDEX==SPECIFIED_MASS_FLUX ) THEN
         ONE_D%ZZ_F(1) = 1._EB
         RETURN
      ENDIF
   ENDIF

   IF ( SF%SPECIES_BC_INDEX==SPECIFIED_MASS_FLUX .AND. ABS(SF%MASS_FLUX(1))<=TWO_EPSILON_EB ) THEN
      ONE_D%ZZ_F(1) = 1._EB
      RETURN
   ENDIF

ENDIF

! Check if suppression by water is to be applied and sum water on surface

IF (PRESENT(WALL_INDEX) .AND. CORRECTOR .AND. SF%E_COEFFICIENT>0._EB .AND. I_WATER>0) THEN
   IF (SPECIES_MIXTURE(I_WATER)%EVAPORATING) THEN
      MPUA_SUM = 0._EB
      DO N=1,N_LAGRANGIAN_CLASSES
         LPC=>LAGRANGIAN_PARTICLE_CLASS(N)
         IF (LPC%Z_INDEX==I_WATER) MPUA_SUM = MPUA_SUM + WC%LP_MPUA(LPC%ARRAY_INDEX)
      ENDDO
      WC%EW = WC%EW + SF%E_COEFFICIENT*MPUA_SUM*DT
   ENDIF
ENDIF

! Get SPECIES_BC_INDEX and adjust for HVAC

SPECIES_BC_INDEX = SF%SPECIES_BC_INDEX

IF (PRESENT(WALL_INDEX)) THEN
   IF (WC%NODE_INDEX > 0) THEN
      IF (-DUCTNODE(WC%NODE_INDEX)%DIR(1)*DUCT_MF(DUCTNODE(WC%NODE_INDEX)%DUCT_INDEX(1))>=0._EB) THEN
         SPECIES_BC_INDEX = SPECIFIED_MASS_FRACTION
      ELSE
         SPECIES_BC_INDEX = SPECIFIED_MASS_FLUX
      ENDIF
   ENDIF
ENDIF

! Apply the different species boundary conditions to non-thermally thick solids

METHOD_OF_MASS_TRANSFER: SELECT CASE(SPECIES_BC_INDEX)

   CASE (INFLOW_OUTFLOW_MASS_FLUX) METHOD_OF_MASS_TRANSFER

      ! OPEN boundary species BC is done in THERMAL_BC under INFLOW_OUTFLOW

   CASE (NO_MASS_FLUX) METHOD_OF_MASS_TRANSFER

      ONE_D%ZZ_F(1:N_TRACKED_SPECIES) = ONE_D%ZZ_G(1:N_TRACKED_SPECIES)

   CASE (SPECIFIED_MASS_FRACTION) METHOD_OF_MASS_TRANSFER

      IF (ABS(ONE_D%T_IGN-T_BEGIN)<SPACING(ONE_D%T_IGN) .AND. ANY(SF%RAMP_INDEX>=1)) THEN
         IF (PREDICTOR) TSI = T + DT
         IF (CORRECTOR) TSI = T
      ELSE
         IF (PREDICTOR) TSI = T + DT - ONE_D%T_IGN
         IF (CORRECTOR) TSI = T      - ONE_D%T_IGN
      ENDIF

      IF (ONE_D%UWS<=0._EB) THEN
         DO N=2,N_TRACKED_SPECIES
            ZZ_GET(N) = SPECIES_MIXTURE(N)%ZZ0 + EVALUATE_RAMP(TSI,SF%TAU(N),SF%RAMP_INDEX(N))* &
                           (SF%MASS_FRACTION(N)-SPECIES_MIXTURE(N)%ZZ0)
         ENDDO
         ZZ_GET(1) = 1._EB-SUM(ZZ_GET(2:N_TRACKED_SPECIES))
         CALL GET_REALIZABLE_MF(ZZ_GET)
         ONE_D%ZZ_F = ZZ_GET
      ELSE
         ONE_D%ZZ_F(1:N_TRACKED_SPECIES) = ONE_D%ZZ_G(1:N_TRACKED_SPECIES)
      ENDIF

   CASE (SPECIFIED_MASS_FLUX) METHOD_OF_MASS_TRANSFER

      ! If the current time is before the "activation" time, T_IGN, apply simple BCs and get out

      IF (T < ONE_D%T_IGN .OR. INITIALIZATION_PHASE) THEN
         ONE_D%ZZ_F(1:N_TRACKED_SPECIES) = ONE_D%ZZ_G(1:N_TRACKED_SPECIES)
         IF (PREDICTOR) ONE_D%UWS = 0._EB
         IF (CORRECTOR) ONE_D%UW  = 0._EB
         ONE_D%MASSFLUX(1:N_TRACKED_SPECIES) = 0._EB
         ONE_D%MASSFLUX_SPEC(1:N_TRACKED_SPECIES) = 0._EB
         ONE_D%MASSFLUX_MATL(1:SF%N_MATL) = 0._EB
         RETURN
      ENDIF

      ! Zero out the running counter of Mass Flux Total (MFT)

      MFT = 0._EB

      ! If the user has specified the burning rate, evaluate the ramp and other related parameters

      SUM_MASSFLUX_LOOP: DO N=1,N_TRACKED_SPECIES
         IF (ABS(SF%MASS_FLUX(N)) > TWO_EPSILON_EB) THEN  ! Use user-specified ramp-up of mass flux
            IF (ABS(ONE_D%T_IGN-T_BEGIN) < SPACING(ONE_D%T_IGN) .AND. SF%RAMP_INDEX(N)>=1) THEN
               IF (PREDICTOR) TSI = T + DT
               IF (CORRECTOR) TSI = T
            ELSE
               IF (PREDICTOR) TSI = T + DT - ONE_D%T_IGN
               IF (CORRECTOR) TSI = T      - ONE_D%T_IGN
            ENDIF
            ONE_D%MASSFLUX(N) = EVALUATE_RAMP(TSI,SF%TAU(N),SF%RAMP_INDEX(N))*SF%MASS_FLUX(N)
            ONE_D%MASSFLUX_SPEC(N) = ONE_D%MASSFLUX(N)
            ONE_D%MASSFLUX(N) = SF%ADJUST_BURN_RATE(N)*ONE_D%MASSFLUX(N)*ONE_D%AREA_ADJUST
         ENDIF
         MFT = MFT + ONE_D%MASSFLUX(N)
      ENDDO SUM_MASSFLUX_LOOP

      ! Apply user-specified mass flux variation

      IF (SF%MASS_FLUX_VAR > TWO_EPSILON_EB) THEN
         ! generate pairs of standard Gaussian random variables
         CALL BOX_MULLER(RN1,RN2)
         TWOMFT = 2._EB*MFT
         MFT = MFT*(1._EB + RN1*SF%MASS_FLUX_VAR)
         MFT = MAX(0._EB,MIN(TWOMFT,MFT))
      ENDIF

      ! Apply water suppression and mass consumption at WALL cell

      IF (PRESENT(WALL_INDEX)) THEN

         ! Apply water suppression coefficient (EW) at a WALL cell

         IF (WC%EW>TWO_EPSILON_EB) THEN
            ONE_D%MASSFLUX(1:N_TRACKED_SPECIES) = ONE_D%MASSFLUX(1:N_TRACKED_SPECIES)*EXP(-WC%EW)
            ONE_D%MASSFLUX_SPEC(1:N_TRACKED_SPECIES) = ONE_D%MASSFLUX_SPEC(1:N_TRACKED_SPECIES)*EXP(-WC%EW)
         ENDIF

         ! Add total consumed mass to various summing arrays

         CONSUME_MASS: IF (CORRECTOR .AND. SF%THERMALLY_THICK .AND. .NOT.SF%THERMALLY_THICK_HT3D) THEN
            OB => OBSTRUCTION(WC%OBST_INDEX)
            DO N=1,N_TRACKED_SPECIES
               OB%MASS = OB%MASS - ONE_D%MASSFLUX_SPEC(N)*DT*ONE_D%AREA
            ENDDO
         ENDIF CONSUME_MASS

      ENDIF

      ! Compute the cell face value of the species mass fraction to get the right mass flux

      RHO_F_PREVIOUS = ONE_D%RHO_F
      IF (N_TRACKED_SPECIES==1) THEN
         ONE_D%RHO_F = PBAR_P(KK,ONE_D%PRESSURE_ZONE)/(RSUM0*ONE_D%TMP_F)
         ONE_D%ZZ_F(1) = 1._EB
         UN = MFT/ONE_D%RHO_F
      ELSE
         DO ITER=1,3
            UN = MFT/ONE_D%RHO_F
            SPECIES_LOOP: DO N=1,N_TRACKED_SPECIES
               ONE_D%RHO_D_F(N) = ONE_D%RHO_D_F(N)*ONE_D%RHO_F/RHO_F_PREVIOUS
               DD = 2._EB*ONE_D%RHO_D_F(N)*ONE_D%RDN
               DENOM = DD + UN*ONE_D%RHO_F
               IF ( ABS(DENOM) > TWO_EPSILON_EB ) THEN
                  ONE_D%ZZ_F(N) = ( ONE_D%MASSFLUX(N) + DD*ONE_D%ZZ_G(N) ) / DENOM
               ELSE
                  ONE_D%ZZ_F(N) = ONE_D%ZZ_G(N)
               ENDIF
            ENDDO SPECIES_LOOP
            CALL GET_REALIZABLE_MF(ONE_D%ZZ_F)
            ZZ_GET(1:N_TRACKED_SPECIES) = ONE_D%ZZ_F(1:N_TRACKED_SPECIES)
            CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,RSUM_F)
            RHO_F_PREVIOUS = ONE_D%RHO_F
            ONE_D%RHO_F = PBAR_P(KK,ONE_D%PRESSURE_ZONE)/(RSUM_F*ONE_D%TMP_F)
         ENDDO
      ENDIF
      IF (PREDICTOR) ONE_D%UWS = -UN
      IF (CORRECTOR) ONE_D%UW  = -UN

END SELECT METHOD_OF_MASS_TRANSFER

END SUBROUTINE CALCULATE_ZZ_F

END SUBROUTINE SPECIES_BC


SUBROUTINE DENSITY_BC

! Compute density at wall from wall temperatures and mass fractions

USE PHYSICAL_FUNCTIONS, ONLY : GET_SPECIFIC_GAS_CONSTANT
REAL(EB) :: ZZ_GET(1:N_TRACKED_SPECIES),RSUM_F,UWP
INTEGER  :: IW,BOUNDARY_TYPE,KK,ICF
REAL(EB), POINTER, DIMENSION(:,:) :: PBAR_P=>NULL()
REAL(EB), POINTER, DIMENSION(:,:,:) :: RHOP=>NULL()

IF (VEG_LEVEL_SET) RETURN

IF (PREDICTOR) THEN
   PBAR_P => PBAR_S
   RHOP => RHOS
ELSE
   PBAR_P => PBAR
   RHOP => RHO
ENDIF

! Loop over all wall cells

WALL_CELL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
   WC => WALL(IW)
   BOUNDARY_TYPE = WC%BOUNDARY_TYPE
   IF (BOUNDARY_TYPE==NULL_BOUNDARY .OR. BOUNDARY_TYPE==INTERPOLATED_BOUNDARY) CYCLE WALL_CELL_LOOP
   ONE_D => WC%ONE_D
   KK = ONE_D%KK
   CALL CALCULATE_RHO_F
   IF (IW<=N_EXTERNAL_WALL_CELLS .AND. BOUNDARY_TYPE/=OPEN_BOUNDARY) &
      RHOP(WC%ONE_D%II,WC%ONE_D%JJ,WC%ONE_D%KK) = 2._EB*ONE_D%RHO_F - ONE_D%RHO_G
ENDDO WALL_CELL_LOOP

! Loop over all cut face cells

CFACE_LOOP: DO ICF=1,N_CFACE_CELLS
   CFA => CFACE(ICF)
   BOUNDARY_TYPE = CFA%BOUNDARY_TYPE
   ONE_D => CFA%ONE_D
   KK = 1
   CALL CALCULATE_RHO_F
ENDDO CFACE_LOOP

CONTAINS

SUBROUTINE CALCULATE_RHO_F

! Compute density, RHO_F, at non-iterpolated boundaries

ZZ_GET(1:N_TRACKED_SPECIES) = MAX(0._EB,ONE_D%ZZ_F(1:N_TRACKED_SPECIES))
CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,RSUM_F)
ONE_D%RHO_F = PBAR_P(KK,ONE_D%PRESSURE_ZONE)/(RSUM_F*ONE_D%TMP_F)

! If the boundary is solid and gas is being drawn in, set surface variables to equal the adjacent gas phase variables

IF (BOUNDARY_TYPE==SOLID_BOUNDARY) THEN
   IF (PREDICTOR) THEN
      UWP = ONE_D%UWS
   ELSE
      UWP = ONE_D%UW
   ENDIF
   IF (UWP>0._EB) THEN
      ONE_D%ZZ_F(1:N_TRACKED_SPECIES) = ONE_D%ZZ_G(1:N_TRACKED_SPECIES)
      ONE_D%RHO_F = ONE_D%RHO_G
   ENDIF
ENDIF

END SUBROUTINE CALCULATE_RHO_F

END SUBROUTINE DENSITY_BC


SUBROUTINE HVAC_BC

! Compute density at wall from wall temperatures and mass fractions

USE HVAC_ROUTINES, ONLY : NODE_AREA_EX,NODE_TMP_EX,DUCT_MF,NODE_ZZ_EX
USE PHYSICAL_FUNCTIONS, ONLY : GET_SPECIFIC_GAS_CONSTANT,GET_AVERAGE_SPECIFIC_HEAT
REAL(EB) :: ZZ_GET(1:N_TRACKED_SPECIES),UN,MFT,RSUM_F,CP_D,CP_G
INTEGER  :: IW,KK,SURF_INDEX,COUNTER,DU
REAL(EB), POINTER, DIMENSION(:,:) :: PBAR_P=>NULL()

IF (PREDICTOR) THEN
   PBAR_P => PBAR_S
ELSE
   PBAR_P => PBAR
ENDIF

! Loop over all internal and external wall cells

WALL_CELL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS

   WC=>WALL(IW)
   IF (WC%NODE_INDEX == 0) CYCLE WALL_CELL_LOOP
   ONE_D => WC%ONE_D
   SURF_INDEX = WC%SURF_INDEX
   SF => SURFACE(SURF_INDEX)
   KK  = ONE_D%KK
   COUNTER = 0

   ! Compute R*Sum(Y_i/W_i) at the wall

   DU=DUCTNODE(WC%NODE_INDEX)%DUCT_INDEX(1)
   MFT = -DUCTNODE(WC%NODE_INDEX)%DIR(1)*DUCT_MF(DU)/NODE_AREA_EX(WC%NODE_INDEX)
   IF (.NOT. ANY(SF%LEAK_PATH>0)) THEN
      IF (DUCTNODE(WC%NODE_INDEX)%DIR(1)*DUCT_MF(DU) > 0._EB) THEN
         IF (SF%THERMAL_BC_INDEX==HVAC_BOUNDARY) THEN
            ONE_D%TMP_F = NODE_TMP_EX(WC%NODE_INDEX)
            ONE_D%HEAT_TRANS_COEF = 0._EB
            ONE_D%QCONF = 0._EB
         ELSE
            IF (DUCT(DU)%LEAK_ENTHALPY) THEN
               ZZ_GET(1:N_TRACKED_SPECIES) = NODE_ZZ_EX(WC%NODE_INDEX,1:N_TRACKED_SPECIES)
               CALL GET_AVERAGE_SPECIFIC_HEAT(ZZ_GET,CP_G,ONE_D%TMP_G)
               CALL GET_AVERAGE_SPECIFIC_HEAT(ZZ_GET,CP_D,NODE_TMP_EX(WC%NODE_INDEX))
               WC%Q_LEAK = -MFT*(CP_D*NODE_TMP_EX(WC%NODE_INDEX)-CP_G*ONE_D%TMP_G)*ONE_D%RDN
            ENDIF
         ENDIF
      ELSE
         IF (SF%THERMAL_BC_INDEX==HVAC_BOUNDARY) THEN
            ONE_D%TMP_F = ONE_D%TMP_G
            ONE_D%HEAT_TRANS_COEF = 0._EB
            ONE_D%QCONF = 0._EB
         ENDIF
      ENDIF
   ENDIF

   IF (MFT >= 0._EB) THEN
      ZZ_GET(1:N_TRACKED_SPECIES) = ONE_D%ZZ_G(1:N_TRACKED_SPECIES)
      CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,RSUM_F)
      ONE_D%RHO_F = PBAR_P(KK,ONE_D%PRESSURE_ZONE)/(RSUM_F*ONE_D%TMP_G)
      UN = MFT/ONE_D%RHO_F
      IF (PREDICTOR) ONE_D%UWS = UN
      IF (CORRECTOR) ONE_D%UW  = UN
   ELSE
      ONE_D%MASSFLUX(1:N_TRACKED_SPECIES) = -NODE_ZZ_EX(WC%NODE_INDEX,1:N_TRACKED_SPECIES)*MFT
   ENDIF

ENDDO WALL_CELL_LOOP

END SUBROUTINE HVAC_BC


SUBROUTINE ZETA_BC

INTEGER :: IW,II,JJ,KK,IIG,JJG,KKG,IOR
REAL(EB) :: UN
LOGICAL :: INFLOW
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZP=>NULL()

IF (PREDICTOR) THEN
   ZZP=>ZZS
ELSE
   ZZP=>ZZ
ENDIF

! Loop through the wall cells, apply zeta boundary conditions

WALL_CELL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS

   WC => WALL(IW)
   ONE_D => WC%ONE_D
   SF => SURFACE(WC%SURF_INDEX)
   II  = WC%ONE_D%II
   JJ  = WC%ONE_D%JJ
   KK  = WC%ONE_D%KK
   IIG = WC%ONE_D%IIG
   JJG = WC%ONE_D%JJG
   KKG = WC%ONE_D%KKG
   IOR = WC%ONE_D%IOR

   BOUNDARY_TYPE_SELECT: SELECT CASE(WC%BOUNDARY_TYPE)

      CASE DEFAULT

         ONE_D%ZZ_F(ZETA_INDEX) = INITIAL_UNMIXED_FRACTION

      CASE(SOLID_BOUNDARY)

         IF (ONE_D%UWS<=0._EB) THEN
            ONE_D%ZZ_F(ZETA_INDEX) = SF%ZETA_FRONT
         ELSE
            ONE_D%ZZ_F(ZETA_INDEX) = ONE_D%ZZ_G(ZETA_INDEX)  ! Use nearest gas cell value if the wall is sucking in gas
         ENDIF

      CASE(OPEN_BOUNDARY)

         INFLOW = .FALSE.
         SELECT CASE(IOR)
            CASE( 1); UN = U(II,JJ,KK)
            CASE(-1); UN = -U(II-1,JJ,KK)
            CASE( 2); UN = V(II,JJ,KK)
            CASE(-2); UN = -V(II,JJ-1,KK)
            CASE( 3); UN = W(II,JJ,KK)
            CASE(-3); UN = -W(II,JJ,KK-1)
         END SELECT
         IF (UN>TWO_EPSILON_EB) INFLOW = .TRUE.

         IF (INFLOW) THEN
            ONE_D%ZZ_F(ZETA_INDEX)=INITIAL_UNMIXED_FRACTION
         ELSE
            ONE_D%ZZ_F(ZETA_INDEX)=ONE_D%ZZ_G(ZETA_INDEX)
         ENDIF
         ZZP(II,JJ,KK,ZETA_INDEX) = ONE_D%ZZ_F(ZETA_INDEX) ! Ghost cell values

      CASE(INTERPOLATED_BOUNDARY)

         ! handled in SPECIES_BC

   END SELECT BOUNDARY_TYPE_SELECT

   ! Only set ghost cell if it is not solid

   IF (IW<=N_EXTERNAL_WALL_CELLS .AND. .NOT.SOLID(CELL_INDEX(II,JJ,KK)) .AND. .NOT.SOLID(CELL_INDEX(IIG,JJG,KKG))) &
       ZZP(II,JJ,KK,ZETA_INDEX) = ONE_D%ZZ_F(ZETA_INDEX)

ENDDO WALL_CELL_LOOP

END SUBROUTINE ZETA_BC


SUBROUTINE SOLID_HEAT_TRANSFER_1D(NM,T,DT_BC,PARTICLE_INDEX,WALL_INDEX,CFACE_INDEX)

! Loop through all the boundary cells that require a 1-D heat transfer calc

USE GEOMETRY_FUNCTIONS
USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
USE COMP_FUNCTIONS, ONLY: SHUTDOWN
REAL(EB), INTENT(IN) :: DT_BC,T
INTEGER, INTENT(IN) :: NM
INTEGER, INTENT(IN), OPTIONAL:: WALL_INDEX,PARTICLE_INDEX,CFACE_INDEX
REAL(EB) :: DUMMY,DTMP,QDXKF,QDXKB,RR,RFACF,RFACB,RFACF2,RFACB2, &
            DXKF,DXKB,QRADINB,RFLUX_UP,RFLUX_DOWN,E_WALLB, &
            VOLSUM, REGRID_MAX, REGRID_SUM,  &
            DXF, DXB,HTCB,Q_WATER_F,Q_WATER_B,TMP_F_OLD, RHO_S0,DT2_BC,TOLERANCE,LAYER_DIVIDE,TMP_BACK,&
            M_DOT_G_PPP_ADJUST(N_TRACKED_SPECIES),M_DOT_G_PPP_ACTUAL(N_TRACKED_SPECIES),&
            M_DOT_S_PPP(MAX_MATERIALS),GEOM_FACTOR,RHO_TEMP(MAX_MATERIALS)
INTEGER :: IIB,JJB,KKB,IWB,NWP,KK,I,NR,NL,N,I_OBST,NS,N_LAYER_CELLS_NEW(MAX_LAYERS),N_CELLS
REAL(EB) :: SMALLEST_CELL_SIZE(MAX_LAYERS),THICKNESS
REAL(EB),ALLOCATABLE,DIMENSION(:) :: TMP_W_NEW
REAL(EB),ALLOCATABLE,DIMENSION(:,:) :: INT_WGT
INTEGER  :: NWP_NEW,I_GRAD,STEPCOUNT,IZERO,SURF_INDEX,PART_INDEX=0
LOGICAL :: REMESH,ITERATE,E_FOUND
CHARACTER(MESSAGE_LENGTH) :: MESSAGE
TYPE(WALL_TYPE), POINTER :: WALL_P

! Copy commonly used derived type variables into local variables.

UNPACK_WALL_PARTICLE: IF (PRESENT(WALL_INDEX)) THEN

   WC => WALL(WALL_INDEX)
   SURF_INDEX = WC%SURF_INDEX
   SF => SURFACE(SURF_INDEX)
   ONE_D => WC%ONE_D
   KK  = ONE_D%KK
   I_OBST = WC%OBST_INDEX
   IWB = WC%BACK_INDEX

   ! Take away energy flux due to water evaporation

   IF (NLP>0) THEN
      Q_WATER_F  = -SUM(WC%LP_CPUA(:))
   ELSE
      Q_WATER_F  = 0._EB
   ENDIF

ELSEIF (PRESENT(CFACE_INDEX)) THEN UNPACK_WALL_PARTICLE

   CFA => CFACE(CFACE_INDEX)
   SURF_INDEX = CFA%SURF_INDEX
   SF => SURFACE(SURF_INDEX)
   ONE_D => CFA%ONE_D
   KK  = ONE_D%KK
   I_OBST = 0
   IWB    = -1
   Q_WATER_F = 0._EB

ELSEIF (PRESENT(PARTICLE_INDEX)) THEN UNPACK_WALL_PARTICLE

   LP => LAGRANGIAN_PARTICLE(PARTICLE_INDEX)
   ONE_D => LP%ONE_D
   SURF_INDEX = LAGRANGIAN_PARTICLE_CLASS(LP%CLASS_INDEX)%SURF_INDEX
   PART_INDEX = PARTICLE_INDEX
   SF => SURFACE(SURF_INDEX)
   KK  = ONE_D%KKG
   I_OBST = 0
   IWB = -1
   Q_WATER_F  = 0._EB

ENDIF UNPACK_WALL_PARTICLE

! Save the old value of the surface temperature to be used at the end of the routine

TMP_F_OLD = ONE_D%TMP_F

! If the fuel has burned away, return

IF (ONE_D%BURNAWAY) THEN
   ONE_D%MASSFLUX(1:N_TRACKED_SPECIES) = 0._EB
   RETURN
ENDIF

! Special case where the gas temperature is fixed by the user

IF (ASSUMED_GAS_TEMPERATURE > 0._EB) ONE_D%TMP_G = TMPA + EVALUATE_RAMP(T-T_BEGIN,DUMMY,I_RAMP_AGT)*(ASSUMED_GAS_TEMPERATURE-TMPA)

! Compute convective heat flux at the surface

DTMP = ONE_D%TMP_G - ONE_D%TMP_F
IF (PRESENT(WALL_INDEX)) THEN
   ONE_D%HEAT_TRANS_COEF = HEAT_TRANSFER_COEFFICIENT(DTMP,SF%H_FIXED,SURF_INDEX,WALL_INDEX=WALL_INDEX)
ELSEIF (PRESENT(CFACE_INDEX)) THEN
   ONE_D%HEAT_TRANS_COEF = HEAT_TRANSFER_COEFFICIENT(DTMP,SF%H_FIXED,SURF_INDEX,CFACE_INDEX=CFACE_INDEX)
ELSEIF (PRESENT(PARTICLE_INDEX)) THEN
   ONE_D%HEAT_TRANS_COEF = HEAT_TRANSFER_COEFFICIENT(DTMP,SF%H_FIXED,SURF_INDEX,PARTICLE_INDEX=PARTICLE_INDEX)
ENDIF
ONE_D%QCONF = ONE_D%HEAT_TRANS_COEF*DTMP

! Exponents for cylindrical or spherical coordinates

SELECT CASE(SF%GEOMETRY)
   CASE(SURF_CARTESIAN)   ; I_GRAD = 1
   CASE(SURF_CYLINDRICAL) ; I_GRAD = 2
   CASE(SURF_SPHERICAL)   ; I_GRAD = 3
END SELECT

! Compute back side emissivity

E_WALLB = SF%EMISSIVITY_BACK
IF (E_WALLB < 0._EB .AND. SF%BACKING /= INSULATED) THEN
   E_WALLB = 0._EB
   VOLSUM = 0._EB
   IF (SF%PYROLYSIS_MODEL==PYROLYSIS_PREDICTED) THEN
      NWP = SUM(ONE_D%N_LAYER_CELLS(1:SF%N_LAYERS))
   ELSE
      NWP = SF%N_CELLS_INI
   ENDIF
   DO N=1,SF%N_MATL
      IF (ONE_D%RHO(NWP,N)<=TWO_EPSILON_EB) CYCLE
      ML => MATERIAL(SF%MATL_INDEX(N))
      VOLSUM = VOLSUM + ONE_D%RHO(NWP,N)/ML%RHO_S
      E_WALLB = E_WALLB + ONE_D%RHO(NWP,N)*ML%EMISSIVITY/ML%RHO_S
   ENDDO
   IF (VOLSUM > 0._EB) E_WALLB = E_WALLB/VOLSUM
ENDIF

! Get heat losses from convection and radiation out of back of surface

LAYER_DIVIDE = SF%LAYER_DIVIDE

SELECT CASE(SF%BACKING)

   CASE(VOID)  ! Non-insulated backing to an ambient void

      IF (SF%TMP_BACK>0._EB) THEN
         TMP_BACK = SF%TMP_BACK
      ELSE
         TMP_BACK = TMP_0(KK)
      ENDIF
      DTMP = TMP_BACK - ONE_D%TMP_B
      HTCB = HEAT_TRANSFER_COEFFICIENT(DTMP,SF%H_FIXED_B,SURF_INDEX=SURF_INDEX)
      QRADINB   =  E_WALLB*SIGMA*TMP_BACK**4
      Q_WATER_B = 0._EB

   CASE(INSULATED)  ! No heat transfer out the back

      HTCB      = 0._EB
      QRADINB   = 0._EB
      E_WALLB   = 0._EB
      Q_WATER_B = 0._EB
      IF (SF%TMP_BACK>0._EB) THEN
         TMP_BACK = SF%TMP_BACK
      ELSE
         TMP_BACK = TMP_0(KK)
      ENDIF

   CASE(EXPOSED)  ! The backside is exposed to gas in current or adjacent mesh.

      Q_WATER_B = 0._EB

      SELECT CASE(IWB)

         CASE(1:) ! Solid backing

            IF (WC%BACK_MESH/=NM .AND. WC%BACK_MESH>0) THEN  ! Back side is in other mesh.
               TMP_BACK = OMESH(WC%BACK_MESH)%EXPOSED_WALL(IWB)%TMP_GAS
               DTMP = TMP_BACK - ONE_D%TMP_B
               HTCB = HEAT_TRANSFER_COEFFICIENT(DTMP,SF%H_FIXED_B,SURF_INDEX=SURF_INDEX)
               QRADINB  = OMESH(WC%BACK_MESH)%EXPOSED_WALL(IWB)%QRADIN
            ELSE  ! Back side is in current mesh.
               WALL_P => WALL(IWB)
               IIB = WALL_P%ONE_D%IIG
               JJB = WALL_P%ONE_D%JJG
               KKB = WALL_P%ONE_D%KKG
               TMP_BACK  = TMP(IIB,JJB,KKB)
               DTMP = TMP_BACK - ONE_D%TMP_B
               HTCB = HEAT_TRANSFER_COEFFICIENT(DTMP,SF%H_FIXED_B,SURF_INDEX=WALL_P%SURF_INDEX,WALL_INDEX=IWB)
               WALL_P%ONE_D%HEAT_TRANS_COEF = HTCB
               QRADINB  = WALL_P%ONE_D%QRADIN
               IF (NLP>0) Q_WATER_B = -SUM(WALL_P%LP_CPUA(:))
            ENDIF

         CASE DEFAULT  ! The back side is an ambient void.

            TMP_BACK = TMP_0(KK)
            DTMP = TMP_BACK - ONE_D%TMP_B
            HTCB = HEAT_TRANSFER_COEFFICIENT(DTMP,SF%H_FIXED_B,SURF_INDEX=SURF_INDEX)
            QRADINB  =  E_WALLB*SIGMA*TMP_BACK**4
            LAYER_DIVIDE = REAL(SF%N_LAYERS+1)

         CASE(-1) ! Particle "backside conditions" are assumed to be from the same gas cell

            TMP_BACK  = ONE_D%TMP_G
            DTMP = TMP_BACK - ONE_D%TMP_B
            HTCB = HEAT_TRANSFER_COEFFICIENT(DTMP,SF%H_FIXED_B,SURF_INDEX,PARTICLE_INDEX=PARTICLE_INDEX)
            QRADINB  = ONE_D%QRADIN

      END SELECT

END SELECT

! Compute grid for reacting nodes

COMPUTE_GRID: IF (SF%PYROLYSIS_MODEL==PYROLYSIS_PREDICTED) THEN
   NWP = SUM(ONE_D%N_LAYER_CELLS(1:SF%N_LAYERS))
   CALL GET_WALL_NODE_WEIGHTS(NWP,SF%N_LAYERS,ONE_D%N_LAYER_CELLS(1:SF%N_LAYERS),ONE_D%LAYER_THICKNESS,SF%GEOMETRY, &
      ONE_D%X(0:NWP),LAYER_DIVIDE,DX_S(1:NWP),RDX_S(0:NWP+1),RDXN_S(0:NWP),DX_WGT_S(0:NWP),DXF,DXB,&
      LAYER_INDEX(0:NWP+1),MF_FRAC(1:NWP),SF%INNER_RADIUS)
ELSE COMPUTE_GRID
   NWP                  = SF%N_CELLS_INI
   DXF                  = SF%DXF
   DXB                  = SF%DXB
   DX_S(1:NWP)          = SF%DX(1:NWP)
   RDX_S(0:NWP+1)       = SF%RDX(0:NWP+1)
   RDXN_S(0:NWP)        = SF%RDXN(0:NWP)
   DX_WGT_S(0:NWP)      = SF%DX_WGT(0:NWP)
   LAYER_INDEX(0:NWP+1) = SF%LAYER_INDEX(0:NWP+1)
   MF_FRAC(1:NWP)       = SF%MF_FRAC(1:NWP)
ENDIF COMPUTE_GRID

! Get total thickness of solid and compute radius for cylindrical and spherical coordinate systems.

THICKNESS = SUM(ONE_D%LAYER_THICKNESS(1:SF%N_LAYERS))

DO I=0,NWP
   R_S(I) = SF%INNER_RADIUS + ONE_D%X(NWP) - ONE_D%X(I)
ENDDO

! Calculate reaction rates based on the solid phase reactions

REMESH = .FALSE.
Q_S = 0._EB
IF (STORE_Q_DOT_PPP_S) Q_DOT_PPP_S = 0._EB

PYROLYSIS_PREDICTED_IF: IF (SF%PYROLYSIS_MODEL==PYROLYSIS_PREDICTED) THEN

   ! Set mass fluxes to 0 and CHANGE_THICKNESS to false.

   ONE_D%MASSFLUX(1:N_TRACKED_SPECIES)      = 0._EB
   ONE_D%MASSFLUX_SPEC(1:N_TRACKED_SPECIES) = 0._EB
   ONE_D%MASSFLUX_MATL(1:SF%N_MATL)         = 0._EB
   ONE_D%CHANGE_THICKNESS                   = .FALSE.

   POINT_LOOP1: DO I=1,NWP

      RHO_S0 = SF%LAYER_DENSITY(LAYER_INDEX(I))
      REGRID_FACTOR(I) = 1._EB
      REGRID_MAX       = 0._EB
      REGRID_SUM       = 0._EB

      ! Compute the pyrolysis reaction terms

      IF (MATERIAL(SF%MATL_INDEX(1))%PYROLYSIS_MODEL/=PYROLYSIS_LIQUID .OR. I<2) THEN

         ! Send the array of component densities, ONE_D%RHO(I,N), into the PYROLYSIS routine

         RHO_TEMP(1:SF%N_MATL) = ONE_D%RHO(I,1:SF%N_MATL)
         CALL PYROLYSIS(SF%N_MATL,SF%MATL_INDEX,SURF_INDEX,ONE_D%IIG,ONE_D%JJG,ONE_D%KKG,ONE_D%TMP(I),ONE_D%TMP_F,&
                        RHO_TEMP(1:SF%N_MATL),RHO_S0,ONE_D%X(I-1),DT_BC,&
                        M_DOT_G_PPP_ADJUST,M_DOT_G_PPP_ACTUAL,M_DOT_S_PPP,Q_S(I),PART_INDEX=PART_INDEX)
         ONE_D%RHO(I,1:SF%N_MATL) = RHO_TEMP(1:SF%N_MATL)

         ! Compute the mass flux of reaction gases at the surface

         GEOM_FACTOR = MF_FRAC(I)*(R_S(I-1)**I_GRAD-R_S(I)**I_GRAD)/ &
                               (I_GRAD*(SF%INNER_RADIUS+SF%THICKNESS)**(I_GRAD-1))
         DO NS = 1,N_TRACKED_SPECIES
            ONE_D%MASSFLUX(NS)      = ONE_D%MASSFLUX(NS)      + M_DOT_G_PPP_ADJUST(NS)*GEOM_FACTOR
            ONE_D%MASSFLUX_SPEC(NS) = ONE_D%MASSFLUX_SPEC(NS) + M_DOT_G_PPP_ACTUAL(NS)*GEOM_FACTOR
         ENDDO
         DO N=1,SF%N_MATL
            ONE_D%MASSFLUX_MATL(N)  = ONE_D%MASSFLUX_MATL(N)  + M_DOT_S_PPP(N)*GEOM_FACTOR
         ENDDO

      ENDIF

      ! Compute regrid factors

      MATERIAL_LOOP1a: DO N=1,SF%N_MATL
         ML  => MATERIAL(SF%MATL_INDEX(N))
         IF (ML%PYROLYSIS_MODEL==PYROLYSIS_LIQUID) THEN
            REMESH = .TRUE.
            IF (I>1) THEN
               REGRID_SUM = 1._EB
               CYCLE MATERIAL_LOOP1a
            ENDIF
         ENDIF
         REGRID_MAX = MAX(REGRID_MAX,ONE_D%RHO(I,N)/ML%RHO_S)
         REGRID_SUM = REGRID_SUM + ONE_D%RHO(I,N)/ML%RHO_S
      ENDDO MATERIAL_LOOP1a

      IF (REGRID_SUM <= 1._EB) REGRID_FACTOR(I) = REGRID_SUM
      IF (REGRID_MAX >= ALMOST_ONE) REGRID_FACTOR(I) = REGRID_MAX

      ! If there is any non-shrinking material, the material matrix will remain, and no shrinking is allowed

      MATERIAL_LOOP1b: DO N=1,SF%N_MATL
         IF (ONE_D%RHO(I,N)<=TWO_EPSILON_EB) CYCLE MATERIAL_LOOP1b
         ML  => MATERIAL(SF%MATL_INDEX(N))
         IF (.NOT. ML%ALLOW_SHRINKING) THEN
            REGRID_FACTOR(I) = MAX(REGRID_FACTOR(I),1._EB)
            EXIT MATERIAL_LOOP1b
         ENDIF
      ENDDO MATERIAL_LOOP1b

      ! If there is any non-swelling material, the material matrix will remain, and no swelling is allowed

      MATERIAL_LOOP1c: DO N=1,SF%N_MATL
         IF (ONE_D%RHO(I,N)<=TWO_EPSILON_EB) CYCLE MATERIAL_LOOP1c
         ML  => MATERIAL(SF%MATL_INDEX(N))
         IF (.NOT. ML%ALLOW_SWELLING) THEN
            REGRID_FACTOR(I) = MIN(REGRID_FACTOR(I),1._EB)
            EXIT MATERIAL_LOOP1c
         ENDIF
      ENDDO MATERIAL_LOOP1c

      ! In points that change thickness, update the density

      IF (ABS(REGRID_FACTOR(I)-1._EB)>=TWO_EPSILON_EB) THEN
         ONE_D%CHANGE_THICKNESS=.TRUE.
         MATERIAL_LOOP1d: DO N=1,SF%N_MATL
            IF(REGRID_FACTOR(I)>TWO_EPSILON_EB) ONE_D%RHO(I,N) = ONE_D%RHO(I,N)/REGRID_FACTOR(I)
         ENDDO MATERIAL_LOOP1d
      ENDIF

   ENDDO POINT_LOOP1

   ! Adjust the MASSFLUX of a wall surface cell to account for non-alignment of the mesh.

   IF (PRESENT(WALL_INDEX)) ONE_D%MASSFLUX(1:N_TRACKED_SPECIES) = ONE_D%MASSFLUX(1:N_TRACKED_SPECIES)*ONE_D%AREA_ADJUST

   ! Compute new coordinates if the solid changes thickness. Save new coordinates in X_S_NEW.

   R_S_NEW(NWP) = 0._EB
   DO I=NWP-1,0,-1
      R_S_NEW(I) = ( R_S_NEW(I+1)**I_GRAD + (R_S(I)**I_GRAD-R_S(I+1)**I_GRAD)*REGRID_FACTOR(I+1) )**(1./REAL(I_GRAD))
   ENDDO

   X_S_NEW(0) = 0._EB
   DO I=1,NWP
      X_S_NEW(I) = R_S_NEW(0) - R_S_NEW(I)
      IF ((X_S_NEW(I)-X_S_NEW(I-1)) < TWO_EPSILON_EB) REMESH = .TRUE.
   ENDDO

   ! If the fuel or water massflux is non-zero, set the ignition time

   IF (ONE_D%T_IGN > T) THEN
      IF (SUM(ONE_D%MASSFLUX(1:N_TRACKED_SPECIES)) > 0._EB) ONE_D%T_IGN = T
   ENDIF

   ! Re-generate grid for a wall changing thickness

   N_LAYER_CELLS_NEW = 0
   SMALLEST_CELL_SIZE = 0._EB

   REMESH_GRID: IF (ONE_D%CHANGE_THICKNESS) THEN
      NWP_NEW = 0
      THICKNESS = 0._EB

      I = 0
      LAYER_LOOP: DO NL=1,SF%N_LAYERS

         ONE_D%LAYER_THICKNESS(NL) = X_S_NEW(I+ONE_D%N_LAYER_CELLS(NL)) - X_S_NEW(I)

         ! Remove very thin layers

         IF (ONE_D%LAYER_THICKNESS(NL) < SF%MINIMUM_LAYER_THICKNESS) THEN
            X_S_NEW(I+ONE_D%N_LAYER_CELLS(NL):NWP) = X_S_NEW(I+ONE_D%N_LAYER_CELLS(NL):NWP)-ONE_D%LAYER_THICKNESS(NL)
            ONE_D%LAYER_THICKNESS(NL) = 0._EB
            N_LAYER_CELLS_NEW(NL) = 0
         ELSE
            CALL GET_N_LAYER_CELLS(SF%MIN_DIFFUSIVITY(NL),ONE_D%LAYER_THICKNESS(NL), &
               SF%STRETCH_FACTOR(NL),SF%CELL_SIZE_FACTOR,SF%N_LAYER_CELLS_MAX(NL),N_LAYER_CELLS_NEW(NL),SMALLEST_CELL_SIZE(NL))
            NWP_NEW = NWP_NEW + N_LAYER_CELLS_NEW(NL)
         ENDIF
         IF ( N_LAYER_CELLS_NEW(NL) /= ONE_D%N_LAYER_CELLS(NL)) REMESH = .TRUE.

         THICKNESS = THICKNESS + ONE_D%LAYER_THICKNESS(NL)
         I = I + ONE_D%N_LAYER_CELLS(NL)
      ENDDO LAYER_LOOP

      ! Check that NWP_NEW has not exceeded the allocated space N_CELLS_MAX
      IF (NWP_NEW > SF%N_CELLS_MAX) THEN
         WRITE(MESSAGE,'(A,I5,A,A)') 'ERROR: N_CELLS_MAX should be at least ',NWP_NEW,' for surface ',TRIM(SF%ID)
         CALL SHUTDOWN(MESSAGE)
      ENDIF

      ! Shrinking wall has gone to zero thickness.
      IF (THICKNESS <=TWO_EPSILON_EB) THEN
         ONE_D%TMP(0:NWP+1) = MAX(TMPMIN,TMP_BACK)
         ONE_D%TMP_F        = MIN(TMPMAX,MAX(TMPMIN,TMP_BACK))
         ONE_D%TMP_B        = MIN(TMPMAX,MAX(TMPMIN,TMP_BACK))
         ONE_D%QCONF        = 0._EB
         ONE_D%MASSFLUX(1:N_TRACKED_SPECIES)  = 0._EB
         ONE_D%MASSFLUX_SPEC(1:N_TRACKED_SPECIES) = 0._EB
         ONE_D%MASSFLUX_MATL(1:SF%N_MATL) = 0._EB
         ONE_D%N_LAYER_CELLS(1:SF%N_LAYERS)     = 0
         ONE_D%BURNAWAY          = .TRUE.
         IF (I_OBST > 0) THEN
            IF (OBSTRUCTION(I_OBST)%CONSUMABLE) OBSTRUCTION(I_OBST)%MASS = -1.
         ENDIF
         RETURN
      ENDIF

      ! Set up new node points following shrinking/swelling

      ONE_D%X(0:NWP) = X_S_NEW(0:NWP)

      X_S_NEW = 0._EB
      IF (REMESH) THEN
         CALL GET_WALL_NODE_COORDINATES(NWP_NEW,SF%N_LAYERS,N_LAYER_CELLS_NEW, &
            SMALLEST_CELL_SIZE(1:SF%N_LAYERS),SF%STRETCH_FACTOR(1:SF%N_LAYERS),X_S_NEW(0:NWP_NEW))
         CALL GET_WALL_NODE_WEIGHTS(NWP_NEW,SF%N_LAYERS,N_LAYER_CELLS_NEW,ONE_D%LAYER_THICKNESS,SF%GEOMETRY, &
            X_S_NEW(0:NWP_NEW),LAYER_DIVIDE,DX_S(1:NWP_NEW),RDX_S(0:NWP_NEW+1),RDXN_S(0:NWP_NEW),&
            DX_WGT_S(0:NWP_NEW),DXF,DXB,LAYER_INDEX(0:NWP_NEW+1),MF_FRAC(1:NWP_NEW),SF%INNER_RADIUS)
         ! Interpolate densities and temperature from old grid to new grid
         ALLOCATE(INT_WGT(NWP_NEW,NWP),STAT=IZERO)
         CALL GET_INTERPOLATION_WEIGHTS(SF%N_LAYERS,NWP,NWP_NEW,ONE_D%N_LAYER_CELLS,N_LAYER_CELLS_NEW, &
                                    ONE_D%X(0:NWP),X_S_NEW(0:NWP_NEW),INT_WGT)
         N_CELLS = MAX(NWP,NWP_NEW)
         CALL INTERPOLATE_WALL_ARRAY(N_CELLS,NWP,NWP_NEW,INT_WGT,ONE_D%TMP(1:N_CELLS))
         ONE_D%TMP(0) = 2*ONE_D%TMP_F-ONE_D%TMP(1) !Make sure surface temperature stays the same
         ONE_D%TMP(NWP_NEW+1) = ONE_D%TMP(NWP+1)
         CALL INTERPOLATE_WALL_ARRAY(N_CELLS,NWP,NWP_NEW,INT_WGT,Q_S(1:N_CELLS))
         DO N=1,SF%N_MATL
            ML  => MATERIAL(SF%MATL_INDEX(N))
            CALL INTERPOLATE_WALL_ARRAY(N_CELLS,NWP,NWP_NEW,INT_WGT,ONE_D%RHO(1:N_CELLS,N))
         ENDDO
         DEALLOCATE(INT_WGT)
         ONE_D%N_LAYER_CELLS(1:SF%N_LAYERS) = N_LAYER_CELLS_NEW(1:SF%N_LAYERS)
         NWP = NWP_NEW
         ONE_D%X(0:NWP) = X_S_NEW(0:NWP)      ! Note: X(NWP+1...) are not set to zero.
      ELSE
         CALL GET_WALL_NODE_WEIGHTS(NWP,SF%N_LAYERS,N_LAYER_CELLS_NEW,ONE_D%LAYER_THICKNESS(1:SF%N_LAYERS),SF%GEOMETRY, &
            ONE_D%X(0:NWP),LAYER_DIVIDE,DX_S(1:NWP),RDX_S(0:NWP+1),RDXN_S(0:NWP),DX_WGT_S(0:NWP),DXF,DXB, &
            LAYER_INDEX(0:NWP+1),MF_FRAC(1:NWP),SF%INNER_RADIUS)
      ENDIF

   ENDIF REMESH_GRID

ELSEIF (SF%PYROLYSIS_MODEL==PYROLYSIS_SPECIFIED) THEN PYROLYSIS_PREDICTED_IF

   ! Take off energy corresponding to specified burning rate

   Q_S(1) = Q_S(1) - ONE_D%MASSFLUX(REACTION(1)%FUEL_SMIX_INDEX)*SF%H_V/DX_S(1)

ENDIF PYROLYSIS_PREDICTED_IF

! Calculate thermal properties

K_S     = 0._EB
RHO_S   = 0._EB
RHOCBAR = 0._EB
ONE_D%EMISSIVITY = 0._EB
E_FOUND = .FALSE.

POINT_LOOP3: DO I=1,NWP
   VOLSUM = 0._EB
   MATERIAL_LOOP3: DO N=1,SF%N_MATL
      IF (ONE_D%RHO(I,N)<=TWO_EPSILON_EB) CYCLE MATERIAL_LOOP3
      ML  => MATERIAL(SF%MATL_INDEX(N))
      VOLSUM = VOLSUM + ONE_D%RHO(I,N)/ML%RHO_S
      IF (ML%K_S>0._EB) THEN
         K_S(I) = K_S(I) + ONE_D%RHO(I,N)*ML%K_S/ML%RHO_S
      ELSE
         NR = -NINT(ML%K_S)
         K_S(I) = K_S(I) + ONE_D%RHO(I,N)*EVALUATE_RAMP(ONE_D%TMP(I),0._EB,NR)/ML%RHO_S
      ENDIF

      IF (ML%C_S>0._EB) THEN
         RHOCBAR(I) = RHOCBAR(I) + ONE_D%RHO(I,N)*ML%C_S
      ELSE
         NR = -NINT(ML%C_S)
         RHOCBAR(I) = RHOCBAR(I) + ONE_D%RHO(I,N)*EVALUATE_RAMP(ONE_D%TMP(I),0._EB,NR)
      ENDIF
      IF (.NOT.E_FOUND) ONE_D%EMISSIVITY = ONE_D%EMISSIVITY + ONE_D%RHO(I,N)*ML%EMISSIVITY/ML%RHO_S
      RHO_S(I) = RHO_S(I) + ONE_D%RHO(I,N)

   ENDDO MATERIAL_LOOP3

   IF (VOLSUM > 0._EB) THEN
      K_S(I) = K_S(I)/VOLSUM
      IF (.NOT.E_FOUND) ONE_D%EMISSIVITY = ONE_D%EMISSIVITY/VOLSUM
   ENDIF
   IF (ONE_D%EMISSIVITY>0._EB) E_FOUND = .TRUE.

   IF (K_S(I)<=TWO_EPSILON_EB)      K_S(I)      = 10000._EB
   IF (RHOCBAR(I)<=TWO_EPSILON_EB)  RHOCBAR(I)  = 0.001_EB

ENDDO POINT_LOOP3

! Calculate average K_S between at grid cell boundaries. Store result in K_S

K_S(0)     = K_S(1)
K_S(NWP+1) = K_S(NWP)
DO I=1,NWP-1
   K_S(I)  = 1._EB / ( DX_WGT_S(I)/K_S(I) + (1._EB-DX_WGT_S(I))/K_S(I+1) )
ENDDO

! Add internal heat source specified by user

IF (SF%SPECIFIED_HEAT_SOURCE) THEN
   DO I=1,NWP
      Q_S(I) = Q_S(I)+SF%INTERNAL_HEAT_SOURCE(LAYER_INDEX(I))
   ENDDO
ENDIF

! Calculate internal radiation for Cartesian geometry only

IF (SF%INTERNAL_RADIATION) THEN
   KAPPA_S = 0._EB
   DO I=1,NWP
      VOLSUM = 0._EB
      DO N=1,SF%N_MATL
         IF (ONE_D%RHO(I,N)<=TWO_EPSILON_EB) CYCLE
         ML  => MATERIAL(SF%MATL_INDEX(N))
         VOLSUM = VOLSUM + ONE_D%RHO(I,N)/ML%RHO_S
         KAPPA_S(I) = KAPPA_S(I) + ONE_D%RHO(I,N)*ML%KAPPA_S/ML%RHO_S
      ENDDO
      IF (VOLSUM>0._EB) KAPPA_S(I) = 2._EB*KAPPA_S(I)/(RDX_S(I)*VOLSUM)    ! kappa = 2*dx*kappa
   ENDDO
   ! solution inwards
   RFLUX_UP = ONE_D%QRADIN + (1._EB-ONE_D%EMISSIVITY)*ONE_D%QRADOUT/(ONE_D%EMISSIVITY+1.0E-10_EB)
   DO I=1,NWP
      RFLUX_DOWN =  ( RFLUX_UP + KAPPA_S(I)*SIGMA*ONE_D%TMP(I)**4 ) / (1._EB + KAPPA_S(I))
      Q_S(I) = Q_S(I) + (RFLUX_UP - RFLUX_DOWN)*RDX_S(I)
      RFLUX_UP = RFLUX_DOWN
   ENDDO
   ! solution outwards
   RFLUX_UP = QRADINB + (1._EB-E_WALLB)*RFLUX_UP
   DO I=NWP,1,-1
      RFLUX_DOWN =  ( RFLUX_UP + KAPPA_S(I)*SIGMA*ONE_D%TMP(I)**4 ) / (1._EB + KAPPA_S(I))
      Q_S(I) = Q_S(I) + (RFLUX_UP - RFLUX_DOWN)*RDX_S(I)
      RFLUX_UP = RFLUX_DOWN
   ENDDO
   ONE_D%QRADOUT = ONE_D%EMISSIVITY*RFLUX_DOWN
ENDIF

! Update the 1-D heat transfer equation

DT2_BC = DT_BC
STEPCOUNT = 1
ALLOCATE(TMP_W_NEW(0:NWP+1),STAT=IZERO)
TMP_W_NEW(0:NWP+1) = ONE_D%TMP(0:NWP+1)
WALL_ITERATE: DO
   ITERATE=.FALSE.
   SUB_TIME: DO N=1,STEPCOUNT
      DXKF   = K_S(0)/DXF
      DXKB   = K_S(NWP)/DXB
      DO I=1,NWP
         BBS(I) = -0.5_EB*DT2_BC*K_S(I-1)*RDXN_S(I-1)*RDX_S(I)/RHOCBAR(I) ! DT_BC->DT2_BC
         AAS(I) = -0.5_EB*DT2_BC*K_S(I)  *RDXN_S(I)  *RDX_S(I)/RHOCBAR(I)
      ENDDO
      DDS(1:NWP) = 1._EB - AAS(1:NWP) - BBS(1:NWP)
      DO I=1,NWP
         CCS(I) = TMP_W_NEW(I) - AAS(I)*(TMP_W_NEW(I+1)-TMP_W_NEW(I)) + BBS(I)*(TMP_W_NEW(I)-TMP_W_NEW(I-1)) &
                  + DT2_BC*Q_S(I)/RHOCBAR(I)
      ENDDO
      IF ( (.NOT.RADIATION) .OR. SF%INTERNAL_RADIATION ) THEN
         RFACF = 0.25_EB*ONE_D%HEAT_TRANS_COEF
         RFACB = 0.25_EB*HTCB
      ELSE
         RFACF = 0.25_EB*ONE_D%HEAT_TRANS_COEF + 2._EB*ONE_D%EMISSIVITY*SIGMA*ONE_D%TMP_F**3
         RFACB = 0.25_EB*HTCB + 2._EB*E_WALLB*SIGMA*ONE_D%TMP_B**3
      ENDIF
      RFACF2 = (DXKF-RFACF)/(DXKF+RFACF)
      RFACB2 = (DXKB-RFACB)/(DXKB+RFACB)
      IF ( (.NOT. RADIATION) .OR. SF%INTERNAL_RADIATION ) THEN
         QDXKF = (ONE_D%HEAT_TRANS_COEF*(ONE_D%TMP_G    - 0.5_EB*ONE_D%TMP_F) + Q_WATER_F)/(DXKF+RFACF)
         QDXKB = (HTCB*                 (      TMP_BACK - 0.5_EB*ONE_D%TMP_B) + Q_WATER_B)/(DXKB+RFACB)
      ELSE
         QDXKF = (ONE_D%HEAT_TRANS_COEF*(ONE_D%TMP_G - 0.5_EB*ONE_D%TMP_F) + ONE_D%QRADIN + &
                 3.*ONE_D%EMISSIVITY*SIGMA*ONE_D%TMP_F**4 + Q_WATER_F) /(DXKF+RFACF)
         QDXKB = (HTCB*(TMP_BACK - 0.5_EB*ONE_D%TMP_B) + QRADINB + 3.*E_WALLB*SIGMA*ONE_D%TMP_B**4 + Q_WATER_B) &
               /(DXKB+RFACB)
      ENDIF
      CCS(1)   = CCS(1)   - BBS(1)  *QDXKF
      CCS(NWP) = CCS(NWP) - AAS(NWP)*QDXKB
      DDT(1:NWP) = DDS(1:NWP)
      DDT(1)   = DDT(1)   + BBS(1)  *RFACF2
      DDT(NWP) = DDT(NWP) + AAS(NWP)*RFACB2
      TRIDIAGONAL_SOLVER_1: DO I=2,NWP
         RR     = BBS(I)/DDT(I-1)
         DDT(I) = DDT(I) - RR*AAS(I-1)
         CCS(I) = CCS(I) - RR*CCS(I-1)
      ENDDO TRIDIAGONAL_SOLVER_1
      CCS(NWP)  = CCS(NWP)/DDT(NWP)
      TRIDIAGONAL_SOLVER_2: DO I=NWP-1,1,-1
         CCS(I) = (CCS(I) - AAS(I)*CCS(I+1))/DDT(I)
      ENDDO TRIDIAGONAL_SOLVER_2
      TMP_W_NEW(1:NWP) = MAX(TMPMIN,CCS(1:NWP))
      TMP_W_NEW(0)     = MAX(TMPMIN,TMP_W_NEW(1)  *RFACF2+QDXKF)
      TMP_W_NEW(NWP+1) = MAX(TMPMIN,TMP_W_NEW(NWP)*RFACB2+QDXKB)
      IF (STEPCOUNT==1) THEN
         TOLERANCE = MAXVAL(ABS((TMP_W_NEW-ONE_D%TMP(0:NWP+1))/ONE_D%TMP(0:NWP+1)), &
            ONE_D%TMP(0:NWP+1)>0._EB) ! returns a negative number, if all TMP_S == 0.
         IF (TOLERANCE<0.0_EB) &
         TOLERANCE = MAXVAL(ABS((TMP_W_NEW-ONE_D%TMP(0:NWP+1))/TMP_W_NEW), &
            TMP_W_NEW>0._EB)
         IF (TOLERANCE > 0.2_EB) THEN
            STEPCOUNT = MIN(200,STEPCOUNT * (INT(TOLERANCE/0.2_EB) + 1))
            ITERATE=.TRUE.
            DT2_BC=DT_BC/REAL(STEPCOUNT)
            TMP_W_NEW = ONE_D%TMP(0:NWP+1)
         ENDIF
      ENDIF
      IF (NWP == 1) THEN
         ONE_D%TMP_F = TMP_W_NEW(1)
         ONE_D%TMP_B = ONE_D%TMP_F
      ELSE
         ONE_D%TMP_F  = 0.5_EB*(TMP_W_NEW(0)+TMP_W_NEW(1))
         ONE_D%TMP_B  = 0.5_EB*(TMP_W_NEW(NWP)+TMP_W_NEW(NWP+1))
      ENDIF
      ONE_D%TMP_F  = MIN(TMPMAX,MAX(TMPMIN,ONE_D%TMP_F))
      ONE_D%TMP_B  = MIN(TMPMAX,MAX(TMPMIN,ONE_D%TMP_B))
   ENDDO SUB_TIME
   IF (.NOT. ITERATE) EXIT WALL_ITERATE
ENDDO WALL_ITERATE

ONE_D%TMP(0:NWP+1) = TMP_W_NEW
DEALLOCATE(TMP_W_NEW)

! If the surface temperature exceeds the ignition temperature, burn it

IF (ONE_D%T_IGN > T ) THEN
   IF (ONE_D%TMP_F >= SF%TMP_IGN) ONE_D%T_IGN = T
ENDIF

! Determine convective heat flux at the wall

ONE_D%QCONF = ONE_D%HEAT_TRANS_COEF * (ONE_D%TMP_G - 0.5_EB * (ONE_D%TMP_F + TMP_F_OLD) )

END SUBROUTINE SOLID_HEAT_TRANSFER_1D


SUBROUTINE PYROLYSIS(N_MATS,MATL_INDEX,SURF_INDEX,IIG,JJG,KKG,TMP_S,TMP_F,RHO_S,RHO_S0,DEPTH,DT_BC,&
                     M_DOT_G_PPP_ADJUST,M_DOT_G_PPP_ACTUAL,M_DOT_S_PPP,Q_DOT_S_PPP,PART_INDEX)

! Calculate the solid phase reaction. Return heat and mass generation rates per unit volume.

! N_MATS = Number of MATerialS
! MATL_INDEX(1:N_MATS) = Indices of the materials from the master material list.
! SURF_INDEX = Index of surface, used only for liquids
! (IIG,JJG,KKG) = Indices of nearest gas phase cell
! TMP_S = Solid temperature (K)
! TMP_F = Solid surface temperature (K)
! RHO_S(1:N_MATS) = Array of component densities (kg/m3)
! RHO_S0 = Original solid density (kg/m3)
! DEPTH = Distance from surface (m)
! DT_BC = Time step used by the solid phase solver (s)
! M_DOT_G_PPP_ADJUST(1:N_TRACKED_SPECIES) = Adjusted mass generation rate per unit volume of the gas species
! M_DOT_G_PPP_ACTUAL(1:N_TRACKED_SPECIES) = Actual mass generation rate per unit volume of the gas species
! M_DOT_S_PPP(1:N_MATS) = Mass generation/depletion rate per unit volume of solid components (kg/m3/s)
! Q_DOT_S_PPP = Heat release rate per unit volume (W/m3)
! PART_INDEX = Optional Lagrangian particle index for vegetation

USE PHYSICAL_FUNCTIONS, ONLY: GET_MASS_FRACTION,GET_MOLECULAR_WEIGHT,GET_VISCOSITY
USE MATH_FUNCTIONS, ONLY: INTERPOLATE1D_UNIFORM
INTEGER, INTENT(IN) :: N_MATS,SURF_INDEX,IIG,JJG,KKG
INTEGER, INTENT(IN), OPTIONAL :: PART_INDEX
REAL(EB), INTENT(IN) :: TMP_S,TMP_F,RHO_S0,DT_BC,DEPTH
REAL(EB), DIMENSION(:) :: RHO_S(N_MATS),ZZ_GET(1:N_TRACKED_SPECIES)
REAL(EB), DIMENSION(:), INTENT(OUT) :: M_DOT_G_PPP_ADJUST(N_TRACKED_SPECIES),M_DOT_G_PPP_ACTUAL(N_TRACKED_SPECIES),&
                                       M_DOT_S_PPP(MAX_MATERIALS)
INTEGER, INTENT(IN), DIMENSION(:) :: MATL_INDEX(N_MATS)
INTEGER :: N,NN,J,NS,SMIX_PTR
TYPE(MATERIAL_TYPE), POINTER :: ML
TYPE(SURFACE_TYPE), POINTER :: SF
REAL(EB) :: DTMP,REACTION_RATE,Y_O2,X_O2,Q_DOT_S_PPP,MW_G,X_G,X_W,D_AIR,H_MASS,RE_L,SHERWOOD,MFLUX,MU_AIR,SC_AIR,U_TANG,&
            SIGMA_BETA,RHO_DOT

Q_DOT_S_PPP = 0._EB
M_DOT_S_PPP = 0._EB
M_DOT_G_PPP_ADJUST = 0._EB
M_DOT_G_PPP_ACTUAL = 0._EB

MATERIAL_LOOP: DO N=1,N_MATS  ! Tech Guide: Sum over the materials, alpha

   IF (RHO_S(N) <= 0._EB) CYCLE MATERIAL_LOOP  ! If component alpha density is zero, go on to the next material.
   ML => MATERIAL(MATL_INDEX(N))

   REACTION_LOOP: DO J=1,ML%N_REACTIONS  ! Tech Guide: Sum over the reactions, beta

      SELECT CASE (ML%PYROLYSIS_MODEL)

         CASE (PYROLYSIS_LIQUID)

            SF => SURFACE(SURF_INDEX)
            SMIX_PTR = MAXLOC(ML%NU_GAS(:,1),1)
            ZZ_GET(1:N_TRACKED_SPECIES) = MAX(0._EB,ZZ(IIG,JJG,KKG,1:N_TRACKED_SPECIES))
            CALL GET_MOLECULAR_WEIGHT(ZZ_GET,MW_G)
            X_G = ZZ_GET(SMIX_PTR)/SPECIES_MIXTURE(SMIX_PTR)%MW*MW_G
            X_W = MIN(1._EB-TWO_EPSILON_EB,EXP(ML%H_R(1)*SPECIES_MIXTURE(SMIX_PTR)%MW/R0*(1._EB/ML%TMP_BOIL-1._EB/TMP_F)))
            IF (SF%HM_FIXED>=0._EB) THEN
               H_MASS = SF%HM_FIXED
            ELSEIF (DNS) THEN
               CALL INTERPOLATE1D_UNIFORM(LBOUND(D_Z(:,SMIX_PTR),1),D_Z(:,SMIX_PTR),TMP(IIG,JJG,KKG),D_AIR)
               H_MASS = 2._EB*D_AIR*(RDX(IIG)*RDY(JJG)*RDZ(KKG))**ONTH
            ELSE
               CALL GET_VISCOSITY(ZZ_GET,MU_AIR,TMP(IIG,JJG,KKG))
               U_TANG   = SQRT(2._EB*KRES(IIG,JJG,KKG))
               RE_L     = MAX(5.E5_EB,RHO(IIG,JJG,KKG)*U_TANG*SF%CONV_LENGTH/MU_AIR)
               SC_AIR   = 0.6_EB     ! NU_AIR/D_AIR (Incropera & DeWitt, Chap 7, External Flow)
               SHERWOOD = 0.037_EB*SC_AIR**ONTH*RE_L**0.8_EB
               H_MASS   = SHERWOOD*MU_AIR/(RHO(IIG,JJG,KKG)*SC*SF%CONV_LENGTH)
            ENDIF
            MFLUX = MAX(0._EB,SPECIES_MIXTURE(SMIX_PTR)%MW/R0/TMP_F*H_MASS*LOG((X_G-1._EB)/(X_W-1._EB)))
            MFLUX = MFLUX * PBAR(KKG,PRESSURE_ZONE(IIG,JJG,KKG))
            RHO_DOT = MIN(MFLUX/DX_S(1),ML%RHO_S/DT_BC)  ! kg/m3/s

         CASE (PYROLYSIS_SOLID)

            ! Reaction rate in 1/s (Tech Guide: r_alpha_beta)
            REACTION_RATE = ML%A(J)*(RHO_S(N)/RHO_S0)**ML%N_S(J)*EXP(-ML%E(J)/(R0*TMP_S))
            ! power term
            DTMP = ML%THR_SIGN(J)*(TMP_S-ML%TMP_THR(J))
            IF (ABS(ML%N_T(J))>=TWO_EPSILON_EB) THEN
               IF (DTMP > 0._EB) THEN
                  REACTION_RATE = REACTION_RATE * DTMP**ML%N_T(J)
               ELSE
                  REACTION_RATE = 0._EB
               ENDIF
            ELSE ! threshold
               IF (DTMP < 0._EB) REACTION_RATE = 0._EB
            ENDIF
            ! Phase change reaction?
            IF (ML%PCR(J)) REACTION_RATE = REACTION_RATE / ((ABS(ML%H_R(J))/1000._EB) * DT_BC)
            ! Oxidation reaction?
            IF ( (ML%N_O2(J)>0._EB) .AND. (O2_INDEX > 0)) THEN
               ! Get oxygen mass fraction
               ZZ_GET(1:N_TRACKED_SPECIES) = MAX(0._EB,ZZ(IIG,JJG,KKG,1:N_TRACKED_SPECIES))
               CALL GET_MASS_FRACTION(ZZ_GET,O2_INDEX,Y_O2)
               ! Calculate oxygen volume fraction in the gas cell
               X_O2 = SPECIES(O2_INDEX)%RCON*Y_O2/RSUM(IIG,JJG,KKG)
               ! Calculate oxygen concentration inside the material, assuming decay function
               X_O2 = X_O2 * EXP(-DEPTH/(TWO_EPSILON_EB+ML%GAS_DIFFUSION_DEPTH(J)))
               REACTION_RATE = REACTION_RATE * X_O2**ML%N_O2(J)
            ENDIF
            RHO_DOT  = MIN(RHO_S0*REACTION_RATE,RHO_S(N)/DT_BC)  ! Tech Guide: rho_s(0)*r_alpha,beta kg/m3/s

         CASE (PYROLYSIS_VEGETATION)

            SF => SURFACE(SURF_INDEX)
            LP => LAGRANGIAN_PARTICLE(PART_INDEX)
            ! Tech Guide: r_alpha,beta (1/s)
            REACTION_RATE = ML%A(J)*(RHO_S(N)/RHO_S0)**ML%N_S(J)*EXP(-ML%E(J)/(R0*TMP_S))
            ! power term
            IF (ABS(ML%N_T(J))>=TWO_EPSILON_EB) REACTION_RATE = REACTION_RATE * TMP_S**ML%N_T(J)
            ! Oxidation reaction?
            IF ( (ML%NU_O2(J)>0._EB) .AND. (O2_INDEX > 0)) THEN
               ! Get oxygen mass fraction
               ZZ_GET(1:N_TRACKED_SPECIES) = MAX(0._EB,ZZ(IIG,JJG,KKG,1:N_TRACKED_SPECIES))
               CALL GET_MASS_FRACTION(ZZ_GET,O2_INDEX,Y_O2)
               SIGMA_BETA = LP%PWT*LP%ONE_D%AREA*RDX(IIG)*RDY(JJG)*RDZ(KKG)
               CALL GET_VISCOSITY(ZZ_GET,MU_AIR,TMP(IIG,JJG,KKG))
               U_TANG = SQRT(2._EB*KRES(IIG,JJG,KKG))
               RE_L   = RHO(IIG,JJG,KKG)*U_TANG*2._EB*SF%THICKNESS/MU_AIR
               REACTION_RATE = REACTION_RATE * &
                               RHO(IIG,JJG,KKG)*Y_O2*SIGMA_BETA*(1._EB+ML%BETA_CHAR(J)*SQRT(RE_L))/(RHO_S0*ML%NU_O2(J))
            ENDIF
            RHO_DOT  = MIN(RHO_S0*REACTION_RATE , RHO_S(N)/DT_BC)  ! Tech Guide: rho_s(0)*r_alpha,beta kg/m3/s

      END SELECT

      RHO_S(N) = MAX( 0._EB , RHO_S(N) - DT_BC*RHO_DOT )  ! Tech Guide: rho_s,alpha_new = rho_s,alpha_old-dt*rho_s(0)*r_alpha,beta
      DO NN=1,N_MATS  ! Loop over other materials, looking for the residue (alpha' represents the other materials)
         ! Tech Guide: rho_s,alpha'_new = rho_s,alpha'_old + rho_s(0)*nu_alpha',alpha,beta*r_alpha,beta
         RHO_S(NN) = RHO_S(NN) + ML%NU_RESIDUE(MATL_INDEX(NN),J)*DT_BC*RHO_DOT
         M_DOT_S_PPP(NN) = M_DOT_S_PPP(NN) + ML%NU_RESIDUE(MATL_INDEX(NN),J)*RHO_DOT  ! (m_dot_alpha')'''
      ENDDO
      Q_DOT_S_PPP = Q_DOT_S_PPP - RHO_DOT * ML%H_R(J)  ! Tech Guide: q_dot_s,c'''
      M_DOT_S_PPP(N) = M_DOT_S_PPP(N) - RHO_DOT  ! m_dot_alpha''' = -rho_s(0) * sum_beta r_alpha,beta
      DO NS=1,N_TRACKED_SPECIES  ! Tech Guide: m_dot_gamma'''
         M_DOT_G_PPP_ADJUST(NS) = M_DOT_G_PPP_ADJUST(NS) + ML%ADJUST_BURN_RATE(NS,J)*ML%NU_GAS(NS,J)*RHO_DOT
         M_DOT_G_PPP_ACTUAL(NS) = M_DOT_G_PPP_ACTUAL(NS) + ML%NU_GAS(NS,J)*RHO_DOT
      ENDDO

   ENDDO REACTION_LOOP

ENDDO MATERIAL_LOOP

END SUBROUTINE PYROLYSIS


REAL(EB) FUNCTION HEAT_TRANSFER_COEFFICIENT(DELTA_TMP,H_FIXED,SURF_INDEX,WALL_INDEX,CFACE_INDEX,PARTICLE_INDEX)

! Compute the convective heat transfer coefficient

USE TURBULENCE, ONLY: LOGLAW_HEAT_FLUX_MODEL,ABL_HEAT_FLUX_MODEL,NATURAL_CONVECTION_MODEL,FORCED_CONVECTION_MODEL,&
                      RAYLEIGH_HEAT_FLUX_MODEL,YUAN_HEAT_FLUX_MODEL
USE PHYSICAL_FUNCTIONS, ONLY: GET_CONDUCTIVITY,GET_VISCOSITY,GET_SPECIFIC_HEAT
REAL(EB), INTENT(IN) :: DELTA_TMP,H_FIXED
INTEGER, INTENT(IN) :: SURF_INDEX
INTEGER, INTENT(IN), OPTIONAL :: WALL_INDEX,PARTICLE_INDEX,CFACE_INDEX
INTEGER  :: ITMP
REAL(EB) :: RE,H_NATURAL,H_FORCED,NUSSELT,FRICTION_VELOCITY,YPLUS,ZSTAR,DN,TMP_FILM,MU_G,K_G,CP_G,ZZ_GET(1:N_TRACKED_SPECIES)
TYPE(SURFACE_TYPE), POINTER :: SFX
TYPE(WALL_TYPE), POINTER :: WCX
TYPE(ONE_D_M_AND_E_XFER_TYPE), POINTER :: ONE_DX

SFX => SURFACE(SURF_INDEX)

! If the user wants a specified HTC, set it and return

IF (H_FIXED >= 0._EB) THEN
   HEAT_TRANSFER_COEFFICIENT = H_FIXED
   RETURN
ENDIF

! Determine if this is a particle or wall cell

IF (PRESENT(PARTICLE_INDEX)) THEN
   ONE_DX => LAGRANGIAN_PARTICLE(PARTICLE_INDEX)%ONE_D
   DN = SFX%CONV_LENGTH
ELSEIF (PRESENT(WALL_INDEX)) THEN
   WCX   => WALL(WALL_INDEX)
   ONE_DX => WALL(WALL_INDEX)%ONE_D
   FRICTION_VELOCITY = WCX%ONE_D%U_TAU
   YPLUS = WCX%ONE_D%Y_PLUS
   DN = 1._EB/ONE_DX%RDN
ELSEIF (PRESENT(CFACE_INDEX)) THEN
   ONE_DX => CFACE(CFACE_INDEX)%ONE_D
   DN = 1._EB/ONE_DX%RDN
ELSE
   HEAT_TRANSFER_COEFFICIENT = SFX%C_VERTICAL*ABS(DELTA_TMP)**ONTH
   RETURN
ENDIF

! If this is a DNS calculation at a solid wall, set HTC and return.

IF ( (DNS .OR. SOLID_PHASE_ONLY) .AND. (PRESENT(WALL_INDEX) .OR. PRESENT(CFACE_INDEX)) ) THEN
   HEAT_TRANSFER_COEFFICIENT = 2._EB * ONE_DX%K_G * ONE_DX%RDN
   RETURN
ENDIF

! Calculate HEAT_TRANSFER_COEFFICIENT

H_NATURAL = 0._EB
H_FORCED  = 0._EB

TMP_FILM = 0.5_EB*(ONE_DX%TMP_G+ONE_DX%TMP_F)
ITMP = MIN(4999,NINT(TMP_FILM))

ZZ_GET(1:N_TRACKED_SPECIES) = ONE_DX%ZZ_G(1:N_TRACKED_SPECIES)

HTC_MODEL_SELECT: SELECT CASE(SFX%HEAT_TRANSFER_MODEL)
   CASE(H_DEFAULT)
      CALL GET_VISCOSITY(ZZ_GET,MU_G,TMP_FILM)
      CALL GET_CONDUCTIVITY(ZZ_GET,K_G,TMP_FILM)
      RE = ONE_DX%RHO_G*ONE_DX%U_TANG*SFX%CONV_LENGTH/MU_G
      CALL FORCED_CONVECTION_MODEL(H_FORCED,RE,K_G,SFX%CONV_LENGTH,SFX%GEOMETRY)
      CALL NATURAL_CONVECTION_MODEL(H_NATURAL,DELTA_TMP,SFX%C_VERTICAL,SFX%C_HORIZONTAL,SFX%GEOMETRY,ONE_DX%IOR,K_G,DN)
   CASE(H_LOGLAW)
      H_NATURAL = 0._EB
      CALL GET_VISCOSITY(ZZ_GET,MU_G,TMP_FILM)
      CALL GET_CONDUCTIVITY(ZZ_GET,K_G,TMP_FILM)
      CALL GET_SPECIFIC_HEAT(ZZ_GET,CP_G,TMP_FILM)
      CALL LOGLAW_HEAT_FLUX_MODEL(H_FORCED,YPLUS,FRICTION_VELOCITY,K_G,ONE_DX%RHO_G,CP_G,MU_G)
   CASE(H_ABL)
      H_NATURAL = 0._EB
      CALL GET_SPECIFIC_HEAT(ZZ_GET,CP_G,TMP_FILM)
      CALL ABL_HEAT_FLUX_MODEL(H_FORCED,FRICTION_VELOCITY,DN,SFX%ROUGHNESS,ONE_DX%TMP_G,ONE_DX%TMP_F,ONE_DX%RHO_G,CP_G)
   CASE(H_RAYLEIGH)
      H_FORCED = 0._EB
      CALL GET_VISCOSITY(ZZ_GET,MU_G,TMP_FILM)
      CALL GET_CONDUCTIVITY(ZZ_GET,K_G,TMP_FILM)
      CALL GET_SPECIFIC_HEAT(ZZ_GET,CP_G,TMP_FILM)
      CALL RAYLEIGH_HEAT_FLUX_MODEL(H_NATURAL,ZSTAR,DN,ONE_DX%TMP_F,ONE_DX%TMP_G,K_G,ONE_DX%RHO_G,CP_G,MU_G)
      IF (PRESENT(WALL_INDEX) .OR. PRESENT(CFACE_INDEX)) ONE_DX%Z_STAR = ZSTAR
   CASE(H_YUAN)
      H_FORCED = 0._EB
      CALL GET_CONDUCTIVITY(ZZ_GET,K_G,TMP_FILM)
      CALL GET_SPECIFIC_HEAT(ZZ_GET,CP_G,TMP_FILM)
      CALL YUAN_HEAT_FLUX_MODEL(H_NATURAL,ZSTAR,DN,ONE_DX%TMP_F,ONE_DX%TMP_G,K_G,ONE_DX%RHO_G,CP_G)
      IF (PRESENT(WALL_INDEX) .OR. PRESENT(CFACE_INDEX)) ONE_DX%Z_STAR = ZSTAR
   CASE(H_CUSTOM)
      CALL GET_VISCOSITY(ZZ_GET,MU_G,TMP_FILM)
      CALL GET_CONDUCTIVITY(ZZ_GET,K_G,TMP_FILM)
      RE = ONE_DX%RHO_G*ONE_DX%U_TANG*SFX%CONV_LENGTH/MU_G
      NUSSELT = MAX(1._EB,SFX%C_FORCED_CONSTANT+SFX%C_FORCED_RE*RE**SFX%C_FORCED_RE_EXP*PR_AIR**SFX%C_FORCED_PR_EXP)
      H_FORCED = NUSSELT*K_G/SFX%CONV_LENGTH
      CALL NATURAL_CONVECTION_MODEL(H_NATURAL,DELTA_TMP,SFX%C_VERTICAL,SFX%C_HORIZONTAL,SFX%GEOMETRY,ONE_DX%IOR,K_G,DN)
END SELECT HTC_MODEL_SELECT

HEAT_TRANSFER_COEFFICIENT = MAX(H_FORCED,H_NATURAL)

END FUNCTION HEAT_TRANSFER_COEFFICIENT


SUBROUTINE TGA_ANALYSIS

! This routine performs a numerical TGA (thermo-gravimetric analysis) at the start of the simulation

USE PHYSICAL_FUNCTIONS, ONLY: SURFACE_DENSITY
USE COMP_FUNCTIONS, ONLY: SHUTDOWN
REAL(EB) :: DT_TGA=0.01_EB,T_TGA,SURF_DEN,SURF_DEN_0,HRR
INTEGER :: N_TGA,I,IW,IP
CHARACTER(80) :: MESSAGE,TCFORM

CALL POINT_TO_MESH(1)

RADIATION = .FALSE.
TGA_HEATING_RATE = TGA_HEATING_RATE/60._EB  ! K/min --> K/s
TGA_FINAL_TEMPERATURE = TGA_FINAL_TEMPERATURE + TMPM  ! C --> K
I_RAMP_AGT = 0
N_TGA = NINT((TGA_FINAL_TEMPERATURE-TMPA)/(TGA_HEATING_RATE*DT_TGA))
T_TGA = 0._EB

IF (TGA_WALL_INDEX>0) THEN
   IW = TGA_WALL_INDEX
   ONE_D => WALL(IW)%ONE_D
ELSEIF (TGA_PARTICLE_INDEX>0) THEN
   IP = TGA_PARTICLE_INDEX
   ONE_D => LAGRANGIAN_PARTICLE(IP)%ONE_D
ELSE
   WRITE(MESSAGE,'(A)') 'ERROR: No wall or particle to which to apply the TGA analysis'
   CALL SHUTDOWN(MESSAGE) ; RETURN
ENDIF

OPEN (LU_TGA,FILE=FN_TGA,FORM='FORMATTED',STATUS='REPLACE')
WRITE(LU_TGA,'(A)') 's,C,g/g,1/s,W/g,W/g'
WRITE(LU_TGA,'(A)') 'Time,Temp,Mass,MLR,MCC,DSC'

SURF_DEN_0 = SURFACE(TGA_SURF_INDEX)%SURFACE_DENSITY
WRITE(TCFORM,'(5A)') "(5(",TRIM(FMT_R),",','),",TRIM(FMT_R),")"

DO I=1,N_TGA
   IF (ONE_D%LAYER_THICKNESS(1)<TWO_EPSILON_EB) EXIT
   T_TGA = I*DT_TGA
   ASSUMED_GAS_TEMPERATURE = TMPA + TGA_HEATING_RATE*T_TGA
   IF (TGA_WALL_INDEX>0) THEN
      CALL SOLID_HEAT_TRANSFER_1D(1,T_TGA,DT_TGA,WALL_INDEX=IW)
      SURF_DEN = SURFACE_DENSITY(1,0,WALL_INDEX=IW)
   ELSE
      CALL SOLID_HEAT_TRANSFER_1D(1,T_TGA,DT_TGA,PARTICLE_INDEX=IP)
      SURF_DEN = SURFACE_DENSITY(1,0,LAGRANGIAN_PARTICLE_INDEX=IP)
   ENDIF
   IF (MOD(I,NINT(1._EB/(TGA_HEATING_RATE*DT_TGA)))==0) THEN
      IF (N_REACTIONS>0) THEN
         HRR = ONE_D%MASSFLUX(REACTION(1)%FUEL_SMIX_INDEX)*0.001*REACTION(1)%HEAT_OF_COMBUSTION/(ONE_D%AREA_ADJUST*SURF_DEN_0)
      ELSE
         HRR = 0._EB
      ENDIF
      WRITE(LU_TGA,TCFORM) REAL(T_TGA,FB), REAL(ONE_D%TMP_F-TMPM,FB), REAL(SURF_DEN/SURF_DEN_0,FB), &
                           REAL(SUM(ONE_D%MASSFLUX_SPEC(1:N_TRACKED_SPECIES))/SURF_DEN_0,FB), &
                           REAL(HRR,FB), REAL(ONE_D%QCONF*0.001_EB/SURF_DEN_0,FB)
   ENDIF
ENDDO

CLOSE(LU_TGA)

END SUBROUTINE TGA_ANALYSIS

END MODULE WALL_ROUTINES
