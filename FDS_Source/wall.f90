MODULE WALL_ROUTINES
 
! Compute the wall boundary conditions
 
USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_POINTERS
 
IMPLICIT NONE
PRIVATE
CHARACTER(255), PARAMETER :: wallid='$Id$'
CHARACTER(255), PARAMETER :: wallrev='$Revision$'
CHARACTER(255), PARAMETER :: walldate='$Date$'

PUBLIC WALL_BC,GET_REV_wall
 

CONTAINS


SUBROUTINE WALL_BC(T,NM)

! This is the main control routine for this module

USE COMP_FUNCTIONS, ONLY: SECOND
REAL(EB) :: TNOW
REAL(EB), INTENT(IN) :: T
INTEGER, INTENT(IN) :: NM

IF (EVACUATION_ONLY(NM)) RETURN

TNOW=SECOND()

CALL POINT_TO_MESH(NM)

CALL THERMAL_BC(T)
IF (SOOT_DEPOSITION .AND. CORRECTOR) CALL CALC_SOOT_DEPOSITION(NM)
CALL SPECIES_BC(T)
CALL DENSITY_BC
IF (HVAC_SOLVE) CALL HVAC_BC

TUSED(6,NM)=TUSED(6,NM)+SECOND()-TNOW
END SUBROUTINE WALL_BC



SUBROUTINE THERMAL_BC(T)

! Thermal boundary conditions for adiabatic, fixed temperature, fixed flux and interpolated boundaries.
! One dimensional heat transfer and pyrolysis is done in PYROLYSIS, which is called at the end of this routine.

USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP 
USE PHYSICAL_FUNCTIONS, ONLY : GET_SPECIFIC_GAS_CONSTANT
REAL(EB) :: DT_BC,T,TSI,TMP_G,RHO_G,DTMP,TMP_OTHER,RAMP_FACTOR,QNET,FDERIV,TMP_EXTERIOR,UN,ARO,UWO,QEXTRA,RSUM_W, &
            YY_G_ALL(MAX_SPECIES),RHO_YY_F(MAX_SPECIES),YY_GET(1:N_SPECIES)
INTEGER  :: IOR,II,JJ,KK,IBC,IIG,JJG,KKG,IW,NOM,IIO,JJO,KKO
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU=>NULL(),VV=>NULL(),WW=>NULL(),RHOP=>NULL(),OM_RHOP=>NULL()
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: YYP=>NULL(),OM_YYP=>NULL()
LOGICAL :: INFLOW
TYPE (SURFACE_TYPE), POINTER :: SF=>NULL()
TYPE (VENTS_TYPE), POINTER :: VT=>NULL()
TYPE (OMESH_TYPE), POINTER :: OM=>NULL()
TYPE (MESH_TYPE), POINTER :: MM=>NULL()
REAL(EB), POINTER, DIMENSION(:,:) :: PBAR_P=>NULL()

IF (PREDICTOR) THEN
   UU => U
   VV => V
   WW => W
   RHOP => RHOS
   YYP  => YYS
   PBAR_P => PBAR_S
ELSE
   UU => US
   VV => VS
   WW => WS
   RHOP => RHO
   YYP  => YY   
   PBAR_P => PBAR
ENDIF
 
! Loop through all boundary cells and apply heat transfer method, except for thermally-thick cells
 
HEAT_FLUX_LOOP: DO IW=1,NWC+NVWC

   IF (BOUNDARY_TYPE(IW)==NULL_BOUNDARY) CYCLE HEAT_FLUX_LOOP

   II  = IJKW(1,IW)
   JJ  = IJKW(2,IW)
   KK  = IJKW(3,IW)
   IIG = IJKW(6,IW)
   JJG = IJKW(7,IW)
   KKG = IJKW(8,IW)
   IOR = IJKW(4,IW)
   IBC = IJKW(5,IW)

   ! Consider special boundary conditions

   IF (BOUNDARY_TYPE(IW)/=SOLID_BOUNDARY) THEN
      IF (BOUNDARY_TYPE(IW)==INTERPOLATED_BOUNDARY) IBC = INTERPOLATED_SURF_INDEX
      IF (BOUNDARY_TYPE(IW)==OPEN_BOUNDARY)         IBC = OPEN_SURF_INDEX
      IF (BOUNDARY_TYPE(IW)==MIRROR_BOUNDARY)       IBC = MIRROR_SURF_INDEX
   ENDIF

   IF (BOUNDARY_TYPE(IW)==POROUS_BOUNDARY) THEN
      IF (IJKW(9,IW)/=0) THEN
         IBC = INTERPOLATED_SURF_INDEX
      ELSE
         CYCLE HEAT_FLUX_LOOP
      ENDIF
   ENDIF

   ! Choose the SURFace type

   SF  => SURFACE(IBC)

   ! Compute surface temperature, TMP_F, and convective heat flux, QCONF, for various boundary conditions

   METHOD_OF_HEAT_TRANSFER: SELECT CASE(SF%THERMAL_BC_INDEX)
 
      CASE (NO_CONVECTION) METHOD_OF_HEAT_TRANSFER

         TMP_F(IW) = TMP(IIG,JJG,KKG)

      CASE (INFLOW_OUTFLOW) METHOD_OF_HEAT_TRANSFER 
 
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
         IF (UN>=0._EB) INFLOW = .TRUE.

         IF (INFLOW) THEN
            TMP_EXTERIOR = TMP_0(KK)
            IF (VENT_INDEX(IW)>0) THEN
               VT => VENTS(VENT_INDEX(IW))
               IF (VT%TMP_EXTERIOR>0._EB) TMP_EXTERIOR = VT%TMP_EXTERIOR
            ENDIF
            TMP_F(IW) = TMP_EXTERIOR
            HEAT_TRANS_COEF(IW) = 0._EB
            QCONF(IW) = 0._EB
         ELSE
            TMP_F(IW) = TMP(IIG,JJG,KKG)
         ENDIF
         TMP(II,JJ,KK) = TMP_F(IW)
 
      CASE (SPECIFIED_TEMPERATURE) METHOD_OF_HEAT_TRANSFER

         TMP_G = TMP(IIG,JJG,KKG)
         IF (TW(IW)==T_BEGIN .AND. SF%RAMP_INDEX(TIME_TEMP)>=1) THEN
            TSI = T
         ELSE
            TSI = T - TW(IW)
         ENDIF
         IF (UW(IW)<=0._EB) THEN
            TMP_F(IW) = TMP_0(KK) + EVALUATE_RAMP(TSI,SF%TAU(TIME_TEMP),SF%RAMP_INDEX(TIME_TEMP))*(SF%TMP_FRONT-TMP_0(KK))
         ELSE
            TMP_F(IW) = TMP_G  ! If gas is being drawn from the domain, set the boundary temperature to the gas temperature
         ENDIF
         DTMP = TMP_G - TMP_F(IW)
         HEAT_TRANS_COEF(IW) = HEAT_TRANSFER_COEFFICIENT(IW,IIG,JJG,KKG,IOR,TMP_G,DTMP,SF%H_FIXED)
         QCONF(IW) = HEAT_TRANS_COEF(IW)*DTMP
         
      CASE (NET_FLUX_BC) METHOD_OF_HEAT_TRANSFER
         
         IF (TW(IW)==T_BEGIN .AND. SF%RAMP_INDEX(TIME_HEAT)>=1) THEN
            TSI = T
         ELSE
            TSI = T - TW(IW)
         ENDIF
         TMP_G = TMP(IIG,JJG,KKG)
         TMP_OTHER = TMP_F(IW)
         RAMP_FACTOR = EVALUATE_RAMP(TSI,SF%TAU(TIME_HEAT),SF%RAMP_INDEX(TIME_HEAT))
         QNET = -RAMP_FACTOR*SF%NET_HEAT_FLUX*AREA_ADJUST(IW)
         ADLOOP: DO
            DTMP = TMP_G - TMP_OTHER
            HEAT_TRANS_COEF(IW) = HEAT_TRANSFER_COEFFICIENT(IW,IIG,JJG,KKG,IOR,TMP_G,DTMP,SF%H_FIXED)
            IF (RADIATION) THEN
               QEXTRA = HEAT_TRANS_COEF(IW)*DTMP + QRADIN(IW) - E_WALL(IW) * SIGMA * TMP_OTHER ** 4 - QNET
               FDERIV = -HEAT_TRANS_COEF(IW) -  4._EB * E_WALL(IW) * SIGMA * TMP_OTHER ** 3
            ELSE
               QEXTRA = HEAT_TRANS_COEF(IW)*DTMP
               FDERIV = -HEAT_TRANS_COEF(IW)
            ENDIF
            IF (FDERIV /= 0._EB) TMP_OTHER = TMP_OTHER - QEXTRA / FDERIV
            IF (ABS(TMP_OTHER - TMP_F(IW)) / TMP_F(IW) < 0.0001) THEN
               TMP_F(IW) = TMP_OTHER
               EXIT ADLOOP
            ELSE
               TMP_F(IW) = TMP_OTHER
               CYCLE ADLOOP
            ENDIF
         ENDDO ADLOOP
         DTMP = TMP_G - TMP_F(IW)
         HEAT_TRANS_COEF(IW) = HEAT_TRANSFER_COEFFICIENT(IW,IIG,JJG,KKG,IOR,TMP_G,DTMP,SF%H_FIXED)
         QCONF(IW) = HEAT_TRANS_COEF(IW)*DTMP

      CASE (CONVECTIVE_FLUX_BC) METHOD_OF_HEAT_TRANSFER
      
         IF (TW(IW)==T_BEGIN .AND. SF%RAMP_INDEX(TIME_HEAT)>=1) THEN
            TSI = T
         ELSE
            TSI = T - TW(IW)
         ENDIF
         RAMP_FACTOR = EVALUATE_RAMP(TSI,SF%TAU(TIME_HEAT),SF%RAMP_INDEX(TIME_HEAT))
         TMP_F(IW) = TMPA + RAMP_FACTOR*(SF%TMP_FRONT-TMPA)
         QCONF(IW) = -RAMP_FACTOR*SF%CONVECTIVE_HEAT_FLUX*AREA_ADJUST(IW)
 
      CASE (INTERPOLATED_BC) METHOD_OF_HEAT_TRANSFER
 
         NOM  =  IJKW(9,IW)
         OM   => OMESH(NOM)
         IF (PREDICTOR) THEN
            OM_RHOP => OM%RHOS
            IF (N_SPECIES>0) OM_YYP => OM%YYS
         ELSE
            OM_RHOP => OM%RHO
            IF (N_SPECIES>0) OM_YYP => OM%YY
         ENDIF
         MM    => MESHES(NOM)
         RHO_G = RHOP(IIG,JJG,KKG)
         RHO_F(IW) = RHO_G  ! Initialize face value of RHO with RHO_G
         IF (N_SPECIES>0) THEN
            YY_G_ALL(1:N_SPECIES) = YYP(IIG,JJG,KKG,1:N_SPECIES)
            RHO_YY_F(1:N_SPECIES) = RHO_G*YY_G_ALL(1:N_SPECIES) ! Initialize face value of RHO_YY with RHO_G*YY_G
         ENDIF

         DO KKO=IJKW(12,IW),IJKW(15,IW)
            DO JJO=IJKW(11,IW),IJKW(14,IW)
               DO IIO=IJKW(10,IW),IJKW(13,IW)
                  SELECT CASE(IOR)
                     CASE( 1)
                        ARO = MIN(1._EB , RDY(JJ)*RDZ(KK)*MM%DY(JJO)*MM%DZ(KKO)) * 2.*DX(II)/(MM%DX(IIO)+DX(II))
                        UWO = -OM%U(IIO,JJO,KKO)
                     CASE(-1)
                        ARO = MIN(1._EB , RDY(JJ)*RDZ(KK)*MM%DY(JJO)*MM%DZ(KKO)) * 2.*DX(II)/(MM%DX(IIO)+DX(II))
                        UWO =  OM%U(IIO-1,JJO,KKO)
                     CASE( 2)
                        ARO = MIN(1._EB , RDX(II)*RDZ(KK)*MM%DX(IIO)*MM%DZ(KKO)) * 2.*DY(JJ)/(MM%DY(JJO)+DY(JJ))
                        UWO = -OM%V(IIO,JJO,KKO)
                     CASE(-2)
                        ARO = MIN(1._EB , RDX(II)*RDZ(KK)*MM%DX(IIO)*MM%DZ(KKO)) * 2.*DY(JJ)/(MM%DY(JJO)+DY(JJ))
                        UWO =  OM%V(IIO,JJO-1,KKO)
                     CASE( 3)
                        ARO = MIN(1._EB , RDX(II)*RDY(JJ)*MM%DX(IIO)*MM%DY(JJO)) * 2.*DZ(KK)/(MM%DZ(KKO)+DZ(KK))
                        UWO = -OM%W(IIO,JJO,KKO)
                     CASE(-3)
                        ARO = MIN(1._EB , RDX(II)*RDY(JJ)*MM%DX(IIO)*MM%DY(JJO)) * 2.*DZ(KK)/(MM%DZ(KKO)+DZ(KK))
                        UWO =  OM%W(IIO,JJO,KKO-1)
                  END SELECT
                  RHO_F(IW)     = RHO_F(IW) + 0.5_EB*ARO*(OM_RHOP(IIO,JJO,KKO)-RHO_G)
                  IF (N_SPECIES>0) RHO_YY_F(1:N_SPECIES) = RHO_YY_F(1:N_SPECIES) + &
                                   0.5_EB*ARO*(OM_RHOP(IIO,JJO,KKO)*OM_YYP(IIO,JJO,KKO,1:N_SPECIES)-RHO_G*YY_G_ALL(1:N_SPECIES))
               ENDDO
            ENDDO
         ENDDO
         RHOP(II,JJ,KK) = 2._EB*RHO_F(IW) - RHO_G
         IF (N_SPECIES==0) THEN
            TMP(II,JJ,KK) = PBAR_P(KK,PRESSURE_ZONE_WALL(IW))/(SPECIES(0)%RCON*RHOP(II,JJ,KK))
         ELSE
            YY_F(IW,1:N_SPECIES)      = RHO_YY_F(1:N_SPECIES)/RHO_F(IW)
            YYP(II,JJ,KK,1:N_SPECIES) = 2._EB*YY_F(IW,1:N_SPECIES) - YY_G_ALL(1:N_SPECIES)
            YY_GET(1:N_SPECIES) = MAX(0._EB,YYP(II,JJ,KK,1:N_SPECIES))
            CALL GET_SPECIFIC_GAS_CONSTANT(YY_GET,RSUM_W)
            TMP(II,JJ,KK) = PBAR_P(KK,PRESSURE_ZONE_WALL(IW))/(RSUM_W*RHOP(II,JJ,KK))
         ENDIF
         QCONF(IW) = 0._EB
         TMP_F(IW) = 0.5*(TMP(IIG,JJG,KKG)+TMP(II,JJ,KK))
 
   END SELECT METHOD_OF_HEAT_TRANSFER
 
ENDDO HEAT_FLUX_LOOP
 
! For thermally-thick boundary conditions, call the routine PYROLYSIS
 
IF (CORRECTOR) THEN
   WALL_COUNTER = WALL_COUNTER + 1
   IF (WALL_COUNTER==WALL_INCREMENT) THEN
      DT_BC    = T - BC_CLOCK
      BC_CLOCK = T
      CALL PYROLYSIS(T,DT_BC)
      WALL_COUNTER = 0
   ENDIF
ENDIF
 
END SUBROUTINE THERMAL_BC
 
 

SUBROUTINE SPECIES_BC(T)

! Compute the species mass fractions at the boundary, YY_F

USE MATH_FUNCTIONS, ONLY : EVALUATE_RAMP
REAL(EB) :: T,YY_EXTERIOR(MAX_SPECIES),YY_G,RHO_G,UN,DD,MFT,TSI,RADIUS,AREA_SCALING,&
            AREA,RVC,RHO_NEW,M_DOT_PPP,CP,MW_RATIO,H_G,DELTA_H_G
INTEGER :: I,IBC,IIG,JJG,KKG,IOR,IW,II,JJ,KK,N,ITMP
TYPE (SURFACE_TYPE), POINTER :: SF=>NULL()
TYPE (VENTS_TYPE), POINTER :: VT=>NULL()
TYPE (PARTICLE_CLASS_TYPE), POINTER :: PC=>NULL()
TYPE (REACTION_TYPE), POINTER :: RN=>NULL()
TYPE (MATERIAL_TYPE), POINTER :: ML=>NULL()
TYPE (DROPLET_TYPE), POINTER :: DR=>NULL()
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU=>NULL(),VV=>NULL(),WW=>NULL(),RHOP=>NULL()
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: YYP=>NULL()
 
IF (PREDICTOR) THEN
   UU => U
   VV => V
   WW => W
   RHOP => RHOS
   YYP => YYS
ELSE
   UU => US
   VV => VS
   WW => WS
   RHOP => RHO
   YYP => YY
ENDIF 
 
! Loop through the wall cells, apply mass boundary conditions

WALL_CELL_LOOP: DO IW=1,NWC+NVWC

   IF (BOUNDARY_TYPE(IW)==NULL_BOUNDARY .OR. BOUNDARY_TYPE(IW)==INTERPOLATED_BOUNDARY) CYCLE WALL_CELL_LOOP
   IF (BOUNDARY_TYPE(IW)==POROUS_BOUNDARY) CYCLE WALL_CELL_LOOP

   IBC = IJKW(5,IW)

   ! Special cases that over-ride the boundary condition index, IBC

   IF (BOUNDARY_TYPE(IW)/=SOLID_BOUNDARY) THEN
      IF (BOUNDARY_TYPE(IW)==OPEN_BOUNDARY)         IBC = OPEN_SURF_INDEX
      IF (BOUNDARY_TYPE(IW)==MIRROR_BOUNDARY)       IBC = MIRROR_SURF_INDEX
   ENDIF

   ! Set the SURFace type

   SF  => SURFACE(IBC)

   ! Special cases

   IF (N_SPECIES==0 .AND. .NOT. SF%SPECIES_BC_INDEX==SPECIFIED_MASS_FLUX)                          CYCLE WALL_CELL_LOOP
   IF (N_SPECIES==0 .AND. SF%SPECIES_BC_INDEX==SPECIFIED_MASS_FLUX .AND. SF%MASS_FLUX(0) == 0._EB) CYCLE WALL_CELL_LOOP

   ! Get the wall indices

   II  = IJKW(1,IW)
   JJ  = IJKW(2,IW)
   KK  = IJKW(3,IW)
   IOR = IJKW(4,IW)
   IIG = IJKW(6,IW)
   JJG = IJKW(7,IW)
   KKG = IJKW(8,IW)
   
   ! Check if suppression by water is to be applied

   IF (CORRECTOR .AND. SF%E_COEFFICIENT>0._EB) THEN
      IF (SUM(WMPUA(IW,:))>0._EB .AND. T>TW(IW)) EW(IW) = EW(IW) + SF%E_COEFFICIENT*SUM(WMPUA(IW,:))*DT
   ENDIF

   ! Apply the different species boundary conditions to non-thermally thick solids

   METHOD_OF_MASS_TRANSFER: SELECT CASE(SF%SPECIES_BC_INDEX)
 
      CASE (NO_MASS_FLUX) METHOD_OF_MASS_TRANSFER 

         IF (.NOT.SOLID(CELL_INDEX(IIG,JJG,KKG))) YY_F(IW,1:N_SPECIES) = YYP(IIG,JJG,KKG,1:N_SPECIES)
         IF (BOUNDARY_TYPE(IW)==OPEN_BOUNDARY) THEN
            YY_EXTERIOR(1:N_SPECIES) = SPECIES(1:N_SPECIES)%YY0
            VT => VENTS(VENT_INDEX(IW))
            DO N=1,N_SPECIES
               IF (VT%MASS_FRACTION(N)>-1._EB) YY_EXTERIOR(N) = VT%MASS_FRACTION(N)
            ENDDO
            SELECT CASE(IOR)
               CASE( 1) 
                  IF (UU(II,JJ,KK)>0._EB)   YY_F(IW,1:N_SPECIES)=YY_EXTERIOR(1:N_SPECIES)
               CASE(-1)
                  IF (UU(II-1,JJ,KK)<0._EB) YY_F(IW,1:N_SPECIES)=YY_EXTERIOR(1:N_SPECIES)
               CASE( 2) 
                  IF (VV(II,JJ,KK)>0._EB)   YY_F(IW,1:N_SPECIES)=YY_EXTERIOR(1:N_SPECIES)
               CASE(-2)
                  IF (VV(II,JJ-1,KK)<0._EB) YY_F(IW,1:N_SPECIES)=YY_EXTERIOR(1:N_SPECIES)
               CASE( 3)
                  IF (WW(II,JJ,KK)>0._EB)   YY_F(IW,1:N_SPECIES)=YY_EXTERIOR(1:N_SPECIES)
               CASE(-3)
                  IF (WW(II,JJ,KK-1)<0._EB) YY_F(IW,1:N_SPECIES)=YY_EXTERIOR(1:N_SPECIES)
            END SELECT
            YYP(II,JJ,KK,1:N_SPECIES)=YY_F(IW,1:N_SPECIES)
         ENDIF
         IF ( SF%LEAK_PATH(1)>-1 .AND. UWS(IW)<0._EB) THEN
            RHO_G = RHOP(IIG,JJG,KKG)
            UN    = -UWS(IW)
            SPECIES_LOOP_1: DO N=1,N_SPECIES
               DD    = 2._EB*RHODW(IW,N)*RDN(IW)
               YY_G  = YYP(IIG,JJG,KKG,N)
               YY_F(IW,N) = ( MASSFLUX(IW,N) + YY_G*DD ) / (DD + UN*RHO_F(IW))
            ENDDO SPECIES_LOOP_1
         ENDIF
 
      CASE (SPECIFIED_MASS_FRACTION) METHOD_OF_MASS_TRANSFER

         IF (TW(IW)==T_BEGIN .AND. ANY(SF%RAMP_INDEX>=1)) THEN
            IF (PREDICTOR) TSI = T + DT
            IF (CORRECTOR) TSI = T
         ELSE
            IF (PREDICTOR) TSI = T + DT - TW(IW)
            IF (CORRECTOR) TSI = T      - TW(IW)
         ENDIF

         DO N=1,N_SPECIES
            YY_F(IW,N) = SPECIES(N)%YY0 + EVALUATE_RAMP(TSI,SF%TAU(N),SF%RAMP_INDEX(N))*(SF%MASS_FRACTION(N)-SPECIES(N)%YY0)
         ENDDO
 
      CASE (SPECIFIED_MASS_FLUX) METHOD_OF_MASS_TRANSFER

         ! If the current time is before the "activation" time, TW, apply simple BCs and get out

         IF (T < TW(IW)) THEN
            MASSFLUX(IW,0) = 0._EB
            IF (N_SPECIES > 0)  THEN
               IF (.NOT.SOLID(CELL_INDEX(IIG,JJG,KKG))) YY_F(IW,1:N_SPECIES) = YYP(IIG,JJG,KKG,1:N_SPECIES)
               IF (PREDICTOR) UWS(IW)   = 0._EB 
               MASSFLUX(IW,1:N_SPECIES) = 0._EB
               MASSFLUX_ACTUAL(IW,1:N_SPECIES) = 0._EB
               ENDIF
            CYCLE WALL_CELL_LOOP
         ENDIF

         ! Zero out the running counter of Mass Flux Total (MFT)

         MFT = 0._EB

         ! If the user has specified the burning rate, evaluate the ramp and other related parameters

         SUM_MASSFLUX_LOOP: DO N=0,N_SPECIES
            IF (SF%MASS_FLUX(N) > 0._EB) THEN  ! Use user-specified ramp-up of mass flux
               IF (TW(IW)==T_BEGIN .AND. SF%RAMP_INDEX(N)>=1) THEN
                  TSI = T
               ELSE
                  TSI = T - TW(IW)
               ENDIF
               MASSFLUX(IW,N) = EVALUATE_RAMP(TSI,SF%TAU(N),SF%RAMP_INDEX(N))*SF%MASS_FLUX(N)
               IF (N==I_FUEL) THEN
                  IF (EW(IW)>0._EB) MASSFLUX(IW,N) = MASSFLUX(IW,N)*EXP(-EW(IW))
                  IF (SF%N_MATL==0) THEN
                     MASSFLUX_ACTUAL(IW,N) = MASSFLUX(IW,N)/AREA_ADJUST(IW)
                  ELSE
                     ML => MATERIAL(SF%MATL_INDEX(1))
                     RN => REACTION(1)
                     MASSFLUX_ACTUAL(IW,N) = MASSFLUX(IW,N)*RN%HEAT_OF_COMBUSTION/ML%HEAT_OF_COMBUSTION(1,I_FUEL)/AREA_ADJUST(IW)
                  ENDIF
               ENDIF
            ENDIF
            MASSFLUX(IW,N) = MASSFLUX(IW,N)*AREA_ADJUST(IW)
            MFT = MFT + MASSFLUX(IW,N)
         ENDDO SUM_MASSFLUX_LOOP

         ! Add total consumed mass to various summing arrays

         CONSUME_MASS: IF (CORRECTOR .AND. SF%THERMALLY_THICK) THEN  
            DO N=1,N_SPECIES
               OBSTRUCTION(OBST_INDEX_W(IW))%MASS = OBSTRUCTION(OBST_INDEX_W(IW))%MASS - MASSFLUX_ACTUAL(IW,N)*DT*AW(IW)
            ENDDO
         ENDIF CONSUME_MASS

         ! Compute the ghost cell value of the species to get the right mass flux
 
         RHO_G = RHOP(IIG,JJG,KKG)
         UN    = MFT/RHO_F(IW)
         IF (PREDICTOR) UWS(IW) = -UN
         SPECIES_LOOP: DO N=1,N_SPECIES
            DD    = 2.*RHODW(IW,N)*RDN(IW)
            YY_G  = YYP(IIG,JJG,KKG,N)
            YY_F(IW,N) = ( MASSFLUX(IW,N) + DD*YY_G ) / (DD + UN*RHO_F(IW))
         ENDDO SPECIES_LOOP
         
   END SELECT METHOD_OF_MASS_TRANSFER

   ! Only set species mass fraction in the ghost cell if it is not solid
    
   IF (IW<=NEWC .AND. N_SPECIES > 0 .AND. .NOT.SOLID(CELL_INDEX(II,JJ,KK)) .AND. .NOT.SOLID(CELL_INDEX(IIG,JJG,KKG))) &
      YYP(II,JJ,KK,1:N_SPECIES) = 2._EB*YY_F(IW,1:N_SPECIES) - YYP(IIG,JJG,KKG,1:N_SPECIES)

ENDDO WALL_CELL_LOOP

! Add gases from virtual particles

IF (VIRTUAL_PARTICLES .AND. CORRECTOR) THEN
   DROPLET_LOOP: DO I=1,NLP
      DR => DROPLET(I)
      PC => PARTICLE_CLASS(DR%CLASS)
      IBC = PC%SURF_INDEX
      IF (IBC<1) CYCLE DROPLET_LOOP
      SF => SURFACE(IBC)
      IW = DR%WALL_INDEX
      II = IJKW(1,IW)
      JJ = IJKW(2,IW)
      KK = IJKW(3,IW)
      IF (SF%SHRINK) THEN
         RADIUS = SUM(WALL(IW)%LAYER_THICKNESS)
      ELSEIF (SF%THERMALLY_THICK) THEN
         RADIUS = SF%THICKNESS
      ELSE
         RADIUS = SF%RADIUS
      ENDIF
      DR%R = RADIUS
      IF (RADIUS==0._EB) CYCLE DROPLET_LOOP

      SELECT CASE(SF%GEOMETRY)
         CASE(SURF_CARTESIAN)
            AREA = 2._EB*SF%LENGTH*SF%WIDTH
            AREA_SCALING = 1._EB
         CASE(SURF_CYLINDRICAL)
            AREA = TWOPI*RADIUS*SF%LENGTH
            AREA_SCALING = (SF%THICKNESS/RADIUS)
         CASE(SURF_SPHERICAL)
            AREA = 4._EB*PI*RADIUS**2
            AREA_SCALING = (SF%THICKNESS/RADIUS)**2
      END SELECT

      ! In PYROLYSIS, all the massfluxes were normalized by a virtual area
      ! based on the INITIAL radius. Here, correct the massflux using the CURRENT radius.

      MASSFLUX(IW,1:N_SPECIES) = MASSFLUX(IW,1:N_SPECIES)*AREA_SCALING
      MASSFLUX_ACTUAL(IW,1:N_SPECIES) = MASSFLUX_ACTUAL(IW,1:N_SPECIES)*AREA_SCALING

      RVC = RDX(II)*RDY(JJ)*RDZ(KK)
      ITMP = MIN(5000,NINT(TMP(II,JJ,KK)))
      CP = Y2CP_C(ITMP)
      H_G = CP*TMP(II,JJ,KK)
      DELTA_H_G = 0._EB
      DO N=1,N_SPECIES
         MW_RATIO = SPECIES(N)%RCON/RSUM(II,JJ,KK)
         M_DOT_PPP = MASSFLUX(IW,N)*AREA*RVC
         D_LAGRANGIAN(II,JJ,KK) =  D_LAGRANGIAN(II,JJ,KK) + DR%PWT*M_DOT_PPP*(MW_RATIO + DELTA_H_G/H_G)/RHO(II,JJ,KK)
         RHO_NEW = RHO(II,JJ,KK) + M_DOT_PPP*DT
         YYP(II,JJ,KK,N) = (RHO(II,JJ,KK)*YYP(II,JJ,KK,N) + M_DOT_PPP*DT)/RHO_NEW
         RHO(II,JJ,KK) = RHO_NEW
      ENDDO
      D_LAGRANGIAN(II,JJ,KK) =  D_LAGRANGIAN(II,JJ,KK) - QCONF(IW)*AREA*RVC/(RHO(II,JJ,KK)*H_G) * DR%PWT
   ENDDO DROPLET_LOOP
ENDIF

END SUBROUTINE SPECIES_BC


 
SUBROUTINE DENSITY_BC
 
! Compute density at wall from wall temperatures and mass fractions 

USE PHYSICAL_FUNCTIONS, ONLY : GET_SPECIFIC_GAS_CONSTANT
REAL(EB) :: YY_GET(1:N_SPECIES)
INTEGER  :: IBC,IIG,JJG,KKG,IW,II,JJ,KK
REAL(EB), POINTER, DIMENSION(:,:) :: PBAR_P=>NULL()
TYPE (SURFACE_TYPE), POINTER :: SF=>NULL()
REAL(EB), POINTER, DIMENSION(:,:,:) :: RHOP=>NULL()
 
IF (PREDICTOR) THEN 
   PBAR_P => PBAR_S
   RHOP => RHOS
ELSE 
   PBAR_P => PBAR
   RHOP => RHO
ENDIF

WALL_CELL_LOOP: DO IW=1,NWC
 
   IF (BOUNDARY_TYPE(IW)==NULL_BOUNDARY .OR. BOUNDARY_TYPE(IW)==POROUS_BOUNDARY) CYCLE WALL_CELL_LOOP

   II  = IJKW(1,IW)
   JJ  = IJKW(2,IW)
   KK  = IJKW(3,IW)
   IIG = IJKW(6,IW)
   JJG = IJKW(7,IW)
   KKG = IJKW(8,IW)
   IBC = IJKW(5,IW)
   IF (BOUNDARY_TYPE(IW)/=SOLID_BOUNDARY) THEN
      IF (BOUNDARY_TYPE(IW)==INTERPOLATED_BOUNDARY) IBC = INTERPOLATED_SURF_INDEX
      IF (BOUNDARY_TYPE(IW)==OPEN_BOUNDARY)         IBC = OPEN_SURF_INDEX
      IF (BOUNDARY_TYPE(IW)==MIRROR_BOUNDARY)       IBC = MIRROR_SURF_INDEX
   ENDIF
   SF => SURFACE(IBC)

   ! Determine ghost cell value of RSUM=R0*Sum(Y_i/M_i) 

   IF (N_SPECIES>0) THEN
      YY_GET = MAX(0._EB,YY_F(IW,:))
      CALL GET_SPECIFIC_GAS_CONSTANT(YY_GET,RSUM_F(IW))
   ENDIF
 
   ! Compute ghost cell density

   IF (BOUNDARY_TYPE(IW)/=INTERPOLATED_BOUNDARY) RHO_F(IW) = PBAR_P(KK,PRESSURE_ZONE_WALL(IW))/(RSUM_F(IW)*TMP_F(IW)) 

   ! Actually set the ghost cell value of density in the ghost cell if it is a solid wall

   IF ( (SOLID(CELL_INDEX(II,JJ,KK)) .AND. .NOT.SOLID(CELL_INDEX(IIG,JJG,KKG))) .OR.  &
         BOUNDARY_TYPE(IW)==OPEN_BOUNDARY .OR. BOUNDARY_TYPE(IW)==INTERPOLATED_BOUNDARY) THEN
      RHOP(II,JJ,KK) = 2._EB*RHO_F(IW) - RHOP(IIG,JJG,KKG)
      IF (N_SPECIES>0) RSUM(II,JJ,KK) = 2._EB*RSUM_F(IW) - RSUM(IIG,JJG,KKG)
   ENDIF
 
ENDDO WALL_CELL_LOOP

END SUBROUTINE DENSITY_BC
 


SUBROUTINE HVAC_BC

! Compute density at wall from wall temperatures and mass fractions 

USE PHYSICAL_FUNCTIONS, ONLY : GET_SPECIFIC_GAS_CONSTANT
REAL(EB) :: YY_GET(1:N_SPECIES),YY_G(1:N_SPECIES),UN,MFT,TMP_G,RHO_G,DD
REAL(EB) :: RHO_0,TMP_0,YY_0(1:N_SPECIES),UW_0,YY_ERR
INTEGER  :: IIG,JJG,KKG,IW,II,JJ,KK,N,IBC,COUNTER
LOGICAL :: ITER = .FALSE.
REAL(EB), POINTER, DIMENSION(:,:) :: PBAR_P=>NULL()
REAL(EB), POINTER, DIMENSION(:,:,:) :: RHOP=>NULL(),UU=>NULL(),VV=>NULL(),WW=>NULL()
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: YYP=>NULL()
TYPE (DUCTNODE_TYPE), POINTER :: DN=>NULL()
TYPE (DUCT_TYPE), POINTER :: DU=>NULL()
TYPE (SURFACE_TYPE), POINTER :: SF=>NULL()


IF (PREDICTOR) THEN 
   UU => U
   VV => V
   WW => W
   RHOP => RHOS
   IF (N_SPECIES > 0) YYP => YYS
   PBAR_P => PBAR_S
ELSE 
   UU => US
   VV => VS
   WW => WS
   RHOP => RHO
   IF (N_SPECIES > 0) YYP => YY
   PBAR_P => PBAR
ENDIF
WALL_CELL_LOOP: DO IW=1,NWC
   IBC = IJKW(5,IW)
   SF => SURFACE(IBC)   
   IF (SF%SPECIES_BC_INDEX/=HVAC_BOUNDARY .AND. SF%THERMAL_BC_INDEX/=HVAC_BOUNDARY) CYCLE WALL_CELL_LOOP
   II  = IJKW(1,IW)
   JJ  = IJKW(2,IW)
   KK  = IJKW(3,IW)
   IIG = IJKW(6,IW)
   JJG = IJKW(7,IW)
   KKG = IJKW(8,IW)
   COUNTER = 0

   ! Determine ghost cell value of RSUM=R0*Sum(Y_i/M_i) 

   IF (N_SPECIES>0) THEN
      YY_G = YYP(IIG,JJG,KKG,:)
      YY_GET = MAX(0._EB,YY_F(IW,:))
      CALL GET_SPECIFIC_GAS_CONSTANT(YY_GET,RSUM_F(IW))
   ENDIF
   TMP_G = TMP(IIG,JJG,KKG)
   IF (VENTS(VENT_INDEX(IW))%NODE_INDEX > 0) THEN
      DN=>DUCTNODE(VENTS(VENT_INDEX(IW))%NODE_INDEX)
      DU=>DUCT(DN%DUCT_INDEX(1))
      MFT = DN%DIR(1)*DU%VEL(NEW)*DU%RHO_D/VENTS(VENT_INDEX(IW))%FDS_AREA*DU%AREA
   ELSE
      MFT = 0._EB
   ENDIF
   RHO_G = RHOP(IIG,JJG,KKG)
   RHO_F(IW) = PBAR_P(KK,PRESSURE_ZONE_WALL(IW))/(RSUM_F(IW)*TMP_F(IW))  
   ITER = .TRUE.
   UN =  -MFT/RHO_F(IW)
   DO WHILE (ITER)
      ITER = .FALSE.
      RHO_0 = RHO_F(IW)
      TMP_0 = TMP_F(IW)
      YY_0  = YY_F(IW,:)
      UW_0 = -UN
      UN    = -MFT/RHO_F(IW)
      IF (PREDICTOR) UWS(IW) = -UN
      HEAT_TRANS_COEF(IW) = 0._EB
      QCONF(IW) = 0._EB      
      IF (UN > 0._EB .AND. VENTS(VENT_INDEX(IW))%NODE_INDEX > 0) THEN
         TMP_F(IW) = DU%TMP_D
      ELSE
         TMP_F(IW) = TMP(IIG,JJG,KKG)
      ENDIF

      IF (N_SPECIES==0) THEN
         MASSFLUX(IW,0) = MFT
      ELSE
            YY_ERR = 0._EB
         IF (UN < 0._EB) THEN
            YY_F(IW,:) = YY_G(:)
         ELSE
            IF (VENTS(VENT_INDEX(IW))%NODE_INDEX > 0) THEN
               MASSFLUX(IW,1:N_SPECIES) = -DN%YY(1:N_SPECIES)*MFT
               MASSFLUX(IW,0) = -MFT - SUM(MASSFLUX(IW,1:N_SPECIES))
            ELSE
               MASSFLUX(IW,:) = 0._EB
            ENDIF
            DO N=1,N_SPECIES
               DD = 2._EB*RHODW(IW,N)*RDN(IW)
               YY_F(IW,N) = ( MASSFLUX(IW,N) + DD*YY_G(N) ) / ( DD + UN*RHO_F(IW) )
               IF (YY_F(IW,N) > 1.E-10_EB) YY_ERR = MAX(YY_ERR,ABS(YY_F(IW,N)-YY_0(N))/YY_0(N))
            ENDDO
         ENDIF
         IF (COUNTER > 5) YY_F(IW,:) = 0.4_EB*YY_F(IW,:)+0.6_EB*YY_0(:)
         YY_GET = MAX(0._EB,YY_F(IW,:))
         CALL GET_SPECIFIC_GAS_CONSTANT(YY_GET,RSUM_F(IW))
         IF (YY_ERR > 1.E-6_EB) ITER = .TRUE.
      ENDIF
      RHO_F(IW) = PBAR_P(KK,PRESSURE_ZONE_WALL(IW))/(RSUM_F(IW)*TMP_F(IW))
      IF (ABS(RHO_0 - RHO_F(IW))/RHO_0 > 1.E-6_EB .OR. &
          ABS(TMP_0 - TMP_F(IW))/TMP_0 > 1.E-6_EB) ITER = .TRUE.
      IF (ABS(UW_0)>0._EB) THEN
          IF(ABS(UW_0 + UN)/ABS(UW_0) > 1.E-6_EB) ITER = .TRUE.
      ENDIF
      COUNTER = COUNTER + 1
      IF (COUNTER > 20) ITER = .FALSE.
   ENDDO

   ! Actually set the ghost cell value of density in the ghost cell if it is a solid wall

   IF ( (SOLID(CELL_INDEX(II,JJ,KK)) .AND. .NOT.SOLID(CELL_INDEX(IIG,JJG,KKG))) .OR.  &
         BOUNDARY_TYPE(IW)==OPEN_BOUNDARY .OR. BOUNDARY_TYPE(IW)==INTERPOLATED_BOUNDARY) THEN
      RHOP(II,JJ,KK) = 2._EB*RHO_F(IW) - RHOP(IIG,JJG,KKG)
      TMP(II,JJ,KK)  = 2._EB*TMP_F(IW) - TMP(IIG,JJG,KKG)
      IF (N_SPECIES>0) RSUM(II,JJ,KK)  = 2._EB*RSUM_F(IW) - RSUM(IIG,JJG,KKG)
      IF (N_SPECIES>0) YYP(II,JJ,KK,:) = 2._EB*YY_F(IW,:) - YYP(IIG,JJG,KKG,:)
   ENDIF
 
ENDDO WALL_CELL_LOOP

END SUBROUTINE HVAC_BC


 
SUBROUTINE PYROLYSIS(T,DT_BC)

! Loop through all the boundary cells that require a 1-D heat transfer calc

USE PHYSICAL_FUNCTIONS, ONLY: GET_MOLECULAR_WEIGHT
USE GEOMETRY_FUNCTIONS
USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
REAL(EB) :: DTMP,QNETF,QDXKF,QDXKB,RR,TMP_G,T,RFACF,RFACB,RFACF2,RFACB2,PPCLAUS,PPSURF,TMP_G_B,DT_BC, &
            DXKF,DXKB,REACTION_RATE,QRADINB,RFLUX_UP,RFLUX_DOWN,E_WALLB, &
            HVRG,Y_MF_G,Y_MF_W,RSUM_W, RSUM_G, MFLUX, MFLUX_S, VOLSUM,YPRSUM, &
            DXF, DXB,HTCB,Q_WATER_F,Q_WATER_B,TMP_F_OLD, DX_GRID, RHO_S0,DT2_BC,TOLERANCE,C_S_ADJUST_UNITS
!REAL(EB) :: YY_S,RSUM_S
INTEGER :: IBC,IIG,JJG,KKG,IIB,JJB,KKB,IWB,NWP,I,J,NR,NNN,NL,II,JJ,KK,IW,IOR,N,I_OBST,NS
REAL(EB) :: SMALLEST_CELL_SIZE(MAX_LAYERS),THICKNESS,YY_GET(1:N_SPECIES)
REAL(EB),ALLOCATABLE,DIMENSION(:) :: TMP_W_NEW
REAL(EB),ALLOCATABLE,DIMENSION(:,:) :: INT_WGT
INTEGER  :: N_LAYER_CELLS_NEW(MAX_LAYERS), NWP_NEW,I_GRAD,STEPCOUNT
LOGICAL :: POINT_SHRINK, RECOMPUTE,ITERATE,E_FOUND
TYPE (WALL_TYPE), POINTER :: WC=>NULL()
TYPE (SURFACE_TYPE), POINTER :: SF=>NULL()
TYPE (MATERIAL_TYPE), POINTER :: ML=>NULL()
 
! Special adjustment of specific heat for steady state applications

C_S_ADJUST_UNITS = 1000._EB/TIME_SHRINK_FACTOR
 
! Loop through the thermally-thick wall cells

WALL_CELL_LOOP: DO IW=1,NWC+NVWC

   IF (BOUNDARY_TYPE(IW)/=SOLID_BOUNDARY .AND. BOUNDARY_TYPE(IW)/=VIRTUAL_BOUNDARY) CYCLE WALL_CELL_LOOP
   IBC = IJKW(5,IW)
   SF  => SURFACE(IBC)
   IF (SF%THERMAL_BC_INDEX /= THERMALLY_THICK) CYCLE WALL_CELL_LOOP
   II  = IJKW(1,IW)
   JJ  = IJKW(2,IW)
   KK  = IJKW(3,IW)
   WC  => WALL(IW)
   IOR = IJKW(4,IW)
   IIG = IJKW(6,IW)
   JJG = IJKW(7,IW)
   KKG = IJKW(8,IW)
   IF (WC%BURNAWAY) THEN
      MASSFLUX(IW,:) = 0._EB
      CYCLE WALL_CELL_LOOP
   ENDIF

   SELECT CASE(SF%GEOMETRY)
   CASE(SURF_CARTESIAN)
      I_GRAD = 0
   CASE(SURF_CYLINDRICAL)
      I_GRAD = 1
   CASE(SURF_SPHERICAL)
      I_GRAD = 2
   END SELECT

   ! Compute convective heat flux at the surface
 
   TMP_G = TMP(IIG,JJG,KKG)
   IF (ASSUMED_GAS_TEMPERATURE > 0._EB) TMP_G = ASSUMED_GAS_TEMPERATURE  ! This is just for diagnostic calcs
   TMP_F_OLD = TMP_F(IW)
   DTMP = TMP_G - TMP_F_OLD
   HEAT_TRANS_COEF(IW) = HEAT_TRANSFER_COEFFICIENT(IW,IIG,JJG,KKG,IOR,TMP_G,DTMP,SF%H_FIXED)
   QCONF(IW) = HEAT_TRANS_COEF(IW)*DTMP
 
   ! Compute back side emissivity
  
   E_WALLB = SF%EMISSIVITY_BACK
   IF (E_WALLB < 0._EB .AND. SF%BACKING /= INSULATED) THEN
      E_WALLB = 0._EB
      VOLSUM = 0._EB
      IF (SF%SHRINK) THEN
         NWP = SUM(WC%N_LAYER_CELLS)
      ELSE
         NWP = SF%N_CELLS
      ENDIF
      DO N=1,SF%N_MATL
         IF (WC%RHO_S(NWP,N)==0._EB) CYCLE
         ML  => MATERIAL(SF%MATL_INDEX(N))
         VOLSUM = VOLSUM + WC%RHO_S(NWP,N)/ML%RHO_S
         E_WALLB = E_WALLB + WC%RHO_S(NWP,N)*ML%EMISSIVITY/ML%RHO_S
      ENDDO
      IF (VOLSUM > 0._EB) E_WALLB = E_WALLB/VOLSUM
   ENDIF

   ! Get heat losses from convection and radiation out of back of surface
  
   SELECT CASE(SF%BACKING)
      CASE(VOID)  ! Non-insulated backing to an ambient void
         DTMP = SF%TMP_BACK - TMP_B(IW)
         HTCB = HEAT_TRANSFER_COEFFICIENT(IW,-1,-1,-1,IOR,SF%TMP_BACK,DTMP,SF%H_FIXED)
         QRADINB   =  E_WALLB*SIGMA*SF%TMP_BACK**4
         Q_WATER_B = 0._EB
         TMP_G_B = SF%TMP_BACK
         
      CASE(INSULATED) 
         HTCB      = 0._EB
         QRADINB   = 0._EB
         E_WALLB   = 0._EB
         Q_WATER_B = 0._EB
         TMP_G_B   = TMPA
                  
      CASE(EXPOSED)  
         IWB = WALL_INDEX_BACK(IW)
         Q_WATER_B = 0._EB
         IF (BOUNDARY_TYPE(IWB)==SOLID_BOUNDARY) THEN
            IIB = IJKW(6,IWB)
            JJB = IJKW(7,IWB)
            KKB = IJKW(8,IWB)
            TMP_G_B  = TMP(IIB,JJB,KKB)
            DTMP = TMP_G_B - TMP_B(IW)
            HTCB = HEAT_TRANSFER_COEFFICIENT(IWB,IIB,JJB,KKB,IOR,TMP_G_B,DTMP,SF%H_FIXED)     
            HEAT_TRANS_COEF(IWB) = HTCB
            QRADINB  = QRADIN(IWB)
            IF (NLP>0) Q_WATER_B = -SUM(WCPUA(IWB,:))
         ELSE
            TMP_G_B  = TMPA
            DTMP = TMP_G_B - TMP_B(IW)
            HTCB = HEAT_TRANSFER_COEFFICIENT(IW,-1,-1,-1,IOR,TMP_G_B,DTMP,SF%H_FIXED)
            QRADINB  =  E_WALLB*SIGMA*TMPA4
         ENDIF
   END SELECT
 
   ! Take away energy flux due to water evaporation
 
   IF (NLP>0) THEN
      Q_WATER_F  = -SUM(WCPUA(IW,:))
   ELSE
      Q_WATER_F  = 0._EB
   ENDIF

   ! Compute grid for shrinking wall nodes

   COMPUTE_GRID: IF (SF%SHRINK) THEN
      NWP = SUM(WC%N_LAYER_CELLS)
      CALL GET_WALL_NODE_WEIGHTS(NWP,SF%N_LAYERS,WC%N_LAYER_CELLS,WC%LAYER_THICKNESS,SF%GEOMETRY, &
         WC%X_S(0:NWP),SF%LAYER_DIVIDE,DX_S(1:NWP),RDX_S(0:NWP+1),RDXN_S(0:NWP),DX_WGT_S(0:NWP),DXF,DXB,&
         LAYER_INDEX(0:NWP+1),MF_FRAC(1:NWP))
   ELSE COMPUTE_GRID
      NWP         = SF%N_CELLS
      DXF         = SF%DXF
      DXB         = SF%DXB
      DX_S(1:NWP)          = SF%DX(1:NWP)
      RDX_S(0:NWP+1)       = SF%RDX(0:NWP+1)
      RDXN_S(0:NWP)        = SF%RDXN(0:NWP)
      DX_WGT_S(0:NWP)      = SF%DX_WGT(0:NWP)
      LAYER_INDEX(0:NWP+1) = SF%LAYER_INDEX(0:NWP+1)
      MF_FRAC(1:NWP)       = SF%MF_FRAC(1:NWP)
   ENDIF COMPUTE_GRID

   ! Calculate reaction rates based on the solid phase reactions
 
   Q_S                  = 0._EB

   PYROLYSIS_MATERIAL_IF: IF (SF%PYROLYSIS_MODEL==PYROLYSIS_MATERIAL) THEN

   MFLUX                = MASSFLUX_ACTUAL(IW,I_FUEL)
   MASSFLUX(IW,:)       = 0._EB
   MASSFLUX_ACTUAL(IW,:)= 0._EB
   POINT_SHRINK         = .FALSE.
   WC%SHRINKING         = .FALSE.
   IF (SF%SHRINK) X_S_NEW(0:NWP)  = WC%X_S(0:NWP)
   
   POINT_LOOP1: DO I=1,NWP

      RHO_S0 = SF%LAYER_DENSITY(LAYER_INDEX(I))
      VOLSUM = 0._EB

      MATERIAL_LOOP1b: DO N=1,SF%N_MATL
         ML  => MATERIAL(SF%MATL_INDEX(N))

         IF (WC%RHO_S(I,N) <= 0._EB) CYCLE MATERIAL_LOOP1b
         IF (ML%PYROLYSIS_MODEL==PYROLYSIS_LIQUID) THEN
            VOLSUM = 1._EB
            CYCLE MATERIAL_LOOP1b
         ENDIF

         REACTION_LOOP: DO J=1,ML%N_REACTIONS
            ! Reaction rate in 1/s
            REACTION_RATE = ML%A(J)*(WC%RHO_S(I,N)/RHO_S0)**ML%N_S(J)*EXP(-ML%E(J)/(R0*WC%TMP_S(I)))
            ! power term
            DTMP = WC%TMP_S(I)-ML%TMP_THR(J)
            IF (ML%N_T(J)/=0._EB) THEN
               IF (DTMP > 0._EB) THEN
                  REACTION_RATE = REACTION_RATE * DTMP**ML%N_T(J)
               ELSE
                  REACTION_RATE = 0._EB
               ENDIF
            ELSE ! threshold
               IF (DTMP < 0._EB) REACTION_RATE = 0._EB
            ENDIF
            ! Reaction rate in kg/(m3s)
            REACTION_RATE = RHO_S0 * REACTION_RATE
            ! Limit reaction rate
            REACTION_RATE = MIN(REACTION_RATE , WC%RHO_S(I,N)/DT_BC)
            ! Compute mdot''_norm = mdot''' * r^(I_GRAD) * \Delta x / R^(I_GRAD)
            MFLUX_S = MF_FRAC(I)*REACTION_RATE/RDX_S(I)/SF%THICKNESS**I_GRAD
            ! Sum up local mass fluxes
            DO NS = 1,N_SPECIES
               MASSFLUX(IW,NS)        = MASSFLUX(IW,NS)        + ML%ADJUST_BURN_RATE(J,NS)*ML%NU_GAS(J,NS)*MFLUX_S
               MASSFLUX_ACTUAL(IW,NS) = MASSFLUX_ACTUAL(IW,NS) +                           ML%NU_GAS(J,NS)*MFLUX_S
            ENDDO
            Q_S(I) = Q_S(I) - REACTION_RATE * ML%H_R(J)
            WC%RHO_S(I,N) = WC%RHO_S(I,N) - DT_BC*REACTION_RATE
            WC%RHO_S(I,N) = MAX(0._EB, WC%RHO_S(I,N))            
            IF (ML%NU_RESIDUE(J) .GT. 0._EB ) THEN
               NNN = SF%RESIDUE_INDEX(N,J)
               WC%RHO_S(I,NNN) = WC%RHO_S(I,NNN) + ML%NU_RESIDUE(J)*DT_BC*REACTION_RATE
            ENDIF
         ENDDO REACTION_LOOP
         VOLSUM = VOLSUM + WC%RHO_S(I,N)/ML%RHO_S
      ENDDO MATERIAL_LOOP1b

      IF (SF%SHRINK) THEN
      POINT_SHRINK = .TRUE.
      MATERIAL_LOOP1a: DO N=1,SF%N_MATL
         IF (WC%RHO_S(I,N).EQ. 0._EB) CYCLE MATERIAL_LOOP1a
         ML  => MATERIAL(SF%MATL_INDEX(N))
         IF (ML%N_REACTIONS.EQ.0) THEN
            POINT_SHRINK = .FALSE.
            EXIT MATERIAL_LOOP1a
         ENDIF
      ENDDO MATERIAL_LOOP1a
      ENDIF

      ! In points that actuall shrink, increase the density to account for filled material

      VOLSUM = MIN(VOLSUM,1._EB)
      IF (POINT_SHRINK) THEN
         IF (VOLSUM<1.0_EB) THEN
            WC%SHRINKING=.TRUE.
            IF (VOLSUM>0.0_EB) THEN
               MATERIAL_LOOP1c: DO N=1,SF%N_MATL
                  WC%RHO_S(I,N) = WC%RHO_S(I,N)/VOLSUM
               ENDDO MATERIAL_LOOP1c
            ENDIF
         ENDIF
      ELSE
         ! In points that do not shrink, do not change the cell size.
         VOLSUM = 1._EB
      ENDIF
      ! Compute new co-ordinates
      IF (SF%SHRINK) X_S_NEW(I) = X_S_NEW(I-1)+(WC%X_S(I)-WC%X_S(I-1))*VOLSUM

   ENDDO POINT_LOOP1

   ! If the fuel or water massflux is non-zero, set the ignition time

   IF (TW(IW)>T) THEN
      IF (SUM(MASSFLUX(IW,:)) > 0._EB) TW(IW) = T
   ENDIF

   ! Special reactions: LIQUID
   ! Liquid evaporation can only take place on the surface (1st cell)

   POINT_SHRINK = .FALSE.
   IF (SF%SHRINK) THEN
      POINT_SHRINK = .TRUE.
      MATERIAL_LOOP2a: DO N=1,SF%N_MATL
         ML  => MATERIAL(SF%MATL_INDEX(N))
         IF (ML%PYROLYSIS_MODEL/=PYROLYSIS_LIQUID .AND. WC%RHO_S(1,N)>0._EB) POINT_SHRINK = .FALSE.
      ENDDO MATERIAL_LOOP2a
      THICKNESS = WC%LAYER_THICKNESS(LAYER_INDEX(1))
   ELSE
      THICKNESS = SF%LAYER_THICKNESS(LAYER_INDEX(1))
   ENDIF

   ! Estimate the previous value of liquid mass fluxes. The possibility of multiple liquids not taken into account.

   MFLUX = MAX(0._EB,MFLUX - MASSFLUX_ACTUAL(IW,I_FUEL)) ! Decrease 
   MATERIAL_LOOP2: DO N=1,SF%N_MATL
      ML  => MATERIAL(SF%MATL_INDEX(N))
      IF (ML%PYROLYSIS_MODEL/=PYROLYSIS_LIQUID) CYCLE MATERIAL_LOOP2
      IF (WC%RHO_S(1,N)==0._EB) CYCLE MATERIAL_LOOP2
      IF (ML%NU_FUEL(1)>0._EB) MFLUX = MFLUX/ML%NU_FUEL(1)
      ! gas phase 
      YY_GET = MAX(0._EB,YY(IIG,JJG,KKG,:))
      Y_MF_G = YY_GET(I_FUEL)
      CALL GET_MOLECULAR_WEIGHT(YY_GET,RSUM_G)
      ! wall values
      YY_GET = MAX(0._EB,YY_F(IW,:))
      Y_MF_W = YY_GET(I_FUEL)
      CALL GET_MOLECULAR_WEIGHT(YY_GET,RSUM_W)
      ! Weighted average of wall and gas values
      ! Alvernative 1
!      YPRSUM  = 0.2*(Y_MF_W/RSUM_W) + 0.8*(Y_MF_G/RSUM_G)
      YPRSUM  = 0.0*(Y_MF_W/RSUM_W) + 1.0*(Y_MF_G/RSUM_G)
      ! Alternative 2
!      YY_S    = 0.2*Y_MF_W+0.8*Y_MF_G
!      RSUM_S  = 0.2*RSUM_W+0.8*RSUM_G
!      YPRSUM  = YY_S/RSUM_S
      HVRG    = REACTION(1)%MW_FUEL*ML%H_R(1)/R0
      PPSURF  = MIN(1._EB,(R0/REACTION(1)%MW_FUEL)*YPRSUM)
      PPCLAUS = MIN(1._EB,EXP(HVRG*(1./ML%TMP_BOIL-1./WC%TMP_S(1))))
      ! Make initial guess
      IF ((MFLUX==0._EB) .AND. (PPSURF<PPCLAUS)) THEN         
         MFLUX = (ML%INIT_VAPOR_FLUX/(R0*TMPA/P_INF))*REACTION(1)%MW_FUEL/ML%ADJUST_BURN_RATE(1,I_FUEL)
      ENDIF
      ! Adjust MFLUX to reach equilibrium vapor pressure
      IF (PPSURF/=PPCLAUS .AND. PPSURF>0._EB) THEN
         MFLUX = MFLUX*MIN(1.02_EB,MAX(0.98_EB,PPCLAUS/PPSURF))
      ENDIF
      IF (MFLUX > 0._EB .AND. TW(IW)>T) TW(IW) = T
      IF (WC%TMP_S(1)>ML%TMP_BOIL) THEN
         ! Net flux guess for liquid evap
         QNETF = Q_WATER_F + QRADIN(IW) - QRADOUT(IW) + QCONF(IW)
         MFLUX = MAX(MFLUX,1.02*(QNETF - 2.*(ML%TMP_BOIL-WC%TMP_S(1))/DXF/K_S(1))/ML%H_R(1))
      ENDIF
      MFLUX = MIN(MFLUX,THICKNESS*ML%RHO_S/DT_BC)
      ! CYLINDRICAL and SPHERICAL scaling not implemented
      DO NS = 1,N_SPECIES
         MASSFLUX(IW,NS)        = MASSFLUX(IW,NS)        + ML%ADJUST_BURN_RATE(1,NS)*ML%NU_GAS(1,NS)*MFLUX
         MASSFLUX_ACTUAL(IW,NS) = MASSFLUX_ACTUAL(IW,NS) +                           ML%NU_GAS(1,NS)*MFLUX
      ENDDO
      Q_S(1) = Q_S(1) - MFLUX*ML%H_R(1)/DX_S(1)  ! no improvement (in cone test) if used updated RDX 

      DX_GRID = DT_BC*MFLUX/ML%RHO_S
      IF (POINT_SHRINK) THEN
         X_S_NEW(1:NWP) = MAX(0._EB,X_S_NEW(1:NWP)-DX_GRID)
         IF (DX_GRID > 0._EB) WC%SHRINKING = .TRUE.
      ENDIF
      
      EXIT MATERIAL_LOOP2     ! Can handle only one LIQUID fuel at the time

   ENDDO MATERIAL_LOOP2

   ! Re-generate grid for shrinking wall

   N_LAYER_CELLS_NEW = 0
   RECOMPUTE_GRID: IF (WC%SHRINKING) THEN
      NWP_NEW = 0
      THICKNESS = 0._EB
      RECOMPUTE = .FALSE.
      I = 0
      LAYER_LOOP: DO NL=1,SF%N_LAYERS
         WC%LAYER_THICKNESS(NL) = X_S_NEW(I+WC%N_LAYER_CELLS(NL))-X_S_NEW(I)
         ! Remove very thin layers
         IF (WC%LAYER_THICKNESS(NL) < 1.0E-6_EB) THEN 
            X_S_NEW(I+WC%N_LAYER_CELLS(NL):NWP) = X_S_NEW(I+WC%N_LAYER_CELLS(NL):NWP)-WC%LAYER_THICKNESS(NL)
            WC%LAYER_THICKNESS(NL) = 0._EB
            N_LAYER_CELLS_NEW(NL)  = 0
         ELSE
            CALL GET_N_LAYER_CELLS(SF%MIN_DIFFUSIVITY(NL),WC%LAYER_THICKNESS(NL), &
               SF%STRETCH_FACTOR,SF%CELL_SIZE_FACTOR,N_LAYER_CELLS_NEW(NL),SMALLEST_CELL_SIZE(NL))
            NWP_NEW = NWP_NEW + N_LAYER_CELLS_NEW(NL)
         ENDIF
         IF ( N_LAYER_CELLS_NEW(NL) .NE. WC%N_LAYER_CELLS(NL)) RECOMPUTE = .TRUE.
         THICKNESS = THICKNESS + WC%LAYER_THICKNESS(NL)
         I = I + WC%N_LAYER_CELLS(NL)
      ENDDO LAYER_LOOP
      DO I = 1,NWP
         IF ( (X_S_NEW(I)-X_S_NEW(I-1)) < SF%REGRID_FACTOR*SMALLEST_CELL_SIZE(LAYER_INDEX(I))) RECOMPUTE = .TRUE.
      ENDDO
      ! Shrinking wall has gone to zero thickness.
      IF (THICKNESS == 0._EB) THEN
         WC%TMP_S(0:NWP+1)    = MAX(TMPMIN,SF%TMP_BACK)
         TMP_F(IW)            = MIN(TMPMAX,MAX(TMPMIN,SF%TMP_BACK))
         TMP_B(IW)            = MIN(TMPMAX,MAX(TMPMIN,SF%TMP_BACK))
         QCONF(IW)            = 0._EB
         MASSFLUX(IW,I_FUEL)  = 0._EB
         MASSFLUX_ACTUAL(IW,I_FUEL) = 0._EB
         IF (I_WATER>0) MASSFLUX(IW,I_WATER)  = 0._EB
         IF (I_WATER>0) MASSFLUX_ACTUAL(IW,I_WATER)  = 0._EB
         WC%N_LAYER_CELLS     = 0
         WC%BURNAWAY          = .TRUE.
         I_OBST               = OBST_INDEX_W(IW)
         IF (I_OBST > 0) THEN
            IF (OBSTRUCTION(I_OBST)%CONSUMABLE) OBSTRUCTION(I_OBST)%MASS = -1.
         ENDIF
         CYCLE WALL_CELL_LOOP
      ENDIF

      WC%X_S(0:NWP) = X_S_NEW(0:NWP)
      X_S_NEW = 0._EB
      IF (RECOMPUTE) THEN
         CALL GET_WALL_NODE_COORDINATES(NWP_NEW,SF%N_LAYERS,N_LAYER_CELLS_NEW, &
            SMALLEST_CELL_SIZE(1:SF%N_LAYERS),SF%STRETCH_FACTOR,X_S_NEW(0:NWP_NEW))
         CALL GET_WALL_NODE_WEIGHTS(NWP_NEW,SF%N_LAYERS,N_LAYER_CELLS_NEW,WC%LAYER_THICKNESS,SF%GEOMETRY, &
            X_S_NEW(0:NWP_NEW),SF%LAYER_DIVIDE,DX_S(1:NWP_NEW),RDX_S(0:NWP_NEW+1),RDXN_S(0:NWP_NEW),&
            DX_WGT_S(0:NWP_NEW),DXF,DXB,LAYER_INDEX(0:NWP_NEW+1),MF_FRAC(1:NWP_NEW))           
         ! Interpolate densities and temperature from old grid to new grid
         ALLOCATE(INT_WGT(NWP_NEW,NWP))
         CALL GET_INTERPOLATION_WEIGHTS(SF%N_LAYERS,NWP,NWP_NEW,WC%N_LAYER_CELLS,N_LAYER_CELLS_NEW, &
                                    WC%X_S(0:NWP),X_S_NEW(0:NWP_NEW),INT_WGT)      
         CALL INTERPOLATE_WALL_ARRAY(NWP,NWP_NEW,INT_WGT,WC%TMP_S(1:NWP))
         WC%TMP_S(NWP_NEW+1) = WC%TMP_S(NWP+1)
         CALL INTERPOLATE_WALL_ARRAY(NWP,NWP_NEW,INT_WGT,Q_S(1:NWP))
         DO N=1,SF%N_MATL
            ML  => MATERIAL(SF%MATL_INDEX(N))
            CALL INTERPOLATE_WALL_ARRAY(NWP,NWP_NEW,INT_WGT,WC%RHO_S(1:NWP,N))
         ENDDO
         DEALLOCATE(INT_WGT)
         WC%N_LAYER_CELLS = N_LAYER_CELLS_NEW(1:SF%N_LAYERS)
         NWP = NWP_NEW
         WC%X_S(0:NWP) = X_S_NEW(0:NWP)      ! Note: WC%X_S(NWP+1...) are not set to zero.
      ELSE      
         CALL GET_WALL_NODE_WEIGHTS(NWP,SF%N_LAYERS,N_LAYER_CELLS_NEW,WC%LAYER_THICKNESS,SF%GEOMETRY, &
            WC%X_S(0:NWP),SF%LAYER_DIVIDE,DX_S(1:NWP),RDX_S(0:NWP+1),RDXN_S(0:NWP),DX_WGT_S(0:NWP),DXF,DXB, &
            LAYER_INDEX(0:NWP+1),MF_FRAC(1:NWP))
      ENDIF
      
   ENDIF RECOMPUTE_GRID

   ELSEIF (SF%PYROLYSIS_MODEL==PYROLYSIS_SPECIFIED) THEN

      ! Take off energy corresponding to specified burning rate
      Q_S(1) = Q_S(1) - MASSFLUX(IW,I_FUEL)*SF%H_V/DX_S(1)

   ENDIF PYROLYSIS_MATERIAL_IF

   ! Calculate thermal properties 
 
   K_S     = 0._EB
   RHO_S   = 0._EB
   RHOCBAR = 0._EB
   E_WALL(IW) = 0._EB
   E_FOUND = .FALSE.
   POINT_LOOP3: DO I=1,NWP
      VOLSUM = 0._EB
      MATERIAL_LOOP3: DO N=1,SF%N_MATL
         IF (WC%RHO_S(I,N)==0._EB) CYCLE MATERIAL_LOOP3
         ML  => MATERIAL(SF%MATL_INDEX(N))
         VOLSUM = VOLSUM + WC%RHO_S(I,N)/ML%RHO_S
         IF (ML%K_S>0._EB) THEN  
            K_S(I) = K_S(I) + WC%RHO_S(I,N)*ML%K_S/ML%RHO_S
         ELSE
            NR     = -NINT(ML%K_S)
            K_S(I) = K_S(I) + WC%RHO_S(I,N)*EVALUATE_RAMP(WC%TMP_S(I),0._EB,NR)/ML%RHO_S
         ENDIF

         IF (ML%C_S>0._EB) THEN
            RHOCBAR(I) = RHOCBAR(I) + WC%RHO_S(I,N)*ML%C_S
         ELSE
            NR     = -NINT(ML%C_S)
            RHOCBAR(I) = RHOCBAR(I) + WC%RHO_S(I,N)*EVALUATE_RAMP(WC%TMP_S(I),0._EB,NR)*C_S_ADJUST_UNITS
         ENDIF
         IF (.NOT.E_FOUND) E_WALL(IW) = E_WALL(IW) + WC%RHO_S(I,N)*ML%EMISSIVITY/ML%RHO_S
         RHO_S(I)   = RHO_S(I) + WC%RHO_S(I,N)
         
      ENDDO MATERIAL_LOOP3
      
      IF (VOLSUM > 0._EB) THEN
         K_S(I) = K_S(I)/VOLSUM
         IF (.NOT.E_FOUND) E_WALL(IW) = E_WALL(IW)/VOLSUM
      ENDIF
      IF (E_WALL(IW)>0._EB) E_FOUND = .TRUE.

      IF (K_S(I)==0._EB)      K_S(I)      = 10000._EB
      IF (RHOCBAR(I)==0._EB)  RHOCBAR(I)  = 0.001_EB

   ENDDO POINT_LOOP3
  
   K_S(0)     = K_S(1)     ! Calculate average K_S between at grid cell boundaries. Store result in K_S
   K_S(NWP+1) = K_S(NWP)
   DO I=1,NWP-1 
      K_S(I)  = 1.0_EB / ( DX_WGT_S(I)/K_S(I) + (1.-DX_WGT_S(I))/K_S(I+1) )
   ENDDO

   ! Calculate internal radiation
      
   IF (SF%INTERNAL_RADIATION) THEN

      KAPPA_S = 0._EB
      DO I=1,NWP
         VOLSUM = 0._EB
         DO N=1,SF%N_MATL
            IF (WC%RHO_S(I,N)==0._EB) CYCLE
            ML  => MATERIAL(SF%MATL_INDEX(N))
            VOLSUM = VOLSUM + WC%RHO_S(I,N)/ML%RHO_S
            KAPPA_S(I) = KAPPA_S(I) + WC%RHO_S(I,N)*ML%KAPPA_S/ML%RHO_S
         ENDDO
     !!  KAPPA_S(I) = 2.*KAPPA_S(I)*DX_S(I)/VOLSUM    ! kappa = 2*dx*kappa
   
         IF (VOLSUM>0._EB) KAPPA_S(I) = 2.*KAPPA_S(I)/(RDX_S(I)*VOLSUM)    ! kappa = 2*dx*kappa or 2*r*dr*kappa
      ENDDO
      DO I=0,NWP
         IF (SF%GEOMETRY==SURF_CYLINDRICAL) THEN
            R_S(I) = SF%THICKNESS-SF%X_S(I)
         ELSE
            R_S(I) = 1._EB
         ENDIF
      ENDDO
      ! solution inwards
      RFLUX_UP = QRADIN(IW) + (1.-E_WALL(IW))*QRADOUT(IW)/(E_WALL(IW)+1.0E-10_EB)
      DO I=1,NWP
      !! RFLUX_DOWN =  ( RFLUX_UP + KAPPA_S(I)*SIGMA*WC%TMP_S(I)**4 ) / (1. + KAPPA_S(I))
         RFLUX_DOWN =  ( R_S(I-1)*RFLUX_UP + KAPPA_S(I)*SIGMA*WC%TMP_S(I)**4 ) / (R_S(I) + KAPPA_S(I))
      !! Q_S(I) = Q_S(I) + (RFLUX_UP - RFLUX_DOWN)/DX_S(I)
         Q_S(I) = Q_S(I) + (R_S(I-1)*RFLUX_UP - R_S(I)*RFLUX_DOWN)*RDX_S(I)
         RFLUX_UP = RFLUX_DOWN
      ENDDO
      IF (SF%BACKING==EXPOSED) THEN
         IF (BOUNDARY_TYPE(IWB)==SOLID_BOUNDARY) QRADOUT(IWB) = E_WALLB*RFLUX_UP
      ENDIF
      ! solution outwards
      RFLUX_UP = QRADINB + (1.-E_WALLB)*RFLUX_UP
      DO I=NWP,1,-1
      !! RFLUX_DOWN =  ( RFLUX_UP + KAPPA_S(I)*SIGMA*WC%TMP_S(I)**4 ) / (1. + KAPPA_S(I))
         RFLUX_DOWN =  ( R_S(I)*RFLUX_UP + KAPPA_S(I)*SIGMA*WC%TMP_S(I)**4 ) / (R_S(I-1) + KAPPA_S(I))
      !! Q_S(I) = Q_S(I) + (RFLUX_UP - RFLUX_DOWN)/DX_S(I)
         Q_S(I) = Q_S(I) + (R_S(I)*RFLUX_UP - R_S(I-1)*RFLUX_DOWN)*RDX_S(I)
         RFLUX_UP = RFLUX_DOWN
      ENDDO
      QRADOUT(IW) = E_WALL(IW)*RFLUX_DOWN
   ENDIF

   ! Update the 1-D heat transfer equation 

   DT2_BC = DT_BC
   STEPCOUNT = 1
   ALLOCATE(TMP_W_NEW(0:NWP+1))
   TMP_W_NEW = WC%TMP_S(0:NWP+1)
   WALL_ITERATE: DO 
      ITERATE=.FALSE.
      SUB_TIME: DO N=1,STEPCOUNT
         DXKF   = K_S(0)/DXF
         DXKB   = K_S(NWP)/DXB

         DO I=1,NWP
            BBS(I) = -0.5_EB*DT2_BC*K_S(I-1)*RDXN_S(I-1)*RDX_S(I)/RHOCBAR(I) ! DT_BC->DT2_BC
            AAS(I) = -0.5_EB*DT2_BC*K_S(I)  *RDXN_S(I)  *RDX_S(I)/RHOCBAR(I)
         ENDDO
         DDS(1:NWP) = 1._EB    - AAS(1:NWP) - BBS(1:NWP)
         DO I=1,NWP
            CCS(I) = TMP_W_NEW(I) - AAS(I)*(TMP_W_NEW(I+1)-TMP_W_NEW(I)) + BBS(I)*(TMP_W_NEW(I)-TMP_W_NEW(I-1)) &
                     + DT2_BC*Q_S(I)/RHOCBAR(I)
         ENDDO
         IF (.NOT. RADIATION .OR. SF%INTERNAL_RADIATION) THEN
            RFACF = 0.25_EB*HEAT_TRANS_COEF(IW)
            RFACB = 0.25_EB*HTCB
         ELSE
            RFACF = 0.25_EB*HEAT_TRANS_COEF(IW)+2._EB*E_WALL(IW)*SIGMA*TMP_F(IW)**3
            RFACB = 0.25_EB*HTCB               +2._EB*E_WALLB*   SIGMA*TMP_B(IW)**3
         ENDIF
         RFACF2 = (DXKF-RFACF)/(DXKF+RFACF)
         RFACB2 = (DXKB-RFACB)/(DXKB+RFACB)
         IF (.NOT. RADIATION .OR. SF%INTERNAL_RADIATION) THEN
            QDXKF = (HEAT_TRANS_COEF(IW)*(TMP_G   - 0.5_EB*TMP_F(IW)) + Q_WATER_F)/(DXKF+RFACF)
            QDXKB = (HTCB*               (TMP_G_B - 0.5_EB*TMP_B(IW)) + Q_WATER_B)/(DXKB+RFACB)
         ELSE
            QDXKF = (HEAT_TRANS_COEF(IW)*(TMP_G   - 0.5_EB*TMP_F(IW)) + QRADIN(IW) + 3.*E_WALL(IW)*SIGMA*TMP_F(IW)**4 + Q_WATER_F) &
                  /(DXKF+RFACF)
            QDXKB = (HTCB*               (TMP_G_B - 0.5_EB*TMP_B(IW)) + QRADINB    + 3.*E_WALLB   *SIGMA*TMP_B(IW)**4 + Q_WATER_B) &
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
            TOLERANCE = MAXVAL(ABS((TMP_W_NEW-WC%TMP_S(0:NWP+1))/WC%TMP_S(0:NWP+1)), &
               WC%TMP_S(0:NWP+1)>0._EB) ! returns a negative number, if all TMP_S == 0.
            IF (TOLERANCE<0.0_EB) &
            TOLERANCE = MAXVAL(ABS((TMP_W_NEW-WC%TMP_S(0:NWP+1))/TMP_W_NEW), &
               TMP_W_NEW>0._EB) 
            IF (TOLERANCE > 0.2_EB) THEN
               STEPCOUNT = MIN(200,STEPCOUNT * (INT(TOLERANCE/0.2_EB) + 1))
               ITERATE=.TRUE.
               DT2_BC=DT_BC/REAL(STEPCOUNT)
               TMP_W_NEW = WC%TMP_S(0:NWP+1)
            ENDIF
         ENDIF
         TMP_F(IW)       = 0.5_EB*(TMP_W_NEW(0)+TMP_W_NEW(1))
         TMP_F(IW)       = MIN(TMPMAX,MAX(TMPMIN,TMP_F(IW)))
         TMP_B(IW)       = 0.5_EB*(TMP_W_NEW(NWP)+TMP_W_NEW(NWP+1))
         TMP_B(IW)       = MIN(TMPMAX,MAX(TMPMIN,TMP_B(IW)))
      ENDDO SUB_TIME
      IF (.NOT. ITERATE) EXIT WALL_ITERATE
   ENDDO WALL_ITERATE
   WC%TMP_S(0:NWP+1) = TMP_W_NEW
   DEALLOCATE(TMP_W_NEW)

   ! If the surface temperature exceeds the ignition temperature, burn it

   IF (TW(IW) > T ) THEN
      IF (TMP_F(IW) >= SF%TMP_IGN) TW(IW) = T
   ENDIF
 
   ! Determine convective heat flux at the wall
 
   QCONF(IW) = HEAT_TRANS_COEF(IW) * (TMP_G - 0.5_EB * (TMP_F(IW) + TMP_F_OLD) )

ENDDO WALL_CELL_LOOP
 
END SUBROUTINE PYROLYSIS
 
 

REAL(EB) FUNCTION HEAT_TRANSFER_COEFFICIENT(IW,IIG,JJG,KKG,IOR,TMP_G,DELTA_TMP,H_FIXED)

USE PHYSICAL_FUNCTIONS, ONLY: GET_SPECIFIC_HEAT,GET_VISCOSITY
USE TURBULENCE, ONLY: WERNER_WENGLE_WALL_MODEL, SURFACE_HEAT_FLUX_MODEL
INTEGER  :: IW,IIG,JJG,KKG,IOR,ITMP,IHTMODEL
INTEGER, PARAMETER :: IHEDDY=1,IHLOGLAW=2,ICHILTON=3
REAL(EB) :: TMP_G,DELTA_TMP,U2=0._EB,V2=0._EB,W2=0._EB,VELCON,H_NATURAL,H_FORCED,H_FIXED,CP,YY_GET(1:N_SPECIES),MU_G
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU=>NULL(),VV=>NULL(),WW=>NULL(),RHOP=>NULL()
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: YYP=>NULL()
TYPE(SURFACE_TYPE), POINTER :: SF
REAL(EB) :: DUMMY,NU,DN

! If convective heat transfer is turned off, just return

IF (.NOT.CONVECTION) THEN
   HEAT_TRANSFER_COEFFICIENT = 0._EB
   RETURN
ENDIF

! If the user wants a specified HTC, set it and return

IF (H_FIXED >= 0._EB) THEN
   HEAT_TRANSFER_COEFFICIENT = H_FIXED
   RETURN
ENDIF

! If this is only a DNS calculation, set HTC and return

IF (DNS) THEN
   HEAT_TRANSFER_COEFFICIENT = 2._EB*KW(IW)*RDN(IW)
   RETURN
ENDIF

! Calculate the HTC for natural convection

SELECT CASE(ABS(IOR))
   CASE(0:2)
      H_NATURAL = HCV*ABS(DELTA_TMP)**ONTH
   CASE(3)
      H_NATURAL = HCH*ABS(DELTA_TMP)**ONTH
END SELECT

! Calculate the HTC for forced convection

H_FORCED = 0._EB

NO_GAS_CELL_IF: IF (IIG>=0) THEN  

   IF (PREDICTOR) THEN
      UU => U
      VV => V
      WW => W
      RHOP => RHOS
      YYP => YYS
   ELSE
      UU => US
      VV => VS
      WW => WS
      RHOP => RHO
      YYP => YY  
   ENDIF

   ! Calculate tangential velocity near the surface

   U2 = 0.25_EB*(UU(IIG,JJG,KKG)+UU(IIG-1,JJG,KKG))**2
   V2 = 0.25_EB*(VV(IIG,JJG,KKG)+VV(IIG,JJG-1,KKG))**2
   W2 = 0.25_EB*(WW(IIG,JJG,KKG)+WW(IIG,JJG,KKG-1))**2
   SELECT CASE(ABS(IOR))
      CASE(1)
         U2 = 0._EB
      CASE(2)
         V2 = 0._EB
      CASE(3)
         W2 = 0._EB
   END SELECT 

   VELCON = (U2+V2+W2)**0.4_EB
   
   ! Extra parameters for experimental models

   IHTMODEL = 0
   IF (H_EDDY .OR. H_LOGLAW .OR. H_CHILTON_COLBURN) THEN
      ITMP = MIN(5000,NINT(TMP_G))
      H_NATURAL = 0._EB
      IF (H_EDDY)            IHTMODEL = IHEDDY
      IF (H_LOGLAW)          IHTMODEL = IHLOGLAW
      IF (H_CHILTON_COLBURN) IHTMODEL = ICHILTON
      SF => SURFACE(IJKW(5,IW))
      IF (N_SPECIES > 0) THEN
         YY_GET = MAX(0._EB,YYP(IIG,JJG,KKG,:))
         CALL GET_VISCOSITY(YY_GET,MU_G,ITMP)
         NU = MU_G/RHOP(IIG,JJG,KKG)
      ELSE
         MU_G = Y2MU_C(ITMP)*SPECIES(0)%MW
         NU = MU_G/RHOP(IIG,JJG,KKG)
      ENDIF
      DN = 1._EB/RDN(IW)
      IF (N_SPECIES > 0) THEN
         CALL GET_SPECIFIC_HEAT(YY_GET,CP,ITMP)
      ELSE
         CP = Y2CP_C(ITMP)
      ENDIF
   ENDIF

   ! Calculate the HTC for forced convection 

   SELECT CASE (IHTMODEL)
      CASE (IHEDDY)
         H_FORCED = MU(IIG,JJG,KKG)*RPR*CP*2._EB*RDN(IW)
      CASE (IHLOGLAW)
         CALL WERNER_WENGLE_WALL_MODEL(DUMMY,U_TAU(IW),SQRT(U2+V2+W2),NU,DN,SF%ROUGHNESS)
         CALL SURFACE_HEAT_FLUX_MODEL(H_FORCED,U_TAU(IW),DN,SF%ROUGHNESS,IOR,RHOP(IIG,JJG,KKG),CP)
      CASE (ICHILTON)
         VELCON = SQRT(U2+V2+W2)
         CALL WERNER_WENGLE_WALL_MODEL(DUMMY,U_TAU(IW),VELCON,NU,DN,SF%ROUGHNESS)
         IF (VELCON > 0._EB) H_FORCED = U_TAU(IW)**2/VELCON*RHO_F(IW)*CP*RPR**TWTH           
      CASE DEFAULT
         H_FORCED  = C_FORCED*VELCON*RHOP(IIG,JJG,KKG)**0.8_EB
   END SELECT
   
ENDIF NO_GAS_CELL_IF

HEAT_TRANSFER_COEFFICIENT = MAX(H_FORCED,H_NATURAL)

END FUNCTION HEAT_TRANSFER_COEFFICIENT



SUBROUTINE GET_REV_wall(MODULE_REV,MODULE_DATE)
INTEGER,INTENT(INOUT) :: MODULE_REV
CHARACTER(255),INTENT(INOUT) :: MODULE_DATE

WRITE(MODULE_DATE,'(A)') wallrev(INDEX(wallrev,':')+1:LEN_TRIM(wallrev)-2)
READ (MODULE_DATE,'(I5)') MODULE_REV
WRITE(MODULE_DATE,'(A)') walldate

END SUBROUTINE GET_REV_wall

 

SUBROUTINE CALC_SOOT_DEPOSITION(NM)

USE PHYSICAL_FUNCTIONS, ONLY: GET_VISCOSITY,GET_CONDUCTIVITY
USE GLOBAL_CONSTANTS, ONLY: EVACUATION_ONLY,SOLID_PHASE_ONLY,SOLID_BOUNDARY,I_PROG_SOOT,N_SPECIES,TUSED
INTEGER, INTENT(IN) :: NM
!REAL(EB), PARAMETER :: CS=1.147_EB,CT=2.18_EB,CM=1.146_EB,KPKG=1._EB,PARTD=2.E-8_EB
REAL(EB) :: U_THERM,U_TURB,TGAS,TWALL,MUGAS,Y_SOOT,RHOG,YIN(1:N_SPECIES),YDEP,K_AIR
INTEGER  :: IIG,JJG,KKG,IW,IOR,ITMP

!REAL(EB), PARAMETER :: A=8.3_EB,B=1._EB/7._EB,Z_PLUS_TURBULENT = 11.81_EB,ALPHA=7.202125273562269_EB 
!REAL(EB), PARAMETER :: BETA=1._EB+B,ETA=(1._EB+B)/A,GAMMA=2._EB/(1._EB+B)

 
IF (EVACUATION_ONLY(NM)) RETURN
IF (SOLID_PHASE_ONLY) RETURN

CALL POINT_TO_MESH(NM)
U_THERM=0._EB 
U_TURB=0._EB 
WALL_CELL_LOOP: DO IW=1,NWC+NVWC
   IF (BOUNDARY_TYPE(IW)/=SOLID_BOUNDARY .OR. UW(IW) < 0._EB) CYCLE WALL_CELL_LOOP
   IOR = IJKW(4,IW)
   IIG = IJKW(6,IW)
   JJG = IJKW(7,IW)
   KKG = IJKW(8,IW)
   YIN = MAX(0._EB,YY(IIG,JJG,KKG,:))
   IF (YIN(I_PROG_SOOT) < 1.E-14_EB) CYCLE WALL_CELL_LOOP
   TGAS = TMP(IIG,JJG,KKG)
   TWALL = TMP_F(IW)
   ITMP = MIN(NINT(0.5_EB*(TGAS+TWALL)),5000)
   RHOG=RHO(IIG,JJG,KKG)
   CALL GET_VISCOSITY(YIN,MUGAS,ITMP)
   CALL GET_CONDUCTIVITY(YIN,K_AIR,ITMP)
   IF (THERMOPHORETIC_DEPOSITION) U_THERM = 0.55_EB*HEAT_TRANS_COEF(IW)*(TGAS-TWALL)*MUGAS/(TGAS*RHOG*K_AIR)
   IF (TURBULENT_DEPOSITION) U_TURB = 0.037_EB*U_TAU(IW)
   IF (U_THERM+U_TURB < 0._EB) CYCLE WALL_CELL_LOOP   
   YIN = YIN * RHOG  
   Y_SOOT = YIN(I_PROG_SOOT)   
   YDEP =Y_SOOT*MIN(1._EB,(U_THERM+U_TURB)*DT*RDN(IW))
   YIN(I_PROG_SOOT) = Y_SOOT - YDEP      
   AWMSOOT(IW)=AWMSOOT(IW)+YDEP/RDN(IW)
   RHO(IIG,JJG,KKG) = RHOG - YDEP
   YY(IIG,JJG,KKG,:) = YIN / RHO(IIG,JJG,KKG)
   
ENDDO WALL_CELL_LOOP

END SUBROUTINE CALC_SOOT_DEPOSITION


END MODULE WALL_ROUTINES
