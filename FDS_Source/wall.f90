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
USE COMP_FUNCTIONS, ONLY: SECOND
REAL(EB) :: TNOW
REAL(EB), INTENT(IN) :: T
INTEGER, INTENT(IN) :: NM

IF (EVACUATION_ONLY(NM)) RETURN

TNOW=SECOND()

CALL POINT_TO_MESH(NM)

CALL THERMAL_BC(T)
CALL SPECIES_BC(T)
CALL DENSITY_BC

TUSED(6,NM)=TUSED(6,NM)+SECOND()-TNOW


CONTAINS
 

SUBROUTINE THERMAL_BC(T)

! Thermal boundary conditions for adiabatic, fixed temperature, fixed flux and interpolated boundaries.
! One dimensional heat transfer and pyrolysis is done in PYROLYSIS, which is called at the end of this routine.

USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP 
USE PHYSICAL_FUNCTIONS, ONLY: GET_AVERAGE_SPECIFIC_HEAT
REAL(EB) :: DT_BC,T,TSI,TMP_G,RHO_G,DTMP,TMP_OTHER,CP_TERM,RHOWAL,RAMP_FACTOR,QNET,FDERIV,TMP_EXTERIOR,UN,ARO,UWO,&
            CP,CP_F,YY_N(1:N_SPECIES),QEXTRA,CP_TERM_F
INTEGER  :: IOR,II,JJ,KK,IBC,IIG,JJG,KKG,IW,NOM,IIO,JJO,KKO,N_INT_CELLS,ITMP,N
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW,RHOP,OM_RHOP
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: YYP
LOGICAL :: INFLOW
TYPE (SURFACE_TYPE), POINTER :: SF
TYPE (VENTS_TYPE), POINTER :: VT
TYPE (OMESH_TYPE), POINTER :: OM
TYPE (MESH_TYPE), POINTER :: MM
REAL(EB), POINTER, DIMENSION(:,:) :: PBAR_P

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
 
HEAT_FLUX_LOOP: DO IW=1,NWC

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

   ! Compute specific heat for certain boundary conditions

   IF (SF%THERMAL_BC_INDEX /= NO_CONVECTION .AND. SF%THERMAL_BC_INDEX /= INTERPOLATED_BC) THEN
      ITMP = MIN(5000,NINT(0.5_EB*(TMP(IIG,JJG,KKG)+TMP_W(IW))))
      IF (N_SPECIES > 0) THEN
         YY_N = 0.5*(YYP(IIG,JJG,KKG,:)+YY_W(IW,:))
         CALL GET_AVERAGE_SPECIFIC_HEAT(YY_N,CP,ITMP)
      ELSE
         CP   = Y2CP_C(ITMP)
      ENDIF
   ENDIF
   METHOD_OF_HEAT_TRANSFER: SELECT CASE(SF%THERMAL_BC_INDEX)
 
      CASE (NO_CONVECTION) METHOD_OF_HEAT_TRANSFER

         TMP_F(IW) = TMP(IIG,JJG,KKG)
         TMP_W(IW) = TMP_F(IW)

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
            TMP_G     = TMP(IIG,JJG,KKG)
            RHOWAL    = 0.5_EB*(RHOP(IIG,JJG,KKG)+RHO_W(IW))
            CP_TERM   = MAX(0._EB,CP*UN*RHOWAL)
            TMP_W(IW) = ( (RDN(IW)*KW(IW)-0.5_EB*CP_TERM)*TMP_G + CP_TERM*TMP_F(IW)-QCONF(IW) )/(0.5_EB*CP_TERM+RDN(IW)*KW(IW))
            TMP_W(IW) = MAX(TMPMIN,TMP_W(IW))
         ELSE
            TMP_F(IW) = TMP(IIG,JJG,KKG)
            TMP_W(IW) = TMP_F(IW)
         ENDIF
         TMP(II,JJ,KK) = TMP_F(IW)
 
      CASE (SPECIFIED_TEMPERATURE) METHOD_OF_HEAT_TRANSFER

         TMP_G = TMP(IIG,JJG,KKG)
         IF (TW(IW)==T_BEGIN .AND. SF%RAMP_INDEX(TIME_TEMP)>=1) THEN
            TSI = T
         ELSE
            TSI = T - TW(IW)
         ENDIF
         TMP_F(IW) = TMP_0(KK) + EVALUATE_RAMP(TSI,SF%TAU(TIME_TEMP),SF%RAMP_INDEX(TIME_TEMP))*(SF%TMP_FRONT-TMP_0(KK))
         DTMP = TMP_G - TMP_F(IW)
         HEAT_TRANS_COEF(IW) = HEAT_TRANSFER_COEFFICIENT(IW,IIG,JJG,KKG,IOR,TMP_G,DTMP)
         QCONF(IW) = HEAT_TRANS_COEF(IW)*DTMP
         RHOWAL    = 0.5_EB*(RHOP(IIG,JJG,KKG)+RHO_W(IW))
         IF (N_SPECIES > 0) THEN
            CALL GET_AVERAGE_SPECIFIC_HEAT(YY_N,CP_F,MIN(5000,NINT(TMP_F(IW))))
         ELSE
            CP_F   = Y2CP_C(MIN(5000,NINT(TMP_F(IW))))
         ENDIF
         CP_TERM   = MAX(0._EB,-UW(IW))*RHOWAL*CP
         CP_TERM_F = MAX(0._EB,-UW(IW))*RHOWAL*CP_F
         TMP_W(IW) = ( (RDN(IW)*KW(IW)-0.5_EB*CP_TERM)*TMP_G + CP_TERM*TMP_F(IW) - QCONF(IW) )/(0.5_EB*CP_TERM+RDN(IW)*KW(IW))
         TMP_W(IW) = MAX(TMPMIN,TMP_W(IW))
         
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
            HEAT_TRANS_COEF(IW) = HEAT_TRANSFER_COEFFICIENT(IW,IIG,JJG,KKG,IOR,TMP_G,DTMP)
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
         HEAT_TRANS_COEF(IW) = HEAT_TRANSFER_COEFFICIENT(IW,IIG,JJG,KKG,IOR,TMP_G,DTMP)
         QCONF(IW) = HEAT_TRANS_COEF(IW)*DTMP
         RHOWAL    = 0.5_EB*(RHOP(IIG,JJG,KKG)+RHO_W(IW))
         CP_TERM   = MAX(0._EB,-CP*UW(IW)*RHOWAL)
         TMP_W(IW) = ( (RDN(IW)*KW(IW)-0.5_EB*CP_TERM)*TMP_G + CP_TERM*TMP_F(IW)-QCONF(IW) )/(0.5_EB*CP_TERM+RDN(IW)*KW(IW))
         TMP_W(IW) = MAX(TMPMIN,TMP_W(IW))

      CASE (CONVECTIVE_FLUX_BC) METHOD_OF_HEAT_TRANSFER
      
         IF (TW(IW)==T_BEGIN .AND. SF%RAMP_INDEX(TIME_HEAT)>=1) THEN
            TSI = T
         ELSE
            TSI = T - TW(IW)
         ENDIF
         RAMP_FACTOR = EVALUATE_RAMP(TSI,SF%TAU(TIME_HEAT),SF%RAMP_INDEX(TIME_HEAT))
         TMP_F(IW) = TMPA + RAMP_FACTOR*(SF%TMP_FRONT-TMPA)
         QCONF(IW) = -RAMP_FACTOR*SF%CONVECTIVE_HEAT_FLUX*AREA_ADJUST(IW)
         RHOWAL    = 0.5_EB*(RHOP(IIG,JJG,KKG)+RHO_W(IW))
         TMP_G     = TMP(IIG,JJG,KKG)
         CP_TERM   = MAX(0._EB,-CP*UW(IW)*RHOWAL)
         TMP_W(IW) = ( (RDN(IW)*KW(IW)-0.5_EB*CP_TERM)*TMP_G + CP_TERM*TMP_F(IW)-QCONF(IW) )/(0.5_EB*CP_TERM+RDN(IW)*KW(IW))
         TMP_W(IW) = MAX(TMPMIN,TMP_W(IW))
 
      CASE (INTERPOLATED_BC) METHOD_OF_HEAT_TRANSFER
 
         NOM   = IJKW(9,IW)
         OM    => OMESH(NOM)
         IF (PREDICTOR) THEN
            OM_RHOP => OM%RHOS
         ELSE
            OM_RHOP => OM%RHO
         ENDIF
         MM    => MESHES(NOM)
         RHO_G = RHOP(IIG,JJG,KKG)
         RHO_W(IW) = RHO_G
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
                  RHO_W(IW) = RHO_W(IW) + ARO*(OM_RHOP(IIO,JJO,KKO)-RHO_G)
               ENDDO
            ENDDO
         ENDDO
         IF (N_SPECIES==0) THEN
            TMP_W(IW) = PBAR_P(KK,PRESSURE_ZONE_WALL(IW))/(SPECIES(0)%RCON*RHO_W(IW))
         ELSE
            TMP_W(IW) = PBAR_P(KK,PRESSURE_ZONE_WALL(IW))/(RSUM_W(IW)*RHO_W(IW))
         ENDIF
         TMP_F(IW)     = TMP_W(IW)
         QCONF(IW)     = 0._EB
         TMP(II,JJ,KK) = TMP_W(IW)
 
   END SELECT METHOD_OF_HEAT_TRANSFER
 
   ! Record wall temperature in the ghost cell
 
   IF (SOLID(CELL_INDEX(II,JJ,KK))) TMP(II,JJ,KK) = MAX(100._EB,MIN(4900._EB,TMP_W(IW)))
 
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

USE MATH_FUNCTIONS, ONLY : EVALUATE_RAMP
REAL(EB) :: T,YY_WALL,DENOM,YY_OTHER(MAX_SPECIES),RHO_YY_OTHER(MAX_SPECIES),YY_EXTERIOR(MAX_SPECIES), &
            YY_G,YY_G_ALL(MAX_SPECIES),RHO_G,UN,DD,EPSB,MFT,TSI,ARO,UWO,RADIUS,AREA,RVC,RHO_NEW,MASS_ADDED
INTEGER :: I,IBC,IIG,JJG,KKG,IOR,IWB,IW,II,JJ,KK,N,IIO,JJO,KKO,N_INT_CELLS,NOM
TYPE (SURFACE_TYPE), POINTER :: SF
TYPE (VENTS_TYPE), POINTER :: VT
TYPE (OMESH_TYPE), POINTER :: OM
TYPE (MESH_TYPE), POINTER :: MM
TYPE (PARTICLE_CLASS_TYPE), POINTER :: PC
TYPE (DROPLET_TYPE), POINTER :: DR
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: YYP,OM_YYP
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW,RHOP,OM_RHOP
 
IF (PREDICTOR) THEN
   UU => U
   VV => V
   WW => W
   RHOP => RHOS
   IF (N_SPECIES > 0) YYP => YYS
ELSE
   UU => US
   VV => VS
   WW => WS
   RHOP => RHO
   IF (N_SPECIES > 0) YYP => YY
ENDIF 
 
! Loop through the wall cells, apply mass boundary conditions

WALL_CELL_LOOP: DO IW=1,NWC

   IF (BOUNDARY_TYPE(IW)==NULL_BOUNDARY) CYCLE WALL_CELL_LOOP

   IBC = IJKW(5,IW)

   ! Special cases that over-ride the boundary condition index, IBC

   IF (BOUNDARY_TYPE(IW)/=SOLID_BOUNDARY) THEN
      IF (BOUNDARY_TYPE(IW)==INTERPOLATED_BOUNDARY) IBC = INTERPOLATED_SURF_INDEX
      IF (BOUNDARY_TYPE(IW)==OPEN_BOUNDARY)         IBC = OPEN_SURF_INDEX
      IF (BOUNDARY_TYPE(IW)==MIRROR_BOUNDARY)       IBC = MIRROR_SURF_INDEX
   ENDIF
   IF (BOUNDARY_TYPE(IW)==POROUS_BOUNDARY) THEN
      IF (IJKW(9,IW)/=0) THEN
         IBC = INTERPOLATED_SURF_INDEX
      ELSE
         CYCLE WALL_CELL_LOOP
      ENDIF
   ENDIF

   ! Set the SURFace type

   SF  => SURFACE(IBC)

   ! Special cases

   IF (N_SPECIES==0 .AND. .NOT. SF%SPECIES_BC_INDEX==SPECIFIED_MASS_FLUX) CYCLE WALL_CELL_LOOP
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

         IF (.NOT.SOLID(CELL_INDEX(IIG,JJG,KKG))) YY_W(IW,1:N_SPECIES) = YYP(IIG,JJG,KKG,1:N_SPECIES)
         IF (BOUNDARY_TYPE(IW)==OPEN_BOUNDARY) THEN
            YY_EXTERIOR(1:N_SPECIES) = SPECIES(1:N_SPECIES)%YY0
            VT => VENTS(VENT_INDEX(IW))
            DO N=1,N_SPECIES
               IF (VT%MASS_FRACTION(N)>-1._EB) YY_EXTERIOR(N) = VT%MASS_FRACTION(N)
            ENDDO
            SELECT CASE(IOR)
               CASE( 1) 
                  IF (UU(II,JJ,KK)>0._EB)   YY_W(IW,1:N_SPECIES)=YY_EXTERIOR(1:N_SPECIES)
               CASE(-1)
                  IF (UU(II-1,JJ,KK)<0._EB) YY_W(IW,1:N_SPECIES)=YY_EXTERIOR(1:N_SPECIES)
               CASE( 2) 
                  IF (VV(II,JJ,KK)>0._EB)   YY_W(IW,1:N_SPECIES)=YY_EXTERIOR(1:N_SPECIES)
               CASE(-2)
                  IF (VV(II,JJ-1,KK)<0._EB) YY_W(IW,1:N_SPECIES)=YY_EXTERIOR(1:N_SPECIES)
               CASE( 3)
                  IF (WW(II,JJ,KK)>0._EB)   YY_W(IW,1:N_SPECIES)=YY_EXTERIOR(1:N_SPECIES)
               CASE(-3)
                  IF (WW(II,JJ,KK-1)<0._EB) YY_W(IW,1:N_SPECIES)=YY_EXTERIOR(1:N_SPECIES)
            END SELECT
            YYP(II,JJ,KK,1:N_SPECIES)=YY_W(IW,1:N_SPECIES)
         ENDIF
         IF ( SF%LEAK_PATH(1)>-1 .AND. UWS(IW)<0._EB) THEN
            RHO_G = RHOP(IIG,JJG,KKG)
            UN = -UWS(IW)
            IF (PREDICTOR) EPSB = -.5_EB*UN**2*DT*RDN(IW)
            IF (CORRECTOR) EPSB =  .5_EB*UN**2*DT*RDN(IW)
            SPECIES_LOOP_1: DO N=1,N_SPECIES
               DD    = RHODW(IW,N)*RDN(IW)
               YY_G  = YYP(IIG,JJG,KKG,N)
               DENOM = DD + (.5_EB*UN+EPSB)*RHO_W(IW)
               YY_W(IW,N) = ( MASSFLUX(IW,N) + YY_G*(DD + (EPSB-.5_EB*UN)*RHO_G) ) / DENOM
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
            IF_FLUX_LIMITER: IF (FLUX_LIMITER>=0) THEN
               YY_W(IW,N) = SPECIES(N)%YY0 + EVALUATE_RAMP(TSI,SF%TAU(N),SF%RAMP_INDEX(N))*(SF%MASS_FRACTION(N)-SPECIES(N)%YY0)
            ELSE IF_FLUX_LIMITER
               IF (UWS(IW)<0._EB .AND. SF%MASS_FRACTION(N) > 0._EB) THEN
                  YY_WALL = SPECIES(N)%YY0 + EVALUATE_RAMP(TSI,SF%TAU(N),SF%RAMP_INDEX(N))*(SF%MASS_FRACTION(N)-SPECIES(N)%YY0)
               ELSEIF (UWS(IW)>0._EB) THEN
                  YY_WALL = YYP(IIG,JJG,KKG,N)
               ELSE  ! This is usually where air is blown into a compartment and we don't want gain or loss of species N
                  RHO_G = RHOP(IIG,JJG,KKG)
                  UN = -UWS(IW)
                  IF (PREDICTOR) EPSB = -.5_EB*UN**2*DT*RDN(IW)
                  IF (CORRECTOR) EPSB =  .5_EB*UN**2*DT*RDN(IW)
                  DD    = RHODW(IW,N)*RDN(IW)
                  YY_G  = YYP(IIG,JJG,KKG,N)
                  DENOM = DD + (.5_EB*UN+EPSB)*RHO_W(IW)
                  YY_WALL = ( YY_G*(DD + (EPSB-.5_EB*UN)*RHO_G) ) / DENOM
               ENDIF
               IF (DNS) YY_W(IW,N) = 2._EB*YY_WALL - YYP(IIG,JJG,KKG,N)
               IF (LES) YY_W(IW,N) =       YY_WALL
            ENDIF IF_FLUX_LIMITER
         ENDDO
 
      CASE (SPECIFIED_MASS_FLUX) METHOD_OF_MASS_TRANSFER

         ! If the current time is before the "activation" time, TW, apply simple BCs and get out

         IF (T < TW(IW) .AND. N_SPECIES > 0) THEN
            IF (.NOT.SOLID(CELL_INDEX(IIG,JJG,KKG))) YY_W(IW,1:N_SPECIES)      = YYP(IIG,JJG,KKG,1:N_SPECIES)
            IF (PREDICTOR) UWS(IW)   = 0._EB 
            MASSFLUX(IW,1:N_SPECIES) = 0._EB
            MASSFLUX_ACTUAL(IW,1:N_SPECIES) = 0._EB
            CYCLE WALL_CELL_LOOP
         ENDIF

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
                  MASSFLUX_ACTUAL(IW,N) = MASSFLUX(IW,N)
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
         UN    = 2._EB*MFT/(RHO_W(IW)+RHO_G)
         IF (PREDICTOR) UWS(IW) = -UN
         IF (PREDICTOR) EPSB = -.5_EB*UN**2*DT*RDN(IW)
         IF (CORRECTOR) EPSB =  .5_EB*UN**2*DT*RDN(IW)
         SPECIES_LOOP: DO N=1,N_SPECIES
            DD    = RHODW(IW,N)*RDN(IW)
            YY_G  = YYP(IIG,JJG,KKG,N)
            DENOM = DD + (.5_EB*UN+EPSB)*RHO_W(IW)
            YY_W(IW,N) = ( MASSFLUX(IW,N) + YY_G*(DD + (EPSB-.5_EB*UN)*RHO_G) ) / DENOM
         ENDDO SPECIES_LOOP

      CASE (INTERPOLATED_BC) METHOD_OF_MASS_TRANSFER
 
         NOM   = IJKW(9,IW)
         OM    => OMESH(NOM)
         IF (PREDICTOR) THEN
            OM_RHOP => OM%RHOS
            OM_YYP  => OM%YYS
         ELSE
            OM_RHOP => OM%RHO
            OM_YYP  => OM%YY
         ENDIF
         MM    => MESHES(NOM)
         RHO_G = RHOP(IIG,JJG,KKG)
         YY_G_ALL(1:N_SPECIES) = YYP(IIG,JJG,KKG,1:N_SPECIES)
         RHO_YY_OTHER(1:N_SPECIES) = RHO_G*YY_G_ALL(1:N_SPECIES) ! Initialize summation of mass fluxes
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
                  RHO_YY_OTHER(1:N_SPECIES) = RHO_YY_OTHER(1:N_SPECIES) + &
                                 ARO*(OM_RHOP(IIO,JJO,KKO)*OM_YYP(IIO,JJO,KKO,1:N_SPECIES)-RHO_G*YY_G_ALL(1:N_SPECIES))
               ENDDO
            ENDDO
         ENDDO
         YY_OTHER(1:N_SPECIES) = RHO_YY_OTHER(1:N_SPECIES)/RHO_W(IW)
         YY_W(IW,1:N_SPECIES)  = YY_OTHER(1:N_SPECIES)
         YYP(II,JJ,KK,1:N_SPECIES) = YY_W(IW,1:N_SPECIES)

   END SELECT METHOD_OF_MASS_TRANSFER
   
!   ! pure periodic bcs for single mesh case
!   IF (FISHPAK_BC(1)==0) THEN
!      RHOP(0,:,:) = RHOP(IBAR,:,:)
!      RHOP(IBP1,:,:) = RHOP(1,:,:)
!   ENDIF
!   IF (FISHPAK_BC(2)==0) THEN
!      RHOP(:,0,:) = RHOP(:,JBAR,:)
!      RHOP(:,JBP1,:) = RHOP(:,1,:)
!   ENDIF
!   IF (FISHPAK_BC(3)==0) THEN
!      RHOP(:,:,0) = RHOP(:,:,KBAR)
!      RHOP(:,:,KBP1) = RHOP(:,:,1)
!   ENDIF
!   DO N=1,N_SPECIES
!      IF (FISHPAK_BC(1)==0) THEN
!         YYP(0,:,:,N) = YYP(IBAR,:,:,N)
!         YYP(IBP1,:,:,N) = YYP(1,:,:,N)
!      ENDIF
!      IF (FISHPAK_BC(2)==0) THEN
!         YYP(:,0,:,N) = YYP(:,JBAR,:,N)
!         YYP(:,JBP1,:,N) = YYP(:,1,:,N)
!      ENDIF
!      IF (FISHPAK_BC(3)==0) THEN
!         YYP(:,:,0,N) = YYP(:,:,KBAR,N)
!         YYP(:,:,KBP1,N) = YYP(:,:,1,N)
!      ENDIF
!   ENDDO

   ! Only set species mass fraction in the ghost cell if it is solid
    
   IF (IW<=NEWC .AND. N_SPECIES > 0 .AND. .NOT.SOLID(CELL_INDEX(II,JJ,KK)) .AND. .NOT.SOLID(CELL_INDEX(IIG,JJG,KKG))) &
      YYP(II,JJ,KK,1:N_SPECIES) = YY_W(IW,1:N_SPECIES)

ENDDO WALL_CELL_LOOP

! Add gases from virtual particles

IF (VIRTUAL_PARTICLES .AND. CORRECTOR .AND. N_SPECIES>0) THEN
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
      ELSE
         RADIUS = MAXVAL(SF%X_S)
      ENDIF
      SELECT CASE(SF%GEOMETRY)
         CASE(SURF_CARTESIAN)
            AREA = 2._EB*PC%LENGTH*PC%WIDTH
         CASE(SURF_CYLINDRICAL)
            AREA = TWOPI*RADIUS*PC%LENGTH
         CASE(SURF_SPHERICAL)
            AREA = 4._EB*PI*RADIUS**2
      END SELECT
      RVC = RDX(II)*RDY(JJ)*RDZ(KK)
      DO N=1,N_SPECIES
         MASS_ADDED = MASSFLUX(IW,N)*AREA*DT
         RHO_NEW = RHO(II,JJ,KK) + MASS_ADDED*RVC
         YYP(II,JJ,KK,N) = (RHO(II,JJ,KK)*YYP(II,JJ,KK,N) + MASS_ADDED*RVC)/RHO_NEW
         RHO(II,JJ,KK) = RHO_NEW
      ENDDO
   ENDDO DROPLET_LOOP
ENDIF

END SUBROUTINE SPECIES_BC
 

 
SUBROUTINE DENSITY_BC
 
! Compute density at wall from wall temperatures and mass fractions 

USE PHYSICAL_FUNCTIONS, ONLY : GET_SPECIFIC_GAS_CONSTANT
REAL(EB) :: WFAC,YY_GET(1:N_SPECIES)
INTEGER  :: IBC,IIG,JJG,KKG,IW,II,JJ,KK, N
REAL(EB), POINTER, DIMENSION(:,:) :: PBAR_P
TYPE (SURFACE_TYPE), POINTER :: SF
REAL(EB), POINTER, DIMENSION(:,:,:) :: RHOP
 
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
      YY_GET = YY_W(IW,:)
      CALL GET_SPECIFIC_GAS_CONSTANT(YY_GET,RSUM_W(IW))
   ENDIF
 
! Compute ghost cell density
   IF (BOUNDARY_TYPE(IW)/=INTERPOLATED_BOUNDARY) RHO_W(IW) = PBAR_P(KK,PRESSURE_ZONE_WALL(IW))/(RSUM_W(IW)*TMP_W(IW))

! Actually set the ghost cell value of density in the ghost cell if it is a solid wall
   IF ( (SOLID(CELL_INDEX(II,JJ,KK)) .AND. .NOT.SOLID(CELL_INDEX(IIG,JJG,KKG))) .OR.  &
         BOUNDARY_TYPE(IW)==OPEN_BOUNDARY .OR. BOUNDARY_TYPE(IW)==INTERPOLATED_BOUNDARY) THEN
      RHOP(II,JJ,KK) = RHO_W(IW)
      IF (N_SPECIES>0) RSUM(II,JJ,KK) = RSUM_W(IW)
   ENDIF
 
ENDDO WALL_CELL_LOOP

END SUBROUTINE DENSITY_BC
 
 
SUBROUTINE PYROLYSIS(T,DT_BC)

! Loop through all the boundary cells that require a 1-D heat transfer calc

USE PHYSICAL_FUNCTIONS, ONLY: GET_MOLECULAR_WEIGHT, GET_AVERAGE_SPECIFIC_HEAT
USE GEOMETRY_FUNCTIONS
USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
REAL(EB) :: DTMP,QNETF,QDXKF,QDXKB,RR,TMP_G,T,RFACF,RFACB,RFACF2,RFACB2, &
            PPCLAUS,PPSURF,TMP_G_B,RHOWAL,CP_TERM,DT_BC, CP,&
            DXKF,DXKB,REACTION_RATE,QCONB,DELTA_RHO_S,QRADINB,RFLUX_UP,RFLUX_DOWN,E_WALLB, &
            HVRG,Y_MF_G,Y_MF_W,YY_S, RSUM_W, RSUM_G, RSUM_S, MFLUX, MFLUX_S, VOLSUM,YPRSUM, &
            DXF, DXB,HTCB,Q_WATER_F,Q_WATER_B,TMP_F_OLD, DX_GRID, RHO_S0,DT2_BC,TOLERANCE,C_S_ADJUST_UNITS
INTEGER :: IBC,IIG,JJG,KKG,IIB,JJB,KKB,IWB,NWP,I,J,NR,NNN,NL,II,JJ,KK,IW,IOR,N,I_OBST,NS,ITMP
REAL(EB) :: SMALLEST_CELL_SIZE(MAX_LAYERS),THICKNESS,LAYER_THICKNESS_NEW(MAX_LAYERS),YY_G_ALL(1:N_SPECIES),YY_GET(1:N_SPECIES)
REAL(EB),ALLOCATABLE,DIMENSION(:) :: TMP_W_NEW
REAL(EB),ALLOCATABLE,DIMENSION(:,:) :: INT_WGT
INTEGER  :: N_LAYER_CELLS_NEW(MAX_LAYERS), NWP_NEW,I_GRAD,STEPCOUNT
LOGICAL :: POINT_SHRINK, RECOMPUTE,ITERATE
TYPE (WALL_TYPE), POINTER :: WC
TYPE (SURFACE_TYPE), POINTER :: SF
TYPE (MATERIAL_TYPE), POINTER :: ML
 
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
   IF (WC%BURNAWAY) CYCLE WALL_CELL_LOOP

   ITMP = MIN(5000,NINT(0.5_EB*(TMP(IIG,JJG,KKG)+TMP_W(IW))))
   IF (N_SPECIES > 0) THEN
      YY_G_ALL(:) = 0.5_EB*(YY(IIG,JJG,KKG,:)+YY_W(IW,:))
      CALL GET_AVERAGE_SPECIFIC_HEAT(YY_G_ALL,CP,ITMP)
   ELSE
      CP = Y2CP_C(ITMP)
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
   HEAT_TRANS_COEF(IW) = HEAT_TRANSFER_COEFFICIENT(IW,IIG,JJG,KKG,IOR,TMP_G,DTMP)
   QCONF(IW) = HEAT_TRANS_COEF(IW)*DTMP
 
   ! Get heat losses from convection and radiation out of back of surface
  
   E_WALLB = MATERIAL(SF%LAYER_MATL_INDEX(SF%N_LAYERS,1))%EMISSIVITY
         
   SELECT CASE(SF%BACKING)
      CASE(VOID)  ! Non-insulated backing to an ambient void
         DTMP = SF%TMP_BACK - TMP_B(IW)
         HTCB = HEAT_TRANSFER_COEFFICIENT(IW,-1,-1,-1,IOR,SF%TMP_BACK,DTMP)
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
            HTCB = HEAT_TRANSFER_COEFFICIENT(IWB,IIB,JJB,KKB,IOR,TMP_G_B,DTMP)     
            HEAT_TRANS_COEF(IWB) = HTCB
            QRADINB  = QRADIN(IWB)
            IF (NLP>0) Q_WATER_B = -SUM(WCPUA(IWB,:))
         ELSE
            TMP_G_B  = TMPA
            DTMP = TMP_G_B - TMP_B(IW)
            HTCB = HEAT_TRANSFER_COEFFICIENT(IW,-1,-1,-1,IOR,TMP_G_B,DTMP)
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
            MFLUX_S = MF_FRAC(I)*REACTION_RATE/RDX_S(I)/SF%THICKNESS**I_GRAD
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
      IF (MASSFLUX(IW,I_FUEL) > 0._EB) TW(IW) = T
      IF (I_WATER > 0) THEN
         IF (MASSFLUX(IW,I_WATER) > 0._EB) TW(IW) = T
      ENDIF
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
      YY_GET = YY(IIG,JJG,KKG,:)
      Y_MF_G = YY_GET(I_FUEL)
      CALL GET_MOLECULAR_WEIGHT(YY_GET,RSUM_G)
      ! wall values
      YY_GET = YY_W(IW,:)
      Y_MF_W = YY_GET(I_FUEL)
      CALL GET_MOLECULAR_WEIGHT(YY_GET,RSUM_W)
      ! Weighted average of wall and gas values
      ! Alvernative 1
      YPRSUM  = 0.2*(Y_MF_W/RSUM_W) + 0.8*(Y_MF_G/RSUM_G)
!      YPRSUM  = 0.0*(Y_MF_W/RSUM_W) + 1.0*(Y_MF_G/RSUM_G)
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
         RHOWAL               = 0.5_EB*(RHO(IIG,JJG,KKG)+RHO_W(IW))
         CP_TERM              = MAX(0._EB,-CP_GAMMA*UW(IW)*RHOWAL)
         QCONF(IW)            = 0._EB
         TMP_W(IW)            = ( (RDN(IW)*KW(IW)-0.5_EB*CP_TERM)*TMP_G+CP_TERM*TMP_F(IW)-QCONF(IW) ) & 
                                /(0.5_EB*CP_TERM+RDN(IW)*KW(IW))
         TMP_W(IW)            = MAX(TMPMIN,TMP_W(IW))
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
         IF (I.EQ.1) E_WALL(IW) = E_WALL(IW) + WC%RHO_S(I,N)*ML%EMISSIVITY/ML%RHO_S
         RHO_S(I)   = RHO_S(I) + WC%RHO_S(I,N)
         
      ENDDO MATERIAL_LOOP3
      
      IF (VOLSUM > 0._EB) THEN
         K_S(I) = K_S(I)/VOLSUM
         IF (I.EQ.1) E_WALL(IW) = E_WALL(IW)/VOLSUM
      ENDIF
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
 
   ! Determine ghost wall temperature
 
   RHOWAL    = 0.5_EB*(RHO(IIG,JJG,KKG)+RHO_W(IW))
   CP_TERM   = MAX(0._EB,-CP_GAMMA*UW(IW)*RHOWAL)
   QCONF(IW) = HEAT_TRANS_COEF(IW) * (TMP_G - 0.5_EB * (TMP_F(IW) + TMP_F_OLD) )
   TMP_W(IW) = ( (RDN(IW)*KW(IW)-0.5_EB*CP_TERM)*TMP_G + CP_TERM*TMP_F(IW)-QCONF(IW) )/(0.5_EB*CP_TERM+RDN(IW)*KW(IW))
   TMP_W(IW) = MAX(TMPMIN,TMP_W(IW))
   IF (SOLID(CELL_INDEX(II,JJ,KK))) TMP(II,JJ,KK) = MAX(100._EB,MIN(4900._EB,TMP_W(IW)))

ENDDO WALL_CELL_LOOP
 
END SUBROUTINE PYROLYSIS
 
 
REAL(EB) FUNCTION HEAT_TRANSFER_COEFFICIENT(IW,IIG,JJG,KKG,IOR,TMP_G,DELTA_TMP)

INTEGER  :: IW,IIG,JJG,KKG,IOR,ITMP
REAL(EB) :: TMP_G,DELTA_TMP,U2,V2,W2,VELCON,H_NATURAL,H_FORCED,H_DNS
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW,RHOP
 
IF (H_FIXED >= 0._EB) THEN
   HEAT_TRANSFER_COEFFICIENT = H_FIXED
   RETURN
ENDIF

IF (PREDICTOR) THEN
   UU => U
   VV => V
   WW => W
   RHOP => RHOS
ELSE
   UU => US
   VV => VS
   WW => WS
   RHOP => RHO
ENDIF 

IF (DNS) THEN
   HEAT_TRANSFER_COEFFICIENT = 2._EB*KW(IW)*RDN(IW)
ELSE
   IF (IIG<0) THEN   ! No gas cell information available
      SELECT CASE(ABS(IOR))
         CASE(0:2)
            H_NATURAL = HCV*ABS(DELTA_TMP)**ONTH
         CASE(3)
            H_NATURAL = HCH*ABS(DELTA_TMP)**ONTH
      END SELECT
      H_FORCED = 0._EB
   ELSE
      SELECT CASE(ABS(IOR))
         CASE(0)
            U2 = 0.25_EB*(UU(IIG,JJG,KKG)+UU(IIG-1,JJG,KKG))**2
            V2 = 0.25_EB*(VV(IIG,JJG,KKG)+VV(IIG,JJG-1,KKG))**2
            W2 = 0.25_EB*(WW(IIG,JJG,KKG)+WW(IIG,JJG,KKG-1))**2
            VELCON = (U2+V2+W2)**0.4_EB
            H_NATURAL = HCV*ABS(DELTA_TMP)**ONTH
         CASE(1)
            V2 = 0.25_EB*(VV(IIG,JJG,KKG)+VV(IIG,JJG-1,KKG))**2
            W2 = 0.25_EB*(WW(IIG,JJG,KKG)+WW(IIG,JJG,KKG-1))**2
            VELCON = (V2+W2)**0.4_EB
            H_NATURAL = HCV*ABS(DELTA_TMP)**ONTH
         CASE(2)
            U2 = 0.25_EB*(UU(IIG,JJG,KKG)+UU(IIG-1,JJG,KKG))**2
            W2 = 0.25_EB*(WW(IIG,JJG,KKG)+WW(IIG,JJG,KKG-1))**2
            VELCON = (U2+W2)**0.4_EB
            H_NATURAL = HCV*ABS(DELTA_TMP)**ONTH
         CASE(3)
            U2 = 0.25_EB*(UU(IIG,JJG,KKG)+UU(IIG-1,JJG,KKG))**2
            V2 = 0.25_EB*(VV(IIG,JJG,KKG)+VV(IIG,JJG-1,KKG))**2
            VELCON = (U2+V2)**0.4_EB
            H_NATURAL = HCH*ABS(DELTA_TMP)**ONTH
      END SELECT
      H_FORCED  = C_FORCED*VELCON*RHOP(IIG,JJG,KKG)**0.8_EB
   ENDIF
   ITMP = MIN(5000,NINT(TMP_G))
   H_DNS = Y2MU_C(ITMP)*CP_GAMMA*RPR*2._EB*RDN(IW)
   HEAT_TRANSFER_COEFFICIENT = MAX(H_DNS,H_FORCED,H_NATURAL)
ENDIF

END FUNCTION HEAT_TRANSFER_COEFFICIENT

END SUBROUTINE WALL_BC

SUBROUTINE GET_REV_wall(MODULE_REV,MODULE_DATE)
INTEGER,INTENT(INOUT) :: MODULE_REV
CHARACTER(255),INTENT(INOUT) :: MODULE_DATE

WRITE(MODULE_DATE,'(A)') wallrev(INDEX(wallrev,':')+1:LEN_TRIM(wallrev)-2)
READ (MODULE_DATE,'(I5)') MODULE_REV
WRITE(MODULE_DATE,'(A)') walldate

END SUBROUTINE GET_REV_wall
 
END MODULE WALL_ROUTINES

MODULE SOOT
 
! Compute soot deposition 
 
USE PRECISION_PARAMETERS
USE MESH_POINTERS

IMPLICIT NONE
PRIVATE
CHARACTER(255), PARAMETER :: sootid='$Id$'
CHARACTER(255), PARAMETER :: sootrev='$Revision$'
CHARACTER(255), PARAMETER :: sootdate='$Date$'

PUBLIC CALC_SOOT_DEPOSITION
 
 
CONTAINS
 
SUBROUTINE CALC_SOOT_DEPOSITION(NM)
USE PHYSICAL_FUNCTIONS, ONLY: GET_VISCOSITY
USE GLOBAL_CONSTANTS, ONLY: EVACUATION_ONLY,SOLID_PHASE_ONLY,SOLID_BOUNDARY,I_PROG_SOOT,N_SPECIES,TUSED
USE COMP_FUNCTIONS, ONLY:SECOND
INTEGER, INTENT(IN) :: NM
REAL(EB), PARAMETER :: CS=1.147_EB,CT=2.18_EB,CM=1.146_EB,KPKG=1._EB,PARTD=2.E-8_EB
REAL(EB) :: U1,U_THERM,U_TURB,TGAS,TWALL,MUGAS,Y_SOOT,RHOG,TNOW,YIN(1:N_SPECIES),YDEP,NU_OVER_DZ
INTEGER  :: II,JJ,KK,IIG,JJG,KKG,IW,IOR,ITMP

REAL(EB), PARAMETER :: A=8.3_EB,B=1._EB/7._EB,Z_PLUS_TURBULENT = 11.81_EB,ALPHA=7.202125273562269_EB 
REAL(EB), PARAMETER :: BETA=1._EB+B,ETA=(1._EB+B)/A,GAMMA=2._EB/(1._EB+B)

 
IF (EVACUATION_ONLY(NM)) RETURN
IF (SOLID_PHASE_ONLY) RETURN

TNOW=SECOND()
CALL POINT_TO_MESH(NM)
 
WALL_CELL_LOOP: DO IW=1,NWC+NVWC
   IF (BOUNDARY_TYPE(IW)/=SOLID_BOUNDARY .OR. UW(IW) < 0._EB) CYCLE WALL_CELL_LOOP
   II  = IJKW(1,IW)
   JJ  = IJKW(2,IW)
   KK  = IJKW(3,IW)
   IOR = IJKW(4,IW)
   IIG = IJKW(6,IW)
   JJG = IJKW(7,IW)
   KKG = IJKW(8,IW)
   YIN = YY(IIG,JJG,KKG,:)
   IF (YIN(I_PROG_SOOT) < 1.E-14_EB) CYCLE WALL_CELL_LOOP
   TGAS = TMP(IIG,JJG,KKG)
   TWALL = TMP_F(IW)
   SELECT CASE(ABS(IOR))
      CASE(1)
         U1 = SQRT(0.25_EB*(V(IIG,JJG-1,KKG)+V(IIG,JJG,KKG))**2+0.25_EB*(W(IIG,JJG,KKG-1)+W(IIG,JJG,KKG))**2)
      CASE(2)
         U1 = SQRT(0.25_EB*(U(IIG-1,JJG,KKG)+U(IIG,JJG,KKG))**2+0.25_EB*(W(IIG,JJG,KKG-1)+W(IIG,JJG,KKG))**2)
      CASE(3)
         U1 = SQRT(0.25_EB*(U(IIG-1,JJG,KKG)+U(IIG,JJG,KKG))**2+0.25_EB*(V(IIG,JJG-1,KKG)+V(IIG,JJG,KKG))**2)
   END SELECT   
   ITMP = MIN(NINT(TGAS),5000)
   CALL GET_VISCOSITY(YIN,MUGAS,ITMP)
   RHOG=RHO(IIG,JJG,KKG)
   NU_OVER_DZ=MUGAS/RHOG*RDN(IW)
   !MFP = MUGAS * SQRT(PI/(2._EB*PBAR(KKG,PRESSURE_ZONE(IIG,JJG,KKG))*RHOG))
   !KN = 2._EB*MFP/PARTD
   !CC = 1+ KN * (1.257_EB+0.4_EB*EXP(-1.1/KN))
   !U_THERM = 2._EB *CC*CS*(KPKG+CT*KN)*(TGAS-TWALL)*RDN(IW)/TGAS*MUGAS/RHOG/((1._EB+3._EB*CM*KN)*(1._EB+2._EB*KPKG+2._EB*CT*KN))
   U_THERM = 0.55_EB*(TGAS-TWALL)*RDN(IW)/TGAS*MUGAS/RHOG
   U_TURB = 0.175_EB*SQRT((ALPHA*NU_OVER_DZ**BETA + ETA*NU_OVER_DZ**B*ABS(U1))**GAMMA) 
   IF (U_THERM+U_TURB < 0._EB) CYCLE WALL_CELL_LOOP
   
   YIN = YIN * RHOG  
   Y_SOOT = YIN(I_PROG_SOOT)   
   YDEP =Y_SOOT*MIN(1._EB,(U_THERM+U_TURB)*DT*RDN(IW))
   YIN(I_PROG_SOOT) = Y_SOOT - YDEP      
   AWMSOOT(IW)=AWMSOOT(IW)+YDEP/RDN(IW)
   RHO(IIG,JJG,KKG) = RHOG - YDEP
   YY(IIG,JJG,KKG,:) = YIN / RHO(IIG,JJG,KKG)
   
ENDDO WALL_CELL_LOOP

TUSED(10,NM)=TUSED(10,NM)+SECOND()-TNOW

END SUBROUTINE CALC_SOOT_DEPOSITION


!---------------------------------------------------------------------------

SUBROUTINE GET_REV_soot(MODULE_REV,MODULE_DATE)
INTEGER,INTENT(INOUT) :: MODULE_REV
CHARACTER(255),INTENT(INOUT) :: MODULE_DATE
WRITE(MODULE_DATE,'(A)') sootrev(INDEX(sootrev,':')+1:LEN_TRIM(sootrev)-2)
READ (MODULE_DATE,'(I5)') MODULE_REV
WRITE(MODULE_DATE,'(A)') sootdate
END SUBROUTINE GET_REV_soot
 
END MODULE SOOT
