MODULE WALL_ROUTINES
 
! Compute the wall boundary conditions
 
USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_POINTERS
 
IMPLICIT NONE
PRIVATE
CHARACTER(255), PARAMETER :: wallid='$Id$'
PUBLIC WALL_BC
 
CONTAINS

SUBROUTINE WALL_BC(T,NM)
USE COMP_FUNCTIONS, ONLY: SECOND
REAL(EB) :: TNOW
REAL(EB), INTENT(IN) :: T
INTEGER, INTENT(IN) :: NM

IF (EVACUATION_ONLY(NM)) RETURN

TNOW=SECOND()

CALL POINT_TO_MESH(NM)

IF (.NOT. ISOTHERMAL) CALL THERMAL_BC(T)
IF (N_SPECIES > 0) CALL SPECIES_BC(T)
CALL DENSITY_BC

TUSED(6,NM)=TUSED(6,NM)+SECOND()-TNOW

CONTAINS
 
SUBROUTINE THERMAL_BC(T)
USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP 
! Thermal boundary conditions for adiabatic, fixed temperature, fixed flux and interpolated boundaries.
! One dimensional heat transfer and pyrolysis is done in PYROLYSIS, which is called at the end of this routine.
 
REAL(EB) :: DT_BC,T,TMP_G,DTMP,TMP_OTHER,CP_TERM,RHOWAL,RAMP_FACTOR,QNET,FDERIV
INTEGER  :: IOR,II,JJ,KK,IBC,IIG,JJG,KKG, IW
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW,RHOP
TYPE (SURFACE_TYPE), POINTER :: SF
 
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
 
! Loop through all boundary cells and apply heat transfer method, except for thermally-thick cells
 
HEAT_FLUX_LOOP: DO IW=1,NWC
   IF (BOUNDARY_TYPE(IW)==NULL_BOUNDARY .OR. BOUNDARY_TYPE(IW)==POROUS_BOUNDARY) CYCLE HEAT_FLUX_LOOP
   II  = IJKW(1,IW)
   JJ  = IJKW(2,IW)
   KK  = IJKW(3,IW)
   IIG = IJKW(6,IW)
   JJG = IJKW(7,IW)
   KKG = IJKW(8,IW)
   IOR = IJKW(4,IW)
   IBC = IJKW(5,IW)
   SF  => SURFACE(IBC)
 
   METHOD_OF_HEAT_TRANSFER: SELECT CASE(SF%THERMAL_BC_INDEX)
 
      CASE (ZERO_GRADIENT) METHOD_OF_HEAT_TRANSFER 
 
         TMP_F(IW) = TMP(IIG,JJG,KKG)
         SELECT CASE(IOR)
            CASE( 1) 
               IF (UU(II,JJ,KK)>=0._EB)   TMP_F(IW) = TMP_0(KK)
            CASE(-1) 
               IF (UU(II-1,JJ,KK)<=0._EB) TMP_F(IW) = TMP_0(KK)
            CASE( 2) 
               IF (VV(II,JJ,KK)>=0._EB)   TMP_F(IW) = TMP_0(KK)
            CASE(-2) 
               IF (VV(II,JJ-1,KK)<=0._EB) TMP_F(IW) = TMP_0(KK)
            CASE( 3) 
               IF (WW(II,JJ,KK)>=0._EB)   TMP_F(IW) = TMP_0(KK)
            CASE(-3) 
               IF (WW(II,JJ,KK-1)<=0._EB) TMP_F(IW) = TMP_0(KK)
            END SELECT
         TMP(II,JJ,KK) = TMP_F(IW)
         TMP_W(IW) = TMP_F(IW)
 
      CASE (SPECIFIED_TEMPERATURE) METHOD_OF_HEAT_TRANSFER
         TMP_G = TMP(IIG,JJG,KKG)
         TMP_F(IW) = TMP_0(KK) + EVALUATE_RAMP(T-TW(IW),SF%TAU(TIME_TEMP),SF%RAMP_INDEX(TIME_TEMP))*(SF%TMP_FRONT-TMP_0(KK))
         DTMP = TMP_G - TMP_F(IW)
         HEAT_TRANS_COEF(IW) = HEAT_TRANSFER_COEFFICIENT(IW,IIG,JJG,KKG,IOR,TMP_G,DTMP)
         QCONF(IW) = HEAT_TRANS_COEF(IW)*DTMP
         RHOWAL    = 0.5_EB*(RHOP(IIG,JJG,KKG)+RHO_W(IW))
         CP_TERM   = MAX(0._EB,-CP_GAMMA*UW(IW)*RHOWAL)
         TMP_W(IW) = ( (RDN(IW)*KW(IW)-0.5_EB*CP_TERM)*TMP_G + CP_TERM*TMP_F(IW)-QCONF(IW) )/(0.5_EB*CP_TERM+RDN(IW)*KW(IW))
         TMP_W(IW) = MAX(TMPMIN,TMP_W(IW))
      CASE (ADIABATIC_INDEX) METHOD_OF_HEAT_TRANSFER
         TMP_G = TMP(IIG,JJG,KKG)
         TMP_OTHER = TMP_F(IW)
         ADLOOP: DO
            DTMP = TMP_G - TMP_OTHER
            HEAT_TRANS_COEF(IW) = HEAT_TRANSFER_COEFFICIENT(IW,IIG,JJG,KKG,IOR,TMP_G,DTMP)
            QNET = HEAT_TRANS_COEF(IW)*DTMP + QRADIN(IW) - E_WALL(IW) * SIGMA * TMP_OTHER ** 4
            FDERIV = -HEAT_TRANS_COEF(IW) -  4._EB * E_WALL(IW) * SIGMA * TMP_OTHER ** 3
            IF (FDERIV /= 0._EB) TMP_OTHER = TMP_OTHER - QNET / FDERIV
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
         CP_TERM   = MAX(0._EB,-CP_GAMMA*UW(IW)*RHOWAL)
         TMP_W(IW) = ( (RDN(IW)*KW(IW)-0.5_EB*CP_TERM)*TMP_G + CP_TERM*TMP_F(IW)-QCONF(IW) )/(0.5_EB*CP_TERM+RDN(IW)*KW(IW))
         TMP_W(IW) = MAX(TMPMIN,TMP_W(IW))

      CASE (SPECIFIED_HEAT_FLUX) METHOD_OF_HEAT_TRANSFER

         RAMP_FACTOR = EVALUATE_RAMP(T-TW(IW),SF%TAU(TIME_HEAT),SF%RAMP_INDEX(TIME_HEAT))
         TMP_F(IW) = TMPA + RAMP_FACTOR*(SF%TMP_FRONT-TMPA)
         QCONF(IW) = -RAMP_FACTOR*SF%CONVECTIVE_HEAT_FLUX*AREA_ADJUST(IW)
         RHOWAL    = 0.5_EB*(RHOP(IIG,JJG,KKG)+RHO_W(IW))
         TMP_G     = TMP(IIG,JJG,KKG)
         CP_TERM   = MAX(0._EB,-CP_GAMMA*UW(IW)*RHOWAL)
         TMP_W(IW) = ( (RDN(IW)*KW(IW)-0.5_EB*CP_TERM)*TMP_G + CP_TERM*TMP_F(IW)-QCONF(IW) )/(0.5_EB*CP_TERM+RDN(IW)*KW(IW))
         TMP_W(IW) = MAX(TMPMIN,TMP_W(IW))
 
      CASE (INTERPOLATED_BC) METHOD_OF_HEAT_TRANSFER
 
         TMP_G = TMP(IIG,JJG,KKG)
         TMP_OTHER =OMESH(IJKW(9,IW))%TMP(IJKW(10,IW),IJKW(11,IW),IJKW(12,IW))
         IF (CELL_VOLUME_RATIO(IW)<0.5 .OR.  CELL_VOLUME_RATIO(IW)>2.0) THEN
            TMP_W(IW) = TMP_G
            SELECT CASE(IOR)
               CASE( 1) 
                  IF (UU(II,JJ,KK)>=0._EB)   TMP_W(IW) = TMP_OTHER
               CASE(-1) 
                  IF (UU(II-1,JJ,KK)<=0._EB) TMP_W(IW) = TMP_OTHER
               CASE( 2) 
                  IF (VV(II,JJ,KK)>=0._EB)   TMP_W(IW) = TMP_OTHER
               CASE(-2) 
                  IF (VV(II,JJ-1,KK)<=0._EB) TMP_W(IW) = TMP_OTHER
               CASE( 3) 
                  IF (WW(II,JJ,KK)>=0._EB)   TMP_W(IW) = TMP_OTHER
               CASE(-3) 
                  IF (WW(II,JJ,KK-1)<=0._EB) TMP_W(IW) = TMP_OTHER
            END SELECT
         ELSE
            TMP_W(IW) = TMP_OTHER
         ENDIF
         TMP_F(IW) = TMP_W(IW)
         QCONF(IW) = KW(IW)*(TMP_G-TMP_W(IW))*RDN(IW)
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
 
REAL(EB) T,YY_WALL,YY_G,DENOM,YY_OTHER(1:20),RHO_G,UN,DD,EPSB,MFT,TSI
INTEGER IBC,IIG,JJG,KKG,IOR,IC,IWB,IW, II, JJ, KK, N
TYPE (SURFACE_TYPE), POINTER :: SF
TYPE (MATERIAL_TYPE), POINTER :: ML
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: YYP
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW,RHOP
 
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

WALL_CELL_LOOP: DO IW=1,NWC
   IF (BOUNDARY_TYPE(IW)==NULL_BOUNDARY .OR. BOUNDARY_TYPE(IW)==POROUS_BOUNDARY) CYCLE WALL_CELL_LOOP
   IBC = IJKW(5,IW)
   SF  => SURFACE(IBC)
   IF (N_SPECIES==0 .AND. SF%MASS_FLUX(0)==0._EB) CYCLE WALL_CELL_LOOP
   II  = IJKW(1,IW)
   JJ  = IJKW(2,IW)
   KK  = IJKW(3,IW)
   IOR = IJKW(4,IW)
   IIG = IJKW(6,IW)
   JJG = IJKW(7,IW)
   KKG = IJKW(8,IW)
   
   METHOD_OF_MASS_TRANSFER: SELECT CASE(SF%SPECIES_BC_INDEX)
 
      CASE (NO_MASS_FLUX) METHOD_OF_MASS_TRANSFER 

         IF (.NOT.SOLID(CELL_INDEX(IIG,JJG,KKG))) YY_W(IW,1:N_SPECIES) = YYP(IIG,JJG,KKG,1:N_SPECIES)
         IF (BOUNDARY_TYPE(IW)==OPEN_BOUNDARY) THEN
            SELECT CASE(IOR)
               CASE( 1) 
                  IF (UU(II,JJ,KK)>0._EB)   YY_W(IW,1:N_SPECIES)=SPECIES(1:N_SPECIES)%YY0
               CASE(-1)
                  IF (UU(II-1,JJ,KK)<0._EB) YY_W(IW,1:N_SPECIES)=SPECIES(1:N_SPECIES)%YY0
               CASE( 2) 
                  IF (VV(II,JJ,KK)>0._EB)   YY_W(IW,1:N_SPECIES)=SPECIES(1:N_SPECIES)%YY0
               CASE(-2)
                  IF (VV(II,JJ-1,KK)<0._EB) YY_W(IW,1:N_SPECIES)=SPECIES(1:N_SPECIES)%YY0
               CASE( 3)
                  IF (WW(II,JJ,KK)>0._EB)   YY_W(IW,1:N_SPECIES)=SPECIES(1:N_SPECIES)%YY0
               CASE(-3)
                  IF (WW(II,JJ,KK-1)<0._EB) YY_W(IW,1:N_SPECIES)=SPECIES(1:N_SPECIES)%YY0
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
         IF (PREDICTOR) TSI = T + DT - TW(IW)
         IF (CORRECTOR) TSI = T      - TW(IW)
         DO N=1,N_SPECIES
            IF (UWS(IW) <= 0._EB) THEN
               YY_WALL = SPECIES(N)%YY0 + EVALUATE_RAMP(TSI,SF%TAU(N),SF%RAMP_INDEX(N))*(SF%MASS_FRACTION(N)-SPECIES(N)%YY0)
            ELSE
               YY_WALL = YYP(IIG,JJG,KKG,N)
            ENDIF
            IF (DNS) YY_W(IW,N) = 2._EB*YY_WALL - YYP(IIG,JJG,KKG,N)
            IF (LES) YY_W(IW,N) =       YY_WALL 
         ENDDO
 
      CASE (SPECIFIED_MASS_FLUX) METHOD_OF_MASS_TRANSFER

         ! If the current time is before the "activation" time, TW, apply simple BCs and get out

         IF (T < TW(IW)) THEN
            IF (.NOT.SOLID(CELL_INDEX(IIG,JJG,KKG))) YY_W(IW,1:N_SPECIES)      = YYP(IIG,JJG,KKG,1:N_SPECIES)
            IF (SOLID(CELL_INDEX(II,JJ,KK)))         YYP(II,JJ,KK,1:N_SPECIES) = YY_W(IW,1:N_SPECIES)
            IF (PREDICTOR) UWS(IW)   = 0._EB 
            MASSFLUX(IW,1:N_SPECIES) = 0._EB
            ACTUAL_BURN_RATE(IW)     = 0._EB
            CYCLE WALL_CELL_LOOP
         ENDIF

         MFT = 0._EB

         ! If the user has specified the burning rate, evaluate the ramp and other related parameters

         SUM_MASSFLUX_LOOP: DO N=0,N_SPECIES
            IF (SF%MASS_FLUX(N) > 0._EB) THEN  ! Use user-specified ramp-up of mass flux
               TSI = T - TW(IW)
               MASSFLUX(IW,N) = EVALUATE_RAMP(TSI,SF%TAU(N),SF%RAMP_INDEX(N))*SF%MASS_FLUX(N)
               IF (N==I_FUEL) THEN
                  IF (EW(IW)>0._EB) MASSFLUX(IW,N) = MASSFLUX(IW,N)*EXP(-EW(IW))
                  ACTUAL_BURN_RATE(IW) = MASSFLUX(IW,N)
               ENDIF
            ENDIF
            MASSFLUX(IW,N) = MASSFLUX(IW,N)*AREA_ADJUST(IW)
            IF (N==I_FUEL) ACTUAL_BURN_RATE(IW) = ACTUAL_BURN_RATE(IW)*AREA_ADJUST(IW)
            MFT = MFT + MASSFLUX(IW,N)
         ENDDO SUM_MASSFLUX_LOOP
 
         ! Add total fuel consumed to various summing arrays

         CONSUME_FUEL: IF (CORRECTOR .AND. SF%THERMALLY_THICK) THEN  
            IF (SF%SURFACE_DENSITY>0._EB .AND. SF%BACKING==EXPOSED) THEN
               IWB = WALL_INDEX_BACK(IW)
               IF (BOUNDARY_TYPE(IWB)==SOLID_BOUNDARY) MASS_LOSS(IWB) = MASS_LOSS(IWB) + MASSFLUX(IW,I_FUEL)*DT
            ENDIF
            MASS_LOSS(IW) = MASS_LOSS(IW) + MASSFLUX(IW,I_FUEL)*DT
            OBSTRUCTION(OBST_INDEX_W(IW))%MASS = OBSTRUCTION(OBST_INDEX_W(IW))%MASS - MASSFLUX(IW,I_FUEL)*DT*AW(IW)
         ENDIF CONSUME_FUEL

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
 
         YY_OTHER(1:N_SPECIES) = OMESH(IJKW(9,IW))%YY(IJKW(10,IW),IJKW(11,IW),IJKW(12,IW),1:N_SPECIES)
         IF (CELL_VOLUME_RATIO(IW)<0.5_EB .OR. CELL_VOLUME_RATIO(IW)>2._EB) THEN
            IF (.NOT.SOLID(CELL_INDEX(IIG,JJG,KKG))) YY_W(IW,1:N_SPECIES) = YYP(IIG,JJG,KKG,1:N_SPECIES) 
            SELECT CASE(IOR)
               CASE( 1)
                  IF (UU(II,JJ,KK)>0._EB)   YY_W(IW,1:N_SPECIES) = YY_OTHER(1:N_SPECIES)
               CASE(-1)
                  IF (UU(II-1,JJ,KK)<0._EB) YY_W(IW,1:N_SPECIES) = YY_OTHER(1:N_SPECIES)
               CASE( 2)
                  IF (VV(II,JJ,KK)>0._EB)   YY_W(IW,1:N_SPECIES) = YY_OTHER(1:N_SPECIES)
               CASE(-2)
                  IF (VV(II,JJ-1,KK)<0._EB) YY_W(IW,1:N_SPECIES) = YY_OTHER(1:N_SPECIES)
               CASE( 3)
                  IF (WW(II,JJ,KK)>0._EB)   YY_W(IW,1:N_SPECIES) = YY_OTHER(1:N_SPECIES)
               CASE(-3)
                  IF (WW(II,JJ,KK-1)<0._EB) YY_W(IW,1:N_SPECIES) = YY_OTHER(1:N_SPECIES)
            END SELECT
         ELSE
            YY_W(IW,1:N_SPECIES) = YY_OTHER(1:N_SPECIES)
         ENDIF
         YYP(II,JJ,KK,1:N_SPECIES) = YY_W(IW,1:N_SPECIES)

   END SELECT METHOD_OF_MASS_TRANSFER

   ! Only set species mass fraction in the ghost cell if it is solid
    
   IF (SOLID(CELL_INDEX(II,JJ,KK)) .AND. .NOT.SOLID(CELL_INDEX(IIG,JJG,KKG))) YYP(II,JJ,KK,1:N_SPECIES) = YY_W(IW,1:N_SPECIES)

ENDDO WALL_CELL_LOOP

END SUBROUTINE SPECIES_BC
 

 
SUBROUTINE DENSITY_BC
 
! Compute density at wall from wall temperatures and mass fractions 

USE PHYSICAL_FUNCTIONS, ONLY : GET_MOLECULAR_WEIGHT
REAL(EB) :: R_SUM_DILUENTS_W,Y_SUM_W,Z_SUM_W, Z_2, WFAC
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
   SF => SURFACE(IBC)

! Determine ghost cell value of RSUM=R0*Sum(Y_i/M_i) for non-mixture fraction case

   IF (.NOT.MIXTURE_FRACTION .AND. N_SPECIES>0) THEN
      RSUM_W(IW) = SPECIES(0)%RCON
      DO N=1,N_SPECIES
         WFAC = SPECIES(N)%RCON - SPECIES(0)%RCON
         RSUM_W(IW) = RSUM_W(IW) + WFAC*YY_W(IW,N)
      ENDDO
   ENDIF
 
! Determine ghost cell value of RSUM=R0*Sum(Y_i/M_i) for mixture fraction case

   IF (MIXTURE_FRACTION) THEN
      IF (CO_PRODUCTION) THEN
         Z_2 = YY_W(IW,I_PROG_CO)
      ELSE
         Z_2 = 0._EB
      ENDIF
      Z_SUM_W  =  0._EB
      Y_SUM_W  =  0._EB
      IF (N_SPEC_DILUENTS > 0) R_SUM_DILUENTS_W = 0._EB
      DO N=1,N_SPECIES
         IF (SPECIES(N)%MODE==MIXTURE_FRACTION_SPECIES) THEN
            Z_SUM_W = Z_SUM_W + YY_W(IW,N)
         ELSEIF (SPECIES(N)%MODE==GAS_SPECIES) THEN
            Y_SUM_W = Y_SUM_W + YY_W(IW,N)
            R_SUM_DILUENTS_W = R_SUM_DILUENTS_W + SPECIES(N)%RCON*YY_W(IW,N)
         ENDIF
      ENDDO
      CALL GET_MOLECULAR_WEIGHT(YY_W(IW,I_FUEL),Z_2,YY_W(IW,I_PROG_F),Y_SUM_W,RSUM_W(IW))
      RSUM_W(IW) = R0/RSUM_W(IW)
      IF (N_SPEC_DILUENTS > 0) RSUM_W(IW) = RSUM_W(IW)*(1._EB-Y_SUM_W) + R_SUM_DILUENTS_W
   ENDIF
 
! Compute ghost cell density

   IF (N_SPECIES==0) THEN
      RHO_W(IW) = PBAR_P(KK,PRESSURE_ZONE_WALL(IW))/(SPECIES(0)%RCON*TMP_W(IW))
   ELSE
      RHO_W(IW) = PBAR_P(KK,PRESSURE_ZONE_WALL(IW))/(RSUM_W(IW)*TMP_W(IW))
   ENDIF
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

USE PHYSICAL_FUNCTIONS, ONLY: GET_MOLECULAR_WEIGHT 
USE GEOMETRY_FUNCTIONS
USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
REAL(EB) :: DTMP,QNETF,QDXKF,QDXKB,RR,TMP_G,T,RFACF,RFACB,RFACF2,RFACB2, &
            PPCLAUS,PPSURF,TMP_G_B,RHOWAL,CP_TERM,TSI,DT_BC, &
            DXKF,DXKB,REACTION_RATE,QCONB,DELTA_RHO_S,QRADINB,RFLUX_UP,RFLUX_DOWN,E_WALLB, &
            HVRG,Z_2,Z_F,Y_MF_G,Y_MF_W,YY_S,Y_SUM_W, RSUM_W, RSUM_G, RSUM_S, MFLUX, VOLSUM, &
            DXF, DXB,HTCB,Q_WATER_F,Q_WATER_B,TMP_F_OLD, DX_GRID, RHO_S0, H_R,DT2_BC,TOLERANCE,C_S_ADJUST_UNITS
INTEGER :: IBC,IIG,JJG,KKG,IIB,JJB,KKB,IWB,IC,NWP,I,J,ITMP,NR,IL,NN,NNN,NL,J1,J2,II,JJ,KK,IW,IOR,N,I_OBST
REAL(EB) :: SMALLEST_CELL_SIZE(MAX_LAYERS),THICKNESS,LAYER_THICKNESS_NEW(MAX_LAYERS)
REAL(EB),ALLOCATABLE,DIMENSION(:) :: TMP_W_NEW
INTEGER  :: N_LAYER_CELLS_NEW(MAX_LAYERS), NWP_NEW,I_GRAD,STEPCOUNT

LOGICAL :: POINT_SHRINK, RECOMPUTE,ITERATE
TYPE (WALL_TYPE), POINTER :: WC
TYPE (SURFACE_TYPE), POINTER :: SF
TYPE (MATERIAL_TYPE), POINTER :: ML
 
! Special adjustment of specific heat for steady state applications

C_S_ADJUST_UNITS = 1000._EB/TIME_SHRINK_FACTOR
 
! Loop through the thermally-thick wall cells

WALL_CELL_LOOP: DO IW=1,NWC
 
   IF (BOUNDARY_TYPE(IW)/=SOLID_BOUNDARY .OR. BOUNDARY_TYPE(IW)==POROUS_BOUNDARY) CYCLE WALL_CELL_LOOP
   IBC = IJKW(5,IW)
   SF  => SURFACE(IBC)
   IF (SF%THERMAL_BC_INDEX /= THERMALLY_THICK) CYCLE WALL_CELL_LOOP
   IF (SUM(WMPUA(IW,:))>0._EB .AND. T>TW(IW)) EW(IW) = EW(IW) + SF%E_COEFFICIENT*SUM(WMPUA(IW,:))*DT_BC
   II  = IJKW(1,IW)
   JJ  = IJKW(2,IW)
   KK  = IJKW(3,IW)
   WC  => WALL(IW)
   IOR = IJKW(4,IW)
   IIG = IJKW(6,IW)
   JJG = IJKW(7,IW)
   KKG = IJKW(8,IW)
   IF (WC%BURNAWAY) CYCLE WALL_CELL_LOOP


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
      THICKNESS = SUM(WC%LAYER_THICKNESS)
      CALL GET_WALL_NODE_WEIGHTS(NWP,SF%N_LAYERS,WC%N_LAYER_CELLS,THICKNESS,SF%GEOMETRY, &
         WC%X_S(0:NWP),DX_S(1:NWP),RDX_S(0:NWP+1),RDXN_S(0:NWP),DX_WGT_S(0:NWP),DXF,DXB,LAYER_INDEX)
   ELSE COMPUTE_GRID
      NWP         = SF%N_CELLS
      DXF         = SF%DXF
      DXB         = SF%DXB
      DX_S(1:NWP)          = SF%DX(1:NWP)
      RDX_S(0:NWP+1)       = SF%RDX(0:NWP+1)
      RDXN_S(0:NWP)        = SF%RDXN(0:NWP)
      DX_WGT_S(0:NWP)      = SF%DX_WGT(0:NWP)
      LAYER_INDEX(0:NWP+1) = SF%LAYER_INDEX(0:NWP+1)
   ENDIF COMPUTE_GRID
 
   ! Calculate reaction rates based on the solid phase reactions
 
   Q_S                  = 0._EB

   PYROLYSIS_MATERIAL_IF: IF (SF%PYROLYSIS_MODEL==PYROLYSIS_MATERIAL) THEN

   MFLUX                = MASSFLUX(IW,I_FUEL)
   MASSFLUX(IW,I_FUEL)  = 0._EB
   ACTUAL_BURN_RATE(IW) = 0._EB
   POINT_SHRINK         = .FALSE.
   WC%SHRINKING         = .FALSE.
   IF (SF%SHRINK) X_S_NEW(0:NWP)  = WC%X_S(0:NWP)
   
   POINT_LOOP1: DO I=1,NWP

      IF (SF%SHRINK) THEN
      POINT_SHRINK = .TRUE.
      MATERIAL_LOOP1a: DO N=1,SF%N_MATL
         IF (WC%RHO_S(I,N).EQ. 0._EB) CYCLE MATERIAL_LOOP1a
         ML  => MATERIAL(SF%MATL_INDEX(N))
         IF (ML%N_REACTIONS.EQ.0 .OR. ANY(ML%NU_RESIDUE.GT.0._EB)) THEN
            POINT_SHRINK = .FALSE.
            EXIT MATERIAL_LOOP1a
         ENDIF
      ENDDO MATERIAL_LOOP1a
      ENDIF
         
      RHO_S0 = SF%LAYER_DENSITY(LAYER_INDEX(I))
      VOLSUM = 0._EB

      MATERIAL_LOOP1b: DO N=1,SF%N_MATL
         ML  => MATERIAL(SF%MATL_INDEX(N))

         IF (WC%RHO_S(I,N) .EQ. 0._EB) CYCLE MATERIAL_LOOP1b
         IF (ML%PYROLYSIS_MODEL==PYROLYSIS_LIQUID) THEN
            VOLSUM = 1._EB
            CYCLE MATERIAL_LOOP1b
         ENDIF

         DO J=1,ML%N_REACTIONS
            ! Reaction rate in 1/s
            REACTION_RATE = ML%A(J)*(WC%RHO_S(I,N)/RHO_S0)**ML%N_S(J)*EXP(-ML%E(J)/(R0*WC%TMP_S(I)))
            IF (ML%N_T(J)/=0.) THEN
               DTMP = WC%TMP_S(I)-ML%TMP_IGN(J)
               IF (DTMP > 0._EB) THEN
                  REACTION_RATE = REACTION_RATE * DTMP**ML%N_T(J)
               ELSE
                  REACTION_RATE = 0._EB
               ENDIF
            ENDIF
            ! Reaction rate in kg/(m3s)
            REACTION_RATE = RHO_S0 * REACTION_RATE
            ! Limit reaction rate
            REACTION_RATE = MIN(REACTION_RATE , WC%RHO_S(I,N)/DT_BC)
            MASSFLUX(IW,I_FUEL)  = MASSFLUX(IW,I_FUEL) +  &
               ML%ADJUST_BURN_RATE*ML%NU_FUEL(J)*REACTION_RATE/RDX_S(I)/SF%THICKNESS**I_GRAD
            ACTUAL_BURN_RATE(IW) = ACTUAL_BURN_RATE(IW) +  &
               ML%NU_FUEL(J)*REACTION_RATE/RDX_S(I)/SF%THICKNESS**I_GRAD
            IF (I_WATER>0) MASSFLUX(IW,I_WATER) = MASSFLUX(IW,I_WATER) + &
               ML%NU_WATER(J)*REACTION_RATE/RDX_S(I)/SF%THICKNESS**I_GRAD
            ITMP = MIN(2000,MAX(0,NINT(WC%TMP_S(I)-TMPA)))
            Q_S(I) = Q_S(I) - REACTION_RATE * ML%Q_ARRAY(J,ITMP)
            WC%RHO_S(I,N) = WC%RHO_S(I,N) - DT_BC*REACTION_RATE
            IF (ML%NU_RESIDUE(J) .GT. 0._EB ) THEN
               NNN = SF%RESIDUE_INDEX(N,J)
               WC%RHO_S(I,NNN) = WC%RHO_S(I,NNN) + ML%NU_RESIDUE(J)*DT_BC*REACTION_RATE
            ENDIF
         ENDDO
         VOLSUM = VOLSUM + WC%RHO_S(I,N)/ML%RHO_S
      ENDDO MATERIAL_LOOP1b

      VOLSUM = MIN(VOLSUM,1._EB)
      IF (SF%SHRINK) X_S_NEW(I) = X_S_NEW(I-1)+(WC%X_S(I)-WC%X_S(I-1))*VOLSUM
      IF (POINT_SHRINK) THEN
         IF (VOLSUM<1.0_EB) THEN
            WC%SHRINKING=.TRUE.
            IF (VOLSUM>0.0_EB) THEN
               MATERIAL_LOOP1c: DO N=1,SF%N_MATL
                  WC%RHO_S(I,N) = WC%RHO_S(I,N)/VOLSUM
               ENDDO MATERIAL_LOOP1c
            ENDIF
         ENDIF
      ENDIF
   ENDDO POINT_LOOP1

   ! If the fuel massflux is non-zero, set the iwgnition time

   IF (TW(IW)>T .AND. MASSFLUX(IW,I_FUEL)>0._EB) TW(IW) = T

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

   MFLUX = MAX(0._EB,MFLUX - MASSFLUX(IW,I_FUEL))
   MATERIAL_LOOP2: DO N=1,SF%N_MATL
      ML  => MATERIAL(SF%MATL_INDEX(N))
      IF (ML%PYROLYSIS_MODEL/=PYROLYSIS_LIQUID) CYCLE MATERIAL_LOOP2
      IF (WC%RHO_S(1,N)==0._EB) CYCLE MATERIAL_LOOP2     
      MFLUX = MFLUX/(ML%NU_FUEL(1)+TINY(MFLUX))
      ! gas phase 
      Y_MF_G = YY(IIG,JJG,KKG,I_FUEL)
      IF (CO_PRODUCTION) THEN
         Z_2 = YY(IIG,JJG,KKG,I_PROG_CO)
      ELSE
         Z_2 = 0._EB
      ENDIF
      CALL GET_MOLECULAR_WEIGHT(Y_MF_G,Z_2,YY(IIG,JJG,KKG,I_PROG_F),  Y_SUM(IIG,JJG,KKG),RSUM_G)
      ! wall values
      Y_MF_W = YY_W(IW,I_FUEL)
      IF (CO_PRODUCTION) THEN
         Z_2 = YY_W(IW,I_PROG_CO)
      ELSE
         Z_2 = 0._EB
      ENDIF
      Y_SUM_W = 0._EB
      DO J=1,N_SPECIES
         IF (SPECIES(J)%MODE==GAS_SPECIES) THEN
            Y_SUM_W = Y_SUM_W + YY_W(IW,J)
         ENDIF
      ENDDO            
      CALL GET_MOLECULAR_WEIGHT(Y_MF_W ,Z_2,YY_W(IW,I_PROG_F),Y_SUM_W,RSUM_W)
!!! Check  this at some point
!      YY_S    = 0.5*(Y_MF_W+Y_MF_G)
!      RSUM_S  = 0.5*(RSUM_W+RSUM_G)
      YY_S    = 0.2*Y_MF_W+0.8*Y_MF_G
      RSUM_S  = 0.2*RSUM_W+0.8*RSUM_G
      HVRG    = REACTION(1)%MW_FUEL*ML%H_R(1)/R0
      PPSURF  = MIN(1._EB,YY_S*R0/(REACTION(1)%MW_FUEL*RSUM_S))
      PPCLAUS = MIN(1._EB,EXP(HVRG*(1./ML%TMP_BOIL(1)-1./WC%TMP_S(1))))
      IF (PPSURF/=PPCLAUS) THEN
         IF (MFLUX==0.AND.PPSURF<PPCLAUS) THEN
            MFLUX = 2.E-5*REACTION(1)%MW_FUEL
            TW(IW)=T
         ELSE
            MFLUX = MFLUX*MIN(1.02_EB,MAX(0.98_EB,PPCLAUS/PPSURF))
         ENDIF
      ENDIF
      IF (WC%TMP_S(1)>ML%TMP_BOIL(1)) THEN
         ! Net flux guess for liquid evap
         QNETF = Q_WATER_F + QRADIN(IW) - QRADOUT(IW) + QCONF(IW)
         MFLUX = MAX(MFLUX,1.02*(QNETF - 2.*(ML%TMP_BOIL(1)-WC%TMP_S(1))/DXF/K_S(1))/ML%H_R(1))
      ENDIF
      MFLUX = MIN(MFLUX,THICKNESS*ML%RHO_S/DT_BC)
      ! CYLINDRICAL and SPHERICAL scaling not implemented
      MASSFLUX(IW,I_FUEL) = MASSFLUX(IW,I_FUEL)  + ML%ADJUST_BURN_RATE*ML%NU_FUEL(1)*MFLUX
      ACTUAL_BURN_RATE(IW)= ACTUAL_BURN_RATE(IW) +                     ML%NU_FUEL(1)*MFLUX
      IF (I_WATER>0) MASSFLUX(IW,I_WATER) = MASSFLUX(IW,I_WATER) + ML%NU_WATER(J)*MFLUX
      DX_GRID = DT_BC*MFLUX/ML%RHO_S
      Q_S(1) = Q_S(1) - MFLUX*ML%H_R(1)/DX_S(1)  ! no improvement (in cone test) if used updated RDX 

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
         IF ( (X_S_NEW(I)-X_S_NEW(I-1)) < 0.5 * SMALLEST_CELL_SIZE(LAYER_INDEX(I))) RECOMPUTE = .TRUE.
      ENDDO
      IF (THICKNESS == 0._EB) THEN
         WC%TMP_S(0:NWP+1)    = MAX(TMPMIN,SF%TMP_BACK)
         TMP_F(IW)            = MIN(TMPMAX,MAX(TMPMIN,SF%TMP_BACK))
         TMP_B(IW)            = MIN(TMPMAX,MAX(TMPMIN,SF%TMP_BACK))
         RHOWAL               = 0.5_EB*(RHO(IIG,JJG,KKG)+RHO_W(IW))
         CP_TERM              = MAX(0._EB,-CP_GAMMA*UW(IW)*RHOWAL)
         QCONF(IW)            = HEAT_TRANS_COEF(IW) * (TMP_G - 0.5_EB * (TMP_F(IW) + TMP_F_OLD) )
         TMP_W(IW)            = ( (RDN(IW)*KW(IW)-0.5_EB*CP_TERM)*TMP_G+CP_TERM*TMP_F(IW)-QCONF(IW) ) & 
                                /(0.5_EB*CP_TERM+RDN(IW)*KW(IW))
         TMP_W(IW)            = MAX(TMPMIN,TMP_W(IW))
         MASSFLUX(IW,I_FUEL)  = 0._EB
         ACTUAL_BURN_RATE(IW) = 0._EB
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
         CALL GET_WALL_NODE_WEIGHTS(NWP_NEW,SF%N_LAYERS,N_LAYER_CELLS_NEW,THICKNESS,SF%GEOMETRY, &
            X_S_NEW(0:NWP_NEW),DX_S(1:NWP_NEW),RDX_S(0:NWP_NEW+1),RDXN_S(0:NWP_NEW),&
            DX_WGT_S(0:NWP_NEW),DXF,DXB,LAYER_INDEX)           
         ! Interpolate densities and temperature from old grid to new grid
         CALL GET_INTERPOLATION_WEIGHTS(SF%N_LAYERS,NWP,NWP_NEW,WC%N_LAYER_CELLS,N_LAYER_CELLS_NEW, &
                                    WC%X_S(0:NWP),X_S_NEW(0:NWP_NEW),INT_WGT(1:NWP_NEW,1:NWP))      
         CALL INTERPOLATE_WALL_ARRAY(NWP,NWP_NEW,INT_WGT(1:NWP_NEW,1:NWP),WC%TMP_S(1:NWP))
         WC%TMP_S(NWP_NEW+1) = WC%TMP_S(NWP+1)
         CALL INTERPOLATE_WALL_ARRAY(NWP,NWP_NEW,INT_WGT(1:NWP_NEW,1:NWP),Q_S(1:NWP))
         DO N=1,SF%N_MATL
            ML  => MATERIAL(SF%MATL_INDEX(N))
            CALL INTERPOLATE_WALL_ARRAY(NWP,NWP_NEW,INT_WGT(1:NWP_NEW,1:NWP),WC%RHO_S(1:NWP,N))
         ENDDO
         WC%N_LAYER_CELLS = N_LAYER_CELLS_NEW(1:SF%N_LAYERS)
         NWP = NWP_NEW
         WC%X_S(0:NWP) = X_S_NEW(0:NWP)      ! Note: WC%X_S(NWP+1...) are not set to zero.
      ELSE      
         CALL GET_WALL_NODE_WEIGHTS(NWP,SF%N_LAYERS,N_LAYER_CELLS_NEW,THICKNESS,SF%GEOMETRY, &
            WC%X_S(0:NWP),DX_S(1:NWP),RDX_S(0:NWP+1),RDXN_S(0:NWP),DX_WGT_S(0:NWP),DXF,DXB,LAYER_INDEX)
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
            ITMP   = NINT(WC%TMP_S(I))
            NR     = -NINT(ML%K_S)
            K_S(I) = K_S(I) + WC%RHO_S(I,N)*EVALUATE_RAMP(WC%TMP_S(I),0._EB,NR)/ML%RHO_S
         ENDIF

         IF (ML%C_S>0._EB) THEN
            RHOCBAR(I) = RHOCBAR(I) + WC%RHO_S(I,N)*ML%C_S
         ELSE
            ITMP   = NINT(WC%TMP_S(I))
            NR     = -NINT(ML%C_S)
            RHOCBAR(I) = RHOCBAR(I) + WC%RHO_S(I,N)*EVALUATE_RAMP(WC%TMP_S(I),0._EB,NR)*C_S_ADJUST_UNITS
         ENDIF
         IF (I.EQ.1) E_WALL(IW) = E_WALL(IW) + WC%RHO_S(I,N)*ML%EMISSIVITY/ML%RHO_S
         RHO_S(I)   = RHO_S(I) + WC%RHO_S(I,N)
         
      ENDDO MATERIAL_LOOP3
      
      K_S(I) = K_S(I)/VOLSUM
      IF (I.EQ.1) E_WALL(IW) = E_WALL(IW)/VOLSUM
   ENDDO POINT_LOOP3
  
   K_S(0)     = K_S(1)     ! Calculate average K_S between at grid cell boundaries. Store result in K_S
   K_S(NWP+1) = K_S(NWP)
   DO I=1,NWP-1 
      K_S(I)  = 1.0_EB / ( DX_WGT_S(I)/K_S(I) + (1.-DX_WGT_S(I))/K_S(I+1) )
   ENDDO

   !RCP_W for part
   IF (N_EVAP_INDICIES>0) THEN
      RCP_W(IW) = RHOCBAR(1)*DX_S(1)
   ENDIF


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
         KAPPA_S(I) = 2.*KAPPA_S(I)*DX_S(I)/VOLSUM    ! kappa = 2*dx*kappa
      ENDDO
      ! solution inwards
      RFLUX_UP = QRADIN(IW) + (1.-E_WALL(IW))*QRADOUT(IW)
      DO I=1,NWP
         RFLUX_DOWN =  ( RFLUX_UP + KAPPA_S(I)*SIGMA*WC%TMP_S(I)**4 ) / (1. + KAPPA_S(I))
         Q_S(I) = Q_S(I) + (RFLUX_UP - RFLUX_DOWN)/DX_S(I)
         RFLUX_UP = RFLUX_DOWN
      ENDDO
      IF (SF%BACKING==EXPOSED) THEN
         IF (BOUNDARY_TYPE(IWB)==SOLID_BOUNDARY) QRADOUT(IWB) = E_WALLB*RFLUX_UP
      ENDIF
      ! solution outwards   
      RFLUX_UP = QRADINB + (1.-E_WALLB)*RFLUX_UP
      DO I=NWP,1,-1
         RFLUX_DOWN =  ( RFLUX_UP + KAPPA_S(I)*SIGMA*WC%TMP_S(I)**4 ) / (1. + KAPPA_S(I))
         Q_S(I) = Q_S(I) + (RFLUX_UP - RFLUX_DOWN)/DX_S(I)
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
            BBS(I) = -0.5_EB*DT_BC*K_S(I-1)*RDXN_S(I-1)*RDX_S(I)/RHOCBAR(I)
            AAS(I) = -0.5_EB*DT_BC*K_S(I)  *RDXN_S(I)  *RDX_S(I)/RHOCBAR(I)
         ENDDO
         DDS(1:NWP) = 1._EB    - AAS(1:NWP) - BBS(1:NWP)
         DO I=1,NWP
            CCS(I) = TMP_W_NEW(I) - AAS(I)*(TMP_W_NEW(I+1)-TMP_W_NEW(I)) + BBS(I)*(TMP_W_NEW(I)-TMP_W_NEW(I-1)) &
                     + DT2_BC*Q_S(I)/RHOCBAR(I)
         ENDDO
         IF (SF%INTERNAL_RADIATION) THEN
            RFACF = 0.25_EB*HEAT_TRANS_COEF(IW)
            RFACB = 0.25_EB*HTCB
         ELSE
            RFACF = 0.25_EB*HEAT_TRANS_COEF(IW)+2._EB*E_WALL(IW)*SIGMA*TMP_F(IW)**3
            RFACB = 0.25_EB*HTCB               +2._EB*E_WALLB*   SIGMA*TMP_B(IW)**3
         ENDIF
         RFACF2 = (DXKF-RFACF)/(DXKF+RFACF)
         RFACB2 = (DXKB-RFACB)/(DXKB+RFACB)
         IF (SF%INTERNAL_RADIATION) THEN
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
      !   WC%TMP_S(1:NWP) = MAX(TMPMIN,CCS(1:NWP))
      !   WC%TMP_S(0)     = MAX(TMPMIN,WC%TMP_S(1)  *RFACF2+QDXKF)
      !   WC%TMP_S(NWP+1) = MAX(TMPMIN,WC%TMP_S(NWP)*RFACB2+QDXKB)
         TMP_W_NEW(1:NWP) = MAX(TMPMIN,CCS(1:NWP))
         TMP_W_NEW(0)     = MAX(TMPMIN,TMP_W_NEW(1)  *RFACF2+QDXKF)
         TMP_W_NEW(NWP+1) = MAX(TMPMIN,TMP_W_NEW(NWP)*RFACB2+QDXKB)
         IF (STEPCOUNT==1) THEN
            TOLERANCE        = MAXVAL(ABS((TMP_W_NEW-WC%TMP_S(0:NWP+1))/WC%TMP_S(0:NWP+1)))
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
         CASE(1:2)
            H_NATURAL = HCV*ABS(DELTA_TMP)**ONTH
         CASE(3)
            H_NATURAL = HCH*ABS(DELTA_TMP)**ONTH
      END SELECT
      H_FORCED = 0._EB
   ELSE
      SELECT CASE(ABS(IOR))
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
   ITMP = 0.1_EB*TMP_G
   H_DNS = SPECIES(0)%MU(ITMP)*CP_GAMMA*RPR*2._EB*RDN(IW)
   HEAT_TRANSFER_COEFFICIENT = MAX(H_DNS,H_FORCED,H_NATURAL)
ENDIF

END FUNCTION HEAT_TRANSFER_COEFFICIENT

END SUBROUTINE WALL_BC
 
END MODULE WALL_ROUTINES

