MODULE PART
 
USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_POINTERS
USE TRAN
USE MEMORY_FUNCTIONS, ONLY: ALLOCATE_STORAGE
USE COMP_FUNCTIONS, ONLY : SECOND 
IMPLICIT NONE
PRIVATE
PUBLIC INSERT_PARTICLES,UPDATE_PARTICLES,REMOVE_PARTICLES,INITIALIZE_PARTICLES,PARTICLE_MOMENTUM_TRANSFER,GET_REV_part
CHARACTER(255), PARAMETER :: partid='$Id$'
CHARACTER(255), PARAMETER :: partrev='$Revision$'
CHARACTER(255), PARAMETER :: partdate='$Date$'

CONTAINS
 

SUBROUTINE INITIALIZE_PARTICLES(NM)

! Insert PARTICLEs into the domain at the start of calculation

USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
USE PHYSICAL_FUNCTIONS, ONLY : PARTICLE_SIZE_DISTRIBUTION 
 
REAL(EB) :: LL,UL,BIN_SIZE,TNOW,DD,DI
INTEGER  :: I,J,IL,IU,ILPC
INTEGER, INTENT(IN) :: NM
TYPE (LAGRANGIAN_PARTICLE_CLASS_TYPE), POINTER :: LPC=>NULL()
TYPE (RAMPS_TYPE), POINTER :: RM=>NULL()

IF (N_LAGRANGIAN_CLASSES==0) RETURN ! Don't waste time if no particles
IF (EVACUATION_ONLY(NM)) RETURN     ! Don't waste time if an evac mesh
 
TNOW=SECOND()
CALL POINT_TO_MESH(NM)

PART_CLASS_LOOP: DO ILPC=1,N_LAGRANGIAN_CLASSES

   LPC=>LAGRANGIAN_PARTICLE_CLASS(ILPC)
   
   ! If particles/PARTICLEs have a size distribution, initialize here
 
   IF_SIZE_DISTRIBUTION: IF (.NOT.LPC%MONODISPERSE .AND. (LPC%DIAMETER > 0._EB) .OR. LPC%CNF_RAMP_INDEX>=0) THEN
      IF (LPC%CNF_RAMP_INDEX<0) THEN
         CALL PARTICLE_SIZE_DISTRIBUTION(LPC%DIAMETER,LPC%R_CNF(:),LPC%CNF(:),LPC%CVF(:),NDC,LPC%GAMMA,LPC%SIGMA,LPC%DISTRIBUTION)
      ELSE
         RM=>RAMPS(LPC%CNF_RAMP_INDEX)
         LPC%DIAMETER = RM%SPAN*0.5_EB ! This isn't actually used, but some routines might expect it to be set.
         DD           = RM%SPAN/NDC
         LPC%R_CNF(0) = RM%T_MIN
         LPC%CNF(0)   = 0._EB
         LPC%CVF(0)   = 0._EB
         DO I=1,NDC
            DI           = RM%T_MIN + I*DD
            LPC%R_CNF(I) = DI
            LPC%CNF(I)   = EVALUATE_RAMP(DI,0._EB,LPC%CNF_RAMP_INDEX) 
            LPC%CVF(I)   = LPC%CVF(I-1) + (DI-0.5*DD)**3*(LPC%CNF(I)-LPC%CNF(I-1))
         ENDDO
         LPC%R_CNF = 1.E-6_EB*0.5_EB*LPC%R_CNF  ! Convert diameter in microns to radius in meters.
         LPC%CNF   = LPC%CNF/LPC%CNF(NDC)
         LPC%CVF   = LPC%CVF/LPC%CVF(NDC)
      ENDIF
      BIN_SIZE = (LPC%R_CNF(NDC)-LPC%R_CNF(0))/REAL(LPC%N_STRATA,EB)
      STRATIFY: DO I=1,LPC%N_STRATA
         LL = LPC%R_CNF(0) + (I-1)*BIN_SIZE
         UL = LPC%R_CNF(0) +  I   *BIN_SIZE
         IL = 1
         LL_LOOP: DO J=1,NDC
            IF (LPC%R_CNF(J)>LL) THEN
               IL = J-1 
               LPC%IL_CNF(I) = J-1
               EXIT LL_LOOP
            ENDIF
         ENDDO LL_LOOP
         IU = NDC
         UL_LOOP: DO J=NDC,1,-1
            IF (LPC%R_CNF(J)<=UL) THEN
               IU = J 
               LPC%IU_CNF(I) = J
               EXIT UL_LOOP
            ENDIF
         ENDDO UL_LOOP
         LPC%W_CNF(I) = LPC%CNF(IU) - LPC%CNF(IL)
      ENDDO STRATIFY
   ENDIF IF_SIZE_DISTRIBUTION

   ! If pacticles/PARTICLEs can break up, compute normalized (median = 1) size distribution for child PARTICLEs

   IF (LPC%BREAKUP .AND. .NOT.LPC%MONODISPERSE) THEN
      IF(LPC%BREAKUP_CNF_RAMP_INDEX<0) THEN
         CALL PARTICLE_SIZE_DISTRIBUTION(1._EB,LPC%BREAKUP_R_CNF(:),LPC%BREAKUP_CNF(:),LPC%BREAKUP_CVF(:),NDC, &
                                        LPC%BREAKUP_GAMMA,LPC%BREAKUP_SIGMA,LPC%BREAKUP_DISTRIBUTION)
      ELSE
         RM=>RAMPS(LPC%BREAKUP_CNF_RAMP_INDEX)
         DD=RM%SPAN/NDC
         DO I=1,NDC
               DI=RM%T_MIN+(I-0.5_EB)*DD
               LPC%BREAKUP_R_CNF(I) = 1E-6_EB*0.5_EB*DI
               LPC%BREAKUP_CNF(I)   = EVALUATE_RAMP(DI,0._EB,LPC%BREAKUP_CNF_RAMP_INDEX) 
         ENDDO
         LPC%BREAKUP_CNF=LPC%BREAKUP_CNF/LPC%BREAKUP_CNF(NDC)
      ENDIF
   ENDIF

ENDDO PART_CLASS_LOOP

TUSED(8,NM)=TUSED(8,NM)+SECOND()-TNOW

END SUBROUTINE INITIALIZE_PARTICLES
 

SUBROUTINE INSERT_PARTICLES(T,NM)

! Insert sprinkler PARTICLEs and lagrangian particles into the domain every time step

USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP,RANDOM_CHOICE
USE GEOMETRY_FUNCTIONS, ONLY: RANDOM_RECTANGLE,RANDOM_CONE,RANDOM_RING
USE TRAN, ONLY: GET_IJK
USE DEVICE_VARIABLES
USE CONTROL_VARIABLES
REAL(EB), INTENT(IN) :: T
INTEGER, INTENT(IN) :: NM
REAL     :: RN,RN2
REAL(EB) :: PHI_RN,FLOW_RATE,THETA_RN,SPHI,CPHI,MASS_SUM,D_PRES_FACTOR, &
            STHETA,CTHETA,PWT0,PARTICLE_SPEED,SHIFT1,SHIFT2,XTMP,YTMP,ZTMP,VLEN, &
            TRIGT1,TRIGT2,TNOW,TSI,PIPE_PRESSURE,X1,X2,Y1,Y2,Z1,Z2, &
            ETA,ETA_MAX,ETA_MIN,XI,YJ,ZK
REAL(EB), PARAMETER :: VENT_OFFSET=0.1
INTEGER :: IP,KS,II,JJ,KK,IC,IL,IU,ILPC,DROP_SUM,IIG,JJG,KKG,IW,IOR,STRATUM,IB
INTEGER :: N,N_INSERT,ILAT,NEW_LP_INDEX
INTEGER, ALLOCATABLE, DIMENSION(:) :: LP_INDEX_LOOKUP
LOGICAL :: INSERT_ANOTHER_BATCH
TYPE (PROPERTY_TYPE), POINTER :: PY=>NULL()
TYPE (TABLES_TYPE), POINTER :: TA=>NULL()
TYPE (DEVICE_TYPE), POINTER :: DV=>NULL()
TYPE (SURFACE_TYPE), POINTER :: SF=>NULL()
TYPE (LAGRANGIAN_PARTICLE_TYPE), POINTER :: LP=>NULL()
TYPE (LAGRANGIAN_PARTICLE_CLASS_TYPE), POINTER :: LPC=>NULL()
TYPE (INITIALIZATION_TYPE), POINTER :: IN=>NULL()
 
IF (EVACUATION_ONLY(NM)) RETURN     ! Don't waste time if an evac mesh
IF (N_LAGRANGIAN_CLASSES==0) RETURN ! Don't waste time if no particles

TNOW=SECOND()
CALL POINT_TO_MESH(NM)

OVERALL_INSERT_LOOP: DO  

   INSERT_ANOTHER_BATCH = .FALSE.

   CALL DEVICE_PARTICLE_INSERT
   CALL WALL_PARTICLE_INSERT
   CALL VOLUME_PARTICLE_INSERT

   ! Reset particle/PARTICLE insertion clocks

   DO N=1,N_INIT
      IN => INITIALIZATION(N)
      IF (IN%SINGLE_INSERTION) CYCLE      
      IF (T >= IN%PARTICLE_INSERT_CLOCK(NM)) IN%PARTICLE_INSERT_CLOCK(NM) = IN%PARTICLE_INSERT_CLOCK(NM) + IN%DT_INSERT
      IF (T >= IN%PARTICLE_INSERT_CLOCK(NM)) INSERT_ANOTHER_BATCH = .TRUE.
   ENDDO

   DO N=1,N_SURF
      SF => SURFACE(N)
      ILPC =  SF%PART_INDEX
      IF (ILPC < 1)    CYCLE
      IF (T >= SF%PARTICLE_INSERT_CLOCK(NM)) SF%PARTICLE_INSERT_CLOCK(NM) = SF%PARTICLE_INSERT_CLOCK(NM) + SF%DT_INSERT
      IF (T >= SF%PARTICLE_INSERT_CLOCK(NM)) INSERT_ANOTHER_BATCH = .TRUE.
   ENDDO

   IF (.NOT.INSERT_ANOTHER_BATCH) EXIT OVERALL_INSERT_LOOP

ENDDO OVERALL_INSERT_LOOP

! Compute initial particle CFL

PART_UVWMAX = 0._EB
IF (PARTICLE_CFL) THEN
   DO IP=1,NLP  
      LP => LAGRANGIAN_PARTICLE(IP)
      LPC => LAGRANGIAN_PARTICLE_CLASS(LP%CLASS_INDEX)
      IF (LPC%MASSLESS_TRACER .OR. LPC%MASSLESS_TARGET) CYCLE
      CALL GET_IJK(LP%X,LP%Y,LP%Z,NM,XI,YJ,ZK,II,JJ,KK)
      PART_UVWMAX = MAX(PART_UVWMAX,MAX(ABS(LP%U)*RDX(II),ABS(LP%V)*RDY(JJ),ABS(LP%W)*RDZ(KK)))
   ENDDO
ENDIF

TUSED(8,NM)=TUSED(8,NM)+SECOND()-TNOW

CONTAINS


SUBROUTINE DEVICE_PARTICLE_INSERT

INTEGER :: I,OI

! Count active sprinklers and nozzles

N_OPEN_NOZZLES = 0
N_ACTUATED_SPRINKLERS = 0

COUNT_OPEN_NOZZLES_LOOP: DO KS=1,N_DEVC ! Loop over all devices, but look for sprinklers or nozzles
   DV => DEVICE(KS)
   PY => PROPERTY(DV%PROP_INDEX)
   IF (.NOT. DV%CURRENT_STATE) CYCLE COUNT_OPEN_NOZZLES_LOOP
   IF (PY%PART_ID == 'null')   CYCLE COUNT_OPEN_NOZZLES_LOOP
   N_OPEN_NOZZLES = N_OPEN_NOZZLES + 1
   IF (PY%QUANTITY=='SPRINKLER LINK TEMPERATURE') N_ACTUATED_SPRINKLERS = N_ACTUATED_SPRINKLERS + 1
ENDDO COUNT_OPEN_NOZZLES_LOOP

! Loop over all devices, but look for sprinklers or nozzles
   
SPRINKLER_INSERT_LOOP: DO KS=1,N_DEVC  

   DV => DEVICE(KS)
   PY => PROPERTY(DV%PROP_INDEX)
   IF (PY%PART_ID == 'null')   CYCLE SPRINKLER_INSERT_LOOP
   IF (DV%MESH/=NM)            CYCLE SPRINKLER_INSERT_LOOP
   IF (.NOT. DV%CURRENT_STATE) CYCLE SPRINKLER_INSERT_LOOP
   LPC=>LAGRANGIAN_PARTICLE_CLASS(PY%PART_INDEX)

   IF (ABS(DV%T_CHANGE-T)<=TWO_EPSILON_EB) THEN 
      DV%T = T
      CYCLE SPRINKLER_INSERT_LOOP
   ENDIF

   N_INSERT = NINT(REAL(PY%PARTICLES_PER_SECOND,EB)*(T-DV%T))

   IF (N_INSERT<1) CYCLE SPRINKLER_INSERT_LOOP

   ! Determine sprinkler/nozzle flow rate

   IF (ABS(T_BEGIN-DV%T_CHANGE)<=TWO_EPSILON_EB .AND. PY%FLOW_RAMP_INDEX>=1) THEN
      TSI = T
   ELSE
      TSI = T - DV%T_CHANGE
   ENDIF

   IF (PY%MASS_FLOW_RATE > 0._EB) THEN
      FLOW_RATE = PY%MASS_FLOW_RATE ! kg/s 
      PIPE_PRESSURE = 1.0_EB
      D_PRES_FACTOR = 1.0_EB
   ELSE
      IF (PY%PRESSURE_RAMP_INDEX>0) THEN
         PIPE_PRESSURE = EVALUATE_RAMP(REAL(DEVC_PIPE_OPERATING(DV%PIPE_INDEX),EB),0._EB,PY%PRESSURE_RAMP_INDEX)
         D_PRES_FACTOR = (PY%OPERATING_PRESSURE/PIPE_PRESSURE)**(1._EB/3._EB)
         FLOW_RATE = PY%K_FACTOR*SQRT(PIPE_PRESSURE)
      ELSE
         PIPE_PRESSURE = PY%OPERATING_PRESSURE
         D_PRES_FACTOR = 1.0_EB
         FLOW_RATE = PY%FLOW_RATE
      ENDIF 
      FLOW_RATE = FLOW_RATE*(LPC%DENSITY/1000._EB)/60._EB  ! kg/s     
   ENDIF
   
   FLOW_RATE = EVALUATE_RAMP(TSI,PY%FLOW_TAU,PY%FLOW_RAMP_INDEX)*FLOW_RATE ! kg/s 

   IF (FLOW_RATE <= 0._EB) THEN
      DV%T = T
      CYCLE SPRINKLER_INSERT_LOOP
   ENDIF

   ! Direction initialization stuff

   OI = DV%ORIENTATION_INDEX
   TRIGT1 = ACOS(-ORIENTATION_VECTOR(3,OI))
   IF (ABS(ORIENTATION_VECTOR(2,OI))<=TWO_EPSILON_EB) THEN
      TRIGT2 = ACOS(1._EB)
   ELSE
      TRIGT2 = ACOS(ABS(ORIENTATION_VECTOR(1,OI))/SQRT(ORIENTATION_VECTOR(1,OI)**2+ORIENTATION_VECTOR(2,OI)**2))
   ENDIF
   
   ! PARTICLE insertion loop
   
   MASS_SUM = 0._EB
   DROP_SUM = 0
   
   ALLOCATE(LP_INDEX_LOOKUP(N_INSERT))
   LP_INDEX_LOOKUP = 0

   PARTICLE_INSERT_LOOP: DO I=1,N_INSERT

      ! Insert a single particle

      IF (NLP+1>MAXIMUM_PARTICLES) THEN
         CALL REMOVE_OLDEST_PARTICLE(NM,PY%PART_INDEX,NLP,NEW_LP_INDEX)
         IF (I>1) LP_INDEX_LOOKUP(I-1) = NEW_LP_INDEX
      ELSE
         NLP = NLP+1
      ENDIF
      LP_INDEX_LOOKUP(I) = NLP
      
      PARTICLE_TAG = PARTICLE_TAG + NMESHES
      CALL ALLOCATE_STORAGE(NM,LAGRANGIAN_PARTICLE_CLASS(PY%PART_INDEX)%SURF_INDEX,LPC_INDEX=PY%PART_INDEX,&
                            LP_INDEX=NLP,TAG=PARTICLE_TAG,NEW_TAG=.TRUE.)
      LAGRANGIAN_PARTICLE => MESHES(NM)%LAGRANGIAN_PARTICLE
      LP=>MESHES(NM)%LAGRANGIAN_PARTICLE(NLP)
   
      ! Set PARTICLE properties
   
      LP%T_INSERT = T                     
      IF (MOD(NLP,LPC%SAMPLING)==0) LP%SHOW = .TRUE.    

      ! Randomly choose particle direction angles, theta and phi
   
      CHOOSE_COORDS: DO
         PICK_PATTERN: IF(PY%SPRAY_PATTERN_INDEX>0) THEN !Use spray pattern table
            TA => TABLES(PY%SPRAY_PATTERN_INDEX)
            CALL RANDOM_NUMBER(RN)
            FIND_ROW: DO II=1,TA%NUMBER_ROWS
               IF (REAL(RN,EB)>PY%TABLE_ROW(II)) CYCLE FIND_ROW
               EXIT FIND_ROW
            END DO FIND_ROW
            CALL RANDOM_NUMBER(RN)
            !THETA_RN = TA%TABLE_DATA(II,1) + RN*(TA%TABLE_DATA(II,2)-TA%TABLE_DATA(II,1))
            ETA_MAX=0.5_EB*(COS(TA%TABLE_DATA(II,1))+1._EB)
            ETA_MIN=0.5_EB*(COS(TA%TABLE_DATA(II,2))+1._EB)
            ETA=ETA_MIN+(ETA_MAX-ETA_MIN)*REAL(RN,EB)
            THETA_RN=ACOS(2._EB*ETA-1._EB)
            CALL RANDOM_NUMBER(RN)
            PHI_RN = TA%TABLE_DATA(II,3) + REAL(RN,EB)*(TA%TABLE_DATA(II,4)-TA%TABLE_DATA(II,3))
            
            IF (PY%PRESSURE_RAMP_INDEX>0) THEN
               PARTICLE_SPEED = PY%V_FACTOR(II)*SQRT(PIPE_PRESSURE)
            ELSE 
               PARTICLE_SPEED = TA%TABLE_DATA(II,5)
            ENDIF
         ELSE PICK_PATTERN !Use conical spray
            CALL RANDOM_CHOICE(PY%SPRAY_LON_CDF(:),PY%SPRAY_LON,NDC2,PHI_RN)
            ILAT=MINLOC(ABS(PY%SPRAY_LON-PHI_RN),1)-1
            CALL RANDOM_CHOICE(PY%SPRAY_LAT_CDF(:,ILAT),PY%SPRAY_LAT,NDC2,THETA_RN)
            IF (PY%PRESSURE_RAMP_INDEX>0) THEN
               PARTICLE_SPEED = PY%V_FACTOR(1)*SQRT(PIPE_PRESSURE)
            ELSE
               PARTICLE_SPEED = PY%PARTICLE_VELOCITY
            ENDIF
         ENDIF PICK_PATTERN
         PHI_RN = PHI_RN + DV%ROTATION  ! Adjust for rotation of head by rotating about z-axis

         !  Adjust for tilt of sprinkler pipe
         SPHI   = SIN(PHI_RN)
         CPHI   = COS(PHI_RN)         
         STHETA = SIN(THETA_RN)
         CTHETA = COS(THETA_RN)
         XTMP   = STHETA*CPHI
         YTMP   = STHETA*SPHI
         ZTMP   = -CTHETA
   
         ! First rotate about y-axis away from x-axis
   
         VLEN   = SQRT(XTMP**2+ZTMP**2)
         SHIFT1 = ACOS(ABS(XTMP)/VLEN)
            SELECT CASE(INT(SIGN(1._EB,ZTMP)))
            CASE (-1)
               IF (XTMP<0) SHIFT1 = PI-SHIFT1
            CASE ( 1)
            SELECT CASE(INT(SIGN(1._EB,XTMP)))
               CASE (-1)
                  SHIFT1 = SHIFT1+PI
               CASE ( 1)
                  SHIFT1 = TWOPI - SHIFT1
            END SELECT
         END SELECT
      
         SHIFT1 = SHIFT1 + TRIGT1
         XTMP = VLEN * COS(SHIFT1)
         ZTMP = -VLEN * SIN(SHIFT1)
   
         ! Second rotate about z-axis away from x-axis
   
         VLEN   = SQRT(XTMP**2+YTMP**2)
         SHIFT1 = ACOS(ABS(XTMP)/VLEN)
         SELECT CASE(INT(SIGN(1._EB,YTMP)))
            CASE ( 1)
               IF (XTMP<0) SHIFT1 = PI-SHIFT1
            CASE (-1)
            SELECT CASE(INT(SIGN(1._EB,XTMP)))
               CASE (-1)
                  SHIFT1 = SHIFT1+PI
               CASE ( 1) 
                  SHIFT1 = TWOPI - SHIFT1
            END SELECT
         END SELECT
   
         SHIFT2 = TRIGT2
         OI = DV%ORIENTATION_INDEX
         SELECT CASE(INT(SIGN(1._EB,ORIENTATION_VECTOR(1,OI))))
            CASE (-1)
               IF (ORIENTATION_VECTOR(2,OI)>0) SHIFT2 = TWOPI - SHIFT2
            CASE ( 1)
            SELECT CASE(INT(SIGN(1._EB,ORIENTATION_VECTOR(2,OI))))
               CASE (-1) 
                  SHIFT2 = PI-SHIFT2
               CASE ( 1)
                  SHIFT2 = SHIFT2+ PI
            END SELECT
         END SELECT
         SHIFT1=SHIFT1+SHIFT2
         XTMP = VLEN * COS(SHIFT1)
         YTMP = VLEN * SIN(SHIFT1)
   
         ! Compute initial position and velocity of PARTICLEs
   
         LP%U = PARTICLE_SPEED*XTMP
         LP%V = PARTICLE_SPEED*YTMP
         LP%W = PARTICLE_SPEED*ZTMP
         LP%X = DV%X + PY%OFFSET*XTMP
         LP%Y = DV%Y + PY%OFFSET*YTMP
         LP%Z = DV%Z + PY%OFFSET*ZTMP
         IF (TWO_D) THEN
            LP%V = 0._EB
            LP%Y = DV%Y
         ENDIF
         IF (LP%X<=XS .OR. LP%X>=XF) CYCLE CHOOSE_COORDS
         IF (LP%Y<=YS .OR. LP%Y>=YF) CYCLE CHOOSE_COORDS
         IF (LP%Z<=ZS .OR. LP%Z>=ZF) CYCLE CHOOSE_COORDS
         CALL GET_IJK(LP%X,LP%Y,LP%Z,NM,XI,YJ,ZK,II,JJ,KK)
         IC = CELL_INDEX(II,JJ,KK)
         LP%ONE_D%IIG = II
         LP%ONE_D%JJG = JJ
         LP%ONE_D%KKG = KK
         IF (.NOT.SOLID(IC)) EXIT CHOOSE_COORDS
      ENDDO CHOOSE_COORDS

      ! Randomly choose PARTICLE size according to Cumulative Distribution Function (CDF)

      CALL MAKE_PARTICLE

      LP => LAGRANGIAN_PARTICLE(NLP)
      SF => SURFACE(LPC%SURF_INDEX)

      ! Adjust particle size to account for pressure dependence of nozzle

      IF (LPC%LIQUID_DROPLET) THEN
         LP%ONE_D%X(1) = LP%ONE_D%X(1)*D_PRES_FACTOR
         LP%ONE_D%LAYER_THICKNESS(1) = LP%ONE_D%X(1)
         LP%MASS = LP%MASS*D_PRES_FACTOR**3
      ENDIF

      ! Sum up mass of liquid being introduced

      MASS_SUM = MASS_SUM + LP%PWT*LP%MASS
      DROP_SUM = DROP_SUM + 1
   ENDDO PARTICLE_INSERT_LOOP

   ! Compute weighting factor for the PARTICLEs just inserted
   
   IF (DROP_SUM > 0) THEN
      PWT0 = FLOW_RATE*(T-DV%T)/MASS_SUM
      DO I=1,N_INSERT
         N = LP_INDEX_LOOKUP(I)
         LAGRANGIAN_PARTICLE(N)%PWT = LAGRANGIAN_PARTICLE(N)%PWT*PWT0      
      ENDDO
   ENDIF

   DEALLOCATE(LP_INDEX_LOOKUP)

   ! Indicate that PARTICLEs from this device have been inserted at this time T

   DV%T = T 
   
ENDDO SPRINKLER_INSERT_LOOP

END SUBROUTINE DEVICE_PARTICLE_INSERT


SUBROUTINE WALL_PARTICLE_INSERT

TYPE(WALL_TYPE), POINTER :: WC=>NULL()
INTEGER :: I

! Loop through all boundary cells and insert particles if appropriate
    
WALL_INSERT_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS

   WC => WALL(IW)
   IF (T < WC%ONE_D%T_IGN)                     CYCLE WALL_INSERT_LOOP
   IF (WALL(IW)%BOUNDARY_TYPE/=SOLID_BOUNDARY) CYCLE WALL_INSERT_LOOP
   SF  => SURFACE(WC%SURF_INDEX)
   ILPC =  SF%PART_INDEX
   IF (ILPC < 1)    CYCLE WALL_INSERT_LOOP
   LPC  => LAGRANGIAN_PARTICLE_CLASS(ILPC)
   IF (LPC%DEVC_INDEX>0) THEN
      IF (.NOT.DEVICE(LPC%DEVC_INDEX)%CURRENT_STATE) CYCLE WALL_INSERT_LOOP
   ENDIF
   IF (LPC%CTRL_INDEX>0) THEN
      IF (.NOT.CONTROL(LPC%CTRL_INDEX)%CURRENT_STATE) CYCLE WALL_INSERT_LOOP
   ENDIF
   IF (T < SF%PARTICLE_INSERT_CLOCK(NM)) CYCLE WALL_INSERT_LOOP
   IF (WC%ONE_D%UW >= -0.0001_EB) CYCLE WALL_INSERT_LOOP
   
   II = WC%ONE_D%II
   JJ = WC%ONE_D%JJ
   KK = WC%ONE_D%KK
   IC = CELL_INDEX(II,JJ,KK)
   IF (.NOT.SOLID(IC)) CYCLE WALL_INSERT_LOOP
   
   IIG = WC%ONE_D%IIG
   JJG = WC%ONE_D%JJG
   KKG = WC%ONE_D%KKG
   IF (NM > 1) THEN
      IF (INTERPOLATED_MESH(IIG,JJG,KKG) > 0) CYCLE WALL_INSERT_LOOP
   ENDIF

   ! Loop over all particles for the IW-th cell

   IOR = WC%ONE_D%IOR
   MASS_SUM = 0._EB

   ALLOCATE(LP_INDEX_LOOKUP(SF%NPPC))
   LP_INDEX_LOOKUP = 0
   PARTICLE_INSERT_LOOP2: DO I=1,SF%NPPC
   
      ! Insert a single droplet at wall cell IW

      IF (NLP+1>MAXIMUM_PARTICLES) THEN
         CALL REMOVE_OLDEST_PARTICLE(NM,ILPC,NLP,NEW_LP_INDEX)
         IF (I>1) LP_INDEX_LOOKUP(I-1)=NEW_LP_INDEX
      ELSE
         NLP = NLP+1
      ENDIF
      LP_INDEX_LOOKUP(I) = NLP
      
      PARTICLE_TAG = PARTICLE_TAG + NMESHES
      CALL ALLOCATE_STORAGE(NM,LPC%SURF_INDEX,LPC_INDEX=ILPC,LP_INDEX=NLP,TAG=PARTICLE_TAG,NEW_TAG=.TRUE.)
      LP=>MESHES(NM)%LAGRANGIAN_PARTICLE(NLP)
      LAGRANGIAN_PARTICLE => MESHES(NM)%LAGRANGIAN_PARTICLE

      ! Assign particle position on the cell face

      CALL RANDOM_NUMBER(RN)
      CALL RANDOM_NUMBER(RN2)

      SELECT CASE (ABS(IOR))
         CASE(1)
            IF (IOR== 1) LP%X = X(II)   + VENT_OFFSET*DX(II+1)
            IF (IOR==-1) LP%X = X(II-1) - VENT_OFFSET*DX(II-1)
            LP%Y = Y(JJ-1) + DY(JJ)*REAL(RN,EB)
            LP%Z = Z(KK-1) + DZ(KK)*REAL(RN2,EB)
         CASE(2)
            IF (IOR== 2) LP%Y = Y(JJ)   + VENT_OFFSET*DY(JJ+1)
            IF (IOR==-2) LP%Y = Y(JJ-1) - VENT_OFFSET*DY(JJ-1)
            LP%X = X(II-1) + DX(II)*REAL(RN,EB)
            LP%Z = Z(KK-1) + DZ(KK)*REAL(RN2,EB)
         CASE(3)
            IF (IOR== 3) LP%Z = Z(KK)   + VENT_OFFSET*DZ(KK+1)
            IF (IOR==-3) LP%Z = Z(KK-1) - VENT_OFFSET*DZ(KK-1)
            LP%X = X(II-1) + DX(II)*REAL(RN,EB)
            LP%Y = Y(JJ-1) + DY(JJ)*REAL(RN2,EB)
      END SELECT

      ! Give particles an initial velocity 
   
      SELECT CASE(IOR) 
         CASE( 1) 
            LP%U = -WALL(IW)%ONE_D%UW
            LP%V = SF%VEL_T(1)
            LP%W = SF%VEL_T(2)
         CASE(-1) 
            LP%U =  WALL(IW)%ONE_D%UW
            LP%V = SF%VEL_T(1)
            LP%W = SF%VEL_T(2)
         CASE( 2)
            LP%U = SF%VEL_T(2)
            LP%V = -WALL(IW)%ONE_D%UW
            LP%W = SF%VEL_T(1)
         CASE(-2)
            LP%U = SF%VEL_T(2)
            LP%V =  WALL(IW)%ONE_D%UW
            LP%W = SF%VEL_T(1)
         CASE( 3)
            LP%U = SF%VEL_T(1)
            LP%V = SF%VEL_T(2)
            LP%W = -WALL(IW)%ONE_D%UW
         CASE(-3)
            LP%U = SF%VEL_T(1)
            LP%V = SF%VEL_T(2)
            LP%W =  WALL(IW)%ONE_D%UW
      END SELECT
   
      LP%ONE_D%IIG = IIG
      LP%ONE_D%JJG = JJG
      LP%ONE_D%KKG = KKG
         
      ! Save the insertion time (TP) and scalar property (SP) for the particle
   
      IF (MOD(NLP,LPC%SAMPLING)==0) LP%SHOW = .TRUE.
      LP%T_INSERT = T

      CALL MAKE_PARTICLE

      LP=>MESHES(NM)%LAGRANGIAN_PARTICLE(NLP)      
      SF=>SURFACE(LPC%SURF_INDEX)
      IF (.NOT. LPC%MASSLESS_TRACER .AND. .NOT.LPC%MASSLESS_TARGET) THEN
         MASS_SUM = MASS_SUM + LP%PWT*LPC%FTPR*LP%ONE_D%X(SF%N_CELLS_INI)**3
      ENDIF

   ENDDO PARTICLE_INSERT_LOOP2

   ! Adjust the particle weighting factors to get the right mass flux

   IF (MASS_SUM > 0._EB) THEN
      SF => SURFACE(WC%SURF_INDEX)
      IF (SF%PARTICLE_MASS_FLUX > 0._EB) THEN
         IF (ABS(WC%ONE_D%T_IGN-T_BEGIN)<=TWO_EPSILON_EB .AND. SF%RAMP_INDEX(TIME_PART)>=1) THEN
            TSI = T
         ELSE
            TSI = T - WC%ONE_D%T_IGN
         ENDIF
         FLOW_RATE = EVALUATE_RAMP(TSI,SF%TAU(TIME_PART),SF%RAMP_INDEX(TIME_PART))*SF%PARTICLE_MASS_FLUX         
         DO I=1,SF%NPPC
            N = LP_INDEX_LOOKUP(I)
            LAGRANGIAN_PARTICLE(N)%PWT = LAGRANGIAN_PARTICLE(N)%PWT*&
                                         FLOW_RATE*WALL(IW)%ONE_D%AREA_ADJUST*WALL(IW)%AW*SF%DT_INSERT/MASS_SUM
         ENDDO
      ENDIF
   ENDIF
   
   DEALLOCATE(LP_INDEX_LOOKUP)
   
ENDDO WALL_INSERT_LOOP

END SUBROUTINE WALL_PARTICLE_INSERT


SUBROUTINE VOLUME_PARTICLE_INSERT

! Loop over all INIT lines and look for particles inserted within a specified volume

INTEGER :: I,ND,N_INSERT,I1,J1,K1,I2,J2,K2,MOD_N_INSERT,N,N_PARTICLES_INSERT
REAL(EB) :: XC1,XC2,YC1,YC2,ZC1,ZC2,X0,Y0,Z0,RR,HH,INSERT_VOLUME,INPUT_VOLUME

VOLUME_INSERT_LOOP: DO IB=1,N_INIT

   IN => INITIALIZATION(IB)

   ! Determine if the INITIALIZATION type involves particles. If not, cycle.

   ILPC = IN%PART_INDEX
   IF (ILPC<1) CYCLE VOLUME_INSERT_LOOP
   IF (IN%SINGLE_INSERTION .AND. IN%ALREADY_INSERTED(NM)) CYCLE VOLUME_INSERT_LOOP

   ! Determine if the particles/PARTICLEs are controlled by devices

   LPC => LAGRANGIAN_PARTICLE_CLASS(ILPC)

   IF (IN%DEVC_INDEX>0) THEN
      IF (.NOT.DEVICE(IN%DEVC_INDEX)%CURRENT_STATE) THEN
         IN%PARTICLE_INSERT_CLOCK(NM) = T
         CYCLE VOLUME_INSERT_LOOP
      ENDIF
   ENDIF
   IF (IN%CTRL_INDEX>0) THEN
      IF (.NOT.CONTROL(IN%CTRL_INDEX)%CURRENT_STATE) THEN
         IN%PARTICLE_INSERT_CLOCK(NM) = T
         CYCLE VOLUME_INSERT_LOOP
      ENDIF
   ENDIF

   ! If it is not time to insert particles/PARTICLEs for this INITIALIZATION block, cycle.
   
   IF (T < IN%PARTICLE_INSERT_CLOCK(NM)) CYCLE VOLUME_INSERT_LOOP

   ! Start processing the INITIALIZATION info

   IF (IN%N_PARTICLES==0 .AND. IN%N_PARTICLES_PER_CELL==0) CYCLE VOLUME_INSERT_LOOP

   ! Cut off parts of the INIT region that are outside the current mesh

   IF (IN%X1>XF .OR. IN%X2<XS .OR. IN%Y1>YF .OR. IN%Y2<YS .OR. IN%Z1>ZF .OR. IN%Z2<ZS) CYCLE VOLUME_INSERT_LOOP
   X1 = MAX(IN%X1,XS) 
   X2 = MIN(IN%X2,XF)
   Y1 = MAX(IN%Y1,YS) 
   Y2 = MIN(IN%Y2,YF)
   Z1 = MAX(IN%Z1,ZS) 
   Z2 = MIN(IN%Z2,ZF)
   
   ! Compute the volume of the INIT region

   SELECT CASE(IN%SHAPE)      
      CASE('BLOCK')
         INSERT_VOLUME = (X2-X1)*(Y2-Y1)*(Z2-Z1)
         INPUT_VOLUME  = (IN%X2-IN%X1)*(IN%Y2-IN%Y1)*(IN%Z2-IN%Z1)
      CASE('CONE')
         INSERT_VOLUME = PI*(0.5_EB*(IN%X2-IN%X1))**2*(IN%Z2-IN%Z1)/3._EB
         INPUT_VOLUME  = INSERT_VOLUME
      CASE('RING')
         INSERT_VOLUME = 0._EB
         INPUT_VOLUME  = 0._EB
   END SELECT

   IF (INSERT_VOLUME<=0._EB .AND. IN%MASS_PER_VOLUME>0._EB) CYCLE VOLUME_INSERT_LOOP

   ! Assign properties to the initial PARTICLEs/particles

   MASS_SUM = 0._EB

   TOTAL_OR_PER_CELL: IF (IN%N_PARTICLES > 0) THEN

      ! If the original INIT volume expands over multiple meshes, insert only a fraction of the specified N_PARTICLES.

      IF (INPUT_VOLUME>TWO_EPSILON_EB) THEN
         N_PARTICLES_INSERT = MAX(1,NINT(IN%N_PARTICLES*INSERT_VOLUME/INPUT_VOLUME))
      ELSE
         N_PARTICLES_INSERT = IN%N_PARTICLES
      ENDIF

      ! Array for speeding up sort-by-age routine
      
      ALLOCATE(LP_INDEX_LOOKUP(N_PARTICLES_INSERT))
      LP_INDEX_LOOKUP = 0
      
      INSERT_PARTICLE_LOOP: DO I=1,N_PARTICLES_INSERT

         IF (NLP+1>MAXIMUM_PARTICLES) THEN
            CALL REMOVE_OLDEST_PARTICLE(NM,ILPC,NLP,NEW_LP_INDEX)
            IF (I>1) LP_INDEX_LOOKUP(I-1) = NEW_LP_INDEX
         ELSE
            NLP = NLP+1
         ENDIF
         LP_INDEX_LOOKUP(I) = NLP
         
         PARTICLE_TAG = PARTICLE_TAG + NMESHES
         CALL ALLOCATE_STORAGE(NM,LAGRANGIAN_PARTICLE_CLASS(ILPC)%SURF_INDEX,LPC_INDEX=ILPC,LP_INDEX=NLP,&
                               TAG=PARTICLE_TAG,NEW_TAG=.TRUE.)
         LAGRANGIAN_PARTICLE => MESHES(NM)%LAGRANGIAN_PARTICLE
         LP=>MESHES(NM)%LAGRANGIAN_PARTICLE(NLP)

         ! Get particle coordinates by randomly choosing within the designated volume

         N = 0
         CHOOSE_XYZ_LOOP:  DO
            N = N + 1
            SELECT CASE(IN%SHAPE)
               CASE('BLOCK') 
                  CALL RANDOM_RECTANGLE(LP%X,LP%Y,LP%Z,X1,X2,Y1,Y2,Z1,Z2)
               CASE('CONE')  
                  X0 = 0.5_EB*(IN%X1+IN%X2)
                  Y0 = 0.5_EB*(IN%Y1+IN%Y2)
                  Z0 = IN%Z1
                  RR = 0.5_EB*(IN%X2-IN%X1)
                  HH = IN%Z2-IN%Z1
                  CALL RANDOM_CONE(NM,LP%X,LP%Y,LP%Z,X0,Y0,Z0,RR,HH)
               CASE('RING')
                  X0 = 0.5_EB*(IN%X1+IN%X2)
                  Y0 = 0.5_EB*(IN%Y1+IN%Y2)
                  RR = 0.5_EB*(IN%X2-IN%X1)
                  LP%Z = IN%Z1
                  CALL RANDOM_RING(NM,LP%X,LP%Y,X0,Y0,RR)
            END SELECT
            CALL GET_IJK(LP%X,LP%Y,LP%Z,NM,XI,YJ,ZK,II,JJ,KK)
            LP%X = LP%X + (I-1)*IN%DX
            LP%Y = LP%Y + (I-1)*IN%DY
            LP%Z = LP%Z + (I-1)*IN%DZ
            IF (.NOT.SOLID(CELL_INDEX(II,JJ,KK))) EXIT CHOOSE_XYZ_LOOP
            ! If cannot find non-solid grid cell, stop searching
            IF (N>10000) EXIT CHOOSE_XYZ_LOOP
         ENDDO CHOOSE_XYZ_LOOP

         CALL VOLUME_INIT_PARTICLE(I)
                  
         IN => INITIALIZATION(IB)
         LP => LAGRANGIAN_PARTICLE(NLP)           

      ENDDO INSERT_PARTICLE_LOOP
      
      N_INSERT = N_PARTICLES_INSERT

   ELSEIF (IN%N_PARTICLES_PER_CELL > 0) THEN TOTAL_OR_PER_CELL

      N_INSERT = 0
      INSERT_VOLUME = 0._EB
      CALL GET_IJK(X1,Y1,Z1,NM,XI,YJ,ZK,I1,J1,K1)
      CALL GET_IJK(X2,Y2,Z2,NM,XI,YJ,ZK,I2,J2,K2)
      I2 = MIN(I2,IBAR)
      J2 = MIN(J2,JBAR)
      K2 = MIN(K2,KBAR)
      N_INSERT = MIN(MAXIMUM_PARTICLES,(I2-I1+1)*(J2-J1+1)*(K2-K1+1)*IN%N_PARTICLES_PER_CELL)
      ALLOCATE(LP_INDEX_LOOKUP(N_INSERT))
      LP_INDEX_LOOKUP = 0
      N_INSERT = 0
      
      DO KK=K1,K2
         DO JJ=J1,J2
            II_LOOP: DO II=I1,I2
               IF (SOLID(CELL_INDEX(II,JJ,KK))) CYCLE II_LOOP
               INSERT_VOLUME = INSERT_VOLUME + DX(II)*DY(JJ)*DZ(KK)               
               INSERT_PARTICLE_LOOP_2: DO I = 1, IN%N_PARTICLES_PER_CELL
                  N_INSERT = N_INSERT + 1
                  IF (N_INSERT > MAXIMUM_PARTICLES) THEN
                     N_INSERT = N_INSERT - 1
                     EXIT INSERT_PARTICLE_LOOP_2
                  ENDIF

                  IF (NLP+1>MAXIMUM_PARTICLES) THEN
                     CALL REMOVE_OLDEST_PARTICLE(NM,ILPC,NLP,NEW_LP_INDEX)
                     IF (N_INSERT>1) LP_INDEX_LOOKUP(N_INSERT-1) = NEW_LP_INDEX
                  ELSE
                     NLP = NLP+1
                  ENDIF
                  LP_INDEX_LOOKUP(N_INSERT) = NLP
                  
                  PARTICLE_TAG = PARTICLE_TAG + NMESHES
                  CALL ALLOCATE_STORAGE(NM,LAGRANGIAN_PARTICLE_CLASS(ILPC)%SURF_INDEX,LPC_INDEX=ILPC,LP_INDEX=NLP,&
                                        TAG=PARTICLE_TAG,NEW_TAG=.TRUE.)
                  LAGRANGIAN_PARTICLE => MESHES(NM)%LAGRANGIAN_PARTICLE
                  LP=>MESHES(NM)%LAGRANGIAN_PARTICLE(NLP)

                  ! Get particle coordinates by randomly choosing within the designated volume

                  XC1 = MAX(X1,X(II-1))
                  YC1 = MAX(Y1,Y(JJ-1))
                  ZC1 = MAX(Z1,Z(KK-1))
                  XC2 = MIN(X2,X(II))
                  YC2 = MIN(Y2,Y(JJ))
                  ZC2 = MIN(Z2,Z(KK))

                  IF (IN%CELL_CENTERED) THEN
                     LP%X = 0.5_EB*(XC1+XC2)
                     LP%Y = 0.5_EB*(YC1+YC2)
                     LP%Z = 0.5_EB*(ZC1+ZC2)
                  ELSE
                     CALL RANDOM_RECTANGLE(LP%X,LP%Y,LP%Z,XC1,XC2,YC1,YC2,ZC1,ZC2)                     
                  ENDIF
                  
                  CALL VOLUME_INIT_PARTICLE(I)
                  
                  IN => INITIALIZATION(IB)
                  LP => LAGRANGIAN_PARTICLE(NLP)
       
               ENDDO INSERT_PARTICLE_LOOP_2
            ENDDO II_LOOP
         ENDDO
      ENDDO

   ENDIF TOTAL_OR_PER_CELL

   ! Adjust particle weighting factor PWT so that desired MASS_PER_VOLUME is achieved

   IF (IN%MASS_PER_TIME>0._EB) THEN
      PWT0 = IN%MASS_PER_TIME*IN%DT_INSERT/MASS_SUM
   ELSEIF (IN%MASS_PER_VOLUME>0._EB) THEN
      PWT0 = IN%MASS_PER_VOLUME*INSERT_VOLUME/MASS_SUM
   ELSE
      PWT0 = IN%PARTICLE_WEIGHT_FACTOR
   ENDIF
   
   DO I=1,MIN(MAXIMUM_PARTICLES,N_INSERT)
      IP = LP_INDEX_LOOKUP(I)
      LAGRANGIAN_PARTICLE(IP)%PWT = LAGRANGIAN_PARTICLE(IP)%PWT*PWT0   
   ENDDO
   
   DEALLOCATE(LP_INDEX_LOOKUP)
   
   IN%ALREADY_INSERTED(NM) = .TRUE.                                             

ENDDO VOLUME_INSERT_LOOP

END SUBROUTINE VOLUME_PARTICLE_INSERT


SUBROUTINE VOLUME_INIT_PARTICLE(I)

! Initialize particle indices and velocity

INTEGER, INTENT(IN) :: I
INTEGER :: ND
IN => INITIALIZATION(IB)
LP => LAGRANGIAN_PARTICLE(NLP)

LP%ONE_D%IIG = II
LP%ONE_D%JJG = JJ
LP%ONE_D%KKG = KK
LP%U = IN%U0
LP%V = IN%V0
LP%W = IN%W0

! If the INITIALIZATION group has an ID, match it with a device

IF (IN%ID/='null') THEN
   DO ND=1,N_DEVC
      DV => DEVICE(ND)
      IF (IN%ID==DV%INIT_ID .AND. I==DV%POINT) THEN
         DV%LP_TAG = PARTICLE_TAG
         DV%PART_INDEX = ILPC
         DV%MESH = NM
         DV%X = LP%X
         DV%Y = LP%Y
         DV%Z = LP%Z
         IF (DV%LINE>0 .AND. DV%LINE_COORD_CODE==123) THEN
            IF (ABS(IN%DX)>TWO_EPSILON_EB .AND. ABS(IN%DY)<TWO_EPSILON_EB .AND. ABS(IN%DZ)<TWO_EPSILON_EB) DV%LINE_COORD_CODE = 1
            IF (ABS(IN%DX)<TWO_EPSILON_EB .AND. ABS(IN%DY)>TWO_EPSILON_EB .AND. ABS(IN%DZ)<TWO_EPSILON_EB) DV%LINE_COORD_CODE = 2
            IF (ABS(IN%DX)<TWO_EPSILON_EB .AND. ABS(IN%DY)<TWO_EPSILON_EB .AND. ABS(IN%DZ)>TWO_EPSILON_EB) DV%LINE_COORD_CODE = 3
         ENDIF
      ENDIF
   ENDDO
ENDIF

! Process particle and set more initial values

CALL MAKE_PARTICLE

LP=>LAGRANGIAN_PARTICLE(NLP)

! Special over-ride of particle diameter for liquid droplets only

IF (IN%DIAMETER>0._EB) THEN
   LP%ONE_D%X(1) = 0.5_EB*IN%DIAMETER
   LP%PWT = 1._EB
   LP%ONE_D%LAYER_THICKNESS(1) = LP%ONE_D%X(1)
   LP%MASS = FOTHPI*LP%ONE_D%RHO(1,1)*LP%ONE_D%X(1)**3
ENDIF

! Save insert time and other miscellaneous attributes

LP%T_INSERT = T                       
IF (MOD(NLP,LPC%SAMPLING)==0) LP%SHOW = .TRUE.    
   
! Get the particle ORIENTATION from the PART line

IF (LPC%N_ORIENTATION>0) THEN
   LP%ORIENTATION_INDEX = LPC%ORIENTATION_INDEX + MOD(I-1,LPC%N_ORIENTATION)
   LP%PWT = LP%PWT*LPC%SOLID_ANGLE(MOD(I-1,LPC%N_ORIENTATION)+1)/(4._EB*PI)
ENDIF     

! Get the particle ORIENTATION from the DEVC line

IF (IN%DEVC_INDEX>0) THEN
   IF (DEVICE(IN%DEVC_INDEX)%QUANTITY=='RADIATIVE HEAT FLUX') LP%ORIENTATION_INDEX = DEVICE(IN%DEVC_INDEX)%ORIENTATION_INDEX
ENDIF

! Sum up the particle masses

MASS_SUM = MASS_SUM + LP%PWT*LP%MASS ! if r=0 the sum will stay 0

END SUBROUTINE VOLUME_INIT_PARTICLE


SUBROUTINE MAKE_PARTICLE

REAL(EB) :: X1,X2,AREA,LENGTH,SCALE_FACTOR,RADIUS,MPUA,LP_VOLUME
INTEGER :: N
TYPE (ONE_D_M_AND_E_XFER_TYPE), POINTER :: ONE_D=>NULL()

SF => SURFACE(LPC%SURF_INDEX)
LP => LAGRANGIAN_PARTICLE(NLP)
ONE_D => LP%ONE_D

IF (LPC%SOLID_PARTICLE) THEN

   IF (LPC%SURF_INDEX==TGA_SURF_INDEX) TGA_PARTICLE_INDEX = NLP

   LP%MASS = 0._EB

   SELECT CASE (LPC%DRAG_LAW)
   
      CASE(SCREEN_DRAG)

         ! Compute special cross-sectional area of screen particle

         AREA = (ABS(ORIENTATION_VECTOR(1,LPC%ORIENTATION_INDEX))*DY(ONE_D%JJG)*DZ(ONE_D%KKG) + &
                 ABS(ORIENTATION_VECTOR(2,LPC%ORIENTATION_INDEX))*DX(ONE_D%IIG)*DZ(ONE_D%KKG) + &
                 ABS(ORIENTATION_VECTOR(3,LPC%ORIENTATION_INDEX))*DX(ONE_D%IIG)*DY(ONE_D%JJG)) * &
                 LAGRANGIAN_PARTICLE_CLASS(LP%CLASS_INDEX)%FREE_AREA_FRACTION
         SELECT CASE (SF%GEOMETRY)
            CASE (SURF_CARTESIAN)
               LP%ONE_D%AREA = AREA
               IF (SF%THERMALLY_THICK) THEN
                  DO N=1,SF%N_LAYERS
                     LP%MASS = LP%MASS + AREA*SF%LAYER_THICKNESS(N)*SF%LAYER_DENSITY(N)
                  END DO
               ENDIF 
            CASE (SURF_CYLINDRICAL)
               LP%ONE_D%AREA = AREA*PI            
               IF (SF%THERMALLY_THICK) THEN
                  X1 = SUM(SF%LAYER_THICKNESS)
                  LENGTH = AREA / (2._EB*X1)
                  DO N=SF%N_LAYERS,1,-1
                     X2 = X1 - SF%LAYER_THICKNESS(N)
                     LP%MASS = LP%MASS + LENGTH*SF%LAYER_DENSITY(N)*PI*(X1**2-X2**2)
                     X1 = X2
                  END DO      
               ENDIF
            CASE (SURF_SPHERICAL)
               LP%ONE_D%AREA = AREA*4._EB
               IF (SF%THERMALLY_THICK) THEN
                  X1 = SUM(SF%LAYER_THICKNESS)               
                  LP%PWT = AREA/(PI*X1**2)
                  DO N=SF%N_LAYERS,1,-1
                     X2 = X1 - SF%LAYER_THICKNESS(N)
                     LP%MASS = LP%MASS + SF%LAYER_DENSITY(N)*FOTHPI*(X1**3-X2**3)
                     X1 = X2
                  END DO                     
               ELSE
                  LP%PWT = AREA/(PI*SF%THICKNESS**2)               
               ENDIF            
         END SELECT
      CASE (POROUS_DRAG)
         MPUA = 0._EB
         LP_VOLUME = LPC%POROUS_VOLUME_FRACTION*DX(II)*DY(JJ)*DZ(KK)
         SELECT CASE (SF%GEOMETRY)
            CASE (SURF_CARTESIAN)
               LP%ONE_D%AREA = LP_VOLUME/SF%THICKNESS            
               IF (SF%THERMALLY_THICK) THEN
                  DO N=1,SF%N_LAYERS
                     LP%MASS = LP%MASS + AREA*SF%LAYER_THICKNESS(N)*SF%LAYER_DENSITY(N)
                  END DO
                  LP%MASS = LP%MASS*LP%ONE_D%AREA
               ENDIF
            CASE (SURF_CYLINDRICAL)
               LP%ONE_D%AREA = 2._EB*LP_VOLUME/SF%THICKNESS
               IF (SF%THERMALLY_THICK) THEN
                  X1 = SUM(SF%LAYER_THICKNESS)
                  LENGTH = LP_VOLUME/(PI*X1**2)
                  DO N=SF%N_LAYERS,1,-1
                     X2 = X1 - SF%LAYER_THICKNESS(N)
                     LP%MASS = LP%MASS + LENGTH*SF%LAYER_DENSITY(N)*PI*(X1**2-X2**2)
                     X1 = X2
                  END DO    
               ENDIF
            CASE (SURF_SPHERICAL)
               LP%ONE_D%AREA = 3._EB*LP_VOLUME/SF%THICKNESS
               LP%PWT = AREA/(4._EB*PI*SF%THICKNESS**2)               
               IF (SF%THERMALLY_THICK) THEN
                  X1 = SUM(SF%LAYER_THICKNESS)               
                  DO N=SF%N_LAYERS,1,-1
                     X2 = X1 - SF%LAYER_THICKNESS(N)
                     LP%MASS = LP%MASS + SF%LAYER_DENSITY(N)*FOTHPI*(X1**3-X2**3)
                     X1 = X2
                  END DO                     
               ENDIF            
         END SELECT         

      CASE DEFAULT
            
         IF (.NOT.LPC%MONODISPERSE) THEN
            CALL PARTICLE_SIZE_WEIGHT(RADIUS,LP%PWT)
            SCALE_FACTOR = RADIUS/SF%THICKNESS
            ONE_D%X(:) = ONE_D%X(:)*SCALE_FACTOR
            ONE_D%LAYER_THICKNESS(:) = ONE_D%LAYER_THICKNESS(:)*SCALE_FACTOR
         ELSE
            SCALE_FACTOR = 1._EB
         ENDIF

         IF (SF%THERMALLY_THICK) THEN
            SELECT CASE (SF%GEOMETRY)
               CASE (SURF_CARTESIAN)
                  DO N=1,SF%N_LAYERS
                     LP%MASS = LP%MASS + 2._EB*SF%LENGTH*SF%WIDTH*SF%LAYER_THICKNESS(N)*SCALE_FACTOR*SF%LAYER_DENSITY(N)
                  ENDDO
               CASE (SURF_CYLINDRICAL)
                  X1 = SUM(SF%LAYER_THICKNESS)*SCALE_FACTOR
                  DO N=SF%N_LAYERS,1,-1
                     X2 = X1 - SF%LAYER_THICKNESS(N)*SCALE_FACTOR
                     LP%MASS = LP%MASS + SF%LENGTH*SF%LAYER_DENSITY(N)*PI*(X1**2-X2**2)
                     X1 = X2
                  ENDDO      
               CASE (SURF_SPHERICAL)
                  X1 = SUM(SF%LAYER_THICKNESS)*SCALE_FACTOR
                  DO N=SF%N_LAYERS,1,-1
                     X2 = X1 - SF%LAYER_THICKNESS(N)*SCALE_FACTOR
                     LP%MASS = LP%MASS + SF%LAYER_DENSITY(N)*FOTHPI*(X1**3-X2**3)
                     X1 = X2
                  ENDDO      
            END SELECT
         ENDIF

   END SELECT   

ELSEIF (LPC%LIQUID_DROPLET) THEN

   LP%ONE_D%RHO(1,1) = LPC%DENSITY
   ONE_D%N_LAYER_CELLS(1:SF%N_LAYERS) = 1
   CALL PARTICLE_SIZE_WEIGHT(ONE_D%X(1),LP%PWT)
   ONE_D%LAYER_THICKNESS(1) = ONE_D%X(1)
   LP%MASS = FOTHPI*ONE_D%RHO(1,1)*ONE_D%X(1)**3

ENDIF

ONE_D%TMP(0:SF%N_CELLS_INI+1) = LPC%TMP_INITIAL
ONE_D%TMP_F = ONE_D%TMP(1)   

! Check if fire spreads radially over this surface type, and if so, set T_IGN appropriately

IF (SF%XYZ(1)>-1.E5_EB) THEN
   ONE_D%T_IGN = T_BEGIN + SQRT((LP%X-SF%XYZ(1))**2 +(LP%Y-SF%XYZ(2))**2 +(LP%Z-SF%XYZ(3))**2)/SF%FIRE_SPREAD_RATE
ELSE
   ONE_D%T_IGN = SF%T_IGN
ENDIF

END SUBROUTINE MAKE_PARTICLE


SUBROUTINE PARTICLE_SIZE_WEIGHT(R,PWT)
REAL(EB), INTENT(OUT):: R,PWT

IF (LPC%MONODISPERSE) THEN
   R   = 0.5_EB*LPC%DIAMETER
   PWT = 1._EB
ELSE
   CALL RANDOM_NUMBER(RN)            
   STRATUM = NINT(REAL(LPC%N_STRATA,EB)*REAL(RN,EB)+0.5_EB)
   IL = LPC%IL_CNF(STRATUM)
   IU = LPC%IU_CNF(STRATUM)
   CALL RANDOM_CHOICE(LPC%CNF(IL:IU),LPC%R_CNF(IL:IU),IU-IL,R)
   PWT = LPC%W_CNF(STRATUM)
   IF (2._EB*R > LPC%MAXIMUM_DIAMETER) THEN
      PWT = PWT*R**3/(0.5_EB*LPC%MAXIMUM_DIAMETER)**3
      R = 0.5_EB*LPC%MAXIMUM_DIAMETER
   ENDIF
   IF (2._EB*R < LPC%MINIMUM_DIAMETER) THEN
      PWT = PWT*R**3/(0.5_EB*LPC%MINIMUM_DIAMETER)**3
      R = 0.5_EB*LPC%MINIMUM_DIAMETER
   ENDIF   
ENDIF

END SUBROUTINE PARTICLE_SIZE_WEIGHT
 

END SUBROUTINE INSERT_PARTICLES
 
 

SUBROUTINE UPDATE_PARTICLES(T,NM)

REAL(EB), INTENT(IN) :: T
INTEGER, INTENT(IN) :: NM
INTEGER :: NOM
REAL(EB) :: TNOW

IF (EVACUATION_ONLY(NM)) RETURN  

! Zero out the number of the PARTICLEs in the "orphanage"; that is, the place to hold PARTICLEs transferring from mesh to mesh

DO NOM=1,NMESHES
   IF (ALLOCATED(MESHES(NM)%OMESH(NOM)%N_PART_ORPHANS)) MESHES(NM)%OMESH(NOM)%N_PART_ORPHANS = 0
ENDDO

! Return if there are no particles in this mesh

IF (MESHES(NM)%NLP==0) THEN
   IF (CORRECTOR .AND. CALC_D_LAGRANGIAN) THEN
      D_LAGRANGIAN = 0._EB
      CALC_D_LAGRANGIAN=.FALSE.
   ENDIF
   RETURN
ENDIF

! Set the CPU timer and point to the current mesh variables

TNOW=SECOND()
CALL POINT_TO_MESH(NM)

! Zero out the contribution by lagrangian particles to divergence

IF (N_LP_ARRAY_INDICES>0 .AND. CORRECTOR) THEN
   D_LAGRANGIAN = 0._EB
   CALC_D_LAGRANGIAN=.FALSE.
ENDIF

! Move the PARTICLEs/particles, then compute mass and energy transfer, then add PARTICLE momentum to gas

IF (CORRECTOR) CALL MOVE_PARTICLES(T,NM)
IF (CORRECTOR) CALL PARTICLE_MASS_ENERGY_TRANSFER(T,NM)

TUSED(8,NM)=TUSED(8,NM)+SECOND()-TNOW

END SUBROUTINE UPDATE_PARTICLES


SUBROUTINE MOVE_PARTICLES(T,NM)

! Momentum transfer from all particles and PARTICLEs
 
USE TRAN, ONLY: GET_IJK 
REAL(EB), INTENT(IN) :: T
REAL     :: RN
REAL(EB) :: SURFACE_PARTICLE_DIAMETER,XI,YJ,ZK,R_D,R_D_0,X_OLD,Y_OLD,Z_OLD,THETA_RN,STEP_FRACTION(-3:3),RDT,THROHALF,B_1
LOGICAL :: HIT_SOLID
INTEGER :: IP,ICN,IIN,JJN,KKN,IW,IWP1,IWM1,IWP2,IWM2,IWP3,IWM3,IOR_OLD,IC,IOR_FIRST,IML,IIG,JJG,KKG
INTEGER, INTENT(IN) :: NM
TYPE (LAGRANGIAN_PARTICLE_TYPE), POINTER :: LP=>NULL()
TYPE (LAGRANGIAN_PARTICLE_CLASS_TYPE), POINTER :: LPC=>NULL()
TYPE (SURFACE_TYPE), POINTER :: SF
REAL(EB), POINTER, DIMENSION(:,:,:) :: NDPC=>NULL()

CALL POINT_TO_MESH(NM)

THROHALF = (0.5_EB)**(1./3.)
B_1 =  1.7321_EB ! SQRT(3)
RDT = 1._EB/DT

SURFACE_PARTICLE_DIAMETER = 0.001_EB  ! All PARTICLEs adjusted to this size when on solid (m)

! Sum up the number of PARTICLEs/particles in each grid cell (NDPC -- Number PARTICLEs Per Cell)

NDPC=>WORK1
NDPC=0._EB

DO IP=1,NLP
   LP => LAGRANGIAN_PARTICLE(IP)
   CALL GET_IJK(LP%X,LP%Y,LP%Z,NM,XI,YJ,ZK,LP%ONE_D%IIG,LP%ONE_D%JJG,LP%ONE_D%KKG)  
   IF (LP%PWT>0._EB .AND. LP%ONE_D%IOR==0) NDPC(LP%ONE_D%IIG,LP%ONE_D%JJG,LP%ONE_D%KKG) = &
                                           NDPC(LP%ONE_D%IIG,LP%ONE_D%JJG,LP%ONE_D%KKG)+LP%PWT
ENDDO

PART_UVWMAX = 0._EB

! Loop through all Lagrangian particles and move them one time step

PARTICLE_LOOP: DO IP=1,NLP  

   ! Assign particle (LP%) and particle class (LPC%) shortcuts

   LP  => LAGRANGIAN_PARTICLE(IP)
   LPC => LAGRANGIAN_PARTICLE_CLASS(LP%CLASS_INDEX)
   SF  => SURFACE(LPC%SURF_INDEX)
   
   ! Determine particle radius

   IF (LPC%LIQUID_DROPLET) THEN
      R_D = LP%ONE_D%X(1)
   ELSEIF (LPC%SOLID_PARTICLE) THEN
      R_D = SUM(LP%ONE_D%LAYER_THICKNESS(1:SF%N_LAYERS))
   ELSEIF (LPC%MASSLESS_TRACER) THEN
      R_D = 0._EB
   ELSEIF (LPC%MASSLESS_TARGET) THEN
      CYCLE PARTICLE_LOOP
   ENDIF

   IF (.NOT. LPC%MASSLESS_TRACER .AND. (R_D<=0._EB .OR. (.NOT. LPC%STATIC .AND. LP%MASS<=TWO_EPSILON_EB))) CYCLE PARTICLE_LOOP

   R_D_0 = R_D

   ! Determine the current coordinates of the particle

   IIG = LP%ONE_D%IIG
   JJG = LP%ONE_D%JJG
   KKG = LP%ONE_D%KKG

   IC = CELL_INDEX(IIG,JJG,KKG)

   X_OLD = LP%X
   Y_OLD = LP%Y
   Z_OLD = LP%Z

   ! Throw out particles that are inside a solid obstruction

   IF (SOLID(IC)) THEN
      LP%X = 1.E6_EB
      CYCLE PARTICLE_LOOP
   ENDIF

   SOLID_GAS_MOVE: IF (LP%ONE_D%IOR/=0) THEN

      CALL MOVE_ON_SOLID

   ELSE SOLID_GAS_MOVE

      CALL MOVE_IN_GAS

      ! If the particle does not move, but does drag, go on to the next particle

      IF (LPC%MASSLESS_TRACER .OR. LP%PWT<=TWO_EPSILON_EB .OR. LPC%STATIC) CYCLE PARTICLE_LOOP   

   ENDIF SOLID_GAS_MOVE

   ! Special case where a particle hits a POROUS_FLOOR
    
   IF (POROUS_FLOOR .AND. LP%Z<ZS .AND. LPC%LIQUID_DROPLET) THEN
      IC = CELL_INDEX(IIG,JJG,1)
      IW = WALL_INDEX(IC,-3)
      IF (WALL(IW)%BOUNDARY_TYPE==SOLID_BOUNDARY .AND. ACCUMULATE_WATER .AND. .NOT.LP%SPLAT) THEN
         WALL(IW)%A_LP_MPUA(LPC%ARRAY_INDEX) = WALL(IW)%A_LP_MPUA(LPC%ARRAY_INDEX) + LP%PWT*LPC%FTPR*R_D**3*WALL(IW)%RAW
         LP%SPLAT = .TRUE.
      ENDIF
      CYCLE PARTICLE_LOOP
   ENDIF

   ! Where is the PARTICLE now? Limit the location by UBOUND and LBOUND due to the possible super fast PARTICLEs

   IIN = MAX(LBOUND(CELLSI,1),MIN(UBOUND(CELLSI,1),FLOOR((LP%X-XS)*RDXINT)))
   JJN = MAX(LBOUND(CELLSJ,1),MIN(UBOUND(CELLSJ,1),FLOOR((LP%Y-YS)*RDYINT)))
   KKN = MAX(LBOUND(CELLSK,1),MIN(UBOUND(CELLSK,1),FLOOR((LP%Z-ZS)*RDZINT)))
   XI  = CELLSI(IIN)
   YJ  = CELLSJ(JJN)
   ZK  = CELLSK(KKN)
   IIN = FLOOR(XI+1._EB)
   JJN = FLOOR(YJ+1._EB)
   KKN = FLOOR(ZK+1._EB)
   IF (IIN<0 .OR. IIN>IBP1) CYCLE PARTICLE_LOOP
   IF (JJN<0 .OR. JJN>JBP1) CYCLE PARTICLE_LOOP
   IF (KKN<0 .OR. KKN>KBP1) CYCLE PARTICLE_LOOP
   ICN = CELL_INDEX(IIN,JJN,KKN)
   IF (IC==0 .OR. ICN==0) CYCLE PARTICLE_LOOP

   IF (LP%X<XS .AND. WALL(WALL_INDEX(IC,-1))%BOUNDARY_TYPE/=SOLID_BOUNDARY) CYCLE PARTICLE_LOOP
   IF (LP%X>XF .AND. WALL(WALL_INDEX(IC, 1))%BOUNDARY_TYPE/=SOLID_BOUNDARY) CYCLE PARTICLE_LOOP
   IF (LP%Y<YS .AND. WALL(WALL_INDEX(IC,-2))%BOUNDARY_TYPE/=SOLID_BOUNDARY) CYCLE PARTICLE_LOOP
   IF (LP%Y>YF .AND. WALL(WALL_INDEX(IC, 2))%BOUNDARY_TYPE/=SOLID_BOUNDARY) CYCLE PARTICLE_LOOP
   IF (LP%Z<ZS .AND. WALL(WALL_INDEX(IC,-3))%BOUNDARY_TYPE/=SOLID_BOUNDARY) CYCLE PARTICLE_LOOP
   IF (LP%Z>ZF .AND. WALL(WALL_INDEX(IC, 3))%BOUNDARY_TYPE/=SOLID_BOUNDARY) CYCLE PARTICLE_LOOP

   ! If PARTICLE hits an obstacle, change its properties

   AIR_TO_SOLID: IF (IIG/=IIN .OR. JJG/=JJN .OR. KKG/=KKN) THEN

      IOR_OLD   = LP%ONE_D%IOR
      HIT_SOLID = .FALSE.

      ! Check if any solid boundaries of original grid cell have been crossed

      IWP1 = WALL_INDEX(IC, 1) 
      IWM1 = WALL_INDEX(IC,-1)
      IWP2 = WALL_INDEX(IC, 2)
      IWM2 = WALL_INDEX(IC,-2)
      IWP3 = WALL_INDEX(IC, 3)
      IWM3 = WALL_INDEX(IC,-3)
      STEP_FRACTION = 1._EB

      IF (KKN>KKG .AND. WALL(IWP3)%BOUNDARY_TYPE==SOLID_BOUNDARY) THEN
         LP%ONE_D%IOR=-3
         HIT_SOLID = .TRUE.
         STEP_FRACTION(LP%ONE_D%IOR) = MAX(0._EB,(Z(KKG)-Z_OLD-0.05_EB*DZ(KKG))/(LP%Z-Z_OLD))
      ENDIF
      IF (KKN<KKG .AND. WALL(IWM3)%BOUNDARY_TYPE==SOLID_BOUNDARY) THEN
         LP%ONE_D%IOR= 3
         HIT_SOLID = .TRUE.
         STEP_FRACTION(LP%ONE_D%IOR) = MAX(0._EB,(Z(KKG-1)-Z_OLD+0.05_EB*DZ(KKG-1))/(LP%Z-Z_OLD))
      ENDIF
      IF (IIN>IIG .AND. WALL(IWP1)%BOUNDARY_TYPE==SOLID_BOUNDARY) THEN
         LP%ONE_D%IOR=-1
         HIT_SOLID = .TRUE.
         STEP_FRACTION(LP%ONE_D%IOR) = MAX(0._EB,(X(IIG)-X_OLD-0.05_EB*DX(IIG))/(LP%X-X_OLD))
      ENDIF
      IF (IIN<IIG .AND. WALL(IWM1)%BOUNDARY_TYPE==SOLID_BOUNDARY) THEN
         LP%ONE_D%IOR= 1
         HIT_SOLID = .TRUE.
         STEP_FRACTION(LP%ONE_D%IOR) = MAX(0._EB,(X(IIG-1)-X_OLD+0.05_EB*DX(IIG-1))/(LP%X-X_OLD))
      ENDIF
      IF (JJN>JJG .AND. WALL(IWP2)%BOUNDARY_TYPE==SOLID_BOUNDARY) THEN
         LP%ONE_D%IOR=-2
         HIT_SOLID = .TRUE.
         STEP_FRACTION(LP%ONE_D%IOR) = MAX(0._EB,(Y(JJG)-Y_OLD-0.05_EB*DY(JJG))/(LP%Y-Y_OLD))
      ENDIF
      IF (JJN<JJG .AND. WALL(IWM2)%BOUNDARY_TYPE==SOLID_BOUNDARY) THEN
         LP%ONE_D%IOR= 2
         HIT_SOLID = .TRUE.
         STEP_FRACTION(LP%ONE_D%IOR) = MAX(0._EB,(Y(JJG-1)-Y_OLD+0.05_EB*DY(JJG-1))/(LP%Y-Y_OLD))
      ENDIF

      ! Remove the particle if it is not allowed on a surface

      IF (LP%ONE_D%IOR/=0 .AND. .NOT.ALLOW_SURFACE_PARTICLES) THEN
         LP%ONE_D%X(1) = 0.9_EB*LPC%KILL_RADIUS
         CYCLE PARTICLE_LOOP
      ENDIF

      ! Get the wall index of the surface

      IML = MINLOC(STEP_FRACTION,DIM=1)
      IOR_FIRST = 0
      SELECT CASE(IML)
         CASE(1)
            IOR_FIRST = -3
         CASE(2)
            IOR_FIRST = -2
         CASE(3)
            IOR_FIRST = -1
         CASE(5)
            IOR_FIRST =  1
         CASE(6)
            IOR_FIRST =  2
         CASE(7)
            IOR_FIRST =  3
      END SELECT
      LP%WALL_INDEX = WALL_INDEX(IC,-IOR_FIRST)

      ! If no solid boundaries of original cell have been crossed, check boundaries of new grid cell
 
      IF (LP%WALL_INDEX==0) THEN
         IWP1 = WALL_INDEX(ICN, 1)
         IWM1 = WALL_INDEX(ICN,-1)
         IWP2 = WALL_INDEX(ICN, 2)
         IWM2 = WALL_INDEX(ICN,-2)
         IWP3 = WALL_INDEX(ICN, 3)
         IWM3 = WALL_INDEX(ICN,-3)
         HIT_SOLID = .FALSE.
         STEP_FRACTION = 1._EB
         IF (KKN>KKG .AND. WALL(IWM3)%BOUNDARY_TYPE==SOLID_BOUNDARY) THEN
            LP%ONE_D%IOR=-3
            HIT_SOLID = .TRUE.
            STEP_FRACTION(LP%ONE_D%IOR) = MAX(0._EB,(Z(KKG)-Z_OLD-0.05_EB*DZ(KKG))/(LP%Z-Z_OLD))
         ENDIF
         IF (KKN<KKG .AND. WALL(IWP3)%BOUNDARY_TYPE==SOLID_BOUNDARY) THEN
            LP%ONE_D%IOR= 3
            HIT_SOLID = .TRUE.
            STEP_FRACTION(LP%ONE_D%IOR) = MAX(0._EB,(Z(KKG-1)-Z_OLD+0.05_EB*DZ(KKG-1))/(LP%Z-Z_OLD))
         ENDIF
         IF (IIN>IIG .AND. WALL(IWM1)%BOUNDARY_TYPE==SOLID_BOUNDARY) THEN
            LP%ONE_D%IOR=-1
            HIT_SOLID = .TRUE.
            STEP_FRACTION(LP%ONE_D%IOR) = MAX(0._EB,(X(IIG)-X_OLD-0.05_EB*DX(IIG))/(LP%X-X_OLD))
         ENDIF
         IF (IIN<IIG .AND. WALL(IWP1)%BOUNDARY_TYPE==SOLID_BOUNDARY) THEN
            LP%ONE_D%IOR= 1
            HIT_SOLID = .TRUE.
            STEP_FRACTION(LP%ONE_D%IOR) = MAX(0._EB,(X(IIG-1)-X_OLD+0.05_EB*DX(IIG-1))/(LP%X-X_OLD))
         ENDIF
         IF (JJN>JJG .AND. WALL(IWM2)%BOUNDARY_TYPE==SOLID_BOUNDARY) THEN
            LP%ONE_D%IOR=-2
            HIT_SOLID = .TRUE.
            STEP_FRACTION(LP%ONE_D%IOR) = MAX(0._EB,(Y(JJG)-Y_OLD-0.05_EB*DY(JJG))/(LP%Y-Y_OLD))
         ENDIF
         IF (JJN<JJG .AND. WALL(IWP2)%BOUNDARY_TYPE==SOLID_BOUNDARY) THEN
            LP%ONE_D%IOR= 2
            HIT_SOLID = .TRUE.
            STEP_FRACTION(LP%ONE_D%IOR) = MAX(0._EB,(Y(JJG-1)-Y_OLD+0.05_EB*DY(JJG-1))/(LP%Y-Y_OLD))
         ENDIF

         IML = MINLOC(STEP_FRACTION,DIM=1)
         IOR_FIRST = 0
         SELECT CASE(IML)
            CASE(1)
               IOR_FIRST = -3
            CASE(2)
               IOR_FIRST = -2
            CASE(3)
               IOR_FIRST = -1
            CASE(5)
               IOR_FIRST =  1
            CASE(6)
               IOR_FIRST =  2
            CASE(7)
               IOR_FIRST =  3
         END SELECT
         LP%WALL_INDEX = WALL_INDEX(ICN,IOR_FIRST)
      ENDIF

      ! Special case where the particle is inside a solid but cannot be processed for some reason. Just return the particle to its
      ! starting position and set its velocity to 0.

      IF (SOLID(ICN) .AND. .NOT.HIT_SOLID) THEN
         LP%X = X_OLD
         LP%Y = Y_OLD
         LP%Z = Z_OLD
         LP%U = 0._EB
         LP%V = 0._EB
         LP%W = 0._EB
         CYCLE PARTICLE_LOOP
      ENDIF

      ! Check if PARTICLE has crossed no solid planes or too many

      IF_HIT_SOLID: IF (HIT_SOLID) THEN

         IF (LP%WALL_INDEX==0) CYCLE PARTICLE_LOOP

         ! Add PARTICLE mass to accumulated liquid array

         IF (ACCUMULATE_WATER .AND. HIT_SOLID .AND. .NOT.LP%SPLAT .AND. LPC%LIQUID_DROPLET) THEN
            WALL(LP%WALL_INDEX)%A_LP_MPUA(LPC%ARRAY_INDEX) = WALL(LP%WALL_INDEX)%A_LP_MPUA(LPC%ARRAY_INDEX)+&
               LP%PWT*LPC%FTPR*R_D**3*WALL(LP%WALL_INDEX)%RAW
            LP%SPLAT = .TRUE.
         ENDIF

         ! Adjust the size of the PARTICLE and weighting factor 

         IF (LPC%LIQUID_DROPLET) THEN
            R_D = MIN(0.5_EB*SURFACE_PARTICLE_DIAMETER,LP%PWT**ONTH*R_D)
            LP%PWT = LP%PWT*(R_D_0/R_D)**3
            LP%ONE_D%X(1) = R_D
            LP%ONE_D%LAYER_THICKNESS(1) = R_D
            LP%MASS = FOTHPI*LP%ONE_D%RHO(1,1)*R_D**3
         ENDIF
         
         ! Move particle to where it almost hits solid

         LP%X = X_OLD + MINVAL(STEP_FRACTION)*(LP%X-X_OLD)
         LP%Y = Y_OLD + MINVAL(STEP_FRACTION)*(LP%Y-Y_OLD)
         LP%Z = Z_OLD + MINVAL(STEP_FRACTION)*(LP%Z-Z_OLD)
         
         CALL GET_IJK(LP%X,LP%Y,LP%Z,NM,XI,YJ,ZK,LP%ONE_D%IIG,LP%ONE_D%JJG,LP%ONE_D%KKG)

         IIG = LP%ONE_D%IIG
         JJG = LP%ONE_D%JJG
         KKG = LP%ONE_D%KKG
         
         ICN = CELL_INDEX(IIG,JJG,KKG)
         IF (IOR_OLD==LP%ONE_D%IOR) CYCLE PARTICLE_LOOP

         ! Check if PARTICLE has not found surface. Simply remove for now. Todo: search algorithm

         IW = WALL_INDEX(ICN, -LP%ONE_D%IOR)
         IF (WALL(IW)%BOUNDARY_TYPE==NULL_BOUNDARY) THEN
            LP%ONE_D%X(1) = 0.9_EB*LPC%KILL_RADIUS
            CYCLE PARTICLE_LOOP
         ENDIF

         ! Choose a direction for the PARTICLEs to move

         DIRECTION: SELECT CASE(LP%ONE_D%IOR)
            CASE (-2:-1,1:2) DIRECTION  
               LP%U = 0._EB
               LP%V = 0._EB
               LP%W = -LPC%VERTICAL_VELOCITY 
            CASE (-3) DIRECTION 
               IF (.NOT.ALLOW_UNDERSIDE_PARTICLES) THEN 
                  LP%U = 0._EB
                  LP%V = 0._EB
                  LP%W = -LPC%VERTICAL_VELOCITY 
                  LP%ONE_D%IOR = 0
               ELSE 
                  CALL RANDOM_NUMBER(RN)
                  THETA_RN = TWOPI*REAL(RN,EB)
                  LP%U = LPC%HORIZONTAL_VELOCITY*COS(THETA_RN)
                  LP%V = LPC%HORIZONTAL_VELOCITY*SIN(THETA_RN)
                  LP%W = 0._EB
               ENDIF
            CASE (3) DIRECTION
               CALL RANDOM_NUMBER(RN)
               THETA_RN = TWOPI*REAL(RN,EB)
               LP%U = LPC%HORIZONTAL_VELOCITY*COS(THETA_RN)
               LP%V = LPC%HORIZONTAL_VELOCITY*SIN(THETA_RN)
               LP%W = 0._EB
         END SELECT DIRECTION

         ! If the particle is solid, as opposed to a liquid droplet, do not make it stick to the wall.

         IF (LPC%SOLID_PARTICLE) LP%ONE_D%IOR = 0

      ENDIF IF_HIT_SOLID

   ENDIF AIR_TO_SOLID 

   ! Check if PARTICLEs that were attached to a solid are still attached after the time update

   IW = WALL_INDEX(ICN, -LP%ONE_D%IOR)

   IF (WALL(IW)%BOUNDARY_TYPE/=SOLID_BOUNDARY) THEN
      SELECT CASE(LP%ONE_D%IOR)
         CASE( 1)
            LP%X = LP%X - 0.2_EB*DX(IIG)
            LP%W = -LP%W
         CASE(-1)
            LP%X = LP%X + 0.2_EB*DX(IIG)
            LP%W = -LP%W
         CASE( 2)
            LP%Y = LP%Y - 0.2_EB*DY(JJG)
            LP%W = -LP%W
         CASE(-2)
            LP%Y = LP%Y + 0.2_EB*DY(JJG)
            LP%W = -LP%W
         CASE( 3) ! Particle has reached the edge of a horizontal surface
            LP%U = -LP%U
            LP%V = -LP%V
            LP%Z =  LP%Z - 0.2_EB*DZ(KKG)
         CASE(-3)
      END SELECT
   ENDIF

   IF (LP%ONE_D%IOR/=0 .AND. WALL(IW)%BOUNDARY_TYPE/=SOLID_BOUNDARY) THEN
      LP%ONE_D%IOR = 0
      LP%WALL_INDEX = 0
   ELSE
      LP%WALL_INDEX = WALL_INDEX(ICN,-LP%ONE_D%IOR)
   ENDIF
   
   CALL GET_IJK(LP%X,LP%Y,LP%Z,NM,XI,YJ,ZK,LP%ONE_D%IIG,LP%ONE_D%JJG,LP%ONE_D%KKG)         

ENDDO PARTICLE_LOOP


! Remove out-of-bounds particles

CALL REMOVE_PARTICLES(T,NM)

CONTAINS


SUBROUTINE MOVE_ON_SOLID

! Move particles attached to solid surfaces

LP%ACCEL_X = 0._EB 
LP%ACCEL_Y = 0._EB 
LP%ACCEL_Z = 0._EB 
LP%X = LP%X + LP%U*DT
LP%Y = LP%Y + LP%V*DT
LP%Z = LP%Z + LP%W*DT

END SUBROUTINE MOVE_ON_SOLID


SUBROUTINE MOVE_IN_GAS

USE PHYSICAL_FUNCTIONS, ONLY : DRAG, GET_VISCOSITY
USE MATH_FUNCTIONS, ONLY : AFILL2, RANDOM_CHOICE, BOX_MULLER
REAL(EB) UBAR,VBAR,WBAR,RVC,UREL,VREL,WREL,QREL,RHO_G,TMP_G,MU_AIR, &
         U_OLD,V_OLD,W_OLD,ZZ_GET(1:N_TRACKED_SPECIES),WAKE_VEL,DROP_DEN,DROP_VOL_FRAC,RE_WAKE,&
         WE_G,T_BU_BAG,T_BU_STRIP,FP_MASS,HALF_DT2,BETA,OBDT,ALPHA,OPA,DTOPA,BDTOA,MPOM,ALBO,SFAC,BREAKUP_RADIUS(0:NDC),&
         DD,DD_X,DD_Y,DD_Z,DW_X,DW_Y,DW_Z,K_TERM(3),Y_TERM(3),C_DRAG,A_DRAG,HAB,PARACOR,QREL2,X_WGT,Y_WGT,Z_WGT
INTEGER IIX,JJY,KKZ

ZZ_GET = 0._EB

CALL GET_IJK(LP%X,LP%Y,LP%Z,NM,XI,YJ,ZK,LP%ONE_D%IIG,LP%ONE_D%JJG,LP%ONE_D%KKG)

U_OLD = LP%U
V_OLD = LP%V
W_OLD = LP%W
   
! Interpolate the nearest velocity components of the gas

IIX  = FLOOR(XI+.5_EB)
JJY  = FLOOR(YJ+.5_EB)
KKZ  = FLOOR(ZK+.5_EB)

X_WGT = XI+.5_EB-IIX
Y_WGT = YJ+.5_EB-JJY
Z_WGT = ZK+.5_EB-KKZ
IF (X_WGT>=0.5_EB .AND. WALL_INDEX(IC,-1)>0) X_WGT = 1._EB
IF (X_WGT< 0.5_EB .AND. WALL_INDEX(IC, 1)>0) X_WGT = 0._EB
IF (Y_WGT>=0.5_EB .AND. WALL_INDEX(IC,-2)>0) Y_WGT = 1._EB
IF (Y_WGT< 0.5_EB .AND. WALL_INDEX(IC, 2)>0) Y_WGT = 0._EB
IF (Z_WGT>=0.5_EB .AND. WALL_INDEX(IC,-3)>0) Z_WGT = 1._EB
IF (Z_WGT< 0.5_EB .AND. WALL_INDEX(IC, 3)>0) Z_WGT = 0._EB
UBAR = AFILL2(U,IIG-1,JJY,KKZ,(LP%X-X(IIG-1))*RDX(IIG),Y_WGT,Z_WGT)
VBAR = AFILL2(V,IIX,JJG-1,KKZ,X_WGT,(LP%Y-Y(JJG-1))*RDY(JJG),Z_WGT)
WBAR = AFILL2(W,IIX,JJY,KKG-1,X_WGT,Y_WGT,(LP%Z-Z(KKG-1))*RDZ(KKG))

! If the particle is massless, just move it and go on to the next particle

IF (LPC%MASSLESS_TRACER .OR. LP%PWT<=TWO_EPSILON_EB) THEN
   IF (LPC%TURBULENT_DISPERSION) THEN
      DD_X = (MU(IIG+1,JJG,KKG) - MU(IIG-1,JJG,KKG))*RSC
      DD_Y = (MU(IIG,JJG+1,KKG) - MU(IIG,JJG-1,KKG))*RSC
      DD_Z = (MU(IIG,JJG,KKG+1) - MU(IIG,JJG,KKG-1))*RSC
      LP%U = UBAR + DD_X*RDX(IIG)/RHO(IIG,JJG,KKG)
      LP%V = VBAR + DD_Y*RDY(JJG)/RHO(IIG,JJG,KKG)
      LP%W = WBAR + DD_Z*RDZ(KKG)/RHO(IIG,JJG,KKG)
      DD   = SQRT(2._EB*MU(IIG,JJG,KKG)/RHO(IIG,JJG,KKG)*RSC*DT)
      ! generate pairs of standard Gaussian random variables
      CALL BOX_MULLER(DW_X,DW_Y)
      CALL BOX_MULLER(DW_Z,DW_X)
      LP%X = LP%X + LP%U*DT + DD*DW_X
      LP%Y = LP%Y + LP%V*DT + DD*DW_Y
      LP%Z = LP%Z + LP%W*DT + DD*DW_Z
   ELSE
      LP%U = UBAR
      LP%V = VBAR
      LP%W = WBAR
      LP%X = LP%X + LP%U*DT
      LP%Y = LP%Y + LP%V*DT
      LP%Z = LP%Z + LP%W*DT
   ENDIF
   RETURN
ENDIF

! Calculate the particle drag coefficient

RVC   = RDX(IIG)*RRN(IIG)*RDY(JJG)*RDZ(KKG)
RHO_G = RHO(IIG,JJG,KKG)
UREL   = LP%U - UBAR
VREL   = LP%V - VBAR
WREL   = LP%W - WBAR
QREL   = SQRT(UREL*UREL + VREL*VREL + WREL*WREL)

DRAG_LAW_SELECT: SELECT CASE (LPC%DRAG_LAW)

   CASE (USER_DRAG)

      C_DRAG = LPC%DRAG_COEFFICIENT(1)

   CASE (SCREEN_DRAG, POROUS_DRAG)  ! Calculate this elsewhere

   CASE DEFAULT

      TMP_G  = MAX(TMPMIN,TMP(IIG,JJG,KKG))
      ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(IIG,JJG,KKG,1:N_TRACKED_SPECIES)
      CALL GET_VISCOSITY(ZZ_GET,MU_AIR,TMP_G)
      LP%RE  = RHO_G*QREL*2._EB*R_D/MU_AIR
      C_DRAG = DRAG(LP%RE,LPC%DRAG_LAW)

      ! Drag reduction model for liquid droplets

      WAKE_VEL=1.0_EB
      IF (LPC%LIQUID_DROPLET) THEN
         DROP_DEN      = AVG_DROP_DEN_ALL(IIG,JJG,KKG) 
         DROP_VOL_FRAC = MIN(1._EB,DROP_DEN/LPC%DENSITY)
         IF (DROP_VOL_FRAC > LPC%DENSE_VOLUME_FRACTION) CALL WAKE_REDUCTION(DROP_VOL_FRAC,LP%RE,C_DRAG,WAKE_VEL)
      ENDIF

      ! Secondary break-up model

      BREAKUP: IF (LPC%BREAKUP) THEN
         ! Use undisturbed wake velocity for breakup calculations
         WAKE_VEL    = WAKE_VEL*QREL
         RE_WAKE     = RHO_G*WAKE_VEL   *2._EB*R_D/MU_AIR
         WE_G        = RHO_G*WAKE_VEL**2*2._EB*R_D/LPC%SURFACE_TENSION
         ! Shape Deformation
         C_DRAG = SHAPE_DEFORMATION(RE_WAKE,WE_G,C_DRAG)
         ! Breakup conditions according to WAVE model by Reitz (1987)
         T_BU_BAG    = T_END-T_BEGIN
         T_BU_STRIP  = T_END-T_BEGIN
         IF (WE_G >= 12.0_EB)               T_BU_BAG   = 1.72_EB*B_1*SQRT(LPC%DENSITY*R_D**3/(2._EB*LPC%SURFACE_TENSION))
         IF (WE_G/SQRT(RE_WAKE) >= 1.0_EB)  T_BU_STRIP = B_1*(R_D/WAKE_VEL)*SQRT(LPC%DENSITY/RHO_G)
         ! PARTICLE age is larger than smallest characteristic breakup time
         AGE_IF: IF ((T-LP%T_INSERT) > MIN(T_BU_BAG,T_BU_STRIP)) THEN
            IF (LPC%MONODISPERSE) THEN
               R_D = THROHALF*R_D
            ELSE
               DO WHILE (R_D >= R_D_0)
                  BREAKUP_RADIUS = LPC%BREAKUP_RATIO*R_D_0*LPC%BREAKUP_R_CNF(:)
                  CALL RANDOM_CHOICE(LPC%BREAKUP_CNF(:),BREAKUP_RADIUS,NDC,R_D)
               END DO
               R_D = MAX(R_D,1.1_EB*LPC%MINIMUM_DIAMETER/2._EB)
            ENDIF
            LP%RE    = RHO_G*QREL*2._EB*R_D/MU_AIR
            C_DRAG   = DRAG(LP%RE,LPC%DRAG_LAW)
            LP%PWT   = LP%PWT*(R_D_0/R_D)**3
            LP%T_INSERT = T
            LP%ONE_D%X(1) = R_D
            LP%ONE_D%LAYER_THICKNESS(1) = R_D
            LP%MASS = FOTHPI*LP%ONE_D%RHO(1,1)*R_D**3
            ! Redo wake reduction and shape deformation for the new drop
            ! Drag reduction, except for particles associated with a SURF line
            WAKE_VEL = 1.0_EB
            IF (LPC%LIQUID_DROPLET) THEN
               DROP_DEN      = AVG_DROP_DEN(IIG,JJG,KKG,LPC%ARRAY_INDEX) 
               DROP_VOL_FRAC = MIN(1._EB,DROP_DEN/LPC%DENSITY)
               IF (DROP_VOL_FRAC > LPC%DENSE_VOLUME_FRACTION) CALL WAKE_REDUCTION(DROP_VOL_FRAC,LP%RE,C_DRAG,WAKE_VEL)
            ENDIF
            ! Change in drag coefficient due to deformation of PARTICLE shape (WE_G > 2)
            WAKE_VEL = WAKE_VEL*QREL
            RE_WAKE  = RHO_G*WAKE_VEL   *2._EB*R_D/MU_AIR
            WE_G     = RHO_G*WAKE_VEL**2*2._EB*R_D/LPC%SURFACE_TENSION
            ! Shape Deformation
            C_DRAG   = SHAPE_DEFORMATION(RE_WAKE,WE_G,C_DRAG)
         ENDIF AGE_IF
      ENDIF BREAKUP

END SELECT DRAG_LAW_SELECT

! Calculate the cross-sectional area of the droplet or particle

IF (LPC%DRAG_LAW/=SCREEN_DRAG .AND. LPC%DRAG_LAW/=POROUS_DRAG) THEN
   IF (LPC%LIQUID_DROPLET) THEN
      A_DRAG = PI*R_D**2
   ELSE
      SELECT CASE(SF%GEOMETRY)
         CASE(SURF_CARTESIAN)
            A_DRAG = SF%LENGTH*SF%WIDTH
         CASE(SURF_CYLINDRICAL)
            A_DRAG = 2._EB*R_D*SF%LENGTH
         CASE(SURF_SPHERICAL)
            A_DRAG = PI*R_D**2
      END SELECT
   ENDIF
ENDIF

! Move the particles unless they are STATIC

PARTICLE_NON_STATIC_IF: IF (.NOT.LPC%STATIC) THEN ! Move airborne, non-stationary particles

   FP_MASS = (RHO_G/RVC)/NDPC(IIG,JJG,KKG) ! fluid parcel mass
   IF (FREEZE_VELOCITY) FP_MASS = 1.E10_EB
   
   HALF_DT2 = 0.5_EB*DT*DT
   BETA  = 0.5_EB*RHO_G*C_DRAG*A_DRAG*(1._EB/LP%MASS+1._EB/FP_MASS)*QREL
   OBDT  = 1._EB+BETA*DT
   ALPHA = FP_MASS/LP%MASS
   OPA   = 1._EB+ALPHA
   DTOPA = DT/OPA
   BDTOA = BETA*DTOPA

   LP%U = ( U_OLD + (U_OLD+ALPHA*UBAR)*BDTOA )/OBDT
   LP%V = ( V_OLD + (V_OLD+ALPHA*VBAR)*BDTOA )/OBDT
   LP%W = ( W_OLD + (W_OLD+ALPHA*WBAR)*BDTOA )/OBDT
            
   IF (BETA>TWO_EPSILON_EB) THEN
      ! fluid momentum source term
      MPOM = LP%PWT*LP%MASS/(RHO_G/RVC)
      LP%ACCEL_X = MPOM*(U_OLD-LP%U)*RDT 
      LP%ACCEL_Y = MPOM*(V_OLD-LP%V)*RDT
      LP%ACCEL_Z = MPOM*(W_OLD-LP%W)*RDT
      ! semi-analytical solution for PARTICLE position
      ALBO = ALPHA*LOG(OBDT)/BETA/OPA
      LP%X = X_OLD + (U_OLD+ALPHA*UBAR)*DTOPA + ALBO*(U_OLD-UBAR) + GVEC(1)*HALF_DT2
      LP%Y = Y_OLD + (V_OLD+ALPHA*VBAR)*DTOPA + ALBO*(V_OLD-VBAR) + GVEC(2)*HALF_DT2
      LP%Z = Z_OLD + (W_OLD+ALPHA*WBAR)*DTOPA + ALBO*(W_OLD-WBAR) + GVEC(3)*HALF_DT2
   ELSE
      ! no drag
      LP%ACCEL_X  = 0._EB
      LP%ACCEL_Y  = 0._EB
      LP%ACCEL_Z  = 0._EB
      LP%X = X_OLD + DT*U_OLD + GVEC(1)*HALF_DT2
      LP%Y = Y_OLD + DT*V_OLD + GVEC(2)*HALF_DT2
      LP%Z = Z_OLD + DT*W_OLD + GVEC(3)*HALF_DT2
   ENDIF
   
   ! gravitational acceleration terms
   
   LP%U = LP%U + GVEC(1)*DT
   LP%V = LP%V + GVEC(2)*DT
   LP%W = LP%W + GVEC(3)*DT
   
   ! 2nd-order terms for the particle relative velocities (experimental)

   IF (SECOND_ORDER_PARTICLE_TRANSPORT .AND. DT>BETA*HALF_DT2) THEN
      HAB = ALPHA*BETA*HALF_DT2/OPA
      QREL2=QREL*QREL
      IF (QREL2<TWO_EPSILON_EB) THEN
         PARACOR = 0._EB
      ELSE
         PARACOR = (UREL*GVEC(1) + VREL*GVEC(2) + WREL*GVEC(3))/QREL2
      ENDIF
      LP%U = LP%U - HAB*(GVEC(1) + UREL*PARACOR)
      LP%V = LP%V - HAB*(GVEC(2) + VREL*PARACOR)
      LP%W = LP%W - HAB*(GVEC(3) + WREL*PARACOR)
   ENDIF

   IF (PARTICLE_CFL) PART_UVWMAX = MAX(PART_UVWMAX,MAX( ABS(LP%U)*RDX(IIG),ABS(LP%V)*RDY(JJG),ABS(LP%W)*RDZ(KKG)))
            
ELSE PARTICLE_NON_STATIC_IF ! Drag calculation for stationary, airborne particles

   SELECT CASE (LPC%DRAG_LAW)
      CASE DEFAULT
         BETA = 0.5_EB*LP%PWT*RVC*C_DRAG*A_DRAG*QREL
         LP%ACCEL_X = -UBAR*BETA
         LP%ACCEL_Y = -VBAR*BETA
         LP%ACCEL_Z = -WBAR*BETA
      CASE (SCREEN_DRAG)
         IF (QREL > 0.015_EB .AND. LPC%FREE_AREA_FRACTION < 1.0_EB ) THEN ! Testing shows below this can have instability  
            TMP_G  = MAX(TMPMIN,TMP(IIG,JJG,KKG))
            ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(IIG,JJG,KKG,1:N_TRACKED_SPECIES)
            CALL GET_VISCOSITY(ZZ_GET,MU_AIR,TMP_G)
            Y_TERM = LPC%DRAG_COEFFICIENT * RHO_G /SQRT(LPC%PERMEABILITY)*QREL*ABS(ORIENTATION_VECTOR(1:3,LPC%ORIENTATION_INDEX))
            K_TERM = MU_AIR/LPC%PERMEABILITY*ABS(ORIENTATION_VECTOR(1:3,LPC%ORIENTATION_INDEX))
            SFAC = 2._EB*MAXVAL(LP%ONE_D%X(0:SF%N_CELLS_MAX))*RVC/RHO_G
            LP%ACCEL_X = -(K_TERM(1)+Y_TERM(1))*UBAR*DY(JJG)*DZ(KKG)*SFAC     
            LP%ACCEL_Y = -(K_TERM(2)+Y_TERM(2))*VBAR*DX(IIG)*DZ(KKG)*SFAC
            LP%ACCEL_Z = -(K_TERM(3)+Y_TERM(3))*WBAR*DX(IIG)*DY(JJG)*SFAC
         ELSE
            LP%ACCEL_X = 0._EB 
            LP%ACCEL_Y = 0._EB
            LP%ACCEL_Z = 0._EB
         ENDIF
      CASE (POROUS_DRAG)
         TMP_G  = MAX(TMPMIN,TMP(IIG,JJG,KKG))
         ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(IIG,JJG,KKG,1:N_TRACKED_SPECIES)
         CALL GET_VISCOSITY(ZZ_GET,MU_AIR,TMP_G)
         Y_TERM = LPC%DRAG_COEFFICIENT * RHO_G /SQRT(LPC%PERMEABILITY)*QREL
         K_TERM = MU_AIR/LPC%PERMEABILITY
         SFAC = 2._EB*RHO_G
         LP%ACCEL_X = -(K_TERM(1)+Y_TERM(1))*UBAR*SFAC
         LP%ACCEL_Y = -(K_TERM(2)+Y_TERM(2))*VBAR*SFAC
         LP%ACCEL_Z = -(K_TERM(3)+Y_TERM(3))*WBAR*SFAC
   END SELECT
ENDIF PARTICLE_NON_STATIC_IF

END SUBROUTINE MOVE_IN_GAS


SUBROUTINE WAKE_REDUCTION(DROP_VOL_FRAC,RE,C_DRAG,WAKE_VEL)

! Compute C_DRAG reduction due to the wake effect (Ramirez, Munoz et al. 2007)

REAL(EB), INTENT(INOUT) :: C_DRAG
REAL(EB), INTENT(IN)  :: DROP_VOL_FRAC,RE
REAL(EB), INTENT(OUT) :: WAKE_VEL 
REAL(EB) :: LODM,RELOD

LODM = (PI/(6._EB*DROP_VOL_FRAC))**(1./3.) - 0.5_EB
RELOD = RE/(16._EB * LODM)
WAKE_VEL = 1._EB - 0.5_EB*C_DRAG*(1._EB - EXP(-RELOD))
WAKE_VEL = MAX(WAKE_VEL,0.15_EB)
C_DRAG = C_DRAG * WAKE_VEL * (1._EB + (RELOD/LODM)*EXP(-RELOD))

END SUBROUTINE WAKE_REDUCTION


REAL(EB) FUNCTION SHAPE_DEFORMATION(RE,WE,C_DRAG)

! SHAPE DEFORMATION Loth, 2008 
! E.Loth, Quasi-steady shape and drag of deformable bubbles and drops, International Journal of Multiphase Flow 34 (2008)

REAL(EB):: RE,WE,C_DRAG,C_DRAGNEW,E
REAL(EB):: DC_DRAGSTAR,fSN,WERE02

IF (WE>2.0_EB) THEN
    WERE02=WE*RE**0.2_EB
    DC_DRAGSTAR=.38E-2_EB*WERE02+3.E-5_EB*WERE02**2+9.E-7_EB*WERE02**3
    fSN=1.0_EB+0.15_EB*RE**0.687_EB
    C_DRAGNEW=1.0_EB/(3.0_EB*RE)*(DC_DRAGSTAR*(8._EB*RE+72._EB-72._EB*fSN)+72._EB*fSN)
    C_DRAGNEW=MIN(8.0_EB/3.0_EB,C_DRAGNEW) ! Bounded from above by drag of a disintegrating drop
    C_DRAGNEW=MAX(C_DRAG,C_DRAGNEW)
    ! Absorb the effect of the larger projected surface area into C_DRAG.
    ! Particle movement routines use projected area of a sphere,
    ! calculate the ratio of projected surface areas of a sphere and an
    ! ellipsoid of the same volume with aspect ratio E.
    E=1._EB-0.75_EB*TANH(0.07_EB*WE)
    SHAPE_DEFORMATION=C_DRAGNEW*E**(-2.0_EB/3.0_EB)
ELSE
    SHAPE_DEFORMATION=C_DRAG
ENDIF    

END FUNCTION SHAPE_DEFORMATION
     
END SUBROUTINE MOVE_PARTICLES


SUBROUTINE PARTICLE_MASS_ENERGY_TRANSFER(T,NM)
    
! Mass and energy transfer between gas and PARTICLEs

USE PHYSICAL_FUNCTIONS, ONLY : GET_MASS_FRACTION,GET_AVERAGE_SPECIFIC_HEAT,GET_MOLECULAR_WEIGHT,GET_SPECIFIC_GAS_CONSTANT,&
                               GET_SPECIFIC_HEAT,GET_MASS_FRACTION_ALL,GET_SENSIBLE_ENTHALPY,GET_VISCOSITY
USE MATH_FUNCTIONS, ONLY: INTERPOLATE1D_UNIFORM
USE COMP_FUNCTIONS, ONLY: SHUTDOWN
USE OUTPUT_DATA, ONLY: M_DOT,Q_DOT
USE TRAN, ONLY: GET_IJK
REAL(EB), POINTER, DIMENSION(:,:,:) :: DROP_DEN=>NULL(),DROP_RAD=>NULL(),DROP_TMP=>NULL(),MVAP_TOT=>NULL(),DROP_AREA=>NULL(),&
                                       DROP_DEN_ALL=>NULL()
REAL(EB), POINTER, DIMENSION(:) :: FILM_THICKNESS=>NULL()
REAL(EB) :: R_DROP,NUSSELT,K_AIR,H_V,H_V_REF, H_L,&
            RVC,WGT,Q_CON_GAS,Q_CON_WALL,Q_RAD,H_HEAT,H_MASS,SH_FAC_GAS,SH_FAC_WALL,NU_FAC_GAS,NU_FAC_WALL, &
            PR_AIR,M_VAP,M_VAP_MAX,MU_AIR,H_SOLID,Q_DOT_RAD,DEN_ADD,AREA_ADD, &
            Y_DROP,Y_GAS,LENGTH,U2,V2,W2,VEL,DENOM,DZ_DTMP_DROP,TMP_DROP_NEW,TMP_WALL,H_WALL, &
            SC_AIR,D_AIR,DHOR,SHERWOOD,X_DROP,M_DROP,RHO_G,MW_RATIO,MW_DROP,FTPR,&
            C_DROP,M_GAS,A_DROP,TMP_G,TMP_DROP,TMP_MELT,TMP_BOIL,MINIMUM_FILM_THICKNESS,RE_L,OMRAF,Q_FRAC,Q_TOT,DT_SUBSTEP, &
            CP,H_NEW,ZZ_GET(1:N_TRACKED_SPECIES),ZZ_GET2(1:N_TRACKED_SPECIES), &
            M_GAS_NEW,MW_GAS,CP2,DELTA_H_G,TMP_G_I,H_G_OLD,H_D_OLD, &
            H_L_REF,TMP_G_NEW,DT_SUM,DCPDT,TMP_WGT,X_EQUIL,Y_EQUIL,Y_ALL(1:N_SPECIES),H_S_B,H_S,C_DROP2,&
            T_BOIL_EFF,P_RATIO,K_LIQUID,CP_LIQUID,MU_LIQUID,NU_LIQUID,BETA_LIQUID,PR_LIQUID,RAYLEIGH,GR
            
INTEGER :: IP,II,JJ,KK,IW,N_LPC,NS,N_SUBSTEPS,ITMP,ITMP2,ITCOUNT,Y_INDEX,Z_INDEX,I_BOIL,I_MELT,I_FUEL
REAL(EB), INTENT(IN) :: T
INTEGER, INTENT(IN) :: NM
LOGICAL :: TEMPITER
CHARACTER(MESSAGE_LENGTH) :: MESSAGE
TYPE (LAGRANGIAN_PARTICLE_TYPE), POINTER :: LP=>NULL()
TYPE (LAGRANGIAN_PARTICLE_CLASS_TYPE), POINTER :: LPC=>NULL()
TYPE (SURFACE_TYPE), POINTER :: SF=>NULL()
TYPE (SPECIES_TYPE), POINTER :: SS=>NULL()

CALL POINT_TO_MESH(NM)

! Initializations

OMRAF  = 1._EB - RUN_AVG_FAC
M_DOT(2,NM) = 0._EB ! Fuel mass loss rate of droplets
M_DOT(4,NM) = 0._EB ! Total mass loss rate of droplets

! Rough estimates

MINIMUM_FILM_THICKNESS = 1.E-5_EB   ! Minimum thickness of liquid film on the surface (m)
H_SOLID                = 300._EB    ! Heat transfer coefficient from solid surface to drop (W/m2/K)

! Empirical coefficients

D_AIR                  = 2.6E-5_EB  ! Water Vapor - Air binary diffusion (m2/s at 25 C, Incropera & DeWitt, Table A.8) 
SC_AIR                 = SC         ! Can be set on MISC line, otherwise dependent on LES/DNS mode
PR_AIR                 = PR         ! Can be set on MISC line, otherwise dependent on LES/DNS mode
SH_FAC_GAS             = 0.6_EB*SC_AIR**ONTH
NU_FAC_GAS             = 0.6_EB*PR_AIR**ONTH        
SH_FAC_WALL            = 0.037_EB*SC_AIR**ONTH
NU_FAC_WALL            = 0.037_EB*PR_AIR**ONTH        

! Liquid properties (temporary place)
MU_LIQUID              = 1040E-6_EB  ! Ns/m2
BETA_LIQUID            = 0.18E-3_EB  ! Volumetric thermal expansion 1/K
K_LIQUID               = 0.60_EB     ! Conductivity W/(m.K)
CP_LIQUID              = 4.19E3_EB   ! J/(kg.K)
PR_LIQUID              = MU_LIQUID*CP_LIQUID/K_LIQUID

! Working arrays

IF (N_LP_ARRAY_INDICES>0) THEN
   MVAP_TOT => WORK7   
   MVAP_TOT = 0._EB
   DO IW = 1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
      WALL(IW)%LP_CPUA  = RUN_AVG_FAC*WALL(IW)%LP_CPUA
      WALL(IW)%LP_MPUA  = RUN_AVG_FAC*WALL(IW)%LP_MPUA
   ENDDO
ENDIF

! Loop over all types of evaporative species

SPECIES_LOOP: DO Z_INDEX = 1,N_TRACKED_SPECIES
   
   ! Initialize quantities common to the evaporation index

   IF (.NOT. SPECIES_MIXTURE(Z_INDEX)%EVAPORATING) CYCLE SPECIES_LOOP
   Y_INDEX = MAXVAL(MAXLOC(SPECIES_MIXTURE(Z_INDEX)%VOLUME_FRACTION))
   SS => SPECIES(Y_INDEX)
   TMP_MELT = SS%TMP_MELT
   TMP_BOIL = SS%TMP_V
   MW_DROP  = SS%MW
   CALL INTERPOLATE1D_UNIFORM(LBOUND(SS%H_V,1),SS%H_V,SS%H_V_REFERENCE_TEMPERATURE,H_V_REF)
   CALL INTERPOLATE1D_UNIFORM(LBOUND(SS%C_P_L_BAR,1),SS%C_P_L_BAR,SS%TMP_MELT,H_L_REF)
   H_L_REF= H_L_REF*TMP_MELT
   I_MELT   = INT(TMP_MELT)
   FILM_THICKNESS => WALL_WORK2
   FILM_THICKNESS =  0._EB

   ! Loop through all PARTICLEs in the class and determine the depth of the liquid film on each surface cell

   FILM_SUMMING_LOOP: DO IP=1,NLP
      LP  => LAGRANGIAN_PARTICLE(IP)
      IF (LP%WALL_INDEX==0) CYCLE FILM_SUMMING_LOOP
      IF (LP%ONE_D%IOR==0)  CYCLE FILM_SUMMING_LOOP
      LPC => LAGRANGIAN_PARTICLE_CLASS(LP%CLASS_INDEX)
      IF (.NOT.LPC%LIQUID_DROPLET) CYCLE FILM_SUMMING_LOOP
      IF (LPC%Z_INDEX/=Z_INDEX) CYCLE FILM_SUMMING_LOOP
      IF (LP%ONE_D%X(1)<=0._EB) CYCLE FILM_SUMMING_LOOP
      IW = LP%WALL_INDEX
      FILM_THICKNESS(IW) = FILM_THICKNESS(IW) + LP%PWT*LP%ONE_D%X(1)**3/WALL(IW)%AW
   ENDDO FILM_SUMMING_LOOP

   FILM_THICKNESS = MAX(MINIMUM_FILM_THICKNESS,FOTH*PI*FILM_THICKNESS) 

   ! Loop through all PARTICLEs within the class and determine mass/energy transfer

   PARTICLE_LOOP: DO IP=1,NLP

      LP  => LAGRANGIAN_PARTICLE(IP)
      LPC => LAGRANGIAN_PARTICLE_CLASS(LP%CLASS_INDEX)      
      IF (LPC%Z_INDEX/=Z_INDEX)     CYCLE PARTICLE_LOOP
      IF (.NOT.LPC%LIQUID_DROPLET)  CYCLE PARTICLE_LOOP
      IF (LP%ONE_D%X(1)<=0._EB)     CYCLE PARTICLE_LOOP

      ! Determine the current coordinates of the particle

      II = LP%ONE_D%IIG
      JJ = LP%ONE_D%JJG
      KK = LP%ONE_D%KKG
      RVC = RDX(II)*RRN(II)*RDY(JJ)*RDZ(KK)
      DHOR = H_V_REF*MW_DROP/R0
      P_RATIO = P_STP/PBAR(0,PRESSURE_ZONE(II,JJ,KK))
      ! Boiling temperature at current background pressure
      T_BOIL_EFF = MAX(0._EB,DHOR*TMP_BOIL/(DHOR-TMP_BOIL*LOG(1._EB/P_RATIO)+TWO_EPSILON_EB))
      I_BOIL   = INT(T_BOIL_EFF)
      ! Determine how many sub-time step iterations are needed and then iterate over the time step.
      ! This is not fully functional. Keep as a placeholder for now.

      N_SUBSTEPS = 1
      DT_SUBSTEP = DT/REAL(N_SUBSTEPS,EB) 
      DT_SUM = 0._EB

      TIME_ITERATION_LOOP: DO WHILE (DT_SUM < DT)

         ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(II,JJ,KK,1:N_TRACKED_SPECIES)
         CALL GET_MASS_FRACTION_ALL(ZZ_GET,Y_ALL)
         MW_GAS = 0._EB
         IF (ABS(Y_ALL(Y_INDEX)-1._EB) > TWO_EPSILON_EB) THEN
            DO NS=1,N_SPECIES
               IF (NS==Y_INDEX) CYCLE
               MW_GAS = MW_GAS + Y_ALL(NS)/SPECIES(NS)%MW
            ENDDO
            IF (MW_GAS<=TWO_EPSILON_EB) THEN
               MW_GAS=SPECIES_MIXTURE(1)%MW
            ELSE
               MW_GAS = (1._EB-Y_ALL(Y_INDEX))/MW_GAS
            ENDIF
         ELSE
            MW_GAS=SPECIES_MIXTURE(1)%MW
         ENDIF

         MW_RATIO = MW_GAS/MW_DROP
      
         ! Initialize PARTICLE thermophysical data

         R_DROP   = LP%ONE_D%X(1)
         FTPR     = FOTHPI * LP%ONE_D%RHO(1,1)
         M_DROP   = FTPR*R_DROP**3
         TMP_DROP = LP%ONE_D%TMP(1)
         CALL INTERPOLATE1D_UNIFORM(LBOUND(SS%H_V,1),SS%H_V,TMP_DROP,H_V)
         IF (H_V < 0._EB) THEN
            WRITE(MESSAGE,'(A,A)') 'Numerical instability in particle energy transport, H_V for ',TRIM(SS%ID)
            CALL SHUTDOWN(MESSAGE)
         ENDIF
         CALL INTERPOLATE1D_UNIFORM(LBOUND(SS%C_P_L,1),SS%C_P_L,TMP_DROP,C_DROP)         
         CALL INTERPOLATE1D_UNIFORM(LBOUND(SS%C_P_L_BAR,1),SS%C_P_L_BAR,TMP_DROP,H_L)                  
         H_L = H_L*TMP_DROP+H_L_REF
         H_D_OLD  = H_L*M_DROP
         WGT      = LP%PWT
         DHOR     = H_V*MW_DROP/R0


         ! Gas conditions

         TMP_G  = TMP(II,JJ,KK)
         RHO_G  = RHO(II,JJ,KK)
         CALL GET_VISCOSITY(ZZ_GET,MU_AIR,TMP_G)
         M_GAS  = RHO_G/RVC  
         M_VAP_MAX = (0.33_EB * M_GAS - MVAP_TOT(II,JJ,KK)) / WGT ! limit to avoid diveregence errors
         K_AIR  = CPOPR*MU_AIR
         IF (Y_INDEX>=0) THEN
            CALL GET_MASS_FRACTION(ZZ_GET,Y_INDEX,Y_GAS)
         ELSE
            Y_GAS = 0._EB
         ENDIF
         U2 = 0.5_EB*(U(II,JJ,KK)+U(II-1,JJ,KK))
         V2 = 0.5_EB*(V(II,JJ,KK)+V(II,JJ-1,KK))
         W2 = 0.5_EB*(W(II,JJ,KK)+W(II,JJ,KK-1))
      
         ! Set variables for heat transfer on solid

         SOLID_OR_GAS_PHASE: IF (LP%ONE_D%IOR/=0 .AND. LP%WALL_INDEX>0) THEN

            IW   = LP%WALL_INDEX
            A_DROP = M_DROP/(FILM_THICKNESS(IW)*LPC%DENSITY)
            TMP_WALL = WALL(IW)%ONE_D%TMP_F
            SELECT CASE(ABS(LP%ONE_D%IOR))
               CASE(1)
                  VEL = SQRT(V2**2+W2**2)
               CASE(2)
                  VEL = SQRT(U2**2+W2**2)
               CASE(3)
                  VEL = SQRT(U2**2+V2**2)
            END SELECT

            IF (.NOT.CONSTANT_H_SOLID) THEN
               LENGTH = 2._EB*R_DROP

               NU_LIQUID = MU_LIQUID / LPC%DENSITY
               !Grashoff number
               GR = MAXVAL(ABS(GVEC))*LENGTH**3*BETA_LIQUID*ABS(TMP_WALL-TMP_DROP)/(NU_LIQUID**2)
               RAYLEIGH = GR*PR_LIQUID
               DIRECTION2: SELECT CASE(LP%ONE_D%IOR)
               CASE (-2:-1,1:2) DIRECTION2
               ! Vertical boundary layers (Churchill, S.W.)
               NUSSELT = 0.68_EB+0.67_EB*RAYLEIGH**(0.25_EB)/((1._EB+(0.492_EB/PR_LIQUID)**(9._EB/16._EB))**(4._EB/9._EB))
               CASE (-3,3) DIRECTION2
               ! Horizontal, unstable boundary layers (top of hot plate)
               ! Raithby, G.D., Hollands, K.G.T. Natural convection. In Rohsenoh, W.M., 
               ! Hartnett, J.P., and Cho, Y.I. (eds.), Handbook of Heat Transfer, chapter
               ! 4. McGraw-Hill, New York, 3rd edition, 1998.
               NUSSELT = 0.560_EB*RAYLEIGH**(0.25_EB)/((1._EB+(0.492_EB/PR_LIQUID)**(9._EB/16._EB))**(4._EB/9._EB))
               END SELECT DIRECTION2
               H_WALL   = NUSSELT*K_LIQUID/LENGTH
               RE_L     = MAX(5.E5_EB,RHO_G*VEL*LENGTH/MU_AIR)
!               NUSSELT  = NU_FAC_WALL*RE_L**0.8_EB
               SHERWOOD = SH_FAC_WALL*RE_L**0.8_EB
            ELSE
               LENGTH   = 1._EB
               RE_L     = MAX(5.E5_EB,RHO_G*VEL*LENGTH/MU_AIR)
               NUSSELT  = NU_FAC_WALL*RE_L**0.8_EB
               SHERWOOD = SH_FAC_WALL*RE_L**0.8_EB
               H_WALL    = H_SOLID
            ENDIF
            H_HEAT   = NUSSELT*K_AIR/LENGTH
            H_MASS   = SHERWOOD*D_AIR/LENGTH
            Q_DOT_RAD = MIN(A_DROP,WALL(IW)%AW*DT/(WGT*DT_SUBSTEP))*WALL(IW)%ONE_D%QRADIN
            WALL(IW)%ONE_D%QRADIN = (WALL(IW)%AW*DT*WALL(IW)%ONE_D%QRADIN - WGT*DT_SUBSTEP*Q_DOT_RAD)/(WALL(IW)%AW*DT)

         ELSE SOLID_OR_GAS_PHASE

            A_DROP   = 4._EB*PI*R_DROP**2
            NUSSELT  = 2._EB + NU_FAC_GAS*SQRT(LP%RE)
            SHERWOOD = 2._EB + SH_FAC_GAS*SQRT(LP%RE)
            H_HEAT   = NUSSELT *K_AIR/(2._EB*R_DROP)
            H_MASS   = SHERWOOD*D_AIR/(2._EB*R_DROP)
            H_WALL   = 0._EB
            TMP_WALL = TMPA
            IF (AVG_DROP_DEN(II,JJ,KK,LPC%ARRAY_INDEX )>0._EB) THEN
               Q_DOT_RAD = (QR_W(II,JJ,KK)/SUM(AVG_DROP_AREA(II,JJ,KK,:)))*(A_DROP/4._EB)              
            ELSE
               Q_DOT_RAD = 0._EB
            ENDIF

         ENDIF SOLID_OR_GAS_PHASE

         ! Compute equilibrium PARTICLE vapor mass fraction, Y_DROP, and its derivative w.r.t. PARTICLE temperature

         X_DROP  = MIN(1._EB,P_RATIO*EXP(DHOR*(1._EB/TMP_BOIL-1._EB/TMP_DROP)))
         Y_DROP  = X_DROP/(MW_RATIO + (1._EB-MW_RATIO)*X_DROP)
         IF (TMP_DROP < T_BOIL_EFF) THEN
            DZ_DTMP_DROP = (MW_RATIO/(X_DROP*(1._EB-MW_RATIO)+MW_RATIO)**2)*DHOR*X_DROP/TMP_DROP**2
         ELSE
            DZ_DTMP_DROP = 0._EB
         ENDIF
         IF (Y_DROP<=Y_GAS) H_MASS = 0._EB

         ! Update the PARTICLE temperature semi_implicitly

         DENOM = 1._EB + (H_HEAT + H_WALL + H_MASS*RHO_G*H_V*DZ_DTMP_DROP)*DT_SUBSTEP*A_DROP/(2._EB*M_DROP*C_DROP) 
         TMP_DROP_NEW = (TMP_DROP + DT_SUBSTEP*( Q_DOT_RAD + &
                        A_DROP*(H_HEAT*(TMP_G   -0.5_EB*TMP_DROP) + H_WALL*(TMP_WALL-0.5_EB*TMP_DROP) -  &
                        H_MASS*RHO_G*H_V*(Y_DROP-0.5_EB*DZ_DTMP_DROP*TMP_DROP-Y_GAS))/(M_DROP*C_DROP)) ) / DENOM

         ! Compute the total amount of heat extracted from the gas, wall and radiative fields

         Q_RAD      = DT_SUBSTEP*Q_DOT_RAD
         Q_CON_GAS  = DT_SUBSTEP*A_DROP*H_HEAT*(TMP_G   -0.5_EB*(TMP_DROP+TMP_DROP_NEW))
         Q_CON_WALL = DT_SUBSTEP*A_DROP*H_WALL*(TMP_WALL-0.5_EB*(TMP_DROP+TMP_DROP_NEW))
         Q_TOT      = Q_RAD+Q_CON_GAS+Q_CON_WALL

         ! Compute the total amount of liquid evaporated

         M_VAP = DT_SUBSTEP*A_DROP*H_MASS*RHO_G*(Y_DROP+0.5_EB*DZ_DTMP_DROP*(TMP_DROP_NEW-TMP_DROP)-Y_GAS) 
         M_VAP = MAX(0._EB,MIN(M_VAP,M_DROP,M_VAP_MAX))
      
         ! Evaporate completely small PARTICLEs

         IF (R_DROP<0.5_EB*LPC%MINIMUM_DIAMETER) THEN
            M_VAP      = M_DROP
            Q_RAD      = Q_RAD * (DT-DT_SUM)/DT_SUBSTEP
            Q_CON_WALL = Q_CON_WALL * (DT-DT_SUM)/DT_SUBSTEP
            IF (Q_RAD + Q_CON_WALL > M_VAP*H_V) THEN
               Q_TOT = Q_RAD + Q_CON_WALL
               Q_RAD = Q_RAD * M_VAP*H_V/Q_TOT
               Q_CON_WALL = Q_CON_WALL * M_VAP*H_V/Q_TOT
            ENDIF
            Q_CON_GAS  = M_VAP*H_V-Q_RAD-Q_CON_WALL
            Q_TOT      = 0._EB
            DT_SUM   = DT
            !IF (Q_TOT>0._EB) THEN
            !   Q_FRAC     = M_VAP*H_V/Q_TOT 
            !   Q_CON_GAS  = Q_CON_GAS*Q_FRAC
            !   Q_CON_WALL = Q_CON_WALL*Q_FRAC
            !   Q_RAD      = Q_RAD*Q_FRAC
            !   Q_TOT      = Q_RAD+Q_CON_GAS+Q_CON_WALL
            !ENDIF
         ENDIF
         IF (M_VAP < M_DROP) THEN
            TMP_DROP_NEW = TMP_DROP + (Q_TOT - M_VAP * H_V)/(C_DROP * (M_DROP - M_VAP))
            ITMP = INT(TMP_DROP)
            ITMP2 = MIN(I_BOIL,MAX(I_MELT,INT(TMP_DROP_NEW)))
            IF (ITMP/=ITMP2) THEN
               C_DROP2 = SUM(SS%C_P_L(MIN(ITMP,ITMP2):MAX(ITMP,ITMP2)))/REAL(ABS(ITMP2-ITMP)+1,EB)
               TMP_DROP_NEW = TMP_DROP + (Q_TOT - M_VAP * H_V)/(C_DROP2 * (M_DROP - M_VAP))
            ENDIF
         ENDIF

         ! If the PARTICLE temperature drops below its freezing point, just reset it

         IF (TMP_DROP_NEW<TMP_MELT) TMP_DROP_NEW = TMP_MELT

         ! If the PARTICLE temperature reaches boiling, use only enough energy from gas to vaporize liquid

         IF (TMP_DROP_NEW>=T_BOIL_EFF) THEN  
            M_VAP  = MIN(M_VAP_MAX,M_DROP,M_VAP + (TMP_DROP_NEW - T_BOIL_EFF)*C_DROP*M_DROP/H_V)
            TMP_DROP_NEW = T_BOIL_EFF
            IF (Q_TOT>0._EB) THEN
               Q_FRAC     = M_VAP*H_V/Q_TOT
               Q_CON_GAS  = Q_CON_GAS*Q_FRAC
               Q_CON_WALL = Q_CON_WALL*Q_FRAC
               Q_RAD      = Q_RAD*Q_FRAC
               Q_TOT      = Q_RAD+Q_CON_GAS+Q_CON_WALL
            ENDIF
         ENDIF
         M_DROP = M_DROP - M_VAP
   
         ! Add fuel evaporation rate to running counter and adjust mass of evaporated fuel to account for different 
         ! Heat of Combustion between fuel PARTICLE and gas

         I_FUEL = 0
         IF (N_REACTIONS>0) I_FUEL = REACTION(1)%FUEL_SMIX_INDEX

         IF (LPC%Z_INDEX==I_FUEL .AND. I_FUEL>0) M_VAP = LPC%ADJUST_EVAPORATION*M_VAP
 
         ! Update gas temperature and determine new subtimestep

         CALL GET_AVERAGE_SPECIFIC_HEAT(ZZ_GET,CP,TMP_G)       
         
         H_G_OLD = M_GAS*CP*TMP_G         
         M_GAS_NEW = M_GAS + WGT*M_VAP
         TMP_G_NEW = TMP_G
         ITMP     = INT(TMP_DROP_NEW)
         TMP_WGT  = TMP_DROP_NEW - AINT(TMP_DROP_NEW)
         H_NEW = H_G_OLD + WGT*(M_VAP*(H_V+H_L) - Q_CON_GAS)

         IF (H_NEW > 0._EB) THEN
            ZZ_GET2 = ZZ_GET * M_GAS/M_GAS_NEW               
            ZZ_GET2(Z_INDEX) = ZZ_GET2(Z_INDEX) + WGT*M_VAP/M_GAS_NEW
            TMP_G_I = TMP_G
            TEMPITER = .TRUE.
            ITCOUNT = 0
            ITERATE_TEMP: DO WHILE (TEMPITER)
               TEMPITER=.FALSE.

               ! Compute approximation of d(cp)/dT

               CALL GET_AVERAGE_SPECIFIC_HEAT(ZZ_GET2,CP2,TMP_G_I)
               IF (TMP_G_I > 1._EB) THEN
                  CALL GET_AVERAGE_SPECIFIC_HEAT(ZZ_GET2,CP,TMP_G_I-1._EB)
                  DCPDT = CP2-CP
               ELSE
                  CALL GET_AVERAGE_SPECIFIC_HEAT(ZZ_GET2,CP,TMP_G_I+1._EB)
                  DCPDT = CP-CP2
               ENDIF

               TMP_G_I = TMP_G_I+(H_NEW-CP2*TMP_G_I*M_GAS_NEW)/(M_GAS_NEW*(CP2+TMP_G_I*DCPDT))

               IF (TMP_G_I < 0._EB) THEN
                  DT_SUBSTEP = DT_SUBSTEP * 0.5_EB
                  N_SUBSTEPS = NINT(DT/DT_SUBSTEP)
                  CYCLE TIME_ITERATION_LOOP
               ENDIF

               ITCOUNT = ITCOUNT + 1
               IF (ABS(TMP_G_NEW-TMP_G_I) > 0.5_EB) TEMPITER = .TRUE.
               IF (ITCOUNT > 10) THEN
                  TMP_G_NEW = 0.5_EB*(TMP_G_I + TMP_G_NEW)
                  EXIT ITERATE_TEMP
               ENDIF               
               TMP_G_NEW = TMP_G_I
            ENDDO ITERATE_TEMP
            TMP_G_NEW = MAX(TMP_G_NEW,TMPMIN)
         ELSE
            DT_SUBSTEP = DT_SUBSTEP * 0.5_EB
            N_SUBSTEPS = NINT(DT/DT_SUBSTEP)
            CYCLE TIME_ITERATION_LOOP
         ENDIF

         CALL INTERPOLATE1D_UNIFORM(LBOUND(SS%H_V,1),SS%H_V,TMP_DROP_NEW,H_V)
         DHOR     = H_V*MW_DROP/R0 
         X_EQUIL  = MIN(1._EB,P_RATIO*EXP(DHOR*(1._EB/TMP_BOIL-1._EB/TMP_DROP_NEW)))
         Y_EQUIL  = X_EQUIL/(MW_RATIO + (1._EB-MW_RATIO)*X_EQUIL)

         ! Limit super-saturation

         IF (Y_GAS < Y_EQUIL) THEN
            CALL GET_MASS_FRACTION(ZZ_GET2,Y_INDEX,Y_GAS)
            IF (Y_GAS/Y_EQUIL > 1.02_EB) THEN
               DT_SUBSTEP = DT_SUBSTEP * 0.5_EB            
               N_SUBSTEPS = NINT(DT/DT_SUBSTEP)
               IF (DT_SUBSTEP <= 0.00001_EB*DT) THEN
                  CALL SHUTDOWN('Numerical instability in particle energy transport, Y_EQUIL')
               ENDIF
               CYCLE TIME_ITERATION_LOOP
            ENDIF
         ENDIF
      
         ! Limit gas temperature change
      
         IF (ABS(TMP_G_NEW/TMP_G - 1._EB) > 0.05_EB) THEN
            DT_SUBSTEP = DT_SUBSTEP * 0.5_EB            
            N_SUBSTEPS = NINT(DT/DT_SUBSTEP)
            IF (DT_SUBSTEP <= 0.00001_EB*DT) THEN
               CALL SHUTDOWN('Numerical instability in particle energy transport, TMP_G')
            ENDIF
            CYCLE TIME_ITERATION_LOOP
         ENDIF

         ! Update gas cell density, temperature, and mass fractions

         RHO(II,JJ,KK) = M_GAS_NEW*RVC
         ZZ(II,JJ,KK,1:N_TRACKED_SPECIES) = ZZ_GET2(1:N_TRACKED_SPECIES)
         CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET2,RSUM(II,JJ,KK))
         TMP(II,JJ,KK) = TMP_G_NEW
    
         ! Compute contribution to the divergence

         CALL GET_SPECIFIC_HEAT(ZZ_GET,CP,TMP_G)
         H_G_OLD = CP * TMP_G * M_GAS
         ZZ_GET = 0._EB
         ZZ_GET(Z_INDEX) = 1._EB
         CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S_B,TMP_DROP)
         CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,TMP_G)
         DELTA_H_G = H_S_B - H_S
         D_LAGRANGIAN(II,JJ,KK) = D_LAGRANGIAN(II,JJ,KK) &
                                + (MW_RATIO*M_VAP/M_GAS + (M_VAP*DELTA_H_G - Q_CON_GAS)/H_G_OLD) * WGT / DT_SUBSTEP
         CALC_D_LAGRANGIAN = .TRUE.

         ! Add energy losses and gains to overall energy budget array

         Q_DOT(7,NM) = Q_DOT(7,NM) - (Q_CON_GAS + Q_CON_WALL + Q_RAD)*WGT/DT_SUBSTEP  ! Q_PART
         Q_DOT(3,NM) = Q_DOT(3,NM) + M_VAP*H_S_B*WGT/DT_SUBSTEP                       ! Q_CONV
         Q_DOT(2,NM) = Q_DOT(2,NM) + Q_RAD*WGT/DT_SUBSTEP                             ! Q_RADI
         Q_DOT(4,NM) = Q_DOT(4,NM) + Q_CON_WALL*WGT/DT_SUBSTEP                        ! Q_COND

         IF (LPC%Z_INDEX==I_FUEL .AND. I_FUEL>0) M_DOT(2,NM) = M_DOT(2,NM) + WGT*M_VAP/DT_SUBSTEP/LPC%ADJUST_EVAPORATION
         M_DOT(4,NM) = M_DOT(4,NM) + WGT*M_VAP/DT_SUBSTEP  ! Total mass loss rate         
         
         ! Keep track of total mass evaporated in cell

         MVAP_TOT(II,JJ,KK) = MVAP_TOT(II,JJ,KK) + WGT*M_VAP
      
         ! Update PARTICLE quantities

         LP%ONE_D%X(1)   = (M_DROP/FTPR)**ONTH
         LP%ONE_D%LAYER_THICKNESS(1) = LP%ONE_D%X(1)
         LP%ONE_D%TMP(1) = TMP_DROP_NEW
         LP%ONE_D%TMP_F  = TMP_DROP_NEW
         LP%MASS = M_DROP

         ! Compute surface cooling

         IF (LP%ONE_D%IOR/=0 .AND. LP%WALL_INDEX>0) &
         WALL(IW)%LP_CPUA(LPC%ARRAY_INDEX) = WALL(IW)%LP_CPUA(LPC%ARRAY_INDEX) + &
                                             OMRAF*WGT*(Q_RAD+Q_CON_WALL)*WALL(IW)%RAW/DT_SUBSTEP

         ! Get out of the loop if the PARTICLE has evaporated completely
         IF (LP%ONE_D%X(1)<=0._EB) CYCLE PARTICLE_LOOP

         DT_SUM = DT_SUM + DT_SUBSTEP
         DT_SUBSTEP = MIN(DT-DT_SUM,DT_SUBSTEP * 1.5_EB)

      ENDDO TIME_ITERATION_LOOP

   ENDDO PARTICLE_LOOP

ENDDO SPECIES_LOOP

TMP = MIN(TMPMAX,MAX(TMPMIN,TMP))

! Second loop is for summing the part quantities

SUM_PART_QUANTITIES: IF (N_LP_ARRAY_INDICES > 0) THEN

   DROP_AREA => WORK1
   DROP_DEN => WORK4
   DROP_RAD => WORK5
   DROP_TMP => WORK6
   DROP_DEN_ALL => WORK8

   PART_CLASS_SUM_LOOP: DO N_LPC = 1,N_LAGRANGIAN_CLASSES

      LPC => LAGRANGIAN_PARTICLE_CLASS(N_LPC)
      IF (LPC%MASSLESS_TRACER .OR. LPC%MASSLESS_TARGET) CYCLE PART_CLASS_SUM_LOOP

      DROP_DEN = 0._EB
      DROP_TMP = 0._EB
      DROP_RAD = 0._EB
      DROP_AREA = 0._EB
      DROP_DEN_ALL = 0._EB

      PARTICLE_LOOP_2: DO IP=1,NLP

         LP => LAGRANGIAN_PARTICLE(IP)
         IF (LP%CLASS_INDEX /= N_LPC) CYCLE PARTICLE_LOOP_2
         II = LP%ONE_D%IIG
         JJ = LP%ONE_D%JJG
         KK = LP%ONE_D%KKG
         RVC = RDX(II)*RRN(II)*RDY(JJ)*RDZ(KK)

         ! Determine the mass of the PARTICLE/particle, depending on whether the particle has a distinct SURFace type.

         IF (LPC%LIQUID_DROPLET) THEN
            R_DROP = LP%ONE_D%X(1)
            A_DROP = PI*R_DROP**2
         ELSE
            SF => SURFACE(LPC%SURF_INDEX)
            R_DROP = MAXVAL(LP%ONE_D%X(0:SF%N_CELLS_MAX))
            SELECT CASE(SF%GEOMETRY)
               CASE(SURF_CARTESIAN)
                  A_DROP = 2._EB*SF%LENGTH*SF%WIDTH
               CASE(SURF_CYLINDRICAL)
                  A_DROP = 2._EB*SF%LENGTH*R_DROP
               CASE(SURF_SPHERICAL)                  
                  A_DROP = PI*R_DROP**2
            END SELECT
         ENDIF

         ! Assign particle or PARTICLE mass to the grid cell if the particle/PARTICLE not on a surface
         
         IF (LP%ONE_D%IOR==0) THEN
            DEN_ADD  =    LP%PWT*LP%MASS * RVC
            AREA_ADD =    LP%PWT*A_DROP * RVC
            DROP_DEN(II,JJ,KK)  = DROP_DEN(II,JJ,KK)  + DEN_ADD
            DROP_TMP(II,JJ,KK)  = DROP_TMP(II,JJ,KK)  + DEN_ADD*LP%ONE_D%TMP_F
            DROP_RAD(II,JJ,KK)  = DROP_RAD(II,JJ,KK)  + AREA_ADD*R_DROP
            DROP_AREA(II,JJ,KK) = DROP_AREA(II,JJ,KK) + AREA_ADD
         ENDIF

         ! Compute surface density

         IF (LP%ONE_D%IOR/=0 .AND. LP%WALL_INDEX>0) THEN
            IW     = LP%WALL_INDEX
            R_DROP = LP%ONE_D%X(1)
            FTPR   = FOTHPI * LP%ONE_D%RHO(1,1)
            M_DROP = FTPR*R_DROP**3
            WALL(IW)%LP_MPUA(LPC%ARRAY_INDEX) = WALL(IW)%LP_MPUA(LPC%ARRAY_INDEX) + OMRAF*LP%PWT*M_DROP*WALL(IW)%RAW
         ENDIF

      ENDDO PARTICLE_LOOP_2

     ! Compute cumulative quantities for PARTICLE "clouds"

      DROP_RAD = DROP_RAD/(DROP_AREA+TWO_EPSILON_EB)
      DROP_TMP = DROP_TMP/(DROP_DEN +TWO_EPSILON_EB)
      AVG_DROP_RAD(:,:,:,LPC%ARRAY_INDEX ) = DROP_RAD
      AVG_DROP_TMP(:,:,:,LPC%ARRAY_INDEX ) = RUN_AVG_FAC*AVG_DROP_TMP(:,:,:,LPC%ARRAY_INDEX ) + OMRAF*DROP_TMP
      IF (LPC%Y_INDEX>0) THEN
         AVG_DROP_TMP(:,:,:,LPC%ARRAY_INDEX ) = MAX(SPECIES(LPC%Y_INDEX)%MELTING_TEMPERATURE,AVG_DROP_TMP(:,:,:,LPC%ARRAY_INDEX ))
      ELSE
         AVG_DROP_TMP(:,:,:,LPC%ARRAY_INDEX ) = MAX(TMPM,AVG_DROP_TMP(:,:,:,LPC%ARRAY_INDEX ))
      ENDIF
      AVG_DROP_DEN(:,:,:,LPC%ARRAY_INDEX ) = RUN_AVG_FAC*AVG_DROP_DEN(:,:,:,LPC%ARRAY_INDEX ) + OMRAF*DROP_DEN
      AVG_DROP_AREA(:,:,:,LPC%ARRAY_INDEX) = RUN_AVG_FAC*AVG_DROP_AREA(:,:,:,LPC%ARRAY_INDEX) + OMRAF*DROP_AREA
      WHERE (AVG_DROP_DEN(:,:,:,LPC%ARRAY_INDEX )<0.0001_EB .AND. ABS(DROP_DEN)<TWO_EPSILON_EB) &
         AVG_DROP_DEN(:,:,:,LPC%ARRAY_INDEX ) = 0.0_EB
    
   ENDDO PART_CLASS_SUM_LOOP

   ! Get total particle density

   DO IP=1,NLP
      LP=>LAGRANGIAN_PARTICLE(IP)       
      LPC=>LAGRANGIAN_PARTICLE_CLASS(LP%CLASS_INDEX)
      IF (LPC%MASSLESS_TRACER .OR. LPC%MASSLESS_TARGET) CYCLE
      II = LP%ONE_D%IIG
      JJ = LP%ONE_D%JJG
      KK = LP%ONE_D%KKG
      RVC = RDX(II)*RRN(II)*RDY(JJ)*RDZ(KK)
      DROP_DEN_ALL(II,JJ,KK) = DROP_DEN(II,JJ,KK) + LP%PWT*LP%MASS*RVC
   ENDDO

   AVG_DROP_DEN_ALL(:,:,:) = RUN_AVG_FAC*AVG_DROP_DEN_ALL(:,:,:) + OMRAF*DROP_DEN_ALL
   
ENDIF SUM_PART_QUANTITIES

! Remove PARTICLEs that have completely evaporated

CALL REMOVE_PARTICLES(T,NM)

END SUBROUTINE PARTICLE_MASS_ENERGY_TRANSFER


SUBROUTINE PARTICLE_MOMENTUM_TRANSFER(NM)

! Add PARTICLE momentum as a force term in momentum equation

USE TRAN, ONLY : GET_IJK
INTEGER, INTENT(IN) :: NM
REAL(EB), POINTER, DIMENSION(:,:,:) :: FVXS=>NULL(),FVYS=>NULL(),FVZS=>NULL(), &
                                       FVXN=>NULL(),FVYN=>NULL(),FVZN=>NULL()
REAL(EB) :: XI,YJ,ZK,SUMX,SUMY,SUMZ,WIDTH,RVC,RK
INTEGER :: II,JJ,KK,IIX,JJY,KKZ,I,J,K,IC,IW,IOR,IIG,JJG,KKG,III,JJJ,KKK,I1,I2,J1,J2,K1,K2
TYPE (LAGRANGIAN_PARTICLE_TYPE), POINTER :: LP=>NULL()
TYPE (LAGRANGIAN_PARTICLE_CLASS_TYPE), POINTER :: LPC=>NULL()
TYPE(WALL_TYPE), POINTER :: WC=>NULL()

CALL POINT_TO_MESH(NM)

FVXS  => WORK1
FVYS  => WORK2
FVZS  => WORK3

FVXS  = 0._EB
FVYS  = 0._EB
FVZS  = 0._EB

ISO_SELECT: SELECT CASE(ISOTROPIC_PARTICLE_FORCE)

CASE DEFAULT

   SUM_MOMENTUM_LOOP: DO I=1,NLP
      LP=>LAGRANGIAN_PARTICLE(I)
      LPC=>LAGRANGIAN_PARTICLE_CLASS(LP%CLASS_INDEX)
      IF (LP%ONE_D%IOR/=0) CYCLE SUM_MOMENTUM_LOOP
      IF (LPC%MASSLESS_TRACER .OR. LPC%MASSLESS_TARGET)    CYCLE SUM_MOMENTUM_LOOP
      CALL GET_IJK(LP%X,LP%Y,LP%Z,NM,XI,YJ,ZK,II,JJ,KK)
      IF (SOLID(CELL_INDEX(II,JJ,KK))) CYCLE SUM_MOMENTUM_LOOP
      IIX = FLOOR(XI+.5_EB)
      JJY = FLOOR(YJ+.5_EB)
      KKZ = FLOOR(ZK+.5_EB)

      IC = CELL_INDEX(IIX,JJ,KK)
      IW = WALL_INDEX(IC,1)
      IF (WALL(IW)%BOUNDARY_TYPE/=SOLID_BOUNDARY) THEN
         FVXS(IIX,JJ,KK) = FVXS(IIX,JJ,KK) - LP%ACCEL_X
      ENDIF
      IC = CELL_INDEX(II,JJY,KK)
      IW = WALL_INDEX(IC,2) 
      IF (WALL(IW)%BOUNDARY_TYPE/=SOLID_BOUNDARY) THEN
         FVYS(II,JJY,KK) = FVYS(II,JJY,KK) - LP%ACCEL_Y
      ENDIF
      IC = CELL_INDEX(II,JJ,KKZ)
      IW = WALL_INDEX(IC,3) 
      IF (WALL(IW)%BOUNDARY_TYPE/=SOLID_BOUNDARY) THEN
         FVZS(II,JJ,KKZ) = FVZS(II,JJ,KKZ) - LP%ACCEL_Z
      ENDIF

      IF (LPC%PERIODIC_X) THEN
         IF (IIX==0)    FVXS(IBAR,JJ,KK) = FVXS(IBAR,JJ,KK) - LP%ACCEL_X
         IF (IIX==IBAR) FVXS(0,JJ,KK)    = FVXS(0,JJ,KK)    - LP%ACCEL_X
      ENDIF
      IF (LPC%PERIODIC_Y) THEN
         IF (JJY==0)    FVYS(II,JBAR,KK) = FVYS(II,JBAR,KK) - LP%ACCEL_Y
         IF (JJY==JBAR) FVYS(II,0,KK)    = FVYS(II,0,KK)    - LP%ACCEL_Y
      ENDIF
      IF (LPC%PERIODIC_Z) THEN
         IF (KKZ==0)    FVZS(II,JJ,KBAR) = FVZS(II,JJ,KBAR) - LP%ACCEL_Z
         IF (KKZ==KBAR) FVZS(II,JJ,0)    = FVZS(II,JJ,0)    - LP%ACCEL_Z
      ENDIF

   ENDDO SUM_MOMENTUM_LOOP

   DO K=0,KBAR
      DO J=0,JBAR
         DO I=0,IBAR
            FVX(I,J,K) = FVX(I,J,K) + FVXS(I,J,K)
            FVY(I,J,K) = FVY(I,J,K) + FVYS(I,J,K)
            FVZ(I,J,K) = FVZ(I,J,K) + FVZS(I,J,K)
         ENDDO
      ENDDO
   ENDDO

CASE(1)

   ! experimental
   SUM_MOMENTUM_LOOP_2: DO I=1,NLP
      LP=>LAGRANGIAN_PARTICLE(I)
      LPC=>LAGRANGIAN_PARTICLE_CLASS(LP%CLASS_INDEX)
      IF (LP%ONE_D%IOR/=0) CYCLE SUM_MOMENTUM_LOOP_2
      IF (LPC%MASSLESS_TRACER .OR. LPC%MASSLESS_TARGET)    CYCLE SUM_MOMENTUM_LOOP_2
      CALL GET_IJK(LP%X,LP%Y,LP%Z,NM,XI,YJ,ZK,II,JJ,KK)
      IF (SOLID(CELL_INDEX(II,JJ,KK))) CYCLE SUM_MOMENTUM_LOOP_2

      FVXS(II,JJ,KK) = FVXS(II,JJ,KK) - LP%ACCEL_X
      FVYS(II,JJ,KK) = FVYS(II,JJ,KK) - LP%ACCEL_Y
      FVZS(II,JJ,KK) = FVZS(II,JJ,KK) - LP%ACCEL_Z
   ENDDO SUM_MOMENTUM_LOOP_2

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
      
      SELECT CASE(WC%BOUNDARY_TYPE)
         CASE(SOLID_BOUNDARY)
            ! apply Dirichlet bc F=0
            SELECT CASE(ABS(IOR))
               CASE(1); FVXS(II,JJ,KK) = -FVXS(IIG,JJG,KKG)
               CASE(2); FVYS(II,JJ,KK) = -FVYS(IIG,JJG,KKG)
               CASE(3); FVZS(II,JJ,KK) = -FVZS(IIG,JJG,KKG)
            END SELECT
         CASE(OPEN_BOUNDARY,MIRROR_BOUNDARY)
            ! apply Neumann bc dF/dx=0
            SELECT CASE(ABS(IOR))
               CASE(1); FVXS(II,JJ,KK) = FVXS(IIG,JJG,KKG)
               CASE(2); FVYS(II,JJ,KK) = FVYS(IIG,JJG,KKG)
               CASE(3); FVZS(II,JJ,KK) = FVZS(IIG,JJG,KKG)
            END SELECT
         CASE(PERIODIC_BOUNDARY)
            SELECT CASE(ABS(IOR))
               CASE( 1); FVXS(0,JJG,KKG)    = FVXS(IBAR,JJG,KKG)
               CASE(-1); FVXS(IBP1,JJG,KKG) = FVXS(1,JJG,KKG)
               CASE( 2); FVYS(IIG,0,KKG)    = FVYS(IIG,JBAR,KKG)
               CASE(-2); FVYS(IIG,JBP1,KKG) = FVYS(IIG,1,KKG)
               CASE( 3); FVZS(IIG,JJG,0)    = FVZS(IIG,JJG,KBAR)
               CASE(-3); FVZS(IIG,JJG,KBP1) = FVZS(IIG,JJG,1)
            END SELECT
      END SELECT

   ENDDO WALL_LOOP

   DO K=0,KBAR
      DO J=0,JBAR
         DO I=0,IBAR
            FVX(I,J,K) = FVX(I,J,K) + 0.5_EB*(FVXS(I,J,K)+FVXS(I+1,J,K))
            FVY(I,J,K) = FVY(I,J,K) + 0.5_EB*(FVYS(I,J,K)+FVYS(I,J+1,K))
            FVZ(I,J,K) = FVZ(I,J,K) + 0.5_EB*(FVZS(I,J,K)+FVZS(I,J,K+1))
         ENDDO
      ENDDO
   ENDDO

CASE(2) ! Kernel smoothing. Experimental

   FVXN  => WORK4
   FVYN  => WORK5
   FVZN  => WORK6

   FVXN  = 0._EB
   FVYN  = 0._EB
   FVZN  = 0._EB

   SUMX=0._EB
   SUMY=0._EB
   SUMZ=0._EB
   SUM_MOMENTUM_LOOP2: DO I=1,NLP
      LP=>LAGRANGIAN_PARTICLE(I)
      LPC=>LAGRANGIAN_PARTICLE_CLASS(LP%CLASS_INDEX)
      IF (LPC%MASSLESS_TRACER .OR. LPC%MASSLESS_TARGET) CYCLE SUM_MOMENTUM_LOOP2
      IF (LP%ONE_D%IOR/=0) CYCLE SUM_MOMENTUM_LOOP2
      CALL GET_IJK(LP%X,LP%Y,LP%Z,NM,XI,YJ,ZK,II,JJ,KK)
      IF (SOLID(CELL_INDEX(II,JJ,KK))) CYCLE SUM_MOMENTUM_LOOP2
      IF (AVG_DROP_DEN_ALL(II,JJ,KK)<TWO_EPSILON_EB) CYCLE SUM_MOMENTUM_LOOP2
      ! LOOP over a small neighborhood of the particle [-1,+1]
      ! Assumming kernel has compact support and half_width DX
      SUMX=SUMX-LP%ACCEL_X
      SUMY=SUMY-LP%ACCEL_Y
      SUMZ=SUMZ-LP%ACCEL_Z
      IIX = FLOOR(XI+.5_EB)
      JJY = FLOOR(YJ+.5_EB)
      KKZ = FLOOR(ZK+.5_EB)
      I1=MAX(IIX-2,0)
      I2=MIN(IIX+2,IBAR)
      J1=MAX(JJ-2,0)
      J2=MIN(JJ+2,JBAR)
      K1=MAX(KK-2,0)
      K2=MIN(KK+2,KBAR)
      WIDTH=0._EB
      
      WIDTH=(DX(II)*DY(JJ)*DZ(KK))**(1._EB/3._EB)
      ! Interpolate force density at (upper) (staggered) cell corners 
      DO III=I1,I2
          DO JJJ=J1,J2
             DO KKK=K1,K2
                   RK=SQRT((XC(III+1)-LP%X)**2+(Y(JJJ)-LP%Y)**2+(Z(KKK)-LP%Z)**2)      
                   FVXN(III,JJJ,KKK) =  FVXN(III,JJJ,KKK) - KERNEL(RK,WIDTH) * LP%ACCEL_X
                   
            ENDDO
          ENDDO
      ENDDO

      I1=MAX(II-2,0)
      I2=MIN(II+2,IBAR)
      J1=MAX(JJY-2,0)
      J2=MIN(JJY+2,JBAR)
      K1=MAX(KK-2,0)
      K2=MIN(KK+2,KBAR)
         DO III=I1,I2
          DO JJJ=J1,J2
             DO KKK=K1,K2                   
                   RK=SQRT((X(III)-LP%X)**2+(YC(JJJ+1)-LP%Y)**2+(Z(KKK)-LP%Z)**2)
                   FVYN(III,JJJ,KKK) =  FVYN(III,JJJ,KKK) - KERNEL(RK,WIDTH) * LP%ACCEL_Y
             ENDDO
          ENDDO
      ENDDO

      I1=MAX(II-2,0)
      I2=MIN(II+2,IBAR)
      J1=MAX(JJ-2,0)
      J2=MIN(JJ+2,JBAR)
      K1=MAX(KKZ-2,0)
      K2=MIN(KKZ+2,KBAR)
      DO III=I1,I2
          DO JJJ=J1,J2
             DO KKK=K1,K2
                IC = CELL_INDEX(III,JJJ,KKK)   

                RK=SQRT((X(III)-LP%X)**2+(Y(JJJ)-LP%Y)**2+(ZC(KKK+1)-LP%Z)**2)
                FVZN(III,JJJ,KKK) =  FVZN(III,JJJ,KKK) - KERNEL(RK,WIDTH) * LP%ACCEL_Z
              ENDDO
           ENDDO
        ENDDO
   ENDDO SUM_MOMENTUM_LOOP2

   ! Approximate volume integrals

   DO III=0,IBAR
      DO JJJ=0,JBAR
         DO KKK=0,KBAR
            I1=MAX(III-1,0)
            I2=III
            J1=MAX(JJJ-1,0)
            J2=JJJ
            K1=MAX(KKK-1,0)
            K2=KKK
            RVC   = RDX(III)*RRN(III)*RDY(JJJ)*RDZ(KKK)
            FVXS(III,JJJ,KKK) = 1._EB/8._EB * (FVXN(I1,J1,K1) + FVXN(I1,J1,K2) + FVXN(I1,J2,K2) + &
                                              FVXN(I2,J2,K2) + FVXN(I2,J2,K1) + FVXN(I2,J1,K1) + &
                                              FVXN(I1,J2,K1) + FVXN(I2,J1,K2)) / RVC
            FVYS(III,JJJ,KKK) = 1._EB/8._EB * (FVYN(I1,J1,K1) + FVYN(I1,J1,K2) + FVYN(I1,J2,K2) + &
                                              FVYN(I2,J2,K2) + FVYN(I2,J2,K1) + FVYN(I2,J1,K1) + &
                                              FVYN(I1,J2,K1) + FVYN(I2,J1,K2)) / RVC
            FVZS(III,JJJ,KKK) = 1._EB/8._EB * (FVZN(I1,J1,K1) + FVZN(I1,J1,K2) + FVZN(I1,J2,K2) + &
                                              FVZN(I2,J2,K2) + FVZN(I2,J2,K1) + FVZN(I2,J1,K1) + &
                                              FVZN(I1,J2,K1) + FVZN(I2,J1,K2)) / RVC
         ENDDO
      ENDDO
   ENDDO
   IF (ABS(SUMX)>TWO_EPSILON_EB) FVXS=FVXS*SUMX/SUM(FVXS)
   IF (ABS(SUMY)>TWO_EPSILON_EB) FVYS=FVYS*SUMY/SUM(FVYS)
   IF (ABS(SUMZ)>TWO_EPSILON_EB) FVZS=FVZS*SUMZ/SUM(FVZS)
   WALL_LOOP1: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
      WC=>WALL(IW)
      IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE WALL_LOOP1
      II  = WC%ONE_D%II
      JJ  = WC%ONE_D%JJ
      KK  = WC%ONE_D%KK
      IOR = WC%ONE_D%IOR
      IIG = WC%ONE_D%IIG
      JJG = WC%ONE_D%JJG
      KKG = WC%ONE_D%KKG
      
      SELECT CASE(WC%BOUNDARY_TYPE)
         CASE(SOLID_BOUNDARY)
            ! apply Dirichlet bc F=0
            SELECT CASE(ABS(IOR))
               CASE(1); FVXS(II,JJ,KK) = -FVXS(IIG,JJG,KKG)
               CASE(2); FVYS(II,JJ,KK) = -FVYS(IIG,JJG,KKG)
               CASE(3); FVZS(II,JJ,KK) = -FVZS(IIG,JJG,KKG)
            END SELECT
         CASE(OPEN_BOUNDARY,MIRROR_BOUNDARY)
            ! apply Neumann bc dF/dx=0
            SELECT CASE(ABS(IOR))
               CASE(1); FVXS(II,JJ,KK) = FVXS(IIG,JJG,KKG)
               CASE(2); FVYS(II,JJ,KK) = FVYS(IIG,JJG,KKG)
               CASE(3); FVZS(II,JJ,KK) = FVZS(IIG,JJG,KKG)
            END SELECT
         CASE(PERIODIC_BOUNDARY)
            SELECT CASE(ABS(IOR))
               CASE( 1); FVXS(0,JJG,KKG)    = FVXS(IBAR,JJG,KKG)
               CASE(-1); FVXS(IBP1,JJG,KKG) = FVXS(1,JJG,KKG)
               CASE( 2); FVYS(IIG,0,KKG)    = FVYS(IIG,JBAR,KKG)
               CASE(-2); FVYS(IIG,JBP1,KKG) = FVYS(IIG,1,KKG)
               CASE( 3); FVZS(IIG,JJG,0)    = FVZS(IIG,JJG,KBAR)
               CASE(-3); FVZS(IIG,JJG,KBP1) = FVZS(IIG,JJG,1)
            END SELECT
      END SELECT
   ENDDO WALL_LOOP1

   DO K=0,KBAR
      DO J=0,JBAR
         DO I=0,IBAR
            FVX(I,J,K) = FVX(I,J,K) + FVXS(I,J,K)
            FVY(I,J,K) = FVY(I,J,K) + FVYS(I,J,K)
            FVZ(I,J,K) = FVZ(I,J,K) + FVZS(I,J,K)
         ENDDO
      ENDDO
   ENDDO

END SELECT ISO_SELECT

CONTAINS

REAL(EB) FUNCTION KERNEL(R,DELTA)
   REAL(EB),INTENT(IN)::R,DELTA
   IF(R/DELTA>2.0_EB) THEN
      KERNEL=0._EB
   ELSE
      KERNEL = 1._EB/(PI * SQRT(PI) * DELTA**3) * EXP(-(R/DELTA)**2)
   ENDIF
END FUNCTION KERNEL

END SUBROUTINE PARTICLE_MOMENTUM_TRANSFER
 
 
SUBROUTINE REMOVE_PARTICLES(T,NM)

! Remove PARTICLEs that have left the current mesh (NM) or are no longer to be tracked
 
INTEGER, INTENT(IN) :: NM
INTEGER :: IKILL,IP
REAL(EB), INTENT(IN) :: T
TYPE (LAGRANGIAN_PARTICLE_TYPE), POINTER :: LP
TYPE (LAGRANGIAN_PARTICLE_CLASS_TYPE), POINTER :: LPC
TYPE (SURFACE_TYPE), POINTER :: SF

IKILL = 0

PARTICLE_LOOP: DO IP=1,NLP

   WEED_LOOP: DO

      IF (IP>NLP-IKILL) EXIT PARTICLE_LOOP

      LP  => MESHES(NM)%LAGRANGIAN_PARTICLE(IP)  
      LPC => LAGRANGIAN_PARTICLE_CLASS(LP%CLASS_INDEX)
      SF  => SURFACE(LPC%SURF_INDEX)
  
      ! Remove particles that are too small

      IF (LPC%SOLID_PARTICLE .AND. SF%THERMALLY_THICK .AND. LP%ONE_D%BURNAWAY) THEN
         CALL PARTICLE_ORPHANAGE
         CYCLE WEED_LOOP
      ENDIF

      IF (LPC%LIQUID_DROPLET .AND. LP%ONE_D%X(1)<=LPC%KILL_RADIUS) THEN
         CALL PARTICLE_ORPHANAGE
         CYCLE WEED_LOOP
      ENDIF      
   
      ! Remove particles that are too old

      IF (T-LP%T_INSERT > LPC%LIFETIME) THEN
         CALL PARTICLE_ORPHANAGE
         CYCLE WEED_LOOP
      ENDIF
   
      ! Remove particles that have left the active mesh

      IF (LP%X>MESHES(NM)%XS .AND. LP%X<MESHES(NM)%XF .AND. LP%Y>MESHES(NM)%YS .AND. LP%Y<MESHES(NM)%YF .AND. &
          LP%Z>MESHES(NM)%ZS .AND. LP%Z<MESHES(NM)%ZF) CYCLE PARTICLE_LOOP

      ! Replace all other particles

      CALL PARTICLE_ORPHANAGE

   ENDDO WEED_LOOP

ENDDO PARTICLE_LOOP
 
NLP = NLP - IKILL

CONTAINS
 
SUBROUTINE PARTICLE_ORPHANAGE

! Determine if the given particle is now in another mesh, and if so, assign it to the ORPHAN array.
 
USE MEMORY_FUNCTIONS, ONLY: REALLOCATE_STORAGE_ARRAYS
INTEGER :: OM,NOM,TAG,N_NEW_STORAGE_SLOTS
TYPE (MESH_TYPE), POINTER :: M=>NULL()
TYPE (OMESH_TYPE), POINTER :: M2=>NULL()

! Check to see if the particle has entered another mesh (NOM)

NOM = 0
SEARCH_LOOP: DO OM=1,NMESHES
   IF (MESHES(NM)%OMESH(OM)%NIC_R==0) CYCLE SEARCH_LOOP
   IF (EVACUATION_ONLY(OM)) CYCLE SEARCH_LOOP
   M=>MESHES(OM)
   IF (LP%X>M%XS .AND. LP%X<M%XF .AND.  LP%Y>M%YS .AND. LP%Y<M%YF .AND.  LP%Z>M%ZS .AND. LP%Z<M%ZF) THEN
      NOM = OM
      EXIT SEARCH_LOOP
   ENDIF
ENDDO SEARCH_LOOP

! Don't remove particles that just cycle in a periodic domain

IF (LPC%PERIODIC_X) THEN
   IF (LP%X>=XF_MAX) THEN
      LP%X = LP%X - (XF_MAX- XS_MIN)
      RETURN
   ENDIF
   IF (LP%X<=XS_MIN) THEN
      LP%X = LP%X + (XF_MAX- XS_MIN)
      RETURN
   ENDIF
ENDIF

IF (LPC%PERIODIC_Y) THEN
   IF (LP%Y>=YF_MAX) THEN
      LP%Y = LP%Y - (YF_MAX- YS_MIN)
      RETURN
   ENDIF
   IF (LP%Y<=YS_MIN) THEN
      LP%Y = LP%Y + (YF_MAX- YS_MIN)
      RETURN
   ENDIF
ENDIF

IF (LPC%PERIODIC_Z) THEN
   IF (LP%Z>=ZF_MAX) THEN
      LP%Z = LP%Z - (ZF_MAX- ZS_MIN)
      RETURN
   ENDIF
   IF (LP%Z<=ZS_MIN) THEN
      LP%Z = LP%Z + (ZF_MAX- ZS_MIN)
      RETURN
   ENDIF
ENDIF

! For particles that have entered another mesh (NOM), add them to the "orphanage"

IF (NOM/=0) THEN
   M2=>MESHES(NM)%OMESH(NOM)
   M2%N_PART_ORPHANS(LP%CLASS_INDEX) = M2%N_PART_ORPHANS(LP%CLASS_INDEX) + 1

   IF (M2%N_PART_ORPHANS(LP%CLASS_INDEX)>M2%ORPHAN_PARTICLE_STORAGE(LP%CLASS_INDEX)%N_STORAGE_SLOTS) THEN
      N_NEW_STORAGE_SLOTS = M2%N_PART_ORPHANS(LP%CLASS_INDEX) - M2%ORPHAN_PARTICLE_STORAGE(LP%CLASS_INDEX)%N_STORAGE_SLOTS
      CALL REALLOCATE_STORAGE_ARRAYS(NM,2,LP%CLASS_INDEX,N_NEW_STORAGE_SLOTS,NOM)
   ENDIF

   M2%ORPHAN_PARTICLE_STORAGE(LP%CLASS_INDEX)%REALS(:,M2%N_PART_ORPHANS(LP%CLASS_INDEX)) = &
      MESHES(NM)%PARTICLE_STORAGE(LP%CLASS_INDEX)%REALS(:,LP%STORAGE_INDEX)
   M2%ORPHAN_PARTICLE_STORAGE(LP%CLASS_INDEX)%INTEGERS(:,M2%N_PART_ORPHANS(LP%CLASS_INDEX)) = &
      MESHES(NM)%PARTICLE_STORAGE(LP%CLASS_INDEX)%INTEGERS(:,LP%STORAGE_INDEX)
   M2%ORPHAN_PARTICLE_STORAGE(LP%CLASS_INDEX)%LOGICALS(:,M2%N_PART_ORPHANS(LP%CLASS_INDEX)) = &
      MESHES(NM)%PARTICLE_STORAGE(LP%CLASS_INDEX)%LOGICALS(:,LP%STORAGE_INDEX)
ENDIF

! Zero out storage for particle that is being removed

MESHES(NM)%PARTICLE_STORAGE(LP%CLASS_INDEX)%NEXT_AVAILABLE_SLOT = &
   MIN(MESHES(NM)%PARTICLE_STORAGE(LP%CLASS_INDEX)%NEXT_AVAILABLE_SLOT,LP%STORAGE_INDEX)
MESHES(NM)%PARTICLE_STORAGE(LP%CLASS_INDEX)%REALS(:,LP%STORAGE_INDEX)    = 0._EB
MESHES(NM)%PARTICLE_STORAGE(LP%CLASS_INDEX)%LOGICALS(:,LP%STORAGE_INDEX) = .FALSE.
MESHES(NM)%PARTICLE_STORAGE(LP%CLASS_INDEX)%INTEGERS(:,LP%STORAGE_INDEX) = 0

! Move particle at the end of the line into the vacated particle

IF (IP<NLP-IKILL) THEN
   MESHES(NM)%LAGRANGIAN_PARTICLE(IP) = MESHES(NM)%LAGRANGIAN_PARTICLE(NLP-IKILL)
   MESHES(NM)%LAGRANGIAN_PARTICLE(IP)%ARRAY_INDEX = IP
ENDIF

IKILL = IKILL + 1

END SUBROUTINE PARTICLE_ORPHANAGE
 
END SUBROUTINE REMOVE_PARTICLES

 
SUBROUTINE REMOVE_OLDEST_PARTICLE(NM,LPC_INDEX,NLP,LP_INDEX)

! Remove the oldest particle of class LPC_INDEX and move particle NLP into its place.

INTEGER, INTENT(IN) :: NM,LPC_INDEX,NLP
INTEGER, INTENT(OUT) :: LP_INDEX
INTEGER :: STORAGE_INDEX,I,TAG,TAG_MIN

! Look for the oldest particle of this class

TAG_MIN = HUGE(NLP)

DO I=1,MESHES(NM)%PARTICLE_STORAGE(LPC_INDEX)%N_STORAGE_SLOTS
   TAG = MESHES(NM)%PARTICLE_STORAGE(LPC_INDEX)%INTEGERS(1,I)
   IF (TAG>0 .AND. TAG<TAG_MIN) THEN
       TAG_MIN = TAG
       STORAGE_INDEX = I
       LP_INDEX = MESHES(NM)%PARTICLE_STORAGE(LPC_INDEX)%INTEGERS(2,I)
   ENDIF
ENDDO

! Zero out the storage for the oldest particle

MESHES(NM)%PARTICLE_STORAGE(LPC_INDEX)%NEXT_AVAILABLE_SLOT = &
   MIN(MESHES(NM)%PARTICLE_STORAGE(LPC_INDEX)%NEXT_AVAILABLE_SLOT,STORAGE_INDEX)
MESHES(NM)%PARTICLE_STORAGE(LPC_INDEX)%REALS(:,STORAGE_INDEX)    = 0._EB
MESHES(NM)%PARTICLE_STORAGE(LPC_INDEX)%LOGICALS(:,STORAGE_INDEX) = .FALSE.
MESHES(NM)%PARTICLE_STORAGE(LPC_INDEX)%INTEGERS(:,STORAGE_INDEX) = 0

! Move the last LAGRANGIAN_PARTICLE into the slot left open after removing the oldest LAGRANGIAN_PARTICLE

MESHES(NM)%LAGRANGIAN_PARTICLE(LP_INDEX) = MESHES(NM)%LAGRANGIAN_PARTICLE(NLP)
MESHES(NM)%LAGRANGIAN_PARTICLE(LP_INDEX)%ARRAY_INDEX = LP_INDEX

END SUBROUTINE REMOVE_OLDEST_PARTICLE


SUBROUTINE GET_REV_part(MODULE_REV,MODULE_DATE)
INTEGER,INTENT(INOUT) :: MODULE_REV
CHARACTER(255),INTENT(INOUT) :: MODULE_DATE

WRITE(MODULE_DATE,'(A)') partrev(INDEX(partrev,':')+2:LEN_TRIM(partrev)-2)
READ (MODULE_DATE,'(I5)') MODULE_REV
WRITE(MODULE_DATE,'(A)') partdate

END SUBROUTINE GET_REV_part
END MODULE PART
