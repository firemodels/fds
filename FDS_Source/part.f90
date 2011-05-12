MODULE PART
 
USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_POINTERS
USE TRAN
IMPLICIT NONE
PRIVATE
PUBLIC INSERT_DROPLETS_AND_PARTICLES,UPDATE_PARTICLES,REMOVE_DROPLETS,INITIALIZE_DROPLETS,GET_REV_part
CHARACTER(255), PARAMETER :: partid='$Id$'
CHARACTER(255), PARAMETER :: partrev='$Revision$'
CHARACTER(255), PARAMETER :: partdate='$Date$'
REAL(EB) :: INSERT_RATE = 0._EB
INTEGER :: INSERT_COUNT = 0
 
CONTAINS
 

SUBROUTINE INITIALIZE_DROPLETS(NM)

! Insert droplets into the domain at the start of calculation

USE COMP_FUNCTIONS, ONLY : SECOND 
USE PHYSICAL_FUNCTIONS, ONLY : DROPLET_SIZE_DISTRIBUTION 
USE MEMORY_FUNCTIONS, ONLY : RE_ALLOCATE_DROPLETS
 
REAL(EB) :: LL,UL,BIN_SIZE,TNOW
INTEGER  :: I,J,IL,IU,IPC
INTEGER, INTENT(IN) :: NM
TYPE (PARTICLE_CLASS_TYPE), POINTER :: PC=>NULL()

IF (N_PART==0) RETURN            ! Don't waste time if no particles
IF (EVACUATION_ONLY(NM)) RETURN  ! Don't waste time if an evac mesh
 
TNOW=SECOND()
CALL POINT_TO_MESH(NM)

PART_CLASS_LOOP: DO IPC=1,N_PART

   PC=>PARTICLE_CLASS(IPC)
   
   ! If particles/droplets have a size distribution, initialize here
 
   IF_SIZE_DISTRIBUTION: IF (.NOT.PC%MONODISPERSE .AND. (PC%DIAMETER > 0._EB)) THEN
      CALL DROPLET_SIZE_DISTRIBUTION(PC%DIAMETER,PC%R_CDF(:),PC%CDF(:),NDC,PC%GAMMA,PC%SIGMA)
      BIN_SIZE = PC%R_CDF(NDC)/REAL(NSTRATA,EB)
      STRATIFY: DO I=1,NSTRATA
         LL = (I-1)*BIN_SIZE
         UL =  I   *BIN_SIZE
         LL_LOOP: DO J=1,NDC
            IF (PC%R_CDF(J)>LL) THEN
               IL = J-1 
               PC%IL_CDF(I) = J-1
               EXIT LL_LOOP
            ENDIF
         ENDDO LL_LOOP
         UL_LOOP: DO J=NDC,1,-1
            IF (PC%R_CDF(J)<=UL) THEN
               IU = J 
               PC%IU_CDF(I) = J
               EXIT UL_LOOP
            ENDIF
         ENDDO UL_LOOP
         PC%W_CDF(I) = PC%CDF(IU) - PC%CDF(IL)
      ENDDO STRATIFY
   ENDIF IF_SIZE_DISTRIBUTION

   ! If pacticles/droplets can break up, compute normalized (median = 1) size distribution for child droplets

   IF (PC%BREAKUP .AND. .NOT.PC%MONODISPERSE) THEN
      CALL DROPLET_SIZE_DISTRIBUTION(1._EB,PC%CHILD_R_CDF(:),PC%CHILD_CDF(:),NDC, &
                                     PC%BREAKUP_CHILD_GAMMA,PC%BREAKUP_CHILD_SIGMA)
   ENDIF

ENDDO PART_CLASS_LOOP

TUSED(8,NM)=TUSED(8,NM)+SECOND()-TNOW
END SUBROUTINE INITIALIZE_DROPLETS
 


SUBROUTINE INSERT_DROPLETS_AND_PARTICLES(T,NM)

! Insert sprinkler droplets and lagrangian particles into the domain every time step

USE COMP_FUNCTIONS, ONLY : SECOND 
USE MEMORY_FUNCTIONS, ONLY : RE_ALLOCATE_DROPLETS
USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
USE GEOMETRY_FUNCTIONS, ONLY: RANDOM_RECTANGLE,RANDOM_CONE
USE DEVICE_VARIABLES
USE CONTROL_VARIABLES
REAL(EB), INTENT(IN) :: T
INTEGER, INTENT(IN) :: NM
REAL(EB) :: PHI_RN,RN,FLOW_RATE,THETA_RN,SPHI,CPHI,MASS_SUM,D_PRES_FACTOR, &
            STHETA,CTHETA,PWT0,DROPLET_SPEED,XI,YJ,ZK,SHIFT1,SHIFT2,XTMP,YTMP,ZTMP,VLEN, &
            TRIGT1,TRIGT2,TNOW,TSI,PIPE_PRESSURE,MASS_PER_TIME,MASS_PER_VOLUME,X1,X2,Y1,Y2,Z1,Z2,BLOCK_VOLUME, &
            ETA,ETA_MAX,ETA_MIN
REAL(EB), PARAMETER :: VENT_OFFSET=0.1
INTEGER :: I,KS,II,JJ,KK,IC,IL,IU,IPC,DROP_SUM,IIG,JJG,KKG,IW,IBC,IOR,STRATUM,IB
INTEGER :: N_INITIAL,N,N_INSERT,ILAT
LOGICAL :: INSERT_ANOTHER_BATCH
TYPE (PROPERTY_TYPE), POINTER :: PY=>NULL()
TYPE (TABLES_TYPE), POINTER :: TA=>NULL()
TYPE (DEVICE_TYPE), POINTER :: DV=>NULL()
TYPE (SURFACE_TYPE), POINTER :: SF=>NULL()
TYPE (DROPLET_TYPE), POINTER :: DR=>NULL()
TYPE (PARTICLE_CLASS_TYPE), POINTER :: PC=>NULL()
TYPE (INITIALIZATION_TYPE), POINTER :: IN=>NULL()
 
IF (EVACUATION_ONLY(NM)) RETURN  ! Don't waste time if an evac mesh
IF (N_PART==0) RETURN            ! Don't waste time if no particles

TNOW=SECOND()
CALL POINT_TO_MESH(NM)

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

! Count the number of droplets/particles inserted

INSERT_RATE = 0._EB
INSERT_COUNT = 0

! Add more than one batch of particles/droplets if FDS time step is large enough

OVERALL_INSERT_LOOP: DO  

   INSERT_ANOTHER_BATCH = .FALSE.

   ! Loop over all devices, but look for sprinklers or nozzles
      
   SPRINKLER_INSERT_LOOP: DO KS=1,N_DEVC  
      DV => DEVICE(KS)
      PY => PROPERTY(DV%PROP_INDEX)
      IF (PY%PART_ID == 'null')   CYCLE SPRINKLER_INSERT_LOOP
      IF (DV%MESH/=NM)            CYCLE SPRINKLER_INSERT_LOOP
      IF (.NOT. DV%CURRENT_STATE) CYCLE SPRINKLER_INSERT_LOOP
      IPC = PY%PART_INDEX
      PC=>PARTICLE_CLASS(IPC)

      IF (ABS(DV%T_CHANGE-T)<=ZERO_P) THEN 
         DV%T = T
         CYCLE SPRINKLER_INSERT_LOOP
      ENDIF
!      IF (T < PY%PARTICLE_INSERT_CLOCK(NM)) CYCLE SPRINKLER_INSERT_LOOP
   
      INSERT_RATE = INSERT_RATE + PY%DROPLETS_PER_SECOND
      N_INSERT = NINT(REAL(PY%DROPLETS_PER_SECOND,EB)*DT)
      INSERT_COUNT = INSERT_COUNT + N_INSERT

      ! Determine sprinkler/nozzle flow rate
   
      IF (ABS(T_BEGIN-DV%T_CHANGE)<=ZERO_P .AND. PY%FLOW_RAMP_INDEX>=1) THEN
         TSI = T
      ELSE
         TSI = T - DV%T_CHANGE
      ENDIF
   
      IF (PY%PRESSURE_RAMP_INDEX>0) THEN
         PIPE_PRESSURE = EVALUATE_RAMP(REAL(N_OPEN_NOZZLES,EB),0._EB,PY%PRESSURE_RAMP_INDEX)
         D_PRES_FACTOR = (PY%OPERATING_PRESSURE/PIPE_PRESSURE)**(1._EB/3._EB)
         FLOW_RATE = PY%K_FACTOR*SQRT(PIPE_PRESSURE)
      ELSE
         PIPE_PRESSURE = PY%OPERATING_PRESSURE
         D_PRES_FACTOR = 1.0_EB
         FLOW_RATE = PY%FLOW_RATE
      ENDIF      
   
      FLOW_RATE = EVALUATE_RAMP(TSI,PY%FLOW_TAU,PY%FLOW_RAMP_INDEX)*FLOW_RATE*(PC%DENSITY/1000._EB)/60._EB  ! kg/s

      IF (FLOW_RATE <= 0._EB) THEN
         DV%T = T
         CYCLE SPRINKLER_INSERT_LOOP
      ENDIF
   
      ! Direction initialization stuff
   
      TRIGT1 = ACOS(-DV%ORIENTATION(3))
      IF (ABS(DV%ORIENTATION(2))<=ZERO_P) THEN
         TRIGT2 = ACOS(1._EB)
      ELSE
         TRIGT2 = ACOS(ABS(DV%ORIENTATION(1))/SQRT(DV%ORIENTATION(1)**2+DV%ORIENTATION(2)**2))
      ENDIF
    
      ! Droplet insertion loop
    
      MASS_SUM = 0._EB
      DROP_SUM = 0
   
      DROPLET_INSERT_LOOP: DO I=1,N_INSERT
         IF (NLP+1>MAXIMUM_DROPLETS) EXIT DROPLET_INSERT_LOOP
         NLP = NLP+1
         PARTICLE_TAG = PARTICLE_TAG + NMESHES
         IF (NLP>NLPDIM) THEN
            CALL RE_ALLOCATE_DROPLETS(1,NM,0,1000)
            DROPLET=>MESHES(NM)%DROPLET
         ENDIF
         DR=>DROPLET(NLP)
    
         ! Set droplet properties
    
         DR%TMP    = PC%TMP_INITIAL         ! Initial temperature
         DR%T      = T                      ! Time of insertion
         DR%IOR    = 0                      ! Orientation index (0 means it is not attached to a solid)
         DR%A_X    = 0._EB                  ! x component of drag on the gas in units of m/s2
         DR%A_Y    = 0._EB                  ! y component of drag on the gas in units of m/s2
         DR%A_Z    = 0._EB                  ! z component of drag on the gas in units of m/s2
         DR%CLASS = IPC                     ! The class of particles
         DR%TAG    = PARTICLE_TAG           ! A unique identifying integer
         IF (MOD(NLP,PC%SAMPLING)==0) THEN  ! Decide whether to show the droplet in Smokeview
            DR%SHOW = .TRUE.    
         ELSE
            DR%SHOW = .FALSE.
         ENDIF
         DR%SPLAT   = .FALSE.
         DR%WALL_INDEX = 0
   
         ! Randomly choose theta and phi
    
         CHOOSE_COORDS: DO
            PICK_PATTERN: IF(PY%SPRAY_PATTERN_INDEX>0) THEN !Use spray pattern table
               TA => TABLES(PY%SPRAY_PATTERN_INDEX)
               CALL RANDOM_NUMBER(RN)
               FIND_ROW: DO II=1,TA%NUMBER_ROWS
                  IF (RN>PY%TABLE_ROW(II)) CYCLE FIND_ROW
                  EXIT FIND_ROW
               END DO FIND_ROW
               CALL RANDOM_NUMBER(RN)
               !THETA_RN = TA%TABLE_DATA(II,1) + RN*(TA%TABLE_DATA(II,2)-TA%TABLE_DATA(II,1))
               ETA_MAX=0.5_EB*(COS(TA%TABLE_DATA(II,1))+1._EB)
               ETA_MIN=0.5_EB*(COS(TA%TABLE_DATA(II,2))+1._EB)
               ETA=ETA_MIN+(ETA_MAX-ETA_MIN)*RN
               THETA_RN=ACOS(2._EB*ETA-1._EB)
               CALL RANDOM_NUMBER(RN)
               PHI_RN = TA%TABLE_DATA(II,3) + RN*(TA%TABLE_DATA(II,4)-TA%TABLE_DATA(II,3))
             
               IF (PY%PRESSURE_RAMP_INDEX>0) THEN
                  DROPLET_SPEED = PY%V_FACTOR(II)*SQRT(PIPE_PRESSURE)
               ELSE 
                  DROPLET_SPEED = TA%TABLE_DATA(II,5)
               ENDIF
            ELSE PICK_PATTERN !Use conical spray
               !CALL RANDOM_NUMBER(RN)
               !THETA_RN = PY%SPRAY_ANGLE(1) + RN*(PY%SPRAY_ANGLE(2)-PY%SPRAY_ANGLE(1))
               CALL RANDOM_CHOICE(PY%SPRAY_LON_CDF(:),PY%SPRAY_LON,NDC2,PHI_RN)
               ILAT=MINLOC(ABS(PY%SPRAY_LON-PHI_RN),1) 
               CALL RANDOM_CHOICE(PY%SPRAY_LAT_CDF(:,ILAT),PY%SPRAY_LAT,NDC2,THETA_RN)
               !CALL RANDOM_NUMBER(RN)
               !PHI_RN = RN*TWOPI
               IF (PY%PRESSURE_RAMP_INDEX>0) THEN
                  DROPLET_SPEED = PY%V_FACTOR(1)*SQRT(PIPE_PRESSURE)
               ELSE
                  DROPLET_SPEED = PY%DROPLET_VELOCITY
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
            SELECT CASE(INT(SIGN(1._EB,DV%ORIENTATION(1))))
               CASE (-1)
                  IF (DV%ORIENTATION(2)>0) SHIFT2 = TWOPI - SHIFT2
               CASE ( 1)
               SELECT CASE(INT(SIGN(1._EB,DV%ORIENTATION(2))))
                  CASE (-1) 
                     SHIFT2 = PI-SHIFT2
                  CASE ( 1)
                     SHIFT2 = SHIFT2+ PI
               END SELECT
            END SELECT
            SHIFT1=SHIFT1+SHIFT2
            XTMP = VLEN * COS(SHIFT1)
            YTMP = VLEN * SIN(SHIFT1)
    
            ! Compute initial position and velocity of droplets
    
            DR%U = DROPLET_SPEED*XTMP
            DR%V = DROPLET_SPEED*YTMP
            DR%W = DROPLET_SPEED*ZTMP
            DR%X = DV%X + PY%OFFSET*XTMP
            DR%Y = DV%Y + PY%OFFSET*YTMP
            DR%Z = DV%Z + PY%OFFSET*ZTMP
            IF (TWO_D) THEN
               DR%V = 0._EB
               DR%Y = DV%Y
            ENDIF
            IF (DR%X<=XS .OR. DR%X>=XF) CYCLE CHOOSE_COORDS
            IF (DR%Y<=YS .OR. DR%Y>=YF) CYCLE CHOOSE_COORDS
            IF (DR%Z<=ZS .OR. DR%Z>=ZF) CYCLE CHOOSE_COORDS
            XI = CELLSI(FLOOR((DR%X-XS)*RDXINT))
            YJ = CELLSJ(FLOOR((DR%Y-YS)*RDYINT))
            ZK = CELLSK(FLOOR((DR%Z-ZS)*RDZINT))
            II = FLOOR(XI+1._EB)
            JJ = FLOOR(YJ+1._EB)
            KK = FLOOR(ZK+1._EB)
            IC = CELL_INDEX(II,JJ,KK)
            IF (.NOT.SOLID(IC)) EXIT CHOOSE_COORDS
    
         ENDDO CHOOSE_COORDS
   
         ! Randomly choose droplet size according to Cumulative Distribution Function (CDF)
    
         IF (PC%MONODISPERSE) THEN
            DR%R   = 0.5_EB*PC%DIAMETER*D_PRES_FACTOR
            DR%PWT = 1._EB
         ELSE
            CALL RANDOM_NUMBER(RN)            
            STRATUM = NINT(REAL(NSTRATA,EB)*RN+0.5_EB)
            IL = PC%IL_CDF(STRATUM)
            IU = PC%IU_CDF(STRATUM)
            CALL RANDOM_CHOICE(PC%CDF(IL:IU), PC%R_CDF(IL:IU), IU-IL,DR%R)
            DR%PWT = PC%W_CDF(STRATUM)
            DR%R = D_PRES_FACTOR*DR%R
            IF (2._EB*DR%R > PC%MAXIMUM_DIAMETER) THEN
               DR%PWT = DR%PWT*DR%R**3/(0.5_EB*PC%MAXIMUM_DIAMETER)**3
               DR%R = 0.5_EB*PC%MAXIMUM_DIAMETER
            ENDIF
         ENDIF
    
         ! Sum up mass of liquid being introduced
   
         MASS_SUM = MASS_SUM + DR%PWT*PC%FTPR*DR%R**3
         DROP_SUM = DROP_SUM + 1
         
      ENDDO DROPLET_INSERT_LOOP
    
      ! Compute weighting factor for the droplets just inserted
      
      IF (DROP_SUM > 0) THEN
         PWT0 = FLOW_RATE*DT/MASS_SUM
         DROPLET(NLP-DROP_SUM+1:NLP)%PWT = DROPLET(NLP-DROP_SUM+1:NLP)%PWT*PWT0
      ENDIF
      
      ! Indicate that droplets from this device have been inserted at this time T
   
      DV%T = T 
   
   ENDDO SPRINKLER_INSERT_LOOP
   
   ! Loop through all boundary cells and insert particles if appropriate
    
   WALL_INSERT_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
    
      IF (T < TW(IW))                        CYCLE WALL_INSERT_LOOP
      IF (BOUNDARY_TYPE(IW)/=SOLID_BOUNDARY) CYCLE WALL_INSERT_LOOP
      IBC =  IJKW(5,IW)
      SF  => SURFACE(IBC)
      IPC =  SF%PART_INDEX
      IF (IPC < 1)    CYCLE WALL_INSERT_LOOP
      PC  => PARTICLE_CLASS(IPC)
      IF (PC%DEVC_INDEX>0) THEN
         IF (.NOT.DEVICE(PC%DEVC_INDEX)%CURRENT_STATE) CYCLE WALL_INSERT_LOOP
      ENDIF
      IF (PC%CTRL_INDEX>0) THEN
         IF (.NOT.CONTROL(PC%CTRL_INDEX)%CURRENT_STATE) CYCLE WALL_INSERT_LOOP
      ENDIF
      IF (T < SF%PARTICLE_INSERT_CLOCK(NM)) CYCLE WALL_INSERT_LOOP
      INSERT_RATE = INSERT_RATE + NPPCW(IW)/SF%DT_INSERT
      INSERT_COUNT = INSERT_COUNT + NPPCW(IW)
      IF (UW(IW) >= -0.0001_EB) CYCLE WALL_INSERT_LOOP
    
      II = IJKW(1,IW)
      JJ = IJKW(2,IW)
      KK = IJKW(3,IW)
      IC = CELL_INDEX(II,JJ,KK)
      IF (.NOT.SOLID(IC)) CYCLE WALL_INSERT_LOOP
    
      IF (NM > 1) THEN
         IIG = IJKW(6,IW)
         JJG = IJKW(7,IW)
         KKG = IJKW(8,IW)
         IF (INTERPOLATED_MESH(IIG,JJG,KKG) > 0) CYCLE WALL_INSERT_LOOP
      ENDIF
    
      ! Loop over all particles for the IW-th cell
   
      IOR = IJKW(4,IW)
      MASS_SUM = 0._EB
      PARTICLE_INSERT_LOOP: DO I=1,NPPCW(IW) 
    
         IF (NLP+1 > MAXIMUM_DROPLETS) EXIT PARTICLE_INSERT_LOOP
    
         NLP = NLP+1
         PARTICLE_TAG = PARTICLE_TAG + NMESHES
         IF (NLP > NLPDIM) THEN
            CALL RE_ALLOCATE_DROPLETS(1,NM,0,1000)
            DROPLET => MESHES(NM)%DROPLET
         ENDIF
         DR=>DROPLET(NLP)
    
         IF (ABS(IOR)==1) THEN
            IF (IOR== 1) DR%X = X(II)   + VENT_OFFSET*DX(II+1)
            IF (IOR==-1) DR%X = X(II-1) - VENT_OFFSET*DX(II-1)
            CALL RANDOM_NUMBER(RN)
            DR%Y = Y(JJ-1) + DY(JJ)*RN
            CALL RANDOM_NUMBER(RN)
            DR%Z = Z(KK-1) + DZ(KK)*RN
         ENDIF
         IF (ABS(IOR)==2) THEN
            IF (IOR== 2) DR%Y = Y(JJ)   + VENT_OFFSET*DY(JJ+1)
            IF (IOR==-2) DR%Y = Y(JJ-1) - VENT_OFFSET*DY(JJ-1)
            CALL RANDOM_NUMBER(RN)
            DR%X = X(II-1) + DX(II)*RN
            CALL RANDOM_NUMBER(RN)
            DR%Z = Z(KK-1) + DZ(KK)*RN
         ENDIF
         IF (ABS(IOR)==3) THEN
            IF (IOR== 3) DR%Z = Z(KK)   + VENT_OFFSET*DZ(KK+1)
            IF (IOR==-3) DR%Z = Z(KK-1) - VENT_OFFSET*DZ(KK-1)
            CALL RANDOM_NUMBER(RN)
            DR%X = X(II-1) + DX(II)*RN
            CALL RANDOM_NUMBER(RN)
            DR%Y = Y(JJ-1) + DY(JJ)*RN
         ENDIF
    
         SELECT CASE(IOR)  ! Give particles an initial velocity (if desired)
            CASE( 1) 
               DR%U = -UW(IW)
            CASE(-1) 
               DR%U =  UW(IW)
            CASE( 2) 
               DR%V = -UW(IW)
            CASE(-2) 
               DR%V =  UW(IW)
            CASE( 3)
               DR%W = -UW(IW)
            CASE(-3) 
               DR%W =  UW(IW)
         END SELECT
    
         ! Save the insertion time (TP) and scalar property (SP) for the particle
    
         DR%A_X    = 0._EB
         DR%A_Y    = 0._EB
         DR%A_Z    = 0._EB
         DR%CLASS  = IPC
         DR%TAG    = PARTICLE_TAG
         IF (MOD(NLP,PC%SAMPLING)==0) THEN
            DR%SHOW = .TRUE.
         ELSE
            DR%SHOW = .FALSE.
         ENDIF
         DR%SPLAT   = .FALSE.
         DR%WALL_INDEX = 0
         DR%R   = 0.0_EB
         DR%TMP = PC%TMP_INITIAL
         DR%T   = T
         DR%PWT = 1._EB
         DR%IOR = 0
         IF (PC%DIAMETER > 0._EB) THEN
            IF (PC%MONODISPERSE) THEN
               DR%R   = 0.5_EB*PC%DIAMETER
               DR%PWT = 1._EB
            ELSE
               CALL RANDOM_NUMBER(RN)            
               STRATUM = NINT(REAL(NSTRATA,EB)*RN+0.5_EB)
               IL = PC%IL_CDF(STRATUM)
               IU = PC%IU_CDF(STRATUM)
               CALL RANDOM_CHOICE(PC%CDF(IL:IU),PC%R_CDF(IL:IU),IU-IL,DR%R)
               DR%PWT = PC%W_CDF(STRATUM)
               IF (2._EB*DR%R > PC%MAXIMUM_DIAMETER) THEN
                  DR%PWT = DR%PWT*DR%R**3/(0.5_EB*PC%MAXIMUM_DIAMETER)**3
                  DR%R = 0.5_EB*PC%MAXIMUM_DIAMETER
               ENDIF
            ENDIF
            MASS_SUM = MASS_SUM + DR%PWT*PC%FTPR*DR%R**3
         ENDIF
   
      ENDDO PARTICLE_INSERT_LOOP

      IF (MASS_SUM > 0._EB) THEN
      IF (SF%PARTICLE_MASS_FLUX > 0._EB) THEN
            IF (ABS(TW(IW)-T_BEGIN)<=ZERO_P .AND. SF%RAMP_INDEX(TIME_PART)>=1) THEN
               TSI = T
            ELSE
               TSI = T - TW(IW)
            ENDIF
            FLOW_RATE = EVALUATE_RAMP(TSI,SF%TAU(TIME_PART),SF%RAMP_INDEX(TIME_PART))*SF%PARTICLE_MASS_FLUX
            DROPLET(NLP-NPPCW(IW)+1:NLP)%PWT = &
            DROPLET(NLP-NPPCW(IW)+1:NLP)%PWT*FLOW_RATE*AREA_ADJUST(IW)*AW(IW)*SF%DT_INSERT/MASS_SUM
         ENDIF
      ENDIF

   ENDDO WALL_INSERT_LOOP

   ! Loop over all INIT lines and look for particles inserted within a specified volume

   VOLUME_INSERT_LOOP: DO IB=1,N_INIT

      IN => INITIALIZATION(IB)

      ! Determine if the INITIALIZATION type involves particles or droplets. If not, cycle.

      IPC = IN%PART_INDEX
      IF (IPC<1) CYCLE VOLUME_INSERT_LOOP
      IF (IN%SINGLE_INSERTION .AND. IN%ALREADY_INSERTED) CYCLE VOLUME_INSERT_LOOP

      ! Determine if the particles/droplets are controlled by devices

      PC => PARTICLE_CLASS(IPC)
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

      ! If it is not time to insert particles/droplets for this INITIALIZATION block, cycle.

      IF (T < IN%PARTICLE_INSERT_CLOCK(NM)) CYCLE VOLUME_INSERT_LOOP
 
      ! Start processing the INITIALIZATION info

      N_INITIAL = IN%NUMBER_INITIAL_DROPLETS

      IF (N_INITIAL==0) CYCLE VOLUME_INSERT_LOOP

      MASS_PER_VOLUME = IN%MASS_PER_VOLUME
      MASS_PER_TIME   = IN%MASS_PER_TIME

      SELECT CASE(IN%SHAPE)
         CASE('BLOCK')
            IF (IN%X1>XF .OR. IN%X2<XS .OR. IN%Y1>YF .OR. IN%Y2<YS .OR. IN%Z1>ZF .OR. IN%Z2<ZS) CYCLE VOLUME_INSERT_LOOP
            X1 = MAX(IN%X1,XS) 
            X2 = MIN(IN%X2,XF)
            Y1 = MAX(IN%Y1,YS) 
            Y2 = MIN(IN%Y2,YF)
            Z1 = MAX(IN%Z1,ZS) 
            Z2 = MIN(IN%Z2,ZF)
            BLOCK_VOLUME = (X2-X1)*(Y2-Y1)*(Z2-Z1)
            IF (BLOCK_VOLUME<=0._EB .AND. (MASS_PER_VOLUME>0._EB .OR. MASS_PER_TIME>0._EB)) CYCLE VOLUME_INSERT_LOOP
         CASE('CONE')
         CASE('POINT')
            IF (IN%X0>XF .OR. IN%X0<XS .OR. IN%Y0>YF .OR. IN%Y0<YS .OR. IN%Z0>ZF .OR. IN%Z0<ZS) CYCLE VOLUME_INSERT_LOOP
      END SELECT
   
      ! Assign properties to the initial droplets/particles
   
      IN%ALREADY_INSERTED = .TRUE.
      MASS_SUM = 0._EB
      INSERT_PARTICLE_LOOP: DO I=1,N_INITIAL
         NLP = NLP + 1
         PARTICLE_TAG = PARTICLE_TAG + NMESHES
         IF (NLP>NLPDIM) THEN
            CALL RE_ALLOCATE_DROPLETS(1,NM,0,1000)
            DROPLET=>MESHES(NM)%DROPLET
         ENDIF
         DR=>DROPLET(NLP)
         
         POINTWISE_IF: IF (IN%POINTWISE_DROPLET_INIT) THEN
            DR%X = IN%X0
            DR%Y = IN%Y0
            DR%Z = IN%Z0
            DR%U = IN%U0
            DR%V = IN%V0
            DR%W = IN%W0
            DR%R = IN%RADIUS
            DR%TMP = IN%TEMPERATURE
            IF (DR%R>0._EB) THEN
               DR%PWT = IN%DT_INSERT*IN%MASS_PER_TIME/(PC%FTPR*DR%R**3)
            ELSE
               DR%PWT = 0._EB
            ENDIF
         ELSE POINTWISE_IF
            BLOCK_OUT_LOOP:  DO
           !!  CALL RANDOM_CONE(DR%X,DR%Y,DR%Z,0.5_EB*(X2+X1),0.5_EB*(Y2+Y1),Z1,0.25_EB,0.5_EB)
               CALL RANDOM_RECTANGLE(DR%X,DR%Y,DR%Z,X1,X2,Y1,Y2,Z1,Z2)
               II = CELLSI(FLOOR((DR%X-XS)*RDXINT)) + 1._EB
               JJ = CELLSJ(FLOOR((DR%Y-YS)*RDYINT)) + 1._EB
               KK = CELLSK(FLOOR((DR%Z-ZS)*RDZINT)) + 1._EB
               IF (.NOT.SOLID(CELL_INDEX(II,JJ,KK))) EXIT BLOCK_OUT_LOOP
            ENDDO BLOCK_OUT_LOOP
            DR%U   = 0._EB                     ! No initial velocity
            DR%V   = 0._EB
            DR%W   = 0._EB
            DR%R   = 0._EB                     ! Radius is zero unless DIAMETER has been specified
            DR%TMP = PC%TMP_INITIAL            ! Initial temperature
            DR%PWT = 1._EB                     ! Weighting factor is one unless changed below
         ENDIF POINTWISE_IF
         
         DR%T   = T                         ! Insertion time is current time
         DR%IOR = 0                         ! Orientation of solid surface (0 means the droplet/particle is not attached)
         DR%CLASS = IPC                     ! Class identifier
         DR%TAG   = PARTICLE_TAG            ! Unique integer tag
         IF (MOD(NLP,PC%SAMPLING)==0) THEN  ! Decide whether to show or not show in Smokeview
            DR%SHOW = .TRUE.    
         ELSE
            DR%SHOW = .FALSE.    
         ENDIF
         DR%SPLAT   = .FALSE.
         DR%WALL_INDEX = 0
    
         IF (PC%DIAMETER>0._EB .AND. .NOT.IN%POINTWISE_DROPLET_INIT) THEN
            IF (PC%MONODISPERSE) THEN
               DR%R   = 0.5_EB*PC%DIAMETER
               DR%PWT = 1._EB
            ELSE
               CALL RANDOM_NUMBER(RN)            
               STRATUM = NINT(REAL(NSTRATA,EB)*RN+0.5_EB)
               IL = PC%IL_CDF(STRATUM)
               IU = PC%IU_CDF(STRATUM)
               CALL RANDOM_CHOICE(PC%CDF(IL:IU),PC%R_CDF(IL:IU),IU-IL,DR%R)
               DR%PWT = PC%W_CDF(STRATUM)
               IF (2._EB*DR%R > PC%MAXIMUM_DIAMETER) THEN
                  DR%PWT = DR%PWT*DR%R**3/(0.5_EB*PC%MAXIMUM_DIAMETER)**3
                  DR%R = 0.5_EB*PC%MAXIMUM_DIAMETER
               ENDIF
            ENDIF
            !MASS_SUM = MASS_SUM + DR%PWT*PC%FTPR*DR%R**3
         ENDIF
         MASS_SUM = MASS_SUM + DR%PWT*PC%FTPR*DR%R**3 ! if r=0 the sum will stay 0
    
         ! Process special particles that are associated with a particular SURFace type
   
         IF (PC%SURF_INDEX>0) THEN
            SF => SURFACE(PC%SURF_INDEX)
            DR%WALL_INDEX = IN%WALL_INDEX_START + I - 1
            IF (SF%THICKNESS>0._EB) THEN
               DR%R = SF%THICKNESS
            ELSEIF (SF%RADIUS>0._EB) THEN
               DR%R = SF%RADIUS
            ENDIF
            DR%PWT = 1._EB
            XW(DR%WALL_INDEX) = DR%X
            YW(DR%WALL_INDEX) = DR%Y
            ZW(DR%WALL_INDEX) = DR%Z
            IJKW(1,DR%WALL_INDEX) = II
            IJKW(2,DR%WALL_INDEX) = JJ
            IJKW(3,DR%WALL_INDEX) = KK
            IJKW(4,DR%WALL_INDEX) = 0
            IJKW(5,DR%WALL_INDEX) = PC%SURF_INDEX
            IJKW(6,DR%WALL_INDEX) = II
            IJKW(7,DR%WALL_INDEX) = JJ
            IJKW(8,DR%WALL_INDEX) = KK
         ENDIF
   
         IF (PC%N_SPLIT>1) THEN
            DR%SPLIT_IOR = MOD(I-1,PC%N_SPLIT)+1
            DR%PWT = DR%PWT/REAL(PC%N_SPLIT,EB)
         ENDIF            
         
      ENDDO INSERT_PARTICLE_LOOP
    
      POINTWISE_IF2: IF (.NOT.IN%POINTWISE_DROPLET_INIT) THEN
       
         ! Adjust particle weighting factor PWT so that desired MASS_PER_VOLUME is achieved
      
         IF (MASS_PER_TIME>0._EB) MASS_PER_VOLUME = MASS_PER_TIME*IN%DT_INSERT/BLOCK_VOLUME

         IF (MASS_PER_VOLUME>0._EB) DROPLET(NLP-N_INITIAL+1:NLP)%PWT = & 
                                    DROPLET(NLP-N_INITIAL+1:NLP)%PWT*MASS_PER_VOLUME*BLOCK_VOLUME/MASS_SUM
                                    
      ENDIF POINTWISE_IF2

   ENDDO VOLUME_INSERT_LOOP

   ! Reset particle/droplet insertion clocks

   DO N=1,N_INIT
      IN => INITIALIZATION(N)
      IF (IN%SINGLE_INSERTION) CYCLE      
      IF (T >= IN%PARTICLE_INSERT_CLOCK(NM)) IN%PARTICLE_INSERT_CLOCK(NM) = IN%PARTICLE_INSERT_CLOCK(NM) + IN%DT_INSERT
      IF (T >= IN%PARTICLE_INSERT_CLOCK(NM)) INSERT_ANOTHER_BATCH = .TRUE.
   ENDDO

   DO N=1,N_SURF
      SF => SURFACE(N)
      IPC =  SF%PART_INDEX
      IF (IPC < 1)    CYCLE
      IF (T >= SF%PARTICLE_INSERT_CLOCK(NM)) SF%PARTICLE_INSERT_CLOCK(NM) = SF%PARTICLE_INSERT_CLOCK(NM) + SF%DT_INSERT
      IF (T >= SF%PARTICLE_INSERT_CLOCK(NM)) INSERT_ANOTHER_BATCH = .TRUE.
   ENDDO

!   DO N=1,N_PROP
!      PY => PROPERTY(N)
!      IF (T >= PY%PARTICLE_INSERT_CLOCK(NM)) PY%PARTICLE_INSERT_CLOCK(NM) = PY%PARTICLE_INSERT_CLOCK(NM) + PY%DT_INSERT
!      IF (T >= PY%PARTICLE_INSERT_CLOCK(NM)) INSERT_ANOTHER_BATCH = .TRUE.
!   ENDDO

   IF (.NOT.INSERT_ANOTHER_BATCH) EXIT OVERALL_INSERT_LOOP

ENDDO OVERALL_INSERT_LOOP
 
TUSED(8,NM)=TUSED(8,NM)+SECOND()-TNOW
END SUBROUTINE INSERT_DROPLETS_AND_PARTICLES
 
 

SUBROUTINE RANDOM_CHOICE(CDF,VAR,NPTS,CHOICE)
 
INTEGER,  INTENT(IN)  :: NPTS
REAL(EB), INTENT(IN)  :: CDF(0:NPTS),VAR(0:NPTS)
REAL(EB), INTENT(OUT) :: CHOICE
INTEGER  :: IT
REAL(EB) :: CFRAC,RN,A,B
 
CALL RANDOM_NUMBER(RN)
A = MINVAL(CDF)
B = MAXVAL(CDF)
RN = A + (B-A)*RN
 
CDF_LOOP: DO IT=1,NPTS
   IF (CDF(IT) > RN) THEN
      CFRAC  = (RN-CDF(IT-1))/(CDF(IT)-CDF(IT-1))
      CHOICE = VAR(IT-1) + (VAR(IT)-VAR(IT-1))*CFRAC  
      EXIT CDF_LOOP
   ENDIF
ENDDO CDF_LOOP
 
END SUBROUTINE RANDOM_CHOICE
 
 

SUBROUTINE UPDATE_PARTICLES(T,NM)

USE COMP_FUNCTIONS, ONLY : SECOND  
REAL(EB), INTENT(IN) :: T
INTEGER, INTENT(IN) :: NM
REAL(EB) :: TNOW

! Return if this is an evacuation mesh
IF (EVACUATION_ONLY(NM)) RETURN  

! Zero out the number of the droplets in the "orphanage"; that is, the place to hold droplets transferring from mesh to mesh

MESHES(NM)%OMESH(1:NMESHES)%N_DROP_ORPHANS = 0

! Return if there are no particles in this mesh

IF (MESHES(NM)%NLP==0)   RETURN

! Set the CPU timer and point to the current mesh variables

TNOW=SECOND()
CALL POINT_TO_MESH(NM)

! Zero out the contribution by lagrangian particles to divergence

IF (N_EVAP_INDICES>0 .AND. .NOT.EVACUATION_ONLY(NM) .AND. CORRECTOR) D_LAGRANGIAN = 0._EB

! Move the droplets/particles, then compute mass and energy transfer, then add droplet momentum to gas

IF (CORRECTOR) CALL MOVE_PARTICLES(T,NM)
IF (CORRECTOR) CALL PARTICLE_MASS_ENERGY_TRANSFER(T,NM)
CALL PARTICLE_MOMENTUM_TRANSFER

TUSED(8,NM)=TUSED(8,NM)+SECOND()-TNOW
END SUBROUTINE UPDATE_PARTICLES



SUBROUTINE MOVE_PARTICLES(T,NM)

! Momentum transfer from all particles and droplets
 
USE COMP_FUNCTIONS, ONLY : SECOND  
USE MATH_FUNCTIONS, ONLY : EVALUATE_RAMP, AFILL2
USE PHYSICAL_FUNCTIONS, ONLY : DRAG, GET_VISCOSITY
 
REAL(EB) :: RHO_G,RVC,RDS,RDC,QREL,UREL,VREL,WREL,TMP_G,RN,THETA_RN, &
            RD,T,C_DRAG,XI,YJ,ZK,MU_AIR,WE_G,DROP_VOL_FRAC,DROP_DEN, &
            UBAR,VBAR,WBAR,GRVT1,GRVT2,GRVT3, &
            UVW,DUMMY=0._EB,X_OLD,Y_OLD,Z_OLD,STEP_FRACTION(-3:3),SURFACE_DROPLET_DIAMETER, &
            T_BU_BAG,T_BU_STRIP,B_1,THROHALF,P_UVWMAX,UVWMAX, &
            ALPHA,BETA,DR_MASS,FP_MASS,OBDT,OPA,BDTOA,U_OLD,V_OLD,W_OLD,MPOM,RDT,&
            WAKE_VEL,RE_WAKE,SFAC
REAL(EB) :: CHILD_RADIUS(0:NDC),ZZ_GET(0:N_TRACKED_SPECIES)
LOGICAL :: HIT_SOLID
INTEGER :: ICN,I,IIN,JJN,KKN,II,JJ,KK,IIX,JJY,KKZ,IW,IWP1,IWM1,IWP2,IWM2,IWP3,IWM3,IOR_OLD,IC,IOR_FIRST,IML
INTEGER, INTENT(IN) :: NM
TYPE (DROPLET_TYPE), POINTER :: DR=>NULL()
TYPE (PARTICLE_CLASS_TYPE), POINTER :: PC=>NULL()
REAL(EB), POINTER, DIMENSION(:,:,:) :: NDPC=>NULL() ! number of droplets per cell

SURFACE_DROPLET_DIAMETER = 0.001_EB  ! All droplets adjusted to this size when on solid (m)
THROHALF = (0.5_EB)**(1./3.)

B_1 =  1.7321_EB ! SQRT(3)

RDT = 1._EB/DT

GRVT1 = -EVALUATE_RAMP(T,DUMMY,I_RAMP_GX)*GVEC(1) 
GRVT2 = -EVALUATE_RAMP(T,DUMMY,I_RAMP_GY)*GVEC(2) 
GRVT3 = -EVALUATE_RAMP(T,DUMMY,I_RAMP_GZ)*GVEC(3) 

UVWMAX = 0._EB
P_UVWMAX = UVWMAX

! Sum up the number of droplets/particles in each grid cell (NDPC -- Number Droplets Per Cell)

NDPC=>WORK1
NDPC=0._EB

DO I=1,NLP
   DR => DROPLET(I)
   XI = CELLSI(FLOOR((DR%X-XS)*RDXINT))
   YJ = CELLSJ(FLOOR((DR%Y-YS)*RDYINT))
   ZK = CELLSK(FLOOR((DR%Z-ZS)*RDZINT))
   II  = FLOOR(XI+1._EB)
   JJ  = FLOOR(YJ+1._EB)
   KK  = FLOOR(ZK+1._EB)
   IF (DR%PWT>0._EB) NDPC(II,JJ,KK)=NDPC(II,JJ,KK)+1._EB
ENDDO

! Loop through all Lagrangian particles and move them one time step

DROPLET_LOOP: DO I=1,NLP  

   ! Assign particle (DR%) and particle class (PC%) shortcuts

   DR => DROPLET(I)
   PC => PARTICLE_CLASS(DR%CLASS)

   ! Determine particle radius

   RD  = DR%R
   IF (.NOT. PC%MASSLESS .AND. RD<=0._EB) CYCLE DROPLET_LOOP
   RDS = RD*RD
   RDC = RD*RDS

   ! Determine the current coordinates of the particle

   XI = CELLSI(FLOOR((DR%X-XS)*RDXINT))
   YJ = CELLSJ(FLOOR((DR%Y-YS)*RDYINT))
   ZK = CELLSK(FLOOR((DR%Z-ZS)*RDZINT))
   II = FLOOR(XI+1._EB)
   JJ = FLOOR(YJ+1._EB)
   KK = FLOOR(ZK+1._EB)
   IC = CELL_INDEX(II,JJ,KK)

   X_OLD = DR%X
   Y_OLD = DR%Y
   Z_OLD = DR%Z
   U_OLD = DR%U
   V_OLD = DR%V
   W_OLD = DR%W

   ! Throw out particles that are inside a solid obstruction

   IF (SOLID(IC)) THEN
      DR%X = 1.E6_EB
      CYCLE DROPLET_LOOP
   ENDIF

   ! Calculate the particle CFL number
    
   UVW = MAX( ABS(DR%U)*RDX(II),ABS(DR%V)*RDY(JJ),ABS(DR%W)*RDZ(KK) )
   P_UVWMAX = MAX(P_UVWMAX,UVW)

   ! Interpolate the nearest velocity components of the gas

   IIX  = FLOOR(XI+.5_EB)
   JJY  = FLOOR(YJ+.5_EB)
   KKZ  = FLOOR(ZK+.5_EB)
   UBAR = AFILL2(U,II-1,JJY,KKZ,(DR%X-X(II-1))*RDX(II),YJ-JJY+.5_EB,ZK-KKZ+.5_EB)
   VBAR = AFILL2(V,IIX,JJ-1,KKZ,XI-IIX+.5_EB,(DR%Y-Y(JJ-1))*RDY(JJ),ZK-KKZ+.5_EB)
   WBAR = AFILL2(W,IIX,JJY,KK-1,XI-IIX+.5_EB,YJ-JJY+.5_EB,(DR%Z-Z(KK-1))*RDZ(KK))
    
   ! If the particle is massless, just move it and go on to the next particle

   IF (PC%MASSLESS .OR. DR%PWT<=ZERO_P) THEN
      DR%U = UBAR
      DR%V = VBAR
      DR%W = WBAR 
      DR%X = DR%X + DR%U*DT
      DR%Y = DR%Y + DR%V*DT
      DR%Z = DR%Z + DR%W*DT
      CYCLE DROPLET_LOOP
   ENDIF

   ! Calculate the particle drag coefficient

   RVC    = RDX(II)*RRN(II)*RDY(JJ)*RDZ(KK)
   UREL   = DR%U - UBAR
   VREL   = DR%V - VBAR
   WREL   = DR%W - WBAR
   QREL   = SQRT(UREL*UREL + VREL*VREL + WREL*WREL)
   TMP_G  = MAX(TMPMIN,TMP(II,JJ,KK))
   RHO_G  = RHO(II,JJ,KK)
   IF (N_TRACKED_SPECIES>0) ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(II,JJ,KK,1:N_TRACKED_SPECIES)
   CALL GET_VISCOSITY(ZZ_GET,MU_AIR,TMP_G)
   DR%RE  = RHO_G*QREL*2._EB*RD/MU_AIR

   IF (PC%DRAG_LAW==USER_DRAG) THEN
      C_DRAG = PC%USER_DRAG_COEFFICIENT
   ELSE
      C_DRAG = DRAG(DR%RE,PC%DRAG_LAW)
   ENDIF

   ! Drag reduction model, except for particles associated with a SURF line

   WAKE_VEL=1.0_EB
   IF (PC%SURF_INDEX<1) THEN
      DROP_DEN      = AVG_DROP_DEN(II,JJ,KK,PC%EVAP_INDEX) 
      DROP_VOL_FRAC = DROP_DEN/PC%DENSITY 
      IF (DROP_VOL_FRAC > PC%DENSE_VOLUME_FRACTION) CALL WAKE_REDUCTION(DROP_VOL_FRAC,DR%RE,C_DRAG,WAKE_VEL)
   ENDIF

  ! Secondary break-up model

   BREAKUP: IF (PC%BREAKUP) THEN
      ! Use undisturbed wake velocity for breakup calculations
      WAKE_VEL    = WAKE_VEL*QREL
      RE_WAKE     = RHO_G*WAKE_VEL   *2._EB*RD/MU_AIR
      WE_G        = RHO_G*WAKE_VEL**2*2._EB*RD/PC%SURFACE_TENSION
      ! Shape Deformation
      C_DRAG = SHAPE_DEFORMATION(RE_WAKE,WE_G,C_DRAG)
      ! Breakup conditions according to WAVE model by Reitz (1987)
      T_BU_BAG    = T_END-T_BEGIN
      T_BU_STRIP  = T_END-T_BEGIN
      IF (WE_G >= 12.0_EB)               T_BU_BAG   = 1.72_EB*B_1*SQRT(PC%DENSITY*RDC/(2._EB*PC%SURFACE_TENSION))
      IF (WE_G/SQRT(RE_WAKE) >= 1.0_EB)  T_BU_STRIP = B_1*(RD/WAKE_VEL)*SQRT(PC%DENSITY/RHO_G)
      ! droplet age is larger than smallest characteristic breakup time
      AGE_IF: IF ((T-DR%T) > MIN(T_BU_BAG,T_BU_STRIP)) THEN
         IF (PC%MONODISPERSE) THEN
            RD    = THROHALF*RD
         ELSE
            DO WHILE (RD >= DR%R)
               CHILD_RADIUS = PC%BREAKUP_CHILD_DIAMETER*DR%R*PC%CHILD_R_CDF(:)
               CALL RANDOM_CHOICE(PC%CHILD_CDF(:),CHILD_RADIUS,NDC,RD)
            END DO
            RD = MAX(RD,1.1_EB*PC%MINIMUM_DIAMETER/2._EB)
         ENDIF
         DR%RE    = RHO_G*QREL*2._EB*RD/MU_AIR
         IF (PC%DRAG_LAW==USER_DRAG) THEN
            C_DRAG = PC%USER_DRAG_COEFFICIENT
         ELSE
            C_DRAG = DRAG(DR%RE,PC%DRAG_LAW)
         ENDIF
         DR%PWT   = DR%PWT*RDC/RD**3
         DR%T     = T
         DR%R     = RD
         RDS      = RD*RD
         RDC      = RD*RDS
         ! Redo wake reduction and shape deformation for the new drop
         ! Drag reduction, except for particles associated with a SURF line
         WAKE_VEL = 1.0_EB
         IF (PC%SURF_INDEX<1) THEN
            DROP_DEN      = AVG_DROP_DEN(II,JJ,KK,PC%EVAP_INDEX) 
            DROP_VOL_FRAC = DROP_DEN/PC%DENSITY 
            IF (DROP_VOL_FRAC > PC%DENSE_VOLUME_FRACTION) CALL WAKE_REDUCTION(DROP_VOL_FRAC,DR%RE,C_DRAG,WAKE_VEL)
         ENDIF
         ! Change in drag coefficient due to deformation of droplet shape (WE_G > 2)
         WAKE_VEL = WAKE_VEL*QREL
         RE_WAKE  = RHO_G*WAKE_VEL   *2._EB*RD/MU_AIR
         WE_G     = RHO_G*WAKE_VEL**2*2._EB*RD/PC%SURFACE_TENSION
         ! Shape Deformation
         C_DRAG   = SHAPE_DEFORMATION(RE_WAKE,WE_G,C_DRAG)
      ENDIF AGE_IF
   ENDIF BREAKUP

   ! Move airborne, non-stationary particles
         
   PARTICLE_NON_STATIC_IF: IF (.NOT.PC%STATIC .AND. DR%IOR==0) THEN

      FP_MASS = (RHO_G/RVC)/NDPC(II,JJ,KK) ! fluid parcel mass
      DR_MASS = PC%DENSITY*FOTH*PI*RDC     ! droplet mass
               
      BETA  = 0.5_EB*RHO_G*C_DRAG*PI*RDS*(1._EB/DR_MASS+1._EB/FP_MASS)*QREL
      OBDT  = 1._EB+BETA*DT
      ALPHA = FP_MASS/DR_MASS
      OPA   = 1._EB+ALPHA
      BDTOA = BETA*DT/OPA
               
      DR%U = ( U_OLD + (U_OLD+ALPHA*UBAR)*BDTOA )/OBDT
      DR%V = ( V_OLD + (V_OLD+ALPHA*VBAR)*BDTOA )/OBDT
      DR%W = ( W_OLD + (W_OLD+ALPHA*WBAR)*BDTOA )/OBDT
               
      IF (BETA>=ZERO_P) THEN
         ! fluid momentum source term
         MPOM = DR%PWT*DR_MASS/(RHO_G/RVC)
         DR%A_X = MPOM*(U_OLD-DR%U)*RDT 
         DR%A_Y = MPOM*(V_OLD-DR%V)*RDT
         DR%A_Z = MPOM*(W_OLD-DR%W)*RDT
         ! gravitational acceleration
         DR%U = DR%U + GVEC(1)*DT
         DR%V = DR%V + GVEC(2)*DT
         DR%W = DR%W + GVEC(3)*DT
         ! analytical solution for droplet position
         DR%X = X_OLD + (U_OLD+ALPHA*UBAR)/OPA*DT + ALPHA/BETA*(U_OLD-UBAR)/OPA*LOG(OBDT) + GVEC(1)*DT*DT
         DR%Y = Y_OLD + (V_OLD+ALPHA*VBAR)/OPA*DT + ALPHA/BETA*(V_OLD-VBAR)/OPA*LOG(OBDT) + GVEC(2)*DT*DT
         DR%Z = Z_OLD + (W_OLD+ALPHA*WBAR)/OPA*DT + ALPHA/BETA*(W_OLD-WBAR)/OPA*LOG(OBDT) + GVEC(3)*DT*DT
      ELSE
         ! no drag
         DR%A_X  = 0._EB
         DR%A_Y  = 0._EB
         DR%A_Z  = 0._EB
         DR%X = X_OLD + (U_OLD + GVEC(1)*DT)*DT
         DR%Y = Y_OLD + (V_OLD + GVEC(2)*DT)*DT
         DR%Z = Z_OLD + (W_OLD + GVEC(3)*DT)*DT
      ENDIF
               
   ENDIF PARTICLE_NON_STATIC_IF 

   ! Drag calculation for stationary, airborne particles

   PARTICLE_STATIC_IF: IF (PC%STATIC .AND. DR%IOR==0 .AND. .NOT. PC%TREE) THEN
      BETA = 0.5_EB*RVC*C_DRAG*(DR%PWT*PI*RDS)*QREL
      OBDT = 1._EB+BETA*DT
      DR%A_X = UBAR*(1._EB/OBDT-1._EB)*RDT 
      DR%A_Y = VBAR*(1._EB/OBDT-1._EB)*RDT
      DR%A_Z = WBAR*(1._EB/OBDT-1._EB)*RDT
   ENDIF PARTICLE_STATIC_IF

   TREE_PARTICLES: IF (DR%IOR==0 .AND. PC%TREE) THEN
      SFAC   = PC%VEG_DRAG_COEFFICIENT*PC%VEG_SV*DR%VEG_PACKING_RATIO*QREL*C_DRAG
      DR%A_X  = SFAC*UREL - DR%VEG_MLR*UBAR/RHO_G
      DR%A_Y  = SFAC*VREL - DR%VEG_MLR*VBAR/RHO_G
      DR%A_Z  = SFAC*WREL - DR%VEG_MLR*WBAR/RHO_G
   ENDIF TREE_PARTICLES

   ! Move particles/droplets attached to solids

   SOLID_IF: IF (DR%IOR/=0) THEN
      DR%A_X = 0._EB 
      DR%A_Y = 0._EB 
      DR%A_Z = 0._EB 
      DR%X = X_OLD + U_OLD*DT
      DR%Y = Y_OLD + V_OLD*DT
      DR%Z = Z_OLD + W_OLD*DT
   ENDIF SOLID_IF

   ! If the particle does not move, but does drag, go on to the next particle

   IF (PC%STATIC) CYCLE DROPLET_LOOP
    
   ! Droplet hits the floor
    
   IF (POROUS_FLOOR .AND. DR%Z<ZS .AND. PC%EVAP_INDEX>0) THEN
      IC = CELL_INDEX(II,JJ,1)
      IW = WALL_INDEX(IC,-3)
      IF (BOUNDARY_TYPE(IW)==SOLID_BOUNDARY .AND. ACCUMULATE_WATER .AND. .NOT.DR%SPLAT) THEN
         AWMPUA(IW,PC%EVAP_INDEX) = AWMPUA(IW,PC%EVAP_INDEX) + DR%PWT*PC%FTPR*RDC*RAW(IW)
         DR%SPLAT = .TRUE.
      ENDIF
      CYCLE DROPLET_LOOP
   ENDIF
       
   IF (.NOT.POROUS_FLOOR .AND. DR%Z<ZS) THEN
      IC = CELL_INDEX(II,JJ,1)
      IW = WALL_INDEX(IC,-3)
      IF (BOUNDARY_TYPE(IW)==SOLID_BOUNDARY) THEN
         IF (ACCUMULATE_WATER .AND. .NOT.DR%SPLAT) THEN
            AWMPUA(IW,PC%EVAP_INDEX) = AWMPUA(IW,PC%EVAP_INDEX) + DR%PWT*PC%FTPR*RDC*RAW(IW)
            DR%SPLAT = .TRUE.
         ENDIF
         DR%Z = ZS + 0.05_EB*DZ(1)
         DR%IOR = 3
         CALL RANDOM_NUMBER(RN)
         THETA_RN = TWOPI*RN
         DR%U = PC%HORIZONTAL_VELOCITY*COS(THETA_RN)
         DR%V = PC%HORIZONTAL_VELOCITY*SIN(THETA_RN)
         DR%W = 0._EB
      ENDIF
   ENDIF
 
   ! Where is the droplet now? Limit the location by UBOUND and LBOUND due to the possible super fast droplets

   IIN = MAX(LBOUND(CELLSI,1),MIN(UBOUND(CELLSI,1),FLOOR((DR%X-XS)*RDXINT)))
   JJN = MAX(LBOUND(CELLSJ,1),MIN(UBOUND(CELLSJ,1),FLOOR((DR%Y-YS)*RDYINT)))
   KKN = MAX(LBOUND(CELLSK,1),MIN(UBOUND(CELLSK,1),FLOOR((DR%Z-ZS)*RDZINT)))
   XI  = CELLSI(IIN)
   YJ  = CELLSJ(JJN)
   ZK  = CELLSK(KKN)
   IIN = FLOOR(XI+1._EB)
   JJN = FLOOR(YJ+1._EB)
   KKN = FLOOR(ZK+1._EB)
   IF (IIN<0 .OR. IIN>IBP1) CYCLE DROPLET_LOOP
   IF (JJN<0 .OR. JJN>JBP1) CYCLE DROPLET_LOOP
   IF (KKN<0 .OR. KKN>KBP1) CYCLE DROPLET_LOOP
   ICN = CELL_INDEX(IIN,JJN,KKN)
   IF (IC==0 .OR. ICN==0) CYCLE DROPLET_LOOP

   IF (DR%X<XS .AND. BOUNDARY_TYPE(WALL_INDEX(IC,-1))/=SOLID_BOUNDARY) CYCLE DROPLET_LOOP
   IF (DR%X>XF .AND. BOUNDARY_TYPE(WALL_INDEX(IC, 1))/=SOLID_BOUNDARY) CYCLE DROPLET_LOOP
   IF (DR%Y<YS .AND. BOUNDARY_TYPE(WALL_INDEX(IC,-2))/=SOLID_BOUNDARY) CYCLE DROPLET_LOOP
   IF (DR%Y>YF .AND. BOUNDARY_TYPE(WALL_INDEX(IC, 2))/=SOLID_BOUNDARY) CYCLE DROPLET_LOOP
   IF (DR%Z<ZS .AND. BOUNDARY_TYPE(WALL_INDEX(IC,-3))/=SOLID_BOUNDARY) CYCLE DROPLET_LOOP
   IF (DR%Z>ZF .AND. BOUNDARY_TYPE(WALL_INDEX(IC, 3))/=SOLID_BOUNDARY) CYCLE DROPLET_LOOP

   ! If droplet hits an obstacle, change its properties

   AIR_TO_SOLID: IF (II/=IIN .OR. JJ/=JJN .OR. KK/=KKN) THEN

      IOR_OLD   = DR%IOR
      HIT_SOLID = .FALSE.

      ! Check if any solid boundaries of original grid cell have been crossed

      IWP1 = WALL_INDEX(IC, 1) 
      IWM1 = WALL_INDEX(IC,-1)
      IWP2 = WALL_INDEX(IC, 2)
      IWM2 = WALL_INDEX(IC,-2)
      IWP3 = WALL_INDEX(IC, 3)
      IWM3 = WALL_INDEX(IC,-3)
      STEP_FRACTION = 1._EB

      IF (KKN>KK .AND. BOUNDARY_TYPE(IWP3)==SOLID_BOUNDARY) THEN
         DR%IOR=-3
         HIT_SOLID = .TRUE.
         STEP_FRACTION(DR%IOR) = MAX(0._EB,(Z(KK)-Z_OLD-0.05_EB*DZ(KK))/(DR%Z-Z_OLD))
      ENDIF
      IF (KKN<KK .AND. BOUNDARY_TYPE(IWM3)==SOLID_BOUNDARY) THEN
         DR%IOR= 3
         HIT_SOLID = .TRUE.
         STEP_FRACTION(DR%IOR) = MAX(0._EB,(Z(KK-1)-Z_OLD+0.05_EB*DZ(KK-1))/(DR%Z-Z_OLD))
      ENDIF
      IF (IIN>II .AND. BOUNDARY_TYPE(IWP1)==SOLID_BOUNDARY) THEN
         DR%IOR=-1
         HIT_SOLID = .TRUE.
         STEP_FRACTION(DR%IOR) = MAX(0._EB,(X(II)-X_OLD-0.05_EB*DX(II))/(DR%X-X_OLD))
      ENDIF
      IF (IIN<II .AND. BOUNDARY_TYPE(IWM1)==SOLID_BOUNDARY) THEN
         DR%IOR= 1
         HIT_SOLID = .TRUE.
         STEP_FRACTION(DR%IOR) = MAX(0._EB,(X(II-1)-X_OLD+0.05_EB*DX(II-1))/(DR%X-X_OLD))
      ENDIF
      IF (JJN>JJ .AND. BOUNDARY_TYPE(IWP2)==SOLID_BOUNDARY) THEN
         DR%IOR=-2
         HIT_SOLID = .TRUE.
         STEP_FRACTION(DR%IOR) = MAX(0._EB,(Y(JJ)-Y_OLD-0.05_EB*DY(JJ))/(DR%Y-Y_OLD))
      ENDIF
      IF (JJN<JJ .AND. BOUNDARY_TYPE(IWM2)==SOLID_BOUNDARY) THEN
         DR%IOR= 2
         HIT_SOLID = .TRUE.
         STEP_FRACTION(DR%IOR) = MAX(0._EB,(Y(JJ-1)-Y_OLD+0.05_EB*DY(JJ-1))/(DR%Y-Y_OLD))
      ENDIF

      IF (DR%IOR/=0 .AND. .NOT.ALLOW_SURFACE_DROPLETS) THEN
         DR%R = 0.9_EB*PC%KILL_RADIUS
         CYCLE DROPLET_LOOP
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
      DR%WALL_INDEX = WALL_INDEX(IC,-IOR_FIRST)

      ! If no solid boundaries of original cell have been crossed, check boundaries of new grid cell
 
      IF (DR%WALL_INDEX==0) THEN
         IWP1 = WALL_INDEX(ICN, 1)
         IWM1 = WALL_INDEX(ICN,-1)
         IWP2 = WALL_INDEX(ICN, 2)
         IWM2 = WALL_INDEX(ICN,-2)
         IWP3 = WALL_INDEX(ICN, 3)
         IWM3 = WALL_INDEX(ICN,-3)
         HIT_SOLID = .FALSE.
         STEP_FRACTION = 1._EB
         IF (KKN>KK .AND. BOUNDARY_TYPE(IWM3)==SOLID_BOUNDARY) THEN
            DR%IOR=-3
            HIT_SOLID = .TRUE.
            STEP_FRACTION(DR%IOR) = MAX(0._EB,(Z(KK)-Z_OLD-0.05_EB*DZ(KK))/(DR%Z-Z_OLD))
         ENDIF
         IF (KKN<KK .AND. BOUNDARY_TYPE(IWP3)==SOLID_BOUNDARY) THEN
            DR%IOR= 3
            HIT_SOLID = .TRUE.
            STEP_FRACTION(DR%IOR) = MAX(0._EB,(Z(KK-1)-Z_OLD+0.05_EB*DZ(KK-1))/(DR%Z-Z_OLD))
         ENDIF
         IF (IIN>II .AND. BOUNDARY_TYPE(IWM1)==SOLID_BOUNDARY) THEN
            DR%IOR=-1
            HIT_SOLID = .TRUE.
            STEP_FRACTION(DR%IOR) = MAX(0._EB,(X(II)-X_OLD-0.05_EB*DX(II))/(DR%X-X_OLD))
         ENDIF
         IF (IIN<II .AND. BOUNDARY_TYPE(IWP1)==SOLID_BOUNDARY) THEN
            DR%IOR= 1
            HIT_SOLID = .TRUE.
            STEP_FRACTION(DR%IOR) = MAX(0._EB,(X(II-1)-X_OLD+0.05_EB*DX(II-1))/(DR%X-X_OLD))
         ENDIF
         IF (JJN>JJ .AND. BOUNDARY_TYPE(IWM2)==SOLID_BOUNDARY) THEN
            DR%IOR=-2
            HIT_SOLID = .TRUE.
            STEP_FRACTION(DR%IOR) = MAX(0._EB,(Y(JJ)-Y_OLD-0.05_EB*DY(JJ))/(DR%Y-Y_OLD))
         ENDIF
         IF (JJN<JJ .AND. BOUNDARY_TYPE(IWP2)==SOLID_BOUNDARY) THEN
            DR%IOR= 2
            HIT_SOLID = .TRUE.
            STEP_FRACTION(DR%IOR) = MAX(0._EB,(Y(JJ-1)-Y_OLD+0.05_EB*DY(JJ-1))/(DR%Y-Y_OLD))
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
         DR%WALL_INDEX = WALL_INDEX(ICN,IOR_FIRST)
      ENDIF

      ! Check if droplet has crossed no solid planes or too many

      IF_HIT_SOLID: IF (HIT_SOLID) THEN

         IF (DR%WALL_INDEX==0) CYCLE DROPLET_LOOP

         ! Add droplet mass to accumulated liquid array

         IF (ACCUMULATE_WATER .AND. HIT_SOLID .AND. .NOT.DR%SPLAT) THEN
            AWMPUA(DR%WALL_INDEX,PC%EVAP_INDEX)=AWMPUA(DR%WALL_INDEX,PC%EVAP_INDEX)+DR%PWT*PC%FTPR*DR%R**3*RAW(DR%WALL_INDEX)
            DR%SPLAT = .TRUE.
         ENDIF

         ! Adjust the size of the droplet and weighting factor 

         DR%R   = MIN(0.5_EB*SURFACE_DROPLET_DIAMETER,(DR%PWT*RDC)**ONTH)
         DR%PWT = DR%PWT*RDC/DR%R**3

         ! Move particle to where it almost hits solid

         DR%X = X_OLD + MINVAL(STEP_FRACTION)*DT*DR%U
         DR%Y = Y_OLD + MINVAL(STEP_FRACTION)*DT*DR%V
         DR%Z = Z_OLD + MINVAL(STEP_FRACTION)*DT*DR%W
         XI  = CELLSI(FLOOR((DR%X-XS)*RDXINT))
         YJ  = CELLSJ(FLOOR((DR%Y-YS)*RDYINT))
         ZK  = CELLSK(FLOOR((DR%Z-ZS)*RDZINT))
         IIN = FLOOR(XI+1._EB)
         JJN = FLOOR(YJ+1._EB)
         KKN = FLOOR(ZK+1._EB)
         ICN = CELL_INDEX(IIN,JJN,KKN)
         IF (IOR_OLD==DR%IOR) CYCLE DROPLET_LOOP

         ! Check if droplet has not found surface. Simply remove for now. Todo: search algorith

         IW = WALL_INDEX(ICN, -DR%IOR)
         IF (BOUNDARY_TYPE(IW)==NULL_BOUNDARY) THEN
            DR%R = 0.9_EB*PC%KILL_RADIUS
            CYCLE DROPLET_LOOP
         ENDIF

         ! Choose a direction for the droplets to move

         DIRECTION: SELECT CASE(DR%IOR)
            CASE (-2:-1,1:2) DIRECTION  
               DR%U = 0._EB
               DR%V = 0._EB
               DR%W = -PC%VERTICAL_VELOCITY 
            CASE (-3) DIRECTION 
               IF (.NOT.ALLOW_UNDERSIDE_DROPLETS) THEN 
               DR%U = 0._EB
               DR%V = 0._EB
               DR%W = -PC%VERTICAL_VELOCITY 
               DR%IOR = 0
               ELSE 
                  CALL RANDOM_NUMBER(RN)
                  THETA_RN = TWOPI*RN
                  DR%U = PC%HORIZONTAL_VELOCITY*COS(THETA_RN)
                  DR%V = PC%HORIZONTAL_VELOCITY*SIN(THETA_RN)
                  DR%W = 0._EB
               ENDIF
            CASE (3) DIRECTION
               CALL RANDOM_NUMBER(RN)
               THETA_RN = TWOPI*RN
               DR%U = PC%HORIZONTAL_VELOCITY*COS(THETA_RN)
               DR%V = PC%HORIZONTAL_VELOCITY*SIN(THETA_RN)
               DR%W = 0._EB
         END SELECT DIRECTION

      ENDIF IF_HIT_SOLID

   ENDIF AIR_TO_SOLID 

   ! Check if droplets that were attached to a solid are still attached after the time update

   IW = WALL_INDEX(ICN, -DR%IOR)

   IF (BOUNDARY_TYPE(IW)/=SOLID_BOUNDARY) THEN
      SELECT CASE(DR%IOR)
         CASE( 1)
            DR%X = DR%X - 0.2_EB*DX(II)
            DR%W = -DR%W
         CASE(-1)
            DR%X = DR%X + 0.2_EB*DX(II)
            DR%W = -DR%W
         CASE( 2)
            DR%Y = DR%Y - 0.2_EB*DY(JJ)
            DR%W = -DR%W
         CASE(-2)
            DR%Y = DR%Y + 0.2_EB*DY(JJ)
            DR%W = -DR%W
         CASE( 3) ! Particle has reached the edge of a horizontal surface
            DR%U = -DR%U
            DR%V = -DR%V
            DR%Z =  DR%Z - 0.2_EB*DZ(KK)
         CASE(-3)
      END SELECT
   ENDIF

   IF (DR%IOR/=0 .AND. BOUNDARY_TYPE(IW)/=SOLID_BOUNDARY) THEN
      DR%IOR = 0
      DR%WALL_INDEX = 0
      IF (RECOUNT_DRIP) DR%SPLAT=.FALSE.
   ELSE
      DR%WALL_INDEX = WALL_INDEX(ICN,-DR%IOR)
   ENDIF

ENDDO DROPLET_LOOP

UVWMAX = MAX(UVWMAX,P_UVWMAX)

PART_CFL = DT*UVWMAX

! Remove out-of-bounds particles
 
CALL REMOVE_DROPLETS(T,NM)

CONTAINS

SUBROUTINE WAKE_REDUCTION(DROP_VOL_FRAC,RE,C_DRAG,WAKE_VEL)

! Compute C_DRAG reduction due to the wake effect (Ramirez, Munoz et al. 2007)

REAL(EB)DROP_VOL_FRAC, RE, C_DRAG
REAL(EB) WAKE_VEL, LODM, RELOD

IF (DROP_VOL_FRAC > 0._EB) THEN
   LODM = (PI/(6._EB*DROP_VOL_FRAC))**(1./3.) - 0.5_EB
   RELOD = RE/(16._EB * LODM)
   WAKE_VEL = 1._EB - 0.5_EB*C_DRAG*(1._EB - EXP(-RELOD))
   WAKE_VEL = MAX(WAKE_VEL,0.15_EB)
   C_DRAG = C_DRAG * WAKE_VEL * (1._EB + (RELOD/LODM)*EXP(-RELOD))
ELSE
   WAKE_VEL=1.0_EB
ENDIF

RETURN
END SUBROUTINE WAKE_REDUCTION

! SHAPE DEFORMATION Loth, 2008 
! E.Loth, Quasi-steady shape and drag of deformable bubbles and drops, International Journal of Multiphase Flow 34 (2008)
REAL(EB) FUNCTION SHAPE_DEFORMATION(RE,WE,C_DRAG)
REAL(EB)::RE,WE,C_DRAG,C_DRAGNEW,E
REAL(EB)::DC_DRAGSTAR,fSN,WERE02

IF(WE>2.0) THEN
    WERE02=WE*RE**0.2
    DC_DRAGSTAR=.38E-2*WERE02+3E-5*WERE02**2+9E-7*WERE02**3
    fSN=1.0+0.15*RE**.0687
    C_DRAGNEW=1.0/(3.0*RE)*(DC_DRAGSTAR*(8*RE+72-72*fSN)+72*fSN)
    C_DRAGNEW=MIN(8.0_EB/3.0_EB,C_DRAGNEW) ! Bounded from above by drag of a disintegrating drop
    C_DRAGNEW=MAX(C_DRAG,C_DRAGNEW)
    ! Absorb the effect of the larger projected surface area into C_DRAG.
    ! Particle movement routines use projected area of a sphere,
    ! calculate the ratio of projected surface areas of a sphere and an
    ! ellipsoid of the same volume with aspect ratio E.
    E=1._EB-0.75_EB*tanh(0.07_EB*WE)
    SHAPE_DEFORMATION=C_DRAGNEW*E**(-2.0/3.0)
ELSE
    SHAPE_DEFORMATION=C_DRAG
ENDIF    

END FUNCTION SHAPE_DEFORMATION
END SUBROUTINE MOVE_PARTICLES



SUBROUTINE PARTICLE_MASS_ENERGY_TRANSFER(T,NM)
    
! Mass and energy transfer between gas and droplets

USE PHYSICAL_FUNCTIONS, ONLY : GET_MASS_FRACTION,GET_AVERAGE_SPECIFIC_HEAT,GET_MOLECULAR_WEIGHT,GET_SPECIFIC_GAS_CONSTANT,&
                               SURFACE_DENSITY,GET_SPECIFIC_HEAT,GET_MASS_FRACTION_ALL
USE COMP_FUNCTIONS, ONLY: SHUTDOWN
USE OUTPUT_DATA, ONLY: FUEL_DROPLET_MLR
REAL(EB), POINTER, DIMENSION(:,:,:) :: DROP_DEN=>NULL(),DROP_RAD=>NULL(),DROP_TMP=>NULL(),MVAP_TOT=>NULL(),DROP_AREA=>NULL()
REAL(EB), POINTER, DIMENSION(:) :: FILM_THICKNESS=>NULL()
REAL(EB) :: R_DROP,NUSSELT,K_AIR,H_V,H_V_REF, H_L,&
            RVC,WGT,Q_CON_GAS,Q_CON_WALL,Q_RAD,H_HEAT,H_MASS,SH_FAC_GAS,SH_FAC_WALL,NU_FAC_GAS,NU_FAC_WALL, &
            T,PR_AIR,M_VAP,M_VAP_MAX,XI,YJ,ZK,RDT,MU_AIR,H_SOLID,Q_DOT_RAD,DEN_ADD,AREA_ADD, &
            Y_DROP,Y_GAS,LENGTH,U2,V2,W2,VEL,DENOM,DZ_DTMP_DROP,TMP_DROP_NEW,TMP_WALL,H_WALL, &
            SC_AIR,D_AIR,DHOR,SHERWOOD,X_DROP,M_DROP,RHO_G,MW_RATIO,MW_DROP,FTPR,&
            C_DROP,M_GAS,A_DROP,TMP_G,TMP_DROP,TMP_MELT,TMP_BOIL,MINIMUM_FILM_THICKNESS,RE_L,OMRAF,Q_FRAC,Q_TOT,DT_SUBSTEP, &
            CP,H_NEW,ZZ_GET(0:N_TRACKED_SPECIES),ZZ_GET2(0:N_TRACKED_SPECIES), &
            M_GAS_NEW,MW_GAS,CP2,VEL_REL,DELTA_H_G,TMP_G_I,H_G_OLD,&
            H_L_REF,TMP_G_NEW,DT_SUM,DCPDT,TMP_WGT,X_EQUIL,Y_EQUIL,Y_ALL(1:N_SPECIES)
INTEGER :: I,II,JJ,KK,IW,N_PC,EVAP_INDEX,N_SUBSTEPS,ITMP,IBC,ITCOUNT,NS
INTEGER, INTENT(IN) :: NM
LOGICAL :: TEMPITER
TYPE (DROPLET_TYPE), POINTER :: DR=>NULL()
TYPE (PARTICLE_CLASS_TYPE), POINTER :: PC=>NULL()
TYPE (SURFACE_TYPE), POINTER :: SF=>NULL()
TYPE (SPECIES_TYPE), POINTER :: SS=>NULL()

! Initializations

RDT    = 1._EB/DT
OMRAF  = 1._EB - RUN_AVG_FAC
FUEL_DROPLET_MLR(NM) = 0._EB

! Rough estimates

MINIMUM_FILM_THICKNESS = 1.E-5_EB   ! Minimum thickness of liquid film on the surface (m)
H_SOLID                = 300._EB    ! Heat transfer coefficient from solid surface to drop (W/m2/K)

! Empirical coefficients

D_AIR                  = 2.6E-5_EB  ! Water Vapor - Air binary diffusion (m2/s at 25 C, Incropera & DeWitt, Table A.8) 
SC_AIR                 = 0.6_EB     ! NU_AIR/D_AIR (Incropera & DeWitt, Chap 7, External Flow)
PR_AIR                 = 0.7_EB     
SH_FAC_GAS             = 0.6_EB*SC_AIR**ONTH
NU_FAC_GAS             = 0.6_EB*PR_AIR**ONTH        
SH_FAC_WALL            = 0.037_EB*SC_AIR**ONTH
NU_FAC_WALL            = 0.037_EB*PR_AIR**ONTH        

! Working arrays

IF (N_EVAP_INDICES>0) THEN
   WCPUA  = RUN_AVG_FAC*WCPUA
   WMPUA  = RUN_AVG_FAC*WMPUA
   DROP_AREA => WORK1
   DROP_DEN => WORK4
   DROP_RAD => WORK5
   DROP_TMP => WORK6
   MVAP_TOT => WORK7
   MVAP_TOT = 0._EB
ENDIF

! Loop over all types of evaporative species
EVAP_INDEX_LOOP: DO EVAP_INDEX = 1,N_EVAP_INDICES

   FILM_THICKNESS => WALL_WORK2
   FILM_THICKNESS =  0._EB

   ! Loop over all particle/droplet classes that have the given evaporative index
   ! First loop is for evaporation and energy transfer
   
   PART_CLASS_LOOP: DO N_PC = 1,N_PART

      PC => PARTICLE_CLASS(N_PC)

      ! If the particle class is associated with a particular SURF, just get its temperature and cycle
      
      IF (PC%SURF_INDEX>0) THEN
         SF => SURFACE(PC%SURF_INDEX)
         DO I=1,NLP
            DR=>DROPLET(I)
            IF (DR%CLASS==N_PC) THEN
               IW = DR%WALL_INDEX
               DR%TMP = TMP_F(IW)
               IF (SF%SHRINK) THEN
                  DR%R = SUM(WALL(IW)%LAYER_THICKNESS)
               ELSE
                  IF (SF%THERMALLY_THICK) THEN
                     DR%R = MAXVAL(SF%X_S)
                  ELSE
                     DR%R = SF%RADIUS
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
         CYCLE PART_CLASS_LOOP
      ENDIF

      ! Check to see if the particles/droplets evaporate
      
      IF (.NOT.PC%EVAPORATE) CYCLE PART_CLASS_LOOP      
      IF (PC%EVAP_INDEX/=EVAP_INDEX) CYCLE PART_CLASS_LOOP
      IF (PC%TREE) CYCLE PART_CLASS_LOOP

      ! Initialize quantities common to the PARTICLE_CLASS 
      SS => SPECIES(PC%Y_INDEX)
      FTPR     = PC%FTPR
      TMP_MELT = SS%TMP_MELT
      TMP_BOIL = SS%TMP_V
      MW_DROP  = SS%MW
      H_V_REF  = SS%H_V(NINT(SS%H_V_REFERENCE_TEMPERATURE))
      H_L_REF  = SS%C_P_L_BAR(NINT(TMP_MELT))*TMP_MELT

      ! Loop through all droplets in the class and determine the depth of the liquid film on each surface cell

      FILM_SUMMING_LOOP: DO I=1,NLP
         DR => DROPLET(I)
         IF (DR%IOR==0)        CYCLE FILM_SUMMING_LOOP
         IF (DR%WALL_INDEX==0) CYCLE FILM_SUMMING_LOOP
         IF (DR%CLASS /= N_PC) CYCLE FILM_SUMMING_LOOP
         IF (DR%R<=0._EB)      CYCLE FILM_SUMMING_LOOP
         IW = DR%WALL_INDEX
         FILM_THICKNESS(IW) = FILM_THICKNESS(IW) + DR%PWT*DR%R**3/AW(IW)
      ENDDO FILM_SUMMING_LOOP
      FILM_THICKNESS = FILM_THICKNESS*FTPR/PC%DENSITY

      FILM_THICKNESS = MAX(MINIMUM_FILM_THICKNESS,FILM_THICKNESS) 

      ! Loop through all droplets within the class and determine mass/energy transfer

      DROPLET_LOOP: DO I=1,NLP

         DR => DROPLET(I)
         IF (DR%CLASS /= N_PC) CYCLE DROPLET_LOOP
         IF (DR%R<=0._EB)      CYCLE DROPLET_LOOP

         ! Determine the current coordinates of the particle

         XI = CELLSI(FLOOR((DR%X-XS)*RDXINT))
         YJ = CELLSJ(FLOOR((DR%Y-YS)*RDYINT))
         ZK = CELLSK(FLOOR((DR%Z-ZS)*RDZINT))
         II  = FLOOR(XI+1._EB)
         JJ  = FLOOR(YJ+1._EB)
         KK  = FLOOR(ZK+1._EB)
         RVC = RDX(II)*RRN(II)*RDY(JJ)*RDZ(KK)
         
         ! Determine how many sub-time step iterations are needed and then iterate over the time step.
         ! This is not fully functional. Keep as a placeholder for now.

         N_SUBSTEPS = 1
         DT_SUBSTEP = DT/REAL(N_SUBSTEPS,EB) 
         DT_SUM = 0._EB

         TIME_ITERATION_LOOP: DO WHILE (DT_SUM < DT)
            ZZ_GET = 0._EB
            IF (N_TRACKED_SPECIES>0) THEN
               ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(II,JJ,KK,1:N_TRACKED_SPECIES)
               CALL GET_MASS_FRACTION_ALL(ZZ_GET,Y_ALL)
               IF (Y_ALL(PC%Y_INDEX) >=1._EB) Y_ALL = SPECIES_MIXTURE(0)%MASS_FRACTION
               MW_GAS = 0._EB
               DO NS=1,N_SPECIES
                  IF (NS==PC%Y_INDEX) CYCLE
                  MW_GAS = MW_GAS + Y_ALL(NS)/SPECIES(NS)%MW
               ENDDO
               MW_GAS = (1._EB-Y_ALL(PC%Y_INDEX))/MW_GAS
            ENDIF
            MW_RATIO = MW_GAS/MW_DROP
            ! Initialize droplet thermophysical data

            R_DROP   = DR%R
            M_DROP   = FTPR*R_DROP**3
            TMP_DROP = DR%TMP
            ITMP     = INT(TMP_DROP)
            TMP_WGT  = TMP_DROP - AINT(TMP_DROP)
            H_V      = SS%H_V(ITMP)+TMP_WGT*(SS%H_V(ITMP+1)-SS%H_V(ITMP))
            C_DROP   = SS%C_P_L(ITMP)+TMP_WGT*(SS%C_P_L(ITMP+1)-SS%C_P_L(ITMP))
            H_L      = (SS%C_P_L_BAR(ITMP)+TMP_WGT*(SS%C_P_L_BAR(ITMP+1)-SS%C_P_L_BAR(ITMP)))*TMP_DROP-H_L_REF
            WGT      = DR%PWT
            DHOR     = H_V*MW_DROP/R0   

            ! Gas conditions

            TMP_G  = TMP(II,JJ,KK)
            RHO_G  = RHO(II,JJ,KK)
            MU_AIR = MU_Z(MIN(5000,NINT(TMP_G)),0)*SPECIES_MIXTURE(0)%MW
            M_GAS  = RHO_G/RVC        
            M_VAP_MAX = (0.33_EB * M_GAS - MVAP_TOT(II,JJ,KK)) / WGT ! limit to avoid diveregence errors
            K_AIR  = CPOPR*MU_AIR
            IF (PC%Y_INDEX>=0) THEN
               CALL GET_MASS_FRACTION(ZZ_GET,PC%Y_INDEX,Y_GAS)
            ELSE
               Y_GAS = 0._EB
            ENDIF
            U2 = 0.5_EB*(U(II,JJ,KK)+U(II-1,JJ,KK))
            V2 = 0.5_EB*(V(II,JJ,KK)+V(II,JJ-1,KK))
            W2 = 0.5_EB*(W(II,JJ,KK)+W(II,JJ,KK-1))
            VEL_REL = SQRT((U2-DR%U)**2+(V2-DR%V)**2+(W2-DR%W)**2)
            ! Set variables for heat transfer on solid

            SOLID_OR_GAS_PHASE: IF (DR%IOR/=0 .AND. DR%WALL_INDEX>0) THEN

               IW   = DR%WALL_INDEX
               A_DROP = M_DROP/(FILM_THICKNESS(IW)*PC%DENSITY)
               TMP_WALL = TMP_F(IW) 
               SELECT CASE(ABS(DR%IOR))
                  CASE(1)
                     VEL = SQRT(V2**2+W2**2)
                  CASE(2)
                     VEL = SQRT(U2**2+W2**2)
                  CASE(3)
                     VEL = SQRT(U2**2+V2**2)
               END SELECT
               LENGTH   = 1._EB
               RE_L     = MAX(5.E5_EB,RHO_G*VEL*LENGTH/MU_AIR)
               NUSSELT  = NU_FAC_WALL*RE_L**0.8_EB
               SHERWOOD = SH_FAC_WALL*RE_L**0.8_EB
               H_HEAT   = NUSSELT*K_AIR/LENGTH
               IF (PC%EVAPORATE) THEN
                  H_MASS = SHERWOOD*D_AIR/LENGTH
               ELSE
                  H_MASS = 0._EB
               ENDIF
               H_WALL    = H_SOLID
               Q_DOT_RAD = A_DROP*QRADIN(IW)

            ELSE SOLID_OR_GAS_PHASE

               A_DROP   = 4._EB*PI*R_DROP**2
               NUSSELT  = 2._EB + NU_FAC_GAS*SQRT(DR%RE)
               SHERWOOD = 2._EB + SH_FAC_GAS*SQRT(DR%RE)
               H_HEAT   = NUSSELT *K_AIR/(2._EB*R_DROP)
               IF (PC%EVAPORATE) THEN
                  H_MASS = SHERWOOD*D_AIR/(2._EB*R_DROP)
               ELSE
                  H_MASS = 0._EB
               ENDIF
               H_WALL   = 0._EB
               TMP_WALL = TMPA
               IF (AVG_DROP_DEN(II,JJ,KK,EVAP_INDEX )>0._EB) THEN
                  Q_DOT_RAD = (QR_W(II,JJ,KK)/SUM(AVG_DROP_AREA(II,JJ,KK,:)))*(A_DROP/4._EB)
               ELSE
                  Q_DOT_RAD = 0._EB
               ENDIF

            ENDIF SOLID_OR_GAS_PHASE
         
            ! Compute equilibrium droplet vapor mass fraction, Y_DROP, and its derivative w.r.t. droplet temperature
    
            IF (PC%EVAPORATE) THEN
               X_DROP  = MIN(1._EB,EXP(DHOR*(1._EB/TMP_BOIL-1._EB/TMP_DROP)))
               Y_DROP  = X_DROP/(MW_RATIO + (1._EB-MW_RATIO)*X_DROP)
               IF (TMP_DROP < TMP_BOIL) THEN
                  DZ_DTMP_DROP = (MW_RATIO/(X_DROP*(1._EB-MW_RATIO)+MW_RATIO)**2)*DHOR*X_DROP/TMP_DROP**2
               ELSE
                  DZ_DTMP_DROP = 0._EB
               ENDIF
               IF (Y_DROP<=Y_GAS) H_MASS = 0._EB
            ELSE
               DZ_DTMP_DROP = 0._EB
               Y_DROP       = 0._EB
            ENDIF
            ! Update the droplet temperature semi_implicitly
    
            DENOM = 1._EB + (H_HEAT + H_WALL + H_MASS*RHO_G*H_V*DZ_DTMP_DROP)*DT_SUBSTEP*A_DROP/(2._EB*M_DROP*C_DROP) 
            TMP_DROP_NEW = ( TMP_DROP + DT_SUBSTEP*( Q_DOT_RAD + &
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
            
            ! Evaporate completely small droplets

            IF (PC%EVAPORATE .AND. R_DROP<0.5_EB*PC%MINIMUM_DIAMETER) THEN
               M_VAP  = M_DROP/N_SUBSTEPS
               IF (Q_TOT>0._EB) THEN
                  Q_FRAC = M_VAP*H_V/Q_TOT 
                  Q_CON_GAS  = Q_CON_GAS*Q_FRAC
                  Q_CON_WALL = Q_CON_WALL*Q_FRAC
                  Q_RAD      = Q_RAD*Q_FRAC
                  Q_TOT  = Q_RAD+Q_CON_GAS+Q_CON_WALL
               ENDIF
            ENDIF
            IF (M_VAP < M_DROP) TMP_DROP_NEW = TMP_DROP + (Q_TOT - M_VAP * H_V)/(C_DROP * (M_DROP - M_VAP))

            ! If the droplet temperature drops below its freezing point, just reset it

            IF (PC%EVAPORATE .AND. TMP_DROP_NEW<TMP_MELT) TMP_DROP_NEW = TMP_MELT

            ! If the droplet temperature reaches boiling, use only enough energy from gas to vaporize liquid

            IF (PC%EVAPORATE .AND. TMP_DROP_NEW>=TMP_BOIL) THEN  
               M_VAP  = MIN(M_VAP_MAX,M_DROP,M_VAP + (TMP_DROP_NEW - TMP_BOIL)*C_DROP*M_DROP/H_V)
               TMP_DROP_NEW = TMP_BOIL
               IF (Q_TOT>0._EB) THEN
                  Q_FRAC = M_VAP*H_V/Q_TOT
                  Q_CON_GAS  = Q_CON_GAS*Q_FRAC
                  Q_CON_WALL = Q_CON_WALL*Q_FRAC
                  Q_RAD      = Q_RAD*Q_FRAC
                  Q_TOT  = Q_RAD+Q_CON_GAS+Q_CON_WALL
               ENDIF
            ENDIF
            M_DROP = M_DROP - M_VAP
        
            ! Compute surface cooling and density

            IF (DR%IOR/=0 .AND. DR%WALL_INDEX>0) THEN
               WCPUA(IW,EVAP_INDEX) = WCPUA(IW,EVAP_INDEX) + OMRAF*WGT*(Q_RAD+Q_CON_WALL)*RAW(IW)/DT_SUBSTEP
               WMPUA(IW,EVAP_INDEX) = WMPUA(IW,EVAP_INDEX) + OMRAF*WGT*M_DROP*RAW(IW)/REAL(N_SUBSTEPS,EB) 
            ENDIF

            ! Add fuel evaporation rate to running counter and adjust mass of evaporated fuel to account for different 
            ! Heat of Combustion between fuel droplet and gas

            IF (N_REACTIONS>0) THEN
               IF (PC%Z_INDEX==REACTION(1)%FUEL_SMIX_INDEX) THEN
                  FUEL_DROPLET_MLR(NM) = FUEL_DROPLET_MLR(NM) + WGT*M_VAP/DT_SUBSTEP
                  M_VAP = PC%ADJUST_EVAPORATION*M_VAP
               ENDIF
            ENDIF

            ! Update gas temperature and determine new subtimestep

            DELTA_H_G = -Q_CON_GAS+M_VAP*(H_V )!- 2011649.919_EB)
            CALL GET_AVERAGE_SPECIFIC_HEAT(ZZ_GET,CP,TMP_G)            
            H_G_OLD = M_GAS*CP*TMP_G
            M_GAS_NEW = M_GAS + WGT*M_VAP
            TMP_G_NEW = TMP_G
            H_NEW = H_G_OLD + DELTA_H_G*WGT!(DELTA_H_G * M_VAP - Q_CON_GAS) * WGT 
            IF (H_NEW > 0._EB) THEN
               ZZ_GET2 = ZZ_GET * M_GAS/M_GAS_NEW               
               IF (PC%Z_INDEX>0) ZZ_GET2(PC%Z_INDEX) = ZZ_GET2(PC%Z_INDEX) + WGT*M_VAP/M_GAS_NEW
               TMP_G_I = TMP_G
               TEMPITER = .TRUE.
               ITCOUNT = 0
               ITERATE_TEMP: DO WHILE (TEMPITER)
                  TEMPITER=.FALSE.
                  CALL GET_AVERAGE_SPECIFIC_HEAT(ZZ_GET2,CP2,TMP_G_I)
                  IF (TMP_G_I < 4999._EB) THEN
                     CALL GET_AVERAGE_SPECIFIC_HEAT(ZZ_GET2,CP,TMP_G_I-1._EB)
                     DCPDT = CP-CP2                     
                  ELSE
                     CALL GET_AVERAGE_SPECIFIC_HEAT(ZZ_GET2,CP,TMP_G_I-1._EB)                  
                     DCPDT = CP2-CP
                  ENDIF
                  ! Compute approximation of d(cp)/dT                  
                  TMP_G_I = TMP_G_I+(H_NEW-CP2*TMP_G_I*M_GAS_NEW)/(M_GAS_NEW*(CP2+TMP_G_I*DCPDT))
                  TMP_G_I = MAX(TMPMIN,TMP_G_I)
                  ITCOUNT = ITCOUNT + 1
                  IF (ABS(TMP_G_NEW-TMP_G_I) > 0.5_EB) TEMPITER = .TRUE.
                  IF (ITCOUNT > 10) THEN
                     TMP_G_NEW = 0.5_EB*(TMP_G_I + TMP_G_NEW)
                     EXIT ITERATE_TEMP
                  ENDIF
                  TMP_G_NEW = TMP_G_I
               ENDDO ITERATE_TEMP
            ELSE
               DT_SUBSTEP = DT_SUBSTEP * 0.5_EB
               N_SUBSTEPS = NINT(DT/DT_SUBSTEP)
               CYCLE TIME_ITERATION_LOOP
            ENDIF
            ITMP     = INT(TMP_DROP_NEW)
            TMP_WGT  = TMP_DROP_NEW - AINT(TMP_DROP_NEW)
            H_V      = SS%H_V(ITMP)+TMP_WGT*(SS%H_V(ITMP+1)-SS%H_V(ITMP))
            DHOR     = H_V*MW_DROP/R0 
            X_EQUIL  = MIN(1._EB,EXP(DHOR*(1._EB/TMP_BOIL-1._EB/MIN(TMP_DROP_NEW,TMP_BOIL))))
            Y_EQUIL = X_EQUIL/(MW_RATIO + (1._EB-MW_RATIO)*X_EQUIL)
            !Limit supersaturation
            IF (Y_GAS < Y_EQUIL) THEN
               CALL GET_MASS_FRACTION(ZZ_GET2,PC%Y_INDEX,Y_GAS)
               IF (Y_GAS/Y_EQUIL > 1.02_EB) THEN
                  DT_SUBSTEP = DT_SUBSTEP * 0.5_EB            
                  N_SUBSTEPS = NINT(DT/DT_SUBSTEP)
                  IF (DT_SUBSTEP <= 0.00001_EB*DT) THEN
                     CALL SHUTDOWN('Numerical instability in particle energy transport')
                  ENDIF
                  CYCLE TIME_ITERATION_LOOP
               ENDIF
            ENDIF
            
            !Limit gas temperature change
            IF (ABS(TMP_G_NEW/TMP_G - 1._EB) > 0.05_EB) THEN
               DT_SUBSTEP = DT_SUBSTEP * 0.5_EB            
               N_SUBSTEPS = NINT(DT/DT_SUBSTEP)
               IF (DT_SUBSTEP <= 0.00001_EB*DT) THEN
                  CALL SHUTDOWN('Numerical instability in particle energy transport')
               ENDIF
               CYCLE TIME_ITERATION_LOOP
            ENDIF

            ! Update gas cell density, temperature, and mass fractions

            RHO(II,JJ,KK) = M_GAS_NEW*RVC

            ZZ(II,JJ,KK,1:N_TRACKED_SPECIES) = ZZ_GET2(1:N_TRACKED_SPECIES)
            CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET2,RSUM(II,JJ,KK))
            TMP(II,JJ,KK) = MIN(TMPMAX,MAX(TMPMIN,TMP_G_NEW))

            ! Compute change in enthalpy between gas and liquid

            CALL GET_SPECIFIC_HEAT(ZZ_GET,CP,TMP_G)
            H_G_OLD = CP * TMP_G * M_GAS
            ZZ_GET = 0._EB
            IF (PC%Z_INDEX>0) ZZ_GET(PC%Z_INDEX) = 1._EB
            DELTA_H_G = H_L + H_V
            CALL GET_AVERAGE_SPECIFIC_HEAT(ZZ_GET,CP2,TMP_G)
            DELTA_H_G = DELTA_H_G - CP2 * TMP_G     

            ! Compute contribution to the divergence

            D_LAGRANGIAN(II,JJ,KK) =  D_LAGRANGIAN(II,JJ,KK) + (MW_RATIO *M_VAP /M_GAS + &
                                      (M_VAP*DELTA_H_G - Q_CON_GAS)/H_G_OLD) * WGT / DT_SUBSTEP 

            ! Keep track of total mass evaporated in cell

            MVAP_TOT(II,JJ,KK) = MVAP_TOT(II,JJ,KK) + WGT*M_VAP
            ! Update droplet quantities

            DR%R   = (M_DROP/FTPR)**ONTH
            DR%TMP = TMP_DROP_NEW

            ! Get out of the loop if the droplet has evaporated completely

            IF (DR%R<=0._EB) CYCLE DROPLET_LOOP
            DT_SUM = DT_SUM + DT_SUBSTEP
            DT_SUBSTEP = MIN(DT-DT_SUM,DT_SUBSTEP * 1.5_EB)

         ENDDO TIME_ITERATION_LOOP
         
      ENDDO DROPLET_LOOP
      
   ENDDO PART_CLASS_LOOP

   ! Second loop is for summing the part quantities

   DROP_DEN = 0._EB
   DROP_TMP = 0._EB
   DROP_RAD = 0._EB
   DROP_AREA = 0._EB

   PART_CLASS_SUM_LOOP: DO N_PC = 1,N_PART

      PC => PARTICLE_CLASS(N_PC)
      IF (PC%EVAP_INDEX/=EVAP_INDEX) CYCLE PART_CLASS_SUM_LOOP

      DROPLET_LOOP_2: DO I=1,NLP

         DR => DROPLET(I)
         IF (DR%CLASS /= N_PC) CYCLE DROPLET_LOOP_2
         IF (DR%R<=0._EB)      CYCLE DROPLET_LOOP_2

         ! Determine the current coordinates of the particle

         XI = CELLSI(FLOOR((DR%X-XS)*RDXINT))
         YJ = CELLSJ(FLOOR((DR%Y-YS)*RDYINT))
         ZK = CELLSK(FLOOR((DR%Z-ZS)*RDZINT))
         II  = FLOOR(XI+1._EB)
         JJ  = FLOOR(YJ+1._EB)
         KK  = FLOOR(ZK+1._EB)
         RVC = RDX(II)*RRN(II)*RDY(JJ)*RDZ(KK)

         ! Determine the mass of the droplet/particle, depending on whether the particle has a distinct SURFace type.

         IBC = PC%SURF_INDEX
         IF (IBC<1) THEN
            M_DROP = PC%FTPR*DR%R**3
         ELSE
            SF => SURFACE(IBC)
            IW = DR%WALL_INDEX
            SELECT CASE(SF%GEOMETRY)
               CASE(SURF_CARTESIAN)
                  M_DROP = 2._EB*SF%LENGTH*SF%WIDTH*DR%R*SURFACE_DENSITY(NM,IW,1)  ! The 1 indicates kg/m3, not kg/m2
               CASE(SURF_CYLINDRICAL)
                  M_DROP = PI*DR%R**2*SF%LENGTH*SURFACE_DENSITY(NM,IW,1)
               CASE(SURF_SPHERICAL)
                  M_DROP = FOTH*PI*DR%R**3*SURFACE_DENSITY(NM,IW,1)
            END SELECT
         ENDIF

         ! Assign particle or droplet mass to the grid cell if the particle/droplet not on a surface

         IF (DR%IOR==0) THEN
            DEN_ADD  =    DR%PWT*M_DROP *RVC
            AREA_ADD = PI*DR%PWT*DR%R**2*RVC
            DROP_DEN(II,JJ,KK)  = DROP_DEN(II,JJ,KK)  + DEN_ADD
            DROP_TMP(II,JJ,KK)  = DROP_TMP(II,JJ,KK)  + DEN_ADD*DR%TMP
            DROP_RAD(II,JJ,KK)  = DROP_RAD(II,JJ,KK)  + AREA_ADD*DR%R
            DROP_AREA(II,JJ,KK) = DROP_AREA(II,JJ,KK) + AREA_ADD
         ENDIF

         ! Compute surface density

         IF (DR%IOR/=0 .AND. DR%WALL_INDEX>0) THEN
            IW = DR%WALL_INDEX
            WMPUA(IW,EVAP_INDEX) = WMPUA(IW,EVAP_INDEX) + OMRAF*DR%PWT*M_DROP*RAW(IW)
         ENDIF

      ENDDO DROPLET_LOOP_2

   ENDDO PART_CLASS_SUM_LOOP

  ! Compute cumulative quantities for droplet "clouds"

   DROP_RAD = DROP_RAD/(DROP_AREA+TINY(1._EB))
   DROP_TMP = DROP_TMP/(DROP_DEN +TINY(1._EB))
   AVG_DROP_RAD(:,:,:,EVAP_INDEX ) = DROP_RAD
   AVG_DROP_TMP(:,:,:,EVAP_INDEX ) = RUN_AVG_FAC*AVG_DROP_TMP(:,:,:,EVAP_INDEX ) + OMRAF*DROP_TMP
   AVG_DROP_TMP(:,:,:,EVAP_INDEX ) = MAX(TMPM,AVG_DROP_TMP(:,:,:,EVAP_INDEX ))
   AVG_DROP_DEN(:,:,:,EVAP_INDEX ) = RUN_AVG_FAC*AVG_DROP_DEN(:,:,:,EVAP_INDEX ) + OMRAF*DROP_DEN
   AVG_DROP_AREA(:,:,:,EVAP_INDEX) = RUN_AVG_FAC*AVG_DROP_AREA(:,:,:,EVAP_INDEX) + OMRAF*DROP_AREA
   WHERE (AVG_DROP_DEN(:,:,:,EVAP_INDEX )<0.0001_EB .AND. ABS(DROP_DEN)<ZERO_P) AVG_DROP_DEN(:,:,:,EVAP_INDEX ) = 0.0_EB

ENDDO EVAP_INDEX_LOOP

! Remove droplets that have completely evaporated
 
CALL REMOVE_DROPLETS(T,NM)
 
END SUBROUTINE PARTICLE_MASS_ENERGY_TRANSFER



SUBROUTINE PARTICLE_MOMENTUM_TRANSFER

! Add droplet momentum as a force term in momentum equation

USE COMP_FUNCTIONS, ONLY : SECOND
REAL(EB), POINTER, DIMENSION(:,:,:) :: FVXS=>NULL(),FVYS=>NULL(),FVZS=>NULL()
REAL(EB) :: XI,YJ,ZK
INTEGER :: II,JJ,KK,IIX,JJY,KKZ,I,J,K,IC,IW
TYPE (DROPLET_TYPE), POINTER :: DR=>NULL()
TYPE (PARTICLE_CLASS_TYPE), POINTER :: PC=>NULL()

FVXS  => WORK1
FVYS  => WORK2
FVZS  => WORK3

FVXS  = 0._EB
FVYS  = 0._EB
FVZS  = 0._EB

SUM_MOMENTUM_LOOP: DO I=1,NLP
   DR=>DROPLET(I)
   PC=>PARTICLE_CLASS(DR%CLASS)
   IF (DR%IOR/=0)   CYCLE SUM_MOMENTUM_LOOP
   IF (PC%MASSLESS) CYCLE SUM_MOMENTUM_LOOP
   XI  = CELLSI(FLOOR((DR%X-XS)*RDXINT))
   YJ  = CELLSJ(FLOOR((DR%Y-YS)*RDYINT))
   ZK  = CELLSK(FLOOR((DR%Z-ZS)*RDZINT))
   II  = FLOOR(XI+1._EB)
   JJ  = FLOOR(YJ+1._EB)
   KK  = FLOOR(ZK+1._EB)
   IF (SOLID(CELL_INDEX(II,JJ,KK))) CYCLE SUM_MOMENTUM_LOOP
   IIX = FLOOR(XI+.5_EB)
   JJY = FLOOR(YJ+.5_EB)
   KKZ = FLOOR(ZK+.5_EB)
   IC = CELL_INDEX(IIX,JJ,KK)
   IW = WALL_INDEX(IC,1)
   IF (BOUNDARY_TYPE(IW)==NULL_BOUNDARY) THEN
      FVXS(IIX,JJ,KK) = FVXS(IIX,JJ,KK) - DR%A_X
   ENDIF
   IC = CELL_INDEX(II,JJY,KK)
   IW = WALL_INDEX(IC,2) 
   IF (BOUNDARY_TYPE(IW)==NULL_BOUNDARY) THEN
      FVYS(II,JJY,KK) = FVYS(II,JJY,KK) - DR%A_Y
   ENDIF
   IC = CELL_INDEX(II,JJ,KKZ)
   IW = WALL_INDEX(IC,3) 
   IF (BOUNDARY_TYPE(IW)==NULL_BOUNDARY) THEN
      FVZS(II,JJ,KKZ) = FVZS(II,JJ,KKZ) - DR%A_Z
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
 
END SUBROUTINE PARTICLE_MOMENTUM_TRANSFER
 
 
 
SUBROUTINE REMOVE_DROPLETS(T,NM)

! Remove droplets that do not lie in any mesh
 
USE MEMORY_FUNCTIONS, ONLY : RE_ALLOCATE_DROPLETS 
INTEGER :: IKILL,I,NM,IPC
REAL(EB) :: T,MAX_DROP_LIFE
TYPE (DROPLET_TYPE), POINTER :: DR=>NULL()
TYPE (PARTICLE_CLASS_TYPE), POINTER :: PC=>NULL()

IF (NLP + INSERT_COUNT > MAXIMUM_DROPLETS) THEN
   MAX_DROP_LIFE = REAL(MAXIMUM_DROPLETS,EB) / INSERT_RATE
ELSE
   MAX_DROP_LIFE = 1.E6
ENDIF
IKILL = 0
DROP_LOOP: DO I=1,NLP
   WEED_LOOP: DO
      DR=>DROPLET(I) 
      IPC=DR%CLASS 
      PC=>PARTICLE_CLASS(IPC)
      IF (I>NLP-IKILL) EXIT DROP_LOOP
      IF (DR%R<=PC%KILL_RADIUS .AND. .NOT.PC%MASSLESS) THEN 
         CALL REPLACE
         CYCLE WEED_LOOP
      ENDIF
      IF ((T-DR%T>PC%LIFETIME) .OR. (T-DR%T>MAX_DROP_LIFE .AND. IKILL <= INSERT_COUNT)) THEN
         CALL REPLACE
         CYCLE WEED_LOOP
      ENDIF
      IF (DR%X>XS .AND. DR%X<XF .AND.DR%Y>YS .AND. DR%Y<YF .AND. DR%Z>ZS .AND. DR%Z<ZF)  CYCLE DROP_LOOP
      CALL REPLACE
   ENDDO WEED_LOOP
ENDDO DROP_LOOP
 
NLP = NLP - IKILL
 
CONTAINS
 
SUBROUTINE REPLACE
 
INTEGER :: OM,NOM
TYPE (MESH_TYPE), POINTER :: M=>NULL()
TYPE (OMESH_TYPE), POINTER :: M2=>NULL()

NOM = 0
SEARCH_LOOP: DO OM=1,NMESHES
   IF (MESHES(NM)%OMESH(OM)%NIC_S==0) CYCLE SEARCH_LOOP
   IF (EVACUATION_ONLY(OM)) CYCLE SEARCH_LOOP
   M=>MESHES(OM)
   IF (DR%X>M%XS .AND. DR%X<M%XF .AND.  DR%Y>M%YS .AND. DR%Y<M%YF .AND.  DR%Z>M%ZS .AND. DR%Z<M%ZF) THEN
      NOM = OM
      EXIT SEARCH_LOOP
   ENDIF
ENDDO SEARCH_LOOP
 
IF (NOM/=0) THEN
   M2=>MESHES(NM)%OMESH(NOM)
   M2%N_DROP_ORPHANS = M2%N_DROP_ORPHANS + 1
   IF (M2%N_DROP_ORPHANS>M2%N_DROP_ORPHANS_DIM) CALL RE_ALLOCATE_DROPLETS(2,NM,NOM,1000)
   M2%DROPLET(M2%N_DROP_ORPHANS) = DROPLET(I)
ENDIF
 
DROPLET(I) = DROPLET(NLP-IKILL)
IKILL = IKILL + 1
 
END SUBROUTINE REPLACE
 
END SUBROUTINE REMOVE_DROPLETS


 
SUBROUTINE GET_REV_part(MODULE_REV,MODULE_DATE)
INTEGER,INTENT(INOUT) :: MODULE_REV
CHARACTER(255),INTENT(INOUT) :: MODULE_DATE

WRITE(MODULE_DATE,'(A)') partrev(INDEX(partrev,':')+1:LEN_TRIM(partrev)-2)
READ (MODULE_DATE,'(I5)') MODULE_REV
WRITE(MODULE_DATE,'(A)') partdate

END SUBROUTINE GET_REV_part
END MODULE PART
