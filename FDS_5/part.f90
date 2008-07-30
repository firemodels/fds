MODULE PART
 
USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_POINTERS
USE TRAN
USE TYPES, ONLY: DROPLET_TYPE, PARTICLE_CLASS_TYPE, PARTICLE_CLASS, WALL_TYPE
IMPLICIT NONE
PRIVATE
PUBLIC INSERT_DROPLETS_AND_PARTICLES,UPDATE_PARTICLES,REMOVE_DROPLETS, &
       INITIALIZE_DROPLETS,GET_REV_part
TYPE (DROPLET_TYPE), POINTER :: DR
TYPE (PARTICLE_CLASS_TYPE), POINTER :: PC
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
 
REAL(EB) RN,MASS_SUM,LL,UL,BIN_SIZE,TNOW
REAL(EB) VOL1, VOL2, X1, X2, Y1, Y2, Z1, Z2
INTEGER I,J,II,JJ,KK,IL,IU,STRATUM,IPC
INTEGER, INTENT(IN) :: NM

IF (N_PART==0) RETURN            ! Don't waste time if no particles
IF (EVACUATION_ONLY(NM)) RETURN  ! Don't waste time if an evac mesh
 
TNOW=SECOND()
CALL POINT_TO_MESH(NM)
 
PART_CLASS_LOOP: DO IPC=1,N_PART
 

   PC=>PARTICLE_CLASS(IPC)

   ! If particles/droplets have a size distribution, initialize here
 
   IF_SIZE_DISTRIBUTION: IF (PC%DIAMETER > 0._EB) THEN
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
 
   ! If there is a specified number of initial droplets/particles, initialize
 
   IF (PC%TREE)         CYCLE PART_CLASS_LOOP
   IF (PC%N_INITIAL==0) CYCLE PART_CLASS_LOOP
 
   IF (PC%X1 == 0.0_EB .AND. PC%X2 == 0.0_EB .AND. PC%Y1 == 0.0_EB .AND. PC%Y2 == 0.0_EB .AND.  &
       PC%Z1 == 0.0_EB .AND. PC%Z2 == 0.0_EB) THEN
      X1 = XS 
      X2 = XF
      Y1 = YS 
      Y2 = YF
      Z1 = ZS 
      Z2 = ZF
      VOL2 = (XF - XS) * (YF - YS) * (ZF - ZS)
      VOL1 = VOL2
   ELSE
      IF (PC%X1>XF .OR. PC%X2<XS .OR. PC%Y1>YF .OR. PC%Y2<YS .OR. PC%Z1>ZF .OR. PC%Z2<ZS) CYCLE PART_CLASS_LOOP
      X1 = MAX(PC%X1,XS) 
      X2 = MIN(PC%X2,XF)
      Y1 = MAX(PC%Y1,YS) 
      Y2 = MIN(PC%Y2,YF)
      Z1 = MAX(PC%Z1,ZS) 
      Z2 = MIN(PC%Z2,ZF)
      VOL2 = (PC%X2 - PC%X1) * (PC%Y2 - PC%Y1) * (PC%Z2 - PC%Z1)
      VOL1 = (X2 - X1) * (Y2 - Y1) * (Z2 - Z1)
   ENDIF

   ! Assign properties to the initial droplets/particles
      
   MASS_SUM = 0._EB
   INITIALIZATION_LOOP: DO I=1,PC%N_INITIAL
      NLP = NLP + 1
      PARTICLE_TAG = PARTICLE_TAG + NMESHES
      IF (NLP>NLPDIM) THEN
         CALL RE_ALLOCATE_DROPLETS(1,NM,0,1000)
         DROPLET=>MESHES(NM)%DROPLET
      ENDIF
      DR=>DROPLET(NLP)
      BLOCK_OUT_LOOP:  DO
         CALL RANDOM_NUMBER(RN)
         DR%X = X1 + RN*(X2-X1)
         CALL RANDOM_NUMBER(RN)
         DR%Y = Y1 + RN*(Y2-Y1)
         CALL RANDOM_NUMBER(RN)
         DR%Z = Z1 + RN*(Z2-Z1)
         II = CELLSI(FLOOR((DR%X-XS)*RDXINT)) + 1._EB
         JJ = CELLSJ(FLOOR((DR%Y-YS)*RDYINT)) + 1._EB
         KK = CELLSK(FLOOR((DR%Z-ZS)*RDZINT)) + 1._EB
         IF (.NOT.SOLID(CELL_INDEX(II,JJ,KK))) EXIT BLOCK_OUT_LOOP
      ENDDO BLOCK_OUT_LOOP
      DR%U   = 0._EB                     ! No initial velocity
      DR%V   = 0._EB
      DR%W   = 0._EB
      DR%TMP = PC%TMP_INITIAL            ! Initial temperature
      DR%T   = 0._EB                     ! Insertion time is 0
      DR%R   = 0._EB                     ! Radius is zero unless DIAMETER has been specified
      DR%PWT = 1._EB                     ! Weighting factor is one unless changed below
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
 
      IF (PC%DIAMETER>0._EB) THEN
         IF (PC%MONODISPERSE) THEN
            DR%R   = 0.5_EB*PC%DIAMETER
            DR%PWT = 1._EB
         ELSE
            STRATUM = MOD(NLP-1,NSTRATA) + 1
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
 
   ENDDO INITIALIZATION_LOOP
 
   ! Adjust particle weighting factor PWT so that desired MASS_PER_VOLUME is achieved

   IF (PC%DIAMETER>0._EB) DROPLET(NLP-PC%N_INITIAL+1:NLP)%PWT = &
      DROPLET(NLP-PC%N_INITIAL+1:NLP)%PWT*PC%MASS_PER_VOLUME*VOL1/MASS_SUM

ENDDO PART_CLASS_LOOP

TUSED(8,NM)=TUSED(8,NM)+SECOND()-TNOW
END SUBROUTINE INITIALIZE_DROPLETS
 

SUBROUTINE INSERT_DROPLETS_AND_PARTICLES(T,NM)

! Insert sprinkler droplets and lagrangian particles into the domain every DT_INSERT seconds

USE COMP_FUNCTIONS, ONLY : SECOND 
USE MEMORY_FUNCTIONS, ONLY : RE_ALLOCATE_DROPLETS
USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
USE DEVICE_VARIABLES
REAL(EB), INTENT(IN) :: T
INTEGER, INTENT(IN) :: NM
REAL(EB) :: PHI_RN,RN,FLOW_RATE,THETA_RN,SPHI,CPHI,MASS_SUM, &
            STHETA,CTHETA,PWT0,DROPLET_SPEED,XI,YJ,ZK,SHIFT1,SHIFT2,XTMP,YTMP,ZTMP,VLEN, &
            TRIGT1,TRIGT2,TNOW,SOLID_ANGLE_FLOWRATE,SOLID_ANGLE_TOTAL_FLOWRATE,TSI
REAL(EB), PARAMETER :: VENT_OFFSET=0.1
INTEGER :: I,KS,II,JJ,KK,IC,IL,IU,STRATUM,IPC,DROP_SUM,IIG,JJG,KKG,IW,IBC,IOR
LOGICAL :: INSERT_ANOTHER_BATCH
TYPE (PROPERTY_TYPE), POINTER :: PY
TYPE (TABLES_TYPE), POINTER :: TA
TYPE (DEVICE_TYPE), POINTER :: DV
TYPE (SURFACE_TYPE), POINTER :: SF
 
IF (EVACUATION_ONLY(NM)) RETURN  ! Don't waste time if an evac mesh
IF (N_PART==0) RETURN  ! Don't waste time if no particles
TNOW=SECOND()
CALL POINT_TO_MESH(NM)

OVERALL_INSERT_LOOP: DO  ! Add more than one batch of particles/droplets if FDS time step is large enough
INSERT_ANOTHER_BATCH = .FALSE.
INSERT_RATE = 0._EB
INSERT_COUNT = 0

SPRINKLER_INSERT_LOOP: DO KS=1,N_DEVC  ! Loop over all devices, but look for sprinklers or nozzles
   DV => DEVICE(KS)
   PY => PROPERTY(DV%PROP_INDEX)
   IF (PY%PART_ID == 'null')   CYCLE SPRINKLER_INSERT_LOOP
   IF (DV%MESH/=NM)            CYCLE SPRINKLER_INSERT_LOOP
   IF (.NOT. DV%CURRENT_STATE) CYCLE SPRINKLER_INSERT_LOOP
   IPC = PROPERTY(DV%PROP_INDEX)%PART_INDEX
   PC=>PARTICLE_CLASS(IPC)
   INSERT_RATE = INSERT_RATE + PC%N_INSERT/PC%DT_INSERT
   INSERT_COUNT = INSERT_COUNT + PC%N_INSERT
   IF (DV%T_CHANGE == T) THEN 
      DV%T = T
      CYCLE SPRINKLER_INSERT_LOOP
   ENDIF
   IF (T < PC%INSERT_CLOCK(NM)) CYCLE SPRINKLER_INSERT_LOOP
   IF (T_BEGIN == DV%T_CHANGE .AND. PY%FLOW_RAMP_INDEX>=1) THEN
      TSI = T
   ELSE
      TSI = T - DV%T_CHANGE
   ENDIF
   FLOW_RATE = EVALUATE_RAMP(TSI,PY%FLOW_TAU,PY%FLOW_RAMP_INDEX)*PY%FLOW_RATE*(PC%DENSITY/1000._EB)/60._EB  ! kg/s
   IF (FLOW_RATE <= 0._EB) THEN
      DV%T = T
      CYCLE SPRINKLER_INSERT_LOOP
   ENDIF

   ! Direction initialization stuff

   TRIGT1 = ACOS(-DV%ORIENTATION(3))
   IF (DV%ORIENTATION(2)==0._EB) THEN
      TRIGT2 = ACOS(1._EB)
   ELSE
      TRIGT2 = ACOS(ABS(DV%ORIENTATION(1))/SQRT(DV%ORIENTATION(1)**2+DV%ORIENTATION(2)**2))
   ENDIF
 
   ! Droplet insertion loop
 
   MASS_SUM = 0._EB
   DROP_SUM = 0
   SOLID_ANGLE_TOTAL_FLOWRATE = 0._EB

   DROPLET_INSERT_LOOP: DO I=1,PC%N_INSERT
      IF (NLP+1>MAXIMUM_DROPLETS) EXIT DROPLET_INSERT_LOOP
      NLP = NLP+1
      PARTICLE_TAG = PARTICLE_TAG + NMESHES
      IF (NLP>NLPDIM) THEN
         CALL RE_ALLOCATE_DROPLETS(1,NM,0,1000)
         DROPLET=>MESHES(NM)%DROPLET
      ENDIF
      DR=>DROPLET(NLP)
 
      ! Set droplet propeties
 
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
            THETA_RN = TA%TABLE_DATA(II,1) + RN*(TA%TABLE_DATA(II,2)-TA%TABLE_DATA(II,1))
            CALL RANDOM_NUMBER(RN)
            PHI_RN = TA%TABLE_DATA(II,3) + RN*(TA%TABLE_DATA(II,4)-TA%TABLE_DATA(II,3))
            SOLID_ANGLE_FLOWRATE = TA%TABLE_DATA(II,6)
            SOLID_ANGLE_TOTAL_FLOWRATE = SOLID_ANGLE_TOTAL_FLOWRATE + TA%TABLE_DATA(II,6)
            PY%DROPLET_VELOCITY  = TA%TABLE_DATA(II,5)
         ELSE PICK_PATTERN !Use conical spray
            CALL RANDOM_NUMBER(RN)
            THETA_RN = PY%SPRAY_ANGLE(1) + RN*(PY%SPRAY_ANGLE(2)-PY%SPRAY_ANGLE(1))
            CALL RANDOM_NUMBER(RN)
            PHI_RN = RN*TWOPI
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
         DROPLET_SPEED = PY%DROPLET_VELOCITY
 
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
         DR%R   = 0.5_EB*PC%DIAMETER
         DR%PWT = 1._EB
      ELSE
         STRATUM = MOD(NLP-1,NSTRATA) + 1
         IL = PC%IL_CDF(STRATUM)
         IU = PC%IU_CDF(STRATUM)
         CALL RANDOM_CHOICE(PC%CDF(IL:IU), PC%R_CDF(IL:IU), IU-IL,DR%R)
         DR%PWT = PC%W_CDF(STRATUM)
         IF (2._EB*DR%R > PC%MAXIMUM_DIAMETER) THEN
            DR%PWT = DR%PWT*DR%R**3/(0.5_EB*PC%MAXIMUM_DIAMETER)**3
            DR%R = 0.5_EB*PC%MAXIMUM_DIAMETER
         ENDIF
      ENDIF
 
      ! Sum up mass of liquid being introduced

      IF(PY%SPRAY_PATTERN_INDEX>0) THEN
         MASS_SUM = MASS_SUM + DR%PWT*PC%FTPR*DR%R**3*SOLID_ANGLE_FLOWRATE
      ELSE
         MASS_SUM = MASS_SUM + DR%PWT*PC%FTPR*DR%R**3
      ENDIF
      DROP_SUM = DROP_SUM + 1
   ENDDO DROPLET_INSERT_LOOP
 
   ! Compute weighting factor for the droplets just inserted

   IF (DROP_SUM > 0) THEN
      IF(PY%SPRAY_PATTERN_INDEX>0) THEN
         PWT0 = FLOW_RATE*PC%DT_INSERT/(MASS_SUM*REAL(PC%N_INSERT,EB)/SOLID_ANGLE_TOTAL_FLOWRATE)
      ELSE
         PWT0 = FLOW_RATE*PC%DT_INSERT/MASS_SUM
      ENDIF
      DROPLET(NLP-DROP_SUM+1:NLP)%PWT = DROPLET(NLP-DROP_SUM+1:NLP)%PWT*PWT0
   ENDIF

   ! Indicate that droplets from this device have been inserted at this time T

   DV%T = T 

ENDDO SPRINKLER_INSERT_LOOP
 
! Loop through all boundary cells and insert particles if appropriate
 
WALL_INSERT_LOOP: DO IW=1,NWC
 
   IF (T < TW(IW))                        CYCLE WALL_INSERT_LOOP
   IF (BOUNDARY_TYPE(IW)/=SOLID_BOUNDARY) CYCLE WALL_INSERT_LOOP
   IBC =  IJKW(5,IW)
   SF  => SURFACE(IBC)
   IPC =  SF%PART_INDEX
   IF (IPC < 1)    CYCLE WALL_INSERT_LOOP
   PC  => PARTICLE_CLASS(IPC)
   IF (T < PC%INSERT_CLOCK(NM)) CYCLE WALL_INSERT_LOOP
   INSERT_RATE = INSERT_RATE + NPPCW(IW)/PC%DT_INSERT
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
 
   IOR = IJKW(4,IW)
   MASS_SUM = 0._EB
   PARTICLE_INSERT_LOOP: DO I=1,NPPCW(IW)  ! Loop over all particles for the IW'th cell
 
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
            STRATUM = MOD(NLP-1,NSTRATA) + 1
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
         DROPLET(NLP-NPPCW(IW)+1:NLP)%PWT = &
         DROPLET(NLP-NPPCW(IW)+1:NLP)%PWT*SF%PARTICLE_MASS_FLUX*AREA_ADJUST(IW)*AW(IW)*PC%DT_INSERT/MASS_SUM
      ENDIF
   ENDIF

ENDDO WALL_INSERT_LOOP

! Reset particle/droplet insertion clock
 
DO IPC=1,N_PART
   PC => PARTICLE_CLASS(IPC)
   IF (T >= PC%INSERT_CLOCK(NM)) PC%INSERT_CLOCK(NM) = PC%INSERT_CLOCK(NM) + PC%DT_INSERT
   IF (T >= PC%INSERT_CLOCK(NM)) INSERT_ANOTHER_BATCH = .TRUE.
ENDDO

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

IF (MESHES(NM)%NLP==0) RETURN
IF (EVACUATION_ONLY(NM)) RETURN  ! Don't waste time if an evac mesh
TNOW=SECOND()
CALL POINT_TO_MESH(NM)

IF (CORRECTOR) CALL MOVE_PARTICLES(T,NM)
IF (CORRECTOR) CALL PARTICLE_MASS_ENERGY_TRANSFER(T,NM)
CALL PARTICLE_MOMENTUM_TRANSFER

TUSED(8,NM)=TUSED(8,NM)+SECOND()-TNOW
END SUBROUTINE UPDATE_PARTICLES



SUBROUTINE MOVE_PARTICLES(T,NM)

! Momentum transfer from all particles and droplets
 
USE COMP_FUNCTIONS, ONLY : SECOND  
USE MATH_FUNCTIONS, ONLY : EVALUATE_RAMP, AFILL2
USE PHYSICAL_FUNCTIONS, ONLY : DRAG
 
REAL(EB) :: RHO_G,RVC,RDS,RDC,QREL,SFAC,UREL,VREL,WREL,TMP_G,RN,THETA_RN, &
            RD,T,C_DRAG,XI,YJ,ZK,MU_AIR, &
            DTSP,DTMIN,UBAR,VBAR,WBAR,BFAC,GRVT1,GRVT2,GRVT3,AUREL,AVREL,AWREL,CONST, &
            UVW,DUMMY=0._EB,X_OLD,Y_OLD,Z_OLD,STEP_FRACTION(-3:3),R_NEW,SURFACE_DROPLET_DIAMETER
LOGICAL :: HIT_SOLID
INTEGER :: ICN,I,IIN,JJN,KKN,II,JJ,KK,IIX,JJY,KKZ,IW,N,NITER,IWP1,IWM1,IWP2,IWM2,IWP3,IWM3,IOR_OLD,IC,IOR_FIRST,IML
INTEGER, INTENT(IN) :: NM

SURFACE_DROPLET_DIAMETER = 0.001_EB  ! All droplets adjusted to this size when on solid (m)

GRVT1 = -EVALUATE_RAMP(T,DUMMY,I_RAMP_GX)*GVEC(1) 
GRVT2 = -EVALUATE_RAMP(T,DUMMY,I_RAMP_GY)*GVEC(2) 
GRVT3 = -EVALUATE_RAMP(T,DUMMY,I_RAMP_GZ)*GVEC(3) 

! Loop through all Lagrangian particles and move them one time step

DROPLET_LOOP: DO I=1,NLP  

   DR => DROPLET(I)
   PC => PARTICLE_CLASS(DR%CLASS)

   ! Determine the current coordinates of the particle

   XI = CELLSI(FLOOR((DR%X-XS)*RDXINT))
   YJ = CELLSJ(FLOOR((DR%Y-YS)*RDYINT))
   ZK = CELLSK(FLOOR((DR%Z-ZS)*RDZINT))
   II  = FLOOR(XI+1._EB)
   JJ  = FLOOR(YJ+1._EB)
   KK  = FLOOR(ZK+1._EB)
   IC  = CELL_INDEX(II,JJ,KK)

   ! Throw out particles that are inside a solid obstruction

   IF (SOLID(IC)) THEN
      DR%X = 1.E6_EB
      CYCLE DROPLET_LOOP
   ENDIF

   ! Interpolate the nearest velocity components of the gas

   IIX  = FLOOR(XI+.5_EB)
   JJY  = FLOOR(YJ+.5_EB)
   KKZ  = FLOOR(ZK+.5_EB)
   UBAR = AFILL2(U,II-1,JJY,KKZ,XI-II+1,YJ-JJY+.5_EB,ZK-KKZ+.5_EB)
   VBAR = AFILL2(V,IIX,JJ-1,KKZ,XI-IIX+.5_EB,YJ-JJ+1,ZK-KKZ+.5_EB)
   WBAR = AFILL2(W,IIX,JJY,KK-1,XI-IIX+.5_EB,YJ-JJY+.5_EB,ZK-KK+1)
    
   ! If the particle is just a massless tracer, just move it and go on to the next particle

   IF (PC%MASSLESS) THEN
      DR%U = UBAR
      DR%V = VBAR
      DR%W = WBAR
      DR%X = DR%X + DR%U*DT
      DR%Y = DR%Y + DR%V*DT
      DR%Z = DR%Z + DR%W*DT
      CYCLE DROPLET_LOOP
   ENDIF

   ! Calculate the particle velocity components and the amount of momentum to transfer to the gas

   RVC = RDX(II)*RDY(JJ)*RDZ(KK)
   RD  = DR%R
   RDS = RD*RD
   RDC = RD*RDS

   UREL  = DR%U - UBAR
   VREL  = DR%V - VBAR
   WREL  = DR%W - WBAR
   QREL  = MAX(1.E-6_EB,SQRT(UREL*UREL + VREL*VREL + WREL*WREL))
   TMP_G = MAX(TMPMIN,TMP(II,JJ,KK))
   RHO_G = RHO(II,JJ,KK)
   MU_AIR = SPECIES(0)%MU(MIN(500,NINT(0.1_EB*TMP_G)))
   DR%RE  = RHO_G*QREL*2._EB*RD/MU_AIR
   DRAG_CALC: IF (DR%IOR==0 .AND. DR%RE>0) THEN
      C_DRAG  = DRAG(DR%RE)
!rm ->
      IF (.NOT. PC%TREE) THEN
!rm <-
      SFAC    = DR%PWT*RDS*PIO2*QREL*C_DRAG
      DR%A_X  = SFAC*UREL*RVC
      DR%A_Y  = SFAC*VREL*RVC
      DR%A_Z  = SFAC*WREL*RVC
!rm ->
      ENDIF
      IF (PC%TREE) THEN
      SFAC    = PC%VEG_DRAG_COEFFICIENT*PC%VEG_SV*DR%VEG_PACKING_RATIO*QREL*C_DRAG
      DR%A_X  = SFAC*UREL
      DR%A_Y  = SFAC*VREL
      DR%A_Z  = SFAC*WREL
      ENDIF
!rm <-

      IF (.NOT.PC%STATIC) THEN
         CONST   = 8._EB*PC%DENSITY*RD/(3._EB*RHO_G*C_DRAG*QREL)
         BFAC    = EXP(-DT/CONST)
         AUREL   = CONST*GRVT1
         AVREL   = CONST*GRVT2
         AWREL   = CONST*GRVT3
         DR%U    = UBAR + (UREL+AUREL)*BFAC - AUREL
         DR%V    = VBAR + (VREL+AVREL)*BFAC - AVREL
         DR%W    = WBAR + (WREL+AWREL)*BFAC - AWREL
      ENDIF
   ELSE
      DR%A_X  = 0._EB
      DR%A_Y  = 0._EB
      DR%A_Z  = 0._EB
   ENDIF DRAG_CALC

   ! If the particle does not move, but does drag, go on to the next particle

   IF (PC%STATIC) CYCLE DROPLET_LOOP
 
   ! Decide how many time steps to use in tracking particle
    
   DTMIN = DT
   UVW = MAX( ABS(DR%U)*RDX(II),ABS(DR%V)*RDY(JJ),ABS(DR%W)*RDZ(KK) )
   IF (UVW/=0._EB) DTMIN = MIN(DTMIN,1._EB/UVW)
    
   NITER = 1
   DTSP  = DT
   NLOOP: DO N=0,3
      IF (DTMIN<DT*0.5_EB**N) THEN
         NITER = 2**(N+1)
         DTSP =  DT*0.5_EB**(N+1)
      ENDIF
   ENDDO NLOOP
    
   ! Iterate over a single time step
    
   SUB_TIME_STEP_ITERATIONS: DO N=1,NITER
    
      ! Get current particle coordinates

      IF (N>1) THEN
         XI = CELLSI(FLOOR((DR%X-XS)*RDXINT))
         YJ = CELLSJ(FLOOR((DR%Y-YS)*RDYINT))
         ZK = CELLSK(FLOOR((DR%Z-ZS)*RDZINT))
         II = FLOOR(XI+1._EB)
         JJ = FLOOR(YJ+1._EB)
         KK = FLOOR(ZK+1._EB)
         IC = CELL_INDEX(II,JJ,KK)
      ENDIF
    
      ! Update droplet position

      X_OLD = DR%X
      Y_OLD = DR%Y
      Z_OLD = DR%Z
      DR%X  = DR%X + DR%U*DTSP
      DR%Y  = DR%Y + DR%V*DTSP
      DR%Z  = DR%Z + DR%W*DTSP
    
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
    
      ! Where is the droplet now?
    
      XI  = CELLSI(FLOOR((DR%X-XS)*RDXINT))
      YJ  = CELLSJ(FLOOR((DR%Y-YS)*RDYINT))
      ZK  = CELLSK(FLOOR((DR%Z-ZS)*RDZINT))
      IIN = FLOOR(XI+1._EB)
      JJN = FLOOR(YJ+1._EB)
      KKN = FLOOR(ZK+1._EB)
      IF (IIN<0 .OR. IIN>IBP1) CYCLE DROPLET_LOOP
      IF (JJN<0 .OR. JJN>JBP1) CYCLE DROPLET_LOOP
      IF (KKN<0 .OR. KKN>KBP1) CYCLE DROPLET_LOOP
      ICN = CELL_INDEX(IIN,JJN,KKN)

      IF (IC==0 .OR. ICN==0) CYCLE SUB_TIME_STEP_ITERATIONS

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
            STEP_FRACTION(DR%IOR) = MAX(0._EB,(Z(KKN)-Z_OLD+0.05_EB*DZ(KKN))/(DR%Z-Z_OLD))
         ENDIF
         IF (IIN>II .AND. BOUNDARY_TYPE(IWP1)==SOLID_BOUNDARY) THEN
            DR%IOR=-1
            HIT_SOLID = .TRUE.
            STEP_FRACTION(DR%IOR) = MAX(0._EB,(X(II)-X_OLD-0.05_EB*DX(II))/(DR%X-X_OLD))
         ENDIF
         IF (IIN<II .AND. BOUNDARY_TYPE(IWM1)==SOLID_BOUNDARY) THEN
            DR%IOR= 1
            HIT_SOLID = .TRUE.
            STEP_FRACTION(DR%IOR) = MAX(0._EB,(X(IIN)-X_OLD+0.05_EB*DX(IIN))/(DR%X-X_OLD))
         ENDIF
         IF (JJN>JJ .AND. BOUNDARY_TYPE(IWP2)==SOLID_BOUNDARY) THEN
            DR%IOR=-2
            HIT_SOLID = .TRUE.
            STEP_FRACTION(DR%IOR) = MAX(0._EB,(Y(JJ)-Y_OLD-0.05_EB*DY(JJ))/(DR%Y-Y_OLD))
         ENDIF
         IF (JJN<JJ .AND. BOUNDARY_TYPE(IWM2)==SOLID_BOUNDARY) THEN
            DR%IOR= 2
            HIT_SOLID = .TRUE.
            STEP_FRACTION(DR%IOR) = MAX(0._EB,(Y(JJN)-Y_OLD+0.05_EB*DY(JJN))/(DR%Y-Y_OLD))
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
         
         ! If no solids boundaries of original cell have been crossed, check boundaries of new grid cell
    
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
               STEP_FRACTION(DR%IOR) = MAX(0._EB,(Z(KKN)-Z_OLD+0.05_EB*DZ(KKN))/(DR%Z-Z_OLD))
            ENDIF
            IF (IIN>II .AND. BOUNDARY_TYPE(IWM1)==SOLID_BOUNDARY) THEN
               DR%IOR=-1
               HIT_SOLID = .TRUE.
               STEP_FRACTION(DR%IOR) = MAX(0._EB,(X(II)-X_OLD-0.05_EB*DX(II))/(DR%X-X_OLD))
            ENDIF
            IF (IIN<II .AND. BOUNDARY_TYPE(IWP1)==SOLID_BOUNDARY) THEN
               DR%IOR= 1
               HIT_SOLID = .TRUE.
               STEP_FRACTION(DR%IOR) = MAX(0._EB,(X(IIN)-X_OLD+0.05_EB*DX(IIN))/(DR%X-X_OLD))
            ENDIF
            IF (JJN>JJ .AND. BOUNDARY_TYPE(IWM2)==SOLID_BOUNDARY) THEN
               DR%IOR=-2
               HIT_SOLID = .TRUE.
               STEP_FRACTION(DR%IOR) = MAX(0._EB,(Y(JJ)-Y_OLD-0.05_EB*DY(JJ))/(DR%Y-Y_OLD))
            ENDIF
            IF (JJN<JJ .AND. BOUNDARY_TYPE(IWP2)==SOLID_BOUNDARY) THEN
               DR%IOR= 2
               HIT_SOLID = .TRUE.
               STEP_FRACTION(DR%IOR) = MAX(0._EB,(Y(JJN)-Y_OLD+0.05_EB*DY(JJN))/(DR%Y-Y_OLD))
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
 
            IF (DR%WALL_INDEX==0) CYCLE SUB_TIME_STEP_ITERATIONS

            ! Add droplet mass to accumulated liquid array

            IF (ACCUMULATE_WATER .AND. HIT_SOLID .AND. .NOT.DR%SPLAT) THEN
               AWMPUA(DR%WALL_INDEX,PC%EVAP_INDEX) = AWMPUA(DR%WALL_INDEX,PC%EVAP_INDEX) + DR%PWT*PC%FTPR*DR%R**3*RAW(DR%WALL_INDEX)
               DR%SPLAT = .TRUE.
            ENDIF

            ! Adjust the size of the droplet and weighting factor 

            R_NEW  = MIN(0.5_EB*SURFACE_DROPLET_DIAMETER,(DR%PWT*DR%R**3)**ONTH)
            DR%PWT = DR%PWT*DR%R**3/R_NEW**3
            DR%R   = R_NEW

            ! Move particle to where it almost hits solid

            DR%X = X_OLD + MINVAL(STEP_FRACTION)*DTSP*DR%U
            DR%Y = Y_OLD + MINVAL(STEP_FRACTION)*DTSP*DR%V
            DR%Z = Z_OLD + MINVAL(STEP_FRACTION)*DTSP*DR%W

            XI  = CELLSI(FLOOR((DR%X-XS)*RDXINT))
            YJ  = CELLSJ(FLOOR((DR%Y-YS)*RDYINT))
            ZK  = CELLSK(FLOOR((DR%Z-ZS)*RDZINT))
            IIN = FLOOR(XI+1._EB)
            JJN = FLOOR(YJ+1._EB)
            KKN = FLOOR(ZK+1._EB)
            ICN = CELL_INDEX(IIN,JJN,KKN)

            IF (IOR_OLD==DR%IOR) CYCLE SUB_TIME_STEP_ITERATIONS

            ! Choose a direction for the droplets to move

            DIRECTION: SELECT CASE(ABS(DR%IOR))
               CASE (1:2) DIRECTION  
                  DR%U = 0._EB
                  DR%V = 0._EB
                  DR%W = -PC%VERTICAL_VELOCITY 
               CASE(3) DIRECTION
                  CALL RANDOM_NUMBER(RN)
                  THETA_RN = TWOPI*RN
                  DR%U = PC%HORIZONTAL_VELOCITY*COS(THETA_RN)
                  DR%V = PC%HORIZONTAL_VELOCITY*SIN(THETA_RN)
                  DR%W = 0._EB
            END SELECT DIRECTION

         ENDIF IF_HIT_SOLID
      ENDIF AIR_TO_SOLID 

      ! Check if droplets that were attached to a solid are still attached after the time update

      IW = 0
      SELECT CASE(DR%IOR)
         CASE( 1)
            IW = WALL_INDEX(ICN,-1)
            IF (BOUNDARY_TYPE(IW)/=SOLID_BOUNDARY) THEN
               DR%X = DR%X - 0.2_EB*DX(II)
               DR%W = -DR%W
            ENDIF
         CASE(-1)
            IW = WALL_INDEX(ICN, 1)
            IF (BOUNDARY_TYPE(IW)/=SOLID_BOUNDARY) THEN
               DR%X = DR%X + 0.2_EB*DX(II)
               DR%W = -DR%W
            ENDIF
         CASE( 2)
            IW = WALL_INDEX(ICN,-2)
            IF (BOUNDARY_TYPE(IW)/=SOLID_BOUNDARY) THEN
               DR%Y = DR%Y - 0.2_EB*DY(JJ)
               DR%W = -DR%W
            ENDIF
         CASE(-2)
            IW = WALL_INDEX(ICN, 2)
            IF (BOUNDARY_TYPE(IW)/=SOLID_BOUNDARY) THEN
               DR%Y = DR%Y + 0.2_EB*DY(JJ)
               DR%W = -DR%W
            ENDIF
         CASE( 3)
            IW = WALL_INDEX(ICN,-3)
            IF (BOUNDARY_TYPE(IW)/=SOLID_BOUNDARY) THEN  ! Particle has reached the edge of a horizontal surface
               DR%U = -DR%U
               DR%V = -DR%V
               DR%Z =  DR%Z - 0.2_EB*DZ(KK)
            ENDIF
         CASE(-3)
            IW = WALL_INDEX(ICN, 3)
      END SELECT

      IF (DR%IOR/=0 .AND. BOUNDARY_TYPE(IW)/=SOLID_BOUNDARY) THEN
         DR%IOR = 0
         DR%WALL_INDEX = 0
      ELSE
         DR%WALL_INDEX = WALL_INDEX(ICN,-DR%IOR)
      ENDIF

   ENDDO SUB_TIME_STEP_ITERATIONS

ENDDO DROPLET_LOOP

! Remove out-of-bounds particles
 
CALL REMOVE_DROPLETS(T,NM)

END SUBROUTINE MOVE_PARTICLES



SUBROUTINE PARTICLE_MASS_ENERGY_TRANSFER(T,NM)
    
! Mass and energy transfer between gas and droplets

USE PHYSICAL_FUNCTIONS, ONLY : GET_MASS_FRACTION

REAL(EB), POINTER, DIMENSION(:,:,:) :: DROP_DEN,DROP_RAD,DROP_TMP,MVAP_TOT
REAL(EB), POINTER, DIMENSION(:) :: FILM_THICKNESS
REAL(EB) :: R_DROP,NUSSELT,K_AIR,H_V, &
            RVC,WGT,OMWGT,Q_CON_GAS,Q_CON_WALL,Q_RAD,H_HEAT,H_MASS,SH_FAC_GAS,SH_FAC_WALL,NU_FAC_GAS,NU_FAC_WALL, &
            T,PR_AIR,M_VAP,M_VAP_MAX,XI,YJ,ZK,RDT,MU_AIR,H_SOLID,Q_DOT_RAD,DEN_ADD, &
            Y_DROP,Y_GAS,LENGTH,U2,V2,W2,VEL,DENOM,DY_DTMP_DROP,TMP_DROP_NEW,TMP_WALL,H_WALL, &
            SC_AIR,D_AIR,DHOR,SHERWOOD,X_DROP,M_DROP,RHO_G,MW_RATIO,MW_DROP,FTPR,&
            C_DROP,M_GAS,A_DROP,TMP_G,TMP_DROP,TMP_MELT,TMP_BOIL,MINIMUM_FILM_THICKNESS,RE_L,OMRAF,Q_FRAC,Q_TOT,DT_SUBSTEP
INTEGER :: I,II,JJ,KK,IW,IGAS,N_PC,EVAP_INDEX,ITER,N_SUBSTEPS
INTEGER, INTENT(IN) :: NM
REAL(EB), PARAMETER :: RUN_AVG_FAC=0.5_EB

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

IF (N_EVAP_INDICIES>0) THEN
   D_VAP  = 0._EB
   WCPUA  = RUN_AVG_FAC*WCPUA
   WMPUA  = RUN_AVG_FAC*WMPUA
   DROP_DEN => WORK4
   DROP_RAD => WORK5
   DROP_TMP => WORK6
   MVAP_TOT => WORK7
   MVAP_TOT = 0._EB
ENDIF

! Loop over all types of evaporative species

EVAP_INDEX_LOOP: DO EVAP_INDEX = 1,N_EVAP_INDICIES

   FILM_THICKNESS => WALL_WORK2
   FILM_THICKNESS =  0._EB

   ! Loop over all particle/droplet classes that have the given evaporative index

   PART_CLASS_LOOP: DO N_PC = 1,N_PART

      PC => PARTICLE_CLASS(N_PC)
      IF (PC%EVAP_INDEX/=EVAP_INDEX) CYCLE PART_CLASS_LOOP
      IF (PC%TREE) CYCLE PART_CLASS_LOOP
      DROP_DEN = 0._EB
      DROP_TMP = 0._EB
      DROP_RAD = 0._EB
      C_DROP   = PC%C_P
      FTPR     = PC%FTPR
      TMP_MELT = PC%TMP_MELT
      TMP_BOIL = PC%TMP_V
      IGAS     = PC%SPEC_INDEX
      MW_DROP  = SPECIES(IGAS)%MW
      MW_RATIO = SPECIES(0)%MW/MW_DROP
      H_V      = PC%H_V
      DHOR     = H_V*MW_DROP/R0

      ! Loop through all droplets in the class and determine the depth of the liquid film on each surface cell

      MASS_SUMMING_LOOP: DO I=1,NLP
         DR => DROPLET(I)
         IF (DR%IOR==0)        CYCLE MASS_SUMMING_LOOP
         IF (DR%WALL_INDEX==0) CYCLE MASS_SUMMING_LOOP
         IF (DR%CLASS /= N_PC) CYCLE MASS_SUMMING_LOOP
         IF (DR%R<=0._EB)      CYCLE MASS_SUMMING_LOOP
         IW = DR%WALL_INDEX
         FILM_THICKNESS(IW) = FILM_THICKNESS(IW) + DR%PWT*DR%R**3/AW(IW)
      ENDDO MASS_SUMMING_LOOP
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
         RVC = RDX(II)*RDY(JJ)*RDZ(KK)

         ! Determine how many sub-time step iterations are needed and then iterate over the time step.
         ! This is not fully functional. Keep as a placeholder for now.

    !!!  N_SUBSTEPS = 0  
    !!!  GET_N: DO ITER=1,100
    !!!     N_SUBSTEPS = N_SUBSTEPS + 1
    !!!     IF (DT/REAL(N_SUBSTEPS,EB)<=10000000.*DR%R) EXIT GET_N
    !!!  ENDDO GET_N

         N_SUBSTEPS = 1
         DT_SUBSTEP = DT/REAL(N_SUBSTEPS,EB) 

         TIME_ITERATION_LOOP: DO ITER=1,N_SUBSTEPS
         
            ! Initialize droplet thermophysical data

            R_DROP   = DR%R
            M_DROP   = FTPR*R_DROP**3
            TMP_DROP = DR%TMP
            WGT      = DR%PWT
               
            ! Gas conditions

            TMP_G  = TMP(II,JJ,KK)
            RHO_G  = RHO(II,JJ,KK)
            MU_AIR = SPECIES(0)%MU(MIN(500,NINT(0.1_EB*TMP_G)))
            M_GAS  = RHO_G/RVC        
            M_VAP_MAX = (0.33_EB * M_GAS - MVAP_TOT(II,JJ,KK)) / WGT!limit to avoid diveregence errors
            K_AIR  = CPOPR*MU_AIR
            IF (IGAS>0) THEN
               Y_GAS = YY(II,JJ,KK,IGAS)
            ELSE
               Y_GAS = 0._EB
            ENDIF

            ! Set variables for heat transfer on solid

            SOLID_OR_GAS_PHASE: IF (DR%IOR/=0 .AND. DR%WALL_INDEX>0) THEN

               IW   = DR%WALL_INDEX
               A_DROP = M_DROP/(FILM_THICKNESS(IW)*PC%DENSITY)
               TMP_WALL = TMP_F(IW) 
               SELECT CASE(ABS(DR%IOR))
                  CASE(1)
                     V2 = 0.25_EB*(V(II,JJ,KK)+V(II,JJ-1,KK))**2
                     W2 = 0.25_EB*(W(II,JJ,KK)+W(II,JJ,KK-1))**2
                     VEL = SQRT(V2+W2)
                  CASE(2)
                     U2 = 0.25_EB*(U(II,JJ,KK)+U(II-1,JJ,KK))**2
                     W2 = 0.25_EB*(W(II,JJ,KK)+W(II,JJ,KK-1))**2
                     VEL = SQRT(U2+W2)
                  CASE(3)
                     U2 = 0.25_EB*(U(II,JJ,KK)+U(II-1,JJ,KK))**2
                     V2 = 0.25_EB*(V(II,JJ,KK)+V(II,JJ-1,KK))**2
                     VEL = SQRT(U2+V2)
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
                  Q_DOT_RAD = QR_W(II,JJ,KK)/SUM(AVG_DROP_DEN(II,JJ,KK,:))*M_DROP
               ELSE
                  Q_DOT_RAD = 0._EB
               ENDIF

            ENDIF SOLID_OR_GAS_PHASE
         
            ! Compute equilibrium droplet vapor mass fraction, Y_DROP, and its derivative w.r.t. droplet temperature
    
            IF (PC%EVAPORATE) THEN
               X_DROP  = MIN(1._EB,EXP(DHOR*(1._EB/TMP_BOIL-1._EB/TMP_DROP)))
               Y_DROP  = X_DROP/(MW_RATIO + (1._EB-MW_RATIO)*X_DROP)
               IF (TMP_DROP < TMP_BOIL) THEN
                  DY_DTMP_DROP = (MW_RATIO/(X_DROP*(1._EB-MW_RATIO)+MW_RATIO)**2)*DHOR*X_DROP/TMP_DROP**2
               ELSE
                  DY_DTMP_DROP = 0._EB
               ENDIF
               IF (Y_DROP<=Y_GAS) H_MASS = 0._EB
            ELSE
               DY_DTMP_DROP = 0._EB
               Y_DROP       = 0._EB
            ENDIF

            ! Update the droplet temperature semi_implicitly
    
            DENOM = 1._EB + (H_HEAT + H_WALL + H_MASS*RHO_G*H_V*DY_DTMP_DROP)*DT_SUBSTEP*A_DROP/(2._EB*M_DROP*C_DROP) 

            TMP_DROP_NEW = ( TMP_DROP + DT_SUBSTEP*( Q_DOT_RAD + &
                             A_DROP*(H_HEAT*(TMP_G   -0.5_EB*TMP_DROP) + H_WALL*(TMP_WALL-0.5_EB*TMP_DROP) -  &
                             H_MASS*RHO_G*H_V*(Y_DROP-0.5_EB*DY_DTMP_DROP*TMP_DROP-Y_GAS))/(M_DROP*C_DROP)) ) / DENOM

            ! Compute the total amount of heat extracted from the gas, wall and radiative fields

            Q_RAD      = DT_SUBSTEP*Q_DOT_RAD
            Q_CON_GAS  = DT_SUBSTEP*A_DROP*H_HEAT*(TMP_G   -0.5_EB*(TMP_DROP+TMP_DROP_NEW))
            Q_CON_WALL = DT_SUBSTEP*A_DROP*H_WALL*(TMP_WALL-0.5_EB*(TMP_DROP+TMP_DROP_NEW))

            ! Compute the total amount of liquid evaporated
  
            M_VAP = DT_SUBSTEP*A_DROP*H_MASS*RHO_G*(Y_DROP+0.5_EB*DY_DTMP_DROP*(TMP_DROP_NEW-TMP_DROP)-Y_GAS) 
            M_VAP = MAX(0._EB,MIN(M_VAP,M_DROP,M_VAP_MAX))
            ! Evaporate completely small droplets

            IF (PC%EVAPORATE .AND. R_DROP<0.5_EB*PC%MINIMUM_DIAMETER) THEN  
               M_VAP  = M_DROP
               Q_TOT  = Q_RAD+Q_CON_GAS+Q_CON_WALL
               IF (Q_TOT>0._EB) THEN
                  Q_FRAC = M_VAP*H_V/Q_TOT 
                  Q_CON_GAS  = Q_CON_GAS*Q_FRAC
                  Q_CON_WALL = Q_CON_WALL*Q_FRAC
                  Q_RAD      = Q_RAD*Q_FRAC
               ENDIF
            ENDIF

            ! If the droplet temperature drops below its freezing point, just reset it

            IF (PC%EVAPORATE .AND. TMP_DROP_NEW<TMP_MELT) TMP_DROP_NEW = TMP_MELT

            ! If the droplet temperature reaches boiling, use only enough energy from gas to vaporize liquid

            IF (PC%EVAPORATE .AND. TMP_DROP_NEW>=TMP_BOIL) THEN  
               M_VAP  = MIN(M_VAP_MAX,M_DROP,M_VAP + (TMP_DROP_NEW - TMP_BOIL)*C_DROP*M_DROP/H_V)
               TMP_DROP_NEW = TMP_BOIL
               Q_TOT  = Q_RAD+Q_CON_GAS+Q_CON_WALL
               IF (Q_TOT>0._EB) THEN
                  Q_FRAC = M_VAP*H_V/Q_TOT
                  Q_CON_GAS  = Q_CON_GAS*Q_FRAC
                  Q_CON_WALL = Q_CON_WALL*Q_FRAC
                  Q_RAD      = Q_RAD*Q_FRAC
               ENDIF
            ENDIF

            ! Update droplet quantities
         
            M_DROP = M_DROP - M_VAP
            DR%R   = (M_DROP/FTPR)**ONTH
            DR%TMP = TMP_DROP_NEW

            ! Compute surface cooling and density

            IF (DR%IOR/=0 .AND. DR%WALL_INDEX>0) THEN
               WCPUA(IW,EVAP_INDEX) = WCPUA(IW,EVAP_INDEX) + OMRAF*WGT*(Q_RAD+Q_CON_WALL)*RAW(IW)/DT_SUBSTEP
               WMPUA(IW,EVAP_INDEX) = WMPUA(IW,EVAP_INDEX) + OMRAF*WGT*M_DROP*RAW(IW)/REAL(N_SUBSTEPS,EB) 
            ENDIF
         
            ! Decrease temperature of the gas cell

            TMP(II,JJ,KK) = (M_GAS*CP_GAMMA*TMP_G - WGT*Q_CON_GAS)/(M_GAS*CP_GAMMA)  
            TMP(II,JJ,KK) = MIN(TMPMAX,MAX(TMPMIN,TMP(II,JJ,KK)))

            ! Compute contribution to the divergence

            D_VAP(II,JJ,KK) = D_VAP(II,JJ,KK) + WGT*RVC/(DT_SUBSTEP*RHO_G)*( M_VAP*MW_RATIO-Q_CON_GAS/(CP_GAMMA*TMP(II,JJ,KK)) )

            ! Add fuel evaporation rate to running counter before adjusting its value

            IF (IGAS>0 .AND. IGAS==I_FUEL) FUEL_DROPLET_MLR(NM) = FUEL_DROPLET_MLR(NM) + WGT*M_VAP/DT_SUBSTEP

            ! Adjust mass of evaporated liquid to account for different Heat of Combustion between fuel droplet and gas

            M_VAP = PC%ADJUST_EVAPORATION*M_VAP

            ! Add vapor or fuel gas to the grid cell

            IF (IGAS>0) YY(II,JJ,KK,IGAS)= MIN(1._EB,(WGT*M_VAP+M_GAS*Y_GAS)/(WGT*M_VAP+M_GAS))

            ! Add new mass from vaporized droplet to the grid cell

            RHO(II,JJ,KK) = RHO(II,JJ,KK) + WGT*M_VAP*RVC
            
            !Track total mass evaporate in cell
            MVAP_TOT(II,JJ,KK) = MVAP_TOT(II,JJ,KK) + WGT*M_VAP

            ! Get out of the loop if the droplet has evaporated completely
            IF (DR%R<=0._EB) CYCLE DROPLET_LOOP

         ENDDO TIME_ITERATION_LOOP

         ! Assign liquid mass to the cell for airborne drops

         IF (DR%IOR==0) THEN
            DEN_ADD = WGT*M_DROP*RVC
            DROP_DEN(II,JJ,KK) = DROP_DEN(II,JJ,KK) + DEN_ADD
            DROP_TMP(II,JJ,KK) = DROP_TMP(II,JJ,KK) + DEN_ADD*DR%TMP
            DROP_RAD(II,JJ,KK) = DROP_RAD(II,JJ,KK) + DEN_ADD*DR%R
         ENDIF
         
      ENDDO DROPLET_LOOP
      ! Calculate cumulative quantities for droplet "clouds"

      WGT   = .5_EB
      OMWGT = 1._EB-WGT
    
      DROP_RAD = DROP_RAD/(DROP_DEN+1.E-10_EB)
      DROP_TMP = DROP_TMP/(DROP_DEN+1.E-10_EB)
    
      AVG_DROP_RAD(:,:,:,EVAP_INDEX ) = (AVG_DROP_DEN(:,:,:,EVAP_INDEX )*AVG_DROP_RAD(:,:,:,EVAP_INDEX )+DROP_DEN*DROP_RAD) &
                                          /(AVG_DROP_DEN(:,:,:,EVAP_INDEX ) + DROP_DEN + 1.0E-10_EB)
      AVG_DROP_TMP(:,:,:,EVAP_INDEX ) = (AVG_DROP_DEN(:,:,:,EVAP_INDEX )*AVG_DROP_TMP(:,:,:,EVAP_INDEX )+DROP_DEN*DROP_TMP) &
                                          /(AVG_DROP_DEN(:,:,:,EVAP_INDEX ) + DROP_DEN + 1.0E-10_EB)
      AVG_DROP_TMP(:,:,:,EVAP_INDEX ) = MAX(TMPM,AVG_DROP_TMP(:,:,:,EVAP_INDEX ))
      AVG_DROP_DEN(:,:,:,EVAP_INDEX ) = WGT*AVG_DROP_DEN(:,:,:,EVAP_INDEX ) + OMWGT*DROP_DEN
      WHERE (AVG_DROP_DEN(:,:,:,EVAP_INDEX )<0.0001_EB .AND. DROP_DEN==0._EB) AVG_DROP_DEN(:,:,:,EVAP_INDEX ) = 0.0_EB

   ENDDO PART_CLASS_LOOP

ENDDO EVAP_INDEX_LOOP

! Remove droplets that have completely evaporated
 
CALL REMOVE_DROPLETS(T,NM)
 
END SUBROUTINE PARTICLE_MASS_ENERGY_TRANSFER



SUBROUTINE PARTICLE_MOMENTUM_TRANSFER

! Add droplet momentum as a force term in momentum equation

USE COMP_FUNCTIONS, ONLY : SECOND
REAL(EB), POINTER, DIMENSION(:,:,:) :: FVXS,FVYS,FVZS
REAL(EB) :: FLUXMIN,FLUXMAX,XI,YJ,ZK
INTEGER :: II,JJ,KK,IIX,JJY,KKZ,I,J,K,IC,IW


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
   IF (BOUNDARY_TYPE(IW)==NULL_BOUNDARY) FVXS(IIX,JJ,KK) = FVXS(IIX,JJ,KK) - DR%A_X
   IC = CELL_INDEX(II,JJY,KK)
   IW = WALL_INDEX(IC,2) 
   IF (BOUNDARY_TYPE(IW)==NULL_BOUNDARY) FVYS(II,JJY,KK) = FVYS(II,JJY,KK) - DR%A_Y
   IC = CELL_INDEX(II,JJ,KKZ)
   IW = WALL_INDEX(IC,3) 
   IF (BOUNDARY_TYPE(IW)==NULL_BOUNDARY) FVZS(II,JJ,KKZ) = FVZS(II,JJ,KKZ) - DR%A_Z
ENDDO SUM_MOMENTUM_LOOP
 
FLUXMAX =  100._EB
FLUXMIN = -FLUXMAX
 
DO K=0,KBAR
   DO J=0,JBAR
      DO I=0,IBAR
         IF (FVXS(I,J,K)>FLUXMAX) FVXS(I,J,K) = FLUXMAX
         IF (FVXS(I,J,K)<FLUXMIN) FVXS(I,J,K) = FLUXMIN
         FVX(I,J,K) = FVX(I,J,K) + FVXS(I,J,K)
         IF (FVYS(I,J,K)>FLUXMAX) FVYS(I,J,K) = FLUXMAX
         IF (FVYS(I,J,K)<FLUXMIN) FVYS(I,J,K) = FLUXMIN
         FVY(I,J,K) = FVY(I,J,K) + FVYS(I,J,K)
         IF (FVZS(I,J,K)>FLUXMAX) FVZS(I,J,K) = FLUXMAX
         IF (FVZS(I,J,K)<FLUXMIN) FVZS(I,J,K) = FLUXMIN
         FVZ(I,J,K) = FVZ(I,J,K) + FVZS(I,J,K)
      ENDDO
   ENDDO
ENDDO
 
END SUBROUTINE PARTICLE_MOMENTUM_TRANSFER
 
 
 
SUBROUTINE REMOVE_DROPLETS(T,NM)
USE MEMORY_FUNCTIONS, ONLY : RE_ALLOCATE_DROPLETS 
! Remove droplets that do not lie in any mesh
 
INTEGER :: IKILL,I,NM,IPC
REAL(EB) :: T,MAX_DROP_LIFE

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
TYPE (MESH_TYPE), POINTER :: M
TYPE (OMESH_TYPE), POINTER :: M2

NOM = 0
SEARCH_LOOP: DO OM=1,NMESHES
   IF (NIC(NM,OM)==0) CYCLE SEARCH_LOOP
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
