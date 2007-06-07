!     Last change:  JEF  25 May 2007   11:07 am
MODULE PART
 
USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_POINTERS
USE TRAN
USE TYPES, ONLY: DROPLET_TYPE, PARTICLE_CLASS_TYPE, PARTICLE_CLASS
IMPLICIT NONE
PRIVATE
PUBLIC INSERT_DROPLETS_AND_PARTICLES,UPDATE_PARTICLES, &
       INITIALIZE_DROPLETS,INITIALIZE_TREES
TYPE (DROPLET_TYPE), POINTER :: DR
TYPE (PARTICLE_CLASS_TYPE), POINTER :: PC
CHARACTER(255), PARAMETER :: partid='$Id: part.f90,v 1.19 2007/04/26 12:50:34 mcgratta Exp $'
 
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

IF (N_PART==0) RETURN !Don't waste time if no particles
IF (EVACUATION_ONLY(NM)) RETURN !Don't waste time if an evac mesh
 
TNOW=SECOND()
CALL POINT_TO_MESH(NM)
 
PART_CLASS_LOOP: DO IPC=1,N_PART
 
   PC=>PARTICLE_CLASS(IPC)
   ! Check if droplets evaporate into a SPECIES
 
   SPECIES_LOOP: DO I=1,N_SPECIES
      IF (PC%SPECIES=='null') CYCLE SPECIES_LOOP
      IF (SPECIES(I)%NAME==PC%SPECIES) THEN
         PC%SPECIES_INDEX = I
         EXIT SPECIES_LOOP
      ENDIF
   ENDDO SPECIES_LOOP
 
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
         CALL RE_ALLOCATE_DROPLETS(1,NM,0)
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
 
      IF (PC%DIAMETER>0._EB) THEN
         STRATUM = MOD(NLP,NSTRATA) + 1
         IL = PC%IL_CDF(STRATUM)
         IU = PC%IU_CDF(STRATUM)
         CALL RANDOM_CHOICE(PC%CDF(IL:IU),PC%R_CDF(IL:IU),IU-IL,DR%R)
         DR%PWT = PC%W_CDF(STRATUM)
         MASS_SUM = MASS_SUM + DR%PWT*PC%FTPR*DR%R**3
      ENDIF
 
   ENDDO INITIALIZATION_LOOP
 
   ! Adjust particle weighting factor PWT so that desired MASS_PER_VOLUME is achieved
   IF (PC%DIAMETER>0._EB) DROPLET(NLP-PC%N_INITIAL+1:NLP)%PWT = &
      DROPLET(NLP-PC%N_INITIAL+1:NLP)%PWT*PC%MASS_PER_VOLUME*VOL1/MASS_SUM
ENDDO PART_CLASS_LOOP
 
TUSED(8,NM)=TUSED(8,NM)+SECOND()-TNOW
END SUBROUTINE INITIALIZE_DROPLETS
 
 

SUBROUTINE INITIALIZE_TREES(NM)
USE MEMORY_FUNCTIONS, ONLY: RE_ALLOCATE_DROPLETS 
REAL(EB) DELTAZ_BIN,CANOPY_LENGTH,CANOPY_VOLUME,COSINE,TANGENT, &
         XC_MAX_ZC_MAX,XC_MAX_ZC_MIN,XC_MIN,XC_MAX, &
         YC_MIN,YC_MAX,ZBIN_VOLUME,ZC_MIN,ZC_MAX,CANOPY_WIDTH
INTEGER NCT,NNLP,NLP_TREE,NZB,NZBINS,P_PER_ZBIN,IPC
REAL(EB) RN,MASS_SUM
INTEGER I,II,JJ,KK,IL,IU,STRATUM
INTEGER, INTENT(IN) :: NM
 
CALL POINT_TO_MESH(NM)
 
TREE_LOOP: DO NCT=1,N_TREES
 
   IF (TREE_MESH(NCT)/=NM) CYCLE TREE_LOOP
   IPC = TREE_PARTICLE_CLASS(NCT)
   PC=>PARTICLE_CLASS(IPC)
 
   ! Build a conical tree
 
   CANOPY_WIDTH  = CANOPY_W(NCT)
   CANOPY_LENGTH = TREE_H(NCT) - CANOPY_B_H(NCT)
   TANGENT = 0.5_EB*CANOPY_W(NCT)/CANOPY_LENGTH
   COSINE = CANOPY_LENGTH/SQRT(CANOPY_LENGTH**2 + (0.5_EB*CANOPY_WIDTH)**2)
   NZBINS = 1._EB/(1._EB-COSINE)
   DELTAZ_BIN = CANOPY_LENGTH/REAL(NZBINS,EB)
   CANOPY_VOLUME = PI*CANOPY_WIDTH**2*CANOPY_LENGTH/12.
 
   NLP_TREE = 0
   ZBIN_LOOP: DO NZB=1,NZBINS
      ZC_MIN = CANOPY_B_H(NCT) + (NZB-1)*DELTAZ_BIN
      ZC_MAX = ZC_MIN + DELTAZ_BIN
      XC_MAX_ZC_MIN = (CANOPY_LENGTH-(ZC_MIN-CANOPY_B_H(NCT)))*TANGENT
      XC_MAX_ZC_MAX = (CANOPY_LENGTH-(ZC_MAX-CANOPY_B_H(NCT)))*TANGENT
      ZBIN_VOLUME = PI*DELTAZ_BIN/3._EB*(XC_MAX_ZC_MIN**2 + XC_MAX_ZC_MAX**2 +XC_MAX_ZC_MIN*XC_MAX_ZC_MAX )
      P_PER_ZBIN  = PC%N_INITIAL*ZBIN_VOLUME/CANOPY_VOLUME
 
      NISP_LOOP1: DO NNLP=1,P_PER_ZBIN
         NLP  = NLP + 1
         PARTICLE_TAG = PARTICLE_TAG + NMESHES
         NLP_TREE = NLP_TREE + 1
         IF (NLP>NLPDIM) THEN
            CALL RE_ALLOCATE_DROPLETS(1,NM,0)
            DROPLET=>MESHES(NM)%DROPLET
            ENDIF
         DR=>DROPLET(NLP)
         DR%TAG = PARTICLE_TAG
         BLK_LOOP:  DO
            CALL RANDOM_NUMBER(RN)
            DR%Z = ZC_MIN + RN*(ZC_MAX-ZC_MIN)
            CALL RANDOM_NUMBER(RN)
            XC_MAX = (CANOPY_LENGTH-(DR%Z-CANOPY_B_H(NCT)))*TANGENT
            XC_MIN = -XC_MAX
            DR%X = XC_MIN + RN*(XC_MAX-XC_MIN)
            CALL RANDOM_NUMBER(RN)
            YC_MAX = SQRT(XC_MAX**2 - DR%X**2)
            YC_MIN = -YC_MAX
            DR%Y = YC_MIN + RN*(YC_MAX-YC_MIN)
            DR%X = X_TREE(NCT) + DR%X
            DR%Y = Y_TREE(NCT) + DR%Y
            DR%Z = Z_TREE(NCT) + DR%Z
            II = CELLSI(FLOOR((DR%X - XS)*RDXINT)) + 1._EB
            JJ = CELLSJ(FLOOR((DR%Y - YS)*RDYINT)) + 1._EB
            KK = CELLSK(FLOOR((DR%Z - ZS)*RDZINT)) + 1._EB
            IF (.NOT.SOLID(CELL_INDEX(II,JJ,KK))) EXIT BLK_LOOP
         ENDDO BLK_LOOP
      ENDDO NISP_LOOP1
   ENDDO ZBIN_LOOP
 
! Define physical characteristics of fuel elements
 
   MASS_SUM = 0._EB
   NISP_LOOP2: DO I=NLP-NLP_TREE+1,NLP
      DR=>DROPLET(I)
      DR%U = 0._EB
      DR%V = 0._EB
      DR%W = 0._EB
      DR%TMP = PC%TMP_INITIAL
      DR%T   = 0._EB
      DR%IOR = 0_EB
      DR%A_X = 0._EB
      DR%A_Y = 0._EB
      DR%A_Z = 0._EB
      DR%CLASS = IPC
      IF (MOD(I,PC%SAMPLING)==0) THEN
         DR%SHOW = .TRUE.
      ELSE
         DR%SHOW = .FALSE.
      ENDIF
      STRATUM = MOD(I,NSTRATA) + 1
      IL = PC%IL_CDF(STRATUM)
      IU = PC%IU_CDF(STRATUM)
      CALL RANDOM_CHOICE(PC%CDF(IL:IU),PC%R_CDF(IL:IU),IU-IL,DR%R)
      DR%PWT = PC%W_CDF(STRATUM)
      MASS_SUM = MASS_SUM + DR%PWT*PC%FTPR*DR%R**3
   ENDDO NISP_LOOP2
 
   ! Adjust particle weighting factor PWT so that desired MASS_PER_VOLUME is achieved
 
   DROPLET(NLP-NLP_TREE+1:NLP)%PWT = DROPLET(NLP-NLP_TREE+1:NLP)%PWT* PC%MASS_PER_VOLUME*CANOPY_VOLUME/MASS_SUM
 
ENDDO TREE_LOOP
 
END SUBROUTINE INITIALIZE_TREES
 
 

SUBROUTINE INSERT_DROPLETS_AND_PARTICLES(T,NM)
USE COMP_FUNCTIONS, ONLY : SECOND 
! Insert sprinkler droplets and lagrangian particles into the domain every DT_INSERT seconds
 
USE MEMORY_FUNCTIONS, ONLY : RE_ALLOCATE_DROPLETS
USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
USE DEVICE_VARIABLES
REAL(EB), INTENT(IN) :: T
INTEGER, INTENT(IN) :: NM
REAL(EB) :: PHI_RN,RN,FLOW_RATE,THETA_RN,SPHI,CPHI,MASS_SUM, &
            STHETA,CTHETA,PWT0,DROPLET_SPEED,XI,YJ,ZK,SHIFT1,SHIFT2,XTMP,YTMP,ZTMP,VLEN, &
            TRIGT1,TRIGT2,TNOW,SOLID_ANGLE_FLOWRATE,SOLID_ANGLE_TOTAL_FLOWRATE
REAL(EB), PARAMETER :: VENT_OFFSET=0.1
INTEGER :: I,KS,II,JJ,KK,IC,IL,IU,STRATUM,IPC,DROP_SUM,IIG,JJG,KKG,IW,IBC,IOR
LOGICAL :: INSERT_ANOTHER_BATCH
TYPE (PROPERTY_TYPE), POINTER :: PY
TYPE (TABLES_TYPE), POINTER :: TA
TYPE (DEVICE_TYPE), POINTER :: DV
TYPE (SURFACE_TYPE), POINTER :: SF
 
IF (N_PART==0) RETURN !Don't waste time if no particles
TNOW=SECOND()
CALL POINT_TO_MESH(NM)

OVERALL_INSERT_LOOP: DO  ! Add more than one batch of particles/droplets if FDS time step is large enough
INSERT_ANOTHER_BATCH = .FALSE.

SPRINKLER_INSERT_LOOP: DO KS=1,N_DEVC  ! Loop over all devices, but look for sprinklers or nozzles
   DV => DEVICE(KS)
   PY => PROPERTY(DV%PROP_INDEX)
   IF (PY%PART_ID == 'null')   CYCLE SPRINKLER_INSERT_LOOP
   IF (DV%MESH/=NM)            CYCLE SPRINKLER_INSERT_LOOP
   IF (.NOT. DV%CURRENT_STATE) CYCLE SPRINKLER_INSERT_LOOP
   IPC = PROPERTY(DV%PROP_INDEX)%PART_INDEX
   PC=>PARTICLE_CLASS(IPC)
   IF (T < PC%INSERT_CLOCK(NM)) CYCLE SPRINKLER_INSERT_LOOP
   FLOW_RATE = EVALUATE_RAMP(T-DV%T_CHANGE,PY%FLOW_TAU,PY%FLOW_RAMP_INDEX)*PY%FLOW_RATE*(PC%DENSITY/1000._EB)/60._EB  ! kg/s
   IF (FLOW_RATE <= 0._EB) CYCLE SPRINKLER_INSERT_LOOP
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
         CALL RE_ALLOCATE_DROPLETS(1,NM,0)
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
 
      STRATUM = MOD(NLP,NSTRATA) + 1
      IL = PC%IL_CDF(STRATUM)
      IU = PC%IU_CDF(STRATUM)
      CALL RANDOM_CHOICE(PC%CDF(IL:IU), PC%R_CDF(IL:IU), IU-IL,DR%R)
      DR%PWT = PC%W_CDF(STRATUM)
 
      ! Sum up mass of water being introduced
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
         PWT0 = FLOW_RATE*MAX(T-DV%T,0.01_EB)/(MASS_SUM*REAL(PC%N_INSERT,EB)/SOLID_ANGLE_TOTAL_FLOWRATE)
      ELSE
         PWT0 = FLOW_RATE*MAX(T-DV%T,0.01_EB)/MASS_SUM
      ENDIF
      DROPLET(NLP-DROP_SUM+1:NLP)%PWT = DROPLET(NLP-DROP_SUM+1:NLP)%PWT*PWT0
   ENDIF
   ! Indicate that droplets from this device have been inserted at this time T

   DV%T = T    
     
ENDDO SPRINKLER_INSERT_LOOP
 

! Loop through all boundary cells and insert particles if appropriate
 
WALL_INSERT_LOOP: DO IW=1,NWC
 
   IF (T < TW(IW))                       CYCLE WALL_INSERT_LOOP
   IF (BOUNDARY_TYPE(IW)==NULL_BOUNDARY) CYCLE WALL_INSERT_LOOP
   IBC =  IJKW(5,IW)
   SF  => SURFACE(IBC)
   IPC =  SF%PART_INDEX
   IF (IPC < 1)    CYCLE WALL_INSERT_LOOP
   PC  => PARTICLE_CLASS(IPC)
   IF (T < PC%INSERT_CLOCK(NM)) CYCLE WALL_INSERT_LOOP
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
         CALL RE_ALLOCATE_DROPLETS(1,NM,0)
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
      DR%R   = 0.0_EB
      DR%TMP = PC%TMP_INITIAL
      DR%T   = T
      DR%PWT = 1._EB
      DR%IOR = 0
 
      IF (PC%DIAMETER > 0._EB) THEN
         STRATUM = MOD(NLP,NSTRATA) + 1
         IL = PC%IL_CDF(STRATUM)
         IU = PC%IU_CDF(STRATUM)
         CALL RANDOM_CHOICE(PC%CDF(IL:IU),PC%R_CDF(IL:IU),IU-IL,DR%R)
         DR%PWT = PC%W_CDF(STRATUM)
         MASS_SUM = MASS_SUM + DR%PWT*PC%FTPR*DR%R**3
      ENDIF

   ENDDO PARTICLE_INSERT_LOOP

   IF (MASS_SUM > 0._EB) THEN
      IF (SF%PARTICLE_MASS_FLUX > 0._EB) THEN
         DROPLET(NLP-NPPCW(IW)+1:NLP)%PWT = DROPLET(NLP-NPPCW(IW)+1:NLP)%PWT*SF%PARTICLE_MASS_FLUX*AW(IW)/MASS_SUM/PC%DT_INSERT
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
IF (EVACUATION_ONLY(NM)) RETURN
TNOW=SECOND()
CALL POINT_TO_MESH(NM)

IF (CORRECTOR) CALL MOVE_PARTICLES(T,NM)
IF (CORRECTOR) CALL PARTICLE_MASS_ENERGY_TRANSFER(T,NM)
CALL PARTICLE_MOMENTUM_TRANSFER

TUSED(8,NM)=TUSED(8,NM)+SECOND()-TNOW
END SUBROUTINE UPDATE_PARTICLES



SUBROUTINE MOVE_PARTICLES(T,NM)

! Momentum and heat transfer from all particles and droplets
 
USE COMP_FUNCTIONS, ONLY : SECOND  
USE MATH_FUNCTIONS, ONLY : EVALUATE_RAMP, AFILL2
USE PHYSICAL_FUNCTIONS, ONLY : DRAG
 
REAL(EB) :: RHO_G,RVC,RDS,RDC,QREL,SFAC,UREL,VREL,WREL,TMP_G,RN,THETA_RN, &
            RD,RRD,T,C_DRAG,XI,YJ,ZK,RDT,MU_AIR, &
            DTSP,DTMIN,UBAR,VBAR,WBAR,BFAC,GRVT1,GRVT2,GRVT3,AUREL,AVREL,AWREL,CONST, &
            UVW,DUMMY=0._EB,X_OLD,Y_OLD,Z_OLD,STEP_FRACTION
LOGICAL :: HIT_SOLID
INTEGER :: ICN,I,IIN,JJN,KKN,II,JJ,KK,IIX,JJY,KKZ,IW,N,NITER, &
           IBC,IC,IWP1,IWM1,IWP2,IWM2,IWP3,IWM3,IOR_OLD
INTEGER, INTENT(IN) :: NM

GRVT1 = -EVALUATE_RAMP(T,DUMMY,I_RAMP_GX)*GVEC(1) 
GRVT2 = -EVALUATE_RAMP(T,DUMMY,I_RAMP_GY)*GVEC(2) 
GRVT3 = -EVALUATE_RAMP(T,DUMMY,I_RAMP_GZ)*GVEC(3) 
RDT   = 1._EB/DT

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
   IIX = FLOOR(XI+.5_EB)
   JJY = FLOOR(YJ+.5_EB)
   KKZ = FLOOR(ZK+.5_EB)
    
   ! Interpolate the nearest velocity components of the gas

   UBAR = AFILL2(U,II-1,JJY,KKZ,XI-II+1,YJ-JJY+.5_EB,ZK-KKZ+.5_EB)
   VBAR = AFILL2(V,IIX,JJ-1,KKZ,XI-IIX+.5_EB,YJ-JJ+1,ZK-KKZ+.5_EB)
   WBAR = AFILL2(W,IIX,JJY,KK-1,XI-IIX+.5_EB,YJ-JJY+.5_EB,ZK-KK+1)
    
   ! If the particle is just a massless tracer, just move it and get out

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
   IC  = CELL_INDEX(II,JJ,KK)
   RD  = DR%R
   RDS = RD*RD
   RDC = RD*RDS
   RRD = 1._EB/RD

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
      SFAC    = DR%PWT*RDS*PIO2*QREL*C_DRAG
      DR%A_X  = SFAC*UREL*RVC
      DR%A_Y  = SFAC*VREL*RVC
      DR%A_Z  = SFAC*WREL*RVC
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

   ! If the particle does not move, but does drag, get out now

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
         IF (BOUNDARY_TYPE(IW)==SOLID_BOUNDARY .AND. ACCUMULATE_WATER) THEN
            AWMPUA(IW,PC%EVAP_INDEX) = AWMPUA(IW,PC%EVAP_INDEX) + DR%PWT*PC%FTPR*RDC*RAW(IW)
         ENDIF
         CYCLE DROPLET_LOOP
      ENDIF
    
      IF (.NOT.POROUS_FLOOR .AND. DR%Z<ZS) THEN
         IC = CELL_INDEX(II,JJ,1)
         IW = WALL_INDEX(IC,-3)
         IF (BOUNDARY_TYPE(IW)==SOLID_BOUNDARY) THEN
            IF (ACCUMULATE_WATER) AWMPUA(IW,PC%EVAP_INDEX) = AWMPUA(IW,PC%EVAP_INDEX) + DR%PWT*PC%FTPR*RDC*RAW(IW)
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

      IF (.NOT.SOLID(ICN)) THEN
         IF (DR%X<XS .OR. DR%X>XF) CYCLE DROPLET_LOOP
         IF (DR%Y<YS .OR. DR%Y>YF) CYCLE DROPLET_LOOP
         IF (DR%Z<ZS .OR. DR%Z>ZF) CYCLE DROPLET_LOOP
      ENDIF
    
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

         IF (KKN>KK .AND. BOUNDARY_TYPE(IWP3)==SOLID_BOUNDARY) THEN
            DR%IOR=-3
            HIT_SOLID = .TRUE.
            STEP_FRACTION = (Z(KK)-Z_OLD-0.05_EB*DZ(KK))/(DR%Z-Z_OLD)
         ENDIF
         IF (KKN<KK .AND. BOUNDARY_TYPE(IWM3)==SOLID_BOUNDARY) THEN
            DR%IOR= 3
            HIT_SOLID = .TRUE.
            STEP_FRACTION = (Z(KKN)-Z_OLD+0.05_EB*DZ(KKN))/(DR%Z-Z_OLD)
         ENDIF
         IF (IIN>II .AND. BOUNDARY_TYPE(IWP1)==SOLID_BOUNDARY) THEN
            DR%IOR=-1
            HIT_SOLID = .TRUE.
            STEP_FRACTION = (X(II)-X_OLD-0.05_EB*DX(II))/(DR%X-X_OLD)
         ENDIF
         IF (IIN<II .AND. BOUNDARY_TYPE(IWM1)==SOLID_BOUNDARY) THEN
            DR%IOR= 1
            HIT_SOLID = .TRUE.
            STEP_FRACTION = (X(IIN)-X_OLD+0.05_EB*DX(IIN))/(DR%X-X_OLD)
         ENDIF
         IF (JJN>JJ .AND. BOUNDARY_TYPE(IWP2)==SOLID_BOUNDARY) THEN
            DR%IOR=-2
            HIT_SOLID = .TRUE.
            STEP_FRACTION = (Y(JJ)-Y_OLD-0.05_EB*DY(JJ))/(DR%Y-Y_OLD)
         ENDIF
         IF (JJN<JJ .AND. BOUNDARY_TYPE(IWM2)==SOLID_BOUNDARY) THEN
            DR%IOR= 2
            HIT_SOLID = .TRUE.
            STEP_FRACTION = (Y(JJN)-Y_OLD+0.05_EB*DY(JJN))/(DR%Y-Y_OLD)
         ENDIF

         DR%WALL_INDEX = WALL_INDEX(IC,-DR%IOR)
         
         ! If no solids boundaries of original cell have been crossed, check boundaries of new grid cell
    
         IF (DR%WALL_INDEX==0) THEN
            IWP1 = WALL_INDEX(ICN, 1)
            IWM1 = WALL_INDEX(ICN,-1)
            IWP2 = WALL_INDEX(ICN, 2)
            IWM2 = WALL_INDEX(ICN,-2)
            IWP3 = WALL_INDEX(ICN, 3)
            IWM3 = WALL_INDEX(ICN,-3)
            HIT_SOLID = .FALSE.
            IF (KKN>KK .AND. BOUNDARY_TYPE(IWM3)==SOLID_BOUNDARY) THEN
               DR%IOR=-3
               HIT_SOLID = .TRUE.
               STEP_FRACTION = (Z(KK)-Z_OLD-0.05_EB*DZ(KK))/(DR%Z-Z_OLD)
            ENDIF
            IF (KKN<KK .AND. BOUNDARY_TYPE(IWP3)==SOLID_BOUNDARY) THEN
               DR%IOR= 3
               HIT_SOLID = .TRUE.
               STEP_FRACTION = (Z(KKN)-Z_OLD+0.05_EB*DZ(KKN))/(DR%Z-Z_OLD)
            ENDIF
            IF (IIN>II .AND. BOUNDARY_TYPE(IWM1)==SOLID_BOUNDARY) THEN
               DR%IOR=-1
               HIT_SOLID = .TRUE.
               STEP_FRACTION = (X(II)-X_OLD-0.05_EB*DX(II))/(DR%X-X_OLD)
            ENDIF
            IF (IIN<II .AND. BOUNDARY_TYPE(IWP1)==SOLID_BOUNDARY) THEN
               DR%IOR= 1
               HIT_SOLID = .TRUE.
               STEP_FRACTION = (X(IIN)-X_OLD+0.05_EB*DX(IIN))/(DR%X-X_OLD)
            ENDIF
            IF (JJN>JJ .AND. BOUNDARY_TYPE(IWM2)==SOLID_BOUNDARY) THEN
               DR%IOR=-2
               HIT_SOLID = .TRUE.
               STEP_FRACTION = (Y(JJ)-Y_OLD-0.05_EB*DY(JJ))/(DR%Y-Y_OLD)
            ENDIF
            IF (JJN<JJ .AND. BOUNDARY_TYPE(IWP2)==SOLID_BOUNDARY) THEN
               DR%IOR= 2
               HIT_SOLID = .TRUE.
               STEP_FRACTION = (Y(JJN)-Y_OLD+0.05_EB*DY(JJN))/(DR%Y-Y_OLD)
            ENDIF
            DR%WALL_INDEX = WALL_INDEX(ICN,DR%IOR)
         ENDIF
   
         ! Check if droplet has crossed no solid planes or too many
  
         IF_HIT_SOLID: IF (HIT_SOLID) THEN
 
            IF (DR%WALL_INDEX==0) CYCLE SUB_TIME_STEP_ITERATIONS

            ! Move particle to where it almost hits solid

            DR%X = X_OLD + STEP_FRACTION*DTSP*DR%U
            DR%Y = Y_OLD + STEP_FRACTION*DTSP*DR%V
            DR%Z = Z_OLD + STEP_FRACTION*DTSP*DR%W

            SELECT CASE(ABS(DR%IOR))
               CASE(1)
                  XI = CELLSI(FLOOR((DR%X-XS)*RDXINT))
                  IIN = FLOOR(XI+1._EB)
               CASE(2)
                  YJ = CELLSJ(FLOOR((DR%Y-YS)*RDYINT))
                  JJN = FLOOR(YJ+1._EB)
               CASE(3)
                  ZK = CELLSK(FLOOR((DR%Z-ZS)*RDZINT))
                  KKN = FLOOR(ZK+1._EB)
            END SELECT

            ICN = CELL_INDEX(IIN,JJN,KKN)

            IF (IOR_OLD==DR%IOR) CYCLE SUB_TIME_STEP_ITERATIONS

            ! Choose a direction for the droplets to move

            DIRECTION: SELECT CASE(ABS(DR%IOR))
               CASE (1:2) DIRECTION  
                  DR%U = 0._EB
                  DR%V = 0._EB
                  DR%W = -PC%VERTICAL_VELOCITY 
               CASE(3) DIRECTION
                  IBC = IJKW(5,DR%WALL_INDEX)
                  CALL RANDOM_NUMBER(RN)
                  THETA_RN = TWOPI*RN
                  DR%U = PC%HORIZONTAL_VELOCITY*COS(THETA_RN)
                  DR%V = PC%HORIZONTAL_VELOCITY*SIN(THETA_RN)
                  DR%W = 0._EB
            END SELECT DIRECTION

         ENDIF IF_HIT_SOLID
      ENDIF AIR_TO_SOLID 

      ! Check if attached droplets are still attached

      IW = 0
      SELECT CASE(DR%IOR)
         CASE( 1)
            IW = WALL_INDEX(ICN,-1)
            IF (BOUNDARY_TYPE(IW)/=SOLID_BOUNDARY) DR%U = -PC%HORIZONTAL_VELOCITY
         CASE(-1)
            IW = WALL_INDEX(ICN, 1)
            IF (BOUNDARY_TYPE(IW)/=SOLID_BOUNDARY) DR%U =  PC%HORIZONTAL_VELOCITY
         CASE( 2)
            IW = WALL_INDEX(ICN,-2)
            IF (BOUNDARY_TYPE(IW)/=SOLID_BOUNDARY) DR%V = -PC%HORIZONTAL_VELOCITY
         CASE(-2)
            IW = WALL_INDEX(ICN, 2)
            IF (BOUNDARY_TYPE(IW)/=SOLID_BOUNDARY) DR%V =  PC%HORIZONTAL_VELOCITY
         CASE( 3)
            IW = WALL_INDEX(ICN,-3)
            IF (BOUNDARY_TYPE(IW)/=SOLID_BOUNDARY) THEN
               DR%U = -DR%U
               DR%V = -DR%V
               DR%Z = DR%Z - 0.2_EB*DZ(KK)
            ENDIF
         CASE(-3)
            IF (.NOT.SOLID(ICN)) DR%IOR = 0
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

REAL(EB), POINTER, DIMENSION(:,:,:) :: DROP_DEN,DROP_RAD,DROP_TMP
REAL(EB), POINTER, DIMENSION(:) :: WMPUAOLD,WCPUAOLD
REAL(EB) :: RD,RDS,RDC,RRD,RD_NEW, &
            NUSSELT,KA,H_V, &
            RVC,WGT,OMWGT, &
            QUSE,HD, &
            SH_FAC,NU_FAC,T,PR,MVAP,XI,YJ,ZK,RDT,MU_AIR, &
            HS,HLF,DTLF,Y_IGAS,QRADD,QWALLC,QWALLR,DEN_ADD, &
            Y_IGASN,Y_EQUIL,Y_GAS,MVAPMIN,MVAPMAX
INTEGER :: I,II,JJ,KK,IW,IGAS,N_PC,ITER,EVAP_INDEX
REAL(EB) :: SC_AIR,D_AIR,DHOR,SHER,XVAP,QVAP, HADT_WALL, HADT_G, &
            MDROP,RHO_G,MW_RATIO,MW_DROP,FTPR,&
            CP_IGAS,CP_DROP,MGAS,ATOP,ABOT,CR2,WALL,FILMD,Z_2,Y_IGAS_Z,Z_F,&
            TMP_GAS_OLD,TMP_GAS_NEW,TMP_DROP_OLD,TMP_DROP_NEW,TMP_MELT,TMP_BOIL,TMP_WALL, &
            DENOM,EGASOLD,ESYSOLD,MGCP,MVCP,MDCPNEW
INTEGER, INTENT(IN) :: NM


D_VAP = 0.

HS     = 3000._EB    ! Heat transfer coef from solid to drop (W/m2/K)
HLF    = 1000._EB    ! Leidenfrost heat transfer coef (W/m2/K)
PR     = 0.7_EB      ! Prandtl number of air
FILMD  = 0.002_EB    ! Thickness of continous water film
NU_FAC = 0.6_EB*PR**ONTH
RDT    = 1._EB/DT
CR2    = 2._EB**ONTH

D_AIR  = 3.33E-5_EB
SC_AIR = 0.6_EB
SH_FAC = 0.6_EB*SC_AIR**ONTH

DROP_DEN => WORK4
DROP_RAD => WORK5
DROP_TMP => WORK6
WMPUAOLD => WALL_WORK1
WCPUAOLD => WALL_WORK2

EVAP_INDEX_LOOP: DO EVAP_INDEX = 1,N_EVAP_INDICIES
   WMPUAOLD = WMPUA(:,EVAP_INDEX)
   WCPUAOLD = WCPUA(:,EVAP_INDEX)
   WMPUA(:,EVAP_INDEX) = 0._EB
   WCPUA(:,EVAP_INDEX) = 0._EB
   PC_LOOP: DO N_PC = 1,N_PART
      PC => PARTICLE_CLASS(N_PC)
      IF (PC%EVAP_INDEX/=EVAP_INDEX) CYCLE PC_LOOP
      DROP_DEN = 0._EB
      DROP_TMP = 0._EB
      DROP_RAD = 0._EB
      CP_DROP =  PC%C_P
      FTPR     = PC%FTPR
      TMP_MELT = PC%TMP_MELT
      TMP_BOIL = PC%TMP_V
      IGAS     = PC%SPECIES_INDEX
      MW_DROP  = SPECIES(IGAS)%MW
      MW_RATIO = SPECIES(0)%MW/MW_DROP
      DTLF     = 2._EB*(TMP_BOIL-TMPM)
      CP_IGAS = PC%GAMMA_VAP*R0/((PC%GAMMA_VAP-1._EB)*MW_DROP)
      ! Heat of Vaporization
   !   IF (IGAS==I_WATER) THEN 
   !      H_V = PC%H_V + 2400._EB*(TMP_BOIL-TMP_DROP_OLD) ! Incropera, 4th, p 846
   !   ELSE
         H_V = PC%H_V
   !   ENDIF

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
         
         !Intiailize droplet thermophysical data
         RD   = DR%R
         RDS  = RD*RD
         RDC  = RD*RDS
         RRD  = 1._EB/RD
         TMP_DROP_OLD = DR%TMP
         TMP_DROP_NEW = TMP_DROP_OLD
         MDROP = FTPR*RDC       ! Mass of drop
         WGT   = DR%PWT
         DHOR  = H_V*MW_DROP/R0
         WALL  = 1._EB
         ABOT  = 0._EB
         ATOP  = 4._EB*PI*RDS
         QUSE  = 0._EB
         QWALLR= 0._EB
         QWALLC= 0._EB
         QVAP  = 0._EB
         MVAP  = 0._EB
         MVAPMIN = 0._EB
         MVAPMAX = MDROP
         RD_NEW = RD
               
         ! Gas conditions
         TMP_GAS_OLD  = TMP(II,JJ,KK)
         TMP_GAS_NEW  = TMP_GAS_OLD
         RHO_G  = RHO(II,JJ,KK)
         MU_AIR = SPECIES(0)%MU(MIN(500,NINT(0.1_EB*TMP_GAS_OLD)))
         RVC    = RDX(II)*RDY(JJ)*RDZ(KK)
         MGAS    = RHO_G/RVC           ! Mass of gas cell
         KA      = CPOPR*MU_AIR
         Y_IGAS     = YY(II,JJ,KK,IGAS)
         Y_IGAS_Z   = 0._EB
         IF (MIXTURE_FRACTION .AND. IGAS==I_WATER) THEN ! Add water vapor from combustion
            IF (CO_PRODUCTION) THEN
               Z_2 = YY(II,JJ,KK,I_PROG_CO)
            ELSE
               Z_2 = 0._EB
            ENDIF
            CALL GET_MASS_FRACTION(YY(II,JJ,KK,I_FUEL),Z_2,YY(II,JJ,KK,I_PROG_F),H2O_INDEX,Y_SUM(II,JJ,KK),Z_F,Y_IGAS_Z)
         ENDIF
         Y_GAS = Y_IGAS + Y_IGAS_Z

         ! Determine radiation absorption for droplet (J)
         QRADD = 0._EB
         IF (AVG_DROP_DEN(II,JJ,KK,EVAP_INDEX )>0._EB) QRADD = QR_W(II,JJ,KK) / SUM(AVG_DROP_DEN(II,JJ,KK,:)) * DT * MDROP

         ! Set variables for heat transfer on solid
         SET_SOLID: IF (DR%IOR/=0 .AND. DR%WALL_INDEX>0) THEN
            IW   = DR%WALL_INDEX
            WALL = 0.5_EB  ! Droplet dimension adjust assuming hemisphere
            RD   = CR2*RD
            RDS  = RD*RD
            RDC  = RD*RDS
            RRD  = 1._EB/RD
            ABOT = PI*RDS
            ATOP = 2._EB*PI*RDS
            IF (WMPUAOLD(IW)/PC%DENSITY>FILMD) THEN
               ABOT = MDROP/WMPUAOLD(IW)
               ATOP = ABOT
            ENDIF
            IF (WGT*ABOT > AW(IW)) ABOT = AW(IW)/WGT
            IF (TMP_WALL-TMP_BOIL>DTLF) THEN  ! Leidenfrost droplet
               HADT_WALL = ABOT*DT*HLF
            ELSE
               HADT_WALL = ABOT*DT*HS
            ENDIF
            TMP_WALL   = TMP_F(IW)
            DR%RE  = RHO_G*PC%HORIZONTAL_VELOCITY*2._EB*RD/MU_AIR
            QRADD  = QRADD + ABOT*DT*QRAD(IW)
            QWALLC = HADT_WALL*(TMP_WALL-TMP_DROP_OLD)
            IF (SURFACE(IJKW(5,IW))%THERMAL_BC_INDEX == THERMALLY_THICK) THEN
              QUSE = AW(IW)*CP_DROP*MDROP*RCP_W(IW)*(TMP_WALL-TMP_DROP_OLD)/ &
                     (AW(IW)*RCP_W(IW)+CP_DROP*MDROP*WGT)
              IF(ABS(QWALLC)>ABS(QUSE)) QWALLC=QUSE
            ENDIF
         ELSE !If not a wall set wall variables just in case
            QWALLC = 0._EB
            ABOT  = ATOP
         ENDIF SET_SOLID  
         
         !Heat Transfer coefficient
         NUSSELT = 2._EB + NU_FAC*SQRT(DR%RE)
         SHER    = 2._EB + SH_FAC*SQRT(DR%RE)
         HD      = 0.5_EB*NUSSELT*KA*RRD
         HADT_G  = HD*ATOP*DT
         
         ! Start the mass/energy transfer calc
         !Compute evaporation
         XVAP  = MIN(1._EB,EXP(DHOR*(1._EB/TMP_BOIL-1._EB/TMP_DROP_OLD)))
         Y_EQUIL  = XVAP/(MW_RATIO + (1._EB-MW_RATIO)*XVAP)
         IF (Y_EQUIL > Y_GAS) THEN
            IF (Y_EQUIL==1._EB) THEN   
               QUSE = QRADD + QWALLC + HADT_G * (TMP_GAS_OLD - TMP_DROP_OLD) + MDROP * CP_DROP * (TMP_DROP_OLD - TMP_MELT)        
               MVAP = QUSE/H_V
            ELSE
               MVAP = SHER*RHO_G*D_AIR*WALL*TWOPI*RD*LOG((Y_GAS-1._EB)/(Y_EQUIL-1._EB))* DT
            ENDIF
            MVAP = MAX(0._EB,MIN(MVAP,MDROP))
         ENDIF
         
         !Misc variables
         MGCP    = MGAS * CP_GAMMA
         MVCP    = MVAP * CP_IGAS * 0.5_EB
         MDCPNEW = (MDROP - MVAP) * CP_DROP
         
         !Compute new drop temperature and energy 
         !Equations from solving for TMP_GAS and TMP_DROP using conservation equations for all mass and for gas mass only
         ESYSOLD = MGCP * TMP_GAS_OLD + WGT * (CP_DROP * MDROP * (TMP_DROP_OLD - TMP_MELT) + QRADD + QWALLC) - &
                  WGT * (MDCPNEW * (-TMP_MELT) + H_V * MVAP + MVCP * (TMP_GAS_OLD - TMP_DROP_OLD))
         EGASOLD = MGCP * TMP_GAS_OLD  -WGT * ((0.5_EB*HADT_G + MVCP) * (TMP_GAS_OLD - TMP_DROP_OLD))

         DENOM = HADT_G * (MGCP + MDCPNEW * WGT) + 2._EB * MDCPNEW * (MGCP + MVCP * WGT)
         TMP_DROP_NEW = (-2._EB * EGASOLD * (MGCP + MVCP * WGT) + ESYSOLD * (2._EB * MGCP + (HADT_G + 2._EB * MVCP) * WGT)) / &
                        (WGT * DENOM)
         TMP_GAS_NEW =  (2._EB * EGASOLD * (MDCPNEW - MVCP) + ESYSOLD * (HADT_G + 2._EB * MVCP)) / DENOM

         !If initial MVAP was zero or new drop temperature is over boiling compute a new guess
         IF (MVAP ==0._EB .OR. TMP_DROP_NEW > TMP_BOIL) THEN
            XVAP  = MIN(1._EB,EXP(DHOR*(1._EB/TMP_BOIL-1._EB/MIN(MAX(TMP_MELT,TMP_DROP_NEW),TMP_BOIL))))
            Y_EQUIL  = XVAP/(MW_RATIO + (1._EB-MW_RATIO)*XVAP)
            IF (Y_EQUIL > Y_GAS) THEN
               IF (Y_EQUIL==1._EB) THEN      
                  QUSE = QRADD + QWALLC + HADT_G * (TMP_GAS_OLD - TMP_DROP_OLD) + MDROP * CP_DROP * (TMP_DROP_OLD - TMP_MELT)
                  MVAP = QUSE/H_V
               ELSE
                  MVAP = SHER*RHO_G*D_AIR*WALL*TWOPI*RD*LOG((Y_GAS-1._EB)/(Y_EQUIL-1._EB))* DT
               ENDIF

               MVAP = MAX(0._EB,MIN(MVAP,MDROP))
               
               MVCP    = MVAP * CP_IGAS * 0.5_EB
               MDCPNEW = (MDROP - MVAP) * CP_DROP

               ESYSOLD = MGCP * TMP_GAS_OLD + WGT * (CP_DROP * MDROP * (TMP_DROP_OLD - TMP_MELT) + QRADD + QWALLC) - &
                        WGT * (MDCPNEW * (-TMP_MELT) + H_V * MVAP + MVCP * (TMP_GAS_OLD - TMP_DROP_OLD))
               EGASOLD = MGCP * TMP_GAS_OLD  -WGT * ((0.5_EB*HADT_G + MVCP) * (TMP_GAS_OLD - TMP_DROP_OLD))

               DENOM = HADT_G * (MGCP + MDCPNEW * WGT) + 2._EB * MDCPNEW * (MGCP + MVCP * WGT)
               TMP_DROP_NEW = (-2._EB * EGASOLD * (MGCP + MVCP * WGT) + ESYSOLD * (2._EB * MGCP + (HADT_G + 2._EB * MVCP) * WGT))&
                              /(WGT * DENOM)
               TMP_GAS_NEW =  (2._EB * EGASOLD * (MDCPNEW - MVCP) + ESYSOLD * (HADT_G + 2._EB * MVCP)) / DENOM
            ENDIF  
         ENDIF

         !Using MVAP guess as max value, iterate to get actual evaporation
         ITERATE_IF: IF (MVAP > 0._EB) THEN
            ITER = 0
            MVAPMAX = MVAP
            MVAPMIN = 0._EB
            ITER_LOOP: DO WHILE (ITER <=10)
               XVAP  = MIN(1._EB,EXP(DHOR*(1._EB/TMP_BOIL-1._EB/MIN(MAX(TMP_MELT,TMP_DROP_NEW),TMP_BOIL))))
               Y_EQUIL  = XVAP/(MW_RATIO + (1._EB-MW_RATIO)*XVAP)
               Y_IGASN = MIN(1._EB,(MVAP*WGT+MGAS*Y_IGAS)/(MVAP*WGT+MGAS))
               Y_GAS = MIN(1._EB,Y_IGASN + Y_IGAS_Z)
               IF (TMP_DROP_OLD < TMP_MELT) THEN
                  ! Don't allow drop to freeze, reduce evaporation
                  MVAPMAX = MVAP
                  MVAP =  0.5_EB * (MVAP + MVAPMIN)
               ELSE
                  IF (TMP_GAS_NEW < TMP_GAS_OLD .AND. TMP_DROP_NEW > TMP_GAS_NEW)  THEN
                     !Don't allow too much gas to drop energy transfer, reduce evaporation
                     MVAPMAX = MVAP
                     MVAP =  0.5_EB * (MVAP + MVAPMIN)
                  ELSE
                     IF (Y_EQUIL < Y_GAS) THEN
                        MVAPMAX = MVAP
                        MVAP =  0.5_EB * (MVAP + MVAPMIN)
                     ELSEIF (Y_EQUIL > Y_GAS)  THEN
                        IF (MVAP == MDROP) EXIT ITER_LOOP !If all the drop is evaporated, we can't evaporate any more
                        MVAPMIN = MVAP
                        MVAP = 0.5_EB * (MVAP + MVAPMAX)
                     ELSE !Got lucky and Y_EQUIL = Y_GAS with good temperatures
                        EXIT ITER_LOOP
                     ENDIF
                  ENDIF
               ENDIF
               ITER = ITER + 1

               MVCP    = MVAP * CP_IGAS * 0.5_EB
               MDCPNEW = (MDROP - MVAP) * CP_DROP
               
               ESYSOLD = MGCP * TMP_GAS_OLD + WGT * (CP_DROP * MDROP * (TMP_DROP_OLD - TMP_MELT) + QRADD + QWALLC) - &
                        WGT * (MDCPNEW * (-TMP_MELT) + H_V * MVAP + MVCP * (TMP_GAS_OLD - TMP_DROP_OLD))
               EGASOLD = MGCP * TMP_GAS_OLD  -WGT * ((0.5_EB*HADT_G + MVCP) * (TMP_GAS_OLD - TMP_DROP_OLD))

               DENOM = HADT_G * (MGCP + MDCPNEW * WGT) + 2._EB * MDCPNEW * (MGCP + MVCP * WGT)
               TMP_DROP_NEW = (-2._EB * EGASOLD * (MGCP + MVCP * WGT) + ESYSOLD * (2._EB * MGCP + (HADT_G + 2._EB * MVCP) * WGT))&
                               /(WGT * DENOM)
               TMP_GAS_NEW =  (2._EB * EGASOLD * (MDCPNEW - MVCP) + ESYSOLD * (HADT_G + 2._EB * MVCP)) / DENOM
            END DO ITER_LOOP
         END IF ITERATE_IF

         IF (TMP_DROP_NEW > TMP_BOIL) TMP_DROP_NEW = TMP_BOIL
         IF (TMP_DROP_NEW < TMP_MELT) TMP_DROP_NEW = TMP_MELT
       
         MDROP = MDROP - MVAP
         MVAP  = MVAP*WGT
         DR%R = (MDROP / FTPR) ** ONTH
        
         IF (WALL<1._EB) THEN
 !           QWALLC = HADT_WALL * (TMP_WALL - TMP_DROP_OLD)      
            QUSE = QWALLR+QWALLC
            RDC       = RD_NEW**3
            WCPUA(IW,EVAP_INDEX ) = WCPUA(IW,EVAP_INDEX ) + WGT*RDT*QUSE*RAW(IW)
            WMPUA(IW,EVAP_INDEX ) = WMPUA(IW,EVAP_INDEX ) + WGT*MDROP*RAW(IW)
         ENDIF

         ! Decrease temperature due to droplet heating/vaporization
         TMP(II,JJ,KK)=TMP_GAS_NEW
         TMP(II,JJ,KK) = MIN(TMPMAX,MAX(TMPMIN,TMP(II,JJ,KK)))

         ! Save water vapor production rate for diverence expression
         D_VAP(II,JJ,KK) = D_VAP(II,JJ,KK) +R0*RDT*( &
               RHO_G*((1-Y_IGAS)/SPECIES(0)%MW+Y_IGAS/MW_DROP)*(TMP(II,JJ,KK)-TMP_GAS_OLD) +RVC*MVAP*TMP(II,JJ,KK)/MW_DROP) 

         ! Adjust mass of evaporated liquid to account for different Heat of Combustion between fuel droplet and gas
         MVAP = PC%ADJUST_EVAPORATION*MVAP

         ! Add water vapor or fuel gas to the cell

         YY(II,JJ,KK,IGAS)= MIN(1._EB,(MVAP+MGAS*Y_IGAS)/(MVAP+MGAS))
         ! Add new mass from vaporized water droplet 

         RHO(II,JJ,KK) = RHO(II,JJ,KK) + MVAP*RVC

         DR%TMP = TMP_DROP_NEW
         IF (DR%R<=0._EB) CYCLE DROPLET_LOOP

         ! Assign water mass to the cell for airborne drops

         IF (DR%IOR==0) THEN
            DEN_ADD = WGT*MDROP*RVC
            DROP_DEN(II,JJ,KK) = DROP_DEN(II,JJ,KK) + DEN_ADD
            DROP_TMP(II,JJ,KK) = DROP_TMP(II,JJ,KK) + DEN_ADD*TMP_DROP_NEW
            DROP_RAD(II,JJ,KK) = DROP_RAD(II,JJ,KK) + DEN_ADD*DR%R
         ENDIF
         
      ENDDO DROPLET_LOOP

      ! Calculate cumulative quantities for droplet "clouds"

      WGT   = .5_EB
      OMWGT = 1._EB-WGT
    
      DROP_RAD = DROP_RAD/(DROP_DEN+1.E-10_EB)
      DROP_TMP = DROP_TMP/(DROP_DEN+1.E-10_EB)
    
      AVG_DROP_RAD(:,:,:,EVAP_INDEX ) = (AVG_DROP_DEN(:,:,:,EVAP_INDEX )*AVG_DROP_RAD(:,:,:,EVAP_INDEX )+ DROP_DEN*DROP_RAD) &
                                          /(AVG_DROP_DEN(:,:,:,EVAP_INDEX ) + DROP_DEN + 1.0E-10_EB)
      AVG_DROP_TMP(:,:,:,EVAP_INDEX ) = (AVG_DROP_DEN(:,:,:,EVAP_INDEX )*AVG_DROP_TMP(:,:,:,EVAP_INDEX )+DROP_DEN*DROP_TMP) &
                                          /(AVG_DROP_DEN(:,:,:,EVAP_INDEX ) + DROP_DEN + 1.0E-10_EB)
      AVG_DROP_TMP(:,:,:,EVAP_INDEX ) = MAX(TMPM,AVG_DROP_TMP(:,:,:,EVAP_INDEX ))
      AVG_DROP_DEN(:,:,:,EVAP_INDEX ) = WGT*AVG_DROP_DEN(:,:,:,EVAP_INDEX ) + OMWGT*DROP_DEN
      WHERE (AVG_DROP_DEN(:,:,:,EVAP_INDEX )<0.0001_EB .AND. DROP_DEN==0._EB) AVG_DROP_DEN(:,:,:,EVAP_INDEX ) = 0.0_EB

      ! Weight the new water mass array

      WMPUA(:,EVAP_INDEX ) = WGT*WMPUAOLD + OMWGT*WMPUA(:,EVAP_INDEX )
      WCPUA(:,EVAP_INDEX ) = WGT*WCPUAOLD + OMWGT*WCPUA(:,EVAP_INDEX)

   ENDDO PC_LOOP
ENDDO EVAP_INDEX_LOOP

! Remove out-of-bounds particles
 
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
REAL(EB) :: T
REAL(EB), PARAMETER :: RDMIN=1.0E-6_EB   ! Minimum droplet size (m)
 
IKILL = 0
DROP_LOOP: DO I=1,NLP
   WEED_LOOP: DO
      DR=>DROPLET(I) 
      IPC=DR%CLASS 
      PC=>PARTICLE_CLASS(IPC)
      IF (I>NLP-IKILL) EXIT DROP_LOOP
      IF (DR%R<=RDMIN .AND. .NOT.PC%MASSLESS) THEN 
         CALL REPLACE
         CYCLE WEED_LOOP
      ENDIF
      IF (T-DR%T>PC%LIFETIME) THEN
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
   IF (M2%N_DROP_ORPHANS>M2%N_DROP_ORPHANS_DIM) CALL RE_ALLOCATE_DROPLETS(2,NM,NOM)
   M2%DROPLET(M2%N_DROP_ORPHANS) = DROPLET(I)
ENDIF
 
DROPLET(I) = DROPLET(NLP-IKILL)
IKILL = IKILL + 1
 
END SUBROUTINE REPLACE
 
END SUBROUTINE REMOVE_DROPLETS
 

END MODULE PART
