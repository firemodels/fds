!> \brief Variables and subroutines pertaining to Lagrangian particles

MODULE PART

USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_POINTERS
USE COMP_FUNCTIONS, ONLY : CURRENT_TIME
USE COMPLEX_GEOMETRY, ONLY: CC_IDCC,CC_IDCF,CC_CGSC,CC_SOLID
IMPLICIT NONE (TYPE,EXTERNAL)

PRIVATE

PUBLIC INSERT_ALL_PARTICLES,MOVE_PARTICLES,PARTICLE_MASS_ENERGY_TRANSFER,REMOVE_PARTICLES,&
       GENERATE_PARTICLE_DISTRIBUTIONS,PARTICLE_MOMENTUM_TRANSFER,GET_RVC

CONTAINS


!> \brief Generate size distribution statistics for all particles classes

SUBROUTINE GENERATE_PARTICLE_DISTRIBUTIONS

USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
USE PHYSICAL_FUNCTIONS, ONLY : PARTICLE_SIZE_DISTRIBUTION

REAL(EB) :: LL,UL,BIN_SIZE,TNOW,DD,DI
INTEGER  :: I,J,IL,IU,ILPC
TYPE (LAGRANGIAN_PARTICLE_CLASS_TYPE), POINTER :: LPC
TYPE (RAMPS_TYPE), POINTER :: RM

IF (N_LAGRANGIAN_CLASSES==0) RETURN ! Don't waste time if no particles

TNOW=CURRENT_TIME()

PART_CLASS_LOOP: DO ILPC=1,N_LAGRANGIAN_CLASSES

   LPC=>LAGRANGIAN_PARTICLE_CLASS(ILPC)

   IF_SIZE_DISTRIBUTION: IF (.NOT.LPC%MONODISPERSE .AND. (LPC%DIAMETER>0._EB) .OR. LPC%CNF_RAMP_INDEX>=0) THEN
      IF (LPC%CNF_RAMP_INDEX<0) THEN
         CALL PARTICLE_SIZE_DISTRIBUTION(LPC%DIAMETER,LPC%R_CNF(:),LPC%CNF(:),LPC%CVF(:),NDC,LPC%GAMMA,LPC%SIGMA,LPC%DISTRIBUTION)
      ELSE
         RM=>RAMPS(LPC%CNF_RAMP_INDEX)
         LPC%DIAMETER = RM%SPAN*0.5_EB ! Not used for MONODISPERSE
         DD           = RM%SPAN/NDC
         LPC%R_CNF(0) = RM%T_MIN
         LPC%CNF(0)   = 0._EB
         LPC%CVF(0)   = 0._EB
         DO I=1,NDC
            DI           = RM%T_MIN + I*DD
            LPC%R_CNF(I) = DI
            LPC%CNF(I)   = EVALUATE_RAMP(DI,LPC%CNF_RAMP_INDEX)
            LPC%CVF(I)   = LPC%CVF(I-1) + (DI-0.5_EB*DD)**3*(LPC%CNF(I)-LPC%CNF(I-1))
         ENDDO
         LPC%R_CNF = 1.E-6_EB*0.5_EB*LPC%R_CNF ! Convert diameter in microns to radius in meters.
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
               LPC%STRATUM_INDEX_LOWER(I) = J-1
               EXIT LL_LOOP
            ENDIF
         ENDDO LL_LOOP
         IU = NDC
         UL_LOOP: DO J=NDC,1,-1
            IF (LPC%R_CNF(J)<=UL) THEN
               IU = J
               LPC%STRATUM_INDEX_UPPER(I) = J
               EXIT UL_LOOP
            ENDIF
         ENDDO UL_LOOP
         LPC%W_CNF(I) = LPC%CNF(IU) - LPC%CNF(IL)
      ENDDO STRATIFY

      ! Compute mean droplet volume for distribution

      LPC%MEAN_DROPLET_VOLUME = 0._EB
      DO I=1,NDC
         LPC%MEAN_DROPLET_VOLUME = LPC%MEAN_DROPLET_VOLUME + ( LPC%CNF(I) - LPC%CNF(I-1) ) * FOTHPI*LPC%R_CNF(I-1)**3
      ENDDO

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
               LPC%BREAKUP_CNF(I)   = EVALUATE_RAMP(DI,LPC%BREAKUP_CNF_RAMP_INDEX)
         ENDDO
         LPC%BREAKUP_CNF=LPC%BREAKUP_CNF/LPC%BREAKUP_CNF(NDC)
      ENDIF
   ENDIF

ENDDO PART_CLASS_LOOP

T_USED(8)=T_USED(8)+CURRENT_TIME()-TNOW
END SUBROUTINE GENERATE_PARTICLE_DISTRIBUTIONS


!> \brief Insert Lagrangian particles into the domain every time step

SUBROUTINE INSERT_ALL_PARTICLES(T,NM)

USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP,RANDOM_CHOICE
USE GEOMETRY_FUNCTIONS, ONLY: RANDOM_RECTANGLE,RANDOM_CONE,RANDOM_RING,CONE_MESH_INTERSECTION_VOLUME,UNIFORM_RING,&
                              RING_MESH_INTERSECTION_ARC
USE COMP_FUNCTIONS, ONLY: GET_FILE_NUMBER, SHUTDOWN
USE TRAN, ONLY: GET_IJK
USE DEVICE_VARIABLES
USE CONTROL_VARIABLES
REAL(EB), INTENT(IN) :: T
INTEGER, INTENT(IN) :: NM
REAL     :: RN,RN2,RN3,INIT_CALLS
REAL(EB) :: PHI_RN,FLOW_RATE,THETA_RN,SPHI,CPHI,MASS_SUM,D_PRES_FACTOR, &
            STHETA,CTHETA,PWT0,PARTICLE_SPEED,SHIFT1,SHIFT2,XTMP,YTMP,ZTMP,VLEN, &
            TRIGT1,TRIGT2,TNOW,TSI,PIPE_PRESSURE,X1,X2,Y1,Y2,Z1,Z2, &
            ETA,ETA_MAX,ETA_MIN,XI,YJ,ZK,X_OFFSET,Y_OFFSET,Z_OFFSET, &
            VXMIN,VXMAX,VYMIN,VYMAX,VZMIN,VZMAX,VDX,VDY,VDZ
REAL(EB), PARAMETER :: VENT_OFFSET=0.5
INTEGER :: IP,KS,II,JJ,KK,IC,IL,IU,ILPC,DROP_SUM,IIG,JJG,KKG,IW,IOR,STRATUM,INIT_INDEX,ICF,ITS,ITERATIONS
INTEGER :: N,N_INSERT,ILAT,NEW_LP_INDEX,III,IOS,LU_VEG_IN,NVOX,VXI
INTEGER, ALLOCATABLE, DIMENSION(:) :: LP_INDEX_LOOKUP
LOGICAL :: INSERT_ANOTHER_BATCH
TYPE (PROPERTY_TYPE), POINTER :: PY
TYPE (TABLES_TYPE), POINTER :: TA
TYPE (DEVICE_TYPE), POINTER :: DV
TYPE (SURFACE_TYPE), POINTER :: SF
TYPE (LAGRANGIAN_PARTICLE_TYPE), POINTER :: LP
TYPE (LAGRANGIAN_PARTICLE_CLASS_TYPE), POINTER :: LPC
TYPE (INITIALIZATION_TYPE), POINTER :: IN
TYPE (WALL_TYPE), POINTER :: WC
TYPE (BOUNDARY_COORD_TYPE), POINTER :: BC
TYPE (BOUNDARY_ONE_D_TYPE), POINTER :: ONE_D
TYPE (BOUNDARY_PROP1_TYPE), POINTER :: B1
CHARACTER(MESSAGE_LENGTH) :: MESSAGE

IF (N_LAGRANGIAN_CLASSES==0) RETURN ! Don't waste time if no particles

TNOW=CURRENT_TIME()
CALL POINT_TO_MESH(NM)

! Insert particles at spray nozzles (INSERT_SPRAY_PARTICLES), VENT surfaces (INSERT_VENT_PARTICLES), specified
! volumes (INSERT_VOLUMETRIC_PARTICLES), or HVAC ducts (INSERT_DUCT_PARTICLES). If the insertions are to be made according
! to a user-specified clock, allow the possibility via the OVERALL_INSERT_LOOP of inserting multiple sets of particles.

OVERALL_INSERT_LOOP: DO

   INSERT_ANOTHER_BATCH = .FALSE.

   ! Insert spray particles/droplets

   CALL INSERT_SPRAY_PARTICLES

   ! Insert particles at a solid boundary or vent.

   CALL INSERT_VENT_PARTICLES

   ! Insert particles within a volume specified by an INIT line.

   DO INIT_INDEX=1,N_INIT
      IN => INITIALIZATION(INIT_INDEX)
      IF (IN%INVOKED_BY_SURF) CYCLE
      X_OFFSET = 0._EB
      Y_OFFSET = 0._EB
      Z_OFFSET = 0._EB
      ! Check if INIT line links to a binary file of bulk density data
      IF (IN%BULK_DENSITY_FILE=='null') THEN
         CALL INSERT_VOLUMETRIC_PARTICLES
      ELSE
         IF (.NOT. IN%ALREADY_INSERTED(NM)) THEN
            ! Build out of blocks
            IN%SHAPE = 'BLOCK'
            LU_VEG_IN = GET_FILE_NUMBER()
            OPEN(UNIT=LU_VEG_IN,FILE=TRIM(IN%BULK_DENSITY_FILE),&
               STATUS='OLD',FORM='UNFORMATTED',ACTION='READ',IOSTAT=IOS)
            IF (IOS==0) THEN
               READ(LU_VEG_IN) VXMIN,VXMAX,VYMIN,VYMAX,VZMIN,VZMAX
               ! Skip if volume containing vegetation is entirely outside the current mesh
               IF (VXMIN>XF .OR. VXMAX<XS .OR. VYMIN>YF .OR. VYMAX<YS .OR. VZMIN>ZF .OR. VZMAX<ZS) THEN
                  IN%ALREADY_INSERTED(NM)=.TRUE.
                  CYCLE
               ENDIF
               ! Voxel resolution
               READ(LU_VEG_IN) VDX,VDY,VDZ
               ! Number of vegetation containing voxels
               READ(LU_VEG_IN) NVOX
               ! Create a pseudo-INIT block for each voxel
               VEG_INSERT_LOOP: DO VXI=1,NVOX
                  CALL INSERT_VOLUMETRIC_PARTICLES
               ENDDO VEG_INSERT_LOOP
               CLOSE(LU_VEG_IN)
               IN%ALREADY_INSERTED(NM) = .TRUE.
            ELSE
               WRITE(MESSAGE,'(A,I0,A,A,A)') 'ERROR: INIT ',INIT_INDEX,', could not read binary bulk density file ', &
                                         TRIM(IN%BULK_DENSITY_FILE),'. Check file exists.'
               CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.); RETURN
            ENDIF
         ENDIF
      ENDIF
   ENDDO

   ! Insert volumetric particles that are invoked by a particular SURF type.

   IF (INIT_INVOKED_BY_SURF) THEN

      DO IW=1,N_INTERNAL_WALL_CELLS+N_EXTERNAL_WALL_CELLS
         WC => MESHES(NM)%WALL(IW)
         BC => MESHES(NM)%BOUNDARY_COORD(WC%BC_INDEX)
         SF => SURFACE(WC%SURF_INDEX)
         IF (BC%IOR==3 .AND. SF%INIT_INDICES(1)>0) THEN
            CALL RANDOM_NUMBER(RN)
            X_OFFSET = BC%X + (RN-0.5_FB)*DX(BC%IIG)
            CALL RANDOM_NUMBER(RN)
            Y_OFFSET = BC%Y + (RN-0.5_FB)*DY(BC%JJG)
            Z_OFFSET = BC%Z
            B1 => MESHES(NM)%BOUNDARY_PROP1(WC%B1_INDEX)
            INIT_CALLS = REAL(SF%INIT_PER_AREA*B1%AREA)
            ITERATIONS = INT(INIT_CALLS)
            CALL RANDOM_NUMBER(RN)
            IF (RN<INIT_CALLS-ITERATIONS) ITERATIONS = ITERATIONS + 1
            DO ITS=1,ITERATIONS
               DO III=1,10
                  IF (SF%INIT_IDS(III)=='null') CYCLE
                  INIT_INDEX = SF%INIT_INDICES(III)
                  CALL INSERT_VOLUMETRIC_PARTICLES
               ENDDO
            ENDDO
         ENDIF
      ENDDO

      DO INIT_INDEX=1,N_INIT
         IF (INITIALIZATION(INIT_INDEX)%INVOKED_BY_SURF) INITIALIZATION(INIT_INDEX)%ALREADY_INSERTED = .TRUE.
      ENDDO

   ENDIF

   ! Insert particles in ducts

   IF (DUCT_HT) CALL INSERT_DUCT_PARTICLES

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
      IF (ILPC < 1 .AND. SF%N_LPC == 0) CYCLE
      IF (T >= SF%PARTICLE_INSERT_CLOCK(NM)) SF%PARTICLE_INSERT_CLOCK(NM) = SF%PARTICLE_INSERT_CLOCK(NM) + SF%DT_INSERT
      IF (T >= SF%PARTICLE_INSERT_CLOCK(NM)) INSERT_ANOTHER_BATCH = .TRUE.
   ENDDO

   IF (.NOT.INSERT_ANOTHER_BATCH) EXIT OVERALL_INSERT_LOOP

ENDDO OVERALL_INSERT_LOOP

! Compute initial particle CFL

IF (PARTICLE_CFL) THEN
   DO IP=1,NLP
      LP => MESHES(NM)%LAGRANGIAN_PARTICLE(IP)
      LPC => LAGRANGIAN_PARTICLE_CLASS(LP%CLASS_INDEX)
      IF (LPC%MASSLESS_TRACER .OR. LPC%MASSLESS_TARGET) CYCLE
      BC => MESHES(NM)%BOUNDARY_COORD(LP%BC_INDEX)
      CALL GET_IJK(BC%X,BC%Y,BC%Z,NM,XI,YJ,ZK,II,JJ,KK)
      PART_UVWMAX = MAX(PART_UVWMAX,MAX(ABS(LP%U)*RDX(II),ABS(LP%V)*RDY(JJ),ABS(LP%W)*RDZ(KK)))
   ENDDO
ENDIF

! If any of the newly inserted particles finds itself in a neighboring mesh, set a flag to indicate that an MPI exchange must
! be done to transfer that particle to the MPI process that controls the neighboring mesh.

DO N=1,N_NEIGHBORING_MESHES
   IF (OMESH(NEIGHBORING_MESH(N))%PARTICLE_SEND_BUFFER%N_ITEMS>0) EXCHANGE_INSERTED_PARTICLES = .TRUE.
ENDDO

T_USED(8)=T_USED(8)+CURRENT_TIME()-TNOW

CONTAINS


!> \brief Insert droplets for sprinklers and nozzles

SUBROUTINE INSERT_SPRAY_PARTICLES

USE MEMORY_FUNCTIONS, ONLY: ALLOCATE_STORAGE
INTEGER :: I,OI

! Loop over all devices, but look for sprinklers or nozzles. Count actuated sprinklers for output purposes.

N_ACTUATED_SPRINKLERS = 0

SPRINKLER_INSERT_LOOP: DO KS=1,N_DEVC

   DV => DEVICE(KS)
   PY => PROPERTY(DV%PROP_INDEX)
   IF (PY%PART_ID == 'null')  CYCLE SPRINKLER_INSERT_LOOP
   IF (.NOT.DV%CURRENT_STATE) CYCLE SPRINKLER_INSERT_LOOP
   IF (PY%QUANTITY=='SPRINKLER LINK TEMPERATURE') N_ACTUATED_SPRINKLERS = N_ACTUATED_SPRINKLERS + 1

   IF (DV%MESH/=NM)           CYCLE SPRINKLER_INSERT_LOOP

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
         PIPE_PRESSURE = EVALUATE_RAMP(REAL(DEVC_PIPE_OPERATING(DV%PIPE_INDEX),EB),PY%PRESSURE_RAMP_INDEX)
         D_PRES_FACTOR = (PY%OPERATING_PRESSURE/PIPE_PRESSURE)**ONTH
         FLOW_RATE = PY%K_FACTOR*SQRT(PIPE_PRESSURE)
      ELSE
         PIPE_PRESSURE = PY%OPERATING_PRESSURE
         D_PRES_FACTOR = 1.0_EB
         FLOW_RATE = PY%FLOW_RATE
      ENDIF
      FLOW_RATE = FLOW_RATE*LPC%DENSITY/1000._EB/60._EB ! convert from L/min to kg/s
   ENDIF

   FLOW_RATE = EVALUATE_RAMP(TSI,PY%FLOW_RAMP_INDEX,TAU=PY%FLOW_TAU)*FLOW_RATE ! kg/s

   IF (FLOW_RATE <= TWO_EPSILON_EB) THEN
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

      CALL ALLOCATE_STORAGE(NM,LP_INDEX=NLP,LPC_INDEX=PY%PART_INDEX,SURF_INDEX=LPC%SURF_INDEX)

      LAGRANGIAN_PARTICLE => MESHES(NM)%LAGRANGIAN_PARTICLE
      LP=>LAGRANGIAN_PARTICLE(NLP)
      LP%CLASS_INDEX = PY%PART_INDEX

      BC=>MESHES(NM)%BOUNDARY_COORD(LP%BC_INDEX)

      ! Set PARTICLE properties

      PARTICLE_TAG = PARTICLE_TAG + NMESHES
      LP%TAG = PARTICLE_TAG
      LP%T_INSERT = T

      ! Randomly choose particle direction angles, theta and phi

      CHOOSE_COORDS: DO
         PICK_PATTERN: IF(PY%SPRAY_PATTERN_INDEX>0) THEN ! Use spray pattern table
            TA => TABLES(PY%SPRAY_PATTERN_INDEX)
            CALL RANDOM_NUMBER(RN)
            FIND_ROW: DO II=1,TA%NUMBER_ROWS
               IF (REAL(RN,EB)>PY%TABLE_ROW(II)) CYCLE FIND_ROW
               EXIT FIND_ROW
            ENDDO FIND_ROW
            CALL RANDOM_NUMBER(RN)
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
         ELSE PICK_PATTERN ! Use conical spray
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
         BC%X = DV%X + PY%OFFSET*XTMP
         BC%Y = DV%Y + PY%OFFSET*YTMP
         BC%Z = DV%Z + PY%OFFSET*ZTMP
         IF (TWO_D) THEN
            LP%V = 0._EB
            BC%Y = DV%Y
         ENDIF
         BC%IOR = 0

         ! If the particle position is outside the current mesh, exit the loop and the particle will be sent to another mesh
         ! or eliminated by the call to REMOVE_PARTICLES at the end of the subroutine.

         IF (BC%X<=XS .OR. BC%X>=XF .OR. BC%Y<=YS .OR. BC%Y>=YF .OR. BC%Z<=ZS .OR. BC%Z>=ZF) THEN
            EXIT CHOOSE_COORDS
         ELSE
            CALL GET_IJK(BC%X,BC%Y,BC%Z,NM,XI,YJ,ZK,II,JJ,KK)
            IC = CELL_INDEX(II,JJ,KK)
            BC%IIG = II; BC%JJG = JJ; BC%KKG = KK
            BC%II  = II; BC%JJ  = JJ; BC%KK  = KK
            IF (.NOT.CELL(IC)%SOLID) EXIT CHOOSE_COORDS
         ENDIF

      ENDDO CHOOSE_COORDS

      ! Randomly choose PARTICLE size according to Cumulative Distribution Function (CDF)

      CALL INITIALIZE_SINGLE_PARTICLE

      ! Adjust particle size to account for pressure dependence of nozzle

      IF (LPC%LIQUID_DROPLET) THEN
         LP%RADIUS = LP%RADIUS*D_PRES_FACTOR
         LP%MASS   = LP%MASS*D_PRES_FACTOR**3
      ENDIF

      ! Sum up mass of liquid being introduced

      MASS_SUM = MASS_SUM + LP%PWT*LP%MASS
      DROP_SUM = DROP_SUM + 1

   ENDDO PARTICLE_INSERT_LOOP

   ! Compute weighting factor for the PARTICLEs just inserted

   IF (DROP_SUM > 0) THEN
      IF (LPC%LIQUID_DROPLET) THEN
         PWT0 = LPC%N_STRATA*FLOW_RATE/(LPC%DENSITY*LPC%MEAN_DROPLET_VOLUME*REAL(PY%PARTICLES_PER_SECOND,EB))/D_PRES_FACTOR**3
      ELSE
         PWT0 = FLOW_RATE*REAL(DROP_SUM,EB)/(MASS_SUM*REAL(PY%PARTICLES_PER_SECOND,EB))
      ENDIF
      DO I=1,N_INSERT
         N = LP_INDEX_LOOKUP(I)
         LAGRANGIAN_PARTICLE(N)%PWT = LAGRANGIAN_PARTICLE(N)%PWT * PWT0
      ENDDO
   ENDIF

   DEALLOCATE(LP_INDEX_LOOKUP)

   ! Indicate that PARTICLEs from this device have been inserted at this time T

   DV%T = T

ENDDO SPRINKLER_INSERT_LOOP

! There might be particles that were inserted outside the current mesh. Either remove them or transfer them to another mesh.

CALL REMOVE_PARTICLES(T,NM)

END SUBROUTINE INSERT_SPRAY_PARTICLES


!> \brief Loop through all boundary cells and insert particles if appropriate

SUBROUTINE INSERT_VENT_PARTICLES

DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
   CALL PARTICLE_FACE_INSERT(WALL_INDEX=IW)
ENDDO

DO ICF=INTERNAL_CFACE_CELLS_LB+1,INTERNAL_CFACE_CELLS_LB+N_INTERNAL_CFACE_CELLS
   CALL PARTICLE_FACE_INSERT(CFACE_INDEX=ICF)
ENDDO

END SUBROUTINE INSERT_VENT_PARTICLES


SUBROUTINE PARTICLE_FACE_INSERT(WALL_INDEX,CFACE_INDEX)

USE COMPLEX_GEOMETRY, ONLY : RANDOM_CFACE_XYZ
USE PHYSICAL_FUNCTIONS, ONLY: GET_PARTICLE_ENTHALPY
USE MEMORY_FUNCTIONS, ONLY: ALLOCATE_STORAGE
INTEGER, INTENT(IN), OPTIONAL :: WALL_INDEX,CFACE_INDEX
INTEGER :: I,N_LPC,N_LPC_MAX,ITER,LPC_INDEX,INSERT_TYPE,NPPC
REAL(EB):: CFA_X,CFA_Y,CFA_Z,RN,VEL_PART,C_S,H_1, H_2,TMP_PART,TMP_GUESS,EMBER_COUNT
LOGICAL :: RETURN_FLAG
TYPE (BOUNDARY_COORD_TYPE), POINTER :: CFA_BC
TYPE (CFACE_TYPE), POINTER :: CFA
TYPE (BOUNDARY_PROP2_TYPE), POINTER :: B2

WALL_OR_CFACE_IF_1: IF (PRESENT(WALL_INDEX)) THEN
   WC => MESHES(NM)%WALL(WALL_INDEX)
   IF (WC%OD_INDEX>0) ONE_D => MESHES(NM)%BOUNDARY_ONE_D(WC%OD_INDEX)
   B1 => MESHES(NM)%BOUNDARY_PROP1(WC%B1_INDEX)
   B2 => MESHES(NM)%BOUNDARY_PROP2(WC%B2_INDEX)
   BC => MESHES(NM)%BOUNDARY_COORD(WC%BC_INDEX)
   SF => SURFACE(WC%SURF_INDEX)
   IF (WC%BOUNDARY_TYPE/=SOLID_BOUNDARY) RETURN
   II = BC%II
   JJ = BC%JJ
   KK = BC%KK
   IC = CELL_INDEX(II,JJ,KK)
   IF (.NOT.CELL(IC)%SOLID) RETURN
   IOR = BC%IOR
ELSEIF (PRESENT(CFACE_INDEX)) THEN
   CFA => MESHES(NM)%CFACE(CFACE_INDEX)
   IF (CFA%OD_INDEX>0) ONE_D => MESHES(NM)%BOUNDARY_ONE_D(CFA%OD_INDEX)
   B1 => MESHES(NM)%BOUNDARY_PROP1(CFA%B1_INDEX)
   B2 => MESHES(NM)%BOUNDARY_PROP2(CFA%B2_INDEX)
   BC => MESHES(NM)%BOUNDARY_COORD(CFA%BC_INDEX)
   SF => SURFACE(CFA%SURF_INDEX)
   IF (CFA%BOUNDARY_TYPE/=SOLID_BOUNDARY) RETURN
ENDIF WALL_OR_CFACE_IF_1

ILPC = SF%PART_INDEX; IF (ILPC < 1 .AND. SF%N_LPC ==0) RETURN

IIG = BC%IIG
JJG = BC%JJG
KKG = BC%KKG

RETURN_FLAG = .FALSE.

IF (T < SF%PARTICLE_INSERT_CLOCK(NM)) RETURN

NPPC=SF%NPPC

INSERT_TYPE_LOOP: DO INSERT_TYPE = 1,2

   SELECT CASE (INSERT_TYPE)
      CASE (1) ! PART_ID on SURF
         IF (ILPC <=0) CYCLE INSERT_TYPE_LOOP
         IF (T < B1%T_IGN) CYCLE INSERT_TYPE_LOOP
         IF (SF%PARTICLE_SURFACE_DENSITY>0._EB .AND. T>T_BEGIN) CYCLE INSERT_TYPE_LOOP

         ! Ember generation
         IF (ANY(SF%EMBER_GENERATION_HEIGHT>=0._EB) .OR. SF%EMBER_YIELD>0._EB) THEN
            ! specify generation only for regions of burning or accumulated embers
            IF (B1%M_DOT_G_PP_ADJUST(REACTION(1)%FUEL_SMIX_INDEX)<TWO_EPSILON_EB .AND. &
               B2%LP_EMPUA(ILPC)<TWO_EPSILON_EB) CYCLE INSERT_TYPE_LOOP
            ! If using an EMBER_YIELD, check enough mass has accumulated or fire has gone out
            IF (SF%EMBER_YIELD>0._EB) THEN
               LPC => LAGRANGIAN_PARTICLE_CLASS(ILPC)
               EMBER_COUNT = B2%LP_EMPUA(ILPC)*B1%AREA_ADJUST*B1%AREA/LPC%INITIAL_MASS
               IF (B2%LP_EMPUA(ILPC)<TWO_EPSILON_EB .OR. (EMBER_COUNT<SF%EMBER_TRACKING_RATIO .AND. &
                   B1%M_DOT_G_PP_ADJUST(REACTION(1)%FUEL_SMIX_INDEX)>TWO_EPSILON_EB)) CYCLE INSERT_TYPE_LOOP
               ! Add multiple particles if neecessary to roughly preserve EMBER_TRACKING_RATIO
               NPPC = MAX(1,NINT(EMBER_COUNT/SF%EMBER_TRACKING_RATIO))
               FLOW_RATE=B2%LP_EMPUA(ILPC)/SF%DT_INSERT
               B2%LP_EMPUA(ILPC) = 0._EB
            ENDIF
         ENDIF

         ! If there is a flow ramp and the value is zero, return.
         IF (SF%PARTICLE_MASS_FLUX > 0._EB) THEN
            IF (ABS(B1%T_IGN-T_BEGIN)<=TWO_EPSILON_EB .AND. SF%RAMP(TIME_PART)%INDEX>=1) THEN
               TSI = T
            ELSE
               TSI = T - B1%T_IGN
            ENDIF
            FLOW_RATE = EVALUATE_RAMP(TSI,SF%RAMP(TIME_PART)%INDEX,TAU=SF%RAMP(TIME_PART)%TAU)*SF%PARTICLE_MASS_FLUX
            IF (FLOW_RATE < TWO_EPSILON_EB) RETURN
         ENDIF
         LPC => LAGRANGIAN_PARTICLE_CLASS(ILPC)

         ! Evalutate if we need to skip ILPC and return if no N_LPC
         IF (LPC%DEVC_INDEX>0) THEN
            IF (.NOT.DEVICE(LPC%DEVC_INDEX)%CURRENT_STATE) RETURN_FLAG = .TRUE.
         ENDIF
         IF (LPC%CTRL_INDEX>0) THEN
            IF (.NOT.CONTROL(LPC%CTRL_INDEX)%CURRENT_STATE) RETURN_FLAG = .TRUE.
         ENDIF
         IF (NM > 1) THEN
            IF (INTERPOLATED_MESH(IIG,JJG,KKG) > 0) RETURN_FLAG = .TRUE.
         ENDIF
         IF (RETURN_FLAG) CYCLE INSERT_TYPE_LOOP
         N_LPC_MAX = 1

      CASE (2) ! PART_ID on MATL
         IF (SF%N_LPC==0) RETURN
         N_LPC_MAX = SF%N_LPC

   END SELECT

   ALLOCATE(LP_INDEX_LOOKUP(NPPC))

   LPC_LOOP: DO N_LPC=1,N_LPC_MAX
      LP_INDEX_LOOKUP = 0
      IF (INSERT_TYPE==1) THEN
         LPC_INDEX = ILPC
      ELSE
         LPC_INDEX = SF%MATL_PART_INDEX(N_LPC)
      ENDIF
      ! Get particle temperature for MATL particle
      IF (INSERT_TYPE==2) THEN

         IF (ONE_D%PART_MASS(N_LPC) < TWO_EPSILON_EB) CYCLE
         LPC => LAGRANGIAN_PARTICLE_CLASS(SF%MATL_PART_INDEX(N_LPC))

         TMP_GUESS = ONE_D%TMP(1)
         H_1 = ONE_D%PART_ENTHALPY(N_LPC)/ONE_D%PART_MASS(N_LPC)
         C_S = H_1/TMP_GUESS
         ITER = 0
         T_SEARCH: DO
            ITER = ITER + 1
            H_2 = GET_PARTICLE_ENTHALPY(SF%MATL_PART_INDEX(N_LPC),TMP_GUESS)
            C_S = H_2/TMP_GUESS
            TMP_PART = TMP_GUESS - (H_2 - H_1)/C_S
            IF (ABS(TMP_PART - TMP_GUESS) < 0.0001_EB) EXIT T_SEARCH
            IF (ITER > 20) THEN
               TMP_PART = 0.5_EB*(TMP_PART+TMP_GUESS)
               EXIT T_SEARCH
            ENDIF
            TMP_GUESS = TMP_PART
         ENDDO T_SEARCH

      ENDIF

      ! Loop over all particles for the WALL_INDEX-th cell

      MASS_SUM = 0._EB

      PARTICLE_INSERT_LOOP: DO I=1,NPPC

        ! Insert a single droplet at wall cell WALL_INDEX or cut-cell face CFACE_INDEXF

         IF (NLP+1>MAXIMUM_PARTICLES) THEN
            CALL REMOVE_OLDEST_PARTICLE(NM,LPC_INDEX,NLP,NEW_LP_INDEX)
            IF (I>1) LP_INDEX_LOOKUP(I-1)=NEW_LP_INDEX
         ELSE
            NLP = NLP+1
         ENDIF
         LP_INDEX_LOOKUP(I) = NLP

         CALL ALLOCATE_STORAGE(NM,LP_INDEX=NLP,LPC_INDEX=LPC_INDEX,SURF_INDEX=LPC%SURF_INDEX)

         LAGRANGIAN_PARTICLE => MESHES(NM)%LAGRANGIAN_PARTICLE
         LP => LAGRANGIAN_PARTICLE(NLP)
         LP%CLASS_INDEX = LPC_INDEX
         PARTICLE_TAG = PARTICLE_TAG + NMESHES
         LP%TAG = PARTICLE_TAG

         ! Reassign pointers after calling ALLOCATE

         BC => MESHES(NM)%BOUNDARY_COORD(LP%BC_INDEX)

         BC%IOR = 0  ! Particle is not stuck to a wall

         IF (PRESENT(WALL_INDEX)) THEN
            WC    => MESHES(NM)%WALL(WALL_INDEX)
            B1    => MESHES(NM)%BOUNDARY_PROP1(WC%B1_INDEX)
            IF (WC%OD_INDEX>0) ONE_D => MESHES(NM)%BOUNDARY_ONE_D(WC%OD_INDEX)
         ELSEIF (PRESENT(CFACE_INDEX)) THEN
            CFA   => MESHES(NM)%CFACE(CFACE_INDEX)
            B1    => MESHES(NM)%BOUNDARY_PROP1(CFA%B1_INDEX)
            CFA_BC => MESHES(NM)%BOUNDARY_COORD(CFA%BC_INDEX)
            IF (CFA%OD_INDEX>0) ONE_D => MESHES(NM)%BOUNDARY_ONE_D(CFA%OD_INDEX)
         ENDIF

         ! Ember flag to be used for outputs

         IF (INSERT_TYPE==1 .AND. (ANY(SF%EMBER_GENERATION_HEIGHT>=0._EB) .OR. SF%EMBER_YIELD>0._EB)) LP%EMBER=.TRUE.

         ! Assign particle position on the cell face

         CALL RANDOM_NUMBER(RN)
         CALL RANDOM_NUMBER(RN2)

         WALL_OR_CFACE_IF_2: IF (PRESENT(WALL_INDEX)) THEN
            SELECT CASE (ABS(IOR))
               CASE(1)
                  IF (IOR== 1) BC%X = X(II)   + VENT_OFFSET*DX(II+1)
                  IF (IOR==-1) BC%X = X(II-1) - VENT_OFFSET*DX(II-1)
                  BC%Y = Y(JJ-1) + DY(JJ)*REAL(RN,EB)
                  BC%Z = Z(KK-1) + DZ(KK)*REAL(RN2,EB)
               CASE(2)
                  IF (IOR== 2) BC%Y = Y(JJ)   + VENT_OFFSET*DY(JJ+1)
                  IF (IOR==-2) BC%Y = Y(JJ-1) - VENT_OFFSET*DY(JJ-1)
                  BC%X = X(II-1) + DX(II)*REAL(RN,EB)
                  BC%Z = Z(KK-1) + DZ(KK)*REAL(RN2,EB)
               CASE(3)
                  IF (IOR== 3) THEN
                     BC%Z = Z(KK)   + VENT_OFFSET*DZ(KK+1)
                     IF (ANY(SF%EMBER_GENERATION_HEIGHT>=0._EB)) THEN
                        CALL RANDOM_NUMBER(RN3)
                        BC%Z = Z(KK) + SF%EMBER_GENERATION_HEIGHT(1) + &
                           (SF%EMBER_GENERATION_HEIGHT(2)-SF%EMBER_GENERATION_HEIGHT(1))*REAL(RN3,EB)
                     ENDIF
                  ENDIF
                  IF (IOR==-3) BC%Z = Z(KK-1) - VENT_OFFSET*DZ(KK-1)
                  BC%X = X(II-1) + DX(II)*REAL(RN,EB)
                  BC%Y = Y(JJ-1) + DY(JJ)*REAL(RN2,EB)
            END SELECT
            ! Give particles an initial velocity
            IF (.NOT.LPC%STATIC) THEN
               IF (SF%VEL_PART >-999999._EB) THEN
                  VEL_PART = SF%VEL_PART
               ELSE
                  VEL_PART = B1%U_NORMAL
               ENDIF
               SELECT CASE(IOR)
                  CASE( 1)
                     LP%U = -VEL_PART
                     LP%V = SF%VEL_T(1)
                     LP%W = SF%VEL_T(2)
                  CASE(-1)
                     LP%U =  VEL_PART
                     LP%V = SF%VEL_T(1)
                     LP%W = SF%VEL_T(2)
                  CASE( 2)
                     LP%U = SF%VEL_T(1)
                     LP%V = -VEL_PART
                     LP%W = SF%VEL_T(2)
                  CASE(-2)
                     LP%U = SF%VEL_T(1)
                     LP%V =  VEL_PART
                     LP%W = SF%VEL_T(2)
                  CASE( 3)
                     LP%U = SF%VEL_T(1)
                     LP%V = SF%VEL_T(2)
                     LP%W = -VEL_PART
                  CASE(-3)
                     LP%U = SF%VEL_T(1)
                     LP%V = SF%VEL_T(2)
                     LP%W =  VEL_PART
               END SELECT
            ENDIF
         ELSEIF (PRESENT(CFACE_INDEX)) THEN
            CALL RANDOM_CFACE_XYZ(NM,CFA,CFA_X,CFA_Y,CFA_Z)
            BC%X = CFA_X + CFA_BC%NVEC(1)*VENT_OFFSET*DX(IIG)
            BC%Y = CFA_Y + CFA_BC%NVEC(2)*VENT_OFFSET*DY(JJG)
            BC%Z = CFA_Z + CFA_BC%NVEC(3)*VENT_OFFSET*DZ(KKG)
            IF (ANY(SF%EMBER_GENERATION_HEIGHT>=0._EB)) THEN
               CALL RANDOM_NUMBER(RN3)
               BC%Z = CFA_Z + SF%EMBER_GENERATION_HEIGHT(1) + &
                      (SF%EMBER_GENERATION_HEIGHT(2)-SF%EMBER_GENERATION_HEIGHT(1))*REAL(RN3,EB)
            ENDIF
            VEL_PART = B1%U_NORMAL
            IF (.NOT.LPC%STATIC .AND. SF%VEL_PART >-999999._EB) VEL_PART = SF%VEL_PART
            LP%U = DOT_PRODUCT(CFA_BC%NVEC,(/-VEL_PART,SF%VEL_T(1),SF%VEL_T(2)/))
            LP%V = DOT_PRODUCT(CFA_BC%NVEC,(/SF%VEL_T(1),-VEL_PART,SF%VEL_T(2)/))
            LP%W = DOT_PRODUCT(CFA_BC%NVEC,(/SF%VEL_T(1),SF%VEL_T(2),-VEL_PART/))

         ENDIF WALL_OR_CFACE_IF_2

         ! Update idicies in case offset puts location in a different cell
         CALL GET_IJK(BC%X,BC%Y,BC%Z,NM,XI,YJ,ZK,BC%IIG,BC%JJG,BC%KKG)
         BC%II = BC%IIG; BC%JJ = BC%JJG; BC%KK = BC%KKG

         ! Save the insertion time (TP) and scalar property (SP) for the particle

         LP%T_INSERT = T

         CALL INITIALIZE_SINGLE_PARTICLE

         IF (INSERT_TYPE==2 .AND. LPC%SOLID_PARTICLE) MESHES(NM)%BOUNDARY_ONE_D(LP%OD_INDEX)%TMP = TMP_PART

         IF (.NOT.LPC%MASSLESS_TRACER .AND. .NOT.LPC%MASSLESS_TARGET) MASS_SUM = MASS_SUM + LP%PWT*LP%MASS

      ENDDO PARTICLE_INSERT_LOOP

      ! Adjust the particle weighting factors to get the right mass flux

      MASS_CHECK: IF (MASS_SUM > 0._EB) THEN
         SELECT CASE (INSERT_TYPE)
            CASE (1) ! PART_ID on SURF
               IF (MASS_SUM > 0._EB) THEN
                  IF (SF%PARTICLE_MASS_FLUX > 0._EB .OR. SF%EMBER_YIELD>0._EB) THEN
                     DO I=1,NPPC
                        LP => LAGRANGIAN_PARTICLE(LP_INDEX_LOOKUP(I))
                        LP%PWT = LP%PWT * FLOW_RATE*B1%AREA_ADJUST*B1%AREA*SF%DT_INSERT/MASS_SUM
                     ENDDO
                  ELSEIF (SF%PARTICLE_SURFACE_DENSITY > 0._EB) THEN
                     DO I=1,NPPC
                        LP => LAGRANGIAN_PARTICLE(LP_INDEX_LOOKUP(I))
                        LP%PWT = LP%PWT * SF%PARTICLE_SURFACE_DENSITY*B1%AREA_ADJUST*B1%AREA/MASS_SUM
                     ENDDO
                  ENDIF
               ENDIF

            CASE (2) ! PART_ID on MATL
               IF (MASS_SUM > 0._EB) THEN
                  DO I=1,NPPC
                     LP => LAGRANGIAN_PARTICLE(LP_INDEX_LOOKUP(I))
                     LP%PWT = LP%PWT * ONE_D%PART_MASS(N_LPC)*B1%AREA_ADJUST*B1%AREA* &
                              MIN(1._EB,SF%DT_INSERT/B1%T_MATL_PART)/MASS_SUM
                  ENDDO
               ENDIF

         END SELECT

      ENDIF MASS_CHECK

   ENDDO LPC_LOOP

   IF(ALLOCATED(LP_INDEX_LOOKUP)) DEALLOCATE(LP_INDEX_LOOKUP)

   IF (INSERT_TYPE==2) THEN
      ! Decrement mass, enthalpy, and aggregation time for current insertion interval
      ONE_D%PART_MASS = MAX(0._EB,1-SF%DT_INSERT/(B1%T_MATL_PART+TINY_EB))*ONE_D%PART_MASS
      ONE_D%PART_ENTHALPY = MAX(0._EB,1-SF%DT_INSERT/(B1%T_MATL_PART+TINY_EB))*ONE_D%PART_ENTHALPY
      B1%T_MATL_PART = MAX(0._EB,B1%T_MATL_PART-SF%DT_INSERT)
   ENDIF

ENDDO INSERT_TYPE_LOOP

END SUBROUTINE PARTICLE_FACE_INSERT


!> \brief Loop over all INIT lines and look for particles inserted within a specified volume

SUBROUTINE INSERT_VOLUMETRIC_PARTICLES

USE MEMORY_FUNCTIONS, ONLY: ALLOCATE_STORAGE
INTEGER :: IIP,N_INSERT,I1,J1,K1,I2,J2,K2,N,N_PARTICLES_INSERT,ND,ICC,IFACE,INDCF,I_RAND,N_INSERT_CELLS
REAL(EB) :: XC1,XC2,YC1,YC2,ZC1,ZC2,X0,Y0,Z0,RR,RRI,HH,INSERT_VOLUME,INPUT_VOLUME,VOLUME_SPLIT_FACTOR,LP_X,LP_Y,LP_Z,RAMP_FACTOR,&
            IN_X1,IN_X2,IN_Y1,IN_Y2,IN_Z1,IN_Z2,IN_X0,IN_Y0,IN_Z0,VCX,VCY,VCZ,MOIST_CONT,FILLED_VOLUME,DIST,DIST_MIN,&
            P_VECTOR(3),P_VECTOR_MIN(3),NVEC_MIN(3)
LOGICAL :: CC_VALID
TYPE (CC_CUTFACE_TYPE), POINTER :: CF

IN => INITIALIZATION(INIT_INDEX)

! Determine if the INITIALIZATION type involves particles. If not, return.

ILPC = IN%PART_INDEX
IF (ILPC<1) RETURN
IF (IN%SINGLE_INSERTION .AND. IN%ALREADY_INSERTED(NM)) RETURN

! If there is a RAMP for MASS_PER_TIME or MASS_PER_VOLUME, evaluate it now and if zero, return.

IF (IN%RAMP_PART_INDEX>0) THEN
   RAMP_FACTOR = EVALUATE_RAMP(T,IN%RAMP_PART_INDEX)
   IF (RAMP_FACTOR<TWO_EPSILON_EB) RETURN
ELSE
   RAMP_FACTOR = 1._EB
ENDIF

! Determine if the particles/PARTICLEs are controlled by devices

LPC => LAGRANGIAN_PARTICLE_CLASS(ILPC)

IF (IN%DEVC_INDEX>0) THEN
   IF (.NOT.DEVICE(IN%DEVC_INDEX)%CURRENT_STATE) THEN
      IN%PARTICLE_INSERT_CLOCK(NM) = T
      RETURN
   ENDIF
ENDIF
IF (IN%CTRL_INDEX>0) THEN
   IF (.NOT.CONTROL(IN%CTRL_INDEX)%CURRENT_STATE) THEN
      IN%PARTICLE_INSERT_CLOCK(NM) = T
      RETURN
   ENDIF
ENDIF

! If it is not time to insert particles for this INITIALIZATION block, cycle.

IF (T < IN%PARTICLE_INSERT_CLOCK(NM)) RETURN

! Start processing the INITIALIZATION info

IF (IN%N_PARTICLES==0 .AND. IN%N_PARTICLES_PER_CELL==0) RETURN

! Adjust INIT values for voxelized bulk density data
IF (IN%BULK_DENSITY_FILE/='null') THEN
   ! Center of voxel
   READ(LU_VEG_IN) VCX,VCY,VCZ
   IN%X1=VCX-VDX/2._EB
   IN%X2=VCX+VDX/2._EB
   IN%Y1=VCY-VDY/2._EB
   IN%Y2=VCY+VDY/2._EB
   IN%Z1=VCZ-VDZ/2._EB
   IN%Z2=VCZ+VDZ/2._EB
   ! Vegetation mass in voxel
   READ(LU_VEG_IN) IN%MASS_PER_VOLUME
   ! Adjust with multiplier
   IN%MASS_PER_VOLUME = MAX(0._EB,IN%BULK_DENSITY_FACTOR)*IN%MASS_PER_VOLUME
   MOIST_CONT=SURFACE(LAGRANGIAN_PARTICLE_CLASS(IN%PART_INDEX)%SURF_INDEX)%MOISTURE_CONTENT(1)
   IF (MOIST_CONT>=0._EB) IN%MASS_PER_VOLUME = IN%MASS_PER_VOLUME*(1._EB+MOIST_CONT)
ENDIF

! Apply coordinate offset if needed

IN_X1 = X_OFFSET + IN%X1
IN_X2 = X_OFFSET + IN%X2
IN_Y1 = Y_OFFSET + IN%Y1
IN_Y2 = Y_OFFSET + IN%Y2
IN_Z1 = Z_OFFSET + IN%Z1
IN_Z2 = Z_OFFSET + IN%Z2
IN_X0 = X_OFFSET + IN%X0
IN_Y0 = X_OFFSET + IN%Y0
IN_Z0 = X_OFFSET + IN%Z0

! If the INIT volume is outside the current mesh, return

IF (IN%SHAPE/='RING' .AND. IN%SHAPE/='LINE') THEN
   IF ((IN_X1-XF)>-50._EB*TWO_EPSILON_EB .OR. (IN_X2-XS)<50._EB*TWO_EPSILON_EB .OR. &
       (IN_Y1-YF)>-50._EB*TWO_EPSILON_EB .OR. (IN_Y2-YS)<50._EB*TWO_EPSILON_EB .OR. &
       (IN_Z1-ZF)>-50._EB*TWO_EPSILON_EB .OR. (IN_Z2-ZS)<50._EB*TWO_EPSILON_EB) RETURN
ELSEIF (IN%SHAPE=='RING') THEN
   IF (RING_MESH_INTERSECTION_ARC(NM,IN%X0,IN%Y0,IN%RADIUS)<TWO_EPSILON_EB) RETURN
ELSE
   IF (IN_X1>XF .OR. IN_X2<XS .OR. IN_Y1>YF .OR. IN_Y2<YS .OR. IN_Z1>ZF .OR. IN_Z2<ZS) RETURN
ENDIF

! Skip mesh that is contained completely within a ring

IF (IN%SHAPE=='RING' .AND. IN_X1<XS .AND. IN_X2>XF .AND. IN_Y1<YS .AND. IN_Y2>YF .AND. IN_Z1<ZS .AND. IN_Z2>ZF) RETURN

! Cut off parts of the INIT region that are outside the current mesh

X1 = MAX(IN_X1,XS)
X2 = MIN(IN_X2,XF)
Y1 = MAX(IN_Y1,YS)
Y2 = MIN(IN_Y2,YF)
Z1 = MAX(IN_Z1,ZS)
Z2 = MIN(IN_Z2,ZF)

! Compute the volume of the INIT region

GEOM_SELECT: SELECT CASE(IN%SHAPE)
   CASE('BLOCK')
      IF (ABS(X2-X1)<TWO_EPSILON_EB) THEN
         IF (ABS(Y2-Y1)<TWO_EPSILON_EB) THEN
            INSERT_VOLUME = Z2 - Z1
            INPUT_VOLUME = IN_Z2 - IN_Z1
         ELSE
            IF(ABS(Z2-Z1)<TWO_EPSILON_EB) THEN
               INSERT_VOLUME = Y2 - Y1
               INPUT_VOLUME = IN_Y2 - IN_Y1
            ELSE
               INSERT_VOLUME = (Y2-Y1)*(Z2-Z1)
               INPUT_VOLUME = (IN_Y2-IN_Y1)*(IN_Z2-IN_Z1)
            ENDIF
         ENDIF
         EXIT GEOM_SELECT
      ENDIF
      IF (ABS(Y2-Y1)<TWO_EPSILON_EB) THEN
         IF (ABS(X2-X1)<TWO_EPSILON_EB) THEN
            INSERT_VOLUME = Z2 - Z1
            INPUT_VOLUME = IN_Z2 - IN_Z1
         ELSE
            IF(ABS(Z2-Z1)<TWO_EPSILON_EB) THEN
               INSERT_VOLUME = X2 - X1
               INPUT_VOLUME = IN_X2 - IN_X1
            ELSE
               INSERT_VOLUME = (X2-X1)*(Z2-Z1)
               INPUT_VOLUME = (IN_X2-IN_X1)*(IN_Z2-IN_Z1)
            ENDIF
         ENDIF
         EXIT GEOM_SELECT
      ENDIF
      IF (ABS(Z2-Z1)<TWO_EPSILON_EB) THEN
         IF (ABS(Y2-Y1)<TWO_EPSILON_EB) THEN
            INSERT_VOLUME = X2 - X1
            INPUT_VOLUME = IN_X2 - IN_X1
         ELSE
            IF(ABS(X2-X1)<TWO_EPSILON_EB) THEN
               INSERT_VOLUME = Y2 - Y1
               INPUT_VOLUME = IN_Y2 - IN_Y1
            ELSE
               INSERT_VOLUME = (X2-X1)*(Y2-Y1)
               INPUT_VOLUME = (IN_X2-IN_X1)*(IN_Y2-IN_Y1)
            ENDIF
         ENDIF
         EXIT GEOM_SELECT
      ENDIF
      IF (CYLINDRICAL) THEN
         INSERT_VOLUME = 0.5_EB*(X2**2-X1**2)*(Y2-Y1)*(Z2-Z1)
         INPUT_VOLUME  = 0.5_EB*(IN_X2**2-IN_X1**2)*(IN_Y2-IN_Y1)*(IN_Z2-IN_Z1)
      ELSE
         INSERT_VOLUME = (X2-X1)*(Y2-Y1)*(Z2-Z1)
         INPUT_VOLUME  = (IN_X2-IN_X1)*(IN_Y2-IN_Y1)*(IN_Z2-IN_Z1)
      ENDIF
   CASE('CONE','CYLINDER')
      X0 = 0.5_EB*(IN_X1+IN_X2)
      Y0 = 0.5_EB*(IN_Y1+IN_Y2)
      Z0 = IN_Z1
      RR = 0.5_EB*(IN_X2-IN_X1)
      RRI = IN%INNER_RADIUS
      HH = IN_Z2-IN_Z1
      IF (IN%SHAPE=='CONE') THEN
         INSERT_VOLUME = CONE_MESH_INTERSECTION_VOLUME(NM,X0,Y0,Z0,RR,RRI,HH,1)
         INPUT_VOLUME  = PI*(RR**2._EB-RRI**2._EB)*HH/3._EB
      ELSE
         INSERT_VOLUME = CONE_MESH_INTERSECTION_VOLUME(NM,X0,Y0,Z0,RR,RRI,HH,0)
         INPUT_VOLUME  = PI*(RR**2._EB-RRI**2._EB)*HH
      ENDIF
   CASE('RING')
      IF (IN%UNIFORM) THEN
         INSERT_VOLUME = 0._EB
         INPUT_VOLUME  = 0._EB
      ELSE
         ! proportion of the circle arc length within mesh (length not volume in this case)
         X0 = 0.5_EB*(IN_X1+IN_X2)
         Y0 = 0.5_EB*(IN_Y1+IN_Y2)
         RR = 0.5_EB*(IN_X2-IN_X1)
         INSERT_VOLUME = RING_MESH_INTERSECTION_ARC(NM,X0,Y0,RR)
         INPUT_VOLUME  = TWOPI*RR
      ENDIF
   CASE('LINE')
      ! proportion of the circle bounding box within mesh
      INSERT_VOLUME = 0._EB
      INPUT_VOLUME  = 0._EB
END SELECT GEOM_SELECT

IF (INSERT_VOLUME<=0._EB .AND. IN%MASS_PER_VOLUME>0._EB) RETURN

! Assign properties to the particles

MASS_SUM = 0._EB
N_INSERT_CELLS = 0

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

   ! Loop through the particles to be inserted, getting their position and then setting up array space.

   N_INSERT = 0

   INSERT_PARTICLE_LOOP: DO IP=1,N_PARTICLES_INSERT

      ! Get particle coordinates by randomly choosing within the designated volume

      N = 0
      CHOOSE_XYZ_LOOP:  DO
         N = N + 1
         SELECT CASE(IN%SHAPE)
            CASE('BLOCK')
               CALL RANDOM_RECTANGLE(LP_X,LP_Y,LP_Z,X1,X2,Y1,Y2,Z1,Z2)
            CASE('CONE','CYLINDER')
               X0 = 0.5_EB*(IN_X1+IN_X2)
               Y0 = 0.5_EB*(IN_Y1+IN_Y2)
               Z0 = IN_Z1
               RR = 0.5_EB*(IN_X2-IN_X1)
               RRI = IN%INNER_RADIUS
               HH = IN_Z2-IN_Z1
               IF (IN%SHAPE=='CONE')     CALL RANDOM_CONE(NM,LP_X,LP_Y,LP_Z,X0,Y0,Z0,RR,RRI,HH,1)
               IF (IN%SHAPE=='CYLINDER') CALL RANDOM_CONE(NM,LP_X,LP_Y,LP_Z,X0,Y0,Z0,RR,RRI,HH,0)
            CASE('RING')
               X0 = 0.5_EB*(IN_X1+IN_X2)
               Y0 = 0.5_EB*(IN_Y1+IN_Y2)
               RR = 0.5_EB*(IN_X2-IN_X1)
               LP_Z = IN_Z1
               IF (IN%UNIFORM) THEN
                  CALL UNIFORM_RING(LP_X,LP_Y,X0,Y0,RR,IP,N_PARTICLES_INSERT)
               ELSE
                  CALL RANDOM_RING(NM,LP_X,LP_Y,X0,Y0,RR)
               ENDIF
            CASE('LINE')
               LP_X = IN_X1 + (IP-1)*IN%DX
               LP_Y = IN_Y1 + (IP-1)*IN%DY
               LP_Z = IN_Z1 + (IP-1)*IN%DZ
               IF (LPC%ID=='RESERVED TARGET PARTICLE') THEN
                  DO ND=1,N_DEVC
                     DV => DEVICE(ND)
                     IF (IN%ID==DV%INIT_ID .AND. IP==DV%POINT) THEN
                        LP_X = DV%X
                        LP_Y = DV%Y
                        LP_Z = DV%Z
                     ENDIF
                  ENDDO
               ENDIF
         END SELECT

         ! Reject particles that are not in the current mesh.

         IF (LP_X<XS .OR. LP_X>XF .OR. LP_Y<YS .OR. LP_Y>YF .OR. LP_Z<ZS .OR. LP_Z>ZF) CYCLE INSERT_PARTICLE_LOOP

         ! Get mesh indices for particle. If the particle is in a solid cell, get another random point. If the particle
         ! is a member of a line of points and this point is SOLID, just skip it.

         CALL GET_IJK(LP_X,LP_Y,LP_Z,NM,XI,YJ,ZK,II,JJ,KK)

         IF (CELL(CELL_INDEX(II,JJ,KK))%SOLID .AND. IN%SHAPE=='LINE') CYCLE INSERT_PARTICLE_LOOP
         ! Check for solid inside GEOM
         CC_VALID = .TRUE.
         IF (CC_IBM) THEN
            IF (CCVAR(II,JJ,KK,CC_CGSC)==CC_SOLID) THEN
               CC_VALID = .FALSE.
            ELSE
               INDCF = CCVAR(II,JJ,KK,CC_IDCF)
               ! If closest CFACE has positive dot-product of normal and centroid-particle vector, location is assumed valid
               IF (INDCF>0) THEN
                  DIST_MIN=HUGE_EB
                  CF => CUT_FACE(INDCF)
                  CC_VALID = .FALSE.
                  CFA_LOOP1: DO IFACE=1,CF%NFACE
                     P_VECTOR = (/LP_X-CF%XYZCEN(IAXIS,IFACE), &
                                 LP_Y-CF%XYZCEN(JAXIS,IFACE), &
                                 LP_Z-CF%XYZCEN(KAXIS,IFACE)/)
                     DIST = NORM2(P_VECTOR)
                     IF (DIST<DIST_MIN) THEN
                        DIST_MIN = DIST
                        P_VECTOR_MIN = P_VECTOR
                        NVEC_MIN = BOUNDARY_COORD(CFACE(CF%CFACE_INDEX(IFACE))%BC_INDEX)%NVEC
                     ENDIF
                  ENDDO CFA_LOOP1
                  IF (DOT_PRODUCT(NVEC_MIN,P_VECTOR_MIN) > TWO_EPSILON_EB) CC_VALID=.TRUE.
               ENDIF
            ENDIF
         ENDIF
         IF (.NOT.CELL(CELL_INDEX(II,JJ,KK))%SOLID .AND. CC_VALID) EXIT CHOOSE_XYZ_LOOP

         ! If cannot find non-solid grid cell, stop searching

         IF (N>N_PARTICLES_INSERT) EXIT INSERT_PARTICLE_LOOP

      ENDDO CHOOSE_XYZ_LOOP

      N_INSERT = N_INSERT + 1

      ! Allocate space for the particle in the appropriate array

      IF (NLP+1>MAXIMUM_PARTICLES) THEN
         CALL REMOVE_OLDEST_PARTICLE(NM,ILPC,NLP,NEW_LP_INDEX)
         IF (N_INSERT>1) LP_INDEX_LOOKUP(N_INSERT-1) = NEW_LP_INDEX
      ELSE
         NLP = NLP+1
      ENDIF
      LP_INDEX_LOOKUP(N_INSERT) = NLP

      CALL ALLOCATE_STORAGE(NM,LP_INDEX=NLP,LPC_INDEX=ILPC,SURF_INDEX=LPC%SURF_INDEX)

      LAGRANGIAN_PARTICLE => MESHES(NM)%LAGRANGIAN_PARTICLE
      LP => LAGRANGIAN_PARTICLE(NLP)
      LP%CLASS_INDEX = ILPC
      PARTICLE_TAG = PARTICLE_TAG + NMESHES
      LP%TAG = PARTICLE_TAG

      BC => MESHES(NM)%BOUNDARY_COORD(LP%BC_INDEX)

      BC%IOR = 0
      BC%X = LP_X
      BC%Y = LP_Y
      BC%Z = LP_Z
      LP%DX = DX(II)
      LP%DY = DY(JJ)
      LP%DZ = DZ(KK)
      LP%INIT_INDEX = INIT_INDEX

      ! Initialize particle properties

      CALL VOLUME_INIT_PARTICLE

   ENDDO INSERT_PARTICLE_LOOP

ELSEIF (IN%N_PARTICLES_PER_CELL > 0) THEN TOTAL_OR_PER_CELL

   N_INSERT = 0
   INSERT_VOLUME = 0._EB
   FILLED_VOLUME = 0._EB
   CALL GET_IJK(MIN(X1+MICRON,X2),MIN(Y1+MICRON,Y2),MIN(Z1+MICRON,Z2),NM,XI,YJ,ZK,I1,J1,K1)
   CALL GET_IJK(MAX(X2-MICRON,X1),MAX(Y2-MICRON,Y1),MAX(Z2-MICRON,Z1),NM,XI,YJ,ZK,I2,J2,K2)
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
            CC_VALID = .TRUE.
            IF (CC_IBM) THEN
               IF (CCVAR(II,JJ,KK,CC_CGSC)==CC_SOLID) CC_VALID = .FALSE.
            ENDIF
            IF (CELL(CELL_INDEX(II,JJ,KK))%SOLID .OR. .NOT.CC_VALID) CYCLE II_LOOP
            IF (IN%SHAPE=='CONE') THEN
               IF (((XC(II)-X0)**2+(YC(JJ)-Y0)**2<(RRI*(1._EB-(ZC(KK)-Z0)/HH))**2) .OR. &
                  ((XC(II)-X0)**2+(YC(JJ)-Y0)**2>(RR*(1._EB-(ZC(KK)-Z0)/HH))**2)) CYCLE II_LOOP
            ENDIF
            IF (IN%SHAPE=='CYLINDER') THEN
               IF (((XC(II)-X0)**2+(YC(JJ)-Y0)**2<RRI**2) .OR. &
                  ((XC(II)-X0)**2+(YC(JJ)-Y0)**2>RR**2)) CYCLE II_LOOP
            ENDIF
            ! If local XB intersects with cutcell, estimate this intersection volume
            VCX = (MIN(X(II),IN_X2)-MAX(X(II-1),IN_X1))
            VCY = (MIN(Y(JJ),IN_Y2)-MAX(Y(JJ-1),IN_Y1))
            VCZ = (MIN(Z(KK),IN_Z2)-MAX(Z(KK-1),IN_Z1))
            IF (CC_IBM) THEN
               ICC = CCVAR(II,JJ,KK,CC_IDCC)
               ! Approximate intersection as min of the two volumes
               IF (ICC>0) THEN
                  INSERT_VOLUME = INSERT_VOLUME + MIN(VCX*VCY*VCZ,SUM(CUT_CELL(ICC)%VOLUME(:)))
                  FILLED_VOLUME = FILLED_VOLUME + SUM(CUT_CELL(ICC)%VOLUME(:))
               ELSE
                  INSERT_VOLUME = INSERT_VOLUME + VCX*VCY*VCZ
                  FILLED_VOLUME = FILLED_VOLUME + DX(II)*DY(JJ)*DZ(KK)
               ENDIF
            ELSE
               INSERT_VOLUME = INSERT_VOLUME + VCX*VCY*VCZ
               FILLED_VOLUME = FILLED_VOLUME + DX(II)*DY(JJ)*DZ(KK)
            ENDIF
            N_INSERT_CELLS = N_INSERT_CELLS + 1

            INSERT_PARTICLE_LOOP_2: DO IP = 1, IN%N_PARTICLES_PER_CELL
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

               CALL ALLOCATE_STORAGE(NM,LP_INDEX=NLP,LPC_INDEX=ILPC,SURF_INDEX=LPC%SURF_INDEX)

               LAGRANGIAN_PARTICLE => MESHES(NM)%LAGRANGIAN_PARTICLE
               LP => LAGRANGIAN_PARTICLE(NLP)
               LP%CLASS_INDEX = ILPC

               BC => MESHES(NM)%BOUNDARY_COORD(LP%BC_INDEX)

               BC%IOR = 0  ! Particle is not stuck to a wall

               PARTICLE_TAG = PARTICLE_TAG + NMESHES
               LP%TAG = PARTICLE_TAG

               ! Get particle coordinates by randomly choosing within the designated volume

               XC1 = MAX(X1,X(II-1))
               YC1 = MAX(Y1,Y(JJ-1))
               ZC1 = MAX(Z1,Z(KK-1))
               XC2 = MIN(X2,X(II))
               YC2 = MIN(Y2,Y(JJ))
               ZC2 = MIN(Z2,Z(KK))
               LP%DX = XC2 - XC1
               LP%DY = YC2 - YC1
               LP%DZ = ZC2 - ZC1

               IF (IN%CELL_CENTERED) THEN
                  BC%X = 0.5_EB*(X(II-1)+X(II))
                  BC%Y = 0.5_EB*(Y(JJ-1)+Y(JJ))
                  BC%Z = 0.5_EB*(Z(KK-1)+Z(KK))
                  ! If particle goes in a cut cell, move it to the centroid
                  IF (CC_IBM) THEN
                     ICC=CCVAR(II,JJ,KK,CC_IDCC)
                     IF (ICC>0) THEN
                        BC%X = CUT_CELL(ICC)%XYZCEN(IAXIS,1)
                        BC%Y = CUT_CELL(ICC)%XYZCEN(JAXIS,1)
                        BC%Z = CUT_CELL(ICC)%XYZCEN(KAXIS,1)
                     ENDIF
                  ENDIF
               ELSE
                  CALL RANDOM_RECTANGLE(BC%X,BC%Y,BC%Z,XC1,XC2,YC1,YC2,ZC1,ZC2)
                  ! Check for particle inside solid part of cut cell
                  IF (CC_IBM) THEN
                     INDCF = CCVAR(II,JJ,KK,CC_IDCF)
                     IF (INDCF>0) THEN
                        ICC = CCVAR(II,JJ,KK,CC_IDCC)
                        CF => CUT_FACE(INDCF)
                        ! Check dot-product of normal and centroid-particle vector for closest CFACE
                        ! Limited number of location guesses before choosing the cut cell centroid
                        RAND_LOCATION_LOOP: DO I_RAND=1,50
                           DIST_MIN=HUGE_EB
                           CFA_LOOP2: DO IFACE=1,CF%NFACE
                              P_VECTOR = (/BC%X-CF%XYZCEN(IAXIS,IFACE), &
                                          BC%Y-CF%XYZCEN(JAXIS,IFACE), &
                                          BC%Z-CF%XYZCEN(KAXIS,IFACE)/)
                              DIST = NORM2(P_VECTOR)
                              IF (DIST<DIST_MIN) THEN
                                 DIST_MIN = DIST
                                 P_VECTOR_MIN = P_VECTOR
                                 NVEC_MIN = BOUNDARY_COORD(CFACE(CF%CFACE_INDEX(IFACE))%BC_INDEX)%NVEC
                              ENDIF
                           ENDDO CFA_LOOP2
                           IF (DOT_PRODUCT(NVEC_MIN,P_VECTOR_MIN) > TWO_EPSILON_EB) EXIT RAND_LOCATION_LOOP
                           CALL RANDOM_RECTANGLE(BC%X,BC%Y,BC%Z,XC1,XC2,YC1,YC2,ZC1,ZC2)
                        ENDDO RAND_LOCATION_LOOP
                        ! No suitable location was found, move to centroid
                        IF (I_RAND==51) THEN
                           BC%X = CUT_CELL(ICC)%XYZCEN(IAXIS,1)
                           BC%Y = CUT_CELL(ICC)%XYZCEN(JAXIS,1)
                           BC%Z = CUT_CELL(ICC)%XYZCEN(KAXIS,1)
                        ENDIF
                     ENDIF
                  ENDIF

               ENDIF

               CALL VOLUME_INIT_PARTICLE

               LP => LAGRANGIAN_PARTICLE(NLP)

            ENDDO INSERT_PARTICLE_LOOP_2
         ENDDO II_LOOP
      ENDDO
   ENDDO

ENDIF TOTAL_OR_PER_CELL

IF (N_INSERT>0) THEN

   ! Adjust particle weighting factor PWT so that desired MASS_PER_VOLUME is achieved

   IF (IN%MASS_PER_TIME>0._EB) THEN
      VOLUME_SPLIT_FACTOR = 1._EB
      IF (INPUT_VOLUME>TWO_EPSILON_EB .AND. INSERT_VOLUME>TWO_EPSILON_EB) VOLUME_SPLIT_FACTOR = INSERT_VOLUME/INPUT_VOLUME
      PWT0 = VOLUME_SPLIT_FACTOR*RAMP_FACTOR*IN%MASS_PER_TIME*IN%DT_INSERT/MASS_SUM
   ELSEIF (IN%MASS_PER_VOLUME>0._EB) THEN
      PWT0 = RAMP_FACTOR*IN%MASS_PER_VOLUME*INSERT_VOLUME/MASS_SUM
   ELSE
      PWT0 = IN%PARTICLE_WEIGHT_FACTOR
   ENDIF

   DO IIP=1,MIN(MAXIMUM_PARTICLES,N_INSERT)
      IP = LP_INDEX_LOOKUP(IIP)
      LP => LAGRANGIAN_PARTICLE(IP)
      ! Scale particles according to contribution to total FILLED_VOLUME (e.g. stretched grid or cut cells)
      IF (IN%MASS_PER_VOLUME>0._EB .AND. IN%N_PARTICLES_PER_CELL > 0) THEN
         LP%PWT = LP%PWT*PWT0*N_INSERT_CELLS/FILLED_VOLUME/LP%RVC
      ELSE
         LP%PWT = LP%PWT*PWT0
      ENDIF
      IF (ANY(IN%PATH_RAMP_INDEX>0)) LP%PATH_PARTICLE=.TRUE.
   ENDDO

ENDIF

DEALLOCATE(LP_INDEX_LOOKUP)

IF (.NOT.IN%INVOKED_BY_SURF .AND. IN%BULK_DENSITY_FILE=='null') IN%ALREADY_INSERTED(NM) = .TRUE.

END SUBROUTINE INSERT_VOLUMETRIC_PARTICLES


!> \brief Place holder for inserting particles for duct heat transfer

SUBROUTINE INSERT_DUCT_PARTICLES
IF (DUCT_HT_INSERTED) RETURN
DUCT_HT_INSERTED = .TRUE.
END SUBROUTINE INSERT_DUCT_PARTICLES


!> \brief Initialize particle indices and velocity

SUBROUTINE VOLUME_INIT_PARTICLE

USE OUTPUT_DATA, ONLY: N_PROF
INTEGER :: ND
TYPE (PROFILE_TYPE), POINTER :: PF

BC%IIG = II; BC%JJG = JJ; BC%KKG = KK
BC%II  = II; BC%JJ  = JJ; BC%KK  = KK
LP%U = IN%U0
LP%V = IN%V0
LP%W = IN%W0

! If the INITIALIZATION group has an ID, match it with a DEVC (device) and/or PROFile

IF (IN%ID/='null') THEN

   DO ND=1,N_DEVC
      DV => DEVICE(ND)
      IF ((IN%ID==DV%INIT_ID .AND. IP==DV%POINT) .OR. DV%LP_TAG==PARTICLE_TAG) THEN
         IF (DV%LP_TAG>0 .AND. DV%LP_TAG/=PARTICLE_TAG) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: INIT_ID on DEVC ',TRIM(DV%ID),' cannot have more than one particle'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.)
         ENDIF
         DV%LP_TAG = PARTICLE_TAG
         DV%PART_CLASS_INDEX = ILPC
         DV%MESH = NM
         DV%X = BC%X
         DV%Y = BC%Y
         DV%Z = BC%Z
         IF (DV%LINE>0 .AND. DV%LINE_COORD_CODE==123) THEN
            IF (ABS(IN%DX)>TWO_EPSILON_EB .AND. ABS(IN%DY)<TWO_EPSILON_EB .AND. ABS(IN%DZ)<TWO_EPSILON_EB) DV%LINE_COORD_CODE = 1
            IF (ABS(IN%DX)<TWO_EPSILON_EB .AND. ABS(IN%DY)>TWO_EPSILON_EB .AND. ABS(IN%DZ)<TWO_EPSILON_EB) DV%LINE_COORD_CODE = 2
            IF (ABS(IN%DX)<TWO_EPSILON_EB .AND. ABS(IN%DY)<TWO_EPSILON_EB .AND. ABS(IN%DZ)>TWO_EPSILON_EB) DV%LINE_COORD_CODE = 3
         ENDIF
      ENDIF
   ENDDO

   DO ND=1,N_PROF
      PF => PROFILE(ND)
      IF (IN%ID==PF%INIT_ID .OR. PF%LP_TAG==PARTICLE_TAG) THEN
         IF (PF%LP_TAG>0 .AND. PF%LP_TAG/=PARTICLE_TAG) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: INIT_ID on PROF ',TRIM(PF%ID),' cannot have more than one particle'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.)
         ENDIF
         PF%LP_TAG = PARTICLE_TAG
         PF%MESH = NM
         PF%X = BC%X
         PF%Y = BC%Y
         PF%Z = BC%Z
      ENDIF
   ENDDO

ENDIF

! Process particle and set more initial values

CALL INITIALIZE_SINGLE_PARTICLE

! Special over-ride of particle diameter for liquid droplets only

IF (IN%DIAMETER>0._EB) THEN
   LP%RADIUS = 0.5_EB*IN%DIAMETER
   LP%PWT = 1._EB
   LP%MASS = FOTHPI*LPC%DENSITY*LP%RADIUS**3
ENDIF

! Save insert time and other miscellaneous attributes

LP%T_INSERT = T

! Get the particle ORIENTATION from the PART line

IF (LPC%N_ORIENTATION>0) LP%ORIENTATION_INDEX = LPC%ORIENTATION_INDEX

! Add property info for 'ADIABATIC SURFACE TEMPERATURE GAS'

IF (IN%ID/='null') THEN
   DO ND=1,N_DEVC
      DV => DEVICE(ND)
      IF (IN%ID==DV%INIT_ID) THEN ! the INIT line is referred to be the DEVC line
         IF (DV%QUANTITY(1)=='RADIATIVE HEAT FLUX' .OR. &
             DV%QUANTITY(1)=='GAUGE HEAT FLUX' .OR. &
             DV%QUANTITY(1)=='RADIANCE' .OR. &
             DV%QUANTITY(1)=='RADIOMETER' .OR. &
             DV%QUANTITY(1)=='ADIABATIC SURFACE TEMPERATURE') THEN
            IF (LPC%ID=='RESERVED TARGET PARTICLE') THEN  ! use the orientation of the DEVC
               LP%ORIENTATION_INDEX = DV%ORIENTATION_INDEX
               LP%PROP_INDEX = DV%PROP_INDEX
               LP%PWT = 1._EB
            ENDIF
            IF (DV%PROP_INDEX>0) THEN
               B1 => MESHES(NM)%BOUNDARY_PROP1(LP%B1_INDEX)
               B1%EMISSIVITY = PROPERTY(DV%PROP_INDEX)%EMISSIVITY
               IF (RADIATION) THEN
                  B1%Q_RAD_OUT  = PROPERTY(DV%PROP_INDEX)%EMISSIVITY*SIGMA*TMPA4
               ELSE
                  B1%Q_RAD_OUT  = 0._EB
               ENDIF
               IF (PROPERTY(DV%PROP_INDEX)%HEAT_TRANSFER_COEFFICIENT>0._EB) &
                  B1%HEAT_TRANS_COEF = PROPERTY(DV%PROP_INDEX)%HEAT_TRANSFER_COEFFICIENT
            ENDIF
         ENDIF
         EXIT
      ENDIF
   ENDDO
ENDIF

! Sum up the particle masses

MASS_SUM = MASS_SUM + LP%PWT*LP%MASS ! if r=0 the sum will stay 0

END SUBROUTINE VOLUME_INIT_PARTICLE


!> \brief Set up the properties of a single, newly inserted particle

SUBROUTINE INITIALIZE_SINGLE_PARTICLE
USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
REAL(EB) :: AREA,SCALE_FACTOR,RADIUS,LP_VOLUME
INTEGER :: N
TYPE(BOUNDARY_ONE_D_TYPE), POINTER :: LP_ONE_D
TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1
TYPE(SURFACE_TYPE), POINTER :: LP_SF

LP%PWT  = 1._EB
LP%RE = 0._EB
LP%C_DRAG = 0._EB

CALL RANDOM_NUMBER(RN)
IF (RN < 1._EB/REAL(LPC%SAMPLING_FACTOR,EB)) THEN
   LP%SHOW = .TRUE.
ELSE
   LP%SHOW = .FALSE.
ENDIF

IF (LPC%MASSLESS_TRACER) RETURN

LP_SF => SURFACE(LPC%SURF_INDEX)
B1 => MESHES(NM)%BOUNDARY_PROP1(LP%B1_INDEX)

IF (LPC%SOLID_PARTICLE) THEN

   IF (LPC%SURF_INDEX==TGA_SURF_INDEX) THEN
      TGA_PARTICLE_INDEX = NLP
      TGA_MESH_INDEX = NM
   ENDIF

   LP%MASS = 0._EB
   SCALE_FACTOR = 1._EB
   SELECT CASE (LP_SF%GEOMETRY)
      CASE(SURF_CARTESIAN)
         LP%RADIUS = LP_SF%THICKNESS
      CASE(SURF_CYLINDRICAL)
         LP%RADIUS = LP_SF%THICKNESS + LP_SF%INNER_RADIUS
      CASE(SURF_SPHERICAL)
         LP%RADIUS = LP_SF%THICKNESS + LP_SF%INNER_RADIUS
   END SELECT

   SELECT CASE (LPC%DRAG_LAW)

      CASE(SCREEN_DRAG)

         ! Compute special cross-sectional area of screen particle

         AREA = (ABS(ORIENTATION_VECTOR(1,LPC%ORIENTATION_INDEX))*DY(BC%JJG)*DZ(BC%KKG) + &
                 ABS(ORIENTATION_VECTOR(2,LPC%ORIENTATION_INDEX))*DX(BC%IIG)*DZ(BC%KKG) + &
                 ABS(ORIENTATION_VECTOR(3,LPC%ORIENTATION_INDEX))*DX(BC%IIG)*DY(BC%JJG)) * &
                 (1._EB-LAGRANGIAN_PARTICLE_CLASS(LP%CLASS_INDEX)%FREE_AREA_FRACTION)

         SELECT CASE (LP_SF%GEOMETRY)
            CASE (SURF_CARTESIAN)
               B1%AREA = 2._EB*AREA
               LP%LENGTH = SQRT(AREA)
               IF (LP_SF%THERMAL_BC_INDEX==THERMALLY_THICK) LP%MASS = LPC%INITIAL_MASS*B1%AREA
            CASE (SURF_CYLINDRICAL)
               B1%AREA = AREA*PI
               LP%LENGTH = 0.5_EB*AREA/LP%RADIUS
               IF (LP_SF%THERMAL_BC_INDEX==THERMALLY_THICK) LP%MASS = LPC%INITIAL_MASS*PI*LP%LENGTH
            CASE (SURF_SPHERICAL)
               B1%AREA = 4._EB*PI*LP%RADIUS**2._EB
               LP%PWT = 4._EB*AREA/B1%AREA
               IF (LP_SF%THERMAL_BC_INDEX==THERMALLY_THICK) LP%MASS = LPC%INITIAL_MASS
         END SELECT

      CASE (POROUS_DRAG)

         LP_VOLUME = LPC%POROUS_VOLUME_FRACTION*LP%DX*LP%DY*LP%DZ

         SELECT CASE (LP_SF%GEOMETRY)
            CASE (SURF_CARTESIAN)
               B1%AREA = LP_VOLUME/LP%RADIUS
               LP%LENGTH = SQRT(0.5_EB*B1%AREA)
               IF (LP_SF%THERMAL_BC_INDEX==THERMALLY_THICK) LP%MASS = LPC%INITIAL_MASS*B1%AREA
            CASE (SURF_CYLINDRICAL)
               B1%AREA = 2._EB*LP_VOLUME/LP%RADIUS
               LP%LENGTH = LP_VOLUME/(PI*LP%RADIUS**2)
               IF (LP_SF%THERMAL_BC_INDEX==THERMALLY_THICK) LP%MASS = LPC%INITIAL_MASS*PI*LP%LENGTH
            CASE (SURF_SPHERICAL)
               B1%AREA = 4._EB*PI*LP%RADIUS**2._EB
               LP%PWT = 3._EB*LP_VOLUME/(LP%RADIUS*B1%AREA)
               IF (LP_SF%THERMAL_BC_INDEX==THERMALLY_THICK) LP%MASS = LPC%INITIAL_MASS
         END SELECT

      CASE DEFAULT

         IF (.NOT.LPC%MONODISPERSE .AND. LPC%SURF_INDEX/=TGA_SURF_INDEX .AND. LP%OD_INDEX>0) THEN
            LP_ONE_D => MESHES(NM)%BOUNDARY_ONE_D(LP%OD_INDEX)
            CALL PARTICLE_SIZE_WEIGHT(RADIUS,LP%PWT)
            SCALE_FACTOR = RADIUS/LP_SF%THICKNESS
            LP_ONE_D%X(:) = LP_ONE_D%X(:)*SCALE_FACTOR
            LP_ONE_D%LAYER_THICKNESS(:) = LP_ONE_D%LAYER_THICKNESS(:)*SCALE_FACTOR
         ENDIF

         LP%RADIUS = LP%RADIUS*SCALE_FACTOR
         SELECT CASE (LP_SF%GEOMETRY)
            CASE (SURF_CARTESIAN)
               B1%AREA = 2._EB*LP_SF%LENGTH*LP_SF%WIDTH
            CASE (SURF_CYLINDRICAL)
               B1%AREA = 2._EB*PI*LP%RADIUS*LP_SF%LENGTH*SCALE_FACTOR
               SCALE_FACTOR=SCALE_FACTOR**2._EB
            CASE (SURF_SPHERICAL)
               B1%AREA = 4._EB*PI*LP%RADIUS**2._EB*SCALE_FACTOR**2._EB
               SCALE_FACTOR=SCALE_FACTOR**3._EB
         END SELECT

         IF (LP_SF%THERMAL_BC_INDEX==THERMALLY_THICK) LP%MASS = LPC%INITIAL_MASS*SCALE_FACTOR

   END SELECT

ELSEIF (LPC%LIQUID_DROPLET) THEN

   CALL PARTICLE_SIZE_WEIGHT(LP%RADIUS,LP%PWT)
   LP%MASS = FOTHPI*LPC%DENSITY*LP%RADIUS**3

ENDIF

! Set the initial inner temperature for all particles

IF (LP%OD_INDEX>0) THEN
   LP_ONE_D => MESHES(NM)%BOUNDARY_ONE_D(LP%OD_INDEX)
   IF (LPC%TMP_INITIAL>0._EB) THEN
      LP_ONE_D%TMP = LPC%TMP_INITIAL
      B1%TMP_F = LP_ONE_D%TMP(1)
      B1%TMP_B = LP_ONE_D%TMP(1)
   ELSEIF (LP_SF%THERMAL_BC_INDEX==THERMALLY_THICK) THEN
      IF (LP_SF%RAMP_T_I_INDEX > 0) THEN
         B1%TMP_F = EVALUATE_RAMP(0._EB,LP_SF%RAMP_T_I_INDEX)
         DO N=1,LP_SF%N_CELLS_INI
            LP_ONE_D%TMP(N) = EVALUATE_RAMP(LP_ONE_D%X(N),LP_SF%RAMP_T_I_INDEX)
         ENDDO
         B1%TMP_B = EVALUATE_RAMP(SUM(SF%LAYER_THICKNESS),LP_SF%RAMP_T_I_INDEX)
         LP_ONE_D%TMP(0) = 2._EB*B1%TMP_F - LP_ONE_D%TMP(0)
         LP_ONE_D%TMP(LP_SF%N_CELLS_INI+1) = 2._EB*B1%TMP_F - LP_ONE_D%TMP(LP_SF%N_CELLS_INI)
      ELSE
         IF(LP_SF%TMP_INNER>0._EB) THEN
            LP_ONE_D%TMP = LP_SF%TMP_INNER
         ELSE
            LP_ONE_D%TMP = TMPA
         ENDIF
         B1%TMP_F = LP_ONE_D%TMP(1)
         B1%TMP_B = LP_ONE_D%TMP(1)
      ENDIF
   ELSE
      LP_ONE_D%TMP = TMPA
      B1%TMP_F = LP_ONE_D%TMP(1)
      B1%TMP_B = LP_ONE_D%TMP(1)
   ENDIF
ELSE
   IF (LPC%TMP_INITIAL>0._EB) THEN
      B1%TMP_F = LPC%TMP_INITIAL
      B1%TMP_B = LPC%TMP_INITIAL
   ELSE
      B1%TMP_F = TMPA
      B1%TMP_B = TMPA
   ENDIF
ENDIF

! Check if fire spreads radially over this surface type, and if so, set T_IGN appropriately

IF (LP_SF%FIRE_SPREAD_RATE>0._EB) THEN
   B1%T_IGN = T_BEGIN + SQRT((BC%X-LP_SF%XYZ(1))**2 + (BC%Y-LP_SF%XYZ(2))**2 + (BC%Z-LP_SF%XYZ(3))**2)/LP_SF%FIRE_SPREAD_RATE
ELSE  ! Assume that the ignition time is relative to the time of insertion, not the time of the simulation starting.
   B1%T_IGN = T + (LP_SF%T_IGN-T_BEGIN)
ENDIF

! Store volume over which gas-solid exchanges occur

IF (BC%X>XS .AND. BC%X<XF .AND. BC%Y>YS .AND. BC%Y<YF .AND. BC%Z>ZS .AND. BC%Z<ZF) CALL GET_RVC(NM,BC%IIG,BC%JJG,BC%KKG,LP%RVC)

END SUBROUTINE INITIALIZE_SINGLE_PARTICLE


!> \brief Determine a particles's size and weight

SUBROUTINE PARTICLE_SIZE_WEIGHT(R,PWT)
REAL(EB), INTENT(OUT):: R,PWT

IF (LPC%MONODISPERSE) THEN
   R   = 0.5_EB*LPC%DIAMETER
   PWT = 1._EB
ELSE
   CALL RANDOM_NUMBER(RN)
   STRATUM = NINT(REAL(LPC%N_STRATA,EB)*REAL(RN,EB)+0.5_EB)
   IL = LPC%STRATUM_INDEX_LOWER(STRATUM)
   IU = LPC%STRATUM_INDEX_UPPER(STRATUM)
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


END SUBROUTINE INSERT_ALL_PARTICLES


!> \brief Update particle position over one gas phase time step

SUBROUTINE MOVE_PARTICLES(T,DT,NM)

USE TRAN, ONLY: GET_IJK
USE COMPLEX_GEOMETRY, ONLY: CC_CGSC,CC_FGSC,CC_IDCF,CC_GASPHASE,CC_SOLID,CC_CUTCFE,POINT_IN_CFACE
USE CC_SCALARS, ONLY: CUTFACE_VELOCITIES
USE MATH_FUNCTIONS, ONLY: CROSS_PRODUCT
USE PHYSICAL_FUNCTIONS, ONLY: EMBER_IGNITION_MODEL
INTEGER :: IFACE,ICF_OLD,ICF_NEW,INDCF,ICF_MIN,IND1,IND2
REAL(EB) :: DIST,DIST_MIN,VEL_VECTOR_1(3),VEL_VECTOR_2(3),P_VECTOR(3),TNOW,&
            DOT_NVECNEW_GVEC,DOT_NVECOLD_GVEC,DOT_NVECOLD_NEW
REAL(EB), INTENT(IN) :: T,DT
INTEGER, INTENT(IN) :: NM
REAL     :: RN
REAL(EB) :: XI,YJ,ZK,R_D,R_D_0,X_OLD,Y_OLD,Z_OLD,X_TRY,Y_TRY,Z_TRY,THETA,THETA_RN,STEP_FRACTION(-3:3),DT_CFL,DT_P,&
            STEP_FRACTION_PREVIOUS,DELTA,PVEC_L
LOGICAL :: HIT_SOLID,CC_CC_GASPHASE,EXTRACT_PARTICLE
INTEGER :: IP,IC_NEW,IIG_OLD,JJG_OLD,KKG_OLD,IIG_TRY,JJG_TRY,KKG_TRY,IW,IC_OLD,IOR_HIT,&
           N_ITER,ITER,I_COORD,J_COORD,IC_TRY,IOR_ORIGINAL,ICC
TYPE (LAGRANGIAN_PARTICLE_TYPE), POINTER :: LP
TYPE (LAGRANGIAN_PARTICLE_CLASS_TYPE), POINTER :: LPC
TYPE (SURFACE_TYPE), POINTER :: SF
TYPE(CFACE_TYPE), POINTER :: CFA_OLD,CFA_NEW
TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC,CFA_OLD_BC,CFA_NEW_BC
REAL(EB), PARAMETER :: ONTHHALF=0.5_EB**ONTH, B_1=1.7321_EB
LOGICAL :: TEST_POS, BOUNCE_CF, IN_CFACE, SLIDE_CF, EXT_CFACE
INTEGER :: DIND, MADD(3,3)
INTEGER, PARAMETER :: EYE3(1:3,1:3)=RESHAPE( (/1,0,0, 0,1,0, 0,0,1 /), (/3,3/) )

! Zero out max particle velocity if CFL number is to be bound by particle speed.

PART_UVWMAX = 0._EB

IF (MESHES(NM)%NLP==0) RETURN

! Set the CPU timer and point to the current mesh variables

TNOW=CURRENT_TIME()
CALL POINT_TO_MESH(NM)

IF (PARTICLE_DRAG) THEN
   FVX_D = 0._EB
   FVY_D = 0._EB
   FVZ_D = 0._EB
ENDIF

IF (CC_IBM) CALL CUTFACE_VELOCITIES(NM,U,V,W,CUTFACES=.TRUE.)

! Loop through all Lagrangian particles and move them one time step

PARTICLE_LOOP: DO IP=1,NLP

   ! Assign particle (LP%) and particle class (LPC%) shortcuts

   LP  => LAGRANGIAN_PARTICLE(IP)
   LPC => LAGRANGIAN_PARTICLE_CLASS(LP%CLASS_INDEX)
   SF  => SURFACE(LPC%SURF_INDEX)
   BC  => BOUNDARY_COORD(LP%BC_INDEX)

   ! Determine the limiting time step (DT_P) to ensure particle does not traverse more than a single grid cell

   DT_CFL = MIN(DX(BC%IIG)/(ABS(LP%U)+TWO_EPSILON_EB),DY(BC%JJG)/(ABS(LP%V)+TWO_EPSILON_EB),DZ(BC%KKG)/(ABS(LP%W)+TWO_EPSILON_EB))
   N_ITER = CEILING(DT/(0.90_EB*DT_CFL))
   DT_P   = DT/REAL(N_ITER,EB)

   ! Initialize drag terms

   LP%ACCEL_X = 0._EB
   LP%ACCEL_Y = 0._EB
   LP%ACCEL_Z = 0._EB

   ! Save value of IOR to determine if the particle hit any surface during the time step

   IOR_ORIGINAL = BC%IOR

   ! Note if Particle was attached to an external CFACE:
   EXT_CFACE = .FALSE.
   IF(LP%CFACE_INDEX==EXTERNAL_CFACE) THEN
      EXT_CFACE = .TRUE.; LP%CFACE_INDEX = 0
   ENDIF

   ! Sub-timesteps

   TIME_STEP_LOOP: DO ITER=1,N_ITER

      ! Determine particle radius

      IF (LPC%MASSLESS_TARGET) THEN
         ! If there is a path ramp, move the target
         IF (LP%PATH_PARTICLE) THEN
            CALL MOVE_IN_GAS
            CALL GET_IJK(BC%X,BC%Y,BC%Z,NM,XI,YJ,ZK,BC%IIG,BC%JJG,BC%KKG)
            BC%II = BC%IIG; BC%JJ = BC%JJG; BC%KK = BC%KKG
         ENDIF
         EXIT TIME_STEP_LOOP
      ENDIF
      IF (LPC%LIQUID_DROPLET .OR. LPC%SOLID_PARTICLE) THEN
         R_D = LP%RADIUS
      ELSEIF (LPC%MASSLESS_TRACER) THEN
         R_D = 0._EB
      ENDIF

      ! Throw out particles that have run out of mass.

      IF (.NOT.LPC%MASSLESS_TRACER .AND. R_D<LPC%KILL_RADIUS) EXIT TIME_STEP_LOOP

      ! Save original particle radius.

      R_D_0 = R_D

      ! Determine the current coordinates of the particle

      CALL GET_IJK(BC%X,BC%Y,BC%Z,NM,XI,YJ,ZK,IIG_OLD,JJG_OLD,KKG_OLD)

      IC_OLD = CELL_INDEX(IIG_OLD,JJG_OLD,KKG_OLD)

      X_OLD = BC%X
      Y_OLD = BC%Y
      Z_OLD = BC%Z

      ! Throw out particles that are inside a solid obstruction unless they are following a path

      IF ((CELL(IC_OLD)%SOLID .OR. CELL(IC_OLD)%EXTERIOR) .AND. .NOT. LP%PATH_PARTICLE ) THEN
         BC%X = 1.E6_EB
         EXIT TIME_STEP_LOOP
      ENDIF

      ! Move the particle one sub-time-step

      SOLID_GAS_MOVE: IF (BC%IOR/=0 .AND. LPC%ADHERE_TO_SOLID) THEN

         CALL MOVE_ON_SOLID

      ELSE SOLID_GAS_MOVE

         CALL MOVE_IN_GAS

         ! If the particle is massless or does not move, go on to the next particle

         IF (LPC%MASSLESS_TRACER .OR. LP%PWT<=TWO_EPSILON_EB .OR. (LPC%STATIC .AND. .NOT.LP%EMBER)) EXIT TIME_STEP_LOOP

      ENDIF SOLID_GAS_MOVE

      ! Determine the cell indices of the new particle location.

      CALL GET_IJK(BC%X,BC%Y,BC%Z,NM,XI,YJ,ZK,BC%IIG,BC%JJG,BC%KKG)
      BC%II = BC%IIG; BC%JJ = BC%JJG; BC%KK = BC%KKG

      ! If the particle is not near a boundary cell, cycle.

      CC_CC_GASPHASE = .TRUE.
      IF (CC_IBM) THEN
         IF (CCVAR(BC%IIG,BC%JJG,BC%KKG,CC_CGSC)/=CC_GASPHASE) &
            CC_CC_GASPHASE = .FALSE.
      ENDIF

      IC_NEW = CELL_INDEX(BC%IIG,BC%JJG,BC%KKG)

      IF ((IC_OLD==0 .OR. IC_NEW==0) .AND. CC_CC_GASPHASE) THEN
         BC%IOR = 0
         LP%CFACE_INDEX = 0
         CYCLE TIME_STEP_LOOP
      ENDIF

      ! The next big loops are to determine if the particle or droplet has hit a solid surface.
      ! HIT_SOLID indicates that it has.

      HIT_SOLID = .FALSE.
      IW = 0

      ! Determine if the particle is near a CFACE, and if so, change its trajectory

      CFACE_SEARCH: IF (CC_IBM) THEN

         ICF_OLD   = LP%CFACE_INDEX; ICF_NEW   = 0
         INDCF     = CCVAR(BC%IIG,BC%JJG,BC%KKG,CC_IDCF)
         BOUNCE_CF = .TRUE.

         IF ( INDCF < 1 ) THEN
            PARTICLE_POSITION_IF: IF (CCVAR(BC%IIG,BC%JJG,BC%KKG,CC_CGSC)==CC_SOLID) THEN
               ! Kinematics of a surface particle moving on Horizontal GEOM surface and passing to CC_SOLID cell.
               ! Bounce back on random direction, maintaining CFACE_INDEX:

               IF (ICF_OLD /= 0 .AND. ABS(LP%W)<TWO_EPSILON_EB) THEN ! Particle already attached to ICF_OLD.
                  CALL RANDOM_NUMBER(RN)
                  DIST_MIN = (1._EB-SIGN(1._EB,LP%V))*PI/2._EB
                  IF(ABS(LP%U) > TWO_EPSILON_EB) DIST_MIN = ATAN2(LP%V,LP%U)
                  THETA_RN = PI*(REAL(RN,EB)+0.5_EB)+DIST_MIN
                  VEL_VECTOR_1(1) = COS(THETA_RN)
                  VEL_VECTOR_1(2) = SIN(THETA_RN)
                  VEL_VECTOR_1(3) = 0._EB
                  BC%X=X_OLD; BC%Y=Y_OLD; BC%Z=Z_OLD
                  LP%U = VEL_VECTOR_1(1)*LPC%HORIZONTAL_VELOCITY
                  LP%V = VEL_VECTOR_1(2)*LPC%HORIZONTAL_VELOCITY
                  LP%W = VEL_VECTOR_1(3)*LPC%VERTICAL_VELOCITY
                  BOUNCE_CF = .FALSE.
               ELSE
                  ! Search for cut-cell in the direction of -GVEC:
                  DIND = MAXLOC(ABS(GVEC(1:3)),DIM=1); MADD(1:3,1:3) = -INT(SIGN(1._EB,GVEC(DIND)))*EYE3
                  INDCF = CCVAR(BC%IIG+MADD(1,DIND),BC%JJG+MADD(2,DIND),BC%KKG+MADD(3,DIND),CC_IDCF) ! One cell up
                  IF(INDCF<1) THEN ! Two cells up
                     MADD(1:3,1:3) = -2*INT(SIGN(1._EB,GVEC(DIND)))*EYE3
                     INDCF = CCVAR(BC%IIG+MADD(1,DIND),BC%JJG+MADD(2,DIND),BC%KKG+MADD(3,DIND),CC_IDCF)
                  ENDIF
               ENDIF
            ELSEIF(CCVAR(BC%IIG,BC%JJG,BC%KKG,CC_CGSC)==CC_CUTCFE .AND. ICF_OLD/=0) THEN PARTICLE_POSITION_IF
               ! Kinematics of Particle moving horizontaly and falling over side walls, or falling of side walls and
               ! moving underneath of object.
               CFA_OLD => CFACE(ICF_OLD)
               CFA_OLD_BC => BOUNDARY_COORD(CFA_OLD%BC_INDEX)
               DOT_NVECOLD_GVEC = DOT_PRODUCT(CFA_OLD_BC%NVEC,GVEC/(NORM2(GVEC)+TWO_EPSILON_EB))
               IF(DOT_NVECOLD_GVEC<-0.99_EB) THEN ! ICF_OLD looking up, search for side walls in the max GVEC component direction.
                  ! Assume Particle has fallen of a side. Search for cut-cell in the direction of GVEC:
                  DIND = MAXLOC(ABS(GVEC(1:3)),DIM=1); MADD(1:3,1:3) = INT(SIGN(1._EB,GVEC(DIND)))*EYE3
                  INDCF = CCVAR(BC%IIG+MADD(1,DIND),BC%JJG+MADD(2,DIND),BC%KKG+MADD(3,DIND),CC_IDCF)
                  IF(INDCF<1) THEN
                     MADD(1:3,1:3) = 2*INT(SIGN(1._EB,GVEC(DIND)))*EYE3
                     INDCF = CCVAR(BC%IIG+MADD(1,DIND),BC%JJG+MADD(2,DIND),BC%KKG+MADD(3,DIND),CC_IDCF)
                  ENDIF
               ELSEIF(ABS(DOT_NVECOLD_GVEC)<0.01_EB) THEN ! ICF_OLD lookinf sideways, search in the X-Y directions for plane walls:
                  ! Search for cut-cell in the direction of max X-Y
                  LOOP_X : DO I_COORD=-1,1
                     DO J_COORD=-1,1
                        MADD(3,3)=ABS(I_COORD)+ABS(J_COORD); IF(ANY(MADD(3,3)==(/0,2/))) CYCLE
                        INDCF = CCVAR(BC%IIG+I_COORD,BC%JJG+J_COORD,BC%KKG,CC_IDCF)
                        IF(INDCF>0) THEN
                           DO IFACE=1,CUT_FACE(INDCF)%NFACE  ! Loop through CFACEs and find the one closest to the particle
                              ICF_NEW = CUT_FACE(INDCF)%CFACE_INDEX(IFACE)
                              CFA_NEW => CFACE(ICF_NEW)
                              CFA_NEW_BC => BOUNDARY_COORD(CFA_NEW%BC_INDEX)
                              DOT_NVECNEW_GVEC = DOT_PRODUCT(CFA_NEW_BC%NVEC,GVEC/(NORM2(GVEC)+TWO_EPSILON_EB))
                              IF(DOT_NVECNEW_GVEC>0.99_EB) THEN ! ICF_NEW looking down, particle movin on side wall.
                                 BOUNCE_CF = .FALSE.
                                 EXIT LOOP_X
                              ENDIF
                           ENDDO
                        ENDIF
                     ENDDO
                  ENDDO LOOP_X
                  IF(BOUNCE_CF) INDCF=-1
                  BOUNCE_CF = .TRUE.
               ENDIF
            ENDIF PARTICLE_POSITION_IF
         ENDIF

         ! Check if there are CFACEs in this cell:
         INDCF_POS : IF ( INDCF > 0 ) THEN  ! Current grid cell has CFACEs
            SLIDE_CF  = .FALSE.
            ! Drop if cell with boundary cut-faces found is outside mesh, no CFACEs defined.
            IF (CELL(CELL_INDEX(CUT_FACE(INDCF)%IJK(1),CUT_FACE(INDCF)%IJK(2),CUT_FACE(INDCF)%IJK(3)))%EXTERIOR) THEN
               IF(ICF_OLD>0) LP%CFACE_INDEX = EXTERNAL_CFACE; CYCLE PARTICLE_LOOP
            ENDIF

            ICF_MIN = 0
            ! First test same slope faces:
            IF(ICF_OLD>0) THEN
               CFA_OLD    => CFACE(ICF_OLD); CFA_OLD_BC => BOUNDARY_COORD(CFA_OLD%BC_INDEX)
               DO IFACE=1,CUT_FACE(INDCF)%NFACE
                  ICF_NEW = CUT_FACE(INDCF)%CFACE_INDEX(IFACE)
                  CFA_NEW => CFACE(ICF_NEW); CFA_NEW_BC => BOUNDARY_COORD(CFA_NEW%BC_INDEX)
                  IF(DOT_PRODUCT(CFA_OLD_BC%NVEC,CFA_NEW_BC%NVEC)<=0.99_EB) CYCLE
                  CALL POINT_IN_CFACE(NM,BC%X,BC%Y,BC%Z,ICF_NEW,IN_CFACE);
                  IF(IN_CFACE) THEN; ICF_MIN = ICF_NEW; EXIT; ENDIF
               ENDDO
            ENDIF
            ! Other Faces:
            IF(ICF_MIN==0) THEN
               DIST_MIN = 1.E6_EB;
               DO IFACE=1,CUT_FACE(INDCF)%NFACE  ! Loop through CFACEs and find the one closest to the particle
                  ICF_NEW = CUT_FACE(INDCF)%CFACE_INDEX(IFACE)
                  CALL POINT_IN_CFACE(NM,BC%X,BC%Y,BC%Z,ICF_NEW,IN_CFACE);
                  IF(IN_CFACE) THEN
                     ICF_MIN = ICF_NEW; EXIT
                  ENDIF
                  CFA_NEW => CFACE(ICF_NEW); CFA_NEW_BC => BOUNDARY_COORD(CFA_NEW%BC_INDEX)
                  DIST = (BC%X-CFA_NEW_BC%X)**2 + (BC%Y-CFA_NEW_BC%Y)**2 + (BC%Z-CFA_NEW_BC%Z)**2
                  IF (DIST<DIST_MIN) THEN
                     DIST_MIN = DIST
                     ICF_MIN = ICF_NEW
                  ENDIF
               ENDDO
            ENDIF
            ICF_NEW = ICF_MIN

            ! Now test case of ICF_NEW same as previous CFACE_INDEX and it's normal opposite to GVEC:
            IF (ICF_OLD>0 .AND. ICF_NEW>0) THEN

               ! Old CFACE:
               CFA_OLD    => CFACE(ICF_OLD)
               CFA_OLD_BC => BOUNDARY_COORD(CFA_OLD%BC_INDEX)
               ! New CFACE:
               CFA_NEW => CFACE(ICF_NEW)
               CFA_NEW_BC => BOUNDARY_COORD(CFA_NEW%BC_INDEX)

               DOT_NVECOLD_GVEC= DOT_PRODUCT(CFA_OLD_BC%NVEC,GVEC/(NORM2(GVEC)+TWO_EPSILON_EB))
               DOT_NVECOLD_NEW = DOT_PRODUCT(CFA_OLD_BC%NVEC,CFA_NEW_BC%NVEC)
               IN_CFACE=.TRUE.
               ! Test for case of CFACE_INDEX switching to an ICF_NEW with similar slope, if so don't do anything.
               IF (.NOT.(ICF_OLD/=ICF_NEW .AND. DOT_NVECOLD_NEW>0.99_EB) ) THEN

                  IF (ICF_OLD/=ICF_NEW .AND. DOT_NVECOLD_NEW<=0.99_EB) THEN
                     ! Case of switching to different ICF_NEW with different slope:
                     DOT_NVECNEW_GVEC = DOT_PRODUCT(CFA_NEW_BC%NVEC,GVEC/(NORM2(GVEC)+TWO_EPSILON_EB))
                     IF(DOT_NVECOLD_GVEC<-0.99_EB) THEN ! Only for slopes less than 8 degrees. CFACE looking up.
                        ! If ICF_NEW is almost vertical pointing in the direction of velocity, set creep velocity.
                        IF(ABS(DOT_NVECNEW_GVEC)<0.01_EB .AND. (CFA_NEW_BC%Z<CFA_OLD_BC%Z)) IN_CFACE=.FALSE.
                     ELSEIF(ABS(DOT_NVECOLD_GVEC)<0.01_EB) THEN ! Side walls almost 90 deg resp to gravity.
                        ! Switching from side wall to ICF_NEW pointing downwards, depends on underside conditions.
                        IF(ABS(DOT_NVECNEW_GVEC)>0.99_EB .AND. (CFA_NEW_BC%Z<CFA_OLD_BC%Z)) IN_CFACE=.FALSE.
                     ELSE
                        LP%CFACE_INDEX = 0
                        BC%IOR = 0
                        SLIDE_CF  = .TRUE.
                     ENDIF
                  ELSE ! ICF_OLD==ICF_NEW.
                     ! Here we test if LP lays outside of CFACE ICF_OLD polygon, being ICF_OLD a 'top' CFACE.
                     ! If so, attach to the most vertical CFACE in surrounding cells (drop in the sides).
                     CALL POINT_IN_CFACE(NM,BC%X,BC%Y,BC%Z,ICF_OLD,IN_CFACE)
                     IF(.NOT.IN_CFACE)THEN
                        ! Select another CFACE in the cell:
                        DIST_MIN = 0.5_EB*MAX(DX(BC%IIG),DY(BC%JJG),DZ(BC%KKG)); ICF_MIN=0
                        DO IFACE=1,CUT_FACE(INDCF)%NFACE  ! Loop through CFACEs and find the one closest to the particle
                           ICF_NEW = CUT_FACE(INDCF)%CFACE_INDEX(IFACE); IF(ICF_OLD==ICF_NEW) CYCLE
                           CFA_NEW => CFACE(ICF_NEW); CFA_NEW_BC => BOUNDARY_COORD(CFA_NEW%BC_INDEX)
                           DIST = SQRT((BC%X-CFA_NEW_BC%X)**2 + (BC%Y-CFA_NEW_BC%Y)**2 + (BC%Z-CFA_NEW_BC%Z)**2)
                           IF (DIST<DIST_MIN) THEN
                              DIST_MIN = DIST
                              ICF_MIN = ICF_NEW
                           ENDIF
                        ENDDO
                        IF (ICF_MIN/=0) THEN ! We found a different CFACE.
                           ICF_NEW = ICF_MIN
                           CFA_NEW => CFACE(ICF_NEW); CFA_NEW_BC => BOUNDARY_COORD(CFA_NEW%BC_INDEX)
                           ! Do not fix velocity if found CFACE is not at lower height than CFACE_INDEX one.
                           IF ((CFA_NEW_BC%Z>(CFA_OLD_BC%Z-TWO_EPSILON_EB))) IN_CFACE = .TRUE.
                        ELSE ! CFACE not found, continue with CFACE_INDEX face if close, else drop.
                           DIST = SQRT((BC%X-CFA_OLD_BC%X)**2 + (BC%Y-CFA_OLD_BC%Y)**2 + (BC%Z-CFA_OLD_BC%Z)**2)
                           IF (DIST>0.2_EB*MAX(DX(BC%IIG),DY(BC%JJG),DZ(BC%KKG)) .AND. &
                               DOT_PRODUCT(CFA_OLD_BC%NVEC,(/LP%U,LP%V,LP%W/))>-TWO_EPSILON_EB) THEN
                              LP%CFACE_INDEX = 0
                              BC%IOR = 0
                              SLIDE_CF  = .TRUE.
                           ELSE
                              ICF_NEW = ICF_OLD
                              IN_CFACE = .TRUE.
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDIF

                  IF(.NOT.IN_CFACE .AND. LP%CFACE_INDEX>0)THEN
                     IF(DOT_NVECOLD_GVEC<-0.99_EB .AND. ABS(DOT_NVECNEW_GVEC)<0.01_EB) THEN !Slope less than 8 deg. CFACE up.
                        ! ICF_NEW found almost vertical pointing in the direction of velocity, set creep velocity.
                        LP%U = 0._EB
                        LP%V = 0._EB
                        LP%W = SIGN(1._EB,GVEC(3))*LPC%VERTICAL_VELOCITY
                        LP%CFACE_INDEX = ICF_NEW
                        BC%IOR = 1
                        HIT_SOLID = .TRUE.
                        SLIDE_CF  = .TRUE.
                     ELSEIF(ABS(DOT_NVECOLD_GVEC)<0.01_EB) THEN ! Side walls almost 90 deg resp to gravity.
                        ! Switching from side wall to ICF_NEW pointing downwards, depends on underside conditions.
                        IF (LPC%SOLID_PARTICLE) THEN
                           LP%CFACE_INDEX = 0
                           BC%IOR = 0
                           SLIDE_CF  = .TRUE.
                        ELSEIF (.NOT.SURFACE(CFACE(ICF_NEW)%SURF_INDEX)%ALLOW_UNDERSIDE_PARTICLES) THEN
                           CFA_NEW => CFACE(ICF_NEW); CFA_NEW_BC => BOUNDARY_COORD(CFA_NEW%BC_INDEX)
                           DOT_NVECOLD_NEW = DOT_PRODUCT(CFA_OLD_BC%NVEC,CFA_NEW_BC%NVEC)
                           IF (DOT_NVECOLD_NEW>0.99_EB) THEN ! ICF_NEW also vertical, particle movin on side wall.
                              LP%CFACE_INDEX = ICF_NEW
                              BC%IOR = 1
                              HIT_SOLID = .TRUE.
                              SLIDE_CF  = .TRUE.
                           ELSE ! ICF_NEW is tilted respect to ICF_OLD, particle falls down slope.
                              CALL RANDOM_NUMBER(RN)
                              BC%X = BC%X + REAL(RN,EB)*(CFA_NEW_BC%X-BC%X)
                              BC%Y = BC%Y + REAL(RN,EB)*(CFA_NEW_BC%Y-BC%Y)
                              BC%Z = BC%Z + REAL(RN,EB)*(CFA_NEW_BC%Z-BC%Z)
                              LP%U = 0._EB
                              LP%V = 0._EB
                              LP%W = -LPC%VERTICAL_VELOCITY
                              LP%CFACE_INDEX = 0
                              BC%IOR = 0
                              SLIDE_CF  = .TRUE.
                           ENDIF
                        ELSE
                           DOT_NVECOLD_NEW = DOT_PRODUCT(CFA_OLD_BC%NVEC,CFA_NEW_BC%NVEC)
                           IF (DOT_NVECOLD_NEW>0.99_EB) THEN ! OLD and NEW CFACEs are aligned.
                              LP%CFACE_INDEX = ICF_NEW
                              BC%IOR = 1
                              HIT_SOLID = .TRUE.
                              SLIDE_CF  = .TRUE.
                           ELSE ! OLD and NEW CFACEs are not aligned, ICF_OLD 90 deg resp to gravity.
                              CFA_NEW => CFACE(ICF_NEW); CFA_NEW_BC => BOUNDARY_COORD(CFA_NEW%BC_INDEX)
                              CALL RANDOM_NUMBER(RN);
                              BC%X = BC%X + REAL(RN,EB)*(CFA_NEW_BC%X-BC%X)
                              BC%Y = BC%Y + REAL(RN,EB)*(CFA_NEW_BC%Y-BC%Y)
                              BC%Z = BC%Z + REAL(RN,EB)*(CFA_NEW_BC%Z-BC%Z)
                              THETA_RN = TWOPI*REAL(RN,EB)
                              CALL CROSS_PRODUCT(VEL_VECTOR_2,CFA_NEW_BC%NVEC,GVEC)
                              IF(NORM2(VEL_VECTOR_2)<TWO_EPSILON_EB) THEN ! Use a random direction in the CFACE plane.
                                 IND1 = CFA_NEW%CUT_FACE_IND1
                                 IND2 = CFA_NEW%CUT_FACE_IND2
                                 VEL_VECTOR_2 = CUT_FACE(IND1)%XYZVERT(1:3,CUT_FACE(IND1)%CFELEM(1,IND2))-&
                                                CUT_FACE(IND1)%XYZVERT(1:3,CUT_FACE(IND1)%CFELEM(2,IND2))
                              ENDIF
                              VEL_VECTOR_2 = VEL_VECTOR_2/NORM2(VEL_VECTOR_2)
                              CALL CROSS_PRODUCT(VEL_VECTOR_1,VEL_VECTOR_2,CFA_NEW_BC%NVEC)
                              ! The surface is tilted; particles go down slope:
                              LP%U = (VEL_VECTOR_1(1)*COS(THETA_RN)+VEL_VECTOR_1(2)*SIN(THETA_RN))*LPC%HORIZONTAL_VELOCITY
                              LP%V = (VEL_VECTOR_2(1)*COS(THETA_RN)+VEL_VECTOR_2(2)*SIN(THETA_RN))*LPC%HORIZONTAL_VELOCITY
                              LP%W = (CFA_NEW_BC%NVEC(1)*COS(THETA_RN)+CFA_NEW_BC%NVEC(2)*SIN(THETA_RN))*LPC%VERTICAL_VELOCITY
                              LP%CFACE_INDEX = ICF_NEW
                              BC%IOR = 1
                              HIT_SOLID = .TRUE.
                              SLIDE_CF  = .TRUE.
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF

            ENDIF


            SLIDE_CF_IF : IF (.NOT.SLIDE_CF) THEN ! Avoid this block in case droplet falls on an obstacle side (filled before).

               CFA_NEW => CFACE(ICF_NEW)
               CFA_NEW_BC => BOUNDARY_COORD(CFA_NEW%BC_INDEX)
               DOT_NVECNEW_GVEC = DOT_PRODUCT(CFA_NEW_BC%NVEC,GVEC/(NORM2(GVEC)+TWO_EPSILON_EB))

               ! If the CFACE normal points up, force the particle to follow the contour. If the normal points down,
               ! put the particle back into the gas phase.

               P_VECTOR = (/BC%X-CFA_NEW_BC%X, BC%Y-CFA_NEW_BC%Y, BC%Z-CFA_NEW_BC%Z/) ! NEW CFACE to particle position vector.
               PVEC_L   = NORM2(P_VECTOR)
               TEST_POS = .FALSE.; IF (ICF_OLD == 0) TEST_POS = DOT_PRODUCT(CFA_NEW_BC%NVEC,P_VECTOR) > TWO_EPSILON_EB
               DIST = TWO_EPSILON_EB; IF (SURFACE(CFACE(ICF_NEW)%SURF_INDEX)%ALLOW_UNDERSIDE_PARTICLES) DIST = 1._EB
               CFACE_ATTACH : IF (DOT_NVECNEW_GVEC>DIST .OR. TEST_POS) THEN

                  ! Normal points down or particle in gas phase.
                  ! If particle is coming from another mesh, attach to ICF_NEW.
                  IF(EXT_CFACE .AND. BC%IOR>0) THEN
                     LP%CFACE_INDEX = ICF_NEW
                     HIT_SOLID = .TRUE.
                  ENDIF
                  ! If solid particle let particle move freely.
                  IF (LPC%SOLID_PARTICLE) THEN
                     LP%CFACE_INDEX = 0
                     BC%IOR = 0
                     HIT_SOLID = .FALSE.
                  ENDIF

               ELSE  CFACE_ATTACH ! normal points up; determine direction for particle to move

                  IF(LPC%ADHERE_TO_SOLID) THEN
                     CALL CROSS_PRODUCT(VEL_VECTOR_1,CFA_NEW_BC%NVEC,GVEC)
                     CALL CROSS_PRODUCT(VEL_VECTOR_2,VEL_VECTOR_1,CFA_NEW_BC%NVEC)
                     CFACE_SLOPE : IF (NORM2(VEL_VECTOR_2) > TWO_EPSILON_EB .AND. ABS(LP%W) > TWO_EPSILON_EB) THEN
                        ! The surface is tilted; particles go down slope:
                        THETA_RN = ACOS(DOT_NVECNEW_GVEC)
                        DIST_MIN= ABS(LPC%HORIZONTAL_VELOCITY*COS(THETA_RN)**2+LPC%VERTICAL_VELOCITY*SIN(THETA_RN)**2)
                        VEL_VECTOR_1 = VEL_VECTOR_2/NORM2(VEL_VECTOR_2)
                        LP%U = VEL_VECTOR_1(1)*DIST_MIN
                        LP%V = VEL_VECTOR_1(2)*DIST_MIN
                        LP%W = VEL_VECTOR_1(3)*DIST_MIN
                     ELSEIF (BC%IOR==0) THEN  CFACE_SLOPE ! surface is flat and particle has no direction,
                                                          ! particle is given random direction
                        CALL RANDOM_NUMBER(RN)
                        THETA_RN = TWOPI*REAL(RN,EB)
                        VEL_VECTOR_1(IAXIS) = COS(THETA_RN)
                        VEL_VECTOR_1(JAXIS) = SIN(THETA_RN)
                        VEL_VECTOR_1(KAXIS) = 0._EB
                        LP%U = VEL_VECTOR_1(IAXIS)*LPC%HORIZONTAL_VELOCITY
                        LP%V = VEL_VECTOR_1(JAXIS)*LPC%HORIZONTAL_VELOCITY
                        LP%W = VEL_VECTOR_1(KAXIS)*LPC%VERTICAL_VELOCITY
                     ELSEIF (ABS(DOT_NVECNEW_GVEC)<0.99_EB .AND. ICF_OLD/=0 .AND. ABS(LP%W)<TWO_EPSILON_EB) THEN CFACE_SLOPE
                        ! Particle moving in horizontal direction and assumed crossing into solid.
                        ! Bounce back on random direction, maintaining CFACE_INDEX:
                        IF (DOT_PRODUCT( (/ LP%U, LP%V /) ,CFA_NEW_BC%NVEC(IAXIS:JAXIS))<-TWO_EPSILON_EB)THEN
                           CALL RANDOM_NUMBER(RN)
                           DIST_MIN = (1._EB-SIGN(1._EB,LP%V))*PI/2._EB
                           IF(ABS(LP%U) > TWO_EPSILON_EB) DIST_MIN = ATAN2(LP%V,LP%U)
                           THETA_RN = PI*(REAL(RN,EB)+0.5_EB)+DIST_MIN
                           VEL_VECTOR_1(1) = COS(THETA_RN)
                           VEL_VECTOR_1(2) = SIN(THETA_RN)
                           VEL_VECTOR_1(3) = 0._EB
                           BC%X=X_OLD; BC%Y=Y_OLD; BC%Z=Z_OLD
                           LP%U = VEL_VECTOR_1(1)*LPC%HORIZONTAL_VELOCITY
                           LP%V = VEL_VECTOR_1(2)*LPC%HORIZONTAL_VELOCITY
                           LP%W = VEL_VECTOR_1(3)*LPC%VERTICAL_VELOCITY
                           EXIT TIME_STEP_LOOP
                        ENDIF
                     ENDIF CFACE_SLOPE
                  ELSE
                      LP%U = 0._EB
                      LP%V = 0._EB
                      LP%W = 0._EB
                  ENDIF

                  ! If the particle is inside the solid, move it to the surface in the normal direction.

                  PVEC_L = NORM2(P_VECTOR)
                  IF (PVEC_L>TWO_EPSILON_EB) THEN
                     THETA = ACOS(DOT_PRODUCT(CFA_NEW_BC%NVEC,P_VECTOR/PVEC_L))
                     IF (THETA>PIO2) THEN
                        DELTA = PVEC_L*SIN(THETA-0.5_EB*PI)+TWO_EPSILON_EB
                        BC%X = BC%X + DELTA*CFA_NEW_BC%NVEC(1)
                        BC%Y = BC%Y + DELTA*CFA_NEW_BC%NVEC(2)
                        BC%Z = BC%Z + DELTA*CFA_NEW_BC%NVEC(3)
                     ENDIF
                  ENDIF
                  LP%CFACE_INDEX = ICF_NEW
                  BC%IOR = 1
                  HIT_SOLID = .TRUE.

                  CALL VENT_PARTICLE_EXTRACTION(HIT_SOLID,EXTRACT_PARTICLE,CFACE_INDEX=ICF_NEW)
                  IF (EXTRACT_PARTICLE) EXIT TIME_STEP_LOOP

               ENDIF CFACE_ATTACH

            ENDIF SLIDE_CF_IF

         ELSEIF (CCVAR(BC%IIG,BC%JJG,BC%KKG,CC_CGSC)/=CC_GASPHASE .AND. BOUNCE_CF) THEN INDCF_POS

            IF(ICF_OLD/=0) THEN
               CFA_OLD => CFACE(ICF_OLD);  CFA_OLD_BC => BOUNDARY_COORD(CFA_OLD%BC_INDEX)
               IF (DOT_PRODUCT(CFA_OLD_BC%NVEC,GVEC/(NORM2(GVEC)+TWO_EPSILON_EB))>0.99_EB) THEN
                  ! Leaving an under side CFACE, set particle to move freely and invert X-Y velocities:
                  LP%U = -2._EB*LP%U
                  LP%V = -2._EB*LP%V
               ENDIF
            ENDIF
            BC%IOR = 0
            LP%CFACE_INDEX = 0

         ENDIF INDCF_POS

      ENDIF CFACE_SEARCH

      ! If the particle crosses a cell boundary, determine its new status and check if it has hit a solid.

      WALL_SEARCH: IF (.NOT.HIT_SOLID .AND. LP%CFACE_INDEX==0 .AND. &
                       (IIG_OLD/=BC%IIG .OR. JJG_OLD/=BC%JJG .OR. KKG_OLD/=BC%KKG)) THEN

         ! Calculate the STEP_FRACTION, which indicates the relative distance between the particles's old and new
         ! position where the particle hits a cell boundary.

         STEP_FRACTION = 1.1_EB

         IF (BC%IIG>IIG_OLD) STEP_FRACTION(-1) = (X(IIG_OLD)  -X_OLD)/(BC%X-X_OLD)
         IF (BC%IIG<IIG_OLD) STEP_FRACTION( 1) = (X(IIG_OLD-1)-X_OLD)/(BC%X-X_OLD)
         IF (BC%JJG>JJG_OLD) STEP_FRACTION(-2) = (Y(JJG_OLD)  -Y_OLD)/(BC%Y-Y_OLD)
         IF (BC%JJG<JJG_OLD) STEP_FRACTION( 2) = (Y(JJG_OLD-1)-Y_OLD)/(BC%Y-Y_OLD)
         IF (BC%KKG>KKG_OLD) STEP_FRACTION(-3) = (Z(KKG_OLD)  -Z_OLD)/(BC%Z-Z_OLD)
         IF (BC%KKG<KKG_OLD) STEP_FRACTION( 3) = (Z(KKG_OLD-1)-Z_OLD)/(BC%Z-Z_OLD)

         ! The minimum value of STEP_FRACTION indicates the relative location along the particle path where it first crosses
         ! a cell boundary. Test this location to see if the cell the particle crosses into is solid. If it is, indicate that
         ! that the particle has HIT_SOLID and EXIT. If not, cycle through the other cell boundary crossings. The particle can
         ! cross at most three cell boundaries in one sub-timestep.

         STEP_FRACTION_PREVIOUS = -1000000._EB
         IIG_TRY = IIG_OLD
         JJG_TRY = JJG_OLD
         KKG_TRY = KKG_OLD
         IW = 0

         TRIAL_LOOP: DO I_COORD=1,3
            IOR_HIT = MINLOC(STEP_FRACTION,DIM=1,MASK=STEP_FRACTION>STEP_FRACTION_PREVIOUS) - 4
            IF (STEP_FRACTION(IOR_HIT)>1._EB) EXIT TRIAL_LOOP
            X_TRY = X_OLD + STEP_FRACTION(IOR_HIT)*(BC%X-X_OLD)
            Y_TRY = Y_OLD + STEP_FRACTION(IOR_HIT)*(BC%Y-Y_OLD)
            Z_TRY = Z_OLD + STEP_FRACTION(IOR_HIT)*(BC%Z-Z_OLD)
            IC_TRY = CELL_INDEX(IIG_TRY,JJG_TRY,KKG_TRY)
            IW = CELL(IC_TRY)%WALL_INDEX(-IOR_HIT)
            IF (WALL(IW)%BOUNDARY_TYPE==SOLID_BOUNDARY) THEN
               LP%WALL_INDEX = IW
               HIT_SOLID = .TRUE.
               CALL VENT_PARTICLE_EXTRACTION(HIT_SOLID,EXTRACT_PARTICLE,WALL_INDEX=IW)
               IF (.NOT. HIT_SOLID) EXIT TRIAL_LOOP
               IF (EXTRACT_PARTICLE) EXIT TIME_STEP_LOOP
               BC%IOR  = IOR_HIT
               BC%X = X_TRY
               BC%Y = Y_TRY
               BC%Z = Z_TRY
               SELECT CASE(IOR_HIT)
                  CASE(-3) ; BC%Z = BC%Z - 0.01*DZ(KKG_TRY)
                  CASE(-2) ; BC%Y = BC%Y - 0.01*DY(JJG_TRY)
                  CASE(-1) ; BC%X = BC%X - 0.01*DX(IIG_TRY)
                  CASE( 1) ; BC%X = BC%X + 0.01*DX(IIG_TRY)
                  CASE( 2) ; BC%Y = BC%Y + 0.01*DY(JJG_TRY)
                  CASE( 3) ; BC%Z = BC%Z + 0.01*DZ(KKG_TRY)
               END SELECT
               BC%IIG = IIG_TRY
               BC%JJG = JJG_TRY
               BC%KKG = KKG_TRY
               IC_NEW = IC_TRY
               EXIT TRIAL_LOOP
            ENDIF
            STEP_FRACTION_PREVIOUS = STEP_FRACTION(IOR_HIT)
            SELECT CASE(IOR_HIT)
               CASE(-3) ; KKG_TRY = KKG_TRY+1
               CASE(-2) ; JJG_TRY = JJG_TRY+1
               CASE(-1) ; IIG_TRY = IIG_TRY+1
               CASE( 1) ; IIG_TRY = IIG_TRY-1
               CASE( 2) ; JJG_TRY = JJG_TRY-1
               CASE( 3) ; KKG_TRY = KKG_TRY-1
            END SELECT
         ENDDO TRIAL_LOOP

         ! If the particle has hit a Cartesian solid, choose a new direction

         IF_HIT_SOLID: IF (HIT_SOLID .AND. IW>0) THEN

            DIRECTION: SELECT CASE(BC%IOR)
               CASE (-2:-1,1:2) DIRECTION
                  LP%U = 0._EB
                  LP%V = 0._EB
                  IF (LPC%ADHERE_TO_SOLID) LP%W = -LPC%VERTICAL_VELOCITY
               CASE (-3) DIRECTION
                  IF (LPC%SOLID_PARTICLE) THEN
                     BC%IOR = 0
                  ELSEIF (.NOT.SURFACE(WALL(IW)%SURF_INDEX)%ALLOW_UNDERSIDE_PARTICLES) THEN
                     LP%U = 0._EB
                     LP%V = 0._EB
                     LP%W = -LPC%VERTICAL_VELOCITY
                     BC%IOR = 0
                  ELSE
                     CALL RANDOM_NUMBER(RN)
                     THETA_RN = TWOPI*REAL(RN,EB)
                     LP%U = LPC%HORIZONTAL_VELOCITY*COS(THETA_RN)
                     LP%V = LPC%HORIZONTAL_VELOCITY*SIN(THETA_RN)
                     LP%W = 0._EB
                  ENDIF
               CASE (3) DIRECTION
                  IF (LPC%ADHERE_TO_SOLID) THEN
                     CALL RANDOM_NUMBER(RN)
                     THETA_RN = TWOPI*REAL(RN,EB)
                     LP%U = LPC%HORIZONTAL_VELOCITY*COS(THETA_RN)
                     LP%V = LPC%HORIZONTAL_VELOCITY*SIN(THETA_RN)
                  ELSE
                     LP%U = 0._EB
                     LP%V = 0._EB
                     IF (LEVEL_SET_MODE>0 .AND. LP%EMBER) CALL EMBER_IGNITION_MODEL(NM,IP)
                  ENDIF
                  LP%W = 0._EB
            END SELECT DIRECTION

         ELSEIF (HIT_SOLID .AND. ICF_NEW>0) THEN

            IF (LEVEL_SET_MODE>0 .AND. LP%EMBER) CALL EMBER_IGNITION_MODEL(NM,IP)

         ENDIF IF_HIT_SOLID

      ENDIF WALL_SEARCH

      ! If the particle has passed outside of its current mesh and it has not
      ! hit anything, schedule it for removal or adoption by another mesh.

      IF (CELL(IC_NEW)%EXTERIOR .AND. .NOT.HIT_SOLID) CYCLE PARTICLE_LOOP

      ! Process the particle if it has hit either a WALL or CFACE

      IF (HIT_SOLID) THEN

         ! Remove the particle if it is not allowed on a surface

         IF (IW>0 .AND. .NOT.SURFACE(WALL(IW)%SURF_INDEX)%ALLOW_SURFACE_PARTICLES) THEN
            LP%RADIUS = 0.9_EB*LPC%KILL_RADIUS
            EXIT TIME_STEP_LOOP
         ENDIF

         ! Adjust the size of the PARTICLE and weighting factor

         IF (LPC%LIQUID_DROPLET .AND. LPC%SURFACE_DIAMETER>0._EB) THEN
            R_D = MIN(0.5_EB*LPC%SURFACE_DIAMETER,LP%PWT**ONTH*R_D)
            LP%PWT = LP%PWT*(R_D_0/R_D)**3
            LP%RADIUS = R_D
            LP%MASS = FOTHPI*LPC%DENSITY*R_D**3
         ENDIF

      ENDIF

      ! If the droplet was attached to a solid WALL (BC%IOR/=0), but now it is not, change its course.
      ! If the droplet was
      ! dripping down a vertical surface (IOR=+-1,2), make it go under the solid and then move upward to (possibly) stick
      ! to the underside or drip off. If the droplet moves off an upward or downward facing horizontal surface (IOR=+-3), reverse
      ! its course and drop it down the side of the solid obstruction.

      IF (LP%CFACE_INDEX==0 .AND. BC%IOR/=0) THEN

         LP%WALL_INDEX = CELL(IC_NEW)%WALL_INDEX(-BC%IOR)

         IF (WALL(LP%WALL_INDEX)%BOUNDARY_TYPE/=SOLID_BOUNDARY) THEN
            IF (LPC%ADHERE_TO_SOLID) THEN
               SELECT CASE(BC%IOR)
                  CASE( 1)
                     BC%X = BC%X - 0.2_EB*DX(BC%IIG)
                     LP%W = SQRT(2._EB*GRAV*DZ(BC%KKG))
                  CASE(-1)
                     BC%X = BC%X + 0.2_EB*DX(BC%IIG)
                     LP%W = SQRT(2._EB*GRAV*DZ(BC%KKG))
                  CASE( 2)
                     BC%Y = BC%Y - 0.2_EB*DY(BC%JJG)
                     LP%W = SQRT(2._EB*GRAV*DZ(BC%KKG))
                  CASE(-2)
                     BC%Y = BC%Y + 0.2_EB*DY(BC%JJG)
                     LP%W = SQRT(2._EB*GRAV*DZ(BC%KKG))
                  CASE(-3,3)
                     LP%U = -2._EB*LP%U
                     LP%V = -2._EB*LP%V
                     BC%Z =  BC%Z - 0.2_EB*DZ(BC%KKG)
               END SELECT
            ENDIF
            BC%IOR = 0
            LP%WALL_INDEX = 0
         ENDIF

      ELSEIF (LP%CFACE_INDEX==0) THEN

         LP%WALL_INDEX = 0  ! The droplet is not stuck to a WALL cell

      ENDIF

      ! Store containing volume at new location
      CALL GET_IJK(BC%X,BC%Y,BC%Z,NM,XI,YJ,ZK,BC%IIG,BC%JJG,BC%KKG)
      BC%II = BC%IIG; BC%JJ = BC%JJG; BC%KK = BC%KKG
      CALL GET_RVC(NM,BC%IIG,BC%JJG,BC%KKG,LP%RVC)

   ENDDO TIME_STEP_LOOP

   ! If the particle is not stuck to a wall, allow it to be counted again.

   IF (BC%IOR==0 .AND. IOR_ORIGINAL==0) LP%SPLAT = .FALSE.

ENDDO PARTICLE_LOOP

IF (CC_IBM) CALL CUTFACE_VELOCITIES(NM,U,V,W,CUTFACES=.FALSE.)

! Remove out-of-bounds particles

CALL REMOVE_PARTICLES(T,NM)

T_USED(8)=T_USED(8)+CURRENT_TIME()-TNOW

CONTAINS


!> \brief Move particles attached to solid surfaces

SUBROUTINE MOVE_ON_SOLID

BC%X = X_OLD + LP%U*DT_P
BC%Y = Y_OLD + LP%V*DT_P
BC%Z = Z_OLD + LP%W*DT_P

END SUBROUTINE MOVE_ON_SOLID


!> \brief Move particles in the gas phase

SUBROUTINE MOVE_IN_GAS

USE PHYSICAL_FUNCTIONS, ONLY : DRAG, GET_VISCOSITY, SURFACE_DENSITY
USE MATH_FUNCTIONS, ONLY : EVALUATE_RAMP, RANDOM_CHOICE, BOX_MULLER
USE SOOT_ROUTINES, ONLY: DROPLET_SCRUBBING
REAL(EB) :: UBAR,VBAR,WBAR,UREL,VREL,WREL,QREL,RHO_G,TMP_G,MU_FILM, &
            U_OLD,V_OLD,W_OLD,ZZ_GET(1:N_TRACKED_SPECIES),WAKE_VEL,DROP_VOL_FRAC,RE_WAKE,&
            WE_G,T_BU_BAG,T_BU_STRIP,MPOM,SFAC,BREAKUP_RADIUS(0:NDC),&
            DD,DD_X,DD_Y,DD_Z,DW_X,DW_Y,DW_Z,K_TERM(3),Y_TERM(3),C_DRAG,A_DRAG,&
            GX_LOC,GY_LOC,GZ_LOC,DRAG_MAX(3)=0._EB,K_SGS,U_P,KN,M_DOT,&
            EMBER_DENSITY,EMBER_VOLUME=0._EB,ACCEL_X,ACCEL_Y,ACCEL_Z,&
            LP_FORCE,FACE_VOLS(2,2,2),VEL_G_INT(3),VOL_WGT(2,2,2),&
            EMBER_PACKING_RATIO,LOCAL_PACKING_RATIO,LPC_GEOM_FACTOR
REAL(EB) :: WGT(2,2,2,3),VEL_G(2,2,2)
REAL(EB), POINTER, DIMENSION(:,:,:) :: FV_D
REAL(EB), SAVE :: BETA
INTEGER :: IIX,JJY,KKZ,IL,JL,KL,AXIS,N_LPC2
LOGICAL :: STUCK=.FALSE.
TYPE(LAGRANGIAN_PARTICLE_CLASS_TYPE), POINTER :: LPC2

! Save current values of particle velocity components

U_OLD = LP%U
V_OLD = LP%V
W_OLD = LP%W

! Interpolate the nearest velocity components of the gas

! Lower cell index for nearest 8 interpolation points
IIX = FLOOR(XI+.5_EB)
JJY = FLOOR(YJ+.5_EB)
KKZ = FLOOR(ZK+.5_EB)
! For cut cells base nearest points on cell centroid
ICC = 0
IF (CC_IBM) THEN
   ICC = CCVAR(IIG_OLD,JJG_OLD,KKG_OLD,CC_IDCC)
   IF (ICC>0) THEN
      IF (BC%X<=CUT_CELL(ICC)%XYZCEN(IAXIS,1) .AND. BC%X>XC(IIG_OLD)) IIX = IIX-1
      IF (BC%X>=CUT_CELL(ICC)%XYZCEN(IAXIS,1) .AND. BC%X<XC(IIG_OLD)) IIX = IIX+1
      IF (BC%Y<=CUT_CELL(ICC)%XYZCEN(JAXIS,1) .AND. BC%Y>YC(JJG_OLD)) JJY = JJY-1
      IF (BC%Y>=CUT_CELL(ICC)%XYZCEN(JAXIS,1) .AND. BC%Y<YC(JJG_OLD)) JJY = JJY+1
      IF (BC%Z<=CUT_CELL(ICC)%XYZCEN(KAXIS,1) .AND. BC%Z>ZC(KKG_OLD)) KKZ = KKZ-1
      IF (BC%Z>=CUT_CELL(ICC)%XYZCEN(KAXIS,1) .AND. BC%Z<ZC(KKG_OLD)) KKZ = KKZ+1
   ENDIF
ENDIF

WGT=0._EB
DO AXIS=IAXIS,KAXIS
   IL = IIX; JL = JJY; KL = KKZ
   IF (AXIS==IAXIS) THEN
      IL = FLOOR(XI)
      VEL_G = U(IL:IL+1,JL:JL+1,KL:KL+1)
   ELSEIF (AXIS==JAXIS) THEN
      JL = FLOOR(YJ)
      VEL_G = V(IL:IL+1,JL:JL+1,KL:KL+1)
   ELSEIF (AXIS==KAXIS) THEN
      KL = FLOOR(ZK)
      VEL_G = W(IL:IL+1,JL:JL+1,KL:KL+1)
   ENDIF
   IF (ICC>0) THEN
      CALL GET_FACE_IDW(AXIS,IL,JL,KL,BC%X,BC%Y,BC%Z,WGT(:,:,:,AXIS))
   ELSE
      CALL GET_FACE_TLW(AXIS,IL,JL,KL,BC%X,BC%Y,BC%Z,WGT(:,:,:,AXIS),VEL_G)
   ENDIF
   VEL_G_INT(AXIS) = SUM(VEL_G*WGT(:,:,:,AXIS))
ENDDO
UBAR = VEL_G_INT(IAXIS)
VBAR = VEL_G_INT(JAXIS)
WBAR = VEL_G_INT(KAXIS)

! If the particle has a path, just follow the path and return

IF (LP%PATH_PARTICLE) THEN
   IF (INITIALIZATION(LP%INIT_INDEX)%PATH_RAMP_INDEX(1) > 0) &
      BC%X = EVALUATE_RAMP(T,INITIALIZATION(LP%INIT_INDEX)%PATH_RAMP_INDEX(1))
   IF (INITIALIZATION(LP%INIT_INDEX)%PATH_RAMP_INDEX(2) > 0) &
      BC%Y = EVALUATE_RAMP(T,INITIALIZATION(LP%INIT_INDEX)%PATH_RAMP_INDEX(2))
   IF (INITIALIZATION(LP%INIT_INDEX)%PATH_RAMP_INDEX(3) > 0) &
      BC%Z = EVALUATE_RAMP(T,INITIALIZATION(LP%INIT_INDEX)%PATH_RAMP_INDEX(3))
   RETURN
ENDIF

! If the particle is massless, just move it and go on to the next particle

TRACER_IF: IF (LPC%MASSLESS_TRACER .OR. LP%PWT<=TWO_EPSILON_EB) THEN
   IF (LPC%TURBULENT_DISPERSION) THEN
      DD_X = RSC_T * (MU(IIG_OLD+1,JJG_OLD,KKG_OLD) - MU(IIG_OLD-1,JJG_OLD,KKG_OLD)) * &
             RDXN(IIG_OLD-1)*RDXN(IIG_OLD)/(RDXN(IIG_OLD-1) + RDXN(IIG_OLD))
      DD_Y = RSC_T * (MU(IIG_OLD,JJG_OLD+1,KKG_OLD) - MU(IIG_OLD,JJG_OLD-1,KKG_OLD)) * &
             RDYN(JJG_OLD-1)*RDYN(JJG_OLD)/(RDYN(JJG_OLD-1) + RDYN(JJG_OLD))
      DD_Z = RSC_T * (MU(IIG_OLD,JJG_OLD,KKG_OLD+1) - MU(IIG_OLD,JJG_OLD,KKG_OLD-1)) * &
             RDZN(KKG_OLD-1)*RDZN(KKG_OLD)/(RDZN(KKG_OLD-1) + RDZN(KKG_OLD))
      LP%U = UBAR + DD_X/RHO(IIG_OLD,JJG_OLD,KKG_OLD)
      LP%V = VBAR + DD_Y/RHO(IIG_OLD,JJG_OLD,KKG_OLD)
      LP%W = WBAR + DD_Z/RHO(IIG_OLD,JJG_OLD,KKG_OLD)
      DD   = SQRT(2._EB*MU(IIG_OLD,JJG_OLD,KKG_OLD)/RHO(IIG_OLD,JJG_OLD,KKG_OLD)*RSC_T*DT_P)
      ! generate pairs of standard Gaussian random variables
      CALL BOX_MULLER(DW_X,DW_Y)
      CALL BOX_MULLER(DW_Z,DW_X)
      BC%X = X_OLD + LP%U*DT_P + DD*DW_X
      BC%Y = Y_OLD + LP%V*DT_P + DD*DW_Y
      BC%Z = Z_OLD + LP%W*DT_P + DD*DW_Z
   ELSE
      LP%U = UBAR
      LP%V = VBAR
      LP%W = WBAR
      BC%X = X_OLD + LP%U*DT_P
      BC%Y = Y_OLD + LP%V*DT_P
      BC%Z = Z_OLD + LP%W*DT_P
   ENDIF
   IF (TWO_D) THEN
      BC%Y=Y_OLD
      LP%V=0._EB
   ENDIF
   RETURN
ENDIF TRACER_IF

! Point to surface properties of the particle

B1  => BOUNDARY_PROP1(LP%B1_INDEX)

! Massive particles undergoing turbulent dispersion (EXPERIMENTAL, under construction)

IF (LPC%TURBULENT_DISPERSION) THEN
   ! Built model from ideas in this ref:
   ! M. Breuer and M. Alletto, Efficient simulation of particle-laden turbulent flows with high mass loadings using LES,
   ! Int. J. Heat and Fluid Flow, 35:2-12, 2012.
   ! The basic idea is to add an isotropic turbulent fluctuation to the cell mean velocity components prior to
   ! computing the drag.
   DELTA = LES_FILTER_WIDTH(IIG_OLD,JJG_OLD,KKG_OLD)
   K_SGS = (MU(IIG_OLD,JJG_OLD,KKG_OLD)/RHO(IIG_OLD,JJG_OLD,KKG_OLD)/C_DEARDORFF/DELTA)**2 ! def of Deardorff eddy viscosity
   U_P = SQRT(TWTH*K_SGS)
   CALL BOX_MULLER(DW_X,DW_Y)
   CALL BOX_MULLER(DW_Z,DW_X)
   UBAR = UBAR + U_P*DW_X
   VBAR = VBAR + U_P*DW_Y
   WBAR = WBAR + U_P*DW_Z
ENDIF

IF (LPC%EMBER_PARTICLE) THEN
   SELECT CASE(SF%GEOMETRY)
      CASE(SURF_CARTESIAN)
         EMBER_VOLUME = SF%LENGTH * SF%WIDTH * 2._EB*R_D
      CASE(SURF_CYLINDRICAL)
         EMBER_VOLUME = SF%LENGTH * PI*R_D**2
      CASE(SURF_SPHERICAL)
         EMBER_VOLUME = FOTHPI * R_D**3
   END SELECT

   ! experimental ember snag model
   IF (LP%EMBER .AND. LPC%EMBER_SNAG_FACTOR>0._EB) THEN
      STUCK=.FALSE.
      ! constrain packing ratio between 0 and 1
      EMBER_PACKING_RATIO = MAX(MIN(1._EB, EMBER_VOLUME * LP%PWT * LP%RVC),0._EB)
      LOCAL_PACKING_RATIO = 0._EB
      LPC2_LOOP: DO N_LPC2=1,N_LAGRANGIAN_CLASSES
         LPC2 => LAGRANGIAN_PARTICLE_CLASS(N_LPC2)
         IF (AVG_DROP_RAD(IIG_OLD,JJG_OLD,KKG_OLD,LPC2%ARRAY_INDEX)<TWO_EPSILON_EB) CYCLE LPC2_LOOP
         SELECT CASE(SURFACE(LPC2%SURF_INDEX)%GEOMETRY)
            CASE(SURF_CARTESIAN)
               LPC_GEOM_FACTOR = 1._EB
            CASE(SURF_CYLINDRICAL)
               LPC_GEOM_FACTOR = 0.5_EB*PI
            CASE(SURF_SPHERICAL)
               LPC_GEOM_FACTOR = FOTH
         END SELECT
         LOCAL_PACKING_RATIO = LOCAL_PACKING_RATIO + AVG_DROP_AREA(IIG_OLD,JJG_OLD,KKG_OLD,LPC2%ARRAY_INDEX)*&
                  AVG_DROP_RAD(IIG_OLD,JJG_OLD,KKG_OLD,LPC2%ARRAY_INDEX)*LPC_GEOM_FACTOR
      ENDDO LPC2_LOOP
      LOCAL_PACKING_RATIO = LOCAL_PACKING_RATIO - EMBER_PACKING_RATIO
      LOCAL_PACKING_RATIO = MAX(MIN(1._EB, LOCAL_PACKING_RATIO),0._EB)
      CALL RANDOM_NUMBER(RN)
      IF (LOCAL_PACKING_RATIO>TWO_EPSILON_EB .AND. &
         RN<(LOCAL_PACKING_RATIO*EMBER_PACKING_RATIO)**LPC%EMBER_SNAG_FACTOR) THEN
         STUCK=.TRUE.
         LP%U = 0._EB
         LP%V = 0._EB
         LP%W = 0._EB
      ENDIF
   ENDIF
ENDIF

! Calculate the particle drag coefficient

RHO_G = RHO(IIG_OLD,JJG_OLD,KKG_OLD)
UREL  = LP%U - UBAR
VREL  = LP%V - VBAR
WREL  = LP%W - WBAR
IF (TWO_D) VREL=0._EB
QREL  = SQRT(UREL*UREL + VREL*VREL + WREL*WREL)
B1%U_TANG = QREL

DRAG_LAW_SELECT: SELECT CASE (LPC%DRAG_LAW)

   CASE (USER_DRAG)

      C_DRAG = LPC%DRAG_COEFFICIENT(1)

   CASE (SCREEN_DRAG, POROUS_DRAG)  ! Calculate this elsewhere

   CASE DEFAULT

      TMP_G  = MAX(TMPMIN,TMP(IIG_OLD,JJG_OLD,KKG_OLD))
      ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(IIG_OLD,JJG_OLD,KKG_OLD,1:N_TRACKED_SPECIES)
      CALL GET_VISCOSITY(ZZ_GET,MU_FILM,TMP_G)
      LP%RE  = RHO_G*QREL*2._EB*R_D/MU_FILM
      KN = 0._EB
      IF (LP%RE<1._EB) KN = MU_FILM*SQRT(0.5_EB*PI/(PBAR(KKG_OLD,PRESSURE_ZONE(IIG_OLD,JJG_OLD,KKG_OLD))*RHO_G))/(2._EB*R_D)
      ! Reynolds number based on hydraulic diameter for circular or square disk
      IF (LPC%DRAG_LAW==DISK_DRAG) LP%RE = RHO_G*QREL*SQRT(B1%AREA/2._EB)/MU_FILM
      C_DRAG = DRAG(LP%RE,LPC%DRAG_LAW,KN)

      ! Primary break-up model

      PRIMARY_BREAKUP_IF: IF ((T-LP%T_INSERT)<LPC%PRIMARY_BREAKUP_TIME) THEN

         C_DRAG = C_DRAG * LPC%PRIMARY_BREAKUP_DRAG_REDUCTION_FACTOR

      ELSE PRIMARY_BREAKUP_IF

         ! Drag reduction model for liquid droplets

         WAKE_VEL=1.0_EB
         IF (LPC%LIQUID_DROPLET) THEN
            DROP_VOL_FRAC = MIN(1._EB,AVG_DROP_DEN(IIG_OLD,JJG_OLD,KKG_OLD,LPC%ARRAY_INDEX)/LPC%DENSITY)
            IF (DROP_VOL_FRAC > LPC%DENSE_VOLUME_FRACTION) CALL WAKE_REDUCTION(DROP_VOL_FRAC,LP%RE,C_DRAG,WAKE_VEL)
         ENDIF

         ! Secondary break-up model

         SECONDARY_BREAKUP_IF: IF (LPC%BREAKUP) THEN
            ! Use undisturbed wake velocity for breakup calculations
            WAKE_VEL    = WAKE_VEL*QREL
            RE_WAKE     = RHO_G*WAKE_VEL   *2._EB*R_D/MU_FILM
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
                  R_D = ONTHHALF*R_D
               ELSE
                  DO WHILE (R_D >= R_D_0)
                     BREAKUP_RADIUS = LPC%BREAKUP_RATIO*R_D_0*LPC%BREAKUP_R_CNF(:)
                     CALL RANDOM_CHOICE(LPC%BREAKUP_CNF(:),BREAKUP_RADIUS,NDC,R_D)
                  END DO
                  R_D = MAX(R_D,1.1_EB*LPC%MINIMUM_DIAMETER/2._EB)
               ENDIF
               LP%RE    = RHO_G*QREL*2._EB*R_D/MU_FILM
               C_DRAG   = DRAG(LP%RE,LPC%DRAG_LAW)
               LP%PWT   = LP%PWT*(R_D_0/R_D)**3
               LP%T_INSERT = T
               LP%RADIUS = R_D
               LP%MASS = FOTHPI*LPC%DENSITY*R_D**3
               ! Redo wake reduction and shape deformation for the new drop
               ! Drag reduction, except for particles associated with a SURF line
               WAKE_VEL = 1.0_EB
               IF (LPC%LIQUID_DROPLET) THEN
                  DROP_VOL_FRAC = MIN(1._EB,AVG_DROP_DEN(IIG_OLD,JJG_OLD,KKG_OLD,LPC%ARRAY_INDEX)/LPC%DENSITY)
                  IF (DROP_VOL_FRAC > LPC%DENSE_VOLUME_FRACTION) CALL WAKE_REDUCTION(DROP_VOL_FRAC,LP%RE,C_DRAG,WAKE_VEL)
               ENDIF
               ! Change in drag coefficient due to deformation of PARTICLE shape (WE_G > 2)
               WAKE_VEL = WAKE_VEL*QREL
               RE_WAKE  = RHO_G*WAKE_VEL   *2._EB*R_D/MU_FILM
               WE_G     = RHO_G*WAKE_VEL**2*2._EB*R_D/LPC%SURFACE_TENSION
               ! Shape Deformation
               C_DRAG   = SHAPE_DEFORMATION(RE_WAKE,WE_G,C_DRAG)
            ENDIF AGE_IF
         ENDIF SECONDARY_BREAKUP_IF

      ENDIF PRIMARY_BREAKUP_IF

END SELECT DRAG_LAW_SELECT

! Calculate the cross-sectional area of the droplet or particle

IF (LPC%DRAG_LAW/=SCREEN_DRAG .AND. LPC%DRAG_LAW/=POROUS_DRAG) THEN
   IF (LPC%LIQUID_DROPLET) THEN
      A_DRAG = PI*R_D**2
   ELSE
      SELECT CASE(SF%GEOMETRY)
         ! Note: LPC%SHAPE_FACTOR=0.25 by default, which is analytically correct for SURF_SPHERICAL,
         !       and accounts for random orientations of other geometries
         CASE(SURF_CARTESIAN)
            A_DRAG = 2._EB*SF%LENGTH*SF%WIDTH*LPC%SHAPE_FACTOR
            ! For disk drag, allow for different area for each particle
            IF (LPC%DRAG_LAW==DISK_DRAG) A_DRAG = B1%AREA/2._EB
         CASE(SURF_CYLINDRICAL)
            A_DRAG = 2._EB*PI*R_D*SF%LENGTH*LPC%SHAPE_FACTOR
         CASE(SURF_SPHERICAL)
            A_DRAG = 4._EB*PI*R_D**2*LPC%SHAPE_FACTOR
      END SELECT
   ENDIF
ENDIF

! Experimental ember generation model

IF (LPC%EMBER_PARTICLE .AND. .NOT.LP%EMBER) THEN
   EMBER_DENSITY = 0._EB
   IF (EMBER_VOLUME>TWO_EPSILON_EB) EMBER_DENSITY = LP%MASS/EMBER_VOLUME
   IF ( EMBER_DENSITY < LPC%EMBER_DENSITY_THRESHOLD .AND. &
        QREL > LPC%EMBER_VELOCITY_THRESHOLD ) THEN
      IF (LPC%TRACK_EMBERS) THEN
         LP%EMBER = .TRUE.
      ELSE
         B1%BURNAWAY = .TRUE.
      ENDIF
   ENDIF
ENDIF

! Move the particles unless they are STATIC

PARTICLE_NON_STATIC_IF: IF (.NOT.LPC%STATIC .OR. (LP%EMBER .AND. .NOT.STUCK)) THEN ! Move airborne, non-stationary particles

   ! Compute gravity components

   IF (SPATIAL_GRAVITY_VARIATION) THEN
      GX_LOC = EVALUATE_RAMP(X_OLD,I_RAMP_GX)*GVEC(1)
      GY_LOC = EVALUATE_RAMP(X_OLD,I_RAMP_GY)*GVEC(2)
      GZ_LOC = EVALUATE_RAMP(X_OLD,I_RAMP_GZ)*GVEC(3)
   ELSEIF (I_RAMP_GX>0 .OR. I_RAMP_GY>0 .OR. I_RAMP_GZ>0) THEN
      GX_LOC = EVALUATE_RAMP(T,I_RAMP_GX)*GVEC(1)
      GY_LOC = EVALUATE_RAMP(T,I_RAMP_GY)*GVEC(2)
      GZ_LOC = EVALUATE_RAMP(T,I_RAMP_GZ)*GVEC(3)
   ELSE
      GX_LOC = GVEC(1)
      GY_LOC = GVEC(2)
      GZ_LOC = GVEC(3)
   ENDIF

   ! Update particle velocity

   BETA = 0.5_EB*RHO_G*C_DRAG*A_DRAG*QREL/LP%MASS

   IF (BETA>TWO_EPSILON_EB) THEN
      LP%U = UBAR + (UREL-GX_LOC/BETA)*EXP(-BETA*DT_P) + GX_LOC/BETA
      LP%V = VBAR + (VREL-GY_LOC/BETA)*EXP(-BETA*DT_P) + GY_LOC/BETA
      LP%W = WBAR + (WREL-GZ_LOC/BETA)*EXP(-BETA*DT_P) + GZ_LOC/BETA
   ELSE
      LP%U = U_OLD + DT_P*GX_LOC
      LP%V = V_OLD + DT_P*GY_LOC
      LP%W = W_OLD + DT_P*GZ_LOC
   ENDIF

   ! Fluid momentum source term

   IF (BETA>TWO_EPSILON_EB) THEN
      MPOM = LP%PWT*LP%RVC/RHO_G
      M_DOT = SUM(B1%M_DOT_G_PP_ACTUAL(1:N_TRACKED_SPECIES))*B1%AREA
      ACCEL_X = MPOM*(LP%MASS*((U_OLD-LP%U)/DT_P+GX_LOC) + M_DOT*UREL)
      ACCEL_Y = MPOM*(LP%MASS*((V_OLD-LP%V)/DT_P+GY_LOC) + M_DOT*VREL)
      ACCEL_Z = MPOM*(LP%MASS*((W_OLD-LP%W)/DT_P+GZ_LOC) + M_DOT*WREL)
   ELSE
      ACCEL_X  = 0._EB
      ACCEL_Y  = 0._EB
      ACCEL_Z  = 0._EB
   ENDIF

   ! Update particle position

   BC%X = X_OLD + 0.5_EB*DT_P*(LP%U+U_OLD)
   BC%Y = Y_OLD + 0.5_EB*DT_P*(LP%V+V_OLD)
   BC%Z = Z_OLD + 0.5_EB*DT_P*(LP%W+W_OLD)

   IF (AEROSOL_SCRUBBING) CALL DROPLET_SCRUBBING(IP,NM,DT,DT_P)

   IF (TWO_D) THEN
      BC%Y = Y_OLD
      LP%V = 0._EB
      ACCEL_Y  = 0._EB
   ENDIF

   IF (PARTICLE_CFL) PART_UVWMAX = MAX(PART_UVWMAX,MAX( ABS(LP%U)*RDX(IIG_OLD),ABS(LP%V)*RDY(JJG_OLD),ABS(LP%W)*RDZ(KKG_OLD)))

ELSE PARTICLE_NON_STATIC_IF ! Drag calculation for stationary, airborne particles

   SELECT CASE (LPC%DRAG_LAW)
      CASE DEFAULT
         M_DOT = SUM(B1%M_DOT_G_PP_ACTUAL(1:N_TRACKED_SPECIES))*B1%AREA
         BETA = LP%PWT*LP%RVC*(0.5_EB*C_DRAG*A_DRAG*QREL + M_DOT/RHO_G)
         ACCEL_X = -UBAR*BETA
         ACCEL_Y = -VBAR*BETA
         ACCEL_Z = -WBAR*BETA
      CASE (SCREEN_DRAG)
         IF (QREL > 0.015_EB .AND. LPC%FREE_AREA_FRACTION < 1.0_EB ) THEN ! Testing shows below this can have instability
            TMP_G  = MAX(TMPMIN,TMP(IIG_OLD,JJG_OLD,KKG_OLD))
            ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(IIG_OLD,JJG_OLD,KKG_OLD,1:N_TRACKED_SPECIES)
            CALL GET_VISCOSITY(ZZ_GET,MU_FILM,TMP_G)
            Y_TERM = LPC%DRAG_COEFFICIENT * RHO_G /SQRT(LPC%PERMEABILITY)*QREL*ABS(ORIENTATION_VECTOR(1:3,LPC%ORIENTATION_INDEX))
            K_TERM = MU_FILM/LPC%PERMEABILITY*ABS(ORIENTATION_VECTOR(1:3,LPC%ORIENTATION_INDEX))
            SFAC = 2._EB*LP%RADIUS*LP%RVC/RHO_G
            ACCEL_X = -(K_TERM(1)+Y_TERM(1))*UBAR*DY(JJG_OLD)*DZ(KKG_OLD)*SFAC
            ACCEL_Y = -(K_TERM(2)+Y_TERM(2))*VBAR*DX(IIG_OLD)*DZ(KKG_OLD)*SFAC
            ACCEL_Z = -(K_TERM(3)+Y_TERM(3))*WBAR*DX(IIG_OLD)*DY(JJG_OLD)*SFAC
         ELSE
            ACCEL_X = 0._EB
            ACCEL_Y = 0._EB
            ACCEL_Z = 0._EB
         ENDIF
      CASE (POROUS_DRAG)
         TMP_G  = MAX(TMPMIN,TMP(IIG_OLD,JJG_OLD,KKG_OLD))
         ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(IIG_OLD,JJG_OLD,KKG_OLD,1:N_TRACKED_SPECIES)
         CALL GET_VISCOSITY(ZZ_GET,MU_FILM,TMP_G)
         Y_TERM = LPC%DRAG_COEFFICIENT * RHO_G /SQRT(LPC%PERMEABILITY)*QREL
         K_TERM = MU_FILM/LPC%PERMEABILITY
         SFAC = 1._EB/RHO_G
         ACCEL_X = -(MIN(DX(IIG_OLD),LP%DX)/DX(IIG_OLD))*(K_TERM(1)+Y_TERM(1))*UBAR*SFAC
         ACCEL_Y = -(MIN(DY(JJG_OLD),LP%DY)/DY(JJG_OLD))*(K_TERM(2)+Y_TERM(2))*VBAR*SFAC
         ACCEL_Z = -(MIN(DZ(KKG_OLD),LP%DZ)/DZ(KKG_OLD))*(K_TERM(3)+Y_TERM(3))*WBAR*SFAC
   END SELECT
   DRAG_MAX(1) = ACCEL_X/(-UBAR-SIGN(1._EB,UBAR)*TWO_EPSILON_EB)
   DRAG_MAX(2) = ACCEL_Y/(-VBAR-SIGN(1._EB,VBAR)*TWO_EPSILON_EB)
   DRAG_MAX(3) = ACCEL_Z/(-WBAR-SIGN(1._EB,WBAR)*TWO_EPSILON_EB)
   IF (TWO_D) THEN
      ACCEL_Y  = 0._EB
      DRAG_MAX(2) = 0._EB
   ENDIF
   IF (ANY(ABS(DRAG_MAX)>PART_UVWMAX)) PART_UVWMAX = MAX(PART_UVWMAX,MAXVAL(DRAG_MAX))

ENDIF PARTICLE_NON_STATIC_IF

! Distribute particle drag onto the mesh

LP%ACCEL_X = LP%ACCEL_X + ACCEL_X
LP%ACCEL_Y = LP%ACCEL_Y + ACCEL_Y
LP%ACCEL_Z = LP%ACCEL_Z + ACCEL_Z

DO AXIS=IAXIS,KAXIS
   IL = IIX; JL = JJY; KL = KKZ
   IF (AXIS == IAXIS) THEN
      LP_FORCE = ACCEL_X/LP%RVC
      IL       = FLOOR(XI)
      FV_D     => FVX_D
   ELSEIF (AXIS == JAXIS) THEN
      LP_FORCE = ACCEL_Y/LP%RVC
      JL       = FLOOR(YJ)
      FV_D     => FVY_D
   ELSEIF (AXIS == KAXIS) THEN
      LP_FORCE = ACCEL_Z/LP%RVC
      KL       = FLOOR(ZK)
      FV_D     => FVZ_D
   ENDIF
   CALL GET_FACE_VOLUMES(AXIS,IL,JL,KL,FACE_VOLS)
   VOL_WGT = FACE_VOLS*WGT(:,:,:,AXIS)
   IF (ANY(VOL_WGT>TWO_EPSILON_EB)) VOL_WGT=WGT(:,:,:,AXIS)/SUM(VOL_WGT)
   FV_D(IL:IL+1, JL:JL+1, KL:KL+1) = FV_D(IL:IL+1, JL:JL+1, KL:KL+1) - LP_FORCE*VOL_WGT
ENDDO

! store C_DRAG for output

SELECT CASE(LPC%DRAG_LAW)
   CASE DEFAULT
      LP%C_DRAG = C_DRAG
   CASE (SCREEN_DRAG, POROUS_DRAG)
      LP%C_DRAG = MAXVAL(LPC%DRAG_COEFFICIENT)
END SELECT

END SUBROUTINE MOVE_IN_GAS


!> \brief Test to see if particles should be removed from solid WALL or CFACE
!> \details If a particle has hit a solid boundary (LP%WALL_INDEX>0 or LP%CFACE_INDEX>0) do the following:
!> If the user has specified AMPUA (Accumulated Mass Per Unit Area) for
!> a particle type, and the particle has not already been counted (LP%SPLAT=F), add its mass to an array.
!> If the solid surface has an inward normal velocity (U_NORMAL>0) and
!> this velocity is greater than a user-specified minimum (PARTICLE_EXTRACTION_VELOCITY) OR if a droplet has
!> hit the floor (BC%Z<ZS) and the floor is POROUS, remove it.

SUBROUTINE VENT_PARTICLE_EXTRACTION(HIT_SOLID,EXTRACT,WALL_INDEX,CFACE_INDEX)

USE GEOMETRY_FUNCTIONS, ONLY: SEARCH_OTHER_MESHES
LOGICAL, INTENT(OUT) :: EXTRACT
LOGICAL, INTENT(INOUT) :: HIT_SOLID
LOGICAL :: SET_EXTRACT
INTEGER, INTENT(IN), OPTIONAL :: WALL_INDEX,CFACE_INDEX
INTEGER :: SURF_INDEX,NOM,IIO,JJO,KKO
TYPE (BOUNDARY_PROP1_TYPE), POINTER :: B1
TYPE (BOUNDARY_PROP2_TYPE), POINTER :: B2
TYPE (BOUNDARY_COORD_TYPE), POINTER :: BCX
TYPE (WALL_TYPE), POINTER :: WC
TYPE (CFACE_TYPE), POINTER :: CFA

EXTRACT = .FALSE.
SET_EXTRACT = .FALSE.

IF (PRESENT(WALL_INDEX)) THEN
   WC => WALL(WALL_INDEX)
   B1 => BOUNDARY_PROP1(WC%B1_INDEX)
   B2 => BOUNDARY_PROP2(WC%B2_INDEX)
   BCX => BOUNDARY_COORD(WC%BC_INDEX)
   SURF_INDEX = WC%SURF_INDEX
ELSEIF (PRESENT(CFACE_INDEX)) THEN
   CFA => CFACE(CFACE_INDEX)
   B1 => BOUNDARY_PROP1(CFA%B1_INDEX)
   B2 => BOUNDARY_PROP2(CFA%B2_INDEX)
   BCX => BOUNDARY_COORD(CFA%BC_INDEX)
   SURF_INDEX = CFA%SURF_INDEX
ENDIF

IF (ACCUMULATE_WATER .AND. .NOT.LP%SPLAT .AND. LPC%ADHERE_TO_SOLID) THEN
   B2%A_LP_MPUA(LPC%ARRAY_INDEX) = B2%A_LP_MPUA(LPC%ARRAY_INDEX) + LP%PWT*LPC%FTPR*R_D**3/B1%AREA
   LP%SPLAT = .TRUE.
ENDIF

IF_EXTRACT: IF ( B1%U_NORMAL>SURFACE(SURF_INDEX)%PARTICLE_EXTRACTION_VELOCITY .OR. &
     (POROUS_FLOOR .AND. BC%Z<ZS .AND. LPC%LIQUID_DROPLET) ) THEN
   CALL SEARCH_OTHER_MESHES(BC%X,BC%Y,BC%Z,NOM,IIO,JJO,KKO)
   IF (POROUS_FLOOR .AND. BC%Z<ZS .AND. LPC%LIQUID_DROPLET .AND. NOM>0) EXIT IF_EXTRACT  ! Do not remove particle still in domain
   IF (B1%NODE_INDEX > 0) THEN
      IF (DUCTNODE(B1%NODE_INDEX)%TRANSPORT_PARTICLES) THEN
         SELECT CASE (BCX%IOR)
            CASE(1)
               IF(BCX%IIG-2 < 1) THEN
                  SET_EXTRACT = .TRUE.
               ELSE
                  IF(CELL(CELL_INDEX(BCX%IIG-2,BCX%JJG,BCX%KKG))%SOLID) THEN
                     SET_EXTRACT = .TRUE.
                  ELSE
                     BC%X = X(BCX%IIG-2)-0.01_EB*DX(BCX%IIG-2)
                     HIT_SOLID = .FALSE.
                  ENDIF
               ENDIF
            CASE(-1)
               IF(BCX%IIG+2 > IBAR) THEN
                  SET_EXTRACT = .TRUE.
               ELSE
                  IF(CELL(CELL_INDEX(BCX%IIG+2,BCX%JJG,BCX%KKG))%SOLID) THEN
                     SET_EXTRACT = .TRUE.
                  ELSE
                     BC%X = X(BCX%IIG+1)+0.01_EB*DX(BCX%IIG+2)
                     HIT_SOLID = .FALSE.
                  ENDIF
               ENDIF
            CASE(2)
               IF(BCX%JJG-2 < 1) THEN
                  SET_EXTRACT = .TRUE.
               ELSE
                  IF(CELL(CELL_INDEX(BCX%IIG,BCX%JJG-2,BCX%KKG))%SOLID) THEN
                     SET_EXTRACT = .TRUE.
                  ELSE
                     BC%Y = Y(BCX%JJG-2)-0.01_EB*DY(BCX%JJG-2)
                     HIT_SOLID = .FALSE.
                  ENDIF
               ENDIF
            CASE(-2)
               IF(BCX%JJG+2 > JBAR) THEN
                  SET_EXTRACT = .TRUE.
               ELSE
                  IF(CELL(CELL_INDEX(BCX%IIG,BCX%JJG+2,BCX%KKG))%SOLID) THEN
                     SET_EXTRACT = .TRUE.
                  ELSE
                     BC%Y = Y(BCX%JJG+1)+0.01_EB*DY(BCX%JJG+2)
                     HIT_SOLID = .FALSE.
                  ENDIF
               ENDIF
            CASE(3)
               IF(BCX%KKG-2 < 1) THEN
                  SET_EXTRACT = .TRUE.
               ELSE
                  IF(CELL(CELL_INDEX(BCX%IIG,BCX%JJG,BCX%KKG-2))%SOLID) THEN
                     SET_EXTRACT = .TRUE.
                  ELSE
                     BC%Z = Z(BCX%KKG-2)-0.01_EB*DZ(BCX%KKG-2)
                     HIT_SOLID = .FALSE.
                  ENDIF
               ENDIF
            CASE(-3)
               IF(BCX%KKG+2 > KBAR) THEN
                  SET_EXTRACT = .TRUE.
               ELSE
                  IF(CELL(CELL_INDEX(BCX%IIG,BCX%JJG,BCX%KKG+2))%SOLID) THEN
                     SET_EXTRACT = .TRUE.
                  ELSE
                     BC%Z = Z(BCX%KKG+1)+0.01_EB*DZ(BCX%KKG+2)
                     HIT_SOLID = .FALSE.
                  ENDIF
               ENDIF
         END SELECT
      ELSE
         SET_EXTRACT = .TRUE.
      ENDIF
   ELSE
      SET_EXTRACT = .TRUE.
   ENDIF
   IF (SET_EXTRACT) THEN
      BC%X=-1.E6_EB
      EXTRACT = .TRUE.
   ENDIF
ENDIF IF_EXTRACT

END SUBROUTINE VENT_PARTICLE_EXTRACTION


!> \brief Compute C_DRAG reduction due to the wake effect (Ramirez, Munoz et al. 2007)

SUBROUTINE WAKE_REDUCTION(DROP_VOL_FRAC,RE,C_DRAG,WAKE_VEL)

REAL(EB), INTENT(INOUT) :: C_DRAG
REAL(EB), INTENT(IN)  :: DROP_VOL_FRAC,RE
REAL(EB), INTENT(OUT) :: WAKE_VEL
REAL(EB) :: LODM,RELOD

LODM     = MAX(TWO_EPSILON_EB,(PI/(6._EB*DROP_VOL_FRAC))**ONTH-0.5_EB)
RELOD    = RE/(16._EB * LODM)
WAKE_VEL = 1._EB - 0.5_EB*C_DRAG*(1._EB - EXP(-RELOD))
WAKE_VEL = MAX(WAKE_VEL,0.15_EB)
C_DRAG   = C_DRAG * WAKE_VEL * (1._EB + (RELOD/LODM)*EXP(-RELOD))

END SUBROUTINE WAKE_REDUCTION


!> \brief Compute shape and drag of deformable droplets
!>
!> \detail E.Loth, Quasi-steady shape and drag of deformable bubbles and drops, International Journal of Multiphase Flow 34 (2008)

REAL(EB) FUNCTION SHAPE_DEFORMATION(RE,WE,C_DRAG)

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
    SHAPE_DEFORMATION=C_DRAGNEW*E**(-TWTH)
ELSE
    SHAPE_DEFORMATION=C_DRAG
ENDIF

END FUNCTION SHAPE_DEFORMATION

!> \brief Get Inverse-Distance Weight (IDW) values for nearest gas faces
!> \param AXIS The axis of the face quantity
!> \param I The lower x index
!> \param J The lower y index
!> \param K The lower z index
!> \param P_X Sample point location in x
!> \param P_Y Sample point location in y
!> \param P_Z Sample point location in z
!> \param IWD 2x2x2 matrix of weights for cartesian faces

SUBROUTINE GET_FACE_IDW(AXIS,I,J,K,P_X,P_Y,P_Z,IDW)

REAL(EB), INTENT(IN) :: P_X,P_Y,P_Z
REAL(EB), INTENT(OUT) :: IDW(2,2,2)
INTEGER, INTENT(IN) :: AXIS,I,J,K
INTEGER :: II,JJ,KK,ICF
REAL(EB) :: DIST,D_WGT,XYZ_INT(3)
REAL(EB), POINTER :: X_F(:),Y_F(:),Z_F(:)

IDW=0._EB

X_F=>XC
Y_F=>YC
Z_F=>ZC
SELECT CASE(AXIS)
   CASE(IAXIS); X_F=>X
   CASE(JAXIS); Y_F=>Y
   CASE(KAXIS); Z_F=>Z
END SELECT

! Get weights from distance to surrounding faces
FACE_LOOP: DO KK=K,K+1
   DO JJ=J,J+1
      DO II=I,I+1
         DIST = 0._EB
         ICF = 0
         IF (CC_IBM) ICF = FCVAR(II,JJ,KK,CC_IDCF,AXIS)
         IF (ICF>0) THEN
            XYZ_INT = CUT_FACE(ICF)%XYZCEN(IAXIS:KAXIS,1)
         ELSE
            XYZ_INT = (/X_F(II),Y_F(JJ),Z_F(KK)/)
         ENDIF
         DIST = NORM2((/P_X,P_Y,P_Z/)-XYZ_INT)
         ! Special case where location is directly on face
         IF (DIST<TWO_EPSILON_EB) THEN
            IDW = 0._EB
            IDW(II-I+1,JJ-J+1,KK-K+1) = 1._EB
            EXIT FACE_LOOP
         ELSE
            D_WGT = 1._EB/DIST**6._EB
         ENDIF
         ! face is solid
         IF(CELL(CELL_INDEX(II,JJ,KK))%WALL_INDEX(AXIS)>0) D_WGT = 0._EB
         IF (CC_IBM .AND. D_WGT>0._EB) THEN
            IF(FCVAR(II,JJ,KK,CC_FGSC,AXIS)==CC_SOLID) D_WGT = 0._EB
         ENDIF
         IDW(II-I+1,JJ-J+1,KK-K+1) = D_WGT
      ENDDO
   ENDDO
ENDDO FACE_LOOP

! Solid face between particle and vector
IF (I<IIG_OLD .AND. CELL(CELL_INDEX(IIG_OLD,JJG_OLD,KKG_OLD))%WALL_INDEX(-1)>0) IDW(1,:,:) = 0._EB
IF (I==IIG_OLD .AND. CELL(CELL_INDEX(IIG_OLD,JJG_OLD,KKG_OLD))%WALL_INDEX(1)>0) IDW(2,:,:) = 0._EB
IF (J<JJG_OLD .AND. CELL(CELL_INDEX(IIG_OLD,JJG_OLD,KKG_OLD))%WALL_INDEX(-2)>0) IDW(:,1,:) = 0._EB
IF (J==JJG_OLD .AND. CELL(CELL_INDEX(IIG_OLD,JJG_OLD,KKG_OLD))%WALL_INDEX(2)>0) IDW(:,2,:) = 0._EB
IF (K<KKG_OLD .AND. CELL(CELL_INDEX(IIG_OLD,JJG_OLD,KKG_OLD))%WALL_INDEX(-3)>0) IDW(:,:,1) = 0._EB
IF (K==KKG_OLD .AND. CELL(CELL_INDEX(IIG_OLD,JJG_OLD,KKG_OLD))%WALL_INDEX(3)>0) IDW(:,:,2) = 0._EB

IF (CC_IBM) THEN
   IF (I<IIG_OLD .AND. FCVAR(IIG_OLD-1,JJG_OLD,KKG_OLD,CC_FGSC,IAXIS)==CC_SOLID) IDW(1,:,:) = 0._EB
   IF (I==IIG_OLD .AND. FCVAR(IIG_OLD,JJG_OLD,KKG_OLD,CC_FGSC,IAXIS)==CC_SOLID)  IDW(2,:,:) = 0._EB
   IF (J<JJG_OLD .AND. FCVAR(IIG_OLD,JJG_OLD-1,KKG_OLD,CC_FGSC,JAXIS)==CC_SOLID) IDW(:,1,:) = 0._EB
   IF (J==JJG_OLD .AND. FCVAR(IIG_OLD,JJG_OLD,KKG_OLD,CC_FGSC,JAXIS)==CC_SOLID)  IDW(:,2,:) = 0._EB
   IF (K<KKG_OLD .AND. FCVAR(IIG_OLD,JJG_OLD,KKG_OLD-1,CC_FGSC,KAXIS)==CC_SOLID) IDW(:,:,1) = 0._EB
   IF (K==KKG_OLD .AND. FCVAR(IIG_OLD,JJG_OLD,KKG_OLD,CC_FGSC,KAXIS)==CC_SOLID)  IDW(:,:,2) = 0._EB
ENDIF

! Normalize
IF (ANY(IDW>TWO_EPSILON_EB)) IDW = IDW/SUM(IDW)

END SUBROUTINE GET_FACE_IDW

!> \brief Get Tri-Linear interpolation Weight (TLW) values for nearest gas faces
!> \param AXIS The axis of the face quantity
!> \param I The lower x index
!> \param J The lower y index
!> \param K The lower z index
!> \param P_X Sample point location in x
!> \param P_Y Sample point location in y
!> \param P_Z Sample point location in z
!> \param TLW 2x2x2 matrix of weights for cartesian faces
!> \param V 2x2x2 array of velocity box corner values

SUBROUTINE GET_FACE_TLW(AXIS,I,J,K,P_X,P_Y,P_Z,TLW,V)

REAL(EB), INTENT(IN) :: P_X,P_Y,P_Z
REAL(EB), INTENT(OUT) :: TLW(0:1,0:1,0:1)
INTEGER, INTENT(IN) :: AXIS,I,J,K
REAL(EB), INTENT(INOUT) :: V(2,2,2)
REAL(EB) :: P,PP,R,RR,S,SS
INTEGER :: IWC(-3:3)
REAL(EB), POINTER :: X_F(:),Y_F(:),Z_F(:)

TLW=0._EB

X_F=>XC
Y_F=>YC
Z_F=>ZC
SELECT CASE(AXIS)
   CASE(IAXIS); X_F=>X
   CASE(JAXIS); Y_F=>Y
   CASE(KAXIS); Z_F=>Z
END SELECT

P = (P_X-X_F(I))/(X_F(I+1)-X_F(I))
R = (P_Y-Y_F(J))/(Y_F(J+1)-Y_F(J))
S = (P_Z-Z_F(K))/(Z_F(K+1)-Z_F(K))

IWC = CELL(IC_OLD)%WALL_INDEX
IF (NEAR_WALL_PARTICLE_INTERPOLATION) THEN

   IF (AXIS/=IAXIS .AND. IIG_OLD> I .AND. IWC(-1)>0) THEN
      IF (WALL(IWC(-1))%BOUNDARY_TYPE==SOLID_BOUNDARY) THEN
         P=(P_X-BOUNDARY_COORD(WALL(IWC(-1))%BC_INDEX)%X)/(X_F(I+1)-BOUNDARY_COORD(WALL(IWC(-1))%BC_INDEX)%X)
         SELECT CASE(AXIS)
            CASE(JAXIS); V(1,:,:)=SURFACE(WALL(IWC(-1))%SURF_INDEX)%VEL_T(1)
            CASE(KAXIS); V(1,:,:)=SURFACE(WALL(IWC(-1))%SURF_INDEX)%VEL_T(2)
         END SELECT
      ENDIF
   ENDIF
   IF (AXIS/=IAXIS .AND. IIG_OLD==I .AND. IWC( 1)>0) THEN
      IF (WALL(IWC( 1))%BOUNDARY_TYPE==SOLID_BOUNDARY) THEN
         P=(P_X-X_F(I))/(BOUNDARY_COORD(WALL(IWC( 1))%BC_INDEX)%X-X_F(I))
         SELECT CASE(AXIS)
            CASE(JAXIS); V(2,:,:)=SURFACE(WALL(IWC( 1))%SURF_INDEX)%VEL_T(1)
            CASE(KAXIS); V(2,:,:)=SURFACE(WALL(IWC( 1))%SURF_INDEX)%VEL_T(2)
         END SELECT
      ENDIF
   ENDIF

   IF (AXIS/=JAXIS .AND. JJG_OLD> J .AND. IWC(-2)>0) THEN
      IF (WALL(IWC(-2))%BOUNDARY_TYPE==SOLID_BOUNDARY) THEN
         R=(P_Y-BOUNDARY_COORD(WALL(IWC(-2))%BC_INDEX)%Y)/(Y_F(J+1)-BOUNDARY_COORD(WALL(IWC(-2))%BC_INDEX)%Y)
         SELECT CASE(AXIS)
            CASE(IAXIS); V(:,1,:)=SURFACE(WALL(IWC(-2))%SURF_INDEX)%VEL_T(1)
            CASE(KAXIS); V(:,1,:)=SURFACE(WALL(IWC(-2))%SURF_INDEX)%VEL_T(2)
         END SELECT
      ENDIF
   ENDIF
   IF (AXIS/=JAXIS .AND. JJG_OLD==J .AND. IWC( 2)>0) THEN
      IF (WALL(IWC( 2))%BOUNDARY_TYPE==SOLID_BOUNDARY) THEN
         R=(P_Y-Y_F(J))/(BOUNDARY_COORD(WALL(IWC( 2))%BC_INDEX)%Y-Y_F(J))
         SELECT CASE(AXIS)
            CASE(IAXIS); V(:,2,:)=SURFACE(WALL(IWC( 2))%SURF_INDEX)%VEL_T(1)
            CASE(KAXIS); V(:,2,:)=SURFACE(WALL(IWC( 2))%SURF_INDEX)%VEL_T(2)
         END SELECT
      ENDIF
   ENDIF

   IF (AXIS/=KAXIS .AND. KKG_OLD> K .AND. IWC(-3)>0) THEN
      IF (WALL(IWC(-3))%BOUNDARY_TYPE==SOLID_BOUNDARY) THEN
         S=(P_Z-BOUNDARY_COORD(WALL(IWC(-3))%BC_INDEX)%Z)/(Z_F(K+1)-BOUNDARY_COORD(WALL(IWC(-3))%BC_INDEX)%Z)
         SELECT CASE(AXIS)
            CASE(IAXIS); V(:,:,1)=SURFACE(WALL(IWC(-3))%SURF_INDEX)%VEL_T(1)
            CASE(JAXIS); V(:,:,1)=SURFACE(WALL(IWC(-3))%SURF_INDEX)%VEL_T(2)
         END SELECT
      ENDIF
   ENDIF
   IF (AXIS/=KAXIS .AND. KKG_OLD==K .AND. IWC( 3)>0) THEN
      IF (WALL(IWC( 3))%BOUNDARY_TYPE==SOLID_BOUNDARY) THEN
         S=(P_Z-Z_F(K))/(BOUNDARY_COORD(WALL(IWC( 3))%BC_INDEX)%Z-Z_F(K))
         SELECT CASE(AXIS)
            CASE(IAXIS); V(:,:,2)=SURFACE(WALL(IWC( 3))%SURF_INDEX)%VEL_T(1)
            CASE(JAXIS); V(:,:,2)=SURFACE(WALL(IWC( 3))%SURF_INDEX)%VEL_T(2)
         END SELECT
      ENDIF
   ENDIF

ELSE
   IF (AXIS/=IAXIS .AND. IIG_OLD> I .AND. IWC(-1)>0) P = 1._EB
   IF (AXIS/=IAXIS .AND. IIG_OLD==I .AND. IWC( 1)>0) P = 0._EB
   IF (AXIS/=JAXIS .AND. JJG_OLD> J .AND. IWC(-2)>0) R = 1._EB
   IF (AXIS/=JAXIS .AND. JJG_OLD==J .AND. IWC( 2)>0) R = 0._EB
   IF (AXIS/=KAXIS .AND. KKG_OLD> K .AND. IWC(-3)>0) S = 1._EB
   IF (AXIS/=KAXIS .AND. KKG_OLD==K .AND. IWC( 3)>0) S = 0._EB
ENDIF

PP = 1._EB-P
RR = 1._EB-R
SS = 1._EB-S

TLW(0,0,0) = PP * RR * SS
TLW(1,0,0) = P  * RR * SS
TLW(0,1,0) = PP * R  * SS
TLW(0,0,1) = PP * RR * S
TLW(1,0,1) = P  * RR * S
TLW(0,1,1) = PP * R  * S
TLW(1,1,0) = P  * R  * SS
TLW(1,1,1) = P  * R  * S

END SUBROUTINE GET_FACE_TLW

!> \brief Return face volumes for distribution of quantity onto faces
!> \param AXIS The axis for the face centered quantity
!> \param I The lower x index
!> \param J The lower y index
!> \param K The lower z index
!> \param FACE_VOLS The face volumes

SUBROUTINE GET_FACE_VOLUMES(AXIS,I,J,K,FACE_VOLS)

INTEGER, INTENT(IN) :: AXIS,I,J,K
INTEGER :: II,JJ,KK,ICF
REAL(EB), INTENT(OUT) :: FACE_VOLS(2,2,2)
REAL(EB), POINTER, DIMENSION(:) :: DX1,DX2,DX3

DX1 => DX
DX2 => DY
DX3 => DZ
SELECT CASE (AXIS)
   CASE(IAXIS); DX1 => DXN
   CASE(JAXIS); DX2 => DYN
   CASE(KAXIS); DX3 => DZN
END SELECT

DO KK=K,K+1
   DO JJ=J,J+1
      DO II=I,I+1
         FACE_VOLS(II-I+1,JJ-J+1,KK-K+1) = DX1(II)*DX2(JJ)*DX3(KK)
         ICF = 0
         IF (CC_IBM) ICF = FCVAR(II,JJ,KK,CC_IDCF,AXIS)
         IF(ICF>0) FACE_VOLS(II-I+1,JJ-J+1,KK-K+1) = FACE_VOLS(II-I+1,JJ-J+1,KK-K+1)*CUT_FACE(ICF)%ALPHA_CF
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE GET_FACE_VOLUMES

END SUBROUTINE MOVE_PARTICLES


!> \brief Compute mass and energy transfer between gas and liquid PARTICLEs
!> \param T Current simulation time (s)
!> \param DT Time step size (s)
!> \param NM Current mesh index

SUBROUTINE PARTICLE_MASS_ENERGY_TRANSFER(T,DT,NM)

USE PHYSICAL_FUNCTIONS, ONLY : GET_FILM_PROPERTIES, GET_MASS_FRACTION,GET_AVERAGE_SPECIFIC_HEAT,&
                               GET_MOLECULAR_WEIGHT,GET_SPECIFIC_HEAT,GET_MASS_FRACTION_ALL,GET_SENSIBLE_ENTHALPY,&
                               GET_MW_RATIO, GET_EQUIL_DATA,DROPLET_H_MASS_H_HEAT_GAS,GET_ENTHALPY,GET_TEMPERATURE
USE MATH_FUNCTIONS, ONLY: INTERPOLATE1D_UNIFORM,F_B
USE COMP_FUNCTIONS, ONLY: SHUTDOWN
USE OUTPUT_DATA, ONLY: M_DOT,Q_DOT
USE TURBULENCE, ONLY: FORCED_CONVECTION_MODEL

REAL(EB), INTENT(IN) :: T,DT
INTEGER, INTENT(IN) :: NM

REAL(EB), POINTER, DIMENSION(:,:,:) :: MVAP_TOT
!< Amount of mass evaporated into a grid cell (kg)
REAL(EB), POINTER, DIMENSION(:,:,:) :: RHO_INTERIM
!< Current gas density (kg/m3)
REAL(EB), POINTER, DIMENSION(:,:,:) :: TMP_INTERIM
!< Current gas temperature (K)
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZ_INTERIM
!< Current gas species mass fractions
REAL(EB) :: C_DROP2 !< Specific heat of particle (J/kg/K)
REAL(EB) :: CP !< Specific heat (J/kg/K)
REAL(EB) :: CP_BAR !< Average specific heat (J/kg/K)
REAL(EB) :: CP_FILM !< Specific heat of the film (J/kg/K) at the film temperature
REAL(EB) :: D_FILM !< Diffusivity into air of the droplet species (m2/s) at the film temperature
REAL(EB) :: DELTA_H_G !< H_S_B - H_S (J)
REAL(EB) :: DH_V_A_DT !< Temperature derivative of H_VA (J/kg/K)
REAL(EB) :: H1 !< Sensible enthalpy (J/kg/K)
REAL(EB) :: H2 !< Sensible enthalpy (J/kg/K)
REAL(EB) :: H_D_OLD !< Particle enthalpy (J) at the start of a sub time step
REAL(EB) :: H_G_OLD !< Gas enthalpy (J) at the start of a sub time step
REAL(EB) :: H_HEAT !< Convection heat transfer coefficient between the particle and the gas (W/m2/K)
REAL(EB) :: H_L !< Enthalpy of the particle (J)
REAL(EB) :: H_MASS !< Mass transfer coefficient (m/s)
REAL(EB) :: H_NEW  !< Total gas and particle enthalpy (J) at the end of a sub time step
REAL(EB) :: H_S !< Sensible enthalpy of the vapor (J) at the gas temperature
REAL(EB) :: H_S_B !< Sensible enthalpy of the vapor (J) at the particle temperature
REAL(EB) :: H_S_G_OLD !< Sensible enthalpy of the gas (J) at the start of a sub time step
REAL(EB) :: H_V !< Heat of vaporization at the particle temperature (J/kg)
REAL(EB) :: H_V2 !< Heat of vaporization at the particle temperature (J/kg)
REAL(EB) :: H_V_A !< Effective heat of vaporization for use in the Clasius-Clapeyron relation (J/kg)
REAL(EB) :: H_V_REF !< Heat of vaporization at the species heat of vaporization reference temperature (J/kg)
REAL(EB) :: H_WALL !< Convection heat transfer coefficient between the particle and a surface (W/m2/K)
REAL(EB) :: K_FILM !< Conducitvity of the film (W/m/K) at the film temperature
REAL(EB) :: M_GAS !< Mass of the gas (kg) in the grid cell
REAL(EB) :: M_GAS_NEW !< Mass of gas (kg) in the grid cell at the end of a sub time step
REAL(EB) :: M_GAS_OLD !< Mass of gas (kg) in a grid cell at the start of a sub time step
REAL(EB) :: MU_FILM !< Viscosity of the film (kg/m/s) at the film temperature
REAL(EB) :: MW_GAS !< Molecular weight (kg/kmol) of the gas
REAL(EB) :: MW_RATIO !< Ratio of average gas molecular weigth to particle species molecular weight
REAL(EB) :: NU_FAC_GAS !< Nusselt number of a particle in the gas
REAL(EB) :: NUSSELT !< Nusselt number
REAL(EB) :: NUSSELT_LIQUID !< Nusselt number for liquid/solid convective heat transfer
REAL(EB) :: VEL_LIQUID !< Droplet/liquid velocity to be used in convective heat transfer coefficient
REAL(EB) :: RE_LIQUID !< Reynolds number used in convective heat transfer coefficient
REAL(EB) :: PR_FILM !< Prandtl number of the film
REAL(EB) :: Q_CON_GAS !< Convective heat transfer between the particle and the gas (J)
REAL(EB) :: Q_CON_SUM !< Sum of convective heat transfer between the particle and the surface over subtimesteps (J)
REAL(EB) :: Q_CON_WALL !< Convective heat transfer between the particle and a surface (J)
REAL(EB) :: Q_DOT_RAD ! Radiant heat transfer rate to particle (J/s)
REAL(EB) :: Q_FRAC !< Heat transfer adjustment factor when particle reaches boiling temperature during a sub time step
REAL(EB) :: Q_RAD !< Net radiation heat transfer to the particle (J)
REAL(EB) :: Q_RAD_SUM !< Sum of radiant heat transfer to the particle over subtimesteps (J)
REAL(EB) :: Q_TOT !< Total heat transfer from convection and radiation to the particle (J)
REAL(EB) :: RHO_FILM !< Density of the film (kg/m3) at the film temperature
REAL(EB) :: RHO_G !< Current gas density (kg/m3)
REAL(EB) :: LENGTH !< Length scale used in computing SH and NU
REAL(EB) :: RE_L !< Particle Reynolds number
REAL(EB) :: SC_FILM !< Schmidt number of the film
REAL(EB) :: SHERWOOD !< Particle Sherwood number
REAL(EB) :: SH_FAC_GAS !< Sherwood number of a particle in the gas
REAL(EB) :: M_VAP !< Mass evaporated (kg) from the particle in the current sub time step
REAL(EB) :: M_VAP_MAX !< Maximum allowable evaporation (kg)
REAL(EB) :: Y_ALL(1:N_SPECIES) !< Mass fraction of all primitive species
REAL(EB) :: X_DROP !< Equilibirum vapor mole fraction at the particle temperature
REAL(EB) :: Y_DROP !< Equilibrium vapor mass fraction at the particle temperature
REAL(EB) :: Y_DROP_A(1) !< Equilibrium vapor mass fraction, array needed for GET_FILM_PROPERTIES
REAL(EB) :: Y_COND !< Fraction of mass associated with any condensed vapor of the particle species
REAL(EB) :: Y_GAS !< Vapor fraction of the particle species
REAL(EB) :: Y_GAS_A(1) !< Vapor fraction of the particle species, array needed for GET_FILM_PROPERTIES
REAL(EB) :: Y_GAS_NEW !< End of sub time step vapor fraction of the particle species
REAL(EB) :: X_EQUIL !< Equilibrium vapor mole fraction
REAL(EB) :: Y_EQUIL !< Equilibrium vapor mass fraction
REAL(EB) :: U2 !< Relative u-velocity (m/s)
REAL(EB) :: V2 !< Relative v-velocity (m/s)
REAL(EB) :: W2 !< Relative w-velocity (m/s)
REAL(EB) :: VEL !< Relative velocity (m/s)
REAL(EB) :: WGT !< LAGRANGIAN_PARTICLE%PWT
REAL(EB) :: A_DROP !< Particle surface area (m2)
REAL(EB) :: C_DROP !< Specific heat of the particle (J/kg/K) at the current temperature
REAL(EB) :: DHOR !< Heat of vaporization divided by the gas constant (K).
REAL(EB) :: FTPR !< 4/3 * PI * Particle density
REAL(EB) :: M_DROP !< Mass of the  particle (kg)
REAL(EB) :: R_DROP !< Particle radius (m)
REAL(EB) :: MW_DROP !< Particle molecular weight (kg/kmol)
REAL(EB) :: T_BOIL_EFF !< Effective boiling temperature (K) of the particle at the current pressuer
REAL(EB) :: TMP_G_OLD !< Gas temperature (K) at the start of sub time step
REAL(EB) :: TMP_G !< Current gas temperature (K)
REAL(EB) :: TMP_G_I !< Gas temperature during search loop (K)
REAL(EB) :: TMP_G_NEW !< Gas temperature (K) at the end of a sub time step
REAL(EB) :: TMP_DROP !< Current particle temperature (K)
REAL(EB) :: TMP_DROP_NEW !< End of sub time step particle temperature (K)
REAL(EB) :: TMP_FILM !< Film temperature (K)
REAL(EB) :: TMP_MELT !< Melting temperaturte of the particle species (K)
REAL(EB) :: TMP_WALL !< Wall surface temperature (K)
REAL(EB) :: TMP_WALL_NEW !< New surface temperature (K)
REAL(EB) :: DT_SUBSTEP !< Current sub time step size (s)
REAL(EB) :: DT_SUM !< Running sum of completed sub time steps (s)
REAL(EB) :: ZZ_GET(1:N_TRACKED_SPECIES)
REAL(EB) :: ZZ_GET2(1:N_TRACKED_SPECIES)
REAL(EB) :: RHOCBAR ! Density of solid surface times specific heat the solid surface (J/m3/K)
REAL(EB) :: MCBAR !< Particle mass time particle specific heat (J/K)
REAL(EB) :: AGHRHO !< Collection of terms used in contstructing A_COL, B_COL, C_COL, and D_COL
REAL(EB) :: DTGOG !< Collection of terms used in contstructing A_COL, B_COL, C_COL, and D_COL
REAL(EB) :: DTGOP !< Collection of terms used in contstructing A_COL, B_COL, C_COL, and D_COL
REAL(EB) :: DTWOW !< Collection of terms used in contstructing A_COL, B_COL, C_COL, and D_COL
REAL(EB) :: DTWOP !< Collection of terms used in contstructing A_COL, B_COL, C_COL, and D_COL
REAL(EB) :: DTOG !< Collection of terms used in contstructing A_COL, B_COL, C_COL, and D_COL
REAL(EB) :: DTOP !< Collection of terms used in contstructing A_COL, B_COL, C_COL, and D_COL
REAL(EB) :: DAHVHLDY !< Collection of terms used in contstructing A_COL, B_COL, C_COL, and D_COL
REAL(EB) :: DADYHV !< Collection of terms used in contstructing A_COL, B_COL, C_COL, and D_COL
REAL(EB) :: DADYDTHVHL !< Collection of terms used in contstructing A_COL, B_COL, C_COL, and D_COL
REAL(EB) :: DADYDTHV !< Collection of terms used in contstructing A_COL, B_COL, C_COL, and D_COL
REAL(EB) :: DYDT !< Temperature derivative of the equilbirum vapor fraciton (1/K)
REAL(EB) :: A_COL(3) !< Gas temperature terms in LHS of solution
REAL(EB) :: B_COL(3) !< Particle temperatre terms in LHS of solution
REAL(EB) :: C_COL(3) !< Wall temperature terms in LHS of solution
REAL(EB) :: D_VEC(3) !< RHS of solution
INTEGER :: IP,II,JJ,KK,IW,ICF,N_LPC,ITMP,ITMP2,Y_INDEX,Z_INDEX,Z_INDEX_A(1),I_BOIL,I_MELT,NMAT
INTEGER :: ARRAY_CASE !< Heat transfer conditions 1=gas only, 2=isothermal wall present, 3=1D wall present
INTEGER(1), ALLOCATABLE, DIMENSION(:) :: PART_WARNING !< Tracks WARNING messages
!< 1 = Particle in gas only, 2 = Particle on constant temperature surface, 3 = Particle on thermally thick surface
CHARACTER(MESSAGE_LENGTH) :: MESSAGE
LOGICAL :: SF_FIXED,SURFACE_DROPLETS
TYPE (LAGRANGIAN_PARTICLE_TYPE), POINTER :: LP
TYPE (LAGRANGIAN_PARTICLE_CLASS_TYPE), POINTER :: LPC
TYPE(BOUNDARY_ONE_D_TYPE), POINTER :: ONE_D
TYPE(BOUNDARY_PROP1_TYPE), POINTER :: LP_B1,B1
TYPE(BOUNDARY_PROP2_TYPE), POINTER :: B2
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC
TYPE (SURFACE_TYPE), POINTER :: SF
TYPE (SPECIES_TYPE), POINTER :: SS
TYPE (WALL_TYPE), POINTER :: WC
TYPE (CFACE_TYPE), POINTER :: CFA
REAL(EB) :: TNOW

! Set the CPU timer and point to the current mesh variables

TNOW=CURRENT_TIME()
CALL POINT_TO_MESH(NM)

IF (NLP > 0) PARTICLES_EXISTED = .TRUE.

IF (.NOT. PARTICLES_EXISTED .OR. ((.NOT. LIQUID_DROPLETS) .AND. (.NOT. SOLID_PARTICLES))) THEN
   ! No particles ever existed in the mesh or only MASSLESS defined, then there is nothing to do.
   T_USED(8)=T_USED(8)+CURRENT_TIME()-TNOW
   RETURN
ENDIF

IF (LIQUID_DROPLETS) THEN

   ALLOCATE(PART_WARNING(NLP))
   PART_WARNING=0

   ! Working arrays

   MVAP_TOT => WORK7
   MVAP_TOT = 0._EB
   DO IW = 1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
      WC => WALL(IW)
      B1 => BOUNDARY_PROP1(WC%B1_INDEX)
      B2 => BOUNDARY_PROP2(WC%B2_INDEX)
      B2%WORK1 = B1%TMP_F
   ENDDO
   DO ICF = INTERNAL_CFACE_CELLS_LB+1,INTERNAL_CFACE_CELLS_LB+N_INTERNAL_CFACE_CELLS
      CFA => CFACE(ICF)
      B1 => BOUNDARY_PROP1(CFA%B1_INDEX)
      B2 => BOUNDARY_PROP2(CFA%B2_INDEX)
      B2%WORK1 = B1%TMP_F
   ENDDO

   RHO_INTERIM => WORK1 ; RHO_INTERIM(0:IBP1,0:JBP1,0:KBP1) = RHO(0:IBP1,0:JBP1,0:KBP1)
   TMP_INTERIM => WORK2 ; TMP_INTERIM(0:IBP1,0:JBP1,0:KBP1) = TMP(0:IBP1,0:JBP1,0:KBP1)
   ZZ_INTERIM  => SWORK1 ; ZZ_INTERIM(0:IBP1,0:JBP1,0:KBP1,1:) = ZZ(0:IBP1,0:JBP1,0:KBP1,1:)

ENDIF

! Keep a running average of surface mass and cooling

DO N_LPC=1,N_LAGRANGIAN_CLASSES
   LPC => LAGRANGIAN_PARTICLE_CLASS(N_LPC)
   IF (LPC%ARRAY_INDEX > 0) THEN
      IF (LPC%LIQUID_DROPLET) THEN
         DO IW = 1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
            WC => WALL(IW)
            B2 => BOUNDARY_PROP2(WC%B2_INDEX)
            B2%LP_TEMP(LPC%ARRAY_INDEX) = 0._EB
            B2%LP_CPUA(LPC%ARRAY_INDEX) = LPC%RUNNING_AVERAGE_FACTOR_WALL*B2%LP_CPUA(LPC%ARRAY_INDEX)
         ENDDO
         DO ICF = INTERNAL_CFACE_CELLS_LB+1,INTERNAL_CFACE_CELLS_LB+N_INTERNAL_CFACE_CELLS
            CFA => CFACE(ICF)
            B2 => BOUNDARY_PROP2(CFA%B2_INDEX)
            B2%LP_TEMP(LPC%ARRAY_INDEX) = 0._EB
            B2%LP_CPUA(LPC%ARRAY_INDEX) = LPC%RUNNING_AVERAGE_FACTOR_WALL*B2%LP_CPUA(LPC%ARRAY_INDEX)
         ENDDO
      ENDIF
      DO IW = 1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
         WC => WALL(IW)
         B2 => BOUNDARY_PROP2(WC%B2_INDEX)
         B2%LP_MPUA(LPC%ARRAY_INDEX) = LPC%RUNNING_AVERAGE_FACTOR_WALL*B2%LP_MPUA(LPC%ARRAY_INDEX)
      ENDDO
      DO ICF = INTERNAL_CFACE_CELLS_LB+1,INTERNAL_CFACE_CELLS_LB+N_INTERNAL_CFACE_CELLS
         CFA => CFACE(ICF)
         B2 => BOUNDARY_PROP2(CFA%B2_INDEX)
         B2%LP_MPUA(LPC%ARRAY_INDEX) = LPC%RUNNING_AVERAGE_FACTOR_WALL*B2%LP_MPUA(LPC%ARRAY_INDEX)
      ENDDO
   ENDIF
ENDDO

IF (.NOT. LIQUID_DROPLETS .AND. SOLID_PARTICLES) THEN
   ! If only solid particles, then no evaporation. Just update averages and remove burned away particles.
   CALL PARTICLE_RUNNING_AVERAGES
   CALL REMOVE_PARTICLES(T,NM)
   T_USED(8)=T_USED(8)+CURRENT_TIME()-TNOW
   RETURN
ENDIF

! Loop over all types of evaporative species

SPECIES_LOOP: DO Z_INDEX = 1,N_TRACKED_SPECIES

   Z_INDEX_A(1) = Z_INDEX
   ! Initialize quantities common to the evaporation index

   IF (.NOT.SPECIES_MIXTURE(Z_INDEX)%EVAPORATING) CYCLE SPECIES_LOOP
   Y_INDEX = MAXVAL(MAXLOC(SPECIES_MIXTURE(Z_INDEX)%VOLUME_FRACTION))
   SS => SPECIES(Y_INDEX)
   TMP_MELT = SS%TMP_MELT
   MW_DROP  = SS%MW
   CALL INTERPOLATE1D_UNIFORM(LBOUND(SS%H_V,1),SS%H_V,SS%H_V_REFERENCE_TEMPERATURE,H_V_REF)
   I_MELT   = INT(TMP_MELT)

   DO IW = 1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
      WC => WALL(IW)
      B2 => BOUNDARY_PROP2(WC%B2_INDEX)
      B2%WORK2 = 0._EB  ! Depth of liquid layer for a given surface cell
      B2%WORK3 = 0._EB  ! Total droplet/surface contact area per surface cell
   ENDDO
   DO ICF = INTERNAL_CFACE_CELLS_LB+1,INTERNAL_CFACE_CELLS_LB+N_INTERNAL_CFACE_CELLS
      CFA => CFACE(ICF)
      B2 => BOUNDARY_PROP2(CFA%B2_INDEX)
      B2%WORK2 = 0._EB  ! Depth of liquid layer for a given surface cell
      B2%WORK3 = 0._EB  ! Total droplet/surface contact area per surface cell
   ENDDO

   ! Loop through all PARTICLEs in the class and determine the depth of the liquid film on each surface cell

   SURFACE_DROPLETS = .FALSE.

   FILM_SUMMING_LOOP: DO IP=1,NLP
      LP  => LAGRANGIAN_PARTICLE(IP)
      LPC => LAGRANGIAN_PARTICLE_CLASS(LP%CLASS_INDEX)
      IF (.NOT.LPC%LIQUID_DROPLET) CYCLE FILM_SUMMING_LOOP
      IF (LPC%Z_INDEX/=Z_INDEX) CYCLE FILM_SUMMING_LOOP
      IF (LP%RADIUS<=LPC%KILL_RADIUS) CYCLE FILM_SUMMING_LOOP
      IF (LP%WALL_INDEX>0) THEN
         WC => WALL(LP%WALL_INDEX)
         B1 => BOUNDARY_PROP1(WC%B1_INDEX)
         B2 => BOUNDARY_PROP2(WC%B2_INDEX)
      ELSEIF (LP%CFACE_INDEX>0) THEN
         CFA => CFACE(LP%CFACE_INDEX)
         B1 => BOUNDARY_PROP1(CFA%B1_INDEX)
         B2 => BOUNDARY_PROP2(CFA%B2_INDEX)
      ELSE
         CYCLE FILM_SUMMING_LOOP
      ENDIF
      B2%WORK2 = B2%WORK2 + LP%PWT*FOTHPI*LP%RADIUS**3/B1%AREA  ! WORK2 stores the depth of liquid film
      B2%WORK3 = B2%WORK3 + LP%PWT*PI*(CR2*LP%RADIUS)**2        ! WORK3 stores total droplet/surface contact area
      SURFACE_DROPLETS = .TRUE.
   ENDDO FILM_SUMMING_LOOP

   ! If the total droplet/surface contact area is less than the wall cell area, do not assume a film layer.

   IF (SURFACE_DROPLETS) THEN
      DO IW = 1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
         WC => WALL(IW)
         B1 => BOUNDARY_PROP1(WC%B1_INDEX)
         B2 => BOUNDARY_PROP2(WC%B2_INDEX)
         IF (B2%WORK3<B1%AREA) B2%WORK2 = 0._EB
         B2%WORK3 = 0._EB  ! WORK3 now represents the sum of absorbed radiation at the WALL cell
      ENDDO
      DO ICF = INTERNAL_CFACE_CELLS_LB+1,INTERNAL_CFACE_CELLS_LB+N_INTERNAL_CFACE_CELLS
         CFA => CFACE(ICF)
         B1 => BOUNDARY_PROP1(CFA%B1_INDEX)
         B2 => BOUNDARY_PROP2(CFA%B2_INDEX)
         IF (B2%WORK3<B1%AREA) B2%WORK2 = 0._EB
         B2%WORK3 = 0._EB  ! WORK3 represents the sum of absorbed radiation at the CFACE
      ENDDO
   ENDIF

   ! Loop through all PARTICLEs within the class and determine mass/energy transfer

   PARTICLE_LOOP: DO IP=1,NLP

      LP  => LAGRANGIAN_PARTICLE(IP)
      LPC => LAGRANGIAN_PARTICLE_CLASS(LP%CLASS_INDEX)
      BC  => BOUNDARY_COORD(LP%BC_INDEX)

      IF (LPC%Z_INDEX/=Z_INDEX)     CYCLE PARTICLE_LOOP
      IF (.NOT.LPC%LIQUID_DROPLET)  CYCLE PARTICLE_LOOP

      LP_B1  => BOUNDARY_PROP1(LP%B1_INDEX)

      SF_FIXED = .FALSE.

      IF (LP%WALL_INDEX>0) THEN
         WC => WALL(LP%WALL_INDEX)
         B1 => BOUNDARY_PROP1(WC%B1_INDEX)
         B2 => BOUNDARY_PROP2(WC%B2_INDEX)
         SF => SURFACE(WC%SURF_INDEX)
         IF (WC%OD_INDEX>0) ONE_D => BOUNDARY_ONE_D(WC%OD_INDEX)
         IF (SF%THERMAL_BC_INDEX == CONVECTIVE_FLUX_BC .OR. SF%THERMAL_BC_INDEX == NET_FLUX_BC) SF_FIXED=.TRUE.
      ELSEIF (LP%CFACE_INDEX>0) THEN
         CFA => CFACE(LP%CFACE_INDEX)
         B1 => BOUNDARY_PROP1(CFA%B1_INDEX)
         B2 => BOUNDARY_PROP2(CFA%B2_INDEX)
         SF => SURFACE(CFA%SURF_INDEX)
         IF (CFA%OD_INDEX>0) ONE_D => BOUNDARY_ONE_D(CFA%OD_INDEX)
         IF (SF%THERMAL_BC_INDEX == CONVECTIVE_FLUX_BC .OR. SF%THERMAL_BC_INDEX == NET_FLUX_BC) SF_FIXED=.TRUE.
      ENDIF

      ! Determine the current coordinates of the particle

      II = BC%IIG
      JJ = BC%JJG
      KK = BC%KKG

      ! Determine how many sub-time step iterations are needed and then iterate over the time step.

      DT_SUBSTEP = DT
      DT_SUM = 0._EB
      WGT    = LP%PWT
      LP_B1%M_DOT_G_PP_ACTUAL = 0._EB
      Q_CON_SUM = 0._EB
      Q_RAD_SUM = 0._EB

      Y_COND = 0._EB

      TIME_ITERATION_LOOP: DO WHILE (DT_SUM < DT)

         KILL_RADIUS_CHECK: IF (LP%RADIUS>LPC%KILL_RADIUS) THEN

            ! Gas conditions

            ZZ_GET(1:N_TRACKED_SPECIES) = ZZ_INTERIM(II,JJ,KK,1:N_TRACKED_SPECIES)
            CALL GET_MASS_FRACTION_ALL(ZZ_GET,Y_ALL)
            Y_GAS = Y_ALL(Y_INDEX)
            IF (SPECIES_MIXTURE(Z_INDEX)%CONDENSATION_SMIX_INDEX > 0) &
               Y_COND = ZZ_GET(SPECIES_MIXTURE(Z_INDEX)%CONDENSATION_SMIX_INDEX)
            Y_GAS = Y_GAS - Y_COND
            TMP_G  = MAX(TMPMIN,TMP_INTERIM(II,JJ,KK))
            RHO_G  = RHO_INTERIM(II,JJ,KK)

            ! Initialize PARTICLE thermophysical data

            R_DROP   = LP%RADIUS
            FTPR     = FOTHPI * LPC%DENSITY
            M_DROP   = FTPR*R_DROP**3
            TMP_DROP = LP_B1%TMP_F
            T_BOIL_EFF = SS%TMP_V
            CALL GET_EQUIL_DATA(MW_DROP,TMP_DROP,PBAR(KK,PRESSURE_ZONE(II,JJ,KK)),H_V,H_V_A,T_BOIL_EFF,X_DROP,SS%H_V)
            I_BOIL   = INT(T_BOIL_EFF)

            IF (H_V < 0._EB) THEN
               WRITE(MESSAGE,'(A,A,A,I0,A,I0,A,I0)') 'Numerical instability in particle energy transport, H_V for ',TRIM(SS%ID), &
                  '. Mesh: ',NM,' Particle: ',IP,' Tag: ',LP%TAG
               CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.)
               DEALLOCATE(PART_WARNING)
               RETURN
            ENDIF

            IF (INT(TMP_DROP) < INT(T_BOIL_EFF)) THEN
               DH_V_A_DT = 0.5_EB*(SS%H_V(INT(TMP_DROP)+1) - SS%H_V(INT(TMP_DROP)))
            ELSE
               DH_V_A_DT = 0.5_EB*(SS%H_V(INT(T_BOIL_EFF)) - SS%H_V(INT(T_BOIL_EFF)-1))
            ENDIF

            M_GAS  = RHO_G/LP%RVC
            IF (DT_SUM <= TWO_EPSILON_EB) THEN
               TMP_G_OLD = TMP_G
               M_GAS_OLD = M_GAS
            ENDIF
            CALL GET_SPECIFIC_HEAT(ZZ_GET,CP,TMP_G)
            CALL GET_ENTHALPY(ZZ_GET,H_G_OLD,TMP_G)
            CP_BAR= H_G_OLD / TMP_G
            H_G_OLD = H_G_OLD * M_GAS
            H_S_G_OLD = CP * TMP_G * M_GAS
            M_VAP_MAX = (0.33_EB * M_GAS - MVAP_TOT(II,JJ,KK)) / WGT ! limit to avoid diveregence errors
            U2 = 0.5_EB*(U(II,JJ,KK)+U(II-1,JJ,KK))
            V2 = 0.5_EB*(V(II,JJ,KK)+V(II,JJ-1,KK))
            W2 = 0.5_EB*(W(II,JJ,KK)+W(II,JJ,KK-1))

            SOLID_OR_GAS_PHASE_1: IF (BC%IOR/=0 .AND. (LP%WALL_INDEX>0 .OR. LP%CFACE_INDEX>0) .AND. .NOT. SF_FIXED) THEN
               SELECT CASE(ABS(BC%IOR))
                  CASE(1)
                     VEL = SQRT(V2**2+W2**2)
                  CASE(2)
                     VEL = SQRT(U2**2+W2**2)
                  CASE(3)
                     VEL = SQRT(U2**2+V2**2)
               END SELECT
               IF (B2%WORK2>0._EB) THEN
                  A_DROP = M_DROP/(B2%WORK2*LPC%DENSITY)  ! WORK2 is the computed liquid film layer depth
               ELSE
                  A_DROP = PI*(CR2*LP%RADIUS)**2          ! If no assumed film layer, use the droplet/surface contact area directly
               ENDIF
               Q_DOT_RAD = MIN(A_DROP,B1%AREA/WGT)*B1%Q_RAD_IN
               B2%WORK3 = B2%WORK3 + Q_DOT_RAD*WGT
               TMP_WALL = MAX(TMPMIN,B2%WORK1)
            ELSE SOLID_OR_GAS_PHASE_1
               VEL = SQRT((U2-LP%U)**2+(V2-LP%V)**2+(W2-LP%W)**2)
               A_DROP   = 4._EB*PI*R_DROP**2
               IF (SUM(AVG_DROP_AREA(II,JJ,KK,:))>0._EB) THEN
                  Q_DOT_RAD = (QR_W(II,JJ,KK)/SUM(AVG_DROP_AREA(II,JJ,KK,:)))*(A_DROP/4._EB)
               ELSE
                  Q_DOT_RAD = 0._EB
               ENDIF
            ENDIF SOLID_OR_GAS_PHASE_1

            ! Determine the ratio of molecular weights between the gas and droplet vapor

            CALL GET_MW_RATIO(Y_INDEX,MW_RATIO,Y_IN=Y_ALL)

            ! Get actual MW of current gas for D_SOURCE

            CALL GET_MOLECULAR_WEIGHT(ZZ_GET,MW_GAS)

            ! Decide whether to evporate the entire droplet or only part of it

            BOIL_ALL: IF (Q_DOT_RAD*DT_SUBSTEP > M_DROP*H_V) THEN  ! boil it all

               M_VAP = M_DROP
               Q_RAD = M_VAP*H_V
               M_VAP = LPC%ADJUST_EVAPORATION*M_VAP

               ! Compute contribution to the divergence

               ZZ_GET = 0._EB
               ZZ_GET(Z_INDEX) = 1._EB
               CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S_B,TMP_DROP)
               CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,TMP_G)
               DELTA_H_G = H_S_B - H_S
               D_SOURCE(II,JJ,KK) = D_SOURCE(II,JJ,KK) + (MW_GAS/MW_DROP*M_VAP/M_GAS + (M_VAP*DELTA_H_G)/H_S_G_OLD) * WGT / DT
               M_DOT_PPP(II,JJ,KK,Z_INDEX) = M_DOT_PPP(II,JJ,KK,Z_INDEX) + M_VAP*LP%RVC*WGT/DT
               LP_B1%M_DOT_G_PP_ACTUAL(Z_INDEX) = LP_B1%M_DOT_G_PP_ACTUAL(Z_INDEX) + M_VAP/(A_DROP*DT)

               ! Add energy losses and gains to overall energy budget array

               Q_DOT(8) = Q_DOT(8) - Q_RAD*WGT/DT  ! Q_PART
               Q_DOT(4) = Q_DOT(4) + M_VAP*H_S_B*WGT/DT                       ! Q_CONV
               Q_DOT(3) = Q_DOT(3) + Q_RAD*WGT/DT                             ! Q_RADI

               IF (LPC%Z_INDEX>0) M_DOT(LPC%Z_INDEX) = M_DOT(LPC%Z_INDEX) + WGT*M_VAP/DT/LPC%ADJUST_EVAPORATION

               ! Keep track of total mass evaporated in cell

               MVAP_TOT(II,JJ,KK) = MVAP_TOT(II,JJ,KK) + WGT*M_VAP

               ! Force removal of droplet

               LP%RADIUS = 0._EB

            ELSE BOIL_ALL  ! boil only some of the liquid

               CALL INTERPOLATE1D_UNIFORM(LBOUND(SS%C_P_L,1),SS%C_P_L,TMP_DROP,C_DROP)
               CALL INTERPOLATE1D_UNIFORM(LBOUND(SS%C_P_L_BAR,1),SS%C_P_L_BAR,TMP_DROP,H_L)
               H_L = H_L*TMP_DROP
               H_D_OLD  = H_L*M_DROP
               ZZ_GET(1:N_TRACKED_SPECIES) = ZZ_INTERIM(II,JJ,KK,1:N_TRACKED_SPECIES)

               ! Compute equilibrium PARTICLE vapor mass fraction, Y_DROP, and its derivative w.r.t. PARTICLE temperature

               Y_DROP  = X_DROP/(MW_RATIO + (1._EB-MW_RATIO)*X_DROP)
               IF (Y_DROP > Y_GAS) LP_B1%B_NUMBER = (Y_DROP - Y_GAS) / MAX(DY_MIN_BLOWING,1._EB-Y_DROP)

               ! Compute temperature deriviative of the vapor mass fraction

               DHOR     = H_V_A*MW_DROP/R0
               DYDT = (MW_RATIO/(X_DROP*(1._EB-MW_RATIO)+MW_RATIO)**2) * &
                      (DHOR*X_DROP/TMP_DROP**2+(1._EB/T_BOIL_EFF-1._EB/TMP_DROP)*DH_V_A_DT*MW_DROP/R0)

               ! Set variables for heat transfer on solid

               Y_GAS_A = Y_GAS
               Y_DROP_A = Y_DROP

               SOLID_OR_GAS_PHASE_2: IF (BC%IOR/=0 .AND. (LP%WALL_INDEX>0 .OR. LP%CFACE_INDEX>0) .AND. .NOT. SF_FIXED) THEN

                  CALL GET_FILM_PROPERTIES(1,PLATE_FILM_FACTOR,Y_DROP_A,Y_GAS_A,Z_INDEX_A,TMP_DROP,TMP_G,ZZ_GET, &
                                           PBAR(KK,PRESSURE_ZONE(II,JJ,KK)),TMP_FILM,MU_FILM,&
                                           K_FILM,CP_FILM,D_FILM,RHO_FILM,PR_FILM,SC_FILM)

                  ! Compute mcbar = rho_w a_w cp_w dx_w for first wall cell for limiting convective heat transfer

                  IF (SF%THERMAL_BC_INDEX==THERMALLY_THICK) THEN
                     RHOCBAR = 0._EB
                     ITMP = MIN(I_MAX_TEMP,NINT(TMP_WALL))
                     DO NMAT=1,SF%N_MATL
                        IF (ONE_D%MATL_COMP(NMAT)%RHO(1)<=TWO_EPSILON_EB) CYCLE
                        RHOCBAR = RHOCBAR + ONE_D%MATL_COMP(NMAT)%RHO(1)*MATERIAL(SF%MATL_INDEX(NMAT))%C_S(ITMP)
                        RHOCBAR = RHOCBAR + LPC%DENSITY*MATERIAL(SF%MATL_INDEX(NMAT))%C_S(ITMP)
                     ENDDO
                     MCBAR = RHOCBAR*B1%AREA*(ONE_D%X(1)-ONE_D%X(0))
                     ARRAY_CASE = 3
                     IF (MCBAR <= TWO_EPSILON_EB) THEN
                        MCBAR = -1._EB
                        ARRAY_CASE = 2
                     ENDIF
                  ELSEIF (SF%THERMAL_BC_INDEX==CONVECTIVE_FLUX_BC .OR. SF%THERMAL_BC_INDEX==NET_FLUX_BC) THEN
                     ARRAY_CASE = 1
                  ELSE
                     MCBAR = -1._EB
                     ARRAY_CASE = 2
                  ENDIF

                  ! Compute heat transfer coefficient between liquid film and solid wall, H_WALL

                  IF (LPC%HEAT_TRANSFER_COEFFICIENT_SOLID<0._EB .AND. ARRAY_CASE>1) THEN
                     LENGTH = 1._EB
                     SELECT CASE(ABS(BC%IOR))
                        CASE(0:2) ;  VEL_LIQUID = LPC%VERTICAL_VELOCITY
                        CASE(3)   ;  VEL_LIQUID = LPC%HORIZONTAL_VELOCITY
                     END SELECT
                     RE_LIQUID = SS%DENSITY_LIQUID*VEL_LIQUID*LENGTH/SS%MU_LIQUID
                     CALL FORCED_CONVECTION_MODEL(NUSSELT_LIQUID,RE_LIQUID,SS%PR_LIQUID**ONTH,SURF_CARTESIAN)
                     H_WALL = MAX(100._EB,NUSSELT_LIQUID*SS%K_LIQUID/LENGTH)
                  ELSEIF (ARRAY_CASE>1) THEN
                     H_WALL = LPC%HEAT_TRANSFER_COEFFICIENT_SOLID
                  ELSE
                     H_WALL = 0._EB
                  ENDIF

                  ! Convective heat transfer coefficient between droplet/film and gas, H_HEAT
                  ! Incropera and Dewitt, Fundamentals of Heat and Mass Transfer, 7th Edition

                  LENGTH   = 1._EB
                  RE_L     = RHO_G*VEL*LENGTH/MU_FILM
                  CALL FORCED_CONVECTION_MODEL(NUSSELT,RE_L,PR_FILM**ONTH,SURF_CARTESIAN)
                  H_HEAT   = NUSSELT*K_FILM/LENGTH

                  ! Compute mass transfer coefficient

                  IF (Y_DROP<=Y_GAS) THEN
                     H_MASS = 0._EB
                  ELSE
                     CALL FORCED_CONVECTION_MODEL(SHERWOOD,RE_L,SC_FILM**ONTH,SURF_CARTESIAN)
                     !M# expressions taken from Sazhin, Prog in Energy and Comb Sci 32 (2006) 162-214
                     SELECT CASE(LPC%EVAP_MODEL)
                        CASE(RM_NO_B)
                           H_MASS   = SHERWOOD*D_FILM/LENGTH
                        CASE(RM_B,RM_LEWIS_B) !Sazhin M0 - M1
                           H_MASS   = SHERWOOD*D_FILM/LENGTH*LOG(1._EB+LP_B1%B_NUMBER)/LP_B1%B_NUMBER
                        CASE(RM_FL_LEWIS_B) !Sazhin M2
                           H_MASS   = SHERWOOD*D_FILM/LENGTH*LOG(1._EB+LP_B1%B_NUMBER)/ &
                                     (LP_B1%B_NUMBER*F_B(LP_B1%B_NUMBER))
                     END SELECT
                  ENDIF

               ELSE SOLID_OR_GAS_PHASE_2

                  CALL GET_FILM_PROPERTIES(1,SPHERE_FILM_FACTOR,Y_DROP_A,Y_GAS_A,Z_INDEX_A,TMP_DROP,TMP_G,ZZ_GET,&
                                           PBAR(KK,PRESSURE_ZONE(II,JJ,KK)),TMP_FILM,MU_FILM,&
                                           K_FILM,CP_FILM,D_FILM,RHO_FILM,PR_FILM,SC_FILM)

                  SH_FAC_GAS = 0.6_EB*SC_FILM**ONTH
                  NU_FAC_GAS = 0.6_EB*PR_FILM**ONTH
                  LENGTH = 2._EB*R_DROP
                  RE_L = RHO_FILM*VEL*LENGTH/MU_FILM
                  CALL DROPLET_H_MASS_H_HEAT_GAS(H_MASS,H_HEAT,D_FILM,K_FILM,CP_FILM,RHO_FILM,LENGTH,Y_DROP,Y_GAS,&
                                                 LP_B1%B_NUMBER,NU_FAC_GAS,SH_FAC_GAS,RE_L,TMP_FILM,ZZ_GET,Z_INDEX,&
                                                 LPC%EVAP_MODEL)
                  H_WALL   = 0._EB
                  TMP_WALL = TMPA
                  ARRAY_CASE = 1

               ENDIF SOLID_OR_GAS_PHASE_2

               ! User-specified heat and mass transfer coefficients over-ride computed ones

               IF (LPC%HEAT_TRANSFER_COEFFICIENT_GAS>=0._EB) H_HEAT=LPC%HEAT_TRANSFER_COEFFICIENT_GAS
               IF (LPC%MASS_TRANSFER_COEFFICIENT>=0._EB)     H_MASS=LPC%MASS_TRANSFER_COEFFICIENT

               ! Build and solve implicit arrays for updating particle, gas, and wall temperatures

               ITMP = INT(TMP_DROP)
               H1 = H_SENS_Z(ITMP,Z_INDEX)+(TMP_DROP-REAL(ITMP,EB))*(H_SENS_Z(ITMP+1,Z_INDEX)-H_SENS_Z(ITMP,Z_INDEX))
               ITMP = INT(TMP_G)
               H2 = H_SENS_Z(ITMP,Z_INDEX)+(TMP_G-REAL(ITMP,EB))*(H_SENS_Z(ITMP+1,Z_INDEX)-H_SENS_Z(ITMP,Z_INDEX))

               AGHRHO = A_DROP*H_MASS*RHO_FILM/(1._EB+0.5_EB*LP%RVC*DT_SUBSTEP*A_DROP*WGT*H_MASS*RHO_FILM*(1._EB-Y_GAS)/RHO_G)

               DTOG = DT_SUBSTEP*WGT/(M_GAS*CP)
               DTGOG = 0.5_EB*DTOG*A_DROP*H_HEAT

               DAHVHLDY = DTOG*AGHRHO*(H1-H2)*(Y_DROP-Y_GAS)
               DADYDTHVHL=0.5_EB*DTOG*AGHRHO*(H1-H2)*DYDT

               DTOP = DT_SUBSTEP/(M_DROP*C_DROP)
               DTGOP = 0.5_EB*DTOP*A_DROP*H_HEAT

               DADYHV = DTOP*AGHRHO*H_V*(Y_DROP-Y_GAS)
               DADYDTHV=0.5_EB*DTOP*AGHRHO*DYDT*H_V

               SELECT CASE (ARRAY_CASE)
                  CASE(1) ! Gas Only
                     A_COL(1) = 1._EB+DTGOG
                     B_COL(1) = -(DTGOG+DADYDTHVHL)
                     A_COL(2) = -DTGOP
                     B_COL(2) = 1._EB+DTGOP+DADYDTHV
                     D_VEC(1) = (1._EB-DTGOG)*TMP_G+(DTGOG-DADYDTHVHL)*TMP_DROP+DAHVHLDY
                     D_VEC(2) = DTGOP*TMP_G+(1-DTGOP+DADYDTHV)*TMP_DROP-DADYHV+DTOP*Q_DOT_RAD
                     TMP_DROP_NEW = -(A_COL(2)*D_VEC(1)-A_COL(1)*D_VEC(2))/(A_COL(1)*B_COL(2)-B_COL(1)*A_COL(2))
                     TMP_G_NEW = (D_VEC(1)-B_COL(1)*TMP_DROP_NEW)/A_COL(1)
                     TMP_WALL_NEW = TMP_WALL
                  CASE(2) ! Const Temp Wall

                     DTWOP = 0.5_EB*DTOP*A_DROP*H_WALL

                     A_COL(1) = 1._EB+DTGOG
                     B_COL(1) = -(DTGOG+DADYDTHVHL)
                     A_COL(2) = -DTGOP
                     B_COL(2) = 1._EB+DTGOP+DTWOP+DADYDTHV
                     D_VEC(1) = (1._EB-DTGOG)*TMP_G+(DTGOG-DADYDTHVHL)*TMP_DROP+DAHVHLDY
                     D_VEC(2) = DTGOP*TMP_G+(1-DTGOP-DTWOP+DADYDTHV)*TMP_DROP+2._EB*DTWOP*TMP_WALL-DADYHV+DTOP*Q_DOT_RAD
                     TMP_DROP_NEW = -(A_COL(2)*D_VEC(1)-A_COL(1)*D_VEC(2))/(A_COL(1)*B_COL(2)-B_COL(1)*A_COL(2))
                     TMP_G_NEW = (D_VEC(1)-B_COL(1)*TMP_DROP_NEW)/A_COL(1)
                     TMP_WALL_NEW = TMP_WALL
                  CASE(3) ! 1D Wall

                     DTWOP = 0.5_EB*DTOP*A_DROP*H_WALL
                     DTWOW = DT_SUBSTEP*A_DROP*WGT*H_WALL/(2._EB*MCBAR)

                     A_COL(1) = 1._EB+DTGOG
                     B_COL(1) = -(DTGOG+DADYDTHVHL)
                     C_COL(1) = 0._EB
                     A_COL(2) = 0._EB
                     B_COL(2) = -DTWOW
                     C_COL(2) = 1._EB+DTWOW
                     A_COL(3) = -DTGOP
                     B_COL(3) = 1._EB+DTGOP+DTWOP+DADYDTHV
                     C_COL(3) = -DTWOP
                     D_VEC(1) = (1._EB-DTGOG)*TMP_G+(DTGOG-DADYDTHVHL)*TMP_DROP+DAHVHLDY
                     D_VEC(2) = DTWOW*TMP_DROP+(1._EB-DTWOW)*TMP_WALL
                     D_VEC(3) = DTGOP*TMP_G+(1-DTGOP-DTWOP+DADYDTHV)*TMP_DROP+DTWOP*TMP_WALL-DADYHV+DTOP*Q_DOT_RAD
                     TMP_DROP_NEW = (A_COL(1)*(D_VEC(2)*C_COL(3)-C_COL(2)*D_VEC(3))+D_VEC(1)*C_COL(2)*A_COL(3))/&
                                    (A_COL(1)*(B_COL(2)*C_COL(3)-C_COL(2)*B_COL(3))+B_COL(1)*C_COL(2)*A_COL(3))
                     TMP_WALL_NEW = (D_VEC(2)-B_COL(2)*TMP_DROP_NEW)/C_COL(2)
                     TMP_G_NEW = (D_VEC(1)-B_COL(1)*TMP_DROP_NEW)/A_COL(1)
               END SELECT
               M_VAP = MAX(0._EB,MIN(M_DROP, DT_SUBSTEP * AGHRHO * (Y_DROP-Y_GAS+0.5_EB*DYDT*(TMP_DROP_NEW-TMP_DROP))))

               ! Compute the total amount of heat extracted from the gas, wall and radiative fields

               Q_RAD      = DT_SUBSTEP*Q_DOT_RAD
               Q_CON_GAS  = DT_SUBSTEP*A_DROP*H_HEAT*0.5_EB*(TMP_G   +TMP_G_NEW   -TMP_DROP-TMP_DROP_NEW)
               Q_CON_WALL = DT_SUBSTEP*A_DROP*H_WALL*0.5_EB*(TMP_WALL+TMP_WALL_NEW-TMP_DROP-TMP_DROP_NEW)
               Q_TOT = Q_RAD+Q_CON_GAS+Q_CON_WALL
               IF (Q_TOT >= M_DROP*H_V) M_VAP = M_DROP

               ! Adjust drop temperature for variable liquid CP

               EVAP_ALL: IF (M_VAP < M_DROP) THEN
                  TMP_DROP_NEW = TMP_DROP + (Q_TOT - M_VAP * H_V)/(C_DROP * (M_DROP - M_VAP))
                  ITMP = NINT(TMP_DROP)
                  ITMP2 = MIN(I_BOIL,MAX(I_MELT,NINT(TMP_DROP_NEW)))
                  IF (ITMP/=ITMP2) THEN
                     C_DROP2 = SUM(SS%C_P_L(MIN(ITMP,ITMP2):MAX(ITMP,ITMP2)))/REAL(ABS(ITMP2-ITMP)+1,EB)
                     TMP_DROP_NEW = TMP_DROP + (Q_TOT - M_VAP * H_V)/(C_DROP2 * (M_DROP - M_VAP))
                  ENDIF

                  IF (TMP_DROP_NEW<TMP_MELT) THEN
                     ! If the PARTICLE temperature drops below its freezing point, just reset it
                     TMP_DROP_NEW = TMP_MELT
                  ELSEIF (TMP_DROP_NEW>T_BOIL_EFF) THEN
                     ! If the PARTICLE temperature reaches boiling, use only enough energy from gas to vaporize liquid
                     ITMP = NINT(TMP_DROP)
                     C_DROP2 = SUM(SS%C_P_L(MIN(ITMP,I_BOIL):MAX(ITMP,ITMP2)))/REAL(ABS(I_BOIL-ITMP)+1,EB)
                     M_VAP  = MIN(M_DROP,(Q_TOT-M_DROP*C_DROP2*(T_BOIL_EFF-TMP_DROP))/(H_V-C_DROP2*(T_BOIL_EFF-TMP_DROP)))
                     IF (M_VAP == M_DROP) THEN
                        Q_FRAC = M_VAP*H_V/Q_TOT
                        Q_RAD = Q_RAD*Q_FRAC
                        Q_CON_GAS = Q_CON_GAS*Q_FRAC
                        Q_CON_WALL = Q_CON_WALL*Q_FRAC
                     ENDIF
                     TMP_DROP_NEW = T_BOIL_EFF
                  ENDIF
               ELSE EVAP_ALL
                  Q_FRAC = M_VAP*H_V/Q_TOT
                  Q_RAD = Q_RAD*Q_FRAC
                  Q_CON_GAS = Q_CON_GAS*Q_FRAC
                  Q_CON_WALL = Q_CON_WALL*Q_FRAC
                  TMP_DROP_NEW = T_BOIL_EFF
               ENDIF EVAP_ALL

               IF (BC%IOR/=0 .AND. (LP%WALL_INDEX>0 .OR. LP%CFACE_INDEX>0) .AND. .NOT. SF_FIXED) &
                  TMP_WALL_NEW=TMP_WALL-Q_CON_WALL/MCBAR
               M_DROP = M_DROP - M_VAP

               ! Add fuel evaporation rate to running counter and adjust mass of evaporated fuel to account for different
               ! Heat of Combustion between fuel PARTICLE and gas

               M_VAP = LPC%ADJUST_EVAPORATION*M_VAP

               M_GAS_NEW = M_GAS + WGT*M_VAP
               Y_GAS_NEW = (Y_GAS*M_GAS+WGT*M_VAP)/M_GAS_NEW
               ZZ_GET2 = M_GAS/M_GAS_NEW*ZZ_GET
               ZZ_GET2(Z_INDEX) = ZZ_GET2(Z_INDEX) + WGT*M_VAP/M_GAS_NEW

               T_BOIL_EFF = SS%TMP_V
               CALL GET_EQUIL_DATA(MW_DROP,TMP_DROP,PBAR(KK,PRESSURE_ZONE(II,JJ,KK)),H_V2,H_V_A,T_BOIL_EFF,X_EQUIL,SS%H_V)
               Y_EQUIL  = X_EQUIL/(MW_RATIO + (1._EB-MW_RATIO)*X_EQUIL)

               ! Limit drop temperature decrease

               IF (TMP_DROP_NEW < 0.9999_EB*TMP_DROP .AND. 1.05_EB*Y_EQUIL < Y_GAS_NEW .AND. M_VAP > 0._EB) THEN
                  IF (Y_GAS_NEW - Y_GAS > 1.E-7_EB) THEN
                     DT_SUBSTEP = DT_SUBSTEP * 0.5_EB
                     IF (DT_SUBSTEP <= 0.00001_EB*DT) THEN
                        DT_SUBSTEP = DT_SUBSTEP * 2.0_EB
                        IF (.NOT. BTEST(PART_WARNING(IP),1)) THEN
                           WRITE(LU_ERR,'(A,I0,A,I0,A,I0)') 'WARNING Y_EQ < Y_G_N. Mesh: ',NM,' Particle: ',IP,' Tag: ',LP%TAG
                           PART_WARNING(IP) = IBSET(PART_WARNING(IP),1)
                        ENDIF
                     ELSE
                        CYCLE TIME_ITERATION_LOOP
                     ENDIF
                  ENDIF
               ENDIF

               ! Limit supersaturation

               IF (Y_GAS < Y_EQUIL .AND. M_VAP > 0._EB) THEN
                  IF (Y_GAS_NEW/Y_EQUIL > 1.01_EB) THEN
                     DT_SUBSTEP = DT_SUBSTEP * 0.5_EB
                     IF (DT_SUBSTEP <= 0.00001_EB*DT) THEN
                        DT_SUBSTEP = DT_SUBSTEP * 2.0_EB
                        IF (.NOT. BTEST(PART_WARNING(IP),2)) THEN
                           WRITE(LU_ERR,'(A,I0,A,I0,A,I0)') 'WARNING Y_G_N > Y_EQ. Mesh: ',NM,'Particle: ',IP,' Tag: ',LP%TAG
                           PART_WARNING(IP) = IBSET(PART_WARNING(IP),2)
                        ENDIF
                     ELSE
                     CYCLE TIME_ITERATION_LOOP
                     ENDIF
                  ENDIF
               ENDIF

               ! Update gas temperature

               CALL INTERPOLATE1D_UNIFORM(LBOUND(SS%C_P_L_BAR,1),SS%C_P_L_BAR,TMP_DROP_NEW,H_L)
               H_NEW = H_G_OLD + (H_D_OLD - M_DROP*TMP_DROP_NEW*H_L + Q_CON_WALL + Q_RAD)*WGT
               TMP_G_I = TMP_G
               TMP_G_NEW = TMP_G
               CALL GET_TEMPERATURE(TMP_G_NEW,H_NEW/M_GAS_NEW,ZZ_GET2)
               IF (TMP_G_NEW < 0._EB) THEN
                  DT_SUBSTEP = DT_SUBSTEP * 0.5_EB
                  IF (DT_SUBSTEP <= 0.00001_EB*DT) THEN
                     DT_SUBSTEP = DT_SUBSTEP * 2.0_EB
                     TMP_G_NEW = 1._EB
                     IF (.NOT. BTEST(PART_WARNING(IP),3)) THEN
                        WRITE(LU_ERR,'(A,I0,A,I0,A,I0)') 'WARNING TMP_G_N < 0. Mesh: ',NM,'Particle: ',IP,' Tag: ',LP%TAG
                        PART_WARNING(IP) = IBSET(PART_WARNING(IP),3)
                     ENDIF
                  ELSE
                     CYCLE TIME_ITERATION_LOOP
                  ENDIF
               ENDIF

               ! Limit gas temperature change

               IF (ABS(TMP_G_NEW/TMP_G - 1._EB) > 0.05_EB) THEN
                  DT_SUBSTEP = DT_SUBSTEP * 0.5_EB
                  IF (DT_SUBSTEP <= 0.00001_EB*DT) THEN
                     DT_SUBSTEP = DT_SUBSTEP * 2.0_EB
                     IF (.NOT. BTEST(PART_WARNING(IP),4)) THEN
                        WRITE(LU_ERR,'(A,I0,A,I0,A,I0)') 'WARNING Delta TMP_G. Mesh: ',NM,'Particle: ',IP,' Tag: ',LP%TAG
                        PART_WARNING(IP) = IBSET(PART_WARNING(IP),4)
                     ENDIF
                  ELSE
                     CYCLE TIME_ITERATION_LOOP
                  ENDIF
               ENDIF

               CALL INTERPOLATE1D_UNIFORM(LBOUND(SS%C_P_L_BAR,1),SS%C_P_L_BAR,TMP_DROP_NEW,H_L)
               H_NEW= H_NEW + M_DROP*WGT*TMP_DROP_NEW*H_L

               IF (TMP_G_NEW < 0.9999_EB*TMP_G .AND. TMP_G_NEW < TMP_DROP_NEW) THEN
                  DT_SUBSTEP = DT_SUBSTEP * 0.5_EB
                  IF (DT_SUBSTEP <= 0.00001_EB*DT) THEN
                     DT_SUBSTEP = DT_SUBSTEP * 2.0_EB
                     IF (.NOT. BTEST(PART_WARNING(IP),5)) THEN
                        WRITE(LU_ERR,'(A,I0,A,I0,A,I0)') 'WARNING TMP_G_N < TMP_D_N. Mesh: ',NM,'Particle: ',IP,' Tag: ',LP%TAG
                        PART_WARNING(IP) = IBSET(PART_WARNING(IP),5)
                     ENDIF
                  ELSE
                     CYCLE TIME_ITERATION_LOOP
                  ENDIF
               ENDIF

               ! Update gas cell density, temperature, and mass fractions

               RHO_INTERIM(II,JJ,KK) = M_GAS_NEW*LP%RVC
               ZZ_INTERIM(II,JJ,KK,1:N_TRACKED_SPECIES) = ZZ_GET2(1:N_TRACKED_SPECIES)
               TMP_INTERIM(II,JJ,KK) = TMP_G_NEW
               IF (BC%IOR/=0 .AND. (LP%WALL_INDEX>0 .OR. LP%CFACE_INDEX>0) .AND. .NOT. SF_FIXED) B2%WORK1 = TMP_WALL_NEW

               ! Compute contribution to the divergence

               ZZ_GET = 0._EB
               ZZ_GET(Z_INDEX) = 1._EB
               CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S_B,TMP_DROP)
               CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,TMP_G)
               DELTA_H_G = H_S_B - H_S
               D_SOURCE(II,JJ,KK) = D_SOURCE(II,JJ,KK) + &
                                    (MW_GAS/MW_DROP*M_VAP/M_GAS + (M_VAP*DELTA_H_G - Q_CON_GAS)/H_S_G_OLD) * WGT / DT
               M_DOT_PPP(II,JJ,KK,Z_INDEX) = M_DOT_PPP(II,JJ,KK,Z_INDEX) + M_VAP*LP%RVC*WGT/DT
               LP_B1%M_DOT_G_PP_ACTUAL(Z_INDEX) = LP_B1%M_DOT_G_PP_ACTUAL(Z_INDEX) + M_VAP/(A_DROP*DT)

               ! Add stability checks

               IF (PARTICLE_CFL) PART_UVWMAX = MAX(PART_UVWMAX,H_MASS,H_HEAT/(RHO_G*CP))

               ! Add energy losses and gains to overall energy budget array

               Q_DOT(8) = Q_DOT(8) - (Q_CON_GAS + Q_CON_WALL + Q_RAD)*WGT/DT  ! Q_PART
               Q_DOT(4) = Q_DOT(4) + M_VAP*H_S_B*WGT/DT                       ! Q_CONV
               Q_DOT(3) = Q_DOT(3) + Q_RAD*WGT/DT                             ! Q_RADI
               Q_DOT(5) = Q_DOT(5) + Q_CON_WALL*WGT/DT                        ! Q_COND

               IF (LPC%Z_INDEX>0) M_DOT(LPC%Z_INDEX) = M_DOT(LPC%Z_INDEX) + WGT*M_VAP/DT/LPC%ADJUST_EVAPORATION

               ! Keep track of total mass evaporated in cell

               MVAP_TOT(II,JJ,KK) = MVAP_TOT(II,JJ,KK) + WGT*M_VAP

               ! Update PARTICLE quantities

               LP%RADIUS    = (M_DROP/FTPR)**ONTH
               LP_B1%TMP_F  = TMP_DROP_NEW
               LP%MASS      = M_DROP

               ! Sum convection and radiation for surfaces
               IF (LP%WALL_INDEX>0 .OR. LP%CFACE_INDEX>0) THEN
                  Q_CON_SUM = Q_CON_SUM + Q_CON_WALL
                  Q_RAD_SUM = Q_RAD_SUM + Q_DOT_RAD*DT_SUBSTEP
               ENDIF

            ENDIF BOIL_ALL

         ENDIF KILL_RADIUS_CHECK

         ! Get out of the loop if the PARTICLE has evaporated completely

         IF (LP%RADIUS<=LPC%KILL_RADIUS) CYCLE PARTICLE_LOOP

         DT_SUM = DT_SUM + DT_SUBSTEP
         DT_SUBSTEP = MIN(DT-DT_SUM,DT_SUBSTEP * 1.5_EB)

      ENDDO TIME_ITERATION_LOOP

      ! Compute surface cooling

      IF (LP%WALL_INDEX>0 .OR. LP%CFACE_INDEX>0) THEN
         R_DROP = LP%RADIUS
         A_DROP = WGT*PI*R_DROP**2
         B2%LP_TEMP(LPC%ARRAY_INDEX) = B2%LP_TEMP(LPC%ARRAY_INDEX) + A_DROP*0.5_EB*(TMP_DROP+TMP_DROP_NEW)
         B2%LP_CPUA(LPC%ARRAY_INDEX) = B2%LP_CPUA(LPC%ARRAY_INDEX) + &
                                       (1._EB-LPC%RUNNING_AVERAGE_FACTOR_WALL)*WGT*Q_CON_WALL/(B1%AREA*DT)
      ENDIF

   ENDDO PARTICLE_LOOP

   ! Subtract off the incoming radiation absorbed by all surface droplets of this type on each wall cell (B2%WORK3).

   IF (SURFACE_DROPLETS) THEN
      DO IW = 1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
         WC => WALL(IW)
         B1 => BOUNDARY_PROP1(WC%B1_INDEX)
         B2 => BOUNDARY_PROP2(WC%B2_INDEX)
         B1%Q_RAD_IN = MAX(0._EB,(B1%AREA*B1%Q_RAD_IN-B2%WORK3)/(B1%AREA))
      ENDDO
      DO ICF = INTERNAL_CFACE_CELLS_LB+1,INTERNAL_CFACE_CELLS_LB+N_INTERNAL_CFACE_CELLS
         CFA => CFACE(ICF)
         B1 => BOUNDARY_PROP1(CFA%B1_INDEX)
         B2 => BOUNDARY_PROP2(CFA%B2_INDEX)
         B1%Q_RAD_IN = MAX(0._EB,(B1%AREA*B1%Q_RAD_IN-B2%WORK3)/(B1%AREA))
      ENDDO
   ENDIF

ENDDO SPECIES_LOOP

DEALLOCATE(PART_WARNING)

! Sum up various quantities used in running averages

CALL PARTICLE_RUNNING_AVERAGES

! Remove PARTICLEs that have completely evaporated

CALL REMOVE_PARTICLES(T,NM)

T_USED(8)=T_USED(8)+CURRENT_TIME()-TNOW

CONTAINS


!> \brief Sum up various quantities used in running averages

SUBROUTINE PARTICLE_RUNNING_AVERAGES

INTEGER :: I,J,K
REAL(EB) :: DEN_ADD,AREA_ADD
REAL(EB), POINTER, DIMENSION(:,:,:) :: DROP_DEN,DROP_RAD,DROP_TMP,DROP_AREA

DROP_AREA => WORK1
DROP_DEN => WORK4
DROP_RAD => WORK5
DROP_TMP => WORK6

PART_CLASS_SUM_LOOP: DO N_LPC = 1,N_LAGRANGIAN_CLASSES

   LPC => LAGRANGIAN_PARTICLE_CLASS(N_LPC)
   IF (LPC%MASSLESS_TRACER .OR. LPC%MASSLESS_TARGET) CYCLE PART_CLASS_SUM_LOOP

   DROP_DEN = 0._EB
   DROP_TMP = 0._EB
   DROP_RAD = 0._EB
   DROP_AREA = 0._EB

   PARTICLE_LOOP_2: DO IP=1,NLP

      LP => LAGRANGIAN_PARTICLE(IP)
      IF (LP%CLASS_INDEX /= N_LPC) CYCLE PARTICLE_LOOP_2
      BC => BOUNDARY_COORD(LP%BC_INDEX)

      ! Determine the mass of the PARTICLE/particle, depending on whether the particle has a distinct SURFace type.

      IF (LPC%LIQUID_DROPLET) THEN
         R_DROP = LP%RADIUS
         A_DROP = PI*R_DROP**2
      ELSE
         SF => SURFACE(LPC%SURF_INDEX)
         R_DROP = LP%RADIUS
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

      IF (BC%IOR==0 .OR. LPC%SOLID_PARTICLE) THEN
         LP_B1 => BOUNDARY_PROP1(LP%B1_INDEX)
         DEN_ADD  =    LP%PWT*LP%MASS*LP%RVC
         AREA_ADD =    LP%PWT*A_DROP*LP%RVC
         DROP_DEN(BC%IIG,BC%JJG,BC%KKG)  = DROP_DEN(BC%IIG,BC%JJG,BC%KKG)  + DEN_ADD
         DROP_TMP(BC%IIG,BC%JJG,BC%KKG)  = DROP_TMP(BC%IIG,BC%JJG,BC%KKG)  + DEN_ADD*LP_B1%TMP_F
         DROP_RAD(BC%IIG,BC%JJG,BC%KKG)  = DROP_RAD(BC%IIG,BC%JJG,BC%KKG)  + AREA_ADD*R_DROP
         DROP_AREA(BC%IIG,BC%JJG,BC%KKG) = DROP_AREA(BC%IIG,BC%JJG,BC%KKG) + AREA_ADD
      ENDIF

      ! Compute Mass Per Unit Area (MPUA) for a liquid droplet stuck to a solid surface

      IF (BC%IOR/=0 .AND. (LP%WALL_INDEX>0 .OR. LP%CFACE_INDEX>0)) THEN
         IF (LP%WALL_INDEX>0) THEN
            WC => WALL(LP%WALL_INDEX)
            B1 => BOUNDARY_PROP1(WC%B1_INDEX)
            B2 => BOUNDARY_PROP2(WC%B2_INDEX)
         ELSE
            CFA => CFACE(LP%CFACE_INDEX)
            B1 => BOUNDARY_PROP1(CFA%B1_INDEX)
            B2 => BOUNDARY_PROP2(CFA%B2_INDEX)
         ENDIF
         IF (LPC%LIQUID_DROPLET) THEN
            R_DROP = LP%RADIUS
            FTPR   = FOTHPI * LPC%DENSITY
            M_DROP = FTPR*R_DROP**3
         ELSE
            M_DROP = LP%MASS
         ENDIF
         B2%LP_MPUA(LPC%ARRAY_INDEX) = B2%LP_MPUA(LPC%ARRAY_INDEX) + &
                                       (1._EB-LPC%RUNNING_AVERAGE_FACTOR_WALL)*LP%PWT*M_DROP/B1%AREA
      ENDIF

   ENDDO PARTICLE_LOOP_2

  ! Compute cumulative quantities for PARTICLE "clouds"

   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            DROP_RAD(I,J,K) = DROP_RAD(I,J,K)/(DROP_AREA(I,J,K)+TWO_EPSILON_EB)
            DROP_TMP(I,J,K) = DROP_TMP(I,J,K)/(DROP_DEN(I,J,K) +TWO_EPSILON_EB)
            AVG_DROP_RAD(I,J,K,LPC%ARRAY_INDEX ) = DROP_RAD(I,J,K)
            AVG_DROP_TMP(I,J,K,LPC%ARRAY_INDEX ) = LPC%RUNNING_AVERAGE_FACTOR*AVG_DROP_TMP(I,J,K,LPC%ARRAY_INDEX ) + &
               (1._EB-LPC%RUNNING_AVERAGE_FACTOR)*DROP_TMP(I,J,K)
            IF (LPC%Y_INDEX>0) THEN
               AVG_DROP_TMP(I,J,K,LPC%ARRAY_INDEX ) = MAX(SPECIES(LPC%Y_INDEX)%TMP_MELT,AVG_DROP_TMP(I,J,K,LPC%ARRAY_INDEX ))
            ELSE
               AVG_DROP_TMP(I,J,K,LPC%ARRAY_INDEX ) = MAX(TMPM,AVG_DROP_TMP(I,J,K,LPC%ARRAY_INDEX ))
            ENDIF
            AVG_DROP_DEN(I,J,K,LPC%ARRAY_INDEX ) = LPC%RUNNING_AVERAGE_FACTOR*AVG_DROP_DEN(I,J,K,LPC%ARRAY_INDEX ) + &
               (1._EB-LPC%RUNNING_AVERAGE_FACTOR)*DROP_DEN(I,J,K)
            AVG_DROP_AREA(I,J,K,LPC%ARRAY_INDEX) = LPC%RUNNING_AVERAGE_FACTOR*AVG_DROP_AREA(I,J,K,LPC%ARRAY_INDEX) + &
               (1._EB-LPC%RUNNING_AVERAGE_FACTOR)*DROP_AREA(I,J,K)
            IF (AVG_DROP_DEN(I,J,K,LPC%ARRAY_INDEX )<0.0001_EB .AND. ABS(DROP_DEN(I,J,K))<TWO_EPSILON_EB) &
               AVG_DROP_DEN(I,J,K,LPC%ARRAY_INDEX ) = 0.0_EB
         ENDDO
      ENDDO
   ENDDO

   ! Compute mean liquid droplet temperature (LP_TEMP) on the walls

   IF (LPC%LIQUID_DROPLET) THEN

      ! Initialize work array for storing droplet total area

      DO IW = 1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
         WC => WALL(IW)
         B2 => BOUNDARY_PROP2(WC%B2_INDEX)
         B2%WORK2 = 0._EB  ! drop area
      ENDDO
      DO ICF = INTERNAL_CFACE_CELLS_LB+1,INTERNAL_CFACE_CELLS_LB+N_INTERNAL_CFACE_CELLS
         CFA => CFACE(ICF)
         B2 => BOUNDARY_PROP2(CFA%B2_INDEX)
         B2%WORK2 = 0._EB  ! drop area
      ENDDO

      PARTICLE_LOOP_3: DO IP=1,NLP

         LP => LAGRANGIAN_PARTICLE(IP)
         IF (LP%CLASS_INDEX /= N_LPC) CYCLE PARTICLE_LOOP_3
         BC => BOUNDARY_COORD(LP%BC_INDEX)

         IF (BC%IOR/=0 .AND. (LP%WALL_INDEX>0 .OR. LP%CFACE_INDEX>0)) THEN
            IF (LP%WALL_INDEX>0) THEN
               WC => WALL(LP%WALL_INDEX)
               B2 => BOUNDARY_PROP2(WC%B2_INDEX)
            ELSE
               CFA => CFACE(LP%CFACE_INDEX)
               B2 => BOUNDARY_PROP2(CFA%B2_INDEX)
            ENDIF
            R_DROP = LP%RADIUS
            A_DROP = LP%PWT*PI*R_DROP**2
            B2%WORK2 = B2%WORK2 + A_DROP
         ENDIF

      ENDDO PARTICLE_LOOP_3

      DO IW = 1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
         WC => WALL(IW)
         B2 => BOUNDARY_PROP2(WC%B2_INDEX)
         IF (B2%WORK2 > 0._EB) THEN
            B2%LP_TEMP(LPC%ARRAY_INDEX) = B2%LP_TEMP(LPC%ARRAY_INDEX)/B2%WORK2
            B2%LP_TEMP(LPC%ARRAY_INDEX) = MAX(SPECIES(LPC%Y_INDEX)%TMP_MELT,B2%LP_TEMP(LPC%ARRAY_INDEX))
         ELSE
            B2%LP_TEMP(LPC%ARRAY_INDEX) = TMPA
         ENDIF
      ENDDO
      DO ICF = INTERNAL_CFACE_CELLS_LB+1,INTERNAL_CFACE_CELLS_LB+N_INTERNAL_CFACE_CELLS
         CFA => CFACE(ICF)
         B2 => BOUNDARY_PROP2(CFA%B2_INDEX)
         IF (B2%WORK2 > 0._EB) THEN
            B2%LP_TEMP(LPC%ARRAY_INDEX) = B2%LP_TEMP(LPC%ARRAY_INDEX)/B2%WORK2
            B2%LP_TEMP(LPC%ARRAY_INDEX) = MAX(SPECIES(LPC%Y_INDEX)%TMP_MELT,B2%LP_TEMP(LPC%ARRAY_INDEX))
         ELSE
            B2%LP_TEMP(LPC%ARRAY_INDEX) = TMPA
         ENDIF
      ENDDO

   ENDIF

ENDDO PART_CLASS_SUM_LOOP

END SUBROUTINE PARTICLE_RUNNING_AVERAGES

END SUBROUTINE PARTICLE_MASS_ENERGY_TRANSFER


!> \brief Add PARTICLE momentum as a force term in momentum equation
!> \param DT Current time step (s)
!> \param NM Current mesh number

SUBROUTINE PARTICLE_MOMENTUM_TRANSFER(DT,NM)

USE CC_SCALARS, ONLY: CUTFACE_VELOCITIES
INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: DT
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW
REAL(EB) :: RDT,UODT,VODT,WODT
INTEGER :: I,J,K

IF (MESHES(NM)%NLP==0) RETURN

CALL POINT_TO_MESH(NM)

RDT = 1._EB/DT

IF (PREDICTOR) THEN
   UU => U
   VV => V
   WW => W
ELSE
   UU => US
   VV => VS
   WW => WS
ENDIF

IF (CC_IBM) CALL CUTFACE_VELOCITIES(NM,UU,VV,WW,CUTFACES=.TRUE.)

! Add summed particle accelerations to the momentum equation. Limit the value to plus/minus abs(u)/dt to prevent a sudden
! change in gas direction.

DO K=0,KBAR
   DO J=0,JBAR
      DO I=0,IBAR
         UODT = ABS(UU(I,J,K)*RDT)
         VODT = ABS(VV(I,J,K)*RDT)
         WODT = ABS(WW(I,J,K)*RDT)
         FVX(I,J,K) = FVX(I,J,K) + MIN(UODT,MAX(-UODT,FVX_D(I,J,K)))
         FVY(I,J,K) = FVY(I,J,K) + MIN(VODT,MAX(-VODT,FVY_D(I,J,K)))
         FVZ(I,J,K) = FVZ(I,J,K) + MIN(WODT,MAX(-WODT,FVZ_D(I,J,K)))
      ENDDO
   ENDDO
ENDDO

IF (CC_IBM) CALL CUTFACE_VELOCITIES(NM,UU,VV,WW,CUTFACES=.FALSE.)

END SUBROUTINE PARTICLE_MOMENTUM_TRANSFER


!> \brief Remove PARTICLEs that have left the current mesh (NM) or are no longer to be tracked
!> \param T Current time (s)
!> \param NM Current mesh number

SUBROUTINE REMOVE_PARTICLES(T,NM)

INTEGER, INTENT(IN) :: NM
INTEGER :: IKILL,IP
REAL(EB), INTENT(IN) :: T
TYPE (LAGRANGIAN_PARTICLE_TYPE), POINTER :: LP
TYPE (LAGRANGIAN_PARTICLE_CLASS_TYPE), POINTER :: LPC
TYPE (SURFACE_TYPE), POINTER :: SF
TYPE (MESH_TYPE), POINTER :: M,M2
TYPE (BOUNDARY_COORD_TYPE), POINTER :: BC
TYPE (BOUNDARY_PROP1_TYPE), POINTER :: B1
TYPE (STORAGE_TYPE), POINTER :: OS

M => MESHES(NM)
IKILL = 0

PARTICLE_LOOP: DO IP=1,NLP

   WEED_LOOP: DO

      IF (IP>NLP-IKILL) EXIT PARTICLE_LOOP

      LP  => M%LAGRANGIAN_PARTICLE(IP)
      LPC => LAGRANGIAN_PARTICLE_CLASS(LP%CLASS_INDEX)
      BC  => M%BOUNDARY_COORD(LP%BC_INDEX)
      SF  => SURFACE(LPC%SURF_INDEX)

      ! Ignore HVAC particles

      IF (LPC%DUCT_PARTICLE) CYCLE PARTICLE_LOOP

      ! Remove particles that are too small

      IF (LPC%SOLID_PARTICLE .AND. SF%THERMAL_BC_INDEX==THERMALLY_THICK) THEN
         B1 => M%BOUNDARY_PROP1(LP%B1_INDEX)
         IF (B1%BURNAWAY) THEN
            CALL ADD_TO_PARTICLE_SEND_BUFFER
            CYCLE WEED_LOOP
         ENDIF
      ENDIF

      IF (LPC%LIQUID_DROPLET) THEN
         IF (LP%RADIUS<=LPC%KILL_RADIUS) THEN
            CALL ADD_TO_PARTICLE_SEND_BUFFER
            CYCLE WEED_LOOP
         ENDIF
      ENDIF

      ! Remove particles that are too old

      IF (T-LP%T_INSERT > LPC%LIFETIME) THEN
         CALL ADD_TO_PARTICLE_SEND_BUFFER
         CYCLE WEED_LOOP
      ENDIF

      ! Remove particles that have left the active mesh

      IF (BC%X>XS .AND. BC%X<XF .AND. &
          BC%Y>YS .AND. BC%Y<YF .AND. &
          BC%Z>ZS .AND. BC%Z<ZF) CYCLE PARTICLE_LOOP

      ! Replace all other particles

      CALL ADD_TO_PARTICLE_SEND_BUFFER

   ENDDO WEED_LOOP

ENDDO PARTICLE_LOOP

NLP = NLP - IKILL

CONTAINS


!> \brief Determine if the given particle is now in another mesh, and if so, pack it into the SEND_BUFFER

SUBROUTINE ADD_TO_PARTICLE_SEND_BUFFER

USE MEMORY_FUNCTIONS, ONLY: REALLOCATE_STORAGE_ARRAYS,NULLIFY_PARTICLE,PACK_PARTICLE
INTEGER :: MM,NOM

! Check to see if the particle has entered another mesh (NOM)

NOM = 0

SEARCH_LOOP: DO MM=1,N_NEIGHBORING_MESHES
   M2=>MESHES(NEIGHBORING_MESH(MM))
   IF (BC%X>=M2%XS .AND. BC%X<=M2%XF .AND.  &
       BC%Y>=M2%YS .AND. BC%Y<=M2%YF .AND.  &
       BC%Z>=M2%ZS .AND. BC%Z<=M2%ZF) THEN
      NOM = NEIGHBORING_MESH(MM)
      EXIT SEARCH_LOOP
   ENDIF
ENDDO SEARCH_LOOP

! Don't remove particles that just cycle in a periodic domain

IF (PERIODIC_DOMAIN_X) THEN
   IF (BC%X>=XF_MAX) THEN
      BC%X = BC%X - (XF_MAX- XS_MIN)
      RETURN
   ENDIF
   IF (BC%X<=XS_MIN) THEN
      BC%X = BC%X + (XF_MAX- XS_MIN)
      RETURN
   ENDIF
ENDIF

IF (PERIODIC_DOMAIN_Y) THEN
   IF (BC%Y>=YF_MAX) THEN
      BC%Y = BC%Y - (YF_MAX- YS_MIN)
      RETURN
   ENDIF
   IF (BC%Y<=YS_MIN) THEN
      BC%Y = BC%Y + (YF_MAX- YS_MIN)
      RETURN
   ENDIF
ENDIF

IF (PERIODIC_DOMAIN_Z) THEN
   IF (BC%Z>=ZF_MAX) THEN
      BC%Z = BC%Z - (ZF_MAX- ZS_MIN)
      RETURN
   ENDIF
   IF (BC%Z<=ZS_MIN) THEN
      BC%Z = BC%Z + (ZF_MAX- ZS_MIN)
      RETURN
   ENDIF
ENDIF

! For particles that have entered another mesh (NOM), add them to the "orphanage"

IF (NOM/=0 .AND. NOM/=NM) THEN

   OS => OMESH(NOM)%PARTICLE_SEND_BUFFER
   IF (OS%N_REALS+LPC%N_REALS>OS%N_REALS_DIM) THEN
      CALL REALLOCATE_STORAGE_ARRAYS(OS,N_REALS=OS%N_REALS_DIM,N_REALS_NEW=OS%N_REALS_DIM+10*LPC%N_REALS)
      OS%N_REALS_DIM = OS%N_REALS_DIM + 10*LPC%N_REALS
   ENDIF
   IF (OS%N_INTEGERS+LPC%N_INTEGERS>OS%N_INTEGERS_DIM) THEN
      CALL REALLOCATE_STORAGE_ARRAYS(OS,N_INTEGERS=OS%N_INTEGERS_DIM,N_INTEGERS_NEW=OS%N_INTEGERS_DIM+10*LPC%N_INTEGERS)
      OS%N_INTEGERS_DIM = OS%N_INTEGERS_DIM + 10*LPC%N_INTEGERS
   ENDIF
   IF (OS%N_LOGICALS+LPC%N_LOGICALS>OS%N_LOGICALS_DIM) THEN
      CALL REALLOCATE_STORAGE_ARRAYS(OS,N_LOGICALS=OS%N_LOGICALS_DIM,N_LOGICALS_NEW=OS%N_LOGICALS_DIM+10*LPC%N_LOGICALS)
      OS%N_LOGICALS_DIM = OS%N_LOGICALS_DIM + 10*LPC%N_LOGICALS
   ENDIF

   OS%N_ITEMS = OS%N_ITEMS + 1
   CALL PACK_PARTICLE(NM,OS,LP,LP%CLASS_INDEX,OS%N_REALS,OS%N_INTEGERS,OS%N_LOGICALS,UNPACK_IT=.FALSE.,COUNT_ONLY=.FALSE.,&
                      CHECK_BOUNDS=.FALSE.)

ENDIF

! Copy last particle into the slot of the vacated particle. Zero out the last particle.

CALL NULLIFY_PARTICLE(NM,IP)

LAGRANGIAN_PARTICLE(IP) = LAGRANGIAN_PARTICLE(NLP-IKILL)

LAGRANGIAN_PARTICLE(NLP-IKILL)%BC_INDEX = 0
LAGRANGIAN_PARTICLE(NLP-IKILL)%OD_INDEX = 0
LAGRANGIAN_PARTICLE(NLP-IKILL)%B1_INDEX = 0
LAGRANGIAN_PARTICLE(NLP-IKILL)%B2_INDEX = 0
LAGRANGIAN_PARTICLE(NLP-IKILL)%BR_INDEX = 0

IKILL = IKILL + 1

END SUBROUTINE ADD_TO_PARTICLE_SEND_BUFFER

END SUBROUTINE REMOVE_PARTICLES


!> \brief Remove the oldest particle of class LPC_INDEX and move particle NLP into its place.
!> \param NM Current mesh number
!> \param LPC_INDEX Particle class index
!> \param NLP Index of the newest Lagrangian particle that is to replace the oldest
!> \param LP_INDEX Index of the oldest particle that is replaced by the newest

SUBROUTINE REMOVE_OLDEST_PARTICLE(NM,LPC_INDEX,NLP,LP_INDEX)

USE MEMORY_FUNCTIONS, ONLY: NULLIFY_PARTICLE
INTEGER, INTENT(IN) :: NM,LPC_INDEX,NLP
INTEGER, INTENT(OUT) :: LP_INDEX
INTEGER :: IP,TAG_MIN
TYPE(LAGRANGIAN_PARTICLE_TYPE), POINTER :: LP

! Look for the oldest particle of this class

TAG_MIN = HUGE(NLP)

DO IP=1,MESHES(NM)%NLP
   LP => MESHES(NM)%LAGRANGIAN_PARTICLE(IP)
   IF (LP%CLASS_INDEX/=LPC_INDEX) CYCLE
   IF (LP%TAG>0 .AND. LP%TAG<TAG_MIN) THEN
       TAG_MIN = LP%TAG
       LP_INDEX = IP
       IF (TAG_MIN<=MESHES(NM)%PARTICLE_LAST(LPC_INDEX)) EXIT
   ENDIF
ENDDO

! Move the last LAGRANGIAN_PARTICLE into the slot left open after removing the oldest LAGRANGIAN_PARTICLE

CALL NULLIFY_PARTICLE(NM,LP_INDEX)

MESHES(NM)%LAGRANGIAN_PARTICLE(LP_INDEX) = MESHES(NM)%LAGRANGIAN_PARTICLE(NLP)
MESHES(NM)%PARTICLE_LAST(LPC_INDEX) = TAG_MIN + NMESHES

LAGRANGIAN_PARTICLE(NLP)%BC_INDEX = 0
LAGRANGIAN_PARTICLE(NLP)%OD_INDEX = 0
LAGRANGIAN_PARTICLE(NLP)%B1_INDEX = 0
LAGRANGIAN_PARTICLE(NLP)%B2_INDEX = 0
LAGRANGIAN_PARTICLE(NLP)%BR_INDEX = 0

END SUBROUTINE REMOVE_OLDEST_PARTICLE

!> \brief Return reciprocal of particle-containing cell volume
!> \param NM Current mesh number
!> \param IIG gas phase cell index
!> \param JJG gas phase cell index
!> \param KKG gas phase cell index
!> \param RVC reciprocal of cell volume

SUBROUTINE GET_RVC(NM,IIG,JJG,KKG,RVC)

INTEGER, INTENT(IN) :: IIG,JJG,KKG,NM
INTEGER :: ICC
REAL(EB), INTENT(OUT) :: RVC

RVC = MESHES(NM)%RDX(IIG)*MESHES(NM)%RRN(IIG)*MESHES(NM)%RDY(JJG)*MESHES(NM)%RDZ(KKG)
IF (CC_IBM) THEN
   ICC = MESHES(NM)%CCVAR(IIG,JJG,KKG,CC_IDCC)
   IF (ICC>0) RVC = 1/SUM(MESHES(NM)%CUT_CELL(ICC)%VOLUME(:))
ENDIF

END SUBROUTINE GET_RVC

END MODULE PART
