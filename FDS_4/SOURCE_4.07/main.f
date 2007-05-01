      PROGRAM FDS  
C
C Fire Dynamics Simulator, Main Program, Single CPU version
C
      USE PREC
      USE VARS
      USE CONS
      USE TRAN
      USE DUMP
      USE READ
      USE INIT
      USE DIVG
      USE PRES
      USE MASS
      USE SPRK
      USE VELO
      USE PYRO
      USE RAD
CEVAC:      USE EVAC
C
      IMPLICIT NONE
C
C     Miscellaneous declarations
C
      LOGICAL  :: EX,DIAGNOSTICS
      INTEGER  :: LO10,NM,IZERO,STOP_CODE=0
      REAL(EB) :: T_MAX,T_MIN
      REAL(EB), ALLOCATABLE, DIMENSION(:) ::  T,DT_SYNC
      INTEGER, ALLOCATABLE, DIMENSION(:) ::  ISTOP
      LOGICAL, ALLOCATABLE, DIMENSION(:) ::  ACTIVE_MESH
      INTEGER NOM,IMIN,IMAX,JMIN,JMAX,KMIN,KMAX,IW,N,NN
      INTEGER, PARAMETER :: N_DROP_ADOPT_MAX=100
      TYPE (MESH_TYPE), POINTER :: M,M4
      TYPE (OMESH_TYPE), POINTER :: M2,M3
      TYPE (HEAT_DETECTOR_TYPE), POINTER :: HD
      TYPE (SMOKE_DETECTOR_TYPE),POINTER :: SD
      TYPE (OBSTRUCTION_TYPE),POINTER :: OB
      TYPE (VENTS_TYPE),POINTER :: VT
      TYPE (SPRINKLER_HEAD_TYPE),POINTER :: SH
      TYPE (SPRINKLER_TYPE),POINTER :: S
C
C     Assign a compilation date (All Nodes)
C
      COMPILE_DATE = 'March 10, 2006'
      VERSION      = 4.07
C
C     Open the input file (All Nodes)
C
C     CALL OPEN_INPUT_FILE
C
C     Read input from CHID.data file (All Nodes)
C
      CALL READ_DATA
CEVAC:C
CEVAC:C     Read input for EVACUATION routine
CEVAC:C
CEVAC:      IF (ANY(EVACUATION_GRID)) CALL READ_EVAC
C
C     Open and write to Smokeview file (Master Node Only)
C
      CALL WRITE_SMOKEVIEW_FILE
C
C     Stop all the processes if this is just a set-up run
C
      IF (SET_UP) CALL SHUTDOWN('Stop FDS, Set-up only')
C
C     Set up Time array (All Nodes)
C
      ALLOCATE(ACTIVE_MESH(NMESHES),STAT=IZERO)
      CALL ChkMemErr('MAIN','ACTIVE_MESH',IZERO)
      ALLOCATE(T(NMESHES),STAT=IZERO)
      CALL ChkMemErr('MAIN','T',IZERO)
      ALLOCATE(DT_SYNC(NMESHES),STAT=IZERO)
      CALL ChkMemErr('MAIN','DT_SYNC',IZERO)
      ALLOCATE(ISTOP(NMESHES),STAT=IZERO)
      CALL ChkMemErr('MAIN','ISTOP',IZERO)
      T     = 0.
      ISTOP = 0
C
C     Initialize global parameters (All Nodes)
C
      CALL INITIALIZE_GLOBAL_VARIABLES
C
C     Initialize radiation (All Nodes)
C
      IF (RADIATION) CALL INIT_RADIATION
C
C     Allocate and initialize mesh-specific variables
C
      DO NM=1,NMESHES
      CALL INITIALIZE_MESH_VARIABLES(NM)
      ENDDO
C
C     Allocate and initialize mesh variable exchange arrays
C
      DO NM=1,NMESHES
      CALL INITIALIZE_MESH_EXCHANGE(NM)
      ENDDO
C
      I_MIN = TRANSPOSE(I_MIN)
      I_MAX = TRANSPOSE(I_MAX)
      J_MIN = TRANSPOSE(J_MIN)
      J_MAX = TRANSPOSE(J_MAX)
      K_MIN = TRANSPOSE(K_MIN)
      K_MAX = TRANSPOSE(K_MAX)
      NIC   = TRANSPOSE(NIC)
C
      DO NM=1,NMESHES
      CALL DOUBLE_CHECK(NM)
      ENDDO
C
C     Initialize miscellaneous routines
C
      IF (PAPERMODEL) CALL INIT_SOLID_PHASE_REACTIONS
C
C     Potentially read data from a previous calculation 
C
      DO NM=1,NMESHES
      IF (RESTART) CALL READ_CORE(T(NM),NM)
      ENDDO
C
C     Initialize output files containing global data (Master Node Only)
C
      CALL INITIALIZE_GLOBAL_DUMPS
C
C     Initialize output files that are mesh-specific
C
      DO NM=1,NMESHES
      CALL INITIALIZE_MESH_DUMPS(NM)
      CALL INITIALIZE_DROPLETS(NM)
      CALL INITIALIZE_TREES(NM)
CEVAC:      IF (ANY(EVACUATION_GRID)) CALL INITIALIZE_EVACUATION(NM)
      ENDDO
C
C     Write out character strings to .smv file
C
      CALL WRITE_STRINGS
C
C     Initialize Mesh Exchange Arrays (All Nodes)
C
      CALL MESH_EXCHANGE(0)
C
C ********************************************************************
C                      MAIN TIMESTEPPING LOOP
C ********************************************************************
CEVAC:C
CEVAC:      IF (ANY(EVACUATION_GRID)) THEN
CEVAC:        DO NM=1,NMESHES
CEVAC:          IF (EVACUATION_GRID(NM))
CEVAC:     .         PARCLK(NM) = -EVAC_DT_FLOWFIELD*EVAC_TIME_ITERATIONS
CEVAC:        ENDDO
CEVAC:        ICYC = -EVAC_TIME_ITERATIONS
CEVAC:      ENDIF
      MAIN_LOOP: DO  
C
      ICYC  = ICYC + 1 
      IF (ANY(EVACUATION_GRID) .AND. ICYC.EQ.0) ICYC = 1
C
C     Check for program stops
C
      INQUIRE(FILE=TRIM(CHID)//'.stop',EXIST=EX)
      IF (EX) ISTOP = 2
C
C     Figure out fastest and slowest meshes
C
      T_MAX = 0.
      T_MIN = 1000000.
      DO NM=1,NMESHES
      T_MIN = MIN(T(NM),T_MIN)
      T_MAX = MAX(T(NM),T_MAX)
      IF (ISTOP(NM).GT.0) STOP_CODE = ISTOP(NM)
      ENDDO
C
      IF (SYNCHRONIZE) THEN
C
      DT_SYNC(1:NMESHES) = MESH(1:NMESHES)%DTNEXT
      DO NM=1,NMESHES
      IF (SYNC_TIME_STEP(NM)) THEN
         MESH(NM)%DTNEXT = MINVAL(DT_SYNC,MASK=SYNC_TIME_STEP)
         T(NM) = MINVAL(T,MASK=SYNC_TIME_STEP)
         ACTIVE_MESH(NM) = .TRUE.
         ELSE
         ACTIVE_MESH(NM) = .FALSE.
         IF (T(NM)+MESH(NM)%DTNEXT.LE.T_MAX) ACTIVE_MESH(NM) = .TRUE.
         IF (STOP_CODE.GT.0) ACTIVE_MESH(NM) = .TRUE.
      ENDIF
      ENDDO
C
      ELSE
C
      ACTIVE_MESH = .FALSE.
      DO NM=1,NMESHES
      IF (T(NM)+MESH(NM)%DTNEXT .LE. T_MAX) ACTIVE_MESH(NM) = .TRUE.
      IF (STOP_CODE.GT.0) ACTIVE_MESH(NM) = .TRUE.
      ENDDO
C
      ENDIF
C
      DIAGNOSTICS = .FALSE.
      LO10 = LOG10(REAL(ABS(ICYC),EB))
      IF (MOD(ICYC,10**LO10).EQ.0 .OR. MOD(ICYC,100).EQ.0 .OR.
     .    T_MIN.GE.TWFIN .OR. STOP_CODE.GT.0) DIAGNOSTICS = .TRUE.
C
C     If no meshes are due to be updated, update them all
C
      IF (ALL(.NOT.ACTIVE_MESH)) ACTIVE_MESH = .TRUE.
CEVAC:      IF (ANY(EVACUATION_GRID).AND.(ICYC.LE.0)) ACTIVE_MESH = .FALSE.
CEVAC:C
CEVAC:C     Do not do EVACuation flows if past the max iteration criteria
CEVAC:C
CEVAC:      IF (ANY(EVACUATION_GRID)) THEN
CEVAC:        EVAC_DT = 1000000.
CEVAC:        IF (ICYC .GT. 0) THEN
CEVAC:          DO NM=1,NMESHES
CEVAC:            IF (.NOT.EVACUATION_ONLY(NM))
CEVAC:     .           EVAC_DT = MIN(EVAC_DT,MESH(NM)%DTNEXT)
CEVAC:          ENDDO
CEVAC:        ENDIF
CEVAC:        DO NM=1,NMESHES
CEVAC:          IF (EVACUATION_ONLY(NM)) THEN
CEVAC:            IF (ICYC .GT. 0) THEN
CEVAC:              ACTIVE_MESH(NM) = .FALSE.
CEVAC:              EVAC_DT = MIN(EVAC_DT, EVAC_DT_STEADY_STATE)
CEVAC:              MESH(NM)%DT     = EVAC_DT
CEVAC:              T(NM)           = T(NM) + MESH(NM)%DT
CEVAC:              IF (EVACUATION_GRID(NM) ) THEN
CEVAC:                CALL EVACUATE_HUMANS(T(NM),NM)
CEVAC:                IF (T(NM).GE.PARCLK(NM)) THEN
CEVAC:                  CALL DUMP_EVAC(PARCLK(NM),NM)
CEVAC:                  DO
CEVAC:                    PARCLK(NM) = PARCLK(NM) + WPAR
CEVAC:                    IF (PARCLK(NM).GE.T(NM)) EXIT
CEVAC:                  ENDDO
CEVAC:                ENDIF
CEVAC:              ENDIF
CEVAC:            ELSE
CEVAC:              ACTIVE_MESH(NM) = .TRUE.
CEVAC:              EVAC_DT = MIN(EVAC_DT, EVAC_DT_FLOWFIELD)
CEVAC:            ENDIF
CEVAC:          ENDIF
CEVAC:        ENDDO
CEVAC:        DO NM=1,NMESHES
CEVAC:          IF (EVACUATION_ONLY(NM)) MESH(NM)%DTNEXT = EVAC_DT
CEVAC:        ENDDO
CEVAC:      ENDIF
C
C=====================================================================
C
      PREDICTOR = .TRUE.
      CORRECTOR = .FALSE.
C
      MULTIBLOCK1: DO NM=1,NMESHES
C
      IF (.NOT.ACTIVE_MESH(NM)) CYCLE MULTIBLOCK1
C
      MESH(NM)%DT = MESH(NM)%DTNEXT
      NTCYC(NM)   = NTCYC(NM) + 1
C
C Insert tracer particles 
C
      IF (T(NM).GE.PINCLK(NM)) CALL INSERT_PARTICLES(T(NM),NM)
C
C Insert sprinkler droplets
C
      IF (T(NM).GE.SPINCLK(NM)) CALL INSERT_DROPLETS(T(NM),NM)
C
C Compute momentum equation spatial differences at the present time step
C
      IF (DNS .AND. ISOTHERMAL .AND. 
     .             (INCOMPRESSIBLE .OR. NSPEC.EQ.0)) THEN
      CALL VELOCITY_FLUX_ISOTHERMAL(T(NM),NM)
      ELSE
      IF (.NOT.CYLINDRICAL) CALL VELOCITY_FLUX(T(NM),NM)
      IF (     CYLINDRICAL) CALL VELOCITY_FLUX_CYLINDRICAL(T(NM),NM)
      ENDIF
C
C Track sprinkler droplets
C
      IF (MESH(NM)%NLP.GT.0) CALL TRACK_DROPLETS(T(NM),NM)
C
C Compute spatial differences of density and species equations
C
      IF (.NOT.ISOTHERMAL .OR. NSPEC.GT.0) CALL MASS_FLUX(NM)
C
C Loop to predict variables at next time step subject to CFL constraint
C
      PREDICTOR_LOOP: DO
C
C Predict density, mass fractions and temperatures at next time step  
C
      IF (.NOT.ISOTHERMAL .OR. NSPEC.GT.0) CALL DENSITY(T(NM),NM)
C
C Compute the divergence and the time derivative of divergence
C
      CALL DIVERGENCE(T(NM),NM)
C
C Solve the Poisson equation for the pressure
C
      CALL PRESSURE_SOLVER(NM)
CEVAC:C
CEVAC:      IF (ANY(EVACUATION_GRID) .AND. EVACUATION_ONLY(NM)) THEN
CEVAC:         PRESSURE_ITERATION_LOOP: DO N=1,EVAC_PRESSURE_ITERATIONS
CEVAC:            CALL NO_FLUX(T(NM))
CEVAC:            CALL PRESSURE_SOLVER(NM)
CEVAC:         ENDDO PRESSURE_ITERATION_LOOP
CEVAC:      ENDIF
C
C Compute velocity estimates at next time step
C
      CALL VELOCITY_PREDICTOR(T(NM),NM,ISTOP(NM))
C
C If time step too small, stop; if time step changed go back
C
      IF (.NOT.MESH(NM)%NEW_TIME_STEP .OR. ISTOP(NM).EQ.1) 
     .   EXIT PREDICTOR_LOOP
C
      ENDDO PREDICTOR_LOOP
C
C Estimate sprinkler link temperatures at the next time step
C
      IF (NSPR.GT.0) CALL CHECK_SPRINKLERS(T(NM),NM)
      IF (NHD .GT.0) CALL CHECK_HEAT_DETECTORS(T(NM),NM)
C
C
C Call smoke detector modeling at next time step.
C
      IF (NSD .GT.0) CALL SMOKE_DETECTORS(T(NM),NM)
C
C Advance the time, start the corrector phase
C
      T(NM) = T(NM) + MESH(NM)%DT
C
      ENDDO MULTIBLOCK1
C
C======================================================================
C
C Exchange information among meshes
C
      CALL MESH_EXCHANGE(1)
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      CORRECTOR = .TRUE.
      PREDICTOR = .FALSE.
      MULTIBLOCK2: DO NM=1,NMESHES
C
      IF (.NOT.ACTIVE_MESH(NM)) CYCLE MULTIBLOCK2
C
C Compute the momentum equation spatial differences
C
      IF (DNS .AND. ISOTHERMAL .AND. 
     .             (INCOMPRESSIBLE .OR. NSPEC.EQ.0)) THEN
      CALL VELOCITY_FLUX_ISOTHERMAL(T(NM),NM)
      ELSE
      IF (.NOT.CYLINDRICAL) CALL VELOCITY_FLUX(T(NM),NM)
      IF (     CYLINDRICAL) CALL VELOCITY_FLUX_CYLINDRICAL(T(NM),NM)
      ENDIF
C
C Correct density, mass fractions and temperature
C
      IF (.NOT.ISOTHERMAL .OR. NSPEC.GT.0) THEN
      CALL MASS_FLUX(NM)
      CALL DENSITY(T(NM),NM)
      ENDIF
C
C Correct sprinkler droplet locations, distribute heat from droplets
C
      IF (MESH(NM)%NLP.GT.0) CALL TRACK_DROPLETS(T(NM),NM)
C
C Compute the divergence and time derivative of divergence
C
      CALL DIVERGENCE(T(NM),NM)
C
C Solve the Poisson equation for pressure
C
      CALL PRESSURE_SOLVER(NM)
CEVAC:C
CEVAC:      IF (ANY(EVACUATION_GRID) .AND. EVACUATION_ONLY(NM)) THEN
CEVAC:         PRESSURE_ITERATION_LOOP2: DO N=1,EVAC_PRESSURE_ITERATIONS
CEVAC:            CALL NO_FLUX(T(NM))
CEVAC:            CALL PRESSURE_SOLVER(NM)
CEVAC:         ENDDO PRESSURE_ITERATION_LOOP2
CEVAC:      ENDIF
C
C Open or close doors
C
      CALL OPEN_AND_CLOSE(T(NM),NM)
C
C Update the velocities 
C
      CALL VELOCITY_CORRECTOR(T(NM),NM)
C
C Check to see if the sprinklers went off
C
      IF (NSPR.GT.0) CALL CHECK_SPRINKLERS(T(NM),NM)
      IF (NHD .GT.0) CALL CHECK_HEAT_DETECTORS(T(NM),NM)

C
C Call smoke detector modeling at next time step.
C
      IF (NSD .GT.0) CALL SMOKE_DETECTORS(T(NM),NM)
C
C Mesh-specific data updates and dumps
C
      CALL UPDATE_HRR(NM)
      IF (NSPEC.GT.0) CALL UPDATE_MASS(NM)
      IF (NTC.GT.0) CALL UPDATE_TC(NM)
CEVAC:C
CEVAC:C EVACuation calculation and dumps for EVACUATION_ONLY
CEVAC:C EVACcuation here, if flow field initialization phase
CEVAC:C
CEVAC:      EVACUATION_DUMP: IF (ANY(EVACUATION_GRID) .AND.
CEVAC:     .     EVACUATION_ONLY(NM)) THEN
CEVAC:        IF (EVACUATION_GRID(NM) ) THEN
CEVAC:          CALL EVACUATE_HUMANS(T(NM)-EVAC_DT_FLOWFIELD*
CEVAC:     .         EVAC_TIME_ITERATIONS,NM)
CEVAC:          IF (T(NM)-EVAC_DT_FLOWFIELD*
CEVAC:     .         EVAC_TIME_ITERATIONS.GE.PARCLK(NM)) THEN
CEVAC:            CALL DUMP_EVAC(PARCLK(NM),NM)
CEVAC:            DO
CEVAC:              PARCLK(NM) = PARCLK(NM) + WPAR
CEVAC:              IF (PARCLK(NM).GE.T(NM)-EVAC_DT_FLOWFIELD*
CEVAC:     .         EVAC_TIME_ITERATIONS) EXIT
CEVAC:            ENDDO
CEVAC:          ENDIF
CEVAC:        ENDIF
CEVAC:C     Dump the EVAC flowfieds for all EVAC meshes
CEVAC:        CALL DUMP_SF(T(NM)-EVAC_DT_FLOWFIELD*
CEVAC:     .         EVAC_TIME_ITERATIONS,NM,0)
CEVAC:C
CEVAC:      ELSE
C
      IF (T(NM).GE.PARCLK(NM).AND.DROPLET_FILE) THEN
         CALL DUMP_PART(PARCLK(NM),NM)
         DO
         PARCLK(NM) = PARCLK(NM) + WPAR
         IF (PARCLK(NM).GE.T(NM)) EXIT
         ENDDO
         ENDIF
C
      IF (T(NM).GE.ISOCLK(NM)) THEN
         CALL DUMP_ISO(ISOCLK(NM),NM)
         DO
         ISOCLK(NM) = ISOCLK(NM) + DTIF
         IF (ISOCLK(NM).GE.T(NM)) EXIT
         ENDDO
         ENDIF
C
      IF (T(NM).GE.SFCLK(NM)) THEN
         CALL DUMP_SF(SFCLK(NM),NM,0)
         IF (SMOKE3D) CALL DUMP_SMOKE3D(SFCLK(NM),NM)
         DO
         SFCLK(NM) = SFCLK(NM) + DTSF
         IF (SFCLK(NM).GE.T(NM)) EXIT
         ENDDO
         ENDIF
C
      IF (T(NM).GE.BFCLK(NM)) THEN
         CALL DUMP_BF(BFCLK(NM),NM)
         DO
         BFCLK(NM) = BFCLK(NM) + DTBF
         IF (BFCLK(NM).GE.T(NM)) EXIT
         ENDDO
         ENDIF
C
      IF (T(NM).GE.PLTCLK(NM) .OR. STOP_CODE.GT.0) THEN
         CALL DUMP_SF(T(NM),NM,1)
         DO
         PLTCLK(NM) = PLTCLK(NM) + WPLT
         IF (PLTCLK(NM).GE.T(NM)) EXIT
         ENDDO
         ENDIF
C
      IF (T(NM).GE.CORCLK(NM) .OR. STOP_CODE.EQ.2) THEN
         CALL DUMP_CORE(T(NM),NM)
         CORCLK(NM) = CORCLK(NM) + DTCORE
         ENDIF
C
CEVAC:      ENDIF EVACUATION_DUMP
C
      IF (DIAGNOSTICS) CALL CHECK_DIVERGENCE(NM)
C
      ENDDO MULTIBLOCK2
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     Exchange information among meshes
C
      CALL MESH_EXCHANGE(2)
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     Write character strings out to the .smv file
C
      CALL WRITE_STRINGS
C
C     Exchange info for diagnostic print out
C
      IF (DIAGNOSTICS) CALL EXCHANGE_DIAGNOSTICS
C
C     Dump out HRR info
C
      IF (T_MIN.GE.HRRCLK .AND.
     .     MINVAL(HRR_COUNT,MASK=.NOT.EVACUATION_ONLY).GT.0.) THEN
      CALL DUMP_HRR(T_MIN)
      HRRCLK = HRRCLK + DTHRR
      HRR_SUM   = 0.
      RHRR_SUM  = 0.
      CHRR_SUM  = 0.
      FHRR_SUM  = 0.
      MLR_SUM   = 0.
      HRR_COUNT = 0.
      ENDIF
CEVAC:C
CEVAC:C     Dump out Evac info
CEVAC:C
CEVAC:      IF (T_MIN.GE.EVACCLK .AND. ANY(EVACUATION_GRID)) THEN
CEVAC:        CALL DUMP_EVAC_CSV(T_MIN)
CEVAC:        EVACCLK = EVACCLK + DTHRR
CEVAC:      ENDIF
C
C     Dump out Mass info
C
      IF (T_MIN.GE.MINTCLK .AND.
     .     MINVAL(MINT_COUNT,MASK=.NOT.EVACUATION_ONLY).GT.0.) THEN
      CALL DUMP_MASS(T_MIN)
      MINTCLK    = MINTCLK + DTMINT
      MINT_SUM   = 0.
      MINT_COUNT = 0.
      ENDIF
C
C     Dump out TC data
C
      IF (T_MIN.GE.TCCLK) THEN
      DO N=1,NTC
      IF (THERMOCOUPLE(N)%COUNT.EQ.0) GOTO 110
      ENDDO
      CALL DUMP_TC(T_MIN)
      TCCLK = TCCLK + DTTC
      THERMOCOUPLE(1:NTC)%VALUE = 0.
      THERMOCOUPLE(1:NTC)%COUNT = 0
      ENDIF
C
  110 CONTINUE
C
C     Dump out Sprinkler and Heat Detector link temps
C
      IF (T_MIN.GT.SPRKCLK) THEN
      CALL DUMP_DETECTORS(T_MIN,1)
      SPRKCLK = SPRKCLK + DT_SPRK
      ENDIF
C
      IF (T_MIN.GT.HEATCLK) THEN
      CALL DUMP_DETECTORS(T_MIN,2)
      HEATCLK = HEATCLK + DT_HEAT
      ENDIF
C
C     Dump out Smoke detector responses
C
      IF (T_MIN.GT.SMOKECLK) THEN
      CALL DUMP_DETECTORS(T_MIN,3)
      SMOKECLK = SMOKECLK + DT_SMOKE
      ENDIF
C
C     Look for actions associated with heat detector activation
C
      HEAT_DETECTOR_LOOP: DO N=1,NHD
C
      HD=>HEAT_DETECTOR(N)
C
      HEAT_DETECTOR_ACT: IF (HD%TMP_L.GE.HD%TMP_ACT) THEN
C
         MESH_LOOP: DO NM=1,NMESHES
            M => MESH(NM)
            DO NN=1,M%NB
            OB=>M%OBSTRUCTION(NN)
            IF (OB%HEAT_INDEX_REMOVE.EQ.N .OR. 
     .          OB%HEAT_INDEX_REMOVE.GT.NHD) THEN
               OB%T_REMOVE = T(NM)
               OB%HEAT_INDEX_REMOVE = 0
               ENDIF
            IF (OB%HEAT_INDEX_CREATE.EQ.N .OR. 
     .          OB%HEAT_INDEX_CREATE.GT.NHD) THEN
               OB%T_CREATE = T(NM)
               OB%HEAT_INDEX_CREATE = 0
               ENDIF
            ENDDO
            DO NN=1,M%NV
            VT=>M%VENTS(NN)
            IF (VT%HEAT_INDEX_ACTIVATE.EQ.N .OR. 
     .          VT%HEAT_INDEX_ACTIVATE.GT.NHD) THEN
               VT%T_OPEN     = T(NM)
               VT%HEAT_INDEX_ACTIVATE = 0
               ENDIF
            IF (VT%HEAT_INDEX_DEACTIVATE.EQ.N .OR. 
     .          VT%HEAT_INDEX_DEACTIVATE.GT.NHD) THEN
               VT%T_CLOSE    = T(NM)
               VT%HEAT_INDEX_DEACTIVATE = 0
               ENDIF
            ENDDO
         ENDDO MESH_LOOP
C
         DO NN=1,NSPR
         SH=>SPRINKLER_HEAD(NN)
         IF (SH%HEAT_INDEX.EQ.N) THEN
            SH%T_ACT      = T_MIN
            SH%HEAT_INDEX = 0
            ENDIF
         ENDDO
C
         M=>MESH(1)
         IF (M%N_STRINGS+2.GT.M%N_STRINGS_MAX) THEN
            CALL RE_ALLOCATE_STRINGS(1)
            ENDIF
         M%N_STRINGS = M%N_STRINGS + 1
         WRITE(M%STRING(M%N_STRINGS),'(A)') 'HEAT_ACT'
         M%N_STRINGS = M%N_STRINGS + 1
         WRITE(M%STRING(M%N_STRINGS),'(I6,F10.2)') N,T_MIN
C
         HD%TMP_ACT = 1000000.
         ENDIF HEAT_DETECTOR_ACT
C
      ENDDO HEAT_DETECTOR_LOOP
C
C
      SPRINKLER_LOOP: DO N=1,NSPR
C
      SH => SPRINKLER_HEAD(N)
      IF (SH%ACT_CODE.NE.0) CYCLE SPRINKLER_LOOP
      S  => SPRINKLER(SH%INDEX)
C
      SPRINKLER_ACT: IF (SH%TMP_L.GE.S%TMP_ACT .OR.
     .                   T_MIN.GT.SH%T_ACT) THEN
         T_ACT_FIRST = MIN(T_ACT_FIRST,T_MIN)
         SH%ACT_CODE = 1
         SH%T        = MAX(T_ACT_FIRST+SYSTEM_DELAY,T_MIN+SH%DELAY)
C
         M=>MESH(1)
         IF (M%N_STRINGS+2.GT.M%N_STRINGS_MAX) 
     .      CALL RE_ALLOCATE_STRINGS(1)
         M%N_STRINGS = M%N_STRINGS + 1
         WRITE(M%STRING(M%N_STRINGS),'(A)') 'SPRK_ACT'
         M%N_STRINGS = M%N_STRINGS + 1
         WRITE(M%STRING(M%N_STRINGS),'(I6,F10.2)') N,T_MIN
C
         ENDIF SPRINKLER_ACT
C
      ENDDO SPRINKLER_LOOP
C
C
      SMOKE_DETECTOR_LOOP: DO N=1,NSD
C
      SD=>SMOKE_DETECTOR(N)
C
      SMOKE_DETECTOR_ACT: IF (SD%YSMOKE_IN.GT.SD%YSMOKE_ACT) THEN
C
         SD%YSMOKE_ACT = 1000000.
C
         M=>MESH(1)
         IF (M%N_STRINGS+2.GT.M%N_STRINGS_MAX) THEN
            CALL RE_ALLOCATE_STRINGS(1)
            ENDIF
         M%N_STRINGS = M%N_STRINGS + 1
         WRITE(M%STRING(M%N_STRINGS),'(A)') 'HEAT_ACT'
         M%N_STRINGS = M%N_STRINGS + 1
         WRITE(M%STRING(M%N_STRINGS),'(I6,F10.2)') N+NHD,T_MIN
C
         ENDIF SMOKE_DETECTOR_ACT
C
      ENDDO SMOKE_DETECTOR_LOOP
C
C     Dump out diagnostics
C
      IF (DIAGNOSTICS) CALL WRITE_DIAGNOSTICS(T)
C
C     Stop the run
C
      IF (T_MIN.GE.TWFIN .OR. STOP_CODE.GT.0) EXIT MAIN_LOOP
C
C     Flush Buffers (All Nodes)
C
      IF (MOD(ICYC,10).EQ.0) CALL FLUSH_BUFFERS
CEVAC:      IF (MOD(ICYC,10).EQ.0 .AND. ANY(EVACUATION_GRID))
CEVAC:     .     CALL FLUSH_BUFFER(LU121_EVAC)
C
      ENDDO MAIN_LOOP
C
C***********************************************************************
C                          END OF TIMESTEP
C***********************************************************************
C
      TUSED(1,1:NMESHES) = SECOND() - TUSED(1,1:NMESHES)
C
      CALL TIMINGS
C
      SELECT CASE(STOP_CODE)
      CASE(0)
            CALL SHUTDOWN('STOP: FDS completed successfully')
      CASE(1)
            CALL SHUTDOWN('STOP: Numerical Instability')
      CASE(2)
            CALL SHUTDOWN('STOP: FDS stopped by user')
      END SELECT
C
C
      CONTAINS
C
C
      SUBROUTINE INITIALIZE_MESH_EXCHANGE(NM)
C
C Create arrays by which info is to exchanged across meshes
C
      INTEGER IMIN,IMAX,JMIN,JMAX,KMIN,KMAX,NOM,IOR,IW
      INTEGER, INTENT(IN) :: NM
      TYPE (MESH_TYPE), POINTER :: M2,M
      LOGICAL FOUND
C
      M=>MESH(NM)
C
      ALLOCATE(M%OMESH(NMESHES))
C
      OTHER_MESH_LOOP: DO NOM=1,NMESHES
C
      IF (NOM.EQ.NM) CYCLE OTHER_MESH_LOOP
C
      M2=>MESH(NOM)
      IMIN=0
       IMAX=M2%IBP1
      JMIN=0
       JMAX=M2%JBP1
      KMIN=0
       KMAX=M2%KBP1
      NIC(NOM,NM) = 0
      FOUND = .FALSE.
      SEARCH_LOOP: DO IW=1,M%NEWC
      IF (M%IJKW(9,IW).NE.NOM) CYCLE SEARCH_LOOP
      NIC(NOM,NM) = NIC(NOM,NM) + 1
      FOUND = .TRUE.
      IOR = M%IJKW(4,IW)
      SELECT CASE(IOR)
      CASE( 1)
            IMIN=MAX(IMIN,M%IJKW(10,IW)-1)
      CASE(-1)
            IMAX=MIN(IMAX,M%IJKW(10,IW))
      CASE( 2)
            JMIN=MAX(JMIN,M%IJKW(11,IW)-1)
      CASE(-2)
            JMAX=MIN(JMAX,M%IJKW(11,IW))
      CASE( 3)
            KMIN=MAX(KMIN,M%IJKW(12,IW)-1)
      CASE(-3)
            KMAX=MIN(KMAX,M%IJKW(12,IW))
      END SELECT
      ENDDO SEARCH_LOOP
C
      IF ( M2%XS.GE.M%XS .AND. M2%XF.LE.M%XF .AND.
     .     M2%YS.GE.M%YS .AND. M2%YF.LE.M%YF .AND.
     .     M2%ZS.GE.M%ZS .AND. M2%ZF.LE.M%ZF ) FOUND = .TRUE.
C
      IF (.NOT.FOUND) CYCLE OTHER_MESH_LOOP
C
      I_MIN(NOM,NM) = IMIN
      I_MAX(NOM,NM) = IMAX
      J_MIN(NOM,NM) = JMIN
      J_MAX(NOM,NM) = JMAX
      K_MIN(NOM,NM) = KMIN
      K_MAX(NOM,NM) = KMAX
C
      ALLOCATE(M%OMESH(NOM)% TMP(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
      ALLOCATE(M%OMESH(NOM)%   H(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
      ALLOCATE(M%OMESH(NOM)%   U(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
      ALLOCATE(M%OMESH(NOM)%   V(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
      ALLOCATE(M%OMESH(NOM)%   W(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
      IF (NSPEC.GT.0) THEN
      ALLOCATE(M%OMESH(NOM)%  YY(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX,NSPEC))
      ALLOCATE(M%OMESH(NOM)% YYS(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX,NSPEC))
      ENDIF
C
C     Wall arrays
C
      ALLOCATE(M%OMESH(NOM)%IV(0:M2%NEWC))
      M%OMESH(NOM)%IV(0:M2%NEWC) = M2%IV(0:M2%NEWC)
      ALLOCATE(M%OMESH(NOM)%IJKW(12,M2%NEWC))
      M%OMESH(NOM)%IJKW(1:12,1:M2%NEWC) = M2%IJKW(1:12,1:M2%NEWC)
C
      ALLOCATE(M%OMESH(NOM)%WALL(0:M2%NEWC))
C
C     Particle and Droplet Orphan Arrays
C
      IF (DROPLET_FILE) THEN
      M%OMESH(NOM)%N_DROP_ORPHANS = 0
      M%OMESH(NOM)%N_DROP_ORPHANS_DIM = 1000
      ALLOCATE(M%OMESH(NOM)%DROPLET(M%OMESH(NOM)%N_DROP_ORPHANS_DIM),
     .         STAT=IZERO)
      CALL ChkMemErr('INIT','DROPLET',IZERO)
      ENDIF
C
      ENDDO OTHER_MESH_LOOP
C
      END SUBROUTINE INITIALIZE_MESH_EXCHANGE
C
C
      SUBROUTINE DOUBLE_CHECK(NM)
C
C Double check exchange pairs
C
      INTEGER NOM
      INTEGER, INTENT(IN) :: NM
      TYPE (MESH_TYPE), POINTER :: M2,M
C
      M=>MESH(NM)
C
      OTHER_MESH_LOOP: DO NOM=1,NMESHES
C
      IF (NOM.EQ.NM) CYCLE OTHER_MESH_LOOP
C
      IF (NIC(NM,NOM).EQ.0 .AND. NIC(NOM,NM).GT.0) THEN
C
      M2=>MESH(NOM)
      ALLOCATE(M%OMESH(NOM)%IJKW(12,M2%NEWC))
      ALLOCATE(M%OMESH(NOM)%IV(0:M2%NEWC))
      ALLOCATE(M%OMESH(NOM)%WALL(0:M2%NEWC))
C
      ENDIF
C
      ENDDO OTHER_MESH_LOOP
C
      END SUBROUTINE DOUBLE_CHECK
C
C
      SUBROUTINE MESH_EXCHANGE(CODE)
C
C Exchange Information between Meshes
C
      INTEGER, INTENT(IN) :: CODE
      INTEGER NM
C
      TNOW = SECOND()
C
      MESH_LOOP: DO NM=1,NMESHES
C
      OTHER_MESH_LOOP: DO NOM=1,NMESHES
C
      IF (CODE.EQ.0 .AND. NIC(NOM,NM).LT.1 .AND.
     .   NIC(NM,NOM).GT.0 .AND. I_MIN(NOM,NM).LT.0) THEN
      IF (RADIATION) THEN
      M =>MESH(NM)
      M2=>MESH(NOM)%OMESH(NM)
         DO IW=1,M%NEWC
         IF (M%IJKW(9,IW).EQ.NOM) THEN
            ALLOCATE(M2%WALL(IW)%ILW(NRA,NSB))
            M2%WALL(IW)%ILW = SIGMA*TMPA4*RPI
            ENDIF
         ENDDO
      ENDIF
      ENDIF
C
      IF (NIC(NOM,NM).EQ.0 .AND. NIC(NM,NOM).EQ.0) CYCLE OTHER_MESH_LOOP
C
      IF (CODE.GT.0) THEN
      IF (.NOT.ACTIVE_MESH(NM) .OR. .NOT.ACTIVE_MESH(NOM)) 
     .     CYCLE OTHER_MESH_LOOP
      ENDIF
C
      IF (DEBUG) THEN
      WRITE(0,*) NOM,' receiving data from ',NM,' code=',CODE
      ENDIF
C
      M =>MESH(NM)
      M2=>MESH(NOM)%OMESH(NM)
      M3=>MESH(NM)%OMESH(NOM)
      M4=>MESH(NOM)
C
      IMIN = I_MIN(NOM,NM)
      IMAX = I_MAX(NOM,NM)
      JMIN = J_MIN(NOM,NM)
      JMAX = J_MAX(NOM,NM)
      KMIN = K_MIN(NOM,NM)
      KMAX = K_MAX(NOM,NM)
C
      INITIALIZE_IF: IF (CODE.EQ.0) THEN
C
      IF (RADIATION) THEN
         DO IW=1,M%NEWC
         IF (M%IJKW(9,IW).EQ.NOM) THEN
            ALLOCATE(M2%WALL(IW)%ILW(NRA,NSB))
            M2%WALL(IW)%ILW = SIGMA*TMPA4*RPI
            ENDIF
         ENDDO
      ENDIF
C
      ENDIF INITIALIZE_IF
C
      PREDICTOR_IF: IF (CODE.EQ.1 .AND. NIC(NOM,NM).GT.0) THEN
C
      M2%TMP(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)=
     . M%TMP(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
      M2%H(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)=
     . M%H(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
      IF (NSPEC.GT.0)
     .M2%YYS(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX,1:NSPEC)=
     . M%YYS(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX,1:NSPEC)
C
      ENDIF PREDICTOR_IF
C
      CORRECTOR_IF: IF (CODE.EQ.0 .OR. CODE.EQ.2) THEN
C
      IF (NIC(NOM,NM).GT.0) THEN
      M2%IV(0:M%NEWC) = M%IV(0:M%NEWC)
      M2%TMP(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)=
     . M%TMP(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
      M2%H(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)=
     . M%H(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
      M2%U(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)=
     . M%U(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
      M2%V(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)=
     . M%V(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
      M2%W(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)=
     . M%W(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
      IF (NSPEC.GT.0)
     .M2%YY(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX,1:NSPEC)=
     . M%YY(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX,1:NSPEC)
      ENDIF
C
C     Radiation Exchange
C
      RADIATION_IF: IF (RADIATION .AND. CODE.EQ.2 .AND. 
     .                  NIC(NOM,NM).GT.0) THEN
C
      DO IW=1,M4%NEWC
      IF (M4%IJKW(9,IW).EQ.NM .AND. M4%IV(IW).EQ.4) THEN
      M4%WALL(IW)%ILW(1:NRA,1:NSB) = M3%WALL(IW)%ILW(1:NRA,1:NSB)
      ENDIF
      ENDDO
C
      ENDIF RADIATION_IF
C
      ENDIF CORRECTOR_IF
C
C     Get Number of Droplet Orphans
C
cc    IF (DROPLET_FILE .AND. NIC(NOM,NM).GT.0) THEN
      IF (DROPLET_FILE) THEN 
C
      M2%N_DROP_ADOPT = MIN(M3%N_DROP_ORPHANS,N_DROP_ADOPT_MAX)
      IF (M4%NLP+M2%N_DROP_ADOPT.GT.M4%NLPDIM) THEN
      CALL RE_ALLOCATE_DROPLETS(1,NOM,0)
      ENDIF
C
      ENDIF
C
C     Sending/Receiving Droplet Buffer Arrays
C
cc    IF_DROPLETS: IF (DROPLET_FILE .AND. NIC(NOM,NM).GT.0) THEN 
      IF_DROPLETS: IF (DROPLET_FILE) THEN 
C
      IF_DROPLETS_SENT: IF (M2%N_DROP_ADOPT.GT.0) THEN
C
      M4%DROPLET(M4%NLP+1:M4%NLP+M2%N_DROP_ADOPT)= 
     .M3%DROPLET(1:M2%N_DROP_ADOPT) 
C
      M4%NLP = M4%NLP + M2%N_DROP_ADOPT
      M3%N_DROP_ORPHANS = 0
C
      ENDIF IF_DROPLETS_SENT
C
      ENDIF IF_DROPLETS
C
      ENDDO OTHER_MESH_LOOP
      ENDDO MESH_LOOP
C
      TUSED(11,:)=TUSED(11,:) + SECOND() - TNOW
      END SUBROUTINE MESH_EXCHANGE
C
C
      SUBROUTINE WRITE_STRINGS
C
      INTEGER :: N,NM
C
C     Write character strings out to the .smv file
C
      MESH_LOOP: DO NM=1,NMESHES
      DO N=1,MESH(NM)%N_STRINGS
      WRITE(LU4,'(A)') TRIM(MESH(NM)%STRING(N))
      ENDDO
      MESH(NM)%N_STRINGS = 0
      ENDDO MESH_LOOP
C
      END SUBROUTINE WRITE_STRINGS
C
C
      SUBROUTINE EXCHANGE_DIAGNOSTICS
C
      INTEGER  :: NM,NECYC,I
      REAL(EB) :: T_SUM
C
      TNOW = SECOND()
C
      MESH_LOOP: DO NM=1,NMESHES

      T_SUM = 0.
      SUM_LOOP: DO I=2,N_TIMERS
      IF (I.EQ.9 .OR. I.EQ.10) CYCLE SUM_LOOP
      T_SUM = T_SUM + TUSED(I,NM)
      ENDDO SUM_LOOP
C
      NECYC          = MAX(1,NTCYC(NM)-NCYC(NM))
      T_PER_STEP(NM) = (T_SUM-T_ACCUM(NM))/REAL(NECYC,EB)
      T_ACCUM(NM)    = T_SUM
      NCYC(NM)       = NTCYC(NM)
C
      ENDDO MESH_LOOP
C
      TUSED(11,:) = TUSED(11,:) + SECOND() - TNOW
      END SUBROUTINE EXCHANGE_DIAGNOSTICS
C
C
      END PROGRAM FDS
