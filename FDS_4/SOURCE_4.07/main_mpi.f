      PROGRAM FDS  
C
C Fire Dynamics Simulator, Main Program, Parallel Version
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
C     USE EVAC
C
      IMPLICIT NONE
C
C     Miscellaneous declarations
C
      LOGICAL  :: EX,DIAGNOSTICS,EXCHANGE_RADIATION=.TRUE.
      INTEGER  :: LO10,NM,IZERO,STOP_CODE=0,DATE_TIME(8),NN
      CHARACTER(10) :: BIG_BEN(3)
      REAL(EB) :: T_MAX,T_MIN
      REAL(EB), ALLOCATABLE, DIMENSION(:) ::  T,TC_GLB,TC_LOC,DT_SYNC
      INTEGER, ALLOCATABLE, DIMENSION(:) ::  ISTOP,COUNT_GLB,COUNT_LOC
      LOGICAL, ALLOCATABLE, DIMENSION(:) ::  ACTIVE_MESH
      INTEGER NOM,IWW,IW
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
C     MPI stuff
C
c     INCLUDE '/usr/local/include/mpif.h'
      INCLUDE 'mpif.h'
      INTEGER :: N,MYID=0,NUMPROCS=1,I,IERR,STATUS(MPI_STATUS_SIZE)
      INTEGER :: RNODE,BUFFER_SIZE,TAG,PNAMELEN
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: TAGS
      INTEGER, ALLOCATABLE, DIMENSION(:) :: REQ
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ARRAY_OF_STATUSES
      INTEGER :: N_REQ
      CHARACTER(MPI_MAX_PROCESSOR_NAME) PNAME
C
C     Initialize MPI
C
      CALL MPI_INIT(IERR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID, IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NUMPROCS, IERR)
      CALL MPI_GET_PROCESSOR_NAME(PNAME, PNAMELEN, IERR)
C
      WRITE(LU0,'(A,I2,A,I2,A,A)') 'Mesh ',MYID+1,' of ',
     .      NUMPROCS,' is alive on ',PNAME(1:PNAMELEN)
C
C     Assign a compilation date (All Nodes)
C
      COMPILE_DATE = 'March 10, 2006'
      VERSION      = 4.07
C
C     Open the input file (All Nodes)
C
      CALL OPEN_INPUT_FILE
C
C     Read input from CHID.data file (All Nodes)
C
      CALL READ_DATA
C
      IF (NMESHES.NE.NUMPROCS) 
     .CALL SHUTDOWN('ERROR: Number of meshes not equal to '//
     .'number of threads')
C
C     Read input for EVACUATION routine
C
C     IF (ANY(EVACUATION_GRID)) CALL READ_EVAC
C
C     Open and write to Smokeview file (Master Node Only)
C
      IF (MYID.EQ.0) CALL WRITE_SMOKEVIEW_FILE
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
      ALLOCATE(COUNT_LOC(NTC),STAT=IZERO)
      CALL ChkMemErr('MAIN','COUNT_LOC',IZERO) ; COUNT_LOC = 0
      ALLOCATE(COUNT_GLB(NTC),STAT=IZERO)
      CALL ChkMemErr('MAIN','COUNT_GLB',IZERO) ; COUNT_GLB = 0
      ALLOCATE(TC_GLB(MAX(NTC,NSPR,NHD,NSD)),STAT=IZERO)
      CALL ChkMemErr('MAIN','TC_GLB',IZERO)
      ALLOCATE(TC_LOC(MAX(NTC,NSPR,NHD,NSD)),STAT=IZERO)
      CALL ChkMemErr('MAIN','TC_LOC',IZERO)
      T     = 0.
      ISTOP = 0
C
C     Create unique tags for all mesh exchanges
C
      ALLOCATE(REQ(NMESHES*NMESHES*10)) ; REQ = MPI_REQUEST_NULL
      ALLOCATE(ARRAY_OF_STATUSES(MPI_STATUS_SIZE,NMESHES*NMESHES*10))
      ALLOCATE(TAGS(NMESHES,NMESHES,0:2))
      TAG = 0
      DO NM=1,NMESHES
      DO NOM=NM,NMESHES
      TAG = TAG+1
      TAGS(NM,NOM,0) = TAG
      TAGS(NOM,NM,0) = TAG
      ENDDO
      ENDDO
      TAGS(:,:,1) = TAGS(:,:,0) + 1000
      TAGS(:,:,2) = TAGS(:,:,0) + 2000
C
C     Initialize global parameters (All Nodes)
C
      CALL INITIALIZE_GLOBAL_VARIABLES
      CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
C
C     Initialize radiation (All Nodes)
C
      IF (RADIATION) CALL INIT_RADIATION
      CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
C
C     Allocate and initialize mesh-specific variables
C
      DO NM=MYID+1,NMESHES,NUMPROCS
      CALL INITIALIZE_MESH_VARIABLES(NM)
      ENDDO
      CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
C
C     Allocate and initialize mesh variable exchange arrays
C
      DO NM=MYID+1,NMESHES,NUMPROCS
      CALL INITIALIZE_MESH_EXCHANGE(NM)
      ENDDO
      CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
C
      CALL MPI_ALLGATHER(I_MIN(1,MYID+1),NMESHES,MPI_INTEGER,I_MIN,
     .                   NMESHES,MPI_INTEGER,MPI_COMM_WORLD,IERR)
      CALL MPI_ALLGATHER(I_MAX(1,MYID+1),NMESHES,MPI_INTEGER,I_MAX,
     .                   NMESHES,MPI_INTEGER,MPI_COMM_WORLD,IERR)
      CALL MPI_ALLGATHER(J_MIN(1,MYID+1),NMESHES,MPI_INTEGER,J_MIN,
     .                   NMESHES,MPI_INTEGER,MPI_COMM_WORLD,IERR)
      CALL MPI_ALLGATHER(J_MAX(1,MYID+1),NMESHES,MPI_INTEGER,J_MAX,
     .                   NMESHES,MPI_INTEGER,MPI_COMM_WORLD,IERR)
      CALL MPI_ALLGATHER(K_MIN(1,MYID+1),NMESHES,MPI_INTEGER,K_MIN,
     .                   NMESHES,MPI_INTEGER,MPI_COMM_WORLD,IERR)
      CALL MPI_ALLGATHER(K_MAX(1,MYID+1),NMESHES,MPI_INTEGER,K_MAX,
     .                   NMESHES,MPI_INTEGER,MPI_COMM_WORLD,IERR)
      CALL MPI_ALLGATHER(NIC(1,MYID+1),  NMESHES,MPI_INTEGER,NIC,
     .                   NMESHES,MPI_INTEGER,MPI_COMM_WORLD,IERR)
C
      I_MIN = TRANSPOSE(I_MIN)
      I_MAX = TRANSPOSE(I_MAX)
      J_MIN = TRANSPOSE(J_MIN)
      J_MAX = TRANSPOSE(J_MAX)
      K_MIN = TRANSPOSE(K_MIN)
      K_MAX = TRANSPOSE(K_MAX)
      NIC   = TRANSPOSE(NIC)
C
      DO NM=MYID+1,NMESHES,NUMPROCS
      CALL DOUBLE_CHECK(NM)
      ENDDO
      CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
C
C     Initialize miscellaneous routines
C
      IF (PAPERMODEL) CALL INIT_SOLID_PHASE_REACTIONS
C
C     Potentially read data from a previous calculation 
C
      DO NM=MYID+1,NMESHES,NUMPROCS
      IF (RESTART) CALL READ_CORE(T(NM),NM)
      ENDDO
      CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
C
C     Initialize output files containing global data (Master Node Only)
C
      IF (MYID.EQ.0) CALL INITIALIZE_GLOBAL_DUMPS
      CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
C
C     Initialize output files that are mesh-specific
C
      DO NM=MYID+1,NMESHES,NUMPROCS
      CALL INITIALIZE_MESH_DUMPS(NM)
      CALL INITIALIZE_DROPLETS(NM)
      CALL INITIALIZE_TREES(NM)
C     IF (ANY(EVACUATION_GRID)) CALL INITIALIZE_EVACUATION(NM)
      CALL POST_RECEIVES(NM,0)
      ENDDO
      CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
C
C     Write out character strings to .smv file
C
      CALL WRITE_STRINGS
C
C     Initialize Mesh Exchange Arrays (All Nodes)
C
      CALL MESH_EXCHANGE(0)
      CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)


c     ACTIVE_MESH=.TRUE.
c     DO NM=MYID+1,NMESHES,NUMPROCS
c     write(0,*) 'mesh ',nm,' call post receives'
c     CALL POST_RECEIVES(NM,2)
c     ENDDO
c     CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
c     write(0,*) 'mesh ',myid+1,' call mesh_exchange'
c     CALL MESH_EXCHANGE(2)
c     write(0,*) 'mesh ',myid+1,' exited mesh_exchange'

C
C ********************************************************************
C                      MAIN TIMESTEPPING LOOP
C ********************************************************************
C
      MAIN_LOOP: DO  
C
      ICYC  = ICYC + 1 
      IF (ANY(EVACUATION_GRID) .AND. ICYC.EQ.0) ICYC = 1
C
      IF (MOD(ICYC,3).EQ.0 .AND. TIMING) THEN
      CALL DATE_AND_TIME(BIG_BEN(1),BIG_BEN(2),BIG_BEN(3),DATE_TIME)
      WRITE(0,'(A,I2,A,I6,A,I2,A,I3.3)') 
     .' Thread ',MYID+1,' starts iteration',ICYC,' at ',
     .DATE_TIME(7),'.',DATE_TIME(8)
      ENDIF
C
C     Radiation Exchange
C
      EXCHANGE_RADIATION = .FALSE.
      IF (RADIATION) THEN
      IF (MOD(ICYC,ANGLE_INCREMENT*TIME_STEP_INCREMENT).EQ.0)
     .   EXCHANGE_RADIATION = .TRUE.
      ENDIF
C
C     Synchronize clocks
C
      CALL MPI_ALLGATHER(T(MYID+1),1,MPI_DOUBLE_PRECISION,T,1,
     .                MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERR)
C
C     Check for program stops
C
      INQUIRE(FILE=TRIM(CHID)//'.stop',EXIST=EX)
      IF (EX) ISTOP = 2
      CALL MPI_ALLGATHER(ISTOP(MYID+1),1,MPI_INTEGER,ISTOP,1,
     .                MPI_INTEGER,MPI_COMM_WORLD,IERR)
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
C     Time step
C
      IF (SYNCHRONIZE) THEN
C
      DT_SYNC(MYID+1) = MESH(MYID+1)%DTNEXT
      CALL MPI_ALLGATHER(DT_SYNC(MYID+1),1,MPI_DOUBLE_PRECISION,
     .     DT_SYNC,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERR)
      IF (SYNC_TIME_STEP(MYID+1)) THEN
         MESH(MYID+1)%DTNEXT = MINVAL(DT_SYNC,MASK=SYNC_TIME_STEP)
         T(MYID+1) = MINVAL(T,MASK=SYNC_TIME_STEP)
         ACTIVE_MESH(MYID+1) = .TRUE.
         ELSE
         ACTIVE_MESH(MYID+1) = .FALSE.
         IF (T(MYID+1)+MESH(MYID+1)%DTNEXT .LE. T_MAX) 
     .      ACTIVE_MESH(MYID+1) = .TRUE.
         IF (STOP_CODE.GT.0) ACTIVE_MESH(MYID+1) = .TRUE.
      ENDIF
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
C     When to dump out diagnostics to the .out file
C
      DIAGNOSTICS = .FALSE.
      LO10 = LOG10(REAL(ICYC,EB))
      IF (MOD(ICYC,10**LO10).EQ.0 .OR. MOD(ICYC,100).EQ.0 .OR.
     .    T_MIN.GE.TWFIN .OR. STOP_CODE.GT.0) DIAGNOSTICS = .TRUE.
C
C     Give every processor the full ACTIVE_MESH array
C
      CALL MPI_ALLGATHER(ACTIVE_MESH(MYID+1), 1, MPI_LOGICAL, 
     .      ACTIVE_MESH,1, MPI_LOGICAL, MPI_COMM_WORLD, IERR)
C
C     If no meshes are due to be updated, update them all
C
      IF (ALL(.NOT.ACTIVE_MESH)) ACTIVE_MESH = .TRUE.
C
C     Do not do EVACuation if past the max iteration criteria
C
c     IF (ANY(EVACUATION_GRID)) THEN
c       EVAC_DT = 1000000.
c       DO NM=1,NMESHES
c         IF (.NOT.EVACUATION_ONLY(NM))
c    .         EVAC_DT = MIN(EVAC_DT,MESH(NM)%DTNEXT)
c       ENDDO
c       DO NM=1,NMESHES
c         IF (EVACUATION_ONLY(NM)) THEN
c           IF (ICYC .GT. EVAC_TIME_ITERATIONS) THEN
c             ACTIVE_MESH(NM) = .FALSE.
c             EVAC_DT = MIN(EVAC_DT, EVAC_DT_STEADY_STATE)
c             MESH(NM)%DT     = EVAC_DT
c             T(NM)           = T(NM) + MESH(NM)%DT
c             IF (EVACUATION_GRID(NM) ) THEN
c               CALL EVACUATE_HUMANS(T(NM),NM)
c               IF (T(NM).GE.PARCLK(NM)) THEN
c                 CALL DUMP_EVAC(PARCLK(NM),NM)
c                 DO
c                   PARCLK(NM) = PARCLK(NM) + WPAR
c                   IF (PARCLK(NM).GE.T(NM)) EXIT
c                 ENDDO
c               ENDIF
c             ENDIF
c           ELSE
c             ACTIVE_MESH(NM) = .TRUE.
c             EVAC_DT = MIN(EVAC_DT, EVAC_DT_FLOWFIELD)
c           ENDIF
c         ENDIF
c       ENDDO
c       DO NM=1,NMESHES
c         IF (EVACUATION_ONLY(NM)) MESH(NM)%DTNEXT = EVAC_DT
c       ENDDO
c     ENDIF
C
      PREDICTOR = .TRUE. ; CORRECTOR = .FALSE.
C
      IF (DEBUG)
     .WRITE(0,*) 'Cycle ',ICYC,' Mesh ',MYID+1,
     .   ' starting',ACTIVE_MESH(MYID+1)
C
      IF (MOD(ICYC,3).EQ.0 .AND. TIMING) THEN
      CALL DATE_AND_TIME(BIG_BEN(1),BIG_BEN(2),BIG_BEN(3),DATE_TIME)
      IF (ACTIVE_MESH(MYID+1))
     .WRITE(0,'(A,I2,A,I2,A,I3.3)')
     .' Thread ',MYID+1,' is active at ',
     .DATE_TIME(7),'.',DATE_TIME(8)
      ENDIF
C
C=====================================================================
C
      MULTIBLOCK1: DO NM=MYID+1,NMESHES,NUMPROCS

      IF (.NOT.ACTIVE_MESH(NM)) CYCLE MULTIBLOCK1
C
      MESH(NM)%DT = MESH(NM)%DTNEXT
      NTCYC(NM)   = NTCYC(NM) + 1
C
C     Set up receiving calls
C
      CALL POST_RECEIVES(NM,1)
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
C Evacuation
C
C     IF (EVACUATION_ONLY(NM)) CALL EVACUATE_HUMANS(T(NM),NM)
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
C
c     IF (ANY(EVACUATION_GRID) .AND. EVACUATION_ONLY(NM)) THEN
c        PRESSURE_ITERATION_LOOP: DO N=1,EVAC_PRESSURE_ITERATIONS
c           CALL NO_FLUX(T(NM))
c           CALL PRESSURE_SOLVER(NM)
c        ENDDO PRESSURE_ITERATION_LOOP
c     ENDIF
C
C Compute velocity estimates at next time step
C
      CALL VELOCITY_PREDICTOR(T(NM),NM,ISTOP(NM))
      IF (ISTOP(NM).EQ.1) THEN
c        OPEN(999,FILE=TRIM(CHID)//'.stop')
c        CLOSE(999)
         ENDIF
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
C Modeling of smoke detector at the next time step
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
      IF (MOD(ICYC,3).EQ.0.AND.TIMING.AND.ACTIVE_MESH(MYID+1)) THEN
      CALL DATE_AND_TIME(BIG_BEN(1),BIG_BEN(2),BIG_BEN(3),DATE_TIME)
      WRITE(0,'(A,I2,A,I2,A,I3.3)') 
     .' Thread ',MYID+1,' enters Mesh Exchange 1 at ',
     .DATE_TIME(7),'.',DATE_TIME(8)
      ENDIF
C
      CALL MESH_EXCHANGE(1)
C
      IF (MOD(ICYC,3).EQ.0.AND.TIMING.AND.ACTIVE_MESH(MYID+1)) THEN
      CALL DATE_AND_TIME(BIG_BEN(1),BIG_BEN(2),BIG_BEN(3),DATE_TIME)
      WRITE(0,'(A,I2,A,I2,A,I3.3)') 
     .' Thread ',MYID+1,' exits  Mesh Exchange 1 at ',
     .DATE_TIME(7),'.',DATE_TIME(8)
      ENDIF
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      CORRECTOR = .TRUE. ; PREDICTOR = .FALSE.
      MULTIBLOCK2: DO NM=MYID+1,NMESHES,NUMPROCS
C
      IF (.NOT.ACTIVE_MESH(NM)) CYCLE MULTIBLOCK2
C
C     Set up receiving arrays
C
      CALL POST_RECEIVES(NM,2)
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
C Evacuation
C
C     IF (EVACUATION_ONLY(NM)) CALL EVACUATE_HUMANS(T(NM),NM)
C
C Compute the divergence and time derivative of divergence
C
      CALL DIVERGENCE(T(NM),NM)
C
C Solve the Poisson equation for pressure
C
      CALL PRESSURE_SOLVER(NM)
C
c     IF (ANY(EVACUATION_GRID) .AND. EVACUATION_ONLY(NM)) THEN
c        PRESSURE_ITERATION_LOOP2: DO N=1,EVAC_PRESSURE_ITERATIONS
c           CALL NO_FLUX(T(NM))
c           CALL PRESSURE_SOLVER(NM)
c        ENDDO PRESSURE_ITERATION_LOOP2
c     ENDIF
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
C Check to see if the smoke detetcor went off 
C
      IF (NSD .GT.0) CALL SMOKE_DETECTORS(T(NM),NM)
C
C Mesh-specific data updates and dumps
C
      CALL UPDATE_HRR(NM)
      IF (NSPEC.GT.0) CALL UPDATE_MASS(NM)
      IF (NTC.GT.0) CALL UPDATE_TC(NM)
C
C EVACuation calculation and dumps for EVACUATION_ONLY
C
c      EVACUATION_DUMP: IF (ANY(EVACUATION_GRID) .AND.
c     .     EVACUATION_ONLY(NM)) THEN
c        IF (EVACUATION_GRID(NM) ) THEN
c          CALL EVACUATE_HUMANS(T(NM),NM)
c          IF (T(NM).GE.PARCLK(NM)) THEN
c            CALL DUMP_EVAC(PARCLK(NM),NM)
c            DO
c              PARCLK(NM) = PARCLK(NM) + WPAR
c              IF (PARCLK(NM).GE.T(NM)) EXIT
c            ENDDO
c          ENDIF
c        ENDIF
C     Dump the EVAC flowfieds for all EVAC meshes
c        CALL DUMP_SF(T(NM),NM,0)
C
c      ELSE
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
c      ENDIF EVACUATION_DUMP
C
      IF (DIAGNOSTICS) CALL CHECK_DIVERGENCE(NM)
C
      ENDDO MULTIBLOCK2
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     Exchange information among meshes
C
      IF (MOD(ICYC,3).EQ.0.AND.TIMING.AND.ACTIVE_MESH(MYID+1)) THEN
      CALL DATE_AND_TIME(BIG_BEN(1),BIG_BEN(2),BIG_BEN(3),DATE_TIME)
      WRITE(0,'(A,I2,A,I2,A,I3.3)') 
     .' Thread ',MYID+1,' enters Mesh Exchange 2 at ',
     .DATE_TIME(7),'.',DATE_TIME(8)
      ENDIF
C
      CALL MESH_EXCHANGE(2)
C
      IF (MOD(ICYC,3).EQ.0.AND.TIMING.AND.ACTIVE_MESH(MYID+1)) THEN
      CALL DATE_AND_TIME(BIG_BEN(1),BIG_BEN(2),BIG_BEN(3),DATE_TIME)
      WRITE(0,'(A,I2,A,I2,A,I3.3)') 
     .' Thread ',MYID+1,' exits  Mesh Exchange 2 at ',
     .DATE_TIME(7),'.',DATE_TIME(8)
      endif
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     Write character strings out to the .smv file
C
      IF (DIAGNOSTICS) CALL WRITE_STRINGS
C
C     Exchange info for diagnostic print out
C
      IF (DIAGNOSTICS) CALL EXCHANGE_DIAGNOSTICS
C
C     Dump out HRR info
C
      IF_DUMP_HRR: IF (T_MIN.GE.HRRCLK) THEN
      CALL MPI_ALLGATHER(HRR_COUNT(MYID+1),1,MPI_DOUBLE_PRECISION,
     .         HRR_COUNT,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERR)
Ccc      IF (MINVAL(HRR_COUNT).GT.0.) THEN
      IF (MINVAL(HRR_COUNT,MASK=.NOT.EVACUATION_ONLY).GT.0.) THEN
      CALL MPI_GATHER(HRR_SUM(MYID+1), 1, MPI_DOUBLE_PRECISION,
     .         HRR_SUM, 1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      CALL MPI_GATHER(RHRR_SUM(MYID+1),1, MPI_DOUBLE_PRECISION,
     .         RHRR_SUM,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      CALL MPI_GATHER(CHRR_SUM(MYID+1),1, MPI_DOUBLE_PRECISION,
     .         CHRR_SUM,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      CALL MPI_GATHER(FHRR_SUM(MYID+1),1, MPI_DOUBLE_PRECISION,
     .         FHRR_SUM,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      CALL MPI_GATHER(MLR_SUM(MYID+1), 1, MPI_DOUBLE_PRECISION,
     .         MLR_SUM, 1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      IF (MYID.EQ.0) CALL DUMP_HRR(T_MIN)
      HRRCLK = HRRCLK + DTHRR
      HRR_SUM   = 0.
      RHRR_SUM  = 0.
      CHRR_SUM  = 0.
      FHRR_SUM  = 0.
      MLR_SUM   = 0.
      HRR_COUNT = 0.
      ENDIF
      ENDIF IF_DUMP_HRR
C
C     Dump out Evac info: EVAC_TODO next lines should be
C     made to work also in the parallel code.
C
c      IF (T_MIN.GE.EVACCLK .AND. ANY(EVACUATION_GRID)) THEN
c        CALL DUMP_EVAC_CSV(T_MIN)
c        EVACCLK = EVACCLK + DTHRR
c      ENDIF
C
C     Dump out Mass info
C
      IF_DUMP_MASS: IF (T_MIN.GE.MINTCLK) THEN
      CALL MPI_ALLGATHER(MINT_COUNT(MYID+1),1,MPI_DOUBLE_PRECISION,
     .     MINT_COUNT,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERR)
Ccc      IF (MINVAL(MINT_COUNT).GT.0.) THEN
      IF (MINVAL(MINT_COUNT,MASK=.NOT.EVACUATION_ONLY).GT.0.) THEN
      CALL MPI_GATHER(MINT_SUM(0,MYID+1),21,MPI_DOUBLE_PRECISION,
     .     MINT_SUM,21,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, IERR)
      IF (MYID.EQ.0) CALL DUMP_MASS(T_MIN)
      MINTCLK    = MINTCLK + DTMINT
      MINT_SUM   = 0.
      MINT_COUNT = 0.
      ENDIF
      ENDIF IF_DUMP_MASS
C
C     Dump out TC data
C
      IF_DUMP_TC: IF (T_MIN.GE.TCCLK .AND. NTC.GT.0) THEN
      COUNT_LOC(1:NTC) = THERMOCOUPLE(1:NTC)%COUNT
      CALL MPI_ALLREDUCE(COUNT_LOC(1),COUNT_GLB(1),NTC,MPI_INTEGER,
     .                MPI_SUM,MPI_COMM_WORLD,IERR)
      IF (MINVAL(COUNT_GLB).GT.0) THEN
      TC_LOC(1:NTC) = THERMOCOUPLE(1:NTC)%VALUE 
      CALL MPI_REDUCE(TC_LOC(1),TC_GLB(1),NTC,MPI_DOUBLE_PRECISION,
     .                MPI_SUM,0,MPI_COMM_WORLD,IERR)
      IF (MYID.EQ.0) THEN
         THERMOCOUPLE(1:NTC)%VALUE = TC_GLB(1:NTC)
         THERMOCOUPLE(1:NTC)%COUNT = COUNT_GLB(1:NTC)
         CALL DUMP_TC(T_MIN)
         ENDIF
      TCCLK = TCCLK + DTTC
      THERMOCOUPLE(1:NTC)%VALUE = 0.
      THERMOCOUPLE(1:NTC)%COUNT = 0
      ENDIF
      ENDIF IF_DUMP_TC
C
C     Dump out Sprinkler and Heat Detector link temps
C
      IF (T_MIN.GT.SPRKCLK) THEN
      TC_LOC(1:NSPR) = SPRINKLER_HEAD(1:NSPR)%TMP_L 
      CALL MPI_ALLREDUCE(TC_LOC(1),TC_GLB(1),NSPR,MPI_DOUBLE_PRECISION,
     .                MPI_MAX,MPI_COMM_WORLD,IERR)
      SPRINKLER_HEAD(1:NSPR)%TMP_L = TC_GLB(1:NSPR)
      IF (MYID.EQ.0) CALL DUMP_DETECTORS(T_MIN,1)
      SPRKCLK = SPRKCLK + DT_SPRK
      ENDIF
C
      IF (T_MIN.GT.HEATCLK) THEN
      TC_LOC(1:NHD) = HEAT_DETECTOR(1:NHD)%TMP_L 
      CALL MPI_ALLREDUCE(TC_LOC(1),TC_GLB(1),NHD,MPI_DOUBLE_PRECISION,
     .                   MPI_MAX,MPI_COMM_WORLD,IERR)
      HEAT_DETECTOR(1:NHD)%TMP_L = TC_GLB(1:NHD)
      IF (MYID.EQ.0) CALL DUMP_DETECTORS(T_MIN,2)
      HEATCLK = HEATCLK + DT_HEAT
      ENDIF
C
      IF (T_MIN.GT.SMOKECLK) THEN
      TC_LOC(1:NSD) = SMOKE_DETECTOR(1:NSD)%YSMOKE_IN
      CALL MPI_ALLREDUCE(TC_LOC(1),TC_GLB(1),NSD,MPI_DOUBLE_PRECISION,
     .                   MPI_MAX,MPI_COMM_WORLD,IERR)
      SMOKE_DETECTOR(1:NSD)%YSMOKE_IN = TC_GLB(1:NSD)
      IF (MYID.EQ.0) CALL DUMP_DETECTORS(T_MIN,3)
      SMOKECLK = SMOKECLK + DT_SMOKE
      ENDIF
C
C     Check for heat detector and sprinkler activation
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
         IF (MYID.EQ.0) THEN
         M=>MESH(1)
         IF (M%N_STRINGS+2.GT.M%N_STRINGS_MAX) THEN
            CALL RE_ALLOCATE_STRINGS(1)
            ENDIF
         M%N_STRINGS = M%N_STRINGS + 1
         WRITE(M%STRING(M%N_STRINGS),'(A)') 'HEAT_ACT'
         M%N_STRINGS = M%N_STRINGS + 1
         WRITE(M%STRING(M%N_STRINGS),'(I3,F8.2)') N,T_MIN
         ENDIF
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
         IF (MYID.EQ.0) THEN
         M=>MESH(1)
         IF (M%N_STRINGS+2.GT.M%N_STRINGS_MAX) 
     .      CALL RE_ALLOCATE_STRINGS(1)
         M%N_STRINGS = M%N_STRINGS + 1
         WRITE(M%STRING(M%N_STRINGS),'(A)') 'SPRK_ACT'
         M%N_STRINGS = M%N_STRINGS + 1
         WRITE(M%STRING(M%N_STRINGS),'(I3,F8.2)') N,T_MIN
         ENDIF
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
         IF (MYID.EQ.0) THEN
         M=>MESH(1)
         IF (M%N_STRINGS+2.GT.M%N_STRINGS_MAX) THEN
            CALL RE_ALLOCATE_STRINGS(1)
            ENDIF
         M%N_STRINGS = M%N_STRINGS + 1
         WRITE(M%STRING(M%N_STRINGS),'(A)') 'HEAT_ACT'
         M%N_STRINGS = M%N_STRINGS + 1
         WRITE(M%STRING(M%N_STRINGS),'(I3,F8.2)') N+NHD,T_MIN
         ENDIF
C
         ENDIF SMOKE_DETECTOR_ACT
C
      ENDDO SMOKE_DETECTOR_LOOP
C
C     Dump out diagnostics
C
      IF (MYID.EQ.0 .AND. DIAGNOSTICS) CALL WRITE_DIAGNOSTICS(T)
C
C     Stop the run
C
      IF (T_MIN.GE.TWFIN .OR. STOP_CODE.GT.0) EXIT MAIN_LOOP
C
C     Flush Buffers (All Nodes)
C
      IF (MOD(ICYC,10).EQ.0) CALL FLUSH_BUFFERS
      IF (MOD(ICYC,10).EQ.0 .AND. ANY(EVACUATION_GRID))
     .     CALL FLUSH_BUFFER(LU121_EVAC)
C
      IF (MOD(ICYC,3).EQ.0 .AND. TIMING) THEN
      CALL DATE_AND_TIME(BIG_BEN(1),BIG_BEN(2),BIG_BEN(3),DATE_TIME)
      WRITE(0,'(A,I2,A,I6,A,I2,A,I3.3)') 
     .' Thread ',MYID+1,' ends iteration',ICYC,' at ',
     .DATE_TIME(7),'.',DATE_TIME(8)
      ENDIF
C
      ENDDO MAIN_LOOP
C
C***********************************************************************
C                          END OF TIMESTEP
C***********************************************************************
C
      TUSED(1,MYID+1) = SECOND() - TUSED(1,MYID+1)
      CALL MPI_GATHER(TUSED(1,MYID+1),N_TIMERS,MPI_DOUBLE_PRECISION,
     .           TUSED,N_TIMERS,MPI_DOUBLE_PRECISION,0,
     .           MPI_COMM_WORLD,IERR)
C
      IF (MYID.EQ.0) CALL TIMINGS
C
      CALL MPI_FINALIZE(IERR)
C
      SELECT CASE(STOP_CODE)
      CASE(0) ; CALL SHUTDOWN('STOP: FDS completed successfully')
      CASE(1) ; CALL SHUTDOWN('STOP: Numerical Instability')
      CASE(2) ; CALL SHUTDOWN('STOP: FDS stopped by user')
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
      ALLOCATE(MESH(NM)%OMESH(NMESHES))
C
      OTHER_MESH_LOOP: DO NOM=1,NMESHES
C
      IF (NOM.EQ.NM) CYCLE OTHER_MESH_LOOP
C
      M2=>MESH(NOM)
      IMIN=0 ; IMAX=M2%IBP1
      JMIN=0 ; JMAX=M2%JBP1
      KMIN=0 ; KMAX=M2%KBP1
      NIC(NOM,NM) = 0
      FOUND = .FALSE.
      SEARCH_LOOP: DO IW=1,M%NEWC
      IF (M%IJKW(9,IW).NE.NOM) CYCLE SEARCH_LOOP
      NIC(NOM,NM) = NIC(NOM,NM) + 1
      FOUND = .TRUE.
      IOR = M%IJKW(4,IW)
      SELECT CASE(IOR)
      CASE( 1) ; IMIN=MAX(IMIN,M%IJKW(10,IW)-1)
      CASE(-1) ; IMAX=MIN(IMAX,M%IJKW(10,IW))
      CASE( 2) ; JMIN=MAX(JMIN,M%IJKW(11,IW)-1)
      CASE(-2) ; JMAX=MIN(JMAX,M%IJKW(11,IW))
      CASE( 3) ; KMIN=MAX(KMIN,M%IJKW(12,IW)-1)
      CASE(-3) ; KMAX=MIN(KMAX,M%IJKW(12,IW))
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
      M%OMESH(NOM)%TMP = TMPA
      ALLOCATE(M%OMESH(NOM)%   H(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
      M%OMESH(NOM)%H = 0.
      ALLOCATE(M%OMESH(NOM)%   U(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
      M%OMESH(NOM)%U = U0
      ALLOCATE(M%OMESH(NOM)%   V(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
      M%OMESH(NOM)%V = V0
      ALLOCATE(M%OMESH(NOM)%   W(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
      M%OMESH(NOM)%W = W0
      IF (NSPEC.GT.0) THEN
      ALLOCATE(M%OMESH(NOM)%  YY(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX,NSPEC))
      ALLOCATE(M%OMESH(NOM)% YYS(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX,NSPEC))
      DO N=1,NSPEC
      M%OMESH(NOM)%YY(:,:,:,N)  = YY0(N)
      M%OMESH(NOM)%YYS(:,:,:,N) = YY0(N)
      ENDDO
      ENDIF
C
C     Wall arrays
C
      ALLOCATE(M%OMESH(NOM)%IJKW(12,M2%NEWC))
      ALLOCATE(M%OMESH(NOM)%IV(0:M2%NEWC))
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
      SUBROUTINE POST_RECEIVES(NM,CODE)
C
      INTEGER, INTENT(IN) :: NM,CODE
C
      N_REQ = 0
C
      OTHER_MESH_LOOP: DO NOM=1,NMESHES
C
      IF (NIC(NM,NOM).EQ.0 .AND. NIC(NOM,NM).EQ.0) CYCLE OTHER_MESH_LOOP
      IF (CODE.GT.0 .AND. .NOT.ACTIVE_MESH(NOM)) CYCLE OTHER_MESH_LOOP
C
      IF (DEBUG) THEN
      WRITE(0,*) NM,' posting receives from ',NOM,' code=',code
      IF (CODE.EQ.0)  
     .WRITE(0,'(A,I2,A,I2,A,I5)') 'NIC(',NM,',',NOM,')=',NIC(NM,NOM)
      ENDIF
C
      M =>MESH(NM)
      M4=>MESH(NOM)
      M3=>MESH(NM)%OMESH(NOM)
C
      RNODE = NOM-1
      TAG   = TAGS(NM,NOM,CODE)
C
      INITIALIZATION_IF: IF (CODE.EQ.0) THEN
C
      IF (NIC(NM,NOM).GT.0) THEN
      ALLOCATE(M3%RPKG1(NIC(NM,NOM)*(3+NSPEC)+1))
      ALLOCATE(M3%RPKG2(NIC(NM,NOM)*(9+NSPEC)+1))
      ALLOCATE(M3%WRPKG((NRA*NSB+1)*NIC(NM,NOM)+1))
      ENDIF
C
      N_REQ = N_REQ+1
      CALL MPI_IRECV(M3%IJKW(1,1),12*M4%NEWC,
     .     MPI_INTEGER,RNODE,TAG,MPI_COMM_WORLD,REQ(N_REQ),IERR)
C
      IF (NIC(NM,NOM).GT.0 .OR. NIC(NOM,NM).GT.0) THEN
      ALLOCATE(M3%R_RDBUF(13*N_DROP_ADOPT_MAX))
      ALLOCATE(M3%R_IDBUF( 2*N_DROP_ADOPT_MAX))
      ALLOCATE(M3%R_LDBUF(   N_DROP_ADOPT_MAX))
      ENDIF
C
      ENDIF INITIALIZATION_IF
C
      PREDICTOR: IF (CODE.EQ.1 .AND. NIC(NM,NOM).GT.0) THEN
      N_REQ = N_REQ+1
      CALL MPI_IRECV(M3%RPKG1(1),NIC(NM,NOM)*(3+NSPEC)+1,
     .     MPI_DOUBLE_PRECISION,RNODE,TAG,MPI_COMM_WORLD,
     .     REQ(N_REQ),IERR)
      ENDIF PREDICTOR
C
      CORRECTOR: IF (CODE.EQ.0 .OR. CODE.EQ.2) THEN
C
      N_REQ = N_REQ+1
      CALL MPI_IRECV(M3%IV(0),M4%NEWC+1,
     .     MPI_INTEGER,RNODE,TAG,MPI_COMM_WORLD,
     .     REQ(N_REQ),IERR)
C
      IF (CODE.EQ.2 .AND. NIC(NM,NOM).GT.0) THEN
      N_REQ=N_REQ+1
      CALL MPI_IRECV(M3%RPKG2(1),NIC(NM,NOM)*(9+NSPEC)+1,
     .     MPI_DOUBLE_PRECISION,RNODE,TAG,MPI_COMM_WORLD,
     .     REQ(N_REQ),IERR)
      ENDIF
      IF (EXCHANGE_RADIATION .AND. NIC(NM,NOM).GT.0 .AND. 
     .    CODE.EQ.2) THEN
      IWW = NIC(NM,NOM)
      N_REQ=N_REQ+1
      CALL MPI_IRECV(M3%WRPKG(1),(NRA*NSB+1)*IWW+1,
     .     MPI_DOUBLE_PRECISION,RNODE,TAG,MPI_COMM_WORLD,
     .     REQ(N_REQ),IERR)
      ENDIF
C
      ENDIF CORRECTOR
C
C     Droplet Orphan Numbers
C
      IF (DROPLET_FILE .AND. 
     .    (NIC(NM,NOM).GT.0 .OR. NIC(NOM,NM).GT.0)) THEN
      N_REQ=N_REQ+1
      CALL MPI_IRECV(M3%N_DROP_ADOPT,
     .     1,MPI_INTEGER,RNODE,TAG,MPI_COMM_WORLD,
     .     REQ(N_REQ),IERR)
      ENDIF
C
C     Droplet Buffer Arrays
C
      IF (DROPLET_FILE .AND. 
     .    (NIC(NM,NOM).GT.0 .OR. NIC(NOM,NM).GT.0)) THEN
      BUFFER_SIZE=13*N_DROP_ADOPT_MAX
      N_REQ=N_REQ+1
      CALL MPI_IRECV(M3%R_RDBUF(1),BUFFER_SIZE,
     .     MPI_DOUBLE_PRECISION,RNODE,TAG,MPI_COMM_WORLD,
     .     REQ(N_REQ),IERR)
      BUFFER_SIZE=2*N_DROP_ADOPT_MAX
      N_REQ=N_REQ+1
      CALL MPI_IRECV(M3%R_IDBUF(1),BUFFER_SIZE,
     .     MPI_INTEGER,RNODE,TAG,MPI_COMM_WORLD,
     .     REQ(N_REQ),IERR)
      BUFFER_SIZE=N_DROP_ADOPT_MAX
      N_REQ=N_REQ+1
      CALL MPI_IRECV(M3%R_LDBUF(1),BUFFER_SIZE,
     .     MPI_LOGICAL,RNODE,TAG,MPI_COMM_WORLD,
     .     REQ(N_REQ),IERR)
      ENDIF
C
      ENDDO OTHER_MESH_LOOP
C
      END SUBROUTINE POST_RECEIVES
C
C
C
      SUBROUTINE MESH_EXCHANGE(CODE)
C
C Exchange Information between Meshes
C
      INTEGER, INTENT(IN) :: CODE
      INTEGER NM,II,JJ,KK,LL,NC,N,NN,SNODE
      INTEGER :: NN1,NN2
C
      TNOW = SECOND()
C
C Send Information to other meshes
C
      NM = MYID+1
C
      SEND_OTHER_MESH_LOOP: DO NOM=1,NMESHES
C
      IF (NIC(NOM,NM).EQ.0 .AND. NIC(NM,NOM).EQ.0) 
     .    CYCLE SEND_OTHER_MESH_LOOP
C
      IF (CODE.GT.0) THEN
      IF (.NOT.ACTIVE_MESH(NM) .OR. .NOT.ACTIVE_MESH(NOM)) 
     .     CYCLE SEND_OTHER_MESH_LOOP
      ENDIF
C
      IF (DEBUG) THEN
      WRITE(0,*) NM,' sending data to ',NOM,' code=',CODE,
     .   ' tag=',TAGS(NM,NOM,CODE)
      ENDIF
C
      M =>MESH(NM)
      M3=>MESH(NM)%OMESH(NOM)
      M4=>MESH(NOM)
C
      SNODE = NOM-1
      RNODE = NM-1
C
      TAG = TAGS(NM,NOM,CODE)
C
      INITIALIZE_SEND_IF: IF (CODE.EQ.0) THEN
C
      IF (NIC(NOM,NM).GT.0) THEN
      ALLOCATE(M3%SPKG1(NIC(NOM,NM)*(3+NSPEC)+1))
      ALLOCATE(M3%SPKG2(NIC(NOM,NM)*(9+NSPEC)+1))
      ENDIF
C
      N_REQ=N_REQ+1
      CALL MPI_ISEND(M%IJKW(1,1),12*M%NEWC,
     .     MPI_INTEGER,SNODE,TAG,MPI_COMM_WORLD,REQ(N_REQ),IERR)
C
      IF (NIC(NOM,NM).GT.0) 
     .   ALLOCATE(M3%WSPKG((NRA*NSB+1)*NIC(NOM,NM)+1))
C
      IF (DROPLET_FILE .AND. 
     .    (NIC(NOM,NM).GT.0 .OR. NIC(NM,NOM).GT.0)) THEN
      ALLOCATE(M3%S_RDBUF(13*N_DROP_ADOPT_MAX))
      ALLOCATE(M3%S_IDBUF( 2*N_DROP_ADOPT_MAX))
      ALLOCATE(M3%S_LDBUF(   N_DROP_ADOPT_MAX))
      ENDIF
C
      ENDIF INITIALIZE_SEND_IF
C
      SEND_PREDICTOR_IF: IF (CODE.EQ.1) THEN
C
      IF (NIC(NOM,NM).GT.0) THEN
      LL = 0
      IWW = 0
      PACK_SPKG1: DO IW=1,M4%NEWC
      IF (M3%IJKW(9,IW).NE.NM .OR. M3%IV(IW).NE.4) CYCLE PACK_SPKG1
      IWW = IWW + 1
      II = M3%IJKW(10,IW)
      JJ = M3%IJKW(11,IW)
      KK = M3%IJKW(12,IW)
      M3%SPKG1(LL+1) = REAL(IW,EB)
      M3%SPKG1(LL+2) = M%TMP(II,JJ,KK)
      M3%SPKG1(LL+3) = M%H(II,JJ,KK)
      IF (NSPEC.GT.0)
     .M3%SPKG1(LL+4:LL+3+NSPEC) = M%YYS(II,JJ,KK,1:NSPEC)
      LL = LL+3+NSPEC
      ENDDO PACK_SPKG1
      M3%SPKG1(IWW*(3+NSPEC)+1) = -999.0_EB
      N_REQ=N_REQ+1
      CALL MPI_ISEND(M3%SPKG1(1),IWW*(3+NSPEC)+1,
     .     MPI_DOUBLE_PRECISION,SNODE,TAG,MPI_COMM_WORLD,
     .     REQ(N_REQ),IERR)
      ENDIF
C
      ENDIF SEND_PREDICTOR_IF
C
      SEND_CORRECTOR_IF: IF (CODE.EQ.0 .OR. CODE.EQ.2) THEN
C
      N_REQ=N_REQ+1
      CALL MPI_ISEND(M%IV(0),M%NEWC+1,
     .     MPI_INTEGER,SNODE,TAG,MPI_COMM_WORLD,
     .     REQ(N_REQ),IERR)
C
      IF (CODE.EQ.2 .AND. NIC(NOM,NM).GT.0) THEN
      LL = 0
      IWW = 0
      PACK_SPKG2: DO IW=1,M4%NEWC
      IF (M3%IJKW(9,IW).NE.NM .OR. M3%IV(IW).NE.4) CYCLE PACK_SPKG2
      IWW = IWW + 1
      II = M3%IJKW(10,IW)
      JJ = M3%IJKW(11,IW)
      KK = M3%IJKW(12,IW)
      M3%SPKG2(LL+1) = REAL(IW,EB)
      M3%SPKG2(LL+2) = M%TMP(II,JJ,KK)
      M3%SPKG2(LL+3) = M%H(II,JJ,KK)
      M3%SPKG2(LL+4) = M%U(II,JJ,KK)
      M3%SPKG2(LL+5) = M%V(II,JJ,KK)
      M3%SPKG2(LL+6) = M%W(II,JJ,KK)
      M3%SPKG2(LL+7) = M%U(II-1,JJ,KK)
      M3%SPKG2(LL+8) = M%V(II,JJ-1,KK)
      M3%SPKG2(LL+9) = M%W(II,JJ,KK-1)
      IF (NSPEC.GT.0)
     .M3%SPKG2(LL+10:LL+9+NSPEC) = M%YY(II,JJ,KK,1:NSPEC)
      LL = LL+9+NSPEC
      ENDDO PACK_SPKG2
      M3%SPKG2(IWW*(9+NSPEC)+1) = -999.0_EB
      N_REQ=N_REQ+1
      CALL MPI_ISEND(M3%SPKG2(1),IWW*(9+NSPEC)+1,
     .     MPI_DOUBLE_PRECISION,SNODE,TAG,MPI_COMM_WORLD,
     .     REQ(N_REQ),IERR)
      ENDIF
C
      SEND_RADIATION: IF ( NIC(NOM,NM).GT.0 .AND. EXCHANGE_RADIATION
     .                     .AND. CODE.EQ.2) THEN
      IWW=0
      LL =0
      PACK_WSPKG: DO IW=1,M4%NEWC
      IF (M3%IJKW(9,IW).NE.NM .OR. M3%IV(IW).NE.4) CYCLE PACK_WSPKG
      IWW = IWW+1
      LL  = LL +1
      M3%WSPKG(LL) = REAL(IW,EB)
      DO NN2=1,NSB
      DO NN1=1,NRA
      LL = LL + 1
      M3%WSPKG(LL) = M3%WALL(IW)%ILW(NN1,NN2)
      ENDDO
      ENDDO
      ENDDO PACK_WSPKG
      M3%WSPKG(LL+1) = -999.0_EB
      N_REQ=N_REQ+1
      CALL MPI_ISEND(M3%WSPKG(1),(NRA*NSB+1)*IWW+1,
     .     MPI_DOUBLE_PRECISION,SNODE,TAG,MPI_COMM_WORLD,
     .     REQ(N_REQ),IERR)
      ENDIF SEND_RADIATION
C
      ENDIF SEND_CORRECTOR_IF
C
C     Get Number of Droplet Orphans
C
      IF (DROPLET_FILE) THEN
C
      N_REQ=N_REQ+1
      CALL MPI_ISEND(M3%N_DROP_ORPHANS,
     .     1,MPI_INTEGER,SNODE,TAG,MPI_COMM_WORLD,REQ(N_REQ),IERR)
C
      ENDIF
C
C     Sending/Receiving Droplet Buffer Arrays
C
      IF_SEND_DROPLETS: IF (DROPLET_FILE) THEN 
C
      NC = 13
      DO N=1,MIN(M3%N_DROP_ORPHANS,N_DROP_ADOPT_MAX)
      M3%S_RDBUF((N-1)*NC+1)  = M3%DROPLET(N)%X
      M3%S_RDBUF((N-1)*NC+2)  = M3%DROPLET(N)%Y
      M3%S_RDBUF((N-1)*NC+3)  = M3%DROPLET(N)%Z
      M3%S_RDBUF((N-1)*NC+4)  = M3%DROPLET(N)%TMP
      M3%S_RDBUF((N-1)*NC+5)  = M3%DROPLET(N)%U
      M3%S_RDBUF((N-1)*NC+6)  = M3%DROPLET(N)%V
      M3%S_RDBUF((N-1)*NC+7)  = M3%DROPLET(N)%W
      M3%S_RDBUF((N-1)*NC+8)  = M3%DROPLET(N)%R
      M3%S_RDBUF((N-1)*NC+9)  = M3%DROPLET(N)%PWT
      M3%S_RDBUF((N-1)*NC+10) = M3%DROPLET(N)%A_X
      M3%S_RDBUF((N-1)*NC+11) = M3%DROPLET(N)%A_Y
      M3%S_RDBUF((N-1)*NC+12) = M3%DROPLET(N)%A_Z
      M3%S_RDBUF((N-1)*NC+13) = M3%DROPLET(N)%T
      ENDDO
      BUFFER_SIZE = NC*MIN(M3%N_DROP_ORPHANS,N_DROP_ADOPT_MAX)
      BUFFER_SIZE = MAX(1,BUFFER_SIZE)
      N_REQ=N_REQ+1
      CALL MPI_ISEND(M3%S_RDBUF(1),BUFFER_SIZE,
     .     MPI_DOUBLE_PRECISION,SNODE,TAG,MPI_COMM_WORLD,
     .     REQ(N_REQ),IERR)
C
      NC = 2
      DO N=1,MIN(M3%N_DROP_ORPHANS,N_DROP_ADOPT_MAX)
      M3%S_IDBUF((N-1)*NC+1) = M3%DROPLET(N)%IOR
      M3%S_IDBUF((N-1)*NC+2) = M3%DROPLET(N)%CLASS
      ENDDO
      BUFFER_SIZE = NC*MIN(M3%N_DROP_ORPHANS,N_DROP_ADOPT_MAX)
      BUFFER_SIZE = MAX(1,BUFFER_SIZE)
      N_REQ=N_REQ+1
      CALL MPI_ISEND(M3%S_IDBUF(1),BUFFER_SIZE,
     .     MPI_INTEGER,SNODE,TAG,MPI_COMM_WORLD,REQ(N_REQ),IERR)
C
      NC = 1
      DO N=1,MIN(M3%N_DROP_ORPHANS,N_DROP_ADOPT_MAX)
      M3%S_LDBUF((N-1)*NC+1) = M3%DROPLET(N)%SHOW
      ENDDO
      BUFFER_SIZE = NC*MIN(M3%N_DROP_ORPHANS,N_DROP_ADOPT_MAX)
      BUFFER_SIZE = MAX(1,BUFFER_SIZE)
      N_REQ=N_REQ+1
      CALL MPI_ISEND(M3%S_LDBUF(1),BUFFER_SIZE,
     .     MPI_LOGICAL,SNODE,TAG,MPI_COMM_WORLD,REQ(N_REQ),IERR)
C
      M3%N_DROP_ORPHANS = 0
C
      ENDIF IF_SEND_DROPLETS
C
      ENDDO SEND_OTHER_MESH_LOOP
C
C Receive Messages (NOM is the receiver, NM is the sender)
C
      CALL MPI_WAITALL(N_REQ,REQ(1:N_REQ),ARRAY_OF_STATUSES,IERR)
C
      NOM = MYID+1
C
      MESH_LOOP: DO NM=1,NMESHES
C
      IF (NIC(NOM,NM).EQ.0 .AND. NIC(NM,NOM).EQ.0) CYCLE MESH_LOOP
C
      IF (CODE.GT.0) THEN
      IF (.NOT.ACTIVE_MESH(NM) .OR. .NOT.ACTIVE_MESH(NOM)) 
     .     CYCLE MESH_LOOP
      ENDIF
C
      IF (DEBUG) THEN
      WRITE(0,*) NOM,' receiving data from ',NM,
     .  ' code=',CODE,' tag=',TAGS(NM,NOM,CODE)
      ENDIF
C
      M =>MESH(NM)
      M2=>MESH(NOM)%OMESH(NM)
      M4=>MESH(NOM)
C
      SNODE = NOM-1
      RNODE = NM-1
C
      TAG = TAGS(NM,NOM,CODE)
C
      INITIALIZE_RECEIVE_IF: IF (CODE.EQ.0) THEN
C
      IF (RADIATION) THEN
         DO IW=1,M%NEWC
         IF (M2%IJKW(9,IW).EQ.NOM) THEN
            ALLOCATE(M2%WALL(IW)%ILW(NRA,NSB))
            M2%WALL(IW)%ILW = SIGMA*TMPA4*RPI
            ENDIF
         ENDDO
      ENDIF
C
      ENDIF INITIALIZE_RECEIVE_IF
C
      RECEIVE_PREDICTOR_IF: IF (CODE.EQ.1) THEN
C
      IF (NIC(NOM,NM).GT.0) THEN
C
      LL = 0
      UNPACK_RPKG1: DO 
      IW = NINT(M2%RPKG1(LL+1))
      IF (IW.EQ.-999) EXIT UNPACK_RPKG1
      II = M4%IJKW(10,IW)
      JJ = M4%IJKW(11,IW)
      KK = M4%IJKW(12,IW)
      M2%TMP(II,JJ,KK) = M2%RPKG1(LL+2)
      M2%H(II,JJ,KK)   = M2%RPKG1(LL+3)
      IF (NSPEC.GT.0) M2%YYS(II,JJ,KK,1:NSPEC)=M2%RPKG1(LL+4:LL+3+NSPEC)
      LL = LL+3+NSPEC
      ENDDO UNPACK_RPKG1
      ENDIF
C
      ENDIF RECEIVE_PREDICTOR_IF
C
      RECEIVE_CORRECTOR_IF: IF (CODE.EQ.0 .OR. CODE.EQ.2) THEN
C
      IF (CODE.EQ.2 .AND. NIC(NOM,NM).GT.0) THEN
      LL = 0
      UNPACK_RPKG2: DO 
      IW = NINT(M2%RPKG2(LL+1))
      IF (IW.EQ.-999) EXIT UNPACK_RPKG2
      II = M4%IJKW(10,IW)
      JJ = M4%IJKW(11,IW)
      KK = M4%IJKW(12,IW)
      M2%TMP(II,JJ,KK) = M2%RPKG2(LL+2)
      M2%H(II,JJ,KK)   = M2%RPKG2(LL+3)
      M2%U(II,JJ,KK)   = M2%RPKG2(LL+4)
      M2%V(II,JJ,KK)   = M2%RPKG2(LL+5)
      M2%W(II,JJ,KK)   = M2%RPKG2(LL+6)
      M2%U(II-1,JJ,KK) = M2%RPKG2(LL+7)
      M2%V(II,JJ-1,KK) = M2%RPKG2(LL+8)
      M2%W(II,JJ,KK-1) = M2%RPKG2(LL+9)
      IF (NSPEC.GT.0) M2%YY(II,JJ,KK,1:NSPEC)=M2%RPKG2(LL+10:LL+9+NSPEC)
      LL = LL+9+NSPEC
      ENDDO UNPACK_RPKG2
      ENDIF
C
      RECEIVE_RADIATION: IF ( NIC(NOM,NM).GT.0 .AND. EXCHANGE_RADIATION 
     .                        .AND. CODE.EQ.2) THEN
      LL =0
      UNPACK_WRPKG: DO 
      LL  = LL + 1
      IW = NINT(M2%WRPKG(LL))
      IF (IW.EQ.-999) EXIT UNPACK_WRPKG
      DO NN2=1,NSB
      DO NN1=1,NRA
      LL = LL + 1
      M4%WALL(IW)%ILW(NN1,NN2) = M2%WRPKG(LL)
      ENDDO
      ENDDO
      ENDDO UNPACK_WRPKG
      ENDIF RECEIVE_RADIATION
C
      ENDIF RECEIVE_CORRECTOR_IF
C
C     Get Number of Droplet Orphans
C
      IF (DROPLET_FILE) THEN
C
      M2%N_DROP_ADOPT = MIN(M2%N_DROP_ADOPT,N_DROP_ADOPT_MAX)
      IF (M4%NLP+M2%N_DROP_ADOPT.GT.M4%NLPDIM) THEN
      CALL RE_ALLOCATE_DROPLETS(1,NOM,0)
      ENDIF
C
      ENDIF
C
C     Sending/Receiving Droplet Buffer Arrays
C
      IF_RECEIVE_DROPLETS: IF (DROPLET_FILE) THEN 
C
      IF_NO_DROPLETS_SENT: IF (M2%N_DROP_ADOPT.EQ.0) THEN
C
      ELSE
C
      NC = 13
      DO N=M4%NLP+1,M4%NLP+M2%N_DROP_ADOPT
      NN = N-M4%NLP-1
      M4%DROPLET(N)%X   = M2%R_RDBUF((NN)*NC+1) 
      M4%DROPLET(N)%Y   = M2%R_RDBUF((NN)*NC+2) 
      M4%DROPLET(N)%Z   = M2%R_RDBUF((NN)*NC+3) 
      M4%DROPLET(N)%TMP = M2%R_RDBUF((NN)*NC+4) 
      M4%DROPLET(N)%U   = M2%R_RDBUF((NN)*NC+5) 
      M4%DROPLET(N)%V   = M2%R_RDBUF((NN)*NC+6) 
      M4%DROPLET(N)%W   = M2%R_RDBUF((NN)*NC+7) 
      M4%DROPLET(N)%R   = M2%R_RDBUF((NN)*NC+8) 
      M4%DROPLET(N)%PWT = M2%R_RDBUF((NN)*NC+9) 
      M4%DROPLET(N)%A_X = M2%R_RDBUF((NN)*NC+10) 
      M4%DROPLET(N)%A_Y = M2%R_RDBUF((NN)*NC+11) 
      M4%DROPLET(N)%A_Z = M2%R_RDBUF((NN)*NC+12) 
      M4%DROPLET(N)%T   = M2%R_RDBUF((NN)*NC+13) 
      ENDDO
C
      NC = 2
      DO N=M4%NLP+1,M4%NLP+M2%N_DROP_ADOPT
      NN = N-M4%NLP-1
      M4%DROPLET(N)%IOR    = M2%R_IDBUF((NN)*NC+1) 
      M4%DROPLET(N)%CLASS  = M2%R_IDBUF((NN)*NC+2) 
      ENDDO
C
      NC = 1
      DO N=M4%NLP+1,M4%NLP+M2%N_DROP_ADOPT
      NN = N-M4%NLP-1
      M4%DROPLET(N)%SHOW = M2%R_LDBUF((NN)*NC+1)  
      ENDDO
C
      M4%NLP = M4%NLP + M2%N_DROP_ADOPT
C
      ENDIF IF_NO_DROPLETS_SENT
C
      ENDIF IF_RECEIVE_DROPLETS
C
      ENDDO MESH_LOOP
C
      TUSED(11,:)=TUSED(11,:) + SECOND() - TNOW
      END SUBROUTINE MESH_EXCHANGE
C
C
      SUBROUTINE WRITE_STRINGS
C
C     Write character strings out to the .smv file
C
      INTEGER :: N,NOM,N_STRINGS_DUM
      CHARACTER(50), ALLOCATABLE, DIMENSION(:) :: STRING_DUM
C
C     All meshes send their STRINGs to node 0
C
      IF (MYID.GT.0) THEN
      CALL MPI_SEND(MESH(MYID+1)%N_STRINGS,1,MPI_INTEGER,
     .              0,1,MPI_COMM_WORLD,IERR)
      IF (MESH(MYID+1)%N_STRINGS.GT.0) 
     .CALL MPI_SEND(MESH(MYID+1)%STRING(1),MESH(MYID+1)%N_STRINGS*50,
     .              MPI_CHARACTER,0,MYID,MPI_COMM_WORLD,IERR)
      ENDIF
C
C     Node 0 receives the STRINGs and writes them to the .smv file
C
      IF (MYID.EQ.0) THEN
         DO N=1,MESH(1)%N_STRINGS
         WRITE(LU4,'(A)') TRIM(MESH(1)%STRING(N))
         ENDDO
         OTHER_MESH_LOOP: DO NOM=2,NMESHES 
         CALL MPI_RECV(N_STRINGS_DUM,1,MPI_INTEGER,
     .                 NOM-1,1,MPI_COMM_WORLD,STATUS,IERR)
         IF (N_STRINGS_DUM.GT.0) THEN
         ALLOCATE(STRING_DUM(N_STRINGS_DUM))
         CALL MPI_RECV(STRING_DUM(1),N_STRINGS_DUM*50,
     .                 MPI_CHARACTER,NOM-1,NOM-1,
     .                 MPI_COMM_WORLD,STATUS,IERR)
         DO N=1,N_STRINGS_DUM
         WRITE(LU4,'(A)') TRIM(STRING_DUM(N))
         ENDDO
         DEALLOCATE(STRING_DUM)
         ENDIF
         ENDDO OTHER_MESH_LOOP
      ENDIF
C
C     All STRING arrays are zeroed out
C
      MESH(MYID+1)%N_STRINGS = 0
C
      END SUBROUTINE WRITE_STRINGS
C
C
      SUBROUTINE EXCHANGE_DIAGNOSTICS
C
      INTEGER  :: NOM,NM,NECYC
      REAL(EB) :: T_SUM
C
      TNOW = SECOND()
C
      NM = MYID+1

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
      CALL MPI_GATHER(T_ACCUM(MYID+1),1,MPI_DOUBLE_PRECISION,
     .           T_ACCUM,1,MPI_DOUBLE_PRECISION,0,
     .           MPI_COMM_WORLD,IERR)
      CALL MPI_GATHER(T_PER_STEP(MYID+1),1,MPI_DOUBLE_PRECISION,
     .           T_PER_STEP,1,MPI_DOUBLE_PRECISION,0,
     .           MPI_COMM_WORLD,IERR)
      CALL MPI_GATHER(NTCYC(MYID+1),1,MPI_INTEGER,
     .           NTCYC,1,MPI_INTEGER,0,
     .           MPI_COMM_WORLD,IERR)
      CALL MPI_GATHER(HRR(MYID+1),1,MPI_DOUBLE_PRECISION,
     .           HRR,1,MPI_DOUBLE_PRECISION,0,
     .           MPI_COMM_WORLD,IERR)
      CALL MPI_GATHER(RHRR(MYID+1),1,MPI_DOUBLE_PRECISION,
     .           RHRR,1,MPI_DOUBLE_PRECISION,0,
     .           MPI_COMM_WORLD,IERR)
C
      IF (MYID.GT.0) THEN
C
         CALL MPI_SEND(MESH(MYID+1)%DT,1,MPI_DOUBLE_PRECISION,
     .                 0,1,MPI_COMM_WORLD,IERR)
         CALL MPI_SEND(MESH(MYID+1)%CFL,1,MPI_DOUBLE_PRECISION,
     .                 0,1,MPI_COMM_WORLD,IERR)
         CALL MPI_SEND(MESH(MYID+1)%DIVMX,1,MPI_DOUBLE_PRECISION,
     .                 0,1,MPI_COMM_WORLD,IERR)
         CALL MPI_SEND(MESH(MYID+1)%DIVMN,1,MPI_DOUBLE_PRECISION,
     .                 0,1,MPI_COMM_WORLD,IERR)
         CALL MPI_SEND(MESH(MYID+1)%RESMAX,1,MPI_DOUBLE_PRECISION,
     .                 0,1,MPI_COMM_WORLD,IERR)
         CALL MPI_SEND(MESH(MYID+1)%POIS_PTB,1,MPI_DOUBLE_PRECISION,
     .                 0,1,MPI_COMM_WORLD,IERR)
         CALL MPI_SEND(MESH(MYID+1)%POIS_ERR,1,MPI_DOUBLE_PRECISION,
     .                 0,1,MPI_COMM_WORLD,IERR)
         CALL MPI_SEND(MESH(MYID+1)%VN,1,MPI_DOUBLE_PRECISION,
     .                 0,1,MPI_COMM_WORLD,IERR)
         CALL MPI_SEND(MESH(MYID+1)%P0,1,MPI_DOUBLE_PRECISION,
     .                 0,1,MPI_COMM_WORLD,IERR)
         CALL MPI_SEND(MESH(MYID+1)%Z_F_EFF,1,MPI_DOUBLE_PRECISION,
     .                 0,1,MPI_COMM_WORLD,IERR)
         CALL MPI_SEND(MESH(MYID+1)%ICFL,1,MPI_INTEGER,
     .                 0,1,MPI_COMM_WORLD,IERR)
         CALL MPI_SEND(MESH(MYID+1)%JCFL,1,MPI_INTEGER,
     .                 0,1,MPI_COMM_WORLD,IERR)
         CALL MPI_SEND(MESH(MYID+1)%KCFL,1,MPI_INTEGER,
     .                 0,1,MPI_COMM_WORLD,IERR)
         CALL MPI_SEND(MESH(MYID+1)%IMX,1,MPI_INTEGER,
     .                 0,1,MPI_COMM_WORLD,IERR)
         CALL MPI_SEND(MESH(MYID+1)%JMX,1,MPI_INTEGER,
     .                 0,1,MPI_COMM_WORLD,IERR)
         CALL MPI_SEND(MESH(MYID+1)%KMX,1,MPI_INTEGER,
     .                 0,1,MPI_COMM_WORLD,IERR)
         CALL MPI_SEND(MESH(MYID+1)%IMN,1,MPI_INTEGER,
     .                 0,1,MPI_COMM_WORLD,IERR)
         CALL MPI_SEND(MESH(MYID+1)%JMN,1,MPI_INTEGER,
     .                 0,1,MPI_COMM_WORLD,IERR)
         CALL MPI_SEND(MESH(MYID+1)%KMN,1,MPI_INTEGER,
     .                 0,1,MPI_COMM_WORLD,IERR)
         CALL MPI_SEND(MESH(MYID+1)%IRM,1,MPI_INTEGER,
     .                 0,1,MPI_COMM_WORLD,IERR)
         CALL MPI_SEND(MESH(MYID+1)%JRM,1,MPI_INTEGER,
     .                 0,1,MPI_COMM_WORLD,IERR)
         CALL MPI_SEND(MESH(MYID+1)%KRM,1,MPI_INTEGER,
     .                 0,1,MPI_COMM_WORLD,IERR)
         CALL MPI_SEND(MESH(MYID+1)%I_VN,1,MPI_INTEGER,
     .                 0,1,MPI_COMM_WORLD,IERR)
         CALL MPI_SEND(MESH(MYID+1)%J_VN,1,MPI_INTEGER,
     .                 0,1,MPI_COMM_WORLD,IERR)
         CALL MPI_SEND(MESH(MYID+1)%K_VN,1,MPI_INTEGER,
     .                 0,1,MPI_COMM_WORLD,IERR)
         CALL MPI_SEND(MESH(MYID+1)%NLP,1,MPI_INTEGER,
     .                 0,1,MPI_COMM_WORLD,IERR)
C
      ELSE 
C
         DO NOM=2,NMESHES
         CALL MPI_RECV(MESH(NOM)%DT,1,MPI_DOUBLE_PRECISION,
     .                 NOM-1,1,MPI_COMM_WORLD,STATUS,IERR)
         CALL MPI_RECV(MESH(NOM)%CFL,1,MPI_DOUBLE_PRECISION,
     .                 NOM-1,1,MPI_COMM_WORLD,STATUS,IERR)
         CALL MPI_RECV(MESH(NOM)%DIVMX,1,MPI_DOUBLE_PRECISION,
     .                 NOM-1,1,MPI_COMM_WORLD,STATUS,IERR)
         CALL MPI_RECV(MESH(NOM)%DIVMN,1,MPI_DOUBLE_PRECISION,
     .                 NOM-1,1,MPI_COMM_WORLD,STATUS,IERR)
         CALL MPI_RECV(MESH(NOM)%RESMAX,1,MPI_DOUBLE_PRECISION,
     .                 NOM-1,1,MPI_COMM_WORLD,STATUS,IERR)
         CALL MPI_RECV(MESH(NOM)%POIS_PTB,1,MPI_DOUBLE_PRECISION,
     .                 NOM-1,1,MPI_COMM_WORLD,STATUS,IERR)
         CALL MPI_RECV(MESH(NOM)%POIS_ERR,1,MPI_DOUBLE_PRECISION,
     .                 NOM-1,1,MPI_COMM_WORLD,STATUS,IERR)
         CALL MPI_RECV(MESH(NOM)%VN,1,MPI_DOUBLE_PRECISION,
     .                 NOM-1,1,MPI_COMM_WORLD,STATUS,IERR)
         CALL MPI_RECV(MESH(NOM)%P0,1,MPI_DOUBLE_PRECISION,
     .                 NOM-1,1,MPI_COMM_WORLD,STATUS,IERR)
         CALL MPI_RECV(MESH(NOM)%Z_F_EFF,1,MPI_DOUBLE_PRECISION,
     .                 NOM-1,1,MPI_COMM_WORLD,STATUS,IERR)
         CALL MPI_RECV(MESH(NOM)%ICFL,1,MPI_INTEGER,
     .                 NOM-1,1,MPI_COMM_WORLD,STATUS,IERR)
         CALL MPI_RECV(MESH(NOM)%JCFL,1,MPI_INTEGER,
     .                 NOM-1,1,MPI_COMM_WORLD,STATUS,IERR)
         CALL MPI_RECV(MESH(NOM)%KCFL,1,MPI_INTEGER,
     .                 NOM-1,1,MPI_COMM_WORLD,STATUS,IERR)
         CALL MPI_RECV(MESH(NOM)%IMX,1,MPI_INTEGER,
     .                 NOM-1,1,MPI_COMM_WORLD,STATUS,IERR)
         CALL MPI_RECV(MESH(NOM)%JMX,1,MPI_INTEGER,
     .                 NOM-1,1,MPI_COMM_WORLD,STATUS,IERR)
         CALL MPI_RECV(MESH(NOM)%KMX,1,MPI_INTEGER,
     .                 NOM-1,1,MPI_COMM_WORLD,STATUS,IERR)
         CALL MPI_RECV(MESH(NOM)%IMN,1,MPI_INTEGER,
     .                 NOM-1,1,MPI_COMM_WORLD,STATUS,IERR)
         CALL MPI_RECV(MESH(NOM)%JMN,1,MPI_INTEGER,
     .                 NOM-1,1,MPI_COMM_WORLD,STATUS,IERR)
         CALL MPI_RECV(MESH(NOM)%KMN,1,MPI_INTEGER,
     .                 NOM-1,1,MPI_COMM_WORLD,STATUS,IERR)
         CALL MPI_RECV(MESH(NOM)%IRM,1,MPI_INTEGER,
     .                 NOM-1,1,MPI_COMM_WORLD,STATUS,IERR)
         CALL MPI_RECV(MESH(NOM)%JRM,1,MPI_INTEGER,
     .                 NOM-1,1,MPI_COMM_WORLD,STATUS,IERR)
         CALL MPI_RECV(MESH(NOM)%KRM,1,MPI_INTEGER,
     .                 NOM-1,1,MPI_COMM_WORLD,STATUS,IERR)
         CALL MPI_RECV(MESH(NOM)%I_VN,1,MPI_INTEGER,
     .                 NOM-1,1,MPI_COMM_WORLD,STATUS,IERR)
         CALL MPI_RECV(MESH(NOM)%J_VN,1,MPI_INTEGER,
     .                 NOM-1,1,MPI_COMM_WORLD,STATUS,IERR)
         CALL MPI_RECV(MESH(NOM)%K_VN,1,MPI_INTEGER,
     .                 NOM-1,1,MPI_COMM_WORLD,STATUS,IERR)
         CALL MPI_RECV(MESH(NOM)%NLP,1,MPI_INTEGER,
     .                 NOM-1,1,MPI_COMM_WORLD,STATUS,IERR)
         ENDDO
C
      ENDIF
C
      TUSED(11,:) = TUSED(11,:) + SECOND() - TNOW
      END SUBROUTINE EXCHANGE_DIAGNOSTICS
C
C
      END PROGRAM FDS
