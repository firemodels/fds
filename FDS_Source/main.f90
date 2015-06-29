PROGRAM FDS  

! Fire Dynamics Simulator, Main Program, Multiple CPU version.

USE PRECISION_PARAMETERS
USE MESH_VARIABLES
USE GLOBAL_CONSTANTS
USE TRAN
USE DUMP
USE READ_INPUT
USE INIT
USE DIVG
USE PRES
USE MASS
USE PART
USE VEGE
USE VELO
USE RAD
USE RADCONS, ONLY : DLN
USE OUTPUT_DATA
USE MEMORY_FUNCTIONS
USE HVAC_ROUTINES
USE COMP_FUNCTIONS, ONLY : SECOND, WALL_CLOCK_TIME
USE DEVICE_VARIABLES
USE WALL_ROUTINES
USE FIRE
USE CONTROL_FUNCTIONS
USE EVAC
USE TURBULENCE, ONLY: NS_ANALYTICAL_SOLUTION,INIT_TURB_ARRAYS,COMPRESSION_WAVE,TWOD_VORTEX_CERFACS,TWOD_VORTEX_UMD, &
                      SYNTHETIC_TURBULENCE,SYNTHETIC_EDDY_SETUP,SANDIA_DAT,GET_REV_turb
USE EMBEDDED_MESH_METHOD, ONLY: SCALARF_EMB,VELOCITY_EMB,RESTRICT_MASS_EMB,RESTRICT_DIV_EMB,SCALAR_GHOST_EMB, &
                                PROJECT_VELOCITY,SORT_MESH_LEVEL,MATCH_VELOCITY_EMB,GET_REV_samr
USE MANUFACTURED_SOLUTIONS, ONLY: SHUNN_MMS_3
USE COMPLEX_GEOMETRY, ONLY: INIT_IBM,GET_REV_geom
USE OPENMP
USE MPI
USE SCRC, ONLY: SCARC_SETUP, SCARC_SOLVER, SCARC_TIMINGS, GET_REV_SCRC
USE BOXTETRA_ROUTINES, ONLY: GET_REV_gsmv
USE SOOT_ROUTINES, ONLY: CALC_AGGLOMERATION

! dummy comment to force svn change - Thu Jun 25 21:56:37 EDT 2015

IMPLICIT NONE

! Miscellaneous declarations

CHARACTER(255), PARAMETER :: mainid='$Id$'
CHARACTER(255), PARAMETER :: mainrev='$Revision$'
CHARACTER(255), PARAMETER :: maindate='$Date$'


LOGICAL  :: EX=.FALSE.,DIAGNOSTICS,EXCHANGE_EVACUATION=.FALSE.,FIRST_PASS,CALL_UPDATE_CONTROLS,CTRL_STOP_STATUS,FOUND
INTEGER  :: LO10,NM,IZERO,CNT,ANG_INC_COUNTER
REAL(EB) :: T_MAX,T_MIN
REAL(EB), ALLOCATABLE, DIMENSION(:) ::  T,TC_GLB,TC_LOC,DT_SYNC,DT_NEXT_SYNC,TI_LOC,TI_GLB, &
                                        DSUM_ALL,PSUM_ALL,USUM_ALL,DSUM_ALL_LOCAL,PSUM_ALL_LOCAL,USUM_ALL_LOCAL
LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: CONNECTED_ZONES_GLOBAL,CONNECTED_ZONES_LOCAL,CONNECTED_MESHES
LOGICAL, ALLOCATABLE, DIMENSION(:) ::  ACTIVE_MESH,STATE_GLB,STATE_LOC
INTEGER :: NOM,IWW,IW !,IERROR
TYPE (MESH_TYPE), POINTER :: M,M4
TYPE (OMESH_TYPE), POINTER :: M2,M3,M5
 
! MPI stuff

INTEGER :: N,I,IERR=0,STATUS(MPI_STATUS_SIZE)
INTEGER :: PNAMELEN=0,TAG_EVAC
INTEGER :: PROVIDED
INTEGER, PARAMETER :: REQUIRED=MPI_THREAD_SINGLE
INTEGER, ALLOCATABLE, DIMENSION(:) :: REQ,REQ1,REQ2,REQ3,REQ4,REQ5,REQ6,REQ7,REQ8,REQ9,COUNTS,DISPLS,&
                                      COUNTS2D,DISPLS2D,COUNTS_TIMERS,DISPLS_TIMERS, &
                                      COUNTS_MASS,DISPLS_MASS,COUNTS_HVAC,DISPLS_HVAC,COUNTS_Q_DOT,DISPLS_Q_DOT, &
                                      COUNTS_M_DOT,DISPLS_M_DOT,COUNTS_HVAC_SPECIES,DISPLS_HVAC_SPECIES
INTEGER :: N_REQ,N_REQ1,N_REQ2,N_REQ3,N_REQ4,N_REQ5,N_REQ6,N_REQ7,N_REQ8,N_REQ9,N_COMMUNICATIONS
CHARACTER(MPI_MAX_PROCESSOR_NAME) :: PNAME
INTEGER, ALLOCATABLE, DIMENSION(:)        :: INTEGER_BUFFER_1
INTEGER, ALLOCATABLE, DIMENSION(:,:)      :: INTEGER_BUFFER_2
REAL(EB), ALLOCATABLE, DIMENSION(:)       :: REAL_BUFFER_1
REAL(EB), ALLOCATABLE, DIMENSION(:,:)     :: REAL_BUFFER_3,REAL_BUFFER_5,REAL_BUFFER_6,REAL_BUFFER_10,REAL_BUFFER_11,REAL_BUFFER_12
REAL(EB), ALLOCATABLE, DIMENSION(:,:,:)   :: REAL_BUFFER_7,REAL_BUFFER_8
REAL(EB), ALLOCATABLE, DIMENSION(:,:,:,:) :: REAL_BUFFER_9
LOGICAL, ALLOCATABLE, DIMENSION(:)        :: LOGICAL_BUFFER_1
 
! Initialize MPI (First executable lines of code)
 
CALL MPI_INIT_THREAD(REQUIRED,PROVIDED,IERR)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID, IERR)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, N_MPI_PROCESSES, IERR)
CALL MPI_GET_PROCESSOR_NAME(PNAME, PNAMELEN, IERR)

IF (PNAME/='null') USE_MPI = .TRUE.
 
! Initialize OpenMP

CALL OPENMP_CHECK

! Check that MPI processes and OpenMP threads are working properly

CALL CHECK_MPI_OPENMP

! Start wall clock timing

INITIALIZATION_PHASE = .TRUE.
WALL_CLOCK_START = WALL_CLOCK_TIME()
 
! Assign a compilation date (All Nodes)

WRITE(VERSION_STRING,'(A)') 'FDS 6.2.0'

CALL GET_INFO (REVISION,REVISION_DATE,COMPILE_DATE)

! Read input from CHID.fds file and stop the code if any errors are found
 
CALL READ_DATA

CALL STOP_CHECK(1)

! Set up send and receive buffer counts and displacements

CALL MPI_INITIALIZATION_CHORES(1)
 
! Open and write to Smokeview and status file (Master Node Only)
 
CALL ASSIGN_FILE_NAMES

DO N=0,N_MPI_PROCESSES-1
   IF (MYID==N) CALL WRITE_SMOKEVIEW_FILE
   IF (N==N_MPI_PROCESSES-1) EXIT
   CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
ENDDO
CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
IF (MYID==0) THEN
   OPEN(LU_SMV,FILE=FN_SMV,FORM='FORMATTED', STATUS='OLD',POSITION='APPEND')
   CALL WRITE_STATUS_FILES
ENDIF

! Stop all the processes if this is just a set-up run
 
IF (SET_UP_ONLY) THEN
   IF (MYID==0) CALL INITIALIZE_DIAGNOSTIC_FILE
   STOP_STATUS = SETUP_ONLY_STOP
   CALL STOP_CHECK(1)
ENDIF

! Allocate various utility arrays
 
CALL MPI_INITIALIZATION_CHORES(2)

! Start the clock

T = T_BEGIN
 
! Initialize global parameters
 
CALL INITIALIZE_GLOBAL_VARIABLES
CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

! Initialize radiation 
 
IF (RADIATION) CALL INIT_RADIATION

! Allocate and initialize mesh-specific variables, and check to see if the code should stop
 
DO NM=1,NMESHES
   IF (PROCESS(NM)==MYID) CALL INITIALIZE_MESH_VARIABLES_1(NM)
ENDDO
CALL STOP_CHECK(1)
 
! Allocate and initialize OMESH arrays to hold "other mesh" data for a given mesh
 
N_COMMUNICATIONS = 0

DO NM=1,NMESHES
   IF (PROCESS(NM)==MYID) CALL INITIALIZE_MESH_EXCHANGE_1(NM)
ENDDO
CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
 
! Allocate "request" arrays to keep track of MPI communications

CALL MPI_INITIALIZATION_CHORES(3)

! Exchange information related to size of OMESH arrays

CALL MPI_INITIALIZATION_CHORES(4)

! Allocate and initialize OMESH arrays to hold "other mesh" data for a given mesh

DO NM=1,NMESHES
   IF (PROCESS(NM)==MYID) CALL INITIALIZE_MESH_EXCHANGE_2(NM)
ENDDO
CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

! Exchange CELL_COUNT, the dimension of various arrays related to obstructions

IF (N_MPI_PROCESSES>1) THEN
   CALL MPI_ALLGATHERV(MPI_IN_PLACE,COUNTS(MYID),MPI_INTEGER,CELL_COUNT,COUNTS,DISPLS,MPI_INTEGER,MPI_COMM_WORLD,IERR)
   IF (IERR/=MPI_SUCCESS) CALL HANDLE_MPI_ERROR('Error in MPI_ALLGATHERV for CELL_COUNT',IERR)
ENDIF

! Initialize persistent MPI sends and receives and allocate buffer arrays.

N_REQ1 = 0
N_REQ2 = 0
N_REQ3 = 0
N_REQ4 = 0
N_REQ5 = 0
N_REQ6 = 0
N_REQ7 = 0
N_REQ8 = 0
N_REQ9 = 0

CALL POST_RECEIVES(0)
CALL MESH_EXCHANGE(0)
CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

! Finish initializing mesh variables

DO NM=1,NMESHES
   IF (PROCESS(NM)==MYID) CALL INITIALIZE_MESH_VARIABLES_2(NM)
ENDDO

! Create arrays and communicators to exchange back wall information across mesh boundaries

CALL INITIALIZE_BACK_WALL_EXCHANGE

CALL STOP_CHECK(1)

! Initialize ScaRC solver

IF (PRES_METHOD == 'SCARC') CALL SCARC_SETUP

! Initialize turb arrays

DO NM=1,NMESHES
   IF (PROCESS(NM)/=MYID .OR. TGA_SURF_INDEX>0) CYCLE
   CALL INIT_TURB_ARRAYS(NM)
ENDDO

! Initialize unstructured geometry

IF (N_FACE>0) THEN
   DO NM=1,NMESHES
      IF (PROCESS(NM)/=MYID) CYCLE
      CALL INIT_IBM(0._EB,NM)
   ENDDO
ENDIF

! Initialize the flow field with random noise to eliminate false symmetries

DO NM=1,NMESHES
   IF (PROCESS(NM)/=MYID .OR. TGA_SURF_INDEX>0) CYCLE
   IF (NOISE) CALL INITIAL_NOISE(NM)
   IF (PERIODIC_TEST==1) CALL NS_ANALYTICAL_SOLUTION(NM,T_BEGIN,RK_STAGE=2)
   IF (PERIODIC_TEST==2) CALL UVW_INIT(NM,UVW_FILE)
   IF (PERIODIC_TEST==3) CALL COMPRESSION_WAVE(NM,0._EB,3)
   IF (PERIODIC_TEST==4) CALL COMPRESSION_WAVE(NM,0._EB,4)
   IF (PERIODIC_TEST==6) CALL TWOD_VORTEX_CERFACS(NM)
   IF (PERIODIC_TEST==7) CALL SHUNN_MMS_3(NM)
   IF (PERIODIC_TEST==8) CALL NS_ANALYTICAL_SOLUTION(NM,T_BEGIN,RK_STAGE=2)
   IF (PERIODIC_TEST==9) CALL SANDIA_DAT(NM,UVW_FILE)
   IF (PERIODIC_TEST==10) CALL TWOD_VORTEX_UMD(NM)
   IF (UVW_RESTART)      CALL UVW_INIT(NM,CSVFINFO(NM)%UVWFILE)
   CALL COMPUTE_VISCOSITY(T_BEGIN,NM) ! needed here for KRES prior to mesh exchange
ENDDO

! Exchange information at mesh boundaries related to the various initialization routines just completed

CALL MESH_EXCHANGE(1)
CALL MESH_EXCHANGE(4)
CALL POST_RECEIVES(6)
CALL MESH_EXCHANGE(6)

! Ensure normal components of velocity match at mesh boundaries and do velocity BCs just in case the flow is not initialized to zero

PREDICTOR = .FALSE.
CORRECTOR = .TRUE.

DO NM=1,NMESHES
   IF (PROCESS(NM)/=MYID .OR. TGA_SURF_INDEX>0) CYCLE
   CALL MATCH_VELOCITY(NM)
   CALL COMPUTE_VISCOSITY(T_BEGIN,NM) ! call again after mesh exchange
   IF (SYNTHETIC_EDDY_METHOD) CALL SYNTHETIC_EDDY_SETUP(NM)
   CALL VELOCITY_BC(T_BEGIN,NM)
   CALL VISCOSITY_BC(NM)
ENDDO

! Iterate surface BCs and radiation in case temperatures are not initialized to ambient

DO I=1,NUMBER_INITIAL_ITERATIONS
   DO NM=1,NMESHES
      IF (PROCESS(NM)/=MYID) CYCLE
      CALL WALL_BC(T_BEGIN,NM)
      IF (RADIATION) CALL COMPUTE_RADIATION(T_BEGIN,NM)
   ENDDO
   DO ANG_INC_COUNTER=1,ANGLE_INCREMENT
      CALL MESH_EXCHANGE(2) ! Exchange radiation intensity at interpolated boundaries
   ENDDO
ENDDO

! Compute divergence just in case the flow field is not initialized to ambient

DO NM=1,NMESHES
   IF (PROCESS(NM)/=MYID) CYCLE
   CALL DIVERGENCE_PART_1(T_BEGIN,NM)
ENDDO

! ! Restrict velocity components from fine mesh to coarse mesh

! IF (EMBEDDED_MESH) THEN
!   DO NM=NMESHES,1,-1
!      DO NOM=NMESHES,1,-1
!         IF (MESHES(NOM)%MESH_LEVEL<=MESHES(NM)%MESH_LEVEL) CYCLE
!         CALL VELOCITY_EMB(NM,NOM,IERROR) ! NM=coarse, NOM=fine
!      ENDDO
!   ENDDO
! ENDIF

! ! Apply normal boundary conditions to embedded meshes
   
! IF (EMBEDDED_MESH) THEN
!    DO NM=1,NMESHES
!       DO NOM=1,NMESHES
!          IF (MESHES(NOM)%MESH_LEVEL<=MESHES(NM)%MESH_LEVEL) CYCLE
!          CALL MATCH_VELOCITY_EMB(NM,NOM,IERROR,T_BEGIN)
!          IF (IERROR==0) CALL PROJECT_VELOCITY(NOM)
!       ENDDO
!    ENDDO
! ENDIF

! Potentially read data from a previous calculation 
 
DO NM=1,NMESHES
   IF (RESTART .AND. PROCESS(NM)==MYID) CALL READ_RESTART(T(NM),NM)
ENDDO
CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
 
! Initialize output files that are mesh-specific
 
DO NM=1,NMESHES
   IF (PROCESS(NM)/=MYID) CYCLE
   IF (TGA_SURF_INDEX<1) CALL INITIALIZE_MESH_DUMPS(NM)
   CALL INITIALIZE_PARTICLES(NM)
   CALL INSERT_PARTICLES(T(NM),NM)
   IF (TGA_SURF_INDEX<1) CALL INITIALIZE_DEVICES(NM)
   IF (TGA_SURF_INDEX<1) CALL INITIALIZE_PROFILES(NM)
ENDDO
CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

! Check for any stop flags at this point in the set up.

CALL STOP_CHECK(1)

! Check to see if only a TGA analysis is to be performed

IF (TGA_SURF_INDEX>0) THEN
   IF (MYID==0) CALL TGA_ANALYSIS
   STOP_STATUS = TGA_ANALYSIS_STOP
   CALL STOP_CHECK(1)
ENDIF

! Initialize output files containing global data (Master Node Only)

IF (MYID==0) CALL INITIALIZE_GLOBAL_DUMPS(T(1))
CALL INIT_EVAC_DUMPS
CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

! Initialize EVACuation routines

IF (ANY(EVACUATION_GRID)) THEN
   CALL INITIALIZE_EVAC
   IF ((.NOT.USE_MPI .OR. N_MPI_PROCESSES==1) .OR. (N_MPI_PROCESSES>1 .AND. MYID==EVAC_PROCESS)) CALL INIT_EVAC_GROUPS
ENDIF

! Initialize HVAC variables
IF (HVAC_SOLVE) THEN
   DO NM=1,NMESHES
      IF (PROCESS(NM)/=MYID) CYCLE  
      CALL HVAC_BC_IN(NM)
   ENDDO
   IF (N_MPI_PROCESSES>1) CALL EXCHANGE_HVAC_BC    
   IF (PROCESS(1)==MYID) THEN
      CALL COLLAPSE_HVAC_BC
      CALL SET_INIT_HVAC
   ENDIF
ENDIF

! Make an initial dump of ambient values

IF (.NOT.RESTART) THEN
   DO NM=1,NMESHES
      IF (PROCESS(NM)/=MYID) CYCLE
      CALL UPDATE_GLOBAL_OUTPUTS(T(NM),NM)      
      CALL DUMP_MESH_OUTPUTS(T(NM),NM)
   ENDDO
ENDIF

! Ensure the time is known to all meshes

IF (N_MPI_PROCESSES>1) THEN
   CALL MPI_ALLGATHERV(MPI_IN_PLACE,COUNTS(MYID),MPI_DOUBLE_PRECISION,T,COUNTS,DISPLS,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERR)
   IF (IERR/=MPI_SUCCESS) CALL HANDLE_MPI_ERROR('Error in MPI_ALLGATHERV for COUNTS',IERR)
ENDIF

! If there are zones and HVAC pass PSUM

IF (HVAC_SOLVE .AND. N_ZONE>0) CALL EXCHANGE_DIVERGENCE_INFO

! Make an initial dump of global output quantities

IF (.NOT.RESTART) THEN
   CALL UPDATE_CONTROLS(T,0._EB,CTRL_STOP_STATUS,.TRUE.)
   CALL DUMP_GLOBAL_OUTPUTS(T(1))
ENDIF

! Check for changes in VENT or OBSTruction control and device status at t=T_BEGIN

IF (.NOT.RESTART) THEN
   DO NM=1,NMESHES
      IF (PROCESS(NM)==MYID) CALL OPEN_AND_CLOSE(T(NM),NM)  
   ENDDO
ENDIF

! Write out character strings to .smv file
 
CALL WRITE_STRINGS

! Check for evacuation initialization stop
 
IF (ANY(EVACUATION_GRID)) CALL STOP_CHECK(1)

! Sprinkler piping calculation

DO CNT=1,N_DEVC
   IF (DEVICE(CNT)%PROP_INDEX > 0 .AND.  DEVICE(CNT)%CURRENT_STATE) THEN
      IF (PROPERTY(DEVICE(CNT)%PROP_INDEX)%PART_INDEX > 0) DEVC_PIPE_OPERATING(DEVICE(CNT)%PIPE_INDEX) = &
         DEVC_PIPE_OPERATING(DEVICE(CNT)%PIPE_INDEX) + 1
   ENDIF
ENDDO

! Start the clock for time stepping

WALL_CLOCK_START_ITERATIONS = WALL_CLOCK_TIME()

! Level Set model for firespread in vegetation (currently uses constant wind: does not need CFD computations).

IF (VEG_LEVEL_SET_UNCOUPLED .OR. VEG_LEVEL_SET_COUPLED) CALL INITIALIZE_LEVEL_SET_FIRESPREAD(1)
IF (VEG_LEVEL_SET_UNCOUPLED) THEN
  CALL LEVEL_SET_FIRESPREAD(T(1),1)
  STOP_STATUS = LEVELSET_STOP
  CALL STOP_CHECK(1)
ENDIF

! This ends the initialization part of the program

INITIALIZATION_PHASE = .FALSE.

!***********************************************************************************************************************************
!                                                   MAIN TIMESTEPPING LOOP
!***********************************************************************************************************************************
 
MAIN_LOOP: DO  
   
   ICYC  = ICYC + 1   ! Time step iterations

   ! Do not print out general diagnostics into .out file every time step

   DIAGNOSTICS = .FALSE.
   EXCHANGE_EVACUATION = .FALSE.

   ! Synchronize clocks
   
   IF (N_MPI_PROCESSES>1) THEN
      CALL MPI_ALLGATHERV(MPI_IN_PLACE,COUNTS(MYID),MPI_DOUBLE_PRECISION,T,COUNTS,DISPLS,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERR)
      IF (IERR/=MPI_SUCCESS) CALL HANDLE_MPI_ERROR('Error in MPI_ALLGATHERV for T',IERR)
   ENDIF
   
   ! Check for program stops
 
   INQUIRE(FILE=FN_STOP,EXIST=EX)
   IF (EX .AND. ICYC>=STOP_AT_ITER) THEN
      STOP_STATUS = USER_STOP
      DIAGNOSTICS = .TRUE.
   ENDIF

   ! Figure out fastest and slowest meshes
 
   T_MAX = MAXVAL(T,MASK=.NOT.EVACUATION_ONLY)
   T_MIN = MINVAL(T,MASK=.NOT.EVACUATION_ONLY)
   IF (ALL(EVACUATION_ONLY)) T_MAX = T_EVAC
   IF (ALL(EVACUATION_ONLY)) T_MIN = T_EVAC
 
   ! Determine time step
 
   DO NM=1,NMESHES
      IF (PROCESS(NM)==MYID) DT_SYNC(NM) = MESHES(NM)%DT_NEXT
   ENDDO

   IF (N_MPI_PROCESSES>1) THEN
      CALL MPI_ALLGATHERV(MPI_IN_PLACE,COUNTS(MYID),MPI_DOUBLE_PRECISION,DT_SYNC,COUNTS,DISPLS, &
                                    MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERR)
      IF (IERR/=MPI_SUCCESS) CALL HANDLE_MPI_ERROR('Error in MPI_ALLGATHERV for DT_SYNC',IERR)
   ENDIF

   DO NM=1,NMESHES
      IF (PROCESS(NM)/=MYID) CYCLE
      MESHES(NM)%DT_NEXT = MINVAL(DT_SYNC,MASK=.NOT.EVACUATION_ONLY)
      T(NM) = MINVAL(T,MASK=.NOT.EVACUATION_ONLY)
      ACTIVE_MESH(NM) = .TRUE.
   ENDDO

   ! Determine when to dump out diagnostics to the .out file

   LO10 = LOG10(REAL(MAX(1,ABS(ICYC)),EB))
   IF (MOD(ICYC,10**LO10)==0 .OR. MOD(ICYC,100)==0 .OR. T_MIN>=T_END) DIAGNOSTICS = .TRUE.

   ! Give every processor the full ACTIVE_MESH array
 
   IF (N_MPI_PROCESSES>1) THEN
      CALL MPI_ALLGATHERV(MPI_IN_PLACE,COUNTS(MYID),MPI_LOGICAL,ACTIVE_MESH,COUNTS,DISPLS,MPI_LOGICAL,MPI_COMM_WORLD,IERR)
      IF (IERR/=MPI_SUCCESS) CALL HANDLE_MPI_ERROR('Error in MPI_ALLGATHERV for ACTIVE_MESH',IERR)
   ENDIF

   ! If no meshes are due to be updated, update them all
 
   IF (ALL(.NOT.ACTIVE_MESH)) ACTIVE_MESH = .TRUE.

   ! If evacuation, set up special time iteration parameters

   CALL EVAC_MAIN_LOOP

   !============================================================================================================================
   !                                          Start of Predictor part of time step
   !============================================================================================================================
 
   PREDICTOR = .TRUE.
   CORRECTOR = .FALSE.

   ! Diagnostic timing calls and initialize energy budget array, Q_DOT

   DO NM=1,NMESHES
      IF (PROCESS(NM)/=MYID) CYCLE
      Q_DOT(:,NM) = 0._EB
   ENDDO

   ! Begin the finite differencing of the PREDICTOR step

   COMPUTE_FINITE_DIFFERENCES_1: DO NM=1,NMESHES
      IF (PROCESS(NM)/=MYID .OR. .NOT.ACTIVE_MESH(NM)) CYCLE COMPUTE_FINITE_DIFFERENCES_1
      MESHES(NM)%DT = MESHES(NM)%DT_NEXT
      NTCYC(NM)     = NTCYC(NM) + 1

      CALL INSERT_PARTICLES(T(NM),NM)
      CALL COMPUTE_VELOCITY_FLUX(T(NM),NM,1)
      CALL MASS_FINITE_DIFFERENCES(NM)
   ENDDO COMPUTE_FINITE_DIFFERENCES_1

!    ! Retrict scalar flux from fine mesh to coarse mesh

!    IF (EMBEDDED_MESH) THEN
!       DO NM=NMESHES,1,-1
!          DO NOM=NMESHES,1,-1
!             IF (MESHES(NOM)%MESH_LEVEL<=MESHES(NM)%MESH_LEVEL) CYCLE
!             CALL SCALARF_EMB(NM,NOM,IERROR)
!          ENDDO
!       ENDDO
!    ENDIF

   ! Estimate quantities at next time step, and decrease/increase time step if necessary based on CFL condition

   FIRST_PASS = .TRUE.
 
   CHANGE_TIME_STEP_LOOP: DO

      ! Predict density and mass fractions at next time step, and then start the divergence calculation
 
      COMPUTE_DENSITY_LOOP: DO NM=1,NMESHES
         IF (PROCESS(NM)/=MYID .OR. .NOT.ACTIVE_MESH(NM)) CYCLE COMPUTE_DENSITY_LOOP
         CALL DENSITY(T(NM),NM)
      ENDDO COMPUTE_DENSITY_LOOP

      ! Restrict mass from fine mesh to coarse mesh (redundant if SCALARF_EMB is done correctly)

!       IF (EMBEDDED_MESH) THEN
!          DO NM=NMESHES,1,-1
!             DO NOM=NMESHES,1,-1
!                IF (MESHES(NOM)%MESH_LEVEL<=MESHES(NM)%MESH_LEVEL) CYCLE
!                CALL RESTRICT_MASS_EMB(NM,NOM,IERROR)
!             ENDDO
!          ENDDO
!       ENDIF
      
      ! Exchange density and species mass fractions in interpolated boundaries

      IF (FIRST_PASS) CALL MESH_EXCHANGE(1)

      ! Do mass and energy boundary conditions, and begin divergence calculation

      COMPUTE_DIVERGENCE_LOOP: DO NM=1,NMESHES
         IF (PROCESS(NM)/=MYID .OR. .NOT.ACTIVE_MESH(NM)) CYCLE COMPUTE_DIVERGENCE_LOOP
         IF (N_FACE>0 .AND. FIRST_PASS) CALL INIT_IBM(T(NM),NM)
         CALL COMPUTE_VELOCITY_FLUX(T(NM),NM,2)
         IF (FIRST_PASS .AND. HVAC_SOLVE) CALL HVAC_BC_IN(NM)
      ENDDO COMPUTE_DIVERGENCE_LOOP

      IF (HVAC_SOLVE) THEN
         IF (FIRST_PASS .AND. N_MPI_PROCESSES>1) CALL EXCHANGE_HVAC_BC
         IF (PROCESS(1)==MYID) CALL HVAC_CALC(T(1),FIRST_PASS)
         IF (N_MPI_PROCESSES>1) CALL EXCHANGE_HVAC_SOLUTION
      ENDIF

      COMPUTE_WALL_BC_LOOP_A: DO NM=1,NMESHES
         IF (PROCESS(NM)/=MYID .OR. .NOT.ACTIVE_MESH(NM)) CYCLE COMPUTE_WALL_BC_LOOP_A
         CALL UPDATE_PARTICLES(T(NM),NM)
         CALL WALL_BC(T(NM),NM)
         CALL PARTICLE_MOMENTUM_TRANSFER(NM)
         CALL DIVERGENCE_PART_1(T(NM),NM)
      ENDDO COMPUTE_WALL_BC_LOOP_A

      ! Apply coarse mesh scalar bc to fine mesh

!       IF (EMBEDDED_MESH) THEN
!          DO NM=1,NMESHES
!             DO NOM=1,NMESHES
!                IF (MESHES(NOM)%MESH_LEVEL<=MESHES(NM)%MESH_LEVEL) CYCLE
!                CALL SCALAR_GHOST_EMB(NM,NOM,IERROR)
!             ENDDO
!          ENDDO
!          DO NM=1,NMESHES
!             IF (MESHES(NM)%MESH_LEVEL==0) CYCLE
!             CALL DIVERGENCE_PART_1(T(NM),NM)
!          ENDDO
!       ENDIF

      ! If there are pressure ZONEs, exchange integrated quantities mesh to mesh for use in the divergence calculation

      IF (N_ZONE>0) CALL EXCHANGE_DIVERGENCE_INFO

      ! Finish the divergence calculation

      FINISH_DIVERGENCE_LOOP: DO NM=1,NMESHES
         IF (PROCESS(NM)/=MYID .OR. .NOT.ACTIVE_MESH(NM)) CYCLE FINISH_DIVERGENCE_LOOP
         CALL DIVERGENCE_PART_2(NM)
      ENDDO FINISH_DIVERGENCE_LOOP

!       ! Restrict fine mesh divergence to coarse mesh

!       IF (EMBEDDED_MESH) THEN
!          DO NM=NMESHES,1,-1
!             DO NOM=NMESHES,1,-1
!                IF (MESHES(NOM)%MESH_LEVEL<=MESHES(NM)%MESH_LEVEL) CYCLE
!                CALL RESTRICT_DIV_EMB(NM,NOM,IERROR)
!             ENDDO
!          ENDDO
!       ENDIF

      ! Solve for the pressure at the current time step

      CALL PRESSURE_ITERATION_SCHEME
      CALL EVAC_PRESSURE_ITERATION_SCHEME

      ! Predict the velocity components at the next time step

      PREDICT_VELOCITY_LOOP: DO NM=1,NMESHES
         IF (PROCESS(NM)/=MYID .OR. .NOT.ACTIVE_MESH(NM)) CYCLE PREDICT_VELOCITY_LOOP
         CALL VELOCITY_PREDICTOR(T(NM)+MESHES(NM)%DT,NM)
      ENDDO PREDICT_VELOCITY_LOOP

      ! Check if there is a numerical instability after updating the velocity field. If there is, finish the time step.

      CALL STOP_CHECK(0)

      IF (STOP_STATUS==INSTABILITY_STOP) THEN
         DIAGNOSTICS = .TRUE.
         EXIT CHANGE_TIME_STEP_LOOP
      ENDIF

      ! Exchange information about the time step status, and if need be, repeat the CHANGE_TIME_STEP_LOOP

      IF (N_MPI_PROCESSES>1) THEN
         CALL MPI_ALLGATHERV(MPI_IN_PLACE,COUNTS(MYID),MPI_LOGICAL,CHANGE_TIME_STEP,COUNTS,DISPLS,MPI_LOGICAL,&
                             MPI_COMM_WORLD,IERR)
         IF (IERR/=MPI_SUCCESS) CALL HANDLE_MPI_ERROR('Error in MPI_ALLGATHERV for CHANGE_TIME_STEP',IERR)
      ENDIF
      IF (ANY(CHANGE_TIME_STEP)) THEN
         CHANGE_TIME_STEP = .TRUE.
         DO NM=1,NMESHES
            IF (EVACUATION_ONLY(NM)) CHANGE_TIME_STEP(NM) = .FALSE.
            IF (PROCESS(NM)/=MYID) CYCLE
            DT_SYNC(NM)      = MESHES(NM)%DT
            DT_NEXT_SYNC(NM) = MESHES(NM)%DT_NEXT
         ENDDO
         IF (N_MPI_PROCESSES>1) THEN
            CALL MPI_ALLGATHERV(MPI_IN_PLACE,COUNTS(MYID),MPI_DOUBLE_PRECISION,DT_SYNC,COUNTS,DISPLS, &
                                MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERR)
            IF (IERR/=MPI_SUCCESS) CALL HANDLE_MPI_ERROR('Error in MPI_ALLGATHERV for DT_SYNC',IERR)
            CALL MPI_ALLGATHERV(MPI_IN_PLACE,COUNTS(MYID),MPI_DOUBLE_PRECISION,DT_NEXT_SYNC,COUNTS,DISPLS, &
                                MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERR)
            IF (IERR/=MPI_SUCCESS) CALL HANDLE_MPI_ERROR('Error in MPI_ALLGATHERV for DT_NEXT_SYNC',IERR)
         ENDIF
         DO NM=1,NMESHES
            IF (PROCESS(NM)/=MYID) CYCLE
            IF (EVACUATION_ONLY(NM)) CYCLE
            MESHES(NM)%DT_NEXT = MINVAL(DT_NEXT_SYNC)
            MESHES(NM)%DT      = MINVAL(DT_SYNC)
         ENDDO
      ENDIF

      IF (.NOT.ANY(CHANGE_TIME_STEP)) EXIT CHANGE_TIME_STEP_LOOP
 
      FIRST_PASS = .FALSE.

   ENDDO CHANGE_TIME_STEP_LOOP

   CHANGE_TIME_STEP = .FALSE.
 
   ! Exchange velocity and pressures at interpolated boundaries

   CALL MESH_EXCHANGE(3)

   ! Force normal components of velocity to match at interpolated boundaries

   DO NM=1,NMESHES
      IF (PROCESS(NM)/=MYID .OR. .NOT.ACTIVE_MESH(NM)) CYCLE
      CALL MATCH_VELOCITY(NM)
   ENDDO

!    ! Restrict fine mesh velocity components to coarse mesh

!    IF (EMBEDDED_MESH) THEN
!       DO NM=NMESHES,1,-1
!          DO NOM=NMESHES,1,-1
!             IF (MESHES(NOM)%MESH_LEVEL<=MESHES(NM)%MESH_LEVEL) CYCLE
!             CALL VELOCITY_EMB(NM,NOM,IERROR)
!          ENDDO
!       ENDDO
!    ENDIF
   
!    ! Apply normal boundary conditions to embedded meshes
   
!    IF (EMBEDDED_MESH) THEN
!       DO NM=1,NMESHES
!          DO NOM=1,NMESHES
!             IF (MESHES(NOM)%MESH_LEVEL<=MESHES(NM)%MESH_LEVEL) CYCLE
!             CALL MATCH_VELOCITY_EMB(NM,NOM,IERROR,T(NM)+MESHES(NM)%DT)
!             IF (IERROR==0) CALL PROJECT_VELOCITY(NOM)
!          ENDDO
!       ENDDO
!    ENDIF

   ! Apply tangential velocity boundary conditions

   VELOCITY_BC_LOOP: DO NM=1,NMESHES
      IF (PROCESS(NM)/=MYID .OR. .NOT.ACTIVE_MESH(NM)) CYCLE VELOCITY_BC_LOOP
      IF (SYNTHETIC_EDDY_METHOD) CALL SYNTHETIC_TURBULENCE(MESHES(NM)%DT,T(NM),NM)
      CALL VELOCITY_BC(T(NM),NM)
   ENDDO VELOCITY_BC_LOOP

   ! Advance the time to start the CORRECTOR step
 
   UPDATE_TIME: DO NM=1,NMESHES
      IF (PROCESS(NM)/=MYID .OR. .NOT.ACTIVE_MESH(NM)) CYCLE UPDATE_TIME
      T(NM) = T(NM) + MESHES(NM)%DT  
   ENDDO UPDATE_TIME

   ! Ensure the time is known to all meshes

   IF (N_MPI_PROCESSES>1) THEN
      CALL MPI_ALLGATHERV(MPI_IN_PLACE,COUNTS(MYID),MPI_DOUBLE_PRECISION,T,COUNTS,DISPLS,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERR)
      IF (IERR/=MPI_SUCCESS) CALL HANDLE_MPI_ERROR('Error in MPI_ALLGATHERV for T',IERR)
   ENDIF
   T_MAX = MAXVAL(T,MASK=.NOT.EVACUATION_ONLY)
   T_MIN = MINVAL(T,MASK=.NOT.EVACUATION_ONLY)
   IF (ALL(EVACUATION_ONLY)) T_MAX = T_EVAC
   IF (ALL(EVACUATION_ONLY)) T_MIN = T_EVAC

   !===============================================================================================================================
   !                                          Start of Corrector part of time step
   !===============================================================================================================================
 
   CORRECTOR = .TRUE.
   PREDICTOR = .FALSE.

!    ! Retrict scalar flux from fine mesh to coarse mesh

!    IF (EMBEDDED_MESH) THEN
!       DO NM=NMESHES,1,-1
!          DO NOM=NMESHES,1,-1
!             IF (MESHES(NOM)%MESH_LEVEL<=MESHES(NM)%MESH_LEVEL) CYCLE
!             CALL SCALARF_EMB(NM,NOM,IERROR)
!          ENDDO
!       ENDDO
!    ENDIF

   ! Finite differences for mass and momentum equations for the second half of the time step

   COMPUTE_FINITE_DIFFERENCES_2: DO NM=1,NMESHES
      IF (PROCESS(NM)/=MYID)    CYCLE COMPUTE_FINITE_DIFFERENCES_2
      CALL OPEN_AND_CLOSE(T(NM),NM)   
      IF (.NOT.ACTIVE_MESH(NM)) CYCLE COMPUTE_FINITE_DIFFERENCES_2
      CALL COMPUTE_VELOCITY_FLUX(T(NM),NM,1)
      CALL MASS_FINITE_DIFFERENCES(NM)
      CALL DENSITY(T(NM),NM)
   ENDDO COMPUTE_FINITE_DIFFERENCES_2

!    ! Restrict mass from fine mesh to coarse mesh (redundant if SCALARF_EMB is done correctly)

!    IF (EMBEDDED_MESH) THEN
!       DO NM=NMESHES,1,-1
!          DO NOM=NMESHES,1,-1
!             IF (MESHES(NOM)%MESH_LEVEL/=MESHES(NM)%MESH_LEVEL+1) CYCLE
!             CALL RESTRICT_MASS_EMB(NM,NOM,IERROR)
!          ENDDO
!       ENDDO
!    ENDIF

   ! Exchange density and mass species

   CALL MESH_EXCHANGE(4)

   ! Apply mass and species boundary conditions, update radiation, particles, and re-compute divergence

   COMPUTE_DIVERGENCE_2: DO NM=1,NMESHES
      IF (PROCESS(NM)/=MYID .OR. .NOT.ACTIVE_MESH(NM)) CYCLE COMPUTE_DIVERGENCE_2
      CALL COMPUTE_VELOCITY_FLUX(T(NM),NM,2)
      IF (AGGLOMERATION) CALL CALC_AGGLOMERATION(NM)
      IF (N_REACTIONS > 0) CALL COMBUSTION (NM)
   ENDDO COMPUTE_DIVERGENCE_2

!    ! Restrict mass from fine mesh to coarse mesh (absolutely necessary to capture combustion source term)

!    IF (EMBEDDED_MESH .AND. N_REACTIONS>0) THEN
!       DO NM=NMESHES,1,-1
!          DO NOM=NMESHES,1,-1
!             IF (MESHES(NOM)%MESH_LEVEL<=MESHES(NM)%MESH_LEVEL) CYCLE
!             CALL RESTRICT_MASS_EMB(NM,NOM,IERROR)
!          ENDDO
!       ENDDO
!    ENDIF

   IF (HVAC_SOLVE .AND. ACTIVE_MESH(1)) CALL HVAC_CALC(T(1),.TRUE.)
 
   COMPUTE_WALL_BC_2A: DO NM=1,NMESHES
      IF (PROCESS(NM)/=MYID .OR. .NOT.ACTIVE_MESH(NM)) CYCLE COMPUTE_WALL_BC_2A
      CALL UPDATE_PARTICLES(T(NM),NM)
      CALL WALL_BC(T(NM),NM)
      CALL PARTICLE_MOMENTUM_TRANSFER(NM)
      CALL BNDRY_VEG_MASS_ENERGY_TRANSFER(T(NM),NM)
      IF (VEG_LEVEL_SET_COUPLED) CALL LEVEL_SET_FIRESPREAD(T(1),1)
      CALL COMPUTE_RADIATION(T(NM),NM)
      CALL DIVERGENCE_PART_1(T(NM),NM)
   ENDDO COMPUTE_WALL_BC_2A

   ! Apply coarse mesh scalar bc to fine mesh
   
!    IF (.FALSE.) THEN
!       DO NM=1,NMESHES
!          DO NOM=1,NMESHES
!             IF (MESHES(NOM)%MESH_LEVEL<=MESHES(NM)%MESH_LEVEL) CYCLE
!             CALL SCALAR_GHOST_EMB(NM,NOM,IERROR)
!          ENDDO
!       ENDDO
!       DO NM=1,NMESHES
!          IF (MESHES(NM)%MESH_LEVEL==0) CYCLE
!          CALL DIVERGENCE_PART_1(T(NM),NM)
!       ENDDO
!    ENDIF

   ! Exchange global pressure zone information

   IF (N_ZONE>0) CALL EXCHANGE_DIVERGENCE_INFO

   ! Finish computing the divergence

   FINISH_DIVERGENCE_LOOP_2: DO NM=1,NMESHES
      IF (PROCESS(NM)/=MYID .OR. .NOT.ACTIVE_MESH(NM)) CYCLE FINISH_DIVERGENCE_LOOP_2
      CALL DIVERGENCE_PART_2(NM)
   ENDDO FINISH_DIVERGENCE_LOOP_2

!    ! Restrict fine mesh divergence to coarse mesh

!    IF (EMBEDDED_MESH) THEN
!       DO NM=NMESHES,1,-1
!          DO NOM=NMESHES,1,-1
!             IF (MESHES(NOM)%MESH_LEVEL<=MESHES(NM)%MESH_LEVEL) CYCLE
!             CALL RESTRICT_DIV_EMB(NM,NOM,IERROR)
!          ENDDO
!       ENDDO
!    ENDIF

   ! Solve the pressure equation

   CALL PRESSURE_ITERATION_SCHEME
   CALL EVAC_PRESSURE_ITERATION_SCHEME

   ! Set up the last big exchange of info

   CALL EVAC_MESH_EXCHANGE(T_EVAC,T_EVAC_SAVE,I_EVAC,ICYC,EXCHANGE_EVACUATION,0)   

   ! Correct the velocity

   CORRECT_VELOCITY_LOOP: DO NM=1,NMESHES
      IF (PROCESS(NM)/=MYID)    CYCLE CORRECT_VELOCITY_LOOP
      IF (.NOT.ACTIVE_MESH(NM)) CYCLE CORRECT_VELOCITY_LOOP
      CALL VELOCITY_CORRECTOR(T(NM),NM)
      IF (DIAGNOSTICS) CALL CHECK_DIVERGENCE(NM)
   ENDDO CORRECT_VELOCITY_LOOP

   ! Exchange the number of particles sent from mesh to mesh

   CALL MESH_EXCHANGE(7)

   ! Exchange velocity, pressure, particles at interpolated boundaries

   CALL POST_RECEIVES(6) 
   CALL MESH_EXCHANGE(6)

   ! Exchange radiation at interpolated boundaries

   DO ANG_INC_COUNTER=1,ANGLE_INCREMENT
      CALL MESH_EXCHANGE(2)
      IF (ICYC>1) EXIT
   ENDDO

   ! Force normal components of velocity to match at interpolated boundaries

   DO NM=1,NMESHES
      IF (PROCESS(NM)/=MYID .OR. .NOT.ACTIVE_MESH(NM)) CYCLE
      CALL MATCH_VELOCITY(NM)
   ENDDO

!    ! Restrict fine mesh velocity components to coarse mesh
   
!    IF (EMBEDDED_MESH) THEN
!       DO NM=NMESHES,1,-1
!          DO NOM=NMESHES,1,-1
!             IF (MESHES(NOM)%MESH_LEVEL<=MESHES(NM)%MESH_LEVEL) CYCLE
!             CALL VELOCITY_EMB(NM,NOM,IERROR)
!          ENDDO
!       ENDDO
!    ENDIF

!    ! Apply normal boundary conditions to embedded meshes
   
!    IF (EMBEDDED_MESH) THEN
!       DO NM=1,NMESHES
!          DO NOM=1,NMESHES
!             IF (MESHES(NOM)%MESH_LEVEL<=MESHES(NM)%MESH_LEVEL) CYCLE
!             CALL MATCH_VELOCITY_EMB(NM,NOM,IERROR,T(NM))
!             IF (IERROR==0) CALL PROJECT_VELOCITY(NOM)
!          ENDDO
!       ENDDO
!    ENDIF

   ! Apply velocity boundary conditions, and update values of HRR, DEVC, etc.

   VELOCITY_BC_LOOP_2: DO NM=1,NMESHES
      IF (PROCESS(NM)/=MYID .OR. .NOT.ACTIVE_MESH(NM)) CYCLE VELOCITY_BC_LOOP_2
      CALL VELOCITY_BC(T(NM),NM)
      CALL UPDATE_GLOBAL_OUTPUTS(T(NM),NM)
   ENDDO VELOCITY_BC_LOOP_2

   ! Check for dumping end of timestep outputs

   CALL_UPDATE_CONTROLS = .FALSE.

   DO NM=1,NMESHES
      IF (PROCESS(NM)==MYID .AND. ACTIVE_MESH(NM)) THEN
         IF (.NOT. CALL_UPDATE_CONTROLS) THEN
            CALL UPDATE_CONTROLS(T,MESHES(NM)%DT,CTRL_STOP_STATUS,.FALSE.)
            IF (CTRL_STOP_STATUS) STOP_STATUS = CTRL_STOP
         ENDIF
         CALL_UPDATE_CONTROLS = .TRUE.
         CALL DUMP_MESH_OUTPUTS(T(NM),NM) 
      ENDIF
   ENDDO

   ! Dump outputs such as HRR, DEVC, etc.
   
   CALL DUMP_GLOBAL_OUTPUTS(T_MIN)

   ! Exchange EVAC information among meshes

   CALL EVAC_EXCHANGE

   ! Dump out diagnostics

   IF (DIAGNOSTICS .OR. T_MIN>=T_END) THEN
      CALL WRITE_STRINGS
      CALL EXCHANGE_DIAGNOSTICS
      IF (MYID==0) CALL WRITE_DIAGNOSTICS(T)
   ENDIF

   ! Flush output file buffers

   IF (T_MIN>=FLUSH_CLOCK .AND. FLUSH_FILE_BUFFERS) THEN
      IF (MYID==0) CALL FLUSH_GLOBAL_BUFFERS
      IF (MYID==MAX(0,EVAC_PROCESS)) CALL FLUSH_EVACUATION_BUFFERS
      DO NM=1,NMESHES
         IF (PROCESS(NM)==MYID) CALL FLUSH_LOCAL_BUFFERS(NM)
      ENDDO
      FLUSH_CLOCK = FLUSH_CLOCK + DT_FLUSH
   ENDIF

   ! Check for abnormal run stop

   CALL STOP_CHECK(1)  ! The argument 1 means that FDS will end unless there is logic associated with the STOP_STATUS

   ! Stop the run normally
   
   IF (T_MIN>=T_END) EXIT MAIN_LOOP

ENDDO MAIN_LOOP
 
!***********************************************************************************************************************************
!                                                     END OF TIME STEPPING LOOP
!***********************************************************************************************************************************
 
! Gather up timings for the meshes

DO NM=1,NMESHES
   IF (PROCESS(NM)==MYID) TUSED(1,NM) = SECOND() - TUSED(1,NM)
ENDDO

IF (N_MPI_PROCESSES>1) THEN
   REAL_BUFFER_3 = TUSED
   CALL MPI_GATHERV(REAL_BUFFER_3(1,DISPLS(MYID)+1),COUNTS_TIMERS(MYID),MPI_DOUBLE_PRECISION, &
                    TUSED,COUNTS_TIMERS,DISPLS_TIMERS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
   IF (IERR/=MPI_SUCCESS) CALL HANDLE_MPI_ERROR('Error in MPI_GATHERV for TUSED',IERR)
ENDIF

IF (MYID==0) CALL TIMINGS

IF (PRES_METHOD == 'SCARC') CALL SCARC_TIMINGS

! Stop the calculation

! Finish unstructured geometry

IF (N_FACE>0) THEN
   DO NM=1,NMESHES
      IF (PROCESS(NM)/=MYID) CYCLE
      CALL INIT_IBM(T(NM),NM)
   ENDDO
ENDIF

CALL END_FDS

! This is the end of program. Supporting routines are listed below.


CONTAINS


SUBROUTINE CHECK_MPI_OPENMP

INTEGER :: THREAD_ID

IF (.NOT.USE_OPENMP .AND. .NOT.USE_MPI) RETURN

! Check the threading support level

IF (USE_MPI .AND. PROVIDED<REQUIRED) THEN
   IF (MYID==0) WRITE(LU_ERR,*) "WARNING:  This MPI implementation provides insufficient threading support."
   !$ CALL OMP_SET_NUM_THREADS(1)
ENDIF

! The multi-threaded section where all threads will say hello

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(THREAD_ID)

THREAD_ID = 0
!$ THREAD_ID = OMP_GET_THREAD_NUM()  

!$OMP CRITICAL
IF (USE_OPENMP .AND. USE_MPI) WRITE(LU_ERR,91) " OpenMP thread ",THREAD_ID," of ",OPENMP_AVAILABLE_THREADS-1,&
   " assigned to MPI process ",MYID," of ",N_MPI_PROCESSES-1," is running on ",PNAME(1:PNAMELEN)
IF (.NOT.USE_OPENMP .AND. USE_MPI) WRITE(LU_ERR,92) " MPI process ",MYID," of ",N_MPI_PROCESSES-1,&
                                                    " is running on ",PNAME(1:PNAMELEN)
IF (USE_OPENMP .AND. .NOT.USE_MPI) WRITE(LU_ERR,93) " OpenMP thread ",THREAD_ID," of ",OPENMP_AVAILABLE_THREADS-1,&
      " is running"
!$OMP END CRITICAL

!$OMP END PARALLEL

91 FORMAT(A,I3,A,I3,A,I3,A,I3,A,A)
92 FORMAT(A,I3,A,I3,A,A)
93 FORMAT(A,I3,A,I3,A)

END SUBROUTINE CHECK_MPI_OPENMP


SUBROUTINE MPI_INITIALIZATION_CHORES(TASK_NUMBER)

INTEGER, INTENT(IN) :: TASK_NUMBER
INTEGER, ALLOCATABLE, DIMENSION(:) :: REQ0
INTEGER :: N_REQ0

SELECT CASE(TASK_NUMBER)

   CASE(1)

      ALLOCATE(CONNECTED_MESHES(NMESHES,NMESHES)) ; CONNECTED_MESHES = .FALSE.

      ! Set up send and receive buffer counts and displacements

      ALLOCATE(REAL_BUFFER_1(NMESHES))
      ALLOCATE(REAL_BUFFER_3(N_TIMERS_DIM,NMESHES))
      ALLOCATE(REAL_BUFFER_5(0:MINT_SPECIES,NMESHES))
      ALLOCATE(REAL_BUFFER_6(N_DUCTNODES,NMESHES))
      ALLOCATE(REAL_BUFFER_7(N_DUCTNODES,N_TOTAL_SCALARS,NMESHES))
      ALLOCATE(REAL_BUFFER_8(0:N_ZONE,0:N_ZONE,NMESHES))
      ALLOCATE(REAL_BUFFER_9(0:N_ZONE,0:N_ZONE,N_TOTAL_SCALARS,NMESHES))
      ALLOCATE(REAL_BUFFER_10(N_DUCTS,NMESHES))
      ALLOCATE(REAL_BUFFER_11(N_Q_DOT,NMESHES))
      ALLOCATE(REAL_BUFFER_12(N_M_DOT,NMESHES))
      ALLOCATE(INTEGER_BUFFER_1(NMESHES))
      ALLOCATE(INTEGER_BUFFER_2(N_DUCTNODES,NMESHES))      
      ALLOCATE(LOGICAL_BUFFER_1(NMESHES))
      
      ALLOCATE(COUNTS(0:N_MPI_PROCESSES-1))
      ALLOCATE(COUNTS2D(0:N_MPI_PROCESSES-1))
      ALLOCATE(COUNTS_TIMERS(0:N_MPI_PROCESSES-1))
      ALLOCATE(COUNTS_HVAC(0:N_MPI_PROCESSES-1))
      ALLOCATE(COUNTS_HVAC_SPECIES(0:N_MPI_PROCESSES-1))
      ALLOCATE(COUNTS_MASS(0:N_MPI_PROCESSES-1))
      ALLOCATE(COUNTS_Q_DOT(0:N_MPI_PROCESSES-1))
      ALLOCATE(COUNTS_M_DOT(0:N_MPI_PROCESSES-1))
      
      ALLOCATE(DISPLS(0:N_MPI_PROCESSES-1))
      ALLOCATE(DISPLS2D(0:N_MPI_PROCESSES-1))
      ALLOCATE(DISPLS_MASS(0:N_MPI_PROCESSES-1))
      ALLOCATE(DISPLS_TIMERS(0:N_MPI_PROCESSES-1))
      ALLOCATE(DISPLS_HVAC(0:N_MPI_PROCESSES-1))
      ALLOCATE(DISPLS_HVAC_SPECIES(0:N_MPI_PROCESSES-1))
      ALLOCATE(DISPLS_Q_DOT(0:N_MPI_PROCESSES-1))
      ALLOCATE(DISPLS_M_DOT(0:N_MPI_PROCESSES-1))
      
      COUNTS    = 0
      DO N=0,N_MPI_PROCESSES-1
         DO NM=1,NMESHES
            IF (PROCESS(NM)==N) COUNTS(N)    = COUNTS(N)    + 1
         ENDDO
      ENDDO
      DISPLS(0)    = 0
      DO N=1,N_MPI_PROCESSES-1
         DISPLS(N)    = COUNTS(N-1)    + DISPLS(N-1)
      ENDDO
      COUNTS2D      = COUNTS*NMESHES
      DISPLS2D      = DISPLS*NMESHES
      COUNTS_TIMERS = COUNTS*N_TIMERS_DIM
      DISPLS_TIMERS = DISPLS*N_TIMERS_DIM
      COUNTS_HVAC   = COUNTS*N_DUCTNODES
      COUNTS_HVAC_SPECIES   = COUNTS*N_DUCTNODES*N_TOTAL_SCALARS
      DISPLS_HVAC   = DISPLS*N_DUCTNODES
      DISPLS_HVAC_SPECIES   = DISPLS*N_DUCTNODES*N_TOTAL_SCALARS
      COUNTS_MASS   = COUNTS*(MINT_SPECIES+1)
      DISPLS_MASS   = DISPLS*(MINT_SPECIES+1)
      COUNTS_Q_DOT  = COUNTS*N_Q_DOT
      DISPLS_Q_DOT  = DISPLS*N_Q_DOT
      COUNTS_M_DOT  = COUNTS*N_M_DOT
      DISPLS_M_DOT  = DISPLS*N_M_DOT

   CASE(2)

      ! Allocate TIME arrays

      ALLOCATE(ACTIVE_MESH(NMESHES),STAT=IZERO)
      CALL ChkMemErr('MAIN','ACTIVE_MESH',IZERO)
      ACTIVE_MESH = .TRUE.
      ALLOCATE(T(NMESHES),STAT=IZERO)
      CALL ChkMemErr('MAIN','T',IZERO)
      ALLOCATE(DT_SYNC(NMESHES),STAT=IZERO)
      CALL ChkMemErr('MAIN','DT_SYNC',IZERO)
      ALLOCATE(DT_NEXT_SYNC(NMESHES),STAT=IZERO)
      CALL ChkMemErr('MAIN','DT_NEXT_SYNC',IZERO)
      
      ! Set up dummy arrays to hold various arrays that must be exchanged among meshes
      
      ALLOCATE(TI_LOC(N_DEVC),STAT=IZERO)
      CALL ChkMemErr('MAIN','TI_LOC',IZERO)
      ALLOCATE(TI_GLB(N_DEVC),STAT=IZERO)
      CALL ChkMemErr('MAIN','TI_GLB',IZERO)
      ALLOCATE(STATE_GLB(N_DEVC),STAT=IZERO)
      CALL ChkMemErr('MAIN','STATE_GLB',IZERO)
      ALLOCATE(STATE_LOC(N_DEVC),STAT=IZERO)
      CALL ChkMemErr('MAIN','STATE_LOC',IZERO)
      ALLOCATE(TC_GLB(N_DEVC),STAT=IZERO)
      CALL ChkMemErr('MAIN','TC_GLB',IZERO)
      ALLOCATE(TC_LOC(N_DEVC),STAT=IZERO)
      CALL ChkMemErr('MAIN','TC_LOC',IZERO)
      
      ! Allocate a few arrays needed to exchange divergence and pressure info among meshes
      
      IF (N_ZONE > 0) THEN
         ALLOCATE(DSUM_ALL(N_ZONE),STAT=IZERO)
         ALLOCATE(PSUM_ALL(N_ZONE),STAT=IZERO)
         ALLOCATE(USUM_ALL(N_ZONE),STAT=IZERO)
         ALLOCATE(CONNECTED_ZONES_GLOBAL(0:N_ZONE,0:N_ZONE),STAT=IZERO)
         ALLOCATE(DSUM_ALL_LOCAL(N_ZONE),STAT=IZERO)
         ALLOCATE(PSUM_ALL_LOCAL(N_ZONE),STAT=IZERO)
         ALLOCATE(USUM_ALL_LOCAL(N_ZONE),STAT=IZERO)
         ALLOCATE(CONNECTED_ZONES_LOCAL(0:N_ZONE,0:N_ZONE),STAT=IZERO)
      ENDIF

   CASE(3)

      ! Allocate "request" arrays to keep track of MPI communications

      ALLOCATE(REQ(N_COMMUNICATIONS*40))
      ALLOCATE(REQ1(N_COMMUNICATIONS*4))
      ALLOCATE(REQ2(N_COMMUNICATIONS*4))
      ALLOCATE(REQ3(N_COMMUNICATIONS*4))
      ALLOCATE(REQ4(N_COMMUNICATIONS*4))
      ALLOCATE(REQ5(N_COMMUNICATIONS*4))
      ALLOCATE(REQ6(N_COMMUNICATIONS*4))
      ALLOCATE(REQ7(N_COMMUNICATIONS*4))
      ALLOCATE(REQ8(N_COMMUNICATIONS*4))
      ALLOCATE(REQ9(N_COMMUNICATIONS*4))

      REQ = MPI_REQUEST_NULL
      REQ1 = MPI_REQUEST_NULL
      REQ2 = MPI_REQUEST_NULL
      REQ3 = MPI_REQUEST_NULL
      REQ4 = MPI_REQUEST_NULL
      REQ5 = MPI_REQUEST_NULL
      REQ6 = MPI_REQUEST_NULL
      REQ7 = MPI_REQUEST_NULL
      REQ8 = MPI_REQUEST_NULL
      REQ9 = MPI_REQUEST_NULL
      
   CASE(4)

      ! Exchange information related to size of OMESH arrays

      CALL MPI_ALLREDUCE(MPI_IN_PLACE,CONNECTED_MESHES(1,1),NMESHES**2,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,IERR)
      IF (IERR/=MPI_SUCCESS) CALL HANDLE_MPI_ERROR('Error in MPI_ALLREDUCE for CONNECTED_MESHES',IERR)

      IF (N_MPI_PROCESSES>1) ALLOCATE(REQ0(NMESHES**2))

      N_REQ0 = 0

      DO NM=1,NMESHES
         IF (EVACUATION_ONLY(NM)) CYCLE
         DO NOM=1,NMESHES
            IF (PROCESS(NOM)/=MYID) CYCLE
            IF (EVACUATION_ONLY(NOM)) CYCLE
            IF (N_MPI_PROCESSES>1 .AND. NM/=NOM .AND. PROCESS(NM)/=MYID .AND. CONNECTED_MESHES(NM,NOM)) THEN
               M2 => MESHES(NOM)%OMESH(NM)
               N_REQ0 = N_REQ0 + 1
               CALL MPI_IRECV(M2%INTEGER_RECV_BUFFER(1),7,MPI_INTEGER,PROCESS(NM),NM,MPI_COMM_WORLD,REQ0(N_REQ0),IERR)
            ENDIF
         ENDDO
      ENDDO

      DO NM=1,NMESHES
         IF (PROCESS(NM)/=MYID) CYCLE
         IF (EVACUATION_ONLY(NM)) CYCLE
         DO NOM=1,NMESHES
            IF (EVACUATION_ONLY(NOM)) CYCLE
            IF (.NOT.CONNECTED_MESHES(NM,NOM)) CYCLE
            M3 => MESHES(NM)%OMESH(NOM)
            IF (N_MPI_PROCESSES>1 .AND. NM/=NOM .AND. PROCESS(NOM)/=MYID .AND. CONNECTED_MESHES(NM,NOM)) THEN
               M3%INTEGER_SEND_BUFFER(1) = M3%I_MIN_R
               M3%INTEGER_SEND_BUFFER(2) = M3%I_MAX_R
               M3%INTEGER_SEND_BUFFER(3) = M3%J_MIN_R
               M3%INTEGER_SEND_BUFFER(4) = M3%J_MAX_R
               M3%INTEGER_SEND_BUFFER(5) = M3%K_MIN_R
               M3%INTEGER_SEND_BUFFER(6) = M3%K_MAX_R
               M3%INTEGER_SEND_BUFFER(7) = M3%NIC_R
               N_REQ0 = N_REQ0 + 1
               CALL MPI_ISEND(M3%INTEGER_SEND_BUFFER(1),7,MPI_INTEGER,PROCESS(NOM),NM,MPI_COMM_WORLD,REQ0(N_REQ0),IERR)
            ELSE
               M2 => MESHES(NOM)%OMESH(NM)
               M2%I_MIN_S = M3%I_MIN_R
               M2%I_MAX_S = M3%I_MAX_R
               M2%J_MIN_S = M3%J_MIN_R
               M2%J_MAX_S = M3%J_MAX_R
               M2%K_MIN_S = M3%K_MIN_R
               M2%K_MAX_S = M3%K_MAX_R
               M2%NIC_S   = M3%NIC_R
            ENDIF
         ENDDO
      ENDDO
 
      IF (N_MPI_PROCESSES>1) CALL MPI_WAITALL(N_REQ0,REQ0(1:N_REQ0),MPI_STATUSES_IGNORE,IERR)

      DO NM=1,NMESHES
         IF (EVACUATION_ONLY(NM)) CYCLE
         DO NOM=1,NMESHES
            IF (PROCESS(NOM)/=MYID) CYCLE
            IF (EVACUATION_ONLY(NOM)) CYCLE
            IF (N_MPI_PROCESSES>1 .AND. NM/=NOM .AND. PROCESS(NM)/=MYID .AND. CONNECTED_MESHES(NM,NOM)) THEN
               M2 => MESHES(NOM)%OMESH(NM)
               M2%I_MIN_S = M2%INTEGER_RECV_BUFFER(1)
               M2%I_MAX_S = M2%INTEGER_RECV_BUFFER(2)
               M2%J_MIN_S = M2%INTEGER_RECV_BUFFER(3)
               M2%J_MAX_S = M2%INTEGER_RECV_BUFFER(4)
               M2%K_MIN_S = M2%INTEGER_RECV_BUFFER(5)
               M2%K_MAX_S = M2%INTEGER_RECV_BUFFER(6)
               M2%NIC_S   = M2%INTEGER_RECV_BUFFER(7)
            ENDIF
         ENDDO
      ENDDO

      IF (N_MPI_PROCESSES>1) THEN 
         DEALLOCATE(REQ0)
         CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
      ENDIF

      ! Exchange IIO_S, etc., the indices of interpolated cells

      DO NM=1,NMESHES
         IF (EVACUATION_ONLY(NM)) CYCLE
         DO NOM=1,NMESHES
            IF (PROCESS(NOM)/=MYID) CYCLE
            IF (EVACUATION_ONLY(NOM)) CYCLE
            IF (MESHES(NOM)%OMESH(NM)%NIC_S>0) THEN
               ALLOCATE(MESHES(NOM)%OMESH(NM)%IIO_S(MESHES(NOM)%OMESH(NM)%NIC_S))
               ALLOCATE(MESHES(NOM)%OMESH(NM)%JJO_S(MESHES(NOM)%OMESH(NM)%NIC_S))
               ALLOCATE(MESHES(NOM)%OMESH(NM)%KKO_S(MESHES(NOM)%OMESH(NM)%NIC_S))
               ALLOCATE(MESHES(NOM)%OMESH(NM)%IOR_S(MESHES(NOM)%OMESH(NM)%NIC_S))
            ENDIF
         ENDDO
      ENDDO

      N_REQ = 0

      DO NM=1,NMESHES
         IF (EVACUATION_ONLY(NM)) CYCLE
         DO NOM=1,NMESHES
            IF (PROCESS(NOM)/=MYID) CYCLE
            IF (EVACUATION_ONLY(NOM)) CYCLE
            M2 => MESHES(NOM)%OMESH(NM)
            IF (N_MPI_PROCESSES>1 .AND. NM/=NOM .AND. PROCESS(NM)/=MYID .AND. M2%NIC_S>0) THEN
               CALL MPI_IRECV(M2%IIO_S(1),M2%NIC_S,MPI_INTEGER,PROCESS(NM),NM,MPI_COMM_WORLD,REQ(N_REQ+1),IERR)
               CALL MPI_IRECV(M2%JJO_S(1),M2%NIC_S,MPI_INTEGER,PROCESS(NM),NM,MPI_COMM_WORLD,REQ(N_REQ+2),IERR)
               CALL MPI_IRECV(M2%KKO_S(1),M2%NIC_S,MPI_INTEGER,PROCESS(NM),NM,MPI_COMM_WORLD,REQ(N_REQ+3),IERR)
               CALL MPI_IRECV(M2%IOR_S(1),M2%NIC_S,MPI_INTEGER,PROCESS(NM),NM,MPI_COMM_WORLD,REQ(N_REQ+4),IERR)
               N_REQ = N_REQ + 4
            ENDIF
         ENDDO
      ENDDO

      DO NM=1,NMESHES
         IF (PROCESS(NM)/=MYID) CYCLE
         IF (EVACUATION_ONLY(NM)) CYCLE
         DO NOM=1,NMESHES
            IF (EVACUATION_ONLY(NOM)) CYCLE
            M3 => MESHES(NM)%OMESH(NOM)
            IF (M3%NIC_R<1) CYCLE
            IF (PROCESS(NOM)/=MYID) THEN
               CALL MPI_ISEND(M3%IIO_R(1),M3%NIC_R,MPI_INTEGER,PROCESS(NOM),NM,MPI_COMM_WORLD,REQ(N_REQ+1),IERR)
               CALL MPI_ISEND(M3%JJO_R(1),M3%NIC_R,MPI_INTEGER,PROCESS(NOM),NM,MPI_COMM_WORLD,REQ(N_REQ+2),IERR)
               CALL MPI_ISEND(M3%KKO_R(1),M3%NIC_R,MPI_INTEGER,PROCESS(NOM),NM,MPI_COMM_WORLD,REQ(N_REQ+3),IERR)
               CALL MPI_ISEND(M3%IOR_R(1),M3%NIC_R,MPI_INTEGER,PROCESS(NOM),NM,MPI_COMM_WORLD,REQ(N_REQ+4),IERR)
               N_REQ = N_REQ + 4
            ELSE
               MESHES(NOM)%OMESH(NM)%IIO_S = M3%IIO_R
               MESHES(NOM)%OMESH(NM)%JJO_S = M3%JJO_R
               MESHES(NOM)%OMESH(NM)%KKO_S = M3%KKO_R
               MESHES(NOM)%OMESH(NM)%IOR_S = M3%IOR_R
            ENDIF
         ENDDO
      ENDDO

      IF (N_REQ>0 .AND. N_MPI_PROCESSES>1) CALL MPI_WAITALL(N_REQ,REQ(1:N_REQ),MPI_STATUSES_IGNORE,IERR)

END SELECT

CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

END SUBROUTINE MPI_INITIALIZATION_CHORES


SUBROUTINE PRESSURE_ITERATION_SCHEME
INTEGER :: NM_MAX

! Iterate calls to pressure solver until velocity tolerance is satisfied

CALL MESH_EXCHANGE(5)

DO NM=1,NMESHES
   IF (EVACUATION_ONLY(NM)) CYCLE
   IF (PROCESS(NM)==MYID .AND. ACTIVE_MESH(NM)) MESHES(NM)%WALL_WORK1 = 0._EB
   IF (PROCESS(NM)==MYID .AND. ACTIVE_MESH(NM)) CALL MATCH_VELOCITY_FLUX(NM)
ENDDO

PRESSURE_ITERATIONS = 0

PRESSURE_ITERATION_LOOP: DO

   PRESSURE_ITERATIONS = PRESSURE_ITERATIONS + 1
   TOTAL_PRESSURE_ITERATIONS = TOTAL_PRESSURE_ITERATIONS + 1

   DO NM=1,NMESHES
      IF (EVACUATION_ONLY(NM)) CYCLE
      IF (PROCESS(NM)==MYID .AND. ACTIVE_MESH(NM)) THEN
         CALL NO_FLUX(NM)
         CALL PRESSURE_SOLVER(T(NM),NM)
      ENDIF
   ENDDO
   IF (PRES_METHOD == 'SCARC') CALL SCARC_SOLVER

   IF (.NOT.ITERATE_PRESSURE) EXIT PRESSURE_ITERATION_LOOP

   CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
   CALL MESH_EXCHANGE(5)

   DO NM=1,NMESHES
      IF (EVACUATION_ONLY(NM)) CYCLE
      IF (PROCESS(NM)==MYID .AND. ACTIVE_MESH(NM)) CALL COMPUTE_VELOCITY_ERROR(NM)
   ENDDO

   IF (N_MPI_PROCESSES>1) THEN
      CALL MPI_ALLGATHERV(MPI_IN_PLACE,COUNTS(MYID),MPI_DOUBLE_PRECISION, &
                          VELOCITY_ERROR_MAX,COUNTS,DISPLS,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERR)
      IF (IERR/=MPI_SUCCESS) CALL HANDLE_MPI_ERROR('Error in MPI_ALLGATHERV for VELOCITY_ERROR_MAX',IERR)
      CALL MPI_ALLGATHERV(MPI_IN_PLACE,COUNTS(MYID),MPI_INTEGER, &
                          VELOCITY_ERROR_MAX_I,COUNTS,DISPLS,MPI_INTEGER,MPI_COMM_WORLD,IERR)
      IF (IERR/=MPI_SUCCESS) CALL HANDLE_MPI_ERROR('Error in MPI_ALLGATHERV for VELOCITY_ERROR_MAX_I',IERR)
      CALL MPI_ALLGATHERV(MPI_IN_PLACE,COUNTS(MYID),MPI_INTEGER, &
                          VELOCITY_ERROR_MAX_J,COUNTS,DISPLS,MPI_INTEGER,MPI_COMM_WORLD,IERR)
      IF (IERR/=MPI_SUCCESS) CALL HANDLE_MPI_ERROR('Error in MPI_ALLGATHERV for VELOCITY_ERROR_MAX_J',IERR)
      CALL MPI_ALLGATHERV(MPI_IN_PLACE,COUNTS(MYID),MPI_INTEGER, &
                          VELOCITY_ERROR_MAX_K,COUNTS,DISPLS,MPI_INTEGER,MPI_COMM_WORLD,IERR)
      IF (IERR/=MPI_SUCCESS) CALL HANDLE_MPI_ERROR('Error in MPI_ALLGATHERV for VELOCITY_ERROR_MAX_K',IERR)
   ENDIF
   IF (VELOCITY_ERROR_FILE .AND. .NOT.ALL(EVACUATION_ONLY)) THEN
      NM_MAX = MAXLOC(VELOCITY_ERROR_MAX,DIM=1)
      WRITE(LU_VELOCITY_ERROR,'(7(I7,A),E16.8)') ICYC,', ',PRESSURE_ITERATIONS,', ',NM_MAX,', ', &
         VELOCITY_ERROR_MAX_I(NM_MAX),', ',VELOCITY_ERROR_MAX_J(NM_MAX),', ',VELOCITY_ERROR_MAX_K(NM_MAX),', ', &
         TOTAL_PRESSURE_ITERATIONS,', ',MAXVAL(VELOCITY_ERROR_MAX)
   ENDIF

   IF (MAXVAL(VELOCITY_ERROR_MAX)<VELOCITY_TOLERANCE .OR. PRESSURE_ITERATIONS>=MAX_PRESSURE_ITERATIONS) &
      EXIT PRESSURE_ITERATION_LOOP

ENDDO PRESSURE_ITERATION_LOOP

END SUBROUTINE PRESSURE_ITERATION_SCHEME


SUBROUTINE STOP_CHECK(END_CODE)

INTEGER, INTENT(IN) :: END_CODE

! Make sure that all MPI processes have the same STOP_STATUS

IF (N_MPI_PROCESSES>1) THEN
   CALL MPI_ALLREDUCE(MPI_IN_PLACE,STOP_STATUS,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERR)
   IF (IERR/=MPI_SUCCESS) CALL HANDLE_MPI_ERROR('Error in MPI_ALLREDUCE for STOP_STATUS',IERR)
   CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
ENDIF

SELECT CASE(STOP_STATUS)
   CASE(NO_STOP) 
      RETURN
   CASE(USER_STOP)
      IF (STOP_AT_ITER==0.AND..NOT.RADIATION_COMPLETED) RETURN
END SELECT

IF (END_CODE==1) CALL END_FDS

END SUBROUTINE STOP_CHECK


SUBROUTINE END_FDS

! End the calculation gracefully, even if there is an error

CHARACTER(255) :: MESSAGE
LOGICAL :: OPN

IF (USE_MPI) THEN
   WRITE(LU_ERR,'(A,I4,A)') 'MPI process ',MYID,' has completed'
   CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
    IF (IERR/=MPI_SUCCESS) CALL HANDLE_MPI_ERROR('Error in last MPI_BARRIER call',IERR)
   CALL MPI_FINALIZE(IERR)
   IF (IERR/=MPI_SUCCESS) CALL HANDLE_MPI_ERROR('Error in MPI_FINALIZE',IERR)
ENDIF

IF (MYID==0) THEN

   SELECT CASE(STOP_STATUS)
      CASE(NO_STOP)
         WRITE(MESSAGE,'(A)') 'STOP: FDS completed successfully'
         IF (STATUS_FILES) CLOSE(LU_NOTREADY,STATUS='DELETE')
      CASE(INSTABILITY_STOP) 
         WRITE(MESSAGE,'(A)') 'STOP: Numerical Instability'
      CASE(USER_STOP) 
         WRITE(MESSAGE,'(A)') 'STOP: FDS stopped by user'
      CASE(SETUP_STOP) 
         WRITE(MESSAGE,'(A)') 'STOP: FDS was improperly set-up'
      CASE(SETUP_ONLY_STOP) 
         WRITE(MESSAGE,'(A)') 'STOP: Set-up only'
      CASE(CTRL_STOP) 
         WRITE(MESSAGE,'(A)') 'STOP: FDS was stopped by KILL control function'
      CASE(TGA_ANALYSIS_STOP) 
         WRITE(MESSAGE,'(A)') 'STOP: TGA analysis only'
      CASE(LEVELSET_STOP) 
         WRITE(MESSAGE,'(A)') 'STOP: Level set analysis only'
      CASE(REALIZABILITY_STOP) 
         WRITE(MESSAGE,'(A)') 'STOP: Unrealizable mass density'
   END SELECT

   WRITE(LU_ERR,'(/A,A,A,A)') TRIM(MESSAGE),' (CHID: ',TRIM(CHID),')'
   INQUIRE(LU_OUTPUT,OPENED=OPN)
   IF (OPN) WRITE(LU_OUTPUT,'(/A,A,A,A)') TRIM(MESSAGE),' (CHID: ',TRIM(CHID),')'

ENDIF

STOP

END SUBROUTINE END_FDS
 
 
SUBROUTINE EXCHANGE_DIVERGENCE_INFO

! Exchange information mesh to mesh needed for divergence integrals
! First, sum DSUM, PSUM and USUM over all meshes controlled by the active process, then reduce over all processes

INTEGER :: IPZ,IOPZ,IOPZ2
REAL(EB) :: TNOW

TNOW = SECOND()

CONNECTED_ZONES_LOCAL = .FALSE.

DO IPZ=1,N_ZONE
   DSUM_ALL_LOCAL(IPZ) = 0._EB
   PSUM_ALL_LOCAL(IPZ) = 0._EB
   USUM_ALL_LOCAL(IPZ) = 0._EB
   IF(P_ZONE(IPZ)%EVACUATION) CYCLE
   DO NM=1,NMESHES
      IF (PROCESS(NM)/=MYID) CYCLE
      IF(EVACUATION_ONLY(NM)) CYCLE
      DSUM_ALL_LOCAL(IPZ) = DSUM_ALL_LOCAL(IPZ) + DSUM(IPZ,NM)
      PSUM_ALL_LOCAL(IPZ) = PSUM_ALL_LOCAL(IPZ) + PSUM(IPZ,NM)
      USUM_ALL_LOCAL(IPZ) = USUM_ALL_LOCAL(IPZ) + USUM(IPZ,NM)
      DO IOPZ=0,N_ZONE
         IF (CONNECTED_ZONES(IPZ,IOPZ,NM)) CONNECTED_ZONES_LOCAL(IPZ,IOPZ) = .TRUE.
      ENDDO
   ENDDO
ENDDO

IF (N_MPI_PROCESSES>1) THEN
   CALL MPI_ALLREDUCE(DSUM_ALL_LOCAL(1),DSUM_ALL(1),N_ZONE,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
   IF (IERR/=MPI_SUCCESS) CALL HANDLE_MPI_ERROR('Error in MPI_ALLREDUCE for DSUM_ALL',IERR)
   CALL MPI_ALLREDUCE(PSUM_ALL_LOCAL(1),PSUM_ALL(1),N_ZONE,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
   IF (IERR/=MPI_SUCCESS) CALL HANDLE_MPI_ERROR('Error in MPI_ALLREDUCE for PSUM_ALL',IERR)
   CALL MPI_ALLREDUCE(USUM_ALL_LOCAL(1),USUM_ALL(1),N_ZONE,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
   IF (IERR/=MPI_SUCCESS) CALL HANDLE_MPI_ERROR('Error in MPI_ALLREDUCE for USUM_ALL',IERR)
   CALL MPI_ALLREDUCE(CONNECTED_ZONES_LOCAL(0,0),CONNECTED_ZONES_GLOBAL(0,0),(N_ZONE+1)**2,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,IERR)
   IF (IERR/=MPI_SUCCESS) CALL HANDLE_MPI_ERROR('Error in MPI_ALLREDUCE for CONNECTED_ZONES',IERR)
ELSE
   DSUM_ALL = DSUM_ALL_LOCAL
   PSUM_ALL = PSUM_ALL_LOCAL
   USUM_ALL = USUM_ALL_LOCAL
   CONNECTED_ZONES_GLOBAL = CONNECTED_ZONES_LOCAL
ENDIF

DO IPZ=1,N_ZONE
   IF(P_ZONE(IPZ)%EVACUATION) CYCLE
   DO NM=1,NMESHES
      IF(EVACUATION_ONLY(NM)) CYCLE
      DSUM(IPZ,NM) = DSUM_ALL(IPZ)
      PSUM(IPZ,NM) = PSUM_ALL(IPZ)
      USUM(IPZ,NM) = USUM_ALL(IPZ)
      CONNECTED_ZONES(IPZ,:,NM) = CONNECTED_ZONES_GLOBAL(IPZ,:)
   ENDDO
ENDDO

! Connect zones to others which are not directly connected

DO NM=1,NMESHES
   IF(EVACUATION_ONLY(NM)) CYCLE 
   DO IPZ=1,N_ZONE
      IF(P_ZONE(IPZ)%EVACUATION) CYCLE
      DO IOPZ=1,N_ZONE
         IF(P_ZONE(IOPZ)%EVACUATION) CYCLE
         IF (IOPZ==IPZ) CYCLE
         IF (CONNECTED_ZONES(IPZ,IOPZ,NM)) THEN
            DO IOPZ2=0,N_ZONE
               IF (IOPZ==IOPZ2) CYCLE
               IF (CONNECTED_ZONES(IOPZ,IOPZ2,NM)) CONNECTED_ZONES(IPZ,IOPZ2,NM) = .TRUE.
               IF (CONNECTED_ZONES(IOPZ,IOPZ2,NM)) CONNECTED_ZONES(IOPZ2,IPZ,NM) = .TRUE.
            ENDDO
         ENDIF
      ENDDO
   ENDDO
ENDDO

TUSED(2,:)=TUSED(2,:) + SECOND() - TNOW
END SUBROUTINE EXCHANGE_DIVERGENCE_INFO


SUBROUTINE INITIALIZE_MESH_EXCHANGE_1(NM)
 
! Create arrays by which info is to exchanged across meshes
 
INTEGER :: IMIN,IMAX,JMIN,JMAX,KMIN,KMAX,NOM,IOR,IW,N,N_STORAGE_SLOTS,IIO,JJO,KKO,NIC_R, &
           IIO_MIN,IIO_MAX,JJO_MIN,JJO_MAX,KKO_MIN,KKO_MAX
INTEGER, INTENT(IN) :: NM
TYPE (MESH_TYPE), POINTER :: M2,M
TYPE (LAGRANGIAN_PARTICLE_CLASS_TYPE), POINTER :: LPC
LOGICAL :: FOUND
 
M=>MESHES(NM)

NOT_EVACUATION_MESH_IF: IF (.NOT.EVACUATION_ONLY(NM)) THEN

ALLOCATE(MESHES(NM)%OMESH(NMESHES))
ALLOCATE(M%OMESH(NM)%IJKW(15,M%N_EXTERNAL_WALL_CELLS))

DO IW=1,M%N_EXTERNAL_WALL_CELLS
   M%OMESH(NM)%IJKW( 1,IW) = M%WALL(IW)%ONE_D%II
   M%OMESH(NM)%IJKW( 2,IW) = M%WALL(IW)%ONE_D%JJ
   M%OMESH(NM)%IJKW( 3,IW) = M%WALL(IW)%ONE_D%KK
   M%OMESH(NM)%IJKW( 4,IW) = M%WALL(IW)%ONE_D%IOR
   M%OMESH(NM)%IJKW( 5,IW) = M%WALL(IW)%SURF_INDEX
   M%OMESH(NM)%IJKW( 6,IW) = M%WALL(IW)%ONE_D%IIG
   M%OMESH(NM)%IJKW( 7,IW) = M%WALL(IW)%ONE_D%JJG
   M%OMESH(NM)%IJKW( 8,IW) = M%WALL(IW)%ONE_D%KKG
   M%OMESH(NM)%IJKW( 9,IW) = M%WALL(IW)%NOM
   M%OMESH(NM)%IJKW(10,IW) = M%WALL(IW)%NOM_IB(1)
   M%OMESH(NM)%IJKW(11,IW) = M%WALL(IW)%NOM_IB(2)
   M%OMESH(NM)%IJKW(12,IW) = M%WALL(IW)%NOM_IB(3)
   M%OMESH(NM)%IJKW(13,IW) = M%WALL(IW)%NOM_IB(4)
   M%OMESH(NM)%IJKW(14,IW) = M%WALL(IW)%NOM_IB(5)
   M%OMESH(NM)%IJKW(15,IW) = M%WALL(IW)%NOM_IB(6)
ENDDO

ALLOCATE(M%OMESH(NM)%BOUNDARY_TYPE(0:M%N_EXTERNAL_WALL_CELLS))
M%OMESH(NM)%BOUNDARY_TYPE(0) = 0
DO IW=1,M%N_EXTERNAL_WALL_CELLS
   M%OMESH(NM)%BOUNDARY_TYPE(IW) = M%WALL(IW)%BOUNDARY_TYPE
ENDDO

END IF NOT_EVACUATION_MESH_IF

OTHER_MESH_LOOP: DO NOM=1,NMESHES
 
   IF (EVACUATION_ONLY(NM)) THEN
      IF (EVACUATION_GRID(NM) .AND. .NOT.EVACUATION_ONLY(NOM)) N_COMMUNICATIONS = N_COMMUNICATIONS + 1
      CYCLE OTHER_MESH_LOOP
   ENDIF
   IF (EVACUATION_ONLY(NOM)) THEN
      IF (EVACUATION_GRID(NOM) .AND. .NOT.EVACUATION_ONLY(NM)) N_COMMUNICATIONS = N_COMMUNICATIONS + 1
      CYCLE OTHER_MESH_LOOP 
   ENDIF
 
   M2=>MESHES(NOM)
   IMIN=0 
   IMAX=M2%IBP1
   JMIN=0 
   JMAX=M2%JBP1
   KMIN=0 
   KMAX=M2%KBP1
   M%OMESH(NOM)%NIC_R = 0
   FOUND = .FALSE.

   SEARCH_LOOP: DO IW=1,M%N_EXTERNAL_WALL_CELLS
      IF (M%WALL(IW)%NOM/=NOM) CYCLE SEARCH_LOOP
      IIO_MIN = M%WALL(IW)%NOM_IB(1)
      JJO_MIN = M%WALL(IW)%NOM_IB(2)
      KKO_MIN = M%WALL(IW)%NOM_IB(3)
      IIO_MAX = M%WALL(IW)%NOM_IB(4)
      JJO_MAX = M%WALL(IW)%NOM_IB(5)
      KKO_MAX = M%WALL(IW)%NOM_IB(6)
      M%WALL(IW)%NIC_MIN = M%OMESH(NOM)%NIC_R + 1
      M%OMESH(NOM)%NIC_R = M%OMESH(NOM)%NIC_R + (IIO_MAX-IIO_MIN+1)*(JJO_MAX-JJO_MIN+1)*(KKO_MAX-KKO_MIN+1)
      M%WALL(IW)%NIC_MAX = M%OMESH(NOM)%NIC_R
      FOUND = .TRUE.
      CONNECTED_MESHES(NM,NOM) = .TRUE.
      CONNECTED_MESHES(NOM,NM) = .TRUE.
      IOR = M%WALL(IW)%ONE_D%IOR
      SELECT CASE(IOR)
         CASE( 1) 
            IMIN=MAX(IMIN,IIO_MIN-1)
         CASE(-1) 
            IMAX=MIN(IMAX,IIO_MAX+1)
         CASE( 2) 
            JMIN=MAX(JMIN,JJO_MIN-1)
         CASE(-2) 
            JMAX=MIN(JMAX,JJO_MAX+1)
         CASE( 3) 
            KMIN=MAX(KMIN,KKO_MIN-1)
         CASE(-3) 
            KMAX=MIN(KMAX,KKO_MAX+1)
      END SELECT
   ENDDO SEARCH_LOOP

   ! Allocate arrays to hold indices of arrays for MPI exchanges

   IF (M%OMESH(NOM)%NIC_R>0) THEN
      ALLOCATE(M%OMESH(NOM)%IIO_R(M%OMESH(NOM)%NIC_R))
      ALLOCATE(M%OMESH(NOM)%JJO_R(M%OMESH(NOM)%NIC_R))
      ALLOCATE(M%OMESH(NOM)%KKO_R(M%OMESH(NOM)%NIC_R))
      ALLOCATE(M%OMESH(NOM)%IOR_R(M%OMESH(NOM)%NIC_R))
      NIC_R = 0
      INDEX_LOOP: DO IW=1,M%N_EXTERNAL_WALL_CELLS
         IF (M%WALL(IW)%NOM/=NOM) CYCLE INDEX_LOOP
         DO KKO=M%WALL(IW)%NOM_IB(3),M%WALL(IW)%NOM_IB(6)
            DO JJO=M%WALL(IW)%NOM_IB(2),M%WALL(IW)%NOM_IB(5)
               DO IIO=M%WALL(IW)%NOM_IB(1),M%WALL(IW)%NOM_IB(4)
                  NIC_R = NIC_R + 1
                  IOR = M%WALL(IW)%ONE_D%IOR
                  M%OMESH(NOM)%IIO_R(NIC_R) = IIO
                  M%OMESH(NOM)%JJO_R(NIC_R) = JJO
                  M%OMESH(NOM)%KKO_R(NIC_R) = KKO
                  M%OMESH(NOM)%IOR_R(NIC_R) = IOR
               ENDDO
            ENDDO
         ENDDO
      ENDDO INDEX_LOOP
   ENDIF

   ! For PERIODIC boundaries with 1 or 2 meshes, we must revert to allocating whole copies of OMESH

   IF (IMIN>IMAX) THEN; IMIN=0; IMAX=M2%IBP1; ENDIF
   IF (JMIN>JMAX) THEN; JMIN=0; JMAX=M2%JBP1; ENDIF
   IF (KMIN>KMAX) THEN; KMIN=0; KMAX=M2%KBP1; ENDIF
 
   ! Embedded meshes

   IF ( NM/=NOM .AND. M2%XS>=M%XS .AND. M2%XF<=M%XF .AND. M2%YS>=M%YS .AND. M2%YF<=M%YF .AND. M2%ZS>=M%ZS .AND. M2%ZF<=M%ZF ) THEN
      FOUND = .TRUE.
      CONNECTED_MESHES(NM,NOM) = .TRUE.
      CONNECTED_MESHES(NOM,NM) = .TRUE.
   ENDIF

   ! Exit the other mesh loop if no neighboring meshes found

   IF (.NOT.FOUND) CYCLE OTHER_MESH_LOOP

   ! Tally the number of communications for this process

   N_COMMUNICATIONS = N_COMMUNICATIONS + 1
 
   ! Save the dimensions of the volume of cells from mesh NOM whose data is received by mesh NM

   M%OMESH(NOM)%I_MIN_R = IMIN
   M%OMESH(NOM)%I_MAX_R = IMAX
   M%OMESH(NOM)%J_MIN_R = JMIN
   M%OMESH(NOM)%J_MAX_R = JMAX
   M%OMESH(NOM)%K_MIN_R = KMIN
   M%OMESH(NOM)%K_MAX_R = KMAX

   ! Allocate the arrays that hold information about the other meshes (OMESH) 
 
   ALLOCATE(M%OMESH(NOM)% RHO(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
   ALLOCATE(M%OMESH(NOM)%RHOS(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
   M%OMESH(NOM)%RHO  = RHOA
   M%OMESH(NOM)%RHOS = RHOA
   ALLOCATE(M%OMESH(NOM)% D(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
   ALLOCATE(M%OMESH(NOM)%DS(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
   M%OMESH(NOM)%D  = 0._EB
   M%OMESH(NOM)%DS = 0._EB
   ALLOCATE(M%OMESH(NOM)%  MU(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
   M%OMESH(NOM)%MU = 0._EB
   ALLOCATE(M%OMESH(NOM)%    H(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
   ALLOCATE(M%OMESH(NOM)%   HS(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
   M%OMESH(NOM)%H  = 0._EB
   M%OMESH(NOM)%HS = 0._EB
   ALLOCATE(M%OMESH(NOM)%   U(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
   ALLOCATE(M%OMESH(NOM)%  US(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
   M%OMESH(NOM)%U  = U0
   M%OMESH(NOM)%US = U0
   ALLOCATE(M%OMESH(NOM)%   V(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
   ALLOCATE(M%OMESH(NOM)%  VS(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
   M%OMESH(NOM)%V  = V0
   M%OMESH(NOM)%VS = V0
   ALLOCATE(M%OMESH(NOM)%   W(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
   ALLOCATE(M%OMESH(NOM)%  WS(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
   M%OMESH(NOM)%W  = W0
   M%OMESH(NOM)%WS = W0
   ALLOCATE(M%OMESH(NOM)% FVX(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
   ALLOCATE(M%OMESH(NOM)% FVY(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
   ALLOCATE(M%OMESH(NOM)% FVZ(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
   M%OMESH(NOM)%FVX = 0._EB
   M%OMESH(NOM)%FVY = 0._EB
   M%OMESH(NOM)%FVZ = 0._EB
   ALLOCATE(M%OMESH(NOM)%KRES(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
   M%OMESH(NOM)%KRES = 0._EB
 
   ALLOCATE(M%OMESH(NOM)%  ZZ(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX,N_TOTAL_SCALARS))
   ALLOCATE(M%OMESH(NOM)% ZZS(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX,N_TOTAL_SCALARS))
   DO N=1,N_TRACKED_SPECIES
      M%OMESH(NOM)%ZZ(:,:,:,N)  = SPECIES_MIXTURE(N)%ZZ0
      M%OMESH(NOM)%ZZS(:,:,:,N) = SPECIES_MIXTURE(N)%ZZ0
   ENDDO
   DO N=N_TRACKED_SPECIES+1,N_TOTAL_SCALARS
      M%OMESH(NOM)%ZZ(:,:,:,N)  = INITIAL_UNMIXED_FRACTION
      M%OMESH(NOM)%ZZS(:,:,:,N) = INITIAL_UNMIXED_FRACTION
   ENDDO
  
   ! Wall arrays
   
   IF (.NOT.ALLOCATED(M%OMESH(NOM)%IJKW)) ALLOCATE(M%OMESH(NOM)%IJKW(15,M2%N_EXTERNAL_WALL_CELLS))
   IF (.NOT.ALLOCATED(M%OMESH(NOM)%BOUNDARY_TYPE)) ALLOCATE(M%OMESH(NOM)%BOUNDARY_TYPE(0:M2%N_EXTERNAL_WALL_CELLS))
   M%OMESH(NOM)%BOUNDARY_TYPE(0)=0

   ! Particle and PARTICLE Orphan Arrays
 
   IF (OMESH_PARTICLES) THEN
      ALLOCATE(M%OMESH(NOM)%N_PART_ORPHANS(N_LAGRANGIAN_CLASSES))
      ALLOCATE(M%OMESH(NOM)%N_PART_ADOPT(N_LAGRANGIAN_CLASSES))
      M%OMESH(NOM)%N_PART_ORPHANS = 0
      M%OMESH(NOM)%N_PART_ADOPT   = 0
      ALLOCATE(M%OMESH(NOM)%ORPHAN_PARTICLE_STORAGE(N_LAGRANGIAN_CLASSES))
      ALLOCATE(M%OMESH(NOM)%ADOPT_PARTICLE_STORAGE(N_LAGRANGIAN_CLASSES))
      DO N=1,N_LAGRANGIAN_CLASSES
         LPC => LAGRANGIAN_PARTICLE_CLASS(N)
         N_STORAGE_SLOTS = 1000
         M%OMESH(NOM)%ORPHAN_PARTICLE_STORAGE(N)%N_STORAGE_SLOTS = N_STORAGE_SLOTS
         M%OMESH(NOM)%ADOPT_PARTICLE_STORAGE(N)%N_STORAGE_SLOTS = N_STORAGE_SLOTS
         ALLOCATE(M%OMESH(NOM)%ORPHAN_PARTICLE_STORAGE(N)%REALS(LPC%N_STORAGE_REALS,N_STORAGE_SLOTS))
         ALLOCATE(M%OMESH(NOM)%ORPHAN_PARTICLE_STORAGE(N)%INTEGERS(LPC%N_STORAGE_INTEGERS,N_STORAGE_SLOTS))
         ALLOCATE(M%OMESH(NOM)%ORPHAN_PARTICLE_STORAGE(N)%LOGICALS(LPC%N_STORAGE_LOGICALS,N_STORAGE_SLOTS))
         ALLOCATE(M%OMESH(NOM)%ADOPT_PARTICLE_STORAGE(N)%REALS(LPC%N_STORAGE_REALS,N_STORAGE_SLOTS))
         ALLOCATE(M%OMESH(NOM)%ADOPT_PARTICLE_STORAGE(N)%INTEGERS(LPC%N_STORAGE_INTEGERS,N_STORAGE_SLOTS))
         ALLOCATE(M%OMESH(NOM)%ADOPT_PARTICLE_STORAGE(N)%LOGICALS(LPC%N_STORAGE_LOGICALS,N_STORAGE_SLOTS))
      ENDDO
   ENDIF

ENDDO OTHER_MESH_LOOP

END SUBROUTINE INITIALIZE_MESH_EXCHANGE_1


SUBROUTINE INITIALIZE_MESH_EXCHANGE_2(NM)

! Create arrays by which info is to exchanged across meshes. In this routine, allocate arrays that involve NIC_R and NIC_S arrays.

INTEGER :: NOM
INTEGER, INTENT(IN) :: NM
TYPE (MESH_TYPE), POINTER :: M

IF (EVACUATION_ONLY(NM)) RETURN

M=>MESHES(NM)

! Allocate arrays to send (IL_S) and receive (IL_R) the radiation intensity (IL) at interpolated boundaries.
! MESHES(NM)%OMESH(NOM)%IL_S are the intensities in mesh NM that are just outside the boundary of mesh NOM. IL_S is populated
! in radi.f90 and then sent to MESHES(NOM)%OMESH(NM)%IL_R in MESH_EXCHANGE. IL_R holds the intensities until they are 
! transferred to the ghost cells of MESHES(NOM)%IL in radi.f90. The IL_S and IL_R arrays are indexed by NIC_S and NIC_R.

DO NOM=1,NMESHES
   IF (M%OMESH(NOM)%NIC_S>0) THEN
      ALLOCATE(M%OMESH(NOM)%IL_S(M%OMESH(NOM)%NIC_S,NUMBER_RADIATION_ANGLES,NUMBER_SPECTRAL_BANDS))
      M%OMESH(NOM)%IL_S = RPI*SIGMA*TMPA4
    ENDIF
   IF (M%OMESH(NOM)%NIC_R>0) THEN
      ALLOCATE(M%OMESH(NOM)%IL_R(M%OMESH(NOM)%NIC_R,NUMBER_RADIATION_ANGLES,NUMBER_SPECTRAL_BANDS))
      M%OMESH(NOM)%IL_R = RPI*SIGMA*TMPA4
   ENDIF
ENDDO

END SUBROUTINE INITIALIZE_MESH_EXCHANGE_2


SUBROUTINE INITIALIZE_BACK_WALL_EXCHANGE

! Bordering meshes swap the number of exposed back wall cells. If the meshes have the same resolution, this number is the same.

CALL POST_RECEIVES(8)
CALL MESH_EXCHANGE(8)
CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

! Create an array of indices that maps the actual wall cell index (IW) to a shorter list of just those that are split across
! mesh boundaries. MESHES(NM)%OMESH(NOM)%EXT_BACK_WALL_CELL_INDEX(I) is the WALL_INDEX (IW) of the Ith listed wall cell in MESH
! NOM. These cells are the back sides of wall cells in MESH NM.

DO NM=1,NMESHES
   IF (PROCESS(NM)/=MYID .OR. EVACUATION_ONLY(NM)) CYCLE
   DO NOM=1,NMESHES
      M3 => MESHES(NM)%OMESH(NOM)
      IF (M3%N_EXT_BACK_WALL_CELLS>0) ALLOCATE(M3%EXT_BACK_WALL_CELL_INDEX(M3%N_EXT_BACK_WALL_CELLS))
   ENDDO
ENDDO
CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

! Exchange BACK_WALL_CELL_INDEX. Note that M3%BACK_WALL_CELL_INDEX is sent from MESH NM to NOM, where it is referred to as
! EXT_BACK_WALL_CELL_INDEX.

CALL POST_RECEIVES(9)
CALL MESH_EXCHANGE(9)
CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

! Allocate arrays for sending and receiving BACK_WALL info. Also, for MESH NM, loop over all wall cells (IW) and if a wall cell
! has a back side in another mesh, NOM, (M%WALL(IW)%BACK_MESH==NOM), loop over the M3%N_EXT_BACK_WALL_CELLS wall cells looking
! for a M3%EXT_BACK_WALL_CELL_INDEX that matches M%WALL(IW)%BACK_INDEX. When found, reassign the M%WALL(IW)%BACK_INDEX to I, the 
! index of the short list of BACK_WALL cells.

MESH_LOOP_1: DO NM=1,NMESHES
   IF (PROCESS(NM)/=MYID .OR. EVACUATION_ONLY(NM)) CYCLE MESH_LOOP_1
   M => MESHES(NM)
   MESH_LOOP_2: DO NOM=1,NMESHES
      IF (NM==NOM .OR. EVACUATION_ONLY(NOM)) CYCLE MESH_LOOP_2
      M3 => M%OMESH(NOM)
      IF (M3%N_EXT_BACK_WALL_CELLS>0) THEN
         IF (M3%N_BACK_WALL_CELLS/=M3%N_EXT_BACK_WALL_CELLS) THEN
            WRITE(0,'(A,I2,A,I2)') 'ERROR: Mismatch of back wall cells between MESH ',NM,', and MESH ',NOM
            STOP_STATUS = SETUP_STOP
            EXIT MESH_LOOP_1
         ENDIF
         ALLOCATE(M3%REAL_SEND_PKG6(M3%N_BACK_WALL_CELLS*2))
         ALLOCATE(M3%REAL_RECV_PKG6(M3%N_EXT_BACK_WALL_CELLS*2))
         ALLOCATE(M3%BACK_WALL(M3%N_EXT_BACK_WALL_CELLS))
         WALL_LOOP: DO IW=1,M%N_EXTERNAL_WALL_CELLS+M%N_INTERNAL_WALL_CELLS
            IF (M%WALL(IW)%BACK_MESH/=NOM) CYCLE WALL_LOOP
            FOUND = .FALSE.
            BACK_WALL_LOOP: DO I=1,M3%N_EXT_BACK_WALL_CELLS
               IF (M%WALL(IW)%BACK_INDEX==M3%EXT_BACK_WALL_CELL_INDEX(I)) THEN
                  M%WALL(IW)%BACK_INDEX = I
                  FOUND = .TRUE.
                  EXIT BACK_WALL_LOOP
               ENDIF
            ENDDO BACK_WALL_LOOP
            IF (.NOT.FOUND) THEN
               WRITE(0,'(A,I2,A,I2)') 'ERROR: Misalignment of obstruction between MESH ',NM,', and MESH ',NOM
               STOP_STATUS = SETUP_STOP
               EXIT MESH_LOOP_1
            ENDIF
         ENDDO WALL_LOOP
      ENDIF
   ENDDO MESH_LOOP_2
ENDDO MESH_LOOP_1

! Check to see if any process has an error. If so, stop the run.

CALL STOP_CHECK(1)

! Set up persistent SEND and RECV calls for BACK_WALL info

CALL POST_RECEIVES(10)
CALL MESH_EXCHANGE(10)
CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

END SUBROUTINE INITIALIZE_BACK_WALL_EXCHANGE


SUBROUTINE POST_RECEIVES(CODE)

! Set up receive buffers for MPI calls.

INTEGER, INTENT(IN) :: CODE
INTEGER :: RNODE,SNODE,IJK_SIZE,N,N_STORAGE_SLOTS,NRA,NRA_MAX,LL,AIC
REAL(EB) :: TNOW
TYPE (LAGRANGIAN_PARTICLE_CLASS_TYPE), POINTER :: LPC

TNOW = SECOND()

! Initialize the number of non-persistent send/receive requests.

N_REQ = 0

! Loop over all receive meshes (NM) and look for the send meshes (NOM). 

MESH_LOOP: DO NM=1,NMESHES

   IF (EVACUATION_ONLY(NM)) CYCLE MESH_LOOP 

   RNODE = PROCESS(NM)
   IF (RNODE/=MYID) CYCLE MESH_LOOP
   M => MESHES(NM)
 
   OTHER_MESH_LOOP: DO NOM=1,NMESHES
    
      IF (EVACUATION_ONLY(NOM)) CYCLE OTHER_MESH_LOOP
      IF (CODE>0 .AND. (.NOT.ACTIVE_MESH(NOM).OR..NOT.ACTIVE_MESH(NM))) CYCLE OTHER_MESH_LOOP
   
      SNODE = PROCESS(NOM)
      IF (RNODE==SNODE) CYCLE OTHER_MESH_LOOP
   
      M4=>MESHES(NOM)
      M3=>MESHES(NM)%OMESH(NOM)
   
      IF (M3%NIC_S==0 .AND. M3%NIC_R==0) CYCLE OTHER_MESH_LOOP
    
      ! Set up receives for one-time exchanges or persistent send/receives.
   
      INITIALIZATION_IF: IF (CODE==0) THEN
   
         IF (.NOT.ALLOCATED(M4%CELL_INDEX)) ALLOCATE(M4%CELL_INDEX(0:M4%IBP1,0:M4%JBP1,0:M4%KBP1))
         IF (.NOT.ALLOCATED(M4%SOLID))      ALLOCATE(M4%SOLID(0:CELL_COUNT(NOM)))
         IF (.NOT.ALLOCATED(M4%WALL_INDEX)) ALLOCATE(M4%WALL_INDEX(0:CELL_COUNT(NOM),-3:3))
   
         N_REQ = MIN(N_REQ+1,SIZE(REQ))
         CALL MPI_IRECV(M4%CELL_INDEX(0,0,0),SIZE(M4%CELL_INDEX),MPI_INTEGER,SNODE,NOM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
         N_REQ = MIN(N_REQ+1,SIZE(REQ))
         CALL MPI_IRECV(M4%SOLID(0),SIZE(M4%SOLID),MPI_INTEGER,SNODE,NOM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
         N_REQ = MIN(N_REQ+1,SIZE(REQ))
         CALL MPI_IRECV(M4%WALL_INDEX(0,-3),SIZE(M4%WALL_INDEX),MPI_INTEGER,SNODE,NOM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
    
         IJK_SIZE = (M3%I_MAX_R-M3%I_MIN_R+1)*(M3%J_MAX_R-M3%J_MIN_R+1)*(M3%K_MAX_R-M3%K_MIN_R+1)
   
         IF (M3%NIC_R>0) THEN

            ! Determine the maximum number of radiation angles that are to be received

            NRA_MAX = 0
            IF (RADIATION) THEN
               DO AIC=1,ANGLE_INCREMENT
                  DO LL=1,M3%NIC_R
                     NRA = 0
                     DO N=NUMBER_RADIATION_ANGLES-AIC+1,1,-ANGLE_INCREMENT
                        IF (DLN(M3%IOR_R(LL),N)>0._EB) NRA = NRA + 1
                     ENDDO
                     NRA_MAX = MAX(NRA_MAX,NRA)
                  ENDDO
               ENDDO
            ENDIF

            ! Allocate the 1-D arrays that hold the big mesh variables that are to be received

            ALLOCATE(M3%REAL_RECV_PKG1(M3%NIC_R*2*(4+N_TOTAL_SCALARS)))
            ALLOCATE(M3%REAL_RECV_PKG2(IJK_SIZE*4))
            ALLOCATE(M3%REAL_RECV_PKG3(M3%NIC_R*2*(4+N_TOTAL_SCALARS)))
            ALLOCATE(M3%REAL_RECV_PKG4(IJK_SIZE*4))
            ALLOCATE(M3%REAL_RECV_PKG5(NRA_MAX*NUMBER_SPECTRAL_BANDS*M3%NIC_R))
            ALLOCATE(M3%REAL_RECV_PKG7(M3%NIC_R*3))
         ENDIF
    
         N_REQ = MIN(N_REQ+1,SIZE(REQ))
         CALL MPI_IRECV(M3%IJKW(1,1),SIZE(M3%IJKW),MPI_INTEGER,SNODE,NOM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
   
         ! Set up persistent receive requests
    
         IF (M3%NIC_R>0) THEN
   
            N_REQ1 = N_REQ1 + 1
            CALL MPI_RECV_INIT(M3%REAL_RECV_PKG1(1),SIZE(M3%REAL_RECV_PKG1),MPI_DOUBLE_PRECISION,SNODE,NOM,MPI_COMM_WORLD,&
                               REQ1(N_REQ1),IERR)
   
            IF (OMESH_PARTICLES) THEN
               N_REQ2 = N_REQ2 + 1
               CALL MPI_RECV_INIT(M3%N_PART_ADOPT,SIZE(M3%N_PART_ADOPT),MPI_INTEGER,SNODE,NOM,MPI_COMM_WORLD,&
                                  REQ2(N_REQ2),IERR)
            ENDIF
   
            N_REQ3 = N_REQ3 + 1
            CALL MPI_RECV_INIT(M3%REAL_RECV_PKG2(1),SIZE(M3%REAL_RECV_PKG2),MPI_DOUBLE_PRECISION,SNODE,NOM,MPI_COMM_WORLD,&
                               REQ3(N_REQ3),IERR)
   
            N_REQ4 = N_REQ4 + 1
            CALL MPI_RECV_INIT(M3%REAL_RECV_PKG3(1),SIZE(M3%REAL_RECV_PKG3),MPI_DOUBLE_PRECISION,SNODE,NOM,MPI_COMM_WORLD,&
                               REQ4(N_REQ4),IERR)
   
            N_REQ5 = N_REQ5 + 1
            CALL MPI_RECV_INIT(M3%REAL_RECV_PKG7(1),SIZE(M3%REAL_RECV_PKG7),MPI_DOUBLE_PRECISION,SNODE,NOM,MPI_COMM_WORLD,&
                               REQ5(N_REQ5),IERR)
   
            N_REQ7 = N_REQ7 + 1
            CALL MPI_RECV_INIT(M3%REAL_RECV_PKG4(1),SIZE(M3%REAL_RECV_PKG4),MPI_DOUBLE_PRECISION,SNODE,NOM,MPI_COMM_WORLD,&
                               REQ7(N_REQ7),IERR)
   
            N_REQ8 = N_REQ8 + 1      
            CALL MPI_RECV_INIT(M3%BOUNDARY_TYPE(0),SIZE(M3%BOUNDARY_TYPE),MPI_INTEGER,SNODE,NOM,MPI_COMM_WORLD,&
                               REQ8(N_REQ8),IERR)
   
            IF (RADIATION) THEN
               N_REQ9 = N_REQ9 + 1
               CALL MPI_RECV_INIT(M3%REAL_RECV_PKG5(1),SIZE(M3%REAL_RECV_PKG5),MPI_DOUBLE_PRECISION,SNODE,NOM,MPI_COMM_WORLD,&
                                  REQ9(N_REQ9),IERR)
            ENDIF
   
         ENDIF
   
      ENDIF INITIALIZATION_IF
   
      ! Exchange BACK_WALL information
   
      IF (CODE==8) THEN
         N_REQ=MIN(N_REQ+1,SIZE(REQ))
         CALL MPI_IRECV(M3%N_EXT_BACK_WALL_CELLS,1,MPI_INTEGER,SNODE,NOM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
      ENDIF
   
      IF (CODE==9 .AND. M3%N_EXT_BACK_WALL_CELLS>0) THEN
            N_REQ=MIN(N_REQ+1,SIZE(REQ))
            CALL MPI_IRECV(M3%EXT_BACK_WALL_CELL_INDEX,SIZE(M3%EXT_BACK_WALL_CELL_INDEX),MPI_INTEGER,SNODE,NOM,MPI_COMM_WORLD,&
                           REQ(N_REQ),IERR)
      ENDIF
   
      IF (CODE==10 .AND. M3%N_EXT_BACK_WALL_CELLS>0) THEN
         N_REQ6 = N_REQ6 + 1     
         CALL MPI_RECV_INIT(M3%REAL_RECV_PKG6(1),SIZE(M3%REAL_RECV_PKG6),MPI_DOUBLE_PRECISION,SNODE,NOM,MPI_COMM_WORLD,&
                            REQ6(N_REQ6),IERR)
      ENDIF
   
      ! PARTICLEs
    
      IF (CODE==6 .AND. OMESH_PARTICLES) THEN
         DO N=1,N_LAGRANGIAN_CLASSES
            IF (M3%N_PART_ADOPT(N)==0) CYCLE
            LPC => LAGRANGIAN_PARTICLE_CLASS(N)
            N_REQ=MIN(N_REQ+1,SIZE(REQ))
            N_STORAGE_SLOTS = M3%ADOPT_PARTICLE_STORAGE(N)%N_STORAGE_SLOTS
            CALL MPI_IRECV(M3%ADOPT_PARTICLE_STORAGE(N)%REALS(1,1),LPC%N_STORAGE_REALS*N_STORAGE_SLOTS, &
                           MPI_DOUBLE_PRECISION,SNODE,NOM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
            N_REQ=MIN(N_REQ+1,SIZE(REQ))
            CALL MPI_IRECV(M3%ADOPT_PARTICLE_STORAGE(N)%INTEGERS(1,1),LPC%N_STORAGE_INTEGERS*N_STORAGE_SLOTS, &
                           MPI_INTEGER,SNODE,NOM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
            N_REQ=MIN(N_REQ+1,SIZE(REQ))
            CALL MPI_IRECV(M3%ADOPT_PARTICLE_STORAGE(N)%LOGICALS(1,1),LPC%N_STORAGE_LOGICALS*N_STORAGE_SLOTS, &
                           MPI_LOGICAL,SNODE,NOM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
         ENDDO
      ENDIF
   
   ENDDO OTHER_MESH_LOOP
   
ENDDO MESH_LOOP

! Receive EVACuation information

DO NOM=1,NMESHES
   SNODE = PROCESS(NOM)
   IF (CODE==6 .AND. EXCHANGE_EVACUATION .AND. MYID==MAX(0,EVAC_PROCESS) .AND. .NOT.EVACUATION_ONLY(NOM)) THEN
      M4=>MESHES(NOM)
      TAG_EVAC = NOM*(MAX(0,EVAC_PROCESS)+1)*CODE*10
      IWW = (M4%IBAR+2)*(M4%JBAR+2)*(M4%KBAR+2)
      N_REQ=MIN(N_REQ+1,SIZE(REQ))
      CALL MPI_IRECV(M4%ZZ(0,0,0,1),IWW*N_TRACKED_SPECIES,MPI_DOUBLE_PRECISION,SNODE,TAG_EVAC,MPI_COMM_WORLD,REQ(N_REQ),IERR)
      N_REQ=MIN(N_REQ+1,SIZE(REQ))
      CALL MPI_IRECV(M4%RHO(0,0,0),IWW,MPI_DOUBLE_PRECISION,SNODE,TAG_EVAC,MPI_COMM_WORLD,REQ(N_REQ),IERR)
      N_REQ=MIN(N_REQ+1,SIZE(REQ))
      CALL MPI_IRECV(M4%RSUM(0,0,0),IWW,MPI_DOUBLE_PRECISION,SNODE,TAG_EVAC,MPI_COMM_WORLD,REQ(N_REQ),IERR)
      N_REQ=MIN(N_REQ+1,SIZE(REQ))
      CALL MPI_IRECV(M4%TMP(0,0,0),IWW,MPI_DOUBLE_PRECISION,SNODE,TAG_EVAC,MPI_COMM_WORLD,REQ(N_REQ),IERR)
      N_REQ=MIN(N_REQ+1,SIZE(REQ))
      CALL MPI_IRECV(M4%UII(0,0,0),IWW,MPI_DOUBLE_PRECISION,SNODE,TAG_EVAC,MPI_COMM_WORLD,REQ(N_REQ),IERR)
      N_REQ=MIN(N_REQ+1,SIZE(REQ))
      CALL MPI_IRECV(M4%CELL_INDEX(0,0,0),IWW,MPI_INTEGER,SNODE,TAG_EVAC,MPI_COMM_WORLD,REQ(N_REQ),IERR)
      IWW = MAXVAL(M4%CELL_INDEX)
      N_REQ=MIN(N_REQ+1,SIZE(REQ))
      CALL MPI_IRECV(M4%SOLID(0),IWW,MPI_LOGICAL,SNODE,TAG_EVAC,MPI_COMM_WORLD,REQ(N_REQ),IERR)
   ENDIF
ENDDO
 
TUSED(11,:)=TUSED(11,:) + SECOND() - TNOW
END SUBROUTINE POST_RECEIVES


SUBROUTINE MESH_EXCHANGE(CODE)
 
! Exchange Information between Meshes
 
REAL(EB) :: TNOW
INTEGER, INTENT(IN) :: CODE
INTEGER :: NM,II,JJ,KK,LL,LLL,N,RNODE,SNODE,IMIN,IMAX,JMIN,JMAX,KMIN,KMAX,IJK_SIZE,N_STORAGE_SLOTS,N_NEW_STORAGE_SLOTS
INTEGER :: NN1,NN2,IPC,CNT,IBC,STORAGE_INDEX_SAVE,II1,II2,JJ1,JJ2,KK1,KK2,NQT2,NN,IOR,NRA,NRA_MAX,AIC
REAL(EB), POINTER, DIMENSION(:,:,:) :: HP,HP2
TYPE (LAGRANGIAN_PARTICLE_TYPE), POINTER :: LP
TYPE (LAGRANGIAN_PARTICLE_CLASS_TYPE), POINTER :: LPC

TNOW = SECOND()

! Special circumstances when doing the radiation exchange (CODE=2)

IF (CODE==2) THEN
   IF (ANY(EVACUATION_ONLY) .AND. N_MPI_PROCESSES>1 .AND. RADIATION) THEN
      IF (.NOT.EXCHANGE_RADIATION .AND. MYID/=EVAC_PROCESS) THEN
         CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
         RETURN
      ENDIF
   ELSE
      IF (.NOT.EXCHANGE_RADIATION .OR. .NOT.RADIATION) RETURN
   ENDIF
ENDIF

SENDING_MESH_LOOP: DO NM=1,NMESHES

   IF (PROCESS(NM)/=MYID)   CYCLE SENDING_MESH_LOOP
   IF (EVACUATION_ONLY(NM)) CYCLE SENDING_MESH_LOOP 

   M =>MESHES(NM)
   M5=>MESHES(NM)%OMESH(NM)

   ! Information about Mesh NM is packed into SEND packages and shipped out to the other meshes (machines) via MPI
 
   RECEIVING_MESH_LOOP: DO NOM=1,NMESHES
 
      SNODE = PROCESS(NOM)
      RNODE = PROCESS(NM)

      IF (EVACUATION_ONLY(NOM)) CYCLE RECEIVING_MESH_LOOP 

      M3=>MESHES(NM)%OMESH(NOM)
      M4=>MESHES(NOM)

      IF (M3%NIC_S==0 .AND. M3%NIC_R==0)  CYCLE RECEIVING_MESH_LOOP
 
      IF (CODE>0) THEN
         IF (.NOT.ACTIVE_MESH(NM) .OR. .NOT.ACTIVE_MESH(NOM))  CYCLE RECEIVING_MESH_LOOP
      ENDIF
 
      IMIN = M3%I_MIN_S
      IMAX = M3%I_MAX_S
      JMIN = M3%J_MIN_S
      JMAX = M3%J_MAX_S
      KMIN = M3%K_MIN_S
      KMAX = M3%K_MAX_S

      IJK_SIZE = (IMAX-IMIN+1)*(JMAX-JMIN+1)*(KMAX-KMIN+1)

      ! Set up sends for one-time exchanges or persistent send/receives.

      INITIALIZE_SEND_IF: IF (CODE==0) THEN

         IF (RNODE/=SNODE) THEN
            N_REQ=MIN(N_REQ+1,SIZE(REQ))
            CALL MPI_ISEND(M%CELL_INDEX(0,0,0),SIZE(M%CELL_INDEX),MPI_INTEGER,SNODE,NM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
            N_REQ=MIN(N_REQ+1,SIZE(REQ))
            CALL MPI_ISEND(M%SOLID(0),SIZE(M%SOLID),MPI_INTEGER,SNODE,NM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
            N_REQ=MIN(N_REQ+1,SIZE(REQ))
            CALL MPI_ISEND(M%WALL_INDEX(0,-3),SIZE(M%WALL_INDEX),MPI_INTEGER,SNODE,NM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
         ENDIF

         IF (M3%NIC_S>0 .AND. RNODE/=SNODE) THEN

            ! Determine the maximum number of radiation angles that are to be sent

            NRA_MAX = 0
            IF (RADIATION) THEN
               DO AIC=1,ANGLE_INCREMENT
                  DO LL=1,M3%NIC_S
                     NRA = 0
                     DO N=NUMBER_RADIATION_ANGLES-AIC+1,1,-ANGLE_INCREMENT
                        IF (DLN(M3%IOR_S(LL),N)>0._EB) NRA = NRA + 1
                     ENDDO
                     NRA_MAX = MAX(NRA_MAX,NRA)
                  ENDDO
               ENDDO
            ENDIF

            ! Allocate 1-D arrays to hold major mesh variables that are to be sent to neighboring meshes

            ALLOCATE(M3%REAL_SEND_PKG1(M3%NIC_S*2*(4+N_TOTAL_SCALARS)))
            ALLOCATE(M3%REAL_SEND_PKG2(IJK_SIZE*4))
            ALLOCATE(M3%REAL_SEND_PKG3(M3%NIC_S*2*(4+N_TOTAL_SCALARS)))
            ALLOCATE(M3%REAL_SEND_PKG4(IJK_SIZE*4))
            ALLOCATE(M3%REAL_SEND_PKG5(NRA_MAX*NUMBER_SPECTRAL_BANDS*M3%NIC_S))
            ALLOCATE(M3%REAL_SEND_PKG7(M3%NIC_S*3))

         ENDIF
 
         IF (RNODE/=SNODE) THEN
            N_REQ=MIN(N_REQ+1,SIZE(REQ))
            CALL MPI_ISEND(M%OMESH(NM)%IJKW(1,1),SIZE(M%OMESH(NM)%IJKW),MPI_INTEGER,SNODE,NM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
         ELSE
            M%OMESH(NOM)%IJKW = M4%OMESH(NOM)%IJKW(:,1:M4%N_EXTERNAL_WALL_CELLS)
         ENDIF
 
         ! Initialize persistent send requests

         IF (M3%NIC_S>0 .AND. RNODE/=SNODE) THEN

            N_REQ1 = N_REQ1 + 1
            CALL MPI_SEND_INIT(M3%REAL_SEND_PKG1(1),SIZE(M3%REAL_SEND_PKG1),MPI_DOUBLE_PRECISION,SNODE,NM,MPI_COMM_WORLD,&
                               REQ1(N_REQ1),IERR)

            IF (OMESH_PARTICLES) THEN
               N_REQ2 = N_REQ2 + 1
               CALL MPI_SEND_INIT(M3%N_PART_ORPHANS,SIZE(M3%N_PART_ORPHANS),MPI_INTEGER,SNODE,NM,MPI_COMM_WORLD,&
                                  REQ2(N_REQ2),IERR)
            ENDIF

            N_REQ3 = N_REQ3 + 1
            CALL MPI_SEND_INIT(M3%REAL_SEND_PKG2(1),SIZE(M3%REAL_SEND_PKG2),MPI_DOUBLE_PRECISION,SNODE,NM,MPI_COMM_WORLD,&
                               REQ3(N_REQ3),IERR)

            N_REQ4 = N_REQ4 + 1
            CALL MPI_SEND_INIT(M3%REAL_SEND_PKG3(1),SIZE(M3%REAL_SEND_PKG3),MPI_DOUBLE_PRECISION,SNODE,NM,MPI_COMM_WORLD,&
                               REQ4(N_REQ4),IERR)

            N_REQ5 = N_REQ5 + 1
            CALL MPI_SEND_INIT(M3%REAL_SEND_PKG7(1),SIZE(M3%REAL_SEND_PKG7),MPI_DOUBLE_PRECISION,SNODE,NM,MPI_COMM_WORLD,&
                               REQ5(N_REQ5),IERR)

            N_REQ7 = N_REQ7 + 1
            CALL MPI_SEND_INIT(M3%REAL_SEND_PKG4(1),SIZE(M3%REAL_SEND_PKG4),MPI_DOUBLE_PRECISION,SNODE,NM,MPI_COMM_WORLD, &
                               REQ7(N_REQ7),IERR)

            N_REQ8 = N_REQ8 + 1
            CALL MPI_SEND_INIT(M5%BOUNDARY_TYPE(0),SIZE(M5%BOUNDARY_TYPE),MPI_INTEGER,SNODE,NM,MPI_COMM_WORLD,&
                               REQ8(N_REQ8),IERR)

            IF (RADIATION) THEN
               N_REQ9 = N_REQ9 + 1
               CALL MPI_SEND_INIT(M3%REAL_SEND_PKG5(1),SIZE(M3%REAL_SEND_PKG5),MPI_DOUBLE_PRECISION,SNODE,NM,MPI_COMM_WORLD,&
                                  REQ9(N_REQ9),IERR)
            ENDIF

         ENDIF

      ENDIF INITIALIZE_SEND_IF

      ! Exchange the number of solid surface cells whose back side is in another mesh

      IF (CODE==8) THEN
         IF (RNODE/=SNODE) THEN
            N_REQ=MIN(N_REQ+1,SIZE(REQ))
            CALL MPI_ISEND(M3%N_BACK_WALL_CELLS,1,MPI_INTEGER,SNODE,NM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
         ELSE
            M2=>MESHES(NOM)%OMESH(NM)
            M2%N_EXT_BACK_WALL_CELLS = M3%N_BACK_WALL_CELLS
         ENDIF
      ENDIF

      IF (CODE==9 .AND. M3%N_BACK_WALL_CELLS>0) THEN
         IF (RNODE/=SNODE) THEN
            N_REQ=MIN(N_REQ+1,SIZE(REQ))
            CALL MPI_ISEND(M3%BACK_WALL_CELL_INDEX,SIZE(M3%BACK_WALL_CELL_INDEX),MPI_INTEGER,SNODE,NM,MPI_COMM_WORLD,&
                           REQ(N_REQ),IERR)
         ELSE
            M2=>MESHES(NOM)%OMESH(NM)
            M2%EXT_BACK_WALL_CELL_INDEX = M3%BACK_WALL_CELL_INDEX
         ENDIF
      ENDIF

      IF (CODE==10 .AND. M3%N_BACK_WALL_CELLS>0) THEN
         IF (RNODE/=SNODE) THEN
            N_REQ6 = N_REQ6 + 1
            CALL MPI_SEND_INIT(M3%REAL_SEND_PKG6(1),SIZE(M3%REAL_SEND_PKG6),MPI_DOUBLE_PRECISION,SNODE,NM,MPI_COMM_WORLD,&
                                REQ6(N_REQ6),IERR)
         ENDIF
      ENDIF

      ! Exchange of density and species mass fractions following the PREDICTOR update

      IF (CODE==1 .AND. M3%NIC_S>0) THEN
         IF (RNODE/=SNODE) THEN
            NQT2 = 2*(4+N_TOTAL_SCALARS)
            PACK_REAL_SEND_PKG1: DO LL=1,M3%NIC_S
               II1 = M3%IIO_S(LL) ; II2 = II1
               JJ1 = M3%JJO_S(LL) ; JJ2 = JJ1
               KK1 = M3%KKO_S(LL) ; KK2 = KK1
               SELECT CASE(M3%IOR_S(LL))
                  CASE(-1) ; II1=M3%IIO_S(LL)   ; II2=II1+1 
                  CASE( 1) ; II1=M3%IIO_S(LL)-1 ; II2=II1+1 
                  CASE(-2) ; JJ1=M3%JJO_S(LL)   ; JJ2=JJ1+1 
                  CASE( 2) ; JJ1=M3%JJO_S(LL)-1 ; JJ2=JJ1+1 
                  CASE(-3) ; KK1=M3%KKO_S(LL)   ; KK2=KK1+1 
                  CASE( 3) ; KK1=M3%KKO_S(LL)-1 ; KK2=KK1+1 
               END SELECT
               M3%REAL_SEND_PKG1(NQT2*(LL-1)+1) = M%RHOS(II1,JJ1,KK1)
               M3%REAL_SEND_PKG1(NQT2*(LL-1)+2) = M%RHOS(II2,JJ2,KK2)
               M3%REAL_SEND_PKG1(NQT2*(LL-1)+3) =   M%MU(II1,JJ1,KK1)
               M3%REAL_SEND_PKG1(NQT2*(LL-1)+4) =   M%MU(II2,JJ2,KK2)
               M3%REAL_SEND_PKG1(NQT2*(LL-1)+5) = M%KRES(II1,JJ1,KK1)
               M3%REAL_SEND_PKG1(NQT2*(LL-1)+6) = M%KRES(II2,JJ2,KK2)
               M3%REAL_SEND_PKG1(NQT2*(LL-1)+7) =    M%D(II1,JJ1,KK1)
               M3%REAL_SEND_PKG1(NQT2*(LL-1)+8) =    M%D(II2,JJ2,KK2)
               DO NN=1,N_TOTAL_SCALARS
                  M3%REAL_SEND_PKG1(NQT2*(LL-1)+8+2*NN-1) = M%ZZS(II1,JJ1,KK1,NN)
                  M3%REAL_SEND_PKG1(NQT2*(LL-1)+8+2*NN  ) = M%ZZS(II2,JJ2,KK2,NN)
               ENDDO
            ENDDO PACK_REAL_SEND_PKG1
         ELSE
            M2=>MESHES(NOM)%OMESH(NM)
            M2%RHOS(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)= M%RHOS(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
            M2%MU(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)  = M%MU(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
            M2%KRES(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)= M%KRES(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
            M2%D(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)   = M%D(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
            M2%ZZS(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX,1:N_TOTAL_SCALARS)= &
                                  M%ZZS(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX,1:N_TOTAL_SCALARS)
         ENDIF
      ENDIF

      ! Exchange velocity/pressure info for ITERATE_PRESSURE

      IF (CODE==5 .AND. M3%NIC_S>0) THEN
         IF (PREDICTOR) HP => M%H
         IF (CORRECTOR) HP => M%HS
         IF (RNODE/=SNODE) THEN
            PACK_REAL_SEND_PKG7: DO LL=1,M3%NIC_S
               SELECT CASE(M3%IOR_S(LL))
                  CASE(-1) ; M3%REAL_SEND_PKG7(3*LL-2) = M%FVX(M3%IIO_S(LL)-1,M3%JJO_S(LL)  ,M3%KKO_S(LL)  )
                             M3%REAL_SEND_PKG7(3*LL-1) =    HP(M3%IIO_S(LL)-1,M3%JJO_S(LL)  ,M3%KKO_S(LL)  )
                             M3%REAL_SEND_PKG7(3*LL  ) =    HP(M3%IIO_S(LL)  ,M3%JJO_S(LL)  ,M3%KKO_S(LL)  )
                  CASE( 1) ; M3%REAL_SEND_PKG7(3*LL-2) = M%FVX(M3%IIO_S(LL)  ,M3%JJO_S(LL)  ,M3%KKO_S(LL)  )
                             M3%REAL_SEND_PKG7(3*LL-1) =    HP(M3%IIO_S(LL)  ,M3%JJO_S(LL)  ,M3%KKO_S(LL)  )
                             M3%REAL_SEND_PKG7(3*LL  ) =    HP(M3%IIO_S(LL)+1,M3%JJO_S(LL)  ,M3%KKO_S(LL)  )
                  CASE(-2) ; M3%REAL_SEND_PKG7(3*LL-2) = M%FVY(M3%IIO_S(LL)  ,M3%JJO_S(LL)-1,M3%KKO_S(LL)  )
                             M3%REAL_SEND_PKG7(3*LL-1) =    HP(M3%IIO_S(LL)  ,M3%JJO_S(LL)-1,M3%KKO_S(LL)  )
                             M3%REAL_SEND_PKG7(3*LL  ) =    HP(M3%IIO_S(LL)  ,M3%JJO_S(LL)  ,M3%KKO_S(LL)  )
                  CASE( 2) ; M3%REAL_SEND_PKG7(3*LL-2) = M%FVY(M3%IIO_S(LL)  ,M3%JJO_S(LL)  ,M3%KKO_S(LL)  )
                             M3%REAL_SEND_PKG7(3*LL-1) =    HP(M3%IIO_S(LL)  ,M3%JJO_S(LL)  ,M3%KKO_S(LL)  )
                             M3%REAL_SEND_PKG7(3*LL  ) =    HP(M3%IIO_S(LL)  ,M3%JJO_S(LL)+1,M3%KKO_S(LL)  )
                  CASE(-3) ; M3%REAL_SEND_PKG7(3*LL-2) = M%FVZ(M3%IIO_S(LL)  ,M3%JJO_S(LL)  ,M3%KKO_S(LL)-1)
                             M3%REAL_SEND_PKG7(3*LL-1) =    HP(M3%IIO_S(LL)  ,M3%JJO_S(LL)  ,M3%KKO_S(LL)-1)
                             M3%REAL_SEND_PKG7(3*LL  ) =    HP(M3%IIO_S(LL)  ,M3%JJO_S(LL)  ,M3%KKO_S(LL)  )
                  CASE( 3) ; M3%REAL_SEND_PKG7(3*LL-2) = M%FVZ(M3%IIO_S(LL)  ,M3%JJO_S(LL)  ,M3%KKO_S(LL)  )
                             M3%REAL_SEND_PKG7(3*LL-1) =    HP(M3%IIO_S(LL)  ,M3%JJO_S(LL)  ,M3%KKO_S(LL)  )
                             M3%REAL_SEND_PKG7(3*LL  ) =    HP(M3%IIO_S(LL)  ,M3%JJO_S(LL)  ,M3%KKO_S(LL)+1)
               END SELECT
            ENDDO PACK_REAL_SEND_PKG7
         ELSE
            M2=>MESHES(NOM)%OMESH(NM)
            IF (PREDICTOR) HP2 => M2%H
            IF (CORRECTOR) HP2 => M2%HS
            M2%FVX(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) = M%FVX(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
            M2%FVY(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) = M%FVY(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
            M2%FVZ(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) = M%FVZ(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)         
            HP2(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)    = HP(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
         ENDIF
      ENDIF

      ! Send pressure information at the end of the PREDICTOR stage of the time step
 
      IF (CODE==3 .AND. M3%NIC_S>0) THEN
         IF (RNODE/=SNODE) THEN
            LL = 0
            DO KK=KMIN,KMAX
               DO JJ=JMIN,JMAX
                  DO II=IMIN,IMAX
                     M3%REAL_SEND_PKG2(LL+1) = M%HS(II,JJ,KK)
                     M3%REAL_SEND_PKG2(LL+2) = M%US(II,JJ,KK)
                     M3%REAL_SEND_PKG2(LL+3) = M%VS(II,JJ,KK)
                     M3%REAL_SEND_PKG2(LL+4) = M%WS(II,JJ,KK)
                     LL = LL+4
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            M2=>MESHES(NOM)%OMESH(NM)
            M2%HS(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) = M%HS(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
            M2%US(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) = M%US(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
            M2%VS(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) = M%VS(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
            M2%WS(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) = M%WS(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
         ENDIF
      ENDIF

      ! Exchange density and mass fraction following CORRECTOR update
 
      IF (CODE==4 .AND. M3%NIC_S>0) THEN
         IF (RNODE/=SNODE) THEN
            NQT2 = 2*(4+N_TOTAL_SCALARS)
            PACK_REAL_SEND_PKG3: DO LL=1,M3%NIC_S
               II1 = M3%IIO_S(LL) ; II2 = II1
               JJ1 = M3%JJO_S(LL) ; JJ2 = JJ1
               KK1 = M3%KKO_S(LL) ; KK2 = KK1
               SELECT CASE(M3%IOR_S(LL))
                  CASE(-1) ; II1=M3%IIO_S(LL)   ; II2=II1+1
                  CASE( 1) ; II1=M3%IIO_S(LL)-1 ; II2=II1+1
                  CASE(-2) ; JJ1=M3%JJO_S(LL)   ; JJ2=JJ1+1
                  CASE( 2) ; JJ1=M3%JJO_S(LL)-1 ; JJ2=JJ1+1
                  CASE(-3) ; KK1=M3%KKO_S(LL)   ; KK2=KK1+1
                  CASE( 3) ; KK1=M3%KKO_S(LL)-1 ; KK2=KK1+1
               END SELECT
               M3%REAL_SEND_PKG3(NQT2*(LL-1)+1) =  M%RHO(II1,JJ1,KK1)
               M3%REAL_SEND_PKG3(NQT2*(LL-1)+2) =  M%RHO(II2,JJ2,KK2)
               M3%REAL_SEND_PKG3(NQT2*(LL-1)+3) =   M%MU(II1,JJ1,KK1)
               M3%REAL_SEND_PKG3(NQT2*(LL-1)+4) =   M%MU(II2,JJ2,KK2)
               M3%REAL_SEND_PKG3(NQT2*(LL-1)+5) = M%KRES(II1,JJ1,KK1)
               M3%REAL_SEND_PKG3(NQT2*(LL-1)+6) = M%KRES(II2,JJ2,KK2)
               M3%REAL_SEND_PKG3(NQT2*(LL-1)+7) =   M%DS(II1,JJ1,KK1)
               M3%REAL_SEND_PKG3(NQT2*(LL-1)+8) =   M%DS(II2,JJ2,KK2)
               DO NN=1,N_TOTAL_SCALARS
                  M3%REAL_SEND_PKG3(NQT2*(LL-1)+8+2*NN-1) = M%ZZ(II1,JJ1,KK1,NN)
                  M3%REAL_SEND_PKG3(NQT2*(LL-1)+8+2*NN  ) = M%ZZ(II2,JJ2,KK2,NN)
               ENDDO
            ENDDO PACK_REAL_SEND_PKG3
         ELSE
            M2=>MESHES(NOM)%OMESH(NM)
            M2%RHO(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) = M%RHO(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
            M2%MU(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)  = M%MU(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
            M2%KRES(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)= M%KRES(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
            M2%DS(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)  = M%DS(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
            M2%ZZ(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX,1:N_TOTAL_SCALARS)= &
                                  M%ZZ(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX,1:N_TOTAL_SCALARS)
         ENDIF
      ENDIF

      ! Exchange BOUNDARY_TYPE following the CORRECTOR stage of the time step 

      IF (CODE==0 .OR. CODE==6) THEN
         IF (RNODE/=SNODE) THEN
            DO IW=1,M%N_EXTERNAL_WALL_CELLS
               M5%BOUNDARY_TYPE(IW) = M%WALL(IW)%BOUNDARY_TYPE
            ENDDO
         ELSE
            M2=>MESHES(NOM)%OMESH(NM)
            M2%BOUNDARY_TYPE(1:M%N_EXTERNAL_WALL_CELLS) = M5%BOUNDARY_TYPE(1:M%N_EXTERNAL_WALL_CELLS)
         ENDIF
      ENDIF

      ! Exchange BACK_WALL information

      IF (CODE==6) THEN
         IF (RNODE/=SNODE) THEN
            LL = 0
            DO II=1,M%OMESH(NOM)%N_BACK_WALL_CELLS
               IW = M%OMESH(NOM)%BACK_WALL_CELL_INDEX(II)
               M3%REAL_SEND_PKG6(LL+1) = M%WALL(IW)%ONE_D%QRADIN
               M3%REAL_SEND_PKG6(LL+2) = M%TMP(M%WALL(IW)%ONE_D%IIG,M%WALL(IW)%ONE_D%JJG,M%WALL(IW)%ONE_D%KKG)
               LL = LL+2
            ENDDO
         ELSE
            M2=>MESHES(NOM)%OMESH(NM)
            DO II=1,M2%N_EXT_BACK_WALL_CELLS
               IW = M2%EXT_BACK_WALL_CELL_INDEX(II)
               M2%BACK_WALL(II)%QRADIN = M%WALL(IW)%ONE_D%QRADIN
               M2%BACK_WALL(II)%TMP_GAS = M%TMP(M%WALL(IW)%ONE_D%IIG,M%WALL(IW)%ONE_D%JJG,M%WALL(IW)%ONE_D%KKG)
            ENDDO
         ENDIF
      ENDIF

      ! Exchange pressure and velocities following CORRECTOR stage of time step
 
      IF (CODE==6 .AND. M3%NIC_S>0) THEN
         IF (RNODE/=SNODE) THEN
            LL = 0
            DO KK=KMIN,KMAX
               DO JJ=JMIN,JMAX
                  DO II=IMIN,IMAX
                     M3%REAL_SEND_PKG4(LL+1) = M%H(II,JJ,KK)
                     M3%REAL_SEND_PKG4(LL+2) = M%U(II,JJ,KK)
                     M3%REAL_SEND_PKG4(LL+3) = M%V(II,JJ,KK)
                     M3%REAL_SEND_PKG4(LL+4) = M%W(II,JJ,KK)
                     LL = LL+4
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            M2=>MESHES(NOM)%OMESH(NM)
            M2%H(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)    = M%H(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
            M2%U(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)    = M%U(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
            M2%V(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)    = M%V(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
            M2%W(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)    = M%W(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
         ENDIF
      ENDIF

      ! Send out radiation info

      SEND_RADIATION: IF (CODE==2 .AND. M3%NIC_S>0) THEN
         IF (RNODE/=SNODE) THEN
            IF (ICYC>1) ANG_INC_COUNTER = M%ANGLE_INC_COUNTER
            LLL = 0
            PACK_REAL_SEND_PKG5: DO LL=1,M3%NIC_S
               IOR = M3%IOR_S(LL)
               DO NN2=1,NUMBER_SPECTRAL_BANDS
                  DO NN1=NUMBER_RADIATION_ANGLES-ANG_INC_COUNTER+1,1,-ANGLE_INCREMENT
                     IF (DLN(IOR,NN1)<=0._EB) CYCLE
                     LLL = LLL + 1
                     M3%REAL_SEND_PKG5(LLL) = M3%IL_S(LL,NN1,NN2)
                  ENDDO
               ENDDO
            ENDDO PACK_REAL_SEND_PKG5
         ELSE
            M2=>MESHES(NOM)%OMESH(NM)
            M2%IL_R = M3%IL_S
         ENDIF
      ENDIF SEND_RADIATION
 
      ! Get Number of PARTICLE Orphans (PARTICLEs that have left other meshes and are waiting to be picked up)

      IF (CODE==7 .AND. OMESH_PARTICLES) THEN
         IF (RNODE==SNODE) THEN
            M2=>MESHES(NOM)%OMESH(NM)
            M2%N_PART_ADOPT = M3%N_PART_ORPHANS
         ENDIF
      ENDIF

      ! Sending/Receiving PARTICLE Buffer Arrays

      IF_SEND_PARTICLES: IF (CODE==6 .AND. OMESH_PARTICLES) THEN 
 
         NODE_CHECK_PARTICLE: IF (SNODE/=RNODE) THEN
            
            DO IPC=1,N_LAGRANGIAN_CLASSES

               IF (M3%N_PART_ORPHANS(IPC)==0) CYCLE

               LPC => LAGRANGIAN_PARTICLE_CLASS(IPC)
               IBC = LPC%SURF_INDEX

               N_STORAGE_SLOTS = M3%ORPHAN_PARTICLE_STORAGE(IPC)%N_STORAGE_SLOTS
               N_REQ=MIN(N_REQ+1,SIZE(REQ))
               CALL MPI_ISEND(M3%ORPHAN_PARTICLE_STORAGE(IPC)%REALS(1,1),LPC%N_STORAGE_REALS*N_STORAGE_SLOTS,MPI_DOUBLE_PRECISION, &
                              SNODE,NM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
               N_REQ=MIN(N_REQ+1,SIZE(REQ))
               CALL MPI_ISEND(M3%ORPHAN_PARTICLE_STORAGE(IPC)%INTEGERS(1,1),LPC%N_STORAGE_INTEGERS*N_STORAGE_SLOTS,MPI_INTEGER, &
                              SNODE,NM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
               N_REQ=MIN(N_REQ+1,SIZE(REQ))
               CALL MPI_ISEND(M3%ORPHAN_PARTICLE_STORAGE(IPC)%LOGICALS(1,1),LPC%N_STORAGE_LOGICALS*N_STORAGE_SLOTS,MPI_LOGICAL, &
                              SNODE,NM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
            ENDDO

         ELSE NODE_CHECK_PARTICLE

            M2 => MESHES(NOM)%OMESH(NM)

            DO IPC=1,N_LAGRANGIAN_CLASSES
               LPC => LAGRANGIAN_PARTICLE_CLASS(IPC)
               M2%ADOPT_PARTICLE_STORAGE(IPC)%REALS    = M3%ORPHAN_PARTICLE_STORAGE(IPC)%REALS
               M2%ADOPT_PARTICLE_STORAGE(IPC)%INTEGERS = M3%ORPHAN_PARTICLE_STORAGE(IPC)%INTEGERS
               M2%ADOPT_PARTICLE_STORAGE(IPC)%LOGICALS = M3%ORPHAN_PARTICLE_STORAGE(IPC)%LOGICALS
            ENDDO

         ENDIF NODE_CHECK_PARTICLE

      ENDIF IF_SEND_PARTICLES
 
   ENDDO RECEIVING_MESH_LOOP

ENDDO SENDING_MESH_LOOP

! Send information needed by EVACuation routine

DO NM=1,NMESHES
   IF (CODE==6 .AND. EXCHANGE_EVACUATION.AND.MYID/=MAX(0,EVAC_PROCESS).AND.PROCESS(NM)==MYID .AND..NOT.EVACUATION_ONLY(NM)) THEN
      M => MESHES(NM)
      TAG_EVAC = NM*(MAX(0,EVAC_PROCESS)+1)*CODE*10
      IWW = (M%IBAR+2)*(M%JBAR+2)*(M%KBAR+2)
      N_REQ=MIN(N_REQ+1,SIZE(REQ))
      CALL MPI_ISEND(M%ZZ(0,0,0,1),IWW*N_TRACKED_SPECIES,MPI_DOUBLE_PRECISION,MAX(0,EVAC_PROCESS),&
           TAG_EVAC,MPI_COMM_WORLD,REQ(N_REQ),IERR)
      N_REQ=MIN(N_REQ+1,SIZE(REQ))
      CALL MPI_ISEND(M%RHO(0,0,0),IWW,MPI_DOUBLE_PRECISION,MAX(0,EVAC_PROCESS),TAG_EVAC,MPI_COMM_WORLD,REQ(N_REQ),IERR)
      N_REQ=MIN(N_REQ+1,SIZE(REQ))
      CALL MPI_ISEND(M%RSUM(0,0,0),IWW,MPI_DOUBLE_PRECISION,MAX(0,EVAC_PROCESS),TAG_EVAC,MPI_COMM_WORLD,REQ(N_REQ),IERR)
      N_REQ=MIN(N_REQ+1,SIZE(REQ))
      CALL MPI_ISEND(M%TMP(0,0,0),IWW,MPI_DOUBLE_PRECISION,MAX(0,EVAC_PROCESS),TAG_EVAC,MPI_COMM_WORLD,REQ(N_REQ),IERR)
      N_REQ=MIN(N_REQ+1,SIZE(REQ))
      CALL MPI_ISEND(M%UII(0,0,0),IWW,MPI_DOUBLE_PRECISION,MAX(0,EVAC_PROCESS),TAG_EVAC,MPI_COMM_WORLD,REQ(N_REQ),IERR)
      N_REQ=MIN(N_REQ+1,SIZE(REQ))
      CALL MPI_ISEND(M%CELL_INDEX(0,0,0),IWW,MPI_INTEGER,MAX(0,EVAC_PROCESS),TAG_EVAC,MPI_COMM_WORLD,REQ(N_REQ),IERR)
      IWW = MAXVAL(M%CELL_INDEX)
      N_REQ=MIN(N_REQ+1,SIZE(REQ))
      CALL MPI_ISEND(M%SOLID(0),IWW,MPI_LOGICAL,MAX(0,EVAC_PROCESS),TAG_EVAC,MPI_COMM_WORLD,REQ(N_REQ),IERR)
   ENDIF
ENDDO


! Halt communications until all processes are ready to receive the data.

IF (N_MPI_PROCESSES>1 .AND. CODE/=1 .AND. CODE/=3 .AND. CODE/=4 .AND. CODE/=5 .AND. N_REQ>0) THEN
   CALL TIMEOUT('REQ',N_REQ,REQ(1:N_REQ))
ENDIF

IF (N_MPI_PROCESSES>1 .AND. CODE==1 .AND. N_REQ1>0) THEN
   CALL MPI_STARTALL(N_REQ1,REQ1(1:N_REQ1),IERR)
   CALL TIMEOUT('REQ1',N_REQ1,REQ1(1:N_REQ1))
ENDIF

IF (N_MPI_PROCESSES>1 .AND. CODE==7 .AND. OMESH_PARTICLES .AND. N_REQ2>0) THEN
   CALL MPI_STARTALL(N_REQ2,REQ2(1:N_REQ2),IERR)
   CALL TIMEOUT('REQ2',N_REQ2,REQ2(1:N_REQ2))
ENDIF

IF (N_MPI_PROCESSES>1 .AND. CODE==3 .AND. N_REQ3>0) THEN
   CALL MPI_STARTALL(N_REQ3,REQ3(1:N_REQ3),IERR)
   CALL TIMEOUT('REQ3',N_REQ3,REQ3(1:N_REQ3))
ENDIF

IF (N_MPI_PROCESSES>1 .AND. CODE==4 .AND. N_REQ4>0) THEN
   CALL MPI_STARTALL(N_REQ4,REQ4(1:N_REQ4),IERR)
   CALL TIMEOUT('REQ4',N_REQ4,REQ4(1:N_REQ4))
ENDIF

IF (N_MPI_PROCESSES>1 .AND. CODE==5 .AND. N_REQ5>0) THEN
   CALL MPI_STARTALL(N_REQ5,REQ5(1:N_REQ5),IERR)
   CALL TIMEOUT('REQ5',N_REQ5,REQ5(1:N_REQ5))
ENDIF

IF (N_MPI_PROCESSES>1 .AND. CODE==6 .AND. N_REQ6>0) THEN
   CALL MPI_STARTALL(N_REQ6,REQ6(1:N_REQ6),IERR)
   CALL TIMEOUT('REQ6',N_REQ6,REQ6(1:N_REQ6))
ENDIF

IF (N_MPI_PROCESSES>1 .AND. CODE==6 .AND. N_REQ7>0) THEN
   CALL MPI_STARTALL(N_REQ7,REQ7(1:N_REQ7),IERR)
   CALL TIMEOUT('REQ7',N_REQ7,REQ7(1:N_REQ7))
ENDIF

IF (N_MPI_PROCESSES>1 .AND. (CODE==0 .OR. CODE==6) .AND. N_REQ8>0) THEN
   CALL MPI_STARTALL(N_REQ8,REQ8(1:N_REQ8),IERR)
   CALL TIMEOUT('REQ8',N_REQ8,REQ8(1:N_REQ8))
ENDIF

IF (N_MPI_PROCESSES>1 .AND. CODE==2 .AND. N_REQ9>0) THEN
   CALL MPI_STARTALL(N_REQ9,REQ9(1:N_REQ9),IERR)
   CALL TIMEOUT('REQ9',N_REQ9,REQ9(1:N_REQ9))
ENDIF


! Receive the information sent above into the appropriate arrays.

SEND_MESH_LOOP: DO NOM=1,NMESHES

IF (EVACUATION_ONLY(NOM)) CYCLE SEND_MESH_LOOP 

SNODE = PROCESS(NOM)
IF (SNODE/=MYID) CYCLE SEND_MESH_LOOP
 
   RECV_MESH_LOOP: DO NM=1,NMESHES
 
      IF (EVACUATION_ONLY(NM)) CYCLE RECV_MESH_LOOP
      IF (CODE>0 .AND. (.NOT.ACTIVE_MESH(NM).OR..NOT.ACTIVE_MESH(NOM))) CYCLE RECV_MESH_LOOP

      RNODE = PROCESS(NM)

      M =>MESHES(NM)
      M2=>MESHES(NOM)%OMESH(NM)
      M4=>MESHES(NOM)

      IF (M2%NIC_S==0 .AND. M2%NIC_R==0) CYCLE RECV_MESH_LOOP
 
      IMIN = M2%I_MIN_R
      IMAX = M2%I_MAX_R
      JMIN = M2%J_MIN_R
      JMAX = M2%J_MAX_R
      KMIN = M2%K_MIN_R
      KMAX = M2%K_MAX_R
     
      ! Unpack densities and species mass fractions following PREDICTOR exchange

      IF (CODE==1 .AND. M2%NIC_R>0 .AND. RNODE/=SNODE) THEN
            NQT2 = 2*(4+N_TOTAL_SCALARS)
            UNPACK_REAL_RECV_PKG1: DO LL=1,M2%NIC_R
               II1 = M2%IIO_R(LL) ; II2 = II1
               JJ1 = M2%JJO_R(LL) ; JJ2 = JJ1
               KK1 = M2%KKO_R(LL) ; KK2 = KK1
               SELECT CASE(M2%IOR_R(LL))
                  CASE(-1) ; II1=M2%IIO_R(LL)   ; II2=II1+1
                  CASE( 1) ; II1=M2%IIO_R(LL)-1 ; II2=II1+1
                  CASE(-2) ; JJ1=M2%JJO_R(LL)   ; JJ2=JJ1+1
                  CASE( 2) ; JJ1=M2%JJO_R(LL)-1 ; JJ2=JJ1+1
                  CASE(-3) ; KK1=M2%KKO_R(LL)   ; KK2=KK1+1
                  CASE( 3) ; KK1=M2%KKO_R(LL)-1 ; KK2=KK1+1
               END SELECT
               M2%RHOS(II1,JJ1,KK1) = M2%REAL_RECV_PKG1(NQT2*(LL-1)+1)
               M2%RHOS(II2,JJ2,KK2) = M2%REAL_RECV_PKG1(NQT2*(LL-1)+2)
                 M2%MU(II1,JJ1,KK1) = M2%REAL_RECV_PKG1(NQT2*(LL-1)+3)
                 M2%MU(II2,JJ2,KK2) = M2%REAL_RECV_PKG1(NQT2*(LL-1)+4)
               M2%KRES(II1,JJ1,KK1) = M2%REAL_RECV_PKG1(NQT2*(LL-1)+5)
               M2%KRES(II2,JJ2,KK2) = M2%REAL_RECV_PKG1(NQT2*(LL-1)+6)
                  M2%D(II1,JJ1,KK1) = M2%REAL_RECV_PKG1(NQT2*(LL-1)+7)
                  M2%D(II2,JJ2,KK2) = M2%REAL_RECV_PKG1(NQT2*(LL-1)+8)
               DO NN=1,N_TOTAL_SCALARS
                  M2%ZZS(II1,JJ1,KK1,NN) = M2%REAL_RECV_PKG1(NQT2*(LL-1)+8+2*NN-1)
                  M2%ZZS(II2,JJ2,KK2,NN) = M2%REAL_RECV_PKG1(NQT2*(LL-1)+8+2*NN  )
               ENDDO
            ENDDO UNPACK_REAL_RECV_PKG1
      ENDIF 
   
      ! Unpack densities and species mass fractions following PREDICTOR exchange
   
      IF (CODE==5 .AND. M2%NIC_R>0 .AND. RNODE/=SNODE) THEN
         IF (PREDICTOR) HP => M2%H
         IF (CORRECTOR) HP => M2%HS
         UNPACK_REAL_RECV_PKG7: DO LL=1,M2%NIC_R
            SELECT CASE(M2%IOR_R(LL))
               CASE(-1) ; M2%FVX(M2%IIO_R(LL)-1,M2%JJO_R(LL)  ,M2%KKO_R(LL)  ) = M2%REAL_RECV_PKG7(3*LL-2)
                              HP(M2%IIO_R(LL)-1,M2%JJO_R(LL)  ,M2%KKO_R(LL)  ) = M2%REAL_RECV_PKG7(3*LL-1)
                              HP(M2%IIO_R(LL)  ,M2%JJO_R(LL)  ,M2%KKO_R(LL)  ) = M2%REAL_RECV_PKG7(3*LL  )
               CASE( 1) ; M2%FVX(M2%IIO_R(LL)  ,M2%JJO_R(LL)  ,M2%KKO_R(LL)  ) = M2%REAL_RECV_PKG7(3*LL-2)
                              HP(M2%IIO_R(LL)  ,M2%JJO_R(LL)  ,M2%KKO_R(LL)  ) = M2%REAL_RECV_PKG7(3*LL-1)
                              HP(M2%IIO_R(LL)+1,M2%JJO_R(LL)  ,M2%KKO_R(LL)  ) = M2%REAL_RECV_PKG7(3*LL  )
               CASE(-2) ; M2%FVY(M2%IIO_R(LL)  ,M2%JJO_R(LL)-1,M2%KKO_R(LL)  ) = M2%REAL_RECV_PKG7(3*LL-2)
                              HP(M2%IIO_R(LL)  ,M2%JJO_R(LL)-1,M2%KKO_R(LL)  ) = M2%REAL_RECV_PKG7(3*LL-1)
                              HP(M2%IIO_R(LL)  ,M2%JJO_R(LL)  ,M2%KKO_R(LL)  ) = M2%REAL_RECV_PKG7(3*LL  )
               CASE( 2) ; M2%FVY(M2%IIO_R(LL)  ,M2%JJO_R(LL)  ,M2%KKO_R(LL)  ) = M2%REAL_RECV_PKG7(3*LL-2)
                              HP(M2%IIO_R(LL)  ,M2%JJO_R(LL)  ,M2%KKO_R(LL)  ) = M2%REAL_RECV_PKG7(3*LL-1)
                              HP(M2%IIO_R(LL)  ,M2%JJO_R(LL)+1,M2%KKO_R(LL)  ) = M2%REAL_RECV_PKG7(3*LL  )
               CASE(-3) ; M2%FVZ(M2%IIO_R(LL)  ,M2%JJO_R(LL)  ,M2%KKO_R(LL)-1) = M2%REAL_RECV_PKG7(3*LL-2)
                              HP(M2%IIO_R(LL)  ,M2%JJO_R(LL)  ,M2%KKO_R(LL)-1) = M2%REAL_RECV_PKG7(3*LL-1)
                              HP(M2%IIO_R(LL)  ,M2%JJO_R(LL)  ,M2%KKO_R(LL)  ) = M2%REAL_RECV_PKG7(3*LL  )
               CASE( 3) ; M2%FVZ(M2%IIO_R(LL)  ,M2%JJO_R(LL)  ,M2%KKO_R(LL)  ) = M2%REAL_RECV_PKG7(3*LL-2)
                              HP(M2%IIO_R(LL)  ,M2%JJO_R(LL)  ,M2%KKO_R(LL)  ) = M2%REAL_RECV_PKG7(3*LL-1)
                              HP(M2%IIO_R(LL)  ,M2%JJO_R(LL)  ,M2%KKO_R(LL)+1) = M2%REAL_RECV_PKG7(3*LL  )
            END SELECT
         ENDDO UNPACK_REAL_RECV_PKG7
      ENDIF
   
      ! Unpack pressure following PREDICTOR stage of time step
    
      IF (CODE==3 .AND. M2%NIC_R>0 .AND. RNODE/=SNODE) THEN
         LL = 0
         DO KK=KMIN,KMAX
            DO JJ=JMIN,JMAX
               DO II=IMIN,IMAX
                  M2%HS(II,JJ,KK)   = M2%REAL_RECV_PKG2(LL+1)
                  M2%US(II,JJ,KK)   = M2%REAL_RECV_PKG2(LL+2)
                  M2%VS(II,JJ,KK)   = M2%REAL_RECV_PKG2(LL+3)
                  M2%WS(II,JJ,KK)   = M2%REAL_RECV_PKG2(LL+4)
                  LL = LL+4
               ENDDO
            ENDDO
         ENDDO
      ENDIF 
   
      ! Unpack density and species mass fractions following CORRECTOR update
   
      IF (CODE==4 .AND. M2%NIC_R>0 .AND. RNODE/=SNODE) THEN
         NQT2 = 2*(4+N_TOTAL_SCALARS)
         UNPACK_REAL_RECV_PKG3: DO LL=1,M2%NIC_R
            II1 = M2%IIO_R(LL) ; II2 = II1
            JJ1 = M2%JJO_R(LL) ; JJ2 = JJ1
            KK1 = M2%KKO_R(LL) ; KK2 = KK1
            SELECT CASE(M2%IOR_R(LL))
               CASE(-1) ; II1=M2%IIO_R(LL)   ; II2=II1+1
               CASE( 1) ; II1=M2%IIO_R(LL)-1 ; II2=II1+1
               CASE(-2) ; JJ1=M2%JJO_R(LL)   ; JJ2=JJ1+1
               CASE( 2) ; JJ1=M2%JJO_R(LL)-1 ; JJ2=JJ1+1
               CASE(-3) ; KK1=M2%KKO_R(LL)   ; KK2=KK1+1
               CASE( 3) ; KK1=M2%KKO_R(LL)-1 ; KK2=KK1+1
            END SELECT
             M2%RHO(II1,JJ1,KK1) = M2%REAL_RECV_PKG3(NQT2*(LL-1)+1)
             M2%RHO(II2,JJ2,KK2) = M2%REAL_RECV_PKG3(NQT2*(LL-1)+2)
              M2%MU(II1,JJ1,KK1) = M2%REAL_RECV_PKG3(NQT2*(LL-1)+3)
              M2%MU(II2,JJ2,KK2) = M2%REAL_RECV_PKG3(NQT2*(LL-1)+4)
            M2%KRES(II1,JJ1,KK1) = M2%REAL_RECV_PKG3(NQT2*(LL-1)+5)
            M2%KRES(II2,JJ2,KK2) = M2%REAL_RECV_PKG3(NQT2*(LL-1)+6)
              M2%DS(II1,JJ1,KK1) = M2%REAL_RECV_PKG3(NQT2*(LL-1)+7)
              M2%DS(II2,JJ2,KK2) = M2%REAL_RECV_PKG3(NQT2*(LL-1)+8)
            DO NN=1,N_TOTAL_SCALARS
               M2%ZZ(II1,JJ1,KK1,NN) = M2%REAL_RECV_PKG3(NQT2*(LL-1)+8+2*NN-1)
               M2%ZZ(II2,JJ2,KK2,NN) = M2%REAL_RECV_PKG3(NQT2*(LL-1)+8+2*NN)
            ENDDO
         ENDDO UNPACK_REAL_RECV_PKG3
      ENDIF
   
      ! Unpack pressure and velocities at the end of the CORRECTOR stage of the time step
    
      IF (CODE==6 .AND. M2%NIC_R>0 .AND. RNODE/=SNODE) THEN
         LL = 0
         DO KK=KMIN,KMAX
            DO JJ=JMIN,JMAX
               DO II=IMIN,IMAX
                  M2%H(II,JJ,KK)    = M2%REAL_RECV_PKG4(LL+1)
                  M2%U(II,JJ,KK)    = M2%REAL_RECV_PKG4(LL+2)
                  M2%V(II,JJ,KK)    = M2%REAL_RECV_PKG4(LL+3)
                  M2%W(II,JJ,KK)    = M2%REAL_RECV_PKG4(LL+4)
                  LL = LL+4
               ENDDO
            ENDDO
         ENDDO
      ENDIF
   
      ! Unpack radiation information at the end of the CORRECTOR stage of the time step
   
      RECEIVE_RADIATION: IF (CODE==2 .AND. M2%NIC_R>0 .AND. RNODE/=SNODE) THEN
         IF (ICYC>1) ANG_INC_COUNTER = M4%ANGLE_INC_COUNTER
         LLL = 0
         UNPACK_REAL_RECV_PKG5: DO LL=1,M2%NIC_R
            IOR = M2%IOR_R(LL)
            DO NN2=1,NUMBER_SPECTRAL_BANDS
               DO NN1=NUMBER_RADIATION_ANGLES-ANG_INC_COUNTER+1,1,-ANGLE_INCREMENT
                  IF (DLN(IOR,NN1)<=0._EB) CYCLE
                  LLL = LLL + 1
                  M2%IL_R(LL,NN1,NN2) = M2%REAL_RECV_PKG5(LLL)
               ENDDO
            ENDDO
         ENDDO UNPACK_REAL_RECV_PKG5
      ENDIF RECEIVE_RADIATION

      ! Unpack back wall information at the end of the CORRECTOR stage of the time step

      RECEIVE_BACK_WALL: IF (CODE==6 .AND. SNODE/=RNODE) THEN
         LL = 0
         DO II=1,M2%N_EXT_BACK_WALL_CELLS
            M2%BACK_WALL(II)%QRADIN  = M2%REAL_RECV_PKG6(LL+1)
            M2%BACK_WALL(II)%TMP_GAS = M2%REAL_RECV_PKG6(LL+2)
            LL = LL+2
         ENDDO
      ENDIF RECEIVE_BACK_WALL

      ! Sending/Receiving PARTICLE Buffer Arrays
    
      IF (CODE==7 .AND. OMESH_PARTICLES) THEN 
         DO IPC=1,N_LAGRANGIAN_CLASSES
            IF (M2%N_PART_ADOPT(IPC)>M2%ADOPT_PARTICLE_STORAGE(IPC)%N_STORAGE_SLOTS) THEN
               N_NEW_STORAGE_SLOTS = M2%N_PART_ADOPT(IPC)-M2%ADOPT_PARTICLE_STORAGE(IPC)%N_STORAGE_SLOTS
               CALL REALLOCATE_STORAGE_ARRAYS(NOM,3,IPC,N_NEW_STORAGE_SLOTS,NM)
            ENDIF
         ENDDO
      ENDIF

      IF_RECEIVE_PARTICLES: IF (CODE==6 .AND. OMESH_PARTICLES) THEN 
   
         DO IPC=1,N_LAGRANGIAN_CLASSES
            IF (M2%N_PART_ADOPT(IPC)==0) CYCLE
            CNT = 0
            DO N=M4%NLP+1,M4%NLP+M2%N_PART_ADOPT(IPC)
               CNT = CNT + 1
               IBC = LAGRANGIAN_PARTICLE_CLASS(IPC)%SURF_INDEX
               CALL ALLOCATE_STORAGE(NOM,IBC,LPC_INDEX=IPC,LP_INDEX=N,TAG=-1)
               LP=>M4%LAGRANGIAN_PARTICLE(N)
               STORAGE_INDEX_SAVE = LP%STORAGE_INDEX
 
               M4%PARTICLE_STORAGE(IPC)%REALS(:,LP%STORAGE_INDEX)    = M2%ADOPT_PARTICLE_STORAGE(IPC)%REALS(:,CNT)
               M4%PARTICLE_STORAGE(IPC)%INTEGERS(:,LP%STORAGE_INDEX) = M2%ADOPT_PARTICLE_STORAGE(IPC)%INTEGERS(:,CNT)

               LP%ARRAY_INDEX = N
               LP%STORAGE_INDEX = STORAGE_INDEX_SAVE

               M4%PARTICLE_STORAGE(IPC)%LOGICALS(:,LP%STORAGE_INDEX) = M2%ADOPT_PARTICLE_STORAGE(IPC)%LOGICALS(:,CNT)
            ENDDO
            M4%NLP = M4%NLP + M2%N_PART_ADOPT(IPC)
         ENDDO
   
      ENDIF IF_RECEIVE_PARTICLES
   
   ENDDO RECV_MESH_LOOP
   
ENDDO SEND_MESH_LOOP

! Ensure that all processes wait until all exchanges are complete

CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
    
TUSED(11,:)=TUSED(11,:) + SECOND() - TNOW
END SUBROUTINE MESH_EXCHANGE


SUBROUTINE TIMEOUT(RNAME,NR,RR)

REAL(EB) :: START_TIME,WAIT_TIME
INTEGER :: NR
INTEGER, DIMENSION(:) :: RR
LOGICAL :: FLAG
CHARACTER(*) :: RNAME

START_TIME = MPI_WTIME()
FLAG = .FALSE.
DO WHILE(.NOT.FLAG)
   CALL MPI_TESTALL(NR,RR(1:NR),FLAG,MPI_STATUSES_IGNORE,IERR)
   WAIT_TIME = MPI_WTIME() - START_TIME
   IF (WAIT_TIME>MPI_TIMEOUT) THEN
      WRITE(LU_ERR,'(A,A,I3,A,A)') TRIM(RNAME),' timed out for MPI process ',MYID,' running on ',PNAME(1:PNAMELEN)
      CALL MPI_ABORT(MPI_COMM_WORLD,0,IERR)
   ENDIF
ENDDO

END SUBROUTINE TIMEOUT


SUBROUTINE HANDLE_MPI_ERROR(ERROR_MESSAGE,ERROR_CODE)

INTEGER :: ERROR_CODE
CHARACTER(*) :: ERROR_MESSAGE

WRITE(LU_ERR,'(A,A,I3,A,I4)') TRIM(ERROR_MESSAGE),', ERROR_CODE=',ERROR_CODE,', MPI_PROCESS=',MYID
CALL MPI_ABORT(MPI_COMM_WORLD,0,IERR)

END SUBROUTINE HANDLE_MPI_ERROR


SUBROUTINE WRITE_STRINGS
 
! Write character strings out to the .smv file
 
INTEGER :: N,NOM,N_STRINGS_DUM
CHARACTER(80), ALLOCATABLE, DIMENSION(:) :: STRING_DUM
 
! All meshes send their STRINGs to node 0
 
DO NM=1,NMESHES
   IF (PROCESS(NM)==MYID .AND. MYID>0) THEN
      CALL MPI_SEND(MESHES(NM)%N_STRINGS,1,MPI_INTEGER,0,1,MPI_COMM_WORLD,IERR)
      IF (MESHES(NM)%N_STRINGS>0) CALL MPI_SEND(MESHES(NM)%STRING(1),MESHES(NM)%N_STRINGS*80,MPI_CHARACTER,0,NM, &
                                                MPI_COMM_WORLD,IERR)
   ENDIF
ENDDO
 
! Node 0 receives the STRINGs and writes them to the .smv file
 
IF (MYID==0) THEN
   DO N=1,MESHES(1)%N_STRINGS
      WRITE(LU_SMV,'(A)') TRIM(MESHES(1)%STRING(N))
   ENDDO
   OTHER_MESH_LOOP: DO NOM=2,NMESHES 
      IF (PROCESS(NOM)>0) THEN
         CALL MPI_RECV(N_STRINGS_DUM,1,MPI_INTEGER,PROCESS(NOM),1,MPI_COMM_WORLD,STATUS,IERR)
         IF (N_STRINGS_DUM>0) THEN
            ALLOCATE(STRING_DUM(N_STRINGS_DUM))
            CALL MPI_RECV(STRING_DUM(1),N_STRINGS_DUM*80,MPI_CHARACTER,PROCESS(NOM),NOM,MPI_COMM_WORLD,STATUS,IERR)
         ENDIF
      ELSE
         N_STRINGS_DUM = MESHES(NOM)%N_STRINGS
         IF (N_STRINGS_DUM>0) THEN
            ALLOCATE(STRING_DUM(N_STRINGS_DUM))
            STRING_DUM(1:N_STRINGS_DUM) = MESHES(NOM)%STRING(1:N_STRINGS_DUM)
         ENDIF
      ENDIF
      DO N=1,N_STRINGS_DUM
         WRITE(LU_SMV,'(A)') TRIM(STRING_DUM(N))
      ENDDO
      IF (ALLOCATED(STRING_DUM)) DEALLOCATE(STRING_DUM)
   ENDDO OTHER_MESH_LOOP
ENDIF
 
! All STRING arrays are zeroed out
 
DO NM=1,NMESHES
   IF (PROCESS(NM)==MYID) MESHES(NM)%N_STRINGS = 0
ENDDO
 
END SUBROUTINE WRITE_STRINGS


SUBROUTINE EXCHANGE_DIAGNOSTICS
 
INTEGER  :: NOM,NECYC,CNT, N_TIMERS_TMP,DISP
REAL(EB) :: T_SUM, TNOW
 
TNOW = SECOND()

DO NM=1,NMESHES
   IF (PROCESS(NM)/=MYID) CYCLE
   T_SUM = 0._EB
   N_TIMERS_TMP = N_TIMERS_FDS
   IF (EVACUATION_GRID(NM)) N_TIMERS_TMP = N_TIMERS_EVAC - 3
   SUM_LOOP: DO I=2,N_TIMERS_TMP
      T_SUM = T_SUM + TUSED(I,NM)
   ENDDO SUM_LOOP
   NECYC          = MAX(1,NTCYC(NM)-NCYC(NM))
   T_PER_STEP(NM) = (T_SUM-T_ACCUM(NM))/REAL(NECYC,EB)
   T_ACCUM(NM)    = T_SUM
   NCYC(NM)       = NTCYC(NM)
ENDDO
 
DISP = DISPLS(MYID)+1
CNT  = COUNTS(MYID)
IF (N_MPI_PROCESSES>1) THEN
   REAL_BUFFER_1 = T
   CALL MPI_GATHERV(REAL_BUFFER_1(DISP),CNT,MPI_DOUBLE_PRECISION,T,COUNTS,DISPLS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
   IF (IERR/=MPI_SUCCESS) CALL HANDLE_MPI_ERROR('Error in MPI_GATHERV for T',IERR)
   REAL_BUFFER_1 = T_ACCUM
   CALL MPI_GATHERV(REAL_BUFFER_1(DISP),CNT,MPI_DOUBLE_PRECISION,T_ACCUM,COUNTS,DISPLS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
   IF (IERR/=MPI_SUCCESS) CALL HANDLE_MPI_ERROR('Error in MPI_GATHERV for T_ACCUM',IERR)
   REAL_BUFFER_1 = T_PER_STEP
   CALL MPI_GATHERV(REAL_BUFFER_1(DISP),CNT,MPI_DOUBLE_PRECISION, &
                    T_PER_STEP,COUNTS,DISPLS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
   IF (IERR/=MPI_SUCCESS) CALL HANDLE_MPI_ERROR('Error in MPI_GATHERV for T_PER_STEP',IERR)
   INTEGER_BUFFER_1 = NTCYC
   CALL MPI_GATHERV(INTEGER_BUFFER_1(DISP),CNT,MPI_INTEGER,NTCYC,COUNTS,DISPLS,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
   IF (IERR/=MPI_SUCCESS) CALL HANDLE_MPI_ERROR('Error in MPI_GATHERV for NTCYC',IERR)
   REAL_BUFFER_11 = Q_DOT
   CALL MPI_GATHERV(REAL_BUFFER_11(1,DISP),CNT*N_Q_DOT,MPI_DOUBLE_PRECISION,Q_DOT,COUNTS_Q_DOT,DISPLS_Q_DOT, &
                    MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
   IF (IERR/=MPI_SUCCESS) CALL HANDLE_MPI_ERROR('Error in MPI_GATHERV for Q_DOT',IERR)
   REAL_BUFFER_12 = M_DOT
   CALL MPI_GATHERV(REAL_BUFFER_12(1,DISP),CNT*N_M_DOT,MPI_DOUBLE_PRECISION,M_DOT,COUNTS_M_DOT,DISPLS_M_DOT, &
                    MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
   IF (IERR/=MPI_SUCCESS) CALL HANDLE_MPI_ERROR('Error in MPI_GATHERV for M_DOT',IERR)
ENDIF
 
! All nodes greater than 0 send various values to node 0

DO NM=1,NMESHES
   IF (PROCESS(NM)/=MYID .OR. MYID==0) CYCLE
   CALL MPI_SEND(MESHES(NM)%DT,                      1,MPI_DOUBLE_PRECISION, 0,1,MPI_COMM_WORLD,IERR)
   CALL MPI_SEND(MESHES(NM)%CFL,                     1,MPI_DOUBLE_PRECISION, 0,1,MPI_COMM_WORLD,IERR)
   CALL MPI_SEND(MESHES(NM)%DIVMX,                   1,MPI_DOUBLE_PRECISION, 0,1,MPI_COMM_WORLD,IERR)
   CALL MPI_SEND(MESHES(NM)%DIVMN,                   1,MPI_DOUBLE_PRECISION, 0,1,MPI_COMM_WORLD,IERR)
   CALL MPI_SEND(MESHES(NM)%RESMAX,                  1,MPI_DOUBLE_PRECISION, 0,1,MPI_COMM_WORLD,IERR)
   CALL MPI_SEND(MESHES(NM)%POIS_PTB,                1,MPI_DOUBLE_PRECISION, 0,1,MPI_COMM_WORLD,IERR)
   CALL MPI_SEND(MESHES(NM)%POIS_ERR,                1,MPI_DOUBLE_PRECISION, 0,1,MPI_COMM_WORLD,IERR)
   CALL MPI_SEND(MESHES(NM)%VN,                      1,MPI_DOUBLE_PRECISION, 0,1,MPI_COMM_WORLD,IERR)
   CALL MPI_SEND(MESHES(NM)%ICFL,                    1,MPI_INTEGER,          0,1,MPI_COMM_WORLD,IERR)
   CALL MPI_SEND(MESHES(NM)%JCFL,                    1,MPI_INTEGER,          0,1,MPI_COMM_WORLD,IERR)
   CALL MPI_SEND(MESHES(NM)%KCFL,                    1,MPI_INTEGER,          0,1,MPI_COMM_WORLD,IERR)
   CALL MPI_SEND(MESHES(NM)%IMX,                     1,MPI_INTEGER,          0,1,MPI_COMM_WORLD,IERR)
   CALL MPI_SEND(MESHES(NM)%JMX,                     1,MPI_INTEGER,          0,1,MPI_COMM_WORLD,IERR)
   CALL MPI_SEND(MESHES(NM)%KMX,                     1,MPI_INTEGER,          0,1,MPI_COMM_WORLD,IERR)
   CALL MPI_SEND(MESHES(NM)%IMN,                     1,MPI_INTEGER,          0,1,MPI_COMM_WORLD,IERR)
   CALL MPI_SEND(MESHES(NM)%JMN,                     1,MPI_INTEGER,          0,1,MPI_COMM_WORLD,IERR)
   CALL MPI_SEND(MESHES(NM)%KMN,                     1,MPI_INTEGER,          0,1,MPI_COMM_WORLD,IERR)
   CALL MPI_SEND(MESHES(NM)%IRM,                     1,MPI_INTEGER,          0,1,MPI_COMM_WORLD,IERR)
   CALL MPI_SEND(MESHES(NM)%JRM,                     1,MPI_INTEGER,          0,1,MPI_COMM_WORLD,IERR)
   CALL MPI_SEND(MESHES(NM)%KRM,                     1,MPI_INTEGER,          0,1,MPI_COMM_WORLD,IERR)
   CALL MPI_SEND(MESHES(NM)%I_VN,                    1,MPI_INTEGER,          0,1,MPI_COMM_WORLD,IERR)
   CALL MPI_SEND(MESHES(NM)%J_VN,                    1,MPI_INTEGER,          0,1,MPI_COMM_WORLD,IERR)
   CALL MPI_SEND(MESHES(NM)%K_VN,                    1,MPI_INTEGER,          0,1,MPI_COMM_WORLD,IERR)
   CALL MPI_SEND(MESHES(NM)%NLP,                     1,MPI_INTEGER,          0,1,MPI_COMM_WORLD,IERR)
ENDDO

! Node 0 receives various values from all other nodes

DO NOM=1,NMESHES
   IF (PROCESS(NOM)==0 .OR. MYID/=0) CYCLE
   CALL MPI_RECV(MESHES(NOM)%DT,                      1,MPI_DOUBLE_PRECISION,PROCESS(NOM),1,MPI_COMM_WORLD,STATUS,IERR)
   CALL MPI_RECV(MESHES(NOM)%CFL,                     1,MPI_DOUBLE_PRECISION,PROCESS(NOM),1,MPI_COMM_WORLD,STATUS,IERR)
   CALL MPI_RECV(MESHES(NOM)%DIVMX,                   1,MPI_DOUBLE_PRECISION,PROCESS(NOM),1,MPI_COMM_WORLD,STATUS,IERR)
   CALL MPI_RECV(MESHES(NOM)%DIVMN,                   1,MPI_DOUBLE_PRECISION,PROCESS(NOM),1,MPI_COMM_WORLD,STATUS,IERR)
   CALL MPI_RECV(MESHES(NOM)%RESMAX,                  1,MPI_DOUBLE_PRECISION,PROCESS(NOM),1,MPI_COMM_WORLD,STATUS,IERR)
   CALL MPI_RECV(MESHES(NOM)%POIS_PTB,                1,MPI_DOUBLE_PRECISION,PROCESS(NOM),1,MPI_COMM_WORLD,STATUS,IERR)
   CALL MPI_RECV(MESHES(NOM)%POIS_ERR,                1,MPI_DOUBLE_PRECISION,PROCESS(NOM),1,MPI_COMM_WORLD,STATUS,IERR)
   CALL MPI_RECV(MESHES(NOM)%VN,                      1,MPI_DOUBLE_PRECISION,PROCESS(NOM),1,MPI_COMM_WORLD,STATUS,IERR)
   CALL MPI_RECV(MESHES(NOM)%ICFL,                    1,MPI_INTEGER,         PROCESS(NOM),1,MPI_COMM_WORLD,STATUS,IERR)
   CALL MPI_RECV(MESHES(NOM)%JCFL,                    1,MPI_INTEGER,         PROCESS(NOM),1,MPI_COMM_WORLD,STATUS,IERR)
   CALL MPI_RECV(MESHES(NOM)%KCFL,                    1,MPI_INTEGER,         PROCESS(NOM),1,MPI_COMM_WORLD,STATUS,IERR)
   CALL MPI_RECV(MESHES(NOM)%IMX,                     1,MPI_INTEGER,         PROCESS(NOM),1,MPI_COMM_WORLD,STATUS,IERR)
   CALL MPI_RECV(MESHES(NOM)%JMX,                     1,MPI_INTEGER,         PROCESS(NOM),1,MPI_COMM_WORLD,STATUS,IERR)
   CALL MPI_RECV(MESHES(NOM)%KMX,                     1,MPI_INTEGER,         PROCESS(NOM),1,MPI_COMM_WORLD,STATUS,IERR)
   CALL MPI_RECV(MESHES(NOM)%IMN,                     1,MPI_INTEGER,         PROCESS(NOM),1,MPI_COMM_WORLD,STATUS,IERR)
   CALL MPI_RECV(MESHES(NOM)%JMN,                     1,MPI_INTEGER,         PROCESS(NOM),1,MPI_COMM_WORLD,STATUS,IERR)
   CALL MPI_RECV(MESHES(NOM)%KMN,                     1,MPI_INTEGER,         PROCESS(NOM),1,MPI_COMM_WORLD,STATUS,IERR)
   CALL MPI_RECV(MESHES(NOM)%IRM,                     1,MPI_INTEGER,         PROCESS(NOM),1,MPI_COMM_WORLD,STATUS,IERR)
   CALL MPI_RECV(MESHES(NOM)%JRM,                     1,MPI_INTEGER,         PROCESS(NOM),1,MPI_COMM_WORLD,STATUS,IERR)
   CALL MPI_RECV(MESHES(NOM)%KRM,                     1,MPI_INTEGER,         PROCESS(NOM),1,MPI_COMM_WORLD,STATUS,IERR)
   CALL MPI_RECV(MESHES(NOM)%I_VN,                    1,MPI_INTEGER,         PROCESS(NOM),1,MPI_COMM_WORLD,STATUS,IERR)
   CALL MPI_RECV(MESHES(NOM)%J_VN,                    1,MPI_INTEGER,         PROCESS(NOM),1,MPI_COMM_WORLD,STATUS,IERR)
   CALL MPI_RECV(MESHES(NOM)%K_VN,                    1,MPI_INTEGER,         PROCESS(NOM),1,MPI_COMM_WORLD,STATUS,IERR)
   CALL MPI_RECV(MESHES(NOM)%NLP,                     1,MPI_INTEGER,         PROCESS(NOM),1,MPI_COMM_WORLD,STATUS,IERR)
ENDDO
 
TUSED(11,:) = TUSED(11,:) + SECOND() - TNOW
END SUBROUTINE EXCHANGE_DIAGNOSTICS


SUBROUTINE DUMP_GLOBAL_OUTPUTS(T)

! Dump HRR data to CHID_hrr.csv, MASS data to CHID_mass.csv, DEVICE data to _devc.csv

REAL(EB) :: T,TNOW
INTEGER :: N,CNT
INTEGER :: NM,DISP

TNOW = SECOND()

! Dump out HRR info  after first "gathering" data to node 0

DISP = DISPLS(MYID)+1
CNT  = COUNTS(MYID)

IF_DUMP_HRR: IF (T>=HRR_CLOCK) THEN
   IF (N_MPI_PROCESSES>1) THEN
      CALL MPI_ALLGATHERV(MPI_IN_PLACE,CNT,MPI_DOUBLE_PRECISION,HRR_TIME_INTERVAL,COUNTS,DISPLS,&
                          MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERR)
      IF (IERR/=MPI_SUCCESS) CALL HANDLE_MPI_ERROR('Error in MPI_ALLGATHERV for HRR_TIME_INTERVAL',IERR)
   ENDIF
   IF (MINVAL(HRR_TIME_INTERVAL,MASK=.NOT.EVACUATION_ONLY)>0._EB) THEN
      IF (N_MPI_PROCESSES>1) THEN
            REAL_BUFFER_11 = Q_DOT_SUM
            CALL MPI_GATHERV(REAL_BUFFER_11(1,DISP),CNT*N_Q_DOT,MPI_DOUBLE_PRECISION, &
                             Q_DOT_SUM,COUNTS_Q_DOT,DISPLS_Q_DOT,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
            IF (IERR/=MPI_SUCCESS) CALL HANDLE_MPI_ERROR('Error in MPI_GATHERV for Q_DOT_SUM',IERR)
            REAL_BUFFER_12 = M_DOT_SUM
            CALL MPI_GATHERV(REAL_BUFFER_12(1,DISP),CNT*N_M_DOT,MPI_DOUBLE_PRECISION, &
                             M_DOT_SUM,COUNTS_M_DOT,DISPLS_M_DOT,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
            IF (IERR/=MPI_SUCCESS) CALL HANDLE_MPI_ERROR('Error in MPI_GATHERV for M_DOT_SUM',IERR)
      ENDIF
      IF (MYID==0) CALL DUMP_HRR(T)
      HRR_CLOCK = HRR_CLOCK + DT_HRR
      Q_DOT_SUM = 0._EB
      M_DOT_SUM = 0._EB
      HRR_TIME_INTERVAL = 0._EB
   ENDIF
ENDIF IF_DUMP_HRR

! Dump unstructured geometry and boundary element info

IF (N_FACE>0 .AND. T>=GEOM_CLOCK) THEN
   IF (MYID==0) CALL DUMP_GEOM(T)
   GEOM_CLOCK = GEOM_CLOCK + DT_GEOM
ENDIF

IF (N_GEOM>0 .AND. T>=BNDC_CLOCK) THEN
   IF (MYID==0) CALL DUMP_BNDC(T)
   BNDC_CLOCK = BNDC_CLOCK + DT_BNDC
ENDIF

IF (N_BNDE>0 .AND. T>=BNDE_CLOCK) THEN
   IF (MYID==0) CALL DUMP_BNDE(T)
   BNDE_CLOCK = BNDE_CLOCK + DT_BNDE
ENDIF

! Dump out Evac info

IF (MYID==MAX(0,EVAC_PROCESS)) CALL EVAC_CSV(T)

! Dump out Mass info after first "gathering" data to node 0

IF_DUMP_MASS: IF (T>=MINT_CLOCK) THEN
   IF (N_MPI_PROCESSES>1) THEN
      CALL MPI_ALLGATHERV(MPI_IN_PLACE,CNT,MPI_DOUBLE_PRECISION,MINT_TIME_INTERVAL,COUNTS,DISPLS,&
                          MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERR)
      IF (IERR/=MPI_SUCCESS) CALL HANDLE_MPI_ERROR('Error in MPI_ALLGATHERV for MINT_TIME_INTERVAL',IERR)
   ENDIF
   IF (MINVAL(MINT_TIME_INTERVAL,MASK=.NOT.EVACUATION_ONLY)>0.) THEN
      IF (N_MPI_PROCESSES>1) THEN
         REAL_BUFFER_5 = MINT_SUM
         CALL MPI_GATHERV(REAL_BUFFER_5(0,DISP),COUNTS_MASS(MYID),MPI_DOUBLE_PRECISION, &
                          MINT_SUM,COUNTS_MASS,DISPLS_MASS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
         IF (IERR/=MPI_SUCCESS) CALL HANDLE_MPI_ERROR('Error in MPI_GATHERV for MINT_SUM',IERR)
      ENDIF
      IF (MYID==0) CALL DUMP_MASS(T)
      MINT_CLOCK    = MINT_CLOCK + DT_MASS
      MINT_SUM   = 0._EB
      MINT_TIME_INTERVAL = 0._EB
   ENDIF
ENDIF IF_DUMP_MASS

! Exchange DEVICE parameters among meshes and dump out DEVICE info after first "gathering" data to node 0
 
IF (N_DEVC>0) THEN
  
   ! Exchange the CURRENT_STATE of each DEViCe

   STATE_LOC(1:N_DEVC) = .FALSE.  ! _LOC is a temporary array that holds the STATE value for the devices on each node
   DO NM=1,NMESHES
      IF (PROCESS(NM)/=MYID) CYCLE
      DO N=1,N_DEVC
         IF (DEVICE(N)%MESH==NM) STATE_LOC(N) = DEVICE(N)%CURRENT_STATE 
      ENDDO
   ENDDO
   IF (N_MPI_PROCESSES>1) THEN
      CALL MPI_ALLREDUCE(STATE_LOC(1),STATE_GLB(1),N_DEVC,MPI_LOGICAL,MPI_LXOR,MPI_COMM_WORLD,IERR)
      IF (IERR/=MPI_SUCCESS) CALL HANDLE_MPI_ERROR('Error in MPI_ALLREDUCE for CURRENT_STATE',IERR)
   ELSE
      STATE_GLB = STATE_LOC
   ENDIF
   DEVICE(1:N_DEVC)%CURRENT_STATE = STATE_GLB(1:N_DEVC)

   DEVC_PIPE_OPERATING = 0
   DO CNT=1,N_DEVC
      IF (DEVICE(CNT)%PROP_INDEX > 0 .AND.  DEVICE(CNT)%CURRENT_STATE) THEN
         IF (PROPERTY(DEVICE(CNT)%PROP_INDEX)%PART_INDEX > 0) DEVC_PIPE_OPERATING(DEVICE(CNT)%PIPE_INDEX) = &
            DEVC_PIPE_OPERATING(DEVICE(CNT)%PIPE_INDEX) + 1
      ENDIF
   ENDDO
   ! Exchange the PRIOR_STATE of each DEViCe

   STATE_LOC(1:N_DEVC) = .FALSE.  ! _LOC is a temporary array that holds the STATE value for the devices on each node
   DO NM=1,NMESHES
      IF (PROCESS(NM)/=MYID) CYCLE
      DO N=1,N_DEVC
         IF (DEVICE(N)%MESH==NM) STATE_LOC(N) = DEVICE(N)%PRIOR_STATE
      ENDDO
   ENDDO
   IF (N_MPI_PROCESSES>1) THEN
      CALL MPI_ALLREDUCE(STATE_LOC(1),STATE_GLB(1),N_DEVC,MPI_LOGICAL,MPI_LXOR,MPI_COMM_WORLD,IERR)
      IF (IERR/=MPI_SUCCESS) CALL HANDLE_MPI_ERROR('Error in MPI_ALLREDUCE for PRIOR_STATE',IERR)
   ELSE
      STATE_GLB = STATE_LOC
   ENDIF
   DEVICE(1:N_DEVC)%PRIOR_STATE = STATE_GLB(1:N_DEVC)

   ! Exchange the INSTANT_VALUE of each DEViCe

   TC_LOC(1:N_DEVC) = 0._EB 
   DO NM=1,NMESHES
      IF (PROCESS(NM)/=MYID) CYCLE
      DO N=1,N_DEVC
         IF (DEVICE(N)%MESH==NM) TC_LOC(N) = DEVICE(N)%INSTANT_VALUE
      ENDDO
   ENDDO
   IF (N_MPI_PROCESSES>1) THEN
      CALL MPI_ALLREDUCE(TC_LOC(1),TC_GLB(1),N_DEVC,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
      IF (IERR/=MPI_SUCCESS) CALL HANDLE_MPI_ERROR('Error in MPI_ALLREDUCE for INSTANT_VALUE',IERR)
   ELSE
      TC_GLB = TC_LOC
   ENDIF
   DEVICE(1:N_DEVC)%INSTANT_VALUE = TC_GLB(1:N_DEVC)

   ! Exchange the SMOOTHED_VALUE of each DEViCe

   TC_LOC(1:N_DEVC) = 0._EB 
   DO NM=1,NMESHES
      IF (PROCESS(NM)/=MYID) CYCLE
      DO N=1,N_DEVC
         IF (DEVICE(N)%MESH==NM) TC_LOC(N) = DEVICE(N)%SMOOTHED_VALUE
      ENDDO
   ENDDO
   IF (N_MPI_PROCESSES>1) THEN
      CALL MPI_ALLREDUCE(TC_LOC(1),TC_GLB(1),N_DEVC,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
      IF (IERR/=MPI_SUCCESS) CALL HANDLE_MPI_ERROR('Error in MPI_ALLREDUCE for SMOOTHED_VALUE',IERR)
   ELSE
      TC_GLB = TC_LOC
   ENDIF
   DEVICE(1:N_DEVC)%SMOOTHED_VALUE = TC_GLB(1:N_DEVC)


   ! Exchange the T_CHANGE of each DEViCe

   TC_LOC(1:N_DEVC) = 0._EB
   DO NM=1,NMESHES
      IF (PROCESS(NM)/=MYID) CYCLE
      DO N=1,N_DEVC
         IF (DEVICE(N)%MESH==NM) TC_LOC(N) = DEVICE(N)%T_CHANGE
      ENDDO
   ENDDO
   IF (N_MPI_PROCESSES>1) THEN
      CALL MPI_ALLREDUCE(TC_LOC(1),TC_GLB(1),N_DEVC,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
      IF (IERR/=MPI_SUCCESS) CALL HANDLE_MPI_ERROR('Error in MPI_ALLREDUCE for T_CHANGE',IERR)
   ELSE
      TC_GLB = TC_LOC
   ENDIF
   DEVICE(1:N_DEVC)%T_CHANGE = TC_GLB(1:N_DEVC)

ENDIF

! Exchange information about Devices that is only needed at print-out time

IF_DUMP_DEVC: IF (T>=DEVC_CLOCK .AND. N_DEVC>0) THEN

   ! Exchange the current COUNT of each DEViCe

   TI_LOC(1:N_DEVC) = DEVICE(1:N_DEVC)%TIME_INTERVAL
   IF (N_MPI_PROCESSES>1) THEN
      CALL MPI_ALLREDUCE(TI_LOC(1),TI_GLB(1),N_DEVC,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
      IF (IERR/=MPI_SUCCESS) CALL HANDLE_MPI_ERROR('Error in MPI_ALLREDUCE for TIME_INTERVAL',IERR)
   ELSE
      TI_GLB = TI_LOC
   ENDIF

   ! Get the current VALUEs of all DEViCes into DEVICE(:)%VALUE on node 0

   IF (MINVAL(TI_GLB)>0._EB) THEN
      TC_LOC(1:N_DEVC) = DEVICE(1:N_DEVC)%VALUE 
      IF (N_MPI_PROCESSES>1) THEN
         CALL MPI_REDUCE(TC_LOC(1),TC_GLB(1),N_DEVC,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
         IF (IERR/=MPI_SUCCESS) CALL HANDLE_MPI_ERROR('Error in MPI_REDUCE for VALUE',IERR)
      ELSE
         TC_GLB = TC_LOC
      ENDIF
      IF (MYID==0) THEN
         DEVICE(1:N_DEVC)%VALUE         = TC_GLB(1:N_DEVC)
         DEVICE(1:N_DEVC)%TIME_INTERVAL = TI_GLB(1:N_DEVC)
         CALL DUMP_DEVICES(T)
      ENDIF
      DEVC_CLOCK = DEVC_CLOCK + DT_DEVC
      DO N=1,N_DEVC
         DEVICE(N)%VALUE = 0._EB
         DEVICE(N)%TIME_INTERVAL = 0._EB
      ENDDO
   ENDIF

ENDIF IF_DUMP_DEVC

! Dump CONTROL info. No gathering required as CONTROL is updated on all meshes

IF (T>=CTRL_CLOCK .AND. N_CTRL>0) THEN
   IF (MYID==0) CALL DUMP_CONTROLS(T)
   CTRL_CLOCK = CTRL_CLOCK + DT_CTRL
ENDIF

TUSED(7,:) = TUSED(7,:) + SECOND() - TNOW
END SUBROUTINE DUMP_GLOBAL_OUTPUTS


SUBROUTINE INITIALIZE_EVAC

! Initialize evacuation meshes
 
DO NM=1,NMESHES
   IF (USE_MPI .AND. N_MPI_PROCESSES>1 .AND. MYID==EVAC_PROCESS .AND. .NOT.EVACUATION_ONLY(NM)) THEN
      M=>MESHES(NM)
      !EVACUATION: SOLID, CELL_INDEX, OBST_INDEX_C, OBSTRUCTION are allocated in READ_OBST for the evac process.
      ALLOCATE(M%ZZ(0:M%IBP1,0:M%JBP1,0:M%KBP1,N_TRACKED_SPECIES),STAT=IZERO)
      CALL ChkMemErr('MAIN','Evac ZZ',IZERO)
      M%ZZ=0._EB
      ALLOCATE(M%RHO(0:M%IBP1,0:M%JBP1,0:M%KBP1),STAT=IZERO)
      CALL ChkMemErr('MAIN','Evac RHO',IZERO)
      M%RHO=RHOA
      ALLOCATE(M%RSUM(0:M%IBP1,0:M%JBP1,0:M%KBP1),STAT=IZERO)
      CALL ChkMemErr('MAIN','Evac RSUM',IZERO)
      M%RSUM=RSUM0
      ALLOCATE(M%TMP(0:M%IBP1,0:M%JBP1,0:M%KBP1),STAT=IZERO)
      CALL ChkMemErr('MAIN','Evac TMP',IZERO)
      M%TMP=TMPA
      ALLOCATE(M%UII(0:M%IBP1,0:M%JBP1,0:M%KBP1),STAT=IZERO)
      CALL ChkMemErr('MAIN','Evac UII',IZERO)
      M%UII=4._EB*SIGMA*TMPA4
   ENDIF
   IF (PROCESS(NM)/=MYID) CYCLE
   IF (EVACUATION_GRID(NM)) PART_CLOCK(NM) = T_EVAC + DT_PART
   IF (MYID/=MAX(0,EVAC_PROCESS)) CYCLE
   IF (ANY(EVACUATION_GRID)) CALL INITIALIZE_EVACUATION(NM)
   IF (EVACUATION_GRID(NM)) CALL DUMP_EVAC(T_EVAC,NM)
ENDDO
IF (ANY(EVACUATION_GRID) .AND. .NOT.RESTART) ICYC = -EVAC_TIME_ITERATIONS

END SUBROUTINE INITIALIZE_EVAC

SUBROUTINE INIT_EVAC_DUMPS

! Initialize evacuation dumps

REAL(EB) :: T_TMP

IF (.NOT.ANY(EVACUATION_ONLY)) RETURN ! No evacuation
 
IF (RESTART) THEN
   T_TMP = MINVAL(T,MASK=.NOT.EVACUATION_ONLY)
   T_EVAC_SAVE = T_TMP
ELSE
   T_EVAC  = - EVAC_DT_FLOWFIELD*EVAC_TIME_ITERATIONS + T_BEGIN
   T_EVAC_SAVE = T_EVAC
   T_TMP = T_EVAC
END IF
IF (.NOT.ANY(EVACUATION_GRID)) RETURN ! No main evacuation meshes
IF ((.NOT.USE_MPI .OR. N_MPI_PROCESSES==1) .OR. (N_MPI_PROCESSES>1 .AND. MYID==EVAC_PROCESS)) &
     CALL INITIALIZE_EVAC_DUMPS(T_TMP,T_EVAC_SAVE)

END SUBROUTINE INIT_EVAC_DUMPS


SUBROUTINE EVAC_CSV(T)
 
! Dump out Evac info

REAL(EB), INTENT(IN) :: T

IF (T>=EVAC_CLOCK .AND. ANY(EVACUATION_GRID)) THEN
   CALL DUMP_EVAC_CSV(T)
   EVAC_CLOCK = EVAC_CLOCK + DT_HRR
ENDIF

END SUBROUTINE EVAC_CSV


SUBROUTINE EVAC_EXCHANGE

LOGICAL EXCHANGE_EVACUATION
INTEGER NM, II, IVENT, I, J, EMESH, JJ, N_END
 
! Fire mesh information ==> Evac meshes
 
IF (.NOT. ANY(EVACUATION_GRID)) RETURN
IF (N_MPI_PROCESSES>1 .AND. MYID /= EVAC_PROCESS) CALL EVAC_MESH_EXCHANGE(T_EVAC,T_EVAC_SAVE,I_EVAC,ICYC,EXCHANGE_EVACUATION,1)
IF ((.NOT.USE_MPI .OR. N_MPI_PROCESSES==1) .OR. (N_MPI_PROCESSES>1 .AND. MYID==EVAC_PROCESS)) &
     CALL EVAC_MESH_EXCHANGE(T_EVAC,T_EVAC_SAVE,I_EVAC,ICYC,EXCHANGE_EVACUATION,2)

! Update evacuation devices

DO NM=1,NMESHES
   IF (ACTIVE_MESH(NM)) CYCLE
   IF (.NOT.EVACUATION_GRID(NM)) CYCLE
   IF (MYID/=MAX(0,EVAC_PROCESS)) CYCLE
   CALL UPDATE_GLOBAL_OUTPUTS(T(NM),NM)      
ENDDO

! Save the evacuation flow fields to the arrays U_EVAC and V_EVAC

N_END = N_EXITS - N_CO_EXITS + N_DOORS
DO NM = 1, NMESHES
   IF (.NOT.ACTIVE_MESH(NM)) CYCLE
   IF (.NOT.EVACUATION_GRID(NM)) CYCLE
   IF (MYID /= MAX(0,EVAC_PROCESS)) CYCLE
   II = EVAC_TIME_ITERATIONS / MAXVAL(EMESH_NFIELDS)
   IF (MOD(ABS(ICYC),II)==0) THEN
      IVENT = (ABS(ICYC))/II + 1
      LOOP_EXITS: DO JJ = 1, N_END
         IF (EMESH_EXITS(JJ)%MAINMESH == NM .AND. EMESH_EXITS(JJ)%I_DOORS_EMESH == IVENT) THEN
            EMESH = EMESH_EXITS(JJ)%EMESH
            DO J = 0, EMESH_IJK(2,EMESH) + 1
               DO I = 0, EMESH_IJK(1,EMESH) + 1
                  EMESH_EXITS(JJ)%U_EVAC(I,J) = MESHES(NM)%U(I,J,1)
                  EMESH_EXITS(JJ)%V_EVAC(I,J) = MESHES(NM)%V(I,J,1)
               END DO
            END DO
            EXIT LOOP_EXITS
         END IF
      END DO LOOP_EXITS
   END IF
   
ENDDO

END SUBROUTINE EVAC_EXCHANGE


SUBROUTINE EVAC_PRESSURE_ITERATION_SCHEME
 
! Evacuation flow field calculation
 
INTEGER :: N
 
COMPUTE_PRESSURE_LOOP: DO NM=1,NMESHES
   IF (.NOT.EVACUATION_ONLY(NM))  CYCLE COMPUTE_PRESSURE_LOOP
   IF (PROCESS(NM)/=MYID)         CYCLE COMPUTE_PRESSURE_LOOP
   IF (.NOT.ACTIVE_MESH(NM))      CYCLE COMPUTE_PRESSURE_LOOP
   PRESSURE_ITERATION_LOOP: DO N=1,EVAC_PRESSURE_ITERATIONS
      CALL NO_FLUX(NM)
      MESHES(NM)%FVZ = 0._EB
      CALL PRESSURE_SOLVER(T(NM),NM)
   ENDDO PRESSURE_ITERATION_LOOP
ENDDO COMPUTE_PRESSURE_LOOP

END SUBROUTINE EVAC_PRESSURE_ITERATION_SCHEME


SUBROUTINE EVAC_MAIN_LOOP

! Call the evacuation routine and adjust the time steps for the evacuation meshes
 
REAL(EB) :: T_FIRE, FIRE_DT, EVAC_DT
INTEGER :: II
 
IF (.NOT. ANY(EVACUATION_GRID)) RETURN
 
IF (ANY(EVACUATION_ONLY) .AND. (ICYC <= 0)) THEN
   ACTIVE_MESH = .FALSE.  ! Be sure that no fire meshes are updated for icyc < 0
ENDIF
EVAC_DT = EVAC_DT_STEADY_STATE
IF (ICYC < 1) EVAC_DT = EVAC_DT_FLOWFIELD
T_FIRE = T_EVAC + EVAC_DT
IF (ICYC > 0) THEN ! Syncrhonize evacuation and fire clocks, if both type meshes present
   IF (.NOT.ALL(EVACUATION_ONLY)) THEN
      T_FIRE = MINVAL(T, MASK=(.NOT.EVACUATION_ONLY) .AND. ACTIVE_MESH)
      DO NM=1,NMESHES
         IF (PROCESS(NM)==MYID) DT_NEXT_SYNC(NM) = MESHES(NM)%DT_NEXT
      ENDDO

      IF (N_MPI_PROCESSES>1) THEN
         CALL MPI_ALLGATHERV(MPI_IN_PLACE,COUNTS(MYID),MPI_DOUBLE_PRECISION, &
                             DT_NEXT_SYNC,COUNTS,DISPLS,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERR)
         IF (IERR/=MPI_SUCCESS) CALL HANDLE_MPI_ERROR('Error in MPI_ALLGATHERV for DT_NEXT_SYNC',IERR)
      ENDIF

      FIRE_DT = MINVAL(DT_NEXT_SYNC, MASK=(.NOT.EVACUATION_ONLY) .AND. ACTIVE_MESH)
      T_FIRE = T_FIRE + FIRE_DT
   ENDIF
ENDIF
EVAC_TIME_STEP_LOOP: DO WHILE (T_EVAC < T_FIRE)
   T_EVAC = T_EVAC + EVAC_DT
   IF ((.NOT.USE_MPI .OR. N_MPI_PROCESSES==1) .OR. (N_MPI_PROCESSES>1 .AND. MYID==EVAC_PROCESS)) CALL PREPARE_TO_EVACUATE(ICYC)
   DO NM = 1, NMESHES
      IF (EVACUATION_ONLY(NM)) THEN
         ACTIVE_MESH(NM)      = .FALSE.
         CHANGE_TIME_STEP(NM) = .FALSE.
         MESHES(NM)%DT        = EVAC_DT
         MESHES(NM)%DT_NEXT   = EVAC_DT
         T(NM)                = T_EVAC
         IF (ICYC <= 1 .AND. .NOT.BTEST(I_EVAC, 2)) THEN
            IF (ICYC <= 0 .AND. EVACUATION_GRID(NM)) THEN
               II = EVAC_TIME_ITERATIONS / MAXVAL(EMESH_NFIELDS)
               IF ((ABS(ICYC)+1) <= EMESH_NFIELDS(EMESH_INDEX(NM))*II) THEN
                  ACTIVE_MESH(NM) = .TRUE.
               END IF
               DIAGNOSTICS = .FALSE.
            END IF
            !
            IF (ICYC <= 0) T(NM) = T_EVAC + EVAC_DT_FLOWFIELD*EVAC_TIME_ITERATIONS - EVAC_DT_FLOWFIELD
         ENDIF
         IF (ICYC <= 1 .AND. BTEST(I_EVAC, 2)) THEN
            IF (ICYC <= 0 .AND. EVACUATION_GRID(NM)) THEN
               DIAGNOSTICS = .FALSE.
            END IF
            IF (ICYC <= 0) T(NM) = T_EVAC + EVAC_DT_FLOWFIELD*EVAC_TIME_ITERATIONS - EVAC_DT_FLOWFIELD
         ENDIF
         IF (.NOT.ACTIVE_MESH(NM)) THEN
            VELOCITY_ERROR_MAX_I(NM) = 1
            VELOCITY_ERROR_MAX_J(NM) = 1
            VELOCITY_ERROR_MAX_K(NM) = 1
            VELOCITY_ERROR_MAX(NM) = 0._EB
            MESHES(NM)%POIS_ERR = 0.0_EB
            MESHES(NM)%POIS_PTB = 0.0_EB
            MESHES(NM)%RESMAX = 0.0_EB
            MESHES(NM)%CFL = 0.0_EB
            MESHES(NM)%ICFL = 0; MESHES(NM)%JCFL = 0; MESHES(NM)%KCFL = 0
            MESHES(NM)%DIVMX = 0.0_EB
            MESHES(NM)%IMX = 0; MESHES(NM)%JMX = 0; MESHES(NM)%KMX = 0
            MESHES(NM)%DIVMN = 0.0_EB
            MESHES(NM)%IMN = 0; MESHES(NM)%JMN = 0; MESHES(NM)%KMN = 0
         END IF
         IF (EVACUATION_GRID(NM) ) THEN
            IF (PROCESS(NM)==MYID .AND. STOP_STATUS==NO_STOP) CALL EVACUATE_HUMANS(T_EVAC,NM,ICYC)
            IF (PROCESS(NM)==MYID .AND. .NOT.ACTIVE_MESH(NM)) NTCYC(NM) = NTCYC(NM) + 1
            IF (T_EVAC >= PART_CLOCK(NM)) THEN
               IF (PROCESS(NM)==MYID) CALL DUMP_EVAC(T_EVAC, NM)
               DO
                  PART_CLOCK(NM) = PART_CLOCK(NM) + DT_PART
                  IF (PART_CLOCK(NM) >= T_EVAC) EXIT
               ENDDO
            ENDIF
         ENDIF
      ENDIF
   ENDDO
   IF (ICYC < 1) EXIT EVAC_TIME_STEP_LOOP
   IF ((.NOT.USE_MPI .OR. N_MPI_PROCESSES==1) .OR. (N_MPI_PROCESSES>1 .AND. MYID==EVAC_PROCESS)) &
        CALL CLEAN_AFTER_EVACUATE(ICYC, I_EVAC)
ENDDO EVAC_TIME_STEP_LOOP
IF (ICYC < 1 .AND. MYID==0) THEN
   ! Write the diagnostic information for the evacuation mesh initialization time steps
   II = EVAC_TIME_ITERATIONS / MAXVAL(EMESH_NFIELDS)
   IF (MOD(ABS(ICYC)+1,II) == 0 .OR. ABS(ICYC)+1 == EVAC_TIME_ITERATIONS) THEN
      WRITE(LU_ERR,'(1X,A,I7,A,F10.3,A)')  'Time Step:',ICYC,',    Evacuation Initialization Time:',T_EVAC,' s'
   END IF
END IF

END SUBROUTINE EVAC_MAIN_LOOP


SUBROUTINE EXCHANGE_HVAC_BC

! Exchange information mesh to mesh needed for performing the HVAC computation

USE HVAC_ROUTINES, ONLY: NODE_H,NODE_P,NODE_RHO,NODE_TMP,NODE_X,NODE_Y,NODE_Z,NODE_ZZ
REAL(EB) :: TNOW

TNOW = SECOND()

REAL_BUFFER_6 = NODE_H
CALL MPI_GATHERV(REAL_BUFFER_6(1,DISPLS(MYID)+1),COUNTS_HVAC(MYID),MPI_DOUBLE_PRECISION, &
                 NODE_H,COUNTS_HVAC,DISPLS_HVAC,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
REAL_BUFFER_6 = NODE_P
CALL MPI_GATHERV(REAL_BUFFER_6(1,DISPLS(MYID)+1),COUNTS_HVAC(MYID),MPI_DOUBLE_PRECISION, &
                 NODE_P,COUNTS_HVAC,DISPLS_HVAC,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
REAL_BUFFER_6 = NODE_RHO
CALL MPI_GATHERV(REAL_BUFFER_6(1,DISPLS(MYID)+1),COUNTS_HVAC(MYID),MPI_DOUBLE_PRECISION, &
                 NODE_RHO,COUNTS_HVAC,DISPLS_HVAC,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
REAL_BUFFER_6 = NODE_TMP
CALL MPI_GATHERV(REAL_BUFFER_6(1,DISPLS(MYID)+1),COUNTS_HVAC(MYID),MPI_DOUBLE_PRECISION, &
                 NODE_TMP,COUNTS_HVAC,DISPLS_HVAC,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
REAL_BUFFER_6 = NODE_X
CALL MPI_GATHERV(REAL_BUFFER_6(1,DISPLS(MYID)+1),COUNTS_HVAC(MYID),MPI_DOUBLE_PRECISION, &
                 NODE_X,COUNTS_HVAC,DISPLS_HVAC,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
REAL_BUFFER_6 = NODE_Y
CALL MPI_GATHERV(REAL_BUFFER_6(1,DISPLS(MYID)+1),COUNTS_HVAC(MYID),MPI_DOUBLE_PRECISION, &
                 NODE_Y,COUNTS_HVAC,DISPLS_HVAC,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
REAL_BUFFER_6 = NODE_Z
CALL MPI_GATHERV(REAL_BUFFER_6(1,DISPLS(MYID)+1),COUNTS_HVAC(MYID),MPI_DOUBLE_PRECISION, &
                 NODE_Z,COUNTS_HVAC,DISPLS_HVAC,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
REAL_BUFFER_6 = NODE_AREA
CALL MPI_GATHERV(REAL_BUFFER_6(1,DISPLS(MYID)+1),COUNTS_HVAC(MYID),MPI_DOUBLE_PRECISION, &
                 NODE_AREA,COUNTS_HVAC,DISPLS_HVAC,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
INTEGER_BUFFER_2 = NODE_ZONE
CALL MPI_GATHERV(INTEGER_BUFFER_2(1,DISPLS(MYID)+1),COUNTS_HVAC(MYID),MPI_INTEGER, &
                 NODE_ZONE,COUNTS_HVAC,DISPLS_HVAC,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
REAL_BUFFER_7 = NODE_ZZ
CALL MPI_GATHERV(REAL_BUFFER_7(1,1,DISPLS(MYID)+1),COUNTS_HVAC_SPECIES(MYID),MPI_DOUBLE_PRECISION, &
                 NODE_ZZ,COUNTS_HVAC_SPECIES,DISPLS_HVAC_SPECIES,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

TUSED(2,:)=TUSED(2,:) + SECOND() - TNOW

END SUBROUTINE EXCHANGE_HVAC_BC


SUBROUTINE EXCHANGE_HVAC_SOLUTION

! Exchange information mesh to mesh needed for performing the HVAC computation

USE HVAC_ROUTINES, ONLY: NODE_AREA_EX,NODE_TMP_EX,NODE_ZZ_EX,DUCT_MF
REAL(EB) :: TNOW

TNOW = SECOND()

REAL_BUFFER_6 = NODE_AREA_EX
CALL MPI_SCATTER(REAL_BUFFER_6(1,1),N_DUCTNODES,MPI_DOUBLE_PRECISION, &
                 NODE_AREA_EX,N_DUCTNODES,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
REAL_BUFFER_6 = NODE_TMP_EX
CALL MPI_SCATTER(REAL_BUFFER_6(1,1),N_DUCTNODES,MPI_DOUBLE_PRECISION, &
                 NODE_TMP_EX,N_DUCTNODES,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
REAL_BUFFER_7 = NODE_ZZ_EX
CALL MPI_SCATTER(REAL_BUFFER_7(1,1,1),N_DUCTNODES*N_TOTAL_SCALARS,MPI_DOUBLE_PRECISION, &
                 NODE_ZZ_EX,N_DUCTNODES*N_TOTAL_SCALARS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
REAL_BUFFER_10 = DUCT_MF
CALL MPI_SCATTER(REAL_BUFFER_10(1,1),N_DUCTS,MPI_DOUBLE_PRECISION, &
                 DUCT_MF,N_DUCTS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

TUSED(2,:)=TUSED(2,:) + SECOND() - TNOW

END SUBROUTINE EXCHANGE_HVAC_SOLUTION

SUBROUTINE GET_INFO (REVISION,REVISION_DATE,COMPILE_DATE)
CHARACTER(LEN=255), INTENT(OUT) :: REVISION, REVISION_DATE, COMPILE_DATE

! Unlike svn, the revisioning system git does not perform keyword substitution.
! To perform this function,  a script named expand_file is called before FDS is
! built that expands the following keywords ($Revision, $RevisionDate and
! $CompileDate) with their proper values. Another script named contract_file is
! called after FDS is built to return these keywords back to their original
! values (so the revisioning system will not think this file has changed).

CHARACTER(255), PARAMETER :: GREVISION='$Revision$'
CHARACTER(255), PARAMETER :: GREVISION_DATE='$RevisionDate: unknown $'
CHARACTER(255), PARAMETER :: GCOMPILE_DATE='$CompileDate: unknown $'

WRITE(REVISION,'(A)')      GREVISION(INDEX(GREVISION,':')+2:LEN_TRIM(GREVISION)-2)
WRITE(REVISION_DATE,'(A)') GREVISION_DATE(INDEX(GREVISION_DATE,':')+2:LEN_TRIM(GREVISION_DATE)-2)
WRITE(COMPILE_DATE,'(A)')  GCOMPILE_DATE(INDEX(GCOMPILE_DATE,':')+2:LEN_TRIM(GCOMPILE_DATE)-2)
RETURN
END SUBROUTINE GET_INFO

END PROGRAM FDS
