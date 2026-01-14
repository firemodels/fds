!> \brief Fire Dynamics Simulator (FDS) is a computational fluid dynamics (CFD) code designed to model
!> fire and other thermal phenomena.

PROGRAM FDS

USE PRECISION_PARAMETERS
USE MESH_VARIABLES
USE GLOBAL_CONSTANTS
USE OUTPUT_CLOCKS
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
USE MISC_FUNCTIONS, ONLY : WRITE_SUMMARY_INFO,VERBOSE_PRINTOUT
USE OUTPUT_DATA
USE MEMORY_FUNCTIONS
USE HVAC_ROUTINES
USE COMP_FUNCTIONS, ONLY : CURRENT_TIME,INITIALIZE_OUTPUT_CLOCKS,READ_EXTERNAL_FILE,GET_FILE_NUMBER
USE DEVICE_VARIABLES
USE WALL_ROUTINES
USE FIRE
USE CHEMCONS, ONLY: WRITE_CVODE_SUBSTEPS
USE CONTROL_FUNCTIONS
USE TURBULENCE, ONLY: NS_ANALYTICAL_SOLUTION,INIT_TURB_ARRAYS,COMPRESSION_WAVE,&
                      TWOD_VORTEX_CERFACS,TWOD_VORTEX_UMD,TWOD_SOBOROT_UMD, &
                      SYNTHETIC_TURBULENCE,SYNTHETIC_EDDY_SETUP,SANDIA_DAT
USE MANUFACTURED_SOLUTIONS, ONLY: SHUNN_MMS_3,SAAD_MMS_1
USE CC_SCALARS,   ONLY: CC_SET_DATA,CC_END_STEP,CC_DENSITY,CC_RHO0W_INTERP,       &
                        CCCOMPUTE_RADIATION,CC_NO_FLUX,CC_COMPUTE_VELOCITY_ERROR, &
                        CC_NO_FLUX,CC_COMPUTE_VELOCITY_ERROR,FINISH_CC,        &
                        INIT_CUTCELL_DATA,MESH_CC_EXCHANGE,ROTATED_CUBE_ANN_SOLN, &
                        CC_RESTORE_UVW_UNLINKED
USE OPENMP_FDS
USE MPI_F08
USE SOOT_ROUTINES, ONLY: CALC_AGGLOMERATION
USE GLOBMAT_SOLVER, ONLY : GLMAT_SOLVER_SETUP, GLMAT_SOLVER, COPY_H_OMESH_TO_MESH, &
                           FINISH_GLMAT_SOLVER,PRESSURE_SOLVER_CHECK_RESIDUALS_U
USE LOCMAT_SOLVER, ONLY : ULMAT_SOLVER,ULMAT_SOLVER_SETUP,FINISH_ULMAT_SOLVER

IMPLICIT NONE (TYPE,EXTERNAL)

! Miscellaneous declarations

LOGICAL  :: EX=.FALSE.,DIAGNOSTICS,CTRL_STOP_STATUS,CHECK_FREEZE_VELOCITY=.TRUE.,EXTERNAL_FAIL
INTEGER  :: LO10,NM,IZERO,ANG_INC_COUNTER
REAL(EB) :: T,DT,TNOW
REAL :: CPUTIME
REAL(EB), ALLOCATABLE, DIMENSION(:) ::  TC_ARRAY,DT_NEW
REAL(EB), ALLOCATABLE, DIMENSION(:,:) ::  TC2_ARRAY
LOGICAL, ALLOCATABLE, DIMENSION(:) ::  STATE_ARRAY
INTEGER :: ITER
TYPE (MESH_TYPE), POINTER :: M,M4
TYPE (OMESH_TYPE), POINTER :: M2,M3

! MPI stuff

TYPE (MPI_STATUS) :: STATUS
INTEGER :: N,I,IERR=0
INTEGER :: PNAMELEN=0
INTEGER :: PROVIDED
INTEGER, PARAMETER :: REQUIRED=MPI_THREAD_FUNNELED
TYPE (MPI_REQUEST), ALLOCATABLE, DIMENSION(:) :: REQ,REQ1,REQ2,REQ3,REQ4,REQ5,REQ7,REQ6,REQ14,REQ15
INTEGER :: N_REQ=0,N_REQ1=0,N_REQ2=0,N_REQ3=0,N_REQ4=0,N_REQ5=0,N_REQ7=0,N_REQ6=0,N_REQ14=0,N_REQ15=0
CHARACTER(MPI_MAX_PROCESSOR_NAME) :: PNAME
LOGICAL, ALLOCATABLE, DIMENSION(:)        :: LOGICAL_BUFFER_EXTERNAL
REAL(EB), ALLOCATABLE, DIMENSION(:)       :: REAL_BUFFER_DUCT,REAL_BUFFER_EXTERNAL
REAL(EB), ALLOCATABLE, DIMENSION(:,:)     :: REAL_BUFFER_10,REAL_BUFFER_20

! Initialize OpenMP

CALL OPENMP_INIT

! Output version info if fds is invoked without any arguments. This must be done before MPI is initialized.

CALL VERSION_INFO

! Initialize MPI

CALL MPI_INIT_THREAD(REQUIRED,PROVIDED,IERR)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, MY_RANK, IERR)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, N_MPI_PROCESSES, IERR)
CALL MPI_GET_PROCESSOR_NAME(PNAME, PNAMELEN, IERR)

! Write out MPI process info to standard error (LU_ERR=0)

IF (MY_RANK==0) WRITE(LU_ERR,'(/A/)') ' Starting FDS ...'

CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
WRITE(LU_ERR,'(A,I6,A,A)') ' MPI Process ',MY_RANK,' started on ',PNAME(1:PNAMELEN)
CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

! Check that MPI processes and OpenMP threads are working properly

CALL CHECK_MPI

! Start wall clock timing

WALL_CLOCK_START = CURRENT_TIME()
CALL CPU_TIME(CPUTIME)
CPU_TIME_START = CPUTIME
ALLOCATE(T_USED(N_TIMERS)) ; T_USED = 0._EB ; T_USED(1) = CURRENT_TIME()

! Assign a compilation date

CALL GET_INFO(REVISION,REVISION_DATE,COMPILE_DATE)

! Read input from CHID.fds file and stop the code if any errors are found

CALL READ_DATA(DT) ; CALL STOP_CHECK(1)

IF (MY_RANK==0) THEN
   CALL WRITE_SUMMARY_INFO(LU_ERR,.TRUE.)
   WRITE(LU_ERR,'(/A,A)')     ' Job TITLE        : ',TRIM(TITLE)
   WRITE(LU_ERR,'(A,A/)')     ' Job ID string    : ',TRIM(CHID)
   IF (VERBOSE) CALL VERBOSE_PRINTOUT('Completed reading of input file')
ENDIF

! Print OpenMP thread status

CALL OPENMP_PRINT_STATUS

! Set up send and receive buffer counts and displacements

CALL MPI_INITIALIZATION_CHORES(1)

! Store max cell aspect ratio

CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,MAX_CELL_ASPECT_RATIO(1:NMESHES),&
                    COUNTS,DISPLS,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERR)
IF (MAXVAL(MAX_CELL_ASPECT_RATIO)>3.99_EB .AND. .NOT.CFL_VELOCITY_NORM_USER_SPECIFIED) CFL_VELOCITY_NORM=1

! Create output file names

CALL ASSIGN_FILE_NAMES

! Start the clock

T = T_BEGIN

! Allocate various utility arrays

CALL MPI_INITIALIZATION_CHORES(2)

! Exchange VENT data

IF (N_VENT_TOTAL > 0) CALL EXCHANGE_VENT_AREA

! Initialize global parameters

CALL INITIALIZE_GLOBAL_VARIABLES

! Initialize radiation

IF (RADIATION) CALL INIT_RADIATION

! Set up persistent MPI communicators for neighborhoods

CALL MPI_INITIALIZATION_CHORES(3)

CALL EXCHANGE_GEOMETRY_INFO

! Allocate and initialized external RAMP and CTRL

IF (READ_EXTERNAL) THEN
   T_EXTERNAL = T_BEGIN + DT_EXTERNAL
   LU_EXTERNAL  = GET_FILE_NUMBER()
   IF (DT_EXTERNAL_HEARTBEAT > 0._EB) LU_EXTERNAL_HEARTBEAT = GET_FILE_NUMBER()
ENDIF

! Set up the background atmosphere and initialize the boundary (WALL) arrays

DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   CALL INITIALIZE_ATMOSPHERE(NM)
   IF (.NOT.SETUP_ONLY .OR. CHECK_MESH_ALIGNMENT) CALL INITIALIZE_WALL_ARRAY(NM)
ENDDO
IF (MY_RANK==0 .AND. VERBOSE) CALL VERBOSE_PRINTOUT('Completed INITIALIZE_WALL_ARRAY')

CALL PROC_HVAC

! Write the Smokeview (.smv) file using parallel MPI writes

CALL WRITE_SMOKEVIEW_FILE

! Stop all the processes if this is just a set-up run

IF (SETUP_ONLY .OR. CHECK_MESH_ALIGNMENT .OR. STOP_STATUS/=0) THEN
   IF (MY_RANK==0) CALL INITIALIZE_DIAGNOSTIC_FILE(DT)
   IF (STOP_STATUS==0) STOP_STATUS = SETUP_ONLY_STOP
   CALL STOP_CHECK(1)
ENDIF

! Allocate and initialize MESH-specific variables

DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   CALL INITIALIZE_MESH_VARIABLES_1(DT,NM)
ENDDO

IF (MY_RANK==0 .AND. VERBOSE) CALL VERBOSE_PRINTOUT('Completed INITIALIZE_MESH_VARIABLES_1')

! MPI process 0 reopens the Smokeview file for additional output

IF (MY_RANK==0) THEN
   OPEN(LU_SMV,FILE=FN_SMV,FORM='FORMATTED', STATUS='OLD',POSITION='APPEND')
   CALL WRITE_STATUS_FILES
ENDIF

! Allocate and initialize OMESH arrays to hold "other mesh" data for a given mesh

DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   CALL INITIALIZE_MESH_EXCHANGE_1(NM)
ENDDO

IF (MY_RANK==0 .AND. VERBOSE) CALL VERBOSE_PRINTOUT('Completed INITIALIZE_MESH_EXCHANGE_1')

! Allocate "request" arrays to keep track of MPI communications

CALL MPI_INITIALIZATION_CHORES(4)

! Exchange information related to size of OMESH arrays

CALL MPI_INITIALIZATION_CHORES(5)

! Initial complex geometry CC setup

IF (CC_IBM) THEN
   CALL CC_SET_DATA(FIRST_CALL=.TRUE.) ! Define Cartesian cell types (used to define pressure zones), cut-cells, cfaces.
   CALL STOP_CHECK(1)
ENDIF

! Initialize PRESSURE_ZONEs

CALL INITIALIZE_PRESSURE_ZONES

CALL STOP_CHECK(1)

IF (MY_RANK==0 .AND. VERBOSE) CALL VERBOSE_PRINTOUT('Completed INITIALIZE_PRESSURE_ZONES')

! Allocate some arrays that are dimensioned with N_ZONE

CALL MPI_INITIALIZATION_CHORES(6)

! Allocate and initialize OMESH arrays to hold "other mesh" data for a given mesh

DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   CALL INITIALIZE_RADIATION_EXCHANGE(NM)
ENDDO

IF (MY_RANK==0 .AND. VERBOSE) CALL VERBOSE_PRINTOUT('Completed INITIALIZE_RADIATION_EXCHANGE')

! Initialize persistent MPI sends and receives and allocate buffer arrays.

CALL POST_RECEIVES(0)
CALL MESH_EXCHANGE(0)

CALL EXCHANGE_GEOMETRY_INFO

! Finish initializing mesh variables

DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   CALL INITIALIZE_MESH_VARIABLES_2(NM)
ENDDO

IF (MY_RANK==0 .AND. VERBOSE) CALL VERBOSE_PRINTOUT('Completed INITIALIZE_MESH_VARIABLES_2')

! Create arrays and communicators to exchange back WALL and THIN_WALL arrays across mesh boundaries.
! In the first call to the subroutine, all the WALL cells that are HT3D need to be exchanged, but once the 3-D noding is 
! done, there is no need to exchange all HT3D cells. The second call reduces the size of the exchange arrays.

CALL INITIALIZE_BACK_WALL_EXCHANGE(1)
CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
CALL INITIALIZE_BACK_WALL_EXCHANGE(2)

IF (MY_RANK==0 .AND. VERBOSE) CALL VERBOSE_PRINTOUT('Completed INITIALIZE_BACK_WALL_EXCHANGE')

CALL STOP_CHECK(1)

! Initialize turb arrays

DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   IF (TGA_SURF_INDEX>0) CYCLE
   CALL INIT_TURB_ARRAYS(NM)
ENDDO

! Final complex geometry CC setup

IF (CC_IBM) THEN
   CALL CC_SET_DATA(FIRST_CALL=.FALSE.) ! Interpolation Stencils, Scalar transport MATVEC data, cface RDNs.
   CALL STOP_CHECK(1)
ENDIF

! Initialize the flow field with random noise to eliminate false symmetries

DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   IF (TGA_SURF_INDEX>0) CYCLE
   CALL INITIAL_NOISE(NM)
   IF (PERIODIC_TEST==1) CALL NS_ANALYTICAL_SOLUTION(NM,T_BEGIN,RK_STAGE=2)
   IF (PERIODIC_TEST==2) CALL UVW_INIT(NM,UVW_FILE)
   IF (PERIODIC_TEST==3) CALL COMPRESSION_WAVE(NM,0._EB,3)
   IF (PERIODIC_TEST==4) CALL COMPRESSION_WAVE(NM,0._EB,4)
   IF (PERIODIC_TEST==6) CALL TWOD_VORTEX_CERFACS(NM)
   IF (PERIODIC_TEST==7) CALL SHUNN_MMS_3(DT,NM)
   IF (PERIODIC_TEST==8) CALL NS_ANALYTICAL_SOLUTION(NM,T_BEGIN,RK_STAGE=2)
   IF (PERIODIC_TEST==9) CALL SANDIA_DAT(NM,UVW_FILE)
   IF (PERIODIC_TEST==10) CALL TWOD_VORTEX_UMD(NM)
   IF (PERIODIC_TEST==11) CALL SAAD_MMS_1(NM)
   IF (PERIODIC_TEST==12) CALL TWOD_SOBOROT_UMD(NM)
   IF (PERIODIC_TEST==13) CALL TWOD_SOBOROT_UMD(NM)
   IF (PERIODIC_TEST==21) CALL ROTATED_CUBE_ANN_SOLN(NM,T_BEGIN) ! No Rotation.
   IF (PERIODIC_TEST==22) CALL ROTATED_CUBE_ANN_SOLN(NM,T_BEGIN) ! 27 deg Rotation.
   IF (PERIODIC_TEST==23) CALL ROTATED_CUBE_ANN_SOLN(NM,T_BEGIN) ! 45 deg Rotation.
   IF (UVW_RESTART)  CALL UVW_INIT(NM,CSVFINFO(NM)%UVWFILE)
   IF (TMP_RESTART)  CALL TMP_INIT(NM,CSVFINFO(NM)%TMPFILE)
   IF (SPEC_RESTART) CALL SPEC_INIT(NM,CSVFINFO(NM)%SPECFILE)
ENDDO

! Init centroid data (i.e. rho,zz) on cut-cells, cut-faces and CFACEs.
IF (CC_IBM) CALL INIT_CUTCELL_DATA(T_BEGIN,DT,FIRST_CALL=.TRUE.)

DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   IF (TGA_SURF_INDEX>0) CYCLE
   CALL COMPUTE_VISCOSITY(NM,APPLY_TO_ESTIMATED_VARIABLES=.FALSE.) ! needed here for KRES prior to mesh exchange
ENDDO
IF (MY_RANK==0 .AND. VERBOSE) CALL VERBOSE_PRINTOUT('Completed COMPUTE_VISCOSITY')

! Exchange information at mesh boundaries related to the various initialization routines just completed

IF (LEVEL_SET_MODE/=1) THEN
   CALL MESH_EXCHANGE(1)
   CALL MESH_EXCHANGE(4)
   CALL MESH_EXCHANGE(6)
ENDIF

! Ensure normal components of velocity match at mesh boundaries and do velocity BCs just in case the flow is not initialized to zero

PREDICTOR = .FALSE.
CORRECTOR = .TRUE.

DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   IF (TGA_SURF_INDEX>0) CYCLE
   CALL MATCH_VELOCITY(NM)
   CALL COMPUTE_VISCOSITY(NM,APPLY_TO_ESTIMATED_VARIABLES=.FALSE.) ! call again after mesh exchange
   IF (SYNTHETIC_EDDY_METHOD) CALL SYNTHETIC_EDDY_SETUP(NM)
   CALL VISCOSITY_BC(NM,APPLY_TO_ESTIMATED_VARIABLES=.FALSE.)
   CALL VELOCITY_BC(T_BEGIN,NM,APPLY_TO_ESTIMATED_VARIABLES=.FALSE.)
ENDDO

IF (MY_RANK==0 .AND. VERBOSE) CALL VERBOSE_PRINTOUT('Completed VELOCITY_BC match')

! Level Set model for firespread in vegetation

IF (LEVEL_SET_MODE>0 .OR. TERRAIN_CASE) THEN
   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL INITIALIZE_LEVEL_SET_FIRESPREAD_1(NM)
   ENDDO
   CALL MESH_EXCHANGE(14)
   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL INITIALIZE_LEVEL_SET_FIRESPREAD_2(NM,1)
   ENDDO
   CALL MESH_EXCHANGE(14)
   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL INITIALIZE_LEVEL_SET_FIRESPREAD_2(NM,2)
   ENDDO
ENDIF

! Check for CC_IBM initialization stop

IF (CC_IBM) THEN
   CALL STOP_CHECK(1)
   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL INITIALIZE_MESH_VARIABLES_3(NM)
   ENDDO
ENDIF

! Potentially read data from a previous calculation

IF (RESTART) THEN
   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL READ_RESTART(T,DT,NM)
   ENDDO
   IF (ABS(T-T_END)<TWO_EPSILON_EB) THEN
      STOP_STATUS = SETUP_STOP
      IF (MY_RANK==0) WRITE(LU_ERR,*) 'ERROR(1051): RESTART initial time equals T_END.'
   ENDIF
   IF (CC_IBM) CALL INIT_CUTCELL_DATA(T,DT,FIRST_CALL=.FALSE.)  ! Init centroid data (rho,zz) on cut-cells and cut-faces.
   CALL STOP_CHECK(1)
ENDIF

! Initialize output clocks

CALL INITIALIZE_OUTPUT_CLOCKS(T)
CALL STOP_CHECK(SETUP_STOP)

! Initialize output files that are mesh-specific

DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   IF (TGA_SURF_INDEX<1) CALL INITIALIZE_MESH_DUMPS(NM)
ENDDO

! Iterate surface BCs and radiation in case temperatures are not initialized to ambient

IF (.NOT.RESTART) THEN

   DO ITER=1,INITIAL_RADIATION_ITERATIONS
      DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
         IF (RADIATION) THEN
            CALL COMPUTE_RADIATION(T_BEGIN,NM,ITER)
            IF (CC_IBM) CALL CCCOMPUTE_RADIATION(T_BEGIN,NM,ITER)
         ENDIF
         CALL WALL_BC(T_BEGIN,DT,NM)
      ENDDO
      IF (RADIATION .AND. EXCHANGE_RADIATION) THEN
         DO ANG_INC_COUNTER=1,ANGLE_INCREMENT
            CALL MESH_EXCHANGE(2)
         ENDDO
      ENDIF
   ENDDO

   IF (MY_RANK==0 .AND. VERBOSE) CALL VERBOSE_PRINTOUT('Completed radiation initialization')
   IF (CC_IBM) CALL CC_RHO0W_INTERP

ENDIF

! Initialize particle distributions

CALL GENERATE_PARTICLE_DISTRIBUTIONS

! Initialize point devices, particles, profiles

DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   CALL INSERT_ALL_PARTICLES(T,NM)
   IF (TGA_SURF_INDEX<1) CALL INITIALIZE_DEVICES(NM)
   IF (TGA_SURF_INDEX<1) CALL INITIALIZE_PROFILES(NM)
ENDDO

IF (MY_RANK==0 .AND. VERBOSE) CALL VERBOSE_PRINTOUT('Completed particle initialization')

! Check for any stop flags at this point in the set up.

CALL STOP_CHECK(1)

! Check to see if only a TGA analysis is to be performed

TGA_IF: IF (TGA_SURF_INDEX>0) THEN
   CALL MPI_ALLREDUCE(MPI_IN_PLACE,TGA_MESH_INDEX,INTEGER_ONE,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,IERR)
   IF (TGA_MESH_INDEX>NMESHES) THEN
      IF (MY_RANK==0) WRITE(LU_ERR,'(3A)') 'ERROR(299): SURF ',TRIM(SURFACE(TGA_SURF_INDEX)%ID),' has no solid or particle for TGA.'
      STOP_STATUS = SETUP_STOP
   ELSE
      IF (MY_RANK==PROCESS(TGA_MESH_INDEX)) CALL TGA_ANALYSIS(TGA_MESH_INDEX)
      IF (MY_RANK==0) CALL INITIALIZE_DIAGNOSTIC_FILE(DT)
      STOP_STATUS = TGA_ANALYSIS_STOP
   ENDIF
   CALL STOP_CHECK(1)
ENDIF TGA_IF

! Initialize output files containing global data (Master Node Only)

IF (MY_RANK==0) THEN
   CALL INITIALIZE_GLOBAL_DUMPS(T,DT)
   IF (VERBOSE) CALL VERBOSE_PRINTOUT('Completed INITIALIZE_GLOBAL_DUMPS')
ENDIF

! Initialize HVAC variables

IF (HVAC_SOLVE) THEN
   IF (MY_RANK==0) CALL INIT_DUCT_NODE
   NODE_PROPERTIES = 0._EB
   NODE_ZONE       = 0
   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL HVAC_BC_IN(NM)
   ENDDO
   CALL STOP_CHECK(1)
   IF (N_MPI_PROCESSES>1) CALL EXCHANGE_HVAC_BC
   IF (MY_RANK==0) THEN
      CALL COLLAPSE_HVAC_BC(T)
      CALL SET_INIT_HVAC
   ENDIF
   CALL STOP_CHECK(1)
ENDIF

! Check for changes in VENT or OBSTruction control and device status at t=T_BEGIN

IF (.NOT.RESTART) CALL CREATE_OR_REMOVE_OBSTRUCTIONS

! Compute divergence just in case the flow field is not initialized to ambient

CALL INITIALIZE_DIVERGENCE_INTEGRALS

DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   CALL DIVERGENCE_PART_1(T_BEGIN,DT,NM)
ENDDO

IF (N_ZONE>0 .AND. .NOT.SOLID_PHASE_ONLY) CALL EXCHANGE_DIVERGENCE_INFO

CALL GLOBAL_MATRIX_REASSIGN(FORCE_REASSIGN=.TRUE.)

IF (MY_RANK==0 .AND. VERBOSE) CALL VERBOSE_PRINTOUT('Completed Poisson initialization')

! Make an initial dump of ambient values

IF (.NOT.RESTART) THEN
   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL UPDATE_GLOBAL_OUTPUTS(T,DT,NM)
      CALL DUMP_MESH_OUTPUTS(T,DT,NM)
   ENDDO
   IF (MY_RANK==0 .AND. VERBOSE) CALL VERBOSE_PRINTOUT('Completed DUMP_MESH_OUTPUTS')
ENDIF

! Make an initial dump of global output quantities

IF (.NOT.RESTART) THEN
   CALL EXCHANGE_GLOBAL_OUTPUTS
   CALL UPDATE_CONTROLS(T,0._EB,CTRL_STOP_STATUS,.TRUE.)
   CALL DUMP_GLOBAL_OUTPUTS
ENDIF

! Write out character strings to .smv file

CALL WRITE_STRINGS

! Sprinkler piping calculation

DO N=1,N_DEVC
   IF (DEVICE(N)%PROP_INDEX > 0 .AND.  DEVICE(N)%CURRENT_STATE) THEN
      IF (PROPERTY(DEVICE(N)%PROP_INDEX)%PART_INDEX > 0) DEVC_PIPE_OPERATING(DEVICE(N)%PIPE_INDEX) = &
         DEVC_PIPE_OPERATING(DEVICE(N)%PIPE_INDEX) + 1
   ENDIF
ENDDO

! Start the clock for time stepping

WALL_CLOCK_START_ITERATIONS = CURRENT_TIME()
T_USED = 0._EB
T_USED(1) = WALL_CLOCK_START_ITERATIONS

! This ends the initialization part of the program

INITIALIZATION_PHASE = .FALSE.

IF (UNFREEZE_TIME > 0._EB) THEN 
   FREEZE_VELOCITY=.TRUE. 
   SOLID_PHASE_ONLY=.TRUE.
   LOCK_TIME_STEP=.TRUE.
ENDIF

IF (MY_RANK==0 .AND. VERBOSE) CALL VERBOSE_PRINTOUT('Starting the time-stepping')

!***********************************************************************************************************************************
!                                                   MAIN TIMESTEPPING LOOP
!***********************************************************************************************************************************

MAIN_LOOP: DO

   ICYC  = ICYC + 1   ! Time step iterations

   ! Do not print out general diagnostics into .out file every time step

   DIAGNOSTICS = .FALSE.

   ! Check for program stops

   IF (MY_RANK==0) INQUIRE(FILE=FN_STOP,EXIST=EX)
   IF (N_MPI_PROCESSES>1) CALL MPI_BCAST(EX,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
   IF (EX .AND. ICYC>=STOP_AT_ITER) THEN
      IF (VERBOSE .AND. STOP_STATUS/=USER_STOP) WRITE(LU_ERR,'(A,I5)') ' STOP file detected, MPI Process =',MY_RANK
      STOP_STATUS = USER_STOP
      DIAGNOSTICS = .TRUE.
   ENDIF

   ! Check to see if the time step can be increased

   IF (ALL(CHANGE_TIME_STEP_INDEX==1)) DT = MINVAL(DT_NEW)

   ! Clip final time step

   IF ((T+DT+DT_END_FILL)>T_END) DT = MAX(T_END-T+TWO_EPSILON_EB,DT_END_MINIMUM)

   ! Determine when to dump out diagnostics to the .out file

   LO10 = INT(LOG10(REAL(MAX(1,ABS(ICYC)),EB)))
   IF (MOD(ICYC,10**LO10)==0 .OR. MOD(ICYC,DIAGNOSTICS_INTERVAL)==0 .OR. (T+DT)>=T_END) DIAGNOSTICS = .TRUE.

   IF ((UNFREEZE_TIME > 0._EB).AND.(T>UNFREEZE_TIME)) THEN 
      FREEZE_VELOCITY=.FALSE.
      SOLID_PHASE_ONLY=.FALSE.
      LOCK_TIME_STEP=.FALSE.
      
      UNFREEZE_DT_LOOP: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
         CALL CHECK_STABILITY(DT,DT_NEW,T,NM)
      ENDDO UNFREEZE_DT_LOOP
      IF (N_MPI_PROCESSES>1) THEN
         CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,DT_NEW(1:NMESHES),&
                             COUNTS,DISPLS,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERR)
      ENDIF
      DT = MINVAL(DT_NEW)
   ENDIF

   !================================================================================================================================
   !                                           Start of Predictor part of time step
   !================================================================================================================================

   PREDICTOR = .TRUE.
   CORRECTOR = .FALSE.

   ! Process externally controlled variables

   IF (READ_EXTERNAL) THEN
      IF (MY_RANK==0 .AND. T > T_EXTERNAL) THEN
         CALL READ_EXTERNAL_FILE(EXTERNAL_FAIL)
         IF (.NOT. EXTERNAL_FAIL) T_EXTERNAL = T + DT_EXTERNAL
      ENDIF
      IF (HEARTBEAT_FAIL) CALL STOP_CHECK(1)
      CALL EXCHANGE_EXTERNAL
   ENDIF

   ! Begin the finite differencing of the PREDICTOR step

   COMPUTE_FINITE_DIFFERENCES_1: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL INSERT_ALL_PARTICLES(T,NM)
      IF (.NOT.SOLID_PHASE_ONLY .AND. .NOT.FREEZE_VELOCITY) CALL COMPUTE_VISCOSITY(NM,APPLY_TO_ESTIMATED_VARIABLES=.FALSE.)
      CALL MASS_FINITE_DIFFERENCES_NEW(NM)
   ENDDO COMPUTE_FINITE_DIFFERENCES_1

   ! Estimate quantities at next time step, and decrease/increase time step if necessary based on CFL condition

   FIRST_PASS = .TRUE.

   CHANGE_TIME_STEP_LOOP: DO

      ! Predict species mass fractions at the next time step.

      COMPUTE_DENSITY_LOOP: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
         IF(CC_IBM .AND. .NOT.FIRST_PASS) CALL CC_RESTORE_UVW_UNLINKED(NM,ASSIGN_UNLINKED_VEL=.TRUE.)
         CALL DENSITY(T,DT,NM)
         IF (LEVEL_SET_MODE>0) CALL LEVEL_SET_FIRESPREAD(T,DT,NM)
      ENDDO COMPUTE_DENSITY_LOOP

      IF (LEVEL_SET_MODE==2 .AND. CHECK_FREEZE_VELOCITY) CALL CHECK_FREEZE_VELOCITY_STATUS

      IF (CC_IBM) CALL CC_DENSITY(T,DT)

      ! Exchange species mass fractions at interpolated boundaries.

      IF (LEVEL_SET_MODE/=1) CALL MESH_EXCHANGE(1)

      ! Exchange level set values, if necessary

      IF (LEVEL_SET_MODE>0) CALL MESH_EXCHANGE(14)

      ! Exchange newly inserted particles, if necessary

      IF (FIRST_PASS .AND. MPI_PARTICLE_EXCHANGE) THEN
         CALL MPI_ALLREDUCE(MPI_IN_PLACE,EXCHANGE_INSERTED_PARTICLES,INTEGER_ONE,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,IERR)
         IF (EXCHANGE_INSERTED_PARTICLES) THEN
            CALL MESH_EXCHANGE(7)
            CALL POST_RECEIVES(11)
            CALL MESH_EXCHANGE(11)
            EXCHANGE_INSERTED_PARTICLES = .FALSE.
         ENDIF
      ENDIF

      ! Calculate convective and diffusive terms of the velocity equation.

      IF (FIRST_PASS .AND. HVAC_SOLVE) THEN
         NODE_PROPERTIES = 0._EB
         NODE_ZONE       = 0
      ENDIF

      COMPUTE_DIVERGENCE_LOOP: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
         IF (.NOT.SOLID_PHASE_ONLY .AND. .NOT.FREEZE_VELOCITY) THEN
            MESHES(NM)%BAROCLINIC_TERMS_ATTACHED = .FALSE.
            CALL VISCOSITY_BC(NM,APPLY_TO_ESTIMATED_VARIABLES=.FALSE.)
            IF (.NOT.CYLINDRICAL) CALL VELOCITY_FLUX(T,DT,NM,APPLY_TO_ESTIMATED_VARIABLES=.FALSE.)
            IF (     CYLINDRICAL) CALL VELOCITY_FLUX_CYLINDRICAL(T,NM,APPLY_TO_ESTIMATED_VARIABLES=.FALSE.)
         ENDIF
         IF (FIRST_PASS .AND. HVAC_SOLVE) CALL HVAC_BC_IN(NM)
      ENDDO COMPUTE_DIVERGENCE_LOOP

      ! HVAC solver

      IF (HVAC_SOLVE) THEN
         IF (FIRST_PASS .AND. N_MPI_PROCESSES>1) CALL EXCHANGE_HVAC_BC
         IF (MY_RANK==0) CALL HVAC_CALC(T,DT,FIRST_PASS)
         IF (N_MPI_PROCESSES>1) CALL EXCHANGE_HVAC_SOLUTION
      ENDIF

      ! Boundary conditions for temperature, species, and density. Start divergence calculation.

      CALL INITIALIZE_DIVERGENCE_INTEGRALS

      COMPUTE_WALL_BC_LOOP_A: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
         CALL WALL_BC(T,DT,NM)
         IF (PARTICLE_DRAG) CALL PARTICLE_MOMENTUM_TRANSFER(DT,NM)
         CALL DIVERGENCE_PART_1(T,DT,NM)
      ENDDO COMPUTE_WALL_BC_LOOP_A

      ! If there are pressure ZONEs, exchange integrated quantities mesh to mesh for use in the divergence calculation

      IF (N_ZONE>0 .AND. .NOT.SOLID_PHASE_ONLY) CALL EXCHANGE_DIVERGENCE_INFO

      ! Finish the divergence calculation

      FINISH_DIVERGENCE_LOOP: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
         CALL DIVERGENCE_PART_2(DT,NM)
      ENDDO FINISH_DIVERGENCE_LOOP

      ! Solve for the pressure at the current time step

      IF (LEVEL_SET_MODE/=1) CALL PRESSURE_ITERATION_SCHEME

      ! Predict the velocity components at the next time step

      CHANGE_TIME_STEP_INDEX = 0
      DT_NEW = DT

      PREDICT_VELOCITY_LOOP: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
         CALL VELOCITY_PREDICTOR(T+DT,DT,DT_NEW,NM)
      ENDDO PREDICT_VELOCITY_LOOP

      ! Check if there is a numerical instability after updating the velocity field. If there is, exit this loop, finish the time
      ! step, and stop the code.

      CALL STOP_CHECK(0)

      IF (STOP_STATUS==INSTABILITY_STOP) THEN
         DIAGNOSTICS = .TRUE.
         SUPPRESS_DIAGNOSTICS = .FALSE.
         EXIT CHANGE_TIME_STEP_LOOP
      ENDIF

      ! Exchange CHANGE_TIME_STEP_INDEX to determine if the time step needs to be decreased (-1) or increased (1). If any mesh
      ! needs to decrease, or all need to increase, exchange the array of new time step values, DT_NEW.

      IF (N_MPI_PROCESSES>1) THEN
         TNOW = CURRENT_TIME()
         CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,CHANGE_TIME_STEP_INDEX(1:NMESHES),&
                             COUNTS,DISPLS,MPI_INTEGER,MPI_COMM_WORLD,IERR)
         IF (ANY(CHANGE_TIME_STEP_INDEX==-1) .OR. ALL(CHANGE_TIME_STEP_INDEX==1)) &
            CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,DT_NEW(1:NMESHES),&
                                COUNTS,DISPLS,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERR)
         T_USED(11) = T_USED(11) + CURRENT_TIME() - TNOW
      ENDIF

      IF (ANY(CHANGE_TIME_STEP_INDEX==-1)) THEN  ! If the time step was reduced, CYCLE CHANGE_TIME_STEP_LOOP
         DT = MINVAL(DT_NEW)
         FIRST_PASS = .FALSE.
      ELSE  ! exit the loop and if the time step is to be increased, this will occur at the next time step.
         EXIT CHANGE_TIME_STEP_LOOP
      ENDIF

   ENDDO CHANGE_TIME_STEP_LOOP

   ! If detailed CFL info needed

   IF (CFL_FILE) CALL WRITE_CFL_FILE

   ! Compute linked velocity arrays. Flux average final velocity to cutfaces.

   IF (CC_IBM) CALL CC_END_STEP(T,DT,DIAGNOSTICS)

   ! Exchange velocity and pressures at interpolated boundaries

   IF (LEVEL_SET_MODE/=1) CALL MESH_EXCHANGE(3)

   ! Force normal components of velocity to match at interpolated boundaries

   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL MATCH_VELOCITY(NM)
   ENDDO

   ! Apply tangential velocity boundary conditions

   VELOCITY_BC_LOOP: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      IF (SYNTHETIC_EDDY_METHOD) CALL SYNTHETIC_TURBULENCE(DT,T,NM)
      CALL VELOCITY_BC(T,NM,APPLY_TO_ESTIMATED_VARIABLES=.TRUE.)
   ENDDO VELOCITY_BC_LOOP

   !================================================================================================================================
   !                                           Start of Corrector part of time step
   !================================================================================================================================

   CORRECTOR = .TRUE.
   PREDICTOR = .FALSE.

   ! Advance the time to start the CORRECTOR step

   T = T + DT

   ! Zero out energy and mass balance arrays

   Q_DOT = 0._EB
   M_DOT = 0._EB

   ! Check for creation or removal of obsructions

   CALL CREATE_OR_REMOVE_OBSTRUCTIONS

   ! Finite differences for mass and momentum equations for the second half of the time step

   COMPUTE_FINITE_DIFFERENCES_2: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      IF (.NOT.SOLID_PHASE_ONLY .AND. .NOT.FREEZE_VELOCITY) CALL COMPUTE_VISCOSITY(NM,APPLY_TO_ESTIMATED_VARIABLES=.TRUE.)
      CALL MASS_FINITE_DIFFERENCES_NEW(NM)
      CALL DENSITY(T,DT,NM)
      IF (LEVEL_SET_MODE>0) CALL LEVEL_SET_FIRESPREAD(T,DT,NM)
   ENDDO COMPUTE_FINITE_DIFFERENCES_2

   IF (LEVEL_SET_MODE==2 .AND. CHECK_FREEZE_VELOCITY) CALL CHECK_FREEZE_VELOCITY_STATUS

   IF (CC_IBM) CALL CC_DENSITY(T,DT)

   ! Exchange species mass fractions.

   IF (LEVEL_SET_MODE/=1) CALL MESH_EXCHANGE(4)
   IF (LEVEL_SET_MODE>0) CALL MESH_EXCHANGE(14)

   ! Apply mass and species boundary conditions, update radiation, particles, and re-compute divergence

   COMPUTE_DIVERGENCE_2_1STPART: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      IF (.NOT.SOLID_PHASE_ONLY .AND. .NOT.FREEZE_VELOCITY) THEN
         MESHES(NM)%BAROCLINIC_TERMS_ATTACHED = .FALSE.
         CALL VISCOSITY_BC(NM,APPLY_TO_ESTIMATED_VARIABLES=.TRUE.)
         IF (.NOT.CYLINDRICAL) CALL VELOCITY_FLUX(T,DT,NM,APPLY_TO_ESTIMATED_VARIABLES=.TRUE.)
         IF (     CYLINDRICAL) CALL VELOCITY_FLUX_CYLINDRICAL(T,NM,APPLY_TO_ESTIMATED_VARIABLES=.TRUE.)
      ENDIF
      IF (N_AGGLOMERATION_SPECIES > 0) CALL CALC_AGGLOMERATION(DT,NM)
   ENDDO COMPUTE_DIVERGENCE_2_1STPART

   IF (N_REACTIONS > 0 .OR. INIT_HRRPUV) CALL COMBUSTION_LOAD_BALANCED(T,DT)

   COMPUTE_DIVERGENCE_2_2NDPART: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      IF (ANY(SPECIES_MIXTURE%CONDENSATION_SMIX_INDEX>0)) CALL CONDENSATION_EVAPORATION(DT,NM)
   ENDDO COMPUTE_DIVERGENCE_2_2NDPART


   IF (HVAC_SOLVE) CALL HVAC_CALC(T,DT,.TRUE.)

   ! Perform droplet mass/energy transfer, then move all particles, then transfer momentum to the gas.
   ! Solid particle mass/energy transfer is handled by WALL_BC.

   COMPUTE_WALL_BC_2A: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL PARTICLE_MASS_ENERGY_TRANSFER(T,DT,NM)
      CALL MOVE_PARTICLES(T,DT,NM)
      IF (PARTICLE_DRAG) CALL PARTICLE_MOMENTUM_TRANSFER(DT,NM)
   ENDDO COMPUTE_WALL_BC_2A

   ! Exchange the number of particles sent from mesh to mesh (7), and if non-zero, exchange particles (11)

   IF (MPI_PARTICLE_EXCHANGE) THEN
      CALL MESH_EXCHANGE(7)
      CALL POST_RECEIVES(11)
      CALL MESH_EXCHANGE(11)
   ENDIF

   ! Update the temperature, species and gas density boundary conditions at all surfaces

   WALL_COUNTER = WALL_COUNTER + 1
   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL WALL_BC(T,DT,NM)
   ENDDO
   IF (WALL_COUNTER==WALL_INCREMENT) WALL_COUNTER = 0

   IF (SOLID_HEAT_TRANSFER_3D .AND. WALL_COUNTER==0) THEN
      CALL MESH_EXCHANGE(6)
      DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
         CALL HT3D_TEMPERATURE_EXCHANGE(NM)
      ENDDO
   ENDIF

   ! Radiation

   DO ITER=1,RADIATION_ITERATIONS
      DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
         CALL COMPUTE_RADIATION(T,NM,ITER)
         IF (CC_IBM) CALL CCCOMPUTE_RADIATION(T,NM,ITER)
      ENDDO
      IF (RADIATION .AND. EXCHANGE_RADIATION .AND. RADIATION_ITERATIONS>1) THEN
         ! Only do an MPI exchange of radiation intensity if multiple iterations are requested.
         DO ANG_INC_COUNTER=1,ANGLE_INCREMENT
            CALL MESH_EXCHANGE(2)
            IF (ICYC>1) EXIT
         ENDDO
      ENDIF
   ENDDO

   ! Start the computation of the divergence term.

   CALL INITIALIZE_DIVERGENCE_INTEGRALS

   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      IF (N_REACTIONS>0 .OR. INIT_HRRPUV) CALL COMBUSTION_BC(NM)
      CALL DIVERGENCE_PART_1(T,DT,NM)
   ENDDO

   ! In most LES fire cases, a correction to the source term in the radiative transport equation is needed.

   IF (RTE_SOURCE_CORRECTION) CALL CALCULATE_RTE_SOURCE_CORRECTION_FACTOR

   ! Exchange global pressure zone information

   IF (N_ZONE>0 .AND. .NOT.SOLID_PHASE_ONLY) CALL EXCHANGE_DIVERGENCE_INFO

   ! Update global pressure matrices after zone connections

   CALL GLOBAL_MATRIX_REASSIGN(FORCE_REASSIGN=.FALSE.)

   ! Exchange mass loss information for OBSTs abutting interpolated boundaries

   IF (EXCHANGE_OBST_MASS) THEN
      CALL MESH_EXCHANGE(15)  ! Exchange number of OBSTs that have mass info to exchange across mesh boundaries
      CALL POST_RECEIVES(16)
      CALL MESH_EXCHANGE(16)  ! Satellite meshes send mass losses to meshes that actually contain the OBSTstruction
      CALL MESH_EXCHANGE(17)  ! Mesh containing the OBSTruction packs up its new mass to be sent to satellite meshes
      CALL POST_RECEIVES(18)
      CALL MESH_EXCHANGE(18)  ! Satellite meshes receive the new mass of OBSTructions that live in neighboring mesh
   ENDIF

   ! Finish computing the divergence

   FINISH_DIVERGENCE_LOOP_2: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL DIVERGENCE_PART_2(DT,NM)
   ENDDO FINISH_DIVERGENCE_LOOP_2

   ! Solve the pressure equation.

   IF (LEVEL_SET_MODE/=1) CALL PRESSURE_ITERATION_SCHEME

   ! Update the  velocity.

   CORRECT_VELOCITY_LOOP: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL VELOCITY_CORRECTOR(T,DT,NM)
      IF (DIAGNOSTICS .OR. STORE_CARTESIAN_DIVERGENCE) CALL CHECK_DIVERGENCE(NM)
   ENDDO CORRECT_VELOCITY_LOOP

   ! Compute linked velocity arrays. Flux average final velocity to cutfaces.

   IF (CC_IBM) CALL CC_END_STEP(T,DT,DIAGNOSTICS)

   ! Exchange velocity, pressure, particles at interpolated boundaries

   IF (LEVEL_SET_MODE/=1) CALL MESH_EXCHANGE(6)

   ! Exchange radiation intensity at interpolated boundaries if only one iteration of the solver is requested.

   IF (RADIATION .AND. EXCHANGE_RADIATION .AND. RADIATION_ITERATIONS==1) THEN
      DO ANG_INC_COUNTER=1,ANGLE_INCREMENT
         CALL MESH_EXCHANGE(2)
         IF (ICYC>1) EXIT
      ENDDO
   ENDIF

   ! Force normal components of velocity to match at interpolated boundaries

   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL MATCH_VELOCITY(NM)
   ENDDO

   ! Apply velocity boundary conditions, and update values of HRR, DEVC, etc.

   VELOCITY_BC_LOOP_2: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL VELOCITY_BC(T,NM,APPLY_TO_ESTIMATED_VARIABLES=.FALSE.)
      CALL UPDATE_GLOBAL_OUTPUTS(T,DT,NM)
   ENDDO VELOCITY_BC_LOOP_2

   ! Share device, HRR, mass data among all processes

   CALL EXCHANGE_GLOBAL_OUTPUTS

   ! Check for dumping end of timestep outputs

   CALL UPDATE_CONTROLS(T,DT,CTRL_STOP_STATUS,.FALSE.)
   IF (CTRL_STOP_STATUS) STOP_STATUS = CTRL_STOP

   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL DUMP_MESH_OUTPUTS(T,DT,NM)
   ENDDO

   ! Dump outputs such as HRR, DEVC, etc.

   CALL DUMP_GLOBAL_OUTPUTS

   ! Dump out diagnostics

   CALL WRITE_STRINGS
   IF (DIAGNOSTICS) THEN
      IF (.NOT.SUPPRESS_DIAGNOSTICS .AND. N_MPI_PROCESSES>1) CALL EXCHANGE_DIAGNOSTICS
      IF (MY_RANK==0) CALL WRITE_DIAGNOSTICS(T,DT)
   ENDIF

   ! Flush output file buffers

   IF (FLUSH_FILE_BUFFERS) THEN
      IF (T>=FLSH_CLOCK(FLSH_COUNTER(1))) THEN
         IF (MY_RANK==0) CALL FLUSH_GLOBAL_BUFFERS
         FLSH_COUNTER(1) = FLSH_COUNTER(1) + 1
      ENDIF
   ENDIF

   ! Dump a restart file if necessary

   IF ( T>=RSRT_CLOCK(RSRT_COUNTER(1)) .OR. STOP_STATUS==USER_STOP) THEN
      CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,RADIATION_COMPLETED(1:NMESHES),&
                          COUNTS,DISPLS,MPI_LOGICAL,MPI_COMM_WORLD,IERR)
      IF (T>=T_END .OR. ALL(RADIATION_COMPLETED) .OR. STOP_STATUS==CTRL_STOP) THEN
         DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
            CALL DUMP_RESTART(T,DT,NM)
         ENDDO
         RSRT_COUNTER(1) = RSRT_COUNTER(1) + 1
      ENDIF
   ENDIF

   ! Check for abnormal run stop

   CALL STOP_CHECK(1)  ! The argument 1 means that FDS will end unless there is logic associated with the STOP_STATUS

   ! Stop the run normally

   IF (MY_RANK==0 .AND. VERBOSE) CALL VERBOSE_PRINTOUT('End of time step')

   IF (T>=T_END .AND. ICYC>0) EXIT MAIN_LOOP

ENDDO MAIN_LOOP

!***********************************************************************************************************************************
!                                                     END OF TIME STEPPING LOOP
!***********************************************************************************************************************************

! Deallocate GLMAT_SOLVER_H variables if needed:

SELECT CASE(PRES_FLAG)
   CASE(ULMAT_FLAG)
      DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
         CALL FINISH_ULMAT_SOLVER(NM)
      ENDDO
   CASE(GLMAT_FLAG,UGLMAT_FLAG); CALL FINISH_GLMAT_SOLVER
END SELECT

! Finish unstructured geometry

IF (CC_IBM) CALL FINISH_CC

! Stop the calculation

CALL END_FDS

! This is the end of program. Supporting routines are listed below.

CONTAINS


!> \brief Check the MPI threading support level

SUBROUTINE CHECK_MPI

IF (PROVIDED<REQUIRED) THEN
   IF (MY_RANK==0) WRITE(LU_ERR,'(A)') ' WARNING:  This MPI implementation provides insufficient threading support.'
   !$ CALL OMP_SET_NUM_THREADS(1)
ENDIF

END SUBROUTINE CHECK_MPI


!> \brief A collection of miscellaneous initializations
!> \param TASK_NUMBER Integer denoting which sets of tasks to do

SUBROUTINE MPI_INITIALIZATION_CHORES(TASK_NUMBER)

INTEGER, INTENT(IN) :: TASK_NUMBER
TYPE (MPI_REQUEST), ALLOCATABLE, DIMENSION(:) :: REQ0,REQ0DUM
TYPE (MPI_GROUP) :: GROUP_WORLD,SUBGROUP
INTEGER :: N_REQ0,SNODE,MEMBERS(0:NMESHES-1),NN,NOM,N_COMMUNICATIONS,NZ
CHARACTER(50) :: DUMMY_STRING

SELECT CASE(TASK_NUMBER)

   CASE(1)

      ! Set up send and receive buffer counts and displacements

      ALLOCATE(REAL_BUFFER_DUCT((2+N_TRACKED_SPECIES)*N_DUCTNODES+N_DUCTS))
      ALLOCATE(REAL_BUFFER_10(10,NMESHES))
      ALLOCATE(REAL_BUFFER_20(20,NMESHES))

      IF (READ_EXTERNAL) THEN
         ALLOCATE(REAL_BUFFER_EXTERNAL(N_RAMP))
         ALLOCATE(LOGICAL_BUFFER_EXTERNAL(N_CTRL))
      ENDIF

      ALLOCATE(COUNTS(0:N_MPI_PROCESSES-1))
      ALLOCATE(COUNTS_10(0:N_MPI_PROCESSES-1))
      ALLOCATE(COUNTS_20(0:N_MPI_PROCESSES-1))
      ALLOCATE(COUNTS_TP(0:N_MPI_PROCESSES-1))

      ALLOCATE(DISPLS(0:N_MPI_PROCESSES-1))
      ALLOCATE(DISPLS_10(0:N_MPI_PROCESSES-1))
      ALLOCATE(DISPLS_20(0:N_MPI_PROCESSES-1))
      ALLOCATE(DISPLS_TP(0:N_MPI_PROCESSES-1))
      ALLOCATE(I_OFFSET(NMESHES))

      COUNTS    = 0
      COUNTS_TP = 0
      DO N=0,N_MPI_PROCESSES-1
         DO NM=1,NMESHES
            IF (PROCESS(NM)==N) THEN
               COUNTS(N)    = COUNTS(N)    + 1
               COUNTS_TP(N) = COUNTS_TP(N) + MESHES(NM)%IBAR
            ENDIF
         ENDDO
      ENDDO

      I_OFFSET(1) = 0
      DO NM=2,NMESHES
         I_OFFSET(NM) = I_OFFSET(NM-1) + MESHES(NM-1)%IBAR  ! Used for TUNNEL_PRECONDITIONER
      ENDDO

      DISPLS(0)    = 0
      DISPLS_TP(0) = 0
      DO N=1,N_MPI_PROCESSES-1
         DISPLS(N)    = COUNTS(N-1)    + DISPLS(N-1)
         DISPLS_TP(N) = COUNTS_TP(N-1) + DISPLS_TP(N-1)
      ENDDO
      COUNTS_10 = COUNTS*10
      DISPLS_10 = DISPLS*10
      COUNTS_20 = COUNTS*20
      DISPLS_20 = DISPLS*20

   CASE(2)

      ! Allocate TIME arrays

      ALLOCATE(DT_NEW(NMESHES),STAT=IZERO) ;  CALL ChkMemErr('MAIN','DT_NEW',IZERO) ; DT_NEW = DT

      ! Set up dummy arrays to hold various arrays that must be exchanged among meshes

      ALLOCATE(STATE_ARRAY(2*N_DEVC),STAT=IZERO) ; CALL ChkMemErr('MAIN','STATE_ARRAY',IZERO)
      ALLOCATE(TC_ARRAY(4*N_DEVC),STAT=IZERO)    ; CALL ChkMemErr('MAIN','TC_ARRAY',IZERO)
      ALLOCATE(TC2_ARRAY(2,N_DEVC),STAT=IZERO)   ; CALL ChkMemErr('MAIN','TC2_ARRAY',IZERO)


   CASE(3)

      ! Allocate arrays that hold geometric information of neighboring meshes.

      DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
         M => MESHES(NM)
         DO N=1,M%N_NEIGHBORING_MESHES
            NOM = M%NEIGHBORING_MESH(N)
            SNODE = PROCESS(NOM)
            IF (SNODE==MY_RANK) CYCLE
            M4 => MESHES(NOM)
            IF (.NOT.ALLOCATED(M4%PRESSURE_ZONE)) THEN
               ALLOCATE(M4%PRESSURE_ZONE(0:M4%IBP1,0:M4%JBP1,0:M4%KBP1))
               M4%PRESSURE_ZONE = -1
            ENDIF
         ENDDO
      ENDDO

      ! Create MPI communicators (sub-groups of MPI_COMM_WORLD) for each mesh
      ! that consists of the MPI processes that control it and its neighbors
      ! MPI_COMM_NEIGHBORS(NM) is the name of the communicator for mesh NM
      ! MPI_COMM_NEIGHBORS_ROOT(NM) is the rank within MPI_COMM_NEIGHBORS(NM) of the process that controls mesh NM

      ALLOCATE(MPI_COMM_NEIGHBORS(NMESHES))
      ALLOCATE(MPI_COMM_NEIGHBORS_ROOT(NMESHES))

      CALL MPI_COMM_GROUP(MPI_COMM_WORLD,GROUP_WORLD,IERR)  ! Get the group handle for MPI_COMM_WORLD

      DO NM=1,NMESHES
         M => MESHES(NM)
         NN = -1
         DO N=1,M%N_NEIGHBORING_MESHES
            NOM = M%NEIGHBORING_MESH(N)
            NN = NN + 1
            IF (NN==0) THEN
               MEMBERS(NN) = PROCESS(NOM)
            ELSEIF (NN>0 .AND. MEMBERS(NN-1)/=PROCESS(NOM)) THEN
               MEMBERS(NN) = PROCESS(NOM)  ! The ranks within MPI_COMM_WORLD to include in the new communicator; i.e. the neighbors
            ELSE
               NN = NN - 1
            ENDIF
            IF (NM==NOM) MPI_COMM_NEIGHBORS_ROOT(NM) = NN
         ENDDO
         CALL MPI_GROUP_INCL(GROUP_WORLD,NN+1,MEMBERS(0:NN),SUBGROUP,IERR)  ! Create the new sub-group of GROUP_WORLD
         CALL MPI_COMM_CREATE(MPI_COMM_WORLD,SUBGROUP,MPI_COMM_NEIGHBORS(NM),IERR) ! Create the new communicator
      ENDDO


   CASE(4)

      ! Allocate "request" arrays to keep track of MPI communications

      N_COMMUNICATIONS = 0
      DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
         M => MESHES(NM)
         N_COMMUNICATIONS = N_COMMUNICATIONS + M%N_NEIGHBORING_MESHES
      ENDDO

      ALLOCATE(REQ(N_COMMUNICATIONS*40))
      ALLOCATE(REQ1(N_COMMUNICATIONS*4))
      ALLOCATE(REQ2(N_COMMUNICATIONS*12))
      ALLOCATE(REQ3(N_COMMUNICATIONS*4))
      ALLOCATE(REQ4(N_COMMUNICATIONS*4))
      ALLOCATE(REQ5(N_COMMUNICATIONS*4))
      ALLOCATE(REQ6(N_COMMUNICATIONS*12))
      ALLOCATE(REQ7(N_COMMUNICATIONS*4))
      ALLOCATE(REQ14(N_COMMUNICATIONS*4))
      ALLOCATE(REQ15(N_COMMUNICATIONS*4))

      REQ = MPI_REQUEST_NULL
      REQ1 = MPI_REQUEST_NULL
      REQ2 = MPI_REQUEST_NULL
      REQ3 = MPI_REQUEST_NULL
      REQ4 = MPI_REQUEST_NULL
      REQ5 = MPI_REQUEST_NULL
      REQ6 = MPI_REQUEST_NULL
      REQ7 = MPI_REQUEST_NULL
      REQ14 = MPI_REQUEST_NULL
      REQ15 = MPI_REQUEST_NULL

   CASE(5)

      IF (N_MPI_PROCESSES>1) THEN

         ALLOCATE(REQ0(NMESHES)) ; N_REQ0 = 0

         DO NM=1,NMESHES
            DO NOM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
               IF (PROCESS(NM)/=MY_RANK .AND. MESHES(NOM)%CONNECTED_MESH(NM)) THEN
                  M2 => MESHES(NOM)%OMESH(NM)
                  N_REQ0 = N_REQ0 + 1
                  IF(N_REQ0>SIZE(REQ0,DIM=1)) THEN
                     ALLOCATE(REQ0DUM(SIZE(REQ0,DIM=1)+NMESHES))
                     REQ0DUM(1:N_REQ0-1) = REQ0(1:N_REQ0-1)
                     CALL MOVE_ALLOC(REQ0DUM,REQ0)
                  ENDIF
                  CALL MPI_IRECV(M2%INTEGER_RECV_BUFFER(1),7,MPI_INTEGER,PROCESS(NM),NM,MPI_COMM_WORLD,REQ0(N_REQ0),IERR)
               ENDIF
            ENDDO
         ENDDO

      ENDIF


      ! DEFINITION NIC_S:   MESHES(NOM)%OMESH(NM)%NIC_S   = MESHES(NM)%OMESH(NOM)%NIC_R
      ! DEFINITION I_MIN_S: MESHES(NOM)%OMESH(NM)%I_MIN_S = MESHES(NM)%OMESH(NOM)%I_MIN_R
      ! DEFINITION I_MAX_S: MESHES(NOM)%OMESH(NM)%I_MAX_S = MESHES(NM)%OMESH(NOM)%I_MAX_R
      ! DEFINITION J_MIN_S: MESHES(NOM)%OMESH(NM)%J_MIN_S = MESHES(NM)%OMESH(NOM)%J_MIN_R
      ! DEFINITION J_MAX_S: MESHES(NOM)%OMESH(NM)%J_MAX_S = MESHES(NM)%OMESH(NOM)%J_MAX_R
      ! DEFINITION K_MIN_S: MESHES(NOM)%OMESH(NM)%K_MIN_S = MESHES(NM)%OMESH(NOM)%K_MIN_R
      ! DEFINITION K_MAX_S: MESHES(NOM)%OMESH(NM)%K_MAX_S = MESHES(NM)%OMESH(NOM)%K_MAX_R

      DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
         DO NOM=1,NMESHES
            IF (.NOT.MESHES(NM)%CONNECTED_MESH(NOM)) CYCLE
            M3 => MESHES(NM)%OMESH(NOM)
            IF (PROCESS(NOM)/=MY_RANK) THEN
               M3%INTEGER_SEND_BUFFER(1) = M3%I_MIN_R
               M3%INTEGER_SEND_BUFFER(2) = M3%I_MAX_R
               M3%INTEGER_SEND_BUFFER(3) = M3%J_MIN_R
               M3%INTEGER_SEND_BUFFER(4) = M3%J_MAX_R
               M3%INTEGER_SEND_BUFFER(5) = M3%K_MIN_R
               M3%INTEGER_SEND_BUFFER(6) = M3%K_MAX_R
               M3%INTEGER_SEND_BUFFER(7) = M3%NIC_R
               N_REQ0 = N_REQ0 + 1
               IF(N_REQ0>SIZE(REQ0,DIM=1)) THEN
                  ALLOCATE(REQ0DUM(SIZE(REQ0,DIM=1)+NMESHES))
                  REQ0DUM(1:N_REQ0-1) = REQ0(1:N_REQ0-1)
                  CALL MOVE_ALLOC(REQ0DUM,REQ0)
               ENDIF
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

      IF (N_MPI_PROCESSES>1) THEN

         CALL TIMEOUT('Mesh connection test',N_REQ0,REQ0(1:N_REQ0))

         DO NM=1,NMESHES
            DO NOM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
               IF (NM/=NOM .AND. PROCESS(NM)/=MY_RANK .AND. MESHES(NOM)%CONNECTED_MESH(NM)) THEN
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

         DEALLOCATE(REQ0)

      ENDIF

      ! Exchange IIO_S, etc., the indices of interpolated cells
      ! DEFINITION IIO_S: MESHES(NOM)%OMESH(NM)%IIO_S = MESHES(NM)%OMESH(NOM)%IIO_R
      ! DEFINITION JJO_S: MESHES(NOM)%OMESH(NM)%JJO_S = MESHES(NM)%OMESH(NOM)%JJO_R
      ! DEFINITION KKO_S: MESHES(NOM)%OMESH(NM)%KKO_S = MESHES(NM)%OMESH(NOM)%KKO_R
      ! DEFINITION IOR_S: MESHES(NOM)%OMESH(NM)%IOR_S = MESHES(NM)%OMESH(NOM)%IOR_R

      DO NM=1,NMESHES
         DO NOM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
            IF (MESHES(NOM)%OMESH(NM)%NIC_S>0) THEN
               M2 => MESHES(NOM)%OMESH(NM)
               ALLOCATE(M2%IIO_S(M2%NIC_S))
               ALLOCATE(M2%JJO_S(M2%NIC_S))
               ALLOCATE(M2%KKO_S(M2%NIC_S))
               ALLOCATE(M2%IOR_S(M2%NIC_S))
            ENDIF
         ENDDO
      ENDDO

      N_REQ = 0

      DO NM=1,NMESHES
         DO NOM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
            M2 => MESHES(NOM)%OMESH(NM)
            IF (N_MPI_PROCESSES>1 .AND. NM/=NOM .AND. PROCESS(NM)/=MY_RANK .AND. M2%NIC_S>0) THEN
               CALL MPI_IRECV(M2%IIO_S(1),M2%NIC_S,MPI_INTEGER,PROCESS(NM),NM,MPI_COMM_WORLD,REQ(N_REQ+1),IERR)
               CALL MPI_IRECV(M2%JJO_S(1),M2%NIC_S,MPI_INTEGER,PROCESS(NM),NM,MPI_COMM_WORLD,REQ(N_REQ+2),IERR)
               CALL MPI_IRECV(M2%KKO_S(1),M2%NIC_S,MPI_INTEGER,PROCESS(NM),NM,MPI_COMM_WORLD,REQ(N_REQ+3),IERR)
               CALL MPI_IRECV(M2%IOR_S(1),M2%NIC_S,MPI_INTEGER,PROCESS(NM),NM,MPI_COMM_WORLD,REQ(N_REQ+4),IERR)
               N_REQ = N_REQ + 4
            ENDIF
         ENDDO
      ENDDO

      DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
         DO NOM=1,NMESHES
            M3 => MESHES(NM)%OMESH(NOM)
            IF (M3%NIC_R<1) CYCLE
            IF (PROCESS(NOM)/=MY_RANK) THEN
               CALL MPI_ISEND(M3%IIO_R(1),M3%NIC_R,MPI_INTEGER,PROCESS(NOM),NM,MPI_COMM_WORLD,REQ(N_REQ+1),IERR)
               CALL MPI_ISEND(M3%JJO_R(1),M3%NIC_R,MPI_INTEGER,PROCESS(NOM),NM,MPI_COMM_WORLD,REQ(N_REQ+2),IERR)
               CALL MPI_ISEND(M3%KKO_R(1),M3%NIC_R,MPI_INTEGER,PROCESS(NOM),NM,MPI_COMM_WORLD,REQ(N_REQ+3),IERR)
               CALL MPI_ISEND(M3%IOR_R(1),M3%NIC_R,MPI_INTEGER,PROCESS(NOM),NM,MPI_COMM_WORLD,REQ(N_REQ+4),IERR)
               N_REQ = N_REQ + 4
            ELSE
               M2 => MESHES(NOM)%OMESH(NM)
               M2%IIO_S = M3%IIO_R
               M2%JJO_S = M3%JJO_R
               M2%KKO_S = M3%KKO_R
               M2%IOR_S = M3%IOR_R
            ENDIF
         ENDDO
      ENDDO

      IF (N_REQ>0 .AND. N_MPI_PROCESSES>1) THEN
         CALL MPI_WAITALL(N_REQ,REQ(1:N_REQ),MPI_STATUSES_IGNORE,IERR)
         N_REQ = 0
      ENDIF

   CASE(6)

      ! CONNECTED_ZONES is a matrix of 0's and 1's such that if two zones are connected, then
      ! CONNECTED_ZONES(NZ1,NZ2)=CONNECTED_ZONES(NZ2,NZ1)=1. The diagonal is always 1, as in a zone is always connected to itself.

      ALLOCATE(CONNECTED_ZONES(0:N_ZONE,0:N_ZONE),STAT=IZERO) ; CALL ChkMemErr('INIT','CONNECTED_ZONES',IZERO) 
      CONNECTED_ZONES = 0
      DO NZ=0,N_ZONE
         CONNECTED_ZONES(NZ,NZ) = 1
      ENDDO

      ! DSUM, PSUM, USUM are summations of different parts of the divergence expression

      ALLOCATE(DSUM(N_ZONE),STAT=IZERO) ; CALL ChkMemErr('MAIN','DSUM',IZERO) ; DSUM = 0._EB
      ALLOCATE(PSUM(N_ZONE),STAT=IZERO) ; CALL ChkMemErr('MAIN','PSUM',IZERO) ; PSUM = 0._EB
      ALLOCATE(USUM(N_ZONE),STAT=IZERO) ; CALL ChkMemErr('MAIN','USUM',IZERO) ; USUM = 0._EB

      ! Determine if consumable OBST masses are to be exchanged

      CALL MPI_ALLREDUCE(MPI_IN_PLACE,EXCHANGE_OBST_MASS,INTEGER_ONE,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,IERR)

END SELECT

IF (MY_RANK==0 .AND. VERBOSE) THEN
   WRITE(DUMMY_STRING,'(A,I0)') 'Completed Initialization Step ',TASK_NUMBER
   CALL VERBOSE_PRINTOUT(DUMMY_STRING)
ENDIF

END SUBROUTINE MPI_INITIALIZATION_CHORES


!> \brief Perform multiple pressure solves until velocity tolerance is satisfied

SUBROUTINE PRESSURE_ITERATION_SCHEME

USE CC_SCALARS, ONLY : GET_LINKED_FV
INTEGER :: NM_MAX_V,NM_MAX_P
REAL(EB) :: TNOW,VELOCITY_ERROR_MAX_OLD,PRESSURE_ERROR_MAX_OLD

PRESSURE_ITERATIONS = 0

IF (BAROCLINIC) THEN
   ITERATE_BAROCLINIC_TERM = .TRUE.
ELSE
   ITERATE_BAROCLINIC_TERM = .FALSE.
ENDIF

IF(CC_IBM) THEN
   ! Here we need an exchange of F for linking:
   CALL MESH_EXCHANGE(5)
   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL GET_LINKED_FV(NM,DO_BAROCLINIC=.FALSE.)
   ENDDO
ENDIF

PRESSURE_ITERATION_LOOP: DO

   PRESSURE_ITERATIONS = PRESSURE_ITERATIONS + 1
   TOTAL_PRESSURE_ITERATIONS = TOTAL_PRESSURE_ITERATIONS + 1

   ! The following loops and exchange always get executed the first pass through the PRESSURE_ITERATION_LOOP.
   ! If we need to iterate the baroclinic torque term, the loop is executed each time.

   IF (ITERATE_BAROCLINIC_TERM .OR. PRESSURE_ITERATIONS==1) THEN
      DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
         IF (BAROCLINIC) CALL BAROCLINIC_CORRECTION(T,NM)
         IF (CC_IBM) CALL CC_NO_FLUX(DT,NM,.TRUE.)
      ENDDO
      CALL MESH_EXCHANGE(5)  ! Exchange FVX, FVY, FVZ
      DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
         CALL MATCH_VELOCITY_FLUX(NM)
      ENDDO
   ENDIF

   ! Compute the right hand side (RHS) and boundary conditions for the Poission equation for pressure.
   ! The WALL_WORK1 array is computed in COMPUTE_VELOCITY_ERROR, but it should
   ! be zero the first time the pressure solver is called.

   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL NO_FLUX(DT,NM)
      IF (CC_IBM) CALL CC_NO_FLUX(DT,NM,.FALSE.) ! set WALL_WORK1 to 0 in cells inside geometries.
      IF (PRESSURE_ITERATIONS==1) MESHES(NM)%WALL_WORK1 = 0._EB
      CALL PRESSURE_SOLVER_COMPUTE_RHS(T,DT,NM)
   ENDDO

   ! Solve the Poission equation using either FFT or ULMAT, GLMAT, or UGLMAT

   SELECT CASE(PRES_FLAG)
      CASE (FFT_FLAG)
         IF (TUNNEL_PRECONDITIONER) CALL TUNNEL_POISSON_SOLVER
         DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
            CALL PRESSURE_SOLVER_FFT(NM)
         ENDDO
      CASE (GLMAT_FLAG,UGLMAT_FLAG)
         CALL GLMAT_SOLVER(T,DT)
         CALL MESH_EXCHANGE(5)
         CALL COPY_H_OMESH_TO_MESH
      CASE (ULMAT_FLAG)
         IF (TUNNEL_PRECONDITIONER) CALL TUNNEL_POISSON_SOLVER
         DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
            CALL ULMAT_SOLVER(NM,T,DT)
         ENDDO
   END SELECT

   ! Check the residuals of the Poisson solution

   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      SELECT CASE(PRES_FLAG)
         CASE DEFAULT
            CALL PRESSURE_SOLVER_CHECK_RESIDUALS(NM)
         CASE (UGLMAT_FLAG,ULMAT_FLAG)
            CALL PRESSURE_SOLVER_CHECK_RESIDUALS_U(NM)
      END SELECT
   ENDDO

   IF (.NOT.ITERATE_PRESSURE) EXIT PRESSURE_ITERATION_LOOP

   ! Exchange both H or HS and FVX, FVY, FVZ and then estimate values of U, V, W (US, VS, WS) at next time step.

   CALL MESH_EXCHANGE(5)

   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL COMPUTE_VELOCITY_ERROR(DT,NM)
      IF (CC_IBM) CALL CC_COMPUTE_VELOCITY_ERROR(DT,NM) ! Inside solids respect to zero velocity.
   ENDDO

   ! Make all MPI processes aware of the maximum velocity error to decide if another pressure iteration is needed.

   IF (N_MPI_PROCESSES>1) THEN
      TNOW = CURRENT_TIME()
      REAL_BUFFER_10(  1,LOWER_MESH_INDEX:UPPER_MESH_INDEX) = VELOCITY_ERROR_MAX(LOWER_MESH_INDEX:UPPER_MESH_INDEX)
      REAL_BUFFER_10(  2,LOWER_MESH_INDEX:UPPER_MESH_INDEX) = PRESSURE_ERROR_MAX(LOWER_MESH_INDEX:UPPER_MESH_INDEX)
      REAL_BUFFER_10(3:5,LOWER_MESH_INDEX:UPPER_MESH_INDEX) = VELOCITY_ERROR_MAX_LOC(1:3,LOWER_MESH_INDEX:UPPER_MESH_INDEX)
      REAL_BUFFER_10(6:8,LOWER_MESH_INDEX:UPPER_MESH_INDEX) = PRESSURE_ERROR_MAX_LOC(1:3,LOWER_MESH_INDEX:UPPER_MESH_INDEX)
      CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,REAL_BUFFER_10(1:10,1:NMESHES),&
                          COUNTS_10,DISPLS_10,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERR)
      VELOCITY_ERROR_MAX(:)         =     REAL_BUFFER_10(1,:)
      PRESSURE_ERROR_MAX(:)         =     REAL_BUFFER_10(2,:)
      VELOCITY_ERROR_MAX_LOC(1:3,:) = INT(REAL_BUFFER_10(3:5,:))
      PRESSURE_ERROR_MAX_LOC(1:3,:) = INT(REAL_BUFFER_10(6:8,:))
      T_USED(11)=T_USED(11) + CURRENT_TIME() - TNOW
   ENDIF

   IF (MY_RANK==0 .AND. VELOCITY_ERROR_FILE) THEN
      NM_MAX_V = MAXLOC(VELOCITY_ERROR_MAX,DIM=1)
      NM_MAX_P = MAXLOC(PRESSURE_ERROR_MAX,DIM=1)
      WRITE(LU_VELOCITY_ERROR,'(E16.8,A,7(I7,A),E16.8,A,4(I7,A),E16.8)') T,',',ICYC,',',PRESSURE_ITERATIONS,',',&
         TOTAL_PRESSURE_ITERATIONS,',',&
         NM_MAX_V,',',VELOCITY_ERROR_MAX_LOC(1,NM_MAX_V),',',VELOCITY_ERROR_MAX_LOC(2,NM_MAX_V),',',&
         VELOCITY_ERROR_MAX_LOC(3,NM_MAX_V),',',MAXVAL(VELOCITY_ERROR_MAX),',',&
         NM_MAX_P,',',PRESSURE_ERROR_MAX_LOC(1,NM_MAX_P),',',PRESSURE_ERROR_MAX_LOC(2,NM_MAX_P),',',&
         PRESSURE_ERROR_MAX_LOC(3,NM_MAX_P),',',MAXVAL(PRESSURE_ERROR_MAX)
   ENDIF

   ! If the VELOCITY_TOLERANCE is satisfied or max/min iterations are hit, exit the loop.

   IF (MAXVAL(PRESSURE_ERROR_MAX)<PRESSURE_TOLERANCE) ITERATE_BAROCLINIC_TERM = .FALSE.

   IF (PREDICTOR .AND. PRESSURE_ITERATIONS>=MAX_PREDICTOR_PRESSURE_ITERATIONS) EXIT PRESSURE_ITERATION_LOOP
   IF (CORRECTOR .AND. PRESSURE_ITERATIONS>=MAX_PRESSURE_ITERATIONS)           EXIT PRESSURE_ITERATION_LOOP

   IF (MAXVAL(PRESSURE_ERROR_MAX)<PRESSURE_TOLERANCE .AND. MAXVAL(VELOCITY_ERROR_MAX)<VELOCITY_TOLERANCE) THEN
      EXIT PRESSURE_ITERATION_LOOP
   ENDIF

   ! Exit the iteration loop if satisfactory progress is not achieved

   IF (SUSPEND_PRESSURE_ITERATIONS .AND. ICYC>10) THEN
      IF (PRESSURE_ITERATIONS>3 .AND.  &
         MAXVAL(VELOCITY_ERROR_MAX)>ITERATION_SUSPEND_FACTOR*VELOCITY_ERROR_MAX_OLD .AND. &
         MAXVAL(PRESSURE_ERROR_MAX)>ITERATION_SUSPEND_FACTOR*PRESSURE_ERROR_MAX_OLD) EXIT PRESSURE_ITERATION_LOOP
      VELOCITY_ERROR_MAX_OLD = MAXVAL(VELOCITY_ERROR_MAX)
      PRESSURE_ERROR_MAX_OLD = MAXVAL(PRESSURE_ERROR_MAX)
   ENDIF

ENDDO PRESSURE_ITERATION_LOOP

END SUBROUTINE PRESSURE_ITERATION_SCHEME


!> \brief Compute a running average of the source correction factor for the radiative transport scheme.

SUBROUTINE CALCULATE_RTE_SOURCE_CORRECTION_FACTOR

REAL(EB), PARAMETER :: WGT=0.9_EB
REAL(EB) :: TNOW

TNOW = CURRENT_TIME()

! Sum up the components of the corrective factor from all the meshes.

IF (N_MPI_PROCESSES>1) THEN
   CALL MPI_ALLREDUCE(MPI_IN_PLACE,RAD_Q_SUM,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
   CALL MPI_ALLREDUCE(MPI_IN_PLACE,KFST4_SUM,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
ENDIF

! Compute the corrective factor for the RTE. Note that the max value of 100 is arbitrary.

IF (KFST4_SUM>TWO_EPSILON_EB) &
   RTE_SOURCE_CORRECTION_FACTOR = WGT*RTE_SOURCE_CORRECTION_FACTOR + (1._EB-WGT)*MIN(C_MAX,MAX(C_MIN,RAD_Q_SUM/KFST4_SUM))

! Reset the components of the corrective factor to zero.

RAD_Q_SUM = 0._EB
KFST4_SUM = 0._EB

T_USED(11)=T_USED(11) + CURRENT_TIME() - TNOW
END SUBROUTINE CALCULATE_RTE_SOURCE_CORRECTION_FACTOR


!> \brief Create or remove obstructions if necessary
!> \details Check if any obstructions are to be created or removed, and if so, call REASSIGN_WALL_CELLS to reset the WALL
!> boundary conditions where obstructions have been removed or created

SUBROUTINE CREATE_OR_REMOVE_OBSTRUCTIONS

DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   CALL OPEN_AND_CLOSE(T,NM)
ENDDO

CALL MPI_ALLREDUCE(MPI_IN_PLACE,OBST_CREATED_OR_REMOVED,INTEGER_ONE,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,IERR)

IF (OBST_CREATED_OR_REMOVED) THEN
   CALL EXCHANGE_GEOMETRY_INFO
   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL REASSIGN_WALL_CELLS(T,NM)
   ENDDO
   UPDATE_ALL_ANGLES = .TRUE.  ! Update all radiation angles the next time radiation is updated
ENDIF

END SUBROUTINE CREATE_OR_REMOVE_OBSTRUCTIONS


!> \brief Reassign global matrices
!> \details Global or local connected zone matrices potentially need to be reassigned after
!> CREATE_OR_REMOVE_OBSTRUCTIONS or CONNECTED_ZONES has been updated.  Zone connections are established
!> in DIVERGENCE_PART_1, so global reassignment must follow.

SUBROUTINE GLOBAL_MATRIX_REASSIGN(FORCE_REASSIGN)

LOGICAL, INTENT(IN) :: FORCE_REASSIGN

IF (OBST_CREATED_OR_REMOVED .OR. FORCE_REASSIGN) THEN
   SELECT CASE (PRES_FLAG)
      CASE (ULMAT_FLAG)
         DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
            IF(ALLOCATED(MESHES(NM)%ZONE_MESH)) CALL FINISH_ULMAT_SOLVER(NM)
            CALL ULMAT_SOLVER_SETUP(NM)
         ENDDO
         CALL STOP_CHECK(1)
      CASE (UGLMAT_FLAG,GLMAT_FLAG)
         IF(ALLOCATED(ZONE_SOLVE)) CALL FINISH_GLMAT_SOLVER
         CALL GLMAT_SOLVER_SETUP(-1) ! Initialize EWC_TYPE, copy wall types to HS
         CALL STOP_CHECK(1)
         CALL MESH_EXCHANGE(3)       ! Exchange guard cell info for IS_WALLT -> HS
         CALL GLMAT_SOLVER_SETUP(0)  ! Process coarse faces, copy updated wall types to HS
         CALL MESH_EXCHANGE(3)       ! Re-exchange guard cell info for IS_WALLT -> HS
         CALL GLMAT_SOLVER_SETUP(1)
         CALL MESH_EXCHANGE(3)       ! Exchange guard cell info for CCVAR(I,J,K,CGSC) -> HS
         CALL GLMAT_SOLVER_SETUP(2)
         CALL MESH_EXCHANGE(3)       ! Exchange guard cell info for CCVAR(I,J,K,UNKH) -> HS
         CALL GLMAT_SOLVER_SETUP(3)
         CALL STOP_CHECK(1)
   END SELECT
   OBST_CREATED_OR_REMOVED = .FALSE.
ENDIF

END SUBROUTINE GLOBAL_MATRIX_REASSIGN


!> \brief Gather all the CFL values and mesh indices to node 0 and then
!> write out to a special file the max value and mesh and indices of the max value.

SUBROUTINE WRITE_CFL_FILE

REAL(EB), DIMENSION(NMESHES) :: CFL_VALUES,VN_VALUES
INTEGER :: NM_CFL_MAX,NM_VN_MAX

DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   REAL_BUFFER_20(1,NM) = MESHES(NM)%CFL
   REAL_BUFFER_20(2,NM) = MESHES(NM)%ICFL
   REAL_BUFFER_20(3,NM) = MESHES(NM)%JCFL
   REAL_BUFFER_20(4,NM) = MESHES(NM)%KCFL
   REAL_BUFFER_20(5,NM) = MESHES(NM)%VN
   REAL_BUFFER_20(6,NM) = MESHES(NM)%I_VN
   REAL_BUFFER_20(7,NM) = MESHES(NM)%J_VN
   REAL_BUFFER_20(8,NM) = MESHES(NM)%K_VN
   REAL_BUFFER_20(9,NM) = MESHES(NM)%US(MESHES(NM)%ICFL-1,MESHES(NM)%JCFL  ,MESHES(NM)%KCFL  )
   REAL_BUFFER_20(10,NM)= MESHES(NM)%US(MESHES(NM)%ICFL  ,MESHES(NM)%JCFL  ,MESHES(NM)%KCFL  )
   REAL_BUFFER_20(11,NM)= MESHES(NM)%VS(MESHES(NM)%ICFL  ,MESHES(NM)%JCFL-1,MESHES(NM)%KCFL  )
   REAL_BUFFER_20(12,NM)= MESHES(NM)%VS(MESHES(NM)%ICFL  ,MESHES(NM)%JCFL  ,MESHES(NM)%KCFL  )
   REAL_BUFFER_20(13,NM)= MESHES(NM)%WS(MESHES(NM)%ICFL  ,MESHES(NM)%JCFL  ,MESHES(NM)%KCFL-1)
   REAL_BUFFER_20(14,NM)= MESHES(NM)%WS(MESHES(NM)%ICFL  ,MESHES(NM)%JCFL  ,MESHES(NM)%KCFL  )
   REAL_BUFFER_20(15,NM)= MESHES(NM)%DS(MESHES(NM)%ICFL  ,MESHES(NM)%JCFL  ,MESHES(NM)%KCFL  )
   REAL_BUFFER_20(16,NM)= MESHES(NM)%MU(MESHES(NM)%ICFL  ,MESHES(NM)%JCFL  ,MESHES(NM)%KCFL  )
   REAL_BUFFER_20(17,NM)= MESHES(NM)%Q( MESHES(NM)%ICFL  ,MESHES(NM)%JCFL  ,MESHES(NM)%KCFL  )*0.001_EB
   REAL_BUFFER_20(18,NM)= MESHES(NM)%MIX_TIME( MESHES(NM)%ICFL  ,MESHES(NM)%JCFL  ,MESHES(NM)%KCFL  )
ENDDO

IF (MY_RANK>0) THEN
   CALL MPI_GATHERV(REAL_BUFFER_20(1,DISPLS(MY_RANK)+1),COUNTS_20(MY_RANK),MPI_DOUBLE_PRECISION,REAL_BUFFER_20(1,1),&
                    COUNTS_20,DISPLS_20,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
ELSE
   CALL MPI_GATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,REAL_BUFFER_20(1,1),COUNTS_20,&
                    DISPLS_20,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
ENDIF

IF (MY_RANK==0) THEN
   CFL_VALUES(:) = REAL_BUFFER_20(1,:)
   VN_VALUES(:)  = REAL_BUFFER_20(5,:)
   NM_CFL_MAX  = MAXLOC(CFL_VALUES,DIM=1)
   NM_VN_MAX   = MAXLOC(VN_VALUES,DIM=1)
   WRITE(LU_CFL,'(I7,A,2(ES12.4,A),F6.3,A,4(I4,A),6(F7.3,A),4(ES12.4,A),F6.3,4(A,I4))') ICYC,',',T,',',DT,',',&
         MAXVAL(CFL_VALUES),',',NM_CFL_MAX,',',&
         NINT(REAL_BUFFER_20(2,NM_CFL_MAX)),',',NINT(REAL_BUFFER_20(3,NM_CFL_MAX)),',',NINT(REAL_BUFFER_20(4,NM_CFL_MAX)),',',&
         REAL_BUFFER_20( 9,NM_CFL_MAX),',',REAL_BUFFER_20(10,NM_CFL_MAX),',',REAL_BUFFER_20(11,NM_CFL_MAX),',',&
         REAL_BUFFER_20(12,NM_CFL_MAX),',',REAL_BUFFER_20(13,NM_CFL_MAX),',',REAL_BUFFER_20(14,NM_CFL_MAX),',',&
         REAL_BUFFER_20(15,NM_CFL_MAX),',',REAL_BUFFER_20(16,NM_CFL_MAX),',',REAL_BUFFER_20(17,NM_CFL_MAX),',',&
         REAL_BUFFER_20(18,NM_CFL_MAX),',',&
         MAXVAL(VN_VALUES),',',NM_VN_MAX,',',&
         NINT(REAL_BUFFER_20(6,NM_VN_MAX)),',',NINT(REAL_BUFFER_20(7,NM_VN_MAX)),',',NINT(REAL_BUFFER_20(8,NM_VN_MAX))
ENDIF

END SUBROUTINE WRITE_CFL_FILE


!> \brief Make sure that all MPI processes have the same STOP_STATUS

SUBROUTINE STOP_CHECK(END_CODE)

INTEGER, INTENT(IN) :: END_CODE
REAL(EB) :: TNOW

IF (N_MPI_PROCESSES>1) THEN
   TNOW = CURRENT_TIME()
   CALL MPI_ALLREDUCE(MPI_IN_PLACE,STOP_STATUS,INTEGER_ONE,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERR)
   T_USED(11)=T_USED(11) + CURRENT_TIME() - TNOW
ENDIF

SELECT CASE(STOP_STATUS)
   CASE(NO_STOP)
      RETURN
   CASE(USER_STOP)
      DIAGNOSTICS = .TRUE.
      IF (STOP_AT_ITER==0 .AND. .NOT.ALL(RADIATION_COMPLETED)) RETURN
END SELECT

IF (END_CODE==1) CALL END_FDS

END SUBROUTINE STOP_CHECK


!> \brief End FDS gracefully, even if there is an error

SUBROUTINE END_FDS

CHARACTER(MESSAGE_LENGTH) :: MESSAGE

IF (STOP_STATUS==NO_STOP .OR. STOP_STATUS==USER_STOP) CALL DUMP_TIMERS

IF (VERBOSE) WRITE(LU_ERR,'(A,I6,A)') ' MPI process ',MY_RANK,' has completed'

IF (MY_RANK==0) THEN

   ! Print out device activation times to the .out file

   IF (STOP_STATUS/=SETUP_STOP) CALL TIMINGS

   ! Print out stop status to .err and .out files

   SELECT CASE(STOP_STATUS)
      CASE(NO_STOP)
         WRITE(MESSAGE,'(A)') 'STOP: FDS completed successfully'
         IF (STATUS_FILES) CLOSE(LU_NOTREADY,STATUS='DELETE')
      CASE(INSTABILITY_STOP)
         WRITE(MESSAGE,'(A)') 'ERROR(374): Numerical Instability - FDS stopped'
      CASE(USER_STOP)
         WRITE(MESSAGE,'(A)') 'STOP: FDS stopped by user'
      CASE(SETUP_STOP)
         WRITE(MESSAGE,'(A)') 'ERROR: FDS was improperly set-up - FDS stopped'
      CASE(SETUP_ONLY_STOP)
         WRITE(MESSAGE,'(A)') 'STOP: Set-up only'
      CASE(CTRL_STOP)
         WRITE(MESSAGE,'(A)') 'STOP: FDS was stopped by KILL control function and completed successfully'
      CASE(TGA_ANALYSIS_STOP)
         WRITE(MESSAGE,'(A)') 'STOP: FDS performed a TGA analysis only and finished successfully'
      CASE(LEVELSET_STOP)
         WRITE(MESSAGE,'(A)') 'STOP: FDS performed a level set analysis only and finished successfully'
      CASE(REALIZABILITY_STOP)
         WRITE(MESSAGE,'(A)') 'ERROR: Unrealizable mass density - FDS stopped'
      CASE(ODE_STOP)
         WRITE(MESSAGE,'(A)') 'ERROR: Combustion ODE solver failure - FDS stopped'
      CASE(HEARTBEAT_STOP)
         WRITE(MESSAGE,'(A)') 'ERROR: External program failure - FDS stopped'
      CASE(CVODE_SUBSTEP_STOP)
         WRITE(MESSAGE,'(A)') 'STOP: CVODE substeps for 1st timestep are completed successfully - FDS stopped'
      CASE DEFAULT
         WRITE(MESSAGE,'(A)') 'null'
   END SELECT

   IF (MESSAGE/='null') THEN
      WRITE(LU_ERR,'(/A,A,A,A)') TRIM(MESSAGE),' (CHID: ',TRIM(CHID),')'
      IF (OUT_FILE_OPENED) WRITE(LU_OUTPUT,'(/A,A,A,A)') TRIM(MESSAGE),' (CHID: ',TRIM(CHID),')'
   ENDIF

ENDIF

! Free MPI send/recv requests

DO I=1,N_REQ1  ; CALL MPI_REQUEST_FREE(REQ1(I) ,IERR) ; ENDDO
DO I=1,N_REQ2  ; CALL MPI_REQUEST_FREE(REQ2(I) ,IERR) ; ENDDO
DO I=1,N_REQ3  ; CALL MPI_REQUEST_FREE(REQ3(I) ,IERR) ; ENDDO
DO I=1,N_REQ4  ; CALL MPI_REQUEST_FREE(REQ4(I) ,IERR) ; ENDDO
DO I=1,N_REQ5  ; CALL MPI_REQUEST_FREE(REQ5(I) ,IERR) ; ENDDO
DO I=1,N_REQ6  ; CALL MPI_REQUEST_FREE(REQ6(I) ,IERR) ; ENDDO
DO I=1,N_REQ7  ; CALL MPI_REQUEST_FREE(REQ7(I) ,IERR) ; ENDDO
DO I=1,N_REQ14 ; CALL MPI_REQUEST_FREE(REQ14(I),IERR) ; ENDDO
IF (EXCHANGE_OBST_MASS) THEN
   DO I=1,N_REQ15 ; CALL MPI_REQUEST_FREE(REQ15(I),IERR) ; ENDDO
ENDIF

! Shutdown MPI

CALL MPI_FINALIZE(IERR)

! Shutdown FDS

STOP

END SUBROUTINE END_FDS


!> \brief Zero out the integrals found in the divergence expression. 

SUBROUTINE INITIALIZE_DIVERGENCE_INTEGRALS

INTEGER :: IPZ

USUM = 0._EB
PSUM = 0._EB
DSUM = 0._EB

IF (OBST_CREATED_OR_REMOVED) THEN
   CONNECTED_ZONES = 0
   DO IPZ=0,N_ZONE
      CONNECTED_ZONES(IPZ,IPZ) = 1
   ENDDO
ENDIF

END SUBROUTINE INITIALIZE_DIVERGENCE_INTEGRALS


!> \brief Exchange information mesh to mesh needed for divergence integrals
!> \details Sum DSUM, PSUM and USUM over all meshes controlled by the active process, then reduce over all processes

SUBROUTINE EXCHANGE_DIVERGENCE_INFO

INTEGER :: IPZ
REAL(EB) :: TNOW

TNOW = CURRENT_TIME()

! Sum up the divergence integrals and zone connection matrix over the MPI processes

IF (N_MPI_PROCESSES>1) THEN
   CALL MPI_ALLREDUCE(MPI_IN_PLACE,DSUM,N_ZONE,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
   CALL MPI_ALLREDUCE(MPI_IN_PLACE,PSUM,N_ZONE,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
   CALL MPI_ALLREDUCE(MPI_IN_PLACE,USUM,N_ZONE,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
ENDIF

! Connect zones to others which are not directly connected

IF (OBST_CREATED_OR_REMOVED) THEN
   IF (N_MPI_PROCESSES>1) CALL MPI_ALLREDUCE(MPI_IN_PLACE,CONNECTED_ZONES,(N_ZONE+1)**2,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERR)
   DO IPZ=1,N_ZONE
      CONNECTED_ZONES = MATMUL(CONNECTED_ZONES,CONNECTED_ZONES)
      CONNECTED_ZONES = MIN(1,CONNECTED_ZONES)
   ENDDO
ENDIF

T_USED(11)=T_USED(11) + CURRENT_TIME() - TNOW
END SUBROUTINE EXCHANGE_DIVERGENCE_INFO


!> \brief Create arrays by which info is to be exchanged across meshes

SUBROUTINE INITIALIZE_MESH_EXCHANGE_1(NM)

INTEGER :: IMIN,IMAX,JMIN,JMAX,KMIN,KMAX,NOM,IOR,IW,N,IIO,JJO,KKO,NIC_R,ILC
INTEGER, INTENT(IN) :: NM
TYPE (MESH_TYPE), POINTER :: M2,M
TYPE (OMESH_TYPE), POINTER :: OM
TYPE (EXTERNAL_WALL_TYPE), POINTER :: EWC
TYPE (WALL_TYPE), POINTER :: WC
TYPE (BOUNDARY_COORD_TYPE), POINTER :: BC
TYPE (LAGRANGIAN_PARTICLE_CLASS_TYPE), POINTER :: LPC
TYPE (LAGRANGIAN_PARTICLE_TYPE), POINTER :: LP
TYPE (STORAGE_TYPE), POINTER :: OS
LOGICAL :: FOUND

M=>MESHES(NM)

ALLOCATE(MESHES(NM)%OMESH(NMESHES))

ALLOCATE(M%CONNECTED_MESH(NMESHES)) ; M%CONNECTED_MESH = .FALSE.

OTHER_MESH_LOOP: DO NOM=1,NMESHES

   OM => M%OMESH(NOM)
   M2 => MESHES(NOM)

   IMIN=0
   IMAX=M2%IBP1
   JMIN=0
   JMAX=M2%JBP1
   KMIN=0
   KMAX=M2%KBP1

   ! DEFINITION NIC_R: (Number of Interpolated Cells -- Receiving) Number of cells in mesh NOM that abut mesh NM.
   ! DEFINITION IIO_R: Array of length NIC_R of I indices of the abutting (outside) cells
   ! DEFINITION JJO_R: Array of length NIC_R of J indices of the abutting (outside) cells
   ! DEFINITION KKO_R: Array of length NIC_R of K indices of the abutting (outside) cells
   ! DEFINITION IOR_R: Array of length NIC_R of orientation of the external boundary cell
   ! DEFINITION NIC_MIN: For external wall cell IW of mesh NM, the indices of abutting cells start with NIC_MIN and end with NIC_MAX
   ! DEFINITION NIC_MAX: For external wall cell IW of mesh NM, the indices of abutting cells start with NIC_MIN and end with NIC_MAX

   OM%NIC_R = 0
   FOUND = .FALSE.

   SEARCH_LOOP: DO IW=1,M%N_EXTERNAL_WALL_CELLS

      EWC => M%EXTERNAL_WALL(IW)
      IF (EWC%NOM/=NOM) CYCLE SEARCH_LOOP
      WC => M%WALL(IW)
      BC => M%BOUNDARY_COORD(WC%BC_INDEX)
      EWC%NIC_MIN = OM%NIC_R + 1
      EWC%NIC = (EWC%IIO_MAX-EWC%IIO_MIN+1)*(EWC%JJO_MAX-EWC%JJO_MIN+1)*(EWC%KKO_MAX-EWC%KKO_MIN+1)
      OM%NIC_R = OM%NIC_R + EWC%NIC
      EWC%NIC_MAX = OM%NIC_R
      FOUND = .TRUE.
      M%CONNECTED_MESH(NOM) = .TRUE.
      SELECT CASE(BC%IOR)
         CASE( 1)
            IMIN=MAX(IMIN,EWC%IIO_MIN-1)
         CASE(-1)
            IMAX=MIN(IMAX,EWC%IIO_MAX+1)
         CASE( 2)
            JMIN=MAX(JMIN,EWC%JJO_MIN-1)
         CASE(-2)
            JMAX=MIN(JMAX,EWC%JJO_MAX+1)
         CASE( 3)
            KMIN=MAX(KMIN,EWC%KKO_MIN-1)
         CASE(-3)
            KMAX=MIN(KMAX,EWC%KKO_MAX+1)
      END SELECT

      SELECT CASE(ABS(BC%IOR))
         CASE(1) ; EWC%AREA_RATIO = &
                   M%DY(BC%JJ)*M%DZ(BC%KK)/((M2%Y(EWC%JJO_MAX)-M2%Y(EWC%JJO_MIN-1))*(M2%Z(EWC%KKO_MAX)-M2%Z(EWC%KKO_MIN-1)))
         CASE(2) ; EWC%AREA_RATIO = &
                   M%DX(BC%II)*M%DZ(BC%KK)/((M2%X(EWC%IIO_MAX)-M2%X(EWC%IIO_MIN-1))*(M2%Z(EWC%KKO_MAX)-M2%Z(EWC%KKO_MIN-1)))
         CASE(3) ; EWC%AREA_RATIO = &
                   M%DX(BC%II)*M%DY(BC%JJ)/((M2%X(EWC%IIO_MAX)-M2%X(EWC%IIO_MIN-1))*(M2%Y(EWC%JJO_MAX)-M2%Y(EWC%JJO_MIN-1)))
      END SELECT
   ENDDO SEARCH_LOOP

   ! Allocate arrays to hold indices of arrays for MPI exchanges

   IF (OM%NIC_R>0) THEN
      ALLOCATE(OM%IIO_R(OM%NIC_R))
      ALLOCATE(OM%JJO_R(OM%NIC_R))
      ALLOCATE(OM%KKO_R(OM%NIC_R))
      ALLOCATE(OM%IOR_R(OM%NIC_R))
      NIC_R = 0
      INDEX_LOOP: DO IW=1,M%N_EXTERNAL_WALL_CELLS
         EWC => M%EXTERNAL_WALL(IW)
         IF (EWC%NOM/=NOM) CYCLE INDEX_LOOP
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  NIC_R = NIC_R + 1
                  IOR = M%BOUNDARY_COORD(M%WALL(IW)%BC_INDEX)%IOR
                  OM%IIO_R(NIC_R) = IIO
                  OM%JJO_R(NIC_R) = JJO
                  OM%KKO_R(NIC_R) = KKO
                  OM%IOR_R(NIC_R) = IOR
               ENDDO
            ENDDO
         ENDDO
      ENDDO INDEX_LOOP

   ENDIF

   ! For a single mesh with PERIODIC boundaries, allocate the entire OMESH array

   IF (IMIN>IMAX .OR. NM==NOM) THEN ; IMIN=0 ; IMAX=M2%IBP1 ; ENDIF
   IF (JMIN>JMAX .OR. NM==NOM) THEN ; JMIN=0 ; JMAX=M2%JBP1 ; ENDIF
   IF (KMIN>KMAX .OR. NM==NOM) THEN ; KMIN=0 ; KMAX=M2%KBP1 ; ENDIF

   ! Embedded meshes. This is the case where mesh NOM is completely inside mesh NM. Mesh NM cannot "see" mesh NOM because mesh NOM
   ! is not connected at the external boundary of mesh NM. The variable CONNECTED_MESH is needed to save this information.

   IF ( NM/=NOM .AND. M2%XS>=M%XS .AND. M2%XF<=M%XF .AND. M2%YS>=M%YS .AND. M2%YF<=M%YF .AND. M2%ZS>=M%ZS .AND. M2%ZF<=M%ZF ) THEN
      FOUND = .TRUE.
      M%CONNECTED_MESH(NOM) = .TRUE.
   ENDIF

   ! Exit the other mesh loop if no neighboring meshes found

   IF (.NOT.FOUND) CYCLE OTHER_MESH_LOOP

   ! Save the dimensions of the volume of cells from mesh NOM whose data is received by mesh NM
   ! DEFINITION I_MIN_R: Starting I index of cell block of mesh NOM whose info is to be received by mesh NM in MPI exchanges.
   ! DEFINITION I_MAX_R: Ending   I index of cell block of mesh NOM whose info is to be received by mesh NM in MPI exchanges.
   ! DEFINITION J_MIN_R: Starting J index of cell block of mesh NOM whose info is to be received by mesh NM in MPI exchanges.
   ! DEFINITION J_MAX_R: Ending   J index of cell block of mesh NOM whose info is to be received by mesh NM in MPI exchanges.
   ! DEFINITION K_MIN_R: Starting K index of cell block of mesh NOM whose info is to be received by mesh NM in MPI exchanges.
   ! DEFINITION K_MAX_R: Ending   K index of cell block of mesh NOM whose info is to be received by mesh NM in MPI exchanges.

   OM%I_MIN_R = IMIN
   OM%I_MAX_R = IMAX
   OM%J_MIN_R = JMIN
   OM%J_MAX_R = JMAX
   OM%K_MIN_R = KMIN
   OM%K_MAX_R = KMAX

   ! Allocate the arrays that hold information about the other meshes (OMESH)

   ALLOCATE(OM% RHO(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
   ALLOCATE(OM%RHOS(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
   OM%RHO  = RHOA
   OM%RHOS = RHOA
   ALLOCATE(OM% D(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
   ALLOCATE(OM%DS(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
   OM%D  = 0._EB
   OM%DS = 0._EB
   ALLOCATE(OM%  MU(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
   OM%MU = 0._EB
   ALLOCATE(OM%    H(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
   ALLOCATE(OM%   HS(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
   OM%H  = 0._EB
   OM%HS = 0._EB
   ALLOCATE(OM%   U(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
   ALLOCATE(OM%  US(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
   OM%U  = U0
   OM%US = U0
   ALLOCATE(OM%   V(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
   ALLOCATE(OM%  VS(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
   OM%V  = V0
   OM%VS = V0
   ALLOCATE(OM%   W(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
   ALLOCATE(OM%  WS(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
   OM%W  = W0
   OM%WS = W0
   ALLOCATE(OM% FVX(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
   ALLOCATE(OM% FVY(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
   ALLOCATE(OM% FVZ(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
   OM%FVX = 0._EB
   OM%FVY = 0._EB
   OM%FVZ = 0._EB
   ALLOCATE(OM%KRES(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
   OM%KRES = 0._EB
   ALLOCATE(OM%Q(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
   OM%Q = 0._EB

   ALLOCATE(OM%  ZZ(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX,N_TOTAL_SCALARS))
   ALLOCATE(OM% ZZS(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX,N_TOTAL_SCALARS))
   DO N=1,N_TRACKED_SPECIES
      OM%ZZ(:,:,:,N)  = SPECIES_MIXTURE(N)%ZZ0
      OM%ZZS(:,:,:,N) = SPECIES_MIXTURE(N)%ZZ0
   ENDDO
   DO N=N_TRACKED_SPECIES+1,N_TOTAL_SCALARS
      OM%ZZ(:,:,:,N)  = INITIAL_UNMIXED_FRACTION
      OM%ZZS(:,:,:,N) = INITIAL_UNMIXED_FRACTION
   ENDDO

   IF (LEVEL_SET_MODE>0 .OR. TERRAIN_CASE) THEN
      ALLOCATE(OM%PHI_LS(IMIN:IMAX,JMIN:JMAX))  ; OM%PHI_LS   = -1._EB
      ALLOCATE(OM%PHI1_LS(IMIN:IMAX,JMIN:JMAX)) ; OM%PHI1_LS  = -1._EB
      ALLOCATE(OM%U_LS(IMIN:IMAX,JMIN:JMAX))    ; OM%U_LS     =  0._EB
      ALLOCATE(OM%V_LS(IMIN:IMAX,JMIN:JMAX))    ; OM%V_LS     =  0._EB
      ALLOCATE(OM%Z_LS(IMIN:IMAX,JMIN:JMAX))    ; OM%Z_LS     =  0._EB
   ENDIF

ENDDO OTHER_MESH_LOOP

! Allocate PARTICLE_STORAGE to hold particle data for a RESTART

DO ILC=1,N_LAGRANGIAN_CLASSES
   LPC => LAGRANGIAN_PARTICLE_CLASS(ILC)
   CALL ALLOCATE_STORAGE(NM,LP_INDEX=1,LPC_INDEX=ILC,SURF_INDEX=LPC%SURF_INDEX)
   LP => M%LAGRANGIAN_PARTICLE(1)
   OS => LPC%PARTICLE_STORAGE
   LPC%N_REALS=0 ; LPC%N_INTEGERS=0 ; LPC%N_LOGICALS=0
   CALL PACK_PARTICLE(NM,OS,LP,ILC,LPC%N_REALS,LPC%N_INTEGERS,LPC%N_LOGICALS,&
                      UNPACK_IT=.FALSE.,COUNT_ONLY=.TRUE.,CHECK_BOUNDS=.FALSE.)
   IF (.NOT.ALLOCATED(LPC%PARTICLE_STORAGE%REALS))    ALLOCATE(LPC%PARTICLE_STORAGE%REALS(LPC%N_REALS))
   IF (.NOT.ALLOCATED(LPC%PARTICLE_STORAGE%INTEGERS)) ALLOCATE(LPC%PARTICLE_STORAGE%INTEGERS(LPC%N_INTEGERS))
   IF (.NOT.ALLOCATED(LPC%PARTICLE_STORAGE%LOGICALS)) ALLOCATE(LPC%PARTICLE_STORAGE%LOGICALS(LPC%N_LOGICALS))
   CALL NULLIFY_PARTICLE(NM,1)
ENDDO

! Allocate SEND_BUFFER to hold particles that leave mesh NM and a RECV_BUFFER for particles entering a NEIGHBORING_MESH

IF (MPI_PARTICLE_EXCHANGE) THEN

   NEIGHBORING_MESH_LOOP: DO N=1,M%N_NEIGHBORING_MESHES
      OM => M%OMESH(M%NEIGHBORING_MESH(N))
      ALLOCATE(OM%PARTICLE_SEND_BUFFER%REALS(100))
      ALLOCATE(OM%PARTICLE_SEND_BUFFER%INTEGERS(100))
      ALLOCATE(OM%PARTICLE_SEND_BUFFER%LOGICALS(100))
      ALLOCATE(OM%PARTICLE_RECV_BUFFER%REALS(100))
      ALLOCATE(OM%PARTICLE_RECV_BUFFER%INTEGERS(100))
      ALLOCATE(OM%PARTICLE_RECV_BUFFER%LOGICALS(100))
   ENDDO NEIGHBORING_MESH_LOOP

ENDIF

END SUBROUTINE INITIALIZE_MESH_EXCHANGE_1


!> \brief Allocate arrays used for MPI exchange of radiation data
!> \param NM Mesh number
!> \details Allocate arrays to send (IL_S) and receive (IL_R) the radiation intensity (IL) at interpolated boundaries.
!> MESHES(NM)%OMESH(NOM)%IL_S are the intensities in mesh NM that are just outside the boundary of mesh NOM. IL_S is populated
!> in radi.f90 and then sent to MESHES(NOM)%OMESH(NM)%IL_R in MESH_EXCHANGE. IL_R holds the intensities until they are
!> transferred to the ghost cells of MESHES(NOM)%IL in radi.f90. The IL_S and IL_R arrays are indexed by NIC_S and NIC_R.

SUBROUTINE INITIALIZE_RADIATION_EXCHANGE(NM)

INTEGER :: NOM
INTEGER, INTENT(IN) :: NM
TYPE (MESH_TYPE), POINTER :: M

M=>MESHES(NM)

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

END SUBROUTINE INITIALIZE_RADIATION_EXCHANGE


!> \brief Bordering meshes tell their neighbors how many exposed back wall cells they expect information for.
!> \param PASS_INDEX An integer with value 1 or 2 indicating whether this is the first or second call

SUBROUTINE INITIALIZE_BACK_WALL_EXCHANGE(PASS_INDEX)

INTEGER, INTENT(IN) :: PASS_INDEX
INTEGER :: NOM,IW,ITW,TRUE_COUNT,II
TYPE(WALL_TYPE), POINTER :: WC
TYPE(THIN_WALL_TYPE), POINTER :: TW
TYPE(STORAGE_TYPE), POINTER :: OS
TYPE(STORAGE_TYPE), ALLOCATABLE :: DUMMY

IF (PASS_INDEX==1) THEN

   ! Create a list of WALL and THIN_WALL indices needed by each mesh and store then in
   ! MESHES(NM)%OMESH(NOM)%[THIN_]WALL_RECV_BUFFER%ITEM_INDEX(1:N_ITEMS)

   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL FIND_WALL_BACK_INDICES(NM)
   ENDDO

   CALL STOP_CHECK(1)

   ! Adjust the thickness and internal noding of HT3D surfaces

   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL ADJUST_HT3D_WALL_CELLS(NM)
   ENDDO

ENDIF

! Current mesh sends to neighboring meshes the number of WALL and THIN_WALL cells that it expects to be SENT

CALL POST_RECEIVES(8)
CALL MESH_EXCHANGE(8)

! Allocate WALL_SEND_BUFFER and THIN_WALL_SEND_BUFFER for each neighboring mesh.
! These derived types hold the information that is to be sent.

DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   DO NOM=1,NMESHES
      IF (NM==NOM) CYCLE
      OS => MESHES(NM)%OMESH(NOM)%WALL_SEND_BUFFER
      IF (ALLOCATED(OS%ITEM_INDEX)) DEALLOCATE(OS%ITEM_INDEX)
      IF (ALLOCATED(OS%SURF_INDEX)) DEALLOCATE(OS%SURF_INDEX)
      IF (OS%N_ITEMS>0) ALLOCATE(OS%ITEM_INDEX(OS%N_ITEMS_DIM))
      IF (OS%N_ITEMS>0) ALLOCATE(OS%SURF_INDEX(OS%N_ITEMS_DIM))
      OS => MESHES(NM)%OMESH(NOM)%THIN_WALL_SEND_BUFFER
      IF (ALLOCATED(OS%ITEM_INDEX)) DEALLOCATE(OS%ITEM_INDEX)
      IF (ALLOCATED(OS%SURF_INDEX)) DEALLOCATE(OS%SURF_INDEX)
      IF (OS%N_ITEMS>0) ALLOCATE(OS%ITEM_INDEX(OS%N_ITEMS_DIM))
      IF (OS%N_ITEMS>0) ALLOCATE(OS%SURF_INDEX(OS%N_ITEMS_DIM))
   ENDDO
ENDDO

! Each mesh sends its neighbors an array of WALL and THIN_WALL indices (ITEM_INDEX) that it expects to be SENT,
! along with the SURF_INDEX for each WALL or THIN_WALL.

CALL POST_RECEIVES(9)
CALL MESH_EXCHANGE(9)

IF (PASS_INDEX==1) THEN

   ! Set up storage arrays for packing WALL and THIN_WALL variables during a RESTART.

   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      M => MESHES(NM)
      OS => M%WALL_STORAGE
      DO IW=1,M%N_WALL_CELLS
         WC => M%WALL(IW)
         OS%N_REALS_DIM    = MAX(OS%N_REALS_DIM,WC%N_REALS)
         OS%N_INTEGERS_DIM = MAX(OS%N_INTEGERS_DIM,WC%N_INTEGERS)
         OS%N_LOGICALS_DIM = MAX(OS%N_LOGICALS_DIM,WC%N_LOGICALS)
      ENDDO
      DO ITW=1,M%N_THIN_WALL_CELLS
         TW => M%THIN_WALL(ITW)
         OS%N_REALS_DIM    = MAX(OS%N_REALS_DIM,TW%N_REALS)
         OS%N_INTEGERS_DIM = MAX(OS%N_INTEGERS_DIM,TW%N_INTEGERS)
         OS%N_LOGICALS_DIM = MAX(OS%N_LOGICALS_DIM,TW%N_LOGICALS)
      ENDDO
      ALLOCATE(OS%REALS(OS%N_REALS_DIM))
      ALLOCATE(OS%INTEGERS(OS%N_INTEGERS_DIM))
      ALLOCATE(OS%LOGICALS(OS%N_LOGICALS_DIM))
   ENDDO

ENDIF

! Allocate arrays to hold real, integer and logical variables for the WALL_SEND_BUFFER

MESH_LOOP_1A: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   M => MESHES(NM)
   MESH_LOOP_1B: DO NOM=1,NMESHES
      IF (PROCESS(NM)==PROCESS(NOM)) CYCLE MESH_LOOP_1B
      M3 => M%OMESH(NOM)
      OS => M3%WALL_SEND_BUFFER
      IF (OS%N_ITEMS==0) CYCLE MESH_LOOP_1B
      DO I=1,OS%N_ITEMS
         WC => M%WALL(OS%ITEM_INDEX(I))
         OS%N_REALS_DIM    = OS%N_REALS_DIM    + WC%N_REALS
         OS%N_INTEGERS_DIM = OS%N_INTEGERS_DIM + WC%N_INTEGERS
         OS%N_LOGICALS_DIM = OS%N_LOGICALS_DIM + WC%N_LOGICALS
      ENDDO
      IF (ALLOCATED(OS%REALS)) DEALLOCATE(OS%REALS)
      IF (ALLOCATED(OS%INTEGERS)) DEALLOCATE(OS%INTEGERS)
      IF (ALLOCATED(OS%LOGICALS)) DEALLOCATE(OS%LOGICALS)
      ALLOCATE(OS%REALS(OS%N_REALS_DIM))
      ALLOCATE(OS%INTEGERS(OS%N_INTEGERS_DIM))
      ALLOCATE(OS%LOGICALS(OS%N_LOGICALS_DIM))
   ENDDO MESH_LOOP_1B
ENDDO MESH_LOOP_1A

! Allocate arrays to hold real, integer and logical variables for the THIN_WALL_SEND_BUFFER

MESH_LOOP_3A: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   M => MESHES(NM)
   MESH_LOOP_3B: DO NOM=1,NMESHES
      IF (PROCESS(NM)==PROCESS(NOM)) CYCLE MESH_LOOP_3B
      M3 => M%OMESH(NOM)
      OS => M3%THIN_WALL_SEND_BUFFER
      IF (OS%N_ITEMS==0) CYCLE MESH_LOOP_3B
      DO I=1,OS%N_ITEMS
         TW => M%THIN_WALL(OS%ITEM_INDEX(I))
         OS%N_REALS_DIM    = OS%N_REALS_DIM    + TW%N_REALS
         OS%N_INTEGERS_DIM = OS%N_INTEGERS_DIM + TW%N_INTEGERS
         OS%N_LOGICALS_DIM = OS%N_LOGICALS_DIM + TW%N_LOGICALS
      ENDDO
      IF (ALLOCATED(OS%REALS)) DEALLOCATE(OS%REALS)
      IF (ALLOCATED(OS%INTEGERS)) DEALLOCATE(OS%INTEGERS)
      IF (ALLOCATED(OS%LOGICALS)) DEALLOCATE(OS%LOGICALS)
      ALLOCATE(OS%REALS(OS%N_REALS_DIM))
      ALLOCATE(OS%INTEGERS(OS%N_INTEGERS_DIM))
      ALLOCATE(OS%LOGICALS(OS%N_LOGICALS_DIM))
   ENDDO MESH_LOOP_3B
ENDDO MESH_LOOP_3A

! Each mesh sends its neighbors the size of the real, integer, and logical buffer arrays

CALL POST_RECEIVES(19)
CALL MESH_EXCHANGE(19)

! Broadcast the required dimension of the WALL arrays for neighboring meshes

IF (N_MPI_PROCESSES>1) THEN
   DO NM=1,NMESHES
      IF (MPI_COMM_NEIGHBORS(NM)==MPI_COMM_NULL) CYCLE
      M => MESHES(NM)
      CALL MPI_BCAST(M%N_WALL_CELLS_DIM,1,MPI_INTEGER,MPI_COMM_NEIGHBORS_ROOT(NM),MPI_COMM_NEIGHBORS(NM),IERR)
   ENDDO
ENDIF

! Each mesh (NM) sets up WALL cells that it expects to receive from its neighbors (NOM)

MESH_LOOP_2A: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   M => MESHES(NM)
   MESH_LOOP_2B: DO NOM=1,NMESHES
      IF (PROCESS(NM)==PROCESS(NOM)) CYCLE MESH_LOOP_2B
      OS => M%OMESH(NOM)%WALL_RECV_BUFFER
      IF (OS%N_ITEMS==0) CYCLE MESH_LOOP_2B
      IF (.NOT.ALLOCATED(MESHES(NOM)%WALL)) ALLOCATE(MESHES(NOM)%WALL(0:MESHES(NOM)%N_WALL_CELLS_DIM))
      IF (ALLOCATED(OS%REALS)) DEALLOCATE(OS%REALS)
      IF (ALLOCATED(OS%INTEGERS)) DEALLOCATE(OS%INTEGERS)
      IF (ALLOCATED(OS%LOGICALS)) DEALLOCATE(OS%LOGICALS)
      ALLOCATE(OS%REALS(OS%N_REALS_DIM))
      ALLOCATE(OS%INTEGERS(OS%N_INTEGERS_DIM))
      ALLOCATE(OS%LOGICALS(OS%N_LOGICALS_DIM))
      DO I=1,OS%N_ITEMS
         CALL ALLOCATE_STORAGE(NOM,WALL_INDEX=OS%ITEM_INDEX(I),SURF_INDEX=OS%SURF_INDEX(I))
      ENDDO
   ENDDO MESH_LOOP_2B
ENDDO MESH_LOOP_2A

! Broadcast the required dimension of the WALL arrays for neighboring meshes

IF (N_MPI_PROCESSES>1) THEN
   DO NM=1,NMESHES
      IF (MPI_COMM_NEIGHBORS(NM)==MPI_COMM_NULL) CYCLE
      M => MESHES(NM)
      CALL MPI_BCAST(M%N_THIN_WALL_CELLS_DIM,1,MPI_INTEGER,MPI_COMM_NEIGHBORS_ROOT(NM),MPI_COMM_NEIGHBORS(NM),IERR)
   ENDDO
ENDIF

! Each mesh (NM) sets up THIN_WALL cells that it expects to receive from its neighbors (NOM)

MESH_LOOP_4A: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   M => MESHES(NM)
   MESH_LOOP_4B: DO NOM=1,NMESHES
      IF (PROCESS(NM)==PROCESS(NOM)) CYCLE MESH_LOOP_4B
      OS => M%OMESH(NOM)%THIN_WALL_RECV_BUFFER
      IF (OS%N_ITEMS==0) CYCLE MESH_LOOP_4B
      IF (.NOT.ALLOCATED(MESHES(NOM)%THIN_WALL)) ALLOCATE(MESHES(NOM)%THIN_WALL(0:MESHES(NOM)%N_THIN_WALL_CELLS_DIM))
      IF (ALLOCATED(OS%REALS)) DEALLOCATE(OS%REALS)
      IF (ALLOCATED(OS%INTEGERS)) DEALLOCATE(OS%INTEGERS)
      IF (ALLOCATED(OS%LOGICALS)) DEALLOCATE(OS%LOGICALS)
      ALLOCATE(OS%REALS(OS%N_REALS_DIM))
      ALLOCATE(OS%INTEGERS(OS%N_INTEGERS_DIM))
      ALLOCATE(OS%LOGICALS(OS%N_LOGICALS_DIM))
      DO I=1,OS%N_ITEMS
         CALL ALLOCATE_STORAGE(NOM,THIN_WALL_INDEX=OS%ITEM_INDEX(I),SURF_INDEX=OS%SURF_INDEX(I))
      ENDDO
   ENDDO MESH_LOOP_4B
ENDDO MESH_LOOP_4A

! Check to see if any process has an error. If so, stop the run.

CALL STOP_CHECK(1)

! Set up persistent SEND and RECV calls for MPI communication of WALL and THIN_WALL buffer arrays

IF (PASS_INDEX==2) THEN  ! Free the previously set permanent send/receive
   DO I=1,N_REQ6  ; CALL MPI_REQUEST_FREE(REQ6(I) ,IERR) ; ENDDO
   N_REQ6 = 0
ENDIF

CALL POST_RECEIVES(10)
CALL MESH_EXCHANGE(10)

IF (PASS_INDEX==1) THEN

   ! Exchange WALL and THIN_WALL cells

   CALL MESH_EXCHANGE(6)

   ! Initialize 3-D solid interpolation arrays

   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL INITIALIZE_HT3D_WALL_CELLS(NM)
   ENDDO

   ! Reduce the size of the send and receive buffer by including only those wall cells needed by other meshes.
   ! The key parameter is M%OM%WALL_RECV_BUFFER%SAVE_FLAG(1:N_ITEMS). If T, include the wall cell on the short list.

   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      M => MESHES(NM)
      DO NOM=1,NMESHES
         IF (NM==NOM) CYCLE
         OS => M%OMESH(NOM)%WALL_RECV_BUFFER
         IF (OS%N_ITEMS==0) CYCLE
         IF (ALLOCATED(DUMMY)) DEALLOCATE(DUMMY) ; ALLOCATE(DUMMY)
         DUMMY%N_ITEMS = OS%N_ITEMS
         ALLOCATE(DUMMY%ITEM_INDEX(OS%N_ITEMS)) ; DUMMY%ITEM_INDEX(1:OS%N_ITEMS) = OS%ITEM_INDEX(1:OS%N_ITEMS)
         ALLOCATE(DUMMY%SURF_INDEX(OS%N_ITEMS)) ; DUMMY%SURF_INDEX(1:OS%N_ITEMS) = OS%SURF_INDEX(1:OS%N_ITEMS)
         ALLOCATE(DUMMY%SAVE_FLAG(OS%N_ITEMS))  ; DUMMY%SAVE_FLAG(1:OS%N_ITEMS)  = OS%SAVE_FLAG(1:OS%N_ITEMS)
         TRUE_COUNT = COUNT(OS%SAVE_FLAG(1:OS%N_ITEMS))
         DEALLOCATE(OS%ITEM_INDEX) ; ALLOCATE(OS%ITEM_INDEX(TRUE_COUNT))
         DEALLOCATE(OS%SURF_INDEX) ; ALLOCATE(OS%SURF_INDEX(TRUE_COUNT))
         DEALLOCATE(OS%SAVE_FLAG)  ; ALLOCATE(OS%SAVE_FLAG(TRUE_COUNT))
         OS%N_ITEMS     = TRUE_COUNT
         OS%N_ITEMS_DIM = TRUE_COUNT
         II = 0
         DO I=1,DUMMY%N_ITEMS
            IF (DUMMY%SAVE_FLAG(I)) THEN
               II = II + 1
               OS%ITEM_INDEX(II) = DUMMY%ITEM_INDEX(I)
               OS%SURF_INDEX(II) = DUMMY%SURF_INDEX(I)
               OS%SAVE_FLAG(II)  = DUMMY%SAVE_FLAG(I)
            ENDIF
         ENDDO
      ENDDO
   ENDDO

ENDIF

END SUBROUTINE INITIALIZE_BACK_WALL_EXCHANGE


!> \brief Initialize the array PRESSURE_ZONE(I,J,K) for all meshes

SUBROUTINE INITIALIZE_PRESSURE_ZONES

USE GEOMETRY_FUNCTIONS, ONLY: ASSIGN_PRESSURE_ZONE,SEARCH_OTHER_MESHES
TYPE(P_ZONE_TYPE), POINTER :: PZ
INTEGER :: N,N_OVERLAP,I,J,K,NOM,IC,IZ,IZZ,IW,IOR,IZO,N_ZONE_ORIG
LOGICAL :: FOUND_UNASSIGNED_CELL
REAL(EB), ALLOCATABLE, DIMENSION(:) :: REAL_BUFFER
INTEGER, ALLOCATABLE, DIMENSION(:) :: INTEGER_BUFFER,INTEGER_BUFFER_4,NEW_ZONE_INDEX
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC

! Ensure that all cells have been assigned a pressure zone, even within solids

IF (NO_PRESSURE_ZONES) THEN
   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      M => MESHES(NM)
      M%PRESSURE_ZONE = 0
   ENDDO
   RETURN
ENDIF

! For each explicitly specified ZONE, populate the array PRESSURE_ZONE(I,J,K) with the ZONE index, N

MESH_LOOP_1: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   ZONE_LOOP: DO N=1,N_ZONE
      PZ => P_ZONE(N)
      CALL SEARCH_OTHER_MESHES(PZ%X,PZ%Y,PZ%Z,NOM,I,J,K)
      IF (NOM/=NM) CYCLE ZONE_LOOP
      CALL ASSIGN_PRESSURE_ZONE(NM,I,J,K,N,N_OVERLAP)
      IF (N_OVERLAP>0) THEN
         WRITE(LU_ERR,'(A,I2,A,I2,A,I4)') 'ERROR(872): ZONE ',N,' overlaps ZONE ',N_OVERLAP,' in MESH ',NM
         STOP_STATUS = SETUP_STOP
         EXIT MESH_LOOP_1
      ENDIF
   ENDDO ZONE_LOOP
ENDDO MESH_LOOP_1

CALL MPI_ALLREDUCE(MPI_IN_PLACE,STOP_STATUS,INTEGER_ONE,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERR)
IF (STOP_STATUS/=0) RETURN

! Propagate PRESSURE_ZONE values that exist at OPEN boundaries and INTERPOLATED boundaries

CALL ZONE_BOUNDARY_EXCHANGE

IF (STOP_STATUS/=0) RETURN

! Starting with the first mesh, look for grid cells that have not yet been assigned a PRESSURE_ZONE value.
! If found (FOUND_UNASSIGNED_CELL=T), exit the mesh loop, tell all processes that a new P_ZONE needs to be created.

MESH_LOOP_2: DO NM=1,NMESHES
   FOUND_UNASSIGNED_CELL = .TRUE.
   SEARCH_LOOP: DO WHILE(FOUND_UNASSIGNED_CELL)
      FOUND_UNASSIGNED_CELL = .FALSE.
      IF (MY_RANK==PROCESS(NM)) THEN
         M => MESHES(NM)
         K_LOOP: DO K=1,M%KBAR
            J_LOOP: DO J=1,M%JBAR
               I_LOOP: DO I=1,M%IBAR
                  IF (CC_IBM) THEN
                     IF (M%CCVAR(I,J,K,1)==1) CYCLE I_LOOP
                  ENDIF
                  IF (M%PRESSURE_ZONE(I,J,K)<0 .AND. .NOT.M%CELL(M%CELL_INDEX(I,J,K))%SOLID) THEN
                     FOUND_UNASSIGNED_CELL = .TRUE.
                     EXIT K_LOOP
                  ENDIF
               ENDDO I_LOOP
            ENDDO J_LOOP
         ENDDO K_LOOP
      ENDIF
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,FOUND_UNASSIGNED_CELL,INTEGER_ONE,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,IERR)
      IF (FOUND_UNASSIGNED_CELL) THEN
         CALL REALLOCATE_P_ZONE(N_ZONE,N_ZONE+1)
         N_ZONE = N_ZONE + 1
         IF (MY_RANK==PROCESS(NM)) CALL ASSIGN_PRESSURE_ZONE(NM,I,J,K,N_ZONE,N_OVERLAP)
         CALL ZONE_BOUNDARY_EXCHANGE
      ENDIF
   ENDDO SEARCH_LOOP
ENDDO MESH_LOOP_2

! Ensure that all cells have been assigned a pressure zone, even within solids.
! Also, compute the volume, number of cells, and cell indices of each pressure zone.

DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   M => MESHES(NM)
   WHERE (M%PRESSURE_ZONE<0) M%PRESSURE_ZONE = 0
   DO K=1,M%KBAR
      DO J=1,M%JBAR
         DO I=1,M%IBAR
            IZ = M%PRESSURE_ZONE(I,J,K)
            IF (IZ==0) CYCLE
            IC = M%CELL_INDEX(I,J,K)
            IF (.NOT.M%CELL(IC)%SOLID) THEN
               P_ZONE(IZ)%VOLUME  = P_ZONE(IZ)%VOLUME  + M%DX(I)*M%RC(I)*M%DY(J)*M%DZ(K)
               P_ZONE(IZ)%N_CELLS = P_ZONE(IZ)%N_CELLS + 1
               P_ZONE(IZ)%CELL_INDICES = (/I,J,K/)
               P_ZONE(IZ)%MESH_INDEX   = NM
            ENDIF
         ENDDO
      ENDDO
   ENDDO
ENDDO

! Compute the volume and number of cells of each pressure zone, as well as a single cell within the pressure zone

IF (N_MPI_PROCESSES>1 .AND. N_ZONE>0) THEN
   ALLOCATE(REAL_BUFFER(N_ZONE))
   ALLOCATE(INTEGER_BUFFER(N_ZONE))     ; INTEGER_BUFFER = 0
   ALLOCATE(INTEGER_BUFFER_4(4*N_ZONE)) ; INTEGER_BUFFER_4 = 0
   ! Determine maximum mesh index for each pressure zone
   DO IZ=1,N_ZONE
      INTEGER_BUFFER(IZ) = P_ZONE(IZ)%MESH_INDEX
   ENDDO
   CALL MPI_ALLREDUCE(MPI_IN_PLACE,INTEGER_BUFFER(1),N_ZONE,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERR)
   ! Put data in linear arrays
   DO IZ=1,N_ZONE
      REAL_BUFFER(IZ) = P_ZONE(IZ)%VOLUME
      IF (PROCESS(INTEGER_BUFFER(IZ))==MY_RANK) THEN  ! Use cell indices in the mesh with largest index
         INTEGER_BUFFER_4(4*(IZ-1)+1) = P_ZONE(IZ)%MESH_INDEX
         INTEGER_BUFFER_4(4*(IZ-1)+2) = P_ZONE(IZ)%CELL_INDICES(1)
         INTEGER_BUFFER_4(4*(IZ-1)+3) = P_ZONE(IZ)%CELL_INDICES(2)
         INTEGER_BUFFER_4(4*(IZ-1)+4) = P_ZONE(IZ)%CELL_INDICES(3)
      ENDIF
      INTEGER_BUFFER(IZ) = P_ZONE(IZ)%N_CELLS
   ENDDO
   CALL MPI_ALLREDUCE(MPI_IN_PLACE,REAL_BUFFER(1),N_ZONE,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
   CALL MPI_ALLREDUCE(MPI_IN_PLACE,INTEGER_BUFFER(1),N_ZONE,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERR)
   CALL MPI_ALLREDUCE(MPI_IN_PLACE,INTEGER_BUFFER_4(1),4*N_ZONE,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERR)
   DO IZ=1,N_ZONE
      P_ZONE(IZ)%VOLUME  = REAL_BUFFER(IZ)
      P_ZONE(IZ)%N_CELLS = INTEGER_BUFFER(IZ)
      P_ZONE(IZ)%MESH_INDEX      = INTEGER_BUFFER_4(4*(IZ-1)+1)
      P_ZONE(IZ)%CELL_INDICES(1) = INTEGER_BUFFER_4(4*(IZ-1)+2)
      P_ZONE(IZ)%CELL_INDICES(2) = INTEGER_BUFFER_4(4*(IZ-1)+3)
      P_ZONE(IZ)%CELL_INDICES(3) = INTEGER_BUFFER_4(4*(IZ-1)+4)
   ENDDO
   DEALLOCATE(REAL_BUFFER)
   DEALLOCATE(INTEGER_BUFFER)
   DEALLOCATE(INTEGER_BUFFER_4)
ENDIF

! Fill ZONEs that are smaller than a user-specified MINIMUM_ZONE_VOLUME

FILL_ZONES: IF (MINIMUM_ZONE_VOLUME>TWO_EPSILON_EB) THEN

   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      M => MESHES(NM)
      DO K=1,M%KBAR
         DO J=1,M%JBAR
            DO I=1,M%IBAR
               IZ = M%PRESSURE_ZONE(I,J,K)
               IF (IZ==0) CYCLE
               IF (P_ZONE(IZ)%VOLUME<MINIMUM_ZONE_VOLUME) THEN
                  IC = M%CELL_INDEX(I,J,K)
                  IF (IC==0) THEN
                     IC = CELL_COUNT(NM)+1
                     CALL REALLOCATE_CELL(NM,IC-1,IC)
                     M%CELL(IC)%I = I
                     M%CELL(IC)%J = J
                     M%CELL(IC)%K = K
                  ENDIF
                  M%CELL(IC)%SOLID = .TRUE.
                  M%PRESSURE_ZONE(I,J,K) = 0
                  DO IOR=-3,3
                     IF (IOR==0) CYCLE
                     IW = M%CELL(IC)%WALL_INDEX(IOR)
                     M%WALL(IW)%BOUNDARY_TYPE = NULL_BOUNDARY
                     IF (IW>0 .AND. IW<=M%N_EXTERNAL_WALL_CELLS) THEN
                        BC => M%BOUNDARY_COORD(M%WALL(IW)%BC_INDEX)
                        M%CELL(M%CELL_INDEX(BC%II,BC%JJ,BC%KK))%SOLID = .TRUE.
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO

   CALL EXCHANGE_GEOMETRY_INFO

   ! Shorten the list of pressure zones

   N_ZONE_ORIG = N_ZONE
   ALLOCATE(NEW_ZONE_INDEX(N_ZONE_ORIG))
   DO IZO=1,N_ZONE_ORIG
      NEW_ZONE_INDEX(IZO) = IZO  ! The original indices of the N_ZONE_ORIG pressure zones
   ENDDO
   IZ = 0
   IZO = 0
   DO
      IZ  = IZ  + 1
      IZO = IZO + 1
      IF (IZ>N_ZONE) EXIT
      IF (P_ZONE(IZ)%VOLUME<MINIMUM_ZONE_VOLUME) THEN
         NEW_ZONE_INDEX(IZO) = 0
         NEW_ZONE_INDEX(IZO+1:N_ZONE_ORIG) = NEW_ZONE_INDEX(IZO+1:N_ZONE_ORIG) - 1
         DO IZZ=IZ,N_ZONE-1
            P_ZONE(IZZ) = P_ZONE(IZZ+1)
         ENDDO
         N_ZONE = N_ZONE - 1
         IZ = IZ - 1
      ENDIF
   ENDDO

   ! Renumber of the ZONEs to account for those that have been removed

   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      M => MESHES(NM)
      DO IZO=1,N_ZONE_ORIG
         WHERE(M%PRESSURE_ZONE==IZO) M%PRESSURE_ZONE = NEW_ZONE_INDEX(IZO)
      ENDDO
   ENDDO

   DEALLOCATE(NEW_ZONE_INDEX)

ENDIF FILL_ZONES

! Warn if the number of pressure zones is too high

IF (N_ZONE>100 .AND. MY_RANK==0) THEN
   WRITE(LU_ERR,'(A,I0,A)') 'WARNING: Number of pressure ZONEs ',N_ZONE,' is large. Consider reducing the number.'
ENDIF

END SUBROUTINE INITIALIZE_PRESSURE_ZONES


!> \brief Propagate PRESSURE_ZONE values from the mesh of origin to neighboring meshes

SUBROUTINE ZONE_BOUNDARY_EXCHANGE

USE GEOMETRY_FUNCTIONS, ONLY: ASSIGN_PRESSURE_ZONE
INTEGER :: N_OVERLAP,IW,II,JJ,KK,IIG,JJG,KKG,IIO,JJO,KKO,NOM,NM
TYPE (MESH_TYPE), POINTER :: OM
TYPE (WALL_TYPE), POINTER :: WC
TYPE (EXTERNAL_WALL_TYPE), POINTER :: EWC
TYPE (BOUNDARY_COORD_TYPE), POINTER :: BC

SETUP_PRESSURE_ZONES_INDEX = 0  ! Flag to indicate that the PRESSURE_ZONE values no longer need to be spread

DO WHILE (ANY(SETUP_PRESSURE_ZONES_INDEX==0))

   ! Broadcast the PRESSURE_ZONE array of each mesh to the MPI processes controlling neighboring meshes

   IF (N_MPI_PROCESSES>1) THEN
      DO NM=1,NMESHES
         IF (MPI_COMM_NEIGHBORS(NM)==MPI_COMM_NULL) CYCLE
         M => MESHES(NM)
         CALL MPI_BCAST(M%PRESSURE_ZONE(0,0,0),SIZE(M%PRESSURE_ZONE),MPI_INTEGER,&
                        MPI_COMM_NEIGHBORS_ROOT(NM),MPI_COMM_NEIGHBORS(NM),IERR)
      ENDDO
   ENDIF

   ! For each mesh, go around the exterior boundary looking for OPEN or INTERPOLATED boundaries.
   ! For the OPEN boundaries, propagate a PRESSURE_ZONE value of 0 into the interior of the mesh.
   ! For the INTERPOLATED boundaries (NOM/=0), get the PRESSURE_ZONE value from the neighboring mesh and propagate it internally.

   MESH_LOOP: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

      M => MESHES(NM)

      SETUP_PRESSURE_ZONES_INDEX(NM) = 1  ! When this flag remains 1 for all meshes, stop iterating

      WALL_LOOP: DO IW=1,M%N_EXTERNAL_WALL_CELLS
         WC  => M%WALL(IW)
         EWC => M%EXTERNAL_WALL(IW)
         NOM = EWC%NOM
         BC  => M%BOUNDARY_COORD(WC%BC_INDEX)
         IIG = BC%IIG
         JJG = BC%JJG
         KKG = BC%KKG
         II  = BC%II
         JJ  = BC%JJ
         KK  = BC%KK
         M%PRESSURE_ZONE(II,JJ,KK) = 0
         IF (WC%BOUNDARY_TYPE==OPEN_BOUNDARY .AND. M%PRESSURE_ZONE(IIG,JJG,KKG)<0) THEN
            CALL ASSIGN_PRESSURE_ZONE(NM,IIG,JJG,KKG,0,N_OVERLAP)
            IF (N_OVERLAP>=0) THEN
               WRITE(LU_ERR,'(A,I2,A,I2,A,I4)') 'ERROR(872): ZONE ',0,' overlaps ZONE ',N_OVERLAP,' in MESH ',NM
               STOP_STATUS = SETUP_STOP
               EXIT MESH_LOOP
            ENDIF
            SETUP_PRESSURE_ZONES_INDEX(NM) = 0
         ENDIF
         IF (NOM>0) THEN
            IIO = EWC%IIO_MIN
            JJO = EWC%JJO_MIN
            KKO = EWC%KKO_MIN
            OM => MESHES(NOM)
            IF (WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .AND. M%PRESSURE_ZONE(IIG,JJG,KKG)>=0 .AND. &
               OM%PRESSURE_ZONE(IIO,JJO,KKO)>=0 .AND. M%PRESSURE_ZONE(IIG,JJG,KKG)/=OM%PRESSURE_ZONE(IIO,JJO,KKO)) THEN
               WRITE(LU_ERR,'(10(A,I0),A)') 'ERROR(873): ZONE ',OM%PRESSURE_ZONE(IIO,JJO,KKO),' meets ZONE ',&
                  M%PRESSURE_ZONE(IIG,JJG,KKG),' at the boundary of MESH ',NOM,' (',IIO,',',JJO,',',KKO,') and MESH ',&
                  NM,' (',IIG,',',JJG,',',KKG,')'
               STOP_STATUS = SETUP_STOP
               EXIT MESH_LOOP
            ENDIF
            IF (WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .AND. M%PRESSURE_ZONE(IIG,JJG,KKG)<0 .AND. &
                OM%PRESSURE_ZONE(IIO,JJO,KKO)>=0) THEN
               CALL ASSIGN_PRESSURE_ZONE(NM,IIG,JJG,KKG,OM%PRESSURE_ZONE(IIO,JJO,KKO),N_OVERLAP)
               IF (N_OVERLAP>0) THEN
                  WRITE(LU_ERR,'(A,I0,A,I0,A,I0)') &
                     'ERROR(872): ZONE ',OM%PRESSURE_ZONE(IIO,JJO,KKO),' overlaps ZONE ',N_OVERLAP,' in MESH ',NM
                  STOP_STATUS = SETUP_STOP
                  EXIT MESH_LOOP
               ENDIF
               SETUP_PRESSURE_ZONES_INDEX(NM) = 0
            ENDIF
            M%PRESSURE_ZONE(II,JJ,KK) = OM%PRESSURE_ZONE(IIO,JJO,KKO)
         ENDIF
      ENDDO WALL_LOOP

      M%PRESSURE_ZONE(0:M%IBP1,       0,       0) = 0
      M%PRESSURE_ZONE(0:M%IBP1,  M%JBP1,       0) = 0
      M%PRESSURE_ZONE(0:M%IBP1,  M%JBP1,  M%KBP1) = 0
      M%PRESSURE_ZONE(0:M%IBP1,       0,  M%KBP1) = 0
      M%PRESSURE_ZONE(       0,0:M%JBP1,       0) = 0
      M%PRESSURE_ZONE(  M%IBP1,0:M%JBP1,       0) = 0
      M%PRESSURE_ZONE(  M%IBP1,0:M%JBP1,  M%KBP1) = 0
      M%PRESSURE_ZONE(       0,0:M%JBP1,  M%KBP1) = 0
      M%PRESSURE_ZONE(       0,       0,0:M%KBP1) = 0
      M%PRESSURE_ZONE(  M%IBP1,       0,0:M%KBP1) = 0
      M%PRESSURE_ZONE(  M%IBP1,  M%JBP1,0:M%KBP1) = 0
      M%PRESSURE_ZONE(       0,  M%JBP1,0:M%KBP1) = 0
   ENDDO MESH_LOOP

   CALL MPI_ALLREDUCE(MPI_IN_PLACE,STOP_STATUS,INTEGER_ONE,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERR)
   IF (STOP_STATUS/=0) RETURN

   CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,SETUP_PRESSURE_ZONES_INDEX(1:NMESHES),&
                       COUNTS,DISPLS,MPI_INTEGER,MPI_COMM_WORLD,IERR)
ENDDO

END SUBROUTINE ZONE_BOUNDARY_EXCHANGE


!> \brief Set up receive buffers for MPI calls.
!> \details This subroutine loops over meshes (NM) that are to receive MPI sends from neighboring meshes (NOM).
!> When CODE=0, initializations are done, in particular the set up of persistent receives.

SUBROUTINE POST_RECEIVES(CODE)

INTEGER, INTENT(IN) :: CODE
INTEGER :: RNODE,SNODE,IJK_SIZE,N,NRA,NRA_MAX,LL,AIC,NN,NOM
REAL(EB) :: TNOW
TYPE (STORAGE_TYPE), POINTER :: OS

TNOW = CURRENT_TIME()

RECEIVING_MESH_LOOP: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   RNODE = PROCESS(NM)
   M => MESHES(NM)

   SENDING_MESH_LOOP: DO NN=1,M%N_NEIGHBORING_MESHES

      NOM = M%NEIGHBORING_MESH(NN)
      SNODE = PROCESS(NOM)
      IF (RNODE==SNODE) CYCLE SENDING_MESH_LOOP

      M3=>MESHES(NM)%OMESH(NOM)

      ! Set up receives for one-time exchanges or persistent send/receives.

      INITIALIZATION_IF: IF (CODE==0) THEN

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

            IJK_SIZE = (M3%I_MAX_R-M3%I_MIN_R+1)*(M3%J_MAX_R-M3%J_MIN_R+1)*(M3%K_MAX_R-M3%K_MIN_R+1)
            ALLOCATE(M3%REAL_RECV_PKG1(M3%NIC_R*(6+2*N_TOTAL_SCALARS)))
            ALLOCATE(M3%REAL_RECV_PKG3(IJK_SIZE*4))
            ALLOCATE(M3%REAL_RECV_PKG5(NRA_MAX*NUMBER_SPECTRAL_BANDS*M3%NIC_R))
            ALLOCATE(M3%REAL_RECV_PKG7(M3%NIC_R*3))
            ALLOCATE(M3%REAL_RECV_PKG8(M3%NIC_R*2))

            IF (LEVEL_SET_MODE>0 .OR. TERRAIN_CASE) ALLOCATE(M3%REAL_RECV_PKG14(5*M3%NIC_R))

         ENDIF

         ! Set up persistent receive requests

         IF (MPI_PARTICLE_EXCHANGE .AND. RNODE/=SNODE) THEN
            OS => M3%PARTICLE_RECV_BUFFER
            N_REQ2=N_REQ2+1 ; CALL MPI_RECV_INIT(OS%N_ITEMS   ,1,MPI_INTEGER,SNODE,NOM,MPI_COMM_WORLD,REQ2(N_REQ2),IERR)
            N_REQ2=N_REQ2+1 ; CALL MPI_RECV_INIT(OS%N_REALS   ,1,MPI_INTEGER,SNODE,NOM,MPI_COMM_WORLD,REQ2(N_REQ2),IERR)
            N_REQ2=N_REQ2+1 ; CALL MPI_RECV_INIT(OS%N_INTEGERS,1,MPI_INTEGER,SNODE,NOM,MPI_COMM_WORLD,REQ2(N_REQ2),IERR)
            N_REQ2=N_REQ2+1 ; CALL MPI_RECV_INIT(OS%N_LOGICALS,1,MPI_INTEGER,SNODE,NOM,MPI_COMM_WORLD,REQ2(N_REQ2),IERR)
         ENDIF

         IF (M3%NIC_R>0) THEN

            N_REQ1 = N_REQ1 + 1
            CALL MPI_RECV_INIT(M3%REAL_RECV_PKG1(1),SIZE(M3%REAL_RECV_PKG1),MPI_DOUBLE_PRECISION,SNODE,NOM,MPI_COMM_WORLD,&
                               REQ1(N_REQ1),IERR)

            IF (EXCHANGE_OBST_MASS) THEN
               N_REQ15 = N_REQ15 + 1
               CALL MPI_RECV_INIT(M3%N_INTERNAL_OBST,INTEGER_ONE,MPI_INTEGER,SNODE,NOM,MPI_COMM_WORLD,REQ15(N_REQ15),IERR)
            ENDIF

            N_REQ3 = N_REQ3 + 1
            CALL MPI_RECV_INIT(M3%REAL_RECV_PKG3(1),SIZE(M3%REAL_RECV_PKG3),MPI_DOUBLE_PRECISION,SNODE,NOM,MPI_COMM_WORLD,&
                               REQ3(N_REQ3),IERR)

            N_REQ7 = N_REQ7 + 1
            CALL MPI_RECV_INIT(M3%REAL_RECV_PKG7(1),SIZE(M3%REAL_RECV_PKG7),MPI_DOUBLE_PRECISION,SNODE,NOM,MPI_COMM_WORLD,&
                               REQ7(N_REQ7),IERR)

            IF (RADIATION) THEN
               N_REQ5 = N_REQ5 + 1
               CALL MPI_RECV_INIT(M3%REAL_RECV_PKG5(1),SIZE(M3%REAL_RECV_PKG5),MPI_DOUBLE_PRECISION,SNODE,NOM,MPI_COMM_WORLD,&
                                  REQ5(N_REQ5),IERR)
            ENDIF

            IF (LEVEL_SET_MODE>0 .OR. TERRAIN_CASE) THEN
               N_REQ14 = N_REQ14 + 1
               CALL MPI_RECV_INIT(M3%REAL_RECV_PKG14(1),SIZE(M3%REAL_RECV_PKG14),MPI_DOUBLE_PRECISION,SNODE,NOM,MPI_COMM_WORLD,&
                                  REQ14(N_REQ14),IERR)
            ENDIF

         ENDIF

      ENDIF INITIALIZATION_IF

      ! Exchange BACK_WALL information

      IF (CODE==8) THEN
         OS => M3%WALL_SEND_BUFFER
         N_REQ=N_REQ+1 ; CALL MPI_IRECV(OS%N_ITEMS    ,1,MPI_INTEGER,SNODE,NOM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
         N_REQ=N_REQ+1 ; CALL MPI_IRECV(OS%N_ITEMS_DIM,1,MPI_INTEGER,SNODE,NOM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
         OS => M3%THIN_WALL_SEND_BUFFER
         N_REQ=N_REQ+1 ; CALL MPI_IRECV(OS%N_ITEMS    ,1,MPI_INTEGER,SNODE,NOM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
         N_REQ=N_REQ+1 ; CALL MPI_IRECV(OS%N_ITEMS_DIM,1,MPI_INTEGER,SNODE,NOM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
      ENDIF

      OS => M3%WALL_SEND_BUFFER
      IF (CODE==9 .AND. OS%N_ITEMS>0) THEN
         N_REQ=N_REQ+1 ; CALL MPI_IRECV(OS%ITEM_INDEX,SIZE(OS%ITEM_INDEX),MPI_INTEGER,SNODE,NOM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
         N_REQ=N_REQ+1 ; CALL MPI_IRECV(OS%SURF_INDEX,SIZE(OS%SURF_INDEX),MPI_INTEGER,SNODE,NOM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
      ENDIF

      OS => M3%THIN_WALL_SEND_BUFFER
      IF (CODE==9 .AND. OS%N_ITEMS>0) THEN
         N_REQ=N_REQ+1 ; CALL MPI_IRECV(OS%ITEM_INDEX,SIZE(OS%ITEM_INDEX),MPI_INTEGER,SNODE,NOM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
         N_REQ=N_REQ+1 ; CALL MPI_IRECV(OS%SURF_INDEX,SIZE(OS%SURF_INDEX),MPI_INTEGER,SNODE,NOM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
      ENDIF

      OS => M3%WALL_RECV_BUFFER
      IF (CODE==19 .AND. OS%N_ITEMS>0) THEN
         N_REQ=N_REQ+1 ; CALL MPI_IRECV(OS%N_REALS_DIM   ,1,MPI_INTEGER,SNODE,NOM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
         N_REQ=N_REQ+1 ; CALL MPI_IRECV(OS%N_INTEGERS_DIM,1,MPI_INTEGER,SNODE,NOM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
         N_REQ=N_REQ+1 ; CALL MPI_IRECV(OS%N_LOGICALS_DIM,1,MPI_INTEGER,SNODE,NOM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
      ENDIF

      OS => M3%THIN_WALL_RECV_BUFFER
      IF (CODE==19 .AND. OS%N_ITEMS>0) THEN
         N_REQ=N_REQ+1 ; CALL MPI_IRECV(OS%N_REALS_DIM   ,1,MPI_INTEGER,SNODE,NOM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
         N_REQ=N_REQ+1 ; CALL MPI_IRECV(OS%N_INTEGERS_DIM,1,MPI_INTEGER,SNODE,NOM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
         N_REQ=N_REQ+1 ; CALL MPI_IRECV(OS%N_LOGICALS_DIM,1,MPI_INTEGER,SNODE,NOM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
      ENDIF

      OS => M3%WALL_RECV_BUFFER
      IF (CODE==10 .AND. OS%N_ITEMS>0) THEN
         N_REQ6=N_REQ6+1 ; CALL MPI_RECV_INIT(OS%REALS,OS%N_REALS_DIM,MPI_DOUBLE_PRECISION,SNODE,NOM,&
                                              MPI_COMM_WORLD,REQ6(N_REQ6),IERR)
         N_REQ6=N_REQ6+1 ; CALL MPI_RECV_INIT(OS%INTEGERS,OS%N_INTEGERS_DIM,MPI_INTEGER,SNODE,NOM,&
                                              MPI_COMM_WORLD,REQ6(N_REQ6),IERR)
         N_REQ6=N_REQ6+1 ; CALL MPI_RECV_INIT(OS%LOGICALS,OS%N_LOGICALS_DIM,MPI_LOGICAL,SNODE,NOM,&
                                              MPI_COMM_WORLD,REQ6(N_REQ6),IERR)
      ENDIF

      OS => M3%THIN_WALL_RECV_BUFFER
      IF (CODE==10 .AND. OS%N_ITEMS>0) THEN
         N_REQ6=N_REQ6+1 ; CALL MPI_RECV_INIT(OS%REALS,OS%N_REALS_DIM,MPI_DOUBLE_PRECISION,SNODE,NOM,&
                                              MPI_COMM_WORLD,REQ6(N_REQ6),IERR)
         N_REQ6=N_REQ6+1 ; CALL MPI_RECV_INIT(OS%INTEGERS,OS%N_INTEGERS_DIM,MPI_INTEGER,SNODE,NOM,&
                                              MPI_COMM_WORLD,REQ6(N_REQ6),IERR)
         N_REQ6=N_REQ6+1 ; CALL MPI_RECV_INIT(OS%LOGICALS,OS%N_LOGICALS_DIM,MPI_LOGICAL,SNODE,NOM,&
                                              MPI_COMM_WORLD,REQ6(N_REQ6),IERR)
      ENDIF

      ! PARTICLEs

      OS => M3%PARTICLE_RECV_BUFFER
      IF (CODE==11 .AND. OS%N_ITEMS>0) THEN
         N_REQ=N_REQ+1 ; CALL MPI_IRECV(OS%REALS,OS%N_REALS,MPI_DOUBLE_PRECISION,SNODE,NOM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
         N_REQ=N_REQ+1 ; CALL MPI_IRECV(OS%INTEGERS,OS%N_INTEGERS,MPI_INTEGER,SNODE,NOM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
         N_REQ=N_REQ+1 ; CALL MPI_IRECV(OS%LOGICALS,OS%N_LOGICALS,MPI_LOGICAL,SNODE,NOM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
      ENDIF

      ! Set up to receive mass loss data for boundary OBSTs

      IF (CODE==16 .AND. M3%N_INTERNAL_OBST>0) THEN
         N_REQ = N_REQ + 1
         CALL MPI_IRECV(M3%REAL_RECV_PKG8(1),2*M3%N_INTERNAL_OBST,MPI_DOUBLE_PRECISION,SNODE,NOM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
      ENDIF
      IF (CODE==18 .AND. M3%N_EXTERNAL_OBST>0) THEN
         N_REQ = N_REQ + 1
         CALL MPI_IRECV(M3%REAL_RECV_PKG8(1),2*M3%N_EXTERNAL_OBST,MPI_DOUBLE_PRECISION,SNODE,NOM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
      ENDIF

   ENDDO SENDING_MESH_LOOP

ENDDO RECEIVING_MESH_LOOP

T_USED(11)=T_USED(11) + CURRENT_TIME() - TNOW
END SUBROUTINE POST_RECEIVES


!> \brief Send and receive the major MPI packages
!> \details For each mesh, NM, controlled by MPI process, SNODE, send data to other meshes, NOM.

SUBROUTINE MESH_EXCHANGE(CODE)

REAL(EB) :: TNOW,XI,YJ,ZK
INTEGER, INTENT(IN) :: CODE
INTEGER :: NM,NOM,I,II,JJ,KK,LL,LLL,N,RNODE,SNODE,IMIN,IMAX,JMIN,JMAX,KMIN,KMAX,IJK_SIZE,RC,IC,LC
INTEGER :: NNN,NN1,NN2,IPC,II1,II2,JJ1,JJ2,KK1,KK2,NQT2,NN,IOR,NRA,NRA_MAX,AIC,OBST_INDEX,IP,IW
CHARACTER(50) :: ERR_STRING
REAL(EB), POINTER, DIMENSION(:,:) :: PHI_LS_P
REAL(EB), POINTER, DIMENSION(:,:,:) :: HP,HP2,RHOP,RHOP2,DP,DP2,UP,UP2,VP,VP2,WP,WP2
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZP,ZZP2
TYPE (LAGRANGIAN_PARTICLE_TYPE), POINTER :: LP
TYPE (BOUNDARY_COORD_TYPE), POINTER :: BC
TYPE (STORAGE_TYPE), POINTER :: OS,OS2
TYPE (WALL_TYPE), POINTER :: WC
TYPE (THIN_WALL_TYPE), POINTER :: TW

IF (CC_IBM) CALL MESH_CC_EXCHANGE(CODE)

TNOW = CURRENT_TIME()

SENDING_MESH_LOOP: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   M =>MESHES(NM)
   SNODE = PROCESS(NM)

   RECEIVING_MESH_LOOP: DO NNN=1,M%N_NEIGHBORING_MESHES

      NOM = M%NEIGHBORING_MESH(NNN)
      M3 => M%OMESH(NOM)
      IF (NOM==NM .AND. M3%NIC_S==0) CYCLE RECEIVING_MESH_LOOP

      RNODE = PROCESS(NOM)

      IF (RNODE==SNODE) THEN
         IMIN = M3%I_MIN_S
         IMAX = M3%I_MAX_S
         JMIN = M3%J_MIN_S
         JMAX = M3%J_MAX_S
         KMIN = M3%K_MIN_S
         KMAX = M3%K_MAX_S
      ENDIF

      ! Set up sends for one-time exchanges or persistent send/receives.

      INITIALIZE_SEND_IF: IF (CODE==0) THEN

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

            IJK_SIZE = (M3%I_MAX_S-M3%I_MIN_S+1)*(M3%J_MAX_S-M3%J_MIN_S+1)*(M3%K_MAX_S-M3%K_MIN_S+1)
            ALLOCATE(M3%REAL_SEND_PKG1(M3%NIC_S*(6+2*N_TOTAL_SCALARS)))
            ALLOCATE(M3%REAL_SEND_PKG3(IJK_SIZE*4))
            ALLOCATE(M3%REAL_SEND_PKG5(NRA_MAX*NUMBER_SPECTRAL_BANDS*M3%NIC_S))
            ALLOCATE(M3%REAL_SEND_PKG7(M3%NIC_S*3))

            IF (LEVEL_SET_MODE>0 .OR. TERRAIN_CASE) ALLOCATE(M3%REAL_SEND_PKG14(5*M3%NIC_S))

         ENDIF

         ! REAL_SEND_PKG8 carries obstruction mass loss information and is used even if RNODE=SNODE

         IF (M3%NIC_S>0) ALLOCATE(M3%REAL_SEND_PKG8(M3%NIC_S*2))

         ! Initialize persistent send requests

         IF (MPI_PARTICLE_EXCHANGE .AND. RNODE/=SNODE) THEN
            OS => M3%PARTICLE_SEND_BUFFER
            N_REQ2=N_REQ2+1 ; CALL MPI_SEND_INIT(OS%N_ITEMS   ,1,MPI_INTEGER,RNODE,NM,MPI_COMM_WORLD,REQ2(N_REQ2),IERR)
            N_REQ2=N_REQ2+1 ; CALL MPI_SEND_INIT(OS%N_REALS   ,1,MPI_INTEGER,RNODE,NM,MPI_COMM_WORLD,REQ2(N_REQ2),IERR)
            N_REQ2=N_REQ2+1 ; CALL MPI_SEND_INIT(OS%N_INTEGERS,1,MPI_INTEGER,RNODE,NM,MPI_COMM_WORLD,REQ2(N_REQ2),IERR)
            N_REQ2=N_REQ2+1 ; CALL MPI_SEND_INIT(OS%N_LOGICALS,1,MPI_INTEGER,RNODE,NM,MPI_COMM_WORLD,REQ2(N_REQ2),IERR)
         ENDIF

         IF (M3%NIC_S>0 .AND. RNODE/=SNODE) THEN

            N_REQ1 = N_REQ1 + 1
            CALL MPI_SEND_INIT(M3%REAL_SEND_PKG1(1),SIZE(M3%REAL_SEND_PKG1),MPI_DOUBLE_PRECISION,RNODE,NM,MPI_COMM_WORLD,&
                               REQ1(N_REQ1),IERR)

            IF (EXCHANGE_OBST_MASS) THEN
               N_REQ15 = N_REQ15 + 1
               CALL MPI_SEND_INIT(M3%N_EXTERNAL_OBST,INTEGER_ONE,MPI_INTEGER,RNODE,NM,MPI_COMM_WORLD,REQ15(N_REQ15),IERR)
            ENDIF

            N_REQ3 = N_REQ3 + 1
            CALL MPI_SEND_INIT(M3%REAL_SEND_PKG3(1),SIZE(M3%REAL_SEND_PKG3),MPI_DOUBLE_PRECISION,RNODE,NM,MPI_COMM_WORLD,&
                               REQ3(N_REQ3),IERR)

            N_REQ7 = N_REQ7 + 1
            CALL MPI_SEND_INIT(M3%REAL_SEND_PKG7(1),SIZE(M3%REAL_SEND_PKG7),MPI_DOUBLE_PRECISION,RNODE,NM,MPI_COMM_WORLD,&
                               REQ7(N_REQ7),IERR)

            IF (RADIATION) THEN
               N_REQ5 = N_REQ5 + 1
               CALL MPI_SEND_INIT(M3%REAL_SEND_PKG5(1),SIZE(M3%REAL_SEND_PKG5),MPI_DOUBLE_PRECISION,RNODE,NM,MPI_COMM_WORLD,&
                                  REQ5(N_REQ5),IERR)
            ENDIF

            IF (LEVEL_SET_MODE>0 .OR. TERRAIN_CASE) THEN
               N_REQ14 = N_REQ14 + 1
               CALL MPI_SEND_INIT(M3%REAL_SEND_PKG14(1),SIZE(M3%REAL_SEND_PKG14),MPI_DOUBLE_PRECISION,RNODE,NM,MPI_COMM_WORLD,&
                                  REQ14(N_REQ14),IERR)
            ENDIF

         ENDIF

      ENDIF INITIALIZE_SEND_IF

      ! Exchange the number of solid surface cells whose back side is in another mesh

      IF (CODE==8) THEN
         IF (RNODE/=SNODE) THEN
            OS => M3%WALL_RECV_BUFFER
            N_REQ=N_REQ+1 ; CALL MPI_ISEND(OS%N_ITEMS    ,1,MPI_INTEGER,RNODE,NM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
            N_REQ=N_REQ+1 ; CALL MPI_ISEND(OS%N_ITEMS_DIM,1,MPI_INTEGER,RNODE,NM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
            OS => M3%THIN_WALL_RECV_BUFFER
            N_REQ=N_REQ+1 ; CALL MPI_ISEND(OS%N_ITEMS    ,1,MPI_INTEGER,RNODE,NM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
            N_REQ=N_REQ+1 ; CALL MPI_ISEND(OS%N_ITEMS_DIM,1,MPI_INTEGER,RNODE,NM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
         ELSE
            M2=>MESHES(NOM)%OMESH(NM)
                 M2%WALL_SEND_BUFFER%N_ITEMS     = M3%WALL_RECV_BUFFER%N_ITEMS
                 M2%WALL_SEND_BUFFER%N_ITEMS_DIM = M3%WALL_RECV_BUFFER%N_ITEMS_DIM
            M2%THIN_WALL_SEND_BUFFER%N_ITEMS     = M3%THIN_WALL_RECV_BUFFER%N_ITEMS
            M2%THIN_WALL_SEND_BUFFER%N_ITEMS_DIM = M3%THIN_WALL_RECV_BUFFER%N_ITEMS_DIM
         ENDIF
      ENDIF

      OS => M3%WALL_RECV_BUFFER
      IF (CODE==9 .AND. OS%N_ITEMS>0 .AND. RNODE/=SNODE) THEN
         N_REQ=N_REQ+1 ; CALL MPI_ISEND(OS%ITEM_INDEX,SIZE(OS%ITEM_INDEX),MPI_INTEGER,RNODE,NM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
         N_REQ=N_REQ+1 ; CALL MPI_ISEND(OS%SURF_INDEX,SIZE(OS%SURF_INDEX),MPI_INTEGER,RNODE,NM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
      ENDIF

      OS => M3%THIN_WALL_RECV_BUFFER
      IF (CODE==9 .AND. OS%N_ITEMS>0 .AND. RNODE/=SNODE) THEN
         N_REQ=N_REQ+1 ; CALL MPI_ISEND(OS%ITEM_INDEX,SIZE(OS%ITEM_INDEX),MPI_INTEGER,RNODE,NM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
         N_REQ=N_REQ+1 ; CALL MPI_ISEND(OS%SURF_INDEX,SIZE(OS%SURF_INDEX),MPI_INTEGER,RNODE,NM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
      ENDIF

      OS => M3%WALL_SEND_BUFFER
      IF (CODE==19 .AND. OS%N_ITEMS>0) THEN
         IF (RNODE/=SNODE) THEN
            N_REQ=N_REQ+1 ; CALL MPI_ISEND(OS%N_REALS_DIM   ,1,MPI_INTEGER,RNODE,NM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
            N_REQ=N_REQ+1 ; CALL MPI_ISEND(OS%N_INTEGERS_DIM,1,MPI_INTEGER,RNODE,NM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
            N_REQ=N_REQ+1 ; CALL MPI_ISEND(OS%N_LOGICALS_DIM,1,MPI_INTEGER,RNODE,NM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
         ELSE
            OS2 => MESHES(NOM)%OMESH(NM)%WALL_RECV_BUFFER
            OS2%N_REALS_DIM    = OS%N_REALS_DIM
            OS2%N_INTEGERS_DIM = OS%N_INTEGERS_DIM
            OS2%N_LOGICALS_DIM = OS%N_LOGICALS_DIM
         ENDIF
      ENDIF

      OS => M3%THIN_WALL_SEND_BUFFER
      IF (CODE==19 .AND. OS%N_ITEMS>0) THEN
         IF (RNODE/=SNODE) THEN
            N_REQ=N_REQ+1 ; CALL MPI_ISEND(OS%N_REALS_DIM   ,1,MPI_INTEGER,RNODE,NM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
            N_REQ=N_REQ+1 ; CALL MPI_ISEND(OS%N_INTEGERS_DIM,1,MPI_INTEGER,RNODE,NM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
            N_REQ=N_REQ+1 ; CALL MPI_ISEND(OS%N_LOGICALS_DIM,1,MPI_INTEGER,RNODE,NM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
         ELSE
            OS2 => MESHES(NOM)%OMESH(NM)%THIN_WALL_RECV_BUFFER
            OS2%N_REALS_DIM    = OS%N_REALS_DIM
            OS2%N_INTEGERS_DIM = OS%N_INTEGERS_DIM
            OS2%N_LOGICALS_DIM = OS%N_LOGICALS_DIM
         ENDIF
      ENDIF

      OS => M3%WALL_SEND_BUFFER
      IF (CODE==10 .AND. OS%N_ITEMS>0 .AND. RNODE/=SNODE) THEN
         N_REQ6=N_REQ6+1; CALL MPI_SEND_INIT(OS%REALS,OS%N_REALS_DIM,MPI_DOUBLE_PRECISION,RNODE,NM,MPI_COMM_WORLD,REQ6(N_REQ6),IERR)
         N_REQ6=N_REQ6+1; CALL MPI_SEND_INIT(OS%INTEGERS,OS%N_INTEGERS_DIM,MPI_INTEGER,RNODE,NM,MPI_COMM_WORLD,REQ6(N_REQ6),IERR)
         N_REQ6=N_REQ6+1; CALL MPI_SEND_INIT(OS%LOGICALS,OS%N_LOGICALS_DIM,MPI_LOGICAL,RNODE,NM,MPI_COMM_WORLD,REQ6(N_REQ6),IERR)
      ENDIF

      OS => M3%THIN_WALL_SEND_BUFFER
      IF (CODE==10 .AND. OS%N_ITEMS>0 .AND. RNODE/=SNODE) THEN
         N_REQ6=N_REQ6+1; CALL MPI_SEND_INIT(OS%REALS,OS%N_REALS_DIM,MPI_DOUBLE_PRECISION,RNODE,NM,MPI_COMM_WORLD,REQ6(N_REQ6),IERR)
         N_REQ6=N_REQ6+1; CALL MPI_SEND_INIT(OS%INTEGERS,OS%N_INTEGERS_DIM,MPI_INTEGER,RNODE,NM,MPI_COMM_WORLD,REQ6(N_REQ6),IERR)
         N_REQ6=N_REQ6+1; CALL MPI_SEND_INIT(OS%LOGICALS,OS%N_LOGICALS_DIM,MPI_LOGICAL,RNODE,NM,MPI_COMM_WORLD,REQ6(N_REQ6),IERR)
      ENDIF

      ! Exchange of density and species mass fractions following the PREDICTOR update (CODE=1) or CORRECTOR (CODE=4) update

      IF ((CODE==1.OR.CODE==4) .AND. M3%NIC_S>0) THEN
         IF (CODE==1) THEN
            RHOP => M%RHOS ; DP => M%D  ; ZZP => M%ZZS
         ELSE
            RHOP => M%RHO  ; DP => M%DS ; ZZP => M%ZZ
         ENDIF
         IF (RNODE/=SNODE) THEN
            NQT2 = 6+2*N_TOTAL_SCALARS
            PACK_REAL_SEND_PKG1: DO LL=1,M3%NIC_S
               II1 = M3%IIO_S(LL) ; II2 = II1
               JJ1 = M3%JJO_S(LL) ; JJ2 = JJ1
               KK1 = M3%KKO_S(LL) ; KK2 = KK1
               SELECT CASE(M3%IOR_S(LL))
                  CASE(-1) ; II2=II1+1
                  CASE( 1) ; II2=II1-1
                  CASE(-2) ; JJ2=JJ1+1
                  CASE( 2) ; JJ2=JJ1-1
                  CASE(-3) ; KK2=KK1+1
                  CASE( 3) ; KK2=KK1-1
               END SELECT
               M3%REAL_SEND_PKG1(NQT2*(LL-1)+1) =   RHOP(II1,JJ1,KK1)
               M3%REAL_SEND_PKG1(NQT2*(LL-1)+2) =   RHOP(II2,JJ2,KK2)
               M3%REAL_SEND_PKG1(NQT2*(LL-1)+3) =   M%MU(II1,JJ1,KK1)
               M3%REAL_SEND_PKG1(NQT2*(LL-1)+4) = M%KRES(II1,JJ1,KK1)
               M3%REAL_SEND_PKG1(NQT2*(LL-1)+5) =     DP(II1,JJ1,KK1)
               M3%REAL_SEND_PKG1(NQT2*(LL-1)+6) =    M%Q(II1,JJ1,KK1)
               DO NN=1,N_TOTAL_SCALARS
                  M3%REAL_SEND_PKG1(NQT2*(LL-1)+6+2*NN-1) = ZZP(II1,JJ1,KK1,NN)
                  M3%REAL_SEND_PKG1(NQT2*(LL-1)+6+2*NN  ) = ZZP(II2,JJ2,KK2,NN)
               ENDDO
            ENDDO PACK_REAL_SEND_PKG1
         ELSE
            M2=>MESHES(NOM)%OMESH(NM)
            IF (CODE==1) THEN
               RHOP2 => M2%RHOS ; DP2 => M2%D  ; ZZP2 => M2%ZZS
            ELSE
               RHOP2 => M2%RHO  ; DP2 => M2%DS ; ZZP2 => M2%ZZ
            ENDIF
            RHOP2(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)   = RHOP(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
            M2%MU(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)   = M%MU(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
            M2%KRES(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) = M%KRES(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
            M2%Q(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)    = M%Q(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
            DP2(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)     = DP(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
            ZZP2(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX,1:N_TOTAL_SCALARS)= ZZP(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX,1:N_TOTAL_SCALARS)
         ENDIF
      ENDIF

      ! Exchange velocity/pressure info for ITERATE_PRESSURE

      IF (CODE==5 .AND. M3%NIC_S>0) THEN
         IF (PREDICTOR) THEN
            HP => M%H
         ELSE
            HP => M%HS
         ENDIF
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
            IF (PREDICTOR) THEN
               HP2 => M2%H
            ELSE
               HP2 => M2%HS
            ENDIF
            M2%FVX(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) = M%FVX(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
            M2%FVY(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) = M%FVY(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
            M2%FVZ(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) = M%FVZ(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
            HP2(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)    = HP(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
         ENDIF
      ENDIF

      ! Send pressure information at the end of the PREDICTOR (CODE=3) or CORRECTOR (CODE=6) stage of the time step

      IF ((CODE==3.OR.CODE==6) .AND. M3%NIC_S>0) THEN
         IF (CODE==3) THEN
            HP => M%HS ; UP => M%US ; VP => M%VS ; WP => M%WS
         ELSE
            HP => M%H  ; UP => M%U  ; VP => M%V  ; WP => M%W
         ENDIF
         IF (RNODE/=SNODE) THEN
            LL = 0
            DO KK=M3%K_MIN_S,M3%K_MAX_S
               DO JJ=M3%J_MIN_S,M3%J_MAX_S
                  DO II=M3%I_MIN_S,M3%I_MAX_S
                     M3%REAL_SEND_PKG3(LL+1) = HP(II,JJ,KK)
                     M3%REAL_SEND_PKG3(LL+2) = UP(II,JJ,KK)
                     M3%REAL_SEND_PKG3(LL+3) = VP(II,JJ,KK)
                     M3%REAL_SEND_PKG3(LL+4) = WP(II,JJ,KK)
                     LL = LL+4
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            M2=>MESHES(NOM)%OMESH(NM)
            IF (CODE==3) THEN
               HP2 => M2%HS ; UP2 => M2%US ; VP2 => M2%VS ; WP2 => M2%WS
            ELSE
               HP2 => M2%H  ; UP2 => M2%U  ; VP2 => M2%V  ; WP2 => M2%W
            ENDIF
            HP2(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) = HP(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
            UP2(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) = UP(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
            VP2(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) = VP(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
            WP2(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) = WP(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
         ENDIF
      ENDIF

      ! Exchange BACK_WALL information

      OS => M3%WALL_SEND_BUFFER
      IF (CODE==6 .AND. RNODE/=SNODE .AND. OS%N_ITEMS>0) THEN
         RC=0 ; IC=0 ; LC=0
         DO I=1,OS%N_ITEMS
            IW = OS%ITEM_INDEX(I)
            WC => M%WALL(IW)
            CALL PACK_WALL(NM,OS,WC,OS%SURF_INDEX(I),RC,IC,LC,UNPACK_IT=.FALSE.,COUNT_ONLY=.FALSE.,CHECK_BOUNDS=.FALSE.)
         ENDDO
      ENDIF

      OS => M3%THIN_WALL_SEND_BUFFER
      IF (CODE==6 .AND. RNODE/=SNODE .AND. OS%N_ITEMS>0) THEN
         RC=0 ; IC=0 ; LC=0
         DO I=1,OS%N_ITEMS
            IW = OS%ITEM_INDEX(I)
            TW => M%THIN_WALL(IW)
            CALL PACK_THIN_WALL(NM,OS,TW,OS%SURF_INDEX(I),RC,IC,LC,UNPACK_IT=.FALSE.,COUNT_ONLY=.FALSE.,CHECK_BOUNDS=.FALSE.)
         ENDDO
      ENDIF

      ! Send out radiation info

      SEND_RADIATION: IF (CODE==2 .AND. M3%NIC_S>0) THEN
         IF (RNODE/=SNODE) THEN
            IF (ICYC>1) THEN
               AIC = M%ANGLE_INC_COUNTER
            ELSE
               AIC = ANG_INC_COUNTER
            ENDIF
            LLL = 0
            PACK_REAL_SEND_PKG5: DO LL=1,M3%NIC_S
               IOR = M3%IOR_S(LL)
               DO NN2=1,NUMBER_SPECTRAL_BANDS
                  DO NN1=NUMBER_RADIATION_ANGLES-AIC+1,1,-ANGLE_INCREMENT
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

      IF (CODE==7 .AND. RNODE==SNODE) THEN
         M2=>MESHES(NOM)%OMESH(NM)
         OS=>M2%PARTICLE_RECV_BUFFER
         OS%N_ITEMS      = M3%PARTICLE_SEND_BUFFER%N_ITEMS
         OS%N_REALS      = M3%PARTICLE_SEND_BUFFER%N_REALS
         OS%N_INTEGERS   = M3%PARTICLE_SEND_BUFFER%N_INTEGERS
         OS%N_LOGICALS   = M3%PARTICLE_SEND_BUFFER%N_LOGICALS
      ENDIF

      ! Sending/Receiving PARTICLE Buffer Arrays

      IF_SEND_PARTICLES: IF (CODE==11 .AND. M3%PARTICLE_SEND_BUFFER%N_ITEMS>0) THEN
         OS => M3%PARTICLE_SEND_BUFFER
         IF (SNODE/=RNODE) THEN
            N_REQ=N_REQ+1 ; CALL MPI_ISEND(OS%REALS   ,OS%N_REALS,MPI_DOUBLE_PRECISION,RNODE,NM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
            N_REQ=N_REQ+1 ; CALL MPI_ISEND(OS%INTEGERS,OS%N_INTEGERS,MPI_INTEGER,RNODE,NM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
            N_REQ=N_REQ+1 ; CALL MPI_ISEND(OS%LOGICALS,OS%N_LOGICALS,MPI_LOGICAL,RNODE,NM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
         ELSE
            M2 => MESHES(NOM)%OMESH(NM)
            M2%PARTICLE_RECV_BUFFER%REALS    = OS%REALS
            M2%PARTICLE_RECV_BUFFER%INTEGERS = OS%INTEGERS
            M2%PARTICLE_RECV_BUFFER%LOGICALS = OS%LOGICALS
         ENDIF
         OS%N_ITEMS = 0
         OS%N_REALS = 0
         OS%N_INTEGERS = 0
         OS%N_LOGICALS = 0
      ENDIF IF_SEND_PARTICLES

      ! Send LEVEL_SET boundary values

      IF (CODE==14 .AND. M3%NIC_S>0) THEN
         IF (PREDICTOR) THEN
            PHI_LS_P => M%PHI1_LS
         ELSE
            PHI_LS_P => M%PHI_LS
         ENDIF
         IF (RNODE/=SNODE) THEN
            NQT2 = 5
            PACK_REAL_SEND_PKG14: DO LL=1,M3%NIC_S
               II1 = M3%IIO_S(LL) ; II2 = II1
               JJ1 = M3%JJO_S(LL) ; JJ2 = JJ1
               SELECT CASE(M3%IOR_S(LL))
                  CASE(-1) ; II2=II1+1
                  CASE( 1) ; II2=II1-1
                  CASE(-2) ; JJ2=JJ1+1
                  CASE( 2) ; JJ2=JJ1-1
               END SELECT
               M3%REAL_SEND_PKG14(NQT2*(LL-1)+1) = PHI_LS_P(II1,JJ1)
               M3%REAL_SEND_PKG14(NQT2*(LL-1)+2) = PHI_LS_P(II2,JJ2)
               M3%REAL_SEND_PKG14(NQT2*(LL-1)+3) = M%U_LS(II1,JJ1)
               M3%REAL_SEND_PKG14(NQT2*(LL-1)+4) = M%V_LS(II1,JJ1)
               M3%REAL_SEND_PKG14(NQT2*(LL-1)+5) = M%Z_LS(II1,JJ1)
            ENDDO PACK_REAL_SEND_PKG14
         ELSE
            M2=>MESHES(NOM)%OMESH(NM)
            IF (PREDICTOR) THEN
               M2%PHI1_LS(IMIN:IMAX,JMIN:JMAX) = PHI_LS_P(IMIN:IMAX,JMIN:JMAX)
            ELSE
               M2%PHI_LS(IMIN:IMAX,JMIN:JMAX)  = PHI_LS_P(IMIN:IMAX,JMIN:JMAX)
            ENDIF
            M2%U_LS(IMIN:IMAX,JMIN:JMAX) = M%U_LS(IMIN:IMAX,JMIN:JMAX)
            M2%V_LS(IMIN:IMAX,JMIN:JMAX) = M%V_LS(IMIN:IMAX,JMIN:JMAX)
            M2%Z_LS(IMIN:IMAX,JMIN:JMAX) = M%Z_LS(IMIN:IMAX,JMIN:JMAX)
         ENDIF
      ENDIF

      ! Send mass losses from OBSTRUCTIONs that only border the current mesh

      IF (CODE==15) THEN
         IF (RNODE==SNODE) THEN
            M2=>MESHES(NOM)%OMESH(NM)
            M2%N_INTERNAL_OBST = M3%N_EXTERNAL_OBST
         ENDIF
      ENDIF

      IF (CODE==16 .AND. M3%N_EXTERNAL_OBST>0) THEN
         IF (SNODE/=RNODE) THEN
            N_REQ = N_REQ + 1
            CALL MPI_ISEND(M3%REAL_SEND_PKG8(1),2*M3%N_EXTERNAL_OBST,MPI_DOUBLE_PRECISION,RNODE,NM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
         ELSE
            M2 => MESHES(NOM)%OMESH(NM)
            M2%REAL_RECV_PKG8 = M3%REAL_SEND_PKG8
         ENDIF
      ENDIF

      IF (CODE==18 .AND. M3%N_INTERNAL_OBST>0) THEN
         IF (SNODE/=RNODE) THEN
            N_REQ = N_REQ + 1
            CALL MPI_ISEND(M3%REAL_SEND_PKG8(1),2*M3%N_INTERNAL_OBST,MPI_DOUBLE_PRECISION,RNODE,NM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
         ELSE
            M2 => MESHES(NOM)%OMESH(NM)
            M2%REAL_RECV_PKG8 = M3%REAL_SEND_PKG8
         ENDIF
      ENDIF

   ENDDO RECEIVING_MESH_LOOP

ENDDO SENDING_MESH_LOOP

! Start the communications. Note that MPI_STARTALL starts the persistent send/receives

IF (N_MPI_PROCESSES>1 .AND. N_REQ>0) THEN
   SELECT CASE(CODE)
      CASE( 8) ; ERR_STRING = 'MPI exchange of WALL information (CODE=8)'
      CASE( 9) ; ERR_STRING = 'MPI exchange of WALL information (CODE=9)'
      CASE(19) ; ERR_STRING = 'MPI exchange of WALL information (CODE=19)'
      CASE(11) ; ERR_STRING = 'MPI exchange of PARTICLE information (CODE=11)'
      CASE(16) ; ERR_STRING = 'MPI exchange of OBST information (CODE=16)'
      CASE(18) ; ERR_STRING = 'MPI exchange of OBST information (CODE=18)'
   END SELECT
   CALL TIMEOUT(ERR_STRING,N_REQ,REQ(1:N_REQ))
   N_REQ = 0
ENDIF

IF (N_MPI_PROCESSES>1 .AND. (CODE==1.OR.CODE==4) .AND. N_REQ1>0) THEN
   CALL MPI_STARTALL(N_REQ1,REQ1(1:N_REQ1),IERR)
   IF (CODE==1) THEN
      CALL TIMEOUT('MPI exchange of gas species densities (CODE=1)',N_REQ1,REQ1(1:N_REQ1))
   ELSE
      CALL TIMEOUT('MPI exchange of gas species densities (CODE=4)',N_REQ1,REQ1(1:N_REQ1))
   ENDIF
ENDIF

IF (N_MPI_PROCESSES>1 .AND. CODE==7 .AND. N_REQ2>0) THEN
   CALL MPI_STARTALL(N_REQ2,REQ2(1:N_REQ2),IERR)
   CALL TIMEOUT('MPI exchange of particle sizes (CODE=7)',N_REQ2,REQ2(1:N_REQ2))
ENDIF

IF (N_MPI_PROCESSES>1 .AND. (CODE==3.OR.CODE==6) .AND. N_REQ3>0) THEN
   CALL MPI_STARTALL(N_REQ3,REQ3(1:N_REQ3),IERR)
   IF (CODE==3) THEN
      CALL TIMEOUT('MPI exchange of velocity and pressure (CODE=3)',N_REQ3,REQ3(1:N_REQ3))
   ELSE
      CALL TIMEOUT('MPI exchange of velocity and pressure (CODE=6)',N_REQ3,REQ3(1:N_REQ3))
   ENDIF
ENDIF

IF (N_MPI_PROCESSES>1 .AND. CODE==5 .AND. N_REQ7>0) THEN
   CALL MPI_STARTALL(N_REQ7,REQ7(1:N_REQ7),IERR)
   CALL TIMEOUT('MPI exchange of pressure (CODE=5)',N_REQ7,REQ7(1:N_REQ7))
ENDIF

IF (N_MPI_PROCESSES>1 .AND. CODE==6 .AND. N_REQ6>0) THEN
   CALL MPI_STARTALL(N_REQ6,REQ6(1:N_REQ6),IERR)
   CALL TIMEOUT('MPI exchange of back wall info (CODE=6)',N_REQ6,REQ6(1:N_REQ6))
ENDIF

IF (N_MPI_PROCESSES>1 .AND. CODE==2 .AND. N_REQ5>0) THEN
   CALL MPI_STARTALL(N_REQ5,REQ5(1:N_REQ5),IERR)
   CALL TIMEOUT('MPI exchange of radiation (CODE=2)',N_REQ5,REQ5(1:N_REQ5))
ENDIF

IF (N_MPI_PROCESSES>1 .AND. CODE==14 .AND. N_REQ14>0) THEN
   CALL MPI_STARTALL(N_REQ14,REQ14(1:N_REQ14),IERR)
   CALL TIMEOUT('MPI exchange of level set values (CODE=14)',N_REQ14,REQ14(1:N_REQ14))
ENDIF

IF (N_MPI_PROCESSES>1 .AND. CODE==15 .AND. N_REQ15>0) THEN
   CALL MPI_STARTALL(N_REQ15,REQ15(1:N_REQ15),IERR)
   CALL TIMEOUT('MPI exchange of obstruction mass (CODE=15)',N_REQ15,REQ15(1:N_REQ15))
ENDIF

! For each mesh, NM, controlled by the current MPI process, RNODE, loop over all
! other meshes, NOM, and load data received into the appropriate arrays.

RECV_MESH_LOOP: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   M => MESHES(NM)
   RNODE = PROCESS(NM)

   SEND_MESH_LOOP: DO NNN=1,M%N_NEIGHBORING_MESHES

      NOM = M%NEIGHBORING_MESH(NNN)
      M2 => M%OMESH(NOM)
      IF (NOM==NM .AND. M2%NIC_R==0) CYCLE SEND_MESH_LOOP

      SNODE = PROCESS(NOM)

      ! Unpack densities and species mass fractions in the PREDICTOR (CODE=1) and CORRECTOR (CODE=4) step

      IF ((CODE==1.OR.CODE==4) .AND. M2%NIC_R>0 .AND. RNODE/=SNODE) THEN
            NQT2 = 6+2*N_TOTAL_SCALARS
            IF (CODE==1) THEN
               RHOP => M2%RHOS ; DP => M2%D  ; ZZP => M2%ZZS
            ELSE
               RHOP => M2%RHO  ; DP => M2%DS ; ZZP => M2%ZZ
            ENDIF
            UNPACK_REAL_RECV_PKG1: DO LL=1,M2%NIC_R
               II1 = M2%IIO_R(LL) ; II2 = II1
               JJ1 = M2%JJO_R(LL) ; JJ2 = JJ1
               KK1 = M2%KKO_R(LL) ; KK2 = KK1
               SELECT CASE(M2%IOR_R(LL))
                  CASE(-1) ; II2=II1+1
                  CASE( 1) ; II2=II1-1
                  CASE(-2) ; JJ2=JJ1+1
                  CASE( 2) ; JJ2=JJ1-1
                  CASE(-3) ; KK2=KK1+1
                  CASE( 3) ; KK2=KK1-1
               END SELECT
                  RHOP(II1,JJ1,KK1) = M2%REAL_RECV_PKG1(NQT2*(LL-1)+1)
                  RHOP(II2,JJ2,KK2) = M2%REAL_RECV_PKG1(NQT2*(LL-1)+2)
                 M2%MU(II1,JJ1,KK1) = M2%REAL_RECV_PKG1(NQT2*(LL-1)+3)
               M2%KRES(II1,JJ1,KK1) = M2%REAL_RECV_PKG1(NQT2*(LL-1)+4)
                    DP(II1,JJ1,KK1) = M2%REAL_RECV_PKG1(NQT2*(LL-1)+5)
                  M2%Q(II1,JJ1,KK1) = M2%REAL_RECV_PKG1(NQT2*(LL-1)+6)
               DO NN=1,N_TOTAL_SCALARS
                     ZZP(II1,JJ1,KK1,NN) = M2%REAL_RECV_PKG1(NQT2*(LL-1)+6+2*NN-1)
                     ZZP(II2,JJ2,KK2,NN) = M2%REAL_RECV_PKG1(NQT2*(LL-1)+6+2*NN  )
               ENDDO
            ENDDO UNPACK_REAL_RECV_PKG1
      ENDIF

      ! Unpack densities and species mass fractions following PREDICTOR exchange

      IF (CODE==5 .AND. M2%NIC_R>0 .AND. RNODE/=SNODE) THEN
         IF (PREDICTOR) THEN
            HP => M2%H
         ELSE
            HP => M2%HS
         ENDIF
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

      IF ((CODE==3.OR.CODE==6) .AND. M2%NIC_R>0 .AND. RNODE/=SNODE) THEN
         IF (CODE==3) THEN
            HP2 => M2%HS ; UP2 => M2%US ; VP2 => M2%VS ; WP2 => M2%WS
         ELSE
            HP2 => M2%H  ; UP2 => M2%U  ; VP2 => M2%V  ; WP2 => M2%W
         ENDIF
         LL = 0
         DO KK=M2%K_MIN_R,M2%K_MAX_R
            DO JJ=M2%J_MIN_R,M2%J_MAX_R
               DO II=M2%I_MIN_R,M2%I_MAX_R
                  HP2(II,JJ,KK) = M2%REAL_RECV_PKG3(LL+1)
                  UP2(II,JJ,KK) = M2%REAL_RECV_PKG3(LL+2)
                  VP2(II,JJ,KK) = M2%REAL_RECV_PKG3(LL+3)
                  WP2(II,JJ,KK) = M2%REAL_RECV_PKG3(LL+4)
                  LL = LL+4
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      ! Unpack radiation information at the end of the CORRECTOR stage of the time step

      RECEIVE_RADIATION: IF (CODE==2 .AND. M2%NIC_R>0 .AND. RNODE/=SNODE) THEN
         IF (ICYC>1) THEN
            AIC = M%ANGLE_INC_COUNTER
         ELSE
            AIC = ANG_INC_COUNTER
         ENDIF
         LLL = 0
         UNPACK_REAL_RECV_PKG5: DO LL=1,M2%NIC_R
            IOR = M2%IOR_R(LL)
            DO NN2=1,NUMBER_SPECTRAL_BANDS
               DO NN1=NUMBER_RADIATION_ANGLES-AIC+1,1,-ANGLE_INCREMENT
                  IF (DLN(IOR,NN1)<=0._EB) CYCLE
                  LLL = LLL + 1
                  M2%IL_R(LL,NN1,NN2) = M2%REAL_RECV_PKG5(LLL)
               ENDDO
            ENDDO
         ENDDO UNPACK_REAL_RECV_PKG5
      ENDIF RECEIVE_RADIATION

      ! Unpack back wall information at the end of the CORRECTOR stage of the time step

      OS => M2%WALL_RECV_BUFFER
      RECEIVE_BACK_WALL: IF (CODE==6 .AND. SNODE/=RNODE .AND. OS%N_ITEMS>0) THEN
         RC=0 ; IC=0 ; LC=0
         DO I=1,OS%N_ITEMS
            IW = OS%ITEM_INDEX(I)
            WC => MESHES(NOM)%WALL(IW)
            CALL PACK_WALL(NOM,OS,WC,OS%SURF_INDEX(I),RC,IC,LC,UNPACK_IT=.TRUE.,COUNT_ONLY=.FALSE.,&
                           CHECK_BOUNDS=INITIALIZATION_PHASE)
         ENDDO
      ENDIF RECEIVE_BACK_WALL

      OS => M2%THIN_WALL_RECV_BUFFER
      RECEIVE_BACK_THIN_WALL: IF (CODE==6 .AND. SNODE/=RNODE .AND. OS%N_ITEMS>0) THEN
         RC=0 ; IC=0 ; LC=0
         DO I=1,OS%N_ITEMS
            IW = OS%ITEM_INDEX(I)
            TW => MESHES(NOM)%THIN_WALL(IW)
            CALL PACK_THIN_WALL(NOM,OS,TW,OS%SURF_INDEX(I),RC,IC,LC,UNPACK_IT=.TRUE.,COUNT_ONLY=.FALSE.,&
                                CHECK_BOUNDS=INITIALIZATION_PHASE)
         ENDDO
      ENDIF RECEIVE_BACK_THIN_WALL

      ! Sending/Receiving PARTICLE Buffer Arrays

      OS => M2%PARTICLE_RECV_BUFFER
      IF (CODE==7 .AND. OS%N_ITEMS>0) THEN
         IF (OS%N_REALS>OS%N_REALS_DIM) THEN
            CALL REALLOCATE_STORAGE_ARRAYS(OS,N_REALS=OS%N_REALS_DIM,N_REALS_NEW=OS%N_REALS)
            OS%N_REALS_DIM = OS%N_REALS
         ENDIF
         IF (OS%N_INTEGERS>OS%N_INTEGERS_DIM) THEN
            CALL REALLOCATE_STORAGE_ARRAYS(OS,N_INTEGERS=OS%N_INTEGERS_DIM,N_INTEGERS_NEW=OS%N_INTEGERS)
            OS%N_INTEGERS_DIM = OS%N_INTEGERS
         ENDIF
         IF (OS%N_LOGICALS>OS%N_LOGICALS_DIM) THEN
            CALL REALLOCATE_STORAGE_ARRAYS(OS,N_LOGICALS=OS%N_LOGICALS_DIM,N_LOGICALS_NEW=OS%N_LOGICALS)
            OS%N_LOGICALS_DIM = OS%N_LOGICALS
         ENDIF
      ENDIF

      OS => M2%PARTICLE_RECV_BUFFER
      IF_RECEIVE_PARTICLES: IF (CODE==11 .AND. OS%N_ITEMS>0) THEN
         RC=0 ; IC=0 ; LC=0
         DO IP=M%NLP+1,M%NLP+OS%N_ITEMS
            IPC = OS%INTEGERS(IC+2)
            CALL ALLOCATE_STORAGE(NM,LP_INDEX=IP,LPC_INDEX=IPC,SURF_INDEX=LAGRANGIAN_PARTICLE_CLASS(IPC)%SURF_INDEX)
            LP => M%LAGRANGIAN_PARTICLE(IP)
            CALL PACK_PARTICLE(NM,OS,LP,IPC,RC,IC,LC,UNPACK_IT=.TRUE.,COUNT_ONLY=.FALSE.,CHECK_BOUNDS=.FALSE.)
            LP%WALL_INDEX  = 0  ! If the droplet was stuck to a wall, unstick it when it arrives in the new mesh
            IF(LP%CFACE_INDEX/=EXTERNAL_CFACE) LP%CFACE_INDEX = 0
            BC=>M%BOUNDARY_COORD(LP%BC_INDEX)
            CALL GET_IJK(BC%X,BC%Y,BC%Z,NM,XI,YJ,ZK,BC%IIG,BC%JJG,BC%KKG)
            BC%II=BC%IIG ; BC%JJ=BC%JJG ; BC%KK=BC%KKG
            CALL GET_RVC(NM,BC%IIG,BC%JJG,BC%KKG,LP%RVC)
            IF (LP%INIT_INDEX>0) THEN
               DO NN=1,N_DEVC
                  IF (DEVICE(NN)%INIT_ID==INITIALIZATION(LP%INIT_INDEX)%ID .AND. DEVICE(NN)%INIT_ID/='null') THEN
                     DEVICE(NN)%LP_TAG = LP%TAG
                     DEVICE(NN)%PART_CLASS_INDEX = IPC
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
         M%NLP = M%NLP + OS%N_ITEMS
         OS%N_ITEMS = 0
      ENDIF IF_RECEIVE_PARTICLES

      IF (CODE==14 .AND. M2%NIC_R>0 .AND. RNODE/=SNODE) THEN
            NQT2 = 5
            IF (PREDICTOR) THEN
               PHI_LS_P => M2%PHI1_LS
            ELSE
               PHI_LS_P => M2%PHI_LS
            ENDIF
            UNPACK_REAL_RECV_PKG14: DO LL=1,M2%NIC_R
               II1 = M2%IIO_R(LL) ; II2 = II1
               JJ1 = M2%JJO_R(LL) ; JJ2 = JJ1
               SELECT CASE(M2%IOR_R(LL))
                  CASE(-1) ; II2=II1+1
                  CASE( 1) ; II2=II1-1
                  CASE(-2) ; JJ2=JJ1+1
                  CASE( 2) ; JJ2=JJ1-1
               END SELECT
               PHI_LS_P(II1,JJ1) = M2%REAL_RECV_PKG14(NQT2*(LL-1)+1)
               PHI_LS_P(II2,JJ2) = M2%REAL_RECV_PKG14(NQT2*(LL-1)+2)
               M2%U_LS(II1,JJ1)  = M2%REAL_RECV_PKG14(NQT2*(LL-1)+3)
               M2%V_LS(II1,JJ1)  = M2%REAL_RECV_PKG14(NQT2*(LL-1)+4)
               M2%Z_LS(II1,JJ1)  = M2%REAL_RECV_PKG14(NQT2*(LL-1)+5)
            ENDDO UNPACK_REAL_RECV_PKG14
      ENDIF

      ! Unpack mass losses from OBST WALL cells in neighboring meshes

      ! CODE 16: MESH that contains the OBST unpacks M2%REAL_RECV_PKG8 which contains the mass lost by exterior cell faces.
      ! CODE 17: MESH that contains the OBST packs up its new masses into REAL_SEND_PKG8
      ! CODE 18: The MESHes containing the exterior cell faces receives the new mass of the OBST that lives in neighboring mesh.

      IF (CODE==16 .AND. M2%N_INTERNAL_OBST>0) THEN
         DO N=1,M2%N_INTERNAL_OBST
            OBST_INDEX = NINT(M2%REAL_RECV_PKG8(2*N-1))
            M%OBSTRUCTION(OBST_INDEX)%MASS = M%OBSTRUCTION(OBST_INDEX)%MASS - M2%REAL_RECV_PKG8(2*N)
         ENDDO
      ENDIF

      IF (CODE==17 .AND. M2%N_INTERNAL_OBST>0) THEN
         DO N=1,M2%N_INTERNAL_OBST
            OBST_INDEX = NINT(M2%REAL_RECV_PKG8(2*N-1))
            M2%REAL_SEND_PKG8(2*N-1) = REAL(OBST_INDEX,EB)
            M2%REAL_SEND_PKG8(2*N)   = M%OBSTRUCTION(OBST_INDEX)%MASS
         ENDDO
      ENDIF

      IF (CODE==18 .AND. M2%N_EXTERNAL_OBST>0) THEN
         M4 => MESHES(NOM)
         DO N=1,M2%N_EXTERNAL_OBST
            OBST_INDEX = NINT(M2%REAL_RECV_PKG8(2*N-1))
            M4%OBSTRUCTION(OBST_INDEX)%MASS = M2%REAL_RECV_PKG8(2*N)
         ENDDO
         M2%N_EXTERNAL_OBST = 0
         M2%REAL_RECV_PKG8 = 0._EB
      ENDIF

   ENDDO SEND_MESH_LOOP

ENDDO RECV_MESH_LOOP

T_USED(11)=T_USED(11) + CURRENT_TIME() - TNOW
END SUBROUTINE MESH_EXCHANGE


!> \brief Probe MPI communications until complete. Write error if the communication does not complete.
!> \param RNAME Name given to the array of MPI_REQUESTs
!> \param NR Number of MPI_REQUESTs
!> \param RR Array of MPI_REQUESTs

SUBROUTINE TIMEOUT(RNAME,NR,RR)

REAL(EB) :: START_TIME,WAIT_TIME
INTEGER, INTENT(IN) :: NR
TYPE (MPI_REQUEST), DIMENSION(:) :: RR
INTEGER :: ERRORCODE
LOGICAL :: FLAG
CHARACTER(*) :: RNAME

IF (.NOT.PROFILING) THEN

   ! Normally, PROFILING=F and this branch continually tests the communication and cancels the requests if too much time elapses.

   START_TIME = MPI_WTIME()
   FLAG = .FALSE.
   DO WHILE(.NOT.FLAG)
      CALL MPI_TESTALL(NR,RR(1:NR),FLAG,MPI_STATUSES_IGNORE,IERR)
      WAIT_TIME = MPI_WTIME() - START_TIME
      IF (WAIT_TIME>MPI_TIMEOUT) THEN
         WRITE(LU_ERR,'(/A,A,A,I0,A,A,A/)') 'ERROR(123): ',TRIM(RNAME),' timed out for MPI process ',MY_RANK,' running on ',&
                                            PNAME(1:PNAMELEN),'. FDS will abort.'
         ERRORCODE = 1
         CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
      ENDIF
   ENDDO

ELSE

   ! If PROFILING=T, do not do MPI_TESTALL because too many calls to this routine swamps the tracing and profiling.

   CALL MPI_WAITALL(NR,RR(1:NR),MPI_STATUSES_IGNORE,IERR)

ENDIF

END SUBROUTINE TIMEOUT


!> \brief Write out the file CHID_cpu.csv containing the CPU time used by the major subroutines for each MPI process.

SUBROUTINE DUMP_TIMERS

INTEGER, PARAMETER :: LINE_LENGTH = 5 + (N_TIMERS+1)*11
REAL(EB) :: T_USED_COPY(N_TIMERS)
CHARACTER(LEN=LINE_LENGTH) :: LINE
CHARACTER(LEN=LINE_LENGTH), DIMENSION(0:N_MPI_PROCESSES-1) :: LINE_ARRAY
CHARACTER(30) :: FRMT

! T_USED_COPY(1) is the time spent in the main routine; i.e. the time not spent in a subroutine.
! T_USED_COPY so not overwriting T_USED(1) when multiple CPU dumps happen.

T_USED_COPY(1) = CURRENT_TIME() - T_USED(1) - SUM(T_USED(2:N_TIMERS))
T_USED_COPY(2:N_TIMERS) = T_USED(2:N_TIMERS)
WRITE(FRMT,'(A,I2.2,A)') '(I5,',N_TIMERS+1,'(",",ES10.3))'
WRITE(LINE,FRMT) MY_RANK,(T_USED_COPY(I),I=1,N_TIMERS),SUM(T_USED_COPY(1:N_TIMERS))

! All MPI processes except root send their timings to the root process. The root process then writes them out to a file.

IF (MY_RANK>0) THEN
   CALL MPI_SEND(LINE,LINE_LENGTH,MPI_CHARACTER,0,MY_RANK,MPI_COMM_WORLD,IERR)
ELSE
   LINE_ARRAY(0) = LINE
   DO N=1,N_MPI_PROCESSES-1
      CALL MPI_RECV(LINE_ARRAY(N),LINE_LENGTH,MPI_CHARACTER,N,N,MPI_COMM_WORLD,STATUS,IERR)
   ENDDO
   FN_CPU = TRIM(CHID)//'_cpu.csv'
   OPEN(LU_CPU,FILE=FN_CPU,STATUS='REPLACE',FORM='FORMATTED')
   WRITE(LU_CPU,'(A)') 'Rank,MAIN,DIVG,MASS,VELO,PRES,WALL,DUMP,PART,RADI,FIRE,COMM,BLNK,HVAC,GEOM,VEGE,CHEM,Total T_USED (s)'
   DO N=0,N_MPI_PROCESSES-1
      WRITE(LU_CPU,'(A)') LINE_ARRAY(N)
   ENDDO
   CLOSE(LU_CPU)
ENDIF

END SUBROUTINE DUMP_TIMERS


!> \brief Write character strings out to the Smokeview (.smv) file

SUBROUTINE WRITE_STRINGS

INTEGER :: N,NOM,N_STRINGS_DUM
CHARACTER(:), ALLOCATABLE :: STRING_BUF
CHARACTER(MESH_STRING_LENGTH), ALLOCATABLE, DIMENSION(:) :: STRING_DUM
REAL(EB) :: TNOW
LOGICAL :: WRITE_ALL

TNOW = CURRENT_TIME()

! All meshes send their STRINGs to node 0

DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   IF (MY_RANK>0) THEN
      CALL MPI_SEND(MESHES(NM)%N_STRINGS,1,MPI_INTEGER,0,1,MPI_COMM_WORLD,IERR)
      IF (MESHES(NM)%N_STRINGS>0) THEN
        ! Pack the strings:
        ALLOCATE(CHARACTER(MESHES(NM)%N_STRINGS * MESH_STRING_LENGTH) :: &
                 STRING_BUF)
        DO N=1,MESHES(NM)%N_STRINGS
          ASSOCIATE ( MST => MESH_STRING_LENGTH )
          STRING_BUF ((N-1)*MST+1:N*MST) = MESHES(NM)%STRING(N)
          END ASSOCIATE
        END DO
        CALL MPI_SEND(STRING_BUF,MESHES(NM)%N_STRINGS*MESH_STRING_LENGTH,MPI_CHARACTER,0,NM,MPI_COMM_WORLD,IERR)
        DEALLOCATE(STRING_BUF)
      ENDIF
   ENDIF
ENDDO

! Node 0 receives the STRINGs and writes them to the .smv file

IF (MY_RANK==0) THEN
   DO N=1,MESHES(1)%N_STRINGS
      WRITE(LU_SMV,'(A)') TRIM(MESHES(1)%STRING(N))
   ENDDO
   OTHER_MESH_LOOP: DO NOM=2,NMESHES
      IF (PROCESS(NOM)>0) THEN
         CALL MPI_RECV(N_STRINGS_DUM,1,MPI_INTEGER,PROCESS(NOM),1,MPI_COMM_WORLD,STATUS,IERR)
         IF (N_STRINGS_DUM>0) THEN
            ALLOCATE(CHARACTER(LEN=N_STRINGS_DUM*MESH_STRING_LENGTH) :: STRING_BUF )
            CALL MPI_RECV(STRING_BUF,N_STRINGS_DUM*MESH_STRING_LENGTH,MPI_CHARACTER,PROCESS(NOM),NOM,MPI_COMM_WORLD,STATUS,IERR)
            ! Unpack the strings:
            ALLOCATE(STRING_DUM(N_STRINGS_DUM))
            ASSOCIATE ( MST => MESH_STRING_LENGTH )
            DO N=1,N_STRINGS_DUM
              STRING_DUM(N) = STRING_BUF((N-1)*MST+1:N*MST)
            END DO
            END ASSOCIATE
            DEALLOCATE (STRING_BUF)
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

   HVAC_IF: IF (HVAC_SOLVE) THEN
      WRITE_ALL=.FALSE.
      IF (ANY(DUCT%FAN_INDEX > 0) .OR. ANY(DUCT%AIRCOIL_INDEX > 0) .OR. ANY(DUCT%DAMPER)) THEN
         IF (ABS(T-T_BEGIN) < TWO_EPSILON_EB) WRITE_ALL = .TRUE.
         DO N=1,N_DUCTS
            IF (DUCT(N)%DAMPER .OR. DUCT(N)%FAN_INDEX > 0 .OR. DUCT(N)%AIRCOIL_INDEX > 0) THEN
               IF (DUCT(N)%DEVC_INDEX > 0) THEN
                  IF (WRITE_ALL .OR. DEVICE(DUCT(N)%DEVC_INDEX)%CURRENT_STATE .NEQV. DEVICE(DUCT(N)%DEVC_INDEX)%PRIOR_STATE) THEN
                     WRITE(LU_SMV,'(A)') 'DUCT_ACT '//TRIM(DUCT(N)%ID)
                     IF (DEVICE(DUCT(N)%DEVC_INDEX)%CURRENT_STATE) THEN
                        WRITE(LU_SMV,'(I6,F10.2,I6)') N,T_BEGIN+(T-T_BEGIN)*TIME_SHRINK_FACTOR,1
                     ELSE
                        WRITE(LU_SMV,'(I6,F10.2,I6)') N,T_BEGIN+(T-T_BEGIN)*TIME_SHRINK_FACTOR,0
                     ENDIF
                  ENDIF
               ELSEIF (DUCT(N)%CTRL_INDEX > 0) THEN
                  IF (WRITE_ALL .OR. CONTROL(DUCT(N)%CTRL_INDEX)%CURRENT_STATE .NEQV. CONTROL(DUCT(N)%CTRL_INDEX)%PRIOR_STATE) THEN
                     WRITE(LU_SMV,'(A)') 'DUCT_ACT '//TRIM(DUCT(N)%ID)
                     IF (CONTROL(DUCT(N)%CTRL_INDEX)%CURRENT_STATE) THEN
                        WRITE(LU_SMV,'(I6,F10.2,I6)') N,T_BEGIN+(T-T_BEGIN)*TIME_SHRINK_FACTOR,1
                     ELSE
                        WRITE(LU_SMV,'(I6,F10.2,I6)') N,T_BEGIN+(T-T_BEGIN)*TIME_SHRINK_FACTOR,0
                     ENDIF
                  ENDIF
               ELSE
                  IF (WRITE_ALL) THEN
                     WRITE(LU_SMV,'(A)') 'DUCT_ACT '//TRIM(DUCT(N)%ID)
                     WRITE(LU_SMV,'(I6,F10.2,I6)') N,T_BEGIN+(T-T_BEGIN)*TIME_SHRINK_FACTOR,1
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
      ENDIF
   ENDIF HVAC_IF

ENDIF

! All STRING arrays are zeroed out

DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   MESHES(NM)%N_STRINGS = 0
ENDDO

T_USED(11) = T_USED(11) + CURRENT_TIME() - TNOW
END SUBROUTINE WRITE_STRINGS


!> \brief MPI exchanges of diagnostic output

SUBROUTINE EXCHANGE_DIAGNOSTICS

INTEGER  :: NOM
REAL(EB) :: TNOW
TYPE :: MESH_STRUCT
    REAL(EB) :: DBLS(7)
    INTEGER :: INTS(16)
END TYPE MESH_STRUCT
INTEGER  :: LENGTH(2)
TYPE (MPI_DATATYPE) :: DTYPES(2),MESH_STRUCT_TYPE
INTEGER (KIND=MPI_ADDRESS_KIND) STRUCT_DISP(2)
TYPE(MESH_STRUCT) :: MESH_SEND, MESH_RECV

TNOW = CURRENT_TIME()

! MPI processes greater than 0 send diagnostic data to MPI process 0

LENGTH(1) = 7
LENGTH(2) = 16

STRUCT_DISP(1) = 0
STRUCT_DISP(2) = (STORAGE_SIZE(1.D0) / 8) * 7

DTYPES(1) = MPI_DOUBLE_PRECISION
DTYPES(2) = MPI_INTEGER

CALL MPI_TYPE_CREATE_STRUCT(2, LENGTH, STRUCT_DISP, DTYPES, MESH_STRUCT_TYPE, IERR)
CALL MPI_TYPE_COMMIT(MESH_STRUCT_TYPE, IERR)

DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   IF (MY_RANK==0) CYCLE
   MESH_SEND%DBLS(1)  = MESHES(NM)%CFL
   MESH_SEND%DBLS(2)  = MESHES(NM)%DIVMX
   MESH_SEND%DBLS(3)  = MESHES(NM)%DIVMN
   MESH_SEND%DBLS(4)  = MESHES(NM)%RESMAX
   MESH_SEND%DBLS(5)  = MESHES(NM)%POIS_PTB
   MESH_SEND%DBLS(6)  = MESHES(NM)%POIS_ERR
   MESH_SEND%DBLS(7)  = MESHES(NM)%VN
   MESH_SEND%INTS(1)  = MESHES(NM)%ICFL
   MESH_SEND%INTS(2)  = MESHES(NM)%JCFL
   MESH_SEND%INTS(3)  = MESHES(NM)%KCFL
   MESH_SEND%INTS(4)  = MESHES(NM)%IMX
   MESH_SEND%INTS(5)  = MESHES(NM)%JMX
   MESH_SEND%INTS(6)  = MESHES(NM)%KMX
   MESH_SEND%INTS(7)  = MESHES(NM)%IMN
   MESH_SEND%INTS(8)  = MESHES(NM)%JMN
   MESH_SEND%INTS(9)  = MESHES(NM)%KMN
   MESH_SEND%INTS(10)  = MESHES(NM)%IRM
   MESH_SEND%INTS(11)  = MESHES(NM)%JRM
   MESH_SEND%INTS(12)  = MESHES(NM)%KRM
   MESH_SEND%INTS(13)  = MESHES(NM)%I_VN
   MESH_SEND%INTS(14)  = MESHES(NM)%J_VN
   MESH_SEND%INTS(15)  = MESHES(NM)%K_VN
   MESH_SEND%INTS(16)  = MESHES(NM)%NLP
   CALL MPI_SSEND(MESH_SEND, 1, MESH_STRUCT_TYPE, 0, 1, MPI_COMM_WORLD, &
       IERR)
ENDDO

! Node 0 receives various values from all other nodes

DO NOM=1,NMESHES
   IF (PROCESS(NOM)==0 .OR. MY_RANK/=0) CYCLE
   CALL MPI_RECV(MESH_RECV, 1, MESH_STRUCT_TYPE, PROCESS(NOM), 1, &
       MPI_COMM_WORLD, STATUS, IERR)

   MESHES(NOM)%CFL       = MESH_RECV%DBLS(1)
   MESHES(NOM)%DIVMX     = MESH_RECV%DBLS(2)
   MESHES(NOM)%DIVMN     = MESH_RECV%DBLS(3)
   MESHES(NOM)%RESMAX    = MESH_RECV%DBLS(4)
   MESHES(NOM)%POIS_PTB  = MESH_RECV%DBLS(5)
   MESHES(NOM)%POIS_ERR  = MESH_RECV%DBLS(6)
   MESHES(NOM)%VN        = MESH_RECV%DBLS(7)
   MESHES(NOM)%ICFL      = MESH_RECV%INTS(1)
   MESHES(NOM)%JCFL      = MESH_RECV%INTS(2)
   MESHES(NOM)%KCFL      = MESH_RECV%INTS(3)
   MESHES(NOM)%IMX       = MESH_RECV%INTS(4)
   MESHES(NOM)%JMX       = MESH_RECV%INTS(5)
   MESHES(NOM)%KMX       = MESH_RECV%INTS(6)
   MESHES(NOM)%IMN       = MESH_RECV%INTS(7)
   MESHES(NOM)%JMN       = MESH_RECV%INTS(8)
   MESHES(NOM)%KMN       = MESH_RECV%INTS(9)
   MESHES(NOM)%IRM       = MESH_RECV%INTS(10)
   MESHES(NOM)%JRM       = MESH_RECV%INTS(11)
   MESHES(NOM)%KRM       = MESH_RECV%INTS(12)
   MESHES(NOM)%I_VN      = MESH_RECV%INTS(13)
   MESHES(NOM)%J_VN      = MESH_RECV%INTS(14)
   MESHES(NOM)%K_VN      = MESH_RECV%INTS(15)
   MESHES(NOM)%NLP       = MESH_RECV%INTS(16)
ENDDO

CALL MPI_TYPE_FREE(MESH_STRUCT_TYPE, IERR)

T_USED(11) = T_USED(11) + CURRENT_TIME() - TNOW
END SUBROUTINE EXCHANGE_DIAGNOSTICS


!> \brief Gather HRR, mass, and device data to node 0

SUBROUTINE EXCHANGE_GLOBAL_OUTPUTS

REAL(EB) :: TNOW
INTEGER :: NN,N,I_STATE,OP_INDEX,NM,DISP,DIM_FAC
TYPE (MPI_OP) :: MPI_OP_INDEX
TYPE(DEVICE_TYPE), POINTER :: DV
TYPE(SUBDEVICE_TYPE), POINTER :: SDV
LOGICAL :: NO_NEED_TO_RECV

TNOW = CURRENT_TIME()

DISP = DISPLS(MY_RANK)+1

! Exchange DEVICE parameters among meshes and dump out DEVICE info after first "gathering" data to node 0

EXCHANGE_DEVICE: IF (N_DEVC>0) THEN

   ! Exchange the CURRENT_STATE and PRIOR_STATE of each DEViCe

   STATE_ARRAY = .FALSE.  ! Temporary array that holds the STATE value for the devices on each node
   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      DO N=1,N_DEVC
         DV => DEVICE(N)
         IF (DV%MESH==NM) THEN
            STATE_ARRAY(N)        = DV%CURRENT_STATE
            STATE_ARRAY(N+N_DEVC) = DV%PRIOR_STATE
         ENDIF
      ENDDO
   ENDDO
   IF (N_MPI_PROCESSES>1) CALL MPI_ALLREDUCE(MPI_IN_PLACE,STATE_ARRAY,2*N_DEVC,MPI_LOGICAL,MPI_LXOR,MPI_COMM_WORLD,IERR)
   DO N=1,N_DEVC
      DV => DEVICE(N)
      DV%CURRENT_STATE = STATE_ARRAY(N)
      DV%PRIOR_STATE   = STATE_ARRAY(N+N_DEVC)
   ENDDO

   ! Dry pipe sprinkler logic

   DEVC_PIPE_OPERATING = 0
   DO N=1,N_DEVC
      IF (DEVICE(N)%PROP_INDEX > 0 .AND.  DEVICE(N)%CURRENT_STATE) THEN
         IF (PROPERTY(DEVICE(N)%PROP_INDEX)%PART_INDEX > 0) DEVC_PIPE_OPERATING(DEVICE(N)%PIPE_INDEX) = &
            DEVC_PIPE_OPERATING(DEVICE(N)%PIPE_INDEX) + 1
      ENDIF
   ENDDO

   ! Each DEViCe has 0 or more SUBDEVICEs. These SUBDEVICEs contain the values of the DEViCe on the meshes controlled
   ! by the copy of the DEViCe controlled by this MPI process. Each MPI process has a copy of every DEViCe, but only
   ! the MPI process that controls the meshes has a copy of the DEViCe for which SUBDEVICEs have been allocated.

   ! In the following loop, for OP_INDEX=1, we add together VALUE_1 and possibly VALUE_2 for all SUBDEVICEs (i.e. meshes)
   ! allocated by the copy of the DEViCe associated with the current MPI process, MY_RANK. Then we do an MPI_ALLREDUCE that
   ! adds the VALUE_1 and VALUE_2 from other SUBDEVICEs allocated by the copies of the DEViCE stored by the other MPI processes.
   ! For OP_INDEX=2 and 3, we take the MIN or MAX of all the VALUEs, along with the MINLOC or MAXLOC.

   OPERATION_LOOP: DO OP_INDEX=1,3

      IF (OP_INDEX==2 .AND. .NOT.MIN_DEVICES_EXIST) CYCLE OPERATION_LOOP
      IF (OP_INDEX==3 .AND. .NOT.MAX_DEVICES_EXIST) CYCLE OPERATION_LOOP
      SELECT CASE(OP_INDEX)
         CASE(1) ; TC_ARRAY  =  0._EB    ; MPI_OP_INDEX = MPI_SUM    ; DIM_FAC = 4
         CASE(2) ; TC2_ARRAY =  1.E10_EB ; MPI_OP_INDEX = MPI_MINLOC ; DIM_FAC = 1
         CASE(3) ; TC2_ARRAY = -1.E10_EB ; MPI_OP_INDEX = MPI_MAXLOC ; DIM_FAC = 1
      END SELECT

      DEVICE_LOOP_1: DO N=1,N_DEVC
         DV => DEVICE(N)
         IF (OP_INDEX==1 .AND. (DV%SPATIAL_STATISTIC(1:3)=='MIN' .OR. DV%SPATIAL_STATISTIC(1:3)=='MAX')) CYCLE DEVICE_LOOP_1
         IF (OP_INDEX==2 .AND.  DV%SPATIAL_STATISTIC(1:3)/='MIN') CYCLE DEVICE_LOOP_1
         IF (OP_INDEX==3 .AND.  DV%SPATIAL_STATISTIC(1:3)/='MAX') CYCLE DEVICE_LOOP_1
         DO NN=1,DV%N_SUBDEVICES
            SDV => DV%SUBDEVICE(NN)
            SELECT CASE(OP_INDEX)
               CASE(1)
                  TC_ARRAY(N)          = TC_ARRAY(N)          + SDV%VALUE_1
                  TC_ARRAY(N+N_DEVC)   = TC_ARRAY(N+N_DEVC)   + SDV%VALUE_2
                  TC_ARRAY(N+2*N_DEVC) = TC_ARRAY(N+2*N_DEVC) + SDV%VALUE_3
                  TC_ARRAY(N+3*N_DEVC) = TC_ARRAY(N+3*N_DEVC) + SDV%VALUE_4
               CASE(2)
                  IF (SDV%VALUE_1<TC2_ARRAY(1,N)) THEN
                     TC2_ARRAY(1,N) = SDV%VALUE_1
                     TC2_ARRAY(2,N) = SDV%VALUE_2
                  ENDIF
               CASE(3)
                  IF (SDV%VALUE_1>TC2_ARRAY(1,N)) THEN
                     TC2_ARRAY(1,N) = SDV%VALUE_1
                     TC2_ARRAY(2,N) = SDV%VALUE_2
                  ENDIF
            END SELECT
         ENDDO
      ENDDO DEVICE_LOOP_1

      ! Perform MPI exchanges to sum or take max/min of device values collected on meshes controlled by different processes

      IF (N_MPI_PROCESSES>1) THEN
         SELECT CASE(OP_INDEX)
            CASE(1) 
               CALL MPI_ALLREDUCE(MPI_IN_PLACE,TC_ARRAY,DIM_FAC*N_DEVC,MPI_DOUBLE_PRECISION,MPI_OP_INDEX,MPI_COMM_WORLD,IERR)
            CASE(2:3) 
               CALL MPI_ALLREDUCE(MPI_IN_PLACE,TC2_ARRAY,N_DEVC,MPI_2DOUBLE_PRECISION,MPI_OP_INDEX,MPI_COMM_WORLD,IERR)
         END SELECT
      ENDIF

      ! Put summed values from the subdevices (SDV) back into the controlling device (DV)

      DEVICE_LOOP_2: DO N=1,N_DEVC

         DV => DEVICE(N)
         IF (OP_INDEX==1 .AND. (DV%SPATIAL_STATISTIC(1:3)=='MIN' .OR. DV%SPATIAL_STATISTIC(1:3)=='MAX')) CYCLE DEVICE_LOOP_2
         IF (OP_INDEX==2 .AND.  DV%SPATIAL_STATISTIC(1:3)/='MIN') CYCLE DEVICE_LOOP_2
         IF (OP_INDEX==3 .AND.  DV%SPATIAL_STATISTIC(1:3)/='MAX') CYCLE DEVICE_LOOP_2
         IF (OP_INDEX==1) THEN
            DV%VALUE_1 = TC_ARRAY(N)
            DV%VALUE_2 = TC_ARRAY(  N_DEVC+N)
            DV%VALUE_3 = TC_ARRAY(2*N_DEVC+N)
            DV%VALUE_4 = TC_ARRAY(3*N_DEVC+N)
         ENDIF
         IF (OP_INDEX>1 .AND.  (DV%SPATIAL_STATISTIC=='MIN'.OR.DV%SPATIAL_STATISTIC=='MAX')) THEN
            DV%VALUE_1 = TC2_ARRAY(1,N)
         ENDIF

         ! Special case for MINLOC or MAXLOC

         IF (OP_INDEX>1 .AND. (DV%SPATIAL_STATISTIC(1:6)=='MINLOC'.OR.DV%SPATIAL_STATISTIC(1:6)=='MAXLOC')) THEN
            NO_NEED_TO_RECV = .FALSE.
            DO NN=1,DV%N_SUBDEVICES
               SDV => DV%SUBDEVICE(NN)
               IF (PROCESS(SDV%MESH)==MY_RANK) THEN
                  IF (SDV%MESH==NINT(TC2_ARRAY(2,N))) THEN
                     DV%VALUE_1 = SDV%VALUE_3
                     IF (MY_RANK>0) THEN
                        CALL MPI_SEND(DV%VALUE_1,1,MPI_DOUBLE_PRECISION,0,999,MPI_COMM_WORLD,IERR)
                     ELSE
                        NO_NEED_TO_RECV = .TRUE.
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
            IF (N_MPI_PROCESSES>1 .AND. MY_RANK==0 .AND. .NOT.NO_NEED_TO_RECV) &
               CALL MPI_RECV(DV%VALUE_1,1,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,STATUS,IERR)
            CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(DV%VALUE_1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
         ENDIF

      ENDDO DEVICE_LOOP_2

   ENDDO OPERATION_LOOP

ENDIF EXCHANGE_DEVICE

! Perform the temporal operations on the device outputs

CALL UPDATE_DEVICES_2(T,DT)

! Check for change in control function output of device

DEVICE_LOOP: DO N=1,N_DEVC

   DV => DEVICE(N)

   LATCHIF: IF (DV%LATCH) THEN
      IF (DV%INITIAL_STATE .EQV. DV%CURRENT_STATE) THEN
         DEVICE_DIRECTION: IF (DV%TRIP_DIRECTION > 0) THEN
            IF (DV%SMOOTHED_VALUE > DV%SETPOINT) DV%CURRENT_STATE = .NOT.DV%INITIAL_STATE
         ELSE DEVICE_DIRECTION
            IF (DV%SMOOTHED_VALUE < DV%SETPOINT) DV%CURRENT_STATE = .NOT.DV%INITIAL_STATE
         ENDIF DEVICE_DIRECTION
      ENDIF
   ELSE LATCHIF
      DEVICE_DIRECTION2: IF (DV%TRIP_DIRECTION > 0) THEN
         IF ((DV%SMOOTHED_VALUE > DV%SETPOINT) .AND. (DV%CURRENT_STATE .EQV.  DV%INITIAL_STATE)) THEN
            DV%CURRENT_STATE = .NOT.DV%INITIAL_STATE
         ELSEIF ((DV%SMOOTHED_VALUE < DV%SETPOINT) .AND. (DV%CURRENT_STATE .NEQV. DV%INITIAL_STATE)) THEN
            DV%CURRENT_STATE = DV%INITIAL_STATE
         ENDIF
      ELSE DEVICE_DIRECTION2
         IF ((DV%SMOOTHED_VALUE < DV%SETPOINT) .AND. (DV%CURRENT_STATE .EQV.  DV%INITIAL_STATE)) THEN
            DV%CURRENT_STATE = .NOT.DV%INITIAL_STATE
         ELSEIF ((DV%SMOOTHED_VALUE > DV%SETPOINT) .AND. (DV%CURRENT_STATE .NEQV. DV%INITIAL_STATE)) THEN
            DV%CURRENT_STATE = DV%INITIAL_STATE
         ENDIF
      ENDIF DEVICE_DIRECTION2
   ENDIF LATCHIF

   ! If a DEViCe changes state, save the Smokeview file strings and time of state change

   IF (DV%CURRENT_STATE.NEQV.DV%PRIOR_STATE) DV%T_CHANGE = T

   IF (PROCESS(DV%MESH)==MY_RANK .AND. &
       ((DV%CURRENT_STATE.NEQV.DV%PRIOR_STATE) .OR. (ABS(T-T_BEGIN)<SPACING(T).AND..NOT.DV%CURRENT_STATE))) THEN
      M=>MESHES(DV%MESH)
      IF (M%N_STRINGS+2>M%N_STRINGS_MAX) CALL RE_ALLOCATE_STRINGS(DV%MESH)
      I_STATE=0
      IF (DV%CURRENT_STATE) I_STATE=1
      M%N_STRINGS = M%N_STRINGS + 1
      WRITE(M%STRING(M%N_STRINGS),'(A,5X,A,1X)') 'DEVICE_ACT',TRIM(DV%ID)
      M%N_STRINGS = M%N_STRINGS + 1
      WRITE(M%STRING(M%N_STRINGS),'(I6,F10.2,I6)') N,T_BEGIN+(T-T_BEGIN)*TIME_SHRINK_FACTOR,I_STATE
   ENDIF

ENDDO DEVICE_LOOP

T_USED(7) = T_USED(7) + CURRENT_TIME() - TNOW
END SUBROUTINE EXCHANGE_GLOBAL_OUTPUTS


!> \brief Dump HRR data to CHID_hrr.csv, MASS data to CHID_mass.csv, DEVICE data to _devc.csv

SUBROUTINE DUMP_GLOBAL_OUTPUTS

REAL(EB) :: TNOW
TYPE(DEVICE_TYPE), POINTER :: DV

TNOW = CURRENT_TIME()

! Dump out HRR info into CHID_hrr.csv

IF (T>=HRR_CLOCK(HRR_COUNTER(1))) THEN
   IF (MY_RANK==0) THEN
      CALL MPI_REDUCE(MPI_IN_PLACE,Q_DOT_SUM,N_Q_DOT,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
      CALL MPI_REDUCE(MPI_IN_PLACE,M_DOT_SUM,N_TRACKED_SPECIES,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
      CALL DUMP_HRR(T,DT)
   ELSE
      CALL MPI_REDUCE(Q_DOT_SUM,Q_DOT_SUM,N_Q_DOT,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
      CALL MPI_REDUCE(M_DOT_SUM,M_DOT_SUM,N_TRACKED_SPECIES,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
   ENDIF
   HRR_COUNTER(1) = HRR_COUNTER(1) + 1
   Q_DOT_SUM = 0._EB
   M_DOT_SUM = 0._EB
   T_LAST_DUMP_HRR = T
ENDIF

! Dump out HRR info into CHID_hrr.csv

IF (MY_RANK==0 .AND. HVAC_SOLVE .AND. (N_DUCT_QUANTITY > 0 .OR. N_NODE_QUANTITY > 0)) THEN
   IF (T>=HVAC_CLOCK(HVAC_COUNTER(1))) THEN
      CALL DUMP_HVAC(T)
      HVAC_COUNTER(1) = HVAC_COUNTER(1) + 1
   ENDIF
ENDIF

! Dump unstructured geometry and boundary element info

IF (N_FACE>0 .AND. T>=GEOM_CLOCK(GEOM_COUNTER(1))) THEN
   IF (MY_RANK==0) THEN
      CALL DUMP_GEOM(T,DO_CFACES=.FALSE.)
   ENDIF
   IF (ABS(T-T_BEGIN)<TWO_EPSILON_EB) CALL DUMP_GEOM(T,DO_CFACES=.TRUE.)
   GEOM_COUNTER(1) = GEOM_COUNTER(1) + 1
ENDIF

! Dump out mass info into CHID_mass.csv

IF (T>=MASS_CLOCK(MASS_COUNTER(1))) THEN
   IF (MY_RANK==0) THEN
      CALL MPI_REDUCE(MPI_IN_PLACE,MASS_DT,1+N_SPECIES+N_TRACKED_SPECIES,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
      CALL DUMP_MASS(T,DT)
   ELSE
      CALL MPI_REDUCE(MASS_DT,MASS_DT,1+N_SPECIES+N_TRACKED_SPECIES,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
   ENDIF
   MASS_COUNTER(1) = MASS_COUNTER(1) + 1
   MASS_DT   = 0._EB
   T_LAST_DUMP_MASS = T
ENDIF

! Dump device info into CHID_devc.csv

IF (T>=DEVC_CLOCK(DEVC_COUNTER(1)) .AND. N_DEVC>0) THEN

   ! Exchange histogram info

   DO N=1,N_DEVC
      DV => DEVICE(N)
      IF (.NOT.PROPERTY(DV%PROP_INDEX)%HISTOGRAM) CYCLE
      IF (PROCESS(DV%MESH)==MY_RANK .AND. MY_RANK>0) CALL MPI_SEND(DV%HISTOGRAM_COUNTS(1), &
          PROPERTY(DV%PROP_INDEX)%HISTOGRAM_NBINS,MPI_DOUBLE_PRECISION,0,DV%MESH,MPI_COMM_WORLD,IERR)
      IF (PROCESS(DV%MESH)>0 .AND. MY_RANK==0) CALL MPI_RECV(DV%HISTOGRAM_COUNTS(1), &
          PROPERTY(DV%PROP_INDEX)%HISTOGRAM_NBINS,MPI_DOUBLE_PRECISION,PROCESS(DV%MESH),DV%MESH,MPI_COMM_WORLD,STATUS,IERR)
   ENDDO

   ! Dump the device output to file

   IF (MINVAL(DEVICE(1:N_DEVC)%TIME_INTERVAL)>0._EB) THEN
      IF (MY_RANK==0) CALL DUMP_DEVICES(T)
      DEVC_COUNTER(1) = DEVC_COUNTER(1) + 1
      DEVICE_LOOP: DO N=1,N_DEVC
         DV => DEVICE(N)
         IF (T>DV%STATISTICS_END) CYCLE
         IF (DV%NO_UPDATE_DEVC_INDEX>0) THEN
            IF (DEVICE(DV%NO_UPDATE_DEVC_INDEX)%CURRENT_STATE) CYCLE DEVICE_LOOP
         ELSEIF (DV%NO_UPDATE_CTRL_INDEX>0) THEN
            IF (CONTROL(DV%NO_UPDATE_CTRL_INDEX)%CURRENT_STATE) CYCLE DEVICE_LOOP
         ENDIF
         IF (DV%TEMPORAL_STATISTIC=='TIME AVERAGE') THEN
            DV%VALUE = 0._EB
            DV%TIME_INTERVAL = 0._EB
         ENDIF
      ENDDO DEVICE_LOOP
   ENDIF

ENDIF

! DUMP CVODE substeps

IF (WRITE_CVODE_SUBSTEPS) CALL DUMP_CVODE_SUBSTEPS()

! Dump CONTROL info. No gathering required as CONTROL is updated on all meshes

IF (T>=CTRL_CLOCK(CTRL_COUNTER(1)) .AND. N_CTRL>0) THEN
   IF (MY_RANK==0) CALL DUMP_CONTROLS(T)
   CTRL_COUNTER(1) = CTRL_COUNTER(1) + 1
ENDIF

! Dump CPU time

IF (T>=CPU_CLOCK(CPU_COUNTER(1))) THEN
   CALL DUMP_TIMERS
   CPU_COUNTER(1) = CPU_COUNTER(1) + 1
ENDIF

! Update event log

IF (MY_RANK==0 .AND. WRITE_DEVC_CTRL) THEN
   IF (N_DEVC > 0) THEN
      IF (MINVAL(ABS(DEVICE%T_CHANGE-T)) < TWO_EPSILON_EB) CALL WRITE_DEVC_CTRL_LOG('DEVC',T)
   ENDIF
   IF (N_CTRL > 0) THEN
      IF (MINVAL(ABS(CONTROL%T_CHANGE-T)) < TWO_EPSILON_EB) CALL WRITE_DEVC_CTRL_LOG('CTRL',T)
   ENDIF
ENDIF

T_USED(7) = T_USED(7) + CURRENT_TIME() - TNOW
END SUBROUTINE DUMP_GLOBAL_OUTPUTS


!> \brief Exchange information needed for performing the HVAC computation

SUBROUTINE EXCHANGE_HVAC_BC

USE HVAC_ROUTINES, ONLY: NODE_PROPERTIES,NODE_ZONE
REAL(EB) :: TNOW

TNOW = CURRENT_TIME()

! Sum up contributions over MPI processes

CALL MPI_ALLREDUCE(MPI_IN_PLACE,NODE_PROPERTIES,N_DUCTNODES*(8+N_TRACKED_SPECIES),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
CALL MPI_ALLREDUCE(MPI_IN_PLACE,NODE_ZONE,N_DUCTNODES,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERR)

T_USED(11)=T_USED(11) + CURRENT_TIME() - TNOW
END SUBROUTINE EXCHANGE_HVAC_BC


!> \brief Exchange information needed for performing the HVAC computation

SUBROUTINE EXCHANGE_HVAC_SOLUTION

USE HVAC_ROUTINES, ONLY: NODE_AREA_EX,NODE_TMP_EX,NODE_ZZ_EX,DUCT_MF
REAL(EB) :: TNOW
INTEGER :: NN

TNOW = CURRENT_TIME()

REAL_BUFFER_DUCT(            1:  N_DUCTNODES) = NODE_AREA_EX(1:N_DUCTNODES)
REAL_BUFFER_DUCT(N_DUCTNODES+1:2*N_DUCTNODES) = NODE_TMP_EX(1:N_DUCTNODES)
NN = 1
DO N=1,N_TRACKED_SPECIES
   NN = NN + 1
   REAL_BUFFER_DUCT(NN*N_DUCTNODES+1:(NN+1)*N_DUCTNODES) = NODE_ZZ_EX(1:N_DUCTNODES,N)
ENDDO
REAL_BUFFER_DUCT((2+N_TRACKED_SPECIES)*N_DUCTNODES+1:(2+N_TRACKED_SPECIES)*N_DUCTNODES+N_DUCTS) = DUCT_MF(1:N_DUCTS)

CALL MPI_BCAST(REAL_BUFFER_DUCT(1),(2+N_TRACKED_SPECIES)*N_DUCTNODES+N_DUCTS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

NODE_AREA_EX(1:N_DUCTNODES) = REAL_BUFFER_DUCT(            1:  N_DUCTNODES)
NODE_TMP_EX(1:N_DUCTNODES)  = REAL_BUFFER_DUCT(N_DUCTNODES+1:2*N_DUCTNODES)
NN = 1
DO N=1,N_TRACKED_SPECIES
   NN = NN + 1
   NODE_ZZ_EX(1:N_DUCTNODES,N) = REAL_BUFFER_DUCT(NN*N_DUCTNODES+1:(NN+1)*N_DUCTNODES)
ENDDO
DUCT_MF(1:N_DUCTS) = REAL_BUFFER_DUCT((2+N_TRACKED_SPECIES)*N_DUCTNODES+1:(2+N_TRACKED_SPECIES)*N_DUCTNODES+N_DUCTS)

T_USED(11)=T_USED(11) + CURRENT_TIME() - TNOW
END SUBROUTINE EXCHANGE_HVAC_SOLUTION


!> \brief Exchange information for externally controlled variables

SUBROUTINE EXCHANGE_EXTERNAL

REAL(EB) :: TNOW

TNOW = CURRENT_TIME()

IF (N_CTRL > 0) THEN
   LOGICAL_BUFFER_EXTERNAL = EXTERNAL_CTRL
   CALL MPI_BCAST(LOGICAL_BUFFER_EXTERNAL(1),N_CTRL,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
   EXTERNAL_CTRL = LOGICAL_BUFFER_EXTERNAL
ENDIF
IF (N_RAMP > 0) THEN
   REAL_BUFFER_EXTERNAL = EXTERNAL_RAMP
   CALL MPI_BCAST(REAL_BUFFER_EXTERNAL(1),N_RAMP,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
   EXTERNAL_RAMP = REAL_BUFFER_EXTERNAL
ENDIF

END SUBROUTINE EXCHANGE_EXTERNAL


!> \brief Check to see if any FREEZE_VELOCITY=T and any SOLID_PHASE_ONLY=T

SUBROUTINE CHECK_FREEZE_VELOCITY_STATUS

CALL MPI_ALLREDUCE(MPI_IN_PLACE,FREEZE_VELOCITY ,INTEGER_ONE,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,IERR)
CALL MPI_ALLREDUCE(MPI_IN_PLACE,SOLID_PHASE_ONLY,INTEGER_ONE,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,IERR)
IF (FREEZE_VELOCITY) CHECK_FREEZE_VELOCITY = .FALSE.

END SUBROUTINE CHECK_FREEZE_VELOCITY_STATUS


!> \brief Sum up local grid-snapped or "FDS" VENT areas so that each process knows the total vent area spanning multiple meshes

SUBROUTINE EXCHANGE_VENT_AREA
INTEGER :: NM,NV

CALL MPI_ALLREDUCE(MPI_IN_PLACE,VENT_TOTAL_AREA,N_VENT_TOTAL,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
   
DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   M => MESHES(NM)
   DO NV=1,M%N_VENT
      M%VENTS(NV)%TOTAL_FDS_AREA = VENT_TOTAL_AREA(M%VENTS(NV)%TOTAL_INDEX)
   ENDDO
ENDDO

END SUBROUTINE EXCHANGE_VENT_AREA


!> \brief Gather revision dates, etc, from subroutines
!> \param REVISION String containing the revision number
!> \param REVISION_DATE String containing the date of the last code revision
!> \param COMPILE_DATE String containing the date of the last code compilation
!> \details Unlike svn, the revisioning system git does not perform keyword substitution.
!> To perform this function,  a script named expand_file is called before FDS is
!> built that expands the following keywords ($Revision, $RevisionDate and
!> $CompileDate) with their proper values. Another script named contract_file is
!> called after FDS is built to return these keywords back to their original
!> values (so the revisioning system will not think this file has changed).

SUBROUTINE GET_INFO (REVISION,REVISION_DATE,COMPILE_DATE)

CHARACTER(LEN=255), INTENT(OUT) :: REVISION, REVISION_DATE, COMPILE_DATE
CHARACTER(255), PARAMETER :: GREVISION='$Revision$'
CHARACTER(255), PARAMETER :: GREVISION_DATE='$RevisionDate: unknown $'
CHARACTER(255), PARAMETER :: GCOMPILE_DATE='$CompileDate: unknown $'

WRITE(REVISION,'(A)')      GREVISION(INDEX(GREVISION,':')+2:LEN_TRIM(GREVISION)-2)
WRITE(REVISION_DATE,'(A)') GREVISION_DATE(INDEX(GREVISION_DATE,':')+2:LEN_TRIM(GREVISION_DATE)-2)
WRITE(COMPILE_DATE,'(A)')  GCOMPILE_DATE(INDEX(GCOMPILE_DATE,':')+2:LEN_TRIM(GCOMPILE_DATE)-2)

END SUBROUTINE GET_INFO

END PROGRAM FDS
