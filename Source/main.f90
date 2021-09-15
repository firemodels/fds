!> \brief FDS is a computational fluid dynamics (CFD) code designed to model
!> fire and other thermal phenomena.

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
USE COMP_FUNCTIONS, ONLY : CURRENT_TIME
USE DEVICE_VARIABLES
USE WALL_ROUTINES
USE FIRE
USE CONTROL_FUNCTIONS
USE EVAC
USE TURBULENCE, ONLY: NS_ANALYTICAL_SOLUTION,INIT_TURB_ARRAYS,COMPRESSION_WAVE,&
                      TWOD_VORTEX_CERFACS,TWOD_VORTEX_UMD,TWOD_SOBOROT_UMD, &
                      SYNTHETIC_TURBULENCE,SYNTHETIC_EDDY_SETUP,SANDIA_DAT
USE MANUFACTURED_SOLUTIONS, ONLY: SHUNN_MMS_3,SAAD_MMS_1
USE CC_SCALARS_IBM,   ONLY: CCIBM_SET_DATA,CCIBM_END_STEP,CCREGION_DENSITY,CCIBM_TARGET_VELOCITY, &
                            CHECK_SPEC_TRANSPORT_CONSERVE,MASS_CONSERVE_INIT,CCIBM_RHO0W_INTERP,  &
                            CCCOMPUTE_RADIATION,CCIBM_NO_FLUX,CCIBM_COMPUTE_VELOCITY_ERROR,       &
                            FINISH_CCIBM,INIT_CUTCELL_DATA,LINEARFIELDS_INTERP_TEST,              &
                            MESH_CC_EXCHANGE,ROTATED_CUBE_ANN_SOLN
USE OPENMP
USE MPI_F08
USE SCRC, ONLY: SCARC_SETUP, SCARC_SOLVER, SCARC_NO_FLUX, SCARC_FINALIZE_INSEPARABLE_POISSON
USE SOOT_ROUTINES, ONLY: CALC_AGGLOMERATION
USE GLOBMAT_SOLVER, ONLY : GLMAT_SOLVER_SETUP_H, GLMAT_SOLVER_H, COPY_H_OMESH_TO_MESH, &
                           FINISH_GLMAT_SOLVER_H,PRESSURE_SOLVER_CHECK_RESIDUALS_U

IMPLICIT NONE (TYPE,EXTERNAL)

! Miscellaneous declarations

LOGICAL  :: EX=.FALSE.,DIAGNOSTICS,EXCHANGE_EVACUATION=.FALSE.,CTRL_STOP_STATUS,CHECK_FREEZE_VELOCITY=.TRUE.
LOGICAL  :: EVACUATION=.FALSE.
INTEGER  :: LO10,NM,IZERO,ANG_INC_COUNTER
REAL(EB) :: T,DT,DT_EVAC,TNOW
REAL :: CPUTIME
REAL(EB), ALLOCATABLE, DIMENSION(:) ::  TC_GLB,TC_LOC,DT_NEW,TI_LOC,TI_GLB, &
                                        DSUM_ALL,PSUM_ALL,USUM_ALL,DSUM_ALL_LOCAL,PSUM_ALL_LOCAL,USUM_ALL_LOCAL
REAL(EB), ALLOCATABLE, DIMENSION(:,:) ::  TC2_GLB,TC2_LOC
LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: CONNECTED_ZONES_GLOBAL,CONNECTED_ZONES_LOCAL
LOGICAL, ALLOCATABLE, DIMENSION(:) ::  STATE_GLB,STATE_LOC
INTEGER :: IWW,IW,ITER
TYPE (MESH_TYPE), POINTER :: M,M4
TYPE (OMESH_TYPE), POINTER :: M2,M3

! MPI stuff

TYPE (MPI_STATUS) :: STATUS
INTEGER :: N,I,IERR=0
INTEGER :: PNAMELEN=0,TAG_EVAC
INTEGER :: PROVIDED
INTEGER, PARAMETER :: REQUIRED=MPI_THREAD_FUNNELED
TYPE (MPI_REQUEST), ALLOCATABLE, DIMENSION(:) :: REQ,REQ1,REQ2,REQ3,REQ4,REQ5,REQ7,REQ6,REQ14,REQ15
INTEGER, ALLOCATABLE, DIMENSION(:) :: COUNTS,DISPLS,COUNTS_MASS,DISPLS_MASS,COUNTS_HVAC,DISPLS_HVAC,&
                                      COUNTS_QM,DISPLS_QM,COUNTS_10,DISPLS_10,COUNTS_20,DISPLS_20
INTEGER :: N_REQ,N_REQ1=0,N_REQ2=0,N_REQ3=0,N_REQ4=0,N_REQ5=0,N_REQ7=0,N_REQ6=0,N_REQ14=0,N_REQ15=0
CHARACTER(MPI_MAX_PROCESSOR_NAME) :: PNAME
REAL(EB), ALLOCATABLE, DIMENSION(:)       :: REAL_BUFFER_DUCT
REAL(EB), ALLOCATABLE, DIMENSION(:,:)     :: REAL_BUFFER_10,REAL_BUFFER_MASS,REAL_BUFFER_HVAC,REAL_BUFFER_QM,REAL_BUFFER_20

! Initialize OpenMP

CALL OPENMP_INIT

! output version info if fds is invoked without any arguments
! (this must be done before MPI is initialized)

CALL VERSION_INFO

! Initialize MPI (First executable lines of code)

CALL MPI_INIT_THREAD(REQUIRED,PROVIDED,IERR)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, MY_RANK, IERR)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, N_MPI_PROCESSES, IERR)
CALL MPI_GET_PROCESSOR_NAME(PNAME, PNAMELEN, IERR)

IF (MY_RANK==0) WRITE(LU_ERR,'(/A/)') ' Starting FDS ...'
CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

WRITE(LU_ERR,'(A,I6,A,A)') ' MPI Process ',MY_RANK,' started on ',PNAME(1:PNAMELEN)
CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

! Check that MPI processes and OpenMP threads are working properly

CALL CHECK_MPI

! Start wall clock timing

WALL_CLOCK_START = CURRENT_TIME()
CALL CPU_TIME(CPUTIME)
CPU_TIME_START = CPUTIME
ALLOCATE(T_USED(N_TIMERS)) ; T_USED = 0._EB ; T_USED(1) = CURRENT_TIME()

! Assign a compilation date (All Nodes)

CALL GET_INFO (REVISION,REVISION_DATE,COMPILE_DATE)

! Read input from CHID.fds file and stop the code if any errors are found

CALL READ_DATA(DT)

CALL STOP_CHECK(1)

IF (EVACUATION_WRITE_FED) EVAC_PROCESS=MAX(0,N_MPI_PROCESSES-1)
IF (ANY(EVACUATION_ONLY)) THEN
   EVACUATION=.TRUE.
   RADIATION=.FALSE.
   HVAC_SOLVE=.FALSE.
ENDIF

! If SOLID_HT3D=T in any mesh, then set SOLID_HT3D=T in all meshes.

CALL MPI_ALLREDUCE(MPI_IN_PLACE,SOLID_HT3D,INTEGER_ONE,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,IERR)

! Setup number of OPENMP threads

CALL OPENMP_SET_THREADS

! Print OPENMP thread status

CALL OPENMP_PRINT_STATUS

! Set up send and receive buffer counts and displacements

CALL MPI_INITIALIZATION_CHORES(1)

! Open and write to Smokeview and status file (Master Node Only)

CALL ASSIGN_FILE_NAMES

DO N=0,N_MPI_PROCESSES-1
   IF (MY_RANK==N) CALL WRITE_SMOKEVIEW_FILE
   IF (N==N_MPI_PROCESSES-1) EXIT
   IF (SHARED_FILE_SYSTEM) CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
ENDDO

! Shut down the run if it is only for checking the set up

IF (SETUP_ONLY .AND. .NOT.CHECK_MESH_ALIGNMENT) STOP_STATUS = SETUP_ONLY_STOP

! Check for errors and shutdown if found

CALL STOP_CHECK(1)

! MPI process 0 reopens the Smokeview file for additional output

IF (MY_RANK==0) THEN
   OPEN(LU_SMV,FILE=FN_SMV,FORM='FORMATTED', STATUS='OLD',POSITION='APPEND')
   CALL WRITE_STATUS_FILES
ENDIF

! Start the clock

T = T_BEGIN

! Allocate various utility arrays

CALL MPI_INITIALIZATION_CHORES(2)

! Initialize global parameters

CALL INITIALIZE_GLOBAL_VARIABLES

! Initialize radiation

IF (RADIATION) CALL INIT_RADIATION

! Set up persistent MPI communicators for neighborhoods

CALL MPI_INITIALIZATION_CHORES(3)

CALL EXCHANGE_GEOMETRY_INFO

! Allocate and initialize mesh-specific variables, and check to see if the code should stop

DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   CALL INITIALIZE_MESH_VARIABLES_1(DT,NM)
ENDDO

! Stop all the processes if this is just a set-up run

IF (CHECK_MESH_ALIGNMENT) THEN
   IF (MY_RANK==0) CALL INITIALIZE_DIAGNOSTIC_FILE(DT)
   STOP_STATUS = SETUP_ONLY_STOP
   IF (MY_RANK==0) WRITE(LU_ERR,'(A)') ' Checking mesh alignment. This could take a few tens of seconds...'
ENDIF
CALL STOP_CHECK(1)

! Allocate and initialize OMESH arrays to hold "other mesh" data for a given mesh

DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   CALL INITIALIZE_MESH_EXCHANGE_1(NM)
ENDDO

! Allocate "request" arrays to keep track of MPI communications

CALL MPI_INITIALIZATION_CHORES(4)

! Exchange information related to size of OMESH arrays

CALL MPI_INITIALIZATION_CHORES(5)

! Initial complex geometry CC setup

IF (CC_IBM) THEN
   CALL CCIBM_SET_DATA(FIRST_CALL=.TRUE.) ! Define Cartesian cell types (used to define pressure zones), cut-cells, cfaces.
   CALL STOP_CHECK(1)
ENDIF

! Initialize PRESSURE_ZONEs

CALL INITIALIZE_PRESSURE_ZONES

CALL STOP_CHECK(1)

! Allocate some arrays that are dimensioned with N_ZONE

CALL MPI_INITIALIZATION_CHORES(6)

! Allocate and initialize OMESH arrays to hold "other mesh" data for a given mesh

DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   CALL INITIALIZE_MESH_EXCHANGE_2(NM)
ENDDO

IF (MY_RANK==0 .AND. VERBOSE) WRITE(LU_ERR,'(A)') ' Completed INITIALIZE_MESH_EXCHANGE_2'

! Initialize persistent MPI sends and receives and allocate buffer arrays.

CALL POST_RECEIVES(0)
CALL MESH_EXCHANGE(0)

CALL EXCHANGE_GEOMETRY_INFO

! Finish initializing mesh variables

DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   CALL INITIALIZE_MESH_VARIABLES_2(NM)
ENDDO
IF (MY_RANK==0 .AND. VERBOSE) WRITE(LU_ERR,'(A)') ' Completed INITIALIZE_MESH_VARIABLES_2'

! Create arrays and communicators to exchange back wall information across mesh boundaries

CALL INITIALIZE_BACK_WALL_EXCHANGE

CALL STOP_CHECK(1)

! Initialize ScaRC solver

IF (PRES_METHOD == 'SCARC' .OR. PRES_METHOD == 'USCARC') THEN
   CALL SCARC_SETUP
   CALL STOP_CHECK(1)
ENDIF

! Initialize turb arrays

DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   IF (TGA_SURF_INDEX>0) CYCLE
   CALL INIT_TURB_ARRAYS(NM)
ENDDO

! Final complex geometry CC setup

IF (CC_IBM) THEN
   CALL CCIBM_SET_DATA(FIRST_CALL=.FALSE.) ! Interpolation Stencils, Scalar transport MATVEC data, cface RDNs.
   CALL STOP_CHECK(1)
ENDIF

! Initialize the flow field with random noise to eliminate false symmetries

DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   IF (TGA_SURF_INDEX>0 .OR. EVACUATION) CYCLE
   IF (NOISE) CALL INITIAL_NOISE(NM)
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
   IF (UVW_RESTART)      CALL UVW_INIT(NM,CSVFINFO(NM)%UVWFILE)
ENDDO

IF (CC_IBM) THEN
   CALL INIT_CUTCELL_DATA(T_BEGIN,DT)  ! Init centroid data (i.e. rho,zz) on cut-cells, cut-faces and CFACEs.
   IF (PERIODIC_TEST==101) CALL LINEARFIELDS_INTERP_TEST
ENDIF

DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   IF (TGA_SURF_INDEX>0 .OR. EVACUATION) CYCLE
   CALL COMPUTE_VISCOSITY(T_BEGIN,NM,APPLY_TO_ESTIMATED_VARIABLES=.FALSE.) ! needed here for KRES prior to mesh exchange
ENDDO

! Exchange information at mesh boundaries related to the various initialization routines just completed

CALL MESH_EXCHANGE(1)
CALL MESH_EXCHANGE(4)
CALL POST_RECEIVES(6)
CALL MESH_EXCHANGE(6)

! Ensure normal components of velocity match at mesh boundaries and do velocity BCs just in case the flow is not initialized to zero

PREDICTOR = .FALSE.
CORRECTOR = .TRUE.

DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   IF (TGA_SURF_INDEX>0 .OR. EVACUATION) CYCLE
   CALL MATCH_VELOCITY(NM)
   CALL COMPUTE_VISCOSITY(T_BEGIN,NM,APPLY_TO_ESTIMATED_VARIABLES=.FALSE.) ! call again after mesh exchange
   IF (SYNTHETIC_EDDY_METHOD) CALL SYNTHETIC_EDDY_SETUP(NM)
   CALL VELOCITY_BC(T_BEGIN,NM,APPLY_TO_ESTIMATED_VARIABLES=.FALSE.)
   CALL VISCOSITY_BC(NM,APPLY_TO_ESTIMATED_VARIABLES=.FALSE.)
ENDDO

! Iterate surface BCs and radiation in case temperatures are not initialized to ambient

DO ITER=1,INITIAL_RADIATION_ITERATIONS
   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL WALL_BC(T_BEGIN,DT,NM)
      IF (RADIATION) THEN
         CALL COMPUTE_RADIATION(T_BEGIN,NM,ITER)
         IF (CC_IBM) CALL CCCOMPUTE_RADIATION(T_BEGIN,NM,ITER)
      ENDIF
   ENDDO
   DO ANG_INC_COUNTER=1,ANGLE_INCREMENT
      CALL MESH_EXCHANGE(2) ! Exchange radiation intensity at interpolated boundaries
   ENDDO
ENDDO
IF (MY_RANK==0 .AND. VERBOSE) WRITE(LU_ERR,'(A)') ' Initialized Radiation'

IF(CHECK_MASS_CONSERVE) CALL MASS_CONSERVE_INIT
IF (CC_IBM .AND. .NOT.COMPUTE_CUTCELLS_ONLY) CALL CCIBM_RHO0W_INTERP

! Compute divergence just in case the flow field is not initialized to ambient

DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   IF (EVACUATION) CYCLE
   CALL DIVERGENCE_PART_1(T_BEGIN,DT,NM)
ENDDO

! Level Set model for firespread in vegetation

IF (LEVEL_SET_MODE>0 .OR. TERRAIN_CASE) THEN
   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL INITIALIZE_LEVEL_SET_FIRESPREAD_1(NM)
   ENDDO
   CALL MESH_EXCHANGE(14)
   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL INITIALIZE_LEVEL_SET_FIRESPREAD_2(NM)
   ENDDO
ENDIF

! Potentially read data from a previous calculation

IF (RESTART) THEN
   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL READ_RESTART(T,DT,NM)
   ENDDO
   IF (CC_IBM) CALL INIT_CUTCELL_DATA(T,DT)  ! Init centroid data (i.e. rho,zz) on cut-cells and cut-faces.
   CALL STOP_CHECK(1)
ENDIF

! Initialize particle distributions

CALL GENERATE_PARTICLE_DISTRIBUTIONS

! Initialize output files that are mesh-specific

DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   CALL INSERT_ALL_PARTICLES(T,NM)
   IF (TGA_SURF_INDEX<1) CALL INITIALIZE_DEVICES(NM)
   IF (TGA_SURF_INDEX<1) CALL INITIALIZE_PROFILES(NM)
   IF (TGA_SURF_INDEX<1) CALL INITIALIZE_MESH_DUMPS(NM)
ENDDO

IF (MY_RANK==0 .AND. VERBOSE) WRITE(LU_ERR,'(A)') ' Inserted particles'

! Check for any stop flags at this point in the set up.

CALL STOP_CHECK(1)

! Check to see if only a TGA analysis is to be performed

IF (TGA_SURF_INDEX>0) THEN
   IF (MY_RANK==0) CALL TGA_ANALYSIS
   STOP_STATUS = TGA_ANALYSIS_STOP
   CALL STOP_CHECK(1)
ENDIF

! Initialize output files containing global data (Master Node Only)

IF (MY_RANK==0) THEN
   CALL INITIALIZE_GLOBAL_DUMPS(T,DT)
   IF (VERBOSE) WRITE(LU_ERR,'(A)') ' Called INITIALIZE_GLOBAL_DUMPS'
ENDIF

! Initialize GLMat solver for H:

IF (GLMAT_SOLVER) THEN
   CALL GLMAT_SOLVER_SETUP_H(1)
   CALL STOP_CHECK(1)
   CALL MESH_EXCHANGE(3) ! Exchange guard cell info for CCVAR(I,J,K,CGSC) -> HS.
   CALL GLMAT_SOLVER_SETUP_H(2)
   CALL MESH_EXCHANGE(3) ! Exchange guard cell info for CCVAR(I,J,K,UNKH) -> HS.
   CALL GLMAT_SOLVER_SETUP_H(3)
   CALL STOP_CHECK(1)
ENDIF
CALL INIT_EVAC_DUMPS

! Initialize EVACuation routines

IF (EVACUATION.OR.EVACUATION_WRITE_FED) THEN
   CALL INITIALIZE_EVAC
   IF (N_MPI_PROCESSES==1 .OR. (N_MPI_PROCESSES>1 .AND. MY_RANK==EVAC_PROCESS)) CALL INIT_EVAC_GROUPS
ENDIF

! Initialize HVAC variables

IF (HVAC_SOLVE) THEN
   IF (MY_RANK==0) CALL INIT_DUCT_NODE
   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL HVAC_BC_IN(NM)
   ENDDO
   IF (N_MPI_PROCESSES>1) CALL EXCHANGE_HVAC_BC
   IF (MY_RANK==0) THEN
      CALL COLLAPSE_HVAC_BC(T)
      CALL SET_INIT_HVAC
   ENDIF
ENDIF

! Make an initial dump of ambient values

IF (.NOT.RESTART) THEN
   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL UPDATE_GLOBAL_OUTPUTS(T,DT,NM)
      CALL DUMP_MESH_OUTPUTS(T,DT,NM)
   ENDDO
   IF (MY_RANK==0 .AND. VERBOSE) WRITE(LU_ERR,'(A)') ' Called DUMP_MESH_OUTPUTS'
ENDIF

! If there are zones and HVAC pass PSUM

IF (HVAC_SOLVE .AND. N_ZONE>0) CALL EXCHANGE_DIVERGENCE_INFO

! Make an initial dump of global output quantities

IF (.NOT.RESTART) THEN
   CALL EXCHANGE_GLOBAL_OUTPUTS
   CALL UPDATE_CONTROLS(T,0._EB,CTRL_STOP_STATUS,.TRUE.)
   CALL DUMP_GLOBAL_OUTPUTS
ENDIF

! Check for changes in VENT or OBSTruction control and device status at t=T_BEGIN

IF (.NOT.RESTART) CALL CREATE_OR_REMOVE_OBSTRUCTIONS

! Write out character strings to .smv file

CALL WRITE_STRINGS

! Check for evacuation initialization stop

IF (EVACUATION) THEN
   CALL STOP_CHECK(1)
   IF (.NOT.RESTART) ICYC = -EVAC_TIME_ITERATIONS
END IF

! Check for CC_IBM initialization stop

IF (CC_IBM) THEN
   IF (COMPUTE_CUTCELLS_ONLY) STOP_STATUS = SETUP_ONLY_STOP
   CALL STOP_CHECK(1)
ENDIF

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

IF (MY_RANK==0 .AND. VERBOSE) WRITE(LU_ERR,'(A)') ' Start the time-stepping loop'

!***********************************************************************************************************************************
!                                                   MAIN TIMESTEPPING LOOP
!***********************************************************************************************************************************

MAIN_LOOP: DO

   ICYC  = ICYC + 1   ! Time step iterations

   ! Do not print out general diagnostics into .out file every time step

   DIAGNOSTICS = .FALSE.
   EXCHANGE_EVACUATION = .FALSE.

   ! Check for program stops

   INQUIRE(FILE=FN_STOP,EXIST=EX)
   IF (EX .AND. ICYC>=STOP_AT_ITER) THEN
      IF (VERBOSE .AND. STOP_STATUS/=USER_STOP) WRITE(LU_ERR,'(A,I5)') ' STOP file detected, MPI Process =',MY_RANK
      STOP_STATUS = USER_STOP
      DIAGNOSTICS = .TRUE.
   ENDIF

   ! Check to see if the time step can be increased

   !IF (ALL(CHANGE_TIME_STEP_INDEX==1)) DT = MINVAL(DT_NEW,MASK=.NOT.EVACUATION_ONLY)
   IF (ALL(CHANGE_TIME_STEP_INDEX==1)) DT = MINVAL(DT_NEW)

   ! Clip final time step

   IF ((T+DT+DT_END_FILL)>T_END) DT = MAX(T_END-T+TWO_EPSILON_EB,DT_END_MINIMUM)

   ! Determine when to dump out diagnostics to the .out file

   LO10 = INT(LOG10(REAL(MAX(1,ABS(ICYC)),EB)))
   IF (MOD(ICYC,10**LO10)==0 .OR. MOD(ICYC,100)==0 .OR. (T+DT)>=T_END) DIAGNOSTICS = .TRUE.

   ! If evacuation, set up special time iteration parameters

   IF (EVACUATION.OR.EVACUATION_WRITE_FED) CALL EVAC_MAIN_LOOP

   !================================================================================================================================
   !                                           Start of Predictor part of time step
   !================================================================================================================================

   PREDICTOR = .TRUE.
   CORRECTOR = .FALSE.

   ! Begin the finite differencing of the PREDICTOR step

   COMPUTE_FINITE_DIFFERENCES_1: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      IF (EVACUATION_SKIP(NM)) CYCLE COMPUTE_FINITE_DIFFERENCES_1
      CALL INSERT_ALL_PARTICLES(T,NM)
      IF (.NOT.SOLID_PHASE_ONLY .AND. .NOT.FREEZE_VELOCITY) CALL COMPUTE_VISCOSITY(T,NM,APPLY_TO_ESTIMATED_VARIABLES=.FALSE.)
      CALL MASS_FINITE_DIFFERENCES(NM)
   ENDDO COMPUTE_FINITE_DIFFERENCES_1

   ! Estimate quantities at next time step, and decrease/increase time step if necessary based on CFL condition

   FIRST_PASS = .TRUE.

   CHANGE_TIME_STEP_LOOP: DO

      ! Predict species mass fractions at the next time step.

      COMPUTE_DENSITY_LOOP: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
         IF (EVACUATION_SKIP(NM)) CYCLE COMPUTE_DENSITY_LOOP
         CALL DENSITY(T,DT,NM)
         IF (LEVEL_SET_MODE>0) CALL LEVEL_SET_FIRESPREAD(T,DT,NM)
      ENDDO COMPUTE_DENSITY_LOOP

      IF (LEVEL_SET_MODE==2 .AND. CHECK_FREEZE_VELOCITY) CALL CHECK_FREEZE_VELOCITY_STATUS

      IF (CC_IBM) CALL CCREGION_DENSITY(T,DT)

      ! Exchange species mass fractions at interpolated boundaries.

      CALL MESH_EXCHANGE(1)

      ! Exchange level set values, if necessary

      IF (LEVEL_SET_MODE>0) CALL MESH_EXCHANGE(14)

      ! Exchange newly inserted particles, if necessary

      CALL MPI_ALLREDUCE(MPI_IN_PLACE,EXCHANGE_INSERTED_PARTICLES,INTEGER_ONE,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,IERR)
      IF (FIRST_PASS .AND. EXCHANGE_INSERTED_PARTICLES) THEN
         CALL MESH_EXCHANGE(7)
         CALL POST_RECEIVES(11)
         CALL MESH_EXCHANGE(11)
      ENDIF

      ! Calculate convective and diffusive terms of the velocity equation.

      COMPUTE_DIVERGENCE_LOOP: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
         IF (EVACUATION_SKIP(NM)) CYCLE COMPUTE_DIVERGENCE_LOOP
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

      COMPUTE_WALL_BC_LOOP_A: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
         IF (EVACUATION_SKIP(NM)) CYCLE COMPUTE_WALL_BC_LOOP_A
         CALL UPDATE_PARTICLES(T,DT,NM)
         CALL WALL_BC(T,DT,NM)
         CALL PARTICLE_MOMENTUM_TRANSFER(DT,NM)
         CALL DIVERGENCE_PART_1(T,DT,NM)
      ENDDO COMPUTE_WALL_BC_LOOP_A

      ! If there are pressure ZONEs, exchange integrated quantities mesh to mesh for use in the divergence calculation

      IF (N_ZONE>0) CALL EXCHANGE_DIVERGENCE_INFO

      ! Finish the divergence calculation

      FINISH_DIVERGENCE_LOOP: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
         IF (EVACUATION_SKIP(NM)) CYCLE FINISH_DIVERGENCE_LOOP
         CALL DIVERGENCE_PART_2(DT,NM)
      ENDDO FINISH_DIVERGENCE_LOOP

      ! Solve for the pressure at the current time step

      CALL PRESSURE_ITERATION_SCHEME
      CALL EVAC_PRESSURE_ITERATION_SCHEME

      ! Predict the velocity components at the next time step

      CHANGE_TIME_STEP_INDEX = 0
      DT_NEW = DT

      PREDICT_VELOCITY_LOOP: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
         IF (EVACUATION_SKIP(NM)) CYCLE PREDICT_VELOCITY_LOOP
         CALL VELOCITY_PREDICTOR(T+DT,DT,DT_NEW,NM)
      ENDDO PREDICT_VELOCITY_LOOP

      ! Check if there is a numerical instability after updating the velocity field. If there is, exit this loop, finish the time
      ! step, and stop the code.

      CALL STOP_CHECK(0)

      IF (STOP_STATUS==INSTABILITY_STOP) THEN
         DIAGNOSTICS = .TRUE.
         EXIT CHANGE_TIME_STEP_LOOP
      ENDIF

      ! Exchange CHANGE_TIME_STEP_INDEX to determine if the time step needs to be decreased (-1) or increased (1). If any mesh
      ! needs to decrease, or all need to increase, exchange the array of new time step values, DT_NEW.

      IF (N_MPI_PROCESSES>1) THEN
         TNOW = CURRENT_TIME()
         CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,CHANGE_TIME_STEP_INDEX,COUNTS,DISPLS,MPI_INTEGER,MPI_COMM_WORLD,IERR)
         IF (ANY(CHANGE_TIME_STEP_INDEX==-1) .OR. ALL(CHANGE_TIME_STEP_INDEX==1)) &
            CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,DT_NEW,COUNTS,DISPLS,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERR)
         T_USED(11) = T_USED(11) + CURRENT_TIME() - TNOW
      ENDIF

      IF (ANY(CHANGE_TIME_STEP_INDEX==-1)) THEN  ! If the time step was reduced, CYCLE CHANGE_TIME_STEP_LOOP
         !DT = MINVAL(DT_NEW,MASK=.NOT.EVACUATION_ONLY)
         DT = MINVAL(DT_NEW)
         FIRST_PASS = .FALSE.
      ELSE  ! exit the loop and if the time step is to be increased, this will occur at the next time step.
         EXIT CHANGE_TIME_STEP_LOOP
      ENDIF

   ENDDO CHANGE_TIME_STEP_LOOP

   ! If detailed CFL info needed

   IF (CFL_FILE) CALL WRITE_CFL_FILE

   ! Exchange velocity and pressures at interpolated boundaries

   CALL MESH_EXCHANGE(3)

   ! Flux average final velocity to cutfaces. Interpolate H to cut-cells from regular fluid cells.

   IF (CC_IBM) CALL CCIBM_END_STEP(T,DT,DIAGNOSTICS)

   ! Force normal components of velocity to match at interpolated boundaries

   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      IF (EVACUATION_SKIP(NM)) CYCLE
      CALL MATCH_VELOCITY(NM)
   ENDDO

   ! Apply tangential velocity boundary conditions

   VELOCITY_BC_LOOP: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      IF (EVACUATION_SKIP(NM)) CYCLE VELOCITY_BC_LOOP
      IF (SYNTHETIC_EDDY_METHOD) CALL SYNTHETIC_TURBULENCE(DT,T,NM)
      CALL VELOCITY_BC(T,NM,APPLY_TO_ESTIMATED_VARIABLES=.TRUE.)
   ENDDO VELOCITY_BC_LOOP

   ! Advance the time to start the CORRECTOR step

   T = T + DT

   !================================================================================================================================
   !                                           Start of Corrector part of time step
   !================================================================================================================================

   CORRECTOR = .TRUE.
   PREDICTOR = .FALSE.

   ! Check for creation or removal of obsructions

   CALL CREATE_OR_REMOVE_OBSTRUCTIONS

   ! Finite differences for mass and momentum equations for the second half of the time step

   COMPUTE_FINITE_DIFFERENCES_2: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      IF (EVACUATION_SKIP(NM)) CYCLE COMPUTE_FINITE_DIFFERENCES_2
      IF (.NOT.SOLID_PHASE_ONLY .AND. .NOT.FREEZE_VELOCITY) CALL COMPUTE_VISCOSITY(T,NM,APPLY_TO_ESTIMATED_VARIABLES=.TRUE.)
      CALL MASS_FINITE_DIFFERENCES(NM)
      CALL DENSITY(T,DT,NM)
      IF (LEVEL_SET_MODE>0) CALL LEVEL_SET_FIRESPREAD(T,DT,NM)
   ENDDO COMPUTE_FINITE_DIFFERENCES_2

   IF (LEVEL_SET_MODE==2 .AND. CHECK_FREEZE_VELOCITY) CALL CHECK_FREEZE_VELOCITY_STATUS

   IF (CC_IBM) CALL CCREGION_DENSITY(T,DT)
   IF (CHECK_MASS_CONSERVE) CALL CHECK_SPEC_TRANSPORT_CONSERVE(T,DT,DIAGNOSTICS)

   ! Exchange species mass fractions.

   CALL MESH_EXCHANGE(4)
   IF (LEVEL_SET_MODE>0) CALL MESH_EXCHANGE(14)

   ! Apply mass and species boundary conditions, update radiation, particles, and re-compute divergence

   COMPUTE_DIVERGENCE_2: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      IF (EVACUATION_SKIP(NM)) CYCLE COMPUTE_DIVERGENCE_2
      IF (.NOT.SOLID_PHASE_ONLY .AND. .NOT.FREEZE_VELOCITY) THEN
         MESHES(NM)%BAROCLINIC_TERMS_ATTACHED = .FALSE.
         CALL VISCOSITY_BC(NM,APPLY_TO_ESTIMATED_VARIABLES=.TRUE.)
         IF (.NOT.CYLINDRICAL) CALL VELOCITY_FLUX(T,DT,NM,APPLY_TO_ESTIMATED_VARIABLES=.TRUE.)
         IF (     CYLINDRICAL) CALL VELOCITY_FLUX_CYLINDRICAL(T,NM,APPLY_TO_ESTIMATED_VARIABLES=.TRUE.)
      ENDIF
      IF (AGGLOMERATION .AND. ANY(AGGLOMERATION_SMIX_INDEX>0)) CALL CALC_AGGLOMERATION(DT,NM)
      IF (N_REACTIONS>0 .OR. INIT_HRRPUV) CALL COMBUSTION(T,DT,NM)
      IF (ANY(SPECIES_MIXTURE%CONDENSATION_SMIX_INDEX>0)) CALL CONDENSATION_EVAPORATION(DT,NM)
   ENDDO COMPUTE_DIVERGENCE_2

   IF (HVAC_SOLVE) CALL HVAC_CALC(T,DT,.TRUE.)

   COMPUTE_WALL_BC_2A: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      IF (EVACUATION_SKIP(NM)) CYCLE COMPUTE_WALL_BC_2A
      IF (N_REACTIONS>0 .OR. INIT_HRRPUV) CALL COMBUSTION_BC(NM)
      CALL UPDATE_PARTICLES(T,DT,NM)
      CALL WALL_BC(T,DT,NM)
      CALL PARTICLE_MOMENTUM_TRANSFER(DT,NM)
   ENDDO COMPUTE_WALL_BC_2A

   DO ITER=1,RADIATION_ITERATIONS
      DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
         IF (EVACUATION_SKIP(NM)) CYCLE
         CALL COMPUTE_RADIATION(T,NM,ITER)
         IF (CC_IBM) CALL CCCOMPUTE_RADIATION(T,NM,ITER)
      ENDDO
      IF (RADIATION_ITERATIONS>1) THEN  ! Only do an MPI exchange of radiation intensity if multiple iterations are requested.
         DO ANG_INC_COUNTER=1,ANGLE_INCREMENT
            CALL MESH_EXCHANGE(2)
            IF (ICYC>1) EXIT
         ENDDO
      ENDIF
   ENDDO

   ! Start the computation of the divergence term.

   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      IF (EVACUATION_SKIP(NM)) CYCLE
      CALL DIVERGENCE_PART_1(T,DT,NM)
   ENDDO

   ! In most LES fire cases, a correction to the source term in the radiative transport equation is needed.

   IF (RTE_SOURCE_CORRECTION) CALL CALCULATE_RTE_SOURCE_CORRECTION_FACTOR

   ! Exchange global pressure zone information

   IF (N_ZONE>0) CALL EXCHANGE_DIVERGENCE_INFO

   ! Exchange mass loss information for OBSTs abutting interpolated boundaries

   CALL MESH_EXCHANGE(15)  ! Exchange number of OBSTs that have mass info to exchange across mesh boundaries
   CALL POST_RECEIVES(16)
   CALL MESH_EXCHANGE(16)  ! Satellite meshes send mass losses to meshes that actually contain the OBSTstruction
   CALL MESH_EXCHANGE(17)  ! Mesh containing the OBSTruction packs up its new mass to be sent to satellite meshes
   CALL POST_RECEIVES(18)
   CALL MESH_EXCHANGE(18)  ! Satellite meshes receive the new mass of OBSTructions that live in neighboring mesh

   ! Finish computing the divergence

   FINISH_DIVERGENCE_LOOP_2: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      IF (EVACUATION_SKIP(NM)) CYCLE FINISH_DIVERGENCE_LOOP_2
      CALL DIVERGENCE_PART_2(DT,NM)
   ENDDO FINISH_DIVERGENCE_LOOP_2

   ! Solve the pressure equation.

   CALL PRESSURE_ITERATION_SCHEME
   CALL EVAC_PRESSURE_ITERATION_SCHEME

   ! Set up the last big exchange of info.

   CALL EVAC_MESH_EXCHANGE(T_EVAC,T_EVAC_SAVE,I_EVAC,ICYC,EXCHANGE_EVACUATION,0)

   ! Update the  velocity.

   CORRECT_VELOCITY_LOOP: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      IF (EVACUATION_SKIP(NM)) CYCLE CORRECT_VELOCITY_LOOP
      CALL VELOCITY_CORRECTOR(T,DT,NM)
      IF (DIAGNOSTICS .AND. .NOT.EVACUATION) CALL CHECK_DIVERGENCE(NM)
   ENDDO CORRECT_VELOCITY_LOOP

   ! Exchange the number of particles sent from mesh to mesh

   CALL MESH_EXCHANGE(7)
   CALL POST_RECEIVES(11)
   CALL MESH_EXCHANGE(11)

   ! Exchange velocity, pressure, particles at interpolated boundaries

   CALL POST_RECEIVES(6)
   CALL MESH_EXCHANGE(6)

   ! Exchange radiation intensity at interpolated boundaries if only one iteration of the solver is requested.

   IF (RADIATION_ITERATIONS==1) THEN
      DO ANG_INC_COUNTER=1,ANGLE_INCREMENT
         CALL MESH_EXCHANGE(2)
         IF (ICYC>1) EXIT
      ENDDO
   ENDIF

   ! Flux average final velocity to cutfaces. Interpolate H to cut-cells from regular fluid cells.

   IF (CC_IBM) CALL CCIBM_END_STEP(T,DT,DIAGNOSTICS)

   ! Force normal components of velocity to match at interpolated boundaries

   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      IF (EVACUATION_SKIP(NM)) CYCLE
      CALL MATCH_VELOCITY(NM)
   ENDDO

   ! Apply velocity boundary conditions, and update values of HRR, DEVC, etc.

   VELOCITY_BC_LOOP_2: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      IF (EVACUATION_SKIP(NM)) CYCLE VELOCITY_BC_LOOP_2
      CALL VELOCITY_BC(T,NM,APPLY_TO_ESTIMATED_VARIABLES=.FALSE.)
      CALL UPDATE_GLOBAL_OUTPUTS(T,DT,NM)
   ENDDO VELOCITY_BC_LOOP_2

   ! Share device, HRR, mass data among all processes

   CALL EXCHANGE_GLOBAL_OUTPUTS

   ! Check for dumping end of timestep outputs

   CALL UPDATE_CONTROLS(T,DT,CTRL_STOP_STATUS,.FALSE.)
   IF (CTRL_STOP_STATUS) STOP_STATUS = CTRL_STOP

   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      IF (EVACUATION_SKIP(NM)) CYCLE
      CALL DUMP_MESH_OUTPUTS(T,DT,NM)
   ENDDO

   ! Dump outputs such as HRR, DEVC, etc.

   CALL DUMP_GLOBAL_OUTPUTS

   ! Exchange EVAC information among meshes

   CALL EVAC_EXCHANGE

   ! Dump out diagnostics

   IF (DIAGNOSTICS) THEN
      CALL WRITE_STRINGS
      IF (.NOT.SUPPRESS_DIAGNOSTICS .AND. N_MPI_PROCESSES>1) CALL EXCHANGE_DIAGNOSTICS
      IF (MY_RANK==0) CALL WRITE_DIAGNOSTICS(T,DT)
   ENDIF

   ! Flush output file buffers

   IF (T>=FLUSH_CLOCK .AND. FLUSH_FILE_BUFFERS) THEN
      IF (MY_RANK==0) CALL FLUSH_GLOBAL_BUFFERS
      IF (MY_RANK==MAX(0,EVAC_PROCESS)) CALL FLUSH_EVACUATION_BUFFERS
      DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
         IF (EVACUATION) CYCLE
         CALL FLUSH_LOCAL_BUFFERS(NM)
      ENDDO
      FLUSH_CLOCK = FLUSH_CLOCK + DT_FLUSH
   ENDIF

   ! Dump a restart file if necessary

   CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,RADIATION_COMPLETED,COUNTS,DISPLS,MPI_LOGICAL,MPI_COMM_WORLD,IERR)
   IF ( (T>=RESTART_CLOCK .OR. STOP_STATUS==USER_STOP) .AND. (T>=T_END .OR. ALL(RADIATION_COMPLETED)) ) THEN
      DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
         IF (EVACUATION_SKIP(NM)) CYCLE
         CALL DUMP_RESTART(T,DT,NM)
      ENDDO
      RESTART_CLOCK = RESTART_CLOCK + DT_RESTART
   ENDIF

   ! Check for abnormal run stop

   CALL STOP_CHECK(1)  ! The argument 1 means that FDS will end unless there is logic associated with the STOP_STATUS

   ! Stop the run normally

   IF (T>=T_END .AND. ICYC>0) EXIT MAIN_LOOP

ENDDO MAIN_LOOP

!***********************************************************************************************************************************
!                                                     END OF TIME STEPPING LOOP
!***********************************************************************************************************************************

! Deallocate GLMAT_SOLVER_H variables if needed:

IF (PRES_METHOD == 'GLMAT') CALL FINISH_GLMAT_SOLVER_H

! Finish unstructured geometry

IF (CC_IBM) CALL FINISH_CCIBM

! Finish the smv file (fire+evacuation results) write for the evacuation MC mode

IF (EVACUATION_MC_MODE .AND. .NOT.EVACUATION_DRILL) CALL FINISH_EVACUATION_SMV

! Stop the calculation

CALL END_FDS

! This is the end of program. Supporting routines are listed below.

CONTAINS


SUBROUTINE CHECK_MPI

! Check the threading support level

IF (PROVIDED<REQUIRED) THEN
   IF (MY_RANK==0) WRITE(LU_ERR,'(A)') ' WARNING:  This MPI implementation provides insufficient threading support.'
   !$ CALL OMP_SET_NUM_THREADS(1)
ENDIF

END SUBROUTINE CHECK_MPI


SUBROUTINE MPI_INITIALIZATION_CHORES(TASK_NUMBER)

INTEGER, INTENT(IN) :: TASK_NUMBER
TYPE (MPI_REQUEST), ALLOCATABLE, DIMENSION(:) :: REQ0
TYPE (MPI_GROUP) :: GROUP_WORLD,SUBGROUP
INTEGER :: N_REQ0,SNODE,MEMBERS(0:NMESHES-1),NN,NOM,N_COMMUNICATIONS

SELECT CASE(TASK_NUMBER)

   CASE(1)

      ! Set up send and receive buffer counts and displacements

      ALLOCATE(REAL_BUFFER_DUCT((2+N_TRACKED_SPECIES)*N_DUCTNODES+N_DUCTS))
      ALLOCATE(REAL_BUFFER_MASS(0:N_SPECIES+N_TRACKED_SPECIES,NMESHES))
      ALLOCATE(REAL_BUFFER_HVAC((9+N_TRACKED_SPECIES)*N_DUCTNODES,NMESHES))
      ALLOCATE(REAL_BUFFER_QM(N_Q_DOT+N_TRACKED_SPECIES,NMESHES))
      ALLOCATE(REAL_BUFFER_10(10,NMESHES))
      ALLOCATE(REAL_BUFFER_20(20,NMESHES))

      ALLOCATE(COUNTS(0:N_MPI_PROCESSES-1))
      ALLOCATE(COUNTS_HVAC(0:N_MPI_PROCESSES-1))
      ALLOCATE(COUNTS_MASS(0:N_MPI_PROCESSES-1))
      ALLOCATE(COUNTS_QM(0:N_MPI_PROCESSES-1))
      ALLOCATE(COUNTS_10(0:N_MPI_PROCESSES-1))
      ALLOCATE(COUNTS_20(0:N_MPI_PROCESSES-1))
      ALLOCATE(COUNTS_TP(0:N_MPI_PROCESSES-1))

      ALLOCATE(DISPLS(0:N_MPI_PROCESSES-1))
      ALLOCATE(DISPLS_MASS(0:N_MPI_PROCESSES-1))
      ALLOCATE(DISPLS_HVAC(0:N_MPI_PROCESSES-1))
      ALLOCATE(DISPLS_QM(0:N_MPI_PROCESSES-1))
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
      COUNTS_HVAC   = COUNTS*((9+N_TRACKED_SPECIES)*N_DUCTNODES)
      DISPLS_HVAC   = DISPLS*((9+N_TRACKED_SPECIES)*N_DUCTNODES)
      COUNTS_MASS   = COUNTS*(N_SPECIES+N_TRACKED_SPECIES+1)
      DISPLS_MASS   = DISPLS*(N_SPECIES+N_TRACKED_SPECIES+1)
      COUNTS_QM = COUNTS*(N_Q_DOT+N_TRACKED_SPECIES)
      DISPLS_QM = DISPLS*(N_Q_DOT+N_TRACKED_SPECIES)
      COUNTS_10 = COUNTS*10
      DISPLS_10 = DISPLS*10
      COUNTS_20 = COUNTS*20
      DISPLS_20 = DISPLS*20

   CASE(2)

      ! Allocate TIME arrays

      ALLOCATE(DT_NEW(NMESHES),STAT=IZERO) ;  CALL ChkMemErr('MAIN','DT_NEW',IZERO) ; DT_NEW = DT

      ! Set up dummy arrays to hold various arrays that must be exchanged among meshes

      ALLOCATE(TI_LOC(N_DEVC),STAT=IZERO)
      CALL ChkMemErr('MAIN','TI_LOC',IZERO)
      ALLOCATE(TI_GLB(N_DEVC),STAT=IZERO)
      CALL ChkMemErr('MAIN','TI_GLB',IZERO)
      ALLOCATE(STATE_GLB(2*N_DEVC),STAT=IZERO)
      CALL ChkMemErr('MAIN','STATE_GLB',IZERO)
      ALLOCATE(STATE_LOC(2*N_DEVC),STAT=IZERO)
      CALL ChkMemErr('MAIN','STATE_LOC',IZERO)
      ALLOCATE(TC_GLB(3*N_DEVC),STAT=IZERO)
      CALL ChkMemErr('MAIN','TC_GLB',IZERO)
      ALLOCATE(TC_LOC(3*N_DEVC),STAT=IZERO)
      CALL ChkMemErr('MAIN','TC_LOC',IZERO)
      ALLOCATE(TC2_GLB(2,N_DEVC),STAT=IZERO)
      CALL ChkMemErr('MAIN','TC2_GLB',IZERO)
      ALLOCATE(TC2_LOC(2,N_DEVC),STAT=IZERO)
      CALL ChkMemErr('MAIN','TC2_LOC',IZERO)


   CASE(3)

      ! Allocate arrays that hold geometric information of neighboring meshes.

      IF (N_MPI_PROCESSES>1) THEN
         CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,CELL_COUNT,COUNTS,DISPLS,MPI_INTEGER,MPI_COMM_WORLD,IERR)
      ENDIF

      DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
         M => MESHES(NM)
         ALLOCATE(M%WALL_INDEX(0:CELL_COUNT(NM),-3:3),STAT=IZERO) ; M%WALL_INDEX = 0
      ENDDO

      DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
         IF (EVACUATION) CYCLE
         M => MESHES(NM)
         DO N=1,M%N_NEIGHBORING_MESHES
            NOM = M%NEIGHBORING_MESH(N)
            SNODE = PROCESS(NOM)
            IF (SNODE==MY_RANK) CYCLE
            M4 => MESHES(NOM)
            IF (.NOT.ALLOCATED(M4%CELL_INDEX))    ALLOCATE(M4%CELL_INDEX(0:M4%IBP1,0:M4%JBP1,0:M4%KBP1))
            IF (.NOT.ALLOCATED(M4%SOLID))         ALLOCATE(M4%SOLID(0:CELL_COUNT(NOM)))
            IF (.NOT.ALLOCATED(M4%WALL_INDEX))    ALLOCATE(M4%WALL_INDEX(0:CELL_COUNT(NOM),-3:3))
            IF (.NOT.ALLOCATED(M4%OBST_INDEX_C))  ALLOCATE(M4%OBST_INDEX_C(0:CELL_COUNT(NOM)))
            IF (.NOT.ALLOCATED(M4%PRESSURE_ZONE)) THEN
               ALLOCATE(M4%PRESSURE_ZONE(0:M4%IBP1,0:M4%JBP1,0:M4%KBP1))
               M4%PRESSURE_ZONE = -1
            ENDIF
         ENDDO
      ENDDO

      ! Create MPI communicators (sub-groups of MPI_COMM_WORLD) for each mesh
      ! that consists of the MPI processes that control it and its neighbors
      ! MPI_COMM_MESH(NM) is the name of the communicator for mesh NM
      ! MPI_COMM_MESH_ROOT(NM) is the rank within MPI_COMM_MESH(NM) of the process that controls mesh NM

      ALLOCATE(MPI_COMM_MESH(NMESHES))
      ALLOCATE(MPI_COMM_MESH_ROOT(NMESHES))

      CALL MPI_COMM_GROUP(MPI_COMM_WORLD,GROUP_WORLD,IERR)  ! Get the group handle for MPI_COMM_WORLD

      DO NM=1,NMESHES
         IF (EVACUATION) THEN
            MPI_COMM_MESH(NM) = MPI_COMM_NULL
            CYCLE
         END IF
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
            IF (NM==NOM) MPI_COMM_MESH_ROOT(NM) = NN
         ENDDO
         CALL MPI_GROUP_INCL(GROUP_WORLD,NN+1,MEMBERS(0:NN),SUBGROUP,IERR)  ! Create the new sub-group of GROUP_WORLD
         CALL MPI_COMM_CREATE(MPI_COMM_WORLD,SUBGROUP,MPI_COMM_MESH(NM),IERR) ! Create the new communicator
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
      ALLOCATE(REQ2(N_COMMUNICATIONS*4))
      ALLOCATE(REQ3(N_COMMUNICATIONS*4))
      ALLOCATE(REQ4(N_COMMUNICATIONS*4))
      ALLOCATE(REQ5(N_COMMUNICATIONS*4))
      ALLOCATE(REQ6(N_COMMUNICATIONS*4))
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

         ALLOCATE(REQ0(NMESHES**2)) ; N_REQ0 = 0

         DO NM=1,NMESHES
            IF (EVACUATION) CYCLE
            DO NOM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
               IF (NM/=NOM .AND. PROCESS(NM)/=MY_RANK .AND. MESHES(NOM)%CONNECTED_MESH(NM)) THEN
                  M2 => MESHES(NOM)%OMESH(NM)
                  N_REQ0 = N_REQ0 + 1
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

      DO NM=1,NMESHES
         IF (PROCESS(NM)/=MY_RANK) CYCLE
         IF (EVACUATION) CYCLE
         DO NOM=1,NMESHES
            IF (.NOT.MESHES(NM)%CONNECTED_MESH(NOM)) CYCLE
            M3 => MESHES(NM)%OMESH(NOM)
            IF (N_MPI_PROCESSES>1 .AND. NM/=NOM .AND. PROCESS(NOM)/=MY_RANK .AND. MESHES(NM)%CONNECTED_MESH(NOM)) THEN
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

      IF (N_MPI_PROCESSES>1) THEN

         CALL MPI_WAITALL(N_REQ0,REQ0(1:N_REQ0),MPI_STATUSES_IGNORE,IERR)

         DO NM=1,NMESHES
            IF (EVACUATION) CYCLE
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
         IF (EVACUATION) CYCLE
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
         IF (EVACUATION) CYCLE
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
         IF (EVACUATION) CYCLE
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

      IF (N_REQ>0 .AND. N_MPI_PROCESSES>1) CALL MPI_WAITALL(N_REQ,REQ(1:N_REQ),MPI_STATUSES_IGNORE,IERR)


      CASE(6)

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

      ALLOCATE(CONNECTED_ZONES(0:N_ZONE,0:N_ZONE,NMESHES),STAT=IZERO)
      CALL ChkMemErr('INIT','CONNECTED_ZONES',IZERO)
      CONNECTED_ZONES = .FALSE.

      ALLOCATE(DSUM(N_ZONE,NMESHES),STAT=IZERO) ; CALL ChkMemErr('MAIN','DSUM',IZERO) ; DSUM = 0._EB
      ALLOCATE(PSUM(N_ZONE,NMESHES),STAT=IZERO) ; CALL ChkMemErr('MAIN','PSUM',IZERO) ; PSUM = 0._EB
      ALLOCATE(USUM(N_ZONE,NMESHES),STAT=IZERO) ; CALL ChkMemErr('MAIN','USUM',IZERO) ; USUM = 0._EB

END SELECT

IF (MY_RANK==0 .AND. VERBOSE) WRITE(LU_ERR,'(A,I2)') ' Completed Initialization Step ',TASK_NUMBER

END SUBROUTINE MPI_INITIALIZATION_CHORES


SUBROUTINE PRESSURE_ITERATION_SCHEME

! Iterate calls to pressure solver until velocity tolerance is satisfied

INTEGER :: NM_MAX_V,NM_MAX_P
REAL(EB) :: TNOW,VELOCITY_ERROR_MAX_OLD,PRESSURE_ERROR_MAX_OLD

PRESSURE_ITERATIONS = 0

IF (BAROCLINIC) THEN
   ITERATE_BAROCLINIC_TERM = .TRUE.
ELSE
   ITERATE_BAROCLINIC_TERM = .FALSE.
ENDIF

PRESSURE_ITERATION_LOOP: DO

   PRESSURE_ITERATIONS = PRESSURE_ITERATIONS + 1
   TOTAL_PRESSURE_ITERATIONS = TOTAL_PRESSURE_ITERATIONS + 1

   ! The following loops and exchange always get executed the first pass through the PRESSURE_ITERATION_LOOP.
   ! If we need to iterate the baroclinic torque term, the loop is executed each time.

   IF (ITERATE_BAROCLINIC_TERM .OR. PRESSURE_ITERATIONS==1 .OR. CC_IBM) THEN
      DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
         IF (EVACUATION) CYCLE
         IF (BAROCLINIC) CALL BAROCLINIC_CORRECTION(T,NM)
         IF (CC_IBM) THEN
            ! Wall model to define target velocities in gas cut faces.
            IF(PRESSURE_ITERATIONS<=1 .AND. .NOT.CC_STRESS_METHOD) CALL CCIBM_TARGET_VELOCITY(DT,NM)
            CALL CCIBM_NO_FLUX(DT,NM,.TRUE.) ! IBM Force, we do it here to get velocity flux matched.
         ENDIF
      ENDDO
      CALL MESH_EXCHANGE(5)  ! Exchange FVX, FVY, FVZ
      DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
         IF (EVACUATION) CYCLE
         CALL MATCH_VELOCITY_FLUX(NM)
      ENDDO
   ENDIF

   ! Compute the right hand side (RHS) and boundary conditions for the Poission equation for pressure.
   ! The WALL_WORK1 array is computed in COMPUTE_VELOCITY_ERROR, but it should
   ! be zero the first time the pressure solver is called.

   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      IF (EVACUATION) CYCLE
      IF (INSEPARABLE_POISSON.AND.(PRES_METHOD=='SCARC'.OR.PRES_METHOD=='USCARC')) THEN
         CALL SCARC_NO_FLUX(DT,NM)
      ELSE
         CALL NO_FLUX(DT,NM)
      ENDIF
      IF (CC_IBM) CALL CCIBM_NO_FLUX(DT,NM,.FALSE.) ! set WALL_WORK1 to 0 in cells inside geometries.
      IF (PRESSURE_ITERATIONS==1) MESHES(NM)%WALL_WORK1 = 0._EB
      CALL PRESSURE_SOLVER_COMPUTE_RHS(T,DT,NM)
   ENDDO

   ! Solve the Poission equation using either FFT, SCARC, or GLMAT

   SELECT CASE(PRES_METHOD)
      CASE ('FFT')
         DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
            IF (EVACUATION) CYCLE
            CALL PRESSURE_SOLVER_FFT(NM)
         ENDDO
      CASE ('SCARC','USCARC')
         CALL SCARC_SOLVER(T,DT)
         IF (INSEPARABLE_POISSON) THEN
            DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
               CALL NO_FLUX(DT,NM)
            ENDDO
            CALL SCARC_FINALIZE_INSEPARABLE_POISSON
         ENDIF
         CALL STOP_CHECK(1)
      CASE ('GLMAT')
         CALL GLMAT_SOLVER_H
         CALL MESH_EXCHANGE(5)
         CALL COPY_H_OMESH_TO_MESH
   END SELECT

   ! Check the residuals of the Poisson solution

   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      IF (EVACUATION) CYCLE
      IF (.NOT.PRES_ON_WHOLE_DOMAIN .AND. (PRES_METHOD=='GLMAT'.OR.PRES_METHOD=='USCARC')) THEN ! unstructured global solvers
         CALL PRESSURE_SOLVER_CHECK_RESIDUALS_U(NM)
      ELSE
         CALL PRESSURE_SOLVER_CHECK_RESIDUALS(NM)
      ENDIF
   ENDDO


   IF (.NOT.ITERATE_PRESSURE) EXIT PRESSURE_ITERATION_LOOP

   ! Exchange both H or HS and FVX, FVY, FVZ and then estimate values of U, V, W (US, VS, WS) at next time step.

   CALL MESH_EXCHANGE(5)

   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      IF (EVACUATION) CYCLE
      CALL COMPUTE_VELOCITY_ERROR(DT,NM)
      IF (CC_IBM) CALL CCIBM_COMPUTE_VELOCITY_ERROR(DT,NM) ! Inside solids respect to zero velocity.
   ENDDO

   ! Make all MPI processes aware of the maximum velocity error to decide if another pressure iteration is needed.

   IF (N_MPI_PROCESSES>1) THEN
      TNOW = CURRENT_TIME()
      REAL_BUFFER_10(  1,LOWER_MESH_INDEX:UPPER_MESH_INDEX) = VELOCITY_ERROR_MAX(LOWER_MESH_INDEX:UPPER_MESH_INDEX)
      REAL_BUFFER_10(  2,LOWER_MESH_INDEX:UPPER_MESH_INDEX) = PRESSURE_ERROR_MAX(LOWER_MESH_INDEX:UPPER_MESH_INDEX)
      REAL_BUFFER_10(3:5,LOWER_MESH_INDEX:UPPER_MESH_INDEX) = VELOCITY_ERROR_MAX_LOC(1:3,LOWER_MESH_INDEX:UPPER_MESH_INDEX)
      REAL_BUFFER_10(6:8,LOWER_MESH_INDEX:UPPER_MESH_INDEX) = PRESSURE_ERROR_MAX_LOC(1:3,LOWER_MESH_INDEX:UPPER_MESH_INDEX)
      CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,REAL_BUFFER_10,COUNTS_10,DISPLS_10,MPI_DOUBLE_PRECISION,&
                          MPI_COMM_WORLD,IERR)
      VELOCITY_ERROR_MAX(:)         =     REAL_BUFFER_10(1,:)
      PRESSURE_ERROR_MAX(:)         =     REAL_BUFFER_10(2,:)
      VELOCITY_ERROR_MAX_LOC(1:3,:) = INT(REAL_BUFFER_10(3:5,:))
      PRESSURE_ERROR_MAX_LOC(1:3,:) = INT(REAL_BUFFER_10(6:8,:))
      T_USED(11)=T_USED(11) + CURRENT_TIME() - TNOW
   ENDIF

   IF (MY_RANK==0 .AND. VELOCITY_ERROR_FILE .AND. .NOT.EVACUATION) THEN
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

   IF ((MAXVAL(PRESSURE_ERROR_MAX)<PRESSURE_TOLERANCE .AND. &
        MAXVAL(VELOCITY_ERROR_MAX)<VELOCITY_TOLERANCE) .OR. PRESSURE_ITERATIONS>=MAX_PRESSURE_ITERATIONS) &
      EXIT PRESSURE_ITERATION_LOOP

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


SUBROUTINE CALCULATE_RTE_SOURCE_CORRECTION_FACTOR

! This routine computes a running average of the source correction factor for the radiative transport scheme.

REAL(EB), PARAMETER :: WGT=0.5_EB
REAL(EB) :: RAD_Q_SUM_ALL,KFST4_SUM_ALL,TNOW

TNOW = CURRENT_TIME()

! Sum up the components of the corrective factor from all the meshes.

IF (N_MPI_PROCESSES>1) THEN
   CALL MPI_ALLREDUCE(RAD_Q_SUM,RAD_Q_SUM_ALL,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
   CALL MPI_ALLREDUCE(KFST4_SUM,KFST4_SUM_ALL,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
ELSE
   RAD_Q_SUM_ALL = RAD_Q_SUM
   KFST4_SUM_ALL = KFST4_SUM
ENDIF

! Compute the corrective factor for the RTE. Note that the max value of 100 is arbitrary.

IF (KFST4_SUM_ALL>TWO_EPSILON_EB) &
   RTE_SOURCE_CORRECTION_FACTOR = WGT*RTE_SOURCE_CORRECTION_FACTOR + (1._EB-WGT)*MIN(C_MAX,MAX(C_MIN,RAD_Q_SUM_ALL/KFST4_SUM_ALL))

! Reset the components of the corrective factor to zero.

RAD_Q_SUM = 0._EB
KFST4_SUM = 0._EB

T_USED(11)=T_USED(11) + CURRENT_TIME() - TNOW
END SUBROUTINE CALCULATE_RTE_SOURCE_CORRECTION_FACTOR


!> \brief Create or remove obstructions if necessary
!>
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
   OBST_CREATED_OR_REMOVED = .FALSE.
ENDIF

END SUBROUTINE CREATE_OR_REMOVE_OBSTRUCTIONS


SUBROUTINE WRITE_CFL_FILE

! This routine gathers all the CFL values and mesh indices to node 0, which then
! writes out the max value and mesh and indices of the max value.

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
   CALL MPI_GATHERV(REAL_BUFFER_20(1,DISPLS(MY_RANK)+1),COUNTS_20(MY_RANK),MPI_DOUBLE_PRECISION,REAL_BUFFER_20,COUNTS_20,DISPLS_20,&
                    MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
ELSE
   CALL MPI_GATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,REAL_BUFFER_20,COUNTS_20,DISPLS_20,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
ENDIF

IF (MY_RANK==0 .AND. .NOT.EVACUATION) THEN
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


SUBROUTINE STOP_CHECK(END_CODE)

INTEGER, INTENT(IN) :: END_CODE
REAL(EB) :: TNOW

! Make sure that all MPI processes have the same STOP_STATUS

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


SUBROUTINE END_FDS

! End the calculation gracefully, even if there is an error

CHARACTER(255) :: MESSAGE

IF (STOP_STATUS==NO_STOP .OR. STOP_STATUS==USER_STOP) CALL DUMP_TIMERS

IF (VERBOSE) WRITE(LU_ERR,'(A,I6,A)') ' MPI process ',MY_RANK,' has completed'

IF (MY_RANK==0) THEN

   ! Print out device activation times to the .out file

   CALL TIMINGS

   ! Print out stop status to .err and .out files

   SELECT CASE(STOP_STATUS)
      CASE(NO_STOP)
         WRITE(MESSAGE,'(A)') 'STOP: FDS completed successfully'
         IF (STATUS_FILES) CLOSE(LU_NOTREADY,STATUS='DELETE')
      CASE(INSTABILITY_STOP)
         WRITE(MESSAGE,'(A)') 'ERROR: Numerical Instability - FDS stopped'
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
      CASE(MPI_TIMEOUT_STOP)
         WRITE(MESSAGE,'(A)') 'ERROR: An MPI exchange timed out - FDS stopped'
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
DO I=1,N_REQ15 ; CALL MPI_REQUEST_FREE(REQ15(I),IERR) ; ENDDO

! Shutdown MPI

CALL MPI_FINALIZE(IERR)

! Shutdown FDS

STOP

END SUBROUTINE END_FDS


!> \brief Broadcast basic geometry info from each mesh, NM, to its neighbors in MPI_COMM_MESH(NM)

SUBROUTINE EXCHANGE_GEOMETRY_INFO

IF (N_MPI_PROCESSES==1) RETURN

DO NM=1,NMESHES
   IF (MPI_COMM_MESH(NM)==MPI_COMM_NULL) CYCLE
   M => MESHES(NM)
   CALL MPI_BCAST(M%CELL_INDEX,SIZE(M%CELL_INDEX),MPI_INTEGER,MPI_COMM_MESH_ROOT(NM),MPI_COMM_MESH(NM),IERR)
   CALL MPI_BCAST(M%SOLID,SIZE(M%SOLID),MPI_INTEGER,MPI_COMM_MESH_ROOT(NM),MPI_COMM_MESH(NM),IERR)
   CALL MPI_BCAST(M%WALL_INDEX,SIZE(M%WALL_INDEX),MPI_INTEGER,MPI_COMM_MESH_ROOT(NM),MPI_COMM_MESH(NM),IERR)
   CALL MPI_BCAST(M%OBST_INDEX_C,SIZE(M%OBST_INDEX_C),MPI_INTEGER,MPI_COMM_MESH_ROOT(NM),MPI_COMM_MESH(NM),IERR)
ENDDO

END SUBROUTINE EXCHANGE_GEOMETRY_INFO


!> \brief Exchange information mesh to mesh needed for divergence integrals
!> \details Sum DSUM, PSUM and USUM over all meshes controlled by the active process, then reduce over all processes

SUBROUTINE EXCHANGE_DIVERGENCE_INFO

INTEGER :: IPZ,IOPZ,IOPZ2
REAL(EB) :: TNOW

TNOW = CURRENT_TIME()

CONNECTED_ZONES_LOCAL = .FALSE.

DO IPZ=1,N_ZONE
   DSUM_ALL_LOCAL(IPZ) = 0._EB
   PSUM_ALL_LOCAL(IPZ) = 0._EB
   USUM_ALL_LOCAL(IPZ) = 0._EB
   IF(P_ZONE(IPZ)%EVACUATION) CYCLE
   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      IF(EVACUATION) CYCLE
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
   CALL MPI_ALLREDUCE(PSUM_ALL_LOCAL(1),PSUM_ALL(1),N_ZONE,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
   CALL MPI_ALLREDUCE(USUM_ALL_LOCAL(1),USUM_ALL(1),N_ZONE,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
   CALL MPI_ALLREDUCE(CONNECTED_ZONES_LOCAL(0,0),CONNECTED_ZONES_GLOBAL(0,0),(N_ZONE+1)**2,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,IERR)
ELSE
   DSUM_ALL = DSUM_ALL_LOCAL
   PSUM_ALL = PSUM_ALL_LOCAL
   USUM_ALL = USUM_ALL_LOCAL
   CONNECTED_ZONES_GLOBAL = CONNECTED_ZONES_LOCAL
ENDIF

DO IPZ=1,N_ZONE
   IF(P_ZONE(IPZ)%EVACUATION) CYCLE
   DO NM=1,NMESHES
      IF(EVACUATION) CYCLE
      DSUM(IPZ,NM) = DSUM_ALL(IPZ)
      PSUM(IPZ,NM) = PSUM_ALL(IPZ)
      USUM(IPZ,NM) = USUM_ALL(IPZ)
      CONNECTED_ZONES(IPZ,:,NM) = CONNECTED_ZONES_GLOBAL(IPZ,:)
   ENDDO
ENDDO

! Connect zones to others which are not directly connected

DO NM=1,NMESHES
   IF(EVACUATION) CYCLE
   DO IPZ=1,N_ZONE
      !IF(P_ZONE(IPZ)%EVACUATION) CYCLE
      DO IOPZ=1,N_ZONE
         !IF(P_ZONE(IOPZ)%EVACUATION) CYCLE
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

T_USED(11)=T_USED(11) + CURRENT_TIME() - TNOW
END SUBROUTINE EXCHANGE_DIVERGENCE_INFO


!> \brief Create arrays by which info is to be exchanged across meshes

SUBROUTINE INITIALIZE_MESH_EXCHANGE_1(NM)

INTEGER :: IMIN,IMAX,JMIN,JMAX,KMIN,KMAX,NOM,IOR,IW,N,N_STORAGE_SLOTS,IIO,JJO,KKO,NIC_R,II,JJ,KK,ILC
INTEGER, INTENT(IN) :: NM
TYPE (MESH_TYPE), POINTER :: M2,M
TYPE (OMESH_TYPE), POINTER :: OM
TYPE (EXTERNAL_WALL_TYPE), POINTER :: EWC
TYPE (WALL_TYPE), POINTER :: WC
TYPE (LAGRANGIAN_PARTICLE_CLASS_TYPE), POINTER :: LPC
LOGICAL :: FOUND

IF (EVACUATION) RETURN

M=>MESHES(NM)

!NOT_EVACUATION_MESH_IF: IF (.NOT.EVACUATION) THEN

ALLOCATE(MESHES(NM)%OMESH(NMESHES))

ALLOCATE(M%CONNECTED_MESH(NMESHES)) ; M%CONNECTED_MESH = .FALSE.

!END IF NOT_EVACUATION_MESH_IF

OTHER_MESH_LOOP: DO NOM=1,NMESHES

   !IF (EVACUATION) CYCLE OTHER_MESH_LOOP

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
      II = WC%BOUNDARY_COORD%II
      JJ = WC%BOUNDARY_COORD%JJ
      KK = WC%BOUNDARY_COORD%KK
      EWC%NIC_MIN = OM%NIC_R + 1
      EWC%NIC = (EWC%IIO_MAX-EWC%IIO_MIN+1)*(EWC%JJO_MAX-EWC%JJO_MIN+1)*(EWC%KKO_MAX-EWC%KKO_MIN+1)
      OM%NIC_R = OM%NIC_R + EWC%NIC
      EWC%NIC_MAX = OM%NIC_R
      FOUND = .TRUE.
      M%CONNECTED_MESH(NOM) = .TRUE.
      IOR = M%WALL(IW)%BOUNDARY_COORD%IOR
      SELECT CASE(IOR)
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

      SELECT CASE(ABS(IOR))
         CASE(1)
            EWC%AREA_RATIO = M%DY(JJ)*M%DZ(KK)/((M2%Y(EWC%JJO_MAX)-M2%Y(EWC%JJO_MIN-1))*(M2%Z(EWC%KKO_MAX)-M2%Z(EWC%KKO_MIN-1)))
         CASE(2)
            EWC%AREA_RATIO = M%DX(II)*M%DZ(KK)/((M2%X(EWC%IIO_MAX)-M2%X(EWC%IIO_MIN-1))*(M2%Z(EWC%KKO_MAX)-M2%Z(EWC%KKO_MIN-1)))
         CASE(3)
            EWC%AREA_RATIO = M%DX(II)*M%DY(JJ)/((M2%X(EWC%IIO_MAX)-M2%X(EWC%IIO_MIN-1))*(M2%Y(EWC%JJO_MAX)-M2%Y(EWC%JJO_MIN-1)))
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
                  IOR = M%WALL(IW)%BOUNDARY_COORD%IOR
                  OM%IIO_R(NIC_R) = IIO
                  OM%JJO_R(NIC_R) = JJO
                  OM%KKO_R(NIC_R) = KKO
                  OM%IOR_R(NIC_R) = IOR
               ENDDO
            ENDDO
         ENDDO
      ENDDO INDEX_LOOP

   ENDIF

   ! For PERIODIC boundaries with 1 or 2 meshes, we must revert to allocating whole copies of OMESH

   IF (IMIN>IMAX) THEN; IMIN=0; IMAX=M2%IBP1; ENDIF
   IF (JMIN>JMAX) THEN; JMIN=0; JMAX=M2%JBP1; ENDIF
   IF (KMIN>KMAX) THEN; KMIN=0; KMAX=M2%KBP1; ENDIF

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

   IF (SOLID_HT3D) THEN
      ALLOCATE(OM%TMP(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
      OM%TMP = TMPA
   ENDIF

   IF (LEVEL_SET_MODE>0 .OR. TERRAIN_CASE) THEN
      ALLOCATE(OM%PHI_LS(IMIN:IMAX,JMIN:JMAX))  ; OM%PHI_LS   = -1._EB
      ALLOCATE(OM%PHI1_LS(IMIN:IMAX,JMIN:JMAX)) ; OM%PHI1_LS  = -1._EB
      ALLOCATE(OM%U_LS(IMIN:IMAX,JMIN:JMAX))    ; OM%U_LS     =  0._EB
      ALLOCATE(OM%V_LS(IMIN:IMAX,JMIN:JMAX))    ; OM%V_LS     =  0._EB
      ALLOCATE(OM%Z_LS(IMIN:IMAX,JMIN:JMAX))    ; OM%Z_LS     =  0._EB
   ENDIF

ENDDO OTHER_MESH_LOOP

! Allocate arrays that hold data for particles that leave mesh NM, enter a
! NEIGHBORING_MESH as ORPHANs and then get ADOPTed by that mesh.

IF (OMESH_PARTICLES) THEN

   NEIGHBORING_MESH_LOOP: DO N=1,M%N_NEIGHBORING_MESHES

      OM => M%OMESH(M%NEIGHBORING_MESH(N))

      ALLOCATE(OM%N_PART_ORPHANS(N_LAGRANGIAN_CLASSES))
      ALLOCATE(OM%N_PART_ADOPT(N_LAGRANGIAN_CLASSES))
      OM%N_PART_ORPHANS = 0
      OM%N_PART_ADOPT   = 0
      ALLOCATE(OM%ORPHAN_PARTICLE_STORAGE(N_LAGRANGIAN_CLASSES))
      ALLOCATE(OM%ADOPT_PARTICLE_STORAGE(N_LAGRANGIAN_CLASSES))
      DO ILC=1,N_LAGRANGIAN_CLASSES
         LPC => LAGRANGIAN_PARTICLE_CLASS(ILC)
         IF (LPC%STATIC .AND. .NOT.LPC%EMBER_PARTICLE) CYCLE
         N_STORAGE_SLOTS = 10
         OM%ORPHAN_PARTICLE_STORAGE(ILC)%N_STORAGE_SLOTS = N_STORAGE_SLOTS
         OM%ADOPT_PARTICLE_STORAGE(ILC)%N_STORAGE_SLOTS = N_STORAGE_SLOTS
         ALLOCATE(OM%ORPHAN_PARTICLE_STORAGE(ILC)%REALS(LPC%N_STORAGE_REALS,N_STORAGE_SLOTS))
         ALLOCATE(OM%ORPHAN_PARTICLE_STORAGE(ILC)%INTEGERS(LPC%N_STORAGE_INTEGERS,N_STORAGE_SLOTS))
         ALLOCATE(OM%ORPHAN_PARTICLE_STORAGE(ILC)%LOGICALS(LPC%N_STORAGE_LOGICALS,N_STORAGE_SLOTS))
         ALLOCATE(OM%ADOPT_PARTICLE_STORAGE(ILC)%REALS(LPC%N_STORAGE_REALS,N_STORAGE_SLOTS))
         ALLOCATE(OM%ADOPT_PARTICLE_STORAGE(ILC)%INTEGERS(LPC%N_STORAGE_INTEGERS,N_STORAGE_SLOTS))
         ALLOCATE(OM%ADOPT_PARTICLE_STORAGE(ILC)%LOGICALS(LPC%N_STORAGE_LOGICALS,N_STORAGE_SLOTS))
      ENDDO

   ENDDO NEIGHBORING_MESH_LOOP

ENDIF

END SUBROUTINE INITIALIZE_MESH_EXCHANGE_1


SUBROUTINE INITIALIZE_MESH_EXCHANGE_2(NM)

! Create arrays by which info is to exchanged across meshes. In this routine, allocate arrays that involve NIC_R and NIC_S arrays.

INTEGER :: NOM
INTEGER, INTENT(IN) :: NM
TYPE (MESH_TYPE), POINTER :: M

IF (EVACUATION) RETURN

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


!> \brief Bordering meshes tell their neighbors how many exposed back wall cells they expect information for.

SUBROUTINE INITIALIZE_BACK_WALL_EXCHANGE

INTEGER :: NOM

CALL POST_RECEIVES(8)
CALL MESH_EXCHANGE(8)

! DEFINITION MESHES(NM)%OMESH(NOM)%N_WALL_CELLS_SEND
! Number of wall cells in Mesh NM for which information must be sent to Mesh NOM.
! DEFINITION MESHES(NM)%OMESH(NOM)%WALL_CELL_INDICES_SEND
! Indices of the wall cells in Mesh NM for which information needs to be sent to Mesh NOM.

DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   IF (EVACUATION) CYCLE
   DO NOM=1,NMESHES
      M3 => MESHES(NM)%OMESH(NOM)
      IF (M3%N_WALL_CELLS_SEND>0) ALLOCATE(M3%WALL_CELL_INDICES_SEND(M3%N_WALL_CELLS_SEND))
   ENDDO
ENDDO

! Mesh NM sends MESHES(NM)%OMESH(NOM)%EXPOSED_WALL_CELL_BACK_INDICES to Mesh NOM where it is received into
! MESHES(NOM)%OMESH(NM)%WALL_CELL_INDICES_SEND

CALL POST_RECEIVES(9)
CALL MESH_EXCHANGE(9)

! Set up arrays to send and receive exposed back wall cell information.

MESH_LOOP_1: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   IF (EVACUATION) CYCLE MESH_LOOP_1
   M => MESHES(NM)
   MESH_LOOP_2: DO NOM=1,NMESHES
      IF (NM==NOM) CYCLE MESH_LOOP_2
      M3 => M%OMESH(NOM)
      IF (M3%N_WALL_CELLS_SEND>0) ALLOCATE(M3%REAL_SEND_PKG6(M3%N_WALL_CELLS_SEND*2))
      IF (M3%N_EXPOSED_WALL_CELLS>0) THEN
         ALLOCATE(M3%REAL_RECV_PKG6(M3%N_EXPOSED_WALL_CELLS*2))
         ALLOCATE(M3%EXPOSED_WALL(M3%N_EXPOSED_WALL_CELLS))
      ENDIF
   ENDDO MESH_LOOP_2
ENDDO MESH_LOOP_1

! Check to see if any process has an error. If so, stop the run.

CALL STOP_CHECK(1)

! Set up persistent SEND and RECV calls for BACK_WALL info

CALL POST_RECEIVES(10)
CALL MESH_EXCHANGE(10)

END SUBROUTINE INITIALIZE_BACK_WALL_EXCHANGE


!> \brief Initialize the array PRESSURE_ZONE(I,J,K) for all meshes

SUBROUTINE INITIALIZE_PRESSURE_ZONES

USE GEOMETRY_FUNCTIONS, ONLY: ASSIGN_PRESSURE_ZONE
TYPE(P_ZONE_TYPE), POINTER :: PZ
INTEGER :: N,N_OVERLAP,I,J,K,N_EVAC_ZONE,N_EVAC_MESH,NM_EVAC

LOGICAL :: FOUND_UNASSIGNED_CELL

! For each explicitly specified ZONE, populate the array PRESSURE_ZONE(I,J,K) with the ZONE index, N

MESH_LOOP_1: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   ZONE_LOOP: DO N=1,N_ZONE
      PZ => P_ZONE(N)
      IF (EVACUATION .OR. PZ%EVACUATION) CYCLE ZONE_LOOP
      CALL ASSIGN_PRESSURE_ZONE(NM,PZ%X,PZ%Y,PZ%Z,N,N_OVERLAP)
      IF (N_OVERLAP>0) THEN
         WRITE(LU_ERR,'(A,I2,A,I2,A,I4)') 'ERROR: ZONE ',N,' overlaps ZONE ',N_OVERLAP,' in MESH ',NM
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
   IF (EVACUATION) CYCLE MESH_LOOP_2
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
                  IF (M%PRESSURE_ZONE(I,J,K)<0 .AND. .NOT.M%SOLID(M%CELL_INDEX(I,J,K))) THEN
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
         IF (MY_RANK==PROCESS(NM)) CALL ASSIGN_PRESSURE_ZONE(NM,M%XC(I),M%YC(J),M%ZC(K),N_ZONE,N_OVERLAP)
         CALL ZONE_BOUNDARY_EXCHANGE
      ENDIF
   ENDDO SEARCH_LOOP
ENDDO MESH_LOOP_2

! Number of evacuation zones = number of main evacuation meshes with flow fields
N_EVAC_ZONE = 0
DO NM=1,NMESHES
   IF (EVACUATION .AND. EVACUATION_SKIP(NM)) N_EVAC_ZONE = N_EVAC_ZONE + 1
ENDDO

EVACUATION_ZONE_IF: IF (N_EVAC_ZONE>0) THEN
   CALL REALLOCATE_P_ZONE(N_ZONE,N_ZONE+N_EVAC_ZONE)
   N_ZONE = N_ZONE + N_EVAC_ZONE

   EVACUATION_ZONE_LOOP: DO N=N_ZONE-N_EVAC_ZONE+1,N_ZONE

      ! Mapping: main evacuation mesh index/p_zone index/mesh index
      ! Each main evacuation mesh has exactly one p_zone

      N_EVAC_MESH = 0
      NM_EVAC = 0
      EVAC_MESH_NUMBER: DO NM=1,NMESHES
         IF (EVACUATION_SKIP(NM)) THEN
            N_EVAC_MESH = N_EVAC_MESH + 1
            NM_EVAC = NM
            IF (N_EVAC_MESH == N-(N_ZONE-N_EVAC_ZONE)) EXIT EVAC_MESH_NUMBER
         END IF
      ENDDO EVAC_MESH_NUMBER

      P_ZONE(N)%MESH_INDEX = NM_EVAC

      P_ZONE(N)%ID = TRIM('EvacZONE_') // TRIM(MESH_NAME(NM_EVAC))
      P_ZONE(N)%EVACUATION = .TRUE.
      P_ZONE(N)%PERIODIC = .FALSE.
   ENDDO EVACUATION_ZONE_LOOP

ENDIF EVACUATION_ZONE_IF

! Ensure that all cells have been assigned a pressure zone, even within solids

DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   M => MESHES(NM)
   WHERE (M%PRESSURE_ZONE<0) M%PRESSURE_ZONE = 0
ENDDO

END SUBROUTINE INITIALIZE_PRESSURE_ZONES


!> \brief Propagate PRESSURE_ZONE values from the mesh of origin to neighboring meshes

SUBROUTINE ZONE_BOUNDARY_EXCHANGE

USE GEOMETRY_FUNCTIONS, ONLY: ASSIGN_PRESSURE_ZONE
INTEGER :: N_OVERLAP,IW,II,JJ,KK,IIG,JJG,KKG,IIO,JJO,KKO,NOM,NM
TYPE (MESH_TYPE), POINTER :: OM
TYPE (WALL_TYPE), POINTER :: WC
TYPE (EXTERNAL_WALL_TYPE), POINTER :: EWC

SETUP_PRESSURE_ZONES_INDEX = 0  ! Flag to indicate that the PRESSURE_ZONE values no longer need to be spread

DO WHILE (ANY(SETUP_PRESSURE_ZONES_INDEX==0))

   ! Broadcast the PRESSURE_ZONE array of each mesh to the MPI processes controlling neighboring meshes

   IF (N_MPI_PROCESSES>1) THEN
      DO NM=1,NMESHES
         IF (MPI_COMM_MESH(NM)==MPI_COMM_NULL) CYCLE
         M => MESHES(NM)
         CALL MPI_BCAST(M%PRESSURE_ZONE,SIZE(M%PRESSURE_ZONE),MPI_INTEGER,MPI_COMM_MESH_ROOT(NM),MPI_COMM_MESH(NM),IERR)
      ENDDO
   ENDIF

   ! For each mesh, go around the exterior boundary looking for OPEN or INTERPOLATED boundaries.
   ! For the OPEN boundaries, propagate a PRESSURE_ZONE value of 0 into the interior of the mesh.
   ! For the INTERPOLATED boundaries (NOM/=0), get the PRESSURE_ZONE value from the neighboring mesh and propagate it internally.

   MESH_LOOP: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

      M => MESHES(NM)

      SETUP_PRESSURE_ZONES_INDEX(NM) = 1  ! When this flag remains 1 for all meshes, stop iterating

      WALL_LOOP: DO IW=1,M%N_EXTERNAL_WALL_CELLS
         IF (EVACUATION) CYCLE WALL_LOOP
         WC  => M%WALL(IW)
         EWC => M%EXTERNAL_WALL(IW)
         NOM = EWC%NOM
         IIG = WC%BOUNDARY_COORD%IIG
         JJG = WC%BOUNDARY_COORD%JJG
         KKG = WC%BOUNDARY_COORD%KKG
         II  = WC%BOUNDARY_COORD%II
         JJ  = WC%BOUNDARY_COORD%JJ
         KK  = WC%BOUNDARY_COORD%KK
         M%PRESSURE_ZONE(II,JJ,KK) = 0
         IF (WC%BOUNDARY_TYPE==OPEN_BOUNDARY .AND. M%PRESSURE_ZONE(IIG,JJG,KKG)<0) THEN
            CALL ASSIGN_PRESSURE_ZONE(NM,M%XC(IIG),M%YC(JJG),M%ZC(KKG),0,N_OVERLAP)
            IF (N_OVERLAP>=0) THEN
               WRITE(LU_ERR,'(A,I2,A,I2,A,I4)') 'ERROR: ZONE ',0,' overlaps ZONE ',N_OVERLAP,' in MESH ',NM
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
               WRITE(LU_ERR,'(A,I2,A,I2,A,I4,A,I4)') 'ERROR: ZONE ',OM%PRESSURE_ZONE(IIO,JJO,KKO),' meets ZONE ',&
                  M%PRESSURE_ZONE(IIG,JJG,KKG),' at the boundary of MESH ',NOM,' and MESH ',NM
               STOP_STATUS = SETUP_STOP
               EXIT MESH_LOOP
            ENDIF
            IF (WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .AND. M%PRESSURE_ZONE(IIG,JJG,KKG)<0 .AND. &
                OM%PRESSURE_ZONE(IIO,JJO,KKO)>=0) THEN
               CALL ASSIGN_PRESSURE_ZONE(NM,M%XC(IIG),M%YC(JJG),M%ZC(KKG),OM%PRESSURE_ZONE(IIO,JJO,KKO),N_OVERLAP)
               IF (N_OVERLAP>0) THEN
                  WRITE(LU_ERR,'(A,I2,A,I2,A,I4)') &
                     'ERROR: ZONE ',OM%PRESSURE_ZONE(IIO,JJO,KKO),' overlaps ZONE ',N_OVERLAP,' in MESH ',NM
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

   CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,SETUP_PRESSURE_ZONES_INDEX,COUNTS,DISPLS,MPI_INTEGER,MPI_COMM_WORLD,IERR)
ENDDO

IF (MY_RANK==0 .AND. VERBOSE) WRITE(LU_ERR,'(A)') ' Completed SETUP_PRESSURE_ZONES'

END SUBROUTINE ZONE_BOUNDARY_EXCHANGE


!> \brief Set up receive buffers for MPI calls.

SUBROUTINE POST_RECEIVES(CODE)

INTEGER, INTENT(IN) :: CODE
INTEGER :: RNODE,SNODE,IJK_SIZE,N,N_STORAGE_SLOTS,NRA,NRA_MAX,LL,AIC,NN,NOM
REAL(EB) :: TNOW
TYPE (LAGRANGIAN_PARTICLE_CLASS_TYPE), POINTER :: LPC

TNOW = CURRENT_TIME()

! Initialize the number of non-persistent send/receive requests.

N_REQ = 0

! Loop over all receive meshes (NM) and look for the send meshes (NOM).

MESH_LOOP: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   IF (EVACUATION) CYCLE MESH_LOOP

   RNODE = PROCESS(NM)
   M => MESHES(NM)

   OTHER_MESH_LOOP: DO NN=1,M%N_NEIGHBORING_MESHES

      NOM = M%NEIGHBORING_MESH(NN)
      M3=>MESHES(NM)%OMESH(NOM)
      IF (NOM==NM .AND. M3%NIC_R==0) CYCLE OTHER_MESH_LOOP

      SNODE = PROCESS(NOM)
      IF (RNODE==SNODE) CYCLE OTHER_MESH_LOOP

      M4=>MESHES(NOM)

      ! Set up receives for one-time exchanges or persistent send/receives.

      INITIALIZATION_IF: IF (CODE==0) THEN

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

            ALLOCATE(M3%REAL_RECV_PKG1(M3%NIC_R*(6+2*N_TOTAL_SCALARS)))
            ALLOCATE(M3%REAL_RECV_PKG3(IJK_SIZE*4))
            ALLOCATE(M3%REAL_RECV_PKG5(NRA_MAX*NUMBER_SPECTRAL_BANDS*M3%NIC_R))
            ALLOCATE(M3%REAL_RECV_PKG7(M3%NIC_R*3))
            ALLOCATE(M3%REAL_RECV_PKG8(M3%NIC_R*2))

            IF (SOLID_HT3D) ALLOCATE(M3%REAL_RECV_PKG4(M3%NIC_R*2))

            IF (LEVEL_SET_MODE>0 .OR. TERRAIN_CASE) ALLOCATE(M3%REAL_RECV_PKG14(4*M3%NIC_R))

         ENDIF

         ! Set up persistent receive requests

         IF (OMESH_PARTICLES) THEN
            N_REQ2 = N_REQ2 + 1
            CALL MPI_RECV_INIT(M3%N_PART_ADOPT,SIZE(M3%N_PART_ADOPT),MPI_INTEGER,SNODE,NOM,MPI_COMM_WORLD,REQ2(N_REQ2),IERR)
         ENDIF

         IF (M3%NIC_R>0) THEN

            N_REQ1 = N_REQ1 + 1
            CALL MPI_RECV_INIT(M3%REAL_RECV_PKG1(1),SIZE(M3%REAL_RECV_PKG1),MPI_DOUBLE_PRECISION,SNODE,NOM,MPI_COMM_WORLD,&
                               REQ1(N_REQ1),IERR)

            N_REQ15 = N_REQ15 + 1
            CALL MPI_RECV_INIT(M3%N_INTERNAL_OBST,INTEGER_ONE,MPI_INTEGER,SNODE,NOM,MPI_COMM_WORLD,REQ15(N_REQ15),IERR)

            N_REQ3 = N_REQ3 + 1
            CALL MPI_RECV_INIT(M3%REAL_RECV_PKG3(1),SIZE(M3%REAL_RECV_PKG3),MPI_DOUBLE_PRECISION,SNODE,NOM,MPI_COMM_WORLD,&
                               REQ3(N_REQ3),IERR)

            IF (SOLID_HT3D) THEN
               N_REQ4 = N_REQ4 + 1
               CALL MPI_RECV_INIT(M3%REAL_RECV_PKG4(1),SIZE(M3%REAL_RECV_PKG4),MPI_DOUBLE_PRECISION,SNODE,NOM,MPI_COMM_WORLD,&
                                  REQ4(N_REQ4),IERR)
            ENDIF

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
         N_REQ=MIN(N_REQ+1,SIZE(REQ))
         CALL MPI_IRECV(M3%N_WALL_CELLS_SEND,1,MPI_INTEGER,SNODE,NOM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
      ENDIF

      IF (CODE==9 .AND. M3%N_WALL_CELLS_SEND>0) THEN
            N_REQ=MIN(N_REQ+1,SIZE(REQ))
            CALL MPI_IRECV(M3%WALL_CELL_INDICES_SEND,SIZE(M3%WALL_CELL_INDICES_SEND),MPI_INTEGER,SNODE,NOM,MPI_COMM_WORLD,&
                           REQ(N_REQ),IERR)
      ENDIF

      IF (CODE==10 .AND. M3%N_EXPOSED_WALL_CELLS>0) THEN
         N_REQ6 = N_REQ6 + 1
         CALL MPI_RECV_INIT(M3%REAL_RECV_PKG6(1),SIZE(M3%REAL_RECV_PKG6),MPI_DOUBLE_PRECISION,SNODE,NOM,MPI_COMM_WORLD,&
                            REQ6(N_REQ6),IERR)
      ENDIF

      ! PARTICLEs

      IF (CODE==11 .AND. OMESH_PARTICLES) THEN
         DO N=1,N_LAGRANGIAN_CLASSES
            IF (M3%N_PART_ADOPT(N)==0) CYCLE
            LPC => LAGRANGIAN_PARTICLE_CLASS(N)
            IF (LPC%STATIC .AND. .NOT.LPC%EMBER_PARTICLE) CYCLE
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

      ! Set up to receive mass loss data for boundary OBSTs

      IF (CODE==16 .AND. M3%N_INTERNAL_OBST>0) THEN
         N_REQ = N_REQ + 1
         CALL MPI_IRECV(M3%REAL_RECV_PKG8(1),2*M3%N_INTERNAL_OBST,MPI_DOUBLE_PRECISION,SNODE,NOM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
      ENDIF
      IF (CODE==18 .AND. M3%N_EXTERNAL_OBST>0) THEN
         N_REQ = N_REQ + 1
         CALL MPI_IRECV(M3%REAL_RECV_PKG8(1),2*M3%N_EXTERNAL_OBST,MPI_DOUBLE_PRECISION,SNODE,NOM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
      ENDIF

   ENDDO OTHER_MESH_LOOP

ENDDO MESH_LOOP

! Receive EVACuation information

OTHER_MESH_LOOP_EVAC: DO NOM=1,NMESHES
   SNODE = PROCESS(NOM)
   RNODE = MAX(0,EVAC_PROCESS)
   IF (RNODE==SNODE) CYCLE OTHER_MESH_LOOP_EVAC
   IF (CODE==6 .AND. EXCHANGE_EVACUATION .AND. MY_RANK==MAX(0,EVAC_PROCESS) .AND. .NOT.EVACUATION) THEN
      M4=>MESHES(NOM)
      TAG_EVAC = NOM*(MAX(0,EVAC_PROCESS+1)+1)*CODE*10
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
ENDDO OTHER_MESH_LOOP_EVAC

T_USED(11)=T_USED(11) + CURRENT_TIME() - TNOW
END SUBROUTINE POST_RECEIVES


SUBROUTINE MESH_EXCHANGE(CODE)

! Exchange Information between Meshes

REAL(EB) :: TNOW
INTEGER, INTENT(IN) :: CODE
INTEGER :: NM,NOM,II,JJ,KK,LL,LLL,N,RNODE,SNODE,IMIN,IMAX,JMIN,JMAX,KMIN,KMAX,IJK_SIZE,N_STORAGE_SLOTS,N_NEW_STORAGE_SLOTS
INTEGER :: NNN,NN1,NN2,IPC,CNT,STORAGE_INDEX_SAVE,II1,II2,JJ1,JJ2,KK1,KK2,NQT2,NN,IOR,NRA,NRA_MAX,AIC,OBST_INDEX
REAL(EB), POINTER, DIMENSION(:,:) :: PHI_LS_P
REAL(EB), POINTER, DIMENSION(:,:,:) :: HP,HP2,RHOP,RHOP2,DP,DP2,UP,UP2,VP,VP2,WP,WP2
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZP,ZZP2
REAL(EB) :: XI,YJ,ZK
TYPE (LAGRANGIAN_PARTICLE_TYPE), POINTER :: LP
TYPE (LAGRANGIAN_PARTICLE_CLASS_TYPE), POINTER :: LPC

IF(CC_IBM) CALL MESH_CC_EXCHANGE(CODE)

TNOW = CURRENT_TIME()

! Special circumstances when doing the radiation exchange (CODE=2)

IF (CODE==2 .AND. (.NOT.EXCHANGE_RADIATION .OR. .NOT.RADIATION)) RETURN

! For each mesh, NM, controlled by MPI process, SNODE, send data to other meshes, NOM.

SENDING_MESH_LOOP: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   IF (EVACUATION) CYCLE SENDING_MESH_LOOP
   M =>MESHES(NM)
   SNODE = PROCESS(NM)

   RECEIVING_MESH_LOOP: DO NNN=1,M%N_NEIGHBORING_MESHES

      NOM = M%NEIGHBORING_MESH(NNN)
      M3 => M%OMESH(NOM)
      IF (NOM==NM .AND. M3%NIC_S==0) CYCLE RECEIVING_MESH_LOOP

      RNODE = PROCESS(NOM)

      IMIN = M3%I_MIN_S
      IMAX = M3%I_MAX_S
      JMIN = M3%J_MIN_S
      JMAX = M3%J_MAX_S
      KMIN = M3%K_MIN_S
      KMAX = M3%K_MAX_S

      IJK_SIZE = (IMAX-IMIN+1)*(JMAX-JMIN+1)*(KMAX-KMIN+1)

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

            ALLOCATE(M3%REAL_SEND_PKG1(M3%NIC_S*(6+2*N_TOTAL_SCALARS)))
            ALLOCATE(M3%REAL_SEND_PKG3(IJK_SIZE*4))
            ALLOCATE(M3%REAL_SEND_PKG5(NRA_MAX*NUMBER_SPECTRAL_BANDS*M3%NIC_S))
            ALLOCATE(M3%REAL_SEND_PKG7(M3%NIC_S*3))

            IF (SOLID_HT3D) ALLOCATE(M3%REAL_SEND_PKG4(M3%NIC_S*2))

            IF (LEVEL_SET_MODE>0 .OR. TERRAIN_CASE) ALLOCATE(M3%REAL_SEND_PKG14(4*M3%NIC_S))

         ENDIF

         ! REAL_SEND_PKG8 carries obstruction mass loss information and is used even if RNODE=SNODE

         IF (M3%NIC_S>0) ALLOCATE(M3%REAL_SEND_PKG8(M3%NIC_S*2))

         ! Initialize persistent send requests

         IF (OMESH_PARTICLES) THEN
            N_REQ2 = N_REQ2 + 1
            CALL MPI_SEND_INIT(M3%N_PART_ORPHANS,SIZE(M3%N_PART_ORPHANS),MPI_INTEGER,RNODE,NM,MPI_COMM_WORLD,REQ2(N_REQ2),IERR)
         ENDIF

         IF (M3%NIC_S>0 .AND. RNODE/=SNODE) THEN

            N_REQ1 = N_REQ1 + 1
            CALL MPI_SEND_INIT(M3%REAL_SEND_PKG1(1),SIZE(M3%REAL_SEND_PKG1),MPI_DOUBLE_PRECISION,RNODE,NM,MPI_COMM_WORLD,&
                               REQ1(N_REQ1),IERR)

            N_REQ15 = N_REQ15 + 1
            CALL MPI_SEND_INIT(M3%N_EXTERNAL_OBST,INTEGER_ONE,MPI_INTEGER,RNODE,NM,MPI_COMM_WORLD,REQ15(N_REQ15),IERR)

            N_REQ3 = N_REQ3 + 1
            CALL MPI_SEND_INIT(M3%REAL_SEND_PKG3(1),SIZE(M3%REAL_SEND_PKG3),MPI_DOUBLE_PRECISION,RNODE,NM,MPI_COMM_WORLD,&
                               REQ3(N_REQ3),IERR)

            IF (SOLID_HT3D) THEN
               N_REQ4 = N_REQ4 + 1
               CALL MPI_SEND_INIT(M3%REAL_SEND_PKG4(1),SIZE(M3%REAL_SEND_PKG4),MPI_DOUBLE_PRECISION,RNODE,NM,MPI_COMM_WORLD,&
                                  REQ4(N_REQ4),IERR)
            ENDIF

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
            N_REQ=MIN(N_REQ+1,SIZE(REQ))
            CALL MPI_ISEND(M3%N_EXPOSED_WALL_CELLS,1,MPI_INTEGER,RNODE,NM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
         ELSE
            M2=>MESHES(NOM)%OMESH(NM)
            M2%N_WALL_CELLS_SEND = M3%N_EXPOSED_WALL_CELLS
         ENDIF
      ENDIF

      IF (CODE==9 .AND. M3%N_EXPOSED_WALL_CELLS>0) THEN
         IF (RNODE/=SNODE) THEN
            N_REQ=MIN(N_REQ+1,SIZE(REQ))
            CALL MPI_ISEND(M3%EXPOSED_WALL_CELL_BACK_INDICES,SIZE(M3%EXPOSED_WALL_CELL_BACK_INDICES),MPI_INTEGER,RNODE,NM,&
                           MPI_COMM_WORLD,REQ(N_REQ),IERR)
         ELSE
            M2=>MESHES(NOM)%OMESH(NM)
            M2%WALL_CELL_INDICES_SEND = M3%EXPOSED_WALL_CELL_BACK_INDICES
         ENDIF
      ENDIF

      IF (CODE==10 .AND. M3%N_WALL_CELLS_SEND>0) THEN
         IF (RNODE/=SNODE) THEN
            N_REQ6 = N_REQ6 + 1
            CALL MPI_SEND_INIT(M3%REAL_SEND_PKG6(1),SIZE(M3%REAL_SEND_PKG6),MPI_DOUBLE_PRECISION,RNODE,NM,MPI_COMM_WORLD,&
                                REQ6(N_REQ6),IERR)
         ENDIF
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
            DO KK=KMIN,KMAX
               DO JJ=JMIN,JMAX
                  DO II=IMIN,IMAX
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

      IF (CODE==6) THEN
         IF (RNODE/=SNODE) THEN
            LL = 0
            DO II=1,M%OMESH(NOM)%N_WALL_CELLS_SEND
               IW = M%OMESH(NOM)%WALL_CELL_INDICES_SEND(II)
               M3%REAL_SEND_PKG6(LL+1) = M%WALL(IW)%ONE_D%Q_RAD_IN
               M3%REAL_SEND_PKG6(LL+2) = M%TMP(M%WALL(IW)%BOUNDARY_COORD%IIG, &
                                               M%WALL(IW)%BOUNDARY_COORD%JJG, &
                                               M%WALL(IW)%BOUNDARY_COORD%KKG)
               LL = LL+2
            ENDDO
         ELSE
            M2=>MESHES(NOM)%OMESH(NM)
            DO II=1,M2%N_EXPOSED_WALL_CELLS
               IW = M2%EXPOSED_WALL_CELL_BACK_INDICES(II)
               M2%EXPOSED_WALL(II)%Q_RAD_IN = M%WALL(IW)%ONE_D%Q_RAD_IN
               M2%EXPOSED_WALL(II)%TMP_GAS  = M%TMP(M%WALL(IW)%BOUNDARY_COORD%IIG, &
                                                    M%WALL(IW)%BOUNDARY_COORD%JJG, &
                                                    M%WALL(IW)%BOUNDARY_COORD%KKG)
            ENDDO
         ENDIF
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

      IF (CODE==7 .AND. OMESH_PARTICLES .AND. RNODE==SNODE) THEN
         M2=>MESHES(NOM)%OMESH(NM)
         M2%N_PART_ADOPT = M3%N_PART_ORPHANS
      ENDIF

      ! Sending/Receiving PARTICLE Buffer Arrays

      IF_SEND_PARTICLES: IF (CODE==11 .AND. OMESH_PARTICLES) THEN

         NODE_CHECK_PARTICLE: IF (SNODE/=RNODE) THEN

            DO IPC=1,N_LAGRANGIAN_CLASSES
               IF (M3%N_PART_ORPHANS(IPC)==0) CYCLE
               LPC => LAGRANGIAN_PARTICLE_CLASS(IPC)
               IF (LPC%STATIC .AND. .NOT.LPC%EMBER_PARTICLE) CYCLE

               N_STORAGE_SLOTS = M3%ORPHAN_PARTICLE_STORAGE(IPC)%N_STORAGE_SLOTS
               N_REQ=MIN(N_REQ+1,SIZE(REQ))
               CALL MPI_ISEND(M3%ORPHAN_PARTICLE_STORAGE(IPC)%REALS(1,1),LPC%N_STORAGE_REALS*N_STORAGE_SLOTS,MPI_DOUBLE_PRECISION, &
                              RNODE,NM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
               N_REQ=MIN(N_REQ+1,SIZE(REQ))
               CALL MPI_ISEND(M3%ORPHAN_PARTICLE_STORAGE(IPC)%INTEGERS(1,1),LPC%N_STORAGE_INTEGERS*N_STORAGE_SLOTS,MPI_INTEGER, &
                              RNODE,NM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
               N_REQ=MIN(N_REQ+1,SIZE(REQ))
               CALL MPI_ISEND(M3%ORPHAN_PARTICLE_STORAGE(IPC)%LOGICALS(1,1),LPC%N_STORAGE_LOGICALS*N_STORAGE_SLOTS,MPI_LOGICAL, &
                              RNODE,NM,MPI_COMM_WORLD,REQ(N_REQ),IERR)
            ENDDO

         ELSE NODE_CHECK_PARTICLE

            M2 => MESHES(NOM)%OMESH(NM)

            DO IPC=1,N_LAGRANGIAN_CLASSES
               LPC => LAGRANGIAN_PARTICLE_CLASS(IPC)
               IF (LPC%STATIC .AND. .NOT.LPC%EMBER_PARTICLE) CYCLE
               M2%ADOPT_PARTICLE_STORAGE(IPC)%REALS    = M3%ORPHAN_PARTICLE_STORAGE(IPC)%REALS
               M2%ADOPT_PARTICLE_STORAGE(IPC)%INTEGERS = M3%ORPHAN_PARTICLE_STORAGE(IPC)%INTEGERS
               M2%ADOPT_PARTICLE_STORAGE(IPC)%LOGICALS = M3%ORPHAN_PARTICLE_STORAGE(IPC)%LOGICALS
            ENDDO

         ENDIF NODE_CHECK_PARTICLE

         M3%N_PART_ORPHANS(1:N_LAGRANGIAN_CLASSES) = 0

      ENDIF IF_SEND_PARTICLES

      ! Send HT3D temperatures

      IF ((CODE==1.OR.CODE==4) .AND. M3%NIC_S>0 .AND. SOLID_HT3D) THEN
         IF (RNODE/=SNODE) THEN
            PACK_REAL_SEND_PKG4: DO LL=1,M3%NIC_S
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
               M3%REAL_SEND_PKG4(2*(LL-1)+1) = M%TMP(II1,JJ1,KK1)
               M3%REAL_SEND_PKG4(2*(LL-1)+2) = M%TMP(II2,JJ2,KK2)
            ENDDO PACK_REAL_SEND_PKG4
         ELSE
            M2=>MESHES(NOM)%OMESH(NM)
            M2%TMP(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)= M%TMP(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
         ENDIF
      ENDIF

      ! Send LEVEL_SET boundary values

      IF (CODE==14 .AND. M3%NIC_S>0) THEN
         IF (PREDICTOR) THEN
            PHI_LS_P => M%PHI1_LS
         ELSE
            PHI_LS_P => M%PHI_LS
         ENDIF
         IF (RNODE/=SNODE) THEN
            NQT2 = 4
            PACK_REAL_SEND_PKG14: DO LL=1,M3%NIC_S
               II1 = M3%IIO_S(LL)
               JJ1 = M3%JJO_S(LL)
               M3%REAL_SEND_PKG14(NQT2*(LL-1)+1) = PHI_LS_P(II1,JJ1)
               M3%REAL_SEND_PKG14(NQT2*(LL-1)+2) = M%U_LS(II1,JJ1)
               M3%REAL_SEND_PKG14(NQT2*(LL-1)+3) = M%V_LS(II1,JJ1)
               M3%REAL_SEND_PKG14(NQT2*(LL-1)+4) = M%Z_LS(II1,JJ1)
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

! Send information needed by EVACuation routine

MESH_LOOP_SEND_EVAC: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   SNODE = PROCESS(NM)
   RNODE = MAX(0,EVAC_PROCESS)
   IF (RNODE==SNODE) CYCLE MESH_LOOP_SEND_EVAC
   IF (CODE==6 .AND. EXCHANGE_EVACUATION.AND.MY_RANK/=MAX(0,EVAC_PROCESS) .AND..NOT.EVACUATION) THEN
      M => MESHES(NM)
      TAG_EVAC = NM*(MAX(0,EVAC_PROCESS+1)+1)*CODE*10
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
ENDDO MESH_LOOP_SEND_EVAC


! Halt communications until all processes are ready to receive the data.

IF (N_MPI_PROCESSES>1 .AND. CODE/=1 .AND. CODE/=3 .AND. CODE/=4 .AND. CODE/=5 .AND. N_REQ>0) THEN
   CALL TIMEOUT('REQ',N_REQ,REQ(1:N_REQ))
ENDIF

IF (N_MPI_PROCESSES>1 .AND. (CODE==1.OR.CODE==4) .AND. N_REQ1>0) THEN
   CALL MPI_STARTALL(N_REQ1,REQ1(1:N_REQ1),IERR)
   CALL TIMEOUT('REQ1',N_REQ1,REQ1(1:N_REQ1))
ENDIF

IF (N_MPI_PROCESSES>1 .AND. CODE==7 .AND. OMESH_PARTICLES .AND. N_REQ2>0) THEN
   CALL MPI_STARTALL(N_REQ2,REQ2(1:N_REQ2),IERR)
   CALL TIMEOUT('REQ2',N_REQ2,REQ2(1:N_REQ2))
ENDIF

IF (N_MPI_PROCESSES>1 .AND. (CODE==3.OR.CODE==6) .AND. N_REQ3>0) THEN
   CALL MPI_STARTALL(N_REQ3,REQ3(1:N_REQ3),IERR)
   CALL TIMEOUT('REQ3',N_REQ3,REQ3(1:N_REQ3))
ENDIF

IF (N_MPI_PROCESSES>1 .AND. (CODE==1.OR.CODE==4) .AND. N_REQ4>0 .AND. SOLID_HT3D) THEN
   CALL MPI_STARTALL(N_REQ4,REQ4(1:N_REQ4),IERR)
   CALL TIMEOUT('REQ4',N_REQ4,REQ4(1:N_REQ4))
ENDIF

IF (N_MPI_PROCESSES>1 .AND. CODE==5 .AND. N_REQ7>0) THEN
   CALL MPI_STARTALL(N_REQ7,REQ7(1:N_REQ7),IERR)
   CALL TIMEOUT('REQ7',N_REQ7,REQ7(1:N_REQ7))
ENDIF

IF (N_MPI_PROCESSES>1 .AND. CODE==6 .AND. N_REQ6>0) THEN
   CALL MPI_STARTALL(N_REQ6,REQ6(1:N_REQ6),IERR)
   CALL TIMEOUT('REQ6',N_REQ6,REQ6(1:N_REQ6))
ENDIF

IF (N_MPI_PROCESSES>1 .AND. CODE==2 .AND. N_REQ5>0) THEN
   CALL MPI_STARTALL(N_REQ5,REQ5(1:N_REQ5),IERR)
   CALL TIMEOUT('REQ5',N_REQ5,REQ5(1:N_REQ5))
ENDIF

IF (N_MPI_PROCESSES>1 .AND. CODE==14 .AND. N_REQ14>0) THEN
   CALL MPI_STARTALL(N_REQ14,REQ14(1:N_REQ14),IERR)
   CALL TIMEOUT('REQ14',N_REQ14,REQ14(1:N_REQ14))
ENDIF

IF (N_MPI_PROCESSES>1 .AND. CODE==15 .AND. N_REQ15>0) THEN
   CALL MPI_STARTALL(N_REQ15,REQ15(1:N_REQ15),IERR)
   CALL TIMEOUT('REQ15',N_REQ15,REQ15(1:N_REQ15))
ENDIF

! For each mesh, NM, controlled by the current MPI process, RNODE, loop over all
! other meshes, NOM, and load data received into the appropriate arrays.

RECV_MESH_LOOP: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   IF (EVACUATION) CYCLE RECV_MESH_LOOP
   M => MESHES(NM)
   RNODE = PROCESS(NM)

   SEND_MESH_LOOP: DO NNN=1,M%N_NEIGHBORING_MESHES

      NOM = M%NEIGHBORING_MESH(NNN)
      M2 => M%OMESH(NOM)
      IF (NOM==NM .AND. M2%NIC_R==0) CYCLE SEND_MESH_LOOP

      SNODE = PROCESS(NOM)

      IMIN = M2%I_MIN_R
      IMAX = M2%I_MAX_R
      JMIN = M2%J_MIN_R
      JMAX = M2%J_MAX_R
      KMIN = M2%K_MIN_R
      KMAX = M2%K_MAX_R

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
         DO KK=KMIN,KMAX
            DO JJ=JMIN,JMAX
               DO II=IMIN,IMAX
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

      RECEIVE_BACK_WALL: IF (CODE==6 .AND. SNODE/=RNODE) THEN
         LL = 0
         DO II=1,M2%N_EXPOSED_WALL_CELLS
            M2%EXPOSED_WALL(II)%Q_RAD_IN = M2%REAL_RECV_PKG6(LL+1)
            M2%EXPOSED_WALL(II)%TMP_GAS  = M2%REAL_RECV_PKG6(LL+2)
            LL = LL+2
         ENDDO
      ENDIF RECEIVE_BACK_WALL

      ! Sending/Receiving PARTICLE Buffer Arrays

      IF (CODE==7 .AND. OMESH_PARTICLES) THEN
         DO IPC=1,N_LAGRANGIAN_CLASSES
            LPC => LAGRANGIAN_PARTICLE_CLASS(IPC)
            IF (LPC%STATIC .AND. .NOT.LPC%EMBER_PARTICLE) CYCLE
            IF (M2%N_PART_ADOPT(IPC)>M2%ADOPT_PARTICLE_STORAGE(IPC)%N_STORAGE_SLOTS) THEN
               N_NEW_STORAGE_SLOTS = M2%N_PART_ADOPT(IPC)-M2%ADOPT_PARTICLE_STORAGE(IPC)%N_STORAGE_SLOTS
               CALL REALLOCATE_STORAGE_ARRAYS(NM,3,IPC,N_NEW_STORAGE_SLOTS,NOM)
            ENDIF
         ENDDO
      ENDIF

      IF_RECEIVE_PARTICLES: IF (CODE==11 .AND. OMESH_PARTICLES) THEN

         DO IPC=1,N_LAGRANGIAN_CLASSES
            IF (M2%N_PART_ADOPT(IPC)==0) CYCLE
            LPC => LAGRANGIAN_PARTICLE_CLASS(IPC)
            IF (LPC%STATIC .AND. .NOT.LPC%EMBER_PARTICLE) CYCLE
            CNT = 0
            DO N=M%NLP+1,M%NLP+M2%N_PART_ADOPT(IPC)
               CNT = CNT + 1
               CALL ALLOCATE_STORAGE(NM,LAGRANGIAN_PARTICLE_CLASS(IPC)%SURF_INDEX,LPC_INDEX=IPC,LP_INDEX=N,TAG=-1)
               LP=>M%LAGRANGIAN_PARTICLE(N)
               STORAGE_INDEX_SAVE = LP%STORAGE_INDEX
               M%PARTICLE_STORAGE(IPC)%REALS(:,STORAGE_INDEX_SAVE)    = M2%ADOPT_PARTICLE_STORAGE(IPC)%REALS(:,CNT)
               M%PARTICLE_STORAGE(IPC)%INTEGERS(:,STORAGE_INDEX_SAVE) = M2%ADOPT_PARTICLE_STORAGE(IPC)%INTEGERS(:,CNT)
               LP%ARRAY_INDEX = N
               LP%STORAGE_INDEX = STORAGE_INDEX_SAVE
               M%PARTICLE_STORAGE(IPC)%LOGICALS(:,STORAGE_INDEX_SAVE) = M2%ADOPT_PARTICLE_STORAGE(IPC)%LOGICALS(:,CNT)
               CALL GET_IJK(LP%BOUNDARY_COORD%X,LP%BOUNDARY_COORD%Y,LP%BOUNDARY_COORD%Z,NM,XI,YJ,ZK,&
                            LP%BOUNDARY_COORD%IIG,LP%BOUNDARY_COORD%JJG,LP%BOUNDARY_COORD%KKG)
               IF (LP%INIT_INDEX>0) THEN
                  DO NN=1,N_DEVC
                     IF (DEVICE(NN)%INIT_ID==INITIALIZATION(LP%INIT_INDEX)%ID) THEN
                        DEVICE(NN)%LP_TAG = LP%TAG
                        DEVICE(NN)%PART_CLASS_INDEX = IPC
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
            M%NLP = M%NLP + M2%N_PART_ADOPT(IPC)
            M2%N_PART_ADOPT(IPC) = 0
         ENDDO

      ENDIF IF_RECEIVE_PARTICLES

      ! Unpack temperature (TMP) only for the case when 3D solid heat conduction is being performed

      IF ((CODE==1.OR.CODE==4) .AND. M2%NIC_R>0 .AND. RNODE/=SNODE .AND.  SOLID_HT3D) THEN
         UNPACK_REAL_RECV_PKG4: DO LL=1,M2%NIC_R
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
            M2%TMP(II1,JJ1,KK1) = M2%REAL_RECV_PKG4(2*(LL-1)+1)
            M2%TMP(II2,JJ2,KK2) = M2%REAL_RECV_PKG4(2*(LL-1)+2)
         ENDDO UNPACK_REAL_RECV_PKG4
      ENDIF

      IF (CODE==14 .AND. M2%NIC_R>0 .AND. RNODE/=SNODE) THEN
            NQT2 = 4
            IF (PREDICTOR) THEN
               PHI_LS_P => M2%PHI1_LS
            ELSE
               PHI_LS_P => M2%PHI_LS
            ENDIF
            UNPACK_REAL_RECV_PKG14: DO LL=1,M2%NIC_R
               II1 = M2%IIO_R(LL)
               JJ1 = M2%JJO_R(LL)
               PHI_LS_P(II1,JJ1) = M2%REAL_RECV_PKG14(NQT2*(LL-1)+1)
               M2%U_LS(II1,JJ1)  = M2%REAL_RECV_PKG14(NQT2*(LL-1)+2)
               M2%V_LS(II1,JJ1)  = M2%REAL_RECV_PKG14(NQT2*(LL-1)+3)
               M2%Z_LS(II1,JJ1)  = M2%REAL_RECV_PKG14(NQT2*(LL-1)+4)
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


SUBROUTINE TIMEOUT(RNAME,NR,RR)

REAL(EB) :: START_TIME,WAIT_TIME
INTEGER, INTENT(IN) :: NR
TYPE (MPI_REQUEST), DIMENSION(:) :: RR
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
         WRITE(LU_ERR,'(A,A,I6,A,A)') TRIM(RNAME),' timed out for MPI process ',MY_RANK,' running on ',PNAME(1:PNAMELEN)
         FLAG = .TRUE.
         STOP_STATUS = MPI_TIMEOUT_STOP
      ENDIF
   ENDDO

ELSE

   ! If PROFILING=T, do not do MPI_TESTALL because too many calls to this routine swamps the tracing and profiling.

   CALL MPI_WAITALL(NR,RR(1:NR),MPI_STATUSES_IGNORE,IERR)

ENDIF

END SUBROUTINE TIMEOUT


SUBROUTINE DUMP_TIMERS

! Write out the file CHID_cpu.csv containing the timing breakdown of each MPI process.

INTEGER, PARAMETER :: LINE_LENGTH = 5 + (N_TIMERS+1)*11
CHARACTER(LEN=LINE_LENGTH) :: LINE
CHARACTER(LEN=LINE_LENGTH), DIMENSION(0:N_MPI_PROCESSES-1) :: LINE_ARRAY
CHARACTER(30) :: FRMT

! T_USED(1) is the time spent in the main routine; i.e. the time not spent in a subroutine.

T_USED(1) = CURRENT_TIME() - T_USED(1) - SUM(T_USED(2:N_TIMERS))
WRITE(FRMT,'(A,I2.2,A)') '(I5,',N_TIMERS+1,'(",",ES10.3))'
WRITE(LINE,FRMT) MY_RANK,(T_USED(I),I=1,N_TIMERS),SUM(T_USED(1:N_TIMERS))

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
   WRITE(LU_CPU,'(A)') 'Rank,MAIN,DIVG,MASS,VELO,PRES,WALL,DUMP,PART,RADI,FIRE,COMM,EVAC,HVAC,GEOM,VEGE,Total T_USED (s)'
   DO N=0,N_MPI_PROCESSES-1
      WRITE(LU_CPU,'(A)') LINE_ARRAY(N)
   ENDDO
   CLOSE(LU_CPU)
ENDIF

END SUBROUTINE DUMP_TIMERS


SUBROUTINE WRITE_STRINGS

! Write character strings out to the .smv file

INTEGER :: N,NOM,N_STRINGS_DUM
CHARACTER(MESH_STRING_LENGTH), ALLOCATABLE, DIMENSION(:) :: STRING_DUM
REAL(EB) :: TNOW

TNOW = CURRENT_TIME()

! All meshes send their STRINGs to node 0

DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   IF (MY_RANK>0) THEN
      CALL MPI_SEND(MESHES(NM)%N_STRINGS,1,MPI_INTEGER,0,1,MPI_COMM_WORLD,IERR)
      IF (MESHES(NM)%N_STRINGS>0) CALL MPI_SEND(MESHES(NM)%STRING(1),MESHES(NM)%N_STRINGS*MESH_STRING_LENGTH,MPI_CHARACTER,0,NM, &
                                                MPI_COMM_WORLD,IERR)
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
            ALLOCATE(STRING_DUM(N_STRINGS_DUM))
            CALL MPI_RECV(STRING_DUM(1),N_STRINGS_DUM*MESH_STRING_LENGTH, &
            MPI_CHARACTER,PROCESS(NOM),NOM,MPI_COMM_WORLD,STATUS,IERR)
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

DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   MESHES(NM)%N_STRINGS = 0
ENDDO

T_USED(11) = T_USED(11) + CURRENT_TIME() - TNOW
END SUBROUTINE WRITE_STRINGS


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

! Gather heat release rates (Q_DOT), mass loss rates (M_DOT), etc., to MPI process 0.

IF (N_MPI_PROCESSES>1) CALL GATHER_Q_AND_M

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


SUBROUTINE EXCHANGE_GLOBAL_OUTPUTS

! Gather HRR, mass, and device data to node 0

USE EVAC, ONLY: N_DOORS, N_EXITS, N_ENTRYS, EVAC_DOORS, EVAC_EXITS, EVAC_ENTRYS, EMESH_INDEX
REAL(EB) :: TNOW
INTEGER :: NN,N,I_STATE,I,OP_INDEX,NM,DISP,DIM_FAC
TYPE (MPI_OP) :: MPI_OP_INDEX
TYPE(DEVICE_TYPE), POINTER :: DV
TYPE(SUBDEVICE_TYPE), POINTER :: SDV
LOGICAL :: NO_NEED_TO_RECV

TNOW = CURRENT_TIME()

IF (EVACUATION .AND. (ICYC<1 .AND. T>T_BEGIN)) RETURN ! No dumps at the evacuation initialization phase

DISP = DISPLS(MY_RANK)+1

! Gather HRR (Q_DOT) and mass loss rate (M_DOT) integrals to node 0

IF (T>=HRR_CLOCK .AND. N_MPI_PROCESSES>1) CALL GATHER_Q_AND_M

! Gather species mass integrals to node 0

IF (T>=MASS_CLOCK .AND. N_MPI_PROCESSES>1) THEN
   REAL_BUFFER_MASS = MASS_DT
   CALL MPI_GATHERV(REAL_BUFFER_MASS(0,DISP),COUNTS_MASS(MY_RANK),MPI_DOUBLE_PRECISION, &
                    MASS_DT,COUNTS_MASS,DISPLS_MASS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
ENDIF

! Exchange DEVICE parameters among meshes and dump out DEVICE info after first "gathering" data to node 0

EXCHANGE_DEVICE: IF (N_DEVC>0) THEN

   ! Exchange the CURRENT_STATE and PRIOR_STATE of each DEViCe

   STATE_LOC = .FALSE.  ! _LOC is a temporary array that holds the STATE value for the devices on each node
   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      DO N=1,N_DEVC
         DV => DEVICE(N)
         IF (DV%MESH==NM) THEN
            STATE_LOC(N)        = DV%CURRENT_STATE
            STATE_LOC(N+N_DEVC) = DV%PRIOR_STATE
         ENDIF
      ENDDO
   ENDDO
   IF (N_MPI_PROCESSES>1) THEN
      CALL MPI_ALLREDUCE(STATE_LOC(1),STATE_GLB(1),2*N_DEVC,MPI_LOGICAL,MPI_LXOR,MPI_COMM_WORLD,IERR)
   ELSE
      STATE_GLB = STATE_LOC
   ENDIF
   DO N=1,N_DEVC
      DV => DEVICE(N)
      DV%CURRENT_STATE = STATE_GLB(N)
      DV%PRIOR_STATE   = STATE_GLB(N+N_DEVC)
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
         CASE(1) ; TC_LOC  =  0._EB    ; MPI_OP_INDEX = MPI_SUM    ; DIM_FAC = 3
         CASE(2) ; TC2_LOC =  1.E10_EB ; MPI_OP_INDEX = MPI_MINLOC ; DIM_FAC = 1
         CASE(3) ; TC2_LOC = -1.E10_EB ; MPI_OP_INDEX = MPI_MAXLOC ; DIM_FAC = 1
      END SELECT
      DEVICE_LOOP_1: DO N=1,N_DEVC
         DV => DEVICE(N)
         IF (OP_INDEX==1 .AND. (DV%SPATIAL_STATISTIC(1:3)=='MIN' .OR. DV%SPATIAL_STATISTIC(1:3)=='MAX')) CYCLE
         IF (OP_INDEX==2 .AND.  DV%SPATIAL_STATISTIC(1:3)/='MIN') CYCLE
         IF (OP_INDEX==3 .AND.  DV%SPATIAL_STATISTIC(1:3)/='MAX') CYCLE
         DO NN=1,DV%N_SUBDEVICES
            SDV => DV%SUBDEVICE(NN)
            SELECT CASE(OP_INDEX)
               CASE(1)
                  TC_LOC(N)            = TC_LOC(N)          + SDV%VALUE_1
                  TC_LOC(N+N_DEVC)     = TC_LOC(N+N_DEVC)   + SDV%VALUE_2
                  TC_LOC(N+2*N_DEVC)   = TC_LOC(N+2*N_DEVC) + SDV%VALUE_3
               CASE(2)
                  IF (SDV%VALUE_1<TC2_LOC(1,N)) THEN
                     TC2_LOC(1,N) = SDV%VALUE_1
                     TC2_LOC(2,N) = SDV%VALUE_2
                  ENDIF
               CASE(3)
                  IF (SDV%VALUE_1>TC2_LOC(1,N)) THEN
                     TC2_LOC(1,N) = SDV%VALUE_1
                     TC2_LOC(2,N) = SDV%VALUE_2
                  ENDIF
            END SELECT
         ENDDO
      ENDDO DEVICE_LOOP_1
      IF (N_MPI_PROCESSES>1) THEN
         SELECT CASE(OP_INDEX)
            CASE(1) ; CALL MPI_ALLREDUCE(TC_LOC(1),TC_GLB(1),DIM_FAC*N_DEVC,MPI_DOUBLE_PRECISION,MPI_OP_INDEX,MPI_COMM_WORLD,IERR)
            CASE(2:3) ; CALL MPI_ALLREDUCE(TC2_LOC,TC2_GLB,N_DEVC,MPI_2DOUBLE_PRECISION,MPI_OP_INDEX,MPI_COMM_WORLD,IERR)
         END SELECT
      ELSE
         SELECT CASE(OP_INDEX)
            CASE(1)   ; TC_GLB = TC_LOC
            CASE(2:3) ; TC2_GLB = TC2_LOC
         END SELECT
      ENDIF
      DEVICE_LOOP_2: DO N=1,N_DEVC
         DV => DEVICE(N)
         IF (OP_INDEX==1 .AND. (DV%SPATIAL_STATISTIC(1:3)=='MIN' .OR. DV%SPATIAL_STATISTIC(1:3)=='MAX')) CYCLE
         IF (OP_INDEX==2 .AND.  DV%SPATIAL_STATISTIC(1:3)/='MIN') CYCLE
         IF (OP_INDEX==3 .AND.  DV%SPATIAL_STATISTIC(1:3)/='MAX') CYCLE
         IF (OP_INDEX==1) THEN
            DV%VALUE_1 = TC_GLB(N)
            DV%VALUE_2 = TC_GLB(  N_DEVC+N)
            DV%VALUE_3 = TC_GLB(2*N_DEVC+N)
         ENDIF
         IF (OP_INDEX>1 .AND.  (DV%SPATIAL_STATISTIC=='MIN'.OR.DV%SPATIAL_STATISTIC=='MAX')) THEN
            DV%VALUE_1 = TC2_GLB(1,N)
         ENDIF
         IF (OP_INDEX>1 .AND. (DV%SPATIAL_STATISTIC(1:6)=='MINLOC'.OR.DV%SPATIAL_STATISTIC(1:6)=='MAXLOC')) THEN
            NO_NEED_TO_RECV = .FALSE.
            DO NN=1,DV%N_SUBDEVICES
               SDV => DV%SUBDEVICE(NN)
               IF (PROCESS(SDV%MESH)==MY_RANK) THEN
                  IF (SDV%MESH==NINT(TC2_GLB(2,N))) THEN
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

! If a door,entr,exit changes state, save the Smokeview file strings and time of state change

EVAC_ONLY: IF (EVACUATION .AND. MY_RANK==MAX(0,EVAC_PROCESS)) THEN
   I=0  ! Counter for evacuation devices, doors+exits+entrys (evss do not change states)
   DO N=1,N_DOORS
      NM = EVAC_DOORS(N)%IMESH
      IF (.NOT.EVAC_DOORS(N)%SHOW .OR. EMESH_INDEX(NM)==0) CYCLE
      I=I+1
      IF (EVAC_DOORS(N)%IMODE>0 .AND. EVAC_DOORS(N)%IMESH==NM) THEN
         EVAC_DOORS(N)%IMODE=-EVAC_DOORS(N)%IMODE   ! +: change status, -: has already changed status
         M=>MESHES(NM)
         IF (M%N_STRINGS+2>M%N_STRINGS_MAX) CALL RE_ALLOCATE_STRINGS(NM)
         I_STATE=ABS(EVAC_DOORS(N)%IMODE)-1
         M%N_STRINGS = M%N_STRINGS + 1
         WRITE(M%STRING(M%N_STRINGS),'(A,5X,A,1X)') 'DEVICE_ACT',TRIM(EVAC_DOORS(N)%ID)
         M%N_STRINGS = M%N_STRINGS + 1
         WRITE(M%STRING(M%N_STRINGS),'(I6,F10.2,I6)') I+N_DEVC,T,I_STATE
      ENDIF
   ENDDO
   DO N=1,N_EXITS
      NM = EVAC_EXITS(N)%IMESH
      IF (EVAC_EXITS(N)%COUNT_ONLY .OR. .NOT.EVAC_EXITS(N)%SHOW .OR. EMESH_INDEX(NM)==0) CYCLE
      I=I+1
      IF (EVAC_EXITS(N)%IMODE>0 .AND. EVAC_EXITS(N)%IMESH==NM) THEN
         EVAC_EXITS(N)%IMODE=-EVAC_EXITS(N)%IMODE   ! +: change status, -: has already changed status
         M=>MESHES(NM)
         IF (M%N_STRINGS+2>M%N_STRINGS_MAX) CALL RE_ALLOCATE_STRINGS(NM)
         I_STATE=ABS(EVAC_EXITS(N)%IMODE)-1
         M%N_STRINGS = M%N_STRINGS + 1
         WRITE(M%STRING(M%N_STRINGS),'(A,5X,A,1X)') 'DEVICE_ACT',TRIM(EVAC_EXITS(N)%ID)
         M%N_STRINGS = M%N_STRINGS + 1
         WRITE(M%STRING(M%N_STRINGS),'(I6,F10.2,I6)') I+N_DEVC,T,I_STATE
      ENDIF
   ENDDO
   DO N=1,N_ENTRYS
      NM = EVAC_ENTRYS(N)%IMESH
      IF (.NOT.EVAC_ENTRYS(N)%SHOW .OR. EMESH_INDEX(NM)==0) CYCLE
      I=I+1
      IF (EVAC_ENTRYS(N)%IMODE>0 .AND. EVAC_ENTRYS(N)%IMESH==NM) THEN
         EVAC_ENTRYS(N)%IMODE=-EVAC_ENTRYS(N)%IMODE   ! +: change status, -: has already changed status
         M=>MESHES(NM)
         IF (M%N_STRINGS+2>M%N_STRINGS_MAX) CALL RE_ALLOCATE_STRINGS(NM)
         I_STATE=ABS(EVAC_ENTRYS(N)%IMODE)-1
         M%N_STRINGS = M%N_STRINGS + 1
         WRITE(M%STRING(M%N_STRINGS),'(A,5X,A,1X)') 'DEVICE_ACT',TRIM(EVAC_ENTRYS(N)%ID)
         M%N_STRINGS = M%N_STRINGS + 1
         WRITE(M%STRING(M%N_STRINGS),'(I6,F10.2,I6)') I+N_DEVC,T,I_STATE
      ENDIF
   ENDDO
ENDIF EVAC_ONLY

T_USED(7) = T_USED(7) + CURRENT_TIME() - TNOW
END SUBROUTINE EXCHANGE_GLOBAL_OUTPUTS


!> \brief Gather heat release and mass generation rates to root MPI process

SUBROUTINE GATHER_Q_AND_M

REAL_BUFFER_QM(1:N_Q_DOT,1:NMESHES)                           = Q_DOT_SUM(1:N_Q_DOT,1:NMESHES)
REAL_BUFFER_QM(N_Q_DOT+1:N_Q_DOT+N_TRACKED_SPECIES,1:NMESHES) = M_DOT_SUM(1:N_TRACKED_SPECIES,1:NMESHES)
IF (MY_RANK>0) THEN
   CALL MPI_GATHERV(REAL_BUFFER_QM(1,DISPLS(MY_RANK)+1),COUNTS_QM(MY_RANK),MPI_DOUBLE_PRECISION, &
                    REAL_BUFFER_QM,COUNTS_QM,DISPLS_QM,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
ELSE
   CALL MPI_GATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,REAL_BUFFER_QM,COUNTS_QM,DISPLS_QM,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
ENDIF
IF (MY_RANK==0) THEN
   Q_DOT_SUM(1:N_Q_DOT,1:NMESHES)           = REAL_BUFFER_QM(1:N_Q_DOT,1:NMESHES)
   M_DOT_SUM(1:N_TRACKED_SPECIES,1:NMESHES) = REAL_BUFFER_QM(N_Q_DOT+1:N_Q_DOT+N_TRACKED_SPECIES,1:NMESHES)
ENDIF

END SUBROUTINE GATHER_Q_AND_M


SUBROUTINE DUMP_GLOBAL_OUTPUTS

! Dump HRR data to CHID_hrr.csv, MASS data to CHID_mass.csv, DEVICE data to _devc.csv

REAL(EB) :: TNOW
TYPE(DEVICE_TYPE), POINTER :: DV

TNOW = CURRENT_TIME()

IF (EVACUATION .AND. (ICYC<1 .AND. T>T_BEGIN)) RETURN ! No dumps at the evacuation initialization phase

! Dump out HRR info into CHID_hrr.csv

IF (T>=HRR_CLOCK) THEN
   IF (MY_RANK==0) CALL DUMP_HRR(T,DT)
   HRR_CLOCK = HRR_CLOCK + DT_HRR
   Q_DOT_SUM = 0._EB
   M_DOT_SUM = 0._EB
   T_LAST_DUMP_HRR = T
ENDIF

! Dump unstructured geometry and boundary element info

IF (N_FACE>0 .AND. T>=GEOM_CLOCK) THEN
   IF (MY_RANK==0) THEN
      CALL DUMP_GEOM(T,DO_CFACES=.FALSE.)
   ENDIF
   IF (ABS(T-T_BEGIN)<TWO_EPSILON_EB) CALL DUMP_GEOM(T,DO_CFACES=.TRUE.)
   GEOM_CLOCK = GEOM_CLOCK + DT_GEOM
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This block is deprecated, but needs to be removed just prior to FDS 7 release
! to avoid RESTART issues in v6
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (N_GEOM>0 .AND. T>=BNDC_CLOCK) THEN
   !IF (MY_RANK==0) CALL DUMP_BNDC(T)
   BNDC_CLOCK = BNDC_CLOCK + DT_BNDC
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Dump out Evac info

IF (MY_RANK==MAX(0,EVAC_PROCESS)) CALL EVAC_CSV(T)

! Dump out mass info into CHID_mass.csv

IF (T>=MASS_CLOCK) THEN
   IF (MY_RANK==0) CALL DUMP_MASS(T,DT)
   MASS_CLOCK = MASS_CLOCK + DT_MASS
   MASS_DT   = 0._EB
   T_LAST_DUMP_MASS = T
ENDIF

! Dump device info into CHID_devc.csv

IF (T>=DEVC_CLOCK .AND. N_DEVC>0) THEN

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
      DEVC_CLOCK = MIN(DEVC_CLOCK + DT_DEVC, T_END)
      DEVICE_LOOP: DO N=1,N_DEVC
         DV => DEVICE(N)
         IF (T>DV%STATISTICS_END) CYCLE
         IF (DV%NO_UPDATE_DEVC_INDEX>0) THEN
            IF (DEVICE(DV%NO_UPDATE_DEVC_INDEX)%CURRENT_STATE) CYCLE DEVICE_LOOP
         ELSEIF (DV%NO_UPDATE_CTRL_INDEX>0) THEN
            IF (CONTROL(DV%NO_UPDATE_CTRL_INDEX)%CURRENT_STATE) CYCLE DEVICE_LOOP
         ENDIF
         DV%VALUE = 0._EB
         DV%TIME_INTERVAL = 0._EB
      ENDDO DEVICE_LOOP
   ENDIF

ENDIF

! Dump CONTROL info. No gathering required as CONTROL is updated on all meshes

IF (T>=CTRL_CLOCK .AND. N_CTRL>0) THEN
   IF (MY_RANK==0) CALL DUMP_CONTROLS(T)
   CTRL_CLOCK = CTRL_CLOCK + DT_CTRL
ENDIF

! Dump CPU time

IF (T>=CPU_CLOCK) THEN
   CALL DUMP_TIMERS
   CPU_CLOCK = CPU_CLOCK + DT_CPU
ENDIF

T_USED(7) = T_USED(7) + CURRENT_TIME() - TNOW
END SUBROUTINE DUMP_GLOBAL_OUTPUTS


SUBROUTINE INITIALIZE_EVAC

! Initialize evacuation meshes

DO NM=1,NMESHES

   IF (EVACUATION_WRITE_FED .AND. N_MPI_PROCESSES>1 .AND. MY_RANK==EVAC_PROCESS .AND. .NOT.EVACUATION) THEN
      ! Just fire meshes present, but write FED file for evacuation => writing process needs info from others
      IF(PROCESS(NM)/=EVAC_PROCESS) THEN
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
   ENDIF

   IF (PROCESS(NM)/=MY_RANK) CYCLE
   IF (EVACUATION) THEN
      IF (EMESH_INDEX(NM)>0) PART_CLOCK(NM) = T_EVAC + DT_PART
   ENDIF
   IF (MY_RANK/=MAX(0,EVAC_PROCESS)) CYCLE
   IF (EVACUATION) THEN
      CALL INITIALIZE_EVACUATION(NM)
      IF (EMESH_INDEX(NM)>0) CALL DUMP_EVAC(T_EVAC,NM)
   ENDIF
ENDDO
IF (EVACUATION .AND. .NOT.RESTART) ICYC = -EVAC_TIME_ITERATIONS
DT_EVAC=DT

END SUBROUTINE INITIALIZE_EVAC

SUBROUTINE INIT_EVAC_DUMPS

! Initialize evacuation dumps

REAL(EB) :: T_TMP

IF (.NOT.EVACUATION) THEN
   ! Only fire meshes (e.g. phase 2 of the new evacuation scheme)
   IF(.NOT.EVACUATION_WRITE_FED) RETURN ! No evacuation FED file write
   ! Write the FED and smoke information for the evacuation calculation
   T_EVAC      = T_BEGIN
   T_EVAC_SAVE = T_BEGIN
   T_TMP       = T_BEGIN
   IF (MY_RANK==EVAC_PROCESS) CALL INITIALIZE_EVAC_DUMPS
   RETURN
ENDIF

IF (RESTART) THEN
   T_TMP = T
   T_EVAC_SAVE = T_TMP
ELSE
   T_EVAC  = - EVAC_DT_FLOWFIELD*EVAC_TIME_ITERATIONS + T_BEGIN
   T_EVAC_SAVE = T_EVAC
   T_TMP = T_EVAC
END IF
IF (.NOT.EVACUATION) RETURN ! No main evacuation meshes
IF (N_MPI_PROCESSES==1 .OR. (N_MPI_PROCESSES>1 .AND. MY_RANK==EVAC_PROCESS))  CALL INITIALIZE_EVAC_DUMPS

END SUBROUTINE INIT_EVAC_DUMPS


SUBROUTINE EVAC_CSV(T)

! Dump out Evac info

REAL(EB), INTENT(IN) :: T

IF (T>=EVAC_CLOCK .AND. EVACUATION) THEN
   CALL DUMP_EVAC_CSV(T)
   EVAC_CLOCK = EVAC_CLOCK + DT_HRR
ENDIF

END SUBROUTINE EVAC_CSV


SUBROUTINE EVAC_EXCHANGE

LOGICAL EXCHANGE_EVACUATION
INTEGER NM, II, IVENT, I, J, EMESH, JJ, N_END

! Fire mesh information ==> Evac meshes

IF (.NOT. EVACUATION .AND. .NOT.EVACUATION_WRITE_FED) RETURN
IF (N_MPI_PROCESSES>1 .AND. MY_RANK /= EVAC_PROCESS) CALL EVAC_MESH_EXCHANGE(T_EVAC,T_EVAC_SAVE,I_EVAC,ICYC,EXCHANGE_EVACUATION,1)
IF (N_MPI_PROCESSES==1 .OR. (N_MPI_PROCESSES>1 .AND. MY_RANK==EVAC_PROCESS)) &
     CALL EVAC_MESH_EXCHANGE(T_EVAC,T_EVAC_SAVE,I_EVAC,ICYC,EXCHANGE_EVACUATION,2)

! Update evacuation devices

DO NM=1,NMESHES
   IF (.NOT.EVACUATION) CYCLE
   IF (EMESH_INDEX(NM)==0.OR.EVACUATION_SKIP(NM)) CYCLE
   IF (MY_RANK/=MAX(0,EVAC_PROCESS)) CYCLE
   CALL UPDATE_GLOBAL_OUTPUTS(T,DT,NM)
ENDDO

! Save the evacuation flow fields to the arrays U_EVAC and V_EVAC

N_END = N_EXITS - N_CO_EXITS + N_DOORS
DO NM = 1, NMESHES
   IF (.NOT.EVACUATION) CYCLE
   IF (EMESH_INDEX(NM)==0.OR.EVACUATION_SKIP(NM)) CYCLE
   IF (MY_RANK /= MAX(0,EVAC_PROCESS)) CYCLE
   II = EVAC_TIME_ITERATIONS / MAXVAL(EMESH_NFIELDS)
   IF (MOD(ABS(ICYC),II)==0) THEN
      IVENT = (ABS(ICYC))/II + 1
      LOOP_EXITS: DO JJ = 1, N_END
         IF (EMESH_EXITS(JJ)%MAINMESH == NM .AND. EMESH_EXITS(JJ)%I_DOORS_EMESH == IVENT) THEN
            EMESH = EMESH_EXITS(JJ)%EMESH
            DO J = 0, EMESH_IJK(2,EMESH) + 1
               DO I = 0, EMESH_IJK(1,EMESH) + 1
                  ! FB is the EFF file precision. Convert to FB here so that EFF read calculation gives
                  ! exactly the same results as EFF write calculation for the agents.
                  IF (MESHES(NM)%PRESSURE_ZONE(I,J,1)>0) THEN
                     EMESH_EXITS(JJ)%U_EVAC(I,J) = REAL(MESHES(NM)%U(I,J,1),FB)
                     EMESH_EXITS(JJ)%V_EVAC(I,J) = REAL(MESHES(NM)%V(I,J,1),FB)
                  ELSE
                     EMESH_EXITS(JJ)%U_EVAC(I,J) = 0.0_FB
                     EMESH_EXITS(JJ)%V_EVAC(I,J) = 0.0_FB
                  END IF
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

COMPUTE_PRESSURE_LOOP: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   IF (EVACUATION_SKIP(NM).OR..NOT.EVACUATION) CYCLE COMPUTE_PRESSURE_LOOP
   PRESSURE_ITERATION_LOOP: DO N=1,EVAC_PRESSURE_ITERATIONS
      CALL NO_FLUX(DT,NM)
      MESHES(NM)%FVZ = 0._EB
      CALL PRESSURE_SOLVER_COMPUTE_RHS(T,DT,NM)
      CALL PRESSURE_SOLVER_FFT(NM)
      CALL PRESSURE_SOLVER_CHECK_RESIDUALS(NM)
   ENDDO PRESSURE_ITERATION_LOOP
ENDDO COMPUTE_PRESSURE_LOOP

END SUBROUTINE EVAC_PRESSURE_ITERATION_SCHEME


SUBROUTINE EVAC_MAIN_LOOP
USE RADCONS, ONLY: TIME_STEP_INCREMENT

! Call the evacuation routine and adjust the time steps for the evacuation meshes

REAL(EB) :: T_FIRE, EVAC_DT, DT_TMP
INTEGER :: II

EVACUATION_SKIP=.FALSE. ! Do not skip the flow calculation

IF (EVACUATION_WRITE_FED.AND. .NOT.EVACUATION) THEN
   ! New fire+evacuation calculation strategy: Phase 2 just fire meshes present
   ! ToDo: Set the mesh exchange options here, if needed (this is just a placeholder right now)
   T_EVAC = T
   EVAC_DT = DT
   RETURN
ENDIF

EVAC_DT = EVAC_DT_STEADY_STATE
T_FIRE = T_EVAC + EVAC_DT
IF (ICYC < 1) EVAC_DT = EVAC_DT_FLOWFIELD
IF (ICYC < 1) DT = EVAC_DT
IF (ICYC < 1) DT_NEW = DT_EVAC
IF (ICYC == 1) DT = DT_EVAC ! Initial fire dt that was read in
IF (ICYC == 1) T  = T_BEGIN ! Initial fire t  that was read in
IF (ICYC == 1) T_EVAC = T_BEGIN - 0.1_EB*MIN(EVAC_DT_FLOWFIELD,EVAC_DT_STEADY_STATE)
IF (ICYC > 0) T_FIRE  = T

DT_TMP = DT
IF ((T+DT)>=T_END) DT_TMP = MAX(MIN(EVAC_DT_STEADY_STATE,T_END-T_EVAC),1.E-10_EB)
IF (ICYC > 0) EVAC_DT = DT_TMP

DO NM = 1, NMESHES
   IF (EVACUATION.AND.EMESH_INDEX(NM)==0) EVACUATION_SKIP(NM) = .TRUE.
   IF (EVACUATION.AND.EMESH_INDEX(NM)>0) DT_NEW(NM) = EVAC_DT
END DO
IF (ICYC <= 0) THEN
   DO NM = 1, NMESHES
      IF (.NOT.EVACUATION) EVACUATION_SKIP(NM) = .TRUE.  ! Be sure that no fire meshes are updated for icyc < 0
   END DO
ENDIF
IF (.NOT.EVACUATION .AND. RADIATION .AND. ICYC > 0) THEN
   DO NM = 1, NMESHES
      IF (.NOT.EVACUATION) CYCLE
      MESHES(NM)%RAD_CALL_COUNTER  = MESHES(NM)%RAD_CALL_COUNTER + 1
      IF (MOD(MESHES(NM)%RAD_CALL_COUNTER,TIME_STEP_INCREMENT)==0 .OR. ICYC==1) THEN
         EXCHANGE_RADIATION = .TRUE.
      ELSE
         EXCHANGE_RADIATION = .FALSE.
      ENDIF
   ENDDO
ENDIF

EVAC_TIME_STEP_LOOP: DO WHILE (T_EVAC < T_FIRE)
   T_EVAC = T_EVAC + EVAC_DT
   IF (N_MPI_PROCESSES==1 .OR. (N_MPI_PROCESSES>1 .AND. MY_RANK==EVAC_PROCESS)) CALL PREPARE_TO_EVACUATE(ICYC)
   DO NM = 1, NMESHES
      IF (EVACUATION) THEN
         EVACUATION_SKIP(NM)  = .TRUE.
         IF (ICYC <= 1 .AND. .NOT.BTEST(I_EVAC, 2)) THEN
            IF (ICYC <= 0 .AND. EMESH_INDEX(NM)>0) THEN
               II = EVAC_TIME_ITERATIONS / MAXVAL(EMESH_NFIELDS)
               IF ((ABS(ICYC)+1) <= EMESH_NFIELDS(EMESH_INDEX(NM))*II) THEN
                  EVACUATION_SKIP(NM) = .FALSE.
               ELSE
                  EVACUATION_SKIP(NM) = .TRUE.
               END IF
               DIAGNOSTICS = .FALSE.
               IF (.NOT.EVACUATION) EVACUATION_SKIP(NM) = .TRUE.
            END IF
            !
            IF (ICYC <= 0) T = T_EVAC + EVAC_DT_FLOWFIELD*EVAC_TIME_ITERATIONS - EVAC_DT_FLOWFIELD
         ENDIF
         IF (ICYC <= 1 .AND. BTEST(I_EVAC, 2)) THEN
            IF (ICYC <= 0 .AND. EMESH_INDEX(NM)>0) THEN
               EVACUATION_SKIP(NM) = .TRUE.
               DIAGNOSTICS = .FALSE.
            END IF
            IF (ICYC <= 0) T = T_EVAC + EVAC_DT_FLOWFIELD*EVAC_TIME_ITERATIONS - EVAC_DT_FLOWFIELD
         ENDIF
         IF (EMESH_INDEX(NM)==0) THEN
            VELOCITY_ERROR_MAX_LOC(:,NM) = 1
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
         IF (EMESH_INDEX(NM)>0) THEN
            IF (PROCESS(NM)==MY_RANK .AND. STOP_STATUS==NO_STOP) CALL EVACUATE_HUMANS(T_EVAC,DT_TMP,NM,ICYC)
            IF (T_EVAC >= PART_CLOCK(NM)) THEN
               IF (PROCESS(NM)==MY_RANK) CALL DUMP_EVAC(T_EVAC, NM)
               DO
                  PART_CLOCK(NM) = PART_CLOCK(NM) + DT_PART
                  IF (PART_CLOCK(NM) >= T_EVAC) EXIT
               ENDDO
            ENDIF
         ENDIF
      ENDIF
   ENDDO
   IF (ICYC < 1) EXIT EVAC_TIME_STEP_LOOP
   IF (N_MPI_PROCESSES==1 .OR. (N_MPI_PROCESSES>1 .AND. MY_RANK==EVAC_PROCESS)) CALL CLEAN_AFTER_EVACUATE(ICYC, I_EVAC)
ENDDO EVAC_TIME_STEP_LOOP
IF (ICYC < 1 .AND. MY_RANK==0) THEN
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
INTEGER :: NN
REAL(EB) :: TNOW

TNOW = CURRENT_TIME()

! Pack HVAC values into REAL_BUFFER_HVAC

REAL_BUFFER_HVAC(              1:  N_DUCTNODES,:) = NODE_H(1:N_DUCTNODES,:)
REAL_BUFFER_HVAC(  N_DUCTNODES+1:2*N_DUCTNODES,:) = NODE_P(1:N_DUCTNODES,:)
REAL_BUFFER_HVAC(2*N_DUCTNODES+1:3*N_DUCTNODES,:) = NODE_RHO(1:N_DUCTNODES,:)
REAL_BUFFER_HVAC(3*N_DUCTNODES+1:4*N_DUCTNODES,:) = NODE_TMP(1:N_DUCTNODES,:)
REAL_BUFFER_HVAC(4*N_DUCTNODES+1:5*N_DUCTNODES,:) = NODE_X(1:N_DUCTNODES,:)
REAL_BUFFER_HVAC(5*N_DUCTNODES+1:6*N_DUCTNODES,:) = NODE_Y(1:N_DUCTNODES,:)
REAL_BUFFER_HVAC(6*N_DUCTNODES+1:7*N_DUCTNODES,:) = NODE_Z(1:N_DUCTNODES,:)
REAL_BUFFER_HVAC(7*N_DUCTNODES+1:8*N_DUCTNODES,:) = NODE_AREA(1:N_DUCTNODES,:)
REAL_BUFFER_HVAC(8*N_DUCTNODES+1:9*N_DUCTNODES,:) = NODE_ZONE(1:N_DUCTNODES,:)
NN = 8
DO N=1,N_TRACKED_SPECIES
   NN = NN + 1
   REAL_BUFFER_HVAC(NN*N_DUCTNODES+1:(NN+1)*N_DUCTNODES,:) = NODE_ZZ(1:N_DUCTNODES,N,:)
ENDDO

IF (MY_RANK>0) THEN
   CALL MPI_GATHERV(REAL_BUFFER_HVAC(1,DISPLS(MY_RANK)+1),COUNTS_HVAC(MY_RANK),MPI_DOUBLE_PRECISION, &
                    REAL_BUFFER_HVAC,COUNTS_HVAC,DISPLS_HVAC,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
ELSE
   CALL MPI_GATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL, &
                    REAL_BUFFER_HVAC,COUNTS_HVAC,DISPLS_HVAC,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
ENDIF

! Unpack HVAC values on MPI process 0

IF (MY_RANK==0) THEN
   NODE_H(1:N_DUCTNODES,:)    =     REAL_BUFFER_HVAC(              1:  N_DUCTNODES,:)
   NODE_P(1:N_DUCTNODES,:)    =     REAL_BUFFER_HVAC(  N_DUCTNODES+1:2*N_DUCTNODES,:)
   NODE_RHO(1:N_DUCTNODES,:)  =     REAL_BUFFER_HVAC(2*N_DUCTNODES+1:3*N_DUCTNODES,:)
   NODE_TMP(1:N_DUCTNODES,:)  =     REAL_BUFFER_HVAC(3*N_DUCTNODES+1:4*N_DUCTNODES,:)
   NODE_X(1:N_DUCTNODES,:)    =     REAL_BUFFER_HVAC(4*N_DUCTNODES+1:5*N_DUCTNODES,:)
   NODE_Y(1:N_DUCTNODES,:)    =     REAL_BUFFER_HVAC(5*N_DUCTNODES+1:6*N_DUCTNODES,:)
   NODE_Z(1:N_DUCTNODES,:)    =     REAL_BUFFER_HVAC(6*N_DUCTNODES+1:7*N_DUCTNODES,:)
   NODE_AREA(1:N_DUCTNODES,:) =     REAL_BUFFER_HVAC(7*N_DUCTNODES+1:8*N_DUCTNODES,:)
   NODE_ZONE(1:N_DUCTNODES,:) = INT(REAL_BUFFER_HVAC(8*N_DUCTNODES+1:9*N_DUCTNODES,:))
   NN = 8
   DO N=1,N_TRACKED_SPECIES
      NN = NN + 1
      NODE_ZZ(1:N_DUCTNODES,N,:) = REAL_BUFFER_HVAC(NN*N_DUCTNODES+1:(NN+1)*N_DUCTNODES,:)
   ENDDO
ENDIF

T_USED(11)=T_USED(11) + CURRENT_TIME() - TNOW
END SUBROUTINE EXCHANGE_HVAC_BC


SUBROUTINE EXCHANGE_HVAC_SOLUTION

! Exchange information mesh to mesh needed for performing the HVAC computation

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

CALL MPI_BCAST(REAL_BUFFER_DUCT,(2+N_TRACKED_SPECIES)*N_DUCTNODES+N_DUCTS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

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


!> \brief Check to see if any FREEZE_VELOCITY=T and any SOLID_PHASE_ONLY=T

SUBROUTINE CHECK_FREEZE_VELOCITY_STATUS

CALL MPI_ALLREDUCE(MPI_IN_PLACE,FREEZE_VELOCITY ,INTEGER_ONE,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,IERR)
CALL MPI_ALLREDUCE(MPI_IN_PLACE,SOLID_PHASE_ONLY,INTEGER_ONE,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,IERR)
IF (FREEZE_VELOCITY) CHECK_FREEZE_VELOCITY = .FALSE.

END SUBROUTINE CHECK_FREEZE_VELOCITY_STATUS


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
