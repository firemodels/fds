!> \brief Compute the HVAC mass and energy transport
!> \details Module contains routines to read the HVAC namelist inputs, intialize the HVAC solver, and solve
!> the flow for each timestep. Note that the HVAC solver is only called for the first MPI process. This requires
!> that HVAC boundary conditions at VENTs be aggregated over MPI processes and that the HVAC solution
!> be shared with all MPI processes.

MODULE HVAC_ROUTINES

USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_POINTERS
USE DEVICE_VARIABLES
USE CONTROL_VARIABLES
USE COMP_FUNCTIONS, ONLY: CURRENT_TIME, CHECKREAD, SHUTDOWN
USE MEMORY_FUNCTIONS, ONLY: ChkMemErr

IMPLICIT NONE (TYPE,EXTERNAL)

PRIVATE

REAL(EB) :: DT_HV !< Size of subtimestep used in HVAC solver
REAL(EB) :: DT_QF = 1._EB !< Size of subtimestep used in HVAC QFAN
REAL(EB) :: DT_MT !< Size of subtimestep used in 1D mass+energy transport solver
CHARACTER(LABEL_LENGTH), ALLOCATABLE, DIMENSION(:,:) :: NODE_DUCT_A  !< Temporary array storing DUCT_ID inputs for NODEs
CHARACTER(LABEL_LENGTH), ALLOCATABLE, DIMENSION(:,:) :: DUCT_NODE_A  !< Temporary array storing NODE_ID inputs for DUCTs
CHARACTER(LABEL_LENGTH), ALLOCATABLE, DIMENSION(:) :: NODE_FILTER_A  !< Temporary array storing FILTER_ID inputs for NODEs
CHARACTER(LABEL_LENGTH), ALLOCATABLE, DIMENSION(:) :: DUCT_FAN_A     !< Temporary array storing FAN_ID inputs for DUCTs
CHARACTER(LABEL_LENGTH), ALLOCATABLE, DIMENSION(:) :: DUCT_AIRCOIL_A !< Temporary array storing AIRCOIL_ID inputs for DUCTS
INTEGER :: LEAK_DUCTS = 0 !< Number of ducts used for leakage
INTEGER, ALLOCATABLE, DIMENSION(:,:):: LEAK_PATH !< Temporary array used to determine ducts to create for leakage paths
CHARACTER(255) :: MESSAGE !< Stores ERROR or WARNING message written to LU_ERR
REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: LHS !< Left hand side of HVAC solution array
REAL(EB), ALLOCATABLE, DIMENSION(:) :: RHS !< Right hand side of HVAC solution array
REAL(EB), ALLOCATABLE, DIMENSION(:) :: DPSTAR !< Array of extrapoloated presure for each ZONE.
INTEGER :: ITER !< Current HVAC solver iterations
INTEGER :: ITER_MAX=10  !< Maximum allowed solver iterations
INTEGER :: ITER_MAX_QFAN=30 !< Maximum allowed solver iterations for QFAN
LOGICAL :: DUCT_NODE_INIT  !< Flag indicating duct nodes have been initialized
LOGICAL :: TRANSPORT_PARTICLES !< Flag indicating particles should be transferred across a localized leak path

PUBLIC HVAC_CALC,READ_HVAC,PROC_HVAC,HVAC_BC_IN,FIND_NETWORKS,COLLAPSE_HVAC_BC,SET_INIT_HVAC,INIT_DUCT_NODE,LEAKAGE_HVAC

REAL(EB),PUBLIC, ALLOCATABLE, DIMENSION(:):: NODE_AREA_EX !< Contain sum of area over all MESHES of all VENTs assigned to each NODE
REAL(EB),PUBLIC, ALLOCATABLE, DIMENSION(:):: NODE_TMP_EX
!< Contains sum of area weighted temperature over all MESHES for all VENTs assigned to each NODE
REAL(EB),PUBLIC, ALLOCATABLE, DIMENSION(:):: DUCT_MF !<Contains mass flow for each duct. Exchanged during MPI exchanges

REAL(EB), PUBLIC, TARGET, ALLOCATABLE, DIMENSION(:,:)  :: NODE_PROPERTIES !< Holding array that is used to do MPI exchanges
REAL(EB), POINTER, DIMENSION(:)   :: NODE_AREA !< Area of each NODE per MESH
REAL(EB), POINTER, DIMENSION(:)   :: NODE_H    !< Area-weighted enthalpy of each NODE
REAL(EB), POINTER, DIMENSION(:)   :: NODE_P    !< Area-weighed pressure of each NODE
REAL(EB), POINTER, DIMENSION(:)   :: NODE_RHO  !< Area-weighted desnity of each NODE
REAL(EB), POINTER, DIMENSION(:)   :: NODE_X    !< Area-weighted x position of each NODE
REAL(EB), POINTER, DIMENSION(:)   :: NODE_Y    !< Area-weighted y position of each NODE
REAL(EB), POINTER, DIMENSION(:)   :: NODE_Z    !< Area-weighted z position of each NODE
REAL(EB), POINTER, DIMENSION(:)   :: NODE_TMP  !< Area-weighted temperature of each NODE
REAL(EB), POINTER, DIMENSION(:,:) :: NODE_ZZ   !< Area-weighted tracked species mass fractions of each NODE

INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: NODE_ZONE !< Array of NODEs belonging to each ZONE

REAL(EB),PUBLIC, ALLOCATABLE, DIMENSION(:,:):: NODE_ZZ_EX !<Area weighted tracked species mass fractions of each NODE per MESH
!< Contains sum of area weighted tracked species mass fractions over all MPI processes for all VENTs assigned to each NODE
REAL(EB),PUBLIC, ALLOCATABLE, DIMENSION(:):: PSUM_TOT !< Contains sum of PSUM for merged pressure zones
INTEGER, PUBLIC :: N_DUCT_QUANTITY=0 !< Number of DUCT output QUANTITY
INTEGER, PUBLIC :: N_NODE_QUANTITY=0 !< Number of NODE output QUANTITY

TYPE(HVAC_QUANTITY_TYPE), PUBLIC, TARGET, ALLOCATABLE, DIMENSION(:) :: DUCT_QUANTITY_ARRAY !< DUCT outputs for .hvac file
TYPE(HVAC_QUANTITY_TYPE), PUBLIC, TARGET, ALLOCATABLE, DIMENSION(:) :: NODE_QUANTITY_ARRAY !< NODE outputs for .hvac file

CONTAINS

!> \brief Reads and processes the HVAC namelist inputs

SUBROUTINE READ_HVAC

USE MISC_FUNCTIONS, ONLY: SEARCH_CONTROLLER,GET_RAMP_INDEX
INTEGER , PARAMETER :: MAX_DUCTS = 20 !< Maximum number of ducts connected to a node
INTEGER :: IOS !< Used for returning the status of a READ statement
INTEGER :: IZERO !< Used for returning the status of an ALLOCATE statement
INTEGER :: N_HVAC_READ !< Counter for number of HVAC inputs that have been read in
INTEGER :: N,NC,ND,NN,NS
INTEGER :: I_AIRCOIL=0 !< AIRCOIL array index
INTEGER :: I_DUCT=0 !< DUCT array index
INTEGER :: I_DUCTNODE=0 !< DUCTNODE array index
INTEGER :: I_FAN=0 !< FAN array index
INTEGER :: I_FILTER=0 !< FILTER array index
INTEGER :: N_CELLS !< Number of cells in a DUCT with HVAC_MASS_TRANSPORT
REAL(EB) :: AREA !< Area (m2) of a DUCT.
REAL(EB) :: DIAMETER !< Diameter (m) of a DUCT.
REAL(EB) :: XYZ(3) !< Position (m) of a DUCTNODE.
REAL(EB) :: LOSS(MAX_DUCTS,MAX_DUCTS) !< Array of flow losses for a DUCT or DUCTNODE.
REAL(EB) :: VOLUME_FLOW !< Fixed volume flow (m3/s) for a FAN or a DUCT.
REAL(EB) :: MAX_FLOW !< Flow at zero pressure (m3/s) for a quadratic FAN.
REAL(EB) :: MAX_PRESSURE !< Stall pressure (Pa) for a quadratic FAN.
REAL(EB) :: ROUGHNESS !< Absolute roughness (m) of a DUCT.
REAL(EB) :: LENGTH !< Length (m) of a DUCT.
REAL(EB) :: TNOW !< Current CPU time (s) used in computing length of time spent in HVAC routines.
REAL(EB) :: TAU_AC !< Time constant (s) for full effectiveness of an AIRCOIL.
REAL(EB) :: TAU_FAN !< Time constant (s) for a FAN to reach its full flow rate.
REAL(EB) :: TAU_VF !< Time constant (s) for fixed flow in a DUCT.
REAL(EB) :: FIXED_Q !< Fixed heat transfer (input as kW) for an AIRCOIL.
REAL(EB) :: CLEAN_LOSS !< Flow loss for a clean filter (e.g. LOADING=0)
REAL(EB) :: COOLANT_MASS_FLOW !< Mass flow (kg/s)of the working fluid for an AIRCOIL.
REAL(EB) :: COOLANT_SPECIFIC_HEAT !< Specific heat (input as kJ/kg/K) of the working fluid for an AIRCOIL.
REAL(EB) :: COOLANT_TEMPERATURE !< Temperature (input as C) of the working fluid for an AIRCOIL.
REAL(EB) :: PERIMETER !< Perimeter (m) of a duct. Computes the hydraulic DIAMETER when specified with AREA.
REAL(EB) :: MASS_FLOW !< Fixed mass flow (kg/s) in a DUCT.
REAL(EB) :: LOADING(MAX_SPECIES) !< Initial species loading (kg) of a FILTER.
REAL(EB) :: EFFICIENCY(MAX_SPECIES) !< Fraction of a species trapped by a FILTER.
REAL(EB) :: LOADING_MULTIPLIER(MAX_SPECIES) !< Multiplier of LOADING used in computing the loss of a FILTER.
REAL(EB) :: LEAK_PRESSURE_EXPONENT !< Exponent for pressure in leakage equation
REAL(EB) :: LEAK_REFERENCE_PRESSURE !< Reference pressure in leakage equation (Pa)
REAL(EB) :: DISCHARGE_COEFFICIENT !< Discharge coefficient leakage equation
REAL(EB) :: WAYPOINTS(30,3) !< Duct waypoints (m)
LOGICAL :: DEBUG !< Flag indicating known values are output to smokeview for debugging
LOGICAL :: ROUND !< Flag indicating DUCT has a round cross-section
LOGICAL :: SQUARE !< Flag indicating DUCT has a square cross-section
LOGICAL :: DAMPER !< Flag indicating that a damper is present in a DUCT.
LOGICAL :: REVERSE !< Flag indicating that a specfied flow or FAN in a DUCT is from the second to the first node.
LOGICAL :: AMBIENT !< Flag indicating a DUCTNODE is connected to the ambient.
LOGICAL :: GEOM !< Flag indicating a DUCTNODE or LEAKAGE VENT_ID is connected to GEOM.
LOGICAL :: GEOM2 !< Flag indicating LEAKAGE VENT2_ID is connected to GEOM.
LOGICAL :: LEAK_ENTHALPY !< Flag indicating that the boundary condition for a LEAKAGE duct should preserve enthalpy.
LOGICAL :: INITIALIZED_HVAC_MASS_TRANSPORT !< Flag indicating DUCTs with N_CELLS>1 have been initiazed.
LOGICAL :: DUCT_QUANTITY_DEFINED=.FALSE. !< Flag indicating a DUCT_QUANTITY list has alreayd been defined
LOGICAL :: NODE_QUANTITY_DEFINED=.FALSE. !< Flag indicating a NODE_QUANTITY list has alreayd been defined
LOGICAL :: DRY(20) !< Flag indicating species output is DRY
CHARACTER(LABEL_LENGTH) :: AIRCOIL_ID !< ID of an AIRCOIL located in a DUCT.
CHARACTER(LABEL_LENGTH) :: CTRL_ID !< Name of a control function controlling a FAN, damper, or AIRCOIL.
CHARACTER(LABEL_LENGTH) :: DEVC_ID !< Name of a device controlling a FAN, damper, or AIRCOIL.
CHARACTER(LABEL_LENGTH) :: DUCT_ID(MAX_DUCTS) !<IDs of DUCTs connected to a DUCTNODE.
CHARACTER(LABEL_LENGTH) :: FAN_ID !< ID of a FAN located in a DUCT.
CHARACTER(LABEL_LENGTH) :: FILTER_ID !< ID of a FILTER located at a DUCTNODE.
CHARACTER(LABEL_LENGTH) :: ID !< Name of an HVAC component
CHARACTER(LABEL_LENGTH) :: NETWORK_ID !< ID of the network.
CHARACTER(LABEL_LENGTH) :: NODE_ID(2) !< IDs of the nodes for each end of a DUCT.
CHARACTER(LABEL_LENGTH) :: QUANTITY(20) !QUANTITY list for .hvac file.
CHARACTER(LABEL_LENGTH) :: QUANTITY_SPEC_ID(20) !<SPEC_ID for QUANTITY list for .hvac file.
CHARACTER(LABEL_LENGTH) :: RAMP_ID !< Name of a RAMP for DUCT flow, FAN curve, or AIRCOIL heat exchange.
CHARACTER(LABEL_LENGTH) :: RAMP_LOSS !< Name of a RAMP for the flow loss of a variable damper where T=damper position and F=loss.
CHARACTER(LABEL_LENGTH) :: SPEC_ID(MAX_SPECIES) !< List of species that are trapped by a filter.
CHARACTER(LABEL_LENGTH) :: TYPE_ID !< Type of HVAC component (e.g. DUCT, FAN, DUCTNODE, etc.)
CHARACTER(LABEL_LENGTH) :: VENT_ID !< Name of a VENT connected to a DUCTNODE or the first node for a LEAKAGE duct
CHARACTER(LABEL_LENGTH) :: VENT2_ID !< VENT connected to the second node for a LEAKAGE duct
TYPE(DUCTNODE_TYPE), POINTER :: DN !< Pointer to a DUCTNODE
TYPE(DUCT_TYPE), POINTER :: DU !< Pointer to a DUCT
TYPE(HVAC_QUANTITY_TYPE), POINTER :: HQT !< Pointer to a DUCT_ or NODE_QUANTITY_ARRAY
NAMELIST /HVAC/ AIRCOIL_ID,AMBIENT,AREA,CLEAN_LOSS,COOLANT_SPECIFIC_HEAT,COOLANT_MASS_FLOW,COOLANT_TEMPERATURE,CTRL_ID,&
                DAMPER,DEBUG,DEVC_ID,DIAMETER,DISCHARGE_COEFFICIENT,DRY,DUCT_ID,&
                EFFICIENCY,FAN_ID,FILTER_ID,FIXED_Q,GEOM,GEOM2,ID,LEAK_ENTHALPY,LEAK_PRESSURE_EXPONENT,LEAK_REFERENCE_PRESSURE,&
                LENGTH,LOADING,LOADING_MULTIPLIER,LOSS,&
                MASS_FLOW,MAX_FLOW,MAX_PRESSURE,N_CELLS,NETWORK_ID,NODE_ID,PERIMETER,QUANTITY,QUANTITY_SPEC_ID,&
                RAMP_ID,RAMP_LOSS,REVERSE,ROUGHNESS,ROUND,SPEC_ID,SQUARE,TAU_AC,TAU_FAN,TAU_VF,TRANSPORT_PARTICLES,&
                TYPE_ID,WAYPOINTS,VENT_ID,VENT2_ID,VOLUME_FLOW,XYZ

TNOW=CURRENT_TIME()

N_HVAC_READ = 0

DUCT_NODE_INIT = .FALSE.

REWIND(LU_INPUT) ; INPUT_FILE_LINE_NUMBER = 0
COUNT_HVAC_LOOP: DO
   CALL CHECKREAD('HVAC',LU_INPUT,IOS)  ; IF (STOP_STATUS==SETUP_STOP) RETURN
   IF (IOS==1) EXIT COUNT_HVAC_LOOP
   DEBUG = .FALSE.
   READ(LU_INPUT,HVAC,END=15,ERR=16,IOSTAT=IOS)
   IF (DEBUG) HVAC_DEBUG = DEBUG
   N_HVAC_READ = N_HVAC_READ + 1
   16 IF (IOS>0) THEN
         WRITE(MESSAGE,'(A,I5,A,I5)') &
            'ERROR(101): Problem with HVAC line number ',N_HVAC_READ+1,', input line number',INPUT_FILE_LINE_NUMBER
         CALL SHUTDOWN(MESSAGE); RETURN
      ENDIF
   IF (TRIM(ID)=='null') THEN
      WRITE(MESSAGE,'(A,I5,A,I5)') &
         'ERROR(501): No ID provided for HVAC line number ',N_HVAC_READ+1,', input line number',INPUT_FILE_LINE_NUMBER
      CALL SHUTDOWN(MESSAGE); RETURN
   ENDIF
   SELECT CASE(TYPE_ID)
      CASE ('DUCT')
         N_DUCTS = N_DUCTS + 1
      CASE ('NODE')
         N_DUCTNODES = N_DUCTNODES + 1
      CASE ('FILTER')
         N_FILTERS = N_FILTERS + 1
      CASE ('FAN')
         N_FANS = N_FANS + 1
      CASE ('AIRCOIL')
         N_AIRCOILS = N_AIRCOILS + 1
      CASE ('LEAK')
         N_DUCTS = N_DUCTS + 1
         N_DUCTNODES = N_DUCTNODES + 2
      CASE ('DUCT QUANTITY LIST')
      CASE ('NODE QUANTITY LIST')
      CASE DEFAULT
         WRITE(MESSAGE,'(A,I5,A,I5)') &
            'ERROR(502): Invalid TYPE_ID provided for HVAC line number ',N_HVAC_READ,', input line number',INPUT_FILE_LINE_NUMBER
         CALL SHUTDOWN(MESSAGE); RETURN
   END SELECT
ENDDO COUNT_HVAC_LOOP
15 CONTINUE

LEAK_DUCTS = 0
IF (ANY(SURFACE%LEAK_PATH(1)>0).OR.ANY(SURFACE%LEAK_PATH(2)>0)) THEN
   ALLOCATE(LEAK_PATH(0:N_ZONE,0:N_ZONE),STAT=IZERO)
   CALL ChkMemErr('HVAC','LEAK_PATH',IZERO)
   LEAK_PATH = 0
   DO NS=1,N_SURF
      IF (SURFACE(NS)%LEAK_PATH(1)>0 .OR. SURFACE(NS)%LEAK_PATH(2)>0) &
         LEAK_PATH(MINVAL(SURFACE(NS)%LEAK_PATH),MAXVAL(SURFACE(NS)%LEAK_PATH))=1
   ENDDO
   LEAK_DUCTS = SUM(LEAK_PATH)
   N_DUCTS = N_DUCTS + LEAK_DUCTS
   N_DUCTNODES = N_DUCTNODES + 2 * LEAK_DUCTS
ENDIF

IF (N_DUCTS > 0) HVAC_SOLVE = .TRUE.

IF ((N_DUCTS > 0 .AND. N_DUCTNODES <= 0) .OR. (N_DUCTS <= 0 .AND. N_DUCTNODES > 0)) THEN
   WRITE(MESSAGE,'(A)') 'ERROR(503): Must have both DUCTs and NODEs in the input file'
   CALL SHUTDOWN(MESSAGE); RETURN
ENDIF

ALLOCATE(DUCT(N_DUCTS),STAT=IZERO)
CALL ChkMemErr('HVAC','DUCT',IZERO)
ALLOCATE(DUCT_NE(N_DUCTS),STAT=IZERO)
CALL ChkMemErr('HVAC','DUCT_NE',IZERO)
DUCT_NE = 0
ALLOCATE(DUCT_DR(N_DUCTS),STAT=IZERO)
CALL ChkMemErr('HVAC','DUCT_DR',IZERO)
DUCT_DR = 0
ALLOCATE(DUCTNODE(N_DUCTNODES),STAT=IZERO)
CALL ChkMemErr('HVAC','DUCTNODE',IZERO)
ALLOCATE(DUCTNODE_NE(N_DUCTNODES),STAT=IZERO)
CALL ChkMemErr('HVAC','DUCTNODE_NE',IZERO)
DUCTNODE_NE = 0
ALLOCATE(DUCTNODE_DR(N_DUCTNODES),STAT=IZERO)
CALL ChkMemErr('HVAC','DUCTNODE_DR',IZERO)
DUCTNODE_DR = 0
ALLOCATE(FILTER(N_FILTERS),STAT=IZERO)
CALL ChkMemErr('HVAC','FILTER',IZERO)
ALLOCATE(FAN(N_FANS),STAT=IZERO)
CALL ChkMemErr('HVAC','FAN',IZERO)
ALLOCATE(AIRCOIL(N_AIRCOILS),STAT=IZERO)
CALL ChkMemErr('HVAC','AIRCOIL',IZERO)

!Temp arrays for input processing
ALLOCATE(DUCT_NODE_A(N_DUCTS,2),STAT=IZERO)
CALL ChkMemErr('HVAC','NODE_DUCTS',IZERO)
ALLOCATE(NODE_DUCT_A(N_DUCTNODES,MAX_DUCTS),STAT=IZERO)
CALL ChkMemErr('HVAC','NODE_DUCTS',IZERO)
ALLOCATE(NODE_FILTER_A(N_DUCTNODES),STAT=IZERO)
CALL ChkMemErr('HVAC','NODE_DUCTS',IZERO)
ALLOCATE(DUCT_FAN_A(N_DUCTS),STAT=IZERO)
CALL ChkMemErr('HVAC','DUCT_FAN',IZERO)
ALLOCATE(DUCT_AIRCOIL_A(N_DUCTS),STAT=IZERO)
CALL ChkMemErr('HVAC','DUCT_AIRCOIL',IZERO)

REWIND(LU_INPUT) ; INPUT_FILE_LINE_NUMBER = 0
DO NN=1,N_HVAC_READ
   CALL SET_HVAC_DEFAULTS
   READ(LU_INPUT,HVAC)
   SELECT CASE (TYPE_ID)
      CASE('DUCT')
         I_DUCT = I_DUCT + 1
         DU=> DUCT(I_DUCT)
         DU%ID   = ID
         IF (TRIM(ID)=='null' ) THEN
            WRITE(MESSAGE,'(A,I5)') 'ERROR(504): Duct has no ID, HVAC line number:',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         IF (DIAMETER < 0._EB .AND. AREA < 0._EB .AND. PERIMETER < 0._EB) THEN
            WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR(505): Duct has no AREA, DIAMETER, or PERIMTER. Duct ID:',TRIM(ID),&
                                        ', HVAC line number:',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         IF (AREA < 0._EB) THEN
            IF (DIAMETER < 0._EB) THEN
               WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR(506): Duct without AREA has no DIAMETER. Duct ID:',TRIM(ID),&
                                           ', HVAC line number:',NN
               CALL SHUTDOWN(MESSAGE); RETURN
            ENDIF
            IF (SQUARE) THEN
               AREA = DIAMETER**2
               PERIMETER = 4._EB*DIAMETER
            ELSEIF (ROUND) THEN
               AREA = 0.25_EB*PI*DIAMETER**2
               PERIMETER = PI*DIAMETER
            ELSE
               IF (PERIMETER < 0._EB) THEN
                  WRITE(MESSAGE,'(A,A,A,I5)') &
                     'ERROR(507): If both ROUND and SQUARE are FALSE, Duct with DIAMETER must also have PERIMETER. Duct ID:',&
                     TRIM(ID),', HVAC line number:',NN
                  CALL SHUTDOWN(MESSAGE); RETURN
               ENDIF
               AREA = 0.25_EB*DIAMETER*PERIMETER
            ENDIF
         ELSEIF (AREA > 0._EB) THEN
            IF (DIAMETER >  0._EB .AND. PERIMETER >  0._EB) THEN
               WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR(508): Duct cannot input both PERIMETER and DIAMETER with AREA. Duct ID:',&
                  TRIM(ID),', HVAC line number:',NN
               CALL SHUTDOWN(MESSAGE); RETURN
            ENDIF
            IF (SQUARE) THEN
               DIAMETER = SQRT(AREA)
               PERIMETER = 4._EB*DIAMETER
            ELSEIF (ROUND) THEN
               DIAMETER = SQRT(4._EB*AREA/PI)
               PERIMETER = PI*DIAMETER
            ELSE
               IF (DIAMETER > 0._EB) THEN
                  PERIMETER = 4._EB*AREA/DIAMETER
               ELSE
                  DIAMETER = 4._EB*AREA/PERIMETER
               ENDIF
            ENDIF
         ENDIF
         DU%AREA_INITIAL = AREA
         DU%AREA = AREA
         DU%DIAMETER = DIAMETER
         DU%LENGTH = LENGTH
         DU%REVERSE = REVERSE
         DU%DP_FAN = 0._EB
         ! Initializing to background/ambient for t=0 outputs
         DU%RHO_D = RHOA
         DU%TMP_D = TMPA
         ALLOCATE(DU%ZZ(N_TRACKED_SPECIES))
         DU%ZZ(1:N_TRACKED_SPECIES) = SPECIES_MIXTURE(1:N_TRACKED_SPECIES)%ZZ0
         ALLOCATE(DU%ZZ_OLD(N_TRACKED_SPECIES))
         DU%ZZ_OLD(1:N_TRACKED_SPECIES) = SPECIES_MIXTURE(1:N_TRACKED_SPECIES)%ZZ0
         DU%LOSS(1:2) = MAX(0._EB,LOSS(1:2,1))
         IF (CTRL_ID /='null' .AND. DEVC_ID /='null') THEN
            WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR(509): Can only specify one of CTRL_ID or DEVC_ID. Duct ID:',TRIM(ID),&
                                        ', HVAC line number:',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         IF (DAMPER .AND. (FAN_ID /='null' .OR. AIRCOIL_ID /='null') .OR. &
             FAN_ID/='null' .AND. (DAMPER .OR. AIRCOIL_ID /='null') .OR. &
             AIRCOIL_ID/='null' .AND. (DAMPER .OR. FAN_ID /='null')) THEN
            WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR(510): Duct can only have one of damper, fan or aircoil. Duct ID:',TRIM(ID),&
                                        ', HVAC line number:',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         IF (FAN_ID/='null' .AND. N_FANS<=0) THEN
            WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR(511): Duct has fan specied but no fans have been defined. Duct ID:',TRIM(ID),&
                                        ', HVAC line number: ',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         DU%DAMPER = DAMPER
         DUCT_AIRCOIL_A(I_DUCT) = AIRCOIL_ID
         DUCT_FAN_A(I_DUCT) = FAN_ID
         IF (CTRL_ID /= 'null' .OR. DEVC_ID /= 'null') &
            CALL SEARCH_CONTROLLER('HVAC',CTRL_ID,DEVC_ID,DU%DEVC_INDEX,DU%CTRL_INDEX,NN)
         IF (DAMPER) THEN
            IF (DU%DEVC_INDEX > 0) THEN
                DU%DAMPER_OPEN = DEVICE(DU%DEVC_INDEX)%INITIAL_STATE
            ELSEIF (DU%CTRL_INDEX > 0) THEN
                DU%DAMPER_OPEN = CONTROL(DU%CTRL_INDEX)%INITIAL_STATE
            ELSE
                DU%DAMPER_OPEN = .TRUE.
            ENDIF
            IF (.NOT. DU%DAMPER_OPEN) DU%AREA = 0._EB
         ELSEIF (FAN_ID /='null') THEN
            IF (DU%DEVC_INDEX > 0) THEN
                DU%FAN_OPERATING = DEVICE(DU%DEVC_INDEX)%INITIAL_STATE
            ELSEIF (DU%CTRL_INDEX > 0) THEN
                DU%FAN_OPERATING = CONTROL(DU%CTRL_INDEX)%INITIAL_STATE
            ELSE
                DU%FAN_OPERATING = .TRUE.
            ENDIF
            IF (DU%FAN_OPERATING) DU%FAN_ON_TIME = T_BEGIN
         ELSEIF (AIRCOIL_ID /='null') THEN
            IF (DU%DEVC_INDEX > 0) THEN
                DU%COIL_OPERATING = DEVICE(DU%DEVC_INDEX)%INITIAL_STATE
            ELSEIF (DU%CTRL_INDEX > 0) THEN
                DU%COIL_OPERATING = CONTROL(DU%CTRL_INDEX)%INITIAL_STATE
            ELSE
                DU%COIL_OPERATING = .TRUE.
            ENDIF
            IF (DU%COIL_OPERATING) DU%COIL_ON_TIME = T_BEGIN
         ENDIF
         DUCT_NODE_A(I_DUCT,:) = NODE_ID
         IF (VOLUME_FLOW < 1.E7_EB .AND. MASS_FLOW < 1.E7_EB) THEN
            WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR(512): Duct has both MASS_FLOW and VOLUME_FLOW defined. Duct ID:',TRIM(DU%ID),&
                                        ', HVAC line number:',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         DU%MASS_FLOW_INITIAL = MASS_FLOW
         DU%VOLUME_FLOW_INITIAL = VOLUME_FLOW
         DU%ROUGHNESS = ROUGHNESS
         DU%TAU = TAU_VF
         IF (TAU_VF > 0._EB) DU%RAMP_INDEX = TANH_RAMP
         IF (TAU_VF < 0._EB) DU%RAMP_INDEX = TSQR_RAMP
         IF (RAMP_ID /= 'null') CALL GET_RAMP_INDEX(RAMP_ID,'DUCT',DU%RAMP_INDEX)
         IF (RAMP_LOSS /= 'null') CALL GET_RAMP_INDEX(RAMP_LOSS,'DUCT',DU%RAMP_LOSS_INDEX)
         IF (N_CELLS > 0) THEN
            HVAC_MASS_TRANSPORT = .TRUE.
            DU%N_CELLS = N_CELLS
         ENDIF

         IF (ANY(WAYPOINTS > -HUGE(1._EB))) THEN
            DU%N_WAYPOINTS = 0
            DO NC=1,30
               IF (ALL(WAYPOINTS(NC,:)>-HUGE(1._EB))) THEN
                  DU%N_WAYPOINTS = DU%N_WAYPOINTS + 1
               ELSE
                  EXIT
               ENDIF
            ENDDO
            ALLOCATE(DU%WAYPOINT_XYZ(DU%N_WAYPOINTS,3))
            DO NC=1,DU%N_WAYPOINTS
               DU%WAYPOINT_XYZ(NC,:) = WAYPOINTS(NC,:)
            ENDDO
         ENDIF

         IF (NETWORK_ID=='null') THEN
            DU%NETWORK_ID='Unassigned'
         ELSE
            DU%NETWORK_ID=NETWORK_ID
         ENDIF

      CASE('NODE')
         I_DUCTNODE = I_DUCTNODE + 1
         NODE_DUCT_A(I_DUCTNODE,:) = DUCT_ID
         NODE_FILTER_A(I_DUCTNODE) = FILTER_ID
         DN => DUCTNODE(I_DUCTNODE)
         DN%ID = ID
         IF (TRIM(ID)=='null' ) THEN
            WRITE(MESSAGE,'(A,I5)') 'ERROR(513): Ductnode has no ID, HVAC line number:',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         IF ((GEOM .OR. GEOM2) .AND. VENT_ID/='null') THEN
            WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR(568): Ductnode with GEOM cannot have a VENT_ID. Ductnode ID:',TRIM(DN%ID),&
                                        ', HVAC line number:',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         IF (.NOT. AMBIENT .AND. .NOT. GEOM .AND. VENT_ID=='null' .AND. N_DUCTS==1) THEN
            WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR(569): Ductnode with one duct needs either AMBIENT, GEOM, or VENT_ID. Ductnode ID:',&
                                        TRIM(DN%ID),', HVAC line number:',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         DN%VENT_ID = VENT_ID
         DN%GEOM = GEOM
         IF (DN%VENT .AND. DN%GEOM) THEN
            WRITE(MESSAGE,'(A,A,A,A)') 'ERROR(570): Problem with ductnode:',TRIM(DN%ID), &
                                          ', cannot assign to both VENT and GEOM.'
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         DN%READ_IN = .TRUE.
         IF (TRIM(VENT_ID)/='null') DN%VENT=.TRUE.
         ! Set temporary node elevation for vents so WIND can establish T and P RAMPs
         IF ((DN%VENT .OR. DN%GEOM).AND. XYZ(3) <-1.E9) XYZ(3) = ZS_MIN
         IF (.NOT. (DN%VENT .OR. DN%GEOM) .AND. XYZ(3)<-1.E9) THEN
            WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR(514): Ambient or internal ductnode requires an elevation, XYZ(3). Ductnode ID:',&
                                        TRIM(DN%ID),', HVAC line number:',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         DN%XYZ      = XYZ
         IF (ALL(DN%XYZ>-1.E9) .OR. DN%VENT .OR. DN%GEOM) DN%SPECIFIED_XYZ = .TRUE.
         IF (.NOT. (DN%VENT .OR. DN%GEOM)) THEN
            IF (DN%XYZ(1) <=-1.E9_EB) DN%XYZ(1) = XS_MIN
            IF (DN%XYZ(2) <=-1.E9_EB) DN%XYZ(2) = YS_MIN
         ENDIF
         DN%AMBIENT  = AMBIENT
         DO ND = 1, MAX_DUCTS
            IF (NODE_DUCT_A(I_DUCTNODE,ND) == 'null') EXIT
            DN%N_DUCTS=ND
         ENDDO
         IF (DN%N_DUCTS == 1 .AND. .NOT. GEOM .AND. .NOT. AMBIENT .AND. VENT_ID=='null') THEN
            WRITE(MESSAGE,'(A,A,A,A,I5)') 'ERROR(515): Non-AMBIENT, non-VENT, and non-GEOM connected ductnode ',&
               '(i.e., internal node) must have >=2 ducts. Ductnode ID:',TRIM(DN%ID),', HVAC line number:',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         IF (DN%N_DUCTS >= 2 .AND. (AMBIENT .OR. GEOM .OR. VENT_ID/='null')) THEN
            WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR(516): AMBIENT, VENT, or GEOM-connected ductnode must have 1 duct. Ductnode ID:',&
                                        TRIM(DN%ID),', HVAC line number:',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         IF (DN%N_DUCTS == 0) THEN
            WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR(517): No ducts specified for ductnode ID:',TRIM(DN%ID),', HVAC line number ',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         IF (DN%N_DUCTS /= 2 .AND. TRIM(FILTER_ID)/='null') THEN
            WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR(518): Ductnode with a filter must have 2 ducts. Ductnode ID:',TRIM(DN%ID),&
                                        ', HVAC line number:',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         ALLOCATE(DN%LOSS_ARRAY(MAX(2,DN%N_DUCTS),MAX(2,DN%N_DUCTS)))
         DN%LOSS_ARRAY = 0._EB
         IF (DN%N_DUCTS >=2) THEN
            DN%LOSS_ARRAY = LOSS(1:DN%N_DUCTS,1:DN%N_DUCTS)
         ELSE
            DN%LOSS_ARRAY(1,2) = LOSS(2,1)
            DN%LOSS_ARRAY(2,1) = LOSS(1,1)
         ENDIF

         IF (NETWORK_ID=='null') THEN
            DN%NETWORK_ID='Unassigned'
         ELSE
            DN%NETWORK_ID=NETWORK_ID
         ENDIF

      CASE('FAN')
         I_FAN = I_FAN + 1
         FAN(I_FAN)%ID = ID
         IF (TRIM(ID)=='null' ) THEN
            WRITE(MESSAGE,'(A,I5)') 'ERROR(519): Fan has no ID, HVAC line number:',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         FAN(I_FAN)%OFF_LOSS = LOSS(1,1)
         IF (FAN(I_FAN)%OFF_LOSS<TWO_EPSILON_EB) FAN(I_FAN)%OFF_LOSS = 1._EB
         FAN(I_FAN)%FAN_RAMP = RAMP_ID
         FAN(I_FAN)%VOL_FLOW = VOLUME_FLOW
         FAN(I_FAN)%MAX_PRES = MAX_PRESSURE
         FAN(I_FAN)%MAX_FLOW = MAX_FLOW
         FAN(I_FAN)%TAU = TAU_FAN
         IF (TAU_FAN > 0._EB) FAN(I_FAN)%SPIN_INDEX = TANH_RAMP
         IF (TAU_FAN < 0._EB) FAN(I_FAN)%SPIN_INDEX = TSQR_RAMP
         IF (RAMP_ID /= 'null') CALL GET_RAMP_INDEX(RAMP_ID,'FAN',FAN(I_FAN)%RAMP_INDEX)
         IF(( (MAX_FLOW<1.E6_EB .OR. MAX_PRESSURE<1.E6_EB) .AND. (VOLUME_FLOW<1.E6_EB .OR. RAMP_ID/='null')))THEN !.OR. &
            WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR(520): FAN can only be one of constant volume, quadratic or ramp. Fan ID:',TRIM(ID),&
                                        ', HVAC line number:',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         IF ((MAX_PRESSURE<1.E6_EB .AND. MAX_FLOW>1.E6_EB) .OR. (MAX_PRESSURE>1.E6_EB .AND. MAX_FLOW<1.E6_EB)) THEN
            WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR(521): IF one of MAX_PRESSURE or MAX_FLOW given, both must be specified. Fan ID:',&
                                        TRIM(ID),', HVAC line number:',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         IF (MAX_PRESSURE <= 0._EB) THEN
            WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR(522): MAX_PRESSURE must be > 0. Fan ID:',TRIM(ID),', HVAC line number:',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         IF (MAX_FLOW <= 0._EB) THEN
            WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR(523): MAX_FLOW must be > 0. Fan ID:',TRIM(ID),', HVAC line number:',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         IF (VOLUME_FLOW < 1.E6_EB) THEN
            FAN(I_FAN)%FAN_TYPE = 1
         ELSEIF(RAMP_ID/='null') THEN
            FAN(I_FAN)%FAN_TYPE = 3
         ELSE
            FAN(I_FAN)%FAN_TYPE = 2
         ENDIF

      CASE('FILTER')
         I_FILTER = I_FILTER + 1
         FILTER(I_FILTER)%ID = ID
         IF (TRIM(ID)=='null' ) THEN
            WRITE(MESSAGE,'(A,I5)') 'ERROR(524): Filter has no ID, HVAC line number:',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         FILTER(I_FILTER)%CLEAN_LOSS = CLEAN_LOSS
         IF (TRIM(RAMP_ID)/='null') CALL GET_RAMP_INDEX(RAMP_ID,'FILTER',FILTER(I_FILTER)%RAMP_INDEX)
         ALLOCATE(FILTER(I_FILTER)%EFFICIENCY(1:N_TRACKED_SPECIES))
         FILTER(I_FILTER)%EFFICIENCY = 0._EB
         ALLOCATE(FILTER(I_FILTER)%MULTIPLIER(1:N_TRACKED_SPECIES))
         FILTER(I_FILTER)%MULTIPLIER = 0._EB
         ALLOCATE(FILTER(I_FILTER)%INITIAL_LOADING(1:N_TRACKED_SPECIES))
         FILTER(I_FILTER)%INITIAL_LOADING = 0._EB
         FILTER(I_FILTER)%LOADING_LOSS = LOSS(1,1)
         SPEC_LOOP2: DO N=1,N_TRACKED_SPECIES
            IF (TRIM(SPEC_ID(N))=='null') EXIT SPEC_LOOP2
            DO NS = 1,N_TRACKED_SPECIES
               IF (TRIM(SPECIES_MIXTURE(NS)%ID)==TRIM(SPEC_ID(N))) THEN
                  FILTER(I_FILTER)%EFFICIENCY(NS)   = EFFICIENCY(N)
                  FILTER(I_FILTER)%MULTIPLIER(NS) = LOADING_MULTIPLIER(N)
                  FILTER(I_FILTER)%INITIAL_LOADING(NS) = LOADING(N)
                  EXIT
               ENDIF
               IF (NS==N_TRACKED_SPECIES) THEN
                  WRITE(MESSAGE,'(A,A,A,A,A)') 'ERROR(525): Problem with filter:',TRIM(ID),'. SPEC ',TRIM(SPEC_ID(N)),' not found'
                  CALL SHUTDOWN(MESSAGE); RETURN
               ENDIF
            ENDDO
         ENDDO SPEC_LOOP2
         FILTER(I_FILTER)%AREA = AREA
      CASE('AIRCOIL')
         I_AIRCOIL = I_AIRCOIL+1
         AIRCOIL(I_AIRCOIL)%COOLANT_SPECIFIC_HEAT   = COOLANT_SPECIFIC_HEAT*1000._EB
         AIRCOIL(I_AIRCOIL)%COOLANT_MASS_FLOW = COOLANT_MASS_FLOW
         AIRCOIL(I_AIRCOIL)%COOLANT_TEMPERATURE = COOLANT_TEMPERATURE+TMPM
         AIRCOIL(I_AIRCOIL)%EFFICIENCY   = EFFICIENCY(1)
         AIRCOIL(I_AIRCOIL)%FIXED_Q      = FIXED_Q*1000._EB
         AIRCOIL(I_AIRCOIL)%ID           = ID
         IF (TRIM(ID)=='null' ) THEN
            WRITE(MESSAGE,'(A,I5)') 'ERROR(526): Aircoil has no ID, HVAC line number:',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         AIRCOIL(I_AIRCOIL)%TAU          = TAU_AC
         IF (TAU_AC > 0._EB) AIRCOIL(I_AIRCOIL)%RAMP_INDEX = TANH_RAMP
         IF (TAU_AC < 0._EB) AIRCOIL(I_AIRCOIL)%RAMP_INDEX = TSQR_RAMP
         IF (RAMP_ID /= 'null') CALL GET_RAMP_INDEX(RAMP_ID,'DUCT',AIRCOIL(I_AIRCOIL)%RAMP_INDEX)
      CASE('LEAK')
         I_DUCTNODE = I_DUCTNODE + 1
         NODE_DUCT_A(I_DUCTNODE,1) = ID
         NODE_FILTER_A(I_DUCTNODE) = 'null'
         DN => DUCTNODE(I_DUCTNODE)
         DN%ID = VENT_ID
         IF (TRIM(ID)=='null' ) THEN
            WRITE(MESSAGE,'(A,I5)') 'ERROR(527): Localized leakage has no ID, HVAC line number:',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         IF (GEOM) THEN
            DN%GEOM = .TRUE.
            IF (TRIM(DN%VENT_ID)=='AMBIENT') THEN
               WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR(571): VENT_ID for leakage cannot be AMBIENT if GEOM is set. Leak ID:',TRIM(ID),&
                  ', HVAC line number:',NN
               CALL SHUTDOWN(MESSAGE); RETURN
            ENDIF
         ELSE
            DN%VENT_ID = VENT_ID
            DN%VENT=.TRUE.
         ENDIF
         DN%NETWORK_ID='LEAK'
         DN%READ_IN = .FALSE.
         DN%TRANSPORT_PARTICLES = TRANSPORT_PARTICLES
         IF (TRIM(DN%VENT_ID)=='null' .AND. .NOT. GEOM) THEN
            WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR(528): Leakage path must have VENT_ID defined. Leak ID:',TRIM(ID),&
               ', HVAC line number:',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         IF (TRIM(DN%VENT_ID)=='AMBIENT') THEN
            WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR(529): Leakage to AMBIENT must have VENT2_ID for the AMBIENT node. Leak ID:',&
               TRIM(ID),', HVAC line number:',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         DN%N_DUCTS=1
         ALLOCATE(DN%LOSS_ARRAY(2,2))
         DN%LOSS_ARRAY = 0._EB

         I_DUCTNODE = I_DUCTNODE + 1
         NODE_DUCT_A(I_DUCTNODE,1) = ID
         NODE_FILTER_A(I_DUCTNODE) = 'null'
         DN => DUCTNODE(I_DUCTNODE)
         DN%ID = VENT2_ID
         IF (GEOM2) THEN
            DN%GEOM = .TRUE.
            IF (TRIM(DN%VENT_ID)=='AMBIENT') THEN
               WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR(572): VENT2_ID for leakage cannot be AMBIENT if GEOM2 is set. Leak ID:',&
                  TRIM(ID),', HVAC line number:',NN
               CALL SHUTDOWN(MESSAGE); RETURN
            ENDIF
         ELSE
            DN%VENT_ID = VENT2_ID
            DN%VENT=.TRUE.
         ENDIF
         DN%VENT=.TRUE.
         DN%NETWORK_ID='LEAK'
         DN%READ_IN = .FALSE.
         DN%TRANSPORT_PARTICLES = TRANSPORT_PARTICLES
         IF (TRIM(VENT2_ID)=='null' .AND. .NOT. GEOM2) THEN
            WRITE(MESSAGE,'(A,A,A,I2)') 'ERROR(530): Leakage path must have VENT2_ID defined. Leak ID:',TRIM(ID),&
                                        ', HVAC line number:',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         ! Set temporary node elevation for vents so WIND can establish T and P RAMPs
         IF (XYZ(3) <-1.E9) XYZ(3) = ZS_MIN
         DN%XYZ      = XYZ
         DN%AMBIENT  = .FALSE.
         IF (TRIM(DN%VENT_ID)=='AMBIENT') THEN
            DN%AMBIENT = .TRUE.
            DN%ID = TRIM(VENT_ID)//' AMB'
            VENT2_ID = DN%ID
            DN%VENT_ID ='null'
            DN%VENT = .FALSE.
         ENDIF
         DN%N_DUCTS=1
         ALLOCATE(DN%LOSS_ARRAY(2,2))
         DN%LOSS_ARRAY = 0._EB

         I_DUCT = I_DUCT + 1
         DU=> DUCT(I_DUCT)
         DU%ID   = ID
         DU%LOCALIZED_LEAKAGE = .TRUE.
         IF (AREA <= 0._EB) THEN
            WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR(531): Leakage has no AREA. Leak ID:',TRIM(ID),', HVAC line number:',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         DU%AREA_INITIAL = AREA
         DU%AREA = AREA
         DU%DIAMETER = -1._EB
         DU%LENGTH = 0.1_EB
         DU%REVERSE = .FALSE.
         ! Initializing to background/ambient for t=0 outputs
         DU%RHO_D = RHOA
         DU%TMP_D = TMPA
         ALLOCATE(DU%ZZ(N_TRACKED_SPECIES))
         DU%ZZ(1:N_TRACKED_SPECIES) = SPECIES_MIXTURE(1:N_TRACKED_SPECIES)%ZZ0
         ALLOCATE(DU%ZZ_OLD(N_TRACKED_SPECIES))
         DU%ZZ_OLD(1:N_TRACKED_SPECIES) = SPECIES_MIXTURE(1:N_TRACKED_SPECIES)%ZZ0
         IF (LOSS(1,1)==0._EB) LOSS(1,1)=1._EB
         DU%LOSS(1:2) = MAX(0._EB,LOSS(1,1))
         DU%DAMPER = .FALSE.
         DUCT_AIRCOIL_A(I_DUCT) = 'null'
         DUCT_FAN_A(I_DUCT) = 'null'
         DUCT_NODE_A(I_DUCT,1) = VENT_ID
         DUCT_NODE_A(I_DUCT,2) = VENT2_ID
         DU%MASS_FLOW_INITIAL = 1.E7_EB
         DU%VOLUME_FLOW_INITIAL = 1.E7_EB
         DU%ROUGHNESS = -1._EB
         DU%TAU = 1._EB
         DU%RAMP_INDEX = TANH_RAMP
         DU%LEAKAGE = .FALSE.
         DU%LEAK_ENTHALPY = LEAK_ENTHALPY
         DU%LEAK_REFERENCE_PRESSURE = LEAK_REFERENCE_PRESSURE
         DU%LEAK_PRESSURE_EXPONENT = LEAK_PRESSURE_EXPONENT
         DU%DISCHARGE_COEFFICIENT = DISCHARGE_COEFFICIENT
         DU%NETWORK_ID='LEAK'
      CASE ('DUCT QUANTITY LIST')
         IF (DUCT_QUANTITY_DEFINED) THEN
            WRITE(MESSAGE,'(A,A,I5)') 'ERROR(532): Can only have one HVAC input with TYPE_ID of DUCT QUANTITY LIST.',&
               ' HVAC line number:',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         DUCT_QUANTITY_DEFINED = .TRUE.
         DO N=1,20
            IF (QUANTITY(N)=='null') EXIT
            N_DUCT_QUANTITY = N_DUCT_QUANTITY + 1
         ENDDO
         ALLOCATE (DUCT_QUANTITY_ARRAY(N_DUCT_QUANTITY))
         DUCT_QUANTITY_ARRAY(1:N_DUCT_QUANTITY)%DRY = DRY(1:N_DUCT_QUANTITY)
         DO N=1, N_DUCT_QUANTITY
            HQT => DUCT_QUANTITY_ARRAY(N)
            CALL GET_QUANTITY_INDEX(HQT%SMOKEVIEW_LABEL,HQT%SMOKEVIEW_BAR_LABEL,HQT%OUTPUT_INDEX,HQT%Y_INDEX,HQT%Z_INDEX,&
               HQT%UNITS,QUANTITY(N),QUANTITY_SPEC_ID(N))
         ENDDO
      CASE ('NODE QUANTITY LIST')
         IF (NODE_QUANTITY_DEFINED) THEN
            WRITE(MESSAGE,'(A,A,I5)') 'ERROR(533): Can only have one HVAC input with TYPE_ID of NODE QUANTITY LIST.',&
               ' HVAC line number:',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         NODE_QUANTITY_DEFINED = .TRUE.
         DO N=1,20
            IF (QUANTITY(N)=='null') EXIT
            N_NODE_QUANTITY = N_NODE_QUANTITY + 1
         ENDDO
         ALLOCATE (NODE_QUANTITY_ARRAY(N_NODE_QUANTITY))
         NODE_QUANTITY_ARRAY(1:N_NODE_QUANTITY)%DRY = DRY(1:N_NODE_QUANTITY)         
         DO N=1, N_NODE_QUANTITY
            HQT => NODE_QUANTITY_ARRAY(N)
            CALL GET_QUANTITY_INDEX(HQT%SMOKEVIEW_LABEL,HQT%SMOKEVIEW_BAR_LABEL,HQT%OUTPUT_INDEX,HQT%Y_INDEX,HQT%Z_INDEX,&
               HQT%UNITS,QUANTITY(N),QUANTITY_SPEC_ID(N))
         ENDDO
   END SELECT
ENDDO

NODE_Z_MIN =  HUGE(EB)
NODE_Z_MAX = -HUGE(EB)
DO N=1,N_DUCTNODES
   NODE_Z_MIN = MIN(NODE_Z_MIN,DUCTNODE(N)%XYZ(3))
   NODE_Z_MAX = MAX(NODE_Z_MAX,DUCTNODE(N)%XYZ(3))
ENDDO

T_USED(13)=T_USED(13)+CURRENT_TIME()-TNOW

RETURN

CONTAINS

!> \brief Sets the default values for the HVAC namelist

SUBROUTINE SET_HVAC_DEFAULTS

AIRCOIL_ID   = 'null'
AMBIENT      = .FALSE.
AREA         = -1._EB
COOLANT_SPECIFIC_HEAT   = 4.186_EB
COOLANT_MASS_FLOW = -1.E10_EB
COOLANT_TEMPERATURE = 20._EB
DISCHARGE_COEFFICIENT = 1._EB
EFFICIENCY   = 1.0_EB
CLEAN_LOSS   = 0._EB
CTRL_ID      = 'null'
DAMPER       = .FALSE.
DEVC_ID      = 'null'
DIAMETER     = -1._EB
DRY          = .FALSE.
DUCT_ID      = 'null'
FAN_ID       = 'null'
FIXED_Q      = -1.E10_EB
FILTER_ID    = 'null'
GEOM         = .FALSE.
GEOM2        = .FALSE.
LEAK_ENTHALPY = .FALSE.
LEAK_PRESSURE_EXPONENT = 0.5_EB
LEAK_REFERENCE_PRESSURE = 4._EB
LENGTH       = -1._EB
LOADING      = 0._EB
LOADING_MULTIPLIER = 1._EB
LOSS         = 0._EB
MASS_FLOW    = 1.E7_EB
MAX_FLOW     = 1.E7_EB
MAX_PRESSURE = 1.E7_EB
N_CELLS      = -999
NETWORK_ID   = 'null'
NODE_ID      = 'null'
QUANTITY = 'null'
QUANTITY_SPEC_ID = 'null'
INITIALIZED_HVAC_MASS_TRANSPORT=.FALSE.
PERIMETER    = -1._EB
RAMP_ID      = 'null'
RAMP_LOSS    = 'null'
REVERSE      = .FALSE.
ROUGHNESS    = -1._EB
ROUND        = .TRUE.
SPEC_ID      = 'null'
SQUARE       = .FALSE.
TYPE_ID      = 'null'
TAU_AC       = TAU_DEFAULT
TAU_FAN      = TAU_DEFAULT
TAU_VF       = TAU_DEFAULT
TRANSPORT_PARTICLES = .FALSE.
VENT_ID      = 'null'
VENT2_ID     = 'null'
VOLUME_FLOW  = 1.E7_EB
WAYPOINTS    = -HUGE(1._EB)
XYZ          = -1.E10_EB

RETURN

END SUBROUTINE SET_HVAC_DEFAULTS

END SUBROUTINE READ_HVAC


!\brief Builds the HVAC network linking together the various types of HVAC inputs

SUBROUTINE PROC_HVAC

USE PHYSICAL_FUNCTIONS, ONLY: GET_ENTHALPY
USE MPI_F08
INTEGER :: N,ND,ND2,NM,NN,NF,NV,IERR,NG,FI(3)
REAL(EB) :: TNOW !< Current CPU time (s) used in computing length of time spent in HVAC routines.
REAL(EB) :: ZZ_GET(1:N_TRACKED_SPECIES) !< Species mass fraction array
REAL(EB) :: CF_AREA(N_DUCTNODES),CF_X(N_DUCTNODES),CF_Y(N_DUCTNODES),CF_Z(N_DUCTNODES),AREA,X,Y,Z,XYZ_F(3,3)
LOGICAL :: FOUND !< Flag indicating search loop has found the HVAC component assocaited with an ID
LOGICAL :: DUMMY = .FALSE. !< Dummy flag
TYPE (LAGRANGIAN_PARTICLE_CLASS_TYPE),DIMENSION(:), POINTER:: TEMPALLOC
TYPE(DUCTNODE_TYPE), POINTER :: DN
TYPE(DUCT_TYPE), POINTER :: DU
TYPE(SURFACE_TYPE), POINTER :: SF
TYPE(GEOMETRY_TYPE), POINTER :: G

TNOW=CURRENT_TIME()

IF (.NOT. HVAC_SOLVE) RETURN

DUCT_LOOP: DO ND = 1, N_DUCTS
   DU => DUCT(ND)
   IF (DU%LEAKAGE) CYCLE DUCT_LOOP
   IF (TRIM(DUCT_NODE_A(ND,1))==TRIM(DUCT_NODE_A(ND,2))) THEN
      WRITE(MESSAGE,'(A,A)') 'ERROR(534): Both nodes have the same ID. Duct ID:',TRIM(DU%ID)
      CALL SHUTDOWN(MESSAGE); RETURN
   ENDIF
   DO N = 1, ND
      IF (N==ND) CYCLE
      IF (TRIM(DU%ID)==TRIM(DUCT(N)%ID)) THEN
         WRITE(MESSAGE,'(A,A)') 'ERROR(535): Two ducts with the same ID. Duct ID:',TRIM(DU%ID)
         CALL SHUTDOWN(MESSAGE); RETURN
      ENDIF
   ENDDO
   NODELOOP_D: DO NN = 1, N_DUCTNODES
      IF(TRIM(DUCTNODE(NN)%ID) == TRIM(DUCT_NODE_A(ND,1))) THEN
         DU%NODE_INDEX(1) = NN
         FOUND = .FALSE.
         DN1_NAME: DO N=1,DUCTNODE(NN)%N_DUCTS
            IF (TRIM(DU%ID)==TRIM(NODE_DUCT_A(NN,N))) THEN
               FOUND = .TRUE.
               EXIT DN1_NAME
            ENDIF
         ENDDO DN1_NAME
         IF (.NOT. FOUND) THEN
            WRITE(MESSAGE,'(A,A,A,A,A)') 'ERROR(536): Duct: ',TRIM(DU%ID),'. Node:',&
                                      TRIM(DUCTNODE(NN)%ID),' does not contain the duct in its list of ducts.'
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
      ENDIF
      IF(TRIM(DUCTNODE(NN)%ID) == TRIM(DUCT_NODE_A(ND,2))) THEN
         DU%NODE_INDEX(2) = NN
         FOUND = .FALSE.
         DN2_NAME: DO N=1,DUCTNODE(NN)%N_DUCTS
            IF (TRIM(DU%ID)==TRIM(NODE_DUCT_A(NN,N))) THEN
               FOUND = .TRUE.
               EXIT DN2_NAME
            ENDIF
         ENDDO DN2_NAME
         IF (.NOT. FOUND) THEN
            WRITE(MESSAGE,'(A,A,A,A,A)') 'ERROR(536): Duct: ',TRIM(DU%ID),'. Node:',&
                                      TRIM(DUCTNODE(NN)%ID),' does not contain the duct in its list of ducts.'
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
      ENDIF
      IF (DU%NODE_INDEX(1) > 0 .AND. DU%NODE_INDEX(2) > 0) EXIT NODELOOP_D
   ENDDO NODELOOP_D
   IF (DU%NODE_INDEX(1) <= 0) THEN
      WRITE(MESSAGE,'(A,I3,A,A)') 'ERROR(537): First ductnode not located for duct:',ND,', Duct ID:',TRIM(DU%ID)
      CALL SHUTDOWN(MESSAGE); RETURN
   ENDIF
   IF (DU%NODE_INDEX(2) <= 0) THEN
      WRITE(MESSAGE,'(A,I3,A,A)') 'ERROR(537): Second ductnode not located for duct:',ND,', Duct ID:',TRIM(DU%ID)
      CALL SHUTDOWN(MESSAGE); RETURN
   ENDIF
   IF (DUCT_FAN_A(ND)/='null') THEN
      DU%DP_FAN = 0._EB
      DO NF = 1, N_FANS
         IF(TRIM(FAN(NF)%ID) == TRIM(DUCT_FAN_A(ND))) THEN
            DU%FAN_INDEX = NF
            EXIT
         ENDIF
      ENDDO
      IF (DU%FAN_INDEX <= 0) THEN
         WRITE(MESSAGE,'(A,I3,A,A)') 'ERROR(538): Fan not located. Duct:',ND,', Duct ID:',TRIM(DU%ID)
         CALL SHUTDOWN(MESSAGE); RETURN
      ENDIF
   ENDIF
   IF (DUCT_AIRCOIL_A(ND)/='null') THEN
      DU%COIL_Q = 0._EB
      DO NF = 1, N_AIRCOILS
         IF(TRIM(AIRCOIL(NF)%ID) == TRIM(DUCT_AIRCOIL_A(ND))) THEN
            DU%AIRCOIL_INDEX = NF
            EXIT
         ENDIF
      ENDDO
      IF (DU%AIRCOIL_INDEX <= 0) THEN
         WRITE(MESSAGE,'(A,I3,A,A)') 'ERROR(539): Aircoil not located. Duct:',ND,', Duct ID:',TRIM(DU%ID)
         CALL SHUTDOWN(MESSAGE); RETURN
      ENDIF
   ENDIF
   IF (DU%SURF_INDEX> 1) THEN
      IF (N_LAGRANGIAN_CLASSES<1) THEN
         N_LAGRANGIAN_CLASSES = 1
         ALLOCATE(LAGRANGIAN_PARTICLE_CLASS(N_LAGRANGIAN_CLASSES))
         DU%LPC_INDEX = 1
      ELSE
         CLASS_LOOP: DO N= 1,N_LAGRANGIAN_CLASSES
            IF (LAGRANGIAN_PARTICLE_CLASS(N)%SURF_INDEX==DU%SURF_INDEX .AND. LAGRANGIAN_PARTICLE_CLASS(N)%DUCT_PARTICLE) THEN
               DU%LPC_INDEX=N
               EXIT CLASS_LOOP
            ENDIF
         ENDDO CLASS_LOOP
         IF (DU%LPC_INDEX < 0) THEN
            ALLOCATE(TEMPALLOC(N_LAGRANGIAN_CLASSES))
            TEMPALLOC = LAGRANGIAN_PARTICLE_CLASS
            DEALLOCATE(LAGRANGIAN_PARTICLE_CLASS)
            ALLOCATE(LAGRANGIAN_PARTICLE_CLASS(N_LAGRANGIAN_CLASSES+1))
            LAGRANGIAN_PARTICLE_CLASS(1:N_LAGRANGIAN_CLASSES) = TEMPALLOC(1:N_LAGRANGIAN_CLASSES)
            N_LAGRANGIAN_CLASSES=N_LAGRANGIAN_CLASSES+1
         ENDIF
         DU%LPC_INDEX = N_LAGRANGIAN_CLASSES
      ENDIF
   ENDIF
      ! Add stuff to initialize class as a duct particle
ENDDO DUCT_LOOP

NODE_LOOP: DO NN = 1, N_DUCTNODES

   DN => DUCTNODE(NN)
   DO N = 1, NN
      IF (N==NN) CYCLE
      IF (TRIM(DN%ID)==TRIM(DUCTNODE(N)%ID)) THEN
         WRITE(MESSAGE,'(A,A)') 'ERROR(540): Two duct nodes with the same ID. Ductnode ID:',TRIM(DN%ID)
         CALL SHUTDOWN(MESSAGE); RETURN
      ENDIF
   ENDDO

   ! Initializes duct node species and RSUM with ambient/background

   ALLOCATE(DN%ZZ(N_TRACKED_SPECIES))
   DN%ZZ(1:N_TRACKED_SPECIES) = SPECIES_MIXTURE(1:N_TRACKED_SPECIES)%ZZ0
   ALLOCATE(DN%ZZ_OLD(N_TRACKED_SPECIES))
   DN%ZZ_OLD(1:N_TRACKED_SPECIES) = SPECIES_MIXTURE(1:N_TRACKED_SPECIES)%ZZ0
   ALLOCATE(DN%ZZ0(N_TRACKED_SPECIES))
   ALLOCATE(DN%ZZ_V(N_TRACKED_SPECIES))
   DN%ZZ_V(1:N_TRACKED_SPECIES) = SPECIES_MIXTURE(1:N_TRACKED_SPECIES)%ZZ0
   ZZ_GET(1:N_TRACKED_SPECIES) = DN%ZZ_V(1:N_TRACKED_SPECIES)
   DN%RSUM   = RSUM0
   ALLOCATE(DN%ZZ_SUM(N_TRACKED_SPECIES))
   DN%ZZ_SUM = 0._EB
   DN%E_SUM = 0._EB

   ! If node is LEAKAGE-related then values are adopted as ambient/background

   IF (DN%LEAKAGE) THEN
      DN%TMP  = TMPA
      DN%RHO  = RHOA
      DN%P    = P_INF
      DN%P_OLD = DN%P
      CALL GET_ENTHALPY(ZZ_GET,DN%CP,TMPA)
      DN%CP = DN%CP / TMPA
      DN%TMP_V  = DN%TMP
      DN%RSUM_V = DN%RSUM
      DN%CP_V   = DN%CP
      DN%RHO_V  = DN%RHO
      CYCLE NODE_LOOP
   ENDIF

   ! If the duct node has a VENT associated with it, find it

   IF (DN%VENT_ID /= 'null') THEN
      FOUND = .FALSE.
      MESH_LOOP: DO NM = 1, NMESHES
         IF (PROCESS(NM)/=MY_RANK)   CYCLE MESH_LOOP  ! Only search meshes controlled by the current MPI process
         NODE_VENT_LOOP:DO NV = 1, MESHES(NM)%N_VENT
            IF(MESHES(NM)%VENTS(NV)%ID == DN%VENT_ID) THEN
               FOUND = .TRUE.
               IF (MESHES(NM)%VENTS(NV)%CTRL_INDEX > 0 .OR. MESHES(NM)%VENTS(NV)%DEVC_INDEX >0) THEN
                  WRITE(MESSAGE,'(A,A)') 'ERROR(541): VENT for ductnode has a DEVC_ID or CTRL_ID, VENT ID:',&
                                          TRIM(MESHES(NM)%VENTS(NV)%ID)
                  CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.)
               ENDIF
               IF (DN%READ_IN .AND. MESHES(NM)%VENTS(NV)%SURF_INDEX /= HVAC_SURF_INDEX) THEN
                  WRITE(MESSAGE,'(A,A)') 'ERROR(542): Ductnode attached to VENT without SURF_ID HVAC for VENT ID:',&
                                          TRIM(MESHES(NM)%VENTS(NV)%ID)
                  CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.)
               ENDIF
               IF (DN%NETWORK_ID=='LEAK') THEN
                  IF (MESHES(NM)%VENTS(NV)%SURF_INDEX == HVAC_SURF_INDEX) THEN
                     WRITE(MESSAGE,'(3A)') 'ERROR(562): VENT ID:',TRIM(MESHES(NM)%VENTS(NV)%ID),&
                        ' used for localized leakage has SURF_ID HVAC'
                     CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.)
                  ENDIF
                  IF (MESHES(NM)%VENTS(NV)%DEVC_INDEX > 0 .OR. MESHES(NM)%VENTS(NV)%CTRL_INDEX > 0) THEN
                     WRITE(MESSAGE,'(3A)') 'ERROR(563): VENT ID:',TRIM(MESHES(NM)%VENTS(NV)%ID),&
                        ' used for localized leakage has a DEVC_ID or CTRL_ID.'
                     CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.)
                  ENDIF
               ENDIF 
               IF (MESHES(NM)%VENTS(NV)%BOUNDARY_TYPE/=HVAC_BOUNDARY) THEN
                  SF => SURFACE(MESHES(NM)%VENTS(NV)%SURF_INDEX)
                  IF (ABS(SF%VEL)>TWO_EPSILON_EB .OR. ABS(SF%VOLUME_FLOW)>TWO_EPSILON_EB .OR. &
                      ABS(SF%MASS_FLUX_TOTAL)>TWO_EPSILON_EB .OR. SF%PYROLYSIS_MODEL/= PYROLYSIS_NONE) THEN
                      WRITE(MESSAGE,'(A,A)') 'ERROR(543):Cannot leak and specify flow or pyrolysis at the same time.  VENT ID:',&
                                             TRIM(MESHES(NM)%VENTS(NV)%ID)
                      CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.)
                  ENDIF
                  IF (ANY(SF%LEAK_PATH>0)) THEN
                      WRITE(MESSAGE,'(A,A,A)') 'ERROR(544):Cannot specify custom leakage and zone leakage with the same surface. ',&
                         ' VENT ID:',TRIM(MESHES(NM)%VENTS(NV)%ID)
                      CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.)
                  ENDIF
               ENDIF
               MESHES(NM)%VENTS(NV)%NODE_INDEX=NN
               ! Sets node to VENT center based on XB. This value will be used for Smokeview visualization. These values
               ! will be modified as needed during the run to reflect the actual area of the VENT that is visible.
               DN%XYZ(1) = 0.5_EB*(MESHES(NM)%VENTS(NV)%X1_ORIG+MESHES(NM)%VENTS(NV)%X2_ORIG)
               DN%XYZ(2) = 0.5_EB*(MESHES(NM)%VENTS(NV)%Y1_ORIG+MESHES(NM)%VENTS(NV)%Y2_ORIG)
               DN%XYZ(3) = 0.5_EB*(MESHES(NM)%VENTS(NV)%Z1_ORIG+MESHES(NM)%VENTS(NV)%Z2_ORIG)
               EXIT NODE_VENT_LOOP
            ENDIF
         ENDDO NODE_VENT_LOOP
         
         IF (.NOT. FOUND) DN%XYZ = -1.E11_EB
      ENDDO MESH_LOOP
      
            ! Check if any MPI process has FOUND the VENT

      IF (N_MPI_PROCESSES>1) CALL MPI_ALLREDUCE(MPI_IN_PLACE,STOP_STATUS,INTEGER_ONE,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERR)

      IF (STOP_STATUS/=0) RETURN

      CALL MPI_ALLREDUCE(MPI_IN_PLACE,FOUND,INTEGER_ONE,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,IERR)
      IF (.NOT. FOUND) THEN
         WRITE(MESSAGE,'(A,A,A,A)') 'ERROR(545): Cannot find VENT_ID: ',TRIM(DN%VENT_ID),' for Ductnode: ',TRIM(DN%ID)
         CALL SHUTDOWN(MESSAGE); RETURN
      ENDIF
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,DN%XYZ(1),INTEGER_ONE,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERR)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,DN%XYZ(2),INTEGER_ONE,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERR)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,DN%XYZ(3),INTEGER_ONE,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERR)
   ENDIF

ENDDO NODE_LOOP

GEOM_IF: IF (N_GEOMETRY > 0 .AND. ANY(SURFACE%NODE_INDEX>0)) THEN
   CF_AREA = 0._EB
   CF_X = 0._EB
   CF_Y = 0._EB
   CF_Z = 0._EB
   NODE_GEOM_LOOP: DO NG=1,N_GEOMETRY
      G=>GEOMETRY(NG)
      IF (.NOT. G%HAVE_NODE) CYCLE
      FACES_LOOP: DO NF=1,G%N_FACES
         SF=>SURFACE(G%SURFS(NF))
         IF (SF%NODE_INDEX==0) CYCLE FACES_LOOP
         FI(1:3)  =  G%FACES((NF-1)*3+1:(NF-1)*3+3)
         XYZ_F(1,1:3) = G%VERTS((FI(1)-1)*3+1:(FI(1)-1)*3+3)
         XYZ_F(2,1:3) = G%VERTS((FI(2)-1)*3+1:(FI(2)-1)*3+3)
         XYZ_F(3,1:3) = G%VERTS((FI(3)-1)*3+1:(FI(3)-1)*3+3)
         AREA = SQRT(((XYZ_F(2,2)-XYZ_F(1,2))*(XYZ_F(3,3)-XYZ_F(1,3))-(XYZ_F(2,3)-XYZ_F(1,3))*(XYZ_F(3,2)-XYZ_F(1,2)))**2 + &
                     ((XYZ_F(2,3)-XYZ_F(1,3))*(XYZ_F(3,1)-XYZ_F(1,1))-(XYZ_F(2,1)-XYZ_F(1,1))*(XYZ_F(3,3)-XYZ_F(1,3)))**2 + &
                     ((XYZ_F(2,1)-XYZ_F(1,1))*(XYZ_F(3,2)-XYZ_F(3,1))-(XYZ_F(2,2)-XYZ_F(1,2))*(XYZ_F(3,1)-XYZ_F(1,1)))**2) * 0.5_EB
         X = SUM(XYZ_F(:,1))*ONTH
         Y = SUM(XYZ_F(:,2))*ONTH
         Z = SUM(XYZ_F(:,3))*ONTH
         CF_AREA(SF%NODE_INDEX) = CF_AREA(SF%NODE_INDEX) + AREA
         CF_X(SF%NODE_INDEX) = CF_X(SF%NODE_INDEX) + X*AREA
         CF_Y(SF%NODE_INDEX) = CF_Y(SF%NODE_INDEX) + Y*AREA
         CF_Z(SF%NODE_INDEX) = CF_Z(SF%NODE_INDEX) + Z*AREA
      ENDDO FACES_LOOP
   ENDDO NODE_GEOM_LOOP 
ENDIF GEOM_IF

NODE_LOOP_2: DO NN=1,N_DUCTNODES
   DN=> DUCTNODE(NN)
   IF (DN%LEAKAGE) CYCLE NODE_LOOP_2
   IF (DN%N_DUCTS==1 .AND. .NOT. (DN%VENT .OR. DN%AMBIENT .OR. DN%GEOM) ) THEN
      WRITE(MESSAGE,'(A,I5,A,A)') 'ERROR(547): Internal ductnode must have at least two attached ducts. Ductnode:',NN,&
                                  ', Ductnode ID:',TRIM(DN%ID)
      CALL SHUTDOWN(MESSAGE); RETURN
   ENDIF
   IF (DN%N_DUCTS> 1 .AND. (DN%AMBIENT .OR. DN%VENT .OR. DN%GEOM) ) THEN
      WRITE(MESSAGE,'(A,I5,A,A)') 'ERROR(548): External ductnode can only have one attached duct. Ductnode:',NN,&
                                  ', Ductnode ID:',TRIM(DN%ID)
      CALL SHUTDOWN(MESSAGE); RETURN
   ENDIF
   IF ((DN%VENT .OR. DN%GEOM) .AND. DN%AMBIENT) THEN
      WRITE(MESSAGE,'(A,I5,A,A)') 'ERROR(546): Ductnode cannot be AMBIENT and have an assigned VENT_ID or GEOM. Ductnode:',NN,&
                                  ', Ductnode ID:',TRIM(DN%ID)
      CALL SHUTDOWN(MESSAGE); RETURN
   ENDIF
   
   IF (DN%GEOM) THEN
      IF (CF_AREA(NN) < TWO_EPSILON_EB) THEN
         WRITE(MESSAGE,'(A,I5,A,A,A)') 'ERROR(573): Ductnode:',NN,', Ductnode ID:',TRIM(DN%ID),&
                                     ' defined with GEOM had no CFACE found.'
         CALL SHUTDOWN(MESSAGE); RETURN
      ENDIF
      DN%XYZ(1) = CF_X(NN)/CF_AREA(NN)
      DN%XYZ(2) = CF_Y(NN)/CF_AREA(NN)
      DN%XYZ(3) = CF_Z(NN)/CF_AREA(NN)
   ENDIF
   
   ALLOCATE(DN%DUCT_INDEX(DN%N_DUCTS))
   ALLOCATE(DN%DIR(DN%N_DUCTS))
   DN%DUCT_INDEX = -1
   DO ND = 1,DN%N_DUCTS
      DO ND2 = 1, N_DUCTS
         IF (NODE_DUCT_A(NN,ND) == DUCT(ND2)%ID) THEN
            DN%DUCT_INDEX(ND) = ND2
            IF (DUCT(ND2)%NODE_INDEX(1)==NN) THEN
               DN%DIR(ND) = -1._EB
            ELSE
               DN%DIR(ND) = 1._EB
            ENDIF
            CYCLE
         ENDIF
      ENDDO
      IF (DN%DUCT_INDEX(ND)==-1) THEN
         WRITE(MESSAGE,'(A,I5,A,I5,A,A)') 'ERROR(549): DUCT ',ND,' not found. Ductnode:',NN,', Ductnode ID:',TRIM(DN%ID)
         CALL SHUTDOWN(MESSAGE); RETURN
      ENDIF
      IF (TRIM(DN%ID) /= TRIM(DUCT_NODE_A(DN%DUCT_INDEX(ND),1)) .AND. TRIM(DN%ID) /= TRIM(DUCT_NODE_A(DN%DUCT_INDEX(ND),2))) THEN
         WRITE(MESSAGE,'(A,A,A,A,A)') 'ERROR(550): Duct: ',TRIM(DUCT(DN%DUCT_INDEX(ND))%ID),' does not contain Ductnode:',&
                                       TRIM(DN%ID), 'in its list of ductnodes.'
         CALL SHUTDOWN(MESSAGE); RETURN
      ENDIF
   ENDDO

   IF(TRIM(NODE_FILTER_A(NN))/='null') THEN
      DO N = 1,N_FILTERS
         IF(TRIM(NODE_FILTER_A(NN))==TRIM(FILTER(N)%ID)) THEN
            DN%FILTER_INDEX = N
            EXIT
         ENDIF
         IF (N==N_FILTERS) THEN
            WRITE(MESSAGE,'(A,A,A,A,A)') 'ERROR(551): Problem with ductnode:',TRIM(DN%ID), &
                                         ',FILTER ',TRIM(NODE_FILTER_A(NN)),' not found'
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
      ENDDO
   ENDIF

   IF (DN%FILTER_INDEX > 0) THEN
      ALLOCATE(DN%FILTER_LOADING(1:N_TRACKED_SPECIES,3))
      DN%FILTER_LOADING(1:N_TRACKED_SPECIES,1) = FILTER(DN%FILTER_INDEX)%INITIAL_LOADING(1:N_TRACKED_SPECIES)
      DN%FILTER_LOADING(1:N_TRACKED_SPECIES,2) = FILTER(DN%FILTER_INDEX)%INITIAL_LOADING(1:N_TRACKED_SPECIES)
      DN%FILTER_LOADING(1:N_TRACKED_SPECIES,3) = FILTER(DN%FILTER_INDEX)%INITIAL_LOADING(1:N_TRACKED_SPECIES)
   ENDIF

ENDDO NODE_LOOP_2

!Temp arrays for input processing
IF (ALLOCATED(DUCT_NODE_A)) DEALLOCATE(DUCT_NODE_A)
IF (ALLOCATED(NODE_DUCT_A)) DEALLOCATE(NODE_DUCT_A)
IF (ALLOCATED(NODE_FILTER_A)) DEALLOCATE(NODE_FILTER_A)
IF (ALLOCATED(DUCT_FAN_A)) DEALLOCATE(DUCT_FAN_A)

ALLOCATE(NODE_PROPERTIES(N_DUCTNODES,8+N_TRACKED_SPECIES)) ; NODE_PROPERTIES = 0._EB
NODE_P    => NODE_PROPERTIES(1:N_DUCTNODES,1)
NODE_TMP  => NODE_PROPERTIES(1:N_DUCTNODES,2)
NODE_RHO  => NODE_PROPERTIES(1:N_DUCTNODES,3)
NODE_H    => NODE_PROPERTIES(1:N_DUCTNODES,4)
NODE_X    => NODE_PROPERTIES(1:N_DUCTNODES,5)
NODE_Y    => NODE_PROPERTIES(1:N_DUCTNODES,6)
NODE_Z    => NODE_PROPERTIES(1:N_DUCTNODES,7)
NODE_AREA => NODE_PROPERTIES(1:N_DUCTNODES,8)
NODE_ZZ   => NODE_PROPERTIES(1:N_DUCTNODES,9:8+N_TRACKED_SPECIES)

ALLOCATE(NODE_ZONE(1:N_DUCTNODES)) ; NODE_ZONE = 0

ALLOCATE(NODE_TMP_EX(1:N_DUCTNODES)) ; NODE_TMP_EX=0._EB
ALLOCATE(NODE_AREA_EX(1:N_DUCTNODES)) ; NODE_AREA_EX = 0._EB
ALLOCATE(NODE_ZZ_EX(1:N_DUCTNODES,1:N_TRACKED_SPECIES)) ; NODE_ZZ_EX = 0._EB

CALL EXAMINE_LOSSES
CALL DETERMINE_FIXED_ELEMENTS(0._EB,DUMMY)
DUCT%AREA_OLD = DUCT%AREA

ALLOCATE(DUCT_MF(1:N_DUCTS))
DUCT_MF = 0._EB

CALL CONNECTION_CHECK

T_USED(13)=T_USED(13)+CURRENT_TIME()-TNOW

RETURN

END SUBROUTINE PROC_HVAC


!< \brief Finds all connected ducts and nodes and hence number of ductruns
SUBROUTINE CONNECTION_CHECK
INTEGER :: NN,NN3,ND
LOGICAL :: NODE_CHECKED(N_DUCTNODES),NODE_CONNECTED(N_DUCTNODES)
TYPE(DUCT_TYPE), POINTER :: DU
TYPE(DUCTNODE_TYPE), POINTER :: DN

NODE_CONNECTED=.FALSE.
NODE_CHECKED=.FALSE.

! Set leakage paths as already checked
DO NN = 1, N_DUCTNODES
   DN=>DUCTNODE(NN)
   IF (DN%LEAKAGE .OR. DUCT(DN%DUCT_INDEX(1))%LOCALIZED_LEAKAGE) NODE_CHECKED(NN) = .TRUE.
END DO

N_CONNECTIVITY_INDICES = 1
NN = 1
L1:DO WHILE (NN <= N_DUCTNODES)
   IF (NODE_CHECKED(NN)) THEN
      NN = NN + 1
      CYCLE L1
   END IF
   DN=>DUCTNODE(NN)
   DL: DO ND = 1, DN%N_DUCTS
      DU=>DUCT(DN%DUCT_INDEX(ND))
      DUCT(DN%DUCT_INDEX(ND))%CONNECTIVITY_INDEX = N_CONNECTIVITY_INDICES
      IF (DU%NODE_INDEX(1)==NN) THEN
         DUCTNODE(DU%NODE_INDEX(2))%CONNECTIVITY_INDEX = N_CONNECTIVITY_INDICES
         NODE_CONNECTED(DU%NODE_INDEX(2)) = .TRUE.
         IF (DUCTNODE(DU%NODE_INDEX(2))%N_DUCTS==1) NODE_CHECKED(DU%NODE_INDEX(2)) = .TRUE.
      ELSE
         DUCTNODE(DU%NODE_INDEX(1))%CONNECTIVITY_INDEX = N_CONNECTIVITY_INDICES
         NODE_CONNECTED(DU%NODE_INDEX(1)) = .TRUE.
         IF (DUCTNODE(DU%NODE_INDEX(1))%N_DUCTS==1) NODE_CHECKED(DU%NODE_INDEX(1)) = .TRUE.
      ENDIF
   END DO DL

   NODE_CHECKED(NN) = .TRUE.
   DUCTNODE(NN)%CONNECTIVITY_INDEX = N_CONNECTIVITY_INDICES
   NODE_CONNECTED(NN) = .TRUE.

   IF (ALL(NODE_CHECKED)) EXIT L1
   L2:DO NN3 = 1, N_DUCTNODES
      IF (NODE_CHECKED(NN3)) THEN
         CYCLE L2
      ELSE IF ((.NOT. NODE_CHECKED(NN3)) .AND. NODE_CONNECTED(NN3)) THEN
         ! If node attached to current node hasn't been checked yet, set the active node to that node
         NN = NN3
         EXIT L2
      ELSE IF (NN3 == N_DUCTNODES) THEN
         ! If there are no nodes part of the current network that haven't been checked, up the network and move to first unchecked
         N_CONNECTIVITY_INDICES = N_CONNECTIVITY_INDICES + 1
         NODE_CONNECTED = .FALSE.
         NN = FINDLOC(NODE_CHECKED,.FALSE.,1)
         EXIT L2
      END IF
   END DO L2
END DO L1

IF (ANY(DUCT%LOCALIZED_LEAKAGE)) THEN
   N_CONNECTIVITY_INDICES = N_CONNECTIVITY_INDICES + 1
   DO ND=1,N_DUCTS
      DU=>DUCT(ND)
      IF (.NOT. DU%LOCALIZED_LEAKAGE) CYCLE
      DU%CONNECTIVITY_INDEX = N_CONNECTIVITY_INDICES
      DUCTNODE(DU%NODE_INDEX(1))%CONNECTIVITY_INDEX = N_CONNECTIVITY_INDICES
      DUCTNODE(DU%NODE_INDEX(2))%CONNECTIVITY_INDEX = N_CONNECTIVITY_INDICES
      IF (DUCTNODE(DU%NODE_INDEX(2))%AMBIENT) DUCTNODE(DU%NODE_INDEX(2))%XYZ = DUCTNODE(DU%NODE_INDEX(1))%XYZ
   ENDDO
ENDIF

END SUBROUTINE CONNECTION_CHECK


!> \brief Initializes duct node properties for non-leakage ducts

SUBROUTINE INIT_DUCT_NODE
USE PHYSICAL_FUNCTIONS, ONLY: GET_ENTHALPY
USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
REAL(EB) :: ZZ_GET(1:N_TRACKED_SPECIES) !< Species mass fraction array.
REAL(EB) :: VF(1:N_TRACKED_SPECIES) !< Species volunme fraction array.
REAL(EB) :: TNOW !< Current CPU time (s) used in computing length of time spent in HVAC routines.
INTEGER :: NN,N
TYPE(DUCTNODE_TYPE), POINTER :: DN
TYPE(INITIALIZATION_TYPE), POINTER :: IN

TNOW = CURRENT_TIME()

DUCT_NODE_INIT = .TRUE.

NODE_LOOP: DO NN = 1, N_DUCTNODES
   DN=> DUCTNODE(NN)
   IF (DN%LEAKAGE) CYCLE NODE_LOOP

   IF (STRATIFICATION .AND. DN%XYZ(3) > -1.E9_EB) THEN
      DN%TMP = TMPA*EVALUATE_RAMP(DN%XYZ(3),I_RAMP_TMP0_Z)
      DN%P   = EVALUATE_RAMP(DN%XYZ(3),I_RAMP_P0_Z)
      DN%RHO = DN%P/(DN%TMP*RSUM0)
   ELSE
      DN%TMP = TMPA
      DN%P   = P_INF
      DN%RHO = RHOA
   ENDIF
   DN%P_OLD = DN%P
   IF (DN%VENT .OR. DN%GEOM) DN%P = -1.E10_EB
   ZZ_GET(1:N_TRACKED_SPECIES) = DN%ZZ_V(1:N_TRACKED_SPECIES)
   CALL GET_ENTHALPY(ZZ_GET,DN%CP,DN%TMP)
   DN%CP     = DN%CP / DN%TMP
   DN%CP_V   = DN%CP
   DN%TMP_V  = DN%TMP
   DN%RSUM_V = DN%RSUM
   DN%RHO_V  = DN%RHO
   DN%ZZ0    = DN%ZZ
   DN%TMP0   = DN%TMP

   INIT_LOOP: DO N = 1,N_INIT
      IF (TRIM(INITIALIZATION(N)%NODE_ID)==TRIM(DN%ID)) THEN
         IN => INITIALIZATION(N)
         IF (IN%VOLUME_FRACTIONS_SPECIFIED) THEN
            VF(2:N_TRACKED_SPECIES) = IN%VOLUME_FRACTION(2:N_TRACKED_SPECIES)
            VF(1) = 1._EB - SUM(VF)
            DN%ZZ0 = VF(1:N_TRACKED_SPECIES)*SPECIES_MIXTURE(1:N_TRACKED_SPECIES)%MW / &
                                          SUM(VF(1:N_TRACKED_SPECIES)*SPECIES_MIXTURE(1:N_TRACKED_SPECIES)%MW)
         ELSEIF (IN%MASS_FRACTIONS_SPECIFIED) THEN
            DN%ZZ0(2:N_TRACKED_SPECIES) = IN%MASS_FRACTION(2:N_TRACKED_SPECIES)
            DN%ZZ0(1) = 1._EB - SUM(DN%ZZ0(2:N_TRACKED_SPECIES))
         ENDIF       
         IF (IN%TEMPERATURE > 0._EB) THEN
            DN%TMP0 = IN%TEMPERATURE
         ENDIF
         EXIT INIT_LOOP
      ENDIF
   ENDDO INIT_LOOP

ENDDO NODE_LOOP

T_USED(13)=T_USED(13)+CURRENT_TIME()-TNOW

END SUBROUTINE INIT_DUCT_NODE

!> \brief Updates the HVAC calculation for a timestep
!>
!> \param T Current time (s)
!> \param DT Current time step (s)
!> \param FIRST_PASS Flag for the first call to the HVAC solver during a timestep, i.e. the first pressure iteration.

SUBROUTINE HVAC_CALC(T,DT,FIRST_PASS)
INTEGER :: NNE,NN,ND
REAL(EB), INTENT(IN) :: T,DT
REAL(EB) :: TNOW !< Current CPU time (s) used in computing length of time spent in HVAC routines.
LOGICAL, INTENT(IN):: FIRST_PASS
LOGICAL :: CHANGE=.TRUE.
LOGICAL, SAVE :: INITIALIZED_HVAC_MASS_TRANSPORT
TYPE(NETWORK_TYPE), POINTER:: NE

TNOW = CURRENT_TIME()

DT_HV = DT
DT_MT = DT

IF (CORRECTOR) THEN
   DUCT%VEL(OLD) = DUCT%VEL(NEW)
   DUCT%DP_FAN(OLD) = DUCT%DP_FAN(NEW)
   DO NN=1,N_DUCTNODES
      IF(DUCTNODE(NN)%FILTER_INDEX > 0) DUCTNODE(NN)%FILTER_LOADING(:,OLD)=DUCTNODE(NN)%FILTER_LOADING(:,NEW)
   ENDDO
   RETURN
ENDIF

FIRST_PASS_IF: IF (FIRST_PASS) THEN
   CALL COLLAPSE_HVAC_BC(T)
   IF (.NOT. INITIALIZED_HVAC_MASS_TRANSPORT .AND. HVAC_MASS_TRANSPORT) CALL SET_INIT_HVAC_MASS_TRANSPORT
   INITIALIZED_HVAC_MASS_TRANSPORT=.TRUE.
   DUCT%CP_D_OLD = DUCT%CP_D
   DUCT%RHO_D_OLD = DUCT%RHO_D
   DUCT%RSUM_D_OLD = DUCT%RSUM_D
   DUCT%DP_FAN(OLD) = DUCT%DP_FAN(NEW)
   WHERE (.NOT. DUCT%FIXED)
      DUCT%VEL(OLD) = DUCT%VEL(NEW)
      DUCT%VEL(PREVIOUS) = DUCT%VEL(NEW)
      DUCT%VEL(GUESS) = DUCT%VEL(NEW)
   END WHERE
   DO ND=1,N_DUCTS
      DUCT(ND)%ZZ_OLD = DUCT(ND)%ZZ
      IF (HVAC_MASS_TRANSPORT .AND. DUCT(ND)%N_CELLS >0) THEN
         DUCT(ND)%CP_C_OLD = DUCT(ND)%CP_C
         DUCT(ND)%RHO_C_OLD = DUCT(ND)%RHO_C
         DUCT(ND)%TMP_C_OLD = DUCT(ND)%TMP_C
         DUCT(ND)%ZZ_C_OLD = DUCT(ND)%ZZ_C
      ENDIF
   ENDDO
   DUCTNODE(:)%RHO_OLD = DUCTNODE(:)%RHO
   DUCTNODE(:)%CP_OLD = DUCTNODE(:)%CP
   DUCTNODE(:)%P_OLD = DUCTNODE(:)%P
   DUCTNODE(:)%TMP_OLD = DUCTNODE(:)%TMP
   DUCTNODE(:)%RSUM_OLD = DUCTNODE(:)%RSUM
   DO NN=1,N_DUCTNODES
      DUCTNODE(NN)%ZZ_OLD = DUCTNODE(NN)%ZZ
   ENDDO
   CALL ADJUST_LEAKAGE_AREA
   CALL FIND_NETWORKS(CHANGE,T) ! calls determined fixed elements (which calls update fan for fixed fans)
   IF (HVAC_QFAN .OR. HVAC_MASS_TRANSPORT) CALL FIND_DUCTRUNS(CHANGE)
ELSE FIRST_PASS_IF !Not FIRST_PASS, reset variables to old solution
   WHERE (.NOT. DUCT%FIXED)
      DUCT(:)%VEL(NEW) = DUCT(:)%VEL(OLD)
      DUCT(:)%VEL(GUESS) = DUCT(:)%VEL(OLD)
      DUCT(:)%VEL(PREVIOUS) = DUCT(:)%VEL(OLD)
   END WHERE
   DUCT(:)%CP_D = DUCT(:)%CP_D_OLD
   DUCT(:)%RHO_D = DUCT(:)%RHO_D_OLD
   DUCT(:)%RHO_D = DUCT(:)%RHO_D_OLD
   DUCT(:)%RSUM_D = DUCT(:)%RSUM_D_OLD
   DO ND=1,N_DUCTS
      DUCT(ND)%ZZ = DUCT(ND)%ZZ_OLD
   ENDDO
   DUCTNODE(:)%RHO = DUCTNODE(:)%RHO_OLD
   DO NN=1,N_DUCTNODES
      DUCTNODE(NN)%ZZ_OLD = DUCTNODE(NN)%ZZ
   ENDDO
   DUCTNODE(:)%CP = DUCTNODE(:)%CP_OLD
   DUCTNODE(:)%P = DUCTNODE(:)%P_OLD
   DUCTNODE(:)%TMP = DUCTNODE(:)%TMP_OLD
   DUCTNODE(:)%RSUM = DUCTNODE(:)%RSUM_OLD
ENDIF FIRST_PASS_IF

ITER = 0

IF (N_ZONE>0) THEN
   ALLOCATE(DPSTAR(1:N_ZONE))
   DPSTAR = 0._EB
   CALL DPSTARCALC
ENDIF

DUCTNODE%P = DUCTNODE%P - P_INF

IF (HVAC_QFAN)  CALL HVAC_QFAN_CALC(T)

! Solution of the HVAC network
NETWORK_LOOP: DO NNE = 1, N_NETWORKS
   NE =>NETWORK(NNE)
   MATRIX_SIZE: IF (NE%N_MATRIX > 0) THEN
      ITER = 0
      ! Reset mass transport for a new iteration
      DUCTNODE%HMT_FILTER = .FALSE.
      ALLOCATE(LHS(NE%N_MATRIX,NE%N_MATRIX))
      ALLOCATE(RHS(NE%N_MATRIX))
      DO WHILE (ITER < ITER_MAX)
         IF(ALLOCATED(DUCTRUN)) DUCTRUN%DT_CFL = DT
         IF (HVAC_MASS_TRANSPORT) THEN
            DO ND=1,N_DUCTS
               IF (DUCT(ND)%N_CELLS <=0) CYCLE
               DUCT(ND)%CP_C = DUCT(ND)%CP_C_OLD
               DUCT(ND)%RHO_C = DUCT(ND)%RHO_C_OLD
               DUCT(ND)%TMP_C = DUCT(ND)%TMP_C_OLD
               DUCT(ND)%ZZ_C = DUCT(ND)%ZZ_C_OLD
            ENDDO
         ENDIF
         LHS = 0._EB
         RHS = 0._EB
         CALL SET_GUESS(NNE,T)
         CALL SET_DONOR(NNE)
         CALL UPDATE_LOSS(T,DT,NNE)
         IF (N_AIRCOILS > 0) CALL COIL_UPDATE(T)
         CALL RHSNODE(NNE)
         CALL RHSDUCT(NNE)
         CALL LHSNODE(NNE)
         CALL LHSDUCT(NNE)
         CALL MATRIX_SOLVE(NNE)
         CALL HVAC_UPDATE(NNE,DT)
         CALL CONVERGENCE_CHECK(NNE)
         ITER = ITER + 1
      ENDDO
      DEALLOCATE(LHS)
      DEALLOCATE(RHS)
   ELSE MATRIX_SIZE
      ! Reset mass transport for a new iteration
      IF (HVAC_MASS_TRANSPORT) THEN
         DO ND=1,N_DUCTS
            IF (DUCT(ND)%N_CELLS <=0) CYCLE
            DUCT(ND)%CP_C = DUCT(ND)%CP_C_OLD
            DUCT(ND)%RHO_C = DUCT(ND)%RHO_C_OLD
            DUCT(ND)%TMP_C = DUCT(ND)%TMP_C_OLD
            DUCT(ND)%ZZ_C = DUCT(ND)%ZZ_C_OLD
         ENDDO
         DUCTRUN%DT_CFL = DT
      ENDIF
      CALL SET_GUESS(NNE,T)
      CALL SET_DONOR(NNE)
      IF (N_AIRCOILS > 0) CALL COIL_UPDATE(T)
      CALL HVAC_UPDATE(NNE,DT)
   ENDIF MATRIX_SIZE
ENDDO NETWORK_LOOP

DUCTNODE%P = DUCTNODE%P + P_INF

!CALL UPDATE_NODE_BC

IF (ALLOCATED(DPSTAR)) DEALLOCATE(DPSTAR)

T_USED(13)=T_USED(13)+CURRENT_TIME()-TNOW

END SUBROUTINE HVAC_CALC

!> \brief Solves the HVAC matrix and extracts the solutions for duct velocity and node pressure.
!>
!> \param NNE Index indicating which HVAC network is being solved

SUBROUTINE MATRIX_SOLVE(NNE)
USE MATH_FUNCTIONS,ONLY : GAUSSJ
INTEGER :: NNE,IERR,ND,NN
TYPE(NETWORK_TYPE), POINTER :: NE
TYPE(DUCT_TYPE), POINTER :: DU
TYPE(DUCTNODE_TYPE), POINTER :: DN

NE =>NETWORK(NNE)
CALL GAUSSJ(LHS,NE%N_MATRIX,NE%N_MATRIX,RHS,1,1,IERR)
DO ND = 1,NE%N_DUCTS
   DU=>DUCT(NE%DUCT_INDEX(ND))
   IF (DU%FIXED .OR. DU%AREA < TWO_EPSILON_EB) CYCLE
   DU%VEL(NEW) = RHS(NE%MATRIX_INDEX(ND))
ENDDO
DO NN = 1,NE%N_DUCTNODES
   DN=>DUCTNODE(NE%NODE_INDEX(NN))
   IF (DN%VENT .OR. DN%AMBIENT .OR. DN%LEAKAGE .OR. DN%GEOM) CYCLE
   DN%P = RHS(NE%MATRIX_INDEX(NE%N_DUCTS+NN))
ENDDO

END SUBROUTINE MATRIX_SOLVE

!> \brief Iterates over the HVAC network updating node and duct quantities
!> \details The routine loops over all ducts and nodes repeatedly. If a node has all of its upstream ducts updated
!> or is a VENT inlet, then the node is updated along with its downstream ducts. The process is repeated until all ducts and
!> nodes have been updated. Following this, ducts with 1D mass transport are then updated.
!>
!> \param DT Current time step (s)
!> \param NNE Index indicating which HVAC network is being solved

SUBROUTINE HVAC_UPDATE(NNE,DT)

USE PHYSICAL_FUNCTIONS, ONLY : GET_SPECIFIC_GAS_CONSTANT,GET_ENTHALPY,GET_TEMPERATURE
REAL(EB),INTENT(IN) :: DT
INTEGER, INTENT(IN) :: NNE
REAL(EB) :: CPTSUM,DU_DX,ETOT,MFLOW,MSUM,MTOT,TGUESS,VFLOW,ZZ_GET(1:N_TRACKED_SPECIES),&
            ZZSUM(1:N_TRACKED_SPECIES),ZZTOT(1:N_TRACKED_SPECIES),HGAS,ESUM_FILTER(1:N_TRACKED_SPECIES)
INTEGER :: NN,ND,NC,NR,NSS,ITMP,N_SUBSTEPS
LOGICAL :: CYCLE_FLAG,SUB_CYCLE_FLAG,D_UPDATED(N_DUCTS),N_UPDATED(N_DUCTNODES),DR_NETWORK(N_DUCTRUNS)
TYPE (DUCTNODE_TYPE), POINTER :: DN
TYPE (DUCT_TYPE), POINTER :: DU
TYPE (DUCTRUN_TYPE), POINTER :: DR
TYPE (NETWORK_TYPE), POINTER :: NE

DR_NETWORK = .FALSE.
NE => NETWORK(NNE)

DO NN = 1,NE%N_DUCTNODES
   DUCTNODE(NE%NODE_INDEX(NN))%ZZ_SUM=0._EB
   DUCTNODE(NE%NODE_INDEX(NN))%E_SUM=0._EB
   NODE_TMP_EX(NE%NODE_INDEX(NN)) = DUCTNODE(NE%NODE_INDEX(NN))%TMP
   NODE_ZZ_EX(NE%NODE_INDEX(NN),1:N_TRACKED_SPECIES) = DUCTNODE(NE%NODE_INDEX(NN))%ZZ(1:N_TRACKED_SPECIES)
ENDDO

HMT_IF: IF (HVAC_MASS_TRANSPORT) THEN
   ! Outputs number of substeps required to maintain CFL for mass transport
   ! Relevant CFL is all mass leaving a duct within one DT.
   DO ND = 1,NE%N_DUCTS
      DU=>DUCT(NE%DUCT_INDEX(ND))
      IF (DU%N_CELLS > 0 .AND. ABS(DU%VEL(NEW)) > 0._EB) THEN
         DR_NETWORK(DU%DUCTRUN) = .TRUE.
         NR = DU%DUCTRUN
         DUCTRUN(NR)%DT_CFL = MIN(DU%LENGTH/ABS(DU%VEL(NEW)),DUCTRUN(NR)%DT_CFL)
      ENDIF
   ENDDO

   ! Update ductruns with mass transport
   UPDATE_DUCTRUNS: DO NR = 1,N_DUCTRUNS
      IF (.NOT. DR_NETWORK(NR)) CYCLE UPDATE_DUCTRUNS
      DR=>DUCTRUN(NR)
      N_SUBSTEPS = CEILING(DT/DR%DT_CFL)
      DT_MT = DT/REAL(N_SUBSTEPS,EB)
      SUB_CYCLE_FLAG = .TRUE.

      SUBSTEP_LOOP: DO NSS = 1, N_SUBSTEPS
         D_UPDATED = .FALSE.
         N_UPDATED = .FALSE.
         DO NN = 1,DR%N_DUCTNODES
            DN=>DUCTNODE(DR%NODE_INDEX(NN))
            IF ((DN%VENT .OR. DN%AMBIENT .OR. DN%GEOM) .AND. &
               DN%DIR(1)*DUCT(DN%DUCT_INDEX(1))%VEL(NEW) <= 0._EB) N_UPDATED(DR%NODE_INDEX(NN)) = .TRUE.
            IF (DN%N_DUCTS==2) THEN
               IF (ABS(DUCT(DN%DUCT_INDEX(1))%VEL(NEW)) < TWO_EPSILON_EB .AND. &
                   ABS(DUCT(DN%DUCT_INDEX(2))%VEL(NEW)) < TWO_EPSILON_EB) N_UPDATED(DR%NODE_INDEX(NN)) = .TRUE.
            ENDIF
         ENDDO

         ITER_LOOP_1: DO
            CYCLE_FLAG = .FALSE.
            DUCT_LOOP_1:DO ND = 1,DR%N_DUCTS
               DU=>DUCT(DR%DUCT_INDEX(ND))
               IF (D_UPDATED(DR%DUCT_INDEX(ND))) CYCLE DUCT_LOOP_1
               CYCLE_FLAG = .TRUE.
               IF (DU%VEL(NEW) > TWO_EPSILON_EB) THEN
                  DN => DUCTNODE(DU%NODE_INDEX(1))
                  NN = DU%NODE_INDEX(1)
               ELSEIF (DU%VEL(NEW) < -TWO_EPSILON_EB) THEN
                  DN => DUCTNODE(DU%NODE_INDEX(2))
                  NN = DU%NODE_INDEX(2)
               ELSE
                  DU%VEL(NEW) = 0._EB
                  IF (N_UPDATED(DU%NODE_INDEX(1)) .AND. N_UPDATED(DU%NODE_INDEX(2))) THEN
                     DU%RHO_D = 0.5_EB*(DUCTNODE(DU%NODE_INDEX(1))%RHO+DUCTNODE(DU%NODE_INDEX(2))%RHO)
                     DU%TMP_D = 0.5_EB*(DUCTNODE(DU%NODE_INDEX(1))%TMP+DUCTNODE(DU%NODE_INDEX(2))%TMP)
                     DU%CP_D = 0.5_EB*(DUCTNODE(DU%NODE_INDEX(1))%CP+DUCTNODE(DU%NODE_INDEX(2))%CP)
                     DU%ZZ(:) = 0.5_EB*(DUCTNODE(DU%NODE_INDEX(1))%ZZ(:)+DUCTNODE(DU%NODE_INDEX(2))%ZZ(:))
                     DUCT_MF(DR%DUCT_INDEX(ND)) = 0._EB
                     D_UPDATED(DR%DUCT_INDEX(ND)) = .TRUE.
                  ENDIF
                  CYCLE DUCT_LOOP_1
               ENDIF
               IF (N_UPDATED(NN)) THEN
                  DU%RHO_D  = DN%RHO
                  DU%TMP_D  = DN%TMP
                  DU%CP_D   = DN%CP
                  DU%ZZ(:) = DN%ZZ(:)
                  DUCT_MF(DR%DUCT_INDEX(ND)) = DU%VEL(NEW)*DU%AREA*DU%RHO_D
                  D_UPDATED(DR%DUCT_INDEX(ND)) = .TRUE.
               ENDIF
            ENDDO DUCT_LOOP_1

            NODE_LOOP_1:DO NN = 1,DR%N_DUCTNODES
               DN=>DUCTNODE(DR%NODE_INDEX(NN))
               IF (N_UPDATED(DR%NODE_INDEX(NN))) CYCLE NODE_LOOP_1
               CYCLE_FLAG = .TRUE.
               MTOT = 0._EB
               ETOT = 0._EB
               ZZTOT = 0._EB
               TGUESS = 0._EB
               CPTSUM = 0
               IF (DN%FILTER_INDEX > 0) THEN
                  DN%HMT_FILTER = .TRUE.
                  DN%FILTER_LOADING(:,3) = 0._EB
                  DN%FILTER_LOADING(:,2) = 0._EB
                  ESUM_FILTER = 0._EB
               ENDIF
               DUCT_NODE_LOOP_1: DO ND = 1,DN%N_DUCTS
                  DU => DUCT(DN%DUCT_INDEX(ND))
                  IF (DU%AREA<=TWO_EPSILON_EB) CYCLE
                  IF (DU%VEL(NEW)*DN%DIR(ND) <= 0._EB) CYCLE
                  IF (.NOT. D_UPDATED(DN%DUCT_INDEX(ND)) .AND. ABS(DU%VEL(NEW)) >= 1.E-5_EB) CYCLE NODE_LOOP_1

                  ! Duct is discretized: we need to find the end of the lump of mass advected within a (sub)step
                  MFLOW = ABS(DU%VEL(NEW)*DU%RHO_D)*DT_MT ! Duct mass flow corrected for any substepping
                  MSUM = 0
                  ZZSUM = 0
                  CPTSUM = 0
                  DU_DX = DU%LENGTH/REAL(DU%N_CELLS,EB)
                  ! Sums cumulative mass, species and energy from the final cell, to locate end of advected lump of mass
                  IF (DU%VEL(NEW) > 0._EB) THEN
                     DO NC = DU%N_CELLS,1,-1
                        IF (MSUM + DU%RHO_C(NC)*DU_DX > MFLOW) THEN
                           DU_DX = MAX(0._EB,(MFLOW - MSUM)/DU%RHO_C(NC))
                           ZZSUM(:) = ZZSUM(:) + DU%RHO_C(NC)*DU%ZZ_C(NC,:)*DU_DX
                           CPTSUM = CPTSUM + DU%RHO_C(NC)*DU%TMP_C(NC)*DU%CP_C(NC)*DU_DX
                           IF (DN%FILTER_INDEX>0) THEN
                              ITMP = MIN(I_MAX_TEMP,NINT(DU%TMP_C(NC)))
                              ESUM_FILTER(:) = ESUM_FILTER(:)+DU%RHO_C(NC)*DU%ZZ_C(NC,:)*DU_DX*CPBAR_Z(ITMP,:)*DU%TMP_C(NC)
                           ENDIF
                           EXIT
                        ELSE
                           MSUM = MSUM + DU%RHO_C(NC)*DU_DX
                           ZZSUM(:) = ZZSUM(:) + DU%RHO_C(NC)*DU%ZZ_C(NC,:)*DU_DX
                           CPTSUM = CPTSUM + DU%RHO_C(NC)*DU%TMP_C(NC)*DU%CP_C(NC)*DU_DX
                           IF (DN%FILTER_INDEX>0) THEN
                              ITMP = MIN(I_MAX_TEMP,NINT(DU%TMP_C(NC)))
                              ESUM_FILTER(:) = ESUM_FILTER(:)+DU%RHO_C(NC)*DU%ZZ_C(NC,:)*DU_DX*CPBAR_Z(ITMP,:)*DU%TMP_C(NC)
                           ENDIF
                        ENDIF
                     ENDDO
                  ELSE
                     DO NC = 1,DU%N_CELLS
                        IF (MSUM + DU%RHO_C(NC)*DU_DX > MFLOW) THEN
                           DU_DX = MAX(0._EB,(MFLOW - MSUM)/DU%RHO_C(NC))
                           ZZSUM(:) = ZZSUM(:) + DU%RHO_C(NC)*DU%ZZ_C(NC,:)*DU_DX
                           CPTSUM = CPTSUM + DU%RHO_C(NC)*DU%TMP_C(NC)*DU%CP_C(NC)*DU_DX
                           IF (DN%FILTER_INDEX>0) THEN
                              ITMP = MIN(I_MAX_TEMP,NINT(DU%TMP_C(NC)))
                              ESUM_FILTER(:) = ESUM_FILTER(:)+DU%RHO_C(NC)*DU%ZZ_C(NC,:)*DU_DX*CPBAR_Z(ITMP,:)*DU%TMP_C(NC)
                           ENDIF
                           EXIT
                        ELSE
                           MSUM = MSUM + DU%RHO_C(NC)*DU_DX
                           ZZSUM(:) = ZZSUM(:) + DU%RHO_C(NC)*DU%ZZ_C(NC,:)*DU_DX
                           CPTSUM = CPTSUM + DU%RHO_C(NC)*DU%TMP_C(NC)*DU%CP_C(NC)*DU_DX
                           IF (DN%FILTER_INDEX>0) THEN
                              ITMP = MIN(I_MAX_TEMP,NINT(DU%TMP_C(NC)))
                              ESUM_FILTER(:) = ESUM_FILTER(:)+DU%RHO_C(NC)*DU%ZZ_C(NC,:)*DU_DX*CPBAR_Z(ITMP,:)*DU%TMP_C(NC)
                           ENDIF
                        ENDIF
                     ENDDO
                  ENDIF

                  ! Cumulative sums of energy, mass, species passing through duct node from ducts
                  ETOT = ETOT + CPTSUM * DU%AREA / DT_MT + DU%COIL_Q
                  MTOT = MTOT + MFLOW * DU%AREA / DT_MT
                  ZZTOT = ZZTOT + ZZSUM * DU%AREA / DT_MT
                  TGUESS = TGUESS + MFLOW * DU%AREA / DT_MT * DU%TMP_D

                  IF (DN%FILTER_INDEX > 0) THEN
                     DN%HMT_FILTER = .TRUE.
                     DN%FILTER_LOADING(:,3) = ZZSUM * DU%AREA / DT_MT * FILTER(DN%FILTER_INDEX)%EFFICIENCY
                     DN%FILTER_LOADING(:,2) = DN%FILTER_LOADING(:,1) + DN%FILTER_LOADING(:,3) * DT_MT
                     MTOT = MTOT - SUM(DN%FILTER_LOADING(:,3))
                     ZZTOT = ZZTOT - DN%FILTER_LOADING(:,3)
                     ITMP = MIN(I_MAX_TEMP,NINT(DU%TMP_D))
                     ESUM_FILTER(:) = ESUM_FILTER(:)*DU%AREA/DT_MT
                     ETOT = ETOT - SUM(ESUM_FILTER*FILTER(DN%FILTER_INDEX)%EFFICIENCY)
                  ENDIF
               ENDDO DUCT_NODE_LOOP_1

               IF (DN%VENT .OR. DN%GEOM) THEN
                  DN%E_SUM = DN%E_SUM + ETOT*DT_MT
                  DN%ZZ_SUM = DN%ZZ_SUM + ZZTOT*DT_MT
               ENDIF

               N_UPDATED(DR%NODE_INDEX(NN)) = .TRUE.
               IF (MTOT==0._EB) CYCLE NODE_LOOP_1
               ZZ_GET = 0._EB
               DN%ZZ(:)  = ZZTOT/MTOT
               ZZ_GET(1:N_TRACKED_SPECIES) = DN%ZZ(1:N_TRACKED_SPECIES)
               CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,DN%RSUM)
               ETOT = ETOT/ MTOT
               DN%TMP = TGUESS/MTOT
               CALL GET_TEMPERATURE(DN%TMP,ETOT,ZZ_GET)
               CALL GET_ENTHALPY(ZZ_GET,HGAS,DN%TMP)
               DN%CP = HGAS/DN%TMP
               DN%RHO = (DN%P+P_INF)/(DN%RSUM*DN%TMP)
            ENDDO NODE_LOOP_1

            IF (.NOT. CYCLE_FLAG) EXIT ITER_LOOP_1
         ENDDO ITER_LOOP_1
         CALL UPDATE_HVAC_MASS_TRANSPORT(DT_MT,NR)

      ENDDO SUBSTEP_LOOP

      DO NN = 1,DR%N_DUCTNODES
         DN=>DUCTNODE(DR%NODE_INDEX(NN))
         IF (.NOT. (DN%VENT .OR. DN%GEOM)) CYCLE
         DU => DUCT(DN%DUCT_INDEX(1))
         IF (DN%DIR(1)*DU%VEL(NEW) <= 0._EB) CYCLE
         MTOT = ABS(DU%VEL(NEW))*DU%AREA*DU%RHO_D*DT
         NODE_ZZ_EX(DR%NODE_INDEX(NN),1:N_TRACKED_SPECIES) = DN%ZZ_SUM(1:N_TRACKED_SPECIES)/MTOT
         ZZ_GET(1:N_TRACKED_SPECIES) = NODE_ZZ_EX(DR%NODE_INDEX(NN),1:N_TRACKED_SPECIES)
         DN%E_SUM = DN%E_SUM/MTOT
         NODE_TMP_EX(DR%NODE_INDEX(NN)) = DN%TMP
         CALL GET_TEMPERATURE(NODE_TMP_EX(DR%NODE_INDEX(NN)),DN%E_SUM,ZZ_GET)
      ENDDO

   ENDDO UPDATE_DUCTRUNS
ENDIF HMT_IF

! Update remainder of network with non-mass transport ducts

D_UPDATED = .FALSE.
N_UPDATED = .FALSE.

DO NN = 1,NE%N_DUCTNODES
   DN=>DUCTNODE(NE%NODE_INDEX(NN))
   IF (DN%DUCTRUN > 0) THEN
      IF (DR_NETWORK(DN%DUCTRUN)) THEN
         N_UPDATED(NE%NODE_INDEX(NN)) = .TRUE.
         CYCLE
      ENDIF
   ENDIF
   IF ((DN%LEAKAGE .OR. DN%VENT .OR. DN%AMBIENT .OR. DN%GEOM) .AND. &
        DN%DIR(1)*DUCT(DN%DUCT_INDEX(1))%VEL(NEW) <= 0._EB) N_UPDATED(NE%NODE_INDEX(NN)) = .TRUE.
   IF (DN%N_DUCTS==2) THEN
      IF (ABS(DUCT(DN%DUCT_INDEX(1))%VEL(NEW)) < TWO_EPSILON_EB .AND. &
          ABS(DUCT(DN%DUCT_INDEX(2))%VEL(NEW)) < TWO_EPSILON_EB) N_UPDATED(NE%NODE_INDEX(NN)) = .TRUE.
   ENDIF
ENDDO

ITER_LOOP_2: DO
   CYCLE_FLAG = .FALSE.
   DUCT_LOOP_2:DO ND = 1,NE%N_DUCTS
      DU=>DUCT(NE%DUCT_INDEX(ND))
      IF (D_UPDATED(NE%DUCT_INDEX(ND))) CYCLE DUCT_LOOP_2
      IF (DU%DUCTRUN > 0) THEN
         IF (DR_NETWORK(DU%DUCTRUN)) THEN
            D_UPDATED(NE%DUCT_INDEX(ND)) = .TRUE.
            CYCLE DUCT_LOOP_2
         ENDIF
      ENDIF
      CYCLE_FLAG = .TRUE.
      IF (DU%VEL(NEW) > TWO_EPSILON_EB) THEN
         DN => DUCTNODE(DU%NODE_INDEX(1))
         NN = DU%NODE_INDEX(1)
      ELSEIF (DU%VEL(NEW) < -TWO_EPSILON_EB) THEN
         DN => DUCTNODE(DU%NODE_INDEX(2))
         NN = DU%NODE_INDEX(2)
      ELSE
         DU%VEL(NEW) = 0._EB
         IF (N_UPDATED(DU%NODE_INDEX(1)) .AND. N_UPDATED(DU%NODE_INDEX(2))) THEN
            DU%RHO_D = 0.5_EB*(DUCTNODE(DU%NODE_INDEX(1))%RHO+DUCTNODE(DU%NODE_INDEX(2))%RHO)
            DU%TMP_D = 0.5_EB*(DUCTNODE(DU%NODE_INDEX(1))%TMP+DUCTNODE(DU%NODE_INDEX(2))%TMP)
            DU%CP_D = 0.5_EB*(DUCTNODE(DU%NODE_INDEX(1))%CP+DUCTNODE(DU%NODE_INDEX(2))%CP)
            DU%ZZ(:) = 0.5_EB*(DUCTNODE(DU%NODE_INDEX(1))%ZZ(:)+DUCTNODE(DU%NODE_INDEX(2))%ZZ(:))
            DUCT_MF(NE%DUCT_INDEX(ND)) = 0._EB
            D_UPDATED(NE%DUCT_INDEX(ND)) = .TRUE.
         ENDIF
         CYCLE DUCT_LOOP_2
      ENDIF
      IF (N_UPDATED(NN)) THEN
         DU%RHO_D  = DN%RHO
         DU%TMP_D  = DN%TMP
         DU%CP_D   = DN%CP
         DU%ZZ(:) = DN%ZZ(:)
         DUCT_MF(NE%DUCT_INDEX(ND)) = DU%VEL(NEW)*DU%AREA*DU%RHO_D
         D_UPDATED(NE%DUCT_INDEX(ND)) = .TRUE.
      ENDIF
   ENDDO DUCT_LOOP_2

   NODE_LOOP_2:DO NN = 1,NE%N_DUCTNODES
      DN=>DUCTNODE(NE%NODE_INDEX(NN))
      IF (N_UPDATED(NE%NODE_INDEX(NN))) CYCLE NODE_LOOP_2
      CYCLE_FLAG = .TRUE.
      MTOT = 0._EB
      ETOT = 0._EB
      ZZTOT = 0._EB
      TGUESS = 0._EB
      CPTSUM = 0
      DUCT_NODE_LOOP_2: DO ND = 1,DN%N_DUCTS
         DU => DUCT(DN%DUCT_INDEX(ND))
         IF (DU%AREA<=TWO_EPSILON_EB) CYCLE
         IF (DU%VEL(NEW)*DN%DIR(ND) < 0._EB) CYCLE
         IF (.NOT. D_UPDATED(DN%DUCT_INDEX(ND)) .AND. ABS(DU%VEL(NEW)) >= 1.E-5_EB) CYCLE NODE_LOOP_2
         VFLOW = ABS(DU%VEL(NEW)*DU%AREA)
         MTOT = MTOT + VFLOW * DU%RHO_D
         ETOT = ETOT + VFLOW * DU%RHO_D * DU%TMP_D * DU%CP_D + DU%COIL_Q
         TGUESS = TGUESS + VFLOW * DU%RHO_D * DU%TMP_D
         ZZTOT = ZZTOT + VFLOW * DU%RHO_D * DU%ZZ
         IF (DN%FILTER_INDEX > 0) THEN
            MTOT = MTOT - SUM(DN%FILTER_LOADING(:,3))
            ZZTOT = ZZTOT - DN%FILTER_LOADING(:,3)
            ITMP = MIN(I_MAX_TEMP,NINT(DU%TMP_D))
            ETOT = ETOT - SUM(CPBAR_Z(ITMP,:)*DU%TMP_D*DN%FILTER_LOADING(:,3))
         ENDIF
      ENDDO DUCT_NODE_LOOP_2
      N_UPDATED(NE%NODE_INDEX(NN)) = .TRUE.
      IF (ABS(MTOT)<=TWO_EPSILON_EB) CYCLE NODE_LOOP_2
      ZZ_GET = 0._EB
      DN%ZZ(:)  = ZZTOT/MTOT
      ZZ_GET(1:N_TRACKED_SPECIES) = DN%ZZ(1:N_TRACKED_SPECIES)
      CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,DN%RSUM)
      ETOT = ETOT/ MTOT
      DN%TMP = TGUESS/MTOT
      CALL GET_TEMPERATURE(DN%TMP,ETOT,ZZ_GET)
      CALL GET_ENTHALPY(ZZ_GET,HGAS,DN%TMP)
      DN%CP = HGAS/DN%TMP
      DN%RHO = (DN%P+P_INF)/(DN%RSUM*DN%TMP)
      NODE_TMP_EX(NE%NODE_INDEX(NN)) = DN%TMP
      NODE_ZZ_EX(NE%NODE_INDEX(NN),:) = DN%ZZ(:)

   ENDDO NODE_LOOP_2

   IF (.NOT. CYCLE_FLAG) EXIT ITER_LOOP_2

ENDDO ITER_LOOP_2

END SUBROUTINE HVAC_UPDATE

!> \brief Builds the right hand side of the HVAC flow matrix for mass conservation at internal nodes
!>
!> \param NETWORK_INDEX Index indicating which HVAC network is being solved

SUBROUTINE RHSNODE(NETWORK_INDEX)
USE GLOBAL_CONSTANTS
INTEGER, INTENT(IN)::NETWORK_INDEX
INTEGER :: NN,ND, ARRAYLOC
TYPE(NETWORK_TYPE), POINTER::NE
TYPE(DUCT_TYPE), POINTER::DU
TYPE(DUCTNODE_TYPE), POINTER::DN

NE => NETWORK(NETWORK_INDEX)
DO NN = 1, NE%N_DUCTNODES
   DN => DUCTNODE(NE%NODE_INDEX(NN))
   IF (DN%VENT .OR. DN%AMBIENT .OR. DN%LEAKAGE .OR. DN%GEOM) CYCLE
   ARRAYLOC = NE%MATRIX_INDEX(NE%N_DUCTS+DUCTNODE_NE(NE%NODE_INDEX(NN)))
   DO ND = 1,DN%N_DUCTS
      DU => DUCT(DN%DUCT_INDEX(ND))
      IF (DU%FIXED) RHS(ARRAYLOC) = RHS(ARRAYLOC) + DN%DIR(ND)*DU%RHO_D*DU%VOLUME_FLOW
   END DO
   IF (DN%FILTER_INDEX > 0) RHS(ARRAYLOC) = RHS(ARRAYLOC) - SUM(DN%FILTER_LOADING(:,3))

ENDDO

END SUBROUTINE RHSNODE

!> \brief Builds the left hand side of the HVAC flow matrix for mass conservation at internal nodes
!>
!> \param NETWORK_INDEX Index indicating which HVAC network is being solved

SUBROUTINE LHSNODE(NETWORK_INDEX)
USE GLOBAL_CONSTANTS
INTEGER, INTENT(IN)::NETWORK_INDEX
INTEGER :: NN,ND, ARRAYLOC1,ARRAYLOC2
TYPE(NETWORK_TYPE), POINTER::NE
TYPE(DUCT_TYPE), POINTER::DU
TYPE(DUCTNODE_TYPE), POINTER::DN

NE => NETWORK(NETWORK_INDEX)
DO NN = 1, NE%N_DUCTNODES
   DN => DUCTNODE(NE%NODE_INDEX(NN))
   IF (DN%VENT .OR. DN%LEAKAGE .OR. DN%AMBIENT .OR. DN%GEOM) CYCLE
   ARRAYLOC1 = NE%MATRIX_INDEX(NE%N_DUCTS+DUCTNODE_NE(NE%NODE_INDEX(NN)))
   DO ND = 1,DN%N_DUCTS
      DU => DUCT(DN%DUCT_INDEX(ND))
      IF (DU%FIXED .OR. DU%AREA <=TWO_EPSILON_EB) CYCLE
      ARRAYLOC2 = NE%MATRIX_INDEX(DUCT_NE(DN%DUCT_INDEX(ND)))
      LHS(ARRAYLOC1,ARRAYLOC2) = -DN%DIR(ND)*DU%RHO_D*DU%AREA
   END DO
ENDDO

END SUBROUTINE LHSNODE

!> \brief Computes the extrapolated pressure for each ZONE

SUBROUTINE DPSTARCALC

USE GLOBAL_CONSTANTS
INTEGER :: NN,IPZ,IOPZ
TYPE(DUCT_TYPE), POINTER::DU
TYPE(DUCTNODE_TYPE), POINTER::DN
TYPE(P_ZONE_TYPE), POINTER::PZ

DPSTAR = 0._EB
DO IPZ = 1,N_ZONE
   PZ => P_ZONE(IPZ)
   IF (PZ%N_DUCTNODES==0) CYCLE
   DPSTAR(IPZ) = P_ZONE(IPZ)%DPSTAR * DT_HV
   PSUM_TOT(IPZ) = PSUM(IPZ)
   DO IOPZ = 1,N_ZONE
      IF (IPZ==IOPZ) CYCLE
      IF (CONNECTED_ZONES(IPZ,IOPZ)>0) PSUM_TOT(IPZ) = PSUM_TOT(IPZ) + PSUM(IOPZ)
   ENDDO
ENDDO

DO IPZ = 1,N_ZONE
   PZ => P_ZONE(IPZ)
   IF (PZ%N_DUCTNODES==0) CYCLE
   DO NN = 1,PZ%N_DUCTNODES
      DN=>DUCTNODE(PZ%NODE_INDEX(NN))
      DU=>DUCT(DN%DUCT_INDEX(1))
      DPSTAR(IPZ) = DPSTAR(IPZ)  - DN%DIR(1) * DU%AREA * DU%VEL(OLD) * DT_HV/PSUM_TOT(IPZ)
      IF (DU%FIXED) DPSTAR(IPZ) = DPSTAR(IPZ)  + DN%DIR(1) * DU%AREA * DU%VEL(NEW) * DT_HV/PSUM_TOT(IPZ)
      DO IOPZ = 1, N_ZONE
         IF (IPZ==IOPZ) CYCLE
         IF (P_ZONE(IOPZ)%N_DUCTNODES==0) CYCLE
         IF (CONNECTED_ZONES(IPZ,IOPZ)>0) THEN
            IF (P_ZONE(IOPZ)%N_DUCTNODES==0) CYCLE
            DPSTAR(IOPZ) = DPSTAR(IOPZ)  - DN%DIR(1) * DU%AREA * DU%VEL(OLD) * DT_HV/PSUM_TOT(IPZ)
            IF (DU%FIXED) DPSTAR(IOPZ) = DPSTAR(IOPZ)  + DN%DIR(1) * DU%AREA * DU%VEL(NEW) * DT_HV/PSUM_TOT(IPZ)
         ENDIF
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE DPSTARCALC

!> \brief Builds the right hand side of the HVAC flow matrix for momentum conservation in a duct
!>
!> \param NETWORK_INDEX Index indicating which HVAC network is being solved

SUBROUTINE RHSDUCT(NETWORK_INDEX)

USE GLOBAL_CONSTANTS
INTEGER, INTENT(IN)::NETWORK_INDEX
INTEGER :: ND, ARRAYLOC,IPZ
REAL(EB) :: HEAD !< DUCT pressure head (Pa)
REAL(EB) :: RGZ !< DUCT elevation head (Pa)
REAL(EB) :: XYZ(3) !< Distance between duct nodes (m)
TYPE(NETWORK_TYPE), POINTER::NE
TYPE(DUCT_TYPE), POINTER::DU
TYPE(DUCTNODE_TYPE), POINTER::DN

NE => NETWORK(NETWORK_INDEX)

DO ND = 1, NE%N_DUCTS
   DU => DUCT(NE%DUCT_INDEX(ND))
   IF (DU%FIXED .OR. DU%AREA < TWO_EPSILON_EB) CYCLE
   HEAD = 0._EB
   RGZ = 0._EB
   ARRAYLOC = NE%MATRIX_INDEX(DUCT_NE(NE%DUCT_INDEX(ND)))
   DN=>DUCTNODE(DU%NODE_INDEX(1))
   IF (DN%AMBIENT) THEN
      HEAD = HEAD + DN%P
   ELSEIF (DN%VENT .OR. DN%LEAKAGE .OR. DN%GEOM) THEN
      HEAD = HEAD + DN%P
      IF ((DN%VENT .OR. DN%GEOM) .AND. HVAC_LOCAL_PRESSURE) &
         HEAD = HEAD - DN%DIR(1) * DU%RHO_D * DU%VEL(OLD) * DU%AREA / NODE_AREA_EX(DU%NODE_INDEX(1)) / (GAMMA * DT_HV)
      IF (N_ZONE > 0) THEN
         IPZ = DN%ZONE_INDEX
         IF (IPZ > 0) HEAD = HEAD + DPSTAR(IPZ)
      ENDIF
   ENDIF
   XYZ = DN%XYZ
   DN=>DUCTNODE(DU%NODE_INDEX(2))
   IF (DN%AMBIENT) THEN
      HEAD = HEAD - DN%P
   ELSEIF (DN%VENT .OR. DN%LEAKAGE .OR. DN%GEOM) THEN
      HEAD = HEAD - DN%P
      IF ((DN%VENT .OR. DN%GEOM) .AND. HVAC_LOCAL_PRESSURE) &
         HEAD = HEAD + DN%DIR(1) * DU%RHO_D * DU%VEL(OLD) * DU%AREA / NODE_AREA_EX(DU%NODE_INDEX(2)) / (GAMMA * DT_HV)
      IF (N_ZONE > 0) THEN
         IPZ = DN%ZONE_INDEX
         IF (IPZ > 0) HEAD = HEAD - DPSTAR(IPZ)
      ENDIF
   ENDIF

   XYZ = DN%XYZ - XYZ
   IF (.NOT. DU%LEAKAGE) THEN
      IF (STRATIFICATION) THEN
         RGZ = (GVEC(1)*XYZ(1)+GVEC(2)*XYZ(2)+GVEC(3)*XYZ(3))*DU%RHO_D
      ELSE
         RGZ = (GVEC(1)*XYZ(1)+GVEC(2)*XYZ(2)+GVEC(3)*XYZ(3))* &
               (DUCTNODE(DU%NODE_INDEX(1))%RHO - DUCTNODE(DU%NODE_INDEX(2))%RHO)
      ENDIF
   ENDIF

   RHS(ARRAYLOC) = DU%VEL(OLD)+DT_HV/DU%LENGTH*((HEAD+RGZ+SUM(DU%DP_FAN)*0.5_EB)/DU%RHO_D + &
                   0.5_EB*DU%TOTAL_LOSS*ABS(DU%VEL(PREVIOUS))*DU%VEL(GUESS))
ENDDO

END SUBROUTINE RHSDUCT

!> \brief Builds the left hand side of the HVAC flow matrix for momentum conservation in a duct
!>
!> \param NETWORK_INDEX Index indicating which HVAC network is being solved

SUBROUTINE LHSDUCT(NETWORK_INDEX)
USE GLOBAL_CONSTANTS
INTEGER, INTENT(IN)::NETWORK_INDEX
INTEGER :: NN,ND, ARRAYLOC1,ARRAYLOC2
TYPE(NETWORK_TYPE), POINTER::NE
TYPE(DUCT_TYPE), POINTER::DU,DU2
TYPE(DUCTNODE_TYPE), POINTER::DN,DN2
TYPE(P_ZONE_TYPE), POINTER::PZ

NE => NETWORK(NETWORK_INDEX)

DUCT_LOOP: DO ND = 1, NE%N_DUCTS
   DU => DUCT(NE%DUCT_INDEX(ND))
   IF (DU%FIXED .OR. DU%AREA < TWO_EPSILON_EB) CYCLE DUCT_LOOP
   ARRAYLOC1 = NE%MATRIX_INDEX(DUCT_NE(NE%DUCT_INDEX(ND)))
   LHS(ARRAYLOC1,ARRAYLOC1) = 1._EB+0.5_EB*DT_HV*DU%TOTAL_LOSS/DU%LENGTH*ABS(DU%VEL(PREVIOUS)+DU%VEL(GUESS))
   DN=>DUCTNODE(DU%NODE_INDEX(1))
   IF (.NOT. (DN%VENT .OR. DN%GEOM)) THEN
      IF (.NOT. DN%AMBIENT .AND. .NOT. DN%LEAKAGE) THEN
         ARRAYLOC2 = NE%MATRIX_INDEX(NE%N_DUCTS+DUCTNODE_NE(DU%NODE_INDEX(1)))
         LHS(ARRAYLOC1,ARRAYLOC2) = -DT_HV/(DU%RHO_D*DU%LENGTH)
      ENDIF
      IF (DN%LEAKAGE .AND. .NOT. DN%AMBIENT) THEN
         !IF (HVAC_LOCAL_PRESSURE) &
         !   LHS(ARRAYLOC1,ARRAYLOC1) = LHS(ARRAYLOC1,ARRAYLOC1)-DN%DIR(1)*DU%AREA/(DU%LENGTH*GAMMA*NODE_AREA_EX(DU%NODE_INDEX(1)))
         PZ => P_ZONE(DN%ZONE_INDEX)
         DO NN=1,PZ%N_DUCTNODES
            DN2=>DUCTNODE(PZ%NODE_INDEX(NN))
            DU2=>DUCT(DN2%DUCT_INDEX(1))
            IF (DU2%AREA < TWO_EPSILON_EB .OR. DU2%FIXED) CYCLE
            ARRAYLOC2 = NE%MATRIX_INDEX(DUCT_NE(DN2%DUCT_INDEX(1)))
            LHS(ARRAYLOC1,ARRAYLOC2) = LHS(ARRAYLOC1,ARRAYLOC2) - DN2%DIR(1)*DU2%AREA*DT_HV**2 / &
                                       (PSUM_TOT(DN%ZONE_INDEX)*DU%RHO_D*DU%LENGTH)
         ENDDO
      ENDIF
   ELSE
      IF (HVAC_LOCAL_PRESSURE) &
         LHS(ARRAYLOC1,ARRAYLOC1) = LHS(ARRAYLOC1,ARRAYLOC1) - DN%DIR(1)*DU%AREA/(DU%LENGTH*GAMMA*NODE_AREA_EX(DU%NODE_INDEX(1)))
      IF (DN%ZONE_INDEX >0) THEN
         PZ => P_ZONE(DN%ZONE_INDEX)
         DO NN=1,PZ%N_DUCTNODES
            DN2=>DUCTNODE(PZ%NODE_INDEX(NN))
            DU2=>DUCT(DN2%DUCT_INDEX(1))
            IF (DU2%AREA < TWO_EPSILON_EB .OR. DU2%FIXED) CYCLE
            ARRAYLOC2 = NE%MATRIX_INDEX(DUCT_NE(DN2%DUCT_INDEX(1)))
            LHS(ARRAYLOC1,ARRAYLOC2) = LHS(ARRAYLOC1,ARRAYLOC2) - DN2%DIR(1)*DU2%AREA*DT_HV**2 / &
                                       (PSUM_TOT(DN%ZONE_INDEX)*DU%RHO_D*DU%LENGTH)
         ENDDO
      ENDIF
   ENDIF
   DN=>DUCTNODE(DU%NODE_INDEX(2))
   IF (.NOT. (DN%VENT .OR. DN%GEOM)) THEN
      IF (.NOT. DN%AMBIENT .AND. .NOT. DN%LEAKAGE) THEN
         ARRAYLOC2 = NE%MATRIX_INDEX(NE%N_DUCTS+DUCTNODE_NE(DU%NODE_INDEX(2)))
         LHS(ARRAYLOC1,ARRAYLOC2) = DT_HV/(DU%RHO_D*DU%LENGTH)
      ENDIF
      IF (DN%LEAKAGE .AND. .NOT. DN%AMBIENT) THEN
         !IF (HVAC_LOCAL_PRESSURE) &
         !   LHS(ARRAYLOC1,ARRAYLOC1) = LHS(ARRAYLOC1,ARRAYLOC1)+DN%DIR(1)*DU%AREA/(DU%LENGTH*GAMMA*NODE_AREA_EX(DU%NODE_INDEX(2)))
         PZ => P_ZONE(DN%ZONE_INDEX)
         DO NN=1,PZ%N_DUCTNODES
            DN2=>DUCTNODE(PZ%NODE_INDEX(NN))
            DU2=>DUCT(DN2%DUCT_INDEX(1))
            IF (DU2%AREA < TWO_EPSILON_EB .OR. DU2%FIXED) CYCLE
            ARRAYLOC2 = NE%MATRIX_INDEX(DUCT_NE(DN2%DUCT_INDEX(1)))
            LHS(ARRAYLOC1,ARRAYLOC2) = LHS(ARRAYLOC1,ARRAYLOC2) + DN2%DIR(1)*DU2%AREA*DT_HV**2 / &
                                       (PSUM_TOT(DN%ZONE_INDEX)*DU%RHO_D*DU%LENGTH)
         ENDDO
      ENDIF
   ELSE
      IF (HVAC_LOCAL_PRESSURE) &
         LHS(ARRAYLOC1,ARRAYLOC1) = LHS(ARRAYLOC1,ARRAYLOC1) + DN%DIR(1)*DU%AREA/(DU%LENGTH*GAMMA*NODE_AREA_EX(DU%NODE_INDEX(2)))
      IF (DN%ZONE_INDEX >0) THEN
         PZ => P_ZONE(DN%ZONE_INDEX)
         DO NN=1,PZ%N_DUCTNODES
            DN2=>DUCTNODE(PZ%NODE_INDEX(NN))
            DU2=>DUCT(DN2%DUCT_INDEX(1))
            IF (DU2%AREA < TWO_EPSILON_EB .OR. DU2%FIXED) CYCLE
            ARRAYLOC2 = NE%MATRIX_INDEX(DUCT_NE(DN2%DUCT_INDEX(1)))
            LHS(ARRAYLOC1,ARRAYLOC2) = LHS(ARRAYLOC1,ARRAYLOC2) + DN2%DIR(1)*DU2%AREA*DT_HV**2 / &
                                       (PSUM_TOT(DN%ZONE_INDEX)*DU%RHO_D*DU%LENGTH)
         ENDDO
      ENDIF
   ENDIF
ENDDO DUCT_LOOP

END SUBROUTINE LHSDUCT

!> \brief Updates the pressure rise imposed by a fan
!>
!> \param T Current time (s)
!> \param ND Index indicating duct where fan is located

SUBROUTINE UPDATE_FAN(T,ND)
USE MATH_FUNCTIONS, ONLY : EVALUATE_RAMP
INTEGER :: FAN_ITER,DR_INDEX,DR_NF,DR_DU
INTEGER, INTENT(IN) :: ND
REAL(EB), INTENT(IN) :: T
REAL(EB) :: DEL_P,VDOT
REAL(EB) :: TSI,FLOW1,FLOW2,FLOW3,ERROR,DR_MAX,DR_MIN,CONST,P_SYS,FLOWGUESS,MAX_PRES,RAMP_FAC,REV_FAC,A,B,C
TYPE(DUCT_TYPE), POINTER::DU
TYPE(FAN_TYPE), POINTER::FA
TYPE(RAMPS_TYPE), POINTER :: RP

DU=> DUCT(ND)
FA=> FAN(DU%FAN_INDEX)
TSI = T - DU%FAN_ON_TIME

QFAN_IF: IF (DU%QFAN_INDEX < 0) THEN
   SELECT CASE (FA%FAN_TYPE)
      CASE(1) ! Constant flow
         DU%DP_FAN(NEW) = 0._EB
         RETURN
      CASE(2) ! Quadratic
         VDOT = (0.25_EB*DU%VEL(NEW)+0.25_EB*DU%VEL(PREVIOUS)+0.5_EB*DU%VEL(OLD))*DU%AREA
         IF (DU%REVERSE) VDOT = -VDOT
         VDOT = MAX(0._EB,MIN(VDOT,FA%MAX_FLOW))
         DEL_P = FA%MAX_PRES - FA%MAX_PRES*(VDOT/FA%MAX_FLOW)**2
         DEL_P = DEL_P*EVALUATE_RAMP(TSI,FA%SPIN_INDEX,TAU=FA%TAU)
      CASE(3) !Fan curve
         VDOT = 0.5_EB*(DU%VEL(NEW)+DU%VEL(OLD))*DU%AREA
         IF (DU%REVERSE) VDOT = -VDOT
         DEL_P = EVALUATE_RAMP(VDOT,FA%RAMP_INDEX)*EVALUATE_RAMP(TSI,FA%SPIN_INDEX,TAU=FA%TAU)
   END SELECT
ELSE QFAN_IF
   DR_INDEX = DU%DUCTRUN
   DR_NF = DU%QFAN_INDEX
   DR_DU = DU%DUCTRUN_M_INDEX
   IF (DU%REVERSE) THEN
      REV_FAC = -1._EB
   ELSE
      REV_FAC = 1._EB
   ENDIF
   DR_MAX = DUCTRUN(DR_INDEX)%VEL(DR_DU,DR_NF,1)*DU%AREA*REV_FAC
   DR_MIN = DUCTRUN(DR_INDEX)%VEL(DR_DU,0,1)*DU%AREA*REV_FAC
   RAMP_FAC = EVALUATE_RAMP(TSI,FA%SPIN_INDEX,TAU=FA%TAU)
   IF (RAMP_FAC < TWO_EPSILON_EB) THEN
      DU%DP_FAN(NEW) = 0._EB
      RETURN
   ENDIF
   ! DEAL WITH REVERSE FLOW IN FAN
   SELECT CASE (FA%FAN_TYPE)
      CASE(2) ! Quadratic
         MAX_PRES = FA%MAX_PRES*RAMP_FAC
         FLOW1 = FA%MAX_FLOW*RAMP_FAC
         IF (DUCTRUN(DR_INDEX)%VEL(DR_DU,0,1)*REV_FAC*DU%AREA > FLOW1) THEN
            DEL_P=0._EB
         ELSEIF (DUCTRUN(DR_INDEX)%VEL(DR_DU,DR_NF,1)*REV_FAC <= 0._EB) THEN
            DEL_P = MAX_PRES
         ELSE !both system and fan curve are quadratic, can get analytic solution
            A = 1._EB/FLOW1**2 + 1._EB/(DR_MAX - DR_MIN)**2
            B = -2._EB*DR_MIN/(DR_MAX - DR_MIN)**2
            C = DR_MIN**2/(DR_MAX - DR_MIN)**2 - 1._EB
            FLOWGUESS = (-B + SQRT(B**2-4._EB*A*C))/(2._EB*A)
            DEL_P = (1._EB-(FLOWGUESS/FLOW1)**2)*MAX_PRES
         ENDIF
      CASE(3) !Fan curve
         RP => RAMPS(FA%RAMP_INDEX)
         MAX_PRES = MAXVAL(RP%DEPENDENT_DATA)*RAMP_FAC
         CONST = MAX_PRES/((DUCTRUN(DR_INDEX)%VEL(DR_DU,DR_NF,1)-DUCTRUN(DR_INDEX)%VEL(DR_DU,0,1))*DU%AREA)**2
         ERROR = 1.E-4_EB
         FLOW1 = MAX(RP%T_MIN*RAMP_FAC,DR_MIN)
         FLOW2 = MIN(DR_MAX,RP%T_MAX*RAMP_FAC)
         FLOWGUESS = MAX(FLOW1,MIN(FLOW2,DU%VEL(OLD)*REV_FAC*DU%AREA))
         IF (DUCTRUN(DR_INDEX)%VEL(DR_DU,0,1)*REV_FAC > RP%T_MAX*RAMP_FAC) THEN
            DEL_P=RP%DEPENDENT_DATA(RP%NUMBER_DATA_POINTS)*RAMP_FAC
         ELSEIF (DUCTRUN(DR_INDEX)%VEL(DR_DU,DR_NF,1)*REV_FAC <= RP%T_MIN) THEN
            DEL_P=RP%DEPENDENT_DATA(1)*RAMP_FAC
         ELSE
            FAN_ITER = 0
            FAN_LOOP2: DO
               FAN_ITER = FAN_ITER + 1
               P_SYS = CONST * (FLOWGUESS - DUCTRUN(DR_INDEX)%VEL(DR_DU,0,1)*REV_FAC*DU%AREA)**2
               DEL_P = EVALUATE_RAMP(FLOWGUESS/RAMP_FAC,FA%RAMP_INDEX)*RAMP_FAC
               IF (P_SYS > DEL_P) THEN
                  FLOW2 = FLOWGUESS
                  FLOW3 = 0.5_EB * (FLOWGUESS + FLOW1)
               ELSEIF (P_SYS < DEL_P) THEN
                  FLOW1 = FLOWGUESS
                  FLOW3 = 0.5_EB * (FLOWGUESS + FLOW2)
               ELSE
                  EXIT FAN_LOOP2
               ENDIF
               IF (ABS(FLOW3 - FLOWGUESS)/(FLOWGUESS+TWO_EPSILON_EB) < ERROR) EXIT FAN_LOOP2
               IF (FAN_ITER==20) THEN
                  FLOWGUESS = 0.5_EB*(FLOW3 + FLOWGUESS)
                  EXIT FAN_LOOP2
               ENDIF
              FLOWGUESS = FLOW3
            ENDDO FAN_LOOP2
            ! Output fan pressure equating to output flow rate
            DEL_P = EVALUATE_RAMP(FLOWGUESS/RAMP_FAC,FA%RAMP_INDEX)*RAMP_FAC
         ENDIF
   END SELECT
ENDIF QFAN_IF

IF (DU%REVERSE) DEL_P=-DEL_P
DU%DP_FAN(NEW) = DEL_P

END SUBROUTINE UPDATE_FAN

!> \brief Averages gas properties at VENTs connected to HVAC system
!>
!> \param NM Current mesh

SUBROUTINE HVAC_BC_IN(NM)

USE PHYSICAL_FUNCTIONS, ONLY: GET_ENTHALPY
USE COMPLEX_GEOMETRY, ONLY: CC_FGSC, CC_SOLID
INTEGER, INTENT(IN) :: NM
INTEGER :: II,JJ,KK,IW,IOR,IZ1,IZ2,ICF
INTEGER, POINTER :: NODE_INDEX,VENT_INDEX
LOGICAL :: WALL_FLAG,ERROR_FLAG
REAL(EB) :: ZZ_GET(1:N_TRACKED_SPECIES),AREA,P_AVE,HGAS
REAL(EB), POINTER, DIMENSION(:,:,:) :: RHOP,UP,VP,WP,HP
REAL(EB), POINTER, DIMENSION(:,:) :: PBARP
REAL(EB) :: TNOW !< Current CPU time (s) used in computing length of time spent in HVAC routines.
TYPE (SURFACE_TYPE), POINTER :: SF
TYPE (WALL_TYPE), POINTER :: WC
TYPE (CFACE_TYPE), POINTER :: CFA
TYPE (BOUNDARY_PROP1_TYPE), POINTER :: B1
TYPE (BOUNDARY_COORD_TYPE), POINTER :: BC

TNOW=CURRENT_TIME()

CALL POINT_TO_MESH(NM)
ERROR_FLAG = .FALSE.
ZZ_GET = 0._EB

IF (PREDICTOR) THEN
   PBARP => PBAR
   RHOP  => RHO
   HP    => H
   UP    => U
   VP    => V
   WP    => W
ELSE
   PBARP => PBAR_S
   RHOP  => RHOS
   HP    => HS
   UP    => US
   VP    => VS
   WP    => WS
ENDIF

WALL_FLAG = .TRUE.
WALL_LOOP: DO IW = 1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
   WC => WALL(IW)
   IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE WALL_LOOP
   B1 => BOUNDARY_PROP1(WC%B1_INDEX)
   BC => BOUNDARY_COORD(WC%BC_INDEX)
   SF => SURFACE(WC%SURF_INDEX)
   NODE_INDEX => B1%NODE_INDEX
   VENT_INDEX => WC%VENT_INDEX
   ! Reset indices to zero in case of adding an OSBT
   NODE_INDEX = 0
   CALL INITIALIZE_HVAC
   IF (ERROR_FLAG) RETURN
ENDDO WALL_LOOP

WALL_FLAG = .FALSE.

CFACE_LOOP: DO ICF=INTERNAL_CFACE_CELLS_LB+1,INTERNAL_CFACE_CELLS_LB+N_INTERNAL_CFACE_CELLS
   CFA => CFACE(ICF)
   B1 => BOUNDARY_PROP1(CFA%B1_INDEX)
   BC => BOUNDARY_COORD(CFA%BC_INDEX)
   SF => SURFACE(CFA%SURF_INDEX)
   NODE_INDEX => B1%NODE_INDEX
   ! Reset indices to zero in case of adding an OSBT
   NODE_INDEX = 0
   CALL INITIALIZE_HVAC
ENDDO CFACE_LOOP

T_USED(13)=T_USED(13)+CURRENT_TIME()-TNOW

CONTAINS

SUBROUTINE INITIALIZE_HVAC

IF (CELL(CELL_INDEX(BC%IIG,BC%JJG,BC%KKG))%SOLID) RETURN
IF (CC_IBM .AND. WALL_FLAG .AND. VENT_INDEX > 0) THEN
   IF (FCVAR(BC%IIG,BC%JJG,BC%KKG,CC_FGSC,ABS(VENTS(VENT_INDEX)%IOR)) == CC_SOLID) RETURN
ENDIF

ZONE_LEAK_IF: IF (ALL(SF%LEAK_PATH < 0)) THEN


   IF (WALL_FLAG) THEN
      IF (VENT_INDEX == 0) RETURN
      IF (VENTS(VENT_INDEX)%NODE_INDEX < 0) RETURN
      NODE_INDEX = VENTS(VENT_INDEX)%NODE_INDEX
      IF (WC%OBST_INDEX > 0) THEN
         IF (OBSTRUCTION(WC%OBST_INDEX)%REMOVABLE) THEN
            ERROR_FLAG = .TRUE.
            WRITE(MESSAGE,'(3A)') 'ERROR(564): VENT: ',TRIM(VENTS(VENT_INDEX)%ID),&
               ', is used for HVAC and attached to a removable OBST.'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.); RETURN
         ENDIF
      ENDIF
   ELSE
      NODE_INDEX = CFA%NODE_INDEX
      IF (NODE_INDEX<=0) RETURN
   ENDIF
   
   IOR = BC%IOR
   II = BC%IIG
   JJ = BC%JJG
   KK = BC%KKG
   AREA = B1%AREA
   IF (B1%PRESSURE_ZONE /= NODE_ZONE(NODE_INDEX)) THEN
      IF (NODE_ZONE(NODE_INDEX) == 0) THEN
         NODE_ZONE(NODE_INDEX) = B1%PRESSURE_ZONE
      ELSE
         ERROR_FLAG = .TRUE.
         WRITE(MESSAGE,'(A,A)') 'ERROR(552): Ductnode must lie with a single pressure zone. Node: ',TRIM(DUCTNODE(NODE_INDEX)%ID)
         CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.); RETURN
      ENDIF
   ENDIF

   NODE_AREA(NODE_INDEX) = NODE_AREA(NODE_INDEX) + AREA
   NODE_RHO(NODE_INDEX) = NODE_RHO(NODE_INDEX) + AREA/RHOP(II,JJ,KK)
   ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(II,JJ,KK,1:N_TRACKED_SPECIES)
   NODE_ZZ(NODE_INDEX,1:N_TRACKED_SPECIES) = NODE_ZZ(NODE_INDEX,1:N_TRACKED_SPECIES) + ZZ_GET(1:N_TRACKED_SPECIES)*AREA
   CALL GET_ENTHALPY(ZZ_GET,HGAS,TMP(II,JJ,KK))
   NODE_TMP(NODE_INDEX) = NODE_TMP(NODE_INDEX) + TMP(II,JJ,KK)*AREA
   NODE_H(NODE_INDEX) = NODE_H(NODE_INDEX) + HGAS * AREA

   SELECT CASE (IOR)
      CASE (1)
         NODE_X(NODE_INDEX) = NODE_X(NODE_INDEX) + X(II-1)*AREA
         NODE_Y(NODE_INDEX) = NODE_Y(NODE_INDEX) + YC(JJ)*AREA
         NODE_Z(NODE_INDEX) = NODE_Z(NODE_INDEX) + ZC(KK)*AREA
      CASE(-1)
         NODE_X(NODE_INDEX) = NODE_X(NODE_INDEX) + X(II)*AREA
         NODE_Y(NODE_INDEX) = NODE_Y(NODE_INDEX) + YC(JJ)*AREA
         NODE_Z(NODE_INDEX) = NODE_Z(NODE_INDEX) + ZC(KK)*AREA
      CASE (2)
         NODE_X(NODE_INDEX) = NODE_X(NODE_INDEX) + XC(II)*AREA
         NODE_Y(NODE_INDEX) = NODE_Y(NODE_INDEX) + Y(JJ-1)*AREA
         NODE_Z(NODE_INDEX) = NODE_Z(NODE_INDEX) + ZC(KK)*AREA
      CASE(-2)
         NODE_X(NODE_INDEX) = NODE_X(NODE_INDEX) + XC(II)*AREA
         NODE_Y(NODE_INDEX) = NODE_Y(NODE_INDEX) + Y(JJ)*AREA
         NODE_Z(NODE_INDEX) = NODE_Z(NODE_INDEX) + ZC(KK)*AREA
      CASE (3)
         NODE_X(NODE_INDEX) = NODE_X(NODE_INDEX) + XC(II)*AREA
         NODE_Y(NODE_INDEX) = NODE_Y(NODE_INDEX) + YC(JJ)*AREA
         NODE_Z(NODE_INDEX) = NODE_Z(NODE_INDEX) + Z(KK-1)*AREA
      CASE (-3)
         NODE_X(NODE_INDEX) = NODE_X(NODE_INDEX) + XC(II)*AREA
         NODE_Y(NODE_INDEX) = NODE_Y(NODE_INDEX) + YC(JJ)*AREA
         NODE_Z(NODE_INDEX) = NODE_Z(NODE_INDEX) + Z(KK)*AREA
   END SELECT
   IF (HVAC_LOCAL_PRESSURE) THEN
      SELECT CASE (IOR)
         CASE DEFAULT
            P_AVE = PBARP(KK,B1%PRESSURE_ZONE)
         CASE (3)
            P_AVE = 0.5_EB*(PBARP(KK-1,B1%PRESSURE_ZONE)+PBARP(KK,B1%PRESSURE_ZONE))
         CASE (-3)
            P_AVE = 0.5_EB*(PBARP(KK,B1%PRESSURE_ZONE)+PBARP(KK+1,B1%PRESSURE_ZONE))
      END SELECT
      NODE_P(NODE_INDEX) = NODE_P(NODE_INDEX) + (P_AVE-P_INF+RHO(II,JJ,KK)*(HP(II,JJ,KK)-KRES(II,JJ,KK)))*AREA
   ELSE
      SELECT CASE (IOR)
         CASE (1)
            NODE_P(NODE_INDEX) = NODE_P(NODE_INDEX) + &
                                       (PBARP(KK,B1%PRESSURE_ZONE)-P_INF-RHO(II,JJ,KK) * &
                                       0.5_EB*(UP(II-1,JJ,KK)+B1%U_NORMAL)**2 * &
                                       SIGN(1._EB,UP(II-1,JJ,KK)+B1%U_NORMAL))* AREA
         CASE(-1)
            NODE_P(NODE_INDEX) = NODE_P(NODE_INDEX) + &
                                       (PBARP(KK,B1%PRESSURE_ZONE)-P_INF+RHO(II,JJ,KK) * &
                                       0.5_EB*(UP(II,JJ,KK)-B1%U_NORMAL)**2  * &
                                       SIGN(1._EB,UP(II,JJ,KK)-B1%U_NORMAL))*AREA
         CASE (2)
            NODE_P(NODE_INDEX) = NODE_P(NODE_INDEX) + &
                                       (PBARP(KK,B1%PRESSURE_ZONE)-P_INF-RHO(II,JJ,KK) * &
                                       0.5_EB*(VP(II,JJ-1,KK)+B1%U_NORMAL)**2 * &
                                       SIGN(1._EB,VP(II,JJ-1,KK)+B1%U_NORMAL))*AREA
         CASE(-2)
            NODE_P(NODE_INDEX) = NODE_P(NODE_INDEX) + &
                                       (PBARP(KK,B1%PRESSURE_ZONE)-P_INF+RHO(II,JJ,KK) * &
                                       0.5_EB*(VP(II,JJ,KK)-B1%U_NORMAL)**2 * &
                                       SIGN(1._EB,VP(II,JJ,KK)-B1%U_NORMAL))*AREA
         CASE (3)
            NODE_P(NODE_INDEX) = NODE_P(NODE_INDEX) + &
                                       (PBARP(KK,B1%PRESSURE_ZONE)-P_INF-RHO(II,JJ,KK) * &
                                       0.5_EB*(WP(II,JJ,KK-1)+B1%U_NORMAL)**2 * &
                                       SIGN(1._EB,WP(II,JJ,KK-1)+B1%U_NORMAL))*AREA
         CASE (-3)
            NODE_P(NODE_INDEX) = NODE_P(NODE_INDEX) + &
                                       (PBARP(KK,B1%PRESSURE_ZONE)-P_INF+RHO(II,JJ,KK) * &
                                       0.5_EB*(WP(II,JJ,KK)-B1%U_NORMAL)**2 * &
                                       SIGN(1._EB,WP(II,JJ,KK)-B1%U_NORMAL))*AREA
      END SELECT
   ENDIF

ELSE ZONE_LEAK_IF

   IF (ALL(SF%LEAK_PATH/=B1%PRESSURE_ZONE)) RETURN
   IF (SF%LEAK_PATH(1) == B1%PRESSURE_ZONE) THEN
      IZ1 = SF%LEAK_PATH(1)
      IZ2 = SF%LEAK_PATH(2)
      NODE_INDEX = DUCT(LEAK_PATH(MIN(IZ1,IZ2),MAX(IZ1,IZ2)))%NODE_INDEX(1)
   ELSE
      IZ1 = SF%LEAK_PATH(2)
      IZ2 = SF%LEAK_PATH(1)
      NODE_INDEX = DUCT(LEAK_PATH(MIN(IZ1,IZ2),MAX(IZ1,IZ2)))%NODE_INDEX(2)
   ENDIF
   IF (NODE_ZONE(NODE_INDEX) == 0) NODE_ZONE(NODE_INDEX) = B1%PRESSURE_ZONE
   IOR = BC%IOR
   II = BC%IIG
   JJ = BC%JJG
   KK = BC%KKG
   AREA = B1%AREA

   NODE_AREA(NODE_INDEX) = NODE_AREA(NODE_INDEX) + AREA
   NODE_RHO(NODE_INDEX) = NODE_RHO(NODE_INDEX) + AREA/RHOP(II,JJ,KK)

   ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(II,JJ,KK,1:N_TRACKED_SPECIES)
   NODE_ZZ(NODE_INDEX,1:N_TRACKED_SPECIES) = NODE_ZZ(NODE_INDEX,1:N_TRACKED_SPECIES) + ZZ_GET(1:N_TRACKED_SPECIES)*AREA
   CALL GET_ENTHALPY(ZZ_GET,HGAS,TMP(II,JJ,KK))
   NODE_TMP(NODE_INDEX) = NODE_TMP(NODE_INDEX) + TMP(II,JJ,KK)*AREA
   NODE_H(NODE_INDEX) = NODE_H(NODE_INDEX) + HGAS * AREA

   SELECT CASE (IOR)
      CASE (3)
         P_AVE = 0.5_EB*(PBARP(KK-1,B1%PRESSURE_ZONE)+PBARP(KK,B1%PRESSURE_ZONE))-P_INF
         NODE_Z(NODE_INDEX) = NODE_Z(NODE_INDEX) + Z(KK-1)*AREA
      CASE (-3)
         P_AVE = 0.5_EB*(PBARP(KK,B1%PRESSURE_ZONE)+PBARP(KK+1,B1%PRESSURE_ZONE))-P_INF
         NODE_Z(NODE_INDEX) = NODE_Z(NODE_INDEX) + Z(KK)*AREA
      CASE DEFAULT
         P_AVE = PBARP(KK,B1%PRESSURE_ZONE)-P_INF
         NODE_Z(NODE_INDEX) = NODE_Z(NODE_INDEX) + ZC(KK)*AREA
   END SELECT
   NODE_P(NODE_INDEX) = NODE_P(NODE_INDEX) + P_AVE*AREA

ENDIF ZONE_LEAK_IF

END SUBROUTINE INITIALIZE_HVAC

END SUBROUTINE HVAC_BC_IN

!> \brief Determines what ducts are "fixed"; i.e. they have been given fixed volume or mass flow (including fixed volume fan)
!>
!> \param T Current time (s)

SUBROUTINE DETERMINE_FIXED_ELEMENTS(T,CHANGE)
USE MATH_FUNCTIONS,ONLY:EVALUATE_RAMP
INTEGER:: NN,ND, COUNTER
REAL(EB), INTENT(IN):: T
LOGICAL, INTENT(INOUT) :: CHANGE
REAL(EB) :: VOLUME_FLOW,TSI
LOGICAL :: ZERO_AREA,FAN_OP_OLD
TYPE(FAN_TYPE), POINTER :: FA
TYPE(DUCT_TYPE), POINTER :: DU
TYPE(DUCTNODE_TYPE), POINTER :: DN

DUCT%FIXED = .FALSE.
DUCT%VOLUME_FLOW = 0._EB

DUCT_LOOP: DO ND=1,N_DUCTS
   DU=>DUCT(ND)

   ! Set duct area
   IF (DU%DAMPER) THEN
      IF (DU%DEVC_INDEX > 0) THEN
         IF (DEVICE(DU%DEVC_INDEX)%CURRENT_STATE .NEQV. DU%DAMPER_OPEN) DU%DAMPER_OPEN = DEVICE(DU%DEVC_INDEX)%CURRENT_STATE
      ELSEIF(DU%CTRL_INDEX > 0) THEN
         IF (CONTROL(DU%CTRL_INDEX)%CURRENT_STATE .NEQV. DU%DAMPER_OPEN) DU%DAMPER_OPEN = CONTROL(DU%CTRL_INDEX)%CURRENT_STATE
      ENDIF
      IF (DU%DAMPER_OPEN) THEN
         DU%AREA = DU%AREA_INITIAL
      ELSE
         DU%AREA = 0._EB
         DU%VEL  = 0._EB
      ENDIF
   ENDIF

   ! If leakage and no vent area then set duct area to zero
   IF (DU%LEAKAGE) THEN
      IF ((NODE_AREA_EX(DU%NODE_INDEX(1)) < TWO_EPSILON_EB .AND. .NOT. DUCTNODE(DU%NODE_INDEX(1))%AMBIENT) &
      .OR.(NODE_AREA_EX(DU%NODE_INDEX(2)) < TWO_EPSILON_EB .AND. .NOT. DUCTNODE(DU%NODE_INDEX(2))%AMBIENT)) DU%AREA=0._EB
   ENDIF

   IF (DU%AREA_OLD <TWO_EPSILON_EB .NEQV. DU%AREA<TWO_EPSILON_EB) CHANGE = .TRUE.
   DU%AREA_OLD = DU%AREA

   ! Set fixed flows
   IF (DU%VOLUME_FLOW_INITIAL<1.E6_EB) THEN
      DU%VOLUME_FLOW = DU%VOLUME_FLOW_INITIAL*EVALUATE_RAMP(T,DU%RAMP_INDEX,TAU=DU%TAU)
      IF(DU%AREA > TWO_EPSILON_EB) DU%VEL(NEW) = DU%VOLUME_FLOW/DU%AREA
      DU%FIXED = .TRUE.
      CYCLE DUCT_LOOP
   ENDIF
   IF (DU%MASS_FLOW_INITIAL<1.E6_EB) THEN
      IF (DU%MASS_FLOW_INITIAL > 0._EB) THEN
         DU%VOLUME_FLOW = DU%MASS_FLOW_INITIAL/DUCTNODE(DU%NODE_INDEX(1))%RHO*EVALUATE_RAMP(T,DU%RAMP_INDEX,TAU=DU%TAU)
      ELSE
         DU%VOLUME_FLOW = DU%MASS_FLOW_INITIAL/DUCTNODE(DU%NODE_INDEX(2))%RHO*EVALUATE_RAMP(T,DU%RAMP_INDEX,TAU=DU%TAU)
      ENDIF
      IF(DU%AREA > TWO_EPSILON_EB) DU%VEL(NEW) = DU%VOLUME_FLOW/DU%AREA
      DU%FIXED = .TRUE.
      CYCLE DUCT_LOOP
   ENDIF

   ! Set fan on and fixed flow fan
   IF (DU%FAN_INDEX>0) THEN
      FAN_OP_OLD = DU%FAN_OPERATING
      IF (DU%DEVC_INDEX > 0) THEN
         DU%FAN_OPERATING = DEVICE(DU%DEVC_INDEX)%CURRENT_STATE
         IF (DU%FAN_OPERATING .AND. (DEVICE(DU%DEVC_INDEX)%CURRENT_STATE .NEQV. DEVICE(DU%DEVC_INDEX)%PRIOR_STATE)) &
            DU%FAN_ON_TIME = T
      ELSEIF (DU%CTRL_INDEX > 0) THEN
         DU%FAN_OPERATING = CONTROL(DU%CTRL_INDEX)%CURRENT_STATE
         IF (DU%FAN_OPERATING .AND. (CONTROL(DU%CTRL_INDEX)%CURRENT_STATE .NEQV. CONTROL(DU%CTRL_INDEX)%PRIOR_STATE)) &
            DU%FAN_ON_TIME = T
      ENDIF
      IF (FAN(DU%FAN_INDEX)%FAN_TYPE==1) THEN
         IF (FAN_OP_OLD .NEQV. DU%FAN_OPERATING) CHANGE = .TRUE.
         IF (DU%FAN_OPERATING) THEN
            DU%FIXED=.TRUE.
            FA => FAN(DU%FAN_INDEX)
            TSI = T - DU%FAN_ON_TIME
            DU%VEL(NEW) = FA%VOL_FLOW/DU%AREA*EVALUATE_RAMP(TSI,FA%SPIN_INDEX,TAU=FA%TAU)
            IF (DU%REVERSE) DU%VEL(NEW)=-DU%VEL(NEW)
            DU%VOLUME_FLOW = DU%VEL(NEW)*DU%AREA
         ELSE
            DU%FIXED = .FALSE.
            DU%VOLUME_FLOW = 0._EB
         ENDIF
      ENDIF
   ENDIF
ENDDO DUCT_LOOP

NODE_LOOP: DO NN = 1, N_DUCTNODES
      DN=>DUCTNODE(NN)
      IF (DN%VENT .OR. DN%AMBIENT .OR. DN%LEAKAGE .OR. DN%GEOM) CYCLE NODE_LOOP
      COUNTER = 0
      VOLUME_FLOW = 0._EB
      ZERO_AREA = .FALSE.
      DO ND = 1,DN%N_DUCTS
         DU=>DUCT(DN%DUCT_INDEX(ND))
         IF (DU%FIXED) THEN
            COUNTER = COUNTER + 1
            VOLUME_FLOW = VOLUME_FLOW + DU%VOLUME_FLOW
         ENDIF
         IF (DU%AREA < TWO_EPSILON_EB) THEN
            COUNTER = COUNTER + 1
            ZERO_AREA = .TRUE.
         ENDIF
      ENDDO
      IF (COUNTER==DN%N_DUCTS) THEN
         IF (ZERO_AREA .AND. ABS(VOLUME_FLOW)<TWO_EPSILON_EB) CYCLE NODE_LOOP
         WRITE(MESSAGE,'(A,A)') 'ERROR(553): Cannot specify fixed flows for all branches of internal ductnode:',TRIM(DN%ID)
         CALL SHUTDOWN(MESSAGE); RETURN
      ENDIF
ENDDO NODE_LOOP

END SUBROUTINE DETERMINE_FIXED_ELEMENTS


!> \brief Determines what HVAC components lies in different, isolated networks (i.e. not sharing a common pressure zone)
!> \param CHANGE Flag to force revaluation of the duct networks. Otherwise done only if there area damper changes.
!> \param T Current time (s)

SUBROUTINE FIND_NETWORKS(CHANGE,T)
INTEGER:: NZ,NN,ND,DUCT_COUNTER(N_DUCTS),NODE_COUNTER(N_DUCTNODES),COUNTER,COUNTER2,ZONE_COUNTER(N_ZONE)
INTEGER, DIMENSION(:), ALLOCATABLE :: NETWORK_DCOUNTER,NETWORK_NCOUNTER,RENUMBER
LOGICAL, INTENT(OUT) :: CHANGE
LOGICAL :: CHANGE2
REAL(EB), INTENT(IN):: T
TYPE(DUCT_TYPE), POINTER :: DU

CHANGE = .FALSE.
IF (.NOT. ALLOCATED(NETWORK)) CHANGE = .TRUE.

IF (N_ZONE > 0) ZONE_COUNTER = 0

CALL DETERMINE_FIXED_ELEMENTS(T,CHANGE)

DO NN = 1, N_DUCTNODES
   NZ = DUCTNODE(NN)%ZONE_INDEX
   IF (NZ>=1) THEN
      NODE_COUNTER(NN) = NZ
      ZONE_COUNTER(NZ) = ZONE_COUNTER(NZ) + 1
   ELSE
      NODE_COUNTER(NN) = NN+N_ZONE
   ENDIF
ENDDO
IF (N_ZONE > 0) THEN
   DO NZ = 1, N_ZONE
      IF (ALLOCATED(P_ZONE(NZ)%NODE_INDEX)) DEALLOCATE(P_ZONE(NZ)%NODE_INDEX)
      ALLOCATE(P_ZONE(NZ)%NODE_INDEX(ZONE_COUNTER(NZ)))
      P_ZONE(NZ)%N_DUCTNODES = ZONE_COUNTER(NZ)
      COUNTER = 1
      DO NN = 1,N_DUCTNODES
         IF (DUCTNODE(NN)%ZONE_INDEX == NZ) THEN
            P_ZONE(NZ)%NODE_INDEX(COUNTER)=NN
            COUNTER = COUNTER + 1
         ENDIF
      ENDDO
   ENDDO
ENDIF


IF (CHANGE) THEN
   IF (ALLOCATED(NETWORK)) DEALLOCATE(NETWORK)

   CHANGE2 = .TRUE.
   DO WHILE (CHANGE2)
      CHANGE2 = .FALSE.
      DO ND = 1, N_DUCTS
         DU => DUCT(ND)
         IF (NODE_COUNTER(DU%NODE_INDEX(1)) /= NODE_COUNTER(DU%NODE_INDEX(2))) THEN
            CHANGE2 = .TRUE.
            COUNTER = MIN(NODE_COUNTER(DU%NODE_INDEX(1)),NODE_COUNTER(DU%NODE_INDEX(2)))
            DUCT_COUNTER(ND) = COUNTER
            NODE_COUNTER(DU%NODE_INDEX(1)) = COUNTER
            NODE_COUNTER(DU%NODE_INDEX(2)) = COUNTER
         ELSE
            DUCT_COUNTER(ND) = NODE_COUNTER(DU%NODE_INDEX(1))
         ENDIF
      ENDDO

      IF (N_ZONE > 0) THEN
         DO NZ = 1, N_ZONE
            COUNTER = 1
            COUNTER2 = 1
            DO NN = 1, P_ZONE(NZ)%N_DUCTNODES
               IF (NN==1) THEN
                  COUNTER = NODE_COUNTER(P_ZONE(NZ)%NODE_INDEX(NN))
                  COUNTER2 = COUNTER
               ELSE
                  IF (COUNTER /= NODE_COUNTER(P_ZONE(NZ)%NODE_INDEX(NN))) &
                     COUNTER2 = MAX(COUNTER2,NODE_COUNTER(P_ZONE(NZ)%NODE_INDEX(NN)))
                     COUNTER  = MIN(COUNTER,NODE_COUNTER(P_ZONE(NZ)%NODE_INDEX(NN)))
               ENDIF
            ENDDO
            IF (COUNTER /= COUNTER2) THEN
               CHANGE2 = .TRUE.
               DO NN = 1, P_ZONE(NZ)%N_DUCTNODES
                  NODE_COUNTER(P_ZONE(NZ)%NODE_INDEX(NN)) = COUNTER
               ENDDO
            ENDIF
         ENDDO
      ENDIF
   END DO

   ALLOCATE(RENUMBER(MAXVAL(NODE_COUNTER)))
   RENUMBER = 0
   N_NETWORKS = 0
   DO NN = 1,N_DUCTNODES
      IF (RENUMBER(NODE_COUNTER(NN)) == 0) THEN
         N_NETWORKS = N_NETWORKS + 1
         RENUMBER(NODE_COUNTER(NN)) = N_NETWORKS
      ENDIF
   ENDDO
   DO NN = 1,N_DUCTNODES
      NODE_COUNTER(NN) = RENUMBER(NODE_COUNTER(NN))
   ENDDO
   DO ND = 1,N_DUCTS
      DUCT_COUNTER(ND) = RENUMBER(DUCT_COUNTER(ND))
   ENDDO
   DEALLOCATE(RENUMBER)

   ALLOCATE(NETWORK(N_NETWORKS))
   NETWORK%N_DUCTS=0
   NETWORK%N_DUCTNODES=0
   ALLOCATE(NETWORK_DCOUNTER(N_NETWORKS))
   NETWORK_DCOUNTER=0
   ALLOCATE(NETWORK_NCOUNTER(N_NETWORKS))
   NETWORK_NCOUNTER=0
   COUNTER = 0
   DO ND = 1, N_DUCTS
      NETWORK(DUCT_COUNTER(ND))%N_DUCTS = NETWORK(DUCT_COUNTER(ND))%N_DUCTS + 1
   ENDDO
   DO NN = 1, N_DUCTNODES
      NETWORK(NODE_COUNTER(NN))%N_DUCTNODES = NETWORK(NODE_COUNTER(NN))%N_DUCTNODES + 1
   ENDDO
   DO NN = 1, N_NETWORKS
      ALLOCATE(NETWORK(NN)%DUCT_INDEX(NETWORK(NN)%N_DUCTS))
      ALLOCATE(NETWORK(NN)%NODE_INDEX(NETWORK(NN)%N_DUCTNODES))
      ALLOCATE(NETWORK(NN)%MATRIX_INDEX(NETWORK(NN)%N_DUCTS+NETWORK(NN)%N_DUCTNODES))
      NETWORK(NN)%MATRIX_INDEX = 0
   ENDDO
   DO ND = 1, N_DUCTS
      NETWORK_DCOUNTER(DUCT_COUNTER(ND)) = NETWORK_DCOUNTER(DUCT_COUNTER(ND)) + 1
      NETWORK(DUCT_COUNTER(ND))%DUCT_INDEX(NETWORK_DCOUNTER(DUCT_COUNTER(ND))) = ND
      DUCT_NE(ND) = NETWORK_DCOUNTER(DUCT_COUNTER(ND))
   ENDDO
   DO NN = 1, N_DUCTNODES
      NETWORK_NCOUNTER(NODE_COUNTER(NN)) = NETWORK_NCOUNTER(NODE_COUNTER(NN)) + 1
      NETWORK(NODE_COUNTER(NN))%NODE_INDEX(NETWORK_NCOUNTER(NODE_COUNTER(NN))) = NN
      DUCTNODE_NE(NN) = NETWORK_NCOUNTER(NODE_COUNTER(NN))
   ENDDO
   DEALLOCATE(NETWORK_DCOUNTER)
   DEALLOCATE(NETWORK_NCOUNTER)
   CALL SETUP_SOLUTION_POINTERS
ENDIF

END SUBROUTINE FIND_NETWORKS

!> \brief For each HVAC network the routine defines which ducts and nodes belowng to which element of the HVAC matrix

SUBROUTINE SETUP_SOLUTION_POINTERS
INTEGER:: NNE,NN,ND,COUNTER
TYPE(DUCT_TYPE), POINTER :: DU
TYPE(DUCTNODE_TYPE), POINTER :: DN
TYPE(NETWORK_TYPE), POINTER :: NE

DO NNE = 1,N_NETWORKS
   COUNTER = 0
   NE => NETWORK(NNE)
   DO ND=1,NE%N_DUCTS
      DU=>DUCT(NE%DUCT_INDEX(ND))
      IF (DU%FIXED .OR. DU%AREA<=TWO_EPSILON_EB) CYCLE
      COUNTER = COUNTER + 1
      NE%MATRIX_INDEX(ND)=COUNTER
   ENDDO
   DO NN=1,NE%N_DUCTNODES
      DN=>DUCTNODE(NE%NODE_INDEX(NN))
      IF (DN%VENT .OR. DN%AMBIENT .OR. DN%LEAKAGE .OR. DN%GEOM) CYCLE
      COUNTER = COUNTER + 1
      NE%MATRIX_INDEX(NE%N_DUCTS+NN)=COUNTER
   ENDDO
   NE%N_MATRIX=COUNTER
ENDDO

END SUBROUTINE SETUP_SOLUTION_POINTERS

!> \brief Determines wall friction loss and assigns node losses to ducts
!>
!> \param T Current time (s)
!> \param DT Current time step (s)
!> \param NNE Index indicating which HVAC network is being solved

SUBROUTINE UPDATE_LOSS(T,DT,NNE)
USE PHYSICAL_FUNCTIONS,ONLY:GET_VISCOSITY
USE MATH_FUNCTIONS,ONLY:EVALUATE_RAMP
REAL(EB) :: FRICTION_FACTOR,LOSS_SUM,ZZ_GET(1:N_TRACKED_SPECIES),VISCOSITY,AREA,VFLOW
INTEGER, INTENT(IN) :: NNE
REAL(EB), INTENT(IN) :: T,DT
INTEGER :: ND,ND2, NN,NUM_OUT
TYPE(DUCT_TYPE), POINTER :: DU,DU2
TYPE(DUCTNODE_TYPE), POINTER :: DN
TYPE(NETWORK_TYPE), POINTER :: NE

NE => NETWORK(NNE)

DO ND = 1, NE%N_DUCTS
   DUCT(NE%DUCT_INDEX(ND))%TOTAL_LOSS=0._EB
ENDDO
NODELOOP : DO NN=1,NE%N_DUCTNODES
   VFLOW = 0._EB
   DN => DUCTNODE(NE%NODE_INDEX(NN))
   IF (DN%LEAKAGE) CYCLE

   ! Add filter loss to the downstream duct if no flow split loss over the two ducts
   IF (DN%FILTER_INDEX > 0) THEN
      IF (.NOT. DN%HMT_FILTER) CALL FILTER_UPDATE(DT,NE%NODE_INDEX(NN))
      IF (DUCT(DN%DUCT_INDEX(1))%AREA < TWO_EPSILON_EB .OR. DUCT(DN%DUCT_INDEX(2))%AREA < TWO_EPSILON_EB) CYCLE
      IF (FILTER(DN%FILTER_INDEX)%AREA > 0._EB) THEN
         AREA = FILTER(DN%FILTER_INDEX)%AREA
      ELSE
         AREA = 0.5_EB*(DUCT(DN%DUCT_INDEX(1))%AREA+DUCT(DN%DUCT_INDEX(2))%AREA)
      ENDIF
      IF(DUCT(DN%DUCT_INDEX(1))%VEL(PREVIOUS)*DN%DIR(1) > 0._EB) THEN
         DUCT(DN%DUCT_INDEX(2))%TOTAL_LOSS = DUCT(DN%DUCT_INDEX(2))%TOTAL_LOSS + &
                                             DN%FILTER_LOSS*(DUCT(DN%DUCT_INDEX(2))%AREA/AREA)**2
      ELSE
         DUCT(DN%DUCT_INDEX(1))%TOTAL_LOSS = DUCT(DN%DUCT_INDEX(1))%TOTAL_LOSS + &
                                             DN%FILTER_LOSS*(DUCT(DN%DUCT_INDEX(1))%AREA/AREA)**2
      ENDIF
   ENDIF

   NODECLASS: IF (DN%VENT .OR. DN%AMBIENT .OR. DN%GEOM) THEN
      ! If node is an external node loss is simply based on inflow or outflow or half loss if no flow
      IF (DUCT(DN%DUCT_INDEX(1))%AREA < TWO_EPSILON_EB .OR. DUCT(DN%DUCT_INDEX(1))%LOCALIZED_LEAKAGE) CYCLE
      IF(DUCT(DN%DUCT_INDEX(1))%VEL(PREVIOUS)*DN%DIR(1) > 1.E-6_EB) THEN
         DUCT(DN%DUCT_INDEX(1))%TOTAL_LOSS = DUCT(DN%DUCT_INDEX(1))%TOTAL_LOSS + DN%LOSS_ARRAY(1,2)
      ELSEIF(DUCT(DN%DUCT_INDEX(1))%VEL(PREVIOUS)*DN%DIR(1) < -1.E-6_EB) THEN
         DUCT(DN%DUCT_INDEX(1))%TOTAL_LOSS = DUCT(DN%DUCT_INDEX(1))%TOTAL_LOSS + DN%LOSS_ARRAY(2,1)
      ELSE
         DUCT(DN%DUCT_INDEX(1))%TOTAL_LOSS = DUCT(DN%DUCT_INDEX(1))%TOTAL_LOSS + 0.5_EB*(DN%LOSS_ARRAY(1,2)+DN%LOSS_ARRAY(2,1))
      ENDIF

   ELSEIF (DN%N_DUCTS==2) THEN NODECLASS
      IF (DUCT(DN%DUCT_INDEX(1))%AREA < TWO_EPSILON_EB .OR. DUCT(DN%DUCT_INDEX(2))%AREA < TWO_EPSILON_EB) CYCLE
      IF(DUCT(DN%DUCT_INDEX(1))%VEL(PREVIOUS)*DN%DIR(1) > 1.E-6_EB) THEN
         DUCT(DN%DUCT_INDEX(2))%TOTAL_LOSS = DUCT(DN%DUCT_INDEX(2))%TOTAL_LOSS + DN%LOSS_ARRAY(1,2)
      ELSEIF (DUCT(DN%DUCT_INDEX(1))%VEL(PREVIOUS)*DN%DIR(1) < -1.E-6_EB) THEN
         DUCT(DN%DUCT_INDEX(1))%TOTAL_LOSS = DUCT(DN%DUCT_INDEX(1))%TOTAL_LOSS + DN%LOSS_ARRAY(2,1)
      ELSE
         DUCT(DN%DUCT_INDEX(1))%TOTAL_LOSS = DUCT(DN%DUCT_INDEX(1))%TOTAL_LOSS + 0.5_EB*DN%LOSS_ARRAY(2,1)
         DUCT(DN%DUCT_INDEX(2))%TOTAL_LOSS = DUCT(DN%DUCT_INDEX(2))%TOTAL_LOSS + 0.5_EB*DN%LOSS_ARRAY(1,2)
      ENDIF

   ELSE NODECLASS
      ! For an internal node each outlet is weights the inlet flows
      NUM_OUT = 0
      DO ND=1,DN%N_DUCTS
         DU => DUCT(DN%DUCT_INDEX(ND))
         IF (DU%AREA < TWO_EPSILON_EB) CYCLE
         IF (DU%VEL(PREVIOUS)*DN%DIR(ND) < 0._EB) NUM_OUT = NUM_OUT + 1
      ENDDO

      NUM_OUT_IF: IF (NUM_OUT==0 .OR. NUM_OUT==DN%N_DUCTS) THEN
         ! If all are inflow or outflow each duct gets the average of its inflowing losses normalized by the number of ducts
         DO ND=1,DN%N_DUCTS
            DU => DUCT(DN%DUCT_INDEX(ND))
            LOSS_SUM = 0._EB
            DO ND2=1,DN%N_DUCTS
               IF (ND2==ND) CYCLE
               LOSS_SUM = LOSS_SUM + DN%LOSS_ARRAY(ND2,ND)
            ENDDO
            DU%TOTAL_LOSS = DU%TOTAL_LOSS + LOSS_SUM / (DN%N_DUCTS-1) / DN%N_DUCTS
         ENDDO

      ELSE NUM_OUT_IF
         ! Weight the outflowing losses based on fraction of inflow volume flow
         VFLOW = 0._EB
         DO ND=1,DN%N_DUCTS
            DU => DUCT(DN%DUCT_INDEX(ND))
            IF (DU%VEL(PREVIOUS)*DN%DIR(ND) > 0._EB) VFLOW = VFLOW + ABS(DU%VEL(PREVIOUS)*DU%AREA)
         ENDDO
         DO ND=1,DN%N_DUCTS
            DU => DUCT(DN%DUCT_INDEX(ND))
            IF (DU%VEL(PREVIOUS)*DN%DIR(ND) > 0._EB) CYCLE
            DO ND2=1,DN%N_DUCTS
               DU2 => DUCT(DN%DUCT_INDEX(ND2))
               IF (ND == ND2 .OR. DU2%VEL(PREVIOUS)*DN%DIR(ND2) <=0._EB) CYCLE
               DU%TOTAL_LOSS = DU%TOTAL_LOSS + DN%LOSS_ARRAY(ND2,ND) * ABS(DU2%VEL(PREVIOUS)*DU2%AREA)/VFLOW
            ENDDO
         ENDDO
      ENDIF NUM_OUT_IF
   ENDIF NODECLASS
ENDDO NODELOOP

DO ND = 1, NE%N_DUCTS
   DU => DUCT(NE%DUCT_INDEX(ND))
   IF (DU%ROUGHNESS > 0._EB) THEN
      ZZ_GET(1:N_TRACKED_SPECIES) = DU%ZZ(1:N_TRACKED_SPECIES)
      CALL GET_VISCOSITY(ZZ_GET,VISCOSITY,DU%TMP_D)
      FRICTION_FACTOR = COMPUTE_FRICTION_FACTOR(DU%RHO_D,VISCOSITY,ABS(DU%VEL(PREVIOUS)),DU%DIAMETER,DU%ROUGHNESS)
   ELSE
      FRICTION_FACTOR = 0._EB
   ENDIF
   IF (DU%VEL(PREVIOUS)>0._EB) THEN
      LOSS_SUM = DU%LOSS(1) * EVALUATE_RAMP(T,DU%RAMP_LOSS_INDEX)
   ELSEIF (DU%VEL(PREVIOUS)<0._EB) THEN
      LOSS_SUM = DU%LOSS(2) * EVALUATE_RAMP(T,DU%RAMP_LOSS_INDEX)
   ELSE
      LOSS_SUM = 0.5_EB*(DU%LOSS(1)+DU%LOSS(2)) * EVALUATE_RAMP(T,DU%RAMP_LOSS_INDEX)
   ENDIF
   DU%TOTAL_LOSS = DU%TOTAL_LOSS + DU%LENGTH/DU%DIAMETER*FRICTION_FACTOR + LOSS_SUM
   IF (DU%FAN_INDEX>0) THEN
      IF(.NOT. DU%FAN_OPERATING) DU%TOTAL_LOSS = DU%TOTAL_LOSS + FAN(DU%FAN_INDEX)%OFF_LOSS
   ENDIF
ENDDO

END SUBROUTINE UPDATE_LOSS

!> \brief Calculates the friction factor for a duct
!>
!> \param RHO Gas density in duct (kg/m3)
!> \param VISCOSITY Gas viscosity in duct (kg/m/s)
!> \param VEL Duct velocity (m/s)
!> \param DIAM Duct diameter (m)
!> \param ROUGHNESS Duct absolute roughness (m)

REAL(EB) FUNCTION COMPUTE_FRICTION_FACTOR(RHO,VISCOSITY,VEL,DIAM,ROUGHNESS)
REAL(EB), INTENT(IN) :: RHO,VISCOSITY,VEL,DIAM,ROUGHNESS
REAL(EB) :: EOD,RE_D

RE_D = MAX(100._EB,RHO*DIAM*VEL/VISCOSITY)
EOD = ROUGHNESS / DIAM
COMPUTE_FRICTION_FACTOR = LOG10(6.9_EB/RE_D+(EOD/3.7_EB)**1.11_EB)
COMPUTE_FRICTION_FACTOR = LOG10(EOD/3.7_EB-4.518_EB/RE_D*COMPUTE_FRICTION_FACTOR)
COMPUTE_FRICTION_FACTOR = (-0.5_EB/COMPUTE_FRICTION_FACTOR)**2

RETURN

END FUNCTION COMPUTE_FRICTION_FACTOR

!> \brief Updates duct velocity previous and guess values and calls UPDATE_FAN
!>
!> \param T Current time (s)
!> \param NNE Index indicating which HVAC network is being solved

SUBROUTINE SET_GUESS(NNE,T)
INTEGER, INTENT(IN) :: NNE
REAL(EB), INTENT(IN):: T
INTEGER :: ND
TYPE(DUCT_TYPE),POINTER :: DU
TYPE(NETWORK_TYPE), POINTER :: NE

NE => NETWORK(NNE)

DO ND = 1,NE%N_DUCTS
   DU => DUCT(NE%DUCT_INDEX(ND))
   IF (DU%FAN_INDEX > 0 .AND. DU%FAN_OPERATING) THEN
      CALL UPDATE_FAN(T,NE%DUCT_INDEX(ND))
   ELSE
      DU%DP_FAN = 0._EB
   ENDIF
   IF (DU%FIXED) THEN
      DU%VEL(PREVIOUS)  = DU%VEL(NEW)
      DU%VEL(GUESS)     = DU%VEL(NEW)
   ELSE
      IF (SIGN(1._EB,DU%VEL(NEW))==SIGN(1._EB,DU%VEL(PREVIOUS))) THEN
         DU%VEL(GUESS) = DU%VEL(NEW)
      ELSE
         DU%VEL(GUESS) = 0._EB
      ENDIF
      DU%VEL(PREVIOUS) = DU%VEL(NEW)
   ENDIF
ENDDO

END SUBROUTINE SET_GUESS

!> \brief sets donor (upstream) values for ducts and nodes
!> \param NNE Index indicating which HVAC network is being solved

SUBROUTINE SET_DONOR(NNE)

INTEGER :: ND,NN
INTEGER, INTENT(IN) :: NNE
REAL(EB) :: RHOLAST,TMPLAST,FVAL,OMFVAL,ITERFRAC
TYPE(DUCT_TYPE), POINTER :: DU
TYPE(DUCTNODE_TYPE), POINTER :: DN
TYPE(NETWORK_TYPE), POINTER :: NE

NE => NETWORK(NNE)
ITERFRAC = REAL(ITER,EB)/REAL(ITER_MAX,EB)
FVAL = MIN(1._EB,MAX(0._EB,(ITERFRAC-ONTH))/ONTH,1._EB)
OMFVAL = 1._EB - FVAL

NODELOOP: DO NN=1,NE%N_DUCTNODES
   DN=>DUCTNODE(NE%NODE_INDEX(NN))
   IF(DN%VENT .OR. DN%AMBIENT .OR. DN%LEAKAGE .OR. DN%GEOM) THEN
      DN%RHO= DN%RHO_V
      DN%TMP= DN%TMP_V
      DN%ZZ = DN%ZZ_V
      DN%RSUM = DN%RSUM_V
      DN%CP = DN%CP_V
   ENDIF
ENDDO NODELOOP

DUCTLOOP: DO ND=1,NE%N_DUCTS
   DU=>DUCT(NE%DUCT_INDEX(ND))
   IF (DU%AREA < TWO_EPSILON_EB) CYCLE DUCTLOOP
   RHOLAST = DU%RHO_D
   TMPLAST = DU%TMP_D

   IF (ABS(DU%VEL(PREVIOUS))>0._EB) THEN
      IF (DU%VEL(PREVIOUS)>0._EB) THEN
         DN=>DUCTNODE(DU%NODE_INDEX(1))
      ELSE
         DN=>DUCTNODE(DU%NODE_INDEX(2))
      ENDIF
      DU%CP_D = DN%CP
      DU%RHO_D = DN%RHO
      DU%TMP_D = DN%TMP
      DU%RSUM_D = DN%RSUM
   ELSE
      DU%CP_D = 0.5_EB*(DUCTNODE(DU%NODE_INDEX(1))%CP+DUCTNODE(DU%NODE_INDEX(2))%CP)
      DU%RHO_D = 0.5_EB*(DUCTNODE(DU%NODE_INDEX(1))%RHO+DUCTNODE(DU%NODE_INDEX(2))%RHO)
      DU%TMP_D = 0.5_EB*(DUCTNODE(DU%NODE_INDEX(1))%TMP+DUCTNODE(DU%NODE_INDEX(2))%TMP)
      DU%RSUM_D = 0.5_EB*(DUCTNODE(DU%NODE_INDEX(1))%RSUM+DUCTNODE(DU%NODE_INDEX(2))%RSUM)
   ENDIF
   IF (ITERFRAC > ONTH) THEN
      DU%RHO_D = FVAL*RHOLAST + OMFVAL*DU%RHO_D
      DU%TMP_D = FVAL*TMPLAST + OMFVAL*DU%TMP_D
   ENDIF

ENDDO DUCTLOOP

END SUBROUTINE SET_DONOR

!> \brief Checks the current iteration for duct velocity convergence and conservation of mass at nodes
!>
!> \param NNE Index indicating which HVAC network is being solved

SUBROUTINE CONVERGENCE_CHECK(NNE)
INTEGER, INTENT(IN) :: NNE
INTEGER :: NN, ND, COUNT
LOGICAL :: CONVERGED
REAL(EB) :: MSUM,MTOT,MFLOW
TYPE(NETWORK_TYPE), POINTER :: NE
TYPE(DUCT_TYPE), POINTER :: DU
TYPE(DUCTNODE_TYPE), POINTER :: DN

NE => NETWORK(NNE)
CONVERGED = .TRUE.
! Check duct velocity convergence
DO ND=1,NE%N_DUCTS
   DU => DUCT(NE%DUCT_INDEX(ND))
   IF (DU%AREA < TWO_EPSILON_EB) CYCLE
   IF (ABS(DU%VEL(PREVIOUS)) < 1.E-5_EB .AND. ABS(DU%VEL(NEW)) < 1.E-5_EB) CYCLE
   IF (DU%VEL(PREVIOUS) < 0._EB .EQV. DU%VEL(NEW) < 0._EB) THEN ! check for flow reversal
         IF (ABS(DU%VEL(PREVIOUS))<=TWO_EPSILON_EB) THEN
            CONVERGED = .FALSE.
            EXIT
         ELSE
            IF (ABS(1._EB-DU%VEL(NEW)/DU%VEL(PREVIOUS)) > 0.05_EB) THEN
            CONVERGED = .FALSE.
            CYCLE
            ENDIF
         ENDIF
   ELSE
      CONVERGED = .FALSE.
   ENDIF
ENDDO

IF (.NOT. CONVERGED) RETURN ! if velocity convergence passes, continue

! Check node mass conservation convergence
DO NN=1,NE%N_DUCTNODES
   DN => DUCTNODE(NE%NODE_INDEX(NN))
   IF (DN%VENT .OR. DN%AMBIENT .OR. DN%LEAKAGE .OR. DN%GEOM) CYCLE
   MSUM = 0._EB
   MTOT = 0._EB
   COUNT = 0
   DO ND=1,DN%N_DUCTS
      DU=>DUCT(DN%DUCT_INDEX(ND))
      IF (ABS(DU%VEL(NEW))<1.E-5_EB) COUNT = COUNT + 1
      MFLOW = DN%DIR(ND)*DU%VEL(NEW)*DU%RHO_D*DU%AREA
      MSUM = MSUM + MFLOW
      MTOT = MTOT + ABS(MFLOW)
   ENDDO
   IF (DN%FILTER_INDEX > 0) THEN
      MFLOW = SUM(DN%FILTER_LOADING(:,3))
      MSUM = MSUM - MFLOW
      MTOT = MTOT + MFLOW
   ENDIF
   IF(ABS(MSUM)< 1.E-6 * MTOT .OR. MTOT < TWO_EPSILON_EB .OR. COUNT >= DN%N_DUCTS-1) CYCLE
   CONVERGED = .FALSE.
ENDDO

IF (CONVERGED) ITER=ITER_MAX+1

END SUBROUTINE CONVERGENCE_CHECK

!> \brief Combines the mesh based arrays of HVAC boundary condtions to determine the HVAC solver boundary conditions
!>
!> \param T Current time (s)

SUBROUTINE COLLAPSE_HVAC_BC(T)

! Takes the MPI gathered mesh array of HVAC boundary conditions and updates the DUCTNODE boundary condition values.
USE PHYSICAL_FUNCTIONS, ONLY : GET_SPECIFIC_GAS_CONSTANT,GET_ENTHALPY,GET_TEMPERATURE
USE MATH_FUNCTIONS, ONLY : EVALUATE_RAMP
REAL(EB), INTENT(IN) :: T
REAL(EB) :: TNOW !< Current CPU time (s) used in computing length of time spent in HVAC routines.
INTEGER:: NN,NS,ND
REAL(EB) :: AREA,RHO_SUM,ZZ_GET(1:N_TRACKED_SPECIES),HGAS,P_OFFSET
TYPE(DUCTNODE_TYPE), POINTER :: DN,DN2

TNOW = CURRENT_TIME()

VENT_CUSTOM_AMBIENT: DO NN=1,N_DUCTNODES
   DN => DUCTNODE(NN)
   IF (.NOT. DN%LEAKAGE) THEN
      DN%ZONE_INDEX = NODE_ZONE(NN)
      IF (CONNECTED_ZONES(0,DN%ZONE_INDEX)>0) DN%ZONE_INDEX = 0
   ENDIF

   INTERNAL_NODE_IF: IF (((DN%VENT .OR. DN%LEAKAGE .OR. DN%GEOM) .AND. .NOT. DN%AMBIENT) .OR. &
                          (DN%AMBIENT .AND. NODE_AREA(NN) > 0._EB)) THEN
      ZZ_GET = 0._EB
      AREA = NODE_AREA(NN)
      IF (AREA<=TWO_EPSILON_EB) THEN
         DUCT(DN%DUCT_INDEX(1))%AREA = 0._EB
         DUCT(DN%DUCT_INDEX(1))%VEL = 0._EB
         CYCLE VENT_CUSTOM_AMBIENT
      ELSE
         DUCT(DN%DUCT_INDEX(1))%AREA = DUCT(DN%DUCT_INDEX(1))%AREA_INITIAL
      ENDIF
      NODE_AREA_EX(NN) = AREA
      DN%XYZ(1) = NODE_X(NN)/AREA
      DN%XYZ(2) = NODE_Y(NN)/AREA
      DN%XYZ(3) = NODE_Z(NN)/AREA
      RHO_SUM = NODE_RHO(NN)

      DO NS=1,N_TRACKED_SPECIES
         DN%ZZ_V(NS) = NODE_ZZ(NN,NS)/AREA
      ENDDO


      ZZ_GET(1:N_TRACKED_SPECIES) = DN%ZZ_V(1:N_TRACKED_SPECIES)
      CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,DN%RSUM_V)

      DN%RHO_V = AREA/RHO_SUM

      !Initialize default values
      IF (DN%P < -1.E9_EB) THEN
         IF (STRATIFICATION) THEN
            DN%TMP = TMPA*EVALUATE_RAMP(DN%XYZ(3),I_RAMP_TMP0_Z)
            DN%P_OLD   = EVALUATE_RAMP(DN%XYZ(3),I_RAMP_P0_Z)
            DN%RHO = DN%P_OLD/(DN%TMP*RSUM0)
         ELSE
            DN%TMP = TMPA
            DN%P   = P_INF
            DN%RHO = RHOA
         ENDIF
         CALL GET_ENTHALPY(ZZ_GET,HGAS,DN%TMP)
         DN%CP = HGAS / DN%TMP
      ENDIF

      IF (DN%LEAKAGE) THEN
         IF(ABS(T-T_BEGIN) < TWO_EPSILON_EB) DN%P_OLD = NODE_P(NN)/AREA+P_INF
      ENDIF

      DN%P = HVAC_PRES_RELAX*(NODE_P(NN)/AREA+P_INF)+(1._EB-HVAC_PRES_RELAX)*DN%P_OLD
      DN%TMP_V = NODE_TMP(NN)/AREA
      HGAS = NODE_H(NN)/AREA
      CALL GET_TEMPERATURE(DN%TMP_V,HGAS,ZZ_GET)
      CALL GET_ENTHALPY(ZZ_GET,HGAS,DN%TMP_V)
      DN%CP_V = HGAS/DN%TMP_V
   ENDIF INTERNAL_NODE_IF
END DO VENT_CUSTOM_AMBIENT

AMBIENT_LEAK: DO NN=1,N_DUCTNODES
   DN => DUCTNODE(NN)
    IF (DN%AMBIENT .AND. NODE_AREA(NN)<=TWO_EPSILON_EB) THEN
      !Initialize ambient nodes outside domain
      IF (DUCT(DN%DUCT_INDEX(1))%NODE_INDEX(1)==NN) THEN
         DN2 => DUCTNODE(DUCT(DN%DUCT_INDEX(1))%NODE_INDEX(2))
      ELSE
         DN2 => DUCTNODE(DUCT(DN%DUCT_INDEX(1))%NODE_INDEX(1))
      ENDIF
      IF (DUCT(DN%DUCT_INDEX(1))%LOCALIZED_LEAKAGE .OR. DN%LEAKAGE .OR. DN%XYZ(3) <-1.E9_EB) DN%XYZ = DN2%XYZ
      IF (DN%XYZ(3) <-1.E9_EB) CYCLE AMBIENT_LEAK
      DN%RSUM   = RSUM0
      IF (STRATIFICATION) THEN
         DN%TMP = TMPA*EVALUATE_RAMP(DN%XYZ(3),I_RAMP_TMP0_Z)
         DN%P   = EVALUATE_RAMP(DN%XYZ(3),I_RAMP_P0_Z)
         DN%RHO = DN%P/(DN%TMP*RSUM0)
      ELSE
         DN%TMP = TMPA
         DN%P   = P_INF
         DN%RHO = RHOA
      ENDIF
      DN%P_OLD = DN%P
      ZZ_GET(1:N_TRACKED_SPECIES) = DN%ZZ_V(1:N_TRACKED_SPECIES)
      CALL GET_ENTHALPY(ZZ_GET,HGAS,DN%TMP)
      DN%TMP_V  = DN%TMP
      DN%RSUM_V = DN%RSUM
      DN%CP     = HGAS/DN%TMP
      DN%CP_V   = DN%CP
      DN%RHO_V  = DN%RHO
   ENDIF
ENDDO AMBIENT_LEAK

SAME_Z_LEAK: DO ND=1,N_DUCTS
   IF (.NOT. DUCT(ND)%LEAKAGE) CYCLE SAME_Z_LEAK
   DN=> DUCTNODE(DUCT(ND)%NODE_INDEX(1))
   DN2=> DUCTNODE(DUCT(ND)%NODE_INDEX(2))
   P_OFFSET = 0.5_EB*(EVALUATE_RAMP(DN%XYZ(3),I_RAMP_P0_Z) - EVALUATE_RAMP(DN2%XYZ(3),I_RAMP_P0_Z))
   DN%P = DN%P - P_OFFSET
   DN2%P = DN2%P + P_OFFSET
ENDDO SAME_Z_LEAK

T_USED(13)=T_USED(13)+CURRENT_TIME()-TNOW

END SUBROUTINE COLLAPSE_HVAC_BC

!> \brief Sets current iteration values for duct nodes connected to a vent to the bondary condition value

SUBROUTINE SET_INIT_HVAC
REAL(EB) :: TNOW !< Current CPU time (s) used in computing length of time spent in HVAC routines.
REAL(EB) :: XYZ1(3), XYZ2(3)
INTEGER:: ND,NN,NW,NC
TYPE(DUCTNODE_TYPE), POINTER :: DN,DN2
TYPE(DUCT_TYPE), POINTER :: DU

TNOW = CURRENT_TIME()

DO NN=1,N_DUCTNODES
   DN=>DUCTNODE(NN)
   IF (.NOT. (DN%VENT .OR. DN%LEAKAGE .OR. DN%GEOM)) CYCLE
   DN%CP = DN%CP_V
   DN%TMP = DN%TMP_V
   DN%RHO = DN%RHO_V
   DN%RSUM = DN%RSUM_V
   DN%ZZ = DN%ZZ_V
ENDDO

DO ND=1,N_DUCTS
   DU => DUCT(ND)
   IF (DU%LENGTH < 0._EB) THEN
      IF (.NOT. DUCTNODE(DU%NODE_INDEX(1))%SPECIFIED_XYZ .OR. .NOT. DUCTNODE(DU%NODE_INDEX(2))%SPECIFIED_XYZ) THEN
         WRITE(MESSAGE,'(A,A)') 'ERROR(554): Duct has no LENGTH and one node lacks an XYZ. Duct: ',TRIM(DU%ID)
         CALL SHUTDOWN(MESSAGE); RETURN
      ENDIF
      DN =>  DUCTNODE(DU%NODE_INDEX(1))
      DN2 => DUCTNODE(DU%NODE_INDEX(2))
      IF (DU%N_WAYPOINTS == 0) THEN
         DU%LENGTH = SQRT((DN%XYZ(1)-DN2%XYZ(1))**2+(DN%XYZ(2)-DN2%XYZ(2))**2+(DN%XYZ(3)-DN2%XYZ(3))**2)
      ELSE
         DU%LENGTH = 0._EB
         DO NW=1,DU%N_WAYPOINTS+1
            IF (NW==1) THEN
               XYZ1 = DN%XYZ
            ELSE
               XYZ1 = XYZ2
            ENDIF
            IF (NW==DU%N_WAYPOINTS+1) THEN
               XYZ2 = DN2%XYZ
            ELSE
               XYZ2(:) = DU%WAYPOINT_XYZ(NW,:)
            ENDIF
            DU%LENGTH = DU%LENGTH + SQRT((XYZ1(1)-XYZ2(1))**2+(XYZ1(2)-XYZ2(2))**2+(XYZ1(3)-XYZ2(3))**2)
         ENDDO
      ENDIF
   ENDIF

   IF (HVAC_MASS_TRANSPORT_CELL_L > 0._EB .AND. .NOT. (DU%LEAKAGE .OR. DU%LOCALIZED_LEAKAGE)) &
      DU%N_CELLS = MAX(1,INT(DU%LENGTH/HVAC_MASS_TRANSPORT_CELL_L))
   IF (DU%N_CELLS > 0) THEN
      DU%DX = DU%LENGTH/DU%N_CELLS
      ALLOCATE(DU%RHO_C(DU%N_CELLS))
      ALLOCATE(DU%TMP_C(DU%N_CELLS))
      ALLOCATE(DU%CP_C(DU%N_CELLS))
      ALLOCATE(DU%ZZ_C(DU%N_CELLS,N_TRACKED_SPECIES))
      ! Initialising as background/ambient here; required for DEVC output at t = 0 s
      DU%RHO_C = RHOA
      DU%TMP_C = TMPA
      DO NC = 1,DU%N_CELLS
         DU%ZZ_C(NC,1:N_TRACKED_SPECIES) = SPECIES_MIXTURE(1:N_TRACKED_SPECIES)%ZZ0
      ENDDO
      ALLOCATE(DU%RHO_C_OLD(DU%N_CELLS))
      ALLOCATE(DU%TMP_C_OLD(DU%N_CELLS))
      ALLOCATE(DU%CP_C_OLD(DU%N_CELLS))
      ALLOCATE(DU%ZZ_C_OLD(DU%N_CELLS,N_TRACKED_SPECIES))
      DO NC = 1,DU%N_CELLS ! Initialising as background here; required for DEVC output at t = 0 s
         DU%ZZ_C_OLD(NC,1:N_TRACKED_SPECIES) = SPECIES_MIXTURE(1:N_TRACKED_SPECIES)%ZZ0
      ENDDO
   ENDIF
ENDDO

IF (N_ZONE>0) ALLOCATE(PSUM_TOT(N_ZONE))

T_USED(13)=T_USED(13)+CURRENT_TIME()-TNOW

END SUBROUTINE SET_INIT_HVAC

!> \brief Initializes cell densities, temperatures, specific heats and species for discretized ducts

SUBROUTINE SET_INIT_HVAC_MASS_TRANSPORT
USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
USE PHYSICAL_FUNCTIONS, ONLY: GET_ENTHALPY,GET_SPECIFIC_GAS_CONSTANT
REAL(EB) :: P1,P2,T1,T2,ZZ1(1:N_TRACKED_SPECIES),ZZ2(1:N_TRACKED_SPECIES),FAC,RSUM
INTEGER :: ND,NC
REAL(EB) :: ZZ_GET(1:N_TRACKED_SPECIES),HGAS
TYPE(DUCT_TYPE), POINTER :: DU
TYPE(DUCTNODE_TYPE), POINTER :: DN

ZZ_GET(1:N_TRACKED_SPECIES) = SPECIES_MIXTURE(1:N_TRACKED_SPECIES)%ZZ0

! Initialize ductnode values and duct cell arrays for duct mass transport
DUCT_LOOP:DO ND = 1, N_DUCTS
   DU => DUCT(ND)
   IF (DU%LEAKAGE .OR. DU%LOCALIZED_LEAKAGE .OR. DU%N_CELLS==0) CYCLE
   DN => DUCTNODE(DU%NODE_INDEX(1))
   P1 = DN%P
   T1 = DN%TMP0
   ZZ1 = DN%ZZ0
   DN => DUCTNODE(DU%NODE_INDEX(2))
   P2 = DN%P
   T2 = DN%TMP0
   ZZ2 = DN%ZZ0
   DO NC=1,DU%N_CELLS
      FAC = (REAL(NC,EB)-0.5_EB)/REAL(DU%N_CELLS)
      DU%ZZ_C(NC,:) = ZZ1(:) + (ZZ2(:)-ZZ1(:)) * FAC
      ZZ_GET = DU%ZZ_C(NC,:)
      DU%TMP_C(NC) = T1 + (T2-T1) * FAC
      CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,RSUM)
      DU%RHO_C(NC) = (P1 + (P2-P1) * FAC)/(RSUM*DU%TMP_C(NC))
      CALL GET_ENTHALPY(ZZ_GET,HGAS,DU%TMP_C(NC))
      DU%CP_C(NC) = HGAS/DU%TMP_C(NC)
   ENDDO
   DU%ZZ_C_OLD = DU%ZZ_C
   DU%TMP_C_OLD = DU%TMP_C
   DU%RHO_C_OLD = DU%RHO_C
   DU%CP_C_OLD = DU%CP_C
ENDDO DUCT_LOOP

END SUBROUTINE SET_INIT_HVAC_MASS_TRANSPORT

!> \brief Initializes leakage ducts and nodes

SUBROUTINE LEAKAGE_HVAC

USE PHYSICAL_FUNCTIONS, ONLY: GET_ENTHALPY
REAL(EB) :: ZZ_GET(1:N_TRACKED_SPECIES),HGAS,TNOW
INTEGER :: I_DUCT,I_DUCTNODE,NZ1,NZ2
TYPE (DUCTNODE_TYPE), POINTER:: DN1,DN2
TYPE (DUCT_TYPE), POINTER:: DU

IF (LEAK_DUCTS==0) RETURN

TNOW=CURRENT_TIME()

I_DUCT = N_DUCTS - LEAK_DUCTS
I_DUCTNODE = N_DUCTNODES - 2 * LEAK_DUCTS

DO NZ1 = 0, N_ZONE
   DO NZ2 = 0, N_ZONE
      IF (LEAK_PATH(NZ1,NZ2) == 1) THEN
         I_DUCT = I_DUCT + 1
         LEAK_PATH(NZ1,NZ2) = I_DUCT
         I_DUCTNODE = I_DUCTNODE + 1
         DU => DUCT(I_DUCT)
         DU%AREA_INITIAL = P_ZONE(NZ2)%LEAK_AREA(NZ1)
         DU%AREA = P_ZONE(NZ2)%LEAK_AREA(NZ1)
         DU%DIAMETER = SQRT(DU%AREA / PI)*2._EB
         DU%LEAKAGE = .TRUE.
         WRITE(DU%ID,'(A,1X,I0,1X,I0)') 'LEAK',NZ1,NZ2
         DU%LENGTH = 0.1_EB
         DU%LOSS(1) = 1._EB
         DU%LOSS(2) = 1._EB
         DU%NODE_INDEX(1) = I_DUCTNODE
         DU%NODE_INDEX(2) = I_DUCTNODE+1
         DU%RHO_D = RHOA
         DU%TMP_D = TMPA
         DU%ROUGHNESS = -1._EB
         DU%VEL = 0._EB
         ALLOCATE(DU%ZZ(N_TRACKED_SPECIES))
         ZZ_GET(1:N_TRACKED_SPECIES) = SPECIES_MIXTURE(1:N_TRACKED_SPECIES)%ZZ0
         DU%ZZ(1:N_TRACKED_SPECIES) = ZZ_GET(1:N_TRACKED_SPECIES)
         CALL GET_ENTHALPY(ZZ_GET,HGAS,TMPA)
         DU%CP_D = HGAS / TMPA
         DN1=>DUCTNODE(I_DUCTNODE)
         IF (NZ1==0) DN1%AMBIENT = .TRUE.
         ALLOCATE(DN1%DUCT_INDEX(1))
         ALLOCATE(DN1%DIR(1))
         ALLOCATE(DN1%LOSS_ARRAY(2,2))
         DN1%LOSS_ARRAY = 0._EB
         DN1%DIR = -1._EB
         DN1%DUCT_INDEX = I_DUCT
         DN1%LEAKAGE = .TRUE.
         DN1%ZONE_INDEX=NZ1
         DN1%N_DUCTS = 1
         DN1%RSUM = RSUM0
         DN1%TMP = TMPA
         DN1%XYZ = (/0._EB,0._EB,0._EB/)
         WRITE(DN1%ID,'(A,1X,I0,1X,I0)') 'LEAK',NZ1,NZ2
         I_DUCTNODE = I_DUCTNODE + 1
         DN2=>DUCTNODE(I_DUCTNODE)
         IF (NZ2==0) DN2%AMBIENT = .TRUE.
         ALLOCATE(DN2%DUCT_INDEX(1))
         ALLOCATE(DN2%DIR(1))
         ALLOCATE(DN2%LOSS_ARRAY(2,2))
         DN2%LOSS_ARRAY = 0._EB
         DN2%DIR = 1._EB
         DN2%DUCT_INDEX = I_DUCT
         DN2%LEAKAGE = .TRUE.
         DN2%ZONE_INDEX=NZ2
         DN2%XYZ = (/0._EB,0._EB,0._EB/)
         DN2%N_DUCTS = 1
         DN2%RSUM = RSUM0
         DN2%TMP = TMPA
         WRITE(DN2%ID,'(A,1X,I0,1X,I0)') 'LEAK',NZ2,NZ1
      ENDIF
   ENDDO
ENDDO

T_USED(13)=T_USED(13)+CURRENT_TIME()-TNOW

END SUBROUTINE LEAKAGE_HVAC

!> \brief Updates filter loading and filter flow loss
!>
!> \param DT Current time step (s)
!> \param NODE_INDEX Index of a node containing a filter

SUBROUTINE FILTER_UPDATE(DT,NODE_INDEX)
!Updates filter loss and loading solution
USE MATH_FUNCTIONS,ONLY:EVALUATE_RAMP
INTEGER,INTENT(IN)::NODE_INDEX
REAL(EB), INTENT(IN) :: DT
REAL(EB) :: TOTAL_LOADING
TYPE(DUCTNODE_TYPE),POINTER::DN,DN2
TYPE(DUCT_TYPE),POINTER::DU
TYPE(FILTER_TYPE),POINTER:: FI

DN => DUCTNODE(NODE_INDEX)
FI => FILTER(DN%FILTER_INDEX)

TOTAL_LOADING = DOT_PRODUCT(FI%MULTIPLIER,DN%FILTER_LOADING(:,1))
IF (FI%RAMP_INDEX > 0) THEN
   DN%FILTER_LOSS = FI%CLEAN_LOSS + EVALUATE_RAMP(TOTAL_LOADING,FI%RAMP_INDEX)
ELSE
   DN%FILTER_LOSS = FI%CLEAN_LOSS + FI%LOADING_LOSS*TOTAL_LOADING
ENDIF

DU=>DUCT(DN%DUCT_INDEX(1))
IF (DU%VEL(PREVIOUS) >= 0._EB .AND. DU%NODE_INDEX(2)==NODE_INDEX) THEN
  DN2 => DUCTNODE(DU%NODE_INDEX(1))
ELSEIF (DU%VEL(PREVIOUS) <= 0._EB .AND. DU%NODE_INDEX(1)==NODE_INDEX) THEN
  DN2 => DUCTNODE(DU%NODE_INDEX(2))
ELSE
   DU=>DUCT(DN%DUCT_INDEX(2))
   IF (DU%VEL(PREVIOUS) >= 0._EB .AND. DU%NODE_INDEX(2)==NODE_INDEX) THEN
   DN2 => DUCTNODE(DU%NODE_INDEX(1))
   ELSE
   DN2 => DUCTNODE(DU%NODE_INDEX(2))
   ENDIF
ENDIF

!Ultimately add in logic for condensible gases
DN%FILTER_LOADING(:,3) = DU%AREA*ABS(DU%VEL(PREVIOUS))*DN2%RHO*DN2%ZZ*FI%EFFICIENCY
DN%FILTER_LOADING(:,2) = DN%FILTER_LOADING(:,1) + DN%FILTER_LOADING(:,3) * DT

END SUBROUTINE FILTER_UPDATE

!> \brief Updates the heat added or removed by aircoils
!>
!> \param T Current time (s)

SUBROUTINE COIL_UPDATE(T)
REAL(EB), INTENT(IN) :: T
INTEGER :: ND
TYPE(DUCT_TYPE),POINTER::DU

COIL_LOOP: DO ND = 1,N_DUCTS
   DU => DUCT(ND)
   IF (DU%AIRCOIL_INDEX < 0) CYCLE COIL_LOOP
   IF (DU%DEVC_INDEX > 0) THEN
      DU%COIL_OPERATING = DEVICE(DU%DEVC_INDEX)%CURRENT_STATE
      IF (DU%COIL_OPERATING .AND. (DEVICE(DU%DEVC_INDEX)%CURRENT_STATE .NEQV. DEVICE(DU%DEVC_INDEX)%PRIOR_STATE)) &
            DU%COIL_ON_TIME = T
   ENDIF
   IF (DU%CTRL_INDEX > 0) THEN
      DU%COIL_OPERATING = CONTROL(DU%CTRL_INDEX)%CURRENT_STATE
      IF (DU%COIL_OPERATING .AND. (CONTROL(DU%CTRL_INDEX)%CURRENT_STATE .NEQV. CONTROL(DU%CTRL_INDEX)%PRIOR_STATE)) &
            DU%COIL_ON_TIME = T
   ENDIF

   DU%COIL_Q = 0._EB
   IF (.NOT. DU%COIL_OPERATING) CYCLE COIL_LOOP
   IF (DU%AREA < TWO_EPSILON_EB) CYCLE COIL_LOOP

   CALL CALC_COIL_Q(DU%AIRCOIL_INDEX,DU%COIL_Q,T-DU%COIL_ON_TIME,DU%VEL(NEW),DU%TMP_D,DU%ZZ,DU%RHO_D,DU%CP_D,DU%AREA)

END DO COIL_LOOP

END SUBROUTINE COIL_UPDATE

!> \brief Calculates heat addition rate of an aircoil
!>
!> \param AC_INDEX Index for a specific aircoil
!> \param COIL_Q Heat addition rate (W)
!> \param T_ON Coil on time (s)
!> \param VEL_D Duct velocity (m/s)
!> \param TMP_D Duct temperature (K)
!> \param ZZ_D Duct tracked species mass fractions (kg/kg)
!> \param RHO_D Duct density (kg/m3)
!> \param VEL_D Duct velocity (m/s)
!> \param CP_D Duct specific heat (J/kg/K)
!> \param AREA_D Duct area (m2)

SUBROUTINE CALC_COIL_Q(AC_INDEX,COIL_Q,TSI,VEL_D,TMP_D,ZZ_D,RHO_D,CP_D,AREA_D)
USE MATH_FUNCTIONS, ONLY : EVALUATE_RAMP
USE PHYSICAL_FUNCTIONS, ONLY : GET_SPECIFIC_HEAT, GET_AVERAGE_SPECIFIC_HEAT, GET_ENTHALPY
INTEGER, INTENT(IN) :: AC_INDEX
REAL(EB), INTENT(IN) :: TSI,VEL_D,TMP_D,ZZ_D(1:N_TRACKED_SPECIES),RHO_D,CP_D,AREA_D
REAL(EB), INTENT(INOUT) :: COIL_Q
REAL(EB) :: TMP_OUT,TMP_GUESS,MDOT_DU,E_IN,MCP_C,CP1,CP2,DCPDT,HGAS
INTEGER :: ITER
TYPE(AIRCOIL_TYPE),POINTER:: AC

AC => AIRCOIL(AC_INDEX)

IF (AC%FIXED_Q > -1.E9_EB) THEN
   COIL_Q = AC%FIXED_Q*EVALUATE_RAMP(TSI,AC%RAMP_INDEX,TAU=AC%TAU)
   MDOT_DU = RHO_D*VEL_D*AREA_D
   CALL GET_SPECIFIC_HEAT(ZZ_D,CP1,TMP_D)
   TMP_GUESS = TMP_D+COIL_Q/(CP1*MDOT_DU)
   IF (INT(TMP_GUESS) >= I_MAX_TEMP) THEN
      COIL_Q = MDOT_DU*CP1*(REAL(I_MAX_TEMP,EB)-TMP_D)
   ELSEIF (TMP_GUESS < TWO_EPSILON_EB) THEN
      COIL_Q = -MDOT_DU*CP1*TMP_D
   ENDIF
ELSE
   ITER = 0
   TMP_GUESS = TMP_D
   MDOT_DU = RHO_D*ABS(VEL_D)*AREA_D
   IF (MDOT_DU < TWO_EPSILON_EB) RETURN
   MCP_C =  AC%COOLANT_MASS_FLOW*AC%COOLANT_SPECIFIC_HEAT
   E_IN = MDOT_DU*TMP_D*CP_D + MCP_C*AC%COOLANT_TEMPERATURE
   DO WHILE (ITER <= 10)
      CALL GET_ENTHALPY(ZZ_D,HGAS,TMP_GUESS)
      CALL GET_AVERAGE_SPECIFIC_HEAT(ZZ_D,CP1,TMP_GUESS)
      IF (TMP_GUESS > 1._EB) THEN
         CALL GET_AVERAGE_SPECIFIC_HEAT(ZZ_D,CP2,TMP_GUESS-1._EB)
         DCPDT = CP1-CP2
      ELSE
         CALL GET_AVERAGE_SPECIFIC_HEAT(ZZ_D,CP2,TMP_GUESS+1._EB)
         DCPDT = CP2-CP1
      ENDIF
      CP1 = HGAS / TMP_GUESS
      TMP_OUT = TMP_GUESS + (E_IN - MDOT_DU*HGAS - MCP_C * TMP_GUESS) / &
                              (MDOT_DU * (DCPDT * TMP_GUESS + CP1) + MCP_C)
      IF (ABS(TMP_OUT-TMP_GUESS) < TWO_EPSILON_EB) EXIT
      IF (ABS(TMP_OUT-TMP_GUESS)/TMP_GUESS < 0.0005) EXIT
      TMP_GUESS = TMP_OUT
      ITER = ITER + 1
   ENDDO

   COIL_Q = AC%COOLANT_MASS_FLOW * AC%COOLANT_SPECIFIC_HEAT*(AC%COOLANT_TEMPERATURE - TMP_OUT)*AC%EFFICIENCY

ENDIF

END SUBROUTINE CALC_COIL_Q

!> \brief Adjusts the leak area based on the zone-to-zone pressure difference

SUBROUTINE ADJUST_LEAKAGE_AREA

INTEGER :: ND
REAL(EB) :: N,DELTA_P,DELTA_P_REF,C_D
TYPE(DUCT_TYPE), POINTER :: DU

DO ND=1,N_DUCTS
   DU => DUCT(ND)
   IF (.NOT.DU%LEAKAGE .AND. .NOT.DU%LOCALIZED_LEAKAGE) CYCLE
   DELTA_P = ABS(DUCTNODE(DU%NODE_INDEX(1))%P-DUCTNODE(DU%NODE_INDEX(2))%P)
   IF (DU%LOCALIZED_LEAKAGE) THEN
      DELTA_P_REF = DU%LEAK_REFERENCE_PRESSURE
      N = DU%LEAK_PRESSURE_EXPONENT
      C_D = DU%DISCHARGE_COEFFICIENT
   ELSE
      DELTA_P_REF = P_ZONE(DUCTNODE(DU%NODE_INDEX(2))%ZONE_INDEX)%LEAK_REFERENCE_PRESSURE(DUCTNODE(DU%NODE_INDEX(1))%ZONE_INDEX)
      N = P_ZONE(DUCTNODE(DU%NODE_INDEX(2))%ZONE_INDEX)%LEAK_PRESSURE_EXPONENT(DUCTNODE(DU%NODE_INDEX(1))%ZONE_INDEX)
      C_D = P_ZONE(DUCTNODE(DU%NODE_INDEX(2))%ZONE_INDEX)%DISCHARGE_COEFFICIENT(DUCTNODE(DU%NODE_INDEX(1))%ZONE_INDEX)
   ENDIF
   DU%AREA = C_D * DU%AREA_INITIAL * (DELTA_P/DELTA_P_REF)**(N-0.5_EB)
ENDDO

END SUBROUTINE ADJUST_LEAKAGE_AREA

!> \brief Updates the 1D mass transport solution in ducts
!>
!> \param DT Current time step (s)

SUBROUTINE UPDATE_HVAC_MASS_TRANSPORT(DT,NR)
USE PHYSICAL_FUNCTIONS,ONLY: GET_TEMPERATURE, GET_ENTHALPY
REAL(EB), INTENT(IN) :: DT
INTEGER, INTENT(IN) :: NR
INTEGER :: N_SUBSTEPS,ND,NS,NC
TYPE(DUCT_TYPE),POINTER :: DU
REAL(EB) :: DT_CFL,DT_DUCT,MASS_FLUX,ZZ_GET(N_TRACKED_SPECIES),HGAS
REAL(EB), ALLOCATABLE, DIMENSION(:) :: CPT_C,CPT_F,RHOCPT_C
REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: RHOZZ_C,ZZ_F

DUCT_LOOP: DO ND = 1,DUCTRUN(NR)%N_DUCTS
   DU => DUCT(DUCTRUN(NR)%DUCT_INDEX(ND))
   IF (DU%N_CELLS == 0 ) CYCLE DUCT_LOOP
   ! Check for zero flow and zero area
   IF (ABS(DU%VEL(NEW))<=TWO_EPSILON_EB .OR. DU%AREA<=TWO_EPSILON_EB) CYCLE DUCT_LOOP

   MASS_FLUX = DU%RHO_D * DU%VEL(NEW)

   ! Set up of CFL and sub time step
   DT_CFL = ABS(DU%DX/(2*DU%VEL(NEW))) ! CFL for Godunov pure upwinding scheme
   N_SUBSTEPS = MAX(1,CEILING(DT/DT_CFL))
   DT_DUCT = DT/REAL(N_SUBSTEPS,EB)

   ! Set upwind face indices and allocate flux arrays
   ALLOCATE(ZZ_F(0:DU%N_CELLS,N_TRACKED_SPECIES))
   ALLOCATE(CPT_F(0:DU%N_CELLS))
   ALLOCATE(CPT_C(DU%N_CELLS))
   ALLOCATE(RHOCPT_C(DU%N_CELLS))
   ALLOCATE(RHOZZ_C(DU%N_CELLS,N_TRACKED_SPECIES))

   SUBSTEP_LOOP: DO NS = 1,N_SUBSTEPS

      ! Populates upwind face variables, accounting for direction of flow (i.e. includes relevant node value as first/last face)
      IF (DU%VEL(NEW)>0._EB) THEN
         ZZ_F(0,:) = DUCTNODE(DU%NODE_INDEX(1))%ZZ(:)
         CPT_F(0) = DUCTNODE(DU%NODE_INDEX(1))%CP*DUCTNODE(DU%NODE_INDEX(1))%TMP
         DO NC = 1,DU%N_CELLS
            ZZ_F(NC,:) = DU%ZZ_C(NC,:) ! Godunov upwinding
            CPT_F(NC) = DU%TMP_C(NC)*DU%CP_C(NC) ! Godunov upwinding
         ENDDO
      ELSE
         ZZ_F(DU%N_CELLS,:) = DUCTNODE(DU%NODE_INDEX(2))%ZZ(:)
         CPT_F(DU%N_CELLS) = DUCTNODE(DU%NODE_INDEX(2))%TMP*DUCTNODE(DU%NODE_INDEX(2))%CP
         DO NC = 0,DU%N_CELLS-1
            ZZ_F(NC,:) = DU%ZZ_C(NC+1,:)
            CPT_F(NC) = DU%TMP_C(NC+1)*DU%CP_C(NC+1)
         ENDDO
      ENDIF

      ! Compute discretized conservation equations using explicit Euler method with Godunov upwinding profile
      DO NC = 1,DU%N_CELLS
         RHOZZ_C(NC,:) = DU%RHO_C(NC)*DU%ZZ_C(NC,:)           - DT_DUCT / DU%DX * MASS_FLUX * ( ZZ_F(NC,:) - ZZ_F(NC-1,:) )
         RHOCPT_C(NC) =  DU%RHO_C(NC)*DU%TMP_C(NC)*DU%CP_C(NC) - DT_DUCT / DU%DX * MASS_FLUX * ( CPT_F(NC) - CPT_F(NC-1) )
      ENDDO

      ! Update cell centre variables
      DU_UPDATE_LOOP: DO NC = 1,DU%N_CELLS
         DU%RHO_C(NC) = SUM(RHOZZ_C(NC,1:N_TRACKED_SPECIES))
         DU%ZZ_C(NC,:) = RHOZZ_C(NC,:)/DU%RHO_C(NC)
         CPT_C(NC) = RHOCPT_C(NC)/DU%RHO_C(NC)
         ZZ_GET = DU%ZZ_C(NC,:) ! Single dimension to be used with GET_AVERAGE_...
         CALL GET_TEMPERATURE(DU%TMP_C(NC),CPT_C(NC),ZZ_GET)
         CALL GET_ENTHALPY(ZZ_GET,HGAS,DU%TMP_C(NC))
         DU%CP_C(NC) = HGAS / DU%TMP_C(NC)
      ENDDO DU_UPDATE_LOOP

   ENDDO SUBSTEP_LOOP

   DEALLOCATE(RHOZZ_C)
   DEALLOCATE(ZZ_F)
   DEALLOCATE(CPT_F)
   DEALLOCATE(CPT_C)
   DEALLOCATE(RHOCPT_C)
   
ENDDO DUCT_LOOP


END SUBROUTINE UPDATE_HVAC_MASS_TRANSPORT

!> \brief Checks duct network for minimal loss defintions

SUBROUTINE EXAMINE_LOSSES
INTEGER :: ND, ND2, ND3, NN, F_COUNTER
LOGICAL :: CHANGE
LOGICAL :: LOSS_D(N_DUCTS,2),LOSS_N(N_DUCTNODES),FIXED_D(N_DUCTS)
TYPE N_LOSS_TEMP
   LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: LOSS
END TYPE N_LOSS_TEMP
TYPE(N_LOSS_TEMP), ALLOCATABLE, DIMENSION(:) :: N_LOSS
TYPE(DUCT_TYPE), POINTER :: DU
TYPE(DUCTNODE_TYPE), POINTER :: DN,DN2

LOSS_D = .FALSE.
LOSS_N = .FALSE.

ALLOCATE (N_LOSS(N_DUCTNODES))

! Allocate array that identifies if a particular direction out of a node has a loss associated with it

DO NN=1,N_DUCTNODES
   DN=>DUCTNODE(NN)
   ! For terminal nodes we have two loss to track: (1,1) into the system and (1,2) out of the system
   IF (DN%N_DUCTS==1) THEN
      ALLOCATE(N_LOSS(NN)%LOSS(1,2))
      N_LOSS(NN)%LOSS = .FALSE.
      IF (DN%LOSS_ARRAY(1,2) > 0._EB) N_LOSS(NN)%LOSS(1,1) = .TRUE.
      IF (DN%LOSS_ARRAY(2,1) > 0._EB) N_LOSS(NN)%LOSS(1,2) = .TRUE.
   ELSE
      ALLOCATE(N_LOSS(NN)%LOSS(DN%N_DUCTS,DN%N_DUCTS))
      N_LOSS(NN)%LOSS = .FALSE.
      DO ND=1,DN%N_DUCTS
         N_LOSS(NN)%LOSS(ND,ND) = .TRUE.
         DO ND2=1,DN%N_DUCTS
            IF (DN%LOSS_ARRAY(ND,ND2) > 0._EB) N_LOSS(NN)%LOSS(ND,ND2) = .TRUE.
         ENDDO
      ENDDO
      IF (DN%FILTER_INDEX > 0) N_LOSS(NN)%LOSS = .TRUE.
   ENDIF
ENDDO

! Set initial LOSS_D array indicating if a duct has a loss associate with forward (node 1 to 2) or reverse flow (node 2 to 1)
DUCTLOOP: DO ND = 1, N_DUCTS
   DU => DUCT(ND)
   ! Fixed flow is equivalent to having losses defined as no pressure solution is required to get the flow.
   IF (DU%VOLUME_FLOW_INITIAL<1.E6_EB .OR. DU%MASS_FLOW_INITIAL<1.E6_EB) THEN
      LOSS_D(ND,:) = .TRUE.
      FIXED_D(ND) = .TRUE.
   ENDIF
   IF (DU%ROUGHNESS > 0._EB) LOSS_D(ND,:) = .TRUE.
   IF (DU%FAN_INDEX > 0) THEN
      IF (FAN(DU%FAN_INDEX)%FAN_TYPE==1) THEN
         FIXED_D(ND) = .TRUE.
         LOSS_D(ND,:) = .TRUE.
      ELSE
         IF (DUCTNODE(DU%NODE_INDEX(1))%N_DUCTS<=2 .AND. DUCTNODE(DU%NODE_INDEX(2))%N_DUCTS<=2) LOSS_D(ND,:) = .TRUE.
      ENDIF
   ENDIF
   IF (DU%LOSS(1)>0._EB) LOSS_D(ND,1) = .TRUE.
   IF (DU%LOSS(2)>0._EB) LOSS_D(ND,2) = .TRUE.
   ! Check the nodes on either side of the duct to see if the node has losses
   DN => DUCTNODE(DU%NODE_INDEX(1))
   IF (DN%N_DUCTS==1) THEN
      IF (N_LOSS(DU%NODE_INDEX(1))%LOSS(1,1)) LOSS_D(ND,1) = .TRUE.
      IF (N_LOSS(DU%NODE_INDEX(1))%LOSS(1,2)) LOSS_D(ND,2) = .TRUE.
   ELSE
      DO ND2 = 1, DN%N_DUCTS
         IF (DN%DUCT_INDEX(ND2) == ND) THEN
            IF (ALL(N_LOSS(DU%NODE_INDEX(1))%LOSS(:,ND2)) .AND. &
                ALL(N_LOSS(DU%NODE_INDEX(1))%LOSS(ND2,:))) THEN
               LOSS_D(ND,1) = .TRUE.
               LOSS_D(ND,2) = .TRUE.
               EXIT
            ENDIF
         ENDIF
      ENDDO
   ENDIF
   DN => DUCTNODE(DU%NODE_INDEX(2))
   IF (DN%N_DUCTS==1) THEN
      IF (N_LOSS(DU%NODE_INDEX(2))%LOSS(1,1)) LOSS_D(ND,2) = .TRUE.
      IF (N_LOSS(DU%NODE_INDEX(2))%LOSS(1,2)) LOSS_D(ND,1) = .TRUE.
   ELSE
      DO ND2 = 1, DN%N_DUCTS
         IF (DN%DUCT_INDEX(ND2) == ND) THEN
            IF (ALL(N_LOSS(DU%NODE_INDEX(2))%LOSS(:,ND2)) .AND. &
                ALL(N_LOSS(DU%NODE_INDEX(2))%LOSS(ND2,:))) THEN
               LOSS_D(ND,1) = .TRUE.
               LOSS_D(ND,2) = .TRUE.
               ENDIF
            EXIT
         ENDIF
      ENDDO
   ENDIF
ENDDO DUCTLOOP
CHANGE = .TRUE.
CHANGELOOP: DO WHILE (CHANGE)
   CHANGE = .FALSE.
   NODELOOP: DO NN=1,N_DUCTNODES
      IF (LOSS_N(NN)) CYCLE NODELOOP
      IF (ALL(N_LOSS(NN)%LOSS)) THEN
         LOSS_N(NN) =.TRUE.
         CHANGE = .TRUE.
         DO ND=1,DUCTNODE(NN)%N_DUCTS
            IF (.NOT. LOSS_D(DUCTNODE(NN)%DUCT_INDEX(ND),1)) THEN
               LOSS_D(DUCTNODE(NN)%DUCT_INDEX(ND),1)= .TRUE.
               CHANGE = .TRUE.
            ENDIF
            IF (.NOT. LOSS_D(DUCTNODE(NN)%DUCT_INDEX(ND),2)) THEN
               LOSS_D(DUCTNODE(NN)%DUCT_INDEX(ND),2)= .TRUE.
               CHANGE = .TRUE.
            ENDIF
         ENDDO
         CYCLE NODELOOP
      ENDIF
      DN => DUCTNODE(NN)
      IF (DN%N_DUCTS==1) THEN
         IF (DUCT(DN%DUCT_INDEX(1))%NODE_INDEX(1)==NN) THEN
            IF (LOSS_D(DN%DUCT_INDEX(1),1) .AND. .NOT. N_LOSS(NN)%LOSS(1,1)) THEN
               N_LOSS(NN)%LOSS(1,1) = .TRUE.
               CHANGE = .TRUE.
            ENDIF
            IF (LOSS_D(DN%DUCT_INDEX(1),2) .AND. .NOT. N_LOSS(NN)%LOSS(1,2)) THEN
               N_LOSS(NN)%LOSS(1,2) = .TRUE.
               CHANGE = .TRUE.
            ENDIF
         ELSE
            IF (LOSS_D(DN%DUCT_INDEX(1),1) .AND. .NOT. N_LOSS(NN)%LOSS(1,2)) THEN
               N_LOSS(NN)%LOSS(1,2) = .TRUE.
               CHANGE = .TRUE.
            ENDIF
            IF (LOSS_D(DN%DUCT_INDEX(1),2) .AND. .NOT. N_LOSS(NN)%LOSS(1,1)) THEN
               N_LOSS(NN)%LOSS(1,1) = .TRUE.
               CHANGE = .TRUE.
            ENDIF
         ENDIF
         IF (ALL(N_LOSS(NN)%LOSS)) THEN
            LOSS_N(NN) = .TRUE.
            CHANGE = .TRUE.
         ENDIF
         CYCLE NODELOOP
      ENDIF
      IF (DN%N_DUCTS==2) THEN
         IF (N_LOSS(NN)%LOSS(1,2)) THEN
            IF (DUCT(DN%DUCT_INDEX(1))%NODE_INDEX(1)==NN) THEN
               IF (.NOT. LOSS_D(DN%DUCT_INDEX(1),2)) THEN
                  LOSS_D(DN%DUCT_INDEX(1),2) = .TRUE.
                  CHANGE = .TRUE.
               ENDIF
            ELSE
               IF (.NOT. LOSS_D(DN%DUCT_INDEX(1),1)) THEN
                  LOSS_D(DN%DUCT_INDEX(1),1) = .TRUE.
                  CHANGE = .TRUE.
               ENDIF
            ENDIF
            IF (DUCT(DN%DUCT_INDEX(2))%NODE_INDEX(1)==NN) THEN
               IF (.NOT. LOSS_D(DN%DUCT_INDEX(2),1)) THEN
                  LOSS_D(DN%DUCT_INDEX(2),1) = .TRUE.
                  CHANGE = .TRUE.
               ENDIF
            ELSE
               IF (.NOT. LOSS_D(DN%DUCT_INDEX(2),2)) THEN
                  LOSS_D(DN%DUCT_INDEX(2),2) = .TRUE.
                  CHANGE = .TRUE.
               ENDIF
            ENDIF
         ENDIF
         IF (N_LOSS(NN)%LOSS(2,1)) THEN
            IF (DUCT(DN%DUCT_INDEX(1))%NODE_INDEX(1)==NN) THEN
               IF (.NOT. LOSS_D(DN%DUCT_INDEX(1),1)) THEN
                  LOSS_D(DN%DUCT_INDEX(1),1) = .TRUE.
                  CHANGE = .TRUE.
               ENDIF
            ELSE
               IF (.NOT. LOSS_D(DN%DUCT_INDEX(1),2)) THEN
                  LOSS_D(DN%DUCT_INDEX(1),2) = .TRUE.
                  CHANGE = .TRUE.
               ENDIF
            ENDIF
            IF (DUCT(DN%DUCT_INDEX(2))%NODE_INDEX(1)==NN) THEN
               IF (.NOT. LOSS_D(DN%DUCT_INDEX(2),2)) THEN
                  LOSS_D(DN%DUCT_INDEX(2),2) = .TRUE.
                  CHANGE = .TRUE.
               ENDIF
            ELSE
               IF (.NOT. LOSS_D(DN%DUCT_INDEX(2),1)) THEN
                  LOSS_D(DN%DUCT_INDEX(2),1) = .TRUE.
                  CHANGE = .TRUE.
               ENDIF
            ENDIF
         ENDIF
         IF (ALL(LOSS_D(DN%DUCT_INDEX(1),:)) .OR. ALL(LOSS_D(DN%DUCT_INDEX(2),:))) THEN
            N_LOSS(NN)%LOSS = .TRUE.
            CHANGE = .TRUE.
         ELSE
            IF (DUCT(DN%DUCT_INDEX(1))%NODE_INDEX(1)==NN) THEN
               IF (LOSS_D(DN%DUCT_INDEX(1),1) .AND. .NOT. N_LOSS(NN)%LOSS(2,1)) THEN
                  N_LOSS(NN)%LOSS(2,1) = .TRUE.
                  CHANGE = .TRUE.
               ENDIF
               IF (LOSS_D(DN%DUCT_INDEX(1),2) .AND. .NOT. N_LOSS(NN)%LOSS(1,2)) THEN
                  N_LOSS(NN)%LOSS(1,2) = .TRUE.
                  CHANGE = .TRUE.
               ENDIF
            ELSE
               IF (LOSS_D(DN%DUCT_INDEX(1),1) .AND. .NOT. N_LOSS(NN)%LOSS(1,2)) THEN
                  N_LOSS(NN)%LOSS(1,2) = .TRUE.
                  CHANGE = .TRUE.
               ENDIF
               IF (LOSS_D(DN%DUCT_INDEX(1),2) .AND. .NOT. N_LOSS(NN)%LOSS(2,1)) THEN
                  N_LOSS(NN)%LOSS(2,1) = .TRUE.
                  CHANGE = .TRUE.
               ENDIF
            ENDIF
            IF (DUCT(DN%DUCT_INDEX(2))%NODE_INDEX(1)==NN) THEN
               IF (LOSS_D(DN%DUCT_INDEX(2),1) .AND. .NOT. N_LOSS(NN)%LOSS(1,2)) THEN
                  N_LOSS(NN)%LOSS(1,2) = .TRUE.
                  CHANGE = .TRUE.
               ENDIF
               IF (LOSS_D(DN%DUCT_INDEX(2),2) .AND. .NOT. N_LOSS(NN)%LOSS(2,1)) THEN
                  N_LOSS(NN)%LOSS(2,1) = .TRUE.
                  CHANGE = .TRUE.
               ENDIF
            ELSE
               IF (LOSS_D(DN%DUCT_INDEX(2),1) .AND. .NOT. N_LOSS(NN)%LOSS(2,1)) THEN
                  N_LOSS(NN)%LOSS(2,1) = .TRUE.
                  CHANGE = .TRUE.
               ENDIF
               IF (LOSS_D(DN%DUCT_INDEX(2),2) .AND. .NOT. N_LOSS(NN)%LOSS(1,2)) THEN
                  N_LOSS(NN)%LOSS(1,2) = .TRUE.
                  CHANGE = .TRUE.
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      F_COUNTER = 0
      DUCT_O_LOOP: DO ND=1,DN%N_DUCTS
         IF (FIXED_D(DN%DUCT_INDEX(ND))) F_COUNTER = F_COUNTER + 1
         IF (ALL(N_LOSS(NN)%LOSS(ND,:)) .AND. ALL(N_LOSS(NN)%LOSS(:,ND))) THEN
            IF (.NOT. LOSS_D(DN%DUCT_INDEX(ND),1)) THEN
               LOSS_D(DN%DUCT_INDEX(ND),1) = .TRUE.
               CHANGE = .TRUE.
            ENDIF
            IF (.NOT. LOSS_D(DN%DUCT_INDEX(ND),2)) THEN
               LOSS_D(DN%DUCT_INDEX(ND),2) = .TRUE.
               CHANGE = .TRUE.
            ENDIF
         ENDIF

         DUCT_I_LOOP: DO ND2=1,DN%N_DUCTS
            IF (DUCT(DN%DUCT_INDEX(ND2))%NODE_INDEX(1)==NN) THEN
               IF (LOSS_D(DN%DUCT_INDEX(ND2),1) .AND. .NOT. N_LOSS(NN)%LOSS(ND,ND2)) THEN
                  N_LOSS(NN)%LOSS(ND,ND2) = .TRUE.
                  CHANGE = .TRUE.
               ENDIF
               IF (.NOT. N_LOSS(NN)%LOSS(ND,ND2)) THEN
                  DN2=>DUCTNODE(DUCT(DN%DUCT_INDEX(ND2))%NODE_INDEX(2))
                  DO ND3=1,DN2%N_DUCTS
                     IF (DN2%DUCT_INDEX(ND3)==DN%DUCT_INDEX(ND2)) THEN
                        IF (ALL(N_LOSS(DUCT(DN2%DUCT_INDEX(ND3))%NODE_INDEX(2))%LOSS(ND3,:))) THEN
                           LOSS_D(DN2%DUCT_INDEX(ND3),1) = .TRUE.
                           CHANGE = .TRUE.
                        ENDIF
                        EXIT
                     ENDIF
                  ENDDO
               ENDIF
            ELSE
               IF (LOSS_D(DN%DUCT_INDEX(ND2),2) .AND. .NOT. N_LOSS(NN)%LOSS(ND,ND2)) THEN
                  N_LOSS(NN)%LOSS(ND,ND2) = .TRUE.
                  CHANGE = .TRUE.
               ENDIF
               IF (.NOT. N_LOSS(NN)%LOSS(ND,ND2)) THEN
                  DN2=>DUCTNODE(DUCT(DN%DUCT_INDEX(ND2))%NODE_INDEX(1))
                  DO ND3=1,DN2%N_DUCTS
                     IF (DN2%DUCT_INDEX(ND3)==DN%DUCT_INDEX(ND2)) THEN
                        IF (ALL(N_LOSS(DUCT(DN2%DUCT_INDEX(ND3))%NODE_INDEX(1))%LOSS(ND3,:))) THEN
                           LOSS_D(DN2%DUCT_INDEX(ND3),2) = .TRUE.
                           CHANGE = .TRUE.
                        ENDIF
                        EXIT
                     ENDIF
                  ENDDO
               ENDIF
            ENDIF
         ENDDO DUCT_I_LOOP
      ENDDO DUCT_O_LOOP
      IF (F_COUNTER >= DN%N_DUCTS - 1) THEN
         N_LOSS(NN)%LOSS = .TRUE.
         IF (.NOT. LOSS_N(NN)) THEN
            LOSS_N(NN) = .TRUE.
            CHANGE = .TRUE.
         ENDIF
         DO ND = 1,DN%N_DUCTS
            IF (.NOT. ALL(LOSS_D(DN%DUCT_INDEX(ND),:))) THEN
               LOSS_D(DN%DUCT_INDEX(ND),:) = .TRUE.
               CHANGE = .TRUE.
            ENDIF
         ENDDO
      ENDIF
   ENDDO NODELOOP
ENDDO CHANGELOOP

DEALLOCATE(N_LOSS)

IF (.NOT. ALL(LOSS_D) .OR. .NOT. ALL(LOSS_N)) THEN
   WRITE(MESSAGE,'(A,I5,A,I5)') 'ERROR(556): Problem with HVAC network. Insufficient LOSS definitions for DUCTs and NODEs'
   CALL SHUTDOWN(MESSAGE); RETURN
ENDIF

END SUBROUTINE EXAMINE_LOSSES

! The routines that follow are for solving system curves when HVAC_QFAN is enabled

!> \brief Obtains the system curves for HVAC fans
!>
!> \param T Current simulation time (s)

SUBROUTINE HVAC_QFAN_CALC(T)
REAL(EB), INTENT(IN) :: T
INTEGER :: NR,NF
TYPE(DUCTRUN_TYPE), POINTER :: DR

! Obtain system curves for HVAC fans
DR_LOOP: DO NR = 1, N_DUCTRUNS
   DR => DUCTRUN(NR)
   IF (DR%N_QFANS ==0) CYCLE
   CALL SET_DONOR_QFAN(NR,0)
   DO NF=1,DR%N_QFANS
      IF (DUCT(DR%FAN_INDEX(NF))%FAN_OPERATING) DR%FAN_OPERATING(NF) = .TRUE.
      CALL SET_DONOR_QFAN(NR,NF)
   ENDDO

   FAN_OP_IF: IF (ANY(DR%FAN_OPERATING)) THEN
      IF (DR%N_DUCTS==1) THEN
         DO NF=0,DR%N_QFANS
            CALL UPDATE_LOSS_QFAN(T,NR,NF)
         ENDDO
         CALL QFAN_SIMPLE(NR,T)
         DO NF=0,DR%N_QFANS
            CALL SET_DONOR_QFAN(NR,NF)
         ENDDO
      ELSE
         ALLOCATE(LHS(DR%N_M_DUCTS+DR%N_M_DUCTNODES,DR%N_M_DUCTS+DR%N_M_DUCTNODES))
         ALLOCATE(RHS(DR%N_M_DUCTS+DR%N_M_DUCTNODES))
         DO NF=0,DR%N_QFANS
            IF (NF/=0) THEN
               IF (.NOT. DR%FAN_OPERATING(NF)) CYCLE
            ENDIF
            ITER = 0
            DO WHILE (ITER < ITER_MAX_QFAN)
               LHS = 0._EB
               RHS = 0._EB
               IF (N_AIRCOILS > 0) CALL COIL_UPDATE_QFAN(T,NR,NF)
               CALL UPDATE_LOSS_QFAN(T,NR,NF)
               CALL RHSNODE_QFAN(NR)
               CALL RHSDUCT_QFAN(NR,NF,T)
               CALL LHSNODE_QFAN(NR,NF)
               CALL LHSDUCT_QFAN(NR,NF)
               CALL MATRIX_SOLVE_QFAN(NR,NF)
               CALL CONVERGENCE_CHECK_QFAN(NR,NF)
               CALL SET_GUESS_QFAN(NR,NF)
               CALL HVAC_UPDATE_QFAN(NR,NF)
               CALL SET_DONOR_QFAN(NR,NF)
               ITER = ITER + 1
            ENDDO
            DUCTRUN(NR)%VEL(:,NF,PREVIOUS) = DUCTRUN(NR)%VEL(:,NF,NEW)
            DUCTRUN(NR)%VEL(:,NF,GUESS) = DUCTRUN(NR)%VEL(:,NF,NEW)
            DUCTRUN(NR)%VEL(:,NF,OLD) = DUCTRUN(NR)%VEL(:,NF,NEW)
            DUCTRUN(NR)%P(:,NF,OLD) = DUCTRUN(NR)%P(:,NF,NEW)
         ENDDO
         DEALLOCATE(LHS)
         DEALLOCATE(RHS)
      ENDIF
      ! Deallocate matrices used for solving steady state system curve
   ENDIF FAN_OP_IF
ENDDO DR_LOOP

END SUBROUTINE HVAC_QFAN_CALC

!> \brief Updates PREVIOUS and GUESS values for QFAN solutiokn
!>
!> \param DUCTRUN_INDEX Index indicating which HVAC ductrun is being solved
!> \param NF Index indicating which ductrun fan is being solved

SUBROUTINE SET_GUESS_QFAN(DUCTRUN_INDEX,NF)
INTEGER, INTENT(IN) :: DUCTRUN_INDEX,NF
INTEGER :: ND
TYPE(DUCTRUN_TYPE), POINTER :: DR

DR => DUCTRUN(DUCTRUN_INDEX)

DO ND = 1,DR%N_M_DUCTS
   IF (SIGN(1._EB,DR%VEL(ND,NF,PREVIOUS)) == SIGN(1._EB,DR%VEL(ND,NF,NEW))) THEN
      DR%VEL(ND,NF,GUESS) = DR%VEL(ND,NF,NEW)
   ELSE
      DR%VEL(ND,NF,GUESS) = 0._EB
   ENDIF
   DR%VEL(ND,NF,PREVIOUS) = DR%VEL(ND,NF,NEW)
ENDDO

END SUBROUTINE SET_GUESS_QFAN

!> \brief Does the analytic solution when a ductrun consists of only one duct.
!>
!> \param DUCTRUN_INDEX Index for the ductrun
!> \param T Current time (s)

SUBROUTINE QFAN_SIMPLE(DUCTRUN_INDEX,T)
INTEGER, INTENT(IN) :: DUCTRUN_INDEX
REAL(EB), INTENT(IN) :: T
REAL(EB) :: P1,P2,FAN_P
TYPE(DUCTRUN_TYPE), POINTER :: DR

DR => DUCTRUN(DUCTRUN_INDEX)

P1 = DUCTNODE(DUCT(DR%DUCT_INDEX(1))%NODE_INDEX(1))%P
P2 = DUCTNODE(DUCT(DR%DUCT_INDEX(1))%NODE_INDEX(2))%P

DR%VEL(1,0,1) = SQRT(2._EB*ABS(P1 - P2)/(DR%LOSS(1,0)*DR%RHO_D(1,0)))
IF (P2 > P1) DR%VEL(1,0,1) = -DR%VEL(1,0,1)

IF (DR%FAN_OPERATING(1)) THEN
   FAN_P = FAN_MAX(DUCT(DR%DUCT_INDEX(1))%FAN_INDEX,T-DUCT(DR%DUCT_INDEX(1))%FAN_ON_TIME)
   IF (DUCT(DR%DUCT_INDEX(1))%REVERSE) FAN_P = -FAN_P
   DR%VEL(1,1,1) = SQRT(2._EB*ABS(P1 - P2 + FAN_P)/(DR%LOSS(1,1)*DR%RHO_D(1,1)))
   IF (P1 - P2 + FAN_P < 0._EB) DR%VEL(1,1,1) = -DR%VEL(1,1,1)
ELSE
   DR%VEL(1,1,1) = DR%VEL(1,0,1)
ENDIF

END SUBROUTINE QFAN_SIMPLE

!> \brief Finds "duct runs"; being ductnodes and ducts directly (via HVAC components) connected to one another

SUBROUTINE FIND_DUCTRUNS(CHANGE)
USE PHYSICAL_FUNCTIONS, ONLY : GET_ENTHALPY
LOGICAL, INTENT(IN) :: CHANGE
INTEGER :: NN,NR,NF,NN3,ND,DUCT_COUNTER(N_DUCTS),NODE_COUNTER(N_DUCTNODES),&
           DR_DUCTS(N_DUCTS),DR_DUCTNODES(N_DUCTS),DUCTRUN_MAP(N_DUCTS),DR_INDEX
LOGICAL :: NODE_CHECKED(N_DUCTNODES),NODE_CONNECTED(N_DUCTNODES),DUCT_FOUND,FAN_PRESENT(N_DUCTS),MT_PRESENT(N_DUCTS)
REAL(EB) :: C0,ZZ_GET(1:N_TRACKED_SPECIES)
TYPE(DUCT_TYPE), POINTER :: DU
TYPE(DUCTNODE_TYPE), POINTER :: DN

IF (.NOT. CHANGE .AND. ALLOCATED(DUCTRUN)) RETURN
IF (ALLOCATED(DUCTRUN)) DEALLOCATE(DUCTRUN)

ZZ_GET = SPECIES_MIXTURE(1:N_TRACKED_SPECIES)%ZZ0
CALL GET_ENTHALPY(ZZ_GET,C0,TMPA)
C0 = C0/TMPA

DUCT_COUNTER=0
NODE_COUNTER=0
FAN_PRESENT=.FALSE.
NODE_CONNECTED=.FALSE.
MT_PRESENT = .FALSE.
NODE_CHECKED=.FALSE.

DUCT%QFAN_INDEX = -1

! Set leakage paths as already checked, they cannot have a fan.
DO NN = 1, N_DUCTNODES
   DN=>DUCTNODE(NN)
   IF (DN%LEAKAGE .OR. DUCT(DN%DUCT_INDEX(1))%LOCALIZED_LEAKAGE) NODE_CHECKED(NN) = .TRUE.
END DO

N_DUCTRUNS = 1
NN = 1

! Finds all connected ducts and nodes and hence number of ductruns
L1:DO WHILE (NN <= N_DUCTNODES)
   IF (NODE_CHECKED(NN)) THEN
      NN = NN + 1
      CYCLE L1
   END IF
   DN=>DUCTNODE(NN)
   DUCT_FOUND = .FALSE.
   DL: DO ND = 1, DN%N_DUCTS
      DU=>DUCT(DN%DUCT_INDEX(ND))
      IF (DU%AREA < TWO_EPSILON_EB) CYCLE DL
      DUCT_FOUND = .TRUE.
      DUCT_COUNTER(DN%DUCT_INDEX(ND)) = N_DUCTRUNS
      IF (DU%N_CELLS > 0) MT_PRESENT(N_DUCTRUNS)=.TRUE.
      IF (DU%NODE_INDEX(1)==NN) THEN
         NODE_COUNTER(DU%NODE_INDEX(2)) = N_DUCTRUNS
         NODE_CONNECTED(DU%NODE_INDEX(2)) = .TRUE.
         IF (DUCTNODE(DU%NODE_INDEX(2))%N_DUCTS==1) NODE_CHECKED(DU%NODE_INDEX(2)) = .TRUE.
      ELSE
         NODE_COUNTER(DU%NODE_INDEX(1)) = N_DUCTRUNS
         NODE_CONNECTED(DU%NODE_INDEX(1)) = .TRUE.
         IF (DUCTNODE(DU%NODE_INDEX(1))%N_DUCTS==1) NODE_CHECKED(DU%NODE_INDEX(1)) = .TRUE.
      ENDIF
      IF (HVAC_QFAN .AND. DU%FAN_INDEX > 0) THEN
         IF (FAN(DU%FAN_INDEX)%FAN_TYPE>1) FAN_PRESENT(N_DUCTRUNS)=.TRUE.
      ENDIF
   END DO DL

   NODE_CHECKED(NN) = .TRUE.

   IF (DUCT_FOUND) THEN
      NODE_COUNTER(NN) = N_DUCTRUNS
      NODE_CONNECTED(NN) = .TRUE.
   ELSE
      ! If no open ducts attached to the node, it doesn't belong to a duct run.
      NN = NN + 1
      CYCLE L1
   ENDIF

   IF (ALL(NODE_CHECKED)) EXIT L1
   L2:DO NN3 = 1, N_DUCTNODES
      IF (NODE_CHECKED(NN3)) THEN
         CYCLE L2
      ELSE IF ((.NOT. NODE_CHECKED(NN3)) .AND. NODE_CONNECTED(NN3)) THEN
         ! If node attached to current node hasn't been checked yet, set the active node to that node
         NN = NN3
         EXIT L2
      ELSE IF (NN3 == N_DUCTNODES) THEN
         ! If there are no nodes part of the current network that haven't been checked, up the network and move to first unchecked
         N_DUCTRUNS = N_DUCTRUNS + 1
         NODE_CONNECTED = .FALSE.
         NN = FINDLOC(NODE_CHECKED,.FALSE.,1)
         EXIT L2
      END IF
   END DO L2
END DO L1

DR_DUCTS = 0
DR_DUCTNODES = 0

! Sums up number of ducts and ductnodes for DUCTRUN arrays.
DO ND = 1, N_DUCTS
   IF (DUCT_COUNTER(ND)==0) CYCLE
   DR_INDEX =DUCT_COUNTER(ND)
   IF (DR_INDEX>0) DR_DUCTS(DR_INDEX) = DR_DUCTS(DR_INDEX) + 1
ENDDO

DO NN = 1, N_DUCTNODES
   IF (NODE_COUNTER(NN)==0) CYCLE
   DR_INDEX =NODE_COUNTER(NN)
   IF (DR_INDEX>0) DR_DUCTNODES(DR_INDEX) = DR_DUCTNODES(DR_INDEX) + 1
ENDDO

DUCTRUN_MAP = 0

ALLOCATE(DUCTRUN(N_DUCTRUNS))
DUCTRUN(1:N_DUCTRUNS)%N_DUCTS = DR_DUCTS(1:N_DUCTRUNS)
DUCTRUN(1:N_DUCTRUNS)%N_DUCTNODES = DR_DUCTNODES(1:N_DUCTRUNS)
DUCTRUN(1:N_DUCTRUNS)%N_QFANS = 0

! Setup indexing arrays

DO NR = 1, N_DUCTRUNS
   ALLOCATE(DUCTRUN(NR)%DUCT_INDEX(DUCTRUN(NR)%N_DUCTS))
   ALLOCATE(DUCTRUN(NR)%NODE_INDEX(DUCTRUN(NR)%N_DUCTNODES))
ENDDO

DUCTRUN_MAP=0
DO ND=1,N_DUCTS
   DR_INDEX = DUCT_COUNTER(ND)
   IF (DR_INDEX > 0) THEN
      DUCTRUN_MAP(DR_INDEX) = DUCTRUN_MAP(DR_INDEX) + 1
      DUCT(ND)%DUCTRUN = DR_INDEX
      DUCT(ND)%DUCTRUN_INDEX = DUCTRUN_MAP(DR_INDEX)
      DUCTRUN(DR_INDEX)%DUCT_INDEX(DUCTRUN_MAP(DR_INDEX)) = ND
      IF (.NOT. DUCT(ND)%FIXED .AND. DUCT(ND)%AREA > TWO_EPSILON_EB) &
         DUCTRUN(DR_INDEX)%N_M_DUCTS = DUCTRUN(DR_INDEX)%N_M_DUCTS + 1
      IF (HVAC_QFAN .AND. DUCT(ND)%FAN_INDEX > 0) THEN
         IF (FAN(DUCT(ND)%FAN_INDEX)%FAN_TYPE > 1) THEN
            DUCTRUN(DR_INDEX)%N_QFANS = DUCTRUN(DR_INDEX)%N_QFANS + 1
            DUCT(ND)%QFAN_INDEX = DUCTRUN(DR_INDEX)%N_QFANS
         ENDIF
      ENDIF
   ENDIF
ENDDO

DUCTRUN_MAP=0
DO NN=1,N_DUCTNODES
   DR_INDEX = NODE_COUNTER(NN)
   IF (DR_INDEX > 0) THEN
      DUCTRUN_MAP(DR_INDEX) = DUCTRUN_MAP(DR_INDEX) + 1
      DUCTNODE(NN)%DUCTRUN = DR_INDEX
      DUCTNODE(NN)%DUCTRUN_INDEX = DUCTRUN_MAP(DR_INDEX)
      DUCTRUN(DR_INDEX)%NODE_INDEX(DUCTRUN_MAP(DR_INDEX)) = NN
      IF (DUCTNODE(NN)%AMBIENT .OR. DUCTNODE(NN)%LEAKAGE .OR. DUCTNODE(NN)%VENT) CYCLE
      DUCTRUN(DR_INDEX)%N_M_DUCTNODES = DUCTRUN(DR_INDEX)%N_M_DUCTNODES + 1
   ENDIF
ENDDO

! Allocates and populates DUCTRUN solution arrays when a QFAN is present
! VEL and P arrays are for solution matrix and use _M_ dimension
DO NR = 1, N_DUCTRUNS
   QFAN_IF: IF (DUCTRUN(NR)%N_QFANS > 0) THEN
      IF (DUCTRUN(NR)%N_M_DUCTS > 0) THEN
         ALLOCATE(DUCTRUN(NR)%DUCT_M_INDEX(DUCTRUN(NR)%N_M_DUCTS))
         DUCTRUN(NR)%DUCT_M_INDEX = -1
         NN3 = 0
         DO ND=1,DUCTRUN(NR)%N_DUCTS
            IF (DUCT(DUCTRUN(NR)%DUCT_INDEX(ND))%FIXED .OR. DUCT(DUCTRUN(NR)%DUCT_INDEX(ND))%AREA < TWO_EPSILON_EB) CYCLE
            NN3 = NN3 + 1
            DUCTRUN(NR)%DUCT_M_INDEX(NN3) = DUCTRUN(NR)%DUCT_INDEX(ND)
            DUCT(DUCTRUN(NR)%DUCT_INDEX(ND))%DUCTRUN_M_INDEX=NN3
         ENDDO
      ENDIF
      IF (DUCTRUN(NR)%N_M_DUCTNODES > 0) THEN
         ALLOCATE(DUCTRUN(NR)%NODE_M_INDEX(DUCTRUN(NR)%N_M_DUCTNODES))
         DUCTRUN(NR)%NODE_M_INDEX = -1
         NN3 = 0
         DO NN=1,DUCTRUN(NR)%N_DUCTNODES
            IF (DUCTNODE(DUCTRUN(NR)%NODE_INDEX(NN))%VENT .OR. &
                DUCTNODE(DUCTRUN(NR)%NODE_INDEX(NN))%AMBIENT .OR. &
                DUCTNODE(DUCTRUN(NR)%NODE_INDEX(NN))%LEAKAGE) CYCLE
            NN3 = NN3 + 1
            DUCTRUN(NR)%NODE_M_INDEX(NN3) = DUCTRUN(NR)%NODE_INDEX(NN)
            DUCTNODE(DUCTRUN(NR)%NODE_INDEX(NN))%DUCTRUN_M_INDEX = NN3
         ENDDO
      ENDIF

      ALLOCATE(DUCTRUN(NR)%FAN_OPERATING(DUCTRUN(NR)%N_QFANS))
      DUCTRUN(NR)%FAN_OPERATING = .FALSE.
      ALLOCATE(DUCTRUN(NR)%FAN_INDEX(DUCTRUN(NR)%N_QFANS))
      DO ND=1,DUCTRUN(NR)%N_DUCTS
         IF (DUCT(DUCTRUN(NR)%DUCT_INDEX(ND))%QFAN_INDEX > 0) &
            DUCTRUN(NR)%FAN_INDEX(DUCT(DUCTRUN(NR)%DUCT_INDEX(ND))%QFAN_INDEX) = DUCTRUN(NR)%DUCT_INDEX(ND)
      ENDDO

      ALLOCATE(DUCTRUN(NR)%RHO_D(DUCTRUN(NR)%N_DUCTS,0:DUCTRUN(NR)%N_QFANS))
      ALLOCATE(DUCTRUN(NR)%TMP_D(DUCTRUN(NR)%N_DUCTS,0:DUCTRUN(NR)%N_QFANS))
      ALLOCATE(DUCTRUN(NR)%CP_D(DUCTRUN(NR)%N_DUCTS,0:DUCTRUN(NR)%N_QFANS))
      ALLOCATE(DUCTRUN(NR)%ZZ_D(DUCTRUN(NR)%N_DUCTS,0:DUCTRUN(NR)%N_QFANS,N_TRACKED_SPECIES))
      ALLOCATE(DUCTRUN(NR)%LOSS(DUCTRUN(NR)%N_DUCTS,0:DUCTRUN(NR)%N_QFANS))
      ALLOCATE(DUCTRUN(NR)%VEL(DUCTRUN(NR)%N_M_DUCTS,0:DUCTRUN(NR)%N_QFANS,4)) ! VEL2(duct,fan,old/new/guess/previous)

      DUCTRUN(NR)%RHO_D = RHOA
      DUCTRUN(NR)%TMP_D = TMPA
      DUCTRUN(NR)%CP_D = C0
      DO ND = 1,DUCTRUN(NR)%N_DUCTS
          DO NF=1,DUCTRUN(NR)%N_QFANS
             DUCTRUN(NR)%ZZ_D(ND,NF,:) = ZZ_GET
          ENDDO
      ENDDO
      DO ND = 1,DUCTRUN(NR)%N_M_DUCTS
          DUCTRUN(NR)%VEL(ND,:,:) = DUCT(DUCTRUN(NR)%DUCT_M_INDEX(ND))%VEL(OLD)
      ENDDO

      ALLOCATE(DUCTRUN(NR)%P(DUCTRUN(NR)%N_M_DUCTNODES,0:DUCTRUN(NR)%N_QFANS,2)) ! P2(node,fan,old/new)
      ALLOCATE(DUCTRUN(NR)%RHO_N(DUCTRUN(NR)%N_DUCTNODES,0:DUCTRUN(NR)%N_QFANS))
      ALLOCATE(DUCTRUN(NR)%TMP_N(DUCTRUN(NR)%N_DUCTNODES,0:DUCTRUN(NR)%N_QFANS))
      ALLOCATE(DUCTRUN(NR)%CP_N(DUCTRUN(NR)%N_DUCTNODES,0:DUCTRUN(NR)%N_QFANS))
      ALLOCATE(DUCTRUN(NR)%ZZ_N(DUCTRUN(NR)%N_DUCTNODES,0:DUCTRUN(NR)%N_QFANS,N_TRACKED_SPECIES))

      DO NN = 1,DUCTRUN(NR)%N_DUCTNODES
         DUCTRUN(NR)%RHO_N(NN,:) = DUCTNODE(DUCTRUN(NR)%NODE_INDEX(NN))%RHO_OLD
         DUCTRUN(NR)%TMP_N(NN,:) = DUCTNODE(DUCTRUN(NR)%NODE_INDEX(NN))%TMP_OLD
         DUCTRUN(NR)%CP_N(NN,:) = DUCTNODE(DUCTRUN(NR)%NODE_INDEX(NN))%CP_OLD
         DO NF = 0,DUCTRUN(NR)%N_QFANS
            DUCTRUN(NR)%ZZ_N(NN,NF,:) = DUCTNODE(DUCTRUN(NR)%NODE_INDEX(NN))%ZZ_OLD(:)
         ENDDO
      ENDDO
      DO NN = 1,DUCTRUN(NR)%N_M_DUCTNODES
         DUCTRUN(NR)%P(NN,:,:) = DUCTNODE(DUCTRUN(NR)%NODE_M_INDEX(NN))%P - P_INF
      ENDDO

      DO NN = 1,DUCTRUN(NR)%N_M_DUCTNODES
      ENDDO
   ELSE
      DUCTRUN(NR)%N_M_DUCTNODES = 0
      DUCTRUN(NR)%N_M_DUCTS = 0
   ENDIF QFAN_IF
ENDDO

DO NR=1,N_DUCTRUNS
   IF (.NOT. MT_PRESENT(NR)) CYCLE
   DO ND=1,DUCTRUN(NR)%N_DUCTS
      IF (DUCT(DUCTRUN(NR)%DUCT_INDEX(ND))%N_CELLS==0) THEN
         WRITE(MESSAGE,'(A,A,A)') &
            'ERROR(556): Duct ',TRIM(DUCT(DUCTRUN(NR)%DUCT_INDEX(ND))%ID),&
            ' has N_CELLS=0 and is in a duct run with ducts that have N_CELLS>0.'
         CALL SHUTDOWN(MESSAGE); RETURN
      ENDIF
   ENDDO
ENDDO

END SUBROUTINE FIND_DUCTRUNS

!> \brief Solves the HVAC matrix for a ductrun used when multiple fans exist within a single run of duct
!>
!> \param DUCTRUN_INDEX Index indicating which HVAC ductrun is being solved
!> \param NF Index indicating which HVAC ductrun fan is being solved

SUBROUTINE MATRIX_SOLVE_QFAN(DUCTRUN_INDEX,NF)
USE MATH_FUNCTIONS,ONLY : GAUSSJ
INTEGER,INTENT(IN) :: DUCTRUN_INDEX,NF
INTEGER :: IERR
TYPE(DUCTRUN_TYPE), POINTER :: DR

DR =>DUCTRUN(DUCTRUN_INDEX)

CALL GAUSSJ(LHS,DR%N_M_DUCTS+DR%N_M_DUCTNODES,DR%N_M_DUCTS+DR%N_M_DUCTNODES,RHS,1,1,IERR)

DR%VEL(1:DR%N_M_DUCTS,NF,NEW) =RHS(1:DR%N_M_DUCTS)
DR%P(1:DR%N_M_DUCTNODES,NF,NEW) = RHS(DR%N_M_DUCTS+1:DR%N_M_DUCTS+DR%N_M_DUCTNODES)

END SUBROUTINE MATRIX_SOLVE_QFAN


!> \brief Builds the right hand side of the HVAC QFAN matrix for mass conservation at internal nodes
!> \param DUCTRUN_INDEX Index indicating which HVAC ductrun is being solved

SUBROUTINE RHSNODE_QFAN(DUCTRUN_INDEX)
USE GLOBAL_CONSTANTS
INTEGER, INTENT(IN)::DUCTRUN_INDEX
INTEGER :: NN,ND
TYPE(DUCTRUN_TYPE), POINTER::DR
TYPE(DUCT_TYPE), POINTER::DU
TYPE(DUCTNODE_TYPE), POINTER::DN

DR => DUCTRUN(DUCTRUN_INDEX)
DO NN = 1, DR%N_M_DUCTNODES
   DN => DUCTNODE(DR%NODE_INDEX(NN))
   DO ND = 1,DN%N_DUCTS
      DU => DUCT(DN%DUCT_INDEX(ND))
      IF (DU%FIXED) RHS(DR%N_M_DUCTS+NN) = RHS(DR%N_M_DUCTS+NN) + DN%DIR(ND)*DU%RHO_D*DU%VOLUME_FLOW
   END DO
ENDDO

END SUBROUTINE RHSNODE_QFAN

!> \brief Builds the left hand side of the HVAC QFAN matrix for mass conservation at internal nodes
!>
!> \param DUCTRUN_INDEX Index indicating which HVAC ductrun is being solved
!> \param NF Index indicating which ductrun fan is being solved

SUBROUTINE LHSNODE_QFAN(DUCTRUN_INDEX,NF)
USE GLOBAL_CONSTANTS
INTEGER, INTENT(IN)::DUCTRUN_INDEX,NF
INTEGER :: NN,ND,ARRAYLOC
TYPE(DUCTRUN_TYPE), POINTER::DR
TYPE(DUCT_TYPE), POINTER::DU
TYPE(DUCTNODE_TYPE), POINTER::DN

DR => DUCTRUN(DUCTRUN_INDEX)
DO NN = 1, DR%N_M_DUCTNODES
   DN => DUCTNODE(DR%NODE_M_INDEX(NN))
   DO ND = 1,DN%N_DUCTS
      DU => DUCT(DN%DUCT_INDEX(ND))
      ARRAYLOC = DU%DUCTRUN_M_INDEX
      IF (ARRAYLOC <=0) CYCLE
      LHS(DR%N_M_DUCTS+NN,ARRAYLOC) = -DN%DIR(ND)*DR%RHO_D(DU%DUCTRUN_INDEX,NF)*DU%AREA
   END DO
ENDDO

END SUBROUTINE LHSNODE_QFAN


!> \brief Builds the right hand side of the HVAC QFAN flow matrix for momentum conservation in a duct
!> \param DUCTRUN_INDEX Index indicating which HVAC ductrun is being solved
!> \param NF Index indicating which ductrun fan is being solved
!> \param T Current time (s)

SUBROUTINE RHSDUCT_QFAN(DUCTRUN_INDEX,NF,T)

USE GLOBAL_CONSTANTS
INTEGER, INTENT(IN) :: DUCTRUN_INDEX,NF
REAL(EB), INTENT(IN) :: T
INTEGER :: ND
REAL(EB) :: HEAD !< DUCT pressure head (Pa)
REAL(EB) :: RGZ !< DUCT elevation head (Pa)
REAL(EB) :: XYZ(3) !< Distance between duct nodes (m)
REAL(EB) :: DP_FAN !< Max pressure of fan (Pa)
TYPE(DUCTRUN_TYPE), POINTER::DR
TYPE(DUCT_TYPE), POINTER::DU
TYPE(DUCTNODE_TYPE), POINTER::DN

DR => DUCTRUN(DUCTRUN_INDEX)

DO ND = 1, DR%N_M_DUCTS
   DU => DUCT(DR%DUCT_M_INDEX(ND))
   HEAD = 0._EB
   RGZ = 0._EB
   DN=>DUCTNODE(DU%NODE_INDEX(1))
   IF (DN%DUCTRUN_M_INDEX <=0) HEAD = DN%P
   IF ((DN%VENT .OR. DN%GEOM) .AND. HVAC_LOCAL_PRESSURE) &
      HEAD = HEAD - DN%DIR(1) * DR%RHO_D(ND,NF) * DU%VEL(OLD) * DU%AREA / NODE_AREA_EX(DU%NODE_INDEX(1)) / (GAMMA * DT_QF)
   XYZ = DN%XYZ
   DN=>DUCTNODE(DU%NODE_INDEX(2))
   IF (DN%DUCTRUN_M_INDEX <=0) HEAD = HEAD - DN%P
   IF ((DN%VENT .OR. DN%GEOM) .AND. HVAC_LOCAL_PRESSURE) &
      HEAD = HEAD + DN%DIR(1) * DR%RHO_D(ND,NF) * DU%VEL(OLD) * DU%AREA / NODE_AREA_EX(DU%NODE_INDEX(2)) / (GAMMA * DT_QF)
   XYZ = DN%XYZ - XYZ
   IF (.NOT. DU%LEAKAGE) THEN
      IF (STRATIFICATION) THEN
         RGZ = (GVEC(1)*XYZ(1)+GVEC(2)*XYZ(2)+GVEC(3)*XYZ(3))*DU%RHO_D
      ELSE
         RGZ = (GVEC(1)*XYZ(1)+GVEC(2)*XYZ(2)+GVEC(3)*XYZ(3))* &
               (DUCTNODE(DU%NODE_INDEX(1))%RHO - DUCTNODE(DU%NODE_INDEX(2))%RHO)
      ENDIF
   ENDIF
   IF (DU%QFAN_INDEX==NF) THEN
      IF ( DU%REVERSE) THEN
         DP_FAN = -FAN_MAX(DU%FAN_INDEX,T-DU%FAN_ON_TIME)
      ELSE
         DP_FAN = FAN_MAX(DU%FAN_INDEX,T-DU%FAN_ON_TIME)
      ENDIF
   ELSE
      DP_FAN = 0._EB
   ENDIF

   RHS(ND) = DR%VEL(ND,NF,OLD)+DT_QF/DU%LENGTH*((HEAD+RGZ+DP_FAN)/DR%RHO_D(ND,NF) + &
                   0.5_EB*DR%LOSS(ND,NF)*ABS(DR%VEL(ND,NF,PREVIOUS))*DR%VEL(ND,NF,GUESS))
ENDDO

END SUBROUTINE RHSDUCT_QFAN

!> \brief Builds the left hand side of the HVAC QFAN flow matrix for momentum conservation in a duct
!>
!> \param DUCTRUN_INDEX Index indicating which HVAC ductrun is being solved
!> \param NF Index indicating which ductrun fan is being solved

SUBROUTINE LHSDUCT_QFAN(DUCTRUN_INDEX,NF)
USE GLOBAL_CONSTANTS
INTEGER, INTENT(IN)::DUCTRUN_INDEX,NF
INTEGER :: ND
TYPE(DUCTRUN_TYPE), POINTER::DR
TYPE(DUCT_TYPE), POINTER::DU
TYPE(DUCTNODE_TYPE), POINTER::DN

DR => DUCTRUN(DUCTRUN_INDEX)

DUCT_LOOP: DO ND = 1, DR%N_M_DUCTS
   DU => DUCT(DR%DUCT_M_INDEX(ND))
   LHS(ND,ND) = 1._EB+0.5_EB*DT_QF*DR%LOSS(ND,NF)/DU%LENGTH*ABS(DR%VEL(ND,NF,PREVIOUS)+DR%VEL(ND,NF,GUESS))
   DN=>DUCTNODE(DU%NODE_INDEX(1))
   IF (DN%DUCTRUN_M_INDEX > 0) LHS(ND,DR%N_M_DUCTS + DN%DUCTRUN_M_INDEX) = -DT_QF/(DR%RHO_D(ND,NF)*DU%LENGTH)
   IF ((DN%VENT .OR. DN%GEOM) .AND. HVAC_LOCAL_PRESSURE) &
      LHS(ND,ND) = LHS(ND,ND)-DN%DIR(1)*DU%AREA/(DU%LENGTH*GAMMA*NODE_AREA_EX(DU%NODE_INDEX(1)))
   DN=>DUCTNODE(DU%NODE_INDEX(2))
   IF (DN%DUCTRUN_M_INDEX > 0) LHS(ND,DR%N_M_DUCTS + DN%DUCTRUN_M_INDEX) =  DT_QF/(DR%RHO_D(ND,NF)*DU%LENGTH)
   IF ((DN%VENT .OR. DN%GEOM) .AND. HVAC_LOCAL_PRESSURE) &
      LHS(ND,ND) = LHS(ND,ND)+DN%DIR(1)*DU%AREA/(DU%LENGTH*GAMMA*NODE_AREA_EX(DU%NODE_INDEX(2)))
ENDDO DUCT_LOOP

END SUBROUTINE LHSDUCT_QFAN

!> \brief Get max fan pressure for a fan
!>
!> \param FAN_INDEX Index of fan type
!> \param TSI Length of time fan has been operating (s)

REAL(EB) FUNCTION FAN_MAX(FAN_INDEX,TSI)
USE MATH_FUNCTIONS, ONLY : EVALUATE_RAMP
INTEGER, INTENT(IN) :: FAN_INDEX
REAL(EB),INTENT(IN) :: TSI

SELECT CASE (FAN(FAN_INDEX)%FAN_TYPE)
   CASE (2)
      FAN_MAX = FAN(FAN_INDEX)%MAX_PRES
   CASE (3)
      FAN_MAX = MAXVAL(RAMPS(FAN(FAN_INDEX)%RAMP_INDEX)%DEPENDENT_DATA)
END SELECT

FAN_MAX = FAN_MAX * EVALUATE_RAMP(TSI,FAN(FAN_INDEX)%SPIN_INDEX,TAU=FAN(FAN_INDEX)%TAU)

END FUNCTION FAN_MAX

!> \brief Checks the current iteration for duct velocity convergence and conservation of mass at nodes
!>
!> \param DUCTRUN_INDEX Index indicating which HVAC ductrun is being solved
!> \param NF Index indicating which HVAC ductrun fan is being solved

SUBROUTINE CONVERGENCE_CHECK_QFAN(DUCTRUN_INDEX,NF)
INTEGER, INTENT(IN) :: DUCTRUN_INDEX,NF
INTEGER :: NN, ND, COUNT
REAL(EB) :: MSUM,MTOT,MFLOW,VEL
LOGICAL :: CONVERGED
TYPE(DUCTRUN_TYPE), POINTER :: DR
TYPE(DUCT_TYPE), POINTER :: DU
TYPE(DUCTNODE_TYPE), POINTER :: DN

DR => DUCTRUN(DUCTRUN_INDEX)
CONVERGED = .TRUE.

! Check duct velocity convergence
DO ND=1,DR%N_M_DUCTS
   IF (ABS(DR%VEL(ND,NF,PREVIOUS)) < 1.E-5_EB .AND. ABS(DR%VEL(ND,NF,NEW)) < 1.E-5_EB) CYCLE
   IF (ABS(DR%VEL(ND,NF,NEW)-DR%VEL(ND,NF,PREVIOUS))/(ABS(DR%VEL(ND,NF,PREVIOUS))+TWO_EPSILON_EB) > 0.01_EB) CONVERGED = .FALSE.
ENDDO

IF (.NOT. CONVERGED) RETURN ! if velocity convergence passes, continue

! Check node mass conservation convergence
DO NN=1,DR%N_M_DUCTNODES
   DN => DUCTNODE(DR%NODE_INDEX(NN))
   IF (DN%VENT .OR. DN%AMBIENT .OR. DN%LEAKAGE .OR. DN%GEOM) CYCLE
   MSUM = 0._EB
   MTOT = 0._EB
   COUNT = 0
   DO ND=1,DN%N_DUCTS
      DU=>DUCT(DN%DUCT_INDEX(ND))
      IF (DU%DUCTRUN_M_INDEX > 0) THEN
         VEL = DR%VEL(DU%DUCTRUN_M_INDEX,NF,NEW)
      ELSE
         VEL = DU%VEL(NEW)
      ENDIF
      IF (ABS(DU%VEL(NEW))<1.E-5_EB) COUNT = COUNT + 1
      MFLOW = DN%DIR(ND)*VEL*DR%RHO_D(DU%DUCTRUN_M_INDEX,NF)*DU%AREA
      MSUM = MSUM + MFLOW
      MTOT = MTOT + ABS(MFLOW)
   ENDDO
   IF(ABS(MSUM)< 1.E-6 * MTOT .OR. MTOT < TWO_EPSILON_EB .OR. COUNT >= DN%N_DUCTS-1) CYCLE
   CONVERGED = .FALSE.
ENDDO

IF (CONVERGED) ITER=ITER_MAX_QFAN+1

END SUBROUTINE CONVERGENCE_CHECK_QFAN

!> \brief Updates the heat added or removed by aircoils in a ductrun
!>
!> \param DUCTRUN_INDEX Index indicating which HVAC ductrun is being solved
!> \param NF Index indicating which HVAC ductrun fan is being solved
!> \param T Current time (s)

SUBROUTINE COIL_UPDATE_QFAN(T,DUCTRUN_INDEX,NF)
USE MATH_FUNCTIONS, ONLY : EVALUATE_RAMP
USE PHYSICAL_FUNCTIONS, ONLY : GET_SPECIFIC_HEAT, GET_AVERAGE_SPECIFIC_HEAT, GET_ENTHALPY
INTEGER, INTENT(IN) :: DUCTRUN_INDEX,NF
REAL(EB), INTENT(IN) :: T
REAL(EB) :: COIL_ON_TIME,DVEL
INTEGER :: ND
TYPE(DUCT_TYPE),POINTER::DU
TYPE(DUCTRUN_TYPE),POINTER :: DR

DR => DUCTRUN(DUCTRUN_INDEX)

COIL_LOOP: DO ND = 1,DR%N_DUCTS
   DU => DUCT(DR%DUCT_INDEX(ND))
   IF (DU%AIRCOIL_INDEX < 0) CYCLE COIL_LOOP
   IF (DU%DEVC_INDEX > 0) THEN
      IF (DEVICE(DU%DEVC_INDEX)%CURRENT_STATE .AND. &
         (DEVICE(DU%DEVC_INDEX)%CURRENT_STATE .NEQV. DEVICE(DU%DEVC_INDEX)%PRIOR_STATE)) COIL_ON_TIME = T
   ENDIF
   IF (DU%CTRL_INDEX > 0) THEN
      IF (CONTROL(DU%CTRL_INDEX)%CURRENT_STATE .AND. &
         (CONTROL(DU%CTRL_INDEX)%CURRENT_STATE .NEQV. CONTROL(DU%CTRL_INDEX)%PRIOR_STATE)) COIL_ON_TIME = T
   ENDIF

   DU%COIL_Q = 0._EB
   IF (.NOT. DU%COIL_OPERATING) CYCLE COIL_LOOP
   IF (DU%AREA < TWO_EPSILON_EB) CYCLE COIL_LOOP
   IF (DU%DUCTRUN_M_INDEX >0) THEN
      DVEL = DR%VEL(ND,NF,OLD)
   ELSE
      DVEL = DU%VEL(NEW)
   ENDIF

   CALL CALC_COIL_Q(DU%AIRCOIL_INDEX,DU%COIL_Q,COIL_ON_TIME,DVEL,DR%TMP_D(ND,NF),DR%ZZ_D(ND,NF,1:N_TRACKED_SPECIES), &
                    DR%RHO_D(ND,NF),DR%CP_D(ND,NF),DU%AREA)

END DO COIL_LOOP

END SUBROUTINE COIL_UPDATE_QFAN


!> \brief Determines wall friction loss and assigns node losses to ducts
!> \param T Current time (s)
!> \param DUCTRUN_INDEX Index indicating which HVAC ductrun is being solved
!> \param NF Index indicating which HVAC ductrun dan is being solved

SUBROUTINE UPDATE_LOSS_QFAN(T,DUCTRUN_INDEX,NF)
USE PHYSICAL_FUNCTIONS,ONLY:GET_VISCOSITY
USE MATH_FUNCTIONS,ONLY:EVALUATE_RAMP
REAL(EB) :: FRICTION_FACTOR,LOSS_SUM,ZZ_GET(1:N_TRACKED_SPECIES),VISCOSITY,VFLOW,AREA,DVEL
INTEGER, INTENT(IN) :: DUCTRUN_INDEX,NF
REAL(EB), INTENT(IN) :: T
INTEGER :: ND,ND2, NN,NUM_OUT,DM_INDEX,D_INDEX
TYPE(DUCT_TYPE), POINTER :: DU,DU2
TYPE(DUCTNODE_TYPE), POINTER :: DN
TYPE(DUCTRUN_TYPE), POINTER :: DR

DR=> DUCTRUN(DUCTRUN_INDEX)
DO ND = 1, DR%N_DUCTS
   DR%LOSS(ND,NF) = 0._EB
ENDDO

NODELOOP : DO NN=1,DR%N_DUCTNODES
   VFLOW = 0._EB
   DN => DUCTNODE(DR%NODE_INDEX(NN))
   IF (DN%LEAKAGE) CYCLE

   ! Add filter loss to the downstream duct if no flow split loss over the two ducts
   IF (DN%FILTER_INDEX > 0) THEN
      IF (DUCT(DN%DUCT_INDEX(1))%AREA < TWO_EPSILON_EB .OR. DUCT(DN%DUCT_INDEX(2))%AREA < TWO_EPSILON_EB) CYCLE
      IF (FILTER(DN%FILTER_INDEX)%AREA > 0._EB) THEN
         AREA = FILTER(DN%FILTER_INDEX)%AREA
      ELSE
         AREA = 0.5_EB*(DUCT(DN%DUCT_INDEX(1))%AREA+DUCT(DN%DUCT_INDEX(2))%AREA)
      ENDIF
      DM_INDEX = DUCT(DN%DUCT_INDEX(1))%DUCTRUN_M_INDEX
      IF (DM_INDEX <=0) THEN
         DVEL = DUCT(DN%DUCT_INDEX(1))%VEL(NEW)
      ELSE
         DVEL = DR%VEL(DM_INDEX,NF,OLD)
      ENDIF
      IF(DVEL*DN%DIR(1) > 0._EB) THEN
         D_INDEX = DUCT(DN%DUCT_INDEX(2))%DUCTRUN_INDEX
         DR%LOSS(D_INDEX,NF) = DR%LOSS(D_INDEX,NF) + DN%FILTER_LOSS*(DUCT(DN%DUCT_INDEX(2))%AREA/AREA)**2
      ELSE
         D_INDEX = DUCT(DN%DUCT_INDEX(1))%DUCTRUN_INDEX
         IF (D_INDEX > 0) DR%LOSS(D_INDEX,NF) = DR%LOSS(D_INDEX,NF) + DN%FILTER_LOSS*(DUCT(DN%DUCT_INDEX(1))%AREA/AREA)**2
      ENDIF
   ENDIF

   NODECLASS: IF(DN%VENT .OR. DN%AMBIENT .OR. DN%GEOM) THEN
      ! If node is an external node loss is simply based on inflow or outflow or half loss if no flow
      IF (DUCT(DN%DUCT_INDEX(1))%AREA < TWO_EPSILON_EB .OR. DUCT(DN%DUCT_INDEX(1))%LOCALIZED_LEAKAGE) CYCLE
      D_INDEX = DUCT(DN%DUCT_INDEX(1))%DUCTRUN_INDEX
      DM_INDEX = DUCT(DN%DUCT_INDEX(1))%DUCTRUN_M_INDEX
      IF (DM_INDEX <=0) THEN
         DVEL = DUCT(DN%DUCT_INDEX(1))%VEL(NEW)
      ELSE
         DVEL = DR%VEL(DM_INDEX,NF,OLD)
      ENDIF
      IF(DVEL*DN%DIR(1) > 1.E-6_EB) THEN
         DR%LOSS(D_INDEX,NF) = DR%LOSS(D_INDEX,NF) + DN%LOSS_ARRAY(1,2)
      ELSEIF(DVEL*DN%DIR(1) < -1.E-6_EB) THEN
         DR%LOSS(D_INDEX,NF) = DR%LOSS(D_INDEX,NF) + DN%LOSS_ARRAY(2,1)
      ELSE
         DR%LOSS(D_INDEX,NF) = DR%LOSS(D_INDEX,NF) + 0.5_EB*(DN%LOSS_ARRAY(1,2)+DN%LOSS_ARRAY(2,1))
      ENDIF

   ELSEIF(DN%N_DUCTS==2) THEN NODECLASS
      D_INDEX = DUCT(DN%DUCT_INDEX(1))%DUCTRUN_INDEX
      IF (D_INDEX <=0) CYCLE NODELOOP
      DM_INDEX = DUCT(DN%DUCT_INDEX(1))%DUCTRUN_M_INDEX
      IF (DM_INDEX <=0) THEN
         DVEL = DUCT(DN%DUCT_INDEX(1))%VEL(NEW)
      ELSE
         DVEL = DR%VEL(DM_INDEX,NF,OLD)
      ENDIF
      IF (DVEL*DN%DIR(1) > 1.E-6_EB) THEN
         DR%LOSS(D_INDEX,NF) = DR%LOSS(D_INDEX,NF) + DN%LOSS_ARRAY(1,2)
      ELSEIF (DVEL*DN%DIR(1) < -1.E-6_EB) THEN
         DR%LOSS(D_INDEX,NF) = DR%LOSS(D_INDEX,NF) + DN%LOSS_ARRAY(2,1)
      ELSE
         DR%LOSS(D_INDEX,NF) = DR%LOSS(D_INDEX,NF) + 0.5_EB*(DN%LOSS_ARRAY(1,2)+DN%LOSS_ARRAY(2,1))
      ENDIF

   ELSE NODECLASS
      ! For an internal node each outlet is weights the inlet flows
      NUM_OUT = 0
      DO ND=1,DN%N_DUCTS
         DU => DUCT(DN%DUCT_INDEX(ND))
         IF (DU%AREA < TWO_EPSILON_EB) CYCLE
         DM_INDEX = DUCT(DN%DUCT_INDEX(ND))%DUCTRUN_M_INDEX
         IF (DM_INDEX <=0) THEN
            DVEL = DUCT(DN%DUCT_INDEX(ND))%VEL(NEW)
         ELSE
            DVEL = DR%VEL(DM_INDEX,NF,OLD)
         ENDIF
         IF (DVEL*DN%DIR(ND) < 0._EB) NUM_OUT = NUM_OUT + 1
      ENDDO

      NUM_OUT_IF: IF (NUM_OUT==0 .OR. NUM_OUT==DN%N_DUCTS) THEN
         ! If all are inflow or outflow each duct gets the average of its inflowing losses normalized by the number of ducts
         DO ND=1,DN%N_DUCTS
            DU => DUCT(DN%DUCT_INDEX(ND))
            LOSS_SUM = 0._EB
            DO ND2=1,DN%N_DUCTS
               IF (ND2==ND) CYCLE
               LOSS_SUM = LOSS_SUM + DN%LOSS_ARRAY(ND2,ND)
            ENDDO
            D_INDEX = DU%DUCTRUN_INDEX
            IF (D_INDEX > 0) DR%LOSS(D_INDEX,NF) = DR%LOSS(D_INDEX,NF) + LOSS_SUM / (DN%N_DUCTS-1) / DN%N_DUCTS
         ENDDO

      ELSE NUM_OUT_IF
         ! Weight the outflowing losses based on fraction of inflow volume flow
         VFLOW = 0._EB
         DO ND=1,DN%N_DUCTS
            DU => DUCT(DN%DUCT_INDEX(ND))
            DM_INDEX = DUCT(DN%DUCT_INDEX(ND))%DUCTRUN_M_INDEX
            IF (DM_INDEX <=0) THEN
               DVEL = DUCT(DN%DUCT_INDEX(ND))%VEL(NEW)
            ELSE
               DVEL = DR%VEL(DM_INDEX,NF,OLD)
            ENDIF
            IF (DVEL*DN%DIR(ND) > 0._EB) VFLOW = VFLOW + ABS(DVEL*DU%AREA)
         ENDDO
         DO ND=1,DN%N_DUCTS
            DU => DUCT(DN%DUCT_INDEX(ND))
            D_INDEX = DU%DUCTRUN_INDEX
            IF (DU%VEL(PREVIOUS)*DN%DIR(ND) > 0._EB .OR. D_INDEX <=0) CYCLE
            DO ND2=1,DN%N_DUCTS
               DU2 => DUCT(DN%DUCT_INDEX(ND2))
               DM_INDEX = DUCT(DN%DUCT_INDEX(ND2))%DUCTRUN_M_INDEX
               IF (DM_INDEX <=0) THEN
                  DVEL = DUCT(DN%DUCT_INDEX(ND2))%VEL(NEW)
               ELSE
                  DVEL = DR%VEL(DM_INDEX,NF,OLD)
               ENDIF
               IF (ND == ND2 .OR. DVEL*DN%DIR(ND2) <=0._EB) CYCLE
               DR%LOSS(D_INDEX,NF) = DR%LOSS(D_INDEX,NF) + DN%LOSS_ARRAY(ND2,ND) * ABS(DVEL*DU2%AREA)/VFLOW
            ENDDO
         ENDDO
      ENDIF NUM_OUT_IF
   ENDIF NODECLASS
ENDDO NODELOOP

DO ND = 1, DR%N_DUCTS
   DU => DUCT(DR%DUCT_INDEX(ND))
   IF (DU%AREA < TWO_EPSILON_EB) CYCLE
   IF (DU%ROUGHNESS > 0._EB) THEN
      DM_INDEX = DU%DUCTRUN_M_INDEX
      IF (DM_INDEX <=0) THEN
         DVEL = DUCT(DN%DUCT_INDEX(ND))%VEL(NEW)
      ELSE
         DVEL = DR%VEL(DM_INDEX,NF,OLD)
      ENDIF
      ZZ_GET(1:N_TRACKED_SPECIES) = DR%ZZ_D(DU%DUCTRUN_INDEX,NF,1:N_TRACKED_SPECIES)
      CALL GET_VISCOSITY(ZZ_GET,VISCOSITY,DR%TMP_D(DU%DUCTRUN_INDEX,NF))
      FRICTION_FACTOR = &
         COMPUTE_FRICTION_FACTOR(DR%RHO_D(DU%DUCTRUN_INDEX,NF),VISCOSITY,ABS(DVEL),DU%DIAMETER,DU%ROUGHNESS)
   ELSE
      FRICTION_FACTOR = 0._EB
   ENDIF
   IF (DVEL > 0._EB) THEN
      LOSS_SUM = DU%LOSS(1) * EVALUATE_RAMP(T,DU%RAMP_LOSS_INDEX)
   ELSEIF (DVEL < 0._EB) THEN
      LOSS_SUM = DU%LOSS(2) * EVALUATE_RAMP(T,DU%RAMP_LOSS_INDEX)
   ELSE
      LOSS_SUM = 0.5_EB*(DU%LOSS(1)+DU%LOSS(2)) * EVALUATE_RAMP(T,DU%RAMP_LOSS_INDEX)
   ENDIF
   DR%LOSS(ND,NF) = DR%LOSS(ND,NF) + DU%LENGTH/DU%DIAMETER*FRICTION_FACTOR + LOSS_SUM
   IF (DU%FAN_INDEX>0) THEN
      IF(.NOT. DU%FAN_OPERATING) DR%LOSS(ND,NF) = DR%LOSS(ND,NF) + FAN(DU%FAN_INDEX)%OFF_LOSS
   ENDIF
ENDDO

END SUBROUTINE UPDATE_LOSS_QFAN

!> \brief sets donor (upstream) values for ducts and nodes
!> \param DUCTRUN_INDEX Index indicating which HVAC ductrun is being solved
!> \param NF Index indicating which HVAC ductrun fan is being solved

SUBROUTINE SET_DONOR_QFAN(DUCTRUN_INDEX,NF)

INTEGER :: ND,NN
INTEGER, INTENT(IN) :: DUCTRUN_INDEX,NF
REAL(EB) :: DVEL
TYPE(DUCT_TYPE), POINTER :: DU
TYPE(DUCTNODE_TYPE), POINTER :: DN
TYPE(DUCTRUN_TYPE), POINTER :: DR

DR => DUCTRUN(DUCTRUN_INDEX)

DUCTLOOP: DO ND=1,DR%N_DUCTS
   DU=>DUCT(DR%DUCT_INDEX(ND))
   IF (DU%AREA < TWO_EPSILON_EB) CYCLE DUCTLOOP
   IF (DR%DUCT_M_INDEX(ND)<=0) THEN
      DVEL = DU%VEL(NEW)
   ELSE
      DVEL = DR%VEL(ND,NF,NEW)
   ENDIF
   IF (ABS(DVEL) > TWO_EPSILON_EB) THEN
      IF (DVEL>0._EB) THEN
         NN = DU%NODE_INDEX(1)
      ELSE
         NN = DU%NODE_INDEX(2)
      ENDIF

      DN=>DUCTNODE(NN)
      IF(DN%VENT .OR. DN%AMBIENT .OR. DN%LEAKAGE .OR. DN%GEOM) THEN
         DR%RHO_D(ND,NF)   = DN%RHO_V
         DR%TMP_D(ND,NF)   = DN%TMP_V
         DR%CP_D(ND,NF)    = DN%CP_V
         DR%ZZ_D(ND,NF,:)  = DN%ZZ_V(:)
      ELSE
         DR%RHO_D(ND,NF)  = DR%RHO_N(DN%DUCTRUN_INDEX,NF)
         DR%TMP_D(ND,NF)  = DR%TMP_N(DN%DUCTRUN_INDEX,NF)
         DR%CP_D(ND,NF)   = DR%CP_N(DN%DUCTRUN_INDEX,NF)
         DR%ZZ_D(ND,NF,:) = DR%ZZ_N(DN%DUCTRUN_INDEX,NF,:)
      ENDIF
   ELSE
      DN=>DUCTNODE(DU%NODE_INDEX(1))
      IF(DN%VENT .OR. DN%AMBIENT .OR. DN%LEAKAGE .OR. DN%GEOM) THEN
         DR%RHO_D(ND,NF)  = 0.5_EB*DN%RHO_V
         DR%TMP_D(ND,NF)  = 0.5_EB*DN%TMP_V
         DR%CP_D(ND,NF)   = 0.5_EB*DN%CP_V
         DR%ZZ_D(ND,NF,:) = 0.5_EB*DN%ZZ_V(:)
      ELSE
         DR%RHO_D(ND,NF)  = 0.5_EB*DR%RHO_N(DN%DUCTRUN_INDEX,NF)
         DR%TMP_D(ND,NF)  = 0.5_EB*DR%TMP_N(DN%DUCTRUN_INDEX,NF)
         DR%CP_D(ND,NF)   = 0.5_EB*DR%CP_N(DN%DUCTRUN_INDEX,NF)
         DR%ZZ_D(ND,NF,:) = 0.5_EB*DR%ZZ_N(DN%DUCTRUN_INDEX,NF,:)
      ENDIF
      DN=>DUCTNODE(DU%NODE_INDEX(2))
      IF(DN%VENT .OR. DN%AMBIENT .OR. DN%LEAKAGE .OR. DN%GEOM) THEN
         DR%RHO_D(ND,NF)   = DR%RHO_D(ND,NF)   + 0.5_EB*DN%RHO_V
         DR%TMP_D(ND,NF)   = DR%TMP_D(ND,NF)   + 0.5_EB*DN%TMP_V
         DR%CP_D(ND,NF)    = DR%CP_D(ND,NF)    + 0.5_EB*DN%CP_V
         DR%ZZ_D(ND,NF,:)  = DR%ZZ_D(ND,NF,:)  + 0.5_EB*DN%ZZ_V(:)
      ELSE
         DN=>DUCTNODE(DU%NODE_INDEX(2))
         DR%RHO_D(ND,NF)   = DR%RHO_D(ND,NF)   + 0.5_EB*DR%RHO_N(DN%DUCTRUN_INDEX,NF)
         DR%TMP_D(ND,NF)   = DR%TMP_D(ND,NF)   + 0.5_EB*DR%TMP_N(DN%DUCTRUN_INDEX,NF)
         DR%CP_D(ND,NF)    = DR%CP_D(ND,NF)    + 0.5_EB*DR%CP_N(DN%DUCTRUN_INDEX,NF)
         DR%ZZ_D(ND,NF,:)  = DR%ZZ_D(ND,NF,:)  + 0.5_EB*DR%ZZ_N(DN%DUCTRUN_INDEX,NF,:)
      ENDIF
   ENDIF
ENDDO DUCTLOOP

END SUBROUTINE SET_DONOR_QFAN

!> \brief Iterates over the HVAC ductrun updating node and duct quantities
!> \details The routine loops over all ducts and nodes repeatedly. If a node has all of its upstream ducts updated
!> or is a VENT inlet, then the node is updated along with its downstream ducts. The process is repeated until all ducts and
!> nodes have been updated.
!>
!> \param DUCTRUNT Current time step (s)
!> \param NNE Index indicating which HVAC network is being solved

SUBROUTINE HVAC_UPDATE_QFAN(DUCTRUN_INDEX,NF)

USE PHYSICAL_FUNCTIONS, ONLY : GET_TEMPERATURE,GET_SPECIFIC_GAS_CONSTANT,GET_ENTHALPY
REAL(EB) :: CPTSUM,ETOT,MTOT,TGUESS,VFLOW,ZZ_GET(1:N_TRACKED_SPECIES),&
            ZZTOT(1:N_TRACKED_SPECIES),HGAS,RSUM,DVEL
INTEGER, INTENT(IN) :: DUCTRUN_INDEX,NF
INTEGER :: NN,ND,ND2
LOGICAL :: CYCLE_FLAG,D_UPDATED(N_DUCTS),N_UPDATED(N_DUCTNODES)
TYPE (DUCTNODE_TYPE), POINTER :: DN,DN2
TYPE (DUCT_TYPE), POINTER :: DU
TYPE (DUCTRUN_TYPE), POINTER :: DR

DR=>DUCTRUN(DUCTRUN_INDEX)

D_UPDATED = .FALSE.
N_UPDATED = .FALSE.

ITER_LOOP: DO
   CYCLE_FLAG = .FALSE.
   DUCT_LOOP:DO ND = 1,DR%N_DUCTS
      IF (D_UPDATED(DR%DUCT_INDEX(ND))) CYCLE DUCT_LOOP
      DU=>DUCT(DR%DUCT_INDEX(ND))
      CYCLE_FLAG = .TRUE.
      IF (DU%DUCTRUN_M_INDEX <= 0) THEN ! Duct is fixed velocity
         DVEL = DU%VEL(NEW)
      ELSE
         DVEL = DR%VEL(DU%DUCTRUN_M_INDEX,NF,NEW)
      ENDIF
      IF (ABS(DVEL) > TWO_EPSILON_EB) THEN
         IF (DVEL > TWO_EPSILON_EB) NN = DU%NODE_INDEX(1)
         IF (DVEL < TWO_EPSILON_EB) NN = DU%NODE_INDEX(2)
         DN => DUCTNODE(NN)
         IF (DN%VENT .OR. DN%AMBIENT .OR. DN%LEAKAGE .OR. DN%GEOM) THEN !VENT node we can update the DR node first then the duct
            N_UPDATED(NN) = .TRUE.
            DR%RHO_N(DN%DUCTRUN_INDEX,NF)  = DN%RHO_V
            DR%TMP_N(DN%DUCTRUN_INDEX,NF)  = DN%TMP_V
            DR%CP_N(DN%DUCTRUN_INDEX,NF)   = DN%CP_V
            DR%ZZ_N(DN%DUCTRUN_INDEX,NF,:) = DN%ZZ_V(:)
         ENDIF
         IF (N_UPDATED(NN)) THEN
            DR%RHO_D(DU%DUCTRUN_INDEX,NF)  = DR%RHO_N(DN%DUCTRUN_INDEX,NF)
            DR%TMP_D(DU%DUCTRUN_INDEX,NF)  = DR%TMP_N(DN%DUCTRUN_INDEX,NF)
            DR%CP_D(DU%DUCTRUN_INDEX,NF)   = DR%CP_N(DN%DUCTRUN_INDEX,NF)
            DR%ZZ_D(DU%DUCTRUN_INDEX,NF,:) = DR%ZZ_N(DN%DUCTRUN_INDEX,NF,:)
            D_UPDATED(DR%DUCT_INDEX(ND)) = .TRUE.
            CYCLE DUCT_LOOP
         ENDIF
      ELSE
         IF (DUCTNODE(DU%NODE_INDEX(1))%VENT .OR. DUCTNODE(DU%NODE_INDEX(1))%AMBIENT .OR. DUCTNODE(DU%NODE_INDEX(1))%LEAKAGE) &
            N_UPDATED(DU%NODE_INDEX(1)) = .TRUE.
         IF (DUCTNODE(DU%NODE_INDEX(2))%VENT .OR. DUCTNODE(DU%NODE_INDEX(2))%AMBIENT .OR. DUCTNODE(DU%NODE_INDEX(2))%LEAKAGE) &
            N_UPDATED(DU%NODE_INDEX(2)) = .TRUE.
         IF (N_UPDATED(DU%NODE_INDEX(1)) .AND. N_UPDATED(DU%NODE_INDEX(2))) THEN
            DN => DUCTNODE(DU%NODE_INDEX(1))
            IF (DN%VENT .OR. DN%AMBIENT .OR. DN%LEAKAGE .OR. DN%GEOM) THEN
               DR%RHO_N(DN%DUCTRUN_INDEX,NF)  = 0.5_EB*DN%RHO_V
               DR%TMP_N(DN%DUCTRUN_INDEX,NF)  = 0.5_EB*DN%TMP_V
               DR%CP_N(DN%DUCTRUN_INDEX,NF)   = 0.5_EB*DN%CP_V
               DR%ZZ_N(DN%DUCTRUN_INDEX,NF,:) = 0.5_EB*DN%ZZ_V(:)
            ELSE
               DR%RHO_D(DU%DUCTRUN_INDEX,NF)  = 0.5_EB*DR%RHO_N(DN%DUCTRUN_INDEX,NF)
               DR%TMP_D(DU%DUCTRUN_INDEX,NF)  = 0.5_EB*DR%TMP_N(DN%DUCTRUN_INDEX,NF)
               DR%CP_D(DU%DUCTRUN_INDEX,NF)   = 0.5_EB*DR%CP_N(DN%DUCTRUN_INDEX,NF)
               DR%ZZ_D(DU%DUCTRUN_INDEX,NF,:) = 0.5_EB*DR%ZZ_N(DN%DUCTRUN_INDEX,NF,:)
            ENDIF
            DN => DUCTNODE(DU%NODE_INDEX(2))
            IF (DN%VENT .OR. DN%AMBIENT .OR. DN%LEAKAGE .OR. DN%GEOM) THEN
               DR%RHO_N(DN%DUCTRUN_INDEX,NF)  = DR%RHO_N(DN%DUCTRUN_INDEX,NF)  + 0.5_EB*DN%RHO_V
               DR%TMP_N(DN%DUCTRUN_INDEX,NF)  = DR%TMP_N(DN%DUCTRUN_INDEX,NF)  + 0.5_EB*DN%TMP_V
               DR%CP_N(DN%DUCTRUN_INDEX,NF)   = DR%CP_N(DN%DUCTRUN_INDEX,NF)   + 0.5_EB*DN%CP_V
               DR%ZZ_N(DN%DUCTRUN_INDEX,NF,:) = DR%ZZ_N(DN%DUCTRUN_INDEX,NF,:) + 0.5_EB*DN%ZZ_V(:)
            ELSE
               DR%RHO_D(DU%DUCTRUN_INDEX,NF)  = DR%RHO_D(DU%DUCTRUN_INDEX,NF)  + 0.5_EB*DR%RHO_N(DN%DUCTRUN_INDEX,NF)
               DR%TMP_D(DU%DUCTRUN_INDEX,NF)  = DR%TMP_D(DU%DUCTRUN_INDEX,NF)  + 0.5_EB*DR%TMP_N(DN%DUCTRUN_INDEX,NF)
               DR%CP_D(DU%DUCTRUN_INDEX,NF)   = DR%CP_D(DU%DUCTRUN_INDEX,NF)   + 0.5_EB*DR%CP_N(DN%DUCTRUN_INDEX,NF)
               DR%ZZ_D(DU%DUCTRUN_INDEX,NF,:) = DR%ZZ_D(DU%DUCTRUN_INDEX,NF,:) + 0.5_EB*DR%ZZ_N(DN%DUCTRUN_INDEX,NF,:)
            ENDIF
            D_UPDATED(DR%DUCT_INDEX(ND)) = .TRUE.
            CYCLE DUCT_LOOP
         ENDIF
      ENDIF
   ENDDO DUCT_LOOP

   NODE_LOOP:DO NN = 1,DR%N_DUCTNODES
      IF(N_UPDATED(DR%NODE_INDEX(NN))) CYCLE NODE_LOOP
      DN=>DUCTNODE(DR%NODE_INDEX(NN))
      CYCLE_FLAG = .TRUE.
      MTOT = 0._EB
      ETOT = 0._EB
      ZZTOT = 0._EB
      TGUESS = 0._EB
      CPTSUM = 0
      DO ND = 1,DN%N_DUCTS
         ND2 = DN%DUCT_INDEX(ND)
         IF (DUCT(ND2)%AREA < TWO_EPSILON_EB) CYCLE
         IF (DUCT(ND2)%DUCTRUN_M_INDEX <= 0) THEN
            DVEL = DUCT(ND2)%VEL(NEW)
         ELSE
            DVEL = DR%VEL(DUCT(ND2)%DUCTRUN_M_INDEX,NF,NEW)
         ENDIF
         IF (ABS(DVEL) <=TWO_EPSILON_EB) CYCLE
         IF (DVEL*DN%DIR(ND) < 0._EB) CYCLE
         IF (.NOT. D_UPDATED(ND2) .AND. ABS(DVEL) >= 1.E-5_EB) CYCLE NODE_LOOP

         VFLOW = ABS(DVEL)*DUCT(ND2)%AREA
         MTOT = MTOT + VFLOW * DR%RHO_D(DUCT(ND2)%DUCTRUN_INDEX,NF)
         ETOT = ETOT + VFLOW * DR%RHO_D(DUCT(ND2)%DUCTRUN_INDEX,NF) * DR%TMP_D(DUCT(ND2)%DUCTRUN_INDEX,NF) * &
            DR%CP_D(DUCT(ND2)%DUCTRUN_INDEX,NF)
         TGUESS = TGUESS + VFLOW * DR%RHO_D(DUCT(ND2)%DUCTRUN_INDEX,NF) * DR%TMP_D(DUCT(ND2)%DUCTRUN_INDEX,NF)
         IF (STRATIFICATION) THEN
            IF (DUCT(ND2)%NODE_INDEX(1)==DR%NODE_INDEX(NN)) THEN
               DN2=>DUCTNODE(DUCT(ND2)%NODE_INDEX(2))
            ELSE
               DN2=>DUCTNODE(DUCT(ND2)%NODE_INDEX(1))
            ENDIF
            ETOT = ETOT + DR%RHO_D(DUCT(ND2)%DUCTRUN_INDEX,NF)*VFLOW*GVEC(3)*(DN%XYZ(3)-DN2%XYZ(3))
         ENDIF
         ZZTOT(:) = ZZTOT(:) + VFLOW * DR%RHO_D(DUCT(ND2)%DUCTRUN_INDEX,NF) * DR%ZZ_D(DUCT(ND2)%DUCTRUN_INDEX,NF,:)
         ETOT = ETOT + DUCT(ND2)%COIL_Q
      ENDDO

      N_UPDATED(DR%NODE_INDEX(NN)) = .TRUE.
      IF (ABS(MTOT)<=TWO_EPSILON_EB) CYCLE NODE_LOOP

      DR%ZZ_N(NN,NF,:)  = ZZTOT(:)/MTOT
      ZZ_GET(:) = DR%ZZ_N(NN,NF,:)
      CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,RSUM)
      ETOT = ETOT/ MTOT
      DR%TMP_N(NN,NF) = TGUESS/MTOT
      CALL GET_TEMPERATURE(DR%TMP_N(NN,NF),ETOT,ZZ_GET)
      CALL GET_ENTHALPY(ZZ_GET,HGAS,DR%TMP_N(NN,NF))
      DR%CP_N(NN,NF) = HGAS/DR%TMP_N(NN,NF)
      IF (DN%DUCTRUN_M_INDEX > 0) THEN
         DR%RHO_N(NN,NF) = (DR%P(DN%DUCTRUN_M_INDEX,NF,2)+P_INF)/(RSUM*DR%TMP_N(NN,NF))
      ELSE
         DR%RHO_N(NN,NF) = (DN%P+P_INF)/(RSUM*DR%TMP_N(NN,NF))
      ENDIF
   ENDDO NODE_LOOP

   IF (.NOT. CYCLE_FLAG) EXIT ITER_LOOP

ENDDO ITER_LOOP

END SUBROUTINE HVAC_UPDATE_QFAN


!> \brief Define the index and other properties of output quantities

SUBROUTINE GET_QUANTITY_INDEX(SMOKEVIEW_LABEL,SMOKEVIEW_BAR_LABEL,OUTPUT_INDEX,Y_INDEX,Z_INDEX,UNITS,QUANTITY,SPEC_ID)
USE MISC_FUNCTIONS, ONLY: GET_SPEC_OR_SMIX_INDEX
USE OUTPUT_DATA
CHARACTER(*), INTENT(INOUT) :: QUANTITY
CHARACTER(*), INTENT(OUT) :: SMOKEVIEW_LABEL,SMOKEVIEW_BAR_LABEL,UNITS
CHARACTER(*), INTENT(IN) :: SPEC_ID
INTEGER, INTENT(OUT) :: OUTPUT_INDEX,Y_INDEX,Z_INDEX

INTEGER :: ND

! Initialize indices

Y_INDEX = -1
Z_INDEX = -1

! Look for the appropriate SPEC or SMIX index

IF (SPEC_ID/='null') THEN
   CALL GET_SPEC_OR_SMIX_INDEX(SPEC_ID,Y_INDEX,Z_INDEX)
   IF (Z_INDEX>=0  .AND. Y_INDEX>=1) Z_INDEX=-999
   IF (Z_INDEX<0 .AND. Y_INDEX<1) THEN
      WRITE(MESSAGE,'(A,A,A,A)')  'ERROR(557): SPEC_ID ',TRIM(SPEC_ID),' is not explicitly specified for QUANTITY ',TRIM(QUANTITY)
      CALL SHUTDOWN(MESSAGE) ; RETURN
   ENDIF
ENDIF

QUANTITY_INDEX_LOOP: DO ND=-N_OUTPUT_QUANTITIES,N_OUTPUT_QUANTITIES

   QUANTITY_IF:IF (QUANTITY==OUTPUT_QUANTITY(ND)%NAME) THEN

      OUTPUT_INDEX = ND

      IF (OUTPUT_QUANTITY(ND)%SPEC_ID_REQUIRED .AND. (Y_INDEX<1 .AND. Z_INDEX<0)) THEN
         IF (SPEC_ID=='null') THEN
            WRITE(MESSAGE,'(3A)')  'ERROR(558): Output QUANTITY ',TRIM(QUANTITY),' requires a SPEC_ID'
         ELSE
            WRITE(MESSAGE,'(5A)')  'ERROR(559): Output QUANTITY ',TRIM(QUANTITY),'. SPEC_ID ',TRIM(SPEC_ID),' not found.'
         ENDIF
         CALL SHUTDOWN(MESSAGE) ; RETURN
      ENDIF

      IF (.NOT. OUTPUT_QUANTITY(ND)%HVAC_SMV) THEN
         WRITE(MESSAGE,'(5A)') 'ERROR(560): HVAC_SMV QUANTITY ',TRIM(QUANTITY), ' not appropriate for .hvac file.'
         CALL SHUTDOWN(MESSAGE) ; RETURN
      ENDIF

      ! Assign Smokeview Label

      IF (Z_INDEX>=0) THEN
         SMOKEVIEW_LABEL = TRIM(SPECIES_MIXTURE(Z_INDEX)%ID)//' '//TRIM(QUANTITY)
         SMOKEVIEW_BAR_LABEL = TRIM(OUTPUT_QUANTITY(ND)%SHORT_NAME)//'_'//TRIM(SPECIES_MIXTURE(Z_INDEX)%ID)
      ELSEIF (Y_INDEX>0) THEN
         SMOKEVIEW_LABEL = TRIM(SPECIES(Y_INDEX)%ID)//' '//TRIM(QUANTITY)
         SMOKEVIEW_BAR_LABEL = TRIM(OUTPUT_QUANTITY(ND)%SHORT_NAME)//'_'//TRIM(SPECIES(Y_INDEX)%FORMULA)
      ELSE
         SMOKEVIEW_LABEL = TRIM(QUANTITY)
         SMOKEVIEW_BAR_LABEL = TRIM(OUTPUT_QUANTITY(ND)%SHORT_NAME)
      ENDIF

      UNITS = TRIM(OUTPUT_QUANTITY(ND)%UNITS)

      RETURN
   ENDIF QUANTITY_IF

ENDDO QUANTITY_INDEX_LOOP

! If no match for desired QUANTITY is found, stop the job

WRITE(MESSAGE,'(3A)') 'ERROR(561): HVAC SMV QUANTITY ',TRIM(QUANTITY), ' not found'
CALL SHUTDOWN(MESSAGE) ; RETURN

END SUBROUTINE GET_QUANTITY_INDEX

END MODULE HVAC_ROUTINES
