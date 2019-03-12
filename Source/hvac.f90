MODULE HVAC_ROUTINES

! Compute the HVAC mass and energy transport

USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_POINTERS
USE DEVICE_VARIABLES
USE CONTROL_VARIABLES
USE COMP_FUNCTIONS, ONLY: CURRENT_TIME, CHECKREAD, SHUTDOWN
USE MEMORY_FUNCTIONS, ONLY: ChkMemErr

IMPLICIT NONE

REAL(EB), ALLOCATABLE, DIMENSION(:):: NODE_AREA_EX,NODE_TMP_EX,DUCT_MF
REAL(EB), ALLOCATABLE, DIMENSION(:,:):: NODE_AREA,NODE_H,NODE_P,NODE_RHO,NODE_X,NODE_Y,NODE_Z,NODE_TMP,NODE_ZZ_EX
REAL(EB), ALLOCATABLE, DIMENSION(:,:,:):: NODE_ZZ
REAL(EB) :: DT_HV,DT_MT
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: NODE_ZONE
CHARACTER(LABEL_LENGTH), ALLOCATABLE, DIMENSION(:,:) :: NODE_DUCT_A,DUCT_NODE_A
CHARACTER(LABEL_LENGTH), ALLOCATABLE, DIMENSION(:) :: NODE_FILTER_A,DUCT_FAN_A,DUCT_AIRCOIL_A
INTEGER :: LEAK_DUCTS = 0
INTEGER, ALLOCATABLE, DIMENSION(:,:):: LEAK_PATH
CHARACTER(255) :: MESSAGE
REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: LHS
REAL(EB), ALLOCATABLE, DIMENSION(:) :: RHS,DPSTAR
INTEGER :: ITER,ITER_MAX=10

PUBLIC HVAC_CALC,READ_HVAC,PROC_HVAC,HVAC_BC_IN,FIND_NETWORKS,COLLAPSE_HVAC_BC,SET_INIT_HVAC

CONTAINS


SUBROUTINE READ_HVAC

! Read and process HVAC networks

USE MATH_FUNCTIONS, ONLY: GET_RAMP_INDEX,GET_TABLE_INDEX
USE MISC_FUNCTIONS, ONLY: SEARCH_CONTROLLER
INTEGER , PARAMETER :: MAX_DUCTS = 20
INTEGER :: IOS,IZERO,N_HVAC_READ,NS,N,NC,ND,NN,I_AIRCOIL=0,I_DUCT=0,I_DUCTNODE=0,I_FAN=0 ,I_FILTER=0,N_CELLS
REAL(EB) :: AREA,DIAMETER,XYZ(3),LOSS(MAX_DUCTS,MAX_DUCTS),VOLUME_FLOW,MAX_FLOW,MAX_PRESSURE,ROUGHNESS,LENGTH,TNOW,TAU_AC,&
            TAU_FAN,TAU_VF,FIXED_Q,CLEAN_LOSS,&
            COOLANT_MASS_FLOW,COOLANT_SPECIFIC_HEAT,COOLANT_TEMPERATURE,PERIMETER,MASS_FLOW
REAL(EB) :: LOADING(MAX_SPECIES),EFFICIENCY(MAX_SPECIES),LOADING_MULTIPLIER(MAX_SPECIES)
LOGICAL :: ROUND, SQUARE, DAMPER, REVERSE, AMBIENT,LEAK_ENTHALPY,INITIALIZED_HVAC_MASS_TRANSPORT,QFAN_BETA=.FALSE.
CHARACTER(LABEL_LENGTH) :: AIRCOIL_ID,CTRL_ID,DEVC_ID,DUCT_ID(MAX_DUCTS),DUCT_INTERP_TYPE,FAN_ID,FILTER_ID,ID,NODE_ID(2),RAMP_ID,&
                           RAMP_LOSS,SPEC_ID(MAX_SPECIES),SURF_ID,TYPE_ID,VENT_ID,VENT2_ID
TYPE(DUCTNODE_TYPE), POINTER :: DN=>NULL()
TYPE(DUCT_TYPE), POINTER :: DU=>NULL()
NAMELIST /HVAC/ AIRCOIL_ID,AMBIENT,AREA,CLEAN_LOSS,COOLANT_SPECIFIC_HEAT,COOLANT_MASS_FLOW,COOLANT_TEMPERATURE,CTRL_ID,&
                DAMPER,DEVC_ID,DIAMETER,DUCT_ID,DUCT_INTERP_TYPE,&
                EFFICIENCY,FAN_ID,FILTER_ID,FIXED_Q,ID,LEAK_ENTHALPY,LENGTH,LOADING,LOADING_MULTIPLIER,LOSS,&
                MASS_FLOW,MAX_FLOW,MAX_PRESSURE,N_CELLS,NODE_ID,PERIMETER,QFAN_BETA,&
                RAMP_ID,RAMP_LOSS,REVERSE,ROUGHNESS,SPEC_ID,SURF_ID,TAU_AC,TAU_FAN,TAU_VF,TYPE_ID,VENT_ID,VENT2_ID,VOLUME_FLOW,XYZ

TNOW=CURRENT_TIME()

N_HVAC_READ = 0

REWIND(LU_INPUT) ; INPUT_FILE_LINE_NUMBER = 0
COUNT_HVAC_LOOP: DO
   CALL CHECKREAD('HVAC',LU_INPUT,IOS)  ; IF (STOP_STATUS==SETUP_STOP) RETURN
   IF (IOS==1) EXIT COUNT_HVAC_LOOP
   READ(LU_INPUT,HVAC,END=15,ERR=16,IOSTAT=IOS)
   N_HVAC_READ = N_HVAC_READ + 1
   16 IF (IOS>0) THEN
         WRITE(MESSAGE,'(A,I5,A,I5)') &
            'ERROR: Problem with HVAC line number ',N_HVAC_READ+1,', input line number',INPUT_FILE_LINE_NUMBER
         CALL SHUTDOWN(MESSAGE); RETURN
      ENDIF
   IF (TRIM(ID)=='null') THEN
      WRITE(MESSAGE,'(A,I5,A,I5)') &
         'ERROR: No ID provided for HVAC line number ',N_HVAC_READ+1,', input line number',INPUT_FILE_LINE_NUMBER
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

IF (N_HVAC_READ > 0) HVAC_SOLVE = .TRUE.

IF ((N_DUCTS > 0 .AND. N_DUCTNODES <= 0) .OR. (N_DUCTS <= 0 .AND. N_DUCTNODES > 0)) THEN
   WRITE(MESSAGE,'(A)') 'ERROR: Must have both DUCTs and DUCTNODEs in the input file'
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
         IF (DIAMETER <= 0._EB .AND. AREA <= 0._EB .AND. PERIMETER <= 0._EB) THEN
            WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR: Duct has no AREA, DIAMETER, or PERIMTER. Duct ID:',TRIM(ID),&
                                        ', HVAC line number:',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         IF (DIAMETER > 0._EB) THEN
            IF (PERIMETER > 0._EB) THEN
               WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR: Duct cannot input both PERIMETER and DIAMETER. Duct ID:',TRIM(ID),&
                                          ', HVAC line number:',NN
               CALL SHUTDOWN(MESSAGE); RETURN
            ENDIF
            AREA = 0.5_EB*PIO2*DIAMETER**2
         ENDIF
         IF (AREA > 0._EB) THEN
            IF (DIAMETER >  0._EB .AND. PERIMETER >  0._EB) THEN
               WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR: Duct cannot input both PERIMETER and DIAMETER with AREA. Duct ID:',TRIM(ID),&
                                          ', HVAC line number:',NN
               CALL SHUTDOWN(MESSAGE); RETURN
            ENDIF
            IF (PERIMETER <= 0._EB) DIAMETER = SQRT(2._EB*AREA/PIO2)
            IF (PERIMETER >  0._EB) DIAMETER = 4._EB*AREA/PERIMETER
         ENDIF
         IF (PERIMETER > 0._EB .AND. AREA <= 0._EB) THEN
            WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR: Duct cannot have PERIMETER without AREA. Duct ID:',TRIM(ID),&
                                        ', HVAC line number:',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         DU%AREA_INITIAL = AREA
         DU%AREA = AREA
         DU%DIAMETER = DIAMETER
         DU%LENGTH = LENGTH
         DU%REVERSE = REVERSE
         ALLOCATE(DU%ZZ(N_TRACKED_SPECIES))
         DU%ZZ(1:N_TRACKED_SPECIES) = SPECIES_MIXTURE(1:N_TRACKED_SPECIES)%ZZ0
         DU%LOSS(1:2) = MAX(0._EB,LOSS(1:2,1))
         IF (CTRL_ID /='null' .AND. DEVC_ID /='null') THEN
            WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR: Can only specify one of CTRL_ID or DEVC_ID. Duct ID:',TRIM(ID),&
                                        ', HVAC line number:',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         IF (DAMPER .AND. (FAN_ID /='null' .OR. AIRCOIL_ID /='null') .OR. &
             FAN_ID/='null' .AND. (DAMPER .OR. AIRCOIL_ID /='null') .OR. &
             AIRCOIL_ID/='null' .AND. (DAMPER .OR. FAN_ID /='null')) THEN
            WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR: Duct can only have one of damper, fan or aircoil. Duct ID:',TRIM(ID),&
                                        ', HVAC line number:',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         IF (FAN_ID/='null' .AND. N_FANS<=0) THEN
            WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR: Duct has fan specied but no fans have been defined. Duct ID:',TRIM(ID),&
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
            IF (.NOT. DU%FAN_OPERATING) DU%DP_FAN = 0._EB
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
            WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR: Duct has both MASS_FLOW and VOLUME_FLOW defined. Duct ID:',TRIM(DU%ID),&
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
         IF (HVAC_MASS_TRANSPORT) THEN
            IF (N_CELLS > 0) THEN
               DU%N_CELLS = N_CELLS
            ELSE
               IF (DU%LENGTH <=0.1) THEN ! If duct is less than 10 cm, adopt two cells, else default to 10 cm cells
                  DU%N_CELLS = 2
               ELSE
                  DU%N_CELLS = NINT(DU%LENGTH*10._EB)
               END IF
            ENDIF
            DU%DX = DU%LENGTH/DU%N_CELLS
            SELECT CASE(DUCT_INTERP_TYPE) ! Reads duct interpolation type user input and sets duct interpolation index
               CASE('NODE1')
                  DU%DUCT_INTERP_TYPE_INDEX=NODE1
               CASE('NODE2')
                  DU%DUCT_INTERP_TYPE_INDEX=NODE2
               CASE('LINEAR_INTERPOLATION')
                  DU%DUCT_INTERP_TYPE_INDEX=LINEAR_INTERPOLATION
               CASE DEFAULT
                  WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR: DUCT_INTERP_TYPE is not correctly specified. Duct ID: ',TRIM(DU%ID)
                  CALL SHUTDOWN(MESSAGE); RETURN
            END SELECT
            ALLOCATE(DU%RHO_C(DU%N_CELLS))
            ALLOCATE(DU%TMP_C(DU%N_CELLS))
            ALLOCATE(DU%CP_C(DU%N_CELLS))
            ALLOCATE(DU%ZZ_C(DU%N_CELLS,N_TRACKED_SPECIES))
            DO NC = 1,DU%N_CELLS ! Initialising as background here; required for DEVC output at t = 0 s
               DU%ZZ_C(NC,1:N_TRACKED_SPECIES) = SPECIES_MIXTURE(1:N_TRACKED_SPECIES)%ZZ0
            ENDDO
         ENDIF

         IF (SURF_ID/='null') THEN
            DUCT_HT = .TRUE.
            SURF_LOOP: DO N= 1,N_SURF
               IF (SURFACE(N)%ID==SURF_ID) THEN
                  DU%SURF_INDEX = N
                  EXIT SURF_LOOP
               ENDIF
            ENDDO SURF_LOOP
            IF (DU%SURF_INDEX== -1) THEN
               WRITE(MESSAGE,'(A,A)') 'ERROR: SURFACE not found for duct ID: ',TRIM(DU%ID)
               CALL SHUTDOWN(MESSAGE); RETURN
            ENDIF
         ENDIF

      CASE('NODE')
         I_DUCTNODE = I_DUCTNODE + 1
         NODE_DUCT_A(I_DUCTNODE,:) = DUCT_ID
         NODE_FILTER_A(I_DUCTNODE) = FILTER_ID
         DN => DUCTNODE(I_DUCTNODE)
         DN%ID = ID
         DN%VENT_ID = VENT_ID
         DN%READ_IN = .TRUE.
         IF (TRIM(VENT_ID)/='null') DN%VENT=.TRUE.
         IF (.NOT. DN%VENT .AND. XYZ(3)<-1.E9) THEN
            WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR: Ambient or internal ductnode requires an elevation, XYZ(3). Ductnode ID:',&
                                        TRIM(DN%ID),', HVAC line number:',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         DN%XYZ      = XYZ
         DN%AMBIENT  = AMBIENT
         DO ND = 1, MAX_DUCTS
            IF (NODE_DUCT_A(I_DUCTNODE,ND) == 'null') EXIT
            DN%N_DUCTS=ND
         ENDDO
         IF (DN%N_DUCTS == 1 .AND. .NOT. AMBIENT .AND. VENT_ID=='null') THEN
            WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR: Non-AMBIENT or non VENT-connected ductnode must have >=2 ducts. Ductnode ID:',&
                                        TRIM(DN%ID),', HVAC line number:',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         IF (DN%N_DUCTS >= 2 .AND. (AMBIENT .OR. VENT_ID/='null')) THEN
            WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR: AMBIENT or VENT-connected ductnode must have 1 duct. Ductnode ID:',&
                                        TRIM(DN%ID),', HVAC line number:',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         IF (DN%N_DUCTS == 0) THEN
            WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR: No ducts specified for ductnode ID:',TRIM(DN%ID),', HVAC line number ',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         IF (DN%N_DUCTS > 2 .AND. TRIM(FILTER_ID)/='null') THEN
            WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR: Ductnode with a filter must have <=2 ducts. Ductnode ID:',TRIM(DN%ID),&
                                        ', HVAC line number:',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         ALLOCATE(DN%LOSS_ARRAY(MAX(2,DN%N_DUCTS),MAX(2,DN%N_DUCTS)))
         DN%LOSS_ARRAY = 0._EB
         IF (DN%N_DUCTS >=2) THEN
            DN%LOSS_ARRAY = LOSS(1:DN%N_DUCTS,1:DN%N_DUCTS)
         ELSE
            DN%LOSS_ARRAY(1,2) = LOSS(1,1)
            DN%LOSS_ARRAY(2,1) = LOSS(2,1)
         ENDIF
         IF (TRIM(FILTER_ID)/='null') THEN
            ALLOCATE(DN%FILTER_LOADING(1:N_TRACKED_SPECIES,3))
            DN%FILTER_LOADING = 0._EB
            SPEC_LOOP1: DO N=1,N_TRACKED_SPECIES
               IF (TRIM(SPEC_ID(N))=='null') EXIT SPEC_LOOP1
               DO NS = 1,N_TRACKED_SPECIES
                  IF (TRIM(SPECIES_MIXTURE(NS)%ID)==TRIM(SPEC_ID(N))) THEN
                     DN%FILTER_LOADING(NS,1) = LOADING(N)
                     EXIT
                  ENDIF
                  IF (NS==N_TRACKED_SPECIES) THEN
                     WRITE(MESSAGE,'(A,A,A,A,A)') 'ERROR: Problem with ductnode:',TRIM(ID),' SPEC ',TRIM(SPEC_ID(N)),' not found'
                     CALL SHUTDOWN(MESSAGE); RETURN
                  ENDIF
               ENDDO
            ENDDO SPEC_LOOP1
         ENDIF

      CASE('FAN')
         I_FAN = I_FAN + 1
         FAN(I_FAN)%ID = ID
         FAN(I_FAN)%OFF_LOSS = LOSS(1,1)
         FAN(I_FAN)%FAN_RAMP = RAMP_ID
         FAN(I_FAN)%VOL_FLOW = VOLUME_FLOW
         FAN(I_FAN)%MAX_PRES = MAX_PRESSURE
         FAN(I_FAN)%MAX_FLOW = MAX_FLOW
         FAN(I_FAN)%TAU = TAU_FAN
         IF (TAU_FAN > 0._EB) FAN(I_FAN)%SPIN_INDEX = TANH_RAMP
         IF (TAU_FAN < 0._EB) FAN(I_FAN)%SPIN_INDEX = TSQR_RAMP
         IF (RAMP_ID /= 'null') CALL GET_RAMP_INDEX(RAMP_ID,'FAN',FAN(I_FAN)%RAMP_INDEX)
         IF(( (MAX_FLOW<1.E6_EB .OR. MAX_PRESSURE<1.E6_EB) .AND. (VOLUME_FLOW<1.E6_EB .OR. RAMP_ID/='null')))THEN !.OR. &
            WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR: FAN can only be one of constant volume, quadratic or ramp. Fan ID:',TRIM(ID),&
                                        ', HVAC line number:',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         IF ((MAX_PRESSURE<1.E6_EB .AND. MAX_FLOW>1.E6_EB) .OR. (MAX_PRESSURE>1.E6_EB .AND. MAX_FLOW<1.E6_EB)) THEN
            WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR: IF one of MAX_PRESSURE or MAX_FLOW given, both must be specified. Fan ID:',&
                                        TRIM(ID),', HVAC line number:',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         IF (MAX_PRESSURE <= 0._EB) THEN
            WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR: MAX_PRESSURE must be > 0. Fan ID:',TRIM(ID),', HVAC line number:',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         IF (MAX_FLOW <= 0._EB) THEN
            WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR: MAX_FLOW must be > 0. Fan ID:',TRIM(ID),', HVAC line number:',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         IF (VOLUME_FLOW < 1.E6_EB) THEN
            FAN(I_FAN)%FAN_TYPE = 1
         ELSEIF(RAMP_ID/='null') THEN
            FAN(I_FAN)%FAN_TYPE = 3
         ELSE
            FAN(I_FAN)%FAN_TYPE = 2
         ENDIF
         IF(FAN(I_FAN)%FAN_TYPE == 2 .AND. QFAN_BETA .EQV. .TRUE.) THEN
            FAN(I_FAN)%FAN_TYPE = 4
            QFAN_BETA_TEST=.TRUE.
         ENDIF

      CASE('FILTER')
         I_FILTER = I_FILTER + 1
         FILTER(I_FILTER)%ID = ID
         FILTER(I_FILTER)%CLEAN_LOSS = CLEAN_LOSS
         IF (TRIM(RAMP_ID)/='null') CALL GET_RAMP_INDEX(RAMP_ID,'FILTER',FILTER(I_FILTER)%RAMP_INDEX)
         ALLOCATE(FILTER(I_FILTER)%EFFICIENCY(1:N_TRACKED_SPECIES))
         FILTER(I_FILTER)%EFFICIENCY = 0._EB
         ALLOCATE(FILTER(I_FILTER)%MULTIPLIER(1:N_TRACKED_SPECIES))
         FILTER(I_FILTER)%MULTIPLIER = 0._EB
         FILTER(I_FILTER)%LOADING_LOSS = LOSS(1,1)
         SPEC_LOOP2: DO N=1,N_TRACKED_SPECIES
            IF (TRIM(SPEC_ID(N))=='null') EXIT SPEC_LOOP2
            DO NS = 1,N_TRACKED_SPECIES
               IF (TRIM(SPECIES_MIXTURE(NS)%ID)==TRIM(SPEC_ID(N))) THEN
                  FILTER(I_FILTER)%EFFICIENCY(NS)   = EFFICIENCY(N)
                  FILTER(I_FILTER)%MULTIPLIER(NS) = LOADING_MULTIPLIER(N)
                  EXIT
               ENDIF
               IF (NS==N_TRACKED_SPECIES) THEN
                  WRITE(MESSAGE,'(A,A,A,A,A)') 'ERROR: Problem with filter:',TRIM(ID),'. SPEC ',TRIM(SPEC_ID(N)),' not found'
                  CALL SHUTDOWN(MESSAGE); RETURN
               ENDIF
            ENDDO
         ENDDO SPEC_LOOP2
      CASE('AIRCOIL')
         I_AIRCOIL = I_AIRCOIL+1
         AIRCOIL(I_AIRCOIL)%COOLANT_SPECIFIC_HEAT   = COOLANT_SPECIFIC_HEAT*1000._EB
         AIRCOIL(I_AIRCOIL)%COOLANT_MASS_FLOW = COOLANT_MASS_FLOW
         AIRCOIL(I_AIRCOIL)%COOLANT_TEMPERATURE = COOLANT_TEMPERATURE+TMPM
         AIRCOIL(I_AIRCOIL)%EFFICIENCY   = EFFICIENCY(1)
         AIRCOIL(I_AIRCOIL)%FIXED_Q      = FIXED_Q*1000._EB
         AIRCOIL(I_AIRCOIL)%ID           = ID
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
         DN%VENT_ID = VENT_ID
         DN%VENT=.TRUE.
         DN%READ_IN = .FALSE.
         IF (TRIM(DN%VENT_ID)=='null') THEN
            WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR: Leakage path must have VENT_ID defined. Leak ID:',TRIM(ID),', HVAC line number:',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         IF (TRIM(DN%VENT_ID)=='AMBIENT') THEN
            WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR: Leakage to AMBIENT must have VENT2_ID for the AMBIENT node. Leak ID:',TRIM(ID),&
                                        ', HVAC line number:',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         DN%XYZ      = XYZ
         DN%AMBIENT  = .FALSE.
         DN%N_DUCTS=1
         ALLOCATE(DN%LOSS_ARRAY(2,2))
         DN%LOSS_ARRAY = 0._EB

         I_DUCTNODE = I_DUCTNODE + 1
         NODE_DUCT_A(I_DUCTNODE,1) = ID
         NODE_FILTER_A(I_DUCTNODE) = 'null'
         DN => DUCTNODE(I_DUCTNODE)
         DN%ID = VENT2_ID
         DN%VENT_ID = VENT2_ID
         DN%VENT=.TRUE.
         DN%READ_IN = .FALSE.
         IF (TRIM(VENT2_ID)=='null') THEN
            WRITE(MESSAGE,'(A,A,A,I2)') 'ERROR: Leakage path must have VENT2_ID defined. Leak ID:',TRIM(ID),&
                                        ', HVAC line number:',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
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
         IF (AREA <= 0._EB) THEN
            WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR: Leakage has no AREA. Leak ID:',TRIM(ID),', HVAC line number:',NN
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
         DU%AREA_INITIAL = AREA
         DU%AREA = AREA
         DU%DIAMETER = -1._EB
         DU%LENGTH = 0.1_EB
         DU%REVERSE = .FALSE.
         ALLOCATE(DU%ZZ(N_TRACKED_SPECIES))
         DU%ZZ(1:N_TRACKED_SPECIES) = SPECIES_MIXTURE(1:N_TRACKED_SPECIES)%ZZ0
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
   END SELECT
ENDDO

T_USED(13)=T_USED(13)+CURRENT_TIME()-TNOW

RETURN

CONTAINS

SUBROUTINE SET_HVAC_DEFAULTS

AIRCOIL_ID   = 'null'
AMBIENT      = .FALSE.
AREA         = -1._EB
COOLANT_SPECIFIC_HEAT   = 4.186_EB
COOLANT_MASS_FLOW = -1.E10_EB
COOLANT_TEMPERATURE = 20._EB
EFFICIENCY   = 1.0_EB
CLEAN_LOSS   = 0._EB
CTRL_ID      = 'null'
DAMPER       = .FALSE.
DEVC_ID      = 'null'
DIAMETER     = -1._EB
DUCT_ID      = 'null'
DUCT_INTERP_TYPE = 'NODE1'
FAN_ID       = 'null'
FIXED_Q      = -1.E10_EB
FILTER_ID    = 'null'
LEAK_ENTHALPY = .FALSE.
LENGTH       = 1._EB
LOADING      = 0._EB
LOADING_MULTIPLIER = 1._EB
LOSS         = 0._EB
MASS_FLOW    = 1.E7_EB
MAX_FLOW     = 1.E7_EB
MAX_PRESSURE = 1.E7_EB
N_CELLS      = -999
NODE_ID      = 'null'
INITIALIZED_HVAC_MASS_TRANSPORT=.FALSE.
PERIMETER    = -1._EB
RAMP_ID      = 'null'
RAMP_LOSS    = 'null'
REVERSE      = .FALSE.
ROUGHNESS    = -1._EB
ROUND        = .TRUE.
SPEC_ID      = 'null'
SQUARE       = .FALSE.
SURF_ID      = 'null'
TYPE_ID      = 'null'
TAU_AC       = TAU_DEFAULT
TAU_FAN      = TAU_DEFAULT
TAU_VF       = TAU_DEFAULT
VENT_ID      = 'null'
VENT2_ID     = 'null'
VOLUME_FLOW  = 1.E7_EB
XYZ          = -1.E10_EB

RETURN

END SUBROUTINE SET_HVAC_DEFAULTS

END SUBROUTINE READ_HVAC


SUBROUTINE PROC_HVAC
USE PHYSICAL_FUNCTIONS, ONLY: GET_ENTHALPY
!INTEGER :: I1,I2,J1,J2,K1,K2,IOR
INTEGER :: N,ND,ND2,NM,NN,NF,NV
REAL(EB) :: TNOW,ZZ_GET(1:N_TRACKED_SPECIES)
LOGICAL :: FOUND
TYPE (LAGRANGIAN_PARTICLE_CLASS_TYPE),DIMENSION(:), POINTER:: TEMPALLOC=>NULL()
TYPE(DUCTNODE_TYPE), POINTER :: DN=>NULL()
TYPE(DUCT_TYPE), POINTER :: DU=>NULL()
TYPE(SURFACE_TYPE), POINTER :: SF=>NULL()
TNOW=CURRENT_TIME()

IF (LEAK_DUCTS > 0) THEN
   HVAC_SOLVE = .TRUE.
   CALL LEAKAGE_HVAC
ENDIF

DUCT_LOOP: DO ND = 1, N_DUCTS
   DU => DUCT(ND)
   IF (DU%LEAKAGE) CYCLE DUCT_LOOP
   IF (TRIM(DUCT_NODE_A(ND,1))==TRIM(DUCT_NODE_A(ND,2))) THEN
      WRITE(MESSAGE,'(A,A)') 'ERROR: Both nodes have the same ID. Duct ID:',TRIM(DU%ID)
      CALL SHUTDOWN(MESSAGE); RETURN
   ENDIF
   DO N = 1, ND
      IF (N==ND) CYCLE
      IF (TRIM(DU%ID)==TRIM(DUCT(N)%ID)) THEN
         WRITE(MESSAGE,'(A,A)') 'ERROR: Two ducts with the same ID. Duct ID:',TRIM(DU%ID)
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
            WRITE(MESSAGE,'(A,A,A,A,A)') 'ERROR: Duct: ',TRIM(DU%ID),'. Node:',&
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
            WRITE(MESSAGE,'(A,A,A,A,A)') 'ERROR: Duct: ',TRIM(DU%ID),'. Node:',&
                                      TRIM(DUCTNODE(NN)%ID),' does not contain the duct in its list of ducts.'
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
      ENDIF
      IF (DU%NODE_INDEX(1) > 0 .AND. DU%NODE_INDEX(2) > 0) EXIT NODELOOP_D
   ENDDO NODELOOP_D
   IF (DU%NODE_INDEX(1) <= 0) THEN
      WRITE(MESSAGE,'(A,I3,A,A)') 'ERROR: First ductnode not located for duct:',ND,', Duct ID:',TRIM(DU%ID)
      CALL SHUTDOWN(MESSAGE); RETURN
   ENDIF
   IF (DU%NODE_INDEX(2) <= 0) THEN
      WRITE(MESSAGE,'(A,I3,A,A)') 'ERROR: Second ductnode not located for duct:',ND,', Duct ID:',TRIM(DU%ID)
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
         WRITE(MESSAGE,'(A,I3,A,A)') 'ERROR: Fan not located. Duct:',ND,', Duct ID:',TRIM(DU%ID)
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
         WRITE(MESSAGE,'(A,I3,A,A)') 'ERROR: Aircoil not located. Duct:',ND,', Duct ID:',TRIM(DU%ID)
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
         WRITE(MESSAGE,'(A,A)') 'ERROR: Two duct nodes with the same ID. Ductnode ID:',TRIM(DN%ID)
         CALL SHUTDOWN(MESSAGE); RETURN
      ENDIF
   ENDDO

   ! Initialises duct node species and RSUM with ambient/background
   ALLOCATE(DN%ZZ(N_TRACKED_SPECIES))
   DN%ZZ(1:N_TRACKED_SPECIES) = SPECIES_MIXTURE(1:N_TRACKED_SPECIES)%ZZ0
   ALLOCATE(DN%ZZ_V(N_TRACKED_SPECIES))
   DN%ZZ_V(1:N_TRACKED_SPECIES) = SPECIES_MIXTURE(1:N_TRACKED_SPECIES)%ZZ0
   ZZ_GET(1:N_TRACKED_SPECIES) = DN%ZZ_V(1:N_TRACKED_SPECIES)
   DN%RSUM   = RSUM0

   ! If node is LEAKAGE related then values are adopted as ambient/background
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
   IF (DN%VENT_ID /= 'null') THEN
      ALLOCATE(DN%IN_MESH(NMESHES))
      DN%IN_MESH=.FALSE.
      FOUND = .FALSE.
      MESH_LOOP: DO NM = 1, NMESHES
         IF (EVACUATION_ONLY(NM)) CYCLE
         NODE_VENT_LOOP:DO NV = 1, MESHES(NM)%N_VENT
            IF(MESHES(NM)%VENTS(NV)%ID == DN%VENT_ID) THEN
               FOUND = .TRUE.
               IF (DN%MESH_INDEX > 0) THEN
                  WRITE(MESSAGE,'(A,A)') 'ERROR: VENT for ductnode is split over more than one mesh for VENT ID:',&
                                          TRIM(MESHES(NM)%VENTS(NV)%ID)
                  CALL SHUTDOWN(MESSAGE); RETURN
               ENDIF
               IF (DN%READ_IN .AND. MESHES(NM)%VENTS(NV)%SURF_INDEX /= HVAC_SURF_INDEX) THEN
                  WRITE(MESSAGE,'(A,A)') 'ERROR: Ductnode attached to VENT without SURF_ID HVAC for VENT ID:',&
                                          TRIM(MESHES(NM)%VENTS(NV)%ID)
                  CALL SHUTDOWN(MESSAGE); RETURN
               ENDIF
               IF (MESHES(NM)%VENTS(NV)%BOUNDARY_TYPE/=HVAC_BOUNDARY) THEN
                  SF => SURFACE(MESHES(NM)%VENTS(NV)%SURF_INDEX)
                  IF (ABS(SF%VEL)>TWO_EPSILON_EB .OR. ABS(SF%VOLUME_FLOW)>TWO_EPSILON_EB .OR. &
                      ABS(SF%MASS_FLUX_TOTAL)>TWO_EPSILON_EB .OR. SF%PYROLYSIS_MODEL/= PYROLYSIS_NONE) THEN
                      WRITE(MESSAGE,'(A,A)') 'Cannot leak and specify flow or pyrolysis at the same time.  VENT ID:',&
                                             TRIM(MESHES(NM)%VENTS(NV)%ID)
                      CALL SHUTDOWN(MESSAGE); RETURN
                  ENDIF
                  IF (ANY(SF%LEAK_PATH>0)) THEN
                      WRITE(MESSAGE,'(A,A)') 'Cannot specify custom leakage and zone leakage with the same surface.  VENT ID:',&
                                             TRIM(MESHES(NM)%VENTS(NV)%ID)
                      CALL SHUTDOWN(MESSAGE); RETURN
                  ENDIF
               ENDIF
               DN%IN_MESH(NM) = .TRUE.
               MESHES(NM)%VENTS(NV)%NODE_INDEX=NN
               EXIT NODE_VENT_LOOP
            ENDIF
         ENDDO NODE_VENT_LOOP
      ENDDO MESH_LOOP
      IF (.NOT. FOUND) THEN
         WRITE(MESSAGE,'(A,A,A,A)') 'Cannot find VENT_ID: ',TRIM(DN%VENT_ID),' for Ductnode: ',TRIM(DN%ID)
         CALL SHUTDOWN(MESSAGE); RETURN
      ENDIF
   ENDIF
   IF (DN%VENT .AND. DN%AMBIENT) THEN
      WRITE(MESSAGE,'(A,I5,A,A)') 'ERROR: DUCTNODE cannot be AMBIENT and have an assigned VENT_ID. Ductnode:',NN,&
                                  ', Ductnode ID:',TRIM(DN%ID)
      CALL SHUTDOWN(MESSAGE); RETURN
   ENDIF
   IF (DN%N_DUCTS == 1 .AND. .NOT. DN%VENT .AND. .NOT. DN%AMBIENT) THEN
      WRITE(MESSAGE,'(A,I5,A,A)') 'ERROR: Internal DUCTNODE must have at least two attached ducts. Ductnode:',NN,&
                                  ', Ductnode ID:',TRIM(DN%ID)
      CALL SHUTDOWN(MESSAGE); RETURN
   ENDIF
   IF (DN%N_DUCTS> 1 .AND. (DN%AMBIENT .OR. DN%VENT) ) THEN
      WRITE(MESSAGE,'(A,I5,A,A)') 'ERROR: External DUCTNODE can only have one attached duct. Ductnode:',NN,&
                                  ', Ductnode ID:',TRIM(DN%ID)
      CALL SHUTDOWN(MESSAGE); RETURN
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
         WRITE(MESSAGE,'(A,I5,A,I5,A,A)') 'ERROR: DUCT ',ND,' not found. Ductnode:',NN,', Ductnode ID:',TRIM(DN%ID)
         CALL SHUTDOWN(MESSAGE); RETURN
      ENDIF
      IF (TRIM(DN%ID) /= TRIM(DUCT_NODE_A(DN%DUCT_INDEX(ND),1)) .AND. TRIM(DN%ID) /= TRIM(DUCT_NODE_A(DN%DUCT_INDEX(ND),2))) THEN
         WRITE(MESSAGE,'(A,A,A,A,A)') 'ERROR: Duct: ',TRIM(DUCT(DN%DUCT_INDEX(ND))%ID),' does not contain Ductnode:',&
                                       TRIM(DN%ID), 'in its list of ductnodes.'
         CALL SHUTDOWN(MESSAGE); RETURN
      ENDIF
   ENDDO

   ! Initialize duct node properties
   IF (STRATIFICATION .AND. DN%XYZ(3) > -1.E9_EB) THEN
      DN%TMP  = TMPA + LAPSE_RATE*DN%XYZ(3)
      IF (ABS(LAPSE_RATE)>TWO_EPSILON_EB) THEN
         DN%P = P_INF*(DN%TMP/TMPA)**(GVEC(3)/RSUM0/LAPSE_RATE)
      ELSE
         DN%P = P_INF*EXP(GVEC(3)*(DN%XYZ(3)-GROUND_LEVEL)/(RSUM0*TMPA))
      ENDIF
      DN%RHO  =  DN%P/(DN%TMP*RSUM0)
   ELSE
      DN%TMP  = TMPA
      DN%P    = P_INF
      DN%RHO  = RHOA
   ENDIF
   DN%P_OLD = DN%P
   IF (DN%VENT) DN%P = -1.E10_EB
   ZZ_GET(1:N_TRACKED_SPECIES) = DN%ZZ_V(1:N_TRACKED_SPECIES)
   CALL GET_ENTHALPY(ZZ_GET,DN%CP,DN%TMP)
   DN%CP_V   = DN%CP / DN%TMP
   DN%TMP_V  = DN%TMP
   DN%RSUM_V = DN%RSUM
   DN%RHO_V  = DN%RHO
   IF(TRIM(NODE_FILTER_A(NN))/='null') THEN
      DO N = 1,N_FILTERS
         IF(TRIM(NODE_FILTER_A(NN))==TRIM(FILTER(N)%ID)) THEN
            DN%FILTER_INDEX = N
            EXIT
         ENDIF
         IF (N==N_FILTERS) THEN
            WRITE(MESSAGE,'(A,A,A,A,A)') 'ERROR: Problem with ductnode:',TRIM(DN%ID), &
                                         ',FILTER ',TRIM(NODE_FILTER_A(NN)),' not found'
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
      ENDDO
   ENDIF
ENDDO NODE_LOOP

!Temp arrays for input processing
IF (ALLOCATED(DUCT_NODE_A)) DEALLOCATE(DUCT_NODE_A)
IF (ALLOCATED(NODE_DUCT_A)) DEALLOCATE(NODE_DUCT_A)
IF (ALLOCATED(NODE_FILTER_A)) DEALLOCATE(NODE_FILTER_A)
IF (ALLOCATED(DUCT_FAN_A)) DEALLOCATE(DUCT_FAN_A)

CALL DETERMINE_FIXED_ELEMENTS(0._EB)

ALLOCATE(NODE_P(1:N_DUCTNODES,1:NMESHES))
ALLOCATE(NODE_TMP(1:N_DUCTNODES,1:NMESHES))
ALLOCATE(NODE_TMP_EX(1:N_DUCTNODES))
NODE_TMP_EX=0._EB
ALLOCATE(NODE_RHO(1:N_DUCTNODES,1:NMESHES))
ALLOCATE(NODE_H(1:N_DUCTNODES,1:NMESHES))
ALLOCATE(NODE_X(1:N_DUCTNODES,1:NMESHES))
ALLOCATE(NODE_Y(1:N_DUCTNODES,1:NMESHES))
ALLOCATE(NODE_Z(1:N_DUCTNODES,1:NMESHES))
ALLOCATE(NODE_ZONE(1:N_DUCTNODES,1:NMESHES))
ALLOCATE(NODE_AREA(1:N_DUCTNODES,1:NMESHES))
ALLOCATE(NODE_AREA_EX(1:N_DUCTNODES))
ALLOCATE(NODE_ZZ(1:N_DUCTNODES,1:N_TRACKED_SPECIES,1:NMESHES))
ALLOCATE(NODE_ZZ_EX(1:N_DUCTNODES,1:N_TRACKED_SPECIES))

NODE_X = 0._EB
NODE_Y = 0._EB
NODE_Z = 0._EB
NODE_ZONE = 0
NODE_H = 0._EB
NODE_P = 0._EB
NODE_TMP= 0._EB
NODE_RHO = 0._EB
NODE_AREA = 0._EB
NODE_AREA_EX = 0._EB
NODE_ZZ = 0._EB
NODE_ZZ_EX = 0._EB

ALLOCATE(DUCT_MF(1:N_DUCTS))
DUCT_MF = 0._EB

T_USED(13)=T_USED(13)+CURRENT_TIME()-TNOW

RETURN

END SUBROUTINE PROC_HVAC

SUBROUTINE HVAC_CALC(T,DT,FIRST_PASS)
! Solve for flows in the HVAC networks
INTEGER :: NNE,NN,NR
REAL(EB), INTENT(IN) :: T,DT
LOGICAL :: CHANGE=.TRUE.
LOGICAL, SAVE :: INITIALIZED_HVAC_MASS_TRANSPORT
LOGICAL, INTENT(IN):: FIRST_PASS
TYPE(NETWORK_TYPE), POINTER:: NE=>NULL()
TYPE(DUCTRUN_TYPE), POINTER :: DR=>NULL()

DT_HV = DT
DT_MT = DT

IF (CORRECTOR) THEN
   DUCT%VEL(OLD) = DUCT%VEL(NEW)
   DUCT%VEL(PREVIOUS) = DUCT%VEL(NEW)
   DUCT%VEL(GUESS) = DUCT%VEL(NEW)
   DUCT%DP_FAN(OLD) = DUCT%DP_FAN(NEW)
   DO NN=1,N_DUCTNODES
      IF(DUCTNODE(NN)%FILTER_INDEX > 0) DUCTNODE(NN)%FILTER_LOADING(:,OLD)=DUCTNODE(NN)%FILTER_LOADING(:,NEW)
      DUCTNODE(NN)%P_OLD = DUCTNODE(NN)%P
   ENDDO
   RETURN
ENDIF

IF (FIRST_PASS) THEN
   CALL COLLAPSE_HVAC_BC
   IF (LEAK_DUCTS > 0) CALL ADJUST_LEAKAGE_AREA
   CALL FIND_NETWORKS(CHANGE,T) ! calls determined fixed elements (which calls update fan for fixed fans)
   IF (.NOT. INITIALIZED_HVAC_MASS_TRANSPORT) CALL FIND_DUCTRUNS ! short term hack to get to run once, requires changing post-BETA
   IF (.NOT. INITIALIZED_HVAC_MASS_TRANSPORT .AND. HVAC_MASS_TRANSPORT) CALL SET_INIT_HVAC_MASS_TRANSPORT
   INITIALIZED_HVAC_MASS_TRANSPORT=.TRUE.
ENDIF

ITER = 0

DO NNE = 1, N_NETWORKS
   NE =>NETWORK(NNE)
   CALL SET_GUESS(NNE,T)
ENDDO

IF (N_ZONE >0) ALLOCATE(DPSTAR(1:N_ZONE))

CALL DPSTARCALC
DUCTNODE%P = DUCTNODE%P - P_INF

IF (QFAN_BETA_TEST) THEN
   ! If a ductrun has a quadratic fan(s) and requires system curve steady state solution: allocate, zero, populate and solve matrices
   ! BETA WARNING: system curve is calculated based on total ductrun, this ignores dampers opening and closing. Need to fix
   DO NR = 1, N_DUCTRUNS
      DR => DUCTRUN(NR)
      IF (DR%N_QFANS > 0) THEN
         ALLOCATE(DR%LHS_SYSTEM_1(DR%N_MATRIX_SYSTEM,DR%N_MATRIX_SYSTEM)) ! LHS of first point of system curve
         ALLOCATE(DR%LHS_SYSTEM_2(DR%N_MATRIX_SYSTEM,DR%N_MATRIX_SYSTEM)) ! LHS of second point of system curve
         ALLOCATE(DR%RHS_SYSTEM_1(DR%N_MATRIX_SYSTEM)) ! RHS of first point of system curve
         ALLOCATE(DR%RHS_SYSTEM_2(DR%N_MATRIX_SYSTEM,DR%N_QFANS)) ! RHS of first point of system curve
         DR%LHS_SYSTEM_1(:,:) = 0._EB
         DR%LHS_SYSTEM_2(:,:) = 0._EB
         DR%RHS_SYSTEM_1(:) = 0._EB
         DR%RHS_SYSTEM_2(:,:) = 0._EB
         CALL RHS_SYSTEM(NR)
         CALL LHS_SYSTEM(NR)
         CALL MATRIX_SYSTEM_SOLVE(NR)
         ! Deallocate matrices used for solving steady state system curve
         DEALLOCATE(DR%LHS_SYSTEM_1)
         DEALLOCATE(DR%LHS_SYSTEM_2)
         DEALLOCATE(DR%RHS_SYSTEM_1)
         DEALLOCATE(DR%RHS_SYSTEM_2)
      ENDIF
   ENDDO
ENDIF

! Solution of the HVAC network
DO NNE = 1, N_NETWORKS
   NE =>NETWORK(NNE)
   IF (NE%N_MATRIX > 0) THEN
      ITER = 0
      ALLOCATE(LHS(NE%N_MATRIX,NE%N_MATRIX))
      ALLOCATE(RHS(NE%N_MATRIX))
      DO WHILE (ITER < ITER_MAX)
         LHS = 0._EB
         RHS = 0._EB
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
         IF (ITER < ITER_MAX) CALL SET_GUESS(NNE,T)
      ENDDO
      DEALLOCATE(LHS)
      DEALLOCATE(RHS)
   ELSE
      CALL SET_DONOR(NNE)
      IF (N_AIRCOILS > 0) CALL COIL_UPDATE(T)
      CALL HVAC_UPDATE(NNE,DT)
   ENDIF
ENDDO

DUCTNODE%P = DUCTNODE%P + P_INF

IF (HVAC_MASS_TRANSPORT) CALL UPDATE_HVAC_MASS_TRANSPORT(DT_MT)

CALL UPDATE_NODE_BC

IF (ALLOCATED(DPSTAR)) DEALLOCATE(DPSTAR)


END SUBROUTINE HVAC_CALC


SUBROUTINE MATRIX_SOLVE(NNE)
USE MATH_FUNCTIONS,ONLY : GAUSSJ
INTEGER :: NNE,IERR,ND,NN
TYPE(NETWORK_TYPE), POINTER :: NE=>NULL()
TYPE(DUCT_TYPE), POINTER :: DU=>NULL()
TYPE(DUCTNODE_TYPE), POINTER :: DN=>NULL()

NE =>NETWORK(NNE)
CALL GAUSSJ(LHS,NE%N_MATRIX,NE%N_MATRIX,RHS,1,1,IERR)
DO ND = 1,NE%N_DUCTS
   DU=>DUCT(NE%DUCT_INDEX(ND))
   IF (DU%FIXED .OR. DU%AREA < TWO_EPSILON_EB) CYCLE
   DU%VEL(NEW) = RHS(NE%MATRIX_INDEX(ND))
ENDDO
DO NN = 1,NE%N_DUCTNODES
   DN=>DUCTNODE(NE%NODE_INDEX(NN))
   IF (DN%FIXED) CYCLE
   DN%P = RHS(NE%MATRIX_INDEX(NE%N_DUCTS+NN))
ENDDO

END SUBROUTINE MATRIX_SOLVE

SUBROUTINE MATRIX_SYSTEM_SOLVE(RN)
USE MATH_FUNCTIONS,ONLY : GAUSSJ
INTEGER :: RN,IERR,ND,FN
TYPE(DUCTRUN_TYPE), POINTER :: DR=>NULL()
TYPE(DUCT_TYPE), POINTER :: DU=>NULL()

DR =>DUCTRUN(RN)
CALL GAUSSJ(DR%LHS_SYSTEM_1,DR%N_MATRIX_SYSTEM,DR%N_MATRIX_SYSTEM,DR%RHS_SYSTEM_1,1,1,IERR)

DO FN = 1, DR%N_QFANS
   CALL GAUSSJ(DR%LHS_SYSTEM_2,DR%N_MATRIX_SYSTEM,DR%N_MATRIX_SYSTEM,DR%RHS_SYSTEM_2(:,FN),1,1,IERR)
ENDDO

DO ND = 1,DR%N_DUCTS
   DU=>DUCT(DR%DUCT_INDEX(ND))
   IF (DU%FIXED .OR. DU%AREA < TWO_EPSILON_EB) CYCLE
   DU%VEL_SYSTEM(1,1,1) = DU%VEL_SYSTEM(1,1,2) ! vel_system(sys#,fan#,old/new)
   DU%VEL_SYSTEM(1,1,2) = DR%RHS_SYSTEM_1(DR%MATRIX_SYSTEM_INDEX(ND))
   DO FN = 1, DR%N_QFANS
      DU%VEL_SYSTEM(2,FN,1) = DU%VEL_SYSTEM(2,FN,2)
      DU%VEL_SYSTEM(2,FN,2) = DR%RHS_SYSTEM_2(DR%MATRIX_SYSTEM_INDEX(ND),FN)
   ENDDO
ENDDO


END SUBROUTINE MATRIX_SYSTEM_SOLVE

SUBROUTINE HVAC_UPDATE(NNE,DT)

!Iterate duct network to update all ducts and nodes
USE COMP_FUNCTIONS, ONLY: CURRENT_TIME
USE PHYSICAL_FUNCTIONS, ONLY : GET_AVERAGE_SPECIFIC_HEAT,GET_SPECIFIC_GAS_CONSTANT,GET_ENTHALPY
REAL(EB) :: CP,CP2,CPTSUM,DCPDT,DU_DX,ETOT,MFLOW,MSUM,MTOT,TGUESS,TNOW,VFLOW,ZZ_GET(1:N_TRACKED_SPECIES),&
            ZZSUM(1:N_TRACKED_SPECIES),ZZTOT(1:N_TRACKED_SPECIES),HGAS,DT_CFL
REAL(EB),INTENT(IN) :: DT
INTEGER, INTENT(IN) :: NNE
INTEGER :: NN,ND,NC,NS,NSS,ITMP,ITCOUNT,N_SUBSTEPS,N_SUBSTEPS_DUCT
LOGICAL :: CYCLE_FLAG,SUB_CYCLE_FLAG
TYPE (DUCTNODE_TYPE), POINTER :: DN=>NULL(),DN2=>NULL()
TYPE (DUCT_TYPE), POINTER :: DU=>NULL()
TYPE (NETWORK_TYPE), POINTER :: NE=>NULL()

TNOW=CURRENT_TIME()
NE => NETWORK(NNE)
DUCT%UPDATED = .FALSE.
DO NN = 1,NE%N_DUCTNODES
   DN=>DUCTNODE(NE%NODE_INDEX(NN))
   DN%UPDATED = .FALSE.
   IF (DN%VENT .OR. DN%LEAKAGE) THEN
      DU => DUCT(DN%DUCT_INDEX(1))
      IF (DU%VEL(NEW)*DN%DIR(1) < -TWO_EPSILON_EB) DN%UPDATED = .TRUE.
   ENDIF
ENDDO

DO ND = 1,NE%N_DUCTS
   DUCT(NE%DUCT_INDEX(ND))%UPDATED = .FALSE.
ENDDO

! Outputs number of substeps required to maintain CFL for mass transport
N_SUBSTEPS = 1
DO ND = 1,NE%N_DUCTS
   DU=>DUCT(NE%DUCT_INDEX(ND))
   IF (DU%N_CELLS > 1 .AND. ABS(DU%VEL(NEW)) > 0._EB) THEN
      DT_CFL = DU%LENGTH/ABS(DU%VEL(NEW)) ! relevant CFL is all mass leaving a duct within one DT
      N_SUBSTEPS_DUCT = MAX(1,CEILING(DT/DT_CFL))
      IF (N_SUBSTEPS_DUCT > N_SUBSTEPS) THEN
         N_SUBSTEPS = N_SUBSTEPS_DUCT
      ENDIF
   ENDIF
ENDDO
DT_MT = DT/REAL(N_SUBSTEPS,EB)
SUB_CYCLE_FLAG = .TRUE.

SUBSTEP_LOOP: DO NSS = 1, N_SUBSTEPS
   ITER_LOOP: DO
      IF (NSS == N_SUBSTEPS) SUB_CYCLE_FLAG = .FALSE. ! Stops the final substep call of HVAC_MASS_... and UPDATE_NODE...
      CYCLE_FLAG = .FALSE.
      DUCT_LOOP:DO ND = 1,NE%N_DUCTS
         DU=>DUCT(NE%DUCT_INDEX(ND))
         IF (DU%UPDATED) CYCLE DUCT_LOOP
         CYCLE_FLAG = .TRUE.
         IF (DU%VEL(NEW) > TWO_EPSILON_EB) THEN
            DN => DUCTNODE(DU%NODE_INDEX(1))
         ELSEIF (DU%VEL(NEW) < -TWO_EPSILON_EB) THEN
            DN => DUCTNODE(DU%NODE_INDEX(2))
         ELSE
            DU%UPDATED = .TRUE.
            DU%VEL(NEW) = 0._EB
            DU%RHO_D = 0.5_EB*(DUCTNODE(DU%NODE_INDEX(1))%RHO+DUCTNODE(DU%NODE_INDEX(2))%RHO)
            DU%TMP_D = 0.5_EB*(DUCTNODE(DU%NODE_INDEX(1))%TMP+DUCTNODE(DU%NODE_INDEX(2))%TMP)
            CYCLE DUCT_LOOP
         ENDIF
         IF (DN%UPDATED) THEN
            DU%RHO_D  = DN%RHO
            DU%TMP_D  = DN%TMP
            DU%CP_D   = DN%CP
            DU%UPDATED = .TRUE.
            DU%ZZ(:) = DN%ZZ(:)
            ZZ_GET = DU%ZZ
         ENDIF
      ENDDO DUCT_LOOP

      NODE_LOOP:DO NN = 1,NE%N_DUCTNODES
         DN=>DUCTNODE(NE%NODE_INDEX(NN))
         IF(DN%UPDATED) CYCLE NODE_LOOP
         CYCLE_FLAG = .TRUE.
         MTOT = 0._EB
         ETOT = 0._EB
         ZZTOT = 0._EB
         TGUESS = 0._EB
         CPTSUM = 0
         DO ND = 1,DN%N_DUCTS
            DU => DUCT(DN%DUCT_INDEX(ND))
            IF (DU%AREA<=TWO_EPSILON_EB) CYCLE
            IF (DU%VEL(NEW)*DN%DIR(ND) <= TWO_EPSILON_EB) CYCLE
            IF (.NOT. DU%UPDATED) CYCLE NODE_LOOP

            MASS_TRANSPORT_IF: IF (DU%N_CELLS==1) THEN
               ! Duct is not discretized
               VFLOW = ABS(DU%VEL(NEW)*DU%AREA)
               MTOT = MTOT + VFLOW * DU%RHO_D
               ETOT = ETOT + VFLOW * DU%RHO_D * DU%TMP_D * DU%CP_D
               TGUESS = TGUESS + VFLOW * DU%RHO_D * DU%TMP_D
               IF (STRATIFICATION) THEN
                  IF (DU%NODE_INDEX(1)==NE%NODE_INDEX(NN)) THEN
                     DN2=>DUCTNODE(DU%NODE_INDEX(2))
                  ELSE
                     DN2=>DUCTNODE(DU%NODE_INDEX(1))
                  ENDIF
                  ETOT = ETOT + DU%RHO_D*VFLOW*GVEC(3)*(DN%XYZ(3)-DN2%XYZ(3))
               ENDIF
               ZZTOT = ZZTOT + VFLOW * DU%RHO_D * DU%ZZ
            ELSE MASS_TRANSPORT_IF
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
                        DU_DX = (MFLOW - MSUM)/DU%RHO_C(NC)
                        ZZSUM(:) = ZZSUM(:) + DU%RHO_C(NC)*DU%ZZ_C(NC,:)*DU_DX
                        CPTSUM = CPTSUM + DU%RHO_C(NC)*DU%TMP_C(NC)*DU%CP_C(NC)*DU_DX
                        EXIT
                     ELSE
                        MSUM = MSUM + DU%RHO_C(NC)*DU_DX
                        ZZSUM(:) = ZZSUM(:) + DU%RHO_C(NC)*DU%ZZ_C(NC,:)*DU_DX
                        CPTSUM = CPTSUM + DU%RHO_C(NC)*DU%TMP_C(NC)*DU%CP_C(NC)*DU_DX
                     ENDIF
                  ENDDO
               ELSE
                  DO NC = 1,DU%N_CELLS
                     IF (MSUM + DU%RHO_C(NC)*DU_DX > MFLOW) THEN
                        DU_DX = (MFLOW - MSUM)/DU%RHO_C(NC)
                        ZZSUM(:) = ZZSUM(:) + DU%RHO_C(NC)*DU%ZZ_C(NC,:)*DU_DX
                        CPTSUM = CPTSUM + DU%RHO_C(NC)*DU%TMP_C(NC)*DU%CP_C(NC)*DU_DX
                        EXIT
                     ELSE
                        MSUM = MSUM + DU%RHO_C(NC)*DU_DX
                        ZZSUM(:) = ZZSUM(:) + DU%RHO_C(NC)*DU%ZZ_C(NC,:)*DU_DX
                        CPTSUM = CPTSUM + DU%RHO_C(NC)*DU%TMP_C(NC)*DU%CP_C(NC)*DU_DX
                     ENDIF
                  ENDDO
               ENDIF
               ! Cumulative sums of energy, mass, species passing through duct node from ducts
               ETOT = ETOT + CPTSUM * DU%AREA / DT_MT
               MTOT = MTOT + MFLOW * DU%AREA / DT_MT
               ZZTOT = ZZTOT + ZZSUM * DU%AREA / DT_MT
               TGUESS = TGUESS + MFLOW * DU%AREA / DT_MT * DU%TMP_D
            ENDIF MASS_TRANSPORT_IF
            ETOT = ETOT + DU%COIL_Q
         ENDDO

         IF (DN%FILTER_INDEX > 0) THEN
            MTOT = MTOT - SUM(DN%FILTER_LOADING(:,3))
            ZZTOT = ZZTOT - DN%FILTER_LOADING(:,3)
            ITMP = MIN(5000,NINT(DU%TMP_D))
            DO NS = 1,N_TRACKED_SPECIES
               ETOT = ETOT - CPBAR_Z(ITMP,NS)*DU%TMP_D*DN%FILTER_LOADING(NS,3)
               !*** CHECK THIS
               IF (STRATIFICATION .AND. .NOT. DN%VENT ) THEN
                  IF (DU%NODE_INDEX(1)==NE%NODE_INDEX(NN)) THEN
                     DN2=>DUCTNODE(DU%NODE_INDEX(2))
                  ELSE
                     DN2=>DUCTNODE(DU%NODE_INDEX(1))
                  ENDIF
                  ETOT = ETOT - DN%FILTER_LOADING(NS,3)*GVEC(3)*(DN%XYZ(3)-DN2%XYZ(3))
               ENDIF
            ENDDO
         ENDIF

         DN%UPDATED = .TRUE.
         IF (ABS(MTOT)<=TWO_EPSILON_EB) CYCLE NODE_LOOP
         ZZ_GET = 0._EB
         DN%ZZ(:)  = ZZTOT/MTOT
         ZZ_GET(1:N_TRACKED_SPECIES) = DN%ZZ(1:N_TRACKED_SPECIES)
         CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,DN%RSUM)
         ETOT = ETOT/ MTOT
         TGUESS = TGUESS/MTOT
         ITCOUNT = 0
         CP_LOOP: DO ! Newton method to find solution of T (and hence cpbar) from enthalpy
            ITCOUNT = ITCOUNT + 1
            CALL GET_ENTHALPY(ZZ_GET,HGAS,TGUESS)
            CALL GET_AVERAGE_SPECIFIC_HEAT(ZZ_GET,CP,TGUESS)
            IF (TGUESS>1._EB) THEN
               CALL GET_AVERAGE_SPECIFIC_HEAT(ZZ_GET,CP2,TGUESS-1._EB)
               DCPDT = CP - CP2
            ELSE
               CALL GET_AVERAGE_SPECIFIC_HEAT(ZZ_GET,CP2,TGUESS+1._EB)
               DCPDT = CP2- CP
            ENDIF
            CP = HGAS / TGUESS
            DN%TMP =TGUESS+(ETOT-HGAS)/(CP+TGUESS*DCPDT)
            IF (ABS(DN%TMP - TGUESS) < TWO_EPSILON_EB) EXIT CP_LOOP
            IF (ABS(DN%TMP - TGUESS)/DN%TMP < 0.0005_EB) EXIT CP_LOOP
            IF (ITCOUNT > 10) THEN
               DN%TMP = 0.5_EB*(DN%TMP+TGUESS)
               EXIT CP_LOOP
            ENDIF
            TGUESS = MAX(0._EB,DN%TMP)
         ENDDO CP_LOOP

         CALL GET_ENTHALPY(ZZ_GET,HGAS,DN%TMP)
         DN%CP = HGAS/DN%TMP
         DN%RHO = (DN%P+P_INF)/(DN%RSUM*DN%TMP)
      ENDDO NODE_LOOP

      IF (.NOT. CYCLE_FLAG) EXIT ITER_LOOP

   ENDDO ITER_LOOP

! If there are substeps remaining, we must advance mass transport and update node BCs to ensure we're conserving mass, energy
   IF (SUB_CYCLE_FLAG) THEN
      CALL UPDATE_HVAC_MASS_TRANSPORT(DT_MT)
      CALL UPDATE_NODE_BC
   ENDIF

ENDDO SUBSTEP_LOOP

T_USED(13)=T_USED(13)+CURRENT_TIME()-TNOW

END SUBROUTINE HVAC_UPDATE


SUBROUTINE RHSNODE(NETWORK_INDEX)
USE GLOBAL_CONSTANTS
INTEGER, INTENT(IN)::NETWORK_INDEX
INTEGER :: NN,ND, ARRAYLOC
TYPE(NETWORK_TYPE), POINTER::NE=>NULL()
TYPE(DUCT_TYPE), POINTER::DU=>NULL()
TYPE(DUCTNODE_TYPE), POINTER::DN=>NULL()

NE => NETWORK(NETWORK_INDEX)
DO NN = 1, NE%N_DUCTNODES
   DN => DUCTNODE(NE%NODE_INDEX(NN))
   IF (DN%FIXED .OR. DN%VENT .OR. DN%AMBIENT .OR. DN%LEAKAGE) CYCLE
   ARRAYLOC = NE%MATRIX_INDEX(NE%N_DUCTS+DUCTNODE_NE(NE%NODE_INDEX(NN)))
   DO ND = 1,DN%N_DUCTS
      DU => DUCT(DN%DUCT_INDEX(ND))
      IF (DU%FIXED) RHS(ARRAYLOC) = RHS(ARRAYLOC) + DN%DIR(ND)*DU%RHO_D*DU%VOLUME_FLOW
   END DO
   IF (DN%FILTER_INDEX > 0) RHS(ARRAYLOC) = RHS(ARRAYLOC) - SUM(DN%FILTER_LOADING(:,3))

ENDDO

END SUBROUTINE RHSNODE


SUBROUTINE LHSNODE(NETWORK_INDEX)
! Populates LHS matrix with nodal conservation data 
USE GLOBAL_CONSTANTS
INTEGER, INTENT(IN)::NETWORK_INDEX
INTEGER :: NN,ND, ARRAYLOC1,ARRAYLOC2
TYPE(NETWORK_TYPE), POINTER::NE=>NULL()
TYPE(DUCT_TYPE), POINTER::DU=>NULL()
TYPE(DUCTNODE_TYPE), POINTER::DN=>NULL()

NE => NETWORK(NETWORK_INDEX)
DO NN = 1, NE%N_DUCTNODES
   DN => DUCTNODE(NE%NODE_INDEX(NN))
   IF (DN%FIXED .OR. DN%VENT .OR. DN%LEAKAGE) CYCLE
   ARRAYLOC1 = NE%MATRIX_INDEX(NE%N_DUCTS+DUCTNODE_NE(NE%NODE_INDEX(NN)))
   DO ND = 1,DN%N_DUCTS
      DU => DUCT(DN%DUCT_INDEX(ND))
      IF (DU%FIXED .OR. DU%AREA <=TWO_EPSILON_EB) CYCLE
      ARRAYLOC2 = NE%MATRIX_INDEX(DUCT_NE(DN%DUCT_INDEX(ND)))
      LHS(ARRAYLOC1,ARRAYLOC2) = -DN%DIR(ND)*DU%RHO_D*DU%AREA
   END DO
ENDDO

END SUBROUTINE LHSNODE


SUBROUTINE DPSTARCALC

USE GLOBAL_CONSTANTS
INTEGER :: NN,IPZ
TYPE(DUCT_TYPE), POINTER::DU=>NULL()
TYPE(DUCTNODE_TYPE), POINTER::DN=>NULL()
TYPE(P_ZONE_TYPE), POINTER::PZ=>NULL()

DO IPZ = 1,N_ZONE
   PZ => P_ZONE(IPZ)
   IF (PZ%N_DUCTNODES==0) CYCLE
   DPSTAR(IPZ) = P_ZONE(IPZ)%DPSTAR * DT_HV * 0.5_EB
   DO NN = 1,PZ%N_DUCTNODES
      DN=>DUCTNODE(PZ%NODE_INDEX(NN))
      DU=>DUCT(DN%DUCT_INDEX(1))
      DPSTAR(IPZ) = DPSTAR(IPZ)  - DN%DIR(1) * DU%AREA * DU%VEL(OLD) * DT_HV/PSUM(IPZ,1) * 0.25_EB
      IF (DU%FIXED) DPSTAR(IPZ) = DPSTAR(IPZ)  + DN%DIR(1) * DU%AREA * DU%VEL(NEW) * DT_HV/PSUM(IPZ,1) *0.25_EB
   ENDDO
ENDDO

END SUBROUTINE DPSTARCALC


SUBROUTINE RHSDUCT(NETWORK_INDEX)

USE GLOBAL_CONSTANTS
INTEGER, INTENT(IN)::NETWORK_INDEX
INTEGER :: ND, ARRAYLOC,IPZ
REAL(EB) :: HEAD,XYZ(3)
TYPE(NETWORK_TYPE), POINTER::NE=>NULL()
TYPE(DUCT_TYPE), POINTER::DU=>NULL()
TYPE(DUCTNODE_TYPE), POINTER::DN=>NULL()

NE => NETWORK(NETWORK_INDEX)

DO ND = 1, NE%N_DUCTS
   DU => DUCT(NE%DUCT_INDEX(ND))
   IF (DU%FIXED .OR. DU%AREA < TWO_EPSILON_EB) CYCLE
   HEAD = 0._EB
   ARRAYLOC = NE%MATRIX_INDEX(DUCT_NE(NE%DUCT_INDEX(ND)))
   DN=>DUCTNODE(DU%NODE_INDEX(1))
   IF (DN%AMBIENT) THEN
      HEAD = HEAD + DN%P
   ELSEIF (DN%VENT .OR. DN%LEAKAGE) THEN
      HEAD = HEAD + DN%P
      IF (N_ZONE > 0) THEN
         IPZ = DN%ZONE_INDEX
         IF (IPZ > 0) HEAD = HEAD + DPSTAR(IPZ)
      ENDIF
   ENDIF
   XYZ = DN%XYZ
   DN=>DUCTNODE(DU%NODE_INDEX(2))
   IF (DN%AMBIENT) THEN
      HEAD = HEAD - DN%P
   ELSEIF (DN%VENT .OR. DN%LEAKAGE) THEN
      HEAD = HEAD - DN%P
      IF (N_ZONE > 0) THEN
         IPZ = DN%ZONE_INDEX
         IF (IPZ > 0) HEAD = HEAD - DPSTAR(IPZ)
      ENDIF
   ENDIF
   XYZ = DN%XYZ - XYZ
   IF(.NOT. DU%LEAKAGE) THEN
     ! IF (STRATIFICATION) THEN
    !     HEAD = HEAD + (DN%P+P_INF)*(1._EB-EXP(-(GVEC(1)*XYZ(1)+GVEC(2)*XYZ(2)+GVEC(3)*XYZ(3))/(DU%RSUM_D*DU%TMP_D)))
    !  ELSE
         HEAD = HEAD + (GVEC(1)*XYZ(1)+GVEC(2)*XYZ(2)+GVEC(3)*XYZ(3))*DU%RHO_D
    !  ENDIF
   ENDIF
   RHS(ARRAYLOC) = DU%VEL(OLD)+DT_HV/DU%LENGTH*((HEAD+SUM(DU%DP_FAN)*0.5_EB)/DU%RHO_D - &
                   0.125_EB*DU%TOTAL_LOSS*ABS(DU%VEL(GUESS)+DU%VEL(OLD))*DU%VEL(OLD))
ENDDO

END SUBROUTINE RHSDUCT


SUBROUTINE LHSDUCT(NETWORK_INDEX)
USE GLOBAL_CONSTANTS
INTEGER, INTENT(IN)::NETWORK_INDEX
INTEGER :: NN,ND, ARRAYLOC1,ARRAYLOC2
TYPE(NETWORK_TYPE), POINTER::NE=>NULL()
TYPE(DUCT_TYPE), POINTER::DU=>NULL(),DU2=>NULL()
TYPE(DUCTNODE_TYPE), POINTER::DN=>NULL(),DN2=>NULL()
TYPE(P_ZONE_TYPE), POINTER::PZ=>NULL()

NE => NETWORK(NETWORK_INDEX)

DUCT_LOOP: DO ND = 1, NE%N_DUCTS
   DU => DUCT(NE%DUCT_INDEX(ND))
   IF (DU%FIXED .OR. DU%AREA < TWO_EPSILON_EB) CYCLE DUCT_LOOP
   ARRAYLOC1 = NE%MATRIX_INDEX(DUCT_NE(NE%DUCT_INDEX(ND)))
   LHS(ARRAYLOC1,ARRAYLOC1) = 1._EB+0.125_EB*DU%TOTAL_LOSS*ABS(DU%VEL(OLD)+DU%VEL(GUESS))*DT_HV/DU%LENGTH
   DN=>DUCTNODE(DU%NODE_INDEX(1))
   IF (.NOT. DN%VENT) THEN
      IF (.NOT. DN%AMBIENT .AND. .NOT. DN%LEAKAGE) THEN
         ARRAYLOC2 = NE%MATRIX_INDEX(NE%N_DUCTS+DUCTNODE_NE(DU%NODE_INDEX(1)))
         LHS(ARRAYLOC1,ARRAYLOC2) = -DT_HV/(DU%RHO_D*DU%LENGTH)
      ENDIF
      IF (DN%LEAKAGE .AND. .NOT. DN%AMBIENT) THEN
         PZ => P_ZONE(DN%ZONE_INDEX)
         DO NN=1,PZ%N_DUCTNODES
            DN2=>DUCTNODE(PZ%NODE_INDEX(NN))
            DU2=>DUCT(DN2%DUCT_INDEX(1))
            IF (DU2%AREA < TWO_EPSILON_EB .OR. DU2%FIXED) CYCLE
            ARRAYLOC2 = NE%MATRIX_INDEX(DUCT_NE(DN2%DUCT_INDEX(1)))
            LHS(ARRAYLOC1,ARRAYLOC2) = LHS(ARRAYLOC1,ARRAYLOC2) - 0.25_EB*DN2%DIR(1)*DU2%AREA*DT_HV**2 / &
                                       (PSUM(DN%ZONE_INDEX,1)*DU%RHO_D*DU%LENGTH)
         ENDDO
      ENDIF
   ELSE
      IF (DN%ZONE_INDEX >0) THEN
         PZ => P_ZONE(DN%ZONE_INDEX)
         DO NN=1,PZ%N_DUCTNODES
            DN2=>DUCTNODE(PZ%NODE_INDEX(NN))
            DU2=>DUCT(DN2%DUCT_INDEX(1))
            IF (DU2%AREA < TWO_EPSILON_EB .OR. DU2%FIXED) CYCLE
            ARRAYLOC2 = NE%MATRIX_INDEX(DUCT_NE(DN2%DUCT_INDEX(1)))
            LHS(ARRAYLOC1,ARRAYLOC2) = LHS(ARRAYLOC1,ARRAYLOC2) - 0.25_EB*DN2%DIR(1)*DU2%AREA*DT_HV**2 / &
                                       (PSUM(DN%ZONE_INDEX,1)*DU%RHO_D*DU%LENGTH)
         ENDDO
      ENDIF
   ENDIF
   DN=>DUCTNODE(DU%NODE_INDEX(2))
   IF (.NOT. DN%VENT) THEN
      IF (.NOT. DN%AMBIENT .AND. .NOT. DN%LEAKAGE) THEN
         ARRAYLOC2 = NE%MATRIX_INDEX(NE%N_DUCTS+DUCTNODE_NE(DU%NODE_INDEX(2)))
         LHS(ARRAYLOC1,ARRAYLOC2) = DT_HV/(DU%RHO_D*DU%LENGTH)
      ENDIF
      IF (DN%LEAKAGE .AND. .NOT. DN%AMBIENT) THEN
         PZ => P_ZONE(DN%ZONE_INDEX)
         DO NN=1,PZ%N_DUCTNODES
            DN2=>DUCTNODE(PZ%NODE_INDEX(NN))
            DU2=>DUCT(DN2%DUCT_INDEX(1))
            IF (DU2%AREA < TWO_EPSILON_EB .OR. DU2%FIXED) CYCLE
            ARRAYLOC2 = NE%MATRIX_INDEX(DUCT_NE(DN2%DUCT_INDEX(1)))
            LHS(ARRAYLOC1,ARRAYLOC2) = LHS(ARRAYLOC1,ARRAYLOC2) + 0.25_EB*DN2%DIR(1)*DU2%AREA*DT_HV**2 / &
                                       (PSUM(DN%ZONE_INDEX,1)*DU%RHO_D*DU%LENGTH)
         ENDDO
      ENDIF
   ELSE
      IF (DN%ZONE_INDEX >0) THEN
         PZ => P_ZONE(DN%ZONE_INDEX)
         DO NN=1,PZ%N_DUCTNODES
            DN2=>DUCTNODE(PZ%NODE_INDEX(NN))
            DU2=>DUCT(DN2%DUCT_INDEX(1))
            IF (DU2%AREA < TWO_EPSILON_EB .OR. DU2%FIXED) CYCLE
            ARRAYLOC2 = NE%MATRIX_INDEX(DUCT_NE(DN2%DUCT_INDEX(1)))
            LHS(ARRAYLOC1,ARRAYLOC2) = LHS(ARRAYLOC1,ARRAYLOC2) + 0.25_EB*DN2%DIR(1)*DU2%AREA*DT_HV**2 / &
                                       (PSUM(DN%ZONE_INDEX,1)*DU%RHO_D*DU%LENGTH)
         ENDDO
      ENDIF
   ENDIF
ENDDO DUCT_LOOP

END SUBROUTINE LHSDUCT


SUBROUTINE RHS_SYSTEM(DUCTRUN_INDEX)
! Populates right hand side of matrix for steady state solution of ductrun to output system curve
USE GLOBAL_CONSTANTS
INTEGER, INTENT(IN)::DUCTRUN_INDEX
INTEGER :: ARRAYLOC, ND, FAN_COUNTER, IPZ
REAL(EB) :: HEAD,FAN_PRES,XYZ(3)
TYPE(DUCTRUN_TYPE), POINTER::DR=>NULL()
TYPE(DUCT_TYPE), POINTER::DU=>NULL()
TYPE(DUCTNODE_TYPE), POINTER::DN=>NULL()

DR => DUCTRUN(DUCTRUN_INDEX)
FAN_COUNTER = 1

DO ND = 1, DR%N_DUCTS
   DU => DUCT(DR%DUCT_INDEX(ND))
   IF (DU%FIXED .OR. DU%AREA < TWO_EPSILON_EB) CYCLE ! duct fixed: has a fixed volume, mass or fan
   HEAD = 0._EB
   ARRAYLOC = DR%MATRIX_SYSTEM_INDEX(DUCT_DR(DR%DUCT_INDEX(ND)))
   DN => DUCTNODE(DU%NODE_INDEX(1)) ! point to node 1, if fixed then need to add p1
   IF (DN%AMBIENT) THEN
      HEAD = DN%P
   ELSEIF(DN%VENT .OR. DN%LEAKAGE) THEN
      HEAD = DN%P
      IF (N_ZONE > 0) THEN
         IPZ = DN%ZONE_INDEX
         IF (IPZ > 0) HEAD = HEAD + DPSTAR(IPZ) ! adds zone pressure
      ENDIF
   ENDIF
   XYZ = DN%XYZ ! get location ready for delta z

   DN => DUCTNODE(DU%NODE_INDEX(2)) ! point to node 2, if fixed then need to subtract p2
   IF (DN%AMBIENT) THEN
      HEAD = HEAD - DN%P
   ELSEIF (DN%VENT .OR. DN%LEAKAGE) THEN
      HEAD = HEAD - DN%P
      IF (N_ZONE > 0) THEN
         IPZ = DN%ZONE_INDEX
         IF (IPZ > 0) HEAD = HEAD - DPSTAR(IPZ)
      ENDIF
   ENDIF

   XYZ = DN%XYZ - XYZ ! get delta z

   ! Add hydrostatic head
   IF(.NOT. DU%LEAKAGE) THEN
      HEAD = HEAD + (GVEC(1)*XYZ(1)+GVEC(2)*XYZ(2)+GVEC(3)*XYZ(3))*DU%RHO_D
   ENDIF

   ! Populate RHS_SYSTEM for each duct equation
   DR%RHS_SYSTEM_1(ARRAYLOC) = HEAD
   DR%RHS_SYSTEM_2(ARRAYLOC,:) = HEAD ! same head for all fan cases

   ! Add pressure from fan, if there are multiple fans, each needs a system curve
   FAN_PRES = 0._EB
   IF (DU%FAN_INDEX > 0) THEN
      FAN_PRES = FAN(DU%FAN_INDEX)%MAX_PRES
      DR%RHS_SYSTEM_2(ND,FAN_COUNTER) = DR%RHS_SYSTEM_2(ND,FAN_COUNTER) + FAN_PRES
      FAN_COUNTER = FAN_COUNTER + 1
   ENDIF
ENDDO

END SUBROUTINE RHS_SYSTEM


SUBROUTINE LHS_SYSTEM(DUCTRUN_INDEX)
! Populates LHS matrix for steady state solution required to output system curve
USE GLOBAL_CONSTANTS
INTEGER, INTENT(IN)::DUCTRUN_INDEX
INTEGER :: NN,ND,FN, ARRAYLOC1,ARRAYLOC2
TYPE(DUCTRUN_TYPE), POINTER::DR=>NULL()
TYPE(DUCT_TYPE), POINTER::DU=>NULL()
TYPE(DUCTNODE_TYPE), POINTER::DN=>NULL()

DR => DUCTRUN(DUCTRUN_INDEX)
! Duct energy equations
DUCT_LOOP: DO ND = 1, DR%N_DUCTS
   DU => DUCT(DR%DUCT_INDEX(ND))
   IF (DU%FIXED .OR. DU%AREA < TWO_EPSILON_EB) CYCLE DUCT_LOOP
   ARRAYLOC1 = DR%MATRIX_SYSTEM_INDEX(DUCT_DR(DR%DUCT_INDEX(ND)))
   ! Velocity coefficient
   DR%LHS_SYSTEM_1(ARRAYLOC1,ARRAYLOC1) = 0.25_EB*DU%RHO_D*DU%TOTAL_LOSS*ABS(DU%VEL_SYSTEM(1,1,1) + DU%VEL_SYSTEM(1,1,2))
   DO FN = 1, DR%N_QFANS
      DR%LHS_SYSTEM_2(ARRAYLOC1,ARRAYLOC1) = 0.25_EB*DU%RHO_D*DU%TOTAL_LOSS*ABS(DU%VEL_SYSTEM(2,FN,1) + DU%VEL_SYSTEM(2,FN,2))
   ENDDO
   ! Pressure coefficients
   ARRAYLOC2 = DR%MATRIX_SYSTEM_INDEX(DR%N_DUCTS+DUCTNODE_DR(DU%NODE_INDEX(1)))
   IF (ARRAYLOC2 > 0) THEN
      DR%LHS_SYSTEM_1(ARRAYLOC1,ARRAYLOC2) = -1._EB
      DR%LHS_SYSTEM_2(ARRAYLOC1,ARRAYLOC2) = -1._EB
   ENDIF
   ARRAYLOC2 = DR%MATRIX_SYSTEM_INDEX(DR%N_DUCTS+DUCTNODE_DR(DU%NODE_INDEX(2)))
   IF (ARRAYLOC2 > 0) THEN
      DR%LHS_SYSTEM_1(ARRAYLOC1,ARRAYLOC2) = 1._EB
      DR%LHS_SYSTEM_2(ARRAYLOC1,ARRAYLOC2) = 1._EB
   ENDIF
ENDDO DUCT_LOOP

! Node conservation
DO NN = 1, DR%N_DUCTNODES
   DN => DUCTNODE(DR%NODE_INDEX(NN))
   IF (DN%FIXED .OR. DN%VENT .OR. DN%LEAKAGE) CYCLE ! only internal nodes
   ARRAYLOC1 = DR%MATRIX_SYSTEM_INDEX(DR%N_DUCTS+DUCTNODE_DR(DR%NODE_INDEX(NN)))
   DO ND = 1,DN%N_DUCTS
      DU => DUCT(DN%DUCT_INDEX(ND))
      IF (DU%FIXED .OR. DU%AREA <=TWO_EPSILON_EB) CYCLE
      ARRAYLOC2 = DR%MATRIX_SYSTEM_INDEX(DUCT_DR(DN%DUCT_INDEX(ND)))
      DR%LHS_SYSTEM_1(ARRAYLOC1,ARRAYLOC2) = -DN%DIR(ND)*DU%RHO_D*DU%AREA
      DR%LHS_SYSTEM_2(ARRAYLOC1,ARRAYLOC2) = -DN%DIR(ND)*DU%RHO_D*DU%AREA
   END DO
ENDDO

END SUBROUTINE LHS_SYSTEM


SUBROUTINE UPDATE_FAN(T,DUCT_INDEX)
USE MATH_FUNCTIONS, ONLY : EVALUATE_RAMP
INTEGER :: FAN_ITER
INTEGER, INTENT(IN) :: DUCT_INDEX
REAL(EB), INTENT(IN) :: T
REAL(EB) :: TSI,FLOW1,FLOW2,FUNC,FLOWGUESS
TYPE(DUCT_TYPE), POINTER::DU=>NULL()
TYPE(FAN_TYPE), POINTER::FA=>NULL()
REAL(EB) :: DEL_P,VDOT

DU=> DUCT(DUCT_INDEX)

FA=> FAN(DU%FAN_INDEX)
TSI = T - DU%FAN_ON_TIME
SELECT CASE (FA%FAN_TYPE)
   CASE(1) ! Constant flow
      DU%VEL(NEW) = FA%VOL_FLOW/DU%AREA*EVALUATE_RAMP(TSI,FA%TAU,FA%SPIN_INDEX)
      IF (DU%REVERSE) DU%VEL(NEW)=-DU%VEL(NEW)
      DU%VOLUME_FLOW = DU%VEL(NEW)*DU%AREA
      RETURN
   CASE(2) ! Quadratic
      VDOT = (0.25_EB*DU%VEL(NEW)+0.25_EB*DU%VEL(PREVIOUS)+0.5_EB*DU%VEL(OLD))*DU%AREA
      IF (DU%REVERSE) VDOT = -VDOT
      VDOT = MAX(0._EB,MIN(VDOT,FA%MAX_FLOW))
      DEL_P = FA%MAX_PRES - FA%MAX_PRES*(VDOT/FA%MAX_FLOW)**2
      DEL_P = DEL_P*EVALUATE_RAMP(TSI,FA%TAU,FA%SPIN_INDEX)
   CASE(3) !Fan curve
      VDOT = 0.5*(DU%VEL(NEW)+DU%VEL(OLD))*DU%AREA
      IF (DU%REVERSE) VDOT = -VDOT
      DEL_P = EVALUATE_RAMP(VDOT,0._EB,FA%RAMP_INDEX)*EVALUATE_RAMP(TSI,FA%TAU,FA%SPIN_INDEX)
   CASE(4) ! System curve-based quadratic fan BETA
      ! Set initial bounds for bisect
      FLOW1 = 0._EB
      FLOW2 = FA%MAX_FLOW
      FAN_ITER = 0
      FAN_LOOP: DO
         FAN_ITER = FAN_ITER + 1
         FLOWGUESS = 0.5*(FLOW1+FLOW2)
         IF (ABS(FLOWGUESS - FLOW1) < TWO_EPSILON_EB) EXIT FAN_LOOP
         FUNC = FA%MAX_PRES/((DU%VEL_SYSTEM(2,DU%QFAN_N,2)-DU%VEL_SYSTEM(1,1,2))*DU%AREA)**2
         FUNC = FUNC * (FLOWGUESS-DU%VEL_SYSTEM(1,1,2)*DU%AREA)**2
         FUNC = FUNC - FA%MAX_PRES + FA%MAX_PRES*(FLOWGUESS/FA%MAX_FLOW)**2
         IF (FUNC > 0) THEN
            FLOW2 = FLOWGUESS
         ELSE
            FLOW1 = FLOWGUESS
         ENDIF
         IF (FAN_ITER > 100) THEN
            FLOWGUESS = 0.5_EB*(FLOW1 + FLOW2)
            EXIT FAN_LOOP
         ENDIF
      ENDDO FAN_LOOP
      ! Output fan pressure equating to output flow rate
      DEL_P = FA%MAX_PRES - FA%MAX_PRES*(FLOWGUESS/FA%MAX_FLOW)**2
      DEL_P = DEL_P*EVALUATE_RAMP(TSI,FA%TAU,FA%SPIN_INDEX)
END SELECT

IF (DU%REVERSE) DEL_P=-DEL_P
DU%DP_FAN(NEW) = DEL_P

END SUBROUTINE UPDATE_FAN


SUBROUTINE HVAC_BC_IN(NM)

! Average gas properties at VENTs connected to HVAC system

USE PHYSICAL_FUNCTIONS, ONLY : GET_SPECIFIC_GAS_CONSTANT,GET_ENTHALPY
INTEGER, INTENT(IN) :: NM
INTEGER :: II,JJ,KK,IW,IOR,IZ1,IZ2
REAL(EB) :: ZZ_GET(1:N_TRACKED_SPECIES),AREA,TNOW,P_AVE,HGAS
REAL(EB), POINTER, DIMENSION(:,:,:) :: RHOP,UP,VP,WP,HP
REAL(EB), POINTER, DIMENSION(:,:) :: PBARP
TYPE (MESH_TYPE),POINTER :: M=>NULL()
TYPE (SURFACE_TYPE), POINTER :: SF=>NULL()
TYPE (WALL_TYPE), POINTER :: WC=>NULL()

IF (EVACUATION_ONLY(NM)) RETURN

TNOW=CURRENT_TIME()

M => MESHES(NM)
CALL POINT_TO_MESH(NM)

NODE_X(:,NM) = 0._EB
NODE_Y(:,NM) = 0._EB
NODE_Z(:,NM) = 0._EB
NODE_ZONE(:,NM) = 0
NODE_H(:,NM) = 0._EB
NODE_P(:,NM) = 0._EB
NODE_TMP(:,NM)= 0._EB
NODE_RHO(:,NM) = 0._EB
NODE_AREA(:,NM) = 0._EB
NODE_ZZ(:,:,NM) = 0._EB
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

WALL_LOOP: DO IW = 1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
   WC=>WALL(IW)
   SF => SURFACE(WC%SURF_INDEX)
   ZONE_LEAK_IF: IF (ALL(SF%LEAK_PATH < 0)) THEN
      IF (WC%VENT_INDEX == 0) CYCLE WALL_LOOP
      IF (VENTS(WC%VENT_INDEX)%NODE_INDEX < 0) CYCLE WALL_LOOP
      WC%NODE_INDEX = VENTS(WC%VENT_INDEX)%NODE_INDEX
      IOR = WC%ONE_D%IOR
      II = WC%ONE_D%IIG
      JJ = WC%ONE_D%JJG
      KK = WC%ONE_D%KKG
      AREA = WC%ONE_D%AREA
      IF (WC%ONE_D%PRESSURE_ZONE /= NODE_ZONE(WC%NODE_INDEX,NM)) THEN
         IF (NODE_ZONE(WC%NODE_INDEX,NM) == 0) THEN
            NODE_ZONE(WC%NODE_INDEX,NM) = WC%ONE_D%PRESSURE_ZONE
         ELSE
            WRITE(MESSAGE,'(A,A)') 'ERROR: DUCTNODE must lie with a single pressure zone. Node: ',TRIM(DUCTNODE(WC%NODE_INDEX)%ID)
            CALL SHUTDOWN(MESSAGE); RETURN
         ENDIF
      ENDIF

      NODE_AREA(WC%NODE_INDEX,NM) = NODE_AREA(WC%NODE_INDEX,NM) + AREA
      NODE_RHO(WC%NODE_INDEX,NM) = NODE_RHO(WC%NODE_INDEX,NM) + AREA/RHOP(II,JJ,KK)
      ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(II,JJ,KK,1:N_TRACKED_SPECIES)
      NODE_ZZ(WC%NODE_INDEX,1:N_TRACKED_SPECIES,NM) = NODE_ZZ(WC%NODE_INDEX,1:N_TRACKED_SPECIES,NM) + &
                                                      ZZ_GET(1:N_TRACKED_SPECIES)*AREA
      CALL GET_ENTHALPY(ZZ_GET,HGAS,TMP(II,JJ,KK))
      NODE_TMP(WC%NODE_INDEX,NM) = NODE_TMP(WC%NODE_INDEX,NM) + TMP(II,JJ,KK)*AREA
      NODE_H(WC%NODE_INDEX,NM) = NODE_H(WC%NODE_INDEX,NM) + HGAS * AREA

      SELECT CASE (IOR)
         CASE (1)
            NODE_X(WC%NODE_INDEX,NM) = NODE_X(WC%NODE_INDEX,NM) + X(II-1)*AREA
            NODE_Y(WC%NODE_INDEX,NM) = NODE_Y(WC%NODE_INDEX,NM) + YC(JJ)*AREA
            NODE_Z(WC%NODE_INDEX,NM) = NODE_Z(WC%NODE_INDEX,NM) + ZC(KK)*AREA
         CASE(-1)
            NODE_X(WC%NODE_INDEX,NM) = NODE_X(WC%NODE_INDEX,NM) + X(II)*AREA
            NODE_Y(WC%NODE_INDEX,NM) = NODE_Y(WC%NODE_INDEX,NM) + YC(JJ)*AREA
            NODE_Z(WC%NODE_INDEX,NM) = NODE_Z(WC%NODE_INDEX,NM) + ZC(KK)*AREA
         CASE (2)
            NODE_X(WC%NODE_INDEX,NM) = NODE_X(WC%NODE_INDEX,NM) + XC(II)*AREA
            NODE_Y(WC%NODE_INDEX,NM) = NODE_Y(WC%NODE_INDEX,NM) + Y(JJ-1)*AREA
            NODE_Z(WC%NODE_INDEX,NM) = NODE_Z(WC%NODE_INDEX,NM) + ZC(KK)*AREA
         CASE(-2)
            NODE_X(WC%NODE_INDEX,NM) = NODE_X(WC%NODE_INDEX,NM) + XC(II)*AREA
            NODE_Y(WC%NODE_INDEX,NM) = NODE_Y(WC%NODE_INDEX,NM) + Y(JJ)*AREA
            NODE_Z(WC%NODE_INDEX,NM) = NODE_Z(WC%NODE_INDEX,NM) + ZC(KK)*AREA
         CASE (3)
            NODE_X(WC%NODE_INDEX,NM) = NODE_X(WC%NODE_INDEX,NM) + XC(II)*AREA
            NODE_Y(WC%NODE_INDEX,NM) = NODE_Y(WC%NODE_INDEX,NM) + YC(JJ)*AREA
            NODE_Z(WC%NODE_INDEX,NM) = NODE_Z(WC%NODE_INDEX,NM) + Z(KK-1)*AREA
         CASE (-3)
            NODE_X(WC%NODE_INDEX,NM) = NODE_X(WC%NODE_INDEX,NM) + XC(II)*AREA
            NODE_Y(WC%NODE_INDEX,NM) = NODE_Y(WC%NODE_INDEX,NM) + YC(JJ)*AREA
            NODE_Z(WC%NODE_INDEX,NM) = NODE_Z(WC%NODE_INDEX,NM) + Z(KK)*AREA
      END SELECT
      IF (HVAC_LOCAL_PRESSURE) THEN
         SELECT CASE (IOR)
            CASE DEFAULT
               P_AVE = PBARP(KK,WC%ONE_D%PRESSURE_ZONE)
            CASE (3)
               P_AVE = 0.5_EB*(PBARP(KK-1,WC%ONE_D%PRESSURE_ZONE)+PBARP(KK,WC%ONE_D%PRESSURE_ZONE))
            CASE (-3)
               P_AVE = 0.5_EB*(PBARP(KK,WC%ONE_D%PRESSURE_ZONE)+PBARP(KK+1,WC%ONE_D%PRESSURE_ZONE))
         END SELECT
         NODE_P(WC%NODE_INDEX,NM) = NODE_P(WC%NODE_INDEX,NM) + (P_AVE+RHO(II,JJ,KK)*(HP(II,JJ,KK)-KRES(II,JJ,KK)))*AREA
      ELSE
         SELECT CASE (IOR)
            CASE (1)
               NODE_P(WC%NODE_INDEX,NM) = NODE_P(WC%NODE_INDEX,NM) + &
                                          (PBARP(KK,WC%ONE_D%PRESSURE_ZONE)-RHO(II,JJ,KK) * &
                                          0.5_EB*(UP(II-1,JJ,KK)+WC%ONE_D%U_NORMAL)**2 * &
                                          SIGN(1._EB,UP(II-1,JJ,KK)+WC%ONE_D%U_NORMAL))* AREA
            CASE(-1)
               NODE_P(WC%NODE_INDEX,NM) = NODE_P(WC%NODE_INDEX,NM) + &
                                          (PBARP(KK,WC%ONE_D%PRESSURE_ZONE)+RHO(II,JJ,KK) * &
                                          0.5_EB*(UP(II,JJ,KK)-WC%ONE_D%U_NORMAL)**2  * &
                                          SIGN(1._EB,UP(II,JJ,KK)-WC%ONE_D%U_NORMAL))*AREA
            CASE (2)
               NODE_P(WC%NODE_INDEX,NM) = NODE_P(WC%NODE_INDEX,NM) + &
                                          (PBARP(KK,WC%ONE_D%PRESSURE_ZONE)-RHO(II,JJ,KK) * &
                                          0.5_EB*(VP(II,JJ-1,KK)+WC%ONE_D%U_NORMAL)**2 * &
                                          SIGN(1._EB,VP(II,JJ-1,KK)+WC%ONE_D%U_NORMAL))*AREA
            CASE(-2)
               NODE_P(WC%NODE_INDEX,NM) = NODE_P(WC%NODE_INDEX,NM) + &
                                          (PBARP(KK,WC%ONE_D%PRESSURE_ZONE)+RHO(II,JJ,KK) * &
                                          0.5_EB*(VP(II,JJ,KK)-WC%ONE_D%U_NORMAL)**2 * &
                                          SIGN(1._EB,VP(II,JJ,KK)-WC%ONE_D%U_NORMAL))*AREA
            CASE (3)
               NODE_P(WC%NODE_INDEX,NM) = NODE_P(WC%NODE_INDEX,NM) + &
                                          (PBARP(KK,WC%ONE_D%PRESSURE_ZONE)-RHO(II,JJ,KK) * &
                                          0.5_EB*(WP(II,JJ,KK-1)+WC%ONE_D%U_NORMAL)**2 * &
                                          SIGN(1._EB,WP(II,JJ,KK-1)+WC%ONE_D%U_NORMAL))*AREA
            CASE (-3)
               NODE_P(WC%NODE_INDEX,NM) = NODE_P(WC%NODE_INDEX,NM) + &
                                          (PBARP(KK,WC%ONE_D%PRESSURE_ZONE)+RHO(II,JJ,KK) * &
                                          0.5_EB*(WP(II,JJ,KK)-WC%ONE_D%U_NORMAL)**2 * &
                                          SIGN(1._EB,WP(II,JJ,KK)-WC%ONE_D%U_NORMAL))*AREA
         END SELECT
      ENDIF
   ELSE ZONE_LEAK_IF
      IF (ALL(SF%LEAK_PATH/=WC%ONE_D%PRESSURE_ZONE)) CYCLE WALL_LOOP
      IF (SF%LEAK_PATH(1) == WC%ONE_D%PRESSURE_ZONE) THEN
         IZ1 = SF%LEAK_PATH(1)
         IZ2 = SF%LEAK_PATH(2)
         WC%NODE_INDEX = DUCT(LEAK_PATH(MIN(IZ1,IZ2),MAX(IZ1,IZ2)))%NODE_INDEX(1)
      ELSE
         IZ1 = SF%LEAK_PATH(2)
         IZ2 = SF%LEAK_PATH(1)
         WC%NODE_INDEX = DUCT(LEAK_PATH(MIN(IZ1,IZ2),MAX(IZ1,IZ2)))%NODE_INDEX(2)
      ENDIF
      IF (NODE_ZONE(WC%NODE_INDEX,NM) == 0) NODE_ZONE(WC%NODE_INDEX,NM) = WC%ONE_D%PRESSURE_ZONE
      IOR = WC%ONE_D%IOR
      II = WC%ONE_D%IIG
      JJ = WC%ONE_D%JJG
      KK = WC%ONE_D%KKG
      AREA = WC%ONE_D%AREA

      NODE_AREA(WC%NODE_INDEX,NM) = NODE_AREA(WC%NODE_INDEX,NM) + AREA
      NODE_RHO(WC%NODE_INDEX,NM) = NODE_RHO(WC%NODE_INDEX,NM) + AREA/RHOP(II,JJ,KK)

      ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(II,JJ,KK,1:N_TRACKED_SPECIES)
      NODE_ZZ(WC%NODE_INDEX,1:N_TRACKED_SPECIES,NM) = NODE_ZZ(WC%NODE_INDEX,1:N_TRACKED_SPECIES,NM) + &
                                                      ZZ_GET(1:N_TRACKED_SPECIES)*AREA
      CALL GET_ENTHALPY(ZZ_GET,HGAS,TMP(II,JJ,KK))
      NODE_TMP(WC%NODE_INDEX,NM) = NODE_TMP(WC%NODE_INDEX,NM) + TMP(II,JJ,KK)*AREA
      NODE_H(WC%NODE_INDEX,NM) = NODE_H(WC%NODE_INDEX,NM) + HGAS * AREA

      SELECT CASE (IOR)
         CASE (3)
            P_AVE = 0.5_EB*(PBARP(KK-1,WC%ONE_D%PRESSURE_ZONE)+PBARP(KK,WC%ONE_D%PRESSURE_ZONE))
            NODE_Z(WC%NODE_INDEX,NM) = NODE_Z(WC%NODE_INDEX,NM) + Z(KK-1)*AREA
         CASE (-3)
            P_AVE = 0.5_EB*(PBARP(KK,WC%ONE_D%PRESSURE_ZONE)+PBARP(KK+1,WC%ONE_D%PRESSURE_ZONE))
            NODE_Z(WC%NODE_INDEX,NM) = NODE_Z(WC%NODE_INDEX,NM) + Z(KK)*AREA
         CASE DEFAULT
            P_AVE = PBARP(KK,WC%ONE_D%PRESSURE_ZONE)
            NODE_Z(WC%NODE_INDEX,NM) = NODE_Z(WC%NODE_INDEX,NM) + ZC(KK)*AREA
      END SELECT
      NODE_P(WC%NODE_INDEX,NM) = NODE_P(WC%NODE_INDEX,NM) + P_AVE*AREA
   ENDIF ZONE_LEAK_IF
ENDDO WALL_LOOP

T_USED(13)=T_USED(13)+CURRENT_TIME()-TNOW


END SUBROUTINE HVAC_BC_IN


SUBROUTINE DETERMINE_FIXED_ELEMENTS(T)
! Determines what ducts are "fixed"; i.e. they have been given fixed volume or mass flow (including fixed volume fan)
USE MATH_FUNCTIONS,ONLY:EVALUATE_RAMP
INTEGER:: NN,ND, COUNTER
REAL(EB), INTENT(IN):: T
LOGICAL :: CHANGE
TYPE(DUCT_TYPE), POINTER :: DU=>NULL()
TYPE(DUCTNODE_TYPE), POINTER :: DN=>NULL()

DUCTNODE%FIXED = .FALSE.
DUCT%FIXED = .FALSE.
DUCT%VOLUME_FLOW = 0._EB
CHANGE = .TRUE.

DUCT_LOOP: DO ND=1,N_DUCTS
   DU=>DUCT(ND)
   IF (DU%VOLUME_FLOW_INITIAL<1.E6_EB) THEN
      DU%VOLUME_FLOW = DU%VOLUME_FLOW_INITIAL*EVALUATE_RAMP(T,DU%TAU,DU%RAMP_INDEX)
      IF(DU%AREA > TWO_EPSILON_EB) DU%VEL(NEW) = DU%VOLUME_FLOW/DU%AREA
      DU%FIXED = .TRUE.
      CYCLE DUCT_LOOP
   ENDIF
   IF (DU%MASS_FLOW_INITIAL<1.E6_EB) THEN
      IF (DU%MASS_FLOW_INITIAL > 0._EB) THEN
         DU%VOLUME_FLOW = DU%MASS_FLOW_INITIAL/DUCTNODE(DU%NODE_INDEX(1))%RHO*EVALUATE_RAMP(T,DU%TAU,DU%RAMP_INDEX)
      ELSE
         DU%VOLUME_FLOW = DU%MASS_FLOW_INITIAL/DUCTNODE(DU%NODE_INDEX(2))%RHO*EVALUATE_RAMP(T,DU%TAU,DU%RAMP_INDEX)
      ENDIF
      IF(DU%AREA > TWO_EPSILON_EB) DU%VEL(NEW) = DU%VOLUME_FLOW/DU%AREA
      DU%FIXED = .TRUE.
      CYCLE DUCT_LOOP
   ENDIF
   IF (DU%FAN_INDEX>0) THEN
      IF (DU%DEVC_INDEX > 0) THEN
         DU%FAN_OPERATING = DEVICE(DU%DEVC_INDEX)%CURRENT_STATE
         IF (DU%FAN_OPERATING .AND. (DEVICE(DU%DEVC_INDEX)%CURRENT_STATE .NEQV. DEVICE(DU%DEVC_INDEX)%PRIOR_STATE)) &
            DU%FAN_ON_TIME = T
      ELSEIF (DU%CTRL_INDEX > 0) THEN
         DU%FAN_OPERATING = CONTROL(DU%CTRL_INDEX)%CURRENT_STATE
         IF (DU%FAN_OPERATING .AND. (CONTROL(DU%CTRL_INDEX)%CURRENT_STATE .NEQV. CONTROL(DU%CTRL_INDEX)%PRIOR_STATE)) &
            DU%FAN_ON_TIME = T
      ENDIF
      IF (DU%FAN_OPERATING .AND. FAN(DU%FAN_INDEX)%FAN_TYPE==1) THEN
         DU%FIXED=.TRUE.
         CALL UPDATE_FAN(T,ND)
      ELSEIF (.NOT. DU%FAN_OPERATING .AND. FAN(DU%FAN_INDEX)%FAN_TYPE==1) THEN
         DU%FIXED = .FALSE.
         DU%VOLUME_FLOW = 0._EB
      ENDIF
   ENDIF
ENDDO DUCT_LOOP

NODE_LOOP: DO NN = 1, N_DUCTNODES
      DN=>DUCTNODE(NN)
      IF (DN%VENT .OR. DN%AMBIENT .OR. DN%LEAKAGE) THEN
         DN%FIXED = .TRUE.
         CHANGE = .TRUE.
         CYCLE NODE_LOOP
      ENDIF
      COUNTER = 0
      DO ND = 1,DN%N_DUCTS
         DU=>DUCT(DN%DUCT_INDEX(ND))
         IF (DU%FIXED .OR. DU%AREA < TWO_EPSILON_EB) COUNTER = COUNTER + 1
      ENDDO
      IF (COUNTER==DN%N_DUCTS) THEN
         WRITE(MESSAGE,'(A,A)') 'ERROR: Cannot specify fixed flows for all branches of internal DUCTNODE:',TRIM(DN%ID)
         CALL SHUTDOWN(MESSAGE); RETURN
      ENDIF
ENDDO NODE_LOOP

END SUBROUTINE DETERMINE_FIXED_ELEMENTS


SUBROUTINE FIND_NETWORKS(CHANGEIN,T)
INTEGER:: NZ,NN,ND,DUCT_COUNTER(N_DUCTS),NODE_COUNTER(N_DUCTNODES),COUNTER,COUNTER2,ZONE_COUNTER(N_ZONE)
INTEGER, DIMENSION(:), ALLOCATABLE :: NETWORK_DCOUNTER,NETWORK_NCOUNTER,RENUMBER
LOGICAL, INTENT(INOUT) :: CHANGEIN
LOGICAL :: CHANGE
REAL(EB), INTENT(IN):: T
TYPE(DUCT_TYPE), POINTER :: DU=>NULL()

CHANGE = CHANGEIN

IF (N_ZONE > 0) ZONE_COUNTER = 0

DO ND = 1, N_DUCTS
   DU => DUCT(ND)
   IF(.NOT. DU%DAMPER) CYCLE
   IF (DU%DEVC_INDEX > 0) THEN
      IF (DEVICE(DU%DEVC_INDEX)%CURRENT_STATE .NEQV. DU%DAMPER_OPEN) THEN
            DU%DAMPER_OPEN = DEVICE(DU%DEVC_INDEX)%CURRENT_STATE
            CHANGE = .TRUE.
      ENDIF
   ELSEIF(DU%CTRL_INDEX > 0) THEN
      IF (CONTROL(DU%CTRL_INDEX)%CURRENT_STATE .NEQV. DU%DAMPER_OPEN) THEN
            DU%DAMPER_OPEN = CONTROL(DU%CTRL_INDEX)%CURRENT_STATE
            CHANGE = .TRUE.
      ENDIF
   ENDIF
   IF (DU%DAMPER_OPEN) THEN
      DU%AREA = DU%AREA_INITIAL
   ELSE
      DU%AREA = 0._EB
      DU%VEL  = 0._EB
   ENDIF
   IF (DU%LEAKAGE) THEN
      IF (NODE_AREA(DU%NODE_INDEX(1),1) < TWO_EPSILON_EB .OR. NODE_AREA(DU%NODE_INDEX(2),1) < TWO_EPSILON_EB) DU%AREA=0._EB
   ENDIF
ENDDO

IF (CHANGE) THEN
   IF (ALLOCATED(NETWORK)) DEALLOCATE(NETWORK)
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

   CHANGE = .TRUE.
   DO WHILE (CHANGE)
      CHANGE = .FALSE.
      DO ND = 1, N_DUCTS
         DU => DUCT(ND)
         IF (NODE_COUNTER(DU%NODE_INDEX(1)) /= NODE_COUNTER(DU%NODE_INDEX(2))) THEN
            CHANGE = .TRUE.
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
               CHANGE = .TRUE.
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
   CALL DETERMINE_FIXED_ELEMENTS(T)
   CALL SETUP_SOLUTION_POINTERS
ENDIF

END SUBROUTINE FIND_NETWORKS


SUBROUTINE FIND_DUCTRUNS
! Finds "duct runs"; being ductnodes and ducts directly (via HVAC components) connected to one another
INTEGER :: NN,NR,NN2,NN3,NN4,ND,DUCT_COUNTER(N_DUCTS),NODE_COUNTER(N_DUCTNODES),&
           NODE_CHECKED(N_DUCTNODES),CHKFLG,NODE_CONNECTED(N_DUCTNODES),N_QFANS
INTEGER, DIMENSION(:), ALLOCATABLE :: DUCTRUN_DCOUNTER,DUCTRUN_NCOUNTER

! Zeroing work arrays and initialising counters
DO NN = 1, N_DUCTNODES
   NODE_CHECKED(NN) = 0
   NODE_CONNECTED(NN) = 0
   NODE_COUNTER(NN) = 0
END DO
DO ND = 1, N_DUCTS
   DUCT_COUNTER(ND) = 0
END DO

NR = 1
NODE_COUNTER(1) = NR
NN = 1

! Finds all connected ducts and nodes and hence number of ductruns
L1:DO
   IF (NODE_CHECKED(NN) == 1) THEN
      NN = NN + 1
      CYCLE L1
   END IF
   NN3 = 1
   DO ND = 1, DUCTNODE(NN)%N_DUCTS
      DUCT_COUNTER(DUCTNODE(NN)%DUCT_INDEX(ND)) = NR
      NODE_COUNTER(DUCT(DUCTNODE(NN)%DUCT_INDEX(ND))%NODE_INDEX(1)) = NR
      NODE_COUNTER(DUCT(DUCTNODE(NN)%DUCT_INDEX(ND))%NODE_INDEX(2)) = NR
      NODE_CONNECTED(DUCT(DUCTNODE(NN)%DUCT_INDEX(ND))%NODE_INDEX(1)) = 1
      NODE_CONNECTED(DUCT(DUCTNODE(NN)%DUCT_INDEX(ND))%NODE_INDEX(2)) = 1
   END DO
   NODE_CHECKED(NN) = 1
   CHKFLG = 0
   DO NN2 = 1, N_DUCTNODES
      CHKFLG = CHKFLG + NODE_CHECKED(NN2)
   END DO
   IF (CHKFLG == N_DUCTNODES) EXIT L1
   L2:DO NN3 = 1, N_DUCTNODES
      IF (NODE_CHECKED(NN3) == 1) THEN
         CYCLE L2
      ELSE IF (NODE_CHECKED(NN3) == 0 .AND. NODE_CONNECTED(NN3) == 1) THEN
         NN = NN3
         EXIT L2
      ELSE IF (NN3 == N_DUCTNODES) THEN
         NR = NR + 1
         DO NN4 = 1, N_DUCTNODES
            NODE_CONNECTED(NN4) = 0
         END DO
         NN = 1
         EXIT L2
      ELSE
         CYCLE L2
      END IF
   END DO L2
END DO L1

N_DUCTRUNS = MAXVAL(NODE_COUNTER)
ALLOCATE(DUCTRUN(N_DUCTRUNS))

! Zeros and then sums up number of ducts and ductnodes for DUCTRUN array
DO NR=1, N_DUCTRUNS
   DUCTRUN(NR)%N_DUCTS=0
   DUCTRUN(NR)%N_DUCTNODES=0
ENDDO
DO ND = 1, N_DUCTS
   DUCTRUN(DUCT_COUNTER(ND))%N_DUCTS = DUCTRUN(DUCT_COUNTER(ND))%N_DUCTS + 1
ENDDO
DO NN = 1, N_DUCTNODES
   DUCTRUN(NODE_COUNTER(NN))%N_DUCTNODES = DUCTRUN(NODE_COUNTER(NN))%N_DUCTNODES + 1
ENDDO

! Allocates work arrays for duct, node and fan indexing
ALLOCATE(DUCTRUN_DCOUNTER(N_DUCTRUNS))
DUCTRUN_DCOUNTER=0
ALLOCATE(DUCTRUN_NCOUNTER(N_DUCTRUNS))
DUCTRUN_NCOUNTER=0

! Allocates and populates DUCTRUN duct and node indexes
DO NR = 1, N_DUCTRUNS
   ALLOCATE(DUCTRUN(NR)%DUCT_INDEX(DUCTRUN(NR)%N_DUCTS))
   ALLOCATE(DUCTRUN(NR)%NODE_INDEX(DUCTRUN(NR)%N_DUCTNODES))
ENDDO
DO ND = 1, N_DUCTS
   DUCTRUN_DCOUNTER(DUCT_COUNTER(ND)) = DUCTRUN_DCOUNTER(DUCT_COUNTER(ND)) + 1
   DUCTRUN(DUCT_COUNTER(ND))%DUCT_INDEX(DUCTRUN_DCOUNTER(DUCT_COUNTER(ND))) = ND
   DUCT_DR(ND) = DUCTRUN_DCOUNTER(DUCT_COUNTER(ND))
ENDDO
DO NN = 1, N_DUCTNODES
   DUCTRUN_NCOUNTER(NODE_COUNTER(NN)) = DUCTRUN_NCOUNTER(NODE_COUNTER(NN)) + 1
   DUCTRUN(NODE_COUNTER(NN))%NODE_INDEX(DUCTRUN_NCOUNTER(NODE_COUNTER(NN))) = NN
   DUCTNODE_DR(NN) = DUCTRUN_NCOUNTER(NODE_COUNTER(NN))
ENDDO

IF (QFAN_BETA_TEST) THEN
   ! Populates number of quadratic fans and indexes them within DUCT
   DO NR = 1, N_DUCTRUNS
      N_QFANS = 0
      DO ND = 1, DUCTRUN(NR)%N_DUCTS
         IF (DUCT(DUCTRUN(NR)%DUCT_INDEX(ND))%FAN_INDEX > 0) THEN
            IF (FAN(DUCT(DUCTRUN(NR)%DUCT_INDEX(ND))%FAN_INDEX)%FAN_TYPE == 4) THEN
               N_QFANS = N_QFANS + 1
               DUCT(DUCTRUN(NR)%DUCT_INDEX(ND))%QFAN_N = N_QFANS
            ENDIF
         ENDIF
      ENDDO
      DUCTRUN(NR)%N_QFANS = N_QFANS
   ENDDO
   ! Allocate and zero solution matrix indices and system curve velocities
   DO NR = 1, N_DUCTRUNS
      ALLOCATE(DUCTRUN(NR)%MATRIX_SYSTEM_INDEX(DUCTRUN(NR)%N_DUCTS+DUCTRUN(NR)%N_DUCTNODES))
      DUCTRUN(NR)%MATRIX_SYSTEM_INDEX = 0
      DO ND = 1, DUCTRUN(NR)%N_DUCTS
         ALLOCATE(DUCT(DUCTRUN(NR)%DUCT_INDEX(ND))%VEL_SYSTEM(2,DUCTRUN(NR)%N_QFANS,2)) ! vel_system(sys#,fan#,old/new)
         DUCT(DUCTRUN(NR)%DUCT_INDEX(ND))%VEL_SYSTEM(1,:,:) = 0._EB
         DUCT(DUCTRUN(NR)%DUCT_INDEX(ND))%VEL_SYSTEM(2,:,:) = TWO_EPSILON_EB
      ENDDO
   ENDDO
CALL SETUP_SOLUTION_SYSTEM_POINTERS
ENDIF

DEALLOCATE(DUCTRUN_DCOUNTER)
DEALLOCATE(DUCTRUN_NCOUNTER)

END SUBROUTINE FIND_DUCTRUNS


SUBROUTINE SETUP_SOLUTION_POINTERS
INTEGER:: NNE,NN,ND,COUNTER
TYPE(DUCT_TYPE), POINTER :: DU=>NULL()
TYPE(DUCTNODE_TYPE), POINTER :: DN=>NULL()
TYPE(NETWORK_TYPE), POINTER :: NE=>NULL()

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
      IF (DN%FIXED .OR. DN%VENT) CYCLE
      COUNTER = COUNTER + 1
      NE%MATRIX_INDEX(NE%N_DUCTS+NN)=COUNTER
   ENDDO
   NE%N_MATRIX=COUNTER
ENDDO

END SUBROUTINE SETUP_SOLUTION_POINTERS


SUBROUTINE SETUP_SOLUTION_SYSTEM_POINTERS
INTEGER:: NR,NN,ND,COUNTER
TYPE(DUCT_TYPE), POINTER :: DU=>NULL()
TYPE(DUCTNODE_TYPE), POINTER :: DN=>NULL()
TYPE(DUCTRUN_TYPE), POINTER :: DR=>NULL()

DO NR = 1,N_DUCTRUNS
   COUNTER = 0
   DR => DUCTRUN(NR)
   DO ND=1,DR%N_DUCTS
      DU=>DUCT(DR%DUCT_INDEX(ND))
      IF (DU%FIXED .OR. DU%AREA<=TWO_EPSILON_EB) CYCLE
      COUNTER = COUNTER + 1
      DR%MATRIX_SYSTEM_INDEX(ND)=COUNTER
   ENDDO
   DO NN=1,DR%N_DUCTNODES
      DN=>DUCTNODE(DR%NODE_INDEX(NN))
      IF (DN%FIXED .OR. DN%VENT) CYCLE
      COUNTER = COUNTER + 1
      DR%MATRIX_SYSTEM_INDEX(DR%N_DUCTS+NN)=COUNTER
   ENDDO
   DR%N_MATRIX_SYSTEM = COUNTER
ENDDO

END SUBROUTINE SETUP_SOLUTION_SYSTEM_POINTERS


SUBROUTINE UPDATE_LOSS(T,DT,NNE)
USE PHYSICAL_FUNCTIONS,ONLY:GET_VISCOSITY
USE MATH_FUNCTIONS,ONLY:EVALUATE_RAMP
REAL(EB) :: FRICTION_FACTOR,LOSS_SUM,ZZ_GET(1:N_TRACKED_SPECIES),VISCOSITY
INTEGER, INTENT(IN) :: NNE
REAL(EB), INTENT(IN) :: T,DT
INTEGER :: ND,ND2, NN,NUM_OUT
TYPE(DUCT_TYPE), POINTER :: DU=>NULL(),DU2=>NULL()
TYPE(DUCTNODE_TYPE), POINTER :: DN=>NULL()
TYPE(NETWORK_TYPE), POINTER :: NE=>NULL()

NE => NETWORK(NNE)

DO ND = 1, NE%N_DUCTS
   DUCT(NE%DUCT_INDEX(ND))%TOTAL_LOSS=0._EB
ENDDO
NODELOOP : DO NN=1,NE%N_DUCTNODES
  DN => DUCTNODE(NE%NODE_INDEX(NN))
  NODECLASS: IF (DN%FILTER_INDEX > 0) THEN
     CALL FILTER_UPDATE(DT,NE%NODE_INDEX(NN))
     IF (DUCT(DN%DUCT_INDEX(1))%AREA < TWO_EPSILON_EB .OR. DUCT(DN%DUCT_INDEX(2))%AREA < TWO_EPSILON_EB) CYCLE
     IF(DUCT(DN%DUCT_INDEX(1))%VEL(GUESS)*DN%DIR(1) > 0._EB) THEN
           DUCT(DN%DUCT_INDEX(1))%TOTAL_LOSS = DUCT(DN%DUCT_INDEX(1))%TOTAL_LOSS + &
              DN%FILTER_LOSS*DUCT(DN%DUCT_INDEX(1))%AREA/DUCT(DN%DUCT_INDEX(2))%AREA
     ELSE
        DUCT(DN%DUCT_INDEX(2))%TOTAL_LOSS = DUCT(DN%DUCT_INDEX(2))%TOTAL_LOSS + DN%FILTER_LOSS
     ENDIF
  ELSEIF(DN%VENT .OR. DN%AMBIENT) THEN NODECLASS
     IF (DUCT(DN%DUCT_INDEX(1))%AREA < TWO_EPSILON_EB) CYCLE
     IF(DUCT(DN%DUCT_INDEX(1))%VEL(GUESS)*DN%DIR(1) < 0._EB) THEN
        DUCT(DN%DUCT_INDEX(1))%TOTAL_LOSS = DUCT(DN%DUCT_INDEX(1))%TOTAL_LOSS + DN%LOSS_ARRAY(1,2)
     ELSE
        DUCT(DN%DUCT_INDEX(1))%TOTAL_LOSS = DUCT(DN%DUCT_INDEX(1))%TOTAL_LOSS + DN%LOSS_ARRAY(2,1)
     ENDIF
  ELSEIF(DN%FILTER_INDEX <=0 .AND. DN%N_DUCTS==2) THEN
     IF (DUCT(DN%DUCT_INDEX(1))%AREA < TWO_EPSILON_EB .OR. DUCT(DN%DUCT_INDEX(2))%AREA < TWO_EPSILON_EB) CYCLE
     IF(ABS(DUCT(DN%DUCT_INDEX(1))%VEL(GUESS)) > 1.E-6_EB) THEN
        IF(DUCT(DN%DUCT_INDEX(1))%VEL(GUESS)*DN%DIR(1) > 0._EB) THEN
           DUCT(DN%DUCT_INDEX(2))%TOTAL_LOSS = DUCT(DN%DUCT_INDEX(2))%TOTAL_LOSS + DN%LOSS_ARRAY(1,2)
        ELSE
           DUCT(DN%DUCT_INDEX(1))%TOTAL_LOSS = DUCT(DN%DUCT_INDEX(1))%TOTAL_LOSS + DN%LOSS_ARRAY(1,2)
        ENDIF
     ELSE
        LOSS_SUM = 0.5_EB*(DN%LOSS_ARRAY(1,2) + DN%LOSS_ARRAY(2,1)*DUCT(DN%DUCT_INDEX(2))%AREA/DUCT(DN%DUCT_INDEX(1))%AREA)
        DUCT(DN%DUCT_INDEX(1))%TOTAL_LOSS = DUCT(DN%DUCT_INDEX(1))%TOTAL_LOSS + &
                                            LOSS_SUM*DUCT(DN%DUCT_INDEX(1))%AREA/DUCT(DN%DUCT_INDEX(2))%AREA
        DUCT(DN%DUCT_INDEX(2))%TOTAL_LOSS = DUCT(DN%DUCT_INDEX(2))%TOTAL_LOSS + LOSS_SUM
     ENDIF
  ELSE NODECLASS
     NUM_OUT = 0
     DO ND=1,DN%N_DUCTS
        DU => DUCT(DN%DUCT_INDEX(ND))
        IF (DU%VEL(GUESS)*DN%DIR(ND) < 0._EB .AND. ABS(DU%VEL(GUESS)) > 1.E-6_EB) NUM_OUT = NUM_OUT + 1
     ENDDO
     IF (NUM_OUT==0) THEN
        DO ND=1,DN%N_DUCTS
           DU => DUCT(DN%DUCT_INDEX(ND))
           DO ND2=1,DN%N_DUCTS
              IF (ND2==ND) CYCLE
              DU%TOTAL_LOSS = DU%TOTAL_LOSS + DN%LOSS_ARRAY(ND2,ND)/DN%N_DUCTS
           ENDDO
        ENDDO
     ELSEIF (NUM_OUT==1) THEN
        DO ND=1,DN%N_DUCTS
           DU => DUCT(DN%DUCT_INDEX(ND))
           IF (DU%VEL(GUESS)*DN%DIR(ND) < 0._EB .AND. ABS(DU%VEL(GUESS)) > 1.E-6_EB) THEN
              NUM_OUT = ND
              EXIT
           ENDIF
        ENDDO
        DO ND=1,DN%N_DUCTS
           IF (ND==NUM_OUT) CYCLE
           DUCT(DN%DUCT_INDEX(ND))%TOTAL_LOSS = DUCT(DN%DUCT_INDEX(ND))%TOTAL_LOSS + &
                                                DN%LOSS_ARRAY(ND,NUM_OUT)*DU%AREA/DUCT(DN%DUCT_INDEX(NUM_OUT))%AREA
        ENDDO
     ELSEIF (NUM_OUT == DN%N_DUCTS - 1) THEN
        DO ND=1,DN%N_DUCTS
           DU => DUCT(DN%DUCT_INDEX(ND))
           IF (DU%VEL(GUESS)*DN%DIR(ND) < 0._EB .AND. ABS(DU%VEL(GUESS)) > 1.E-6_EB) THEN
              CYCLE
           ELSE
              NUM_OUT = ND
              EXIT
           ENDIF
        ENDDO
        DO ND=1,DN%N_DUCTS
           IF (ND==NUM_OUT) CYCLE
           DUCT(DN%DUCT_INDEX(ND))%TOTAL_LOSS = DUCT(DN%DUCT_INDEX(ND))%TOTAL_LOSS + DN%LOSS_ARRAY(NUM_OUT,ND)
        ENDDO
     ELSE
         LOSS_SUM = 0._EB
         DO ND=1,DN%N_DUCTS
            DU => DUCT(DN%DUCT_INDEX(ND))
            IF(DU%VEL(GUESS)*DN%DIR(ND) > 0._EB)  LOSS_SUM = LOSS_SUM + DU%VEL(GUESS)*DN%DIR(ND)*DU%AREA
         ENDDO
         DO ND=1,DN%N_DUCTS
            DU => DUCT(DN%DUCT_INDEX(ND))
            IF (DU%VEL(GUESS)*DN%DIR(ND) < 0._EB .AND. ABS(DU%VEL(GUESS)) > 1.E-6_EB) THEN
               DO ND2=1,DN%N_DUCTS
                  DU2 => DUCT(DN%DUCT_INDEX(ND2))
                  IF (DU2%VEL(GUESS)*DN%DIR(ND2) > 0._EB) DU%TOTAL_LOSS = DU%TOTAL_LOSS + &
                                                          DU2%VEL(GUESS)*DN%DIR(ND2)*DU2%AREA*DN%LOSS_ARRAY(ND2,ND)/LOSS_SUM
               ENDDO
            ENDIF
         ENDDO
     ENDIF
  ENDIF NODECLASS
ENDDO NODELOOP

DO ND = 1, NE%N_DUCTS
   DU => DUCT(NE%DUCT_INDEX(ND))
   IF (DU%ROUGHNESS > 0._EB) THEN
      ZZ_GET(1:N_TRACKED_SPECIES) = DU%ZZ(1:N_TRACKED_SPECIES)
      CALL GET_VISCOSITY(ZZ_GET,VISCOSITY,DU%TMP_D)
      FRICTION_FACTOR = COMPUTE_FRICTION_FACTOR(DU%RHO_D,VISCOSITY,ABS(DU%VEL(GUESS)),DU%DIAMETER,DU%ROUGHNESS)
   ELSE
      FRICTION_FACTOR = 0._EB
   ENDIF
   IF (DU%VEL(GUESS)>0._EB) THEN
      LOSS_SUM = DU%LOSS(1) * EVALUATE_RAMP(T,0._EB,DU%RAMP_LOSS_INDEX)
   ELSEIF (DU%VEL(GUESS)<0._EB) THEN
      LOSS_SUM = DU%LOSS(2) * EVALUATE_RAMP(T,0._EB,DU%RAMP_LOSS_INDEX)
   ELSE
      LOSS_SUM = 0.5_EB*(DU%LOSS(1)+DU%LOSS(2)) * EVALUATE_RAMP(T,0._EB,DU%RAMP_LOSS_INDEX)
   ENDIF
   DU%TOTAL_LOSS = DU%TOTAL_LOSS + DU%LENGTH/DU%DIAMETER*FRICTION_FACTOR + LOSS_SUM
   IF (DU%FAN_INDEX>0) THEN
      IF(.NOT. DU%FAN_OPERATING) DU%TOTAL_LOSS = DU%TOTAL_LOSS + FAN(DU%FAN_INDEX)%OFF_LOSS
   ENDIF
ENDDO

END SUBROUTINE UPDATE_LOSS


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


SUBROUTINE SET_GUESS(NNE,T)
INTEGER, INTENT(IN) :: NNE
REAL(EB), INTENT(IN):: T
INTEGER :: ND
TYPE(DUCT_TYPE),POINTER :: DU=>NULL()
TYPE(NETWORK_TYPE), POINTER :: NE=>NULL()

NE => NETWORK(NNE)

DO ND = 1,NE%N_DUCTS
   DU => DUCT(NE%DUCT_INDEX(ND))
   IF (DU%FAN_INDEX > 0 .AND. DU%FAN_OPERATING) THEN
      CALL UPDATE_FAN(T,NE%DUCT_INDEX(ND))
   ELSEIF (DU%FAN_INDEX > 0 .AND. .NOT. DU%FAN_OPERATING) THEN
      DU%DP_FAN = 0._EB
   ENDIF
   IF (DU%FAN_INDEX > 0 .AND. .NOT. DU%FAN_OPERATING) DU%DP_FAN = 0._EB
   IF (DU%FIXED) THEN
      DU%VEL(PREVIOUS)  = DU%VEL(GUESS)
      DU%VEL(GUESS)     = DU%VEL(NEW)
      CYCLE
   ENDIF
   DU%VEL(PREVIOUS)  = DU%VEL(GUESS)
   DU%VEL(GUESS)     = DU%VEL(NEW)
ENDDO

END SUBROUTINE SET_GUESS


SUBROUTINE SET_DONOR(NNE)
! Sets relevant boundary conditions for DUCTs and DUCTNODEs
USE MATH_FUNCTIONS, ONLY : EVALUATE_RAMP
INTEGER :: ND,NN
INTEGER, INTENT(IN) :: NNE
REAL(EB) :: RHOLAST,TMPLAST,FVAL,OMFVAL,ITERFRAC
TYPE(DUCT_TYPE), POINTER :: DU=>NULL()
TYPE(DUCTNODE_TYPE), POINTER :: DN=>NULL()
TYPE(NETWORK_TYPE), POINTER :: NE=>NULL()

NE => NETWORK(NNE)
ITERFRAC = REAL(ITER,EB)/REAL(ITER_MAX,EB)
FVAL = MIN(1._EB,MAX(0._EB,(ITERFRAC-ONTH))/ONTH,1._EB)
OMFVAL = 1._EB - FVAL

NODELOOP: DO NN=1,NE%N_DUCTNODES
   DN=>DUCTNODE(NE%NODE_INDEX(NN))
   IF(DN%VENT .OR. DN%AMBIENT .OR. DN%LEAKAGE) THEN
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
   IF (DU%FIXED .AND. DU%FAN_INDEX < 0) THEN ! fixed volume flow fan set by user
      DU%VEL(PREVIOUS) = DU%VEL(NEW)
      DU%VEL(GUESS) = DU%VEL(NEW)
   ENDIF
   IF (.NOT. DU%FIXED) THEN ! no user set volume or mass flow
      DU%VEL(PREVIOUS) = DU%VEL(GUESS)
      DU%VEL(GUESS) = DU%VEL(NEW)
   ENDIF
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


SUBROUTINE CONVERGENCE_CHECK(NNE)
INTEGER, INTENT(IN) :: NNE
INTEGER :: NN, ND
LOGICAL :: CONVERGED
REAL(EB) :: MSUM,MTOT,MFLOW,VEL
TYPE(NETWORK_TYPE), POINTER :: NE=>NULL()
TYPE(DUCT_TYPE), POINTER :: DU=>NULL()
TYPE(DUCTNODE_TYPE), POINTER :: DN=>NULL()

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
   IF (DN%FIXED) CYCLE
   MSUM = 0._EB
   MTOT = 0._EB
   DO ND=1,DN%N_DUCTS
      DU=>DUCT(DN%DUCT_INDEX(ND))
      IF (ABS(DU%VEL(NEW))<1.E-7_EB) THEN
         VEL = 0._EB
      ELSE
         VEL = DU%VEL(NEW)
      ENDIF
      MFLOW = DN%DIR(ND)*VEL*DU%RHO_D*DU%AREA
      MSUM = MSUM + MFLOW
      MTOT = MTOT + ABS(MFLOW)
   ENDDO
   IF (DN%FILTER_INDEX > 0) THEN
      MFLOW = SUM(DN%FILTER_LOADING(:,3))
      MSUM = MSUM - MFLOW
      MTOT = MTOT + MFLOW
   ENDIF
   IF(ABS(MSUM)< 1.E-6 * MTOT .OR. MTOT < TWO_EPSILON_EB) CYCLE
   CONVERGED = .FALSE.
ENDDO

IF (CONVERGED) ITER=ITER_MAX+1

END SUBROUTINE CONVERGENCE_CHECK


SUBROUTINE COLLAPSE_HVAC_BC

! Takes the MPI gathered mesh array of HVAC boundary conditions and updates the DUCTNODE boundary condition values.
USE PHYSICAL_FUNCTIONS, ONLY : GET_SPECIFIC_GAS_CONSTANT,GET_ENTHALPY,GET_AVERAGE_SPECIFIC_HEAT
INTEGER:: NN,NS,ITCOUNT,ZONE_TEST(NMESHES)
REAL(EB) :: AREA,RHO_SUM,CPBAR,H_G,TMP_SUM,TMP_NEW,ZZ_GET(1:N_TRACKED_SPECIES),CPBAR2,DCPDT,HGAS
TYPE(DUCTNODE_TYPE), POINTER :: DN=>NULL(),DN2=>NULL()

VENT_CUSTOM_AMBIENT: DO NN=1,N_DUCTNODES
   DN => DUCTNODE(NN)
   IF (.NOT. DN%LEAKAGE) DN%ZONE_INDEX = MAXVAL(NODE_ZONE(NN,:))
   ZONE_TEST=-DN%ZONE_INDEX-1
   WHERE (NODE_ZONE(NN,:) < DN%ZONE_INDEX) ZONE_TEST=NODE_ZONE(NN,:)-DN%ZONE_INDEX
   IF (MAXVAL(ZONE_TEST)>-DN%ZONE_INDEX) THEN
      WRITE(MESSAGE,'(A,A)') 'ERROR: DUCTNODE must lie with a single pressure zone. Node: ',TRIM(DUCTNODE(NN)%ID)
      CALL SHUTDOWN(MESSAGE); RETURN
   ENDIF

   INTERNAL_NODE_IF: IF (((DN%VENT .OR. DN%LEAKAGE) .AND. .NOT. DN%AMBIENT) .OR. &
                          (DN%AMBIENT .AND. SUM(NODE_AREA(NN,:)) > 0._EB)) THEN
      ZZ_GET = 0._EB
      AREA = SUM(NODE_AREA(NN,:))
      IF (AREA<=TWO_EPSILON_EB) THEN
         DUCT(DN%DUCT_INDEX(1))%AREA = 0._EB
         DUCT(DN%DUCT_INDEX(1))%VEL = 0._EB
         CYCLE VENT_CUSTOM_AMBIENT
      ELSE
         DUCT(DN%DUCT_INDEX(1))%AREA = DUCT(DN%DUCT_INDEX(1))%AREA_INITIAL
      ENDIF
      NODE_AREA_EX(NN) = AREA
      DN%XYZ(1) = SUM(NODE_X(NN,:))/AREA
      DN%XYZ(2) = SUM(NODE_Y(NN,:))/AREA
      DN%XYZ(3) = SUM(NODE_Z(NN,:))/AREA
      RHO_SUM = SUM(NODE_RHO(NN,:))

      DO NS=1,N_TRACKED_SPECIES
         DN%ZZ_V(NS) = SUM(NODE_ZZ(NN,NS,:))/AREA
      ENDDO


      ZZ_GET(1:N_TRACKED_SPECIES) = DN%ZZ_V(1:N_TRACKED_SPECIES)
      CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,DN%RSUM_V)

      DN%RHO_V = AREA/RHO_SUM

      !Initialize default values
      IF (DN%P < -1.E9_EB) THEN
         IF (STRATIFICATION) THEN
            DN%TMP  = TMPA + LAPSE_RATE*DN%XYZ(3)
            IF (ABS(LAPSE_RATE)>TWO_EPSILON_EB) THEN
               DN%P_OLD = P_INF*(DN%TMP/TMPA)**(GVEC(3)/RSUM0/LAPSE_RATE)
            ELSE
               DN%P_OLD = P_INF*EXP(GVEC(3)*(DN%XYZ(3)-GROUND_LEVEL)/(RSUM0*TMPA))
            ENDIF
            DN%RHO   =  DN%P_OLD/(TMPA*RSUM0)
         ELSE
            DN%P_OLD      = P_INF
            DN%RHO   =  DN%P/(DN%TMP*RSUM0)
         ENDIF
         CALL GET_ENTHALPY(ZZ_GET,HGAS,DN%TMP)
         DN%CP = HGAS / DN%TMP
      ENDIF

      DN%P = HVAC_PRES_RELAX*(SUM(NODE_P(NN,:))/AREA)+(1._EB-HVAC_PRES_RELAX)*DN%P_OLD
      TMP_SUM = SUM(NODE_TMP(NN,:))/AREA
      H_G = SUM(NODE_H(NN,:))/AREA
      ITCOUNT = 0

      DO
         ITCOUNT = ITCOUNT + 1
         CALL GET_ENTHALPY(ZZ_GET,HGAS,TMP_SUM)
         CALL GET_AVERAGE_SPECIFIC_HEAT(ZZ_GET,CPBAR,TMP_SUM)
         IF (TMP_SUM > 1._EB) THEN
            CALL GET_AVERAGE_SPECIFIC_HEAT(ZZ_GET,CPBAR2,TMP_SUM-1._EB)
            DCPDT = CPBAR - CPBAR2
         ELSE
            CALL GET_AVERAGE_SPECIFIC_HEAT(ZZ_GET,CPBAR2,TMP_SUM+1._EB)
            DCPDT = CPBAR2 - CPBAR
         ENDIF
         CPBAR = HGAS / TMP_SUM
         TMP_NEW = TMP_SUM + (H_G - HGAS)/(CPBAR+TMP_SUM*DCPDT)
         IF (ABS(TMP_SUM - TMP_NEW) < SPACING(TMP_NEW) .OR. ABS(TMP_SUM - TMP_NEW)/TMP_NEW < 0.0005_EB) EXIT
         IF (ITCOUNT>10) THEN
            TMP_NEW = 0.5_EB*(TMP_SUM+TMP_NEW)
            EXIT
         ENDIF
         TMP_SUM = TMP_NEW
      ENDDO
      DN%TMP_V = TMP_NEW
      CALL GET_ENTHALPY(ZZ_GET,HGAS,TMP_NEW)
      DN%CP_V = HGAS/TMP_NEW
   ENDIF INTERNAL_NODE_IF
END DO VENT_CUSTOM_AMBIENT

AMBIENT_LEAK: DO NN=1,N_DUCTNODES
   DN => DUCTNODE(NN)
    IF (DN%AMBIENT .AND. SUM(NODE_AREA(NN,:))<=TWO_EPSILON_EB) THEN
      !Initialize ambient nodes outside domain
      IF (DUCT(DN%DUCT_INDEX(1))%NODE_INDEX(1)==NN) THEN
         DN2 => DUCTNODE(DUCT(DN%DUCT_INDEX(1))%NODE_INDEX(2))
      ELSE
         DN2 => DUCTNODE(DUCT(DN%DUCT_INDEX(1))%NODE_INDEX(1))
      ENDIF
      IF (DN%LEAKAGE .OR. DN%XYZ(3) <-1.E9_EB) DN%XYZ = DN2%XYZ
      IF (DN%XYZ(3) <-1.E9_EB) CYCLE AMBIENT_LEAK
      DN%RSUM   = RSUM0
      IF (STRATIFICATION) THEN
         DN%TMP  = TMPA + LAPSE_RATE*DN%XYZ(3)
         IF (ABS(LAPSE_RATE)>TWO_EPSILON_EB) THEN
            DN%P = P_INF*(DN%TMP/TMPA)**(GVEC(3)/RSUM0/LAPSE_RATE)
         ELSE
            DN%P = P_INF*EXP(GVEC(3)*(DN%XYZ(3)-GROUND_LEVEL)/(RSUM0*TMPA))
         ENDIF
         DN%RHO   =  DN%P/(DN%TMP*RSUM0)
      ELSE
         DN%TMP  = TMPA
         DN%P      = P_INF
         DN%RHO    = RHOA
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

END SUBROUTINE COLLAPSE_HVAC_BC


SUBROUTINE SET_INIT_HVAC
INTEGER:: NN
TYPE(DUCTNODE_TYPE), POINTER :: DN=>NULL()

DO NN=1,N_DUCTNODES
   DN=>DUCTNODE(NN)
   IF (.NOT. DN%VENT .AND. .NOT. DN%LEAKAGE) CYCLE
   DN%CP = DN%CP_V
   DN%TMP = DN%TMP_V
   DN%RHO = DN%RHO_V
   DN%RSUM = DN%RSUM_V
   DN%ZZ = DN%ZZ_V
ENDDO

END SUBROUTINE SET_INIT_HVAC


SUBROUTINE SET_INIT_HVAC_MASS_TRANSPORT
! Initialises cell densities, temperatures, specific heats and species' for discretised ducts
USE PHYSICAL_FUNCTIONS, ONLY: GET_ENTHALPY
USE COMP_FUNCTIONS, ONLY: CURRENT_TIME
INTEGER :: ND,NN,NR
REAL(EB) :: DRHO,DTMP,DZZ(1:N_TRACKED_SPECIES), ZZ_GET(1:N_TRACKED_SPECIES),HGAS
TYPE(DUCT_TYPE), POINTER :: DU=>NULL()
TYPE(DUCTRUN_TYPE), POINTER :: DR=>NULL()

! Propagates duct interpolation method through the ductruns, checking for errors in input
DO NR = 1, N_DUCTRUNS
   DR => DUCTRUN(NR)
   DO ND = 1, DR%N_DUCTS
      DU => DUCT(DR%DUCT_INDEX(ND))
      IF (DU%DUCT_INTERP_TYPE_INDEX /= DUCT(DR%DUCT_INDEX(1))%DUCT_INTERP_TYPE_INDEX) THEN
         WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR: DUCT_INTERP_TYPE must be consistent through ductruns. Duct ID: ',TRIM(DU%ID)
         CALL SHUTDOWN(MESSAGE); RETURN
      END IF
      DU%DUCT_INTERP_TYPE_INDEX = DUCT(DR%DUCT_INDEX(1))%DUCT_INTERP_TYPE_INDEX
   END DO
END DO

! Initialize ductnode values and duct cell arrays for duct mass transport
DUCTRUN_LOOP:DO NR = 1, N_DUCTRUNS
   DR => DUCTRUN(NR)
   DUCT_LOOP:DO ND = 1, DR%N_DUCTS
      DU => DUCT(DR%DUCT_INDEX(ND))
      IF (DU%LEAKAGE) CYCLE
      IF (DU%N_CELLS==1) CYCLE
      SELECT CASE (DU%DUCT_INTERP_TYPE_INDEX)
         CASE (NODE1) ! duct cells and ductnodes adopt values from node 1 of lowest duct_index duct in ductrun
            DU%RHO_C(:) = DUCTNODE(DUCT(DR%DUCT_INDEX(1))%NODE_INDEX(1))%RHO_V
            DU%TMP_C(:) = DUCTNODE(DUCT(DR%DUCT_INDEX(1))%NODE_INDEX(1))%TMP
            DU%CP_C(:) = DUCTNODE(DUCT(DR%DUCT_INDEX(1))%NODE_INDEX(1))%CP
            DO NN = 1, DU%N_CELLS
               DU%ZZ_C(NN,1:N_TRACKED_SPECIES) = DUCTNODE(DUCT(DR%DUCT_INDEX(1))%NODE_INDEX(1))%ZZ
            ENDDO
            DO NN = 1, 2
               DUCTNODE(DU%NODE_INDEX(NN))%RHO = DUCTNODE(DUCT(DR%DUCT_INDEX(1))%NODE_INDEX(1))%RHO_V
               DUCTNODE(DU%NODE_INDEX(NN))%TMP = DUCTNODE(DUCT(DR%DUCT_INDEX(1))%NODE_INDEX(1))%TMP
               DUCTNODE(DU%NODE_INDEX(NN))%CP = DUCTNODE(DUCT(DR%DUCT_INDEX(1))%NODE_INDEX(1))%CP
               DUCTNODE(DU%NODE_INDEX(NN))%ZZ = DUCTNODE(DUCT(DR%DUCT_INDEX(1))%NODE_INDEX(1))%ZZ
            END DO
         CASE (NODE2) ! duct cells and ductnodes adopt values from node 2 of highest duct_index duct in ductrun
            DU%RHO_C = DUCTNODE(DUCT(DR%DUCT_INDEX(DR%N_DUCTS))%NODE_INDEX(2))%RHO_V
            DU%TMP_C = DUCTNODE(DUCT(DR%DUCT_INDEX(DR%N_DUCTS))%NODE_INDEX(2))%TMP
            DU%CP_C = DUCTNODE(DUCT(DR%DUCT_INDEX(DR%N_DUCTS))%NODE_INDEX(2))%CP
            DO NN = 1, DU%N_CELLS
               DU%ZZ_C(NN,1:N_TRACKED_SPECIES) = DUCTNODE(DUCT(DR%DUCT_INDEX(DR%N_DUCTS))%NODE_INDEX(2))%ZZ
            ENDDO
            DO NN = 1, 2
               DUCTNODE(DU%NODE_INDEX(NN))%RHO = DUCTNODE(DUCT(DR%DUCT_INDEX(DR%N_DUCTS))%NODE_INDEX(2))%RHO_V
               DUCTNODE(DU%NODE_INDEX(NN))%TMP = DUCTNODE(DUCT(DR%DUCT_INDEX(DR%N_DUCTS))%NODE_INDEX(2))%TMP
               DUCTNODE(DU%NODE_INDEX(NN))%CP = DUCTNODE(DUCT(DR%DUCT_INDEX(DR%N_DUCTS))%NODE_INDEX(2))%CP
               DUCTNODE(DU%NODE_INDEX(NN))%ZZ = DUCTNODE(DUCT(DR%DUCT_INDEX(DR%N_DUCTS))%NODE_INDEX(2))%ZZ
            END DO
         CASE (LINEAR_INTERPOLATION) ! linear interp' between node 1 and 2 (for verification case, has problem w/ complex networks)
            DRHO = (DUCTNODE(DU%NODE_INDEX(2))%RHO - DUCTNODE(DU%NODE_INDEX(1))%RHO) / DU%N_CELLS
            DTMP = (DUCTNODE(DU%NODE_INDEX(2))%TMP - DUCTNODE(DU%NODE_INDEX(1))%TMP) / DU%N_CELLS
            DZZ(1:N_TRACKED_SPECIES) = (DUCTNODE(DU%NODE_INDEX(2))%ZZ - DUCTNODE(DU%NODE_INDEX(1))%ZZ) / DU%N_CELLS
            DO NN = 1, DU%N_CELLS
               DU%RHO_C(NN) = DUCTNODE(DU%NODE_INDEX(1))%RHO + DRHO*(REAL(NN,EB) - 0.5_EB)
               DU%TMP_C(NN) = DUCTNODE(DU%NODE_INDEX(1))%TMP + DTMP*(REAL(NN,EB) - 0.5_EB)
               DU%ZZ_C(NN,1:N_TRACKED_SPECIES) = DUCTNODE(DU%NODE_INDEX(1))%ZZ + DZZ(:)*(REAL(NN,EB) - 0.5_EB)
               ZZ_GET = DU%ZZ_C(NN,1:N_TRACKED_SPECIES)
               CALL GET_ENTHALPY(ZZ_GET,HGAS,DU%TMP_C(NN))
               DU%CP_C(NN) = HGAS / DU%TMP_C(NN)
            ENDDO
         CASE DEFAULT
            WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR: DUCT_INTERP_TYPE is not correctly specified. Duct ID: ',TRIM(DU%ID)
            CALL SHUTDOWN(MESSAGE); RETURN
      END SELECT
   ENDDO DUCT_LOOP
ENDDO DUCTRUN_LOOP

END SUBROUTINE SET_INIT_HVAC_MASS_TRANSPORT


SUBROUTINE UPDATE_NODE_BC
!Takes the MPI gathered mesh array of HVAC boundary conditions and updates the DUCTNODE boundary condition values.
INTEGER:: NN, NS, ND
TYPE(DUCTNODE_TYPE), POINTER :: DN=>NULL()

DO NN=1,N_DUCTNODES
   DN => DUCTNODE(NN)
   NODE_TMP_EX(NN) = DN%TMP
   DO NS=1,N_TRACKED_SPECIES
      NODE_ZZ_EX(NN,NS) = DN%ZZ(NS)
   ENDDO
END DO

DO ND=1,N_DUCTS
   DUCT_MF(ND) = DUCT(ND)%VEL(NEW)*DUCT(ND)%AREA*DUCT(ND)%RHO_D
ENDDO

END SUBROUTINE UPDATE_NODE_BC


SUBROUTINE LEAKAGE_HVAC

USE PHYSICAL_FUNCTIONS, ONLY: GET_ENTHALPY
REAL(EB) :: ZZ_GET(1:N_TRACKED_SPECIES),HGAS
INTEGER :: I_DUCT,I_DUCTNODE,NZ1,NZ2
TYPE (DUCTNODE_TYPE), POINTER:: DN1=>NULL(),DN2=>NULL()
TYPE (DUCT_TYPE), POINTER:: DU=>NULL()

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
         DN1%MESH_INDEX = 1
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
         DN2%MESH_INDEX = 1
         DN2%N_DUCTS = 1
         DN2%RSUM = RSUM0
         DN2%TMP = TMPA
         WRITE(DN2%ID,'(A,1X,I0,1X,I0)') 'LEAK',NZ2,NZ1
      ENDIF
   ENDDO
ENDDO

END SUBROUTINE LEAKAGE_HVAC


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
   DN%FILTER_LOSS = FI%CLEAN_LOSS + EVALUATE_RAMP(TOTAL_LOADING,0._EB,FI%RAMP_INDEX)
ELSE
   DN%FILTER_LOSS = FI%CLEAN_LOSS + FI%LOADING_LOSS*TOTAL_LOADING
ENDIF

DU=>DUCT(DN%DUCT_INDEX(1))
IF (DU%VEL(GUESS) >= 0._EB .AND. DU%NODE_INDEX(2)==NODE_INDEX) THEN
  DN2 => DUCTNODE(DU%NODE_INDEX(1))
ELSEIF (DU%VEL(GUESS) <= 0._EB .AND. DU%NODE_INDEX(1)==NODE_INDEX) THEN
  DN2 => DUCTNODE(DU%NODE_INDEX(2))
ELSE
   DU=>DUCT(DN%DUCT_INDEX(2))
   IF (DU%VEL(GUESS) >= 0._EB .AND. DU%NODE_INDEX(2)==NODE_INDEX) THEN
   DN2 => DUCTNODE(DU%NODE_INDEX(1))
   ELSE
   DN2 => DUCTNODE(DU%NODE_INDEX(2))
   ENDIF
ENDIF

!Ultimately add in logic for condensible gases
DN%FILTER_LOADING(:,3) = DU%AREA*ABS(DU%VEL(GUESS))*DN2%RHO*DN2%ZZ*FI%EFFICIENCY
DN%FILTER_LOADING(:,2) = DN%FILTER_LOADING(:,1) + DN%FILTER_LOADING(:,3) * DT

END SUBROUTINE FILTER_UPDATE


SUBROUTINE COIL_UPDATE(T)
USE MATH_FUNCTIONS, ONLY : EVALUATE_RAMP
USE PHYSICAL_FUNCTIONS, ONLY : GET_AVERAGE_SPECIFIC_HEAT, GET_ENTHALPY
REAL(EB), INTENT(IN) :: T
REAL(EB) :: TMP_IN,TMP_OUT,TMP_GUESS,MDOT_DU,E_IN,MCP_C,CP1,CP2,DCPDT,TSI,ZZ_GET(1:N_TRACKED_SPECIES), HGAS
INTEGER :: ND,ITER
TYPE(DUCT_TYPE),POINTER::DU
TYPE(AIRCOIL_TYPE),POINTER:: AC

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

   AC => AIRCOIL(DU%AIRCOIL_INDEX)
   IF (AC%FIXED_Q > -1.E9_EB) THEN
      TSI = T - DU%COIL_ON_TIME
      DU%COIL_Q = AC%FIXED_Q*EVALUATE_RAMP(TSI,AC%TAU,AC%RAMP_INDEX)
   ELSE
      ITER = 0
      TMP_IN = DU%TMP_D
      TMP_GUESS = TMP_IN
      MDOT_DU = DU%RHO_D*ABS(DU%VEL(GUESS))*DU%AREA
      MCP_C =  AC%COOLANT_MASS_FLOW*AC%COOLANT_SPECIFIC_HEAT
      E_IN = MDOT_DU*TMP_IN*DU%CP_D + MCP_C*AC%COOLANT_TEMPERATURE
      ZZ_GET = DU%ZZ
      DO WHILE (ITER <= 10)
         CALL GET_ENTHALPY(ZZ_GET,HGAS,TMP_GUESS)
         CALL GET_AVERAGE_SPECIFIC_HEAT(ZZ_GET,CP1,TMP_GUESS)
         IF (TMP_GUESS > 1._EB) THEN
            CALL GET_AVERAGE_SPECIFIC_HEAT(ZZ_GET,CP2,TMP_GUESS-1._EB)
            DCPDT = CP1-CP2
         ELSE
            CALL GET_AVERAGE_SPECIFIC_HEAT(ZZ_GET,CP2,TMP_GUESS+1._EB)
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

      DU%COIL_Q = AC%COOLANT_MASS_FLOW * AC%COOLANT_SPECIFIC_HEAT*(AC%COOLANT_TEMPERATURE - TMP_OUT)*AC%EFFICIENCY
   ENDIF
END DO COIL_LOOP

END SUBROUTINE COIL_UPDATE


SUBROUTINE ADJUST_LEAKAGE_AREA
INTEGER :: ND
TYPE(DUCT_TYPE),POINTER::DU

DO ND=1,N_DUCTS
   DU => DUCT(ND)
   IF (.NOT. DU%LEAKAGE) CYCLE
   DU%AREA = DU%AREA_INITIAL * (ABS(DUCTNODE(DU%NODE_INDEX(1))%P-DUCTNODE(DU%NODE_INDEX(2))%P) / &
             P_ZONE(DUCTNODE(DU%NODE_INDEX(2))%ZONE_INDEX)%LEAK_REFERENCE_PRESSURE(DUCTNODE(DU%NODE_INDEX(1))%ZONE_INDEX)) ** &
             (P_ZONE(DUCTNODE(DU%NODE_INDEX(2))%ZONE_INDEX)%LEAK_PRESSURE_EXPONENT(DUCTNODE(DU%NODE_INDEX(1))%ZONE_INDEX)-0.5_EB)
ENDDO

END SUBROUTINE ADJUST_LEAKAGE_AREA


SUBROUTINE UPDATE_HVAC_MASS_TRANSPORT(DT)
USE PHYSICAL_FUNCTIONS,ONLY: GET_AVERAGE_SPECIFIC_HEAT, GET_ENTHALPY
REAL(EB), INTENT(IN) :: DT
INTEGER :: N_SUBSTEPS,ND,NS,NC,ITCOUNT
TYPE(DUCT_TYPE),POINTER :: DU=>NULL()
REAL(EB) :: CP,CP2,DCPDT,DT_CFL,DT_DUCT,MASS_FLUX,TGUESS,ZZ_GET(N_TRACKED_SPECIES),HGAS
REAL(EB), ALLOCATABLE, DIMENSION(:) :: CPT_C,CPT_F,RHOCPT_C
REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: RHOZZ_C,ZZ_F

DUCT_LOOP: DO ND = 1,N_DUCTS
   DU => DUCT(ND)
   IF (DU%N_CELLS == 1 ) CYCLE DUCT_LOOP

   ! Check for zero flow and zero area
   IF (ABS(DU%VEL(NEW))<=TWO_EPSILON_EB .OR. DU%AREA<=TWO_EPSILON_EB) CYCLE DUCT_LOOP

   MASS_FLUX = DU%RHO_D * DU%VEL(NEW)

   ! Set up of CFL and sub time step
   DT_CFL = DU%DX/(2*DU%VEL(NEW)) ! CFL for Godunov pure upwinding scheme
   N_SUBSTEPS = MAX(1,CEILING(DT/DT_CFL))
   DT_DUCT = DT/REAL(N_SUBSTEPS,EB)

   SUBSTEP_LOOP: DO NS = 1,N_SUBSTEPS
      ! Set upwind face indices and allocate flux arrays
      ALLOCATE(ZZ_F(0:DU%N_CELLS,N_TRACKED_SPECIES))
      ALLOCATE(CPT_F(0:DU%N_CELLS))
      ALLOCATE(CPT_C(DU%N_CELLS))
      ALLOCATE(RHOCPT_C(DU%N_CELLS))
      ALLOCATE(RHOZZ_C(DU%N_CELLS,N_TRACKED_SPECIES))

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
         TGUESS = DU%TMP_C(NC)
         ITCOUNT = 0
         CP_LOOP: DO ! Uses Newton method to iterate to find solution of TMP_C from enthalpy
            ITCOUNT = ITCOUNT + 1
            CALL GET_ENTHALPY(ZZ_GET,HGAS,TGUESS)
            CALL GET_AVERAGE_SPECIFIC_HEAT(ZZ_GET,CP,TGUESS)
            IF (TGUESS>1._EB) THEN
               CALL GET_AVERAGE_SPECIFIC_HEAT(ZZ_GET,CP2,TGUESS-1._EB)
               DCPDT = CP - CP2
            ELSE
               CALL GET_AVERAGE_SPECIFIC_HEAT(ZZ_GET,CP2,TGUESS+1._EB)
               DCPDT = CP2 - CP
            ENDIF
            CP = HGAS/TGUESS
            DU%TMP_C(NC) = TGUESS + ( CPT_C(NC) - HGAS ) / ( CP + TGUESS * DCPDT )
            IF (ABS(DU%TMP_C(NC) - TGUESS) < TWO_EPSILON_EB) EXIT CP_LOOP
            IF ((DU%TMP_C(NC) - TGUESS)/DU%TMP_C(NC) < 0.0005_EB) EXIT CP_LOOP
            IF (ITCOUNT > 10) THEN
               DU%TMP_C(NC) = 0.5_EB*(DU%TMP_C(NC)+TGUESS)
               EXIT CP_LOOP
            ENDIF
            TGUESS = DU%TMP_C(NC)
         ENDDO CP_LOOP
         CALL GET_ENTHALPY(ZZ_GET,HGAS,DU%TMP_C(NC))
         DU%CP_C(NC) = HGAS / DU%TMP_C(NC)
      ENDDO DU_UPDATE_LOOP

      DEALLOCATE(RHOZZ_C)
      DEALLOCATE(ZZ_F)
      DEALLOCATE(CPT_F)
      DEALLOCATE(CPT_C)
      DEALLOCATE(RHOCPT_C)

   ENDDO SUBSTEP_LOOP

ENDDO DUCT_LOOP


END SUBROUTINE UPDATE_HVAC_MASS_TRANSPORT

END MODULE HVAC_ROUTINES
