!> \brief Derived variable types for devices

MODULE DEVICE_VARIABLES

USE PRECISION_PARAMETERS
IMPLICIT NONE

INTEGER, ALLOCATABLE, DIMENSION(:) :: DEVC_PIPE_OPERATING

!> \brief Derived type storing inputs for the PROP namelist group 

TYPE PROPERTY_TYPE
   REAL(EB) :: DENSITY,DIAMETER,EMISSIVITY,HEAT_TRANSFER_COEFFICIENT,SPECIFIC_HEAT,RTI, &
               ACTIVATION_TEMPERATURE,ACTIVATION_OBSCURATION, &
               ALPHA_E,ALPHA_C,BETA_E,BETA_C,CHARACTERISTIC_VELOCITY,PARTICLE_VELOCITY,MASS_FLOW_RATE,FLOW_RATE,FLOW_TAU, &
               GAUGE_EMISSIVITY,GAUGE_TEMPERATURE,INITIAL_TEMPERATURE,K_FACTOR,C_FACTOR,OPERATING_PRESSURE,OFFSET,&
               SPRAY_ANGLE(2,2),P0=0._EB,PX(3)=0._EB,PXX(3,3)=0._EB
   INTEGER  :: PDPA_M=0,PDPA_N=0,N_SMOKEVIEW_PARAMETERS=0,N_SMOKEVIEW_IDS=0,N_INSERT,I_VEL=0,PARTICLES_PER_SECOND
   LOGICAL  :: PDPA_INTEGRATE=.TRUE.,PDPA_NORMALIZE=.TRUE.,HISTOGRAM_NORMALIZE=.TRUE.,HISTOGRAM=.FALSE., &
               HISTOGRAM_CUMULATIVE=.FALSE.
   REAL(EB) :: PDPA_START=0._EB,PDPA_END=1.E6_EB,PDPA_RADIUS=0.1_EB
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: TABLE_ROW, V_FACTOR
   INTEGER  :: PART_INDEX=-1,FLOW_RAMP_INDEX,SPRAY_PATTERN_INDEX,Z_INDEX=-999,Y_INDEX=-999,PRESSURE_RAMP_INDEX
   CHARACTER(LABEL_LENGTH) :: SMOKEVIEW_ID(SMOKEVIEW_OBJECTS_DIMENSION),PART_ID,ID,QUANTITY,TABLE_ID,SPEC_ID='null', &
                    SMOKEVIEW_PARAMETERS(SMOKEVIEW_OBJECTS_DIMENSION)='null'
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: SPRAY_LON_CDF,SPRAY_LON,SPRAY_LAT
   REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: SPRAY_LAT_CDF
   REAL(EB) :: HISTOGRAM_LIMITS(2) = 0._EB
   INTEGER :: HISTOGRAM_NBINS=10,FED_ACTIVITY=2
END TYPE PROPERTY_TYPE

!> \brief Derived type storing location and value informaton for DEVICE\%SUBDEVICE

TYPE SUBDEVICE_TYPE
   !> !\{
   !> Intermediate value used for computing device SPATIAL_STATISTIC or TEMPORAL_STATISTIC
   REAL(EB) :: VALUE_1=0._EB,VALUE_2=0._EB,VALUE_3=0._EB
   !> !\}
   !> !\{
   !> Subdevice point, line, or bounding box physical coordinate (m)
   REAL(EB) :: X1,X2,Y1,Y2,Z1,Z2
   !> !\}
   INTEGER :: MESH !< Subdevice mesh location
   !> !\{
   !> Subdevice point, line, or bounding box grid index
   INTEGER :: I1=-1,I2=-1,J1=-1,J2=-1,K1=-1,K2=-1
   !> !\}
   INTEGER :: N_PATH=0 !< Number of grid cells along subdevice path for TRANSMISSION or PATH OBSCURATION
   !> !\{
   !> Grid index for a grid cell along subdevice path for TRANSMISSION or PATH OBSCURATION
   INTEGER, ALLOCATABLE, DIMENSION(:) :: I_PATH,J_PATH,K_PATH
   ! !\}
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: D_PATH 
   !<Segment length in a grid cell along subdevice path for TRANSMISSION or PATH OBSCURATION
END TYPE SUBDEVICE_TYPE

!> \brief Derived type for a measurement device (DEVC)

TYPE DEVICE_TYPE
   TYPE(SUBDEVICE_TYPE), ALLOCATABLE, DIMENSION(:) :: SUBDEVICE !<Array of subdevices
   REAL(EB) :: T !< Used to track time shit for a DEVC that is part of an ASPIRATION detector.
   !> !\{
   !> Physical coordinate of a point DEVC
   REAL(EB) :: X,Y,Z
   !> !\}
   !> !\{
   !> Physical coordinates for DEVC spanning a line, plane, or volume
   REAL(EB) :: X1,X2,Y1,Y2,Z1,Z2
   !> !\}
   REAL(EB) :: INITIAL_VALUE=-1.E10_EB,INSTANT_VALUE,VALUE=0._EB,SMOOTHED_VALUE=-1.E10_EB, &
               DEPTH,TMP_L,Y_C,OBSCURATION,DELAY,ROTATION,SMOOTHING_FACTOR=0._EB,VALUE_1,VALUE_2,VALUE_3,&
               SETPOINT, T_CHANGE=1.E6_EB,BYPASS_FLOWRATE,DT,TOTAL_FLOWRATE,FLOWRATE, &
               CONVERSION_ADDEND=0._EB,CONVERSION_FACTOR=1._EB,COORD_FACTOR=1._EB, &
               TIME_INTERVAL=0._EB,QUANTITY_RANGE(2),STATISTICS_START=-1.E20_EB,STATISTICS_END=1.E20_EB,&
               OVEC(3),DFVEC(3),CELL_L=-1._EB,RMS_VALUE=0._EB,RMS_VALUE2=0._EB,&
               COV_VALUE=0._EB,AVERAGE_VALUE=0._EB,AVERAGE_VALUE2=0._EB,&
               PDPA_NUMER=0._EB,PDPA_DENUM=0._EB,Z_INT,TMP_UP,TMP_LOW,TIME_PERIOD=-1._EB, &
               VALUE_1_PREVIOUS=0._EB,VALUE_2_PREVIOUS=0._EB,VALUE_3_PREVIOUS=0._EB
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: TIME_ARRAY,TIME_MAX_VALUE,TIME_MIN_VALUE
   REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: YY_SOOT
   REAL(EB), POINTER, DIMENSION(:) :: T_E,Y_E
   INTEGER  :: N_SUBDEVICES=0,IOR,IOR_ASSUMED=0,WALL_INDEX=-1,ORDINAL,MESH,&
               I_DEPTH,N_T_E,PROP_INDEX,ORIENTATION_INDEX=0,TRIP_DIRECTION,CTRL_INDEX=-1,N_INPUTS,SURF_INDEX=-1,&
               Z_INDEX=-999,Y_INDEX=-999,MATL_INDEX=-999,PART_CLASS_INDEX=0,REAC_INDEX=0,VELO_INDEX=0,&
               NO_UPDATE_DEVC_INDEX=-1,NO_UPDATE_CTRL_INDEX=-1,CFACE_INDEX=-1,FED_ACTIVITY=2,&
               DUCT_INDEX=-1,NODE_INDEX(2)=-1,POINT=1,LINE=0,LINE_COORD_CODE=123,PIPE_INDEX=1,LP_TAG=0,ORIENTATION_NUMBER=0,&
               DUCT_CELL_INDEX=-1,LOWEST_MESH=0,N_INTERVALS=-1,N_QUANTITY=1
   INTEGER, ALLOCATABLE, DIMENSION(:) :: DEVC_INDEX,SUBDEVICE_INDEX,QUANTITY_INDEX,I,J,K
   CHARACTER(LABEL_LENGTH), ALLOCATABLE, DIMENSION(:) :: QUANTITY
   CHARACTER(LABEL_LENGTH) :: ID,PROP_ID,CTRL_ID,DEVC_ID,SPATIAL_STATISTIC='null',TEMPORAL_STATISTIC='null',&
                    SURF_ID,PART_ID='null',SPEC_ID='null',MATL_ID='null',&
                    SMOKEVIEW_BAR_LABEL,UNITS='null',XYZ_UNITS='m',DUCT_ID='null',NODE_ID(2)='null',MOVE_ID='null',&
                    D_ID='null',R_ID='null',X_ID='null',Y_ID='null',Z_ID='null',&
                    INIT_ID='null',NO_UPDATE_DEVC_ID='null',NO_UPDATE_CTRL_ID='null',REAC_ID='null'
   CHARACTER(LABEL_LENGTH) :: SMOKEVIEW_LABEL
   LOGICAL :: INITIAL_STATE,CURRENT_STATE,LATCH,PRIOR_STATE,DRY=.FALSE.,HIDE_COORDINATES=.FALSE., &
              EVACUATION=.FALSE.,RELATIVE=.FALSE.,OUTPUT=.TRUE.,ABSOLUTE_VALUE=.FALSE.,USE_PREVIOUS_VALUE=.FALSE.
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: HISTOGRAM_COUNTS
END TYPE DEVICE_TYPE

! Device arrays

INTEGER :: N_PROP !< Number of PROP inputs
INTEGER :: N_DEVC ! <Number of DEVC inputs 
INTEGER :: N_DEVC_TIME !<Number of DEVC with DV\%QUANTITY='TIME'
INTEGER :: N_DEVC_LINE !<Number of DEVC with DV\%POINTS > 0
INTEGER :: MAX_DEVC_LINE_POINTS=1
INTEGER :: N_HISTOGRAM=0
INTEGER :: MAX_HISTOGRAM_NBINS=0
CHARACTER(80), ALLOCATABLE, DIMENSION(:) :: TIME_DEVC_LABEL !<Label for each DEVC written to header of CHID_devc.csv
CHARACTER(80), ALLOCATABLE, DIMENSION(:) :: TIME_DEVC_UNITS !<Units for each DEVC written to header of CHID_devc.csv
CHARACTER(80), ALLOCATABLE, DIMENSION(:) :: LINE_DEVC_LABEL !<Label for each line DEVC written to header of CHID_line.csv
CHARACTER(80), ALLOCATABLE, DIMENSION(:) :: LINE_DEVC_UNITS !<Units for each line DEVC written to header of CHID_line.csv
CHARACTER(80), ALLOCATABLE, DIMENSION(:) :: HISTOGRAM_LABEL !<Label for each histogram written to header of CHID_hist.csv
CHARACTER(80), ALLOCATABLE, DIMENSION(:) :: HISTOGRAM_UNITS !<Units for each histogram written to header of CHID_hist.csv
CHARACTER(80), ALLOCATABLE, DIMENSION(:,:) :: LINE_DEVC_VALUE !Array of line DEVC values written to CHID_line.csv
CHARACTER(80), ALLOCATABLE, DIMENSION(:,:) :: HISTOGRAM_VALUE !Array of histogram values written to CHID_hist.csv
REAL(EB), ALLOCATABLE, DIMENSION(:) :: TIME_DEVC_VALUE !<Array of DEVC values written to CHID_devc.csv
TYPE (PROPERTY_TYPE),  DIMENSION(:), ALLOCATABLE, TARGET :: PROPERTY !< Array of all defined PROP inputs
TYPE (DEVICE_TYPE),  DIMENSION(:), ALLOCATABLE, TARGET :: DEVICE !< Array of all defined DEVC inputs

END MODULE DEVICE_VARIABLES


!> \brief Derived types for evaluating control functions

MODULE CONTROL_VARIABLES

USE PRECISION_PARAMETERS

IMPLICIT NONE

!> !\{
!> Parameter defining the type of control function for CONTROL\%CONTROL_INDEX
INTEGER, PARAMETER :: AND_GATE=1, OR_GATE=2, XOR_GATE=3, X_OF_N_GATE=4, TIME_DELAY=5, DEADBAND=6, CYCLING=7, &
                      CUSTOM=8,KILL=9,CORE_DUMP=10,CF_SUM=11,CF_SUBTRACT=12,CF_MULTIPLY=13,CF_DIVIDE=14,CF_POWER=15,CF_PID=16,&
                      CF_EXP=17,CF_LOG=18,CF_SIN=19,CF_COS=20,CF_ASIN=21,CF_ACOS=22,CF_MIN=23,CF_MAX=24, CF_ABS=25
!> !\}
!> !\{
!> Parameter used to define the type of input for CONTROL\%INPUT_TYPE
INTEGER, PARAMETER :: DEVICE_INPUT=1,CONTROL_INPUT=2,CONSTANT_INPUT=3
!> !\} 

INTEGER :: N_CTRL = 0 !< Length of CONTROL
INTEGER :: N_CTRL_FILES = 0 !< Number of CHID_ctrl.csv output files

!> \brief Derived type for a control function (CTRL)

TYPE CONTROL_TYPE
   LOGICAL :: INITIAL_STATE=.FALSE. !< Initial logical state of the control function assigned at the start of the calculation
   LOGICAL :: CURRENT_STATE=.FALSE. !< Current timestep logical state of the control function
   LOGICAL :: PRIOR_STATE=.FALSE.   !< Prior timestep logical state of the control function
   LOGICAL :: LATCH=.TRUE.          !< Control function can only change state once
   LOGICAL :: UPDATED=.FALSE.       !< Control function has been updatead in the current timestep
   LOGICAL :: EVACUATION=.FALSE.    !< Control function is used by EVAC
   INTEGER :: CONTROL_INDEX=0       !< Indicates the type of control function
   INTEGER :: N_INPUTS=0            !< Number of inputs to the control function
   INTEGER :: RAMP_INDEX=0          !< Index of a RAMP used for the CUSTOM control function
   INTEGER :: MESH=1                !< Mesh location of the control function
   INTEGER :: N=1                   
   !< Number of inputs that must be true to change the INITIAL_STATE of an AT LEAST or ONLY functions
   INTEGER :: ON_BOUND=0            
   !< Negative indicates a DEADBAND INITIAL_STATE trips at the lower bound. Positive trips at the upper bound.
   INTEGER :: TRIP_DIRECTION        
   !< Negative indicates INITIAL_STATE trips when dropping below the SETPOINT. Positive is increase above SETPOINT.
   INTEGER, ALLOCATABLE, DIMENSION (:) :: INPUT
   !< Array of indices containing the device and control function inputs
   INTEGER, ALLOCATABLE, DIMENSION (:) :: INPUT_TYPE
   !< Array of inidicating if a specific input to a control function is a device and a control functon
   REAL(EB) :: SETPOINT(2)=1.E30_EB
   !<Setpoint for a control function. For a DEADBAND function contains the lower and upper bounds of the DEADBAND
   REAL(EB) :: DELAY=0._EB !<Delay time (s) for a TIME_DELAY function
   REAL(EB) :: T_CHANGE=1000000._EB !<Time the control function changed state
   REAL(EB) :: CONSTANT=-9.E30_EB !<Value assigned to CONSTANT on a CTRL input
   REAL(EB) :: INSTANT_VALUE !<Current mathematical value of the control function
   REAL(EB) :: PROPORTIONAL_GAIN !<Proportional gain for a PID control function
   REAL(EB) :: INTEGRAL_GAIN !<Intregal gain for a PID control function
   REAL(EB) :: DIFFERENTIAL_GAIN !<Differential gain for a PID control function
   REAL(EB) :: PREVIOUS_VALUE=-9.E30_EB !<Prior timestep value used in computing the time derivative term for a PID control function
   REAL(EB) :: INTEGRAL=0._EB !<Current value of the time integral term for a PID control function
   REAL(EB) :: TARGET_VALUE = -9.E30_EB !<Desired target used to compute error value for input to a PID control function
   CHARACTER(LABEL_LENGTH) :: ID='null' !<Name of control function
   CHARACTER(LABEL_LENGTH) :: RAMP_ID='null' !<Name of RAMP used for a CUSTOM control function
   CHARACTER(LABEL_LENGTH) :: INPUT_ID(40)='null' !<Array of DEVC\%ID or CTRL\%ID 
END TYPE CONTROL_TYPE

TYPE (CONTROL_TYPE),  DIMENSION(:), ALLOCATABLE, TARGET :: CONTROL !< Array of all defined CTRL inputs

END MODULE CONTROL_VARIABLES
