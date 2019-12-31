!
! This module contains the FDS+Evac human movement algorithm and
! related subprograms.  Some of the types are defined in the
! type.f90 module file.  The statistic (cumulative distributions) are in
! the ieva.f90, which includes a module dcdflib.f by netlib (http://www.netlib.org/).
! The DCDFLIB has been written in the public domain.  Some of the DCDFLIB code has
! appeared in an ACM publication and it is subject to their algorithms policy,
! see the comments at the start of the DCDFLIB in ieva.f90.
!
! Author: Timo Korhonen, VTT Technical Research Centre of Finland, 2007-2017
!
!!!!!!!!!!!!!!
!
MODULE EVAC
  !
  USE PRECISION_PARAMETERS
  USE MESH_VARIABLES
  USE GLOBAL_CONSTANTS
  USE TRAN
  USE TYPES
  USE STAT
  USE COMP_FUNCTIONS
  USE MATH_FUNCTIONS
  USE MEMORY_FUNCTIONS
  USE MESH_POINTERS
  USE PHYSICAL_FUNCTIONS, ONLY : GET_MASS_FRACTION, FED
  USE DCDFLIB, ONLY : DCDFLIB_Gamma => Gamma_ieva
  USE DEVICE_VARIABLES
  USE CONTROL_VARIABLES

  IMPLICIT NONE
  PRIVATE

  ! Public subprograms (called from the main program or read or dump)

  PUBLIC EVACUATE_HUMANS, INITIALIZE_EVACUATION, INIT_EVAC_GROUPS
  PUBLIC READ_EVAC, DUMP_EVAC, DUMP_EVAC_CSV, PREPARE_TO_EVACUATE, CLEAN_AFTER_EVACUATE
  PUBLIC EVAC_MESH_EXCHANGE, INITIALIZE_EVAC_DUMPS
  ! Public variables (needed in the main program or dump):
  PUBLIC N_DOORS, N_EXITS, N_ENTRYS, N_SSTANDS, EVAC_DOORS, EVAC_EXITS, EVAC_ENTRYS, EVAC_SSTANDS, &
       EVAC_EXIT_TYPE, EVAC_DOOR_TYPE, EVAC_ENTR_TYPE, EVAC_SSTAND_TYPE, NPC_EVAC, N_HOLES, &
       EVACUATION_TYPE, EVAC_HOLE_TYPE, EVAC_EVACS, EVAC_HOLES, N_CO_EXITS, N_DOOR_MESHES, N_STRS
  PUBLIC EVAC_EMESH_EXITS_TYPE, EMESH_EXITS, EMESH_ID, EMESH_IJK, EMESH_XB, EMESH_NM, EMESH_NFIELDS, &
       EMESH_INDEX, HUMAN_SMOKE_HEIGHT, EVAC_DELTA_SEE, EVAC_EMESH_STAIRS_TYPE, EMESH_STAIRS
  PUBLIC NO_EVAC_MESHES, INPUT_EVAC_GRIDS
  !
  CHARACTER(255):: EVAC_VERSION = '2.5.2'

  INTEGER :: INPUT_EVAC_GRIDS
  LOGICAL :: NO_EVAC_MESHES
  !
  ! This is a group of persons, who are initialized together, i.e., they have same mass, speed,
  ! etc distributions and they are all put in the given rectangle.
  ! (&EVAC lines)
  TYPE EVACUATION_TYPE
     REAL(EB) :: X1=0._EB,X2=0._EB,Y1=0._EB,Y2=0._EB,Z1=0._EB,Z2=0._EB,T_START=0._EB, Angle=0._EB
     REAL(EB) :: Tpre_mean=0._EB, Tpre_para=0._EB, Tpre_para2=0._EB, Tpre_low=0._EB, Tpre_high=0._EB
     REAL(EB) :: Tdet_mean=0._EB, Tdet_para=0._EB, Tdet_para2=0._EB, Tdet_low=0._EB, Tdet_high=0._EB
     REAL(EB) :: TIME_FALL_DOWN=0._EB, TARGET_X=0._EB, TARGET_Y=0._EB, DELTA_X=0._EB, DELTA_Y=0._EB
     REAL(EB) :: T_START_FED=0._EB
     CHARACTER(60) :: CLASS_NAME='null', ID='null', AVATAR_TYPE_NAME='null'
     CHARACTER(LABEL_LENGTH) :: GRID_NAME='null', PROP_ID='null', CROWBAR_INPUT_FILE='null'
     LOGICAL :: EVACFILE=.FALSE., After_Tpre=.FALSE., No_Persons=.FALSE., SHOW=.TRUE.
     INTEGER :: N_INITIAL=0,SAMPLING=0, IPC=0, IMESH=0, AVATAR_TYPE_INDEX=0
     INTEGER :: I_PRE_DIST=0, I_DET_DIST=0, GUARD_MEN=0
     INTEGER :: GN_MIN=0, GN_MAX=0
     INTEGER :: N_VENT_FFIELDS=0, Avatar_Color_Index=0, I_AGENT_TYPE=2, I_EVAC_THIS_MESH=1
     INTEGER, DIMENSION(3) :: RGB=-1, AVATAR_RGB=-1
     INTEGER, POINTER, DIMENSION(:) :: I_DOOR_NODES =>NULL()
     INTEGER, POINTER, DIMENSION(:) :: I_VENT_FFIELDS =>NULL()
     REAL(EB), POINTER, DIMENSION(:) :: P_VENT_FFIELDS =>NULL()
  END TYPE EVACUATION_TYPE
  !
  ! An evacuation hole, i.e., a rectangle where humans should not be put.
  ! This makes the &EVAC lines easier to define.
  ! (&EVHO lines)
  TYPE EVAC_HOLE_TYPE
     REAL(EB) :: X1=0._EB,X2=0._EB,Y1=0._EB,Y2=0._EB,Z1=0._EB,Z2=0._EB, TIME_FALL_DOWN=-1.0_EB
     CHARACTER(60) :: ID='null', PERS_ID='null', EVAC_ID='null'
     CHARACTER(LABEL_LENGTH) :: GRID_NAME='null'
     INTEGER, DIMENSION(3) :: RGB=-1
     LOGICAL :: SHOW=.TRUE.
     INTEGER :: IMESH=0, I_EVHO_THIS_MESH=1
  END TYPE EVAC_HOLE_TYPE
  !
  ! A spectator stand. IOR: which x,y line is the bottom line of the stand.
  ! ior=+1 x=x2, ior=-1 x=x1, ior=+2 y=y2, ior=-2 y=y1
  ! H is the height of the stand, S is the length along the incline.
  ! (&EVSS lines)
  TYPE EVAC_SSTAND_TYPE
     REAL(EB) :: X1=0._EB,X2=0._EB,Y1=0._EB,Y2=0._EB,Z1=0._EB,Z2=0._EB, H=0._EB, H0=0._EB, S=0._EB
     REAL(EB) :: Esc_SpeedUp=0._EB, Esc_SpeedDn=0._EB
     REAL(EB) :: FAC_V0_UP=1._EB, FAC_V0_DOWN=1._EB, FAC_V0_HORI=1._EB
     REAL(EB) :: cos_x=1._EB, cos_y=1._EB, sin_x=0._EB, sin_y=0._EB
     CHARACTER(60) :: ID='null'
     CHARACTER(LABEL_LENGTH) :: GRID_NAME='null'
     CHARACTER(LABEL_LENGTH) :: VENT_FFIELD='null'
     INTEGER, DIMENSION(3) :: RGB=-1
     INTEGER :: IMESH=0, IOR=0, I_VENT_FFIELD=0
     REAL(EB) :: UBAR0=0._EB, VBAR0=0._EB
     LOGICAL :: Use_v0=.FALSE., SHOW=.TRUE., COUNT_ONLY=.FALSE.
     REAL(EB), DIMENSION(3) :: ORIENTATION=0.0_EB
  END TYPE EVAC_SSTAND_TYPE
  !
  ! Humans belong to some small group (1 to about 5 persons).  This type
  ! collects the common properties of the group.
  TYPE GROUP_TYPE
     REAL(EB) :: GROUP_X=0._EB, GROUP_Y=0._EB, MAX_DIST_CENTER=0._EB, LIMIT_COMP=0._EB
     REAL(EB) :: GROUP_EFF=0._EB, RADIUS_COMPLETE_0=0._EB, RADIUS_COMPLETE_1=0._EB
     REAL(EB) :: Speed=0._EB, IntDose=0._EB, Tpre=0._EB, Tdoor=0._EB, Tdet=0._EB
     INTEGER :: GROUP_SIZE=0, GROUP_ID=0, COMPLETE=0, IEL=0, Avatar_Color_Index=0
     INTEGER, DIMENSION(3) :: AVATAR_RGB=-1
     INTEGER, POINTER, DIMENSION(:) :: GROUP_I_FFIELDS =>NULL()
  END TYPE GROUP_TYPE

  TYPE KNOWN_DOOR_TYPE
     INTEGER :: N_nodes=0, I_Target=0
     INTEGER, POINTER, DIMENSION(:) :: I_nodes =>NULL()
  END TYPE KNOWN_DOOR_TYPE
  !
  ! This defines a class of persons, e.g. soccer fan.
  ! (&PERS lines)
  TYPE EVAC_PERS_TYPE
     REAL(EB) :: D_mean=0._EB, D_para=0._EB, D_para2=0._EB, D_low=0._EB, D_high=0._EB
     REAL(EB) :: V_mean=0._EB, V_para=0._EB, V_para2=0._EB, V_low=0._EB, V_high=0._EB
     REAL(EB) :: Tau_mean=0._EB, Tau_para=0._EB, Tau_para2=0._EB, Tau_low=0._EB, Tau_high=0._EB
     REAL(EB) :: Tpre_mean=0._EB, Tpre_para=0._EB, Tpre_para2=0._EB, Tpre_low=0._EB, Tpre_high=0._EB
     REAL(EB) :: Tdet_mean=0._EB, Tdet_para=0._EB, Tdet_para2=0._EB, Tdet_low=0._EB, Tdet_high=0._EB
     REAL(EB) :: A=0._EB,B=0._EB,Lambda=0._EB,C_Young=0._EB,Gamma=0._EB,Kappa=0._EB, m_agent=0._EB
     REAL(EB) :: r_torso=0._EB,r_shoulder=0._EB,d_shoulder=0._EB,m_iner=0._EB, Tau_iner=0._EB
     REAL(EB) :: FAC_V0_UP=-1._EB, FAC_V0_DOWN=-1._EB, FAC_V0_HORI=-1._EB, MAXIMUM_V0_FACTOR=-1.0_EB
     !Issue1547: Added MAXIMUM_V0_FACTOR to person class type (structured type).
     CHARACTER(60) :: ID='null'
     INTEGER :: I_DIA_DIST=0, I_VEL_DIST=0, I_PRE_DIST=0, I_DET_DIST=0, I_TAU_DIST=0
     INTEGER :: Avatar_Color_Index=0
     INTEGER, DIMENSION(3) :: RGB=-1, AVATAR_RGB=-1
  END TYPE EVAC_PERS_TYPE
  !
  ! Exit door type: this just count the number of persons
  ! T_first: first person's exit time (saved for output)
  ! CHECK_FLOW: If true then the flow can not exceed Flow_max
  ! (&EXIT lines)
  TYPE EVAC_EXIT_TYPE
     REAL(EB) :: T_first=0._EB, T_last=0._EB, Flow_max=0._EB, Width=0._EB, Height=2.0_EB
     REAL(EB) :: X1=0._EB, X2=0._EB, Y1=0._EB, Y2=0._EB, Z1=0._EB, Z2=0._EB, &
          X=0._EB, Y=0._EB, Z=0._EB, Xsmoke=0._EB, Ysmoke=0._EB, Zsmoke=0._EB, &
          TIME_OPEN=0._EB, TIME_CLOSE=0._EB, R_NTARGET=0._EB, XYZ_Z=0._EB
     INTEGER :: IOR=0, ICOUNT=0, IMESH=0, INODE=0, IMODE=1
     INTEGER, DIMENSION(50) :: NTARGET=0
     REAL(EB) :: FED_CO_CO2_O2=0._EB, SOOT_DENS=0._EB, TMP_G=0._EB, RADFLUX=0._EB
     INTEGER :: II=0, JJ=0, KK=0, FED_MESH=0
     LOGICAL :: CHECK_FLOW=.FALSE., COUNT_ONLY=.FALSE., SHOW=.TRUE., COUNT_DENSITY=.FALSE., KNOWN_DOOR=.FALSE.
     LOGICAL :: LOCKED_WHEN_CLOSED=.FALSE., TARGET_WHEN_CLOSED=.FALSE., WAIT_AT_XYZ=.TRUE.
     INTEGER :: STR_INDX=0, STR_SUB_INDX=0
     CHARACTER(60) :: ID='null', PERS_ID='null', EVAC_ID='null'
     CHARACTER(60) :: TO_NODE='null'
     CHARACTER(LABEL_LENGTH) :: GRID_NAME='null'
     CHARACTER(LABEL_LENGTH) :: VENT_FFIELD='null'
     INTEGER :: I_VENT_FFIELD=0, Avatar_Color_Index=0, I_EMESH_EXITS=0
     INTEGER, DIMENSION(3) :: RGB=-1
     REAL(EB), DIMENSION(3) :: ORIENTATION=0.0_EB
  END TYPE EVAC_EXIT_TYPE
  !
  ! Like exit, but door will always put the persons to some
  ! other node. (Thus no count_only option.)
  ! (&DOOR lines)
  TYPE EVAC_DOOR_TYPE
     REAL(EB) :: T_first=0._EB, T_last=0._EB, Flow_max=0._EB, Width=0._EB, Height=2.0_EB
     REAL(EB) :: X1=0._EB, X2=0._EB, Y1=0._EB, Y2=0._EB, Z1=0._EB, Z2=0._EB, &
          X=0._EB, Y=0._EB, Z=0._EB, Xsmoke=0._EB, Ysmoke=0._EB, Zsmoke=0._EB, &
          TIME_OPEN=0._EB, TIME_CLOSE=0._EB, R_NTARGET=0._EB, XYZ_Z=0._EB
     INTEGER :: IOR=0, ICOUNT=0, INODE=0, INODE2=0, IMESH=0, IMESH2=0, IMODE=1, ICORR_TO=0
     INTEGER, DIMENSION(50) :: NTARGET=0
     INTEGER :: STR_INDX=0, STR_SUB_INDX=0
     REAL(EB) :: FED_CO_CO2_O2=0._EB, SOOT_DENS=0._EB, TMP_G=0._EB, RADFLUX=0._EB
     INTEGER :: II=0, JJ=0, KK=0, FED_MESH=0
     LOGICAL :: CHECK_FLOW=.FALSE., EXIT_SIGN=.FALSE., KEEP_XY=.FALSE., SHOW=.TRUE., KNOWN_DOOR=.FALSE.
     LOGICAL :: LOCKED_WHEN_CLOSED=.FALSE., TARGET_WHEN_CLOSED=.FALSE., WAIT_AT_XYZ=.TRUE.
     CHARACTER(60) :: ID='null'
     CHARACTER(60) :: TO_NODE='null'
     CHARACTER(LABEL_LENGTH) :: GRID_NAME='null'
     CHARACTER(LABEL_LENGTH) :: VENT_FFIELD='null'
     INTEGER :: I_VENT_FFIELD=0, Avatar_Color_Index=0, I_EMESH_EXITS=0
     INTEGER, DIMENSION(3) :: RGB=-1
     REAL(EB), DIMENSION(3) :: ORIENTATION=0.0_EB
  END TYPE EVAC_DOOR_TYPE
  !
  ! Like door, but corr will model stairs (or corridors).
  ! The parameters, like velocity as function of density etc.
  ! define if it is corridor or stairway
  ! (&CORR lines)
  TYPE EVAC_CORR_TYPE
     REAL(EB) :: T_first=0._EB, T_last=0._EB, Flow_max=0._EB, Width1=0._EB, Width2=0._EB
     REAL(EB) :: X1=0._EB,X2=0._EB,Y1=0._EB,Y2=0._EB,Z1=0._EB,Z2=0._EB, Width=0._EB
     REAL(EB) :: Eff_Width=0._EB, Eff_Length=0._EB, Eff_Area=0._EB, Fac_Speed=0._EB
     REAL(EB) :: Travel_Time=0._EB, Time_Wait_To_Fill=0._EB, Time_Last_In=0._EB, Time_Stop=0._EB, Time_Start=0._EB
     INTEGER :: Elevator_Mode=1
     ! Note: Corridor may have 2 different points, where smoke etc. data is saved.
     REAL(EB), DIMENSION(2) :: FED_CO_CO2_O2=0._EB, SOOT_DENS=0._EB, TMP_G=0._EB, RADFLUX=0._EB
     INTEGER :: FED_MESH=0, FED_MESH2=0
     INTEGER, DIMENSION(2) :: II=0, JJ=0, KK=0
     INTEGER :: IOR=0, ICOUNT=0, INODE=0, INODE2=0, IMESH=0, IMESH2=0, IDOOR_FROM=0
     INTEGER :: MAX_HUMANS_INSIDE=0, n_inside=0
     LOGICAL :: CHECK_FLOW=.FALSE., ELEVATOR=.FALSE.
     INTEGER, DIMENSION(3) :: RGB=-1
     CHARACTER(60) :: ID='null'
     CHARACTER(60) :: TO_NODE='null'
     CHARACTER(LABEL_LENGTH) :: GRID_NAME='null'
     TYPE (CORR_LL_TYPE), POINTER :: First =>NULL()
  END TYPE EVAC_CORR_TYPE
  !
  ! STRS is a construct to build a staircase. STRS consists of stairs and landings.
  TYPE EVAC_STRS_TYPE
     REAL(EB) :: XB(6)
     REAL(EB), POINTER, DIMENSION(:,:)   :: XB_NODE =>NULL(), XB_CORE =>NULL()
     REAL(EB) :: FAC_V0_HORI=1._EB, FAC_V0_DOWN=1._EB, FAC_V0_UP=1._EB
     INTEGER :: ICOUNT=0, INODE=0, INODE2=0, IMESH=0, IMESH2=0, N_CORES = 0
     INTEGER :: N_LANDINGS, N_NODES, N_NODES_OUT, N_NODES_IN
     INTEGER, POINTER, DIMENSION(:) :: NODE_IOR =>NULL(), NODE_TYPE =>NULL(), NODES_IN =>NULL()
     INTEGER, POINTER, DIMENSION(:) :: NODES_OUT =>NULL(), I_CORE =>NULL()
     CHARACTER(60) :: ID
     CHARACTER(LABEL_LENGTH) :: MESH_ID
     LOGICAL RIGHT_HANDED
  END TYPE EVAC_STRS_TYPE
  !
  ! This produces more humans on the floor specified by the coordinates.
  ! Person types ('soccer_fan' etc) are also defined here for these persons.
  ! (&ENTR lines)
  TYPE EVAC_ENTR_TYPE
     REAL(EB) :: T_first=0._EB, T_last=0._EB, Flow=0._EB, Width=0._EB, T_Start=0._EB, T_Stop=0._EB
     REAL(EB) :: X1=0._EB,X2=0._EB,Y1=0._EB,Y2=0._EB,Z1=0._EB,Z2=0._EB, Height=2.0_EB, Z=0.0_EB
     REAL(EB) :: CB_TimeLastRead=0._EB, CB_TimeNextRead=0._EB, Area=0._EB
     INTEGER :: IOR=0, ICOUNT=0, IPC=0, IMESH=0, INODE=0, IMODE=-1, TO_INODE=0, N_Initial=0, Max_Humans=-1, &
          STR_INDX=0, STR_SUB_INDX=0, AVATAR_TYPE_INDEX=0, CB_N_Agents, CB_N_CrowbarAgents
     CHARACTER(60) :: CLASS_NAME='null', ID='null', AVATAR_TYPE_NAME='null'
     CHARACTER(60) :: TO_NODE='null'
     CHARACTER(LABEL_LENGTH) :: GRID_NAME='null', Max_Humans_Ramp, PROP_ID='null', CROWBAR_INPUT_FILE='null'
     LOGICAL :: After_Tpre=.FALSE., No_Persons=.FALSE., SHOW=.TRUE., CROWBAR_READ_IN=.FALSE.
     INTEGER :: N_VENT_FFIELDS=0, Avatar_Color_Index=0, I_AGENT_TYPE=2
     INTEGER, POINTER, DIMENSION(:) :: I_DOOR_NODES =>NULL()
     INTEGER, POINTER, DIMENSION(:) :: I_VENT_FFIELDS =>NULL()
     REAL(EB), POINTER, DIMENSION(:) :: P_VENT_FFIELDS =>NULL()
     INTEGER, DIMENSION(3) :: RGB=-1, AVATAR_RGB=-1
     REAL(EB), DIMENSION(3) :: ORIENTATION=0.0_EB
  END TYPE EVAC_ENTR_TYPE
  !
  TYPE EVAC_NODE_TYPE
     INTEGER :: Node_Index=0, IMESH=0
     CHARACTER(60) :: ID='null', Node_Type='null'
     CHARACTER(LABEL_LENGTH) :: GRID_NAME='null'
  END TYPE EVAC_NODE_TYPE
  !
  ! Linked list, needed for the corridors
  TYPE CORR_LL_TYPE
     TYPE (HUMAN_TYPE) :: HUMAN
     REAL(EB) :: T_in=0._EB, T_out=0._EB
     LOGICAL :: From1_To2=.FALSE.
     INTEGER :: Index=0
     TYPE (CORR_LL_TYPE), POINTER :: Next =>NULL()
  END TYPE CORR_LL_TYPE
  !
  ! Save device data for the evacuation calculation
  TYPE EVAC_DEVC_TYPE
     REAL(EB) :: T_Change=0._EB
     LOGICAL :: CURRENT=.FALSE., PRIOR=.FALSE., USE_NOW=.FALSE.
     INTEGER :: I_Devc=0, I_Devc_Evac=0, I_Type=0
     CHARACTER(LABEL_LENGTH) :: DEVC_ID='null'
  END TYPE EVAC_DEVC_TYPE
  !
  ! The EDEV namelists, which contain information how the agents
  ! act to the device and/or control information.
  TYPE EVAC_EDEV_TYPE
     CHARACTER(60) :: ID='null', PERS_ID='null', EVAC_ID='null'
     CHARACTER(LABEL_LENGTH) :: MESH_ID='null'
     CHARACTER(LABEL_LENGTH), POINTER, DIMENSION(:) :: INPUT_ID =>NULL()
     INTEGER, POINTER, DIMENSION(:) :: INPUT_DEVC_INDEX =>NULL()
     INTEGER :: I_PRE_DIST=0, IMESH=0, N_INPUTS=0
     REAL(EB) :: Tpre_mean=0._EB, Tpre_para=0._EB, Tpre_para2=0._EB, Tpre_low=0._EB, Tpre_high=0._EB
     REAL(EB) :: TIME_DELAY=0._EB, PROB=0._EB
     LOGICAL :: GLOBAL=.TRUE.
  END TYPE EVAC_EDEV_TYPE
  !
  !
  ! Next holds door information for groups
  TYPE (KNOWN_DOOR_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: Group_Known_Doors
  ! Next holds door information for lonely humans (group_id=0)
  TYPE (KNOWN_DOOR_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: Human_Known_Doors
  INTEGER :: ilh, ilh_dim

  ! Holds the list of the different human groups, i33 is a running index
  ! for the groups, i33_dim is last index, i.e., the dimension of the array.
  TYPE (GROUP_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: Group_List
  INTEGER :: i33, i33_dim

  ! Holds the information of the nodes
  TYPE (EVAC_NODE_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: Evac_Node_List

  ! Holds the information of the EVAC-lines.
  TYPE (EVACUATION_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: EVAC_EVACS

  ! Holds the information of the EVHO-lines.
  TYPE (EVAC_HOLE_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: EVAC_HOLES

  ! Holds the information of the EVSS-lines.
  TYPE (EVAC_SSTAND_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: EVAC_SSTANDS

  ! Holds the information of the EXIT-lines.
  TYPE (EVAC_EXIT_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: EVAC_EXITS

  ! Holds the information of the DOOR-lines.
  TYPE (EVAC_DOOR_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: EVAC_DOORS

  ! Holds the information of the ENTR-lines.
  TYPE (EVAC_ENTR_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: EVAC_ENTRYS

  ! Holds the information of the CORR-lines.
  TYPE (EVAC_CORR_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: EVAC_CORRS

  ! Holds the information on the STRS-lines.
  TYPE (EVAC_STRS_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: EVAC_STRS

  ! Holds the information of the PERS-lines.
  TYPE (EVAC_PERS_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: EVAC_PERSON_CLASSES

  ! Holds the information of the EDEV-lines.
  TYPE (EVAC_EDEV_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: EVAC_EDEV

  ! Holds the information of the devices used in the evacuation calculation.
  TYPE (EVAC_DEVC_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: EVAC_DEVICES
  !
  ! Next are needed for the Gaussian random numbers
  INTEGER GaussFlag
  REAL(EB) GaussSet1, GaussSet2
  INTEGER GTrunFlag
  REAL(EB) GTrunSet1, GTrunSet2
  !
  INTEGER :: NPC_EVAC, NPC_PERS, N_EXITS, N_DOORS, N_ENTRYS, &
       N_CORRS, N_EGRIDS, N_NODES, N_HOLES, N_SSTANDS, N_STRS, N_CO_EXITS, N_DEVC_EVAC, N_EDEV, N_DOOR_MESHES
  INTEGER :: NPPS
  INTEGER :: ILABEL_last, I_FED_FILE_FORMAT=-3
  ! I_FED_FILE_FORMAT: -1 Mesh, Corr, Door, and Exit data on the FED file
  ! I_FED_FILE_FORMAT: -2 Mesh and Corr data on the FED file
  ! I_FED_FILE_FORMAT: -3 Mesh, Corr, and Devc data on the FED file
  CHARACTER(MESSAGE_LENGTH) :: MESSAGE
  REAL(FB) :: EVAC_Z_MIN, EVAC_Z_MAX

  ! Main evacuation mesh information that is read before READ_MESH in read.f90.
  ! These are needed to automate the generation of the additional door flow meshes.
  !
  ! Type for exit+door information needed in READ_MESH to define the additional door flow meshes.
  TYPE EVAC_EMESH_EXITS_TYPE
     REAL(EB), DIMENSION(6) :: XB=0.0_EB
     REAL(EB), POINTER, DIMENSION(:,:) :: U_EVAC =>NULL(), V_EVAC =>NULL()
     INTEGER :: EMESH=0, NEXIT=0, IOR=0, IMESH=0, I_OBST=0, I_VENT=0, MAINMESH=0, I_DOORS_EMESH=0
     INTEGER :: IBAR=0, JBAR=0, KBAR=0
     INTEGER, DIMENSION(3) :: RGB=-1
     LOGICAL :: IS_EXIT=.FALSE., DEFINE_MESH=.TRUE.
     CHARACTER(LABEL_LENGTH) :: ID='null'
  END TYPE EVAC_EMESH_EXITS_TYPE
  ! Type for strs information needed in READ_MESH to define the main evacuation meshes for stairs.
  TYPE EVAC_EMESH_STAIRS_TYPE
     REAL(EB), DIMENSION(6) :: XB=0.0_EB, XB_CORE=0.0_EB
     REAL(EB) :: EVAC_Z_OFFSET
     INTEGER :: EMESH=0, NSTRS=0, IMESH=0, N_CORES=0, N_LANDINGS=0, IBAR=0, JBAR=0, KBAR=0, I_OBST=0, I_HOLE=0
     INTEGER, DIMENSION(3) :: RGB=-1
     LOGICAL :: DEFINE_MESH=.TRUE.
     CHARACTER(LABEL_LENGTH) :: ID='null'
  END TYPE EVAC_EMESH_STAIRS_TYPE
  CHARACTER(LABEL_LENGTH), DIMENSION(:), ALLOCATABLE :: EMESH_ID
  INTEGER, DIMENSION(:), ALLOCATABLE :: EMESH_NM
  INTEGER, DIMENSION(:), ALLOCATABLE :: EMESH_NFIELDS
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: EMESH_IJK
  REAL(EB), DIMENSION(:,:), ALLOCATABLE :: EMESH_XB
  INTEGER, DIMENSION(:), ALLOCATABLE :: EMESH_INDEX
  ! Holds the VENT information of the EXIT/DOOR lines for READ_MESH etc.
  TYPE (EVAC_EMESH_EXITS_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: EMESH_EXITS
  TYPE (EVAC_EMESH_STAIRS_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: EMESH_STAIRS

  !
  !
  !
  LOGICAL :: NOT_RANDOM, L_FALLING_MODEL=.FALSE., DISCARD_SMOKE_INFO=.FALSE.

  INTEGER :: I_FRIC_SW, COLOR_METHOD, COLOR_METHOD_TMP, I_AVATAR_COLOR, MAX_HUMANS_DIM, SMOKE_KS_SPEED_FUNCTION, &
             FED_ACTIVITY, I_HERDING_TYPE
  REAL(EB) :: EVAC_MASS_EXTINCTION_COEFF
  REAL(EB) ::  FAC_A_WALL, FAC_B_WALL, LAMBDA_WALL, &
       NOISEME, NOISETH, NOISECM, RADIUS_COMPLETE_0, &
       RADIUS_COMPLETE_1, GROUP_EFF, FED_DOOR_CRIT, &
       TDET_SMOKE_DENS, DENS_INIT, EVAC_DT_MAX, GROUP_DENS, &
       FC_DAMPING, EVAC_DT_MIN, V_MAX, V_ANGULAR_MAX, V_ANGULAR, &
       SMOKE_MIN_SPEED, SMOKE_MIN_SPEED_FACTOR, SMOKE_MIN_SPEED_VISIBILITY, TAU_CHANGE_DOOR, &
       HUMAN_SMOKE_HEIGHT, TAU_CHANGE_V0, THETA_SECTOR, CONST_DF, FAC_DF, &
       CONST_CF, FAC_CF, FAC_1_WALL, FAC_2_WALL, FAC_V0_DIR, FAC_V0_NOCF, FAC_NOCF, &
       CF_MIN_A, CF_FAC_A_WALL, CF_MIN_TAU, CF_MIN_TAU_INER, CF_FAC_TAUS, FAC_DOOR_QUEUE, FAC_DOOR_ALPHA,&
       FAC_DOOR_WAIT, CF_MIN_B, FAC_DOOR_OLD, FAC_DOOR_OLD2, R_HERDING, W0_HERDING, WR_HERDING, &
       DOT_HERDING, EVAC_DELTA_SEE, C_HAWK, R_HAWK_DOVE, A_HAWK_T_FACTOR, HERDING_TAU_FACTOR, &
       ALPHA_HAWK, DELTA_HAWK, EPSILON_HAWK, THETA_HAWK, A_HAWK_T_START, T_STOP_HD_GAME, &
       F_MIN_FALL, F_MAX_FALL, D_OVERLAP_FALL, TAU_FALL_DOWN, A_FAC_FALLEN, TIME_FALL_DOWN, PROB_FALL_DOWN, &
       T_ASET_HAWK, T_0_HAWK, T_ASET_TFAC_HAWK, MAX_INITIAL_OVERLAP, TIME_INIT_NERVOUSNESS, &
       SMOKE_SPEED_ALPHA, SMOKE_SPEED_BETA, CROWBAR_DT_READ, NASH_CLOSE_ENOUGH, HERDING_TPRE_PROB
  INTEGER, DIMENSION(3) :: DEAD_RGB
  !
  REAL(EB), DIMENSION(:), ALLOCATABLE :: Tsteps
  !
  REAL(EB), DIMENSION(:), ALLOCATABLE :: K_ave_Door
  REAL(EB), DIMENSION(:), ALLOCATABLE :: FED_max_Door
  LOGICAL, DIMENSION(:), ALLOCATABLE :: Is_Known_Door, Is_Visible_Door
  LOGICAL, DIMENSION(:,:), ALLOCATABLE :: Is_BeeLine_Door, Is_BeeLineXB_Door
  INTEGER, DIMENSION(:), ALLOCATABLE :: Color_Tmp
  !
  INTEGER :: n_dead=0, icyc_old=0, n_change_doors=0, n_change_trials=0
  REAL(EB) :: fed_max_alive, fed_max
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: N_HawkDoveCount
  !
  ! Stairs constants
  INTEGER :: STRS_LANDING_TYPE=1, STRS_STAIR_TYPE=2

 ! Human constants
  INTEGER :: HUMAN_NO_TARGET=0,&
             HUMAN_TARGET_UNSPECIFIED = 1, &
             HUMAN_ANOTHER_MESH_TARGET=2,&
             HUMAN_STRS_TARGET=4
         ! hr%ior = 0: not entering a door
         !          1: moving but target not specified
         !          2: put to an another mesh (target is door/entry)
         !          3: target is corridor (remove from floor)
         !          4: not used (floor node...)
         !          5: target is exit (remove from floor)
         !         -1: can not move to the target node
         !         -2: move to door/entry on the same floor
CONTAINS
  !
  SUBROUTINE READ_EVAC(IMODE)
    USE OUTPUT_DATA, ONLY: COLOR2RGB
    IMPLICIT NONE
    !
    ! IMODE = 1:
    ! Read some evacuation input from the main evacuation meshes:
    ! Read the ijk and XB info of the main evacuation meshes.
    ! Read the door/exit information and check which main evacuation mesh.
    ! Calculate the number of additional door flow meshes and give names to them (doormesh_ID), where
    ! ID is the ID of the door.  Make the 'outflow' VENTs also (xb, ior, mesh_id, mesh_#).  The
    ! additional door flow meshes are numbered starting from the last mesh read in.
    !
    ! IMODE = 2:
    ! Read the evacuation input.
    !
    ! Passed variables
    INTEGER, INTENT(IN) :: IMODE
    !
    ! Local variables
    INTEGER :: NUMBER_INITIAL_PERSONS, SAMPLING_FACTOR, IPC, n_tmp, GN_MIN, GN_MAX, N_RAMP_INI
    REAL(EB) :: DTSAM
    LOGICAL :: EVACFILE, DEFINE_MESH

    REAL(EB) :: DUMMY
    REAL(EB) :: XB(6), XB1(6), XB2(6)
    REAL(EB), DIMENSION(3) :: XYZ, XYZ_SMOKE
    INTEGER :: IOS, IZERO, N, I, J, K, IOR
    CHARACTER(LABEL_LENGTH) QUANTITY, MAX_HUMANS_RAMP, INPUT_ID(40)
    CHARACTER(60) FYI,ID,PERS_ID,TO_NODE,EVAC_ID, DEFAULT_PROPERTIES, AVATAR_TYPE
    CHARACTER(LABEL_LENGTH) FLOW_FIELD_ID
    INTEGER :: DIAMETER_DIST,VELOCITY_DIST,PRE_EVAC_DIST,DET_EVAC_DIST,TAU_EVAC_DIST, GUARD_MEN_IN
    REAL(EB) :: VEL_MEAN,VEL_PARA,VEL_PARA2,VEL_LOW,VEL_HIGH, &
         DIA_MEAN,DIA_PARA,DIA_PARA2,DIA_LOW,DIA_HIGH, &
         PRE_MEAN,PRE_PARA,PRE_PARA2,PRE_LOW,PRE_HIGH, &
         DET_MEAN,DET_PARA,DET_PARA2,DET_LOW,DET_HIGH, &
         TAU_MEAN,TAU_PARA,TAU_PARA2,TAU_LOW,TAU_HIGH, &
         FCONST_A,FCONST_B,L_NON_SP,C_YOUNG,GAMMA,KAPPA,ANGLE, &
         D_TORSO_MEAN,D_SHOULDER_MEAN, TAU_ROT, M_INERTIA, TARGET_X, TARGET_Y, &
         DELTA_X, DELTA_Y, MAXIMUM_V0_FACTOR, TIME_START_FED, MASS_OF_AGENT
    !Issue1547: The MAXIMUM_V0_FACTOR should be difined as a variable in Fortran, because it is used
    !Issue1547: in a namelist (PERS-namelist) to read an user input, if given.
    INTEGER :: MAX_HUMANS_INSIDE, n_max_in_corrs, COLOR_INDEX, MAX_HUMANS, AGENT_TYPE
    REAL(EB) :: MAX_FLOW, WIDTH, TIME_START, TIME_STOP, WIDTH1, &
         WIDTH2, EFF_WIDTH, EFF_LENGTH, FAC_SPEED, TIME_OPEN, TIME_CLOSE
    REAL(EB) :: UBAR0, VBAR0, TIME_DELAY, TRAVEL_TIME
    LOGICAL :: CHECK_FLOW, COUNT_ONLY, AFTER_REACTION_TIME, EXIT_SIGN, KEEP_XY, USE_V0, SHOW, COUNT_DENSITY, GLOBAL
    LOGICAL :: OUTPUT_SPEED, OUTPUT_MOTIVE_FORCE, OUTPUT_FED, OUTPUT_OMEGA, OUTPUT_DENSITY, OUTPUT_ANGLE, &
         OUTPUT_CONTACT_FORCE, OUTPUT_TOTAL_FORCE, OUTPUT_MOTIVE_ANGLE, OUTPUT_ACCELERATION, OUTPUT_NERVOUSNESS, &
         LOCKED_WHEN_CLOSED, TARGET_WHEN_CLOSED, WAIT_AT_XYZ, ELEVATOR, EVAC_FDS6, KNOWN_DOOR
    !Issue1547: Added new output keyword for the PERS namelist, here is the definition
    !Issue1547: of the variable (OUTPUT_NERVOUSNES) that is needed by Fortran.
    INTEGER, DIMENSION(3) :: RGB, AVATAR_RGB
    CHARACTER(LABEL_LENGTH) :: VENT_FFIELD, MESH_ID, EVAC_MESH, PROP_ID, CROWBAR_INPUT_FILE
    REAL(EB) :: FAC_V0_UP, FAC_V0_DOWN, FAC_V0_HORI, HEIGHT, HEIGHT0, ESC_SPEED, PROB
    CHARACTER(25) :: COLOR, DEAD_COLOR, AVATAR_COLOR

    ! Stairs variables
    REAL(EB) :: XB_CORE(4), XB_CORES(500,6), XB_LANDINGS(500,6), XB_STAIRS(500,8)
    REAL(EB) VERTICAL_LANDING_SEPARATION, STR_Length, STR_Height
    INTEGER N_LANDINGS, NL, NODES_TMP(500)
    LOGICAL :: RIGHT_HANDED, LEFT_HANDED

    CHARACTER(60), DIMENSION(51) :: KNOWN_DOOR_NAMES
    REAL(EB), DIMENSION(51) :: KNOWN_DOOR_PROBS

    INTEGER, DIMENSION(:), ALLOCATABLE :: TMP_AVATAR_TYPE_INDEX
    CHARACTER(60), DIMENSION(:), ALLOCATABLE :: TMP_AVATAR_TYPE_NAME, TMP_AVATAR_TYPE_PROP

    INTEGER :: ii,jj,kk

    INTEGER :: IJK(3)
    REAL(EB) :: EVAC_Z_OFFSET
    INTEGER :: N_CORES
    INTEGER :: N_AVATAR_TYPE
    LOGICAL :: CROWBAR_DUMP

    NAMELIST /EXIT/ ID, XB, IOR, FLOW_FIELD_ID, CHECK_FLOW, &
         MAX_FLOW, FYI, COUNT_ONLY, WIDTH, XYZ, VENT_FFIELD, COUNT_DENSITY, &
         MESH_ID, COLOR_INDEX, XYZ_SMOKE, EVAC_ID, PERS_ID, KNOWN_DOOR, &
         TIME_OPEN, TIME_CLOSE, EVAC_MESH, RGB, COLOR, SHOW, HEIGHT, LOCKED_WHEN_CLOSED, TARGET_WHEN_CLOSED, WAIT_AT_XYZ
    NAMELIST /DOOR/ ID, XB, IOR, FLOW_FIELD_ID, CHECK_FLOW, &
         MAX_FLOW, TO_NODE, FYI, WIDTH, XYZ, VENT_FFIELD, KNOWN_DOOR, &
         EXIT_SIGN, MESH_ID, COLOR_INDEX, XYZ_SMOKE, KEEP_XY, &
         TIME_OPEN, TIME_CLOSE, EVAC_MESH, RGB, COLOR, SHOW, HEIGHT, LOCKED_WHEN_CLOSED, TARGET_WHEN_CLOSED, WAIT_AT_XYZ
    NAMELIST /ENTR/ ID, XB, IOR, FLOW_FIELD_ID, MAX_FLOW, &
         FYI, WIDTH, QUANTITY, PERS_ID, TIME_START, &
         TIME_STOP, AFTER_REACTION_TIME, &
         KNOWN_DOOR_NAMES, KNOWN_DOOR_PROBS, &
         MESH_ID, COLOR_INDEX, EVAC_MESH, RGB, COLOR, &
         AVATAR_COLOR, AVATAR_RGB, MAX_HUMANS, MAX_HUMANS_RAMP, SHOW, HEIGHT, AGENT_TYPE, AVATAR_TYPE, PROP_ID, &
         CROWBAR_INPUT_FILE
    NAMELIST /CORR/ ID, XB, IOR, FLOW_FIELD_ID, CHECK_FLOW, &
         MAX_FLOW, TO_NODE, FYI, WIDTH, WIDTH1, WIDTH2, &
         EFF_WIDTH, EFF_LENGTH, MAX_HUMANS_INSIDE, FAC_SPEED, &
         XB1, XB2, RGB, COLOR, ELEVATOR, TRAVEL_TIME
    NAMELIST /STRS/ ID, XB, XB_CORE, XB_CORES, TO_NODE, RIGHT_HANDED, LEFT_HANDED, MESH_ID, &
         N_LANDINGS, XB_LANDINGS, VERTICAL_LANDING_SEPARATION, &
         FAC_V0_UP, FAC_V0_DOWN, FAC_V0_HORI, RGB, COLOR, IJK, EVAC_Z_OFFSET
    NAMELIST /EVAC/ NUMBER_INITIAL_PERSONS, QUANTITY, FYI, &
         ID, DTSAM, XB, FLOW_FIELD_ID, PERS_ID, &
         TIME_START, TIME_STOP, IOR, MAX_FLOW, WIDTH, ANGLE, &
         AFTER_REACTION_TIME, GN_MIN, GN_MAX, &
         KNOWN_DOOR_NAMES, KNOWN_DOOR_PROBS, MESH_ID, &
         COLOR_INDEX, EVAC_MESH, RGB, COLOR, &
         AVATAR_COLOR, AVATAR_RGB, SHOW, PRE_EVAC_DIST, DET_EVAC_DIST, &
         PRE_MEAN,PRE_PARA,PRE_PARA2,PRE_LOW,PRE_HIGH, &
         DET_MEAN,DET_PARA,DET_PARA2,DET_LOW,DET_HIGH, AGENT_TYPE, AVATAR_TYPE, PROP_ID, TIME_FALL_DOWN, &
         GUARD_MEN_IN, TARGET_X, TARGET_Y, DELTA_X, DELTA_Y, TIME_START_FED, CROWBAR_INPUT_FILE
    NAMELIST /EVHO/ FYI, ID, XB, EVAC_ID, PERS_ID, MESH_ID, EVAC_MESH, RGB, COLOR, SHOW, TIME_FALL_DOWN

    NAMELIST /EVSS/ FYI, ID, XB, MESH_ID, HEIGHT, HEIGHT0, IOR, &
         FAC_V0_UP, FAC_V0_DOWN, FAC_V0_HORI, ESC_SPEED, EVAC_MESH, RGB, COLOR, &
         UBAR0, VBAR0, USE_V0, SHOW, VENT_FFIELD, COUNT_ONLY

    NAMELIST /PERS/ FYI,ID,DIAMETER_DIST,VELOCITY_DIST, &
         PRE_EVAC_DIST,DET_EVAC_DIST,TAU_EVAC_DIST, &
         VEL_MEAN,VEL_PARA,VEL_PARA2,VEL_LOW,VEL_HIGH, &
         DIA_MEAN,DIA_PARA,DIA_PARA2,DIA_LOW,DIA_HIGH, &
         PRE_MEAN,PRE_PARA,PRE_PARA2,PRE_LOW,PRE_HIGH, &
         DET_MEAN,DET_PARA,DET_PARA2,DET_LOW,DET_HIGH, &
         TAU_MEAN,TAU_PARA,TAU_PARA2,TAU_LOW,TAU_HIGH, &
         FCONST_A,FCONST_B,L_NON_SP, C_YOUNG,GAMMA,KAPPA, GROUP_DENS, &
         FAC_A_WALL, FAC_B_WALL, LAMBDA_WALL, NOISEME, NOISETH, NOISECM, &
         I_FRIC_SW, GROUP_EFF, RADIUS_COMPLETE_0, &
         RADIUS_COMPLETE_1, DEFAULT_PROPERTIES, &
         NOT_RANDOM, FED_DOOR_CRIT, COLOR_METHOD, &
         TDET_SMOKE_DENS, DENS_INIT, EVAC_DT_MAX, EVAC_DT_MIN, &
         D_TORSO_MEAN, D_SHOULDER_MEAN, TAU_ROT, M_INERTIA, &
         FC_DAMPING, V_MAX, V_ANGULAR_MAX, V_ANGULAR, &
         OUTPUT_SPEED, OUTPUT_MOTIVE_FORCE, OUTPUT_FED, OUTPUT_OMEGA, OUTPUT_DENSITY, &
         OUTPUT_MOTIVE_ANGLE, OUTPUT_ANGLE, OUTPUT_CONTACT_FORCE, OUTPUT_TOTAL_FORCE, &
         OUTPUT_ACCELERATION, OUTPUT_NERVOUSNESS, COLOR_INDEX, DEAD_RGB, DEAD_COLOR, &
         SMOKE_MIN_SPEED, SMOKE_MIN_SPEED_FACTOR, SMOKE_MIN_SPEED_VISIBILITY, &
         TAU_CHANGE_DOOR, RGB, COLOR, AVATAR_COLOR, AVATAR_RGB, HUMAN_SMOKE_HEIGHT, &
         TAU_CHANGE_V0, THETA_SECTOR, CONST_DF, FAC_DF, CONST_CF, FAC_CF, &
         FAC_1_WALL, FAC_2_WALL, FAC_V0_DIR, FAC_V0_NOCF, FAC_NOCF, &
         CF_MIN_A, CF_FAC_A_WALL, CF_MIN_TAU, CF_MIN_TAU_INER, CF_FAC_TAUS, &
         FAC_DOOR_QUEUE, FAC_DOOR_ALPHA, FAC_DOOR_WAIT, CF_MIN_B, NASH_CLOSE_ENOUGH, &
         FAC_V0_UP, FAC_V0_DOWN, FAC_V0_HORI, FAC_DOOR_OLD, FAC_DOOR_OLD2, &
         R_HERDING, W0_HERDING, WR_HERDING, I_HERDING_TYPE, DOT_HERDING, EVAC_FDS6, EVAC_DELTA_SEE, &
         MAX_HUMANS_DIM, C_HAWK, R_HAWK_DOVE, A_HAWK_T_FACTOR, HERDING_TAU_FACTOR, &
         ALPHA_HAWK, DELTA_HAWK, EPSILON_HAWK, THETA_HAWK, T_STOP_HD_GAME, A_HAWK_T_START, &
         F_MIN_FALL, F_MAX_FALL, D_OVERLAP_FALL, TAU_FALL_DOWN, A_FAC_FALLEN, PROB_FALL_DOWN, &
         T_ASET_HAWK, T_0_HAWK, T_ASET_TFAC_HAWK, &
         MAXIMUM_V0_FACTOR, MAX_INITIAL_OVERLAP, TIME_INIT_NERVOUSNESS, &
         SMOKE_SPEED_ALPHA, SMOKE_SPEED_BETA, SMOKE_KS_SPEED_FUNCTION, FED_ACTIVITY, &
         CROWBAR_DT_READ, MASS_OF_AGENT, DISCARD_SMOKE_INFO, HERDING_TPRE_PROB
    !Issue1547: Added new output keyword for the PERS namelist, here the new output
    !Issue1547: keyword OUTPUT_NERVOUSNES is added to the namelist. Also the user input
    !Issue1547: for the social force MAXIMUM_V0_FACTOR is added to the namelist.
    !
    NAMELIST /EDEV/ FYI, ID, TIME_DELAY, GLOBAL, EVAC_ID, PERS_ID, MESH_ID, INPUT_ID, &
         PRE_EVAC_DIST, PRE_MEAN, PRE_PARA, PRE_PARA2, PRE_LOW, PRE_HIGH, PROB
    !
    IF_IMODE_1: IF (IMODE == 1) THEN
       N_DOOR_MESHES = 0
       IF (NO_EVACUATION) THEN
          N_EVAC = 0
          RETURN ! skip evacuation calculation
       END IF

       CALL COUNT_EVAC_NODES(1) ; IF (STOP_STATUS==SETUP_STOP) RETURN

       ALLOCATE(EMESH_NM(     MAX(1,INPUT_EVAC_GRIDS)), STAT=IZERO)
       CALL ChkMemErr('READ_EVAC','EMESH_NM',IZERO)
       ALLOCATE(EMESH_NFIELDS(MAX(1,INPUT_EVAC_GRIDS+N_STRS)), STAT=IZERO)
       CALL ChkMemErr('READ_EVAC','EMESH_NFIELDS',IZERO)
       EMESH_NFIELDS = 0
       ! Array for the exit+door information (not for the count only exits)
       ALLOCATE(EMESH_EXITS(MAX(1,N_EXITS-N_CO_EXITS+N_DOORS)), STAT=IZERO)
       CALL ChkMemErr('READ_EVAC','EMESH_EXITS',IZERO)
       ! Array for the strs information
       ALLOCATE(EMESH_STAIRS(MAX(1,N_STRS)), STAT=IZERO)
       CALL ChkMemErr('READ_EVAC','EMESH_STAIRS',IZERO)

       ! Read the STRS line
       READ_STRS_LOOP: DO N = 1, N_STRS
          !
          ! Set MESH defaults
          IJK(1)= 10
          IJK(2)= 10
          IJK(3)= 1
          ID                          = 'null'
          XB                          = 0._EB
          XB_CORE                     = 0._EB
          XB_CORES                    = -9999999999999._EB
          XB_LANDINGS                 = 0._EB
          MESH_ID                     = 'null'
          N_LANDINGS                  = 0
          VERTICAL_LANDING_SEPARATION = 0._EB
          EVAC_Z_OFFSET               = 1.0_EB
          !
          CALL CHECKREAD('STRS',LU_INPUT,IOS)  ; IF (STOP_STATUS==SETUP_STOP) RETURN
          IF (IOS == 1) THEN
             EXIT READ_STRS_LOOP
          END IF
          READ(LU_INPUT,STRS,END=32,IOSTAT=IOS)

          ! Reorder XB coordinates if necessary
          CALL CHECK_XB(XB)

          ! Colors, integer RGB(3), e.g., (23,255,0)
          IF (ANY(RGB < 0) .AND. COLOR=='null') COLOR = 'BLACK'
          IF (COLOR /= 'null') CALL COLOR2RGB(RGB,COLOR)

          ! Count number of cores
          N_CORES = 0
          DO I = 1,500
             IF (ANY(XB_CORES(I,:)>-9999999999998._EB)) N_CORES = N_CORES + 1
          ENDDO

          IF (N_LANDINGS>500) THEN
             WRITE(MESSAGE,'(A,A,A)') 'ERROR: STRS ',TRIM(ID),' N_LANDINGS > 500'
             CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
          END IF
          !N_LANDINGS = N_LANDINGS
          !N_NODES = 2*N_LANDINGS - 1

          IF (N_CORES == 0) THEN
             N_CORES = 1
             XB_CORES(1,1:4) = XB_CORE(1:4)
             XB_CORES(1,5) = XB(5)
             XB_CORES(1,6) = XB(6)
             IF (ALL(XB_CORE==0._EB)) THEN
                WRITE(MESSAGE,'(3A)') 'ERROR: STRS object ', TRIM(ID), ' has no XB_CORE defined.'
                CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
             ENDIF
          ELSE
             ! For now, just the first one is used (this is under construction)
             XB_CORES(1,1:6) = XB_CORES(1,1:6)
             !   XB_CORE(1:N_CORES,1:6) = XB_CORES(1:N_CORES,1:6)
          ENDIF
          DEFINE_MESH = .TRUE.
          IF (TRIM(MESH_ID)/='null') DEFINE_MESH = .FALSE.

          EMESH_STAIRS(N)%IBAR = IJK(1)
          EMESH_STAIRS(N)%JBAR = IJK(2)
          EMESH_STAIRS(N)%KBAR = IJK(3) ! kbar = 1 is forced in read_mesh (for now)
          EMESH_STAIRS(N)%EVAC_Z_OFFSET = EVAC_Z_OFFSET  ! Not yet used properly in read_mesh???

          EMESH_STAIRS(N)%N_CORES = N_CORES
          EMESH_STAIRS(N)%N_LANDINGS = N_LANDINGS
          EMESH_STAIRS(N)%NSTRS = N
          EMESH_STAIRS(N)%RGB = RGB
          EMESH_STAIRS(N)%ID = TRIM(ID)
          EMESH_STAIRS(N)%DEFINE_MESH = DEFINE_MESH
          EMESH_STAIRS(N)%XB = XB
          EMESH_STAIRS(N)%XB_CORE(1:6) = XB_CORES(1,1:6)

       END DO READ_STRS_LOOP
32     REWIND(LU_INPUT)


       N_TMP = 0    ! N_TMP is exit+door index without count only exits
       READ_EXIT_LOOP_0: DO N = 1, N_EXITS
          DEFINE_MESH   = .TRUE.
          ID            = 'null'
          RGB           = -1
          COLOR         = 'null'
          XB            = 0.0_EB
          IOR           = 0
          MESH_ID       = 'null'
          EVAC_MESH     = 'null'
          VENT_FFIELD   = 'null'
          COUNT_ONLY    = .FALSE.
          COUNT_DENSITY = .FALSE.
          COLOR_INDEX   = -1

          CALL CHECKREAD('EXIT', LU_INPUT, IOS)  ; IF (STOP_STATUS==SETUP_STOP) RETURN
          IF (IOS == 1) THEN
             EXIT READ_EXIT_LOOP_0
          END IF
          READ(LU_INPUT, EXIT, END=26, IOSTAT=IOS)
          IF (VENT_FFIELD /= 'null') DEFINE_MESH = .FALSE. ! Old input, user gives the door flow meshes

          ! Old input used COLOR_INDEX, next lines are needed for that
          IF (MYID==MAX(0,EVAC_PROCESS) .AND. COLOR_INDEX/=-1) WRITE (LU_ERR,'(A,A)') &
               ' WARNING: keyword COLOR_INDEX is replaced by COLOR at EXIT line ', TRIM(ID)
          IF (COLOR_INDEX == 1) COLOR = 'BLACK'
          IF (COLOR_INDEX == 2) COLOR = 'YELLOW'
          IF (COLOR_INDEX == 3) COLOR = 'BLUE'
          IF (COLOR_INDEX == 4) COLOR = 'RED'
          IF (COLOR_INDEX == 5) COLOR = 'GREEN'
          IF (COLOR_INDEX == 6) COLOR = 'MAGENTA'
          IF (COLOR_INDEX == 7) COLOR = 'CYAN'

          IF (COUNT_DENSITY) THEN
             COUNT_ONLY = .TRUE.
             IOR = +3
          END IF

          IF (COUNT_ONLY) CYCLE READ_EXIT_LOOP_0
          N_TMP = N_TMP + 1   ! N_TMP is exit+door index without count only exits
          ! Colors, integer RGB(3), e.g., (23,255,0)
          IF (ANY(RGB < 0) .AND. COLOR=='null') COLOR = 'FOREST GREEN'
          IF (COLOR /= 'null') CALL COLOR2RGB(RGB,COLOR)

          IF (EVAC_MESH /= 'null') THEN
             MESH_ID = TRIM(EVAC_MESH)
             IF (MYID==MAX(0,EVAC_PROCESS)) WRITE (LU_ERR,'(A,A)') &
                  ' WARNING: keyword EVAC_MESH is replaced by MESH_ID at EXIT line ', TRIM(ID)
          END IF

          ! Reorder XB coordinates if necessary
          CALL CHECK_XB(XB)

          ! Check which evacuation mesh, INPUT_EVAC_GRIDS is the number of main evacuation meshes.
          II = 0
          PEX_MeshLoop_0: DO I = 1, INPUT_EVAC_GRIDS
             IF (Is_Within_Bounds(XB(1), XB(2), XB(3), XB(4), XB(5), XB(6), &
                  EMESH_XB(1,I), EMESH_XB(2,I), EMESH_XB(3,I), EMESH_XB(4,I), EMESH_XB(5,I), EMESH_XB(6,I), &
                  0._EB, 0._EB, 0._EB)) THEN
                IF (TRIM(MESH_ID) == 'null' .OR. TRIM(MESH_ID) == TRIM(EMESH_ID(I))) THEN
                   II = II + 1
                   ! N_TMP is exit+door index without count only exits
                   ! N is the exit index, I is the main evacuation mesh index
                   EMESH_EXITS(N_TMP)%EMESH = I
                END IF
             END IF
          END DO PEX_MeshLoop_0
          PEX_StrsLoop_0: DO I = 1, N_STRS
             IF (Is_Within_Bounds(XB(1), XB(2), XB(3), XB(4), XB(5), XB(6), &
                  EMESH_STAIRS(I)%XB(1), EMESH_STAIRS(I)%XB(2), EMESH_STAIRS(I)%XB(3), &
                  EMESH_STAIRS(I)%XB(4), EMESH_STAIRS(I)%XB(5), EMESH_STAIRS(I)%XB(6), &
                  0._EB, 0._EB, 0._EB)) THEN
                IF (TRIM(MESH_ID) == 'null' .OR. TRIM(MESH_ID) == ('Emesh_' // TRIM(EMESH_STAIRS(I)%ID))) THEN
                   II = II + 1
                   ! N_TMP is exit+door index without count only exits
                   ! N is the exit index, I is the main evacuation mesh index
                   EMESH_EXITS(N_TMP)%EMESH = 0
                   DEFINE_MESH = .FALSE.
                END IF
             END IF
          END DO PEX_StrsLoop_0

          IF (II == 0 .AND. MYID==MAX(0,EVAC_PROCESS)) THEN
             WRITE(MESSAGE,'(A,A,A)') 'ERROR: EXIT ',TRIM(ID), ' problem with IMESH, no mesh found'
             CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
          END IF
          IF (II > 1 .AND. MYID==MAX(0,EVAC_PROCESS)) THEN
             WRITE(MESSAGE,'(A,A,A)') 'ERROR: EXIT ',TRIM(ID), ' not an unique mesh found '
             CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
          END IF
          SELECT CASE (IOR)
          CASE (-1,+1)
             IF ( (XB(2)-XB(1)) >= TWO_EPSILON_EB ) THEN
                IF (MYID==MAX(0,EVAC_PROCESS)) THEN
                   WRITE(MESSAGE,'(A,A,A)') 'ERROR: EXIT ',TRIM(ID), ' IOR=+-1 but not a vertical plane '
                   CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
                END IF
             END IF
          CASE (-2,+2)
             IF ( (XB(4)-XB(3)) >= TWO_EPSILON_EB ) THEN
                IF (MYID==MAX(0,EVAC_PROCESS)) THEN
                   WRITE(MESSAGE,'(A,A,A)') 'ERROR: EXIT ',TRIM(ID), ' IOR=+-2 but not a vertical plane '
                   CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
                END IF
             END IF
          CASE Default
             IF (MYID==MAX(0,EVAC_PROCESS)) THEN
                WRITE(MESSAGE,'(A,A,A)') 'ERROR: EXIT ',TRIM(ID), ' No IOR given or not -2,-1,+1,+2 '
                CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
             END IF
          END SELECT
          ! Doors/Exits are always vertical
          IF (EMESH_EXITS(N_TMP)%EMESH>0) THEN
             XB(5) = EMESH_XB(5,EMESH_EXITS(N_TMP)%EMESH)
             XB(6) = EMESH_XB(6,EMESH_EXITS(N_TMP)%EMESH)
          END IF
          !
          ! Save the information needed in READ_MESH, etc.
          EMESH_EXITS(N_TMP)%NEXIT       = N     ! N is the exit line index (including count only ones)
          EMESH_EXITS(N_TMP)%IOR         = IOR
          EMESH_EXITS(N_TMP)%IS_EXIT     = .TRUE.
          EMESH_EXITS(N_TMP)%DEFINE_MESH = DEFINE_MESH  ! Count_only exits are not included in 1...N_TMP list
          EMESH_EXITS(N_TMP)%ID          = TRIM(ID)
          EMESH_EXITS(N_TMP)%XB(:)       = XB(:)
          EMESH_EXITS(N_TMP)%RGB(:)      = RGB(:)
          EMESH_EXITS(N_TMP)%IMESH       = EMESH_EXITS(N_TMP)%EMESH  ! Updated in READ_MESH, the mesh index
          EMESH_EXITS(N_TMP)%MAINMESH    = EMESH_EXITS(N_TMP)%EMESH  ! Updated in READ_MESH, the main mesh index
          ! EMESH_EXITS(N_TMP)%EMESH is the main evacuation mesh index: EMESH_ID, EMESH_IJK, EMESH_XB
       END DO READ_EXIT_LOOP_0
26     CONTINUE
       REWIND(LU_INPUT)

       READ_DOOR_LOOP_0: DO N = 1, N_DOORS
          DEFINE_MESH   = .TRUE.
          ID            = 'null'
          RGB           = -1
          COLOR         = 'null'
          XB            = 0.0_EB
          IOR           = 0
          MESH_ID       = 'null'
          EVAC_MESH     = 'null'
          TO_NODE       = 'null'
          VENT_FFIELD   = 'null'
          COLOR_INDEX   = -1

          CALL CHECKREAD('DOOR', LU_INPUT, IOS)  ; IF (STOP_STATUS==SETUP_STOP) RETURN
          IF (IOS == 1) THEN
             EXIT READ_DOOR_LOOP_0
          END IF
          READ(LU_INPUT, DOOR, END=27, IOSTAT=IOS)
          IF (VENT_FFIELD /= 'null') DEFINE_MESH = .FALSE. ! Old input, user gives the door flow meshes
          IF (TO_NODE     == 'null') DEFINE_MESH = .FALSE. ! This is more or less like an entry.

          ! Old input used COLOR_INDEX, next lines are needed for that
          IF (MYID==MAX(0,EVAC_PROCESS) .AND. COLOR_INDEX/=-1) WRITE (LU_ERR,'(A,A)') &
               ' WARNING: keyword COLOR_INDEX is replaced by COLOR at DOOR line ', TRIM(ID)
          IF (COLOR_INDEX == 1) COLOR = 'BLACK'
          IF (COLOR_INDEX == 2) COLOR = 'YELLOW'
          IF (COLOR_INDEX == 3) COLOR = 'BLUE'
          IF (COLOR_INDEX == 4) COLOR = 'RED'
          IF (COLOR_INDEX == 5) COLOR = 'GREEN'
          IF (COLOR_INDEX == 6) COLOR = 'MAGENTA'
          IF (COLOR_INDEX == 7) COLOR = 'CYAN'

          N_TMP = N_TMP + 1   ! N_TMP is exit+door index without count only exits
          ! Colors, integer RGB(3), e.g., (23,255,0)
          IF (ANY(RGB < 0) .AND. COLOR=='null') COLOR = 'FOREST GREEN'
          IF (COLOR /= 'null') CALL COLOR2RGB(RGB,COLOR)

          IF (EVAC_MESH /= 'null') THEN
             MESH_ID = TRIM(EVAC_MESH)
             IF (MYID==MAX(0,EVAC_PROCESS)) WRITE (LU_ERR,'(A,A)') &
                  ' WARNING: keyword EVAC_MESH is replaced by MESH_ID at DOOR line ', TRIM(ID)
          END IF

          ! Reorder XB coordinates if necessary
          CALL CHECK_XB(XB)

          ! Check which evacuation mesh, INPUT_EVAC_GRIDS is the number of main evacuation meshes.
          II = 0
          PDX_MeshLoop_0: DO I = 1, INPUT_EVAC_GRIDS
             IF (Is_Within_Bounds(XB(1), XB(2), XB(3), XB(4), XB(5), XB(6), &
                  EMESH_XB(1,I), EMESH_XB(2,I), EMESH_XB(3,I), EMESH_XB(4,I), EMESH_XB(5,I), EMESH_XB(6,I), &
                  0._EB, 0._EB, 0._EB)) THEN
                IF (TRIM(MESH_ID) == 'null' .OR. TRIM(MESH_ID) == TRIM(EMESH_ID(I))) THEN
                   II = II + 1
                   ! N_TMP is exit+door index without count only exits
                   ! N is the door index, I is the main evacuation mesh index
                   EMESH_EXITS(N_TMP)%EMESH = I
                END IF
             END IF
          END DO PDX_MeshLoop_0
          PDX_StrsLoop_0: DO I = 1, N_STRS
             IF (Is_Within_Bounds(XB(1), XB(2), XB(3), XB(4), XB(5), XB(6), &
                  EMESH_STAIRS(I)%XB(1), EMESH_STAIRS(I)%XB(2), EMESH_STAIRS(I)%XB(3), &
                  EMESH_STAIRS(I)%XB(4), EMESH_STAIRS(I)%XB(5), EMESH_STAIRS(I)%XB(6), &
                  0._EB, 0._EB, 0._EB)) THEN
                IF (TRIM(MESH_ID) == 'null' .OR. TRIM(MESH_ID) == ('Emesh_' // TRIM(EMESH_STAIRS(I)%ID))) THEN
                   II = II + 1
                   ! N_TMP is exit+door index without count only exits
                   ! N is the exit index, I is the main evacuation mesh index
                   DEFINE_MESH = .FALSE.
                   EMESH_EXITS(N_TMP)%EMESH = 0
                END IF
             END IF
          END DO PDX_StrsLoop_0

          IF (II == 0 .AND. MYID==MAX(0,EVAC_PROCESS)) THEN
             WRITE(MESSAGE,'(A,A,A)') 'ERROR: DOOR ',TRIM(ID), ' problem with IMESH, no mesh found'
             CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
          END IF
          IF (II > 1 .AND. MYID==MAX(0,EVAC_PROCESS)) THEN
             WRITE(MESSAGE,'(A,A,A)') 'ERROR: DOOR ',TRIM(ID), ' not an unique mesh found '
             CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
          END IF
          SELECT CASE (IOR)
          CASE (-1,+1)
             IF ( (XB(2)-XB(1)) >= TWO_EPSILON_EB ) THEN
                IF (MYID==MAX(0,EVAC_PROCESS)) THEN
                   WRITE(MESSAGE,'(A,A,A)') 'ERROR: DOOR ',TRIM(ID), ' IOR=+-1 but not a vertical plane '
                   CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
                END IF
             END IF
          CASE (-2,+2)
             IF ( (XB(4)-XB(3)) >= TWO_EPSILON_EB ) THEN
                IF (MYID==MAX(0,EVAC_PROCESS)) THEN
                   WRITE(MESSAGE,'(A,A,A)') 'ERROR: DOOR ',TRIM(ID), ' IOR=+-2 but not a vertical plane '
                   CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
                END IF
             END IF
          CASE Default
             IF (MYID==MAX(0,EVAC_PROCESS)) THEN
                WRITE(MESSAGE,'(A,A,A)') 'ERROR: DOOR ',TRIM(ID), ' No IOR given or not -2,-1,+1,+2 '
                CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
             END IF
          END SELECT
          ! Doors/Exits are always vertical
          IF (EMESH_EXITS(N_TMP)%EMESH>0) THEN
             XB(5) = EMESH_XB(5,EMESH_EXITS(N_TMP)%EMESH)
             XB(6) = EMESH_XB(6,EMESH_EXITS(N_TMP)%EMESH)
          END IF
          !
          ! Save the information needed in READ_MESH, etc.
          ! Default, for "no flow field" exits
          EMESH_EXITS(N_TMP)%NEXIT       = N     ! N is the door line index (including all doors)
          EMESH_EXITS(N_TMP)%IOR         = IOR
          EMESH_EXITS(N_TMP)%IS_EXIT     = .FALSE.
          EMESH_EXITS(N_TMP)%DEFINE_MESH = DEFINE_MESH
          EMESH_EXITS(N_TMP)%ID          = TRIM(ID)
          EMESH_EXITS(N_TMP)%XB(:)       = XB(:)
          EMESH_EXITS(N_TMP)%RGB(:)      = RGB(:)
          EMESH_EXITS(N_TMP)%IMESH       = EMESH_EXITS(N_TMP)%EMESH  ! Updated in READ_MESH, the flow mesh index
          EMESH_EXITS(N_TMP)%MAINMESH    = EMESH_EXITS(N_TMP)%EMESH  ! Updated in READ_MESH, the main mesh index
          ! EMESH_EXITS(N_TMP)%EMESH is the main evacuation mesh index: EMESH_ID, EMESH_IJK, EMESH_XB
       END DO READ_DOOR_LOOP_0
27     CONTINUE
       REWIND(LU_INPUT)

       RETURN  ! imode=1: Initialization call to READ_EVAC
    END IF IF_IMODE_1

    ! Lines below are only for imode=2, i.e., after the READ_MESH in read.f90.
    IF (.NOT. ANY(EVACUATION_ONLY)) THEN
       N_EVAC = 0
       IF (MYID==MAX(0,EVAC_PROCESS)) THEN
          IF (ANY(EVACUATION_ONLY)) THEN
             WRITE(MESSAGE,'(A,A,A)') 'ERROR: No main evacuation meshes defined.'
             CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
          END IF
       END IF
       RETURN
    END IF

    NPPS = 30000 ! Number Persons Per Set (dump to a file)
    !
    EVAC_CLOCK = T_BEGIN    ! clock for the CHID_evac.csv file
    EVAC_N_QUANTITIES = 0

    N_RAMP_INI = N_RAMP
    I33 = 0
    ILH = 0

    IF (MYID==MAX(0,EVAC_PROCESS)) THEN
       ALLOCATE(Tsteps(NMESHES),STAT=IZERO)
       CALL ChkMemErr('READ_EVAC','Tsteps',IZERO)
       Tsteps(:) = EVAC_DT_FLOWFIELD
       IF (ABS(TIME_SHRINK_FACTOR-1.0_EB) >= TWO_EPSILON_EB ) THEN
          CALL SHUTDOWN('ERROR: Evac is not ready for TIME_SHRINK_FACTOR',PROCESS_0_ONLY=.FALSE.)
          RETURN
       END IF
    END IF
    !
    ! I_EVAC: 'binary' index:
    ! xxxxx = (0,1)*16 + (0,1)*8 + (0,1)*4 + (0,1)*2 + (0,1)
    ! 0. bit (xxxx1): save flow fields (xxxx0) do not save
    ! 1. bit (xxx1x): save soot + fed  (xxx0x) do not save
    ! 2. bit (xx1xx): read flow fields (xx0xx) do not read
    ! 3. bit (x1xxx): read soot + fed  (x0xxx) do not read
    ! 4. bit (1xxxx): fire calculation (0xxxx) no fire calculation
    ! These are just initialization. Later it is checked if files exists.
    ! If EFF/FED file does not exists, then it is calculated (and saved).
    I_EVAC = 16*0 + 8*0 + 4*0 + 2*0 + 1*0 ! do not save soot,fed files
    IF (.NOT. ALL(EVACUATION_ONLY) ) THEN
       ! Note: If EVACUATION_DRILL=true there are no fire meshes
       ! Note: If NO_EVACUATION=true there are no evacuation meshes
       ! There are fire grids ==> save fed and evac flow fields
       ! Simple chemistry need always REAC line, non-simple chemistry does not need this
       I_EVAC = 16*1 + 8*0 + 4*0 + 2*1 + 1*1
       WRITE(LU_ERR,    '(A,A)')  ' FDS+Evac pressure method : ', TRIM(PRES_METHOD)
       IF (TRIM(PRES_METHOD) .NE. 'FFT') THEN
          WRITE(LU_ERR,'(A)') ' FDS+Evac WARNING: FDS+Evac was developed for FFT solver'
          WRITE(LU_ERR,'(A)') '                   Verify results for other pressure solvers'
          WRITE(LU_ERR,'(A)') '                   The FED file is only saved for later use'
          I_EVAC = 16*1 + 8*0 + 4*0 + 2*1 + 1*0  ! I_EVAC=16 => do only fire calculation, save fed info
       END IF
    ELSE
       ! There are no fire meshes
       WRITE(LU_ERR,'(A,A)')  ' FDS+Evac pressure method : ', TRIM(PRES_METHOD)
       IF (TRIM(PRES_METHOD) .NE. 'FFT') THEN
          WRITE(LU_ERR,    '(A)') ' FDS+Evac ERROR: FDS+Evac was developed for FFT solver'
          WRITE(MESSAGE, '(A,A)') ' ERROR: Evacuation only and pressure solver is: ',TRIM(PRES_METHOD)
          CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
       END IF
       IF (EVACUATION_MC_MODE) THEN
          ! MC-mode: Try to read EFF file if exists on the hard disk
          IF (EVACUATION_DRILL) THEN
             ! There are no fire grids ==> try to read evac flow fields if MC mode
             I_EVAC = 16*0 + 8*0 + 4*1 + 2*0 + 1*0
          ELSE
             ! There are no fire grids ==> try to read fed and evac flow fields if MC mode
             I_EVAC = 16*0 + 8*1 + 4*1 + 2*0 + 1*0
          END IF
       ELSE
          IF (EVACUATION_DRILL) THEN
             ! There are no fire grids ==> recalculate evac flow fields if not MC mode
             I_EVAC = 16*0 + 8*0 + 4*0 + 2*0 + 1*1
          ELSE
             ! There are no fire grids ==> try to read fed and recalculate evac flow fields if not MC mode
             I_EVAC = 16*0 + 8*1 + 4*0 + 2*0 + 1*1
          END IF
       END IF
    END IF
    !
    ! Every human has an identification number, ILABEL_last is the last used number,
    ! so next human will have an identification number, which is ILABEL_last + 1
    ILABEL_last = 0

    CALL COUNT_EVAC_NODES(2) ; IF (STOP_STATUS==SETUP_STOP) RETURN
    N_AVATAR_TYPE = 1 ! The default type index is one and the default name is 'Human'.
    ALLOCATE(TMP_AVATAR_TYPE_NAME(N_ENTRYS+NPC_EVAC),STAT=IZERO) ! etries first then evac lines
    CALL ChkMemErr('READ_EVAC','TMP_AVATAR_TYPE_NAME',IZERO)
    ALLOCATE(TMP_AVATAR_TYPE_INDEX(N_ENTRYS+NPC_EVAC),STAT=IZERO)
    CALL ChkMemErr('READ_EVAC','TMP_AVATAR_TYPE_INDEX',IZERO)
    ALLOCATE(TMP_AVATAR_TYPE_PROP(N_ENTRYS+NPC_EVAC),STAT=IZERO) ! etries first then evac lines
    CALL ChkMemErr('READ_EVAC','TMP_AVATAR_TYPE_PROP',IZERO)
    CALL READ_PERS ; IF (STOP_STATUS==SETUP_STOP) RETURN
    CALL READ_STRS ; IF (STOP_STATUS==SETUP_STOP) RETURN
    CALL READ_EXIT ; IF (STOP_STATUS==SETUP_STOP) RETURN
    CALL READ_DOOR ; IF (STOP_STATUS==SETUP_STOP) RETURN
    CALL READ_CORR ; IF (STOP_STATUS==SETUP_STOP) RETURN
    CALL READ_ENTRIES ; IF (STOP_STATUS==SETUP_STOP) RETURN
    CALL COLLECT_NODE_INFO
    CALL READ_EVAC_LINES ; IF (STOP_STATUS==SETUP_STOP) RETURN
    CALL READ_EVHO ; IF (STOP_STATUS==SETUP_STOP) RETURN
    CALL READ_EVSS ; IF (STOP_STATUS==SETUP_STOP) RETURN
    CALL READ_EDEV ; IF (STOP_STATUS==SETUP_STOP) RETURN

    DEALLOCATE(TMP_AVATAR_TYPE_PROP)
    DEALLOCATE(TMP_AVATAR_TYPE_INDEX)
    DEALLOCATE(TMP_AVATAR_TYPE_NAME)
    IF (MYID /= MAX(0,EVAC_PROCESS)) RETURN

    CALL CHECK_EVAC_NODES ; IF (STOP_STATUS==SETUP_STOP) RETURN

  CONTAINS

    SUBROUTINE COUNT_EVAC_NODES(IMODE)
      IMPLICIT NONE
      !
      ! Passed variables
      INTEGER, INTENT(IN) :: IMODE

      IMODE_IF: IF (IMODE == 1) THEN
         !
         ! Determine total number of PERS lines in the input file
         !
         EVAC_AVATAR_NCOLOR = 0  ! Dimension of Avatar color table
         I_AVATAR_COLOR     = 0  ! Counter for avatar colors
         COLOR_METHOD       = -1 ! Default is standard human colors in Smokeview
         DEAD_RGB           = (/  0,255,255/) ! cyan
         NPC_PERS = 0
         EVAC_DELTA_SEE     = 0.29_EB    ! Eye range (plus minus mid smoke level)
         HUMAN_SMOKE_HEIGHT   = 1.60_EB  ! Nose above floor level
         MAX_HUMANS_DIM = 10000 ! Maximum number of agents in an evacuation mesh
         REWIND(LU_INPUT) ; INPUT_FILE_LINE_NUMBER = 0
         COUNT_PERS_LOOP: DO
            CALL CHECKREAD('PERS',LU_INPUT,IOS)  ; IF (STOP_STATUS==SETUP_STOP) RETURN
            IF (IOS == 1) THEN
               EXIT COUNT_PERS_LOOP
            END IF
            READ(LU_INPUT,NML=PERS,END=221,ERR=222,IOSTAT=IOS)
            NPC_PERS = NPC_PERS + 1
222         IF (IOS > 0) THEN
               WRITE(MESSAGE,'(A,I4,A,I6)') 'ERROR: Problem with PERS ',NPC_PERS+1,', line number',INPUT_FILE_LINE_NUMBER
               CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
            END IF
         END DO COUNT_PERS_LOOP
221      REWIND(LU_INPUT) ; INPUT_FILE_LINE_NUMBER = 0
         IF (COLOR_METHOD == 3) EVAC_AVATAR_NCOLOR = NPC_PERS + 1
         COLOR_METHOD_TMP = COLOR_METHOD
         !
         ! Determine total number of EVAC lines in the input file
         !
         NPC_EVAC = 0
         COUNT_EVAC_LOOP: DO
            NUMBER_INITIAL_PERSONS = 0
            CALL CHECKREAD('EVAC',LU_INPUT,IOS)  ; IF (STOP_STATUS==SETUP_STOP) RETURN
            IF (IOS == 1) THEN
               EXIT COUNT_EVAC_LOOP
            END IF
            READ(LU_INPUT,NML=EVAC,END=219,ERR=220,IOSTAT=IOS)
            NPC_EVAC = NPC_EVAC + 1
            IF (COLOR_METHOD == 0 .AND. NUMBER_INITIAL_PERSONS > 0) THEN
               EVAC_AVATAR_NCOLOR = EVAC_AVATAR_NCOLOR + 1
            END IF
            !
220         IF (IOS > 0) THEN
               WRITE(MESSAGE,'(A,I4,A,I6)') 'ERROR: Problem with EVAC ',NPC_EVAC+1,', line number',INPUT_FILE_LINE_NUMBER
               CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
            END IF
         END DO COUNT_EVAC_LOOP
219      REWIND(LU_INPUT) ; INPUT_FILE_LINE_NUMBER = 0
         !
         ! Determine total number of EXIT lines in the input file
         !
         N_EXITS = 0
         N_CO_EXITS = 0
         COUNT_EXITS_LOOP: DO
            COUNT_ONLY = .FALSE.
            COUNT_DENSITY = .FALSE.
            CALL CHECKREAD('EXIT',LU_INPUT,IOS)  ; IF (STOP_STATUS==SETUP_STOP) RETURN
            IF (IOS == 1) THEN
               EXIT COUNT_EXITS_LOOP
            END IF
            READ(LU_INPUT,NML=EXIT,END=223,ERR=224,IOSTAT=IOS)
            N_EXITS = N_EXITS + 1
            IF (COLOR_METHOD == 4 .AND. .NOT.COUNT_ONLY) THEN
               EVAC_AVATAR_NCOLOR = EVAC_AVATAR_NCOLOR + 1
            END IF
            IF (COUNT_DENSITY) COUNT_ONLY=.TRUE.
            IF (COUNT_ONLY) N_CO_EXITS = N_CO_EXITS + 1
224         IF (IOS > 0) THEN
               WRITE(MESSAGE,'(A,I4,A,I6)') 'ERROR: Problem with EXIT ',N_EXITS+1,', line number',INPUT_FILE_LINE_NUMBER
               CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
            END IF
         END DO COUNT_EXITS_LOOP
223      REWIND(LU_INPUT) ; INPUT_FILE_LINE_NUMBER = 0
         !
         ! Determine total number of DOOR lines in the input file
         !
         N_DOORS = 0
         COUNT_DOORS_LOOP: DO
            CALL CHECKREAD('DOOR',LU_INPUT,IOS)  ; IF (STOP_STATUS==SETUP_STOP) RETURN
            IF (IOS == 1) THEN
               EXIT COUNT_DOORS_LOOP
            END IF
            READ(LU_INPUT,NML=DOOR,END=225,ERR=226,IOSTAT=IOS)
            N_DOORS = N_DOORS + 1
            IF (COLOR_METHOD == 4) THEN
               EVAC_AVATAR_NCOLOR = EVAC_AVATAR_NCOLOR + 1
            END IF
226         IF (IOS > 0) THEN
               WRITE(MESSAGE,'(A,I4,A,I6)') 'ERROR: Problem with DOOR ',N_DOORS + 1,', line number',INPUT_FILE_LINE_NUMBER
               CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
            END IF
         END DO COUNT_DOORS_LOOP
225      REWIND(LU_INPUT) ; INPUT_FILE_LINE_NUMBER = 0
         !
         ! Determine total number of ENTR lines in the input file
         !
         N_ENTRYS = 0
         COUNT_ENTRYS_LOOP: DO
            MAX_FLOW        = 0.0_EB
            MAX_HUMANS_RAMP = 'null'
            CALL CHECKREAD('ENTR',LU_INPUT,IOS)  ; IF (STOP_STATUS==SETUP_STOP) RETURN
            IF (IOS == 1) THEN
               EXIT COUNT_ENTRYS_LOOP
            END IF
            READ(LU_INPUT,NML=ENTR,END=227,ERR=228,IOSTAT=IOS)
            N_ENTRYS = N_ENTRYS + 1
            IF (COLOR_METHOD == 0 .AND. (MAX_FLOW >= TWO_EPSILON_EB .OR. Trim(MAX_HUMANS_RAMP)/='null')) THEN
               EVAC_AVATAR_NCOLOR = EVAC_AVATAR_NCOLOR + 1
            END IF
228         IF (IOS > 0) THEN
               WRITE(MESSAGE,'(A,I4,A,I6)') 'ERROR: Problem with ENTR ',N_ENTRYS + 1,', line number',INPUT_FILE_LINE_NUMBER
               CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
            END IF
         END DO COUNT_ENTRYS_LOOP
227      REWIND(LU_INPUT) ; INPUT_FILE_LINE_NUMBER = 0
         !
         ! Determine total number of CORR lines in the input file
         !
         N_CORRS = 0
         COUNT_CORRS_LOOP: DO
            CALL CHECKREAD('CORR',LU_INPUT,IOS)  ; IF (STOP_STATUS==SETUP_STOP) RETURN
            IF (IOS == 1) THEN
               EXIT COUNT_CORRS_LOOP
            END IF
            READ(LU_INPUT,NML=CORR,END=229,ERR=230,IOSTAT=IOS)
            N_CORRS = N_CORRS + 1
230         IF (IOS > 0) THEN
               WRITE(MESSAGE,'(A,I4,A,I6)') 'ERROR: Problem with CORR ',N_CORRS + 1,', line number',INPUT_FILE_LINE_NUMBER
               CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
            END IF
         END DO COUNT_CORRS_LOOP
229      REWIND(LU_INPUT) ; INPUT_FILE_LINE_NUMBER = 0
         !
         ! Determine total number of EVHO lines in the input file
         !
         N_HOLES = 0
         COUNT_EVHO_LOOP: DO
            CALL CHECKREAD('EVHO',LU_INPUT,IOS)  ; IF (STOP_STATUS==SETUP_STOP) RETURN
            IF (IOS == 1) THEN
               EXIT COUNT_EVHO_LOOP
            END IF
            READ(LU_INPUT,NML=EVHO,END=231,ERR=232,IOSTAT=IOS)
            N_HOLES = N_HOLES + 1
232         IF (IOS > 0) THEN
               WRITE(MESSAGE,'(A,I4,A,I6)') 'ERROR: Problem with EVHO ',N_HOLES + 1,', line number',INPUT_FILE_LINE_NUMBER
               CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
            END IF
         END DO COUNT_EVHO_LOOP
231      REWIND(LU_INPUT) ; INPUT_FILE_LINE_NUMBER = 0
         !
         ! Determine total number of EVSS lines in the input file
         !
         N_SSTANDS = 0
         COUNT_EVSS_LOOP: DO
            CALL CHECKREAD('EVSS',LU_INPUT,IOS)  ; IF (STOP_STATUS==SETUP_STOP) RETURN
            IF (IOS == 1) THEN
               EXIT COUNT_EVSS_LOOP
            END IF
            READ(LU_INPUT,NML=EVSS,END=233,ERR=234,IOSTAT=IOS)
            N_SSTANDS = N_SSTANDS + 1
234         IF (IOS > 0) THEN
               WRITE(MESSAGE,'(A,I4,A,I6)') 'ERROR: Problem with EVSS ',N_SSTANDS + 1,', line number',INPUT_FILE_LINE_NUMBER
               CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
            END IF
         END DO COUNT_EVSS_LOOP
233      REWIND(LU_INPUT) ; INPUT_FILE_LINE_NUMBER = 0
         !
         ! Determine total number of STRS lines in the input file
         !
         N_STRS = 0
         COUNT_STRS_LOOP: DO
            CALL CHECKREAD('STRS',LU_INPUT,IOS)  ; IF (STOP_STATUS==SETUP_STOP) RETURN
            IF (IOS == 1) THEN
               EXIT COUNT_STRS_LOOP
            END IF
            READ(LU_INPUT,NML=STRS,END=235,ERR=236,IOSTAT=IOS)
            N_STRS = N_STRS + 1
236         IF (IOS > 0) THEN
               WRITE(MESSAGE,'(A,I4,A,I6)') 'ERROR: Problem with STRS ',N_STRS + 1,', line number',INPUT_FILE_LINE_NUMBER
               CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
            END IF
         END DO COUNT_STRS_LOOP
235      REWIND(LU_INPUT) ; INPUT_FILE_LINE_NUMBER = 0
         !
         ! Determine total number of EDEV lines in the input file
         !
         N_EDEV = 0
         COUNT_EDEV_LOOP: DO
            CALL CHECKREAD('EDEV',LU_INPUT,IOS)  ; IF (STOP_STATUS==SETUP_STOP) RETURN
            IF (IOS == 1) THEN
               EXIT COUNT_EDEV_LOOP
            END IF
            READ(LU_INPUT,NML=EDEV,END=237,ERR=238,IOSTAT=IOS)
            N_EDEV = N_EDEV + 1
238         IF (IOS > 0) THEN
               WRITE(MESSAGE,'(A,I4,A,I6)') 'ERROR: Problem with EDEV ',N_EDEV + 1,', line number',INPUT_FILE_LINE_NUMBER
               CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
            END IF
         END DO COUNT_EDEV_LOOP
237      REWIND(LU_INPUT) ; INPUT_FILE_LINE_NUMBER = 0

         SELECT CASE (COLOR_METHOD)
         CASE (-1)
            EVAC_AVATAR_NCOLOR = 1
         CASE (0,3,4)
            EVAC_AVATAR_NCOLOR = EVAC_AVATAR_NCOLOR + 1
         CASE (1,2,5)
            EVAC_AVATAR_NCOLOR = 7
         CASE Default
            EVAC_AVATAR_NCOLOR = 1
         END SELECT

         ! Allocate avatar color array for Smokeview file write
         EVAC_AVATAR_NCOLOR = MAX(1,EVAC_AVATAR_NCOLOR)
         ALLOCATE(EVAC_AVATAR_RGB(3,EVAC_AVATAR_NCOLOR),STAT=IZERO)
         CALL ChkMemErr('READ_EVAC','EVAC_AVATAR_RGB',IZERO)
         EVAC_AVATAR_RGB = 0

         SELECT CASE (COLOR_METHOD)
         CASE (-1)
            EVAC_AVATAR_RGB(1:3,1) = (/ 39, 64,139/)  ! ROYAL BLUE 4
         CASE (3)
            EVAC_AVATAR_RGB(1:3,1) = (/ 39, 64,139/)  ! ROYAL BLUE 4
            EVAC_AVATAR_RGB(1:3,EVAC_AVATAR_NCOLOR) = DEAD_RGB
         CASE (0,4)
            EVAC_AVATAR_RGB(1:3,EVAC_AVATAR_NCOLOR) = DEAD_RGB
         CASE (1,2,5)
            EVAC_AVATAR_RGB(1:3,1) = (/  0,  0,  0/)  ! black
            EVAC_AVATAR_RGB(1:3,2) = (/255,255,  0/)  ! yellow
            EVAC_AVATAR_RGB(1:3,3) = (/  0,  0,255/)  ! blue
            EVAC_AVATAR_RGB(1:3,4) = (/255,  0,  0/)  ! red
            EVAC_AVATAR_RGB(1:3,5) = (/  0,255,  0/)  ! green
            EVAC_AVATAR_RGB(1:3,6) = (/255,  0,255/)  ! magenta
            EVAC_AVATAR_RGB(1:3,7) = DEAD_RGB
         CASE Default
            EVAC_AVATAR_RGB(1:3,1) = (/ 39, 64,139/)  ! ROYAL BLUE 4
         END SELECT

      ELSE ! IMODE_IF = 2
         ALLOCATE(EMESH_INDEX(MAX(1,NMESHES)), STAT=IZERO)
         CALL ChkMemErr('READ_EVAC','EMESH_INDEX',IZERO)
         EMESH_INDEX = 0
         N_EGRIDS = 0
         DO N = 1, NMESHES
            IF (EVACUATION_ONLY(N) .AND. EVACUATION_SKIP(N)) THEN
               N_EGRIDS = N_EGRIDS + 1
               EMESH_INDEX(N) = N_EGRIDS
            END IF
         END  DO

         ! Allocate quantities for EVAC, PERS, EXIT types and initialize some variables
         !
         EVAC_PROC_IF: IF (MYID==MAX(0,EVAC_PROCESS)) THEN
            IF (NPC_EVAC > 0 ) THEN
               ALLOCATE(EVAC_EVACS(NPC_EVAC),STAT=IZERO)
               CALL ChkMemErr('READ','EVAC_EVACS',IZERO)
               EVAC_EVACS(:)%I_EVAC_THIS_MESH = 1
            END IF

            IF (N_HOLES > 0 ) THEN
               ALLOCATE(EVAC_HOLES(N_HOLES),STAT=IZERO)
               CALL ChkMemErr('READ','EVAC_HOLES',IZERO)
               EVAC_HOLES(:)%I_EVHO_THIS_MESH = 1
            END IF

            IF (N_SSTANDS > 0 ) THEN
               ALLOCATE(EVAC_SSTANDS(N_SSTANDS),STAT=IZERO)
               CALL ChkMemErr('READ','EVAC_SSTANDS',IZERO)
            END IF

            IF (N_EXITS > 0 ) THEN
               ALLOCATE(EVAC_EXITS(N_EXITS),STAT=IZERO)
               CALL ChkMemErr('READ','EVAC_EXITS',IZERO)
            END IF

            IF (N_DOORS > 0 ) THEN
               ALLOCATE(EVAC_DOORS(N_DOORS),STAT=IZERO)
               CALL ChkMemErr('READ','EVAC_DOORS',IZERO)
            END IF

            IF (N_ENTRYS > 0 ) THEN
               ALLOCATE(EVAC_ENTRYS(N_ENTRYS),STAT=IZERO)
               CALL ChkMemErr('READ','EVAC_ENTRYS',IZERO)
            END IF

            IF (N_CORRS > 0 ) THEN
               ALLOCATE(EVAC_CORRS(N_CORRS),STAT=IZERO)
               CALL ChkMemErr('READ','EVAC_CORRS',IZERO)
            END IF

            IF (N_STRS > 0 ) THEN
               ALLOCATE(EVAC_STRS(N_STRS),STAT=IZERO)
               CALL ChkMemErr('READ','EVAC_STRS',IZERO)
            END IF

            IF (N_EDEV > 0 ) THEN
               ALLOCATE(EVAC_EDEV(N_EDEV),STAT=IZERO)
               CALL ChkMemErr('READ','EVAC_EDEV',IZERO)
            END IF

            ALLOCATE(EVAC_PERSON_CLASSES(0:NPC_PERS),STAT=IZERO)
            CALL ChkMemErr('READ','EVAC_PERSON_CLASSES',IZERO)

            IF (NPC_EVAC > 0) THEN
               EVAC_EVACS(1:NPC_EVAC)%GRID_NAME   = 'null'
               EVAC_EVACS(1:NPC_EVAC)%CLASS_NAME  = 'null'
               EVAC_EVACS(1:NPC_EVAC)%IMESH       = 0
               EVAC_EVACS(1:NPC_EVAC)%ID          = 'null'
            END IF
            IF (N_HOLES > 0) THEN
               EVAC_HOLES(1:N_HOLES)%GRID_NAME   = 'null'
               EVAC_HOLES(1:N_HOLES)%PERS_ID     = 'null'
               EVAC_HOLES(1:N_HOLES)%EVAC_ID     = 'null'
               EVAC_HOLES(1:N_HOLES)%IMESH       = 0
               EVAC_HOLES(1:N_HOLES)%ID          = 'null'
            END IF

            EVAC_PERSON_CLASSES(0:NPC_PERS)%ID = 'null'

            IF (N_EXITS > 0) THEN
               EVAC_EXITS(1:N_EXITS)%ID        = 'null'
               EVAC_EXITS(1:N_EXITS)%TO_NODE   = 'null'
               EVAC_EXITS(1:N_EXITS)%GRID_NAME = 'null'
               EVAC_EXITS(1:N_EXITS)%IMESH     = 0
            END IF

            IF (N_DOORS > 0) THEN
               EVAC_DOORS(1:N_DOORS)%ID        = 'null'
               EVAC_DOORS(1:N_DOORS)%TO_NODE   = 'null'
               EVAC_DOORS(1:N_DOORS)%GRID_NAME = 'null'
               EVAC_DOORS(1:N_DOORS)%IMESH     = 0
               EVAC_DOORS(1:N_DOORS)%IMESH2    = 0
            END IF

            IF (N_CORRS > 0) THEN
               EVAC_CORRS(1:N_CORRS)%ID        = 'null'
               EVAC_CORRS(1:N_CORRS)%TO_NODE   = 'null'
               EVAC_CORRS(1:N_CORRS)%GRID_NAME = 'null'
               EVAC_CORRS(1:N_CORRS)%IMESH     = 0
               EVAC_CORRS(1:N_CORRS)%IMESH2    = 0
            END IF

            IF (N_ENTRYS > 0) THEN
               EVAC_ENTRYS(1:N_ENTRYS)%ID          = 'null'
               EVAC_ENTRYS(1:N_ENTRYS)%TO_NODE     = 'null'
               EVAC_ENTRYS(1:N_ENTRYS)%GRID_NAME   = 'null'
               EVAC_ENTRYS(1:N_ENTRYS)%CLASS_NAME  = 'null'
               EVAC_ENTRYS(1:N_ENTRYS)%IMESH       = 0
            END IF

         END IF EVAC_PROC_IF

         EVAC_PROC_IF_2: IF (MYID==MAX(0,EVAC_PROCESS)) THEN
            N_NODES = N_ENTRYS + N_EXITS + N_DOORS + N_CORRS + N_EGRIDS + N_STRS
            IF (N_NODES > 0) THEN
               ALLOCATE(EVAC_Node_List(1:N_NODES),STAT=IZERO)
               CALL ChkMemErr('READ','EVAC_NODE_LIST',IZERO)
            END IF
         END IF EVAC_PROC_IF_2

      END IF IMODE_IF

    END SUBROUTINE COUNT_EVAC_NODES

    SUBROUTINE READ_PERS
      IMPLICIT NONE
      !
      ! Local variables
      TYPE (EVAC_PERS_TYPE), POINTER :: PCP=>NULL()
      !
      ! Next parameters are same for all humans. the last values read in from 'pers' lines are valid.
      FAC_A_WALL  = 1.0_EB
      FAC_B_WALL  = 0.5_EB
      LAMBDA_WALL = 0.2_EB
      NOISEME     = 0.0_EB
      NOISETH     = 0.01_EB
      NOISECM     = 3.0_EB
      NOT_RANDOM  = .FALSE.
      EVAC_FDS6   = .TRUE.  ! The default and only one, fds5 style is not supported anymore
      I_FRIC_SW   = 1
      V_MAX             = 20.0_EB ! m/s
      V_ANGULAR_MAX     = 8.0_EB  ! rps
      V_ANGULAR         = 2.0_EB  ! rps
      FC_DAMPING        = 500.0_EB ! N/(m/s)
      ! Z_smoke = XB_z - EVACUATION_Z_OFFSET(NM) + HUMAN_SMOKE_HEIGHT, i.e. position of the nose/eyes
      ! above the floor.  The smoke and gas densities are taken from this level (FED calculation and
      ! visible doors etc.)

      GROUP_EFF         = 0.0_EB
      RADIUS_COMPLETE_0 = 0.2_EB
      RADIUS_COMPLETE_1 = 0.5_EB

      ! TDET_SMOKE_DENS: Smoke is detected, when its density is larger than, e.g. 1 mg/m3
      !                  Default is no detection due to smoke.
      TDET_SMOKE_DENS = -999.9_EB  !
      FED_DOOR_CRIT   = -100.0_EB ! Which doors are 'smoke free' Evac >= 2.2.2
      GROUP_DENS      = 0.0_EB
      SMOKE_MIN_SPEED = -99999._EB
      SMOKE_MIN_SPEED_FACTOR = -99999._EB
      SMOKE_MIN_SPEED_VISIBILITY = 0.0_EB
      SMOKE_KS_SPEED_FUNCTION = 1 ! Linear function c(Ks) = (1 + beta*Ks/alpha)
      FED_ACTIVITY = 2  ! 2:light work, 1: at rest, 3: heavy work
      TAU_CHANGE_DOOR = 1.0_EB
      DENS_INIT           = -99999.0_EB
      MAX_INITIAL_OVERLAP = -99999.0_EB
      EVAC_DT_MAX     = 0.01_EB
      EVAC_DT_MIN     = 0.001_EB
      !Issue1547: Some time is needed, because the agents are initially at rest, so the initial
      !Issue1547: acceleration time is not counted. Time units are seconds, so 60 seconds.
      TIME_INIT_NERVOUSNESS    = 60.0_EB

      CROWBAR_DT_READ = 1.0_EB

      ! Linear speed vs smoke density function, default is the Lund fit.
      ! c(Ks) = ( 1 + (BETA*Ks)/ALPHA )  [0,1] interval, c(Ks=0)=1, c(Ks=inf)=0
      SMOKE_SPEED_ALPHA = 0.706_EB   ! Lund 2003, report 3126 (Frantzich & Nilsson)
      SMOKE_SPEED_BETA  = -0.057_EB  ! Lund 2003, report 3126 (Frantzich & Nilsson)

      ! Next parameters are for the counterflow (CF)
      ! Evac 2.2.0: Counterflow treatment is the default
      TAU_CHANGE_V0   = 0.1_EB   !CF: How often direction is updated?
      THETA_SECTOR    = -40.0_EB !CF: The angle of the first sector
      CONST_DF        = 2.0_EB   !CF: prefer agents going in the same direction
      FAC_DF          = 1.0_EB   !CF: prefer agents going in the same direction
      CONST_CF        = 1.0_EB   !CF: dislike agents going in the opposite direction
      FAC_CF          = 2.0_EB   !CF: dislike agents going in the opposite direction
      FAC_1_WALL      = 5.0_EB   !CF: direction is towards a wall
      FAC_2_WALL      = 10.0_EB  !CF: direction leads too close to a wall
      FAC_V0_DIR      = 1.0_EB   !CF: Prefer right and straight ahead (if CF)
      FAC_V0_NOCF     = 1.0_EB   !CF: prefer v0, if no counterflow
      FAC_NOCF        = 2.0_EB   !CF: prefer v0, if no counterflow
      CF_MIN_A        = 0.5_EB   !CF: decrease social force
      CF_MIN_B        = 0.3_EB   !CF: decrease social force range
      CF_FAC_A_WALL   = 1.0_EB   !CF: decrease social force
      CF_MIN_TAU      = 0.10_EB  !CF: increase motive force
      CF_MIN_TAU_INER = 0.05_EB  !CF: increase motive force, rotation
      CF_FAC_TAUS     = 0.25_EB  !CF: increase motive force, trans+rot

      ! Evac 2.2.0: Use the waiting time in the door selection is the default
      FAC_DOOR_QUEUE  = 1.3_EB  ! Door selection algorithm: persons/m/s
      FAC_DOOR_ALPHA  = 0.5_EB  ! Door selection algorithm: alpha*t_walk + (1-alpha)*t_queue
      FAC_DOOR_WAIT   = 0.9_EB  ! Door selection algorithm: patience factor
      FAC_DOOR_OLD    = 0.1_EB  ! The present door is considered "smoke free" longer than others
      FAC_DOOR_OLD2   = 0.9_EB  ! The present door is considered "not too much smoke" longer than others
      NASH_CLOSE_ENOUGH = 0.05_EB ! If less than 5 percent changes targer door in the initialization => converged

      R_HERDING       =  5.0_EB  ! Herding agents: How far to look, the radius
      W0_HERDING      = 10.0_EB  ! Herding agents: Weight at the distance r=0 (linear function)
      WR_HERDING      =  0.0_EB  ! Herding agents: Weight at the distance r=R_HERDING
      DOT_HERDING     = -0.2_EB  ! Herding agents: dot product parameter (do not count agents that are heading towards)
      HERDING_TAU_FACTOR=0.1 ! If no herding agents close, how probable it is to use the information
                             ! from the agents further away, a probability 0.0-1.0 range.
                             ! If =1.0 then this is done everytime, when the agents chooses a door
                             ! and this is done every TAU_CHANGE_DOOR second on the average.
      I_HERDING_TYPE  = 0    ! Herding agents: >1 do not move if no door (0 default ffield)
      HERDING_TPRE_PROB = 0.0_EB
      ! Hawk - Dove game: Just an academic exercise at this moment, works only for a simple geometries.
      ! Hawk - Dove game parameters: If C_HAWK < 0 then no game is played, i.e., normal FDS+Evac
      ! Hawk - Dove game: Only rational agents play (AGENT_TYPE=1 on EVAC/ENTR namelist)
      ! U_i(Hawk) = p_i*c/(2*Delta*alpha*T_i) - (1-p_i), U_i(Dove) = p_i
      C_HAWK       = -2.0_EB ! Ctilde_Hawk = c/(2*alpha*Delta): 0-1 always hawks, >1 more and more doves
      ALPHA_HAWK   =  1.0_EB !
      DELTA_HAWK   =  1.0_EB !
      EPSILON_HAWK =  0.0_EB ! 2*Alpha*Delta*T_i < Epsilon then dove, i.e., hawk should be at least better than epsilon
      THETA_HAWK   =  1.0_EB ! random force factor for hawks: Gaussian theta = theta*theta_hawk
      R_HAWK_DOVE  =  0.5_EB ! cut-off radius, play against the agents which are inside the skin-skin distance
      A_HAWK_T_FACTOR=0.0_EB ! linear growth factor of the alpha parameter (the threath parameter)
                             ! alpha = alpha_0*(1+A_HAWK_T_FACTOR*(t-t_start))
      A_HAWK_T_START=T_BEGIN ! linear growth factor of the alpha parameter (the threath parameter)
                             ! alpha is constant upto this time, then a=a0*(1+a_fac*(t-t_start))
      T_STOP_HD_GAME = T_END ! Do not play after the given time
      ! New parameters for the Hawk-Dove game:
      !    C_HAWK = T_0/(T_i-T_ASET+T_0), when T_i > T_ASET-T_0
      !    C_HAWK = very large number,    when T_i =< T_ASET-T_0
      !    If EPSILON_HAWK > 0 then
      T_ASET_HAWK      = -1.0_EB ! seconds
      T_0_HAWK         =  1.0_EB ! seconds, how much before T_ASET cost function starts to increase
      T_ASET_TFAC_HAWK = 0.0_EB  ! seconds per seconds, how much T_ASET changes

      ! Paramters of the fall down algorithm
      F_MIN_FALL           = 2000.0_EB
      F_MAX_FALL           = 9000.0_EB
      D_OVERLAP_FALL       = 0.05_EB
      TAU_FALL_DOWN        = -1.0_EB
      A_FAC_FALLEN         = 1.0_EB
      PROB_FALL_DOWN       = -1.0_EB

      OUTPUT_SPEED         = .FALSE.
      OUTPUT_DENSITY       = .FALSE.
      OUTPUT_MOTIVE_FORCE  = .FALSE.
      OUTPUT_MOTIVE_ANGLE  = .FALSE.
      OUTPUT_FED           = .FALSE.
      OUTPUT_OMEGA         = .FALSE.
      OUTPUT_ANGLE         = .FALSE.
      OUTPUT_ACCELERATION  = .FALSE.
      OUTPUT_CONTACT_FORCE = .FALSE.
      OUTPUT_TOTAL_FORCE   = .FALSE.
      OUTPUT_NERVOUSNESS   = .FALSE.
      DISCARD_SMOKE_INFO   = .FALSE. ! Set true in input file for debug purposes only
      !Issue1547: This subroutine reads the PERS namelists. All the namelist entries should
      !Issue1547: have some default values that are used if no value is given in the input.
      !Issue1547: This way you will not have any uninitialized variables in Fortran.
      DEAD_COLOR           = 'null'
      !
      ! Read the PERS lines (no read for default n=0 case)
      !
      READ_PERS_LOOP: DO N = 0, NPC_PERS
         !
         ID           = 'null'
         RGB          = -1
         COLOR        = 'null'
         AVATAR_RGB   = -1
         AVATAR_COLOR = 'null'
         ! Default: No distributions for human properties
         DEFAULT_PROPERTIES = 'null'
         DIAMETER_DIST = -1
         VELOCITY_DIST = -1
         TAU_EVAC_DIST = -1
         PRE_EVAC_DIST = 0
         DET_EVAC_DIST = 0
         VEL_PARA = 0.0_EB
         DIA_PARA = 0.0_EB
         PRE_PARA = 0.0_EB
         DET_PARA = 0.0_EB
         TAU_PARA = 0.0_EB
         VEL_PARA2 = 0.0_EB
         DIA_PARA2 = 0.0_EB
         PRE_PARA2 = 0.0_EB
         DET_PARA2 = 0.0_EB
         TAU_PARA2 = 0.0_EB
         VEL_LOW = 0.0_EB
         DIA_LOW = 0.0_EB
         PRE_LOW = 0.0_EB
         DET_LOW = T_BEGIN
         TAU_LOW = 0.0_EB
         VEL_HIGH = 999.0_EB
         DIA_HIGH = 999.0_EB
         PRE_HIGH = HUGE(PRE_HIGH)
         DET_HIGH = HUGE(PRE_HIGH)
         TAU_HIGH = 999.0_EB
         ! Default values for persons
         VEL_MEAN = 1.25_EB
         DIA_MEAN = -10.0_EB
         PRE_MEAN = 10.0_EB
         DET_MEAN = T_BEGIN
         TAU_MEAN = 1.0_EB
         FCONST_A = 2000.0_EB
         FCONST_B = 0.08_EB
         ! Evac 2.2.0: The default social force anisotropy parameter lambda changed
         L_NON_SP = 0.3_EB ! Evac 2.2.0
         C_YOUNG  = 120000.0_EB
         GAMMA    = 16000.0_EB
         KAPPA    = 240000.0_EB
         ! Rotational freedom constants
         D_TORSO_MEAN = 0.30_EB
         D_SHOULDER_MEAN = 0.19_EB
         TAU_ROT   = 0.2_EB
         M_INERTIA = -4.0_EB ! Default male mass of inertia
         MASS_OF_AGENT = -80.0_EB ! Default male mass

         ! If not given on PERS line, use those given on EVSS lines
         FAC_V0_UP   = -1.0_EB
         FAC_V0_DOWN = -1.0_EB
         FAC_V0_HORI = -1.0_EB

         !Issue1547: Nervousness dependet social force parameter (desired velocity is changed)
         !Issue1547: This is read for each PERS line separately. Above the OUTPUT_NERVOUSNESS
         !Issue1547: is outside the read PERS loop, so just the last value read in will be
         !Issue1547: used. The read in value of MAXIMUM_V0_FACTOR is saved in the corresponding
         !Issue1547: structured data type (later PCP => EVAC_PERSON_CLASSES(N) is a pointer).
         MAXIMUM_V0_FACTOR    = -1.0_EB ! If less than zero, do not use nervousness algorithm
         !
         ! No read for default values
         IF (N > 0) THEN
            CALL CHECKREAD('PERS',LU_INPUT,IOS)  ; IF (STOP_STATUS==SETUP_STOP) RETURN
            IF (IOS == 1) THEN
               EXIT READ_PERS_LOOP
            END IF
            READ(LU_INPUT,PERS,END=24,IOSTAT=IOS)

            IF (.NOT.EVAC_FDS6) THEN
               WRITE(MESSAGE,'(A,A,A)') 'ERROR: PERS ',TRIM(ID),' EVAC_FDS6=F is not anymore supported.'
               CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
            END IF
            ! IMO MSC.1/Circ. 1238, ANNEX 2, Table 3.4, uniform distributions for velocity
            SELECT CASE (TRIM(DEFAULT_PROPERTIES))
            CASE ('IMO_Male<30','imo_male<30','IMO_MALE<30','IMO_MaleCrew','imo_malecrew','IMO_MALECREW')
               IF (VELOCITY_DIST < 0) THEN
                  VELOCITY_DIST = 1
                  VEL_LOW  = 1.11_EB
                  VEL_HIGH = 1.85_EB
                  VEL_MEAN = 0.5_EB*(VEL_LOW+VEL_HIGH)
               END IF
               ! Stairs down,up 0.76-1.26 0.5-0.84
               !  down about 0.683, up 0.452
               FAC_V0_UP   = 0.452_EB
               FAC_V0_DOWN = 0.683_EB
               FAC_V0_HORI = 1.0_EB
               DEFAULT_PROPERTIES='Male' ! Male body size etc.
            CASE ('IMO_Male30-50','imo_male30-50','IMO_MALE30-50')
               IF (VELOCITY_DIST < 0) THEN
                  VELOCITY_DIST = 1
                  VEL_LOW  = 0.97_EB
                  VEL_HIGH = 1.62_EB
                  VEL_MEAN = 0.5_EB*(VEL_LOW+VEL_HIGH)
               END IF
               ! Stairs down,up 0.64-1.07 0.47-0.79
               !  down about 0.660, up 0.485
               FAC_V0_UP   = 0.485_EB
               FAC_V0_DOWN = 0.660_EB
               FAC_V0_HORI = 1.0_EB
               DEFAULT_PROPERTIES='Male' ! Male body size etc.
            CASE ('IMO_Male>50','imo_male>50','IMO_MALE>50')
               IF (VELOCITY_DIST < 0) THEN
                  VELOCITY_DIST = 1
                  VEL_LOW  = 0.84_EB
                  VEL_HIGH = 1.40_EB
                  VEL_MEAN = 0.5_EB*(VEL_LOW+VEL_HIGH)
               END IF
               ! Stairs down,up 0.5-0.84 0.38-0.64
               !  down about 0.598, up 0.454
               FAC_V0_UP   = 0.454_EB
               FAC_V0_DOWN = 0.598_EB
               FAC_V0_HORI = 1.0_EB
               DEFAULT_PROPERTIES='Male' ! Male body size etc.
            CASE ('IMO_Male>50MI1','imo_male>50mi1','IMO_MALE>50MI1')
               IF (VELOCITY_DIST < 0) THEN
                  VELOCITY_DIST = 1
                  VEL_LOW  = 0.64_EB
                  VEL_HIGH = 1.06_EB
                  VEL_MEAN = 0.5_EB*(VEL_LOW+VEL_HIGH)
               END IF
               ! Stairs down,up 0.38-0.64 0.29-0.49
               !  down about 0.598, up 0.458
               FAC_V0_UP   = 0.458_EB
               FAC_V0_DOWN = 0.598_EB
               FAC_V0_HORI = 1.0_EB
               DEFAULT_PROPERTIES='Male' ! Male body size etc.
            CASE ('IMO_Male>50MI2','imo_male>50mi2','IMO_MALE>50MI2')
               IF (VELOCITY_DIST < 0) THEN
                  VELOCITY_DIST = 1
                  VEL_LOW  = 0.55_EB
                  VEL_HIGH = 0.91_EB
                  VEL_MEAN = 0.5_EB*(VEL_LOW+VEL_HIGH)
               END IF
               ! Stairs down,up 0.33-0.55 0.25-0.41
               !  down about 0.602, up 0.452
               FAC_V0_UP   = 0.452_EB
               FAC_V0_DOWN = 0.602_EB
               FAC_V0_HORI = 1.0_EB
               DEFAULT_PROPERTIES='Male' ! Male body size etc.
            CASE ('IMO_Female<30','imo_female<30','IMO_FEMALE<30','IMO_FemaleCrew','imo_femalecrew','IMO_FEMALECREW')
               IF (VELOCITY_DIST < 0) THEN
                  VELOCITY_DIST = 1
                  VEL_LOW  = 0.93_EB
                  VEL_HIGH = 1.55_EB
                  VEL_MEAN = 0.5_EB*(VEL_LOW+VEL_HIGH)
               END IF
               ! Stairs down,up 0.56-0.94 0.47-0.79
               !  down about 0.604, up 0.507
               FAC_V0_UP   = 0.507_EB
               FAC_V0_DOWN = 0.604_EB
               FAC_V0_HORI = 1.0_EB
               DEFAULT_PROPERTIES='Female' ! Female body size etc.
            CASE ('IMO_Female30-50','imo_female30-50','IMO_FEMALE30-50')
               IF (VELOCITY_DIST < 0) THEN
                  VELOCITY_DIST = 1
                  VEL_LOW  = 0.71_EB
                  VEL_HIGH = 1.19_EB
                  VEL_MEAN = 0.5_EB*(VEL_LOW+VEL_HIGH)
               END IF
               ! Stairs down,up 0.49-0.81 0.44-0.74
               !  down about 0.685, up 0.621
               FAC_V0_UP   = 0.621_EB
               FAC_V0_DOWN = 0.685_EB
               FAC_V0_HORI = 1.0_EB
               DEFAULT_PROPERTIES='Female' ! Female body size etc.
            CASE ('IMO_Female>50','imo_female>50','IMO_FEMALE>50')
               IF (VELOCITY_DIST < 0) THEN
                  VELOCITY_DIST = 1
                  VEL_LOW  = 0.56_EB
                  VEL_HIGH = 0.94_EB
                  VEL_MEAN = 0.5_EB*(VEL_LOW+VEL_HIGH)
               END IF
               ! Stairs down,up 0.45-0.75 0.37-0.61
               !  down about 0.800, up 0.654
               FAC_V0_UP   = 0.654_EB
               FAC_V0_DOWN = 0.800_EB
               FAC_V0_HORI = 1.0_EB
               DEFAULT_PROPERTIES='Female' ! Female body size etc.
            CASE ('IMO_Female>50MI1','imo_female>50mi1','IMO_FEMALE>50MI1')
               IF (VELOCITY_DIST < 0) THEN
                  VELOCITY_DIST = 1
                  VEL_LOW  = 0.43_EB
                  VEL_HIGH = 0.71_EB
                  VEL_MEAN = 0.5_EB*(VEL_LOW+VEL_HIGH)
               END IF
               ! Stairs down,up 0.34-0.56 0.28-0.46
               !  down about 0.790, up 0.649
               FAC_V0_UP   = 0.649_EB
               FAC_V0_DOWN = 0.790_EB
               FAC_V0_HORI = 1.0_EB
               DEFAULT_PROPERTIES='Female' ! Female body size etc.
            CASE ('IMO_Female>50MI2','imo_female>50mi2','IMO_FEMALE>50MI2')
               IF (VELOCITY_DIST < 0) THEN
                  VELOCITY_DIST = 1
                  VEL_LOW  = 0.37_EB
                  VEL_HIGH = 0.61_EB
                  VEL_MEAN = 0.5_EB*(VEL_LOW+VEL_HIGH)
               END IF
               ! Stairs down,up 0.29-0.49 0.23-0.39
               !  down about 0.790, up 0.631
               FAC_V0_UP   = 0.631_EB
               FAC_V0_DOWN = 0.790_EB
               FAC_V0_HORI = 1.0_EB
               DEFAULT_PROPERTIES='Female' ! Female body size etc.
            CASE Default
            END SELECT

            ! Check if some default human group is given.
            SELECT CASE (TRIM(DEFAULT_PROPERTIES))
            CASE ('Adult','adult','ADULT')
               IF (VELOCITY_DIST < 0) THEN
                  VELOCITY_DIST = 1
                  VEL_MEAN = 1.25_EB
                  VEL_PARA = 0.30_EB
                  VEL_LOW  = 0.95_EB
                  VEL_HIGH = 1.55_EB
               END IF
               IF (DIAMETER_DIST < 0) THEN
                  DIAMETER_DIST = 1
                  DIA_MEAN = 0.51_EB
                  DIA_PARA = 0.07_EB
                  DIA_LOW  = 0.44_EB
                  DIA_HIGH = 0.58_EB
                  D_TORSO_MEAN    = 0.30_EB
                  D_SHOULDER_MEAN = 0.19_EB
               END IF
               IF (TAU_EVAC_DIST < 0) THEN
                  TAU_EVAC_DIST = 1
                  TAU_MEAN = 1.00_EB
                  TAU_PARA = 0.10_EB
                  TAU_LOW  = 0.80_EB
                  TAU_HIGH = 1.20_EB
               END IF
            CASE ('Male','male','MALE')
               IF (VELOCITY_DIST < 0) THEN
                  VELOCITY_DIST = 1
                  VEL_MEAN = 1.35_EB
                  VEL_PARA = 0.20_EB
                  VEL_LOW  = 1.15_EB
                  VEL_HIGH = 1.55_EB
               END IF
               IF (DIAMETER_DIST < 0) THEN
                  DIAMETER_DIST = 1
                  DIA_MEAN = 0.54_EB
                  DIA_PARA = 0.04_EB
                  DIA_LOW  = 0.50_EB
                  DIA_HIGH = 0.58_EB
                  D_TORSO_MEAN    = 0.32_EB
                  D_SHOULDER_MEAN = 0.20_EB
               END IF
               IF (TAU_EVAC_DIST < 0) THEN
                  TAU_EVAC_DIST = 1
                  TAU_MEAN = 1.00_EB
                  TAU_PARA = 0.10_EB
                  TAU_LOW  = 0.80_EB
                  TAU_HIGH = 1.20_EB
               END IF
            CASE ('Female','female','FEMALE')
               IF (VELOCITY_DIST < 0) THEN
                  VELOCITY_DIST = 1
                  VEL_MEAN = 1.15_EB
                  VEL_PARA = 0.20_EB
                  VEL_LOW  = 0.95_EB
                  VEL_HIGH = 1.35_EB
               END IF
               IF (DIAMETER_DIST < 0) THEN
                  DIAMETER_DIST = 1
                  DIA_MEAN = 0.48_EB
                  DIA_PARA = 0.04_EB
                  DIA_LOW  = 0.44_EB
                  DIA_HIGH = 0.52_EB
                  D_TORSO_MEAN    = 0.28_EB
                  D_SHOULDER_MEAN = 0.18_EB
               END IF
               IF (TAU_EVAC_DIST < 0) THEN
                  TAU_EVAC_DIST = 1
                  TAU_MEAN = 1.00_EB
                  TAU_PARA = 0.10_EB
                  TAU_LOW  = 0.80_EB
                  TAU_HIGH = 1.20_EB
               END IF
            CASE ('Child','child','CHILD')
               IF (VELOCITY_DIST < 0) THEN
                  VELOCITY_DIST = 1
                  VEL_MEAN = 0.90_EB
                  VEL_PARA = 0.30_EB
                  VEL_LOW  = 0.60_EB
                  VEL_HIGH = 1.20_EB
               END IF
               IF (DIAMETER_DIST < 0) THEN
                  DIAMETER_DIST = 1
                  DIA_MEAN = 0.42_EB
                  DIA_PARA = 0.03_EB
                  DIA_LOW  = 0.39_EB
                  DIA_HIGH = 0.45_EB
                  D_TORSO_MEAN    = 0.24_EB
                  D_SHOULDER_MEAN = 0.14_EB
               END IF
               IF (TAU_EVAC_DIST < 0) THEN
                  TAU_EVAC_DIST = 1
                  TAU_MEAN = 1.00_EB
                  TAU_PARA = 0.10_EB
                  TAU_LOW  = 0.80_EB
                  TAU_HIGH = 1.20_EB
               END IF
            CASE ('Elderly','elderly','ELDERLY')
               IF (VELOCITY_DIST < 0) THEN
                  VELOCITY_DIST = 1
                  VEL_MEAN = 0.80_EB
                  VEL_PARA = 0.30_EB
                  VEL_LOW  = 0.50_EB
                  VEL_HIGH = 1.10_EB
               END IF
               IF (DIAMETER_DIST < 0) THEN
                  DIAMETER_DIST = 1
                  DIA_MEAN = 0.50_EB
                  DIA_PARA = 0.04_EB
                  DIA_LOW  = 0.46_EB
                  DIA_HIGH = 0.54_EB
                  D_TORSO_MEAN    = 0.30_EB
                  D_SHOULDER_MEAN = 0.18_EB
               END IF
               IF (TAU_EVAC_DIST < 0) THEN
                  TAU_EVAC_DIST = 1
                  TAU_MEAN = 1.00_EB
                  TAU_PARA = 0.10_EB
                  TAU_LOW  = 0.80_EB
                  TAU_HIGH = 1.20_EB
               END IF
            CASE ('null')
               ! Do nothing, use the defaults
               IF (MYID == MAX(0,EVAC_PROCESS)) THEN
                  WRITE (LU_ERR,'(A,A)') ' WARNING: PERS ',TRIM(ID),' no DEFAULT_PROPERTIES given'
               END IF
            CASE Default
               IF (MYID == MAX(0,EVAC_PROCESS)) THEN
                  WRITE(MESSAGE,'(A,A,A)') 'ERROR: PERS ',TRIM(ID),' problem with DEFAULT_PROPERTIES'
                  CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
               END IF
            END SELECT

         END IF

         IF (PRE_MEAN < 0._EB .OR. PRE_LOW < 0._EB) THEN
            IF (MYID == MAX(0,EVAC_PROCESS)) THEN
               WRITE(MESSAGE,'(A,A,A)') 'ERROR: PERS ',TRIM(ID), ' PRE-evacuation time should positive.'
               CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
            END IF
         END IF

         IF (FED_ACTIVITY < 1 .OR. FED_ACTIVITY > 3) THEN
            IF (MYID == MAX(0,EVAC_PROCESS)) THEN
               WRITE(MESSAGE,'(A,A,A,I3)') 'ERROR: PERS ',TRIM(ID), ' not a valid FED activity switch: ',FED_ACTIVITY
               CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
            END IF
         END IF

         DIAMETER_DIST = MAX(0,DIAMETER_DIST)
         VELOCITY_DIST = MAX(0,VELOCITY_DIST)
         TAU_EVAC_DIST = MAX(0,TAU_EVAC_DIST)
         IF (DIA_MEAN < 0.0_EB) THEN
            SELECT CASE (DIAMETER_DIST)
            CASE (1)  ! Uniform
               DIA_MEAN = 0.5_EB*(DIA_HIGH+DIA_LOW)
            CASE (3)  ! Gamma: mean = alpha*beta
               DIA_MEAN = DIA_PARA*DIA_PARA2
            CASE (6)  ! Beta: mean alpha/(alpha+beta)
               DIA_MEAN = DIA_PARA/(DIA_PARA*DIA_PARA2)
            CASE (8)  ! Weibull (Exp, alpha=1):
               ! mean = (1/lambda)*GammaFunc(1 + (1/alpha) )
               ! median = (1/lambda)* (ln(2))^(1/alpha)
               ! exp: mean = (1/lambda), median = ln(2)/lambda
               DIA_MEAN = (1.0_EB/DIA_PARA2)*DCDFLIB_Gamma(1.0_EB+(1.0_EB/DIA_PARA))
            CASE (9)  ! Gumbel
               ! mean = gamma / alpha, gamma = 0.5772156649015328606_EB
               ! median = -(1/alpha)*ln(ln(2))
               DIA_MEAN = 0.5772156649015328606_EB / DIA_PARA
            CASE Default
               CONTINUE
            END SELECT
         END IF

         !
         ! Avatar colors, integer RGB(3), e.g., (23,255,0)
         IF (DEAD_COLOR /= 'null') THEN
            CALL COLOR2RGB(DEAD_RGB,DEAD_COLOR)
            EVAC_AVATAR_RGB(1:3,EVAC_AVATAR_NCOLOR) = DEAD_RGB
         END IF
         IF (ANY(AVATAR_RGB < 0) .AND. AVATAR_COLOR=='null') AVATAR_COLOR = 'ROYAL BLUE 4'
         IF (AVATAR_COLOR /= 'null') CALL COLOR2RGB(AVATAR_RGB,AVATAR_COLOR)
         IF (COLOR_METHOD_TMP == 3) THEN
            i_avatar_color = i_avatar_color + 1
            EVAC_AVATAR_RGB(1:3,i_avatar_color) = AVATAR_RGB
         END IF
         !
         IF (MYID /= MAX(0,EVAC_PROCESS)) CYCLE READ_PERS_LOOP
         !
         PCP=>EVAC_PERSON_CLASSES(N)
         !
         IF (TRIM(ID) /= 'null') THEN
            DO I = 1, N-1
               IF (TRIM(ID) == TRIM(EVAC_PERSON_CLASSES(I)%ID)) THEN
                  WRITE(MESSAGE,'(A,I4,A,I4,A,A)') 'ERROR: PERS lines',I,' and',N,', ID strings are not unique: ',TRIM(ID)
                  CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN

               END IF
            END DO
         END IF
         ! Colors, integer RGB(3), e.g., (23,255,0)
         IF (ANY(RGB < 0) .AND. COLOR=='null') COLOR = 'ROYAL BLUE 4'
         IF (COLOR /= 'null') CALL COLOR2RGB(RGB,COLOR)
         PCP%RGB = RGB
         PCP%AVATAR_RGB = AVATAR_RGB
         IF (COLOR_METHOD_TMP == 3) PCP%Avatar_Color_Index = i_avatar_color
         PCP%ID = ID
         !
         PCP%D_mean = DIA_MEAN
         PCP%I_DIA_DIST  = DIAMETER_DIST
         PCP%D_low  = DIA_LOW
         PCP%D_high = DIA_HIGH
         PCP%D_para = DIA_PARA
         PCP%D_para2 = DIA_PARA2
         !
         PCP%V_mean = VEL_MEAN
         PCP%I_VEL_DIST  = VELOCITY_DIST
         PCP%V_low  = VEL_LOW
         PCP%V_high = VEL_HIGH
         PCP%V_para = VEL_PARA
         PCP%V_para2 = VEL_PARA2
         !
         PCP%Tau_mean = TAU_MEAN
         PCP%I_TAU_DIST  = TAU_EVAC_DIST
         PCP%Tau_low  = TAU_LOW
         PCP%Tau_high = TAU_HIGH
         PCP%Tau_para = TAU_PARA
         PCP%Tau_para2 = TAU_PARA2
         !
         PCP%Tpre_mean  = PRE_MEAN
         PCP%I_PRE_DIST = PRE_EVAC_DIST
         PCP%Tpre_low   = PRE_LOW
         PCP%Tpre_high  = PRE_HIGH
         PCP%Tpre_para  = PRE_PARA
         PCP%Tpre_para2  = PRE_PARA2
         !
         PCP%Tdet_mean  = DET_MEAN
         PCP%I_DET_DIST = DET_EVAC_DIST
         PCP%Tdet_low   = DET_LOW
         PCP%Tdet_high  = DET_HIGH
         PCP%Tdet_para  = DET_PARA
         PCP%Tdet_para2  = DET_PARA2
         !
         PCP%A       = FCONST_A
         PCP%B       = FCONST_B
         PCP%Lambda  = L_NON_SP
         PCP%C_Young = C_YOUNG
         PCP%Gamma   = GAMMA
         PCP%Kappa   = KAPPA
         !
         PCP%r_torso    = 0.5_EB*D_TORSO_MEAN
         PCP%r_shoulder = 0.5_EB*D_SHOULDER_MEAN
         PCP%Tau_iner   = TAU_ROT  ! s
         IF (M_INERTIA < 0.0_EB ) THEN
            PCP%m_iner = -M_INERTIA*(0.25_EB*PCP%D_mean**2+PCP%r_torso**2)**2 / (0.27_EB**2+0.16_EB**2)**2
         ELSE
            PCP%m_iner = M_INERTIA  ! kg m2
         END IF
         PCP%m_agent = ABS(MASS_OF_AGENT) ! kg, default male body size mass

         PCP%FAC_V0_UP = FAC_V0_UP
         PCP%FAC_V0_DOWN = FAC_V0_DOWN
         PCP%FAC_V0_HORI = FAC_V0_HORI

         !Issue1547: Save the read in value in the structured type of the person types.
         PCP%MAXIMUM_V0_FACTOR = MAXIMUM_V0_FACTOR
         !
      END DO READ_PERS_LOOP
24    REWIND(LU_INPUT)
      COLOR_METHOD = COLOR_METHOD_TMP
      ! Check if falling down model is used or not
      IF (TAU_FALL_DOWN >= 0.0_EB) THEN
         L_FALLING_MODEL=.TRUE.
      ELSE
         L_FALLING_MODEL=.FALSE.
      END IF
      IF (T_ASET_HAWK > 0.0_EB) THEN
         C_HAWK = 1000000000000000.0_EB ! Default, all doves if t_i < t_aset-t_0
      END IF

      IF (MAX_INITIAL_OVERLAP <= -99998.0_EB .AND. DENS_INIT <= -99998.0_EB) THEN
         ! No inputs given, use the defaults
         DENS_INIT = 1.0_EB
         MAX_INITIAL_OVERLAP = 0.08_EB
      END IF
      IF (MAX_INITIAL_OVERLAP <= -99998.0_EB .AND. DENS_INIT > -99998.0_EB) THEN
         ! Just the density given => old input format behavour
         IF (DENS_INIT > 2.0_EB) THEN
            MAX_INITIAL_OVERLAP = 0.0_EB
         ELSE
            MAX_INITIAL_OVERLAP = 0.08_EB
         END IF
         IF (DENS_INIT > 10.0_EB)  MAX_INITIAL_OVERLAP = -0.1_EB
      END IF
      ! The default is at least 0.08 m empty space between the agents
      IF (MAX_INITIAL_OVERLAP < -99998.0_EB) MAX_INITIAL_OVERLAP = 0.08_EB
      DENS_INIT = MAX(DENS_INIT,1.0_EB)

      R_HERDING = MAX(0.1_EB,R_HERDING)  ! Avoid divisions by zero
      IF (GROUP_DENS <= 0.01_EB) GROUP_DENS = 0.25_EB
      IF (GROUP_DENS > 3.50_EB) GROUP_DENS = 3.50_EB
      DENS_INIT = MAX(GROUP_DENS,DENS_INIT)
      IF (TDET_SMOKE_DENS < 0.0_EB) TDET_SMOKE_DENS = HUGE(TDET_SMOKE_DENS)

      IF(SMOKE_MIN_SPEED < 0.0_EB .AND. SMOKE_MIN_SPEED_FACTOR < 0.0_EB) THEN
         SMOKE_MIN_SPEED_FACTOR = 0.1_EB ! Default min speed = v0_i*v_s_min_fac
      END IF
      IF(SMOKE_MIN_SPEED > -99998.0_EB .AND. SMOKE_MIN_SPEED_FACTOR > -99998.0_EB) THEN
         WRITE(MESSAGE,'(A,I3,A)') 'ERROR: READ_EVAC both SMOKE_MIN_SPEED and SMOKE_MIN_SPEED_FACTOR are given'
         CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
      END IF
      IF(SMOKE_MIN_SPEED_VISIBILITY < 0.01_EB) THEN
         SMOKE_MIN_SPEED_VISIBILITY = 0.01_EB  ! No divisions by zero
      END IF

      SELECT CASE (COLOR_METHOD)
      CASE (-1:5)
         CONTINUE
      CASE (7)
         COLOR_METHOD = -1
         IF (MYID==MAX(0,EVAC_PROCESS)) WRITE (LU_ERR,'(A)') &
              ' WARNING: COLOR_METHOD=7 is not defined anymore, the default (-1) is used.'
      CASE Default
         WRITE(MESSAGE,'(A,I3,A)') 'ERROR: READ_EVAC COLOR METHOD',COLOR_METHOD, ' is not defined'
         CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
      END SELECT

      IF (COLOR_METHOD >= 0) EVAC_N_QUANTITIES = EVAC_N_QUANTITIES + 1
      IF (OUTPUT_MOTIVE_FORCE) EVAC_N_QUANTITIES = EVAC_N_QUANTITIES + 1
      IF (OUTPUT_MOTIVE_ANGLE) EVAC_N_QUANTITIES = EVAC_N_QUANTITIES + 1
      IF (OUTPUT_FED) EVAC_N_QUANTITIES = EVAC_N_QUANTITIES + 1
      IF (OUTPUT_SPEED) EVAC_N_QUANTITIES = EVAC_N_QUANTITIES + 1
      IF (OUTPUT_DENSITY) EVAC_N_QUANTITIES = EVAC_N_QUANTITIES + 1
      IF (OUTPUT_OMEGA) EVAC_N_QUANTITIES = EVAC_N_QUANTITIES + 1
      IF (OUTPUT_ACCELERATION) EVAC_N_QUANTITIES = EVAC_N_QUANTITIES + 1
      IF (OUTPUT_CONTACT_FORCE) EVAC_N_QUANTITIES = EVAC_N_QUANTITIES + 1
      IF (OUTPUT_TOTAL_FORCE) EVAC_N_QUANTITIES = EVAC_N_QUANTITIES + 1
      !Issue1547: Add the counter of the special quantities that are stored in the .prt5 files.
      IF (OUTPUT_NERVOUSNESS) EVAC_N_QUANTITIES = EVAC_N_QUANTITIES + 1

      IF (EVAC_N_QUANTITIES > 0) THEN
         ALLOCATE(EVAC_QUANTITIES_INDEX(EVAC_N_QUANTITIES),STAT=IZERO)
         CALL ChkMemErr('READ','EVAC_QUANTITIES_INDEX',IZERO)

         n = 1
         IF (COLOR_METHOD >= 0) THEN
            EVAC_QUANTITIES_INDEX(n)=247  ! HUMAN_COLOR
            n = n + 1
         END IF
         IF (OUTPUT_MOTIVE_FORCE) THEN
            EVAC_QUANTITIES_INDEX(n)=240
            n = n + 1
         END IF
         IF (OUTPUT_FED) THEN
            EVAC_QUANTITIES_INDEX(n)=241
            n = n + 1
         END IF
         IF (OUTPUT_SPEED) THEN
            EVAC_QUANTITIES_INDEX(n)=242
            n = n + 1
         END IF
         IF (OUTPUT_OMEGA) THEN
            EVAC_QUANTITIES_INDEX(n)=243
            n = n + 1
         END IF
         IF (OUTPUT_ACCELERATION) THEN
            EVAC_QUANTITIES_INDEX(n)=244
            n = n + 1
         END IF
         IF (OUTPUT_CONTACT_FORCE) THEN
            EVAC_QUANTITIES_INDEX(n)=245
            n = n + 1
         END IF
         IF (OUTPUT_TOTAL_FORCE) THEN
            EVAC_QUANTITIES_INDEX(n)=246
            n = n + 1
         END IF
         IF (OUTPUT_MOTIVE_ANGLE) THEN
            EVAC_QUANTITIES_INDEX(n)=248
            n = n + 1
         END IF
         IF (OUTPUT_DENSITY) THEN
            EVAC_QUANTITIES_INDEX(n)=249
            n = n + 1
         END IF
         !Issue1547: The output quantities are numbered in "data.f90", see it also (search for NERVOUSNESS)
         !Issue1547: Here "the value" of OUTPUT_NERVOUSNES is saved for later use, the variable
         !Issue1547: OUTPUT_NERVOUSNES is defined only for this subroutine.
         IF (OUTPUT_NERVOUSNESS) THEN
            EVAC_QUANTITIES_INDEX(n)=250
            n = n + 1
         END IF

         IF ( n-1 /= EVAC_N_QUANTITIES ) THEN
            WRITE(MESSAGE,'(A,2I4,A)') 'ERROR: Evac output quantities ',EVAC_N_QUANTITIES,n-1, ' Some bug in the program.'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END IF
      END IF

    END SUBROUTINE READ_PERS

    SUBROUTINE READ_EXIT
      IMPLICIT NONE
      !
      ! Read the EXIT lines
      !
      ! Local variables
      INTEGER nm, i1, i2, j1, j2, i, ii, iii, N_TMP
      LOGICAL L_TMP
      TYPE (EVAC_EXIT_TYPE), POINTER :: PEX=>NULL()
      TYPE (EVAC_STRS_TYPE),  POINTER :: STRP=>NULL()
      TYPE (MESH_TYPE), POINTER :: M=>NULL()
      !
      N_TMP = 0
      READ_EXIT_LOOP: DO N = 1, N_EXITS
         !
         ID            = 'null'
         RGB           = -1
         COLOR         = 'null'
         XB            = 0.0_EB
         IOR           = 0
         FLOW_FIELD_ID = 'null'
         VENT_FFIELD   = 'null'
         EVAC_ID       = 'null'
         PERS_ID       = 'null'
         MESH_ID       = 'null'
         EVAC_MESH     = 'null'
         CHECK_FLOW    = .FALSE.
         COUNT_ONLY    = .FALSE.
         COUNT_DENSITY = .FALSE.
         DEFINE_MESH   = .TRUE.
         SHOW          = .TRUE.
         MAX_FLOW      = 0.0_EB
         WIDTH         = 0.0_EB
         HEIGHT        = 2.0_EB
         TIME_OPEN     = -HUGE(TIME_OPEN)
         TIME_CLOSE    = HUGE(TIME_CLOSE)
         XYZ(:)        = HUGE(XYZ)
         XYZ_SMOKE(:)  = HUGE(XYZ_SMOKE)
         COLOR_INDEX   = -1
         LOCKED_WHEN_CLOSED = .FALSE.
         TARGET_WHEN_CLOSED = .FALSE.
         WAIT_AT_XYZ        = .FALSE.
         KNOWN_DOOR         = .FALSE.
         !
         CALL CHECKREAD('EXIT',LU_INPUT,IOS)  ; IF (STOP_STATUS==SETUP_STOP) RETURN
         IF (IOS == 1) THEN
            EXIT READ_EXIT_LOOP
         END IF
         READ(LU_INPUT,EXIT,END=26,IOSTAT=IOS)
         !
         IF (VENT_FFIELD /= 'null') DEFINE_MESH = .FALSE. ! Old input, user gives the door flow meshes
         IF (COUNT_ONLY) DEFINE_MESH = .FALSE.
         ! Old input used COLOR_INDEX, next lines are needed for that
         IF (MYID==MAX(0,EVAC_PROCESS) .AND. COLOR_INDEX/=-1) WRITE (LU_ERR,'(A,A)') &
              ' WARNING: keyword COLOR_INDEX is replaced by COLOR at EXIT line ',TRIM(ID)
         IF (COLOR_INDEX == 1) COLOR = 'BLACK'
         IF (COLOR_INDEX == 2) COLOR = 'YELLOW'
         IF (COLOR_INDEX == 3) COLOR = 'BLUE'
         IF (COLOR_INDEX == 4) COLOR = 'RED'
         IF (COLOR_INDEX == 5) COLOR = 'GREEN'
         IF (COLOR_INDEX == 6) COLOR = 'MAGENTA'
         IF (COLOR_INDEX == 7) COLOR = 'CYAN'

         IF (COUNT_DENSITY) THEN
            COUNT_ONLY = .TRUE.
            IOR = +3
         END IF
         ! Colors, integer RGB(3), e.g., (23,255,0)
         IF (ANY(RGB < 0) .AND. COLOR=='null') COLOR = 'FOREST GREEN'
         IF (COLOR /= 'null') CALL COLOR2RGB(RGB,COLOR)
         IF (COLOR_METHOD == 4 .AND. .NOT.COUNT_ONLY) THEN
            i_avatar_color = i_avatar_color + 1
            EVAC_AVATAR_RGB(1:3,i_avatar_color) = RGB
         END IF

         IF (MYID /= MAX(0,EVAC_PROCESS)) CYCLE READ_EXIT_LOOP

         PEX=>EVAC_EXITS(N)

         IF (.NOT. COUNT_ONLY) N_TMP = N_TMP + 1 ! real exits + door index (not count_only ones)

         EVAC_EXITS(N)%I_EMESH_EXITS = 0  ! Default, for "no flow field" exits
         PEX%RGB = RGB
         IF (COLOR_METHOD == 4 .AND. .NOT.COUNT_ONLY) PEX%Avatar_Color_Index = i_avatar_color

         IF (TRIM(ID) /= 'null') THEN
            DO I = 1, N-1
               IF (TRIM(ID) == TRIM(EVAC_EXITS(I)%ID)) THEN
                  WRITE(MESSAGE,'(A,I4,A,I4,A,A)') 'ERROR: EXIT lines',I,' and',N,', ID strings are not unique: ',TRIM(ID)
                  CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
               END IF
            END DO
         END IF
         IF (EVAC_MESH /= 'null') THEN
            MESH_ID = TRIM(EVAC_MESH)
            IF (MYID==MAX(0,EVAC_PROCESS)) WRITE (LU_ERR,'(A,A)') &
                 ' WARNING: keyword EVAC_MESH is replaced by MESH_ID at EXIT line ', TRIM(ID)
         END IF

         ! Check that the exit is properly specified

         ! Reorder XB coordinates if necessary
         CALL CHECK_XB(XB)
         !
         ! Check which evacuation mesh
         II = 0
         PEX_MeshLoop: DO I = 1, NMESHES
            IF (EVACUATION_ONLY(I) .AND. EVACUATION_SKIP(I)) THEN
               IF (Is_Within_Bounds(XB(1),XB(2),XB(3),XB(4),XB(5),XB(6),&
                    MESHES(i)%XS,MESHES(i)%XF,MESHES(i)%YS,MESHES(i)%YF,MESHES(i)%ZS,MESHES(i)%ZF, 0._EB, 0._EB, 0._EB)) THEN
                  IF (TRIM(MESH_ID) == 'null' .OR. TRIM(MESH_ID) == TRIM(MESH_NAME(i))) THEN
                     II = II + 1
                     PEX%IMESH = i
                  END IF
               END IF
            END IF
         END DO PEX_MeshLoop
         PEX_StrsLoop: DO I = 1, N_STRS
            IF (Is_Within_Bounds(XB(1), XB(2), XB(3), XB(4), XB(5), XB(6), &
                 EMESH_STAIRS(I)%XB(1), EMESH_STAIRS(I)%XB(2), EMESH_STAIRS(I)%XB(3), &
                 EMESH_STAIRS(I)%XB(4), EMESH_STAIRS(I)%XB(5), EMESH_STAIRS(I)%XB(6), &
                 0._EB, 0._EB, 0._EB)) THEN
               IF (TRIM(MESH_ID) == 'null' .OR. TRIM(MESH_ID) == ('Emesh_' // TRIM(EMESH_STAIRS(I)%ID))) THEN
                  DEFINE_MESH = .FALSE.
               END IF
            END IF
         END DO PEX_StrsLoop

         IF (PEX%IMESH == 0) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: EXIT ',TRIM(ID), ' problem with IMESH, no mesh found'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END IF
         IF (ii > 1) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: EXIT ',TRIM(ID), ' not an unique mesh found '
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END IF

         IF (DEFINE_MESH) EVAC_EXITS(N)%I_EMESH_EXITS = N_TMP
         NM = PEX%IMESH
         M => MESHES(PEX%IMESH)

         IF (XB(1)/=XB(2) .AND. XB(3)/=XB(4) .AND. .NOT. ABS(IOR)==3) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: EXIT ',TRIM(ID),' must be a vertical plane'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         ENDIF
         IF (XB(1)==XB(2) .AND. XB(3)==XB(4) .AND. ABS(IOR)==3) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: Plane counter EXIT ',TRIM(ID),' must be a horizontal plane'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         ENDIF

         ! User input
         PEX%X1 = XB(1)
         PEX%X2 = XB(2)
         PEX%Y1 = XB(3)
         PEX%Y2 = XB(4)
         PEX%Z1 = XB(5)
         PEX%Z2 = XB(6)

         ! Move user input to mesh cell boundaries
         XB(1) = MAX(XB(1),M%XS)
         XB(2) = MIN(XB(2),M%XF)
         XB(3) = MAX(XB(3),M%YS)
         XB(4) = MIN(XB(4),M%YF)
         XB(5) = MAX(XB(5),M%ZS)
         XB(6) = MIN(XB(6),M%ZF)

         I1 = NINT( GINV(XB(1)-M%XS,1,nm)*M%RDXI )
         I2 = NINT( GINV(XB(2)-M%XS,1,nm)*M%RDXI )
         J1 = NINT( GINV(XB(3)-M%YS,2,nm)*M%RDETA)
         J2 = NINT( GINV(XB(4)-M%YS,2,nm)*M%RDETA)

         XB(1) = M%X(I1)
         XB(2) = M%X(I2)
         XB(3) = M%Y(J1)
         XB(4) = M%Y(J2)

         IF ( ABS(XB(1)-PEX%X1)>1.E-4_EB .OR. ABS(XB(2)-PEX%X2)>1.E-4_EB .OR. &
              ABS(XB(3)-PEX%Y1)>1.E-4_EB .OR. ABS(XB(4)-PEX%Y2)>1.E-4_EB ) THEN
            WRITE(LU_ERR,fmt='(a,a,a,a)') ' WARNING: Exit ',TRIM(ID),' XB adjusted to mesh ',TRIM(MESH_NAME(nm))
            WRITE(LU_ERR,fmt='(a,6f12.4)') ' Old XB:', PEX%X1,PEX%X2,PEX%Y1,PEX%Y2,PEX%Z1,PEX%Z2
            WRITE(LU_ERR,fmt='(a,6f12.4)') ' New XB:', XB(1:6)
         END IF

         II = (I1+I2)/2
         JJ = (J1+J2)/2
         SELECT CASE (IOR)
         CASE (-1)
            II = II + 1
         CASE (-2)
            JJ = JJ + 1
         END SELECT
         IF (M%SOLID(M%CELL_INDEX(II,JJ,1)) .AND. .NOT. COUNT_ONLY) THEN
            WRITE(LU_ERR,fmt='(a,a,a)') ' WARNING: Exit ',TRIM(ID),' problem with XB, mid point facing solid'
         END IF

         ! Coordinates are lined up with the mesh.
         PEX%X1 = XB(1)
         PEX%X2 = XB(2)
         PEX%Y1 = XB(3)
         PEX%Y2 = XB(4)
         PEX%Z1 = XB(5)
         PEX%Z2 = XB(6)
         !
         ! Find the implicitely generated mesh behind this exit/door
         I = PEX%I_EMESH_EXITS
         IF (VENT_FFIELD == 'null' .AND. I > 0) THEN
            IF (EMESH_EXITS(I)%NEXIT == N .AND. EMESH_EXITS(I)%IS_EXIT) THEN
               VENT_FFIELD = TRIM(MESH_NAME(EMESH_EXITS(I)%IMESH))
            END IF
         END IF

         PEX%IOR = IOR
         PEX%ID  = TRIM(ID)
         PEX%GRID_NAME  = TRIM(FLOW_FIELD_ID)
         PEX%CHECK_FLOW = CHECK_FLOW
         PEX%VENT_FFIELD= VENT_FFIELD
         PEX%INODE      = 0
         PEX%T_first    = T_BEGIN
         PEX%T_last     = T_BEGIN
         PEX%ICOUNT     = 0
         PEX%Flow_max   = 0.0_EB
         PEX%TIME_OPEN  = TIME_OPEN
         PEX%TIME_CLOSE = TIME_CLOSE
         PEX%LOCKED_WHEN_CLOSED = LOCKED_WHEN_CLOSED
         PEX%TARGET_WHEN_CLOSED = TARGET_WHEN_CLOSED
         PEX%WAIT_AT_XYZ        = WAIT_AT_XYZ
         PEX%IMODE      = -1 ! Exit is open by default
         IF (TIME_OPEN > TIME_CLOSE) THEN
            PEX%IMODE = -1
         ELSE
            IF (TIME_OPEN > T_BEGIN) PEX%IMODE = +2
         END IF
         IF (CHECK_FLOW) PEX%Flow_max   = MAX_FLOW
         PEX%COUNT_ONLY = .FALSE.
         IF (COUNT_ONLY) PEX%COUNT_ONLY = .TRUE.
         PEX%COUNT_DENSITY = COUNT_DENSITY
         PEX%EVAC_ID = EVAC_ID
         PEX%PERS_ID = PERS_ID
         IF (COUNT_ONLY  ) SHOW = .FALSE.
         PEX%SHOW    = SHOW
         PEX%HEIGHT = HEIGHT
         PEX%KNOWN_DOOR = KNOWN_DOOR


         SELECT CASE (IOR)
         CASE (-1,+1)
            IF (WIDTH <= 0.0_EB) THEN
               PEX%Width = XB(4) - XB(3)
            ELSE
               PEX%Width = WIDTH
            END IF
            PEX%ORIENTATION(1)=REAL(SIGN(1,IOR),EB)
         CASE (-2,+2)
            IF (WIDTH <= 0.0_EB) THEN
               PEX%Width = XB(2) - XB(1)
            ELSE
               PEX%Width = WIDTH
            END IF
            PEX%ORIENTATION(2)=REAL(SIGN(1,IOR),EB)
         CASE (-3,+3)
            IF ( (XB(4)-XB(3)) <= 0.0_EB .OR. (XB(2)-XB(1)) <= 0.0_EB) THEN
               WRITE(MESSAGE,'(A,A,A)') 'ERROR: EXIT ',TRIM(ID),' IOR=+-3 but not a horizontal plane'
               CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
            END IF
            PEX%ORIENTATION(3)=REAL(SIGN(1,IOR),EB)
            PEX%Width = ABS(XB(2)-XB(1)) * ABS(XB(4)-XB(3))  ! Area of the plane counter
         CASE (0)
            IF ( (XB(4)-XB(3)) <= 0.0_EB .OR. (XB(2)-XB(1)) <= 0.0_EB) THEN
               WRITE(MESSAGE,'(A,A,A)') 'ERROR: EXIT ',TRIM(ID),' no IOR but not a horizontal plane'
               CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
            END IF
         CASE Default
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: EXIT ',TRIM(ID),' problem with IOR'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END SELECT

         L_TMP=.FALSE.
         DO i = 1, NMESHES
            IF (.NOT. EVACUATION_ONLY(i)) CYCLE
            IF (TRIM(FLOW_FIELD_ID)==TRIM(MESH_NAME(i))) THEN
               L_TMP=.TRUE.
               EXIT
            END IF
         END DO
         IF (.NOT.(TRIM(FLOW_FIELD_ID)=='null' .OR. L_TMP)) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: EXIT ',TRIM(ID),' problem with FLOW_FIELD_ID'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END IF
         L_TMP=.FALSE.
         DO i = 1, NMESHES
            IF (.NOT. EVACUATION_ONLY(i)) CYCLE
            IF (TRIM(VENT_FFIELD)==TRIM(MESH_NAME(i))) THEN
               L_TMP=.TRUE.
               EXIT
            END IF
         END DO
         IF (.NOT.(TRIM(VENT_FFIELD)=='null' .OR. L_TMP)) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: EXIT ',TRIM(ID),' problem with VENT_FFIELD'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END IF
         !
         ! Check which vent field. If VENT_FFIELD is not found, use the main evac grid.
         PEX%I_VENT_FFIELD = 0
         PEX_Mesh2Loop: DO i = 1, NMESHES
            IF ( EVACUATION_ONLY(I) .AND. (TRIM(MESH_NAME(i)) == TRIM(PEX%VENT_FFIELD)) ) THEN
               IF ( (PEX%Z1 >= MESHES(i)%ZS .AND. PEX%Z2 <= MESHES(i)%ZF).AND. &
                    (PEX%Y1 >= MESHES(i)%YS .AND. PEX%Y2 <= MESHES(i)%YF).AND. &
                    (PEX%X1 >= MESHES(i)%XS .AND. PEX%X2 <= MESHES(i)%XF)) THEN
                  PEX%I_VENT_FFIELD = i
                  EXIT PEX_Mesh2Loop
               END IF
            END IF
         END DO PEX_Mesh2Loop
         ! If no vent field is given, then use the main evac grid.
         IF (PEX%I_VENT_FFIELD == 0) THEN
            PEX%I_VENT_FFIELD = PEX%IMESH
            PEX%VENT_FFIELD = TRIM(MESH_NAME(PEX%IMESH))
         END IF

         PEX%FED_MESH = 0
         IF (XYZ(1) < HUGE(XYZ)) THEN
            PEX%X = XYZ(1)
            PEX%Y = XYZ(2)
            PEX%Z = 0.5_EB*(XB(5)+XB(6))
         ELSE
            PEX%X = 0.5_EB*(XB(1)+XB(2))
            PEX%Y = 0.5_EB*(XB(3)+XB(4))
            PEX%Z = 0.5_EB*(XB(5)+XB(6))
            SELECT CASE (IOR)
            CASE (-1,+1)
               PEX%X = PEX%X - IOR*0.20_EB
            CASE (-2,+2)
               PEX%Y = PEX%Y - IOR*0.10_EB
            END SELECT
         END IF
         !
         ! Check which evacuation floor.  Now there may be overlapping meshes.
         ii  = 0
         iii = 0
         PEX_Mesh3Loop: DO i = 1, NMESHES
            IF (EVACUATION_ONLY(I) .AND. EVACUATION_SKIP(I)) THEN
               IF ( (PEX%Z >= MESHES(i)%ZS .AND. PEX%Z <= MESHES(i)%ZF).AND. &
                    (PEX%Y >= MESHES(i)%YS .AND. PEX%Y <= MESHES(i)%YF).AND. &
                    (PEX%X >= MESHES(i)%XS .AND. PEX%X <= MESHES(i)%XF)) THEN
                  IF (TRIM(MESH_ID) == 'null' .OR. TRIM(MESH_ID) == TRIM(MESH_NAME(i))) THEN
                     ii  = ii + 1
                     iii = i
                  END IF
               END IF
            END IF
         END DO PEX_Mesh3Loop
         IF (ii == 0 .AND. .NOT.COUNT_ONLY .AND. SHOW) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: EXIT line ',TRIM(PEX%ID), ' problem with XYZ, no mesh found'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END IF
         IF (ii > 1) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: EXIT line ',TRIM(PEX%ID), ' problem with XYZ, not an unique mesh found'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END IF
         IF (iii /= PEX%IMESH .AND. .NOT.COUNT_ONLY .AND. SHOW) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: EXIT line ',TRIM(PEX%ID), ' problem with XYZ, wrong mesh found'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END IF
         M => MESHES(PEX%IMESH)
         II = FLOOR(M%CELLSI(FLOOR((PEX%X-M%XS)*M%RDXINT))+ 1.0_EB)
         JJ = FLOOR(M%CELLSJ(FLOOR((PEX%Y-M%YS)*M%RDYINT))+ 1.0_EB)
         KK = FLOOR(M%CELLSK(FLOOR((PEX%Z-M%ZS)*M%RDZINT))+ 1.0_EB)
         IF (M%SOLID(M%CELL_INDEX(II,JJ,KK)) .AND. .NOT.COUNT_ONLY .AND. SHOW) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: EXIT line ',TRIM(PEX%ID), ' problem with XYZ, inside solid'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END IF

         ! PEX%Z is used to plot the door on the correct height in Smokeview.
         PEX%Z = PEX%Z + 0.5_EB*PEX%Height - EVACUATION_Z_OFFSET(PEX%IMESH)
         PEX%XYZ_Z = PEX%Z

         IF (XYZ_SMOKE(1) < HUGE(XYZ_SMOKE)) THEN
            PEX%Xsmoke = XYZ_SMOKE(1)
            PEX%Ysmoke = XYZ_SMOKE(2)
            PEX%Zsmoke = XYZ_SMOKE(3)
         ELSE
            PEX%Xsmoke = PEX%X
            PEX%Ysmoke = PEX%Y
            PEX%Zsmoke = 0.5_EB*(XB(5)+XB(6)) - EVACUATION_Z_OFFSET(PEX%IMESH) + HUMAN_SMOKE_HEIGHT
         END IF

         ! Check if exit is in Stairs
         PEX%STR_INDX = 0
         PEX%STR_SUB_INDX = 0
         CheckExitStrLoop: DO i = 1, N_STRS
            IF (EVAC_STRS(i)%IMESH==PEX%IMESH ) THEN
               STRP=>EVAC_STRS(i)
               PEX%STR_INDX = i
               DO j = 1,STRP%N_NODES
                  IF ( Is_Within_Bounds(PEX%X1,PEX%X2,PEX%Y1,PEX%Y2,PEX%Z1,PEX%Z2, &
                       STRP%XB_NODE(j,1), STRP%XB_NODE(j,2), STRP%XB_NODE(j,3),STRP%XB_NODE(j,4), &
                       STRP%XB_NODE(j,5), STRP%XB_NODE(j,6), 0._EB, 0._EB, 0._EB)) THEN
                     PEX%STR_SUB_INDX = j
                     EXIT CheckExitStrLoop
                  END IF
               END DO
            END IF
         END DO CheckExitStrLoop
         !
         ! Check, which fire grid and i,j,k (xyz)
         PEX_SmokeLoop: DO i = 1, NMESHES
            IF (.NOT. EVACUATION_ONLY(I)) THEN
               IF ( Is_Within_Bounds(PEX%Xsmoke,PEX%Xsmoke,PEX%Ysmoke,PEX%Ysmoke,PEX%Zsmoke,PEX%Zsmoke,&
                    MESHES(i)%XS,MESHES(i)%XF,MESHES(i)%YS,MESHES(i)%YF,MESHES(i)%ZS,MESHES(i)%ZF,0._EB,0._EB,0._EB)) THEN
                  PEX%FED_MESH = i
                  EXIT PEX_SmokeLoop
               END IF
            END IF
            !     No mesh found
            PEX%FED_MESH = -1
         END DO PEX_SmokeLoop
         !   No mesh found
         IF (PEX%FED_MESH == 0) PEX%FED_MESH = -1

         IF (PEX%FED_MESH > 0) THEN
            M => MESHES(PEX%FED_MESH)
            II = FLOOR(M%CELLSI(FLOOR((PEX%Xsmoke-M%XS)*M%RDXINT))+ 1.0_EB)
            JJ = FLOOR(M%CELLSJ(FLOOR((PEX%Ysmoke-M%YS)*M%RDYINT))+ 1.0_EB)
            KK = FLOOR(M%CELLSK(FLOOR((PEX%Zsmoke-M%ZS)*M%RDZINT))+ 1.0_EB)
            IF ( M%SOLID(M%CELL_INDEX(II,JJ,KK)) ) THEN
               PEX%FED_MESH = -1   ! no smoke at a solid object
               PEX%II = 0
               PEX%JJ = 0
               PEX%KK = 0
            ELSE
               PEX%II = II
               PEX%JJ = JJ
               PEX%KK = KK
            END IF
         ELSE
            PEX%II = 0
            PEX%JJ = 0
            PEX%KK = 0
         END IF
         !
      END DO READ_EXIT_LOOP
26    REWIND(LU_INPUT)

    END SUBROUTINE READ_EXIT

    SUBROUTINE READ_DOOR
      IMPLICIT NONE
      !
      ! Read the DOOR lines
      !
      ! Local variables
      INTEGER nm, i1, i2, j1, j2, i, ii, iii, N_TMP
      LOGICAL L_TMP
      TYPE (EVAC_DOOR_TYPE), POINTER :: PDX=>NULL()
      TYPE (MESH_TYPE), POINTER :: M=>NULL()
      !
      N_TMP = N_EXITS - N_CO_EXITS
      READ_DOOR_LOOP: DO N = 1, N_DOORS
         !
         ID            = 'null'
         RGB           = -1
         COLOR         = 'null'
         XB            = 0.0_EB
         IOR           = 0
         FLOW_FIELD_ID = 'null'
         VENT_FFIELD   = 'null'
         MESH_ID       = 'null'
         EVAC_MESH     = 'null'
         TO_NODE       = 'null'
         CHECK_FLOW    = .FALSE.
         DEFINE_MESH   = .TRUE.
         EXIT_SIGN     = .TRUE.
         SHOW          = .TRUE.
         MAX_FLOW      = 0.0_EB
         WIDTH         = 0.0_EB
         HEIGHT        = 2.0_EB
         TIME_OPEN     = -HUGE(TIME_OPEN)
         TIME_CLOSE    = HUGE(TIME_CLOSE)
         XYZ(:)        = HUGE(XYZ)
         XYZ_SMOKE(:)  = HUGE(XYZ_SMOKE)
         COLOR_INDEX   = -1
         KEEP_XY       = .TRUE.
         LOCKED_WHEN_CLOSED = .FALSE.
         TARGET_WHEN_CLOSED = .FALSE.
         WAIT_AT_XYZ        = .FALSE.
         KNOWN_DOOR         = .FALSE.
         !
         CALL CHECKREAD('DOOR',LU_INPUT,IOS)  ; IF (STOP_STATUS==SETUP_STOP) RETURN
         IF (IOS == 1) THEN
            EXIT READ_DOOR_LOOP
         END IF
         READ(LU_INPUT,DOOR,END=27,IOSTAT=IOS)
         !
         IF (VENT_FFIELD /= 'null') DEFINE_MESH = .FALSE. ! Old input, user gives the door flow meshes

         IF (TO_NODE == 'null') EXIT_SIGN = .FALSE. ! This is more or less like an entry.
         IF (TO_NODE == 'null') KEEP_XY = .FALSE. ! This is more or less like an entry.
         IF (TO_NODE == 'null') DEFINE_MESH = .FALSE. ! This is more or less like an entry.


         ! Old input used COLOR_INDEX, next lines are needed for that
         IF (MYID==MAX(0,EVAC_PROCESS) .AND. COLOR_INDEX/=-1) WRITE (LU_ERR,'(A,A)') &
              ' WARNING: keyword COLOR_INDEX is replaced by COLOR at DOOR line ',TRIM(ID)
         IF (COLOR_INDEX == 1) COLOR = 'BLACK'
         IF (COLOR_INDEX == 2) COLOR = 'YELLOW'
         IF (COLOR_INDEX == 3) COLOR = 'BLUE'
         IF (COLOR_INDEX == 4) COLOR = 'RED'
         IF (COLOR_INDEX == 5) COLOR = 'GREEN'
         IF (COLOR_INDEX == 6) COLOR = 'MAGENTA'
         IF (COLOR_INDEX == 7) COLOR = 'CYAN'
         !
         ! Colors, integer RGB(3), e.g., (23,255,0)
         IF (ANY(RGB < 0) .AND. COLOR=='null') COLOR = 'FOREST GREEN'
         IF (COLOR /= 'null') CALL COLOR2RGB(RGB,COLOR)
         IF (COLOR_METHOD == 4) THEN
            i_avatar_color = i_avatar_color + 1
            EVAC_AVATAR_RGB(1:3,i_avatar_color) = RGB
         END IF

         IF (MYID /= MAX(0,EVAC_PROCESS)) CYCLE READ_DOOR_LOOP

         PDX=>EVAC_DOORS(N)

         N_TMP = N_TMP + 1 ! real exits + door index (not count_only ones)

         PDX%RGB = RGB
         IF (COLOR_METHOD == 4) PDX%Avatar_Color_Index = i_avatar_color
         IF (TRIM(ID) /= 'null') THEN
            DO I = 1, N-1
               IF (TRIM(ID) == TRIM(EVAC_DOORS(I)%ID)) THEN
                  WRITE(MESSAGE,'(A,I4,A,I4,A,A)') 'ERROR: DOOR lines',I,' and',N,', ID strings are not unique: ',TRIM(ID)
                  CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN

               END IF
            END DO
         END IF

         IF (EVAC_MESH /= 'null') THEN
            MESH_ID = EVAC_MESH
            IF (MYID==MAX(0,EVAC_PROCESS)) WRITE (LU_ERR,'(A,A)') &
                 ' WARNING: keyword EVAC_MESH is replaced by MESH_ID at DOOR line ', TRIM(ID)
         END IF

         IF (TO_NODE=='null') THEN
            EXIT_SIGN = .FALSE.
         END IF

         ! Check that the door is properly specified

         DO I=1,5,2
            IF (XB(I) > XB(I+1)) THEN
               DUMMY   = XB(I)
               XB(I)   = XB(I+1)
               XB(I+1) = DUMMY
            END IF
         END DO
         !
         ! Check which evacuation floor. Now there may be overlapping meshes.
         ii = 0
         PDX_MeshLoop: DO i = 1, NMESHES
            IF (EVACUATION_ONLY(I) .AND. EVACUATION_SKIP(I)) THEN
               IF (Is_Within_Bounds(XB(1),XB(2),XB(3),XB(4),XB(5),XB(6),&
                    MESHES(i)%XS,MESHES(i)%XF,MESHES(i)%YS,MESHES(i)%YF,MESHES(i)%ZS,MESHES(i)%ZF, 0._EB, 0._EB, 0._EB)) THEN
                  IF (TRIM(MESH_ID) == 'null' .OR. TRIM(MESH_ID) == TRIM(MESH_NAME(i))) THEN
                     ii = ii + 1
                     PDX%IMESH = i
                  END IF
                  IF (TRIM(TO_NODE) == TRIM(MESH_NAME(i))) THEN
                     PDX%IMESH2 = i
                  END IF
               END IF
            END IF
         END DO PDX_MeshLoop
         PDX_StrsLoopB: DO I = 1, N_STRS
            IF (Is_Within_Bounds(XB(1), XB(2), XB(3), XB(4), XB(5), XB(6), &
                 EMESH_STAIRS(I)%XB(1), EMESH_STAIRS(I)%XB(2), EMESH_STAIRS(I)%XB(3), &
                 EMESH_STAIRS(I)%XB(4), EMESH_STAIRS(I)%XB(5), EMESH_STAIRS(I)%XB(6), &
                 0._EB, 0._EB, 0._EB)) THEN
               IF (TRIM(MESH_ID) == 'null' .OR. TRIM(MESH_ID) == ('Emesh_' // TRIM(EMESH_STAIRS(I)%ID))) THEN
                  DEFINE_MESH = .FALSE.
               END IF
            END IF
         END DO PDX_StrsLoopB

         IF (DEFINE_MESH) EVAC_DOORS(N)%I_EMESH_EXITS = N_TMP
         PDX_StrsLoop: DO I = 1, N_STRS
            IF (TRIM(PDX%TO_NODE) == TRIM(EVAC_STRS(I)%ID)) THEN
               PDX%IMESH2 = EVAC_STRS(I)%IMESH
               EXIT PDX_StrsLoop
            END IF
         END DO PDX_StrsLoop
         IF (PDX%IMESH == 0) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: DOOR ',TRIM(ID), ' problem with IMESH, no mesh found'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END IF
         IF (ii > 1) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: DOOR ',TRIM(ID), ' not an unique mesh found '
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END IF
         ! Use the main evacuation mesh flow field if none is given
         IF (TRIM(FLOW_FIELD_ID) == 'null') THEN
            FLOW_FIELD_ID = TRIM(MESH_NAME(PDX%IMESH))
         END IF

         nm = PDX%IMESH
         M => MESHES(PDX%IMESH)

         IF (XB(1)/=XB(2) .AND. XB(3)/=XB(4)) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: DOOR ',TRIM(ID),' must be a plane'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         ENDIF

         ! User input
         PDX%X1 = XB(1)
         PDX%X2 = XB(2)
         PDX%Y1 = XB(3)
         PDX%Y2 = XB(4)
         PDX%Z1 = XB(5)
         PDX%Z2 = XB(6)

         ! Move user input to mesh cell boundaries
         XB(1) = MAX(XB(1),M%XS)
         XB(2) = MIN(XB(2),M%XF)
         XB(3) = MAX(XB(3),M%YS)
         XB(4) = MIN(XB(4),M%YF)
         XB(5) = MAX(XB(5),M%ZS)
         XB(6) = MIN(XB(6),M%ZF)

         I1 = NINT( GINV(XB(1)-M%XS,1,nm)*M%RDXI )
         I2 = NINT( GINV(XB(2)-M%XS,1,nm)*M%RDXI )
         J1 = NINT( GINV(XB(3)-M%YS,2,nm)*M%RDETA)
         J2 = NINT( GINV(XB(4)-M%YS,2,nm)*M%RDETA)

         XB(1) = M%X(I1)
         XB(2) = M%X(I2)
         XB(3) = M%Y(J1)
         XB(4) = M%Y(J2)

         IF ( ABS(XB(1)-PDX%X1)>1.E-4_EB .OR. ABS(XB(2)-PDX%X2)>1.E-4_EB .OR. &
              ABS(XB(3)-PDX%Y1)>1.E-4_EB .OR. ABS(XB(4)-PDX%Y2)>1.E-4_EB ) THEN
            WRITE(LU_ERR,fmt='(a,a,a,a)') ' WARNING: Door ',TRIM(ID),' XB adjusted to mesh ',TRIM(MESH_NAME(nm))
            WRITE(LU_ERR,fmt='(a,6f12.4)') ' Old XB:', PDX%X1,PDX%X2,PDX%Y1,PDX%Y2,PDX%Z1,PDX%Z2
            WRITE(LU_ERR,fmt='(a,6f12.4)') ' New XB:', XB(1:6)
         END IF

         II = (I1+I2)/2
         JJ = (J1+J2)/2
         SELECT CASE (IOR)
         CASE (-1)
            II = II + 1
         CASE (-2)
            JJ = JJ + 1
         END SELECT
         IF (M%SOLID(M%CELL_INDEX(II,JJ,1))) THEN
            WRITE(LU_ERR,fmt='(a,a,a)') ' WARNING: Door ',TRIM(ID),' problem with XB, mid point facing solid'
         END IF

         ! Coordinates are lined up with the mesh.
         PDX%X1 = XB(1)
         PDX%X2 = XB(2)
         PDX%Y1 = XB(3)
         PDX%Y2 = XB(4)
         PDX%Z1 = XB(5)
         PDX%Z2 = XB(6)
         !
         ! Find the implicitely generated mesh behind this exit/door
         I = PDX%I_EMESH_EXITS
         IF (VENT_FFIELD == 'null' .AND. I > 0) THEN
            IF (EMESH_EXITS(I)%NEXIT == N .AND. .NOT.EMESH_EXITS(I)%IS_EXIT) THEN
               VENT_FFIELD = TRIM(MESH_NAME(EMESH_EXITS(I)%IMESH))
            END IF
         END IF

         IF (TO_NODE == 'null') THEN
            TO_NODE = TRIM(ID)
         END IF

         PDX%IOR        = IOR
         PDX%ID         = Trim(ID)
         PDX%GRID_NAME  = Trim(FLOW_FIELD_ID)
         PDX%VENT_FFIELD= Trim(VENT_FFIELD)
         PDX%CHECK_FLOW = CHECK_FLOW
         PDX%EXIT_SIGN  = EXIT_SIGN
         PDX%KEEP_XY    = KEEP_XY
         PDX%SHOW       = SHOW
         PDX%TO_NODE    = Trim(TO_NODE)
         PDX%INODE      = 0
         PDX%INODE2     = 0
         PDX%T_first    = T_BEGIN
         PDX%T_last     = T_BEGIN
         PDX%ICOUNT     = 0
         PDX%Flow_max   = 0.0_EB
         PDX%TIME_OPEN  = TIME_OPEN
         PDX%TIME_CLOSE = TIME_CLOSE
         PDX%LOCKED_WHEN_CLOSED = LOCKED_WHEN_CLOSED
         PDX%TARGET_WHEN_CLOSED = TARGET_WHEN_CLOSED
         PDX%WAIT_AT_XYZ        = WAIT_AT_XYZ
         PDX%ICORR_TO   = 0  ! The corr index if target is an elevator corr, otherwise zero
         PDX%IMODE      = -1 ! Door is open by default
         PDX%KNOWN_DOOR = KNOWN_DOOR
         IF (TIME_OPEN > TIME_CLOSE) THEN
            PDX%IMODE = -1
         ELSE
            IF (TIME_OPEN > T_BEGIN) PDX%IMODE = +2
         END IF
         IF (CHECK_FLOW) PDX%Flow_max   = MAX_FLOW
         PDX%HEIGHT = HEIGHT


         SELECT CASE (IOR)
         CASE (-1,+1)
            IF (WIDTH <= 0.0_EB) THEN
               PDX%Width = XB(4) - XB(3)
            ELSE
               PDX%Width = WIDTH
            END IF
            PDX%ORIENTATION(1)=REAL(SIGN(1,IOR),EB)
         CASE (-2,+2)
            IF (WIDTH <= 0.0_EB) THEN
               PDX%Width = XB(2) - XB(1)
            ELSE
               PDX%Width = WIDTH
            END IF
            PDX%ORIENTATION(2)=REAL(SIGN(1,IOR),EB)
         CASE (-3)
            IF ( (XB(4)-XB(3)) <= 0.0_EB .OR. (XB(2)-XB(1)) <= 0.0_EB) THEN
               WRITE(MESSAGE,'(A,A,A)') 'ERROR: DOOR ',TRIM(ID),' IOR=-3 but not 3-dim object'
               CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
            END IF
            PDX%ORIENTATION(3)=REAL(SIGN(1,IOR),EB)
         CASE (0)
            IF ( (XB(4)-XB(3)) <= 0.0_EB .OR. (XB(2)-XB(1)) <= 0.0_EB) THEN
               WRITE(MESSAGE,'(A,A,A)') 'ERROR: DOOR ',TRIM(ID),' no IOR but not 3-dim object'
               CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
            END IF
         CASE Default
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: DOOR ',TRIM(ID),' problem with IOR'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END SELECT
         !

         L_TMP=.FALSE.
         DO i = 1, NMESHES
            IF (.NOT. EVACUATION_ONLY(i)) CYCLE
            IF (TRIM(FLOW_FIELD_ID)==TRIM(MESH_NAME(i))) THEN
               L_TMP=.TRUE.
               EXIT
            END IF
         END DO
         IF (.NOT.(TRIM(FLOW_FIELD_ID)=='null' .OR. L_TMP)) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: DOOR ',TRIM(ID),' problem with FLOW_FIELD_ID'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END IF
         L_TMP=.FALSE.
         DO i = 1, NMESHES
            IF (.NOT. EVACUATION_ONLY(i)) CYCLE
            IF (TRIM(VENT_FFIELD)==TRIM(MESH_NAME(i))) THEN
               L_TMP=.TRUE.
               EXIT
            END IF
         END DO
         IF (.NOT.(TRIM(VENT_FFIELD)=='null' .OR. L_TMP)) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: DOOR ',TRIM(ID),' problem with VENT_FFIELD'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END IF

         ! Check which vent field. If VENT_FFIELD is not found, use the main evac grid.
         PDX%I_VENT_FFIELD = 0
         PDX_Mesh2Loop: DO i = 1, NMESHES
            IF ( EVACUATION_ONLY(I) .AND. (TRIM(MESH_NAME(i)) == TRIM(PDX%VENT_FFIELD)) ) THEN
               IF ( (PDX%Z1 >= MESHES(i)%ZS .AND. PDX%Z2 <= MESHES(i)%ZF).AND. &
                    (PDX%Y1 >= MESHES(i)%YS .AND. PDX%Y2 <= MESHES(i)%YF).AND. &
                    (PDX%X1 >= MESHES(i)%XS .AND. PDX%X2 <= MESHES(i)%XF)) THEN
                  PDX%I_VENT_FFIELD = i
                  EXIT PDX_Mesh2Loop
               END IF
            END IF
         END DO PDX_Mesh2Loop
         ! If no vent field is given, then use the main evac grid.
         IF (PDX%I_VENT_FFIELD == 0) THEN
            PDX%I_VENT_FFIELD = PDX%IMESH
            PDX%VENT_FFIELD = TRIM(MESH_NAME(PDX%IMESH))
         END IF

         PDX%FED_MESH = 0
         IF (XYZ(1) < HUGE(XYZ)) THEN
            PDX%X = XYZ(1)
            PDX%Y = XYZ(2)
            PDX%Z = 0.5_EB*(XB(5)+XB(6))
         ELSE
            PDX%X = 0.5_EB*(XB(1)+XB(2))
            PDX%Y = 0.5_EB*(XB(3)+XB(4))
            PDX%Z = 0.5_EB*(XB(5)+XB(6))
            SELECT CASE (IOR)
            CASE (-1,+1)
               PDX%X = PDX%X - IOR*0.20_EB
            CASE (-2,+2)
               PDX%Y = PDX%Y - IOR*0.10_EB
            END SELECT
         END IF
         !
         ! Check which evacuation floor
         ii = 0
         iii = 0
         PDX_Mesh3Loop: DO i = 1, NMESHES
            IF (EVACUATION_ONLY(I) .AND. EVACUATION_SKIP(I)) THEN
               IF ( (PDX%Z >= MESHES(i)%ZS .AND. PDX%Z <= MESHES(i)%ZF).AND. &
                    (PDX%Y >= MESHES(i)%YS .AND. PDX%Y <= MESHES(i)%YF).AND. &
                    (PDX%X >= MESHES(i)%XS .AND. PDX%X <= MESHES(i)%XF)) THEN
                  IF (TRIM(MESH_ID) == 'null' .OR. TRIM(MESH_ID) == TRIM(MESH_NAME(i))) THEN
                     ii  = ii + 1
                     iii = i
                  END IF
               END IF
            END IF
         END DO PDX_Mesh3Loop
         IF (ii == 0) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: DOOR line ',TRIM(PDX%ID), ' problem with XYZ, no mesh found'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END IF
         IF (ii > PDX%IMESH) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: DOOR line ',TRIM(PDX%ID), ' problem with XYZ, not an unique mesh found'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END IF
         IF (iii /= PDX%IMESH .AND. SHOW) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: DOOR line ',TRIM(PDX%ID), ' problem with XYZ, wrong mesh found'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END IF
         M => MESHES(PDX%IMESH)
         II = FLOOR(M%CELLSI(FLOOR((PDX%X-M%XS)*M%RDXINT))+ 1.0_EB)
         JJ = FLOOR(M%CELLSJ(FLOOR((PDX%Y-M%YS)*M%RDYINT))+ 1.0_EB)
         KK = FLOOR(M%CELLSK(FLOOR((PDX%Z-M%ZS)*M%RDZINT))+ 1.0_EB)
         IF (M%SOLID(M%CELL_INDEX(II,JJ,KK))) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: DOOR line ',TRIM(PDX%ID), ' problem with XYZ, inside solid'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END IF

         ! PDX%Z is used to plot the door on the correct height in Smokeview.
         PDX%Z = PDX%Z + 0.5_EB*PDX%Height - EVACUATION_Z_OFFSET(PDX%IMESH)
         PDX%XYZ_Z = PDX%Z

         IF (XYZ_SMOKE(1) < HUGE(XYZ_SMOKE)) THEN
            PDX%Xsmoke = XYZ_SMOKE(1)
            PDX%Ysmoke = XYZ_SMOKE(2)
            PDX%Zsmoke = XYZ_SMOKE(3)
         ELSE
            PDX%Xsmoke = PDX%X
            PDX%Ysmoke = PDX%Y
            PDX%Zsmoke = 0.5_EB*(XB(5)+XB(6)) - EVACUATION_Z_OFFSET(PDX%IMESH) + HUMAN_SMOKE_HEIGHT
         END IF

         ! Check, which fire grid and i,j,k (xyz)
         PDX_SmokeLoop: DO i = 1, NMESHES
            IF (.NOT. EVACUATION_ONLY(I)) THEN
               IF ( (PDX%Zsmoke >= MESHES(i)%ZS .AND. PDX%Zsmoke <= MESHES(i)%ZF) .AND. &
                    (PDX%Ysmoke >= MESHES(i)%YS .AND. PDX%Ysmoke <= MESHES(i)%YF) .AND. &
                    (PDX%Xsmoke >= MESHES(i)%XS .AND. PDX%Xsmoke <= MESHES(i)%XF)) THEN
                  PDX%FED_MESH = i
                  EXIT PDX_SmokeLoop
               END IF
            END IF
            !     No mesh found
            PDX%FED_MESH = -1
         END DO PDX_SmokeLoop
         !   No mesh found
         IF (PDX%FED_MESH == 0) PDX%FED_MESH = -1

         IF (PDX%FED_MESH > 0) THEN
            M => MESHES(PDX%FED_MESH)
            II = FLOOR(M%CELLSI(FLOOR((PDX%Xsmoke-M%XS)*M%RDXINT))+ 1.0_EB)
            JJ = FLOOR(M%CELLSJ(FLOOR((PDX%Ysmoke-M%YS)*M%RDYINT))+ 1.0_EB)
            KK = FLOOR(M%CELLSK(FLOOR((PDX%Zsmoke-M%ZS)*M%RDZINT))+ 1.0_EB)
            IF ( M%SOLID(M%CELL_INDEX(II,JJ,KK)) ) THEN
               PDX%FED_MESH = -1   ! no smoke at a solid object
               PDX%II = 0
               PDX%JJ = 0
               PDX%KK = 0
            ELSE
               PDX%II = II
               PDX%JJ = JJ
               PDX%KK = KK
            END IF
         ELSE
            PDX%II = 0
            PDX%JJ = 0
            PDX%KK = 0
         END IF
         !
      END DO READ_DOOR_LOOP
27    REWIND(LU_INPUT)

    END SUBROUTINE READ_DOOR

    SUBROUTINE READ_CORR
      IMPLICIT NONE
      !
      ! Local variables
      TYPE (EVAC_CORR_TYPE), POINTER :: PCX=>NULL()
      TYPE (MESH_TYPE), POINTER :: M=>NULL()
      !
      ! Read the CORR line
      !
      n_max_in_corrs = 0
      READ_CORR_LOOP: DO N = 1, N_CORRS
         IF (MYID /= MAX(0,EVAC_PROCESS)) CYCLE READ_CORR_LOOP
         PCX=>EVAC_CORRS(N)
         !
         ID            = 'null'
         RGB           = -1
         COLOR         = 'null'
         XB            = HUGE(XB)
         XB1           = HUGE(XB1)
         XB2           = HUGE(XB2)
         IOR           = 0
         FLOW_FIELD_ID = 'null'
         TO_NODE       = 'null'
         CHECK_FLOW    = .FALSE.
         ELEVATOR      = .FALSE.
         MAX_FLOW      = 0.0_EB
         WIDTH         = 0.0_EB
         WIDTH1        = 0.0_EB
         WIDTH2        = 0.0_EB
         FAC_SPEED     = 0.0_EB
         EFF_WIDTH     = 0.0_EB
         EFF_LENGTH    = 0.0_EB
         TRAVEL_TIME   = -.1_EB
         MAX_HUMANS_INSIDE = 0
         !
         CALL CHECKREAD('CORR',LU_INPUT,IOS)  ; IF (STOP_STATUS==SETUP_STOP) RETURN
         IF (IOS == 1) THEN
            EXIT READ_CORR_LOOP
         END IF
         READ(LU_INPUT,CORR,END=29,IOSTAT=IOS)
         !
         IF (TRIM(ID) /= 'null') THEN
            DO I = 1, N-1
               IF (TRIM(ID) == TRIM(EVAC_CORRS(I)%ID)) THEN
                  WRITE(MESSAGE,'(A,I4,A,I4,A,A)') 'ERROR: CORR lines',I,' and',N,', ID strings are not unique: ',TRIM(ID)
                  CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
               END IF
            END DO
         END IF
         !
         DO I=1,5,2
            IF (XB(I) > XB(I+1)) THEN
               DUMMY   = XB(I)
               XB(I)   = XB(I+1)
               XB(I+1) = DUMMY
            END IF
         END DO
         DO I=1,5,2
            IF (XB1(I) > XB1(I+1)) THEN
               DUMMY   = XB1(I)
               XB1(I)   = XB1(I+1)
               XB1(I+1) = DUMMY
            END IF
         END DO
         DO I=1,5,2
            IF (XB2(I) > XB2(I+1)) THEN
               DUMMY   = XB2(I)
               XB2(I)   = XB2(I+1)
               XB2(I+1) = DUMMY
            END IF
         END DO
         !
         ! Position, where smoke etc. is saved. If both XB and XB1 are given, use XB1
         IF ( XB(1) < HUGE(XB) ) THEN
            PCX%FED_MESH = 0
            PCX%X1 = 0.5_EB*( XB(1) +  XB(2))
            PCX%Y1 = 0.5_EB*( XB(3) +  XB(4))
            PCX%Z1 = 0.5_EB*( XB(5) +  XB(6))
         ELSE
            PCX%FED_MESH = -1
            PCX%X1 = 0.0_EB
            PCX%Y1 = 0.0_EB
            PCX%Z1 = 0.0_EB
         END IF
         IF ( XB1(1) < HUGE(XB1) ) THEN
            PCX%FED_MESH = 0
            PCX%X1 = 0.5_EB*( XB1(1) +  XB1(2))
            PCX%Y1 = 0.5_EB*( XB1(3) +  XB1(4))
            PCX%Z1 = 0.5_EB*( XB1(5) +  XB1(6))
         ELSE IF (XB(1) == HUGE(XB) ) THEN
            PCX%FED_MESH = -1
            PCX%X1 = 0.0_EB
            PCX%Y1 = 0.0_EB
            PCX%Z1 = 0.0_EB
         END IF
         IF ( XB2(1) < HUGE(XB2) ) THEN
            PCX%FED_MESH2 = 0
            PCX%X2 = 0.5_EB*(XB2(1) + XB2(2))
            PCX%Y2 = 0.5_EB*(XB2(3) + XB2(4))
            PCX%Z2 = 0.5_EB*(XB2(5) + XB2(6))
         ELSE
            PCX%FED_MESH2 = -1
            PCX%X2 = 0.0_EB
            PCX%Y2 = 0.0_EB
            PCX%Z2 = 0.0_EB
         END IF

         PCX%IOR        = IOR
         PCX%ID         = ID
         PCX%GRID_NAME  = FLOW_FIELD_ID
         PCX%CHECK_FLOW = CHECK_FLOW
         PCX%TO_NODE    = TO_NODE
         PCX%INODE      = 0
         PCX%INODE2     = 0
         PCX%T_first    = T_BEGIN
         PCX%T_last     = T_BEGIN
         PCX%ICOUNT     = 0
         PCX%IDOOR_FROM = 0 ! "Top end door", if exits (could also be corr==>coor, not door==>corr)

         PCX%MAX_HUMANS_INSIDE = 0
         IF (MAX_HUMANS_INSIDE > 0 ) THEN
            PCX%MAX_HUMANS_INSIDE = MAX_HUMANS_INSIDE
         ELSE
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: CORR ',TRIM(ID),' MAX_HUMANS_INSIDE <= 0'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END IF

         IF (ELEVATOR) THEN
            PCX%ELEVATOR = .TRUE.
            IF (TRAVEL_TIME < 0) THEN
               WRITE(MESSAGE,'(A,A,A)') 'ERROR: CORR ',TRIM(ID),' TRAVEL_TIME < 0 for an elevator'
               CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
            END IF
            PCX%Travel_Time = TRAVEL_TIME  ! up-down and down-up times supposed to be the same
            PCX%Time_Wait_To_Fill = 5.0_EB ! Wait 5 s if there still empty places in the elevator
            PCX%Time_Last_In = T_BEGIN     ! Time stamp of the last person entered the elevator car
            PCX%Elevator_Mode = 2          ! 0 moving down, 1 down, 2 up, 3 moving up
            PCX%Time_Stop = 0.0_EB         ! Time when the car stopped (and doors open)
            PCX%Time_Start = 0.0_EB        ! Time when the car closed the doors
            FAC_SPEED = 1.0_EB
         END IF


         IF (FAC_SPEED < 0) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: CORR ',TRIM(ID),' FAC_SPEED < 0'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         ELSE
            IF (FAC_SPEED == 0.0_EB) FAC_SPEED = 0.6_EB
            PCX%Fac_Speed = FAC_SPEED
         END IF

         PCX%Flow_max   = 0.0_EB
         IF (CHECK_FLOW) PCX%Flow_max   = MAX_FLOW

         PCX%Width = MAX( ABS(XB(4)-XB(3)) , ABS(XB(2)-XB(1)) )

         PCX%Eff_Width = 0.0_EB
         IF (EFF_WIDTH >= TWO_EPSILON_EB ) THEN
            PCX%Eff_Width = EFF_WIDTH
         ELSE
            PCX%Eff_Width = PCX%Width
         END IF

         PCX%Eff_Length = 0.0_EB
         IF (EFF_LENGTH >= TWO_EPSILON_EB ) THEN
            PCX%Eff_Length = EFF_LENGTH
         ELSE
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: CORR ',TRIM(PCX%ID),' EFF_LENGTH <= 0'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END IF
         PCX%Eff_Area = PCX%Eff_Length*PCX%Eff_Width

         PCX%Width1 = WIDTH1
         PCX%Width2 = WIDTH2
         IF (WIDTH1*WIDTH2 <= 0.0_EB) THEN
            PCX%Width1 = PCX%Width
            PCX%Width2 = PCX%Width
         END IF
         !
         PCX_MeshLoop: DO i = 1, NMESHES
            IF (.NOT. EVACUATION_ONLY(I) .AND. PCX%FED_MESH >= 0) THEN
               IF ( (PCX%Z1 >= MESHES(i)%ZS .AND. PCX%Z1 <= MESHES(i)%ZF).AND. &
                    (PCX%Y1 >= MESHES(i)%YS .AND. PCX%Y1 <= MESHES(i)%YF).AND. &
                    (PCX%X1 >= MESHES(i)%XS .AND. PCX%X1 <= MESHES(i)%XF)) THEN
                  PCX%FED_MESH = i
                  EXIT PCX_MeshLoop
               END IF
            END IF
            !     No mesh found
            PCX%FED_MESH = -1
         END DO PCX_MeshLoop
         !   No mesh found
         IF (PCX%FED_MESH == 0) PCX%FED_MESH = -1

         IF (PCX%FED_MESH > 0) THEN
            M => MESHES(PCX%FED_MESH)
            II = FLOOR( M%CELLSI(FLOOR((PCX%X1-M%XS)*M%RDXINT)) + 1.0_EB  )
            JJ = FLOOR( M%CELLSJ(FLOOR((PCX%Y1-M%YS)*M%RDYINT)) + 1.0_EB  )
            KK = FLOOR( M%CELLSK(FLOOR((PCX%Z1-M%ZS)*M%RDZINT)) + 1.0_EB  )
            IF ( M%SOLID(M%CELL_INDEX(II,JJ,KK)) ) THEN
               PCX%FED_MESH = -1   ! no smoke at a solid object
               PCX%II(1) = 0
               PCX%JJ(1) = 0
               PCX%KK(1) = 0
            ELSE
               PCX%II(1) = II
               PCX%JJ(1) = JJ
               PCX%KK(1) = KK
            END IF
         ELSE
            PCX%II(1) = 0
            PCX%JJ(1) = 0
            PCX%KK(1) = 0
         END IF

         PCX_MeshLoop2: DO i = 1, NMESHES
            IF (.NOT. EVACUATION_ONLY(I) .AND. PCX%FED_MESH2 >= 0) THEN
               IF ( (PCX%Z2 >= MESHES(i)%ZS .AND. PCX%Z2 <= MESHES(i)%ZF).AND. &
                    (PCX%Y2 >= MESHES(i)%YS .AND. PCX%Y2 <= MESHES(i)%YF).AND. &
                    (PCX%X2 >= MESHES(i)%XS .AND. PCX%X2 <= MESHES(i)%XF)) THEN
                  PCX%FED_MESH2 = i
                  EXIT PCX_MeshLoop2
               END IF
            END IF
         END DO PCX_MeshLoop2
         !   No mesh found
         IF (PCX%FED_MESH2 == 0) PCX%FED_MESH2 = -1

         IF (PCX%FED_MESH2 > 0) THEN
            M => MESHES(PCX%FED_MESH2)
            II = FLOOR( M%CELLSI(FLOOR((PCX%X2-M%XS)*M%RDXINT)) + 1.0_EB  )
            JJ = FLOOR( M%CELLSJ(FLOOR((PCX%Y2-M%YS)*M%RDYINT)) + 1.0_EB  )
            KK = FLOOR( M%CELLSK(FLOOR((PCX%Z2-M%ZS)*M%RDZINT)) + 1.0_EB  )
            IF ( M%SOLID(M%CELL_INDEX(II,JJ,KK)) ) THEN
               PCX%FED_MESH2 = -1   ! no smoke at a solid object
               PCX%II(2) = 0
               PCX%JJ(2) = 0
               PCX%KK(2) = 0
            ELSE
               PCX%II(2) = II
               PCX%JJ(2) = JJ
               PCX%KK(2) = KK
            END IF
         ELSE
            PCX%II(2) = 0
            PCX%JJ(2) = 0
            PCX%KK(2) = 0
         END IF
         !
         n_max_in_corrs = MAX(n_max_in_corrs,PCX%MAX_HUMANS_INSIDE)

         ! Initialize the linked lists of persons who are inside corridors.
         PCX%n_inside = 0
         NULLIFY(PCX%First)

         ! Colors, integer RGB(3), e.g., (23,255,0)
         IF (ANY(RGB < 0) .AND. COLOR=='null') COLOR = 'BLACK'
         IF (COLOR /= 'null') CALL COLOR2RGB(RGB,COLOR)
         PCX%RGB = RGB
         !
      END DO READ_CORR_LOOP
29    REWIND(LU_INPUT)

    END SUBROUTINE READ_CORR

    SUBROUTINE READ_STRS
      IMPLICIT NONE
      !
      ! Local variables
      REAL(EB) Z_TMP
      TYPE (EVAC_STRS_TYPE), POINTER :: STRP=>NULL()

      ! Read the STRS line
      READ_STRS_LOOP: DO N = 1,N_STRS
         IF (MYID /= MAX(0,EVAC_PROCESS)) CYCLE READ_STRS_LOOP
         STRP=>EVAC_STRS(N)
         !
         ID                          = 'null'
         XB                          = 0._EB
         XB_CORE                     = 0._EB
         XB_CORES                    = 0._EB
         XB_LANDINGS                 = 0._EB
         RIGHT_HANDED                = .TRUE.
         LEFT_HANDED                 = .FALSE.
         MESH_ID                     = 'null'
         N_LANDINGS                  = 0
         VERTICAL_LANDING_SEPARATION = 0._EB
         FAC_V0_UP                   = 1.0_EB
         FAC_V0_DOWN                 = 1.0_EB
         FAC_V0_HORI                 = 1.0_EB
         !
         CALL CHECKREAD('STRS',LU_INPUT,IOS)  ; IF (STOP_STATUS==SETUP_STOP) RETURN
         IF (IOS == 1) THEN
            EXIT READ_STRS_LOOP
         END IF
         READ(LU_INPUT,STRS,END=32,IOSTAT=IOS)
         !
         DO I=1,5,2
            IF (XB(I) > XB(I+1)) THEN
               DUMMY   = XB(I)
               XB(I)   = XB(I+1)
               XB(I+1) = DUMMY
            END IF
         END DO
         !
         STRP%ID          = ID
         STRP%XB          = XB
         ii = 0
         IF (TRIM(MESH_ID) == 'null') THEN
            MESH_ID = 'Emesh_' // TRIM(ID)
         END IF
         STRP_MeshLoop: DO I = 1, NMESHES
            IF (.NOT. EVACUATION_ONLY(I)) CYCLE
            IF (.NOT. EVACUATION_SKIP(I)) CYCLE
            IF (TRIM(MESH_ID) == 'null' .OR. TRIM(MESH_ID)==TRIM(MESH_NAME(I))) THEN
               ii = ii + 1
               STRP%IMESH = I
               EXIT STRP_MeshLoop
            END IF
         END DO STRP_MeshLoop
         IF (ii == 0) THEN
            WRITE(MESSAGE,'(A,A,A,A,A)') 'ERROR: STRS ',TRIM(STRP%ID),' no mesh ',TRIM(MESH_ID),' found'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END IF
         IF (ii > 1) THEN
            WRITE(MESSAGE,'(A,A,A,A,A)') 'ERROR: STRS ',TRIM(STRP%ID),' no unique mesh ',TRIM(MESH_ID),' found'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END IF

         ! Count number of cores
         STRP%N_CORES = 0
         DO I = 1,500
            IF (ANY(XB_CORES(I,:)/=0._EB)) STRP%N_CORES = STRP%N_CORES + 1
         ENDDO

         STRP%INODE          = 0
         IF (LEFT_HANDED) RIGHT_HANDED = .FALSE.
         STRP%RIGHT_HANDED = RIGHT_HANDED
         STRP%MESH_ID     = MESH_ID
         STRP%FAC_V0_UP   = FAC_V0_UP
         STRP%FAC_V0_DOWN = FAC_V0_DOWN
         STRP%FAC_V0_HORI = FAC_V0_HORI

         IF (N_LANDINGS>500) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: STRS ',TRIM(STRP%ID),' N_LANDINGS > 500'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END IF
         STRP%N_LANDINGS = N_LANDINGS
         STRP%N_NODES = 2*N_LANDINGS - 1

         ! Allocate landing and stair geometry
         ALLOCATE(STRP%XB_NODE(1:STRP%N_NODES,1:8), STAT=IZERO)
         CALL ChkMemErr('Read_Evac','STRP%XB_NODE',IZERO)
         ALLOCATE(STRP%NODE_IOR(1:STRP%N_NODES), STAT=IZERO)
         CALL ChkMemErr('Read_Evac','STRP%NODE_IOR',IZERO)
         ALLOCATE(STRP%NODE_TYPE(1:STRP%N_NODES), STAT=IZERO)
         CALL ChkMemErr('Read_Evac','STRP%NODE_TYPE',IZERO)
         STRP%XB_NODE = 0._EB
         STRP%NODE_IOR = 0

         ALLOCATE(STRP%XB_CORE(1:MAX(1,STRP%N_CORES),1:6), STAT=IZERO)
         CALL ChkMemErr('Read_Evac','STRP%XB_CORE',IZERO)
         IF (STRP%N_CORES == 0) THEN
            STRP%N_CORES = 1
            STRP%XB_CORE(1,1:4) = XB_CORE(1:4)
            STRP%XB_CORE(1,5:6) = XB(5:6)
            IF (ALL(XB_CORE==0._EB)) THEN
               WRITE(MESSAGE,'(3A)') 'ERROR: STRS object ', TRIM(ID), ' has no XB_CORE defined.'
               CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
            ENDIF

         ELSE
            STRP%XB_CORE(1:STRP%N_CORES,1:6) = XB_CORES(1:STRP%N_CORES,1:6)
         ENDIF
         ALLOCATE(STRP%I_CORE(1:STRP%N_NODES), STAT=IZERO)
         CALL ChkMemErr('Read_Evac','STRP%I_CORE',IZERO)

         ! Count and copy explicitly given landing geometries
         NL = 0
         DO I = 1,N_LANDINGS
            IF (ANY(XB_LANDINGS(I,:)/=0._EB)) NL = NL + 1
         END DO
         DO I = NL+1,N_LANDINGS
            XB_LANDINGS(I,1:4) = XB_LANDINGS(I-2,1:4)
            XB_LANDINGS(I,5:6) = XB_LANDINGS(I-1,5:6)+VERTICAL_LANDING_SEPARATION
         END DO
         ! Compute stair geometry
         DO I = 1,N_LANDINGS-1
            XB_STAIRS(I,5) = 0.5_EB*(XB_LANDINGS(I,5)+XB_LANDINGS(I,6))
            XB_STAIRS(I,6) = 0.5_EB*(XB_LANDINGS(I+1,5)+XB_LANDINGS(I+1,6))
            STR_Height = XB_STAIRS(I,6)-XB_STAIRS(I,5)
            Z_TMP = 0.5_EB * (XB_STAIRS(I,6)+XB_STAIRS(I,5))
            ! Choose core
            DO J = 1,STRP%N_CORES
               IF ((Z_TMP >= STRP%XB_CORE(J,5)) .AND. (Z_TMP <= STRP%XB_CORE(J,6))) THEN
                  EXIT
               ENDIF
            ENDDO
            J = MIN(J,STRP%N_CORES)
            IF (XB_LANDINGS(I+1,1)>XB_LANDINGS(I,2)) THEN ! From -x to +x
               XB_STAIRS(I,1) = XB_LANDINGS(I,2)
               XB_STAIRS(I,2) = XB_LANDINGS(I+1,1)
               STR_Length = XB_STAIRS(I,2)-XB_STAIRS(I,1)
               IF (RIGHT_HANDED) THEN
                  XB_STAIRS(I,3) = STRP%XB(3)
                  XB_STAIRS(I,4) = STRP%XB_CORE(J,3)
               ELSE
                  XB_STAIRS(I,3) = STRP%XB_CORE(J,4)
                  XB_STAIRS(I,4) = STRP%XB(4)
               END IF
               STRP%NODE_IOR(2*I) = +1
               IF (STR_Length >= TWO_EPSILON_EB) XB_STAIRS(I,7) = COS(ATAN(STR_Height/STR_Length))
               XB_STAIRS(I,8) = 1._EB
            ELSE IF (XB_LANDINGS(I+1,2)<XB_LANDINGS(I,1)) THEN ! From +x to -x
               XB_STAIRS(I,1) = XB_LANDINGS(I+1,2)
               XB_STAIRS(I,2) = XB_LANDINGS(I,1)
               STR_Length = XB_STAIRS(I,2)-XB_STAIRS(I,1)
               IF (RIGHT_HANDED) THEN
                  XB_STAIRS(I,3) = STRP%XB_CORE(J,4)
                  XB_STAIRS(I,4) = STRP%XB(4)
               ELSE
                  XB_STAIRS(I,3) = STRP%XB(3)
                  XB_STAIRS(I,4) = STRP%XB_CORE(J,3)
               END IF
               STRP%NODE_IOR(2*I) = -1
               IF (STR_Length >= TWO_EPSILON_EB) XB_STAIRS(I,7) = COS(ATAN(STR_Height/STR_Length))
               XB_STAIRS(I,8) = 1._EB
            END IF
            IF (XB_LANDINGS(I+1,3)>XB_LANDINGS(I,4)) THEN ! From -y to +y
               IF (RIGHT_HANDED) THEN
                  XB_STAIRS(I,1) = STRP%XB_CORE(J,2)
                  XB_STAIRS(I,2) = STRP%XB(2)
               ELSE
                  XB_STAIRS(I,1) = STRP%XB(1)
                  XB_STAIRS(I,2) = STRP%XB_CORE(J,1)
               END IF
               XB_STAIRS(I,3) = XB_LANDINGS(I,4)
               XB_STAIRS(I,4) = XB_LANDINGS(I+1,3)
               STR_Length = XB_STAIRS(I,4)-XB_STAIRS(I,3)
               STRP%NODE_IOR(2*I) = +2
               XB_STAIRS(I,7) = 1._EB
               IF (STR_Length >= TWO_EPSILON_EB) XB_STAIRS(I,8) = COS(ATAN(STR_Height/STR_Length))
            ELSE IF (XB_LANDINGS(I+1,4)<XB_LANDINGS(I,3)) THEN ! From +y to -y
               IF (RIGHT_HANDED) THEN
                  XB_STAIRS(I,1) = STRP%XB(1)
                  XB_STAIRS(I,2) = STRP%XB_CORE(J,1)
               ELSE
                  XB_STAIRS(I,1) = STRP%XB_CORE(J,2)
                  XB_STAIRS(I,2) = STRP%XB(2)
               END IF
               XB_STAIRS(I,3) = XB_LANDINGS(I+1,4)
               XB_STAIRS(I,4) = XB_LANDINGS(I,3)
               STR_Length = XB_STAIRS(I,4)-XB_STAIRS(I,3)
               STRP%NODE_IOR(2*I) = -2
               XB_STAIRS(I,7) = 1._EB
               IF (STR_Length >= TWO_EPSILON_EB) XB_STAIRS(I,8) = COS(ATAN(STR_Height/STR_Length))
            END IF
         END DO

         ! Collect sub-node coordinates
         J = 1
         DO I = 1, N_LANDINGS
            STRP%XB_NODE(J,1:6)     = XB_LANDINGS(I,:)
            STRP%NODE_TYPE(J)       = STRS_LANDING_TYPE
            ! Choose core
            Z_TMP = 0.5_EB * (STRP%XB_NODE(J,5)+STRP%XB_NODE(J,6))
            DO K = 1,STRP%N_CORES
               IF ((Z_TMP >= STRP%XB_CORE(K,5)) .AND. (Z_TMP <= STRP%XB_CORE(K,6))) THEN
                  EXIT
               ENDIF
            ENDDO
            K = MIN(K,STRP%N_CORES)
            STRP%I_CORE(J) = K
            IF (I<N_LANDINGS) THEN
               STRP%XB_NODE(J+1,1:8) = XB_STAIRS(I,:)
               STRP%NODE_TYPE(J+1)   = STRS_STAIR_TYPE
               ! Choose core
               Z_TMP = 0.5_EB * (STRP%XB_NODE(J+1,5)+STRP%XB_NODE(J+1,6))
               DO K = 1,STRP%N_CORES
                  IF ((Z_TMP >= STRP%XB_CORE(K,5)) .AND. (Z_TMP <= STRP%XB_CORE(K,6))) THEN
                     EXIT
                  ENDIF
               ENDDO
               K = MIN(K,STRP%N_CORES)
               STRP%I_CORE(J+1) = K
            END IF
            J = J + 2
         END DO

      END DO READ_STRS_LOOP
32    REWIND(LU_INPUT)

    END SUBROUTINE READ_STRS

    LOGICAL FUNCTION Is_Within_Bounds(P1x1,P1x2,P1y1,P1y2,P1z1,P1z2,&
         P2x1,P2x2,P2y1,P2y2,P2z1,P2z2,xtol,ytol,ztol)
      IMPLICIT NONE
      !
      REAL(EB), INTENT(IN) :: P1x1,P1x2,P1y1,P1y2,P1z1,P1z2
      REAL(EB), INTENT(IN) :: P2x1,P2x2,P2y1,P2y2,P2z1,P2z2,xtol,ytol,ztol
      !Returns .TRUE. if P2 is within the bounds of P1 with tolerances.
      Is_Within_Bounds = .FALSE.
      IF ( P1x1 >= P2x1-xtol .AND. P1x2 <= P2x2+xtol .AND. &
           P1y1 >= P2y1-ytol .AND. P1y2 <= P2y2+ytol .AND. &
           P1z1 >= P2z1-ztol .AND. P1z2 <= P2z2+ztol ) Is_Within_Bounds = .TRUE.
      RETURN
    END FUNCTION Is_Within_Bounds

    SUBROUTINE COLLECT_NODE_INFO
      IMPLICIT NONE
      !
      ! Now exits, doors, corrs and strs are already read in
      IF (N_NODES > 0 .AND. MYID==MAX(0,EVAC_PROCESS)) THEN
         n_tmp = 0
         DO N = 1, NMESHES
            IF (EVACUATION_ONLY(N).AND.EVACUATION_SKIP(N)) THEN
               n_tmp = n_tmp + 1
               EVAC_Node_List(n_tmp)%Node_Index = n_tmp
               EVAC_Node_List(n_tmp)%Node_Type  = 'Floor'
               EVAC_Node_List(n_tmp)%ID         = MESH_NAME(n)
               EVAC_Node_List(n_tmp)%GRID_NAME  = MESH_NAME(n)
               EVAC_Node_List(n_tmp)%IMESH      = n
            END IF
         END DO
         DO N = 1, N_ENTRYS
            n_tmp = n_tmp + 1
            EVAC_ENTRYS(N)%INODE             = n_tmp
            EVAC_Node_List(n_tmp)%Node_Index = n
            EVAC_Node_List(n_tmp)%Node_Type  = 'Entry'
            EVAC_Node_List(n_tmp)%ID         = EVAC_ENTRYS(n)%ID
            EVAC_Node_List(n_tmp)%IMESH      = EVAC_ENTRYS(n)%IMESH
         END DO
         DO N = 1, N_DOORS
            n_tmp = n_tmp + 1
            EVAC_DOORS(n)%INODE              = n_tmp
            EVAC_Node_List(n_tmp)%Node_Index = n
            EVAC_Node_List(n_tmp)%Node_Type  = 'Door'
            EVAC_Node_List(n_tmp)%ID         = EVAC_DOORS(n)%ID
            EVAC_Node_List(n_tmp)%IMESH      = EVAC_DOORS(n)%IMESH
         END DO
         DO N = 1, N_EXITS
            n_tmp = n_tmp + 1
            EVAC_EXITS(N)%INODE              = n_tmp
            EVAC_Node_List(n_tmp)%Node_Index = n
            EVAC_Node_List(n_tmp)%Node_Type  = 'Exit'
            EVAC_Node_List(n_tmp)%ID         = EVAC_EXITS(n)%ID
            EVAC_Node_List(n_tmp)%IMESH      = EVAC_EXITS(n)%IMESH
         END DO
         DO N = 1, N_CORRS
            n_tmp = n_tmp + 1
            EVAC_CORRS(N)%INODE              = n_tmp
            EVAC_Node_List(n_tmp)%Node_Index = n
            EVAC_Node_List(n_tmp)%Node_Type  = 'Corr'
            EVAC_Node_List(n_tmp)%ID         = EVAC_CORRS(n)%ID
            CorrDoorLoop: DO I = 1, N_DOORS
               IF (TRIM(EVAC_DOORS(I)%TO_NODE)==TRIM(EVAC_CORRS(N)%ID)) THEN
                  EVAC_CORRS(N)%IDOOR_FROM = I
                  IF (EVAC_CORRS(N)%ELEVATOR) EVAC_DOORS(I)%ICORR_TO = N
                  EXIT CorrDoorLoop
               END IF
            END DO CorrDoorLoop
         END DO
         DO N = 1, N_STRS
            n_tmp = n_tmp + 1
            EVAC_STRS(n)%INODE               = n_tmp
            EVAC_Node_List(n_tmp)%Node_Index = n
            EVAC_Node_List(n_tmp)%Node_Type  = 'Stairs'
            EVAC_Node_List(n_tmp)%ID         = EVAC_STRS(n)%ID
         END DO

         ! Check that door/corr/entry/exit have unique names
         DO N = 1, n_nodes - 1
            DO I = n + 1, n_nodes
               IF (TRIM(EVAC_Node_List(n)%ID) == TRIM(EVAC_Node_List(i)%ID)) THEN
                  WRITE(MESSAGE,'(8A)') 'ERROR: ', TRIM(EVAC_Node_List(n)%Node_Type), ': ', &
                       TRIM(EVAC_Node_List(n)%ID), ' has same ID as ', &
                       TRIM(EVAC_Node_List(i)%Node_Type), ': ', TRIM(EVAC_Node_List(i)%ID)
                  CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
               END IF
            END DO
         END DO

         DO N = 1, N_ENTRYS
            DO I = 1, EVAC_ENTRYS(n)%N_VENT_FFIELDS
               IF (EVAC_ENTRYS(n)%I_DOOR_NODES(i) < 0) THEN
                  EVAC_ENTRYS(n)%I_DOOR_NODES(i) = EVAC_EXITS(ABS(EVAC_ENTRYS(n)%I_DOOR_NODES(i)))%INODE
               ELSE IF (EVAC_ENTRYS(n)%I_DOOR_NODES(i) > 0) THEN
                  EVAC_ENTRYS(n)%I_DOOR_NODES(i) = EVAC_DOORS(ABS(EVAC_ENTRYS(n)%I_DOOR_NODES(i)))%INODE
               END IF
            END DO
         END DO

      END IF

    END SUBROUTINE COLLECT_NODE_INFO

    SUBROUTINE READ_ENTRIES
      IMPLICIT NONE
      !
      ! Read the ENTR lines
      !
      ! Local variables
      INTEGER :: nm, i1, i2, j1, j2, NR
      LOGICAL :: L_TMP
      TYPE (EVAC_ENTR_TYPE), POINTER :: PNX=>NULL()
      TYPE (EVAC_PERS_TYPE), POINTER :: PCP=>NULL()
      TYPE (EVAC_STRS_TYPE), POINTER :: STRP=>NULL()
      TYPE (MESH_TYPE), POINTER :: M=>NULL()
      !
      CROWBAR_DUMP = .FALSE.
      READ_ENTR_LOOP: DO N = 1, N_ENTRYS
         !
         ID            = 'null'
         RGB           = -1
         COLOR         = 'null'
         AVATAR_RGB    = -1
         AVATAR_COLOR  = 'null'
         XB            = 0.0_EB
         IOR           = 0
         FLOW_FIELD_ID = 'null'
         MESH_ID       = 'null'
         EVAC_MESH     = 'null'
         TO_NODE       = 'null'
         PERS_ID       = 'null'
         QUANTITY      = 'null'
         AVATAR_TYPE   = 'null'
         PROP_ID       = 'null'
         MAX_FLOW      = 0.0_EB
         WIDTH         = 0.0_EB
         HEIGHT        = 2.0_EB
         AFTER_REACTION_TIME = .FALSE.
         SHOW          = .TRUE.
         TIME_START          = -HUGE(TIME_START)
         TIME_STOP           =  HUGE(TIME_STOP)
         MAX_HUMANS    = -1
         MAX_HUMANS_RAMP         = 'null'
         KNOWN_DOOR_NAMES         = 'null'
         KNOWN_DOOR_PROBS         = 1.0_EB
         AGENT_TYPE               = 2  ! Default is "known door" agent
         CROWBAR_INPUT_FILE = 'null'
         !
         !
         CALL CHECKREAD('ENTR',LU_INPUT,IOS)  ; IF (STOP_STATUS==SETUP_STOP) RETURN
         IF (IOS == 1) THEN
            EXIT READ_ENTR_LOOP
         END IF
         READ(LU_INPUT,ENTR,END=28,IOSTAT=IOS)
         !
         ! Old input used QUANTITY, next lines are needed for that
         IF (MYID==MAX(0,EVAC_PROCESS) .AND. QUANTITY /= 'null') WRITE (LU_ERR,'(A,A)') &
              ' WARNING: keyword QUANTITY is replaced by AVATAR_COLOR at ENTR line ',TRIM(ID)
         IF (QUANTITY == 'BLACK')   AVATAR_COLOR = 'BLACK'
         IF (QUANTITY == 'YELLOW')  AVATAR_COLOR = 'YELLOW'
         IF (QUANTITY == 'BLUE')    AVATAR_COLOR = 'BLUE'
         IF (QUANTITY == 'RED')     AVATAR_COLOR = 'RED'
         IF (QUANTITY == 'GREEN')   AVATAR_COLOR = 'GREEN'
         IF (QUANTITY == 'MAGENTA') AVATAR_COLOR = 'MAGENTA'
         IF (QUANTITY == 'CYAN')    AVATAR_COLOR = 'CYAN'
         !
         IF (MYID==MAX(0,EVAC_PROCESS)) THEN
            IF(TRIM(CROWBAR_INPUT_FILE)/='null' .AND. &
                 ((MAX_FLOW > TWO_EPSILON_EB) .OR. TRIM(MAX_HUMANS_RAMP)/='null')) THEN
               WRITE(MESSAGE,'(A,A,A)') 'ERROR: ENTR line ',TRIM(ID),&
                    ' problem with CROWBAR_INPUT_FILE and MAX_FLOW or MAX_HUMANS_RAMP'
               CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
            END IF
         END IF
         !
         ! Colors, integer RGB(3), e.g., (23,255,0)
         IF (ANY(RGB < 0) .AND. COLOR=='null') COLOR = 'SKY BLUE'
         IF (COLOR /= 'null') CALL COLOR2RGB(RGB,COLOR)
         IF (ANY(AVATAR_RGB < 0) .AND. AVATAR_COLOR=='null') AVATAR_COLOR = 'ROYAL BLUE 4'
         IF (AVATAR_COLOR /= 'null') CALL COLOR2RGB(AVATAR_RGB,AVATAR_COLOR)
         IF (COLOR_METHOD == 0 .AND. (MAX_FLOW >= TWO_EPSILON_EB .OR. &
              Trim(MAX_HUMANS_RAMP)/='null') .OR. TRIM(CROWBAR_INPUT_FILE)/='null') THEN
            i_avatar_color = i_avatar_color + 1
            EVAC_AVATAR_RGB(1:3,i_avatar_color) = AVATAR_RGB
         END IF

         IF (TRIM(CROWBAR_INPUT_FILE)/='null') CROWBAR_DUMP=.TRUE.

         TMP_AVATAR_TYPE_INDEX(N) = 0       ! Defaults for entries that do not generate new agents
         TMP_AVATAR_TYPE_NAME(N)  = 'null'
         TMP_AVATAR_TYPE_PROP(N)  = TRIM(PROP_ID)
         IF (MAX_FLOW >= TWO_EPSILON_EB .OR. TRIM(MAX_HUMANS_RAMP)/='null' &
              .OR. TRIM(CROWBAR_INPUT_FILE)/='null') THEN
            IF (TRIM(AVATAR_TYPE) == 'null' .OR. TRIM(AVATAR_TYPE) == 'Human') THEN
               TMP_AVATAR_TYPE_INDEX(N) = 1
               TMP_AVATAR_TYPE_NAME(N) = TRIM('Human')
            ELSE
               TMP_AVATAR_TYPE_NAME(N)  = TRIM(AVATAR_TYPE)
               TMP_AVATAR_TYPE_INDEX(N) = 0
               Avatar_Type_Loop1: DO I = 1, N-1
                  IF (TRIM(TMP_AVATAR_TYPE_NAME(N)) == TRIM(TMP_AVATAR_TYPE_NAME(I))) THEN
                     TMP_AVATAR_TYPE_INDEX(N) = TMP_AVATAR_TYPE_INDEX(I)
                     EXIT Avatar_Type_Loop1
                  END IF
               END DO Avatar_Type_Loop1
               IF (TMP_AVATAR_TYPE_INDEX(N) == 0) THEN  ! New type found
                  N_AVATAR_TYPE = N_AVATAR_TYPE + 1
                  TMP_AVATAR_TYPE_INDEX(N) = N_AVATAR_TYPE
               END IF
            END IF
         END IF

         IF (MYID /= MAX(0,EVAC_PROCESS)) CYCLE READ_ENTR_LOOP

         ! Restart is not working for evacuation part, but fire+evacuation calculation could
         ! be restarted to produce a good CHID_evac.fed file for later evacuation simulations.
         ! So, the restarted fire+evacuation calculation does not have any agents in it.
         IF (RESTART) THEN
            MAX_HUMANS = 0
            MAX_FLOW   = 0.0_EB
            MAX_HUMANS_RAMP = 'null'
            CROWBAR_INPUT_FILE = 'null'
         END IF

         IF (MAX_HUMANS < 0) MAX_HUMANS = HUGE(MAX_HUMANS)
         IF (MAX_FLOW <= 0.0_EB .AND. TRIM(CROWBAR_INPUT_FILE)=='null') MAX_HUMANS = 0
         IF (Trim(MAX_HUMANS_RAMP)/='null') THEN
            CALL GET_RAMP_INDEX(Trim(MAX_HUMANS_RAMP),'TIME',NR)
            MAX_HUMANS = -NR
         ENDIF

         PNX=>EVAC_ENTRYS(N)

         PNX%AVATAR_TYPE_NAME  = TRIM(TMP_AVATAR_TYPE_NAME(N))
         PNX%AVATAR_TYPE_INDEX = TMP_AVATAR_TYPE_INDEX(N)
         PNX%RGB = RGB
         PNX%AVATAR_RGB =AVATAR_RGB
         IF (TRIM(ID) /= 'null') THEN
            DO I = 1, N-1
               IF (TRIM(ID) == TRIM(EVAC_ENTRYS(I)%ID)) THEN
                  WRITE(MESSAGE,'(A,I4,A,I4,A,A)') 'ERROR: ENTR lines',I,' and',N,', ID strings are not unique: ',TRIM(ID)
                  CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN

               END IF
            END DO
         END IF
         IF (COLOR_METHOD == 0 .AND. (MAX_FLOW >= TWO_EPSILON_EB .OR. &
              TRIM(MAX_HUMANS_RAMP)/='null' .OR. TRIM(CROWBAR_INPUT_FILE)/='null')) &
              PNX%Avatar_Color_Index = i_avatar_color

         IF (EVAC_MESH /= 'null') THEN
            MESH_ID = EVAC_MESH
            IF (MYID==MAX(0,EVAC_PROCESS)) WRITE (LU_ERR,'(A,A)') &
                 ' WARNING: keyword EVAC_MESH is replaced by MESH_ID at ENTR line ', TRIM(ID)
         END IF

         IF (TRIM(KNOWN_DOOR_NAMES(51)) /= 'null') THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: ENTR line ',TRIM(ID), ' problem with KNOWN_DOOR_NAMES'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END IF
         IF (TRIM(KNOWN_DOOR_NAMES(1)) == 'null') THEN
            i = 0 ! no doors given
         ELSE
            i = 50 ! known door names given
            DO WHILE ( TRIM(KNOWN_DOOR_NAMES(i)) == 'null' .AND. i > 0)
               i = i-1
            END DO
         END IF
         PNX%N_VENT_FFIELDS = i
         ALLOCATE(PNX%I_DOOR_NODES(0:i),STAT=IZERO)
         CALL ChkMemErr('Read_Entries','PNX%I_DOOR_NODES',IZERO)
         ALLOCATE(PNX%I_VENT_FFIELDS(0:i),STAT=IZERO)
         CALL ChkMemErr('Read_Entries','PNX%I_VENT_FFIELDS',IZERO)
         ALLOCATE(PNX%P_VENT_FFIELDS(0:i),STAT=IZERO)
         CALL ChkMemErr('Read_Entries','PNX%P_VENT_FFIELDS',IZERO)
         !

         PNX%TO_NODE    = TO_NODE
         PNX%T_first    = T_BEGIN
         PNX%T_last     = T_BEGIN
         PNX%ICOUNT     = 0
         PNX%Flow       = MAX_FLOW
         PNX%T_Start    = TIME_START
         PNX%T_Stop     = TIME_STOP
         PNX%Max_Humans = MAX_HUMANS
         PNX%Max_Humans_Ramp  = MAX_HUMANS_RAMP
         PNX%HEIGHT = HEIGHT
         PNX%SHOW   = SHOW
         PNX%I_AGENT_TYPE = AGENT_TYPE
         PNX%CB_N_Agents = 0 ! Agents in the XB area
         PNX%CB_N_CrowbarAgents = 0 ! Crowbar read in agents in the XB area

         ! Check that the entry is properly specified

         DO I=1,5,2
            IF (XB(I) > XB(I+1)) THEN
               DUMMY   = XB(I)
               XB(I)   = XB(I+1)
               XB(I+1) = DUMMY
            END IF
         END DO

         L_TMP=.FALSE.
         DO i = 1, NMESHES
            IF (.NOT. EVACUATION_ONLY(i)) CYCLE
            IF (TRIM(FLOW_FIELD_ID)==TRIM(MESH_NAME(i))) THEN
               L_TMP=.TRUE.
               EXIT
            END IF
         END DO
         IF (.NOT.(TRIM(FLOW_FIELD_ID)=='null' .OR. L_TMP)) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: ENTR ',TRIM(ID),' problem with FLOW_FIELD_ID'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END IF
         !
         ! Check which evacuation floor
         ii = 0
         n_tmp = 0
         PNX_MeshLoop: DO i = 1, NMESHES
            IF (EVACUATION_ONLY(I) .AND. EVACUATION_SKIP(I)) THEN
               n_tmp = n_tmp + 1
               IF (Is_Within_Bounds(XB(1),XB(2),XB(3),XB(4),XB(5),XB(6),&
                    MESHES(i)%XS,MESHES(i)%XF,MESHES(i)%YS,MESHES(i)%YF,MESHES(i)%ZS,MESHES(i)%ZF, 0._EB, 0._EB, 0._EB)) THEN
                  IF (TRIM(MESH_ID) == 'null' .OR. TRIM(MESH_ID) == TRIM(MESH_NAME(i))) THEN
                     ii = ii + 1
                     PNX%IMESH = i
                     PNX%TO_INODE = n_tmp
                     PNX%TO_NODE  = MESH_NAME(i)
                  END IF
               END IF
            END IF
         END DO PNX_MeshLoop
         IF (PNX%IMESH == 0) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: ENTR ',TRIM(ID), ' problem with IMESH, no mesh found'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END IF
         IF (ii > 1) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: ENTR ',TRIM(ID), ' not an unique mesh found '
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END IF

         nm = PNX%IMESH
         M => MESHES(PNX%IMESH)

         IF (ABS(IOR)<3 .AND. XB(1)/=XB(2) .AND. XB(3)/=XB(4)) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: ENTR ',TRIM(ID),' must be a plane'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         ENDIF
         IF (ABS(IOR)==3 .AND. ABS(XB(1)-XB(2))*ABS(XB(3)-XB(4)) < TWO_EPSILON_EB) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: ENTR ',TRIM(ID),' must have a positive area'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         ENDIF

         ! User input
         PNX%X1 = XB(1)
         PNX%X2 = XB(2)
         PNX%Y1 = XB(3)
         PNX%Y2 = XB(4)
         PNX%Z1 = XB(5)
         PNX%Z2 = XB(6)

         ! Move user input to mesh cell boundaries
         XB(1) = MAX(XB(1),M%XS)
         XB(2) = MIN(XB(2),M%XF)
         XB(3) = MAX(XB(3),M%YS)
         XB(4) = MIN(XB(4),M%YF)
         XB(5) = MAX(XB(5),M%ZS)
         XB(6) = MIN(XB(6),M%ZF)

         IF (ABS(IOR)<3) THEN

            I1 = NINT( GINV(XB(1)-M%XS,1,nm)*M%RDXI )
            I2 = NINT( GINV(XB(2)-M%XS,1,nm)*M%RDXI )
            J1 = NINT( GINV(XB(3)-M%YS,2,nm)*M%RDETA)
            J2 = NINT( GINV(XB(4)-M%YS,2,nm)*M%RDETA)

            XB(1) = M%X(I1)
            XB(2) = M%X(I2)
            XB(3) = M%Y(J1)
            XB(4) = M%Y(J2)

            IF ( ABS(XB(1)-PNX%X1)>1.E-4_EB .OR. ABS(XB(2)-PNX%X2)>1.E-4_EB .OR. &
                 ABS(XB(3)-PNX%Y1)>1.E-4_EB .OR. ABS(XB(4)-PNX%Y2)>1.E-4_EB ) THEN
               WRITE(LU_ERR,fmt='(a,a,a,a)') ' WARNING: Entr ',TRIM(ID),' XB adjusted to mesh ',TRIM(MESH_NAME(nm))
               WRITE(LU_ERR,fmt='(a,6f12.4)') ' Old XB:', PNX%X1,PNX%X2,PNX%Y1,PNX%Y2,PNX%Z1,PNX%Z2
               WRITE(LU_ERR,fmt='(a,6f12.4)') ' New XB:', XB(1:6)
            END IF
            II = (I1+I2)/2
            JJ = (J1+J2)/2
            SELECT CASE (IOR)
            CASE (+1)
               II = II + 1
            CASE (+2)
               JJ = JJ + 1
            END SELECT
            IF (M%SOLID(M%CELL_INDEX(II,JJ,1)) .AND. ABS(IOR)<3) THEN
               WRITE(LU_ERR,fmt='(a,a,a)') ' WARNING: Entr ',TRIM(ID),' problem with XB, mid point facing solid'
            END IF

            ! Coordinates are lined up with the mesh.
            PNX%X1 = XB(1)
            PNX%X2 = XB(2)
            PNX%Y1 = XB(3)
            PNX%Y2 = XB(4)
            PNX%Z1 = XB(5)
            PNX%Z2 = XB(6)
         END IF
         !
         PNX%IOR        = IOR
         PNX%ID         = ID
         PNX%CLASS_NAME = PERS_ID
         PNX%PROP_ID    = TRIM(PROP_ID)
         PNX%CROWBAR_INPUT_FILE = TRIM(CROWBAR_INPUT_FILE)
         PNX%CROWBAR_READ_IN = .FALSE.

         ! Be sure that the crowbar entry is checked at the first initialization time step
!         PNX%CB_TimeLastRead = - EVAC_DT_FLOWFIELD*EVAC_TIME_ITERATIONS + T_BEGIN - 2.0_EB
         PNX%CB_TimeNextRead = - EVAC_DT_FLOWFIELD*EVAC_TIME_ITERATIONS + T_BEGIN - 1.0_EB
         PNX%CB_TimeLastRead = HUGE(PNX%CB_TimeLastRead)
!         PNX%CB_TimeNextRead = HUGE(PNX%CB_TimeNextRead)
         IF(TRIM(CROWBAR_INPUT_FILE)/='null') PNX%CROWBAR_READ_IN = .TRUE.

         ! PNX%Z is used to plot the door on the correct height in Smokeview.
         PNX%Z = 0.5_EB*(PNX%Z1+PNX%Z2) + 0.5_EB*PNX%Height - EVACUATION_Z_OFFSET(PNX%IMESH)

         PNX%IPC = 0
         DO ipc= 1, npc_pers
            pcp => evac_person_classes(ipc)
            IF ( pcp%id == PERS_ID ) PNX%IPC = IPC
         END DO

         SELECT CASE (IOR)
         CASE (-1,+1)
            IF (WIDTH <= 0.0_EB) THEN
               PNX%Width = XB(4) - XB(3)
            ELSE
               PNX%Width = WIDTH
            END IF
            PNX%ORIENTATION(1)=-REAL(SIGN(1,IOR),EB)
         CASE (-2,+2)
            IF (WIDTH <= 0.0_EB) THEN
               PNX%Width = XB(2) - XB(1)
            ELSE
               PNX%Width = WIDTH
            END IF
            PNX%ORIENTATION(2)=-REAL(SIGN(1,IOR),EB)
         CASE (3)
            IF ( (XB(4)-XB(3)) <= TWO_EPSILON_EB .OR. (XB(2)-XB(1)) <= TWO_EPSILON_EB) THEN
               WRITE(MESSAGE,'(A,A,A)') 'ERROR: ENTR',TRIM(ID),' IOR=3 but not 3-dim object'
               CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
            END IF
            PNX%ORIENTATION(3)=-REAL(SIGN(1,IOR),EB)
            PNX%AREA = ABS(XB(4)-XB(3))*ABS(XB(2)-XB(1))  ! XB area used for Crowbar project
         CASE (0)
            IF ( (XB(4)-XB(3)) <= 0.0_EB .OR. (XB(2)-XB(1)) <= 0.0_EB) THEN
               WRITE(MESSAGE,'(A,A,A)') 'ERROR: ENTR',TRIM(ID),' no IOR but not 3-dim object'
               CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
            END IF
         CASE Default
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: ENTR',TRIM(ID),' problem with IOR'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END SELECT

         ! Check if entry leads to Stairs
         PNX%STR_INDX = 0
         PNX%STR_SUB_INDX = 0
         CheckEntrStrLoop: DO i = 1, N_STRS
            IF (TRIM(EVAC_STRS(i)%MESH_ID)==PNX%TO_NODE) THEN
               STRP=>EVAC_STRS(i)
               PNX%STR_INDX = i
               DO j = 1,STRP%N_NODES
                  IF ( PNX%Z1 >= STRP%XB_NODE(j,5) .AND. PNX%Z2 <= STRP%XB_NODE(j,6) .AND. &
                       PNX%Y1 >= STRP%XB_NODE(j,3) .AND. PNX%Y2 <= STRP%XB_NODE(j,4) .AND. &
                       PNX%X1 >= STRP%XB_NODE(j,1) .AND. PNX%X2 <= STRP%XB_NODE(j,2) ) THEN
                     PNX%STR_SUB_INDX = j
                     EXIT CheckEntrStrLoop
                  END IF
               END DO
            END IF
         END DO CheckEntrStrLoop

         ! Use the main_evac_grid flow field if none is given
         IF (TRIM(FLOW_FIELD_ID) == 'null') THEN
            PNX%GRID_NAME  = TRIM(PNX%TO_NODE)
         ELSE
            PNX%GRID_NAME  = TRIM(FLOW_FIELD_ID)
         END IF

         DO i = 1, PNX%N_VENT_FFIELDS
            ! P = 0 or 1 for entrys.
            IF ( .NOT.( ABS(KNOWN_DOOR_PROBS(i)-1.0_EB) < 0.0001_EB .OR. ABS(KNOWN_DOOR_PROBS(i)) < 0.0001_EB ) )  THEN
               WRITE(MESSAGE,'(A,A,A,f12.6,A)') 'ERROR: ENTR line ',TRIM(PNX%ID), &
                    ' problem with probability, ', KNOWN_DOOR_PROBS(i),' it should be zero or one.'
               CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
            END IF
            PNX%P_VENT_FFIELDS(i) = MAX(0.0_EB,KNOWN_DOOR_PROBS(i))
            PNX%I_VENT_FFIELDS(i) = 0
            PNX%I_DOOR_NODES(i) = 0
            DO j = 1, N_EXITS
               IF ( TRIM(EVAC_EXITS(j)%ID) == TRIM(KNOWN_DOOR_NAMES(i)) ) THEN
                  PNX%I_VENT_FFIELDS(i) = EVAC_EXITS(j)%I_VENT_FFIELD
                  PNX%I_DOOR_NODES(i)   = -j
               END IF
            END DO
            DO j = 1, N_DOORS
               IF ( TRIM(EVAC_DOORS(j)%ID) == TRIM(KNOWN_DOOR_NAMES(i)) ) THEN
                  PNX%I_VENT_FFIELDS(i) = EVAC_DOORS(j)%I_VENT_FFIELD
                  PNX%I_DOOR_NODES(i)   = +j
               END IF
            END DO
            IF ( PNX%I_VENT_FFIELDS(i)*PNX%I_DOOR_NODES(i) == 0 ) THEN
               WRITE(MESSAGE,'(A,A,A,A,A)') 'ERROR: ENTR line ',TRIM(PNX%ID), &
                    ' problem with door/exit names, ', TRIM(KNOWN_DOOR_NAMES(i)),' not found'
               CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
            END IF
         END DO
         !
         ! No known doors given, use the flow_field_id value
         !
         PNX%P_VENT_FFIELDS(0) = 1.0_EB
         PNX%I_VENT_FFIELDS(0) = 0
         PNX%I_DOOR_NODES(0) = 0
         PNX_Mesh2Loop: DO i = 1, NMESHES
            IF ( EVACUATION_ONLY(I) .AND. TRIM(PNX%GRID_NAME) == TRIM(MESH_NAME(i)) ) THEN
               PNX%I_VENT_FFIELDS(0) = i
               EXIT PNX_Mesh2Loop
            END IF
         END DO PNX_Mesh2Loop
         IF ( PNX%I_VENT_FFIELDS(0) == 0 ) THEN
            WRITE(MESSAGE,'(A,A,A,A,A)') 'ERROR: ENTR line ',TRIM(PNX%ID),&
                 ' problem with flow field name, ', TRIM(PNX%GRID_NAME),' not found'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END IF
         !
      END DO READ_ENTR_LOOP

28    REWIND(LU_INPUT)

    END SUBROUTINE READ_ENTRIES

    SUBROUTINE READ_EVAC_LINES
      IMPLICIT NONE
      !
      ! Read the EVAC lines
      !
      ! Local variables
      LOGICAL :: L_TMP
      TYPE (EVACUATION_TYPE), POINTER :: HPT=>NULL()
      TYPE (EVAC_PERS_TYPE),  POINTER :: PCP=>NULL()

      READ_EVAC_LOOP: DO N = 1, NPC_EVAC
         !
         ID                       = 'null'
         RGB                      = -1
         COLOR                    = 'null'
         AVATAR_RGB               = -1
         AVATAR_COLOR             = 'null'
         QUANTITY                 = 'null'
         FLOW_FIELD_ID            = 'null'
         MESH_ID                  = 'null'
         EVAC_MESH                = 'null'
         PERS_ID                  = 'null'
         AVATAR_TYPE              = 'null'
         SAMPLING_FACTOR          = 1
         NUMBER_INITIAL_PERSONS   = 0
         XB                       = 0.0_EB
         ANGLE                    = -1000.0_EB
         EVACFILE                 = .FALSE.
         AFTER_REACTION_TIME      = .FALSE.
         SHOW                     = .TRUE.
         TIME_START               = -99.0_EB
         TIME_STOP                = -99.0_EB
         TIME_FALL_DOWN           = HUGE(TIME_FALL_DOWN)
         TIME_START_FED           = T_BEGIN
         GN_MIN                   = 1
         GN_MAX                   = 1
         PRE_EVAC_DIST = -1  ! If Tpre given on EVAC namelist, override PERS
         DET_EVAC_DIST = -1  ! If Tdet given on EVAC namelist, override PERS
         PRE_PARA      = 0.0_EB
         DET_PARA      = 0.0_EB
         PRE_PARA2     = 0.0_EB
         DET_PARA2     = 0.0_EB
         PRE_LOW       = 0.0_EB
         DET_LOW       = T_BEGIN
         PRE_HIGH      = HUGE(PRE_HIGH)
         DET_HIGH      = HUGE(PRE_HIGH)
         PRE_MEAN      = 10.0_EB
         DET_MEAN      = T_BEGIN
         GUARD_MEN_IN  = 0
         TARGET_X      = 0.0_EB
         TARGET_Y      = 0.0_EB
         DELTA_X       = -1.0_EB ! By default no ordered input
         DELTA_Y       = -1.0_EB ! By default no ordered input

         KNOWN_DOOR_NAMES         = 'null'
         KNOWN_DOOR_PROBS         = 1.0_EB
         AGENT_TYPE               = 2  ! Default is "known door" agent
         PROP_ID                  = 'null'
         CROWBAR_INPUT_FILE = 'null'
         !
         CALL CHECKREAD('EVAC',LU_INPUT,IOS)  ; IF (STOP_STATUS==SETUP_STOP) RETURN
         IF (IOS == 1) THEN
            EXIT READ_EVAC_LOOP
         END IF
         READ(LU_INPUT,EVAC,END=25,IOSTAT=IOS)
         !
         ! Old input used QUANTITY, next lines are needed for that
         IF (MYID==MAX(0,EVAC_PROCESS) .AND. QUANTITY /= 'null') WRITE (LU_ERR,'(A,A)') &
              ' WARNING: keyword QUANTITY is replaced by AVATAR_COLOR at EVAC line ',TRIM(ID)
         IF (MYID==MAX(0,EVAC_PROCESS) .AND. TRIM(CROWBAR_INPUT_FILE)/='null' .AND. NUMBER_INITIAL_PERSONS < 1) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: EVAC line ',TRIM(ID),&
                 ' problem with CROWBAR_INPUT_FILE and NUMBER_INITIAL_PERSONS'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END IF
         IF (QUANTITY == 'BLACK')   AVATAR_COLOR = 'BLACK'
         IF (QUANTITY == 'YELLOW')  AVATAR_COLOR = 'YELLOW'
         IF (QUANTITY == 'BLUE')    AVATAR_COLOR = 'BLUE'
         IF (QUANTITY == 'RED')     AVATAR_COLOR = 'RED'
         IF (QUANTITY == 'GREEN')   AVATAR_COLOR = 'GREEN'
         IF (QUANTITY == 'MAGENTA') AVATAR_COLOR = 'MAGENTA'
         IF (QUANTITY == 'CYAN')    AVATAR_COLOR = 'CYAN'
         !
         ! Colors, integer RGB(3), e.g., (23,255,0)
         IF (ANY(RGB < 0) .AND. COLOR=='null') COLOR = 'BLACK'
         IF (COLOR /= 'null') CALL COLOR2RGB(RGB,COLOR)
         IF (ANY(AVATAR_RGB < 0) .AND. AVATAR_COLOR=='null') AVATAR_COLOR = 'ROYAL BLUE 4'
         IF (AVATAR_COLOR /= 'null') CALL COLOR2RGB(AVATAR_RGB,AVATAR_COLOR)
         IF (COLOR_METHOD == 0 .AND. NUMBER_INITIAL_PERSONS > 0) THEN
            I_AVATAR_COLOR = I_AVATAR_COLOR + 1
            EVAC_AVATAR_RGB(1:3,I_AVATAR_COLOR) = AVATAR_RGB
         END IF

         ! Restart is not working for evacuation part, but fire+evacuation calculation could
         ! be restarted to produce a good CHID_evac.fed file for later evacuation simulations.
         ! So, the restarted fire+evacuation calculation does not have any agents in it.
         IF (RESTART) THEN
            NUMBER_INITIAL_PERSONS = 0
            CROWBAR_INPUT_FILE = 'null'
         END IF

         TMP_AVATAR_TYPE_PROP(N+N_ENTRYS) = TRIM(PROP_ID)
         IF (TRIM(AVATAR_TYPE) == 'null' .OR. TRIM(AVATAR_TYPE) == 'Human') THEN
            TMP_AVATAR_TYPE_INDEX(N+N_ENTRYS) = 1
            TMP_AVATAR_TYPE_NAME(N+N_ENTRYS) = TRIM('Human')
         ELSE
            TMP_AVATAR_TYPE_NAME(N+N_ENTRYS)  = TRIM(AVATAR_TYPE)
            TMP_AVATAR_TYPE_INDEX(N+N_ENTRYS) = 0
            Avatar_Type_Loop: DO I = 1+N_ENTRYS, N-1+N_ENTRYS
               IF (TRIM(TMP_AVATAR_TYPE_NAME(N+N_ENTRYS)) == TRIM(TMP_AVATAR_TYPE_NAME(I))) THEN
                  TMP_AVATAR_TYPE_INDEX(N+N_ENTRYS) = TMP_AVATAR_TYPE_INDEX(I)
                  EXIT Avatar_Type_Loop
               END IF
            END DO Avatar_Type_Loop
            Avatar_Type_Loop2: DO I = 1, N_ENTRYS
               IF (TRIM(TMP_AVATAR_TYPE_NAME(N+N_ENTRYS)) == TRIM(TMP_AVATAR_TYPE_NAME(I)) .AND. &
                    TMP_AVATAR_TYPE_INDEX(N+N_ENTRYS) == 0) THEN
                  TMP_AVATAR_TYPE_INDEX(N+N_ENTRYS) = TMP_AVATAR_TYPE_INDEX(I)
                  EXIT Avatar_Type_Loop2
               END IF
            END DO Avatar_Type_Loop2
            IF (TMP_AVATAR_TYPE_INDEX(N+N_ENTRYS) == 0) THEN  ! New type found
               N_AVATAR_TYPE = N_AVATAR_TYPE + 1
               TMP_AVATAR_TYPE_INDEX(N+N_ENTRYS) = N_AVATAR_TYPE
            END IF
         END IF

         IF (MYID /= MAX(0,EVAC_PROCESS)) CYCLE READ_EVAC_LOOP

         HPT=>EVAC_EVACS(N)

         HPT%AVATAR_TYPE_NAME  = TRIM(TMP_AVATAR_TYPE_NAME(N+N_ENTRYS))
         HPT%AVATAR_TYPE_INDEX = TMP_AVATAR_TYPE_INDEX(N+N_ENTRYS)
         HPT%RGB = RGB
         HPT%AVATAR_RGB = AVATAR_RGB
         IF (COLOR_METHOD == 0 .AND. NUMBER_INITIAL_PERSONS > 0) HPT%Avatar_Color_Index = i_avatar_color
         IF (EVAC_MESH /= 'null') THEN
            MESH_ID = EVAC_MESH
            IF (MYID==MAX(0,EVAC_PROCESS)) WRITE (LU_ERR,'(A,A)') &
                 ' WARNING: keyword EVAC_MESH is replaced by MESH_ID at EVAC line ', TRIM(ID)
         END IF

         IF (TRIM(ID) /= 'null') THEN
            DO I = 1, N-1
               IF (TRIM(ID) == TRIM(EVAC_EVACS(I)%ID)) THEN
                  WRITE(MESSAGE,'(A,I4,A,I4,A,A)') 'ERROR: EVAC lines',I,' and',N,', ID strings are not unique: ',TRIM(ID)
                  CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN

               END IF
            END DO
         END IF
         IF (TRIM(PERS_ID) == 'null') THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: EVAC line ',TRIM(ID),' no PERS_ID given'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         ELSE
            ii = 1
            DO i = 1,NPC_PERS
               IF (TRIM(EVAC_PERSON_CLASSES(i)%ID) == TRIM(PERS_ID)) CYCLE
               ii = ii + 1
            END DO
            IF (ii > NPC_PERS) THEN
               WRITE(MESSAGE,'(A,A,A)') 'ERROR: EVAC line ',TRIM(ID), ' problem with PERS_ID'
               CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
            END IF
         END IF

         IF (TRIM(KNOWN_DOOR_NAMES(51)) /= 'null') THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: EVAC line ',TRIM(ID), ' problem with KNOWN_DOOR_NAMES'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END IF
         IF (TRIM(KNOWN_DOOR_NAMES(1)) == 'null') THEN
            i = 0 ! no doors given
         ELSE
            i = 50 ! known door names given
            DO WHILE ( TRIM(KNOWN_DOOR_NAMES(i)) == 'null' .AND. i > 1)
               i = i-1
            END DO
         END IF
         HPT%N_VENT_FFIELDS = i
         ALLOCATE(HPT%I_DOOR_NODES(0:i),STAT=IZERO)
         CALL ChkMemErr('Read_Evac','HPT%I_DOOR_NODES',IZERO)
         ALLOCATE(HPT%I_VENT_FFIELDS(0:i),STAT=IZERO)
         CALL ChkMemErr('Read_Evac','HPT%I_VENT_FFIELDS',IZERO)
         ALLOCATE(HPT%P_VENT_FFIELDS(0:i),STAT=IZERO)
         CALL ChkMemErr('Read_Evac','HPT%P_VENT_FFIELDS',IZERO)
         !
         IF (NUMBER_INITIAL_PERSONS > 0 .AND. RESTART) NUMBER_INITIAL_PERSONS =  0
         !
         IF (NUMBER_INITIAL_PERSONS > 0) EVACFILE = .TRUE.
         !
         HPT%CLASS_NAME = PERS_ID
         HPT%T_START    = TIME_START
         HPT%T_START_FED = TIME_START_FED
         HPT%I_AGENT_TYPE = AGENT_TYPE
         HPT%PROP_ID    = TRIM(PROP_ID)
         HPT%TIME_FALL_DOWN = TIME_FALL_DOWN

         HPT%GN_MIN = GN_MIN
         HPT%GN_MAX = GN_MAX

         ! input in degrees, internal units are radians, [0,2pi)
         ! If no angle is given then use random angle
         IF (ANGLE > -999.9_EB) THEN
            ANGLE = Pi*ANGLE/(180.0_EB)
            DO WHILE (ANGLE >= 2.0_EB*Pi)
               ANGLE = ANGLE - 2.0_EB*Pi
            END DO
            DO WHILE (ANGLE < 0.0_EB)
               ANGLE = ANGLE + 2.0_EB*Pi
            END DO
         END IF
         HPT%Angle = ANGLE

         HPT%IPC = 0
         DO ipc= 1, npc_pers
            pcp => evac_person_classes(ipc)
            IF ( TRIM(pcp%id) == TRIM(PERS_ID) ) HPT%IPC = IPC
         END DO
         !
         HPT%SAMPLING = SAMPLING_FACTOR
         !
         IF (NUMBER_INITIAL_PERSONS > 0) THEN
            DO I=1,5,2
               IF (XB(I) > XB(I+1)) THEN
                  DUMMY   = XB(I)
                  XB(I)   = XB(I+1)
                  XB(I+1) = DUMMY
               END IF
            END DO
         END IF
         !
         HPT%X1 = XB(1)
         HPT%X2 = XB(2)
         HPT%Y1 = XB(3)
         HPT%Y2 = XB(4)
         HPT%Z1 = XB(5)
         HPT%Z2 = XB(6)
         HPT%N_INITIAL = NUMBER_INITIAL_PERSONS
         HPT%EVACFILE = EVACFILE
         HPT%IMESH = 0
         HPT%ID = TRIM(ID)
         HPT%SHOW   = SHOW
         HPT%I_PRE_DIST  = PRE_EVAC_DIST
         HPT%Tpre_mean   = PRE_MEAN
         HPT%Tpre_low    = PRE_LOW
         HPT%Tpre_high   = PRE_HIGH
         HPT%Tpre_para   = PRE_PARA
         HPT%Tpre_para2  = PRE_PARA2
         HPT%I_DET_DIST  = DET_EVAC_DIST
         HPT%Tdet_mean   = DET_MEAN
         HPT%Tdet_low    = DET_LOW
         HPT%Tdet_high   = DET_HIGH
         HPT%Tdet_para   = DET_PARA
         HPT%Tdet_para2  = DET_PARA2
         HPT%GUARD_MEN   = GUARD_MEN_IN
         HPT%TARGET_X    = TARGET_X
         HPT%TARGET_Y    = TARGET_Y
         HPT%CROWBAR_INPUT_FILE = TRIM(CROWBAR_INPUT_FILE)
         IF (DELTA_X>0.0_EB .AND. DELTA_Y>0.0_EB) THEN
            HPT%DELTA_X = DELTA_X ; HPT%DELTA_Y = DELTA_Y
         ELSE
            HPT%DELTA_X = -1.0_EB ; HPT%DELTA_Y = -1.0_EB
         END IF

         L_TMP=.FALSE.
         DO i = 1, NMESHES
            IF (.NOT. EVACUATION_ONLY(i)) CYCLE
            IF (TRIM(FLOW_FIELD_ID)==TRIM(MESH_NAME(i))) THEN
               L_TMP=.TRUE.
               EXIT
            END IF
         END DO
         IF (.NOT.(TRIM(FLOW_FIELD_ID)=='null' .OR. L_TMP)) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: EVAC ',TRIM(ID),' problem with FLOW_FIELD_ID'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END IF

         ! Check which evacuation floor
         ii = 0
         HP_MeshLoop: DO i = 1, nmeshes
            IF (EVACUATION_ONLY(I) .AND. EVACUATION_SKIP(I)) THEN
               IF ( Is_Within_Bounds(HPT%X1,HPT%X2,HPT%Y1,HPT%Y2,HPT%Z1,HPT%Z2,&
                    MESHES(i)%XS,MESHES(i)%XF,MESHES(i)%YS,MESHES(i)%YF,MESHES(i)%ZS,MESHES(i)%ZF, 0._EB, 0._EB, 0._EB)) THEN
                  IF (TRIM(MESH_ID) == 'null' .OR. TRIM(MESH_ID) == TRIM(MESH_NAME(i))) THEN
                     ii = ii + 1
                     HPT%IMESH = i
                  END IF
               END IF
            END IF
         END DO HP_MeshLoop
         IF (HPT%IMESH == 0) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: EVAC line ',TRIM(ID),' problem with IMESH, no mesh found'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END IF
         IF (ii > 1) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: EVAC line ',TRIM(ID), ' not an unique mesh found '
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END IF

         ! Use the main_evac_grid flow field if none is given
         IF (TRIM(FLOW_FIELD_ID) == 'null') THEN
            HPT%GRID_NAME  = TRIM(MESH_NAME(HPT%IMESH))
         ELSE
            HPT%GRID_NAME  = TRIM(FLOW_FIELD_ID)
         END IF
         !
         DO i = 1, HPT%N_VENT_FFIELDS
            HPT%P_VENT_FFIELDS(i) = MAX(0.0_EB,KNOWN_DOOR_PROBS(i))
            HPT%I_VENT_FFIELDS(i) = 0
            HPT%I_DOOR_NODES(i) = 0
            DO j = 1, N_EXITS
               IF ( TRIM(EVAC_EXITS(j)%ID) == TRIM(KNOWN_DOOR_NAMES(i)) ) THEN
                  HPT%I_VENT_FFIELDS(i) = EVAC_EXITS(j)%I_VENT_FFIELD
                  HPT%I_DOOR_NODES(i)   = EVAC_EXITS(j)%INODE
               END IF
            END DO
            DO j = 1, N_DOORS
               IF ( TRIM(EVAC_DOORS(j)%ID) == TRIM(KNOWN_DOOR_NAMES(i)) ) THEN
                  HPT%I_VENT_FFIELDS(i) = EVAC_DOORS(j)%I_VENT_FFIELD
                  HPT%I_DOOR_NODES(i)   = EVAC_DOORS(j)%INODE
               END IF
            END DO
            IF ( HPT%I_VENT_FFIELDS(i)*HPT%I_DOOR_NODES(i) == 0 ) THEN
               WRITE(MESSAGE,'(A,A,A,A,A)') 'ERROR: EVAC line ',TRIM(HPT%ID), &
                    ' problem with door/exit names, ', TRIM(KNOWN_DOOR_NAMES(i)),' not found'
               CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
            END IF
         END DO
         !
         ! No known doors given, use the flow_field_id value
         !
         HPT%P_VENT_FFIELDS(0) = 1.0_EB
         HPT%I_VENT_FFIELDS(0) = 0
         HPT%I_DOOR_NODES(0) = 0
         HP_Mesh2Loop: DO i = 1, NMESHES
            IF ( EVACUATION_ONLY(I) .AND. TRIM(HPT%GRID_NAME) == TRIM(MESH_NAME(I)) ) THEN
               HPT%I_VENT_FFIELDS(0) = i
               EXIT HP_Mesh2Loop
            END IF
         END DO HP_Mesh2Loop
         IF ( HPT%I_VENT_FFIELDS(0) == 0 ) THEN
            WRITE(MESSAGE,'(A,A,A,A,A)') 'ERROR: EVAC line ',TRIM(HPT%ID), &
                 ' problem with flow field name, ', TRIM(HPT%GRID_NAME),' not found'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END IF
         !
      END DO READ_EVAC_LOOP
25    REWIND(LU_INPUT)


      IF (CROWBAR_DUMP) THEN
         ! Dump also the read in positions
         N_AVATAR_TYPE = N_AVATAR_TYPE + 1
      END IF

      N_EVAC = MAX(1,N_AVATAR_TYPE)
      ALLOCATE(EVAC_CLASS_NAME(N_EVAC),STAT=IZERO)
      CALL ChkMemErr('READ_EVAC','EVAC_CLASS_NAME',IZERO)
      ALLOCATE(EVAC_CLASS_RGB(3,N_EVAC),STAT=IZERO)
      CALL ChkMemErr('READ_EVAC','EVAC_CLASS_RGB',IZERO)
      ALLOCATE(EVAC_CLASS_PROP(N_EVAC),STAT=IZERO)
      CALL ChkMemErr('READ_EVAC','EVAC_CLASS_PROP',IZERO)
      EVAC_CLASS_NAME(1) = TRIM('Human')
      EVAC_CLASS_PROP = TRIM('null')
      DO N = 1, N_ENTRYS + NPC_EVAC
         IF (TMP_AVATAR_TYPE_INDEX(N) > 0) THEN
            I = TMP_AVATAR_TYPE_INDEX(N)
            EVAC_CLASS_NAME(I) = TRIM(TMP_AVATAR_TYPE_NAME(N))
            IF (TRIM(TMP_AVATAR_TYPE_PROP(N)) /= 'null') EVAC_CLASS_PROP(I) = TRIM(TMP_AVATAR_TYPE_PROP(N))
         END IF
      END DO

      DO N = 1, N_EVAC
         EVAC_CLASS_RGB(1:3,N) = (/ 39, 64,139/)  ! ROYAL BLUE 4
      END DO
      IF (CROWBAR_DUMP) THEN
         ! Dump also the read in positions
         EVAC_CLASS_NAME(N_EVAC) = TRIM('Crowbar')
         EVAC_CLASS_PROP(N_EVAC) = TRIM('Crowbar_props')
         EVAC_CLASS_RGB(1:3,N_EVAC) = (/ 0, 0, 0/)  ! Crowbar is BLACK
      END IF

    END SUBROUTINE READ_EVAC_LINES

    SUBROUTINE READ_EDEV
      IMPLICIT NONE
      !
      ! Read the EDEV lines
      !
      ! Local variables
      TYPE (EVAC_EDEV_TYPE), POINTER :: EDV=>NULL()

      READ_EDEV_LOOP: DO N = 1, N_EDEV
         IF (MYID /= MAX(0,EVAC_PROCESS)) CYCLE READ_EDEV_LOOP
         EDV => EVAC_EDEV(N)
         !
         ID            = 'null'
         EVAC_ID       = 'null'
         PERS_ID       = 'null'
         MESH_ID       = 'null'
         TIME_DELAY    = 0.0_EB
         GLOBAL        = .TRUE.
         PROB          = 0.0_EB
         PRE_EVAC_DIST = -1 ! If pre dist given on EDEV, override evac/pers values
         PRE_MEAN      = 0.0_EB
         PRE_PARA      = 0.0_EB
         PRE_PARA2     = 0.0_EB
         PRE_LOW       = 0.0_EB
         PRE_HIGH      = HUGE(PRE_HIGH)
         INPUT_ID      = 'null'
         !
         CALL CHECKREAD('EDEV',LU_INPUT,IOS)  ; IF (STOP_STATUS==SETUP_STOP) RETURN
         IF (IOS == 1) THEN
            EXIT READ_EDEV_LOOP
         END IF
         READ(LU_INPUT,EDEV,END=30,IOSTAT=IOS)
         !
         EDV%ID         = TRIM(ID)
         EDV%EVAC_ID    = TRIM(EVAC_ID)
         EDV%PERS_ID    = TRIM(PERS_ID)
         EDV%GLOBAL     = GLOBAL
         EDV%TIME_DELAY = TIME_DELAY
         EDV%PROB       = MIN(1.0_EB,MAX(0.0_EB,PROB))
         EDV%I_pre_dist = PRE_EVAC_DIST
         EDV%Tpre_mean  = PRE_MEAN
         EDV%Tpre_para  = PRE_PARA
         EDV%Tpre_para2 = PRE_PARA2
         EDV%Tpre_low   = PRE_LOW
         EDV%Tpre_high  = PRE_HIGH

         ! Check which main evacuation mesh if MESH_ID given
         EDV%IMESH     = 0
         EDV%MESH_ID   = TRIM(MESH_ID)
         IF (TRIM(MESH_ID) /= 'null') THEN
            EDV_MeshLoop: DO I = 1, NMESHES
               IF (EVACUATION_ONLY(I) .AND. EVACUATION_SKIP(I) .AND. &
                    TRIM(MESH_ID) == TRIM(MESH_NAME(i))) THEN
                  EDV%IMESH = I
                  EXIT EDV_MeshLoop
               END IF
            END DO EDV_MeshLoop
         END IF
         IF (TRIM(MESH_ID) /= 'null' .AND. EDV%IMESH == 0) THEN
            IF (TRIM(ID)=='null') THEN
               WRITE(MESSAGE,'(A,I5,A,A,A)')  'ERROR: EDEV ',N,', ',TRIM(MESH_ID),' is not a main evacuation mesh'
            ELSE
               WRITE(MESSAGE,'(A,A,A,A,A)')  'ERROR: EDEV ',TRIM(ID),', ',TRIM(MESH_ID),' is not a main evacuation mesh'
            END IF
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         ENDIF
         IF (TRIM(MESH_ID) == 'null' .AND. .NOT.GLOBAL) THEN
            WRITE(MESSAGE,'(A,I4,A)')  'ERROR: EDEV ',N,' is not global, it must have MESH_ID'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         ENDIF

         ! Check the INPUT_IDs
         EDV%N_INPUTS  = 0
         INPUT_COUNT: DO
            IF (EDV%N_INPUTS==40) THEN
               EXIT INPUT_COUNT
            END IF
            IF (INPUT_ID(EDV%N_INPUTS+1)=='null') THEN
               EXIT INPUT_COUNT
            END IF
            EDV%N_INPUTS  = EDV%N_INPUTS + 1
         END DO INPUT_COUNT
         IF (EDV%N_INPUTS==0) THEN
            IF (TRIM(ID)=='null') THEN
               WRITE(MESSAGE,'(A,I5,A)')  'ERROR: EDEV ',N,' must have at least one input'
            ELSE
               WRITE(MESSAGE,'(A,A,A)')  'ERROR: EDEV ',TRIM(ID),' must have at least one input'
            END IF
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         ENDIF
         ALLOCATE(EDV%INPUT_ID(EDV%N_INPUTS),STAT=IZERO)
         CALL ChkMemErr('READ_EVAC','INPUT_ID',IZERO)
         ALLOCATE(EDV%INPUT_DEVC_INDEX(EDV%N_INPUTS),STAT=IZERO)
         CALL ChkMemErr('READ_EVAC','INPUT_DEVC_INDEX',IZERO)
         EDV%INPUT_ID(1:EDV%N_INPUTS) = INPUT_ID(1:EDV%N_INPUTS)

      END DO READ_EDEV_LOOP
30    REWIND(LU_INPUT)

    END SUBROUTINE READ_EDEV

    SUBROUTINE READ_EVHO
      IMPLICIT NONE
      !
      ! Read the EVHO lines
      !
      ! Local variables
      TYPE (EVAC_HOLE_TYPE),  POINTER :: EHX=>NULL()

      READ_EVHO_LOOP: DO N = 1, N_HOLES
         IF (MYID /= MAX(0,EVAC_PROCESS)) CYCLE READ_EVHO_LOOP
         EHX=>EVAC_HOLES(N)
         !
         ID            = 'null'
         RGB           = -1
         COLOR         = 'null'
         XB            = 0.0_EB
         EVAC_ID       = 'null'
         PERS_ID       = 'null'
         MESH_ID       = 'null'
         EVAC_MESH     = 'null'
         SHOW          = .TRUE.
         TIME_FALL_DOWN = -1.0_EB
         !
         CALL CHECKREAD('EVHO',LU_INPUT,IOS)  ; IF (STOP_STATUS==SETUP_STOP) RETURN
         IF (IOS == 1) THEN
            EXIT READ_EVHO_LOOP
         END IF
         READ(LU_INPUT,EVHO,END=30,IOSTAT=IOS)
         !
         DO I=1,5,2
            IF (XB(I) > XB(I+1)) THEN
               DUMMY   = XB(I)
               XB(I)   = XB(I+1)
               XB(I+1) = DUMMY
            END IF
         END DO
         IF (EVAC_MESH /= 'null') THEN
            MESH_ID = EVAC_MESH
            IF (MYID==MAX(0,EVAC_PROCESS)) WRITE (LU_ERR,'(A,A)') &
                 ' WARNING: keyword EVAC_MESH is replaced by MESH_ID at EVHO line ', TRIM(ID)
         END IF

         !
         EHX%X1 = XB(1)
         EHX%X2 = XB(2)
         EHX%Y1 = XB(3)
         EHX%Y2 = XB(4)
         EHX%Z1 = XB(5)
         EHX%Z2 = XB(6)
         EHX%ID        = ID
         EHX%EVAC_ID   = EVAC_ID
         EHX%PERS_ID   = PERS_ID
         EHX%SHOW      = SHOW
         EHX%TIME_FALL_DOWN = TIME_FALL_DOWN
         !
         ! Check which evacuation floor
         ii = 0
         EHX_MeshLoop: DO i = 1, NMESHES
            IF (EVACUATION_ONLY(I) .AND. EVACUATION_SKIP(I)) THEN
               IF ( (EHX%Z1 >= MESHES(i)%ZS .AND. EHX%Z2 <= MESHES(i)%ZF).AND. &
                    (EHX%Y1 >= MESHES(i)%YS .AND. EHX%Y2 <= MESHES(i)%YF).AND. &
                    (EHX%X1 >= MESHES(i)%XS .AND. EHX%X2 <= MESHES(i)%XF)) THEN
                  IF (TRIM(MESH_ID) == 'null' .OR. TRIM(MESH_ID) == TRIM(MESH_NAME(i))) THEN
                     ii = ii + 1
                     EHX%IMESH = i
                     EHX%GRID_NAME  = MESH_NAME(i)
                  END IF
               END IF
            END IF
         END DO EHX_MeshLoop
         IF (EHX%IMESH == 0) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: EVHO line ',TRIM(EHX%ID), ' problem with IMESH, no mesh found'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.)
         END IF
         IF (ii > 1) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: EVHO line ',TRIM(EHX%ID), ' not an unique mesh found '
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END IF
         IF (TRIM(ID) /= 'null') THEN
            DO I = 1, N-1
               IF (TRIM(ID) == TRIM(EVAC_HOLES(I)%ID)) THEN
                  WRITE(MESSAGE,'(A,I4,A,I4,A,A)') 'ERROR: EVHO lines',I,' and',N,', ID strings are not unique: ',TRIM(ID)
                  CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN

               END IF
            END DO
         END IF

         ! Colors, integer RGB(3), e.g., (23,255,0)
         IF (ANY(RGB < 0) .AND. COLOR=='null') COLOR = 'BLACK'
         IF (COLOR /= 'null') CALL COLOR2RGB(RGB,COLOR)
         EHX%RGB = RGB
      END DO READ_EVHO_LOOP
30    REWIND(LU_INPUT)

    END SUBROUTINE READ_EVHO

    SUBROUTINE READ_EVSS
      IMPLICIT NONE
      !
      ! Read the EVSS lines
      !
      ! Local variables
      TYPE (EVAC_SSTAND_TYPE), POINTER :: ESS=>NULL()
      REAL(EB) :: X, Y, Z
      LOGICAL L_TMP

      READ_EVSS_LOOP: DO N = 1, N_SSTANDS
         IF (MYID /= MAX(0,EVAC_PROCESS)) CYCLE READ_EVSS_LOOP
         ESS => EVAC_SSTANDS(N)
         !
         ID            = 'null'
         RGB           = -1
         COLOR         = 'null'
         XB            = 0.0_EB
         MESH_ID       = 'null'
         EVAC_MESH     = 'null'
         VENT_FFIELD   = 'null'
         IOR           = 0
         HEIGHT        = 0.0_EB
         HEIGHT0       = 0.0_EB
         FAC_V0_UP     = 1.0_EB
         FAC_V0_DOWN   = 1.0_EB
         FAC_V0_HORI   = 1.0_EB
         ESC_SPEED     = 0.0_EB
         UBAR0         = 0.0_EB
         VBAR0         = 0.0_EB
         USE_V0        = .FALSE.
         SHOW          = .TRUE.
         !
         CALL CHECKREAD('EVSS',LU_INPUT,IOS)  ; IF (STOP_STATUS==SETUP_STOP) RETURN
         IF (IOS == 1) THEN
            EXIT READ_EVSS_LOOP
         END IF
         READ(LU_INPUT,EVSS,END=31,IOSTAT=IOS)
         !
         DO I=1,5,2
            IF (XB(I) > XB(I+1)) THEN
               DUMMY   = XB(I)
               XB(I)   = XB(I+1)
               XB(I+1) = DUMMY
            END IF
         END DO
         IF (EVAC_MESH /= 'null') THEN
            MESH_ID = EVAC_MESH
            IF (MYID==MAX(0,EVAC_PROCESS)) WRITE (LU_ERR,'(A,A)') &
                 ' WARNING: keyword EVAC_MESH is replaced by MESH_ID at EVSS line ', TRIM(ID)
         END IF
         !
         ESS%X1 = XB(1)
         ESS%X2 = XB(2)
         ESS%Y1 = XB(3)
         ESS%Y2 = XB(4)
         ESS%Z1 = XB(5)
         ESS%Z2 = XB(6)
         ESS%ID          = ID
         ESS%H           = HEIGHT
         ESS%H0          = HEIGHT0
         ESS%FAC_V0_HORI = FAC_V0_HORI
         IF ( (ESS%H - ESS%H0) < 0.0_EB) THEN
            ESS%FAC_V0_UP   = FAC_V0_DOWN
            ESS%FAC_V0_DOWN = FAC_V0_UP
            ESS%Esc_SpeedDn   = MAX(0.0_EB,+ESC_SPEED)
            ESS%Esc_SpeedUp   = MAX(0.0_EB,-ESC_SPEED)
         ELSE
            ESS%FAC_V0_UP   = FAC_V0_UP
            ESS%FAC_V0_DOWN = FAC_V0_DOWN
            ESS%Esc_SpeedUp   = MAX(0.0_EB,+ESC_SPEED)
            ESS%Esc_SpeedDn   = MAX(0.0_EB,-ESC_SPEED)
         END IF
         ESS%IOR    = IOR
         ESS%UBAR0  = UBAR0
         ESS%VBAR0  = VBAR0
         ESS%Use_v0 = USE_V0
         ESS%SHOW   = SHOW
         ESS%VENT_FFIELD= Trim(VENT_FFIELD)
         !
         ! Check which evacuation floor
         ii = 0
         ESS_MeshLoop: DO i = 1, NMESHES
            IF (EVACUATION_ONLY(I) .AND. EVACUATION_SKIP(I)) THEN
               IF ( (ESS%Z1 >= MESHES(i)%ZS .AND. ESS%Z2 <= MESHES(i)%ZF).AND. &
                    (ESS%Y1 >= MESHES(i)%YS .AND. ESS%Y2 <= MESHES(i)%YF).AND. &
                    (ESS%X1 >= MESHES(i)%XS .AND. ESS%X2 <= MESHES(i)%XF)) THEN
                  IF (TRIM(MESH_ID) == 'null' .OR. TRIM(MESH_ID) == TRIM(MESH_NAME(i))) THEN
                     ii = ii + 1
                     ESS%IMESH = i
                     ESS%GRID_NAME  = MESH_NAME(i)
                  END IF
               END IF
            END IF
         END DO ESS_MeshLoop
         IF (ESS%IMESH == 0) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: EVSS line ',TRIM(ESS%ID), ' problem with IMESH, no mesh found'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END IF
         IF (ii > 1) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: EVSS line ',TRIM(ESS%ID), ' not an unique mesh found '
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END IF
         IF (TRIM(ID) /= 'null') THEN
            DO I = 1, N-1
               IF (TRIM(ID) == TRIM(EVAC_SSTANDS(I)%ID)) THEN
                  WRITE(MESSAGE,'(A,I4,A,I4,A,A)') 'ERROR: EVSS lines',I,' and',N,', ID strings are not unique: ',TRIM(ID)
                  CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN

               END IF
            END DO
         END IF
         L_TMP=.FALSE.
         DO I = 1, NMESHES
            IF (.NOT. EVACUATION_ONLY(I)) CYCLE
            IF (TRIM(VENT_FFIELD)==TRIM(MESH_NAME(I))) THEN
               L_TMP=.TRUE.
               EXIT
            END IF
         END DO
         IF (.NOT.(TRIM(VENT_FFIELD)=='null' .OR. L_TMP)) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: EVSS ',TRIM(ID),' problem with VENT_FFIELD'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END IF
         ! Check which vent field. If VENT_FFIELD is not found, use the main evac grid.
         ESS%I_VENT_FFIELD = 0
         ESS_Mesh2Loop: DO I = 1, NMESHES
            IF ( EVACUATION_ONLY(I) .AND. (TRIM(MESH_NAME(I)) == TRIM(VENT_FFIELD)) ) THEN
               IF ( (ESS%Z1 >= MESHES(I)%ZS .AND. ESS%Z2 <= MESHES(I)%ZF).AND. &
                    (ESS%Y1 >= MESHES(I)%YS .AND. ESS%Y2 <= MESHES(I)%YF).AND. &
                    (ESS%X1 >= MESHES(I)%XS .AND. ESS%X2 <= MESHES(I)%XF)) THEN
                  ESS%I_VENT_FFIELD = I
                  EXIT ESS_Mesh2Loop
               END IF
            END IF
         END DO ESS_Mesh2Loop
         ! If no vent field is given, then use the main evac grid.
         IF (ESS%I_VENT_FFIELD == 0) THEN
            ESS%I_VENT_FFIELD = ESS%IMESH
            ESS%VENT_FFIELD = TRIM(MESH_NAME(ESS%IMESH))
         END IF

         IF (ABS(ESS%H0-ESS%H) < 1.0E-3) THEN
            ESS%H = ESS%H0
            IOR=1
            ESS%IOR    = IOR
         END IF

         SELECT CASE (IOR)
         CASE(-1,+1)
            ESS%S = SQRT((ESS%X2-ESS%X1)**2 + (ESS%H-ESS%H0)**2)
            ESS%COS_X = ABS(ESS%X2-ESS%X1)/ SQRT((ESS%X2-ESS%X1)**2 + (ESS%H-ESS%H0)**2)
            ESS%COS_Y = 1.0_EB
            ESS%SIN_X = ABS(ESS%H-ESS%H0)/ SQRT((ESS%X2-ESS%X1)**2 + (ESS%H-ESS%H0)**2)
            ESS%SIN_Y = 0.0_EB
            ESS%ORIENTATION(1) = IOR*(ESS%H-ESS%H0)/ SQRT((ESS%X2-ESS%X1)**2 + (ESS%H-ESS%H0)**2)
            ESS%ORIENTATION(2) = 0.0_EB
            ESS%ORIENTATION(3) = ESS%COS_X
         CASE(-2,+2)
            ESS%S = SQRT((ESS%Y2-ESS%Y1)**2 + (ESS%H-ESS%H0)**2)
            ESS%COS_X = 1.0_EB
            ESS%COS_Y = ABS(ESS%Y2-ESS%Y1)/ SQRT((ESS%Y2-ESS%Y1)**2 + (ESS%H-ESS%H0)**2)
            ESS%SIN_X = 0.0_EB
            ESS%SIN_Y = ABS(ESS%H-ESS%H0)/ SQRT((ESS%Y2-ESS%Y1)**2 + (ESS%H-ESS%H0)**2)
            ESS%ORIENTATION(1) = 0.0_EB
            ESS%ORIENTATION(2) = 0.5_EB*IOR*(ESS%H-ESS%H0)/ SQRT((ESS%Y2-ESS%Y1)**2 + (ESS%H-ESS%H0)**2)
            ESS%ORIENTATION(3) = ESS%COS_Y
         CASE Default
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: EVSS ',TRIM(ESS%ID),' problem with IOR'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END SELECT

         ! Colors, integer RGB(3), e.g., (23,255,0)
         IF (ANY(RGB < 0) .AND. COLOR=='null') COLOR = 'AZURE 2'  ! RGB = 193 205 205
         IF (COLOR /= 'null') CALL COLOR2RGB(RGB,COLOR)
         ESS%RGB = RGB

         ! Check if doors, exits and entrys are on the inclines.
         ! If they are, change the 'z' so that Smokeview plots correctly.
         DO I = 1, N_EXITS
            IF (EVAC_EXITS(I)%IMESH /= ESS%IMESH) CYCLE
            X = EVAC_EXITS(I)%X
            Y = EVAC_EXITS(I)%Y
            Z = 0.0_EB
            IF ( (ESS%X1 <= X .AND. ESS%X2 >= X) .AND. (ESS%Y1 <= Y .AND. ESS%Y2 >= Y)) THEN
               SELECT CASE (ESS%IOR)
               CASE(-1)
                  Z = ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%X1-X)/ABS(ESS%X1-ESS%X2)
               CASE(+1)
                  Z = ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%X2-X)/ABS(ESS%X1-ESS%X2)
               CASE(-2)
                  Z = ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%Y1-Y)/ABS(ESS%Y1-ESS%Y2)
               CASE(+2)
                  Z = ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%Y2-Y)/ABS(ESS%Y1-ESS%Y2)
               END SELECT
               EVAC_EXITS(I)%XYZ_Z = EVAC_EXITS(I)%XYZ_Z + Z
            END IF

            X = 0.5_EB*(EVAC_EXITS(I)%X1+EVAC_EXITS(I)%X2)
            Y = 0.5_EB*(EVAC_EXITS(I)%Y1+EVAC_EXITS(I)%Y2)
            Z = 0.0_EB
            IF ( (ESS%X1 <= X .AND. ESS%X2 >= X) .AND. (ESS%Y1 <= Y .AND. ESS%Y2 >= Y)) THEN
               SELECT CASE (ESS%IOR)
               CASE(-1)
                  Z = ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%X1-X)/ABS(ESS%X1-ESS%X2)
               CASE(+1)
                  Z = ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%X2-X)/ABS(ESS%X1-ESS%X2)
               CASE(-2)
                  Z = ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%Y1-Y)/ABS(ESS%Y1-ESS%Y2)
               CASE(+2)
                  Z = ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%Y2-Y)/ABS(ESS%Y1-ESS%Y2)
               END SELECT
               EVAC_EXITS(I)%Z = EVAC_EXITS(I)%Z + Z
            END IF
         END DO

         DO I = 1, N_DOORS
            IF (EVAC_DOORS(I)%IMESH /= ESS%IMESH) CYCLE
            X = EVAC_DOORS(I)%X
            Y = EVAC_DOORS(I)%Y
            Z = 0.0_EB
            IF ( (ESS%X1 <= X .AND. ESS%X2 >= X) .AND. (ESS%Y1 <= Y .AND. ESS%Y2 >= Y)) THEN
               SELECT CASE (ESS%IOR)
               CASE(-1)
                  Z = ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%X1-X)/ABS(ESS%X1-ESS%X2)
               CASE(+1)
                  Z = ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%X2-X)/ABS(ESS%X1-ESS%X2)
               CASE(-2)
                  Z = ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%Y1-Y)/ABS(ESS%Y1-ESS%Y2)
               CASE(+2)
                  Z = ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%Y2-Y)/ABS(ESS%Y1-ESS%Y2)
               END SELECT
               EVAC_DOORS(I)%XYZ_Z = EVAC_DOORS(I)%XYZ_Z + Z
            END IF
            X = 0.5_EB*(EVAC_DOORS(I)%X1+EVAC_DOORS(I)%X2)
            Y = 0.5_EB*(EVAC_DOORS(I)%Y1+EVAC_DOORS(I)%Y2)
            Z = 0.0_EB
            IF ( (ESS%X1 <= X .AND. ESS%X2 >= X) .AND. (ESS%Y1 <= Y .AND. ESS%Y2 >= Y)) THEN
               SELECT CASE (ESS%IOR)
               CASE(-1)
                  Z = ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%X1-X)/ABS(ESS%X1-ESS%X2)
               CASE(+1)
                  Z = ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%X2-X)/ABS(ESS%X1-ESS%X2)
               CASE(-2)
                  Z = ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%Y1-Y)/ABS(ESS%Y1-ESS%Y2)
               CASE(+2)
                  Z = ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%Y2-Y)/ABS(ESS%Y1-ESS%Y2)
               END SELECT
               EVAC_DOORS(I)%Z = EVAC_DOORS(I)%Z + Z
            END IF
         END DO

         DO I = 1, N_ENTRYS
            IF (EVAC_ENTRYS(I)%IMESH /= ESS%IMESH) CYCLE
            X = 0.5_EB*(EVAC_ENTRYS(I)%X1+EVAC_ENTRYS(I)%X2)
            Y = 0.5_EB*(EVAC_ENTRYS(I)%Y1+EVAC_ENTRYS(I)%Y2)
            Z = 0.0_EB
            IF ( (ESS%X1 <= X .AND. ESS%X2 >= X) .AND. (ESS%Y1 <= Y .AND. ESS%Y2 >= Y)) THEN
               SELECT CASE (ESS%IOR)
               CASE(-1)
                  Z = ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%X1-X)/ABS(ESS%X1-ESS%X2)
               CASE(+1)
                  Z = ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%X2-X)/ABS(ESS%X1-ESS%X2)
               CASE(-2)
                  Z = ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%Y1-Y)/ABS(ESS%Y1-ESS%Y2)
               CASE(+2)
                  Z = ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%Y2-Y)/ABS(ESS%Y1-ESS%Y2)
               END SELECT
               EVAC_ENTRYS(I)%Z = EVAC_ENTRYS(I)%Z + Z
            END IF
         END DO

      END DO READ_EVSS_LOOP
31    REWIND(LU_INPUT)

    END SUBROUTINE READ_EVSS

    SUBROUTINE CHECK_EVAC_NODES
      IMPLICIT NONE
      !
      ! Local variables
      REAL(EB) xtol,ytol
      LOGICAL :: L_TMP
      INTEGER, DIMENSION(:), ALLOCATABLE :: N_MESH_TMP
      TYPE (EVAC_ENTR_TYPE), POINTER :: PNX=>NULL()
      TYPE (EVAC_DOOR_TYPE), POINTER :: PDX=>NULL()
      TYPE (EVAC_STRS_TYPE), POINTER :: STRP=>NULL()
      !
      ! Set the IMESH and IMESH2 for corridors
      DO n = 1, n_corrs
         Nodeloop2: DO i = 1, n_nodes
            IF (EVAC_Node_List(i)%ID == EVAC_CORRS(n)%TO_NODE) THEN
               EVAC_CORRS(n)%INODE2 = i
               IF ( TRIM(EVAC_Node_List(i)%Node_Type) == 'Exit' .OR. &
                    TRIM(EVAC_Node_List(i)%Node_Type) == 'Door' .OR. &
                    TRIM(EVAC_Node_List(i)%Node_Type) == 'Entry' ) THEN
                  EVAC_CORRS(n)%IMESH  = EVAC_Node_List(i)%IMESH
                  EVAC_CORRS(n)%IMESH2 = EVAC_Node_List(i)%IMESH
               ELSE
                  EVAC_CORRS(n)%IMESH  = 0
                  EVAC_CORRS(n)%IMESH2 = n_egrids
               END IF
               EVAC_Node_List(evac_corrs(n)%INODE)%IMESH = EVAC_CORRS(n)%IMESH2
               EXIT Nodeloop2
            END IF
         END DO Nodeloop2
         IF (EVAC_CORRS(n)%INODE2 == 0 .OR. EVAC_CORRS(n)%IMESH2 == 0) THEN
            WRITE(MESSAGE,'(A,A,A,A)') 'ERROR: CORR ',Trim(EVAC_CORRS(n)%ID),' problem with TO_NODE ',&
                 Trim(EVAC_CORRS(n)%TO_NODE)
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END IF
      END DO

      ! Set the z offset counters for the EVHO and EVAC rectangles shown in SMV
      ALLOCATE(N_MESH_TMP(NMESHES))
      CALL ChkMemErr('CHECK_EVAC_NODES','N_MESH_TMP',IZERO)
      N_MESH_TMP = 0 ! Counter array for the EVAC and EVHO lines on the meshes
      DO N = 1, NPC_EVAC
         IF (.NOT.EVAC_EVACS(N)%SHOW) CYCLE
         ! IMESH: The mesh index, including fire meshes also, [1,NMESHES]
         N_MESH_TMP(EVAC_EVACS(N)%IMESH) = N_MESH_TMP(EVAC_EVACS(N)%IMESH) + 1
         EVAC_EVACS(N)%I_EVAC_THIS_MESH  = N_MESH_TMP(EVAC_EVACS(N)%IMESH)
      END DO
      DO N = 1, N_HOLES
         IF (.NOT.EVAC_HOLES(N)%SHOW) CYCLE
         N_MESH_TMP(EVAC_HOLES(N)%IMESH) = N_MESH_TMP(EVAC_HOLES(N)%IMESH) + 1
         EVAC_HOLES(N)%I_EVHO_THIS_MESH  = N_MESH_TMP(EVAC_HOLES(N)%IMESH)
      END DO
      DEALLOCATE(N_MESH_TMP)
      !
      DO n = 1, N_DOORS
         NodeLoop: DO i = 1, n_nodes
            IF (TRIM(EVAC_Node_List(i)%ID) == TRIM(EVAC_DOORS(n)%TO_NODE)) THEN
               EVAC_DOORS(n)%INODE2 = i
               EVAC_DOORS(n)%IMESH2 = EVAC_Node_List(i)%IMESH
               EXIT NodeLoop
            END IF
         END DO NodeLoop
         PDX_StrsLoop: DO I = 1, N_STRS
            ! N_EGRIDS: n_egrids_inputs + n_strs (other are additional meshes, evac_humans=.false.)
            IF (TRIM(EVAC_DOORS(n)%TO_NODE) == TRIM(EVAC_STRS(I)%ID)) THEN
               EVAC_DOORS(n)%IMESH2 = EVAC_STRS(I)%IMESH
               EVAC_DOORS(n)%INODE2 = N_EGRIDS - N_STRS + I
               EXIT PDX_StrsLoop
            END IF
         END DO PDX_StrsLoop
         IF (EVAC_DOORS(n)%INODE2 == 0 .OR. EVAC_DOORS(n)%IMESH2 == 0) THEN
            WRITE(MESSAGE,'(A,A,A,A)') 'ERROR: DOOR ',TRIM(EVAC_DOORS(n)%ID),' problem with TO_NODE ',&
                 TRIM(EVAC_DOORS(n)%TO_NODE)
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END IF
      END DO


      DO n = 1, N_DOORS
         IF (EVAC_DOORS(n)%KEEP_XY) THEN
            i = EVAC_DOORS(n)%INODE2
            IF (TRIM(EVAC_Node_List(i)%Node_Type) == 'Door') THEN
               PDX => EVAC_DOORS(EVAC_Node_List(i)%Node_Index)
               IF ((EVAC_DOORS(n)%IOR /= -PDX%IOR) .OR. ABS(EVAC_DOORS(n)%Width-PDX%Width) > 0.1_EB ) THEN
                  WRITE(MESSAGE,'(A,A,A)') 'ERROR: DOOR ',TRIM(EVAC_DOORS(n)%ID),' KEEP_XY Problem'
                  CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
               END IF
            END IF
            IF (TRIM(EVAC_Node_List(i)%Node_Type) == 'Entry') THEN
               PNX => EVAC_ENTRYS(EVAC_Node_List(i)%Node_Index)
               IF ((EVAC_DOORS(n)%IOR /= PNX%IOR) .OR. ABS(EVAC_DOORS(n)%Width-PNX%Width) > 0.1_EB ) THEN
                  WRITE(MESSAGE,'(A,A,A)') 'ERROR: DOOR ',TRIM(EVAC_DOORS(n)%ID),' KEEP_XY Problem'
                  CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
               END IF
            END IF
         END IF
      END DO

      DO N = 1,N_DOORS
         PDX => EVAC_DOORS(N)
         ! Check if door leads to Stairs
         PDX%STR_INDX = 0
         PDX%STR_SUB_INDX = 0
         IF (ABS(PDX%IOR)>1) THEN
            xtol=0.05_EB
            ytol=0.2_EB
         ELSE
            xtol=0.2_EB
            ytol=0.05_EB
         END IF
         CheckDoorStrLoop: DO I = 1, N_STRS
            STRP=>EVAC_STRS(I)
            IF (STRP%IMESH==PDX%IMESH .OR. STRP%IMESH==PDX%IMESH2) THEN
               PDX%STR_INDX = I
               DO J = 1,STRP%N_NODES
                  IF ( Is_Within_Bounds(PDX%X1,PDX%X2,PDX%Y1,PDX%Y2,PDX%Z1,PDX%Z2, &
                       STRP%XB_NODE(j,1), STRP%XB_NODE(j,2), STRP%XB_NODE(j,3),STRP%XB_NODE(j,4), &
                       STRP%XB_NODE(j,5), STRP%XB_NODE(j,6), xtol, ytol, 0.1_EB)) THEN
                     PDX%STR_SUB_INDX = j
                     EXIT CheckDoorStrLoop
                  END IF
               END DO
            END IF
         END DO CheckDoorStrLoop
         IF (MYID==MAX(0,EVAC_PROCESS) .AND. PDX%STR_INDX/=0 .AND. PDX%STR_SUB_INDX==0) THEN
            WRITE(MESSAGE,'(A,A,A,A,A)') 'ERROR: DOOR line ', TRIM(PDX%ID), ' is in stairs ', &
                 TRIM(EVAC_STRS(PDX%STR_INDX)%ID), ', no node found'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END IF
      END DO

      ! Create list of incoming nodes for Stairs
      DO N = 1,N_STRS
         STRP => EVAC_STRS(N)
         STRP%N_NODES_IN = 0
         STRP%N_NODES_OUT = 0
         NODES_TMP = 0
         DO I = 1,N_NODES
            SELECT CASE (EVAC_NODE_List(I)%Node_type)
            CASE ('Door')
               J = EVAC_NODE_List(I)%Node_index
               IF (EVAC_DOORS(J)%IMESH2 == STRP%IMESH) THEN
                  STRP%N_NODES_IN = STRP%N_NODES_IN + 1
                  NODES_TMP(STRP%N_NODES_IN) = I
               END IF
            CASE ('Entry')
               J = EVAC_NODE_List(I)%Node_index
               IF (EVAC_ENTRYS(J)%IMESH == STRP%IMESH) THEN
                  STRP%N_NODES_IN = STRP%N_NODES_IN + 1
                  NODES_TMP(STRP%N_NODES_IN) = I
               END IF
            END SELECT
         END DO
         ALLOCATE(STRP%NODES_IN(1:STRP%N_NODES_IN),STAT=IZERO)
         CALL ChkMemErr('Read_Evac','STRP%NODES_IN',IZERO)
         STRP%NODES_IN = NODES_TMP(1:STRP%N_NODES_IN)
         ! Create List of outgoing nodes for stairs
         NODES_TMP = 0
         DO I = 1, N_NODES
            IF (EVAC_NODE_List(I)%Node_type == 'Floor') CYCLE
            IF (EVAC_NODE_List(I)%Node_type == 'Entry') CYCLE
            IF (EVAC_NODE_List(I)%Node_type == 'Stair') CYCLE
            ! count only exits should not be counted here
            IF (EVAC_NODE_List(I)%Node_type == 'Exit') THEN
               IF (EVAC_EXITS(EVAC_NODE_List(I)%Node_Index)%COUNT_ONLY) CYCLE
            ENDIF
            IF (EVAC_NODE_List(I)%IMESH == STRP%IMESH) THEN
               STRP%N_NODES_OUT = STRP%N_NODES_OUT + 1
               NODES_TMP(STRP%N_NODES_OUT) = I
            END IF
         END DO
         ALLOCATE(STRP%NODES_OUT(1:STRP%N_NODES_OUT),STAT=IZERO)
         CALL ChkMemErr('Read_Evac','STRP%NODES_OUT',IZERO)
         STRP%NODES_OUT = NODES_TMP(1:STRP%N_NODES_OUT)
      END DO

      DO N = 1,N_EXITS
         L_TMP = .TRUE.
         DO i = 1, NPC_PERS
            IF (TRIM(EVAC_PERSON_CLASSES(i)%ID) == TRIM(EVAC_EXITS(N)%PERS_ID)) L_TMP=.FALSE.
         END DO
         IF (L_TMP .AND. .NOT.(TRIM(EVAC_EXITS(N)%PERS_ID))=='null') THEN
            WRITE(MESSAGE,'(4A)') 'ERROR: EXIT ',TRIM(EVAC_EXITS(N)%ID),' problem with PERS_ID ',TRIM(EVAC_EXITS(N)%PERS_ID)
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END IF
         L_TMP = .TRUE.
         DO i = 1, NPC_EVAC
            IF (TRIM(EVAC_EVACS(i)%ID) == TRIM(EVAC_EXITS(N)%EVAC_ID)) L_TMP=.FALSE.
         END DO
         IF (L_TMP .AND. .NOT.(TRIM(EVAC_EXITS(N)%EVAC_ID))=='null') THEN
            WRITE(MESSAGE,'(4A)') 'ERROR: EXIT ',TRIM(EVAC_EXITS(N)%ID),' problem with EVAC_ID ',TRIM(EVAC_EXITS(n)%EVAC_ID)
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END IF
      END DO

    END SUBROUTINE CHECK_EVAC_NODES

  END SUBROUTINE READ_EVAC

  SUBROUTINE INITIALIZE_EVAC_DUMPS(Tin,T_SAVE)
    IMPLICIT NONE
    !
    ! Passed variables
    REAL(EB), INTENT(IN) :: Tin
    REAL(EB), INTENT(INOUT) :: T_SAVE
    ! Local variables
    CHARACTER(50) tcform
    CHARACTER(LABEL_LENGTH) DEVC_ID
    INTEGER n_cols, i, j,k, nm, izero, j_ntargets, j_density, n_devc_read, IFIELD, NFIELDS, &
         JJ, JJ_NOW, N_END
    LOGICAL L_fed_read, L_fed_save, L_eff_read, L_eff_save, L_status, CURRENT_STATE, PRIOR_STATE
    INTEGER(4) n_egrids_tmp, ibar_tmp, jbar_tmp, kbar_tmp, &
         ntmp1, ntmp2, ntmp3, ntmp4, ntmp5, ntmp6, IOS, N
    REAL(FB) u_tmp, v_tmp
    INTEGER(4) N_TMP
    REAL(FB) TMPOUT1, TMPOUT2, TMPOUT3, TMPOUT4, T_TMP, DT_TMP
    REAL(FB) TMPOUT5, TMPOUT6, TMPOUT7, TMPOUT8
    CHARACTER(60), ALLOCATABLE, DIMENSION(:) :: CTEMP
    INTEGER :: size_rnd
    INTEGER, DIMENSION(8) :: t_rnd
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed_rnd
    CHARACTER(LABEL_LENGTH), DIMENSION(3) :: FED_ACT_NAME
    DATA FED_ACT_NAME /'at rest or at sleep', 'light work', 'heavy work'/
    !
    TYPE (MESH_TYPE), POINTER :: MFF =>NULL()
    TYPE (DEVICE_TYPE), POINTER :: DV =>NULL()
    TYPE (CONTROL_TYPE), POINTER :: CV =>NULL()
    TYPE (EVAC_EDEV_TYPE), POINTER :: EDV=>NULL()
    !
    ! Logical unit numbers
    ! LU_EVACCSV: CHID_evac.csv, number of persons
    ! LU_EVACEFF: CHID_evac.eff, evacflow fields, binary
    ! LU_EVACFED: CHID_evac.fed, FED and soot, time dependent, binary
    !      Format: 1. row: n_egrids >=0  (Old Format, version 1.10)
    !              1a. row: n < 0 (New Format)
    !              1b. row: n_egrids,4,n_corrs=0,4 (New Format, version 1.11)
    !      New Format, version 2.2.2: n=-2 ==> no fed infor for exits and doors

    IF (.NOT.ANY(EVACUATION_ONLY)) RETURN

    ! Evacuation files

    LU_EVACCSV = GET_FILE_NUMBER()
    FN_EVACCSV = TRIM(CHID)//'_evac.csv'
    LU_EVACEFF = GET_FILE_NUMBER()
    FN_EVACEFF = TRIM(CHID)//'_evac.eff'
    LU_EVACFED = GET_FILE_NUMBER()
    FN_EVACFED = TRIM(CHID)//'_evac.fed'
    LU_EVACOUT = GET_FILE_NUMBER()
    FN_EVACOUT = TRIM(CHID)//'_evac.out'
    LU_EVAC_CB = GET_FILE_NUMBER()

    L_fed_read = BTEST(I_EVAC,3)
    L_fed_save = BTEST(I_EVAC,1)
    L_eff_read = BTEST(I_EVAC,2)
    L_eff_save = BTEST(I_EVAC,0)

    !

    ! Open evacuation output file

    IF (APPEND) THEN
       OPEN (LU_EVACOUT,file=FN_EVACOUT,form='formatted',status='old', position='append')
    ELSE
       OPEN (LU_EVACOUT,file=FN_EVACOUT,form='formatted', status='replace')
    ENDIF

    ! Write program info

    WRITE(LU_EVACOUT,'(/A)')   ' FDS+Evac Evacuation Module'
    WRITE(LU_EVACOUT,'(A,A)')  ' FDS+Evac Version         : ', TRIM(EVAC_VERSION)
    WRITE(LU_EVACOUT,FMT='(/A,I2)')  ' FDS+Evac Color_Method    :', COLOR_METHOD
    IF (L_fed_read .OR. L_fed_save) THEN
       IF (Fed_Door_Crit >= 0) THEN
          WRITE(LU_EVACOUT,FMT='(A,F14.8)') ' FDS+Evac Fed_Door_Crit   :', FED_DOOR_CRIT
       ELSE
          ! Visibility S = 3/K, K is extinction coeff.
          WRITE(LU_EVACOUT,FMT='(A,F14.8,A)') ' FDS+Evac Vis_Door_Crit   :', ABS(FED_DOOR_CRIT), ' m'
       END IF
       WRITE(LU_EVACOUT,FMT='(A,A)') ' FDS+Evac: FED activity level is ', TRIM(FED_ACT_NAME(FED_ACTIVITY))
    ELSE
       WRITE(LU_EVACOUT,FMT='(A)') ' FDS+Evac: Evacuation drill mode '
    END IF
    IF (NOT_RANDOM ) WRITE(LU_EVACOUT,FMT='(A)') ' FDS+Evac Random seed is not used.'
    CALL RANDOM_SEED(size=size_rnd)
    ALLOCATE(seed_rnd(size_rnd),STAT=IZERO)
    CALL ChkMemErr('READ_EVAC','seed_rnd',IZERO)
    IF (.NOT. NOT_RANDOM) THEN    ! Initialize the generator randomly
       CALL DATE_AND_TIME(values = t_rnd)
       seed_rnd = 31*t_rnd(7) + 29*t_rnd(8)
    ELSE
       ! Do not use a random seed, use a constant seed
       seed_rnd = 2819
    END IF
    CALL RANDOM_SEED(put=seed_rnd)
    DEALLOCATE(seed_rnd)
    ! Fds uses a constant key to initialize the generator, but it is better for humans to use a random key.
    IF (Fed_Door_Crit < 0) THEN
       ! Visibility S = 3/K, K is extinction coeff.
       FED_DOOR_CRIT = 3.0_EB/FED_DOOR_CRIT ! Extinction coeff (1/m)
    END IF
    ! Next information should be saved in the FED-file so that evacuation calculations
    ! reading smoke from disk are using correct value.
    IF (SOOT_INDEX > 0) THEN
       EVAC_MASS_EXTINCTION_COEFF = SPECIES(SOOT_INDEX)%MASS_EXTINCTION_COEFFICIENT
       WRITE(LU_EVACOUT,FMT='(A,f12.2)') ' FDS+Evac: Soot mass extinction coeff is used: ',EVAC_MASS_EXTINCTION_COEFF
    ELSE
       EVAC_MASS_EXTINCTION_COEFF = 8700.0_EB  ! m2/kg
       WRITE(LU_EVACOUT,FMT='(A,f12.2)') ' FDS+Evac: Default mass extinction coeff is used: ',EVAC_MASS_EXTINCTION_COEFF
    END IF
    ! Write information about the falling down algorithm
    IF (L_FALLING_MODEL) THEN
       WRITE(LU_EVACOUT,FMT='(A)') ' FDS+Evac: Fallind down algoritm is used'
       WRITE(LU_EVACOUT,FMT='(A,f12.2)') '           TAU_FALL_DOWN  = ',TAU_FALL_DOWN
       WRITE(LU_EVACOUT,FMT='(A,f12.2)') '           D_OVERLAP_FALL = ',D_OVERLAP_FALL
       WRITE(LU_EVACOUT,FMT='(A,f12.2)') '           A_FAC_FALLEN   = ',A_FAC_FALLEN
       WRITE(LU_EVACOUT,FMT='(A,f12.2)') '           PROB_FALL_DOWN = ',PROB_FALL_DOWN
       WRITE(LU_EVACOUT,FMT='(A,f12.2)') '           F_MIN_FALL     = ',F_MIN_FALL
       WRITE(LU_EVACOUT,FMT='(A,f12.2)') '           F_MAX_FALL     = ',F_MAX_FALL
    END IF

    n_cols = n_egrids + n_corrs + N_EXITS + N_DOORS + 1 + N_EXITS - n_co_exits + N_DOORS
    ! Initialize the FED counters:
    icyc_old = -1
    n_dead = -1
    fed_max_alive = 0.0_EB
    fed_max = 0.0_EB
    !
    APPEND_IF: IF (APPEND) THEN
       OPEN (LU_EVACCSV,file=FN_EVACCSV,form='formatted',status='old', position='append')
       !
       FED_SAVE_RESTART: IF (L_fed_save) THEN
          OPEN (LU_EVACFED,file=FN_EVACFED,form='unformatted', status='old',position='rewind')
          READ (LU_EVACFED,IOSTAT=IOS) ntmp1
          IF (ios/=0) THEN
             WRITE(MESSAGE,'(A)') 'ERROR: Init Evac Dumps: FED READ ERROR'
             CLOSE (LU_EVACFED)
             CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
          END IF
          I_FED_FILE_FORMAT = ntmp1
          READ (LU_EVACFED,IOSTAT=IOS) n_egrids_tmp, ntmp2, ntmp3, ntmp4, ntmp5, ntmp6
          IF (ios/=0) THEN
             WRITE(MESSAGE,'(A)') 'ERROR: Init Evac Dumps: FED READ ERROR, Restart failed'
             CLOSE (LU_EVACFED)
             CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
          END IF
          ! Do not read old format. Do not read new format, if there the numbers are not: n_egrids, 4, n_corrs, 8
          IF (I_FED_FILE_FORMAT==-3) THEN  ! version 2.2.2 file format, no doors and exits
             IF ( ntmp2 /= 4 .OR. ntmp3 /= n_corrs .OR. ntmp1 >= 0 .OR. ntmp4 /= 8  .OR. &
                  ntmp5 /= 0 .OR. ntmp6 /= 4) THEN
                WRITE (MESSAGE,FMT='(a,a,a)') ' FDS+Evac Error in FED File: ', TRIM(FN_EVACFED), ', Restart failed'
                CLOSE (LU_EVACFED)
                CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
             END IF
          ELSE
             WRITE (MESSAGE,FMT='(a,a,a)') ' FDS+Evac Error in FED File: ', TRIM(FN_EVACFED), ', Old FED file format'
             CLOSE (LU_EVACFED)
             CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
          END IF
          IF (n_egrids_tmp /= n_egrids) THEN
             WRITE(MESSAGE,'(A,2I4,A)') 'ERROR: Init Evac Dumps: FED ',n_egrids_tmp, n_egrids, ', Restart failed'
             CLOSE (LU_EVACFED)
             CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
          END IF
          ! Next loop is for evacuation devices (like heat detectors)
          READ (LU_EVACFED,IOSTAT=IOS) ntmp1
          IF (IOS/=0) THEN
             WRITE(MESSAGE,'(A)') 'ERROR: Init Evac Dumps, Restart: FED read error1 for DEVC'
             CLOSE (LU_EVACFED)
             CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
          END IF
          ALLOCATE(EVAC_DEVICES(MAX(1,ntmp1)), STAT = IZERO)
          CALL ChkMemErr('INITIALIZE_EVAC_DUMPS','EVAC_DEVICES', IZERO)
          WRITE(LU_EVACOUT,'(A,I4)') ' FDS+Evac Number of evacuation devices written to the FED file: ',ntmp1
          N_DEVC_EVAC = 0
          DEVC_LOOP: DO I = 1, N_DEVC
             DV => DEVICE(I)
             IF (.NOT. DV%EVACUATION) CYCLE DEVC_LOOP
             IF (TRIM(DV%QUANTITY)=='TIME' .AND. DV%SETPOINT<=T_BEGIN) CYCLE DEVC_LOOP
             N_DEVC_EVAC = N_DEVC_EVAC + 1
             EVAC_DEVICES(N_DEVC_EVAC)%DEVC_ID     = DV%ID
             EVAC_DEVICES(N_DEVC_EVAC)%T_Change    = 0.0_EB
             EVAC_DEVICES(N_DEVC_EVAC)%CURRENT     = .FALSE.
             EVAC_DEVICES(N_DEVC_EVAC)%PRIOR       = .FALSE.
             EVAC_DEVICES(N_DEVC_EVAC)%I_Devc      = I
             EVAC_DEVICES(N_DEVC_EVAC)%I_Type      = DEVICE_INPUT  ! 1
             EVAC_DEVICES(N_DEVC_EVAC)%I_Devc_Evac = N_DEVC_EVAC
          END DO DEVC_LOOP
          CTRL_LOOP: DO I = 1, N_CTRL
             CV => CONTROL(I)
             IF (.NOT. CV%EVACUATION) CYCLE CTRL_LOOP
             N_DEVC_EVAC = N_DEVC_EVAC + 1
             EVAC_DEVICES(N_DEVC_EVAC)%DEVC_ID     = CV%ID
             EVAC_DEVICES(N_DEVC_EVAC)%T_Change    = 0.0_EB
             EVAC_DEVICES(N_DEVC_EVAC)%CURRENT     = .FALSE.
             EVAC_DEVICES(N_DEVC_EVAC)%PRIOR       = .FALSE.
             EVAC_DEVICES(N_DEVC_EVAC)%I_Devc      = I
             EVAC_DEVICES(N_DEVC_EVAC)%I_Type      = CONTROL_INPUT  ! 2
             EVAC_DEVICES(N_DEVC_EVAC)%I_Devc_Evac = N_DEVC_EVAC
          END DO CTRL_LOOP
          IF (ntmp1 /= N_DEVC_EVAC) THEN
             WRITE(MESSAGE,'(A)') 'ERROR: Init Evac Dumps, Restart: FED read error2 for DEVC'
             CLOSE (LU_EVACFED)
             CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
          END IF

          ! Position the FED file at the correct position, i.e., at the restart point.
          ! The FED file might have some time points after the restart time, because
          ! FED file is written every 2 seconds.
          ! LU_EVACFED: CHID_EVAC.FED, FED AND SOOT, TIME DEPENDENT, BINARY
          ! FILE FORMAT: 1A. ROW: N < 0 (NEW FORMAT)
          !              1B. ROW: N_EGRIDS,4,N_CORRS,8 (NEW FORMATS)
          !              1C. ROW: N_DEVC_EVAC (NEW FORMAT: -3)
          !                 2. ROW: T AND DT
          !                    3. ROW: IBAR,JBAR,KBAR, N_QUANTITIES
          !                       4. ROW: ONWARDS DATA
          !                    GOTO 3. (MESHES)
          !                       N. ROW: CORR DATA (8 REAL NUMBERS)
          !                       N+1. ROW: NEXT CORR DATA...
          !                 GOTO 2. (TIME POINTS)
          T_TMP  = REAL(T_BEGIN,FB)
          DT_TMP = 0.0_FB
          T_SAVE = 0.0_FB
          TIME_LOOP: DO WHILE (REAL(T_SAVE,EB) < Tin)
             IOS = 0
             READ (LU_EVACFED,END=324,IOSTAT=IOS) T_TMP, DT_TMP
             T_SAVE = T_TMP + DT_TMP ! Next time point in the file
             IF (IOS/=0) THEN
                WRITE(MESSAGE,'(A)') 'ERROR: Init Evac Dumps, Restart: FED read error'
                CLOSE (LU_EVACFED)
                CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
             END IF
             MESH_LOOP: DO NM=1,NMESHES
                IF ( .NOT.(EMESH_INDEX(NM)>0 .AND. EVACUATION_ONLY(NM)) ) CYCLE
                CALL POINT_TO_MESH(NM)
                READ (LU_EVACFED,IOSTAT=IOS) IBAR_TMP, JBAR_TMP, KBAR_TMP, N_TMP
                IF (IOS/=0) THEN
                   WRITE(MESSAGE,'(A)') 'ERROR: Init Evac Dumps, Restart: FED read error'
                   CLOSE (LU_EVACFED)
                   CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
                END IF
                IF (IBAR_TMP /= IBAR .OR. JBAR_TMP /= JBAR .OR. N_TMP < 4 ) THEN
                   WRITE(MESSAGE,'(A)') 'ERROR: Init Evac Dumps, Restart: Problems to read the FED file'
                   CLOSE (LU_EVACFED)
                   CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
                END IF
                DO I = 1, IBAR
                   DO J= 1, JBAR
                      READ (LU_EVACFED,IOSTAT=IOS) TMPOUT1, TMPOUT2, TMPOUT3, TMPOUT4
                      IF (IOS/=0) THEN
                         WRITE(MESSAGE,'(A)') 'ERROR: Init Evac Dumps, Restart: FED read error'
                         CLOSE (LU_EVACFED)
                         CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
                      END IF
                   END DO     ! J=1,JBAR
                END DO       ! I=1,IBAR
             END DO MESH_LOOP
             CORR_LOOP: DO I = 1, N_CORRS
                READ (LU_EVACFED,IOSTAT=IOS) TMPOUT1, TMPOUT2, TMPOUT3, TMPOUT4, TMPOUT5, TMPOUT6, TMPOUT7, TMPOUT8
                IF (IOS/=0) THEN
                   WRITE(MESSAGE,'(A)') 'ERROR: Init Evac Dumps, Restart: FED read error'
                   CLOSE (LU_EVACFED)
                   CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
                END IF
             END DO CORR_LOOP
             READ (LU_EVACFED,IOSTAT=IOS) ntmp1
             IF (IOS/=0) THEN
                WRITE(MESSAGE,'(A)') 'ERROR: Init Evac Dumps, Restart: FED read error3 for DEVC'
                CLOSE (LU_EVACFED)
                CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
             END IF
             n_devc_read = ntmp1
             DO I = 1, n_devc_read
                READ (LU_EVACFED,IOSTAT=IOS) ntmp3, DEVC_ID
                IF (IOS/=0) THEN
                   WRITE(MESSAGE,'(A)') 'ERROR: Init Evac Dumps, Restart: FED read error4 for DEVC'
                   CLOSE (LU_EVACFED)
                   CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
                END IF
                READ (LU_EVACFED,IOSTAT=IOS) ntmp1, ntmp2, CURRENT_STATE, PRIOR_STATE, tmpout1
                IF (IOS /= 0 .OR. ntmp1 > N_DEVC_EVAC) THEN
                   WRITE(MESSAGE,'(A)') 'ERROR: Init Evac Dumps, Restart: FED read error5 for DEVC'
                   CLOSE (LU_EVACFED)
                   CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
                END IF
                j = ntmp1
                IF (ntmp2 /= EVAC_DEVICES(j)%I_Devc .OR. TRIM(DEVC_ID) /= TRIM(EVAC_DEVICES(j)%DEVC_ID) .OR. &
                     ntmp3 /= EVAC_DEVICES(j)%I_Type) THEN
                   WRITE(MESSAGE,'(A)') 'ERROR: Init Evac Dumps, Restart: FED read error6 for DEVC'
                   CLOSE (LU_EVACFED)
                   CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
                END IF
                EVAC_DEVICES(j)%T_Change = tmpout1
                EVAC_DEVICES(j)%CURRENT = CURRENT_STATE
                EVAC_DEVICES(j)%PRIOR = PRIOR_STATE
             END DO
          END DO TIME_LOOP
324       CONTINUE

          WRITE (LU_EVACOUT,fmt='(a,a,a)') ' FDS+Evac FED File: ', TRIM(FN_EVACFED), ' is calculated and used'
       END IF FED_SAVE_RESTART
       !
       FED_READ_RESTART: IF (L_fed_read) THEN
          INQUIRE (file=FN_EVACFED,exist=L_status)
          IF (.NOT. L_status) THEN
             WRITE (LU_EVACOUT,fmt='(a,a,a)') ' FDS+Evac No FED File: ', TRIM(FN_EVACFED), ', FED and soot not used'
             l_fed_read = .FALSE.
             l_fed_save = .FALSE.
             I_EVAC = IBCLR(I_EVAC,3)  ! do not read FED
             I_EVAC = IBCLR(I_EVAC,1)  ! do not save FED
          ELSE
             CALL SHUTDOWN('ERROR: Evac Dumps: FED, no restart yet',PROCESS_0_ONLY=.FALSE.) ; RETURN
             OPEN (LU_EVACFED,file=FN_EVACFED,form='unformatted', status='old')
             READ (LU_EVACFED,Iostat=ios) n_egrids_tmp
             IF (ios/=0) THEN
                WRITE(MESSAGE,'(A)') 'ERROR: Init Evac Dumps: FED READ ERROR'
                CLOSE (LU_EVACFED)
                CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
             END IF
             IF (n_egrids_tmp /= n_egrids) THEN
                WRITE(MESSAGE,'(A,2I4)') 'ERROR: Init Evac Dumps: FED ',n_egrids_tmp, n_egrids
                CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
             END IF
          END IF
       END IF FED_READ_RESTART
       EFF_READ_RESTART: IF (L_eff_read) THEN
          INQUIRE (file=FN_EVACEFF,exist=L_status)
          IF (L_status) THEN
             WRITE (LU_EVACOUT,fmt='(a,a,a/)') ' FDS+Evac EFF File: ', TRIM(FN_EVACEFF), ' is used'
             l_eff_save = .FALSE.
             I_EVAC = IBCLR(I_EVAC,0)  ! do not save EFF
             OPEN (LU_EVACEFF,file=FN_EVACEFF,form='unformatted', status='old')
             READ (LU_EVACEFF,Iostat=ios) n_egrids_tmp
             IF (ios/=0) THEN
                WRITE(MESSAGE,'(A)') 'ERROR: Init Evac Dumps, Restart: EFF READ ERROR'
                CLOSE (LU_EVACEFF)
                CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
             END IF
             IF (N_EGRIDS_TMP /= SUM(EMESH_NFIELDS)) THEN
                WRITE(MESSAGE,'(A,2I4)') 'ERROR: Init Evac Dumps, Restart: EFF ',N_EGRIDS_TMP, SUM(EMESH_NFIELDS)
                CLOSE (LU_EVACEFF)
                CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
             END IF
          ELSE
             WRITE(MESSAGE,'(A,2I4)') 'ERROR: Init Evac Dumps: EFF, no restart yet'
             CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
          END IF
       END IF EFF_READ_RESTART
       !
       IF ( l_fed_read .OR. l_fed_save ) n_dead = 0
       !
       ! Restart: read always from the hard drive
       ios = 3
       INQUIRE (file=FN_EVACEFF,exist=L_status)
       IF (L_status) THEN
          l_eff_save = .FALSE.
          l_eff_read = .TRUE.
          I_EVAC = IBCLR(I_EVAC,0) ! do not save EFF
          I_EVAC = IBSET(I_EVAC,2) ! read EFF
          OPEN (LU_EVACEFF,file=FN_EVACEFF,form='unformatted', status='old')
          READ (LU_EVACEFF,IOSTAT=IOS) N_EGRIDS_TMP
          IF (IOS/=0) THEN
             IOS = 1
             WRITE(LU_EVACOUT,'(A)') ' WARNING: Init Evac Dumps: EFF READ ERROR'
             WRITE(LU_EVACOUT,'(A)') ' WARNING: EFF file is not read in'
             CLOSE (LU_EVACEFF)
          END IF
          IF (N_EGRIDS_TMP /=  SUM(EMESH_NFIELDS) .AND. IOS < 1) THEN
             IOS = 2
             WRITE(LU_EVACOUT,'(A,2I4)') ' WARNING: Init Evac Dumps: EFF READ ERROR ',n_egrids_tmp, SUM(EMESH_NFIELDS)
             WRITE(LU_EVACOUT,'(A)')     ' WARNING: EFF file is not read in'
             CLOSE (LU_EVACEFF)
          END IF
       END IF
       IF (IOS /= 0) THEN
         WRITE(MESSAGE,'(A)') 'ERROR: Restart problem: EFF READ ERROR'
         CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
       END IF

    ELSE ! not a restart, rewrite files or just read in
       !
       FED_SAVE: IF (L_fed_save) THEN
          l_fed_read = .FALSE.
          I_EVAC = IBCLR(I_EVAC,3)  ! do not read FED
          OPEN (LU_EVACFED,file=FN_EVACFED,form='unformatted', status='replace')
          ! First line: <0 new format
          !             -1: second line: #mesh #reals #corrs #reals #doors+exits #nreals
          !              (#reals: fed,soot,temp,radflux,...)
          ! First line: >0: nmeshes, fed and soot saved/read for meshes
          !
          ! New format -1: fed information is save for doors and all exits (also for count_only)
          ! New format -2: fed information is not saved for doors.
          ntmp1 = I_FED_FILE_FORMAT
          ntmp2 = 4
          ntmp3 = N_CORRS
          ! Corrs: save for both XB1 and XB2 (if only XB, XB2 is then zeros)
          ntmp4 = 8
          ntmp5 = 0
          ntmp6 = 4
          n_egrids_tmp = n_egrids
          WRITE (LU_EVACFED) ntmp1
          WRITE (LU_EVACFED) n_egrids_tmp, ntmp2, ntmp3, ntmp4, ntmp5, ntmp6
          N_DEVC_EVAC = 0
          DEVC_LOOP_WRITE: DO I = 1, N_DEVC
             DV => DEVICE(I)
             IF (.NOT. DV%EVACUATION) CYCLE DEVC_LOOP_WRITE
             IF (TRIM(DV%QUANTITY)=='TIME' .AND. DV%SETPOINT<=T_BEGIN) CYCLE DEVC_LOOP_WRITE
             N_DEVC_EVAC = N_DEVC_EVAC + 1
          END DO DEVC_LOOP_WRITE
          CTRL_LOOP_WRITE: DO I = 1, N_CTRL
             CV => CONTROL(I)
             IF (.NOT. CV%EVACUATION) CYCLE CTRL_LOOP_WRITE
             N_DEVC_EVAC = N_DEVC_EVAC + 1
          END DO CTRL_LOOP_WRITE
          ntmp1 = N_DEVC_EVAC
          WRITE (LU_EVACFED) ntmp1
          WRITE(LU_EVACOUT,'(A,I4)') ' FDS+Evac Number of evacuation devices written to the FED file: ',N_DEVC_EVAC
          ALLOCATE(EVAC_DEVICES(MAX(1,N_DEVC_EVAC)), STAT = IZERO)
          CALL ChkMemErr('INITIALIZE_EVAC_DUMPS','EVAC_DEVICES', IZERO)
          N_DEVC_EVAC = 0
          DEVC_LOOP_WRITE_2: DO I = 1, N_DEVC
             DV => DEVICE(I)
             IF (.NOT. DV%EVACUATION) CYCLE DEVC_LOOP_WRITE_2
             IF (TRIM(DV%QUANTITY)=='TIME' .AND. DV%SETPOINT<=T_BEGIN) CYCLE DEVC_LOOP_WRITE_2
             N_DEVC_EVAC = N_DEVC_EVAC + 1
             EVAC_DEVICES(N_DEVC_EVAC)%DEVC_ID     = DV%ID
             EVAC_DEVICES(N_DEVC_EVAC)%T_Change    = DV%T_CHANGE
             EVAC_DEVICES(N_DEVC_EVAC)%CURRENT     = DV%CURRENT_STATE
             EVAC_DEVICES(N_DEVC_EVAC)%PRIOR       = DV%PRIOR_STATE
             EVAC_DEVICES(N_DEVC_EVAC)%I_Devc      = I
             EVAC_DEVICES(N_DEVC_EVAC)%I_Type      = DEVICE_INPUT  ! 1
             EVAC_DEVICES(N_DEVC_EVAC)%I_Devc_Evac = N_DEVC_EVAC
          END DO DEVC_LOOP_WRITE_2
          CTRL_LOOP_WRITE_2: DO I = 1, N_CTRL
             CV => CONTROL(I)
             IF (.NOT. DV%EVACUATION) CYCLE CTRL_LOOP_WRITE_2
             N_DEVC_EVAC = N_DEVC_EVAC + 1
             EVAC_DEVICES(N_DEVC_EVAC)%DEVC_ID     = CV%ID
             EVAC_DEVICES(N_DEVC_EVAC)%T_Change    = CV%T_CHANGE
             EVAC_DEVICES(N_DEVC_EVAC)%CURRENT     = CV%CURRENT_STATE
             EVAC_DEVICES(N_DEVC_EVAC)%PRIOR       = CV%PRIOR_STATE
             EVAC_DEVICES(N_DEVC_EVAC)%I_Devc      = I
             EVAC_DEVICES(N_DEVC_EVAC)%I_Type      = CONTROL_INPUT  ! 2
             EVAC_DEVICES(N_DEVC_EVAC)%I_Devc_Evac = N_DEVC_EVAC
          END DO CTRL_LOOP_WRITE_2
          ! Check the inputs for the EDEV namelist array EVAC_EDEV(1:N_EDEV)
          EDEV_LOOP: DO I = 1, N_EDEV
             EDV => EVAC_EDEV(I)
             EDV%INPUT_DEVC_INDEX(:) = 0
             INPUTS_LOOP: DO J = 1, EDV%N_INPUTS
                DEVICES_LOOP: DO K = 1, N_DEVC_EVAC
                   IF (TRIM(EDV%INPUT_ID(J)) == TRIM(EVAC_DEVICES(K)%DEVC_ID)) THEN
                      EDV%INPUT_DEVC_INDEX(J) = K
                      CYCLE INPUTS_LOOP
                   END IF
                END DO DEVICES_LOOP
                WRITE(MESSAGE,'(A,I5,A,A,A)')  'ERROR: EDEV ',I,' input ',TRIM(EDV%INPUT_ID(J)),' is not found'
                CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
             END DO INPUTS_LOOP
          END DO EDEV_LOOP

          WRITE (LU_EVACOUT,fmt='(a,a,a)') ' FDS+Evac FED File: ', TRIM(FN_EVACFED), ' is calculated and used'
       END IF FED_SAVE
       !
       FED_READ: IF (L_fed_read) THEN
          INQUIRE (file=FN_EVACFED,exist=L_status)
          NO_FED_AT_DISK: IF (.NOT. L_status) THEN
             WRITE (LU_EVACOUT,fmt='(a,a,a)') ' FDS+Evac No FED File: ', TRIM(FN_EVACFED), ', FED and soot not used'
             l_fed_read = .FALSE.
             l_fed_save = .FALSE.
             I_EVAC = IBCLR(I_EVAC,3)  ! do not read FED
             I_EVAC = IBCLR(I_EVAC,1)  ! do not save FED
          ELSE
             OPEN (LU_EVACFED,file=FN_EVACFED,form='unformatted', status='old')
             READ (LU_EVACFED,Iostat=ios) ntmp1
             IF (ios/=0) THEN
                WRITE(MESSAGE,'(A)') 'ERROR: Init Evac Dumps: FED READ ERROR'
                CLOSE (LU_EVACFED)
                CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
             END IF
             I_FED_FILE_FORMAT = ntmp1
             IF (I_FED_FILE_FORMAT/=-3) THEN
                WRITE(MESSAGE,'(A)') 'ERROR: Init Evac Dumps: FED READ ERROR, WRONG FORAT FOR FED FILE'
                CLOSE (LU_EVACFED)
                CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
             END IF
             ! New format (version 1.11)
             READ (LU_EVACFED,Iostat=ios) n_egrids_tmp, ntmp2, ntmp3, ntmp4, ntmp5, ntmp6
             IF (ios/=0) THEN
                WRITE(MESSAGE,'(A)') 'ERROR: Init Evac Dumps: FED READ ERROR'
                CLOSE (LU_EVACFED)
                CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
             END IF

             ! Do not read old format. Do not read new format, if there the numbers are not: n_egrids, 4, n_corrs, 8
             IF ( ntmp2 /= 4 .OR. ntmp3 /= n_corrs .OR. ntmp1 >= 0 .OR. ntmp4 /= 8  .OR. &
                  ntmp5 /= 0 .OR. ntmp6 /= 4) THEN
                WRITE (LU_EVACOUT,fmt='(a,a,a)') ' FDS+Evac Error in FED File: ', TRIM(FN_EVACFED), ', FED and soot not used'
                l_fed_read = .FALSE.
                l_fed_save = .FALSE.
                I_EVAC = IBCLR(I_EVAC,3) ! do not read FED
                I_EVAC = IBCLR(I_EVAC,1) ! do not save FED
                CLOSE (LU_EVACFED)
             END IF
             IF (n_egrids_tmp /= n_egrids) THEN
                WRITE(MESSAGE,'(A,2I4)') 'ERROR: Init Evac Dumps: FED ',n_egrids_tmp, n_egrids
                CLOSE (LU_EVACFED)
                CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
             END IF
             ! Next is for evacuation devices (like heat detectors)
             READ (LU_EVACFED,IOSTAT=IOS) ntmp1
             IF (IOS/=0) THEN
                WRITE(MESSAGE,'(A)') 'ERROR: Init Evac Dumps: FED read error1 for DEVC'
                CLOSE (LU_EVACFED)
                CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
             END IF
             N_DEVC_EVAC = ntmp1
             WRITE(LU_EVACOUT,'(A,I4)') ' FDS+Evac Number of evacuation devices in the FED file: ',N_DEVC_EVAC
             ALLOCATE(EVAC_DEVICES(MAX(1,N_DEVC_EVAC)), STAT = IZERO)
             CALL ChkMemErr('INITIALIZE_EVAC_DUMPS','EVAC_DEVICES', IZERO)

             IF (l_fed_read .OR. l_fed_save) THEN
                WRITE (LU_EVACOUT,fmt='(a,a,a)') ' FDS+Evac FED File: ', TRIM(FN_EVACFED), ' is used'
             END IF
          END IF NO_FED_AT_DISK
       END IF FED_READ
       !
       ! Number of evac flow fields is same as the number of all evac grids.
       EFF_READ: IF (L_eff_read) THEN
          ios = 3
          INQUIRE (file=FN_EVACEFF,exist=L_status)
          IF (L_status) THEN
             l_eff_save = .FALSE.
             l_eff_read = .TRUE.
             I_EVAC = IBCLR(I_EVAC,0) ! do not save EFF
             I_EVAC = IBSET(I_EVAC,2) ! read EFF
             OPEN (LU_EVACEFF,file=FN_EVACEFF,form='unformatted', status='old')
             READ (LU_EVACEFF,IOSTAT=IOS) N_EGRIDS_TMP
             IF (IOS/=0) THEN
                IOS = 1
                WRITE(LU_EVACOUT,'(A)') ' WARNING: Init Evac Dumps: EFF READ ERROR'
                WRITE(LU_EVACOUT,'(A)') ' WARNING: EFF file is not read in'
                CLOSE (LU_EVACEFF)
             END IF
             IF (N_EGRIDS_TMP /= SUM(EMESH_NFIELDS) .AND. IOS < 1) THEN
                IOS = 2
                WRITE(LU_EVACOUT,'(A,2I4)') ' WARNING: Init Evac Dumps: EFF READ ERROR ',n_egrids_tmp, SUM(EMESH_NFIELDS)
                WRITE(LU_EVACOUT,'(A)')     ' WARNING: EFF file is not read in'
                CLOSE (LU_EVACEFF)
             END IF
          END IF
          IF (IOS /= 0) THEN
             ! Read error ==> recalculate EFF
             L_EFF_SAVE = .TRUE.
             L_EFF_READ = .FALSE.
             I_EVAC = IBCLR(I_EVAC,2) ! do not read EFF
             I_EVAC = IBSET(I_EVAC,0) ! save EFF
          END IF
       END IF EFF_READ
       !
       ALLOCATE(CTEMP(MAX(1,N_EXITS)), STAT = IZERO)
       CALL ChkMemErr('INITIALIZE_EVAC_DUMPS','CTEMP', IZERO)
       j = 0
       DO i = 1, N_EXITS
          IF (.NOT. EVAC_EXITS(i)%COUNT_ONLY) THEN
             j = j + 1
             CTEMP(j) = TRIM(EVAC_EXITS(i)%ID)
          END IF
       END DO
       j_ntargets = j
       DO i = 1, N_EXITS
          IF (EVAC_EXITS(i)%COUNT_DENSITY .AND. EVAC_EXITS(i)%COUNT_ONLY) THEN
             j = j + 1
             CTEMP(j) = TRIM(EVAC_EXITS(i)%ID)
          END IF
       END DO
       j_density = j

       IF ( l_fed_read .OR. l_fed_save .OR. L_FALLING_MODEL) THEN
          ! Write the 'fed' columns
          n_dead = 0
          OPEN (LU_EVACCSV,file=FN_EVACCSV,form='formatted',status='replace')
          ! June 2009: Changed the .csv file format to the fds5 style
          ! first row: units (or variable class)
          ! second row: variable name
          ! third row-: data
          WRITE (tcform,'(a,i4.4,a)') "(",n_cols+3+(j_density-j_ntargets),"(a,','),a)"
          WRITE (LU_EVACCSV,tcform) 's','AgentsInside', &
               ('AgentsInsideMesh', i=1,n_egrids), &
               ('AgentsInsideCorr', i=1,n_corrs), &
               ('ExitCounter', i=1,N_EXITS), &
               ('DoorCounter', i=1,N_DOORS), &
               ('TargetExitCounter', i=1,N_EXITS-n_co_exits), &
               ('TargetDoorCounter', i=1,N_DOORS), &
               ('DensityCounter', i=1,j_density-j_ntargets), &
               'Agents','FED_Index','FED_Index'
          WRITE (LU_EVACCSV,tcform) 'EVAC_Time','AllAgents', &
               (TRIM(EVAC_Node_List(i)%GRID_NAME), i=1,n_egrids), &
               (TRIM(EVAC_CORRS(i)%ID), i=1,n_corrs), &
               (TRIM(EVAC_EXITS(i)%ID), i=1,N_EXITS), &
               (TRIM(EVAC_DOORS(i)%ID), i=1,N_DOORS), &
               (TRIM(CTEMP(i)), i=1,N_EXITS-n_co_exits), &
               (TRIM(EVAC_DOORS(i)%ID), i=1,N_DOORS), &
               (TRIM(CTEMP(i)), i=j_ntargets+1,j_density), &
               'Number_of_Deads','FED_max','FED_max_alive'
       ELSE
          ! Do not write the 'fed' columns
          OPEN (LU_EVACCSV,file=FN_EVACCSV,form='formatted',status='replace')
          ! June 2009: Changed the .csv file format to the fds5 style
          ! first row: units (or variable class)
          ! second row: variable name
          ! third row-: data
          IF (C_HAWK >= 0.0_EB) THEN
             WRITE (tcform,'(a,i4.4,a,i4.4,a)') "(",n_cols+(j_density-j_ntargets),"(a,','),a,", n_egrids*4, "(',',a))"
             WRITE (LU_EVACCSV,tcform) 's','AgentsInside', &
                  ('AgentsInsideMesh', i=1,n_egrids), &
                  ('AgentsInsideCorr', i=1,n_corrs), &
                  ('ExitCounter', i=1,N_EXITS), &
                  ('DoorCounter', i=1,N_DOORS), &
                  ('TargetExitCounter', i=1,N_EXITS-n_co_exits), &
                  ('TargetDoorCounter', i=1,N_DOORS), &
                  ('DensityCounter', i=1,j_density-j_ntargets), &
                  (('HawkDoveCounter',j=1,4), i=1,n_egrids)
             WRITE (LU_EVACCSV,tcform) 'EVAC_Time','AllAgents', &
                  (TRIM(EVAC_Node_List(i)%GRID_NAME), i=1,n_egrids), &
                  (TRIM(EVAC_CORRS(i)%ID), i=1,n_corrs), &
                  (TRIM(EVAC_EXITS(i)%ID), i=1,N_EXITS), &
                  (TRIM(EVAC_DOORS(i)%ID), i=1,N_DOORS), &
                  (TRIM(CTEMP(i)), i=1,N_EXITS-n_co_exits), &
                  (TRIM(EVAC_DOORS(i)%ID), i=1,N_DOORS), &
                  (TRIM(CTEMP(i)), i=j_ntargets+1,j_density), &
                  ('Hawk1','Dove1','Hawk2','Dove2', i=1,n_egrids)
          ELSE
             WRITE (tcform,'(a,i4.4,a)') "(",n_cols+(j_density-j_ntargets),"(a,','),a)"
             WRITE (LU_EVACCSV,tcform) 's','AgentsInside', &
                  ('AgentsInsideMesh', i=1,n_egrids), &
                  ('AgentsInsideCorr', i=1,n_corrs), &
                  ('ExitCounter', i=1,N_EXITS), &
                  ('DoorCounter', i=1,N_DOORS), &
                  ('TargetExitCounter', i=1,N_EXITS-n_co_exits), &
                  ('TargetDoorCounter', i=1,N_DOORS), &
                  ('DensityCounter', i=1,j_density-j_ntargets)
             WRITE (LU_EVACCSV,tcform) 'EVAC_Time','AllAgents', &
                  (TRIM(EVAC_Node_List(i)%GRID_NAME), i=1,n_egrids), &
                  (TRIM(EVAC_CORRS(i)%ID), i=1,n_corrs), &
                  (TRIM(EVAC_EXITS(i)%ID), i=1,N_EXITS), &
                  (TRIM(EVAC_DOORS(i)%ID), i=1,N_DOORS), &
                  (TRIM(CTEMP(i)), i=1,N_EXITS-n_co_exits), &
                  (TRIM(EVAC_DOORS(i)%ID), i=1,N_DOORS), &
                  (TRIM(CTEMP(i)), i=j_ntargets+1,j_density)
          END IF
       END IF
       DEALLOCATE(CTEMP)
    END IF APPEND_IF
    !
    ! Read the evac flow fields from the disk, if they exist.
    IF (L_eff_read .OR. APPEND) THEN
       ios = 0
       N_END = N_EXITS - N_CO_EXITS + N_DOORS
       ReadEffLoop: DO NM = 1, NMESHES
          EVAC_ONLY_NM: IF (EVACUATION_ONLY(NM)) THEN
             !IF (.NOT. EVACUATION_GRID(NM)) CYCLE ReadEffLoop
             IF (EMESH_INDEX(NM)==0) CYCLE ReadEffLoop
             MFF=>MESHES(NM)
             NFIELDS = EMESH_NFIELDS(EMESH_INDEX(NM))
             DO IFIELD = 1, NFIELDS

                LOOP_EXITS: DO JJ = 1, N_END
                   IF (EMESH_EXITS(JJ)%MAINMESH == NM .AND. EMESH_EXITS(JJ)%I_DOORS_EMESH == IFIELD) THEN
                      JJ_NOW = JJ
                      EXIT LOOP_EXITS
                   END IF
                END DO LOOP_EXITS

                READ (LU_EVACEFF,Iostat=ios) ibar_tmp, jbar_tmp, kbar_tmp
                IF (IOS /= 0) THEN
                   WRITE (LU_EVACOUT,'(A)') ' WARNING: Init Evac Dumps: EFF READ ERROR'
                   WRITE (LU_EVACOUT,'(A)') ' WARNING: EFF file is not read in'
                   CLOSE (LU_EVACEFF)
                   EXIT ReadEffLoop
                END IF
                IF ( MFF%IBAR /= ibar_tmp .OR. MFF%JBAR /= jbar_tmp .OR. MFF%KBAR /= kbar_tmp ) THEN
                   ios = 2
                   WRITE (LU_EVACOUT,'(A)') ' WARNING: Init Evac Dumps: EFF READ ERROR'
                   WRITE (LU_EVACOUT,'(A)') ' WARNING: EFF file is not read in'
                   CLOSE (LU_EVACEFF)
                   EXIT ReadEffLoop
                END IF
                DO  i = 0, MFF%IBAR+1
                   DO j= 0, MFF%JBAR+1
                      READ (LU_EVACEFF,Iostat=ios) u_tmp, v_tmp
                      IF (ios /= 0) THEN
                         WRITE (LU_EVACOUT,'(A)') ' WARNING: Init Evac Dumps: EFF READ ERROR'
                         WRITE (LU_EVACOUT,'(A)') ' WARNING: EFF file is not read in'
                         CLOSE (LU_EVACEFF)
                         EXIT ReadEffLoop
                      END IF
                      EMESH_EXITS(JJ_NOW)%U_EVAC(I,J) = REAL(u_tmp,EB)
                      EMESH_EXITS(JJ_NOW)%V_EVAC(I,J) = REAL(v_tmp,EB)
                      MFF%U(i,j,:) = REAL(u_tmp,EB)
                      MFF%V(i,j,:) = REAL(v_tmp,EB)
                      MFF%W(i,j,:) = 0.0_EB
                      MFF%KRES(i,j,:) = 0.5_EB*SQRT(REAL(u_tmp,EB)**2 + REAL(v_tmp,EB)**2)
                   END DO
                END DO
                MFF%U_EDGE_Y(:)=-1.E6_EB
                MFF%U_EDGE_Z(:)=-1.E6_EB
                MFF%V_EDGE_X(:)=-1.E6_EB
                MFF%V_EDGE_Z(:)=-1.E6_EB
             END DO
          END IF EVAC_ONLY_NM
       END DO ReadEffLoop
       IF (ios /= 0) THEN
          ! Read error ==> recalculate EFF
          l_eff_save = .TRUE.
          l_eff_read = .FALSE.
          I_EVAC = IBCLR(I_EVAC,2) ! do not read EFF
          I_EVAC = IBSET(I_EVAC,0) ! save EFF
       END IF
    END IF
    !
    IF (L_eff_read) THEN
       WRITE (LU_EVACOUT,fmt='(a,a,a/)') ' FDS+Evac EFF File: ', TRIM(FN_EVACEFF), ' is read in and used'
    END IF
    IF (L_eff_save .AND. .NOT. APPEND) THEN
       l_eff_read = .FALSE.
       I_EVAC = IBCLR(I_EVAC,2)  ! do not read EFF
       OPEN (LU_EVACEFF,file=FN_EVACEFF,form='unformatted', status='replace')
       n_egrids_tmp = SUM(EMESH_NFIELDS)
       WRITE (LU_EVACEFF) n_egrids_tmp
       WRITE (LU_EVACOUT,fmt='(a,a,a/)') ' FDS+Evac EFF File: ', TRIM(FN_EVACEFF), ' is (re)calculated and used'
    END IF

    EVAC_Z_MIN =  HUGE(EVAC_Z_MIN)
    EVAC_Z_MAX = -HUGE(EVAC_Z_MIN)
    DO nm = 1, NMESHES
       MFF=>MESHES(nm)
       EVAC_Z_MIN = MIN(EVAC_Z_MIN,REAL(MFF%ZS,FB))
       EVAC_Z_MAX = MAX(EVAC_Z_MAX,REAL(MFF%ZF,FB))
    END DO

    ! Write STRS properties
    IF (MYID==MAX(0,EVAC_PROCESS)) THEN
       DO N = 1, N_STRS
          WRITE (LU_EVACOUT,'(A,A)')      '  Stair ',TRIM(EVAC_STRS(N)%ID)
          WRITE (LU_EVACOUT,'(A,6F10.3)') '   Co-ordinates: ',EVAC_STRS(N)%XB(1:6)
          WRITE (LU_EVACOUT,'(A)')          '   Node coordinates'
          DO NM = 1,EVAC_STRS(N)%N_NODES
             IF (EVAC_STRS(N)%NODE_TYPE(NM) == STRS_LANDING_TYPE) THEN
                WRITE (LU_EVACOUT,'(A,I3,8F8.3)') '   Landing ',NM,(EVAC_STRS(N)%XB_NODE(NM,I),I=1,8)
             ELSE
                WRITE (LU_EVACOUT,'(A,I3,8F8.3)') '   Stair   ',NM,(EVAC_STRS(N)%XB_NODE(NM,I),I=1,8)
             ENDIF
          ENDDO
          WRITE (LU_EVACOUT,'(A)')          '   Nodes in '
          DO NM = 1, EVAC_STRS(N)%N_NODES_IN
             WRITE (LU_EVACOUT,'(I5,A,A)')          NM, ' ', Trim(EVAC_NODE_List(EVAC_STRS(N)%NODES_IN(NM))%ID)
          ENDDO
          WRITE (LU_EVACOUT,'(A)')          '   Nodes out '
          DO NM = 1, EVAC_STRS(N)%N_NODES_OUT
             WRITE (LU_EVACOUT,'(I5,A,A)')          NM, ' ', Trim(EVAC_NODE_List(EVAC_STRS(N)%NODES_OUT(NM))%ID)
          ENDDO
       ENDDO
    END IF

  END SUBROUTINE INITIALIZE_EVAC_DUMPS

!
  SUBROUTINE INITIALIZE_EVACUATION(NM)
    IMPLICIT NONE
    !
    ! Insert humans into the domain at the start of calculation
    !
    ! Passed variables
    INTEGER, INTENT(IN) :: NM
    !
    ! Local variables
    REAL :: RN_REAL
    REAL(EB) RN, simoDX, simoDY, TNOW
    REAL(EB) VOL1, VOL2, X1, X2, Y1, Y2, Z1, Z2, dist, d_max, G_mean, G_sd, G_high, G_low, &
         x1_old, y1_old, CB_TIME, ANGLE, EVEL
    INTEGER i,j,k,ii,jj,kk,ipc, izero, n_tmp, ie, nom, I_OBST, NM_SEE
    INTEGER i11, i22, group_size
    LOGICAL pp_see_group, is_solid, HPT_ORDERED, CROWBAR_READ_IN
    INTEGER jjj, iii, i44, CB_LINE_TYPE, CB_N_AGENTS, CB_ID_CAMERA
    INTEGER(4) IOS
    REAL(EB) x11, y11, group_x_sum, group_y_sum, group_x_center, group_y_center, dens_fac
    INTEGER :: i_endless_loop, istat, I_DX, J_DY
    REAL(EB), DIMENSION(6) :: y_tmp, x_tmp, r_tmp
    REAL(EB), DIMENSION(4) :: d_xy
    LOGICAL, DIMENSION(4) :: FoundWall_xy
    INTEGER, DIMENSION(:), ALLOCATABLE :: CB_I_AGENT
    REAL(EB), DIMENSION(:), ALLOCATABLE :: CB_H_AGENT
    REAL(EB), DIMENSION(:,:), ALLOCATABLE :: CB_XYZ_AGENT,CB_V_AGENT
    !
    TYPE (MESH_TYPE), POINTER :: M =>NULL()
    TYPE (EVAC_SSTAND_TYPE),POINTER :: ESS=>NULL()
    TYPE (EVACUATION_TYPE), POINTER :: HPT=>NULL()
    TYPE (EVAC_PERS_TYPE),  POINTER :: PCP=>NULL()
    TYPE (EVAC_HOLE_TYPE),  POINTER :: EHX=>NULL()
    TYPE (HUMAN_TYPE), POINTER :: HR=>NULL(), HRE=>NULL()
    !
    !IF ( .NOT.(EVACUATION_ONLY(NM) .AND. EVACUATION_GRID(NM)) ) RETURN
    IF (.NOT.(EVACUATION_ONLY(NM))) RETURN
    IF (EMESH_INDEX(NM)==0) RETURN
    ! Next means that only EVAC_PROCESS is doing something
    IF (MYID /= PROCESS(NM)) RETURN

    TNOW = CURRENT_TIME()
    !
    ! Gaussian random numbers, initialize (only once during
    ! the whole calculation is needed). We are now in the
    ! initialization routine, so this is called before the
    ! use of GaussFlag in the evacuation movement part.
    GTrunFlag = 0
    GaussFlag = 0
    !
    MESHES(NM)%N_HUMANS = 0
    MESHES(NM)%N_HUMANS_DIM = MAX_HUMANS_DIM
    ALLOCATE(MESHES(NM)%HUMAN(MESHES(NM)%N_HUMANS_DIM),STAT=IZERO)
    CALL ChkMemErr('INIT_EVACUATION','HUMAN',IZERO)
    !
    ! HUMAN_GRID: (x,y,z): Center of the grid cells
    !             SOOT_DENS: Smoke density at the center of the grid cells (mg/m3)
    !             FED_CO_CO2_O2: FED index
    !             IMESH: Which fire mesh, if any
    !             II,JJ,KK: Fire mesh indexes
    ALLOCATE(MESHES(NM)%HUMAN_GRID(MESHES(NM)%IBAR,MESHES(NM)%JBAR), STAT=IZERO)
    CALL ChkMemErr('INIT_EVACUATION','HUMAN_GRID',IZERO)
    !
    CALL POINT_TO_MESH(NM)
    NM_SEE = NMESHES - N_EGRIDS + EMESH_INDEX(NM)
    !
    ! Initialise HUMAN_GRID
    !
    FED_I_LOOP: DO i = 1,IBAR
       FED_J_LOOP: DO j= 1,JBAR
          x1 = XS + (i-1)*DXI + 0.5_EB*DXI
          y1 = YS + (j-1)*DETA + 0.5_EB*DETA
          ! z1 is here the head height where smoke information is taken, so it should be about the head height
          z1 = 0.5_EB*(ZF+ZS) - EVACUATION_Z_OFFSET(NM) + HUMAN_SMOKE_HEIGHT
          SS_Loop: DO k = 1, N_SSTANDS
             ! Inclines, i.e., EVSS, one should take the smoke at the correct height.
             ESS => EVAC_SSTANDS(k)
             IF (ESS%IMESH == nm .AND. (ESS%X1 <= x1 .AND. ESS%X2 >= x1) .AND. &
                  (ESS%Y1 <= y1 .AND. ESS%Y2 >= y1) ) THEN
                SELECT CASE (ESS%IOR)
                CASE(-1)
                   z1 = z1 + ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%X1-x1)/ABS(ESS%X1-ESS%X2)
                CASE(+1)
                   z1 = z1 + ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%X2-x1)/ABS(ESS%X1-ESS%X2)
                CASE(-2)
                   z1 = z1 + ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%Y1-y1)/ABS(ESS%Y1-ESS%Y2)
                CASE(+2)
                   z1 = z1 + ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%Y2-y1)/ABS(ESS%Y1-ESS%Y2)
                END SELECT
                EXIT SS_Loop
             END IF
          END DO SS_Loop
          HUMAN_GRID(i,j)%X = x1
          HUMAN_GRID(i,j)%Y = y1
          HUMAN_GRID(i,j)%Z = z1
          HUMAN_GRID(i,j)%SOOT_DENS     = 0.0_EB
          HUMAN_GRID(i,j)%FED_CO_CO2_O2 = 0.0_EB
          HUMAN_GRID(i,j)%TMP_G         = 0.0_EB
          HUMAN_GRID(i,j)%RADFLUX        = 0.0_EB
          HUMAN_GRID(i,j)%IMESH         = 0
          HUMAN_GRID(i,j)%II = i
          HUMAN_GRID(i,j)%JJ = j
          HUMAN_GRID(i,j)%KK = 1

          ! If there are no fire grids, skip this intialization part
          IF (.NOT. BTEST(I_EVAC,4) ) CYCLE FED_J_LOOP

          MESH_LOOP: DO NOM = 1, NMESHES
             IF (.NOT. EVACUATION_ONLY(NOM)) THEN
                M => MESHES(NOM)
                IF ( X1 >= M%XS .AND. X1 < M%XF .AND. &
                     Y1 >= M%YS .AND. Y1 < M%YF .AND. &
                     Z1 >= M%ZS .AND. Z1 < M%ZF) THEN
                   II = FLOOR( M%CELLSI(FLOOR((X1-M%XS)*M%RDXINT)) + 1.0_EB  )
                   JJ = FLOOR( M%CELLSJ(FLOOR((Y1-M%YS)*M%RDYINT)) + 1.0_EB  )
                   KK = FLOOR( M%CELLSK(FLOOR((Z1-M%ZS)*M%RDZINT)) + 1.0_EB  )
                   I_OBST = M%OBST_INDEX_C(M%CELL_INDEX(II,JJ,KK))
                   IF ( M%SOLID(M%CELL_INDEX(II,JJ,KK)) .AND. .NOT.M%OBSTRUCTION(I_OBST)%HIDDEN ) THEN
                      HUMAN_GRID(i,j)%IMESH = 0 ! No smoke inside OBSTs
                   ELSE
                      HUMAN_GRID(i,j)%II = II
                      HUMAN_GRID(i,j)%JJ = JJ
                      HUMAN_GRID(i,j)%KK = KK
                      HUMAN_GRID(i,j)%IMESH = NOM
                   END IF
                   EXIT MESH_LOOP
                END IF
             END IF
             ! No fire mesh is found
             HUMAN_GRID(i,j)%IMESH = 0
          END DO MESH_LOOP
          I_OBST = OBST_INDEX_C(CELL_INDEX(I,J,1))
          IF (.NOT. (SOLID(CELL_INDEX(i,j,1)) .AND. .NOT. OBSTRUCTION(I_OBST)%HIDDEN)) THEN
             HUMAN_GRID(i,j)%IMESH = HUMAN_GRID(i,j)%IMESH
          ELSE
             ! This grid cell is solid ==> No humans in this cell
             ! Zero, if fire mesh obst or no fire mesh at all
             ! Negative (-nom), if fire mesh gas cell but evac mesh solid.
             HUMAN_GRID(i,j)%IMESH = -HUMAN_GRID(i,j)%IMESH
          END IF
       END DO FED_J_LOOP
    END DO FED_I_LOOP

    EVAC_CLASS_LOOP: DO IPC=1,NPC_EVAC
       !
       HPT=>EVAC_EVACS(IPC)
       !
       ! Check the mesh
       IF (HPT%IMESH /= NM) CYCLE EVAC_CLASS_LOOP
       !
       ! If there is an initial number of humans, initialize
       IF (HPT%N_INITIAL == 0) CYCLE EVAC_CLASS_LOOP
       !
       IF (HPT%X1 == 0.0_EB .AND. HPT%X2 == 0.0_EB .AND. &
            HPT%Y1 == 0.0_EB .AND. HPT%Y2 == 0.0_EB .AND. &
            HPT%Z1 == 0.0_EB .AND. HPT%Z2 == 0.0_EB ) THEN
          X1 = XS ; X2 = XF
          Y1 = YS ; Y2 = YF
          Z1 = ZS ; Z2 = ZF
          VOL2 = (XF - XS) * (YF - YS)
          VOL1 = VOL2
       ELSE
          IF (HPT%X1 > XF .OR. HPT%X2 < XS .OR. &
               HPT%Y1 > YF .OR. HPT%Y2 < YS .OR. &
               HPT%Z1 > ZF .OR. HPT%Z2 < ZS) THEN
             CYCLE EVAC_CLASS_LOOP
          END IF
          X1 = MAX(HPT%X1,XS) ; X2 = MIN(HPT%X2, XF)
          Y1 = MAX(HPT%Y1,YS) ; Y2 = MIN(HPT%Y2, YF)
          Z1 = MAX(HPT%Z1,ZS) ; Z2 = MIN(HPT%Z2, ZF)
          VOL2 = (HPT%X2 - HPT%X1) * (HPT%Y2 - HPT%Y1)
          VOL1 = (X2 - X1) * (Y2 - Y1) * (Z2 - Z1)
       END IF

       !CROWBAR project: read in positions and velocities of the agents
       IF (TRIM(HPT%CROWBAR_INPUT_FILE)/=TRIM('null')) THEN
          CROWBAR_READ_IN = .TRUE. ! Read the initial position and vels from a file
       ELSE
          CROWBAR_READ_IN = .FALSE.
       END IF
       IF (CROWBAR_READ_IN) THEN
          HPT%GN_MIN = 1; HPT%GN_MAX = 1 !
          IOS=0
          CB_TIME=0.0_EB
          CB_ID_CAMERA=0
          OPEN (LU_EVAC_CB,file=TRIM(HPT%CROWBAR_INPUT_FILE),form='formatted', status='old')
          READ (LU_EVAC_CB,*) CB_LINE_TYPE, CB_TIME, CB_N_AGENTS, CB_ID_CAMERA
          IF (CB_LINE_TYPE/=0 .OR. CB_N_AGENTS<0) THEN
             WRITE(MESSAGE,'(A)') 'ERROR: Initialize Evacuation:: CROWBAR INPUT FILE READ ERROR 1'
             CLOSE (LU_EVAC_CB)
             CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
          END IF
          WRITE(LU_EVACOUT,'(A,A)') 'EVAC: Crowbar input file is: ',TRIM(HPT%CROWBAR_INPUT_FILE)
          WRITE(LU_EVACOUT,'(I2,F8.2,I6,I16)') CB_LINE_TYPE, CB_TIME, CB_N_AGENTS, CB_ID_CAMERA
          IF (CB_N_AGENTS>HPT%N_INITIAL) THEN
             WRITE(MESSAGE,'(A)') 'ERROR: Initialize Evacuation:: CROWBAR N_AGENTS'
             CLOSE (LU_EVAC_CB)
             CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
          END IF
          HPT%N_INITIAL = CB_N_AGENTS ! Use the number that is going to be read in
          ALLOCATE(CB_I_AGENT(CB_N_AGENTS),STAT=IZERO)
          CALL ChkMemErr('INIT_EVACUATION','CB_I_AGENT',IZERO)
          ALLOCATE(CB_XYZ_AGENT(CB_N_AGENTS,3),STAT=IZERO)
          CALL ChkMemErr('INIT_EVACUATION','CB_XYZ_AGENT',IZERO)
          ALLOCATE(CB_V_AGENT(CB_N_AGENTS,3),STAT=IZERO)
          CALL ChkMemErr('INIT_EVACUATION','CB_V_AGENT',IZERO)
          ALLOCATE(CB_H_AGENT(CB_N_AGENTS),STAT=IZERO)
          CALL ChkMemErr('INIT_EVACUATION','CB_H_AGENT',IZERO)
          CB_I_AGENT = 0
          CB_XYZ_AGENT = 0.0_EB; CB_H_AGENT = 0.0_EB; CB_V_AGENT = 0.0_EB
          DO I = 1, CB_N_AGENTS
             CB_LINE_TYPE = -1
             READ (LU_EVAC_CB,*) CB_LINE_TYPE, CB_I_AGENT(I),&
                  CB_XYZ_AGENT(I,1),CB_XYZ_AGENT(I,2),CB_XYZ_AGENT(I,3),CB_H_AGENT(I),&
                  CB_V_AGENT(I,1),CB_V_AGENT(I,2),CB_V_AGENT(I,3)
             IF (CB_LINE_TYPE/=1) THEN
                WRITE(MESSAGE,'(A)') 'ERROR: Initialize Evacuation:: CROWBAR INPUT FILE READ ERROR 2'
                CLOSE (LU_EVAC_CB)
                CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
             END IF
          END DO
          CB_XYZ_AGENT = CB_XYZ_AGENT/1000.0_EB ! mm to m transformation
          CB_H_AGENT   = CB_H_AGENT  /1000.0_EB ! mm to m transformation
          CB_V_AGENT   = CB_V_AGENT  /1000.0_EB ! mm/s to m/s transformation
          WRITE(LU_EVACOUT,'(A)') 'i_agent x_i     y_i     z_i   h_i     vx_i    vy_i    vz_i'
          DO I = 1, CB_N_AGENTS
             WRITE(LU_EVACOUT,'(I4,3F8.2,F6.2,3F8.2)')CB_I_AGENT(I),&
                  CB_XYZ_AGENT(I,1),CB_XYZ_AGENT(I,2),CB_XYZ_AGENT(I,3),CB_H_AGENT(I),&
                  CB_V_AGENT(I,1),CB_V_AGENT(I,2),CB_V_AGENT(I,3)
          END DO
          CLOSE(LU_EVAC_CB)
       END IF

       IF (HPT%DELTA_X>0.0_EB .AND. HPT%DELTA_Y>0.0_EB) THEN
          HPT_ORDERED=.TRUE. ! Ordered input, delta_x, delta_y, XB-midpoint used
       ELSE
          HPT_ORDERED=.FALSE. ! Random placement
       END IF
       !
       ! Check which evacuation floor node  (=1,...,n_egrids)
       n_tmp = 0
       HP_MeshLoop: DO i = 1, NMESHES
          !IF (EVACUATION_ONLY(I) .AND. EVACUATION_GRID(I)) THEN
          IF (EVACUATION_ONLY(I) .AND. EMESH_INDEX(I)>0) THEN
             n_tmp = n_tmp +1
             IF (HPT%IMESH == i) THEN
                EXIT HP_MeshLoop
             END IF
          END IF
       END DO HP_MeshLoop
       IF (n_tmp < 1 .OR. n_tmp > n_egrids) THEN
          WRITE(MESSAGE,'(A,A,A,I4)') 'ERROR: INIT_EVAC: EVAC ',TRIM(HPT%ID),' problem evac node, INODE= ',n_tmp
          CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
       END IF
       !
       PCP => EVAC_PERSON_CLASSES(HPT%IPC)
       !
       ! i11: counter, humans per EVAC-line
       ! i22: group member index
       ! i33: group index (=0: no group, i.e., single humans)
       ! ilh: lonely human index
       !
       i11 = 0
       I_DX = 0; J_DY = 0

       INITIALIZATION_LOOP: DO I=1,HPT%N_INITIAL
          !
          CALL RANDOM_NUMBER(RN_REAL)
          RN = REAL(RN_REAL,EB)
          group_size = HPT%GN_MIN - 1 +  INT((HPT%GN_MAX-HPT%GN_MIN+1)*RN+0.5_EB)
          group_size = MAX(MIN(group_size, HPT%GN_MAX),HPT%GN_MIN)

          IF ( i11+group_size > HPT%N_INITIAL ) THEN
             group_size = HPT%N_INITIAL - i11
          END IF

          IF (i11 >= HPT%N_INITIAL) THEN
             EXIT INITIALIZATION_LOOP
          END IF
          IF (HPT_ORDERED .AND. I_DX < -0.5_EB*(X2-X1)/HPT%DELTA_X) THEN
             EXIT INITIALIZATION_LOOP
          END IF

          i22 = 0
          i_endless_loop = 0
          group_X_sum = 0
          group_Y_sum = 0
          GROUP_SIZE_LOOP: DO
             ! i22: Counter, humans on this group (group index i33)
             i22 = i22 + 1
             IF (i22 > group_size) THEN
                EXIT GROUP_SIZE_LOOP
             END IF
             i11 = i11 + 1

             IF (i22 == 1) THEN
                ! One member groups are not counted as group
                IF (group_size > 1) i33 = i33 + 1
                IF (group_size == 1) ilh = ilh + 1
             END IF
             N_HUMANS = N_HUMANS + 1
             !
             IF (N_HUMANS > N_HUMANS_DIM) THEN
                CALL SHUTDOWN('ERROR: Init Humans: no re-allocation yet',PROCESS_0_ONLY=.FALSE.) ; RETURN
                CALL RE_ALLOCATE_HUMANS(1,NM)
                HUMAN=>MESHES(NM)%HUMAN
             END IF
             !
             HR => HUMAN(N_HUMANS)
             HR%IPC = HPT%IPC  ! PERS-line index
             HR%IEL = IPC      ! EVAC-line index
             CALL CLASS_PROPERTIES(HR,PCP,HR%IEL) ; IF (STOP_STATUS>NO_STOP) RETURN
             HR%I_Target = 0
             HR%I_DoorAlgo = HPT%I_AGENT_TYPE
             CALL RANDOM_NUMBER(RN_REAL)
             RN = REAL(RN_REAL,EB)
             HR%F_FallDown = F_MIN_FALL + RN*(F_MAX_FALL - F_MIN_FALL)
             HR%T_FallenDown = HPT%TIME_FALL_DOWN
             ! HR%T_FallenDown = HUGE(HR%T_FallenDown)
             HR%Angle_FallenDown = 0.0_EB
             HR%SizeFac_FallenDown = 0.0_EB
             CALL RANDOM_NUMBER(RN_REAL)
             RN = REAL(RN_REAL,EB)
             HR%T_CheckFallDown = T_BEGIN + RN

             !
             BLK_LOOP:  DO
                IF (i22 == 1) THEN
                   ! First member of the group
                   IF (HPT_ORDERED) THEN
                      IF (I_DX < -FLOOR(0.5_EB*(X2-X1)/HPT%DELTA_X)) THEN
                         N_HUMANS = N_HUMANS - 1
                         CYCLE INITIALIZATION_LOOP
                      END IF
                      HR%X = I_DX*HPT%DELTA_X + 0.5*(X1+X2)
                      HR%Y = J_DY*HPT%DELTA_Y + 0.5*(Y1+Y2)
                      x1_old = HR%X ; y1_old = HR%Y
                      J_DY = J_DY + SIGN(1,J_DY) ! if j_dy=0 => sign is positive
                      IF (J_DY > 0.5_EB*(Y2-Y1)/HPT%DELTA_Y) J_DY = -1
                      IF (J_DY < -0.5_EB*(Y2-Y1)/HPT%DELTA_Y) THEN
                         J_DY = 0
                         I_DX = I_DX + SIGN(1,I_DX)
                      END IF
                      IF (I_DX >  FLOOR(0.5_EB*(X2-X1)/HPT%DELTA_X)) I_DX = -1
                   ELSEIF (CROWBAR_READ_IN) THEN
                      HR%X = CB_XYZ_AGENT(i11,1) ! x_i
                      HR%Y = CB_XYZ_AGENT(i11,2) ! y_i
                      ! These are overwritten at the initial phase (time < t_det+t_pre)
                      ! so set these once again (evac-line should have t_pre=t_det=0.0
                      ! for the initial agents that are read in. (Or some other specific time)
                      HR%U_CB = CB_V_AGENT(i11,1) ! vx_i, x component of the velocity
                      HR%V_CB = CB_V_AGENT(i11,2) ! vy_i, x component of the velocity
                      HR%CROWBAR_READ_IN = .TRUE.
                      x1_old = HR%X
                      y1_old = HR%Y
                   ELSE
                      CALL RANDOM_NUMBER(RN_REAL)
                      RN = REAL(RN_REAL,EB)
                      HR%X = X1 + RN*(X2-X1)
                      x1_old = HR%X
                      CALL RANDOM_NUMBER(RN_REAL)
                      RN = REAL(RN_REAL,EB)
                      HR%Y = Y1 + RN*(Y2-Y1)
                      y1_old = HR%Y
                   END IF
                ELSE
                   ! Next members of the group are put around the first member.
                   G_mean = 0.0_EB
                   G_sd   = 4.0_EB  ! units (m) std.dev.
                   G_high = 6.0_EB  ! units (m) cut-off
                   G_low  = -6.0_EB ! units (m) cut-off
                   ! First the angle, then the radial distance
                   CALL RANDOM_NUMBER(RN_REAL)
                   RN = REAL(RN_REAL,EB)
                   simoDX = SIN(2.0_EB*Pi*rn)
                   simoDY = COS(2.0_EB*Pi*rn)
                   G_mean = (2.0_EB/3.0_EB)*SQRT(group_size/(Pi*GROUP_DENS))
                   G_sd   = G_mean  ! units (m) std.dev.
                   G_high = MAX(3.0_EB,3.0_EB/GROUP_DENS)* G_mean ! units (m) cut-off
                   G_low  = 0.25_EB      ! units (m) cut-off
                   GTrunFlag=0
                   rn = GaussTrun(G_mean,G_sd,G_low,G_high)
                   simoDX = rn*simoDX
                   simoDY = rn*simoDY
                   HR%X = MIN(X2,MAX(X1, x1_old + simoDX))
                   HR%Y = MIN(Y2,MAX(Y1, y1_old + simoDY))
                END IF
                HR%Z = Z1 + 0.5_EB*(Z2-Z1) ! The mid level of the evacuation mesh
                IF (CROWBAR_READ_IN) THEN
                   ! Set the body angle according to the velocity that is read in.
                   ! If velocity close to zero => use random angle
                   EVEL = SQRT(HR%U_CB**2 + HR%V_CB**2)
                   IF (EVEL > 0.1_EB) THEN
                      IF (HR%V_CB >= 0.0_EB) THEN
                         ANGLE = ACOS(HR%U_CB/EVEL)
                      ELSE
                         ANGLE = 2.0_EB*PI - ACOS(HR%U_CB/EVEL)
                      END IF
                      HR%Angle = ANGLE
                   ELSE
                      CALL RANDOM_NUMBER(RN_REAL)
                      RN = REAL(RN_REAL,EB)
                      HR%Angle = 2.0_EB*Pi*rn
                   END IF
                ELSEIF (HPT%Angle > -999.9_EB) THEN
                   HR%Angle = HPT%Angle
                ELSE
                   CALL RANDOM_NUMBER(RN_REAL)
                   RN = REAL(RN_REAL,EB)
                   HR%Angle = 2.0_EB*Pi*rn
                END IF

                DO WHILE (HR%Angle >= 2.0_EB*Pi)
                   HR%Angle = HR%Angle - 2.0_EB*Pi
                END DO
                DO WHILE (HR%Angle < 0.0_EB)
                   HR%Angle = HR%Angle + 2.0_EB*Pi
                END DO

                !Check, that a person is not put on top of some other person
                IF (HPT_ORDERED) THEN
                   ! High density is wanted
                   d_max = MIN(0.0_EB, MAX_INITIAL_OVERLAP)
                ELSE
                   d_max = MAX_INITIAL_OVERLAP
                END IF
                dens_fac = MAX(1.0_EB,DENS_INIT)

                VOL2 = MAX(VOL2,0.08_EB)
                IF (i_endless_loop >= 10*INT(dens_fac*(16.0_EB*MAX(1.0_EB,LOG10(2.5_EB*VOL2))) / &
                     MAX(1.0_EB,LOG10((2.5_EB*VOL2)/MAX(1.0_EB,2.5_EB*VOL2-1.0_EB)))) ) THEN
                   WRITE (LU_EVACOUT,fmt='(A,A,A,I4,A,I6)') ' ERROR: Initialize_Humans, EVAC line ', &
                        TRIM(EVAC_EVACS(IPC)%ID), ', Mesh ', NM, ', i_human ', n_humans
                   WRITE (LU_EVACOUT,fmt='(a)') '      x       y       z     Rd      Rt      Rs      ds  '
                   WRITE (LU_EVACOUT,fmt='(3f8.2,4f8.4)') HR%X, HR%Y, HR%Z, &
                        2.0_EB*HR%Radius, HR%r_torso, HR%r_shoulder, HR%d_shoulder
                   STOP_STATUS = SETUP_STOP  ! Stop code: FDS improperly set-up
                   HR%SHOW = .TRUE.
                   HR%COLOR_INDEX = EVAC_AVATAR_NCOLOR  ! Cyan
                   EXIT INITIALIZATION_LOOP
                END IF

                ! EVHO (evacuation hole) checking
                EH_Loop: DO ie = 1, n_holes
                   EHX => EVAC_HOLES(ie)
                   IF (EHX%IMESH /= NM) CYCLE EH_Loop
                   IF ( TRIM(EHX%EVAC_ID) /= 'null' .AND. TRIM(EHX%EVAC_ID) /= TRIM(HPT%ID)) CYCLE EH_Loop
                   IF ( TRIM(EHX%PERS_ID) /= 'null' .AND. TRIM(EHX%PERS_ID) /= TRIM(HPT%CLASS_NAME)) CYCLE EH_Loop
                   IF ( (EHX%Y1 <= HR%Y .AND. EHX%Y2 >= HR%Y) .AND. (EHX%X1 <= HR%X .AND. EHX%X2 >= HR%X) ) THEN
                      ! User should not give too large EVHO:s
                      i_endless_loop = i_endless_loop + 1
                      CYCLE BLK_LOOP
                   END IF
                END DO EH_Loop

                ! Put a human in a mesh, if there is enough empty space, i.e., check the OBSTs.
                Is_Solid = .FALSE.
                KK = 1

                r_tmp(1) = HR%r_shoulder ! right circle
                r_tmp(2) = HR%r_torso    ! center circle
                r_tmp(3) = HR%r_shoulder ! left circle
                y_tmp(1) = HR%Y - COS(HR%angle)*HR%d_shoulder ! right
                x_tmp(1) = HR%X + SIN(HR%angle)*HR%d_shoulder
                y_tmp(2) = HR%Y ! torso
                x_tmp(2) = HR%X
                y_tmp(3) = HR%Y + COS(HR%angle)*HR%d_shoulder ! left
                x_tmp(3) = HR%X - SIN(HR%angle)*HR%d_shoulder

                ! Issue 944 bug fix
                IF (Min(x_tmp(1),x_tmp(2),x_tmp(3)) <= XS) Is_Solid = .TRUE.
                IF (Min(y_tmp(1),y_tmp(2),y_tmp(3)) <= YS) Is_Solid = .TRUE.
                IF (Max(x_tmp(1),x_tmp(2),x_tmp(3)) >= XF) Is_Solid = .TRUE.
                IF (Max(y_tmp(1),y_tmp(2),y_tmp(3)) >= YF) Is_Solid = .TRUE.

                IF (.NOT.Is_Solid) THEN
                   II = FLOOR( CELLSI(FLOOR((x_tmp(1)-XS)*RDXINT)) + 1.0_EB )
                   JJ = FLOOR( CELLSJ(FLOOR((y_tmp(1)-YS)*RDYINT)) + 1.0_EB )
                   I_OBST = OBST_INDEX_C(CELL_INDEX(II,JJ,KK))
                   Is_Solid = (Is_Solid .OR. (SOLID(CELL_INDEX(II,JJ,KK)) .AND. .NOT. OBSTRUCTION(I_OBST)%HIDDEN))
                   II = FLOOR( CELLSI(FLOOR((x_tmp(3)-XS)*RDXINT)) + 1.0_EB )
                   JJ = FLOOR( CELLSJ(FLOOR((y_tmp(3)-YS)*RDYINT)) + 1.0_EB )
                   I_OBST = OBST_INDEX_C(CELL_INDEX(II,JJ,KK))
                   Is_Solid = (Is_Solid .OR. (SOLID(CELL_INDEX(II,JJ,KK)) .AND. .NOT. OBSTRUCTION(I_OBST)%HIDDEN))
                   II = FLOOR( CELLSI(FLOOR((x_tmp(2)-XS)*RDXINT)) + 1.0_EB )
                   JJ = FLOOR( CELLSJ(FLOOR((y_tmp(2)-YS)*RDYINT)) + 1.0_EB )
                   I_OBST = OBST_INDEX_C(CELL_INDEX(II,JJ,KK))
                   Is_Solid = (Is_Solid .OR. (SOLID(CELL_INDEX(II,JJ,KK)) .AND. .NOT. OBSTRUCTION(I_OBST)%HIDDEN))
                END IF

                IF (.NOT.Is_Solid) THEN
                   ! Check the distances to walls (not to wall corners)
                   ! Skip wall force ior = 0: Check all walls (do not put too close to doors/exits)
                   CALL Find_walls(nm, x_tmp(1), y_tmp(1), r_tmp(1), d_max, 0, d_xy, FoundWall_xy, istat)
                   IF (istat/=0) Is_Solid = .TRUE.
                   CALL Find_walls(nm, x_tmp(2), y_tmp(2), r_tmp(2), d_max, 0, d_xy, FoundWall_xy, istat)
                   IF (istat/=0) Is_Solid = .TRUE.
                   CALL Find_walls(nm, x_tmp(3), y_tmp(3), r_tmp(3), d_max, 0, d_xy, FoundWall_xy, istat)
                   IF (istat/=0) Is_Solid = .TRUE.
                   IF (STOP_STATUS==SETUP_STOP) RETURN
                END IF

                IF (.NOT.Is_Solid) THEN
                   ! Check that the agent is not too close to other agents, who already introduced in the calculation.
                   P2PLoop: DO ie = 1, n_humans - 1
                      HRE => HUMAN(ie)
                      r_tmp(4) = HRE%r_shoulder ! right circle
                      r_tmp(5) = HRE%r_torso    ! center circle
                      r_tmp(6) = HRE%r_shoulder ! left circle
                      y_tmp(4) = HRE%Y - COS(HRE%angle)*HRE%d_shoulder ! right circle
                      x_tmp(4) = HRE%X + SIN(HRE%angle)*HRE%d_shoulder
                      y_tmp(5) = HRE%Y ! center circle
                      x_tmp(5) = HRE%X
                      y_tmp(6) = HRE%Y + COS(HRE%angle)*HRE%d_shoulder ! left circle
                      x_tmp(6) = HRE%X - SIN(HRE%angle)*HRE%d_shoulder
                      DO iii = 1, 3
                         DO jjj = 4, 6
                            DIST = SQRT((x_tmp(jjj)-x_tmp(iii))**2 + (y_tmp(jjj)-y_tmp(iii))**2) - (r_tmp(jjj)+r_tmp(iii))
                            IF ( DIST < d_max) THEN
                               i_endless_loop = i_endless_loop + 1
                               CYCLE BLK_LOOP
                            END IF
                         END DO
                      END DO
                   END DO P2PLoop

                   IF (i22 > 1) THEN
                      ! Check if the new member will see the first member of the group.
                      X11 = HR%X
                      Y11 = HR%Y

                      PP_see_group = See_each_other(nm_see, x11, y11, x1_old, y1_old)

                   ELSE
                      PP_see_group = .TRUE.
                   END IF

                   IF ( .NOT. PP_see_group ) THEN
                      i_endless_loop = i_endless_loop + 1
                      CYCLE BLK_LOOP
                   END IF

                   ! Coordinates are OK, exit the loop
                   ! Only the centres of the circles are tested agains solids.
                   EXIT BLK_LOOP
                ELSE
                   i_endless_loop = i_endless_loop + 1
                END IF            ! Not inside a solid object
             END DO BLK_LOOP
             !
             ! The coordinates of the new person has passed the tests.
             !
             ILABEL_last = ILABEL_last + 1
             HR%ILABEL = ILABEL_last
             IF (group_size > 1 ) THEN
                ! Is a member of a group
                HR%GROUP_ID = i33
             ELSE
                ! Is a lonely soul
                HR%GROUP_ID = -ilh
             END IF
             HR%Commitment = 0.0_EB
             HR%SHOW = .TRUE.


             SELECT CASE (COLOR_METHOD)
             CASE (-1)
                HR%COLOR_INDEX = 1
             CASE (0)
                HR%COLOR_INDEX = HPT%Avatar_Color_Index
             CASE (1)
                HR%COLOR_INDEX = MOD(group_size-1,6) + 1
             CASE (2)
                IF (HR%GROUP_ID > 0 ) THEN
                   HR%COLOR_INDEX = MOD(HR%GROUP_ID,6) + 1
                ELSE
                   HR%COLOR_INDEX = 1 ! lonely human
                END IF
             CASE (3)
                HR%COLOR_INDEX = evac_person_classes(HPT%IPC)%Avatar_Color_Index
             CASE (4)
                HR%COLOR_INDEX = 1
             CASE (5)
                HR%COLOR_INDEX = 1
             CASE Default
                WRITE(MESSAGE,'(A,I3,A)') 'ERROR: READ_EVAC COLOR METHOD',COLOR_METHOD, ' is not defined'
                CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
             END SELECT

             HR%IMESH       = HPT%IMESH
             HR%INODE       = n_tmp
             HR%NODE_NAME   = TRIM(EVAC_Node_List(n_tmp)%ID)
             HR%I_FFIELD  = 0
             HR%IOR       = HUMAN_NO_TARGET
             IF (.NOT.CROWBAR_READ_IN) THEN
                HR%U         = 0.0_EB
                HR%V         = 0.0_EB
             END IF
             HR%W         = 0.0_EB
             HR%F_Y       = 0.0_EB
             HR%F_X       = 0.0_EB
             HR%v0_fac    = 1.0_EB
             HR%SumForces = 0.0_EB
             HR%SumForces2 = 0.0_EB
             HR%IntDose   = 0.0_EB
             HR%Eta       = 0.0_EB
             HR%Ksi       = 0.0_EB
             HR%NewRnd    = .TRUE.
             HR%STR_SUB_INDX = 0
             HR%SKIP_WALL_FORCE_IOR = 0

             group_X_sum = group_X_sum + HR%X
             group_Y_sum = group_Y_sum + HR%Y
             IF (group_size == 1) THEN
                EXIT  GROUP_SIZE_LOOP
             END IF
             IF (i22 == group_size .AND. group_size > 1) THEN
                group_X_center = group_X_sum / MAX(1,group_size)
                group_Y_center = group_Y_sum / MAX(1,group_size)
                !
                II = FLOOR( CELLSI(FLOOR((group_X_center-XS)*RDXINT)) + 1.0_EB )
                JJ = FLOOR( CELLSJ(FLOOR((group_Y_center-YS)*RDYINT)) + 1.0_EB )
                KK = 1
                x1_old = group_X_center
                y1_old = group_Y_center
                PP_see_group = .TRUE.
                DO i44 = 1,group_size
                   HR => HUMAN( N_HUMANS - group_size + i44 )
                   X11 = HR%X
                   Y11 = HR%Y
                   PP_see_group = See_each_other(nm_see, x11, y11, x1_old, y1_old)

                   IF ( .NOT. PP_see_group ) THEN
                      i_endless_loop = i_endless_loop + 1
                      i22 = 0  ! Start at the beginning of the group
                      i33 = i33 - 1  ! group index
                      i11 = i11 - group_size ! human index (evac-line)
                      N_HUMANS = N_HUMANS - group_size ! total # of humans
                      ILABEL_last = ILABEL_last - group_size ! labels of humans
                      CYCLE GROUP_SIZE_LOOP
                   END IF
                END DO

                EXIT GROUP_SIZE_LOOP

             END IF ! group_size>1 and i22=group_size

          END DO GROUP_SIZE_LOOP ! i22, group_size loop
          !
       END DO INITIALIZATION_LOOP  ! i, # persons in a evac-line
       !
       IF (CROWBAR_READ_IN) THEN
          DEALLOCATE(CB_H_AGENT)
          DEALLOCATE(CB_V_AGENT)
          DEALLOCATE(CB_XYZ_AGENT)
          DEALLOCATE(CB_I_AGENT)
       END IF
    END DO EVAC_CLASS_LOOP ! ipc, number of evac-lines

    WRITE (LU_EVACOUT,fmt='(a,f8.2,a,i0,a,i0/)') ' EVAC: Time ', 0.0_EB,' mesh ',nm,' number of humans ',n_humans
    T_USED(12)=T_USED(12)+CURRENT_TIME()-TNOW
    !
  END SUBROUTINE INITIALIZE_EVACUATION

  !
  SUBROUTINE INIT_EVAC_GROUPS
    IMPLICIT NONE
    !
    ! Initialize group lists, known doors, etc
    !
    ! Passed variables
    !
    ! Local variables
    INTEGER I,J, IZERO, nom, j1, ii, i_target_old, i_change_old, i_tmp, i_tmp2, I_MODE
    INTEGER :: i_egrid, i_target, color_index, i_new_ffield
    REAL(EB) :: evel, TNOW
    !
    TYPE (MESH_TYPE), POINTER :: M=>NULL()
    TYPE (EVAC_SSTAND_TYPE), POINTER :: ESS=>NULL()
    TYPE (HUMAN_TYPE), POINTER :: HR=>NULL()
    !
    IF (.NOT.ANY(EVACUATION_ONLY)) RETURN
    IF (STOP_STATUS > 0) RETURN

    TNOW=CURRENT_TIME()

    !
    ilh_dim = ilh           ! lonely humans dimension
    i33_dim = i33
    ALLOCATE(Group_List(0:i33_dim),STAT=IZERO)
    CALL ChkMemErr('Initialize_Evacuation','Group_List',IZERO)
    Group_List(:)%Tdoor   = HUGE(Group_List(0)%Tdoor)
    Group_List(:)%COMPLETE = 0
    !
    DO i = 0,i33_dim
       ALLOCATE(Group_List(i)%GROUP_I_FFIELDS(n_egrids),STAT=izero)
       Group_List(i)%GROUP_I_FFIELDS(:) = 0
    END DO

    ALLOCATE(Group_Known_Doors(1:i33_dim),STAT=IZERO)
    CALL ChkMemErr('Initialize_Evacuation', 'Group_Known_Doors',IZERO)
    ALLOCATE(Human_Known_Doors(1:ilh_dim),STAT=IZERO)
    CALL ChkMemErr('Initialize_Evacuation', 'Human_Known_Doors',IZERO)
    !
    ! These arrays are common for the whole module EVAC.  They are only allocated once.
    ALLOCATE(Is_Known_Door(MAX(1,N_DOORS+N_EXITS)),STAT=IZERO)
    CALL ChkMemErr('Initialize_Evacuation','Is_Known_Door',IZERO)
    ALLOCATE(Is_Visible_Door(MAX(1,N_DOORS+N_EXITS)),STAT=IZERO)
    CALL ChkMemErr('Initialize_Evacuation','Is_Visible_Door',IZERO)
    ALLOCATE(Is_BeeLine_Door(MAX(1,N_DOORS+N_EXITS),3),STAT=IZERO)
    CALL ChkMemErr('Initialize_Evacuation','Is_BeeLine_Door',IZERO)
    ALLOCATE(Is_BeeLineXB_Door(MAX(1,N_DOORS+N_EXITS),3),STAT=IZERO)
    CALL ChkMemErr('Initialize_Evacuation','Is_BeeLineXB_Door',IZERO)
    ALLOCATE(FED_max_Door(MAX(1,N_DOORS+N_EXITS)),STAT=IZERO)
    CALL ChkMemErr('Initialize_Evacuation','FED_max_Door',IZERO)
    ALLOCATE(K_ave_Door(MAX(1,N_DOORS+N_EXITS)),STAT=IZERO)
    CALL ChkMemErr('Initialize_Evacuation','K_ave_Door',IZERO)
    ALLOCATE(Color_Tmp(MAX(1,i33_dim)),STAT=IZERO)
    CALL ChkMemErr('Initialize_Evacuation','Color_Tmp',IZERO)

    ! Counters for Hawk - Dove game
    IF (C_HAWK >= 0.0_EB) THEN
       ALLOCATE(N_HawkDoveCount(4,MAX(1,n_egrids)),STAT=IZERO)
       CALL ChkMemErr('Initialize_Evacuation','N_HawkDoveCount',IZERO)
       N_HawkDoveCount = 0
    END IF

    i_egrid = 0
    IF (N_DOORS >0) EVAC_DOORS(:)%R_NTARGET = 5.0_EB
    IF (N_EXITS >0) EVAC_EXITS(:)%R_NTARGET = 5.0_EB
    DO NOM = 1, NMESHES
       N_CHANGE_DOORS  = 0 ! Count the initialization Nash equilibrium iterations
       N_CHANGE_TRIALS = 0 ! Count the initialization Nash equilibrium iterations
       I_CHANGE_OLD    = -1
       IF ( .NOT.(EVACUATION_ONLY(NOM) .AND. EMESH_INDEX(NOM)>0) ) CYCLE
       M => MESHES(NOM)
       GROUP_LIST(:)%GROUP_SIZE  = 0
       GROUP_LIST(:)%GROUP_X = 0.0_EB
       GROUP_LIST(:)%GROUP_Y = 0.0_EB
       GROUP_LIST(:)%SPEED   = 0.0_EB
       I_EGRID = I_EGRID+1
       DO I = 1, M%N_HUMANS
          HR => M%HUMAN(I)
          J = MAX(0,HR%GROUP_ID)
          J1 = -MIN(0,HR%GROUP_ID)
          IF (j > 0) THEN    ! The agent is a member of a group
             Group_Known_Doors(j)%N_nodes = -1
          ELSE    ! The agent is a lonely one
             Human_Known_Doors(j1)%N_nodes = -1
          END IF
          GROUP_LIST(J)%GROUP_SIZE = GROUP_LIST(J)%GROUP_SIZE + 1
          GROUP_LIST(J)%GROUP_X    = GROUP_LIST(J)%GROUP_X + HR%X
          GROUP_LIST(J)%GROUP_Y    = GROUP_LIST(J)%GROUP_Y + HR%Y
          GROUP_LIST(J)%SPEED      = GROUP_LIST(J)%SPEED + HR%SPEED
          DO II = 1,N_DOORS
             IF (NOM /= EVAC_DOORS(II)%IMESH) CYCLE
             EVEL = SQRT((HR%X-EVAC_DOORS(II)%X)**2 + (HR%Y-EVAC_DOORS(II)%Y)**2)
             EVAC_DOORS(II)%R_NTARGET = MAX(EVEL,EVAC_DOORS(II)%R_NTARGET)
          END DO
          DO II = 1,N_EXITS
             IF (NOM /= EVAC_EXITS(II)%IMESH) CYCLE
             EVEL = SQRT((HR%X-EVAC_EXITS(II)%X)**2 + (HR%Y-EVAC_EXITS(II)%Y)**2)
             EVAC_EXITS(II)%R_NTARGET = MAX(EVEL,EVAC_EXITS(II)%R_NTARGET)
          END DO
       END DO
       GROUP_LIST(1:)%GROUP_X = GROUP_LIST(1:)%GROUP_X / MAX(1,GROUP_LIST(1:)%GROUP_SIZE)
       GROUP_LIST(1:)%GROUP_Y = GROUP_LIST(1:)%GROUP_Y / MAX(1,GROUP_LIST(1:)%GROUP_SIZE)
       GROUP_LIST(1:)%SPEED   = GROUP_LIST(1:)%SPEED / MAX(1,GROUP_LIST(1:)%GROUP_SIZE)
       I_TMP = 0
       I_TMP2 = -1
       IF (M%N_HUMANS < 1) CYCLE
       I_MODE = 0
       DO WHILE (.NOT.(I_CHANGE_OLD == N_CHANGE_DOORS) .AND. .NOT.(I_TMP2 == I_TMP)) ! Iterate Nash equilibrium
          I_CHANGE_OLD = N_CHANGE_DOORS
          I_TMP = 0
          DO I = 1, M%N_HUMANS
             HR => M%HUMAN(I)
             HR%I_Door_Mode = 0 ! Default, no target door yet
             ! GROUP_ID > 0: +GROUP_ID
             ! GROUP_ID < 0: -HUMAN_ID (lonely humans)
             J  =  MAX(0,HR%GROUP_ID)
             J1 = -MIN(0,HR%GROUP_ID)
             ! Test, if this group has already a ffield (on this floor)
             ! Lonely humans have J=0 and GROUP_I_FFIELDS is always 0.
             IF (GROUP_LIST(J)%GROUP_I_FFIELDS(I_EGRID) == 0) THEN
                I_TARGET_OLD = HR%I_TARGET
                N_CHANGE_TRIALS = N_CHANGE_TRIALS + 1
                CALL CHANGE_TARGET_DOOR(NOM, NOM, I, J, J1, I_EGRID, I_MODE, HR%X, HR%Y, I_TARGET, COLOR_INDEX, I_NEW_FFIELD, HR)
                IF (STOP_STATUS>NO_STOP) RETURN
                IF (ABS(I_TARGET_OLD) /= ABS(I_TARGET)) THEN
                   N_CHANGE_DOORS = N_CHANGE_DOORS + 1
                   I_TMP = I
                   IF (I_TARGET > 0) THEN
                      IF (I_TARGET <= N_DOORS ) THEN
                         EVEL = SQRT((HR%X-0.5_EB*(EVAC_DOORS(I_TARGET)%X1+EVAC_DOORS(I_TARGET)%X2))**2 + &
                              (HR%Y-0.5_EB*(EVAC_DOORS(I_TARGET)%Y1+EVAC_DOORS(I_TARGET)%Y2))**2)
                         EVEL = 50.0_EB*EVEL/EVAC_DOORS(I_TARGET)%R_NTARGET + 1.0_EB
                         II = MIN(50,MAX(1,INT(EVEL)))
                         EVAC_DOORS(I_TARGET)%NTARGET(II:50) = EVAC_DOORS(I_TARGET)%NTARGET(II:50) + 1
                      ELSE
                         EVEL = SQRT((HR%X-0.5_EB*(EVAC_EXITS(I_TARGET-N_DOORS)%X1+EVAC_EXITS(I_TARGET-N_DOORS)%X2))**2 + &
                              (HR%Y-0.5_EB*(EVAC_EXITS(I_TARGET-N_DOORS)%Y1+EVAC_EXITS(I_TARGET-N_DOORS)%Y2))**2)
                         EVEL = 50.0_EB*EVEL/EVAC_EXITS(I_TARGET-N_DOORS)%R_NTARGET + 1.0_EB
                         II = MIN(50,MAX(1,INT(EVEL)))
                         EVAC_EXITS(I_TARGET-N_DOORS)%NTARGET(II:50) = EVAC_EXITS(I_TARGET-N_DOORS)%NTARGET(II:50) + 1
                      END IF
                   ELSE  ! I_TARGET < 0 means non-visible target door
                      N_CHANGE_DOORS = N_CHANGE_DOORS - 1  ! Do not iterate non-visible doors
                   END IF
                   IF (I_TARGET_OLD > 0) THEN
                      IF (I_TARGET_OLD <= N_DOORS ) THEN
                         EVEL = SQRT((HR%X-0.5_EB*(EVAC_DOORS(I_TARGET_OLD)%X1+EVAC_DOORS(I_TARGET_OLD)%X2))**2 + &
                              (HR%Y-0.5_EB*(EVAC_DOORS(I_TARGET_OLD)%Y1+EVAC_DOORS(I_TARGET_OLD)%Y2))**2)
                         EVEL = 50.0_EB*EVEL/EVAC_DOORS(I_TARGET_OLD)%R_NTARGET + 1.0_EB
                         II = MIN(50,MAX(1,INT(EVEL)))
                         EVAC_DOORS(I_TARGET_OLD)%NTARGET(II:50) = EVAC_DOORS(I_TARGET_OLD)%NTARGET(II:50) - 1
                      ELSE
                         EVEL = SQRT((HR%X-0.5_EB*(EVAC_EXITS(I_TARGET_OLD-N_DOORS)%X1+EVAC_EXITS(I_TARGET_OLD-N_DOORS)%X2))**2+&
                              (HR%Y-0.5_EB*(EVAC_EXITS(I_TARGET_OLD-N_DOORS)%Y1+EVAC_EXITS(I_TARGET_OLD-N_DOORS)%Y2))**2)
                         EVEL = 50.0_EB*EVEL/EVAC_EXITS(I_TARGET_OLD-N_DOORS)%R_NTARGET + 1.0_EB
                         II = MIN(50,MAX(1,INT(EVEL)))
                         EVAC_EXITS(I_TARGET_OLD-N_DOORS)%NTARGET(II:50) = EVAC_EXITS(I_TARGET_OLD-N_DOORS)%NTARGET(II:50) - 1
                      END IF
                   END IF
                END IF
             ELSE              ! This group has already a flow field
                ! This group has already tried to change the field
                IF (COLOR_METHOD == 5 .AND. J > 0) HR%COLOR_INDEX = COLOR_TMP(J)
                IF (COLOR_METHOD == 4 .AND. J > 0) HR%COLOR_INDEX = COLOR_TMP(J)
                HR%I_FFIELD    = GROUP_LIST(J)%GROUP_I_FFIELDS(I_EGRID)
                HR%FFIELD_NAME = TRIM(MESH_NAME(HR%I_FFIELD))
                HR%I_TARGET = GROUP_KNOWN_DOORS(J)%I_TARGET
             END IF           ! First member of a group or a lonely soul
          END DO              ! 1, N_HUMANS
          ! Stop Nash iterations if close enough, NASH_CLOSE_ENOUGH
          IF (REAL(N_CHANGE_DOORS-I_CHANGE_OLD,EB)/REAL(MAX(1,M%N_HUMANS),EB) <= NASH_CLOSE_ENOUGH) THEN
             WRITE(LU_EVACOUT,FMT='(A,2F10.4)') ' INIT: Nash Changes per agent close enough: ', &
                  REAL(N_CHANGE_DOORS-I_CHANGE_OLD,EB)/REAL(MAX(1,M%N_HUMANS),EB), NASH_CLOSE_ENOUGH
             I_TMP2 = I_TMP
          ELSE
             I_TMP2 = -1
          END IF
          ! Maximum number of Nash iterations is twenty per agent
          IF (N_CHANGE_TRIALS/MAX(1,M%N_HUMANS) > 20) THEN
             I_TMP2 = I_TMP
             WRITE(LU_EVACOUT,FMT='(A,2F10.4)') ' INIT: Exit Nash iteration after 20 ones: ', &
                  REAL(N_CHANGE_DOORS-I_CHANGE_OLD,EB)/REAL(MAX(1,M%N_HUMANS),EB), NASH_CLOSE_ENOUGH
          END IF
          IF (ABS(FAC_DOOR_QUEUE) < TWO_EPSILON_EB) I_CHANGE_OLD = N_CHANGE_DOORS  ! DO NOT ITERATE THE NASH EQUILIBRIUM
          I_MODE = 1  ! change_target_door initializations have now been done
       END DO         ! Nash iterations
       IF (ABS(FAC_DOOR_QUEUE) >= TWO_EPSILON_EB) WRITE(LU_EVACOUT,FMT='(A,F14.2,A,I8)') &
            ' INIT: Changes per agent ', REAL(N_CHANGE_DOORS,EB)/REAL(M%N_HUMANS,EB), &
            ', Nash iterations', N_CHANGE_TRIALS/M%N_HUMANS
    END DO                  ! 1, NMESHES

    WRITE (LU_EVACOUT,FMT='(/A)') ' EVAC: Initial positions of the agents'
    WRITE (LU_EVACOUT,FMT='(A,A)') ' Agent      X       Y       Z    Tpre    Tdet  ', &
         ' Dia    V0   Tau   I_gr I_ff'

    ! Initialize the GROUP_I_FFIELDS
    I_EGRID = 0
    DO NOM = 1, NMESHES
       IF ( .NOT.(EVACUATION_ONLY(NOM) .AND. EMESH_INDEX(NOM)>0) ) CYCLE
       M => MESHES(NOM)
       I_EGRID = I_EGRID + 1
       DO I = 1, M%N_HUMANS
          HR => M%HUMAN(I)
          IF (C_HAWK >= 0.0_EB) THEN
             ! Hawk - Dove game counters, initial values
             ! (1,emesh) hawks; (2,emesh) doves; (3,emesh) dummy hawks; (4,emesh) dummy doves
             IF (HR%A <= 500.0_EB .AND. HR%I_DoorAlgo==1) N_HawkDoveCount(1,I_EGRID) = N_HawkDoveCount(1,I_EGRID) + 1
             IF (HR%A >  500.0_EB .AND. HR%I_DoorAlgo==1) N_HawkDoveCount(2,I_EGRID) = N_HawkDoveCount(2,I_EGRID) + 1
             IF (HR%A <= 500.0_EB .AND. HR%I_DoorAlgo==2) N_HawkDoveCount(3,I_EGRID) = N_HawkDoveCount(3,I_EGRID) + 1
             IF (HR%A >  500.0_EB .AND. HR%I_DoorAlgo==2) N_HawkDoveCount(4,I_EGRID) = N_HawkDoveCount(4,I_EGRID) + 1
          END IF

          ! Check spectator stands, correct the z-coordinate
          HR%Z = 0.5_EB*(M%ZS+M%ZF)  ! The agent is not on any incline
          SS_LOOP: DO J = 1, N_SSTANDS
             ESS => EVAC_SSTANDS(J)
             IF (ESS%IMESH == NOM .AND. (ESS%X1 <= HR%X .AND. ESS%X2 >= HR%X) .AND. &
                  (ESS%Y1 <= HR%Y .AND. ESS%Y2 >= HR%Y) ) THEN
                SELECT CASE (ESS%IOR)
                CASE(-1)
                   HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%X1-HR%X)/ABS(ESS%X1-ESS%X2)
                CASE(+1)
                   HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%X2-HR%X)/ABS(ESS%X1-ESS%X2)
                CASE(-2)
                   HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%Y1-HR%Y)/ABS(ESS%Y1-ESS%Y2)
                CASE(+2)
                   HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%Y2-HR%Y)/ABS(ESS%Y1-ESS%Y2)
                END SELECT
                EXIT SS_LOOP
             END IF
          END DO SS_LOOP

          J = MAX(0,HR%GROUP_ID)
          GROUP_LIST(J)%IEL = HR%IEL

          WRITE (LU_EVACOUT,FMT='(I6,5F8.2,3F6.2,I6,I4,I4)') HR%ILABEL, &
               HR%X, HR%Y, HR%Z, HR%TPRE, HR%TDET,2.0_EB*HR%RADIUS, &
               HR%SPEED, HR%TAU, HR%GROUP_ID, HR%I_FFIELD, HR%COLOR_INDEX
       END DO
    END DO
    WRITE (LU_EVACOUT,FMT='(/)')

    T_USED(12)=T_USED(12)+CURRENT_TIME()-TNOW
  END SUBROUTINE INIT_EVAC_GROUPS
!
  SUBROUTINE EVAC_MESH_EXCHANGE(T,T_SAVE,I_MODE, ICYC, EXCHANGE_EVACUATION, MODE)
    IMPLICIT NONE
    !
    ! Passed variables
    REAL(EB), INTENT(IN) :: T
    INTEGER, INTENT(IN) :: I_MODE, ICYC, MODE
    LOGICAL, INTENT(OUT) :: EXCHANGE_EVACUATION
    REAL(EB), INTENT(INOUT) :: T_SAVE
    !
    ! Local variables
    INTEGER :: NM, NOM, I, J, K, I1, J1, K1
    INTEGER :: IOS, IZERO, N_DEVC_WRITE
    LOGICAL L_USE_FED, L_FED_READ, L_FED_SAVE, CURRENT_STATE, PRIOR_STATE
    REAL(EB) DT_SAVE, TNOW, T_CHANGE
    INTEGER(4) IBAR_TMP, JBAR_TMP, KBAR_TMP, N_TMP, N_DEVC_WRITE_TMP, N_TMP2
    REAL(FB) TMPOUT1, TMPOUT2, TMPOUT3, TMPOUT4, T_TMP, DT_TMP
    REAL(FB) TMPOUT5, TMPOUT6, TMPOUT7, TMPOUT8
    CHARACTER(LABEL_LENGTH) DEVC_ID
    REAL(EB), ALLOCATABLE, DIMENSION(:) :: ZZ_GET
    TYPE (DEVICE_TYPE), POINTER :: DV =>NULL()
    TYPE (CONTROL_TYPE), POINTER :: CV =>NULL()
    TYPE (EVAC_EDEV_TYPE), POINTER :: EDV=>NULL()

    EXCHANGE_EVACUATION=.FALSE.
    !
    IF (.NOT. ANY(EVACUATION_ONLY)) RETURN
    IF (ICYC < 1) RETURN
    !
    ! I_MODE: 'binary' index:
    ! XXXXX = (0,1)*16 + (0,1)*8 + (0,1)*4 + (0,1)*2 + (0,1)
    ! 0. BIT (XXXX1): Save flow fields (XXXX0) Do not save
    ! 1. BIT (XXX1X): Save soot + FED  (XXX0X) Do not save
    ! 2. BIT (XX1XX): Read flow fields (XX0XX) Do not read
    ! 3. BIT (X1XXX): Read soot + FED  (X0XXX) Do not read
    ! 4. BIT (1XXXX): Fire calculation (0XXXX) No fire calculation
    ! BTEST(INTEGER,BIT) bit=0 least significant bit
    ! IBCLR(INTEGER,BIT) clears a bit
    ! IBSET(INTEGER,BIT) sets a bit
    !
    ! LU_EVACFED: CHID_EVAC.FED, FED AND SOOT, TIME DEPENDENT, BINARY
    ! FILE FORMAT: 1. ROW: N_EGRIDS >=0  (OLD FORMAT)
    !              1A. ROW: N < 0 (NEW FORMAT)
    !              1B. ROW: N_EGRIDS,4,N_CORRS,8 (NEW FORMATS)
    !              1C. ROW: N_DEVC_EVAC (NEW FORMAT: -3)
    !                 2. ROW: T AND DT
    !                    3. ROW: IBAR,JBAR,KBAR, N_QUANTITIES
    !                       4. ROW: ONWARDS DATA
    !                    GOTO 3. (MESHES)
    !                       N. ROW: CORR DATA (8 REAL NUMBERS)
    !                       N+1. ROW: NEXT CORR DATA...
    !                 GOTO 2. (TIME POINTS)
    !
    ! Update interval (seconds) fire ==> evac information

    TNOW = CURRENT_TIME()
    DT_SAVE = 2.0_EB
    IOS = 0
    L_USE_FED  = .FALSE.
    L_FED_READ = BTEST(I_MODE,3)
    L_FED_SAVE = BTEST(I_MODE,1)
    !
    ! Change information: Fire meshes ==> evac meshes
    IF ( T >= T_BEGIN .AND. REAL(T,FB) >= REAL(T_SAVE,FB) ) THEN
       IF (MODE > 0) T_SAVE = T + DT_SAVE
       L_USE_FED = .TRUE.
    END IF

    L_USE_FED = L_USE_FED .AND. (L_FED_READ .OR. L_FED_SAVE)

    IF (MODE < 1 .AND. L_USE_FED) EXCHANGE_EVACUATION = .TRUE.
    IF (MODE < 2) RETURN

    IF (L_FED_SAVE) THEN
       N_DEVC_WRITE = 0
       DEVC_LOOP_0: DO I = 1, N_DEVC
          DV => DEVICE(I)
          IF (.NOT. DV%EVACUATION) CYCLE DEVC_LOOP_0
          IF (TRIM(DV%QUANTITY)=='TIME' .AND. DV%SETPOINT<=T_BEGIN) CYCLE DEVC_LOOP_0
          N_DEVC_WRITE = N_DEVC_WRITE + 1
          IF (DV%CURRENT_STATE .NEQV. DV%PRIOR_STATE) THEN
             EVAC_DEVICES(N_DEVC_WRITE)%T_Change = DV%T_CHANGE
             EVAC_DEVICES(N_DEVC_WRITE)%CURRENT  = DV%CURRENT_STATE
             EVAC_DEVICES(N_DEVC_WRITE)%PRIOR    = DV%PRIOR_STATE
             EVAC_DEVICES(N_DEVC_WRITE)%USE_NOW  = .FALSE.
          END IF
       END DO DEVC_LOOP_0
       CTRL_LOOP_0: DO I = 1, N_CTRL
          CV => CONTROL(I)
          IF (.NOT. CV%EVACUATION) CYCLE CTRL_LOOP_0
          N_DEVC_WRITE = N_DEVC_WRITE + 1
          IF (CONTROL(I)%CURRENT_STATE .NEQV. CONTROL(I)%PRIOR_STATE) THEN
             EVAC_DEVICES(N_DEVC_WRITE)%T_Change = CV%T_CHANGE
             EVAC_DEVICES(N_DEVC_WRITE)%CURRENT  = CV%CURRENT_STATE
             EVAC_DEVICES(N_DEVC_WRITE)%PRIOR    = CV%PRIOR_STATE
             EVAC_DEVICES(N_DEVC_WRITE)%USE_NOW  = .FALSE.
          END IF
       END DO CTRL_LOOP_0
    END IF

    IF (L_USE_FED) THEN  ! Update at this time step (do this at DT_SAVE intervals)
       IF (L_FED_SAVE) THEN
          ALLOCATE(ZZ_GET(1:N_TRACKED_SPECIES),STAT=IZERO)
          CALL CHKMEMERR('EVAC_MESH_EXCHANGE', 'ZZ_GET',IZERO)
          WRITE (LU_EVACFED) REAL(T,FB), REAL(DT_SAVE,FB)
       ELSE
          READ (LU_EVACFED,END=324,IOSTAT=IOS) T_TMP, DT_TMP
          T_SAVE = T_TMP + DT_TMP ! Next time point in the file
       END IF

       ! Next loop interpolates fire mesh (soot+fed) into human_grids and
       ! saves it to the disk, or it reads fed+soot from the disk.
       MESH_LOOP: DO NM=1,NMESHES
          IF ( .NOT.(EVACUATION_ONLY(NM) .AND. EMESH_INDEX(NM)>0) ) CYCLE
          !
          CALL POINT_TO_MESH(NM)
          IF (L_FED_SAVE) THEN
             IBAR_TMP = IBAR
             JBAR_TMP = JBAR
             KBAR_TMP = 1
             N_TMP    = 4  ! New format (version 1.11)
             WRITE (LU_EVACFED) IBAR_TMP, JBAR_TMP, KBAR_TMP, N_TMP
          ELSE
             READ (LU_EVACFED,IOSTAT=IOS) IBAR_TMP, JBAR_TMP, KBAR_TMP, N_TMP
             IF (IOS/=0) THEN
                WRITE(MESSAGE,'(A)') 'ERROR: EVAC_MESH_EXCHANGE: FED read error'
                CLOSE (LU_EVACFED)
                CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
             END IF
             IF (IBAR_TMP /= IBAR .OR. JBAR_TMP /= JBAR .OR. N_TMP < 4 ) THEN
                CLOSE (LU_EVACFED)
                CALL SHUTDOWN('ERROR: Problems to read the FED file',PROCESS_0_ONLY=.FALSE.) ; RETURN
             END IF

          END IF

          DO I = 1, IBAR
             DO J= 1, JBAR

                IF (L_FED_SAVE) THEN
                   IF ( ABS(HUMAN_GRID(I,J)%IMESH) > 0 ) THEN
                      ! IMESH > 0, i.e. fire grid found
                      I1 = HUMAN_GRID(I,J)%II
                      J1 = HUMAN_GRID(I,J)%JJ
                      K1 = HUMAN_GRID(I,J)%KK
                      NOM = ABS(HUMAN_GRID(I,J)%IMESH)
                      CALL GET_FIRE_CONDITIONS(NOM,I1,J1,K1,&
                           HUMAN_GRID(I,J)%FED_CO_CO2_O2,HUMAN_GRID(I,J)%SOOT_DENS,&
                           HUMAN_GRID(I,J)%TMP_G, HUMAN_GRID(I,J)%RADFLUX, ZZ_GET, FED_ACTIVITY)
                   END IF
                   ! Save FED, SOOT, TEMP(C), and RADFLUX
                   WRITE (LU_EVACFED) &
                        REAL(HUMAN_GRID(I,J)%FED_CO_CO2_O2,FB), &
                        REAL(HUMAN_GRID(I,J)%SOOT_DENS,FB), &
                        REAL(HUMAN_GRID(I,J)%TMP_G,FB), &
                        REAL(HUMAN_GRID(I,J)%RADFLUX,FB)
                ELSE ! Read FED from a file
                   ! Read FED, SOOT, TEMP(C), and RADFLUX
                   READ (LU_EVACFED,IOSTAT=IOS) TMPOUT1, TMPOUT2, TMPOUT3, TMPOUT4
                   IF (IOS/=0) THEN
                      WRITE(MESSAGE,'(A)') 'ERROR: EVAC_MESH_EXCHANGE: FED read error'
                      CLOSE (LU_EVACFED)
                      CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
                   END IF
                   HUMAN_GRID(I,J)%FED_CO_CO2_O2 = TMPOUT1
                   HUMAN_GRID(I,J)%SOOT_DENS = TMPOUT2
                   HUMAN_GRID(I,J)%TMP_G = TMPOUT3
                   HUMAN_GRID(I,J)%RADFLUX = TMPOUT4
                END IF   ! calculate and save FED
             END DO     ! J=1,JBAR
          END DO       ! I=1,IBAR
          IF (DISCARD_SMOKE_INFO)  HUMAN_GRID(:,:)%FED_CO_CO2_O2 = 0.0_EB
          IF (DISCARD_SMOKE_INFO)  HUMAN_GRID(:,:)%TMP_G         = 0.0_EB
          IF (DISCARD_SMOKE_INFO)  HUMAN_GRID(:,:)%SOOT_DENS     = 0.0_EB
          IF (DISCARD_SMOKE_INFO)  HUMAN_GRID(:,:)%RADFLUX       = 0.0_EB

       END DO MESH_LOOP

       ! Next loop interpolates fire mesh (soot+fed) into human_grids and
       ! saves it to the disk, or it reads fed+soot from the disk.
       CORR_LOOP: DO I = 1, N_CORRS
          !
          IF (L_FED_SAVE) THEN

             IF ( EVAC_CORRS(I)%FED_MESH > 0 .AND. .NOT.DISCARD_SMOKE_INFO) THEN
                ! Here the fire properties are saved to the arrays.
                I1 = EVAC_CORRS(I)%II(1)
                J1 = EVAC_CORRS(I)%JJ(1)
                K1 = EVAC_CORRS(I)%KK(1)
                NOM = EVAC_CORRS(I)%FED_MESH
                CALL GET_FIRE_CONDITIONS(NOM,I1,J1,K1,&
                     EVAC_CORRS(I)%FED_CO_CO2_O2(1),EVAC_CORRS(I)%SOOT_DENS(1),&
                     EVAC_CORRS(I)%TMP_G(1), EVAC_CORRS(I)%RADFLUX(1), ZZ_GET, FED_ACTIVITY)
             ELSE
                ! No FED_MESH found
                EVAC_CORRS(I)%FED_CO_CO2_O2(1) = 0.0_EB
                EVAC_CORRS(I)%SOOT_DENS(1) = 0.0_EB
                EVAC_CORRS(I)%TMP_G(1) = 0.0_EB
                EVAC_CORRS(I)%RADFLUX(1) = 0.0_EB
             END IF                ! FED_MESH > 0, i.e. fire grid found

             IF ( EVAC_CORRS(I)%FED_MESH2 > 0 .AND. .NOT.DISCARD_SMOKE_INFO) THEN
                I1 = EVAC_CORRS(I)%II(2)
                J1 = EVAC_CORRS(I)%JJ(2)
                K1 = EVAC_CORRS(I)%KK(2)
                NOM = EVAC_CORRS(I)%FED_MESH2
                CALL GET_FIRE_CONDITIONS(NOM,I1,J1,K1,&
                     EVAC_CORRS(I)%FED_CO_CO2_O2(2),EVAC_CORRS(I)%SOOT_DENS(2),&
                     EVAC_CORRS(I)%TMP_G(2), EVAC_CORRS(I)%RADFLUX(2), ZZ_GET, FED_ACTIVITY)
             ELSE
                ! No FED_MESH2 found
                EVAC_CORRS(I)%FED_CO_CO2_O2(2) = 0.0_EB
                EVAC_CORRS(I)%SOOT_DENS(2) = 0.0_EB
                EVAC_CORRS(I)%TMP_G(2) = 0.0_EB
                EVAC_CORRS(I)%RADFLUX(2) = 0.0_EB
             END IF                ! FED_MESH2 > 0, i.e. fire grid found

             ! Save FED, SOOT, TEMP(C), and RADFLUX
             WRITE (LU_EVACFED) &
                  REAL(EVAC_CORRS(I)%FED_CO_CO2_O2(1),FB), &
                  REAL(EVAC_CORRS(I)%SOOT_DENS(1),FB), &
                  REAL(EVAC_CORRS(I)%TMP_G(1),FB), &
                  REAL(EVAC_CORRS(I)%RADFLUX(1),FB), &
                  REAL(EVAC_CORRS(I)%FED_CO_CO2_O2(2),FB), &
                  REAL(EVAC_CORRS(I)%SOOT_DENS(2),FB), &
                  REAL(EVAC_CORRS(I)%TMP_G(2),FB), &
                  REAL(EVAC_CORRS(I)%RADFLUX(2),FB)
          ELSE                    ! Read FED from a file
             ! Read FED, SOOT, TEMP(C), and RADFLUX
             READ (LU_EVACFED,IOSTAT=IOS) TMPOUT1, TMPOUT2, TMPOUT3, TMPOUT4, TMPOUT5, TMPOUT6, TMPOUT7, TMPOUT8
             IF (IOS/=0) THEN
                WRITE(MESSAGE,'(A)') 'ERROR: EVAC_MESH_EXCHANGE: FED read error'
                CLOSE (LU_EVACFED)
                CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
             END IF
             IF (DISCARD_SMOKE_INFO) THEN
                TMPOUT1 = 0.0_EB; TMPOUT2 = 0.0_EB; TMPOUT3 = 0.0_EB; TMPOUT4 = 0.0_EB
                TMPOUT5 = 0.0_EB; TMPOUT6 = 0.0_EB; TMPOUT7 = 0.0_EB; TMPOUT8 = 0.0_EB
             END IF
             EVAC_CORRS(I)%FED_CO_CO2_O2(1) = TMPOUT1
             EVAC_CORRS(I)%SOOT_DENS(1) = TMPOUT2
             EVAC_CORRS(I)%TMP_G(1) = TMPOUT3
             EVAC_CORRS(I)%RADFLUX(1) = TMPOUT4
             EVAC_CORRS(I)%FED_CO_CO2_O2(2) = TMPOUT5
             EVAC_CORRS(I)%SOOT_DENS(2) = TMPOUT6
             EVAC_CORRS(I)%TMP_G(2) = TMPOUT7
             EVAC_CORRS(I)%RADFLUX(2) = TMPOUT8

          END IF                  ! Calculate and save FED
       END DO CORR_LOOP

       EDEV_FED_SAVE: IF (L_FED_SAVE) THEN
          ! Next loop is for evacuation devices (like heat detectors)
          N_DEVC_WRITE_TMP = 0
          DEVC_LOOP: DO I = 1, N_DEVC_EVAC
             IF (EVAC_DEVICES(I)%CURRENT .NEQV. EVAC_DEVICES(I)%PRIOR .OR. ICYC==1) THEN
                N_DEVC_WRITE_TMP = N_DEVC_WRITE_TMP + 1
             END IF
          END DO DEVC_LOOP
          WRITE (LU_EVACFED) N_DEVC_WRITE_TMP
          DEVC_LOOP_2: DO I = 1, N_DEVC_EVAC
             IF (EVAC_DEVICES(I)%CURRENT .NEQV. EVAC_DEVICES(I)%PRIOR .OR. ICYC==1) THEN
                N_DEVC_WRITE_TMP = EVAC_DEVICES(I)%I_Devc
                N_TMP = I
                N_TMP2 = EVAC_DEVICES(I)%I_Type
                WRITE (LU_EVACFED) N_TMP2, EVAC_DEVICES(I)%DEVC_ID
                WRITE (LU_EVACFED) N_TMP, N_DEVC_WRITE_TMP, EVAC_DEVICES(I)%CURRENT, &
                     EVAC_DEVICES(I)%PRIOR, REAL(EVAC_DEVICES(I)%T_Change,FB)
                EVAC_DEVICES(I)%PRIOR = EVAC_DEVICES(I)%CURRENT
                EVAC_DEVICES(I)%USE_NOW = .TRUE.
             END IF
          END DO DEVC_LOOP_2
       ELSE  ! Read evacuation devices info from the FED file
          READ (LU_EVACFED,IOSTAT=IOS) N_DEVC_WRITE_TMP
          IF (IOS/=0 .OR. N_DEVC_WRITE_TMP > N_DEVC_EVAC) THEN
             WRITE(MESSAGE,'(A)') 'ERROR: EVAC_MESH_EXCHANGE: FED read error1 for DEVC'
             CLOSE (LU_EVACFED)
             CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
          END IF
          N_DEVC_WRITE = N_DEVC_WRITE_TMP
          DEVC_LOOP_3: DO I = 1, N_DEVC_WRITE
             READ (LU_EVACFED,IOSTAT=IOS) N_TMP2, DEVC_ID
             IF (IOS/=0) THEN
                WRITE(MESSAGE,'(A)') 'ERROR: EVAC_MESH_EXCHANGE: FED read error2 for DEVC'
                CLOSE (LU_EVACFED)
                CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
             END IF
             READ (LU_EVACFED,IOSTAT=IOS) N_DEVC_WRITE_TMP, N_TMP, CURRENT_STATE, PRIOR_STATE, T_TMP
             IF (IOS /= 0 .OR. N_DEVC_WRITE_TMP > N_DEVC_EVAC) THEN
                WRITE(MESSAGE,'(A)') 'ERROR: EVAC_MESH_EXCHANGE: FED read error3 for DEVC'
                CLOSE (LU_EVACFED)
                CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
             END IF
             j = N_DEVC_WRITE_TMP
             IF (ICYC==1) THEN
                EVAC_DEVICES(I)%DEVC_ID = TRIM(DEVC_ID)
                EVAC_DEVICES(I)%I_Devc = N_TMP
                EVAC_DEVICES(I)%I_Type = N_TMP2
                EVAC_DEVICES(I)%I_Devc_Evac = j
             END IF
             IF (N_TMP /= EVAC_DEVICES(j)%I_Devc .OR. TRIM(DEVC_ID) /= EVAC_DEVICES(j)%DEVC_ID .OR. &
                  N_TMP2 /= EVAC_DEVICES(j)%I_Type) THEN
                WRITE(MESSAGE,'(A)') 'ERROR: EVAC_MESH_EXCHANGE: FED read error4 for DEVC'
                CLOSE (LU_EVACFED)
                CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
             END IF
             T_CHANGE = T_TMP
             EVAC_DEVICES(j)%T_Change = T_CHANGE
             EVAC_DEVICES(j)%CURRENT = CURRENT_STATE
             EVAC_DEVICES(j)%PRIOR = PRIOR_STATE
             EVAC_DEVICES(I)%USE_NOW = .TRUE.
          END DO DEVC_LOOP_3
          IF (ICYC==1) THEN
             ! Check the inputs for the EDEV namelist array EVAC_EDEV(1:N_EDEV)
             EDEV_LOOP: DO I = 1, N_EDEV
                EDV => EVAC_EDEV(I)
                EDV%INPUT_DEVC_INDEX(:) = 0
                INPUTS_LOOP: DO J = 1, EDV%N_INPUTS
                   DEVICES_LOOP: DO K = 1, N_DEVC_EVAC
                      IF (TRIM(EDV%INPUT_ID(J)) == TRIM(EVAC_DEVICES(K)%DEVC_ID)) THEN
                         EDV%INPUT_DEVC_INDEX(J) = K
                         CYCLE INPUTS_LOOP
                      END IF
                   END DO DEVICES_LOOP
                   WRITE(MESSAGE,'(A,I5,A,A,A)')  'ERROR: EDEV ',I,' input ',TRIM(EDV%INPUT_ID(J)),' is not found'
                   CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
                END DO INPUTS_LOOP
             END DO EDEV_LOOP
          END IF
       END IF EDEV_FED_SAVE

       IF (L_FED_SAVE) DEALLOCATE(ZZ_GET)

    END IF                    ! L_USE_FED

324 CONTINUE
    IF (IOS < 0) THEN
       WRITE (LU_EVACOUT,FMT='(A,F12.4,A)') 'FED FILE EOF: Time ', T_SAVE-DT_SAVE, ' not found'
       WRITE (LU_EVACOUT,FMT='(A)') 'FED file EOF: Use previous values'
       T_SAVE = 1.0E15
    END IF

    T_USED(12) = T_USED(12) + CURRENT_TIME() - TNOW
  END SUBROUTINE EVAC_MESH_EXCHANGE
!
  SUBROUTINE PREPARE_TO_EVACUATE(ICYC)
    IMPLICIT NONE
    !
    ! Do the mesh independent initializations for the subroutine EVACUATE_HUMANS.
    !
    ! Passed variables
    INTEGER, INTENT(IN) :: ICYC
    !
    ! Local variables
    LOGICAL L_EFF_READ, L_EFF_SAVE
    INTEGER(4) IBAR_TMP, JBAR_TMP, KBAR_TMP
    INTEGER NM_TIM, I, J, IFIELD, NFIELDS, N_END, JJ, JJ_NOW
    REAL(FB) U_TMP, V_TMP
    !
    TYPE (MESH_TYPE), POINTER :: MFF=>NULL()

    IF (.NOT.ANY(EVACUATION_ONLY)) RETURN

    L_EFF_READ = BTEST(I_EVAC,2)
    L_EFF_SAVE = BTEST(I_EVAC,0)
    IF (ICYC == 1 .AND. L_EFF_SAVE) THEN
       N_END = N_EXITS - N_CO_EXITS + N_DOORS
       WRITE_EFF_LOOP: DO NM_TIM = 1, NMESHES
          EVAC_ONLY_NM: IF (EVACUATION_ONLY(NM_TIM)) THEN
             IF (.NOT. EMESH_INDEX(NM_TIM)>0) CYCLE WRITE_EFF_LOOP
             MFF=>MESHES(NM_TIM)
             IBAR_TMP = MFF%IBAR
             JBAR_TMP = MFF%JBAR
             KBAR_TMP = 1
             NFIELDS = EMESH_NFIELDS(EMESH_INDEX(NM_TIM))
             DO IFIELD = 1, NFIELDS

                LOOP_EXITS: DO JJ = 1, N_END
                   IF (EMESH_EXITS(JJ)%MAINMESH == NM_TIM .AND. EMESH_EXITS(JJ)%I_DOORS_EMESH == IFIELD) THEN
                      JJ_NOW = JJ
                      EXIT LOOP_EXITS
                   END IF
                END DO LOOP_EXITS

                WRITE (LU_EVACEFF) IBAR_TMP, JBAR_TMP, KBAR_TMP
                DO  I = 0, MFF%IBAR+1
                   DO J= 0, MFF%JBAR+1
                      U_TMP = REAL(EMESH_EXITS(JJ_NOW)%U_EVAC(I,J),FB)
                      V_TMP = REAL(EMESH_EXITS(JJ_NOW)%V_EVAC(I,J),FB)
                      WRITE (LU_EVACEFF) U_TMP, V_TMP
                   END DO
                END DO
             END DO
          END IF EVAC_ONLY_NM
       END DO WRITE_EFF_LOOP
       ! Clear the save bit, save is done only once.
       L_EFF_SAVE = .FALSE.
       I_EVAC = IBCLR(I_EVAC,0)
    END IF
    !
    ! Initialize counters only once for each time step.
    IF ( ICYC >= 0 .AND. ICYC_OLD < ICYC ) THEN
       ICYC_OLD = ICYC
       FED_MAX_ALIVE = 0.0_EB
       FED_MAX       = 0.0_EB
    END IF

  END SUBROUTINE PREPARE_TO_EVACUATE
!
  SUBROUTINE CLEAN_AFTER_EVACUATE(ICYC,I_MODE)
    IMPLICIT NONE
    !
    ! Do the mesh independent clean up for the subroutine EVACUATE_HUMANS.
    !
    ! Passed variables
    INTEGER, INTENT(IN) :: ICYC, I_MODE
    !
    ! Local variables
    !
    INTRINSIC :: BTEST
    IF (.NOT.ANY(EVACUATION_ONLY)) RETURN
    IF (ICYC < 1) RETURN
    ! Check if FED is used
    IF (BTEST(I_MODE,3) .OR. BTEST(I_MODE,1)) EVAC_DEVICES(:)%USE_NOW = .FALSE.

  END SUBROUTINE CLEAN_AFTER_EVACUATE
!
  SUBROUTINE EVACUATE_HUMANS(TIN,DT,NM,ICYC)
    IMPLICIT NONE
    !
    ! Calculates the forces on humans and moves them.
    ! Uses a modified Velocity-Verlet algorithm.
    ! The modification is that the motive force part of
    ! the force is done using the dissipative self-consistent
    ! Velocity-Verlet algorithm by vattulainen et al.
    !
    ! INPUTS:
    !   TIN: End time, this routine makes the move TIN-DT ==> TIN (DT=FDS DT)
    !   NM: Mesh index
    !   ICYC: Index of the fds fire time step
    !
    ! Passed variables
    REAL(EB), INTENT(IN) :: TIN,DT
    INTEGER, INTENT(IN) :: NM,ICYC
    !
    ! Local variables
    INTEGER, PARAMETER :: N_SECTORS = 2
    REAL(EB) :: DTSP,UBAR,VBAR,X1,Y1,XI,YJ,ZK, WSPA=0._EB, WSPB=0._EB
    INTEGER :: ICN,I,J,IIN,JJN,KKN,II,JJ,KK,IIX,JJY,KKZ,ICX, ICY, N, J1, I_OBST, I_OBSTX, I_OBSTY
    INTEGER  :: IE, TIM_IC, TIM_IW, TIM_IWX, TIM_IWY, TIM_IW2, TIM_IC2, SURF_INDEX, NM_SEE
    REAL(EB) :: P2P_DIST, P2P_DIST_MAX, P2P_U, P2P_V, EVEL, TIM_DIST, EVEL2, MAX_V0_FAC
    !Issue1547: MAX_V0_FAC is declared as a real variable.
    REAL(EB), DIMENSION(4) :: D_XY
    LOGICAL, DIMENSION(4) :: FOUNDWALL_XY
    INTEGER :: ISTAT, STRS_INDX, I_TARGET, COLOR_INDEX, I_NEW_FFIELD, I_TARGET_OLD
    !
    !
    REAL(EB) ::  SCAL_PROD_OVER_RSQR, U_NEW, V_NEW, VMAX_TIMO, COSPHIFAC, &
         SPEED_MAX, DELTA_MIN, DT_SUM, C_YEFF, LAMBDAW, B_WALL, A_WALL, &
         T, CONTACT_F, SOCIAL_F, SMOKE_BETA, SMOKE_ALPHA, SMOKE_SPEED_FAC, tmp1, tmp2, CONTACT_FX, CONTACT_FY
    INTEGER :: IIE, JJE, IIO, JJO, III, JJJ, I_EGRID, I_TMP
    REAL(EB) :: X_NOW, Y_NOW, D_HUMANS, D_WALLS, DTSP_NEW, FAC_TIM, DT_GROUP_DOOR, X11, Y11, SPEED, TPRE
    LOGICAL PP_SEE_EACH, L_FIRST_PASS, USE_FED, PP_SEE_DOOR
    LOGICAL L_EFF_READ, L_EFF_SAVE, L_DEAD
    REAL(EB) :: COS_X, COS_Y, SPEED_XM, SPEED_XP, SPEED_YM, SPEED_YP, HR_Z, HR_A, HR_B, HR_TAU, HR_TAU_INER
    !
    REAL :: RN_REAL
    REAL(EB) :: RN, RNCF
    REAL(EB) :: GAME, GATH, GACM, TARGET_X, TARGET_Y
    !
    INTEGER :: I_OLD_FFIELD, IZERO, IMODE_OLD, I_STRS_DOOR, NM_STRS_INDEX, GUARD_MEN
    CHARACTER(LABEL_LENGTH) :: NAME_OLD_FFIELD
    !
    REAL(EB) :: P2P_TORQUE, FC_X, FC_Y, OMEGA_NEW, ANGLE, A1, TC_Z, FC_X1, FC_Y1
    REAL(EB) :: OMEGA_MAX, OMEGA_0, FAC_V0_UP, FAC_V0_DOWN, FAC_V0_HORI, TAU_FAC, SPEED_X, SPEED_Y
    REAL(EB) :: SMOKE_MIN_SPEED_TMP
    REAL(EB), DIMENSION(6) :: Y_TMP, X_TMP, R_TMP, V_TMP, U_TMP
    !
    ! Next are for the block list of the double agent-agent loop (speed up)
    INTEGER :: MAX_HUMANS_CELL, I_DX, J_DY, IE_MAX, BL_MAX
    REAL(EB) :: DX_MIN, DY_MIN
    INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: BLOCK_GRID
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: BLOCK_GRID_N
    INTEGER, DIMENSION(:), ALLOCATABLE :: BLOCK_LIST
    !
    REAL(EB), DIMENSION(:), ALLOCATABLE :: HERDING_LIST_DOORS
    INTEGER, DIMENSION(10) :: HERDING_LIST_IHUMAN
    REAL(EB), DIMENSION(10) :: HERDING_LIST_P2PDIST
    INTEGER :: HERDING_LIST_N
    REAL(EB) :: HERDING_LIST_P2PMAX, R_HERD_HR, DOT_HERD_HR, Other_TPRE
    !
    REAL(EB) :: D_HUMANS_MIN, D_WALLS_MIN
    REAL(EB) :: TNOW
    !
    LOGICAL :: NM_STRS_MESH
    REAL(EB), DIMENSION(N_SECTORS+1) :: SUM_SUUNTA, U_THETA, V_THETA, COS_THETA, SIN_THETA, THETAS
    REAL(EB) :: SUM_SUUNTA_MAX, THETA_STEP, THETA_START, VR_2R, ANGLE_HRE, &
         ANGLE_HR, V_HR, V_HRE, P2P_SUUNTA_MAX, ANGLE_OLD, HR_X, HR_Y, D_NEW, D_SHIFT, COMMITMENT
    INTEGER :: I_SUUNTA_MAX, N_SUUNTA_BACK, N_SUUNTA_BACKCF, I_COUNT_DENSITY, I_COUNT_DENSITYR, I_COUNT_DENSITYL
    INTEGER, DIMENSION(N_SECTORS+1) :: N_SUUNTA, N_SUUNTACF
    REAL(EB) :: U_i_Hawk, U_I_Dove, x_o, y_o, T_tmp, T_tmp1, Width, p_i, HR_GAMMA, HR_KAPPA
    INTEGER :: N_queue, N_i_Hawk, N_I_Dove, N_i_Hawk2, N_I_Dove2
    LOGICAL :: L_DO_GAME, L_FALLEN_DOWN, L_HRE_FALLEN_DOWN, L_FALLDOWN_THIS_TIME_STEP
    !Issue1547: Some variables that are needed
    REAL(EB) :: Eta_Nervous, Speed_Nervous, HR_Speed

    TYPE (MESH_TYPE),        POINTER :: MFF=>NULL()
    TYPE (EVAC_STRS_TYPE),   POINTER :: STRP=>NULL()
    TYPE (EVAC_SSTAND_TYPE), POINTER :: ESS=>NULL()
    TYPE (HUMAN_TYPE),       POINTER :: HR=>NULL(), HRE=>NULL()
    TYPE (EVACUATION_TYPE),  POINTER :: HPT =>NULL()
    TYPE (EVAC_ENTR_TYPE),   POINTER :: PNX =>NULL()
    TYPE (EVAC_EDEV_TYPE),   POINTER :: EDV=>NULL()
    TYPE (EVAC_HOLE_TYPE),   POINTER :: EHX=>NULL()
    !
    INTRINSIC :: BTEST
    !
    IF (.NOT.(EVACUATION_ONLY(NM) .AND. EMESH_INDEX(NM)>0)) RETURN
    TNOW=CURRENT_TIME()
    ! Check if FED is used
    USE_FED = .FALSE.
    IF (BTEST(I_EVAC,3) .OR. BTEST(I_EVAC,1)) USE_FED = .TRUE.
    !
    ! Maximun speed of the agents.
    VMAX_TIMO      = V_MAX
    ! Maximum angular speed of the agents.
    OMEGA_MAX      = V_ANGULAR_MAX*2.0_EB*PI ! 8 rounds per second in radians
    ! The target angular speed of the agents (like v0 for translational motion)
    OMEGA_0        = V_ANGULAR*2.0_EB*PI     ! 2 rounds per second in radians
    !
    DT_SUM         = 0.0_EB
    ! DT_GROUP_DOOR: How often (on the average) an agent (or group) tries to change the door.
    DT_GROUP_DOOR  = TAU_CHANGE_DOOR
    !
    CALL POINT_TO_MESH(NM)
    NM_SEE = NMESHES - N_EGRIDS + EMESH_INDEX(NM)
    !
    ! Find the smallest DELTA_X and/or DELTA_Y.  Note: Evac meshes should not be stretched ones.
    DX_MIN = MINVAL(DX)
    DY_MIN = MINVAL(DY)
    DELTA_MIN = MIN(DY_MIN, DX_MIN)
    !
    ! Read/write the EFF-file?
    L_EFF_READ = BTEST(I_EVAC,2)
    L_EFF_SAVE = BTEST(I_EVAC,0)
    !
    ! Find the egrid index of this main evac mesh. [1, # main evac meshes]
    I_EGRID = 0
    DO I = 1, N_EGRIDS
       IF (EVAC_NODE_LIST(I)%IMESH == NM) THEN
          I_EGRID = EVAC_NODE_LIST(I)%NODE_INDEX
       END IF
    END DO
    IF (I_EGRID == 0) THEN
       WRITE(LU_ERR,'(A,I6)') 'ERROR: EVACUATE_HUMANS, No mesh found ',NM
       WRITE(LU_EVACOUT,'(A,I6)') 'ERROR: EVACUATE_HUMANS, No mesh found ',NM
       STOP_STATUS = EVACUATION_STOP ; RETURN
    END IF
    NM_STRS_MESH = .FALSE.
    NM_STRS_INDEX = 0
    CHECK_STRS_LOOPM0: DO I = 1, N_STRS
       IF (EVAC_STRS(I)%IMESH == NM) THEN
          NM_STRS_MESH = .TRUE.
          NM_STRS_INDEX = I
          EXIT CHECK_STRS_LOOPM0
       END IF
    END DO CHECK_STRS_LOOPM0
    !
    ! Blocks are used to speed the double loops over the agent-agent interactions.
    ALLOCATE(BLOCK_GRID_N(1:IBAR,1:JBAR),STAT=IZERO)
    CALL CHKMEMERR('EVACUATE_HUMANS','BLOCK_GRID_N',IZERO)

    ALLOCATE(HERDING_LIST_DOORS(MAX(1,N_DOORS+N_EXITS)),STAT=IZERO)
    CALL CHKMEMERR('EVACUATE_HUMANS','HERDING_LIST_DOORS',IZERO)
    ! Initialize some counters etc. for this main evac mesh.
    L_FIRST_PASS = .TRUE.
    HUMAN_TIME_LOOP: DO WHILE ( DT_SUM < DT )
       ! DT is the fds flow calculation time step.
       ! Sometimes agent time step is smaller than the fire time step, so syncronize the agent clock correctly.
       ! TSTEPS: Time steps of the agent movement algorithm for different meshes, which were calculated at the
       !         end of the previous main time step.
       ! DTSP: The present time step for the agent movement algorithm, which may be smaller than the fds main time step DT.
       DTSP = MIN( (DT-DT_SUM), TSTEPS(NM) )
       DT_SUM = DT_SUM + DTSP
       T = TIN - DT + DT_SUM     ! Current time for the agents

       DO I = 1, N_DOORS
          IF (EVAC_DOORS(I)%IMESH == NM .AND. EVAC_DOORS(I)%ICORR_TO == 0) THEN
             IF (EVAC_DOORS(I)%TIME_OPEN > EVAC_DOORS(I)%TIME_CLOSE) THEN
                IMODE_OLD = EVAC_DOORS(I)%IMODE
                IF (T >  EVAC_DOORS(I)%TIME_CLOSE .AND. T < EVAC_DOORS(I)%TIME_OPEN .AND. IMODE_OLD == -1) &
                     EVAC_DOORS(I)%IMODE=+2
                IF (T >= EVAC_DOORS(I)%TIME_OPEN .AND. IMODE_OLD == -2) EVAC_DOORS(I)%IMODE=+1
             ELSE
                IMODE_OLD = EVAC_DOORS(I)%IMODE
                IF (T >  EVAC_DOORS(I)%TIME_OPEN .AND. T < EVAC_DOORS(I)%TIME_CLOSE .AND. IMODE_OLD == -2) &
                     EVAC_DOORS(I)%IMODE=+1
                IF (T >= EVAC_DOORS(I)%TIME_CLOSE .AND. IMODE_OLD==-1) EVAC_DOORS(I)%IMODE=+2
             END IF
          END IF
          IF (EVAC_DOORS(I)%IMESH == NM .AND. EVAC_DOORS(I)%ICORR_TO > 0) THEN
             ! Target of the door is an elevator, the corr index is doors%ICORR_TO
             IF (EVAC_DOORS(I)%TIME_OPEN > EVAC_DOORS(I)%TIME_CLOSE) THEN
                ! Cycle: T_BEGIN open, close at TIME_CLOSE, open TIME_OPEN and stay open
                IMODE_OLD = EVAC_DOORS(I)%IMODE
                IF ((T > EVAC_DOORS(I)%TIME_CLOSE .AND. T <= EVAC_DOORS(I)%TIME_OPEN) .OR. &
                     EVAC_CORRS(EVAC_DOORS(I)%ICORR_TO)%Elevator_Mode /= 2) THEN
                   ! The elevator is not in use by the user input or it is not at the floor
                   IF (IMODE_OLD == -1) EVAC_DOORS(I)%IMODE=+2 ! close the door
                   IF (IMODE_OLD == +1) EVAC_DOORS(I)%IMODE=-2 ! do not open, let stay closed
                ELSE
                   ! The elevator is in use by the user input and it is at the floor
                   IF (EVAC_CORRS(EVAC_DOORS(I)%ICORR_TO)%Elevator_Mode == 2) THEN
                      IF (IMODE_OLD == -2) EVAC_DOORS(I)%IMODE=+1 ! open the door
                      IF (IMODE_OLD == +2) EVAC_DOORS(I)%IMODE=-1 ! stay open
                   END IF
                END IF
             ELSE
                ! Cycle: T_BEGIN closed, open at TIME_OPEN, close TIME_CLOSE and stay closed
                IMODE_OLD = EVAC_DOORS(I)%IMODE
                IF (T < EVAC_DOORS(I)%TIME_OPEN .OR. T > EVAC_DOORS(I)%TIME_CLOSE .OR. &
                     EVAC_CORRS(EVAC_DOORS(I)%ICORR_TO)%Elevator_Mode /= 2) THEN
                   ! The elevator is not in use by the user input or it is not at the floor
                   IF (IMODE_OLD == -1) EVAC_DOORS(I)%IMODE=+2 ! close the door
                   IF (IMODE_OLD == +1) EVAC_DOORS(I)%IMODE=-2 ! do not open, let stay closed
                ELSE
                   ! The elevator is in use by the user input and it is at the floor
                   IF (EVAC_CORRS(EVAC_DOORS(I)%ICORR_TO)%Elevator_Mode == 2) THEN
                      IF (IMODE_OLD == -2) EVAC_DOORS(I)%IMODE=+1 ! open the door
                      IF (IMODE_OLD == +2) EVAC_DOORS(I)%IMODE=-1 ! stay open
                   END IF
                END IF
             END IF
          END IF
       END DO
       DO I = 1, N_EXITS
          IF (EVAC_EXITS(I)%IMESH == NM .AND. .NOT.EVAC_EXITS(I)%COUNT_ONLY) THEN
             IF (EVAC_EXITS(I)%TIME_OPEN > EVAC_EXITS(I)%TIME_CLOSE) THEN
                IMODE_OLD = EVAC_EXITS(I)%IMODE
                IF (T >  EVAC_EXITS(I)%TIME_CLOSE .AND. T < EVAC_EXITS(I)%TIME_OPEN .AND. IMODE_OLD==-1) &
                     EVAC_EXITS(I)%IMODE=+2
                IF (T >= EVAC_EXITS(I)%TIME_OPEN .AND. IMODE_OLD==-2) EVAC_EXITS(I)%IMODE=+1
             ELSE
                IMODE_OLD = EVAC_EXITS(I)%IMODE
                IF (T >  EVAC_EXITS(I)%TIME_OPEN .AND. T < EVAC_EXITS(I)%TIME_CLOSE .AND. IMODE_OLD==-2) &
                     EVAC_EXITS(I)%IMODE=+1
                IF (T >= EVAC_EXITS(I)%TIME_CLOSE .AND. IMODE_OLD==-1) EVAC_EXITS(I)%IMODE=+2
             END IF
          END IF
       END DO
       !

       DO I = 1, N_DOORS
          IF (EVAC_DOORS(I)%IMESH == NM) THEN
             EVAC_DOORS(I)%NTARGET(1:50) = 0
          END IF
       END DO
       DO I = 1, N_EXITS
          IF (EVAC_EXITS(I)%IMESH == NM) THEN
             EVAC_EXITS(I)%NTARGET(1:50) = 0
          END IF
       END DO
       DO I = 1, N_ENTRYS
          IF (EVAC_ENTRYS(I)%IMESH == NM .AND. EVAC_ENTRYS(I)%CROWBAR_READ_IN) THEN
             EVAC_ENTRYS(I)%CB_N_Agents = 0
             EVAC_ENTRYS(I)%CB_N_CrowbarAgents = 0
          END IF
       END DO

       ! ================================================
       ! Initialize group arrays for this main evac mesh.
       ! ================================================
       GROUP_LIST(:)%GROUP_SIZE      = 0
       GROUP_LIST(:)%GROUP_X         = 0.0_EB
       GROUP_LIST(:)%GROUP_Y         = 0.0_EB
       GROUP_LIST(:)%MAX_DIST_CENTER = 0.0_EB
       GROUP_LIST(:)%SPEED           = 0.0_EB
       GROUP_LIST(:)%INTDOSE         = 0.0_EB
       GROUP_LIST(:)%TPRE            = 0.0_EB
       GROUP_LIST(:)%TDET            = HUGE(GROUP_LIST(:)%TDET)
       !
       DO J = 0, I33_DIM
          GROUP_LIST(J)%GROUP_I_FFIELDS(I_EGRID) = 0
       END DO
       IF (N_DOORS >0) EVAC_DOORS(:)%R_NTARGET = 5.0_EB
       IF (N_EXITS >0) EVAC_EXITS(:)%R_NTARGET = 5.0_EB
       DO I = 1, N_HUMANS
          HR=>HUMAN(I)
          J = MAX(0,HR%GROUP_ID)  ! Group index of the agent
          GROUP_LIST(J)%GROUP_SIZE = GROUP_LIST(J)%GROUP_SIZE + 1
          GROUP_LIST(J)%GROUP_X    = GROUP_LIST(J)%GROUP_X + HR%X
          GROUP_LIST(J)%GROUP_Y    = GROUP_LIST(J)%GROUP_Y + HR%Y
          GROUP_LIST(J)%SPEED      = GROUP_LIST(J)%SPEED + HR%SPEED
          GROUP_LIST(J)%INTDOSE    = GROUP_LIST(J)%INTDOSE + HR%INTDOSE
          GROUP_LIST(J)%TPRE       = MAX(GROUP_LIST(J)%TPRE,HR%TPRE)
          GROUP_LIST(J)%TDET       = MIN(GROUP_LIST(J)%TDET,HR%TDET)
          I_TMP = HR%I_TARGET
          IF (I_TMP > 0 .AND. I_TMP <= N_DOORS ) THEN
             IF ((NM_STRS_MESH .AND. HR%STR_SUB_INDX==EVAC_DOORS(I_TMP)%STR_SUB_INDX) .OR. .NOT.NM_STRS_MESH) THEN
                EVEL = SQRT((HR%X- EVAC_DOORS(I_TMP)%X)**2 + (HR%Y-EVAC_DOORS(I_TMP)%Y)**2)
                EVAC_DOORS(I_TMP)%R_NTARGET = MAX(EVEL,EVAC_DOORS(I_TMP)%R_NTARGET)
             END IF
          END IF
          I_TMP = I_TMP - N_DOORS
          IF (I_TMP > 0 .AND. I_TMP <= N_EXITS ) THEN
             IF ((NM_STRS_MESH .AND. HR%STR_SUB_INDX==EVAC_EXITS(I_TMP)%STR_SUB_INDX) .OR. .NOT.NM_STRS_MESH) THEN
                EVEL = SQRT((HR%X-EVAC_EXITS(I_TMP)%X)**2 + (HR%Y-EVAC_EXITS(I_TMP)%Y)**2)
                EVAC_EXITS(I_TMP)%R_NTARGET = MAX(EVEL,EVAC_EXITS(I_TMP)%R_NTARGET)
             END IF
          END IF
       END DO

       GROUP_LIST(1:)%GROUP_X = GROUP_LIST(1:)%GROUP_X / MAX(1,GROUP_LIST(1:)%GROUP_SIZE)
       GROUP_LIST(1:)%GROUP_Y = GROUP_LIST(1:)%GROUP_Y / MAX(1,GROUP_LIST(1:)%GROUP_SIZE)
       GROUP_LIST(1:)%SPEED   = GROUP_LIST(1:)%SPEED   / MAX(1,GROUP_LIST(1:)%GROUP_SIZE)
       GROUP_LIST(1:)%INTDOSE = GROUP_LIST(1:)%INTDOSE / MAX(1,GROUP_LIST(1:)%GROUP_SIZE)
       GROUP_LIST(:)%MAX_DIST_CENTER = 0.0_EB

       ! GROUP_LIST intitialization and count the number of agents going towards different doors/exits.
       DO I = 1, N_HUMANS
          HR=>HUMAN(I)
          ! GROUP_ID > 0: +GROUP_ID
          ! GROUP_ID < 0: -HUMAN_ID (lonely agents)
          J  =  MAX(0,HR%GROUP_ID)  ! group index > 0
          J1 = -MIN(0,HR%GROUP_ID)  ! lonely agent index > 0
          GROUP_LIST(J)%MAX_DIST_CENTER = MAX(GROUP_LIST(J)%MAX_DIST_CENTER, &
               SQRT((HR%X - GROUP_LIST(J)%GROUP_X)**2 + (HR%Y - GROUP_LIST(J)%GROUP_Y)**2))
          I_TMP = HR%I_TARGET
          IF (I_TMP > 0 .AND. I_TMP <= N_DOORS ) THEN
             IF ((NM_STRS_MESH .AND. HR%STR_SUB_INDX==EVAC_DOORS(I_TMP)%STR_SUB_INDX) .OR. .NOT.NM_STRS_MESH) THEN
                EVEL = SQRT((HR%X-0.5_EB*(EVAC_DOORS(I_TMP)%X1+EVAC_DOORS(I_TMP)%X2))**2 + &
                     (HR%Y-0.5_EB*(EVAC_DOORS(I_TMP)%Y1+EVAC_DOORS(I_TMP)%Y2))**2)
                EVEL = 50.0_EB*EVEL/EVAC_DOORS(I_TMP)%R_NTARGET + 1.0_EB
                IE = MIN(50,MAX(1,INT(EVEL)))
                EVAC_DOORS(I_TMP)%NTARGET(IE:50) = EVAC_DOORS(I_TMP)%NTARGET(IE:50) + 1
             END IF
          END IF
          I_TMP = I_TMP - N_DOORS
          IF (I_TMP > 0 .AND. I_TMP <= N_EXITS ) THEN
             IF ((NM_STRS_MESH .AND. HR%STR_SUB_INDX==EVAC_EXITS(I_TMP)%STR_SUB_INDX) .OR. .NOT.NM_STRS_MESH) THEN
                EVEL = SQRT((HR%X-0.5_EB*(EVAC_EXITS(I_TMP)%X1+EVAC_EXITS(I_TMP)%X2))**2 + &
                     (HR%Y-0.5_EB*(EVAC_EXITS(I_TMP)%Y1+EVAC_EXITS(I_TMP)%Y2))**2)
                EVEL = 50.0_EB*EVEL/EVAC_EXITS(I_TMP)%R_NTARGET + 1.0_EB
                IE = MIN(50,MAX(1,INT(EVEL)))
                EVAC_EXITS(I_TMP)%NTARGET(IE:50) = EVAC_EXITS(I_TMP)%NTARGET(IE:50) + 1
             END IF
          END IF
       END DO

       ! J=0, i.e., lonely humans: GROUP_LIST is not used , but it is safe to initialize.
       IF (N_HUMANS > 0) THEN
          GROUP_LIST(0)%GROUP_SIZE = 1
          GROUP_LIST(0)%COMPLETE   = 1
          GROUP_LIST(0)%GROUP_X    = 0.5_EB*(XS+XF)
          GROUP_LIST(0)%GROUP_Y    = 0.5_EB*(YS+YF)
          GROUP_LIST(0)%SPEED      = 1.0_EB
          GROUP_LIST(0)%INTDOSE    = 0.0_EB
          GROUP_LIST(0)%MAX_DIST_CENTER = 0.0_EB
       END IF

       ! Note: Group_size is the number of group members on this main evac mesh.
       ! Check if the groups are already together or not.
       DO J = 1, I33_DIM
          GROUP_LIST(J)%LIMIT_COMP = RADIUS_COMPLETE_0 + RADIUS_COMPLETE_1*GROUP_LIST(J)%GROUP_SIZE
          IF ( ((GROUP_LIST(J)%MAX_DIST_CENTER <=  GROUP_LIST(J)%LIMIT_COMP) .OR. &
               (GROUP_LIST(J)%COMPLETE == 1)) .AND. GROUP_LIST(J)%GROUP_SIZE > 0) THEN
             ! If COMPLETE=1 already, it stays at 1.
             ! If many floors, check only floors where there are group members.
             ! Group may become complete before fire is detected. Do not count these as complete.
             IF (T > GROUP_LIST(J)%TDET) THEN
                IF (GROUP_LIST(J)%COMPLETE == 0) THEN
                   ! TDOOR: Saves the time point when the group started to move towards the exit.
                   GROUP_LIST(J)%TDOOR = MAX(T,GROUP_LIST(J)%TDET)
                END IF
                GROUP_LIST(J)%COMPLETE = 1
             END IF
          END IF
       END DO

       Door_Mode_Loop: DO I = 1, N_HUMANS
          ! Check if the agent is moving towards a door or other target.
          ! Other agents can only follow agents moving towards doors ("herding mode").
          ! I_DoorAlgo=0: Use the main evacuation mesh flow field ("stupid agents"), do not follow these.
          ! I_Door_Mode:  0) no target door or not (yet) moving
          !               1) target door + moving by using flow fields
          !               2) target door + moving + bee line towards the exit (do not use flow fields) XYZ visible
          !               3) target door + moving + bee line towards the exit (do not use flow fields) XB visible
          !               4) target door + moving + bee line towards the exit (do not use flow fields) XYZ+XB visible
          !              <0) came out form a door or an entry, but does not have a target door
          HR=>HUMAN(I)
          IF (HR%I_DoorAlgo == 0) CYCLE Door_Mode_Loop
          J  =  MAX(0,HR%GROUP_ID)  ! group index > 0 (=0 for lonely agents)
          J1 = -MIN(0,HR%GROUP_ID)  ! lonely agent index > 0
          IF (J1 > 0 .AND. (HR%I_Target /= 0 .AND. HR%I_DoorAlgo > 0) .AND. T > HR%TPRE+HR%TDET) &
               HR%I_Door_Mode = MAX(1,HR%I_Door_Mode)
          IF (J > 0 .AND. (HR%I_Target /= 0 .AND. HR%I_DoorAlgo > 0) .AND. T > GROUP_LIST(J)%TPRE + GROUP_LIST(J)%TDOOR) &
               HR%I_Door_Mode = MAX(1,HR%I_Door_Mode)
       END DO Door_Mode_Loop

       ! ========================================================
       ! Change target door?
       ! ========================================================
       IF (T >= TWO_EPSILON_EB .AND. .NOT.NM_STRS_MESH) THEN
          CHANGE_DOOR_LOOP: DO I = 1, N_HUMANS
             HR => HUMAN(I)

             I_OLD_FFIELD    = HR%I_FFIELD
             NAME_OLD_FFIELD = TRIM(MESH_NAME(I_OLD_FFIELD))
             J  =  MAX(0,HR%GROUP_ID)   ! Group index
             J1 = -MIN(0,HR%GROUP_ID)   ! Lonely agent index

             ! If the group is not yet together, they are not moving towards the exit.
             IF (GROUP_LIST(J)%COMPLETE == 0) CYCLE CHANGE_DOOR_LOOP

             ! Agents start to move towards the exit after the pre-evacuation (reaction) time.

             IF (GROUP_LIST(J)%GROUP_I_FFIELDS(I_EGRID) == 0) THEN
                ! Test if this group/agent should update the exit door (exp. distribution)
                IF (ABS(HR%I_TARGET) > 0) THEN
                   ! Check, if the target door is closed.
                   IF (ABS(HR%I_TARGET) <= N_DOORS) THEN
                      IF (ABS(EVAC_DOORS(ABS(HR%I_TARGET))%IMODE)==2 .AND. &
                           .NOT.EVAC_DOORS(ABS(HR%I_TARGET))%TARGET_WHEN_CLOSED) THEN
                         HR%I_TARGET = 0
                         !RN = 2.0_EB ! update door, present one is closed
                      END IF
                   ELSE
                      IF (ABS(EVAC_EXITS(ABS(HR%I_TARGET)-N_DOORS)%IMODE)==2 .AND. &
                           .NOT.EVAC_EXITS(ABS(HR%I_TARGET)-N_DOORS)%TARGET_WHEN_CLOSED) THEN
                         HR%I_TARGET = 0
                         !RN = 2.0_EB ! update door, present one is closed
                      END IF
                   END IF
                END IF
                IF (J == 0 .AND. T < HR%TPRE+HR%TDET .AND. HR%I_Target/=0) CYCLE CHANGE_DOOR_LOOP   ! LONELY AGENTS
                IF (J > 0 .AND. T < GROUP_LIST(J)%TPRE + GROUP_LIST(J)%TDOOR) CYCLE CHANGE_DOOR_LOOP ! GROUPS
                CALL RANDOM_NUMBER(RN_REAL)
                RN = REAL(RN_REAL,EB)
                IF (RN > EXP(-DTSP/DT_GROUP_DOOR) ) THEN
                   I_TARGET_OLD = HR%I_TARGET
                   N_CHANGE_TRIALS = N_CHANGE_TRIALS + 1
                   CALL CHANGE_TARGET_DOOR(NM,NM,I,J,J1,I_EGRID,1,HR%X,HR%Y,I_TARGET,COLOR_INDEX,I_NEW_FFIELD,HR)
                   IF (STOP_STATUS>NO_STOP) RETURN
                   IF (ABS(I_TARGET_OLD) > 0) THEN
                      IF (.NOT.Is_Visible_Door(ABS(I_TARGET_OLD)) .AND. .NOT.Is_Known_Door(ABS(I_TARGET_OLD))) THEN
                         I_TARGET_OLD = 0
                      END IF
                   END IF
                   IF (I_TARGET /= 0) THEN  ! One should put this inside change_target_door
                      IF (Is_BeeLine_Door(ABS(I_TARGET),2))   HR%SeeDoorXYZ1=.TRUE.
                      IF (Is_BeeLine_Door(ABS(I_TARGET),3))   HR%SeeDoorXYZ2=.TRUE.
                      IF (Is_BeeLineXB_Door(ABS(I_TARGET),2)) HR%SeeDoorXB1=.TRUE.
                      IF (Is_BeeLineXB_Door(ABS(I_TARGET),3)) HR%SeeDoorXB2=.TRUE.
                      IF (Is_BeeLine_Door(ABS(I_TARGET),1) .OR. Is_BeeLineXB_Door(ABS(I_TARGET),1)) THEN
                         IF (HR%D_Walls > 4.0_EB .AND. HR%DENSITY < 2.0_EB) THEN
                            HR%I_Door_Mode = HR%I_Door_Mode
                         ELSE
                            IF (ABS(I_TARGET) > N_DOORS) THEN
                               N = ABS(I_TARGET) - N_DOORS
                               IF (ABS(EVAC_EXITS(N)%IMODE)==2 .AND. EVAC_EXITS(N)%TARGET_WHEN_CLOSED .AND. &
                                    EVAC_EXITS(N)%WAIT_AT_XYZ) THEN
                                  HR%I_Door_Mode = HR%I_Door_Mode
                               ELSE
                                  HR%I_Door_Mode = 1
                               END IF
                            ELSE
                               N = ABS(I_TARGET)
                               IF (ABS(EVAC_DOORS(N)%IMODE)==2 .AND. EVAC_DOORS(N)%TARGET_WHEN_CLOSED .AND. &
                                    EVAC_DOORS(N)%WAIT_AT_XYZ) THEN
                                  HR%I_Door_Mode = HR%I_Door_Mode
                               ELSE
                                  HR%I_Door_Mode = 1
                               END IF
                            END IF
                         END IF
                         ! HR%I_Door_Mode = 2   ! Use bee line towards the door
                      ELSE
                         HR%I_Door_Mode = 1   ! Use flow field towards the door
                      END IF
                   END IF

                   ! HR%I_DoorAlgo: 1 rational agents, 2 known doors, 3 herding, 0 main evac mesh ff
                   ! Herding behaviour: Count and weight by the distance the target doors of the other
                   ! agents around this agent.  The five (or less) nearest neighbor agents are used.
                   ! I_HERDING_TYPE 0: default herding, 1: keep the first choice, 2: do not move,
                   !                3: do not move + keep the first choice
                   Herding_Type_Agent: IF (HR%I_DoorAlgo>=3 .OR. (HR%I_Target==0 .AND. HR%I_DoorAlgo>0)) THEN
                      IF (HR%I_DoorAlgo==3 .AND. HR%I_Target/=0 .AND. (I_HERDING_TYPE==1 .OR. I_HERDING_TYPE==3)) THEN
                         ! These herding agents keep the first choise
                         I_TARGET = I_TARGET_OLD
                         IF (ABS(I_TARGET) > 0) HR%I_Door_Mode = 1
                         IF (ABS(I_TARGET) > 0 .AND. Is_BeeLine_Door(MAX(1,ABS(I_TARGET)),1)) HR%I_Door_Mode = 2
                         IF (ABS(I_TARGET) > 0 .AND. Is_BeeLineXB_Door(MAX(1,ABS(I_TARGET)),1)) HR%I_Door_Mode = HR%I_Door_Mode + 2
                         IF (HR%D_Walls <= 4.0_EB .OR. HR%DENSITY >= 2.0_EB) HR%I_Door_Mode = 1
                         ! IF (HR%D_Walls <= 4.0_EB .AND. HR%DENSITY >= 2.0_EB) HR%I_Door_Mode = 1
                      ELSE
                         HERDING_LIST_DOORS = 0.0_EB    ! Real array, weights by the p2p distance (linear function)
                         HERDING_LIST_IHUMAN = 0        ! indices of the ten closest other agents
                         ! Now 5 neighbors are checked (for future the array has maximum of 10 places)
                         HERDING_LIST_P2PDIST(1:5) = HUGE(HERDING_LIST_P2PMAX) ! person to person distances (squares)
                         HERDING_LIST_P2PDIST(6:10) = 0.0_EB ! person to person distances (squares)
                         HERDING_LIST_N = 0 ! number of neighbors (max 10, see dimension)
                         HERDING_LIST_P2PMAX = HUGE(HERDING_LIST_P2PMAX)
                         TIM_DIST = HUGE(TIM_DIST)
                         R_HERD_HR = R_HERDING
                         DOT_HERD_HR = DOT_HERDING
                         IF (HR%I_DoorAlgo == 4) THEN
                            R_HERD_HR = MIN(2.5_EB, R_HERDING)
                            DOT_HERD_HR = 0.2_EB
                         END IF
                         Other_Agent_Loop: DO IE = 1, N_HUMANS
                            IF (I==IE) CYCLE Other_Agent_Loop
                            HRE => HUMAN(IE)
                            P2P_DIST = (HRE%X-HR%X)**2 + (HRE%Y-HR%Y)**2
                            IF (P2P_DIST > R_HERD_HR**2 .AND. HRE%I_Target==0) CYCLE Other_Agent_Loop
                            IF (P2P_DIST > R_HERD_HR**2 .AND. HRE%I_Door_Mode <= 0) CYCLE Other_Agent_Loop
                            ! Check, that the persons are seeing each other, i.e., there are no walls between.
                            PP_SEE_EACH = SEE_EACH_OTHER(NM_SEE, HR%X, HR%Y, HRE%X, HRE%Y)
                            IF (.NOT. PP_SEE_EACH) CYCLE Other_Agent_Loop
                            IF (ABS(HRE%I_Target) > 0 .AND. P2P_DIST < TIM_DIST) THEN
                               IF (ABS(HRE%I_TARGET) <= N_DOORS) THEN
                                  x11 = EVAC_DOORS(ABS(HRE%I_TARGET))%X
                                  y11 = EVAC_DOORS(ABS(HRE%I_TARGET))%Y
                               ELSE
                                  x11 = EVAC_EXITS(ABS(HRE%I_TARGET)-N_DOORS)%X
                                  y11 = EVAC_EXITS(ABS(HRE%I_TARGET)-N_DOORS)%Y
                               END IF
                               PP_SEE_DOOR = SEE_DOOR(NM, HR%X, HR%Y, x11, y11, tmp1, tmp2)
                               IF (PP_SEE_DOOR) THEN
                                  EVEL  = SQRT((HR%X-x11)**2  + (HR%Y-y11)**2)  ! L2 norm
                                  EVEL2 = SQRT((HRE%X-x11)**2 + (HRE%Y-y11)**2) ! L2 norm
                               ELSE
                                  EVEL  = ABS(HR%X-x11)  + ABS(HR%Y-y11)  ! L1 norm
                                  EVEL2 = ABS(HRE%X-x11) + ABS(HRE%Y-y11) ! L1 norm
                               END IF
                               x_now =((HR%X-HRE%X)*HRE%UBAR + (HR%Y-HRE%Y)*HRE%VBAR) / ( &
                                 MAX(0.01_EB,SQRT(HRE%UBAR**2+HRE%VBAR**2))*SQRT((HRE%X-HR%X)**2+(HRE%Y-HR%Y)**2) )
                               ! The other agent should be closer to the door if the door is visible.
                               ! If it is not visible the other agent should go away.
                               IF ((EVEL2 < 0.5_EB*EVEL .AND. PP_SEE_DOOR) .OR. &
                                    ( (x_now < 0.0_EB .OR. EVEL2 < 0.9_EB*EVEL) .AND. .NOT.PP_SEE_DOOR)) THEN
                                  TIM_DIST = P2P_DIST
                                  IF (PP_SEE_DOOR) HR%I_Target2 = ABS(HRE%I_Target)
                                  IF (.NOT.PP_SEE_DOOR) HR%I_Target2 = -ABS(HRE%I_Target)
                               END IF
                            END IF
                            IF (P2P_DIST > R_HERD_HR**2) CYCLE Other_Agent_Loop
                            ! If not yet moving phase, do not follow others
                            IF (J == 0 .AND. T < HR%TPRE+HR%TDET) CYCLE Other_Agent_Loop
                            IF (J > 0 .AND. T < GROUP_LIST(J)%TPRE + GROUP_LIST(J)%TDOOR) CYCLE Other_Agent_Loop
                            IF (P2P_DIST <= HERDING_LIST_P2PMAX) THEN
                               ! Dot product used to check if the other agent is walking towards the agent
                               IF (HERDING_LIST_N < 5) THEN
                                  ! List is not yet full, all agents are ok
                                  HERDING_LIST_N = HERDING_LIST_N + 1
                                  HERDING_LIST_IHUMAN(HERDING_LIST_N) = +IE
                                  HERDING_LIST_P2PDIST(HERDING_LIST_N) = P2P_DIST
                               ELSE
                                  ! There are already five agents in the list, find the furthest one
                                  I_TMP = MAXVAL(MAXLOC(HERDING_LIST_P2PDIST))
                                  HERDING_LIST_IHUMAN(I_TMP) = +IE
                                  HERDING_LIST_P2PDIST(I_TMP) = P2P_DIST
                               END IF
                               HERDING_LIST_P2PMAX = MAXVAL(HERDING_LIST_P2PDIST)
                            END IF
                         END DO Other_Agent_Loop
                         Other_Agent_Loop_2: DO IE = 1, HERDING_LIST_N
                            HRE => HUMAN(HERDING_LIST_IHUMAN(IE))
                            IF (HRE%I_Door_Mode <= 0) CYCLE Other_Agent_Loop_2
                            IF (ABS(HRE%I_Target)==0) CYCLE Other_Agent_Loop_2
                            IF (T < HRE%TPRE+HRE%TDET) CYCLE Other_Agent_Loop_2
                            EVEL = ((HR%X-HRE%X)*HRE%UBAR + (HR%Y-HRE%Y)*HRE%VBAR) / ( &
                                 MAX(0.01_EB,SQRT(HRE%UBAR**2+HRE%VBAR**2))*SQRT((HRE%X-HR%X)**2+(HRE%Y-HR%Y)**2) )
                            IF (EVEL > DOT_HERDING) CYCLE Other_Agent_Loop_2
                            ! Linear weight function
                            P2P_DIST = SQRT((HRE%X-HR%X)**2 + (HRE%Y-HR%Y)**2)
                            HERDING_LIST_DOORS(ABS(HRE%I_Target)) = HERDING_LIST_DOORS(ABS(HRE%I_Target)) + &
                                 W0_HERDING -((W0_HERDING-WR_HERDING)/R_HERD_HR)*P2P_DIST
                         END DO Other_Agent_Loop_2
                         
                         Other_TPRE = 0.0_EB
                         Other_Agent_Loop_3: DO IE = 1, HERDING_LIST_N
                            HRE => HUMAN(HERDING_LIST_IHUMAN(IE))
                            Other_TPRE = Other_TPRE + HRE%TPRE
                         END DO Other_Agent_Loop_3
                         Other_TPRE = Other_TPRE/REAL(MAX(1,HERDING_LIST_N))
                         
                         DO II = 1, N_DOORS+N_EXITS
                            IF (HERDING_LIST_DOORS(II)>0.0_EB) THEN
                               ! Make it symmetrical with respect the doors.
                               ! The order of the doors is now random if they are equal.
                               CALL RANDOM_NUMBER(RN_REAL)
                               RN = REAL(RN_REAL,EB)
                               HERDING_LIST_DOORS(II) = HERDING_LIST_DOORS(II) + 0.0001_EB*RN
                            END IF
                         END DO
                         I_TMP = 0
                         EVEL = 0.0_EB
                         DO II = 1, N_DOORS+N_EXITS
                            IF (HERDING_LIST_DOORS(II)>EVEL) THEN
                               IF (Is_BeeLine_Door(II,1)) THEN
                                  I_TMP = II
                               ELSE
                                  I_TMP = -II ! Set the door non-visible status
                               END IF
                               EVEL = HERDING_LIST_DOORS(II)
                            END IF
                         END DO
                         IF (I_TMP /= 0 .AND. HR%I_Door_Mode == 0) THEN
                            ! I_Door_Mode: 0: default, has not found a door by the algorithm or is not moving (Tpre+Tdet)
                            !              1: has selected the door by the door selection algorithm and is moving
                            !             <0: came out of a door/entr (-inode) and did not find a new target on this floor
                            ! Found a door using herding algorithm, start to move after one second.
                            ! Merge one's TPRE with others' TPRE
                            IF (HERDING_TPRE_PROB > TWO_EPSILON_EB) THEN
                               !HR%TPRE = 0.5_EB*HR%TPRE + 0.5_EB*Other_TPRE
                               HR%TPRE = (1.0_EB-HERDING_TPRE_PROB)*HR%TPRE + HERDING_TPRE_PROB*Other_TPRE
                            ELSE
                               HR%TPRE = DT_GROUP_DOOR ! Original FDS+Evac
                            END IF
                            HR%TDET = MIN(T,HR%TDET)
                         END IF
                         IF (HR%I_DoorAlgo == 4) THEN
                            ! These agents follow others if that door is better than their own choice.
                            ! How: the i_target2 is set to this new door, where others are going and
                            ! it is set as a known door in the door selection subroutine.
                            IF (I_TMP /= 0) THEN
                               HR%I_Target2 = I_TMP ! Follow others, make the door known in door selection
                               IF (I_TARGET /= 0) I_TMP = I_TARGET ! No door, follow others right now
                            END IF
                            IF (I_TMP == 0 .AND. I_TARGET /= 0) THEN
                               ! If there is no one to follow, use the normal door.
                               ! If no door at all, use the i_target2 method
                               I_TMP = I_TARGET
                               HR%I_Target2 = 0
                            END IF
                         END IF
                         IF (I_TARGET_OLD==0 .AND. I_TMP==0 .AND. ABS(HR%I_Target2) > 0 .AND. T > HR%TPRE+HR%TDET) THEN
                            CALL RANDOM_NUMBER(RN_REAL)
                            RN = REAL(RN_REAL,EB)
                            IF (RN <= HERDING_TAU_FACTOR) I_TMP = HR%I_Target2
                         END IF
                         IF (ABS(I_TMP) > 0) THEN  ! Found neighbors
                            I_TARGET = I_TMP
                            IF (ABS(I_TMP) > N_DOORS) THEN
                               I_New_FField = EVAC_EXITS(ABS(I_TMP)-N_DOORS)%I_VENT_FFIELD
                            ELSE
                               I_New_FField = EVAC_DOORS(ABS(I_TMP))%I_VENT_FFIELD
                            END IF
                            HR%I_FFIELD = I_New_FField
                            HR%FFIELD_NAME = TRIM(MESH_NAME(I_New_FField))
                            HR%I_Door_Mode = 1
                            IF (Is_BeeLine_Door(ABS(I_TMP),1)) HR%I_Door_Mode = 2
                            IF (Is_BeeLineXB_Door(ABS(I_TMP),1)) HR%I_Door_Mode = HR%I_Door_Mode + 2
                            IF (HR%D_Walls <= 4.0_EB .OR. HR%DENSITY >= 2.0_EB) HR%I_Door_Mode = 1
                            ! IF (HR%D_Walls <= 4.0_EB .AND. HR%DENSITY >= 2.0_EB) HR%I_Door_Mode = 1
                         ELSE
                            IF (ABS(I_TARGET_OLD)>0) THEN  ! Go to the old target door
                               I_TARGET = I_TARGET_OLD
                               I_New_FField = I_OLD_FFIELD
                               HR%I_Door_Mode = 1
                               IF (Is_BeeLine_Door(ABS(I_TARGET),1)) HR%I_Door_Mode = 2
                               IF (Is_BeeLineXB_Door(ABS(I_TARGET),1)) HR%I_Door_Mode = HR%I_Door_Mode + 2
                               IF (HR%D_Walls <= 4.0_EB .OR. HR%DENSITY >= 2.0_EB) HR%I_Door_Mode = 1
                               ! IF (HR%D_Walls <= 4.0_EB .AND. HR%DENSITY >= 2.0_EB) HR%I_Door_Mode = 1
                            ELSE  ! no target door set, set the default flow field
                               I_TARGET = 0
                               IF (HR%IEL > 0) THEN
                                  ! The agent is from some evac line
                                  HPT => EVAC_EVACS(HR%IEL)
                                  IF (HPT%IMESH == NM) THEN
                                     I_NEW_FFIELD    = HPT%I_VENT_FFIELDS(0)
                                     HR%FFIELD_NAME = TRIM(MESH_NAME(I_NEW_FFIELD))
                                  ELSE
                                     I_NEW_FFIELD    = NM
                                     HR%FFIELD_NAME = TRIM(MESH_NAME(NM))
                                  END IF
                               ELSE
                                  ! The agent is from some entr line
                                  PNX => EVAC_ENTRYS(ABS(HR%IEL))
                                  IF (PNX%IMESH == NM) THEN
                                     I_NEW_FFIELD    = PNX%I_VENT_FFIELDS(0)
                                     HR%FFIELD_NAME = TRIM(MESH_NAME(I_NEW_FFIELD))
                                  ELSE
                                     I_NEW_FFIELD    = NM
                                     HR%FFIELD_NAME = TRIM(MESH_NAME(NM))
                                  END IF
                               END IF
                               IF (HR%I_Door_Mode < 0) THEN
                                  IF (EVAC_NODE_LIST(ABS(HR%I_Door_Mode))%NODE_TYPE == 'Door') THEN
                                     HR%FFIELD_NAME = &
                                          TRIM(EVAC_DOORS(EVAC_NODE_LIST(ABS(HR%I_Door_Mode))%NODE_INDEX)%GRID_NAME)
                                     II_LOOP1: DO II = 1, NMESHES
                                        IF (EVACUATION_ONLY(II) .AND. TRIM(HR%FFIELD_NAME) == TRIM(MESH_NAME(II))) THEN
                                           I_NEW_FFIELD = II
                                           EXIT II_LOOP1
                                        END IF
                                     END DO II_LOOP1
                                  END IF
                                  IF (EVAC_NODE_LIST(ABS(HR%I_Door_Mode))%NODE_TYPE == 'Entry') THEN
                                     HR%FFIELD_NAME = &
                                          TRIM(EVAC_ENTRYS(EVAC_NODE_LIST(ABS(HR%I_Door_Mode))%NODE_INDEX)%GRID_NAME)
                                     II_LOOP2: DO II = 1, NMESHES
                                        IF (EVACUATION_ONLY(II) .AND. TRIM(HR%FFIELD_NAME) == TRIM(MESH_NAME(II))) THEN
                                           I_NEW_FFIELD = II
                                           EXIT II_LOOP2
                                        END IF
                                     END DO II_LOOP2
                                  END IF
                               END IF
                            END IF
                            HR%I_FFIELD = I_New_FField
                            HR%FFIELD_NAME = TRIM(MESH_NAME(I_New_FField))
                         END IF
                         COLOR_INDEX = HR%COLOR_INDEX
                         IF (COLOR_METHOD == 4) THEN
                            COLOR_INDEX = EVAC_AVATAR_NCOLOR ! default, cyan
                            IF (ABS(I_TARGET)>0 .AND. ABS(I_TARGET)<=N_DOORS ) &
                                 COLOR_INDEX = EVAC_DOORS(ABS(I_TARGET))%Avatar_Color_Index
                            IF (ABS(I_TARGET)>N_DOORS .AND. ABS(I_TARGET)<=N_DOORS+N_EXITS) &
                                 COLOR_INDEX = EVAC_EXITS(ABS(I_TARGET)-N_DOORS)%Avatar_Color_Index
                         END IF
                         IF (COLOR_METHOD == 5) COLOR_INDEX = EVAC_AVATAR_NCOLOR ! default, cyan
                         IF (COLOR_METHOD == 5 .AND. I_TARGET /= 0) COLOR_INDEX = HR%COLOR_INDEX
                         HR%I_TARGET = I_TARGET
                         HR%COLOR_INDEX = COLOR_INDEX
                         I_TARGET_OLD = I_TARGET  ! No Nash iterations
                      END IF
                   END IF Herding_Type_Agent

                   IF (ABS(I_TARGET_OLD) /= ABS(I_TARGET)) THEN
                      N_CHANGE_DOORS = N_CHANGE_DOORS + 1
                      IF (I_TARGET > 0) THEN  ! The door is visible, if i_target > 0, if < 0 not visible
                         IF (I_TARGET <= N_DOORS ) THEN
                            EVEL = SQRT((HR%X-0.5_EB*(EVAC_DOORS(I_TARGET)%X1+EVAC_DOORS(I_TARGET)%X2))**2 + &
                                 (HR%Y-0.5_EB*(EVAC_DOORS(I_TARGET)%Y1+EVAC_DOORS(I_TARGET)%Y2))**2)
                            EVEL = 50.0_EB*EVEL/EVAC_DOORS(I_TARGET)%R_NTARGET + 1.0_EB
                            II = MIN(50,MAX(1,INT(EVEL)))
                            EVAC_DOORS(I_TARGET)%NTARGET(II:50) = EVAC_DOORS(I_TARGET)%NTARGET(II:50) + 1
                         ELSE
                            EVEL = SQRT((HR%X-0.5_EB*(EVAC_EXITS(I_TARGET-N_DOORS)%X1+EVAC_EXITS(I_TARGET-N_DOORS)%X2))**2 + &
                                 (HR%Y-0.5_EB*(EVAC_EXITS(I_TARGET-N_DOORS)%Y1+EVAC_EXITS(I_TARGET-N_DOORS)%Y2))**2)
                            EVEL = 50.0_EB*EVEL/EVAC_EXITS(I_TARGET-N_DOORS)%R_NTARGET + 1.0_EB
                            II = MIN(50,MAX(1,INT(EVEL)))
                            EVAC_EXITS(I_TARGET-N_DOORS)%NTARGET(II:50) = EVAC_EXITS(I_TARGET-N_DOORS)%NTARGET(II:50) + 1
                         END IF
                      END IF
                      IF (I_TARGET_OLD > 0) THEN
                         IF (I_TARGET_OLD <= N_DOORS ) THEN
                            EVEL = SQRT((HR%X-0.5_EB*(EVAC_DOORS(I_TARGET_OLD)%X1+EVAC_DOORS(I_TARGET_OLD)%X2))**2+&
                                 (HR%Y-0.5_EB*(EVAC_DOORS(I_TARGET_OLD)%Y1+EVAC_DOORS(I_TARGET_OLD)%Y2))**2)
                            EVEL = 50.0_EB*EVEL/EVAC_DOORS(I_TARGET_OLD)%R_NTARGET + 1.0_EB
                            II = MIN(50,MAX(1,INT(EVEL)))
                            EVAC_DOORS(I_TARGET_OLD)%NTARGET(II:50) = EVAC_DOORS(I_TARGET_OLD)%NTARGET(II:50) - 1
                         ELSE
                            EVEL = SQRT( &
                                 (HR%X-0.5_EB*(EVAC_EXITS(I_TARGET_OLD-N_DOORS)%X1 + EVAC_EXITS(I_TARGET_OLD-N_DOORS)%X2))**2 + &
                                 (HR%Y-0.5_EB*(EVAC_EXITS(I_TARGET_OLD-N_DOORS)%Y1 + EVAC_EXITS(I_TARGET_OLD-N_DOORS)%Y2))**2 )
                            EVEL = 50.0_EB*EVEL/EVAC_EXITS(I_TARGET_OLD-N_DOORS)%R_NTARGET + 1.0_EB
                            II = MIN(50,MAX(1,INT(EVEL)))
                            EVAC_EXITS(I_TARGET_OLD-N_DOORS)%NTARGET(II:50) = &
                                 EVAC_EXITS(I_TARGET_OLD-N_DOORS)%NTARGET(II:50) - 1
                         END IF
                      END IF
                   END IF
                ELSE ! Do not update door flow field at this time step
                   IF (J > 0) GROUP_KNOWN_DOORS(J)%I_TARGET = HR%I_TARGET
                   IF (COLOR_METHOD == 5 .AND. J > 0) COLOR_TMP(J) = HR%COLOR_INDEX
                   IF (COLOR_METHOD == 4 .AND. J > 0) COLOR_TMP(J) = HR%COLOR_INDEX
                   IF (J > 0) GROUP_LIST(J)%GROUP_I_FFIELDS(I_EGRID) = HR%I_FFIELD
                END IF  ! Update door flow field: is RN large enough?
             ELSE
                ! This group has already tried to change the field.
                IF (COLOR_METHOD == 5 .AND. J > 0) HR%COLOR_INDEX = COLOR_TMP(J)
                IF (COLOR_METHOD == 4 .AND. J > 0) HR%COLOR_INDEX = COLOR_TMP(J)
                HR%I_FFIELD    = GROUP_LIST(J)%GROUP_I_FFIELDS(I_EGRID)
                HR%FFIELD_NAME = TRIM(MESH_NAME(HR%I_FFIELD))
                IF (J > 0) HR%I_TARGET = GROUP_KNOWN_DOORS(J)%I_TARGET
             END IF  ! GROUP_I_FIELD=0
          END DO CHANGE_DOOR_LOOP  ! Loop over agents
       END IF  ! T > 0

       ! ========================================================
       ! Move loop: The first step of the VV algorithm
       ! ========================================================
       ! ========================================================
       ! Step 1: Here the humans are moved to the new positions using the present velocities.
       ! After the first step the forces at the new postitions are calculated and the velocities are updated.
       ! (The 'dissipative' self-driving force contribution to the velocities is updated self-consistently,
       ! but the other terms are not.)
       ! ========================================================
       D_HUMANS_MIN = HUGE(D_HUMANS_MIN)
       D_WALLS_MIN  = HUGE(D_WALLS_MIN)
       IF (C_HAWK >= 0.0_EB) N_HawkDoveCount(:,I_EGRID) = 0

       EVAC_MOVE_LOOP: DO I=1, N_HUMANS
          HR=>HUMAN(I)
          IF (HR%IEL > 0) THEN
             ! The agent is from some evac line
             GUARD_MEN = EVAC_EVACS(HR%IEL)%GUARD_MEN
             TARGET_X = EVAC_EVACS(HR%IEL)%TARGET_X
             TARGET_Y = EVAC_EVACS(HR%IEL)%TARGET_Y
          ELSE
             ! The agent is from some entr line, not yet as an input
             ! GUARD_MEN = EVAC_ENTRYS(ABS(HR%IEL))%GUARD_MEN
             GUARD_MEN = 0
             TARGET_X = 0.5_EB*(XS+XF)
             TARGET_Y = 0.5_EB*(YS+YF)
          END IF
          ! HR%Z is the real z-coordinate of the agent (inclines, stairs,etc),
          ! HR_Z is z-coordinate of the main evac mesh
          HR_Z = 0.5_EB*(ZS+ZF)
          LAMBDAW = LAMBDA_WALL
          A_WALL  = FAC_A_WALL*HR%A
          B_WALL  = FAC_B_WALL*HR%B

          GAME    = NOISEME
          GATH    = NOISETH
          GACM    = NOISECM

          L_FALLEN_DOWN = .FALSE.
          IF (HR%T_FallenDown < T) THEN
             ! This agent has already fallen down
             L_FALLEN_DOWN = .TRUE.
             GATH = 0.0_EB
             A_WALL = 0.0_EB
          END IF
          HoleFallLoop: DO J = 1, N_HOLES
             EHX => EVAC_HOLES(J)
             IF (EHX%TIME_FALL_DOWN<0.0_EB .OR. EHX%IMESH/=NM) CYCLE HoleFallLoop
             IF (EHX%TIME_FALL_DOWN<T .AND. &
                  EHX%X1<=HR%X .AND. EHX%X2>=HR%X .AND. EHX%Y1<=HR%Y .AND. EHX%Y2>=HR%Y) THEN
                L_FALLEN_DOWN = .TRUE.
                N_DEAD = N_DEAD+1
                HR%T_FallenDown = T
                GATH = 0.0_EB
                A_WALL = 0.0_EB
                EHX%TIME_FALL_DOWN = -1.0_EB
                CYCLE HoleFallLoop
             END IF
          END DO HoleFallLoop
          L_DEAD  = .FALSE.
          IF ( HR%INTDOSE >= 1.0_EB  ) THEN
             L_DEAD = .TRUE.
             ! No random force for a dead person.
             GATH = 0.0_EB
             ! No psychological force terms for a dead person.
             A_WALL = 0.0_EB
             IF (HR%TPRE /= HUGE(HR%TPRE)) THEN
                IF (L_FALLEN_DOWN) N_DEAD = N_DEAD-1 ! Do not duoble count
                N_DEAD = N_DEAD+1
                WRITE (LU_EVACOUT,FMT='(A,I6,A,F8.2,A,I6)') ' EVAC: Agent n:o', &
                     HR%ILABEL, ' dead at ', T, ' s, number of casualties ', N_DEAD
             END IF
             IF (.NOT. L_FALLEN_DOWN) THEN
                HR%T_FallenDown = T
                L_FALLEN_DOWN = .TRUE.
             END IF
             HR%TPRE = HUGE(HR%TPRE)
             HR%TDET = HUGE(HR%TDET)
             HR%COLOR_INDEX = EVAC_AVATAR_NCOLOR
          ELSE
             FED_MAX_ALIVE = MAX(FED_MAX_ALIVE,HR%INTDOSE)
          END IF
          IF (L_FALLEN_DOWN) HR%COLOR_INDEX = EVAC_AVATAR_NCOLOR
          FED_MAX = MAX(FED_MAX,HR%INTDOSE)  ! Dead or alive
          HR_TAU      = HR%TAU
          HR_TAU_INER = HR%TAU_INER
          ! MAX(0,HR%GROUP_ID)  Group index
          IF (MAX(0,HR%GROUP_ID)==0) THEN
             TPRE = HR%TPRE + HR%TDET ! Lonely agent
          ELSE
             TPRE = HR%TDET ! Member of a group
          END IF

          !Issue1547: Nervousness related parameter is taken from the agent's person type.
          !Issue1547: Later just the MAX_V0_FAC can be used inside this EVAC_MOVE_LOOP.
          HR_Speed = HR%SPEED
          MAX_V0_FAC = EVAC_PERSON_CLASSES(HR%IPC)%MAXIMUM_V0_FACTOR
          IF (MAX_V0_FAC>0.0_EB .AND. T > TPRE) THEN
             !Issue1547: Nervousness related parameter is taken from the agent's person type.
             !Issue1547: Later just the MAX_V0_FAC can be used inside this EVAC_FORCE_LOOP.
             !Issue1547: Nervousness related calculation: average speed is saved as HR%Speed_ave= vi_ave/vi_0
             !Issue1547: HR%Speed_ave is almost the nervousness parameter (1 - vi_ave/vi_0)
             Eta_Nervous = MIN(MAX(1.0_EB -  HR%Speed_ave, 0.0_EB),1.0_EB) ! To be sure that 0 <= Eta <= 1
             Speed_Nervous = (1.0_EB-Eta_Nervous)*HR%Speed + Eta_Nervous*MAX_V0_FAC
             HR_Speed = Speed_Nervous
          END IF

          IF (T < T_BEGIN .OR. T < TPRE) THEN
             HR_TAU = MAX(CF_MIN_TAU, 0.1_EB*HR%TAU)
             HR_TAU_INER = MAX(CF_MIN_TAU_INER, 0.1_EB*HR%TAU_INER)
          END IF

          ! Counterflow: Increase motivation to go ahead
          IF (HR%COMMITMENT > 0.01_EB .AND. T > T_BEGIN) THEN
             HR_TAU      = MAX(CF_MIN_TAU, HR%COMMITMENT*CF_FAC_TAUS*HR_TAU + (1.0_EB-HR%COMMITMENT)*HR_TAU)
             HR_TAU_INER = MAX(CF_MIN_TAU_INER, HR%COMMITMENT*CF_FAC_TAUS*HR_TAU_INER + (1.0_EB-HR%COMMITMENT)*HR_TAU_INER)
          END IF
          IF (C_HAWK >= 0.0_EB)THEN
             IF (HR%A <= 500.0_EB) THEN
                ! Hawk - dove game: more random noise for hawks
                GATH = GATH*THETA_HAWK
                ! DoseCrit1 = hr%a   modified, no sharp transition, HRA_dove  => HRA_hawk in one second
                IF (HR%I_DoorAlgo==1) HR%DoseCrit1 = MAX(HR%A, HR%DoseCrit1-9.0_EB*HR%A*DTSP)
                IF (HR%I_DoorAlgo==1) N_HawkDoveCount(1,I_EGRID) = N_HawkDoveCount(1,I_EGRID) + 1
                IF (HR%I_DoorAlgo==2) N_HawkDoveCount(3,I_EGRID) = N_HawkDoveCount(3,I_EGRID) + 1
             ELSE
                ! DoseCrit1 = hr%a   modified, no sharp transition, HRA_hawk  => HRA_dove in one second
                IF (HR%I_DoorAlgo==1) HR%DoseCrit1 = MIN(HR%A, HR%DoseCrit1+0.9_EB*HR%A*DTSP)
                IF (HR%I_DoorAlgo==1) N_HawkDoveCount(2,I_EGRID) = N_HawkDoveCount(2,I_EGRID) + 1
                IF (HR%I_DoorAlgo==2) N_HawkDoveCount(4,I_EGRID) = N_HawkDoveCount(4,I_EGRID) + 1
             END IF
          END IF
          !
          ! In which grid cell is the agent?
          !
          XI = CELLSI(FLOOR((HR%X-XS)*RDXINT))
          YJ = CELLSJ(FLOOR((HR%Y-YS)*RDYINT))
          ZK = CELLSK(FLOOR((HR_Z-ZS)*RDZINT))
          II = FLOOR(XI+1.0_EB)
          JJ = FLOOR(YJ+1.0_EB)
          KK = 1
          IIX = FLOOR(XI+0.5_EB)
          JJY = FLOOR(YJ+0.5_EB)
          KKZ = FLOOR(ZK+0.5_EB)
          !
          HR%W = 0.0_EB
          SMOKE_SPEED_FAC = 1.0_EB
          IF (.NOT.L_DEAD .AND. T >= HR%TDET) HR%DETECT1 = IBSET(HR%DETECT1,0)  ! Detected by the T_det distribution, bit 0
          Not_Init_Phase_Use_FED: IF (T > T_BEGIN .AND. ICYC > 1 .AND. USE_FED) THEN
             ! Calculate purser's fractional effective dose (FED)
             ! Note: Purser uses minutes, here DT is in seconds
             ! FED_DOSE = FED_LCO*FED_VCO2 + FED_LO
             IF (HR%IEL > 0) THEN  ! From an evac line
                IF (T >= EVAC_EVACS(HR%IEL)%T_START_FED) THEN
                   HR%INTDOSE =DTSP*HUMAN_GRID(II,JJ)%FED_CO_CO2_O2 + HR%INTDOSE
                END IF
             ELSE ! From an entr line
                HR%INTDOSE = DTSP*HUMAN_GRID(II,JJ)%FED_CO_CO2_O2 + HR%INTDOSE
             END IF
             ! Smoke density vs speed
             ! Lund 2003, report 3126 (Frantzich & Nilsson)
             ! V0(K) = V0*( 1 + (BETA*K)/ALPHA ), Where [K]=1/M EXT.COEFF
             ! BETA=-0.057 m2/s, SIGMA=0.015, [-0.087,-0.028]
             ! ALPHA=0.706 m/s, SIGMA=0.069, [0.565,0.847]
             ! [SOOT_DENS] = mg/m3
             ! [K] = 1/m
             ! [MASS_EXTINCTION_COEFFICIENT] = m2/kg
             ! K = MASS_EXTINCTION_COEFFICIENT*SOOT_DENS*1.0E-6
             ! VISIBILITY = 3/K  [m]
             SMOKE_BETA  = SMOKE_SPEED_BETA
             SMOKE_ALPHA = SMOKE_SPEED_ALPHA
             IF (SMOKE_MIN_SPEED_FACTOR > -99998.0_EB) THEN
                SMOKE_MIN_SPEED_TMP = SMOKE_MIN_SPEED_FACTOR
             ELSE
                SMOKE_MIN_SPEED_TMP = SMOKE_MIN_SPEED/HR%SPEED
             END IF

             SELECT CASE(SMOKE_KS_SPEED_FUNCTION)
             CASE (1)
                ! Linear smoke density vs speed function, c(Ks) = 1 + beta*Ks/alpha
                SMOKE_SPEED_FAC = 1.0_EB + (SMOKE_BETA/SMOKE_ALPHA)* &
                     EVAC_MASS_EXTINCTION_COEFF*1.0E-6_EB*HUMAN_GRID(II,JJ)%SOOT_DENS
             CASE DEFAULT
                WRITE(LU_ERR,'(A,I3)') 'ERROR: PERS line has a non-valid SMOKE_KS_SPEED_FUNCTION,',SMOKE_KS_SPEED_FUNCTION
                WRITE(LU_EVACOUT,'(A,I3)') 'ERROR: PERS line has a non-valid SMOKE_KS_SPEED_FUNCTION,',SMOKE_KS_SPEED_FUNCTION
                STOP_STATUS = EVACUATION_STOP ; RETURN
             END SELECT
             ! Next if is just an obsolote feature (Poland request), it is not used
             ! if default SMOKE_MIN_SPEED_VISIBILITY is given.
             IF (EVAC_MASS_EXTINCTION_COEFF*1.0E-6_EB*HUMAN_GRID(II,JJ)%SOOT_DENS > &
                  3.0_EB/SMOKE_MIN_SPEED_VISIBILITY) THEN  ! 1/m > 300 1/m
                SMOKE_SPEED_FAC = MIN(SMOKE_SPEED_FAC, SMOKE_SPEED_FAC*( 2.0_EB - &
                     (EVAC_MASS_EXTINCTION_COEFF*1.0E-6_EB*HUMAN_GRID(II,JJ)%SOOT_DENS - &
                     3.0_EB/SMOKE_MIN_SPEED_VISIBILITY ) / (3.0_EB/SMOKE_MIN_SPEED_VISIBILITY) ) )
             END IF
             SMOKE_SPEED_FAC = MAX(SMOKE_SPEED_FAC, SMOKE_MIN_SPEED_TMP)

             ! Check the smoke density for the detection by smoke
             IF (.NOT.L_DEAD .AND. HUMAN_GRID(II,JJ)%SOOT_DENS > TDET_SMOKE_DENS) THEN
                HR%TDET = MIN(HR%TDET,T)
                HR%DETECT1 = IBSET(HR%DETECT1,1)  ! Detected by the local smoke consentration, bit 1
             END IF

             ! Check the evacuation devices for the detection
             IF (L_FIRST_PASS) THEN
                EDEV_LOOP: DO J = 1, N_EDEV
                   EDV => EVAC_EDEV(J)
                   INPUTS_LOOP: DO N = 1, EDV%N_INPUTS
                      KK = EDV%INPUT_DEVC_INDEX(N)
                      IF (EVAC_DEVICES(KK)%T_Change <= T .AND. EVAC_DEVICES(KK)%CURRENT .AND. &
                           EVAC_DEVICES(KK)%USE_NOW) THEN

                         IF (EDV%i_pre_dist > -1) THEN
                            CALL TPRE_GENERATION(EDV%i_pre_dist,EDV%Tpre_low,EDV%Tpre_high,EDV%Tpre_mean,EDV%Tpre_para, &
                                 EDV%Tpre_para2,TPRE)
                            IF (STOP_STATUS>NO_STOP) RETURN
                         ELSE
                            TPRE = HR%TPRE ! Use evac/pers namelist reaction times
                         END IF
                         HR%DETECT1 = IBSET(HR%DETECT1,2)  ! Detected by some device, bit 2
                         IF (T+TPRE < HR%TDET+HR%TPRE) THEN
                            WRITE(LU_EVACOUT,FMT='(A,I6,A,A,A,F8.2,A,F6.2,A)') ' Agent n:o ', HR%ILABEL, ' detection by ', &
                                 TRIM(EVAC_DEVICES(KK)%DEVC_ID), ' at', T, ' s, new Tpre is ', TPRE, ' s'
                            HR%TDET = T
                            HR%TPRE = TPRE
                         ELSE
                            WRITE(LU_EVACOUT,FMT='(A,I6,A,A,A,F8.2,A)') ' Agent n:o ', HR%ILABEL, ' detection by ', &
                                 TRIM(EVAC_DEVICES(KK)%DEVC_ID), ' at', T, ' s'
                         END IF
                      END IF
                   END DO INPUTS_LOOP
                END DO EDEV_LOOP
             END IF
          END IF Not_Init_Phase_Use_FED
          HR%V0_FAC = SMOKE_SPEED_FAC

          ! ========================================================
          ! Calculate persons prefered walking direction V0
          ! ========================================================
          N = NM_STRS_INDEX
          IF (TAU_CHANGE_V0 > 1.0E-12_EB) THEN
             ! Stairs mesh: next subroutine call is needed to set up the targets etc correctly,
             ! but UBAR,VBAR are not needed.
             IF (NM_STRS_MESH .OR. HR%CROWBAR_UPDATE_V0 .OR. HR%CROWBAR_READ_IN) &
                  CALL FIND_PREFERRED_DIRECTION(I, N, T,T_BEGIN, L_DEAD, NM_STRS_MESH, &
                  II, JJ, XI, YJ, UBAR, VBAR, HR_TAU, TPRE, NM, I_STRS_DOOR, HR_SPEED)
             IF (STOP_STATUS>NO_STOP) RETURN
             ! Crowbar (read in camera data) needs tau and speed updates
             ! Collision avoidance (incl. counterflow), do not update v0 on every time step.
             UBAR = HR%UBAR; VBAR = HR%VBAR
          ELSE
             ! Update v0 on every time step, no collision avoidance.
             IF(GUARD_MEN == 1) THEN
!!$                CALL GUARD(NM, XI, YJ, WSPA, WSPB, TARGET_X, TARGET_Y)
                UBAR = WSPA
                VBAR = WSPB
                EVEL = SQRT(UBAR**2 + VBAR**2)
                UBAR = UBAR/EVEL ; VBAR = VBAR/EVEL
             ELSE
                CALL FIND_PREFERRED_DIRECTION(I, N, T, T_BEGIN, L_DEAD, NM_STRS_MESH, &
                     II, JJ, XI, YJ, UBAR, VBAR, HR_TAU, TPRE, NM, I_STRS_DOOR, HR_SPEED)
                IF (STOP_STATUS>NO_STOP) RETURN
             END IF
          END IF
          ! ========================================================
          ! Prefered walking direction V0 is now (UBAR,VBAR)
          ! ========================================================

          ! Inclines: Velocities are along the incline
          !           coordinates are projected on the (x,y) plane
          COS_X = 1.0_EB
          COS_Y = 1.0_EB
          SPEED_XM = HR_SPEED
          SPEED_XP = HR_SPEED
          SPEED_YM = HR_SPEED
          SPEED_YP = HR_SPEED
          IF (NM_STRS_MESH) THEN
             STRS_INDX = N
             STRP=>EVAC_STRS(N)
          END IF
          !
          ! Check if an agent is on a spectator stand.
          IF (.NOT. NM_STRS_MESH) HR%Z = 0.5_EB*(ZS+ZF)  ! The agent is not on any incline
          SS_LOOP1: DO J = 1, N_SSTANDS
             ESS => EVAC_SSTANDS(J)
             IF (ESS%IMESH /= NM) CYCLE SS_LOOP1
             IF (ESS%X1 > HR%X) CYCLE SS_LOOP1
             IF (ESS%X2 < HR%X) CYCLE SS_LOOP1
             IF (ESS%Y1 > HR%Y) CYCLE SS_LOOP1
             IF (ESS%Y2 < HR%Y) CYCLE SS_LOOP1
             IF (ESS%USE_V0 .AND. .NOT.(HR%I_FFIELD == ESS%I_VENT_FFIELD .OR. ESS%I_VENT_FFIELD == 0)) CYCLE SS_LOOP1

             ! Next are here just for a test case.  The evacuation flow fields
             ! in front of the doors/exits are not optimal.  Close to a door it
             ! would be better to use an ad-hoc movement direction, something
             ! like is done in the STRS staircase model.
             IF (.NOT.L_DEAD .AND. ESS%USE_V0 .AND. TAU_CHANGE_V0 <= 1.0E-12_EB) THEN
                ! Update v0 on every time step, if no collision avoidance.
                EVEL = SQRT(ESS%UBAR0**2 + ESS%VBAR0**2)
                IF (HR%I_FFIELD == ESS%I_VENT_FFIELD .OR. ESS%I_VENT_FFIELD == 0 .AND. EVEL > 1.0E-12_EB) THEN
                   UBAR = ESS%UBAR0/EVEL
                   VBAR = ESS%VBAR0/EVEL

                END IF
             END IF
             COS_X = ESS%COS_X
             COS_Y = ESS%COS_Y
             FAC_V0_UP   = ESS%FAC_V0_UP
             FAC_V0_DOWN = ESS%FAC_V0_DOWN
             FAC_V0_HORI = ESS%FAC_V0_HORI
             IF (EVAC_PERSON_CLASSES(HR%IPC)%FAC_V0_HORI >= TWO_EPSILON_EB) THEN
                FAC_V0_HORI = EVAC_PERSON_CLASSES(HR%IPC)%FAC_V0_HORI
             END IF
             IF (EVAC_PERSON_CLASSES(HR%IPC)%FAC_V0_UP >= TWO_EPSILON_EB) THEN
                IF ((ESS%H - ESS%H0) <= -TWO_EPSILON_EB) THEN
                   FAC_V0_DOWN = EVAC_PERSON_CLASSES(HR%IPC)%FAC_V0_UP
                ELSE
                   FAC_V0_UP = EVAC_PERSON_CLASSES(HR%IPC)%FAC_V0_UP
                END IF
             END IF
             IF (EVAC_PERSON_CLASSES(HR%IPC)%FAC_V0_DOWN >= TWO_EPSILON_EB) THEN
                IF ((ESS%H - ESS%H0) <= -TWO_EPSILON_EB) THEN
                   FAC_V0_UP = EVAC_PERSON_CLASSES(HR%IPC)%FAC_V0_DOWN
                ELSE
                   FAC_V0_DOWN = EVAC_PERSON_CLASSES(HR%IPC)%FAC_V0_DOWN
                END IF
             END IF
             SELECT CASE (ESS%IOR)
             CASE(-1)
                SPEED_XM = COS_X*HR_SPEED*FAC_V0_DOWN
                SPEED_XP = COS_X*HR_SPEED*FAC_V0_UP
                SPEED_YM = HR_SPEED*FAC_V0_HORI
                SPEED_YP = HR_SPEED*FAC_V0_HORI
                HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%X1-HR%X)/ABS(ESS%X1-ESS%X2)
             CASE(+1)
                SPEED_XM = COS_X*HR_SPEED*FAC_V0_UP
                SPEED_XP = COS_X*HR_SPEED*FAC_V0_DOWN
                SPEED_YM = HR_SPEED*FAC_V0_HORI
                SPEED_YP = HR_SPEED*FAC_V0_HORI
                HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%X2-HR%X)/ABS(ESS%X1-ESS%X2)
             CASE(-2)
                SPEED_XM = HR_SPEED*FAC_V0_HORI
                SPEED_XP = HR_SPEED*FAC_V0_HORI
                SPEED_YM = COS_Y*HR_SPEED*FAC_V0_DOWN
                SPEED_YP = COS_Y*HR_SPEED*FAC_V0_UP
                HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%Y1-HR%Y)/ABS(ESS%Y1-ESS%Y2)
             CASE(+2)
                SPEED_XM = HR_SPEED*FAC_V0_HORI
                SPEED_XP = HR_SPEED*FAC_V0_HORI
                SPEED_YM = COS_Y*HR_SPEED*FAC_V0_UP
                SPEED_YP = COS_Y*HR_SPEED*FAC_V0_DOWN
                HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%Y2-HR%Y)/ABS(ESS%Y1-ESS%Y2)
             END SELECT
          END DO SS_LOOP1

          ! Set height and speed for an agent in stairs
          IF (NM_STRS_MESH) THEN
             CALL GETSTAIRSPEEDANDZ(SPEED_XM, SPEED_XP, SPEED_YM, SPEED_YP, HR_SPEED, STRP, HR)
             IF (STOP_STATUS>NO_STOP) RETURN
          END IF

          ! ========================================================
          ! Step 2:  The new velocities are calculated using the old forces.
          ! ========================================================
          ! (F_X,F_Y): Wall forces + agent-agent forces (social + contact forces)
          ! Random forces and self-driving force are treated separately.
          U_NEW = HR%U + 0.5_EB*HR%F_X*DTSP/HR%MASS
          V_NEW = HR%V + 0.5_EB*HR%F_Y*DTSP/HR%MASS
          ! Rotational motion:
          OMEGA_NEW = HR%OMEGA + 0.5_EB*DTSP*HR%TORQUE/HR%M_INER

          ! Add self-driving force and torque
          EVEL = SQRT(UBAR**2 + VBAR**2)
          IF (EVEL > 0.001_EB) THEN
             !Inclines: U,V are the (x,y) plane projection velocities

             SPEED_X = SPEED_XP*(0.5_EB + SIGN(0.5_EB,UBAR)) + SPEED_XM*(0.5_EB - SIGN(0.5_EB,UBAR))
             SPEED_X = SPEED_X*HR%V0_FAC*(UBAR/EVEL)
             SPEED_Y = SPEED_YP*(0.5_EB + SIGN(0.5_EB,VBAR)) + SPEED_YM*(0.5_EB - SIGN(0.5_EB,VBAR))
             SPEED_Y = SPEED_Y*HR%V0_FAC*(VBAR/EVEL)
             TAU_FAC = MAX(0.1_EB,SQRT(SPEED_X**2 + SPEED_Y**2)/MAX(0.1_EB,HR_SPEED))
             HR_TAU  = HR_TAU*TAU_FAC
             U_NEW = U_NEW + 0.5_EB*(DTSP/HR_TAU)*(SPEED_X - HR%U)
             HR%UBAR = UBAR
             V_NEW = V_NEW + 0.5_EB*(DTSP/HR_TAU)*(SPEED_Y - HR%V)
             HR%VBAR = VBAR

             !Issue1547: Nervousness related calculation: average speed as vi_ave/vi_0
             !Issue1547: HR%Speed_ave is almost the nervousness parameter (1 - vi_ave/vi_0)
             !Issue1547: The vi_0 above is the current desired velocity of an agent.
             !Issue1547: The Speed_ave is declared in "type.f90" as a compenent in the human_type.
             !Issue1547: Pointer HR => HUMAN(I), I is the evac_move_loop agent index.
             !Issue1547: The agent movement algorithm consist of two loops over all agents
             !Issue1547: at a mesh, but it is enough to calculate the nervousness value here
             !Issue1547: at the move loop only. The later loop is the force calculation loop.
             !Issue1547: But the change of the current v0_i for an agent should be done in
             !Issue1547: both loops.  Note that the agent may have different desired velocity
             !Issue1547: due to: Inclines (HR%Speed is along the inclines, HR%U,HR%V are plane
             !Issue1547: projected velocities, FAC_V0_DOWN etc also may be given), Smoke (smoke
             !Issue1547: reduces desired velocity)
             IF (MAX_V0_FAC>0.0_EB .AND. T > TPRE) THEN
                EVEL2 = SQRT(SPEED_X**2 + SPEED_Y**2) ! Initial desired velocity projected to (x,y) plane
                HR%Speed_ave = (HR%T_Speed_ave*HR%Speed_ave + DTSP*(SQRT(HR%U**2+HR%V**2)/EVEL2)) / &
                     (HR%T_Speed_ave+DTSP)
                HR%T_Speed_ave = HR%T_Speed_ave + DTSP
             END IF

             IF (VBAR >= 0.0_EB) THEN
                ANGLE = ACOS(UBAR/EVEL)
             ELSE
                ANGLE = 2.0_EB*PI - ACOS(UBAR/EVEL)
             END IF

             ! Collison avoidance
             ANGLE = ANGLE + HR%ANGLE_OLD
             DO WHILE (ANGLE >= 2.0_EB*PI)
                ANGLE = ANGLE - 2.0_EB*PI
             END DO
             DO WHILE (ANGLE < 0.0_EB)
                ANGLE = ANGLE + 2.0_EB*PI
             END DO

             IF (ANGLE == 2.0_EB*PI) ANGLE = 0.0_EB  ! Angle is [0,2pi)

             ! Rotational motion: J(DW/DT) = (J/T_INER)*( ((ANGLE-ANGLE_0/PI))W_0 - W )
             IF (ABS(ANGLE-HR%ANGLE) <= PI ) THEN
                ! Zero is not crossed.
                OMEGA_NEW = OMEGA_NEW + 0.5_EB*(DTSP/HR_TAU_INER)*( (ANGLE-HR%ANGLE)*(OMEGA_0/PI) - HR%OMEGA)
             ELSE
                ! Zero is crossed
                OMEGA_NEW = OMEGA_NEW + 0.5_EB*(DTSP/HR_TAU_INER)* &
                     ( (2.0_EB*PI-ABS(ANGLE-HR%ANGLE))*SIGN(1.0_EB , HR%ANGLE-ANGLE)*(OMEGA_0/PI) - HR%OMEGA)
             END IF
          ELSE  ! No V0
             U_NEW = U_NEW + 0.5_EB*(DTSP/HR_TAU)* (- HR%U)
             V_NEW = V_NEW + 0.5_EB*(DTSP/HR_TAU)* (- HR%V)
             ! Slow rotation down if no direction available, i.e., target OMEGA_0 is zero.
             OMEGA_NEW = OMEGA_NEW + 0.5_EB*(DTSP/HR_TAU_INER)*(-HR%OMEGA)
             HR%UBAR = 0.0_EB
             HR%VBAR = 0.0_EB
          END IF
          !
          ! check, if a new random force is needed on next time step.
          ! Poisson distribution, i.e., the agents do not have memory.
          ! P[NO CHANGE DURING DT] = EXP(-DTSP/HR%TAU)
          IF (GATH >= TWO_EPSILON_EB .AND. T > T_BEGIN) THEN
             CALL RANDOM_NUMBER(RN_REAL)
             RN = REAL(RN_REAL,EB)
             IF (RN > EXP(-DTSP/(0.2_EB*HR_TAU))) HR%NEWRND = .TRUE.
             IF (HR%NEWRND) THEN
                HR%NEWRND = .FALSE.
                IF (EVEL<0.001_EB) GATH=GATH/100._EB  ! No prefered direction, stand still
                GATH = GATH*(HR_SPEED/1.3_EB)**2
                ! GATH is found to be more or less nice for speeds about 1.3 m/s
                HR%KSI = (GAUSSRAND(GAME, GATH, GACM))
                CALL RANDOM_NUMBER(RN_REAL)
                RN = REAL(RN_REAL,EB)
                HR%ETA = 2.0_EB*PI*RN
                EVEL = SQRT(HR%U**2 + HR%V**2)
                ! Random noice variance: GATH = MAX( GATH, GATH*(100.0_EB-110.0_EB*ABS(EVEL/HR%SPEED)) )
                ! Above GATH is found to be more or less nice for speeds about 1.3 m/s
                ! Scale the random force by the current target speed V0 (Note: SQRT(GATH) = STD.DEV.)
                ! NOTE: Speed reduction due to smoke is treated elsewhere.
                HR%KSI = ABS(HR%KSI)
                HR%KSI = MAX(HR%KSI, HR%KSI*10.0_EB*(1.0_EB-MIN(0.9_EB,EVEL/HR_SPEED)))
             END IF
          END IF

          ! Add random force term
          IF (GATH >= TWO_EPSILON_EB .AND. T > T_BEGIN) THEN
             U_NEW = U_NEW + 0.5_EB*DTSP*HR%V0_FAC*HR%MASS*HR%KSI*COS(HR%ETA)/HR%MASS
             V_NEW = V_NEW + 0.5_EB*DTSP*HR%V0_FAC*HR%MASS*HR%KSI*SIN(HR%ETA)/HR%MASS
             OMEGA_NEW = OMEGA_NEW + 0.5_EB*DTSP* 1.0_EB*SIGN(HR%KSI,HR%ETA-PI)
          END IF

          ! Check that velocities are not too large, i.e., unphysical (less than 10 m/s for humans)
          EVEL = SQRT(U_NEW**2 + V_NEW**2)
          IF ( EVEL > VMAX_TIMO ) THEN
             U_NEW = U_NEW*(VMAX_TIMO/EVEL)
             V_NEW = V_NEW*(VMAX_TIMO/EVEL)
          END IF
          ! Check that angular velocity is not too large
          IF ( ABS(OMEGA_NEW) > OMEGA_MAX ) THEN
             OMEGA_NEW = SIGN(OMEGA_MAX,OMEGA_NEW)
          END IF
          IF (TAU_FALL_DOWN >= 0.0_EB .AND. L_FALLEN_DOWN) THEN
            T_tmp = MIN(TAU_FALL_DOWN, T - HR%T_FallenDown) ! 0 - tau_fall_down
            U_NEW = (TAU_FALL_DOWN - T_tmp)*U_NEW/TAU_FALL_DOWN
            V_NEW = (TAU_FALL_DOWN - T_tmp)*V_NEW/TAU_FALL_DOWN
            OMEGA_NEW = (TAU_FALL_DOWN - T_tmp)*OMEGA_NEW/TAU_FALL_DOWN
          END IF

          ! ========================================================
          ! Step 3:  The new coordinates are calculated using the new velocities.
          ! ========================================================
          X1 = HR%X + U_NEW*DTSP
          Y1 = HR%Y + V_NEW*DTSP
          HR%U = U_NEW
          HR%V = V_NEW

          ! The new body angle, should be on the interval [0,2pi)
          A1 = HR%ANGLE + OMEGA_NEW*DTSP
          DO WHILE (A1 >= 2.0_EB*PI)
             A1 = A1 - 2.0_EB*PI
          END DO
          DO WHILE (A1 < 0.0_EB)
             A1 = A1 + 2.0_EB*PI
          END DO
          HR%OMEGA = OMEGA_NEW
          !
          ! Check, if human is on an escalator and change the coordinates.
          SS_LOOP1B: DO J = 1, N_SSTANDS
             ESS => EVAC_SSTANDS(J)
             IF (ESS%IMESH == NM .AND. (ESS%X1 <= HR%X .AND. ESS%X2 >= HR%X) .AND. (ESS%Y1 <= HR%Y .AND. ESS%Y2 >= HR%Y) ) THEN
                COS_X = ESS%COS_X
                COS_Y = ESS%COS_Y
                SELECT CASE (ESS%IOR)
                CASE(-1)
                   X1 = X1 + COS_X*(ESS%ESC_SPEEDUP-ESS%ESC_SPEEDDN)*DTSP
                CASE(+1)
                   X1 = X1 - COS_X*(ESS%ESC_SPEEDUP-ESS%ESC_SPEEDDN)*DTSP
                CASE(-2)
                   Y1 = Y1 + COS_Y*(ESS%ESC_SPEEDUP-ESS%ESC_SPEEDDN)*DTSP
                CASE(+2)
                   Y1 = Y1 - COS_Y*(ESS%ESC_SPEEDUP-ESS%ESC_SPEEDDN)*DTSP
                END SELECT
                EXIT SS_LOOP1B
             END IF
          END DO SS_LOOP1B
          !
          ! In which grid cell is the agent, the new coordinates (T + DT)?
          XI  = CELLSI(FLOOR((X1-XS)*RDXINT))
          YJ  = CELLSJ(FLOOR((Y1-YS)*RDYINT))
          ZK  = CELLSK(FLOOR((HR_Z-ZS)*RDZINT))
          IIN = FLOOR(XI+1.0_EB)
          JJN = FLOOR(YJ+1.0_EB)
          KKN = 1
          ICN = CELL_INDEX(IIN,JJN,KKN)
          ICX = CELL_INDEX(IIN,JJ ,KKN)
          ICY = CELL_INDEX(II ,JJN,KKN)
          I_OBST  = OBST_INDEX_C(CELL_INDEX(IIN,JJN,KKN))
          I_OBSTX = OBST_INDEX_C(CELL_INDEX(IIN,JJ ,KKN))
          I_OBSTY = OBST_INDEX_C(CELL_INDEX(II ,JJN,KKN))
          HR%X_OLD = HR%X
          HR%Y_OLD = HR%Y

          ! Check, if the agent moves inside a solid object ==> might be an open
          ! vent or a 'sucking vent' used to calculate the flow fields.
          ! This is just to be fail safe.  If the user input is correct, this
          ! should never happen.
          IF ( SOLID(ICN) .AND. .NOT.OBSTRUCTION(I_OBST)%HIDDEN) THEN
             IF ( (SOLID(ICX).AND. .NOT.OBSTRUCTION(I_OBSTX)%HIDDEN) .AND. .NOT. &
                  (SOLID(ICY).AND. .NOT.OBSTRUCTION(I_OBSTY)%HIDDEN) ) THEN
                IF ( II < IIN ) THEN
                   TIM_IC = CELL_INDEX(IIN,JJN,KK)
                   CALL GET_IW(IIN,JJN,KK,-1,TIM_IW)
                   SURF_INDEX =WALL(TIM_IW)%SURF_INDEX
                   IF (SURFACE(SURF_INDEX)%VEL >= TWO_EPSILON_EB .OR. WALL(TIM_IW)%BOUNDARY_TYPE==OPEN_BOUNDARY) THEN
                      HR%X = X1
                      HR%Y = HR%Y
                      Y1 = HR%Y
                   END IF
                ELSE
                   TIM_IC = CELL_INDEX(IIN,JJN,KK)
                   CALL GET_IW(IIN,JJN,KK,+1,TIM_IW)
                   SURF_INDEX = WALL(TIM_IW)%SURF_INDEX
                   IF (SURFACE(SURF_INDEX)%VEL >= TWO_EPSILON_EB .OR. WALL(TIM_IW)%BOUNDARY_TYPE==OPEN_BOUNDARY) THEN
                      HR%X = X1
                      HR%Y = HR%Y
                      Y1 = HR%Y
                   END IF
                END IF
             ELSE IF ( (SOLID(ICY).AND. .NOT.OBSTRUCTION(I_OBSTY)%HIDDEN) .AND. .NOT. &
                  (SOLID(ICX).AND. .NOT.OBSTRUCTION(I_OBSTX)%HIDDEN) ) THEN
                IF ( JJ < JJN ) THEN
                   TIM_IC = CELL_INDEX(IIN,JJN,KK)
                   CALL GET_IW(IIN,JJN,KK,-2,TIM_IW)
                   SURF_INDEX = WALL(TIM_IW)%SURF_INDEX
                   IF (SURFACE(SURF_INDEX)%VEL > TWO_EPSILON_EB .OR. WALL(TIM_IW)%BOUNDARY_TYPE==OPEN_BOUNDARY) THEN
                      HR%X = HR%X
                      X1 = HR%X
                      HR%Y = Y1
                   END IF
                ELSE
                   TIM_IC = CELL_INDEX(IIN,JJN,KK)
                   CALL GET_IW(IIN,JJN,KK,+2,TIM_IW)
                   SURF_INDEX = WALL(TIM_IW)%SURF_INDEX
                   IF (SURFACE(SURF_INDEX)%VEL > TWO_EPSILON_EB .OR. WALL(TIM_IW)%BOUNDARY_TYPE==OPEN_BOUNDARY) THEN
                      HR%X = HR%X
                      X1 = HR%X
                      HR%Y = Y1
                   END IF
                END IF
             ELSE
                WRITE(MESSAGE,'(A,I4,A,2F8.2)') 'ERROR: EVACUATE_HUMANS, Solid icx and icy, mesh ', NM, ' pos ',X1,Y1
             END IF
          ELSE
             ! Target cell is not a solid ==> move
             HR%X = X1
             HR%Y = Y1
          END IF
          HR%X = X1
          HR%Y = Y1
          HR%ANGLE = A1
          !
          IF (HR%IEL < 0) THEN ! From an entry line
             PNX => EVAC_ENTRYS(ABS(HR%IEL))
             IF (PNX%CROWBAR_READ_IN .AND. PNX%IMESH == NM) THEN
                IF (Is_XY_Within_Bounds(HR%X,HR%Y,PNX%X1,PNX%X2,PNX%Y1,PNX%Y2,0._EB,0._EB)) THEN
                   PNX%CB_N_Agents = PNX%CB_N_Agents + 1
                   IF (HR%CROWBAR_READ_IN) PNX%CB_N_CrowbarAgents = PNX%CB_N_CrowbarAgents + 1
                END IF
             END IF
          END IF

       END DO EVAC_MOVE_LOOP
       ! ========================================================
       ! Move loop ends here
       ! ========================================================

       ! ========================================================
       ! Check if persons are leaving this mesh via doors/exits and put these persons to the target.
       ! ========================================================
       IF (N_HUMANS > 0) CALL CHECK_EXITS(T,NM)
       IF (STOP_STATUS>NO_STOP) RETURN
       IF (N_HUMANS > 0) CALL CHECK_DOORS(T,NM)
       IF (STOP_STATUS>NO_STOP) RETURN

       ! ========================================================
       ! Check if persons are entering this mesh via a corr.
       ! ========================================================
       CALL CHECK_CORRS(T,NM,DTSP)
       IF (STOP_STATUS>NO_STOP) RETURN

       ! ========================================================
       ! Add persons from entrys (specified flow rates)
       ! ========================================================
       DO I = 1, N_ENTRYS
          IF (T > T_BEGIN ) THEN
             CALL ENTRY_HUMAN(I,T,NM,ISTAT) ; IF (STOP_STATUS>NO_STOP) RETURN
          END IF
          CALL ENTRY_CROWBAR_HUMAN(I,T,NM,ISTAT) ; IF (STOP_STATUS>NO_STOP) RETURN
       END DO

       ! ========================================================
       ! Remove out-of-bounds persons (outside the grid)
       ! ========================================================
       IF (N_HUMANS > 0) CALL REMOVE_OUT_OF_GRIDS
       IF (STOP_STATUS>NO_STOP) RETURN

       IF ( ICYC >= 0) THEN
          DTSP_NEW = EVAC_DT_STEADY_STATE
       ELSE
          DTSP_NEW = EVAC_DT_FLOWFIELD  ! Initialization phase
       END IF
       DTSP_NEW = MAX(DTSP_NEW, EVAC_DT_MIN)
       DTSP_NEW = MIN(DTSP_NEW, EVAC_DT_MAX)

       SPEED_MAX  = 0.0_EB

       ! ================================================
       ! Prepare to calculate the new forces, initialize different variables and arrays for the step 3 of
       ! the SC-VV algorithm.
       ! ================================================

       ! ================================================
       ! Initialize group arrays for this main evac mesh.
       ! ================================================
       GROUP_LIST(:)%GROUP_SIZE  = 0
       GROUP_LIST(:)%GROUP_X = 0.0_EB
       GROUP_LIST(:)%GROUP_Y = 0.0_EB
       GROUP_LIST(:)%MAX_DIST_CENTER = 0.0_EB
       GROUP_LIST(:)%TPRE    = 0.0_EB
       GROUP_LIST(:)%TDET    = HUGE(GROUP_LIST(:)%TDET)

       DO J = 0, I33_DIM
          GROUP_LIST(J)%GROUP_I_FFIELDS(I_EGRID) = 0
       END DO
       DO I = 1, N_HUMANS
          HR=>HUMAN(I)
          J = MAX(0,HR%GROUP_ID)
          GROUP_LIST(J)%GROUP_SIZE = GROUP_LIST(J)%GROUP_SIZE + 1
          GROUP_LIST(J)%GROUP_X    = GROUP_LIST(J)%GROUP_X + HR%X
          GROUP_LIST(J)%GROUP_Y    = GROUP_LIST(J)%GROUP_Y + HR%Y
          GROUP_LIST(J)%GROUP_I_FFIELDS(I_EGRID) = HR%I_FFIELD
          GROUP_LIST(J)%TPRE       = MAX(GROUP_LIST(J)%TPRE,HR%TPRE)
          GROUP_LIST(J)%TDET       = MIN(GROUP_LIST(J)%TDET,HR%TDET)
       END DO
       !
       GROUP_LIST(1:)%GROUP_X = GROUP_LIST(1:)%GROUP_X / MAX(1,GROUP_LIST(1:)%GROUP_SIZE)
       GROUP_LIST(1:)%GROUP_Y = GROUP_LIST(1:)%GROUP_Y / MAX(1,GROUP_LIST(1:)%GROUP_SIZE)

       ! BLOCK_GRID_N(I,J): How many agents are in this grid cell
       BLOCK_GRID_N = 0
       GROUP_LIST(:)%MAX_DIST_CENTER = 0.0_EB
       DO I = 1, N_HUMANS
          HR => HUMAN(I)
          ! Which cell, new coordinates:
          IIN = FLOOR( CELLSI(FLOOR((HR%X-XS)*RDXINT)) + 1.0_EB)
          JJN = FLOOR( CELLSJ(FLOOR((HR%Y-YS)*RDYINT)) + 1.0_EB)
          BLOCK_GRID_N(IIN,JJN) = BLOCK_GRID_N(IIN,JJN) + 1
          J = MAX(0,HR%GROUP_ID)
          GROUP_LIST(J)%MAX_DIST_CENTER = MAX(GROUP_LIST(J)%MAX_DIST_CENTER, &
               SQRT((HR%X - GROUP_LIST(J)%GROUP_X)**2 + (HR%Y - GROUP_LIST(J)%GROUP_Y)**2))
       END DO
       MAX_HUMANS_CELL = MAX(1,MAXVAL(BLOCK_GRID_N))

       IF (N_HUMANS > 0) THEN
          GROUP_LIST(0)%GROUP_SIZE = 1
          GROUP_LIST(0)%GROUP_X    = 0.5_EB*(XS+XF)
          GROUP_LIST(0)%GROUP_Y    = 0.5_EB*(YS+YF)
          GROUP_LIST(0)%SPEED      = 1.0_EB
          GROUP_LIST(0)%INTDOSE    = 0.0_EB
          GROUP_LIST(0)%MAX_DIST_CENTER = 0.0_EB
          GROUP_LIST(0)%COMPLETE   = 1
       END IF
       ! Check if the groups are already gathered together or not.
       DO J = 1, I33_DIM
          GROUP_LIST(J)%LIMIT_COMP = RADIUS_COMPLETE_0 + RADIUS_COMPLETE_1*GROUP_LIST(J)%GROUP_SIZE
          IF ( ((GROUP_LIST(J)%MAX_DIST_CENTER <=  GROUP_LIST(J)%LIMIT_COMP) .OR. &
               (GROUP_LIST(J)%COMPLETE == 1)) .AND. GROUP_LIST(J)%GROUP_SIZE > 0 ) THEN
             IF (T > GROUP_LIST(J)%TDET) THEN
                IF (GROUP_LIST(J)%COMPLETE == 0) THEN
                   GROUP_LIST(J)%TDOOR = MAX(T,GROUP_LIST(J)%TDET)
                END IF
                GROUP_LIST(J)%COMPLETE = 1
             END IF
          END IF
       END DO

       ! ========================================================
       ! Use blocks to speed up double loops in the force loop.
       ! ========================================================
       ALLOCATE(BLOCK_GRID(1:IBAR,1:JBAR,MAX_HUMANS_CELL),STAT=IZERO)
       CALL CHKMEMERR('EVACUATE_HUMANS','BLOCK_GRID',IZERO)
       BLOCK_GRID = 0
       ! BLOCK_GRID_N(:,:) How many agents in this cell.
       ! BLOCK_GRID(:,:,1-N) Agent indeces

       BLOCK_GRID_N = 0
       DO I = 1, N_HUMANS
          HR => HUMAN(I)
          ! Which grid (block grid) cell, new coordinates:
          IIN = FLOOR( CELLSI(FLOOR((HR%X-XS)*RDXINT)) + 1.0_EB )
          JJN = FLOOR( CELLSJ(FLOOR((HR%Y-YS)*RDYINT)) + 1.0_EB )
          BLOCK_GRID_N(IIN,JJN) = BLOCK_GRID_N(IIN,JJN) + 1
          BLOCK_GRID(IIN,JJN,BLOCK_GRID_N(IIN,JJN)) = I
       END DO
       I_DX = 2*(INT((2.0_EB*0.3_EB+5.0_EB)/DX_MIN) + 1) + 1
       J_DY = 2*(INT((2.0_EB*0.3_EB+5.0_EB)/DY_MIN) + 1) + 1
       I_DX = MAX(1,I_DX)
       J_DY = MAX(1,J_DY)
       BL_MAX = I_DX*J_DY*MAX_HUMANS_CELL
       ! Next list will contain the agents that should be looped over for the present agent.
       ALLOCATE(BLOCK_LIST(BL_MAX),STAT=IZERO)
       CALL CHKMEMERR('EVACUATE_HUMANS','BLOCK_LIST',IZERO)

       ! ========================================================
       ! Step (3) of SC-VV starts here: Calculate new forces
       ! ========================================================
       EVAC_FORCE_LOOP: DO I=1,N_HUMANS
          HR => HUMAN(I)
          IF (HR%IEL > 0) THEN
             ! The agent is from some evac line
             GUARD_MEN = EVAC_EVACS(HR%IEL)%GUARD_MEN
             TARGET_X = EVAC_EVACS(HR%IEL)%TARGET_X
             TARGET_Y = EVAC_EVACS(HR%IEL)%TARGET_Y
          ELSE
             ! The agent is from some entr line, not yet as an input
             ! GUARD_MEN = EVAC_ENTRYS(ABS(HR%IEL))%GUARD_MEN
             GUARD_MEN = 0
             TARGET_X = 0.5_EB*(XS+XF)
             TARGET_Y = 0.5_EB*(YS+YF)
          END IF
          IF (MAX(0,HR%GROUP_ID)==0) THEN
             TPRE = HR%TPRE + HR%TDET ! Lonely agent
          ELSE
             TPRE = HR%TDET ! Member of a group
          END IF
          IF (HR%CROWBAR_READ_IN .AND. T <= TPRE+0.01*EVAC_DT_MIN .AND. ((T+DTSP) > TPRE)) THEN
             HR%U = HR%U_CB; HR%V = HR%V_CB
             EVEL = MAX(0.01_EB,SQRT(HR%U_CB**2+HR%V_CB**2))
             HR%UBAR = HR%U_CB/EVEL ; HR%VBAR = HR%V_CB/EVEL
             UBAR = HR%UBAR ; VBAR = HR%VBAR
          END IF
          J  =  MAX(0,HR%GROUP_ID)   ! Group index
          J1 = -MIN(0,HR%GROUP_ID)   ! Lonely agent index
          ! HR%Z is the real z-coordinate of the agent (inclines, stairs,etc),
          ! HR_Z is z-coordinate of the main evac mesh
          HR_Z = 0.5_EB*(ZS+ZF)
          ! Calculate the position and velocities of the shoulder cirles
          Y_TMP(1) = HR%Y - COS(HR%ANGLE)*HR%D_SHOULDER ! RIGHT
          X_TMP(1) = HR%X + SIN(HR%ANGLE)*HR%D_SHOULDER ! RIGHT
          Y_TMP(2) = HR%Y ! TORSO
          X_TMP(2) = HR%X ! TORSO
          Y_TMP(3) = HR%Y + COS(HR%ANGLE)*HR%D_SHOULDER ! LEFT
          X_TMP(3) = HR%X - SIN(HR%ANGLE)*HR%D_SHOULDER ! LEFT
          R_TMP(1) = HR%R_SHOULDER
          R_TMP(2) = HR%R_TORSO
          R_TMP(3) = HR%R_SHOULDER
          U_TMP(1) = HR%U + COS(HR%ANGLE)*HR%OMEGA*HR%D_SHOULDER ! RIGHT
          V_TMP(1) = HR%V + SIN(HR%ANGLE)*HR%OMEGA*HR%D_SHOULDER ! RIGHT
          U_TMP(2) = HR%U ! TORSO
          V_TMP(2) = HR%V ! TORSO
          U_TMP(3) = HR%U - COS(HR%ANGLE)*HR%OMEGA*HR%D_SHOULDER ! LEFT
          V_TMP(3) = HR%V - SIN(HR%ANGLE)*HR%OMEGA*HR%D_SHOULDER ! LEFT
          HR_A = HR%A
          HR_B = HR%B
          HR_KAPPA = HR%KAPPA
          HR_GAMMA = HR%GAMMA
          IF (T <= T_BEGIN) THEN
             HR_KAPPA = 0.0_EB ; HR_GAMMA = 0.0_EB
          END IF
          IF (HR%CROWBAR_READ_IN) THEN
             ! Crowbar project: Use camera data to steer the agents
             HR_A = 0.5_EB*HR_A
             HR_B = 0.5_EB*HR_B
             HR_KAPPA = HR%KAPPA
             HR_GAMMA = HR%GAMMA
          END IF

          CONTACT_F  = 0.0_EB
          CONTACT_FX = 0.0_EB
          CONTACT_FY = 0.0_EB
          SOCIAL_F  = 0.0_EB
          LAMBDAW = LAMBDA_WALL
          A_WALL  = FAC_A_WALL*HR_A
          B_WALL  = FAC_B_WALL*HR_B
          GAME    = NOISEME
          GATH    = NOISETH
          GACM    = NOISECM
          I_COUNT_DENSITY  = 0
          I_COUNT_DENSITYR = 0
          I_COUNT_DENSITYL = 0
          D_HUMANS = HUGE(D_HUMANS)
          D_WALLS  = HUGE(D_WALLS)
          L_FALLEN_DOWN = .FALSE.
          IF (HR%T_FallenDown <= T) THEN
             L_FALLEN_DOWN = .TRUE.
             GATH = 0.0_EB
             A_WALL = 0.0_EB
          END IF
          L_DEAD  = .FALSE.
          IF (HR%INTDOSE >= 1.0_EB) THEN
             L_DEAD = .TRUE.
             ! No random force for a dead person.
             GATH = 0.0_EB
             ! No psychological force terms for a dead person.
             A_WALL = 0.0_EB
             IF (HR%TPRE /= HUGE(HR%TPRE)) THEN
                IF (L_FALLEN_DOWN) N_DEAD = N_DEAD-1 ! Do not duoble count
                N_DEAD = N_DEAD+1
                WRITE (LU_EVACOUT,FMT='(A,I6,A,F8.2,A,I6)') ' EVAC: Agent n:o', HR%ILABEL, &
                     ' dead at ', T, ' s, number of casualties ', N_DEAD
             END IF
             IF (.NOT. L_FALLEN_DOWN) THEN
                HR%T_FallenDown = T
                L_FALLEN_DOWN = .TRUE.
             END IF
             HR%TDET = HUGE(HR%TDET)
             HR%TPRE = HUGE(HR%TPRE)
             HR%COLOR_INDEX = EVAC_AVATAR_NCOLOR
          END IF
          IF (L_FALLEN_DOWN) HR%COLOR_INDEX = EVAC_AVATAR_NCOLOR
          HR_TAU = HR%TAU
          HR_TAU_INER = HR%TAU_INER
          IF (MAX(0,HR%GROUP_ID)==0) THEN
             TPRE = HR%TPRE + HR%TDET ! Lonely agent
          ELSE
             TPRE = HR%TDET ! Member of a group
          END IF
          IF (T < T_BEGIN .OR. T < TPRE) THEN
             HR_TAU = MAX(CF_MIN_TAU, 0.1_EB*HR%TAU)
             HR_TAU_INER = MAX(CF_MIN_TAU_INER, 0.1_EB*HR%TAU_INER)
          END IF
          ! =======================================================
          ! Speed dependent social force
          ! =======================================================
          IF (C_HAWK >= 0.0_EB .AND. HR%I_DoorAlgo==1) THEN
             HR_A = HR%DoseCrit1
          END IF
          HR_A =  HR_A*MAX(0.5_EB,(SQRT(HR%U**2+HR%V**2)/HR%SPEED))
          A_WALL = MIN(A_WALL, FAC_A_WALL*HR_A)

          HR_Speed = HR%SPEED
          MAX_V0_FAC = EVAC_PERSON_CLASSES(HR%IPC)%MAXIMUM_V0_FACTOR
          IF (MAX_V0_FAC>0.0_EB .AND. T > TPRE) THEN
             !Issue1547: Nervousness related parameter is taken from the agent's person type.
             !Issue1547: Later just the MAX_V0_FAC can be used inside this EVAC_FORCE_LOOP.
             !Issue1547: Nervousness related speed change, read the comments in the evac_move_loop above.
             Eta_Nervous = MIN(MAX(1.0_EB -  HR%Speed_ave, 0.0_EB),1.0_EB) ! To be sure that 0 <= Eta <= 1
             Speed_Nervous = (1.0_EB-Eta_Nervous)*HR%Speed + Eta_Nervous*MAX_V0_FAC
             HR_Speed = Speed_Nervous
          END IF

          ! Counterflow: increase motivation to go ahead and decrease social force
          IF (HR%COMMITMENT > 0.01_EB) THEN
             EVEL = MIN(1.0_EB,SQRT(HR%U**2+HR%V**2)/HR_SPEED)
             EVEL = HR%COMMITMENT*EVEL + (1.0_EB-HR%COMMITMENT)*1.0_EB
             HR_TAU      = MAX(CF_MIN_TAU, &
                  HR%COMMITMENT*CF_FAC_TAUS*HR_TAU + (1.0_EB-HR%COMMITMENT)*HR_TAU)
             HR_TAU_INER = MAX(CF_MIN_TAU_INER, &
                  HR%COMMITMENT*CF_FAC_TAUS*HR_TAU_INER + (1.0_EB-HR%COMMITMENT)*HR_TAU_INER)
             A_WALL = MIN(A_WALL, CF_FAC_A_WALL*FAC_A_WALL*HR_A)
             HR_A =  HR_A*MAX(CF_MIN_A,EVEL)
             HR_B =  HR_B*MAX(CF_MIN_B,EVEL)
          END IF
          ! Hawk - dove game: more random noise for hawks
          IF (C_HAWK >= 0.0_EB .AND. HR%I_DoorAlgo==1 .AND. HR%A <= 500.0_EB) THEN
             GATH = GATH*THETA_HAWK
          END IF
          L_FALLDOWN_THIS_TIME_STEP = .FALSE.
          IF (.NOT. L_FALLEN_DOWN .AND. TAU_FALL_DOWN >= 0.0_EB .AND. T-HR%T_CheckFallDown > 0.1_EB) THEN
             ! Check later if this agent falls down, do not do it at every time step
             HR%T_CheckFallDown = T ! Do the check at this time step, save the check time
             CALL RANDOM_NUMBER(RN_REAL)
             RN = REAL(RN_REAL,EB)
             ! P[falls down during Delta_t] = 1 - exp(-lambda*Delta_t), now Delta_t = 0.1 s
             IF (1.0_EB - EXP(-PROB_FALL_DOWN*0.1_EB) > RN) L_FALLDOWN_THIS_TIME_STEP = .TRUE.
          END IF

          ! Psychological force: Cut-off when acceleration below 0.0001 m/s**2
          ! P2P_DIST_MAX = HR%B*LOG(HR%A/0.0001_EB)
          P2P_DIST_MAX = HR_B*LOG(HR_A/0.0001_EB)
          P2P_DIST_MAX = MIN(P2P_DIST_MAX, 5.0_EB)  ! 5.0 m is the maximum range of pp-force
          IF ( HR%SUMFORCES2 > 0.1_EB ) THEN
             ! If large pressure then short range forces only (speed up)
             ! P2P_DIST_MAX = MIN( P2P_DIST_MAX, -HR%B*LOG(HR%SUMFORCES2/(100.0_EB*HR%A)) )
             P2P_DIST_MAX = MIN( P2P_DIST_MAX, -HR_B*LOG(HR%SUMFORCES2/(100.0_EB*HR_A)) )
          END IF
          P2P_DIST_MAX = MAX(P2P_DIST_MAX, 3.0_EB*HR_B)
          ! Next is the max distance for the collision avoidance, counterflow, etc.

          ! In which grid cell is the agent, new coordinates:
          XI = CELLSI(FLOOR((HR%X-XS)*RDXINT))
          YJ = CELLSJ(FLOOR((HR%Y-YS)*RDYINT))
          ZK = CELLSK(FLOOR((HR_Z-ZS)*RDZINT))
          IIN  = FLOOR(XI+1.0_EB)
          JJN  = FLOOR(YJ+1.0_EB)
          KKN = 1
          IIX = FLOOR(XI+0.5_EB)
          JJY = FLOOR(YJ+0.5_EB)
          KKZ = 1
          ICN = CELL_INDEX(IIN,JJN,KKN)
          X1 = HR%X
          Y1 = HR%Y
          HR%W = 0.0_EB

          ! ========================================================
          ! Calculate persons prefered walking direction
          ! ========================================================
          N = NM_STRS_INDEX
          IF(GUARD_MEN == 1) THEN
!!$             CALL GUARD(NM, XI, YJ, WSPA, WSPB, TARGET_X, TARGET_Y)
             UBAR = WSPA
             VBAR = WSPB
             EVEL = SQRT(UBAR**2 + VBAR**2)
             UBAR = UBAR/EVEL ; VBAR = VBAR/EVEL
          ELSE
             CALL FIND_PREFERRED_DIRECTION(I, N, T+DTSP_NEW, T_BEGIN, L_DEAD, NM_STRS_MESH, &
                  IIN, JJN, XI, YJ, UBAR, VBAR, HR_TAU, TPRE, NM, I_STRS_DOOR, HR_SPEED)
             IF (STOP_STATUS>NO_STOP) RETURN
          END IF
          EVEL = UBAR**2 + VBAR**2

          ! (UBAR,VBAR) is now the direction of the flow field of the evacuation mesh
          ! leading towards the chosen door (or the main evacuation field).  It has
          ! no counterflow corrections yet.
          ! ========================================================
          ! The prefered walking direction V0 is (UBAR,VBAR)
          ! ========================================================

          ! =======================================================
          ! Update the block_grid array search ranges
          ! =======================================================
          !
          ! Collision avoidance (incl. counterflow)
          ! Do not do collision avoidance on every time step, do it every 0.1 s on the average by default.
          ! No need for collision avoidance if the target velocity is zero.
          IF (TAU_CHANGE_V0 >= TWO_EPSILON_EB .AND. EVEL >= TWO_EPSILON_EB) THEN
             CALL RANDOM_NUMBER(RN_REAL)
             RNCF = REAL(RN_REAL,EB)
          ELSE
             RNCF = -1.0_EB ; ANGLE_OLD = 0.0_EB ; COMMITMENT = 0.0_EB
          END IF
          IF ( RNCF > EXP(-DTSP/TAU_CHANGE_V0) ) THEN
             ! Collision avoidance needs some range, at least 3.0 m
             P2P_SUUNTA_MAX = MAX(P2P_DIST_MAX, 3.0_EB)
          ELSE
             P2P_SUUNTA_MAX = P2P_DIST_MAX
          END IF

          ! Initialize hawk - dove game utility functions for the current agent (HR%)
          L_DO_GAME = .FALSE.
          U_i_Hawk  = 0.0_EB ; U_i_Dove  = 0.0_EB
          N_i_Hawk  = 0      ; N_i_Dove  = 0
          N_i_Hawk2 = 0      ; N_i_Dove2 = 0
          Hawk_Dove_Game: IF (C_HAWK >= 0.0_EB .AND. T > HR%TPRE+HR%TDET .AND. HR%I_Target>0 .AND. HR%I_DoorAlgo==1) THEN
             ! Hawk - Dove game, only for visible target doors/exits (i_target>0)
             CALL RANDOM_NUMBER(RN_REAL)
             RN = REAL(RN_REAL,EB)
             IF (L_DEAD .OR. L_FALLEN_DOWN .OR. T > T_STOP_HD_GAME) RN = -1.0_EB ! Do not play
             Play_This_Time_Step: IF (RN > EXP(-DTSP/DT_GROUP_DOOR)) THEN
                L_DO_GAME = .TRUE.
                I_TMP = ABS(HR%I_TARGET)
                P2P_SUUNTA_MAX = MAX(P2P_DIST_MAX, P2P_SUUNTA_MAX, 2.0_EB*R_HAWK_DOVE+2.0_EB*HR%RADIUS)
                IF (I_TMP <= N_DOORS) THEN
                   x_o = 0.5_EB*(EVAC_DOORS(I_TMP)%X1 + EVAC_DOORS(I_TMP)%X2)
                   y_o = 0.5_EB*(EVAC_DOORS(I_TMP)%Y1 + EVAC_DOORS(I_TMP)%Y2)
                   T_tmp1 = 50.0_EB*SQRT((X1-x_o)**2 + (Y1-y_o)**2)/EVAC_DOORS(I_TMP)%R_NTARGET + 1.0_EB
                   ii = MIN(50,MAX(1,INT(T_tmp1)-1))
                   N_queue = EVAC_DOORS(I_TMP)%NTARGET(ii)
                   Width = EVAC_DOORS(I_TMP)%Width
                ELSE  ! The target is an exit
                   I_TMP = I_TMP - N_DOORS ! The exit index
                   x_o = 0.5_EB*(EVAC_EXITS(I_TMP)%X1 + EVAC_EXITS(I_TMP)%X2)
                   y_o = 0.5_EB*(EVAC_EXITS(I_TMP)%Y1 + EVAC_EXITS(I_TMP)%Y2)
                   T_tmp1 = 50.0_EB*SQRT((X1-x_o)**2 + (Y1-y_o)**2)/EVAC_EXITS(I_TMP)%R_NTARGET + 1.0_EB
                   ii = MIN(50,MAX(1,INT(T_tmp1)-1))
                   N_queue = EVAC_EXITS(I_TMP)%NTARGET(ii)
                   Width = EVAC_EXITS(I_TMP)%Width
                END IF
                IF (ABS(FAC_DOOR_QUEUE) >= TWO_EPSILON_EB) THEN
                   ! Queueing time is included
                   T_tmp  = SQRT((x_o-X1)**2 + (y_o-Y1)**2)
                   T_tmp1 = REAL(N_queue,EB)/(ABS(FAC_DOOR_QUEUE)*Width)
                   ! T_tmp1 is the estimated queueing time, maximum is 3 p/m2 density and a semi circle
                   IF (FAC_DOOR_QUEUE < -0.001_EB) THEN
                      ! Estimated time is the longer one of the queueing and walking times
                      T_tmp = MAX((T_tmp/HR_Speed), T_tmp1)
                   ELSE
                      ! Estimated time is the queueing time plus the walking time
                      ! Door selection algorithm: alpha*t_walk + (1-alpha)*t_queue
                      T_tmp = (FAC_DOOR_ALPHA*(T_tmp/HR_Speed) +  (1.0_EB-FAC_DOOR_ALPHA)*T_tmp1) / &
                           MAX(1.0_EB-FAC_DOOR_ALPHA,FAC_DOOR_ALPHA)
                   END IF
                ELSE
                   ! No queueing time, just the walking time
                   T_tmp = SQRT((x_o-X1)**2 + (y_o-Y1)**2) / HR_Speed
                END IF
             END IF Play_This_Time_Step
          END IF Hawk_Dove_Game

          ! Speed up the dead agent loop, only contact forces are needed.
          IF (L_DEAD .OR. L_FALLEN_DOWN) P2P_DIST_MAX = 0.0_EB
          IF (L_DEAD .OR. L_FALLEN_DOWN) P2P_SUUNTA_MAX = 0.0_EB

          BLOCK_LIST = 0
          I_DX = INT((2.0_EB*0.3_EB+MAX(P2P_SUUNTA_MAX,P2P_DIST_MAX))/DX(IIN)) + 1
          J_DY = INT((2.0_EB*0.3_EB+MAX(P2P_SUUNTA_MAX,P2P_DIST_MAX))/DY(JJN)) + 1
          IIO = MAX(1,IIN-I_DX)
          IIE = MIN(IBAR,IIN+I_DX)
          JJO = MAX(1,JJN-J_DY)
          JJE = MIN(JBAR,JJN+J_DY)
          IE_MAX = 0
          DO IIX = IIO, IIE
             DO JJY = JJO, JJE
                BG_LOOP: DO J = 1, BLOCK_GRID_N(IIX,JJY)
                   IE_MAX = IE_MAX + 1
                   BLOCK_LIST(MIN(IE_MAX,BL_MAX)) = BLOCK_GRID(IIX,JJY,J)
                END DO BG_LOOP
             END DO
          END DO
          IF (IE_MAX > BL_MAX ) THEN
             WRITE(LU_ERR,'(A,2I6)') 'ERROR: EVACUATE_HUMANS, ie_max, bl_max ', IE_MAX, BL_MAX
             WRITE(LU_EVACOUT,'(A,2I6)') 'ERROR: EVACUATE_HUMANS, ie_max, bl_max ', IE_MAX, BL_MAX
             STOP_STATUS = EVACUATION_STOP ; RETURN
          END IF

          ! =======================================================
          ! Inclines: Velocities are along the incline
          !           coordinates are projected on the (x,y) plane
          ! =======================================================
          COS_X = 1.0_EB
          COS_Y = 1.0_EB
          SPEED_XM = HR_SPEED
          SPEED_XP = HR_SPEED
          SPEED_YM = HR_SPEED
          SPEED_YP = HR_SPEED
          IF (NM_STRS_MESH) THEN
             STRS_INDX = N
             STRP=>EVAC_STRS(N)
          END IF
          ! Check if an agent is on a spectator stand.
          IF (.NOT. NM_STRS_MESH) HR%Z = 0.5_EB*(ZS+ZF)  ! The agent is not on any incline
          SS_LOOP2: DO J = 1, N_SSTANDS
             ESS => EVAC_SSTANDS(J)
             IF (ESS%IMESH == NM .AND. (ESS%X1 <= HR%X .AND. ESS%X2 >= HR%X) .AND. (ESS%Y1 <= HR%Y .AND. ESS%Y2 >= HR%Y) ) THEN
                IF (ESS%USE_V0 .AND. .NOT.(HR%I_FFIELD == ESS%I_VENT_FFIELD .OR. ESS%I_VENT_FFIELD == 0)) CYCLE SS_LOOP2
                ! Next are here just for a test case.  The evacuation flow fields
                ! in front of the doors/exits are not optimal.  Close to a door it
                ! would be better to use an ad-hoc movement direction, something
                ! like is done in the STRS staircase model.
                IF (.NOT.(L_DEAD .OR. L_FALLEN_DOWN) .AND. ESS%USE_V0) THEN
                   EVEL = SQRT(ESS%UBAR0**2 + ESS%VBAR0**2)
                   IF (HR%I_FFIELD == ESS%I_VENT_FFIELD .OR. ESS%I_VENT_FFIELD == 0 .AND. EVEL > 1.0E-12_EB) THEN
                      UBAR = ESS%UBAR0/EVEL
                      VBAR = ESS%VBAR0/EVEL
                   END IF
                END IF
                COS_X = ESS%COS_X
                COS_Y = ESS%COS_Y
                FAC_V0_UP   = ESS%FAC_V0_UP
                FAC_V0_DOWN = ESS%FAC_V0_DOWN
                FAC_V0_HORI = ESS%FAC_V0_HORI
                IF (EVAC_PERSON_CLASSES(HR%IPC)%FAC_V0_HORI >= TWO_EPSILON_EB) THEN
                   FAC_V0_HORI = EVAC_PERSON_CLASSES(HR%IPC)%FAC_V0_UP
                END IF
                IF (EVAC_PERSON_CLASSES(HR%IPC)%FAC_V0_UP >= TWO_EPSILON_EB) THEN
                   IF ((ESS%H - ESS%H0) <= -TWO_EPSILON_EB) THEN
                      FAC_V0_DOWN = EVAC_PERSON_CLASSES(HR%IPC)%FAC_V0_UP
                   ELSE
                      FAC_V0_UP = EVAC_PERSON_CLASSES(HR%IPC)%FAC_V0_UP
                   END IF
                END IF
                IF (EVAC_PERSON_CLASSES(HR%IPC)%FAC_V0_DOWN >= TWO_EPSILON_EB) THEN
                   IF ((ESS%H - ESS%H0) <= -TWO_EPSILON_EB) THEN
                      FAC_V0_UP = EVAC_PERSON_CLASSES(HR%IPC)%FAC_V0_DOWN
                   ELSE
                      FAC_V0_DOWN = EVAC_PERSON_CLASSES(HR%IPC)%FAC_V0_DOWN
                   END IF
                END IF
                SELECT CASE (ESS%IOR)
                CASE(-1)
                   SPEED_XM = COS_X*HR_SPEED*FAC_V0_DOWN
                   SPEED_XP = COS_X*HR_SPEED*FAC_V0_UP
                   SPEED_YM = HR_SPEED*FAC_V0_HORI
                   SPEED_YP = HR_SPEED*FAC_V0_HORI
                   HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%X1-HR%X)/ABS(ESS%X1-ESS%X2)
                CASE(+1)
                   SPEED_XM = COS_X*HR_SPEED*FAC_V0_UP
                   SPEED_XP = COS_X*HR_SPEED*FAC_V0_DOWN
                   SPEED_YM = HR_SPEED*FAC_V0_HORI
                   SPEED_YP = HR_SPEED*FAC_V0_HORI
                   HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%X2-HR%X)/ABS(ESS%X1-ESS%X2)
                CASE(-2)
                   SPEED_XM = HR_SPEED*FAC_V0_HORI
                   SPEED_XP = HR_SPEED*FAC_V0_HORI
                   SPEED_YM = COS_Y*HR_SPEED*FAC_V0_DOWN
                   SPEED_YP = COS_Y*HR_SPEED*FAC_V0_UP
                   HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%Y1-HR%Y)/ABS(ESS%Y1-ESS%Y2)
                CASE(+2)
                   SPEED_XM = HR_SPEED*FAC_V0_HORI
                   SPEED_XP = HR_SPEED*FAC_V0_HORI
                   SPEED_YM = COS_Y*HR_SPEED*FAC_V0_UP
                   SPEED_YP = COS_Y*HR_SPEED*FAC_V0_DOWN
                   HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%Y2-HR%Y)/ABS(ESS%Y1-ESS%Y2)
                END SELECT
             END IF
          END DO SS_LOOP2
          ! =======================================================
          ! Inclines ends
          ! =======================================================

          ! =======================================================
          ! Set height and speed for an agent in stairs
          ! =======================================================
          IF (NM_STRS_MESH) THEN
            CALL GETSTAIRSPEEDANDZ(SPEED_XM, SPEED_XP, SPEED_YM, SPEED_YP, HR_SPEED, STRP, HR)
            IF (STOP_STATUS>NO_STOP) RETURN
          END IF

          ! ========================================================
          ! Agent-agent interaction forces
          ! ========================================================
          ! Look for other agents, use blocks to speed up the loop.
          P2P_U      = 0.0_EB
          P2P_V      = 0.0_EB
          P2P_TORQUE = 0.0_EB
          ! ========================================================
          ! Collision avoidance (incl. counterflow)
          ! ========================================================
          CHANGE_V0_RNCF0: IF ( RNCF > EXP(-DTSP/TAU_CHANGE_V0) ) THEN

             V_HR = MAX(0.1_EB,MIN(1.0_EB,SQRT(HR%U**2 + HR%V**2)/HR_SPEED))
             THETA_START = -ABS(THETA_SECTOR)
             THETA_STEP = 2.0_EB*ABS(THETA_START)/REAL(N_SECTORS-1,EB)
             EVEL = MAX(0.0_EB,MIN(1.0_EB,(HR%U**2+HR%V**2)/HR_SPEED))
             THETA_START = THETA_START - 0.5_EB*(1.0_EB-EVEL)*MAX(0.0_EB,(90.0_EB+THETA_START-0.5_EB*THETA_STEP))
             THETA_STEP = 2.0_EB*ABS(THETA_START)/REAL(N_SECTORS-1,EB)
             IF (HR%UBAR**2 + HR%VBAR**2 < 0.1_EB) THEN
                HR%UBAR = UBAR; HR%VBAR = VBAR
             END IF
             SUM_SUUNTA = 0.0_EB
             DO III = 1, N_SECTORS
                THETAS(III)    = -ABS(THETA_START) + (III-1)*THETA_STEP ! Degrees
                COS_THETA(III) = COS(PI*THETAS(III)/180._EB) ! Radians in fortran functions
                SIN_THETA(III) = SIN(PI*THETAS(III)/180._EB)
                U_THETA(III) = COS_THETA(III)*UBAR - SIN_THETA(III)*VBAR
                V_THETA(III) = SIN_THETA(III)*UBAR + COS_THETA(III)*VBAR
                CALL RANDOM_NUMBER(RN_REAL)
                RN = REAL(RN_REAL,EB)
                SUM_SUUNTA(III) = 0.000001_EB*RN ! left/right symmetry trick. Now the order of the sectors do not matter.
             END DO
             SUM_SUUNTA(N_SECTORS+1) = 0.1_EB ! prefer straight ahead a little bit
             COS_THETA(N_SECTORS+1) = 1.0_EB ! V0 direction
             SIN_THETA(N_SECTORS+1) = 0.0_EB
             THETAS(N_SECTORS+1)    = 0.0_EB
             U_THETA(N_SECTORS+1)   = UBAR
             V_THETA(N_SECTORS+1)   = VBAR
             ANGLE_OLD       = 0.0_EB
             N_SUUNTA        = 0 ! How many agents per sector
             N_SUUNTACF      = 0 ! How many agents per sector (counterflow)
             N_SUUNTA_BACK   = 0 ! How many agents "behind"
             N_SUUNTA_BACKCF = 0 ! How many agents "behind" (counterflow)
          END IF CHANGE_V0_RNCF0

          ! The force loop over the other agents starts here
          P2PLOOP: DO IE = 1, IE_MAX
             ! I: Index of the current agent (outer loop), IE: Index of the other agent (inner loop)
             IF (BLOCK_LIST(IE) == I) CYCLE P2PLOOP  ! No self interaction
             HRE => HUMAN(BLOCK_LIST(IE))
             L_HRE_FALLEN_DOWN = .FALSE.
             IF (HRE%T_FallenDown < T) L_HRE_FALLEN_DOWN = .TRUE.
             ! In stairs, only consider humans at the same, next and previous sub nodes (landings, stairs)
             IF (NM_STRS_MESH .AND. ABS(HRE%STR_SUB_INDX - HR%STR_SUB_INDX)>1) CYCLE P2PLOOP

             P2P_DIST = ((HRE%X-X1)**2 + (HRE%Y-Y1)**2)
             Hawk_Dove_Game1: IF (L_DO_GAME .AND. P2P_DIST-(HR%RADIUS+HRE%RADIUS)**2 < R_HAWK_DOVE**2) THEN
                ! Hawk - Dove game
                IF (HRE%A > 500.0_EB) THEN
                   N_i_Dove = N_i_Dove + 1
                ELSE
                   N_i_Hawk = N_i_Hawk + 1
                END IF
             END IF Hawk_Dove_Game1
             Hawk_Dove_Game1B: IF (L_DO_GAME .AND. P2P_DIST-(HR%RADIUS+HRE%RADIUS)**2 < 4.0_EB*R_HAWK_DOVE**2) THEN
                ! Hawk - Dove game
                IF (HRE%A > 500.0_EB) THEN
                   N_i_Dove2 = N_i_Dove2 + 1
                ELSE
                   N_i_Hawk2 = N_i_Hawk2 + 1
                END IF
             END IF Hawk_Dove_Game1B
             IF (P2P_DIST < (MIN(MAX(P2P_DIST_MAX,P2P_SUUNTA_MAX),1.5_EB)+HR%RADIUS+HRE%RADIUS)**2) THEN
                I_COUNT_DENSITY = I_COUNT_DENSITY + 1
                EVEL = SQRT(HR%UBAR**2 + HR%VBAR**2)
                IF (EVEL >= TWO_EPSILON_EB) THEN
                   EVEL = MAX(0.001_EB,SQRT(HR%UBAR**2 + HR%VBAR**2))
                   IF (HR%VBAR >= 0.0_EB) THEN
                      ANGLE_HR = ACOS(HR%UBAR/EVEL)
                   ELSE
                      ANGLE_HR = 2.0_EB*PI - ACOS(HR%UBAR/EVEL)
                   END IF
                   IF (ANGLE_HR == 2.0_EB*PI) ANGLE_HR = 0.0_EB  ! Agent HR angle is [0,2PI)
                ELSE
                   ANGLE_HR = HR%ANGLE
                END IF

                D_NEW = SQRT((HRE%X-HR%X)**2 + (HRE%Y-HR%Y)**2)
                IF ((HRE%Y-HR%Y) >= TWO_EPSILON_EB) THEN
                   ANGLE_HRE = ACOS((HRE%X-HR%X)/D_NEW)
                ELSE
                   ANGLE_HRE = 2.0_EB*PI - ACOS((HRE%X-HR%X)/D_NEW)
                END IF
                IF (ANGLE_HRE == 2.0_EB*PI) ANGLE_HRE = 0.0_EB  ! Agent HRE angle is [0,2PI)
                ANGLE_HRE = ANGLE_HR - ANGLE_HRE
                IF (ANGLE_HRE >  PI) ANGLE_HRE = ANGLE_HRE - 2.0_EB*PI
                IF (ANGLE_HRE < -PI) ANGLE_HRE = ANGLE_HRE + 2.0_EB*PI
                ! If HRE is on the right hand side of HR then the angle is negative,
                ! i.e., positive direction is anti-clockwise.
                ANGLE_HRE = -180.0_EB*ANGLE_HRE/PI  ! Degrees
                IF (ANGLE_HRE >= -90.0_EB .AND. ANGLE_HRE <=  0.0_EB) I_COUNT_DENSITYR = I_COUNT_DENSITYR + 1
                IF (ANGLE_HRE >=   0.0_EB .AND. ANGLE_HRE <= 90.0_EB) I_COUNT_DENSITYL = I_COUNT_DENSITYL + 1
             END IF
             IF (P2P_DIST > (MAX(P2P_DIST_MAX,P2P_SUUNTA_MAX)+HR%RADIUS+HRE%RADIUS)**2) CYCLE P2PLOOP
             P2P_DIST = SQRT(P2P_DIST)
             !
             ! Check, that the persons are seeing each other, i.e., there are no walls between.
             PP_SEE_EACH = SEE_EACH_OTHER(NM, X1, Y1, HRE%X, HRE%Y)
             IF (.NOT. PP_SEE_EACH) CYCLE P2PLOOP

             ! Collision avoidance, counterflow, etc.
             CHANGE_V0_RNCF1: IF ( RNCF > EXP(-DTSP/TAU_CHANGE_V0) ) THEN
                EVEL = MAX(0.001_EB,SQRT(UBAR**2 + VBAR**2))
                IF (VBAR >= 0.0_EB) THEN
                   ANGLE_HR = ACOS(UBAR/EVEL)
                ELSE
                   ANGLE_HR = 2.0_EB*PI - ACOS(UBAR/EVEL)
                END IF
                IF (ANGLE_HR == 2.0_EB*PI) ANGLE_HR = 0.0_EB  ! Agent HR angle is [0,2PI)

                D_SHIFT = 2.0_EB*(0.5_EB-MIN(1.0_EB,SQRT(HR%U**2+HR%V**2)/HR_SPEED))*HR%RADIUS
                HR_X = HR%X - UBAR*D_SHIFT
                HR_Y = HR%Y - VBAR*D_SHIFT
                D_NEW = SQRT((HRE%X-HR_X)**2 + (HRE%Y-HR_Y)**2)
                IF (D_SHIFT < 0.0_EB .OR. D_NEW > HR%R_TORSO) THEN
                   IF ((HRE%Y-HR_Y) >= 0.0_EB) THEN
                      ANGLE_HRE = ACOS((HRE%X-HR_X)/D_NEW)
                   ELSE
                      ANGLE_HRE = 2.0_EB*PI - ACOS((HRE%X-HR_X)/D_NEW)
                   END IF
                   IF (ANGLE_HRE == 2.0_EB*PI) ANGLE_HRE = 0.0_EB  ! Agent HRE angle is [0,2PI)
                   ANGLE_HRE = ANGLE_HR - ANGLE_HRE

                   IF (ANGLE_HRE >  PI) ANGLE_HRE = ANGLE_HRE - 2.0_EB*PI
                   IF (ANGLE_HRE < -PI) ANGLE_HRE = ANGLE_HRE + 2.0_EB*PI
                   ! Agent hre is at this angle when measured from HR
                   ! If HRE is on the right hand side of HR then the angle is negative,
                   ! i.e., positive direction is anti-clockwise (as usual).
                   ANGLE_HRE = -180.0_EB*ANGLE_HRE/PI  ! Degrees
                   ANGLE_HR  = -180.0_EB*ANGLE_HR /PI  ! Degrees

                   V_HRE = (HRE%X-HR%X)*UBAR + (HRE%Y-HR%Y)*VBAR
                   V_HRE = V_HRE/SQRT((HRE%X-HR%X)**2+(HRE%Y-HR%Y)**2) ! COS R_HRE VS V0
                   V_HRE = MAX(0.5_EB,V_HRE)  ! Look further ahead than sideways
                   EVEL = MAX(0.2_EB,MIN(1.0_EB,SQRT(HR%U**2+HR%V**2)/HR_SPEED))
                   TIM_DIST = MAX(0.2_EB,MIN(EVEL,V_HRE))*P2P_SUUNTA_MAX + HR%RADIUS + HRE%RADIUS
                   EVEL = MAX(0.5_EB,MIN(1.0_EB,SQRT(HR%U**2+HR%V**2)/HR_SPEED))
                   V_HRE = MAX(0.5_EB,MIN(EVEL,V_HRE))*P2P_SUUNTA_MAX ! Look further ahead than sideways
                   IF (P2P_DIST < (V_HRE+HR%RADIUS+HRE%RADIUS) ) THEN
                      VR_2R = HRE%UBAR*UBAR + HRE%VBAR*VBAR ! Counterflow or not?
                      ! Which sector if any
                      III = INT((ANGLE_HRE-(THETA_START-1.5_EB*THETA_STEP))/THETA_STEP)
                      IF (III > 0 .AND. III < N_SECTORS+1) THEN
                         ! (UBAR,VBAR) are unit vectors
                         IF (P2P_DIST < TIM_DIST) N_SUUNTA(III) = N_SUUNTA(III) + 1
                         IF (P2P_DIST < TIM_DIST .AND. VR_2R <= -0.2_EB) N_SUUNTACF(III) = N_SUUNTACF(III) + 1
                         IF (VR_2R >= TWO_EPSILON_EB) THEN ! Same direction
                            V_HRE = HRE%U*UBAR + HRE%V*VBAR  ! HRE speed along the v0 direction
                            V_HR  = MAX(0.0_EB,HR%U*UBAR + HR%V*VBAR)  ! HR speed along the v0 direction
                            V_HRE = CONST_DF + FAC_DF*(MIN(HR%V0_FAC*HR_SPEED,V_HRE) - &
                                 MIN(HR%V0_FAC*HR_SPEED,V_HR))
                         ELSE ! Counterflow
                            V_HRE = HRE%U*UBAR + HRE%V*VBAR  ! HRE SPEED ALONG THE V0 DIRECTION
                            V_HRE = -1.0_EB*(CONST_CF + FAC_CF*MAX(0.0_EB,-V_HRE))
                         END IF
                         SUM_SUUNTA(III) = SUM_SUUNTA(III) + V_HRE/MAX(0.2_EB,(P2P_DIST-HR%RADIUS-HRE%RADIUS))
                         IF (ANGLE_HRE > -0.5_EB*THETA_STEP .AND. ANGLE_HRE < 0.5_EB*THETA_STEP) THEN
                            ! The "additional" sector pointing on the v0 direction
                            N_SUUNTA(N_SECTORS+1) = N_SUUNTA(N_SECTORS+1) + 1
                            IF (VR_2R <= -0.2_EB) N_SUUNTACF(N_SECTORS+1) = N_SUUNTACF(N_SECTORS+1) + 1
                            SUM_SUUNTA(N_SECTORS+1) = SUM_SUUNTA(N_SECTORS+1) + &
                                 V_HRE/MAX(0.2_EB,(P2P_DIST-HR%RADIUS-HRE%RADIUS))
                         END IF
                      END IF
                   END IF
                END IF
             END IF CHANGE_V0_RNCF1

             IF ( P2P_DIST > P2P_DIST_MAX+HR%RADIUS+HRE%RADIUS ) CYCLE P2PLOOP
             !
             ! ========================================================
             ! Calculate the combination of spring constant for the two agents
             ! ========================================================
             C_YEFF = (2.0_EB*HR%C_YOUNG*2.0_EB*HRE%C_YOUNG)/(2.0_EB*HR%C_YOUNG+2.0_EB*HRE%C_YOUNG)
             !
             ! ========================================================
             ! Angle dependent social force:
             ! ========================================================
             IF ( (HR%U**2 +HR%V**2) >= TWO_EPSILON_EB ) THEN
                COSPHIFAC = ( (HRE%X-X1)*HR%U + (HRE%Y-Y1)*HR%V ) &
                     / ( SQRT((HRE%X-X1)**2 + (HRE%Y-Y1)**2)*SQRT(HR%U**2 +HR%V**2) )
                COSPHIFAC = HR%LAMBDA + 0.5_EB*(1.0_EB-HR%LAMBDA)*(1.0_EB+COSPHIFAC)
             ELSE
                COSPHIFAC = 1.0_EB
             END IF

             ! Calculate the position and velocities of the shoulder cirles for HRE
             R_TMP(4) = HRE%R_SHOULDER ! Right circle
             R_TMP(5) = HRE%R_TORSO    ! Center circle
             R_TMP(6) = HRE%R_SHOULDER ! Left circle
             Y_TMP(4) = HRE%Y - COS(HRE%ANGLE)*HRE%D_SHOULDER ! Right circle
             X_TMP(4) = HRE%X + SIN(HRE%ANGLE)*HRE%D_SHOULDER ! Right circle
             Y_TMP(5) = HRE%Y ! Center circle
             X_TMP(5) = HRE%X ! Center circle
             Y_TMP(6) = HRE%Y + COS(HRE%ANGLE)*HRE%D_SHOULDER ! Left circle
             X_TMP(6) = HRE%X - SIN(HRE%ANGLE)*HRE%D_SHOULDER ! Left circle

             ! ========================================================
             ! Add psychological (social) force term
             ! ========================================================
             IF (.NOT. (L_DEAD .OR. L_FALLEN_DOWN)) THEN
                FC_X = 0.0_EB
                FC_Y = 0.0_EB
                IF (L_HRE_FALLEN_DOWN) COSPHIFAC = A_FAC_FALLEN * COSPHIFAC
                ! Use the closest circles to calculate the psychological force
                DO III = 1, 3
                   DO JJJ = 4, 6
                      TIM_DIST = MAX(0.001_EB,SQRT((X_TMP(III)-X_TMP(JJJ))**2 + (Y_TMP(III)-Y_TMP(JJJ))**2))
                      ! Next is |vector1|*|vector2|
                      EVEL = SQRT((X_TMP(JJJ)-X_TMP(III))**2+(Y_TMP(JJJ)-Y_TMP(III))**2)* SQRT(U_TMP(III)**2+V_TMP(III)**2)
                      IF (EVEL >= TWO_EPSILON_EB) EVEL = ((X_TMP(JJJ)-X_TMP(III))*U_TMP(III) + &
                           (Y_TMP(JJJ)-Y_TMP(III))*V_TMP(III)) / EVEL   ! COS THETA (SCAL_PROD/(LENGHT1*LENGTH2)
                      IF (EVEL >= TWO_EPSILON_EB) THEN
                         D_HUMANS = MIN( (TIM_DIST-(R_TMP(III)+R_TMP(JJJ))) /EVEL, D_HUMANS)
                      ELSE
                         D_HUMANS = MIN( (TIM_DIST-(R_TMP(III)+R_TMP(JJJ))) /0.01_EB , D_HUMANS)
                      END IF
                      IF (TAU_FALL_DOWN >= 0.0_EB .AND. TIM_DIST <= R_TMP(III)+R_TMP(JJJ)-D_OVERLAP_FALL) THEN
                         ! Circles are touching each others more than D_OVERLAP_FALL => Fall down
                         IF (L_HRE_FALLEN_DOWN) COSPHIFAC = 0.0_EB ! On top of a fallen agent, no social force anymore
                         IF (T > HR%TDET .AND. L_FALLDOWN_THIS_TIME_STEP .AND. L_HRE_FALLEN_DOWN .AND. .NOT.L_FALLEN_DOWN) THEN
                            L_FALLEN_DOWN = .TRUE.
                            HR%T_FallenDown = T
                            HR%Angle_FallenDown = HR%angle
                            HR%SizeFac_FallenDown = 0.0_EB
                            N_DEAD = N_DEAD+1
                         END IF
                      END IF

                      FC_X1 = (X_TMP(III)-X_TMP(JJJ))*HR_A*COSPHIFAC* &
                           EXP( -(TIM_DIST-( R_TMP(III)+R_TMP(JJJ) ))/HR_B )/TIM_DIST
                      FC_Y1 = (Y_TMP(III)-Y_TMP(JJJ))*HR_A*COSPHIFAC* &
                           EXP( -(TIM_DIST-( R_TMP(III)+R_TMP(JJJ) ))/HR_B )/TIM_DIST
                      IF ( (FC_X1**2+FC_Y1**2) > (FC_X**2+FC_Y**2) ) THEN
                         FC_X = FC_X1
                         FC_Y = FC_Y1
                      END IF
                   END DO
                END DO
                IF (L_HRE_FALLEN_DOWN) THEN ! Maximum force is skin-skin value, i.e., the pre-exp factor
                   EVEL = SQRT(FC_X**2+FC_Y**2)
                   IF (EVEL > HR_A*COSPHIFAC ) THEN
                      FC_X = HR_A*COSPHIFAC * FC_X/EVEL
                      FC_Y = HR_A*COSPHIFAC * FC_Y/EVEL
                   END IF
                END IF
                P2P_U = P2P_U + FC_X
                P2P_V = P2P_V + FC_Y
                SOCIAL_F = SOCIAL_F + SQRT(FC_X**2 + FC_Y**2)
                TC_Z = 0.0_EB
                ! Calculate the torque due to the social force. use the closest circles.
                DO JJJ = 4, 6
                   ! First the right shoulder
                   TIM_DIST = SQRT( (X_TMP(JJJ)-X_TMP(1))**2 + (Y_TMP(JJJ)-Y_TMP(1))**2 )
                   FC_X = (X_TMP(1)-X_TMP(JJJ)) * HR_A*COSPHIFAC*EXP( -(TIM_DIST - (R_TMP(1)+R_TMP(JJJ)))/HR_B )/TIM_DIST
                   FC_Y = (Y_TMP(1)-Y_TMP(JJJ)) * HR_A*COSPHIFAC*EXP( -(TIM_DIST - (R_TMP(1)+R_TMP(JJJ)))/HR_B )/TIM_DIST
                   IF ( ABS(FC_Y*(X_TMP(1)-HR%X) - FC_X*(Y_TMP(1)-HR%Y)) > ABS(TC_Z) ) THEN
                      TC_Z = FC_Y*(X_TMP(1)-HR%X) - FC_X*(Y_TMP(1)-HR%Y)
                   END IF
                END DO
                P2P_TORQUE = P2P_TORQUE + TC_Z
                TC_Z = 0.0_EB
                DO JJJ = 4, 6
                   ! Then the left shoulder
                   TIM_DIST = SQRT( (X_TMP(JJJ)-X_TMP(3))**2 + (Y_TMP(JJJ)-Y_TMP(3))**2 )
                   FC_X = (X_TMP(3)-X_TMP(JJJ)) * HR_A*COSPHIFAC*EXP( -(TIM_DIST - (R_TMP(3)+R_TMP(JJJ)))/HR_B )/TIM_DIST
                   FC_Y = (Y_TMP(3)-Y_TMP(JJJ)) * HR_A*COSPHIFAC*EXP( -(TIM_DIST - (R_TMP(3)+R_TMP(JJJ)))/HR_B )/TIM_DIST
                   IF ( ABS(FC_Y*(X_TMP(3)-HR%X) - FC_X*(Y_TMP(3)-HR%Y)) > ABS(TC_Z) ) THEN
                      TC_Z = FC_Y*(X_TMP(3)-HR%X) - FC_X*(Y_TMP(3)-HR%Y)
                   END IF
                END DO
                P2P_TORQUE = P2P_TORQUE + TC_Z
             END IF  ! Dead or alive?

             ! ========================================================
             ! Add contact force terms
             ! ========================================================
             ! No contact forces for fallen down agent due to other agents (wall contact forces are calculated)
             IF ( (P2P_DIST <= (HR%RADIUS+HRE%RADIUS)) .AND. .NOT.L_FALLEN_DOWN) THEN
                ! Calculate the velocities of the shoulder cirles
                V_TMP(4) = HRE%V + SIN(HRE%ANGLE)*HRE%OMEGA*HRE%D_SHOULDER
                V_TMP(5) = HRE%V
                V_TMP(6) = HRE%V - SIN(HRE%ANGLE)*HRE%OMEGA*HRE%D_SHOULDER
                U_TMP(4) = HRE%U + COS(HRE%ANGLE)*HRE%OMEGA*HRE%D_SHOULDER
                U_TMP(5) = HRE%U
                U_TMP(6) = HRE%U - COS(HRE%ANGLE)*HRE%OMEGA*HRE%D_SHOULDER

                DO III = 1, 3
                   DO JJJ = 4, 6
                      TIM_DIST = SQRT((X_TMP(III)-X_TMP(JJJ))**2 + (Y_TMP(III)-Y_TMP(JJJ))**2)
                      ! Next is |vector1|*|vector2|
                      EVEL = SQRT((X_TMP(JJJ)-X_TMP(III))**2+(Y_TMP(JJJ)-Y_TMP(III))**2)* SQRT(U_TMP(III)**2+V_TMP(III)**2)
                      IF (EVEL >= TWO_EPSILON_EB) EVEL = ((X_TMP(JJJ)-X_TMP(III))*U_TMP(III) + &
                           (Y_TMP(JJJ)-Y_TMP(III))*V_TMP(III)) / EVEL   ! COS THETA (SCAL_PROD/(LENGHT1*LENGTH2)
                      IF (EVEL >= TWO_EPSILON_EB) THEN
                         D_HUMANS = MIN( (TIM_DIST-(R_TMP(III)+R_TMP(JJJ))) /EVEL, D_HUMANS)
                      ELSE
                         D_HUMANS = MIN( (TIM_DIST-(R_TMP(III)+R_TMP(JJJ))) /0.01_EB , D_HUMANS)
                      END IF
!!$                      IF (TAU_FALL_DOWN >= 0.0_EB .AND. TIM_DIST <= R_TMP(III)+R_TMP(JJJ)-D_OVERLAP_FALL) THEN
!!$                         ! Circles are touching each others more than D_OVERLAP_FALL => Fall down
!!$                         !P_FallDownNow = HR%T_CheckFallDown, pers-rivin inputti PROB_FALL_DOWN
!!$                         IF (T > HR%TDET .AND. L_FALLDOWN_THIS_TIME_STEP .AND. L_HRE_FALLEN_DOWN) THEN
!!$                            L_FALLEN_DOWN = .TRUE.
!!$                            HR%T_FallenDown = T
!!$                            HR%Angle_FallenDown = HR%angle
!!$                            HR%SizeFac_FallenDown = 0.0_EB
!!$                         END IF
!!$                      END IF
                      ! No contact forces for standing-fallen down
                      IF (TIM_DIST <= R_TMP(III)+R_TMP(JJJ) .AND. .NOT.L_HRE_FALLEN_DOWN) THEN
                         ! Circles are touching each others
                         FC_X =(X_TMP(III)-X_TMP(JJJ))*C_YEFF*((R_TMP(III)+R_TMP(JJJ))-TIM_DIST)/TIM_DIST
                         FC_Y =(Y_TMP(III)-Y_TMP(JJJ))*C_YEFF*((R_TMP(III)+R_TMP(JJJ))-TIM_DIST)/TIM_DIST
                         FC_X = FC_X - FC_DAMPING*(U_TMP(III)-U_TMP(JJJ))*(X_TMP(III)-X_TMP(JJJ))/TIM_DIST
                         FC_Y = FC_Y - FC_DAMPING*(V_TMP(III)-V_TMP(JJJ))*(Y_TMP(III)-Y_TMP(JJJ))/TIM_DIST
                         CONTACT_F = CONTACT_F + SQRT(FC_X**2 + FC_Y**2)
                         CONTACT_FX = CONTACT_FX + FC_X
                         CONTACT_FY = CONTACT_FY + FC_Y
                         P2P_U = P2P_U + FC_X
                         P2P_V = P2P_V + FC_Y
                         P2P_TORQUE = P2P_TORQUE + FC_Y*(X_TMP(III)-HR%X) - FC_X*(Y_TMP(III)-HR%Y)
                         SCAL_PROD_OVER_RSQR = ((Y_TMP(III)-Y_TMP(JJJ))*(U_TMP(III)-U_TMP(JJJ)) - &
                              (X_TMP(III)-X_TMP(JJJ))*(V_TMP(III)-V_TMP(JJJ))) / (TIM_DIST**2)
                         IF (I_FRIC_SW >= 1 ) THEN  ! This is the default
                            FC_X = - HR_KAPPA*((R_TMP(III)+R_TMP(JJJ))-TIM_DIST)* &
                                 ( (Y_TMP(III)-Y_TMP(JJJ)) * SCAL_PROD_OVER_RSQR )
                            FC_Y = - HR_KAPPA*((R_TMP(III)+R_TMP(JJJ))-TIM_DIST)* &
                                 (-(X_TMP(III)-X_TMP(JJJ)) * SCAL_PROD_OVER_RSQR )
                            P2P_U = P2P_U + FC_X
                            P2P_V = P2P_V + FC_Y
                            P2P_TORQUE = P2P_TORQUE + FC_Y*( (X_TMP(III) + &
                                 (R_TMP(III)/R_TMP(JJJ))*(X_TMP(JJJ)-X_TMP(III)) ) - HR%X )
                            P2P_TORQUE = P2P_TORQUE - FC_X*( (Y_TMP(III) + &
                                 (R_TMP(III)/R_TMP(JJJ))*(Y_TMP(JJJ)-Y_TMP(III)) ) - HR%Y )
                         ELSE
                            FC_X = -HR_GAMMA*( (Y_TMP(III)-Y_TMP(JJJ))*SCAL_PROD_OVER_RSQR)
                            FC_Y = -HR_GAMMA*(-(X_TMP(III)-X_TMP(JJJ))*SCAL_PROD_OVER_RSQR)
                            P2P_U = P2P_U + FC_X
                            P2P_V = P2P_V + FC_Y
                            P2P_TORQUE = P2P_TORQUE + FC_Y*( (X_TMP(III) + &
                                 (R_TMP(III)/R_TMP(JJJ))*(X_TMP(JJJ)-X_TMP(III)) ) - HR%X )
                            P2P_TORQUE = P2P_TORQUE - FC_X*( (Y_TMP(III) + &
                                 (R_TMP(III)/R_TMP(JJJ))*(Y_TMP(JJJ)-Y_TMP(III)) ) - HR%Y )
                         END IF
                         CONTACT_FX = CONTACT_FX + FC_X
                         CONTACT_FY = CONTACT_FY + FC_Y
                      END IF
                   END DO
                END DO
             END IF  ! Contact forces?
          END DO P2PLOOP
          IF (MAX(P2P_DIST_MAX,P2P_SUUNTA_MAX)>0.1_EB) THEN
             HR%DENSITY  = I_COUNT_DENSITY  / (PI*(MIN(MAX(P2P_DIST_MAX,P2P_SUUNTA_MAX),1.5_EB)+HR%RADIUS+HR%RADIUS)**2)
             HR%DENSITYR = I_COUNT_DENSITYR / (0.25_EB*PI*(MIN(MAX(P2P_DIST_MAX,P2P_SUUNTA_MAX),1.5_EB)+HR%RADIUS+HR%RADIUS)**2)
             HR%DENSITYL = I_COUNT_DENSITYL / (0.25_EB*PI*(MIN(MAX(P2P_DIST_MAX,P2P_SUUNTA_MAX),1.5_EB)+HR%RADIUS+HR%RADIUS)**2)
          ELSE
             HR%DENSITY  = 0.0_EB
             HR%DENSITYR = 0.0_EB
             HR%DENSITYL = 0.0_EB
          END IF
          Hawk_Dove_Game2: IF (L_DO_GAME) THEN
             ! Hawk - Dove game
             IF (N_i_Dove + N_i_Hawk == 0) THEN
                N_i_Dove = N_i_Dove2 ; N_i_Hawk = N_i_Hawk2
             END IF
             IF (N_i_Dove + N_i_Hawk > 0) THEN
                p_i = REAL(N_i_Hawk,EB)/REAL(N_i_Dove + N_i_Hawk,EB)
                ! U_i_Hawk = (p_i*C_Hawk/T_tmp) - (1.0_EB-p_i)
                ! Time dependent (linearly growing) threath parameter alpha, in the final version
                ! the threath is estimated from the smoke consentration etc, but for the development
                ! of the game, a simple linear function is introduced.
                ! New Hawk-Dove game parameters:
                IF (T_Tmp > (T_ASET_HAWK+T_ASET_TFAC_HAWK*MAX(0.0_EB,T-A_HAWK_T_START))-T_0_HAWK) THEN
                   !U_i_Hawk = p_i * T_0_HAWK/(T_tmp-T_ASET_HAWK+T_0_HAWK) - (1.0_EB-p_i)
                   U_i_Hawk = p_i * T_0_HAWK/(T_tmp-(T_ASET_HAWK+T_ASET_TFAC_HAWK*MAX(0.0_EB,T-A_HAWK_T_START))+T_0_HAWK)&
                        - (1.0_EB-p_i)
                ELSE
                   IF (EPSILON_HAWK <= 0.0_EB) THEN
                      U_i_Hawk = 10000000000000.0_EB ! all doves close to the door
                   ELSE
                      U_i_Hawk = p_i*EPSILON_HAWK - (1.0_EB-p_i)
                   END IF
                END IF
                U_i_Dove = p_i
!!$                IF (HR%A <= 500.0_EB .AND. (U_i_Dove < U_i_Hawk .OR. 2.0_EB*Delta_Hawk*Alpha_Hawk* &
!!$                     (1.0_EB+A_HAWK_T_FACTOR*MAX(0.0_EB,T-A_HAWK_T_START))*T_tmp < Epsilon_Hawk)) THEN   ! Was a hawk
                IF (HR%A <= 500.0_EB .AND. (U_i_Dove < U_i_Hawk)) THEN   ! Was a hawk
                   ! DoseCrit1 = hr%a   modified, no sharp transition
                   ! DoseCrit2 = hr%tau modified, now sharp transition
                   HR%DoseCrit1 = MAX(HR%A,  HR%DoseCrit1) ! Starting value for hawk => dove transition
                   HR%DoseCrit2 = MAX(HR%TAU,HR%DoseCrit2) ! Starting value for hawk => dove transition
                   HR%A = 10.0_EB*HR%A ; HR%TAU = 2.0_EB*HR%TAU
                   HR%DoseCrit2 = HR%TAU ! sharp transition
                END IF
!!$                IF (HR%A >  500.0_EB .AND. U_i_Hawk < U_i_Dove .AND. 2.0_EB*Delta_Hawk*Alpha_Hawk* &
!!$                     (1.0_EB+A_HAWK_T_FACTOR*MAX(0.0_EB,T-A_HAWK_T_START))*T_tmp >= Epsilon_Hawk) THEN   ! Was a dove
                IF (HR%A >  500.0_EB .AND. (U_i_Hawk < U_i_Dove) ) THEN   ! Was a dove
                   ! DoseCrit1 = hr%a   modified, no sharp transition
                   ! DoseCrit2 = hr%tau modified, now sharp transition
                   HR%DoseCrit1 = MIN(HR%A,  HR%DoseCrit1) ! Starting value for hawk => dove transition
                   HR%DoseCrit2 = MIN(HR%TAU,HR%DoseCrit2) ! Starting value for hawk => dove transition
                   HR%A =  0.1_EB*HR%A ; HR%TAU = 0.5_EB*HR%TAU
                   HR%DoseCrit2 = HR%TAU ! sharp transition
                END IF
             ELSE
                ! No neighbors, set to a dove if it is a hawk
                ! It is better looking if the status is not changed.u_i_
                !IF (HR%A <= 500.0_EB) THEN
                !   HR%A = 10.0_EB*HR%A ; HR%TAU = 2.0_EB*HR%TAU
                !END IF
             END IF
          END IF Hawk_Dove_Game2
          ! ========================================================
          ! Person-person interaction forces ends here
          ! ========================================================

          ! ========================================================
          ! The person-wall forces
          ! ========================================================
          !
          ! Walls are looked for the body circle
          CALL FIND_WALLS(NM, X1, Y1, HR%RADIUS, P2P_DIST_MAX, HR%SKIP_WALL_FORCE_IOR, D_XY, FOUNDWALL_XY, ISTAT)
          IF (STOP_STATUS>NO_STOP) THEN
             STOP_STATUS=EVACUATION_STOP
             RETURN
          END IF
          ! ========================================================
          ! Collision avoidance, counterflow, etc.
          ! ========================================================
          CHANGE_V0_RNCF2: IF ( RNCF > EXP(-DTSP/TAU_CHANGE_V0) ) THEN
             V_HR  = MAX(0.1_EB,SQRT(HR%U**2+HR%V**2)/HR_SPEED)
             TIM_DIST = 0.0_EB
             DO III = 1, N_SECTORS
                ! Avoid walls, do not take a direction where there is a wall closer than
                ! D_PERP = 0.6 m (perpendicular).
                IF (ABS(SIN_THETA(III)) >= TWO_EPSILON_EB) THEN
                   X11 = HR%X + U_THETA(III)*MIN(P2P_SUUNTA_MAX, 0.6_EB/ABS(SIN_THETA(III)))
                   Y11 = HR%Y + V_THETA(III)*MIN(P2P_SUUNTA_MAX, 0.6_EB/ABS(SIN_THETA(III)))
                ELSE ! Straight ahead
                   X11 = HR%X + U_THETA(III)*P2P_SUUNTA_MAX
                   Y11 = HR%Y + V_THETA(III)*P2P_SUUNTA_MAX
                END IF
                P2P_DIST = SQRT((HR%X-X11)**2 + (HR%Y-Y11)**2) - HR%RADIUS
                X11 = MAX(XS,MIN(XF,X11))
                Y11 = MAX(YS,MIN(YF,Y11))
                PP_SEE_EACH = SEE_EACH_OTHER(NM, HR%X, HR%Y, X11, Y11)
                IF(.NOT.PP_SEE_EACH) THEN
                   VR_2R = -FAC_1_WALL*P2P_SUUNTA_MAX*V_HR
                   SUM_SUUNTA(III) = SUM_SUUNTA(III) + VR_2R/MAX(0.2_EB,P2P_DIST)
                END IF
             END DO

             IF (FOUNDWALL_XY(1)) TIM_DIST = TIM_DIST - 3.0_EB/MAX(0.3_EB,ABS(HR%X-D_XY(1))-HR%RADIUS)
             IF (FOUNDWALL_XY(2)) TIM_DIST = TIM_DIST - 3.0_EB/MAX(0.3_EB,ABS(HR%X-D_XY(2))-HR%RADIUS)
             IF (FOUNDWALL_XY(3)) TIM_DIST = TIM_DIST - 3.0_EB/MAX(0.3_EB,ABS(HR%Y-D_XY(3))-HR%RADIUS)
             IF (FOUNDWALL_XY(4)) TIM_DIST = TIM_DIST - 3.0_EB/MAX(0.3_EB,ABS(HR%Y-D_XY(4))-HR%RADIUS)

             IF (FOUNDWALL_XY(1) .AND. ABS(HR%X-D_XY(1))-HR%RADIUS < 0.1_EB) THEN
                DO III = 1, N_SECTORS
                   IF(U_THETA(III) < -0.10_EB) THEN
                      VR_2R = U_THETA(III)/MAX(0.1_EB, ABS(HR%X-HR%RADIUS-D_XY(1)))
                      SUM_SUUNTA(III) = SUM_SUUNTA(III) - ABS(P2P_SUUNTA_MAX*FAC_2_WALL*VR_2R)
                   END IF
                END DO
             END IF
             IF (FOUNDWALL_XY(2) .AND. ABS(HR%X-D_XY(2))-HR%RADIUS < 0.1_EB) THEN
                DO III = 1, N_SECTORS
                   IF(U_THETA(III) > +0.10_EB) THEN
                      VR_2R = -U_THETA(III)/MAX(0.1_EB, ABS(HR%X+HR%RADIUS-D_XY(2)))
                      SUM_SUUNTA(III) = SUM_SUUNTA(III) - ABS(P2P_SUUNTA_MAX*FAC_2_WALL*VR_2R)
                   END IF
                END DO
             END IF
             IF (FOUNDWALL_XY(3) .AND. ABS(HR%Y-D_XY(3))-HR%RADIUS < 0.1_EB) THEN
                DO III = 1, N_SECTORS
                   IF(V_THETA(III) < -0.10_EB) THEN
                      VR_2R = V_THETA(III)/MAX(0.1_EB, ABS(HR%Y-HR%RADIUS-D_XY(3)))
                      SUM_SUUNTA(III) = SUM_SUUNTA(III) - ABS(P2P_SUUNTA_MAX*FAC_2_WALL*VR_2R)
                   END IF
                END DO
             END IF
             IF (FOUNDWALL_XY(4) .AND. ABS(HR%Y-D_XY(4))-HR%RADIUS < 0.1_EB) THEN
                DO III = 1, N_SECTORS
                   IF(V_THETA(III) > +0.10_EB) THEN
                      VR_2R = -V_THETA(III)/MAX(0.1_EB, ABS(HR%Y+HR%RADIUS-D_XY(4)))
                      SUM_SUUNTA(III) = SUM_SUUNTA(III) - ABS(P2P_SUUNTA_MAX*FAC_2_WALL*VR_2R)
                   END IF
                END DO
             END IF
             SUM_SUUNTA(N_SECTORS+1) = SUM_SUUNTA(N_SECTORS+1) + ABS(FAC_V0_DIR)*V_HR
             DO III = 1, N_SECTORS
                IF (N_SUUNTACF(N_SECTORS+1) < 1) THEN
                   ! No counterflow: prefer left (and straight ahead)
                ELSE
                   ! If there is at least one counterflow agent in the front sector
                   ! Counterflow: prefer right (and straight ahead)
                   SUM_SUUNTA(III) = SUM_SUUNTA(III) - SIGN(1.0_EB,THETAS(III))*FAC_V0_DIR
                END IF
             END DO

             ! If there is no agents in the front sector, do not use counterflow algorithm
             IF (N_SUUNTA(N_SECTORS+1) < 1) SUM_SUUNTA(N_SECTORS+1) = SUM_SUUNTA(N_SECTORS+1) + 50.0_EB  ! Empty space ahead
             IF ((N_SUUNTACF(N_SECTORS+1)) < 1) THEN
                ! No counterflow in the front sector, prefer v0 direction, i.e., "stay on line"
                SUM_SUUNTA(N_SECTORS+1) = SUM_SUUNTA(N_SECTORS+1) + &
                     (FAC_NOCF + FAC_V0_NOCF*V_HR)*(N_SUUNTA(N_SECTORS+1))
             END IF
             I_SUUNTA_MAX = N_SECTORS + 1
             SUM_SUUNTA_MAX = -HUGE(SUM_SUUNTA_MAX)
             DO III = 1, N_SECTORS+1
                IF (SUM_SUUNTA(III) > SUM_SUUNTA_MAX) THEN
                   UBAR = U_THETA(III)
                   VBAR = V_THETA(III)
                   SUM_SUUNTA_MAX = SUM_SUUNTA(III)
                   I_SUUNTA_MAX = III
                END IF
             END DO
             ! If counterflow then try to pass sideways and decrease TAU and TAU_INER.
             ! Use this if there are more counterflow agents than downstrean agents.
             ANGLE_OLD = 0.0_EB
             COMMITMENT = 0.0_EB
             IF (V_HR < 0.3_EB .AND. (N_SUUNTA(N_SECTORS+1)-N_SUUNTACF(N_SECTORS+1)) < &
                  N_SUUNTACF(N_SECTORS+1)-0) THEN
                COMMITMENT = MAX(0.5_EB,REAL(SUM(N_SUUNTACF(1:N_SECTORS+1)),EB)/MAX(1,SUM(N_SUUNTA(1:N_SECTORS+1))))
                ANGLE_OLD = -SIGN(1.0_EB,THETAS(I_SUUNTA_MAX))*MAX(0.1_EB,ABS(THETAS(I_SUUNTA_MAX)))*PI/180.0_EB
                ANGLE_OLD = ANGLE_OLD + 85.0_EB*PI/180.0_EB
             ELSE IF (V_HR < 0.3_EB .AND. SUM(N_SUUNTACF(1:N_SECTORS)) >= 1 .AND. TIM_DIST < -15.0_EB) THEN
                COMMITMENT = MAX(0.5_EB,REAL(SUM(N_SUUNTACF(1:N_SECTORS)),EB)/MAX(1,SUM(N_SUUNTA(1:N_SECTORS))))
                ANGLE_OLD = -SIGN(1.0_EB,THETAS(I_SUUNTA_MAX))*MAX(0.1_EB,ABS(THETAS(I_SUUNTA_MAX)))*PI/180.0_EB
                ANGLE_OLD = ANGLE_OLD + 85.0_EB*PI/180.0_EB
             ELSE
                ANGLE_OLD = 0.0_EB
                COMMITMENT = 0.0_EB
             END IF
          ELSE ! CHANGE_V0_RNCF2
             ! Do not change direction during this time step, use the previous direction
             EVEL = UBAR**2 + VBAR**2
             IF (TAU_CHANGE_V0 > 1.0E-12_EB .AND. EVEL >= TWO_EPSILON_EB) THEN
                UBAR = HR%UBAR
                VBAR = HR%VBAR
                ANGLE_OLD = HR%ANGLE_OLD
                COMMITMENT = HR%COMMITMENT
             ELSE
                HR%UBAR = UBAR
                HR%VBAR = VBAR
                ANGLE_OLD = 0.0_EB
                COMMITMENT = 0.0_EB
             END IF
          END IF CHANGE_V0_RNCF2 ! COLLISION AVOIDANCE ENDS
          HR%ANGLE_OLD = ANGLE_OLD
          HR%COMMITMENT = COMMITMENT

          IF (.NOT.L_FALLEN_DOWN .AND. .NOT.L_DEAD) &
               CALL WALL_SOCIALFORCES(X_TMP, Y_TMP, R_TMP, P2P_DIST_MAX, D_XY, P2P_U, P2P_V, SOCIAL_F, FOUNDWALL_XY)

          CALL WALL_CONTACTFORCES(X_TMP(1), Y_TMP(1), R_TMP(1), U_TMP(1), V_TMP(1), D_XY, &
               P2P_U, P2P_V, P2P_TORQUE, CONTACT_F, D_WALLS, FOUNDWALL_XY, CONTACT_FX, CONTACT_FY)
          CALL WALL_CONTACTFORCES(X_TMP(2), Y_TMP(2), R_TMP(2), U_TMP(2), V_TMP(2), D_XY, &
               P2P_U, P2P_V, P2P_TORQUE, CONTACT_F, D_WALLS, FOUNDWALL_XY, CONTACT_FX, CONTACT_FY)
          CALL WALL_CONTACTFORCES(X_TMP(3), Y_TMP(3), R_TMP(3), U_TMP(3), V_TMP(3), D_XY, &
               P2P_U, P2P_V, P2P_TORQUE, CONTACT_F, D_WALLS, FOUNDWALL_XY, CONTACT_FX, CONTACT_FY)
          IF (STOP_STATUS>NO_STOP) THEN
             STOP_STATUS=EVACUATION_STOP
             RETURN
          END IF

          ! Add forces from the door case
          CALL DOOR_FORCES(X_TMP, Y_TMP, R_TMP, U_TMP, V_TMP, P2P_DIST_MAX, D_XY,&
               P2P_U, P2P_V, SOCIAL_F, CONTACT_F, P2P_TORQUE, FOUNDWALL_XY, CONTACT_FX, CONTACT_FY)
          IF (STOP_STATUS>NO_STOP) THEN
             STOP_STATUS=EVACUATION_STOP
             RETURN
          END IF

          IF (NM_STRS_MESH .AND. ABS(HR%SKIP_WALL_FORCE_IOR)>0 .AND. I_STRS_DOOR>0) THEN
             IF (I_STRS_DOOR > N_DOORS) THEN
                X11 = EVAC_EXITS(I_STRS_DOOR-N_DOORS)%X1
                Y11 = EVAC_EXITS(I_STRS_DOOR-N_DOORS)%Y1
             ELSE
                X11 = EVAC_DOORS(I_STRS_DOOR)%X1
                Y11 = EVAC_DOORS(I_STRS_DOOR)%Y1
             END IF
             ! x1, y1 agent,  x11, y11, coordiantes of the corner
             CALL CORNER_FORCES(X1, Y1, X11, Y11, P2P_DIST_MAX, P2P_U, P2P_V, SOCIAL_F, &
                  CONTACT_F, P2P_TORQUE, D_WALLS, X_TMP, Y_TMP, R_TMP, U_TMP, V_TMP, ISTAT, CONTACT_FX, CONTACT_FY)
             IF (I_STRS_DOOR > N_DOORS) THEN
                X11 = EVAC_EXITS(I_STRS_DOOR-N_DOORS)%X2
                Y11 = EVAC_EXITS(I_STRS_DOOR-N_DOORS)%Y2
             ELSE
                X11 = EVAC_DOORS(I_STRS_DOOR)%X2
                Y11 = EVAC_DOORS(I_STRS_DOOR)%Y2
             END IF
             ! x1, y1 agent,  x11, y11, coordiantes of the corner
             CALL CORNER_FORCES(X1, Y1, X11, Y11, P2P_DIST_MAX, P2P_U, P2P_V, SOCIAL_F, &
                  CONTACT_F, P2P_TORQUE, D_WALLS, X_TMP, Y_TMP, R_TMP, U_TMP, V_TMP, ISTAT, CONTACT_FX, CONTACT_FY)
             IF (STOP_STATUS>NO_STOP) THEN
                STOP_STATUS=EVACUATION_STOP
                RETURN
             END IF
          END IF

          ! Add wall corner - person forces
          ! Top right corner (X > X_HUMAN, Y > Y_HUMAN)
          X_NOW = -DX(IIN+1)
          LOOP_PX: DO II = IIN, IBAR, +1
             X_NOW = X_NOW + DX(II)
             IF (X_NOW-HR%RADIUS > P2P_DIST_MAX) EXIT LOOP_PX
             Y_NOW = -DY(JJN-1)
             LOOP_PXPY: DO JJ = JJN, JBAR
                Y_NOW = Y_NOW + DY(JJ)
                IF (SQRT(X_NOW**2 + Y_NOW**2)-HR%RADIUS > P2P_DIST_MAX) EXIT LOOP_PXPY
                TIM_IC  = CELL_INDEX(II,JJ,KKN)   ! Present
                TIM_IWX = WALL_INDEX(TIM_IC, +1)  ! Right
                IF (TIM_IWX>0) THEN
                   I_OBST = WALL(TIM_IWX)%OBST_INDEX
                   IF (OBSTRUCTION(I_OBST)%HIDDEN) TIM_IWX = 0
                END IF
                TIM_IWY = WALL_INDEX(TIM_IC, +2)  ! Up
                IF (TIM_IWY>0) THEN
                   I_OBST = WALL(TIM_IWY)%OBST_INDEX
                   IF (OBSTRUCTION(I_OBST)%HIDDEN) TIM_IWY = 0
                END IF
                TIM_IC  = CELL_INDEX(II,JJ+1,KKN) ! One cell up
                TIM_IW  = WALL_INDEX(TIM_IC, +1)  ! Up and right
                IF (TIM_IW>0) THEN
                   I_OBST = WALL(TIM_IW)%OBST_INDEX
                   IF (OBSTRUCTION(I_OBST)%HIDDEN) TIM_IW = 0
                END IF
                TIM_IC2 = CELL_INDEX(II+1,JJ,KKN) ! One cell right
                TIM_IW2 = WALL_INDEX(TIM_IC2, +2) ! Right and up
                IF (TIM_IW2>0) THEN
                   I_OBST = WALL(TIM_IW2)%OBST_INDEX
                   IF (OBSTRUCTION(I_OBST)%HIDDEN) TIM_IW2 = 0
                END IF
                IF (TIM_IWY /= 0) EXIT LOOP_PXPY  !
                IF ( (TIM_IWX==0).AND.(TIM_IWY==0).AND.(TIM_IW/=0 .OR. TIM_IW2/=0) ) THEN
                   IF (TIM_IW/=0) THEN
                      ! First y-direction then x-direction
                      X11 = WALL(TIM_IW )%X                 ! Corner point x
                      Y11 = WALL(TIM_IW )%Y-0.5_EB*DY(JJ+1) ! Corner point y
                   ELSE
                      ! First x-direction then y-direction
                      X11 = WALL(TIM_IW2)%X-0.5_EB*DX(II+1) ! Corner point x
                      Y11 = WALL(TIM_IW2)%Y                 ! Corner point y
                   END IF

                   CALL CORNER_FORCES(X1, Y1, X11, Y11, P2P_DIST_MAX, P2P_U, P2P_V, SOCIAL_F, &
                        CONTACT_F, P2P_TORQUE, D_WALLS, X_TMP, Y_TMP, R_TMP, U_TMP, V_TMP, ISTAT, CONTACT_FX, CONTACT_FY)
                   EXIT LOOP_PXPY
                END IF
             END DO LOOP_PXPY
          END DO LOOP_PX

          ! Top left corner (X < X_HUMAN, Y > Y_HUMAN)
          X_NOW = -DX(IIN-1)
          LOOP_MX: DO II = IIN, 1, -1
             X_NOW = X_NOW + DX(II)
             IF (X_NOW-HR%RADIUS > P2P_DIST_MAX) EXIT LOOP_MX
             Y_NOW = -DY(JJN)
             LOOP_MXPY: DO JJ = JJN, JBAR
                Y_NOW = Y_NOW + DY(JJ)
                IF (SQRT(X_NOW**2 + Y_NOW**2)-HR%RADIUS > P2P_DIST_MAX) EXIT LOOP_MXPY
                TIM_IC  = CELL_INDEX(II,JJ,KKN)   ! Present
                TIM_IWX = WALL_INDEX(TIM_IC, -1)  ! Left
                IF (TIM_IWX>0) THEN
                   I_OBST = WALL(TIM_IWX)%OBST_INDEX
                   IF (OBSTRUCTION(I_OBST)%HIDDEN) TIM_IWX = 0
                END IF
                TIM_IWY = WALL_INDEX(TIM_IC, +2)  ! Up
                IF (TIM_IWY>0) THEN
                   I_OBST = WALL(TIM_IWY)%OBST_INDEX
                   IF (OBSTRUCTION(I_OBST)%HIDDEN) TIM_IWY = 0
                END IF
                TIM_IC  = CELL_INDEX(II,JJ+1,KKN) ! One cell up
                TIM_IW  = WALL_INDEX(TIM_IC, -1)  ! Up and left
                IF (TIM_IW>0) THEN
                   I_OBST = WALL(TIM_IW)%OBST_INDEX
                   IF (OBSTRUCTION(I_OBST)%HIDDEN) TIM_IW = 0
                END IF
                TIM_IC2 = CELL_INDEX(II-1,JJ,KKN) ! One cell left
                TIM_IW2 = WALL_INDEX(TIM_IC2, +2) ! Left and up
                IF (TIM_IW2>0) THEN
                   I_OBST = WALL(TIM_IW2)%OBST_INDEX
                   IF (OBSTRUCTION(I_OBST)%HIDDEN) TIM_IW2 = 0
                END IF
                IF (TIM_IWY /= 0) EXIT LOOP_MXPY
                IF ( (TIM_IWX==0).AND.(TIM_IWY==0).AND.(TIM_IW/=0 .OR. TIM_IW2/=0) ) THEN
                   IF (TIM_IW/=0) THEN
                      ! First y-direction then x-direction
                      X11 = WALL(TIM_IW )%X                 ! Corner point x
                      Y11 = WALL(TIM_IW )%Y-0.5_EB*DY(JJ+1) ! Corner point y
                   ELSE
                      ! First x-direction then y-direction
                      X11 = WALL(TIM_IW2)%X+0.5_EB*DX(II-1) ! Corner point x
                      Y11 = WALL(TIM_IW2)%Y                 ! Corner point y
                   END IF

                   CALL CORNER_FORCES(X1, Y1, X11, Y11, P2P_DIST_MAX, P2P_U, P2P_V, SOCIAL_F, &
                        CONTACT_F, P2P_TORQUE, D_WALLS, X_TMP, Y_TMP, R_TMP, U_TMP, V_TMP, ISTAT, CONTACT_FX, CONTACT_FY)
                   EXIT LOOP_MXPY
                END IF
             END DO LOOP_MXPY
          END DO LOOP_MX

          ! Bottom right corner (X > X_HUMAN, Y < Y_HUMAN)
          X_NOW = -DX(IIN+1)
          LOOP_PY: DO II = IIN, IBAR, +1
             X_NOW = X_NOW + DX(II)
             IF (X_NOW-HR%RADIUS > P2P_DIST_MAX) EXIT LOOP_PY
             Y_NOW = -DY(JJN-1)
             LOOP_PXMY: DO JJ = JJN, 1, -1
                Y_NOW = Y_NOW + DY(JJ)
                IF (SQRT(X_NOW**2 + Y_NOW**2)-HR%RADIUS > P2P_DIST_MAX) EXIT LOOP_PXMY
                TIM_IC  = CELL_INDEX(II,JJ,KKN)   ! Present
                TIM_IWX = WALL_INDEX(TIM_IC, +1)  ! Right
                IF (TIM_IWX>0) THEN
                   I_OBST = WALL(TIM_IWX)%OBST_INDEX
                   IF (OBSTRUCTION(I_OBST)%HIDDEN) TIM_IWX = 0
                END IF
                TIM_IWY = WALL_INDEX(TIM_IC, -2)  ! Down
                IF (TIM_IWY>0) THEN
                   I_OBST = WALL(TIM_IWY)%OBST_INDEX
                   IF (OBSTRUCTION(I_OBST)%HIDDEN) TIM_IWY = 0
                END IF
                TIM_IC  = CELL_INDEX(II,JJ-1,KKN) ! One cell down
                TIM_IW  = WALL_INDEX(TIM_IC, +1)  ! Down and right
                IF (TIM_IW>0) THEN
                   I_OBST = WALL(TIM_IW)%OBST_INDEX
                   IF (OBSTRUCTION(I_OBST)%HIDDEN) TIM_IW = 0
                END IF
                TIM_IC2 = CELL_INDEX(II+1,JJ,KKN) ! One cell right
                TIM_IW2 = WALL_INDEX(TIM_IC2, -2) ! Right and down
                IF (TIM_IW2>0) THEN
                   I_OBST = WALL(TIM_IW2)%OBST_INDEX
                   IF (OBSTRUCTION(I_OBST)%HIDDEN) TIM_IW2 = 0
                END IF
                IF (TIM_IWY /= 0) EXIT LOOP_PXMY
                IF ( (TIM_IWX==0).AND.(TIM_IWY==0).AND.(TIM_IW/=0 .OR. TIM_IW2/=0) ) THEN
                   IF (TIM_IW/=0) THEN
                      ! First y-direction then x-direction
                      X11 = WALL(TIM_IW )%X                 ! Corner point x
                      Y11 = WALL(TIM_IW )%Y+0.5_EB*DY(JJ-1) ! Corner point y
                   ELSE
                      ! First x-direction then y-direction
                      X11 = WALL(TIM_IW2)%X-0.5_EB*DX(II+1) ! Corner point x
                      Y11 = WALL(TIM_IW2)%Y                 ! Corner point y
                   END IF

                   CALL CORNER_FORCES(X1, Y1, X11, Y11, P2P_DIST_MAX, P2P_U, P2P_V, SOCIAL_F, &
                        CONTACT_F, P2P_TORQUE, D_WALLS, X_TMP, Y_TMP, R_TMP, U_TMP, V_TMP, ISTAT, CONTACT_FX, CONTACT_FY)
                   EXIT LOOP_PXMY
                END IF
             END DO LOOP_PXMY
          END DO LOOP_PY

          ! Bottom left corner (X < X_HUMAN, Y < Y_HUMAN)
          X_NOW = -DX(IIN-1)
          LOOP_MY: DO II = IIN, 1, -1
             X_NOW = X_NOW + DX(II)
             IF (X_NOW-HR%RADIUS > P2P_DIST_MAX) EXIT LOOP_MY
             Y_NOW = -DY(JJN-1)
             LOOP_MXMY: DO JJ = JJN, 1, -1
                Y_NOW = Y_NOW + DY(JJ)
                IF (SQRT(X_NOW**2+Y_NOW**2)-HR%RADIUS > P2P_DIST_MAX) EXIT LOOP_MXMY
                TIM_IC  = CELL_INDEX(II,JJ,KKN)   ! Present
                TIM_IWX = WALL_INDEX(TIM_IC, -1)  ! Left
                IF (TIM_IWX>0) THEN
                   I_OBST = WALL(TIM_IWX)%OBST_INDEX
                   IF (OBSTRUCTION(I_OBST)%HIDDEN) TIM_IWX = 0
                END IF
                TIM_IWY = WALL_INDEX(TIM_IC, -2)  ! Down
                IF (TIM_IWY>0) THEN
                   I_OBST = WALL(TIM_IWY)%OBST_INDEX
                   IF (OBSTRUCTION(I_OBST)%HIDDEN) TIM_IWY = 0
                END IF
                TIM_IC  = CELL_INDEX(II,JJ-1,KKN) ! Once cell down
                TIM_IW  = WALL_INDEX(TIM_IC, -1)  ! Down and left
                IF (TIM_IW>0) THEN
                   I_OBST = WALL(TIM_IW)%OBST_INDEX
                   IF (OBSTRUCTION(I_OBST)%HIDDEN) TIM_IW = 0
                END IF
                TIM_IC2 = CELL_INDEX(II-1,JJ,KKN) ! One cell left
                TIM_IW2 = WALL_INDEX(TIM_IC2, -2) ! Left and down
                IF (TIM_IW2>0) THEN
                   I_OBST = WALL(TIM_IW2)%OBST_INDEX
                   IF (OBSTRUCTION(I_OBST)%HIDDEN) TIM_IW2 = 0
                END IF
                IF (TIM_IWY /= 0) EXIT LOOP_MXMY
                IF ( (TIM_IWX==0).AND.(TIM_IWY==0).AND.(TIM_IW/=0 .OR. TIM_IW2/=0) ) THEN
                   IF (TIM_IW/=0) THEN
                      ! First y-direction then x-direction
                      X11 = WALL(TIM_IW )%X                 ! Corner point x
                      Y11 = WALL(TIM_IW )%Y+0.5_EB*DY(JJ-1) ! Corner point y
                   ELSE
                      ! First x-direction then y-direction
                      X11 = WALL(TIM_IW2)%X+0.5_EB*DX(II-1) ! Corner point x
                      Y11 = WALL(TIM_IW2)%Y                 ! Corner point y
                   END IF

                   CALL CORNER_FORCES(X1, Y1, X11, Y11, P2P_DIST_MAX, P2P_U, P2P_V, SOCIAL_F, &
                        CONTACT_F, P2P_TORQUE, D_WALLS, X_TMP, Y_TMP, R_TMP, U_TMP, V_TMP, ISTAT, CONTACT_FX, CONTACT_FY)
                   EXIT LOOP_MXMY
                END IF
             END DO LOOP_MXMY
          END DO LOOP_MY
          IF (STOP_STATUS>NO_STOP) THEN
             STOP_STATUS=EVACUATION_STOP
             RETURN
          END IF

          ! ========================================================
          ! The person-wall forces ends here
          ! ========================================================

          ! Save the forces for the loop evac_move_loop (next time step)
          HR%F_X = P2P_U
          HR%F_Y = P2P_V
          HR%SUMFORCES  = CONTACT_F
          HR%SUMFORCES2 = SOCIAL_F + CONTACT_F
          HR%TORQUE = P2P_TORQUE

          ! Fall down algorithm: Check the resultant contact force
          IF (TAU_FALL_DOWN >= 0.0_EB .AND. T > HR%TDET) THEN
             IF (.NOT.L_FALLEN_DOWN .AND. SQRT(CONTACT_FX**2+CONTACT_FY**2) > HR%F_FallDown) THEN
                L_FALLEN_DOWN = .TRUE.
                N_DEAD = N_DEAD+1
                HR%T_FallenDown = T
                HR%Angle_FallenDown = HR%angle
                HR%SizeFac_FallenDown = 0.0_EB
             END IF
          END IF

          IF ( T <= T_BEGIN ) THEN
             IF ( ABS(P2P_U)/HR%MASS > 550.0_EB ) P2P_U =550.0_EB*HR%MASS*P2P_U/ABS(P2P_U)
             IF ( ABS(P2P_V)/HR%MASS > 550.0_EB ) P2P_V =550.0_EB*HR%MASS*P2P_V/ABS(P2P_V)
             HR%F_X = P2P_U
             HR%F_Y = P2P_V
          END IF

          ! Add the social and contact force terms
          U_NEW     = HR%U + 0.5_EB*HR%F_X*DTSP/HR%MASS
          V_NEW     = HR%V + 0.5_EB*HR%F_Y*DTSP/HR%MASS
          OMEGA_NEW = HR%OMEGA + 0.5_EB*DTSP*HR%TORQUE/HR%M_INER

          ! Add the effect of the random force
          IF (GATH >= TWO_EPSILON_EB .AND. T > T_BEGIN ) THEN
             U_NEW = U_NEW + 0.5_EB*DTSP*HR%V0_FAC*HR%MASS*HR%KSI*COS(HR%ETA)/HR%MASS
             V_NEW = V_NEW + 0.5_EB*DTSP*HR%V0_FAC*HR%MASS*HR%KSI*SIN(HR%ETA)/HR%MASS
             P2P_U = P2P_U + HR%V0_FAC*HR%MASS*HR%KSI*COS(HR%ETA)
             P2P_V = P2P_V + HR%V0_FAC*HR%MASS*HR%KSI*SIN(HR%ETA)
             OMEGA_NEW = OMEGA_NEW + 0.5_EB*DTSP*1.0_EB*SIGN(HR%KSI,HR%ETA-PI)
             P2P_TORQUE = P2P_TORQUE + 1.0_EB*SIGN(HR%KSI,HR%ETA-PI)*HR%M_INER
          ELSE
             HR%KSI = 0.0_EB
             HR%ETA = 0.0_EB
          END IF
          ! Now the step (4a) of SC-VV by VATTULAINEN is ended.

          ! Add self-propelling force terms, self-consistent VV
          ! (First time step towards the exit door)
          IF ( T <= TPRE ) THEN
             HR_TAU = MAX(CF_MIN_TAU, 0.1_EB*HR%TAU)
             IF ( (T+DTSP_NEW) > TPRE) THEN
                EVEL = SQRT(UBAR**2 + VBAR**2)
                IF (EVEL >= TWO_EPSILON_EB) THEN
                   SPEED = SPEED_YP*(0.5_EB + SIGN(0.5_EB,VBAR)) + SPEED_YM*(0.5_EB - SIGN(0.5_EB,VBAR))
                   SPEED = SPEED*HR%V0_FAC
                   TAU_FAC = (SPEED*(VBAR/EVEL))**2
                   SPEED = SPEED_XP*(0.5_EB + SIGN(0.5_EB,UBAR)) + SPEED_XM*(0.5_EB - SIGN(0.5_EB,UBAR))
                   SPEED = SPEED*HR%V0_FAC
                   TAU_FAC = MAX(0.1_EB,SQRT(TAU_FAC + (SPEED*(UBAR/EVEL))**2)/MAX(0.1_EB,HR_SPEED))
                   HR_TAU=HR_TAU*TAU_FAC
                   P2P_U = P2P_U + (HR%MASS/HR_TAU)*SPEED*(UBAR/EVEL)
                   SPEED = SPEED_YP*(0.5_EB + SIGN(0.5_EB,VBAR)) + SPEED_YM*(0.5_EB - SIGN(0.5_EB,VBAR))
                   SPEED = SPEED*HR%V0_FAC
                   P2P_V = P2P_V + (HR%MASS/HR_TAU)*SPEED*(VBAR/EVEL)
                END IF
             END IF
             IF (.NOT. HR%CROWBAR_READ_IN) THEN
                UBAR = 0.0_EB
                VBAR = 0.0_EB
             END IF
          END IF

          EVEL = SQRT(UBAR**2 + VBAR**2)
          IF (T > TPRE .AND. EVEL >= TWO_EPSILON_EB) THEN
             SPEED = SPEED_YP*(0.5_EB + SIGN(0.5_EB,VBAR)) + SPEED_YM*(0.5_EB - SIGN(0.5_EB,VBAR))
             SPEED = SPEED*HR%V0_FAC
             TAU_FAC = (SPEED*(VBAR/EVEL))**2
             SPEED = SPEED_XP*(0.5_EB + SIGN(0.5_EB,UBAR)) + SPEED_XM*(0.5_EB - SIGN(0.5_EB,UBAR))
             SPEED = SPEED*HR%V0_FAC
             TAU_FAC = MAX(0.1_EB,SQRT(TAU_FAC + (SPEED*(UBAR/EVEL))**2)/MAX(0.1_EB,HR_SPEED))
             HR_TAU=HR_TAU*TAU_FAC
          END IF

          FAC_TIM =  1.0_EB + (DTSP/(2.0_EB*HR_TAU))

          ! Add self-propelling force terms, self-consistent VV
          EVEL = SQRT(UBAR**2 + VBAR**2)
          IF (EVEL >= TWO_EPSILON_EB) THEN
             SPEED = SPEED_XP*(0.5_EB + SIGN(0.5_EB,UBAR)) + SPEED_XM*(0.5_EB - SIGN(0.5_EB,UBAR))
             SPEED = SPEED*HR%V0_FAC
             SPEED_X = SPEED
             U_NEW = (U_NEW + 0.5_EB*(DTSP/HR_TAU)*SPEED*(UBAR/EVEL)) / FAC_TIM
             HR%UBAR = UBAR
             P2P_U = P2P_U + (HR%MASS/HR_TAU)*(SPEED*(UBAR/EVEL) - HR%U)

             SPEED = SPEED_YP*(0.5_EB + SIGN(0.5_EB,VBAR)) + SPEED_YM*(0.5_EB - SIGN(0.5_EB,VBAR))
             SPEED = SPEED*HR%V0_FAC
             SPEED_Y = SPEED
             V_NEW = (V_NEW + 0.5_EB*(DTSP/HR_TAU)*SPEED*(VBAR/EVEL)) / FAC_TIM
             HR%VBAR = VBAR
             P2P_V = P2P_V + (HR%MASS/HR_TAU)*(SPEED*(VBAR/EVEL) - HR%V)

             IF (VBAR >= 0.0_EB) THEN
                ANGLE = ACOS(UBAR/EVEL)
             ELSE
                ANGLE = 2.0_EB*PI - ACOS(UBAR/EVEL)
             END IF

             ! Collision avoidance has HR%ANGLE_OLD /= 0.0
             ANGLE = ANGLE + HR%ANGLE_OLD
             DO WHILE (ANGLE >= 2.0_EB*PI)
                ANGLE = ANGLE - 2.0_EB*PI
             END DO
             DO WHILE (ANGLE < 0.0_EB)
                ANGLE = ANGLE + 2.0_EB*PI
             END DO

             IF (ANGLE == 2.0_EB*PI) ANGLE = 0.0_EB  ! ANGLE IS [0,2PI)

             ! Rotational motion: J(DW/DT) = (J/T_INER)*( ((ANGLE-ANGLE_0)/PI)*W_0 - W)
             IF (ABS(ANGLE-HR%ANGLE) <= PI ) THEN
                ! Zero is not crossed.
                OMEGA_NEW = OMEGA_NEW + 0.5_EB*(DTSP/HR_TAU_INER)*((ANGLE-HR%ANGLE)*(OMEGA_0/PI) - HR%OMEGA)
                P2P_TORQUE = P2P_TORQUE + (HR%M_INER/HR_TAU_INER)*((ANGLE-HR%ANGLE)*(OMEGA_0/PI) - HR%OMEGA)
             ELSE
                ! Zero is crossed
                OMEGA_NEW = OMEGA_NEW + 0.5_EB*(DTSP/HR_TAU_INER)* &
                     ( (2.0_EB*PI-ABS(ANGLE-HR%ANGLE))*SIGN(1.0_EB , HR%ANGLE-ANGLE)*(OMEGA_0/PI) - HR%OMEGA)
                P2P_TORQUE = P2P_TORQUE + (HR%M_INER/HR_TAU_INER)* &
                     ( (2.0_EB*PI-ABS(ANGLE-HR%ANGLE))*SIGN(1.0_EB , HR%ANGLE-ANGLE)*(OMEGA_0/PI) - HR%OMEGA)
             END IF
          ELSE  ! No target direction
             HR%UBAR = 0.0_EB
             HR%VBAR = 0.0_EB
             U_NEW = U_NEW / FAC_TIM
             V_NEW = V_NEW / FAC_TIM
             ! Slow rotation down if no direction available, i.e., target omega_0 is zero.
             OMEGA_NEW = OMEGA_NEW + 0.5_EB*(DTSP/HR_TAU_INER)*(-HR%OMEGA)
          END IF
          ! Check that velocities are not too large, i.e., less than 10 m/s for humans
          EVEL = SQRT(U_NEW**2 + V_NEW**2)
          IF ( EVEL > VMAX_TIMO ) THEN
             U_NEW = U_NEW*(VMAX_TIMO/EVEL)
             V_NEW = V_NEW*(VMAX_TIMO/EVEL)
          END IF
          ! Check that angular velocity is not too large
          IF ( ABS(OMEGA_NEW) > OMEGA_MAX ) THEN
             OMEGA_NEW = SIGN(OMEGA_MAX,OMEGA_NEW)
          END IF
          IF (TAU_FALL_DOWN >= 0.0_EB .AND. L_FALLEN_DOWN) THEN
            T_tmp = MIN(TAU_FALL_DOWN, T - HR%T_FallenDown) ! 0 - tau_fall_down
            U_NEW = (TAU_FALL_DOWN - T_tmp)*U_NEW/TAU_FALL_DOWN
            V_NEW = (TAU_FALL_DOWN - T_tmp)*V_NEW/TAU_FALL_DOWN
            OMEGA_NEW = (TAU_FALL_DOWN - T_tmp)*OMEGA_NEW/TAU_FALL_DOWN
          END IF

          HR%U = U_NEW
          HR%V = V_NEW
          HR%OMEGA = OMEGA_NEW

          ! ========================================================
          ! Decide the time step for the next human movement loop
          ! ========================================================
          ! Distances to closest walls and other agents
          D_HUMANS = MAX(D_HUMANS,0.0005_EB)        ! This agent
          HR%D_Walls = MIN(D_WALLS,10.0_EB)
          D_WALLS  = MAX(D_WALLS, 0.0005_EB)        ! This agent
          D_HUMANS_MIN = MIN(D_HUMANS_MIN,D_HUMANS) ! Among all agents
          D_WALLS_MIN  = MIN(D_WALLS_MIN, D_WALLS)  ! Among all agents

          IF ( T > T_BEGIN ) THEN
             ! Time step, do not move too close to other agents, walls, or 0.5*grid spacing.
             ! DELTA_MIN is minimum of DX and DY for the mesh.  minimum movement is 0.0001 m.
             DT_LOOP: DO
                U_TMP(2) = HR%U + 0.5_EB*DTSP_NEW*P2P_U/HR%MASS
                V_TMP(2) = HR%V + 0.5_EB*DTSP_NEW*P2P_V/HR%MASS
                OMEGA_NEW = HR%OMEGA + 0.5_EB*DTSP*P2P_TORQUE/HR%M_INER
                U_TMP(1) = U_TMP(2) + COS(HR%ANGLE)*OMEGA_NEW*HR%D_SHOULDER
                U_TMP(3) = U_TMP(2) - COS(HR%ANGLE)*OMEGA_NEW*HR%D_SHOULDER
                V_TMP(1) = V_TMP(2) + SIN(HR%ANGLE)*OMEGA_NEW*HR%D_SHOULDER
                V_TMP(3) = V_TMP(2) - SIN(HR%ANGLE)*OMEGA_NEW*HR%D_SHOULDER
                IF ( MAX(U_TMP(1)**2+V_TMP(1)**2, U_TMP(2)**2+V_TMP(2)**2, U_TMP(3)**2+V_TMP(3)**2)*DTSP_NEW**2 > &
                     (MIN(0.2_EB*D_HUMANS, 0.2_EB*D_WALLS, 0.5_EB*DELTA_MIN))**2 ) THEN
                   DTSP_NEW = DTSP_NEW*0.8_EB
                   CYCLE DT_LOOP
                END IF
                EXIT DT_LOOP
             END DO DT_LOOP

             DTSP_NEW = MAX(DTSP_NEW, EVAC_DT_MIN)
             DTSP_NEW = MIN(DTSP_NEW, EVAC_DT_MAX)
          ELSE  ! Initialization phase
             HR%U = 0.5_EB*HR%U ; HR%V = 0.5_EB*HR%V
             EVEL = SQRT(HR%U**2 + HR%V**2)
             IF ( EVEL > 0.1_EB*VMAX_TIMO ) THEN
                HR%U = HR%U*(0.1_EB*VMAX_TIMO/EVEL)
                HR%V = HR%V*(0.1_EB*VMAX_TIMO/EVEL)
             END IF
             DTSP_NEW = MAX(DTSP_NEW, EVAC_DT_MIN)
             DTSP_NEW = MIN(DTSP_NEW, EVAC_DT_MAX)
             IF (ICYC==0) THEN
                HR%U = 0.0_EB; HR%V = 0.0_EB; HR%OMEGA = 0.0_EB
             END IF
          END IF

       END DO EVAC_FORCE_LOOP
       DEALLOCATE(BLOCK_LIST)
       DEALLOCATE(BLOCK_GRID)
       ! ========================================================
       ! Force loop ends here
       ! ========================================================

       ! Save the maximum time steps allowed by the SC-VV algorithm
       TSTEPS(NM) = DTSP_NEW
       IF (L_FIRST_PASS) L_FIRST_PASS = .FALSE.

    END DO HUMAN_TIME_LOOP
    DEALLOCATE(BLOCK_GRID_N)
    ! ========================================================
    ! Human time loop ends here (TIN ==> T = TIN + DT)
    ! ========================================================

    ! ========================================================
    ! Evacuation routine ends here
    ! ========================================================
    T_USED(12)=T_USED(12)+CURRENT_TIME()-TNOW

  CONTAINS

!!$    SUBROUTINE GUARD(NM, XI, YJ, WSPA, WSPB, TARGET_X, TARGET_Y)
!!$      ! Dane potrzebne na wejxxciu:
!!$      !   mapWidth
!!$      !   mapHeight
!!$      !   tileSize
!!$      !   numberPeopl
!!$      !   startX,Y
!!$      !   targetX,Y
!!$      !   walkability
!!$      IMPLICIT NONE
!!$      INTEGER, INTENT(IN) :: NM
!!$      REAL(EB), INTENT(IN) :: XI, YJ, TARGET_X, TARGET_Y
!!$      REAL(EB), INTENT(OUT) :: WSPA, WSPB
!!$      INTEGER :: ic, iw, xt, yt
!!$      TYPE (MESH_TYPE), POINTER :: M =>NULL()
!!$
!!$      INTEGER :: mapWidth, mapHeight, tileSize, numberPeople = 2, &
!!$           onClosedList = 10, &
!!$           notStarted = 0, &
!!$           found = 1, nonexistent = 2, &
!!$           walkable = 0, unwalkable = 1
!!$      ! Create needed arrays
!!$      INTEGER, DIMENSION(:,:), ALLOCATABLE :: walkability, pathBank
!!$      INTEGER, DIMENSION(:), ALLOCATABLE :: openList, openX, openY, Fcost, &
!!$           Hcost, pathLength, pathLocation, pathStatus, xPath, yPath
!!$      INTEGER, DIMENSION(:,:), ALLOCATABLE :: whichList, parentX, parentY, &
!!$           Gcost
!!$
!!$      ! Local pathfinding
!!$      INTEGER :: onOpenList=0, parentXval=0, parentYval=0, a=0, b=0, n=0, &
!!$           u=0, v=0, temp=0, corner=0, numberOfOpenListItems=0, addedGCost=0, &
!!$           tempGcost=0, path=0, tempx, pathX, pathY, newOpenListItemID=0
!!$      INTEGER :: startX, startY, x, y
!!$      INTEGER :: targetX, targetY, pathfinderID ! in
!!$
!!$      M => MESHES(NM)
!!$
!!$      mapWidth = NINT((M%XF - M%XS)*100)
!!$      mapHeight = NINT((M%YF - M%YS)*100)
!!$      tileSize = mapWidth/M%IBAR
!!$
!!$      ALLOCATE(walkability(0:M%IBAR,0:M%JBAR)) ! Dimensions of mesh
!!$      ALLOCATE(openList(0:mapWidth*mapHeight+2))
!!$      ALLOCATE(whichList(0:mapWidth+1,0:mapHeight+1))
!!$      ALLOCATE(openX(0:mapWidth*mapHeight+2))
!!$      ALLOCATE(openY(0:mapWidth*mapHeight+2))
!!$      ALLOCATE(parentX(0:mapWidth+1,0:mapHeight+1))
!!$      ALLOCATE(parentY(0:mapWidth+1,0:mapHeight+1))
!!$      ALLOCATE(Fcost(0:mapWidth*mapHeight+2))
!!$      ALLOCATE(Gcost(0:mapWidth+1,0:mapHeight+1))
!!$      ALLOCATE(Hcost(0:mapWidth*mapHeight+2))
!!$      ALLOCATE(pathLength(0:numberPeople+1))
!!$      ALLOCATE(pathLocation(0:numberPeople+1))
!!$
!!$      ! Path reading variables
!!$      ALLOCATE(pathStatus(0:numberPeople+1))
!!$      ALLOCATE(xPath(0:numberPeople+1))
!!$      ALLOCATE(yPath(0:numberPeople+1))
!!$
!!$      ! Mam juxx wymiar mapy
!!$      ! WRITE(LU_ERR,*) M%XS
!!$      ! WRITE(LU_ERR,*) M%XF
!!$      ! WRITE(LU_ERR,*) M%YS
!!$      ! WRITE(LU_ERR,*) M%YF
!!$
!!$      ! Oraz na ile jest podzielona elementxxw
!!$      ! WRITE(LU_ERR,*) M%IBAR
!!$      ! WRITE(LU_ERR,*) M%JBAR
!!$
!!$      ! Teraz podzielxx jxx na komxxrki i sprawdzxx czy sxx przeszkody
!!$      ! Ewidentnie robimy tylko raz
!!$      DO xt = 1, M%IBAR
!!$         DO yt = 1, M%JBAR
!!$            ic  = M%CELL_INDEX(xt,yt,1)
!!$            iw  = M%OBST_INDEX_C(ic)
!!$            walkability(xt-1,yt-1) = iw
!!$         END DO
!!$      END DO
!!$
!!$      ! 1. Convert location data (in pixels) to coordinates in the walkability array.
!!$      !startX = startingX/tileSize
!!$      !startY = startingY/tileSize
!!$      !targetX = targetX/tileSize
!!$      !targetY = targetY/tileSize
!!$      startX = FLOOR(XI)
!!$      startY = FLOOR(YJ)
!!$
!!$      targetX = FLOOR(M%CELLSI(FLOOR((TARGET_X-M%XS)*M%RDXINT))+1.0_EB)
!!$      targetY = FLOOR(M%CELLSJ(FLOOR((TARGET_Y-M%YS)*M%RDYINT))+1.0_EB)
!!$      !write(lu_err,*)'**** targetX, targetY ',targetX, targetY,target_X, target_Y
!!$      !targetX = 19
!!$      !targetY = 19
!!$
!!$      ! 2.Quick Path Checks: Under the some circumstances no path needs to
!!$      ! be generated ...
!!$      ! If starting location and target are in the same location...
!!$      IF((startX==targetX .AND. startY==targetY .AND. pathLocation(pathfinderID)>0) .OR.(startX==targetX .AND. &
!!$           startY==targetY .AND. pathLocation(pathfinderID)==0)) THEN
!!$         !path = found
!!$      ELSE
!!$         !IF(startX==targetX .AND. startY==targetY .AND. pathLocation(pathfinderID)==0) THEN
!!$         !path = nonexistent
!!$         !END IF
!!$
!!$         ! If target square is unwalkable, return that it's a nonexistent path.
!!$         !if (walkability[targetX][targetY] == unwalkable)
!!$         !goto noPath;
!!$
!!$         ! 3.Reset some variables that need to be cleared
!!$         IF(onClosedList>1000000) THEN ! reset whichList occasionally
!!$            DO x=0, mapWidth-1
!!$               DO y=0, mapHeight-1
!!$                  whichList(x,y) = 0
!!$               END DO
!!$            END DO
!!$            onClosedList = 10
!!$         END IF
!!$         ! changing the values of onOpenList and onClosed list is faster than redimming whichList() array
!!$         onClosedList = onClosedList+2
!!$         onOpenList = onClosedList-1
!!$         pathLength(pathfinderID) = notStarted ! i.e, = 0
!!$         pathLocation(pathfinderID) = notStarted ! i.e, = 0
!!$         Gcost(startX,startY) = 0 ! reset starting square's G value to 0
!!$
!!$         ! 4.Add the starting location to the open list of squares to be checked.
!!$         numberOfOpenListItems = 1
!!$         ! assign it as the top (and currently only) item in the open list, which is maintained
!!$         ! as a binary heap (explained below)
!!$         openList(1) = 1
!!$         openX(1) = startX
!!$         openY(1) = startY
!!$
!!$         ! 5.Do the following until a path is found or deemed nonexistent.
!!$         DO
!!$
!!$            ! 6.If the open list is not empty, take the first cell off of the list.
!!$            ! This is the lowest F cost cell on the open list.
!!$            IF(numberOfOpenListItems /= 0) THEN
!!$
!!$               ! 7. Pop the first item off the open list.
!!$               parentXval = openX(openList(1))
!!$               parentYval = openY(openList(1)) ! record cell coordinates of the item
!!$               whichList(parentXval,parentYval) = onClosedList ! add the item to the closed list
!!$
!!$               ! Open List = Binary Heap: Delete this item from the open list, which
!!$               ! is maintained as a binary heap. For more information on binary heaps, see:
!!$               ! http://www.policyalmanac.org/games/binaryHeaps.htm
!!$               numberOfOpenListItems = numberOfOpenListItems - 1 ! reduce number of open list items by 1
!!$
!!$               ! Delete the top item in binary heap and reorder the heap, with the lowest F cost item rising to the top.
!!$               openList(1) = openList(numberOfOpenListItems+1) ! move the last item in the heap up to slot #1
!!$               v = 1
!!$
!!$               ! Repeat the following until the new item in slot #1 sinks to its proper spot in the heap.
!!$               DO
!!$                  u = v
!!$                  IF(2*u+1<=numberOfOpenListItems) THEN ! if both children exist
!!$                     ! Check if the F cost of the parent is greater than each child.
!!$                     ! Select the lowest of the two children.
!!$                     IF(Fcost(openList(u)) >= Fcost(openList(2*u))) THEN
!!$                        v = 2*u;
!!$                     END IF
!!$                     IF(Fcost(openList(v)) >= Fcost(openList(2*u+1))) THEN
!!$                        v = 2*u+1
!!$                     END IF
!!$                  ELSE
!!$                     IF(2*u<=numberOfOpenListItems) THEN ! if only child #1 exists
!!$                        ! Check if the F cost of the parent is greater than child #1
!!$                        IF(Fcost(openList(u)) >= Fcost(openList(2*u))) THEN
!!$                           v = 2*u;
!!$                        END IF
!!$                     END IF
!!$                  END IF
!!$
!!$                  IF(u /= v) THEN ! if parent's F is > one of its children, swap them
!!$                     temp = openList(u)
!!$                     openList(u) = openList(v)
!!$                     openList(v) = temp
!!$                  ELSE
!!$                     EXIT ! otherwise, exit loop
!!$                  END IF
!!$               END DO
!!$               !while (!KeyDown(27))//reorder the binary heap
!!$
!!$               ! 7.Check the adjacent squares. (Its "children" -- these path children
!!$               ! are similar, conceptually, to the binary heap children mentioned
!!$               ! above, but don't confuse them. They are different. Path children
!!$               ! are portrayed in Demo 1 with grey pointers pointing toward
!!$               ! their parents.) Add these adjacent child squares to the open list
!!$               ! for later consideration if appropriate (see various if statements
!!$               ! below).
!!$               DO b=parentYval-1, parentYval+1
!!$                  DO a=parentXval-1, parentXval+1
!!$                     ! If not off the map (do this first to avoid array out-of-bounds errors)
!!$                     IF(a /= -1 .AND. b /= -1 .AND. a /= mapWidth/tileSize .AND. b /= mapHeight/tileSize) THEN
!!$                        ! If not already on the closed list (items on the closed list have
!!$                        ! already been considered and can now be ignored).
!!$                        IF(whichList(a,b) /= onClosedList) THEN
!!$                           ! If not a wall/obstacle square.
!!$                           IF(walkability(a,b) < unwalkable) THEN
!!$                              ! Don't cut across corners
!!$                              corner = walkable
!!$                              IF(a==parentXval-1) THEN
!!$                                 IF(b == parentYval-1) THEN
!!$                                    IF(walkability(parentXval-1,parentYval) >= unwalkable .OR. &
!!$                                         walkability(parentXval,parentYval-1) >= unwalkable) THEN
!!$                                       corner = unwalkable
!!$                                    END IF
!!$                                 ELSE IF(b == parentYval+1) THEN
!!$                                    IF(walkability(parentXval,parentYval+1) >= unwalkable .OR. &
!!$                                         walkability(parentXval-1,parentYval) >= unwalkable) THEN
!!$                                       corner = unwalkable
!!$                                    END IF
!!$                                 END IF
!!$                              ELSE IF(a == parentXval+1) THEN
!!$                                 IF(b == parentYval-1) THEN
!!$                                    IF(walkability(parentXval,parentYval-1) >= unwalkable .OR. &
!!$                                         walkability(parentXval+1,parentYval) >= unwalkable) THEN
!!$                                       corner = unwalkable
!!$                                    END IF
!!$                                 ELSE IF(b == parentYval+1) THEN
!!$                                    IF(walkability(parentXval+1,parentYval) >= unwalkable .OR. &
!!$                                         walkability(parentXval,parentYval+1) >= unwalkable) THEN
!!$                                       corner = unwalkable
!!$                                    END IF
!!$                                 END IF
!!$                              END IF
!!$                              IF(corner == walkable) THEN
!!$                                 ! If not already on the open list, add it to the open list.
!!$                                 IF(whichList(a,b) /= onOpenList) THEN
!!$                                    ! Create a new open list item in the binary heap.
!!$                                    newOpenListItemID = newOpenListItemID + 1 ! each new item has a unique ID #
!!$                                    n = numberOfOpenListItems+1
!!$                                    ! place the new open list item (actually, its ID#) at the bottom of the heap
!!$                                    openList(n) = newOpenListItemID
!!$                                    openX(newOpenListItemID) = a
!!$                                    openY(newOpenListItemID) = b ! record the x and y coordinates of the new item
!!$
!!$                                    ! Figure out its G cost
!!$                                    IF(abs(a-parentXval) == 1 .AND. abs(b-parentYval) == 1) THEN
!!$                                       addedGCost = 14 ! cost of going to diagonal squares
!!$                                    ELSE
!!$                                       addedGCost = 10 !cost of going to non-diagonal squares
!!$                                    END IF
!!$                                    Gcost(a,b) = Gcost(parentXval,parentYval) + addedGCost
!!$                                    ! Figure out its H and F costs and parent
!!$                                    Hcost(openList(n)) = 10*(abs(a - targetX) + abs(b - targetY))
!!$                                    Fcost(openList(n)) = Gcost(a,b) + Hcost(openList(n))
!!$                                    parentX(a,b) = parentXval
!!$                                    parentY(a,b) = parentYval
!!$                                    ! Move the new open list item to the proper place in the binary heap.
!!$                                    ! Starting at the bottom, successively compare to parent items,
!!$                                    ! swapping as needed until the item finds its place in the heap
!!$                                    ! or bubbles all the way to the top (if it has the lowest F cost).
!!$                                    DO WHILE(n /= 1) ! While item hasn't bubbled to the top (n=1)
!!$                                       ! Check if child's F cost is < parent's F cost. If so, swap them.
!!$                                       IF(Fcost(openList(n)) <= Fcost(openList(n/2))) THEN
!!$                                          temp = openList(n/2)
!!$                                          openList(n/2) = openList(n)
!!$                                          openList(n) = temp
!!$                                          n = n/2
!!$                                       ELSE
!!$                                          EXIT
!!$                                       END IF
!!$                                    END DO
!!$                                    numberOfOpenListItems = numberOfOpenListItems+1 ! add one to the number of items in the heap
!!$                                    ! Change whichList to show that the new item is on the open list.
!!$                                    whichList(a,b) = onOpenList
!!$
!!$                                    ! 8.If adjacent cell is already on the open list, check to see if this
!!$                                    ! path to that cell from the starting location is a better one.
!!$                                    ! If so, change the parent of the cell and its G and F costs.
!!$                                 ELSE ! If whichList(a,b) = onOpenList
!!$                                    ! Figure out the G cost of this possible new path
!!$                                    IF(abs(a-parentXval) == 1 .AND. abs(b-parentYval) == 1) THEN
!!$                                       addedGCost = 14 ! cost of going to diagonal tiles
!!$                                    ELSE
!!$                                       addedGCost = 10 ! cost of going to non-diagonal tiles
!!$                                    END IF
!!$                                    tempGcost = Gcost(parentXval,parentYval) + addedGCost
!!$                                    ! If this path is shorter (G cost is lower) then change
!!$                                    ! the parent cell, G cost and F cost.
!!$                                    IF(tempGcost < Gcost(a,b)) THEN ! if G cost is less,
!!$                                       parentX(a,b) = parentXval ! change the square's parent
!!$                                       parentY(a,b) = parentYval
!!$                                       Gcost(a,b) = tempGcost ! change the G cost
!!$                                       ! Because changing the G cost also changes the F cost, if
!!$                                       ! the item is on the open list we need to change the item's
!!$                                       ! recorded F cost and its position on the open list to make
!!$                                       ! sure that we maintain a properly ordered open list.
!!$                                       DO x = 1, numberOfOpenListItems !//look for the item in the heap
!!$                                          IF(openX(openList(x)) == a .AND. openY(openList(x)) == b) THEN ! item found
!!$                                             Fcost(openList(x)) = Gcost(a,b) + Hcost(openList(x)) ! change the F cost
!!$                                             ! See if changing the F score bubbles the item up from it's current location
!!$                                             ! in the heap
!!$                                             n = x;
!!$                                             DO WHILE (n /= 1) ! While item hasn't bubbled to the top (n=1)
!!$                                                ! Check if child is < parent. If so, swap them.
!!$                                                IF(Fcost(openList(n)) < Fcost(openList(n/2))) THEN
!!$                                                   temp = openList(n/2)
!!$                                                   openList(n/2) = openList(n)
!!$                                                   openList(n) = temp
!!$                                                   n = n/2
!!$                                                ELSE
!!$                                                   EXIT
!!$                                                END IF
!!$                                             END DO
!!$                                             EXIT ! exit for x = loop
!!$                                          END IF ! If openX(openList(x)) = a
!!$                                       END DO ! For x = 1 To numberOfOpenListItems
!!$                                    END IF ! If tempGcost < Gcost(a,b)
!!$                                 END IF ! else If whichList(a,b) = onOpenList
!!$                              END IF ! If not cutting a corner
!!$                           END IF ! If not a wall/obstacle square.
!!$                        END IF ! If not already on the closed list
!!$                     END IF ! If not off the map
!!$                  END DO ! for (a = parentXval-1; a <= parentXval+1; a++){
!!$               END DO ! for (b = parentYval-1; b <= parentYval+1; b++){
!!$            ELSE ! if (numberOfOpenListItems != 0)
!!$               ! 9.If open list is empty then there is no path.
!!$               path = nonexistent
!!$               EXIT
!!$            END IF
!!$            ! If target is added to open list then path has been found.
!!$            IF(whichList(targetX,targetY) == onOpenList) THEN
!!$               path = found
!!$               EXIT
!!$            END IF
!!$         END DO
!!$
!!$         ! 10.Save the path if it exists.
!!$         IF(path == found) THEN
!!$            ! a.Working backwards from the target to the starting location by checking
!!$            ! each cell's parent, figure out the length of the path.
!!$            pathX = targetX
!!$            pathY = targetY
!!$            DO WHILE(pathX /= startX .OR. pathY /= startY)
!!$               ! Look up the parent of the current cell.
!!$               tempx = parentX(pathX,pathY)
!!$               pathY = parentY(pathX,pathY)
!!$               pathX = tempx
!!$               ! Figure out the path length
!!$               pathLength(pathfinderID) = pathLength(pathfinderID) + 1;
!!$            END DO
!!$            !b.Resize the data bank to the right size in bytes
!!$            ALLOCATE(pathBank(0:pathLength(pathfinderID),0:1))
!!$
!!$            !c. Now copy the path information over to the databank. Since we are
!!$            ! working backwards from the target to the start location, we copy
!!$            ! the information to the data bank in reverse order. The result is
!!$            ! a properly ordered set of path data, from the first step to the
!!$            ! last.
!!$            pathX = targetX
!!$            pathY = targetY
!!$            pathBank(0,0) = startX
!!$            pathBank(0,1) = startY
!!$            !cellPosition = pathLength(pathfinderID)*2 ! start at the end
!!$            DO WHILE(pathX /= startX .OR. pathY /= startY)
!!$               !cellPosition = cellPosition - 2 ! work backwards 2 integers
!!$               !pathBank(pathfinderID,cellPosition) = pathX
!!$               !pathBank(pathfinderID,cellPosition+1) = pathY
!!$               ! d.Look up the parent of the current cell.
!!$               tempx = parentX(pathX,pathY)
!!$               pathY = parentY(pathX,pathY)
!!$               pathX = tempx
!!$               pathBank(pathLength(pathfinderID),0) = pathX
!!$               pathBank(pathLength(pathfinderID),1) = pathY
!!$               pathLength(pathfinderID) = pathLength(pathfinderID) - 1
!!$               ! WRITE(LU_ERR,*) 'pathX=',pathX
!!$               ! WRITE(LU_ERR,*) 'pathY=',pathY
!!$               ! e.If we have reached the starting square, exit the loop.
!!$            END DO
!!$            ! WRITE(LU_ERR,*) 'koniec----------'
!!$
!!$            ! 11.Read the first path step into xPath/yPath arrays
!!$            !ReadPath(pathfinderID,startingX,startingY,1);
!!$         END IF
!!$         !return path;
!!$
!!$
!!$
!!$         ! 13.If there is no path to the selected target, set the pathfinder's
!!$         ! xPath and yPath equal to its current location and return that the
!!$         ! path is nonexistent.
!!$         ! noPath:
!!$         ! xPath[pathfinderID] = startingX;
!!$         ! yPath[pathfinderID] = startingY;
!!$         ! return nonexistent;
!!$
!!$         !WRITE(LU_ERR,*) 'Xend=',pathBank(2,0)
!!$         !WRITE(LU_ERR,*) 'Xsta=',startX
!!$         !WRITE(LU_ERR,*) 'Yend=',pathBank(2,1)
!!$         !WRITE(LU_ERR,*) 'Ysta=',startY
!!$         !WRITE(LU_ERR,*) '====================='
!!$
!!$         WSPA = pathBank(2,0) - startX
!!$         WSPB = pathBank(2,1) - startY
!!$
!!$         IF(path == found) THEN
!!$            DEALLOCATE(pathBank)
!!$         END IF
!!$
!!$      END IF ! if startX==targetX
!!$
!!$      DEALLOCATE(walkability)
!!$      DEALLOCATE(openList)
!!$      DEALLOCATE(whichList)
!!$      DEALLOCATE(openX)
!!$      DEALLOCATE(openY)
!!$      DEALLOCATE(parentX)
!!$      DEALLOCATE(parentY)
!!$      DEALLOCATE(Fcost)
!!$      DEALLOCATE(Gcost)
!!$      DEALLOCATE(Hcost)
!!$      DEALLOCATE(pathLength)
!!$      DEALLOCATE(pathLocation)
!!$
!!$      ! Path reading variables
!!$      DEALLOCATE(pathStatus)
!!$      DEALLOCATE(xPath)
!!$      DEALLOCATE(yPath)
!!$
!!$    END SUBROUTINE GUARD

    SUBROUTINE FIND_PREFERRED_DIRECTION(I, NOUT, T, T_BEGIN, L_DEAD, NM_STRS_MESH, &
         II, JJ, XI, YJ, UBAR, VBAR, HR_TAU, TPRE, NM, I_STRS_DOOR, HR_SPEED)
      IMPLICIT NONE
      !
      ! Calculate the prefered walking direction
      !
      ! Inputs:
      !   I: The index of the agent (on this mesh)
      !   T: Time
      !   T_BEGIN: The starting time of the simulation
      !   L_DEAD: Is the agent dead or not
      !   II,JJ: The grid cell indices of the agent
      !   XI,YJ: The grid cell coordinates of the agent for the velocity
      !   NM: The main evacuation mesh index
      !   NOUT: The index of the stairs mesh (if the agent is on stairs)
      !   NM_STRS_MESH: True, if the mesh is a stair mesh
      ! Outputs:
      !   UBAR,VBAR: The prefered walking direction
      !   HR_TAU: Translational motion tau
      !   TPRE: Pre-movement time
      !   I_STRS_DOOR: In STRS mesh the target door/exit index
      !
      !   HR%I_Door_Mode: 1: Use flow field, 2: use bee line to the door
      !
      ! Passed variables
      INTEGER, INTENT(IN) :: II, JJ, I, NM, NOUT
      REAL(EB), INTENT(IN) :: XI, YJ, T, T_BEGIN
      LOGICAL, INTENT(IN) :: L_DEAD, NM_STRS_MESH
      INTEGER, INTENT(OUT) :: I_STRS_DOOR
      REAL(EB), INTENT(INOUT) :: HR_TAU
      REAL(EB), INTENT(OUT) :: UBAR, VBAR, TPRE, HR_SPEED
      !
      ! Local variables
      INTEGER :: N,NM_NOW, STRS_INDX, DOOR_IOR, KKZ, J
      REAL(EB) :: X_TARGET, Y_TARGET, DOOR_WIDTH, DOOR_DIST, EVEL, X1, X2, Y1, Y2, D_TMP, X_XYZ, Y_XYZ, D_X, D_Y
      REAL(EB) :: alpha, Dv0, Dv0_x, Dv0_y, D_DrV0, Dr_length, Dr_x, Dr_y
      LOGICAL :: STRAIGHT_LINE_TO_TARGET, Is_Known_Door_tmp, Is_XB_Visible, Is_XYZ_Visible, Is_InFront
      LOGICAL :: V0_IS_SET_ZERO
      TYPE (MESH_TYPE), POINTER :: MFF=>NULL()
      TYPE (HUMAN_TYPE), POINTER :: HR=>NULL()
      TYPE (EVAC_STRS_TYPE), POINTER :: STRP=>NULL()
      INTEGER :: I_TARGET_TMP, INODE, I_TARGET_OLD, III, N_QUEUE, JJ_NOW, N_END
      REAL(EB) :: DIST_TO_DOOR, DIST_TO_DOOR_TMP, X_NODE, Y_NODE, WIDTH, T_TMP1
      REAL(EB) :: ANGLE_HR, ANGLE_XYZ, ANGLE_XB, D_NEW, RN

      HR=>HUMAN(I)
      KKZ = 1
      ! HR%I_FFIELD is the flow field mesh index for this agent
      NM_NOW = MAX(1,HR%I_FFIELD)
      IF (.NOT.EVACUATION_ONLY(NM_NOW) .OR. HR%I_FFIELD < 1) THEN
         WRITE(LU_ERR,*)'*** ',I,HR%IPC,HR%IEL,HR%I_TARGET,HR%I_FFIELD,HR%ILABEL
         WRITE(LU_ERR,'(A,I6,A)') 'ERROR: EVACUATE_HUMANS, mesh ', HR%I_FFIELD, ' is not an evacuation flow field.'
         WRITE(LU_EVACOUT,'(A,I6,A)') 'ERROR: EVACUATE_HUMANS, mesh ', HR%I_FFIELD, ' is not an evacuation flow field.'
         STOP_STATUS = EVACUATION_STOP ; RETURN
      END IF
      IF (T <= T_BEGIN) THEN
         ! Initialization phase, i.e., flow field calculation
         UBAR = 0.0_EB ; VBAR = 0.0_EB ; HR_TAU = HR%TAU
         RETURN
      END IF
      !
      ! Determine if the mesh is a stairs-mesh
      IF (NM_STRS_MESH) THEN
         N = NOUT
         NM_NOW = NM
         STRS_INDX = N
         STRP=>EVAC_STRS(N)
      END IF

      MFF=>MESHES(NM_NOW)  ! Pointer to the flow field mesh

      ! Check if going straight line to target
      STRAIGHT_LINE_TO_TARGET = .FALSE.
      I_STRS_DOOR = 0
      HR%SKIP_WALL_FORCE_IOR = 0
      V0_IS_SET_ZERO = .FALSE.
      If_StrsMesh: IF (NM_STRS_MESH) THEN
         IF (HR%I_TARGET == 0) THEN
            CALL FIND_TARGET_NODE_IN_STRS(STRP,HR)
         END IF
         N = ABS(HR%I_TARGET)
         IF (N>N_DOORS) THEN
            N = N - N_DOORS
            IF (EVAC_EXITS(N)%STR_SUB_INDX == HR%STR_SUB_INDX) THEN
               STRAIGHT_LINE_TO_TARGET = .TRUE.
            END IF
         ELSE IF (N>0) THEN
            IF (EVAC_DOORS(N)%STR_SUB_INDX == HR%STR_SUB_INDX) THEN
               STRAIGHT_LINE_TO_TARGET = .TRUE.
            END IF
         END IF
         ! Test if this agent should update the exit door (exp. distribution)
         RN = 2.0_EB
         IF (HR%I_TARGET > 0 .AND. STRAIGHT_LINE_TO_TARGET) THEN
            ! The agent is at the final landing and has already chosen the closest door once
            CALL RANDOM_NUMBER(RN_REAL)
            RN = REAL(RN_REAL,EB)
         END IF
         If_StraightLineNewDoor: IF (STRAIGHT_LINE_TO_TARGET .AND. RN > EXP(-DTSP/TAU_CHANGE_DOOR)) THEN
            ! Find the target among the nodes leading out of the strs (doors and exits)
            I_TARGET_OLD = HR%I_TARGET
            DIST_TO_DOOR = HUGE(DIST_TO_DOOR)

            NODES_OUT_LOOP: DO J = 1, STRP%N_NODES_OUT
               INODE = STRP%NODES_OUT(J)

               SELECT CASE(EVAC_NODE_LIST(INODE)%NODE_TYPE)
               CASE('Door')
                  IF (EVAC_DOORS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%STR_SUB_INDX /= HR%STR_SUB_INDX) CYCLE NODES_OUT_LOOP

                  ! Group index=0: the agent is from an entry line (no evac line)
                  ! J  =  MAX(0,HR%GROUP_ID)   ! Group index
                  ! J1 = -MIN(0,HR%GROUP_ID)   ! Lonely agent index
                  Is_Known_Door_tmp = .FALSE.
                  IF (EVAC_DOORS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%KNOWN_DOOR .AND. .NOT.HR%I_DoorAlgo==3) &
                       Is_Known_Door_tmp = .TRUE.
                  IF ( HR%GROUP_ID /= 0 ) THEN
                     IF ( HR%GROUP_ID < 0 ) THEN
                        ! A lonely soul
                        DO iii = 1, Human_Known_Doors(-HR%GROUP_ID)%N_nodes
                           IF (EVAC_DOORS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%INODE == &
                                ABS(Human_Known_Doors(-HR%GROUP_ID)%I_nodes(iii))) Is_Known_Door_tmp = .TRUE.
                        END DO
                     ELSE
                        ! A member of a group
                        DO iii = 1, Group_Known_Doors(HR%GROUP_ID)%N_nodes
                           IF (EVAC_DOORS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%INODE == &
                                ABS(Group_Known_Doors(HR%GROUP_ID)%I_nodes(iii))) Is_Known_Door_tmp = .TRUE.
                        END DO
                     END IF
                  ELSE
                     ! The agent is from an entry line. P_door= 0.0 or 1.0, known doors
                     ! i_door_nodes: <0 exit number, >0 door number
                     DO iii = 1, EVAC_ENTRYS(ABS(HR%IEL))%N_VENT_FFIELDS
                        IF (EVAC_DOORS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%INODE == EVAC_ENTRYS(ABS(HR%IEL))%I_DOOR_NODES(iii) &
                             .AND. EVAC_ENTRYS(ABS(HR%IEL))%P_VENT_FFIELDS(III)>0.5_EB ) Is_Known_Door_tmp = .TRUE.
                     END DO
                  END IF

                  ! IF (.NOT. EVAC_DOORS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%EXIT_SIGN) CYCLE NODES_IN_LOOP
                  ! Put here a check: If a known door for this agent, then exit_sign = T or F are good.
                  ! If not known door, then exit_sign=T only counted as a possible door.
                  X_NODE = 0.5_EB*(EVAC_DOORS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%X1 + &
                       EVAC_DOORS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%X2)
                  Y_NODE = 0.5_EB*(EVAC_DOORS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%Y1 + &
                       EVAC_DOORS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%Y2)
                  I_TARGET_TMP = EVAC_NODE_LIST(INODE)%NODE_INDEX
                  DIST_TO_DOOR_TMP = 50.0_EB*SQRT( (X_NODE-HR%X)**2 + (Y_NODE-HR%Y)**2 ) / &
                       EVAC_DOORS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%R_NTARGET + 1.0_EB
                  III = MIN(50,MAX(1,INT(DIST_TO_DOOR_TMP)-1))
                  WIDTH   = EVAC_DOORS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%WIDTH
                  N_QUEUE = EVAC_DOORS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%NTARGET(III)
               CASE('Exit')
                  ! Next check is actually not needed, because count_only exits are not counted as STR nodes.
                  IF (EVAC_EXITS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%COUNT_ONLY) CYCLE NODES_OUT_LOOP
                  IF (EVAC_EXITS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%STR_SUB_INDX /= HR%STR_SUB_INDX) CYCLE NODES_OUT_LOOP

                  ! Group index=0: the agent is from an entry line (no evac line)
                  ! J  =  MAX(0,HR%GROUP_ID)   ! Group index
                  ! J1 = -MIN(0,HR%GROUP_ID)   ! Lonely agent index
                  Is_Known_Door_tmp = .FALSE.
                  IF (EVAC_EXITS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%KNOWN_DOOR .AND. .NOT.HR%I_DoorAlgo==3) &
                       Is_Known_Door_tmp = .TRUE.
                  IF ( HR%GROUP_ID /= 0 ) THEN
                     IF ( HR%GROUP_ID < 0 ) THEN
                        ! A lonely soul
                        DO iii = 1, Human_Known_Doors(-HR%GROUP_ID)%N_nodes
                           IF (EVAC_EXITS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%INODE == &
                                ABS(Human_Known_Doors(-HR%GROUP_ID)%I_nodes(iii))) Is_Known_Door_tmp = .TRUE.
                        END DO
                     ELSE
                        ! A member of a group
                        DO iii = 1, Group_Known_Doors(HR%GROUP_ID)%N_nodes
                           IF (EVAC_EXITS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%INODE == &
                                ABS(Group_Known_Doors(HR%GROUP_ID)%I_nodes(iii))) Is_Known_Door_tmp = .TRUE.
                        END DO
                     END IF
                  ELSE
                     ! The agent is from an entry line. P_door= 0.0 or 1.0, known doors
                     ! i_door_nodes: <0 exit number, >0 door number
                     DO iii = 1, EVAC_ENTRYS(ABS(HR%IEL))%N_VENT_FFIELDS
                        IF (EVAC_EXITS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%INODE == EVAC_ENTRYS(ABS(HR%IEL))%I_DOOR_NODES(iii) &
                             .AND. EVAC_ENTRYS(ABS(HR%IEL))%P_VENT_FFIELDS(III)>0.5_EB ) Is_Known_Door_tmp = .TRUE.
                     END DO
                  END IF

                  X_NODE = 0.5_EB*(EVAC_EXITS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%X1 + &
                       EVAC_EXITS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%X2)
                  Y_NODE = 0.5_EB*(EVAC_EXITS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%Y1 + &
                       EVAC_EXITS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%Y2)
                  I_TARGET_TMP = EVAC_NODE_LIST(INODE)%NODE_INDEX + N_DOORS
                  DIST_TO_DOOR_TMP = 50.0_EB*SQRT( (X_NODE-HR%X)**2 + (Y_NODE-HR%Y)**2 ) / &
                       EVAC_EXITS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%R_NTARGET + 1.0_EB
                  III = MIN(50,MAX(1,INT(DIST_TO_DOOR_TMP)-1))
                  WIDTH   = EVAC_EXITS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%WIDTH
                  N_QUEUE = EVAC_EXITS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%NTARGET(III)
               CASE DEFAULT
                  WRITE(LU_ERR,*) 'ERROR (DEBUG): Unknown node type in Find_prefered_direction nodes_out_loop'
               END SELECT
               ! It is assumed that exit sign is on either side of the door
               DIST_TO_DOOR_TMP = 0.0_EB
               IF (.NOT.Is_Known_Door_tmp) DIST_TO_DOOR_TMP = DIST_TO_DOOR_TMP + 1000.0_EB
               IF (Trim(EVAC_NODE_LIST(INODE)%NODE_TYPE)=='Door') THEN
                  IF (.NOT.EVAC_DOORS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%EXIT_SIGN .AND. .NOT.Is_Known_Door_tmp) THEN
                     DIST_TO_DOOR_TMP = DIST_TO_DOOR_TMP + 1000.0_EB ! Penalty for floor ==> door ==> stairs doors
                  END IF
               END IF
               IF (ABS(FAC_DOOR_QUEUE) >= TWO_EPSILON_EB) THEN
                  DIST_TO_DOOR_TMP = DIST_TO_DOOR_TMP + SQRT( (X_NODE-HR%X)**2 + (Y_NODE-HR%Y)**2 )
                  T_TMP1 = MIN(1.5_EB*PI*DIST_TO_DOOR_TMP**2/(ABS(FAC_DOOR_QUEUE)*WIDTH), &
                       REAL(N_QUEUE,EB)/(ABS(FAC_DOOR_QUEUE)*WIDTH))
                  DIST_TO_DOOR_TMP = DIST_TO_DOOR_TMP + (DIST_TO_DOOR_TMP/HR%SPEED) +  T_TMP1
                  IF (I_TARGET_TMP == I_TARGET_OLD) DIST_TO_DOOR_TMP = DIST_TO_DOOR_TMP*FAC_DOOR_WAIT
                  IF (DIST_TO_DOOR_TMP < DIST_TO_DOOR) THEN
                     DIST_TO_DOOR = MAX(0.0_EB,DIST_TO_DOOR_TMP)
                     HR%I_TARGET = ABS(I_TARGET_TMP)  ! i_target > 0: A good target door is found
                  END IF
               ELSE
                  DIST_TO_DOOR_TMP = DIST_TO_DOOR_TMP + SQRT( (X_NODE-HR%X)**2 + (Y_NODE-HR%Y)**2 )
                  IF (I_TARGET_TMP == I_TARGET_OLD) DIST_TO_DOOR_TMP = DIST_TO_DOOR_TMP*FAC_DOOR_WAIT
                  IF (DIST_TO_DOOR_TMP < DIST_TO_DOOR) THEN
                     DIST_TO_DOOR = MAX(0.0_EB,DIST_TO_DOOR_TMP)
                     HR%I_TARGET = ABS(I_TARGET_TMP)  ! i_target > 0: A good target door is found
                  END IF
               END IF
            END DO NODES_OUT_LOOP
            DIST_TO_DOOR = HUGE(DIST_TO_DOOR)
            ! i_target > 0: A good target door is already found in NODES_OUT
            ! Prefer NODES_OUT doors over NODES_IN doors
            If_Nodes_In_Loop: IF (HR%I_TARGET <= 0) THEN
               NODES_IN_LOOP: DO J = 1, STRP%N_NODES_IN
                  INODE = STRP%NODES_IN(J)
                  SELECT CASE(EVAC_NODE_LIST(INODE)%NODE_TYPE)
                  CASE('Door')
                     IF (EVAC_DOORS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%STR_SUB_INDX /= HR%STR_SUB_INDX) CYCLE NODES_IN_LOOP
                     X_NODE = 0.5_EB*(EVAC_DOORS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%X1 + &
                          EVAC_DOORS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%X2)
                     Y_NODE = 0.5_EB*(EVAC_DOORS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%Y1 + &
                          EVAC_DOORS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%Y2)
                     I_TARGET_TMP = EVAC_NODE_LIST(INODE)%NODE_INDEX
                     DIST_TO_DOOR_TMP = 50.0_EB*SQRT( (X_NODE-HR%X)**2 + (Y_NODE-HR%Y)**2 ) / &
                          EVAC_DOORS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%R_NTARGET + 1.0_EB
                     III = MIN(50,MAX(1,INT(DIST_TO_DOOR_TMP)-1))
                     WIDTH   = EVAC_DOORS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%WIDTH
                     N_QUEUE = EVAC_DOORS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%NTARGET(III)
                     !CASE('Exit')
                  CASE DEFAULT
                     WRITE(LU_ERR,*) 'ERROR (DEBUG): Unknown node type in Find_prefered_direction nodes_in_loop'
                  END SELECT
                  ! It is assumed that exit sign is on either side of the door
                  IF (EVAC_DOORS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%EXIT_SIGN) THEN
                     DIST_TO_DOOR_TMP = 1000.0_EB ! Penalty for floor ==> door ==> stairs doors
                  ELSE
                     DIST_TO_DOOR_TMP = 0.0_EB
                  END IF
                  IF (ABS(FAC_DOOR_QUEUE) >= TWO_EPSILON_EB) THEN
                     DIST_TO_DOOR_TMP = DIST_TO_DOOR_TMP + SQRT( (X_NODE-HR%X)**2 + (Y_NODE-HR%Y)**2 )
                     T_TMP1 = MIN(1.5_EB*PI*DIST_TO_DOOR_TMP**2/(ABS(FAC_DOOR_QUEUE)*WIDTH), &
                          REAL(N_QUEUE,EB)/(ABS(FAC_DOOR_QUEUE)*WIDTH))
                     DIST_TO_DOOR_TMP = (DIST_TO_DOOR_TMP/HR%SPEED) +  T_TMP1
                     IF (I_TARGET_TMP == I_TARGET_OLD) DIST_TO_DOOR_TMP = DIST_TO_DOOR_TMP*FAC_DOOR_WAIT
                     IF (DIST_TO_DOOR_TMP < DIST_TO_DOOR) THEN
                        DIST_TO_DOOR = MAX(0.0_EB,DIST_TO_DOOR_TMP)
                        HR%I_TARGET = ABS(I_TARGET_TMP)  ! i_target > 0: A good target door is found
                     END IF
                  ELSE
                     DIST_TO_DOOR_TMP = DIST_TO_DOOR_TMP + SQRT( (X_NODE-HR%X)**2 + (Y_NODE-HR%Y)**2 )
                     IF (I_TARGET_TMP == I_TARGET_OLD) DIST_TO_DOOR_TMP = DIST_TO_DOOR_TMP*FAC_DOOR_WAIT
                     !IF (.NOT.Is_Known_Door_tmp) DIST_TO_DOOR_TMP = DIST_TO_DOOR_TMP + 100000.0_EB
                     IF (DIST_TO_DOOR_TMP < DIST_TO_DOOR) THEN
                        DIST_TO_DOOR = MAX(0.0_EB,DIST_TO_DOOR_TMP)
                        HR%I_TARGET = ABS(I_TARGET_TMP)  ! i_target > 0: A good target door is found
                     END IF
                  END IF
               END DO NODES_IN_LOOP
            END IF If_Nodes_In_Loop
         END IF If_StraightLineNewDoor

         N = ABS(HR%I_TARGET)
         I_STRS_DOOR = N ! Which door is the target is needed in door post forces
         IF (N>N_DOORS) THEN
            N = N - N_DOORS
            IF (EVAC_EXITS(N)%STR_SUB_INDX == HR%STR_SUB_INDX) THEN
               STRAIGHT_LINE_TO_TARGET = .TRUE.
               HR%I_TARGET = ABS(HR%I_TARGET)
               X_TARGET = (EVAC_EXITS(N)%X1 + EVAC_EXITS(N)%X2)/2.0_EB
               Y_TARGET = (EVAC_EXITS(N)%Y1 + EVAC_EXITS(N)%Y2)/2.0_EB
               DOOR_WIDTH = EVAC_EXITS(N)%WIDTH
               DOOR_IOR   = EVAC_EXITS(N)%IOR
               X1=EVAC_EXITS(N)%X1 ; X2=EVAC_EXITS(N)%X2
               Y1=EVAC_EXITS(N)%Y1 ; Y2=EVAC_EXITS(N)%Y2
            END IF
         ELSE IF (N>0) THEN
            IF (EVAC_DOORS(N)%STR_SUB_INDX == HR%STR_SUB_INDX) THEN
               STRAIGHT_LINE_TO_TARGET = .TRUE.
               HR%I_TARGET = ABS(HR%I_TARGET)
               X_TARGET = (EVAC_DOORS(N)%X1 + EVAC_DOORS(N)%X2)/2.0_EB
               Y_TARGET = (EVAC_DOORS(N)%Y1 + EVAC_DOORS(N)%Y2)/2.0_EB
               DOOR_WIDTH = EVAC_DOORS(N)%WIDTH
               DOOR_IOR   = EVAC_DOORS(N)%IOR
               X1=EVAC_DOORS(N)%X1 ; X2=EVAC_DOORS(N)%X2
               Y1=EVAC_DOORS(N)%Y1 ; Y2=EVAC_DOORS(N)%Y2
            END IF
         END IF
      END IF If_StrsMesh
      If_Strs_Mesh2: IF (STRAIGHT_LINE_TO_TARGET) THEN
         UBAR = X_TARGET-HR%X
         VBAR = Y_TARGET-HR%Y
         DOOR_DIST = SQRT((X_TARGET-HR%X)**2+(Y_TARGET-HR%Y)**2)
         SELECT CASE(DOOR_IOR)
         CASE(-1,+1)
            IF ((HR%Y > Y1) .AND. (HR%Y < Y2)) THEN
               IF (UBAR==0.0_EB) UBAR = HR%UBAR  ! used old one if at the door line
               HR%SKIP_WALL_FORCE_IOR = NINT(SIGN(1.0_EB,UBAR))
               IF ((HR%Y > Y1+HR%Radius) .AND. (HR%Y < Y2-HR%Radius)) THEN
                  UBAR = SIGN(1.0_EB,UBAR)
                  VBAR = 0._EB
               END IF
            END IF
         CASE(-2,+2)
            IF ((HR%X > X1) .AND. (HR%X < X2)) THEN
               IF (VBAR==0.0_EB) VBAR = HR%VBAR  ! used old one if at the door line
               HR%SKIP_WALL_FORCE_IOR = NINT(SIGN(1.0_EB,VBAR)*2)
               IF ((HR%X > X1+HR%Radius) .AND. (HR%X < X2-HR%Radius)) THEN
                  UBAR = 0._EB
                  VBAR = SIGN(1.0_EB,VBAR)
               END IF
            END IF
         END SELECT
      ELSE IF (NM_STRS_MESH) THEN ! If_Strs_Mesh2: In stairs but not straight line to target
         CALL STRS_U_AND_V(STRP,HR%STR_SUB_INDX,HR%X,HR%Y,HR%STRS_DIRECTION,UBAR,VBAR)
      ELSE ! If_Strs_Mesh2: Is not a Strs Mesh

         IF (HR%CROWBAR_READ_IN .AND. HR%T_LastRead_CB < T-5.0_EB*CROWBAR_DT_READ) THEN
            HR%CROWBAR_READ_IN = .FALSE. ! No more camera information available for this agent
         END IF
         IF (HR%CROWBAR_READ_IN .AND. .NOT. HR%CROWBAR_UPDATE_V0) RETURN
         IF (.NOT.HR%CROWBAR_READ_IN .AND. HR%CROWBAR_UPDATE_V0) THEN
            ! Read in agents that have lately "disappered" from the camera data at or close to the entr XB boundary
            ! Continue with the current UBAR,VBAR and use the Crowbar modified tau and speed
            UBAR = HR%UBAR; VBAR = HR%VBAR
            ! Update target velocity using camera information
            EVEL = MAX(0.01_EB,SQRT(HR%U_CB**2+HR%V_CB**2))
            HR_SPEED = MAX(0.05_EB*HR_SPEED,MIN(HR_SPEED,1.0_EB*EVEL))
            ! Slow speed => tau should be scaled also (motive force strenght scales as v0/tau)
            HR_TAU = MAX(0.1_EB,MIN(1.0_EB,HR_SPEED/HR%SPEED))*HR%TAU
            IF (HR%T_LastRead_CB < T-15.0_EB*CROWBAR_DT_READ) HR%CROWBAR_UPDATE_V0 = .FALSE.
            RETURN
         END IF
         ! IF (HR%CROWBAR_READ_IN .AND. HR%CROWBAR_UPDATE_V0 .AND. HR%T_LastRead_CB >= T-5.0_EB*CROWBAR_DT_READ) THEN
         IF (HR%CROWBAR_READ_IN .AND. HR%CROWBAR_UPDATE_V0) THEN
            ! Crowbar project: Use the read in camera information
            IF (T >= TPRE) THEN
               ! Update target velocity using camera information
               EVEL = MAX(0.01_EB,SQRT(HR%U_CB**2+HR%V_CB**2))
               HR_SPEED = MIN(HR_SPEED,1.0_EB*EVEL)
               UBAR = HR%U_CB/EVEL ; VBAR = HR%V_CB/EVEL
               ! Slow speed => tau should be scaled also (motive force strenght scales as v0/tau)
               HR_TAU = MAX(0.1_EB,MIN(1.0_EB,HR_SPEED/HR%SPEED))*HR%TAU
               ! Force strong motive force so that the agents follow the read in trajectories
               HR_TAU = MAX(0.05_EB,0.25*HR_TAU)
               ! X2_CB,Y2_CB: Agent's real (HR%X,HR%Y) postition at the time of last read in
               Dr_x = (HR%X2_CB-HR%X_CB) ; Dr_y = (HR%Y2_CB-HR%Y_CB)
               Dr_length = SQRT(Dr_x**2 + Dr_y**2)
               IF (EVEL > 0.1_EB .AND. Dr_length > 0.01_EB) THEN
                  ! Read in velocity is not too close to zero
                  D_DrV0 = (Dr_x*UBAR + Dr_y*VBAR) ! distance along the v0 direction
                  ! Agent should be slowed down to get it closer to camera position
                  ! If more than 1.0 m => max slow down that is now half speed
                  ! Agent should be accelerated to get it closer to camera position
                  ! If more than 1.0 m => max acc that is now 1.8*speed
                  ! y(x) = 1 - 2.0*x, where x is d_drv0
                  HR_SPEED = HR_SPEED*MIN(1.8_EB,MAX(0.2_EB,1.0_EB-5.0_EB*D_DrV0))

                  Dv0_x = D_DrV0*UBAR - Dr_x  ! Dv0 is a vector pointing normally towards the v0 direction
                  Dv0_y = D_DrV0*VBAR - Dr_y
                  Dv0 = SQRT(Dv0_x**2 + Dv0_y**2)
                  IF (Dv0 > 0.01_EB) THEN
                     Dv0_x = Dv0_x/Dv0 ; Dv0_y = Dv0_y/Dv0
                  ELSE
                     Dv0_x = 0.0_EB ; Dv0_y = 0.0_EB
                  END IF
                  alpha = Dv0/HR_SPEED ! normal distance in seconds
                  alpha = alpha/(2.0_EB*CROWBAR_DT_READ) ! maximum change in direction, if more than dt_read time difference
                  alpha = MIN(0.5_EB,alpha) ! no more than 45 deg change in the v0 direction
                  UBAR = (1.0_EB-alpha)*UBAR + alpha*Dv0_x ; VBAR = (1.0_EB-alpha)*VBAR + alpha*Dv0_y
                  EVEL = MAX(0.01_EB,SQRT(HR%U_CB**2+HR%V_CB**2))
                  UBAR = UBAR/EVEL ; VBAR = VBAR/EVEL  ! Should be a unit vector
               ELSEIF (EVEL <= 0.1_EB .AND. Dr_length > 1.0_EB*CROWBAR_DT_READ) THEN
                  ! Agent in the camera is not moving much and xy_camera differs xy_i more than
                  ! 1.0 m/s * dt_read metres then go towards the read in xy position.
                  Dr_x = (HR%X_CB-HR%X2_CB) ; Dr_y = (HR%Y_CB-HR%Y2_CB) ! Points towards the read in xy
                  Dr_length = SQRT(Dr_x**2 + Dr_y**2)
                  UBAR = Dr_x/MAX(0.01_EB,Dr_length); VBAR = Dr_y/MAX(0.01_EB,Dr_length)
                  HR_SPEED = HR_SPEED*MIN(1.8_EB,MAX(0.01_EB,1.0_EB-0.5_EB*Dr_length))
               END IF
            ELSE
               ! Initialization phase, i.e., initial positions from the camera data
               UBAR = 0.0_EB ; VBAR = 0.0_EB
            END IF
            RETURN
         END IF

         If_BeeLine: IF (HR%I_Door_Mode >= 2 .AND. HR%I_TARGET /= 0) THEN
            ! Use the flow field by default
            N = ABS(HR%I_TARGET)
            N_END = N_EXITS - N_CO_EXITS + N_DOORS
            IF (N > N_DOORS) THEN
               N = N - N_DOORS  ! The exit line number including also count only exits
               JJ_NOW = EVAC_EXITS(N)%I_EMESH_EXITS
               IF (TRIM(EVAC_EXITS(N)%ID) /= TRIM(EMESH_EXITS(JJ_NOW)%ID)) THEN
                  write(lu_err,*)'*** ERROR IN FIND_PREFERRED_DIRECTION'
                  write(lu_err,*)'*** ',TRIM(EVAC_EXITS(N)%ID), ' is not ',TRIM(EMESH_EXITS(JJ_NOW)%ID)
               END IF
            ELSE
               JJ_NOW = EVAC_DOORS(N)%I_EMESH_EXITS
               IF (TRIM(EVAC_DOORS(N)%ID) /= TRIM(EMESH_EXITS(JJ_NOW)%ID)) THEN
                  write(lu_err,*)'*** ERROR IN FIND_PREFERRED_DIRECTION'
                  write(lu_err,*)'*** ',TRIM(EVAC_DOORS(N)%ID), ' is not ',TRIM(EMESH_EXITS(JJ_NOW)%ID)
               END IF
            END IF
            ! JJ_NOW = ABS(HR%I_TARGET)  ! 1,....,n_doors+n_exits
            UBAR = (1.0_EB-(XI-II+1))*EMESH_EXITS(JJ_NOW)%U_EVAC(II-1,JJ) + (XI-II+1)*EMESH_EXITS(JJ_NOW)%U_EVAC(II,JJ)
            VBAR = (1.0_EB-(YJ-JJ+1))*EMESH_EXITS(JJ_NOW)%V_EVAC(II,JJ-1) + (YJ-JJ+1)*EMESH_EXITS(JJ_NOW)%V_EVAC(II,JJ)
            Is_InFront = .FALSE. ; Is_XB_Visible = .FALSE. ; Is_XYZ_Visible = .FALSE.
            IF (HR%I_Door_Mode > 2) Is_XB_Visible = .TRUE.
            IF (HR%I_Door_Mode == 2 .OR. HR%I_Door_Mode == 4) Is_XYZ_Visible = .TRUE.

            N = ABS(HR%I_TARGET)
            IF (N>N_DOORS) THEN
               N = N - N_DOORS
               IF (ABS(EVAC_EXITS(N)%IMODE)==2 .AND. EVAC_EXITS(N)%TARGET_WHEN_CLOSED .AND. &
                    EVAC_EXITS(N)%WAIT_AT_XYZ) THEN
                  X_XYZ = EVAC_EXITS(N)%X ; Y_XYZ = EVAC_EXITS(N)%Y
                  DOOR_IOR   = EVAC_EXITS(N)%IOR
                  DOOR_WIDTH = EVAC_EXITS(N)%WIDTH
                  IF (ABS(DOOR_IOR)==1) THEN
                     IF (DOOR_IOR*(HR%X-X_XYZ) >= TWO_EPSILON_EB) THEN
                        UBAR = 0.0_EB ; VBAR = 0.0_EB
                        V0_IS_SET_ZERO = .TRUE.
                        Is_InFront = .FALSE. ; Is_XB_Visible = .FALSE. ; Is_XYZ_Visible = .FALSE.
                     END IF
                  ELSE
                     IF (DOOR_IOR*(HR%Y-Y_XYZ) >= TWO_EPSILON_EB) THEN
                        UBAR = 0.0_EB ; VBAR = 0.0_EB
                        V0_IS_SET_ZERO = .TRUE.
                        Is_InFront = .FALSE. ; Is_XB_Visible = .FALSE. ; Is_XYZ_Visible = .FALSE.
                     END IF
                  END IF
               END IF
            ELSE
               IF (ABS(EVAC_DOORS(N)%IMODE)==2 .AND. EVAC_DOORS(N)%TARGET_WHEN_CLOSED .AND. &
                    EVAC_DOORS(N)%WAIT_AT_XYZ) THEN
                  X_XYZ = EVAC_DOORS(N)%X ; Y_XYZ = EVAC_DOORS(N)%Y
                  DOOR_IOR   = EVAC_DOORS(N)%IOR
                  DOOR_WIDTH = EVAC_DOORS(N)%WIDTH
                  IF (ABS(DOOR_IOR)==1) THEN
                     IF (DOOR_IOR*(HR%X-X_XYZ) >= TWO_EPSILON_EB) THEN
                        UBAR = 0.0_EB ; VBAR = 0.0_EB
                        V0_IS_SET_ZERO = .TRUE.
                        Is_InFront = .FALSE. ; Is_XB_Visible = .FALSE. ; Is_XYZ_Visible = .FALSE.
                     END IF
                  ELSE
                     IF (DOOR_IOR*(HR%Y-Y_XYZ) >= TWO_EPSILON_EB) THEN
                        UBAR = 0.0_EB ; VBAR = 0.0_EB
                        V0_IS_SET_ZERO = .TRUE.
                        Is_InFront = .FALSE. ; Is_XB_Visible = .FALSE. ; Is_XYZ_Visible = .FALSE.
                     END IF
                  END IF
               END IF
            END IF

            N = ABS(HR%I_TARGET)
            IF (Is_XB_Visible .OR. Is_XYZ_Visible) THEN
               IF (N>N_DOORS) THEN
                  N = N - N_DOORS
                  X_XYZ = EVAC_EXITS(N)%X ; Y_XYZ = EVAC_EXITS(N)%Y
                  X1 = EVAC_EXITS(N)%X1 ; X2 = EVAC_EXITS(N)%X2
                  Y1 = EVAC_EXITS(N)%Y1 ; Y2 = EVAC_EXITS(N)%Y2
                  DOOR_IOR = EVAC_EXITS(N)%IOR
                  DOOR_WIDTH = EVAC_EXITS(N)%WIDTH
               ELSE
                  X_XYZ = EVAC_DOORS(N)%X ; Y_XYZ = EVAC_DOORS(N)%Y
                  X1 = EVAC_DOORS(N)%X1 ; X2 = EVAC_DOORS(N)%X2
                  Y1 = EVAC_DOORS(N)%Y1 ; Y2 = EVAC_DOORS(N)%Y2
                  DOOR_IOR = EVAC_DOORS(N)%IOR
                  DOOR_WIDTH = EVAC_DOORS(N)%WIDTH
               END IF

               If_IOR_Is_1: IF (ABS(DOOR_IOR)==1)THEN
                  Y1 = Y1+MIN(0.3_EB,0.5_EB*DOOR_WIDTH)
                  Y2 = Y2-MIN(0.3_EB,0.5_EB*DOOR_WIDTH)
                  Y_TARGET = Y_XYZ
                  D_X = ABS(HR%X-X_XYZ) ; D_Y = ABS(HR%Y-Y_XYZ)
                  EVEL = ABS(ATAN(D_Y/MAX(D_X,0.01_EB)))/(0.5_EB*PI)  ! 0-1 scale for 0-90 degrees
                  IF (HR%Y > Y_XYZ - 0.5_EB*DOOR_WIDTH + MIN(0.3_EB,0.5_EB*DOOR_WIDTH) .AND. &
                       HR%Y < Y_XYZ + 0.5_EB*DOOR_WIDTH - MIN(0.3_EB,0.5_EB*DOOR_WIDTH)) Y_TARGET = HR%Y
                  IF (HR%Y < Y_XYZ - 0.5_EB*DOOR_WIDTH + MIN(0.3_EB,0.5_EB*DOOR_WIDTH) .AND. HR%SeeDoorXYZ1) &
                       Y_TARGET = Y_XYZ - EVEL*(0.5_EB*DOOR_WIDTH - MIN(0.3_EB,0.5_EB*DOOR_WIDTH))
                  IF (HR%Y > Y_XYZ + 0.5_EB*DOOR_WIDTH - MIN(0.3_EB,0.5_EB*DOOR_WIDTH) .AND. HR%SeeDoorXYZ2) &
                       Y_TARGET = Y_XYZ + EVEL*(0.5_EB*DOOR_WIDTH - MIN(0.3_EB,0.5_EB*DOOR_WIDTH))
                  Y_XYZ = Y_TARGET
                  D_TMP = DOOR_IOR*(X_XYZ-HR%X) ! if < 0 then between door and xyz
                  X_TARGET = 0.5_EB*(X1+X2) ; Y_TARGET = 0.5_EB*(Y1+Y2)
                  D_X = ABS(HR%X-0.5_EB*(X1+X2)) ; D_Y = ABS(HR%Y-0.5_EB*(Y1+Y2))
                  EVEL = ABS(ATAN(D_Y/MAX(D_X,0.01_EB)))/(0.5_EB*PI)  ! 0-1 scale for 0-90 degrees
                  IF (HR%Y > Y1 .AND. HR%Y < Y2) Is_InFront = .TRUE.
                  IF (HR%Y < Y1 .AND. HR%SeeDoorXB1) &
                       Y_TARGET = 0.5_EB*(Y1+Y2) - EVEL*(0.5_EB*DOOR_WIDTH - MIN(0.3_EB,0.5_EB*DOOR_WIDTH))
                  IF (HR%Y > Y2 .AND. HR%SeeDoorXB2) &
                       Y_TARGET = 0.5_EB*(Y1+Y2) + EVEL*(0.5_EB*DOOR_WIDTH - MIN(0.3_EB,0.5_EB*DOOR_WIDTH))
               ELSE ! abs(ior)=2 branch, If_IOR_Is_1
                  X1 = X1+MIN(0.3_EB,0.5_EB*DOOR_WIDTH) ; X2 = X2-MIN(0.3_EB,0.5_EB*DOOR_WIDTH)
                  X_TARGET = X_XYZ
                  D_X = ABS(HR%X-X_XYZ) ; D_Y = ABS(HR%Y-Y_XYZ)
                  EVEL = ABS(ATAN(D_X/MAX(D_Y,0.01_EB)))/(0.5_EB*PI)  ! 0-1 scale for 0-90 degrees
                  IF (HR%X > X_XYZ - 0.5_EB*DOOR_WIDTH + MIN(0.3_EB,0.5_EB*DOOR_WIDTH) .AND. &
                       HR%X < X_XYZ + 0.5_EB*DOOR_WIDTH - MIN(0.3_EB,0.5_EB*DOOR_WIDTH)) X_TARGET = HR%X
                  IF (HR%X < X_XYZ - 0.5_EB*DOOR_WIDTH + MIN(0.3_EB,0.5_EB*DOOR_WIDTH) .AND. HR%SeeDoorXYZ1) &
                       X_TARGET = X_XYZ - EVEL*(0.5_EB*DOOR_WIDTH - MIN(0.3_EB,0.5_EB*DOOR_WIDTH))
                  IF (HR%X > X_XYZ + 0.5_EB*DOOR_WIDTH - MIN(0.3_EB,0.5_EB*DOOR_WIDTH) .AND. HR%SeeDoorXYZ2) &
                       X_TARGET = X_XYZ + EVEL*(0.5_EB*DOOR_WIDTH - MIN(0.3_EB,0.5_EB*DOOR_WIDTH))
                  X_XYZ = X_TARGET
                  D_TMP = 0.5_EB*DOOR_IOR*(Y_XYZ-HR%Y) ! if < 0 then between door and xyz
                  Y_TARGET = 0.5_EB*(Y1+Y2) ; X_TARGET = 0.5_EB*(X1+X2)
                  D_X = ABS(HR%X-0.5_EB*(X1+X2)) ; D_Y = ABS(HR%Y-0.5_EB*(Y1+Y2))
                  EVEL = ABS(ATAN(D_X/MAX(D_Y,0.01_EB)))/(0.5_EB*PI)  ! 0-1 scale for 0-90 degrees
                  IF (HR%X > X1 .AND. HR%X < X2) Is_InFront = .TRUE.
                  IF (HR%X < X1 .AND. HR%SeeDoorXB1) &
                       X_TARGET = 0.5_EB*(X1+X2) - EVEL*(0.5_EB*DOOR_WIDTH - MIN(0.3_EB,0.5_EB*DOOR_WIDTH))
                  IF (HR%X > X2 .AND. HR%SeeDoorXB2) &
                       X_TARGET = 0.5_EB*(X1+X2) + EVEL*(0.5_EB*DOOR_WIDTH - MIN(0.3_EB,0.5_EB*DOOR_WIDTH))
               END IF If_IOR_Is_1

               IF (D_TMP <= 0.0_EB) Is_XYZ_Visible = .FALSE.  ! use door XB only
               DOOR_DIST = SQRT((X_XYZ-HR%X)**2+(Y_XYZ-HR%Y)**2)
               EVEL = MAX(DOOR_WIDTH,SQRT((X_XYZ-0.5_EB*(X1+X2))**2+(Y_XYZ-0.5_EB*(Y1+Y2))**2))
               IF (DOOR_DIST <= EVEL .AND. Is_XB_Visible) Is_XYZ_Visible = .FALSE.  ! use door XB only
               EVEL = 1.5_EB*DOOR_WIDTH
               IF (DOOR_DIST <= EVEL .AND. .NOT.Is_XB_Visible) Is_XYZ_Visible = .FALSE.  ! use flow field

               ! Calculate the angle (left or right) to the door relative to the v0 direction
               EVEL = SQRT(UBAR**2 + VBAR**2)
               IF (EVEL >= TWO_EPSILON_EB) THEN
                  EVEL = MAX(0.001_EB,SQRT(UBAR**2 + VBAR**2))
                  IF (VBAR >= 0.0_EB) THEN
                     ANGLE_HR = ACOS(UBAR/EVEL)
                  ELSE
                     ANGLE_HR = 2.0_EB*PI - ACOS(UBAR/EVEL)
                  END IF
                  IF (ANGLE_HR == 2.0_EB*PI) ANGLE_HR = 0.0_EB  ! Agent HR angle is [0,2PI)
               ELSE
                  ANGLE_HR = HR%ANGLE
               END IF

               D_NEW = SQRT((X_XYZ-HR%X)**2 + (Y_XYZ-HR%Y)**2)
               IF ((Y_XYZ-HR%Y) >= 0.0_EB) THEN
                  ANGLE_XYZ = ACOS((X_XYZ-HR%X)/D_NEW)
               ELSE
                  ANGLE_XYZ = 2.0_EB*PI - ACOS((X_XYZ-HR%X)/D_NEW)
               END IF
               IF (ANGLE_XYZ == 2.0_EB*PI) ANGLE_XYZ = 0.0_EB  ! Agent HRE angle is [0,2PI)
               ANGLE_XYZ = ANGLE_HR - ANGLE_XYZ
               IF (ANGLE_XYZ >= PI) ANGLE_XYZ = 2.0_EB*PI - ANGLE_XYZ
               IF (ANGLE_XYZ <= -PI) ANGLE_XYZ = 2.0_EB*PI + ANGLE_XYZ
               ANGLE_XYZ = -180.0_EB*ANGLE_XYZ/PI  ! Degrees, -180 - +180
               ! If XYZ is on the right then the angle is negative, i.e., positive direction is anti-clockwise.
               IF (ANGLE_XYZ < -20.0_EB .AND. HR%DensityR > 1.0_EB .AND. Is_XYZ_Visible) THEN
                  Is_XYZ_Visible = .FALSE.
               END IF
               IF (ANGLE_XYZ >  20.0_EB .AND. HR%DensityL > 1.0_EB .AND. Is_XYZ_Visible) THEN
                  Is_XYZ_Visible = .FALSE.
               END IF

               D_NEW = SQRT((0.5_EB*(X1+X2)-HR%X)**2 + (0.5_EB*(Y1+Y2)-HR%Y)**2)
               IF ((0.5_EB*(Y1+Y2)-HR%Y) >= 0.0_EB) THEN
                  ANGLE_XB = ACOS((0.5_EB*(X1+X2)-HR%X)/D_NEW)
               ELSE
                  ANGLE_XB = 2.0_EB*PI - ACOS((0.5_EB*(X1+X2)-HR%X)/D_NEW)
               END IF
               IF (ANGLE_XB == 2.0_EB*PI) ANGLE_XB = 0.0_EB  ! Agent HRE angle is [0,2PI)
               ANGLE_XB = ANGLE_HR - ANGLE_XB
               IF (ANGLE_XB >= PI) ANGLE_XB = 2.0_EB*PI - ANGLE_XB
               IF (ANGLE_XB <= -PI) ANGLE_XB = 2.0_EB*PI + ANGLE_XB
               ANGLE_XB = -180.0_EB*ANGLE_XB/PI  ! Degrees, -180 - +180
               ! If XB is on the right then the angle is negative, i.e., positive direction is anti-clockwise.
               DOOR_DIST = SQRT((0.5_EB*(X1+X1)-HR%X)**2+(0.5_EB*(Y1+Y1)-HR%Y)**2)
               IF (ANGLE_XB < -20.0_EB .AND. HR%DensityR > 1.0_EB .AND. Is_XB_Visible .AND. DOOR_DIST > DOOR_WIDTH) THEN
                  Is_XB_Visible = .FALSE.
               END IF
               IF (ANGLE_XB >  20.0_EB .AND. HR%DensityL > 1.0_EB .AND. Is_XB_Visible .AND. DOOR_DIST > DOOR_WIDTH) THEN
                  Is_XB_Visible = .FALSE.
               END IF
               !IF (ANGLE_XB < -0.0_EB .AND. HR%DensityR > 1.0_EB .AND. HR%DensityL < 0.5_EB) Is_XB_Visible = .FALSE.
               !IF (ANGLE_XB >  0.0_EB .AND. HR%DensityL > 1.0_EB .AND. HR%DensityR < 0.5_EB) Is_XB_Visible = .FALSE.

               IF (D_TMP > 0 .AND. Is_XYZ_Visible) THEN
                  IF (Is_InFront .AND. Is_XB_Visible) THEN
                     ! perpendicular to the door plane, i.e., straight
                     UBAR = REAL(MOD(DOOR_IOR,2),EB) ; VBAR = REAL(DOOR_IOR/2,EB)
                  ELSE
                     ! Aim towards XYZ
                     UBAR = X_XYZ-HR%X ; VBAR = Y_XYZ-HR%Y
                  END IF
               ELSEIF (Is_XB_Visible) THEN
                  ! D_TMP>0 and xyz not visible .OR. D_TMP<0 for all cases (d_tmp<0: between xb and xyz)
                  ! Nothing visible is corrected below by overwriting ubar,vbar
                  IF (Is_InFront) THEN
                     ! perpendicular to the door plane, i.e., straight
                     UBAR = REAL(MOD(DOOR_IOR,2),EB) ; VBAR = REAL(DOOR_IOR/2,EB)
                  ELSE
                     ! Aim towards XB, towards the "door posts", if not in front
                     UBAR = X_TARGET-HR%X ; VBAR = Y_TARGET-HR%Y
                  END IF
               END IF
            END IF
         ELSE  ! If_BeeLine (or not fds6)
            ! Normal case, use the flow field
            N = ABS(HR%I_TARGET)
            N_END = N_EXITS - N_CO_EXITS + N_DOORS
            IF (N==0) THEN
               DEFAULT_TARGET: DO JJ_NOW = 1, N_END
                  IF (EMESH_EXITS(JJ_NOW)%MAINMESH == NM .AND. EMESH_EXITS(JJ_NOW)%DEFINE_MESH) THEN
                     IF (HR%I_Door_Mode < 0 .AND. .NOT.EMESH_EXITS(JJ_NOW)%IS_EXIT) THEN
                        IF (EVAC_NODE_LIST(ABS(HR%I_Door_Mode))%NODE_TYPE == 'Door') THEN
                           IF (Trim(EVAC_NODE_LIST(ABS(HR%I_Door_Mode))%ID) == TRIM(EMESH_EXITS(JJ_NOW)%ID)) &
                                CYCLE DEFAULT_TARGET
                        END IF
                     END IF
                     N = JJ_NOW
                     EXIT DEFAULT_TARGET
                  END IF
               END DO DEFAULT_TARGET
               IF (EMESH_EXITS(N)%IS_EXIT) THEN
                  N = EMESH_EXITS(N)%NEXIT + N_DOORS ! I_Target index
               ELSE
                  N = EMESH_EXITS(N)%NEXIT  ! I_Target index
               END IF
            END IF

            IF (N > N_DOORS) THEN
               N = N - N_DOORS  ! The exit line number including also count only exits
               JJ_NOW = EVAC_EXITS(N)%I_EMESH_EXITS
               IF (TRIM(EVAC_EXITS(N)%ID) /= TRIM(EMESH_EXITS(JJ_NOW)%ID)) THEN
                  WRITE(LU_ERR,*)'*** ERROR IN FINE_PREFERRED_DIRECTION'
                  WRITE(LU_ERR,*)'*** ',TRIM(EVAC_DOORS(N)%ID), ' is not ',TRIM(EMESH_EXITS(JJ_NOW)%ID)
               END IF
            ELSEIF (N > 0) THEN
               JJ_NOW = EVAC_DOORS(N)%I_EMESH_EXITS
               IF (TRIM(EVAC_DOORS(N)%ID) /= TRIM(EMESH_EXITS(JJ_NOW)%ID)) THEN
                  WRITE(LU_ERR,*)'*** ERROR IN FINE_PREFERRED_DIRECTION'
                  WRITE(LU_ERR,*)'*** ',TRIM(EVAC_DOORS(N)%ID), ' is not ',TRIM(EMESH_EXITS(JJ_NOW)%ID)
               END IF
            END IF
            IF (EMESH_EXITS(JJ_NOW)%MAINMESH /= NM) THEN
               WRITE(LU_ERR,*)'*** find direction error: ',JJ_NOW,NM, &
                    EMESH_EXITS(JJ_NOW)%MAINMESH,' ',TRIM(EMESH_EXITS(JJ_NOW)%ID)
            END IF
            UBAR = (1.0_EB-(XI-II+1))*EMESH_EXITS(JJ_NOW)%U_EVAC(II-1,JJ) + (XI-II+1)*EMESH_EXITS(JJ_NOW)%U_EVAC(II,JJ)
            VBAR = (1.0_EB-(YJ-JJ+1))*EMESH_EXITS(JJ_NOW)%V_EVAC(II,JJ-1) + (YJ-JJ+1)*EMESH_EXITS(JJ_NOW)%V_EVAC(II,JJ)
         END IF If_BeeLine
      END IF If_Strs_Mesh2

      EVEL = SQRT(UBAR**2 + VBAR**2)  ! (UBAR,VBAR) is an unit vector
      IF (EVEL >= TWO_EPSILON_EB) THEN
         UBAR = UBAR/EVEL ; VBAR = VBAR/EVEL
      ELSE
         ! If no v0 found for the current position of the agent, use the previous value of the agent
         IF (.NOT.V0_IS_SET_ZERO) THEN
            UBAR = HR%UBAR ; VBAR = HR%VBAR
         END IF
      END IF
      IF (L_DEAD .OR. L_FALLEN_DOWN) THEN
         UBAR = 0.0_EB; VBAR = 0.0_EB
      END IF
      IF (I_HERDING_TYPE>1 .AND. HR%I_Target==0 .AND. HR%I_DoorAlgo>0 .AND. .NOT.V0_IS_SET_ZERO) THEN
         ! Herding type agent without any target door: do not move
         ! Agent without any target door: do not move
         DIST_TO_DOOR_TMP = 0.0_EB
         IF (HR%I_Door_Mode < 0) THEN
            IF (EVAC_NODE_LIST(ABS(HR%I_Door_Mode))%NODE_TYPE == 'Door') THEN
               DIST_TO_DOOR_TMP = SQRT((HR%X - 0.5_EB*(EVAC_DOORS(EVAC_NODE_LIST(ABS(HR%I_Door_Mode))%NODE_INDEX)%X1 + &
                    EVAC_DOORS(EVAC_NODE_LIST(ABS(HR%I_Door_Mode))%NODE_INDEX)%X2))**2 + &
                    (HR%Y - 0.5_EB*(EVAC_DOORS(EVAC_NODE_LIST(ABS(HR%I_Door_Mode))%NODE_INDEX)%Y1 + &
                    EVAC_DOORS(EVAC_NODE_LIST(ABS(HR%I_Door_Mode))%NODE_INDEX)%Y2))**2)
            END IF
            IF (EVAC_NODE_LIST(ABS(HR%I_Door_Mode))%NODE_TYPE == 'Entry') THEN
               DIST_TO_DOOR_TMP = SQRT((HR%X - 0.5_EB*(EVAC_ENTRYS(EVAC_NODE_LIST(ABS(HR%I_Door_Mode))%NODE_INDEX)%X1 + &
                    EVAC_ENTRYS(EVAC_NODE_LIST(ABS(HR%I_Door_Mode))%NODE_INDEX)%X2))**2 + &
                    (HR%Y - 0.5_EB*(EVAC_ENTRYS(EVAC_NODE_LIST(ABS(HR%I_Door_Mode))%NODE_INDEX)%Y1 + &
                    EVAC_ENTRYS(EVAC_NODE_LIST(ABS(HR%I_Door_Mode))%NODE_INDEX)%Y2))**2)
            END IF
            IF (DIST_TO_DOOR_TMP > 3.0_EB) HR%I_Door_Mode = 0
            ! Herding agents that do not have any target will stop after 3 m from the
            ! door/entr they use to came to this floor.
         ELSE
            UBAR = 0.0_EB ; VBAR = 0.0_EB
         END IF
      END IF
      J = MAX(0,HR%GROUP_ID)
      IF (J == 0 ) THEN
         GROUP_LIST(0)%TPRE = HR%TPRE + HR%TDET
         TPRE = HR%TPRE + HR%TDET
      ELSE
         TPRE = HR%TDET
      END IF
      IF (L_DEAD) TPRE = HUGE(TPRE)

      ! Direction to the centre of the group (if any)
      IF (GROUP_LIST(J)%GROUP_SIZE >= 2) THEN
         HR%UBAR_CENTER = (GROUP_LIST(J)%GROUP_X - HR%X)
         HR%VBAR_CENTER = (GROUP_LIST(J)%GROUP_Y - HR%Y)
         EVEL = SQRT(HR%UBAR_CENTER**2 + HR%VBAR_CENTER**2)
         IF ( EVEL >= TWO_EPSILON_EB .AND. .NOT. L_DEAD ) THEN
            HR%UBAR_CENTER = HR%UBAR_CENTER / EVEL
            HR%VBAR_CENTER = HR%VBAR_CENTER / EVEL
         ELSE
            HR%UBAR_CENTER = 0.0_EB ; HR%VBAR_CENTER = 0.0_EB
         END IF
      ELSE
         ! Only one person in the group
         HR%UBAR_CENTER = 0.0_EB ; HR%VBAR_CENTER = 0.0_EB
      END IF

      IF ( J > 0 ) THEN
         ! Group is already gathered together, but not yet moving towards the door
         IF (GROUP_LIST(J)%COMPLETE == 1 .AND. T <= GROUP_LIST(J)%TPRE+GROUP_LIST(J)%TDOOR) THEN
            UBAR = 0.0_EB ; VBAR = 0.0_EB
         END IF

         EVEL = UBAR**2 + VBAR**2
         IF (GROUP_LIST(J)%COMPLETE == 1 .AND. EVEL >= TWO_EPSILON_EB) THEN
            ! The group is already gathered together
            UBAR = (1-GROUP_EFF)*UBAR + GROUP_EFF*HR%UBAR_CENTER/ &
                 SQRT(((1-GROUP_EFF)*UBAR + GROUP_EFF*HR%UBAR_CENTER)**2+ &
                 ((1-GROUP_EFF)*VBAR + GROUP_EFF*HR%VBAR_CENTER)**2)
            VBAR = (1-GROUP_EFF)*VBAR + GROUP_EFF*HR%VBAR_CENTER/ &
                 SQRT(((1-GROUP_EFF)*UBAR + GROUP_EFF*HR%UBAR_CENTER)**2+ &
                 ((1-GROUP_EFF)*VBAR + GROUP_EFF*HR%VBAR_CENTER)**2)
         ELSE
            ! The group is still in the gathering phase
            UBAR = HR%UBAR_CENTER ; VBAR = HR%VBAR_CENTER
         END IF
      END IF

      EVEL = UBAR**2 + VBAR**2  ! (UBAR,VBAR) is an unit vector
      IF  (T <= TPRE .OR. EVEL < 0.01_EB) THEN
         ! No movement yet or no direction
         UBAR = 0.0_EB ; VBAR = 0.0_EB
      END IF
      RETURN
    END SUBROUTINE FIND_PREFERRED_DIRECTION


    SUBROUTINE GETSTAIRSPEEDANDZ(SPEED_XM, SPEED_XP, SPEED_YM, SPEED_YP, HP_SPEED, SP, HP)
      IMPLICIT NONE
      !
      ! Passed variables
      REAL(EB), INTENT(OUT) :: SPEED_XM, SPEED_XP, SPEED_YM, SPEED_YP
      REAL(EB), INTENT(IN)  :: HP_SPEED
      TYPE (EVAC_STRS_TYPE), POINTER ::  SP
      TYPE (HUMAN_TYPE), POINTER :: HP
      !
      ! Local variables
      REAL(EB) COS_X, COS_Y
      INTEGER J1, J2, J

      IF (HP%STR_SUB_INDX > 0) THEN
         J1 = MAX(HP%STR_SUB_INDX-1,1)
         J2 = MIN(HP%STR_SUB_INDX+1,STRP%N_NODES)
      ELSE
         J1 = 1
         J2 = STRP%N_NODES
      END IF

      LOOPSUBINDX: DO J = J1,J2
         IF (HP%X < SP%XB_NODE(J,1)) CYCLE
         IF (HP%X > SP%XB_NODE(J,2)) CYCLE
         IF (HP%Y < SP%XB_NODE(J,3)) CYCLE
         IF (HP%Y > SP%XB_NODE(J,4)) CYCLE
         IF (HP%Z < SP%XB_NODE(J,5)) CYCLE
         IF (HP%Z > SP%XB_NODE(J,6)) CYCLE
         ! Ok, human is inside the subnode J
         HP%STR_SUB_INDX = J
         IF (SP%NODE_TYPE(J)==STRS_STAIR_TYPE) THEN
            COS_X = SP%XB_NODE(J,7)
            COS_Y = SP%XB_NODE(J,8)
            SELECT CASE (SP%NODE_IOR(J))
            CASE(-1)
               HP%Z = SP%XB_NODE(J,5) + (SP%XB_NODE(J,6)-SP%XB_NODE(J,5))* &
                    ABS(SP%XB_NODE(J,2)-HP%X)/ABS(SP%XB_NODE(J,2)-SP%XB_NODE(J,1))
            CASE(+1)
               HP%Z = SP%XB_NODE(J,5) + (SP%XB_NODE(J,6)-SP%XB_NODE(J,5))* &
                    ABS(SP%XB_NODE(J,1)-HP%X)/ABS(SP%XB_NODE(J,2)-SP%XB_NODE(J,1))
            CASE(-2)
               HP%Z = SP%XB_NODE(J,5) + (SP%XB_NODE(J,6)-SP%XB_NODE(J,5))* &
                    ABS(SP%XB_NODE(J,4)-HP%Y)/ABS(SP%XB_NODE(J,4)-SP%XB_NODE(J,3))
            CASE(+2)
               HP%Z = SP%XB_NODE(J,5) + (SP%XB_NODE(J,6)-SP%XB_NODE(J,5))* &
                    ABS(SP%XB_NODE(J,3)-HP%Y)/ABS(SP%XB_NODE(J,4)-SP%XB_NODE(J,3))
            END SELECT
            !            SELECT CASE (SP%NODE_IOR(J)*HP%STRS_DIRECTION)
            SELECT CASE (SP%NODE_IOR(J))
            CASE(-1)
               SPEED_XM = COS_X* HP_SPEED* SP%FAC_V0_UP
               SPEED_XP = COS_X* HP_SPEED* SP%FAC_V0_DOWN
               SPEED_YM =        HP_SPEED* SP%FAC_V0_HORI
               SPEED_YP =        HP_SPEED* SP%FAC_V0_HORI
            CASE(+1)
               SPEED_XM = COS_X* HP_SPEED* SP%FAC_V0_DOWN
               SPEED_XP = COS_X* HP_SPEED* SP%FAC_V0_UP
               SPEED_YM =        HP_SPEED* SP%FAC_V0_HORI
               SPEED_YP =        HP_SPEED* SP%FAC_V0_HORI
            CASE(-2)
               SPEED_XM =        HP_SPEED* SP%FAC_V0_HORI
               SPEED_XP =        HP_SPEED* SP%FAC_V0_HORI
               SPEED_YM = COS_Y* HP_SPEED* SP%FAC_V0_UP
               SPEED_YP = COS_Y* HP_SPEED* SP%FAC_V0_DOWN
            CASE(+2)
               SPEED_XM =        HP_SPEED* SP%FAC_V0_HORI
               SPEED_XP =        HP_SPEED* SP%FAC_V0_HORI
               SPEED_YM = COS_Y* HP_SPEED* SP%FAC_V0_DOWN
               SPEED_YP = COS_Y* HP_SPEED* SP%FAC_V0_UP
            END SELECT
         ELSE
            HP%Z = 0.5_EB*(SP%XB_NODE(J,5)+SP%XB_NODE(J,6))
         END IF
         EXIT LOOPSUBINDX
      END DO LOOPSUBINDX
    END SUBROUTINE GETSTAIRSPEEDANDZ

    SUBROUTINE FIND_TARGET_NODE_IN_STRS(SP,HP)
      IMPLICIT NONE
      !
      ! This subroutine sets
      !   HP%I_TARGET: the target door/exit index (now <0, not visible status is returned)
      !   HP%STRS_DIRECTION: +1 or -1
      !
      ! Passed variables
      TYPE (EVAC_STRS_TYPE), POINTER :: SP
      TYPE (HUMAN_TYPE), POINTER :: HP
      !
      ! Local variables
      LOGICAL ISKNOWNDOOR, FINALTARGETFOUND, SOMEKNOWNDOORS
      INTEGER :: I, ID, FINAL_NODE, IG, IN, INODE, STR_SUB_INDX
      REAL(EB) :: Z_NODE, Z_FINAL, DZ_NODE, DZ_FINAL, Z_FINAL_UNKNOWN,DZ_TMP1, DZ_TMP2, DZ_NODE_ACTUAL
      REAL(EB) :: DIST_TO_DOOR, X_NODE, Y_NODE, DIST_TO_DOOR_TMP
      TYPE (EVAC_ENTR_TYPE), POINTER :: PNX =>NULL()

      FINALTARGETFOUND = .FALSE.
      SOMEKNOWNDOORS = .FALSE.
      IG = ABS(HP%GROUP_ID)

      Z_FINAL_UNKNOWN = 0._EB
      DZ_NODE_ACTUAL = 0._EB

      DO I = 1, N_NODES
         IF (EVAC_NODE_LIST(I)%NODE_TYPE=='Door') THEN
            IF (EVAC_DOORS(EVAC_NODE_LIST(I)%NODE_INDEX)%KNOWN_DOOR) SOMEKNOWNDOORS=.TRUE.
         END IF
         IF (EVAC_NODE_LIST(I)%NODE_TYPE=='Exit') THEN
            IF (EVAC_EXITS(EVAC_NODE_LIST(I)%NODE_INDEX)%KNOWN_DOOR) SOMEKNOWNDOORS=.TRUE.
         END IF
      END DO
      DO I = 1, N_NODES
         ISKNOWNDOOR = .FALSE.
         IF (EVAC_NODE_LIST(I)%NODE_TYPE=='Door') THEN
            IF (EVAC_DOORS(EVAC_NODE_LIST(I)%NODE_INDEX)%KNOWN_DOOR) ISKNOWNDOOR=.TRUE.
         END IF
         IF (EVAC_NODE_LIST(I)%NODE_TYPE=='Exit') THEN
            IF (EVAC_EXITS(EVAC_NODE_LIST(I)%NODE_INDEX)%KNOWN_DOOR) ISKNOWNDOOR=.TRUE.
         END IF
         IF (IG == 0 ) THEN
            ! Human came from entry
            PNX => EVAC_ENTRYS(ABS(HR%IEL))
            IF (PNX%N_VENT_FFIELDS == 0 .AND. .NOT.SOMEKNOWNDOORS) THEN
               ISKNOWNDOOR = .TRUE.
            ELSE
               DO ID = 1, PNX%N_VENT_FFIELDS
                  IF (I == PNX%I_DOOR_NODES(ID)) THEN
                     ISKNOWNDOOR = .TRUE.
                     EXIT
                  END IF
               END DO
            ENDIF
         ELSE
            ! Check the group know doors
            IF (HUMAN_KNOWN_DOORS(IG)%N_NODES == 0 .AND. .NOT.SOMEKNOWNDOORS) THEN
               ISKNOWNDOOR = .TRUE.
            ELSE
               IN = EVAC_NODE_LIST(I)%NODE_INDEX
               DO ID = 1, HUMAN_KNOWN_DOORS(IG)%N_NODES
                  IF (I == ABS(HUMAN_KNOWN_DOORS(IG)%I_NODES(ID))) THEN
                     ISKNOWNDOOR = .TRUE.
                     EXIT
                  END IF
               END DO
            END IF
         END IF
         IF (EVAC_NODE_LIST(I)%NODE_TYPE=='Exit') THEN
            Z_FINAL_UNKNOWN = 0.5_EB*(EVAC_EXITS(EVAC_NODE_LIST(I)%NODE_INDEX)%Z1+EVAC_EXITS(EVAC_NODE_LIST(I)%NODE_INDEX)%Z2)
            IF (ISKNOWNDOOR) THEN
               FINAL_NODE = I
               Z_FINAL = Z_FINAL_UNKNOWN
               FINALTARGETFOUND = .TRUE.
               EXIT
            END IF
         END IF
      END DO
      IF (.NOT.FINALTARGETFOUND) Z_FINAL = Z_FINAL_UNKNOWN
      DZ_FINAL = Z_FINAL - HP%Z

      ! Find the target among the nodes leading out of the strs (doors and exits)
      FINALTARGETFOUND = .FALSE.
      DZ_TMP2 = HUGE(1.0_EB)
      DIST_TO_DOOR = HUGE(DIST_TO_DOOR)
      FINDTARGETNODELOOP: DO I = 1,SP%N_NODES_OUT
         INODE = SP%NODES_OUT(I)
         DZ_NODE = -1._EB * SIGN(1.0_EB,DZ_FINAL) ! Initialize dz_node to different direction than final target
         SELECT CASE(EVAC_NODE_LIST(INODE)%NODE_TYPE)
         CASE('Door')
            Z_NODE = EVAC_DOORS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%Z1
            X_NODE = 0.5_EB*(EVAC_DOORS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%X1+EVAC_DOORS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%X2)
            Y_NODE = 0.5_EB*(EVAC_DOORS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%Y1+EVAC_DOORS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%Y2)
         CASE('Exit')
            Z_NODE = EVAC_EXITS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%Z1
            X_NODE = 0.5_EB*(EVAC_EXITS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%X1+EVAC_EXITS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%X2)
            Y_NODE = 0.5_EB*(EVAC_EXITS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%Y1+EVAC_EXITS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%Y2)
         CASE DEFAULT
            WRITE(LU_ERR,*) 'ERROR (DEBUG): Unknown node type in find_target_node_in_strs'
         END SELECT
         DZ_NODE = Z_NODE - HP%Z
         IF (SIGN(1.0_EB,DZ_NODE)==SIGN(1.0_EB,DZ_FINAL)) THEN
            DZ_TMP1 = ABS(DZ_FINAL)-ABS(DZ_NODE)
            IF ( DZ_TMP1 < DZ_TMP2 ) THEN
               FINALTARGETFOUND = .TRUE.
               HP%I_TARGET = -EVAC_NODE_LIST(SP%NODES_OUT(I))%NODE_INDEX
               DZ_TMP2 = DZ_TMP1
               DZ_NODE_ACTUAL = DZ_NODE
            END IF
         END IF
      END DO FINDTARGETNODELOOP

      ! Find the target among the nodes leading in to the strs, i.e. doors
      IF (.NOT.FINALTARGETFOUND) THEN
         DIST_TO_DOOR = HUGE(DIST_TO_DOOR)
         FINDTARGETNODELOOP2: DO I = 1,SP%N_NODES_IN
            INODE = SP%NODES_IN(I)
            DZ_NODE = -1._EB * SIGN(1.0_EB,DZ_FINAL) ! Initialize dz_node to different direction than final target
            SELECT CASE(EVAC_NODE_LIST(INODE)%NODE_TYPE)
            CASE('Door')
               Z_NODE = 0.5_EB*(EVAC_DOORS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%Z1+EVAC_DOORS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%Z2)
               X_NODE = 0.5_EB*(EVAC_DOORS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%X1+EVAC_DOORS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%X2)
               Y_NODE = 0.5_EB*(EVAC_DOORS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%Y1+EVAC_DOORS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%Y2)
               STR_SUB_INDX = EVAC_DOORS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%STR_SUB_INDX
            CASE('Exit')
               !            Z_NODE = EVAC_EXITS(EVAC_NODE_LIST(SP%NODES_IN(I))%NODE_INDEX)%Z1
            CASE('Entry')
               !            Z_NODE = EVAC_ENTRYS(EVAC_NODE_LIST(SP%NODES_IN(I))%NODE_INDEX)%Z1
            CASE DEFAULT
               WRITE(LU_ERR,*) 'ERROR (DEBUG): Unknown node type in find_target_node_in_strs'
            END SELECT
            DZ_NODE = Z_NODE - HP%Z
            IF (SIGN(1.0_EB,DZ_NODE)==SIGN(1.0_EB,DZ_FINAL)) THEN
               DZ_TMP1 = ABS(DZ_FINAL-DZ_NODE)
               IF (HP%STR_SUB_INDX==STR_SUB_INDX) THEN
                  DIST_TO_DOOR_TMP = (X_NODE-HP%X)**2 + (Y_NODE-HP%Y)**2
               ELSE
                  DIST_TO_DOOR_TMP = -1.0_EB
               END IF
               IF ( DZ_TMP1 < DZ_TMP2) THEN
                     FINALTARGETFOUND = .TRUE.
                     HP%I_TARGET = -EVAC_NODE_LIST(INODE)%NODE_INDEX
                     DZ_TMP2 = DZ_TMP1
                     DZ_NODE_ACTUAL = DZ_NODE
               END IF
            END IF
         END DO FINDTARGETNODELOOP2
      END IF
      IF (DZ_NODE_ACTUAL >= 0._EB) THEN
         HP%STRS_DIRECTION = 1
      ELSE
         HP%STRS_DIRECTION = -1
      END IF
    END SUBROUTINE FIND_TARGET_NODE_IN_STRS

    SUBROUTINE STRS_U_AND_V(STRP,I_NODE,X,Y,DIRECTION,UBAR,VBAR)
      IMPLICIT NONE
      !
      ! Get preferred direction in STRS
      !
      ! Passed variables
      TYPE (EVAC_STRS_TYPE), POINTER::  STRP
      INTEGER, INTENT(IN) :: I_NODE, DIRECTION
      REAL(EB), INTENT(IN) :: X, Y
      REAL(EB), INTENT(OUT) :: UBAR, VBAR
      !
      ! Local variables
      INTEGER I_CORE

      I_CORE = STRP%I_CORE(I_NODE)

      UBAR = 0._EB
      VBAR = 0._EB

      IF (STRP%RIGHT_HANDED) THEN
         IF (DIRECTION > 0) THEN
            IF ((X <= STRP%XB_CORE(I_CORE,2)) .AND. (Y <= STRP%XB_CORE(I_CORE,3))) UBAR =  1._EB
            IF ((X >= STRP%XB_CORE(I_CORE,2)) .AND. (Y <= STRP%XB_CORE(I_CORE,4))) VBAR =  1._EB
            IF ((X >= STRP%XB_CORE(I_CORE,1)) .AND. (Y >= STRP%XB_CORE(I_CORE,4))) UBAR = -1._EB
            IF ((X <= STRP%XB_CORE(I_CORE,1)) .AND. (Y >= STRP%XB_CORE(I_CORE,3))) VBAR = -1._EB
         ELSE
            IF ((X <= STRP%XB_CORE(I_CORE,1)) .AND. (Y <= STRP%XB_CORE(I_CORE,4))) VBAR =  1._EB
            IF ((X <= STRP%XB_CORE(I_CORE,2)) .AND. (Y >= STRP%XB_CORE(I_CORE,4))) UBAR =  1._EB
            IF ((X >= STRP%XB_CORE(I_CORE,2)) .AND. (Y >= STRP%XB_CORE(I_CORE,3))) VBAR = -1._EB
            IF ((X >= STRP%XB_CORE(I_CORE,1)) .AND. (Y <= STRP%XB_CORE(I_CORE,3))) UBAR = -1._EB
         ENDIF
      ELSE
         IF (DIRECTION < 0) THEN
            IF ((X <= STRP%XB_CORE(I_CORE,2)) .AND. (Y <= STRP%XB_CORE(I_CORE,3))) UBAR =  1._EB
            IF ((X >= STRP%XB_CORE(I_CORE,2)) .AND. (Y <= STRP%XB_CORE(I_CORE,4))) VBAR =  1._EB
            IF ((X >= STRP%XB_CORE(I_CORE,1)) .AND. (Y >= STRP%XB_CORE(I_CORE,4))) UBAR = -1._EB
            IF ((X <= STRP%XB_CORE(I_CORE,1)) .AND. (Y >= STRP%XB_CORE(I_CORE,3))) VBAR = -1._EB
         ELSE
            IF ((X <= STRP%XB_CORE(I_CORE,1)) .AND. (Y <= STRP%XB_CORE(I_CORE,4))) VBAR =  1._EB
            IF ((X <= STRP%XB_CORE(I_CORE,2)) .AND. (Y >= STRP%XB_CORE(I_CORE,4))) UBAR =  1._EB
            IF ((X >= STRP%XB_CORE(I_CORE,2)) .AND. (Y >= STRP%XB_CORE(I_CORE,3))) VBAR = -1._EB
            IF ((X >= STRP%XB_CORE(I_CORE,1)) .AND. (Y <= STRP%XB_CORE(I_CORE,3))) UBAR = -1._EB
         ENDIF
      ENDIF
    END SUBROUTINE STRS_U_AND_V

    SUBROUTINE GET_IW(IIIN,JJIN,KKIN,IOR,IW)
      IMPLICIT NONE
      !
      ! Passed variables
      INTEGER, INTENT(IN) :: IIIN, JJIN, KKIN, IOR
      INTEGER, INTENT(OUT) :: IW
      !
      ! Local variables
      INTEGER :: II, JJ, KK, IC, I_OBST
      !
      II = IIIN
      JJ = JJIN
      KK = KKIN

      IC  = CELL_INDEX(II,JJ,KK)
      I_OBST = OBST_INDEX_C(IC)

      IF (SOLID(IC) .AND. .NOT.OBSTRUCTION(I_OBST)%HIDDEN) THEN
         SELECT CASE(IOR)
         CASE(-1)
            IF (II>0)      II = II-1
         CASE( 1)
            IF (II<IBAR+1) II = II+1
         CASE(-2)
            IF (JJ>0)      JJ = JJ-1
         CASE( 2)
            IF (JJ<JBAR+1) JJ = JJ+1
         CASE(-3)
            IF (KK>0)      KK = KK-1
         CASE( 3)
            IF (KK<KBAR+1) KK = KK+1
         END SELECT
      END IF

      IC  = CELL_INDEX(II,JJ,KK)
      IW  = WALL_INDEX(IC,-IOR)

      IF (IW<=0) THEN
         SELECT CASE(IOR)
         CASE(-1)
            IF (II>0)      IC = CELL_INDEX(II-1,JJ,KK)
         CASE( 1)
            IF (II<IBAR+1) IC = CELL_INDEX(II+1,JJ,KK)
         CASE(-2)
            IF (JJ>0)      IC = CELL_INDEX(II,JJ-1,KK)
         CASE( 2)
            IF (JJ<JBAR+1) IC = CELL_INDEX(II,JJ+1,KK)
         CASE(-3)
            IF (KK>0)      IC = CELL_INDEX(II,JJ,KK-1)
         CASE( 3)
            IF (KK<KBAR+1) IC = CELL_INDEX(II,JJ,KK+1)
         END SELECT
         IW = WALL_INDEX(IC,-IOR)
      END IF

    END SUBROUTINE GET_IW

    !
    SUBROUTINE CHECK_EXITS(T,NM)
      IMPLICIT NONE
      !
      ! Remove persons if they are found at an exit or just count if COUNT_ONLY=T.
      !
      ! Passed variables
      REAL(EB), INTENT(IN) :: T
      INTEGER, INTENT(IN) :: NM
      !
      ! Local variables
      REAL(EB) X_OLD, Y_OLD, PEXX1, PEXX2, PEXY1, PEXY2
      INTEGER :: IE,I,N_TMP
      LOGICAL :: L_INIT_IOR, CLOSED
      TYPE (EVAC_EXIT_TYPE), POINTER :: PEX =>NULL()
      TYPE (HUMAN_TYPE), POINTER :: HR =>NULL()
      !
      L_INIT_IOR = .TRUE.
      ! HUMAN(:)%IOR = 0
      PEXLOOP: DO IE = 1, N_EXITS
         PEX=>EVAC_EXITS(IE)
         IF (PEX%IMESH /= NM ) CYCLE PEXLOOP
         IF(PEX%LOCKED_WHEN_CLOSED .AND. ABS(PEX%IMODE) == 2) THEN
            CLOSED = .TRUE.
         ELSE
            CLOSED = .FALSE.
         END IF
         IF (L_INIT_IOR) THEN
            HUMAN(:)%IOR = 0 ; L_INIT_IOR = .FALSE.
         END IF
         SELECT CASE (PEX%IOR)
         CASE (-1,+1)
            PEXX1 = PEX%X1
            PEXX2 = PEX%X2
            PEXY1 = 0.5_EB*(PEX%Y1+PEX%Y2)-0.5_EB*PEX%WIDTH
            PEXY2 = 0.5_EB*(PEX%Y1+PEX%Y2)+0.5_EB*PEX%WIDTH
         CASE (-2,+2)
            PEXX1 = 0.5_EB*(PEX%X1+PEX%X2)-0.5_EB*PEX%WIDTH
            PEXX2 = 0.5_EB*(PEX%X1+PEX%X2)+0.5_EB*PEX%WIDTH
            PEXY1 = PEX%Y1
            PEXY2 = PEX%Y2
         CASE (-3,+3)  ! Count persons on a horizontal plane given by XB
            PEXX1 = PEX%X1
            PEXX2 = PEX%X2
            PEXY1 = PEX%Y1
            PEXY2 = PEX%Y2
            PEX%ICOUNT = 0
         END SELECT
         HUMLOOP: DO I = 1, N_HUMANS
            HR=>HUMAN(I)

            ! Check only one true, i.e., not a count only, exit.
            IF ( HR%IOR > 0 ) CYCLE HUMLOOP

            ! Counters can distinguish PERS and EVAC/ENTR namelists.
            IF (PEX%COUNT_ONLY .AND. Trim(PEX%PERS_ID) /= 'null') THEN
               IF (Trim(PEX%PERS_ID) /= Trim(EVAC_PERSON_CLASSES(HR%IPC)%ID)) Cycle HUMLOOP
            END IF
            IF (PEX%COUNT_ONLY .AND. Trim(PEX%EVAC_ID) /= 'null') THEN
               IF (HR%IEL > 0) THEN
                  IF (Trim(PEX%EVAC_ID) /= Trim(EVAC_EVACS(HR%IEL)%ID)) Cycle HUMLOOP
               ELSE
                  IF (Trim(PEX%EVAC_ID) /= Trim(EVAC_ENTRYS(ABS(HR%IEL))%ID)) Cycle HUMLOOP
               END IF
            END IF

            X_OLD = HR%X_OLD
            Y_OLD = HR%Y_OLD
            SELECT CASE (PEX%IOR)
            CASE (+1)
               IF ((HR%X >= PEX%X1 .AND. X_OLD < PEX%X1) .AND. (HR%Y >= PEX%Y1 .AND. HR%Y <= PEX%Y2) ) THEN
                  IF (CLOSED) THEN
                     HR%IOR = 0  ! Do not remove the person
                     HR%X = HR%X_OLD
                     HR%Y = HR%Y_OLD
                     HR%U = 0.0_EB
                     HR%V = 0.0_EB
                  ELSE
                     IF (PEX%COUNT_ONLY) HR%IOR = 2
                     IF (.NOT. PEX%COUNT_ONLY) HR%IOR = 1
                     PEX%T_LAST=T
                     PEX%ICOUNT = PEX%ICOUNT + 1
                     IF (PEX%T_FIRST <= T_BEGIN) PEX%T_FIRST = T
                  END IF
               END IF
            CASE (-1)
               IF ((HR%X <= PEX%X2 .AND. X_OLD > PEX%X2) .AND. (HR%Y >= PEX%Y1 .AND. HR%Y <= PEX%Y2) ) THEN
                  IF (CLOSED) THEN
                     HR%IOR = 0  ! Do not remove the person
                     HR%X = HR%X_OLD
                     HR%Y = HR%Y_OLD
                     HR%U = 0.0_EB
                     HR%V = 0.0_EB
                  ELSE
                     IF (PEX%COUNT_ONLY) HR%IOR = 2
                     IF (.NOT. PEX%COUNT_ONLY) HR%IOR = 1
                     PEX%T_LAST=T
                     PEX%ICOUNT = PEX%ICOUNT + 1
                     IF (PEX%T_FIRST <= T_BEGIN) PEX%T_FIRST = T
                  END IF
               END IF
            CASE (+2)
               IF ((HR%Y >= PEX%Y1 .AND. Y_OLD < PEX%Y1) .AND. (HR%X >= PEX%X1 .AND. HR%X <= PEX%X2) ) THEN
                  IF (CLOSED) THEN
                     HR%IOR = 0  ! Do not remove the person
                     HR%X = HR%X_OLD
                     HR%Y = HR%Y_OLD
                     HR%U = 0.0_EB
                     HR%V = 0.0_EB
                  ELSE
                     IF (PEX%COUNT_ONLY) HR%IOR = 2
                     IF (.NOT. PEX%COUNT_ONLY) HR%IOR = 1
                     PEX%T_LAST=T
                     PEX%ICOUNT = PEX%ICOUNT + 1
                     IF (PEX%T_FIRST <= T_BEGIN) PEX%T_FIRST = T
                  END IF
               END IF
            CASE (-2)
               IF ((HR%Y <= PEX%Y2 .AND. Y_OLD > PEX%Y2) .AND. (HR%X >= PEX%X1 .AND. HR%X <= PEX%X2) ) THEN
                  IF (CLOSED) THEN
                     HR%IOR = 0  ! Do not remove the person
                     HR%X = HR%X_OLD
                     HR%Y = HR%Y_OLD
                     HR%U = 0.0_EB
                     HR%V = 0.0_EB
                  ELSE
                     IF (PEX%COUNT_ONLY) HR%IOR = 2
                     IF (.NOT. PEX%COUNT_ONLY) HR%IOR = 1
                     PEX%T_LAST=T
                     PEX%ICOUNT = PEX%ICOUNT + 1
                     IF (PEX%T_FIRST <= T_BEGIN) PEX%T_FIRST = T
                  END IF
               END IF
            CASE (-3,+3)
               ! These are always count only exits, i.e., just counters.
               IF (HR%Y >= PEX%Y1 .AND. HR%Y <= PEX%Y2 .AND. HR%X >= PEX%X1 .AND. HR%X <= PEX%X2) THEN
                  HR%IOR = 2
                  PEX%T_LAST=T
                  PEX%ICOUNT = PEX%ICOUNT + 1
                  IF (PEX%T_FIRST <= T_BEGIN) PEX%T_FIRST = T
               END IF
            END SELECT
            IF (HR%IOR > 0) THEN
               IF (C_HAWK >= 0.0_EB) THEN
                  IF (.NOT. PEX%COUNT_ONLY) WRITE (LU_EVACOUT,FMT='(A,I6,A,F8.2,A,A,A,F8.1,A,3I3)') &
                       ' Agent n:o', HR%ILABEL, ' out at ', T, ' s, exit ', TRIM(PEX%ID), &
                       ', A=', HR%A, ', AgentType=', HR%I_DoorAlgo, HR%COLOR_INDEX, HR%I_Target
                  IF (PEX%COUNT_ONLY .AND. ABS(PEX%IOR)<3) WRITE (LU_EVACOUT,FMT='(A,I6,A,F8.2,A,A,A,F8.1,A,3I3)') &
                       ' Agent n:o', HR%ILABEL, ' counted at ', T, ' s, counter ', TRIM(PEX%ID), &
                       ', A=', HR%A, ', AgentType=', HR%I_DoorAlgo, HR%COLOR_INDEX, HR%I_Target
                  IF (HR%IOR == 2) HR%IOR = HUMAN_NO_TARGET
               ELSE
                  IF (.NOT. PEX%COUNT_ONLY) WRITE (LU_EVACOUT,FMT='(A,I6,A,F8.2,A,A,A,F8.4,A,I4,2I3)') &
                       ' Agent n:o', HR%ILABEL, ' out at ', T, ' s, exit ', TRIM(PEX%ID), &
                       ', FED=', HR%INTDOSE, ', Color_i=', HR%COLOR_INDEX, HR%I_Target, HR%I_DoorAlgo
                  IF (PEX%COUNT_ONLY .AND. ABS(PEX%IOR)<3) WRITE (LU_EVACOUT,FMT='(A,I6,A,F8.2,A,A,A,F8.4,A,I4,2I3)') &
                       ' Agent n:o', HR%ILABEL, ' counted at ', T, ' s, counter ', TRIM(PEX%ID), &
                       ', FED=', HR%INTDOSE, ', Color_i=', HR%COLOR_INDEX, HR%I_Target, HR%I_DoorAlgo
                  IF (HR%IOR == 2) HR%IOR = HUMAN_NO_TARGET
               END IF
            END IF
         END DO HUMLOOP
      END DO PEXLOOP
      N_TMP = N_HUMANS
      REMOVE_LOOP: DO I = N_TMP, 1, -1
         HR=>HUMAN(I)
         ! Remove person, if wanted.
         IF (HR%IOR >= 1) THEN
            HR%IOR = HUMAN_NO_TARGET
            CALL REMOVE_PERSON(I)
         END IF
      END DO REMOVE_LOOP
      !
    END SUBROUTINE CHECK_EXITS
    !
    SUBROUTINE CHECK_DOORS(T,NM)
      IMPLICIT NONE
      !
      ! Replace persons if they are found at a door.
      !
      ! Passed variables
      REAL(EB), INTENT(IN) :: T
      INTEGER, INTENT(IN) :: NM
      !
      ! Local variables
      REAL(EB) X_OLD, Y_OLD, XX, YY, ZZ, PDXX1, PDXX2, PDXY1, PDXY2, V, ANGLE
      INTEGER :: IE,I,N_TMP, ISTAT, IOR_NEW, INODE2, IMESH2, N, IOR
      INTEGER :: NEW_FFIELD_I, COLOR_INDEX, I_TARGET, INODE, STR_INDX, STR_SUB_INDX
      CHARACTER(LABEL_LENGTH) :: NEW_FFIELD_NAME
      LOGICAL :: KEEP_XY, UPSTREAM, NO_TO_NODE, L_INIT_IOR, CLOSED
      TYPE (EVAC_DOOR_TYPE), POINTER :: PDX =>NULL()
      TYPE (HUMAN_TYPE), POINTER :: HR =>NULL()
      !
      L_INIT_IOR = .TRUE.
      PDXLOOP: DO IE = 1, N_DOORS
         NO_TO_NODE = .FALSE.
         KEEP_XY = .FALSE.
         PDX=>EVAC_DOORS(IE)
         IF (Trim(PDX%TO_NODE)=='null') NO_TO_NODE = .TRUE.
         ! Note: IMESH2 is not good for corr targets
         IF (NO_TO_NODE) CYCLE PDXLOOP
         IF (PDX%IMESH /= NM .AND. PDX%IMESH2 /= NM) CYCLE PDXLOOP
         IF (PDX%IMESH /= NM .AND. .NOT.NM_STRS_MESH) CYCLE PDXLOOP
         IF (TRIM(PDX%TO_NODE) == TRIM(PDX%ID)) CYCLE PDXLOOP ! An entry mode for a door.
         IF(PDX%LOCKED_WHEN_CLOSED .AND. ABS(PDX%IMODE) == 2) THEN
            CLOSED = .TRUE.
         ELSE
            CLOSED = .FALSE.
         END IF
         IF (L_INIT_IOR) THEN
            HUMAN(:)%IOR = 0 ; L_INIT_IOR = .FALSE.
         END IF
         KEEP_XY = PDX%KEEP_XY
         SELECT CASE (PDX%IOR)
         CASE (-1,+1)
            PDXX1 = PDX%X1
            PDXX2 = PDX%X2
            PDXY1 = 0.5_EB*(PDX%Y1+PDX%Y2)-0.5_EB*PDX%WIDTH
            PDXY2 = 0.5_EB*(PDX%Y1+PDX%Y2)+0.5_EB*PDX%WIDTH
         CASE (-2,+2)
            PDXX1 = 0.5_EB*(PDX%X1+PDX%X2)-0.5_EB*PDX%WIDTH
            PDXX2 = 0.5_EB*(PDX%X1+PDX%X2)+0.5_EB*PDX%WIDTH
            PDXY1 = PDX%Y1
            PDXY2 = PDX%Y2
         END SELECT
         HUMLOOP: DO I = 1, N_HUMANS
            HR=>HUMAN(I)
            NEW_FFIELD_NAME = TRIM(HR%FFIELD_NAME)
            NEW_FFIELD_I = HR%I_FFIELD
            ! Check only one door... (Yes, this can happen...)
            IF ( HR%IOR /= HUMAN_NO_TARGET) CYCLE HUMLOOP
            ! Coming out upstream?
            UPSTREAM = .FALSE.
            IF (PDX%IMESH2 == NM .AND. NM_STRS_MESH) THEN
               ! Here, check for correct height
               IF ((HR%Z > PDX%Z2) .OR. (HR%Z < PDX%Z1)) CYCLE HUMLOOP
               UPSTREAM = .TRUE.
            END IF
            X_OLD = HR%X_OLD
            Y_OLD = HR%Y_OLD
            SELECT CASE (ABS(PDX%IOR))
            CASE (1)
               IF (SIGN(1.0_EB,HR%X-PDX%X1)/=SIGN(1.0_EB,X_OLD-PDX%X1) .AND. (HR%Y >= PDX%Y1 .AND. HR%Y <= PDX%Y2) ) THEN
                  HR%IOR = HUMAN_TARGET_UNSPECIFIED
               END IF
            CASE (2)
               IF (SIGN(1.0_EB,HR%Y-PDX%Y1)/=SIGN(1.0_EB,Y_OLD-PDX%Y1) .AND. (HR%X >= PDX%X1 .AND. HR%X <= PDX%X2) ) THEN
                  HR%IOR = HUMAN_TARGET_UNSPECIFIED
               END IF
            END SELECT
            IF ( HR%IOR == HUMAN_TARGET_UNSPECIFIED ) THEN
               ISTAT = 0
               IF (UPSTREAM) THEN
                  INODE  = PDX%INODE2
                  INODE2 = PDX%INODE
               ELSE
                  INODE  = PDX%INODE
                  INODE2 = PDX%INODE2
               END IF
               IF (CLOSED) THEN
                  ISTAT = 1
               ELSE
                  CALL CHECK_TARGET_NODE(INODE,INODE2,ISTAT,XX,YY,ZZ,IOR_NEW, IMESH2,T,NEW_FFIELD_NAME,NEW_FFIELD_I, &
                       COLOR_INDEX, I_TARGET, KEEP_XY, ANGLE, STR_INDX, STR_SUB_INDX, HR)
               END IF
               IF (ISTAT == 0 ) THEN
                  ! Put person to a new node, i.e., target node has empty space
                  HR%X = XX
                  HR%Y = YY
                  HR%Z = ZZ
                  HR%ANGLE = ANGLE
                  HR%I_Target2 = 0
                  IF (STR_INDX>0) HR%STR_SUB_INDX = STR_SUB_INDX
                  IF (KEEP_XY .AND. IOR_NEW == HUMAN_TARGET_UNSPECIFIED) THEN
                     IF (EVAC_NODE_LIST(INODE2)%NODE_TYPE == 'Door') THEN
                        XX = 0.5_EB*(EVAC_DOORS(EVAC_NODE_LIST(INODE2)%NODE_INDEX)%X1 + &
                             EVAC_DOORS(EVAC_NODE_LIST(INODE2)%NODE_INDEX)%X2 - (PDX%X1+PDX%X2))
                        YY = 0.5_EB*(EVAC_DOORS(EVAC_NODE_LIST(INODE2)%NODE_INDEX)%Y1 + &
                             EVAC_DOORS(EVAC_NODE_LIST(INODE2)%NODE_INDEX)%Y2 - (PDX%Y1+PDX%Y2))
                     ELSE
                        XX = 0.5_EB*(EVAC_ENTRYS(EVAC_NODE_LIST(INODE2)%NODE_INDEX)%X1 + &
                             EVAC_ENTRYS(EVAC_NODE_LIST(INODE2)%NODE_INDEX)%X2 - (PDX%X1+PDX%X2))
                        YY = 0.5_EB*(EVAC_ENTRYS(EVAC_NODE_LIST(INODE2)%NODE_INDEX)%Y1 + &
                             EVAC_ENTRYS(EVAC_NODE_LIST(INODE2)%NODE_INDEX)%Y2 - (PDX%Y1+PDX%Y2))
                     END IF
                     IF (ABS(PDX%IOR) == 2) HR%X = HR%X + XX
                     IF (ABS(PDX%IOR) == 1) HR%Y = HR%Y + YY
                  END IF
                  V = SQRT( HR%U**2 + HR%V**2 )

                  HR%IOR = IOR_NEW
                  IF (IOR_NEW == HUMAN_TARGET_UNSPECIFIED .OR. IOR_NEW == HUMAN_STRS_TARGET) THEN
                     HR%I_Door_Mode = 0  ! Default, main evacuation flow field agent
                     IF (I_Target/=0 .AND. HR%I_DoorAlgo>0) HR%I_Door_Mode = 1  ! Found a target door/exit
                     IF (I_Target==0 .AND. HR%I_DoorAlgo>0) HR%I_Door_Mode = -INODE2
                     IF (EVAC_NODE_LIST(INODE2)%NODE_TYPE == 'Door') THEN
                        IOR = - EVAC_DOORS(EVAC_NODE_LIST(INODE2)%NODE_INDEX)%IOR
                     END IF
                     IF (EVAC_NODE_LIST(INODE2)%NODE_TYPE == 'Entry') THEN
                        IOR = EVAC_ENTRYS(EVAC_NODE_LIST(INODE2)%NODE_INDEX)%IOR
                     END IF
                     IF (EVAC_NODE_LIST(INODE2)%NODE_TYPE == 'Floor') THEN
                        ! For STRS the first node should be DOOR
                        IOR = EVAC_DOORS(EVAC_NODE_LIST(INODE)%NODE_INDEX)%IOR
                     END IF

                     ! DOOR/ENTRY or STRS-mesh
                     IF (HR%IMESH /= IMESH2 ) THEN
                        ! Put the person to a new floor
                        HR%IOR         = HUMAN_ANOTHER_MESH_TARGET
                        HR%INODE = 0
                        DO N = 1, N_EGRIDS
                           IF (EVAC_NODE_LIST(N)%IMESH == IMESH2) HR%INODE = N
                        END DO
                        HR%NODE_NAME   = TRIM(MESH_NAME(IMESH2))
                        HR%FFIELD_NAME = TRIM(NEW_FFIELD_NAME)
                        HR%I_FFIELD = NEW_FFIELD_I
                        HR%SKIP_WALL_FORCE_IOR = 0
                     ELSE
                        HR%IOR = -2   ! Same floor (DOOR/ENTRY)
                        HR%FFIELD_NAME = TRIM(NEW_FFIELD_NAME)
                        HR%I_FFIELD = NEW_FFIELD_I
                     END IF
                     HR%I_TARGET = I_TARGET
                     IF (COLOR_METHOD == 4) THEN
                        IF (ABS(I_TARGET)>0 .AND. ABS(I_TARGET)<=N_DOORS) &
                             HR%COLOR_INDEX = EVAC_DOORS(ABS(I_TARGET))%Avatar_Color_Index
                        IF (ABS(I_TARGET)>N_DOORS .AND. ABS(I_TARGET)<=N_DOORS+N_EXITS) &
                             HR%COLOR_INDEX = EVAC_EXITS(ABS(I_TARGET)-N_DOORS)%Avatar_Color_Index
                        IF (I_TARGET==0) HR%COLOR_INDEX = EVAC_AVATAR_NCOLOR ! default, cyan
                     END IF
                     IF (COLOR_METHOD==5 .AND. (HR%I_DoorAlgo==3 .OR. &
                          (I_Target==0 .AND. HR%I_DoorAlgo>0))) HR%COLOR_INDEX = EVAC_AVATAR_NCOLOR ! default, cyan
                     HR%X_OLD = HR%X
                     HR%Y_OLD = HR%Y
                     ! IOR is the direction where the human is ejected.
                     IF ( ABS(IOR) == 1 ) THEN
                        HR%U = V*IOR
                        HR%V = 0.0_EB
                        IF (.NOT.KEEP_XY) THEN
                           ! 180 or 0 degrees, i.e., pi or 0 radians
                           HR%ANGLE = (0.5_EB-(IOR/2.0_EB))*PI
                        END IF
                     END IF
                     IF ( ABS(IOR) == 2 ) THEN
                        HR%U = 0.0_EB
                        HR%V = 0.5_EB*V*IOR
                        IF (.NOT.KEEP_XY) THEN
                           ! 270 or 90 degrees, i.e., 3pi/2 or pi/2 radians
                           HR%ANGLE = (1.0_EB-(IOR/4.0_EB))*PI
                        END IF
                     END IF

                     IF (HR%IMESH /= IMESH2 ) THEN
                        HR%IMESH       = IMESH2
                        IF (MESHES(IMESH2)%N_HUMANS+1 > MESHES(IMESH2)%N_HUMANS_DIM) THEN
                           ! Re-allocation is not yet checked.
                           WRITE(LU_ERR,*) 'ERROR: HUMANS: NO RE-ALLOCATION YET'
                           WRITE(LU_EVACOUT,*) 'ERROR: HUMANS: NO RE-ALLOCATION YET'
                           STOP_STATUS = EVACUATION_STOP ; RETURN
                        END IF
                        MESHES(IMESH2)%N_HUMANS = MESHES(IMESH2)%N_HUMANS + 1
                        MESHES(IMESH2)%HUMAN(MESHES(IMESH2)%N_HUMANS) = HUMAN(I)
                        MESHES(IMESH2)%HUMAN(MESHES(IMESH2)%N_HUMANS)%IOR = 0
                     END IF

                  END IF            ! Target is DOOR or ENTRY

                  PDX%T_LAST=T
                  PDX%ICOUNT = PDX%ICOUNT + 1
                  IF (PDX%T_FIRST <= T_BEGIN) PDX%T_FIRST = T

                  IF (C_HAWK >= 0.0_EB) THEN
                     WRITE (LU_EVACOUT,FMT='(A,I6,A,F8.2,A,A,A,F8.1,A,3I3)') &
                          ' Agent n:o', HR%ILABEL, ' out at ', T, ' s, door ', TRIM(PDX%ID), &
                          ' A=', HR%A, ', AgentType=', HR%I_DoorAlgo, HR%COLOR_INDEX, HR%I_Target
                  ELSE
                     WRITE (LU_EVACOUT,FMT='(A,I6,A,F8.2,A,A,A,F8.4,A,I4,2I3)') &
                          ' Agent n:o', HR%ILABEL, ' out at ', T, ' s, door ', TRIM(PDX%ID), &
                          ' FED ', HR%INTDOSE, ', Color_i=', HR%COLOR_INDEX, HR%I_Target, HR%I_DoorAlgo
                  END IF

               ELSE    ! ISTAT = 1 ==> DO NOT MOVE TO NODE
                  ! Can not move to the next node, so do not allow to move inside
                  ! the DOOR ==> keep old position and put velocity to zero.
                  HR%X = HR%X_OLD
                  HR%Y = HR%Y_OLD
                  HR%U = 0.0_EB
                  HR%V = 0.0_EB
                  HR%IOR = -1  ! Can not move to target node
               END IF
            END IF
         END DO HUMLOOP
      END DO PDXLOOP

      ! Remove humans from this floor and put them some other floor
      N_TMP = N_HUMANS
      REMOVE_LOOP: DO I = N_TMP, 1, -1
         HR=>HUMAN(I)
         ! HR%IOR = 2: Put to an another floor (target is DOOR/ENTRY)
         !          0: Not entering a DOOR
         !          3: Target is corridor (remove from floor)
         !          4: Not used (floor node...)
         !          5: Target is EXIT (Remove from floor)
         !         -1: Can not move to the target node
         !         -2: Move to DOOR/ENTRY on the same floor
         !
         ! Remove person from the floor
         IF (HR%IOR >= 1) THEN
            HR%IOR = HUMAN_NO_TARGET
            CALL REMOVE_PERSON(I)
         END IF
      END DO REMOVE_LOOP
      !
    END SUBROUTINE CHECK_DOORS
    !
    SUBROUTINE CHECK_CORRS(T,NM,DTSP)
      IMPLICIT NONE
      !
      ! Entry persons from the corridors.
      !
      ! Passed variables
      REAL(EB), INTENT(IN) :: T,DTSP
      INTEGER, INTENT(IN) :: NM
      !
      ! Local variables
      REAL(EB) XX, YY, ZZ, V, X_INT, ANGLE
      INTEGER :: IE, ISTAT, IOR_NEW, INODE2, IMESH2, N, IOR
      INTEGER :: NEW_FFIELD_I, COLOR_INDEX, I_TARGET, INODE, STR_INDX, STR_SUB_INDX
      CHARACTER(LABEL_LENGTH) :: NEW_FFIELD_NAME
      LOGICAL :: KEEP_XY
      TYPE (EVAC_CORR_TYPE),  POINTER :: PCX=>NULL()
      TYPE (CORR_LL_TYPE), POINTER :: NOW_LL=>NULL(), TMP_LL=>NULL(), NEXT_LL=>NULL(), PREV_LL=>NULL()
      TYPE (HUMAN_TYPE), POINTER :: HR =>NULL()
      !
      KEEP_XY = .FALSE.
      PCXLOOP: DO IE = 1, N_CORRS
         PCX=>EVAC_CORRS(IE)
         IF (PCX%IMESH2 /= NM) CYCLE PCXLOOP
         Elevator_If: IF (PCX%ELEVATOR) THEN
            IF (PCX%Elevator_Mode == 2 .AND. PCX%Time_Start < T .AND. PCX%N_INSIDE > 0)  THEN
               ! The car is at the upstairs, there is someone in the car, it is time to go down
               PCX%Time_Stop = PCX%Time_Start + PCX%Travel_Time
               PCX%Elevator_Mode = 0  ! The car is moving down
               IF (PCX%IDOOR_FROM > 0) THEN
                  ! Close the door if it is not yet closed.
                  IF (EVAC_DOORS(PCX%IDOOR_FROM)%IMODE > -2) EVAC_DOORS(PCX%IDOOR_FROM)%IMODE = +2
               END IF
               NULLIFY(PREV_LL)
               NULLIFY(NOW_LL)
               NULLIFY(NEXT_LL)
               IF (ASSOCIATED(PCX%FIRST)) THEN
                  NOW_LL => PCX%FIRST
                  IF (ASSOCIATED(PCX%FIRST%NEXT)) NEXT_LL => PCX%FIRST%NEXT
               END IF
               HUMLOOP_ELE: DO
                  IF (.NOT. ASSOCIATED(NOW_LL)) EXIT HUMLOOP_ELE
                  NOW_LL%T_out = PCX%Time_Stop ! The time when the car is at downstairs
                  PREV_LL => NOW_LL  ! Save the previous element
                  IF (ASSOCIATED(NOW_LL%NEXT)) THEN
                     NOW_LL => NOW_LL%NEXT ! Jump to the next element in LL
                  ELSE
                     NULLIFY(NOW_LL)  ! End the loop
                  END IF
               END DO HUMLOOP_ELE
            END IF
            IF (PCX%Elevator_Mode == 0 .AND. PCX%Time_Stop < T) THEN
               ! The car is moving down
               PCX%Elevator_Mode = 1  ! The car is at downstairs
            END IF
            IF (PCX%Elevator_Mode == 1 .AND. PCX%N_INSIDE <= 0) THEN
               ! The car is at the downstairs and it is already empty
               PCX%Elevator_Mode = 3  ! The car is moving up
               PCX%Time_Start = T
               PCX%Time_Stop  = PCX%Time_Start + PCX%Travel_Time
            END IF
            IF (PCX%Elevator_Mode == 3 .AND. PCX%Time_Stop < T) THEN
               ! The car is moving up and coming to the upper floor
               PCX%Elevator_Mode = 2  ! The car is at upstairs
               IF (PCX%IDOOR_FROM > 0) THEN
                  ! Open the door if it is not yet open.
                  IF (EVAC_DOORS(PCX%IDOOR_FROM)%IMODE /= -1) EVAC_DOORS(PCX%IDOOR_FROM)%IMODE = +1
               END IF
            END IF
         END IF Elevator_If

         IF (PCX%N_INSIDE <= 0) CYCLE PCXLOOP
         NULLIFY(PREV_LL)
         NULLIFY(NOW_LL)
         NULLIFY(NEXT_LL)
         IF (ASSOCIATED(PCX%FIRST)) THEN
            NOW_LL => PCX%FIRST
            IF (ASSOCIATED(PCX%FIRST%NEXT)) NEXT_LL => PCX%FIRST%NEXT
         END IF
         HUMLOOP: DO
            IF (.NOT. ASSOCIATED(NOW_LL)) EXIT HUMLOOP
            ISTAT = 0
            INODE = PCX%INODE
            INODE2 = PCX%INODE2
            HR => NOW_LL%HUMAN
            IF ( HR%INTDOSE >= 1.0_EB  ) THEN
               IF (HR%TPRE /= HUGE(HR%TPRE)) THEN
                  N_DEAD = N_DEAD+1
                  HR%TPRE = HUGE(HR%TPRE)
                  HR%TDET = HUGE(HR%TDET)
                  HR%TAU  = HR%TAU
                  HR%MASS = HR%MASS
               END IF
            ELSE
               FED_MAX_ALIVE = MAX(FED_MAX_ALIVE,HR%INTDOSE)
            END IF
            FED_MAX = MAX(FED_MAX,HR%INTDOSE)

            ! calculate Purser's fractional effective dose (FED)
            IF (T > T_BEGIN) THEN
               IF ( PCX%FED_MESH2 > 0 ) THEN
                  X_INT = MIN(1.0_EB,MAX(0.0_EB,(NOW_LL%T_OUT-T)) / (NOW_LL%T_OUT - NOW_LL%T_IN))
               ELSE
                  X_INT = 1.0_EB
               END IF
               HR%INTDOSE = DTSP*( (1.0_EB-X_INT)*PCX%FED_CO_CO2_O2(2) + X_INT*PCX%FED_CO_CO2_O2(1) ) + HR%INTDOSE
            END IF

            IF ( (NOW_LL%T_OUT) <= T) THEN
               CALL CHECK_TARGET_NODE(INODE,INODE2,ISTAT,XX,YY,ZZ,IOR_NEW, &
                    IMESH2,T,NEW_FFIELD_NAME,NEW_FFIELD_I,COLOR_INDEX, &
                    I_TARGET,KEEP_XY,ANGLE,STR_INDX,STR_SUB_INDX, HR)

               IF (ISTAT == 0 ) THEN
                  ! Put person to a new node, i.e., target node has empty space
                  HR%X = XX
                  HR%Y = YY
                  HR%Z = ZZ
                  HR%ANGLE = ANGLE

                  V = SQRT( HR%U**2 + HR%V**2 )
                  IF (V > 0.01_EB) V = PCX%FAC_SPEED*HR%V0_FAC*HR%SPEED
                  HR%IOR = 1        ! Remove from corridor
                  ! Put the person to a floor
                  IF (IOR_NEW == 1) THEN ! DOOR or ENTRY
                     HR%IMESH = IMESH2
                     HR%INODE = 0
                     DO N = 1, N_EGRIDS
                        IF (EVAC_NODE_LIST(N)%IMESH == IMESH2) HR%INODE = N
                     END DO
                     HR%NODE_NAME   = TRIM(MESH_NAME(IMESH2))
                     HR%FFIELD_NAME = TRIM(NEW_FFIELD_NAME)
                     HR%I_FFIELD = NEW_FFIELD_I
                     HR%I_TARGET = I_TARGET

                     IF (EVAC_NODE_LIST(INODE2)%NODE_TYPE == 'Door') THEN
                        IOR = - EVAC_DOORS( EVAC_NODE_LIST(INODE2)%NODE_INDEX)%IOR
                     END IF
                     IF (EVAC_NODE_LIST(INODE2)%NODE_TYPE == 'Entry') THEN
                        IOR = EVAC_ENTRYS( EVAC_NODE_LIST(INODE2)%NODE_INDEX)%IOR
                     END IF
                     IF ( ABS(IOR) == 1 ) THEN
                        HR%U = V*IOR
                        HR%V = 0.0_EB
                     END IF
                     IF ( ABS(IOR) == 2 ) THEN
                        HR%U = 0.0_EB
                        HR%V = 0.5_EB*V*IOR
                     END IF
                     HR%X_OLD = HR%X
                     HR%Y_OLD = HR%Y

                     IF (MESHES(IMESH2)%N_HUMANS+1 > MESHES(IMESH2)%N_HUMANS_DIM) THEN
                        ! Re-allocation is not yet checked.
                        WRITE(LU_ERR,*) 'ERROR: HUMANS: NO RE-ALLOCATION YET'
                        WRITE(LU_EVACOUT,*) 'ERROR: HUMANS: NO RE-ALLOCATION YET'
                        STOP_STATUS = EVACUATION_STOP ; RETURN
                        CALL RE_ALLOCATE_HUMANS(1,IMESH2)
                     END IF
                     MESHES(IMESH2)%N_HUMANS = MESHES(IMESH2)%N_HUMANS + 1
                     MESHES(IMESH2)%HUMAN(MESHES(IMESH2)%N_HUMANS) = HR
                     MESHES(IMESH2)%HUMAN(MESHES(IMESH2)%N_HUMANS)%IOR = HUMAN_NO_TARGET
                     WRITE (LU_EVACOUT,FMT='(A,I6,A,F8.2,A,A,A,F8.4)') ' Agent n:o', &
                          HR%ILABEL, ' out at ', T, ' s, corr ', TRIM(PCX%ID), ' FED ', HR%INTDOSE
                  END IF            ! Target is door or entry

                  IF (IOR_NEW == 3) THEN ! CORR
                     WRITE (LU_EVACOUT,FMT='(A,I6,A,F8.2,A,A,F8.4)') ' Agent n:o', &
                          HR%ILABEL, ' change corr ', T, ' s, corr ', TRIM(PCX%ID), HR%INTDOSE
                  END IF

                  IF (IOR_NEW == 5) THEN ! EXIT
                     WRITE (LU_EVACOUT,FMT='(A,I6,A,F8.2,A,A,F8.4)') ' Agent n:o', &
                          HR%ILABEL, ' exits ', T, ' s, corr ', TRIM(PCX%ID), HR%INTDOSE
                  END IF
               ELSE
                  ! Can not move to the next node, so do not allow to move inside
                  ! the DOOR ==> keep old position and put velocity to zero.
                  HR%X = HR%X_OLD
                  HR%Y = HR%Y_OLD
                  HR%I_TARGET = 0
                  HR%U = 0.0_EB
                  HR%V = 0.0_EB
                  HR%IOR = HUMAN_NO_TARGET
               END IF         ! ISTAT == 0

            ELSE             ! T_OUT > T, i.e., still in corridor
               HR%X = HR%X_OLD
               HR%Y = HR%Y_OLD
               HR%IOR = HUMAN_NO_TARGET
               HR%I_TARGET = 0
            END IF

            ! Here person is removed from the corridor (linked list)
            IF (HR%IOR > 0) THEN
               PCX%N_INSIDE = PCX%N_INSIDE - 1
               IF (ASSOCIATED(NOW_LL%NEXT)) THEN
                  TMP_LL => NOW_LL%NEXT
                  NULLIFY(HR)       ! HR points to NOW_LL%HUMAN
                  NULLIFY(NOW_LL%NEXT) ! Remove the pointer
                  DEALLOCATE(NOW_LL) ! Free memory
                  NULLIFY(NOW_LL)   ! Remove the pointer
                  NOW_LL => TMP_LL  ! Jump to the next element in LL
                  IF (.NOT. ASSOCIATED(PREV_LL)) THEN
                     PCX%FIRST => NOW_LL   ! First element removed
                  ELSE
                     PREV_LL%NEXT => NOW_LL ! Deleted element skipped
                  END IF
               ELSE                ! Remove the last element
                  NULLIFY(HR)       ! HR points to NOW_LL%HUMAN
                  NULLIFY(NOW_LL%NEXT) ! End list
                  DEALLOCATE(NOW_LL) ! Free memory
                  NULLIFY(NOW_LL)   ! End loop
                  IF (ASSOCIATED(PREV_LL)) THEN
                     NULLIFY(PREV_LL%NEXT) ! End the list
                  ELSE
                     NULLIFY(PCX%FIRST)   ! Empty list
                  END IF
               END IF
            ELSE
               PREV_LL => NOW_LL  ! Save the previous element
               IF (ASSOCIATED(NOW_LL%NEXT)) THEN
                  NOW_LL => NOW_LL%NEXT ! Jump to the next element in LL
               ELSE
                  NULLIFY(NOW_LL)  ! End the loop
               END IF
            END IF
         END DO HUMLOOP
      END DO PCXLOOP
      !
    END SUBROUTINE CHECK_CORRS
    !
    SUBROUTINE CHECK_TARGET_NODE(INODE,INODE2,ISTAT,XX,YY,ZZ,IOR_NEW,IMESH2,T,NEW_FFIELD_NAME,NEW_FFIELD_I, &
         COLOR_INDEX,I_TARGET,KEEP_XY,ANGLE,STR_INDX, STR_SUB_INDX, HR)
      IMPLICIT NONE
      !
      ! Check, that the target node is free.
      ! If target node is DOOR/ENTRY, try to put person to the floor.
      !
      ! IOR_NEW = 1: target is a DOOR or ENTRY
      !           3: Target is a corridor
      !           4: Target is a floor
      !           5: Target is an EXIT
      ! ISTAT   = 1: Target is not free (do not move there)
      !           0: Target is free
      !
      ! Passed variables
      INTEGER, INTENT(IN) :: INODE, INODE2
      INTEGER, INTENT(OUT) :: ISTAT, IOR_NEW, IMESH2, COLOR_INDEX, I_TARGET
      REAL(EB), INTENT(OUT) :: XX, YY, ZZ, ANGLE
      INTEGER, INTENT(OUT) :: STR_INDX, STR_SUB_INDX
      REAL(EB), INTENT(IN) :: T
      INTEGER, INTENT(INOUT) :: NEW_FFIELD_I
      LOGICAL, INTENT(IN) :: KEEP_XY
      CHARACTER(LABEL_LENGTH), INTENT(INOUT) :: NEW_FFIELD_NAME
      TYPE (HUMAN_TYPE), POINTER :: HR
      !
      ! Local variables
      REAL :: RN_REAL
      REAL(EB) RN, X1, X2, Y1, Y2, Z1, Z2, D_MAX, DIST, WIDTH, &
           XX1,YY1
      INTEGER  II, JJ, KK, IOR, IRNMAX, IRN, IE, IZERO, J1
      REAL(EB), DIMENSION(6) :: R_TMP, X_TMP, Y_TMP
      INTEGER :: I_TMP, III, JJJ, I_OBST
      LOGICAL :: KEEP_XY2, NM_STRS_MESH

      TYPE (CORR_LL_TYPE), POINTER :: TMPCURRENT =>NULL()
      TYPE (EVAC_STRS_TYPE), POINTER :: STRP =>NULL()
      TYPE (EVAC_DOOR_TYPE), POINTER :: PDX2 =>NULL()
      TYPE (EVACUATION_TYPE), POINTER :: HPT =>NULL()
      TYPE (HUMAN_TYPE), POINTER :: HRE =>NULL()
      TYPE (EVAC_ENTR_TYPE), POINTER :: PNX =>NULL(), PNX2 =>NULL()
      TYPE (EVAC_CORR_TYPE), POINTER :: PCX2 =>NULL()
      !
      XX = 0.0_EB ; YY = 0.0_EB ; ZZ = 0.0_EB
      I_TARGET = 0
      STR_INDX = 0
      STR_SUB_INDX = 0
      NM_STRS_MESH = .FALSE.

      ! Where is the person going to?
      SELECTTARGETTYPE: SELECT CASE (EVAC_NODE_LIST(INODE2)%NODE_TYPE)
      CASE ('Door','Entry') SELECTTARGETTYPE
         IOR_NEW = 1
         IF (EVAC_NODE_LIST(INODE2)%NODE_TYPE == 'Door') THEN
            PDX2 => EVAC_DOORS(EVAC_NODE_LIST(INODE2)%NODE_INDEX)
            IMESH2  = PDX2%IMESH
            NEW_FFIELD_NAME = TRIM(PDX2%GRID_NAME)
            IRN_LOOP1: DO IRN = 1, NMESHES
               IF ( EVACUATION_ONLY(IRN) .AND. TRIM(PDX2%GRID_NAME) == TRIM(MESH_NAME(irn)) ) THEN
                  new_ffield_i = irn
                  EXIT Irn_Loop1
               END IF
            END DO Irn_Loop1
            X1  = PDX2%X1
            X2  = PDX2%X2
            Y1  = PDX2%Y1
            Y2  = PDX2%Y2
            Z1  = PDX2%Z1
            Z2  = PDX2%Z2
            ior = -PDX2%IOR        ! Now entering the room
            Width = PDX2%Width
         ELSE
            PNX2 => EVAC_ENTRYS(EVAC_Node_List(INODE2)%Node_Index)
            imesh2  = PNX2%IMESH
            new_ffield_name = TRIM(PNX2%GRID_NAME)
            Irn_Loop2: DO irn = 1, NMESHES
               IF ( EVACUATION_ONLY(IRN) .AND. TRIM(PNX2%GRID_NAME) == TRIM(MESH_NAME(irn)) ) THEN
                  new_ffield_i = irn
                  EXIT Irn_Loop2
               END IF
            END DO Irn_Loop2
            X1  = PNX2%X1
            X2  = PNX2%X2
            Y1  = PNX2%Y1
            Y2  = PNX2%Y2
            Z1  = PNX2%Z1
            Z2  = PNX2%Z2
            ior = PNX2%IOR
            Width = PNX2%Width
            STR_INDX = PNX2%STR_INDX
            STR_SUB_INDX = PNX2%STR_SUB_INDX
            IF (STR_INDX>0) THEN
               Z1 = EVAC_STRS(STR_INDX)%XB_NODE(STR_SUB_INDX,5)
               Z2 = EVAC_STRS(STR_INDX)%XB_NODE(STR_SUB_INDX,6)
            END IF
         END IF
         IF (ABS(ior) == 1) irnmax = INT(Width*4.0_EB)
         IF (ABS(ior) == 2) irnmax = INT(Width*4.0_EB)
         IF (ABS(ior) == 3 .OR. ior == 0) irnmax = INT((x2-x1)*4.0_EB)*INT((y2-y1)*4.0_EB)
         irnmax = MAX(irnmax,5)
         istat = 1
         MFF=>MESHES(imesh2)
         DO N = 1, N_STRS
            IF (EVAC_STRS(N)%IMESH == IMESH2) THEN
               STRS_Indx = N
               STRP=>EVAC_STRS(N)
               NM_STRS_MESH = .TRUE.
               EXIT
            END IF
         END DO
         !
         irn = 0
         angle = 0.0_EB
         CheckPPForce: DO WHILE (irn < irnmax)
            SELECT CASE (ior)
            CASE(-1,1)
               ! 180 or 0 degrees, i.e., pi or 0 radians
               angle = (0.5_EB-(ior/2.0_EB))*Pi
               CALL RANDOM_NUMBER(RN_REAL)
               RN = REAL(RN_REAL,EB)
               IF (keep_xy) THEN
                  yy = HR%Y + (irn/MAX(irn,1))*(rn-0.5_EB)*0.25_EB*HR%Radius
                  xx = x1 + ior*(1.0_EB*HR%B + HR%Radius)
                  angle = HR%Angle
               ELSE
                  yy = 0.5_EB*(y1+y2) + (rn-0.5_EB)*MAX(0.0_EB,Width-2.0_EB*HR%Radius-2.0_EB*HR%B)
                  xx = x1 + ior*(1.0_EB*HR%B + HR%r_torso)
               END IF
            CASE(-2,2)
               ! 270 or 90 degrees, i.e., 3pi/2 or pi/2 radians
               angle = (1.0_EB-(ior/4.0_EB))*Pi
               CALL RANDOM_NUMBER(RN_REAL)
               RN = REAL(RN_REAL,EB)
               IF (keep_xy) THEN
                  xx = HR%X + (irn/MAX(irn,1))*(rn-0.5_EB)*0.25_EB*HR%Radius
                  yy = y1 + (ior/ABS(ior))*(1.0_EB*HR%B + HR%Radius)
                  angle = HR%Angle
               ELSE
                  xx = 0.5_EB*(x1+x2) + (rn-0.5_EB)*MAX(0.0_EB,Width-2.0_EB*HR%Radius-2.0_EB*HR%B)
                  yy = y1 + (ior/ABS(ior))*(1.0_EB*HR%B + HR%r_torso)
               END IF
            CASE(0,3)
               CALL RANDOM_NUMBER(RN_REAL)
               RN = REAL(RN_REAL,EB)
               yy = 0.5_EB*(y1+y2) + (rn-0.5_EB)*MAX(0.0_EB,(y2-y1)-2.0_EB*HR%Radius-2.0_EB*HR%B)
               CALL RANDOM_NUMBER(RN_REAL)
               RN = REAL(RN_REAL,EB)
               xx = 0.5_EB*(x1+x2) + (rn-0.5_EB)*MAX(0.0_EB,(x2-x1)-2.0_EB*HR%Radius-2.0_EB*HR%B)
            END SELECT
            zz = 0.5_EB*(z1+z2)
            II = FLOOR( MFF%CELLSI(FLOOR((xx-MFF%XS)*MFF%RDXINT)) + 1.0_EB )
            JJ = FLOOR( MFF%CELLSJ(FLOOR((yy-MFF%YS)*MFF%RDYINT)) + 1.0_EB )
            KK = 1
            I_OBST = MFF%OBST_INDEX_C(MFF%CELL_INDEX(II,JJ,KK))

            irn = irn + 1

            IF (MFF%SOLID(MFF%CELL_INDEX(II,JJ,KK)) .AND. .NOT.MFF%OBSTRUCTION(I_OBST)%HIDDEN) CYCLE CheckPPForce
            IF ( ABS(ior) == 2 .AND. .NOT. keep_xy ) THEN
               xx1 = xx - HR%Radius - 1.0_EB*HR%B
               II = FLOOR(MFF%CELLSI(FLOOR((xx1-MFF%XS)*MFF%RDXINT))+1.0_EB)
               I_OBST = MFF%OBST_INDEX_C(MFF%CELL_INDEX(II,JJ,KK))
               IF (MFF%SOLID(MFF%CELL_INDEX(II,JJ,KK)) .AND. .NOT.MFF%OBSTRUCTION(I_OBST)%HIDDEN) CYCLE CheckPPForce
               xx1 = xx + HR%Radius + 1.0_EB*HR%B
               II = FLOOR(MFF%CELLSI(FLOOR((xx1-MFF%XS)*MFF%RDXINT))+1.0_EB)
               I_OBST = MFF%OBST_INDEX_C(MFF%CELL_INDEX(II,JJ,KK))
               IF (MFF%SOLID(MFF%CELL_INDEX(II,JJ,KK)) .AND. .NOT.MFF%OBSTRUCTION(I_OBST)%HIDDEN) CYCLE CheckPPForce
            END IF
            IF ( ABS(ior) == 1 .AND. .NOT. keep_xy ) THEN
               yy1 = yy - HR%Radius - 1.0_EB*HR%B
               JJ = FLOOR(MFF%CELLSJ(FLOOR((yy1-MFF%YS)*MFF%RDYINT))+1.0_EB)
               I_OBST = MFF%OBST_INDEX_C(MFF%CELL_INDEX(II,JJ,KK))
               IF (MFF%SOLID(MFF%CELL_INDEX(II,JJ,KK)) .AND. .NOT.MFF%OBSTRUCTION(I_OBST)%HIDDEN) CYCLE CheckPPForce
               yy1 = yy + HR%Radius + 1.0_EB*HR%B
               JJ = FLOOR(MFF%CELLSJ(FLOOR((yy1-MFF%YS)*MFF%RDYINT))+1.0_EB)
               I_OBST = MFF%OBST_INDEX_C(MFF%CELL_INDEX(II,JJ,KK))
               IF (MFF%SOLID(MFF%CELL_INDEX(II,JJ,KK)) .AND. .NOT.MFF%OBSTRUCTION(I_OBST)%HIDDEN) CYCLE CheckPPForce
            END IF

            !
            d_max =  1.0_EB*HR%B
            r_tmp(1) = HR%r_shoulder ! right circle
            r_tmp(2) = HR%r_torso     ! center circle
            r_tmp(3) = HR%r_shoulder ! left circle
            y_tmp(1) = yy - COS(angle)*HR%d_shoulder ! right
            x_tmp(1) = xx + SIN(angle)*HR%d_shoulder
            y_tmp(2) = yy ! torso
            x_tmp(2) = xx
            y_tmp(3) = yy + COS(angle)*HR%d_shoulder ! left
            x_tmp(3) = xx - SIN(angle)*HR%d_shoulder
            P2PLoop: DO ie = 1, MFF%N_HUMANS
               HRE=>MFF%HUMAN(IE)
               IF (NM_STRS_MESH .AND. STR_SUB_INDX /= HRE%STR_SUB_INDX) CYCLE P2PLoop
               r_tmp(4) = HRE%r_shoulder ! right circle
               r_tmp(5) = HRE%r_torso     ! center circle
               r_tmp(6) = HRE%r_shoulder ! left circle
               y_tmp(4) = HRE%Y - COS(HRE%angle)*HRE%d_shoulder ! right circle
               x_tmp(4) = HRE%X + SIN(HRE%angle)*HRE%d_shoulder
               y_tmp(5) = HRE%Y ! center circle
               x_tmp(5) = HRE%X
               y_tmp(6) = HRE%Y + COS(HRE%angle)*HRE%d_shoulder ! left circle
               x_tmp(6) = HRE%X - SIN(HRE%angle)*HRE%d_shoulder
               DO iii = 1, 3
                  DO jjj = 4, 6
                     DIST = SQRT((x_tmp(jjj)-x_tmp(iii))**2 + (y_tmp(jjj)-y_tmp(iii))**2) - (r_tmp(jjj)+r_tmp(iii))
                     IF ( DIST < d_max) CYCLE CheckPPForce
                  END DO
               END DO
            END DO P2PLoop
            istat = 0    ! target is free
            ior_new = 1  ! target is a door or an entry
            EXIT CheckPPForce
         END DO CheckPPForce

         IF (istat == 0) THEN
            j  =  MAX(0,HR%GROUP_ID)
            j1 = -MIN(0,HR%GROUP_ID)
            IF (j > 0) THEN
               Group_Known_Doors(j)%I_Target = 0
               P2PLoop2: DO ie = 1, MFF%N_HUMANS
                  HRE=>MFF%HUMAN(IE)
                  IF (HRE%GROUP_ID == j .AND. (HRE%Ilabel /= HR%Ilabel)) THEN
                     Group_Known_Doors(j)%I_Target = HRE%I_Target
                     EXIT P2PLoop2
                  END IF
               END DO P2PLoop2
            END IF
            IF (HR%IEL > 0 ) THEN
               HPT => EVAC_EVACS(HR%IEL)
            ELSE
               PNX => EVAC_ENTRYS(ABS(HR%IEL))
            END IF

            i_tmp = 0
            DO ie = 1, n_egrids
               IF (EVAC_Node_List(ie)%IMESH == imesh2 ) THEN
                  i_tmp = ie
               END IF
            END DO
            IF (i_tmp == 0 .OR. i_tmp > n_egrids) THEN
               WRITE(LU_ERR,'(A,I6)') 'ERROR: Check_Target_Node, no imesh2 found ',imesh2
               WRITE(LU_EVACOUT,'(A,I6)') 'ERROR: Check_Target_Node, no imesh2 found ',imesh2
               STOP_STATUS = EVACUATION_STOP ; RETURN
            END IF

            IF (j > 0 .AND. Group_List(j)%GROUP_I_FFIELDS(i_tmp) > 0) THEN
               ! There are already group members on this floor, use the same field
               new_ffield_i = Group_List(j)%GROUP_I_FFIELDS(i_tmp)
               new_ffield_name = TRIM(MESH_NAME(new_ffield_i))
               IF (j > 0) I_Target = Group_Known_Doors(j)%I_Target
            ELSE

               HR%I_Target = 0
               CALL Change_Target_Door(imesh2, imesh2, 1, j, j1, 0, 2, xx, yy, I_Target, color_index, new_ffield_i, HR)

               new_ffield_name = TRIM(MESH_NAME(new_ffield_i))
               IF ( j > 0 ) THEN
                  Group_List(j)%GROUP_I_FFIELDS(i_tmp) = new_ffield_i
               END IF
            END IF  ! first member of a group?

         END IF ! istat=0, i.e., put human to a new node

      CASE ('Corr') SelectTargetType
         PCX2 => EVAC_CORRS(EVAC_Node_List(INODE2)%Node_Index)
         IF ((.NOT.PCX2%ELEVATOR .AND. PCX2%n_inside < PCX2%MAX_HUMANS_INSIDE) .OR. &
              (PCX2%ELEVATOR .AND. PCX2%Elevator_Mode == 2 .AND. PCX2%n_inside < PCX2%MAX_HUMANS_INSIDE)) THEN
            xx = 0.0_EB ; yy = 0.0_EB ; zz = 0.0_EB ! not used for corridors
            ior_new = 3     ! target is a corridor
            imesh2 = HR%IMESH     ! do not put to an another floor yet
            HR%INODE=evac_corrs(EVAC_Node_List(INODE2)%Node_Index)%INODE
            istat = 0       ! target is free
            PCX2%n_inside = PCX2%n_inside + 1
            ALLOCATE(TmpCurrent,STAT=IZERO)
            CALL ChkMemErr('LINK','Check_Target_Node',IZERO)
            TmpCurrent%HUMAN  = HR    ! Save the human data
            TmpCurrent%T_in   = T
            TmpCurrent%T_out  = T + PCX2%Eff_Length/(PCX2%Fac_Speed*HR%v0_fac*HR%Speed)
            IF (PCX2%ELEVATOR) THEN
               PCX2%Time_Last_In = T
               PCX2%Time_Start   = T + PCX2%Time_Wait_To_Fill
               IF (PCX2%n_inside == PCX2%MAX_HUMANS_INSIDE) THEN
                  PCX2%Time_Start = T
                  PCX2%Elevator_Mode = 0 ! Start moving down
                  PCX2%Time_Stop = T + PCX2%Travel_Time
               END IF
               TmpCurrent%T_out  = PCX2%Time_Start + PCX2%Travel_Time
            END IF
            TmpCurrent%From1_To2 = .TRUE.
            TmpCurrent%Index     = PCX2%n_inside
            TmpCurrent%Next      => PCX2%First
            PCX2%First            => TmpCurrent
         ELSE
            ior_new = 3     ! target is a corridor
            istat = 1       ! do not enter the corridor
            imesh2 = HR%IMESH
         END IF
      CASE ('Floor') SelectTargetType
         ! Next is used for meshes
         !
         NM_STRS_MESH = .FALSE.
         ior_new = 4       ! target is a mesh
         istat = 1         ! Default: no space in the mesh
         SELECT CASE (EVAC_NODE_List(INODE)%Node_Type)
         CASE ('Door')
            PDX2 => EVAC_DOORS(EVAC_Node_List(INODE)%Node_Index)
            imesh2 = EVAC_Node_List(INODE2)%IMESH
            MFF=>MESHES(imesh2)
            DO N = 1, N_STRS
               IF (EVAC_STRS(N)%IMESH == imesh2) THEN
                  STRS_Indx = N
                  STRP=>EVAC_STRS(N)
                  NM_STRS_MESH = .TRUE.
                  EXIT
               END IF
            END DO
            X1  = MIN(MAX(MFF%XS,PDX2%X1),MFF%XF)
            X2  = MIN(MAX(MFF%XS,PDX2%X2),MFF%XF)
            Y1  = MIN(MAX(MFF%YS,PDX2%Y1),MFF%YF)
            Y2  = MIN(MAX(MFF%YS,PDX2%Y2),MFF%YF)
            Z1  = HR%Z   ! use agent's z
            Z2  = HR%Z   ! use agent's z
            zz  = HR%Z   ! use agent's z
            ior = PDX2%IOR        ! direction is mesh1==>door1==>mesh2
            Width = PDX2%Width
            keep_xy2 = .TRUE. ! do not move horizontally
            IF (PDX2%STR_INDX > 0) THEN
               STR_INDX = PDX2%STR_INDX
               STR_SUB_INDX = PDX2%STR_SUB_INDX
            END IF
         CASE ('Entry')
            ! Todo:  Case ('Entry') and geometry part from CHECK_ENTRY
            PNX2 => EVAC_ENTRYS(EVAC_Node_List(INODE2)%Node_Index)
            IF (PNX2%STR_INDX > 0) THEN
               STR_INDX = PNX2%STR_INDX
               STR_SUB_INDX = PNX2%STR_SUB_INDX
            END IF
         CASE Default
            WRITE(MESSAGE,'(A)') 'ERROR (debug): Check_Target_Node and STRS: Node 1 is not a door.'
         END SELECT
         IF (ABS(ior) == 1) irnmax = INT(Width*4.0_EB)
         IF (ABS(ior) == 2) irnmax = INT(Width*4.0_EB)
         IF (ABS(ior) == 3 .OR. ior == 0) irnmax = INT((x2-x1)*4.0_EB)*INT((y2-y1)*4.0_EB)
         irnmax = MAX(irnmax,5)
         !
         irn = 0
         angle = 0.0_EB
         CheckPPForceSTRS: DO WHILE (irn < irnmax)
            SELECT CASE (ior)
            CASE(-1,1)
               ! 180 or 0 degrees, i.e., pi or 0 radians
               angle = (0.5_EB-(ior/2.0_EB))*Pi
               CALL RANDOM_NUMBER(RN_REAL)
               RN = REAL(RN_REAL,EB)
               IF (keep_xy2) THEN
                  yy = HR%Y + (irn/MAX(irn,1))*(rn-0.5_EB)*0.25_EB*HR%Radius
                  xx = x1 + ior*(1.0_EB*HR%B + HR%Radius)
                  angle = HR%Angle
               ELSE
                  yy = 0.5_EB*(y1+y2) + (rn-0.5_EB)* MAX(0.0_EB,Width-2.0_EB*HR%Radius-2.0_EB*HR%B)
                  xx = x1 + ior*(1.0_EB*HR%B + HR%r_torso)
               END IF
            CASE(-2,2)
               ! 270 or 90 degrees, i.e., 3pi/2 or pi/2 radians
               angle = (1.0_EB-(ior/4.0_EB))*Pi
               CALL RANDOM_NUMBER(RN_REAL)
               RN = REAL(RN_REAL,EB)
               IF (keep_xy2) THEN
                  xx = HR%X + (irn/MAX(irn,1))*(rn-0.5_EB)*0.25_EB*HR%Radius
                  yy = y1 + (ior/ABS(ior))*(1.0_EB*HR%B + HR%Radius)
                  angle = HR%Angle
               ELSE
                  xx = 0.5_EB*(x1+x2) + (rn-0.5_EB)* MAX(0.0_EB,Width-2.0_EB*HR%Radius-2.0_EB*HR%B)
                  yy = y1 + (ior/ABS(ior))*(1.0_EB*HR%B + HR%r_torso)
               END IF
            CASE(0,3)
               CALL RANDOM_NUMBER(RN_REAL)
               RN = REAL(RN_REAL,EB)
               yy = 0.5_EB*(y1+y2) + (rn-0.5_EB)* MAX(0.0_EB,(y2-y1)-2.0_EB*HR%Radius-2.0_EB*HR%B)
               CALL RANDOM_NUMBER(RN_REAL)
               RN = REAL(RN_REAL,EB)
               xx = 0.5_EB*(x1+x2) + (rn-0.5_EB)* MAX(0.0_EB,(x2-x1)-2.0_EB*HR%Radius-2.0_EB*HR%B)
            END SELECT
            II = FLOOR( MFF%CELLSI(FLOOR((xx-MFF%XS)*MFF%RDXINT)) + 1.0_EB )
            JJ = FLOOR( MFF%CELLSJ(FLOOR((yy-MFF%YS)*MFF%RDYINT)) + 1.0_EB )
            KK = 1

            irn = irn + 1

            I_OBST = MFF%OBST_INDEX_C(MFF%CELL_INDEX(II,JJ,KK))
            IF (MFF%SOLID(MFF%CELL_INDEX(II,JJ,KK)) .AND. .NOT.MFF%OBSTRUCTION(I_OBST)%HIDDEN) CYCLE CheckPPForceSTRS

            !
            d_max =  1.0_EB*HR%B
            r_tmp(1) = HR%r_shoulder ! right circle
            r_tmp(2) = HR%r_torso     ! center circle
            r_tmp(3) = HR%r_shoulder ! left circle
            y_tmp(1) = yy - COS(angle)*HR%d_shoulder ! right
            x_tmp(1) = xx + SIN(angle)*HR%d_shoulder
            y_tmp(2) = yy ! torso
            x_tmp(2) = xx
            y_tmp(3) = yy + COS(angle)*HR%d_shoulder ! left
            x_tmp(3) = xx - SIN(angle)*HR%d_shoulder
            P2PLoopSTRS: DO ie = 1, MFF%N_HUMANS
               HRE=>MFF%HUMAN(IE)
               IF (ABS(STR_SUB_INDX-HRE%STR_SUB_INDX)>1) CYCLE P2PLoopSTRS
               r_tmp(4) = HRE%r_shoulder ! right circle
               r_tmp(5) = HRE%r_torso     ! center circle
               r_tmp(6) = HRE%r_shoulder ! left circle
               y_tmp(4) = HRE%Y - COS(HRE%angle)*HRE%d_shoulder ! right circle
               x_tmp(4) = HRE%X + SIN(HRE%angle)*HRE%d_shoulder
               y_tmp(5) = HRE%Y ! center circle
               x_tmp(5) = HRE%X
               y_tmp(6) = HRE%Y + COS(HRE%angle)*HRE%d_shoulder ! left circle
               x_tmp(6) = HRE%X - SIN(HRE%angle)*HRE%d_shoulder
               DO iii = 1, 3
                  DO jjj = 4, 6
                     DIST = SQRT((x_tmp(jjj)-x_tmp(iii))**2 + (y_tmp(jjj)-y_tmp(iii))**2) - (r_tmp(jjj)+r_tmp(iii))
                     IF ( DIST < d_max) CYCLE CheckPPForceSTRS
                  END DO
               END DO
            END DO P2PLoopSTRS
            istat = 0    ! target is free
            EXIT CheckPPForceSTRS
         END DO CheckPPForceSTRS
         IF (istat == 0) THEN
            IF (NM_STRS_MESH) CALL Find_Target_Node_In_Strs(STRP,HR)
            I_target = HR%I_Target
            new_ffield_i    = imesh2
            new_ffield_name = TRIM(MESH_NAME(new_ffield_i))
         END IF ! istat=0, i.e., put human to a new node
      CASE ('Exit') SelectTargetType
         ior_new = 5  ! target is an exit
         istat = 0    ! remove from the floor (exit is always free)
         imesh2 = HR%IMESH
      CASE Default SelectTargetType
         ior_new = 6  ! target is not defined
         istat = 1    ! do not remove from the floor
         imesh2 = HR%IMESH
      END SELECT SelectTargetType
      !
    END SUBROUTINE Check_Target_Node
    !
    SUBROUTINE REMOVE_PERSON(I)
      IMPLICIT NONE
      !
      ! Remove a person
      !
      ! Passed variables
      INTEGER, INTENT(IN) :: I
      !
      ! Local variables
      !
      HUMAN(I) = HUMAN(N_HUMANS)
      N_HUMANS = N_HUMANS - 1
      !
    END SUBROUTINE REMOVE_PERSON
    !
    SUBROUTINE REMOVE_OUT_OF_GRIDS
      IMPLICIT NONE
      !
      ! Remove humans that do not lie in any mesh
      !
      INTEGER :: IKILL, I
      !
      IKILL = 0
      DROP_LOOP: DO I=1,N_HUMANS
         !
         HR=>HUMAN(I)
         IF (I > N_HUMANS-IKILL) EXIT DROP_LOOP
         IF (HR%X > XS .AND. HR%X < XF .AND. HR%Y > YS .AND. HR%Y < YF) CYCLE DROP_LOOP
         WRITE(LU_ERR,'(A,I6,A,3F8.2,I3)') 'WARNING: Agent n:o ',HR%ILABEL,' removed at coord: ', &
              HR%X,HR%Y,HR%Z,HR%SKIP_WALL_FORCE_IOR
         HUMAN(I) = HUMAN(N_HUMANS-IKILL)
         IKILL = IKILL + 1
      END DO DROP_LOOP
      N_HUMANS = N_HUMANS - IKILL

!!$      IKILL = 0
!!$      CROWBAR_DROP_LOOP: DO I=1,N_HUMANS
!!$         HR=>HUMAN(I)
!!$         IF (.NOT.HR%CROWBAR_READ_IN) CYCLE CROWBAR_DROP_LOOP
!!$         IF (I > N_HUMANS-IKILL) EXIT CROWBAR_DROP_LOOP
!!$         IF (HR%IEL>=0) CYCLE CROWBAR_DROP_LOOP ! not from an entry
!!$         IF (HR%CROWBAR_READ_IN .AND. HR%T_LastRead_CB > T-1.5_EB*CROWBAR_DT_READ) CYCLE CROWBAR_DROP_LOOP
!!$         ! IF (HR%CROWBAR_READ_IN .AND. HR%T_LastRead_CB > T-2.5_EB*CROWBAR_DT_READ) CYCLE CROWBAR_DROP_LOOP
!!$         ! Now just those crowbar agents that have not lately been read in (updated)
!!$
!!$         ! Is the agent (x,y) outside or close to the ENTR bounds?
!!$         DIST_CB = HUGE(DIST_CB)
!!$         X1 = EVAC_ENTRYS(ABS(HR%IEL))%X1; X2 = EVAC_ENTRYS(ABS(HR%IEL))%X2;
!!$         Y1 = EVAC_ENTRYS(ABS(HR%IEL))%Y1; Y2 = EVAC_ENTRYS(ABS(HR%IEL))%Y2;
!!$         DIST_CB = MIN(DIST_CB,HR%X-X1) ! The -x face, how much inside of XB
!!$         DIST_CB = MIN(DIST_CB,X2-HR%X) ! The +x face, how much inside of XB
!!$         DIST_CB = MIN(DIST_CB,HR%Y-Y1) ! The -y face, how much inside of XB
!!$         DIST_CB = MIN(DIST_CB,Y2-HR%Y) ! The +y face, how much inside of XB
!!$         ! Agent speed = 1.0 m/s, dt=0.2s => more than 0.1 m out of XB => do not "remove"
!!$         IF (DIST_CB < 0.5_EB*CROWBAR_DT_READ) THEN
!!$            HR%CROWBAR_READ_IN = .FALSE.  ! Out of the camera area
!!$            HR%CROWBAR_UPDATE_V0 = .TRUE. ! Use still the last camera information for speed for a while
!!$            CYCLE CROWBAR_DROP_LOOP ! Is not in the XB area anymore
!!$         END IF
!!$
!!$         ! Is the agent (x,y)_readin outside or close to the ENTR bounds?
!!$         DIST_CB = HUGE(DIST_CB)
!!$         X1 = EVAC_ENTRYS(ABS(HR%IEL))%X1; X2 = EVAC_ENTRYS(ABS(HR%IEL))%X2;
!!$         Y1 = EVAC_ENTRYS(ABS(HR%IEL))%Y1; Y2 = EVAC_ENTRYS(ABS(HR%IEL))%Y2;
!!$         DIST_CB = MIN(DIST_CB,HR%X_CB-X1) ! The -x face, how much inside of XB
!!$         DIST_CB = MIN(DIST_CB,X2-HR%X_CB) ! The +x face, how much inside of XB
!!$         DIST_CB = MIN(DIST_CB,HR%Y_CB-Y1) ! The -y face, how much inside of XB
!!$         DIST_CB = MIN(DIST_CB,Y2-HR%Y_CB) ! The +y face, how much inside of XB
!!$         ! Agent speed = 1.0 m/s, dt=0.2s => more than 0.1 m out of XB => do not remove
!!$         IF (DIST_CB < 0.5_EB*CROWBAR_DT_READ) THEN
!!$            HR%CROWBAR_READ_IN = .FALSE.  ! Out of the camera area
!!$            HR%CROWBAR_UPDATE_V0 = .TRUE. ! Use still the last camera information for speed for a while
!!$            CYCLE CROWBAR_DROP_LOOP ! Is not in the XB area anymore
!!$         END IF
!!$      END DO CROWBAR_DROP_LOOP
!!$      N_HUMANS = N_HUMANS - IKILL

    END SUBROUTINE REMOVE_OUT_OF_GRIDS
    !
    SUBROUTINE ENTRY_HUMAN(I_entry, Tin, NM, istat)
      IMPLICIT NONE
      !
      ! Insert humans into the domain every 1/Flow seconds.
      !
      ! Passed variables
      REAL(EB), INTENT(IN) :: Tin
      INTEGER,  INTENT(IN) :: I_entry, NM
      INTEGER, INTENT(OUT) :: istat
      !
      ! Local variables
      REAL :: RN_REAL
      REAL(EB) RN, x1, x2, y1, y2, z1, z2, d_max, dist, xx, yy, zz, xx1, yy1
      INTEGER  I,II, JJ, KK, ior, irnmax, irn, ie, NR
      REAL(EB), DIMENSION(6) ::y_tmp, x_tmp, r_tmp

      TYPE (EVAC_ENTR_TYPE), POINTER :: PNX =>NULL()
      TYPE (MESH_TYPE), POINTER :: MFF =>NULL()
      TYPE (EVAC_PERS_TYPE), POINTER :: PCP =>NULL()
      TYPE (HUMAN_TYPE), POINTER :: HR=>NULL(), HRE =>NULL()
      !
      istat = 1
      PNX => EVAC_ENTRYS(I_entry)
      IF (PNX%CROWBAR_READ_IN) RETURN
      IF (PNX%IMESH /= NM ) RETURN
      IF (PNX%T_Start > Tin) RETURN
      IF (PNX%T_Stop < Tin) RETURN
      IF (PNX%Max_Humans < 0) THEN
         NR     = -PNX%Max_Humans
         IF (PNX%ICOUNT >= INT(EVALUATE_RAMP(Tin,0._EB,NR))) RETURN
      ELSE
         IF (PNX%ICOUNT >= PNX%Max_Humans) RETURN
         IF ( (Tin-PNX%T_last) < (1.0_EB/Max(0.0001_EB,PNX%Flow)) ) RETURN
      ENDIF
      MFF => MESHES(NM)
      X1  = PNX%X1
      X2  = PNX%X2
      Y1  = PNX%Y1
      Y2  = PNX%Y2
      Z1  = PNX%Z1
      Z2  = PNX%Z2
      ior = PNX%IOR
      !
      IF (N_HUMANS+1 > N_HUMANS_DIM) THEN
         ! Re-allocation is not yet checked.
         WRITE(LU_ERR,*) 'ERROR: INSERT HUMANS: NO RE-ALLOCATION YET'
         WRITE(LU_EVACOUT,*) 'ERROR: INSERT HUMANS: NO RE-ALLOCATION YET'
         STOP_STATUS = EVACUATION_STOP ; RETURN
         HUMAN=>MESHES(NM)%HUMAN
      END IF
      !
      PCP => EVAC_PERSON_CLASSES(PNX%IPC)
      HR  => MESHES(NM)%HUMAN(N_HUMANS+1)
      CALL CLASS_PROPERTIES(HR,PCP,-I_entry)
      HR%Tpre = 0.0_EB
      HR%Tdet = T_BEGIN
      HR%IPC  = PNX%IPC
      HR%IEL  = -I_entry
      HR%GROUP_ID = 0
      HR%I_Target = 0
      HR%I_DoorAlgo = PNX%I_AGENT_TYPE
      HR%I_Door_Mode = 0 ! Default, no target door yet
      CALL RANDOM_NUMBER(RN_REAL)
      RN = REAL(RN_REAL,EB)
      HR%F_FallDown = F_MIN_FALL + RN*(F_MAX_FALL - F_MIN_FALL)
      HR%T_FallenDown = HUGE(HR%T_FallenDown)
      HR%Angle_FallenDown = 0.0_EB
      HR%SizeFac_FallenDown = 0.0_EB
      HR%T_CheckFallDown = Tin
      !
      IF (ABS(ior) == 1) irnmax = INT(PNX%Width*4.0_EB)
      IF (ABS(ior) == 2) irnmax = INT(PNX%Width*4.0_EB)
      IF (ABS(ior) == 3 .OR. ior == 0) irnmax = INT((x2-x1)*4.0_EB)*INT((y2-y1)*4.0_EB)
      irnmax = MAX(irnmax,5)
      !
      irn = 0
      CheckPPForce: DO WHILE (irn < irnmax)
         SELECT CASE (ior)
         CASE(-1,1)
            CALL RANDOM_NUMBER(RN_REAL)
            RN = REAL(RN_REAL,EB)
            yy = 0.5_EB*(y1+y2) + (rn-0.5_EB)*MAX(0.0_EB,PNX%Width-2.0_EB*HR%Radius-2.0_EB*HR%B)
            xx = x1 + ior*5.0_EB*HR%B
            HR%Angle = (1-ior)*Pi/2.0_EB  ! ior=1: 0,  ior=-1: pi
         CASE(-2,2)
            CALL RANDOM_NUMBER(RN_REAL)
            RN = REAL(RN_REAL,EB)
            xx = 0.5_EB*(x1+x2) + (rn-0.5_EB)*MAX(0.0_EB,PNX%Width-2.0_EB*HR%Radius-2.0_EB*HR%B)
            yy = y1 + (ior/ABS(ior))*5.0_EB*HR%B
            HR%Angle = Pi/2.0_EB + (2-ior)*Pi/4.0_EB  ! ior=2: (3/2)pi,  ior=-2: pi/2
         CASE(0,3)
            CALL RANDOM_NUMBER(RN_REAL)
            RN = REAL(RN_REAL,EB)
            yy = 0.5_EB*(y1+y2) + (rn-0.5_EB)*MAX(0.0_EB,(y2-y1)-2.0_EB*HR%Radius-2.0_EB*HR%B)
            CALL RANDOM_NUMBER(RN_REAL)
            RN = REAL(RN_REAL,EB)
            xx = 0.5_EB*(x1+x2) + (rn-0.5_EB)*MAX(0.0_EB,(x2-x1)-2.0_EB*HR%Radius-2.0_EB*HR%B)
            CALL RANDOM_NUMBER(RN_REAL)
            RN = REAL(RN_REAL,EB)
            HR%Angle = 2.0_EB*Pi*rn
         END SELECT
         DO WHILE (HR%Angle >= 2.0_EB*Pi)
            HR%Angle = HR%Angle - 2.0_EB*Pi
         END DO
         DO WHILE (HR%Angle < 0.0_EB)
            HR%Angle = HR%Angle + 2.0_EB*Pi
         END DO
         zz  = 0.5_EB*(z1+z2)
         II = FLOOR( CELLSI(FLOOR((xx-XS)*RDXINT)) + 1.0_EB )
         JJ = FLOOR( CELLSJ(FLOOR((yy-YS)*RDYINT)) + 1.0_EB )
         KK = 1

         irn = irn + 1

         I_OBST = OBST_INDEX_C(CELL_INDEX(II,JJ,KK))
         IF (SOLID(CELL_INDEX(II,JJ,KK)) .AND. .NOT.OBSTRUCTION(I_OBST)%HIDDEN) CYCLE CheckPPForce
         IF ( ABS(ior) == 2 ) THEN
            xx1 = xx - HR%Radius - 1.0_EB*HR%B
            II = FLOOR(CELLSI(FLOOR((xx1-XS)*RDXINT))+1.0_EB)
            I_OBST = OBST_INDEX_C(CELL_INDEX(II,JJ,KK))
            IF (SOLID(CELL_INDEX(II,JJ,KK)) .AND. .NOT.OBSTRUCTION(I_OBST)%HIDDEN) CYCLE CheckPPForce
            xx1 = xx + HR%Radius + 1.0_EB*HR%B
            II = FLOOR(CELLSI(FLOOR((xx1-XS)*RDXINT))+1.0_EB)
            I_OBST = OBST_INDEX_C(CELL_INDEX(II,JJ,KK))
            IF (SOLID(CELL_INDEX(II,JJ,KK)) .AND. .NOT.OBSTRUCTION(I_OBST)%HIDDEN) CYCLE CheckPPForce
         END IF
         IF ( ABS(ior) == 1 ) THEN
            yy1 = yy - HR%Radius - 1.0_EB*HR%B
            JJ = FLOOR(CELLSJ(FLOOR((yy1-YS)*RDYINT))+1.0_EB)
            I_OBST = OBST_INDEX_C(CELL_INDEX(II,JJ,KK))
            IF (SOLID(CELL_INDEX(II,JJ,KK)) .AND. .NOT.OBSTRUCTION(I_OBST)%HIDDEN) CYCLE CheckPPForce
            yy1 = yy + HR%Radius + 1.0_EB*HR%B
            JJ = FLOOR(CELLSJ(FLOOR((yy1-YS)*RDYINT))+1.0_EB)
            I_OBST = OBST_INDEX_C(CELL_INDEX(II,JJ,KK))
            IF (SOLID(CELL_INDEX(II,JJ,KK)) .AND. .NOT.OBSTRUCTION(I_OBST)%HIDDEN) CYCLE CheckPPForce
         END IF
         !
         d_max = 1.0_EB*HR%B
         r_tmp(1) = HR%r_shoulder ! right circle
         r_tmp(2) = HR%r_torso     ! center circle
         r_tmp(3) = HR%r_shoulder ! left circle
         y_tmp(1) = yy - COS(HR%angle)*HR%d_shoulder ! right
         x_tmp(1) = xx + SIN(HR%angle)*HR%d_shoulder
         y_tmp(2) = yy ! torso
         x_tmp(2) = xx
         y_tmp(3) = yy + COS(HR%angle)*HR%d_shoulder ! left
         x_tmp(3) = xx - SIN(HR%angle)*HR%d_shoulder
         P2PLoop: DO ie = 1, n_humans
            HRE=>HUMAN(IE)
            r_tmp(4) = HRE%r_shoulder ! right circle
            r_tmp(5) = HRE%r_torso     ! center circle
            r_tmp(6) = HRE%r_shoulder ! left circle
            y_tmp(4) = HRE%Y - COS(HRE%angle)*HRE%d_shoulder ! right circle
            x_tmp(4) = HRE%X + SIN(HRE%angle)*HRE%d_shoulder
            y_tmp(5) = HRE%Y ! center circle
            x_tmp(5) = HRE%X
            y_tmp(6) = HRE%Y + COS(HRE%angle)*HRE%d_shoulder ! left circle
            x_tmp(6) = HRE%X - SIN(HRE%angle)*HRE%d_shoulder
            DO iii = 1, 3
               DO jjj = 4, 6
                  DIST = SQRT((x_tmp(jjj)-x_tmp(iii))**2 + (y_tmp(jjj)-y_tmp(iii))**2) - (r_tmp(jjj)+r_tmp(iii))
                  IF (DIST < d_max) CYCLE CheckPPForce
               END DO
            END DO
         END DO P2PLoop
         istat = 0
         EXIT CheckPPForce
      END DO CheckPPForce

      IF (istat == 0 ) THEN
         N_HUMANS = N_HUMANS + 1
         PNX%T_last = Tin
         PNX%ICOUNT = PNX%ICOUNT + 1
         IF (PNX%T_first <= T_BEGIN) PNX%T_first = Tin
         HR%X = xx
         HR%Y = yy
         HR%Z = zz
         !
         ILABEL_last = ILABEL_last + 1
         HR%ILABEL = ILABEL_last
         HR%SHOW = .TRUE.
         HR%COLOR_INDEX = 1
         SELECT CASE (COLOR_METHOD)
         CASE (-1)
            HR%COLOR_INDEX = 1
         CASE (0)
            HR%COLOR_INDEX = PNX%Avatar_Color_Index
         CASE (1,2)
            HR%COLOR_INDEX = 1    ! lonely human
         CASE (3)
            HR%COLOR_INDEX = evac_person_classes(PNX%IPC)%Avatar_Color_Index
         CASE (4)
            HR%COLOR_INDEX = 1
         CASE (5)
            ! Correct color is put, where the flow fields are chosen.
            HR%COLOR_INDEX = 1
         CASE Default
            WRITE(MESSAGE,'(A,I3,A)') 'ERROR: ENTRY_HUMAN COLOR METHOD',COLOR_METHOD, ' is not defined'
            CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
         END SELECT
         HR%FFIELD_NAME = TRIM(PNX%GRID_NAME)
         HR%I_FFIELD    = 0
         Mesh2Loop: DO i = 1, NMESHES
            IF ( EVACUATION_ONLY(I) .AND. TRIM(HR%FFIELD_NAME) == TRIM(MESH_NAME(i)) ) THEN
               HR%I_FFIELD = i
               EXIT Mesh2Loop
            END IF
         END DO Mesh2Loop
         IF ( HR%I_FFIELD == 0 ) THEN
            WRITE(LU_ERR,'(A,A,A,A)') 'ERROR: ENTR line ',TRIM(PNX%ID), ' problem with flow field name, ', &
                 TRIM(PNX%GRID_NAME),' not found'
            WRITE(LU_EVACOUT,'(A,A,A,A)') 'ERROR: ENTR line ',TRIM(PNX%ID), ' problem with flow field name, ', &
                 TRIM(PNX%GRID_NAME),' not found'
            STOP_STATUS = EVACUATION_STOP ; RETURN
         END IF
         HR%IMESH       = PNX%IMESH
         HR%INODE       = PNX%TO_INODE
         HR%NODE_NAME   = TRIM(PNX%TO_NODE)
         HR%U = 0.0_EB ; HR%V= 0.0_EB ; HR%W = 0.0_EB
         HR%IOR       = 0
         HR%F_Y       = 0.0_EB
         HR%F_X       = 0.0_EB
         HR%v0_fac    = 1.0_EB
         HR%SumForces = 0.0_EB
         HR%SumForces2 = 0.0_EB
         HR%IntDose   = 0.0_EB
         HR%Eta       = 0.0_EB
         HR%Ksi       = 0.0_EB
         HR%NewRnd    = .TRUE.
         HR%STR_SUB_INDX = PNX%STR_SUB_INDX
         HR%SKIP_WALL_FORCE_IOR = 0
         WRITE (LU_EVACOUT,fmt='(a,i6,a,f8.2,a,3a)') ' EVAC: Person n:o', HR%ILABEL, ' inserted ', Tin, &
              ' s, entry ', TRIM(PNX%ID),' ffield ', TRIM(HR%FFIELD_NAME)
         WRITE (LU_EVACOUT,fmt='(a,a)') ' person     x       y       z    Tpre    Tdet  ', &
         ' dia    v0   tau   i_gr i_ff'
         WRITE (LU_EVACOUT,fmt='(i6,5f8.2,3f6.2,i6,i4,i4)') HR%ILABEL, &
              HR%X, HR%Y, HR%Z, HR%Tpre, HR%Tdet,2.0_EB*HR%Radius, &
              HR%Speed, HR%Tau, HR%GROUP_ID, HR%i_ffield, HR%COLOR_INDEX
      END IF
      !
    END SUBROUTINE ENTRY_HUMAN

    SUBROUTINE ENTRY_CROWBAR_HUMAN(I_entry, Tin, NM, istat)
      IMPLICIT NONE
      !
      ! Insert humans into the domain if needed.
      ! Crowbar project: use camera information
      !
      ! Now the camera file for the entry is reread every time.
      ! The file information is not read in and saved for later use.
      ! To minimize the file reads, also the next time stamp is read
      ! and saved for later use. Also the last read time stamp is
      ! saved, so one can figure out, when one should read once again.
      ! I.e., It Tin >= next_time_stamp
      ! Entry type variables: CB_TimeLastRead, CB_TimeNextRead
      !
      ! This routine is also called for times T < T_BEGIN, so the initialization
      ! could be done.
      !
      ! Passed variables
      REAL(EB), INTENT(IN) :: Tin  ! Current time in the internal evacuation clock
      INTEGER,  INTENT(IN) :: I_entry, NM ! Index of the ENTR, Index of the mesh
      INTEGER, INTENT(OUT) :: istat ! 0: Agents are generated, 1: No new agents
      !
      ! Local variables
      REAL :: RN_REAL
      REAL(EB) :: RN, x1, x2, y1, y2, z1, z2, d_max, dist, xx, yy, zz, xx1, yy1
      REAL(EB) :: CB_TIME, EVEL, CB_TIME_OLD, CB_TIME_NOW
      INTEGER :: I, J, II, JJ, KK, ior, irnmax, irn, ie, CB_LINE_TYPE, CB_N_AGENTS, CB_ID_CAMERA, &
           N_A_old, N_CBA_old, CB_N_AGENTS_PREV_DT
      LOGICAL :: FIRST_PASS, TITLE_WROTED, SOLID_OBSERVATION, CROWBAR_READ_IN_TMP
      REAL(EB), DIMENSION(6) :: y_tmp, x_tmp, r_tmp
      INTEGER, DIMENSION(:), ALLOCATABLE :: CB_I_AGENT, CB_I_AGENT_PREV_DT
      REAL(EB), DIMENSION(:), ALLOCATABLE :: CB_H_AGENT
      REAL(EB), DIMENSION(:,:), ALLOCATABLE :: CB_XYZ_AGENT,CB_V_AGENT, CB_XYZ_AGENT_PREV_DT

      TYPE (EVAC_ENTR_TYPE), POINTER :: PNX =>NULL()
      TYPE (MESH_TYPE), POINTER :: MFF =>NULL()
      TYPE (EVAC_PERS_TYPE), POINTER :: PCP =>NULL()
      TYPE (HUMAN_TYPE), POINTER :: HR=>NULL(), HRE =>NULL()

      PNX => EVAC_ENTRYS(I_entry)

      IF (.NOT.PNX%CROWBAR_READ_IN) RETURN
      IF (PNX%IMESH /= NM ) RETURN
      IF (PNX%ICOUNT >= PNX%Max_Humans) RETURN ! This is not needed for Crowbar entries
      IF (Tin < PNX%CB_TimeNextRead-TWO_EPSILON_EB) RETURN ! Check the time for Crowbar entries

      MFF => MESHES(NM)
      FIRST_PASS = .FALSE.
      CB_TIME = 0.0_EB
      CB_ID_CAMERA = 0

      OPEN (LU_EVAC_CB,file=TRIM(PNX%CROWBAR_INPUT_FILE),form='formatted', status='old')
      CB_N_AGENTS = 0
      READ (LU_EVAC_CB,*) CB_LINE_TYPE, CB_TIME, CB_N_AGENTS, CB_ID_CAMERA
      IF (CB_LINE_TYPE/=0 .OR. CB_N_AGENTS<0) THEN
         CLOSE (LU_EVAC_CB)
          WRITE(LU_ERR,'(A)') 'ERROR, Entry_Crowbar_Human: CROWBAR INPUT FILE READ ERROR 1'
         STOP_STATUS = EVACUATION_STOP ; RETURN
      END IF
      IF (PNX%CB_TimeLastRead > Tin) THEN
         ! This is the first call to this cb_entry
         FIRST_PASS = .TRUE.
         PNX%CB_TimeLastRead = Tin
         WRITE(LU_EVACOUT,'(A,A)')   ' FDS+Evac: Crowbar input entry file is: ',TRIM(PNX%CROWBAR_INPUT_FILE)
         WRITE(LU_EVACOUT,'(A,I16)') '                          camera ID is: ',CB_ID_CAMERA
      END IF
      CB_N_AGENTS_PREV_DT = CB_N_AGENTS

      ! DO WHILE (PNX%CB_TimeLastRead >= CB_TIME-TWO_EPSILON_EB)
      DO WHILE (PNX%CB_TimeLastRead >= CB_TIME+0.5_EB*CROWBAR_DT_READ)
         ! Skip the times that are already read in (and used), but not the last read
         DO I = 1, CB_N_AGENTS_PREV_DT
            CB_LINE_TYPE = -1
            READ (LU_EVAC_CB,*) CB_LINE_TYPE
            IF (CB_LINE_TYPE/=1) THEN
               WRITE(LU_ERR,'(A)') 'ERROR, Entry_Crowbar_Human: CROWBAR INPUT FILE READ ERROR 2'
               CLOSE (LU_EVAC_CB)
               STOP_STATUS = EVACUATION_STOP ; RETURN
            END IF
         END DO
         CB_N_AGENTS_PREV_DT = 0
         READ (LU_EVAC_CB,*) CB_LINE_TYPE, CB_TIME, CB_N_AGENTS_PREV_DT, CB_ID_CAMERA
         IF (CB_LINE_TYPE/=0 .OR. CB_N_AGENTS_PREV_DT<0) THEN
            WRITE(LU_ERR,'(A)') 'ERROR, Entry_Crowbar_Human: CROWBAR INPUT FILE READ ERROR 1'
            CLOSE (LU_EVAC_CB)
            STOP_STATUS = EVACUATION_STOP ; RETURN
         END IF
         CB_TIME_OLD = CB_TIME
      END DO

      IfFirstPass: IF (.NOT.FIRST_PASS) THEN
         ! Read the previous time information for a camera
         ALLOCATE(CB_I_AGENT_PREV_DT(CB_N_AGENTS_PREV_DT),STAT=IZERO)
         CALL ChkMemErr('INIT_EVACUATION','CB_I_AGENT_PREV_DT',IZERO)
         ALLOCATE(CB_XYZ_AGENT_PREV_DT(CB_N_AGENTS_PREV_DT,3),STAT=IZERO)
         CALL ChkMemErr('INIT_EVACUATION','CB_XYZ_AGENT_PREV_DT',IZERO)
         CB_I_AGENT_PREV_DT = 0
         CB_XYZ_AGENT_PREV_DT = 0.0_EB
         DO I = 1, CB_N_AGENTS_PREV_DT
            CB_LINE_TYPE = -1
            READ (LU_EVAC_CB,*) CB_LINE_TYPE, CB_I_AGENT_PREV_DT(I),&
                 CB_XYZ_AGENT_PREV_DT(I,1),CB_XYZ_AGENT_PREV_DT(I,2),CB_XYZ_AGENT_PREV_DT(I,3)
            IF (CB_LINE_TYPE/=1) THEN
               WRITE(LU_ERR,'(A)') 'ERROR, Entry_Crowbar_Human: CROWBAR INPUT FILE READ ERROR 2A'
               CLOSE (LU_EVAC_CB)
               STOP_STATUS = EVACUATION_STOP ; RETURN
            END IF
         END DO
         CB_XYZ_AGENT_PREV_DT = CB_XYZ_AGENT_PREV_DT/1000.0_EB ! mm to m transformation

         CB_N_AGENTS = 0
         READ (LU_EVAC_CB,*) CB_LINE_TYPE, CB_TIME, CB_N_AGENTS, CB_ID_CAMERA
         IF (CB_LINE_TYPE/=0 .OR. CB_N_AGENTS<0) THEN
            WRITE(LU_ERR,'(A)') 'ERROR, Entry_Crowbar_Human: CROWBAR INPUT FILE READ ERROR 1B'
            CLOSE (LU_EVAC_CB)
            STOP_STATUS = EVACUATION_STOP ; RETURN
         END IF
      END IF IfFirstPass

      ! Read the current time information for a camera
      ALLOCATE(CB_I_AGENT(CB_N_AGENTS),STAT=IZERO)
      CALL ChkMemErr('INIT_EVACUATION','CB_I_AGENT',IZERO)
      ALLOCATE(CB_XYZ_AGENT(CB_N_AGENTS,3),STAT=IZERO)
      CALL ChkMemErr('INIT_EVACUATION','CB_XYZ_AGENT',IZERO)
      ALLOCATE(CB_V_AGENT(CB_N_AGENTS,3),STAT=IZERO)
      CALL ChkMemErr('INIT_EVACUATION','CB_V_AGENT',IZERO)
      ALLOCATE(CB_H_AGENT(CB_N_AGENTS),STAT=IZERO)
      CALL ChkMemErr('INIT_EVACUATION','CB_H_AGENT',IZERO)
      CB_I_AGENT = 0
      CB_XYZ_AGENT = 0.0_EB; CB_H_AGENT = 0.0_EB; CB_V_AGENT = 0.0_EB
      DO I = 1, CB_N_AGENTS
         CB_LINE_TYPE = -1
         READ (LU_EVAC_CB,*) CB_LINE_TYPE, CB_I_AGENT(I),&
              CB_XYZ_AGENT(I,1),CB_XYZ_AGENT(I,2),CB_XYZ_AGENT(I,3),CB_H_AGENT(I),&
              CB_V_AGENT(I,1),CB_V_AGENT(I,2),CB_V_AGENT(I,3)
         IF (CB_LINE_TYPE/=1) THEN
            WRITE(LU_ERR,'(A)') 'ERROR, Entry_Crowbar_Human: CROWBAR INPUT FILE READ ERROR 2'
            CLOSE (LU_EVAC_CB)
            STOP_STATUS = EVACUATION_STOP ; RETURN
         END IF
      END DO
      CB_XYZ_AGENT = CB_XYZ_AGENT/1000.0_EB ! mm to m transformation
      CB_H_AGENT   = CB_H_AGENT  /1000.0_EB ! mm to m transformation
      CB_V_AGENT   = CB_V_AGENT  /1000.0_EB ! mm/s to m/s transformation
      CB_V_AGENT   = MIN(CB_V_AGENT,5.0_EB) ! speed componets [-5 m/s, 5 m/s] interval
      CB_V_AGENT   = MAX(CB_V_AGENT,-5.0_EB)

      ! Check and read the next time stamp on the input file
      PNX%CB_TimeNextRead = HUGE(PNX%CB_TimeNextRead) ! If the file ends, do not read anymore
      PNX%CB_TimeLastRead = CB_TIME
      CB_TIME_NOW = CB_TIME
      READ (LU_EVAC_CB,*,END=100) CB_LINE_TYPE, CB_TIME
      IF (CB_LINE_TYPE/=0) THEN
         WRITE(LU_ERR,'(A)') 'ERROR, Entry_Crowbar_Human: CROWBAR INPUT FILE READ ERROR 3'
         CLOSE (LU_EVAC_CB)
         STOP_STATUS = EVACUATION_STOP ; RETURN
      END IF
      PNX%CB_TimeNextRead = CB_TIME ! Next time stamp in the camera file
100   CONTINUE
      CLOSE(LU_EVAC_CB)

      ! The XB of the ENTR namelist:
      X1  = PNX%X1
      X2  = PNX%X2
      Y1  = PNX%Y1
      Y2  = PNX%Y2
      Z1  = PNX%Z1
      Z2  = PNX%Z2
      ior = PNX%IOR

      PCP => EVAC_PERSON_CLASSES(PNX%IPC)

      TITLE_WROTED = .FALSE.
      N_A_old   = PNX%CB_N_Agents  ! CB_N_(Crowbar)Agents is updated in the EVAC_MOVE_LOOP for this time step
      N_CBA_old = PNX%CB_N_CrowbarAgents

      CB_Agents: DO J = 1, CB_N_AGENTS
         istat = 1 ! Default is "no space found", i.e., do not place this agent
         SOLID_OBSERVATION = .TRUE.
         IF (.NOT.FIRST_PASS) THEN
            SOLID_OBSERVATION = .FALSE.
            OldTimeAgents: DO I = 1, CB_N_AGENTS_PREV_DT
               ! Check if the new agent ID in the camera is a solid observation, i.e.,
               ! it exists at two consecutive time steps (the current and the previous one)
               IF (CB_I_AGENT(J)==CB_I_AGENT_PREV_DT(I)) THEN
                  SOLID_OBSERVATION = .TRUE.
                  CB_V_AGENT(J,1) = (CB_XYZ_AGENT(J,1)-CB_XYZ_AGENT_PREV_DT(I,1)) / MAX(EVAC_DT_MIN,(CB_TIME_NOW-CB_TIME_OLD))
                  CB_V_AGENT(J,2) = (CB_XYZ_AGENT(J,2)-CB_XYZ_AGENT_PREV_DT(I,2)) / MAX(EVAC_DT_MIN,(CB_TIME_NOW-CB_TIME_OLD))
                  CB_V_AGENT(J,:)   = MIN(CB_V_AGENT(J,:),5.0_EB) ! speed componets [-5 m/s, 5 m/s] interval
                  CB_V_AGENT(J,:)   = MAX(CB_V_AGENT(J,:),-5.0_EB)
                  EXIT OldTimeAgents
               END IF
            END DO OldTimeAgents
         END IF

         P2PLoop2: DO I = 1, N_HUMANS
            HRE=>HUMAN(I)
            ! HR%IEL: Which camera, i.e, ENTR line, (iel<0 for entries, >0 for evac namelists)
            IF (HRE%CROWBAR_READ_IN .AND. HRE%ID_CB==CB_I_AGENT(J) .AND. HRE%IEL==-I_entry) THEN
               ! The current agent read in at this time step exists already for this camera.
               ! Update the current location and velocity information: (x,y), v0
               HRE%X_CB = CB_XYZ_AGENT(J,1)
               HRE%Y_CB = CB_XYZ_AGENT(J,2)
               HRE%X2_CB = HRE%X
               HRE%Y2_CB = HRE%Y
               HRE%U_CB = CB_V_AGENT(J,1)
               HRE%V_CB = CB_V_AGENT(J,2)
               HRE%CROWBAR_READ_IN   = .TRUE.
               HRE%CROWBAR_UPDATE_V0 = .TRUE.
               HRE%T_LastRead_CB = Tin
               CYCLE CB_Agents ! The read in agent is already present in the calculation
            END IF

            ! Check the actual and the last read in positions, close to either one is enough
            ! The last read in position should be a relatively recent one, dt is about 0.1 s
            dist = SQRT((HRE%X-CB_XYZ_AGENT(J,1))**2 + (HRE%Y-CB_XYZ_AGENT(J,2))**2)
            IF (HRE%T_LastRead_CB > T-MIN(1.0_EB,5.0_EB*CROWBAR_DT_READ)) THEN
               ! Just for agents that are read in (camera information), not for ordinaty agents.
               ! Ordinary agents have T_LastRead_CB < t_begin.
               dist = MIN(dist,SQRT((HRE%X_CB-CB_XYZ_AGENT(J,1))**2 + (HRE%Y_CB-CB_XYZ_AGENT(J,2))**2))
            END IF
            ! If read at this or previous time step then this other agent is not considered for this "new"
            ! read in agent position. I.e., this HRE agent has already its own read in information.
            CROWBAR_READ_IN_TMP = .FALSE.
            IF (HRE%CROWBAR_READ_IN .AND. HRE%T_LastRead_CB > T-1.8_EB*CROWBAR_DT_READ) CROWBAR_READ_IN_TMP = .TRUE.
            IF (dist < 1.5_EB*HRE%Radius .AND. .NOT.CROWBAR_READ_IN_TMP) THEN
               ! There is a "normal" agent close at the read in position => change this to a Crowbar agent.
               write(lu_err,*)'*** Crowbar time,old ID, new ID ',Tin,HRE%ID_CB,CB_I_AGENT(J)
               HRE%ID_CB = CB_I_AGENT(J)
               HRE%X_CB  = CB_XYZ_AGENT(J,1)
               HRE%Y_CB  = CB_XYZ_AGENT(J,2)
               HRE%X2_CB = HRE%X
               HRE%Y2_CB = HRE%Y
               HRE%U_CB  = CB_V_AGENT(J,1)
               HRE%V_CB  = CB_V_AGENT(J,2)
               HRE%CROWBAR_READ_IN   = .TRUE.
               HRE%CROWBAR_UPDATE_V0 = .TRUE.
               HRE%T_LastRead_CB     = Tin
               CYCLE CB_Agents ! The read in agent is allocated to an existing non-cb agent
            END IF

         END DO P2PLoop2
         IF (.NOT.SOLID_OBSERVATION) CYCLE CB_Agents ! Not yet a solid observation

         IF (N_HUMANS+1 > N_HUMANS_DIM) THEN
            ! Re-allocation is not yet checked.
            WRITE(LU_ERR,*) 'ERROR: Insert CB Humans: no re-allocation yet'
            STOP_STATUS = EVACUATION_STOP ; RETURN
            CALL RE_ALLOCATE_HUMANS(1,NM)
            HUMAN=>MESHES(NM)%HUMAN
         END IF
         HR  => MESHES(NM)%HUMAN(N_HUMANS+1)
         CALL CLASS_PROPERTIES(HR,PCP,-I_entry)
         IF (STOP_STATUS==SETUP_STOP) RETURN
         HR%Tpre = 0.0_EB   ! Start to move directly, ENTR-agents
         HR%Tdet = T_BEGIN
         HR%IPC  = PNX%IPC  ! Index of PERS namelist
         HR%IEL  = -I_entry ! Evac lines > 0, Entr lines < 0
         HR%GROUP_ID = 0
         HR%I_Target = 0    ! Default, no target door yet
         HR%I_Door_Mode = 0 ! Default, no target door yet
         HR%I_DoorAlgo = PNX%I_AGENT_TYPE

         ! Hawk-Dove game parameters (ongoing research/develepment work)
         CALL RANDOM_NUMBER(RN_REAL)
         RN = REAL(RN_REAL,EB)
         HR%F_FallDown = F_MIN_FALL + RN*(F_MAX_FALL - F_MIN_FALL)
         HR%T_FallenDown = HUGE(HR%T_FallenDown)
         HR%Angle_FallenDown = 0.0_EB
         HR%SizeFac_FallenDown = 0.0_EB
         HR%T_CheckFallDown = Tin

         ! Agent positions and velocities are read from camera file
         HR%X    = CB_XYZ_AGENT(J,1)
         HR%Y    = CB_XYZ_AGENT(J,2)
         HR%X_CB = CB_XYZ_AGENT(J,1)
         HR%Y_CB = CB_XYZ_AGENT(J,2)
         HR%X2_CB = HR%X
         HR%Y2_CB = HR%Y
         HR%U    = CB_V_AGENT(J,1)
         HR%V    = CB_V_AGENT(J,2)
         HR%U_CB = CB_V_AGENT(J,1)
         HR%V_CB = CB_V_AGENT(J,2)
         HR%ID_CB = CB_I_AGENT(J)
         HR%CROWBAR_READ_IN = .TRUE.
         HR%CROWBAR_UPDATE_V0 = .FALSE.
         xx = CB_XYZ_AGENT(J,1)
         yy = CB_XYZ_AGENT(J,2)

         ! Set the body angle according to the velocity that is read in.
         ! If velocity close to zero => use random angle
         EVEL = SQRT(HR%U_CB**2 + HR%V_CB**2)
         IF (EVEL > 0.1_EB) THEN
            IF (HR%V_CB >= 0.0_EB) THEN
               HR%ANGLE = ACOS(HR%U_CB/EVEL)
            ELSE
               HR%ANGLE = 2.0_EB*PI - ACOS(HR%U_CB/EVEL)
            END IF
         ELSE
            CALL RANDOM_NUMBER(RN_REAL)
            RN = REAL(RN_REAL,EB)
            HR%Angle = 2.0_EB*Pi*rn
         END IF

         DO WHILE (HR%Angle >= 2.0_EB*Pi)
            HR%Angle = HR%Angle - 2.0_EB*Pi
         END DO
         DO WHILE (HR%Angle < 0.0_EB)
            HR%Angle = HR%Angle + 2.0_EB*Pi
         END DO
         !
         IF (ABS(ior) == 1) irnmax = INT(PNX%Width*4.0_EB)
         IF (ABS(ior) == 2) irnmax = INT(PNX%Width*4.0_EB)
         IF (ABS(ior) == 3 .OR. ior == 0) irnmax = INT((x2-x1)*4.0_EB)*INT((y2-y1)*4.0_EB)
         irnmax = MAX(irnmax,5)
         !
         irn = 0
         CheckPPForce: DO WHILE (irn < irnmax)

            zz  = 0.5_EB*(z1+z2)
            II = FLOOR( CELLSI(FLOOR((xx-XS)*RDXINT)) + 1.0_EB )
            JJ = FLOOR( CELLSJ(FLOOR((yy-YS)*RDYINT)) + 1.0_EB )
            KK = 1

            irn = irn + 1 ! Add counter of placement trials for this agent

            I_OBST = OBST_INDEX_C(CELL_INDEX(II,JJ,KK))
            IF (SOLID(CELL_INDEX(II,JJ,KK)) .AND. .NOT.OBSTRUCTION(I_OBST)%HIDDEN) CYCLE CheckPPForce
            IF ( ABS(ior) == 2 ) THEN
               xx1 = xx - HR%Radius - 1.0_EB*HR%B
               II = FLOOR(CELLSI(FLOOR((xx1-XS)*RDXINT))+1.0_EB)
               I_OBST = OBST_INDEX_C(CELL_INDEX(II,JJ,KK))
               IF (SOLID(CELL_INDEX(II,JJ,KK)) .AND. .NOT.OBSTRUCTION(I_OBST)%HIDDEN) CYCLE CheckPPForce
               xx1 = xx + HR%Radius + 1.0_EB*HR%B
               II = FLOOR(CELLSI(FLOOR((xx1-XS)*RDXINT))+1.0_EB)
               I_OBST = OBST_INDEX_C(CELL_INDEX(II,JJ,KK))
               IF (SOLID(CELL_INDEX(II,JJ,KK)) .AND. .NOT.OBSTRUCTION(I_OBST)%HIDDEN) CYCLE CheckPPForce
            END IF
            IF ( ABS(ior) == 1 ) THEN
               yy1 = yy - HR%Radius - 1.0_EB*HR%B
               JJ = FLOOR(CELLSJ(FLOOR((yy1-YS)*RDYINT))+1.0_EB)
               I_OBST = OBST_INDEX_C(CELL_INDEX(II,JJ,KK))
               IF (SOLID(CELL_INDEX(II,JJ,KK)) .AND. .NOT.OBSTRUCTION(I_OBST)%HIDDEN) CYCLE CheckPPForce
               yy1 = yy + HR%Radius + 1.0_EB*HR%B
               JJ = FLOOR(CELLSJ(FLOOR((yy1-YS)*RDYINT))+1.0_EB)
               I_OBST = OBST_INDEX_C(CELL_INDEX(II,JJ,KK))
               IF (SOLID(CELL_INDEX(II,JJ,KK)) .AND. .NOT.OBSTRUCTION(I_OBST)%HIDDEN) CYCLE CheckPPForce
            END IF
            !
            d_max = 1.0_EB*HR%B
            r_tmp(1) = HR%r_shoulder ! right circle
            r_tmp(2) = HR%r_torso     ! center circle
            r_tmp(3) = HR%r_shoulder ! left circle
            y_tmp(1) = yy - COS(HR%angle)*HR%d_shoulder ! right
            x_tmp(1) = xx + SIN(HR%angle)*HR%d_shoulder
            y_tmp(2) = yy ! torso
            x_tmp(2) = xx
            y_tmp(3) = yy + COS(HR%angle)*HR%d_shoulder ! left
            x_tmp(3) = xx - SIN(HR%angle)*HR%d_shoulder
            P2PLoop: DO ie = 1, n_humans
               HRE=>HUMAN(IE)
               ! HR%IEL: Which camera, i.e, ENTR line, (iel<0 for entries, >0 for evac namelists)
!!$               IF (HRE%CROWBAR_READ_IN .AND. HRE%ID_CB==HR%ID_CB .AND. HRE%IEL==HR%IEL) THEN
!!$                  ! The current agent read in at this time step exists already for this camera.
!!$                  ! Update the current location and velocity information: (x,y), v0
!!$                  HRE%U_CB = (CB_XYZ_AGENT(J,1)-HRE%X_CB)/MAX(EVAC_DT_MIN,Tin-HRE%T_LastRead_CB)
!!$                  HRE%V_CB = (CB_XYZ_AGENT(J,2)-HRE%Y_CB)/MAX(EVAC_DT_MIN,Tin-HRE%T_LastRead_CB)
!!$
!!$                  HRE%X_CB = CB_XYZ_AGENT(J,1)
!!$                  HRE%Y_CB = CB_XYZ_AGENT(J,2)
!!$                  HRE%X2_CB = HRE%X
!!$                  HRE%Y2_CB = HRE%Y
!!$                  HRE%U_CB = CB_V_AGENT(J,1)
!!$                  HRE%V_CB = CB_V_AGENT(J,2)
!!$                  HRE%CROWBAR_UPDATE_V0 = .TRUE.
!!$                  HRE%T_LastRead_CB = Tin
!!$                  EXIT CheckPPForce
!!$               END IF
               r_tmp(4) = HRE%r_shoulder ! right circle
               r_tmp(5) = HRE%r_torso     ! center circle
               r_tmp(6) = HRE%r_shoulder ! left circle
               y_tmp(4) = HRE%Y - COS(HRE%angle)*HRE%d_shoulder ! right circle
               x_tmp(4) = HRE%X + SIN(HRE%angle)*HRE%d_shoulder
               y_tmp(5) = HRE%Y ! center circle
               x_tmp(5) = HRE%X
               y_tmp(6) = HRE%Y + COS(HRE%angle)*HRE%d_shoulder ! left circle
               x_tmp(6) = HRE%X - SIN(HRE%angle)*HRE%d_shoulder
               DO iii = 1, 3
                  DO jjj = 4, 6
                     DIST = SQRT((x_tmp(jjj)-x_tmp(iii))**2 + (y_tmp(jjj)-y_tmp(iii))**2) - (r_tmp(jjj)+r_tmp(iii))
                     IF (DIST < d_max) CYCLE CheckPPForce
                  END DO
               END DO
            END DO P2PLoop
            istat = 0
            EXIT CheckPPForce
         END DO CheckPPForce

         IF (istat == 0 ) THEN ! An empty space for the agent is found
            N_HUMANS = N_HUMANS + 1
            PNX%T_last = Tin
            PNX%ICOUNT = PNX%ICOUNT + 1
            IF (PNX%T_first <= T_BEGIN) PNX%T_first = Tin
            HR%X = xx
            HR%Y = yy
            HR%Z = zz
            !
            ILABEL_last = ILABEL_last + 1
            HR%ILABEL = ILABEL_last
            HR%SHOW = .TRUE.
            HR%COLOR_INDEX = 1
            SELECT CASE (COLOR_METHOD)
            CASE (-1)
               HR%COLOR_INDEX = 1
            CASE (0)
               HR%COLOR_INDEX = PNX%Avatar_Color_Index
            CASE (1,2)
               HR%COLOR_INDEX = 1    ! lonely human
            CASE (3)
               HR%COLOR_INDEX = evac_person_classes(PNX%IPC)%Avatar_Color_Index
            CASE (4)
               HR%COLOR_INDEX = 1
            CASE (5)
               ! Correct color is put, where the flow fields are chosen.
               HR%COLOR_INDEX = 1
            CASE Default
               WRITE(LU_ERR,'(A,I3,A)') 'ERROR: ENTRY_HUMAN COLOR METHOD',COLOR_METHOD, ' is not defined'
               STOP_STATUS = EVACUATION_STOP ; RETURN
            END SELECT
            HR%FFIELD_NAME = TRIM(PNX%GRID_NAME)
            HR%I_FFIELD    = 0
            Mesh2Loop: DO i = 1, NMESHES
               IF ( EVACUATION_ONLY(I) .AND. TRIM(HR%FFIELD_NAME) == TRIM(MESH_NAME(i)) ) THEN
                  HR%I_FFIELD = i
                  EXIT Mesh2Loop
               END IF
            END DO Mesh2Loop
            IF ( HR%I_FFIELD == 0 ) THEN
               WRITE(LU_ERR,'(A,A,A,A)') 'ERROR: ENTR line ',TRIM(PNX%ID), ' problem with flow field name, ', &
                    TRIM(PNX%GRID_NAME),' not found'
               STOP_STATUS = EVACUATION_STOP ; RETURN
            END IF
            HR%IMESH       = PNX%IMESH
            HR%INODE       = PNX%TO_INODE
            HR%NODE_NAME   = TRIM(PNX%TO_NODE)
            HR%UBAR = HR%U_CB ; HR%VBAR= HR%V_CB
            HR%U = 0.0_EB ; HR%V= 0.0_EB
            HR%Tdet = MAX(EVAC_DT_MIN,Tin,PNX%CB_TimeLastRead)   ! Start to move directly, ENTR-agents
            !HR%Tdet = MAX(EVAC_DT_MIN,PNX%CB_TimeLastRead)   ! Start to move directly, ENTR-agents
            HR%Tpre = 0.0_EB
            HR%T_LastRead_CB = PNX%CB_TimeLastRead
            HR%CROWBAR_UPDATE_V0 = .TRUE.
            HR%W          = 0.0_EB
            HR%IOR        = 0
            HR%F_Y        = 0.0_EB
            HR%F_X        = 0.0_EB
            HR%v0_fac     = 1.0_EB
            HR%SumForces  = 0.0_EB
            HR%SumForces2 = 0.0_EB
            HR%IntDose    = 0.0_EB
            HR%Eta        = 0.0_EB
            HR%Ksi        = 0.0_EB
            HR%NewRnd     = .TRUE.
            HR%STR_SUB_INDX = PNX%STR_SUB_INDX
            HR%SKIP_WALL_FORCE_IOR = 0
            IF (.NOT.TITLE_WROTED) WRITE(LU_EVACOUT,'(A)') &
                 'i_agent x_i     y_i     z_i   h_i     vx_i    vy_i    vz_i  dia    v0   tau '
            TITLE_WROTED = .TRUE.
            WRITE(LU_EVACOUT,'(I4,3F8.2,F6.2,3F8.2,3F6.2)')CB_I_AGENT(J),&
                 HR%X,HR%Y,HR%Z,CB_H_AGENT(J),&
                 HR%U_CB,HR%V_CB,0.0_EB, 2.0_EB*HR%Radius, HR%Speed, HR%Tau
         END IF
      END DO CB_Agents


      DEALLOCATE(CB_H_AGENT)
      DEALLOCATE(CB_V_AGENT)
      DEALLOCATE(CB_XYZ_AGENT)
      DEALLOCATE(CB_I_AGENT)
      IF (.NOT.FIRST_PASS) THEN
         DEALLOCATE(CB_XYZ_AGENT_PREV_DT)
         DEALLOCATE(CB_I_AGENT_PREV_DT)
      END IF
      !
    END SUBROUTINE ENTRY_CROWBAR_HUMAN

    LOGICAL FUNCTION Is_XY_Within_Bounds(P1x,P1y,P2x1,P2x2,P2y1,P2y2,xtol,ytol)
      IMPLICIT NONE
      !
      REAL(EB), INTENT(IN) :: P1x,P1y
      REAL(EB), INTENT(IN) :: P2x1,P2x2,P2y1,P2y2,xtol,ytol
      !Returns .TRUE. if P2 is within the bounds of P1 with tolerances.
      Is_XY_Within_Bounds = .FALSE.
      IF (P1x >= P2x1-xtol .AND. P1x <= P2x2+xtol .AND. P1y >= P2y1-ytol .AND. P1y <= P2y2+ytol) &
           Is_XY_Within_Bounds = .TRUE.
      RETURN
    END FUNCTION Is_XY_Within_Bounds

    SUBROUTINE Corner_Forces(x1, y1, x11, y11, p2p_dist_max, P2P_U, P2P_V, Social_F, &
         Contact_F, P2P_Torque, d_walls, x_tmp, y_tmp, r_tmp, u_tmp, v_tmp, istat, CONTACT_FX, CONTACT_FY)
      IMPLICIT NONE

      ! Corner point - agent social and contact forces and torques
      !
      ! Inputs:  x1,y1: coordinates of the agent
      !          x11,y11: coordiantes of the corner
      !          x/y_tmp: coordinates of the centres of the spheres (1-3 used)
      !                   1: right, 2: centre, 3: left circle
      !          r_tmp: radii of the spheres
      !          u/v_tmp: velocities of the centres of the spheres (1-3 used)
      !          p2p_dist_max: cutoff distance of the social force
      ! Outputs: istat: 0 no errors
      !          P2P_U/V: Forces x/y directions
      !          P2P_Torque: Torque
      !          Social_F: Radial social line force
      !          Contact_F: Radial contact line force
      !          d_walls: Shortest distance to walls
      !
      !
      ! Passed variables
      REAL(EB), INTENT(IN) :: x1, y1, x11, y11, p2p_dist_max
      REAL(EB), DIMENSION(6), INTENT(IN) :: x_tmp, y_tmp, r_tmp, u_tmp, v_tmp
      REAL(EB), INTENT(INOUT) :: P2P_U, P2P_V, Social_F, Contact_F, P2P_Torque, d_walls, CONTACT_FX, CONTACT_FY
      INTEGER, INTENT(OUT) :: istat
      !
      ! Local variables
      REAL(EB) :: dist, CosPhiFac, u, v, Fc_x, Fc_y, FricFac, k_fric, evel
      INTEGER :: iii

      istat = 0
      IF (I_Fric_sw >= 1 ) THEN
         FricFac = 0.0_EB
         k_fric = HR%Kappa
      ELSE
         FricFac = 1.0_EB
         k_fric = HR%Gamma
      END IF

      dist = SQRT((y11-y1)**2 + (x11-x1)**2)
      IF (dist-HR%Radius > p2p_dist_max) THEN
         istat = 1
         RETURN
      END IF
      u = u_tmp(2) ; v = v_tmp(2)

      IF ( (u**2+v**2) >= TWO_EPSILON_EB ) THEN
         CosPhiFac = ((x11-X1)*u + (y11-Y1)*v) / (dist*SQRT(u**2+v**2))
         CosPhiFac = LambdaW + 0.5_EB*(1.0_EB-LambdaW)*(1.0_EB+CosPhiFac)
      ELSE
         CosPhiFac = 1.0_EB
      END IF
      IF (L_DEAD .OR. L_FALLEN_DOWN) CosPhiFac = 0.0_EB
      P2P_U = P2P_U + (X1-x11)*A_Wall*CosPhiFac*EXP(-(dist-HR%Radius)/B_Wall) / dist
      P2P_V = P2P_V + (Y1-y11)*A_Wall*CosPhiFac*EXP(-(dist-HR%Radius)/B_Wall) / dist
      Social_F = Social_F + ABS(A_Wall*CosPhiFac*EXP(-(dist-HR%Radius)/B_Wall))

      DO iii = 1, 3
         dist = SQRT( (y11-y_tmp(iii))**2 + (x11-x_tmp(iii))**2 )

         ! Next is |vector1|*|vector2|
         evel = SQRT((x11-x_tmp(iii))**2+(y11-y_tmp(iii))**2)*SQRT(u_tmp(iii)**2+v_tmp(iii)**2)
         IF (evel >= TWO_EPSILON_EB) evel = ((x11-x_tmp(iii))*u_tmp(iii) + &
              (y11-y_tmp(iii))*v_tmp(iii)) / evel   ! cos theta (scal_prod/(lenght1*length2)
         IF (evel >= TWO_EPSILON_EB) THEN
            d_walls = MIN( (dist-r_tmp(iii))/evel, d_walls)
         ELSE
            d_walls = MIN( (dist-r_tmp(iii))/0.01_EB, d_walls)
         END IF

         ! d_walls = Min(dist-r_tmp(iii),d_walls)

         IF (dist <= r_tmp(iii)) THEN
            Fc_x = (x_tmp(iii)-x11)* 2.0_EB*HR%C_Young*(r_tmp(iii)-dist)/dist
            Fc_y = (y_tmp(iii)-y11)* 2.0_EB*HR%C_Young*(r_tmp(iii)-dist)/dist
            !Only radial contact forces, i.e., pressure calculation
            Contact_F = Contact_F + SQRT(Fc_x**2 + Fc_y**2)
            ! Tangential contact force:
            ! I_FRIC_SW >= 1 is the default and it is now the only one supported
               Fc_x = Fc_x - k_fric*( (1.0_EB-FricFac)*(r_tmp(iii)-dist)+FricFac ) *(y_tmp(iii)-y11)* &
                    ( (y_tmp(iii)-y11)*u_tmp(iii) - (x_tmp(iii)-x11)*v_tmp(iii) )/dist**2
               Fc_y = Fc_y + k_fric*( (1.0_EB-FricFac)*(r_tmp(iii)-dist)+FricFac ) *(x_tmp(iii)-x11)* &
                    ( (y_tmp(iii)-y11)*u_tmp(iii) - (x_tmp(iii)-x11)*v_tmp(iii) )/dist**2
            P2P_Torque = P2P_Torque + Fc_y*(x11-x_tmp(iii)) - Fc_x*(y11-y_tmp(iii))
            P2P_U = P2P_U + Fc_x
            P2P_V = P2P_V + Fc_y
            Contact_FX = Contact_FX + Fc_x
            Contact_FY = Contact_FY + Fc_y
         END IF
      END DO

    END SUBROUTINE Corner_Forces

    SUBROUTINE Door_Forces(x_tmp, y_tmp, r_tmp, u_tmp, v_tmp, p2p_dist_max, d_xy,&
         P2P_U, P2P_V, Social_F, Contact_F, P2P_Torque, FoundWall_xy, CONTACT_FX, CONTACT_FY)
      IMPLICIT NONE
      !
      ! This routine adds forces from the door posts. (VENTs with VEL>0)
      !
      ! Inputs:  FoundWall_xy(1-4): True for solid walls (False for outflow VENTs)
      !          x/y_tmp: coordinates of the centres of the spheres (1-3 used)
      !                   1: right, 2: centre, 3: left circle
      !          r_tmp: radii of the spheres
      !          u/v_tmp: velocities of the centres of the spheres (1-3 used)
      !          p2p_dist_max: cutoff distance of the social force
      ! Outputs:
      !          P2P_U/V: Forces x/y directions
      !          P2P_Torque: Torque
      !          Social_F: Radial social line force
      !          Contact_F: Radial contact line force
      !          d_walls: Shortest distance to walls
      !
      ! Passed variables
      REAL(EB), INTENT(IN) :: p2p_dist_max
      REAL(EB), DIMENSION(4), INTENT(IN) :: d_xy
      LOGICAL, DIMENSION(4), INTENT(IN) :: FoundWall_xy
      REAL(EB), DIMENSION(6), INTENT(IN) :: x_tmp, y_tmp, r_tmp, u_tmp, v_tmp
      REAL(EB), INTENT(INOUT) :: P2P_U, P2P_V, Social_F, Contact_F, P2P_Torque, CONTACT_FX, CONTACT_FY
      !
      ! Local variables
      REAL(EB) :: dist1, dist2

      ! Check if there are doors (vents with vel >0)
      DO ii = 1, N_VENT
         IF (ABS(VENTS(ii)%IOR)>2 .OR. SURFACE(VENTS(ii)%SURF_INDEX)%VEL<=0) CYCLE
         dist1 = SQRT((VENTS(ii)%x1-x_tmp(2))**2 + (VENTS(ii)%y1-y_tmp(2))**2) ! door - agent centre distance
         dist2 = SQRT((VENTS(ii)%x2-x_tmp(2))**2 + (VENTS(ii)%y2-y_tmp(2))**2)
         IF ( (MIN(dist1,dist2)-HR%Radius) > P2P_DIST_MAX) CYCLE


         IF (VENTS(ii)%IOR== -1 .AND. FoundWall_xy(2)) CYCLE
         IF (VENTS(ii)%IOR== +1 .AND. FoundWall_xy(1)) CYCLE
         IF (VENTS(ii)%IOR== -2 .AND. FoundWall_xy(4)) CYCLE
         IF (VENTS(ii)%IOR== +2 .AND. FoundWall_xy(3)) CYCLE
         SELECT CASE(VENTS(ii)%IOR)
         CASE(-1)  ! wall at +x direction
            IF (.NOT.FoundWall_xy(2) .AND. ABS(VENTS(ii)%x1-d_xy(2))<0.01_EB .AND. &
                 (VENTS(ii)%y1<y_tmp(2).AND.VENTS(ii)%y2>y_tmp(2))) THEN
               ! There is a outflow vent here
               x11 = VENTS(ii)%x1
               y11 = VENTS(ii)%y1
               CALL Corner_Forces(x1, y1, x11, y11, p2p_dist_max, P2P_U, P2P_V, Social_F, &
                    Contact_F, P2P_Torque, d_walls, x_tmp, y_tmp, r_tmp, u_tmp, v_tmp, istat, CONTACT_FX, CONTACT_FY)
               x11 = VENTS(ii)%x2
               y11 = VENTS(ii)%y2
               CALL Corner_Forces(x1, y1, x11, y11, p2p_dist_max, P2P_U, P2P_V, Social_F, &
                    Contact_F, P2P_Torque, d_walls, x_tmp, y_tmp, r_tmp, u_tmp, v_tmp, istat, CONTACT_FX, CONTACT_FY)
            END IF
         CASE(+1)  ! wall at -x direction
            IF (.NOT.FoundWall_xy(1) .AND. ABS(VENTS(ii)%x1-d_xy(1))<0.01_EB .AND. &
                 (VENTS(ii)%y1<y_tmp(2).AND.VENTS(ii)%y2>y_tmp(2))) THEN
               ! There is a outflow vent here
               x11 = VENTS(ii)%x1
               y11 = VENTS(ii)%y1
               CALL Corner_Forces(x1, y1, x11, y11, p2p_dist_max, P2P_U, P2P_V, Social_F, &
                    Contact_F, P2P_Torque, d_walls, x_tmp, y_tmp, r_tmp, u_tmp, v_tmp, istat, CONTACT_FX, CONTACT_FY)
               x11 = VENTS(ii)%x2
               y11 = VENTS(ii)%y2
               CALL Corner_Forces(x1, y1, x11, y11, p2p_dist_max, P2P_U, P2P_V, Social_F, &
                    Contact_F, P2P_Torque, d_walls, x_tmp, y_tmp, r_tmp, u_tmp, v_tmp, istat, CONTACT_FX, CONTACT_FY)
            END IF
         CASE(-2)  ! wall at +y direction
            IF (.NOT.FoundWall_xy(4) .AND. ABS(VENTS(ii)%y1-d_xy(4))<0.01_EB .AND. &
                 (VENTS(ii)%x1<x_tmp(2).AND.VENTS(ii)%x2>x_tmp(2))) THEN
               ! There is a outflow vent here
               x11 = VENTS(ii)%x1
               y11 = VENTS(ii)%y1
               CALL Corner_Forces(x1, y1, x11, y11, p2p_dist_max, P2P_U, P2P_V, Social_F, &
                    Contact_F, P2P_Torque, d_walls, x_tmp, y_tmp, r_tmp, u_tmp, v_tmp, istat, CONTACT_FX, CONTACT_FY)
               x11 = VENTS(ii)%x2
               y11 = VENTS(ii)%y2
               CALL Corner_Forces(x1, y1, x11, y11, p2p_dist_max, P2P_U, P2P_V, Social_F, &
                    Contact_F, P2P_Torque, d_walls, x_tmp, y_tmp, r_tmp, u_tmp, v_tmp, istat, CONTACT_FX, CONTACT_FY)
            END IF
         CASE(+2)  ! wall at -y direction
            IF (.NOT.FoundWall_xy(3) .AND. ABS(VENTS(ii)%y1-d_xy(3))<0.01_EB .AND. &
                 (VENTS(ii)%x1<x_tmp(2).AND.VENTS(ii)%x2>x_tmp(2))) THEN
               ! There is a outflow vent here
               x11 = VENTS(ii)%x1
               y11 = VENTS(ii)%y1
               CALL Corner_Forces(x1, y1, x11, y11, p2p_dist_max, P2P_U, P2P_V, Social_F, &
                    Contact_F, P2P_Torque, d_walls, x_tmp, y_tmp, r_tmp, u_tmp, v_tmp, istat, CONTACT_FX, CONTACT_FY)
               x11 = VENTS(ii)%x2
               y11 = VENTS(ii)%y2
               CALL Corner_Forces(x1, y1, x11, y11, p2p_dist_max, P2P_U, P2P_V, Social_F, &
                    Contact_F, P2P_Torque, d_walls, x_tmp, y_tmp, r_tmp, u_tmp, v_tmp, istat, CONTACT_FX, CONTACT_FY)
            END IF
         END SELECT
      END DO

    END SUBROUTINE Door_Forces

    SUBROUTINE Wall_SocialForces(x_tmp, y_tmp, r_tmp, p2p_dist_max, d_xy, P2P_U, P2P_V, Social_F, FoundWall_xy)
      IMPLICIT NONE
      !
      ! wall - agent social forces
      !
      ! Inputs:  FoundWall_xy(1-4): True for solid walls (False for outflow VENTs)
      !          x/y_tmp: coordinates of the centres of the spheres (1-3 used)
      !                   1: right, 2: centre, 3: left circle
      !          r_tmp: radii of the spheres
      !          p2p_dist_max: cutoff distance of the social force
      ! Outputs:
      !          P2P_U/V: Forces x/y directions
      !          Social_F: Radial social line force
      !
      !
      ! Passed variables
      REAL(EB), INTENT(IN) :: p2p_dist_max
      REAL(EB), DIMENSION(6), INTENT(IN) :: x_tmp, y_tmp, r_tmp
      REAL(EB), DIMENSION(4), INTENT(IN) :: d_xy
      LOGICAL, DIMENSION(4), INTENT(IN) :: FoundWall_xy
      REAL(EB), INTENT(INOUT) :: P2P_U, P2P_V, Social_F
      !
      ! Local variables
      INTEGER :: is, idir, iii
      REAL(EB) :: CosPhiFac, dist, F_soc, F_tmp

      ! -x direction
      is   = -1
      idir =  1
      F_soc = 0.0_EB
      DO iii = 1,3
         dist = ABS(d_xy(idir) - x_tmp(iii)) ! wall - agent centre distance
         IF (dist-r_tmp(iii) <= P2P_DIST_MAX .AND. FoundWall_xy(idir)) THEN
            IF ( (HR%U**2+HR%V**2) >= TWO_EPSILON_EB ) THEN
               CosPhiFac = (is*HR%U)/SQRT(HR%U**2+HR%V**2)
               CosPhiFac = LambdaW + 0.5_EB*(1.0_EB-LambdaW)*(1.0_EB+CosPhiFac)
            ELSE
               CosPhiFac = 1.0_EB
            END IF
            F_tmp = - is*A_Wall*CosPhiFac*EXP( -(dist-r_tmp(iii))/B_Wall )
            IF (ABS(F_tmp) > ABS(F_soc)) F_soc = F_tmp
         END IF
      END DO
      P2P_U = P2P_U + F_soc
      Social_F = Social_F + ABS(F_soc)

      is   = +1
      idir =  2
      F_soc = 0.0_EB
      DO iii = 1,3
         dist = ABS(d_xy(idir) - x_tmp(iii)) ! wall - agent centre distance
         IF (dist-r_tmp(iii) <= P2P_DIST_MAX .AND. FoundWall_xy(idir)) THEN
            IF ( (HR%U**2+HR%V**2) >= TWO_EPSILON_EB ) THEN
               CosPhiFac = (is*HR%U)/SQRT(HR%U**2+HR%V**2)
               CosPhiFac = LambdaW + 0.5_EB*(1.0_EB-LambdaW)*(1.0_EB+CosPhiFac)
            ELSE
               CosPhiFac = 1.0_EB
            END IF
            F_tmp = - is*A_Wall*CosPhiFac*EXP( -(dist-r_tmp(iii))/B_Wall )
            IF (ABS(F_tmp) > ABS(F_soc)) F_soc = F_tmp
         END IF
      END DO
      P2P_U = P2P_U + F_soc
      Social_F = Social_F + ABS(F_soc)

      is   = -1
      idir =  3
      F_soc = 0.0_EB
      DO iii = 1,3
         dist = ABS(d_xy(idir) - y_tmp(iii)) ! wall - agent centre distance
         IF (dist-r_tmp(iii) <= P2P_DIST_MAX .AND. FoundWall_xy(idir)) THEN
            IF ( (HR%U**2+HR%V**2) >= TWO_EPSILON_EB ) THEN
               CosPhiFac = (is*HR%V)/SQRT(HR%U**2+HR%V**2)
               CosPhiFac = LambdaW + 0.5_EB*(1.0_EB-LambdaW)*(1.0_EB+CosPhiFac)
            ELSE
               CosPhiFac = 1.0_EB
            END IF
            F_tmp = - is*A_Wall*CosPhiFac*EXP( -(dist-r_tmp(iii))/B_Wall )
            IF (ABS(F_tmp) > ABS(F_soc)) F_soc = F_tmp
         END IF
      END DO
      P2P_V = P2P_V + F_soc
      Social_F = Social_F + ABS(F_soc)

      is   = +1
      idir =  4
      F_soc = 0.0_EB
      DO iii = 1,3
         dist = ABS(d_xy(idir) - y_tmp(iii)) ! wall - agent centre distance
         IF (dist-r_tmp(iii) <= P2P_DIST_MAX .AND. FoundWall_xy(idir)) THEN
            IF ( (HR%U**2+HR%V**2) >= TWO_EPSILON_EB ) THEN
               CosPhiFac = (is*HR%V)/SQRT(HR%U**2+HR%V**2)
               CosPhiFac = LambdaW + 0.5_EB*(1.0_EB-LambdaW)*(1.0_EB+CosPhiFac)
            ELSE
               CosPhiFac = 1.0_EB
            END IF
            F_tmp = - is*A_Wall*CosPhiFac*EXP( -(dist-r_tmp(iii))/B_Wall )
            IF (ABS(F_tmp) > ABS(F_soc)) F_soc = F_tmp
         END IF
      END DO
      P2P_V = P2P_V + F_soc
      Social_F = Social_F + ABS(F_soc)

    END SUBROUTINE Wall_SocialForces

    SUBROUTINE Wall_ContactForces(x_tmp, y_tmp, r_tmp, u_tmp, v_tmp, d_xy,&
         P2P_U, P2P_V, P2P_Torque, Contact_F, d_walls, FoundWall_xy, CONTACT_FX, CONTACT_FY)
      IMPLICIT NONE
      !
      ! wall - agent contact forces
      !
      ! Inputs:  FoundWall_xy(1-4): True for solid walls (False for outflow VENTs)
      !          x/y_tmp: coordinates of the centres of the spheres (1-3 used)
      !                   1: right, 2: centre, 3: left circle
      !          r_tmp: radii of the spheres
      !          u/v_tmp: velocities of the centres of the spheres (1-3 used)
      !          p2p_dist_max: cutoff distance of the social force
      ! Outputs:
      !          P2P_U/V: Forces x/y directions
      !          P2P_Torque: Torque
      !          Social_F: Radial social line force
      !          Contact_F: Radial contact line force
      !          d_walls: Shortest distance to walls
      !
      ! Passed variables
      REAL(EB), INTENT(IN) :: x_tmp, y_tmp, r_tmp, u_tmp, v_tmp
      REAL(EB), DIMENSION(4), INTENT(IN) :: d_xy
      LOGICAL, DIMENSION(4), INTENT(IN) :: FoundWall_xy
      REAL(EB), INTENT(INOUT) :: P2P_U, P2P_V, P2P_Torque, Contact_F, d_walls, CONTACT_FX, CONTACT_FY
      !
      ! Local variables
      INTEGER :: is, idir
      REAL(EB) :: Fc_y, Fc_x, dist, evel

      ! -x direction
      is   = -1
      idir =  1
      dist = ABS(d_xy(idir) - x_tmp)

      IF (FoundWall_xy(idir)) THEN
         ! Next is |vector1|*|vector2|
         evel = ABS(d_xy(idir)-x_tmp)*SQRT(u_tmp**2+v_tmp**2)
         ! Next is cos(theta) = (scal_prod/(lenght1*length2)
         IF (evel >= TWO_EPSILON_EB) evel = ((d_xy(idir)-x_tmp)*u_tmp) / evel
         IF (evel > 0.01_EB) THEN
            d_walls = MIN( (dist-r_tmp)/evel, d_walls)
         ELSE
            d_walls = MIN( (dist-r_tmp)/0.01_EB, d_walls)
         END IF
      END IF

      IF (dist-r_tmp <= 0.0_EB .AND. FoundWall_xy(idir)) THEN
         Fc_x = -is*2.0_EB*HR%C_Young*(r_tmp-dist)
         Fc_y = 0.0_EB
         Contact_F = Contact_F + ABS(Fc_x)
         IF (I_Fric_sw >= 1 ) THEN
            Fc_y = Fc_y - HR%Kappa*(r_tmp-dist)*v_tmp
         ELSE
            Fc_y = Fc_y - HR%Gamma*v_tmp
         END IF
         Contact_FX = Contact_FX + Fc_x
         Contact_FY = Contact_FY + Fc_y
         P2P_Torque = P2P_Torque + Fc_y*(d_xy(idir)-HR%X) - Fc_x*(y_tmp-HR%Y)
         P2P_U = P2P_U + Fc_x
         P2P_V = P2P_V + Fc_y
      END IF

      ! +x direction
      is   = +1
      idir =  2
      dist = ABS(d_xy(idir) - x_tmp)

      IF (FoundWall_xy(idir)) THEN
         ! Next is |vector1|*|vector2|
         evel = ABS(d_xy(idir)-x_tmp)*SQRT(u_tmp**2+v_tmp**2)
         ! Next is cos(theta) = (scal_prod/(lenght1*length2)
         IF (evel >= TWO_EPSILON_EB) evel = ((d_xy(idir)-x_tmp)*u_tmp) / evel
         IF (evel > 0.01_EB) THEN
            d_walls = MIN( (dist-r_tmp)/evel, d_walls)
         ELSE
            d_walls = MIN( (dist-r_tmp)/0.01_EB, d_walls)
         END IF
      END IF

      IF (dist-r_tmp <= 0.0_EB .AND. FoundWall_xy(idir)) THEN
         Fc_x = -is*2.0_EB*HR%C_Young*(r_tmp-dist)
         Fc_y = 0.0_EB
         Contact_F = Contact_F + ABS(Fc_x)
         IF (I_Fric_sw >= 1 ) THEN
            Fc_y = Fc_y - HR%Kappa*(r_tmp-dist)*v_tmp
         ELSE
            Fc_y = Fc_y - HR%Gamma*v_tmp
         END IF
         Contact_FX = Contact_FX + Fc_x
         Contact_FY = Contact_FY + Fc_y
         P2P_Torque = P2P_Torque + Fc_y*(d_xy(idir)-HR%X) - Fc_x*(y_tmp-HR%Y)
         P2P_U = P2P_U + Fc_x
         P2P_V = P2P_V + Fc_y
      END IF

      ! -y direction
      is   = -1
      idir =  3
      dist = ABS(d_xy(idir) - y_tmp)

      IF (FoundWall_xy(idir)) THEN
         ! Next is |vector1|*|vector2|
         evel = ABS(d_xy(idir)-y_tmp)*SQRT(u_tmp**2+v_tmp**2)
         ! Next is cos(theta) = (scal_prod/(lenght1*length2)
         IF (evel >= TWO_EPSILON_EB) evel = ((d_xy(idir)-y_tmp)*v_tmp) / evel
         IF (evel > 0.01_EB) THEN
            d_walls = MIN( (dist-r_tmp)/evel, d_walls)
         ELSE
            d_walls = MIN( (dist-r_tmp)/0.01_EB, d_walls)
         END IF
      END IF

      IF (dist-r_tmp <= 0.0_EB .AND. FoundWall_xy(idir)) THEN
         Fc_y = -is*2.0_EB*HR%C_Young*(r_tmp-dist)
         Fc_x = 0.0_EB
         Contact_F = Contact_F + ABS(Fc_y)
         IF (I_Fric_sw >= 1 ) THEN
            Fc_x = Fc_x - HR%Kappa*(r_tmp-dist)*u_tmp
         ELSE
            Fc_x = Fc_x - HR%Gamma*u_tmp
         END IF
         Contact_FX = Contact_FX + Fc_x
         Contact_FY = Contact_FY + Fc_y
         P2P_Torque = P2P_Torque + Fc_y*(x_tmp-HR%X) - Fc_x*(d_xy(idir)-HR%Y)
         P2P_U = P2P_U + Fc_x
         P2P_V = P2P_V + Fc_y
      END IF

      ! +y direction
      is   = +1
      idir =  4
      dist = ABS(d_xy(idir) - y_tmp)

      IF (FoundWall_xy(idir)) THEN
         ! Next is |vector1|*|vector2|
         evel = ABS(d_xy(idir)-y_tmp)*SQRT(u_tmp**2+v_tmp**2)
         ! Next is cos(theta) = (scal_prod/(lenght1*length2)
         IF (evel >= TWO_EPSILON_EB) evel = ((d_xy(idir)-y_tmp)*v_tmp) / evel
         IF (evel > 0.01_EB) THEN
            d_walls = MIN( (dist-r_tmp)/evel, d_walls)
         ELSE
            d_walls = MIN( (dist-r_tmp)/0.01_EB, d_walls)
         END IF
      END IF

      IF (dist-r_tmp <= 0.0_EB .AND. FoundWall_xy(idir)) THEN
         Fc_y = -is*2.0_EB*HR%C_Young*(r_tmp-dist)
         Fc_x = 0.0_EB
         Contact_F = Contact_F + ABS(Fc_y)
         IF (I_Fric_sw >= 1 ) THEN
            Fc_x = Fc_x - HR%Kappa*(r_tmp-dist)*u_tmp
         ELSE
            Fc_x = Fc_x - HR%Gamma*u_tmp
         END IF
         Contact_FX = Contact_FX + Fc_x
         Contact_FY = Contact_FY + Fc_y
         P2P_Torque = P2P_Torque + Fc_y*(x_tmp-HR%X) - Fc_x*(d_xy(idir)-HR%Y)
         P2P_U = P2P_U + Fc_x
         P2P_V = P2P_V + Fc_y
      END IF

    END SUBROUTINE Wall_ContactForces
    !
  END SUBROUTINE EVACUATE_HUMANS
! ============================================================
! EVACUATE_HUMANS ENDS HERE.
! ============================================================
!
! ============================================================
! NEXT ARE EVAC MODULE SUBPROGRMAS
! ============================================================
!
  SUBROUTINE CLASS_PROPERTIES(HR,PCP,IEL)
    IMPLICIT NONE
    !
    ! Passed variables
    INTEGER :: IEL  ! >0: EVAC line index, <0: ENTR line index=Abs(iel)
    TYPE (HUMAN_TYPE), POINTER:: HR
    TYPE (EVAC_PERS_TYPE), POINTER:: PCP
    !
    ! Local variables
    ! How many rnd numbers per one call to the rnd routine
    INTEGER, PARAMETER  :: n_rnd=1, n_max_par=4
    INTEGER  :: n_par, RandomType, I_det_dist, I_pre_dist
    REAL(EB) :: Tdet_low, Tdet_high, Tdet_mean, Tdet_para, Tdet_para2
    REAL(EB) :: Tpre_low, Tpre_high, Tpre_mean, Tpre_para, Tpre_para2
    ! No more than 4 numbers needed to specify the distributions
    REAL(EB) :: RandomPara(n_max_par)
    REAL(EB) :: rnd_vec(n_rnd)
    ! 1: uniform (TESTED: OK)
    ! 2: Truncated normal (TESTED: OK)
    ! 3: gamma  (TESTED: OK)
    ! 4: normal (TESTED: OK)
    ! 5: lognormal (TESTED: OK)
    ! 6: beta (TESTED: OK)
    ! 7: Triangular (TESTED: OK)
    ! 8: Weibull (TESTED: OK) (alpha=1: Exponential)
    ! 9: Gumbel (TESTED: OK)

    I_det_dist = PCP%I_DET_DIST
    Tdet_low   = PCP%Tdet_low
    Tdet_high  = PCP%Tdet_high
    Tdet_mean  = PCP%Tdet_mean
    Tdet_para  = PCP%Tdet_para
    Tdet_para2 = PCP%Tdet_para2
    I_pre_dist = PCP%I_PRE_DIST
    Tpre_low   = PCP%Tpre_low
    Tpre_high  = PCP%Tpre_high
    Tpre_mean  = PCP%Tpre_mean
    Tpre_para  = PCP%Tpre_para
    Tpre_para2 = PCP%Tpre_para2
    IF (IEL>0) THEN
       IF (EVAC_EVACS(IEL)%I_DET_DIST>-1) THEN
          I_det_dist = EVAC_EVACS(IEL)%I_DET_DIST
          Tdet_low   = EVAC_EVACS(IEL)%Tdet_low
          Tdet_high  = EVAC_EVACS(IEL)%Tdet_high
          Tdet_mean  = EVAC_EVACS(IEL)%Tdet_mean
          Tdet_para  = EVAC_EVACS(IEL)%Tdet_para
          Tdet_para2 = EVAC_EVACS(IEL)%Tdet_para2
       END IF
       IF (EVAC_EVACS(IEL)%I_PRE_DIST>-1) THEN
          I_pre_dist = EVAC_EVACS(IEL)%I_PRE_DIST
          Tpre_low   = EVAC_EVACS(IEL)%Tpre_low
          Tpre_high  = EVAC_EVACS(IEL)%Tpre_high
          Tpre_mean  = EVAC_EVACS(IEL)%Tpre_mean
          Tpre_para  = EVAC_EVACS(IEL)%Tpre_para
          Tpre_para2 = EVAC_EVACS(IEL)%Tpre_para2
       END IF
    END IF

    SELECT CASE(PCP%I_VEL_DIST)
    CASE(-1)
       CALL SHUTDOWN('ERROR: Class_Properties: -1',PROCESS_0_ONLY=.FALSE.) ; RETURN
    CASE(0)
       HR%Speed  = PCP%V_mean
    CASE(1)   ! Uniform
       ! Parameters: (ave,min,max) ave not used
       n_par = 3
       Randomtype = 1
       RandomPara(1) = 0.5_EB*(PCP%V_high+PCP%V_low)
       RandomPara(2) = PCP%V_low
       RandomPara(3) = PCP%V_high
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Speed = rnd_vec(1)
    CASE(2)   ! Truncated Normal
       ! Parameters: (ave,sigma,min,max)
       n_par = 4
       Randomtype = 2
       RandomPara(1) = PCP%V_mean
       RandomPara(2) = PCP%V_para
       RandomPara(3) = PCP%V_low
       RandomPara(4) = PCP%V_high
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Speed = rnd_vec(1)
    CASE(3)   ! Gamma
       ! Parameters: (ave,alpha,beta) ave not used
       n_par = 3
       Randomtype = 3
       RandomPara(1) = PCP%V_mean
       RandomPara(2) = PCP%V_para
       RandomPara(3) = PCP%V_para2
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Speed = rnd_vec(1)
    CASE(4)   ! Normal
       ! Parameters: (ave,sigma)
       n_par = 2
       Randomtype = 4
       RandomPara(1) = PCP%V_mean
       RandomPara(2) = PCP%V_para
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Speed = rnd_vec(1)
    CASE(5)   ! LogNormal
       ! mean and variance of log(x) should be given
       ! Parameters: (ave,sigma) of ln(x)
       n_par = 4
       Randomtype = 5
       RandomPara(1) = PCP%V_mean  ! mu of ln(x)
       RandomPara(2) = PCP%V_para  ! sigma of ln(x)
       RandomPara(3) = PCP%V_high  ! high end cutoff of x
       RandomPara(4) = PCP%V_para2 ! shift of x
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Speed = rnd_vec(1)
    CASE(6)   ! Beta
       ! Parameters: (ave,a,b) ave not used
       n_par = 3
       Randomtype = 6
       RandomPara(1) = PCP%V_mean
       RandomPara(2) = PCP%V_para
       RandomPara(3) = PCP%V_para2
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Speed = rnd_vec(1)
    CASE(7)   ! Triangular
       ! Parameters: (peak,min,max)
       n_par = 3
       Randomtype = 7
       RandomPara(1) = PCP%V_mean
       RandomPara(2) = PCP%V_low
       RandomPara(3) = PCP%V_high
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Speed = rnd_vec(1)
    CASE(8)   ! Weibull  (alpha=1: Exponential)
       ! Parameters: (ave,alpha,lambda) ave not used
       n_par = 3
       Randomtype = 8
       RandomPara(1) = PCP%V_mean
       RandomPara(2) = PCP%V_para
       RandomPara(3) = PCP%V_para2
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Speed = rnd_vec(1)
    CASE(9)   ! Gumbel
       ! Parameters: (ave,alpha) ave not used
       n_par = 2
       Randomtype = 9
       RandomPara(1) = PCP%V_mean
       RandomPara(2) = PCP%V_para
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Speed = rnd_vec(1)
    CASE Default
       CALL SHUTDOWN('ERROR: Class_Properties I_VEL_DIST',PROCESS_0_ONLY=.FALSE.) ; RETURN
    END SELECT
    HR%Speed = Max(HR%Speed, 0.0_EB)
    
    SELECT CASE(PCP%I_DIA_DIST)
    CASE(-1)
       CALL SHUTDOWN('ERROR: Class_Properties: -1',PROCESS_0_ONLY=.FALSE.) ; RETURN
    CASE(0)
       HR%Radius  = 0.5_EB*PCP%D_mean
    CASE(1)   ! Uniform
       ! Parameters: (ave,min,max) ave not used
       n_par = 3
       Randomtype = 1
       RandomPara(1) = PCP%D_mean
       RandomPara(2) = PCP%D_low
       RandomPara(3) = PCP%D_high
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Radius = 0.5_EB*rnd_vec(1)
    CASE(2)   ! Truncated Normal
       n_par = 4
       Randomtype = 2
       RandomPara(1) = PCP%D_mean
       RandomPara(2) = PCP%D_para
       RandomPara(3) = PCP%D_low
       RandomPara(4) = PCP%D_high
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Radius = 0.5_EB*rnd_vec(1)
    CASE(3)   ! Gamma
       ! Parameters: (ave,alpha,beta) ave not used
       n_par = 3
       Randomtype = 3
       RandomPara(1) = PCP%D_mean
       RandomPara(2) = PCP%D_para
       RandomPara(3) = PCP%D_para2
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Radius = 0.5_EB*rnd_vec(1)
    CASE(4)   ! Normal
       ! Parameters: (ave,sigma)
       n_par = 2
       Randomtype = 4
       RandomPara(1) = PCP%D_mean
       RandomPara(2) = PCP%D_para
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Radius = 0.5_EB*rnd_vec(1)
    CASE(5)   ! LogNormal
       ! mean and variance of log(x) should be given
       ! Parameters: (ave,sigma) of ln(x)
       n_par = 4
       Randomtype = 5
       RandomPara(1) = PCP%D_mean  ! mu of ln(x)
       RandomPara(2) = PCP%D_para  ! sigma of ln(x)
       RandomPara(3) = PCP%D_high  ! high end cutoff
       RandomPara(4) = PCP%D_para2 ! shift
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Radius = 0.5_EB*rnd_vec(1)
    CASE(6)   ! Beta
       ! Parameters: (ave,a,b) ave not used
       n_par = 3
       Randomtype = 6
       RandomPara(1) = PCP%D_mean
       RandomPara(2) = PCP%D_para
       RandomPara(3) = PCP%D_para2
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Radius = 0.5_EB*rnd_vec(1)
    CASE(7)   ! Triangular
       ! Parameters: (peak,min,max)
       n_par = 3
       Randomtype = 7
       RandomPara(1) = PCP%D_mean
       RandomPara(2) = PCP%D_low
       RandomPara(3) = PCP%D_high
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Radius = 0.5_EB*rnd_vec(1)
    CASE(8)   ! Weibull  (alpha=1: Exponential)
       ! Parameters: (ave,alpha,lambda)
       n_par = 3
       Randomtype = 8
       RandomPara(1) = PCP%D_mean
       RandomPara(2) = PCP%D_para
       RandomPara(3) = PCP%D_para2
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Radius = 0.5_EB*rnd_vec(1)
    CASE(9)   ! Gumbel
       ! Parameters: (ave,alpha)
       n_par = 2
       Randomtype = 9
       RandomPara(1) = PCP%D_mean
       RandomPara(2) = PCP%D_para
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Radius = 0.5_EB*rnd_vec(1)
    CASE Default
       CALL SHUTDOWN('ERROR: Class_Properties I_DIA_DIST',PROCESS_0_ONLY=.FALSE.) ; RETURN
    END SELECT
    HR%Radius = Max(HR%Radius, 0.05_EB) ! 5cm minimum radius
    HR%Mass   = PCP%m_agent*(HR%Radius/0.27_EB)**2

    SELECT CASE(PCP%I_TAU_DIST)
    CASE(-1)
       CALL SHUTDOWN('ERROR: Class_Properties: -1',PROCESS_0_ONLY=.FALSE.) ; RETURN
    CASE(0)
       HR%Tau  = PCP%Tau_mean
    CASE(1)   ! Uniform
       ! Parameters: (ave,min,max) ave not used
       n_par = 3
       Randomtype = 1
       RandomPara(1) = 0.5_EB*(PCP%Tau_high+PCP%Tau_low)
       PCP%Tau_mean = 0.5_EB*(PCP%Tau_high+PCP%Tau_low)
       RandomPara(2) = PCP%Tau_low
       RandomPara(3) = PCP%Tau_high
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tau = rnd_vec(1)
    CASE(2)   ! Truncated Normal
       ! Parameters: (ave,sigma,min,max)
       n_par = 4
       Randomtype = 2
       RandomPara(1) = PCP%Tau_mean
       RandomPara(2) = PCP%Tau_para
       RandomPara(3) = PCP%Tau_low
       RandomPara(4) = PCP%Tau_high
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tau = rnd_vec(1)
    CASE(3)   ! Gamma
       ! Parameters: (ave,alpha,beta) ave not used
       n_par = 3
       Randomtype = 3
       RandomPara(1) = PCP%Tau_mean
       RandomPara(2) = PCP%Tau_para
       RandomPara(3) = PCP%Tau_para2
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tau = rnd_vec(1)
    CASE(4)   ! Normal
       ! Parameters: (ave,sigma)
       n_par = 2
       Randomtype = 4
       RandomPara(1) = PCP%Tau_mean
       RandomPara(2) = PCP%Tau_para
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tau = rnd_vec(1)
    CASE(5)   ! LogNormal
       ! mean and variance of log(x) should be given
       ! Parameters: (ave,sigma) of ln(x)
       n_par = 4
       Randomtype = 5
       RandomPara(1) = PCP%Tau_mean
       RandomPara(2) = PCP%Tau_para
       RandomPara(3) = PCP%Tau_high  ! high end cutoff
       RandomPara(4) = PCP%Tau_para2 ! shift
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tau = rnd_vec(1)
    CASE(6)   ! Beta
       ! Parameters: (ave,a,b) ave not used
       n_par = 3
       Randomtype = 6
       RandomPara(1) = PCP%Tau_mean
       RandomPara(2) = PCP%Tau_para
       RandomPara(3) = PCP%Tau_para2
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tau = rnd_vec(1)
    CASE(7)   ! Triangular
       ! Parameters: (peak,min,max)
       n_par = 3
       Randomtype = 7
       RandomPara(1) = PCP%Tau_mean
       RandomPara(2) = PCP%Tau_low
       RandomPara(3) = PCP%Tau_high
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tau = rnd_vec(1)
    CASE(8)   ! Weibull  (alpha=1: Exponential)
       ! Parameters: (ave,alpha,lambda)
       n_par = 3
       Randomtype = 8
       RandomPara(1) = PCP%Tau_mean
       RandomPara(2) = PCP%Tau_para
       RandomPara(3) = PCP%Tau_para2
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tau = rnd_vec(1)
    CASE(9)   ! Gumbel
       ! Parameters: (ave,alpha)
       n_par = 2
       Randomtype = 9
       RandomPara(1) = PCP%Tau_mean
       RandomPara(2) = PCP%Tau_para
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tau = rnd_vec(1)
    CASE Default
       CALL SHUTDOWN('ERROR: Class_Properties I_TAU_DIST',PROCESS_0_ONLY=.FALSE.) ; RETURN
    END SELECT
    HR%Tau = Max(HR%Tau, 0.01_EB) ! 0.01s minimum motive force parameter

    SELECT CASE(I_DET_DIST)
    CASE(-1)
       CALL SHUTDOWN('ERROR: Class_Properties: -1',PROCESS_0_ONLY=.FALSE.) ; RETURN
    CASE(0)
       HR%Tdet  = Tdet_mean
    CASE(1)   ! Uniform
       ! Parameters: (ave,min,max) ave not used
       n_par = 3
       Randomtype = 1
       RandomPara(1) = 0.5_EB*(Tdet_high+Tdet_low)
       RandomPara(2) = Tdet_low
       RandomPara(3) = Tdet_high
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tdet = rnd_vec(1)
    CASE(2)   ! Truncated Normal
       ! Parameters: (ave,sigma,min,max)
       n_par = 4
       Randomtype = 2
       RandomPara(1) = Tdet_mean
       RandomPara(2) = Tdet_para
       RandomPara(3) = Tdet_low
       RandomPara(4) = Tdet_high
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tdet = rnd_vec(1)
    CASE(3)   ! Gamma
       ! Parameters: (ave,alpha,beta) ave not used
       n_par = 3
       Randomtype = 3
       RandomPara(1) = Tdet_mean
       RandomPara(2) = Tdet_para
       RandomPara(3) = Tdet_para2
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tdet = rnd_vec(1)
    CASE(4)   ! Normal
       ! Parameters: (ave,sigma)
       n_par = 2
       Randomtype = 4
       RandomPara(1) = Tdet_mean
       RandomPara(2) = Tdet_para
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tdet = rnd_vec(1)
    CASE(5)   ! LogNormal
       ! mean and variance of log(x) should be given
       ! Parameters: (ave,sigma) of ln(x)
       n_par = 4
       Randomtype = 5
       RandomPara(1) = Tdet_mean
       RandomPara(2) = Tdet_para
       RandomPara(3) = Tdet_high  ! high end cutoff
       RandomPara(4) = Tdet_para2 ! shift
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tdet = rnd_vec(1)
    CASE(6)   ! Beta
       ! Parameters: (ave,a,b) ave not used
       n_par = 3
       Randomtype = 6
       RandomPara(1) = Tdet_mean
       RandomPara(2) = Tdet_para
       RandomPara(3) = Tdet_para2
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tdet = rnd_vec(1)
    CASE(7)   ! Triangular
       ! Parameters: (peak,min,max)
       n_par = 3
       Randomtype = 7
       RandomPara(1) = Tdet_mean
       RandomPara(2) = Tdet_low
       RandomPara(3) = Tdet_high
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tdet = rnd_vec(1)
    CASE(8)   ! Weibull  (alpha=1: Exponential)
       ! Parameters: (ave,alpha,lambda)
       n_par = 3
       Randomtype = 8
       RandomPara(1) = Tdet_mean
       RandomPara(2) = Tdet_para
       RandomPara(3) = Tdet_para2
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tdet = rnd_vec(1)
    CASE(9)   ! Gumbel
       ! Parameters: (ave,alpha)
       n_par = 2
       Randomtype = 9
       RandomPara(1) = Tdet_mean
       RandomPara(2) = Tdet_para
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tdet = rnd_vec(1)
    CASE Default
       CALL SHUTDOWN('ERROR: Class_Properties I_DET_DIST',PROCESS_0_ONLY=.FALSE.) ; RETURN
    END SELECT

    SELECT CASE(I_PRE_DIST)
    CASE(-1)
       CALL SHUTDOWN('ERROR: Class_Properties: -1',PROCESS_0_ONLY=.FALSE.) ; RETURN
    CASE(0)
       HR%Tpre  = MAX(0._EB,Tpre_mean)
    CASE(1)   ! Uniform
       ! Parameters: (ave,min,max) ave not used
       n_par = 3
       Randomtype = 1
       RandomPara(1) = 0.5_EB*(Tpre_high+Tpre_low)
       RandomPara(2) = Tpre_low
       RandomPara(3) = Tpre_high
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tpre = MAX(0._EB,rnd_vec(1))
    CASE(2)   ! Truncated Normal
       ! Parameters: (ave,sigma,min,max)
       n_par = 4
       Randomtype = 2
       RandomPara(1) = Tpre_mean
       RandomPara(2) = Tpre_para
       RandomPara(3) = Tpre_low
       RandomPara(4) = Tpre_high
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tpre = MAX(0._EB,rnd_vec(1))
    CASE(3)   ! Gamma
       ! Parameters: (ave,alpha,beta) ave not used
       n_par = 3
       Randomtype = 3
       RandomPara(1) = Tpre_mean
       RandomPara(2) = Tpre_para
       RandomPara(3) = Tpre_para2
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tpre = MAX(0._EB,rnd_vec(1))
    CASE(4)   ! Normal
       ! Parameters: (ave,sigma)
       n_par = 2
       Randomtype = 4
       RandomPara(1) = Tpre_mean
       RandomPara(2) = Tpre_para
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tpre = MAX(0._EB,rnd_vec(1))
    CASE(5)   ! LogNormal
       ! mean and variance of log(x) should be given
       ! Parameters: (ave,sigma) of ln(x)
       n_par = 4
       Randomtype = 5
       RandomPara(1) = Tpre_mean
       RandomPara(2) = Tpre_para
       RandomPara(3) = Tpre_high  ! high end cutoff
       RandomPara(4) = Tpre_para2 ! shift
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tpre = MAX(0._EB,rnd_vec(1))
    CASE(6)   ! Beta
       ! Parameters: (ave,a,b) ave not used
       n_par = 3
       Randomtype = 6
       RandomPara(1) = Tpre_mean
       RandomPara(2) = Tpre_para
       RandomPara(3) = Tpre_para2
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tpre = MAX(0._EB,rnd_vec(1))
    CASE(7)   ! Triangular
       ! Parameters: (peak,min,max)
       n_par = 3
       Randomtype = 7
       RandomPara(1) = Tpre_mean
       RandomPara(2) = Tpre_low
       RandomPara(3) = Tpre_high
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tpre = MAX(0._EB,rnd_vec(1))
    CASE(8)   ! Weibull  (alpha=1: Exponential)
       ! Parameters: (ave,alpha,lambda)
       n_par = 3
       Randomtype = 8
       RandomPara(1) = Tpre_mean
       RandomPara(2) = Tpre_para
       RandomPara(3) = Tpre_para2
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tpre = MAX(0._EB,rnd_vec(1))
    CASE(9)   ! Gumbel
       ! Parameters: (ave,alpha)
       n_par = 2
       Randomtype = 9
       RandomPara(1) = Tpre_mean
       RandomPara(2) = Tpre_para
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tpre = MAX(0._EB,rnd_vec(1))
    CASE Default
       CALL SHUTDOWN('ERROR: Class_Properties I_PRE_DIST',PROCESS_0_ONLY=.FALSE.) ; RETURN
    END SELECT
    !
    ! Constants for the 'psychological' potential
    HR%A      = PCP%A
    HR%B      = PCP%B
    HR%Lambda = PCP%Lambda
    HR%C_Young = PCP%C_Young
    HR%Gamma   = PCP%Gamma
    HR%Kappa   = PCP%Kappa
    !
    HR%r_torso    = HR%Radius*(2.0_EB*PCP%r_torso/PCP%D_mean)
    HR%r_shoulder = HR%Radius*(2.0_EB*PCP%r_shoulder/PCP%D_mean)
    HR%d_shoulder = HR%Radius-HR%r_shoulder
    HR%m_iner     = PCP%m_iner*(HR%r_torso**2)**2/(PCP%r_torso**2)**2
    HR%Tau_iner   = PCP%Tau_iner*(HR%Tau/PCP%Tau_mean) ! s
    HR%Torque = 0.0_EB
    HR%angle_old = 0.0_EB  ! rad
    HR%Omega = 0.0_EB  ! rad/s
    !Issue1547: Nervousness parameter initialization, put one minute non-nervous time as default.
    !Issue1547: Some time is needed, because the agents are initially at rest, so the initial
    !Issue1547: acceleration time is not counted. Time units are seconds, so 60 seconds.
    HR%Speed_ave = 1.0_EB
    ! HR%T_Speed_ave = 60.0_EB
    HR%T_Speed_ave = TIME_INIT_NERVOUSNESS
    HR%ID_CB = -1 ! Crowbar project, ids start at zero
    !
  END SUBROUTINE CLASS_PROPERTIES

  SUBROUTINE TPRE_GENERATION(i_pre_dist,Tpre_low,Tpre_high,Tpre_mean,Tpre_para,Tpre_para2,TPRE_OUT)
    IMPLICIT NONE
    !
    ! Passed variables
    INTEGER, INTENT(IN) :: i_pre_dist
    REAL(EB), INTENT(IN) :: Tpre_low,Tpre_high,Tpre_mean,Tpre_para,Tpre_para2
    REAL(EB), INTENT(OUT) :: TPRE_OUT
    !
    ! Local variables
    ! How many rnd numbers per one call to the rnd routine
    INTEGER, PARAMETER  :: n_rnd=1, n_max_par=4
    INTEGER  :: n_par, RandomType
    ! No more than 4 numbers needed to specify the distributions
    REAL(EB) :: RandomPara(n_max_par)
    REAL(EB) :: rnd_vec(n_rnd)
    ! 1: uniform (TESTED: OK)
    ! 2: Truncated normal (TESTED: OK)
    ! 3: gamma  (TESTED: OK)
    ! 4: normal (TESTED: OK)
    ! 5: lognormal (TESTED: OK)
    ! 6: beta (TESTED: OK)
    ! 7: Triangular (TESTED: OK)
    ! 8: Weibull (TESTED: OK) (alpha=1: Exponential)
    ! 9: Gumbel (TESTED: OK)

    SELECT CASE(I_PRE_DIST)
    CASE(-1)
       CALL SHUTDOWN('ERROR: Tpre_Generation: -1',PROCESS_0_ONLY=.FALSE.) ; RETURN
    CASE(0)
       TPRE_OUT = MAX(0._EB,Tpre_mean)
    CASE(1)   ! Uniform
       ! Parameters: (ave,min,max) ave not used
       n_par = 3
       Randomtype = 1
       RandomPara(1) = 0.5_EB*(Tpre_high+Tpre_low)
       RandomPara(2) = Tpre_low
       RandomPara(3) = Tpre_high
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       TPRE_OUT = MAX(0._EB,rnd_vec(1))
    CASE(2)   ! Truncated Normal
       ! Parameters: (ave,sigma,min,max)
       n_par = 4
       Randomtype = 2
       RandomPara(1) = Tpre_mean
       RandomPara(2) = Tpre_para
       RandomPara(3) = Tpre_low
       RandomPara(4) = Tpre_high
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       TPRE_OUT = MAX(0._EB,rnd_vec(1))
    CASE(3)   ! Gamma
       ! Parameters: (ave,alpha,beta) ave not used
       n_par = 3
       Randomtype = 3
       RandomPara(1) = Tpre_mean
       RandomPara(2) = Tpre_para
       RandomPara(3) = Tpre_para2
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       TPRE_OUT = MAX(0._EB,rnd_vec(1))
    CASE(4)   ! Normal
       ! Parameters: (ave,sigma)
       n_par = 2
       Randomtype = 4
       RandomPara(1) = Tpre_mean
       RandomPara(2) = Tpre_para
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       TPRE_OUT = MAX(0._EB,rnd_vec(1))
    CASE(5)   ! LogNormal
       ! mean and variance of log(x) should be given
       ! Parameters: (ave,sigma) of ln(x)
       n_par = 4
       Randomtype = 5
       RandomPara(1) = Tpre_mean
       RandomPara(2) = Tpre_para
       RandomPara(3) = Tpre_high  ! high end cutoff
       RandomPara(4) = Tpre_para2 ! shift
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       TPRE_OUT = MAX(0._EB,rnd_vec(1))
    CASE(6)   ! Beta
       ! Parameters: (ave,a,b) ave not used
       n_par = 3
       Randomtype = 6
       RandomPara(1) = Tpre_mean
       RandomPara(2) = Tpre_para
       RandomPara(3) = Tpre_para2
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       TPRE_OUT = MAX(0._EB,rnd_vec(1))
    CASE(7)   ! Triangular
       ! Parameters: (peak,min,max)
       n_par = 3
       Randomtype = 7
       RandomPara(1) = Tpre_mean
       RandomPara(2) = Tpre_low
       RandomPara(3) = Tpre_high
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       TPRE_OUT = MAX(0._EB,rnd_vec(1))
    CASE(8)   ! Weibull  (alpha=1: Exponential)
       ! Parameters: (ave,alpha,lambda)
       n_par = 3
       Randomtype = 8
       RandomPara(1) = Tpre_mean
       RandomPara(2) = Tpre_para
       RandomPara(3) = Tpre_para2
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       TPRE_OUT = MAX(0._EB,rnd_vec(1))
    CASE(9)   ! Gumbel
       ! Parameters: (ave,alpha)
       n_par = 2
       Randomtype = 9
       RandomPara(1) = Tpre_mean
       RandomPara(2) = Tpre_para
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       TPRE_OUT = MAX(0._EB,rnd_vec(1))
    CASE Default
       CALL SHUTDOWN('ERROR: Tpre_Generation I_PRE_DIST',PROCESS_0_ONLY=.FALSE.) ; RETURN
    END SELECT
  END SUBROUTINE TPRE_GENERATION
!
  SUBROUTINE RE_ALLOCATE_HUMANS(CODE,NM)
    IMPLICIT NONE
    !
    ! Passed variables
    INTEGER, INTENT(IN) :: CODE,NM
    !
    ! Local variables
    INTEGER IZERO
    TYPE (HUMAN_TYPE), ALLOCATABLE, DIMENSION(:) :: DUMMY
    TYPE (MESH_TYPE), POINTER :: M =>NULL()
    !
    IF (.NOT.ANY(EVACUATION_ONLY)) RETURN
    IF (.NOT.EVACUATION_ONLY(NM)) RETURN

    SELECT CASE(CODE)
       !
    CASE(1)
       !
       M=>MESHES(NM)
       ALLOCATE(DUMMY(1:M%N_HUMANS_DIM),STAT=IZERO)
       DUMMY = M%HUMAN
       !
       DEALLOCATE(M%HUMAN)
       ALLOCATE(M%HUMAN(M%N_HUMANS_DIM+1000),STAT=IZERO)
       M%HUMAN(1:M%N_HUMANS_DIM) = DUMMY(1:M%N_HUMANS_DIM)
       M%N_HUMANS_DIM = M%N_HUMANS_DIM+1000
       !
    END SELECT
    !
    DEALLOCATE(DUMMY)
    !
  END SUBROUTINE RE_ALLOCATE_HUMANS
!
  SUBROUTINE DUMP_EVAC(T,NM)
    IMPLICIT NONE
    !
    ! Passed variables
    INTEGER, INTENT(IN) :: NM
    REAL(EB), INTENT(IN) :: T
    !
    ! Local variables
    INTEGER :: NPP,NPLIM,i,izero,nn,n
    LOGICAL :: CROWBAR_DUMP
    REAL(EB) :: TNOW, EVEL, angle_hr
    REAL(FB) :: PART_MIN, PART_MAX
    REAL(FB), ALLOCATABLE, DIMENSION(:) :: XP,YP,ZP
    REAL(FB), ALLOCATABLE, DIMENSION(:,:) :: QP, AP
    INTEGER, ALLOCATABLE, DIMENSION(:) :: TA
    TYPE (HUMAN_TYPE), POINTER :: HR =>NULL()
    INTEGER, PARAMETER :: PART_BOUNDFILE_VERSION=1
    !
    IF (.NOT.ANY(EVACUATION_ONLY)) RETURN
    IF (.NOT.(EVACUATION_ONLY(NM) .AND. EMESH_INDEX(NM)>0)) RETURN
    TNOW=CURRENT_TIME()
    !
    CALL POINT_TO_MESH(NM)

    ! Write the current time to the prt5 file, then start looping through the particle classes
    WRITE(LU_PART(NM)) REAL(T,FB)

    WRITE(LU_PART(NM+NMESHES),'(ES13.6,1X,I4,1X,I4)')T, N_EVAC, PART_BOUNDFILE_VERSION

    HUMAN_CLASS_LOOP: DO N = 1, N_EVAC
       ! Count the number of humans to dump out

       CROWBAR_DUMP = .FALSE.
       IF (N==N_EVAC .AND. TRIM(EVAC_CLASS_NAME(N_EVAC))=='Crowbar') CROWBAR_DUMP = .TRUE.

       NPLIM = 0
       HUMLOOP: DO I = 1, N_HUMANS
          HR => HUMAN(I)
          ! Crowbar: Dump the read in positions
          IF (CROWBAR_DUMP .AND. N<N_EVAC) CYCLE HUMLOOP  ! The last avatar type is Crowbar
          IF (CROWBAR_DUMP .AND. .NOT.HR%CROWBAR_READ_IN) CYCLE HUMLOOP ! Just the read in positions are dumped
          IF (HR%COLOR_INDEX < 1) HR%COLOR_INDEX = 1
          IF (HR%COLOR_INDEX > EVAC_AVATAR_NCOLOR) HR%COLOR_INDEX = EVAC_AVATAR_NCOLOR
          IF (HR%SHOW) THEN
             IF (CROWBAR_DUMP) THEN
                IF (HR%CROWBAR_READ_IN) NPLIM = NPLIM + 1
             ELSE
                IF (HR%IEL > 0) THEN
                   IF (EVAC_EVACS(HR%IEL)%AVATAR_TYPE_INDEX == N) NPLIM = NPLIM + 1
                ELSE
                   IF (EVAC_ENTRYS(ABS(HR%IEL))%AVATAR_TYPE_INDEX == N) NPLIM = NPLIM + 1
                END IF
             END IF
          END IF
       END DO HUMLOOP
       NPLIM = MIN(NPPS,NPLIM)

       ALLOCATE(TA(NPLIM),STAT=IZERO)
       CALL ChkMemErr('DUMP_EVAC','TA',IZERO)
       ALLOCATE(XP(NPLIM),STAT=IZERO)
       CALL ChkMemErr('DUMP_EVAC','XP',IZERO)
       ALLOCATE(YP(NPLIM),STAT=IZERO)
       CALL ChkMemErr('DUMP_EVAC','YP',IZERO)
       ALLOCATE(ZP(NPLIM),STAT=IZERO)
       CALL ChkMemErr('DUMP_EVAC','ZP',IZERO)
       ! body angle, semi major axis, semi minor axis
       ALLOCATE(AP(NPLIM,4),STAT=IZERO)
       CALL ChkMemErr('DUMP_EVAC','AP',IZERO)
       IF (EVAC_N_QUANTITIES > 0) THEN
          ALLOCATE(QP(NPLIM,EVAC_N_QUANTITIES),STAT=IZERO)
          CALL ChkMemErr('DUMP_EVAC','QP',IZERO)
       END IF

       ! Load human coordinates into single precision array

       NPP = 0
       PLOOP: DO I = 1, N_HUMANS
          HR => HUMAN(I)
          IF (.NOT. HR%SHOW ) CYCLE PLOOP
          IF (CROWBAR_DUMP .AND. N<N_EVAC) CYCLE PLOOP  ! The last avatar type is Crowbar
          IF (CROWBAR_DUMP .AND. .NOT.HR%CROWBAR_READ_IN) CYCLE PLOOP ! Just the read in positions are dumped
          IF (.NOT.CROWBAR_DUMP) THEN
             IF (HR%IEL > 0) THEN
                IF (.NOT.(EVAC_EVACS(HR%IEL)%AVATAR_TYPE_INDEX == N)) CYCLE PLOOP
             ELSE
                IF (.NOT.(EVAC_ENTRYS(ABS(HR%IEL))%AVATAR_TYPE_INDEX == N)) CYCLE PLOOP
             END IF
          END IF

          NPP = NPP + 1
          TA(NPP) = HR%ILABEL
          XP(NPP) = REAL(HR%X,FB)
          YP(NPP) = REAL(HR%Y,FB)

          ZP(NPP) = MIN( MAX(REAL(HR%Z,FB), EVAC_Z_MIN), EVAC_Z_MAX)

          AP(NPP,1) = 180.0_FB*REAL(HR%Angle/Pi,FB)
          AP(NPP,2) =   2.0_FB*REAL(HR%Radius,FB) ! diameter
          AP(NPP,3) =   2.0_FB*REAL(HR%r_torso,FB) ! diameter
          ! Height of a human scaled by radius, default male 1.80 m
          AP(NPP,4) =  1.80_FB*REAL(HR%Radius/0.27_EB,FB)
          
          IF (CROWBAR_DUMP) THEN
             EVEL = SQRT(HR%U_CB**2 + HR%V_CB**2)
             IF (EVEL >= TWO_EPSILON_EB) THEN
                IF (HR%V_CB >= 0.0_EB) THEN
                   angle_hr = ACOS(HR%U_CB/EVEL)
                ELSE
                   angle_hr = 2.0_EB*Pi - ACOS(HR%U_CB/EVEL)
                END IF
                IF (angle_hr == 2.0_EB*Pi) angle_hr = 0.0_EB  ! agent HR angle is [0,2Pi)
             ELSE
                angle_hr = HR%ANGLE
             END IF
             XP(NPP) = REAL(HR%X_CB,FB)
             YP(NPP) = REAL(HR%Y_CB,FB)
             ! Shift the crowbar avatars a little bit upwards
             ZP(NPP) = MIN(MAX(REAL(HR%Z+0.1_EB,FB), EVAC_Z_MIN), EVAC_Z_MAX)
             AP(NPP,1) = REAL(180.0_EB*angle_hr/Pi,FB) ! read in velocity direction
             AP(NPP,2) =   0.1_FB ! diameter
             AP(NPP,3) =   0.05_FB ! torso diameter
             AP(NPP,4) =  1.00_FB ! height
          END IF

          DO NN=1,EVAC_N_QUANTITIES
             SELECT CASE(EVAC_QUANTITIES_INDEX(NN))
             CASE(240)  ! MOTIVE_ACCELERATION, Unimpeded walking Speed / tau
                EVEL = SQRT(HR%UBAR**2 + HR%VBAR**2)
                IF (EVEL >= TWO_EPSILON_EB) THEN
                   QP(NPP,NN) = REAL(HR%Speed/HR%Tau,FB)
                ELSE
                   QP(NPP,NN) = REAL(EVEL,FB)
                END IF
             CASE(241)  ! FED_DOSE, Fractional Effective Dose
                QP(NPP,NN) = REAL(HR%IntDose,FB)
                IF (CROWBAR_DUMP .AND. HR%CROWBAR_READ_IN) QP(NPP,NN) = REAL(HR%ID_CB,FB)
             CASE(242)  ! SPEED, Human speed
                QP(NPP,NN) = REAL(SQRT(HR%U**2 + HR%V**2),FB)
             CASE(243)  ! ANGULAR_SPEED, Human Angular Velocity
                QP(NPP,NN) = REAL(HR%Omega,FB)
             CASE(244)  ! ACCELERATION, Human acc.
                QP(NPP,NN) = REAL(SQRT( ((HR%UBAR-HR%U)/HR%Tau + (HR%F_X/HR%Mass))**2 + &
                     ((HR%VBAR-HR%V)/HR%Tau + (HR%F_Y/HR%Mass))**2 ),FB)
             CASE(245)  ! CONTACT_LINEFORCE, Human Pressure: contact forces
                QP(NPP,NN) = REAL(HR%SumForces ,FB)
             CASE(246)  ! TOTAL_LINEFORCE, Human Pressure2: contact + social
                QP(NPP,NN) = REAL(HR%SumForces2 ,FB)
             CASE(247)  ! COLOR, Human color index
                QP(NPP,NN) = REAL(HR%COLOR_INDEX - 1,FB)
             CASE(248)  ! MOTIVE_ANGLE,
                EVEL = SQRT(HR%UBAR**2 + HR%VBAR**2)
                IF (EVEL >= TWO_EPSILON_EB) THEN
                   IF (HR%VBAR >= 0.0_EB) THEN
                      angle_hr = ACOS(HR%UBAR/EVEL)
                   ELSE
                      angle_hr = 2.0_EB*Pi - ACOS(HR%UBAR/EVEL)
                   END IF
                   IF (angle_hr == 2.0_EB*Pi) angle_hr = 0.0_EB  ! agent HR angle is [0,2Pi)
                ELSE
                   angle_hr = 0.0_EB
                END IF
                QP(NPP,NN) = REAL(angle_hr*180.0_EB/Pi,FB)
             CASE(249)  ! DENSITY, agent density, 1/m2
                QP(NPP,NN) = REAL(HR%Density ,FB)
             CASE(250)  !Issue1547: Nervousness parameter
                !Issue1547: Speed_ave is average speed versus v0_i (v_ave/v0_i)
                !Issue1547: Nervousness = 1 - (v_ave/v0_i)
                QP(NPP,NN) = REAL(MAX(1.0_EB - HR%Speed_ave,0.0_EB), FB)
             END SELECT
          END DO

          IF (NPP>=NPPS) EXIT PLOOP
       END DO PLOOP
       !
       ! Dump human data into the .prt5 file
       !
       WRITE(LU_PART(NM)) NPLIM
       WRITE(LU_PART(NM)) (XP(I),I=1,NPLIM),(YP(I),I=1,NPLIM),(ZP(I),I=1,NPLIM), &
            (AP(I,1),I=1,NPLIM),(AP(I,2),I=1,NPLIM),(AP(I,3),I=1,NPLIM),(AP(I,4),I=1,NPLIM)
       WRITE(LU_PART(NM)) (TA(I),I=1,NPLIM)
       IF (EVAC_N_QUANTITIES > 0) THEN
          WRITE(LU_PART(NM)) ((QP(I,NN),I=1,NPLIM),NN=1,EVAC_N_QUANTITIES)
       END IF
       !
       WRITE(LU_PART(NM+NMESHES),'(I4,1X,I7)') EVAC_N_QUANTITIES, NPLIM
       IF (EVAC_N_QUANTITIES > 0) THEN
          DO NN = 1, EVAC_N_QUANTITIES
             IF (NPLIM > 0) THEN
                PART_MAX = QP(1,NN)
                PART_MIN = PART_MAX
                DO I = 2, NPLIM
                   PART_MIN = MIN(QP(I,NN),PART_MIN)
                   PART_MAX = MAX(QP(I,NN),PART_MAX)
                END DO
             ELSE
                PART_MIN = 1.0_FB
                PART_MAX = 0.0_FB
             ENDIF
             WRITE(LU_PART(NM+NMESHES),'(5X,ES13.6,1X,ES13.6)')PART_MIN, PART_MAX
          ENDDO
          DEALLOCATE(QP)
       END IF
       DEALLOCATE(AP)
       DEALLOCATE(ZP)
       DEALLOCATE(YP)
       DEALLOCATE(XP)
       DEALLOCATE(TA)

    END DO HUMAN_CLASS_LOOP

    !
    T_USED(12) = T_USED(12) + CURRENT_TIME() - TNOW
  END SUBROUTINE DUMP_EVAC
!
  FUNCTION GaussRand( gmean, gtheta, gcutmult )
    IMPLICIT NONE
    !
    ! Random numbers from the Gaussian distribution
    !
    REAL(EB) GaussRand
    !  generates a random number (x) with
    !  P(x) = exp[- (x-gmean)^2 / (2*gtheta)], if x is in
    !         [gmean - gcutmult*sqrt(gtheta), gmean + gcutmult*sqrt(gtheta)]
    !       = 0                              , if not
    !
    ! Passed variables
    REAL(EB) gmean, gtheta, gcutmult

    ! Local variables
    REAL :: RN_REAL
    REAL(EB) v1,v2,rsq,fac, RN

    IF ( (GaussFlag == 1)  .AND. (ABS(GaussSet2-gmean) <= gcutmult*SQRT(gtheta)) ) THEN
       GaussFlag = 0
       GaussRand = GaussSet2
    ELSE
       GaussFlag = 0
       DO
          CALL RANDOM_NUMBER(RN_REAL)
          RN = REAL(RN_REAL,EB)
          v1 = 1.0_EB - 2.0_EB*rn
          CALL RANDOM_NUMBER(RN_REAL)
          RN = REAL(RN_REAL,EB)
          v2 = 1.0_EB - 2.0_EB*rn
          rsq = v1*v1 + v2*v2
          DO WHILE (rsq >= 1.0 )
             CALL RANDOM_NUMBER(RN_REAL)
             RN = REAL(RN_REAL,EB)
             v1 = 1.0_EB - 2.0_EB*rn
             CALL RANDOM_NUMBER(RN_REAL)
             RN = REAL(RN_REAL,EB)
             v2 = 1.0_EB - 2.0_EB*rn
             rsq=v1*v1+v2*v2
          END DO
          fac = SQRT(-2.0_EB*gtheta*LOG(rsq)/rsq)
          GaussSet1 = v1*fac + gmean
          GaussSet2 = v2*fac + gmean

          IF ( (ABS(GaussSet1-gmean) <= gcutmult*SQRT(gtheta)) .OR. (ABS(GaussSet2-gmean) <= gcutmult*SQRT(gtheta)) ) THEN
             EXIT
          END IF

       END DO

       IF (ABS(GaussSet1-gmean) <= gcutmult*SQRT(gtheta)) THEN
          GaussFlag = 1
          GaussRand = GaussSet1
       ELSE
          GaussFlag = 0
          GaussRand = GaussSet2
       END IF
    END IF

  END FUNCTION GaussRand
!
  FUNCTION GaussTrun( gmean, gsigma, glow, ghigh )
    IMPLICIT NONE
    !
    ! Random numbers from the Gaussian distribution
    !
    REAL(EB) GaussTrun
    !  generates a random number (x) with
    !  P(x) = exp[- (x-gmean)^2 / (2*gsigma^2)], if x is in
    !         [gmean - gcutmult*gsigma, gmean + gcutmult*gsigma]
    !       = 0                              , if not
    !
    ! Passed variables
    REAL(EB) gmean, gsigma, glow, ghigh

    ! Local variables
    REAL :: RN_REAL
    REAL(EB) v1, v2, rsq, fac, RN

    IF ( (GTrunFlag == 1)  .AND. (GTrunSet2 >= glow) .AND. (GTrunSet2 <= ghigh) ) THEN
       GTrunFlag = 0
       GaussTrun = GTrunSet2
    ELSE
       GTrunFlag = 0
       DO
          CALL RANDOM_NUMBER(RN_REAL)
          RN = REAL(RN_REAL,EB)
          v1 = 1.0_EB - 2.0_EB*rn
          CALL RANDOM_NUMBER(RN_REAL)
          RN = REAL(RN_REAL,EB)
          v2 = 1.0_EB - 2.0_EB*rn
          rsq = v1*v1 + v2*v2
          DO WHILE (rsq >= 1.0 )
             CALL RANDOM_NUMBER(RN_REAL)
             RN = REAL(RN_REAL,EB)
             v1 = 1.0_EB - 2.0_EB*rn
             CALL RANDOM_NUMBER(RN_REAL)
             RN = REAL(RN_REAL,EB)
             v2 = 1.0_EB - 2.0_EB*rn
             rsq = v1*v1 + v2*v2
          END DO
          fac = SQRT(-2.0_EB*(gsigma**2)*LOG(rsq)/rsq)
          GTrunSet1 = v1*fac + gmean
          GTrunSet2 = v2*fac + gmean

          IF ( ((GTrunSet1 >= glow) .AND. (GTrunSet1 <= ghigh)) .OR. ((GTrunSet2 >= glow) .AND. (GTrunSet2 <= ghigh)) ) THEN
             EXIT
          END IF

       END DO

       IF ( (GTrunSet1 >= glow) .AND. (GTrunSet1 <= ghigh) ) THEN
          GTrunFlag = 1
          GaussTrun = GTrunSet1
       ELSE
          GTrunFlag = 0
          GaussTrun = GTrunSet2
       END IF
    END IF

  END FUNCTION GaussTrun
!
  SUBROUTINE DUMP_EVAC_CSV(Tin)
    IMPLICIT NONE
    !
    ! Dump agent data to CHID_evac.csv
    ! This subroutine is called from the main program.
    !
    ! Passed variables
    REAL(EB), INTENT(IN) :: Tin
    !
    ! Local variables
    CHARACTER(80) tcform
    INTEGER n_cols, n_tot_humans, i, ii, izero, ii_ntargets, ii_density
    REAL(FB), ALLOCATABLE, DIMENSION(:) :: ITEMP
    !
    IF (.NOT.ANY(EVACUATION_ONLY)) RETURN
    !
    ALLOCATE(ITEMP(MAX(1,N_EXITS+N_DOORS)), STAT = IZERO)
    CALL ChkMemErr('DUMP_EVAC_CSV','ITEMP', IZERO)
    !
    ! Output first the floors then the corridors
    n_cols = n_egrids + n_corrs + N_EXITS + N_DOORS + 1 + N_EXITS - n_co_exits + N_DOORS
    n_tot_humans = 0
    DO i = 1, n_egrids
       n_tot_humans = n_tot_humans + MESHES(EVAC_Node_List(i)%IMESH)%N_HUMANS
    END DO
    DO i = 1, n_corrs
       n_tot_humans = n_tot_humans + EVAC_CORRS(i)%n_inside
    END DO
    !
    ii = 0
    DO i = 1, N_EXITS
       IF (.NOT.EVAC_EXITS(i)%COUNT_ONLY) THEN
          ii = ii + 1
          ITEMP(ii) = REAL(EVAC_EXITS(i)%NTARGET(50),FB)
       END IF
    END DO
    DO i = 1, N_DOORS
       ii = ii + 1
       ITEMP(ii) = REAL(EVAC_DOORS(i)%NTARGET(50),FB)
    END DO
    ii_ntargets = ii
    DO i = 1, N_EXITS
       IF (EVAC_EXITS(i)%COUNT_ONLY .AND. EVAC_EXITS(i)%COUNT_DENSITY) THEN
          ii = ii + 1
          ITEMP(ii) = REAL(EVAC_EXITS(i)%ICOUNT/MAX(0.0001_EB,EVAC_EXITS(i)%Width),FB)
       END IF
    END DO
    ii_density = ii
    !
    IF (n_dead >= 0) THEN
       ! Write the 'fed' columns
       IF (ii_density > ii_ntargets) THEN
          WRITE(tcform,'(a,i4.4,a,a,i4.4,a,a)') "(ES13.5E3,", n_cols, "(',',i8)", "," , &
               ii_density-ii_ntargets, "(',',ES13.5E3)", ",',',i8,',',ES13.5E3,',',ES13.5E3)"
          WRITE (LU_EVACCSV,fmt=tcform) Tin, n_tot_humans, &
               (MESHES(EVAC_Node_List(i)%IMESH)%N_HUMANS, i=1,n_egrids), &
               (EVAC_CORRS(i)%n_inside, i = 1,n_corrs), &
               (EVAC_EXITS(i)%ICOUNT, i = 1,N_EXITS), &
               (EVAC_DOORS(i)%ICOUNT, i = 1,N_DOORS), &
               (NINT(ITEMP(i)), i = 1,N_EXITS-n_co_exits+N_DOORS), &
               (ITEMP(i), i = ii_ntargets+1,ii_density), &
               n_dead, fed_max, fed_max_alive
       ELSE
          WRITE(tcform,'(a,i4.4,a,a)') "(ES13.5E3,",n_cols+1, &
               "(',',i8)", ",',',ES13.5E3,',',ES13.5E3)"
          WRITE (LU_EVACCSV,fmt=tcform) Tin, n_tot_humans, &
               (MESHES(EVAC_Node_List(i)%IMESH)%N_HUMANS, i=1,n_egrids), &
               (EVAC_CORRS(i)%n_inside, i = 1,n_corrs), &
               (EVAC_EXITS(i)%ICOUNT, i = 1,N_EXITS), &
               (EVAC_DOORS(i)%ICOUNT, i = 1,N_DOORS), &
               (NINT(ITEMP(i)), i = 1,N_EXITS-n_co_exits+N_DOORS), &
               n_dead, fed_max, fed_max_alive
       END IF
    ELSE
       ! Do not write the 'fed' columns
       IF (ii_density > ii_ntargets) THEN
          IF (C_HAWK >= 0.0_EB) THEN
             WRITE(tcform,'(a,i4.4,a,a,i4.4,a,i4.4,a)') "(ES13.5E3,", n_cols, "(',',i8)", "," , &
                  ii_density-ii_ntargets, "(',',ES13.5E3),", n_egrids*4, "(',',i8))"
             WRITE (LU_EVACCSV,fmt=tcform) Tin, n_tot_humans, &
                  (MESHES(EVAC_Node_List(i)%IMESH)%N_HUMANS, i=1,n_egrids), &
                  (EVAC_CORRS(i)%n_inside, i = 1,n_corrs), &
                  (EVAC_EXITS(i)%ICOUNT, i = 1,N_EXITS), &
                  (EVAC_DOORS(i)%ICOUNT, i = 1,N_DOORS), &
                  (NINT(ITEMP(i)), i = 1,N_EXITS-n_co_exits+N_DOORS), &
                  (ITEMP(i), i = ii_ntargets+1,ii_density), &
                  ((N_HawkDoveCount(ii,i),ii=1,4),i=1,n_egrids)
          ELSE
             WRITE(tcform,'(a,i4.4,a,a,i4.4,a)') "(ES13.5E3,", n_cols, "(',',i8)", "," , &
                  ii_density-ii_ntargets, "(',',ES13.5E3))"
             WRITE (LU_EVACCSV,fmt=tcform) Tin, n_tot_humans, &
                  (MESHES(EVAC_Node_List(i)%IMESH)%N_HUMANS, i=1,n_egrids), &
                  (EVAC_CORRS(i)%n_inside, i = 1,n_corrs), &
                  (EVAC_EXITS(i)%ICOUNT, i = 1,N_EXITS), &
                  (EVAC_DOORS(i)%ICOUNT, i = 1,N_DOORS), &
                  (NINT(ITEMP(i)), i = 1,N_EXITS-n_co_exits+N_DOORS), &
                  (ITEMP(i), i = ii_ntargets+1,ii_density)
          END IF
       ELSE
          IF (C_HAWK >= 0.0_EB) THEN
             WRITE(tcform,'(a,i4.4,a,i4.4,a)') "(ES13.5E3,",n_cols, "(',',i8),", n_egrids*4,"(',',i8))"
             WRITE (LU_EVACCSV,fmt=tcform) Tin, n_tot_humans, &
                  (MESHES(EVAC_Node_List(i)%IMESH)%N_HUMANS, i=1,n_egrids), &
                  (EVAC_CORRS(i)%n_inside, i = 1,n_corrs), &
                  (EVAC_EXITS(i)%ICOUNT, i = 1,N_EXITS), &
                  (EVAC_DOORS(i)%ICOUNT, i = 1,N_DOORS), &
                  (NINT(ITEMP(i)), i = 1,N_EXITS-n_co_exits+N_DOORS), &
                  ((N_HawkDoveCount(ii,i),ii=1,4), i=1,n_egrids)
          ELSE
             WRITE(tcform,'(a,i4.4,a)') "(ES13.5E3,",n_cols, "(',',i8),i8)"
             WRITE (LU_EVACCSV,fmt=tcform) Tin, n_tot_humans, &
                  (MESHES(EVAC_Node_List(i)%IMESH)%N_HUMANS, i=1,n_egrids), &
                  (EVAC_CORRS(i)%n_inside, i = 1,n_corrs), &
                  (EVAC_EXITS(i)%ICOUNT, i = 1,N_EXITS), &
                  (EVAC_DOORS(i)%ICOUNT, i = 1,N_DOORS), &
                  (NINT(ITEMP(i)), i = 1,N_EXITS-n_co_exits+N_DOORS)
          END IF
       END IF
    END IF
    DEALLOCATE(ITEMP)
    !
  END SUBROUTINE DUMP_EVAC_CSV

  LOGICAL FUNCTION See_each_other(nm, r1_x, r1_y, r2_x, r2_y)
    ! This function returns true, if the two points have a line-of-sight.
    ! This function does not use smoke information, i.e., it just sees if
    ! there are obstacles between the two points.
    ! Inputs:  nm: mesh index, r1 an r2 should belong to the same mesh
    !          (r1_x,r1_y): co-ordinates of the first agent
    !          (r2_x,r2_y): co-ordinates of the second agent
    ! NOTE: This works for thick and thin OBSTs.
    !
    ! Passed variables
    INTEGER, INTENT(IN) :: nm
    REAL(EB), INTENT(IN) :: r1_x, r1_y, r2_x, r2_y
    !
    ! Local variables
    INTEGER :: i, j, isx, isy, i_r1, i_r2, j_r1, j_r2
    INTEGER :: i_old, j_old, ic, ic2, iw, iw1, iw2
    REAL(EB) :: x, y
    TYPE (MESH_TYPE), POINTER :: M =>NULL()

    M => MESHES(NM)
    See_each_other = .TRUE.  ! Default

    isx = INT(SIGN(1.0_EB,r2_x - r1_x)) ! loop increment +1 or -1
    isy = INT(SIGN(1.0_EB,r2_y - r1_y)) ! loop increment +1 or -1
    i_r1 = FLOOR(M%CELLSI(FLOOR((r1_x-M%XS)*M%RDXINT)) + 1.0_EB) ! II start
    i_r2 = FLOOR(M%CELLSI(FLOOR((r2_x-M%XS)*M%RDXINT)) + 1.0_EB) ! II end
    j_r1 = FLOOR(M%CELLSJ(FLOOR((r1_y-M%YS)*M%RDYINT)) + 1.0_EB) ! JJ start
    j_r2 = FLOOR(M%CELLSJ(FLOOR((r2_y-M%YS)*M%RDYINT)) + 1.0_EB) ! JJ end
    i_r1 = MAX(1,MIN(i_r1,M%IBAR)) ! To be sure that indices are always ok.
    i_r2 = MAX(1,MIN(i_r2,M%IBAR))
    j_r1 = MAX(1,MIN(j_r1,M%JBAR))
    j_r2 = MAX(1,MIN(j_r2,M%JBAR))

    ! Same cell: sees always each other
    IF (ABS(i_r2-i_r1)+ABS(j_r2-j_r1) < 1) RETURN

    ! Choose the main direction:
    IF (ABS(r2_x-r1_x)/M%DX(1) < ABS(r2_y-r1_y)/M%DY(1)) THEN
       ! Now y is the main direction
       i_old = i_r1 ; j_old = j_r1
       y = 0.0_EB
       MainLoopY: DO j = j_r1+isy, j_r2, isy
          y = y + isy*M%DY(j)
          x = MAX(M%XS,MIN(M%XF,r1_x + y*(r2_x - r1_x)/(r2_y - r1_y)))
          i = FLOOR(M%CELLSI(FLOOR((x-M%XS)*M%RDXINT)) + 1.0_EB)
          i = isx*MIN(isx*i_r2,isx*i) ! i in interval i_r1...i_r2
          ic  = M%CELL_INDEX(i_old,j_old,1)
          ic2 = M%CELL_INDEX(i    ,j_old,1) ! side cell
          iw  = M%WALL_INDEX(ic, isy*2) ! main direction
          iw1 = M%WALL_INDEX(ic ,isx*1) ! sideways
          iw2 = M%WALL_INDEX(ic2,isy*2) ! side + main direction
          IF (IW >0) THEN
             IF (M%OBSTRUCTION(M%WALL(IW )%OBST_INDEX)%HIDDEN .AND. M%WALL(IW )%OBST_INDEX>0) iw  = 0
          END IF
          IF (IW1>0) THEN
             IF (M%OBSTRUCTION(M%WALL(IW1)%OBST_INDEX)%HIDDEN .AND. M%WALL(IW1)%OBST_INDEX>0) iw1 = 0
          END IF
          IF (IW2>0) THEN
             IF (M%OBSTRUCTION(M%WALL(IW2)%OBST_INDEX)%HIDDEN .AND. M%WALL(IW2)%OBST_INDEX>0) iw2 = 0
          END IF
          ! iw is zero, if there is no solid boundary
          ! from (i,j)==>(i,jnew):    iw and iw2 are zero, iw1 does not matter
          !                           ic=ic2 ==> iw=iw2
          ! from (i,j)==>(inew,jnew): iw1 and iw2 are zero, iw does not matter
          IF ((i==i_old .AND. iw/=0) .OR. (i/=i_old .AND. (iw1/=0 .OR. iw2/=0))) THEN
             See_each_other = .FALSE.
             EXIT MainLoopY
          END IF
          i_old = i ; j_old = j
       END DO MainLoopY
    ELSE
       ! Now x is the main direction
       i_old = i_r1 ; j_old = j_r1
       x = 0.0_EB
       MainLoopX: DO i = i_r1+isx, i_r2, isx
          x = x + isx*M%DX(i)
          y = MAX(M%YS,MIN(M%YF,r1_y + x*(r2_y - r1_y)/(r2_x - r1_x)))
          j = FLOOR(M%CELLSJ(FLOOR((y-M%YS)*M%RDYINT)) + 1.0_EB)
          j = isy*MIN(isy*j_r2,isy*j) ! j in interval j_r1...j_r2
          ic  = M%CELL_INDEX(i_old,j_old,1)
          ic2 = M%CELL_INDEX(i_old,j    ,1) ! side cell
          iw  = M%WALL_INDEX(ic, isx*1) ! main direction
          iw1 = M%WALL_INDEX(ic ,isy*2) ! sideways
          iw2 = M%WALL_INDEX(ic2,isx*1) ! side + main direction
          IF (IW >0) THEN
             IF (M%OBSTRUCTION(M%WALL(IW )%OBST_INDEX)%HIDDEN .AND. M%WALL(IW )%OBST_INDEX>0) iw  = 0
          END IF
          IF (IW1>0) THEN
             IF (M%OBSTRUCTION(M%WALL(IW1)%OBST_INDEX)%HIDDEN .AND. M%WALL(IW1)%OBST_INDEX>0) iw1 = 0
          END IF
          IF (IW2>0) THEN
             IF (M%OBSTRUCTION(M%WALL(IW2)%OBST_INDEX)%HIDDEN .AND. M%WALL(IW2)%OBST_INDEX>0) iw2 = 0
          END IF
          ! iw is zero, if there is no solid boundary
          ! from (i,j)==>(inew,j):    iw and iw2 are zero, iw1 does not matter
          !                           ic=ic2 ==> iw=iw2
          ! from (i,j)==>(inew,jnew): iw1 and iw2 are zero, iw does not matter
          IF ((j==j_old .AND. iw/=0) .OR. (j/=j_old .AND. (iw1/=0 .OR. iw2/=0))) THEN
             See_each_other = .FALSE.
             EXIT MainLoopX
          END IF
          i_old = i ; j_old = j
       END DO MainLoopX
    END IF
  END FUNCTION See_each_other
  !
  LOGICAL FUNCTION See_door(nm, r1_x, r1_y, r2_x, r2_y, ave_K, max_fed)
    ! This function returns true, if the two points have a line-of-sight.
    ! This function does use smoke information, i.e., it sees if
    ! there are obstacles and/or too much smoke between the two points.
    ! Inputs:  nm: mesh index, r1 an r2 should belong to the same mesh
    !          (r1_x,r1_y): co-ordinates of the agent
    !          (r2_x,r2_y): co-ordinates of the door
    ! Outputs: ave_K: average extinction coefficient of the path
    !          max_fed: maximum level of FED at the path
    !
    ! Passed variables
    INTEGER, INTENT(IN) :: nm
    REAL(EB), INTENT(IN) :: r1_x, r1_y, r2_x, r2_y
    REAL(EB), INTENT(OUT) :: ave_K, max_fed
    !
    ! Local variables
    INTEGER :: i, j, isx, isy, i_r1, i_r2, j_r1, j_r2
    INTEGER :: i_old, j_old, ic, ic2, iw, iw1, iw2, NM_SEE
    REAL(EB) :: x, y
    TYPE (MESH_TYPE), POINTER :: M =>NULL(), MFF =>NULL()

    M   => MESHES(NM)
    NM_SEE = NMESHES - N_EGRIDS + EMESH_INDEX(NM)
    MFF => MESHES(NM_SEE)

    See_door = .TRUE.  ! Default
    ave_K = 0.0_EB
    max_fed = 0.0_EB

    isx = INT(SIGN(1.0_EB,r2_x - r1_x)) ! loop increment +1 or -1
    isy = INT(SIGN(1.0_EB,r2_y - r1_y)) ! loop increment +1 or -1
    i_r1 = FLOOR(MFF%CELLSI(FLOOR((r1_x-MFF%XS)*MFF%RDXINT)) + 1.0_EB) ! II start
    i_r2 = FLOOR(MFF%CELLSI(FLOOR((r2_x-MFF%XS)*MFF%RDXINT)) + 1.0_EB) ! II end
    j_r1 = FLOOR(MFF%CELLSJ(FLOOR((r1_y-MFF%YS)*MFF%RDYINT)) + 1.0_EB) ! JJ start
    j_r2 = FLOOR(MFF%CELLSJ(FLOOR((r2_y-MFF%YS)*MFF%RDYINT)) + 1.0_EB) ! JJ end
    i_r1 = MAX(1,MIN(i_r1,MFF%IBAR)) ! To be sure that indices are always ok.
    i_r2 = MAX(1,MIN(i_r2,MFF%IBAR))
    j_r1 = MAX(1,MIN(j_r1,MFF%JBAR))
    j_r2 = MAX(1,MIN(j_r2,MFF%JBAR))
    ! Same cell: sees always each other
    IF (ABS(i_r2-i_r1)+ABS(j_r2-j_r1) < 1) THEN
       ave_K = EVAC_MASS_EXTINCTION_COEFF*1.0E-6_EB*M%HUMAN_GRID(i_r1,j_r1)%SOOT_DENS
       max_fed = M%HUMAN_GRID(i_r1,j_r1)%FED_CO_CO2_O2
       RETURN
    END IF

    ! Choose the main direction:
    IF (ABS(r2_x-r1_x)/MFF%DX(1) < ABS(r2_y-r1_y)/MFF%DY(1)) THEN
       ! Now y is the main direction

       ! Add the last cell (r1) to the average
       ave_K = EVAC_MASS_EXTINCTION_COEFF*1.0E-6_EB*M%HUMAN_GRID(i_r1,j_r1)%SOOT_DENS/(ABS(j_r1-j_r2)+1)
       max_fed =  M%HUMAN_GRID(i_r1,j_r1)%FED_CO_CO2_O2

       i_old = i_r1 ; j_old = j_r1
       y = 0.0_EB
       MainLoopY: DO j = j_r1+isy, j_r2, isy
          y = y + isy*MFF%DY(j)
          x = MAX(MFF%XS,MIN(MFF%XF,r1_x + y*(r2_x - r1_x)/(r2_y - r1_y)))
          i = FLOOR(MFF%CELLSI(FLOOR((x-MFF%XS)*MFF%RDXINT)) + 1.0_EB)
          i = isx*MIN(isx*i_r2,isx*i) ! i in interval j_r1...j_r2
          ic  = MFF%CELL_INDEX(i_old,j_old,1)
          ic2 = MFF%CELL_INDEX(i    ,j_old,1) ! side cell
          iw  = MFF%WALL_INDEX(ic, isy*2) ! main direction
          iw1 = MFF%WALL_INDEX(ic ,isx*1) ! sideways
          iw2 = MFF%WALL_INDEX(ic2,isy*2) ! side + main direction
          IF (IW >0) THEN
             IF (MFF%OBSTRUCTION(MFF%WALL(IW )%OBST_INDEX)%HIDDEN .AND. MFF%WALL(IW )%OBST_INDEX>0) iw  = 0
          END IF
          IF (IW1>0) THEN
             IF (MFF%OBSTRUCTION(MFF%WALL(IW1)%OBST_INDEX)%HIDDEN .AND. MFF%WALL(IW1)%OBST_INDEX>0) iw1 = 0
          END IF
          IF (IW2>0) THEN
             IF (MFF%OBSTRUCTION(MFF%WALL(IW2)%OBST_INDEX)%HIDDEN .AND. MFF%WALL(IW2)%OBST_INDEX>0) iw2 = 0
          END IF
          ! iw is zero, if there is no solid boundary
          ! from (i,j)==>(i,jnew):    iw and iw2 are zero, iw1 does not matter
          !                           ic=ic2 ==> iw=iw2
          ! from (i,j)==>(inew,jnew): iw1 and iw2 are zero, iw does not matter
          IF ((i==i_old .AND. iw/=0) .OR. (i/=i_old .AND. (iw1/=0 .OR. iw2/=0))) THEN
             See_door = .FALSE.
             ave_K = EVAC_MASS_EXTINCTION_COEFF*1.0E-6_EB*M%HUMAN_GRID(i_r1,j_r1)%SOOT_DENS
             max_fed = M%HUMAN_GRID(i_r1,j_r1)%FED_CO_CO2_O2
             EXIT MainLoopY
          END IF
          ave_K = ave_K +  EVAC_MASS_EXTINCTION_COEFF*1.0E-6_EB*M%HUMAN_GRID(i,j)%SOOT_DENS/(ABS(j_r1-j_r2)+1)
          max_fed = MAX(max_fed, M%HUMAN_GRID(i,j)%FED_CO_CO2_O2)
          i_old = i ; j_old = j
       END DO MainLoopY

    ELSE
       ! Now x is the main direction
       ! Add the first cell (r1) to the average
       ave_K = EVAC_MASS_EXTINCTION_COEFF*1.0E-6_EB*M%HUMAN_GRID(i_r1,j_r1)%SOOT_DENS/(ABS(i_r1-i_r2)+1)
       max_fed =  M%HUMAN_GRID(i_r1,j_r1)%FED_CO_CO2_O2
       i_old = i_r1 ; j_old = j_r1
       x = 0.0_EB
       MainLoopX: DO i = i_r1+isx, i_r2, isx
          x = x + isx*MFF%DX(i)
          y = MAX(MFF%YS,MIN(MFF%YF,r1_y + x*(r2_y - r1_y)/(r2_x - r1_x)))
          j = FLOOR(MFF%CELLSJ(FLOOR((y-MFF%YS)*MFF%RDYINT)) + 1.0_EB)
          j = isy*MIN(isy*j_r2,isy*j) ! j in interval j_r1...j_r2
          ic  = MFF%CELL_INDEX(i_old,j_old,1)
          ic2 = MFF%CELL_INDEX(i_old,j    ,1) ! side cell
          iw  = MFF%WALL_INDEX(ic, isx*1) ! main direction
          iw1 = MFF%WALL_INDEX(ic ,isy*2) ! sideways
          iw2 = MFF%WALL_INDEX(ic2,isx*1) ! side + main direction
          IF (IW >0) THEN
             IF (MFF%OBSTRUCTION(MFF%WALL(IW )%OBST_INDEX)%HIDDEN .AND. MFF%WALL(IW )%OBST_INDEX>0) iw  = 0
          END IF
          IF (IW1>0) THEN
             IF (MFF%OBSTRUCTION(MFF%WALL(IW1)%OBST_INDEX)%HIDDEN .AND. MFF%WALL(IW1)%OBST_INDEX>0) iw1 = 0
          END IF
          IF (IW2>0) THEN
             IF (MFF%OBSTRUCTION(MFF%WALL(IW2)%OBST_INDEX)%HIDDEN .AND. MFF%WALL(IW2)%OBST_INDEX>0) iw2 = 0
          END IF
          ! iw is zero, if there is no solid boundary
          ! from (i,j)==>(inew,j):    iw and iw2 are zero, iw1 does not matter
          !                           ic=ic2 ==> iw=iw2
          ! from (i,j)==>(inew,jnew): iw1 and iw2 are zero, iw does not matter
          IF ((j==j_old .AND. iw/=0) .OR. (j/=j_old .AND. (iw1/=0 .OR. iw2/=0))) THEN
             See_door = .FALSE.
             ave_K = EVAC_MASS_EXTINCTION_COEFF*1.0E-6_EB*M%HUMAN_GRID(i_r1,j_r1)%SOOT_DENS
             max_fed = M%HUMAN_GRID(i_r1,j_r1)%FED_CO_CO2_O2
             EXIT MainLoopX
          END IF
          ave_K = ave_K + EVAC_MASS_EXTINCTION_COEFF*1.0E-6_EB*M%HUMAN_GRID(i,j)%SOOT_DENS/(ABS(i_r1-i_r2)+1)
          max_fed = MAX(max_fed, M%HUMAN_GRID(i,j)%FED_CO_CO2_O2)
          i_old = i ; j_old = j
       END DO MainLoopX
    END IF
  END FUNCTION See_door

  SUBROUTINE Find_walls(nm, r1_x, r1_y, r_circle, d_cutoff, Skip_Wall_Force_Ior, d_xy, FoundWall_xy, istat)
    IMPLICIT NONE
    !
    ! This subroutine checks if the circle is inside a solid (totally or partially inside).
    ! Only the positions of the closest walls in the four main directions are returned.
    ! If a wall is further away than d_cutoff from the arc of the circle then
    ! the corresponding edge of the mesh is returned.
    ! Inputs:  nm: mesh index
    !          (r1_x,r1_y): co-ordinates of the circle
    !          r_circle: radius of the circle
    !          d_cutoff: cut-off distance (e.g., P2P_DIST_MAX)
    !          Skip_Wall_Force_Ior: skip wall force in that direction
    ! Outputs: d_xy(1): wall position, -x direction
    !          d_xy(2): wall position, +x direction
    !          d_xy(3): wall position, -y direction
    !          d_xy(4): wall position, +y direction
    !          FoundWall_xy(1-4): True if found an obst (VEL>0 vent: False)
    !          istat: error or not?
    !            =-1: circle center is inside a solid ==> error
    !            =0 : circle is not inside a solid (only social forces)
    !            =1 : circle is partially inside a solid (contact forces also)
    ! NOTE: This works for thick and thin OBSTs.
    !
    ! Passed variables
    INTEGER, INTENT(IN) :: nm, Skip_Wall_Force_Ior
    REAL(EB), INTENT(IN) :: r1_x, r1_y, r_circle, d_cutoff
    REAL(EB), DIMENSION(4), INTENT(OUT) :: d_xy
    LOGICAL, DIMENSION(4), INTENT(OUT) :: FoundWall_xy
    INTEGER, INTENT(OUT) :: istat
    !
    ! Local variables
    INTEGER :: ii, jj, iw, ic, SURF_INDEX, is, i_end, iin, jjn, kkn, I_OBST
    REAL(EB) :: dx, dy, d_mx, d_px, d_my, d_py
    TYPE (MESH_TYPE), POINTER :: M =>NULL()

    M => MESHES(NM)
    istat = 0            ! Default
    FoundWall_xy = .FALSE.
    d_mx = M%XS
    d_px = M%XF
    d_my = M%YS
    d_py = M%YF

    ! Where is the circle
    IIN = FLOOR(M%CELLSI(FLOOR((r1_x-M%XS)*M%RDXINT))+1.0_EB)
    JJN = FLOOR(M%CELLSJ(FLOOR((r1_y-M%YS)*M%RDYINT))+1.0_EB)
    KKN = 1

    ! Find the closest wall at the -x direction

    is = -1    ! minus or plus direction
    IF (Skip_Wall_Force_Ior==-1) THEN
       d_mx = d_mx + is*(2.0_EB*d_cutoff)
    ELSE
       dx = M%DX(iin) ! no stretched meshes for evac
       i_end = INT((d_cutoff+r_circle)/dx) + 1
       i_end = iin + is*i_end
       i_end = MIN(M%IBAR,MAX(1,i_end))
       ii = iin
       ic = M%cell_index(ii,jjn,kkn)  ! cell index
       I_OBST = M%OBST_INDEX_C(IC)
       IF (M%Solid(ic) .AND. .NOT.M%OBSTRUCTION(I_OBST)%HIDDEN) THEN
          istat = -1
          RETURN
       END IF
       iw = M%wall_index(ic, is*1)      ! wall index
       IF (IW>0) THEN
          IF (M%OBSTRUCTION(M%WALL(IW)%OBST_INDEX)%HIDDEN .AND. M%WALL(IW)%OBST_INDEX>0) iw = 0
       END IF
       DO WHILE (iw==0 .AND. ii/=i_end)
          ii = ii + is
          ic = M%cell_index(ii,jjn,kkn)  ! cell index
          iw = M%wall_index(ic, is*1)      ! wall index
          IF (IW>0) THEN
             IF (M%OBSTRUCTION(M%WALL(IW)%OBST_INDEX)%HIDDEN .AND. M%WALL(IW)%OBST_INDEX>0) iw = 0
          END IF
       END DO

       IF (iw /= 0) THEN
          FoundWall_xy(1) = .TRUE.
          SURF_INDEX = M%WALL(IW)%SURF_INDEX ! Boundary condition index
          ! There is a 'door', i.e., outflow-boundary (or open boundary)
          ! so no wall forces ==> exit this loop
          d_mx = M%WALL(IW)%X
          I_OBST = M%OBST_INDEX_C(IC)
          IF (M%Solid(ic) .AND. .NOT.M%OBSTRUCTION(I_OBST)%HIDDEN) THEN
             WRITE(MESSAGE,'(A,I4,2I6)') 'ERROR: Find_Walls ',nm, ii,jjn
             CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
          END IF
          IF (SURFACE(SURF_INDEX)%VEL >= TWO_EPSILON_EB .OR. M%WALL(IW)%BOUNDARY_TYPE==OPEN_BOUNDARY) THEN
             !d_mx = d_mx + is*(2.0_EB*d_cutoff)
             FoundWall_xy(1) = .FALSE.
          END IF
          IF (ABS(d_mx-r1_x) < r_circle) istat = 1
       END IF
    END IF

    ! Find the closest wall at the +x direction

    is = +1    ! minus or plus direction
    IF (Skip_Wall_Force_Ior==1) THEN
       d_px = d_px + is*(2.0_EB*d_cutoff)
    ELSE
       dx = M%DX(iin) ! no stretched meshes for evac
       i_end = INT((d_cutoff+r_circle)/dx) + 1
       i_end = iin + is*i_end
       i_end = MIN(M%IBAR,MAX(1,i_end))
       ii = iin
       ic = M%cell_index(ii,jjn,kkn)  ! cell index
       I_OBST = M%OBST_INDEX_C(IC)
       IF (M%Solid(ic) .AND. .NOT.M%OBSTRUCTION(I_OBST)%HIDDEN) THEN
          istat = -1
          RETURN
       END IF
       iw = M%wall_index(ic, is*1)      ! wall index
       IF (IW>0) THEN
          IF (M%OBSTRUCTION(M%WALL(IW)%OBST_INDEX)%HIDDEN .AND. M%WALL(IW)%OBST_INDEX>0) iw = 0
       END IF
       DO WHILE (iw==0 .AND. ii/=i_end)
          ii = ii + is
          ic = M%cell_index(ii,jjn,kkn)  ! cell index
          iw = M%wall_index(ic, is*1)      ! wall index
          IF (IW>0) THEN
             IF (M%OBSTRUCTION(M%WALL(IW)%OBST_INDEX)%HIDDEN .AND. M%WALL(IW)%OBST_INDEX>0) iw = 0
          END IF
       END DO

       IF (iw /= 0) THEN
          FoundWall_xy(2) = .TRUE.
          SURF_INDEX = M%WALL(IW)%SURF_INDEX ! Boundary condition index
          ! There is a 'door', i.e., outflow-boundary (or open boundary)
          ! so no wall forces ==> exit this loop
          d_px = M%WALL(IW)%X
          I_OBST = M%OBST_INDEX_C(IC)
          IF (M%Solid(ic) .AND. .NOT.M%OBSTRUCTION(I_OBST)%HIDDEN) THEN
             WRITE(MESSAGE,'(A,I4,2I6)') 'ERROR: Find_Walls ',nm, ii,jjn
             CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
          END IF
          IF (SURFACE(SURF_INDEX)%VEL > TWO_EPSILON_EB .OR. M%WALL(IW)%BOUNDARY_TYPE==OPEN_BOUNDARY) THEN
             !d_px = d_px + is*(2.0_EB*d_cutoff)
             FoundWall_xy(2) = .FALSE.
          END IF
          IF (ABS(d_px-r1_x) < r_circle) istat = 1
       END IF
    END IF

    ! Find the closest wall at the -y direction

    is = -1    ! minus or plus direction
    IF (Skip_Wall_Force_Ior==-2) THEN
       d_my = d_my + is*(2.0_EB*d_cutoff)
    ELSE
       dy = M%DY(jjn) ! no stretched meshes for evac
       i_end = INT((d_cutoff+r_circle)/dy) + 1
       i_end = jjn + is*i_end
       i_end = MIN(M%JBAR,MAX(1,i_end))
       jj = jjn
       ic = M%cell_index(iin,jj,kkn)  ! cell index
       I_OBST = M%OBST_INDEX_C(IC)
       IF (M%Solid(ic) .AND. .NOT.M%OBSTRUCTION(I_OBST)%HIDDEN) THEN
          istat = -1
          RETURN
       END IF
       iw = M%wall_index(ic, is*2)      ! wall index
       IF (IW>0) THEN
          IF (M%OBSTRUCTION(M%WALL(IW)%OBST_INDEX)%HIDDEN .AND. M%WALL(IW)%OBST_INDEX>0) iw = 0
       END IF
       DO WHILE (iw==0 .AND. jj/=i_end)
          jj = jj + is
          ic = M%cell_index(iin,jj,kkn)  ! cell index
          iw = M%wall_index(ic, is*2)      ! wall index
          IF (IW>0) THEN
             IF (M%OBSTRUCTION(M%WALL(IW)%OBST_INDEX)%HIDDEN .AND. M%WALL(IW)%OBST_INDEX>0) iw = 0
          END IF
       END DO

       IF (iw /= 0) THEN
          FoundWall_xy(3) = .TRUE.
          SURF_INDEX = M%WALL(IW)%SURF_INDEX ! Boundary condition index
          ! There is a 'door', i.e., outflow-boundary (or open boundary)
          ! so no wall forces ==> exit this loop
          d_my = M%WALL(IW)%Y
          I_OBST = M%OBST_INDEX_C(IC)
          IF (M%Solid(ic) .AND. .NOT.M%OBSTRUCTION(I_OBST)%HIDDEN) THEN
             WRITE(MESSAGE,'(A,I4,2I6)') 'ERROR: Find_Walls ',nm, ii,jjn
             CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
          END IF
          IF (SURFACE(SURF_INDEX)%VEL >= TWO_EPSILON_EB .OR. M%WALL(IW)%BOUNDARY_TYPE==OPEN_BOUNDARY) THEN
             !d_my = d_my + is*(2.0_EB*d_cutoff)
             FoundWall_xy(3) = .FALSE.
          END IF
          IF (ABS(d_my-r1_y) < r_circle) istat = 1
       END IF
    END IF

    ! Find the closest wall at the +y direction

    is = +1    ! minus or plus direction
    IF (Skip_Wall_Force_Ior==2) THEN
       d_py = d_py + is*(2.0_EB*d_cutoff)
    ELSE
       dy = M%DY(jjn) ! no stretched meshes for evac
       i_end = INT((d_cutoff+r_circle)/dy) + 1
       i_end = jjn + is*i_end
       i_end = MIN(M%JBAR,MAX(1,i_end))
       jj = jjn
       ic = M%cell_index(iin,jj,kkn)  ! cell index
       I_OBST = M%OBST_INDEX_C(IC)
       IF (M%Solid(ic) .AND. .NOT.M%OBSTRUCTION(I_OBST)%HIDDEN) THEN
          istat = -1
          RETURN
       END IF
       iw = M%wall_index(ic, is*2)      ! wall index
       IF (IW>0) THEN
          IF (M%OBSTRUCTION(M%WALL(IW)%OBST_INDEX)%HIDDEN .AND. M%WALL(IW)%OBST_INDEX>0) iw = 0
       END IF
       DO WHILE (iw==0 .AND. jj/=i_end)
          jj = jj + is
          ic = M%cell_index(iin,jj,kkn)  ! cell index
          iw = M%wall_index(ic, is*2)      ! wall index
          IF (IW>0) THEN
             IF (M%OBSTRUCTION(M%WALL(IW)%OBST_INDEX)%HIDDEN .AND. M%WALL(IW)%OBST_INDEX>0) iw = 0
          END IF
       END DO

       IF (iw /= 0) THEN
          FoundWall_xy(4) = .TRUE.
          SURF_INDEX = M%WALL(IW)%SURF_INDEX  ! Boundary condition index
          ! There is a 'door', i.e., outflow-boundary (or open boundary)
          ! so no wall forces ==> exit this loop
          d_py = M%WALL(IW)%Y
          I_OBST = M%OBST_INDEX_C(IC)
          IF (M%Solid(ic) .AND. .NOT.M%OBSTRUCTION(I_OBST)%HIDDEN) THEN
             WRITE(MESSAGE,'(A,I4,2I6)') 'ERROR: Find_Walls ',nm, ii,jjn
             CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.) ; RETURN
          END IF
          IF (SURFACE(SURF_INDEX)%VEL > TWO_EPSILON_EB .OR. M%WALL(IW)%BOUNDARY_TYPE==OPEN_BOUNDARY) THEN
             !d_py = d_py + is*(2.0_EB*d_cutoff)
             FoundWall_xy(4) = .FALSE.
          END IF
          IF (ABS(d_py-r1_y) < r_circle) istat = 1
       END IF
    END IF
    d_xy(1) = d_mx
    d_xy(2) = d_px
    d_xy(3) = d_my
    d_xy(4) = d_py
    RETURN
  END SUBROUTINE Find_walls

  SUBROUTINE GET_FIRE_CONDITIONS(NOM,I,J,K,fed_indx,soot_dens,gas_temp,rad_flux, ZZ_GET, FED_ACTIVITY)
    IMPLICIT NONE
    !
    ! Passed variables
    INTEGER, INTENT(IN) :: I, J, K, NOM, FED_ACTIVITY
    REAL(EB), INTENT(OUT) :: fed_indx, soot_dens, gas_temp, rad_flux
    REAL(EB), INTENT(INOUT) :: ZZ_GET(1:N_TRACKED_SPECIES)
    !
    ! Local variables
    REAL(EB) :: Y_MF_INT

    soot_dens = 0._EB ; fed_indx = 0._EB ; gas_temp = 0._EB ; rad_flux = 0._EB
    ! Mass fraction array ==> soot density (mg/m3)
    ! Next is for soot (mg/m3)
    ZZ_GET(1:N_TRACKED_SPECIES) = MESHES(nom)%ZZ(I,J,K,1:N_TRACKED_SPECIES)
    IF (SOOT_INDEX > 0) THEN
       CALL GET_MASS_FRACTION(ZZ_GET,SOOT_INDEX,Y_MF_INT)
       soot_dens = Y_MF_INT*MESHES(nom)%RHO(I,J,K)*1.E6_EB
    ELSE
       soot_dens = 0._EB
    ENDIF
    ! Calculate Purser's fractional effective dose (FED)
    fed_indx = FED(ZZ_GET,MESHES(nom)%RSUM(I,J,K),FED_ACTIVITY)
    ! Gas temperature, ind=5, C
    gas_temp  = MESHES(nom)%TMP(I,J,K)
    ! Rad flux, ind=18, kW/m2 (no -sigma*Tamb^4 term)
    rad_flux = MAX(MESHES(nom)%UII(I,J,K)/4.0_EB,SIGMA*TMPA4)

  END SUBROUTINE GET_FIRE_CONDITIONS

  SUBROUTINE Change_Target_Door(nm, nm2, ie, j, j1, i_egrid, imode, xx, yy, I_Target, I_Color, I_Field, HR)
    IMPLICIT NONE
    !
    ! Door selection algorithm
    !
    ! Inputs:
    !   nm: mesh index
    !   nm2: mesh2 index
    !   ie: the index of the agent on this mesh  NOT USED
    !   j:  group index
    !   j1: lonely agent index
    !   i_egrid: main evacuation mesh index
    !   imode: 0= initialization, 1= evacuate_humans, 2= check_target_node
    !   xx: x co-ordinate of the agent
    !   yy: y co-ordinate of the agent
    !
    ! Outputs:
    !   I_Target: Target door/exit node index
    !   I_Color: color index
    !   I_Field: new flow field index
    !
    ! Passed variables
    INTEGER, INTENT(IN) :: nm, nm2, ie, j, j1, i_egrid, imode
    REAL(EB), INTENT(IN) :: xx, yy
    INTEGER, INTENT(OUT)  :: I_Target, I_Color
    INTEGER, INTENT(INOUT)  :: I_Field
    TYPE (HUMAN_TYPE), POINTER :: HR
    !
    ! Local variables
    REAL :: RN_REAL
    REAL(EB) :: L2_min, max_fed, ave_K, L2_tmp, RN, max_fed2, ave_K2, L2_tmp2
    REAL(EB) :: x1_old, y1_old, Speed, X11, Y11, x_o, y_o, XBx, XBy
    REAL(EB) :: X1, Y1, X2, Y2, DOOR_WIDTH, X_XYZ, Y_XYZ
    INTEGER :: i_old_ffield, i_tmp, i_new_ffield, IEL, color_index, DOOR_IOR
    INTEGER :: i, i_o, izero, nm_tmp, I_Agent_Type, I_Old_Target
    CHARACTER(LABEL_LENGTH) :: name_old_ffield, name_new_ffield
    LOGICAL :: PP_see_door, PP_see_doorXB, PP_correct_side
    REAL(EB) :: T_tmp, T_tmp1, Width
    INTEGER :: N_queue, II, JJ, JJ_NOW
    TYPE (EVACUATION_TYPE), POINTER :: HPT =>NULL()
    TYPE (EVAC_ENTR_TYPE),  POINTER :: PNX =>NULL()


    IF (imode < 2) THEN
       nm_tmp = nm
    ELSE
       ! This is for check_target_node
       nm_tmp = nm2
    END IF

    I_Agent_Type = HR%I_DoorAlgo
    i_old_ffield = HR%I_FFIELD
    I_Old_Target = HR%I_Target
    IF (i_old_ffield > 0) THEN
       name_old_ffield = TRIM(MESH_NAME(i_old_ffield))
    ELSE
       name_old_ffield = TRIM(MESH_NAME(nm_tmp))
    END IF
    IF (imode==1 .AND. I_Agent_Type==3 .AND. HR%I_Target/=0 .AND. (I_HERDING_TYPE==1 .OR. I_HERDING_TYPE==3)) THEN
       ! These herding agents keep the first choise
       I_Target = HR%I_Target
       I_Color = HR%COLOR_INDEX
       RETURN
    END IF
    IF (HR%IEL > 0 ) THEN
       ! Agent HR originates from an evac line
       HPT => EVAC_EVACS(HR%IEL)
    ELSE
       ! Agent HR originates from an entr line (and is a lonely agent)
       PNX => EVAC_ENTRYS(ABS(HR%IEL))
    END IF

    ! Only those doors are possible which are in the same main evac mesh.
    K_ave_door        = 0.0000001_EB*ABS(FED_DOOR_CRIT)
    FED_max_Door      = 0.0_EB
    Is_Known_Door     = .FALSE.
    Is_Visible_Door   = .FALSE.
    Is_BeeLine_Door   = .FALSE.
    Is_BeeLineXB_Door = .FALSE.
    ! Group index=0: the agent is from an entry line (no evac line)
    IF ( HR%GROUP_ID /= 0 ) THEN
       IF ( HR%GROUP_ID < 0 ) THEN
          ! A lonely soul
          x1_old = xx
          y1_old = yy
          IEL    = HR%IEL
          Speed  = HR%Speed
          DO i = 1, N_DOORS
             IF ( EVAC_DOORS(i)%IMESH == nm_tmp .AND. .NOT.(TRIM(EVAC_DOORS(i)%TO_NODE)==TRIM(EVAC_DOORS(i)%ID))) THEN
                Is_Visible_Door(i) = .TRUE.
                IF (EVAC_DOORS(i)%KNOWN_DOOR) Is_Known_Door(i) = .TRUE.
                IF (imode == 0) CYCLE   ! Initialization call
                DO i_tmp = 1, Human_Known_Doors(j1)%N_nodes
                   IF (EVAC_DOORS(i)%INODE == ABS(Human_Known_Doors(j1)%I_nodes(i_tmp))) Is_Known_Door(i) = .TRUE.
                END DO
             END IF
          END DO
          DO i = 1, N_EXITS
             IF ( EVAC_EXITS(i)%IMESH == nm_tmp .AND. .NOT. EVAC_EXITS(i)%COUNT_ONLY ) THEN
                Is_Visible_Door(N_DOORS+i) = .TRUE.
                IF (EVAC_EXITS(i)%KNOWN_DOOR) Is_Known_Door(N_DOORS+i) = .TRUE.
                IF (imode == 0) CYCLE   ! Initialization call
                DO i_tmp = 1, Human_Known_Doors(j1)%N_nodes
                   IF (EVAC_EXITS(i)%INODE == ABS(Human_Known_Doors(j1)%I_nodes(i_tmp))) Is_Known_Door(N_DOORS+i) = .TRUE.
                END DO
             END IF
          END DO
       ELSE
          ! A member of a group
          x1_old = Group_List(j)%GROUP_X
          y1_old = Group_List(j)%GROUP_Y
          IEL    = Group_List(j)%IEL
          Speed  = Group_List(j)%Speed
          DO i = 1, N_DOORS
             IF ( EVAC_DOORS(i)%IMESH == nm_tmp .AND. .NOT.(TRIM(EVAC_DOORS(i)%TO_NODE)==TRIM(EVAC_DOORS(i)%ID))) THEN
                Is_Visible_Door(i) = .TRUE.
                IF (EVAC_DOORS(i)%KNOWN_DOOR) Is_Known_Door(i) = .TRUE.
                IF (imode == 0) CYCLE   ! Initialization call
                DO i_tmp = 1, Group_Known_Doors(j)%N_nodes
                   IF (EVAC_DOORS(i)%INODE == ABS(Group_Known_Doors(j)%I_nodes(i_tmp))) Is_Known_Door(i) = .TRUE.
                END DO
             END IF
          END DO
          DO i = 1, N_EXITS
             IF ( EVAC_EXITS(i)%IMESH == nm_tmp .AND. .NOT. EVAC_EXITS(i)%COUNT_ONLY ) THEN
                Is_Visible_Door(N_DOORS+i) = .TRUE.
                IF (EVAC_EXITS(i)%KNOWN_DOOR) Is_Known_Door(N_DOORS+i) = .TRUE.
                IF (imode == 0) CYCLE   ! Initialization call
                DO i_tmp = 1, Group_Known_Doors(j)%N_nodes
                   IF (EVAC_EXITS(i)%INODE == ABS(Group_Known_Doors(j)%I_nodes(i_tmp))) Is_Known_Door(N_DOORS+i) = .TRUE.
                END DO
             END IF
          END DO
       END IF

       IF (imode == 0) THEN
          ! Initialization phase: Draw the known doors/exits.
          Known_Doors_Loop: DO i = 1, HPT%N_VENT_FFIELDS
             i_tmp = 1
             IF (TRIM(EVAC_Node_List(HPT%I_DOOR_NODES(i))%Node_Type) == 'Door') THEN
                i_tmp = EVAC_Node_List(HPT%I_DOOR_NODES(i))%Node_Index
                IF (TRIM(EVAC_DOORS(i_tmp)%ID) == TRIM(EVAC_DOORS(i_tmp)%TO_NODE)) THEN
                   Is_Known_Door(i_tmp) = .FALSE.  ! An entry type door
                   CYCLE Known_Doors_Loop
                END IF
             END IF
             IF (TRIM(EVAC_Node_List(HPT%I_DOOR_NODES(i))%Node_Type) == 'Exit' ) THEN
                i_tmp = N_DOORS + EVAC_Node_List(HPT%I_DOOR_NODES(i))%Node_Index
             END IF
             IF (HPT%P_VENT_FFIELDS(i) < 1.0_EB) THEN
                CALL RANDOM_NUMBER(RN_REAL)
                RN = REAL(RN_REAL,EB)
                IF ( RN < HPT%P_VENT_FFIELDS(i) ) THEN
                   Is_Known_Door(i_tmp) = .TRUE.
                ELSE
                   Is_Known_Door(i_tmp) = .FALSE.
                END IF
             ELSE
                Is_Known_Door(i_tmp) = .TRUE.
             END IF
          END DO Known_Doors_Loop

          ! Save the random known door information
          i_tmp = COUNT(Is_Known_Door)
          IF (j > 0) THEN    ! The agent is a member of a group
             Group_Known_Doors(j)%N_nodes = i_tmp
             IF (COUNT(Is_Known_Door) > 0) THEN
                ALLOCATE(Group_Known_Doors(j)%I_nodes(i_tmp), STAT=IZERO)
                CALL ChkMemErr('Change_Target_Door', 'Group_Known_Doors',IZERO)
                i_tmp = 0
                DO i = 1, N_DOORS
                   IF (Is_Known_Door(i)) THEN
                      i_tmp = i_tmp + 1
                      Group_Known_Doors(j)%I_nodes(i_tmp) = EVAC_DOORS(i)%INODE
                   END IF
                END DO
                DO i = 1, N_EXITS
                   IF (Is_Known_Door(i+N_DOORS)) THEN
                      i_tmp = i_tmp + 1
                      Group_Known_Doors(j)%I_nodes(i_tmp) = EVAC_EXITS(i)%INODE
                   END IF
                END DO
             END IF        ! there are known doors for this group
          ELSE    ! The agent is a lonely one
             Human_Known_Doors(j1)%N_nodes = i_tmp
             IF (COUNT(Is_Known_Door) > 0) THEN
                ALLOCATE(Human_Known_Doors(j1)%I_nodes(i_tmp),STAT=IZERO)
                CALL ChkMemErr('Change_Target_Door', 'Human_Known_Doors',IZERO)
                i_tmp = 0
                DO i = 1, N_DOORS
                   IF (Is_Known_Door(i)) THEN
                      i_tmp = i_tmp + 1
                      Human_Known_Doors(j1)%I_nodes(i_tmp) = EVAC_DOORS(i)%INODE
                   END IF
                END DO
                DO i = 1, N_EXITS
                   IF (Is_Known_Door(i+N_DOORS)) THEN
                      i_tmp = i_tmp + 1
                      Human_Known_Doors(j1)%I_nodes(i_tmp) = EVAC_EXITS(i)%INODE
                   END IF
                END DO
             END IF   ! Is there any known doors for this group?
          END IF    ! Is the agent a member of a group or not?

          ! Check that the door is in correct main evac mesh.
          DO i = 1, N_DOORS
             IF ( EVAC_DOORS(i)%IMESH /= HPT%imesh ) THEN
                Is_Known_Door(i) = .FALSE.
             END IF
          END DO
          DO i = 1, N_EXITS
             IF ( EVAC_EXITS(i)%IMESH /= HPT%imesh .OR. EVAC_EXITS(i)%COUNT_ONLY ) THEN
                Is_Known_Door(N_DOORS+i) = .FALSE.
             END IF
          END DO
          HR%I_Door_Mode = 0 ! Initialization
       END IF  ! imode=0, Initialization call

    ELSE
       ! The agent is from an entry line. P_door= 0.0 or 1.0, known doors
       x1_old = xx
       y1_old = yy
       IEL    = ABS(HR%IEL)
       Speed  = HR%Speed
       DO i = 1, N_DOORS
          IF ( EVAC_DOORS(i)%IMESH == nm_tmp .AND. .NOT.(TRIM(EVAC_DOORS(i)%TO_NODE)==TRIM(EVAC_DOORS(i)%ID))) THEN
             Is_Visible_Door(i) = .TRUE.
          END IF
       END DO
       DO i = 1, N_EXITS
          IF ( EVAC_EXITS(i)%IMESH == nm_tmp .AND. .NOT. EVAC_EXITS(i)%COUNT_ONLY ) THEN
             Is_Visible_Door(N_DOORS+i) = .TRUE.
          END IF
       END DO
       Known_Doors_Loop2: DO i = 1, PNX%N_VENT_FFIELDS
          ! P = 0 or 1 for entrys.
          ! Check that the door/exit is in the correct main evac grid.
          i_tmp = 1
          IF (TRIM(EVAC_Node_List(PNX%I_DOOR_NODES(i))%Node_Type) == 'Door') THEN
             i_tmp = EVAC_Node_List(PNX%I_DOOR_NODES(i))%Node_Index
             IF (TRIM(EVAC_DOORS(i_tmp)%ID) == TRIM(EVAC_DOORS(i_tmp)%TO_NODE)) THEN
                Is_Known_Door(i_tmp) = .FALSE.  ! An entry type door
                CYCLE Known_Doors_Loop2
             END IF
          END IF
          IF (TRIM(EVAC_Node_List(PNX%I_DOOR_NODES(i))%Node_Type) == 'Exit' ) THEN
             i_tmp = N_DOORS + EVAC_Node_List(PNX%I_DOOR_NODES(i))%Node_Index
          END IF
          IF ( Is_Visible_Door(i_tmp) ) THEN
             ! Door/exit is on this main evac mesh.
             IF (PNX%P_VENT_FFIELDS(i) < 0.5_EB) THEN
                Is_Known_Door(i_tmp) = .FALSE.
             ELSE
                Is_Known_Door(i_tmp) = .TRUE.
             END IF
          END IF
       END DO Known_Doors_Loop2
       DO i = 1, N_DOORS
          IF (Is_Visible_Door(i) .AND. EVAC_DOORS(i)%KNOWN_DOOR) Is_Known_Door(i) = .TRUE.
       END DO
       DO i = 1, N_EXITS
          IF (Is_Visible_Door(N_DOORS+i) .AND. EVAC_EXITS(i)%KNOWN_DOOR)&
               Is_Known_Door(N_DOORS+i) = .TRUE.
       END DO
    END IF   ! Is the agent from an entr or from an evac line
    ! Now Is_Visible_Door means that on the same floor.
    ! Now Is_Know_Door means: known + correct floor

    ! Check that the known/visible doors are in the correct evacuation zone.
    ! The U_EVAC and V_EVAC is set to zero outside the zone of the door/exit.
    ! Agent is at the cell II,JJ, check the flow field at this cell.
    ! This check is done for actual evacuation movement. The initial selection of the
    ! doors might not be nice, if one has many evacuation zones in one evacuation mesh.
    ! (At the initialization phase, the u_evac,v_evac information is not up to date.)
    IF (ICYC>0) THEN
       II = FLOOR(CELLSI(FLOOR((HR%X-XS)*RDXINT)) + 1.0_EB)
       JJ = FLOOR(CELLSJ(FLOOR((HR%Y-YS)*RDYINT)) + 1.0_EB)
       DO I = 1, N_DOORS
          IF (.NOT.(Is_Visible_Door(i))) CYCLE
          JJ_NOW = EVAC_DOORS(I)%I_EMESH_EXITS
          IF (ABS(EMESH_EXITS(JJ_NOW)%U_EVAC(II,JJ)**2+EMESH_EXITS(JJ_NOW)%V_EVAC(II,JJ)**2)==0.0_EB) THEN
             Is_Known_Door(I) = .FALSE.
             Is_Visible_Door(I) = .FALSE.
          END IF
       END DO
       DO I = 1, N_EXITS
          IF(EVAC_EXITS(I)%COUNT_ONLY) CYCLE
          IF (.NOT.(Is_Visible_Door(N_DOORS+i))) CYCLE
          JJ_NOW = EVAC_EXITS(I)%I_EMESH_EXITS
          IF (ABS(EMESH_EXITS(JJ_NOW)%U_EVAC(II,JJ)**2+EMESH_EXITS(JJ_NOW)%V_EVAC(II,JJ)**2)==0.0_EB) THEN
             Is_Known_Door(N_DOORS+I) = .FALSE.
             Is_Visible_Door(N_DOORS+I) = .FALSE.
          END IF
       END DO
    END IF

    ! Agent types: 1 rational agents, 2 known doors, 3 herding, 0 main evac ff
    IF (I_Agent_Type == 0) THEN  ! main evac ff (or evac/entr line ff)
       Is_Known_Door     = .FALSE.
       Is_Visible_Door   = .FALSE.
       Is_BeeLine_Door   = .FALSE.
       Is_BeeLineXB_Door = .FALSE.
    END IF
    ! Agent type 4 follow others if they can and find the door to better than their own choise
    IF (I_Agent_Type == 4 .AND. HR%I_Target2 /= 0) THEN
       IF (Is_Visible_Door(ABS(HR%I_Target2))) Is_Known_Door(ABS(HR%I_Target2)) = .TRUE.
    END IF

    DO i = 1, N_DOORS
       IF (ABS(EVAC_DOORS(I)%IMODE)==2 .AND. .NOT.EVAC_DOORS(I)%TARGET_WHEN_CLOSED) THEN
          Is_Visible_Door(i) = .FALSE.
          Is_Known_Door(i) = .FALSE.
       END IF
    END DO
    DO i = 1, N_EXITS
       IF ( (ABS(EVAC_EXITS(i)%IMODE)==2) .AND. .NOT. EVAC_EXITS(i)%COUNT_ONLY .AND. &
            .NOT.EVAC_EXITS(I)%TARGET_WHEN_CLOSED) THEN
          Is_Visible_Door(N_DOORS+i) = .FALSE.
          Is_Known_Door(N_DOORS+i) = .FALSE.
       END IF
    END DO

    ! Find the visible doors.
    DO i = 1, N_DOORS + N_EXITS
       PP_correct_side = .TRUE.  ! If the door xb is seen on the correct side of the door
       IF ( Is_Visible_Door(i) ) THEN

          IF (EVAC_Node_List(n_egrids+N_ENTRYS+i)%Node_Type == 'Door' ) THEN
             X_XYZ = EVAC_DOORS(i)%X ; Y_XYZ = EVAC_DOORS(i)%Y
             X1 = EVAC_DOORS(i)%X1 ; X2 = EVAC_DOORS(i)%X2
             Y1 = EVAC_DOORS(i)%Y1 ; Y2 = EVAC_DOORS(i)%Y2
             DOOR_IOR = EVAC_DOORS(i)%IOR
             DOOR_WIDTH = EVAC_DOORS(i)%WIDTH
          ELSE
             X_XYZ = EVAC_EXITS(i-N_DOORS)%X ; Y_XYZ = EVAC_EXITS(i-N_DOORS)%Y
             X1 = EVAC_EXITS(i-N_DOORS)%X1 ; X2 = EVAC_EXITS(i-N_DOORS)%X2
             Y1 = EVAC_EXITS(i-N_DOORS)%Y1 ; Y2 = EVAC_EXITS(i-N_DOORS)%Y2
             DOOR_IOR = EVAC_EXITS(i-N_DOORS)%IOR
             DOOR_WIDTH = EVAC_EXITS(i-N_DOORS)%WIDTH
          END IF
          XBx = 0.5_EB*(X1+X2) ; XBy = 0.5_EB*(Y1+Y2)  ! Mid point of XB
          ! Find the door post positions
          IF (ABS(DOOR_IOR)==1)THEN
             Y1 = Y1+MIN(0.3_EB,0.5_EB*DOOR_WIDTH) ; Y2 = Y2-MIN(0.3_EB,0.5_EB*DOOR_WIDTH)
          ELSE
             X1 = X1+MIN(0.3_EB,0.5_EB*DOOR_WIDTH) ; X2 = X2-MIN(0.3_EB,0.5_EB*DOOR_WIDTH)
          END IF
          ! Check if the agent is on the correct side of the door (thin obsts)
          IF (ABS(DOOR_IOR)==1)THEN
             IF ( DOOR_IOR*(x1_old-XBx)>0.0_EB ) PP_correct_side = .FALSE.
          ELSE
             IF ( DOOR_IOR*(y1_old-XBy)>0.0_EB ) PP_correct_side = .FALSE.
          END IF

          PP_see_door = See_each_other(nm_tmp, x1_old, y1_old, x1, y1)
          Is_BeeLineXB_Door(i,2) = PP_see_door .AND. PP_correct_side
          PP_see_door = See_each_other(nm_tmp, x1_old, y1_old, x2, y2)
          Is_BeeLineXB_Door(i,3) = PP_see_door .AND. PP_correct_side
          IF (ABS(DOOR_IOR)==1)THEN
             Y1 = Y_XYZ - 0.5_EB*DOOR_WIDTH + MIN(0.3_EB,0.5_EB*DOOR_WIDTH)
             Y2 = Y_XYZ + 0.5_EB*DOOR_WIDTH - MIN(0.3_EB,0.5_EB*DOOR_WIDTH)
             X1 = X_XYZ ; X2 = X_XYZ
          ELSE
             X1 = X_XYZ - 0.5_EB*DOOR_WIDTH + MIN(0.3_EB,0.5_EB*DOOR_WIDTH)
             X2 = X_XYZ + 0.5_EB*DOOR_WIDTH - MIN(0.3_EB,0.5_EB*DOOR_WIDTH)
             Y1 = Y_XYZ ; Y2 = Y_XYZ
          END IF
          PP_see_door = See_each_other(nm_tmp, x1_old, y1_old, x1, y1)
          Is_BeeLine_Door(i,2) = PP_see_door
          PP_see_door = See_each_other(nm_tmp, x1_old, y1_old, x2, y2)
          Is_BeeLine_Door(i,3) = PP_see_door

          IF (EVAC_Node_List(n_egrids+N_ENTRYS+i)%Node_Type == 'Door' ) THEN
             X11 = 0.5_EB*(EVAC_DOORS(i)%X1 + EVAC_DOORS(i)%X2)
             Y11 = 0.5_EB*(EVAC_DOORS(i)%Y1 + EVAC_DOORS(i)%Y2)
          ELSE
             X11 = 0.5_EB*(EVAC_EXITS(i-N_DOORS)%X1 + EVAC_EXITS(i-N_DOORS)%X2)
             Y11 = 0.5_EB*(EVAC_EXITS(i-N_DOORS)%Y1 + EVAC_EXITS(i-N_DOORS)%Y2)
          END IF
          ! Groups: the first member (x1_old,y1_old) of the group is used.
          PP_see_door = See_each_other(nm_tmp, x1_old, y1_old, x11, y11)
          Is_BeeLineXB_Door(i,1) = PP_see_door .AND. PP_correct_side
          IF (EVAC_Node_List(n_egrids+N_ENTRYS+i)%Node_Type == 'Door' ) THEN
             X11 = EVAC_DOORS(i)%X
             Y11 = EVAC_DOORS(i)%Y
          ELSE
             X11 = EVAC_EXITS(i-N_DOORS)%X
             Y11 = EVAC_EXITS(i-N_DOORS)%Y
          END IF
          PP_see_door = See_each_other(nm_tmp, x1_old, y1_old, x11, y11)
          Is_BeeLine_Door(i,1) = PP_see_door
          PP_see_door   = See_door(nm_tmp, x1_old, y1_old, x11, y11, ave_K, max_fed)
          PP_see_doorXB = See_door(nm_tmp, x1_old, y1_old, XBx, XBy, ave_K2, max_fed2) .AND. PP_correct_side
          IF (PP_see_doorXB .AND. .NOT.PP_see_door) THEN
             max_fed = max_fed2 ; ave_K = ave_K2
          END IF
          IF (PP_see_doorXB .AND. PP_see_door) THEN
             max_fed = MIN(max_fed, max_fed2) ; ave_K = MIN(ave_K, ave_K2)
          END IF
          PP_see_door = PP_see_door .OR. (PP_see_doorXB .AND. PP_correct_side)
             
          FED_max_Door(i) = max_fed
          K_ave_Door(i) = MAX(ave_K,0.5_EB*ABS(FED_DOOR_CRIT)) ! no divisions by zero

          ! Note: a DOOR is not counted as visible door, if it does not have an
          ! EXIT_SIGN, unless it is already been a target door for this agent/group.
          IF (PP_see_door) THEN
             IF (EVAC_Node_List(n_egrids+N_ENTRYS+i)%Node_Type == 'Door') THEN
                IF (.NOT. EVAC_DOORS(i)%EXIT_SIGN .AND. .NOT. ABS(HR%I_Target) == i .AND. &
                    .NOT. Is_Known_Door(i)) THEN
                     ! no exit sign, not the current target door, not known
                   Is_Visible_Door(i) = .FALSE.
                END IF
             END IF
          ELSE
             Is_Visible_Door(i) = .FALSE.
          END IF
          ! Rational Agents know all visible doors (and the given known doors)
          IF (I_Agent_Type==1 .AND. Is_Visible_Door(i)) Is_Known_Door(i) = .TRUE.
          ! Herding agents do not use visible doors if they are not known doors
          IF (I_Agent_Type==3 .AND. .NOT.Is_Known_Door(i)) Is_Visible_Door(i) = .FALSE.
          ! I_Agent_Type==4: These agents follow others, but if there is nobody to follow
          ! then they behave like the agents I_Agent_Type==2, i.e., normal door choise.
       END IF ! correct main evac mesh
    END DO ! all doors and exits

    ! Now: Is_Visible_Door correct floor + is visible
    ! Now: Is_Known_Door correct floor + known door

    ! If the current target door is visible then it is considered to be known.
    ! If there would be other known doors then this door would not be chosen if
    ! it was not known door.  And if there are no known doors at all then the
    ! one that have already be chosen is better known that the others.
    ! Note: I_Target < 0: not visible, >0: visible

    ! If the target door has been visible before, it is considered to be visible still.
    IF (HR%I_Target > 0) Is_Visible_Door(HR%I_Target) = .TRUE.

    IF (imode == 1) THEN
       DO i = 1, N_DOORS + N_EXITS
          IF (ABS(HR%I_Target) == i .AND. Is_Visible_Door(i)) Is_Known_Door(i) = .TRUE.
          IF (HR%GROUP_ID < 0) THEN ! This does not work for groups yet nor ENTR agents
             ! If the door is known to have too much smoke, set it unknown and not visible
             ! If the door is known to have some smoke, set it unknown
             Do_Known_doors0: DO ii=1,Human_Known_Doors(j1)%N_nodes
                IF (ABS(Human_Known_Doors(j1)%I_nodes(ii)) == n_egrids+N_ENTRYS+i) THEN
                   IF (Human_Known_Doors(j1)%I_nodes(ii) == 0) Is_Visible_Door(i) = .FALSE.
                   IF (Human_Known_Doors(j1)%I_nodes(ii) <= 0) Is_Known_Door(i) = .FALSE.
                END IF
             END DO Do_Known_doors0
          END IF
       END DO
    END IF

    ! If no visible nor known doors then the present target door is the only option so it
    ! is set to be known.
    IF (ABS(HR%I_Target)>0 .AND. .NOT.(ANY(Is_Visible_Door) .OR. ANY(Is_Known_Door))) &
         Is_Known_Door(ABS(HR%I_Target)) = .TRUE.


    IF (ANY(Is_Known_Door) .OR. ANY(Is_Visible_Door)) THEN
       i_tmp   = 0
       L2_min = HUGE(L2_min)
       DO i = 1, N_DOORS + N_EXITS
          IF ( Is_Known_Door(i) .AND. Is_Visible_Door(i) ) THEN
             x_o = 0.0_EB
             y_o = 0.0_EB
             N_queue = 0
             IF (TRIM(EVAC_Node_List(n_egrids+N_ENTRYS+i)%Node_Type) == 'Door' ) THEN
                x_o = 0.5_EB*(EVAC_DOORS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%X1 + &
                     EVAC_DOORS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%X2)
                y_o = 0.5_EB*(EVAC_DOORS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%Y1 + &
                     EVAC_DOORS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%Y2)
                i_o = EVAC_DOORS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%I_VENT_FFIELD
                T_tmp1 = 50.0_EB*SQRT((x1_old-x_o)**2 + (y1_old-y_o)**2)/ &
                     EVAC_DOORS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%R_NTARGET + 1.0_EB
                ii = MIN(50,MAX(1,INT(T_tmp1)-1))
                Width = EVAC_DOORS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%Width
                N_queue = EVAC_DOORS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%NTARGET(ii)
             ELSE      ! 'Exit'
                x_o = 0.5_EB*(EVAC_EXITS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%X1 + &
                     EVAC_EXITS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%X2)
                y_o = 0.5_EB*(EVAC_EXITS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%Y1 + &
                     EVAC_EXITS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%Y2)
                i_o = EVAC_EXITS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%I_VENT_FFIELD
                T_tmp1 = 50.0_EB*SQRT((x1_old-x_o)**2 + (y1_old-y_o)**2)/ &
                     EVAC_EXITS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%R_NTARGET + 1.0_EB
                ii = MIN(50,MAX(1,INT(T_tmp1)-1))
                Width = EVAC_EXITS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%Width
                N_queue = EVAC_EXITS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%NTARGET(ii)
             END IF
             IF (FED_DOOR_CRIT >= TWO_EPSILON_EB) THEN
                L2_tmp = FED_max_Door(i) * SQRT((x1_old-x_o)**2 + (y1_old-y_o)**2)/Speed
             ELSE
                L2_tmp = K_ave_Door(i)
             END IF
             IF (i_o == i_old_ffield) L2_tmp = FAC_DOOR_OLD*L2_tmp
             IF (ABS(FAC_DOOR_QUEUE) >= TWO_EPSILON_EB) THEN
                T_tmp  = SQRT((x_o-x1_old)**2 + (y_o-y1_old)**2)
                T_tmp1 = MIN(1.5_EB*Pi*T_tmp**2/(ABS(FAC_DOOR_QUEUE)*Width), REAL(N_queue,EB)/(ABS(FAC_DOOR_QUEUE)*Width))
                IF (FAC_DOOR_QUEUE < -0.001_EB) THEN
                   T_tmp = MAX((T_tmp/Speed), T_tmp1)
                ELSE
                   T_tmp = FAC_DOOR_ALPHA*(T_tmp/Speed) +  (1.0_EB-FAC_DOOR_ALPHA)*T_tmp1
                END IF
                IF (i_o == i_old_ffield) T_tmp = T_tmp*FAC_DOOR_WAIT
                IF ( T_tmp < L2_min .AND. L2_tmp < ABS(FED_DOOR_CRIT) ) THEN
                   L2_min = MAX(0.0_EB,T_tmp)
                   i_tmp = i
                END IF
             ELSE
                T_tmp  = SQRT((x_o-x1_old)**2 + (y_o-y1_old)**2)
                IF (i_o == i_old_ffield) T_tmp = T_tmp*FAC_DOOR_WAIT
                IF (T_tmp < L2_min .AND. L2_tmp < ABS(FED_DOOR_CRIT)) THEN
                   L2_min = MAX(0.0_EB,T_tmp)
                   i_tmp = i
                END IF
             END IF
             IF (i == ABS(I_Old_Target)) THEN
                IF (FED_DOOR_CRIT >= TWO_EPSILON_EB) THEN
                   IF (j > 0) THEN
                      L2_tmp = Group_List(j)%IntDose + FED_max_Door(i) * SQRT((x1_old-x_o)**2+(y1_old-y_o)**2)/Speed
                   ELSE
                      L2_tmp = HR%IntDose + FED_max_Door(i) * SQRT((x1_old-x_o)**2+(y1_old-y_o)**2)/Speed
                   END IF
                   L2_tmp2 = FED_max_Door(i) * SQRT((x1_old-x_o)**2 + (y1_old-y_o)**2)/Speed  ! no smoke criteria
                ELSE
                   ! Check that visibility > 0.5*distance to the door
                   L2_tmp = SQRT((x1_old-x_o)**2+(y1_old-y_o)**2)*0.5_EB/(3.0_EB/K_ave_Door(i))
                   L2_tmp2 = K_ave_Door(i)  ! no smoke criteria
                END IF
                L2_tmp  = FAC_DOOR_OLD2*L2_tmp
                L2_tmp2 = FAC_DOOR_OLD *L2_tmp2
                IF (HR%GROUP_ID < 0 .AND. L2_tmp2 >= ABS(FED_DOOR_CRIT)) THEN ! This is not working for groups yet nor entr agents
                   ! The present (known) target door has some smoke, record this information for later use.
                   ! If there is too much smoke => set i_nodes = 0 for this door, is not anymore known door
                   ! If there is some smoke => set i_nodes = -ABS(i_nodes), i.e., the door node number is negative
                   Do_Known_doors1: DO ii=1,Human_Known_Doors(j1)%N_nodes
                      IF (ABS(Human_Known_Doors(j1)%I_nodes(ii)) == n_egrids+N_ENTRYS+ABS(I_Old_Target)) THEN
                         IF (L2_tmp >= 1.0_EB) THEN
                            Human_Known_Doors(j1)%I_nodes(ii) = 0 ! no door, too much smoke
                         ELSE
                            Human_Known_Doors(j1)%I_nodes(ii) = -(n_egrids+N_ENTRYS+ABS(I_Old_Target)) ! some smoke
                         END IF
                         EXIT Do_Known_doors1
                      END IF
                   END DO Do_Known_doors1
                END IF
             END IF
          END IF
       END DO
       IF (i_tmp > 0 ) THEN
          ! Known and visible door, no smoke
          IF (EVAC_Node_List(n_egrids+N_ENTRYS+i_tmp)%Node_Type == 'Door' ) THEN
             name_new_ffield = TRIM(EVAC_DOORS(i_tmp)%VENT_FFIELD)
             i_new_ffield = EVAC_DOORS(i_tmp)%I_VENT_FFIELD
          ELSE        ! 'Exit'
             name_new_ffield = TRIM(EVAC_EXITS(i_tmp-N_DOORS)%VENT_FFIELD)
             i_new_ffield = EVAC_EXITS(i_tmp-N_DOORS)%I_VENT_FFIELD
          END IF
          color_index = 1
          I_Target = i_tmp
       ELSE
          ! No visible known doors available, try non-visible known doors
          i_tmp   = 0
          L2_min = HUGE(L2_min)
          DO i = 1, N_DOORS + N_EXITS
             IF ( Is_Known_Door(i) .AND. .NOT. Is_Visible_Door(i) ) THEN
                x_o = 0.0_EB
                y_o = 0.0_EB
                IF (EVAC_Node_List(n_egrids+N_ENTRYS+i)%Node_Type == 'Door' ) THEN
                   x_o = 0.5_EB*(EVAC_DOORS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%X1 + &
                        EVAC_DOORS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%X2)
                   y_o = 0.5_EB*(EVAC_DOORS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%Y1 + &
                        EVAC_DOORS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%Y2)
                   i_o = EVAC_DOORS( EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%I_VENT_FFIELD
                ELSE    ! 'Exit'
                   x_o = 0.5_EB*(EVAC_EXITS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%X1 + &
                        EVAC_EXITS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%X2)
                   y_o = 0.5_EB*(EVAC_EXITS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%Y1 + &
                        EVAC_EXITS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%Y2)
                   i_o = EVAC_EXITS( EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%I_VENT_FFIELD
                END IF
                IF (FED_DOOR_CRIT >= TWO_EPSILON_EB) THEN
                   ! L2_tmp = FED_max_Door(i)*SQRT((x1_old-x_o)**2 + (y1_old-y_o)**2)/Speed
                   ! Now L1 norm for non-visible doors (21.9.2010)
                   L2_tmp = FED_max_Door(i)*(ABS(x1_old-x_o) + ABS(y1_old-y_o))/Speed
                ELSE
                   L2_tmp = K_ave_Door(i)
                END IF
                IF (i_o == i_old_ffield) L2_tmp = FAC_DOOR_OLD*L2_tmp
                ! T_tmp  = SQRT((x_o-x1_old)**2 + (y_o-y1_old)**2)
                ! Now L1 norm for non-visible doors (21.9.2010)
                T_tmp  = ABS(x_o-x1_old) + ABS(y_o-y1_old)
                IF (i_o == i_old_ffield) T_tmp = T_tmp*FAC_DOOR_WAIT
                IF (T_tmp < L2_min .AND. L2_tmp < ABS(FED_DOOR_CRIT)) THEN
                   L2_min = MAX(0.0_EB,T_tmp)
                   i_tmp = i
                END IF
             END IF
          END DO
          IF (i_tmp > 0 ) THEN
             !  Non-visible known door, no smoke is found
             IF (EVAC_Node_List( n_egrids+N_ENTRYS+i_tmp)%Node_Type == 'Door' ) THEN
                name_new_ffield = TRIM(EVAC_DOORS(i_tmp)%VENT_FFIELD)
                i_new_ffield = EVAC_DOORS(i_tmp)%I_VENT_FFIELD
             ELSE      ! 'Exit'
                name_new_ffield = TRIM(EVAC_EXITS(i_tmp-N_DOORS)%VENT_FFIELD)
                i_new_ffield = EVAC_EXITS(i_tmp-N_DOORS)%I_VENT_FFIELD
             END IF
             color_index = 2
             I_Target = -i_tmp
          ELSE
             ! Known doors with no smoke have not been found.
             ! Try visible, not known door with no smoke.
             i_tmp   = 0
             L2_min = HUGE(L2_min)
             DO i = 1, N_DOORS + N_EXITS
                IF (Is_Visible_Door(i)) THEN
                   x_o = 0.0_EB
                   y_o = 0.0_EB
                   IF (EVAC_Node_List(n_egrids+N_ENTRYS+i)%Node_Type == 'Door' ) THEN
                      x_o = 0.5_EB*(EVAC_DOORS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%X1 + &
                           EVAC_DOORS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%X2)
                      y_o = 0.5_EB*(EVAC_DOORS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%Y1 + &
                           EVAC_DOORS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%Y2)
                      i_o = EVAC_DOORS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%I_VENT_FFIELD
                   ELSE  ! 'Exit'
                      x_o = 0.5_EB*(EVAC_EXITS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%X1 + &
                           EVAC_EXITS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%X2)
                      y_o = 0.5_EB*(EVAC_EXITS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%Y1 + &
                           EVAC_EXITS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%Y2)
                      i_o = EVAC_EXITS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%I_VENT_FFIELD
                   END IF
                   IF (FED_DOOR_CRIT >= TWO_EPSILON_EB) THEN
                      L2_tmp = FED_max_Door(i)*SQRT((x1_old-x_o)**2 + (y1_old-y_o)**2)/Speed
                   ELSE
                      L2_tmp = K_ave_Door(i)
                   END IF
                   IF (i_o == i_old_ffield) L2_tmp = FAC_DOOR_OLD*L2_tmp
                   T_tmp  = SQRT((x_o-x1_old)**2 + (y_o-y1_old)**2)
                   IF (i_o == i_old_ffield) T_tmp = T_tmp*FAC_DOOR_WAIT
                   IF (T_tmp < L2_min .AND. L2_tmp < ABS(FED_DOOR_CRIT)) THEN
                      L2_min = MAX(0.0_EB,T_tmp)
                      i_tmp = i
                   END IF
                END IF
             END DO
             IF (i_tmp > 0 ) THEN
                ! No smoke, visible door (not known) is found
                IF (EVAC_Node_List( n_egrids+N_ENTRYS+i_tmp)%Node_Type == 'Door' ) THEN
                   name_new_ffield = TRIM(EVAC_DOORS(i_tmp)%VENT_FFIELD)
                   i_new_ffield = EVAC_DOORS(i_tmp)%I_VENT_FFIELD
                ELSE    ! 'Exit'
                   name_new_ffield = TRIM( EVAC_EXITS(i_tmp-N_DOORS)%VENT_FFIELD)
                   i_new_ffield = EVAC_EXITS(i_tmp-N_DOORS)%I_VENT_FFIELD
                END IF
                color_index = 3
                I_Target = i_tmp
             ELSE
                ! Now we have some smoke and some visible or known doors
                i_tmp   = 0
                L2_min = HUGE(L2_min)
                DO i = 1, N_DOORS + N_EXITS
                   IF ( Is_Known_Door(i) .OR. Is_Visible_Door(i) ) THEN
                      x_o = 0.0_EB
                      y_o = 0.0_EB
                      IF (EVAC_Node_List(n_egrids+N_ENTRYS+i)%Node_Type == 'Door' ) THEN
                         x_o = 0.5_EB*(EVAC_DOORS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%X1 + &
                              EVAC_DOORS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%X2)
                         y_o = 0.5_EB*(EVAC_DOORS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%Y1 + &
                              EVAC_DOORS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%Y2)
                         i_o = EVAC_DOORS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%I_VENT_FFIELD
                      ELSE ! 'Exit'
                         x_o = 0.5_EB*(EVAC_EXITS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%X1 + &
                              EVAC_EXITS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%X2)
                         y_o = 0.5_EB*(EVAC_EXITS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%Y1 + &
                              EVAC_EXITS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%Y2)
                         i_o = EVAC_EXITS(EVAC_Node_List(i+n_egrids+N_ENTRYS)%Node_Index)%I_VENT_FFIELD
                      END IF
                      IF (FED_DOOR_CRIT >= TWO_EPSILON_EB) THEN
                         IF (j > 0) THEN
                            IF (Is_Visible_Door(i)) THEN
                               L2_tmp = Group_List(j)%IntDose + FED_max_Door(i) * SQRT((x1_old-x_o)**2+(y1_old-y_o)**2)/Speed
                            ELSE
                               ! Use L1 norm for non-visible doors
                               L2_tmp = Group_List(j)%IntDose + FED_max_Door(i) * (ABS(x1_old-x_o)+ABS(y1_old-y_o))/Speed
                            END IF
                         ELSE
                            IF (Is_Visible_Door(i)) THEN
                               L2_tmp = HR%IntDose + FED_max_Door(i) * SQRT((x1_old-x_o)**2+(y1_old-y_o)**2)/Speed
                            ELSE
                               ! Use L1 norm for non-visible doors
                               L2_tmp = HR%IntDose + FED_max_Door(i) * (ABS(x1_old-x_o)+ABS(y1_old-y_o))/Speed
                            END IF
                         END IF
                      ELSE
                         ! Check that visibility > 0.5*distance to the door
                         IF (Is_Visible_Door(i)) THEN
                            L2_tmp = SQRT((x1_old-x_o)**2+(y1_old-y_o)**2)*0.5_EB/(3.0_EB/K_ave_Door(i))
                         ELSE
                            L2_tmp = (ABS(x1_old-x_o)+ABS(y1_old-y_o))*0.5_EB/(3.0_EB/K_ave_Door(i))
                         END IF
                      END IF
                      IF (L2_tmp >= 1.0_EB) THEN
                         Is_Visible_Door(i) = .FALSE.  ! Too much smoke
                         Is_Known_Door(i)   = .FALSE.  ! Too much smoke
                      END IF
                      IF (i_o == i_old_ffield) L2_tmp = FAC_DOOR_OLD2*L2_tmp
                      IF (L2_tmp < L2_min) THEN
                         L2_min = L2_tmp
                         i_tmp = i
                      END IF
                   END IF
                END DO

                IF (i_tmp > 0 .AND. L2_min < 1.0_EB) THEN
                   ! Not too much smoke, (known and/or visible doors)
                   IF (EVAC_Node_List(n_egrids+N_ENTRYS+i_tmp)%Node_Type == 'Door' ) THEN
                      name_new_ffield = TRIM(EVAC_DOORS(i_tmp)%VENT_FFIELD)
                      i_new_ffield = EVAC_DOORS(i_tmp)%I_VENT_FFIELD
                   ELSE  ! 'Exit'
                      name_new_ffield = TRIM(EVAC_EXITS(i_tmp-N_DOORS)%VENT_FFIELD)
                      i_new_ffield = EVAC_EXITS(i_tmp-N_DOORS)%I_VENT_FFIELD
                   END IF
                   IF (Is_Known_Door(i_tmp) .AND. Is_Visible_Door(i_tmp)) color_index = 4
                   IF (Is_Known_Door(i_tmp) .AND. .NOT. Is_Visible_Door(i_tmp)) color_index = 5
                   IF (.NOT. Is_Known_Door(i_tmp) .AND. Is_Visible_Door(i_tmp)) color_index = 6
                   I_Target = i_tmp
                   IF (.NOT.Is_Visible_Door(i_tmp)) I_Target = -i_tmp
                ELSE    ! no match
                   ! No door found (or too much smoke), use the main evac mesh ffield by default.
                   ! Note: Now there is a herding algorithm and it is also used for "normal" agents
                   !       if they do not have any door.  If I_HERDING_TYPE > 1 then the "lost"
                   !       agents are not moving (v0=0), otherwise they are following the main
                   !       evacuation mesh flow field.  Exception: If agents came from a door/entr
                   !       to the room then they continue about 3 m away from the door and stop there.
                   i_tmp = 0
                   I_Target = i_tmp
                   IF (imode /= 2) THEN
                      i_new_ffield    = nm_tmp
                      name_new_ffield = TRIM(MESH_NAME(i_new_ffield))
                   ELSE ! Came out of a door or an entry
                      i_new_ffield = I_Field
                      name_new_ffield = TRIM(MESH_NAME(i_new_ffield))
                   END IF
                   IF (imode == 0) THEN  ! Initialization call, use evac line info
                      name_new_ffield = TRIM(HPT%GRID_NAME)
                      i_new_ffield    = HPT%I_VENT_FFIELDS(0)
                   END IF
                   color_index = EVAC_AVATAR_NCOLOR ! default, cyan
                END IF  ! case 4
             END IF    ! case 3
          END IF      ! case 2
       END IF        ! case 1
       IF (Color_Method == 4 ) THEN
          color_index = EVAC_AVATAR_NCOLOR ! default, cyan
          IF (i_tmp > 0 .AND. i_tmp <= N_DOORS ) color_index = EVAC_DOORS(i_tmp)%Avatar_Color_Index
          IF (i_tmp > N_DOORS .AND. i_tmp <= N_DOORS + N_EXITS) &
               color_index = EVAC_EXITS(i_tmp-N_DOORS)%Avatar_Color_Index
       END IF
       IF (imode == 2) THEN   ! check_target_node calls
          I_Target = i_tmp
          IF (i_tmp > 0 .AND. .NOT. Is_Visible_Door(MAX(1,i_tmp)) ) THEN
             ! I_Target >0: visible, <0: not visible, =0: No door found
             I_Target = -i_tmp
          END IF
          I_Color  = color_index
          I_Field  = i_new_ffield
          RETURN
       END IF
       I_Field = i_new_ffield

    ELSE ! No known/visible door
       i_tmp = 0 ! no door found
       color_index = EVAC_AVATAR_NCOLOR ! default, cyan
       I_Target = i_tmp
       IF (imode == 2) THEN   ! check_target_node calls
          I_Target = 0
          i_new_ffield = I_Field
          name_new_ffield = TRIM(MESH_NAME(i_new_ffield))
          I_Color  = color_index
          I_Field  = i_new_ffield
          RETURN
       END IF
       IF (HR%IEL > 0 ) THEN
          ! The agent is from some evac line
          IF (HPT%IMESH == nm_tmp) THEN
             i_new_ffield    = HPT%I_VENT_FFIELDS(0)
             name_new_ffield = TRIM(Mesh_Name(i_new_ffield))
          ELSE
             name_new_ffield = TRIM(MESH_NAME(nm_tmp))
             i_new_ffield    = nm_tmp
          END IF
          IF (Color_Method == 4) color_index = EVAC_AVATAR_NCOLOR
       ELSE
          ! The agent is from some entr line
          IF (PNX%IMESH == nm_tmp) THEN
             i_new_ffield    = PNX%I_VENT_FFIELDS(0)
             name_new_ffield = TRIM(Mesh_Name(i_new_ffield))
          ELSE
             name_new_ffield = TRIM(MESH_NAME(nm_tmp))
             i_new_ffield    = nm_tmp
          END IF
          IF (Color_Method == 4) color_index = EVAC_AVATAR_NCOLOR
       END IF
       I_Field = i_new_ffield
    END IF ! Any known or visible door
    HR%I_Target = i_tmp
    IF (ABS(i_tmp) > 0) HR%I_Door_Mode = 1
    IF (ABS(i_tmp) > 0 .AND. Is_BeeLine_Door(MAX(1,ABS(i_tmp)),1)) HR%I_Door_Mode = 2
    IF (ABS(i_tmp) > 0 .AND. Is_BeeLineXB_Door(MAX(1,ABS(i_tmp)),1)) HR%I_Door_Mode = HR%I_Door_Mode + 2
!    IF (ABS(i_tmp) > 0 .AND. .NOT.Is_BeeLine_Door(MAX(1,ABS(i_tmp))) .AND. &
!        .NOT.Is_BeeLineXB_Door(MAX(1,ABS(i_tmp)))) HR%I_Door_Mode = 1
    IF (i_tmp == 0 .AND. .NOT.HR%I_Door_Mode < 0) HR%I_Door_Mode = 0
    IF (imode < 2) I_Target = i_tmp
    IF (i_tmp > 0 .AND. .NOT. Is_Visible_Door(MAX(1,i_tmp))) THEN
       ! I_Target > 0: visible, < 0: not visible
       HR%I_Target = -i_tmp
       IF (imode < 2) I_Target = -i_tmp
    END IF
    IF (j > 0) Group_Known_Doors(j)%I_Target = HR%I_Target
    IF ( (i_new_ffield /= i_old_ffield) .OR. (imode == 0) .OR. (I_Target /= I_Old_Target)) THEN
       ! Change the field of this group/agent.
       IF ( j == 0 ) THEN
          ! Group index=0: 'lost souls' or lonely agents
          HR%I_FFIELD    = i_new_ffield
          HR%FFIELD_NAME = TRIM(name_new_ffield)
          IF (COLOR_METHOD == 5) HR%COLOR_INDEX = color_index
          IF (COLOR_METHOD == 4) HR%COLOR_INDEX = color_index
          IF (ABS(FAC_DOOR_QUEUE) >= TWO_EPSILON_EB) RETURN
          IF (imode > 0) THEN
             WRITE (LU_EVACOUT,fmt='(a,i5,a,a,a,a)') ' EVAC: Human ',ie,', new ffield: ', &
                  TRIM(name_new_ffield), ', old ffield: ',TRIM(name_old_ffield)
          ELSE
             WRITE (LU_EVACOUT,fmt='(a,i5,a,a)') ' EVAC: Human ',ie,', flow field: ', TRIM(name_new_ffield)
          END IF
       ELSE
          Group_Known_Doors(j)%I_Target = HR%I_Target
          Group_List(j)%GROUP_I_FFIELDS(i_egrid) = i_new_ffield
          HR%I_FFIELD = i_new_ffield
          HR%FFIELD_NAME = TRIM(name_new_ffield)
          IF (COLOR_METHOD == 5) HR%COLOR_INDEX = color_index
          IF (COLOR_METHOD == 4) HR%COLOR_INDEX = color_index
          Color_Tmp(j) = color_index
          IF (ABS(FAC_DOOR_QUEUE) >= TWO_EPSILON_EB) RETURN
          IF (imode > 0) THEN
             WRITE (LU_EVACOUT,fmt='(a,i5,a,a,a,a)') ' EVAC: Group ',j,', new ffield: ', &
                  TRIM(name_new_ffield), ', old ffield: ', TRIM(name_old_ffield)
          ELSE
             WRITE (LU_EVACOUT,fmt='(a,i5,a,a)') ' EVAC: Group ',ie,', flow field: ', TRIM(name_new_ffield)
          END IF
       END IF
    ELSE ! The new door is the same as the old target door.
       IF (COLOR_METHOD == 5) HR%COLOR_INDEX = color_index
       IF (COLOR_METHOD == 5 .AND. j > 0) Color_Tmp(j) = HR%COLOR_INDEX
       IF (COLOR_METHOD == 4) HR%COLOR_INDEX = color_index
       IF (COLOR_METHOD == 4 .AND. j > 0) Color_Tmp(j) = HR%COLOR_INDEX
    END IF

  END SUBROUTINE Change_Target_Door

END MODULE EVAC
