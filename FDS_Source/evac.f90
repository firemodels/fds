!
!
! This module contains the FDS+Evac human movement algorithm and
! related subprograms.  Some of the types are defined in the
! type.f90 module file.  The statistic (cumulative distributions) are in
! the ieva.f90, which includes a module dcdflib.f by netlib (http://www.netlib.org/).
! The DCDFLIB has been written in the public domain.  Some of the DCDFLIB code has
! appeared in an ACM publication and it is subject to their algorithms policy,
! see the comments at the start of the DCDFLIB in ieva.f90.
!
! Author: Timo Korhonen, VTT Technical Research Centre of Finland, 2007-2009
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
!  Use MESH_POINTERS, ONLY: DT,IJKW,BOUNDARY_TYPE,XW,YW,WALL_INDEX,POINT_TO_MESH
!  Use EVAC_MESH_POINTERS
  USE PHYSICAL_FUNCTIONS, ONLY : GET_MASS_FRACTION, FED
  USE DCDFLIB, ONLY :  DCDFLIB_Gamma => Gamma
  !
  IMPLICIT NONE
  !
  CHARACTER(255), PARAMETER :: evacid='$Id$'
  CHARACTER(255), PARAMETER :: evacrev='$Revision$'
  CHARACTER(255), PARAMETER :: evacdate='$Date$'
  !
  PRIVATE
  ! Public subprograms (called from the main program or read or dump)
  PUBLIC EVACUATE_HUMANS, INITIALIZE_EVACUATION, INIT_EVAC_GROUPS
  PUBLIC READ_EVAC, DUMP_EVAC, DUMP_EVAC_CSV, PREPARE_TO_EVACUATE
  PUBLIC EVAC_MESH_EXCHANGE, INITIALIZE_EVAC_DUMPS, GET_REV_EVAC
  ! Public variables (needed in the main program):
  ! Public variables (needed in the dump routine):
  PUBLIC N_DOORS, N_EXITS, N_ENTRYS, N_SSTANDS, EVAC_DOORS, EVAC_EXITS, EVAC_ENTRYS, EVAC_SSTANDS, & 
       EVAC_EXIT_TYPE, EVAC_DOOR_TYPE, EVAC_ENTR_TYPE, EVAC_SSTAND_TYPE
  !
  CHARACTER(255):: EVAC_VERSION = '2.1.2'
  CHARACTER(255) :: EVAC_COMPILE_DATE
  INTEGER :: EVAC_MODULE_REV
  !
  ! This is a group of persons, who are initialized together,
  ! i.e., they have same mass, speed, etc distributions and
  ! they are all put in the given rectangle.
  ! (&EVAC lines)
  TYPE EVACUATION_TYPE
     REAL(EB) :: X1=0._EB,X2=0._EB,Y1=0._EB,Y2=0._EB,Z1=0._EB,Z2=0._EB,T_START=0._EB, Angle=0._EB
     CHARACTER(60) :: CLASS_NAME='null', ID='null'
     CHARACTER(30) :: GRID_NAME='null'
     LOGICAL :: EVACFILE=.FALSE., After_Tpre=.FALSE., No_Persons=.FALSE., SHOW=.TRUE.
     INTEGER :: N_INITIAL=0,SAMPLING=0, IPC=0, IMESH=0
     INTEGER :: GN_MIN=0, GN_MAX=0
     INTEGER :: N_VENT_FFIELDS=0, Avatar_Color_Index=0
     INTEGER, DIMENSION(3) :: RGB=-1, AVATAR_RGB=-1
     INTEGER, POINTER, DIMENSION(:) :: I_DOOR_NODES =>NULL()
     INTEGER, POINTER, DIMENSION(:) :: I_VENT_FFIELDS =>NULL()
     REAL(EB), POINTER, DIMENSION(:) :: P_VENT_FFIELDS =>NULL()
  END TYPE EVACUATION_TYPE
  !
  ! An evacuatio hole, i.e., a rectangle where humans should
  ! not be put.  This makes the &EVAC lines easier to define.
  ! (&EVHO lines)
  TYPE EVAC_HOLE_TYPE
     REAL(EB) :: X1=0._EB,X2=0._EB,Y1=0._EB,Y2=0._EB,Z1=0._EB,Z2=0._EB
     CHARACTER(60) :: ID='null', PERS_ID='null', EVAC_ID='null'
     CHARACTER(30) :: GRID_NAME='null'
     INTEGER, DIMENSION(3) :: RGB=-1
     LOGICAL :: SHOW=.TRUE.
     INTEGER :: IMESH=0
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
     CHARACTER(26) :: GRID_NAME='null'
     INTEGER, DIMENSION(3) :: RGB=-1
     INTEGER :: IMESH=0, IOR=0
     REAL(EB) :: UBAR0=0._EB, VBAR0=0._EB
     LOGICAL :: Use_v0=.FALSE., SHOW=.TRUE.
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
     REAL(EB) :: A=0._EB,B=0._EB,Lambda=0._EB,C_Young=0._EB,Gamma=0._EB,Kappa=0._EB
     REAL(EB) :: r_torso=0._EB,r_shoulder=0._EB,d_shoulder=0._EB,m_iner=0._EB, Tau_iner=0._EB
     REAL(EB) :: FAC_V0_UP=-1._EB, FAC_V0_DOWN=-1._EB, FAC_V0_HORI=-1._EB
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
          TIME_OPEN=0._EB, TIME_CLOSE=0._EB, R_NTARGET=0._EB
     INTEGER :: IOR=0, ICOUNT=0, IMESH=0, INODE=0, IMODE=1
     INTEGER, DIMENSION(50) :: NTARGET=0
     REAL(EB) :: FED_CO_CO2_O2=0._EB, SOOT_DENS=0._EB, TMP_G=0._EB, RADFLUX=0._EB
     INTEGER :: II=0, JJ=0, KK=0, FED_MESH=0
     LOGICAL :: CHECK_FLOW=.FALSE., COUNT_ONLY=.FALSE., SHOW=.TRUE.
     INTEGER :: STR_INDX=0, STR_SUB_INDX=0
     CHARACTER(60) :: ID='null'
     CHARACTER(60) :: TO_NODE='null'
     CHARACTER(30) :: GRID_NAME='null'
     CHARACTER(26) :: VENT_FFIELD='null'
     INTEGER :: I_VENT_FFIELD=0, Avatar_Color_Index=0
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
          TIME_OPEN=0._EB, TIME_CLOSE=0._EB, R_NTARGET=0._EB
     INTEGER :: IOR=0, ICOUNT=0, INODE=0, INODE2=0, IMESH=0, IMESH2=0, IMODE=1
     INTEGER, DIMENSION(50) :: NTARGET=0
     INTEGER :: STR_INDX=0, STR_SUB_INDX=0
     REAL(EB) :: FED_CO_CO2_O2=0._EB, SOOT_DENS=0._EB, TMP_G=0._EB, RADFLUX=0._EB
     INTEGER :: II=0, JJ=0, KK=0, FED_MESH=0
     LOGICAL :: CHECK_FLOW=.FALSE., EXIT_SIGN=.FALSE., KEEP_XY=.FALSE., SHOW=.TRUE.
     CHARACTER(60) :: ID='null'
     CHARACTER(60) :: TO_NODE='null'
     CHARACTER(30) :: GRID_NAME='null'
     CHARACTER(26) :: VENT_FFIELD='null'
     INTEGER :: I_VENT_FFIELD=0, Avatar_Color_Index=0
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
     ! Note: Corridor may have 2 different points, where smoke etc. data
     ! is saved.
     REAL(EB), DIMENSION(2) :: FED_CO_CO2_O2=0._EB, SOOT_DENS=0._EB, TMP_G=0._EB, RADFLUX=0._EB
     INTEGER :: FED_MESH=0, FED_MESH2=0
     INTEGER, DIMENSION(2) :: II=0, JJ=0, KK=0
     INTEGER :: IOR=0, ICOUNT=0, INODE=0, INODE2=0, IMESH=0, IMESH2=0
     INTEGER :: MAX_HUMANS_INSIDE=0, n_inside=0
     LOGICAL :: CHECK_FLOW=.FALSE.
     INTEGER, DIMENSION(3) :: RGB=-1
     CHARACTER(60) :: ID='null'
     CHARACTER(60) :: TO_NODE='null'
     CHARACTER(30) :: GRID_NAME='null'
     TYPE (CORR_LL_Type), POINTER :: First =>NULL()
  END TYPE EVAC_CORR_TYPE
  ! 
  ! STRS is a construct to build a staircase. STRS consists of stairs and
  ! landings. 
  TYPE EVAC_STRS_TYPE
     REAL(EB) :: XB(6)
     REAL(EB), POINTER, DIMENSION(:,:)   :: XB_NODE =>NULL(), XB_CORE =>NULL()
     REAL(EB) :: FAC_V0_HORI=1._EB, FAC_V0_DOWN=1._EB, FAC_V0_UP=1._EB
     INTEGER :: ICOUNT=0, INODE=0, INODE2=0, IMESH=0, IMESH2=0, N_CORES = 0
     INTEGER :: N_LANDINGS, N_NODES, N_NODES_OUT, N_NODES_IN
     INTEGER, POINTER, DIMENSION(:) :: NODE_IOR =>NULL(), NODE_TYPE =>NULL(), NODES_IN =>NULL()
     INTEGER, POINTER, DIMENSION(:) :: NODES_OUT =>NULL(), I_CORE =>NULL()
     CHARACTER(60) :: ID
     CHARACTER(24) :: MESH_ID
     LOGICAL RIGHT_HANDED
  END TYPE EVAC_STRS_TYPE
  !
  ! This produces more humans on the floor specified by the
  ! coordinates. the person type ('soccer_fan' etc) are also
  ! defined here for these persons.
  ! (&ENTR lines)
  TYPE EVAC_ENTR_TYPE
     REAL(EB) :: T_first=0._EB, T_last=0._EB, Flow=0._EB, Width=0._EB, T_Start=0._EB, T_Stop=0._EB
     REAL(EB) :: X1=0._EB,X2=0._EB,Y1=0._EB,Y2=0._EB,Z1=0._EB,Z2=0._EB, Height=2.0_EB, Z=0.0_EB
     INTEGER :: IOR=0, ICOUNT=0, IPC=0, IMESH=0, INODE=0, IMODE=-1, &
          TO_INODE=0, N_Initial=0, Max_Humans=-1, &
          STR_INDX=0, STR_SUB_INDX=0
     CHARACTER(60) :: CLASS_NAME='null', ID='null'
     CHARACTER(60) :: TO_NODE='null'
     CHARACTER(30) :: GRID_NAME='null', Max_Humans_Ramp
     LOGICAL :: After_Tpre=.FALSE., No_Persons=.FALSE., SHOW=.TRUE.
     INTEGER :: N_VENT_FFIELDS=0, Avatar_Color_Index=0
     INTEGER, POINTER, DIMENSION(:) :: I_DOOR_NODES =>NULL()
     INTEGER, POINTER, DIMENSION(:) :: I_VENT_FFIELDS =>NULL()
     REAL(EB), POINTER, DIMENSION(:) :: P_VENT_FFIELDS =>NULL()
     INTEGER, DIMENSION(3) :: RGB=-1, AVATAR_RGB=-1
     REAL(EB), DIMENSION(3) :: ORIENTATION=0.0_EB
  END TYPE EVAC_ENTR_TYPE
  !
  ! coordinates. the person type ('soccer_fan' etc) are also
  ! defined here for these persons.
  TYPE EVAC_NODE_TYPE
     INTEGER :: Node_Index=0, IMESH=0
     CHARACTER(60) :: ID='null', Node_Type='null'
     CHARACTER(30) :: GRID_NAME='null'
  END TYPE EVAC_NODE_TYPE
  !
  ! Linked list, needed for the corridors
  TYPE CORR_LL_TYPE
     TYPE (HUMAN_TYPE) :: HUMAN
     REAL(EB) :: T_in=0._EB, T_out=0._EB
     LOGICAL :: From1_To2=.FALSE.
     INTEGER :: Index=0
     TYPE (CORR_LL_Type), POINTER :: Next =>NULL()
  END TYPE CORR_LL_TYPE
  !
  ! Pointers to the allocatable arrays so one can use these as
  ! shorthands to the array elements.

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
  TYPE (EVAC_NODE_Type), DIMENSION(:), ALLOCATABLE, TARGET :: Evac_Node_List

  ! Holds the information of the EVAC-lines.
  TYPE (EVACUATION_Type), DIMENSION(:), ALLOCATABLE, TARGET :: EVACUATION

  ! Holds the information of the EVHO-lines.
  TYPE (EVAC_HOLE_Type), DIMENSION(:), ALLOCATABLE, TARGET :: EVAC_HOLES

  ! Holds the information of the EVSS-lines.
  TYPE (EVAC_SSTAND_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: EVAC_SSTANDS

  ! Holds the information of the EXIT-lines.
  TYPE (EVAC_EXIT_Type), DIMENSION(:), ALLOCATABLE, TARGET :: EVAC_EXITS

  ! Holds the information of the DOOR-lines.
  TYPE (EVAC_DOOR_Type), DIMENSION(:), ALLOCATABLE, TARGET :: EVAC_DOORS

  ! Holds the information of the ENTR-lines.
  TYPE (EVAC_ENTR_Type), DIMENSION(:), ALLOCATABLE, TARGET :: EVAC_ENTRYS

  ! Holds the information of the CORR-lines.
  TYPE (EVAC_CORR_Type), DIMENSION(:), ALLOCATABLE, TARGET :: EVAC_CORRS

  ! Holds the information on the STRS-lines.
  TYPE (EVAC_STRS_Type), DIMENSION(:), ALLOCATABLE, TARGET :: EVAC_STRS

  ! Holds the information of the PERS-lines.
  TYPE (EVAC_PERS_Type), DIMENSION(:), ALLOCATABLE, TARGET :: EVAC_PERSON_CLASSES
  !
  ! Next are needed for the Gaussian random numbers
  INTEGER GaussFlag
  REAL(EB) GaussSet1, GaussSet2
  INTEGER GTrunFlag
  REAL(EB) GTrunSet1, GTrunSet2
  !
  INTEGER :: NPC_EVAC, NPC_PERS, N_EXITS, N_DOORS, N_ENTRYS, &
       N_CORRS, N_EGRIDS, N_NODES, N_HOLES, N_SSTANDS, N_STRS, N_CO_EXITS
  INTEGER :: NPPS
  INTEGER :: ILABEL_last
  CHARACTER(100) :: MESSAGE
  REAL(FB) :: EVAC_Z_MIN, EVAC_Z_MAX
  !
  REAL(EB), DIMENSION(:,:), ALLOCATABLE :: TT_Evac, FF_Evac
  INTEGER, DIMENSION(:), ALLOCATABLE :: NTT_Evac
  !
  !
  LOGICAL :: NOT_RANDOM
  INTEGER :: I_FRIC_SW, COLOR_METHOD
  REAL(EB) ::  FAC_A_WALL, FAC_B_WALL, LAMBDA_WALL, &
       NOISEME, NOISETH, NOISECM, RADIUS_COMPLETE_0, &
       RADIUS_COMPLETE_1, GROUP_EFF, FED_DOOR_CRIT, &
       TDET_SMOKE_DENS, DENS_INIT, EVAC_DT_MAX, GROUP_DENS, &
       FC_DAMPING, EVAC_DT_MIN, V_MAX, V_ANGULAR_MAX, V_ANGULAR, &
       SMOKE_MIN_SPEED, SMOKE_MIN_SPEED_VISIBILITY, TAU_CHANGE_DOOR, &
       HUMAN_SMOKE_HEIGHT, TAU_CHANGE_V0, THETA_SECTOR, CONST_DF, FAC_DF, &
       CONST_CF, FAC_CF, FAC_1_WALL, FAC_2_WALL, FAC_V0_DIR, FAC_V0_NOCF, FAC_NOCF, &
       CF_MIN_A, CF_FAC_A_WALL, CF_MIN_TAU, CF_MIN_TAU_INER, CF_FAC_TAUS, FAC_DOOR_QUEUE, &
       FAC_DOOR_WAIT, CF_MIN_B
  INTEGER, DIMENSION(3) :: DEAD_RGB
  !
  REAL(EB), DIMENSION(:), ALLOCATABLE :: Tsteps
  !
  REAL(EB), DIMENSION(:), ALLOCATABLE :: K_ave_Door
  REAL(EB), DIMENSION(:), ALLOCATABLE :: FED_max_Door
  LOGICAL, DIMENSION(:), ALLOCATABLE :: Is_Known_Door, Is_Visible_Door
  INTEGER, DIMENSION(:), ALLOCATABLE :: Color_Tmp
  !
  INTEGER :: n_dead=0, icyc_old=0, n_change_doors=0, n_change_trials=0
  REAL(EB) :: fed_max_alive, fed_max
  !
  ! Stairs constants
  INTEGER :: STRS_LANDING_TYPE=1, STRS_STAIR_TYPE=2

 ! Human constants
  INTEGER :: HUMAN_SAME_MESH_TARGET=-2,HUMAN_IMPOSSIBLE_TARGET=-1, &
             HUMAN_NO_TARGET=0,&
             HUMAN_TARGET_UNSPECIFIED = 1, &
             HUMAN_ANOTHER_MESH_TARGET=2,&
             HUMAN_CORRIDOR_TARGET=3,&
             HUMAN_STRS_TARGET=4, &
             HUMAN_EXIT_TARGET=5
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
  SUBROUTINE READ_EVAC
    IMPLICIT NONE
    !
    ! Local variables
    INTEGER :: NUMBER_INITIAL_PERSONS, COLOR_METHOD_TMP, &
         SAMPLING_FACTOR, IPC, n_tmp, GN_MIN, GN_MAX, N_RAMP_INI
    REAL(EB) :: DTSAM
    LOGICAL :: EVACFILE

    REAL(EB) :: DUMMY
    REAL(EB) :: XB(6), XB1(6), XB2(6)
    REAL(EB), DIMENSION(3) :: XYZ, XYZ_SMOKE
    INTEGER :: IOS, IZERO, N, I, J, K, IOR
    CHARACTER(30) QUANTITY, MAX_HUMANS_RAMP
    CHARACTER(60) FYI,ID,PERS_ID,TO_NODE,EVAC_ID, DEFAULT_PROPERTIES
    CHARACTER(26) FLOW_FIELD_ID
    INTEGER :: DIAMETER_DIST,VELOCITY_DIST, &
         PRE_EVAC_DIST,DET_EVAC_DIST,TAU_EVAC_DIST
    REAL(EB) :: VEL_MEAN,VEL_PARA,VEL_PARA2,VEL_LOW,VEL_HIGH, &
         DIA_MEAN,DIA_PARA,DIA_PARA2,DIA_LOW,DIA_HIGH, &
         PRE_MEAN,PRE_PARA,PRE_PARA2,PRE_LOW,PRE_HIGH, &
         DET_MEAN,DET_PARA,DET_PARA2,DET_LOW,DET_HIGH, &
         TAU_MEAN,TAU_PARA,TAU_PARA2,TAU_LOW,TAU_HIGH, &
         FCONST_A,FCONST_B,L_NON_SP, &
         C_YOUNG,GAMMA,KAPPA, ANGLE, &
         D_TORSO_MEAN,D_SHOULDER_MEAN, TAU_ROT, M_INERTIA
    INTEGER :: MAX_HUMANS_INSIDE, n_max_in_corrs, COLOR_INDEX, i_avatar_color, MAX_HUMANS
    REAL(EB) :: MAX_FLOW, WIDTH, TIME_START, TIME_STOP, WIDTH1, &
         WIDTH2, EFF_WIDTH, EFF_LENGTH, FAC_SPEED, TIME_OPEN, TIME_CLOSE
    REAL(EB) :: UBAR0, VBAR0
    LOGICAL :: CHECK_FLOW, COUNT_ONLY, AFTER_REACTION_TIME, EXIT_SIGN, KEEP_XY, USE_V0, SHOW
    LOGICAL :: OUTPUT_SPEED, OUTPUT_MOTIVE_FORCE, OUTPUT_FED, OUTPUT_OMEGA, &
         OUTPUT_ANGLE, OUTPUT_CONTACT_FORCE, OUTPUT_TOTAL_FORCE, OUTPUT_MOTIVE_ANGLE
    INTEGER, DIMENSION(3) :: RGB, AVATAR_RGB
    CHARACTER(26) :: VENT_FFIELD, MESH_ID, EVAC_MESH
    REAL(EB) :: FAC_V0_UP, FAC_V0_DOWN, FAC_V0_HORI, HEIGHT, HEIGHT0, ESC_SPEED
    CHARACTER(25) :: COLOR, DEAD_COLOR, AVATAR_COLOR

    ! Stairs variables
    REAL(EB) :: XB_CORE(4), XB_CORES(500,6), XB_LANDINGS(500,6), XB_STAIRS(500,8)
    REAL(EB) VERTICAL_LANDING_SEPARATION, STR_Length, STR_Height
    INTEGER N_LANDINGS, NL, NODES_TMP(500)
    LOGICAL :: RIGHT_HANDED, LEFT_HANDED

    CHARACTER(26), DIMENSION(51) :: KNOWN_DOOR_NAMES
    REAL(EB), DIMENSION(51) :: KNOWN_DOOR_PROBS

    INTEGER :: ii,jj,kk

    INTEGER :: size_rnd
    INTEGER, DIMENSION(8) :: t_rnd
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed_rnd

    NAMELIST /EXIT/ ID, XB, IOR, FLOW_FIELD_ID, CHECK_FLOW, &
         MAX_FLOW, FYI, COUNT_ONLY, WIDTH, XYZ, VENT_FFIELD, &
         MESH_ID, COLOR_INDEX, XYZ_SMOKE, &
         TIME_OPEN, TIME_CLOSE, EVAC_MESH, RGB, COLOR, SHOW, HEIGHT
    NAMELIST /DOOR/ ID, XB, IOR, FLOW_FIELD_ID, CHECK_FLOW, &
         MAX_FLOW, TO_NODE, FYI, WIDTH, XYZ, VENT_FFIELD, &
         EXIT_SIGN, MESH_ID, COLOR_INDEX, XYZ_SMOKE, KEEP_XY, &
         TIME_OPEN, TIME_CLOSE, EVAC_MESH, RGB, COLOR, SHOW, HEIGHT
    NAMELIST /ENTR/ ID, XB, IOR, FLOW_FIELD_ID, MAX_FLOW, &
         FYI, WIDTH, QUANTITY, PERS_ID, TIME_START, &
         TIME_STOP, AFTER_REACTION_TIME, &
         KNOWN_DOOR_NAMES, KNOWN_DOOR_PROBS, &
         MESH_ID, COLOR_INDEX, EVAC_MESH, RGB, COLOR, &
         AVATAR_COLOR, AVATAR_RGB, MAX_HUMANS, MAX_HUMANS_RAMP, SHOW, HEIGHT
    NAMELIST /CORR/ ID, XB, IOR, FLOW_FIELD_ID, CHECK_FLOW, &
         MAX_FLOW, TO_NODE, FYI, WIDTH, WIDTH1, WIDTH2, &
         EFF_WIDTH, EFF_LENGTH, MAX_HUMANS_INSIDE, FAC_SPEED, &
         XB1, XB2, RGB, COLOR
    NAMELIST /STRS/ ID, XB, XB_CORE, XB_CORES, TO_NODE, RIGHT_HANDED, LEFT_HANDED, MESH_ID, &
         N_LANDINGS, XB_LANDINGS, VERTICAL_LANDING_SEPARATION, &
         FAC_V0_UP, FAC_V0_DOWN, FAC_V0_HORI
    NAMELIST /EVAC/ NUMBER_INITIAL_PERSONS, QUANTITY, FYI, &
         ID, DTSAM, XB, FLOW_FIELD_ID, PERS_ID, &
         TIME_START, TIME_STOP, IOR, MAX_FLOW, WIDTH, ANGLE, &
         AFTER_REACTION_TIME, GN_MIN, GN_MAX, &
         KNOWN_DOOR_NAMES, KNOWN_DOOR_PROBS, MESH_ID, &
         COLOR_INDEX, EVAC_MESH, RGB, COLOR, &
         AVATAR_COLOR, AVATAR_RGB, SHOW
    NAMELIST /EVHO/ FYI, ID, XB, EVAC_ID, PERS_ID, MESH_ID, EVAC_MESH, RGB, COLOR, SHOW

    NAMELIST /EVSS/ FYI, ID, XB, MESH_ID, HEIGHT, HEIGHT0, IOR, &
         FAC_V0_UP, FAC_V0_DOWN, FAC_V0_HORI, ESC_SPEED, EVAC_MESH, RGB, COLOR, &
         UBAR0, VBAR0, USE_V0, SHOW

    NAMELIST /PERS/ FYI,ID,DIAMETER_DIST,VELOCITY_DIST, &
         PRE_EVAC_DIST,DET_EVAC_DIST,TAU_EVAC_DIST, &
         VEL_MEAN,VEL_PARA,VEL_PARA2,VEL_LOW,VEL_HIGH, &
         DIA_MEAN,DIA_PARA,DIA_PARA2,DIA_LOW,DIA_HIGH, &
         PRE_MEAN,PRE_PARA,PRE_PARA2,PRE_LOW,PRE_HIGH, &
         DET_MEAN,DET_PARA,DET_PARA2,DET_LOW,DET_HIGH, &
         TAU_MEAN,TAU_PARA,TAU_PARA2,TAU_LOW,TAU_HIGH, &
         FCONST_A,FCONST_B,L_NON_SP, &
         C_YOUNG,GAMMA,KAPPA, GROUP_DENS, &
         FAC_A_WALL, FAC_B_WALL, LAMBDA_WALL, &
         NOISEME, NOISETH, NOISECM, &
         I_FRIC_SW, GROUP_EFF, RADIUS_COMPLETE_0, &
         RADIUS_COMPLETE_1, DEFAULT_PROPERTIES, &
         NOT_RANDOM, FED_DOOR_CRIT, COLOR_METHOD, &
         TDET_SMOKE_DENS, DENS_INIT, EVAC_DT_MAX, EVAC_DT_MIN, &
         D_TORSO_MEAN, D_SHOULDER_MEAN, TAU_ROT, M_INERTIA, &
         FC_DAMPING, V_MAX, V_ANGULAR_MAX, V_ANGULAR, &
         OUTPUT_SPEED, OUTPUT_MOTIVE_FORCE, OUTPUT_FED, OUTPUT_OMEGA,&
         OUTPUT_MOTIVE_ANGLE, OUTPUT_ANGLE, OUTPUT_CONTACT_FORCE, OUTPUT_TOTAL_FORCE, &
         COLOR_INDEX, DEAD_RGB, DEAD_COLOR, &
         SMOKE_MIN_SPEED, SMOKE_MIN_SPEED_VISIBILITY, &
         TAU_CHANGE_DOOR, RGB, COLOR, &
         AVATAR_COLOR, AVATAR_RGB, HUMAN_SMOKE_HEIGHT, &
         TAU_CHANGE_V0, THETA_SECTOR, CONST_DF, FAC_DF, CONST_CF, FAC_CF, &
         FAC_1_WALL, FAC_2_WALL, FAC_V0_DIR, FAC_V0_NOCF, FAC_NOCF, &
         CF_MIN_A, CF_FAC_A_WALL, CF_MIN_TAU, CF_MIN_TAU_INER, CF_FAC_TAUS, &
         FAC_DOOR_QUEUE, FAC_DOOR_WAIT, CF_MIN_B, &
         FAC_V0_UP, FAC_V0_DOWN, FAC_V0_HORI
    !
    IF (.NOT. ANY(EVACUATION_GRID)) THEN
       N_EVAC = 0
       RETURN
    END IF

    NPPS = 30000 ! Number Persons Per Set (dump to a file)
    !
    EVAC_DT = EVAC_DT_FLOWFIELD     ! Initialize the clock
    EVAC_CLOCK = T_BEGIN    ! clock for the CHID_evac.csv file
    EVAC_N_QUANTITIES = 0

    N_RAMP_INI = N_RAMP
    i33 = 0
    ilh = 0

    IF (MYID==MAX(0,EVAC_PROCESS)) THEN
       ALLOCATE(Tsteps(NMESHES),STAT=IZERO)
       CALL ChkMemErr('READ','Tsteps',IZERO) 
       Tsteps(:) = EVAC_DT_FLOWFIELD
       IF (ABS(TIME_SHRINK_FACTOR-1.0_EB) > 0.000000000001_EB ) CALL SHUTDOWN('ERROR: Evac is not ready for TIME_SHRINK_FACTOR')
    END IF
    !
    ! I_EVAC: 'binary' index:
    ! These are just initialization. Later it is checked if files
    ! exists. If EFF/FED file does not exists, then it is calculated 
    ! (and saved).
    I_EVAC = 16*0 + 8*0 + 4*0 + 2*0 + 1*0 ! do not save soot,fed files
    IF (.NOT. ALL(EVACUATION_ONLY) ) THEN
       ! There are fire grids ==> save fed and evac flow fields
       I_EVAC = 16*1 + 8*0 + 4*0 + 2*1 + 1*1
       IF (N_REACTIONS == 0) I_EVAC = 16*1 + 8*0 + 4*0 + 2*0 + 1*1
    ELSE
       IF (EVACUATION_MC_MODE) THEN
          ! There are no fire grids ==> try to read fed and evac flow fields if MC mode
          I_EVAC = 16*0 + 8*1 + 4*1 + 2*0 + 1*0
       ELSE
          ! There are no fire grids ==> try to read fed and recalculate evac flow fields if not MC mode
          I_EVAC = 16*0 + 8*1 + 4*0 + 2*0 + 1*1
       END IF
    END IF
    !
    ! Every human has an identification number, ILABEL_last is
    ! the last used number, so next human will have an
    ! identification number, which is ILABEL_last + 1
    ILABEL_last = 0


    CALL COUNT_EVAC_NODES
    ! Write (lu_err,*) 'Evac: Counted evacuation nodes'
    CALL READ_PERS
    ! Write (lu_err,*) 'Evac: Read person classes'
    CALL READ_STRS
    ! Write (lu_err,*) 'Evac: Read stairs'
    CALL READ_EXIT
    ! Write (lu_err,*) 'Evac: Read exits'
    CALL READ_DOOR    
    ! Write (lu_err,*) 'Evac: Read doors'
    CALL READ_CORR
    ! Write (lu_err,*) 'Evac: Read corridors'
    CALL READ_ENTRIES
    ! Write (lu_err,*) 'Evac: Read entries'
    CALL COLLECT_NODE_INFO
    ! Write (lu_err,*) 'Evac: Collected node info '
    CALL READ_EVAC_LINES
    ! Write (lu_err,*) 'Evac: Read evac namelists'
    CALL READ_EVHO
    ! Write (lu_err,*) 'Evac: Read evho namelists'
    CALL READ_EVSS
    ! Write (lu_err,*) 'Evac: Read inclines'

    IF (MYID /= MAX(0,EVAC_PROCESS)) RETURN

    CALL CHECK_EVAC_NODES

  CONTAINS

    SUBROUTINE COUNT_EVAC_NODES
      IMPLICIT NONE
      !
      ! Determine total number of PERS lines in the input file
      !
      EVAC_AVATAR_NCOLOR = 0  ! Dimension of Avatar color table
      i_avatar_color     = 0  ! Counter for avatar colors
      COLOR_METHOD       = -1 ! Default is standard human colors in Smokeview
      DEAD_RGB           = (/  0,255,255/) ! cyan
      NPC_PERS = 0
      COUNT_PERS_LOOP: DO
         CALL CHECKREAD('PERS',LU_INPUT,IOS) 
         IF (IOS == 1) THEN
            EXIT COUNT_PERS_LOOP
         END IF
         READ(LU_INPUT,NML=PERS,END=221,ERR=222,IOSTAT=IOS)
         NPC_PERS = NPC_PERS + 1
222      IF (IOS > 0) CALL SHUTDOWN('ERROR: Problem with PERS line')
      END DO COUNT_PERS_LOOP
221   REWIND(LU_INPUT)
      IF (COLOR_METHOD == 3) EVAC_AVATAR_NCOLOR = NPC_PERS + 1
      COLOR_METHOD_TMP = COLOR_METHOD
      !
      ! Determine total number of EVAC lines in the input file
      !
      NPC_EVAC = 0
      COUNT_EVAC_LOOP: DO
         NUMBER_INITIAL_PERSONS = 0
         CALL CHECKREAD('EVAC',LU_INPUT,IOS) 
         IF (IOS == 1) THEN
            EXIT COUNT_EVAC_LOOP
         END IF
         READ(LU_INPUT,NML=EVAC,END=219,ERR=220,IOSTAT=IOS)
         NPC_EVAC = NPC_EVAC + 1
         IF (COLOR_METHOD == 0 .AND. NUMBER_INITIAL_PERSONS > 0) THEN
            EVAC_AVATAR_NCOLOR = EVAC_AVATAR_NCOLOR + 1
         END IF
         !
220      IF (IOS > 0) CALL SHUTDOWN('ERROR: Problem with EVAC line')
      END DO COUNT_EVAC_LOOP
219   REWIND(LU_INPUT)
      !
      ! Determine total number of EXIT lines in the input file
      !
      N_EXITS = 0
      N_CO_EXITS = 0
      COUNT_EXITS_LOOP: DO
         COUNT_ONLY = .FALSE.
         CALL CHECKREAD('EXIT',LU_INPUT,IOS) 
         IF (IOS == 1) THEN
            EXIT COUNT_EXITS_LOOP
         END IF
         READ(LU_INPUT,NML=EXIT,END=223,ERR=224,IOSTAT=IOS)
         N_EXITS = N_EXITS + 1
         IF (COLOR_METHOD == 4 .AND. .NOT.COUNT_ONLY) THEN
            EVAC_AVATAR_NCOLOR = EVAC_AVATAR_NCOLOR + 1
         END IF
         IF (COUNT_ONLY) N_CO_EXITS = N_CO_EXITS + 1
224      IF (IOS > 0) CALL SHUTDOWN('ERROR: Problem with EXIT line')
      END DO COUNT_EXITS_LOOP
223   REWIND(LU_INPUT)
      !
      ! Determine total number of DOOR lines in the input file
      !
      N_DOORS = 0
      COUNT_DOORS_LOOP: DO
         CALL CHECKREAD('DOOR',LU_INPUT,IOS) 
         IF (IOS == 1) THEN
            EXIT COUNT_DOORS_LOOP
         END IF
         READ(LU_INPUT,NML=DOOR,END=225,ERR=226,IOSTAT=IOS)
         N_DOORS = N_DOORS + 1
         IF (COLOR_METHOD == 4) THEN
            EVAC_AVATAR_NCOLOR = EVAC_AVATAR_NCOLOR + 1
         END IF
226      IF (IOS > 0) CALL SHUTDOWN('ERROR: Problem with DOOR line')
      END DO COUNT_DOORS_LOOP
225   REWIND(LU_INPUT)
      !
      ! Determine total number of ENTR lines in the input file
      !
      N_ENTRYS = 0
      COUNT_ENTRYS_LOOP: DO
         MAX_FLOW     = 0.0_EB
         CALL CHECKREAD('ENTR',LU_INPUT,IOS) 
         IF (IOS == 1) THEN
            EXIT COUNT_ENTRYS_LOOP
         END IF
         READ(LU_INPUT,NML=ENTR,END=227,ERR=228,IOSTAT=IOS)
         N_ENTRYS = N_ENTRYS + 1
         IF (COLOR_METHOD == 0 .AND. MAX_FLOW > 0) THEN
            EVAC_AVATAR_NCOLOR = EVAC_AVATAR_NCOLOR + 1
         END IF
228      IF (IOS > 0) CALL SHUTDOWN('ERROR: Problem with ENTR line')
      END DO COUNT_ENTRYS_LOOP
227   REWIND(LU_INPUT)
      !
      ! Determine total number of CORR lines in the input file
      !
      N_CORRS = 0
      COUNT_CORRS_LOOP: DO
         CALL CHECKREAD('CORR',LU_INPUT,IOS) 
         IF (IOS == 1) THEN
            EXIT COUNT_CORRS_LOOP
         END IF
         READ(LU_INPUT,NML=CORR,END=229,ERR=230,IOSTAT=IOS)
         N_CORRS = N_CORRS + 1
230      IF (IOS > 0) CALL SHUTDOWN('ERROR: Problem with CORR line')
      END DO COUNT_CORRS_LOOP
229   REWIND(LU_INPUT)
      !
      ! Determine total number of EVHO lines in the input file
      !
      N_HOLES = 0
      COUNT_EVHO_LOOP: DO
         CALL CHECKREAD('EVHO',LU_INPUT,IOS) 
         IF (IOS == 1) THEN
            EXIT COUNT_EVHO_LOOP
         END IF
         READ(LU_INPUT,NML=EVHO,END=231,ERR=232,IOSTAT=IOS)
         N_HOLES = N_HOLES + 1
232      IF (IOS > 0) CALL SHUTDOWN('ERROR: Problem with EVHO line')
      END DO COUNT_EVHO_LOOP
231   REWIND(LU_INPUT)
      !
      ! Determine total number of EVSS lines in the input file
      !
      N_SSTANDS = 0
      COUNT_EVSS_LOOP: DO
         CALL CHECKREAD('EVSS',LU_INPUT,IOS) 
         IF (IOS == 1) THEN
            EXIT COUNT_EVSS_LOOP
         END IF
         READ(LU_INPUT,NML=EVSS,END=233,ERR=234,IOSTAT=IOS)
         N_SSTANDS = N_SSTANDS + 1
234      IF (IOS > 0) CALL SHUTDOWN('ERROR: Problem with EVSS line')
      END DO COUNT_EVSS_LOOP
233   REWIND(LU_INPUT)
      !
      ! Determine total number of STRS lines in the input file
      !
      N_STRS = 0
      COUNT_STRS_LOOP: DO
         CALL CHECKREAD('STRS',LU_INPUT,IOS) 
         IF (IOS == 1) THEN
            EXIT COUNT_STRS_LOOP
         END IF
         READ(LU_INPUT,NML=STRS,END=235,ERR=236,IOSTAT=IOS)
         N_STRS = N_STRS + 1
236      IF (IOS > 0) CALL SHUTDOWN('ERROR: Problem with STRS line')
      END DO COUNT_STRS_LOOP
235   REWIND(LU_INPUT)

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
      !
      ! Allocate quantities for EVAC, PERS, EXIT types
      !
      EVAC_PROC_IF: IF (MYID==MAX(0,EVAC_PROCESS)) THEN
         IF (npc_evac > 0 ) THEN
            ALLOCATE(EVACUATION(NPC_EVAC),STAT=IZERO)
            CALL ChkMemErr('READ','EVACUATION',IZERO)
         !Else
         !   Allocate(EVACUATION(1),STAT=IZERO)
         !   Call ChkMemErr('READ','EVACUATION',IZERO)
         END IF

         IF (n_holes > 0 ) THEN
            ALLOCATE(EVAC_HOLES(N_HOLES),STAT=IZERO)
            CALL ChkMemErr('READ','EVAC_HOLES',IZERO)
         !Else 
         !   Allocate(EVAC_HOLES(1),STAT=IZERO)
         !   Call ChkMemErr('READ','EVAC_HOLES',IZERO)
         END IF

         IF (n_sstands > 0 ) THEN
            ALLOCATE(EVAC_SSTANDS(N_SSTANDS),STAT=IZERO)
            CALL ChkMemErr('READ','EVAC_SSTANDS',IZERO)
         !Else
         !   Allocate(EVAC_SSTANDS(1),STAT=IZERO)
         !   Call ChkMemErr('READ','EVAC_SSTANDS',IZERO)
         END IF

         IF (n_exits > 0 ) THEN
            ALLOCATE(EVAC_EXITS(N_EXITS),STAT=IZERO)
            CALL ChkMemErr('READ','EVAC_EXITS',IZERO) 
         !Else
         !   Allocate(EVAC_EXITS(1),STAT=IZERO)
         !   Call ChkMemErr('READ','EVAC_EXITS',IZERO) 
         END IF

         IF (n_doors > 0 ) THEN
            ALLOCATE(EVAC_DOORS(N_DOORS),STAT=IZERO)
            CALL ChkMemErr('READ','EVAC_DOORS',IZERO) 
         !Else
         !   Allocate(EVAC_DOORS(1),STAT=IZERO)
         !   Call ChkMemErr('READ','EVAC_DOORS',IZERO) 
         END IF

         IF (n_entrys > 0 ) THEN
            ALLOCATE(EVAC_ENTRYS(N_ENTRYS),STAT=IZERO)
            CALL ChkMemErr('READ','EVAC_ENTRYS',IZERO)
         !Else
         !   Allocate(EVAC_ENTRYS(1),STAT=IZERO)
         !   Call ChkMemErr('READ','EVAC_ENTRYS',IZERO)
         END IF

         IF (n_corrs > 0 ) THEN
            ALLOCATE(EVAC_CORRS(N_CORRS),STAT=IZERO)
            CALL ChkMemErr('READ','EVAC_CORRS',IZERO)
         !Else
         !   Allocate(EVAC_CORRS(1),STAT=IZERO)
         !   Call ChkMemErr('READ','EVAC_CORRS',IZERO)
         END IF

         IF (N_STRS > 0 ) THEN
            ALLOCATE(EVAC_STRS(N_STRS),STAT=IZERO)
            CALL ChkMemErr('READ','EVAC_STRS',IZERO)
         END IF

         ALLOCATE(EVAC_PERSON_CLASSES(0:NPC_PERS),STAT=IZERO)
         CALL ChkMemErr('READ','EVAC_PERSON_CLASSES',IZERO) 

         n_egrids = 0
         DO n = 1, nmeshes
            IF (evacuation_only(n) .AND. evacuation_grid(n) ) THEN
               n_egrids = n_egrids + 1
            END IF
         END  DO

         n_nodes = n_entrys + n_exits + n_doors + n_corrs + n_egrids + n_strs
         IF (n_nodes > 0 ) THEN
            ALLOCATE(EVAC_Node_List(1:n_nodes),STAT=IZERO)
            CALL ChkMemErr('READ','EVAC_NODE_LIST',IZERO) 
         END IF

         IF (npc_evac > 0 ) THEN
            !          EVACUATION(1:NPC_EVAC)%COLOR_INDEX = 1
            EVACUATION(1:NPC_EVAC)%GRID_NAME   = 'null'
            EVACUATION(1:NPC_EVAC)%CLASS_NAME  = 'null'
            EVACUATION(1:NPC_EVAC)%IMESH       = 0
            EVACUATION(1:NPC_EVAC)%ID          = 'null'
         END IF
         IF (n_holes > 0 ) THEN
            EVAC_HOLES(1:N_HOLES)%GRID_NAME   = 'null'
            EVAC_HOLES(1:N_HOLES)%PERS_ID     = 'null'
            EVAC_HOLES(1:N_HOLES)%EVAC_ID     = 'null'
            EVAC_HOLES(1:N_HOLES)%IMESH       = 0
            EVAC_HOLES(1:N_HOLES)%ID          = 'null'
         END IF

         EVAC_PERSON_CLASSES(0:NPC_PERS)%ID = 'null'

         IF (n_exits > 0 ) THEN
            EVAC_EXITS(1:N_EXITS)%ID        = 'null'
            EVAC_EXITS(1:N_EXITS)%TO_NODE   = 'null'
            EVAC_EXITS(1:N_EXITS)%GRID_NAME = 'null'
            EVAC_EXITS(1:N_EXITS)%IMESH     = 0
            !          EVAC_EXITS(1:N_EXITS)%COLOR_INDEX = 1
         END IF

         IF (n_doors > 0 ) THEN
            EVAC_DOORS(1:N_DOORS)%ID        = 'null'
            EVAC_DOORS(1:N_DOORS)%TO_NODE   = 'null'
            EVAC_DOORS(1:N_DOORS)%GRID_NAME = 'null'
            EVAC_DOORS(1:N_DOORS)%IMESH     = 0
            EVAC_DOORS(1:N_DOORS)%IMESH2    = 0
            !          EVAC_DOORS(1:N_DOORS)%COLOR_INDEX = 1
         END IF

         IF (n_corrs > 0 ) THEN
            EVAC_CORRS(1:N_CORRS)%ID        = 'null'
            EVAC_CORRS(1:N_CORRS)%TO_NODE   = 'null'
            EVAC_CORRS(1:N_CORRS)%GRID_NAME = 'null'
            EVAC_CORRS(1:N_CORRS)%IMESH     = 0
            EVAC_CORRS(1:N_CORRS)%IMESH2    = 0
         END IF

         IF (n_entrys > 0 ) THEN
            EVAC_ENTRYS(1:N_ENTRYS)%ID          = 'null'
            EVAC_ENTRYS(1:N_ENTRYS)%TO_NODE     = 'null'
            EVAC_ENTRYS(1:N_ENTRYS)%GRID_NAME   = 'null'
            EVAC_ENTRYS(1:N_ENTRYS)%CLASS_NAME  = 'null'
            EVAC_ENTRYS(1:N_ENTRYS)%IMESH       = 0
            !          EVAC_ENTRYS(1:N_ENTRYS)%COLOR_INDEX = 1
         END IF

      END IF EVAC_PROC_IF

    END SUBROUTINE COUNT_EVAC_NODES

    SUBROUTINE READ_PERS
      IMPLICIT NONE
      !
      ! Local variables
      TYPE (EVAC_PERS_Type), POINTER :: PCP=>NULL()
      !
      ! NEXT PARAMETERS ARE SAME FOR ALL HUMANS. THE LAST
      ! VALUES READ IN FROM 'PERS' LINES ARE VALID.
      FAC_A_WALL  = 1.0_EB
      FAC_B_WALL  = 0.5_EB
      LAMBDA_WALL = 0.2_EB
      NOISEME     = 0.0_EB
      NOISETH     = 0.01_EB
      NOISECM     = 3.0_EB
      NOT_RANDOM = .FALSE.
      I_FRIC_SW   = 1
      V_MAX             = 20.0_EB ! m/s
      V_ANGULAR_MAX     = 8.0_EB  ! rps
      V_ANGULAR         = 2.0_EB  ! rps
      FC_DAMPING        = 500.0_EB ! N/(m/s)
      ! Z_smoke = XB_z - EVACUATION_Z_OFFSET(NM) + HUMAN_SMOKE_HEIGHT, i.e. position
      ! of the nose/eyes above the floor.  The smoke and gas densities are
      ! taken from this level (FED calculation and visible doors etc.)
      HUMAN_SMOKE_HEIGHT   = 1.60_EB  ! Nose above floor level

      GROUP_EFF         = 0.0_EB
      RADIUS_COMPLETE_0 = 0.2_EB
      RADIUS_COMPLETE_1 = 0.5_EB

      ! TDET_SMOKE_DENS: Smoke is detected, when its density is larger than, e.g. 1 mg/m3
      !                  Default is no detection due to smoke.
      TDET_SMOKE_DENS = -999.9_EB  ! 
      FED_DOOR_CRIT = 0.000001_EB ! Which doors are 'smoke free'
      GROUP_DENS      = 0.0_EB
      SMOKE_MIN_SPEED = 0.1_EB
      SMOKE_MIN_SPEED_VISIBILITY = 0.0_EB
      TAU_CHANGE_DOOR = 1.0_EB
      DENS_INIT       = 0.0_EB
      EVAC_DT_MAX     = 0.01_EB
      EVAC_DT_MIN     = 0.001_EB

      ! Next parameters are for the counterflow (CF)
      TAU_CHANGE_V0   = -0.1_EB  !CF: How often direction is updated?
      THETA_SECTOR    = -40.0_EB !CF: The angle of the first sector
      CONST_DF        = 0.5_EB   !CF: prefer agents going in the same direction
      FAC_DF          = 1.0_EB   !CF: prefer agents going in the same direction
      CONST_CF        = 1.0_EB   !CF: dislike agents going in the opposite direction
      FAC_CF          = 2.0_EB   !CF: dislike agents going in the opposite direction
      FAC_1_WALL      = 5.0_EB   !CF: direction is towards a wall
      FAC_2_WALL      = 10.0_EB  !CF: direction leads too close to a wall
      FAC_V0_DIR      = 1.0_EB   !CF: v0*cos term for all sectors
      FAC_V0_NOCF     = 1.0_EB   !CF: prefer v0, if no counterflow
      FAC_NOCF        = 0.5_EB   !CF: prefer v0, if no counterflow
      CF_MIN_A        = 0.5_EB   !CF: decrease social force
      CF_MIN_B        = 0.3_EB   !CF: decrease social force range
      CF_FAC_A_WALL   = 1.0_EB   !CF: decrease social force
      CF_MIN_TAU      = 0.10_EB  !CF: increase motive force
      CF_MIN_TAU_INER = 0.05_EB  !CF: increase motive force, rotation
      CF_FAC_TAUS     = 0.25_EB  !CF: increase motive force, trans+rot
      FAC_DOOR_QUEUE  = 0.0_EB   ! Door selection algorithm: persons/m/s
      FAC_DOOR_WAIT   = 1.0_EB   ! Door selection algorithm: patience factor
      
      OUTPUT_SPEED         = .FALSE.
      OUTPUT_MOTIVE_FORCE  = .FALSE.
      OUTPUT_MOTIVE_ANGLE  = .FALSE.
      OUTPUT_FED           = .FALSE.
      OUTPUT_OMEGA         = .FALSE.
      OUTPUT_ANGLE         = .FALSE.
      OUTPUT_CONTACT_FORCE = .FALSE.
      OUTPUT_TOTAL_FORCE   = .FALSE.
      DEAD_COLOR = 'null'
      ! 
      ! Read the PERS lines (no read for default n=0 case)
      !
      READ_PERS_LOOP: DO N=0,NPC_PERS
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
         L_NON_SP = 0.5_EB
         C_YOUNG  = 120000.0_EB
         GAMMA    = 16000.0_EB
         KAPPA    = 40000.0_EB
         ! Rotational freedom constants
         D_TORSO_MEAN = 0.30_EB
         D_SHOULDER_MEAN = 0.19_EB
         TAU_ROT   = 0.2_EB
         M_INERTIA = -4.0_EB

         ! If not given on PERS line, use those given on EVSS lines
         FAC_V0_UP   = -1.0_EB
         FAC_V0_DOWN = -1.0_EB
         FAC_V0_HORI = -1.0_EB
         !
         ! No read for default values
         IF ( N > 0 ) THEN
            CALL CHECKREAD('PERS',LU_INPUT,IOS)
            IF (IOS == 1) THEN
               EXIT READ_PERS_LOOP
            END IF
            READ(LU_INPUT,PERS,END=24,IOSTAT=IOS)

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
            CASE Default
               WRITE(MESSAGE,'(A,A,A)') 'ERROR: PERS ',TRIM(ID),' problem with DEFAULT_PROPERTIES'
               CALL SHUTDOWN(MESSAGE)
            END SELECT

         END IF

         IF (PRE_MEAN < 0._EB .OR. PRE_LOW < 0._EB) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: PERS ',TRIM(ID), ' PRE-evacuation time should positive.'
            CALL SHUTDOWN(MESSAGE)
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
            DIA_MEAN = 0.51_EB
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

         PCP%FAC_V0_UP = FAC_V0_UP
         PCP%FAC_V0_DOWN = FAC_V0_DOWN
         PCP%FAC_V0_HORI = FAC_V0_HORI
         !
      END DO READ_PERS_LOOP
24    REWIND(LU_INPUT)
      COLOR_METHOD = COLOR_METHOD_TMP

      IF (GROUP_DENS .LE. 0.01_EB) GROUP_DENS = 0.25_EB
      IF (GROUP_DENS .GT. 3.50_EB) GROUP_DENS = 3.50_EB
      DENS_INIT = MAX(GROUP_DENS,DENS_INIT)
      IF (TDET_SMOKE_DENS < 0.0_EB) TDET_SMOKE_DENS = HUGE(TDET_SMOKE_DENS)

      IF(SMOKE_MIN_SPEED < 0.0_EB ) SMOKE_MIN_SPEED = 0.0_EB
      IF(SMOKE_MIN_SPEED_VISIBILITY < 0.01_EB) THEN
         SMOKE_MIN_SPEED_VISIBILITY = 0.01_EB  ! No divisions by zero
      END IF

      IF (.NOT. NOT_RANDOM .AND. MYID==MAX(0,EVAC_PROCESS)) THEN    ! Initialize the generator randomly
         CALL RANDOM_SEED(size=size_rnd)
         ALLOCATE(seed_rnd(size_rnd),STAT=IZERO)
         CALL ChkMemErr('READ_EVAC','seed_rnd',IZERO)
         CALL DATE_AND_TIME(values = t_rnd)
         seed_rnd = 100.0_EB*t_rnd(7) + t_rnd(8)/10.0_EB
         CALL RANDOM_SEED(put=seed_rnd)
         DEALLOCATE(seed_rnd)
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
         CALL SHUTDOWN(MESSAGE)
      END SELECT

      IF (COLOR_METHOD >= 0) EVAC_N_QUANTITIES = EVAC_N_QUANTITIES + 1
      IF (OUTPUT_MOTIVE_FORCE) EVAC_N_QUANTITIES = EVAC_N_QUANTITIES + 1
      IF (OUTPUT_MOTIVE_ANGLE) EVAC_N_QUANTITIES = EVAC_N_QUANTITIES + 1
      IF (OUTPUT_FED) EVAC_N_QUANTITIES = EVAC_N_QUANTITIES + 1
      IF (OUTPUT_SPEED) EVAC_N_QUANTITIES = EVAC_N_QUANTITIES + 1
      IF (OUTPUT_OMEGA) EVAC_N_QUANTITIES = EVAC_N_QUANTITIES + 1
      IF (OUTPUT_ANGLE) EVAC_N_QUANTITIES = EVAC_N_QUANTITIES + 1
      IF (OUTPUT_CONTACT_FORCE) EVAC_N_QUANTITIES = EVAC_N_QUANTITIES + 1
      IF (OUTPUT_TOTAL_FORCE) EVAC_N_QUANTITIES = EVAC_N_QUANTITIES + 1

      IF (EVAC_N_QUANTITIES > 0) THEN
         ALLOCATE(EVAC_QUANTITIES_INDEX(EVAC_N_QUANTITIES),STAT=IZERO)
         CALL ChkMemErr('READ','EVAC_QUANTITIES_INDEX',IZERO)

         OUTPUT_QUANTITY(240)%NAME       = 'HUMAN_MOTIVE_ACCELERATION'
         OUTPUT_QUANTITY(240)%UNITS      = 'm/s2'
         OUTPUT_QUANTITY(240)%SHORT_NAME = 'motiveAcc'

         OUTPUT_QUANTITY(241)%NAME       = 'HUMAN_FED_DOSE'
         OUTPUT_QUANTITY(241)%UNITS      = '  '
         OUTPUT_QUANTITY(241)%SHORT_NAME = 'FED'

         OUTPUT_QUANTITY(242)%NAME       = 'HUMAN_SPEED'
         OUTPUT_QUANTITY(242)%UNITS      = 'm/s'
         OUTPUT_QUANTITY(242)%SHORT_NAME = 'speed'

         OUTPUT_QUANTITY(243)%NAME       = 'HUMAN_ANGULAR_SPEED'
         OUTPUT_QUANTITY(243)%UNITS      = 'rad/s'
         OUTPUT_QUANTITY(243)%SHORT_NAME = 'omega'

         OUTPUT_QUANTITY(244)%NAME       = 'HUMAN_ACCELERATION'
         OUTPUT_QUANTITY(244)%UNITS      = 'm/s2'
         OUTPUT_QUANTITY(244)%SHORT_NAME = 'acc'

         OUTPUT_QUANTITY(245)%NAME       = 'HUMAN_CONTACT_LINEFORCE'
         OUTPUT_QUANTITY(245)%UNITS      = 'N/m'
         OUTPUT_QUANTITY(245)%SHORT_NAME = 'force_c'

         OUTPUT_QUANTITY(246)%NAME       = 'HUMAN_TOTAL_LINEFORCE'
         OUTPUT_QUANTITY(246)%UNITS      = 'N/m'
         OUTPUT_QUANTITY(246)%SHORT_NAME = 'force_tot'

         OUTPUT_QUANTITY(247)%NAME       = 'HUMAN_COLOR'
         OUTPUT_QUANTITY(247)%UNITS      = '  '
         OUTPUT_QUANTITY(247)%SHORT_NAME = 'color'

         OUTPUT_QUANTITY(248)%NAME       = 'HUMAN_MOTIVE_ANGLE'
         OUTPUT_QUANTITY(248)%UNITS      = 'rad'
         OUTPUT_QUANTITY(248)%SHORT_NAME = 'mAngle'

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
         IF (OUTPUT_ANGLE) THEN
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

         IF ( n-1 .NE. EVAC_N_QUANTITIES ) THEN
            WRITE(MESSAGE,'(A,2I4,A)') 'ERROR: Evac output quantities ',EVAC_N_QUANTITIES,n-1, ' Some bug in the program.'
            CALL SHUTDOWN(MESSAGE)
         END IF
      END IF

    END SUBROUTINE READ_PERS

    SUBROUTINE READ_EXIT
      IMPLICIT NONE
      !
      ! Read the EXIT lines
      !
      ! Local variables
      INTEGER nm, i1, i2, j1, j2, i, ii, iii
      TYPE (EVAC_EXIT_Type), POINTER :: PEX=>NULL()
      TYPE (EVAC_STRS_Type),  POINTER :: STRP=>NULL()
      TYPE (MESH_TYPE), POINTER :: M=>NULL()
      !
      READ_EXIT_LOOP: DO N = 1, N_EXITS
         !
         ID            = 'null'
         RGB   = -1
         COLOR = 'null'
         XB            = 0.0_EB
         IOR           = 0
         FLOW_FIELD_ID = 'null'
         VENT_FFIELD   = 'null'
         MESH_ID       = 'null'
         EVAC_MESH     = 'null'
         CHECK_FLOW    = .FALSE.
         COUNT_ONLY    = .FALSE.
         SHOW          = .TRUE.
         MAX_FLOW      = 0.0_EB
         WIDTH         = 0.0_EB
         HEIGHT        = 2.0_EB
         TIME_OPEN     = -HUGE(TIME_OPEN)
         TIME_CLOSE    = HUGE(TIME_CLOSE)
         XYZ(:)        = HUGE(XYZ)
         XYZ_SMOKE(:)  = HUGE(XYZ_SMOKE)
         COLOR_INDEX   = -1
         !
         CALL CHECKREAD('EXIT',LU_INPUT,IOS)
         IF (IOS == 1) THEN
            EXIT READ_EXIT_LOOP
         END IF
         READ(LU_INPUT,EXIT,END=26,IOSTAT=IOS)
         !
         ! Old input used COLOR_INDEX, next lines are needed for that
         IF (MYID==MAX(0,EVAC_PROCESS) .AND. COLOR_INDEX.NE.-1) WRITE (LU_ERR,'(A,A)') &
              ' WARNING: keyword COLOR_INDEX is replaced by COLOR at EXIT line ',TRIM(ID)
         IF (COLOR_INDEX == 1) COLOR = 'BLACK'  
         IF (COLOR_INDEX == 2) COLOR = 'YELLOW' 
         IF (COLOR_INDEX == 3) COLOR = 'BLUE'   
         IF (COLOR_INDEX == 4) COLOR = 'RED'    
         IF (COLOR_INDEX == 5) COLOR = 'GREEN'  
         IF (COLOR_INDEX == 6) COLOR = 'MAGENTA'
         IF (COLOR_INDEX == 7) COLOR = 'CYAN'   

         ! Colors, integer RGB(3), e.g., (23,255,0)
         IF (ANY(RGB < 0) .AND. COLOR=='null') COLOR = 'FOREST GREEN'
         IF (COLOR /= 'null') CALL COLOR2RGB(RGB,COLOR)
         IF (COLOR_METHOD == 4 .AND. .NOT.COUNT_ONLY) THEN
            i_avatar_color = i_avatar_color + 1
            EVAC_AVATAR_RGB(1:3,i_avatar_color) = RGB
         END IF

         IF (MYID /= MAX(0,EVAC_PROCESS)) CYCLE READ_EXIT_LOOP

         PEX=>EVAC_EXITS(N)

         PEX%RGB = RGB
         IF (COLOR_METHOD == 4 .AND. .NOT.COUNT_ONLY) PEX%Avatar_Color_Index = i_avatar_color

         IF (EVAC_MESH /= 'null') THEN
            MESH_ID = EVAC_MESH
            IF (MYID==MAX(0,EVAC_PROCESS)) WRITE (LU_ERR,'(A,A)') &
                 ' WARNING: keyword EVAC_MESH is replaced by MESH_ID at EXIT line ', TRIM(ID)
         END IF

         ! Check that the exit is properly specified
 
         DO I=1,5,2
            IF (XB(I) > XB(I+1)) THEN
               DUMMY   = XB(I)
               XB(I)   = XB(I+1)
               XB(I+1) = DUMMY
            END IF
         END DO
         ! 
         ! Check which evacuation mesh
         ii = 0
         PEX_MeshLoop: DO i = 1, nmeshes
            IF (evacuation_only(i) .AND. evacuation_grid(i)) THEN
               IF (Is_Within_Bounds(XB(1),XB(2),XB(3),XB(4),XB(5),XB(6),&
                  Meshes(i)%XS,Meshes(i)%XF,Meshes(i)%YS,Meshes(i)%YF,Meshes(i)%ZS,Meshes(i)%ZF, 0._EB, 0._EB, 0._EB)) THEN
                  IF (TRIM(MESH_ID) == 'null' .OR. TRIM(MESH_ID) == TRIM(MESH_NAME(i))) THEN
                     ii = ii + 1
                     PEX%IMESH = i
                  END IF
               END IF
            END IF
         END DO PEX_MeshLoop
         IF (PEX%IMESH == 0) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: EXIT ',TRIM(ID), ' problem with IMESH, no mesh found'
            CALL SHUTDOWN(MESSAGE)
         END IF
         IF (ii > 1) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: EXIT ',TRIM(ID), ' not an unique mesh found '
            CALL SHUTDOWN(MESSAGE)
         END IF
 
         nm = PEX%IMESH
         M => MESHES(PEX%IMESH)
 
         IF (XB(1)/=XB(2) .AND. XB(3)/=XB(4)) THEN
            WRITE(MESSAGE,'(A,I4,A)') 'ERROR: EXIT ',TRIM(ID),' must be a plane'
            CALL SHUTDOWN(MESSAGE)
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
            WRITE(LU_ERR,fmt='(a,6f12.4)') 'Old XB:', PEX%X1,PEX%X2,PEX%Y1,PEX%Y2,PEX%Z1,PEX%Z2
            WRITE(LU_ERR,fmt='(a,6f12.4)') 'New XB:', XB(1:6)
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
         PEX%IMODE      = -1 ! Exit is open by default
         IF (TIME_OPEN>TIME_CLOSE) THEN
            PEX%IMODE=-1
         ELSE
            IF (TIME_OPEN>T_BEGIN) PEX%IMODE=+2
         END IF
         IF (CHECK_FLOW) PEX%Flow_max   = MAX_FLOW
         PEX%COUNT_ONLY = .FALSE.
         IF (COUNT_ONLY) PEX%COUNT_ONLY = .TRUE.
         PEX%SHOW = SHOW
         IF (COUNT_ONLY  ) PEX%SHOW = .FALSE.
         PEX%HEIGHT = HEIGHT

         !       PEX%COLOR_INDEX = Mod(Max(0,COLOR_INDEX-1),7) + 1 ! 1-7 always

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
         CASE (-3)
            IF ( (XB(4)-XB(3)) <= 0.0_EB .OR. (XB(2)-XB(1)) <= 0.0_EB) THEN
               WRITE(MESSAGE,'(A,I4,A)') 'ERROR: EXIT',N,' IOR=-3 but not 3-dim object'
               CALL SHUTDOWN(MESSAGE)
            END IF
            PEX%ORIENTATION(3)=REAL(SIGN(1,IOR),EB)
         CASE (0)
            IF ( (XB(4)-XB(3)) <= 0.0_EB .OR. (XB(2)-XB(1)) <= 0.0_EB) THEN
               WRITE(MESSAGE,'(A,I4,A)') 'ERROR: EXIT',N,' no IOR but not 3-dim object'
               CALL SHUTDOWN(MESSAGE)
            END IF
         CASE Default
            WRITE(MESSAGE,'(A,I4,A)') 'ERROR: EXIT',N,' problem with IOR'
            CALL SHUTDOWN(MESSAGE)
         END SELECT
         ! 
         ! Check which vent field. If VENT_FFIELD is not found, use
         ! the main evac grid.
         PEX%I_VENT_FFIELD = 0
         PEX_Mesh2Loop: DO i = 1, nmeshes
            IF ( evacuation_only(i) .AND. (TRIM(MESH_NAME(i)) == TRIM(PEX%VENT_FFIELD)) ) THEN
               IF ( (PEX%Z1 >= Meshes(i)%ZS .AND. PEX%Z2 <= Meshes(i)%ZF).AND. &
                    (PEX%Y1 >= Meshes(i)%YS .AND. PEX%Y2 <= Meshes(i)%YF).AND. &
                    (PEX%X1 >= Meshes(i)%XS .AND. PEX%X2 <= Meshes(i)%XF)) THEN
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
         ! Check which evacuation floor
         ! Now there may be overlapping meshes.
         ii  = 0
         iii = 0
         PEX_Mesh3Loop: DO i = 1, nmeshes
            IF (evacuation_only(i) .AND. evacuation_grid(i)) THEN
               IF ( (PEX%Z >= Meshes(i)%ZS .AND. PEX%Z <= Meshes(i)%ZF).AND. &
                    (PEX%Y >= Meshes(i)%YS .AND. PEX%Y <= Meshes(i)%YF).AND. &
                    (PEX%X >= Meshes(i)%XS .AND. PEX%X <= Meshes(i)%XF)) THEN
                  IF (TRIM(MESH_ID) == 'null' .OR. TRIM(MESH_ID) == TRIM(MESH_NAME(i))) THEN
                     ii  = ii + 1
                     iii = i
                  END IF
               END IF
            END IF
         END DO PEX_Mesh3Loop
         IF (ii == 0 .AND. .NOT.COUNT_ONLY .AND. SHOW) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: EXIT line ',TRIM(PEX%ID), ' problem with XYZ, no mesh found'
            CALL SHUTDOWN(MESSAGE)
         END IF
         IF (ii > 1) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: EXIT line ',TRIM(PEX%ID), ' problem with XYZ, not an unique mesh found'
            CALL SHUTDOWN(MESSAGE)
         END IF
         IF (iii /= PEX%IMESH .AND. .NOT.COUNT_ONLY .AND. SHOW) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: EXIT line ',TRIM(PEX%ID), ' problem with XYZ, wrong mesh found'
            CALL SHUTDOWN(MESSAGE)
         END IF
         M => MESHES(PEX%IMESH)
         II = FLOOR(M%CELLSI(FLOOR((PEX%X-M%XS)*M%RDXINT))+ 1.0_EB)
         JJ = FLOOR(M%CELLSJ(FLOOR((PEX%Y-M%YS)*M%RDYINT))+ 1.0_EB)
         KK = FLOOR(M%CELLSK(FLOOR((PEX%Z-M%ZS)*M%RDZINT))+ 1.0_EB)
         IF (M%SOLID(M%CELL_INDEX(II,JJ,KK)) .AND. .NOT.COUNT_ONLY .AND. SHOW) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: EXIT line ',TRIM(PEX%ID), ' problem with XYZ, inside solid'
            CALL SHUTDOWN(MESSAGE)
  
         END IF

         ! PEX%Z is used to plot the door on the correct height in Smokeview.
         PEX%Z = PEX%Z + 0.5_EB*PEX%Height - EVACUATION_Z_OFFSET(PEX%IMESH)

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
         PEX_SmokeLoop: DO i = 1, nmeshes
            IF (.NOT. evacuation_only(i)) THEN
               IF ( Is_Within_Bounds(PEX%Xsmoke,PEX%Xsmoke,PEX%Ysmoke,PEX%Ysmoke,PEX%Zsmoke,PEX%Zsmoke,&
                    Meshes(i)%XS,Meshes(i)%XF,Meshes(i)%YS,Meshes(i)%YF,Meshes(i)%ZS,Meshes(i)%ZF,0._EB,0._EB,0._EB)) THEN
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
      TYPE (EVAC_DOOR_Type), POINTER :: PDX=>NULL()
      TYPE (MESH_TYPE), POINTER :: M=>NULL()
      INTEGER nm, i1, i2, j1, j2, i, ii, iii
      !
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
         EXIT_SIGN     = .FALSE.
         SHOW          = .TRUE.
         MAX_FLOW      = 0.0_EB
         WIDTH         = 0.0_EB
         HEIGHT        = 2.0_EB
         TIME_OPEN     = -HUGE(TIME_OPEN)
         TIME_CLOSE    = HUGE(TIME_CLOSE)
         XYZ(:)        = HUGE(XYZ)
         XYZ_SMOKE(:)  = HUGE(XYZ_SMOKE)
         COLOR_INDEX   = -1
         KEEP_XY       = .FALSE.
         !
         CALL CHECKREAD('DOOR',LU_INPUT,IOS)
         IF (IOS == 1) THEN
            EXIT READ_DOOR_LOOP
         END IF
         READ(LU_INPUT,DOOR,END=27,IOSTAT=IOS)
         !
         ! Old input used COLOR_INDEX, next lines are needed for that
         IF (MYID==MAX(0,EVAC_PROCESS) .AND. COLOR_INDEX.NE.-1) WRITE (LU_ERR,'(A,A)') &
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

         PDX%RGB = RGB
         IF (COLOR_METHOD == 4) PDX%Avatar_Color_Index = i_avatar_color

         IF (EVAC_MESH /= 'null') THEN
            MESH_ID = EVAC_MESH
            IF (MYID==MAX(0,EVAC_PROCESS)) WRITE (LU_ERR,'(A,A)') &
                 ' WARNING: keyword EVAC_MESH is replaced by MESH_ID at DOOR line ', TRIM(ID)
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
         ! Check which evacuation floor
         ! Now there may be overlapping meshes.
         ii = 0
         PDX_MeshLoop: DO i = 1, nmeshes
            IF (evacuation_only(i) .AND. evacuation_grid(i)) THEN
               IF (Is_Within_Bounds(XB(1),XB(2),XB(3),XB(4),XB(5),XB(6),&
                  Meshes(i)%XS,Meshes(i)%XF,Meshes(i)%YS,Meshes(i)%YF,Meshes(i)%ZS,Meshes(i)%ZF, 0._EB, 0._EB, 0._EB)) THEN
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
         IF (PDX%IMESH == 0) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: DOOR ',TRIM(ID), ' problem with IMESH, no mesh found'
            CALL SHUTDOWN(MESSAGE)
         END IF
         IF (ii > 1) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: DOOR ',TRIM(ID), ' not an unique mesh found '
            CALL SHUTDOWN(MESSAGE)
         END IF

         nm = PDX%IMESH
         M => MESHES(PDX%IMESH)
 
         IF (XB(1)/=XB(2) .AND. XB(3)/=XB(4)) THEN
            WRITE(MESSAGE,'(A,I4,A)') 'ERROR: DOOR',TRIM(ID),' must be a plane'
            CALL SHUTDOWN(MESSAGE)
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
            WRITE(LU_ERR,fmt='(a,6f12.4)') 'Old XB:', PDX%X1,PDX%X2,PDX%Y1,PDX%Y2,PDX%Z1,PDX%Z2
            WRITE(LU_ERR,fmt='(a,6f12.4)') 'New XB:', XB(1:6)
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
         PDX%IOR        = IOR
         PDX%ID         = ID
         PDX%GRID_NAME  = FLOW_FIELD_ID
         PDX%VENT_FFIELD= VENT_FFIELD
         PDX%CHECK_FLOW = CHECK_FLOW
         PDX%EXIT_SIGN  = EXIT_SIGN
         PDX%KEEP_XY    = KEEP_XY
         PDX%SHOW       = SHOW
         PDX%TO_NODE    = TO_NODE
         PDX%INODE      = 0
         PDX%INODE2     = 0
         PDX%T_first    = T_BEGIN
         PDX%T_last     = T_BEGIN
         PDX%ICOUNT     = 0
         PDX%Flow_max   = 0.0_EB
         PDX%TIME_OPEN  = TIME_OPEN
         PDX%TIME_CLOSE = TIME_CLOSE
         PDX%IMODE      = -1 ! Door is open by default
         IF (TIME_OPEN>TIME_CLOSE) THEN
            PDX%IMODE=-1
         ELSE
            IF (TIME_OPEN>T_BEGIN) PDX%IMODE=+2
         END IF
         IF (CHECK_FLOW) PDX%Flow_max   = MAX_FLOW
         PDX%HEIGHT = HEIGHT

         !       PDX%COLOR_INDEX = Mod(Max(0,COLOR_INDEX-1),7) ! 1-7 always

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
               WRITE(MESSAGE,'(A,I4,A)') 'ERROR: DOOR',N,' IOR=-3 but not 3-dim object'
               CALL SHUTDOWN(MESSAGE)
            END IF
            PDX%ORIENTATION(3)=REAL(SIGN(1,IOR),EB)
         CASE (0)
            IF ( (XB(4)-XB(3)) <= 0.0_EB .OR. (XB(2)-XB(1)) <= 0.0_EB) THEN
               WRITE(MESSAGE,'(A,I4,A)') 'ERROR: DOOR',N,' no IOR but not 3-dim object'
               CALL SHUTDOWN(MESSAGE)
            END IF
         CASE Default
            WRITE(MESSAGE,'(A,I4,A)') 'ERROR: DOOR',N,' problem with IOR'
            CALL SHUTDOWN(MESSAGE)
         END SELECT
         ! 
         ! Check which vent field. If VENT_FFIELD is not found, use
         ! the main evac grid.
         PDX%I_VENT_FFIELD = 0
         PDX_Mesh2Loop: DO i = 1, nmeshes
            IF ( evacuation_only(i) .AND. (TRIM(MESH_NAME(i)) == TRIM(PDX%VENT_FFIELD)) ) THEN
               IF ( (PDX%Z1 >= Meshes(i)%ZS .AND. PDX%Z2 <= Meshes(i)%ZF).AND. &
                    (PDX%Y1 >= Meshes(i)%YS .AND. PDX%Y2 <= Meshes(i)%YF).AND. &
                    (PDX%X1 >= Meshes(i)%XS .AND. PDX%X2 <= Meshes(i)%XF)) THEN
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

         IF (XYZ_SMOKE(1) < HUGE(XYZ_SMOKE)) THEN
            PDX%Xsmoke = XYZ_SMOKE(1)
            PDX%Ysmoke = XYZ_SMOKE(2)
            PDX%Zsmoke = XYZ_SMOKE(3)
         ELSE
            PDX%Xsmoke = PDX%X
            PDX%Ysmoke = PDX%Y
            PDX%Zsmoke = 0.5_EB*(XB(5)+XB(6)) - EVACUATION_Z_OFFSET(PDX%IMESH) + HUMAN_SMOKE_HEIGHT
         END IF

         ! Check which evacuation floor
         ii = 0
         iii = 0
         PDX_Mesh3Loop: DO i = 1, nmeshes
            IF (evacuation_only(i) .AND. evacuation_grid(i)) THEN
               IF ( (PDX%Z >= Meshes(i)%ZS .AND. PDX%Z <= Meshes(i)%ZF).AND. &
                    (PDX%Y >= Meshes(i)%YS .AND. PDX%Y <= Meshes(i)%YF).AND. &
                    (PDX%X >= Meshes(i)%XS .AND. PDX%X <= Meshes(i)%XF)) THEN
                  IF (TRIM(MESH_ID) == 'null' .OR. TRIM(MESH_ID) == TRIM(MESH_NAME(i))) THEN
                     ii  = ii + 1
                     iii = i
                  END IF
               END IF
            END IF
         END DO PDX_Mesh3Loop
         IF (ii == 0) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: DOOR line ',TRIM(PDX%ID), ' problem with XYZ, no mesh found'
            CALL SHUTDOWN(MESSAGE)
         END IF
         IF (ii > PDX%IMESH) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: DOOR line ',TRIM(PDX%ID), ' problem with XYZ, not an unique mesh found'
            CALL SHUTDOWN(MESSAGE)
         END IF
         IF (iii /= PDX%IMESH .AND. SHOW) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: DOOR line ',TRIM(PDX%ID), ' problem with XYZ, wrong mesh found'
            CALL SHUTDOWN(MESSAGE)
         END IF
         M => MESHES(PDX%IMESH)
         II = FLOOR(M%CELLSI(FLOOR((PDX%X-M%XS)*M%RDXINT))+ 1.0_EB)
         JJ = FLOOR(M%CELLSJ(FLOOR((PDX%Y-M%YS)*M%RDYINT))+ 1.0_EB)
         KK = FLOOR(M%CELLSK(FLOOR((PDX%Z-M%ZS)*M%RDZINT))+ 1.0_EB)
         IF (M%SOLID(M%CELL_INDEX(II,JJ,KK))) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: DOOR line ',TRIM(PDX%ID), ' problem with XYZ, inside solid'
            CALL SHUTDOWN(MESSAGE)
         END IF
         ! PDX%Z is used to plot the door on the correct height in Smokeview.
         PDX%Z = PDX%Z + 0.5_EB*PDX%Height - EVACUATION_Z_OFFSET(PDX%IMESH)

         ! Check, which fire grid and i,j,k (xyz)
         PDX_SmokeLoop: DO i = 1, nmeshes
            IF (.NOT. evacuation_only(i)) THEN
               IF ( (PDX%Zsmoke >= Meshes(i)%ZS .AND. PDX%Zsmoke <= Meshes(i)%ZF) .AND. &
                    (PDX%Ysmoke >= Meshes(i)%YS .AND. PDX%Ysmoke <= Meshes(i)%YF) .AND. &
                    (PDX%Xsmoke >= Meshes(i)%XS .AND. PDX%Xsmoke <= Meshes(i)%XF)) THEN
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
      TYPE (EVAC_CORR_Type), POINTER :: PCX=>NULL()
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
         MAX_FLOW      = 0.0_EB
         WIDTH         = 0.0_EB
         WIDTH1        = 0.0_EB
         WIDTH2        = 0.0_EB
         FAC_SPEED     = 0.0_EB
         EFF_WIDTH     = 0.0_EB
         EFF_LENGTH    = 0.0_EB
         MAX_HUMANS_INSIDE = 0
         !
         CALL CHECKREAD('CORR',LU_INPUT,IOS)
         IF (IOS == 1) THEN
            EXIT READ_CORR_LOOP
         END IF
         READ(LU_INPUT,CORR,END=29,IOSTAT=IOS)
         !
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
         ! Position, where smoke etc. is saved.
         ! If both XB and XB1 are given, use XB1
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

         PCX%MAX_HUMANS_INSIDE = 0
         IF (MAX_HUMANS_INSIDE > 0 ) THEN
            PCX%MAX_HUMANS_INSIDE = MAX_HUMANS_INSIDE
         ELSE
            WRITE(MESSAGE,'(A,I4,A)') 'ERROR: CORR',N,' MAX_HUMANS_INSIDE <= 0'
            CALL SHUTDOWN(MESSAGE)
         END IF

         IF (FAC_SPEED < 0 ) THEN
            WRITE(MESSAGE,'(A,I4,A)') 'ERROR: CORR',N,' FAC_SPEED < 0'
            CALL SHUTDOWN(MESSAGE)
         ELSE
            IF (FAC_SPEED == 0.0_EB) FAC_SPEED = 0.6_EB
            PCX%Fac_Speed = FAC_SPEED
         END IF

         PCX%Flow_max   = 0.0_EB
         IF (CHECK_FLOW) PCX%Flow_max   = MAX_FLOW

         PCX%Width = MAX( ABS(XB(4)-XB(3)) , ABS(XB(2)-XB(1)) )

         PCX%Eff_Width = 0.0_EB
         IF (EFF_WIDTH > 0.0_EB ) THEN
            PCX%Eff_Width = EFF_WIDTH
         ELSE
            PCX%Eff_Width = PCX%Width
         END IF

         PCX%Eff_Length = 0.0_EB
         IF (EFF_LENGTH > 0.0_EB ) THEN
            PCX%Eff_Length = EFF_LENGTH
         ELSE
            WRITE(MESSAGE,'(A,I4,A)') 'ERROR: CORR',N,' EFF_LENGTH <= 0'
            CALL SHUTDOWN(MESSAGE)
         END IF
         PCX%Eff_Area = PCX%Eff_Length*PCX%Eff_Width

         PCX%Width1 = WIDTH1
         PCX%Width2 = WIDTH2
         IF (WIDTH1*WIDTH2 <= 0.0_EB) THEN
            PCX%Width1 = PCX%Width
            PCX%Width2 = PCX%Width
         END IF
         ! 
         PCX_MeshLoop: DO i = 1, nmeshes
            IF (.NOT. evacuation_only(i) .AND. PCX%FED_MESH >= 0) THEN
               IF ( (PCX%Z1 >= Meshes(i)%ZS .AND. PCX%Z1 <= Meshes(i)%ZF).AND. &
                    (PCX%Y1 >= Meshes(i)%YS .AND. PCX%Y1 <= Meshes(i)%YF).AND. &
                    (PCX%X1 >= Meshes(i)%XS .AND. PCX%X1 <= Meshes(i)%XF)) THEN
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

         PCX_MeshLoop2: DO i = 1, nmeshes
            IF (.NOT. evacuation_only(i) .AND. PCX%FED_MESH2 >= 0) THEN
               IF ( (PCX%Z2 >= Meshes(i)%ZS .AND. PCX%Z2 <= Meshes(i)%ZF).AND. &
                    (PCX%Y2 >= Meshes(i)%YS .AND. PCX%Y2 <= Meshes(i)%YF).AND. &
                    (PCX%X2 >= Meshes(i)%XS .AND. PCX%X2 <= Meshes(i)%XF)) THEN
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
      TYPE (EVAC_STRS_Type), POINTER :: STRP=>NULL()

      ! Read the STRS line
      READ_STRS_LOOP: DO N = 1,N_STRS
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
         CALL CHECKREAD('STRS',LU_INPUT,IOS)
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
         STRP_MeshLoop: DO I = 1, NMESHES
            IF (.NOT. evacuation_only(I)) CYCLE
            IF (.NOT. evacuation_grid(I)) CYCLE
            IF (TRIM(MESH_ID) == 'null' .OR. TRIM(MESH_ID)==TRIM(MESH_NAME(I))) THEN
               ii = ii + 1
               STRP%IMESH = I
               EXIT STRP_MeshLoop
            END IF
         END DO STRP_MeshLoop

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
            WRITE(MESSAGE,'(A,I4,A)') 'ERROR: STRS',N,' N_LANDINGS > 500'
            CALL SHUTDOWN(MESSAGE)
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
                CALL SHUTDOWN(MESSAGE)
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
               IF ((Z_TMP >= STRP%XB_CORE(J,5)) .AND. (Z_TMP <= STRP%XB_CORE(J,6))) EXIT
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
               IF (STR_Length > 0._EB) XB_STAIRS(I,7) = COS(ATAN(STR_Height/STR_Length))
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
               IF (STR_Length > 0._EB) XB_STAIRS(I,7) = COS(ATAN(STR_Height/STR_Length))
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
               IF (STR_Length > 0._EB) XB_STAIRS(I,8) = COS(ATAN(STR_Height/STR_Length))
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
               IF (STR_Length > 0._EB) XB_STAIRS(I,8) = COS(ATAN(STR_Height/STR_Length))
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
               IF ((Z_TMP >= STRP%XB_CORE(K,5)) .AND. (Z_TMP <= STRP%XB_CORE(K,6))) EXIT
            ENDDO
            K = MIN(K,STRP%N_CORES)
            STRP%I_CORE(J) = K
            IF (I<N_LANDINGS) THEN
               STRP%XB_NODE(J+1,1:8) = XB_STAIRS(I,:)
               STRP%NODE_TYPE(J+1)   = STRS_STAIR_TYPE
               ! Choose core
               Z_TMP = 0.5_EB * (STRP%XB_NODE(J+1,5)+STRP%XB_NODE(J+1,6))
               DO K = 1,STRP%N_CORES
                  IF ((Z_TMP >= STRP%XB_CORE(K,5)) .AND. (Z_TMP <= STRP%XB_CORE(K,6))) EXIT
               ENDDO
               K = MIN(K,STRP%N_CORES)
               STRP%I_CORE(J+1) = K
            END IF
            J = J + 2
         END DO

         ! Compute velocity fields
!         Do NM = 1, NMESHES
!            If (MESH_ID == MESH_NAME(NM)) Then
!               M=>MESHES(NM)
!               Allocate(STRP%U_UP(1:M%IBAR,1:M%JBAR),STAT=IZERO)
!               Call ChkMemErr('Read_Evac','STRP%U_UP',IZERO) 
!               Allocate(STRP%V_UP(1:M%IBAR,1:M%JBAR),STAT=IZERO)
!               Call ChkMemErr('Read_Evac','STRP%V_UP',IZERO) 
!               Allocate(STRP%U_DOWN(1:M%IBAR,1:M%JBAR),STAT=IZERO)
!               Call ChkMemErr('Read_Evac','STRP%U_DOWN',IZERO) 
!               Allocate(STRP%V_DOWN(1:M%IBAR,1:M%JBAR),STAT=IZERO)
!               Call ChkMemErr('Read_Evac','STRP%V_DOWN',IZERO) 
!               STRP%U_UP = 0._EB
!               STRP%V_UP = 0._EB
!               STRP%U_DOWN  = 0._EB
!               STRP%V_DOWN  = 0._EB

!               II_C1 = Floor(M%CELLSI(Floor((STRP%XB_CORE(1)-M%XS)*M%RDXINT)) + 1.0_EB)
!               II_C2 = Floor(M%CELLSI(Floor((STRP%XB_CORE(2)-M%XS)*M%RDXINT)) + 1.0_EB)
!               JJ_C1 = Floor(M%CELLSJ(Floor((STRP%XB_CORE(3)-M%YS)*M%RDYINT)) + 1.0_EB)
!               JJ_C2 = Floor(M%CELLSJ(Floor((STRP%XB_CORE(4)-M%YS)*M%RDYINT)) + 1.0_EB)

!               If (RIGHT_HANDED) Then
!                  STRP%U_UP(1    :II_C2,  1    :JJ_C1 ) = 1._EB
!                  STRP%V_UP(II_C2:M%IBAR, 1    :JJ_C2 ) = 1._EB             
!                  STRP%U_UP(II_C1:M%IBAR, JJ_C2:M%JBAR) = -1._EB
!                  STRP%V_UP(1    :II_C1,  JJ_C1:M%JBAR) = -1._EB

!                  STRP%V_DOWN(1    :II_C1,  1    :JJ_C2 )  = 1._EB
!                  STRP%U_DOWN(1    :II_C2,  JJ_C2:M%JBAR)  = 1._EB             
!                  STRP%V_DOWN(II_C2:M%IBAR, JJ_C1:M%JBAR)  = -1._EB
!                  STRP%U_DOWN(II_C1:M%IBAR, 1    :JJ_C1 )  = -1._EB
!               Else
!                  STRP%U_DOWN(1    :II_C2,  1    :JJ_C1 ) = 1._EB
!                  STRP%V_DOWN(II_C2:M%IBAR, 1    :JJ_C2 ) = 1._EB             
!                  STRP%U_DOWN(II_C1:M%IBAR, JJ_C2:M%JBAR) = -1._EB
!                  STRP%V_DOWN(1    :II_C1,  JJ_C1:M%JBAR) = -1._EB

!                  STRP%V_UP(1    :II_C1,  1    :JJ_C2 )  = 1._EB
!                  STRP%U_UP(1    :II_C2,  JJ_C2:M%JBAR)  = 1._EB             
!                  STRP%V_UP(II_C2:M%IBAR, JJ_C1:M%JBAR)  = -1._EB
!                  STRP%U_UP(II_C1:M%IBAR, 1    :JJ_C1 )  = -1._EB
!               End If
!            End If
!         End Do

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
      IF (n_nodes > 0 .AND. MYID==MAX(0,EVAC_PROCESS)) THEN
         n_tmp = 0
         DO n = 1, nmeshes
            IF (evacuation_only(n).AND.evacuation_grid(n)) THEN
               n_tmp = n_tmp + 1
               EVAC_Node_List(n_tmp)%Node_Index = n_tmp
               EVAC_Node_List(n_tmp)%Node_Type  = 'Floor'
               EVAC_Node_List(n_tmp)%ID         = MESH_NAME(n)
               EVAC_Node_List(n_tmp)%GRID_NAME  = MESH_NAME(n)
               EVAC_Node_List(n_tmp)%IMESH      = n
            END IF
         END DO
         DO n = 1, n_entrys
            n_tmp = n_tmp + 1
            evac_entrys(n)%INODE             = n_tmp 
            EVAC_Node_List(n_tmp)%Node_Index = n
            EVAC_Node_List(n_tmp)%Node_Type  = 'Entry'
            EVAC_Node_List(n_tmp)%ID         = EVAC_ENTRYS(n)%ID
            EVAC_Node_List(n_tmp)%IMESH      = EVAC_ENTRYS(n)%IMESH
         END DO
         DO n = 1, n_doors
            n_tmp = n_tmp + 1
            evac_doors(n)%INODE              = n_tmp 
            EVAC_Node_List(n_tmp)%Node_Index = n
            EVAC_Node_List(n_tmp)%Node_Type  = 'Door'
            EVAC_Node_List(n_tmp)%ID         = EVAC_DOORS(n)%ID
            EVAC_Node_List(n_tmp)%IMESH      = EVAC_DOORS(n)%IMESH
         END DO
         DO n = 1, n_exits
            n_tmp = n_tmp + 1
            evac_exits(n)%INODE              = n_tmp 
            EVAC_Node_List(n_tmp)%Node_Index = n
            EVAC_Node_List(n_tmp)%Node_Type  = 'Exit'
            EVAC_Node_List(n_tmp)%ID         = EVAC_EXITS(n)%ID
            EVAC_Node_List(n_tmp)%IMESH      = EVAC_EXITS(n)%IMESH
         END DO
         DO n = 1, n_corrs
            n_tmp = n_tmp + 1
            evac_corrs(n)%INODE              = n_tmp 
            EVAC_Node_List(n_tmp)%Node_Index = n
            EVAC_Node_List(n_tmp)%Node_Type  = 'Corr'
            EVAC_Node_List(n_tmp)%ID         = EVAC_CORRS(n)%ID
         END DO
         DO n = 1, n_strs
            n_tmp = n_tmp + 1
            EVAC_STRS(n)%INODE               = n_tmp
            EVAC_Node_List(n_tmp)%Node_Index = n
            EVAC_Node_List(n_tmp)%Node_Type  = 'Stairs'
            EVAC_Node_List(n_tmp)%ID         = EVAC_STRS(n)%ID
         END DO

         ! Check that door/corr/entry/exit have unique names
         DO n = 1, n_nodes - 1
            DO i = n + 1, n_nodes
               IF (TRIM(EVAC_Node_List(n)%ID) == TRIM(EVAC_Node_List(i)%ID)) THEN
                  WRITE(MESSAGE,'(8A)') 'ERROR: ', TRIM(EVAC_Node_List(n)%Node_Type), ': ', &
                       TRIM(EVAC_Node_List(n)%ID), ' has same ID as ', &
                       TRIM(EVAC_Node_List(i)%Node_Type), ': ', TRIM(EVAC_Node_List(i)%ID)
                  CALL SHUTDOWN(MESSAGE)
               END IF
            END DO
         END DO

         ! BUG fix, 17.11.2008
         DO n = 1, n_entrys
            DO i = 1, evac_entrys(n)%N_VENT_FFIELDS
               IF (evac_entrys(n)%I_DOOR_NODES(i) < 0) THEN
                  evac_entrys(n)%I_DOOR_NODES(i) = EVAC_EXITS(ABS(evac_entrys(n)%I_DOOR_NODES(i)))%INODE
               ELSE IF (evac_entrys(n)%I_DOOR_NODES(i) > 0) THEN
                  evac_entrys(n)%I_DOOR_NODES(i) = EVAC_DOORS(ABS(evac_entrys(n)%I_DOOR_NODES(i)))%INODE
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
      INTEGER nm, i1, i2, j1, j2, NR
      TYPE (EVAC_ENTR_Type), POINTER :: PNX=>NULL()
      TYPE (EVAC_PERS_Type), POINTER :: PCP=>NULL()
      TYPE (EVAC_STRS_Type), POINTER :: STRP=>NULL()
      TYPE (MESH_TYPE), POINTER :: M=>NULL()
      !
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
         !
         !
         CALL CHECKREAD('ENTR',LU_INPUT,IOS)
         IF (IOS == 1) THEN
            EXIT READ_ENTR_LOOP
         END IF
         READ(LU_INPUT,ENTR,END=28,IOSTAT=IOS)
         ! 
         ! Old input used QUANTITY, next lines are needed for that
         IF (MYID==MAX(0,EVAC_PROCESS) .AND. QUANTITY .NE. 'null') WRITE (LU_ERR,'(A,A)') &
              ' WARNING: keyword QUANTITY is replaced by AVATAR_COLOR at ENTR line ',TRIM(ID)
         IF (QUANTITY == 'BLACK')   AVATAR_COLOR = 'BLACK'  
         IF (QUANTITY == 'YELLOW')  AVATAR_COLOR = 'YELLOW' 
         IF (QUANTITY == 'BLUE')    AVATAR_COLOR = 'BLUE'   
         IF (QUANTITY == 'RED')     AVATAR_COLOR = 'RED'    
         IF (QUANTITY == 'GREEN')   AVATAR_COLOR = 'GREEN'  
         IF (QUANTITY == 'MAGENTA') AVATAR_COLOR = 'MAGENTA'
         IF (QUANTITY == 'CYAN')    AVATAR_COLOR = 'CYAN'   
         !
         ! Colors, integer RGB(3), e.g., (23,255,0)
         IF (ANY(RGB < 0) .AND. COLOR=='null') COLOR = 'SKY BLUE'
         IF (COLOR /= 'null') CALL COLOR2RGB(RGB,COLOR)
         IF (ANY(AVATAR_RGB < 0) .AND. AVATAR_COLOR=='null') AVATAR_COLOR = 'ROYAL BLUE 4'
         IF (AVATAR_COLOR /= 'null') CALL COLOR2RGB(AVATAR_RGB,AVATAR_COLOR)
         IF (COLOR_METHOD == 0 .AND. MAX_FLOW > 0.0_EB) THEN
            i_avatar_color = i_avatar_color + 1
            EVAC_AVATAR_RGB(1:3,i_avatar_color) = AVATAR_RGB
         END IF
         IF (MYID /= MAX(0,EVAC_PROCESS)) CYCLE READ_ENTR_LOOP

         IF (MAX_HUMANS < 0) MAX_HUMANS = HUGE(MAX_HUMANS)
         IF (MAX_FLOW <= 0.0_EB) MAX_HUMANS = 0
         IF (MAX_HUMANS_RAMP/='null') THEN
            CALL GET_RAMP_INDEX(MAX_HUMANS_RAMP,'TIME',NR)
            MAX_HUMANS = -NR
         ENDIF 

         PNX=>EVAC_ENTRYS(N)

         PNX%RGB = RGB
         PNX%AVATAR_RGB =AVATAR_RGB
         IF (COLOR_METHOD == 0 .AND. MAX_FLOW > 0.0_EB) PNX%Avatar_Color_Index = i_avatar_color

         IF (EVAC_MESH /= 'null') THEN
            MESH_ID = EVAC_MESH
            IF (MYID==MAX(0,EVAC_PROCESS)) WRITE (LU_ERR,'(A,A)') &
                 ' WARNING: keyword EVAC_MESH is replaced by MESH_ID at ENTR line ', TRIM(ID)
         END IF

         IF (TRIM(KNOWN_DOOR_NAMES(51)) /= 'null') THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: ENTR line ',TRIM(ID), ' problem with KNOWN_DOOR_NAMES'
            CALL SHUTDOWN(MESSAGE)
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
         CALL ChkMemErr('Read_Evac','PNX%I_DOOR_NODES',IZERO) 
         ALLOCATE(PNX%I_VENT_FFIELDS(0:i),STAT=IZERO)
         CALL ChkMemErr('Read_Evac','PNX%I_VENT_FFIELDS',IZERO) 
         ALLOCATE(PNX%P_VENT_FFIELDS(0:i),STAT=IZERO)
         CALL ChkMemErr('Read_Evac','PNX%P_VENT_FFIELDS',IZERO) 
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

         ! Check that the entry is properly specified

         DO I=1,5,2
            IF (XB(I) > XB(I+1)) THEN
               DUMMY   = XB(I)
               XB(I)   = XB(I+1)
               XB(I+1) = DUMMY
            END IF
         END DO
         ! 
         ! Check which evacuation floor
         ii = 0
         n_tmp = 0
         PNX_MeshLoop: DO i = 1, nmeshes
            IF (evacuation_only(i) .AND. evacuation_grid(i)) THEN
               n_tmp = n_tmp + 1
               IF (Is_Within_Bounds(XB(1),XB(2),XB(3),XB(4),XB(5),XB(6),&
                  Meshes(i)%XS,Meshes(i)%XF,Meshes(i)%YS,Meshes(i)%YF,Meshes(i)%ZS,Meshes(i)%ZF, 0._EB, 0._EB, 0._EB)) THEN
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
            CALL SHUTDOWN(MESSAGE)
         END IF
         IF (ii > 1) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: ENTR ',TRIM(ID), ' not an unique mesh found '
            CALL SHUTDOWN(MESSAGE)
         END IF

         nm = PNX%IMESH
         M => MESHES(PNX%IMESH)
 
         IF (XB(1)/=XB(2) .AND. XB(3)/=XB(4)) THEN
            WRITE(MESSAGE,'(A,I4,A)') 'ERROR: ENTR',N,' must be a plane'
            CALL SHUTDOWN(MESSAGE)
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
            WRITE(LU_ERR,fmt='(a,6f12.4)') 'Old XB:', PNX%X1,PNX%X2,PNX%Y1,PNX%Y2,PNX%Z1,PNX%Z2
            WRITE(LU_ERR,fmt='(a,6f12.4)') 'New XB:', XB(1:6)
         END IF

         II = (I1+I2)/2
         JJ = (J1+J2)/2
         SELECT CASE (IOR)
         CASE (+1)
            II = II + 1
         CASE (+2)
            JJ = JJ + 1
         END SELECT
         IF (M%SOLID(M%CELL_INDEX(II,JJ,1))) THEN
            WRITE(LU_ERR,fmt='(a,a,a)') ' WARNING: Entr ',TRIM(ID),' problem with XB, mid point facing solid'
         END IF

         ! Coordinates are lined up with the mesh.
         PNX%X1 = XB(1)
         PNX%X2 = XB(2)
         PNX%Y1 = XB(3)
         PNX%Y2 = XB(4)
         PNX%Z1 = XB(5)
         PNX%Z2 = XB(6)
         ! 
         PNX%IOR        = IOR
         PNX%ID         = ID
         PNX%CLASS_NAME = PERS_ID

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
            IF ( (XB(4)-XB(3)) <= 0.0_EB .OR. (XB(2)-XB(1)) <= 0.0_EB) THEN
               WRITE(MESSAGE,'(A,I4,A)') 'ERROR: ENTR',N,' IOR=3 but not 3-dim object'
               CALL SHUTDOWN(MESSAGE)
            END IF
            PNX%ORIENTATION(3)=-REAL(SIGN(1,IOR),EB)
         CASE (0)
            IF ( (XB(4)-XB(3)) <= 0.0_EB .OR. (XB(2)-XB(1)) <= 0.0_EB) THEN
               WRITE(MESSAGE,'(A,I4,A)') 'ERROR: ENTR',N,' no IOR but not 3-dim object'
               CALL SHUTDOWN(MESSAGE)
            END IF
         CASE Default
            WRITE(MESSAGE,'(A,I4,A)') 'ERROR: ENTR',N,' problem with IOR'
            CALL SHUTDOWN(MESSAGE)
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
               CALL SHUTDOWN(MESSAGE)
            END IF
            PNX%P_VENT_FFIELDS(i) = MAX(0.0_EB,KNOWN_DOOR_PROBS(i))
            PNX%I_VENT_FFIELDS(i) = 0
            PNX%I_DOOR_NODES(i) = 0
            DO j = 1, n_exits
               IF ( TRIM(EVAC_EXITS(j)%ID) == TRIM(KNOWN_DOOR_NAMES(i)) ) THEN
                  PNX%I_VENT_FFIELDS(i) = EVAC_EXITS(j)%I_VENT_FFIELD
                  ! BUG fix, 17.11.2008
                  ! PNX%I_DOOR_NODES(i)   = EVAC_EXITS(j)%INODE
                  PNX%I_DOOR_NODES(i)   = -j
               END IF
            END DO
            DO j = 1, n_doors
               IF ( TRIM(EVAC_DOORS(j)%ID) == TRIM(KNOWN_DOOR_NAMES(i)) ) THEN
                  PNX%I_VENT_FFIELDS(i) = EVAC_DOORS(j)%I_VENT_FFIELD
                  ! BUG fix, 17.11.2008
                  ! PNX%I_DOOR_NODES(i)   = EVAC_DOORS(j)%INODE
                  PNX%I_DOOR_NODES(i)   = +j
               END IF
            END DO
            IF ( PNX%I_VENT_FFIELDS(i)*PNX%I_DOOR_NODES(i) == 0 ) THEN
               WRITE(MESSAGE,'(A,A,A,A,A)') 'ERROR: ENTR line ',TRIM(PNX%ID), &
                    ' problem with door/exit names, ', TRIM(KNOWN_DOOR_NAMES(i)),' not found'
               CALL SHUTDOWN(MESSAGE)
            END IF
         END DO
         !
         ! No known doors given, use the flow_field_id value
         ! 
         PNX%P_VENT_FFIELDS(0) = 1.0_EB
         PNX%I_VENT_FFIELDS(0) = 0
         PNX%I_DOOR_NODES(0) = 0
         PNX_Mesh2Loop: DO i = 1, nmeshes
            IF ( evacuation_only(i) .AND. TRIM(PNX%GRID_NAME) == TRIM(MESH_NAME(i)) ) THEN
               PNX%I_VENT_FFIELDS(0) = i
               EXIT PNX_Mesh2Loop
            END IF
         END DO PNX_Mesh2Loop
         IF ( PNX%I_VENT_FFIELDS(0) == 0 ) THEN
            WRITE(MESSAGE,'(A,A,A,A,A)') 'ERROR: ENTR line ',TRIM(PNX%ID),&
                 ' problem with flow field name, ', TRIM(PNX%GRID_NAME),' not found'
            CALL SHUTDOWN(MESSAGE)
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
      TYPE (EVACUATION_Type), POINTER :: HPT=>NULL()
      TYPE (EVAC_PERS_Type),  POINTER :: PCP=>NULL()

      READ_EVAC_LOOP: DO N=1,NPC_EVAC
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
         SAMPLING_FACTOR          = 1      
         NUMBER_INITIAL_PERSONS   = 0
         XB                       = 0.0_EB
         ANGLE                    = -1000.0_EB
         EVACFILE                 = .FALSE.
         AFTER_REACTION_TIME      = .FALSE.
         SHOW                     = .TRUE.
         TIME_START               = -99.0_EB
         TIME_STOP                = -99.0_EB
         GN_MIN                   = 1
         GN_MAX                   = 1      

         KNOWN_DOOR_NAMES         = 'null'
         KNOWN_DOOR_PROBS         = 1.0_EB
         !
         CALL CHECKREAD('EVAC',LU_INPUT,IOS)
         IF (IOS == 1) THEN
            EXIT READ_EVAC_LOOP
         END IF
         READ(LU_INPUT,EVAC,END=25,IOSTAT=IOS)
         ! 
         ! Old input used QUANTITY, next lines are needed for that
         IF (MYID==MAX(0,EVAC_PROCESS) .AND. QUANTITY .NE. 'null') WRITE (LU_ERR,'(A,A)') &
              ' WARNING: keyword QUANTITY is replaced by AVATAR_COLOR at EVAC line ',TRIM(ID)
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
            i_avatar_color = i_avatar_color + 1
            EVAC_AVATAR_RGB(1:3,i_avatar_color) = AVATAR_RGB
         END IF

         IF (MYID /= MAX(0,EVAC_PROCESS)) CYCLE READ_EVAC_LOOP

         HPT=>EVACUATION(N)

         HPT%RGB = RGB
         HPT%AVATAR_RGB = AVATAR_RGB
         IF (COLOR_METHOD == 0 .AND. NUMBER_INITIAL_PERSONS > 0) HPT%Avatar_Color_Index = i_avatar_color
         IF (EVAC_MESH /= 'null') THEN
            MESH_ID = EVAC_MESH
            IF (MYID==MAX(0,EVAC_PROCESS)) WRITE (LU_ERR,'(A,A)') &
                 ' WARNING: keyword EVAC_MESH is replaced by MESH_ID at EVAC line ', TRIM(ID)
         END IF
         IF (TRIM(PERS_ID) == 'null') THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: EVAC line ',TRIM(ID),' no PERS_ID given'
            CALL SHUTDOWN(MESSAGE)
         ELSE
            ii = 1
            DO i = 1,NPC_PERS
               IF (TRIM(EVAC_PERSON_CLASSES(i)%ID) == TRIM(PERS_ID)) CYCLE
               ii = ii + 1
            END DO
            IF (ii > NPC_PERS) THEN
               WRITE(MESSAGE,'(A,A,A)') 'ERROR: EVAC line ',TRIM(ID), ' prblem with PERS_ID'
               CALL SHUTDOWN(MESSAGE)
            END IF
         END IF

         IF (TRIM(KNOWN_DOOR_NAMES(51)) /= 'null') THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: EVAC line ',TRIM(ID), ' problem with KNOWN_DOOR_NAMES'
            CALL SHUTDOWN(MESSAGE)
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

         ! Check which evacuation floor
         ii = 0
         HP_MeshLoop: DO i = 1, nmeshes
            IF (evacuation_only(i) .AND. evacuation_grid(i)) THEN
               IF ( Is_Within_Bounds(HPT%X1,HPT%X2,HPT%Y1,HPT%Y2,HPT%Z1,HPT%Z2,&
                    Meshes(i)%XS,Meshes(i)%XF,Meshes(i)%YS,Meshes(i)%YF,Meshes(i)%ZS,Meshes(i)%ZF, 0._EB, 0._EB, 0._EB)) THEN
                  IF (TRIM(MESH_ID) == 'null' .OR. TRIM(MESH_ID) == TRIM(MESH_NAME(i))) THEN
                     ii = ii + 1
                     HPT%IMESH = i
                  END IF
               END IF
            END IF
         END DO HP_MeshLoop
         IF (HPT%IMESH == 0) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: EVAC line ',TRIM(ID),' problem with IMESH, no mesh found'
            CALL SHUTDOWN(MESSAGE)
         END IF
         IF (ii > 1) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: EVAC line ',TRIM(ID), ' not an unique mesh found '
            CALL SHUTDOWN(MESSAGE)
         END IF

         ! Use the main_evac_grid flow field if none is given
         IF (TRIM(FLOW_FIELD_ID) == 'null') THEN
            HPT%GRID_NAME  = TRIM(MESH_NAME(HPT%IMESH))
         ELSE
            HPT%GRID_NAME  = FLOW_FIELD_ID
         END IF
         !
         DO i = 1, HPT%N_VENT_FFIELDS
            HPT%P_VENT_FFIELDS(i) = MAX(0.0_EB,KNOWN_DOOR_PROBS(i))
            HPT%I_VENT_FFIELDS(i) = 0
            HPT%I_DOOR_NODES(i) = 0
            DO j = 1, n_exits
               IF ( TRIM(EVAC_EXITS(j)%ID) == TRIM(KNOWN_DOOR_NAMES(i)) ) THEN
                  HPT%I_VENT_FFIELDS(i) = EVAC_EXITS(j)%I_VENT_FFIELD
                  HPT%I_DOOR_NODES(i)   = EVAC_EXITS(j)%INODE
               END IF
            END DO
            DO j = 1, n_doors
               IF ( TRIM(EVAC_DOORS(j)%ID) == TRIM(KNOWN_DOOR_NAMES(i)) ) THEN
                  HPT%I_VENT_FFIELDS(i) = EVAC_DOORS(j)%I_VENT_FFIELD
                  HPT%I_DOOR_NODES(i)   = EVAC_DOORS(j)%INODE
               END IF
            END DO
            IF ( HPT%I_VENT_FFIELDS(i)*HPT%I_DOOR_NODES(i) == 0 ) THEN
               WRITE(MESSAGE,'(A,A,A,A,A)') 'ERROR: EVAC line ',TRIM(HPT%ID), &
                    ' problem with door/exit names, ', TRIM(KNOWN_DOOR_NAMES(i)),' not found'
               CALL SHUTDOWN(MESSAGE)
            END IF
         END DO
         !
         ! No known doors given, use the flow_field_id value
         ! 
         HPT%P_VENT_FFIELDS(0) = 1.0_EB
         HPT%I_VENT_FFIELDS(0) = 0
         HPT%I_DOOR_NODES(0) = 0
         HP_Mesh2Loop: DO i = 1, nmeshes
            IF ( evacuation_only(i) .AND. TRIM(HPT%GRID_NAME) == TRIM(MESH_NAME(i)) ) THEN
               HPT%I_VENT_FFIELDS(0) = i
               EXIT HP_Mesh2Loop
            END IF
         END DO HP_Mesh2Loop
         IF ( HPT%I_VENT_FFIELDS(0) == 0 ) THEN
            WRITE(MESSAGE,'(A,A,A,A,A)') 'ERROR: EVAC line ',TRIM(HPT%ID), &
                 ' problem with flow field name, ', TRIM(HPT%GRID_NAME),' not found'
            CALL SHUTDOWN(MESSAGE)
         END IF
         !
      END DO READ_EVAC_LOOP
25    REWIND(LU_INPUT)

      N_EVAC = 1
      ALLOCATE(EVAC_CLASS_NAME(N_EVAC),STAT=IZERO)
      CALL ChkMemErr('READ_EVAC','EVAC_CLASS_NAME',IZERO)
      ALLOCATE(EVAC_CLASS_RGB(3,N_EVAC),STAT=IZERO)
      CALL ChkMemErr('READ_EVAC','EVAC_CLASS_RGB',IZERO)
      EVAC_CLASS_NAME(1) = 'Human'
      DO N = 1, N_EVAC
         EVAC_CLASS_RGB(1:3,N) = (/ 39, 64,139/)  ! ROYAL BLUE 4
      END DO
      ! Default color table for agents

    END SUBROUTINE READ_EVAC_LINES

    SUBROUTINE READ_EVHO
      IMPLICIT NONE
      !
      ! Read the EVHO lines
      !
      ! Local variables
      TYPE (EVAC_HOLE_Type),  POINTER :: EHX=>NULL()

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
         !
         CALL CHECKREAD('EVHO',LU_INPUT,IOS)
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
         ! 
         ! Check which evacuation floor
         ii = 0
         EHX_MeshLoop: DO i = 1, nmeshes
            IF (evacuation_only(i) .AND. evacuation_grid(i)) THEN
               IF ( (EHX%Z1 >= Meshes(i)%ZS .AND. EHX%Z2 <= Meshes(i)%ZF).AND. &
                    (EHX%Y1 >= Meshes(i)%YS .AND. EHX%Y2 <= Meshes(i)%YF).AND. &
                    (EHX%X1 >= Meshes(i)%XS .AND. EHX%X2 <= Meshes(i)%XF)) THEN
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
            CALL SHUTDOWN(MESSAGE)
         END IF
         IF (ii > 1) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: EVHO line ',TRIM(EHX%ID), ' not an unique mesh found '
            CALL SHUTDOWN(MESSAGE)
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
      TYPE (EVAC_SSTAND_Type), POINTER :: ESS=>NULL()
      REAL(EB) :: X, Y, Z

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
         CALL CHECKREAD('EVSS',LU_INPUT,IOS)
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
         ! 
         ! Check which evacuation floor
         ii = 0
         ESS_MeshLoop: DO i = 1, nmeshes
            IF (evacuation_only(i) .AND. evacuation_grid(i)) THEN
               IF ( (ESS%Z1 >= Meshes(i)%ZS .AND. ESS%Z2 <= Meshes(i)%ZF).AND. &
                    (ESS%Y1 >= Meshes(i)%YS .AND. ESS%Y2 <= Meshes(i)%YF).AND. &
                    (ESS%X1 >= Meshes(i)%XS .AND. ESS%X2 <= Meshes(i)%XF)) THEN
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
            CALL SHUTDOWN(MESSAGE)
         END IF
         IF (ii > 1) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: EVSS line ',TRIM(ESS%ID), ' not an unique mesh found '
            CALL SHUTDOWN(MESSAGE)
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
            WRITE(MESSAGE,'(A,I4,A)') 'ERROR: EVSS',N,' problem with IOR'
            CALL SHUTDOWN(MESSAGE)
         END SELECT

         ! Colors, integer RGB(3), e.g., (23,255,0)
         ! If (Any(RGB < 0) .And. COLOR=='null') COLOR = 'BLACK'
         IF (ANY(RGB < 0) .AND. COLOR=='null') COLOR = 'AZURE 2'  ! RGB = 193 205 205
         IF (COLOR /= 'null') CALL COLOR2RGB(RGB,COLOR)
         ESS%RGB = RGB

         ! Check if doors, exits and entrys are on the inclines.
         ! If they are, change the 'z' so that Smokeview plots correctly.
         DO I = 1, N_EXITS
            IF (EVAC_EXITS(I)%IMESH /= ESS%IMESH) CYCLE
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
      TYPE (EVAC_ENTR_Type), POINTER :: PNX=>NULL()
      TYPE (EVAC_DOOR_Type), POINTER :: PDX=>NULL()
      TYPE (EVAC_STRS_Type), POINTER :: STRP=>NULL()
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
            WRITE(MESSAGE,'(A,I4,A)') 'ERROR: CORR',n,' problem with TO_NODE, loop2'
            CALL SHUTDOWN(MESSAGE)
         END IF
      END DO

      !
      DO n = 1, n_doors
         NodeLoop: DO i = 1, n_nodes
            IF (EVAC_Node_List(i)%ID == EVAC_DOORS(n)%TO_NODE) THEN
               EVAC_DOORS(n)%INODE2 = i
               EVAC_DOORS(n)%IMESH2 = EVAC_Node_List(i)%IMESH
               EXIT NodeLoop
            END IF
         END DO NodeLoop
         IF (EVAC_DOORS(n)%INODE2 == 0 .OR. EVAC_DOORS(n)%IMESH2 == 0) THEN
            WRITE(MESSAGE,'(A,I4,A)') 'ERROR: DOOR',n,' problem with TO_NODE'
            CALL SHUTDOWN(MESSAGE)
         END IF
      END DO


      DO n = 1, n_doors
         IF (EVAC_DOORS(n)%KEEP_XY) THEN
            i = EVAC_DOORS(n)%INODE2
            IF (TRIM(EVAC_Node_List(i)%Node_Type) == 'Door') THEN
               PDX => EVAC_DOORS(EVAC_Node_List(i)%Node_Index)
               IF ((EVAC_DOORS(n)%IOR /= -PDX%IOR) .OR. ABS(EVAC_DOORS(n)%Width-PDX%Width) > 0.1_EB ) THEN
                  WRITE(MESSAGE,'(A,I4,A)') 'ERROR: DOOR',N,' KEEP_XY Problem'
                  CALL SHUTDOWN(MESSAGE)
               END IF
            END IF
            IF (TRIM(EVAC_Node_List(i)%Node_Type) == 'Entry') THEN
               PNX => EVAC_ENTRYS(EVAC_Node_List(i)%Node_Index)
               IF ((EVAC_DOORS(n)%IOR /= PNX%IOR) .OR. ABS(EVAC_DOORS(n)%Width-PNX%Width) > 0.1_EB ) THEN
                  WRITE(MESSAGE,'(A,I4,A)') 'ERROR: DOOR',N,' KEEP_XY Problem'
                  CALL SHUTDOWN(MESSAGE)
               END IF
            END IF
         END IF
      END DO

   DO N = 1,N_DOORS
      PDX => EVAC_DOORS(N)
      ! Check if door leads to Stairs
       PDX%STR_INDX = 0
       PDX%STR_SUB_INDX = 0
       CheckDoorStrLoop: DO i = 1, N_STRS
         STRP=>EVAC_STRS(i)
         IF (STRP%IMESH==PDX%IMESH .OR. STRP%IMESH==PDX%IMESH2) THEN
            PDX%STR_INDX = i
            DO j = 1,STRP%N_NODES
               IF ( Is_Within_Bounds(PDX%X1,PDX%X2,PDX%Y1,PDX%Y2,PDX%Z1,PDX%Z2, &
                    STRP%XB_NODE(j,1), STRP%XB_NODE(j,2), STRP%XB_NODE(j,3),STRP%XB_NODE(j,4), &
                    STRP%XB_NODE(j,5), STRP%XB_NODE(j,6), 0._EB, 0._EB, 0._EB)) THEN
                  PDX%STR_SUB_INDX = j
                  EXIT CheckDoorStrLoop
               END IF
            END DO
         END IF
       END DO CheckDoorStrLoop
    END DO

    ! Create list  of incoming nodes for STairs
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
      DO I = NMESHES+1,N_NODES
         IF (EVAC_NODE_List(I)%Node_type == 'Entry') CYCLE
         IF (EVAC_NODE_List(I)%Node_type == 'Stair') CYCLE
         IF (EVAC_NODE_List(I)%IMESH == STRP%IMESH) THEN
            STRP%N_NODES_OUT = STRP%N_NODES_OUT + 1
            NODES_TMP(STRP%N_NODES_OUT) = I
         END IF
      END DO
      ALLOCATE(STRP%NODES_OUT(1:STRP%N_NODES_OUT),STAT=IZERO)
      CALL ChkMemErr('Read_Evac','STRP%NODES_OUT',IZERO) 
      STRP%NODES_OUT = NODES_TMP(1:STRP%N_NODES_OUT)
    END DO

    END SUBROUTINE CHECK_EVAC_NODES

  END SUBROUTINE READ_EVAC

  SUBROUTINE Initialize_Evac_Dumps
    IMPLICIT NONE
    !
    ! Local variables
    CHARACTER(50) tcform
    INTEGER n_cols, i, j, nm, izero
    LOGICAL L_fed_read, L_fed_save, L_eff_read, L_eff_save, L_status
    INTEGER(4) n_egrids_tmp, ibar_tmp, jbar_tmp, kbar_tmp, &
         ntmp1, ntmp2, ntmp3, ntmp4, ntmp5, ntmp6, ios, N
    REAL(FB) u_tmp, v_tmp
    CHARACTER(60), ALLOCATABLE, DIMENSION(:) :: CTEMP
    !
    TYPE (MESH_TYPE), POINTER :: MFF =>NULL()
    !

    ! Logical unit numbers
    ! LU_EVACCSV: CHID_evac.csv, number of persons
    ! LU_EVACEFF: CHID_evac.eff, evacflow fields, binary
    ! LU_EVACFED: CHID_evac.fed, FED and soot, time dependent, binary
    !      Format: 1. row: n_egrids >=0  (Old Format, version 1.10)
    !              1a. row: n < 0 (New Format)
    !              1b. row: n_egrids,4,n_corrs=0,4 (New Format, version 1.11)

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

    !

    ! Open evacuation output file

    IF (APPEND) THEN
       OPEN (LU_EVACOUT,file=FN_EVACOUT,form='formatted',status='old', position='append')
    ELSE
       OPEN (LU_EVACOUT,file=FN_EVACOUT,form='formatted', status='replace')
    ENDIF

    ! Write program info

    WRITE(EVAC_COMPILE_DATE,'(A)') evacrev(INDEX(evacrev,':')+1:LEN_TRIM(evacrev)-2)
    READ (EVAC_COMPILE_DATE,'(I5)') EVAC_MODULE_REV
    WRITE(EVAC_COMPILE_DATE,'(A)') evacdate
    CALL GET_REV_evac(EVAC_MODULE_REV,EVAC_COMPILE_DATE)
    !
    WRITE(LU_EVACOUT,'(/A)')          ' FDS+Evac Evacuation Module'
    WRITE(LU_EVACOUT,'(/A,A)')        ' FDS+Evac Compilation Date: ', &
         TRIM(EVAC_COMPILE_DATE(INDEX(EVAC_COMPILE_DATE,'(')+1:INDEX(EVAC_COMPILE_DATE,')')-1))
    WRITE(LU_EVACOUT,'(A,A)')  ' FDS+Evac Version         : ', TRIM(EVAC_VERSION)
    WRITE(LU_EVACOUT,'(A,i0)')  ' FDS+Evac SVN Revision No.: ', EVAC_MODULE_REV

    WRITE(LU_EVACOUT,fmt='(/a,i2)')  ' FDS+Evac Color_Method    :', COLOR_METHOD
    IF (Fed_Door_Crit >= 0) THEN
       WRITE(LU_EVACOUT,fmt='(a,f14.8)') ' FDS+Evac Fed_Door_Crit   :', FED_DOOR_CRIT
    ELSE
       ! Visibility S = 3/K, K is extinction coeff.
       WRITE(LU_EVACOUT,fmt='(a,f14.8,a)') ' FDS+Evac Vis_Door_Crit   :', ABS(FED_DOOR_CRIT), ' m'
    END IF
    IF (NOT_RANDOM ) WRITE(LU_EVACOUT,fmt='(a)') ' FDS+Evac Random seed is not used.'
    IF (Fed_Door_Crit < 0) THEN
       ! Visibility S = 3/K, K is extinction coeff.
       FED_DOOR_CRIT = 3.0_EB/FED_DOOR_CRIT ! Extinction coeff (1/m)
    END IF
    !
    L_fed_read = BTEST(I_EVAC,3)
    L_fed_save = BTEST(I_EVAC,1)
    L_eff_read = BTEST(I_EVAC,2)
    L_eff_save = BTEST(I_EVAC,0)

    n_cols = n_egrids + n_corrs + n_exits + n_doors + 1 + n_exits - n_co_exits + n_doors
    ! Initialize the FED counters:
    icyc_old = -1
    n_dead = -1
    fed_max_alive = 0.0_EB
    fed_max = 0.0_EB
    !
    IF (append) THEN
       OPEN (LU_EVACCSV,file=FN_EVACCSV,form='formatted',status='old', position='append')
       !
       IF (L_fed_save) THEN
          OPEN (LU_EVACFED,file=FN_EVACFED,form='unformatted', status='old',position='append')
       END IF
       ! 
       IF (L_fed_read) THEN
          INQUIRE (file=FN_EVACFED,exist=L_status)
          IF (.NOT. L_status) THEN
             WRITE (LU_EVACOUT,fmt='(a,a,a)') ' FDS+Evac No FED File: ', TRIM(FN_EVACFED), ', FED and soot not used'
             l_fed_read = .FALSE.
             l_fed_save = .FALSE.
             I_EVAC = IBCLR(I_EVAC,3)  ! do not read FED
             I_EVAC = IBCLR(I_EVAC,1)  ! do not save FED
          ELSE
             CALL SHUTDOWN('ERROR: Evac Dumps: FED, no restart yet')
             OPEN (LU_EVACFED,file=FN_EVACFED,form='unformatted', status='old')
             READ (LU_EVACFED,Iostat=ios) n_egrids_tmp
             IF (ios.NE.0) THEN
                WRITE(MESSAGE,'(A)') 'ERROR: Init Evac Dumps: FED READ ERROR'
                CLOSE (LU_EVACFED)
                CALL SHUTDOWN(MESSAGE)
             END IF
             IF (n_egrids_tmp /= n_egrids) THEN
                WRITE(MESSAGE,'(A,2I4)') 'ERROR: Init Evac Dumps: FED ',n_egrids_tmp, n_egrids
                CALL SHUTDOWN(MESSAGE)
             END IF
          END IF
       END IF
       IF (L_eff_read) THEN
          INQUIRE (file=FN_EVACEFF,exist=L_status)
          IF (L_status) THEN
             WRITE (LU_EVACOUT,fmt='(a,a,a/)') ' FDS+Evac EFF File: ', TRIM(FN_EVACEFF), ' is used'
             l_eff_save = .FALSE.
             I_EVAC = IBCLR(I_EVAC,0)  ! do not save EFF
             OPEN (LU_EVACEFF,file=FN_EVACEFF,form='unformatted', status='old')
             READ (LU_EVACEFF,Iostat=ios) n_egrids_tmp
             IF (ios.NE.0) THEN
                WRITE(MESSAGE,'(A)') 'ERROR: Init Evac Dumps: EFF READ ERROR'
                CLOSE (LU_EVACEFF)
                CALL SHUTDOWN(MESSAGE)
             END IF
             IF (n_egrids_tmp /= COUNT(EVACUATION_ONLY)) THEN
                WRITE(MESSAGE,'(A,2I4)') 'ERROR: Init Evac Dumps: EFF ',n_egrids_tmp, COUNT(EVACUATION_ONLY)
                CLOSE (LU_EVACEFF)
                CALL SHUTDOWN(MESSAGE)
             END IF
          ELSE
             WRITE(MESSAGE,'(A,2I4)') 'ERROR: Init Evac Dumps: EFF, no restart yet'
             CALL SHUTDOWN(MESSAGE)
          END IF
       END IF
       !
       IF ( l_fed_read .OR. l_fed_save ) n_dead = 0
       !

    ELSE                      ! replace files
       !
       IF (L_fed_save) THEN
          l_fed_read = .FALSE.
          I_EVAC = IBCLR(I_EVAC,3)  ! do not read FED
          OPEN (LU_EVACFED,file=FN_EVACFED,form='unformatted', status='replace')
          ! First line: <0 new format
          !             -1: second line: #mesh #reals #corrs #reals #doors+exits #nreals
          !              (#reals: fed,soot,temp,radflux,...)
          ! First line: >0: nmeshes, fed and soot saved/read for meshes
          ntmp1 = -1
          ntmp2 = 4
          ntmp3 = N_CORRS
          ! Corrs: save for both XB1 and XB2 (if only XB, XB2 is then zeros)
          ntmp4 = 8
          ntmp5 = N_DOORS+N_EXITS
          ntmp6 = 4
          n_egrids_tmp = n_egrids
          WRITE (LU_EVACFED) ntmp1
          WRITE (LU_EVACFED) n_egrids_tmp, ntmp2, ntmp3, ntmp4, ntmp5, ntmp6
          WRITE (LU_EVACOUT,fmt='(a,a,a)') ' FDS+Evac FED File: ', TRIM(FN_EVACFED), ' is calculated and used'
       END IF
       ! 
       IF (L_fed_read) THEN
          INQUIRE (file=FN_EVACFED,exist=L_status)
          IF (.NOT. L_status) THEN
             WRITE (LU_EVACOUT,fmt='(a,a,a)') ' FDS+Evac No FED File: ', TRIM(FN_EVACFED), ', FED and soot not used'
             l_fed_read = .FALSE.
             l_fed_save = .FALSE.
             I_EVAC = IBCLR(I_EVAC,3)  ! do not read FED
             I_EVAC = IBCLR(I_EVAC,1)  ! do not save FED
          ELSE
             OPEN (LU_EVACFED,file=FN_EVACFED,form='unformatted', status='old')
             READ (LU_EVACFED,Iostat=ios) ntmp1
             IF (ios.NE.0) THEN
                WRITE(MESSAGE,'(A)') 'ERROR: Init Evac Dumps: FED READ ERROR'
                CLOSE (LU_EVACFED)
                CALL SHUTDOWN(MESSAGE)
             END IF
             IF ( ntmp1 >= 0 ) THEN
                ! Old format (version 1.10)
                n_egrids_tmp = ntmp1
             ELSE
                ! New format (version 1.11)
                READ (LU_EVACFED,Iostat=ios) n_egrids_tmp, ntmp2, ntmp3, ntmp4, ntmp5, ntmp6
                IF (ios.NE.0) THEN
                   WRITE(MESSAGE,'(A)') 'ERROR: Init Evac Dumps: FED READ ERROR'
                   CLOSE (LU_EVACFED)
                   CALL SHUTDOWN(MESSAGE)
                END IF
             END IF

             ! Do not read old format. Do not read new format, if there the numbers
             ! are not: n_egrids, 4, n_corrs, 8 
             IF ( ntmp2 /= 4 .OR. ntmp3 /= n_corrs .OR. ntmp1 >= 0 .OR. ntmp4 /= 8  .OR. &
                  ntmp5 /= n_doors+n_exits .OR. ntmp6 /= 4) THEN
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
                CALL SHUTDOWN(MESSAGE)
             END IF

             IF (l_fed_read .OR. l_fed_save) THEN
                WRITE (LU_EVACOUT,fmt='(a,a,a)') ' FDS+Evac FED File: ', TRIM(FN_EVACFED), ' is used'
             END IF
          END IF
       END IF
       ! 
       ! Number of evac flow fields is same as the number of all evac grids.
       IF (L_eff_read) THEN
          ios = 3
          INQUIRE (file=FN_EVACEFF,exist=L_status)
          IF (L_status) THEN
             l_eff_save = .FALSE.
             l_eff_read = .TRUE.
             I_EVAC = IBCLR(I_EVAC,0) ! do not save EFF
             I_EVAC = IBSET(I_EVAC,2) ! read EFF
             OPEN (LU_EVACEFF,file=FN_EVACEFF,form='unformatted', status='old')
             READ (LU_EVACEFF,Iostat=ios) n_egrids_tmp
             IF (ios.NE.0) THEN
                ios = 1
                WRITE(LU_EVACOUT,'(A)') ' WARNING: Init Evac Dumps: EFF READ ERROR'
                WRITE(LU_EVACOUT,'(A)') ' WARNING: EFF file is not read in'
                CLOSE (LU_EVACEFF)
             END IF
             IF (n_egrids_tmp /= COUNT(EVACUATION_ONLY) .AND. ios < 1) THEN
                ios = 2
                WRITE(LU_EVACOUT,'(A,2I4)') ' WARNING: Init Evac Dumps: EFF READ ERROR ',n_egrids_tmp, COUNT(EVACUATION_ONLY)
                WRITE(LU_EVACOUT,'(A)')     ' WARNING: EFF file is not read in'
                CLOSE (LU_EVACEFF)
             END IF
          END IF
          IF (ios .NE. 0) THEN
             ! Read error ==> recalculate EFF
             l_eff_save = .TRUE.
             l_eff_read = .FALSE.
             I_EVAC = IBCLR(I_EVAC,2) ! do not read EFF
             I_EVAC = IBSET(I_EVAC,0) ! save EFF
          END IF
       END IF  ! L_eff_read
       !
       ALLOCATE(CTEMP(MAX(1,n_exits-n_co_exits)), STAT = IZERO)
       CALL ChkMemErr('Initialize_Evac_Dumps','CTEMP', IZERO)
       j = 0
       DO i = 1, n_exits
          IF (.NOT. EVAC_EXITS(i)%COUNT_ONLY) THEN
             j = j + 1
             CTEMP(j) = TRIM(EVAC_EXITS(i)%ID)
          END IF
       END DO
    
       IF ( l_fed_read .OR. l_fed_save ) THEN
          ! Write the 'fed' columns
          n_dead = 0
          OPEN (LU_EVACCSV,file=FN_EVACCSV,form='formatted',status='replace')
          ! June 2009: Changed the .csv file format to the fds5 style
          ! first row: units (or variable class)
          ! second row: variable name
          ! third row-: data
          ! Write (LU_EVACCSV,*) n_cols+3
          WRITE (tcform,'(a,i4.4,a)') "(",n_cols+3,"(a,','),a)"
          ! Write (LU_EVACCSV,tcform) 'Time','Humans', &
          !      ('Floor', i=1,n_egrids), &
          !      ('Corridor', i=1,n_corrs), &
          !      ('Exit', i=1,n_exits), &
          !      ('Door', i=1,n_doors), &
          !      ('Exit', i=1,n_exits-n_co_exits), &
          !      ('Door', i=1,n_doors), &
          !      'Fed','Fed','Fed'
          WRITE (LU_EVACCSV,tcform) 's','AgentsInside', &
               ('AgentsInsideMesh', i=1,n_egrids), &
               ('AgentsInsideCorr', i=1,n_corrs), &
               ('ExitCounter', i=1,n_exits), &
               ('DoorCounter', i=1,n_doors), &
               ('TargetExitCounter', i=1,n_exits-n_co_exits), &
               ('TargetDoorCounter', i=1,n_doors), &
               'Agents','FED_Index','FED_Index'
          WRITE (LU_EVACCSV,tcform) 'EVAC_Time','AllAgents', &
               (TRIM(EVAC_Node_List(i)%GRID_NAME), i=1,n_egrids), &
               (TRIM(EVAC_CORRS(i)%ID), i=1,n_corrs), &
               (TRIM(EVAC_EXITS(i)%ID), i=1,n_exits), &
               (TRIM(EVAC_DOORS(i)%ID), i=1,n_doors), &
               (TRIM(CTEMP(i)), i=1,n_exits-n_co_exits), &
               (TRIM(EVAC_DOORS(i)%ID), i=1,n_doors), &
               'Number_of_Deads','FED_max','FED_max_alive'
       ELSE
          ! Do not write the 'fed' columns
          OPEN (LU_EVACCSV,file=FN_EVACCSV,form='formatted',status='replace')
          ! June 2009: Changed the .csv file format to the fds5 style
          ! first row: units (or variable class)
          ! second row: variable name
          ! third row-: data
          ! Write (LU_EVACCSV,*) n_cols
          WRITE (tcform,'(a,i4.4,a)') "(",n_cols,"(a,','),a)"
          ! Write (LU_EVACCSV,tcform) 'Time','Humans', &
          !      ('Floor', i=1,n_egrids), &
          !      ('Corridor', i=1,n_corrs), &
          !      ('Exit', i=1,n_exits), &
          !      ('Door', i=1,n_doors), &
          !      ('Exit', i=1,n_exits-n_co_exits), &
          !      ('Door', i=1,n_doors)
          WRITE (LU_EVACCSV,tcform) 's','AgentsInside', &
               ('AgentsInsideMesh', i=1,n_egrids), &
               ('AgentsInsideCorr', i=1,n_corrs), &
               ('ExitCounter', i=1,n_exits), &
               ('DoorCounter', i=1,n_doors), &
               ('TargetExitCounter', i=1,n_exits-n_co_exits), &
               ('TargetDoorCounter', i=1,n_doors)
          WRITE (LU_EVACCSV,tcform) 'EVAC_Time','AllAgents', &
               (TRIM(EVAC_Node_List(i)%GRID_NAME), i=1,n_egrids), &
               (TRIM(EVAC_CORRS(i)%ID), i=1,n_corrs), &
               (TRIM(EVAC_EXITS(i)%ID), i=1,n_exits), &
               (TRIM(EVAC_DOORS(i)%ID), i=1,n_doors), &
               (TRIM(CTEMP(i)), i=1,n_exits-n_co_exits), &
               (TRIM(EVAC_DOORS(i)%ID), i=1,n_doors)
       END IF
       DEALLOCATE(CTEMP)
    END IF                  ! if-append-else
    ! 
    ! Read the evac flow fields from the disk, if they exist.
    IF ( L_eff_read ) THEN
       ios = 0
       ReadEffLoop: DO nm = 1, NMESHES
          IF (EVACUATION_ONLY(NM)) THEN
             MFF=>MESHES(nm)
             READ (LU_EVACEFF,Iostat=ios) ibar_tmp, jbar_tmp, kbar_tmp
             IF (ios .NE. 0) THEN
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
                   IF (ios .NE. 0) THEN
                      WRITE (LU_EVACOUT,'(A)') ' WARNING: Init Evac Dumps: EFF READ ERROR'
                      WRITE (LU_EVACOUT,'(A)') ' WARNING: EFF file is not read in'
                      CLOSE (LU_EVACEFF)
                      EXIT ReadEffLoop
                   END IF
                   MFF%U(i,j,:) = u_tmp
                   MFF%V(i,j,:) = v_tmp
                   MFF%W(i,j,:) = 0.0_EB
                END DO
             END DO
             MFF%UVW_GHOST(:,1)=-1.E6_EB
             MFF%UVW_GHOST(:,2)=-1.E6_EB
          END IF
       END DO ReadEffLoop
       IF (ios .NE. 0) THEN
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
    IF (L_eff_save .AND. .NOT. append) THEN
       l_eff_read = .FALSE.
       I_EVAC = IBCLR(I_EVAC,2)  ! do not read EFF
       OPEN (LU_EVACEFF,file=FN_EVACEFF,form='unformatted', status='replace')
       n_egrids_tmp = COUNT(EVACUATION_ONLY)
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

    ! write STRS properties
    IF (MYID==MAX(0,EVAC_PROCESS)) THEN
       DO N = 1, N_STRS
          WRITE (LU_EVACOUT,'(A,A)')      '  Stair ',TRIM(EVAC_STRS(N)%ID)
          WRITE (LU_EVACOUT,'(A,6F10.3)') '   Co-ordinates: ',EVAC_STRS(N)%XB(1:6)
          WRITE (LU_EVACOUT,'(A)')          '   Node coordinates'
          DO NM = 1,EVAC_STRS(N)%N_NODES
             IF (EVAC_STRS(N)%NODE_TYPE(NM) == STRS_LANDING_TYPE) THEN
                WRITE (LU_EVACOUT,'(A,I3,8F8.3)') '   Landing ',NM,EVAC_STRS(N)%XB_NODE(NM,1:8)
             ELSE
                WRITE (LU_EVACOUT,'(A,I3,8F8.3)') '   Stair   ',NM,EVAC_STRS(N)%XB_NODE(NM,1:8)
             ENDIF
          ENDDO
          WRITE (LU_EVACOUT,'(A)')          '   Nodes in '
          DO NM = 1, EVAC_STRS(N)%N_NODES_IN
             WRITE (LU_EVACOUT,'(I5,A,A)')          NM, ' ', EVAC_NODE_List(EVAC_STRS(N)%NODES_IN(NM))%ID
          ENDDO
          WRITE (LU_EVACOUT,'(A)')          '   Nodes out '
          DO NM = 1, EVAC_STRS(N)%N_NODES_OUT
             WRITE (LU_EVACOUT,'(I5,A,A)')          NM, ' ', EVAC_NODE_List(EVAC_STRS(N)%NODES_OUT(NM))%ID
          ENDDO
       ENDDO
   END IF

  END SUBROUTINE Initialize_Evac_Dumps
      
!
  SUBROUTINE INITIALIZE_EVACUATION(NM,ISTOP)
    IMPLICIT NONE
    !
    ! Insert humans into the domain at the start of calculation
    !
    ! Passed variables
    INTEGER, INTENT(IN) :: NM
    INTEGER, INTENT(OUT) :: ISTOP
    !
    ! Local variables
    REAL(EB) RN, simoDX, simoDY, TNOW
    REAL(EB) VOL1, VOL2, X1, X2, Y1, Y2, Z1, Z2, &
         dist, d_max, G_mean, G_sd, G_high, G_low, x1_old, y1_old
    INTEGER i,j,k,ii,jj,kk,ipc, izero, n_tmp, ie, nom, I_OBST
    INTEGER i11, i22, group_size
    LOGICAL pp_see_group, is_solid
    INTEGER jjj, iii, i44
    REAL(EB) x11, y11, group_x_sum, group_y_sum, group_x_center, group_y_center, dens_fac
    INTEGER :: i_endless_loop, istat
    REAL(EB), DIMENSION(6) :: y_tmp, x_tmp, r_tmp
    REAL(EB), DIMENSION(4) :: d_xy
    LOGICAL, DIMENSION(4) :: FoundWall_xy
    ! 
    TYPE (MESH_TYPE), POINTER :: M =>NULL()
    TYPE (EVAC_SSTAND_Type),POINTER :: ESS=>NULL()
    TYPE (EVACUATION_Type), POINTER :: HPT=>NULL()
    TYPE (EVAC_PERS_Type),  POINTER :: PCP=>NULL()
    TYPE (EVAC_HOLE_Type),  POINTER :: EHX=>NULL()
    TYPE (HUMAN_TYPE), POINTER :: HR=>NULL(), HRE=>NULL()
    !
    TNOW = SECOND()

    IF ( .NOT.(EVACUATION_ONLY(NM) .AND. EVACUATION_GRID(NM)) ) RETURN
    ! Next means that only EVAC_PROCESS is doing something
    IF (MYID /= PROCESS(NM)) RETURN
    !
    ! Gaussian random numbers, initialize (only once during
    ! the whole calculation is needed). We are now in the
    ! initialization routine, so this is called before the
    ! use of GaussFlag in the evacuation movement part.
    GTrunFlag = 0
    GaussFlag = 0
    !
    MESHES(NM)%N_HUMANS = 0
    MESHES(NM)%N_HUMANS_DIM = 10000
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
    !
    ! Initialise HUMAN_GRID
    !
    FED_I_LOOP: DO i = 1,IBAR
       FED_J_LOOP: DO j= 1,JBAR
          x1 = XS + (i-1)*DXI + 0.5_EB*DXI
          y1 = YS + (j-1)*DETA + 0.5_EB*DETA
          ! z1 is here the head height where smoke information is taken, so it
          ! should be about the head height
          z1 = 0.5_EB*(ZF+ZS) - EVACUATION_Z_OFFSET(NM) + HUMAN_SMOKE_HEIGHT
          SS_Loop: DO k = 1, n_sstands
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
             END IF
             EXIT SS_Loop
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
                IF ( X1 >= M%XS .AND. X1 <= M%XF .AND. &
                     Y1 >= M%YS .AND. Y1 <= M%YF .AND. &
                     Z1 >= M%ZS .AND. Z1 <= M%ZF) THEN
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
       HPT=>EVACUATION(IPC)
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
       ! 
       ! Check which evacuation floor node  (=1,...,n_egrids)
       n_tmp = 0
       HP_MeshLoop: DO i = 1, nmeshes
          IF (evacuation_only(i) .AND. evacuation_grid(i)) THEN
             n_tmp = n_tmp +1
             IF (HPT%IMESH == i) THEN
                EXIT HP_MeshLoop
             END IF
          END IF
       END DO HP_MeshLoop
       IF (n_tmp < 1 .OR. n_tmp > n_egrids) THEN
          WRITE(MESSAGE,'(A,A,A,I4)') 'ERROR: INIT_EVAC ',TRIM(HPT%ID),' problem evac node, INODE= ',n_tmp
          CALL SHUTDOWN(MESSAGE)
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

       INITIALIZATION_LOOP: DO I=1,HPT%N_INITIAL
          !
          CALL RANDOM_NUMBER(RN)
          group_size = HPT%GN_MIN - 1 +  INT((HPT%GN_MAX-HPT%GN_MIN+1)*RN+0.5_EB)
          group_size = MAX(MIN(group_size, HPT%GN_MAX),HPT%GN_MIN)

          IF ( i11+group_size > HPT%N_INITIAL ) THEN
             group_size = HPT%N_INITIAL - i11
          END IF

          IF ( i11 >= HPT%N_INITIAL ) EXIT INITIALIZATION_LOOP

          i22 = 0
          i_endless_loop = 0
          group_X_sum = 0
          group_Y_sum = 0
          GROUP_SIZE_LOOP: DO 
             ! i22: Counter, humans on this group (group index i33)
             i22 = i22 + 1
             IF (i22 > group_size) EXIT GROUP_SIZE_LOOP
             i11 = i11 + 1

             IF (i22 == 1) THEN
                ! One member groups are not counted as group
                IF (group_size > 1) i33 = i33 + 1
                IF (group_size == 1) ilh = ilh + 1
             END IF
             N_HUMANS = N_HUMANS + 1
             !
             IF (N_HUMANS > N_HUMANS_DIM) THEN
                CALL SHUTDOWN('ERROR: Init Humans: no re-allocation yet')
                CALL RE_ALLOCATE_HUMANS(1,NM)
                HUMAN=>MESHES(NM)%HUMAN
             END IF
             !
             HR => HUMAN(N_HUMANS)
             CALL CLASS_PROPERTIES(HR,PCP)
             HR%IPC = HPT%IPC  ! PERS-line index
             HR%IEL = IPC      ! EVAC-line index
             HR%I_Target = 0

             !
             BLK_LOOP:  DO
                IF (i22 == 1) THEN
                   ! First member of the group
                   CALL RANDOM_NUMBER(RN)
                   HR%X = X1 + RN*(X2-X1)
                   x1_old = HR%X
                   CALL RANDOM_NUMBER(RN)
                   HR%Y = Y1 + RN*(Y2-Y1)
                   y1_old = HR%Y
                ELSE
                   ! Next members of the group are put around the first member.
                   G_mean = 0.0_EB
                   G_sd   = 4.0_EB  ! units (m) std.dev.
                   G_high = 6.0_EB  ! units (m) cut-off
                   G_low  = -6.0_EB ! units (m) cut-off
                   ! First the angle, then the radial distance
                   CALL RANDOM_NUMBER(rn)
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
                HR%Z = Z1 + 0.5_EB*(Z2-Z1)
                IF (HPT%Angle > -999.9_EB) THEN
                   HR%Angle = HPT%Angle
                ELSE
                   CALL RANDOM_NUMBER(RN)
                   HR%Angle = 2.0_EB*Pi*rn
                END IF
                DO WHILE (HR%Angle >= 2.0_EB*Pi)
                   HR%Angle = HR%Angle - 2.0_EB*Pi
                END DO
                DO WHILE (HR%Angle < 0.0_EB)
                   HR%Angle = HR%Angle + 2.0_EB*Pi
                END DO

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

                !Check, that a person is not put on top of some other person
                IF (DENS_INIT > 2.0_EB) THEN
                   ! High density is wanted
                   d_max = 0.0_EB
                ELSE
                   d_max = 1.0_EB*HR%B
                END IF
                dens_fac = MAX(1.0_EB,DENS_INIT)

                IF (i_endless_loop >= 10*INT(dens_fac*(16.0_EB*MAX(1.0_EB,LOG10(2.5_EB*VOL2))) / &
                     MAX(1.0_EB,LOG10((2.5_EB*VOL2)/(2.5_EB*VOL2-1)))) ) THEN
                   WRITE (LU_EVACOUT,fmt='(A,I4,A,I4,A,I6)') ' ERROR: Initialize_Humans, EVAC line ', &
                        IPC, ', Mesh ', NM, ', i_human ', n_humans
                   WRITE (LU_EVACOUT,fmt='(a)') '      x       y       z     Rd      Rt      Rs      ds  '
                   WRITE (LU_EVACOUT,fmt='(3f8.2,4f8.4)') HR%X, HR%Y, HR%Z, &
                        2.0_EB*HR%Radius, HR%r_torso, HR%r_shoulder, HR%d_shoulder
                   ISTOP = 3  ! Stop code: FDS improperly set-up
                   HR%SHOW = .TRUE.    
                   HR%COLOR_INDEX = 7  ! Cyan
                   EXIT INITIALIZATION_LOOP
                END IF

                ! Put a human in a mesh, if there is enough empty space, i.e.,
                ! check the OBSTs.
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
                II = FLOOR( CELLSI(FLOOR((x_tmp(1)-XS)*RDXINT)) + 1.0_EB )
                JJ = FLOOR( CELLSJ(FLOOR((y_tmp(1)-YS)*RDYINT)) + 1.0_EB )
                I_OBST = OBST_INDEX_C(CELL_INDEX(II,JJ,KK))
                Is_Solid = (Is_Solid .OR. (SOLID(CELL_INDEX(II,JJ,KK)) .AND. .NOT. OBSTRUCTION(I_OBST)%HIDDEN))
                II = FLOOR( CELLSI(FLOOR((x_tmp(3)-XS)*RDXINT)) + 1.0_EB )
                JJ = FLOOR( CELLSJ(FLOOR((y_tmp(3)-YS)*RDYINT)) + 1.0_EB )
                Is_Solid = (Is_Solid .OR. (SOLID(CELL_INDEX(II,JJ,KK)) .AND. .NOT. OBSTRUCTION(I_OBST)%HIDDEN))
                II = FLOOR( CELLSI(FLOOR((x_tmp(2)-XS)*RDXINT)) + 1.0_EB )
                JJ = FLOOR( CELLSJ(FLOOR((y_tmp(2)-YS)*RDYINT)) + 1.0_EB )
                Is_Solid = (Is_Solid .OR. (SOLID(CELL_INDEX(II,JJ,KK)) .AND. .NOT. OBSTRUCTION(I_OBST)%HIDDEN))

                IF (.NOT.Is_Solid) THEN
                   ! Check the distances to walls (not to wall corners)
                   ! Skip wall force ior = 0: Check all walls (do not put too
                   ! close to doors/exits)
                   CALL Find_walls(nm, x_tmp(1), y_tmp(1), r_tmp(1), d_max, 0, d_xy, FoundWall_xy, istat)
                   IF (istat/=0) Is_Solid = .TRUE.
                   CALL Find_walls(nm, x_tmp(2), y_tmp(2), r_tmp(2), d_max, 0, d_xy, FoundWall_xy, istat)
                   IF (istat/=0) Is_Solid = .TRUE.
                   CALL Find_walls(nm, x_tmp(3), y_tmp(3), r_tmp(3), d_max, 0, d_xy, FoundWall_xy, istat)
                   IF (istat/=0) Is_Solid = .TRUE.
                END IF

                IF (.NOT.Is_Solid) THEN
                   ! Check that the agent is not too close to other agents, who
                   ! already introduced in the calculation.
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

                      PP_see_group = See_each_other(nm, x11, y11, x1_old, y1_old)
                      
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
                CALL SHUTDOWN(MESSAGE)
             END SELECT

             HR%IMESH       = HPT%IMESH
             HR%INODE       = n_tmp
             HR%NODE_NAME   = TRIM(EVAC_Node_List(n_tmp)%ID)
             HR%I_FFIELD  = 0
             HR%IOR       = HUMAN_NO_TARGET
             HR%U         = 0.0_EB
             HR%V         = 0.0_EB
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
             IF ( group_size == 1 ) EXIT  GROUP_SIZE_LOOP
             IF ( i22 == group_size .AND. group_size > 1 ) THEN
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
                   PP_see_group = See_each_other(nm, x11, y11, x1_old, y1_old)

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
    END DO EVAC_CLASS_LOOP ! ipc, number of evac-lines

    WRITE (LU_EVACOUT,fmt='(a,f8.2,a,i0,a,i0/)') ' EVAC: Time ', 0.0_EB,' mesh ',nm,' number of humans ',n_humans
    ! Write (LU_OUTPUT,fmt='(a,f8.2,a,i0,a,i0/)') ' EVAC: Time ', &
    !      0.0_EB,' mesh ',nm,' number of humans ',n_humans
    !
    TUSED(12,NM)=TUSED(12,NM)+SECOND()-TNOW
    !
  END SUBROUTINE INITIALIZE_EVACUATION

  !
  SUBROUTINE INIT_EVAC_GROUPS
    IMPLICIT NONE
    !
    ! Initialize group lists, known doors, etc
    !
    ! Local variables
    INTEGER I,J, IZERO, nom, j1, ii, i_target_old, i_change_old, i_tmp, i_tmp2
    INTEGER :: i_egrid, i_target, color_index, i_new_ffield
    REAL(EB) :: evel
    ! 
    TYPE (MESH_TYPE), POINTER :: M=>NULL()
    TYPE (EVAC_SSTAND_Type), POINTER :: ESS=>NULL()
    TYPE (HUMAN_TYPE), POINTER :: HR=>NULL()
    !
    IF (.NOT.ANY(EVACUATION_GRID)) RETURN
    IF (PROCESS_STOP_STATUS > 0) RETURN

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
    ALLOCATE(Is_Known_Door(MAX(1,n_doors+n_exits)),STAT=IZERO)
    CALL ChkMemErr('Initialize_Evacuation','Is_Known_Door',IZERO) 
    ALLOCATE(Is_Visible_Door(MAX(1,n_doors+n_exits)),STAT=IZERO)
    CALL ChkMemErr('Initialize_Evacuation','Is_Visible_Door',IZERO)
    ALLOCATE(FED_max_Door(MAX(1,n_doors+n_exits)),STAT=IZERO)
    CALL ChkMemErr('Initialize_Evacuation','FED_max_Door',IZERO) 
    ALLOCATE(K_ave_Door(MAX(1,n_doors+n_exits)),STAT=IZERO)
    CALL ChkMemErr('Initialize_Evacuation','K_ave_Door',IZERO) 
    ALLOCATE(Color_Tmp(MAX(1,i33_dim)),STAT=IZERO)
    CALL ChkMemErr('Initialize_Evacuation','Color_Tmp',IZERO) 

    i_egrid = 0
    IF (n_doors >0) EVAC_DOORS(:)%R_NTARGET = 5.0_EB
    IF (n_exits >0) EVAC_EXITS(:)%R_NTARGET = 5.0_EB
    DO nom = 1, NMESHES
       n_change_doors  = 0 ! Count the initialization Nash Equilibrium iterations
       n_change_trials = 0 ! Count the initialization Nash Equilibrium iterations
       i_change_old    = -1
       IF ( .NOT.(EVACUATION_ONLY(NoM) .AND. EVACUATION_GRID(NoM)) ) CYCLE
       M => MESHES(NOM)
       Group_List(:)%GROUP_SIZE  = 0
       Group_List(:)%GROUP_X = 0.0_EB
       Group_List(:)%GROUP_Y = 0.0_EB
       Group_List(:)%Speed   = 0.0_EB
       i_egrid = i_egrid+1
       DO i = 1, M%N_HUMANS
          HR => M%HUMAN(I)
          j = MAX(0,HR%GROUP_ID)
          Group_List(j)%GROUP_SIZE = Group_List(j)%GROUP_SIZE + 1
          Group_List(j)%GROUP_X    = Group_List(j)%GROUP_X + HR%X
          Group_List(j)%GROUP_Y    = Group_List(j)%GROUP_Y + HR%Y
          Group_List(j)%Speed      = Group_List(j)%Speed + HR%Speed
          DO ii = 1,n_doors
             IF (nom .NE. EVAC_DOORS(ii)%Imesh) CYCLE
             evel = SQRT((HR%X-EVAC_DOORS(ii)%X)**2 + (HR%Y-EVAC_DOORS(ii)%Y)**2)
             EVAC_DOORS(ii)%R_NTARGET = MAX(evel,EVAC_DOORS(ii)%R_NTARGET)
          END DO
          DO ii = 1,n_exits
             IF (nom .NE. EVAC_EXITS(ii)%Imesh) CYCLE
             evel = SQRT((HR%X-EVAC_EXITS(ii)%X)**2 + (HR%Y-EVAC_EXITS(ii)%Y)**2)
             EVAC_EXITS(ii)%R_NTARGET = MAX(evel,EVAC_EXITS(ii)%R_NTARGET)
          END DO
       END DO
       Group_List(1:)%GROUP_X = Group_List(1:)%GROUP_X / MAX(1,Group_List(1:)%GROUP_SIZE)
       Group_List(1:)%GROUP_Y = Group_List(1:)%GROUP_Y / MAX(1,Group_List(1:)%GROUP_SIZE)
       Group_List(1:)%Speed   = Group_List(1:)%Speed / MAX(1,Group_List(1:)%GROUP_SIZE)
       i_tmp = 0
       i_tmp2 = -1
       IF (M%N_HUMANS < 1) CYCLE
       DO WHILE (.NOT.(i_change_old == n_change_doors) .AND. .NOT.(i_tmp2 == i_tmp)) ! Iterate Nash equilibrium
          i_change_old = n_change_doors
          i_tmp = 0
          DO i = 1, M%N_HUMANS
             HR => M%HUMAN(I)
             ! group_id > 0: +group_id
             ! group_id < 0: -human_id (lonely humans)
             j  =  MAX(0,HR%GROUP_ID)
             j1 = -MIN(0,HR%GROUP_ID)
             ! Test, if this group has already a ffield (on this floor)
             ! Lonely humans have j=0 and group_i_ffields is always 0.
             IF (Group_List(j)%GROUP_I_FFIELDS(i_egrid) == 0) THEN
                i_target_old = HR%I_Target
                n_change_trials = n_change_trials + 1
                CALL Change_Target_Door(nom, nom, i, j, j1, i_egrid, 0, HR%X, HR%Y, i_target, color_index, i_new_ffield, HR)
                IF (ABS(i_target_old) .NE. ABS(i_target)) THEN
                   n_change_doors = n_change_doors + 1
                   i_tmp = i
                   IF (i_target > 0) THEN
                      IF (i_target <= n_doors ) THEN
                         evel = SQRT((HR%X-EVAC_DOORS(i_target)%X)**2 + (HR%Y-EVAC_DOORS(i_target)%Y)**2)
                         evel = 50.0_EB*evel/EVAC_DOORS(i_target)%R_NTARGET + 1.0_EB
                         ii = MIN(50,MAX(1,INT(evel)))
                         EVAC_DOORS(i_target)%NTARGET(ii:50) = EVAC_DOORS(i_target)%NTARGET(ii:50) + 1
                      ELSE
                         evel = SQRT((HR%X-EVAC_EXITS(i_target-n_doors)%X)**2 + (HR%Y-EVAC_EXITS(i_target-n_doors)%Y)**2)
                         evel = 50.0_EB*evel/EVAC_EXITS(i_target-n_doors)%R_NTARGET + 1.0_EB
                         ii = MIN(50,MAX(1,INT(evel)))
                         EVAC_EXITS(i_target-n_doors)%NTARGET(ii:50) = EVAC_EXITS(i_target-n_doors)%NTARGET(ii:50) + 1
                      END IF
                   ELSE  ! i_target < 0 means non-visible target door
                      n_change_doors = n_change_doors - 1  ! do not iterate non-visible doors
                   END IF
                   IF (i_target_old > 0) THEN
                      IF (i_target_old <= n_doors ) THEN
                         evel = SQRT((HR%X-EVAC_DOORS(i_target_old)%X)**2 + (HR%Y-EVAC_DOORS(i_target_old)%Y)**2)
                         evel = 50.0_EB*evel/EVAC_DOORS(i_target_old)%R_NTARGET + 1.0_EB
                         ii = MIN(50,MAX(1,INT(evel)))
                         EVAC_DOORS(i_target_old)%NTARGET(ii:50) = EVAC_DOORS(i_target_old)%NTARGET(ii:50) - 1
                      ELSE
                         evel = SQRT((HR%X-EVAC_EXITS(i_target_old-n_doors)%X)**2 + &
                              (HR%Y-EVAC_EXITS(i_target_old-n_doors)%Y)**2)
                         evel = 50.0_EB*evel/EVAC_EXITS(i_target_old-n_doors)%R_NTARGET + 1.0_EB
                         ii = MIN(50,MAX(1,INT(evel)))
                         EVAC_EXITS(i_target_old-n_doors)%NTARGET(ii:50) = &
                              EVAC_EXITS(i_target_old-n_doors)%NTARGET(ii:50) - 1
                      END IF
                   END IF
                END IF
             ELSE              ! this group has already a flow field
                ! This group has already tried to change the field
                IF (COLOR_METHOD == 5 .AND. j > 0) HR%COLOR_INDEX = Color_Tmp(j)
                IF (COLOR_METHOD == 4 .AND. j > 0) HR%COLOR_INDEX = Color_Tmp(j)
                HR%I_FFIELD    = Group_List(j)%GROUP_I_FFIELDS(i_egrid)
                HR%FFIELD_NAME = TRIM(MESH_NAME(HR%I_FFIELD))
                HR%I_Target = Group_Known_Doors(j)%I_Target
             END IF            ! first member of a group or a lonely soul
          END DO              ! 1, n_humans
          IF (n_change_doors-i_change_old == 1) THEN
             WRITE(LU_EVACOUT,fmt='(a,2i10)') ' Init: Door Changes i_tmp ',  i_tmp, i_tmp2
             i_tmp2 = i_tmp
          ELSE
             i_tmp2 = -1
          END IF
          IF (n_change_doors/MAX(1,M%N_HUMANS) > 10*M%N_HUMANS) i_tmp2 = i_tmp
          IF (FAC_DOOR_QUEUE <= 0.001_EB) i_change_old = n_change_doors  ! Do not iterate the Nash equilibrium
       END DO         ! Nash iterations
       IF (FAC_DOOR_QUEUE > 0.001_EB) WRITE(LU_EVACOUT,fmt='(a,f14.2,a,i8)') &
            ' Init: Changes per agent ', REAL(n_change_doors,EB)/REAL(M%N_HUMANS,EB), &
            ', Nash iterations', n_change_trials/M%N_HUMANS
    END DO                  ! 1, nmeshes

    WRITE (LU_EVACOUT,fmt='(/a)') ' EVAC: Initial positions of the humans'
    WRITE (LU_EVACOUT,fmt='(a,a)') ' person     x       y       z    Tpre    Tdet  ', &
         ' dia    v0   tau   i_gr i_ff'

    ! Initialize the group_i_ffields
    i_egrid = 0
    DO nom = 1, NMESHES
       IF ( .NOT.(EVACUATION_ONLY(NoM) .AND. EVACUATION_GRID(NoM)) ) CYCLE
       M => MESHES(NOM)
       i_egrid = i_egrid+1
       DO i = 1, M%N_HUMANS
          HR => M%HUMAN(I)

          ! Check spectator stands, correct the z-coordiante
          SS_Loop: DO j = 1, n_sstands
             ESS => EVAC_SSTANDS(j)
             IF (ESS%IMESH == nom .AND. (ESS%X1 <= HR%X .AND. ESS%X2 >= HR%X) .AND. &
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
                EXIT SS_Loop
             END IF
          END DO SS_Loop

          j = MAX(0,HR%GROUP_ID)
          Group_List(j)%IEL = HR%IEL

          WRITE (LU_EVACOUT,fmt='(i6,5f8.2,3f6.2,i6,i4,i4)') HR%ILABEL, &
               HR%X, HR%Y, HR%Z, HR%Tpre, HR%Tdet,2.0_EB*HR%Radius, &
               HR%Speed, HR%Tau, HR%GROUP_ID, HR%i_ffield, HR%COLOR_INDEX
       END DO
    END DO
    WRITE (LU_EVACOUT,fmt='(/)')
    !n_change_doors  = 0 ! Do not count the initialization Nash Equilibrium iterations
    !n_change_trials = 0 ! Do not count the initialization Nash Equilibrium iterations

  END SUBROUTINE INIT_EVAC_GROUPS
!
  SUBROUTINE EVAC_MESH_EXCHANGE(T,T_Save,I_mode, ICYC, EXCHANGE_EVACUATION, MODE)
    IMPLICIT NONE
    !
    ! Passed variables
    REAL(EB), INTENT(IN) :: T
    INTEGER, INTENT(IN) :: I_mode, ICYC, MODE
    LOGICAL, INTENT(OUT) :: EXCHANGE_EVACUATION
    REAL(EB), INTENT(INOUT) :: T_Save
    !
    ! Local variables
    INTEGER :: nm, nom, i, j, i1, j1, k1
    INTEGER :: ios, IZERO
    LOGICAL L_use_fed, L_fed_read, L_fed_save
    REAL(EB) DT_Save
    INTEGER(4) ibar_tmp, jbar_tmp, kbar_tmp, n_tmp
    REAL(FB) tmpout1, tmpout2, tmpout3, tmpout4, t_tmp, dt_tmp
    REAL(FB) tmpout5, tmpout6, tmpout7, tmpout8
    REAL(EB), ALLOCATABLE, DIMENSION(:) :: YY_GET
    !
    IF (.NOT. ANY(EVACUATION_GRID)) RETURN
    IF (ICYC < 1) RETURN
    !
    ! I_mode: 'binary' index:
    ! xxxxx = (0,1)*16 + (0,1)*8 + (0,1)*4 + (0,1)*2 + (0,1)
    ! 0. bit (xxxx1): save flow fields (xxxx0) do not save
    ! 1. bit (xxx1x): save soot + fed  (xxx0x) do not save
    ! 2. bit (xx1xx): read flow fields (xx0xx) do not read
    ! 3. bit (x1xxx): read soot + fed  (x0xxx) do not read
    ! 4. bit (1xxxx): fire calculation (0xxxx) no fire calculation
    ! BTEST(integer,bit) bit=0 least significant bit
    ! IBCLR(integer,bit) clears a bit
    ! IBSET(integer,bit) sets a bit
    !
    ! LU_EVACFED: CHID_evac.fed, FED and soot, time dependent, binary
    ! File format: 1. row: n_egrids >=0  (Old Format)
    !              1a. row: n < 0 (New Format)
    !              1b. row: n_egrids,4,n_corrs,8 (New Format)
    !                 2. row: T and DT
    !                    3. row: ibar,jbar,kbar, n_quantities
    !                       4. row: onwards data
    !                    goto 3. (meshes)
    !                       N. row: corr data (8 real numbers)
    !                       N+1. row: next corr data...
    !                 goto 2. (time points)
    !
    ! Update interval (seconds) fire ==> evac information
    DT_Save = 2.0_EB
    ios = 0

    L_use_fed  = .FALSE.
    L_fed_read = BTEST(I_mode,3)
    L_fed_save = BTEST(I_mode,1)

    !
    ! Change information: Fire meshes ==> Evac meshes
    IF ( T >= T_BEGIN .AND. REAL(T,FB) >= REAL(T_Save,FB) ) THEN
       IF (MODE > 0) T_Save = T + DT_Save
       L_use_fed = .TRUE.
    END IF

    L_use_fed = L_use_fed .AND. (L_fed_read .OR. L_fed_save)

    IF (MODE < 1 .AND. L_use_fed) EXCHANGE_EVACUATION = .TRUE.
    IF (MODE < 2) RETURN

    IF (L_use_fed) THEN
       IF (L_fed_save) THEN
          ALLOCATE(YY_GET(1:MAX(1,N_SPECIES)),STAT=IZERO)
          CALL ChkMemErr('EVAC_MESH_EXCHANGE', 'YY_GET',IZERO) 
          WRITE (LU_EVACFED) REAL(T,FB), REAL(DT_Save,FB)
       ELSE
          READ (LU_EVACFED,END=324,Iostat=ios) t_tmp, dt_tmp
          T_Save = t_tmp + dt_tmp ! next time point in the file
       END IF

       ! Next loop interpolates fire mesh (soot+fed) into human_grids and
       ! saves it to the disk. Or it reads FED+soot from the disk.
       MESH_LOOP: DO NM=1,NMESHES
          IF ( .NOT.(EVACUATION_GRID(NM) .AND. EVACUATION_ONLY(NM)) ) CYCLE
          !
          CALL POINT_TO_MESH(NM)
          !             Call POINT_TO_EVAC_MESH(NM)
          IF (L_fed_save) THEN
             ibar_tmp = IBAR
             jbar_tmp = JBAR
             kbar_tmp = 1
             n_tmp    = 4  ! New format (version 1.11)
             WRITE (LU_EVACFED) ibar_tmp, jbar_tmp, kbar_tmp, n_tmp
          ELSE
             READ (LU_EVACFED,Iostat=ios) ibar_tmp, jbar_tmp, kbar_tmp, n_tmp
             IF (ios.NE.0) THEN
                WRITE(MESSAGE,'(A)') 'ERROR: Evac Mesh Exchange: FED READ ERROR'
                CLOSE (LU_EVACFED)
                CALL SHUTDOWN(MESSAGE)
             END IF
             IF (ibar_tmp /= IBAR .OR. jbar_tmp /= JBAR .OR. n_tmp < 4 ) THEN
                CLOSE (LU_EVACFED)
                CALL SHUTDOWN('ERROR: Problems to read the FED file')
             END IF

          END IF

          DO i = 1, IBAR
             DO j= 1, JBAR

                IF (L_fed_save) THEN
                   IF ( ABS(HUMAN_GRID(i,j)%IMESH) > 0 ) THEN
                      ! imesh > 0, i.e. fire grid found
                      i1 = HUMAN_GRID(i,j)%II
                      j1 = HUMAN_GRID(i,j)%JJ
                      k1 = HUMAN_GRID(i,j)%KK
                      nom = ABS(HUMAN_GRID(i,j)%IMESH)
                      CALL GET_FIRE_CONDITIONS(nom,i1,j1,k1,&
                           HUMAN_GRID(i,j)%FED_CO_CO2_O2,HUMAN_GRID(i,j)%SOOT_DENS,&
                           HUMAN_GRID(i,j)%TMP_G, HUMAN_GRID(i,j)%RADFLUX, YY_GET)
                   END IF
                   ! Save Fed, Soot, Temp(C), and Radflux
                   WRITE (LU_EVACFED) &
                        REAL(HUMAN_GRID(i,j)%FED_CO_CO2_O2,FB), &
                        REAL(HUMAN_GRID(i,j)%SOOT_DENS,FB), &
                        REAL(HUMAN_GRID(i,j)%TMP_G,FB), &
                        REAL(HUMAN_GRID(i,j)%RADFLUX,FB)
                ELSE ! Read FED from a file
                   ! Read Fed, Soot, Temp(C), and Radflux
                   READ (LU_EVACFED,Iostat=ios) tmpout1, tmpout2, tmpout3, tmpout4
                   IF (ios.NE.0) THEN
                      WRITE(MESSAGE,'(A)') 'ERROR: Evac Mesh Exchange: FED READ ERROR'
                      CLOSE (LU_EVACFED)
                      CALL SHUTDOWN(MESSAGE)
                   END IF
                   HUMAN_GRID(i,j)%FED_CO_CO2_O2 = tmpout1
                   HUMAN_GRID(i,j)%SOOT_DENS = tmpout2
                   HUMAN_GRID(i,j)%TMP_G = tmpout3
                   HUMAN_GRID(i,j)%RADFLUX = tmpout4
                END IF   ! Calculate and save FED

             END DO     ! j=1,JBAR
          END DO       ! i=1,IBAR

       END DO MESH_LOOP

       ! Next loop interpolates fire mesh (soot+fed) into human_grids and
       ! saves it to the disk. Or it reads FED+soot from the disk.
       CORR_LOOP: DO i = 1, N_CORRS
          !
          IF (L_fed_save) THEN

             IF ( EVAC_CORRS(i)%FED_MESH > 0 ) THEN
                ! Here the fire properties are saved to the arrays.
                i1 = EVAC_CORRS(i)%II(1)
                j1 = EVAC_CORRS(i)%JJ(1)
                k1 = EVAC_CORRS(i)%KK(1)
                nom = EVAC_CORRS(i)%FED_MESH
                CALL GET_FIRE_CONDITIONS(nom,i1,j1,k1,&
                         EVAC_CORRS(i)%FED_CO_CO2_O2(1),EVAC_CORRS(i)%SOOT_DENS(1),&
                         EVAC_CORRS(i)%TMP_G(1), EVAC_CORRS(i)%RADFLUX(1), YY_GET)
             ELSE
                ! No fed_mesh found
                EVAC_CORRS(i)%FED_CO_CO2_O2(1) = 0.0_EB
                EVAC_CORRS(i)%SOOT_DENS(1) = 0.0_EB
                EVAC_CORRS(i)%TMP_G(1) = 0.0_EB
                EVAC_CORRS(i)%RADFLUX(1) = 0.0_EB
             END IF                ! fed_mesh > 0, i.e. fire grid found

             IF ( EVAC_CORRS(i)%FED_MESH2 > 0 ) THEN
                i1 = EVAC_CORRS(i)%II(2)
                j1 = EVAC_CORRS(i)%JJ(2)
                k1 = EVAC_CORRS(i)%KK(2)
                nom = EVAC_CORRS(i)%FED_MESH2
                CALL GET_FIRE_CONDITIONS(nom,i1,j1,k1,&
                         EVAC_CORRS(i)%FED_CO_CO2_O2(2),EVAC_CORRS(i)%SOOT_DENS(2),&
                         EVAC_CORRS(i)%TMP_G(2), EVAC_CORRS(i)%RADFLUX(2), YY_GET)
             ELSE
                ! No fed_mesh2 found
                EVAC_CORRS(i)%FED_CO_CO2_O2(2) = 0.0_EB
                EVAC_CORRS(i)%SOOT_DENS(2) = 0.0_EB
                EVAC_CORRS(i)%TMP_G(2) = 0.0_EB
                EVAC_CORRS(i)%RADFLUX(2) = 0.0_EB
             END IF                ! fed_mesh2 > 0, i.e. fire grid found

             ! Save Fed, Soot, Temp(C), and Radflux
             WRITE (LU_EVACFED) &
                  REAL(EVAC_CORRS(i)%FED_CO_CO2_O2(1),FB), &
                  REAL(EVAC_CORRS(i)%SOOT_DENS(1),FB), &
                  REAL(EVAC_CORRS(i)%TMP_G(1),FB), &
                  REAL(EVAC_CORRS(i)%RADFLUX(1),FB), &
                  REAL(EVAC_CORRS(i)%FED_CO_CO2_O2(2),FB), &
                  REAL(EVAC_CORRS(i)%SOOT_DENS(2),FB), &
                  REAL(EVAC_CORRS(i)%TMP_G(2),FB), &
                  REAL(EVAC_CORRS(i)%RADFLUX(2),FB)
          ELSE                    ! Read FED from a file
             ! Read Fed, Soot, Temp(C), and Radflux
             READ (LU_EVACFED,Iostat=ios) tmpout1, tmpout2, tmpout3, tmpout4, tmpout5, tmpout6, tmpout7, tmpout8
             IF (ios.NE.0) THEN
                WRITE(MESSAGE,'(A)') 'ERROR: Evac Mesh Exchange: FED READ ERROR'
                CLOSE (LU_EVACEFF)
                CALL SHUTDOWN(MESSAGE)
             END IF
             EVAC_CORRS(i)%FED_CO_CO2_O2(1) = tmpout1
             EVAC_CORRS(i)%SOOT_DENS(1) = tmpout2
             EVAC_CORRS(i)%TMP_G(1) = tmpout3
             EVAC_CORRS(i)%RADFLUX(1) = tmpout4
             EVAC_CORRS(i)%FED_CO_CO2_O2(2) = tmpout5
             EVAC_CORRS(i)%SOOT_DENS(2) = tmpout6
             EVAC_CORRS(i)%TMP_G(2) = tmpout7
             EVAC_CORRS(i)%RADFLUX(2) = tmpout8

          END IF                  ! Calculate and save FED
       END DO CORR_LOOP

       DOOR_LOOP: DO i = 1, N_DOORS
          !
          IF (L_fed_save) THEN
             IF ( EVAC_DOORS(i)%FED_MESH > 0 ) THEN
                i1 = EVAC_DOORS(i)%II
                j1 = EVAC_DOORS(i)%JJ
                k1 = EVAC_DOORS(i)%KK
                nom = EVAC_DOORS(i)%FED_MESH
                CALL GET_FIRE_CONDITIONS(nom,i1,j1,k1,&
                              EVAC_DOORS(i)%FED_CO_CO2_O2,EVAC_DOORS(i)%SOOT_DENS,&
                              EVAC_DOORS(i)%TMP_G, EVAC_DOORS(i)%RADFLUX, YY_GET)
             ELSE
                ! No fed_mesh found
                EVAC_DOORS(i)%FED_CO_CO2_O2 = 0.0_EB
                EVAC_DOORS(i)%SOOT_DENS = 0.0_EB
                EVAC_DOORS(i)%TMP_G = 0.0_EB
                EVAC_DOORS(i)%RADFLUX = 0.0_EB
             END IF                ! fed_mesh > 0, i.e. fire grid found

             ! Save Fed, Soot, Temp(C), and Radflux
             WRITE (LU_EVACFED) &
                  REAL(EVAC_DOORS(i)%FED_CO_CO2_O2,FB), &
                  REAL(EVAC_DOORS(i)%SOOT_DENS,FB), &
                  REAL(EVAC_DOORS(i)%TMP_G,FB), &
                  REAL(EVAC_DOORS(i)%RADFLUX,FB)
          ELSE                    ! Read FED from a file
             ! Read Fed, Soot, Temp(C), and Radflux
             READ (LU_EVACFED,Iostat=ios) tmpout1, tmpout2, tmpout3, tmpout4
             IF (ios.NE.0) THEN
                WRITE(MESSAGE,'(A)') 'ERROR: Evac Mesh Exchange: FED READ ERROR'
                CLOSE (LU_EVACFED)
                CALL SHUTDOWN(MESSAGE)
             END IF
             EVAC_DOORS(i)%FED_CO_CO2_O2 = tmpout1
             EVAC_DOORS(i)%SOOT_DENS = tmpout2
             EVAC_DOORS(i)%TMP_G = tmpout3
             EVAC_DOORS(i)%RADFLUX = tmpout4
          END IF                  ! Calculate and save FED
       END DO DOOR_LOOP

       EXIT_LOOP: DO i = 1, N_EXITS
          !
          ! Do not save/read data for counters.
          IF (EVAC_EXITS(i)%COUNT_ONLY) CYCLE EXIT_LOOP
          IF (L_fed_save) THEN
             IF ( EVAC_EXITS(i)%FED_MESH > 0 ) THEN
                i1 = EVAC_EXITS(i)%II
                j1 = EVAC_EXITS(i)%JJ
                k1 = EVAC_EXITS(i)%KK
                nom = EVAC_EXITS(i)%FED_MESH
                CALL GET_FIRE_CONDITIONS(nom,i1,j1,k1,&
                         EVAC_EXITS(i)%FED_CO_CO2_O2,EVAC_EXITS(i)%SOOT_DENS,&
                         EVAC_EXITS(i)%TMP_G, EVAC_EXITS(i)%RADFLUX, YY_GET)
             ELSE
                ! No fed_mesh found
                EVAC_EXITS(i)%FED_CO_CO2_O2 = 0.0_EB
                EVAC_EXITS(i)%SOOT_DENS = 0.0_EB
                EVAC_EXITS(i)%TMP_G = 0.0_EB
                EVAC_EXITS(i)%RADFLUX = 0.0_EB
             END IF                ! fed_mesh > 0, i.e. fire grid found

             ! Save Fed, Soot, Temp(C), and Radflux
             WRITE (LU_EVACFED) &
                  REAL(EVAC_EXITS(i)%FED_CO_CO2_O2,FB), &
                  REAL(EVAC_EXITS(i)%SOOT_DENS,FB), &
                  REAL(EVAC_EXITS(i)%TMP_G,FB), &
                  REAL(EVAC_EXITS(i)%RADFLUX,FB)
          ELSE                    ! Read FED from a file
             ! Read Fed, Soot, Temp(C), and Radflux
             READ (LU_EVACFED,Iostat=ios) tmpout1, tmpout2, tmpout3, tmpout4
             IF (ios.NE.0) THEN
                WRITE(MESSAGE,'(A)') 'ERROR: Evac Mesh Exchange: FED READ ERROR'
                CLOSE (LU_EVACFED)
                CALL SHUTDOWN(MESSAGE)
             END IF
             EVAC_EXITS(i)%FED_CO_CO2_O2 = tmpout1
             EVAC_EXITS(i)%SOOT_DENS = tmpout2
             EVAC_EXITS(i)%TMP_G = tmpout3
             EVAC_EXITS(i)%RADFLUX = tmpout4
          END IF                  ! Calculate and save FED
       END DO EXIT_LOOP

       IF (L_fed_save) DEALLOCATE(YY_GET)

    END IF                    ! l_use_fed

324 CONTINUE
    IF (ios < 0) THEN 
       WRITE (LU_EVACOUT,fmt='(a,f12.4,a)') 'FED file EOF: time ', T_Save-DT_Save, ' not found'
       WRITE (LU_EVACOUT,fmt='(a)') 'FED file EOF: use previous values'
       T_Save = 1.0E15
    END IF

  END SUBROUTINE EVAC_MESH_EXCHANGE
!
  SUBROUTINE PREPARE_TO_EVACUATE(ICYC)
    IMPLICIT NONE
    !
    ! Do the mesh independent initializations for the 
    ! subroutine evacuate_humans.
    !
    ! Passed variables
    INTEGER, INTENT(IN) :: ICYC
    !
    ! Local variables
    LOGICAL L_eff_read, L_eff_save
    INTEGER(4) ibar_tmp, jbar_tmp, kbar_tmp
    INTEGER nm_tim, i, j
    ! 
    TYPE (MESH_TYPE), POINTER :: MFF=>NULL()

    IF (.NOT.ANY(EVACUATION_GRID)) RETURN

    L_eff_read = BTEST(I_EVAC,2)
    L_eff_save = BTEST(I_EVAC,0)
    IF ( ICYC == 0 .AND. L_eff_save ) THEN
       DO nm_tim = 1, NMESHES
          IF (EVACUATION_ONLY(nm_tim)) THEN
             MFF=>MESHES(nm_tim)
             ibar_tmp = MFF%IBAR
             jbar_tmp = MFF%JBAR
             kbar_tmp = 1
             WRITE (LU_EVACEFF) ibar_tmp, jbar_tmp, kbar_tmp
             DO  i = 0, MFF%IBAR+1
                DO j= 0, MFF%JBAR+1
                   WRITE (LU_EVACEFF) REAL(MFF%U(i,j,1),FB), REAL(MFF%V(i,j,1),FB)
                END DO
             END DO
          END IF
       END DO
       ! Clear the save bit, save is done only once.
       L_eff_save = .FALSE.
       I_EVAC = IBCLR(I_EVAC,0)
    END IF
    !
    ! Initialize counters only once for each time step.
    IF ( ICYC >= 0 .AND. icyc_old < ICYC ) THEN
       icyc_old = ICYC
       fed_max_alive = 0.0_EB
       fed_max       = 0.0_EB
    END IF

  END SUBROUTINE PREPARE_TO_EVACUATE
!
  SUBROUTINE EVACUATE_HUMANS(Tin,NM,ICYC)
    IMPLICIT NONE
    !
    ! Calculates the forces on humans and moves them.
    ! Uses a modified velocity-Verlet algorithm.
    ! The modification is that the motive force part of
    ! the force is done using the dissipative self-consistent
    ! velocity-Verlet algorithm by Vattulainen et al.
    !
    ! Inputs:
    !   Tin: End time, this routine makes the move Tin-DT ==> Tin (DT=FDS dt)
    !   NM: Mesh index
    !   ICYC: index of the FDS fire time step
    !
    ! Passed variables
    REAL(EB), INTENT(IN) :: Tin
    INTEGER, INTENT(IN) :: NM,ICYC
    !
    ! Local variables
    INTEGER, PARAMETER :: n_sectors = 2
    REAL(EB) DTSP,UBAR,VBAR,X1,Y1,XI,YJ,ZK
    INTEGER ICN,I,J,IIN,JJN,KKN,II,JJ,KK,IIX,JJY,KKZ,ICX, ICY, N, J1, I_OBST, I_OBSTX, I_OBSTY
    INTEGER  IE, tim_ic, tim_iw, NM_now, tim_iwx, tim_iwy, tim_iw2, tim_ic2, ibc
    REAL(EB) :: P2P_DIST, P2P_DIST_MAX, P2P_U, P2P_V, EVEL, tim_dist
    REAL(EB), DIMENSION(4) :: d_xy
    LOGICAL, DIMENSION(4) :: FoundWall_xy
    INTEGER :: istat, STRS_INDX, i_target, color_index, i_new_ffield, i_target_old
    !
    !
    REAL(EB) ::  scal_prod_over_rsqr, U_new, V_new, Vmax_timo, CosPhiFac, &
         Speed_max, Delta_min, Dt_sum, C_Yeff, LambdaW, B_Wall, A_Wall, &
         T, Contact_F, Social_F, smoke_beta, smoke_alpha, smoke_speed_fac
    INTEGER :: iie, jje, iio, jjo, iii, jjj, i_egrid, i_tmp
    REAL(EB) :: x_now, y_now, d_humans, d_walls, DTSP_new, &
         fac_tim, dt_group_door, x11, y11, Speed, Tpre
    LOGICAL PP_see_each
    LOGICAL L_eff_read, L_eff_save, L_Dead
    REAL(EB) :: cos_x, cos_y, speed_xm, speed_xp, speed_ym, speed_yp, hr_z, hr_a, hr_b, hr_tau, hr_tau_iner
    !
    REAL(EB) :: rn, rnCF
    REAL(EB) :: GaMe, GaTh, GaCM
    !
    INTEGER :: i_old_ffield, IZERO, imode_old
    CHARACTER(26) :: name_old_ffield
    !
    REAL(EB) :: P2P_Torque, Fc_x, Fc_y, Omega_new, angle, A1, Tc_z, Fc_x1, Fc_y1
    REAL(EB) :: Omega_max, Omega_0, FAC_V0_UP, FAC_V0_DOWN, FAC_V0_HORI
    REAL(EB), DIMENSION(6) :: y_tmp, x_tmp, r_tmp, v_tmp, u_tmp
    !
    ! Next are for the block list of the double agent-agent loop (speed up)
    INTEGER :: Max_Humans_Cell, i_dx, j_dy, ie_max, bl_max
    REAL(EB) :: dx_min, dy_min
    INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: BLOCK_GRID
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: BLOCK_GRID_N
    INTEGER, DIMENSION(:), ALLOCATABLE :: BLOCK_LIST
    !
    REAL(EB) :: d_humans_min, d_walls_min
    REAL(EB) :: TNOW, tnow13, tnow14, tnow15
    !
    LOGICAL :: NM_STRS_MESH
    REAL(EB), DIMENSION(n_sectors+1) :: Sum_suunta, u_theta, v_theta, cos_theta, sin_theta, thetas
    REAL(EB) :: Sum_suunta_max, theta_step, theta_start, vr_2r, angle_hre, &
         angle_hr, v_hr, v_hre, P2P_Suunta_MAX, angle_old, hr_x, hr_y, d_new, d_shift, commitment
    INTEGER :: i_suunta_max, N_suunta_back, N_suunta_backCF
    INTEGER, DIMENSION(n_sectors+1) :: N_suunta, N_suuntaCF

    TYPE (MESH_TYPE), POINTER :: MFF=>NULL()
    TYPE (EVAC_STRS_Type), POINTER :: STRP=>NULL()
    TYPE (EVAC_SSTAND_Type), POINTER :: ESS=>NULL()
    TYPE (HUMAN_TYPE), POINTER :: HR=>NULL(), HRE=>NULL()
    !
    TNOW=SECOND()
    IF ( .NOT.(EVACUATION_ONLY(NM) .AND. EVACUATION_GRID(NM)) ) RETURN
    !
    ! Maximun speed of the agents.  
    Vmax_timo      = V_MAX
    ! Maximum angular speed of the agents.
    Omega_max      = V_ANGULAR_MAX*2.0_EB*Pi ! 8 rounds per second in radians
    ! The target angular speed of the agents (like v0 for translational motion)
    Omega_0        = V_ANGULAR*2.0_EB*Pi     ! 2 rounds per second in radians
    !
    Dt_sum         = 0.0_EB
    ! dt_group_door: how often (on the average) an agent (or group) tries to change the door.
    dt_group_door  = TAU_CHANGE_DOOR
    !
    CALL POINT_TO_MESH(NM)
    !
    ! Find the smallest delta_x and/or delta_y.  Note: evac meshes should not be stretched ones.
    dx_min = MINVAL(DX)
    dy_min = MINVAL(DY)
    Delta_min = MIN(dy_min, dx_min)
    !
    ! Read/Write the EFF-file?
    L_eff_read = BTEST(I_EVAC,2)
    L_eff_save = BTEST(I_EVAC,0)
    !
    ! Find the egrid index of this main evac mesh. [1, # main evac meshes]
    i_egrid = 0
    DO i = 1, n_egrids
       IF (EVAC_Node_List(i)%IMESH == NM) THEN
          i_egrid = EVAC_Node_List(i)%Node_Index
       END IF
    END DO
    IF (i_egrid == 0) THEN
       WRITE(MESSAGE,'(A,I6)') 'ERROR: Evacuate_Humans, no mesh found ',NM
       CALL SHUTDOWN(MESSAGE)
    END IF
    !
    ! Blocks are used to speed the double loops over the agent-agent interactions.
    ALLOCATE(BLOCK_GRID_N(1:IBAR,1:JBAR),STAT=IZERO)
    CALL ChkMemErr('EVACUATE_HUMANS','BLOCK_GRID_N',IZERO)
    ! Initialize some counters etc. for this main evac mesh.
    DO i = 1, n_doors
       IF (EVAC_DOORS(i)%IMESH == nm) THEN
          IF (EVAC_DOORS(i)%TIME_OPEN>EVAC_DOORS(i)%TIME_CLOSE) THEN
             imode_old=EVAC_DOORS(i)%IMODE
             IF (Tin>EVAC_DOORS(i)%TIME_CLOSE .AND. Tin<EVAC_DOORS(i)%TIME_OPEN .AND. imode_old==-1) EVAC_DOORS(i)%IMODE=+2
             IF (Tin>=EVAC_DOORS(i)%TIME_OPEN .AND. imode_old==-2) EVAC_DOORS(i)%IMODE=+1
          ELSE
             imode_old=EVAC_DOORS(i)%IMODE
             IF (Tin>EVAC_DOORS(i)%TIME_OPEN .AND. Tin<EVAC_DOORS(i)%TIME_CLOSE .AND. imode_old==-2) EVAC_DOORS(i)%IMODE=+1
             IF (Tin>=EVAC_DOORS(i)%TIME_CLOSE .AND. imode_old==-1) EVAC_DOORS(i)%IMODE=+2
          END IF
       END IF
    END DO
    DO i = 1, n_exits
       IF (EVAC_EXITS(i)%IMESH == nm .AND. .NOT.EVAC_EXITS(i)%COUNT_ONLY) THEN
          IF (EVAC_EXITS(i)%TIME_OPEN>EVAC_EXITS(i)%TIME_CLOSE) THEN
             imode_old=EVAC_EXITS(i)%IMODE
             IF (Tin>EVAC_EXITS(i)%TIME_CLOSE .AND. Tin<EVAC_EXITS(i)%TIME_OPEN .AND. imode_old==-1) EVAC_EXITS(i)%IMODE=+2
             IF (Tin>=EVAC_EXITS(i)%TIME_OPEN .AND. imode_old==-2) EVAC_EXITS(i)%IMODE=+1
          ELSE
             imode_old=EVAC_EXITS(i)%IMODE
             IF (Tin>EVAC_EXITS(i)%TIME_OPEN .AND. Tin<EVAC_EXITS(i)%TIME_CLOSE .AND. imode_old==-2) EVAC_EXITS(i)%IMODE=+1
             IF (Tin>=EVAC_EXITS(i)%TIME_CLOSE .AND. imode_old==-1) EVAC_EXITS(i)%IMODE=+2
          END IF
       END IF
    END DO
    ! 
    HUMAN_TIME_LOOP: DO WHILE ( Dt_sum < DT )
       ! DT is the fds flow calculation time step.
       ! Sometimes agent time step is smaller than the fire time step, so
       ! syncronize the agent clock correctly.
       ! Tsteps: Time steps of the agent movement algorithm for different meshes,
       ! which were calculated at the end of the previous main time step.
       ! DTSP: The present time step for the agent movement algorithm, which may be
       ! smaller than the FDS main time step DT.
       DTSP = MIN( (DT-Dt_sum), Tsteps(nm) )

       Dt_sum = Dt_sum + DTSP
       T = Tin - DT + Dt_sum     ! Current time for the agents

       DO i = 1, n_doors
          IF (EVAC_DOORS(i)%IMESH == nm) THEN
             EVAC_DOORS(i)%NTARGET(1:50) = 0
          END IF
       END DO
       DO i = 1, n_exits
          IF (EVAC_EXITS(i)%IMESH == nm) THEN
             EVAC_EXITS(i)%NTARGET(1:50) = 0
          END IF
       END DO

       ! ================================================
       ! Initialize group arrays for this main evac mesh.
       ! ================================================
       Group_List(:)%GROUP_SIZE      = 0
       Group_List(:)%GROUP_X         = 0.0_EB
       Group_List(:)%GROUP_Y         = 0.0_EB
       Group_List(:)%MAX_DIST_CENTER = 0.0_EB
       Group_List(:)%Speed           = 0.0_EB
       Group_List(:)%IntDose         = 0.0_EB
       Group_List(:)%Tpre            = 0.0_EB
       Group_List(:)%Tdet            = HUGE(Group_List(:)%Tdet)
       !
       DO j = 0, i33_dim
          Group_List(j)%GROUP_I_FFIELDS(i_egrid) = 0
       END DO
       IF (n_doors >0) EVAC_DOORS(:)%R_NTARGET = 5.0_EB
       IF (n_exits >0) EVAC_EXITS(:)%R_NTARGET = 5.0_EB
       DO i = 1, N_HUMANS
          HR=>HUMAN(I)
          j = MAX(0,HR%GROUP_ID)  ! Group index of the agent
          Group_List(j)%GROUP_SIZE = Group_List(j)%GROUP_SIZE + 1
          Group_List(j)%GROUP_X    = Group_List(j)%GROUP_X + HR%X
          Group_List(j)%GROUP_Y    = Group_List(j)%GROUP_Y + HR%Y
          Group_List(j)%Speed      = Group_List(j)%Speed + HR%Speed
          Group_List(j)%IntDose    = Group_List(j)%IntDose + HR%IntDose
          Group_List(j)%Tpre       = MAX(Group_List(j)%Tpre,HR%Tpre)
          Group_List(j)%Tdet       = MIN(Group_List(j)%Tdet,HR%Tdet)
          i_tmp = HR%I_TARGET
          IF (i_tmp > 0 .AND. i_tmp <= n_doors ) THEN
             evel = SQRT((HR%X- EVAC_DOORS(i_tmp)%X)**2 + (HR%Y-EVAC_DOORS(i_tmp)%Y)**2)
             EVAC_DOORS(i_tmp)%R_NTARGET = MAX(evel,EVAC_DOORS(i_tmp)%R_NTARGET)
          END IF
          i_tmp = i_tmp - n_doors
          IF (i_tmp > 0 .AND. i_tmp <= n_exits ) THEN
             evel = SQRT((HR%X-EVAC_EXITS(i_tmp)%X)**2 + (HR%Y-EVAC_EXITS(i_tmp)%Y)**2)
             EVAC_EXITS(i_tmp)%R_NTARGET = MAX(evel,EVAC_EXITS(i_tmp)%R_NTARGET)
          END IF
       END DO

       Group_List(1:)%GROUP_X = Group_List(1:)%GROUP_X / MAX(1,Group_List(1:)%GROUP_SIZE)
       Group_List(1:)%GROUP_Y = Group_List(1:)%GROUP_Y / MAX(1,Group_List(1:)%GROUP_SIZE)
       Group_List(1:)%Speed   = Group_List(1:)%Speed   / MAX(1,Group_List(1:)%GROUP_SIZE)
       Group_List(1:)%IntDose = Group_List(1:)%IntDose / MAX(1,Group_List(1:)%GROUP_SIZE)
       Group_List(:)%MAX_DIST_CENTER = 0.0_EB

       ! Group_List intitialization and count the number of agents going towards different doors/exits.
       DO i = 1, N_HUMANS
          HR=>HUMAN(I)
          ! group_id > 0: +group_id
          ! group_id < 0: -human_id (lonely agents)
          j  =  MAX(0,HR%GROUP_ID)  ! group index > 0
          j1 = -MIN(0,HR%GROUP_ID)  ! lonely agent index > 0
          Group_List(j)%MAX_DIST_CENTER = MAX(Group_List(j)%MAX_DIST_CENTER, &
               SQRT((HR%X - Group_List(j)%GROUP_X)**2 + (HR%Y - Group_List(j)%GROUP_Y)**2))
          i_tmp = HR%I_TARGET
          IF (i_tmp > 0 .AND. i_tmp <= n_doors ) THEN
             evel = SQRT((HR%X- EVAC_DOORS(i_tmp)%X)**2 + (HR%Y-EVAC_DOORS(i_tmp)%Y)**2)
             evel = 50.0_EB*evel/EVAC_DOORS(i_tmp)%R_NTARGET + 1.0_EB
             ie = MIN(50,MAX(1,INT(evel)))
             EVAC_DOORS(i_tmp)%NTARGET(ie:50) = EVAC_DOORS(i_tmp)%NTARGET(ie:50) + 1
          END IF
          i_tmp = i_tmp - n_doors
          IF (i_tmp > 0 .AND. i_tmp <= n_exits ) THEN
             evel = SQRT((HR%X-EVAC_EXITS(i_tmp)%X)**2 + (HR%Y-EVAC_EXITS(i_tmp)%Y)**2)
             evel = 50.0_EB*evel/EVAC_EXITS(i_tmp)%R_NTARGET + 1.0_EB
             ie = MIN(50,MAX(1,INT(evel)))
             EVAC_EXITS(i_tmp)%NTARGET(ie:50) = EVAC_EXITS(i_tmp)%NTARGET(ie:50) + 1
          END IF
       END DO

       ! j=0, i.e., lonely humans: Group_List is not used , but it is safe to initialize.
       IF (N_HUMANS > 0) THEN
          Group_List(0)%GROUP_SIZE = 1
          Group_List(0)%COMPLETE   = 1
          Group_List(0)%GROUP_X    = 0.5_EB*(XS+XF)
          Group_List(0)%GROUP_Y    = 0.5_EB*(YS+YF)
          Group_List(0)%Speed      = 1.0_EB
          Group_List(0)%IntDose    = 0.0_EB
          Group_List(0)%MAX_DIST_CENTER = 0.0_EB
       END IF

       ! Note: group_size is the number of group members on this main evac mesh.
       ! Check if the groups are already together or not.
       DO j = 1, i33_dim
          Group_List(j)%LIMIT_COMP = RADIUS_COMPLETE_0 + RADIUS_COMPLETE_1*Group_List(j)%GROUP_SIZE
          IF ( ((Group_List(j)%MAX_DIST_CENTER <=  Group_List(j)%LIMIT_COMP) .OR. &
               (Group_List(j)%COMPLETE == 1)) .AND. Group_List(j)%GROUP_SIZE > 0) THEN
             ! If complete=1 already, it stays at 1.
             ! If many floors, check only floors where there are group members.
             ! Group may become complete before fire is detected. Do not count these as complete.
             IF (T > Group_List(j)%Tdet) THEN
                IF (Group_List(j)%COMPLETE == 0) THEN
                   ! Tdoor: saves the time point when the group started to move towards the exit.
                   Group_List(j)%Tdoor = MAX(T,Group_List(j)%Tdet)
                END IF
                Group_List(j)%COMPLETE = 1
             END IF
          END IF
       END DO

       ! ========================================================
       ! CHANGE TARGET DOOR?
       ! ========================================================
       IF (T > 0.0_EB) THEN
          Change_Door_Loop: DO ie = 1, N_HUMANS
             HR => HUMAN(ie)
             ! Check and cycle if in stairs
             Check_Strs_Loop0: DO N = 1, N_STRS
                IF (EVAC_STRS(N)%IMESH == HR%IMESH) CYCLE Change_Door_Loop
             END DO Check_Strs_Loop0

             i_old_ffield    = HR%I_FFIELD
             name_old_ffield = TRIM(MESH_NAME(i_old_ffield))
             j  =  MAX(0,HR%GROUP_ID)   ! group index
             j1 = -MIN(0,HR%GROUP_ID)   ! lonely agent index

             ! If the group is not yet together, they are not moving towards the exit.
             IF (Group_List(j)%COMPLETE == 0) CYCLE Change_Door_Loop

             ! Agents start to move towards the exit after the pre-evacuation (reaction) time.
             IF (j == 0 .AND. T < HR%Tpre+HR%Tdet) CYCLE Change_Door_Loop   ! lonely agents
             IF (j > 0 .AND. T < Group_List(j)%Tpre + Group_List(j)%Tdoor) CYCLE Change_Door_Loop ! groups

             IF (Group_List(j)%GROUP_I_FFIELDS(i_egrid) == 0) THEN
                ! Test if this group/agent should update the exit door (exp. distribution)
                CALL RANDOM_NUMBER(rn)
                IF (rn > EXP(-DTSP/dt_group_door) ) THEN
                   i_target_old = HR%I_Target 
                   n_change_trials = n_change_trials + 1
                   CALL Change_Target_Door(nm, nm, ie, j, j1, i_egrid, 1, HR%X, HR%Y, i_target, color_index, i_new_ffield, HR)
                   IF (ABS(i_target_old) .NE. ABS(i_target)) THEN
                      n_change_doors = n_change_doors + 1
                      IF (i_target > 0) THEN
                         IF (i_target <= n_doors ) THEN
                            evel = SQRT((HR%X-EVAC_DOORS(i_target)%X)**2 + (HR%Y-EVAC_DOORS(i_target)%Y)**2)
                            evel = 50.0_EB*evel/EVAC_DOORS(i_target)%R_NTARGET + 1.0_EB
                            ii = MIN(50,MAX(1,INT(evel)))
                            EVAC_DOORS(i_target)%NTARGET(ii:50) = EVAC_DOORS(i_target)%NTARGET(ii:50) + 1
                         ELSE
                            evel = SQRT((HR%X-EVAC_EXITS(i_target-n_doors)%X)**2 + (HR%Y-EVAC_EXITS(i_target-n_doors)%Y)**2)
                            evel = 50.0_EB*evel/EVAC_EXITS(i_target-n_doors)%R_NTARGET + 1.0_EB
                            ii = MIN(50,MAX(1,INT(evel)))
                            EVAC_EXITS(i_target-n_doors)%NTARGET(ii:50) = EVAC_EXITS(i_target-n_doors)%NTARGET(ii:50) + 1
                         END IF
                      END IF
                      IF (i_target_old > 0) THEN
                         IF (i_target_old <= n_doors ) THEN
                            evel = SQRT((HR%X-EVAC_DOORS(i_target_old)%X)**2 + (HR%Y-EVAC_DOORS(i_target_old)%Y)**2)
                            evel = 50.0_EB*evel/EVAC_DOORS(i_target_old)%R_NTARGET + 1.0_EB
                            ii = MIN(50,MAX(1,INT(evel)))
                            EVAC_DOORS(i_target_old)%NTARGET(ii:50) = EVAC_DOORS(i_target_old)%NTARGET(ii:50) - 1
                         ELSE
                            evel = SQRT((HR%X-EVAC_EXITS(i_target_old-n_doors)%X)**2 + &
                                 (HR%Y-EVAC_EXITS(i_target_old-n_doors)%Y)**2)
                            evel = 50.0_EB*evel/EVAC_EXITS(i_target_old-n_doors)%R_NTARGET + 1.0_EB
                            ii = MIN(50,MAX(1,INT(evel)))
                            EVAC_EXITS(i_target_old-n_doors)%NTARGET(ii:50) = &
                                 EVAC_EXITS(i_target_old-n_doors)%NTARGET(ii:50) - 1
                         END IF
                      END IF
                   END IF
                ELSE ! Do not update door flow field at this time step
                   IF (j > 0) Group_Known_Doors(j)%I_Target = HR%I_Target 
                   IF (COLOR_METHOD == 5 .AND. j > 0) Color_Tmp(j) = HR%COLOR_INDEX
                   IF (COLOR_METHOD == 4 .AND. j > 0) Color_Tmp(j) = HR%COLOR_INDEX
                   IF (j > 0) Group_List(j)%GROUP_I_FFIELDS(i_egrid) = HR%I_FFIELD 
                END IF  ! update door flow field: is rn large enough?
             ELSE
                ! This group has already tried to change the field.
                IF (COLOR_METHOD == 5 .AND. j > 0) HR%COLOR_INDEX = Color_Tmp(j)
                IF (COLOR_METHOD == 4 .AND. j > 0) HR%COLOR_INDEX = Color_Tmp(j)
                HR%I_FFIELD    = Group_List(j)%GROUP_I_FFIELDS(i_egrid)
                HR%FFIELD_NAME = TRIM(MESH_NAME(HR%I_FFIELD))
                IF (j > 0) HR%I_Target = Group_Known_Doors(j)%I_Target
             END IF  ! group_i_field=0
          END DO Change_Door_Loop  ! loop over agents 
       END IF  ! t > 0

       ! ========================================================
       ! MOVE LOOP: The first step of the VV algorithm
       ! ========================================================
       ! ========================================================
       ! Step 1: Here the humans are moved to the new positions using the 
       ! present velocities.  
       ! After the first step the forces at the new postitions 
       ! are calculated and the velocities are updated.  
       ! (The 'dissipative' self-driving force contribution to the velocities 
       ! is updated self-consistently, but the other terms are not.)
       ! ========================================================
       TNOW15=SECOND()
       d_humans_min = HUGE(d_humans_min)
       d_walls_min  = HUGE(d_walls_min)

       EVAC_MOVE_LOOP: DO I=1, N_HUMANS
          HR=>HUMAN(I)
          ! hr%z is the real z-coordinate of the agent (inclines, stairs,etc),
          ! hr_z is z-coordinate of the main evac mesh
          hr_z = 0.5_EB*(ZS+ZF)
          LambdaW = LAMBDA_WALL
          A_Wall  = FAC_A_WALL*HR%A
          B_Wall  = FAC_B_WALL*HR%B
          GaMe    = NOISEME
          GaTh    = NOISETH
          GaCM    = NOISECM
          EVEL = SQRT(HR%U**2 + HR%V**2)
          L_Dead  = .FALSE.
          IF ( HR%IntDose >= 1.0_EB  ) THEN
             L_Dead = .TRUE.
             ! No random force for a dead person.
             GaTh = 0.0_EB
             ! No psychological force terms for a dead person.
             A_Wall = 0.0_EB
             IF (HR%Tpre /= HUGE(HR%Tpre)) THEN
                n_dead = n_dead+1
                WRITE (LU_EVACOUT,fmt='(a,i6,a,f8.2,a,i6)') ' EVAC: Person n:o', &
                     HR%ILABEL, ' dead at ', T, ' s, number of casualties ', n_dead
             END IF
             HR%Tpre = HUGE(HR%Tpre)
             HR%Tdet = HUGE(HR%Tdet)
             HR%COLOR_INDEX = EVAC_AVATAR_NCOLOR
          ELSE
             fed_max_alive = MAX(fed_max_alive,HR%IntDose)
          END IF
          fed_max = MAX(fed_max,HR%IntDose)  ! dead or alive
          hr_tau      = HR%Tau
          hr_tau_iner = HR%Tau_Iner
          ! Counterflow: increase motivation to go ahead
          IF (HR%Commitment > 0.01_EB) THEN
             hr_tau      = MAX(CF_MIN_TAU, &
                  HR%Commitment*CF_FAC_TAUS*hr_tau + (1.0_EB-HR%Commitment)*hr_tau)
             hr_tau_iner = MAX(CF_MIN_TAU_INER, &
                  HR%Commitment*CF_FAC_TAUS*hr_tau_iner + (1.0_EB-HR%Commitment)*hr_tau_iner)
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
          HR%W = 0.0_EB

          ! Check the smoke density for the detection by smoke
          IF (.NOT.L_Dead .AND. HUMAN_GRID(ii,jj)%SOOT_DENS > TDET_SMOKE_DENS) THEN
             HR%Tdet = MIN(HR%Tdet,T)
          END IF

          smoke_speed_fac = 1.0_EB
          IF (T > T_BEGIN) THEN
             ! Calculate Purser's fractional effective dose (FED)
             ! Note: Purser uses minutes, here dt is in seconds
             ! fed_dose = fed_lco*fed_vco2 + fed_lo
             HR%IntDose = DTSP*HUMAN_GRID(II,JJ)%FED_CO_CO2_O2 + HR%IntDose
             ! Smoke density vs speed
             ! Lund 2003, Report 3126 (Frantzich & Nilsson)
             ! v0(K) = v0*( 1 + (beta*K)/alpha ), where [K]=1/m ext.coeff
             ! beta=-0.057 m2/s, sigma=0.015, [-0.087,-0.028]
             ! alpha=0.706 m/s, sigma=0.069, [0.565,0.847]
             ! [SOOT_DENS] = mg/m3
             ! [K] = 1/m
             ! [MASS_EXTINCTION_COEFFICIENT] = m2/kg
             ! K = MASS_EXTINCTION_COEFFICIENT*SOOT_DENS*1.0E-6
             ! Visibility = 3/K  [m]
             smoke_beta  = -0.057_EB
             smoke_alpha = 0.706_EB
             smoke_speed_fac = 1.0_EB + (smoke_beta/smoke_alpha)* &
                  MASS_EXTINCTION_COEFFICIENT*1.0E-6_EB*HUMAN_GRID(II,JJ)%SOOT_DENS
             IF (MASS_EXTINCTION_COEFFICIENT*1.0E-6_EB*HUMAN_GRID(II,JJ)%SOOT_DENS > 3.0_EB/SMOKE_MIN_SPEED_VISIBILITY) THEN
                smoke_speed_fac = MIN(smoke_speed_fac, smoke_speed_fac*( 2.0_EB - &
                     ( MASS_EXTINCTION_COEFFICIENT*1.0E-6_EB*HUMAN_GRID(II,JJ)%SOOT_DENS - & 
                     3.0_EB/SMOKE_MIN_SPEED_VISIBILITY ) / (3.0_EB/SMOKE_MIN_SPEED_VISIBILITY) ) )
             END IF
             smoke_speed_fac = MAX(smoke_speed_fac, SMOKE_MIN_SPEED)
          END IF
          HR%v0_fac = smoke_speed_fac

          ! ========================================================
          ! Calculate persons prefered walking direction v0
          ! ========================================================
          NM_STRS_MESH = .FALSE.
          StrsMeshLoop: DO N = 1, N_STRS
             IF (EVAC_STRS(N)%IMESH==NM_now) THEN     
                NM_STRS_MESH = .TRUE.
                EXIT StrsMeshLoop
             END IF
          END DO StrsMeshLoop
          IF (TAU_CHANGE_V0 > 1.0E-12_EB) THEN
             ! Collision avoidance (incl. counterflow), do not update v0 on every time step.
             UBAR = HR%UBAR; VBAR = HR%VBAR
          ELSE
             CALL Find_Prefered_Direction(I, N, T, T_BEGIN, L_Dead, NM_STRS_MESH, &
                  II, JJ, IIX, JJY, XI, YJ, ZK, UBAR, VBAR, hr_tau, Tpre)
          END IF
          ! ========================================================
          ! Prefered walking direction v0 is now (UBAR,VBAR)
          ! ========================================================

          ! Inclines: Velocities are along the incline
          !           Coordinates are projected on the (x,y) plane
          cos_x = 1.0_EB
          cos_y = 1.0_EB
          speed_xm = HR%Speed
          speed_xp = HR%Speed
          speed_ym = HR%Speed
          speed_yp = HR%Speed
          IF (NM_STRS_MESH) THEN
             STRS_Indx = N
             STRP=>EVAC_STRS(N)     
          END IF
          !
          ! Check if an agent is on a spectator stand.
          SS_Loop1: DO j = 1, n_sstands
             ESS => EVAC_SSTANDS(j)
             IF (ESS%IMESH /= NM) CYCLE SS_Loop1
             IF (ESS%X1 > HR%X) CYCLE SS_Loop1
             IF (ESS%X2 < HR%X) CYCLE SS_Loop1
             IF (ESS%Y1 > HR%Y) CYCLE SS_Loop1
             IF (ESS%Y2 < HR%Y) CYCLE SS_Loop1
             cos_x = ESS%cos_x
             cos_y = ESS%cos_y
             FAC_V0_UP   = ESS%FAC_V0_UP
             FAC_V0_DOWN = ESS%FAC_V0_DOWN
             FAC_V0_HORI = ESS%FAC_V0_HORI
             IF (EVAC_PERSON_CLASSES(HR%IPC)%FAC_V0_HORI > 0.0_EB) THEN
                FAC_V0_HORI = EVAC_PERSON_CLASSES(HR%IPC)%FAC_V0_UP
             END IF
             IF (EVAC_PERSON_CLASSES(HR%IPC)%FAC_V0_UP > 0.0_EB) THEN
                IF ((ESS%H - ESS%H0) < 0.0_EB) THEN
                   FAC_V0_DOWN = EVAC_PERSON_CLASSES(HR%IPC)%FAC_V0_UP
                ELSE
                   FAC_V0_UP = EVAC_PERSON_CLASSES(HR%IPC)%FAC_V0_UP
                END IF
             END IF
             IF (EVAC_PERSON_CLASSES(HR%IPC)%FAC_V0_DOWN > 0.0_EB) THEN
                IF ((ESS%H - ESS%H0) < 0.0_EB) THEN
                   FAC_V0_UP = EVAC_PERSON_CLASSES(HR%IPC)%FAC_V0_DOWN 
                ELSE
                   FAC_V0_DOWN = EVAC_PERSON_CLASSES(HR%IPC)%FAC_V0_DOWN 
                END IF
             END IF
             SELECT CASE (ESS%IOR)
             CASE(-1)
                speed_xm = cos_x*HR%Speed*FAC_V0_DOWN
                speed_xp = cos_x*HR%Speed*FAC_V0_UP
                speed_ym = HR%Speed*FAC_V0_HORI
                speed_yp = HR%Speed*FAC_V0_HORI
                HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%X1-HR%X)/ABS(ESS%X1-ESS%X2)
             CASE(+1)
                speed_xm = cos_x*HR%Speed*FAC_V0_UP
                speed_xp = cos_x*HR%Speed*FAC_V0_DOWN
                speed_ym = HR%Speed*FAC_V0_HORI
                speed_yp = HR%Speed*FAC_V0_HORI
                HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%X2-HR%X)/ABS(ESS%X1-ESS%X2)
             CASE(-2)
                speed_xm = HR%Speed*FAC_V0_HORI
                speed_xp = HR%Speed*FAC_V0_HORI
                speed_ym = cos_y*HR%Speed*FAC_V0_DOWN
                speed_yp = cos_y*HR%Speed*FAC_V0_UP
                HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%Y1-HR%Y)/ABS(ESS%Y1-ESS%Y2)
             CASE(+2)
                speed_xm = HR%Speed*FAC_V0_HORI
                speed_xp = HR%Speed*FAC_V0_HORI
                speed_ym = cos_y*HR%Speed*FAC_V0_UP
                speed_yp = cos_y*HR%Speed*FAC_V0_DOWN
                HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%Y2-HR%Y)/ABS(ESS%Y1-ESS%Y2)
             END SELECT
             EXIT SS_Loop1
             HR%Z = 0.5_EB*(ZS+ZF)  ! The agent is not on any incline
          END DO SS_Loop1

          ! Set height and speed for an agent in stairs
          IF (NM_STRS_MESH) THEN 
             CALL GetStairSpeedAndZ(speed_xm, speed_xp, speed_ym, speed_yp, STRP, HR)
          END IF

          ! ========================================================
          ! Step 2:  The new velocities are calculated using the old forces.
          ! ========================================================
          ! (F_X,F_Y): wall forces + agent-agent forces (social + contact forces)
          ! Random forces and self-driving force are treated separately.
          U_new = HR%U + 0.5_EB*HR%F_X*DTSP/HR%Mass
          V_new = HR%V + 0.5_EB*HR%F_Y*DTSP/HR%Mass
          ! Rotational motion:
          Omega_new = HR%Omega + 0.5_EB*DTSP*HR%Torque/HR%M_iner

          ! Add self-driving force and torque
          EVEL = SQRT(UBAR**2 + VBAR**2)
          IF (EVEL > 0.0_EB) THEN
             !Inclines: U,V are the (x,y) plane projection velocities
             speed = speed_xp*(0.5_EB + SIGN(0.5_EB,UBAR)) + speed_xm*(0.5_EB - SIGN(0.5_EB,UBAR)) 
             speed = speed*HR%v0_fac
             U_new = U_new + 0.5_EB*(DTSP/hr_tau)*(speed*(UBAR/EVEL) - HR%U)
             HR%UBAR = UBAR

             speed = speed_yp*(0.5_EB + SIGN(0.5_EB,VBAR)) + speed_ym*(0.5_EB - SIGN(0.5_EB,VBAR)) 
             speed = speed*HR%v0_fac
             V_new = V_new + 0.5_EB*(DTSP/hr_tau)*(speed*(VBAR/EVEL) - HR%V)
             HR%VBAR = VBAR

             IF (VBAR >= 0.0_EB) THEN
                angle = ACOS(UBAR/EVEL)
             ELSE
                angle = 2.0_EB*Pi - ACOS(UBAR/EVEL)
             END IF

             ! Collison avoidance
             angle = angle + HR%angle_old
             DO WHILE (angle >= 2.0_EB*Pi)
                angle = angle - 2.0_EB*Pi
             END DO
             DO WHILE (angle < 0.0_EB)
                angle = angle + 2.0_EB*Pi
             END DO

             IF (angle == 2.0_EB*Pi) angle = 0.0_EB  ! angle is [0,2Pi)

             ! Rotational motion: J(dw/dt) = (J/t_iner)*( ((angle-angle_0/pi))w_0 - w )
             IF (ABS(angle-HR%angle) <= Pi ) THEN
                ! zero is not crossed.
                Omega_new = Omega_new + 0.5_EB*(DTSP/HR_tau_iner)*( (angle-HR%angle)*(Omega_0/Pi) - HR%Omega)
             ELSE
                ! zero is crossed
                Omega_new = Omega_new + 0.5_EB*(DTSP/HR_tau_iner)* &
                     ( (2.0_EB*Pi-ABS(angle-HR%angle))*SIGN(1.0_EB , HR%angle-angle)*(Omega_0/Pi) - HR%Omega)
             END IF
          ELSE  ! No v0
             U_new = U_new + 0.5_EB*(DTSP/hr_tau)* (- HR%U)
             V_new = V_new + 0.5_EB*(DTSP/hr_tau)* (- HR%V)
             ! Slow rotation down if no direction available, i.e., target Omega_0 is zero.
             Omega_new = Omega_new + 0.5_EB*(DTSP/HR_tau_iner)*(-HR%Omega)
             HR%UBAR = 0.0_EB
             HR%VBAR = 0.0_EB
          END IF
          !
          ! Check, if a new random force is needed on next time step.
          ! Poisson distribution, i.e., the agents do not have memory.
          ! P[No change during dt] = Exp(-DTSP/HR%Tau)
          IF ( GaTh > 0.0_EB .AND. T > T_BEGIN ) THEN
             CALL RANDOM_NUMBER(rn)
             IF ( rn > EXP(-DTSP/(0.2_EB*hr_tau)) ) HR%NewRnd = .TRUE.
             IF ( HR%NewRnd ) THEN
                HR%NewRnd = .FALSE.
                GaTh = GaTh*(HR%Speed/1.3_EB)**2
                ! GaTh is found to be more or less nice for speeds about 1.3 m/s
                HR%ksi = (GaussRand(GaMe, GaTh, GaCM))
                CALL RANDOM_NUMBER(rn)
                HR%eta = 2.0_EB*Pi*rn
                EVEL = SQRT(HR%U**2 + HR%V**2)
                ! Random noice variance: GaTh = Max( GaTh, GaTh*(100.0_EB-110.0_EB*Abs(EVEL/HR%Speed)) )
                ! Above GaTh is found to be more or less nice for speeds about 1.3 m/s
                ! Scale the random force by the current target speed v0 (note: Sqrt(GaTh) = std.dev.)
                ! Note: Speed reduction due to smoke is treated elsewhere.
                HR%ksi = ABS(HR%ksi)
                HR%ksi = MAX(HR%ksi, HR%ksi*10.0_EB*(1.0_EB-MIN(0.9_EB,EVEL/HR%Speed)))
             END IF
          END IF

          ! Add random force term
          IF ( GaTh > 0.0_EB .AND. T > T_BEGIN ) THEN
             U_new = U_new + 0.5_EB*DTSP*HR%v0_fac*HR%Mass*HR%ksi*COS(HR%eta)/HR%Mass
             V_new = V_new + 0.5_EB*DTSP*HR%v0_fac*HR%Mass*HR%ksi*SIN(HR%eta)/HR%Mass
             Omega_new = Omega_new + 0.5_EB*DTSP* 1.0_EB*SIGN(HR%ksi,HR%eta-Pi)
          END IF

          ! Check that velocities are not too large, i.e., unphysical (less than 10 m/s for humans)
          EVEL = SQRT(U_new**2 + V_new**2)
          IF ( EVEL > Vmax_timo ) THEN
             U_new = U_new*(Vmax_timo/EVEL)
             V_new = V_new*(Vmax_timo/EVEL)
          END IF
          ! Check that angular velocity is not too large
          IF ( ABS(Omega_new) > Omega_max ) THEN
             Omega_new = SIGN(Omega_max,Omega_new)
          END IF

          ! ========================================================
          ! Step 3:  The new coordinates are calculated using the new velocities.
          ! ========================================================
          X1 = HR%X + U_new*DTSP
          Y1 = HR%Y + V_new*DTSP
          HR%U = U_new
          HR%V = V_new
          ! The new body angle, should be on the interval [0,2Pi)
          A1 = HR%Angle + Omega_new*DTSP
          DO WHILE (A1 >= 2.0_EB*Pi)
             A1 = A1 - 2.0_EB*Pi
          END DO
          DO WHILE (A1 < 0.0_EB)
             A1 = A1 + 2.0_EB*Pi
          END DO
          HR%Omega = Omega_new
          !
          ! Check, if human is on an escalator and change the coordinates.
          SS_Loop1b: DO j = 1, n_sstands
             ESS => EVAC_SSTANDS(j)
             IF (ESS%IMESH == nm .AND. (ESS%X1 <= HR%X .AND. ESS%X2 >= HR%X) .AND. (ESS%Y1 <= HR%Y .AND. ESS%Y2 >= HR%Y) ) THEN
                cos_x = ESS%cos_x
                cos_y = ESS%cos_y
                SELECT CASE (ESS%IOR)
                CASE(-1)
                   X1 = X1 - cos_x*(ESS%Esc_SpeedUp-ESS%Esc_SpeedDn)*DTSP
                CASE(+1)
                   X1 = X1 - cos_x*(ESS%Esc_SpeedUp-ESS%Esc_SpeedDn)*DTSP
                CASE(-2)
                   Y1 = Y1 + cos_y*(ESS%Esc_SpeedUp-ESS%Esc_SpeedDn)*DTSP
                CASE(+2)
                   Y1 = Y1 - cos_y*(ESS%Esc_SpeedUp-ESS%Esc_SpeedDn)*DTSP
                END SELECT
                EXIT SS_Loop1b
             END IF
          END DO SS_Loop1b
          !
          ! In which grid cell is the agent, the new coordinates (t + dt)?
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
          HR%X_old = HR%X
          HR%Y_old = HR%Y

          ! Check, if the agent moves inside a solid object ==> might be an open
          ! vent or a 'sucking vent' used to calculate the flow fields.
          ! This is just to be fail safe.  If the user input is correct, this
          ! should never happen.
          IF ( SOLID(ICN) .AND. .NOT.OBSTRUCTION(I_OBST)%HIDDEN) THEN
             IF ( (Solid(ICX).AND. .NOT.OBSTRUCTION(I_OBSTX)%HIDDEN) .AND. .NOT. &
                  (Solid(ICY).AND. .NOT.OBSTRUCTION(I_OBSTY)%HIDDEN) ) THEN
                IF ( ii < iin ) THEN
                   tim_ic = cell_index(iin,jjn,kk)
                   CALL Get_iw(iin,jjn,kk,-1,tim_iw)
                   ibc = ijkw(5,tim_iw)
                   IF (SURFACE(IBC)%VEL> 0.0_EB .OR. BOUNDARY_TYPE(tim_iw)==OPEN_BOUNDARY) THEN
                      HR%X = X1 
                      HR%Y = HR%Y
                   END IF
                ELSE
                   tim_ic = cell_index(iin,jjn,kk)
                   CALL Get_iw(iin,jjn,kk,+1,tim_iw)
                   ibc = ijkw(5,tim_iw)
                   IF (SURFACE(IBC)%VEL> 0.0_EB .OR. BOUNDARY_TYPE(tim_iw)==OPEN_BOUNDARY) THEN
                      HR%X = X1 
                      HR%Y = HR%Y
                   END IF
                END IF
             ELSE IF ( (Solid(ICY).AND. .NOT.OBSTRUCTION(I_OBSTY)%HIDDEN) .AND. .NOT. &
                  (Solid(ICX).AND. .NOT.OBSTRUCTION(I_OBSTX)%HIDDEN) ) THEN
                IF ( jj < jjn ) THEN
                   tim_ic = cell_index(iin,jjn,kk)
                   CALL Get_iw(iin,jjn,kk,-2,tim_iw)
                   ibc = ijkw(5,tim_iw)
                   IF (SURFACE(IBC)%VEL> 0.0_EB .OR. BOUNDARY_TYPE(tim_iw)==OPEN_BOUNDARY) THEN
                      HR%X = HR%X
                      HR%Y = Y1
                   END IF
                ELSE
                   tim_ic = cell_index(iin,jjn,kk)
                   CALL Get_iw(iin,jjn,kk,+2,tim_iw)
                   ibc = ijkw(5,tim_iw)
                   IF (SURFACE(IBC)%VEL> 0.0_EB .OR. BOUNDARY_TYPE(tim_iw)==OPEN_BOUNDARY) THEN
                      HR%X = HR%X
                      HR%Y = Y1 
                   END IF
                END IF
             ELSE
                WRITE(MESSAGE,'(A,I4,A,2F8.2)') 'ERROR: Evacuate_Humans, Solid ICX and ICY, mesh ', nm, ' pos ',X1,Y1
             END IF
          ELSE
             ! Target cell is not a solid ==> move
             HR%X = X1
             HR%Y = Y1
          END IF
          HR%X = X1
          HR%Y = Y1
          HR%Angle = A1
          !
       END DO EVAC_MOVE_LOOP
       ! ========================================================
       ! MOVE LOOP ENDS HERE
       ! ========================================================

       ! ========================================================
       ! Check if persons are leaving this mesh via doors/exits
       ! and put these persons to the target.
       ! ========================================================
       IF (N_HUMANS > 0) CALL CHECK_DOORS(T,NM)
       IF (N_HUMANS > 0) CALL CHECK_EXITS(T,NM)

       ! ========================================================
       ! Check if persons are entering this mesh via a corr.
       ! ========================================================
       CALL CHECK_CORRS(T,NM,DTSP)

       ! ========================================================
       ! Add persons from entrys (specified flow rates)
       ! ========================================================
       IF (T > T_BEGIN ) THEN
          DO i = 1, N_ENTRYS
             CALL ENTRY_HUMAN(i,T,NM,istat)
          END DO
       END IF

       ! ========================================================
       ! Remove out-of-bounds persons (outside the grid)
       ! ========================================================
       IF (N_HUMANS > 0) CALL REMOVE_OUT_OF_GRIDS(T,NM)

       IF ( ICYC >= 0) THEN
          DTSP_new = EVAC_DT_STEADY_STATE
       ELSE
          DTSP_new = EVAC_DT_FLOWFIELD  ! Initialization phase
       END IF
       Speed_max  = 0.0_EB
       TUSED(15,NM)=TUSED(15,NM)+SECOND()-TNOW15  ! CPU timing

       ! ================================================
       ! Prepare to calculate the new forces, initialize
       ! different variables and arrys for the step 3 of 
       ! the SC-VV algorithm.
       ! ================================================

       ! ================================================
       ! Initialize group arrays for this main evac mesh.
       ! ================================================
       Group_List(:)%GROUP_SIZE  = 0
       Group_List(:)%GROUP_X = 0.0_EB
       Group_List(:)%GROUP_Y = 0.0_EB
       Group_List(:)%MAX_DIST_CENTER = 0.0_EB
       Group_List(:)%Tpre    = 0.0_EB
       Group_List(:)%Tdet    = HUGE(Group_List(:)%Tdet)

       DO j = 0, i33_dim
          Group_List(j)%GROUP_I_FFIELDS(i_egrid) = 0
       END DO
       DO i = 1, N_HUMANS
          HR=>HUMAN(I)
          j = MAX(0,HR%GROUP_ID)
          Group_List(j)%GROUP_SIZE = Group_List(j)%GROUP_SIZE + 1
          Group_List(j)%GROUP_X    = Group_List(j)%GROUP_X + HR%X
          Group_List(j)%GROUP_Y    = Group_List(j)%GROUP_Y + HR%Y
          Group_List(j)%GROUP_I_FFIELDS(i_egrid) = HR%i_ffield
          Group_List(j)%Tpre       = MAX(Group_List(j)%Tpre,HR%Tpre)
          Group_List(j)%Tdet       = MIN(Group_List(j)%Tdet,HR%Tdet)
       END DO
       !
       Group_List(1:)%GROUP_X = Group_List(1:)%GROUP_X / MAX(1,Group_List(1:)%GROUP_SIZE)
       Group_List(1:)%GROUP_Y = Group_List(1:)%GROUP_Y / MAX(1,Group_List(1:)%GROUP_SIZE)

       ! BLOCK_GRID_N(i,j): How many agents are in this grid cell
       BLOCK_GRID_N = 0
       Group_List(:)%MAX_DIST_CENTER = 0.0_EB
       DO i = 1, N_HUMANS
          HR => HUMAN(i)
          ! Which cell, new coordinates:
          IIN = FLOOR( CELLSI(FLOOR((HR%X-XS)*RDXINT)) + 1.0_EB)
          JJN = FLOOR( CELLSJ(FLOOR((HR%Y-YS)*RDYINT)) + 1.0_EB)
          BLOCK_GRID_N(IIN,JJN) = BLOCK_GRID_N(IIN,JJN) + 1
          j = MAX(0,HR%GROUP_ID)
          Group_List(j)%MAX_DIST_CENTER = MAX(Group_List(j)%MAX_DIST_CENTER, &
               SQRT((HR%X - Group_List(j)%GROUP_X)**2 + (HR%Y - Group_List(j)%GROUP_Y)**2))
       END DO
       Max_Humans_Cell = MAX(1,MAXVAL(BLOCK_GRID_N))

       IF (N_HUMANS > 0) THEN
          Group_List(0)%GROUP_SIZE = 1
          Group_List(0)%GROUP_X    = 0.5_EB*(XS+XF)
          Group_List(0)%GROUP_Y    = 0.5_EB*(YS+YF)
          Group_List(0)%Speed      = 1.0_EB
          Group_List(0)%IntDose    = 0.0_EB
          Group_List(0)%MAX_DIST_CENTER = 0.0_EB
          Group_List(0)%COMPLETE   = 1
       END IF
       ! Check if the groups are already gathered together or not.
       DO j = 1, i33_dim
          Group_List(j)%LIMIT_COMP = RADIUS_COMPLETE_0 + RADIUS_COMPLETE_1*Group_List(j)%GROUP_SIZE
          IF ( ((Group_List(j)%MAX_DIST_CENTER <=  Group_List(j)%LIMIT_COMP) .OR. &
               (Group_List(j)%COMPLETE == 1)) .AND. Group_List(j)%GROUP_SIZE > 0 ) THEN
             IF (T > Group_List(j)%Tdet) THEN
                IF (Group_List(j)%COMPLETE == 0) THEN
                   Group_List(j)%Tdoor = MAX(T,Group_List(j)%Tdet)
                END IF
                Group_List(j)%COMPLETE = 1
             END IF
          END IF
       END DO

       ! ========================================================
       ! Use blocks to speed up double loops in the force loop.
       ! ========================================================
       ALLOCATE(BLOCK_GRID(1:IBAR,1:JBAR,Max_Humans_Cell),STAT=IZERO)
       CALL ChkMemErr('EVACUATE_HUMANS','BLOCK_GRID',IZERO)
       BLOCK_GRID = 0
       ! Block_Grid_N(:,:) how many agents in this cell.
       ! Block_Grid(:,:,1-n) agent indeces

       BLOCK_GRID_N = 0
       DO i = 1, N_HUMANS
          HR => HUMAN(i)
          ! Which grid (block grid) cell, new coordinates:
          IIN = FLOOR( CELLSI(FLOOR((HR%X-XS)*RDXINT)) + 1.0_EB )
          JJN = FLOOR( CELLSJ(FLOOR((HR%Y-YS)*RDYINT)) + 1.0_EB )
          BLOCK_GRID_N(IIN,JJN) = BLOCK_GRID_N(IIN,JJN) + 1
          BLOCK_GRID(IIN,JJN,BLOCK_GRID_N(IIN,JJN)) = i
       END DO
       i_dx = 2*(INT((2.0_EB*0.3_EB+5.0_EB)/dx_min) + 1) + 1
       j_dy = 2*(INT((2.0_EB*0.3_EB+5.0_EB)/dy_min) + 1) + 1
       i_dx = MAX(1,i_dx)
       j_dy = MAX(1,j_dy)
       bl_max = i_dx*j_dy*Max_Humans_Cell
       ! Next list will contain the agents that should be looped over for the present agent.
       ALLOCATE(BLOCK_LIST(bl_max),STAT=IZERO)
       CALL ChkMemErr('EVACUATE_HUMANS','BLOCK_LIST',IZERO)

       ! ========================================================
       ! Step (3) of SC-VV starts here: Calculate new forces
       ! ========================================================
       TNOW13=SECOND()
       EVAC_FORCE_LOOP: DO I=1,N_HUMANS  
          HR => HUMAN(i)
          j  =  MAX(0,HR%GROUP_ID)   ! group index
          j1 = -MIN(0,HR%GROUP_ID)   ! lonely agent index
          ! hr%z is the real z-coordinate of the agent (inclines, stairs,etc),
          ! hr_z is z-coordinate of the main evac mesh
          hr_z = 0.5_EB*(ZS+ZF)
          ! Calculate the position and velocities of the shoulder cirles
          y_tmp(1) = HR%Y - COS(HR%angle)*HR%d_shoulder ! right
          x_tmp(1) = HR%X + SIN(HR%angle)*HR%d_shoulder ! right
          y_tmp(2) = HR%Y ! torso
          x_tmp(2) = HR%X ! torso
          y_tmp(3) = HR%Y + COS(HR%angle)*HR%d_shoulder ! left
          x_tmp(3) = HR%X - SIN(HR%angle)*HR%d_shoulder ! left
          r_tmp(1) = HR%r_shoulder
          r_tmp(2) = HR%r_torso
          r_tmp(3) = HR%r_shoulder
          u_tmp(1) = HR%U + COS(HR%angle)*HR%Omega*HR%d_shoulder ! right
          v_tmp(1) = HR%V + SIN(HR%angle)*HR%Omega*HR%d_shoulder ! right
          u_tmp(2) = HR%U ! torso
          v_tmp(2) = HR%V ! torso
          u_tmp(3) = HR%U - COS(HR%angle)*HR%Omega*HR%d_shoulder ! left
          v_tmp(3) = HR%V - SIN(HR%angle)*HR%Omega*HR%d_shoulder ! left
          hr_a = HR%A
          hr_b = HR%B

          Contact_F = 0.0_EB
          Social_F  = 0.0_EB
          LambdaW = LAMBDA_WALL
          A_Wall  = FAC_A_WALL*HR%A
          B_Wall  = FAC_B_WALL*HR%B
          GaMe    = NOISEME
          GaTh    = NOISETH
          GaCM    = NOISECM
          d_humans = HUGE(d_humans)
          d_walls  = HUGE(d_walls)
          L_Dead  = .FALSE.
          IF (HR%IntDose >= 1.0_EB) THEN
             L_Dead = .TRUE.
             ! No random force for a dead person.
             GaTh = 0.0_EB
             ! No psychological force terms for a dead person.
             A_Wall = 0.0_EB
             IF (HR%Tpre /= HUGE(HR%Tpre)) THEN
                n_dead = n_dead+1
                WRITE (LU_EVACOUT,fmt='(a,i6,a,f8.2,a,i6)') ' EVAC: Person n:o', HR%ILABEL, &
                     ' dead at ', T, ' s, number of casualties ', n_dead
             END IF
             HR%Tdet = HUGE(HR%Tdet)
             HR%Tpre = HUGE(HR%Tpre)
             HR%COLOR_INDEX = EVAC_AVATAR_NCOLOR
          END IF
          hr_tau = HR%Tau
          hr_tau_iner = HR%Tau_Iner
          ! =======================================================
          ! Speed dependent social force
          ! =======================================================
          hr_a =  HR%A*MAX(0.5_EB,(SQRT(HR%U**2+HR%V**2)/HR%Speed))
          A_Wall = MIN(A_Wall, FAC_A_WALL*hr_a)

          ! Counterflow: increase motivation to go ahead and decrease social force
          IF (HR%Commitment > 0.01_EB) THEN
             evel = MIN(1.0_EB,SQRT(HR%U**2+HR%V**2)/HR%Speed)
             evel = HR%Commitment*evel + (1.0_EB-HR%Commitment)*1.0_EB
             hr_tau      = MAX(CF_MIN_TAU, &
                  HR%Commitment*CF_FAC_TAUS*hr_tau + (1.0_EB-HR%Commitment)*hr_tau)
             hr_tau_iner = MAX(CF_MIN_TAU_INER, &
                  HR%Commitment*CF_FAC_TAUS*hr_tau_iner + (1.0_EB-HR%Commitment)*hr_tau_iner)
             hr_a =  HR%A*MAX(CF_MIN_A,evel)
             hr_b =  HR%B*MAX(CF_MIN_B,evel)
             A_Wall = MIN(A_Wall, CF_FAC_A_WALL*FAC_A_WALL*hr_a)
          END IF
          !
          ! Psychological force: cut-off when acceleration below 0.0001 m/s**2
          P2P_DIST_MAX = HR%B*LOG(HR%A/0.0001_EB)
          P2P_DIST_MAX = MIN(P2P_DIST_MAX, 5.0_EB)  ! 5.0 m is the maximum range of PP-force
          IF ( HR%SumForces2 > 0.1_EB ) THEN
             ! If large pressure then short range forces only (speed up)
             P2P_DIST_MAX = MIN( P2P_DIST_MAX, -HR%B*LOG(HR%SumForces2/(100.0_EB*HR%A)) )
          END IF
          P2P_DIST_MAX = MAX(P2P_DIST_MAX, 3.0_EB*HR%B)
          ! Next is the max distance for the collision avoidance, counterflow, etc.
          P2P_Suunta_MAX = MAX(P2P_DIST_MAX, 3.0_EB)

          ! Speed up the dead agent loop, only contact forces are needed.
          IF (L_Dead) P2P_DIST_MAX = 0.0_EB
          IF (L_Dead) P2P_Suunta_MAX = 0.0_EB
          ! Check if counterflow algorithm is not used.
          IF (TAU_CHANGE_V0 < 1.E-12_EB) P2P_Suunta_MAX = P2P_DIST_MAX

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
          CALL Find_Prefered_Direction(I, N, T+DTSP_new, T_BEGIN, L_Dead, NM_STRS_MESH, &
               IIN, JJN, IIX, JJY, XI, YJ, ZK, UBAR, VBAR, hr_tau, Tpre)
          ! ========================================================
          ! The prefered walking direction v0 is (UBAR,VBAR)
          ! ========================================================

          ! =======================================================
          ! Update the Block_Grid array search ranges
          ! =======================================================
          BLOCK_LIST = 0
          i_dx = INT((2.0_EB*0.3_EB+MAX(P2P_Suunta_MAX,P2P_DIST_MAX))/DX(IIN)) + 1
          j_dy = INT((2.0_EB*0.3_EB+MAX(P2P_Suunta_MAX,P2P_DIST_MAX))/DY(JJN)) + 1
          iio = MAX(1,IIN-i_dx)
          iie = MIN(IBAR,IIN+i_dx)
          jjo = MAX(1,JJN-j_dy)
          jje = MIN(JBAR,JJN+j_dy)
          ie_max = 0
          DO iix = iio, iie
             DO jjy = jjo, jje
                BG_Loop: DO j = 1, BLOCK_GRID_N(iix,jjy)
                   ie_max = ie_max + 1
                   BLOCK_LIST(MIN(ie_max,bl_max)) = BLOCK_GRID(iix,jjy,j)
                END DO BG_Loop
             END DO
          END DO
          IF (ie_max > bl_max ) THEN
             WRITE(MESSAGE,'(A,2I6)') 'ERROR: Evacuate_Humans, ie_max, bl_max ', ie_max, bl_max
             CALL SHUTDOWN(MESSAGE)
          END IF

          ! =======================================================
          ! Inclines: Velocities are along the incline
          !           Coordinates are projected on the (x,y) plane
          ! =======================================================
          cos_x = 1.0_EB
          cos_y = 1.0_EB
          speed_xm = HR%Speed
          speed_xp = HR%Speed
          speed_ym = HR%Speed
          speed_yp = HR%Speed
          IF (NM_STRS_MESH) THEN
             STRS_Indx = N
             STRP=>EVAC_STRS(N)     
          END IF
          ! Check if an agent is on a spectator stand.
          SS_Loop2: DO j = 1, n_sstands
             ESS => EVAC_SSTANDS(j)
             IF (ESS%IMESH == nm .AND. (ESS%X1 <= HR%X .AND. ESS%X2 >= HR%X) .AND. (ESS%Y1 <= HR%Y .AND. ESS%Y2 >= HR%Y) ) THEN
                cos_x = ESS%cos_x
                cos_y = ESS%cos_y
                FAC_V0_UP   = ESS%FAC_V0_UP
                FAC_V0_DOWN = ESS%FAC_V0_DOWN
                FAC_V0_HORI = ESS%FAC_V0_HORI
                IF (EVAC_PERSON_CLASSES(HR%IPC)%FAC_V0_HORI > 0.0_EB) THEN
                   FAC_V0_HORI = EVAC_PERSON_CLASSES(HR%IPC)%FAC_V0_UP
                END IF
                IF (EVAC_PERSON_CLASSES(HR%IPC)%FAC_V0_UP > 0.0_EB) THEN
                   IF ((ESS%H - ESS%H0) < 0.0_EB) THEN
                      FAC_V0_DOWN = EVAC_PERSON_CLASSES(HR%IPC)%FAC_V0_UP
                   ELSE
                      FAC_V0_UP = EVAC_PERSON_CLASSES(HR%IPC)%FAC_V0_UP
                   END IF
                END IF
                IF (EVAC_PERSON_CLASSES(HR%IPC)%FAC_V0_DOWN > 0.0_EB) THEN
                   IF ((ESS%H - ESS%H0) < 0.0_EB) THEN
                      FAC_V0_UP = EVAC_PERSON_CLASSES(HR%IPC)%FAC_V0_DOWN 
                   ELSE
                      FAC_V0_DOWN = EVAC_PERSON_CLASSES(HR%IPC)%FAC_V0_DOWN 
                   END IF
                END IF
                SELECT CASE (ESS%IOR)
                CASE(-1)
                   speed_xm = cos_x*HR%Speed*FAC_V0_DOWN
                   speed_xp = cos_x*HR%Speed*FAC_V0_UP
                   speed_ym = HR%Speed*FAC_V0_HORI
                   speed_yp = HR%Speed*FAC_V0_HORI
                   HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%X1-HR%X)/ABS(ESS%X1-ESS%X2)
                CASE(+1)
                   speed_xm = cos_x*HR%Speed*FAC_V0_UP
                   speed_xp = cos_x*HR%Speed*FAC_V0_DOWN
                   speed_ym = HR%Speed*FAC_V0_HORI
                   speed_yp = HR%Speed*FAC_V0_HORI
                   HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%X2-HR%X)/ABS(ESS%X1-ESS%X2)
                CASE(-2)
                   speed_xm = HR%Speed*FAC_V0_HORI
                   speed_xp = HR%Speed*FAC_V0_HORI
                   speed_ym = cos_y*HR%Speed*FAC_V0_DOWN
                   speed_yp = cos_y*HR%Speed*FAC_V0_UP
                   HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%Y1-HR%Y)/ABS(ESS%Y1-ESS%Y2)
                CASE(+2)
                   speed_xm = HR%Speed*FAC_V0_HORI
                   speed_xp = HR%Speed*FAC_V0_HORI
                   speed_ym = cos_y*HR%Speed*FAC_V0_UP
                   speed_yp = cos_y*HR%Speed*FAC_V0_DOWN
                   HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + (ESS%H-ESS%H0)*ABS(ESS%Y2-HR%Y)/ABS(ESS%Y1-ESS%Y2)
                END SELECT
                EXIT SS_Loop2
             END IF
             HR%Z = 0.5_EB*(ZS+ZF)  ! The agent is not on any incline
          END DO SS_Loop2
          ! =======================================================
          ! Inclines Ends
          ! =======================================================

          ! =======================================================
          ! Set height and speed for an agent in stairs
          ! =======================================================
          IF (NM_STRS_MESH) THEN 
            CALL GetStairSpeedAndZ(speed_xm, speed_xp, speed_ym, speed_yp,STRP,HR)
          END IF

          ! ========================================================
          ! AGENT-AGENT INTERACTION FORCES
          ! ========================================================
          ! Look for other agents, use blocks to speed up the loop.
          P2P_U      = 0.0_EB
          P2P_V      = 0.0_EB
          P2P_Torque = 0.0_EB
          TNOW14=SECOND()   ! person-person force loop timing
          ! ========================================================
          ! Collision avoidance (incl. counterflow)
          ! ========================================================
          ! Do not do this on every time step, do it every 0.1 s on the average by default.
          IF (TAU_CHANGE_V0 > 1.0E-12_EB) THEN
             CALL RANDOM_NUMBER(rnCF)
          ELSE
             rnCF = -1.0_EB
             angle_old = 0.0_EB
             commitment = 0.0_EB
          END IF
          Change_v0_rnCF0: IF ( rnCF > EXP(-DTSP/TAU_CHANGE_V0) ) THEN

             v_hr = MAX(0.1_EB,MIN(1.0_EB,SQRT(HR%U**2 + HR%V**2)/HR%Speed))
             theta_start = -ABS(THETA_SECTOR)
             theta_step = 2.0_EB*ABS(theta_start)/REAL(n_sectors-1,EB)
             evel = MAX(0.0_EB,MIN(1.0_EB,(HR%U**2+HR%V**2)/HR%Speed))
             theta_start = theta_start - 0.5_EB*(1.0_EB-evel)*MAX(0.0_EB,(90.0_EB+theta_start-0.5_EB*theta_step))
             theta_step = 2.0_EB*ABS(theta_start)/REAL(n_sectors-1,EB)
             IF(HR%UBAR**2+HR%VBAR**2 < 0.1_EB) THEN
                HR%UBAR = UBAR; HR%VBAR = VBAR
             END IF
             Sum_suunta = 0.0_EB
             DO iii = 1, n_sectors
                thetas(iii)    = -ABS(theta_start) + (iii-1)*theta_step ! Degrees
                cos_theta(iii) = COS(Pi*thetas(iii)/180._EB) ! Radians in Fortran functions
                sin_theta(iii) = SIN(Pi*thetas(iii)/180._EB)
                u_theta(iii) = cos_theta(iii)*UBAR - sin_theta(iii)*VBAR
                v_theta(iii) = sin_theta(iii)*UBAR + cos_theta(iii)*VBAR
             END DO
             cos_theta(n_sectors+1) = 1.0_EB ! v0 direction
             sin_theta(n_sectors+1) = 0.0_EB
             thetas(n_sectors+1)    = 0.0_EB
             u_theta(n_sectors+1)   = UBAR
             v_theta(n_sectors+1)   = VBAR
             angle_old       = 0.0_EB
             N_suunta        = 0 ! How many agents per sector
             N_suuntaCF      = 0 ! How many agents per sector (counterflow)
             N_suunta_back   = 0 ! How many agents "behind"
             N_suunta_backCF = 0 ! How many agents "behind" (counterflow)
          END IF Change_v0_rnCF0

          ! The force loop over the other agents starts here
          P2PLOOP: DO IE = 1, ie_max
             IF (BLOCK_LIST(ie) == I) CYCLE P2PLOOP
             HRE => HUMAN(BLOCK_LIST(ie))
             ! In stairs, only consider humans at the same, next and previous sub nodes (landings, stairs)
             IF (ABS(HRE%STR_SUB_INDX - HR%STR_SUB_INDX)>1) CYCLE P2PLOOP
             P2P_DIST = ((HRE%X-X1)**2 + (HRE%Y-Y1)**2)
             IF (P2P_DIST > (MAX(P2P_Suunta_Max,P2P_DIST_MAX)+HR%Radius+HRE%Radius)**2) CYCLE P2PLOOP
             P2P_DIST = SQRT(P2P_DIST)
             !
             ! Check, that the persons are seeing each other, i.e., there are no walls between.
             PP_see_each = See_each_other(nm, x1, y1, HRE%X, HRE%Y)
             IF (.NOT. PP_see_each) CYCLE P2PLOOP

             ! Collision avoidance, counterflow, etc.
             Change_v0_rnCF1: IF ( rnCF > EXP(-DTSP/TAU_CHANGE_V0) ) THEN
                EVEL = SQRT(UBAR**2 + VBAR**2)
                IF (VBAR >= 0.0_EB) THEN
                   angle_hr = ACOS(UBAR/EVEL)
                ELSE
                   angle_hr = 2.0_EB*Pi - ACOS(UBAR/EVEL)
                END IF
                IF (angle_hr == 2.0_EB*Pi) angle_hr = 0.0_EB  ! agent HR angle is [0,2Pi)
                
                d_shift = 2.0_EB*(0.5_EB-MIN(1.0_EB,SQRT(HR%U**2+HR%V**2)/HR%Speed))*HR%Radius
                hr_x = HR%X - UBAR*d_shift
                hr_y = HR%Y - VBAR*d_shift
                d_new = SQRT((HRE%X-hr_x)**2 + (HRE%Y-hr_y)**2)
                IF (d_shift < 0.0_EB .OR. d_new > HR%r_torso) THEN
                   IF ((HRE%Y-HR_Y) >= 0.0_EB) THEN
                      angle_hre = ACOS((HRE%X-HR_X)/d_new)
                   ELSE
                      angle_hre = 2.0_EB*Pi - ACOS((HRE%X-HR_X)/d_new)
                   END IF
                   IF (angle_hre == 2.0_EB*Pi) angle_hre = 0.0_EB  ! agent HRE angle is [0,2Pi)
                   angle_hre = angle_hr - angle_hre
                   
                   IF (angle_hre >= Pi) angle_hre = 2.0_EB*Pi - angle_hre
                   IF (angle_hre <= -Pi) angle_hre = 2.0_EB*Pi + angle_hre
                   ! Agent HRE is at this angle when measured from HR
                   ! If HRE is on the right hand side of HR then the angle is negative,
                   ! i.e., positive direction is anti-clockwise (as usual).
                   angle_hre = -180.0_EB*angle_hre/Pi  ! degrees
                   
                   v_hre = (HRE%X-HR%X)*UBAR + (HRE%Y-HR%Y)*VBAR
                   v_hre = v_hre/SQRT((HRE%X-HR%X)**2+(HRE%Y-HR%Y)**2) ! cos r_HRE vs v0
                   v_hre = MAX(0.5_EB,v_hre)  ! look further ahead than sideways
                   evel = MAX(0.2_EB,MIN(1.0_EB,SQRT(HR%U**2+HR%V**2)/HR%Speed))
                   tim_dist = MAX(0.2_EB,MIN(evel,v_hre))*P2P_Suunta_MAX + HR%Radius + HRE%Radius
                   evel = MAX(0.5_EB,MIN(1.0_EB,SQRT(HR%U**2+HR%V**2)/HR%Speed))
                   v_hre = MAX(0.5_EB,MIN(evel,v_hre))*P2P_Suunta_MAX ! look further ahead than sideways
                   IF (P2P_DIST < (v_hre+HR%Radius+HRE%Radius) ) THEN
                      vr_2r = HRE%UBAR*UBAR+HRE%VBAR*VBAR ! counterflow or not?
                      ! Which sector if any
                      iii = INT((angle_hre-(theta_start-1.5_EB*theta_step))/theta_step)
                      IF (iii > 0 .AND. iii < n_sectors+1) THEN
                         ! (UBAR,VBAR) are unit vectors
                         IF (P2P_DIST < tim_dist) N_suunta(iii) = N_suunta(iii) + 1
                         IF (P2P_DIST < tim_dist .AND. vr_2r <= -0.2_EB) N_suuntaCF(iii) = N_suuntaCF(iii) + 1
                         IF (vr_2r > 0.0_EB) THEN ! Same direction
                            v_hre = HRE%U*UBAR + HRE%V*VBAR  ! HRE speed along the v0 direction
                            v_hr  = MAX(0.0_EB,HR%U*UBAR + HR%V*VBAR)  ! HR speed along the v0 direction
                            v_hre = CONST_DF + FAC_DF*(MIN(HR%v0_fac*HR%Speed,v_hre) - &
                                 MIN(HR%v0_fac*HR%Speed,v_hr))
                         ELSE ! Counterflow
                            v_hre = HRE%U*UBAR + HRE%V*VBAR  ! HRE speed along the v0 direction
                            v_hre = -1.0_EB*(CONST_CF + FAC_CF*MAX(0.0_EB,-v_hre))
                         END IF
                         Sum_suunta(iii) = Sum_suunta(iii) + v_hre/MAX(0.2_EB,(P2P_DIST-HR%Radius-HRE%Radius))
                         IF (angle_hre > -0.5_EB*theta_step .AND. angle_hre < 0.5_EB*theta_step) THEN
                            ! The "additional" sector pointing on the v0 direction
                            N_suunta(n_sectors+1) = N_suunta(n_sectors+1) + 1
                            IF (vr_2r <= -0.2_EB) N_suuntaCF(n_sectors+1) = N_suuntaCF(n_sectors+1) + 1
                            Sum_suunta(n_sectors+1) = Sum_suunta(n_sectors+1) + &
                                 v_hre/MAX(0.2_EB,(P2P_DIST-HR%Radius-HRE%Radius))
                         END IF
                      END IF
                   END IF
                END IF
             END IF Change_v0_rnCF1
             IF ( P2P_DIST > (P2P_DIST_MAX+HR%Radius+HRE%Radius)**2 ) CYCLE P2PLOOP
             ! 
             ! ========================================================
             ! Calculate the combination of spring constant for the two agents
             ! ========================================================
             C_Yeff = (2.0_EB*HR%C_Young*2.0_EB*HRE%C_Young)/(2.0_EB*HR%C_Young+2.0_EB*HRE%C_Young)
             !
             ! ========================================================
             ! Angle dependent social force:
             ! ========================================================
             IF ( (HR%U**2 +HR%V**2) > 0.0_EB ) THEN
                CosPhiFac = ( (HRE%X-X1)*HR%U + (HRE%Y-Y1)*HR%V ) &
                     / ( SQRT((HRE%X-X1)**2 + (HRE%Y-Y1)**2)*SQRT(HR%U**2 +HR%V**2) )
                CosPhiFac = HR%Lambda + 0.5_EB*(1.0_EB-HR%Lambda)*(1.0_EB+CosPhiFac)
             ELSE
                CosPhiFac = 1.0_EB
             END IF

             ! Calculate the position and velocities of the shoulder cirles for HRE
             r_tmp(4) = HRE%r_shoulder ! right circle
             r_tmp(5) = HRE%r_torso    ! center circle
             r_tmp(6) = HRE%r_shoulder ! left circle
             y_tmp(4) = HRE%Y - COS(HRE%angle)*HRE%d_shoulder ! right circle
             x_tmp(4) = HRE%X + SIN(HRE%angle)*HRE%d_shoulder ! right circle
             y_tmp(5) = HRE%Y ! center circle
             x_tmp(5) = HRE%X ! center circle
             y_tmp(6) = HRE%Y + COS(HRE%angle)*HRE%d_shoulder ! left circle
             x_tmp(6) = HRE%X - SIN(HRE%angle)*HRE%d_shoulder ! left circle

             ! ========================================================
             ! Add psychological (social) force term
             ! ========================================================
             IF (.NOT. L_Dead) THEN
                Fc_x = 0.0_EB
                Fc_y = 0.0_EB
                ! Use the closest circles to calculate the psychological force
                DO iii = 1, 3
                   DO jjj = 4, 6
                      tim_dist = SQRT((x_tmp(iii)-x_tmp(jjj))**2 + (y_tmp(iii)-y_tmp(jjj))**2)
                      ! d_humans = Min( tim_dist-(r_tmp(iii)+r_tmp(jjj)) , d_humans )
                      ! Next is |vector1|*|vector2|
                      evel = SQRT((x_tmp(jjj)-x_tmp(iii))**2+(y_tmp(jjj)-y_tmp(iii))**2)* SQRT(u_tmp(iii)**2+v_tmp(iii)**2)
                      IF (evel > 0.0_EB) evel = ((x_tmp(jjj)-x_tmp(iii))*u_tmp(iii) + &
                           (y_tmp(jjj)-y_tmp(iii))*v_tmp(iii)) / evel   ! cos theta (scal_prod/(lenght1*length2)
                      IF (evel > 0.01_EB) THEN
                         d_humans = MIN( (tim_dist-(r_tmp(iii)+r_tmp(jjj))) /evel, d_humans)
                      ELSE
                         d_humans = MIN( (tim_dist-(r_tmp(iii)+r_tmp(jjj))) /0.01_EB , d_humans)
                      END IF

                      Fc_x1 = (x_tmp(iii)-x_tmp(jjj))*hr_a*CosPhiFac* &
                           EXP( -(tim_dist-( r_tmp(iii)+r_tmp(jjj) ))/hr_b )/tim_dist 
                      Fc_y1 = (y_tmp(iii)-y_tmp(jjj))*hr_a*CosPhiFac* &
                           EXP( -(tim_dist-( r_tmp(iii)+r_tmp(jjj) ))/hr_b )/tim_dist 
                      IF ( (Fc_x1**2+Fc_y1**2) > (Fc_x**2+Fc_y**2) ) THEN
                         Fc_x = Fc_x1
                         Fc_y = Fc_y1
                      END IF
                   END DO
                END DO
                P2P_U = P2P_U + Fc_x
                P2P_V = P2P_V + Fc_y
                Social_F = Social_F + SQRT(Fc_x**2 + Fc_y**2)
                Tc_z = 0.0_EB
                ! Calculate the torque due to the social force. Use the closest circles.
                DO jjj = 4, 6
                   ! First the right shoulder
                   tim_dist = SQRT( (x_tmp(jjj)-x_tmp(1))**2 + (y_tmp(jjj)-y_tmp(1))**2 )
                   Fc_x = (x_tmp(1)-x_tmp(jjj)) * hr_a*CosPhiFac*EXP( -(tim_dist - (r_tmp(1)+r_tmp(jjj)))/hr_b )/tim_dist
                   Fc_y = (y_tmp(1)-y_tmp(jjj)) * hr_a*CosPhiFac*EXP( -(tim_dist - (r_tmp(1)+r_tmp(jjj)))/hr_b )/tim_dist
                   IF ( ABS(Fc_y*(x_tmp(1)-HR%X) - Fc_x*(y_tmp(1)-HR%Y)) > ABS(Tc_z) ) THEN
                      Tc_z = Fc_y*(x_tmp(1)-HR%X) - Fc_x*(y_tmp(1)-HR%Y)
                   END IF
                END DO
                P2P_Torque = P2P_Torque + Tc_z
                Tc_z = 0.0_EB
                DO jjj = 4, 6
                   ! Then the left shoulder
                   tim_dist = SQRT( (x_tmp(jjj)-x_tmp(3))**2 + (y_tmp(jjj)-y_tmp(3))**2 )
                   Fc_x = (x_tmp(3)-x_tmp(jjj)) * hr_a*CosPhiFac*EXP( -(tim_dist - (r_tmp(3)+r_tmp(jjj)))/hr_b )/tim_dist
                   Fc_y = (y_tmp(3)-y_tmp(jjj)) * hr_a*CosPhiFac*EXP( -(tim_dist - (r_tmp(3)+r_tmp(jjj)))/hr_b )/tim_dist
                   IF ( ABS(Fc_y*(x_tmp(3)-HR%X) - Fc_x*(y_tmp(3)-HR%Y)) > ABS(Tc_z) ) THEN
                      Tc_z = Fc_y*(x_tmp(3)-HR%X) - Fc_x*(y_tmp(3)-HR%Y)
                   END IF
                END DO
                P2P_Torque = P2P_Torque + Tc_z
             END IF  ! dead or alive?

             ! ========================================================
             ! Add contact force terms
             ! ========================================================
             IF ( P2P_DIST <= (HR%Radius+HRE%Radius) ) THEN
                ! Calculate the velocities of the shoulder cirles
                v_tmp(4) = HRE%V + SIN(HRE%angle)*HRE%Omega*HRE%d_shoulder
                v_tmp(5) = HRE%V
                v_tmp(6) = HRE%V - SIN(HRE%angle)*HRE%Omega*HRE%d_shoulder
                u_tmp(4) = HRE%U + COS(HRE%angle)*HRE%Omega*HRE%d_shoulder
                u_tmp(5) = HRE%U
                u_tmp(6) = HRE%U - COS(HRE%angle)*HRE%Omega*HRE%d_shoulder

                DO iii = 1, 3
                   DO jjj = 4, 6
                      tim_dist = SQRT((x_tmp(iii)-x_tmp(jjj))**2 + (y_tmp(iii)-y_tmp(jjj))**2)
                      ! d_humans = Min( tim_dist-(r_tmp(iii)+r_tmp(jjj)) , d_humans )
                      ! Next is |vector1|*|vector2|
                      evel = SQRT((x_tmp(jjj)-x_tmp(iii))**2+(y_tmp(jjj)-y_tmp(iii))**2)* SQRT(u_tmp(iii)**2+v_tmp(iii)**2)
                      IF (evel > 0.0_EB) evel = ((x_tmp(jjj)-x_tmp(iii))*u_tmp(iii) + &
                           (y_tmp(jjj)-y_tmp(iii))*v_tmp(iii)) / evel   ! cos theta (scal_prod/(lenght1*length2)
                      IF (evel > 0.01_EB) THEN
                         d_humans = MIN( (tim_dist-(r_tmp(iii)+r_tmp(jjj))) /evel, d_humans)
                      ELSE
                         d_humans = MIN( (tim_dist-(r_tmp(iii)+r_tmp(jjj))) /0.01_EB , d_humans)
                      END IF
                      IF (tim_dist <= r_tmp(iii)+r_tmp(jjj) ) THEN
                         ! Circles are touching each others
                         Fc_x =(x_tmp(iii)-x_tmp(jjj))*C_Yeff*((r_tmp(iii)+r_tmp(jjj))-tim_dist)/tim_dist
                         Fc_y =(y_tmp(iii)-y_tmp(jjj))*C_Yeff*((r_tmp(iii)+r_tmp(jjj))-tim_dist)/tim_dist
                         Fc_x = Fc_x - FC_DAMPING*(u_tmp(iii)-u_tmp(jjj))*(x_tmp(iii)-x_tmp(jjj))/tim_dist
                         Fc_y = Fc_y - FC_DAMPING*(v_tmp(iii)-v_tmp(jjj))*(y_tmp(iii)-y_tmp(jjj))/tim_dist
                         Contact_F = Contact_F + SQRT(Fc_x**2 + Fc_y**2)
                         P2P_U = P2P_U + Fc_x
                         P2P_V = P2P_V + Fc_y
                         P2P_Torque = P2P_Torque + Fc_y*(x_tmp(iii)-HR%X) - Fc_x*(y_tmp(iii)-HR%Y)
                         scal_prod_over_rsqr = ((y_tmp(iii)-y_tmp(jjj))*(u_tmp(iii)-u_tmp(jjj)) - &
                              (x_tmp(iii)-x_tmp(jjj))*(v_tmp(iii)-v_tmp(jjj))) / (tim_dist**2)
                         IF (I_Fric_sw >= 1 ) THEN  ! This is the default
                            Fc_x = - HR%Kappa*((r_tmp(iii)+r_tmp(jjj))-tim_dist)* &
                                 ( (y_tmp(iii)-y_tmp(jjj)) * scal_prod_over_rsqr )
                            Fc_y = - HR%Kappa*((r_tmp(iii)+r_tmp(jjj))-tim_dist)* &
                                 (-(x_tmp(iii)-x_tmp(jjj)) * scal_prod_over_rsqr )
                            P2P_U = P2P_U + Fc_x
                            P2P_V = P2P_V + Fc_y
                            P2P_Torque = P2P_Torque + Fc_y*( (x_tmp(iii) + &
                                 (r_tmp(iii)/r_tmp(jjj))*(x_tmp(jjj)-x_tmp(iii)) ) - HR%X ) 
                            P2P_Torque = P2P_Torque - Fc_x*( (y_tmp(iii) + &
                                 (r_tmp(iii)/r_tmp(jjj))*(y_tmp(jjj)-y_tmp(iii)) ) - HR%Y ) 
                         ELSE
                            Fc_x = -HR%Gamma*( (y_tmp(iii)-y_tmp(jjj))*scal_prod_over_rsqr)
                            Fc_y = -HR%Gamma*(-(x_tmp(iii)-x_tmp(jjj))*scal_prod_over_rsqr)
                            P2P_U = P2P_U + Fc_x
                            P2P_V = P2P_V + Fc_y
                            P2P_Torque = P2P_Torque + Fc_y*( (x_tmp(iii) + &
                                 (r_tmp(iii)/r_tmp(jjj))*(x_tmp(jjj)-x_tmp(iii)) ) - HR%X ) 
                            P2P_Torque = P2P_Torque - Fc_x*( (y_tmp(iii) + &
                                 (r_tmp(iii)/r_tmp(jjj))*(y_tmp(jjj)-y_tmp(iii)) ) - HR%Y ) 
                         END IF
                      END IF
                   END DO
                END DO
             END IF  ! contact forces?
          END DO P2PLOOP
          TUSED(14,NM)=TUSED(14,NM)+SECOND()-TNOW14
          ! ========================================================
          ! PERSON-PERSON INTERACTION FORCES ENDS HERE
          ! ========================================================

          ! ========================================================
          ! THE PERSON-WALL FORCES
          ! ========================================================
          !
          ! Walls are looked for the body circle
          CALL Find_walls(nm, x1, y1, HR%Radius, P2P_DIST_MAX, HR%SKIP_WALL_FORCE_IOR, d_xy, FoundWall_xy, istat)

          ! ========================================================
          ! Collision avoidance, counterflow, etc.
          ! ========================================================
          Change_v0_rnCF2: IF ( rnCF > EXP(-DTSP/TAU_CHANGE_V0) ) THEN
             v_hr  = MAX(0.1_EB,SQRT(HR%U**2+HR%V**2)/HR%Speed)
             tim_dist = 0.0_EB
             DO iii = 1, n_sectors
                ! Awoid walls, do not take a direction where there is a wall closer than
                ! d_perp = 0.6 m (perpendicular).
                IF (ABS(sin_theta(iii)) > 0.0001_EB) THEN
                   x11 = HR%X + u_theta(iii)*MIN(P2P_Suunta_MAX, 0.6_EB/ABS(sin_theta(iii)))
                   y11 = HR%Y + v_theta(iii)*MIN(P2P_Suunta_MAX, 0.6_EB/ABS(sin_theta(iii)))
                ELSE ! straight ahead
                   x11 = HR%X + u_theta(iii)*P2P_Suunta_MAX
                   y11 = HR%Y + v_theta(iii)*P2P_Suunta_MAX
                END IF
                P2P_DIST = SQRT((HR%X-x11)**2 + (HR%Y-y11)**2) - HR%Radius
                x11 = MAX(XS,MIN(XF,x11))
                y11 = MAX(YS,MIN(YF,y11))
                PP_see_each = See_each_other(nm, HR%X, HR%Y, x11, y11)
                IF(.NOT.PP_see_each) THEN
                   vr_2r = -FAC_1_WALL*P2P_Suunta_MAX*v_hr
                   Sum_suunta(iii) = Sum_suunta(iii) + vr_2r/MAX(0.2_EB,P2P_DIST)
                END IF
             END DO

             IF (FoundWall_xy(1)) tim_dist = tim_dist - 3.0_EB/MAX(0.3_EB,ABS(HR%X-d_xy(1))-HR%Radius)
             IF (FoundWall_xy(2)) tim_dist = tim_dist - 3.0_EB/MAX(0.3_EB,ABS(HR%X-d_xy(2))-HR%Radius)
             IF (FoundWall_xy(3)) tim_dist = tim_dist - 3.0_EB/MAX(0.3_EB,ABS(HR%Y-d_xy(3))-HR%Radius)
             IF (FoundWall_xy(4)) tim_dist = tim_dist - 3.0_EB/MAX(0.3_EB,ABS(HR%Y-d_xy(4))-HR%Radius)
!!$             If (FoundWall_xy(1)) tim_dist = tim_dist - 3.0_EB/Max(0.3_EB,Abs(HR%X-d_xy(1))-HR%Radius)
!!$             If (FoundWall_xy(2)) tim_dist = tim_dist - 3.0_EB/Max(0.3_EB,Abs(HR%X-d_xy(2))-HR%Radius)
!!$             If (FoundWall_xy(3)) tim_dist = tim_dist - 3.0_EB/Max(0.3_EB,Abs(HR%Y-d_xy(3))-HR%Radius)
!!$             If (FoundWall_xy(4)) tim_dist = tim_dist - 3.0_EB/Max(0.3_EB,Abs(HR%Y-d_xy(4))-HR%Radius)

             IF (FoundWall_xy(1) .AND. ABS(HR%X-d_xy(1))-HR%Radius < 0.1_EB) THEN
                DO iii = 1, n_sectors
                   IF(u_theta(iii) < -0.10_EB) THEN
                      vr_2r = u_theta(iii)/MAX(0.1_EB, ABS(HR%X-HR%Radius-d_xy(1)))
                      Sum_suunta(iii) = Sum_suunta(iii) - ABS(P2P_Suunta_MAX*FAC_2_WALL*vr_2r)
                   END IF
                END DO
             END IF
             IF (FoundWall_xy(2) .AND. ABS(HR%X-d_xy(2))-HR%Radius < 0.1_EB) THEN
                DO iii = 1, n_sectors
                   IF(u_theta(iii) > +0.10_EB) THEN
                      vr_2r = -u_theta(iii)/MAX(0.1_EB, ABS(HR%X+HR%Radius-d_xy(2)))
                      Sum_suunta(iii) = Sum_suunta(iii) - ABS(P2P_Suunta_MAX*FAC_2_WALL*vr_2r)
                   END IF
                END DO
             END IF
             IF (FoundWall_xy(3) .AND. ABS(HR%Y-d_xy(3))-HR%Radius < 0.1_EB) THEN
                DO iii = 1, n_sectors
                   IF(v_theta(iii) < -0.10_EB) THEN
                      vr_2r = v_theta(iii)/MAX(0.1_EB, ABS(HR%Y-HR%Radius-d_xy(3)))
                      Sum_suunta(iii) = Sum_suunta(iii) - ABS(P2P_Suunta_MAX*FAC_2_WALL*vr_2r)
                   END IF
                END DO
             END IF
             IF (FoundWall_xy(4) .AND. ABS(HR%Y-d_xy(4))-HR%Radius < 0.1_EB) THEN
                DO iii = 1, n_sectors
                   IF(v_theta(iii) > +0.10_EB) THEN
                      vr_2r = -v_theta(iii)/MAX(0.1_EB, ABS(HR%Y+HR%Radius-d_xy(4)))
                      Sum_suunta(iii) = Sum_suunta(iii) - ABS(P2P_Suunta_MAX*FAC_2_WALL*vr_2r)
                   END IF
                END DO
             END IF
!!$             Do iii = 1, n_sectors+1
!!$                ! Prefer flow field direction
!!$                Sum_suunta(iii) = Sum_suunta(iii) + FAC_V0_DIR*v_hr*((UBAR*u_theta(iii)+VBAR*v_theta(iii)))
!!$             End Do
             Sum_suunta(n_sectors+1) = Sum_suunta(n_sectors+1) + ABS(FAC_V0_DIR)*v_hr
             DO iii = 1, n_sectors
                IF (N_suuntaCF(n_sectors+1) < 1) THEN
                   ! No counterflow: Prefer left (and straight ahead)
                   ! Sum_suunta(iii) = Sum_suunta(iii) + Sign(1.0_EB,thetas(iii))* &
                   !     FAC_V0_DIR*v_hr
                ELSE
                   ! Counterflow: Prefer right (and straight ahead)
                   Sum_suunta(iii) = Sum_suunta(iii) - SIGN(1.0_EB,thetas(iii))*FAC_V0_DIR
                END IF
             END DO

             IF (N_suunta(n_sectors+1) < 1) Sum_suunta(n_sectors+1) = Sum_suunta(n_sectors+1) + 50.0_EB  ! Empty space ahead
             IF ((N_suuntaCF(n_sectors+1)) < 1) THEN
                ! No counterflow, prefer v0 direction, i.e., "stay on line"
                Sum_suunta(n_sectors+1) = Sum_suunta(n_sectors+1) + &
                     (FAC_NOCF + FAC_V0_NOCF*v_hr)*(N_suunta(n_sectors+1))
             END IF
             i_suunta_max = n_sectors + 1
             Sum_suunta_max = -HUGE(Sum_suunta_max)
             DO iii = 1, n_sectors+1
                IF (Sum_suunta(iii) > Sum_suunta_max) THEN
                   UBAR = u_theta(iii)
                   VBAR = v_theta(iii)
                   Sum_suunta_max = Sum_suunta(iii)
                   i_suunta_max = iii
                END IF
             END DO
             ! If counterflow then try to pass sideways and decrease tau and tau_iner.
             ! Use this if there are more counterflow agents than downstrean agents.
             angle_old = 0.0_EB
             commitment = 0.0_EB
             IF (v_hr < 0.3_EB .AND. (N_suunta(n_sectors+1)-N_suuntaCF(n_sectors+1)) < &
                  N_suuntaCF(n_sectors+1)-0) THEN
                commitment = MAX(0.5_EB,REAL(SUM(N_suuntaCF(1:n_sectors+1)),EB)/MAX(1,SUM(N_suunta(1:n_sectors+1))))
                angle_old = -SIGN(1.0_EB,thetas(i_suunta_max))*MAX(0.1_EB,ABS(thetas(i_suunta_max)))*Pi/180.0_EB 
                angle_old = angle_old + 85.0_EB*Pi/180.0_EB
             ELSE IF (v_hr < 0.3_EB .AND. SUM(N_suuntaCF(1:n_sectors)) >= 1 .AND. tim_dist < -15.0_EB) THEN
                commitment = MAX(0.5_EB,REAL(SUM(N_suuntaCF(1:n_sectors)),EB)/MAX(1,SUM(N_suunta(1:n_sectors))))
                angle_old = -SIGN(1.0_EB,thetas(i_suunta_max))*MAX(0.1_EB,ABS(thetas(i_suunta_max)))*Pi/180.0_EB 
                angle_old = angle_old + 85.0_EB*Pi/180.0_EB
             ELSE
                angle_old = 0.0_EB
                commitment = 0.0_EB
             END IF
          ELSE ! Change_v0_rnCF2
             ! Do not change direction during this time step, use the previous direction
             IF (TAU_CHANGE_V0 > 1.0E-12_EB) THEN
                UBAR = HR%UBAR
                VBAR = HR%VBAR
                angle_old = HR%angle_old
                commitment = HR%Commitment
             ELSE
                HR%UBAR = UBAR
                HR%VBAR = VBAR
                angle_old = 0.0_EB
                commitment = 0.0_EB
             END IF
          END IF Change_v0_rnCF2 ! Collision avoidance ends
          HR%angle_old = angle_old
          HR%Commitment = commitment

          CALL Wall_SocialForces(nm, x_tmp, y_tmp, r_tmp, P2P_DIST_MAX, d_xy, P2P_U, P2P_V, Social_F, FoundWall_xy)

          CALL Wall_ContactForces(nm, x_tmp(1), y_tmp(1), r_tmp(1), u_tmp(1), v_tmp(1), d_xy, &
               P2P_U, P2P_V, P2P_Torque, Contact_F, d_walls, FoundWall_xy)
          CALL Wall_ContactForces(nm, x_tmp(2), y_tmp(2), r_tmp(2), u_tmp(2), v_tmp(2), d_xy, &
               P2P_U, P2P_V, P2P_Torque, Contact_F, d_walls, FoundWall_xy)
          CALL Wall_ContactForces(nm, x_tmp(3), y_tmp(3), r_tmp(3), u_tmp(3), v_tmp(3), d_xy, &
               P2P_U, P2P_V, P2P_Torque, Contact_F, d_walls, FoundWall_xy)

          ! Add forces from the door case
          CALL Door_Forces(nm, x_tmp, y_tmp, r_tmp, u_tmp, v_tmp, p2p_dist_max, d_xy,&
               P2P_U, P2P_V, Social_F, Contact_F, P2P_Torque, FoundWall_xy)

          ! Add wall corner - person forces
          ! top right corner (x > x_human, y > y_human)
          x_now = -DX(iin+1)
          Loop_px: DO ii = iin, IBAR, +1
             x_now = x_now + DX(ii)
             IF (x_now-HR%Radius > P2P_DIST_MAX) EXIT Loop_px
             y_now = -DY(jjn-1)
             Loop_pxpy: DO jj = jjn, JBAR
                y_now = y_now + DY(jj)
                IF (SQRT(x_now**2 + y_now**2)-HR%Radius > P2P_DIST_MAX) EXIT Loop_pxpy
                tim_ic  = cell_index(ii,jj,kkn)   ! present
                tim_iwx = wall_index(tim_ic, +1)  ! right
                IF (tim_iwx>0) THEN
                   I_OBST = OBST_INDEX_W(tim_IWX)
                   IF (OBSTRUCTION(I_OBST)%HIDDEN) tim_iwx = 0
                END IF
                tim_iwy = wall_index(tim_ic, +2)  ! up
                IF (tim_iwy>0) THEN
                   I_OBST = OBST_INDEX_W(tim_IWY)
                   IF (OBSTRUCTION(I_OBST)%HIDDEN) tim_iwy = 0
                END IF
                tim_ic  = cell_index(ii,jj+1,kkn) ! one cell up
                tim_iw  = wall_index(tim_ic, +1)  ! up and right
                IF (tim_iw>0) THEN
                   I_OBST = OBST_INDEX_W(tim_IW)
                   IF (OBSTRUCTION(I_OBST)%HIDDEN) tim_iw = 0
                END IF
                tim_ic2 = cell_index(ii+1,jj,kkn) ! one cell right
                tim_iw2 = wall_index(tim_ic2, +2) ! right and up
                IF (tim_iw2>0) THEN
                   I_OBST = OBST_INDEX_W(tim_IW2)
                   IF (OBSTRUCTION(I_OBST)%HIDDEN) tim_iw2 = 0
                END IF
                IF (tim_iwy /= 0) EXIT Loop_pxpy  ! 
                IF ( (tim_iwx==0).AND.(tim_iwy==0).AND.(tim_iw/=0 .OR. tim_iw2/=0) ) THEN
                   IF (tim_iw/=0) THEN
                      ! first y-direction then x-direction
                      x11 = xw(tim_iw )                 ! corner point x
                      y11 = yw(tim_iw )-0.5_EB*DY(jj+1) ! corner point y
                   ELSE
                      ! first x-direction then y-direction
                      x11 = xw(tim_iw2)-0.5_EB*DX(ii+1) ! corner point x
                      y11 = yw(tim_iw2)                 ! corner point y
                   END IF

                   CALL Corner_Forces(x1, y1, x11, y11, p2p_dist_max, P2P_U, P2P_V, Social_F, &
                        Contact_F, P2P_Torque, d_walls, x_tmp, y_tmp, r_tmp, u_tmp, v_tmp, istat)
                   EXIT Loop_pxpy
                END IF
             END DO Loop_pxpy
          END DO Loop_px

          ! top left corner (x < x_human, y > y_human)
          x_now = -DX(iin-1)
          Loop_mx: DO ii = iin, 1, -1
             x_now = x_now + DX(ii)
             IF (x_now-HR%Radius > P2P_DIST_MAX) EXIT Loop_mx
             y_now = -DY(jjn)
             Loop_mxpy: DO jj = jjn, JBAR
                y_now = y_now + DY(jj)
                IF (SQRT(x_now**2 + y_now**2)-HR%Radius > P2P_DIST_MAX) EXIT Loop_mxpy
                tim_ic  = cell_index(ii,jj,kkn)   ! present
                tim_iwx = wall_index(tim_ic, -1)  ! left
                IF (tim_iwx>0) THEN
                   I_OBST = OBST_INDEX_W(tim_IWX)
                   IF (OBSTRUCTION(I_OBST)%HIDDEN) tim_iwx = 0
                END IF
                tim_iwy = wall_index(tim_ic, +2)  ! up
                IF (tim_iwy>0) THEN
                   I_OBST = OBST_INDEX_W(tim_IWY)
                   IF (OBSTRUCTION(I_OBST)%HIDDEN) tim_iwy = 0
                END IF
                tim_ic  = cell_index(ii,jj+1,kkn) ! one cell up
                tim_iw  = wall_index(tim_ic, -1)  ! up and left 
                IF (tim_iw>0) THEN
                   I_OBST = OBST_INDEX_W(tim_IW)
                   IF (OBSTRUCTION(I_OBST)%HIDDEN) tim_iw = 0
                END IF
                tim_ic2 = cell_index(ii-1,jj,kkn) ! one cell left
                tim_iw2 = wall_index(tim_ic2, +2) ! left and up
                IF (tim_iw2>0) THEN
                   I_OBST = OBST_INDEX_W(tim_IW2)
                   IF (OBSTRUCTION(I_OBST)%HIDDEN) tim_iw2 = 0
                END IF
                IF (tim_iwy /= 0) EXIT Loop_mxpy
                IF ( (tim_iwx==0).AND.(tim_iwy==0).AND.(tim_iw/=0 .OR. tim_iw2/=0) ) THEN
                   IF (tim_iw/=0) THEN
                      ! first y-direction then x-direction
                      x11 = xw(tim_iw )                 ! corner point x
                      y11 = yw(tim_iw )-0.5_EB*DY(jj+1) ! corner point y
                   ELSE
                      ! first x-direction then y-direction
                      x11 = xw(tim_iw2)+0.5_EB*DX(ii-1) ! corner point x
                      y11 = yw(tim_iw2)                 ! corner point y
                   END IF

                   CALL Corner_Forces(x1, y1, x11, y11, p2p_dist_max, P2P_U, P2P_V, Social_F, &
                        Contact_F, P2P_Torque, d_walls, x_tmp, y_tmp, r_tmp, u_tmp, v_tmp, istat)
                   EXIT Loop_mxpy
                END IF
             END DO Loop_mxpy
          END DO Loop_mx

          ! bottom right corner (x > x_human, y < y_human)
          x_now = -DX(iin+1)
          Loop_py: DO ii = iin, IBAR, +1
             x_now = x_now + DX(ii)
             IF (x_now-HR%Radius > P2P_DIST_MAX) EXIT Loop_py
             y_now = -DY(jjn-1)
             Loop_pxmy: DO jj = jjn, 1, -1
                y_now = y_now + DY(jj)
                IF (SQRT(x_now**2 + y_now**2)-HR%Radius > P2P_DIST_MAX) EXIT Loop_pxmy
                tim_ic  = cell_index(ii,jj,kkn)   ! present
                tim_iwx = wall_index(tim_ic, +1)  ! right
                IF (tim_iwx>0) THEN
                   I_OBST = OBST_INDEX_W(tim_IWX)
                   IF (OBSTRUCTION(I_OBST)%HIDDEN) tim_iwx = 0
                END IF
                tim_iwy = wall_index(tim_ic, -2)  ! down
                IF (tim_iwy>0) THEN
                   I_OBST = OBST_INDEX_W(tim_IWY)
                   IF (OBSTRUCTION(I_OBST)%HIDDEN) tim_iwy = 0
                END IF
                tim_ic  = cell_index(ii,jj-1,kkn) ! one cell down
                tim_iw  = wall_index(tim_ic, +1)  ! down and right
                IF (tim_iw>0) THEN
                   I_OBST = OBST_INDEX_W(tim_IW)
                   IF (OBSTRUCTION(I_OBST)%HIDDEN) tim_iw = 0
                END IF
                tim_ic2 = cell_index(ii+1,jj,kkn) ! one cell right
                tim_iw2 = wall_index(tim_ic2, -2) ! right and down
                IF (tim_iw2>0) THEN
                   I_OBST = OBST_INDEX_W(tim_IW2)
                   IF (OBSTRUCTION(I_OBST)%HIDDEN) tim_iw2 = 0
                END IF
                IF (tim_iwy /= 0) EXIT Loop_pxmy
                IF ( (tim_iwx==0).AND.(tim_iwy==0).AND.(tim_iw/=0 .OR. tim_iw2/=0) ) THEN
                   IF (tim_iw/=0) THEN
                      ! first y-direction then x-direction
                      x11 = xw(tim_iw )                 ! corner point x
                      y11 = yw(tim_iw )+0.5_EB*DY(jj-1) ! corner point y
                   ELSE
                      ! first x-direction then y-direction
                      x11 = xw(tim_iw2)-0.5_EB*DX(ii+1) ! corner point x
                      y11 = yw(tim_iw2)                 ! corner point y
                   END IF

                   CALL Corner_Forces(x1, y1, x11, y11, p2p_dist_max, P2P_U, P2P_V, Social_F, &
                        Contact_F, P2P_Torque, d_walls, x_tmp, y_tmp, r_tmp, u_tmp, v_tmp, istat)
                   EXIT Loop_pxmy
                END IF
             END DO Loop_pxmy
          END DO Loop_py

          ! bottom left corner (x < x_human, y < y_human)
          x_now = -DX(iin-1)
          Loop_my: DO ii = iin, 1, -1
             x_now = x_now + DX(ii)
             IF (x_now-HR%Radius > P2P_DIST_MAX) EXIT Loop_my
             y_now = -DY(jjn-1)
             Loop_mxmy: DO jj = jjn, 1, -1
                y_now = y_now + DY(jj)
                IF (SQRT(x_now**2+y_now**2)-HR%Radius > P2P_DIST_MAX) EXIT Loop_mxmy
                tim_ic  = cell_index(ii,jj,kkn)   ! present
                tim_iwx = wall_index(tim_ic, -1)  ! left
                IF (tim_iwx>0) THEN
                   I_OBST = OBST_INDEX_W(tim_IWX)
                   IF (OBSTRUCTION(I_OBST)%HIDDEN) tim_iwx = 0
                END IF
                tim_iwy = wall_index(tim_ic, -2)  ! down
                IF (tim_iwy>0) THEN
                   I_OBST = OBST_INDEX_W(tim_IWY)
                   IF (OBSTRUCTION(I_OBST)%HIDDEN) tim_iwy = 0
                END IF
                tim_ic  = cell_index(ii,jj-1,kkn) ! once cell down
                tim_iw  = wall_index(tim_ic, -1)  ! down and left
                IF (tim_iw>0) THEN
                   I_OBST = OBST_INDEX_W(tim_IW)
                   IF (OBSTRUCTION(I_OBST)%HIDDEN) tim_iw = 0
                END IF
                tim_ic2 = cell_index(ii-1,jj,kkn) ! one cell left
                tim_iw2 = wall_index(tim_ic2, -2) ! left and down
                IF (tim_iw2>0) THEN
                   I_OBST = OBST_INDEX_W(tim_IW2)
                   IF (OBSTRUCTION(I_OBST)%HIDDEN) tim_iw2 = 0
                END IF
                IF (tim_iwy /= 0) EXIT Loop_mxmy
                IF ( (tim_iwx==0).AND.(tim_iwy==0).AND.(tim_iw/=0 .OR. tim_iw2/=0) ) THEN
                   IF (tim_iw/=0) THEN
                      ! first y-direction then x-direction
                      x11 = xw(tim_iw )                 ! corner point x
                      y11 = yw(tim_iw )+0.5_EB*DY(jj-1) ! corner point y
                   ELSE
                      ! first x-direction then y-direction
                      x11 = xw(tim_iw2)+0.5_EB*DX(ii-1) ! corner point x
                      y11 = yw(tim_iw2)                 ! corner point y
                   END IF

                   CALL Corner_Forces(x1, y1, x11, y11, p2p_dist_max, P2P_U, P2P_V, Social_F, &
                        Contact_F, P2P_Torque, d_walls, x_tmp, y_tmp, r_tmp, u_tmp, v_tmp, istat)
                   EXIT Loop_mxmy
                END IF
             END DO Loop_mxmy
          END DO Loop_my

          ! ========================================================
          ! THE PERSON-WALL FORCES ENDS HERE
          ! ========================================================

          ! Save the forces for the loop EVAC_MOVE_LOOP (next time step)
          HR%F_X = P2P_U
          HR%F_Y = P2P_V
          HR%SumForces  = Contact_F
          HR%SumForces2 = Social_F + Contact_F
          HR%Torque = P2P_Torque

          IF ( T <= T_BEGIN ) THEN
             IF ( ABS(P2P_U)/HR%Mass > 550.0_EB ) P2P_U =550.0_EB*HR%Mass*P2P_U/ABS(P2P_U)
             IF ( ABS(P2P_V)/HR%Mass > 550.0_EB ) P2P_V =550.0_EB*HR%Mass*P2P_V/ABS(P2P_V)
             HR%F_X = P2P_U
             HR%F_Y = P2P_V
          END IF

          ! Add the social and contact force terms
          U_new     = HR%U + 0.5_EB*HR%F_X*DTSP/HR%Mass
          V_new     = HR%V + 0.5_EB*HR%F_Y*DTSP/HR%Mass
          Omega_new = HR%Omega + 0.5_EB*DTSP*HR%Torque/HR%M_iner

          ! Add the effect of the random force
          IF (GaTh > 0.0_EB .AND. T > T_BEGIN ) THEN
             U_new = U_new + 0.5_EB*DTSP*HR%v0_fac*HR%Mass*HR%ksi*COS(HR%eta)/HR%Mass
             V_new = V_new + 0.5_EB*DTSP*HR%v0_fac*HR%Mass*HR%ksi*SIN(HR%eta)/HR%Mass
             P2P_U = P2P_U + HR%v0_fac*HR%Mass*HR%ksi*COS(HR%eta)
             P2P_V = P2P_V + HR%v0_fac*HR%Mass*HR%ksi*SIN(HR%eta)
             Omega_new = Omega_new + 0.5_EB*DTSP*1.0_EB*SIGN(HR%ksi,HR%eta-Pi)
             P2P_Torque = P2P_Torque + 1.0_EB*SIGN(HR%ksi,HR%eta-Pi)*HR%M_iner
          ELSE
             HR%ksi = 0.0_EB
             HR%eta = 0.0_EB
          END IF
          ! Now the step (4a) of SC-VV by Vattulainen is ended.

          ! Add self-propelling force terms, self-consistent VV 
          ! (first time step towards the exit door)
          fac_tim =  1.0_EB + (DTSP/(2.0_EB*hr_tau))
          IF ( T <= Tpre ) THEN
             IF ( (T+DTSP_new) > Tpre) THEN
                EVEL = SQRT(UBAR**2 + VBAR**2)
                IF (EVEL > 0.0_EB) THEN
                   speed = speed_xp*(0.5_EB + SIGN(0.5_EB,UBAR)) + speed_xm*(0.5_EB - SIGN(0.5_EB,UBAR)) 
                   speed = speed*HR%v0_fac
                   P2P_U = P2P_U + (HR%Mass/hr_tau)*speed*(UBAR/EVEL)
                   speed = speed_yp*(0.5_EB + SIGN(0.5_EB,VBAR)) + speed_ym*(0.5_EB - SIGN(0.5_EB,VBAR)) 
                   speed = speed*HR%v0_fac
                   P2P_V = P2P_V + (HR%Mass/hr_tau)*speed*(VBAR/EVEL)
                END IF
             END IF
             UBAR = 0.0_EB
             VBAR = 0.0_EB
          END IF

          ! Add self-propelling force terms, self-consistent VV
          EVEL = SQRT(UBAR**2 + VBAR**2)
          IF (EVEL > 0.0_EB) THEN
             speed = speed_xp*(0.5_EB + SIGN(0.5_EB,UBAR)) + speed_xm*(0.5_EB - SIGN(0.5_EB,UBAR)) 
             speed = speed*HR%v0_fac
             U_new = (U_new + 0.5_EB*(DTSP/hr_tau)*speed*(UBAR/EVEL)) / fac_tim
             HR%UBAR = UBAR
             P2P_U = P2P_U + (HR%Mass/hr_tau)*(speed*(UBAR/EVEL) - HR%U)

             speed = speed_yp*(0.5_EB + SIGN(0.5_EB,VBAR)) + speed_ym*(0.5_EB - SIGN(0.5_EB,VBAR)) 
             speed = speed*HR%v0_fac
             V_new = (V_new + 0.5_EB*(DTSP/hr_tau)*speed*(VBAR/EVEL)) / fac_tim
             HR%VBAR = VBAR
             P2P_V = P2P_V + (HR%Mass/hr_tau)*(speed*(VBAR/EVEL) - HR%V)

             IF (VBAR >= 0.0_EB) THEN
                angle = ACOS(UBAR/EVEL)
             ELSE
                angle = 2.0_EB*Pi - ACOS(UBAR/EVEL)
             END IF

             ! Collision avoidance has HR%angle_old .ne. 0.0
             angle = angle + HR%angle_old
             DO WHILE (angle >= 2.0_EB*Pi)
                angle = angle - 2.0_EB*Pi
             END DO
             DO WHILE (angle < 0.0_EB)
                angle = angle + 2.0_EB*Pi
             END DO

             IF (angle == 2.0_EB*Pi) angle = 0.0_EB  ! angle is [0,2Pi)

             ! Rotational motion: J(dw/dt) = (J/t_iner)*( ((angle-angle_0)/pi)*w_0 - w)
             IF (ABS(angle-HR%angle) <= Pi ) THEN
                ! zero is not crossed.
                Omega_new = Omega_new + 0.5_EB*(DTSP/HR_tau_iner)*((angle-HR%angle)*(Omega_0/Pi) - HR%Omega)
                P2P_Torque = P2P_Torque + (HR%M_iner/HR_tau_iner)*((angle-HR%angle)*(Omega_0/Pi) - HR%Omega)
             ELSE
                ! zero is crossed
                Omega_new = Omega_new + 0.5_EB*(DTSP/HR_Tau_iner)* &
                     ( (2.0_EB*Pi-ABS(angle-HR%angle))*SIGN(1.0_EB , HR%angle-angle)*(Omega_0/Pi) - HR%Omega)
                P2P_Torque = P2P_Torque + (HR%M_iner/HR_Tau_iner)* &
                     ( (2.0_EB*Pi-ABS(angle-HR%angle))*SIGN(1.0_EB , HR%angle-angle)*(Omega_0/Pi) - HR%Omega)
             END IF
          ELSE  ! No target direction
             HR%UBAR = 0.0_EB
             HR%VBAR = 0.0_EB
             U_new = U_new / fac_tim
             V_new = V_new / fac_tim
             ! Slow rotation down if no direction available, i.e., target Omega_0 is zero.
             Omega_new = Omega_new + 0.5_EB*(DTSP/HR_Tau_iner)*(-HR%Omega)
          END IF

          ! Check that velocities are not too large, i.e., less than 10 m/s for humans
          EVEL = SQRT(U_new**2 + V_new**2)
          IF ( EVEL > Vmax_timo ) THEN
             U_new = U_new*(Vmax_timo/EVEL)
             V_new = V_new*(Vmax_timo/EVEL)
          END IF
          ! Check that angular velocity is not too large
          IF ( ABS(Omega_new) > Omega_max ) THEN
             Omega_new = SIGN(Omega_max,Omega_new)
          END IF

          HR%U = U_new
          HR%V = V_new
          HR%Omega = Omega_new

          ! ========================================================
          ! Decide the time step for the next human movement loop
          ! ========================================================
          ! Distances to closest walls and other agents
          d_humans = MAX(d_humans,0.0005_EB)        ! this agent
          d_walls  = MAX(d_walls, 0.0005_EB)        ! this agent
          d_humans_min = MIN(d_humans_min,d_humans) ! among all agents
          d_walls_min  = MIN(d_walls_min, d_walls)  ! among all agents

          IF ( T > T_BEGIN ) THEN
             ! Time step, do not move too close to other agents, walls, or 0.5*grid spacing.
             ! Delta_min is minimum of dx and dy for the mesh.  Minimum movement is 0.0001 m.
             dt_Loop: DO
                u_tmp(2) = HR%U + 0.5_EB*DTSP_new*P2P_U/HR%Mass
                v_tmp(2) = HR%V + 0.5_EB*DTSP_new*P2P_V/HR%Mass
                Omega_new = HR%Omega + 0.5_EB*DTSP*P2P_Torque/HR%M_iner
                u_tmp(1) = u_tmp(2) + COS(HR%angle)*Omega_new*HR%d_shoulder
                u_tmp(3) = u_tmp(2) - COS(HR%angle)*Omega_new*HR%d_shoulder
                v_tmp(1) = v_tmp(2) + SIN(HR%angle)*Omega_new*HR%d_shoulder
                v_tmp(3) = v_tmp(2) - SIN(HR%angle)*Omega_new*HR%d_shoulder
                IF ( MAX(u_tmp(1)**2+v_tmp(1)**2, u_tmp(2)**2+v_tmp(2)**2, u_tmp(3)**2+v_tmp(3)**2)*DTSP_new**2 > &
                     (MIN(0.2_EB*d_humans, 0.2_EB*d_walls, 0.5_EB*Delta_min))**2 ) THEN
                   DTSP_new = DTSP_new*0.8_EB
                   CYCLE dt_Loop
                END IF
                EXIT dt_Loop
             END DO dt_Loop

             DTSP_new = MAX(DTSP_new, EVAC_DT_MIN)
             DTSP_new = MIN(DTSP_new, EVAC_DT_MAX)
          ELSE  ! Initialization phase
             HR%U     = 0.0_EB
             HR%V     = 0.0_EB
             HR%Omega = 0.0_EB
          END IF

       END DO EVAC_FORCE_LOOP
       TUSED(13,NM)=TUSED(13,NM)+SECOND()-TNOW13
       DEALLOCATE(BLOCK_LIST)
       DEALLOCATE(BLOCK_GRID)
       ! ========================================================
       ! FORCE LOOP ENDS HERE
       ! ========================================================

       ! Save the maximum time steps allowed by the SC-VV algorithm
       Tsteps(NM) = DTSP_new

    END DO HUMAN_TIME_LOOP
    DEALLOCATE(BLOCK_GRID_N)
    ! ========================================================
    ! HUMAN TIME LOOP ENDS HERE (Tin ==> T = Tin + DT)
    ! ========================================================

    ! ========================================================
    ! Evacuation routine ends here
    ! ========================================================
    TUSED(12,NM)=TUSED(12,NM)+SECOND()-TNOW

  CONTAINS

    SUBROUTINE Find_Prefered_Direction(I, Nout, T, T_BEGIN, L_Dead, NM_STRS_MESH, &
         II, JJ, IIX, JJY, XI, YJ, ZK, UBAR, VBAR, hr_tau, Tpre)
      IMPLICIT NONE
      !
      ! Calculate the prefered walking direction
      !
      ! Inputs:
      !   I: The index of the agent (on this mesh)
      !   T: Time
      !   T_BEGIN: The starting time of the simulation
      !   L_Dead: Is the agent dead or not
      !   II,JJ: The grid cell indices of the agent
      !   IIX,JJY: The grid cell indices of the agent for the velocity
      !   XI,YJ,ZK: The grid cell coordinates of the agent for the velocity
      ! Outputs:
      !   Nout: The index of the stairs mesh (if the agent is on stairs)
      !   NM_STRS_MESH: True, if the mesh is a stair mesh
      !   UBAR,VBAR: The prefered walking direction
      !   hr_tau: translational motion Tau
      !   Tpre: Pre-movement time
      !
      ! Passed variables
      INTEGER, INTENT(IN) :: II, JJ, IIX, JJY, I
      REAL(EB), INTENT(IN) :: XI, YJ, ZK, T, T_BEGIN
      LOGICAL, INTENT(IN) :: L_Dead
      INTEGER, INTENT(OUT) :: Nout
      LOGICAL, INTENT(INOUT) :: NM_STRS_MESH
      REAL(EB), INTENT(INOUT) :: hr_tau
      REAL(EB), INTENT(OUT) :: UBAR, VBAR, Tpre
      !
      ! Local variables
      INTEGER :: N,MN_now, STRS_Indx, door_ior, KKZ, j, i1, j1
      REAL(EB) :: x_target, y_target, door_width, Door_Dist, EVEL
      LOGICAL :: NM_STRS_MESHS, straight_Line_To_Target
      TYPE (MESH_TYPE), POINTER :: MFF=>NULL()
      TYPE (HUMAN_TYPE), POINTER :: HR=>NULL()
      TYPE (EVAC_STRS_Type), POINTER :: STRP=>NULL()

      HR=>HUMAN(I)
      KKZ = 1
      ! HR%I_FFIELD is the flow field mesh index for this agent
      NM_now = MAX(1,HR%I_FFIELD)
      IF (.NOT.EVACUATION_ONLY(NM_now)) THEN
         WRITE(lu_err,*)'*** ',i,HR%IPC,HR%IEL,HR%I_Target,HR%I_FFIELD,HR%ILABEL
         WRITE(MESSAGE,'(A,I6,A)') 'ERROR: Evacuate_Humans, mesh ', NM_now, ' is not an evacuation flow field.'
         CALL SHUTDOWN(MESSAGE)
      END IF
      ! 
      ! Determine if the mesh is a stairs-mesh
      NM_STRS_MESH = .FALSE.
      StrsMeshLoop: DO N = 1, N_STRS
         IF (EVAC_STRS(N)%IMESH==NM_now) THEN     
            NM_STRS_MESH = .TRUE.
            STRS_Indx = N
            STRP=>EVAC_STRS(N)
            Nout = N
            EXIT StrsMeshLoop
         END IF
      END DO StrsMeshLoop

      MFF=>MESHES(NM_now)  ! Pointer to the flow field mesh

      ! Check if going straight line to target
      Straight_Line_To_Target = .FALSE.
      IF (NM_STRS_MESH) THEN
         IF (HR%I_Target == 0) THEN
            !               write(*,*) 'Finding target within move loop'
            CALL Find_Target_Node_In_Strs(STRP,HR)
         END IF
         N = HR%I_Target
         IF (N>N_DOORS) THEN
            N = N - N_DOORS
            IF (EVAC_EXITS(N)%STR_SUB_INDX == HR%STR_SUB_INDX) THEN
               Straight_Line_To_Target = .TRUE.
               x_target = (EVAC_EXITS(N)%X1 + EVAC_EXITS(N)%X2)/2.0_EB
               y_target = (EVAC_EXITS(N)%Y1 + EVAC_EXITS(N)%Y2)/2.0_EB
               door_width = EVAC_EXITS(N)%Width
               door_ior   = EVAC_EXITS(N)%IOR
            END IF
         ELSE IF (N>0) THEN
            IF (EVAC_DOORS(N)%STR_SUB_INDX == HR%STR_SUB_INDX) THEN
               Straight_Line_To_Target = .TRUE.
               x_target = (EVAC_DOORS(N)%X1 + EVAC_DOORS(N)%X2)/2.0_EB
               y_target = (EVAC_DOORS(N)%Y1 + EVAC_DOORS(N)%Y2)/2.0_EB
               door_width = EVAC_DOORS(N)%Width
               door_ior   = EVAC_DOORS(N)%IOR
            END IF
         END IF
      END IF
      IF (Straight_Line_To_Target) THEN
         UBAR = x_target-HR%X
         VBAR = y_target-HR%Y
         Door_Dist = SQRT((x_target-HR%X)**2+(y_target-HR%Y)**2)
         IF ( Door_Dist < 0.5_EB*Door_width ) THEN
            SELECT CASE(door_ior)
            CASE(-1,+1)
               UBAR = SIGN(1.0_EB,UBAR)
               VBAR = 0._EB
               HR%SKIP_WALL_FORCE_IOR = NINT(UBAR)
            CASE(-2,+2)
               UBAR = 0._EB
               VBAR = SIGN(1.0_EB,VBAR)
               HR%SKIP_WALL_FORCE_IOR = NINT(VBAR)*2
            END SELECT
         END IF
      ELSE IF (NM_STRS_MESH) THEN 
         CALL STRS_U_AND_V(STRP,HR%STR_SUB_INDX,HR%X,HR%Y,HR%STRS_Direction,UBAR,VBAR)
!         If (HR%STRS_Direction > 0) Then
!            UBAR = STRP%U_UP(II,JJ)
!            VBAR = STRP%V_UP(II,JJ)
!         Else
!            UBAR = STRP%U_DOWN(II,JJ)
!            VBAR = STRP%V_DOWN(II,JJ)
!         End If
!!$      Else If (II < IBAR/4 .And. Abs(JJ-JBAR/2)<3 ) Then 
!!$         UBAR = 0.0_EB
!!$         VBAR = Sign(1.0_EB,HR%VBAR)
!!$      ! Else If (MFF%SOLID(MFF%CELL_INDEX(II,JJ,1))) Then 
!!$         ! There might be additional OBSTs at the door flow fields, which are not
!!$         ! at the main evacuation mesh.  Try to find one grid cell away a field.
!!$         !If (mn_now == HR%IMESH) Then
!!$         !   UBAR = HR%UBAR
!!$         !   VBAR = HR%VBAR
!!$         !Else
!!$         !   i1 = II - Nint(Sign(1.0_EB,HR%UBAR))
!!$         !   j1 = JJ - Nint(Sign(1.0_EB,HR%VBAR))
!!$         !   i1 = Max(1,Min(i1,MFF%IBAR))
!!$         !   j1 = Max(1,Min(j1,MFF%JBAR))
!!$         !   If (.Not. MFF%SOLID(MFF%CELL_INDEX(II,j1,1))) Then
!!$         !      i1 = II
!!$         !   Else If (.Not. MFF%SOLID(MFF%CELL_INDEX(i1,JJ,1))) Then
!!$         !      j1 = JJ
!!$         !   End If
!!$         !   UBAR = MFF%U(i1,j1,1)
!!$         !   VBAR = MFF%V(i1,j1,1)
!!$         !End If
      ELSE
         ! For thin walls and large dx,dy one do not get nice
         ! interpolation by using AFILL2. AFILL2 does not notice that
         ! there are thin walls and, thus, sometimes takes values from
         ! the other side of the thin wall in order to interpolate.
         ! UBAR = AFILL2(MFF%U,II-1,JJY,KKZ,XI-II+1,YJ-JJY+.5_EB,ZK-KKZ+.5_EB)
         ! VBAR = AFILL2(MFF%V,IIX,JJ-1,KKZ,XI-IIX+.5_EB,YJ-JJ+1,ZK-KKZ+.5_EB)
         UBAR = (1.0_EB-(XI-II+1))*MFF%U(II-1,JJ,1) + (XI-II+1)*MFF%U(II,JJ,1)
         VBAR = (1.0_EB-(YJ-JJ+1))*MFF%V(II,JJ-1,1) + (YJ-JJ+1)*MFF%V(II,JJ,1)
      END IF

      EVEL = SQRT(UBAR**2 + VBAR**2)  ! (UBAR,VBAR) is an unit vector
      IF (EVEL > 0.0_EB) THEN
         UBAR = UBAR/EVEL
         VBAR = VBAR/EVEL
      ELSE
         ! No v0 found for the current location of the agent, use previous value
         UBAR = HR%UBAR
         VBAR = HR%VBAR
      END IF
      IF (L_Dead) THEN
         UBAR = 0.0_EB
         VBAR = 0.0_EB
      END IF

      j = MAX(0,HR%GROUP_ID)
      IF (j == 0 ) THEN
         Group_List(0)%Tpre = HR%Tpre + HR%Tdet
         Tpre = HR%Tpre + HR%Tdet
      ELSE
         Tpre = HR%Tdet
      END IF
      IF (L_Dead) Tpre = HUGE(Tpre)

      ! Direction to the centre of the group (if any)
      IF (Group_List(j)%GROUP_SIZE >= 2) THEN
         HR%UBAR_Center = (Group_List(j)%GROUP_X - HR%X)
         HR%VBAR_Center = (Group_List(j)%GROUP_Y - HR%Y)
         EVEL = SQRT(HR%UBAR_Center**2 + HR%VBAR_Center**2)
         IF ( EVEL > 0.0_EB .AND. .NOT. L_Dead ) THEN
            HR%UBAR_Center = HR%UBAR_Center / EVEL
            HR%VBAR_Center = HR%VBAR_Center / EVEL
         ELSE
            HR%UBAR_Center = 0.0_EB
            HR%VBAR_Center = 0.0_EB
         END IF
      ELSE
         HR%UBAR_Center = 0.0_EB ! only one person in the group
         HR%VBAR_Center = 0.0_EB ! only one person in the group
      END IF

      IF ( j > 0 ) THEN
         ! Group is already gathered together, but not yet moving towards the door
         IF (Group_List(j)%COMPLETE == 1 .AND. T <= Group_List(j)%Tpre+Group_List(j)%Tdoor) THEN
            UBAR = 0.0_EB
            VBAR = 0.0_EB
         END IF

         EVEL = UBAR**2 + VBAR**2
         IF (Group_List(j)%COMPLETE == 1 .AND. EVEL > 0.0_EB) THEN
            ! The group is already gathered together
            UBAR = (1-GROUP_EFF)*UBAR + GROUP_EFF*HR%UBAR_Center/ &
                 SQRT(((1-GROUP_EFF)*UBAR + GROUP_EFF*HR%UBAR_Center)**2+ &
                 ((1-GROUP_EFF)*VBAR + GROUP_EFF*HR%VBAR_Center)**2)
            VBAR = (1-GROUP_EFF)*VBAR + GROUP_EFF*HR%VBAR_Center/ &
                 SQRT(((1-GROUP_EFF)*UBAR + GROUP_EFF*HR%UBAR_Center)**2+ &
                 ((1-GROUP_EFF)*VBAR + GROUP_EFF*HR%VBAR_Center)**2)
         ELSE
            ! The group is still in the gathering phase
            UBAR = HR%UBAR_Center
            VBAR = HR%VBAR_Center
         END IF
      END IF

      IF ( T <= Tpre ) THEN
         ! No movement yet
         UBAR = 0.0_EB
         VBAR = 0.0_EB
         hr_tau = MAX(0.1_EB,HR%Tau/10.0_EB)
      END IF
      IF ( T <= T_BEGIN ) THEN
         ! Initialization phase, i.e., flow field calculation
         UBAR = 0.0_EB
         VBAR = 0.0_EB
         hr_tau = HR%Tau
      END IF
      RETURN
    END SUBROUTINE Find_Prefered_Direction


    SUBROUTINE GetStairSpeedAndZ(speed_xm,speed_xp, speed_ym,speed_yp,SP,HP)
      IMPLICIT NONE
      !
      ! Passed variables
      REAL(EB), INTENT(OUT) :: speed_xm, speed_xp, speed_ym, speed_yp
      TYPE (EVAC_STRS_TYPE), POINTER ::  SP
      TYPE (HUMAN_TYPE), POINTER :: HP
      !
      ! Local variables
      REAL(EB) cos_x, cos_y
      INTEGER J1, J2, J

      IF (HP%STR_SUB_INDX > 0) THEN
         J1 = MAX(HP%STR_SUB_INDX-1,1)
         J2 = MIN(HP%STR_SUB_INDX+1,STRP%N_NODES)
      ELSE
         J1 = 1
         J2 = STRP%N_NODES
      END IF

      LoopSubIndx: DO J = J1,J2
         IF (HP%X < SP%XB_NODE(J,1)) CYCLE                
         IF (HP%X > SP%XB_NODE(J,2)) CYCLE                
         IF (HP%Y < SP%XB_NODE(J,3)) CYCLE                
         IF (HP%Y > SP%XB_NODE(J,4)) CYCLE                
         IF (HP%Z < SP%XB_NODE(J,5)) CYCLE                
         IF (HP%Z > SP%XB_NODE(J,6)) CYCLE
         ! Ok, human is inside the subnode J
         HP%STR_SUB_INDX = J      
         IF (SP%NODE_TYPE(J)==STRS_STAIR_TYPE) THEN
            cos_x = SP%XB_NODE(J,7)
            cos_y = SP%XB_NODE(J,8)
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
            !            Select Case (SP%NODE_IOR(J)*HP%STRS_direction)
            SELECT CASE (SP%NODE_IOR(J))
            CASE(-1)
               speed_xm = cos_x* HP%Speed* SP%FAC_V0_UP
               speed_xp = cos_x* HP%Speed* SP%FAC_V0_DOWN
               speed_ym =        HP%Speed* SP%FAC_V0_HORI
               speed_yp =        HP%Speed* SP%FAC_V0_HORI
            CASE(+1)
               speed_xm = cos_x* HP%Speed* SP%FAC_V0_DOWN
               speed_xp = cos_x* HP%Speed* SP%FAC_V0_UP
               speed_ym =        HP%Speed* SP%FAC_V0_HORI
               speed_yp =        HP%Speed* SP%FAC_V0_HORI
            CASE(-2)
               speed_xm =        HP%Speed* SP%FAC_V0_HORI
               speed_xp =        HP%Speed* SP%FAC_V0_HORI
               speed_ym = cos_y* HP%Speed* SP%FAC_V0_UP
               speed_yp = cos_y* HP%Speed* SP%FAC_V0_DOWN
            CASE(+2)
               speed_xm =        HP%Speed* SP%FAC_V0_HORI
               speed_xp =        HP%Speed* SP%FAC_V0_HORI
               speed_ym = cos_y* HP%Speed* SP%FAC_V0_DOWN
               speed_yp = cos_y* HP%Speed* SP%FAC_V0_UP
            END SELECT
         ELSE
            HP%Z = 0.5_EB*(SP%XB_NODE(J,5)+SP%XB_NODE(J,6))      
         END IF
         EXIT LoopSubIndx
      END DO LoopSubIndx
    END SUBROUTINE GetStairSpeedAndZ

    SUBROUTINE Find_Target_Node_In_Strs(SP,HP)
      IMPLICIT NONE
      !
      ! Passed variables
      TYPE (EVAC_STRS_TYPE), POINTER :: SP
      TYPE (HUMAN_TYPE), POINTER :: HP
      !
      ! Local variables
      LOGICAL IsKnownDoor,FinalTargetFound
      INTEGER :: I_Target = 0, I, Id, Final_node, IG, IN, Inode
      REAL(EB) z_node, z_final, dz_node, dz_final, z_final_unknown,dz_tmp1, dz_tmp2, dz_node_actual
      TYPE (EVAC_ENTR_TYPE), POINTER :: PNX =>NULL()

      FinalTargetFound = .FALSE.
      IG = ABS(HP%GROUP_ID)

      z_final_unknown = 0._EB
      dz_node_actual = 0._EB

      DO I = 1, n_nodes
         IsKnownDoor = .FALSE.
         IF (IG == 0 ) THEN
            ! Human came from entry
            PNX => EVAC_ENTRYS(ABS(HR%IEL))
            IF (PNX%N_VENT_FFIELDS == 0) THEN
               IsKnownDoor = .TRUE.
            ELSE
               DO Id = 1, PNX%N_VENT_FFIELDS
                  IF (I == EVAC_Node_List(PNX%I_DOOR_NODES(II))%Node_Index) THEN
                     IsKnownDoor = .TRUE.
                     EXIT
                  END IF
               END DO
            ENDIF
         ELSE            
            ! Check the group know doors
            IF (Human_Known_Doors(IG)%N_nodes == 0) THEN
               IsKnownDoor = .TRUE.
            ELSE
               IN = EVAC_Node_List(I)%Node_index
               DO Id = 1, Human_Known_Doors(IG)%N_nodes
                  IF (I == Human_Known_Doors(IG)%I_nodes(Id)) THEN
                     IsKnownDoor = .TRUE.
                     EXIT
                  END IF
               END DO
            END IF
         END IF
         IF (EVAC_NODE_List(I)%Node_type=='Exit') THEN
            z_final_unknown = EVAC_EXITS(EVAC_NODE_List(I)%Node_index)%Z1
            IF (IsKnownDoor) THEN
               Final_node = I
               z_final = z_final_unknown
               FinalTargetFound = .TRUE.
               EXIT
            END IF
         END IF
      END DO
      IF (.NOT.FinalTargetFound) z_final = z_final_unknown
      dz_final = z_final - HP%Z

      ! Find the target among the nodes leading out of the STRS (doors and exits)
      FinalTargetFound = .FALSE.
      dz_tmp2 = HUGE(1.0_EB)
      FindTargetNodeLoop: DO I = 1,SP%N_NODES_OUT
         Inode = SP%NODES_OUT(I)
         dz_node = -1._EB * SIGN(1.0_EB,dz_final) ! initialize dz_node to different direction than final target
         SELECT CASE(EVAC_NODE_List(Inode)%Node_type)
         CASE('Door')
            z_node = EVAC_DOORS(EVAC_Node_List(Inode)%Node_index)%Z1
         CASE('Exit')
            z_node = EVAC_EXITS(EVAC_Node_List(Inode)%Node_index)%Z1
         CASE default
            WRITE(LU_ERR,*) 'ERROR (debug): unknown node type in Find_Target_Node_In_Strs'
         END SELECT
         dz_node = z_node - HP%Z
         IF (SIGN(1.0_EB,dz_node)==SIGN(1.0_EB,dz_final)) THEN
            dz_tmp1 = ABS(dz_final)-ABS(dz_node)
            IF ( dz_tmp1 < dz_tmp2) THEN
               FinalTargetFound = .TRUE.
               HP%I_Target = EVAC_Node_List(SP%NODES_OUT(I))%Node_index
               dz_tmp2 = dz_tmp1
               dz_node_actual = dz_node
            END IF
         END IF
      END DO FindTargetNodeLoop

      ! Find the target among the nodes leading in to the STRS, i.e. DOORS
      IF (.NOT.FinalTargetFound) THEN
      FindTargetNodeLoop2: DO I = 1,SP%N_NODES_IN
         Inode = SP%NODES_IN(I)
         dz_node = -1._EB * SIGN(1.0_EB,dz_final) ! initialize dz_node to different direction than final target
         SELECT CASE(EVAC_NODE_List(Inode)%Node_type)
         CASE('Door')
            z_node = EVAC_DOORS(EVAC_Node_List(Inode)%Node_index)%Z1
         CASE('Exit')
!            z_node = EVAC_EXITS(EVAC_Node_List(SP%NODES_IN(I))%Node_index)%Z1
         CASE('Entry')
!            z_node = EVAC_ENTRYS(EVAC_Node_list(SP%NODES_IN(I))%Node_index)%Z1
         CASE default
            WRITE(LU_ERR,*) 'ERROR (debug): unknown node type in Find_Target_Node_In_Strs'
         END SELECT
         dz_node = z_node - HP%Z
         IF (SIGN(1.0_EB,dz_node)==SIGN(1.0_EB,dz_final)) THEN
            dz_tmp1 = ABS(dz_final-dz_node)
            IF ( dz_tmp1 < dz_tmp2) THEN
               FinalTargetFound = .TRUE.
               HP%I_Target = EVAC_Node_List(Inode)%Node_index
               dz_tmp2 = dz_tmp1
               dz_node_actual = dz_node
            END IF
         END IF
      END DO FindTargetNodeLoop2
      END IF
      IF (dz_node_actual >= 0._EB) THEN
         HP%STRS_Direction = 1
      ELSE
         HP%STRS_Direction = -1   
      END IF
    END SUBROUTINE Find_Target_Node_In_Strs

    SUBROUTINE STRS_U_AND_V(STRP,I_NODE,X,Y,Direction,UBAR,VBAR)
      IMPLICIT NONE
      !
      ! Get preferred direction in STRS
      !
      ! Passed variables
      TYPE (EVAC_STRS_TYPE), POINTER::  STRP
      INTEGER, INTENT(IN) :: I_NODE, Direction
      REAL(EB), INTENT(IN) :: X, Y
      REAL(EB), INTENT(OUT) :: UBAR, VBAR
      !
      ! Local variables
      INTEGER I_CORE

      I_CORE = STRP%I_CORE(I_NODE)

      UBAR = 0._EB
      VBAR = 0._EB

      IF (STRP%RIGHT_HANDED) THEN
         IF (Direction > 0) THEN
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
         IF (Direction < 0) THEN
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

    SUBROUTINE GET_IW(IIin,JJin,KKin,IOR,IW)
      IMPLICIT NONE
      !
      ! Passed variables
      INTEGER, INTENT(IN) :: IIin, JJin, KKin, IOR
      INTEGER, INTENT(OUT) :: IW
      !
      ! Local variables
      INTEGER :: ii, jj, kk, ic, I_OBST
      !
      ii = IIin
      jj = JJin
      kk = KKin

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
      ! Remove persons if they are found at an exit.
      !
      ! Passed variables
      REAL(EB), INTENT(IN) :: T
      INTEGER, INTENT(IN) :: NM
      !
      ! Local variables
      REAL(EB) x_old, y_old, pexx1, pexx2, pexy1, pexy2
      INTEGER :: ie,i,n_tmp
      TYPE (EVAC_EXIT_Type), POINTER :: PEX =>NULL()
      TYPE (HUMAN_TYPE), POINTER :: HR =>NULL()
      !
      HUMAN(:)%IOR = 0
      PexLoop: DO ie = 1, n_exits
         PEX=>EVAC_EXITS(ie)
         IF (PEX%IMESH /= NM ) CYCLE PexLoop
         SELECT CASE (PEX%IOR)
         CASE (-1,+1)
            pexx1 = PEX%X1
            pexx2 = PEX%X2
            pexy1 = 0.5_EB*(PEX%Y1+PEX%Y2)-0.5_EB*PEX%Width
            pexy2 = 0.5_EB*(PEX%Y1+PEX%Y2)+0.5_EB*PEX%Width
         CASE (-2,+2)
            pexx1 = 0.5_EB*(PEX%X1+PEX%X2)-0.5_EB*PEX%Width
            pexx2 = 0.5_EB*(PEX%X1+PEX%X2)+0.5_EB*PEX%Width
            pexy1 = PEX%Y1
            pexy2 = PEX%Y2
         END SELECT
         HumLoop: DO i = 1, n_humans
            HR=>HUMAN(I)
            IF ( HR%IOR > 0 ) CYCLE HumLoop
            x_old = HR%X_old
            y_old = HR%Y_old
            SELECT CASE (PEX%IOR)
            CASE (+1)
               IF ((HR%X >= pex%x1 .AND. x_old < pex%x1) .AND. (HR%Y >= pex%y1 .AND. HR%Y <= pex%y2) ) THEN
                  IF (PEX%COUNT_ONLY) HR%IOR = 2
                  IF (.NOT. PEX%COUNT_ONLY) HR%IOR = 1
                  PEX%T_last=T
                  PEX%ICOUNT = PEX%ICOUNT + 1
                  IF (PEX%T_first <= T_BEGIN) PEX%T_first = T
               END IF
            CASE (-1)
               IF ((HR%X <= pex%x2 .AND. x_old > pex%x2) .AND. (HR%Y >= pex%y1 .AND. HR%Y <= pex%y2) ) THEN
                  IF (PEX%COUNT_ONLY) HR%IOR = 2
                  IF (.NOT. PEX%COUNT_ONLY) HR%IOR = 1
                  PEX%T_last=T
                  PEX%ICOUNT = PEX%ICOUNT + 1
                  IF (PEX%T_first <= T_BEGIN) PEX%T_first = T
               END IF
            CASE (+2)
               IF ((HR%Y >= pex%y1 .AND. y_old < pex%y1) .AND. (HR%X >= pex%x1 .AND. HR%X <= pex%x2) ) THEN
                  IF (PEX%COUNT_ONLY) HR%IOR = 2
                  IF (.NOT. PEX%COUNT_ONLY) HR%IOR = 1
                  PEX%T_last=T
                  PEX%ICOUNT = PEX%ICOUNT + 1
                  IF (PEX%T_first <= T_BEGIN) PEX%T_first = T
               END IF
            CASE (-2)
               IF ((HR%Y <= pex%y2 .AND. y_old > pex%y2) .AND. (HR%X >= pex%x1 .AND. HR%X <= pex%x2) ) THEN
                  IF (PEX%COUNT_ONLY) HR%IOR = 2
                  IF (.NOT. PEX%COUNT_ONLY) HR%IOR = 1
                  PEX%T_last=T
                  PEX%ICOUNT = PEX%ICOUNT + 1
                  IF (PEX%T_first <= T_BEGIN) PEX%T_first = T
               END IF
            END SELECT
            IF (HR%IOR > 0) THEN
               IF (.NOT. PEX%COUNT_ONLY) WRITE (LU_EVACOUT,fmt='(a,i6,a,f8.2,a,a,a,f8.4,a,i4)') &
                    ' Agent n:o', HR%ILABEL, ' out at ', T, ' s, exit ', TRIM(PEX%ID), &
                    ', FED=', HR%IntDose, ', color_i=', HR%COLOR_INDEX
               IF (      PEX%COUNT_ONLY) WRITE (LU_EVACOUT,fmt='(a,i6,a,f8.2,a,a,a,f8.4,a,i4)') &
                    ' Agent n:o', HR%ILABEL, ' counted at ', T, ' s, counter ', TRIM(PEX%ID), &
                    ', FED=', HR%IntDose, ', color_i=', HR%COLOR_INDEX
               IF (HR%IOR == 2) HR%IOR = HUMAN_NO_TARGET
            END IF
         END DO HumLoop
      END DO PexLoop
      n_tmp = n_humans
      Remove_loop: DO i = n_tmp, 1, -1
         HR=>HUMAN(I)
         ! Remove person, if wanted.
         IF (HR%IOR >= 1) THEN
            HR%IOR = HUMAN_NO_TARGET
            CALL Remove_Person(I)
         END IF
      END DO Remove_loop
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
      REAL(EB) x_old, y_old, xx, yy, zz, pdxx1, pdxx2, pdxy1, pdxy2, v, angle
      INTEGER :: ie,i,n_tmp, istat, ior_new, inode2, imesh2, n, ior
      INTEGER :: new_ffield_i, color_index, i_target, inode, STR_INDX, STR_SUB_INDX
      CHARACTER(60) :: TO_NODE
      CHARACTER(26) :: new_ffield_name
      LOGICAL :: keep_xy, upstream
      TYPE (EVAC_DOOR_Type), POINTER :: PDX =>NULL()
      TYPE (HUMAN_TYPE), POINTER :: HR =>NULL()
      !
      keep_xy = .FALSE.
      HUMAN(:)%IOR = HUMAN_NO_TARGET
      PdxLoop: DO ie = 1, n_doors
         PDX=>EVAC_DOORS(ie)
         ! Note: IMESH2 is not good for CORR targets
         IF (PDX%IMESH /= NM .AND. PDX%IMESH2 /= NM) CYCLE PdxLoop
         IF (PDX%IMESH /= NM .AND. .NOT.NM_STRS_MESH) CYCLE PdxLoop
         keep_xy = PDX%KEEP_XY
         SELECT CASE (PDX%IOR)
         CASE (-1,+1)
            pdxx1 = PDX%X1
            pdxx2 = PDX%X2
            pdxy1 = 0.5_EB*(PDX%Y1+PDX%Y2)-0.5_EB*PDX%Width
            pdxy2 = 0.5_EB*(PDX%Y1+PDX%Y2)+0.5_EB*PDX%Width
         CASE (-2,+2)
            pdxx1 = 0.5_EB*(PDX%X1+PDX%X2)-0.5_EB*PDX%Width
            pdxx2 = 0.5_EB*(PDX%X1+PDX%X2)+0.5_EB*PDX%Width
            pdxy1 = PDX%Y1
            pdxy2 = PDX%Y2
         END SELECT
         HumLoop: DO i = 1, n_humans
            HR=>HUMAN(I)
            new_ffield_name = TRIM(HR%FFIELD_NAME)
            new_ffield_i = HR%I_FFIELD
            ! Check only one door... (Yes, this can happen...)
            IF ( HR%IOR /= HUMAN_NO_TARGET) CYCLE HumLoop
            ! Coming out upstream?
            upstream = .FALSE.
            IF (PDX%IMESH2 == NM .AND. NM_STRS_MESH) THEN
               ! here, check for correct height
               IF ((HR%Z > PDX%Z2) .OR. (HR%Z < PDX%Z1)) CYCLE HumLoop
               upstream = .TRUE.
            END IF
            x_old = HR%X_old
            y_old = HR%Y_old
            SELECT CASE (ABS(PDX%IOR))
            CASE (1)
               IF (SIGN(1.0_EB,HR%X-pdx%x1)/=SIGN(1.0_EB,x_old-pdx%x1) .AND. (HR%Y >= pdx%y1 .AND. HR%Y <= pdx%y2) ) THEN
                  HR%IOR = HUMAN_TARGET_UNSPECIFIED
               END IF
            CASE (2)
               IF (SIGN(1.0_EB,HR%Y-pdx%y1)/=SIGN(1.0_EB,y_old-pdx%y1) .AND. (HR%X >= pdx%x1 .AND. HR%X <= pdx%x2) ) THEN
                  HR%IOR = HUMAN_TARGET_UNSPECIFIED
               END IF
            END SELECT
            IF ( HR%IOR == HUMAN_TARGET_UNSPECIFIED ) THEN
               istat = 0
               IF (upstream) THEN
                  inode  = PDX%INODE2
                  inode2 = PDX%INODE
               ELSE
                  inode  = PDX%INODE
                  inode2 = PDX%INODE2
               END IF
               CALL Check_Target_Node(inode,inode2,istat,xx,yy,zz,ior_new, &
                    imesh2,T,new_ffield_name,new_ffield_i,color_index,&
                    i_target, keep_xy, angle, STR_INDX, STR_SUB_INDX, HR)
               IF (istat == 0 ) THEN
                  ! Put person to a new node, i.e., target node has empty space
                  HR%X = xx
                  HR%Y = yy
                  HR%Z = zz
                  HR%Angle = angle
                  IF (STR_INDX>0) HR%STR_SUB_INDX = STR_SUB_INDX
                  IF (keep_xy .AND. ior_new == HUMAN_TARGET_UNSPECIFIED) THEN
                     IF (EVAC_Node_List(INODE2)%Node_Type == 'Door') THEN
                        xx = 0.5_EB*(EVAC_DOORS(EVAC_Node_List(INODE2)%Node_Index)%X1 + &
                             EVAC_DOORS(EVAC_Node_List(INODE2)%Node_Index)%X2 - (PDX%X1+PDX%X2))
                        yy = 0.5_EB*(EVAC_DOORS(EVAC_Node_List(INODE2)%Node_Index)%Y1 + &
                             EVAC_DOORS(EVAC_Node_List(INODE2)%Node_Index)%Y2 - (PDX%Y1+PDX%Y2))
                     ELSE
                        xx = 0.5_EB*(EVAC_ENTRYS(EVAC_Node_List(INODE2)%Node_Index)%X1 + &
                             EVAC_ENTRYS(EVAC_Node_List(INODE2)%Node_Index)%X2 - (PDX%X1+PDX%X2))
                        yy = 0.5_EB*(EVAC_ENTRYS(EVAC_Node_List(INODE2)%Node_Index)%Y1 + &
                             EVAC_ENTRYS(EVAC_Node_List(INODE2)%Node_Index)%Y2 - (PDX%Y1+PDX%Y2))
                     END IF
                     IF (ABS(PDX%ior) == 2) HR%X = HR%X + xx  
                     IF (ABS(PDX%ior) == 1) HR%Y = HR%Y + yy 
                  END IF
                  v = SQRT( HR%U**2 + HR%V**2 )

                  HR%IOR = ior_new
                  IF (ior_new == HUMAN_TARGET_UNSPECIFIED .OR. ior_new == HUMAN_STRS_TARGET) THEN
                     ! door/entry or STRS-mesh
                     IF (HR%IMESH /= imesh2 ) THEN
                        ! Put the person to a new floor
                        HR%IOR         = HUMAN_ANOTHER_MESH_TARGET
                        HR%INODE = 0
                        DO n = 1, n_egrids
                           IF (EVAC_Node_List(n)%IMESH == imesh2) HR%INODE = n
                        END DO
                        HR%NODE_NAME   = TRIM(MESH_NAME(imesh2))
                        HR%FFIELD_NAME = TRIM(new_ffield_name)
                        HR%I_FFIELD = new_ffield_i
                        HR%SKIP_WALL_FORCE_IOR = 0
                     ELSE
                        HR%IOR = -2   ! same floor (door/entry)
                        HR%FFIELD_NAME = TRIM(new_ffield_name)
                        HR%I_FFIELD = new_ffield_i
                     END IF
                     HR%I_Target = I_Target
                     HR%X_old = HR%X
                     HR%Y_old = HR%Y
                     ! ior is the direction where the human is ejected.
                     IF (EVAC_Node_List(INODE2)%Node_Type == 'Door') THEN
                        ior = - EVAC_DOORS(EVAC_Node_List(INODE2)%Node_Index)%IOR
                     END IF
                     IF (EVAC_Node_List(INODE2)%Node_Type == 'Entry') THEN
                        ior = EVAC_ENTRYS(EVAC_Node_List(INODE2)%Node_Index)%IOR
                     END IF
                     IF (EVAC_Node_List(INODE2)%Node_Type == 'Floor') THEN
                        ! For STRS the first node should be DOOR
                        ior = EVAC_DOORS(EVAC_Node_List(INODE)%Node_Index)%IOR
                     END IF
                     IF ( ABS(ior) == 1 ) THEN
                        HR%U = v*ior
                        HR%V = 0.0_EB
                        IF (.NOT.keep_xy) THEN
                           ! 180 or 0 degrees, i.e., pi or 0 radians
                           HR%Angle = (0.5_EB-(ior/2.0_EB))*Pi
                        END IF
                     END IF
                     IF ( ABS(ior) == 2 ) THEN
                        HR%U = 0.0_EB
                        HR%V = 0.5_EB*v*ior
                        IF (.NOT.keep_xy) THEN
                           ! 270 or 90 degrees, i.e., 3pi/2 or pi/2 radians
                           HR%Angle = (1.0_EB-(ior/4.0_EB))*Pi
                        END IF
                     END IF

                     IF (HR%IMESH /= imesh2 ) THEN
                        HR%IMESH       = imesh2
                        IF (MESHES(imesh2)%N_HUMANS+1 > MESHES(imesh2)%N_HUMANS_DIM) THEN
                           ! Re-allocation is not yet checked.
                           CALL SHUTDOWN('ERROR: Humans: no re-allocation yet')
                           CALL RE_ALLOCATE_HUMANS(1,imesh2)
                        END IF
                        MESHES(imesh2)%N_HUMANS = MESHES(imesh2)%N_HUMANS + 1
                        MESHES(imesh2)%HUMAN(MESHES(imesh2)%N_HUMANS) = HUMAN(I)
                        MESHES(imesh2)%HUMAN(MESHES(imesh2)%N_HUMANS)%IOR = 0
                     END IF

                  END IF            ! target is door or entry

                  PDX%T_last=T
                  PDX%ICOUNT = PDX%ICOUNT + 1
                  IF (PDX%T_first <= T_BEGIN) PDX%T_first = T

                  WRITE (LU_EVACOUT,fmt='(a,i6,a,f8.2,a,a,a,f8.4)') ' EVAC: Person n:o', HR%ILABEL, ' out at ', T, &
                       ' s, door ', TRIM(PDX%ID), ' fed ', HR%IntDose

               ELSE    ! istat = 1 ==> do not move to node
                  ! Can not move to the next node, so do not allow to move inside
                  ! the door ==> keep old position and put velocity to zero.
                  HR%X = HR%X_old
                  HR%Y = HR%Y_old
                  HR%U = 0.0_EB
                  HR%V = 0.0_EB
                  HR%IOR = -1  ! can not move to target node
               END IF
            END IF
         END DO HumLoop
      END DO PdxLoop

      ! Remove humans from this floor and put them some other floor
      n_tmp = n_humans
      Remove_loop: DO i = n_tmp, 1, -1
         HR=>HUMAN(I)
         ! hr%ior = 2: put to an another floor (target is door/entry)
         !          0: not entering a door
         !          3: target is corridor (remove from floor)
         !          4: not used (floor node...)
         !          5: target is exit (remove from floor)
         !         -1: can not move to the target node
         !         -2: move to door/entry on the same floor
         !
         ! Remove person from the floor
         IF (HR%IOR >= 1) THEN
            HR%IOR = HUMAN_NO_TARGET
            CALL Remove_Person(I)
         END IF
      END DO Remove_loop
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
      REAL(EB) x_old, y_old, xx, yy, zz, pcxx1, pcxx2, pcxy1, pcxy2, v, x_int, angle
      INTEGER :: ie,i,n_tmp, istat, ior_new, inode2, imesh2, n, ior
      INTEGER :: new_ffield_i, color_index, i_target, inode, STR_INDX, STR_SUB_INDX
      CHARACTER(60) :: TO_NODE
      CHARACTER(26) :: new_ffield_name
      LOGICAL :: keep_xy
      TYPE (EVAC_CORR_Type),  POINTER :: PCX=>NULL()
      TYPE (CORR_LL_TYPE), POINTER :: Now_LL=>NULL(), Tmp_LL=>NULL(), Next_LL=>NULL(), Prev_LL=>NULL()
      TYPE (HUMAN_TYPE), POINTER :: HR =>NULL()
      !
      keep_xy = .FALSE.
      PcxLoop: DO ie = 1, n_corrs
         PCX=>EVAC_CORRS(ie)
         IF (PCX%IMESH2 /= NM) CYCLE PcxLoop
         IF (PCX%n_inside <= 0) CYCLE PcxLoop
         NULLIFY(Prev_LL)
         NULLIFY(Now_LL)
         NULLIFY(Next_LL)
         IF (ASSOCIATED(PCX%First)) THEN
            Now_LL => PCX%First
            IF (ASSOCIATED(PCX%First%Next)) Next_LL => PCX%First%Next
         END IF
         HumLoop: DO
            IF (.NOT. ASSOCIATED(Now_LL)) EXIT HumLoop
            istat = 0
            inode = PCX%INODE
            inode2 = PCX%INODE2
            HR => Now_LL%HUMAN
            IF ( HR%IntDose >= 1.0_EB  ) THEN
               IF (HR%Tpre /= HUGE(HR%Tpre)) THEN
                  n_dead = n_dead+1
                  HR%Tpre = HUGE(HR%Tpre)
                  HR%Tdet = HUGE(HR%Tdet)
                  HR%Tau  = HR%Tau
                  HR%Mass = HR%Mass
               END IF
            ELSE
               fed_max_alive = MAX(fed_max_alive,HR%IntDose)
            END IF
            fed_max = MAX(fed_max,HR%IntDose)

            ! Calculate Purser's fractional effective dose (FED)
            IF (T > T_BEGIN) THEN
               IF ( PCX%FED_MESH2 > 0 ) THEN
                  x_int = MIN(1.0_EB,MAX(0.0_EB,(Now_LL%T_out-T)) / (Now_LL%T_out - Now_LL%T_in))
               ELSE
                  x_int = 1.0_EB
               END IF
               HR%IntDose = DTSP*( (1.0_EB-x_int)*PCX%FED_CO_CO2_O2(2) + x_int*PCX%FED_CO_CO2_O2(1) ) + HR%IntDose
            END IF

            IF ( (Now_LL%T_out) <= T) THEN
               CALL Check_Target_Node(inode,inode2,istat,xx,yy,zz,ior_new, &
                    imesh2,T,new_ffield_name,new_ffield_i,color_index, &
                    i_target,keep_xy,angle,STR_INDX,STR_SUB_INDX, HR)

               IF (istat == 0 ) THEN
                  ! Put person to a new node, i.e., target node has empty space
                  HR%X = xx
                  HR%Y = yy
                  HR%Z = zz
                  HR%Angle = angle

                  v = SQRT( HR%U**2 + HR%V**2 )
                  IF (v > 0.01_EB) v = PCX%Fac_Speed*HR%v0_fac*HR%Speed
                  HR%IOR = 1        ! remove from corridor
                  ! Put the person to a floor
                  IF (ior_new == 1) THEN ! door or entry
                     HR%IMESH = imesh2
                     HR%INODE = 0
                     DO n = 1, n_egrids
                        IF (EVAC_Node_List(n)%IMESH == imesh2) HR%INODE = n
                     END DO
                     HR%NODE_NAME   = TRIM(MESH_NAME(imesh2))
                     HR%FFIELD_NAME = TRIM(new_ffield_name)
                     HR%I_FFIELD = new_ffield_i
                     HR%I_Target = I_Target

                     IF (EVAC_Node_List(INODE2)%Node_Type == 'Door') THEN
                        ior = - EVAC_DOORS( EVAC_Node_List(INODE2)%Node_Index)%IOR
                     END IF
                     IF (EVAC_Node_List(INODE2)%Node_Type == 'Entry') THEN
                        ior = EVAC_ENTRYS( EVAC_Node_List(INODE2)%Node_Index)%IOR
                     END IF
                     IF ( ABS(ior) == 1 ) THEN
                        HR%U = v*ior
                        HR%V = 0.0_EB
                     END IF
                     IF ( ABS(ior) == 2 ) THEN
                        HR%U = 0.0_EB
                        HR%V = 0.5_EB*v*ior
                     END IF
                     HR%X_old = HR%X
                     HR%Y_old = HR%Y

                     IF (MESHES(imesh2)%N_HUMANS+1 > MESHES(imesh2)%N_HUMANS_DIM) THEN
                        ! Re-allocation is not yet checked.
                        CALL SHUTDOWN('ERROR: Humans: no re-allocation yet')
                        CALL RE_ALLOCATE_HUMANS(1,imesh2)
                     END IF
                     MESHES(imesh2)%N_HUMANS = MESHES(imesh2)%N_HUMANS + 1
                     MESHES(imesh2)%HUMAN(MESHES(imesh2)%N_HUMANS) = HR
                     WRITE (LU_EVACOUT,fmt='(a,i6,a,f8.2,a,a,a,f8.4)') ' EVAC: Person n:o', &
                          HR%ILABEL, ' out at ', T, ' s, corr ', TRIM(PCX%ID), ' fed ', HR%IntDose
                  END IF            ! target is door or entry

                  IF (ior_new == 3) THEN ! corridor
                     WRITE (LU_EVACOUT,fmt='(a,i6,a,f8.2,a,a,f8.4)') ' EVAC: Person n:o', &
                          HR%ILABEL, ' change corr ', T, ' s, corr ', TRIM(PCX%ID), HR%IntDose
                  END IF

                  IF (ior_new == 5) THEN ! exit
                     WRITE (LU_EVACOUT,fmt='(a,i6,a,f8.2,a,a,f8.4)') ' EVAC: Person n:o', &
                          HR%ILABEL, ' exits ', T, ' s, corr ', TRIM(PCX%ID), HR%IntDose
                  END IF
               ELSE
                  ! Can not move to the next node, so do not allow to move inside
                  ! the door ==> keep old position and put velocity to zero.
                  HR%X = HR%X_old
                  HR%Y = HR%Y_old
                  HR%I_Target = 0
                  HR%U = 0.0_EB
                  HR%V = 0.0_EB
                  HR%IOR = HUMAN_NO_TARGET
               END IF         ! istat == 0

            ELSE             ! T_out > T, i.e., still in corridor
               HR%X = HR%X_old
               HR%Y = HR%Y_old
               HR%IOR = HUMAN_NO_TARGET
               HR%I_Target = 0
            END IF

            ! Here person is removed from the corridor (linked list)
            IF (HR%IOR > 0) THEN
               PCX%n_inside = PCX%n_inside - 1
               IF (ASSOCIATED(Now_LL%Next)) THEN
                  Tmp_LL => Now_LL%Next 
                  NULLIFY(HR)       ! HR points to Now_LL%Human
                  NULLIFY(Now_LL%Next) ! remove the pointer
                  DEALLOCATE(Now_LL) ! free memory
                  NULLIFY(Now_LL)   ! remove the pointer
                  Now_LL => Tmp_LL  ! Jump to the next element in LL
                  IF (.NOT. ASSOCIATED(Prev_LL)) THEN
                     PCX%First => Now_LL   ! first element removed
                  ELSE
                     Prev_LL%Next => Now_LL ! deleted element skipped
                  END IF
               ELSE                ! remove the last element
                  NULLIFY(HR)       ! HR points to Now_LL%Human
                  NULLIFY(Now_LL%Next) ! end list
                  DEALLOCATE(Now_LL) ! free memory
                  NULLIFY(Now_LL)   ! end loop
                  IF (ASSOCIATED(Prev_LL)) THEN
                     NULLIFY(Prev_LL%Next) ! end the list
                  ELSE
                     NULLIFY(PCX%First)   ! empty list
                  END IF
               END IF
            ELSE
               Prev_LL => Now_LL  ! save the previous element
               IF (ASSOCIATED(Now_LL%Next)) THEN
                  Now_LL => Now_LL%Next ! Jump to the next element in LL
               ELSE
                  NULLIFY(Now_LL)  ! end the loop
               END IF
            END IF
         END DO HumLoop
      END DO PcxLoop
      !
    END SUBROUTINE CHECK_CORRS
    !
    SUBROUTINE Check_Target_Node(inode,inode2,istat,xx,yy,zz,ior_new, &
         imesh2,T,new_ffield_name,new_ffield_i,color_index,i_target,keep_xy,angle,&
         STR_INDX, STR_SUB_INDX, HR)
      IMPLICIT NONE
      !
      ! Check, that the target node is free.
      ! If target node is door/entry, try to put person to the
      ! floor. 
      !
      ! ior_new = 1: target is a door or entry
      !           3: target is a corridor
      !           4: target is a floor
      !           5: target is an exit
      ! istat   = 1: target is not free (do not move there)
      !           0: target is free
      !
      ! Passed variables
      INTEGER, INTENT(IN) :: inode, inode2
      INTEGER, INTENT(OUT) :: istat, ior_new, imesh2, color_index, i_target
      REAL(EB), INTENT(OUT) :: xx, yy, zz, angle
      INTEGER, INTENT(OUT) :: STR_INDX, STR_SUB_INDX
      REAL(EB), INTENT(IN) :: T
      INTEGER, INTENT(INOUT) :: new_ffield_i
      LOGICAL, INTENT(IN) :: keep_xy
      CHARACTER(26), INTENT(INOUT) :: new_ffield_name
      TYPE (HUMAN_TYPE), POINTER :: HR
      !
      ! Local variables
      REAL(EB) RN, x1, x2, y1, y2, z1, z2, d_max, dist, Width, &
           xx1,yy1, max_fed, ave_K
      INTEGER  II, JJ, KK, ior, irnmax, irn, ie, izero, j1
      REAL(EB), DIMENSION(6) :: r_tmp, x_tmp, y_tmp
      INTEGER :: i_tmp, i_tim, iii, jjj, I_OBST
      LOGICAL :: PP_see_door, keep_xy2, NM_STRS_MESH

      TYPE (CORR_LL_TYPE), POINTER :: TmpCurrent =>NULL(), TmpLoop =>NULL()
      TYPE (EVAC_STRS_TYPE), POINTER :: STRP =>NULL()
      TYPE (EVAC_DOOR_Type), POINTER :: PDX2 =>NULL()
      TYPE (MESH_TYPE), POINTER :: MMF =>NULL()
      TYPE (EVACUATION_Type), POINTER :: HPT =>NULL()
      TYPE (HUMAN_TYPE), POINTER :: HRE =>NULL()
      TYPE (EVAC_ENTR_Type), POINTER :: PNX =>NULL(), PNX2 =>NULL()
      TYPE (EVAC_CORR_Type), POINTER :: PCX2 =>NULL()
      !
      xx = 0.0_EB ; yy = 0.0_EB ; zz = 0.0_EB 
      I_Target = 0
      STR_INDX = 0
      STR_SUB_INDX = 0

      ! Where is the person going to?
      SelectTargetType: SELECT CASE (EVAC_Node_List(INODE2)%Node_Type)
      CASE ('Door','Entry') SelectTargetType
         ior_new = 1
         IF (EVAC_Node_List(INODE2)%Node_Type == 'Door') THEN
            PDX2 => EVAC_DOORS(EVAC_Node_List(INODE2)%Node_Index)
            imesh2  = PDX2%IMESH
            new_ffield_name = TRIM(PDX2%GRID_NAME)
            Irn_Loop1: DO irn = 1, nmeshes
               IF ( evacuation_only(irn) .AND. TRIM(PDX2%GRID_NAME) == TRIM(MESH_NAME(irn)) ) THEN
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
            ior = -PDX2%IOR        ! now entering the room
            Width = PDX2%Width
         ELSE
            PNX2 => EVAC_ENTRYS(EVAC_Node_List(INODE2)%Node_Index) 
            imesh2  = PNX2%IMESH
            new_ffield_name = TRIM(PNX2%GRID_NAME)
            Irn_Loop2: DO irn = 1, nmeshes
               IF ( evacuation_only(irn) .AND. TRIM(PNX2%GRID_NAME) == TRIM(MESH_NAME(irn)) ) THEN
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
         MFF=>Meshes(imesh2)
         !
         irn = 0
         angle = 0.0_EB
         CheckPPForce: DO WHILE (irn < irnmax)
            SELECT CASE (ior)
            CASE(-1,1)
               ! 180 or 0 degrees, i.e., pi or 0 radians
               angle = (0.5_EB-(ior/2.0_EB))*Pi
               CALL RANDOM_NUMBER(rn)
               IF (keep_xy) THEN
                  yy = HR%Y + (irn/MAX(irn,1))*(rn-0.5_EB)*0.25_EB*HR%Radius
                  xx = x1 + ior*(1.0_EB*HR%B + HR%Radius)
                  angle = HR%Angle
               ELSE
                  yy = 0.5_EB*(y1+y2) + (rn-0.5_EB)*MAX(0.0_EB,Width-2.0_EB*HR%Radius-2.0_EB*HR%B)
                  xx = x1 + ior*(1.0_EB*HR%B + HR%r_torso)
               END IF
            CASE(-2,2)
               ! 270 or90 degrees, i.e., 3pi/2 or pi/2 radians
               angle = (1.0_EB-(ior/4.0_EB))*Pi
               CALL RANDOM_NUMBER(rn)
               IF (keep_xy) THEN
                  xx = HR%X + (irn/MAX(irn,1))*(rn-0.5_EB)*0.25_EB*HR%Radius
                  yy = y1 + (ior/ABS(ior))*(1.0_EB*HR%B + HR%Radius)
                  angle = HR%Angle
               ELSE
                  xx = 0.5_EB*(x1+x2) + (rn-0.5_EB)*MAX(0.0_EB,Width-2.0_EB*HR%Radius-2.0_EB*HR%B)
                  yy = y1 + (ior/ABS(ior))*(1.0_EB*HR%B + HR%r_torso)
               END IF
            CASE(0,3)
               CALL RANDOM_NUMBER(rn)
               yy = 0.5_EB*(y1+y2) + (rn-0.5_EB)*MAX(0.0_EB,(y2-y1)-2.0_EB*HR%Radius-2.0_EB*HR%B)
               CALL RANDOM_NUMBER(rn)
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
               IF (STR_SUB_INDX /= HRE%STR_SUB_INDX) CYCLE P2PLoop
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
               HPT => EVACUATION(HR%IEL)
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
               WRITE(MESSAGE,'(A,I6)') 'ERROR: Check_Target_Node, no imesh2 found ',imesh2
               CALL SHUTDOWN(MESSAGE)
            END IF

            IF (j > 0 .AND. Group_List(j)%GROUP_I_FFIELDS(i_tmp) > 0) THEN
               ! There are already group members on this floor, use the same field
               new_ffield_i = Group_List(j)%GROUP_I_FFIELDS(i_tmp)
               new_ffield_name = TRIM(MESH_NAME(new_ffield_i))
               IF (j > 0) I_Target = Group_Known_Doors(j)%I_Target
            ELSE

               CALL Change_Target_Door(imesh2, imesh2, 1, j, j1, 0, 2, HR%X, HR%Y, I_Target, color_index, new_ffield_i, HR)

               new_ffield_name = TRIM(MESH_NAME(new_ffield_i))
               IF ( j > 0 ) THEN
                  Group_List(j)%GROUP_I_FFIELDS(i_tmp) = new_ffield_i
               END IF
            END IF  ! first member of a group?

         END IF ! istat=0, i.e., put human to a new node

      CASE ('Corr') SelectTargetType
         PCX2 => EVAC_CORRS(EVAC_Node_List(INODE2)%Node_Index)
         IF (PCX2%n_inside < PCX2%MAX_HUMANS_INSIDE) THEN
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
            MFF=>Meshes(imesh2)
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
               CALL RANDOM_NUMBER(rn)
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
               CALL RANDOM_NUMBER(rn)
               IF (keep_xy2) THEN
                  xx = HR%X + (irn/MAX(irn,1))*(rn-0.5_EB)*0.25_EB*HR%Radius
                  yy = y1 + (ior/ABS(ior))*(1.0_EB*HR%B + HR%Radius)
                  angle = HR%Angle
               ELSE
                  xx = 0.5_EB*(x1+x2) + (rn-0.5_EB)* MAX(0.0_EB,Width-2.0_EB*HR%Radius-2.0_EB*HR%B)
                  yy = y1 + (ior/ABS(ior))*(1.0_EB*HR%B + HR%r_torso)
               END IF
            CASE(0,3)
               CALL RANDOM_NUMBER(rn)
               yy = 0.5_EB*(y1+y2) + (rn-0.5_EB)* MAX(0.0_EB,(y2-y1)-2.0_EB*HR%Radius-2.0_EB*HR%B)
               CALL RANDOM_NUMBER(rn)
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
    SUBROUTINE REMOVE_OUT_OF_GRIDS(T,NM)
      IMPLICIT NONE
      !
      ! Remove humans that do not lie in any mesh
      !
      ! Passed variables
      INTEGER, INTENT(IN) :: NM
      REAL(EB), INTENT(IN) :: T
      !
      ! Local variables
      INTEGER :: IKILL, I
      !
      IKILL = 0
      DROP_LOOP: DO I=1,N_HUMANS
         !
         HR=>HUMAN(I)
         IF (I > N_HUMANS-IKILL) EXIT DROP_LOOP
         IF (HR%X > XS .AND. HR%X < XF .AND. HR%Y > YS .AND. HR%Y < YF) CYCLE DROP_LOOP
         WRITE(LU_ERR,'(A, 3F8.2, I3)') 'WARNING: Person removed at coord: ', HR%X,HR%Y,HR%Z,HR%SKIP_WALL_FORCE_IOR
         HUMAN(I) = HUMAN(N_HUMANS-IKILL)
         IKILL = IKILL + 1
      END DO DROP_LOOP
      N_HUMANS = N_HUMANS - IKILL
      !
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
      REAL(EB) RN, x1, x2, y1, y2, z1, z2, d_max, dist, xx, yy, zz, xx1, yy1
      INTEGER  II, JJ, KK, ior, irnmax, irn, ie, NR
      REAL(EB), DIMENSION(6) ::y_tmp, x_tmp, r_tmp

      TYPE (EVAC_ENTR_Type), POINTER :: PNX =>NULL()
      TYPE (MESH_TYPE), POINTER :: MFF =>NULL()
      TYPE (EVAC_PERS_Type), POINTER :: PCP =>NULL()
      TYPE (HUMAN_TYPE), POINTER :: HR=>NULL(), HRE =>NULL()
      !
      istat = 1
      PNX => EVAC_ENTRYS(I_entry)
      IF (PNX%IMESH /= NM ) RETURN
      IF (PNX%Flow <= 0.0_EB ) RETURN
      IF (PNX%T_Start > Tin) RETURN
      IF (PNX%T_Stop < Tin) RETURN
      IF (PNX%Max_Humans < 0) THEN
         NR     = -PNX%Max_Humans  
         IF (PNX%ICOUNT >= INT(EVALUATE_RAMP(Tin,0._EB,NR))) RETURN
      ELSE
         IF (PNX%ICOUNT >= PNX%Max_Humans) RETURN
      ENDIF
      MFF => MESHES(NM)
      IF ( (Tin-PNX%T_last) < (1.0_EB/PNX%Flow) ) RETURN
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
         CALL SHUTDOWN('ERROR: Insert Humans: no re-allocation yet')
         CALL RE_ALLOCATE_HUMANS(1,NM)
         HUMAN=>MESHES(NM)%HUMAN
      END IF
      !
      PCP => EVAC_PERSON_CLASSES(PNX%IPC)
      HR  => MESHES(NM)%HUMAN(N_HUMANS+1)
      CALL CLASS_PROPERTIES(HR,PCP)
      HR%Tpre = 0.0_EB
      HR%Tdet = T_BEGIN
      HR%IPC  = PNX%IPC
      HR%IEL  = -I_entry
      HR%GROUP_ID = 0
      HR%I_Target = 0
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
            CALL RANDOM_NUMBER(rn)
            yy = 0.5_EB*(y1+y2) + (rn-0.5_EB)*MAX(0.0_EB,PNX%Width-2.0_EB*HR%Radius-2.0_EB*HR%B)
            xx = x1 + ior*5.0_EB*HR%B
            HR%Angle = (1-ior)*Pi/2.0_EB  ! ior=1: 0,  ior=-1: pi
         CASE(-2,2)
            CALL RANDOM_NUMBER(rn)
            xx = 0.5_EB*(x1+x2) + (rn-0.5_EB)*MAX(0.0_EB,PNX%Width-2.0_EB*HR%Radius-2.0_EB*HR%B)
            yy = y1 + (ior/ABS(ior))*5.0_EB*HR%B
            HR%Angle = Pi/2.0_EB + (2-ior)*Pi/4.0_EB  ! ior=2: (3/2)pi,  ior=-2: pi/2
         CASE(0,3)
            CALL RANDOM_NUMBER(rn)
            yy = 0.5_EB*(y1+y2) + (rn-0.5_EB)*MAX(0.0_EB,(y2-y1)-2.0_EB*HR%Radius-2.0_EB*HR%B)
            CALL RANDOM_NUMBER(rn)
            xx = 0.5_EB*(x1+x2) + (rn-0.5_EB)*MAX(0.0_EB,(x2-x1)-2.0_EB*HR%Radius-2.0_EB*HR%B)
            CALL RANDOM_NUMBER(rn)
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
            CALL SHUTDOWN(MESSAGE)
         END SELECT
         HR%FFIELD_NAME = TRIM(PNX%GRID_NAME)
         HR%I_FFIELD    = 0
         Mesh2Loop: DO i = 1, nmeshes
            IF ( evacuation_only(i) .AND. TRIM(HR%FFIELD_NAME) == TRIM(MESH_NAME(i)) ) THEN
               HR%I_FFIELD = i
               EXIT Mesh2Loop
            END IF
         END DO Mesh2Loop
         IF ( HR%I_FFIELD == 0 ) THEN
            WRITE(MESSAGE,'(A,A,A,A)') 'ERROR: ENTR line ',TRIM(PNX%ID), ' problem with flow field name, ', &
                 TRIM(PNX%GRID_NAME),' not found'
            CALL SHUTDOWN(MESSAGE)
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

    SUBROUTINE Corner_Forces(x1, y1, x11, y11, p2p_dist_max, P2P_U, P2P_V, Social_F, &
         Contact_F, P2P_Torque, d_walls, x_tmp, y_tmp, r_tmp, u_tmp, v_tmp, istat)
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
      REAL(EB), INTENT(INOUT) :: P2P_U, P2P_V, Social_F, Contact_F, P2P_Torque, d_walls
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

      IF ( (u**2+v**2) > 0.0_EB ) THEN
         CosPhiFac = ((x11-X1)*u + (y11-Y1)*v) / (dist*SQRT(u**2+v**2))
         CosPhiFac = LambdaW + 0.5_EB*(1.0_EB-LambdaW)*(1.0_EB+CosPhiFac)
      ELSE
         CosPhiFac = 1.0_EB
      END IF
      P2P_U = P2P_U + (X1-x11)*A_Wall*CosPhiFac*EXP(-(dist-HR%Radius)/B_Wall) / dist
      P2P_V = P2P_V + (Y1-y11)*A_Wall*CosPhiFac*EXP(-(dist-HR%Radius)/B_Wall) / dist
      Social_F = Social_F + ABS(A_Wall*CosPhiFac*EXP(-(dist-HR%Radius)/B_Wall))

      DO iii = 1, 3
         dist = SQRT( (y11-y_tmp(iii))**2 + (x11-x_tmp(iii))**2 )

         ! Next is |vector1|*|vector2|
         evel = SQRT((x11-x_tmp(iii))**2+(y11-y_tmp(iii))**2)*SQRT(u_tmp(iii)**2+v_tmp(iii)**2)
         IF (evel > 0.0_EB) evel = ((x11-x_tmp(iii))*u_tmp(iii) + &
              (y11-y_tmp(iii))*v_tmp(iii)) / evel   ! cos theta (scal_prod/(lenght1*length2)
         IF (evel > 0.01_EB) THEN
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
               Fc_x = Fc_x - k_fric*( (1.0_EB-FricFac)*(r_tmp(iii)-dist)+FricFac ) *(y_tmp(iii)-y11)* &
                    ( (y_tmp(iii)-y11)*u_tmp(iii) - (x_tmp(iii)-x11)*v_tmp(iii) )/dist**2
               Fc_y = Fc_y + k_fric*( (1.0_EB-FricFac)*(r_tmp(iii)-dist)+FricFac ) *(x_tmp(iii)-x11)* &
                    ( (y_tmp(iii)-y11)*u_tmp(iii) - (x_tmp(iii)-x11)*v_tmp(iii) )/dist**2
            P2P_Torque = P2P_Torque + Fc_y*(x11-x_tmp(iii)) - Fc_x*(y11-y_tmp(iii))
            P2P_U = P2P_U + Fc_x
            P2P_V = P2P_V + Fc_y
         END IF
      END DO

    END SUBROUTINE Corner_Forces

    SUBROUTINE Door_Forces(nm, x_tmp, y_tmp, r_tmp, u_tmp, v_tmp, p2p_dist_max, d_xy,&
         P2P_U, P2P_V, Social_F, Contact_F, P2P_Torque, FoundWall_xy)
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
      INTEGER, INTENT(IN) :: nm
      REAL(EB), INTENT(IN) :: p2p_dist_max
      REAL(EB), DIMENSION(4), INTENT(IN) :: d_xy
      LOGICAL, DIMENSION(4), INTENT(IN) :: FoundWall_xy
      REAL(EB), DIMENSION(6), INTENT(IN) :: x_tmp, y_tmp, r_tmp, u_tmp, v_tmp
      REAL(EB), INTENT(INOUT) :: P2P_U, P2P_V, Social_F, Contact_F, P2P_Torque
      !
      ! Local variables
      INTEGER :: is, idir, iin, jjn, istat
      REAL(EB) :: CosPhiFac, dist, dist1, dist2

      
      ! Check if there are doors (vents with vel >0)
      DO ii = 1, N_VENT
         IF (ABS(VENTS(ii)%IOR)>2 .OR. SURFACE(VENTS(ii)%IBC)%VEL<=0) CYCLE
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
                    Contact_F, P2P_Torque, d_walls, x_tmp, y_tmp, r_tmp, u_tmp, v_tmp, istat)
               x11 = VENTS(ii)%x2
               y11 = VENTS(ii)%y2
               CALL Corner_Forces(x1, y1, x11, y11, p2p_dist_max, P2P_U, P2P_V, Social_F, &
                    Contact_F, P2P_Torque, d_walls, x_tmp, y_tmp, r_tmp, u_tmp, v_tmp, istat)
            END IF
         CASE(+1)  ! wall at -x direction
            IF (.NOT.FoundWall_xy(1) .AND. ABS(VENTS(ii)%x1-d_xy(1))<0.01_EB .AND. &
                 (VENTS(ii)%y1<y_tmp(2).AND.VENTS(ii)%y2>y_tmp(2))) THEN
               ! There is a outflow vent here
               x11 = VENTS(ii)%x1
               y11 = VENTS(ii)%y1
               CALL Corner_Forces(x1, y1, x11, y11, p2p_dist_max, P2P_U, P2P_V, Social_F, &
                    Contact_F, P2P_Torque, d_walls, x_tmp, y_tmp, r_tmp, u_tmp, v_tmp, istat)
               x11 = VENTS(ii)%x2
               y11 = VENTS(ii)%y2
               CALL Corner_Forces(x1, y1, x11, y11, p2p_dist_max, P2P_U, P2P_V, Social_F, &
                    Contact_F, P2P_Torque, d_walls, x_tmp, y_tmp, r_tmp, u_tmp, v_tmp, istat)
            END IF
         CASE(-2)  ! wall at +y direction
            IF (.NOT.FoundWall_xy(4) .AND. ABS(VENTS(ii)%y1-d_xy(4))<0.01_EB .AND. &
                 (VENTS(ii)%x1<x_tmp(2).AND.VENTS(ii)%x2>x_tmp(2))) THEN
               ! There is a outflow vent here
               x11 = VENTS(ii)%x1
               y11 = VENTS(ii)%y1
               CALL Corner_Forces(x1, y1, x11, y11, p2p_dist_max, P2P_U, P2P_V, Social_F, &
                    Contact_F, P2P_Torque, d_walls, x_tmp, y_tmp, r_tmp, u_tmp, v_tmp, istat)
               x11 = VENTS(ii)%x2
               y11 = VENTS(ii)%y2
               CALL Corner_Forces(x1, y1, x11, y11, p2p_dist_max, P2P_U, P2P_V, Social_F, &
                    Contact_F, P2P_Torque, d_walls, x_tmp, y_tmp, r_tmp, u_tmp, v_tmp, istat)
            END IF
         CASE(+2)  ! wall at -y direction
            IF (.NOT.FoundWall_xy(3) .AND. ABS(VENTS(ii)%y1-d_xy(3))<0.01_EB .AND. &
                 (VENTS(ii)%x1<x_tmp(2).AND.VENTS(ii)%x2>x_tmp(2))) THEN
               ! There is a outflow vent here
               x11 = VENTS(ii)%x1
               y11 = VENTS(ii)%y1
               CALL Corner_Forces(x1, y1, x11, y11, p2p_dist_max, P2P_U, P2P_V, Social_F, &
                    Contact_F, P2P_Torque, d_walls, x_tmp, y_tmp, r_tmp, u_tmp, v_tmp, istat)
               x11 = VENTS(ii)%x2
               y11 = VENTS(ii)%y2
               CALL Corner_Forces(x1, y1, x11, y11, p2p_dist_max, P2P_U, P2P_V, Social_F, &
                    Contact_F, P2P_Torque, d_walls, x_tmp, y_tmp, r_tmp, u_tmp, v_tmp, istat)
            END IF
         END SELECT
      END DO

    END SUBROUTINE Door_Forces

    SUBROUTINE Wall_SocialForces(nm, x_tmp, y_tmp, r_tmp, p2p_dist_max, d_xy, P2P_U, P2P_V, Social_F, FoundWall_xy)
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
      INTEGER, INTENT(IN) :: nm
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
            IF ( (HR%U**2+HR%V**2) > 0.0_EB ) THEN
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
            IF ( (HR%U**2+HR%V**2) > 0.0_EB ) THEN
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
            IF ( (HR%U**2+HR%V**2) > 0.0_EB ) THEN
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
            IF ( (HR%U**2+HR%V**2) > 0.0_EB ) THEN
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
!!$    Subroutine Wall_SocialForces(nm, x_tmp, y_tmp, r_tmp, p2p_dist_max, d_xy, P2P_U, P2P_V, Social_F, FoundWall_xy)
!!$      Implicit None
!!$
!!$      ! wall - agent social forces
!!$      !
!!$      ! Inputs:  FoundWall_xy(1-4): True for solid walls (False for outflow VENTs)
!!$      !          x/y_tmp: coordinates of the centres of the spheres (1-3 used)
!!$      !                   1: right, 2: centre, 3: left circle
!!$      !          r_tmp: radii of the spheres
!!$      !          p2p_dist_max: cutoff distance of the social force
!!$      ! Outputs: 
!!$      !          P2P_U/V: Forces x/y directions
!!$      !          Social_F: Radial social line force
!!$      !
!!$
!!$      Integer, Intent(IN) :: nm
!!$      Real(EB), Intent(IN) :: x_tmp, y_tmp, r_tmp, p2p_dist_max
!!$      !Real(EB), Dimension(6), Intent(IN) :: x_tmp, y_tmp, r_tmp
!!$      Real(EB), Dimension(4), Intent(IN) :: d_xy
!!$      Logical, Dimension(4), Intent(IN) :: FoundWall_xy
!!$      Real(EB), Intent(INOUT) :: P2P_U, P2P_V, Social_F
!!$      !
!!$      Integer :: is, idir
!!$      Real(EB) :: CosPhiFac, dist
!!$
!!$      ! -x direction
!!$      is   = -1
!!$      idir =  1
!!$      dist = Abs(d_xy(idir) - x_tmp) ! wall - agent centre distance
!!$      If (dist-r_tmp <= P2P_DIST_MAX .And. FoundWall_xy(idir)) Then
!!$         If ( (HR%U**2+HR%V**2) > 0.0_EB ) Then
!!$            CosPhiFac = (is*HR%U)/Sqrt(HR%U**2+HR%V**2)
!!$            CosPhiFac = LambdaW + 0.5_EB*(1.0_EB-LambdaW)*(1.0_EB+CosPhiFac)
!!$         Else
!!$            CosPhiFac = 1.0_EB
!!$         End If
!!$         P2P_U = P2P_U - is*A_Wall*CosPhiFac*Exp( -(dist-r_tmp)/B_Wall )
!!$         Social_F = Social_F + Abs(A_Wall*CosPhiFac*Exp( -(dist-r_tmp)/B_Wall ))
!!$      End If
!!$
!!$      is   = +1
!!$      idir =  2
!!$      dist = Abs(d_xy(idir) - x_tmp) ! wall - agent centre distance
!!$      If (dist-r_tmp <= P2P_DIST_MAX .And. FoundWall_xy(idir)) Then
!!$         If ( (HR%U**2+HR%V**2) > 0.0_EB ) Then
!!$            CosPhiFac = (is*HR%U)/Sqrt(HR%U**2+HR%V**2)
!!$            CosPhiFac = LambdaW + 0.5_EB*(1.0_EB-LambdaW)*(1.0_EB+CosPhiFac)
!!$         Else
!!$            CosPhiFac = 1.0_EB
!!$         End If
!!$         P2P_U = P2P_U - is*A_Wall*CosPhiFac*Exp( -(dist-r_tmp)/B_Wall )
!!$         Social_F = Social_F + Abs(A_Wall*CosPhiFac*Exp( -(dist-r_tmp)/B_Wall ))
!!$      End If
!!$
!!$      is   = -1
!!$      idir =  3
!!$      dist = Abs(d_xy(idir) - y_tmp) ! wall - agent centre distance
!!$      If (dist-r_tmp <= P2P_DIST_MAX .And. FoundWall_xy(idir)) Then
!!$         If ( (HR%U**2+HR%V**2) > 0.0_EB ) Then
!!$            CosPhiFac = (is*HR%V)/Sqrt(HR%U**2+HR%V**2)
!!$            CosPhiFac = LambdaW + 0.5_EB*(1.0_EB-LambdaW)*(1.0_EB+CosPhiFac)
!!$         Else
!!$            CosPhiFac = 1.0_EB
!!$         End If
!!$         P2P_V = P2P_V - is*A_Wall*CosPhiFac*Exp( -(dist-r_tmp)/B_Wall )
!!$         Social_F = Social_F + Abs(A_Wall*CosPhiFac*Exp( -(dist-r_tmp)/B_Wall ))
!!$      End If
!!$
!!$      is   = +1
!!$      idir =  4
!!$      dist = Abs(d_xy(idir) - y_tmp) ! wall - agent centre distance
!!$      If (dist-r_tmp <= P2P_DIST_MAX .And. FoundWall_xy(idir)) Then
!!$         If ( (HR%U**2+HR%V**2) > 0.0_EB ) Then
!!$            CosPhiFac = (is*HR%V)/Sqrt(HR%U**2+HR%V**2)
!!$            CosPhiFac = LambdaW + 0.5_EB*(1.0_EB-LambdaW)*(1.0_EB+CosPhiFac)
!!$         Else
!!$            CosPhiFac = 1.0_EB
!!$         End If
!!$         P2P_V = P2P_V - is*A_Wall*CosPhiFac*Exp( -(dist-r_tmp)/B_Wall )
!!$         Social_F = Social_F + Abs(A_Wall*CosPhiFac*Exp( -(dist-r_tmp)/B_Wall ))
!!$      End If
!!$
!!$    End Subroutine Wall_SocialForces

    SUBROUTINE Wall_ContactForces(nm, x_tmp, y_tmp, r_tmp, u_tmp, v_tmp, d_xy,&
         P2P_U, P2P_V, P2P_Torque, Contact_F, d_walls, FoundWall_xy)
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
      INTEGER, INTENT(IN) :: nm
      REAL(EB), INTENT(IN) :: x_tmp, y_tmp, r_tmp, u_tmp, v_tmp
      REAL(EB), DIMENSION(4), INTENT(IN) :: d_xy
      LOGICAL, DIMENSION(4), INTENT(IN) :: FoundWall_xy
      REAL(EB), INTENT(INOUT) :: P2P_U, P2P_V, P2P_Torque, Contact_F, d_walls
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
         IF (evel > 0.0_EB) evel = ((d_xy(idir)-x_tmp)*u_tmp) / evel
         IF (evel > 0.01_EB) THEN
            d_walls = MIN( (dist-r_tmp)/evel, d_walls)
         ELSE
            d_walls = MIN( (dist-r_tmp)/0.01_EB, d_walls)
         END IF
         ! d_walls = Min(dist-r_tmp,d_walls) ! wall - circle centre distance
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
         !write(lu_evacout,*)'*** -x Fc',Fc_x,Fc_y,x_tmp,y_tmp,u_tmp,v_tmp,icyc
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
         IF (evel > 0.0_EB) evel = ((d_xy(idir)-x_tmp)*u_tmp) / evel
         IF (evel > 0.01_EB) THEN
            d_walls = MIN( (dist-r_tmp)/evel, d_walls)
         ELSE
            d_walls = MIN( (dist-r_tmp)/0.01_EB, d_walls)
         END IF
         ! d_walls = Min(dist-r_tmp,d_walls) ! wall - circle centre distance
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
         !write(lu_evacout,*)'*** +x Fc',Fc_x,Fc_y,x_tmp,y_tmp,u_tmp,v_tmp,icyc
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
         IF (evel > 0.0_EB) evel = ((d_xy(idir)-y_tmp)*v_tmp) / evel
         IF (evel > 0.01_EB) THEN
            d_walls = MIN( (dist-r_tmp)/evel, d_walls)
         ELSE
            d_walls = MIN( (dist-r_tmp)/0.01_EB, d_walls)
         END IF
         !d_walls = Min(dist-r_tmp,d_walls) ! wall - circle centre distance
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
         !write(lu_evacout,*)'*** -y Fc',Fc_x,Fc_y,x_tmp,y_tmp,u_tmp,v_tmp,icyc
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
         IF (evel > 0.0_EB) evel = ((d_xy(idir)-y_tmp)*v_tmp) / evel
         IF (evel > 0.01_EB) THEN
            d_walls = MIN( (dist-r_tmp)/evel, d_walls)
         ELSE
            d_walls = MIN( (dist-r_tmp)/0.01_EB, d_walls)
         END IF
         !d_walls = Min(dist-r_tmp,d_walls) ! wall - circle centre distance
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
         !write(lu_evacout,*)'*** +y Fc',Fc_x,Fc_y,x_tmp,y_tmp,u_tmp,v_tmp,icyc
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
  SUBROUTINE CLASS_PROPERTIES(HR,PCP)
    IMPLICIT NONE
    !
    ! Passed variables
    TYPE (HUMAN_TYPE), POINTER:: HR
    TYPE (EVAC_PERS_Type), POINTER:: PCP
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

    SELECT CASE(PCP%I_VEL_DIST)
    CASE(-1)
       CALL SHUTDOWN('ERROR: Class_Properties: -1')
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
       CALL SHUTDOWN('ERROR: Class_Properties I_VEL_DIST')
    END SELECT

    SELECT CASE(PCP%I_DIA_DIST)
    CASE(-1)
       CALL SHUTDOWN('ERROR: Class_Properties: -1')
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
       CALL SHUTDOWN('ERROR: Class_Properties I_DIA_DIST')
    END SELECT
    HR%Mass   = 80.0_EB*(HR%Radius/0.27_EB)**2 

    SELECT CASE(PCP%I_TAU_DIST)
    CASE(-1)
       CALL SHUTDOWN('ERROR: Class_Properties: -1')
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
       CALL SHUTDOWN('ERROR: Class_Properties I_TAU_DIST')
    END SELECT

    SELECT CASE(PCP%I_DET_DIST)
    CASE(-1)
       CALL SHUTDOWN('ERROR: Class_Properties: -1')
    CASE(0)
       HR%Tdet  = PCP%Tdet_mean
    CASE(1)   ! Uniform
       ! Parameters: (ave,min,max) ave not used
       n_par = 3
       Randomtype = 1
       RandomPara(1) = 0.5_EB*(PCP%Tdet_high+PCP%Tdet_low)
       RandomPara(2) = PCP%Tdet_low
       RandomPara(3) = PCP%Tdet_high
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tdet = rnd_vec(1)
    CASE(2)   ! Truncated Normal
       ! Parameters: (ave,sigma,min,max)
       n_par = 4
       Randomtype = 2
       RandomPara(1) = PCP%Tdet_mean
       RandomPara(2) = PCP%Tdet_para
       RandomPara(3) = PCP%Tdet_low
       RandomPara(4) = PCP%Tdet_high
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tdet = rnd_vec(1)
    CASE(3)   ! Gamma
       ! Parameters: (ave,alpha,beta) ave not used
       n_par = 3
       Randomtype = 3
       RandomPara(1) = PCP%Tdet_mean
       RandomPara(2) = PCP%Tdet_para
       RandomPara(3) = PCP%Tdet_para2
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tdet = rnd_vec(1)
    CASE(4)   ! Normal
       ! Parameters: (ave,sigma)
       n_par = 2
       Randomtype = 4
       RandomPara(1) = PCP%Tdet_mean
       RandomPara(2) = PCP%Tdet_para
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tdet = rnd_vec(1)
    CASE(5)   ! LogNormal
       ! mean and variance of log(x) should be given
       ! Parameters: (ave,sigma) of ln(x)
       n_par = 4
       Randomtype = 5
       RandomPara(1) = PCP%Tdet_mean
       RandomPara(2) = PCP%Tdet_para
       RandomPara(3) = PCP%Tdet_high  ! high end cutoff
       RandomPara(4) = PCP%Tdet_para2 ! shift
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tdet = rnd_vec(1)
    CASE(6)   ! Beta
       ! Parameters: (ave,a,b) ave not used
       n_par = 3
       Randomtype = 6
       RandomPara(1) = PCP%Tdet_mean
       RandomPara(2) = PCP%Tdet_para
       RandomPara(3) = PCP%Tdet_para2
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tdet = rnd_vec(1)
    CASE(7)   ! Triangular
       ! Parameters: (peak,min,max)
       n_par = 3
       Randomtype = 7
       RandomPara(1) = PCP%Tdet_mean
       RandomPara(2) = PCP%Tdet_low
       RandomPara(3) = PCP%Tdet_high
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tdet = rnd_vec(1)
    CASE(8)   ! Weibull  (alpha=1: Exponential)
       ! Parameters: (ave,alpha,lambda)
       n_par = 3
       Randomtype = 8
       RandomPara(1) = PCP%Tdet_mean
       RandomPara(2) = PCP%Tdet_para
       RandomPara(3) = PCP%Tdet_para2
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tdet = rnd_vec(1)
    CASE(9)   ! Gumbel
       ! Parameters: (ave,alpha)
       n_par = 2
       Randomtype = 9
       RandomPara(1) = PCP%Tdet_mean
       RandomPara(2) = PCP%Tdet_para
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tdet = rnd_vec(1)
    CASE Default
       CALL SHUTDOWN('ERROR: Class_Properties I_DET_DIST')
    END SELECT

    SELECT CASE(PCP%I_PRE_DIST)
    CASE(-1)
       CALL SHUTDOWN('ERROR: Class_Properties: -1')
    CASE(0)
       HR%Tpre  = MAX(0._EB,PCP%Tpre_mean)
    CASE(1)   ! Uniform
       ! Parameters: (ave,min,max) ave not used
       n_par = 3
       Randomtype = 1
       RandomPara(1) = 0.5_EB*(PCP%Tpre_high+PCP%Tpre_low)
       RandomPara(2) = PCP%Tpre_low
       RandomPara(3) = PCP%Tpre_high
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tpre = MAX(0._EB,rnd_vec(1))
    CASE(2)   ! Truncated Normal
       ! Parameters: (ave,sigma,min,max)
       n_par = 4
       Randomtype = 2
       RandomPara(1) = PCP%Tpre_mean
       RandomPara(2) = PCP%Tpre_para
       RandomPara(3) = PCP%Tpre_low
       RandomPara(4) = PCP%Tpre_high
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tpre = MAX(0._EB,rnd_vec(1))
    CASE(3)   ! Gamma
       ! Parameters: (ave,alpha,beta) ave not used
       n_par = 3
       Randomtype = 3
       RandomPara(1) = PCP%Tpre_mean
       RandomPara(2) = PCP%Tpre_para
       RandomPara(3) = PCP%Tpre_para2
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tpre = MAX(0._EB,rnd_vec(1))
    CASE(4)   ! Normal
       ! Parameters: (ave,sigma)
       n_par = 2
       Randomtype = 4
       RandomPara(1) = PCP%Tpre_mean
       RandomPara(2) = PCP%Tpre_para
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tpre = MAX(0._EB,rnd_vec(1))
    CASE(5)   ! LogNormal
       ! mean and variance of log(x) should be given
       ! Parameters: (ave,sigma) of ln(x)
       n_par = 4
       Randomtype = 5
       RandomPara(1) = PCP%Tpre_mean
       RandomPara(2) = PCP%Tpre_para
       RandomPara(3) = PCP%Tpre_high  ! high end cutoff
       RandomPara(4) = PCP%Tpre_para2 ! shift
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tpre = MAX(0._EB,rnd_vec(1))
    CASE(6)   ! Beta
       ! Parameters: (ave,a,b) ave not used
       n_par = 3
       Randomtype = 6
       RandomPara(1) = PCP%Tpre_mean
       RandomPara(2) = PCP%Tpre_para
       RandomPara(3) = PCP%Tpre_para2
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tpre = MAX(0._EB,rnd_vec(1))
    CASE(7)   ! Triangular
       ! Parameters: (peak,min,max)
       n_par = 3
       Randomtype = 7
       RandomPara(1) = PCP%Tpre_mean
       RandomPara(2) = PCP%Tpre_low
       RandomPara(3) = PCP%Tpre_high
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tpre = MAX(0._EB,rnd_vec(1))
    CASE(8)   ! Weibull  (alpha=1: Exponential)
       ! Parameters: (ave,alpha,lambda)
       n_par = 3
       Randomtype = 8
       RandomPara(1) = PCP%Tpre_mean
       RandomPara(2) = PCP%Tpre_para
       RandomPara(3) = PCP%Tpre_para2
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tpre = MAX(0._EB,rnd_vec(1))
    CASE(9)   ! Gumbel
       ! Parameters: (ave,alpha)
       n_par = 2
       Randomtype = 9
       RandomPara(1) = PCP%Tpre_mean
       RandomPara(2) = PCP%Tpre_para
       CALL RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tpre = MAX(0._EB,rnd_vec(1))
    CASE Default
       CALL SHUTDOWN('ERROR: Class_Properties I_PRE_DIST')
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
    !
  END SUBROUTINE CLASS_PROPERTIES
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
    IF (.NOT.ANY(EVACUATION_GRID)) RETURN
    IF ( .NOT.(EVACUATION_ONLY(NM) .AND. EVACUATION_GRID(NM)) ) RETURN

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
    REAL(EB) :: TNOW, EVEL, angle_hr
    REAL(FB), ALLOCATABLE, DIMENSION(:) :: XP,YP,ZP
    REAL(FB), ALLOCATABLE, DIMENSION(:,:) :: QP, AP
    INTEGER, ALLOCATABLE, DIMENSION(:) :: TA
    TYPE (HUMAN_TYPE), POINTER :: HR =>NULL()
    !
    TNOW=SECOND() 
    !
    IF (.NOT.ANY(EVACUATION_GRID)) RETURN
    IF (.NOT.(EVACUATION_ONLY(NM) .AND. EVACUATION_GRID(NM))) RETURN
    CALL POINT_TO_MESH(NM)
    !       Call POINT_TO_EVAC_MESH(NM)

    ! Write the current time to the prt5 file, then start
    ! looping through the particle classes
    WRITE(LU_PART(NM)) REAL(T,FB)

    HUMAN_CLASS_LOOP: DO N=1,N_EVAC
       ! Count the number of humans to dump out
       !
       NPLIM = 0
       DO I=1,N_HUMANS
          HR=>HUMAN(I)
          IF (HR%COLOR_INDEX < 1) HR%COLOR_INDEX = 1
          IF (HR%COLOR_INDEX > EVAC_AVATAR_NCOLOR) HR%COLOR_INDEX = EVAC_AVATAR_NCOLOR
          IF (HR%SHOW .AND. N_EVAC == 1) NPLIM = NPLIM + 1
       END DO
       NPLIM = MIN(NPPS,NPLIM)
       !
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
       !
       ! Load human coordinates into single precision array
       !
       NPP = 0
       PLOOP: DO I=1,N_HUMANS
          HR=>HUMAN(I)
          IF (.NOT. HR%SHOW ) CYCLE PLOOP
          NPP = NPP + 1
          TA(NPP) = HR%ILABEL
          XP(NPP) = REAL(HR%X,FB)
          YP(NPP) = REAL(HR%Y,FB)
          ZP(NPP) = MIN( MAX(REAL(HR%Z,FB), EVAC_Z_MIN), EVAC_Z_MAX)

          AP(NPP,1) = 180.0_FB*REAL(HR%Angle/Pi,FB)
          AP(NPP,2) =   2.0_FB*REAL(HR%Radius,FB)
          AP(NPP,3) =   2.0_FB*REAL(HR%r_torso,FB)
          ! Height of a human scaled by radius, default male 1.80 m
          AP(NPP,4) =  1.80_FB*REAL(HR%Radius/0.27_EB,FB)

          DO NN=1,EVAC_N_QUANTITIES
             SELECT CASE(EVAC_QUANTITIES_INDEX(NN))
             CASE(240)  ! MOTIVE_ACCELERATION, Unimpeded walking Speed / tau
                EVEL = SQRT(HR%UBAR**2 + HR%VBAR**2)
                IF (EVEL > 0.0_EB) THEN
                   QP(NPP,NN) = REAL(HR%Speed/HR%Tau,FB)
                ELSE
                   QP(NPP,NN) = REAL(EVEL,FB)
                END IF
             CASE(241)  ! FED_DOSE, Fractional Effective Dose
                QP(NPP,NN) = REAL(HR%IntDose,FB)
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
                IF (EVEL > 0.0_EB) THEN
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
       IF (EVAC_N_QUANTITIES > 0) THEN
          DEALLOCATE(QP)
       END IF
       DEALLOCATE(AP)
       DEALLOCATE(ZP)
       DEALLOCATE(YP)
       DEALLOCATE(XP)
       DEALLOCATE(TA)

    END DO HUMAN_CLASS_LOOP

    !
    TUSED(7,NM) = TUSED(7,NM) + SECOND() - TNOW
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
    REAL(EB) v1,v2,rsq,fac, rn

    IF ( (GaussFlag == 1)  .AND. (ABS(GaussSet2-gmean) <= gcutmult*SQRT(gtheta)) ) THEN
       GaussFlag = 0
       GaussRand = GaussSet2
    ELSE
       GaussFlag = 0
       DO 
          CALL RANDOM_NUMBER(rn)
          v1 = 1.0_EB - 2.0_EB*rn
          CALL RANDOM_NUMBER(rn)
          v2 = 1.0_EB - 2.0_EB*rn
          rsq = v1*v1 + v2*v2
          DO WHILE (rsq >= 1.0 )
             CALL RANDOM_NUMBER(rn)
             v1 = 1.0_EB - 2.0_EB*rn
             CALL RANDOM_NUMBER(rn)
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
    REAL(EB) v1, v2, rsq, fac, rn

    IF ( (GTrunFlag == 1)  .AND. (GTrunSet2 >= glow) .AND. (GTrunSet2 <= ghigh) ) THEN
       GTrunFlag = 0
       GaussTrun = GTrunSet2
    ELSE
       GTrunFlag = 0
       DO 
          CALL RANDOM_NUMBER(rn)
          v1 = 1.0_EB - 2.0_EB*rn
          CALL RANDOM_NUMBER(rn)
          v2 = 1.0_EB - 2.0_EB*rn
          rsq = v1*v1 + v2*v2
          DO WHILE (rsq >= 1.0 )
             CALL RANDOM_NUMBER(rn)
             v1 = 1.0_EB - 2.0_EB*rn
             CALL RANDOM_NUMBER(rn)
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
    ! Dump human data to CHID_evac.csv
    !
    ! Passed variables
    REAL(EB), INTENT(IN) :: Tin
    !
    ! Local variables
    CHARACTER(50) tcform
    INTEGER n_cols, n_tot_humans, i, ii, izero
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ITEMP
    !
    IF (.NOT.ANY(EVACUATION_GRID)) RETURN
    !
    ALLOCATE(ITEMP(MAX(1,n_exits-n_co_exits+n_doors)), STAT = IZERO)
    CALL ChkMemErr('DUMP_EVAC_CSV','ITEMP', IZERO)
    !
    ! Output first the floors then the corridors
    n_cols = n_egrids + n_corrs + n_exits + n_doors + 1 + n_exits - n_co_exits + n_doors
    n_tot_humans = 0
    DO i = 1, n_egrids
       n_tot_humans = n_tot_humans + MESHES(EVAC_Node_List(i)%IMESH)%N_HUMANS
    END DO
    DO i = 1, n_corrs
       n_tot_humans = n_tot_humans + EVAC_CORRS(i)%n_inside
    END DO
    !
    ii = 0
    DO i = 1, n_exits
       IF (.NOT.EVAC_EXITS(i)%COUNT_ONLY) THEN
          ii = ii + 1
          ITEMP(ii) = EVAC_EXITS(i)%NTARGET(50)
       END IF
    END DO
    DO i = 1, n_doors
       ii = ii + 1
       ITEMP(ii) = EVAC_DOORS(i)%NTARGET(50)
    END DO
    !
    IF (n_dead >= 0) THEN
       ! Write the 'fed' columns
       WRITE(tcform,'(a,i4.4,a,a)') "(ES13.5E3,",n_cols+1, &
            "(',',i8)", ",',',ES13.5E3,',',ES13.5E3)"
       WRITE (LU_EVACCSV,fmt=tcform) Tin, n_tot_humans, &
            (MESHES(EVAC_Node_List(i)%IMESH)%N_HUMANS, i=1,n_egrids), &
            (EVAC_CORRS(i)%n_inside, i = 1,n_corrs), &
            (EVAC_EXITS(i)%ICOUNT, i = 1,n_exits), &
            (EVAC_DOORS(i)%ICOUNT, i = 1,n_doors), &
            (ITEMP(i), i = 1,n_exits-n_co_exits+n_doors), &
            n_dead, fed_max, fed_max_alive
    ELSE
       ! Do not write the 'fed' columns
       WRITE(tcform,'(a,i4.4,a)') "(ES13.5E3,",n_cols, "(',',i8),i8)"
       WRITE (LU_EVACCSV,fmt=tcform) Tin, n_tot_humans, &
            (MESHES(EVAC_Node_List(i)%IMESH)%N_HUMANS, i=1,n_egrids), &
            (EVAC_CORRS(i)%n_inside, i = 1,n_corrs), &
            (EVAC_EXITS(i)%ICOUNT, i = 1,n_exits), &
            (EVAC_DOORS(i)%ICOUNT, i = 1,n_doors), &
            (ITEMP(i), i = 1,n_exits-n_co_exits+n_doors)
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
    IF (ABS(i_r2-i_r1)+ABS(j_r2-j_r1) .LT. 1) RETURN

    ! Choose the main direction:
    IF (ABS(i_r2-i_r1) .LT. ABS(j_r2-j_r1)) THEN
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
          IF (iw >0) THEN
             IF (M%OBSTRUCTION(M%OBST_INDEX_W(iw ))%HIDDEN .AND. M%OBST_INDEX_W(iw )>0) iw  = 0
          END IF
          IF (iw1>0) THEN
             IF (M%OBSTRUCTION(M%OBST_INDEX_W(iw1))%HIDDEN .AND. M%OBST_INDEX_W(iw1)>0) iw1 = 0
          END IF
          IF (iw2>0) THEN
             IF (M%OBSTRUCTION(M%OBST_INDEX_W(iw2))%HIDDEN .AND. M%OBST_INDEX_W(iw2)>0) iw2 = 0
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
          IF (iw >0) THEN
             IF (M%OBSTRUCTION(M%OBST_INDEX_W(iw ))%HIDDEN .AND. M%OBST_INDEX_W(iw )>0) iw  = 0
          END IF
          IF (iw1>0) THEN
             IF (M%OBSTRUCTION(M%OBST_INDEX_W(iw1))%HIDDEN .AND. M%OBST_INDEX_W(iw1)>0) iw1 = 0
          END IF
          IF (iw2>0) THEN
             IF (M%OBSTRUCTION(M%OBST_INDEX_W(iw2))%HIDDEN .AND. M%OBST_INDEX_W(iw2)>0) iw2 = 0
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
  LOGICAL FUNCTION See_door(nm, idoor, itarget, r1_x, r1_y, r2_x, r2_y, ave_K, max_fed)
    ! This function returns true, if the two points have a line-of-sight.
    ! This function does use smoke information, i.e., it sees if
    ! there are obstacles and/or too much smoke between the two points.
    ! NOTE: This works only for thick OBSTs (at least one grid cell thick)
    ! Inputs:  nm: mesh index, r1 an r2 should belong to the same mesh
    !          idoor: index of the  door/exit
    !          itarget: the current target door of the agent
    !          (r1_x,r1_y): co-ordinates of the agent
    !          (r2_x,r2_y): co-ordinates of the door
    ! Outputs: ave_K: average extinction coefficient of the path
    !          max_fed: maximum level of FED at the path
    !
    ! Passed variables
    INTEGER, INTENT(IN) :: nm, idoor, itarget
    REAL(EB), INTENT(IN) :: r1_x, r1_y, r2_x, r2_y
    REAL(EB), INTENT(OUT) :: ave_K, max_fed
    !
    ! Local variables
    INTEGER :: i, j, isx, isy, i_r1, i_r2, j_r1, j_r2
    INTEGER :: i_old, j_old, ic, ic2, iw, iw1, iw2
    REAL(EB) :: x, y
    TYPE (MESH_TYPE), POINTER :: M =>NULL()

    M => MESHES(NM)
    See_door = .TRUE.  ! Default
    ave_K = 0.0_EB
    max_fed = 0.0_EB

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
    IF (ABS(i_r2-i_r1)+ABS(j_r2-j_r1) .LT. 1) THEN
       ave_K = MASS_EXTINCTION_COEFFICIENT*1.0E-6_EB*M%HUMAN_GRID(i_r1,j_r1)%SOOT_DENS
       max_fed = M%HUMAN_GRID(i_r1,j_r1)%FED_CO_CO2_O2
       RETURN
    END IF

    ! Choose the main direction:
    IF (ABS(i_r2-i_r1) .LT. ABS(j_r2-j_r1)) THEN
       ! Now y is the main direction

       ! Add the last cell (r1) to the average
       ave_K = Mass_extinction_coefficient*1.0E-6_EB*M%HUMAN_GRID(i_r1,j_r1)%SOOT_DENS/(ABS(j_r1-j_r2)+1)
       max_fed =  M%HUMAN_GRID(i_r1,j_r1)%FED_CO_CO2_O2

       i_old = i_r1 ; j_old = j_r1
       y = 0.0_EB
       MainLoopY: DO j = j_r1+isy, j_r2, isy
          y = y + isy*M%DY(j)
          x = MAX(M%XS,MIN(M%XF,r1_x + y*(r2_x - r1_x)/(r2_y - r1_y)))
          i = FLOOR(M%CELLSI(FLOOR((x-M%XS)*M%RDXINT)) + 1.0_EB)
          i = isx*MIN(isx*i_r2,isx*i) ! i in interval j_r1...j_r2
          ic  = M%CELL_INDEX(i_old,j_old,1)
          ic2 = M%CELL_INDEX(i    ,j_old,1) ! side cell
          iw  = M%WALL_INDEX(ic, isy*2) ! main direction
          iw1 = M%WALL_INDEX(ic ,isx*1) ! sideways
          iw2 = M%WALL_INDEX(ic2,isy*2) ! side + main direction
          IF (iw >0) THEN
             IF (M%OBSTRUCTION(M%OBST_INDEX_W(iw ))%HIDDEN .AND. M%OBST_INDEX_W(iw )>0) iw  = 0
          END IF
          IF (iw1>0) THEN
             IF (M%OBSTRUCTION(M%OBST_INDEX_W(iw1))%HIDDEN .AND. M%OBST_INDEX_W(iw1)>0) iw1 = 0
          END IF
          IF (iw2>0) THEN
             IF (M%OBSTRUCTION(M%OBST_INDEX_W(iw2))%HIDDEN .AND. M%OBST_INDEX_W(iw2)>0) iw2 = 0
          END IF
          ! iw is zero, if there is no solid boundary
          ! from (i,j)==>(i,jnew):    iw and iw2 are zero, iw1 does not matter
          !                           ic=ic2 ==> iw=iw2
          ! from (i,j)==>(inew,jnew): iw1 and iw2 are zero, iw does not matter
          IF ((i==i_old .AND. iw/=0) .OR. (i/=i_old .AND. (iw1/=0 .OR. iw2/=0)) .AND. (idoor .NE. itarget)) THEN
             See_door = .FALSE.
             ave_K = MASS_EXTINCTION_COEFFICIENT*1.0E-6_EB*M%HUMAN_GRID(i_r1,j_r1)%SOOT_DENS
             max_fed = M%HUMAN_GRID(i_r1,j_r1)%FED_CO_CO2_O2
             EXIT MainLoopY
          END IF
          ave_K = ave_K + MASS_EXTINCTION_COEFFICIENT*1.0E-6_EB*M%HUMAN_GRID(i,j)%SOOT_DENS/(ABS(j_r1-j_r2)+1)
          max_fed = MAX(max_fed, M%HUMAN_GRID(i,j)%FED_CO_CO2_O2)
          i_old = i ; j_old = j
       END DO MainLoopY

    ELSE
       ! Now x is the main direction
       ! Add the first cell (r1) to the average
       ave_K = Mass_extinction_coefficient*1.0E-6_EB*M%HUMAN_GRID(i_r1,j_r1)%SOOT_DENS/(ABS(i_r1-i_r2)+1)
       max_fed =  M%HUMAN_GRID(i_r1,j_r1)%FED_CO_CO2_O2
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
          IF (iw >0) THEN
             IF (M%OBSTRUCTION(M%OBST_INDEX_W(iw ))%HIDDEN .AND. M%OBST_INDEX_W(iw )>0) iw  = 0
          END IF
          IF (iw1>0) THEN
             IF (M%OBSTRUCTION(M%OBST_INDEX_W(iw1))%HIDDEN .AND. M%OBST_INDEX_W(iw1)>0) iw1 = 0
          END IF
          IF (iw2>0) THEN
             IF (M%OBSTRUCTION(M%OBST_INDEX_W(iw2))%HIDDEN .AND. M%OBST_INDEX_W(iw2)>0) iw2 = 0
          END IF
          ! iw is zero, if there is no solid boundary
          ! from (i,j)==>(inew,j):    iw and iw2 are zero, iw1 does not matter
          !                           ic=ic2 ==> iw=iw2
          ! from (i,j)==>(inew,jnew): iw1 and iw2 are zero, iw does not matter
          IF ((i==i_old .AND. iw/=0) .OR. (i/=i_old .AND. (iw1/=0 .OR. iw2/=0)) .AND. (idoor .NE. itarget)) THEN
             See_door = .FALSE.
             ave_K = MASS_EXTINCTION_COEFFICIENT*1.0E-6_EB*M%HUMAN_GRID(i_r1,j_r1)%SOOT_DENS
             max_fed = M%HUMAN_GRID(i_r1,j_r1)%FED_CO_CO2_O2
             EXIT MainLoopX
          END IF
          ave_K = ave_K + MASS_EXTINCTION_COEFFICIENT*1.0E-6_EB*M%HUMAN_GRID(i,j)%SOOT_DENS/(ABS(i_r1-i_r2)+1)
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
    INTEGER :: ii, jj, iw, ic, ibc, is, i_end, iin, jjn, kkn, I_OBST
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
       IF (iw>0) THEN
          IF (M%OBSTRUCTION(M%OBST_INDEX_W(iw))%HIDDEN .AND. M%OBST_INDEX_W(iw)>0) iw = 0
       END IF
       DO WHILE (iw==0 .AND. ii/=i_end)
          ii = ii + is
          ic = M%cell_index(ii,jjn,kkn)  ! cell index
          iw = M%wall_index(ic, is*1)      ! wall index
          IF (iw>0) THEN
             IF (M%OBSTRUCTION(M%OBST_INDEX_W(iw))%HIDDEN .AND. M%OBST_INDEX_W(iw)>0) iw = 0
          END IF
       END DO

       IF (iw /= 0) THEN
          FoundWall_xy(1) = .TRUE.
          ibc = M%IJKW(5,iw)         ! Boundary condition index
          ! There is a 'door', i.e., outflow-boundary (or open boundary)
          ! so no wall forces ==> exit this loop
          d_mx = M%xw(iw)
          I_OBST = M%OBST_INDEX_C(IC)
          IF (M%Solid(ic) .AND. .NOT.M%OBSTRUCTION(I_OBST)%HIDDEN) THEN
             WRITE(MESSAGE,'(A,I4,2I6)') 'ERROR: Find_Walls ',nm, ii,jjn
             CALL SHUTDOWN(MESSAGE)
          END IF
          IF (SURFACE(ibc)%VEL> 0.0_EB .OR. M%BOUNDARY_TYPE(iw)==OPEN_BOUNDARY) THEN
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
       IF (iw>0) THEN
          IF (M%OBSTRUCTION(M%OBST_INDEX_W(iw))%HIDDEN .AND. M%OBST_INDEX_W(iw)>0) iw = 0
       END IF
       DO WHILE (iw==0 .AND. ii/=i_end)
          ii = ii + is
          ic = M%cell_index(ii,jjn,kkn)  ! cell index
          iw = M%wall_index(ic, is*1)      ! wall index
          IF (iw>0) THEN
             IF (M%OBSTRUCTION(M%OBST_INDEX_W(iw))%HIDDEN .AND. M%OBST_INDEX_W(iw)>0) iw = 0
          END IF
       END DO

       IF (iw /= 0) THEN
          FoundWall_xy(2) = .TRUE.
          ibc = M%IJKW(5,iw)         ! Boundary condition index
          ! There is a 'door', i.e., outflow-boundary (or open boundary)
          ! so no wall forces ==> exit this loop
          d_px = M%xw(iw)
          I_OBST = M%OBST_INDEX_C(IC)
          IF (M%Solid(ic) .AND. .NOT.M%OBSTRUCTION(I_OBST)%HIDDEN) THEN
             WRITE(MESSAGE,'(A,I4,2I6)') 'ERROR: Find_Walls ',nm, ii,jjn
             CALL SHUTDOWN(MESSAGE)
          END IF
          IF (SURFACE(ibc)%VEL> 0.0_EB .OR. M%BOUNDARY_TYPE(iw)==OPEN_BOUNDARY) THEN
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
       IF (iw>0) THEN
          IF (M%OBSTRUCTION(M%OBST_INDEX_W(iw))%HIDDEN .AND. M%OBST_INDEX_W(iw)>0) iw = 0
       END IF
       DO WHILE (iw==0 .AND. jj/=i_end)
          jj = jj + is
          ic = M%cell_index(iin,jj,kkn)  ! cell index
          iw = M%wall_index(ic, is*2)      ! wall index
          IF (iw>0) THEN
             IF (M%OBSTRUCTION(M%OBST_INDEX_W(iw))%HIDDEN .AND. M%OBST_INDEX_W(iw)>0) iw = 0
          END IF
       END DO

       IF (iw /= 0) THEN
          FoundWall_xy(3) = .TRUE.
          ibc = M%IJKW(5,iw)         ! Boundary condition index
          ! There is a 'door', i.e., outflow-boundary (or open boundary)
          ! so no wall forces ==> exit this loop
          d_my = M%yw(iw)
          I_OBST = M%OBST_INDEX_C(IC)
          IF (M%Solid(ic) .AND. .NOT.M%OBSTRUCTION(I_OBST)%HIDDEN) THEN
             WRITE(MESSAGE,'(A,I4,2I6)') 'ERROR: Find_Walls ',nm, ii,jjn
             CALL SHUTDOWN(MESSAGE)
          END IF
          IF (SURFACE(ibc)%VEL> 0.0_EB .OR. M%BOUNDARY_TYPE(iw)==OPEN_BOUNDARY) THEN
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
       IF (iw>0) THEN
          IF (M%OBSTRUCTION(M%OBST_INDEX_W(iw))%HIDDEN .AND. M%OBST_INDEX_W(iw)>0) iw = 0
       END IF
       DO WHILE (iw==0 .AND. jj/=i_end)
          jj = jj + is
          ic = M%cell_index(iin,jj,kkn)  ! cell index
          iw = M%wall_index(ic, is*2)      ! wall index
          IF (iw>0) THEN
             IF (M%OBSTRUCTION(M%OBST_INDEX_W(iw))%HIDDEN .AND. M%OBST_INDEX_W(iw)>0) iw = 0
          END IF
       END DO

       IF (iw /= 0) THEN
          FoundWall_xy(4) = .TRUE.
          ibc = M%IJKW(5,iw)         ! Boundary condition index
          ! There is a 'door', i.e., outflow-boundary (or open boundary)
          ! so no wall forces ==> exit this loop
          d_py = M%yw(iw)
          I_OBST = M%OBST_INDEX_C(IC)
          IF (M%Solid(ic) .AND. .NOT.M%OBSTRUCTION(I_OBST)%HIDDEN) THEN
             WRITE(MESSAGE,'(A,I4,2I6)') 'ERROR: Find_Walls ',nm, ii,jjn
             CALL SHUTDOWN(MESSAGE)
          END IF
          IF (SURFACE(ibc)%VEL> 0.0_EB .OR. M%BOUNDARY_TYPE(iw)==OPEN_BOUNDARY) THEN
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

  SUBROUTINE GET_FIRE_CONDITIONS(NOM,I,J,K,fed_indx,soot_dens,gas_temp,rad_flux, YY_GET)
    IMPLICIT NONE
    !
    ! Passed variables
    INTEGER, INTENT(IN) :: I, J, K, NOM
    REAL(EB), INTENT(OUT) :: fed_indx, soot_dens, gas_temp, rad_flux
    REAL(EB), INTENT(INOUT) :: YY_GET(1:N_SPECIES)
    !
    ! Local variables
    REAL(EB) :: Y_MF_INT

    ! Mass fraction array ==> soot density (mg/m3)
    ! Next is for soot (mg/m3)
    YY_GET(:) = MESHES(nom)%YY(I,J,K,:)
    CALL GET_MASS_FRACTION(YY_GET,SOOT_INDEX,Y_MF_INT)
    soot_dens = Y_MF_INT*MESHES(nom)%RHO(I,J,K)*1.E6_EB
    ! Calculate Purser's fractional effective dose (FED)
    fed_indx = FED(YY_GET,MESHES(nom)%RSUM(I,J,K))
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
    INTEGER, INTENT(OUT)  :: I_Target, I_Color, I_Field
    TYPE (HUMAN_TYPE), POINTER :: HR
    !
    ! Local variables
    REAL(EB) :: L2_min, max_fed, ave_K, L2_tmp, rn
    REAL(EB) :: x1_old, y1_old, Speed, X11, Y11, x_o, y_o
    INTEGER :: i_old_ffield, i_tmp, i_new_ffield, IEL, color_index
    INTEGER :: i, i_o, izero, nm_tmp
    CHARACTER(26) :: name_old_ffield, name_new_ffield
    LOGICAL :: PP_see_door
    REAL(EB) :: T_tmp, T_tmp1, Width
    INTEGER :: N_queue, ii
    TYPE (EVACUATION_Type), POINTER :: HPT =>NULL()
    TYPE (EVAC_ENTR_Type),  POINTER :: PNX =>NULL()


    IF (imode < 2) THEN
       nm_tmp = nm
    ELSE
       ! This is for check_target_node
       nm_tmp = nm2
    END IF

    i_old_ffield = HR%I_FFIELD
    IF (i_old_ffield > 0) THEN
       name_old_ffield = TRIM(MESH_NAME(i_old_ffield))
    ELSE
       name_old_ffield = TRIM(MESH_NAME(nm_tmp))
    END IF
    IF (HR%IEL > 0 ) THEN
       ! Agent HR originates from an evac line
       HPT => EVACUATION(HR%IEL)
    ELSE
       ! Agent HR originates from an entr line (and is a lonely agent)
       PNX => EVAC_ENTRYS(ABS(HR%IEL))
    END IF

    ! Only those doors are possible which are in the same main evac mesh.
    K_ave_Door(:)      = 0.0_EB
    FED_max_Door(:)    = 0.0_EB
    Is_Known_Door(:)   = .FALSE.
    Is_Visible_Door(:) = .FALSE.
    ! Group index=0: the agent is from an entry line (no evac line)
    IF ( HR%GROUP_ID /= 0 ) THEN
       IF ( HR%GROUP_ID < 0 ) THEN
          ! A lonely soul
          x1_old = xx
          y1_old = yy
          IEL    = HR%IEL
          Speed  = HR%Speed
          DO i = 1, n_doors
             IF ( EVAC_DOORS(i)%IMESH == nm_tmp) THEN
                Is_Visible_Door(i) = .TRUE.
                IF (imode == 0) CYCLE   ! Initialization call
                DO i_tmp = 1, Human_Known_Doors(j1)%N_nodes
                   IF (EVAC_DOORS(i)%INODE == Human_Known_Doors(j1)%I_nodes(i_tmp)) Is_Known_Door(i) = .TRUE.
                END DO
             END IF
          END DO
          DO i = 1, n_exits
             IF ( EVAC_EXITS(i)%IMESH == nm_tmp .AND. .NOT. EVAC_EXITS(i)%COUNT_ONLY ) THEN
                Is_Visible_Door(n_doors+i) = .TRUE.
                IF (imode == 0) CYCLE   ! Initialization call
                DO i_tmp = 1, Human_Known_Doors(j1)%N_nodes
                   IF (EVAC_EXITS(i)%INODE == Human_Known_Doors(j1)%I_nodes(i_tmp)) Is_Known_Door(n_doors+i) = .TRUE.
                END DO
             END IF
          END DO
       ELSE
          ! A member of a group
          x1_old = Group_List(j)%GROUP_X
          y1_old = Group_List(j)%GROUP_Y
          IEL    = Group_List(j)%IEL
          Speed  = Group_List(j)%Speed
          DO i = 1, n_doors
             IF ( EVAC_DOORS(i)%IMESH == nm_tmp) THEN
                Is_Visible_Door(i) = .TRUE.
                IF (imode == 0) CYCLE   ! Initialization call
                DO i_tmp = 1, Group_Known_Doors(j)%N_nodes
                   IF (EVAC_DOORS(i)%INODE == Group_Known_Doors(j)%I_nodes(i_tmp)) Is_Known_Door(i) = .TRUE.
                END DO
             END IF
          END DO
          DO i = 1, n_exits
             IF ( EVAC_EXITS(i)%IMESH == nm_tmp .AND. .NOT. EVAC_EXITS(i)%COUNT_ONLY ) THEN
                Is_Visible_Door(n_doors+i) = .TRUE.
                IF (imode == 0) CYCLE   ! Initialization call
                DO i_tmp = 1, Group_Known_Doors(j)%N_nodes
                   IF (EVAC_EXITS(i)%INODE == Group_Known_Doors(j)%I_nodes(i_tmp)) Is_Known_Door(n_doors+i) = .TRUE.
                END DO
             END IF
          END DO
       END IF

       IF (imode == 0) THEN
          ! Initialization phase: Draw the known doors/exits.
          DO i = 1, HPT%N_VENT_FFIELDS 
             i_tmp = 1
             IF (TRIM(EVAC_Node_List(HPT%I_DOOR_NODES(i))%Node_Type) == 'Door') THEN
                i_tmp = EVAC_Node_List(HPT%I_DOOR_NODES(i))%Node_Index 
             END IF
             IF (TRIM(EVAC_Node_List(HPT%I_DOOR_NODES(i))%Node_Type) == 'Exit' ) THEN
                i_tmp = n_doors + EVAC_Node_List(HPT%I_DOOR_NODES(i))%Node_Index 
             END IF
             IF (HPT%P_VENT_FFIELDS(i) < 1.0_EB) THEN
                CALL RANDOM_NUMBER(RN)
                IF ( RN < HPT%P_VENT_FFIELDS(i) ) THEN
                   Is_Known_Door(i_tmp) = .TRUE.
                ELSE
                   Is_Known_Door(i_tmp) = .FALSE.
                END IF
             ELSE
                Is_Known_Door(i_tmp) = .TRUE.
             END IF
          END DO

          ! Save the random known door information
          i_tmp = COUNT(Is_Known_Door)
          IF (j > 0) THEN    ! The agent is a member of a group
             Group_Known_Doors(j)%N_nodes = i_tmp
             IF (COUNT(Is_Known_Door) > 0) THEN
                ALLOCATE(Group_Known_Doors(j)%I_nodes(i_tmp), STAT=IZERO)
                CALL ChkMemErr('Change_Target_Door', 'Group_Known_Doors',IZERO) 
                i_tmp = 0
                DO i = 1, n_doors
                   IF (Is_Known_Door(i)) THEN
                      i_tmp = i_tmp + 1
                      Group_Known_Doors(j)%I_nodes(i_tmp) = EVAC_DOORS(i)%INODE  
                   END IF
                END DO
                DO i = 1, n_exits
                   IF (Is_Known_Door(i+n_doors)) THEN
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
                DO i = 1, n_doors
                   IF (Is_Known_Door(i)) THEN
                      i_tmp = i_tmp + 1
                      Human_Known_Doors(j1)%I_nodes(i_tmp) = EVAC_DOORS(i)%INODE  
                   END IF
                END DO
                DO i = 1, n_exits
                   IF (Is_Known_Door(i+n_doors)) THEN
                      i_tmp = i_tmp + 1
                      Human_Known_Doors(j1)%I_nodes(i_tmp) = EVAC_EXITS(i)%INODE  
                   END IF
                END DO
             END IF   ! Is there any known doors for this group?
          END IF    ! Is the agent a member of a group or not?
          
          ! Check that the door is in correct main evac mesh.
          DO i = 1, n_doors
             IF ( EVAC_DOORS(i)%IMESH /= HPT%imesh ) THEN
                Is_Known_Door(i) = .FALSE.
             END IF
          END DO
          DO i = 1, n_exits
             IF ( EVAC_EXITS(i)%IMESH /= HPT%imesh .OR. EVAC_EXITS(i)%COUNT_ONLY ) THEN
                Is_Known_Door(n_doors+i) = .FALSE.
             END IF
          END DO
       END IF  ! imode=0, Initialization call

    ELSE
       ! The agent is from an entry line. P_door= 0.0 or 1.0, known doors
       x1_old = xx
       y1_old = yy
       IEL    = ABS(HR%IEL)
       Speed  = HR%Speed
       DO i = 1, n_doors
          IF ( EVAC_DOORS(i)%IMESH == nm_tmp) THEN
             Is_Visible_Door(i) = .TRUE.
          END IF
       END DO
       DO i = 1, n_exits
          IF ( EVAC_EXITS(i)%IMESH == nm_tmp .AND. .NOT. EVAC_EXITS(i)%COUNT_ONLY ) THEN
             Is_Visible_Door(n_doors+i) = .TRUE.
          END IF
       END DO
       DO i = 1, PNX%N_VENT_FFIELDS 
          ! Check that the door/exit is in the correct main evac grid.
          i_tmp = 1
          IF (TRIM(EVAC_Node_List(PNX%I_DOOR_NODES(i))%Node_Type) == 'Door') THEN
             i_tmp = EVAC_Node_List(PNX%I_DOOR_NODES(i))%Node_Index 
          END IF
          IF (TRIM(EVAC_Node_List(PNX%I_DOOR_NODES(i))%Node_Type) == 'Exit' ) THEN
             i_tmp = n_doors + EVAC_Node_List(PNX%I_DOOR_NODES(i))%Node_Index 
          END IF
          IF ( Is_Visible_Door(i_tmp) ) THEN
             ! Door/exit is on this main evac mesh.
             IF (PNX%P_VENT_FFIELDS(i) < 0.5_EB) THEN
                Is_Known_Door(i_tmp) = .FALSE.
             ELSE
                Is_Known_Door(i_tmp) = .TRUE.
             END IF
          END IF
       END DO
    END IF   ! Is the agent from an entr or from an evac line

    ! Find the visible doors.
    DO i = 1, n_doors + n_exits
       IF ( Is_Visible_Door(i) ) THEN
          IF (EVAC_Node_List(n_egrids+n_entrys+i)%Node_Type == 'Door' ) THEN
             X11 = EVAC_DOORS(i)%X 
             Y11 = EVAC_DOORS(i)%Y
          ELSE
             X11 = EVAC_EXITS(i-n_doors)%X 
             Y11 = EVAC_EXITS(i-n_doors)%Y
          END IF
          ! Groups: the first member (x1_old,y1_old) of the group is used.
          PP_see_door = See_door(nm_tmp, i, HR%I_Target, x1_old, y1_old, x11, y11, ave_K, max_fed)
          FED_max_Door(i) = max_fed
          K_ave_Door(i) = ave_K 

          ! Note: a DOOR is not counted as visible door, if it does not have an
          ! EXIT_SIGN, unless it is already been a target door for this agent/group.
          IF (PP_see_door) THEN
             IF (EVAC_Node_List(n_egrids+n_entrys+i)%Node_Type == 'Door') THEN
                IF (.NOT. EVAC_DOORS(i)%EXIT_SIGN .AND. .NOT. HR%I_Target == i) THEN
                   Is_Visible_Door(i) = .FALSE.
                END IF
             END IF
          ELSE
             Is_Visible_Door(i) = .FALSE.
          END IF
       END IF ! correct main evac mesh
    END DO ! all doors and exits


    ! Note: I_Target < 0: not visible, >0: visible
    IF (ANY(Is_Visible_Door) .AND. imode == 1) THEN
       DO i = 1, n_doors + n_exits
          IF (ABS(HR%I_Target) == i .AND. Is_Visible_Door(i)) Is_Known_Door(i) = .TRUE.
       END DO
    END IF

    DO i = 1, n_doors
       ! If ( EVAC_DOORS(i)%TIME_OPEN > T .Or. EVAC_DOORS(i)%TIME_CLOSE < T) Then
       IF ( ABS(EVAC_DOORS(i)%IMODE)==2) THEN
          Is_Visible_Door(i) = .FALSE.
          Is_Known_Door(i) = .FALSE.
       END IF
    END DO
    DO i = 1, n_exits
       ! If ( (EVAC_EXITS(i)%TIME_OPEN > T .Or. EVAC_EXITS(i)%TIME_CLOSE < T) .And. .Not. EVAC_EXITS(i)%COUNT_ONLY ) Then
       IF ( (ABS(EVAC_EXITS(i)%IMODE)==2) .AND. .NOT. EVAC_EXITS(i)%COUNT_ONLY ) THEN
          Is_Visible_Door(n_doors+i) = .FALSE.
          Is_Known_Door(n_doors+i) = .FALSE.
       END IF
    END DO

    IF (ANY(Is_Known_Door) .OR. ANY(Is_Visible_Door)) THEN
       i_tmp   = 0
       L2_min = HUGE(L2_min)
       DO i = 1, n_doors + n_exits
          IF ( Is_Known_Door(i) .AND. Is_Visible_Door(i) ) THEN
             x_o = 0.0_EB
             y_o = 0.0_EB
             N_queue = 0
             IF (TRIM(EVAC_Node_List(n_egrids+n_entrys+i)%Node_Type) == 'Door' ) THEN
                x_o = EVAC_DOORS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%X
                y_o = EVAC_DOORS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%Y
                i_o = EVAC_DOORS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%I_VENT_FFIELD
                T_tmp1 = 50.0_EB*SQRT((x1_old-x_o)**2 + (y1_old-y_o)**2)/ &
                     EVAC_DOORS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%R_NTARGET + 1.0_EB
                ii = MIN(50,MAX(1,INT(T_tmp1)-1))
                Width = EVAC_DOORS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%Width
                N_queue = EVAC_DOORS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%NTARGET(ii)
             ELSE      ! 'Exit'
                x_o = EVAC_EXITS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%X
                y_o = EVAC_EXITS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%Y
                i_o = EVAC_EXITS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%I_VENT_FFIELD
                T_tmp1 = 50.0_EB*SQRT((x1_old-x_o)**2 + (y1_old-y_o)**2)/ &
                     EVAC_EXITS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%R_NTARGET + 1.0_EB
                ii = MIN(50,MAX(1,INT(T_tmp1)-1))
                Width = EVAC_EXITS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%Width
                N_queue = EVAC_EXITS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%NTARGET(ii)
             END IF
             IF (FED_DOOR_CRIT > 0.0_EB) THEN
                L2_tmp = FED_max_Door(i) * SQRT((x1_old-x_o)**2 + (y1_old-y_o)**2)/Speed
             ELSE
                L2_tmp = K_ave_Door(i)
             END IF
             IF (i_o == i_old_ffield) L2_tmp = 0.1_EB*L2_tmp
             IF (FAC_DOOR_QUEUE > 0.001_EB) THEN
                T_tmp  = SQRT((x_o-x1_old)**2 + (y_o-y1_old)**2)
                T_tmp1 = MIN(1.5_EB*Pi*T_tmp**2/(FAC_DOOR_QUEUE*Width), REAL(N_queue,EB)/(FAC_DOOR_QUEUE*Width))
                T_tmp = (T_tmp/Speed) +  T_tmp1
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
          END IF
       END DO
       IF (i_tmp > 0 ) THEN
          ! Known and visible door, no smoke
          IF (EVAC_Node_List(n_egrids+n_entrys+i_tmp)%Node_Type == 'Door' ) THEN
             name_new_ffield = TRIM(EVAC_DOORS(i_tmp)%VENT_FFIELD)
             i_new_ffield = EVAC_DOORS(i_tmp)%I_VENT_FFIELD
          ELSE        ! 'Exit'
             name_new_ffield = TRIM(EVAC_EXITS(i_tmp-n_doors)%VENT_FFIELD)
             i_new_ffield = EVAC_EXITS(i_tmp-n_doors)%I_VENT_FFIELD
          END IF
          color_index = 1
       ELSE
          ! No visible known doors available, try non-visible known doors
          i_tmp   = 0
          L2_min = HUGE(L2_min)
          DO i = 1, n_doors + n_exits
             IF ( Is_Known_Door(i) .AND. .NOT. Is_Visible_Door(i) ) THEN
                x_o = 0.0_EB
                y_o = 0.0_EB
                IF (EVAC_Node_List(n_egrids+n_entrys+i)%Node_Type == 'Door' ) THEN
                   x_o = EVAC_DOORS( EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%X
                   y_o = EVAC_DOORS( EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%Y
                   i_o = EVAC_DOORS( EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%I_VENT_FFIELD
                ELSE    ! 'Exit'
                   x_o = EVAC_EXITS( EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%X
                   y_o = EVAC_EXITS( EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%Y
                   i_o = EVAC_EXITS( EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%I_VENT_FFIELD
                END IF
                IF (FED_DOOR_CRIT > 0.0_EB) THEN
                   L2_tmp = FED_max_Door(i)*SQRT((x1_old-x_o)**2 + (y1_old-y_o)**2)/Speed
                ELSE
                   l2_tmp = K_ave_Door(i)
                END IF
                IF (i_o == i_old_ffield) L2_tmp = 0.1_EB*L2_tmp
                T_tmp  = SQRT((x_o-x1_old)**2 + (y_o-y1_old)**2)
                IF (i_o == i_old_ffield) T_tmp = T_tmp*FAC_DOOR_WAIT
                IF (T_tmp < L2_min .AND. L2_tmp < ABS(FED_DOOR_CRIT)) THEN
                   L2_min = MAX(0.0_EB,T_tmp)
                   i_tmp = i
                END IF
             END IF
          END DO
          IF (i_tmp > 0 ) THEN
             !  Non-visible known door, no smoke is found
             IF (EVAC_Node_List( n_egrids+n_entrys+i_tmp)%Node_Type == 'Door' ) THEN
                name_new_ffield = TRIM(EVAC_DOORS(i_tmp)%VENT_FFIELD)
                i_new_ffield = EVAC_DOORS(i_tmp)%I_VENT_FFIELD
             ELSE      ! 'Exit'
                name_new_ffield = TRIM(EVAC_EXITS(i_tmp-n_doors)%VENT_FFIELD)
                i_new_ffield = EVAC_EXITS(i_tmp-n_doors)%I_VENT_FFIELD
             END IF
             color_index = 2
          ELSE
             ! Known doors with no smoke have not been found.
             ! Try visible, not known door with no smoke.
             i_tmp   = 0
             L2_min = HUGE(L2_min)
             DO i = 1, n_doors + n_exits
                IF (Is_Visible_Door(i)) THEN
                   x_o = 0.0_EB
                   y_o = 0.0_EB
                   IF (EVAC_Node_List(n_egrids+n_entrys+i)%Node_Type == 'Door' ) THEN
                      x_o = EVAC_DOORS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%X
                      y_o = EVAC_DOORS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%Y
                      i_o = EVAC_DOORS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%I_VENT_FFIELD
                   ELSE  ! 'Exit'
                      x_o = EVAC_EXITS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%X
                      y_o = EVAC_EXITS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%Y
                      i_o = EVAC_EXITS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%I_VENT_FFIELD
                   END IF
                   IF (FED_DOOR_CRIT > 0.0_EB) THEN
                      L2_tmp = FED_max_Door(i)*SQRT((x1_old-x_o)**2 + (y1_old-y_o)**2)/Speed
                   ELSE
                      l2_tmp = K_ave_Door(i)
                   END IF
                   IF (i_o == i_old_ffield) L2_tmp = 0.1_EB*L2_tmp
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
                IF (EVAC_Node_List( n_egrids+n_entrys+i_tmp)%Node_Type == 'Door' ) THEN
                   name_new_ffield = TRIM(EVAC_DOORS(i_tmp)%VENT_FFIELD)
                   i_new_ffield = EVAC_DOORS(i_tmp)%I_VENT_FFIELD
                ELSE    ! 'Exit'
                   name_new_ffield = TRIM( EVAC_EXITS(i_tmp-n_doors)%VENT_FFIELD)
                   i_new_ffield = EVAC_EXITS(i_tmp-n_doors)%I_VENT_FFIELD
                END IF
                color_index = 3
             ELSE
                ! Now we have some smoke and some visible or known doors
                i_tmp   = 0
                L2_min = HUGE(L2_min)
                DO i = 1, n_doors + n_exits
                   IF ( Is_Known_Door(i) .OR. Is_Visible_Door(i) ) THEN
                      x_o = 0.0_EB
                      y_o = 0.0_EB
                      IF (EVAC_Node_List(n_egrids+n_entrys+i)%Node_Type == 'Door' ) THEN
                         x_o = EVAC_DOORS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%X
                         y_o = EVAC_DOORS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%Y
                         i_o = EVAC_DOORS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%I_VENT_FFIELD
                      ELSE ! 'Exit'
                         x_o = EVAC_EXITS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%X
                         y_o = EVAC_EXITS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%Y
                         i_o = EVAC_EXITS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%I_VENT_FFIELD
                      END IF
                      IF (FED_DOOR_CRIT > 0.0_EB) THEN
                         IF (j > 0) THEN
                            L2_tmp = Group_List(j)%IntDose + FED_max_Door(i) * SQRT((x1_old-x_o)**2+(y1_old-y_o)**2)/Speed
                         ELSE
                            L2_tmp = HR%IntDose + FED_max_Door(i) * SQRT((x1_old-x_o)**2+(y1_old-y_o)**2)/Speed
                         END IF
                      ELSE
                         l2_tmp = (SQRT((x1_old-x_o)**2+(y1_old-y_o)**2)*0.5_EB)/(3.0_EB/K_ave_Door(i))
                      END IF
                      IF (i_o == i_old_ffield) L2_tmp = 0.9_EB*L2_tmp
                      IF (L2_tmp < L2_min) THEN
                         L2_min = L2_tmp
                         i_tmp = i
                      END IF
                   END IF
                END DO

                IF (i_tmp > 0 .AND. L2_min < 1.0_EB) THEN
                   ! Not too much smoke, (known and/or visible doors)
                   IF (EVAC_Node_List(n_egrids+n_entrys+i_tmp)%Node_Type == 'Door' ) THEN
                      name_new_ffield = TRIM(EVAC_DOORS(i_tmp)%VENT_FFIELD)
                      i_new_ffield = EVAC_DOORS(i_tmp)%I_VENT_FFIELD
                   ELSE  ! 'Exit'
                      name_new_ffield = TRIM(EVAC_EXITS(i_tmp-n_doors)%VENT_FFIELD)
                      i_new_ffield = EVAC_EXITS(i_tmp-n_doors)%I_VENT_FFIELD
                   END IF
                   IF (Is_Known_Door(i_tmp) .AND. Is_Visible_Door(i_tmp)) color_index = 4
                   IF (Is_Known_Door(i_tmp) .AND. .NOT. Is_Visible_Door(i_tmp)) color_index = 5
                   IF (.NOT. Is_Known_Door(i_tmp) .AND. Is_Visible_Door(i_tmp)) color_index = 6
                ELSE    ! no match 
                   ! No door found (or too much smoke), use the main evac grid ffield
                   ! Note: This should be changed in later versions of the program.
                   i_tmp = 0
                   name_new_ffield = TRIM(MESH_NAME(nm_tmp))
                   i_new_ffield    = nm_tmp
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
          IF (i_tmp > 0 .AND. i_tmp <= n_doors ) color_index = EVAC_DOORS(i_tmp)%Avatar_Color_Index
          IF (i_tmp > n_doors .AND. i_tmp <= n_doors + n_exits) &
               color_index = EVAC_EXITS(i_tmp-n_doors)%Avatar_Color_Index
       END IF
       IF (imode == 2) THEN   ! check_target_node calls
          I_Target = i_tmp
          IF (i_tmp > 0 .AND. .NOT. Is_Visible_Door(MAX(1,i_tmp)) ) THEN
             ! I_Target >0: visible, <0: not visible
             I_Target = -i_tmp
          END IF
          I_Color  = color_index
          I_Field  = i_new_ffield
          RETURN
       END IF

    ELSE ! No known/visible door
       i_tmp = 0 ! no door found
       color_index = EVAC_AVATAR_NCOLOR ! default, cyan
       IF (imode == 2) THEN   ! check_target_node calls
          I_Target = 0
          I_Color  = color_index
          I_Field  = nm_tmp
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
    END IF ! Any known or visible door
    HR%I_Target = i_tmp
    IF (imode < 2) I_Target = i_tmp
    IF (i_tmp > 0 .AND. .NOT. Is_Visible_Door(MAX(1,i_tmp))) THEN
       ! I_Target > 0: visible, < 0: not visible
       HR%I_Target = -i_tmp
       IF (imode < 2) I_Target = -i_tmp
    END IF
    IF (j > 0) Group_Known_Doors(j)%I_Target = HR%I_Target
    IF ( (i_new_ffield /= i_old_ffield) .OR. (imode == 0) ) THEN
       ! Change the field of this group/agent.
       IF ( j == 0 ) THEN
          ! Group index=0: 'lost souls' or lonely agents
          HR%I_FFIELD    = i_new_ffield
          HR%FFIELD_NAME = TRIM(name_new_ffield)
          IF (COLOR_METHOD == 5) HR%COLOR_INDEX = color_index
          IF (COLOR_METHOD == 4) HR%COLOR_INDEX = color_index
          IF (FAC_DOOR_QUEUE > 0.001_EB) RETURN
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
          IF (FAC_DOOR_QUEUE > 0.001_EB) RETURN
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
  !
  SUBROUTINE GET_REV_evac(MODULE_REV,MODULE_DATE)
    !
    ! Passed variables
    INTEGER,INTENT(INOUT) :: MODULE_REV
    CHARACTER(255),INTENT(INOUT) :: MODULE_DATE
    !
    ! Local variables
    !

    WRITE(MODULE_DATE,'(A)') evacrev(INDEX(evacrev,':')+1:LEN_TRIM(evacrev)-2)
    READ (MODULE_DATE,'(I5)') MODULE_REV
    WRITE(MODULE_DATE,'(A)') evacdate

  END SUBROUTINE GET_REV_evac
  
END MODULE EVAC
