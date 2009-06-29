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
Module EVAC
  !
  Use PRECISION_PARAMETERS
  Use MESH_VARIABLES
  Use GLOBAL_CONSTANTS
  Use TRAN
  Use TYPES
  Use STAT
  Use COMP_FUNCTIONS
  Use MATH_FUNCTIONS
  Use MEMORY_FUNCTIONS
  Use MESH_POINTERS
!  Use MESH_POINTERS, ONLY: DT,IJKW,BOUNDARY_TYPE,XW,YW,WALL_INDEX,POINT_TO_MESH
!  Use EVAC_MESH_POINTERS
  Use PHYSICAL_FUNCTIONS, Only : GET_MASS_FRACTION, FED
  Use DCDFLIB, Only :  DCDFLIB_Gamma => Gamma
  !
  Implicit None
  !
  Character(255), Parameter :: evacid='$Id$'
  Character(255), Parameter :: evacrev='$Revision$'
  Character(255), Parameter :: evacdate='$Date$'
  !
  Private
  ! Public subprograms (called from the main program or read or dump)
  Public EVACUATE_HUMANS, INITIALIZE_EVACUATION, INIT_EVAC_GROUPS
  Public READ_EVAC, DUMP_EVAC, DUMP_EVAC_CSV, PREPARE_TO_EVACUATE
  Public EVAC_MESH_EXCHANGE, INITIALIZE_EVAC_DUMPS, GET_REV_EVAC
  ! Public variables (needed in the main program):
  ! Public variables (needed in the dump routine):
  Public N_DOORS, N_EXITS, N_ENTRYS, N_SSTANDS, EVAC_DOORS, EVAC_EXITS, EVAC_ENTRYS, EVAC_SSTANDS, & 
       EVAC_EXIT_TYPE, EVAC_DOOR_TYPE, EVAC_ENTR_TYPE, EVAC_SSTAND_TYPE
  !
  Character(255):: EVAC_VERSION = '2.1.2'
  Character(255) :: EVAC_COMPILE_DATE
  Integer :: EVAC_MODULE_REV
  !
  ! This is a group of persons, who are initialized together,
  ! i.e., they have same mass, speed, etc distributions and
  ! they are all put in the given rectangle.
  ! (&EVAC lines)
  Type EVACUATION_TYPE
     Real(EB) :: X1=0._EB,X2=0._EB,Y1=0._EB,Y2=0._EB,Z1=0._EB,Z2=0._EB,T_START=0._EB, Angle=0._EB
     Character(60) :: CLASS_NAME='null', ID='null'
     Character(30) :: GRID_NAME='null'
     Logical :: EVACFILE=.False., After_Tpre=.False., No_Persons=.False., SHOW=.True.
     Integer :: N_INITIAL=0,SAMPLING=0, IPC=0, IMESH=0
     Integer :: GN_MIN=0, GN_MAX=0
     Integer :: N_VENT_FFIELDS=0, Avatar_Color_Index=0
     Integer, Dimension(3) :: RGB=-1, AVATAR_RGB=-1
     Integer, Pointer, Dimension(:) :: I_DOOR_NODES =>NULL()
     Integer, Pointer, Dimension(:) :: I_VENT_FFIELDS =>NULL()
     Real(EB), Pointer, Dimension(:) :: P_VENT_FFIELDS =>NULL()
  End Type EVACUATION_TYPE
  !
  ! An evacuatio hole, i.e., a rectangle where humans should
  ! not be put.  This makes the &EVAC lines easier to define.
  ! (&EVHO lines)
  Type EVAC_HOLE_TYPE
     Real(EB) :: X1=0._EB,X2=0._EB,Y1=0._EB,Y2=0._EB,Z1=0._EB,Z2=0._EB
     Character(60) :: ID='null', PERS_ID='null', EVAC_ID='null'
     Character(30) :: GRID_NAME='null'
     Integer, Dimension(3) :: RGB=-1
     Logical :: SHOW=.True.
     Integer :: IMESH=0
  End Type EVAC_HOLE_TYPE
  !
  ! A spectator stand. IOR: which x,y line is the bottom line of the stand.
  ! ior=+1 x=x2, ior=-1 x=x1, ior=+2 y=y2, ior=-2 y=y1
  ! H is the height of the stand, S is the length along the incline.
  ! (&EVSS lines)
  Type EVAC_SSTAND_TYPE
     Real(EB) :: X1=0._EB,X2=0._EB,Y1=0._EB,Y2=0._EB,Z1=0._EB,Z2=0._EB, H=0._EB, H0=0._EB, S=0._EB
     Real(EB) :: Esc_SpeedUp=0._EB, Esc_SpeedDn=0._EB
     Real(EB) :: FAC_V0_UP=1._EB, FAC_V0_DOWN=1._EB, FAC_V0_HORI=1._EB
     Real(EB) :: cos_x=1._EB, cos_y=1._EB, sin_x=0._EB, sin_y=0._EB
     Character(60) :: ID='null'
     Character(26) :: GRID_NAME='null'
     Integer, Dimension(3) :: RGB=-1
     Integer :: IMESH=0, IOR=0
     Real(EB) :: UBAR0=0._EB, VBAR0=0._EB
     Logical :: Use_v0=.False., SHOW=.True.
     Real(EB), Dimension(3) :: ORIENTATION=0.0_EB
  End Type EVAC_SSTAND_TYPE
  !
  ! Humans belong to some small group (1 to about 5 persons).  This type
  ! collects the common properties of the group.
  Type GROUP_TYPE
     Real(EB) :: GROUP_X=0._EB, GROUP_Y=0._EB, MAX_DIST_CENTER=0._EB, LIMIT_COMP=0._EB
     Real(EB) :: GROUP_EFF=0._EB, RADIUS_COMPLETE_0=0._EB, RADIUS_COMPLETE_1=0._EB
     Real(EB) :: Speed=0._EB, IntDose=0._EB, Tpre=0._EB, Tdoor=0._EB, Tdet=0._EB
     Integer :: GROUP_SIZE=0, GROUP_ID=0, COMPLETE=0, IEL=0, Avatar_Color_Index=0
     Integer, Dimension(3) :: AVATAR_RGB=-1
     Integer, Pointer, Dimension(:) :: GROUP_I_FFIELDS =>NULL()
  End Type GROUP_TYPE
  
  Type KNOWN_DOOR_TYPE
     Integer :: N_nodes=0, I_Target=0
     Integer, Pointer, Dimension(:) :: I_nodes =>NULL()
  End Type KNOWN_DOOR_TYPE
  !
  ! This defines a class of persons, e.g. soccer fan.
  ! (&PERS lines)
  Type EVAC_PERS_TYPE
     Real(EB) :: D_mean=0._EB, D_para=0._EB, D_para2=0._EB, D_low=0._EB, D_high=0._EB
     Real(EB) :: V_mean=0._EB, V_para=0._EB, V_para2=0._EB, V_low=0._EB, V_high=0._EB
     Real(EB) :: Tau_mean=0._EB, Tau_para=0._EB, Tau_para2=0._EB, Tau_low=0._EB, Tau_high=0._EB
     Real(EB) :: Tpre_mean=0._EB, Tpre_para=0._EB, Tpre_para2=0._EB, Tpre_low=0._EB, Tpre_high=0._EB
     Real(EB) :: Tdet_mean=0._EB, Tdet_para=0._EB, Tdet_para2=0._EB, Tdet_low=0._EB, Tdet_high=0._EB
     Real(EB) :: A=0._EB,B=0._EB,Lambda=0._EB,C_Young=0._EB,Gamma=0._EB,Kappa=0._EB
     Real(EB) :: r_torso=0._EB,r_shoulder=0._EB,d_shoulder=0._EB,m_iner=0._EB, Tau_iner=0._EB
     Real(EB) :: FAC_V0_UP=-1._EB, FAC_V0_DOWN=-1._EB
     Character(60) :: ID='null'
     Integer :: I_DIA_DIST=0, I_VEL_DIST=0, I_PRE_DIST=0, I_DET_DIST=0, I_TAU_DIST=0
     Integer :: Avatar_Color_Index=0
     Integer, Dimension(3) :: RGB=-1, AVATAR_RGB=-1
  End Type EVAC_PERS_TYPE
  !
  ! Exit door type: this just count the number of persons
  ! T_first: first person's exit time (saved for output)
  ! CHECK_FLOW: If true then the flow can not exceed Flow_max
  ! (&EXIT lines)
  Type EVAC_EXIT_TYPE
     Real(EB) :: T_first=0._EB, T_last=0._EB, Flow_max=0._EB, Width=0._EB, Height=2.0_EB
     Real(EB) :: X1=0._EB, X2=0._EB, Y1=0._EB, Y2=0._EB, Z1=0._EB, Z2=0._EB, &
          X=0._EB, Y=0._EB, Z=0._EB, Xsmoke=0._EB, Ysmoke=0._EB, Zsmoke=0._EB, &
          TIME_OPEN=0._EB, TIME_CLOSE=0._EB, R_NTARGET=0._EB
     Integer :: IOR=0, ICOUNT=0, IMESH=0, INODE=0, IMODE=1
     Integer, Dimension(50) :: NTARGET=0
     Real(EB) :: FED_CO_CO2_O2=0._EB, SOOT_DENS=0._EB, TMP_G=0._EB, RADFLUX=0._EB
     Integer :: II=0, JJ=0, KK=0, FED_MESH=0
     Logical :: CHECK_FLOW=.False., COUNT_ONLY=.False., SHOW=.True.
     Integer :: STR_INDX=0, STR_SUB_INDX=0
     Character(60) :: ID='null'
     Character(60) :: TO_NODE='null'
     Character(30) :: GRID_NAME='null'
     Character(26) :: VENT_FFIELD='null'
     Integer :: I_VENT_FFIELD=0, Avatar_Color_Index=0
     Integer, Dimension(3) :: RGB=-1
     Real(EB), Dimension(3) :: ORIENTATION=0.0_EB
  End Type EVAC_EXIT_TYPE
  !
  ! Like exit, but door will always put the persons to some
  ! other node. (Thus no count_only option.)
  ! (&DOOR lines)
  Type EVAC_DOOR_TYPE
     Real(EB) :: T_first=0._EB, T_last=0._EB, Flow_max=0._EB, Width=0._EB, Height=2.0_EB
     Real(EB) :: X1=0._EB, X2=0._EB, Y1=0._EB, Y2=0._EB, Z1=0._EB, Z2=0._EB, &
          X=0._EB, Y=0._EB, Z=0._EB, Xsmoke=0._EB, Ysmoke=0._EB, Zsmoke=0._EB, &
          TIME_OPEN=0._EB, TIME_CLOSE=0._EB, R_NTARGET=0._EB
     Integer :: IOR=0, ICOUNT=0, INODE=0, INODE2=0, IMESH=0, IMESH2=0, IMODE=1
     Integer, Dimension(50) :: NTARGET=0
     Integer :: STR_INDX=0, STR_SUB_INDX=0
     Real(EB) :: FED_CO_CO2_O2=0._EB, SOOT_DENS=0._EB, TMP_G=0._EB, RADFLUX=0._EB
     Integer :: II=0, JJ=0, KK=0, FED_MESH=0
     Logical :: CHECK_FLOW=.False., EXIT_SIGN=.False., KEEP_XY=.False., SHOW=.True.
     Character(60) :: ID='null'
     Character(60) :: TO_NODE='null'
     Character(30) :: GRID_NAME='null'
     Character(26) :: VENT_FFIELD='null'
     Integer :: I_VENT_FFIELD=0, Avatar_Color_Index=0
     Integer, Dimension(3) :: RGB=-1
     Real(EB), Dimension(3) :: ORIENTATION=0.0_EB
  End Type EVAC_DOOR_TYPE
  !
  ! Like door, but corr will model stairs (or corridors). 
  ! The parameters, like velocity as function of density etc.
  ! define if it is corridor or stairway
  ! (&CORR lines)
  Type EVAC_CORR_TYPE
     Real(EB) :: T_first=0._EB, T_last=0._EB, Flow_max=0._EB, Width1=0._EB, Width2=0._EB
     Real(EB) :: X1=0._EB,X2=0._EB,Y1=0._EB,Y2=0._EB,Z1=0._EB,Z2=0._EB, Width=0._EB
     Real(EB) :: Eff_Width=0._EB, Eff_Length=0._EB, Eff_Area=0._EB, Fac_Speed=0._EB
     ! Note: Corridor may have 2 different points, where smoke etc. data
     ! is saved.
     Real(EB), Dimension(2) :: FED_CO_CO2_O2=0._EB, SOOT_DENS=0._EB, TMP_G=0._EB, RADFLUX=0._EB
     Integer :: FED_MESH=0, FED_MESH2=0
     Integer, Dimension(2) :: II=0, JJ=0, KK=0
     Integer :: IOR=0, ICOUNT=0, INODE=0, INODE2=0, IMESH=0, IMESH2=0
     Integer :: MAX_HUMANS_INSIDE=0, n_inside=0
     Logical :: CHECK_FLOW=.False.
     Integer, Dimension(3) :: RGB=-1
     Character(60) :: ID='null'
     Character(60) :: TO_NODE='null'
     Character(30) :: GRID_NAME='null'
     Type (CORR_LL_Type), Pointer :: First =>NULL()
  End Type EVAC_CORR_TYPE
  ! 
  ! STRS is a construct to build a staircase. STRS consists of stairs and
  ! landings. 
  Type EVAC_STRS_TYPE
     Real(EB) :: XB(6)
     Real(EB), Pointer, Dimension(:,:)   :: XB_NODE =>NULL(), XB_CORE =>NULL()
     Real(EB) :: FAC_V0_HORI=1._EB, FAC_V0_DOWN=1._EB, FAC_V0_UP=1._EB
     Integer :: ICOUNT=0, INODE=0, INODE2=0, IMESH=0, IMESH2=0, N_CORES = 0
     Integer :: N_LANDINGS, N_NODES, N_NODES_OUT, N_NODES_IN
     Integer, Pointer, Dimension(:) :: NODE_IOR =>NULL(), NODE_TYPE =>NULL(), NODES_IN =>NULL()
     Integer, Pointer, Dimension(:) :: NODES_OUT =>NULL(), I_CORE =>NULL()
     Character(60) :: ID
     Character(24) :: MESH_ID
     Logical RIGHT_HANDED
  End Type EVAC_STRS_TYPE
  !
  ! This produces more humans on the floor specified by the
  ! coordinates. the person type ('soccer_fan' etc) are also
  ! defined here for these persons.
  ! (&ENTR lines)
  Type EVAC_ENTR_TYPE
     Real(EB) :: T_first=0._EB, T_last=0._EB, Flow=0._EB, Width=0._EB, T_Start=0._EB, T_Stop=0._EB
     Real(EB) :: X1=0._EB,X2=0._EB,Y1=0._EB,Y2=0._EB,Z1=0._EB,Z2=0._EB, Height=2.0_EB, Z=0.0_EB
     Integer :: IOR=0, ICOUNT=0, IPC=0, IMESH=0, INODE=0, IMODE=-1, &
          TO_INODE=0, N_Initial=0, Max_Humans=-1, &
          STR_INDX=0, STR_SUB_INDX=0
     Character(60) :: CLASS_NAME='null', ID='null'
     Character(60) :: TO_NODE='null'
     Character(30) :: GRID_NAME='null', Max_Humans_Ramp
     Logical :: After_Tpre=.False., No_Persons=.False., SHOW=.True.
     Integer :: N_VENT_FFIELDS=0, Avatar_Color_Index=0
     Integer, Pointer, Dimension(:) :: I_DOOR_NODES =>NULL()
     Integer, Pointer, Dimension(:) :: I_VENT_FFIELDS =>NULL()
     Real(EB), Pointer, Dimension(:) :: P_VENT_FFIELDS =>NULL()
     Integer, Dimension(3) :: RGB=-1, AVATAR_RGB=-1
     Real(EB), Dimension(3) :: ORIENTATION=0.0_EB
  End Type EVAC_ENTR_TYPE
  !
  ! coordinates. the person type ('soccer_fan' etc) are also
  ! defined here for these persons.
  Type EVAC_NODE_TYPE
     Integer :: Node_Index=0, IMESH=0
     Character(60) :: ID='null', Node_Type='null'
     Character(30) :: GRID_NAME='null'
  End Type EVAC_NODE_TYPE
  !
  ! Linked list, needed for the corridors
  Type CORR_LL_TYPE
     Type (HUMAN_TYPE) :: HUMAN
     Real(EB) :: T_in=0._EB, T_out=0._EB
     Logical :: From1_To2=.False.
     Integer :: Index=0
     Type (CORR_LL_Type), Pointer :: Next =>NULL()
  End Type CORR_LL_TYPE
  !
  ! Pointers to the allocatable arrays so one can use these as
  ! shorthands to the array elements.

  !
  ! Next holds door information for groups
  Type (KNOWN_DOOR_TYPE), Dimension(:), Allocatable, Target :: Group_Known_Doors
  ! Next holds door information for lonely humans (group_id=0)
  Type (KNOWN_DOOR_TYPE), Dimension(:), Allocatable, Target :: Human_Known_Doors
  Integer :: ilh, ilh_dim
  
  ! Holds the list of the different human groups, i33 is a running index 
  ! for the groups, i33_dim is last index, i.e., the dimension of the array.
  Type (GROUP_TYPE), Dimension(:), Allocatable, Target :: Group_List
  Integer :: i33, i33_dim

  ! Holds the information of the nodes
  Type (EVAC_NODE_Type), Dimension(:), Allocatable, Target :: Evac_Node_List

  ! Holds the information of the EVAC-lines.
  Type (EVACUATION_Type), Dimension(:), Allocatable, Target :: EVACUATION

  ! Holds the information of the EVHO-lines.
  Type (EVAC_HOLE_Type), Dimension(:), Allocatable, Target :: EVAC_HOLES

  ! Holds the information of the EVSS-lines.
  Type (EVAC_SSTAND_TYPE), Dimension(:), Allocatable, Target :: EVAC_SSTANDS

  ! Holds the information of the EXIT-lines.
  Type (EVAC_EXIT_Type), Dimension(:), Allocatable, Target :: EVAC_EXITS

  ! Holds the information of the DOOR-lines.
  Type (EVAC_DOOR_Type), Dimension(:), Allocatable, Target :: EVAC_DOORS

  ! Holds the information of the ENTR-lines.
  Type (EVAC_ENTR_Type), Dimension(:), Allocatable, Target :: EVAC_ENTRYS

  ! Holds the information of the CORR-lines.
  Type (EVAC_CORR_Type), Dimension(:), Allocatable, Target :: EVAC_CORRS

  ! Holds the information on the STRS-lines.
  Type (EVAC_STRS_Type), Dimension(:), Allocatable, Target :: EVAC_STRS

  ! Holds the information of the PERS-lines.
  Type (EVAC_PERS_Type), Dimension(:), Allocatable, Target :: EVAC_PERSON_CLASSES
  !
  ! Next are needed for the Gaussian random numbers
  Integer GaussFlag
  Real(EB) GaussSet1, GaussSet2
  Integer GTrunFlag
  Real(EB) GTrunSet1, GTrunSet2
  !
  Integer :: NPC_EVAC, NPC_PERS, N_EXITS, N_DOORS, N_ENTRYS, &
       N_CORRS, N_EGRIDS, N_NODES, N_HOLES, N_SSTANDS, N_STRS, N_CO_EXITS
  Integer :: NPPS
  Integer :: ILABEL_last
  Character(100) :: MESSAGE
  Real(FB) :: EVAC_Z_MIN, EVAC_Z_MAX
  !
  Real(EB), Dimension(:,:), Allocatable :: TT_Evac, FF_Evac
  Integer, Dimension(:), Allocatable :: NTT_Evac
  !
  !
  Logical :: NOT_RANDOM
  Integer :: I_FRIC_SW, COLOR_METHOD
  Real(EB) ::  FAC_A_WALL, FAC_B_WALL, LAMBDA_WALL, &
       NOISEME, NOISETH, NOISECM, RADIUS_COMPLETE_0, &
       RADIUS_COMPLETE_1, GROUP_EFF, FED_DOOR_CRIT, &
       TDET_SMOKE_DENS, DENS_INIT, EVAC_DT_MAX, GROUP_DENS, &
       FC_DAMPING, EVAC_DT_MIN, V_MAX, V_ANGULAR_MAX, V_ANGULAR, &
       SMOKE_MIN_SPEED, SMOKE_MIN_SPEED_VISIBILITY, TAU_CHANGE_DOOR, &
       HUMAN_SMOKE_HEIGHT, TAU_CHANGE_V0, THETA_SECTOR, CONST_DF, FAC_DF, &
       CONST_CF, FAC_CF, FAC_1_WALL, FAC_2_WALL, FAC_V0_DIR, FAC_V0_NOCF, FAC_NOCF, &
       CF_MIN_A, CF_FAC_A_WALL, CF_MIN_TAU, CF_MIN_TAU_INER, CF_FAC_TAUS, FAC_DOOR_QUEUE, &
       FAC_DOOR_WAIT, CF_MIN_B
  Integer, Dimension(3) :: DEAD_RGB
  !
  Real(EB), Dimension(:), Allocatable :: Tsteps
  !
  Real(EB), Dimension(:), Allocatable :: K_ave_Door
  Real(EB), Dimension(:), Allocatable :: FED_max_Door
  Logical, Dimension(:), Allocatable :: Is_Known_Door, Is_Visible_Door
  Integer, Dimension(:), Allocatable :: Color_Tmp
  !
  Integer :: n_dead=0, icyc_old=0, n_change_doors=0, n_change_trials=0
  Real(EB) :: fed_max_alive, fed_max
  !
  ! Stairs constants
  Integer :: STRS_LANDING_TYPE=1, STRS_STAIR_TYPE=2

 ! Human constants
  Integer :: HUMAN_SAME_MESH_TARGET=-2,HUMAN_IMPOSSIBLE_TARGET=-1, &
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
Contains
  !
  Subroutine READ_EVAC
    Implicit None
    !
    ! Local variables
    Integer :: NUMBER_INITIAL_PERSONS, COLOR_METHOD_TMP, &
         SAMPLING_FACTOR, IPC, n_tmp, GN_MIN, GN_MAX, N_RAMP_INI
    Real(EB) :: DTSAM
    Logical :: EVACFILE

    Real(EB) :: DUMMY
    Real(EB) :: XB(6), XB1(6), XB2(6)
    Real(EB), Dimension(3) :: XYZ, XYZ_SMOKE
    Integer :: IOS, IZERO, N, I, J, K, IOR
    Character(30) QUANTITY, MAX_HUMANS_RAMP
    Character(60) FYI,ID,PERS_ID,TO_NODE,EVAC_ID, DEFAULT_PROPERTIES
    Character(26) FLOW_FIELD_ID
    Integer :: DIAMETER_DIST,VELOCITY_DIST, &
         PRE_EVAC_DIST,DET_EVAC_DIST,TAU_EVAC_DIST
    Real(EB) :: VEL_MEAN,VEL_PARA,VEL_PARA2,VEL_LOW,VEL_HIGH, &
         DIA_MEAN,DIA_PARA,DIA_PARA2,DIA_LOW,DIA_HIGH, &
         PRE_MEAN,PRE_PARA,PRE_PARA2,PRE_LOW,PRE_HIGH, &
         DET_MEAN,DET_PARA,DET_PARA2,DET_LOW,DET_HIGH, &
         TAU_MEAN,TAU_PARA,TAU_PARA2,TAU_LOW,TAU_HIGH, &
         FCONST_A,FCONST_B,L_NON_SP, &
         C_YOUNG,GAMMA,KAPPA, ANGLE, &
         D_TORSO_MEAN,D_SHOULDER_MEAN, TAU_ROT, M_INERTIA
    Integer :: MAX_HUMANS_INSIDE, n_max_in_corrs, COLOR_INDEX, i_avatar_color, MAX_HUMANS
    Real(EB) :: MAX_FLOW, WIDTH, TIME_START, TIME_STOP, WIDTH1, &
         WIDTH2, EFF_WIDTH, EFF_LENGTH, FAC_SPEED, TIME_OPEN, TIME_CLOSE
    Real(EB) :: UBAR0, VBAR0
    Logical :: CHECK_FLOW, COUNT_ONLY, AFTER_REACTION_TIME, EXIT_SIGN, KEEP_XY, USE_V0, SHOW
    Logical :: OUTPUT_SPEED, OUTPUT_MOTIVE_FORCE, OUTPUT_FED, OUTPUT_OMEGA, &
         OUTPUT_ANGLE, OUTPUT_CONTACT_FORCE, OUTPUT_TOTAL_FORCE, OUTPUT_MOTIVE_ANGLE
    Integer, Dimension(3) :: RGB, AVATAR_RGB
    Character(26) :: VENT_FFIELD, MESH_ID, EVAC_MESH
    Real(EB) :: FAC_V0_UP, FAC_V0_DOWN, FAC_V0_HORI, HEIGHT, HEIGHT0, ESC_SPEED
    Character(25) :: COLOR, DEAD_COLOR, AVATAR_COLOR

    ! Stairs variables
    Real(EB) :: XB_CORE(4), XB_CORES(500,6), XB_LANDINGS(500,6), XB_STAIRS(500,8)
    Real(EB) VERTICAL_LANDING_SEPARATION, STR_Length, STR_Height
    Integer N_LANDINGS, NL, NODES_TMP(500)
    Logical :: RIGHT_HANDED, LEFT_HANDED

    Character(26), Dimension(51) :: KNOWN_DOOR_NAMES
    Real(EB), Dimension(51) :: KNOWN_DOOR_PROBS

    Integer :: ii,jj,kk

    Integer :: size_rnd
    Integer, Dimension(8) :: t_rnd
    Integer, Dimension(:), Allocatable :: seed_rnd

    Namelist /EXIT/ ID, XB, IOR, FLOW_FIELD_ID, CHECK_FLOW, &
         MAX_FLOW, FYI, COUNT_ONLY, WIDTH, XYZ, VENT_FFIELD, &
         MESH_ID, COLOR_INDEX, XYZ_SMOKE, &
         TIME_OPEN, TIME_CLOSE, EVAC_MESH, RGB, COLOR, SHOW, HEIGHT
    Namelist /DOOR/ ID, XB, IOR, FLOW_FIELD_ID, CHECK_FLOW, &
         MAX_FLOW, TO_NODE, FYI, WIDTH, XYZ, VENT_FFIELD, &
         EXIT_SIGN, MESH_ID, COLOR_INDEX, XYZ_SMOKE, KEEP_XY, &
         TIME_OPEN, TIME_CLOSE, EVAC_MESH, RGB, COLOR, SHOW, HEIGHT
    Namelist /ENTR/ ID, XB, IOR, FLOW_FIELD_ID, MAX_FLOW, &
         FYI, WIDTH, QUANTITY, PERS_ID, TIME_START, &
         TIME_STOP, AFTER_REACTION_TIME, &
         KNOWN_DOOR_NAMES, KNOWN_DOOR_PROBS, &
         MESH_ID, COLOR_INDEX, EVAC_MESH, RGB, COLOR, &
         AVATAR_COLOR, AVATAR_RGB, MAX_HUMANS, MAX_HUMANS_RAMP, SHOW, HEIGHT
    Namelist /CORR/ ID, XB, IOR, FLOW_FIELD_ID, CHECK_FLOW, &
         MAX_FLOW, TO_NODE, FYI, WIDTH, WIDTH1, WIDTH2, &
         EFF_WIDTH, EFF_LENGTH, MAX_HUMANS_INSIDE, FAC_SPEED, &
         XB1, XB2, RGB, COLOR
    Namelist /STRS/ ID, XB, XB_CORE, XB_CORES, TO_NODE, RIGHT_HANDED, LEFT_HANDED, MESH_ID, &
         N_LANDINGS, XB_LANDINGS, VERTICAL_LANDING_SEPARATION, &
         FAC_V0_UP, FAC_V0_DOWN, FAC_V0_HORI
    Namelist /EVAC/ NUMBER_INITIAL_PERSONS, QUANTITY, FYI, &
         ID, DTSAM, XB, FLOW_FIELD_ID, PERS_ID, &
         TIME_START, TIME_STOP, IOR, MAX_FLOW, WIDTH, ANGLE, &
         AFTER_REACTION_TIME, GN_MIN, GN_MAX, &
         KNOWN_DOOR_NAMES, KNOWN_DOOR_PROBS, MESH_ID, &
         COLOR_INDEX, EVAC_MESH, RGB, COLOR, &
         AVATAR_COLOR, AVATAR_RGB, SHOW
    Namelist /EVHO/ FYI, ID, XB, EVAC_ID, PERS_ID, MESH_ID, EVAC_MESH, RGB, COLOR, SHOW

    Namelist /EVSS/ FYI, ID, XB, MESH_ID, HEIGHT, HEIGHT0, IOR, &
         FAC_V0_UP, FAC_V0_DOWN, FAC_V0_HORI, ESC_SPEED, EVAC_MESH, RGB, COLOR, &
         UBAR0, VBAR0, USE_V0, SHOW

    Namelist /PERS/ FYI,ID,DIAMETER_DIST,VELOCITY_DIST, &
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
    If (.Not. ANY(EVACUATION_GRID)) Then
       N_EVAC = 0
       Return
    End If

    NPPS = 30000 ! Number Persons Per Set (dump to a file)
    !
    EVAC_DT = EVAC_DT_FLOWFIELD     ! Initialize the clock
    EVAC_CLOCK = T_BEGIN    ! clock for the CHID_evac.csv file
    EVAC_N_QUANTITIES = 0

    N_RAMP_INI = N_RAMP
    i33 = 0
    ilh = 0

    If (MYID==Max(0,EVAC_PROCESS)) Then
       Allocate(Tsteps(NMESHES),STAT=IZERO)
       Call ChkMemErr('READ','Tsteps',IZERO) 
       Tsteps(:) = EVAC_DT_FLOWFIELD
       If (Abs(TIME_SHRINK_FACTOR-1.0_EB) > 0.000000000001_EB ) Call SHUTDOWN('ERROR: Evac is not ready for TIME_SHRINK_FACTOR')
    End If
    !
    ! I_EVAC: 'binary' index:
    ! These are just initialization. Later it is checked if files
    ! exists. If EFF/FED file does not exists, then it is calculated 
    ! (and saved).
    I_EVAC = 16*0 + 8*0 + 4*0 + 2*0 + 1*0 ! do not save soot,fed files
    If (.Not. All(EVACUATION_ONLY) ) Then
       ! There are fire grids ==> save fed and evac flow fields
       I_EVAC = 16*1 + 8*0 + 4*0 + 2*1 + 1*1
       If (N_REACTIONS == 0) I_EVAC = 16*1 + 8*0 + 4*0 + 2*0 + 1*1
    Else
       If (EVACUATION_MC_MODE) Then
          ! There are no fire grids ==> try to read fed and evac flow fields if MC mode
          I_EVAC = 16*0 + 8*1 + 4*1 + 2*0 + 1*0
       Else
          ! There are no fire grids ==> try to read fed and recalculate evac flow fields if not MC mode
          I_EVAC = 16*0 + 8*1 + 4*0 + 2*0 + 1*1
       End If
    End If
    !
    ! Every human has an identification number, ILABEL_last is
    ! the last used number, so next human will have an
    ! identification number, which is ILABEL_last + 1
    ILABEL_last = 0


    Call COUNT_EVAC_NODES
    ! Write (lu_err,*) 'Evac: Counted evacuation nodes'
    Call READ_PERS
    ! Write (lu_err,*) 'Evac: Read person classes'
    Call READ_STRS
    ! Write (lu_err,*) 'Evac: Read stairs'
    Call READ_EXIT
    ! Write (lu_err,*) 'Evac: Read exits'
    Call READ_DOOR    
    ! Write (lu_err,*) 'Evac: Read doors'
    Call READ_CORR
    ! Write (lu_err,*) 'Evac: Read corridors'
    Call READ_ENTRIES
    ! Write (lu_err,*) 'Evac: Read entries'
    Call COLLECT_NODE_INFO
    ! Write (lu_err,*) 'Evac: Collected node info '
    Call READ_EVAC_LINES
    ! Write (lu_err,*) 'Evac: Read evac namelists'
    Call READ_EVHO
    ! Write (lu_err,*) 'Evac: Read evho namelists'
    Call READ_EVSS
    ! Write (lu_err,*) 'Evac: Read inclines'

    If (MYID /= Max(0,EVAC_PROCESS)) Return

    Call CHECK_EVAC_NODES

  Contains

    Subroutine COUNT_EVAC_NODES
      Implicit None
      !
      ! Determine total number of PERS lines in the input file
      !
      EVAC_AVATAR_NCOLOR = 0  ! Dimension of Avatar color table
      i_avatar_color     = 0  ! Counter for avatar colors
      COLOR_METHOD       = -1 ! Default is standard human colors in Smokeview
      DEAD_RGB           = (/  0,255,255/) ! cyan
      NPC_PERS = 0
      COUNT_PERS_LOOP: Do
         Call CHECKREAD('PERS',LU_INPUT,IOS) 
         If (IOS == 1) Then
            Exit COUNT_PERS_LOOP
         End If
         Read(LU_INPUT,NML=PERS,End=221,ERR=222,IOSTAT=IOS)
         NPC_PERS = NPC_PERS + 1
222      If (IOS > 0) Call SHUTDOWN('ERROR: Problem with PERS line')
      End Do COUNT_PERS_LOOP
221   Rewind(LU_INPUT)
      If (COLOR_METHOD == 3) EVAC_AVATAR_NCOLOR = NPC_PERS + 1
      COLOR_METHOD_TMP = COLOR_METHOD
      !
      ! Determine total number of EVAC lines in the input file
      !
      NPC_EVAC = 0
      COUNT_EVAC_LOOP: Do
         NUMBER_INITIAL_PERSONS = 0
         Call CHECKREAD('EVAC',LU_INPUT,IOS) 
         If (IOS == 1) Then
            Exit COUNT_EVAC_LOOP
         End If
         Read(LU_INPUT,NML=EVAC,End=219,ERR=220,IOSTAT=IOS)
         NPC_EVAC = NPC_EVAC + 1
         If (COLOR_METHOD == 0 .And. NUMBER_INITIAL_PERSONS > 0) Then
            EVAC_AVATAR_NCOLOR = EVAC_AVATAR_NCOLOR + 1
         End If
         !
220      If (IOS > 0) Call SHUTDOWN('ERROR: Problem with EVAC line')
      End Do COUNT_EVAC_LOOP
219   Rewind(LU_INPUT)
      !
      ! Determine total number of EXIT lines in the input file
      !
      N_EXITS = 0
      N_CO_EXITS = 0
      COUNT_EXITS_LOOP: Do
         COUNT_ONLY = .False.
         Call CHECKREAD('EXIT',LU_INPUT,IOS) 
         If (IOS == 1) Then
            Exit COUNT_EXITS_LOOP
         End If
         Read(LU_INPUT,NML=Exit,End=223,ERR=224,IOSTAT=IOS)
         N_EXITS = N_EXITS + 1
         If (COLOR_METHOD == 4 .And. .Not.COUNT_ONLY) Then
            EVAC_AVATAR_NCOLOR = EVAC_AVATAR_NCOLOR + 1
         End If
         If (COUNT_ONLY) N_CO_EXITS = N_CO_EXITS + 1
224      If (IOS > 0) Call SHUTDOWN('ERROR: Problem with EXIT line')
      End Do COUNT_EXITS_LOOP
223   Rewind(LU_INPUT)
      !
      ! Determine total number of DOOR lines in the input file
      !
      N_DOORS = 0
      COUNT_DOORS_LOOP: Do
         Call CHECKREAD('DOOR',LU_INPUT,IOS) 
         If (IOS == 1) Then
            Exit COUNT_DOORS_LOOP
         End If
         Read(LU_INPUT,NML=DOOR,End=225,ERR=226,IOSTAT=IOS)
         N_DOORS = N_DOORS + 1
         If (COLOR_METHOD == 4) Then
            EVAC_AVATAR_NCOLOR = EVAC_AVATAR_NCOLOR + 1
         End If
226      If (IOS > 0) Call SHUTDOWN('ERROR: Problem with DOOR line')
      End Do COUNT_DOORS_LOOP
225   Rewind(LU_INPUT)
      !
      ! Determine total number of ENTR lines in the input file
      !
      N_ENTRYS = 0
      COUNT_ENTRYS_LOOP: Do
         MAX_FLOW     = 0.0_EB
         Call CHECKREAD('ENTR',LU_INPUT,IOS) 
         If (IOS == 1) Then
            Exit COUNT_ENTRYS_LOOP
         End If
         Read(LU_INPUT,NML=ENTR,End=227,ERR=228,IOSTAT=IOS)
         N_ENTRYS = N_ENTRYS + 1
         If (COLOR_METHOD == 0 .And. MAX_FLOW > 0) Then
            EVAC_AVATAR_NCOLOR = EVAC_AVATAR_NCOLOR + 1
         End If
228      If (IOS > 0) Call SHUTDOWN('ERROR: Problem with ENTR line')
      End Do COUNT_ENTRYS_LOOP
227   Rewind(LU_INPUT)
      !
      ! Determine total number of CORR lines in the input file
      !
      N_CORRS = 0
      COUNT_CORRS_LOOP: Do
         Call CHECKREAD('CORR',LU_INPUT,IOS) 
         If (IOS == 1) Then
            Exit COUNT_CORRS_LOOP
         End If
         Read(LU_INPUT,NML=CORR,End=229,ERR=230,IOSTAT=IOS)
         N_CORRS = N_CORRS + 1
230      If (IOS > 0) Call SHUTDOWN('ERROR: Problem with CORR line')
      End Do COUNT_CORRS_LOOP
229   Rewind(LU_INPUT)
      !
      ! Determine total number of EVHO lines in the input file
      !
      N_HOLES = 0
      COUNT_EVHO_LOOP: Do
         Call CHECKREAD('EVHO',LU_INPUT,IOS) 
         If (IOS == 1) Then
            Exit COUNT_EVHO_LOOP
         End If
         Read(LU_INPUT,NML=EVHO,End=231,ERR=232,IOSTAT=IOS)
         N_HOLES = N_HOLES + 1
232      If (IOS > 0) Call SHUTDOWN('ERROR: Problem with EVHO line')
      End Do COUNT_EVHO_LOOP
231   Rewind(LU_INPUT)
      !
      ! Determine total number of EVSS lines in the input file
      !
      N_SSTANDS = 0
      COUNT_EVSS_LOOP: Do
         Call CHECKREAD('EVSS',LU_INPUT,IOS) 
         If (IOS == 1) Then
            Exit COUNT_EVSS_LOOP
         End If
         Read(LU_INPUT,NML=EVSS,End=233,ERR=234,IOSTAT=IOS)
         N_SSTANDS = N_SSTANDS + 1
234      If (IOS > 0) Call SHUTDOWN('ERROR: Problem with EVSS line')
      End Do COUNT_EVSS_LOOP
233   Rewind(LU_INPUT)
      !
      ! Determine total number of STRS lines in the input file
      !
      N_STRS = 0
      COUNT_STRS_LOOP: Do
         Call CHECKREAD('STRS',LU_INPUT,IOS) 
         If (IOS == 1) Then
            Exit COUNT_STRS_LOOP
         End If
         Read(LU_INPUT,NML=STRS,End=235,ERR=236,IOSTAT=IOS)
         N_STRS = N_STRS + 1
236      If (IOS > 0) Call SHUTDOWN('ERROR: Problem with STRS line')
      End Do COUNT_STRS_LOOP
235   Rewind(LU_INPUT)

      Select Case (COLOR_METHOD)
      Case (-1)
         EVAC_AVATAR_NCOLOR = 1
      Case (0,3,4)
         EVAC_AVATAR_NCOLOR = EVAC_AVATAR_NCOLOR + 1
      Case (1,2,5)
         EVAC_AVATAR_NCOLOR = 7
      Case Default
         EVAC_AVATAR_NCOLOR = 1
      End Select

      ! Allocate avatar color array for Smokeview file write
      EVAC_AVATAR_NCOLOR = Max(1,EVAC_AVATAR_NCOLOR)
      Allocate(EVAC_AVATAR_RGB(3,EVAC_AVATAR_NCOLOR),STAT=IZERO)
      Call ChkMemErr('READ_EVAC','EVAC_AVATAR_RGB',IZERO)

      Select Case (COLOR_METHOD)
      Case (-1)
         EVAC_AVATAR_RGB(1:3,1) = (/ 39, 64,139/)  ! ROYAL BLUE 4
      Case (3)
         EVAC_AVATAR_RGB(1:3,1) = (/ 39, 64,139/)  ! ROYAL BLUE 4
         EVAC_AVATAR_RGB(1:3,EVAC_AVATAR_NCOLOR) = DEAD_RGB
      Case (0,4)
         EVAC_AVATAR_RGB(1:3,EVAC_AVATAR_NCOLOR) = DEAD_RGB
      Case (1,2,5)
         EVAC_AVATAR_RGB(1:3,1) = (/  0,  0,  0/)  ! black
         EVAC_AVATAR_RGB(1:3,2) = (/255,255,  0/)  ! yellow
         EVAC_AVATAR_RGB(1:3,3) = (/  0,  0,255/)  ! blue
         EVAC_AVATAR_RGB(1:3,4) = (/255,  0,  0/)  ! red
         EVAC_AVATAR_RGB(1:3,5) = (/  0,255,  0/)  ! green
         EVAC_AVATAR_RGB(1:3,6) = (/255,  0,255/)  ! magenta
         EVAC_AVATAR_RGB(1:3,7) = DEAD_RGB
      Case Default
         EVAC_AVATAR_RGB(1:3,1) = (/ 39, 64,139/)  ! ROYAL BLUE 4
      End Select
      !
      ! Allocate quantities for EVAC, PERS, EXIT types
      !
      EVAC_PROC_IF: If (MYID==Max(0,EVAC_PROCESS)) Then
         If (npc_evac > 0 ) Then
            Allocate(EVACUATION(NPC_EVAC),STAT=IZERO)
            Call ChkMemErr('READ','EVACUATION',IZERO)
         !Else
         !   Allocate(EVACUATION(1),STAT=IZERO)
         !   Call ChkMemErr('READ','EVACUATION',IZERO)
         End If

         If (n_holes > 0 ) Then
            Allocate(EVAC_HOLES(N_HOLES),STAT=IZERO)
            Call ChkMemErr('READ','EVAC_HOLES',IZERO)
         !Else 
         !   Allocate(EVAC_HOLES(1),STAT=IZERO)
         !   Call ChkMemErr('READ','EVAC_HOLES',IZERO)
         End If

         If (n_sstands > 0 ) Then
            Allocate(EVAC_SSTANDS(N_SSTANDS),STAT=IZERO)
            Call ChkMemErr('READ','EVAC_SSTANDS',IZERO)
         !Else
         !   Allocate(EVAC_SSTANDS(1),STAT=IZERO)
         !   Call ChkMemErr('READ','EVAC_SSTANDS',IZERO)
         End If

         If (n_exits > 0 ) Then
            Allocate(EVAC_EXITS(N_EXITS),STAT=IZERO)
            Call ChkMemErr('READ','EVAC_EXITS',IZERO) 
         !Else
         !   Allocate(EVAC_EXITS(1),STAT=IZERO)
         !   Call ChkMemErr('READ','EVAC_EXITS',IZERO) 
         End If

         If (n_doors > 0 ) Then
            Allocate(EVAC_DOORS(N_DOORS),STAT=IZERO)
            Call ChkMemErr('READ','EVAC_DOORS',IZERO) 
         !Else
         !   Allocate(EVAC_DOORS(1),STAT=IZERO)
         !   Call ChkMemErr('READ','EVAC_DOORS',IZERO) 
         End If

         If (n_entrys > 0 ) Then
            Allocate(EVAC_ENTRYS(N_ENTRYS),STAT=IZERO)
            Call ChkMemErr('READ','EVAC_ENTRYS',IZERO)
         !Else
         !   Allocate(EVAC_ENTRYS(1),STAT=IZERO)
         !   Call ChkMemErr('READ','EVAC_ENTRYS',IZERO)
         End If

         If (n_corrs > 0 ) Then
            Allocate(EVAC_CORRS(N_CORRS),STAT=IZERO)
            Call ChkMemErr('READ','EVAC_CORRS',IZERO)
         !Else
         !   Allocate(EVAC_CORRS(1),STAT=IZERO)
         !   Call ChkMemErr('READ','EVAC_CORRS',IZERO)
         End If

         If (N_STRS > 0 ) Then
            Allocate(EVAC_STRS(N_STRS),STAT=IZERO)
            Call ChkMemErr('READ','EVAC_STRS',IZERO)
         End If

         Allocate(EVAC_PERSON_CLASSES(0:NPC_PERS),STAT=IZERO)
         Call ChkMemErr('READ','EVAC_PERSON_CLASSES',IZERO) 

         n_egrids = 0
         Do n = 1, nmeshes
            If (evacuation_only(n) .And. evacuation_grid(n) ) Then
               n_egrids = n_egrids + 1
            End If
         End  Do

         n_nodes = n_entrys + n_exits + n_doors + n_corrs + n_egrids + n_strs
         If (n_nodes > 0 ) Then
            Allocate(EVAC_Node_List(1:n_nodes),STAT=IZERO)
            Call ChkMemErr('READ','EVAC_NODE_LIST',IZERO) 
         End If

         If (npc_evac > 0 ) Then
            !          EVACUATION(1:NPC_EVAC)%COLOR_INDEX = 1
            EVACUATION(1:NPC_EVAC)%GRID_NAME   = 'null'
            EVACUATION(1:NPC_EVAC)%CLASS_NAME  = 'null'
            EVACUATION(1:NPC_EVAC)%IMESH       = 0
            EVACUATION(1:NPC_EVAC)%ID          = 'null'
         End If
         If (n_holes > 0 ) Then
            EVAC_HOLES(1:N_HOLES)%GRID_NAME   = 'null'
            EVAC_HOLES(1:N_HOLES)%PERS_ID     = 'null'
            EVAC_HOLES(1:N_HOLES)%EVAC_ID     = 'null'
            EVAC_HOLES(1:N_HOLES)%IMESH       = 0
            EVAC_HOLES(1:N_HOLES)%ID          = 'null'
         End If

         EVAC_PERSON_CLASSES(0:NPC_PERS)%ID = 'null'

         If (n_exits > 0 ) Then
            EVAC_EXITS(1:N_EXITS)%ID        = 'null'
            EVAC_EXITS(1:N_EXITS)%TO_NODE   = 'null'
            EVAC_EXITS(1:N_EXITS)%GRID_NAME = 'null'
            EVAC_EXITS(1:N_EXITS)%IMESH     = 0
            !          EVAC_EXITS(1:N_EXITS)%COLOR_INDEX = 1
         End If

         If (n_doors > 0 ) Then
            EVAC_DOORS(1:N_DOORS)%ID        = 'null'
            EVAC_DOORS(1:N_DOORS)%TO_NODE   = 'null'
            EVAC_DOORS(1:N_DOORS)%GRID_NAME = 'null'
            EVAC_DOORS(1:N_DOORS)%IMESH     = 0
            EVAC_DOORS(1:N_DOORS)%IMESH2    = 0
            !          EVAC_DOORS(1:N_DOORS)%COLOR_INDEX = 1
         End If

         If (n_corrs > 0 ) Then
            EVAC_CORRS(1:N_CORRS)%ID        = 'null'
            EVAC_CORRS(1:N_CORRS)%TO_NODE   = 'null'
            EVAC_CORRS(1:N_CORRS)%GRID_NAME = 'null'
            EVAC_CORRS(1:N_CORRS)%IMESH     = 0
            EVAC_CORRS(1:N_CORRS)%IMESH2    = 0
         End If

         If (n_entrys > 0 ) Then
            EVAC_ENTRYS(1:N_ENTRYS)%ID          = 'null'
            EVAC_ENTRYS(1:N_ENTRYS)%TO_NODE     = 'null'
            EVAC_ENTRYS(1:N_ENTRYS)%GRID_NAME   = 'null'
            EVAC_ENTRYS(1:N_ENTRYS)%CLASS_NAME  = 'null'
            EVAC_ENTRYS(1:N_ENTRYS)%IMESH       = 0
            !          EVAC_ENTRYS(1:N_ENTRYS)%COLOR_INDEX = 1
         End If

      End If EVAC_PROC_IF

    End Subroutine COUNT_EVAC_NODES

    Subroutine READ_PERS
      Implicit None
      !
      ! Local variables
      Type (EVAC_PERS_Type), Pointer :: PCP=>NULL()
      !
      ! NEXT PARAMETERS ARE SAME FOR ALL HUMANS. THE LAST
      ! VALUES READ IN FROM 'PERS' LINES ARE VALID.
      FAC_A_WALL  = 1.0_EB
      FAC_B_WALL  = 0.5_EB
      LAMBDA_WALL = 0.2_EB
      NOISEME     = 0.0_EB
      NOISETH     = 0.01_EB
      NOISECM     = 3.0_EB
      NOT_RANDOM = .False.
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
      
      OUTPUT_SPEED         = .False.
      OUTPUT_MOTIVE_FORCE  = .False.
      OUTPUT_MOTIVE_ANGLE  = .False.
      OUTPUT_FED           = .False.
      OUTPUT_OMEGA         = .False.
      OUTPUT_ANGLE         = .False.
      OUTPUT_CONTACT_FORCE = .False.
      OUTPUT_TOTAL_FORCE   = .False.
      DEAD_COLOR = 'null'
      ! 
      ! Read the PERS lines (no read for default n=0 case)
      !
      READ_PERS_LOOP: Do N=0,NPC_PERS
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
         PRE_HIGH = Huge(PRE_HIGH)
         DET_HIGH = Huge(PRE_HIGH)
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
         !
         ! No read for default values
         If ( N > 0 ) Then
            Call CHECKREAD('PERS',LU_INPUT,IOS)
            If (IOS == 1) Then
               Exit READ_PERS_LOOP
            End If
            Read(LU_INPUT,PERS,End=24,IOSTAT=IOS)

            ! Check if some default human group is given.
            Select Case (Trim(DEFAULT_PROPERTIES))
            Case ('Adult','adult','ADULT')
               If (VELOCITY_DIST < 0) Then
                  VELOCITY_DIST = 1
                  VEL_MEAN = 1.25_EB
                  VEL_PARA = 0.30_EB
                  VEL_LOW  = 0.95_EB
                  VEL_HIGH = 1.55_EB
               End If
               If (DIAMETER_DIST < 0) Then
                  DIAMETER_DIST = 1
                  DIA_MEAN = 0.51_EB
                  DIA_PARA = 0.07_EB
                  DIA_LOW  = 0.44_EB
                  DIA_HIGH = 0.58_EB
                  D_TORSO_MEAN    = 0.30_EB
                  D_SHOULDER_MEAN = 0.19_EB
               End If
               If (TAU_EVAC_DIST < 0) Then
                  TAU_EVAC_DIST = 1
                  TAU_MEAN = 1.00_EB
                  TAU_PARA = 0.10_EB
                  TAU_LOW  = 0.80_EB
                  TAU_HIGH = 1.20_EB
               End If
            Case ('Male','male','MALE')
               If (VELOCITY_DIST < 0) Then
                  VELOCITY_DIST = 1
                  VEL_MEAN = 1.35_EB
                  VEL_PARA = 0.20_EB
                  VEL_LOW  = 1.15_EB
                  VEL_HIGH = 1.55_EB
               End If
               If (DIAMETER_DIST < 0) Then
                  DIAMETER_DIST = 1
                  DIA_MEAN = 0.54_EB
                  DIA_PARA = 0.04_EB
                  DIA_LOW  = 0.50_EB
                  DIA_HIGH = 0.58_EB
                  D_TORSO_MEAN    = 0.32_EB
                  D_SHOULDER_MEAN = 0.20_EB
               End If
               If (TAU_EVAC_DIST < 0) Then
                  TAU_EVAC_DIST = 1
                  TAU_MEAN = 1.00_EB
                  TAU_PARA = 0.10_EB
                  TAU_LOW  = 0.80_EB
                  TAU_HIGH = 1.20_EB
               End If
            Case ('Female','female','FEMALE')
               If (VELOCITY_DIST < 0) Then
                  VELOCITY_DIST = 1
                  VEL_MEAN = 1.15_EB
                  VEL_PARA = 0.20_EB
                  VEL_LOW  = 0.95_EB
                  VEL_HIGH = 1.35_EB
               End If
               If (DIAMETER_DIST < 0) Then
                  DIAMETER_DIST = 1
                  DIA_MEAN = 0.48_EB
                  DIA_PARA = 0.04_EB
                  DIA_LOW  = 0.44_EB
                  DIA_HIGH = 0.52_EB
                  D_TORSO_MEAN    = 0.28_EB
                  D_SHOULDER_MEAN = 0.18_EB
               End If
               If (TAU_EVAC_DIST < 0) Then
                  TAU_EVAC_DIST = 1
                  TAU_MEAN = 1.00_EB
                  TAU_PARA = 0.10_EB
                  TAU_LOW  = 0.80_EB
                  TAU_HIGH = 1.20_EB
               End If
            Case ('Child','child','CHILD')
               If (VELOCITY_DIST < 0) Then
                  VELOCITY_DIST = 1
                  VEL_MEAN = 0.90_EB
                  VEL_PARA = 0.30_EB
                  VEL_LOW  = 0.60_EB
                  VEL_HIGH = 1.20_EB
               End If
               If (DIAMETER_DIST < 0) Then
                  DIAMETER_DIST = 1
                  DIA_MEAN = 0.42_EB
                  DIA_PARA = 0.03_EB
                  DIA_LOW  = 0.39_EB
                  DIA_HIGH = 0.45_EB
                  D_TORSO_MEAN    = 0.24_EB
                  D_SHOULDER_MEAN = 0.14_EB
               End If
               If (TAU_EVAC_DIST < 0) Then
                  TAU_EVAC_DIST = 1
                  TAU_MEAN = 1.00_EB
                  TAU_PARA = 0.10_EB
                  TAU_LOW  = 0.80_EB
                  TAU_HIGH = 1.20_EB
               End If
            Case ('Elderly','elderly','ELDERLY')
               If (VELOCITY_DIST < 0) Then
                  VELOCITY_DIST = 1
                  VEL_MEAN = 0.80_EB
                  VEL_PARA = 0.30_EB
                  VEL_LOW  = 0.50_EB
                  VEL_HIGH = 1.10_EB
               End If
               If (DIAMETER_DIST < 0) Then
                  DIAMETER_DIST = 1
                  DIA_MEAN = 0.50_EB
                  DIA_PARA = 0.04_EB
                  DIA_LOW  = 0.46_EB
                  DIA_HIGH = 0.54_EB
                  D_TORSO_MEAN    = 0.30_EB
                  D_SHOULDER_MEAN = 0.18_EB
               End If
               If (TAU_EVAC_DIST < 0) Then
                  TAU_EVAC_DIST = 1
                  TAU_MEAN = 1.00_EB
                  TAU_PARA = 0.10_EB
                  TAU_LOW  = 0.80_EB
                  TAU_HIGH = 1.20_EB
               End If
            Case ('null')
               ! Do nothing, use the defaults
            Case Default
               Write(MESSAGE,'(A,A,A)') 'ERROR: PERS ',Trim(ID),' problem with DEFAULT_PROPERTIES'
               Call SHUTDOWN(MESSAGE)
            End Select

         End If

         If (PRE_MEAN < 0._EB .Or. PRE_LOW < 0._EB) Then
            Write(MESSAGE,'(A,A,A)') 'ERROR: PERS ',Trim(ID), ' PRE-evacuation time should positive.'
            Call SHUTDOWN(MESSAGE)
         End If
         
         DIAMETER_DIST = Max(0,DIAMETER_DIST)
         VELOCITY_DIST = Max(0,VELOCITY_DIST)
         TAU_EVAC_DIST = Max(0,TAU_EVAC_DIST)
         If (DIA_MEAN < 0.0_EB) Then
            Select Case (DIAMETER_DIST)
               Case (1)  ! Uniform
                  DIA_MEAN = 0.5_EB*(DIA_HIGH+DIA_LOW)
               Case (3)  ! Gamma: mean = alpha*beta
                  DIA_MEAN = DIA_PARA*DIA_PARA2
               Case (6)  ! Beta: mean alpha/(alpha+beta)
                  DIA_MEAN = DIA_PARA/(DIA_PARA*DIA_PARA2)
               Case (8)  ! Weibull (Exp, alpha=1): 
                  ! mean = (1/lambda)*GammaFunc(1 + (1/alpha) )
                  ! median = (1/lambda)* (ln(2))^(1/alpha)
                  ! exp: mean = (1/lambda), median = ln(2)/lambda
                  DIA_MEAN = (1.0_EB/DIA_PARA2)*DCDFLIB_Gamma(1.0_EB+(1.0_EB/DIA_PARA))
               Case (9)  ! Gumbel
                  ! mean = gamma / alpha, gamma = 0.5772156649015328606_EB
                  ! median = -(1/alpha)*ln(ln(2))
                  DIA_MEAN = 0.5772156649015328606_EB / DIA_PARA
               Case Default
                  Continue
            End Select
            DIA_MEAN = 0.51_EB
         End If

         !
         ! Avatar colors, integer RGB(3), e.g., (23,255,0)
         If (DEAD_COLOR /= 'null') Then
            Call COLOR2RGB(DEAD_RGB,DEAD_COLOR)
            EVAC_AVATAR_RGB(1:3,EVAC_AVATAR_NCOLOR) = DEAD_RGB
         End If
         If (Any(AVATAR_RGB < 0) .And. AVATAR_COLOR=='null') AVATAR_COLOR = 'ROYAL BLUE 4'
         If (AVATAR_COLOR /= 'null') Call COLOR2RGB(AVATAR_RGB,AVATAR_COLOR)
         If (COLOR_METHOD_TMP == 3) Then
            i_avatar_color = i_avatar_color + 1
            EVAC_AVATAR_RGB(1:3,i_avatar_color) = AVATAR_RGB
         End If
         !
         If (MYID /= Max(0,EVAC_PROCESS)) Cycle READ_PERS_LOOP
         !
         PCP=>EVAC_PERSON_CLASSES(N)
         !
         ! Colors, integer RGB(3), e.g., (23,255,0)
         If (Any(RGB < 0) .And. COLOR=='null') COLOR = 'ROYAL BLUE 4'
         If (COLOR /= 'null') Call COLOR2RGB(RGB,COLOR)
         PCP%RGB = RGB
         PCP%AVATAR_RGB = AVATAR_RGB
         If (COLOR_METHOD_TMP == 3) PCP%Avatar_Color_Index = i_avatar_color
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
         If (M_INERTIA < 0.0_EB ) Then
            PCP%m_iner = -M_INERTIA*(0.25_EB*PCP%D_mean**2+PCP%r_torso**2)**2 / (0.27_EB**2+0.16_EB**2)**2
         Else
            PCP%m_iner = M_INERTIA  ! kg m2
         End If

         PCP%FAC_V0_UP = FAC_V0_UP
         PCP%FAC_V0_DOWN = FAC_V0_DOWN
         !
      End Do READ_PERS_LOOP
24    Rewind(LU_INPUT)
      COLOR_METHOD = COLOR_METHOD_TMP

      If (GROUP_DENS .Le. 0.01_EB) GROUP_DENS = 0.25_EB
      If (GROUP_DENS .Gt. 3.50_EB) GROUP_DENS = 3.50_EB
      DENS_INIT = Max(GROUP_DENS,DENS_INIT)
      If (TDET_SMOKE_DENS < 0.0_EB) TDET_SMOKE_DENS = Huge(TDET_SMOKE_DENS)

      If(SMOKE_MIN_SPEED < 0.0_EB ) SMOKE_MIN_SPEED = 0.0_EB
      If(SMOKE_MIN_SPEED_VISIBILITY < 0.01_EB) Then
         SMOKE_MIN_SPEED_VISIBILITY = 0.01_EB  ! No divisions by zero
      End If

      If (.Not. NOT_RANDOM .And. MYID==Max(0,EVAC_PROCESS)) Then    ! Initialize the generator randomly
         Call Random_seed(size=size_rnd)
         Allocate(seed_rnd(size_rnd),STAT=IZERO)
         Call ChkMemErr('READ_EVAC','seed_rnd',IZERO)
         Call Date_and_time(values = t_rnd)
         seed_rnd = 100.0_EB*t_rnd(7) + t_rnd(8)/10.0_EB
         Call Random_seed(put=seed_rnd)
         Deallocate(seed_rnd)
      End If

      Select Case (COLOR_METHOD)
      Case (-1:5)
         Continue
      Case (7)
         COLOR_METHOD = -1
         If (MYID==Max(0,EVAC_PROCESS)) Write (LU_ERR,'(A)') &
              ' WARNING: COLOR_METHOD=7 is not defined anymore, the default (-1) is used.'
      Case Default
         Write(MESSAGE,'(A,I3,A)') 'ERROR: READ_EVAC COLOR METHOD',COLOR_METHOD, ' is not defined'
         Call SHUTDOWN(MESSAGE)
      End Select

      If (COLOR_METHOD >= 0) EVAC_N_QUANTITIES = EVAC_N_QUANTITIES + 1
      If (OUTPUT_MOTIVE_FORCE) EVAC_N_QUANTITIES = EVAC_N_QUANTITIES + 1
      If (OUTPUT_MOTIVE_ANGLE) EVAC_N_QUANTITIES = EVAC_N_QUANTITIES + 1
      If (OUTPUT_FED) EVAC_N_QUANTITIES = EVAC_N_QUANTITIES + 1
      If (OUTPUT_SPEED) EVAC_N_QUANTITIES = EVAC_N_QUANTITIES + 1
      If (OUTPUT_OMEGA) EVAC_N_QUANTITIES = EVAC_N_QUANTITIES + 1
      If (OUTPUT_ANGLE) EVAC_N_QUANTITIES = EVAC_N_QUANTITIES + 1
      If (OUTPUT_CONTACT_FORCE) EVAC_N_QUANTITIES = EVAC_N_QUANTITIES + 1
      If (OUTPUT_TOTAL_FORCE) EVAC_N_QUANTITIES = EVAC_N_QUANTITIES + 1

      If (EVAC_N_QUANTITIES > 0) Then
         Allocate(EVAC_QUANTITIES_INDEX(EVAC_N_QUANTITIES),STAT=IZERO)
         Call ChkMemErr('READ','EVAC_QUANTITIES_INDEX',IZERO)

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
         If (COLOR_METHOD >= 0) Then
            EVAC_QUANTITIES_INDEX(n)=247  ! HUMAN_COLOR
            n = n + 1
         End If
         If (OUTPUT_MOTIVE_FORCE) Then
            EVAC_QUANTITIES_INDEX(n)=240
            n = n + 1
         End If
         If (OUTPUT_FED) Then
            EVAC_QUANTITIES_INDEX(n)=241
            n = n + 1
         End If
         If (OUTPUT_SPEED) Then
            EVAC_QUANTITIES_INDEX(n)=242
            n = n + 1
         End If
         If (OUTPUT_OMEGA) Then
            EVAC_QUANTITIES_INDEX(n)=243
            n = n + 1
         End If
         If (OUTPUT_ANGLE) Then
            EVAC_QUANTITIES_INDEX(n)=244
            n = n + 1
         End If
         If (OUTPUT_CONTACT_FORCE) Then
            EVAC_QUANTITIES_INDEX(n)=245
            n = n + 1
         End If
         If (OUTPUT_TOTAL_FORCE) Then
            EVAC_QUANTITIES_INDEX(n)=246
            n = n + 1
         End If
         If (OUTPUT_MOTIVE_ANGLE) Then
            EVAC_QUANTITIES_INDEX(n)=248
            n = n + 1
         End If

         If ( n-1 .Ne. EVAC_N_QUANTITIES ) Then
            Write(MESSAGE,'(A,2I4,A)') 'ERROR: Evac output quantities ',EVAC_N_QUANTITIES,n-1, ' Some bug in the program.'
            Call SHUTDOWN(MESSAGE)
         End If
      End If

    End Subroutine READ_PERS

    Subroutine READ_EXIT
      Implicit None
      !
      ! Read the EXIT lines
      !
      ! Local variables
      Integer nm, i1, i2, j1, j2
      Type (EVAC_EXIT_Type), Pointer :: PEX=>NULL()
      Type (EVAC_STRS_Type),  Pointer :: STRP=>NULL()
      Type (MESH_TYPE), Pointer :: M=>NULL()
      !
      READ_EXIT_LOOP: Do N = 1, N_EXITS
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
         CHECK_FLOW    = .False.
         COUNT_ONLY    = .False.
         SHOW          = .True.
         MAX_FLOW      = 0.0_EB
         WIDTH         = 0.0_EB
         HEIGHT        = 2.0_EB
         TIME_OPEN     = -Huge(TIME_OPEN)
         TIME_CLOSE    = Huge(TIME_CLOSE)
         XYZ(:)        = Huge(XYZ)
         XYZ_SMOKE(:)  = Huge(XYZ_SMOKE)
         COLOR_INDEX   = -1
         !
         Call CHECKREAD('EXIT',LU_INPUT,IOS)
         If (IOS == 1) Then
            Exit READ_EXIT_LOOP
         End If
         Read(LU_INPUT,EXIT,End=26,IOSTAT=IOS)
         !
         ! Old input used COLOR_INDEX, next lines are needed for that
         If (MYID==Max(0,EVAC_PROCESS) .And. COLOR_INDEX.Ne.-1) Write (LU_ERR,'(A,A)') &
              ' WARNING: keyword COLOR_INDEX is replaced by COLOR at EXIT line ',Trim(ID)
         If (COLOR_INDEX == 1) COLOR = 'BLACK'  
         If (COLOR_INDEX == 2) COLOR = 'YELLOW' 
         If (COLOR_INDEX == 3) COLOR = 'BLUE'   
         If (COLOR_INDEX == 4) COLOR = 'RED'    
         If (COLOR_INDEX == 5) COLOR = 'GREEN'  
         If (COLOR_INDEX == 6) COLOR = 'MAGENTA'
         If (COLOR_INDEX == 7) COLOR = 'CYAN'   

         ! Colors, integer RGB(3), e.g., (23,255,0)
         If (Any(RGB < 0) .And. COLOR=='null') COLOR = 'FOREST GREEN'
         If (COLOR /= 'null') Call COLOR2RGB(RGB,COLOR)
         If (COLOR_METHOD == 4 .And. .Not.COUNT_ONLY) Then
            i_avatar_color = i_avatar_color + 1
            EVAC_AVATAR_RGB(1:3,i_avatar_color) = RGB
         End If

         If (MYID /= Max(0,EVAC_PROCESS)) Cycle READ_EXIT_LOOP

         PEX=>EVAC_EXITS(N)

         PEX%RGB = RGB
         If (COLOR_METHOD == 4 .And. .Not.COUNT_ONLY) PEX%Avatar_Color_Index = i_avatar_color

         If (EVAC_MESH /= 'null') Then
            MESH_ID = EVAC_MESH
            If (MYID==Max(0,EVAC_PROCESS)) Write (LU_ERR,'(A,A)') &
                 ' WARNING: keyword EVAC_MESH is replaced by MESH_ID at EXIT line ', Trim(ID)
         End If

         ! Check that the exit is properly specified
 
         Do I=1,5,2
            If (XB(I) > XB(I+1)) Then
               DUMMY   = XB(I)
               XB(I)   = XB(I+1)
               XB(I+1) = DUMMY
            End If
         End Do
         ! 
         ! Check which evacuation mesh
         ii = 0
         PEX_MeshLoop: Do i = 1, nmeshes
            If (evacuation_only(i) .And. evacuation_grid(i)) Then
               If (Is_Within_Bounds(XB(1),XB(2),XB(3),XB(4),XB(5),XB(6),&
                  Meshes(i)%XS,Meshes(i)%XF,Meshes(i)%YS,Meshes(i)%YF,Meshes(i)%ZS,Meshes(i)%ZF, 0._EB, 0._EB, 0._EB)) Then
                  If (Trim(MESH_ID) == 'null' .Or. Trim(MESH_ID) == Trim(MESH_NAME(i))) Then
                     ii = ii + 1
                     PEX%IMESH = i
                     !cc             Exit PEX_MeshLoop
                  End If
               End If
            End If
         End Do PEX_MeshLoop
         If (PEX%IMESH == 0) Then
            Write(MESSAGE,'(A,A,A)') 'ERROR: EXIT line ',Trim(ID), ' problem with IMESH, no mesh found'
            Call SHUTDOWN(MESSAGE)
         End If
         If (ii > 1) Then
            Write(MESSAGE,'(A,A,A)') 'ERROR: EXIT line ',Trim(ID), ' not an unique mesh found '
            Call SHUTDOWN(MESSAGE)
         End If
 
         nm = PEX%IMESH
         XB(5) = Meshes(nm)%ZS 
         XB(6) = Meshes(nm)%ZF 
 
         If (XB(1)/=XB(2) .And. XB(3)/=XB(4)) Then
            Write(MESSAGE,'(A,I4,A)') 'ERROR: EXIT',N,' must be a plane'
            Call SHUTDOWN(MESSAGE)
         Endif

         ! User input
         PEX%X1 = XB(1)
         PEX%X2 = XB(2)
         PEX%Y1 = XB(3)
         PEX%Y2 = XB(4)
         PEX%Z1 = XB(5)
         PEX%Z2 = XB(6)

         ! Move user input to mesh cell boundaries
         XB(1) = Max(XB(1),Meshes(nm)%XS)
         XB(2) = Min(XB(2),Meshes(nm)%XF)
         XB(3) = Max(XB(3),Meshes(nm)%YS)
         XB(4) = Min(XB(4),Meshes(nm)%YF)

         I1 = Nint( GINV(XB(1)-Meshes(nm)%XS,1,nm)*Meshes(nm)%RDXI ) 
         I2 = Nint( GINV(XB(2)-Meshes(nm)%XS,1,nm)*Meshes(nm)%RDXI )
         J1 = Nint( GINV(XB(3)-Meshes(nm)%YS,2,nm)*Meshes(nm)%RDETA) 
         J2 = Nint( GINV(XB(4)-Meshes(nm)%YS,2,nm)*Meshes(nm)%RDETA)

         XB(1) = Meshes(nm)%X(I1)
         XB(2) = Meshes(nm)%X(I2)
         XB(3) = Meshes(nm)%Y(J1)
         XB(4) = Meshes(nm)%Y(J2)
         If ( Abs(XB(1)-PEX%X1)>1.E-4_EB .Or. Abs(XB(2)-PEX%X2)>1.E-4_EB .Or. &
              Abs(XB(3)-PEX%Y1)>1.E-4_EB .Or. Abs(XB(4)-PEX%Y2)>1.E-4_EB ) Then
            Write(LU_ERR,fmt='(a,a,a,a)') ' WARNING: Exit line ',Trim(ID),' XB adjusted to mesh ',Trim(MESH_NAME(nm))
            Write(LU_ERR,fmt='(a,6f12.4)') 'Old XB:', PEX%X1,PEX%X2,PEX%Y1,PEX%Y2,PEX%Z1,PEX%Z2
            Write(LU_ERR,fmt='(a,6f12.4)') 'New XB:', XB(1:6)
         End If

         ! Coordinates are lined up with the mesh.
         PEX%X1 = XB(1)
         PEX%X2 = XB(2)
         PEX%Y1 = XB(3)
         PEX%Y2 = XB(4)
         PEX%Z1 = XB(5)
         PEX%Z2 = XB(6)
         !
         PEX%IOR = IOR
         PEX%ID  = Trim(ID)
         PEX%GRID_NAME  = Trim(FLOW_FIELD_ID)
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
         If (TIME_OPEN>TIME_CLOSE) Then
            PEX%IMODE=-1
         Else
            If (TIME_OPEN>T_BEGIN) PEX%IMODE=+2
         End If
         If (CHECK_FLOW) PEX%Flow_max   = MAX_FLOW
         PEX%COUNT_ONLY = .False.
         If (COUNT_ONLY) PEX%COUNT_ONLY = .True.
         PEX%SHOW = SHOW
         If (COUNT_ONLY  ) PEX%SHOW = .False.
         PEX%HEIGHT = HEIGHT

         !       PEX%COLOR_INDEX = Mod(Max(0,COLOR_INDEX-1),7) + 1 ! 1-7 always

         Select Case (IOR)
         Case (-1,+1)
            If (WIDTH <= 0.0_EB) Then
               PEX%Width = XB(4) - XB(3)
            Else
               PEX%Width = WIDTH
            End If
            PEX%ORIENTATION(1)=Real(Sign(1,IOR),EB)
         Case (-2,+2)
            If (WIDTH <= 0.0_EB) Then
               PEX%Width = XB(2) - XB(1)
            Else
               PEX%Width = WIDTH
            End If
            PEX%ORIENTATION(2)=Real(Sign(1,IOR),EB)
         Case (-3)
            If ( (XB(4)-XB(3)) <= 0.0_EB .Or. (XB(2)-XB(1)) <= 0.0_EB) Then
               Write(MESSAGE,'(A,I4,A)') 'ERROR: EXIT',N,' IOR=-3 but not 3-dim object'
               Call SHUTDOWN(MESSAGE)
            End If
            PEX%ORIENTATION(3)=Real(Sign(1,IOR),EB)
         Case (0)
            If ( (XB(4)-XB(3)) <= 0.0_EB .Or. (XB(2)-XB(1)) <= 0.0_EB) Then
               Write(MESSAGE,'(A,I4,A)') 'ERROR: EXIT',N,' no IOR but not 3-dim object'
               Call SHUTDOWN(MESSAGE)
            End If
         Case Default
            Write(MESSAGE,'(A,I4,A)') 'ERROR: EXIT',N,' problem with IOR'
            Call SHUTDOWN(MESSAGE)
         End Select
         ! 
         ! Check which vent field. If VENT_FFIELD is not found, use
         ! the main evac grid.
         PEX%I_VENT_FFIELD = 0
         PEX_Mesh2Loop: Do i = 1, nmeshes
            If ( evacuation_only(i) .And. (Trim(MESH_NAME(i)) == Trim(PEX%VENT_FFIELD)) ) Then
               If ( (PEX%Z1 >= Meshes(i)%ZS .And. PEX%Z2 <= Meshes(i)%ZF).And. &
                    (PEX%Y1 >= Meshes(i)%YS .And. PEX%Y2 <= Meshes(i)%YF).And. &
                    (PEX%X1 >= Meshes(i)%XS .And. PEX%X2 <= Meshes(i)%XF)) Then
                  PEX%I_VENT_FFIELD = i
                  Exit PEX_Mesh2Loop
               End If
            End If
         End Do PEX_Mesh2Loop
         ! If no vent field is given, then use the main evac grid.
         If (PEX%I_VENT_FFIELD == 0) Then
            PEX%I_VENT_FFIELD = PEX%IMESH
            PEX%VENT_FFIELD = Trim(MESH_NAME(PEX%IMESH))
         End If

         PEX%FED_MESH = 0
         If (XYZ(1) < Huge(XYZ)) Then
            PEX%X = XYZ(1)
            PEX%Y = XYZ(2)
            PEX%Z = 0.5_EB*(XB(5)+XB(6))
         Else
            PEX%X = 0.5_EB*(XB(1)+XB(2))
            PEX%Y = 0.5_EB*(XB(3)+XB(4))
            PEX%Z = 0.5_EB*(XB(5)+XB(6))
         End If
         ! 
         ! Check which evacuation floor
         ! Now there may be overlapping meshes.
         ii = 0
         PEX_Mesh3Loop: Do i = 1, nmeshes
            If (evacuation_only(i) .And. evacuation_grid(i)) Then
               If ( (PEX%Z >= Meshes(i)%ZS .And. PEX%Z <= Meshes(i)%ZF).And. &
                    (PEX%Y >= Meshes(i)%YS .And. PEX%Y <= Meshes(i)%YF).And. &
                    (PEX%X >= Meshes(i)%XS .And. PEX%X <= Meshes(i)%XF)) Then
                  If (PEX%IMESH == i ) ii = ii + 1
               End If
            End If
         End Do PEX_Mesh3Loop
         If (ii == 0) Then
            Write(MESSAGE,'(A,A,A)') 'ERROR: EXIT line ',Trim(PEX%ID), ' problem with XYZ, no mesh found'
            Call SHUTDOWN(MESSAGE)
         End If

         ! PEX%Z is used to plot the door on the correct height in Smokeview.
         PEX%Z = PEX%Z + 0.5_EB*PEX%Height - EVACUATION_Z_OFFSET(PEX%IMESH)

         If (XYZ_SMOKE(1) < Huge(XYZ_SMOKE)) Then
            PEX%Xsmoke = XYZ_SMOKE(1)
            PEX%Ysmoke = XYZ_SMOKE(2)
            PEX%Zsmoke = XYZ_SMOKE(3)
         Else
            PEX%Xsmoke = PEX%X
            PEX%Ysmoke = PEX%Y
            PEX%Zsmoke = 0.5_EB*(XB(5)+XB(6)) - EVACUATION_Z_OFFSET(PEX%IMESH) + HUMAN_SMOKE_HEIGHT
         End If

         ! Check if exit is in Stairs
         PEX%STR_INDX = 0
         PEX%STR_SUB_INDX = 0
         CheckExitStrLoop: Do i = 1, N_STRS
            If (EVAC_STRS(i)%IMESH==PEX%IMESH ) Then
               STRP=>EVAC_STRS(i)
               PEX%STR_INDX = i
               Do j = 1,STRP%N_NODES
                  If ( Is_Within_Bounds(PEX%X1,PEX%X2,PEX%Y1,PEX%Y2,PEX%Z1,PEX%Z2, &
                       STRP%XB_NODE(j,1), STRP%XB_NODE(j,2), STRP%XB_NODE(j,3),STRP%XB_NODE(j,4), &
                       STRP%XB_NODE(j,5), STRP%XB_NODE(j,6), 0._EB, 0._EB, 0._EB)) Then
                     PEX%STR_SUB_INDX = j
                     Exit CheckExitStrLoop
                  End If
               End Do
            End If
         End Do CheckExitStrLoop
         ! 
         ! Check, which fire grid and i,j,k (xyz)
         PEX_SmokeLoop: Do i = 1, nmeshes
            If (.Not. evacuation_only(i)) Then
               If ( Is_Within_Bounds(PEX%Xsmoke,PEX%Xsmoke,PEX%Ysmoke,PEX%Ysmoke,PEX%Zsmoke,PEX%Zsmoke,&
                    Meshes(i)%XS,Meshes(i)%XF,Meshes(i)%YS,Meshes(i)%YF,Meshes(i)%ZS,Meshes(i)%ZF,0._EB,0._EB,0._EB)) Then
                  PEX%FED_MESH = i
                  Exit PEX_SmokeLoop
               End If
            End If
            !     No mesh found
            PEX%FED_MESH = -1
         End Do PEX_SmokeLoop
         !   No mesh found
         If (PEX%FED_MESH == 0) PEX%FED_MESH = -1

         If (PEX%FED_MESH > 0) Then 
            M => MESHES(PEX%FED_MESH)
            II = Floor(M%CELLSI(Floor((PEX%Xsmoke-M%XS)*M%RDXINT))+ 1.0_EB)
            JJ = Floor(M%CELLSJ(Floor((PEX%Ysmoke-M%YS)*M%RDYINT))+ 1.0_EB)
            KK = Floor(M%CELLSK(Floor((PEX%Zsmoke-M%ZS)*M%RDZINT))+ 1.0_EB)
            If ( M%SOLID(M%CELL_INDEX(II,JJ,KK)) ) Then
               PEX%FED_MESH = -1   ! no smoke at a solid object
               PEX%II = 0
               PEX%JJ = 0
               PEX%KK = 0
            Else
               PEX%II = II
               PEX%JJ = JJ
               PEX%KK = KK
            End If
         Else
            PEX%II = 0
            PEX%JJ = 0
            PEX%KK = 0
         End If
         !
      End Do READ_EXIT_LOOP
26    Rewind(LU_INPUT)

    End Subroutine READ_EXIT

    Subroutine READ_DOOR    
      Implicit None
      !
      ! Read the DOOR lines
      !
      ! Local variables
      Type (EVAC_DOOR_Type), Pointer :: PDX=>NULL()
      Type (MESH_TYPE), Pointer :: M=>NULL()
      Integer nm, i1, i2, j1, j2
      !
      READ_DOOR_LOOP: Do N = 1, N_DOORS
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
         CHECK_FLOW    = .False.
         EXIT_SIGN     = .False.
         SHOW          = .True.
         MAX_FLOW      = 0.0_EB
         WIDTH         = 0.0_EB
         HEIGHT        = 2.0_EB
         TIME_OPEN     = -Huge(TIME_OPEN)
         TIME_CLOSE    = Huge(TIME_CLOSE)
         XYZ(:)        = Huge(XYZ)
         XYZ_SMOKE(:)  = Huge(XYZ_SMOKE)
         COLOR_INDEX   = -1
         KEEP_XY       = .False.
         !
         Call CHECKREAD('DOOR',LU_INPUT,IOS)
         If (IOS == 1) Then
            Exit READ_DOOR_LOOP
         End If
         Read(LU_INPUT,DOOR,End=27,IOSTAT=IOS)
         !
         ! Old input used COLOR_INDEX, next lines are needed for that
         If (MYID==Max(0,EVAC_PROCESS) .And. COLOR_INDEX.Ne.-1) Write (LU_ERR,'(A,A)') &
              ' WARNING: keyword COLOR_INDEX is replaced by COLOR at DOOR line ',Trim(ID)
         If (COLOR_INDEX == 1) COLOR = 'BLACK'  
         If (COLOR_INDEX == 2) COLOR = 'YELLOW' 
         If (COLOR_INDEX == 3) COLOR = 'BLUE'   
         If (COLOR_INDEX == 4) COLOR = 'RED'    
         If (COLOR_INDEX == 5) COLOR = 'GREEN'  
         If (COLOR_INDEX == 6) COLOR = 'MAGENTA'
         If (COLOR_INDEX == 7) COLOR = 'CYAN'   
         !
         ! Colors, integer RGB(3), e.g., (23,255,0)
         If (Any(RGB < 0) .And. COLOR=='null') COLOR = 'FOREST GREEN'
         If (COLOR /= 'null') Call COLOR2RGB(RGB,COLOR)
         If (COLOR_METHOD == 4) Then
            i_avatar_color = i_avatar_color + 1
            EVAC_AVATAR_RGB(1:3,i_avatar_color) = RGB
         End If

         If (MYID /= Max(0,EVAC_PROCESS)) Cycle READ_DOOR_LOOP

         PDX=>EVAC_DOORS(N)

         PDX%RGB = RGB
         If (COLOR_METHOD == 4) PDX%Avatar_Color_Index = i_avatar_color

         If (EVAC_MESH /= 'null') Then
            MESH_ID = EVAC_MESH
            If (MYID==Max(0,EVAC_PROCESS)) Write (LU_ERR,'(A,A)') &
                 ' WARNING: keyword EVAC_MESH is replaced by MESH_ID at DOOR line ', Trim(ID)
         End If

         ! Check that the door is properly specified

         Do I=1,5,2
            If (XB(I) > XB(I+1)) Then
               DUMMY   = XB(I)
               XB(I)   = XB(I+1)
               XB(I+1) = DUMMY
            End If
         End Do
         ! 
         ! Check which evacuation floor
         ! Now there may be overlapping meshes.
         ii = 0
         PDX_MeshLoop: Do i = 1, nmeshes
            If (evacuation_only(i) .And. evacuation_grid(i)) Then
               If (Is_Within_Bounds(XB(1),XB(2),XB(3),XB(4),XB(5),XB(6),&
                  Meshes(i)%XS,Meshes(i)%XF,Meshes(i)%YS,Meshes(i)%YF,Meshes(i)%ZS,Meshes(i)%ZF, 0._EB, 0._EB, 0._EB)) Then
                  If (Trim(MESH_ID) == 'null' .Or. Trim(MESH_ID) == Trim(MESH_NAME(i))) Then
                     ii = ii + 1
                     PDX%IMESH = i
                  End If
                  If (Trim(TO_NODE) == Trim(MESH_NAME(i))) Then
                     PDX%IMESH2 = i
                  End If
               End If
            End If
         End Do PDX_MeshLoop
         If (PDX%IMESH == 0) Then
            Write(MESSAGE,'(A,A,A)') 'ERROR: DOOR line ',Trim(ID), ' problem with IMESH, no mesh found'
            Call SHUTDOWN(MESSAGE)
         End If
         If (ii > 1) Then
            Write(MESSAGE,'(A,A,A)') 'ERROR: DOOR line ',Trim(ID), ' not an unique mesh found '
            Call SHUTDOWN(MESSAGE)
         End If

         nm = PDX%IMESH
 
         If (XB(1)/=XB(2) .And. XB(3)/=XB(4)) Then
            Write(MESSAGE,'(A,I4,A)') 'ERROR: DOOR',N,' must be a plane'
            Call SHUTDOWN(MESSAGE)
         Endif

         ! User input
         PDX%X1 = XB(1)
         PDX%X2 = XB(2)
         PDX%Y1 = XB(3)
         PDX%Y2 = XB(4)
         PDX%Z1 = XB(5)
         PDX%Z2 = XB(6)

         ! Move user input to mesh cell boundaries
         XB(1) = Max(XB(1),Meshes(nm)%XS)
         XB(2) = Min(XB(2),Meshes(nm)%XF)
         XB(3) = Max(XB(3),Meshes(nm)%YS)
         XB(4) = Min(XB(4),Meshes(nm)%YF)
         XB(5) = Max(XB(5),Meshes(nm)%ZS)
         XB(6) = Min(XB(6),Meshes(nm)%ZF)

         I1 = Nint( GINV(XB(1)-Meshes(nm)%XS,1,nm)*Meshes(nm)%RDXI ) 
         I2 = Nint( GINV(XB(2)-Meshes(nm)%XS,1,nm)*Meshes(nm)%RDXI )
         J1 = Nint( GINV(XB(3)-Meshes(nm)%YS,2,nm)*Meshes(nm)%RDETA) 
         J2 = Nint( GINV(XB(4)-Meshes(nm)%YS,2,nm)*Meshes(nm)%RDETA)

         XB(1) = Meshes(nm)%X(I1)
         XB(2) = Meshes(nm)%X(I2)
         XB(3) = Meshes(nm)%Y(J1)
         XB(4) = Meshes(nm)%Y(J2)
         If ( Abs(XB(1)-PDX%X1)>1.E-4_EB .Or. Abs(XB(2)-PDX%X2)>1.E-4_EB .Or. &
              Abs(XB(3)-PDX%Y1)>1.E-4_EB .Or. Abs(XB(4)-PDX%Y2)>1.E-4_EB ) Then
            Write(LU_ERR,fmt='(a,a,a,a)') ' WARNING: Door line ',Trim(ID),' XB adjusted to mesh ',Trim(MESH_NAME(nm))
            Write(LU_ERR,fmt='(a,6f12.4)') 'Old XB:', PDX%X1,PDX%X2,PDX%Y1,PDX%Y2,PDX%Z1,PDX%Z2
            Write(LU_ERR,fmt='(a,6f12.4)') 'New XB:', XB(1:6)
         End If

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
         If (TIME_OPEN>TIME_CLOSE) Then
            PDX%IMODE=-1
         Else
            If (TIME_OPEN>T_BEGIN) PDX%IMODE=+2
         End If
         If (CHECK_FLOW) PDX%Flow_max   = MAX_FLOW
         PDX%HEIGHT = HEIGHT

         !       PDX%COLOR_INDEX = Mod(Max(0,COLOR_INDEX-1),7) ! 1-7 always

         Select Case (IOR)
         Case (-1,+1)
            If (WIDTH <= 0.0_EB) Then
               PDX%Width = XB(4) - XB(3)
            Else
               PDX%Width = WIDTH
            End If
            PDX%ORIENTATION(1)=Real(Sign(1,IOR),EB)
         Case (-2,+2)
            If (WIDTH <= 0.0_EB) Then
               PDX%Width = XB(2) - XB(1)
            Else
               PDX%Width = WIDTH
            End If
            PDX%ORIENTATION(2)=Real(Sign(1,IOR),EB)
         Case (-3)
            If ( (XB(4)-XB(3)) <= 0.0_EB .Or. (XB(2)-XB(1)) <= 0.0_EB) Then
               Write(MESSAGE,'(A,I4,A)') 'ERROR: DOOR',N,' IOR=-3 but not 3-dim object'
               Call SHUTDOWN(MESSAGE)
            End If
            PDX%ORIENTATION(3)=Real(Sign(1,IOR),EB)
         Case (0)
            If ( (XB(4)-XB(3)) <= 0.0_EB .Or. (XB(2)-XB(1)) <= 0.0_EB) Then
               Write(MESSAGE,'(A,I4,A)') 'ERROR: DOOR',N,' no IOR but not 3-dim object'
               Call SHUTDOWN(MESSAGE)
            End If
         Case Default
            Write(MESSAGE,'(A,I4,A)') 'ERROR: DOOR',N,' problem with IOR'
            Call SHUTDOWN(MESSAGE)
         End Select
         ! 
         ! Check which vent field. If VENT_FFIELD is not found, use
         ! the main evac grid.
         PDX%I_VENT_FFIELD = 0
         PDX_Mesh2Loop: Do i = 1, nmeshes
            If ( evacuation_only(i) .And. (Trim(MESH_NAME(i)) == Trim(PDX%VENT_FFIELD)) ) Then
               If ( (PDX%Z1 >= Meshes(i)%ZS .And. PDX%Z2 <= Meshes(i)%ZF).And. &
                    (PDX%Y1 >= Meshes(i)%YS .And. PDX%Y2 <= Meshes(i)%YF).And. &
                    (PDX%X1 >= Meshes(i)%XS .And. PDX%X2 <= Meshes(i)%XF)) Then
                  PDX%I_VENT_FFIELD = i
                  Exit PDX_Mesh2Loop
               End If
            End If
         End Do PDX_Mesh2Loop
         ! If no vent field is given, then use the main evac grid.
         If (PDX%I_VENT_FFIELD == 0) Then
            PDX%I_VENT_FFIELD = PDX%IMESH
            PDX%VENT_FFIELD = Trim(MESH_NAME(PDX%IMESH))
         End If

         PDX%FED_MESH = 0
         If (XYZ(1) < Huge(XYZ)) Then
            PDX%X = XYZ(1)
            PDX%Y = XYZ(2)
            PDX%Z = 0.5_EB*(XB(5)+XB(6))
         Else
            PDX%X = 0.5_EB*(XB(1)+XB(2))
            PDX%Y = 0.5_EB*(XB(3)+XB(4))
            PDX%Z = 0.5_EB*(XB(5)+XB(6))
         End If

         If (XYZ_SMOKE(1) < Huge(XYZ_SMOKE)) Then
            PDX%Xsmoke = XYZ_SMOKE(1)
            PDX%Ysmoke = XYZ_SMOKE(2)
            PDX%Zsmoke = XYZ_SMOKE(3)
         Else
            PDX%Xsmoke = PDX%X
            PDX%Ysmoke = PDX%Y
            PDX%Zsmoke = 0.5_EB*(XB(5)+XB(6)) - EVACUATION_Z_OFFSET(PDX%IMESH) + HUMAN_SMOKE_HEIGHT
         End If

         ! Check which evacuation floor
         ii = 0
         PDX_Mesh3Loop: Do i = 1, nmeshes
            If (evacuation_only(i) .And. evacuation_grid(i)) Then
               If ( (PDX%Z >= Meshes(i)%ZS .And. PDX%Z <= Meshes(i)%ZF).And. &
                    (PDX%Y >= Meshes(i)%YS .And. PDX%Y <= Meshes(i)%YF).And. &
                    (PDX%X >= Meshes(i)%XS .And. PDX%X <= Meshes(i)%XF)) Then
                  If (PDX%IMESH == i ) ii = ii + 1
               End If
            End If
         End Do PDX_Mesh3Loop
         If (ii == 0) Then
            Write(MESSAGE,'(A,A,A)') 'ERROR: DOOR line ',Trim(PDX%ID), ' problem with XYZ, no mesh found'
            Call SHUTDOWN(MESSAGE)
         End If

         ! PDX%Z is used to plot the door on the correct height in Smokeview.
         PDX%Z = PDX%Z + 0.5_EB*PDX%Height - EVACUATION_Z_OFFSET(PDX%IMESH)

         ! Check, which fire grid and i,j,k (xyz)
         PDX_SmokeLoop: Do i = 1, nmeshes
            If (.Not. evacuation_only(i)) Then
               If ( (PDX%Zsmoke >= Meshes(i)%ZS .And. PDX%Zsmoke <= Meshes(i)%ZF) .And. &
                    (PDX%Ysmoke >= Meshes(i)%YS .And. PDX%Ysmoke <= Meshes(i)%YF) .And. &
                    (PDX%Xsmoke >= Meshes(i)%XS .And. PDX%Xsmoke <= Meshes(i)%XF)) Then
                  PDX%FED_MESH = i
                  Exit PDX_SmokeLoop
               End If
            End If
            !     No mesh found
            PDX%FED_MESH = -1
         End Do PDX_SmokeLoop
         !   No mesh found
         If (PDX%FED_MESH == 0) PDX%FED_MESH = -1

         If (PDX%FED_MESH > 0) Then 
            M => MESHES(PDX%FED_MESH)
            II = Floor(M%CELLSI(Floor((PDX%Xsmoke-M%XS)*M%RDXINT))+ 1.0_EB)
            JJ = Floor(M%CELLSJ(Floor((PDX%Ysmoke-M%YS)*M%RDYINT))+ 1.0_EB)
            KK = Floor(M%CELLSK(Floor((PDX%Zsmoke-M%ZS)*M%RDZINT))+ 1.0_EB)
            If ( M%SOLID(M%CELL_INDEX(II,JJ,KK)) ) Then
               PDX%FED_MESH = -1   ! no smoke at a solid object
               PDX%II = 0
               PDX%JJ = 0
               PDX%KK = 0
            Else
               PDX%II = II
               PDX%JJ = JJ
               PDX%KK = KK
            End If
         Else
            PDX%II = 0
            PDX%JJ = 0
            PDX%KK = 0
         End If
         ! 
      End Do READ_DOOR_LOOP
27    Rewind(LU_INPUT)

    End Subroutine READ_DOOR

    Subroutine READ_CORR
      Implicit None
      !
      ! Local variables
      Type (EVAC_CORR_Type), Pointer :: PCX=>NULL()
      Type (MESH_TYPE), Pointer :: M=>NULL()
      !
      ! Read the CORR line
      !
      n_max_in_corrs = 0
      READ_CORR_LOOP: Do N = 1, N_CORRS
         If (MYID /= Max(0,EVAC_PROCESS)) Cycle READ_CORR_LOOP
         PCX=>EVAC_CORRS(N)
         !
         ID            = 'null'
         RGB           = -1
         COLOR         = 'null'
         XB            = Huge(XB)
         XB1           = Huge(XB1)
         XB2           = Huge(XB2)
         IOR           = 0
         FLOW_FIELD_ID = 'null'
         TO_NODE       = 'null'
         CHECK_FLOW    = .False.
         MAX_FLOW      = 0.0_EB
         WIDTH         = 0.0_EB
         WIDTH1        = 0.0_EB
         WIDTH2        = 0.0_EB
         FAC_SPEED     = 0.0_EB
         EFF_WIDTH     = 0.0_EB
         EFF_LENGTH    = 0.0_EB
         MAX_HUMANS_INSIDE = 0
         !
         Call CHECKREAD('CORR',LU_INPUT,IOS)
         If (IOS == 1) Then
            Exit READ_CORR_LOOP
         End If
         Read(LU_INPUT,CORR,End=29,IOSTAT=IOS)
         !
         !
         Do I=1,5,2
            If (XB(I) > XB(I+1)) Then
               DUMMY   = XB(I)
               XB(I)   = XB(I+1)
               XB(I+1) = DUMMY
            End If
         End Do
         Do I=1,5,2
            If (XB1(I) > XB1(I+1)) Then
               DUMMY   = XB1(I)
               XB1(I)   = XB1(I+1)
               XB1(I+1) = DUMMY
            End If
         End Do
         Do I=1,5,2
            If (XB2(I) > XB2(I+1)) Then
               DUMMY   = XB2(I)
               XB2(I)   = XB2(I+1)
               XB2(I+1) = DUMMY
            End If
         End Do
         !
         ! Position, where smoke etc. is saved.
         ! If both XB and XB1 are given, use XB1
         If ( XB(1) < Huge(XB) ) Then
            PCX%FED_MESH = 0
            PCX%X1 = 0.5_EB*( XB(1) +  XB(2))
            PCX%Y1 = 0.5_EB*( XB(3) +  XB(4))
            PCX%Z1 = 0.5_EB*( XB(5) +  XB(6))
         Else
            PCX%FED_MESH = -1
            PCX%X1 = 0.0_EB
            PCX%Y1 = 0.0_EB
            PCX%Z1 = 0.0_EB
         End If
         If ( XB1(1) < Huge(XB1) ) Then
            PCX%FED_MESH = 0
            PCX%X1 = 0.5_EB*( XB1(1) +  XB1(2))
            PCX%Y1 = 0.5_EB*( XB1(3) +  XB1(4))
            PCX%Z1 = 0.5_EB*( XB1(5) +  XB1(6))
         Else If (XB(1) == Huge(XB) ) Then
            PCX%FED_MESH = -1
            PCX%X1 = 0.0_EB
            PCX%Y1 = 0.0_EB
            PCX%Z1 = 0.0_EB
         End If
         If ( XB2(1) < Huge(XB2) ) Then
            PCX%FED_MESH2 = 0
            PCX%X2 = 0.5_EB*(XB2(1) + XB2(2))
            PCX%Y2 = 0.5_EB*(XB2(3) + XB2(4))
            PCX%Z2 = 0.5_EB*(XB2(5) + XB2(6))
         Else
            PCX%FED_MESH2 = -1
            PCX%X2 = 0.0_EB
            PCX%Y2 = 0.0_EB
            PCX%Z2 = 0.0_EB
         End If

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
         If (MAX_HUMANS_INSIDE > 0 ) Then
            PCX%MAX_HUMANS_INSIDE = MAX_HUMANS_INSIDE
         Else
            Write(MESSAGE,'(A,I4,A)') 'ERROR: CORR',N,' MAX_HUMANS_INSIDE <= 0'
            Call SHUTDOWN(MESSAGE)
         End If

         If (FAC_SPEED < 0 ) Then
            Write(MESSAGE,'(A,I4,A)') 'ERROR: CORR',N,' FAC_SPEED < 0'
            Call SHUTDOWN(MESSAGE)
         Else
            If (FAC_SPEED == 0.0_EB) FAC_SPEED = 0.6_EB
            PCX%Fac_Speed = FAC_SPEED
         End If

         PCX%Flow_max   = 0.0_EB
         If (CHECK_FLOW) PCX%Flow_max   = MAX_FLOW

         PCX%Width = Max( Abs(XB(4)-XB(3)) , Abs(XB(2)-XB(1)) )

         PCX%Eff_Width = 0.0_EB
         If (EFF_WIDTH > 0.0_EB ) Then
            PCX%Eff_Width = EFF_WIDTH
         Else
            PCX%Eff_Width = PCX%Width
         End If

         PCX%Eff_Length = 0.0_EB
         If (EFF_LENGTH > 0.0_EB ) Then
            PCX%Eff_Length = EFF_LENGTH
         Else
            Write(MESSAGE,'(A,I4,A)') 'ERROR: CORR',N,' EFF_LENGTH <= 0'
            Call SHUTDOWN(MESSAGE)
         End If
         PCX%Eff_Area = PCX%Eff_Length*PCX%Eff_Width

         PCX%Width1 = WIDTH1
         PCX%Width2 = WIDTH2
         If (WIDTH1*WIDTH2 <= 0.0_EB) Then
            PCX%Width1 = PCX%Width
            PCX%Width2 = PCX%Width
         End If
         ! 
         PCX_MeshLoop: Do i = 1, nmeshes
            If (.Not. evacuation_only(i) .And. PCX%FED_MESH >= 0) Then
               If ( (PCX%Z1 >= Meshes(i)%ZS .And. PCX%Z1 <= Meshes(i)%ZF).And. &
                    (PCX%Y1 >= Meshes(i)%YS .And. PCX%Y1 <= Meshes(i)%YF).And. &
                    (PCX%X1 >= Meshes(i)%XS .And. PCX%X1 <= Meshes(i)%XF)) Then
                  PCX%FED_MESH = i
                  Exit PCX_MeshLoop
               End If
            End If
            !     No mesh found
            PCX%FED_MESH = -1
         End Do PCX_MeshLoop
         !   No mesh found
         If (PCX%FED_MESH == 0) PCX%FED_MESH = -1

         If (PCX%FED_MESH > 0) Then 
            M => MESHES(PCX%FED_MESH)
            II = Floor( M%CELLSI(Floor((PCX%X1-M%XS)*M%RDXINT)) + 1.0_EB  )
            JJ = Floor( M%CELLSJ(Floor((PCX%Y1-M%YS)*M%RDYINT)) + 1.0_EB  )
            KK = Floor( M%CELLSK(Floor((PCX%Z1-M%ZS)*M%RDZINT)) + 1.0_EB  )
            If ( M%SOLID(M%CELL_INDEX(II,JJ,KK)) ) Then
               PCX%FED_MESH = -1   ! no smoke at a solid object
               PCX%II(1) = 0
               PCX%JJ(1) = 0
               PCX%KK(1) = 0
            Else
               PCX%II(1) = II
               PCX%JJ(1) = JJ
               PCX%KK(1) = KK
            End If
         Else
            PCX%II(1) = 0
            PCX%JJ(1) = 0
            PCX%KK(1) = 0
         End If

         PCX_MeshLoop2: Do i = 1, nmeshes
            If (.Not. evacuation_only(i) .And. PCX%FED_MESH2 >= 0) Then
               If ( (PCX%Z2 >= Meshes(i)%ZS .And. PCX%Z2 <= Meshes(i)%ZF).And. &
                    (PCX%Y2 >= Meshes(i)%YS .And. PCX%Y2 <= Meshes(i)%YF).And. &
                    (PCX%X2 >= Meshes(i)%XS .And. PCX%X2 <= Meshes(i)%XF)) Then
                  PCX%FED_MESH2 = i
                  Exit PCX_MeshLoop2
               End If
            End If
         End Do PCX_MeshLoop2
         !   No mesh found
         If (PCX%FED_MESH2 == 0) PCX%FED_MESH2 = -1

         If (PCX%FED_MESH2 > 0) Then 
            M => MESHES(PCX%FED_MESH2)
            II = Floor( M%CELLSI(Floor((PCX%X2-M%XS)*M%RDXINT)) + 1.0_EB  )
            JJ = Floor( M%CELLSJ(Floor((PCX%Y2-M%YS)*M%RDYINT)) + 1.0_EB  )
            KK = Floor( M%CELLSK(Floor((PCX%Z2-M%ZS)*M%RDZINT)) + 1.0_EB  )
            If ( M%SOLID(M%CELL_INDEX(II,JJ,KK)) ) Then
               PCX%FED_MESH2 = -1   ! no smoke at a solid object
               PCX%II(2) = 0
               PCX%JJ(2) = 0
               PCX%KK(2) = 0
            Else
               PCX%II(2) = II
               PCX%JJ(2) = JJ
               PCX%KK(2) = KK
            End If
         Else
            PCX%II(2) = 0
            PCX%JJ(2) = 0
            PCX%KK(2) = 0
         End If
         ! 
         n_max_in_corrs = Max(n_max_in_corrs,PCX%MAX_HUMANS_INSIDE)

         ! Initialize the linked lists of persons who are inside corridors.
         PCX%n_inside = 0
         Nullify(PCX%First)

         ! Colors, integer RGB(3), e.g., (23,255,0)
         If (Any(RGB < 0) .And. COLOR=='null') COLOR = 'BLACK'
         If (COLOR /= 'null') Call COLOR2RGB(RGB,COLOR)
         PCX%RGB = RGB
         !
      End Do READ_CORR_LOOP
29    Rewind(LU_INPUT)

    End Subroutine READ_CORR

    Subroutine READ_STRS
      Implicit None
      !
      ! Local variables
      REAL(EB) Z_TMP
      Type (EVAC_STRS_Type), Pointer :: STRP=>NULL()

      ! Read the STRS line
      READ_STRS_LOOP: Do N = 1,N_STRS
         STRP=>EVAC_STRS(N)
         !
         ID                          = 'null'
         XB                          = 0._EB
         XB_CORE                     = 0._EB
         XB_CORES                    = 0._EB
         XB_LANDINGS                 = 0._EB
         RIGHT_HANDED                = .True.
         LEFT_HANDED                 = .False.
         MESH_ID                     = 'null'
         N_LANDINGS                  = 0
         VERTICAL_LANDING_SEPARATION = 0._EB
         FAC_V0_UP                   = 1.0_EB
         FAC_V0_DOWN                 = 1.0_EB
         FAC_V0_HORI                 = 1.0_EB
         !
         Call CHECKREAD('STRS',LU_INPUT,IOS)
         If (IOS == 1) Then
            Exit READ_STRS_LOOP
         End If
         Read(LU_INPUT,STRS,End=32,IOSTAT=IOS)
         !
         Do I=1,5,2
            If (XB(I) > XB(I+1)) Then
               DUMMY   = XB(I)
               XB(I)   = XB(I+1)
               XB(I+1) = DUMMY
            End If
         End Do
         !
         STRP%ID          = ID
         STRP%XB          = XB
         ii = 0
         STRP_MeshLoop: Do I = 1, NMESHES
            If (.Not. evacuation_only(I)) Cycle
            If (.Not. evacuation_grid(I)) Cycle
            If (Trim(MESH_ID) == 'null' .Or. Trim(MESH_ID)==Trim(MESH_NAME(I))) Then
               ii = ii + 1
               STRP%IMESH = I
               Exit STRP_MeshLoop
            End If
         End Do STRP_MeshLoop

         ! Count number of cores
         STRP%N_CORES = 0
         DO I = 1,500
            If (ANY(XB_CORES(I,:)/=0._EB)) STRP%N_CORES = STRP%N_CORES + 1
         ENDDO

         STRP%INODE          = 0
         If (LEFT_HANDED) RIGHT_HANDED = .FALSE.
         STRP%RIGHT_HANDED = RIGHT_HANDED
         STRP%MESH_ID     = MESH_ID
         STRP%FAC_V0_UP   = FAC_V0_UP
         STRP%FAC_V0_DOWN = FAC_V0_DOWN
         STRP%FAC_V0_HORI = FAC_V0_HORI

         If (N_LANDINGS>500) Then
            Write(MESSAGE,'(A,I4,A)') 'ERROR: STRS',N,' N_LANDINGS > 500'
            Call SHUTDOWN(MESSAGE)
         End If
         STRP%N_LANDINGS = N_LANDINGS
         STRP%N_NODES = 2*N_LANDINGS - 1

         ! Allocate landing and stair geometry
         Allocate(STRP%XB_NODE(1:STRP%N_NODES,1:8), STAT=IZERO)
         Call ChkMemErr('Read_Evac','STRP%XB_NODE',IZERO)
         Allocate(STRP%NODE_IOR(1:STRP%N_NODES), STAT=IZERO)
         Call ChkMemErr('Read_Evac','STRP%NODE_IOR',IZERO)
         Allocate(STRP%NODE_TYPE(1:STRP%N_NODES), STAT=IZERO)
         Call ChkMemErr('Read_Evac','STRP%NODE_TYPE',IZERO)
         STRP%XB_NODE = 0._EB 
         STRP%NODE_IOR = 0

         Allocate(STRP%XB_CORE(1:MAX(1,STRP%N_CORES),1:6), STAT=IZERO)
         Call ChkMemErr('Read_Evac','STRP%XB_CORE',IZERO)
         If (STRP%N_CORES == 0) Then
            STRP%N_CORES = 1
            STRP%XB_CORE(1,1:4) = XB_CORE(1:4)
            STRP%XB_CORE(1,5:6) = XB(5:6)
            If (All(XB_CORE==0._EB)) Then
                Write(MESSAGE,'(3A)') 'ERROR: STRS object ', Trim(ID), ' has no XB_CORE defined.'
                Call SHUTDOWN(MESSAGE)
            Endif
         
         Else
            STRP%XB_CORE(1:STRP%N_CORES,1:6) = XB_CORES(1:STRP%N_CORES,1:6)
         Endif
         Allocate(STRP%I_CORE(1:STRP%N_NODES), STAT=IZERO)
         Call ChkMemErr('Read_Evac','STRP%I_CORE',IZERO)

         ! Count and copy explicitly given landing geometries
         NL = 0
         Do I = 1,N_LANDINGS
            If (Any(XB_LANDINGS(I,:)/=0._EB)) NL = NL + 1
         End Do
         Do I = NL+1,N_LANDINGS
            XB_LANDINGS(I,1:4) = XB_LANDINGS(I-2,1:4)
            XB_LANDINGS(I,5:6) = XB_LANDINGS(I-1,5:6)+VERTICAL_LANDING_SEPARATION
         End Do
         ! Compute stair geometry
         Do I = 1,N_LANDINGS-1
            XB_STAIRS(I,5) = 0.5_EB*(XB_LANDINGS(I,5)+XB_LANDINGS(I,6))
            XB_STAIRS(I,6) = 0.5_EB*(XB_LANDINGS(I+1,5)+XB_LANDINGS(I+1,6))
            STR_Height = XB_STAIRS(I,6)-XB_STAIRS(I,5)
            Z_TMP = 0.5_EB * (XB_STAIRS(I,6)+XB_STAIRS(I,5))
            ! Choose core
            Do J = 1,STRP%N_CORES
               If ((Z_TMP >= STRP%XB_CORE(J,5)) .AND. (Z_TMP <= STRP%XB_CORE(J,6))) Exit
            Enddo
            J = Min(J,STRP%N_CORES)
            If (XB_LANDINGS(I+1,1)>XB_LANDINGS(I,2)) Then ! From -x to +x
               XB_STAIRS(I,1) = XB_LANDINGS(I,2)
               XB_STAIRS(I,2) = XB_LANDINGS(I+1,1)
               STR_Length = XB_STAIRS(I,2)-XB_STAIRS(I,1)
               If (RIGHT_HANDED) Then 
                  XB_STAIRS(I,3) = STRP%XB(3)   
                  XB_STAIRS(I,4) = STRP%XB_CORE(J,3)
               Else
                  XB_STAIRS(I,3) = STRP%XB_CORE(J,4)   
                  XB_STAIRS(I,4) = STRP%XB(4)
               End If
               STRP%NODE_IOR(2*I) = +1
               If (STR_Length > 0._EB) XB_STAIRS(I,7) = Cos(Atan(STR_Height/STR_Length))
               XB_STAIRS(I,8) = 1._EB
            Else If (XB_LANDINGS(I+1,2)<XB_LANDINGS(I,1)) Then ! From +x to -x
               XB_STAIRS(I,1) = XB_LANDINGS(I+1,2)
               XB_STAIRS(I,2) = XB_LANDINGS(I,1)
               STR_Length = XB_STAIRS(I,2)-XB_STAIRS(I,1)
               If (RIGHT_HANDED) Then
                  XB_STAIRS(I,3) = STRP%XB_CORE(J,4)
                  XB_STAIRS(I,4) = STRP%XB(4)
               Else
                  XB_STAIRS(I,3) = STRP%XB(3)
                  XB_STAIRS(I,4) = STRP%XB_CORE(J,3)
               End If
               STRP%NODE_IOR(2*I) = -1
               If (STR_Length > 0._EB) XB_STAIRS(I,7) = Cos(Atan(STR_Height/STR_Length))
               XB_STAIRS(I,8) = 1._EB
            End If
            If (XB_LANDINGS(I+1,3)>XB_LANDINGS(I,4)) Then ! From -y to +y
               If (RIGHT_HANDED) Then
                  XB_STAIRS(I,1) = STRP%XB_CORE(J,2)
                  XB_STAIRS(I,2) = STRP%XB(2)
               Else
                  XB_STAIRS(I,1) = STRP%XB(1)
                  XB_STAIRS(I,2) = STRP%XB_CORE(J,1)
               End If
               XB_STAIRS(I,3) = XB_LANDINGS(I,4)
               XB_STAIRS(I,4) = XB_LANDINGS(I+1,3)
               STR_Length = XB_STAIRS(I,4)-XB_STAIRS(I,3)
               STRP%NODE_IOR(2*I) = +2
               XB_STAIRS(I,7) = 1._EB
               If (STR_Length > 0._EB) XB_STAIRS(I,8) = Cos(Atan(STR_Height/STR_Length))
            Else If (XB_LANDINGS(I+1,4)<XB_LANDINGS(I,3)) Then ! From +y to -y
               If (RIGHT_HANDED) Then
                  XB_STAIRS(I,1) = STRP%XB(1)
                  XB_STAIRS(I,2) = STRP%XB_CORE(J,1)
               Else
                  XB_STAIRS(I,1) = STRP%XB_CORE(J,2)
                  XB_STAIRS(I,2) = STRP%XB(2)
               End If
               XB_STAIRS(I,3) = XB_LANDINGS(I+1,4)
               XB_STAIRS(I,4) = XB_LANDINGS(I,3)
               STR_Length = XB_STAIRS(I,4)-XB_STAIRS(I,3)
               STRP%NODE_IOR(2*I) = -2
               XB_STAIRS(I,7) = 1._EB
               If (STR_Length > 0._EB) XB_STAIRS(I,8) = Cos(Atan(STR_Height/STR_Length))
            End If
         End Do

         ! Collect sub-node coordinates
         J = 1
         Do I = 1, N_LANDINGS
            STRP%XB_NODE(J,1:6)     = XB_LANDINGS(I,:)
            STRP%NODE_TYPE(J)       = STRS_LANDING_TYPE
            ! Choose core
            Z_TMP = 0.5_EB * (STRP%XB_NODE(J,5)+STRP%XB_NODE(J,6))
            Do K = 1,STRP%N_CORES
               If ((Z_TMP >= STRP%XB_CORE(K,5)) .AND. (Z_TMP <= STRP%XB_CORE(K,6))) Exit
            Enddo
            K = Min(K,STRP%N_CORES)
            STRP%I_CORE(J) = K
            If (I<N_LANDINGS) Then
               STRP%XB_NODE(J+1,1:8) = XB_STAIRS(I,:)
               STRP%NODE_TYPE(J+1)   = STRS_STAIR_TYPE
               ! Choose core
               Z_TMP = 0.5_EB * (STRP%XB_NODE(J+1,5)+STRP%XB_NODE(J+1,6))
               Do K = 1,STRP%N_CORES
                  If ((Z_TMP >= STRP%XB_CORE(K,5)) .AND. (Z_TMP <= STRP%XB_CORE(K,6))) Exit
               Enddo
               K = Min(K,STRP%N_CORES)
               STRP%I_CORE(J+1) = K
            End If
            J = J + 2
         End Do

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

      End Do READ_STRS_LOOP
32    Rewind(LU_INPUT)

    End Subroutine READ_STRS

    Logical Function Is_Within_Bounds(P1x1,P1x2,P1y1,P1y2,P1z1,P1z2,& 
         P2x1,P2x2,P2y1,P2y2,P2z1,P2z2,xtol,ytol,ztol)
      Implicit None
      !
      Real(EB), Intent(IN) :: P1x1,P1x2,P1y1,P1y2,P1z1,P1z2
      Real(EB), Intent(IN) :: P2x1,P2x2,P2y1,P2y2,P2z1,P2z2,xtol,ytol,ztol
      !Returns .TRUE. if P2 is within the bounds of P1 with tolerances.
      Is_Within_Bounds = .False.
      If ( P1x1 >= P2x1-xtol .And. P1x2 <= P2x2+xtol .And. &
           P1y1 >= P2y1-ytol .And. P1y2 <= P2y2+ytol .And. &
           P1z1 >= P2z1-ztol .And. P1z2 <= P2z2+ztol ) Is_Within_Bounds = .True.
      Return
    End Function Is_Within_Bounds

    Subroutine COLLECT_NODE_INFO
      Implicit None
      !
      ! Now exits, doors, corrs and strs are already read in
      If (n_nodes > 0 .And. MYID==Max(0,EVAC_PROCESS)) Then
         n_tmp = 0
         Do n = 1, nmeshes
            If (evacuation_only(n).And.evacuation_grid(n)) Then
               n_tmp = n_tmp + 1
               EVAC_Node_List(n_tmp)%Node_Index = n_tmp
               EVAC_Node_List(n_tmp)%Node_Type  = 'Floor'
               EVAC_Node_List(n_tmp)%ID         = MESH_NAME(n)
               EVAC_Node_List(n_tmp)%GRID_NAME  = MESH_NAME(n)
               EVAC_Node_List(n_tmp)%IMESH      = n
            End If
         End Do
         Do n = 1, n_entrys
            n_tmp = n_tmp + 1
            evac_entrys(n)%INODE             = n_tmp 
            EVAC_Node_List(n_tmp)%Node_Index = n
            EVAC_Node_List(n_tmp)%Node_Type  = 'Entry'
            EVAC_Node_List(n_tmp)%ID         = EVAC_ENTRYS(n)%ID
            EVAC_Node_List(n_tmp)%IMESH      = EVAC_ENTRYS(n)%IMESH
         End Do
         Do n = 1, n_doors
            n_tmp = n_tmp + 1
            evac_doors(n)%INODE              = n_tmp 
            EVAC_Node_List(n_tmp)%Node_Index = n
            EVAC_Node_List(n_tmp)%Node_Type  = 'Door'
            EVAC_Node_List(n_tmp)%ID         = EVAC_DOORS(n)%ID
            EVAC_Node_List(n_tmp)%IMESH      = EVAC_DOORS(n)%IMESH
         End Do
         Do n = 1, n_exits
            n_tmp = n_tmp + 1
            evac_exits(n)%INODE              = n_tmp 
            EVAC_Node_List(n_tmp)%Node_Index = n
            EVAC_Node_List(n_tmp)%Node_Type  = 'Exit'
            EVAC_Node_List(n_tmp)%ID         = EVAC_EXITS(n)%ID
            EVAC_Node_List(n_tmp)%IMESH      = EVAC_EXITS(n)%IMESH
         End Do
         Do n = 1, n_corrs
            n_tmp = n_tmp + 1
            evac_corrs(n)%INODE              = n_tmp 
            EVAC_Node_List(n_tmp)%Node_Index = n
            EVAC_Node_List(n_tmp)%Node_Type  = 'Corr'
            EVAC_Node_List(n_tmp)%ID         = EVAC_CORRS(n)%ID
         End Do
         Do n = 1, n_strs
            n_tmp = n_tmp + 1
            EVAC_STRS(n)%INODE               = n_tmp
            EVAC_Node_List(n_tmp)%Node_Index = n
            EVAC_Node_List(n_tmp)%Node_Type  = 'Stairs'
            EVAC_Node_List(n_tmp)%ID         = EVAC_STRS(n)%ID
         End Do

         ! Check that door/corr/entry/exit have unique names
         Do n = 1, n_nodes - 1
            Do i = n + 1, n_nodes
               If (Trim(EVAC_Node_List(n)%ID) == Trim(EVAC_Node_List(i)%ID)) Then
                  Write(MESSAGE,'(8A)') 'ERROR: ', Trim(EVAC_Node_List(n)%Node_Type), ': ', &
                       Trim(EVAC_Node_List(n)%ID), ' has same ID as ', &
                       Trim(EVAC_Node_List(i)%Node_Type), ': ', Trim(EVAC_Node_List(i)%ID)
                  Call SHUTDOWN(MESSAGE)
               End If
            End Do
         End Do

         ! BUG fix, 17.11.2008
         Do n = 1, n_entrys
            Do i = 1, evac_entrys(n)%N_VENT_FFIELDS
               If (evac_entrys(n)%I_DOOR_NODES(i) < 0) Then
                  evac_entrys(n)%I_DOOR_NODES(i) = EVAC_EXITS(Abs(evac_entrys(n)%I_DOOR_NODES(i)))%INODE
               Else If (evac_entrys(n)%I_DOOR_NODES(i) > 0) Then
                  evac_entrys(n)%I_DOOR_NODES(i) = EVAC_DOORS(Abs(evac_entrys(n)%I_DOOR_NODES(i)))%INODE
               End If
            End Do
         End Do

      End If

    End Subroutine COLLECT_NODE_INFO

    Subroutine READ_ENTRIES
      Implicit None
      !
      ! Read the ENTR lines
      !
      ! Local variables
      Integer nm, i1, i2, j1, j2, NR
      Type (EVAC_ENTR_Type), Pointer :: PNX=>NULL()
      Type (EVAC_PERS_Type), Pointer :: PCP=>NULL()
      Type (EVAC_STRS_Type), Pointer :: STRP=>NULL()
      !
      READ_ENTR_LOOP: Do N = 1, N_ENTRYS
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
         AFTER_REACTION_TIME = .False.
         SHOW          = .True.
         TIME_START          = -Huge(TIME_START)
         TIME_STOP           =  Huge(TIME_STOP)
         MAX_HUMANS    = -1
         MAX_HUMANS_RAMP         = 'null'
         KNOWN_DOOR_NAMES         = 'null'
         KNOWN_DOOR_PROBS         = 1.0_EB
         !
         !
         Call CHECKREAD('ENTR',LU_INPUT,IOS)
         If (IOS == 1) Then
            Exit READ_ENTR_LOOP
         End If
         Read(LU_INPUT,ENTR,End=28,IOSTAT=IOS)
         ! 
         ! Old input used QUANTITY, next lines are needed for that
         If (MYID==Max(0,EVAC_PROCESS) .And. QUANTITY .Ne. 'null') Write (LU_ERR,'(A,A)') &
              ' WARNING: keyword QUANTITY is replaced by AVATAR_COLOR at ENTR line ',Trim(ID)
         If (QUANTITY == 'BLACK')   AVATAR_COLOR = 'BLACK'  
         If (QUANTITY == 'YELLOW')  AVATAR_COLOR = 'YELLOW' 
         If (QUANTITY == 'BLUE')    AVATAR_COLOR = 'BLUE'   
         If (QUANTITY == 'RED')     AVATAR_COLOR = 'RED'    
         If (QUANTITY == 'GREEN')   AVATAR_COLOR = 'GREEN'  
         If (QUANTITY == 'MAGENTA') AVATAR_COLOR = 'MAGENTA'
         If (QUANTITY == 'CYAN')    AVATAR_COLOR = 'CYAN'   
         !
         ! Colors, integer RGB(3), e.g., (23,255,0)
         If (Any(RGB < 0) .And. COLOR=='null') COLOR = 'SKY BLUE'
         If (COLOR /= 'null') Call COLOR2RGB(RGB,COLOR)
         If (Any(AVATAR_RGB < 0) .And. AVATAR_COLOR=='null') AVATAR_COLOR = 'ROYAL BLUE 4'
         If (AVATAR_COLOR /= 'null') Call COLOR2RGB(AVATAR_RGB,AVATAR_COLOR)
         If (COLOR_METHOD == 0 .And. MAX_FLOW > 0.0_EB) Then
            i_avatar_color = i_avatar_color + 1
            EVAC_AVATAR_RGB(1:3,i_avatar_color) = AVATAR_RGB
         End If
         If (MYID /= Max(0,EVAC_PROCESS)) Cycle READ_ENTR_LOOP

         If (MAX_HUMANS < 0) MAX_HUMANS = Huge(MAX_HUMANS)
         If (MAX_FLOW <= 0.0_EB) MAX_HUMANS = 0
         IF (MAX_HUMANS_RAMP/='null') THEN
            CALL GET_RAMP_INDEX(MAX_HUMANS_RAMP,'TIME',NR)
            MAX_HUMANS = -NR
         ENDIF 

         PNX=>EVAC_ENTRYS(N)

         PNX%RGB = RGB
         PNX%AVATAR_RGB =AVATAR_RGB
         If (COLOR_METHOD == 0 .And. MAX_FLOW > 0.0_EB) PNX%Avatar_Color_Index = i_avatar_color

         If (EVAC_MESH /= 'null') Then
            MESH_ID = EVAC_MESH
            If (MYID==Max(0,EVAC_PROCESS)) Write (LU_ERR,'(A,A)') &
                 ' WARNING: keyword EVAC_MESH is replaced by MESH_ID at ENTR line ', Trim(ID)
         End If

         If (Trim(KNOWN_DOOR_NAMES(51)) /= 'null') Then
            Write(MESSAGE,'(A,A,A)') 'ERROR: ENTR line ',Trim(ID), ' problem with KNOWN_DOOR_NAMES'
            Call SHUTDOWN(MESSAGE)
         End If
         If (Trim(KNOWN_DOOR_NAMES(1)) == 'null') Then
            i = 0 ! no doors given
         Else
            i = 50 ! known door names given
            Do While ( Trim(KNOWN_DOOR_NAMES(i)) == 'null' .And. i > 0)
               i = i-1
            End Do
         End If
         PNX%N_VENT_FFIELDS = i
         Allocate(PNX%I_DOOR_NODES(0:i),STAT=IZERO)
         Call ChkMemErr('Read_Evac','PNX%I_DOOR_NODES',IZERO) 
         Allocate(PNX%I_VENT_FFIELDS(0:i),STAT=IZERO)
         Call ChkMemErr('Read_Evac','PNX%I_VENT_FFIELDS',IZERO) 
         Allocate(PNX%P_VENT_FFIELDS(0:i),STAT=IZERO)
         Call ChkMemErr('Read_Evac','PNX%P_VENT_FFIELDS',IZERO) 
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

         Do I=1,5,2
            If (XB(I) > XB(I+1)) Then
               DUMMY   = XB(I)
               XB(I)   = XB(I+1)
               XB(I+1) = DUMMY
            End If
         End Do
         ! 
         ! Check which evacuation floor
         ii = 0
         n_tmp = 0
         PNX_MeshLoop: Do i = 1, nmeshes
            If (evacuation_only(i) .And. evacuation_grid(i)) Then
               n_tmp = n_tmp + 1
               If (Is_Within_Bounds(XB(1),XB(2),XB(3),XB(4),XB(5),XB(6),&
                  Meshes(i)%XS,Meshes(i)%XF,Meshes(i)%YS,Meshes(i)%YF,Meshes(i)%ZS,Meshes(i)%ZF, 0._EB, 0._EB, 0._EB)) Then
                  If (Trim(MESH_ID) == 'null' .Or. Trim(MESH_ID) == Trim(MESH_NAME(i))) Then
                     ii = ii + 1
                     PNX%IMESH = i
                     PNX%TO_INODE = n_tmp
                     PNX%TO_NODE  = MESH_NAME(i)
                  End If
               End If
            End If
         End Do PNX_MeshLoop
         If (PNX%IMESH == 0) Then
            Write(MESSAGE,'(A,A,A)') 'ERROR: ENTR line ',Trim(ID), ' problem with IMESH, no mesh found'
            Call SHUTDOWN(MESSAGE)
         End If
         If (ii > 1) Then
            Write(MESSAGE,'(A,A,A)') 'ERROR: ENTR line ',Trim(ID), ' not an unique mesh found '
            Call SHUTDOWN(MESSAGE)
         End If

         nm = PNX%IMESH
 
         If (XB(1)/=XB(2) .And. XB(3)/=XB(4)) Then
            Write(MESSAGE,'(A,I4,A)') 'ERROR: ENTR',N,' must be a plane'
            Call SHUTDOWN(MESSAGE)
         Endif

         ! User input
         PNX%X1 = XB(1)
         PNX%X2 = XB(2)
         PNX%Y1 = XB(3)
         PNX%Y2 = XB(4)

         ! Move user input to mesh cell boundaries
         XB(1) = Max(XB(1),Meshes(nm)%XS)
         XB(2) = Min(XB(2),Meshes(nm)%XF)
         XB(3) = Max(XB(3),Meshes(nm)%YS)
         XB(4) = Min(XB(4),Meshes(nm)%YF)
         XB(5) = Max(XB(5),Meshes(nm)%ZS)
         XB(6) = Min(XB(6),Meshes(nm)%ZF)

         I1 = Nint( GINV(XB(1)-Meshes(nm)%XS,1,nm)*Meshes(nm)%RDXI ) 
         I2 = Nint( GINV(XB(2)-Meshes(nm)%XS,1,nm)*Meshes(nm)%RDXI )
         J1 = Nint( GINV(XB(3)-Meshes(nm)%YS,2,nm)*Meshes(nm)%RDETA) 
         J2 = Nint( GINV(XB(4)-Meshes(nm)%YS,2,nm)*Meshes(nm)%RDETA)

         XB(1) = Meshes(nm)%X(I1)
         XB(2) = Meshes(nm)%X(I2)
         XB(3) = Meshes(nm)%Y(J1)
         XB(4) = Meshes(nm)%Y(J2)
         If ( Abs(XB(1)-PNX%X1)>1.E-4_EB .Or. Abs(XB(2)-PNX%X2)>1.E-4_EB .Or. &
              Abs(XB(3)-PNX%Y1)>1.E-4_EB .Or. Abs(XB(4)-PNX%Y2)>1.E-4_EB ) Then
            Write(LU_ERR,fmt='(a,a,a,a)') ' WARNING: Entr line ',Trim(ID),' XB adjusted to mesh ',Trim(MESH_NAME(nm))
            Write(LU_ERR,fmt='(a,6f12.4)') 'Old XB:', PNX%X1,PNX%X2,PNX%Y1,PNX%Y2,PNX%Z1,PNX%Z2
            Write(LU_ERR,fmt='(a,6f12.4)') 'New XB:', XB(1:6)
         End If

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
         Do ipc= 1, npc_pers
            pcp => evac_person_classes(ipc)
            If ( pcp%id == PERS_ID ) PNX%IPC = IPC
         End Do

         Select Case (IOR)
         Case (-1,+1)
            If (WIDTH <= 0.0_EB) Then
               PNX%Width = XB(4) - XB(3)
            Else
               PNX%Width = WIDTH
            End If
            PNX%ORIENTATION(1)=-Real(Sign(1,IOR),EB)
         Case (-2,+2)
            If (WIDTH <= 0.0_EB) Then
               PNX%Width = XB(2) - XB(1)
            Else
               PNX%Width = WIDTH
            End If
            PNX%ORIENTATION(2)=-Real(Sign(1,IOR),EB)
         Case (3)
            If ( (XB(4)-XB(3)) <= 0.0_EB .Or. (XB(2)-XB(1)) <= 0.0_EB) Then
               Write(MESSAGE,'(A,I4,A)') 'ERROR: ENTR',N,' IOR=3 but not 3-dim object'
               Call SHUTDOWN(MESSAGE)
            End If
            PNX%ORIENTATION(3)=-Real(Sign(1,IOR),EB)
         Case (0)
            If ( (XB(4)-XB(3)) <= 0.0_EB .Or. (XB(2)-XB(1)) <= 0.0_EB) Then
               Write(MESSAGE,'(A,I4,A)') 'ERROR: ENTR',N,' no IOR but not 3-dim object'
               Call SHUTDOWN(MESSAGE)
            End If
         Case Default
            Write(MESSAGE,'(A,I4,A)') 'ERROR: ENTR',N,' problem with IOR'
            Call SHUTDOWN(MESSAGE)
         End Select

         ! Check if entry leads to Stairs
         PNX%STR_INDX = 0
         PNX%STR_SUB_INDX = 0
         CheckEntrStrLoop: Do i = 1, N_STRS
            If (Trim(EVAC_STRS(i)%MESH_ID)==PNX%TO_NODE) Then     
               STRP=>EVAC_STRS(i)
               PNX%STR_INDX = i
               Do j = 1,STRP%N_NODES
                  If ( PNX%Z1 >= STRP%XB_NODE(j,5) .And. PNX%Z2 <= STRP%XB_NODE(j,6) .And. &
                       PNX%Y1 >= STRP%XB_NODE(j,3) .And. PNX%Y2 <= STRP%XB_NODE(j,4) .And. &
                       PNX%X1 >= STRP%XB_NODE(j,1) .And. PNX%X2 <= STRP%XB_NODE(j,2) ) Then
                     PNX%STR_SUB_INDX = j
                     Exit CheckEntrStrLoop
                  End If
               End Do
            End If
         End Do CheckEntrStrLoop

         ! Use the main_evac_grid flow field if none is given
         If (Trim(FLOW_FIELD_ID) == 'null') Then
            PNX%GRID_NAME  = Trim(PNX%TO_NODE)
         Else
            PNX%GRID_NAME  = Trim(FLOW_FIELD_ID)
         End If

         Do i = 1, PNX%N_VENT_FFIELDS
            ! P = 0 or 1 for entrys.
            If ( .Not.( Abs(KNOWN_DOOR_PROBS(i)-1.0_EB) < 0.0001_EB .Or. Abs(KNOWN_DOOR_PROBS(i)) < 0.0001_EB ) )  Then
               Write(MESSAGE,'(A,A,A,f12.6,A)') 'ERROR: ENTR line ',Trim(PNX%ID), &
                    ' problem with probability, ', KNOWN_DOOR_PROBS(i),' it should be zero or one.'
               Call SHUTDOWN(MESSAGE)
            End If
            PNX%P_VENT_FFIELDS(i) = Max(0.0_EB,KNOWN_DOOR_PROBS(i))
            PNX%I_VENT_FFIELDS(i) = 0
            PNX%I_DOOR_NODES(i) = 0
            Do j = 1, n_exits
               If ( Trim(EVAC_EXITS(j)%ID) == Trim(KNOWN_DOOR_NAMES(i)) ) Then
                  PNX%I_VENT_FFIELDS(i) = EVAC_EXITS(j)%I_VENT_FFIELD
                  ! BUG fix, 17.11.2008
                  ! PNX%I_DOOR_NODES(i)   = EVAC_EXITS(j)%INODE
                  PNX%I_DOOR_NODES(i)   = -j
               End If
            End Do
            Do j = 1, n_doors
               If ( Trim(EVAC_DOORS(j)%ID) == Trim(KNOWN_DOOR_NAMES(i)) ) Then
                  PNX%I_VENT_FFIELDS(i) = EVAC_DOORS(j)%I_VENT_FFIELD
                  ! BUG fix, 17.11.2008
                  ! PNX%I_DOOR_NODES(i)   = EVAC_DOORS(j)%INODE
                  PNX%I_DOOR_NODES(i)   = +j
               End If
            End Do
            If ( PNX%I_VENT_FFIELDS(i)*PNX%I_DOOR_NODES(i) == 0 ) Then
               Write(MESSAGE,'(A,A,A,A,A)') 'ERROR: ENTR line ',Trim(PNX%ID), &
                    ' problem with door/exit names, ', Trim(KNOWN_DOOR_NAMES(i)),' not found'
               Call SHUTDOWN(MESSAGE)
            End If
         End Do
         !
         ! No known doors given, use the flow_field_id value
         ! 
         PNX%P_VENT_FFIELDS(0) = 1.0_EB
         PNX%I_VENT_FFIELDS(0) = 0
         PNX%I_DOOR_NODES(0) = 0
         PNX_Mesh2Loop: Do i = 1, nmeshes
            If ( evacuation_only(i) .And. Trim(PNX%GRID_NAME) == Trim(MESH_NAME(i)) ) Then
               PNX%I_VENT_FFIELDS(0) = i
               Exit PNX_Mesh2Loop
            End If
         End Do PNX_Mesh2Loop
         If ( PNX%I_VENT_FFIELDS(0) == 0 ) Then
            Write(MESSAGE,'(A,A,A,A,A)') 'ERROR: ENTR line ',Trim(PNX%ID),&
                 ' problem with flow field name, ', Trim(PNX%GRID_NAME),' not found'
            Call SHUTDOWN(MESSAGE)
         End If
         ! 
      End Do READ_ENTR_LOOP
28    Rewind(LU_INPUT)

    End Subroutine READ_ENTRIES

    Subroutine READ_EVAC_LINES
      Implicit None
      !
      ! Read the EVAC lines
      ! 
      ! Local variables
      Type (EVACUATION_Type), Pointer :: HPT=>NULL()
      Type (EVAC_PERS_Type),  Pointer :: PCP=>NULL()

      READ_EVAC_LOOP: Do N=1,NPC_EVAC
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
         EVACFILE                 = .False.
         AFTER_REACTION_TIME      = .False.
         SHOW                     = .True.
         TIME_START               = -99.0_EB
         TIME_STOP                = -99.0_EB
         GN_MIN                   = 1
         GN_MAX                   = 1      

         KNOWN_DOOR_NAMES         = 'null'
         KNOWN_DOOR_PROBS         = 1.0_EB
         !
         Call CHECKREAD('EVAC',LU_INPUT,IOS)
         If (IOS == 1) Then
            Exit READ_EVAC_LOOP
         End If
         Read(LU_INPUT,EVAC,End=25,IOSTAT=IOS)
         ! 
         ! Old input used QUANTITY, next lines are needed for that
         If (MYID==Max(0,EVAC_PROCESS) .And. QUANTITY .Ne. 'null') Write (LU_ERR,'(A,A)') &
              ' WARNING: keyword QUANTITY is replaced by AVATAR_COLOR at EVAC line ',Trim(ID)
         If (QUANTITY == 'BLACK')   AVATAR_COLOR = 'BLACK'  
         If (QUANTITY == 'YELLOW')  AVATAR_COLOR = 'YELLOW' 
         If (QUANTITY == 'BLUE')    AVATAR_COLOR = 'BLUE'   
         If (QUANTITY == 'RED')     AVATAR_COLOR = 'RED'    
         If (QUANTITY == 'GREEN')   AVATAR_COLOR = 'GREEN'  
         If (QUANTITY == 'MAGENTA') AVATAR_COLOR = 'MAGENTA'
         If (QUANTITY == 'CYAN')    AVATAR_COLOR = 'CYAN'   
         !
         ! Colors, integer RGB(3), e.g., (23,255,0)
         If (Any(RGB < 0) .And. COLOR=='null') COLOR = 'BLACK'
         If (COLOR /= 'null') Call COLOR2RGB(RGB,COLOR)
         If (Any(AVATAR_RGB < 0) .And. AVATAR_COLOR=='null') AVATAR_COLOR = 'ROYAL BLUE 4'
         If (AVATAR_COLOR /= 'null') Call COLOR2RGB(AVATAR_RGB,AVATAR_COLOR)
         If (COLOR_METHOD == 0 .And. NUMBER_INITIAL_PERSONS > 0) Then
            i_avatar_color = i_avatar_color + 1
            EVAC_AVATAR_RGB(1:3,i_avatar_color) = AVATAR_RGB
         End If

         If (MYID /= Max(0,EVAC_PROCESS)) Cycle READ_EVAC_LOOP

         HPT=>EVACUATION(N)

         HPT%RGB = RGB
         HPT%AVATAR_RGB = AVATAR_RGB
         If (COLOR_METHOD == 0 .And. NUMBER_INITIAL_PERSONS > 0) HPT%Avatar_Color_Index = i_avatar_color
         If (EVAC_MESH /= 'null') Then
            MESH_ID = EVAC_MESH
            If (MYID==Max(0,EVAC_PROCESS)) Write (LU_ERR,'(A,A)') &
                 ' WARNING: keyword EVAC_MESH is replaced by MESH_ID at EVAC line ', Trim(ID)
         End If
         If (Trim(PERS_ID) == 'null') Then
            Write(MESSAGE,'(A,A,A)') 'ERROR: EVAC line ',Trim(ID),' no PERS_ID given'
            Call SHUTDOWN(MESSAGE)
         Else
            ii = 1
            Do i = 1,NPC_PERS
               If (Trim(EVAC_PERSON_CLASSES(i)%ID) == Trim(PERS_ID)) Cycle
               ii = ii + 1
            End Do
            If (ii > NPC_PERS) Then
               Write(MESSAGE,'(A,A,A)') 'ERROR: EVAC line ',Trim(ID), ' prblem with PERS_ID'
               Call SHUTDOWN(MESSAGE)
            End If
         End If

         If (Trim(KNOWN_DOOR_NAMES(51)) /= 'null') Then
            Write(MESSAGE,'(A,A,A)') 'ERROR: EVAC line ',Trim(ID), ' problem with KNOWN_DOOR_NAMES'
            Call SHUTDOWN(MESSAGE)
         End If
         If (Trim(KNOWN_DOOR_NAMES(1)) == 'null') Then
            i = 0 ! no doors given
         Else
            i = 50 ! known door names given
            Do While ( Trim(KNOWN_DOOR_NAMES(i)) == 'null' .And. i > 1)
               i = i-1
            End Do
         End If
         HPT%N_VENT_FFIELDS = i
         Allocate(HPT%I_DOOR_NODES(0:i),STAT=IZERO)
         Call ChkMemErr('Read_Evac','HPT%I_DOOR_NODES',IZERO) 
         Allocate(HPT%I_VENT_FFIELDS(0:i),STAT=IZERO)
         Call ChkMemErr('Read_Evac','HPT%I_VENT_FFIELDS',IZERO) 
         Allocate(HPT%P_VENT_FFIELDS(0:i),STAT=IZERO)
         Call ChkMemErr('Read_Evac','HPT%P_VENT_FFIELDS',IZERO) 
         ! 
         If (NUMBER_INITIAL_PERSONS > 0 .And. RESTART) NUMBER_INITIAL_PERSONS =  0
         ! 
         If (NUMBER_INITIAL_PERSONS > 0) EVACFILE = .True.
         !
         HPT%CLASS_NAME = PERS_ID
         HPT%T_START    = TIME_START

         HPT%GN_MIN = GN_MIN
         HPT%GN_MAX = GN_MAX

         ! input in degrees, internal units are radians, [0,2pi)
         ! If no angle is given then use random angle
         If (ANGLE > -999.9_EB) Then
            ANGLE = Pi*ANGLE/(180.0_EB)
            Do While (ANGLE >= 2.0_EB*Pi)
               ANGLE = ANGLE - 2.0_EB*Pi
            End Do
            Do While (ANGLE < 0.0_EB)
               ANGLE = ANGLE + 2.0_EB*Pi
            End Do
         End If
         HPT%Angle = ANGLE

         HPT%IPC = 0
         Do ipc= 1, npc_pers
            pcp => evac_person_classes(ipc)
            If ( Trim(pcp%id) == Trim(PERS_ID) ) HPT%IPC = IPC
         End Do
         ! 
         HPT%SAMPLING = SAMPLING_FACTOR
         !
         If (NUMBER_INITIAL_PERSONS > 0) Then
            Do I=1,5,2
               If (XB(I) > XB(I+1)) Then
                  DUMMY   = XB(I)
                  XB(I)   = XB(I+1)
                  XB(I+1) = DUMMY
               End If
            End Do
         End If
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
         HPT%ID = Trim(ID)
         HPT%SHOW   = SHOW

         ! Check which evacuation floor
         ii = 0
         HP_MeshLoop: Do i = 1, nmeshes
            If (evacuation_only(i) .And. evacuation_grid(i)) Then
               If ( Is_Within_Bounds(HPT%X1,HPT%X2,HPT%Y1,HPT%Y2,HPT%Z1,HPT%Z2,&
                    Meshes(i)%XS,Meshes(i)%XF,Meshes(i)%YS,Meshes(i)%YF,Meshes(i)%ZS,Meshes(i)%ZF, 0._EB, 0._EB, 0._EB)) Then
                  If (Trim(MESH_ID) == 'null' .Or. Trim(MESH_ID) == Trim(MESH_NAME(i))) Then
                     ii = ii + 1
                     HPT%IMESH = i
                  End If
               End If
            End If
         End Do HP_MeshLoop
         If (HPT%IMESH == 0) Then
            Write(MESSAGE,'(A,A,A)') 'ERROR: EVAC line ',Trim(ID),' problem with IMESH, no mesh found'
            Call SHUTDOWN(MESSAGE)
         End If
         If (ii > 1) Then
            Write(MESSAGE,'(A,A,A)') 'ERROR: EVAC line ',Trim(ID), ' not an unique mesh found '
            Call SHUTDOWN(MESSAGE)
         End If

         ! Use the main_evac_grid flow field if none is given
         If (Trim(FLOW_FIELD_ID) == 'null') Then
            HPT%GRID_NAME  = Trim(MESH_NAME(HPT%IMESH))
         Else
            HPT%GRID_NAME  = FLOW_FIELD_ID
         End If
         !
         Do i = 1, HPT%N_VENT_FFIELDS
            HPT%P_VENT_FFIELDS(i) = Max(0.0_EB,KNOWN_DOOR_PROBS(i))
            HPT%I_VENT_FFIELDS(i) = 0
            HPT%I_DOOR_NODES(i) = 0
            Do j = 1, n_exits
               If ( Trim(EVAC_EXITS(j)%ID) == Trim(KNOWN_DOOR_NAMES(i)) ) Then
                  HPT%I_VENT_FFIELDS(i) = EVAC_EXITS(j)%I_VENT_FFIELD
                  HPT%I_DOOR_NODES(i)   = EVAC_EXITS(j)%INODE
               End If
            End Do
            Do j = 1, n_doors
               If ( Trim(EVAC_DOORS(j)%ID) == Trim(KNOWN_DOOR_NAMES(i)) ) Then
                  HPT%I_VENT_FFIELDS(i) = EVAC_DOORS(j)%I_VENT_FFIELD
                  HPT%I_DOOR_NODES(i)   = EVAC_DOORS(j)%INODE
               End If
            End Do
            If ( HPT%I_VENT_FFIELDS(i)*HPT%I_DOOR_NODES(i) == 0 ) Then
               Write(MESSAGE,'(A,A,A,A,A)') 'ERROR: EVAC line ',Trim(HPT%ID), &
                    ' problem with door/exit names, ', Trim(KNOWN_DOOR_NAMES(i)),' not found'
               Call SHUTDOWN(MESSAGE)
            End If
         End Do
         !
         ! No known doors given, use the flow_field_id value
         ! 
         HPT%P_VENT_FFIELDS(0) = 1.0_EB
         HPT%I_VENT_FFIELDS(0) = 0
         HPT%I_DOOR_NODES(0) = 0
         HP_Mesh2Loop: Do i = 1, nmeshes
            If ( evacuation_only(i) .And. Trim(HPT%GRID_NAME) == Trim(MESH_NAME(i)) ) Then
               HPT%I_VENT_FFIELDS(0) = i
               Exit HP_Mesh2Loop
            End If
         End Do HP_Mesh2Loop
         If ( HPT%I_VENT_FFIELDS(0) == 0 ) Then
            Write(MESSAGE,'(A,A,A,A,A)') 'ERROR: EVAC line ',Trim(HPT%ID), &
                 ' problem with flow field name, ', Trim(HPT%GRID_NAME),' not found'
            Call SHUTDOWN(MESSAGE)
         End If
         !
      End Do READ_EVAC_LOOP
25    Rewind(LU_INPUT)

      N_EVAC = 1
      Allocate(EVAC_CLASS_NAME(N_EVAC),STAT=IZERO)
      Call ChkMemErr('READ_EVAC','EVAC_CLASS_NAME',IZERO)
      Allocate(EVAC_CLASS_RGB(3,N_EVAC),STAT=IZERO)
      Call ChkMemErr('READ_EVAC','EVAC_CLASS_RGB',IZERO)
      EVAC_CLASS_NAME(1) = 'Human'
      Do N = 1, N_EVAC
         EVAC_CLASS_RGB(1:3,N) = (/ 39, 64,139/)  ! ROYAL BLUE 4
      End Do
      ! Default color table for agents

    End Subroutine READ_EVAC_LINES

    Subroutine READ_EVHO
      Implicit None
      !
      ! Read the EVHO lines
      !
      ! Local variables
      Type (EVAC_HOLE_Type),  Pointer :: EHX=>NULL()

      READ_EVHO_LOOP: Do N = 1, N_HOLES
         If (MYID /= Max(0,EVAC_PROCESS)) Cycle READ_EVHO_LOOP
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
         SHOW          = .True.
         !
         Call CHECKREAD('EVHO',LU_INPUT,IOS)
         If (IOS == 1) Then
            Exit READ_EVHO_LOOP
         End If
         Read(LU_INPUT,EVHO,End=30,IOSTAT=IOS)
         !
         Do I=1,5,2
            If (XB(I) > XB(I+1)) Then
               DUMMY   = XB(I)
               XB(I)   = XB(I+1)
               XB(I+1) = DUMMY
            End If
         End Do
         If (EVAC_MESH /= 'null') Then
            MESH_ID = EVAC_MESH
            If (MYID==Max(0,EVAC_PROCESS)) Write (LU_ERR,'(A,A)') &
                 ' WARNING: keyword EVAC_MESH is replaced by MESH_ID at EVHO line ', Trim(ID)
         End If

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
         EHX_MeshLoop: Do i = 1, nmeshes
            If (evacuation_only(i) .And. evacuation_grid(i)) Then
               If ( (EHX%Z1 >= Meshes(i)%ZS .And. EHX%Z2 <= Meshes(i)%ZF).And. &
                    (EHX%Y1 >= Meshes(i)%YS .And. EHX%Y2 <= Meshes(i)%YF).And. &
                    (EHX%X1 >= Meshes(i)%XS .And. EHX%X2 <= Meshes(i)%XF)) Then
                  If (Trim(MESH_ID) == 'null' .Or. Trim(MESH_ID) == Trim(MESH_NAME(i))) Then
                     ii = ii + 1
                     EHX%IMESH = i
                     EHX%GRID_NAME  = MESH_NAME(i)
                  End If
               End If
            End If
         End Do EHX_MeshLoop
         If (EHX%IMESH == 0) Then
            Write(MESSAGE,'(A,A,A)') 'ERROR: EVHO line ',Trim(EHX%ID), ' problem with IMESH, no mesh found'
            Call SHUTDOWN(MESSAGE)
         End If
         If (ii > 1) Then
            Write(MESSAGE,'(A,A,A)') 'ERROR: EVHO line ',Trim(EHX%ID), ' not an unique mesh found '
            Call SHUTDOWN(MESSAGE)
         End If

         ! Colors, integer RGB(3), e.g., (23,255,0)
         If (Any(RGB < 0) .And. COLOR=='null') COLOR = 'BLACK'
         If (COLOR /= 'null') Call COLOR2RGB(RGB,COLOR)
         EHX%RGB = RGB
      End Do READ_EVHO_LOOP
30    Rewind(LU_INPUT)

    End Subroutine READ_EVHO

    Subroutine READ_EVSS
      Implicit None
      !
      ! Read the EVSS lines
      !
      ! Local variables
      Type (EVAC_SSTAND_Type), Pointer :: ESS=>NULL()
      Real(EB) :: X, Y, Z

      READ_EVSS_LOOP: Do N = 1, N_SSTANDS
         If (MYID /= Max(0,EVAC_PROCESS)) Cycle READ_EVSS_LOOP
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
         USE_V0        = .False.
         SHOW          = .True.
         !
         Call CHECKREAD('EVSS',LU_INPUT,IOS)
         If (IOS == 1) Then
            Exit READ_EVSS_LOOP
         End If
         Read(LU_INPUT,EVSS,End=31,IOSTAT=IOS)
         !
         Do I=1,5,2
            If (XB(I) > XB(I+1)) Then
               DUMMY   = XB(I)
               XB(I)   = XB(I+1)
               XB(I+1) = DUMMY
            End If
         End Do
         If (EVAC_MESH /= 'null') Then
            MESH_ID = EVAC_MESH
            If (MYID==Max(0,EVAC_PROCESS)) Write (LU_ERR,'(A,A)') &
                 ' WARNING: keyword EVAC_MESH is replaced by MESH_ID at EVSS line ', Trim(ID)
         End If
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
         If ( (ESS%H - ESS%H0) < 0.0_EB) Then
            ESS%FAC_V0_UP   = FAC_V0_DOWN
            ESS%FAC_V0_DOWN = FAC_V0_UP
            ESS%Esc_SpeedDn   = Max(0.0_EB,+ESC_SPEED)
            ESS%Esc_SpeedUp   = Max(0.0_EB,-ESC_SPEED)
         Else
            ESS%FAC_V0_UP   = FAC_V0_UP
            ESS%FAC_V0_DOWN = FAC_V0_DOWN
            ESS%Esc_SpeedUp   = Max(0.0_EB,+ESC_SPEED)
            ESS%Esc_SpeedDn   = Max(0.0_EB,-ESC_SPEED)
         End If
         ESS%IOR    = IOR
         ESS%UBAR0  = UBAR0
         ESS%VBAR0  = VBAR0
         ESS%Use_v0 = USE_V0
         ESS%SHOW   = SHOW
         ! 
         ! Check which evacuation floor
         ii = 0
         ESS_MeshLoop: Do i = 1, nmeshes
            If (evacuation_only(i) .And. evacuation_grid(i)) Then
               If ( (ESS%Z1 >= Meshes(i)%ZS .And. ESS%Z2 <= Meshes(i)%ZF).And. &
                    (ESS%Y1 >= Meshes(i)%YS .And. ESS%Y2 <= Meshes(i)%YF).And. &
                    (ESS%X1 >= Meshes(i)%XS .And. ESS%X2 <= Meshes(i)%XF)) Then
                  If (Trim(MESH_ID) == 'null' .Or. Trim(MESH_ID) == Trim(MESH_NAME(i))) Then
                     ii = ii + 1
                     ESS%IMESH = i
                     ESS%GRID_NAME  = MESH_NAME(i)
                  End If
               End If
            End If
         End Do ESS_MeshLoop
         If (ESS%IMESH == 0) Then
            Write(MESSAGE,'(A,A,A)') 'ERROR: EVSS line ',Trim(ESS%ID), ' problem with IMESH, no mesh found'
            Call SHUTDOWN(MESSAGE)
         End If
         If (ii > 1) Then
            Write(MESSAGE,'(A,A,A)') 'ERROR: EVSS line ',Trim(ESS%ID), ' not an unique mesh found '
            Call SHUTDOWN(MESSAGE)
         End If

         If (Abs(ESS%H0-ESS%H) < 1.0E-3) Then
            ESS%H = ESS%H0
            IOR=1
            ESS%IOR    = IOR
         End If

         Select Case (IOR)
         Case(-1,+1)
            ESS%S = Sqrt((ESS%X2-ESS%X1)**2 + (ESS%H-ESS%H0)**2)
            ESS%COS_X = Abs(ESS%X2-ESS%X1)/ Sqrt((ESS%X2-ESS%X1)**2 + (ESS%H-ESS%H0)**2)
            ESS%COS_Y = 1.0_EB
            ESS%SIN_X = Abs(ESS%H-ESS%H0)/ Sqrt((ESS%X2-ESS%X1)**2 + (ESS%H-ESS%H0)**2)
            ESS%SIN_Y = 0.0_EB
            ESS%ORIENTATION(1) = IOR*(ESS%H-ESS%H0)/ Sqrt((ESS%X2-ESS%X1)**2 + (ESS%H-ESS%H0)**2)
            ESS%ORIENTATION(2) = 0.0_EB
            ESS%ORIENTATION(3) = ESS%COS_X
         Case(-2,+2)
            ESS%S = Sqrt((ESS%Y2-ESS%Y1)**2 + (ESS%H-ESS%H0)**2)
            ESS%COS_X = 1.0_EB
            ESS%COS_Y = Abs(ESS%Y2-ESS%Y1)/ Sqrt((ESS%Y2-ESS%Y1)**2 + (ESS%H-ESS%H0)**2)
            ESS%SIN_X = 0.0_EB
            ESS%SIN_Y = Abs(ESS%H-ESS%H0)/ Sqrt((ESS%Y2-ESS%Y1)**2 + (ESS%H-ESS%H0)**2)
            ESS%ORIENTATION(1) = 0.0_EB
            ESS%ORIENTATION(2) = 0.5_EB*IOR*(ESS%H-ESS%H0)/ Sqrt((ESS%Y2-ESS%Y1)**2 + (ESS%H-ESS%H0)**2)
            ESS%ORIENTATION(3) = ESS%COS_Y
         Case Default
            Write(MESSAGE,'(A,I4,A)') 'ERROR: EVSS',N,' problem with IOR'
            Call SHUTDOWN(MESSAGE)
         End Select

         ! Colors, integer RGB(3), e.g., (23,255,0)
         ! If (Any(RGB < 0) .And. COLOR=='null') COLOR = 'BLACK'
         If (Any(RGB < 0) .And. COLOR=='null') COLOR = 'AZURE 2'  ! RGB = 193 205 205
         If (COLOR /= 'null') Call COLOR2RGB(RGB,COLOR)
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

      End Do READ_EVSS_LOOP
31    Rewind(LU_INPUT)

    End Subroutine READ_EVSS

    Subroutine CHECK_EVAC_NODES
      Implicit None
      !
      ! Local variables
      Type (EVAC_ENTR_Type), Pointer :: PNX=>NULL()
      Type (EVAC_DOOR_Type), Pointer :: PDX=>NULL()
      Type (EVAC_STRS_Type), Pointer :: STRP=>NULL()
      !
      ! Set the IMESH and IMESH2 for corridors
      Do n = 1, n_corrs
         Nodeloop2: Do i = 1, n_nodes
            If (EVAC_Node_List(i)%ID == EVAC_CORRS(n)%TO_NODE) Then
               EVAC_CORRS(n)%INODE2 = i
               If ( Trim(EVAC_Node_List(i)%Node_Type) == 'Exit' .Or. &
                    Trim(EVAC_Node_List(i)%Node_Type) == 'Door' .Or. &
                    Trim(EVAC_Node_List(i)%Node_Type) == 'Entry' ) Then
                  EVAC_CORRS(n)%IMESH  = EVAC_Node_List(i)%IMESH
                  EVAC_CORRS(n)%IMESH2 = EVAC_Node_List(i)%IMESH
               Else
                  EVAC_CORRS(n)%IMESH  = 0
                  EVAC_CORRS(n)%IMESH2 = n_egrids
               End If
               EVAC_Node_List(evac_corrs(n)%INODE)%IMESH = EVAC_CORRS(n)%IMESH2
               Exit Nodeloop2
            End If
         End Do Nodeloop2
         If (EVAC_CORRS(n)%INODE2 == 0 .Or. EVAC_CORRS(n)%IMESH2 == 0) Then
            Write(MESSAGE,'(A,I4,A)') 'ERROR: CORR',n,' problem with TO_NODE, loop2'
            Call SHUTDOWN(MESSAGE)
         End If
      End Do

      !
      Do n = 1, n_doors
         NodeLoop: Do i = 1, n_nodes
            If (EVAC_Node_List(i)%ID == EVAC_DOORS(n)%TO_NODE) Then
               EVAC_DOORS(n)%INODE2 = i
               EVAC_DOORS(n)%IMESH2 = EVAC_Node_List(i)%IMESH
               Exit NodeLoop
            End If
         End Do NodeLoop
         If (EVAC_DOORS(n)%INODE2 == 0 .Or. EVAC_DOORS(n)%IMESH2 == 0) Then
            Write(MESSAGE,'(A,I4,A)') 'ERROR: DOOR',n,' problem with TO_NODE'
            Call SHUTDOWN(MESSAGE)
         End If
      End Do


      Do n = 1, n_doors
         If (EVAC_DOORS(n)%KEEP_XY) Then
            i = EVAC_DOORS(n)%INODE2
            If (Trim(EVAC_Node_List(i)%Node_Type) == 'Door') Then
               PDX => EVAC_DOORS(EVAC_Node_List(i)%Node_Index)
               If ((EVAC_DOORS(n)%IOR /= -PDX%IOR) .Or. Abs(EVAC_DOORS(n)%Width-PDX%Width) > 0.1_EB ) Then
                  Write(MESSAGE,'(A,I4,A)') 'ERROR: DOOR',N,' KEEP_XY Problem'
                  Call SHUTDOWN(MESSAGE)
               End If
            End If
            If (Trim(EVAC_Node_List(i)%Node_Type) == 'Entry') Then
               PNX => EVAC_ENTRYS(EVAC_Node_List(i)%Node_Index)
               If ((EVAC_DOORS(n)%IOR /= PNX%IOR) .Or. Abs(EVAC_DOORS(n)%Width-PNX%Width) > 0.1_EB ) Then
                  Write(MESSAGE,'(A,I4,A)') 'ERROR: DOOR',N,' KEEP_XY Problem'
                  Call SHUTDOWN(MESSAGE)
               End If
            End If
         End If
      End Do

   Do N = 1,N_DOORS
      PDX => EVAC_DOORS(N)
      ! Check if door leads to Stairs
       PDX%STR_INDX = 0
       PDX%STR_SUB_INDX = 0
       CheckDoorStrLoop: Do i = 1, N_STRS
         STRP=>EVAC_STRS(i)
         If (STRP%IMESH==PDX%IMESH .OR. STRP%IMESH==PDX%IMESH2) Then
            PDX%STR_INDX = i
            Do j = 1,STRP%N_NODES
               If ( Is_Within_Bounds(PDX%X1,PDX%X2,PDX%Y1,PDX%Y2,PDX%Z1,PDX%Z2, &
                    STRP%XB_NODE(j,1), STRP%XB_NODE(j,2), STRP%XB_NODE(j,3),STRP%XB_NODE(j,4), &
                    STRP%XB_NODE(j,5), STRP%XB_NODE(j,6), 0._EB, 0._EB, 0._EB)) Then
                  PDX%STR_SUB_INDX = j
                  Exit CheckDoorStrLoop
               End If
            End Do
         End If
       End Do CheckDoorStrLoop
    End Do

    ! Create list  of incoming nodes for STairs
    Do N = 1,N_STRS
      STRP => EVAC_STRS(N)
      STRP%N_NODES_IN = 0
      STRP%N_NODES_OUT = 0
      NODES_TMP = 0
      Do I = 1,N_NODES
         Select Case (EVAC_NODE_List(I)%Node_type)
         Case ('Door')
            J = EVAC_NODE_List(I)%Node_index
            If (EVAC_DOORS(J)%IMESH2 == STRP%IMESH) Then
               STRP%N_NODES_IN = STRP%N_NODES_IN + 1
               NODES_TMP(STRP%N_NODES_IN) = I
            End If
         Case ('Entry')
            J = EVAC_NODE_List(I)%Node_index
            If (EVAC_ENTRYS(J)%IMESH == STRP%IMESH) Then
               STRP%N_NODES_IN = STRP%N_NODES_IN + 1
               NODES_TMP(STRP%N_NODES_IN) = I
            End If
         End Select
      End Do
      Allocate(STRP%NODES_IN(1:STRP%N_NODES_IN),STAT=IZERO)
      Call ChkMemErr('Read_Evac','STRP%NODES_IN',IZERO) 
      STRP%NODES_IN = NODES_TMP(1:STRP%N_NODES_IN)
      ! Create List of outgoing nodes for stairs
      NODES_TMP = 0
      Do I = NMESHES+1,N_NODES
         If (EVAC_NODE_List(I)%Node_type == 'Entry') Cycle
         If (EVAC_NODE_List(I)%Node_type == 'Stair') Cycle
         If (EVAC_NODE_List(I)%IMESH == STRP%IMESH) Then
            STRP%N_NODES_OUT = STRP%N_NODES_OUT + 1
            NODES_TMP(STRP%N_NODES_OUT) = I
         End If
      End Do
      Allocate(STRP%NODES_OUT(1:STRP%N_NODES_OUT),STAT=IZERO)
      Call ChkMemErr('Read_Evac','STRP%NODES_OUT',IZERO) 
      STRP%NODES_OUT = NODES_TMP(1:STRP%N_NODES_OUT)
    End Do

    End Subroutine CHECK_EVAC_NODES

  End Subroutine READ_EVAC

  Subroutine Initialize_Evac_Dumps
    Implicit None
    !
    ! Local variables
    Character(50) tcform
    Integer n_cols, i, j, nm, izero
    Logical L_fed_read, L_fed_save, L_eff_read, L_eff_save, L_status
    Integer(4) n_egrids_tmp, ibar_tmp, jbar_tmp, kbar_tmp, &
         ntmp1, ntmp2, ntmp3, ntmp4, ntmp5, ntmp6, ios, N
    Real(FB) u_tmp, v_tmp
    Character(60), Allocatable, Dimension(:) :: CTEMP
    !
    Type (MESH_TYPE), Pointer :: MFF =>NULL()
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

    Write(EVAC_COMPILE_DATE,'(A)') evacrev(Index(evacrev,':')+1:Len_trim(evacrev)-2)
    Read (EVAC_COMPILE_DATE,'(I5)') EVAC_MODULE_REV
    Write(EVAC_COMPILE_DATE,'(A)') evacdate
    Call GET_REV_evac(EVAC_MODULE_REV,EVAC_COMPILE_DATE)
    !
    Write(LU_EVACOUT,'(/A)')          ' FDS+Evac Evacuation Module'
    Write(LU_EVACOUT,'(/A,A)')        ' FDS+Evac Compilation Date: ', &
         Trim(EVAC_COMPILE_DATE(Index(EVAC_COMPILE_DATE,'(')+1:Index(EVAC_COMPILE_DATE,')')-1))
    Write(LU_EVACOUT,'(A,A)')  ' FDS+Evac Version         : ', Trim(EVAC_VERSION)
    Write(LU_EVACOUT,'(A,i0)')  ' FDS+Evac SVN Revision No.: ', EVAC_MODULE_REV

    Write(LU_EVACOUT,fmt='(/a,i2)')  ' FDS+Evac Color_Method    :', COLOR_METHOD
    If (Fed_Door_Crit >= 0) Then
       Write(LU_EVACOUT,fmt='(a,f14.8)') ' FDS+Evac Fed_Door_Crit   :', FED_DOOR_CRIT
    Else
       ! Visibility S = 3/K, K is extinction coeff.
       Write(LU_EVACOUT,fmt='(a,f14.8,a)') ' FDS+Evac Vis_Door_Crit   :', Abs(FED_DOOR_CRIT), ' m'
    End If
    If (NOT_RANDOM ) Write(LU_EVACOUT,fmt='(a)') ' FDS+Evac Random seed is not used.'
    If (Fed_Door_Crit < 0) Then
       ! Visibility S = 3/K, K is extinction coeff.
       FED_DOOR_CRIT = 3.0_EB/FED_DOOR_CRIT ! Extinction coeff (1/m)
    End If
    !
    L_fed_read = Btest(I_EVAC,3)
    L_fed_save = Btest(I_EVAC,1)
    L_eff_read = Btest(I_EVAC,2)
    L_eff_save = Btest(I_EVAC,0)

    n_cols = n_egrids + n_corrs + n_exits + n_doors + 1 + n_exits - n_co_exits + n_doors
    ! Initialize the FED counters:
    icyc_old = -1
    n_dead = -1
    fed_max_alive = 0.0_EB
    fed_max = 0.0_EB
    !
    If (append) Then
       Open (LU_EVACCSV,file=FN_EVACCSV,form='formatted',status='old', position='append')
       !
       If (L_fed_save) Then
          Open (LU_EVACFED,file=FN_EVACFED,form='unformatted', status='old',position='append')
       End If
       ! 
       If (L_fed_read) Then
          Inquire (file=FN_EVACFED,exist=L_status)
          If (.Not. L_status) Then
             Write (LU_EVACOUT,fmt='(a,a,a)') ' FDS+Evac No FED File: ', Trim(FN_EVACFED), ', FED and soot not used'
             l_fed_read = .False.
             l_fed_save = .False.
             I_EVAC = Ibclr(I_EVAC,3)  ! do not read FED
             I_EVAC = Ibclr(I_EVAC,1)  ! do not save FED
          Else
             Call SHUTDOWN('ERROR: Evac Dumps: FED, no restart yet')
             Open (LU_EVACFED,file=FN_EVACFED,form='unformatted', status='old')
             Read (LU_EVACFED,Iostat=ios) n_egrids_tmp
             If (ios.Ne.0) Then
                Write(MESSAGE,'(A)') 'ERROR: Init Evac Dumps: FED READ ERROR'
                Close (LU_EVACFED)
                Call SHUTDOWN(MESSAGE)
             End If
             If (n_egrids_tmp /= n_egrids) Then
                Write(MESSAGE,'(A,2I4)') 'ERROR: Init Evac Dumps: FED ',n_egrids_tmp, n_egrids
                Call SHUTDOWN(MESSAGE)
             End If
          End If
       End If
       If (L_eff_read) Then
          Inquire (file=FN_EVACEFF,exist=L_status)
          If (L_status) Then
             Write (LU_EVACOUT,fmt='(a,a,a/)') ' FDS+Evac EFF File: ', Trim(FN_EVACEFF), ' is used'
             l_eff_save = .False.
             I_EVAC = Ibclr(I_EVAC,0)  ! do not save EFF
             Open (LU_EVACEFF,file=FN_EVACEFF,form='unformatted', status='old')
             Read (LU_EVACEFF,Iostat=ios) n_egrids_tmp
             If (ios.Ne.0) Then
                Write(MESSAGE,'(A)') 'ERROR: Init Evac Dumps: EFF READ ERROR'
                Close (LU_EVACEFF)
                Call SHUTDOWN(MESSAGE)
             End If
             If (n_egrids_tmp /= Count(EVACUATION_ONLY)) Then
                Write(MESSAGE,'(A,2I4)') 'ERROR: Init Evac Dumps: EFF ',n_egrids_tmp, Count(EVACUATION_ONLY)
                Close (LU_EVACEFF)
                Call SHUTDOWN(MESSAGE)
             End If
          Else
             Write(MESSAGE,'(A,2I4)') 'ERROR: Init Evac Dumps: EFF, no restart yet'
             Call SHUTDOWN(MESSAGE)
          End If
       End If
       !
       If ( l_fed_read .Or. l_fed_save ) n_dead = 0
       !

    Else                      ! replace files
       !
       If (L_fed_save) Then
          l_fed_read = .False.
          I_EVAC = Ibclr(I_EVAC,3)  ! do not read FED
          Open (LU_EVACFED,file=FN_EVACFED,form='unformatted', status='replace')
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
          Write (LU_EVACFED) ntmp1
          Write (LU_EVACFED) n_egrids_tmp, ntmp2, ntmp3, ntmp4, ntmp5, ntmp6
          Write (LU_EVACOUT,fmt='(a,a,a)') ' FDS+Evac FED File: ', Trim(FN_EVACFED), ' is calculated and used'
       End If
       ! 
       If (L_fed_read) Then
          Inquire (file=FN_EVACFED,exist=L_status)
          If (.Not. L_status) Then
             Write (LU_EVACOUT,fmt='(a,a,a)') ' FDS+Evac No FED File: ', Trim(FN_EVACFED), ', FED and soot not used'
             l_fed_read = .False.
             l_fed_save = .False.
             I_EVAC = Ibclr(I_EVAC,3)  ! do not read FED
             I_EVAC = Ibclr(I_EVAC,1)  ! do not save FED
          Else
             Open (LU_EVACFED,file=FN_EVACFED,form='unformatted', status='old')
             Read (LU_EVACFED,Iostat=ios) ntmp1
             If (ios.Ne.0) Then
                Write(MESSAGE,'(A)') 'ERROR: Init Evac Dumps: FED READ ERROR'
                Close (LU_EVACFED)
                Call SHUTDOWN(MESSAGE)
             End If
             If ( ntmp1 >= 0 ) Then
                ! Old format (version 1.10)
                n_egrids_tmp = ntmp1
             Else
                ! New format (version 1.11)
                Read (LU_EVACFED,Iostat=ios) n_egrids_tmp, ntmp2, ntmp3, ntmp4, ntmp5, ntmp6
                If (ios.Ne.0) Then
                   Write(MESSAGE,'(A)') 'ERROR: Init Evac Dumps: FED READ ERROR'
                   Close (LU_EVACFED)
                   Call SHUTDOWN(MESSAGE)
                End If
             End If

             ! Do not read old format. Do not read new format, if there the numbers
             ! are not: n_egrids, 4, n_corrs, 8 
             If ( ntmp2 /= 4 .Or. ntmp3 /= n_corrs .Or. ntmp1 >= 0 .Or. ntmp4 /= 8  .Or. &
                  ntmp5 /= n_doors+n_exits .Or. ntmp6 /= 4) Then
                Write (LU_EVACOUT,fmt='(a,a,a)') ' FDS+Evac Error in FED File: ', Trim(FN_EVACFED), ', FED and soot not used'
                l_fed_read = .False.
                l_fed_save = .False.
                I_EVAC = Ibclr(I_EVAC,3) ! do not read FED
                I_EVAC = Ibclr(I_EVAC,1) ! do not save FED
                Close (LU_EVACFED)
             End If
             If (n_egrids_tmp /= n_egrids) Then
                Write(MESSAGE,'(A,2I4)') 'ERROR: Init Evac Dumps: FED ',n_egrids_tmp, n_egrids
                Close (LU_EVACFED)
                Call SHUTDOWN(MESSAGE)
             End If

             If (l_fed_read .Or. l_fed_save) Then
                Write (LU_EVACOUT,fmt='(a,a,a)') ' FDS+Evac FED File: ', Trim(FN_EVACFED), ' is used'
             End If
          End If
       End If
       ! 
       ! Number of evac flow fields is same as the number of all evac grids.
       If (L_eff_read) Then
          ios = 3
          Inquire (file=FN_EVACEFF,exist=L_status)
          If (L_status) Then
             l_eff_save = .False.
             l_eff_read = .True.
             I_EVAC = Ibclr(I_EVAC,0) ! do not save EFF
             I_EVAC = Ibset(I_EVAC,2) ! read EFF
             Open (LU_EVACEFF,file=FN_EVACEFF,form='unformatted', status='old')
             Read (LU_EVACEFF,Iostat=ios) n_egrids_tmp
             If (ios.Ne.0) Then
                ios = 1
                Write(LU_EVACOUT,'(A)') ' WARNING: Init Evac Dumps: EFF READ ERROR'
                Write(LU_EVACOUT,'(A)') ' WARNING: EFF file is not read in'
                Close (LU_EVACEFF)
             End If
             If (n_egrids_tmp /= Count(EVACUATION_ONLY) .And. ios < 1) Then
                ios = 2
                Write(LU_EVACOUT,'(A,2I4)') ' WARNING: Init Evac Dumps: EFF READ ERROR ',n_egrids_tmp, Count(EVACUATION_ONLY)
                Write(LU_EVACOUT,'(A)')     ' WARNING: EFF file is not read in'
                Close (LU_EVACEFF)
             End If
          End If
          If (ios .Ne. 0) Then
             ! Read error ==> recalculate EFF
             l_eff_save = .True.
             l_eff_read = .False.
             I_EVAC = Ibclr(I_EVAC,2) ! do not read EFF
             I_EVAC = Ibset(I_EVAC,0) ! save EFF
          End If
       End If  ! L_eff_read
       !
       Allocate(CTEMP(Max(1,n_exits-n_co_exits)), STAT = IZERO)
       Call ChkMemErr('Initialize_Evac_Dumps','CTEMP', IZERO)
       j = 0
       Do i = 1, n_exits
          If (.Not. EVAC_EXITS(i)%COUNT_ONLY) Then
             j = j + 1
             CTEMP(j) = Trim(EVAC_EXITS(i)%ID)
          End If
       End Do
    
       If ( l_fed_read .Or. l_fed_save ) Then
          ! Write the 'fed' columns
          n_dead = 0
          Open (LU_EVACCSV,file=FN_EVACCSV,form='formatted',status='replace')
          ! June 2009: Changed the .csv file format to the fds5 style
          ! first row: units (or variable class)
          ! second row: variable name
          ! third row-: data
          ! Write (LU_EVACCSV,*) n_cols+3
          Write (tcform,'(a,i4.4,a)') "(",n_cols+3,"(a,','),a)"
          ! Write (LU_EVACCSV,tcform) 'Time','Humans', &
          !      ('Floor', i=1,n_egrids), &
          !      ('Corridor', i=1,n_corrs), &
          !      ('Exit', i=1,n_exits), &
          !      ('Door', i=1,n_doors), &
          !      ('Exit', i=1,n_exits-n_co_exits), &
          !      ('Door', i=1,n_doors), &
          !      'Fed','Fed','Fed'
          Write (LU_EVACCSV,tcform) 's','AgentsInside', &
               ('AgentsInsideMesh', i=1,n_egrids), &
               ('AgentsInsideCorr', i=1,n_corrs), &
               ('ExitCounter', i=1,n_exits), &
               ('DoorCounter', i=1,n_doors), &
               ('TargetExitCounter', i=1,n_exits-n_co_exits), &
               ('TargetDoorCounter', i=1,n_doors), &
               'Agents','FED_Index','FED_Index'
          Write (LU_EVACCSV,tcform) 'EVAC_Time','AllAgents', &
               (Trim(EVAC_Node_List(i)%GRID_NAME), i=1,n_egrids), &
               (Trim(EVAC_CORRS(i)%ID), i=1,n_corrs), &
               (Trim(EVAC_EXITS(i)%ID), i=1,n_exits), &
               (Trim(EVAC_DOORS(i)%ID), i=1,n_doors), &
               (Trim(CTEMP(i)), i=1,n_exits-n_co_exits), &
               (Trim(EVAC_DOORS(i)%ID), i=1,n_doors), &
               'Number_of_Deads','FED_max','FED_max_alive'
       Else
          ! Do not write the 'fed' columns
          Open (LU_EVACCSV,file=FN_EVACCSV,form='formatted',status='replace')
          ! June 2009: Changed the .csv file format to the fds5 style
          ! first row: units (or variable class)
          ! second row: variable name
          ! third row-: data
          ! Write (LU_EVACCSV,*) n_cols
          Write (tcform,'(a,i4.4,a)') "(",n_cols,"(a,','),a)"
          ! Write (LU_EVACCSV,tcform) 'Time','Humans', &
          !      ('Floor', i=1,n_egrids), &
          !      ('Corridor', i=1,n_corrs), &
          !      ('Exit', i=1,n_exits), &
          !      ('Door', i=1,n_doors), &
          !      ('Exit', i=1,n_exits-n_co_exits), &
          !      ('Door', i=1,n_doors)
          Write (LU_EVACCSV,tcform) 's','AgentsInside', &
               ('AgentsInsideMesh', i=1,n_egrids), &
               ('AgentsInsideCorr', i=1,n_corrs), &
               ('ExitCounter', i=1,n_exits), &
               ('DoorCounter', i=1,n_doors), &
               ('TargetExitCounter', i=1,n_exits-n_co_exits), &
               ('TargetDoorCounter', i=1,n_doors)
          Write (LU_EVACCSV,tcform) 'EVAC_Time','AllAgents', &
               (Trim(EVAC_Node_List(i)%GRID_NAME), i=1,n_egrids), &
               (Trim(EVAC_CORRS(i)%ID), i=1,n_corrs), &
               (Trim(EVAC_EXITS(i)%ID), i=1,n_exits), &
               (Trim(EVAC_DOORS(i)%ID), i=1,n_doors), &
               (Trim(CTEMP(i)), i=1,n_exits-n_co_exits), &
               (Trim(EVAC_DOORS(i)%ID), i=1,n_doors)
       End If
       Deallocate(CTEMP)
    End If                  ! if-append-else
    ! 
    ! Read the evac flow fields from the disk, if they exist.
    If ( L_eff_read ) Then
       ios = 0
       ReadEffLoop: Do nm = 1, NMESHES
          If (EVACUATION_ONLY(NM)) Then
             MFF=>MESHES(nm)
             Read (LU_EVACEFF,Iostat=ios) ibar_tmp, jbar_tmp, kbar_tmp
             If (ios .Ne. 0) Then
                Write (LU_EVACOUT,'(A)') ' WARNING: Init Evac Dumps: EFF READ ERROR'
                Write (LU_EVACOUT,'(A)') ' WARNING: EFF file is not read in'
                Close (LU_EVACEFF)
                Exit ReadEffLoop
             End If
             If ( MFF%IBAR /= ibar_tmp .Or. MFF%JBAR /= jbar_tmp .Or. MFF%KBAR /= kbar_tmp ) Then
                ios = 2
                Write (LU_EVACOUT,'(A)') ' WARNING: Init Evac Dumps: EFF READ ERROR'
                Write (LU_EVACOUT,'(A)') ' WARNING: EFF file is not read in'
                Close (LU_EVACEFF)
                Exit ReadEffLoop
             End If
             Do  i = 0, MFF%IBAR+1
                Do j= 0, MFF%JBAR+1
                   Read (LU_EVACEFF,Iostat=ios) u_tmp, v_tmp
                   If (ios .Ne. 0) Then
                      Write (LU_EVACOUT,'(A)') ' WARNING: Init Evac Dumps: EFF READ ERROR'
                      Write (LU_EVACOUT,'(A)') ' WARNING: EFF file is not read in'
                      Close (LU_EVACEFF)
                      Exit ReadEffLoop
                   End If
                   MFF%U(i,j,:) = u_tmp
                   MFF%V(i,j,:) = v_tmp
                   MFF%W(i,j,:) = 0.0_EB
                End Do
             End Do
             MFF%UVW_GHOST(:,1)=-1.E6_EB
             MFF%UVW_GHOST(:,2)=-1.E6_EB
          End If
       End Do ReadEffLoop
       If (ios .Ne. 0) Then
          ! Read error ==> recalculate EFF
          l_eff_save = .True.
          l_eff_read = .False.
          I_EVAC = Ibclr(I_EVAC,2) ! do not read EFF
          I_EVAC = Ibset(I_EVAC,0) ! save EFF
       End If
    End If
    !
    If (L_eff_read) Then
       Write (LU_EVACOUT,fmt='(a,a,a/)') ' FDS+Evac EFF File: ', Trim(FN_EVACEFF), ' is read in and used'
    End If
    If (L_eff_save .And. .Not. append) Then
       l_eff_read = .False.
       I_EVAC = Ibclr(I_EVAC,2)  ! do not read EFF
       Open (LU_EVACEFF,file=FN_EVACEFF,form='unformatted', status='replace')
       n_egrids_tmp = Count(EVACUATION_ONLY)
       Write (LU_EVACEFF) n_egrids_tmp
       Write (LU_EVACOUT,fmt='(a,a,a/)') ' FDS+Evac EFF File: ', Trim(FN_EVACEFF), ' is (re)calculated and used'
    End If

    EVAC_Z_MIN =  Huge(EVAC_Z_MIN)
    EVAC_Z_MAX = -Huge(EVAC_Z_MIN)
    Do nm = 1, NMESHES
       MFF=>MESHES(nm)
       EVAC_Z_MIN = Min(EVAC_Z_MIN,Real(MFF%ZS,FB))
       EVAC_Z_MAX = Max(EVAC_Z_MAX,Real(MFF%ZF,FB))
    End Do

    ! write STRS properties
    If (MYID==Max(0,EVAC_PROCESS)) Then
       Do N = 1, N_STRS
          Write (LU_EVACOUT,'(A,A)')      '  Stair ',Trim(EVAC_STRS(N)%ID)
          Write (LU_EVACOUT,'(A,6F10.3)') '   Co-ordinates: ',EVAC_STRS(N)%XB(1:6)
          Write (LU_EVACOUT,'(A)')          '   Node coordinates'
          Do NM = 1,EVAC_STRS(N)%N_NODES
             If (EVAC_STRS(N)%NODE_TYPE(NM) == STRS_LANDING_TYPE) Then
                Write (LU_EVACOUT,'(A,I3,8F8.3)') '   Landing ',NM,EVAC_STRS(N)%XB_NODE(NM,1:8)
             Else
                Write (LU_EVACOUT,'(A,I3,8F8.3)') '   Stair   ',NM,EVAC_STRS(N)%XB_NODE(NM,1:8)
             Endif
          ENDDO
          Write (LU_EVACOUT,'(A)')          '   Nodes in '
          DO NM = 1, EVAC_STRS(N)%N_NODES_IN
             Write (LU_EVACOUT,'(I5,A,A)')          NM, ' ', EVAC_NODE_List(EVAC_STRS(N)%NODES_IN(NM))%ID
          ENDDO
          Write (LU_EVACOUT,'(A)')          '   Nodes out '
          DO NM = 1, EVAC_STRS(N)%N_NODES_OUT
             Write (LU_EVACOUT,'(I5,A,A)')          NM, ' ', EVAC_NODE_List(EVAC_STRS(N)%NODES_OUT(NM))%ID
          ENDDO
       Enddo
   End if

  End Subroutine Initialize_Evac_Dumps
      
!
  Subroutine INITIALIZE_EVACUATION(NM,ISTOP)
    Implicit None
    !
    ! Insert humans into the domain at the start of calculation
    !
    ! Passed variables
    Integer, Intent(IN) :: NM
    Integer, Intent(OUT) :: ISTOP
    !
    ! Local variables
    Real(EB) RN, simoDX, simoDY, TNOW
    Real(EB) VOL1, VOL2, X1, X2, Y1, Y2, Z1, Z2, &
         dist, d_max, G_mean, G_sd, G_high, G_low, x1_old, y1_old
    Integer i,j,k,ii,jj,kk,ipc, izero, n_tmp, ie, nom, I_OBST
    Integer i11, i22, group_size
    Logical pp_see_group, is_solid
    Integer jjj, iii, i44
    Real(EB) x11, y11, group_x_sum, group_y_sum, group_x_center, group_y_center, dens_fac
    Integer :: i_endless_loop, istat
    Real(EB), Dimension(6) :: y_tmp, x_tmp, r_tmp
    Real(EB), Dimension(4) :: d_xy
    Logical, Dimension(4) :: FoundWall_xy
    ! 
    Type (MESH_TYPE), Pointer :: M =>NULL()
    Type (EVAC_SSTAND_Type),Pointer :: ESS=>NULL()
    Type (EVACUATION_Type), Pointer :: HPT=>NULL()
    Type (EVAC_PERS_Type),  Pointer :: PCP=>NULL()
    Type (EVAC_HOLE_Type),  Pointer :: EHX=>NULL()
    Type (HUMAN_TYPE), Pointer :: HR=>NULL(), HRE=>NULL()
    !
    TNOW = SECOND()

    If ( .Not.(EVACUATION_ONLY(NM) .And. EVACUATION_GRID(NM)) ) Return
    ! Next means that only EVAC_PROCESS is doing something
    If (MYID /= PROCESS(NM)) Return
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
    Allocate(MESHES(NM)%HUMAN(MESHES(NM)%N_HUMANS_DIM),STAT=IZERO)
    Call ChkMemErr('INIT_EVACUATION','HUMAN',IZERO)
    !
    ! HUMAN_GRID: (x,y,z): Center of the grid cells
    !             SOOT_DENS: Smoke density at the center of the grid cells (mg/m3)
    !             FED_CO_CO2_O2: FED index
    !             IMESH: Which fire mesh, if any
    !             II,JJ,KK: Fire mesh indexes
    Allocate(MESHES(NM)%HUMAN_GRID(MESHES(NM)%IBAR,MESHES(NM)%JBAR), STAT=IZERO)
    Call ChkMemErr('INIT_EVACUATION','HUMAN_GRID',IZERO)
    !
    Call POINT_TO_MESH(NM)
    !
    ! Initialise HUMAN_GRID
    !
    FED_I_LOOP: Do i = 1,IBAR
       FED_J_LOOP: Do j= 1,JBAR
          x1 = XS + (i-1)*DXI + 0.5_EB*DXI
          y1 = YS + (j-1)*DETA + 0.5_EB*DETA
          ! z1 is here the head height where smoke information is taken, so it
          ! should be about the head height
          z1 = 0.5_EB*(ZF+ZS) - EVACUATION_Z_OFFSET(NM) + HUMAN_SMOKE_HEIGHT
          SS_Loop: Do k = 1, n_sstands
             ! Inclines, i.e., EVSS, one should take the smoke at the correct height.
             ESS => EVAC_SSTANDS(k)
             If (ESS%IMESH == nm .And. (ESS%X1 <= x1 .And. ESS%X2 >= x1) .And. &
                  (ESS%Y1 <= y1 .And. ESS%Y2 >= y1) ) Then
                Select Case (ESS%IOR)
                Case(-1)
                   z1 = z1 + ESS%H0 + (ESS%H-ESS%H0)*Abs(ESS%X1-x1)/Abs(ESS%X1-ESS%X2)
                Case(+1)
                   z1 = z1 + ESS%H0 + (ESS%H-ESS%H0)*Abs(ESS%X2-x1)/Abs(ESS%X1-ESS%X2)
                Case(-2)
                   z1 = z1 + ESS%H0 + (ESS%H-ESS%H0)*Abs(ESS%Y1-y1)/Abs(ESS%Y1-ESS%Y2)
                Case(+2)
                   z1 = z1 + ESS%H0 + (ESS%H-ESS%H0)*Abs(ESS%Y2-y1)/Abs(ESS%Y1-ESS%Y2)
                End Select
             End If
             Exit SS_Loop
          End Do SS_Loop
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
          If (.Not. Btest(I_EVAC,4) ) Cycle FED_J_LOOP

          MESH_LOOP: Do NOM = 1, NMESHES
             If (.Not. EVACUATION_ONLY(NOM)) Then
                M => MESHES(NOM)
                If ( X1 >= M%XS .And. X1 <= M%XF .And. &
                     Y1 >= M%YS .And. Y1 <= M%YF .And. &
                     Z1 >= M%ZS .And. Z1 <= M%ZF) Then
                   II = Floor( M%CELLSI(Floor((X1-M%XS)*M%RDXINT)) + 1.0_EB  )
                   JJ = Floor( M%CELLSJ(Floor((Y1-M%YS)*M%RDYINT)) + 1.0_EB  )
                   KK = Floor( M%CELLSK(Floor((Z1-M%ZS)*M%RDZINT)) + 1.0_EB  )
                   I_OBST = OBST_INDEX_C(M%CELL_INDEX(II,JJ,KK))
                   If ( M%SOLID(M%CELL_INDEX(II,JJ,KK)) .AND. .NOT.M%OBSTRUCTION(I_OBST)%HIDDEN ) Then
                      HUMAN_GRID(i,j)%IMESH = 0 ! No smoke inside OBSTs
                   Else
                      HUMAN_GRID(i,j)%II = II
                      HUMAN_GRID(i,j)%JJ = JJ
                      HUMAN_GRID(i,j)%KK = KK
                      HUMAN_GRID(i,j)%IMESH = NOM
                   End If
                   Exit MESH_LOOP
                End If
             End If
             ! No fire mesh is found
             HUMAN_GRID(i,j)%IMESH = 0
          End Do MESH_LOOP
          I_OBST = OBST_INDEX_C(CELL_INDEX(I,J,1))
          If (.Not. (SOLID(CELL_INDEX(i,j,1)) .AND. .NOT. OBSTRUCTION(I_OBST)%HIDDEN) ) Then
             HUMAN_GRID(i,j)%IMESH = HUMAN_GRID(i,j)%IMESH
          Else
             ! This grid cell is solid ==> No humans in this cell
             ! Zero, if fire mesh obst or no fire mesh at all
             ! Negative (-nom), if fire mesh gas cell but evac mesh solid.
             HUMAN_GRID(i,j)%IMESH = -HUMAN_GRID(i,j)%IMESH
          End If
       End Do FED_J_LOOP
    End Do FED_I_LOOP

    EVAC_CLASS_LOOP: Do IPC=1,NPC_EVAC
       !
       HPT=>EVACUATION(IPC)
       !
       ! Check the mesh
       If (HPT%IMESH /= NM) Cycle EVAC_CLASS_LOOP
       !
       ! If there is an initial number of humans, initialize
       If (HPT%N_INITIAL == 0) Cycle EVAC_CLASS_LOOP
       ! 
       If (HPT%X1 == 0.0_EB .And. HPT%X2 == 0.0_EB .And. &
            HPT%Y1 == 0.0_EB .And. HPT%Y2 == 0.0_EB .And. &
            HPT%Z1 == 0.0_EB .And. HPT%Z2 == 0.0_EB ) Then
          X1 = XS ; X2 = XF
          Y1 = YS ; Y2 = YF
          Z1 = ZS ; Z2 = ZF
          VOL2 = (XF - XS) * (YF - YS)
          VOL1 = VOL2
       Else
          If (HPT%X1 > XF .Or. HPT%X2 < XS .Or. &
               HPT%Y1 > YF .Or. HPT%Y2 < YS .Or. &
               HPT%Z1 > ZF .Or. HPT%Z2 < ZS) Then
             Cycle EVAC_CLASS_LOOP
          End If
          X1 = Max(HPT%X1,XS) ; X2 = Min(HPT%X2, XF)
          Y1 = Max(HPT%Y1,YS) ; Y2 = Min(HPT%Y2, YF)
          Z1 = Max(HPT%Z1,ZS) ; Z2 = Min(HPT%Z2, ZF)
          VOL2 = (HPT%X2 - HPT%X1) * (HPT%Y2 - HPT%Y1)
          VOL1 = (X2 - X1) * (Y2 - Y1) * (Z2 - Z1)
       End If
       ! 
       ! Check which evacuation floor node  (=1,...,n_egrids)
       n_tmp = 0
       HP_MeshLoop: Do i = 1, nmeshes
          If (evacuation_only(i) .And. evacuation_grid(i)) Then
             n_tmp = n_tmp +1
             If (HPT%IMESH == i) Then
                Exit HP_MeshLoop
             End If
          End If
       End Do HP_MeshLoop
       If (n_tmp < 1 .Or. n_tmp > n_egrids) Then
          Write(MESSAGE,'(A,A,A,I4)') 'ERROR: INIT_EVAC ',Trim(HPT%ID),' problem evac node, INODE= ',n_tmp
          Call SHUTDOWN(MESSAGE)
       End If
       !
       PCP => EVAC_PERSON_CLASSES(HPT%IPC)
       !
       ! i11: counter, humans per EVAC-line
       ! i22: group member index
       ! i33: group index (=0: no group, i.e., single humans)
       ! ilh: lonely human index
       !
       i11 = 0

       INITIALIZATION_LOOP: Do I=1,HPT%N_INITIAL
          !
          Call Random_number(RN)
          group_size = HPT%GN_MIN - 1 +  Int((HPT%GN_MAX-HPT%GN_MIN+1)*RN+0.5_EB)
          group_size = Max(Min(group_size, HPT%GN_MAX),HPT%GN_MIN)

          If ( i11+group_size > HPT%N_INITIAL ) Then
             group_size = HPT%N_INITIAL - i11
          End If

          If ( i11 >= HPT%N_INITIAL ) Exit INITIALIZATION_LOOP

          i22 = 0
          i_endless_loop = 0
          group_X_sum = 0
          group_Y_sum = 0
          GROUP_SIZE_LOOP: Do 
             ! i22: Counter, humans on this group (group index i33)
             i22 = i22 + 1
             If (i22 > group_size) Exit GROUP_SIZE_LOOP
             i11 = i11 + 1

             If (i22 == 1) Then
                ! One member groups are not counted as group
                If (group_size > 1) i33 = i33 + 1
                If (group_size == 1) ilh = ilh + 1
             End If
             N_HUMANS = N_HUMANS + 1
             !
             If (N_HUMANS > N_HUMANS_DIM) Then
                Call SHUTDOWN('ERROR: Init Humans: no re-allocation yet')
                Call RE_ALLOCATE_HUMANS(1,NM)
                HUMAN=>MESHES(NM)%HUMAN
             End If
             !
             HR => HUMAN(N_HUMANS)
             Call CLASS_PROPERTIES(HR,PCP)
             HR%IPC = HPT%IPC  ! PERS-line index
             HR%IEL = IPC      ! EVAC-line index
             HR%I_Target = 0

             !
             BLK_LOOP:  Do
                If (i22 == 1) Then
                   ! First member of the group
                   Call Random_number(RN)
                   HR%X = X1 + RN*(X2-X1)
                   x1_old = HR%X
                   Call Random_number(RN)
                   HR%Y = Y1 + RN*(Y2-Y1)
                   y1_old = HR%Y
                Else
                   ! Next members of the group are put around the first member.
                   G_mean = 0.0_EB
                   G_sd   = 4.0_EB  ! units (m) std.dev.
                   G_high = 6.0_EB  ! units (m) cut-off
                   G_low  = -6.0_EB ! units (m) cut-off
                   ! First the angle, then the radial distance
                   Call Random_number(rn)
                   simoDX = Sin(2.0_EB*Pi*rn)
                   simoDY = Cos(2.0_EB*Pi*rn)
                   G_mean = (2.0_EB/3.0_EB)*Sqrt(group_size/(Pi*GROUP_DENS))
                   G_sd   = G_mean  ! units (m) std.dev.
                   G_high = Max(3.0_EB,3.0_EB/GROUP_DENS)* G_mean ! units (m) cut-off
                   G_low  = 0.25_EB      ! units (m) cut-off
                   GTrunFlag=0
                   rn = GaussTrun(G_mean,G_sd,G_low,G_high)
                   simoDX = rn*simoDX
                   simoDY = rn*simoDY
                   HR%X = Min(X2,Max(X1, x1_old + simoDX))
                   HR%Y = Min(Y2,Max(Y1, y1_old + simoDY))
                End If
                HR%Z = Z1 + 0.5_EB*(Z2-Z1)
                If (HPT%Angle > -999.9_EB) Then
                   HR%Angle = HPT%Angle
                Else
                   Call Random_number(RN)
                   HR%Angle = 2.0_EB*Pi*rn
                End If
                Do While (HR%Angle >= 2.0_EB*Pi)
                   HR%Angle = HR%Angle - 2.0_EB*Pi
                End Do
                Do While (HR%Angle < 0.0_EB)
                   HR%Angle = HR%Angle + 2.0_EB*Pi
                End Do

                ! EVHO (evacuation hole) checking
                EH_Loop: Do ie = 1, n_holes
                   EHX => EVAC_HOLES(ie)
                   If (EHX%IMESH /= NM) Cycle EH_Loop
                   If ( Trim(EHX%EVAC_ID) /= 'null' .And. Trim(EHX%EVAC_ID) /= Trim(HPT%ID)) Cycle EH_Loop
                   If ( Trim(EHX%PERS_ID) /= 'null' .And. Trim(EHX%PERS_ID) /= Trim(HPT%CLASS_NAME)) Cycle EH_Loop
                   If ( (EHX%Y1 <= HR%Y .And. EHX%Y2 >= HR%Y) .And. (EHX%X1 <= HR%X .And. EHX%X2 >= HR%X) ) Then
                      ! User should not give too large EVHO:s
                      i_endless_loop = i_endless_loop + 1
                      Cycle BLK_LOOP
                   End If
                End Do EH_Loop

                !Check, that a person is not put on top of some other person
                If (DENS_INIT > 2.0_EB) Then
                   ! High density is wanted
                   d_max = 0.0_EB
                Else
                   d_max = 1.0_EB*HR%B
                End If
                dens_fac = Max(1.0_EB,DENS_INIT)

                If (i_endless_loop >= 10*Int(dens_fac*(16.0_EB*Max(1.0_EB,Log10(2.5_EB*VOL2))) / &
                     Max(1.0_EB,Log10((2.5_EB*VOL2)/(2.5_EB*VOL2-1)))) ) Then
                   Write (LU_EVACOUT,fmt='(A,I4,A,I4,A,I6)') ' ERROR: Initialize_Humans, EVAC line ', &
                        IPC, ', Mesh ', NM, ', i_human ', n_humans
                   Write (LU_EVACOUT,fmt='(a)') '      x       y       z     Rd      Rt      Rs      ds  '
                   Write (LU_EVACOUT,fmt='(3f8.2,4f8.4)') HR%X, HR%Y, HR%Z, &
                        2.0_EB*HR%Radius, HR%r_torso, HR%r_shoulder, HR%d_shoulder
                   ISTOP = 3  ! Stop code: FDS improperly set-up
                   HR%SHOW = .True.    
                   HR%COLOR_INDEX = 7  ! Cyan
                   Exit INITIALIZATION_LOOP
                End If

                ! Put a human in a mesh, if there is enough empty space, i.e.,
                ! check the OBSTs.
                Is_Solid = .False.
                KK = 1

                r_tmp(1) = HR%r_shoulder ! right circle
                r_tmp(2) = HR%r_torso    ! center circle
                r_tmp(3) = HR%r_shoulder ! left circle
                y_tmp(1) = HR%Y - Cos(HR%angle)*HR%d_shoulder ! right
                x_tmp(1) = HR%X + Sin(HR%angle)*HR%d_shoulder
                y_tmp(2) = HR%Y ! torso
                x_tmp(2) = HR%X
                y_tmp(3) = HR%Y + Cos(HR%angle)*HR%d_shoulder ! left
                x_tmp(3) = HR%X - Sin(HR%angle)*HR%d_shoulder
                II = Floor( CELLSI(Floor((x_tmp(1)-XS)*RDXINT)) + 1.0_EB )
                JJ = Floor( CELLSJ(Floor((y_tmp(1)-YS)*RDYINT)) + 1.0_EB )
                I_OBST = OBST_INDEX_C(CELL_INDEX(II,JJ,KK))
                Is_Solid = (Is_Solid .Or. (SOLID(CELL_INDEX(II,JJ,KK)) .AND. .NOT. OBSTRUCTION(I_OBST)%HIDDEN))
                II = Floor( CELLSI(Floor((x_tmp(3)-XS)*RDXINT)) + 1.0_EB )
                JJ = Floor( CELLSJ(Floor((y_tmp(3)-YS)*RDYINT)) + 1.0_EB )
                Is_Solid = (Is_Solid .Or. (SOLID(CELL_INDEX(II,JJ,KK)) .AND. .NOT. OBSTRUCTION(I_OBST)%HIDDEN))
                II = Floor( CELLSI(Floor((x_tmp(2)-XS)*RDXINT)) + 1.0_EB )
                JJ = Floor( CELLSJ(Floor((y_tmp(2)-YS)*RDYINT)) + 1.0_EB )
                Is_Solid = (Is_Solid .Or. (SOLID(CELL_INDEX(II,JJ,KK)) .AND. .NOT. OBSTRUCTION(I_OBST)%HIDDEN))

                If (.Not.Is_Solid) Then
                   ! Check the distances to walls (not to wall corners)
                   ! Skip wall force ior = 0: Check all walls (do not put too
                   ! close to doors/exits)
                   Call Find_walls(nm, x_tmp(1), y_tmp(1), r_tmp(1), d_max, 0, d_xy, FoundWall_xy, istat)
                   If (istat/=0) Is_Solid = .True.
                   Call Find_walls(nm, x_tmp(2), y_tmp(2), r_tmp(2), d_max, 0, d_xy, FoundWall_xy, istat)
                   If (istat/=0) Is_Solid = .True.
                   Call Find_walls(nm, x_tmp(3), y_tmp(3), r_tmp(3), d_max, 0, d_xy, FoundWall_xy, istat)
                   If (istat/=0) Is_Solid = .True.
                End If

                If (.Not.Is_Solid) Then
                   ! Check that the agent is not too close to other agents, who
                   ! already introduced in the calculation.
                   P2PLoop: Do ie = 1, n_humans - 1
                      HRE => HUMAN(ie)
                      r_tmp(4) = HRE%r_shoulder ! right circle
                      r_tmp(5) = HRE%r_torso    ! center circle
                      r_tmp(6) = HRE%r_shoulder ! left circle
                      y_tmp(4) = HRE%Y - Cos(HRE%angle)*HRE%d_shoulder ! right circle
                      x_tmp(4) = HRE%X + Sin(HRE%angle)*HRE%d_shoulder
                      y_tmp(5) = HRE%Y ! center circle
                      x_tmp(5) = HRE%X
                      y_tmp(6) = HRE%Y + Cos(HRE%angle)*HRE%d_shoulder ! left circle
                      x_tmp(6) = HRE%X - Sin(HRE%angle)*HRE%d_shoulder
                      Do iii = 1, 3
                         Do jjj = 4, 6
                            DIST = Sqrt((x_tmp(jjj)-x_tmp(iii))**2 + (y_tmp(jjj)-y_tmp(iii))**2) - (r_tmp(jjj)+r_tmp(iii))
                            If ( DIST < d_max) Then
                               i_endless_loop = i_endless_loop + 1
                               Cycle BLK_LOOP
                            End If
                         End Do
                      End Do
                   End Do P2PLoop

                   If (i22 > 1) Then
                      ! Check if the new member will see the first member of the group.
                      X11 = HR%X
                      Y11 = HR%Y

                      PP_see_group = See_each_other(nm, x11, y11, x1_old, y1_old)
                      
                   Else
                      PP_see_group = .True.
                   End If

                   If ( .Not. PP_see_group ) Then 
                      i_endless_loop = i_endless_loop + 1
                      Cycle BLK_LOOP
                   End If

                   ! Coordinates are OK, exit the loop
                   ! Only the centres of the circles are tested agains solids.
                   Exit BLK_LOOP
                Else
                   i_endless_loop = i_endless_loop + 1
                End If            ! Not inside a solid object
             End Do BLK_LOOP
             ! 
             ! The coordinates of the new person has passed the tests.
             !
             ILABEL_last = ILABEL_last + 1
             HR%ILABEL = ILABEL_last
             If (group_size > 1 ) Then
                ! Is a member of a group
                HR%GROUP_ID = i33
             Else
                ! Is a lonely soul
                HR%GROUP_ID = -ilh
             End If
             HR%Commitment = 0.0_EB
             HR%SHOW = .True.    


             Select Case (COLOR_METHOD)
             Case (-1)
                HR%COLOR_INDEX = 1
             Case (0)
                HR%COLOR_INDEX = HPT%Avatar_Color_Index
             Case (1)
                HR%COLOR_INDEX = Mod(group_size-1,6) + 1
             Case (2)
                If (HR%GROUP_ID > 0 ) Then
                   HR%COLOR_INDEX = Mod(HR%GROUP_ID,6) + 1
                Else
                   HR%COLOR_INDEX = 1 ! lonely human
                End If
             Case (3)
                HR%COLOR_INDEX = evac_person_classes(HPT%IPC)%Avatar_Color_Index
             Case (4)
                HR%COLOR_INDEX = 1
             Case (5)
                HR%COLOR_INDEX = 1
             Case Default
                Write(MESSAGE,'(A,I3,A)') 'ERROR: READ_EVAC COLOR METHOD',COLOR_METHOD, ' is not defined'
                Call SHUTDOWN(MESSAGE)
             End Select

             HR%IMESH       = HPT%IMESH
             HR%INODE       = n_tmp
             HR%NODE_NAME   = Trim(EVAC_Node_List(n_tmp)%ID)
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
             HR%NewRnd    = .True.
             HR%STR_SUB_INDX = 0
             HR%SKIP_WALL_FORCE_IOR = 0



             group_X_sum = group_X_sum + HR%X
             group_Y_sum = group_Y_sum + HR%Y
             If ( group_size == 1 ) Exit  GROUP_SIZE_LOOP
             If ( i22 == group_size .And. group_size > 1 ) Then
                group_X_center = group_X_sum / Max(1,group_size)
                group_Y_center = group_Y_sum / Max(1,group_size)
                !
                II = Floor( CELLSI(Floor((group_X_center-XS)*RDXINT)) + 1.0_EB )
                JJ = Floor( CELLSJ(Floor((group_Y_center-YS)*RDYINT)) + 1.0_EB )
                KK = 1
                x1_old = group_X_center
                y1_old = group_Y_center
                PP_see_group = .True.
                Do i44 = 1,group_size
                   HR => HUMAN( N_HUMANS - group_size + i44 )
                   X11 = HR%X
                   Y11 = HR%Y
                   PP_see_group = See_each_other(nm, x11, y11, x1_old, y1_old)

                   If ( .Not. PP_see_group ) Then 
                      i_endless_loop = i_endless_loop + 1
                      i22 = 0  ! Start at the beginning of the group
                      i33 = i33 - 1  ! group index
                      i11 = i11 - group_size ! human index (evac-line)
                      N_HUMANS = N_HUMANS - group_size ! total # of humans
                      ILABEL_last = ILABEL_last - group_size ! labels of humans
                      Cycle GROUP_SIZE_LOOP
                   End If
                End Do

                Exit GROUP_SIZE_LOOP

             End If ! group_size>1 and i22=group_size

          End Do GROUP_SIZE_LOOP ! i22, group_size loop
          !
       End Do INITIALIZATION_LOOP  ! i, # persons in a evac-line
       ! 
    End Do EVAC_CLASS_LOOP ! ipc, number of evac-lines

    Write (LU_EVACOUT,fmt='(a,f8.2,a,i0,a,i0/)') ' EVAC: Time ', 0.0_EB,' mesh ',nm,' number of humans ',n_humans
    ! Write (LU_OUTPUT,fmt='(a,f8.2,a,i0,a,i0/)') ' EVAC: Time ', &
    !      0.0_EB,' mesh ',nm,' number of humans ',n_humans
    !
    TUSED(12,NM)=TUSED(12,NM)+SECOND()-TNOW
    !
  End Subroutine INITIALIZE_EVACUATION

  !
  Subroutine INIT_EVAC_GROUPS
    Implicit None
    !
    ! Initialize group lists, known doors, etc
    !
    ! Local variables
    Integer I,J, IZERO, nom, j1, ii, i_target_old, i_change_old, i_tmp, i_tmp2
    Integer :: i_egrid, i_target, color_index, i_new_ffield
    Real(EB) :: evel
    ! 
    Type (MESH_TYPE), Pointer :: M=>NULL()
    Type (EVAC_SSTAND_Type), Pointer :: ESS=>NULL()
    Type (HUMAN_TYPE), Pointer :: HR=>NULL()
    !
    If (.Not.Any(EVACUATION_GRID)) Return
    If (PROCESS_STOP_STATUS > 0) Return

    !
    ilh_dim = ilh           ! lonely humans dimension
    i33_dim = i33
    Allocate(Group_List(0:i33_dim),STAT=IZERO)
    Call ChkMemErr('Initialize_Evacuation','Group_List',IZERO) 
    Group_List(:)%Tdoor   = Huge(Group_List(0)%Tdoor)
    Group_List(:)%COMPLETE = 0
    !
    Do i = 0,i33_dim
       Allocate(Group_List(i)%GROUP_I_FFIELDS(n_egrids),STAT=izero)
       Group_List(i)%GROUP_I_FFIELDS(:) = 0
    End Do

    Allocate(Group_Known_Doors(1:i33_dim),STAT=IZERO)
    Call ChkMemErr('Initialize_Evacuation', 'Group_Known_Doors',IZERO) 
    Allocate(Human_Known_Doors(1:ilh_dim),STAT=IZERO)
    Call ChkMemErr('Initialize_Evacuation', 'Human_Known_Doors',IZERO) 
    !
    ! These arrays are common for the whole module EVAC.  They are only allocated once.
    Allocate(Is_Known_Door(Max(1,n_doors+n_exits)),STAT=IZERO)
    Call ChkMemErr('Initialize_Evacuation','Is_Known_Door',IZERO) 
    Allocate(Is_Visible_Door(Max(1,n_doors+n_exits)),STAT=IZERO)
    Call ChkMemErr('Initialize_Evacuation','Is_Visible_Door',IZERO)
    Allocate(FED_max_Door(Max(1,n_doors+n_exits)),STAT=IZERO)
    Call ChkMemErr('Initialize_Evacuation','FED_max_Door',IZERO) 
    Allocate(K_ave_Door(Max(1,n_doors+n_exits)),STAT=IZERO)
    Call ChkMemErr('Initialize_Evacuation','K_ave_Door',IZERO) 
    Allocate(Color_Tmp(Max(1,i33_dim)),STAT=IZERO)
    Call ChkMemErr('Initialize_Evacuation','Color_Tmp',IZERO) 

    i_egrid = 0
    If (n_doors >0) EVAC_DOORS(:)%R_NTARGET = 5.0_EB
    If (n_exits >0) EVAC_EXITS(:)%R_NTARGET = 5.0_EB
    Do nom = 1, NMESHES
       n_change_doors  = 0 ! Count the initialization Nash Equilibrium iterations
       n_change_trials = 0 ! Count the initialization Nash Equilibrium iterations
       i_change_old    = -1
       If ( .Not.(EVACUATION_ONLY(NoM) .And. EVACUATION_GRID(NoM)) ) Cycle
       M => MESHES(NOM)
       Group_List(:)%GROUP_SIZE  = 0
       Group_List(:)%GROUP_X = 0.0_EB
       Group_List(:)%GROUP_Y = 0.0_EB
       Group_List(:)%Speed   = 0.0_EB
       i_egrid = i_egrid+1
       Do i = 1, M%N_HUMANS
          HR => M%HUMAN(I)
          j = Max(0,HR%GROUP_ID)
          Group_List(j)%GROUP_SIZE = Group_List(j)%GROUP_SIZE + 1
          Group_List(j)%GROUP_X    = Group_List(j)%GROUP_X + HR%X
          Group_List(j)%GROUP_Y    = Group_List(j)%GROUP_Y + HR%Y
          Group_List(j)%Speed      = Group_List(j)%Speed + HR%Speed
          Do ii = 1,n_doors
             If (nom .Ne. EVAC_DOORS(ii)%Imesh) Cycle
             evel = Sqrt((HR%X-EVAC_DOORS(ii)%X)**2 + (HR%Y-EVAC_DOORS(ii)%Y)**2)
             EVAC_DOORS(ii)%R_NTARGET = Max(evel,EVAC_DOORS(ii)%R_NTARGET)
          End Do
          Do ii = 1,n_exits
             If (nom .Ne. EVAC_EXITS(ii)%Imesh) Cycle
             evel = Sqrt((HR%X-EVAC_EXITS(ii)%X)**2 + (HR%Y-EVAC_EXITS(ii)%Y)**2)
             EVAC_EXITS(ii)%R_NTARGET = Max(evel,EVAC_EXITS(ii)%R_NTARGET)
          End Do
       End Do
       Group_List(1:)%GROUP_X = Group_List(1:)%GROUP_X / Max(1,Group_List(1:)%GROUP_SIZE)
       Group_List(1:)%GROUP_Y = Group_List(1:)%GROUP_Y / Max(1,Group_List(1:)%GROUP_SIZE)
       Group_List(1:)%Speed   = Group_List(1:)%Speed / Max(1,Group_List(1:)%GROUP_SIZE)
       i_tmp = 0
       i_tmp2 = -1
       If (M%N_HUMANS < 1) Cycle
       Do While (.Not.(i_change_old == n_change_doors) .And. .Not.(i_tmp2 == i_tmp)) ! Iterate Nash equilibrium
          i_change_old = n_change_doors
          i_tmp = 0
          Do i = 1, M%N_HUMANS
             HR => M%HUMAN(I)
             ! group_id > 0: +group_id
             ! group_id < 0: -human_id (lonely humans)
             j  =  Max(0,HR%GROUP_ID)
             j1 = -Min(0,HR%GROUP_ID)
             ! Test, if this group has already a ffield (on this floor)
             ! Lonely humans have j=0 and group_i_ffields is always 0.
             If (Group_List(j)%GROUP_I_FFIELDS(i_egrid) == 0) Then
                i_target_old = HR%I_Target
                n_change_trials = n_change_trials + 1
                Call Change_Target_Door(nom, nom, i, j, j1, i_egrid, 0, HR%X, HR%Y, i_target, color_index, i_new_ffield, HR)
                If (Abs(i_target_old) .Ne. Abs(i_target)) Then
                   n_change_doors = n_change_doors + 1
                   i_tmp = i
                   If (i_target > 0) Then
                      If (i_target <= n_doors ) Then
                         evel = Sqrt((HR%X-EVAC_DOORS(i_target)%X)**2 + (HR%Y-EVAC_DOORS(i_target)%Y)**2)
                         evel = 50.0_EB*evel/EVAC_DOORS(i_target)%R_NTARGET + 1.0_EB
                         ii = Min(50,Max(1,Int(evel)))
                         EVAC_DOORS(i_target)%NTARGET(ii:50) = EVAC_DOORS(i_target)%NTARGET(ii:50) + 1
                      Else
                         evel = Sqrt((HR%X-EVAC_EXITS(i_target-n_doors)%X)**2 + (HR%Y-EVAC_EXITS(i_target-n_doors)%Y)**2)
                         evel = 50.0_EB*evel/EVAC_EXITS(i_target-n_doors)%R_NTARGET + 1.0_EB
                         ii = Min(50,Max(1,Int(evel)))
                         EVAC_EXITS(i_target-n_doors)%NTARGET(ii:50) = EVAC_EXITS(i_target-n_doors)%NTARGET(ii:50) + 1
                      End If
                   Else  ! i_target < 0 means non-visible target door
                      n_change_doors = n_change_doors - 1  ! do not iterate non-visible doors
                   End If
                   If (i_target_old > 0) Then
                      If (i_target_old <= n_doors ) Then
                         evel = Sqrt((HR%X-EVAC_DOORS(i_target_old)%X)**2 + (HR%Y-EVAC_DOORS(i_target_old)%Y)**2)
                         evel = 50.0_EB*evel/EVAC_DOORS(i_target_old)%R_NTARGET + 1.0_EB
                         ii = Min(50,Max(1,Int(evel)))
                         EVAC_DOORS(i_target_old)%NTARGET(ii:50) = EVAC_DOORS(i_target_old)%NTARGET(ii:50) - 1
                      Else
                         evel = Sqrt((HR%X-EVAC_EXITS(i_target_old-n_doors)%X)**2 + &
                              (HR%Y-EVAC_EXITS(i_target_old-n_doors)%Y)**2)
                         evel = 50.0_EB*evel/EVAC_EXITS(i_target_old-n_doors)%R_NTARGET + 1.0_EB
                         ii = Min(50,Max(1,Int(evel)))
                         EVAC_EXITS(i_target_old-n_doors)%NTARGET(ii:50) = &
                              EVAC_EXITS(i_target_old-n_doors)%NTARGET(ii:50) - 1
                      End If
                   End If
                End If
             Else              ! this group has already a flow field
                ! This group has already tried to change the field
                If (COLOR_METHOD == 5 .And. j > 0) HR%COLOR_INDEX = Color_Tmp(j)
                If (COLOR_METHOD == 4 .And. j > 0) HR%COLOR_INDEX = Color_Tmp(j)
                HR%I_FFIELD    = Group_List(j)%GROUP_I_FFIELDS(i_egrid)
                HR%FFIELD_NAME = Trim(MESH_NAME(HR%I_FFIELD))
                HR%I_Target = Group_Known_Doors(j)%I_Target
             End If            ! first member of a group or a lonely soul
          End Do              ! 1, n_humans
          If (n_change_doors-i_change_old == 1) Then
             Write(LU_EVACOUT,fmt='(a,2i10)') ' Init: Door Changes i_tmp ',  i_tmp, i_tmp2
             i_tmp2 = i_tmp
          Else
             i_tmp2 = -1
          End If
          If (n_change_doors/Max(1,M%N_HUMANS) > 10*M%N_HUMANS) i_tmp2 = i_tmp
          If (FAC_DOOR_QUEUE <= 0.001_EB) i_change_old = n_change_doors  ! Do not iterate the Nash equilibrium
       End Do         ! Nash iterations
       If (FAC_DOOR_QUEUE > 0.001_EB) Write(LU_EVACOUT,fmt='(a,f14.2,a,i8)') &
            ' Init: Changes per agent ', Real(n_change_doors,EB)/Real(M%N_HUMANS,EB), &
            ', Nash iterations', n_change_trials/M%N_HUMANS
    End Do                  ! 1, nmeshes

    Write (LU_EVACOUT,fmt='(/a)') ' EVAC: Initial positions of the humans'
    Write (LU_EVACOUT,fmt='(a,a)') ' person     x       y       z    Tpre    Tdet  ', &
         ' dia    v0   tau   i_gr i_ff'

    ! Initialize the group_i_ffields
    i_egrid = 0
    Do nom = 1, NMESHES
       If ( .Not.(EVACUATION_ONLY(NoM) .And. EVACUATION_GRID(NoM)) ) Cycle
       M => MESHES(NOM)
       i_egrid = i_egrid+1
       Do i = 1, M%N_HUMANS
          HR => M%HUMAN(I)

          ! Check spectator stands, correct the z-coordiante
          SS_Loop: Do j = 1, n_sstands
             ESS => EVAC_SSTANDS(j)
             If (ESS%IMESH == nom .And. (ESS%X1 <= HR%X .And. ESS%X2 >= HR%X) .And. &
                  (ESS%Y1 <= HR%Y .And. ESS%Y2 >= HR%Y) ) Then
                Select Case (ESS%IOR)
                Case(-1)
                   HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + (ESS%H-ESS%H0)*Abs(ESS%X1-HR%X)/Abs(ESS%X1-ESS%X2)
                Case(+1)
                   HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + (ESS%H-ESS%H0)*Abs(ESS%X2-HR%X)/Abs(ESS%X1-ESS%X2)
                Case(-2)
                   HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + (ESS%H-ESS%H0)*Abs(ESS%Y1-HR%Y)/Abs(ESS%Y1-ESS%Y2)
                Case(+2)
                   HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + (ESS%H-ESS%H0)*Abs(ESS%Y2-HR%Y)/Abs(ESS%Y1-ESS%Y2)
                End Select
                Exit SS_Loop
             End If
          End Do SS_Loop

          j = Max(0,HR%GROUP_ID)
          Group_List(j)%IEL = HR%IEL

          Write (LU_EVACOUT,fmt='(i6,5f8.2,3f6.2,i6,i4,i4)') HR%ILABEL, &
               HR%X, HR%Y, HR%Z, HR%Tpre, HR%Tdet,2.0_EB*HR%Radius, &
               HR%Speed, HR%Tau, HR%GROUP_ID, HR%i_ffield, HR%COLOR_INDEX
       End Do
    End Do
    Write (LU_EVACOUT,fmt='(/)')
    !n_change_doors  = 0 ! Do not count the initialization Nash Equilibrium iterations
    !n_change_trials = 0 ! Do not count the initialization Nash Equilibrium iterations

  End Subroutine INIT_EVAC_GROUPS
!
  Subroutine EVAC_MESH_EXCHANGE(T,T_Save,I_mode, ICYC, EXCHANGE_EVACUATION, MODE)
    Implicit None
    !
    ! Passed variables
    Real(EB), Intent(IN) :: T
    Integer, Intent(IN) :: I_mode, ICYC, MODE
    Logical, Intent(OUT) :: EXCHANGE_EVACUATION
    Real(EB), Intent(INOUT) :: T_Save
    !
    ! Local variables
    Integer :: nm, nom, i, j, i1, j1, k1
    Integer :: ios, IZERO
    Logical L_use_fed, L_fed_read, L_fed_save
    Real(EB) DT_Save
    Integer(4) ibar_tmp, jbar_tmp, kbar_tmp, n_tmp
    Real(FB) tmpout1, tmpout2, tmpout3, tmpout4, t_tmp, dt_tmp
    Real(FB) tmpout5, tmpout6, tmpout7, tmpout8
    REAL(EB), Allocatable, Dimension(:) :: YY_GET
    !
    If (.Not. Any(EVACUATION_GRID)) Return
    If (ICYC < 1) Return
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

    L_use_fed  = .False.
    L_fed_read = Btest(I_mode,3)
    L_fed_save = Btest(I_mode,1)

    !
    ! Change information: Fire meshes ==> Evac meshes
    If ( T >= T_BEGIN .And. Real(T,FB) >= Real(T_Save,FB) ) Then
       If (MODE > 0) T_Save = T + DT_Save
       L_use_fed = .True.
    End If

    L_use_fed = L_use_fed .And. (L_fed_read .Or. L_fed_save)

    If (MODE < 1 .And. L_use_fed) EXCHANGE_EVACUATION = .True.
    If (MODE < 2) Return

    If (L_use_fed) Then
       If (L_fed_save) Then
          Allocate(YY_GET(1:Max(1,N_SPECIES)),STAT=IZERO)
          Call ChkMemErr('EVAC_MESH_EXCHANGE', 'YY_GET',IZERO) 
          Write (LU_EVACFED) Real(T,FB), Real(DT_Save,FB)
       Else
          Read (LU_EVACFED,End=324,Iostat=ios) t_tmp, dt_tmp
          T_Save = t_tmp + dt_tmp ! next time point in the file
       End If

       ! Next loop interpolates fire mesh (soot+fed) into human_grids and
       ! saves it to the disk. Or it reads FED+soot from the disk.
       MESH_LOOP: Do NM=1,NMESHES
          If ( .Not.(EVACUATION_GRID(NM) .And. EVACUATION_ONLY(NM)) ) Cycle
          !
          Call POINT_TO_MESH(NM)
          !             Call POINT_TO_EVAC_MESH(NM)
          If (L_fed_save) Then
             ibar_tmp = IBAR
             jbar_tmp = JBAR
             kbar_tmp = 1
             n_tmp    = 4  ! New format (version 1.11)
             Write (LU_EVACFED) ibar_tmp, jbar_tmp, kbar_tmp, n_tmp
          Else
             Read (LU_EVACFED,Iostat=ios) ibar_tmp, jbar_tmp, kbar_tmp, n_tmp
             If (ios.Ne.0) Then
                Write(MESSAGE,'(A)') 'ERROR: Evac Mesh Exchange: FED READ ERROR'
                Close (LU_EVACFED)
                Call SHUTDOWN(MESSAGE)
             End If
             If (ibar_tmp /= IBAR .Or. jbar_tmp /= JBAR .Or. n_tmp < 4 ) Then
                Close (LU_EVACFED)
                Call SHUTDOWN('ERROR: Problems to read the FED file')
             End If

          End If

          Do i = 1, IBAR
             Do j= 1, JBAR

                If (L_fed_save) Then
                   If ( Abs(HUMAN_GRID(i,j)%IMESH) > 0 ) Then
                      ! imesh > 0, i.e. fire grid found
                      i1 = HUMAN_GRID(i,j)%II
                      j1 = HUMAN_GRID(i,j)%JJ
                      k1 = HUMAN_GRID(i,j)%KK
                      nom = Abs(HUMAN_GRID(i,j)%IMESH)
                      CALL GET_FIRE_CONDITIONS(nom,i1,j1,k1,&
                           HUMAN_GRID(i,j)%FED_CO_CO2_O2,HUMAN_GRID(i,j)%SOOT_DENS,&
                           HUMAN_GRID(i,j)%TMP_G, HUMAN_GRID(i,j)%RADFLUX, YY_GET)
                   End If
                   ! Save Fed, Soot, Temp(C), and Radflux
                   Write (LU_EVACFED) &
                        Real(HUMAN_GRID(i,j)%FED_CO_CO2_O2,FB), &
                        Real(HUMAN_GRID(i,j)%SOOT_DENS,FB), &
                        Real(HUMAN_GRID(i,j)%TMP_G,FB), &
                        Real(HUMAN_GRID(i,j)%RADFLUX,FB)
                Else ! Read FED from a file
                   ! Read Fed, Soot, Temp(C), and Radflux
                   Read (LU_EVACFED,Iostat=ios) tmpout1, tmpout2, tmpout3, tmpout4
                   If (ios.Ne.0) Then
                      Write(MESSAGE,'(A)') 'ERROR: Evac Mesh Exchange: FED READ ERROR'
                      Close (LU_EVACFED)
                      Call SHUTDOWN(MESSAGE)
                   End If
                   HUMAN_GRID(i,j)%FED_CO_CO2_O2 = tmpout1
                   HUMAN_GRID(i,j)%SOOT_DENS = tmpout2
                   HUMAN_GRID(i,j)%TMP_G = tmpout3
                   HUMAN_GRID(i,j)%RADFLUX = tmpout4
                End If   ! Calculate and save FED

             End Do     ! j=1,JBAR
          End Do       ! i=1,IBAR

       End Do MESH_LOOP

       ! Next loop interpolates fire mesh (soot+fed) into human_grids and
       ! saves it to the disk. Or it reads FED+soot from the disk.
       CORR_LOOP: Do i = 1, N_CORRS
          !
          If (L_fed_save) Then

             If ( EVAC_CORRS(i)%FED_MESH > 0 ) Then
                ! Here the fire properties are saved to the arrays.
                i1 = EVAC_CORRS(i)%II(1)
                j1 = EVAC_CORRS(i)%JJ(1)
                k1 = EVAC_CORRS(i)%KK(1)
                nom = EVAC_CORRS(i)%FED_MESH
                CALL GET_FIRE_CONDITIONS(nom,i1,j1,k1,&
                         EVAC_CORRS(i)%FED_CO_CO2_O2(1),EVAC_CORRS(i)%SOOT_DENS(1),&
                         EVAC_CORRS(i)%TMP_G(1), EVAC_CORRS(i)%RADFLUX(1), YY_GET)
             Else
                ! No fed_mesh found
                EVAC_CORRS(i)%FED_CO_CO2_O2(1) = 0.0_EB
                EVAC_CORRS(i)%SOOT_DENS(1) = 0.0_EB
                EVAC_CORRS(i)%TMP_G(1) = 0.0_EB
                EVAC_CORRS(i)%RADFLUX(1) = 0.0_EB
             End If                ! fed_mesh > 0, i.e. fire grid found

             If ( EVAC_CORRS(i)%FED_MESH2 > 0 ) Then
                i1 = EVAC_CORRS(i)%II(2)
                j1 = EVAC_CORRS(i)%JJ(2)
                k1 = EVAC_CORRS(i)%KK(2)
                nom = EVAC_CORRS(i)%FED_MESH2
                CALL GET_FIRE_CONDITIONS(nom,i1,j1,k1,&
                         EVAC_CORRS(i)%FED_CO_CO2_O2(2),EVAC_CORRS(i)%SOOT_DENS(2),&
                         EVAC_CORRS(i)%TMP_G(2), EVAC_CORRS(i)%RADFLUX(2), YY_GET)
             Else
                ! No fed_mesh2 found
                EVAC_CORRS(i)%FED_CO_CO2_O2(2) = 0.0_EB
                EVAC_CORRS(i)%SOOT_DENS(2) = 0.0_EB
                EVAC_CORRS(i)%TMP_G(2) = 0.0_EB
                EVAC_CORRS(i)%RADFLUX(2) = 0.0_EB
             End If                ! fed_mesh2 > 0, i.e. fire grid found

             ! Save Fed, Soot, Temp(C), and Radflux
             Write (LU_EVACFED) &
                  Real(EVAC_CORRS(i)%FED_CO_CO2_O2(1),FB), &
                  Real(EVAC_CORRS(i)%SOOT_DENS(1),FB), &
                  Real(EVAC_CORRS(i)%TMP_G(1),FB), &
                  Real(EVAC_CORRS(i)%RADFLUX(1),FB), &
                  Real(EVAC_CORRS(i)%FED_CO_CO2_O2(2),FB), &
                  Real(EVAC_CORRS(i)%SOOT_DENS(2),FB), &
                  Real(EVAC_CORRS(i)%TMP_G(2),FB), &
                  Real(EVAC_CORRS(i)%RADFLUX(2),FB)
          Else                    ! Read FED from a file
             ! Read Fed, Soot, Temp(C), and Radflux
             Read (LU_EVACFED,Iostat=ios) tmpout1, tmpout2, tmpout3, tmpout4, tmpout5, tmpout6, tmpout7, tmpout8
             If (ios.Ne.0) Then
                Write(MESSAGE,'(A)') 'ERROR: Evac Mesh Exchange: FED READ ERROR'
                Close (LU_EVACEFF)
                Call SHUTDOWN(MESSAGE)
             End If
             EVAC_CORRS(i)%FED_CO_CO2_O2(1) = tmpout1
             EVAC_CORRS(i)%SOOT_DENS(1) = tmpout2
             EVAC_CORRS(i)%TMP_G(1) = tmpout3
             EVAC_CORRS(i)%RADFLUX(1) = tmpout4
             EVAC_CORRS(i)%FED_CO_CO2_O2(2) = tmpout5
             EVAC_CORRS(i)%SOOT_DENS(2) = tmpout6
             EVAC_CORRS(i)%TMP_G(2) = tmpout7
             EVAC_CORRS(i)%RADFLUX(2) = tmpout8

          End If                  ! Calculate and save FED
       End Do CORR_LOOP

       DOOR_LOOP: Do i = 1, N_DOORS
          !
          If (L_fed_save) Then
             If ( EVAC_DOORS(i)%FED_MESH > 0 ) Then
                i1 = EVAC_DOORS(i)%II
                j1 = EVAC_DOORS(i)%JJ
                k1 = EVAC_DOORS(i)%KK
                nom = EVAC_DOORS(i)%FED_MESH
                CALL GET_FIRE_CONDITIONS(nom,i1,j1,k1,&
                              EVAC_DOORS(i)%FED_CO_CO2_O2,EVAC_DOORS(i)%SOOT_DENS,&
                              EVAC_DOORS(i)%TMP_G, EVAC_DOORS(i)%RADFLUX, YY_GET)
             Else
                ! No fed_mesh found
                EVAC_DOORS(i)%FED_CO_CO2_O2 = 0.0_EB
                EVAC_DOORS(i)%SOOT_DENS = 0.0_EB
                EVAC_DOORS(i)%TMP_G = 0.0_EB
                EVAC_DOORS(i)%RADFLUX = 0.0_EB
             End If                ! fed_mesh > 0, i.e. fire grid found

             ! Save Fed, Soot, Temp(C), and Radflux
             Write (LU_EVACFED) &
                  Real(EVAC_DOORS(i)%FED_CO_CO2_O2,FB), &
                  Real(EVAC_DOORS(i)%SOOT_DENS,FB), &
                  Real(EVAC_DOORS(i)%TMP_G,FB), &
                  Real(EVAC_DOORS(i)%RADFLUX,FB)
          Else                    ! Read FED from a file
             ! Read Fed, Soot, Temp(C), and Radflux
             Read (LU_EVACFED,Iostat=ios) tmpout1, tmpout2, tmpout3, tmpout4
             If (ios.Ne.0) Then
                Write(MESSAGE,'(A)') 'ERROR: Evac Mesh Exchange: FED READ ERROR'
                Close (LU_EVACFED)
                Call SHUTDOWN(MESSAGE)
             End If
             EVAC_DOORS(i)%FED_CO_CO2_O2 = tmpout1
             EVAC_DOORS(i)%SOOT_DENS = tmpout2
             EVAC_DOORS(i)%TMP_G = tmpout3
             EVAC_DOORS(i)%RADFLUX = tmpout4
          End If                  ! Calculate and save FED
       End Do DOOR_LOOP

       EXIT_LOOP: Do i = 1, N_EXITS
          !
          ! Do not save/read data for counters.
          If (EVAC_EXITS(i)%COUNT_ONLY) Cycle EXIT_LOOP
          If (L_fed_save) Then
             If ( EVAC_EXITS(i)%FED_MESH > 0 ) Then
                i1 = EVAC_EXITS(i)%II
                j1 = EVAC_EXITS(i)%JJ
                k1 = EVAC_EXITS(i)%KK
                nom = EVAC_EXITS(i)%FED_MESH
                CALL GET_FIRE_CONDITIONS(nom,i1,j1,k1,&
                         EVAC_EXITS(i)%FED_CO_CO2_O2,EVAC_EXITS(i)%SOOT_DENS,&
                         EVAC_EXITS(i)%TMP_G, EVAC_EXITS(i)%RADFLUX, YY_GET)
             Else
                ! No fed_mesh found
                EVAC_EXITS(i)%FED_CO_CO2_O2 = 0.0_EB
                EVAC_EXITS(i)%SOOT_DENS = 0.0_EB
                EVAC_EXITS(i)%TMP_G = 0.0_EB
                EVAC_EXITS(i)%RADFLUX = 0.0_EB
             End If                ! fed_mesh > 0, i.e. fire grid found

             ! Save Fed, Soot, Temp(C), and Radflux
             Write (LU_EVACFED) &
                  Real(EVAC_EXITS(i)%FED_CO_CO2_O2,FB), &
                  Real(EVAC_EXITS(i)%SOOT_DENS,FB), &
                  Real(EVAC_EXITS(i)%TMP_G,FB), &
                  Real(EVAC_EXITS(i)%RADFLUX,FB)
          Else                    ! Read FED from a file
             ! Read Fed, Soot, Temp(C), and Radflux
             Read (LU_EVACFED,Iostat=ios) tmpout1, tmpout2, tmpout3, tmpout4
             If (ios.Ne.0) Then
                Write(MESSAGE,'(A)') 'ERROR: Evac Mesh Exchange: FED READ ERROR'
                Close (LU_EVACFED)
                Call SHUTDOWN(MESSAGE)
             End If
             EVAC_EXITS(i)%FED_CO_CO2_O2 = tmpout1
             EVAC_EXITS(i)%SOOT_DENS = tmpout2
             EVAC_EXITS(i)%TMP_G = tmpout3
             EVAC_EXITS(i)%RADFLUX = tmpout4
          End If                  ! Calculate and save FED
       End Do EXIT_LOOP

       If (L_fed_save) Deallocate(YY_GET)

    End If                    ! l_use_fed

324 Continue
    If (ios < 0) Then 
       Write (LU_EVACOUT,fmt='(a,f12.4,a)') 'FED file EOF: time ', T_Save-DT_Save, ' not found'
       Write (LU_EVACOUT,fmt='(a)') 'FED file EOF: use previous values'
       T_Save = 1.0E15
    End If

  End Subroutine EVAC_MESH_EXCHANGE
!
  Subroutine PREPARE_TO_EVACUATE(ICYC)
    Implicit None
    !
    ! Do the mesh independent initializations for the 
    ! subroutine evacuate_humans.
    !
    ! Passed variables
    Integer, Intent(IN) :: ICYC
    !
    ! Local variables
    Logical L_eff_read, L_eff_save
    Integer(4) ibar_tmp, jbar_tmp, kbar_tmp
    Integer nm_tim, i, j
    ! 
    Type (MESH_TYPE), Pointer :: MFF=>NULL()

    If (.Not.Any(EVACUATION_GRID)) Return

    L_eff_read = Btest(I_EVAC,2)
    L_eff_save = Btest(I_EVAC,0)
    If ( ICYC == 0 .And. L_eff_save ) Then
       Do nm_tim = 1, NMESHES
          If (EVACUATION_ONLY(nm_tim)) Then
             MFF=>MESHES(nm_tim)
             ibar_tmp = MFF%IBAR
             jbar_tmp = MFF%JBAR
             kbar_tmp = 1
             Write (LU_EVACEFF) ibar_tmp, jbar_tmp, kbar_tmp
             Do  i = 0, MFF%IBAR+1
                Do j= 0, MFF%JBAR+1
                   Write (LU_EVACEFF) Real(MFF%U(i,j,1),FB), Real(MFF%V(i,j,1),FB)
                End Do
             End Do
          End If
       End Do
       ! Clear the save bit, save is done only once.
       L_eff_save = .False.
       I_EVAC = Ibclr(I_EVAC,0)
    End If
    !
    ! Initialize counters only once for each time step.
    If ( ICYC >= 0 .And. icyc_old < ICYC ) Then
       icyc_old = ICYC
       fed_max_alive = 0.0_EB
       fed_max       = 0.0_EB
    End If

  End Subroutine PREPARE_TO_EVACUATE
!
  Subroutine EVACUATE_HUMANS(Tin,NM,ICYC)
    Implicit None
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
    Real(EB), Intent(IN) :: Tin
    Integer, Intent(IN) :: NM,ICYC
    !
    ! Local variables
    Integer, Parameter :: n_sectors = 2
    Real(EB) DTSP,UBAR,VBAR,X1,Y1,XI,YJ,ZK
    Integer ICN,I,J,IIN,JJN,KKN,II,JJ,KK,IIX,JJY,KKZ,ICX, ICY, N, J1, I_OBST, I_OBSTX, I_OBSTY
    Integer  IE, tim_ic, tim_iw, NM_now, tim_iwx, tim_iwy, tim_iw2, tim_ic2, ibc
    Real(EB) :: P2P_DIST, P2P_DIST_MAX, P2P_U, P2P_V, EVEL, tim_dist
    Real(EB), Dimension(4) :: d_xy
    Logical, Dimension(4) :: FoundWall_xy
    Integer :: istat, STRS_INDX, i_target, color_index, i_new_ffield, i_target_old
    !
    !
    Real(EB) ::  scal_prod_over_rsqr, U_new, V_new, Vmax_timo, CosPhiFac, &
         Speed_max, Delta_min, Dt_sum, C_Yeff, LambdaW, B_Wall, A_Wall, &
         T, Contact_F, Social_F, smoke_beta, smoke_alpha, smoke_speed_fac
    Integer :: iie, jje, iio, jjo, iii, jjj, i_egrid, i_tmp
    Real(EB) :: x_now, y_now, d_humans, d_walls, DTSP_new, &
         fac_tim, dt_group_door, x11, y11, Speed, Tpre
    Logical PP_see_each
    Logical L_eff_read, L_eff_save, L_Dead
    Real(EB) :: cos_x, cos_y, speed_xm, speed_xp, speed_ym, speed_yp, hr_z, hr_a, hr_b, hr_tau, hr_tau_iner
    !
    Real(EB) :: rn, rnCF
    Real(EB) :: GaMe, GaTh, GaCM
    !
    Integer :: i_old_ffield, IZERO, imode_old
    Character(26) :: name_old_ffield
    !
    Real(EB) :: P2P_Torque, Fc_x, Fc_y, Omega_new, angle, A1, Tc_z, Fc_x1, Fc_y1
    Real(EB) :: Omega_max, Omega_0, FAC_V0_UP, FAC_V0_DOWN
    Real(EB), Dimension(6) :: y_tmp, x_tmp, r_tmp, v_tmp, u_tmp
    !
    ! Next are for the block list of the double agent-agent loop (speed up)
    Integer :: Max_Humans_Cell, i_dx, j_dy, ie_max, bl_max
    Real(EB) :: dx_min, dy_min
    Integer, Dimension(:,:,:), Allocatable :: BLOCK_GRID
    Integer, Dimension(:,:), Allocatable :: BLOCK_GRID_N
    Integer, Dimension(:), Allocatable :: BLOCK_LIST
    !
    Real(EB) :: d_humans_min, d_walls_min
    Real(EB) :: TNOW, tnow13, tnow14, tnow15
    !
    Logical :: NM_STRS_MESH
    Real(EB), Dimension(n_sectors+1) :: Sum_suunta, u_theta, v_theta, cos_theta, sin_theta, thetas
    Real(EB) :: Sum_suunta_max, theta_step, theta_start, vr_2r, angle_hre, &
         angle_hr, v_hr, v_hre, P2P_Suunta_MAX, angle_old, hr_x, hr_y, d_new, d_shift, commitment
    Integer :: i_suunta_max, N_suunta_back, N_suunta_backCF
    Integer, Dimension(n_sectors+1) :: N_suunta, N_suuntaCF

    Type (MESH_TYPE), Pointer :: MFF=>NULL()
    Type (EVAC_STRS_Type), Pointer :: STRP=>NULL()
    Type (EVAC_SSTAND_Type), Pointer :: ESS=>NULL()
    Type (HUMAN_TYPE), Pointer :: HR=>NULL(), HRE=>NULL()
    !
    TNOW=SECOND()
    If ( .Not.(EVACUATION_ONLY(NM) .And. EVACUATION_GRID(NM)) ) Return
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
    Call POINT_TO_MESH(NM)
    !
    ! Find the smallest delta_x and/or delta_y.  Note: evac meshes should not be stretched ones.
    dx_min = Minval(DX)
    dy_min = Minval(DY)
    Delta_min = Min(dy_min, dx_min)
    !
    ! Read/Write the EFF-file?
    L_eff_read = Btest(I_EVAC,2)
    L_eff_save = Btest(I_EVAC,0)
    !
    ! Find the egrid index of this main evac mesh. [1, # main evac meshes]
    i_egrid = 0
    Do i = 1, n_egrids
       If (EVAC_Node_List(i)%IMESH == NM) Then
          i_egrid = EVAC_Node_List(i)%Node_Index
       End If
    End Do
    If (i_egrid == 0) Then
       Write(MESSAGE,'(A,I6)') 'ERROR: Evacuate_Humans, no mesh found ',NM
       Call SHUTDOWN(MESSAGE)
    End If
    !
    ! Blocks are used to speed the double loops over the agent-agent interactions.
    Allocate(BLOCK_GRID_N(1:IBAR,1:JBAR),STAT=IZERO)
    Call ChkMemErr('EVACUATE_HUMANS','BLOCK_GRID_N',IZERO)
    ! Initialize some counters etc. for this main evac mesh.
    Do i = 1, n_doors
       If (EVAC_DOORS(i)%IMESH == nm) Then
          If (EVAC_DOORS(i)%TIME_OPEN>EVAC_DOORS(i)%TIME_CLOSE) Then
             imode_old=EVAC_DOORS(i)%IMODE
             If (Tin>EVAC_DOORS(i)%TIME_CLOSE .And. Tin<EVAC_DOORS(i)%TIME_OPEN .And. imode_old==-1) EVAC_DOORS(i)%IMODE=+2
             If (Tin>=EVAC_DOORS(i)%TIME_OPEN .And. imode_old==-2) EVAC_DOORS(i)%IMODE=+1
          Else
             imode_old=EVAC_DOORS(i)%IMODE
             If (Tin>EVAC_DOORS(i)%TIME_OPEN .And. Tin<EVAC_DOORS(i)%TIME_CLOSE .And. imode_old==-2) EVAC_DOORS(i)%IMODE=+1
             If (Tin>=EVAC_DOORS(i)%TIME_CLOSE .And. imode_old==-1) EVAC_DOORS(i)%IMODE=+2
          End If
       End If
    End Do
    Do i = 1, n_exits
       If (EVAC_EXITS(i)%IMESH == nm .And. .Not.EVAC_EXITS(i)%COUNT_ONLY) Then
          If (EVAC_EXITS(i)%TIME_OPEN>EVAC_EXITS(i)%TIME_CLOSE) Then
             imode_old=EVAC_EXITS(i)%IMODE
             If (Tin>EVAC_EXITS(i)%TIME_CLOSE .And. Tin<EVAC_EXITS(i)%TIME_OPEN .And. imode_old==-1) EVAC_EXITS(i)%IMODE=+2
             If (Tin>=EVAC_EXITS(i)%TIME_OPEN .And. imode_old==-2) EVAC_EXITS(i)%IMODE=+1
          Else
             imode_old=EVAC_EXITS(i)%IMODE
             If (Tin>EVAC_EXITS(i)%TIME_OPEN .And. Tin<EVAC_EXITS(i)%TIME_CLOSE .And. imode_old==-2) EVAC_EXITS(i)%IMODE=+1
             If (Tin>=EVAC_EXITS(i)%TIME_CLOSE .And. imode_old==-1) EVAC_EXITS(i)%IMODE=+2
          End If
       End If
    End Do
    ! 
    HUMAN_TIME_LOOP: Do While ( Dt_sum < DT )
       ! DT is the fds flow calculation time step.
       ! Sometimes agent time step is smaller than the fire time step, so
       ! syncronize the agent clock correctly.
       ! Tsteps: Time steps of the agent movement algorithm for different meshes,
       ! which were calculated at the end of the previous main time step.
       ! DTSP: The present time step for the agent movement algorithm, which may be
       ! smaller than the FDS main time step DT.
       DTSP = Min( (DT-Dt_sum), Tsteps(nm) )

       Dt_sum = Dt_sum + DTSP
       T = Tin - DT + Dt_sum     ! Current time for the agents

       Do i = 1, n_doors
          If (EVAC_DOORS(i)%IMESH == nm) Then
             EVAC_DOORS(i)%NTARGET(1:50) = 0
          End If
       End Do
       Do i = 1, n_exits
          If (EVAC_EXITS(i)%IMESH == nm) Then
             EVAC_EXITS(i)%NTARGET(1:50) = 0
          End If
       End Do

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
       Group_List(:)%Tdet            = Huge(Group_List(:)%Tdet)
       !
       Do j = 0, i33_dim
          Group_List(j)%GROUP_I_FFIELDS(i_egrid) = 0
       End Do
       If (n_doors >0) EVAC_DOORS(:)%R_NTARGET = 5.0_EB
       If (n_exits >0) EVAC_EXITS(:)%R_NTARGET = 5.0_EB
       Do i = 1, N_HUMANS
          HR=>HUMAN(I)
          j = Max(0,HR%GROUP_ID)  ! Group index of the agent
          Group_List(j)%GROUP_SIZE = Group_List(j)%GROUP_SIZE + 1
          Group_List(j)%GROUP_X    = Group_List(j)%GROUP_X + HR%X
          Group_List(j)%GROUP_Y    = Group_List(j)%GROUP_Y + HR%Y
          Group_List(j)%Speed      = Group_List(j)%Speed + HR%Speed
          Group_List(j)%IntDose    = Group_List(j)%IntDose + HR%IntDose
          Group_List(j)%Tpre       = Max(Group_List(j)%Tpre,HR%Tpre)
          Group_List(j)%Tdet       = Min(Group_List(j)%Tdet,HR%Tdet)
          i_tmp = HR%I_TARGET
          If (i_tmp > 0 .And. i_tmp <= n_doors ) Then
             evel = Sqrt((HR%X- EVAC_DOORS(i_tmp)%X)**2 + (HR%Y-EVAC_DOORS(i_tmp)%Y)**2)
             EVAC_DOORS(i_tmp)%R_NTARGET = Max(evel,EVAC_DOORS(i_tmp)%R_NTARGET)
          End If
          i_tmp = i_tmp - n_doors
          If (i_tmp > 0 .And. i_tmp <= n_exits ) Then
             evel = Sqrt((HR%X-EVAC_EXITS(i_tmp)%X)**2 + (HR%Y-EVAC_EXITS(i_tmp)%Y)**2)
             EVAC_EXITS(i_tmp)%R_NTARGET = Max(evel,EVAC_EXITS(i_tmp)%R_NTARGET)
          End If
       End Do

       Group_List(1:)%GROUP_X = Group_List(1:)%GROUP_X / Max(1,Group_List(1:)%GROUP_SIZE)
       Group_List(1:)%GROUP_Y = Group_List(1:)%GROUP_Y / Max(1,Group_List(1:)%GROUP_SIZE)
       Group_List(1:)%Speed   = Group_List(1:)%Speed   / Max(1,Group_List(1:)%GROUP_SIZE)
       Group_List(1:)%IntDose = Group_List(1:)%IntDose / Max(1,Group_List(1:)%GROUP_SIZE)
       Group_List(:)%MAX_DIST_CENTER = 0.0_EB

       ! Group_List intitialization and count the number of agents going towards different doors/exits.
       Do i = 1, N_HUMANS
          HR=>HUMAN(I)
          ! group_id > 0: +group_id
          ! group_id < 0: -human_id (lonely agents)
          j  =  Max(0,HR%GROUP_ID)  ! group index > 0
          j1 = -Min(0,HR%GROUP_ID)  ! lonely agent index > 0
          Group_List(j)%MAX_DIST_CENTER = Max(Group_List(j)%MAX_DIST_CENTER, &
               Sqrt((HR%X - Group_List(j)%GROUP_X)**2 + (HR%Y - Group_List(j)%GROUP_Y)**2))
          i_tmp = HR%I_TARGET
          If (i_tmp > 0 .And. i_tmp <= n_doors ) Then
             evel = Sqrt((HR%X- EVAC_DOORS(i_tmp)%X)**2 + (HR%Y-EVAC_DOORS(i_tmp)%Y)**2)
             evel = 50.0_EB*evel/EVAC_DOORS(i_tmp)%R_NTARGET + 1.0_EB
             ie = Min(50,Max(1,Int(evel)))
             EVAC_DOORS(i_tmp)%NTARGET(ie:50) = EVAC_DOORS(i_tmp)%NTARGET(ie:50) + 1
          End If
          i_tmp = i_tmp - n_doors
          If (i_tmp > 0 .And. i_tmp <= n_exits ) Then
             evel = Sqrt((HR%X-EVAC_EXITS(i_tmp)%X)**2 + (HR%Y-EVAC_EXITS(i_tmp)%Y)**2)
             evel = 50.0_EB*evel/EVAC_EXITS(i_tmp)%R_NTARGET + 1.0_EB
             ie = Min(50,Max(1,Int(evel)))
             EVAC_EXITS(i_tmp)%NTARGET(ie:50) = EVAC_EXITS(i_tmp)%NTARGET(ie:50) + 1
          End If
       End Do

       ! j=0, i.e., lonely humans: Group_List is not used , but it is safe to initialize.
       If (N_HUMANS > 0) Then
          Group_List(0)%GROUP_SIZE = 1
          Group_List(0)%COMPLETE   = 1
          Group_List(0)%GROUP_X    = 0.5_EB*(XS+XF)
          Group_List(0)%GROUP_Y    = 0.5_EB*(YS+YF)
          Group_List(0)%Speed      = 1.0_EB
          Group_List(0)%IntDose    = 0.0_EB
          Group_List(0)%MAX_DIST_CENTER = 0.0_EB
       End If

       ! Note: group_size is the number of group members on this main evac mesh.
       ! Check if the groups are already together or not.
       Do j = 1, i33_dim
          Group_List(j)%LIMIT_COMP = RADIUS_COMPLETE_0 + RADIUS_COMPLETE_1*Group_List(j)%GROUP_SIZE
          If ( ((Group_List(j)%MAX_DIST_CENTER <=  Group_List(j)%LIMIT_COMP) .Or. &
               (Group_List(j)%COMPLETE == 1)) .And. Group_List(j)%GROUP_SIZE > 0) Then
             ! If complete=1 already, it stays at 1.
             ! If many floors, check only floors where there are group members.
             ! Group may become complete before fire is detected. Do not count these as complete.
             If (T > Group_List(j)%Tdet) Then
                If (Group_List(j)%COMPLETE == 0) Then
                   ! Tdoor: saves the time point when the group started to move towards the exit.
                   Group_List(j)%Tdoor = Max(T,Group_List(j)%Tdet)
                End If
                Group_List(j)%COMPLETE = 1
             End If
          End If
       End Do

       ! ========================================================
       ! CHANGE TARGET DOOR?
       ! ========================================================
       If (T > 0.0_EB) Then
          Change_Door_Loop: Do ie = 1, N_HUMANS
             HR => HUMAN(ie)
             ! Check and cycle if in stairs
             Check_Strs_Loop0: Do N = 1, N_STRS
                IF (EVAC_STRS(N)%IMESH == HR%IMESH) Cycle Change_Door_Loop
             End Do Check_Strs_Loop0

             i_old_ffield    = HR%I_FFIELD
             name_old_ffield = Trim(MESH_NAME(i_old_ffield))
             j  =  Max(0,HR%GROUP_ID)   ! group index
             j1 = -Min(0,HR%GROUP_ID)   ! lonely agent index

             ! If the group is not yet together, they are not moving towards the exit.
             If (Group_List(j)%COMPLETE == 0) Cycle Change_Door_Loop

             ! Agents start to move towards the exit after the pre-evacuation (reaction) time.
             If (j == 0 .And. T < HR%Tpre+HR%Tdet) Cycle Change_Door_Loop   ! lonely agents
             If (j > 0 .And. T < Group_List(j)%Tpre + Group_List(j)%Tdoor) Cycle Change_Door_Loop ! groups

             If (Group_List(j)%GROUP_I_FFIELDS(i_egrid) == 0) Then
                ! Test if this group/agent should update the exit door (exp. distribution)
                Call Random_number(rn)
                If (rn > Exp(-DTSP/dt_group_door) ) Then
                   i_target_old = HR%I_Target 
                   n_change_trials = n_change_trials + 1
                   Call Change_Target_Door(nm, nm, ie, j, j1, i_egrid, 1, HR%X, HR%Y, i_target, color_index, i_new_ffield, HR)
                   If (Abs(i_target_old) .Ne. Abs(i_target)) Then
                      n_change_doors = n_change_doors + 1
                      If (i_target > 0) Then
                         If (i_target <= n_doors ) Then
                            evel = Sqrt((HR%X-EVAC_DOORS(i_target)%X)**2 + (HR%Y-EVAC_DOORS(i_target)%Y)**2)
                            evel = 50.0_EB*evel/EVAC_DOORS(i_target)%R_NTARGET + 1.0_EB
                            ii = Min(50,Max(1,Int(evel)))
                            EVAC_DOORS(i_target)%NTARGET(ii:50) = EVAC_DOORS(i_target)%NTARGET(ii:50) + 1
                         Else
                            evel = Sqrt((HR%X-EVAC_EXITS(i_target-n_doors)%X)**2 + (HR%Y-EVAC_EXITS(i_target-n_doors)%Y)**2)
                            evel = 50.0_EB*evel/EVAC_EXITS(i_target-n_doors)%R_NTARGET + 1.0_EB
                            ii = Min(50,Max(1,Int(evel)))
                            EVAC_EXITS(i_target-n_doors)%NTARGET(ii:50) = EVAC_EXITS(i_target-n_doors)%NTARGET(ii:50) + 1
                         End If
                      End If
                      If (i_target_old > 0) Then
                         If (i_target_old <= n_doors ) Then
                            evel = Sqrt((HR%X-EVAC_DOORS(i_target_old)%X)**2 + (HR%Y-EVAC_DOORS(i_target_old)%Y)**2)
                            evel = 50.0_EB*evel/EVAC_DOORS(i_target_old)%R_NTARGET + 1.0_EB
                            ii = Min(50,Max(1,Int(evel)))
                            EVAC_DOORS(i_target_old)%NTARGET(ii:50) = EVAC_DOORS(i_target_old)%NTARGET(ii:50) - 1
                         Else
                            evel = Sqrt((HR%X-EVAC_EXITS(i_target_old-n_doors)%X)**2 + &
                                 (HR%Y-EVAC_EXITS(i_target_old-n_doors)%Y)**2)
                            evel = 50.0_EB*evel/EVAC_EXITS(i_target_old-n_doors)%R_NTARGET + 1.0_EB
                            ii = Min(50,Max(1,Int(evel)))
                            EVAC_EXITS(i_target_old-n_doors)%NTARGET(ii:50) = &
                                 EVAC_EXITS(i_target_old-n_doors)%NTARGET(ii:50) - 1
                         End If
                      End If
                   End If
                Else ! Do not update door flow field at this time step
                   If (j > 0) Group_Known_Doors(j)%I_Target = HR%I_Target 
                   If (COLOR_METHOD == 5 .And. j > 0) Color_Tmp(j) = HR%COLOR_INDEX
                   If (COLOR_METHOD == 4 .And. j > 0) Color_Tmp(j) = HR%COLOR_INDEX
                   If (j > 0) Group_List(j)%GROUP_I_FFIELDS(i_egrid) = HR%I_FFIELD 
                End If  ! update door flow field: is rn large enough?
             Else
                ! This group has already tried to change the field.
                If (COLOR_METHOD == 5 .And. j > 0) HR%COLOR_INDEX = Color_Tmp(j)
                If (COLOR_METHOD == 4 .And. j > 0) HR%COLOR_INDEX = Color_Tmp(j)
                HR%I_FFIELD    = Group_List(j)%GROUP_I_FFIELDS(i_egrid)
                HR%FFIELD_NAME = Trim(MESH_NAME(HR%I_FFIELD))
                If (j > 0) HR%I_Target = Group_Known_Doors(j)%I_Target
             End If  ! group_i_field=0
          End Do Change_Door_Loop  ! loop over agents 
       End If  ! t > 0

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
       d_humans_min = Huge(d_humans_min)
       d_walls_min  = Huge(d_walls_min)

       EVAC_MOVE_LOOP: Do I=1, N_HUMANS
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
          EVEL = Sqrt(HR%U**2 + HR%V**2)
          L_Dead  = .False.
          If ( HR%IntDose >= 1.0_EB  ) Then
             L_Dead = .True.
             ! No random force for a dead person.
             GaTh = 0.0_EB
             ! No psychological force terms for a dead person.
             A_Wall = 0.0_EB
             If (HR%Tpre /= Huge(HR%Tpre)) Then
                n_dead = n_dead+1
                Write (LU_EVACOUT,fmt='(a,i6,a,f8.2,a,i6)') ' EVAC: Person n:o', &
                     HR%ILABEL, ' dead at ', T, ' s, number of casualties ', n_dead
             End If
             HR%Tpre = Huge(HR%Tpre)
             HR%Tdet = Huge(HR%Tdet)
             HR%COLOR_INDEX = EVAC_AVATAR_NCOLOR
          Else
             fed_max_alive = Max(fed_max_alive,HR%IntDose)
          End If
          fed_max = Max(fed_max,HR%IntDose)  ! dead or alive
          hr_tau      = HR%Tau
          hr_tau_iner = HR%Tau_Iner
          ! Counterflow: increase motivation to go ahead
          If (HR%Commitment > 0.01_EB) Then
             hr_tau      = Max(CF_MIN_TAU, &
                  HR%Commitment*CF_FAC_TAUS*hr_tau + (1.0_EB-HR%Commitment)*hr_tau)
             hr_tau_iner = Max(CF_MIN_TAU_INER, &
                  HR%Commitment*CF_FAC_TAUS*hr_tau_iner + (1.0_EB-HR%Commitment)*hr_tau_iner)
          End If
          !
          ! In which grid cell is the agent?
          !
          XI = CELLSI(Floor((HR%X-XS)*RDXINT))
          YJ = CELLSJ(Floor((HR%Y-YS)*RDYINT))
          ZK = CELLSK(Floor((HR_Z-ZS)*RDZINT))
          II = Floor(XI+1.0_EB)
          JJ = Floor(YJ+1.0_EB)
          KK = 1
          IIX = Floor(XI+0.5_EB)
          JJY = Floor(YJ+0.5_EB)
          KKZ = Floor(ZK+0.5_EB)
          HR%W = 0.0_EB

          ! Check the smoke density for the detection by smoke
          If (.Not.L_Dead .And. HUMAN_GRID(ii,jj)%SOOT_DENS > TDET_SMOKE_DENS) Then
             HR%Tdet = Min(HR%Tdet,T)
          End If

          smoke_speed_fac = 1.0_EB
          If (T > T_BEGIN) Then
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
             If (MASS_EXTINCTION_COEFFICIENT*1.0E-6_EB*HUMAN_GRID(II,JJ)%SOOT_DENS > 3.0_EB/SMOKE_MIN_SPEED_VISIBILITY) Then
                smoke_speed_fac = Min(smoke_speed_fac, smoke_speed_fac*( 2.0_EB - &
                     ( MASS_EXTINCTION_COEFFICIENT*1.0E-6_EB*HUMAN_GRID(II,JJ)%SOOT_DENS - & 
                     3.0_EB/SMOKE_MIN_SPEED_VISIBILITY ) / (3.0_EB/SMOKE_MIN_SPEED_VISIBILITY) ) )
             End If
             smoke_speed_fac = Max(smoke_speed_fac, SMOKE_MIN_SPEED)
          End If
          HR%v0_fac = smoke_speed_fac

          ! ========================================================
          ! Calculate persons prefered walking direction v0
          ! ========================================================
          NM_STRS_MESH = .False.
          StrsMeshLoop: Do N = 1, N_STRS
             If (EVAC_STRS(N)%IMESH==NM_now) Then     
                NM_STRS_MESH = .True.
                Exit StrsMeshLoop
             End If
          End Do StrsMeshLoop
          If (TAU_CHANGE_V0 > 1.0E-12_EB) Then
             ! Collision avoidance (incl. counterflow), do not update v0 on every time step.
             UBAR = HR%UBAR; VBAR = HR%VBAR
          Else
             Call Find_Prefered_Direction(I, N, T, T_BEGIN, L_Dead, NM_STRS_MESH, &
                  II, JJ, IIX, JJY, XI, YJ, ZK, UBAR, VBAR, hr_tau, Tpre)
          End If
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
          If (NM_STRS_MESH) Then
             STRS_Indx = N
             STRP=>EVAC_STRS(N)     
          End If
          !
          ! Check if an agent is on a spectator stand.
          SS_Loop1: Do j = 1, n_sstands
             ESS => EVAC_SSTANDS(j)
             If (ESS%IMESH /= NM) Cycle SS_Loop1
             If (ESS%X1 > HR%X) Cycle SS_Loop1
             If (ESS%X2 < HR%X) Cycle SS_Loop1
             If (ESS%Y1 > HR%Y) Cycle SS_Loop1
             If (ESS%Y2 < HR%Y) Cycle SS_Loop1
             cos_x = ESS%cos_x
             cos_y = ESS%cos_y
             FAC_V0_UP   = ESS%FAC_V0_UP
             FAC_V0_DOWN = ESS%FAC_V0_DOWN
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
             Select Case (ESS%IOR)
             Case(-1)
                speed_xm = cos_x*HR%Speed*FAC_V0_DOWN
                speed_xp = cos_x*HR%Speed*FAC_V0_UP
                speed_ym = HR%Speed*ESS%FAC_V0_HORI
                speed_yp = HR%Speed*ESS%FAC_V0_HORI
                HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + (ESS%H-ESS%H0)*Abs(ESS%X1-HR%X)/Abs(ESS%X1-ESS%X2)
             Case(+1)
                speed_xm = cos_x*HR%Speed*FAC_V0_UP
                speed_xp = cos_x*HR%Speed*FAC_V0_DOWN
                speed_ym = HR%Speed*ESS%FAC_V0_HORI
                speed_yp = HR%Speed*ESS%FAC_V0_HORI
                HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + (ESS%H-ESS%H0)*Abs(ESS%X2-HR%X)/Abs(ESS%X1-ESS%X2)
             Case(-2)
                speed_xm = HR%Speed*ESS%FAC_V0_HORI
                speed_xp = HR%Speed*ESS%FAC_V0_HORI
                speed_ym = cos_y*HR%Speed*FAC_V0_DOWN
                speed_yp = cos_y*HR%Speed*FAC_V0_UP
                HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + (ESS%H-ESS%H0)*Abs(ESS%Y1-HR%Y)/Abs(ESS%Y1-ESS%Y2)
             Case(+2)
                speed_xm = HR%Speed*ESS%FAC_V0_HORI
                speed_xp = HR%Speed*ESS%FAC_V0_HORI
                speed_ym = cos_y*HR%Speed*FAC_V0_UP
                speed_yp = cos_y*HR%Speed*FAC_V0_DOWN
                HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + (ESS%H-ESS%H0)*Abs(ESS%Y2-HR%Y)/Abs(ESS%Y1-ESS%Y2)
             End Select
             Exit SS_Loop1
             HR%Z = 0.5_EB*(ZS+ZF)  ! The agent is not on any incline
          End Do SS_Loop1

          ! Set height and speed for an agent in stairs
          If (NM_STRS_MESH) Then 
             Call GetStairSpeedAndZ(speed_xm, speed_xp, speed_ym, speed_yp, STRP, HR)
          End If

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
          EVEL = Sqrt(UBAR**2 + VBAR**2)
          If (EVEL > 0.0_EB) Then
             !Inclines: U,V are the (x,y) plane projection velocities
             speed = speed_xp*(0.5_EB + Sign(0.5_EB,UBAR)) + speed_xm*(0.5_EB - Sign(0.5_EB,UBAR)) 
             speed = speed*HR%v0_fac
             U_new = U_new + 0.5_EB*(DTSP/hr_tau)*(speed*(UBAR/EVEL) - HR%U)
             HR%UBAR = UBAR

             speed = speed_yp*(0.5_EB + Sign(0.5_EB,VBAR)) + speed_ym*(0.5_EB - Sign(0.5_EB,VBAR)) 
             speed = speed*HR%v0_fac
             V_new = V_new + 0.5_EB*(DTSP/hr_tau)*(speed*(VBAR/EVEL) - HR%V)
             HR%VBAR = VBAR

             If (VBAR >= 0.0_EB) Then
                angle = Acos(UBAR/EVEL)
             Else
                angle = 2.0_EB*Pi - Acos(UBAR/EVEL)
             End If

             ! Collison avoidance
             angle = angle + HR%angle_old
             Do While (angle >= 2.0_EB*Pi)
                angle = angle - 2.0_EB*Pi
             End Do
             Do While (angle < 0.0_EB)
                angle = angle + 2.0_EB*Pi
             End Do

             If (angle == 2.0_EB*Pi) angle = 0.0_EB  ! angle is [0,2Pi)

             ! Rotational motion: J(dw/dt) = (J/t_iner)*( ((angle-angle_0/pi))w_0 - w )
             If (Abs(angle-HR%angle) <= Pi ) Then
                ! zero is not crossed.
                Omega_new = Omega_new + 0.5_EB*(DTSP/HR_tau_iner)*( (angle-HR%angle)*(Omega_0/Pi) - HR%Omega)
             Else
                ! zero is crossed
                Omega_new = Omega_new + 0.5_EB*(DTSP/HR_tau_iner)* &
                     ( (2.0_EB*Pi-Abs(angle-HR%angle))*Sign(1.0_EB , HR%angle-angle)*(Omega_0/Pi) - HR%Omega)
             End If
          Else  ! No v0
             U_new = U_new + 0.5_EB*(DTSP/hr_tau)* (- HR%U)
             V_new = V_new + 0.5_EB*(DTSP/hr_tau)* (- HR%V)
             ! Slow rotation down if no direction available, i.e., target Omega_0 is zero.
             Omega_new = Omega_new + 0.5_EB*(DTSP/HR_tau_iner)*(-HR%Omega)
             HR%UBAR = 0.0_EB
             HR%VBAR = 0.0_EB
          End If
          !
          ! Check, if a new random force is needed on next time step.
          ! Poisson distribution, i.e., the agents do not have memory.
          ! P[No change during dt] = Exp(-DTSP/HR%Tau)
          If ( GaTh > 0.0_EB .And. T > T_BEGIN ) Then
             Call Random_number(rn)
             If ( rn > Exp(-DTSP/(0.2_EB*hr_tau)) ) HR%NewRnd = .True.
             If ( HR%NewRnd ) Then
                HR%NewRnd = .False.
                GaTh = GaTh*(HR%Speed/1.3_EB)**2
                ! GaTh is found to be more or less nice for speeds about 1.3 m/s
                HR%ksi = (GaussRand(GaMe, GaTh, GaCM))
                Call Random_number(rn)
                HR%eta = 2.0_EB*Pi*rn
                EVEL = Sqrt(HR%U**2 + HR%V**2)
                ! Random noice variance: GaTh = Max( GaTh, GaTh*(100.0_EB-110.0_EB*Abs(EVEL/HR%Speed)) )
                ! Above GaTh is found to be more or less nice for speeds about 1.3 m/s
                ! Scale the random force by the current target speed v0 (note: Sqrt(GaTh) = std.dev.)
                ! Note: Speed reduction due to smoke is treated elsewhere.
                HR%ksi = Abs(HR%ksi)
                HR%ksi = Max(HR%ksi, HR%ksi*10.0_EB*(1.0_EB-Min(0.9_EB,EVEL/HR%Speed)))
             End If
          End If

          ! Add random force term
          If ( GaTh > 0.0_EB .And. T > T_BEGIN ) Then
             U_new = U_new + 0.5_EB*DTSP*HR%v0_fac*HR%Mass*HR%ksi*Cos(HR%eta)/HR%Mass
             V_new = V_new + 0.5_EB*DTSP*HR%v0_fac*HR%Mass*HR%ksi*Sin(HR%eta)/HR%Mass
             Omega_new = Omega_new + 0.5_EB*DTSP* 1.0_EB*Sign(HR%ksi,HR%eta-Pi)
          End If

          ! Check that velocities are not too large, i.e., unphysical (less than 10 m/s for humans)
          EVEL = Sqrt(U_new**2 + V_new**2)
          If ( EVEL > Vmax_timo ) Then
             U_new = U_new*(Vmax_timo/EVEL)
             V_new = V_new*(Vmax_timo/EVEL)
          End If
          ! Check that angular velocity is not too large
          If ( Abs(Omega_new) > Omega_max ) Then
             Omega_new = Sign(Omega_max,Omega_new)
          End If

          ! ========================================================
          ! Step 3:  The new coordinates are calculated using the new velocities.
          ! ========================================================
          X1 = HR%X + U_new*DTSP
          Y1 = HR%Y + V_new*DTSP
          HR%U = U_new
          HR%V = V_new
          ! The new body angle, should be on the interval [0,2Pi)
          A1 = HR%Angle + Omega_new*DTSP
          Do While (A1 >= 2.0_EB*Pi)
             A1 = A1 - 2.0_EB*Pi
          End Do
          Do While (A1 < 0.0_EB)
             A1 = A1 + 2.0_EB*Pi
          End Do
          HR%Omega = Omega_new
          !
          ! Check, if human is on an escalator and change the coordinates.
          SS_Loop1b: Do j = 1, n_sstands
             ESS => EVAC_SSTANDS(j)
             If (ESS%IMESH == nm .And. (ESS%X1 <= HR%X .And. ESS%X2 >= HR%X) .And. (ESS%Y1 <= HR%Y .And. ESS%Y2 >= HR%Y) ) Then
                cos_x = ESS%cos_x
                cos_y = ESS%cos_y
                Select Case (ESS%IOR)
                Case(-1)
                   X1 = X1 - cos_x*(ESS%Esc_SpeedUp-ESS%Esc_SpeedDn)*DTSP
                Case(+1)
                   X1 = X1 - cos_x*(ESS%Esc_SpeedUp-ESS%Esc_SpeedDn)*DTSP
                Case(-2)
                   Y1 = Y1 + cos_y*(ESS%Esc_SpeedUp-ESS%Esc_SpeedDn)*DTSP
                Case(+2)
                   Y1 = Y1 - cos_y*(ESS%Esc_SpeedUp-ESS%Esc_SpeedDn)*DTSP
                End Select
                Exit SS_Loop1b
             End If
          End Do SS_Loop1b
          !
          ! In which grid cell is the agent, the new coordinates (t + dt)?
          XI  = CELLSI(Floor((X1-XS)*RDXINT))
          YJ  = CELLSJ(Floor((Y1-YS)*RDYINT))
          ZK  = CELLSK(Floor((HR_Z-ZS)*RDZINT))
          IIN = Floor(XI+1.0_EB)
          JJN = Floor(YJ+1.0_EB)
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
          If ( SOLID(ICN) .AND. .NOT.OBSTRUCTION(I_OBST)%HIDDEN) Then
             If ( (Solid(ICX).AND. .NOT.OBSTRUCTION(I_OBSTX)%HIDDEN) .And. .Not. &
                  (Solid(ICY).AND. .NOT.OBSTRUCTION(I_OBSTY)%HIDDEN) ) Then
                If ( ii < iin ) Then
                   tim_ic = cell_index(iin,jjn,kk)
                   Call Get_iw(iin,jjn,kk,-1,tim_iw)
                   ibc = ijkw(5,tim_iw)
                   If (SURFACE(IBC)%VEL> 0.0_EB .Or. BOUNDARY_TYPE(tim_iw)==OPEN_BOUNDARY) Then
                      HR%X = X1 
                      HR%Y = HR%Y
                   End If
                Else
                   tim_ic = cell_index(iin,jjn,kk)
                   Call Get_iw(iin,jjn,kk,+1,tim_iw)
                   ibc = ijkw(5,tim_iw)
                   If (SURFACE(IBC)%VEL> 0.0_EB .Or. BOUNDARY_TYPE(tim_iw)==OPEN_BOUNDARY) Then
                      HR%X = X1 
                      HR%Y = HR%Y
                   End If
                End If
             Else If ( (Solid(ICY).AND. .NOT.OBSTRUCTION(I_OBSTY)%HIDDEN) .And. .Not. &
                  (Solid(ICX).AND. .NOT.OBSTRUCTION(I_OBSTX)%HIDDEN) ) Then
                If ( jj < jjn ) Then
                   tim_ic = cell_index(iin,jjn,kk)
                   Call Get_iw(iin,jjn,kk,-2,tim_iw)
                   ibc = ijkw(5,tim_iw)
                   If (SURFACE(IBC)%VEL> 0.0_EB .Or. BOUNDARY_TYPE(tim_iw)==OPEN_BOUNDARY) Then
                      HR%X = HR%X
                      HR%Y = Y1
                   End If
                Else
                   tim_ic = cell_index(iin,jjn,kk)
                   Call Get_iw(iin,jjn,kk,+2,tim_iw)
                   ibc = ijkw(5,tim_iw)
                   If (SURFACE(IBC)%VEL> 0.0_EB .Or. BOUNDARY_TYPE(tim_iw)==OPEN_BOUNDARY) Then
                      HR%X = HR%X
                      HR%Y = Y1 
                   End If
                End If
             Else
                Write(MESSAGE,'(A,I4,A,2F8.2)') 'ERROR: Evacuate_Humans, Solid ICX and ICY, mesh ', nm, ' pos ',X1,Y1
             End If
          Else
             ! Target cell is not a solid ==> move
             HR%X = X1
             HR%Y = Y1
          End If
          HR%X = X1
          HR%Y = Y1
          HR%Angle = A1
          !
       End Do EVAC_MOVE_LOOP
       ! ========================================================
       ! MOVE LOOP ENDS HERE
       ! ========================================================

       ! ========================================================
       ! Check if persons are leaving this mesh via doors/exits
       ! and put these persons to the target.
       ! ========================================================
       If (N_HUMANS > 0) Call CHECK_DOORS(T,NM)
       If (N_HUMANS > 0) Call CHECK_EXITS(T,NM)

       ! ========================================================
       ! Check if persons are entering this mesh via a corr.
       ! ========================================================
       Call CHECK_CORRS(T,NM,DTSP)

       ! ========================================================
       ! Add persons from entrys (specified flow rates)
       ! ========================================================
       If (T > T_BEGIN ) Then
          Do i = 1, N_ENTRYS
             Call ENTRY_HUMAN(i,T,NM,istat)
          End Do
       End If

       ! ========================================================
       ! Remove out-of-bounds persons (outside the grid)
       ! ========================================================
       If (N_HUMANS > 0) Call REMOVE_OUT_OF_GRIDS(T,NM)

       If ( ICYC >= 0) Then
          DTSP_new = EVAC_DT_STEADY_STATE
       Else
          DTSP_new = EVAC_DT_FLOWFIELD  ! Initialization phase
       End If
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
       Group_List(:)%Tdet    = Huge(Group_List(:)%Tdet)

       Do j = 0, i33_dim
          Group_List(j)%GROUP_I_FFIELDS(i_egrid) = 0
       End Do
       Do i = 1, N_HUMANS
          HR=>HUMAN(I)
          j = Max(0,HR%GROUP_ID)
          Group_List(j)%GROUP_SIZE = Group_List(j)%GROUP_SIZE + 1
          Group_List(j)%GROUP_X    = Group_List(j)%GROUP_X + HR%X
          Group_List(j)%GROUP_Y    = Group_List(j)%GROUP_Y + HR%Y
          Group_List(j)%GROUP_I_FFIELDS(i_egrid) = HR%i_ffield
          Group_List(j)%Tpre       = Max(Group_List(j)%Tpre,HR%Tpre)
          Group_List(j)%Tdet       = Min(Group_List(j)%Tdet,HR%Tdet)
       End Do
       !
       Group_List(1:)%GROUP_X = Group_List(1:)%GROUP_X / Max(1,Group_List(1:)%GROUP_SIZE)
       Group_List(1:)%GROUP_Y = Group_List(1:)%GROUP_Y / Max(1,Group_List(1:)%GROUP_SIZE)

       ! BLOCK_GRID_N(i,j): How many agents are in this grid cell
       BLOCK_GRID_N = 0
       Group_List(:)%MAX_DIST_CENTER = 0.0_EB
       Do i = 1, N_HUMANS
          HR => HUMAN(i)
          ! Which cell, new coordinates:
          IIN = Floor( CELLSI(Floor((HR%X-XS)*RDXINT)) + 1.0_EB)
          JJN = Floor( CELLSJ(Floor((HR%Y-YS)*RDYINT)) + 1.0_EB)
          BLOCK_GRID_N(IIN,JJN) = BLOCK_GRID_N(IIN,JJN) + 1
          j = Max(0,HR%GROUP_ID)
          Group_List(j)%MAX_DIST_CENTER = Max(Group_List(j)%MAX_DIST_CENTER, &
               Sqrt((HR%X - Group_List(j)%GROUP_X)**2 + (HR%Y - Group_List(j)%GROUP_Y)**2))
       End Do
       Max_Humans_Cell = Max(1,Maxval(BLOCK_GRID_N))

       If (N_HUMANS > 0) Then
          Group_List(0)%GROUP_SIZE = 1
          Group_List(0)%GROUP_X    = 0.5_EB*(XS+XF)
          Group_List(0)%GROUP_Y    = 0.5_EB*(YS+YF)
          Group_List(0)%Speed      = 1.0_EB
          Group_List(0)%IntDose    = 0.0_EB
          Group_List(0)%MAX_DIST_CENTER = 0.0_EB
          Group_List(0)%COMPLETE   = 1
       End If
       ! Check if the groups are already gathered together or not.
       Do j = 1, i33_dim
          Group_List(j)%LIMIT_COMP = RADIUS_COMPLETE_0 + RADIUS_COMPLETE_1*Group_List(j)%GROUP_SIZE
          If ( ((Group_List(j)%MAX_DIST_CENTER <=  Group_List(j)%LIMIT_COMP) .Or. &
               (Group_List(j)%COMPLETE == 1)) .And. Group_List(j)%GROUP_SIZE > 0 ) Then
             If (T > Group_List(j)%Tdet) Then
                If (Group_List(j)%COMPLETE == 0) Then
                   Group_List(j)%Tdoor = Max(T,Group_List(j)%Tdet)
                End If
                Group_List(j)%COMPLETE = 1
             End If
          End If
       End Do

       ! ========================================================
       ! Use blocks to speed up double loops in the force loop.
       ! ========================================================
       Allocate(BLOCK_GRID(1:IBAR,1:JBAR,Max_Humans_Cell),STAT=IZERO)
       Call ChkMemErr('EVACUATE_HUMANS','BLOCK_GRID',IZERO)
       BLOCK_GRID = 0
       ! Block_Grid_N(:,:) how many agents in this cell.
       ! Block_Grid(:,:,1-n) agent indeces

       BLOCK_GRID_N = 0
       Do i = 1, N_HUMANS
          HR => HUMAN(i)
          ! Which grid (block grid) cell, new coordinates:
          IIN = Floor( CELLSI(Floor((HR%X-XS)*RDXINT)) + 1.0_EB )
          JJN = Floor( CELLSJ(Floor((HR%Y-YS)*RDYINT)) + 1.0_EB )
          BLOCK_GRID_N(IIN,JJN) = BLOCK_GRID_N(IIN,JJN) + 1
          BLOCK_GRID(IIN,JJN,BLOCK_GRID_N(IIN,JJN)) = i
       End Do
       i_dx = 2*(Int((2.0_EB*0.3_EB+5.0_EB)/dx_min) + 1) + 1
       j_dy = 2*(Int((2.0_EB*0.3_EB+5.0_EB)/dy_min) + 1) + 1
       i_dx = Max(1,i_dx)
       j_dy = Max(1,j_dy)
       bl_max = i_dx*j_dy*Max_Humans_Cell
       ! Next list will contain the agents that should be looped over for the present agent.
       Allocate(BLOCK_LIST(bl_max),STAT=IZERO)
       Call ChkMemErr('EVACUATE_HUMANS','BLOCK_LIST',IZERO)

       ! ========================================================
       ! Step (3) of SC-VV starts here: Calculate new forces
       ! ========================================================
       TNOW13=SECOND()
       EVAC_FORCE_LOOP: Do I=1,N_HUMANS  
          HR => HUMAN(i)
          j  =  Max(0,HR%GROUP_ID)   ! group index
          j1 = -Min(0,HR%GROUP_ID)   ! lonely agent index
          ! hr%z is the real z-coordinate of the agent (inclines, stairs,etc),
          ! hr_z is z-coordinate of the main evac mesh
          hr_z = 0.5_EB*(ZS+ZF)
          ! Calculate the position and velocities of the shoulder cirles
          y_tmp(1) = HR%Y - Cos(HR%angle)*HR%d_shoulder ! right
          x_tmp(1) = HR%X + Sin(HR%angle)*HR%d_shoulder ! right
          y_tmp(2) = HR%Y ! torso
          x_tmp(2) = HR%X ! torso
          y_tmp(3) = HR%Y + Cos(HR%angle)*HR%d_shoulder ! left
          x_tmp(3) = HR%X - Sin(HR%angle)*HR%d_shoulder ! left
          r_tmp(1) = HR%r_shoulder
          r_tmp(2) = HR%r_torso
          r_tmp(3) = HR%r_shoulder
          u_tmp(1) = HR%U + Cos(HR%angle)*HR%Omega*HR%d_shoulder ! right
          v_tmp(1) = HR%V + Sin(HR%angle)*HR%Omega*HR%d_shoulder ! right
          u_tmp(2) = HR%U ! torso
          v_tmp(2) = HR%V ! torso
          u_tmp(3) = HR%U - Cos(HR%angle)*HR%Omega*HR%d_shoulder ! left
          v_tmp(3) = HR%V - Sin(HR%angle)*HR%Omega*HR%d_shoulder ! left
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
          d_humans = Huge(d_humans)
          d_walls  = Huge(d_walls)
          L_Dead  = .False.
          If (HR%IntDose >= 1.0_EB) Then
             L_Dead = .True.
             ! No random force for a dead person.
             GaTh = 0.0_EB
             ! No psychological force terms for a dead person.
             A_Wall = 0.0_EB
             If (HR%Tpre /= Huge(HR%Tpre)) Then
                n_dead = n_dead+1
                Write (LU_EVACOUT,fmt='(a,i6,a,f8.2,a,i6)') ' EVAC: Person n:o', HR%ILABEL, &
                     ' dead at ', T, ' s, number of casualties ', n_dead
             End If
             HR%Tdet = Huge(HR%Tdet)
             HR%Tpre = Huge(HR%Tpre)
             HR%COLOR_INDEX = EVAC_AVATAR_NCOLOR
          End If
          hr_tau = HR%Tau
          hr_tau_iner = HR%Tau_Iner
          ! =======================================================
          ! Speed dependent social force
          ! =======================================================
          hr_a =  HR%A*Max(0.5_EB,(Sqrt(HR%U**2+HR%V**2)/HR%Speed))
          A_Wall = Min(A_Wall, FAC_A_WALL*hr_a)

          ! Counterflow: increase motivation to go ahead and decrease social force
          If (HR%Commitment > 0.01_EB) Then
             evel = Min(1.0_EB,Sqrt(HR%U**2+HR%V**2)/HR%Speed)
             evel = HR%Commitment*evel + (1.0_EB-HR%Commitment)*1.0_EB
             hr_tau      = Max(CF_MIN_TAU, &
                  HR%Commitment*CF_FAC_TAUS*hr_tau + (1.0_EB-HR%Commitment)*hr_tau)
             hr_tau_iner = Max(CF_MIN_TAU_INER, &
                  HR%Commitment*CF_FAC_TAUS*hr_tau_iner + (1.0_EB-HR%Commitment)*hr_tau_iner)
             hr_a =  HR%A*Max(CF_MIN_A,evel)
             hr_b =  HR%B*Max(CF_MIN_B,evel)
             A_Wall = Min(A_Wall, CF_FAC_A_WALL*FAC_A_WALL*hr_a)
          End If
          !
          ! Psychological force: cut-off when acceleration below 0.0001 m/s**2
          P2P_DIST_MAX = HR%B*Log(HR%A/0.0001_EB)
          P2P_DIST_MAX = Min(P2P_DIST_MAX, 5.0_EB)  ! 5.0 m is the maximum range of PP-force
          If ( HR%SumForces2 > 0.1_EB ) Then
             ! If large pressure then short range forces only (speed up)
             P2P_DIST_MAX = Min( P2P_DIST_MAX, -HR%B*Log(HR%SumForces2/(100.0_EB*HR%A)) )
          End If
          P2P_DIST_MAX = Max(P2P_DIST_MAX, 3.0_EB*HR%B)
          ! Next is the max distance for the collision avoidance, counterflow, etc.
          P2P_Suunta_MAX = Max(P2P_DIST_MAX, 3.0_EB)

          ! Speed up the dead agent loop, only contact forces are needed.
          If (L_Dead) P2P_DIST_MAX = 0.0_EB
          If (L_Dead) P2P_Suunta_MAX = 0.0_EB
          ! Check if counterflow algorithm is not used.
          If (TAU_CHANGE_V0 < 1.E-12_EB) P2P_Suunta_MAX = P2P_DIST_MAX

          ! In which grid cell is the agent, new coordinates:
          XI = CELLSI(Floor((HR%X-XS)*RDXINT))
          YJ = CELLSJ(Floor((HR%Y-YS)*RDYINT))
          ZK = CELLSK(Floor((HR_Z-ZS)*RDZINT))
          IIN  = Floor(XI+1.0_EB)
          JJN  = Floor(YJ+1.0_EB)
          KKN = 1
          IIX = Floor(XI+0.5_EB)
          JJY = Floor(YJ+0.5_EB)
          KKZ = 1
          ICN = CELL_INDEX(IIN,JJN,KKN)
          X1 = HR%X 
          Y1 = HR%Y 
          HR%W = 0.0_EB

          ! ========================================================
          ! Calculate persons prefered walking direction
          ! ========================================================
          Call Find_Prefered_Direction(I, N, T+DTSP_new, T_BEGIN, L_Dead, NM_STRS_MESH, &
               IIN, JJN, IIX, JJY, XI, YJ, ZK, UBAR, VBAR, hr_tau, Tpre)
          ! ========================================================
          ! The prefered walking direction v0 is (UBAR,VBAR)
          ! ========================================================

          ! =======================================================
          ! Update the Block_Grid array search ranges
          ! =======================================================
          BLOCK_LIST = 0
          i_dx = Int((2.0_EB*0.3_EB+Max(P2P_Suunta_MAX,P2P_DIST_MAX))/DX(IIN)) + 1
          j_dy = Int((2.0_EB*0.3_EB+Max(P2P_Suunta_MAX,P2P_DIST_MAX))/DY(JJN)) + 1
          iio = Max(1,IIN-i_dx)
          iie = Min(IBAR,IIN+i_dx)
          jjo = Max(1,JJN-j_dy)
          jje = Min(JBAR,JJN+j_dy)
          ie_max = 0
          Do iix = iio, iie
             Do jjy = jjo, jje
                BG_Loop: Do j = 1, BLOCK_GRID_N(iix,jjy)
                   ie_max = ie_max + 1
                   BLOCK_LIST(Min(ie_max,bl_max)) = BLOCK_GRID(iix,jjy,j)
                End Do BG_Loop
             End Do
          End Do
          If (ie_max > bl_max ) Then
             Write(MESSAGE,'(A,2I6)') 'ERROR: Evacuate_Humans, ie_max, bl_max ', ie_max, bl_max
             Call SHUTDOWN(MESSAGE)
          End If

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
          If (NM_STRS_MESH) Then
             STRS_Indx = N
             STRP=>EVAC_STRS(N)     
          End If
          ! Check if an agent is on a spectator stand.
          SS_Loop2: Do j = 1, n_sstands
             ESS => EVAC_SSTANDS(j)
             If (ESS%IMESH == nm .And. (ESS%X1 <= HR%X .And. ESS%X2 >= HR%X) .And. (ESS%Y1 <= HR%Y .And. ESS%Y2 >= HR%Y) ) Then
                cos_x = ESS%cos_x
                cos_y = ESS%cos_y
                FAC_V0_UP   = ESS%FAC_V0_UP
                FAC_V0_DOWN = ESS%FAC_V0_DOWN
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
                Select Case (ESS%IOR)
                Case(-1)
                   speed_xm = cos_x*HR%Speed*FAC_V0_DOWN
                   speed_xp = cos_x*HR%Speed*FAC_V0_UP
                   speed_ym = HR%Speed*ESS%FAC_V0_HORI
                   speed_yp = HR%Speed*ESS%FAC_V0_HORI
                   HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + (ESS%H-ESS%H0)*Abs(ESS%X1-HR%X)/Abs(ESS%X1-ESS%X2)
                Case(+1)
                   speed_xm = cos_x*HR%Speed*FAC_V0_UP
                   speed_xp = cos_x*HR%Speed*FAC_V0_DOWN
                   speed_ym = HR%Speed*ESS%FAC_V0_HORI
                   speed_yp = HR%Speed*ESS%FAC_V0_HORI
                   HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + (ESS%H-ESS%H0)*Abs(ESS%X2-HR%X)/Abs(ESS%X1-ESS%X2)
                Case(-2)
                   speed_xm = HR%Speed*ESS%FAC_V0_HORI
                   speed_xp = HR%Speed*ESS%FAC_V0_HORI
                   speed_ym = cos_y*HR%Speed*FAC_V0_DOWN
                   speed_yp = cos_y*HR%Speed*FAC_V0_UP
                   HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + (ESS%H-ESS%H0)*Abs(ESS%Y1-HR%Y)/Abs(ESS%Y1-ESS%Y2)
                Case(+2)
                   speed_xm = HR%Speed*ESS%FAC_V0_HORI
                   speed_xp = HR%Speed*ESS%FAC_V0_HORI
                   speed_ym = cos_y*HR%Speed*FAC_V0_UP
                   speed_yp = cos_y*HR%Speed*FAC_V0_DOWN
                   HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + (ESS%H-ESS%H0)*Abs(ESS%Y2-HR%Y)/Abs(ESS%Y1-ESS%Y2)
                End Select
                Exit SS_Loop2
             End If
             HR%Z = 0.5_EB*(ZS+ZF)  ! The agent is not on any incline
          End Do SS_Loop2
          ! =======================================================
          ! Inclines Ends
          ! =======================================================

          ! =======================================================
          ! Set height and speed for an agent in stairs
          ! =======================================================
          If (NM_STRS_MESH) Then 
            Call GetStairSpeedAndZ(speed_xm, speed_xp, speed_ym, speed_yp,STRP,HR)
          End If

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
          If (TAU_CHANGE_V0 > 1.0E-12_EB) Then
             Call Random_number(rnCF)
          Else
             rnCF = -1.0_EB
             angle_old = 0.0_EB
             commitment = 0.0_EB
          End If
          Change_v0_rnCF0: If ( rnCF > Exp(-DTSP/TAU_CHANGE_V0) ) Then

             v_hr = Max(0.1_EB,Min(1.0_EB,Sqrt(HR%U**2 + HR%V**2)/HR%Speed))
             theta_start = -Abs(THETA_SECTOR)
             theta_step = 2.0_EB*Abs(theta_start)/Real(n_sectors-1,EB)
             evel = Max(0.0_EB,Min(1.0_EB,(HR%U**2+HR%V**2)/HR%Speed))
             theta_start = theta_start - 0.5_EB*(1.0_EB-evel)*Max(0.0_EB,(90.0_EB+theta_start-0.5_EB*theta_step))
             theta_step = 2.0_EB*Abs(theta_start)/Real(n_sectors-1,EB)
             If(HR%UBAR**2+HR%VBAR**2 < 0.1_EB) Then
                HR%UBAR = UBAR; HR%VBAR = VBAR
             End If
             Sum_suunta = 0.0_EB
             Do iii = 1, n_sectors
                thetas(iii)    = -Abs(theta_start) + (iii-1)*theta_step ! Degrees
                cos_theta(iii) = Cos(Pi*thetas(iii)/180._EB) ! Radians in Fortran functions
                sin_theta(iii) = Sin(Pi*thetas(iii)/180._EB)
                u_theta(iii) = cos_theta(iii)*UBAR - sin_theta(iii)*VBAR
                v_theta(iii) = sin_theta(iii)*UBAR + cos_theta(iii)*VBAR
             End Do
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
          End If Change_v0_rnCF0

          ! The force loop over the other agents starts here
          P2PLOOP: Do IE = 1, ie_max
             If (BLOCK_LIST(ie) == I) Cycle P2PLOOP
             HRE => HUMAN(BLOCK_LIST(ie))
             ! In stairs, only consider humans at the same, next and previous sub nodes (landings, stairs)
             If (ABS(HRE%STR_SUB_INDX - HR%STR_SUB_INDX)>1) Cycle P2PLOOP
             P2P_DIST = ((HRE%X-X1)**2 + (HRE%Y-Y1)**2)
             If (P2P_DIST > (Max(P2P_Suunta_Max,P2P_DIST_MAX)+HR%Radius+HRE%Radius)**2) Cycle P2PLOOP
             P2P_DIST = Sqrt(P2P_DIST)
             !
             ! Check, that the persons are seeing each other, i.e., there are no walls between.
             PP_see_each = See_each_other(nm, x1, y1, HRE%X, HRE%Y)
             If (.Not. PP_see_each) Cycle P2PLOOP

             ! Collision avoidance, counterflow, etc.
             Change_v0_rnCF1: If ( rnCF > Exp(-DTSP/TAU_CHANGE_V0) ) Then
                EVEL = Sqrt(UBAR**2 + VBAR**2)
                If (VBAR >= 0.0_EB) Then
                   angle_hr = Acos(UBAR/EVEL)
                Else
                   angle_hr = 2.0_EB*Pi - Acos(UBAR/EVEL)
                End If
                If (angle_hr == 2.0_EB*Pi) angle_hr = 0.0_EB  ! agent HR angle is [0,2Pi)
                
                d_shift = 2.0_EB*(0.5_EB-Min(1.0_EB,Sqrt(HR%U**2+HR%V**2)/HR%Speed))*HR%Radius
                hr_x = HR%X - UBAR*d_shift
                hr_y = HR%Y - VBAR*d_shift
                d_new = Sqrt((HRE%X-hr_x)**2 + (HRE%Y-hr_y)**2)
                If (d_shift < 0.0_EB .Or. d_new > HR%r_torso) Then
                   If ((HRE%Y-HR_Y) >= 0.0_EB) Then
                      angle_hre = Acos((HRE%X-HR_X)/d_new)
                   Else
                      angle_hre = 2.0_EB*Pi - Acos((HRE%X-HR_X)/d_new)
                   End If
                   If (angle_hre == 2.0_EB*Pi) angle_hre = 0.0_EB  ! agent HRE angle is [0,2Pi)
                   angle_hre = angle_hr - angle_hre
                   
                   If (angle_hre >= Pi) angle_hre = 2.0_EB*Pi - angle_hre
                   If (angle_hre <= -Pi) angle_hre = 2.0_EB*Pi + angle_hre
                   ! Agent HRE is at this angle when measured from HR
                   ! If HRE is on the right hand side of HR then the angle is negative,
                   ! i.e., positive direction is anti-clockwise (as usual).
                   angle_hre = -180.0_EB*angle_hre/Pi  ! degrees
                   
                   v_hre = (HRE%X-HR%X)*UBAR + (HRE%Y-HR%Y)*VBAR
                   v_hre = v_hre/Sqrt((HRE%X-HR%X)**2+(HRE%Y-HR%Y)**2) ! cos r_HRE vs v0
                   v_hre = Max(0.5_EB,v_hre)  ! look further ahead than sideways
                   evel = Max(0.2_EB,Min(1.0_EB,Sqrt(HR%U**2+HR%V**2)/HR%Speed))
                   tim_dist = Max(0.2_EB,Min(evel,v_hre))*P2P_Suunta_MAX + HR%Radius + HRE%Radius
                   evel = Max(0.5_EB,Min(1.0_EB,Sqrt(HR%U**2+HR%V**2)/HR%Speed))
                   v_hre = Max(0.5_EB,Min(evel,v_hre))*P2P_Suunta_MAX ! look further ahead than sideways
                   If (P2P_DIST < (v_hre+HR%Radius+HRE%Radius) ) Then
                      vr_2r = HRE%UBAR*UBAR+HRE%VBAR*VBAR ! counterflow or not?
                      ! Which sector if any
                      iii = Int((angle_hre-(theta_start-1.5_EB*theta_step))/theta_step)
                      If (iii > 0 .And. iii < n_sectors+1) Then
                         ! (UBAR,VBAR) are unit vectors
                         If (P2P_DIST < tim_dist) N_suunta(iii) = N_suunta(iii) + 1
                         If (P2P_DIST < tim_dist .And. vr_2r <= -0.2_EB) N_suuntaCF(iii) = N_suuntaCF(iii) + 1
                         If (vr_2r > 0.0_EB) Then ! Same direction
                            v_hre = HRE%U*UBAR + HRE%V*VBAR  ! HRE speed along the v0 direction
                            v_hr  = Max(0.0_EB,HR%U*UBAR + HR%V*VBAR)  ! HR speed along the v0 direction
                            v_hre = CONST_DF + FAC_DF*(Min(HR%v0_fac*HR%Speed,v_hre) - &
                                 Min(HR%v0_fac*HR%Speed,v_hr))
                         Else ! Counterflow
                            v_hre = HRE%U*UBAR + HRE%V*VBAR  ! HRE speed along the v0 direction
                            v_hre = -1.0_EB*(CONST_CF + FAC_CF*Max(0.0_EB,-v_hre))
                         End If
                         Sum_suunta(iii) = Sum_suunta(iii) + v_hre/Max(0.2_EB,(P2P_DIST-HR%Radius-HRE%Radius))
                         If (angle_hre > -0.5_EB*theta_step .And. angle_hre < 0.5_EB*theta_step) Then
                            ! The "additional" sector pointing on the v0 direction
                            N_suunta(n_sectors+1) = N_suunta(n_sectors+1) + 1
                            If (vr_2r <= -0.2_EB) N_suuntaCF(n_sectors+1) = N_suuntaCF(n_sectors+1) + 1
                            Sum_suunta(n_sectors+1) = Sum_suunta(n_sectors+1) + &
                                 v_hre/Max(0.2_EB,(P2P_DIST-HR%Radius-HRE%Radius))
                         End If
                      End If
                   End If
                End If
             End If Change_v0_rnCF1
             If ( P2P_DIST > (P2P_DIST_MAX+HR%Radius+HRE%Radius)**2 ) Cycle P2PLOOP
             ! 
             ! ========================================================
             ! Calculate the combination of spring constant for the two agents
             ! ========================================================
             C_Yeff = (2.0_EB*HR%C_Young*2.0_EB*HRE%C_Young)/(2.0_EB*HR%C_Young+2.0_EB*HRE%C_Young)
             !
             ! ========================================================
             ! Angle dependent social force:
             ! ========================================================
             If ( (HR%U**2 +HR%V**2) > 0.0_EB ) Then
                CosPhiFac = ( (HRE%X-X1)*HR%U + (HRE%Y-Y1)*HR%V ) &
                     / ( Sqrt((HRE%X-X1)**2 + (HRE%Y-Y1)**2)*Sqrt(HR%U**2 +HR%V**2) )
                CosPhiFac = HR%Lambda + 0.5_EB*(1.0_EB-HR%Lambda)*(1.0_EB+CosPhiFac)
             Else
                CosPhiFac = 1.0_EB
             End If

             ! Calculate the position and velocities of the shoulder cirles for HRE
             r_tmp(4) = HRE%r_shoulder ! right circle
             r_tmp(5) = HRE%r_torso    ! center circle
             r_tmp(6) = HRE%r_shoulder ! left circle
             y_tmp(4) = HRE%Y - Cos(HRE%angle)*HRE%d_shoulder ! right circle
             x_tmp(4) = HRE%X + Sin(HRE%angle)*HRE%d_shoulder ! right circle
             y_tmp(5) = HRE%Y ! center circle
             x_tmp(5) = HRE%X ! center circle
             y_tmp(6) = HRE%Y + Cos(HRE%angle)*HRE%d_shoulder ! left circle
             x_tmp(6) = HRE%X - Sin(HRE%angle)*HRE%d_shoulder ! left circle

             ! ========================================================
             ! Add psychological (social) force term
             ! ========================================================
             If (.Not. L_Dead) Then
                Fc_x = 0.0_EB
                Fc_y = 0.0_EB
                ! Use the closest circles to calculate the psychological force
                Do iii = 1, 3
                   Do jjj = 4, 6
                      tim_dist = Sqrt((x_tmp(iii)-x_tmp(jjj))**2 + (y_tmp(iii)-y_tmp(jjj))**2)
                      ! d_humans = Min( tim_dist-(r_tmp(iii)+r_tmp(jjj)) , d_humans )
                      ! Next is |vector1|*|vector2|
                      evel = Sqrt((x_tmp(jjj)-x_tmp(iii))**2+(y_tmp(jjj)-y_tmp(iii))**2)* Sqrt(u_tmp(iii)**2+v_tmp(iii)**2)
                      If (evel > 0.0_EB) evel = ((x_tmp(jjj)-x_tmp(iii))*u_tmp(iii) + &
                           (y_tmp(jjj)-y_tmp(iii))*v_tmp(iii)) / evel   ! cos theta (scal_prod/(lenght1*length2)
                      If (evel > 0.01_EB) Then
                         d_humans = Min( (tim_dist-(r_tmp(iii)+r_tmp(jjj))) /evel, d_humans)
                      Else
                         d_humans = Min( (tim_dist-(r_tmp(iii)+r_tmp(jjj))) /0.01_EB , d_humans)
                      End If

                      Fc_x1 = (x_tmp(iii)-x_tmp(jjj))*hr_a*CosPhiFac* &
                           Exp( -(tim_dist-( r_tmp(iii)+r_tmp(jjj) ))/hr_b )/tim_dist 
                      Fc_y1 = (y_tmp(iii)-y_tmp(jjj))*hr_a*CosPhiFac* &
                           Exp( -(tim_dist-( r_tmp(iii)+r_tmp(jjj) ))/hr_b )/tim_dist 
                      If ( (Fc_x1**2+Fc_y1**2) > (Fc_x**2+Fc_y**2) ) Then
                         Fc_x = Fc_x1
                         Fc_y = Fc_y1
                      End If
                   End Do
                End Do
                P2P_U = P2P_U + Fc_x
                P2P_V = P2P_V + Fc_y
                Social_F = Social_F + Sqrt(Fc_x**2 + Fc_y**2)
                Tc_z = 0.0_EB
                ! Calculate the torque due to the social force. Use the closest circles.
                Do jjj = 4, 6
                   ! First the right shoulder
                   tim_dist = Sqrt( (x_tmp(jjj)-x_tmp(1))**2 + (y_tmp(jjj)-y_tmp(1))**2 )
                   Fc_x = (x_tmp(1)-x_tmp(jjj)) * hr_a*CosPhiFac*Exp( -(tim_dist - (r_tmp(1)+r_tmp(jjj)))/hr_b )/tim_dist
                   Fc_y = (y_tmp(1)-y_tmp(jjj)) * hr_a*CosPhiFac*Exp( -(tim_dist - (r_tmp(1)+r_tmp(jjj)))/hr_b )/tim_dist
                   If ( Abs(Fc_y*(x_tmp(1)-HR%X) - Fc_x*(y_tmp(1)-HR%Y)) > Abs(Tc_z) ) Then
                      Tc_z = Fc_y*(x_tmp(1)-HR%X) - Fc_x*(y_tmp(1)-HR%Y)
                   End If
                End Do
                P2P_Torque = P2P_Torque + Tc_z
                Tc_z = 0.0_EB
                Do jjj = 4, 6
                   ! Then the left shoulder
                   tim_dist = Sqrt( (x_tmp(jjj)-x_tmp(3))**2 + (y_tmp(jjj)-y_tmp(3))**2 )
                   Fc_x = (x_tmp(3)-x_tmp(jjj)) * hr_a*CosPhiFac*Exp( -(tim_dist - (r_tmp(3)+r_tmp(jjj)))/hr_b )/tim_dist
                   Fc_y = (y_tmp(3)-y_tmp(jjj)) * hr_a*CosPhiFac*Exp( -(tim_dist - (r_tmp(3)+r_tmp(jjj)))/hr_b )/tim_dist
                   If ( Abs(Fc_y*(x_tmp(3)-HR%X) - Fc_x*(y_tmp(3)-HR%Y)) > Abs(Tc_z) ) Then
                      Tc_z = Fc_y*(x_tmp(3)-HR%X) - Fc_x*(y_tmp(3)-HR%Y)
                   End If
                End Do
                P2P_Torque = P2P_Torque + Tc_z
             End If  ! dead or alive?

             ! ========================================================
             ! Add contact force terms
             ! ========================================================
             If ( P2P_DIST <= (HR%Radius+HRE%Radius) ) Then
                ! Calculate the velocities of the shoulder cirles
                v_tmp(4) = HRE%V + Sin(HRE%angle)*HRE%Omega*HRE%d_shoulder
                v_tmp(5) = HRE%V
                v_tmp(6) = HRE%V - Sin(HRE%angle)*HRE%Omega*HRE%d_shoulder
                u_tmp(4) = HRE%U + Cos(HRE%angle)*HRE%Omega*HRE%d_shoulder
                u_tmp(5) = HRE%U
                u_tmp(6) = HRE%U - Cos(HRE%angle)*HRE%Omega*HRE%d_shoulder

                Do iii = 1, 3
                   Do jjj = 4, 6
                      tim_dist = Sqrt((x_tmp(iii)-x_tmp(jjj))**2 + (y_tmp(iii)-y_tmp(jjj))**2)
                      ! d_humans = Min( tim_dist-(r_tmp(iii)+r_tmp(jjj)) , d_humans )
                      ! Next is |vector1|*|vector2|
                      evel = Sqrt((x_tmp(jjj)-x_tmp(iii))**2+(y_tmp(jjj)-y_tmp(iii))**2)* Sqrt(u_tmp(iii)**2+v_tmp(iii)**2)
                      If (evel > 0.0_EB) evel = ((x_tmp(jjj)-x_tmp(iii))*u_tmp(iii) + &
                           (y_tmp(jjj)-y_tmp(iii))*v_tmp(iii)) / evel   ! cos theta (scal_prod/(lenght1*length2)
                      If (evel > 0.01_EB) Then
                         d_humans = Min( (tim_dist-(r_tmp(iii)+r_tmp(jjj))) /evel, d_humans)
                      Else
                         d_humans = Min( (tim_dist-(r_tmp(iii)+r_tmp(jjj))) /0.01_EB , d_humans)
                      End If
                      If (tim_dist <= r_tmp(iii)+r_tmp(jjj) ) Then
                         ! Circles are touching each others
                         Fc_x =(x_tmp(iii)-x_tmp(jjj))*C_Yeff*((r_tmp(iii)+r_tmp(jjj))-tim_dist)/tim_dist
                         Fc_y =(y_tmp(iii)-y_tmp(jjj))*C_Yeff*((r_tmp(iii)+r_tmp(jjj))-tim_dist)/tim_dist
                         Fc_x = Fc_x - FC_DAMPING*(u_tmp(iii)-u_tmp(jjj))*(x_tmp(iii)-x_tmp(jjj))/tim_dist
                         Fc_y = Fc_y - FC_DAMPING*(v_tmp(iii)-v_tmp(jjj))*(y_tmp(iii)-y_tmp(jjj))/tim_dist
                         Contact_F = Contact_F + Sqrt(Fc_x**2 + Fc_y**2)
                         P2P_U = P2P_U + Fc_x
                         P2P_V = P2P_V + Fc_y
                         P2P_Torque = P2P_Torque + Fc_y*(x_tmp(iii)-HR%X) - Fc_x*(y_tmp(iii)-HR%Y)
                         scal_prod_over_rsqr = ((y_tmp(iii)-y_tmp(jjj))*(u_tmp(iii)-u_tmp(jjj)) - &
                              (x_tmp(iii)-x_tmp(jjj))*(v_tmp(iii)-v_tmp(jjj))) / (tim_dist**2)
                         If (I_Fric_sw >= 1 ) Then  ! This is the default
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
                         Else
                            Fc_x = -HR%Gamma*( (y_tmp(iii)-y_tmp(jjj))*scal_prod_over_rsqr)
                            Fc_y = -HR%Gamma*(-(x_tmp(iii)-x_tmp(jjj))*scal_prod_over_rsqr)
                            P2P_U = P2P_U + Fc_x
                            P2P_V = P2P_V + Fc_y
                            P2P_Torque = P2P_Torque + Fc_y*( (x_tmp(iii) + &
                                 (r_tmp(iii)/r_tmp(jjj))*(x_tmp(jjj)-x_tmp(iii)) ) - HR%X ) 
                            P2P_Torque = P2P_Torque - Fc_x*( (y_tmp(iii) + &
                                 (r_tmp(iii)/r_tmp(jjj))*(y_tmp(jjj)-y_tmp(iii)) ) - HR%Y ) 
                         End If
                      End If
                   End Do
                End Do
             End If  ! contact forces?
          End Do P2PLOOP
          TUSED(14,NM)=TUSED(14,NM)+SECOND()-TNOW14
          ! ========================================================
          ! PERSON-PERSON INTERACTION FORCES ENDS HERE
          ! ========================================================

          ! ========================================================
          ! THE PERSON-WALL FORCES
          ! ========================================================
          !
          ! Walls are looked for the body circle
          Call Find_walls(nm, x1, y1, HR%Radius, P2P_DIST_MAX, HR%SKIP_WALL_FORCE_IOR, d_xy, FoundWall_xy, istat)

          ! ========================================================
          ! Collision avoidance, counterflow, etc.
          ! ========================================================
          Change_v0_rnCF2: If ( rnCF > Exp(-DTSP/TAU_CHANGE_V0) ) Then
             v_hr  = Max(0.1_EB,Sqrt(HR%U**2+HR%V**2)/HR%Speed)
             tim_dist = 0.0_EB
             Do iii = 1, n_sectors
                ! Awoid walls, do not take a direction where there is a wall closer than
                ! d_perp = 0.6 m (perpendicular).
                If (Abs(sin_theta(iii)) > 0.0001_EB) Then
                   x11 = HR%X + u_theta(iii)*Min(P2P_Suunta_MAX, 0.6_EB/Abs(sin_theta(iii)))
                   y11 = HR%Y + v_theta(iii)*Min(P2P_Suunta_MAX, 0.6_EB/Abs(sin_theta(iii)))
                Else ! straight ahead
                   x11 = HR%X + u_theta(iii)*P2P_Suunta_MAX
                   y11 = HR%Y + v_theta(iii)*P2P_Suunta_MAX
                End If
                P2P_DIST = Sqrt((HR%X-x11)**2 + (HR%Y-y11)**2) - HR%Radius
                x11 = Max(XS,Min(XF,x11))
                y11 = Max(YS,Min(YF,y11))
                PP_see_each = See_each_other(nm, HR%X, HR%Y, x11, y11)
                If(.Not.PP_see_each) Then
                   vr_2r = -FAC_1_WALL*P2P_Suunta_MAX*v_hr
                   Sum_suunta(iii) = Sum_suunta(iii) + vr_2r/Max(0.2_EB,P2P_DIST)
                End If
             End Do

             If (FoundWall_xy(1)) tim_dist = tim_dist - 3.0_EB/Max(0.3_EB,Abs(HR%X-d_xy(1))-HR%Radius)
             If (FoundWall_xy(2)) tim_dist = tim_dist - 3.0_EB/Max(0.3_EB,Abs(HR%X-d_xy(2))-HR%Radius)
             If (FoundWall_xy(3)) tim_dist = tim_dist - 3.0_EB/Max(0.3_EB,Abs(HR%Y-d_xy(3))-HR%Radius)
             If (FoundWall_xy(4)) tim_dist = tim_dist - 3.0_EB/Max(0.3_EB,Abs(HR%Y-d_xy(4))-HR%Radius)
!!$             If (FoundWall_xy(1)) tim_dist = tim_dist - 3.0_EB/Max(0.3_EB,Abs(HR%X-d_xy(1))-HR%Radius)
!!$             If (FoundWall_xy(2)) tim_dist = tim_dist - 3.0_EB/Max(0.3_EB,Abs(HR%X-d_xy(2))-HR%Radius)
!!$             If (FoundWall_xy(3)) tim_dist = tim_dist - 3.0_EB/Max(0.3_EB,Abs(HR%Y-d_xy(3))-HR%Radius)
!!$             If (FoundWall_xy(4)) tim_dist = tim_dist - 3.0_EB/Max(0.3_EB,Abs(HR%Y-d_xy(4))-HR%Radius)

             If (FoundWall_xy(1) .And. Abs(HR%X-d_xy(1))-HR%Radius < 0.1_EB) Then
                Do iii = 1, n_sectors
                   If(u_theta(iii) < -0.10_EB) Then
                      vr_2r = u_theta(iii)/Max(0.1_EB, Abs(HR%X-HR%Radius-d_xy(1)))
                      Sum_suunta(iii) = Sum_suunta(iii) - Abs(P2P_Suunta_MAX*FAC_2_WALL*vr_2r)
                   End If
                End Do
             End If
             If (FoundWall_xy(2) .And. Abs(HR%X-d_xy(2))-HR%Radius < 0.1_EB) Then
                Do iii = 1, n_sectors
                   If(u_theta(iii) > +0.10_EB) Then
                      vr_2r = -u_theta(iii)/Max(0.1_EB, Abs(HR%X+HR%Radius-d_xy(2)))
                      Sum_suunta(iii) = Sum_suunta(iii) - Abs(P2P_Suunta_MAX*FAC_2_WALL*vr_2r)
                   End If
                End Do
             End If
             If (FoundWall_xy(3) .And. Abs(HR%Y-d_xy(3))-HR%Radius < 0.1_EB) Then
                Do iii = 1, n_sectors
                   If(v_theta(iii) < -0.10_EB) Then
                      vr_2r = v_theta(iii)/Max(0.1_EB, Abs(HR%Y-HR%Radius-d_xy(3)))
                      Sum_suunta(iii) = Sum_suunta(iii) - Abs(P2P_Suunta_MAX*FAC_2_WALL*vr_2r)
                   End If
                End Do
             End If
             If (FoundWall_xy(4) .And. Abs(HR%Y-d_xy(4))-HR%Radius < 0.1_EB) Then
                Do iii = 1, n_sectors
                   If(v_theta(iii) > +0.10_EB) Then
                      vr_2r = -v_theta(iii)/Max(0.1_EB, Abs(HR%Y+HR%Radius-d_xy(4)))
                      Sum_suunta(iii) = Sum_suunta(iii) - Abs(P2P_Suunta_MAX*FAC_2_WALL*vr_2r)
                   End If
                End Do
             End If
!!$             Do iii = 1, n_sectors+1
!!$                ! Prefer flow field direction
!!$                Sum_suunta(iii) = Sum_suunta(iii) + FAC_V0_DIR*v_hr*((UBAR*u_theta(iii)+VBAR*v_theta(iii)))
!!$             End Do
             Sum_suunta(n_sectors+1) = Sum_suunta(n_sectors+1) + Abs(FAC_V0_DIR)*v_hr
             Do iii = 1, n_sectors
                If (N_suuntaCF(n_sectors+1) < 1) Then
                   ! No counterflow: Prefer left (and straight ahead)
                   ! Sum_suunta(iii) = Sum_suunta(iii) + Sign(1.0_EB,thetas(iii))* &
                   !     FAC_V0_DIR*v_hr
                Else
                   ! Counterflow: Prefer right (and straight ahead)
                   Sum_suunta(iii) = Sum_suunta(iii) - Sign(1.0_EB,thetas(iii))*FAC_V0_DIR
                End If
             End Do

             If (N_suunta(n_sectors+1) < 1) Sum_suunta(n_sectors+1) = Sum_suunta(n_sectors+1) + 50.0_EB  ! Empty space ahead
             If ((N_suuntaCF(n_sectors+1)) < 1) Then
                ! No counterflow, prefer v0 direction, i.e., "stay on line"
                Sum_suunta(n_sectors+1) = Sum_suunta(n_sectors+1) + &
                     (FAC_NOCF + FAC_V0_NOCF*v_hr)*(N_suunta(n_sectors+1))
             End If
             i_suunta_max = n_sectors + 1
             Sum_suunta_max = -Huge(Sum_suunta_max)
             Do iii = 1, n_sectors+1
                If (Sum_suunta(iii) > Sum_suunta_max) Then
                   UBAR = u_theta(iii)
                   VBAR = v_theta(iii)
                   Sum_suunta_max = Sum_suunta(iii)
                   i_suunta_max = iii
                End If
             End Do
             ! If counterflow then try to pass sideways and decrease tau and tau_iner.
             ! Use this if there are more counterflow agents than downstrean agents.
             angle_old = 0.0_EB
             commitment = 0.0_EB
             If (v_hr < 0.3_EB .And. (N_suunta(n_sectors+1)-N_suuntaCF(n_sectors+1)) < &
                  N_suuntaCF(n_sectors+1)-0) Then
                commitment = Max(0.5_EB,Real(Sum(N_suuntaCF(1:n_sectors+1)),EB)/Max(1,Sum(N_suunta(1:n_sectors+1))))
                angle_old = -Sign(1.0_EB,thetas(i_suunta_max))*Max(0.1_EB,Abs(thetas(i_suunta_max)))*Pi/180.0_EB 
                angle_old = angle_old + 85.0_EB*Pi/180.0_EB
             Else If (v_hr < 0.3_EB .And. Sum(N_suuntaCF(1:n_sectors)) >= 1 .And. tim_dist < -15.0_EB) Then
                commitment = Max(0.5_EB,Real(Sum(N_suuntaCF(1:n_sectors)),EB)/Max(1,Sum(N_suunta(1:n_sectors))))
                angle_old = -Sign(1.0_EB,thetas(i_suunta_max))*Max(0.1_EB,Abs(thetas(i_suunta_max)))*Pi/180.0_EB 
                angle_old = angle_old + 85.0_EB*Pi/180.0_EB
             Else
                angle_old = 0.0_EB
                commitment = 0.0_EB
             End If
          Else ! Change_v0_rnCF2
             ! Do not change direction during this time step, use the previous direction
             If (TAU_CHANGE_V0 > 1.0E-12_EB) Then
                UBAR = HR%UBAR
                VBAR = HR%VBAR
                angle_old = HR%angle_old
                commitment = HR%Commitment
             Else
                HR%UBAR = UBAR
                HR%VBAR = VBAR
                angle_old = 0.0_EB
                commitment = 0.0_EB
             End If
          End If Change_v0_rnCF2 ! Collision avoidance ends
          HR%angle_old = angle_old
          HR%Commitment = commitment

          Call Wall_SocialForces(nm, x_tmp, y_tmp, r_tmp, P2P_DIST_MAX, d_xy, P2P_U, P2P_V, Social_F, FoundWall_xy)

          Call Wall_ContactForces(nm, x_tmp(1), y_tmp(1), r_tmp(1), u_tmp(1), v_tmp(1), d_xy, &
               P2P_U, P2P_V, P2P_Torque, Contact_F, d_walls, FoundWall_xy)
          Call Wall_ContactForces(nm, x_tmp(2), y_tmp(2), r_tmp(2), u_tmp(2), v_tmp(2), d_xy, &
               P2P_U, P2P_V, P2P_Torque, Contact_F, d_walls, FoundWall_xy)
          Call Wall_ContactForces(nm, x_tmp(3), y_tmp(3), r_tmp(3), u_tmp(3), v_tmp(3), d_xy, &
               P2P_U, P2P_V, P2P_Torque, Contact_F, d_walls, FoundWall_xy)

          ! Add forces from the door case
          Call Door_Forces(nm, x_tmp, y_tmp, r_tmp, u_tmp, v_tmp, p2p_dist_max, d_xy,&
               P2P_U, P2P_V, Social_F, Contact_F, P2P_Torque, FoundWall_xy)

          ! Add wall corner - person forces
          ! top right corner (x > x_human, y > y_human)
          x_now = -DX(iin+1)
          Loop_px: Do ii = iin, IBAR, +1
             x_now = x_now + DX(ii)
             If (x_now-HR%Radius > P2P_DIST_MAX) Exit Loop_px
             y_now = -DY(jjn-1)
             Loop_pxpy: Do jj = jjn, JBAR
                y_now = y_now + DY(jj)
                If (Sqrt(x_now**2 + y_now**2)-HR%Radius > P2P_DIST_MAX) Exit Loop_pxpy
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
                If (tim_iwy /= 0) Exit Loop_pxpy  ! 
                If ( (tim_iwx==0).And.(tim_iwy==0).And.(tim_iw/=0 .Or. tim_iw2/=0) ) Then
                   If (tim_iw/=0) Then
                      ! first y-direction then x-direction
                      x11 = xw(tim_iw )                 ! corner point x
                      y11 = yw(tim_iw )-0.5_EB*DY(jj+1) ! corner point y
                   Else
                      ! first x-direction then y-direction
                      x11 = xw(tim_iw2)-0.5_EB*DX(ii+1) ! corner point x
                      y11 = yw(tim_iw2)                 ! corner point y
                   End If

                   Call Corner_Forces(x1, y1, x11, y11, p2p_dist_max, P2P_U, P2P_V, Social_F, &
                        Contact_F, P2P_Torque, d_walls, x_tmp, y_tmp, r_tmp, u_tmp, v_tmp, istat)
                   Exit Loop_pxpy
                End If
             End Do Loop_pxpy
          End Do Loop_px

          ! top left corner (x < x_human, y > y_human)
          x_now = -DX(iin-1)
          Loop_mx: Do ii = iin, 1, -1
             x_now = x_now + DX(ii)
             If (x_now-HR%Radius > P2P_DIST_MAX) Exit Loop_mx
             y_now = -DY(jjn)
             Loop_mxpy: Do jj = jjn, JBAR
                y_now = y_now + DY(jj)
                If (Sqrt(x_now**2 + y_now**2)-HR%Radius > P2P_DIST_MAX) Exit Loop_mxpy
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
                If (tim_iwy /= 0) Exit Loop_mxpy
                If ( (tim_iwx==0).And.(tim_iwy==0).And.(tim_iw/=0 .Or. tim_iw2/=0) ) Then
                   If (tim_iw/=0) Then
                      ! first y-direction then x-direction
                      x11 = xw(tim_iw )                 ! corner point x
                      y11 = yw(tim_iw )-0.5_EB*DY(jj+1) ! corner point y
                   Else
                      ! first x-direction then y-direction
                      x11 = xw(tim_iw2)+0.5_EB*DX(ii-1) ! corner point x
                      y11 = yw(tim_iw2)                 ! corner point y
                   End If

                   Call Corner_Forces(x1, y1, x11, y11, p2p_dist_max, P2P_U, P2P_V, Social_F, &
                        Contact_F, P2P_Torque, d_walls, x_tmp, y_tmp, r_tmp, u_tmp, v_tmp, istat)
                   Exit Loop_mxpy
                End If
             End Do Loop_mxpy
          End Do Loop_mx

          ! bottom right corner (x > x_human, y < y_human)
          x_now = -DX(iin+1)
          Loop_py: Do ii = iin, IBAR, +1
             x_now = x_now + DX(ii)
             If (x_now-HR%Radius > P2P_DIST_MAX) Exit Loop_py
             y_now = -DY(jjn-1)
             Loop_pxmy: Do jj = jjn, 1, -1
                y_now = y_now + DY(jj)
                If (Sqrt(x_now**2 + y_now**2)-HR%Radius > P2P_DIST_MAX) Exit Loop_pxmy
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
                I_OBST = OBST_INDEX_W(tim_IW)
                IF (tim_iw>0) THEN
                   IF (OBSTRUCTION(I_OBST)%HIDDEN) tim_iw = 0
                END IF
                tim_ic2 = cell_index(ii+1,jj,kkn) ! one cell right
                tim_iw2 = wall_index(tim_ic2, -2) ! right and down
                IF (tim_iw2>0) THEN
                   I_OBST = OBST_INDEX_W(tim_IW2)
                   IF (OBSTRUCTION(I_OBST)%HIDDEN) tim_iw2 = 0
                END IF
                If (tim_iwy /= 0) Exit Loop_pxmy
                If ( (tim_iwx==0).And.(tim_iwy==0).And.(tim_iw/=0 .Or. tim_iw2/=0) ) Then
                   If (tim_iw/=0) Then
                      ! first y-direction then x-direction
                      x11 = xw(tim_iw )                 ! corner point x
                      y11 = yw(tim_iw )+0.5_EB*DY(jj-1) ! corner point y
                   Else
                      ! first x-direction then y-direction
                      x11 = xw(tim_iw2)-0.5_EB*DX(ii+1) ! corner point x
                      y11 = yw(tim_iw2)                 ! corner point y
                   End If

                   Call Corner_Forces(x1, y1, x11, y11, p2p_dist_max, P2P_U, P2P_V, Social_F, &
                        Contact_F, P2P_Torque, d_walls, x_tmp, y_tmp, r_tmp, u_tmp, v_tmp, istat)
                   Exit Loop_pxmy
                End If
             End Do Loop_pxmy
          End Do Loop_py

          ! bottom left corner (x < x_human, y < y_human)
          x_now = -DX(iin-1)
          Loop_my: Do ii = iin, 1, -1
             x_now = x_now + DX(ii)
             If (x_now-HR%Radius > P2P_DIST_MAX) Exit Loop_my
             y_now = -DY(jjn-1)
             Loop_mxmy: Do jj = jjn, 1, -1
                y_now = y_now + DY(jj)
                If (Sqrt(x_now**2+y_now**2)-HR%Radius > P2P_DIST_MAX) Exit Loop_mxmy
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
                If (tim_iwy /= 0) Exit Loop_mxmy
                If ( (tim_iwx==0).And.(tim_iwy==0).And.(tim_iw/=0 .Or. tim_iw2/=0) ) Then
                   If (tim_iw/=0) Then
                      ! first y-direction then x-direction
                      x11 = xw(tim_iw )                 ! corner point x
                      y11 = yw(tim_iw )+0.5_EB*DY(jj-1) ! corner point y
                   Else
                      ! first x-direction then y-direction
                      x11 = xw(tim_iw2)+0.5_EB*DX(ii-1) ! corner point x
                      y11 = yw(tim_iw2)                 ! corner point y
                   End If

                   Call Corner_Forces(x1, y1, x11, y11, p2p_dist_max, P2P_U, P2P_V, Social_F, &
                        Contact_F, P2P_Torque, d_walls, x_tmp, y_tmp, r_tmp, u_tmp, v_tmp, istat)
                   Exit Loop_mxmy
                End If
             End Do Loop_mxmy
          End Do Loop_my

          ! ========================================================
          ! THE PERSON-WALL FORCES ENDS HERE
          ! ========================================================

          ! Save the forces for the loop EVAC_MOVE_LOOP (next time step)
          HR%F_X = P2P_U
          HR%F_Y = P2P_V
          HR%SumForces  = Contact_F
          HR%SumForces2 = Social_F + Contact_F
          HR%Torque = P2P_Torque

          If ( T <= T_BEGIN ) Then
             If ( Abs(P2P_U)/HR%Mass > 550.0_EB ) P2P_U =550.0_EB*HR%Mass*P2P_U/Abs(P2P_U)
             If ( Abs(P2P_V)/HR%Mass > 550.0_EB ) P2P_V =550.0_EB*HR%Mass*P2P_V/Abs(P2P_V)
             HR%F_X = P2P_U
             HR%F_Y = P2P_V
          End If

          ! Add the social and contact force terms
          U_new     = HR%U + 0.5_EB*HR%F_X*DTSP/HR%Mass
          V_new     = HR%V + 0.5_EB*HR%F_Y*DTSP/HR%Mass
          Omega_new = HR%Omega + 0.5_EB*DTSP*HR%Torque/HR%M_iner

          ! Add the effect of the random force
          If (GaTh > 0.0_EB .And. T > T_BEGIN ) Then
             U_new = U_new + 0.5_EB*DTSP*HR%v0_fac*HR%Mass*HR%ksi*Cos(HR%eta)/HR%Mass
             V_new = V_new + 0.5_EB*DTSP*HR%v0_fac*HR%Mass*HR%ksi*Sin(HR%eta)/HR%Mass
             P2P_U = P2P_U + HR%v0_fac*HR%Mass*HR%ksi*Cos(HR%eta)
             P2P_V = P2P_V + HR%v0_fac*HR%Mass*HR%ksi*Sin(HR%eta)
             Omega_new = Omega_new + 0.5_EB*DTSP*1.0_EB*Sign(HR%ksi,HR%eta-Pi)
             P2P_Torque = P2P_Torque + 1.0_EB*Sign(HR%ksi,HR%eta-Pi)*HR%M_iner
          Else
             HR%ksi = 0.0_EB
             HR%eta = 0.0_EB
          End If
          ! Now the step (4a) of SC-VV by Vattulainen is ended.

          ! Add self-propelling force terms, self-consistent VV 
          ! (first time step towards the exit door)
          fac_tim =  1.0_EB + (DTSP/(2.0_EB*hr_tau))
          If ( T <= Tpre ) Then
             If ( (T+DTSP_new) > Tpre) Then
                EVEL = Sqrt(UBAR**2 + VBAR**2)
                If (EVEL > 0.0_EB) Then
                   speed = speed_xp*(0.5_EB + Sign(0.5_EB,UBAR)) + speed_xm*(0.5_EB - Sign(0.5_EB,UBAR)) 
                   speed = speed*HR%v0_fac
                   P2P_U = P2P_U + (HR%Mass/hr_tau)*speed*(UBAR/EVEL)
                   speed = speed_yp*(0.5_EB + Sign(0.5_EB,VBAR)) + speed_ym*(0.5_EB - Sign(0.5_EB,VBAR)) 
                   speed = speed*HR%v0_fac
                   P2P_V = P2P_V + (HR%Mass/hr_tau)*speed*(VBAR/EVEL)
                End If
             End If
             UBAR = 0.0_EB
             VBAR = 0.0_EB
          End If

          ! Add self-propelling force terms, self-consistent VV
          EVEL = Sqrt(UBAR**2 + VBAR**2)
          If (EVEL > 0.0_EB) Then
             speed = speed_xp*(0.5_EB + Sign(0.5_EB,UBAR)) + speed_xm*(0.5_EB - Sign(0.5_EB,UBAR)) 
             speed = speed*HR%v0_fac
             U_new = (U_new + 0.5_EB*(DTSP/hr_tau)*speed*(UBAR/EVEL)) / fac_tim
             HR%UBAR = UBAR
             P2P_U = P2P_U + (HR%Mass/hr_tau)*(speed*(UBAR/EVEL) - HR%U)

             speed = speed_yp*(0.5_EB + Sign(0.5_EB,VBAR)) + speed_ym*(0.5_EB - Sign(0.5_EB,VBAR)) 
             speed = speed*HR%v0_fac
             V_new = (V_new + 0.5_EB*(DTSP/hr_tau)*speed*(VBAR/EVEL)) / fac_tim
             HR%VBAR = VBAR
             P2P_V = P2P_V + (HR%Mass/hr_tau)*(speed*(VBAR/EVEL) - HR%V)

             If (VBAR >= 0.0_EB) Then
                angle = Acos(UBAR/EVEL)
             Else
                angle = 2.0_EB*Pi - Acos(UBAR/EVEL)
             End If

             ! Collision avoidance has HR%angle_old .ne. 0.0
             angle = angle + HR%angle_old
             Do While (angle >= 2.0_EB*Pi)
                angle = angle - 2.0_EB*Pi
             End Do
             Do While (angle < 0.0_EB)
                angle = angle + 2.0_EB*Pi
             End Do

             If (angle == 2.0_EB*Pi) angle = 0.0_EB  ! angle is [0,2Pi)

             ! Rotational motion: J(dw/dt) = (J/t_iner)*( ((angle-angle_0)/pi)*w_0 - w)
             If (Abs(angle-HR%angle) <= Pi ) Then
                ! zero is not crossed.
                Omega_new = Omega_new + 0.5_EB*(DTSP/HR_tau_iner)*((angle-HR%angle)*(Omega_0/Pi) - HR%Omega)
                P2P_Torque = P2P_Torque + (HR%M_iner/HR_tau_iner)*((angle-HR%angle)*(Omega_0/Pi) - HR%Omega)
             Else
                ! zero is crossed
                Omega_new = Omega_new + 0.5_EB*(DTSP/HR_Tau_iner)* &
                     ( (2.0_EB*Pi-Abs(angle-HR%angle))*Sign(1.0_EB , HR%angle-angle)*(Omega_0/Pi) - HR%Omega)
                P2P_Torque = P2P_Torque + (HR%M_iner/HR_Tau_iner)* &
                     ( (2.0_EB*Pi-Abs(angle-HR%angle))*Sign(1.0_EB , HR%angle-angle)*(Omega_0/Pi) - HR%Omega)
             End If
          Else  ! No target direction
             HR%UBAR = 0.0_EB
             HR%VBAR = 0.0_EB
             U_new = U_new / fac_tim
             V_new = V_new / fac_tim
             ! Slow rotation down if no direction available, i.e., target Omega_0 is zero.
             Omega_new = Omega_new + 0.5_EB*(DTSP/HR_Tau_iner)*(-HR%Omega)
          End If

          ! Check that velocities are not too large, i.e., less than 10 m/s for humans
          EVEL = Sqrt(U_new**2 + V_new**2)
          If ( EVEL > Vmax_timo ) Then
             U_new = U_new*(Vmax_timo/EVEL)
             V_new = V_new*(Vmax_timo/EVEL)
          End If
          ! Check that angular velocity is not too large
          If ( Abs(Omega_new) > Omega_max ) Then
             Omega_new = Sign(Omega_max,Omega_new)
          End If

          HR%U = U_new
          HR%V = V_new
          HR%Omega = Omega_new

          ! ========================================================
          ! Decide the time step for the next human movement loop
          ! ========================================================
          ! Distances to closest walls and other agents
          d_humans = Max(d_humans,0.0005_EB)        ! this agent
          d_walls  = Max(d_walls, 0.0005_EB)        ! this agent
          d_humans_min = Min(d_humans_min,d_humans) ! among all agents
          d_walls_min  = Min(d_walls_min, d_walls)  ! among all agents

          If ( T > T_BEGIN ) Then
             ! Time step, do not move too close to other agents, walls, or 0.5*grid spacing.
             ! Delta_min is minimum of dx and dy for the mesh.  Minimum movement is 0.0001 m.
             dt_Loop: Do
                u_tmp(2) = HR%U + 0.5_EB*DTSP_new*P2P_U/HR%Mass
                v_tmp(2) = HR%V + 0.5_EB*DTSP_new*P2P_V/HR%Mass
                Omega_new = HR%Omega + 0.5_EB*DTSP*P2P_Torque/HR%M_iner
                u_tmp(1) = u_tmp(2) + Cos(HR%angle)*Omega_new*HR%d_shoulder
                u_tmp(3) = u_tmp(2) - Cos(HR%angle)*Omega_new*HR%d_shoulder
                v_tmp(1) = v_tmp(2) + Sin(HR%angle)*Omega_new*HR%d_shoulder
                v_tmp(3) = v_tmp(2) - Sin(HR%angle)*Omega_new*HR%d_shoulder
                If ( Max(u_tmp(1)**2+v_tmp(1)**2, u_tmp(2)**2+v_tmp(2)**2, u_tmp(3)**2+v_tmp(3)**2)*DTSP_new**2 > &
                     (Min(0.2_EB*d_humans, 0.2_EB*d_walls, 0.5_EB*Delta_min))**2 ) Then
                   DTSP_new = DTSP_new*0.8_EB
                   Cycle dt_Loop
                End If
                Exit dt_Loop
             End Do dt_Loop

             DTSP_new = Max(DTSP_new, EVAC_DT_MIN)
             DTSP_new = Min(DTSP_new, EVAC_DT_MAX)
          Else  ! Initialization phase
             HR%U     = 0.0_EB
             HR%V     = 0.0_EB
             HR%Omega = 0.0_EB
          End If

       End Do EVAC_FORCE_LOOP
       TUSED(13,NM)=TUSED(13,NM)+SECOND()-TNOW13
       Deallocate(BLOCK_LIST)
       Deallocate(BLOCK_GRID)
       ! ========================================================
       ! FORCE LOOP ENDS HERE
       ! ========================================================

       ! Save the maximum time steps allowed by the SC-VV algorithm
       Tsteps(NM) = DTSP_new

    End Do HUMAN_TIME_LOOP
    Deallocate(BLOCK_GRID_N)
    ! ========================================================
    ! HUMAN TIME LOOP ENDS HERE (Tin ==> T = Tin + DT)
    ! ========================================================

    ! ========================================================
    ! Evacuation routine ends here
    ! ========================================================
    TUSED(12,NM)=TUSED(12,NM)+SECOND()-TNOW

  Contains

    Subroutine Find_Prefered_Direction(I, Nout, T, T_BEGIN, L_Dead, NM_STRS_MESH, &
         II, JJ, IIX, JJY, XI, YJ, ZK, UBAR, VBAR, hr_tau, Tpre)
      Implicit None
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
      Integer, Intent(IN) :: II, JJ, IIX, JJY, I
      Real(EB), Intent(IN) :: XI, YJ, ZK, T, T_BEGIN
      Logical, Intent(IN) :: L_Dead
      Integer, Intent(OUT) :: Nout
      Logical, Intent(INOUT) :: NM_STRS_MESH
      Real(EB), Intent(INOUT) :: hr_tau
      Real(EB), Intent(OUT) :: UBAR, VBAR, Tpre
      !
      ! Local variables
      Integer :: N,MN_now, STRS_Indx, door_ior, KKZ, j, i1, j1
      Real(EB) :: x_target, y_target, door_width, Door_Dist, EVEL
      Logical :: NM_STRS_MESHS, straight_Line_To_Target
      Type (MESH_TYPE), Pointer :: MFF=>NULL()
      Type (HUMAN_TYPE), Pointer :: HR=>NULL()
      Type (EVAC_STRS_Type), Pointer :: STRP=>NULL()

      HR=>HUMAN(I)
      KKZ = 1
      ! HR%I_FFIELD is the flow field mesh index for this agent
      NM_now = Max(1,HR%I_FFIELD)
      If (.Not.EVACUATION_ONLY(NM_now)) Then
         write(lu_err,*)'*** ',i,HR%IPC,HR%IEL,HR%I_Target,HR%I_FFIELD,HR%ILABEL
         Write(MESSAGE,'(A,I6,A)') 'ERROR: Evacuate_Humans, mesh ', NM_now, ' is not an evacuation flow field.'
         Call SHUTDOWN(MESSAGE)
      End If
      ! 
      ! Determine if the mesh is a stairs-mesh
      NM_STRS_MESH = .False.
      StrsMeshLoop: Do N = 1, N_STRS
         If (EVAC_STRS(N)%IMESH==NM_now) Then     
            NM_STRS_MESH = .True.
            STRS_Indx = N
            STRP=>EVAC_STRS(N)
            Nout = N
            Exit StrsMeshLoop
         End If
      End Do StrsMeshLoop

      MFF=>MESHES(NM_now)  ! Pointer to the flow field mesh

      ! Check if going straight line to target
      Straight_Line_To_Target = .False.
      If (NM_STRS_MESH) Then
         If (HR%I_Target == 0) Then
            !               write(*,*) 'Finding target within move loop'
            Call Find_Target_Node_In_Strs(STRP,HR)
         End If
         N = HR%I_Target
         If (N>N_DOORS) Then
            N = N - N_DOORS
            If (EVAC_EXITS(N)%STR_SUB_INDX == HR%STR_SUB_INDX) Then
               Straight_Line_To_Target = .TRUE.
               x_target = (EVAC_EXITS(N)%X1 + EVAC_EXITS(N)%X2)/2.0_EB
               y_target = (EVAC_EXITS(N)%Y1 + EVAC_EXITS(N)%Y2)/2.0_EB
               door_width = EVAC_EXITS(N)%Width
               door_ior   = EVAC_EXITS(N)%IOR
            End If
         Else If (N>0) Then
            If (EVAC_DOORS(N)%STR_SUB_INDX == HR%STR_SUB_INDX) Then
               Straight_Line_To_Target = .TRUE.
               x_target = (EVAC_DOORS(N)%X1 + EVAC_DOORS(N)%X2)/2.0_EB
               y_target = (EVAC_DOORS(N)%Y1 + EVAC_DOORS(N)%Y2)/2.0_EB
               door_width = EVAC_DOORS(N)%Width
               door_ior   = EVAC_DOORS(N)%IOR
            End If
         End If
      End If
      If (Straight_Line_To_Target) Then
         UBAR = x_target-HR%X
         VBAR = y_target-HR%Y
         Door_Dist = SQRT((x_target-HR%X)**2+(y_target-HR%Y)**2)
         If ( Door_Dist < 0.5_EB*Door_width ) Then
            Select Case(door_ior)
            Case(-1,+1)
               UBAR = SIGN(1.0_EB,UBAR)
               VBAR = 0._EB
               HR%SKIP_WALL_FORCE_IOR = NINT(UBAR)
            Case(-2,+2)
               UBAR = 0._EB
               VBAR = SIGN(1.0_EB,VBAR)
               HR%SKIP_WALL_FORCE_IOR = NINT(VBAR)*2
            End Select
         End If
      Else If (NM_STRS_MESH) Then 
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
      Else
         ! For thin walls and large dx,dy one do not get nice
         ! interpolation by using AFILL2. AFILL2 does not notice that
         ! there are thin walls and, thus, sometimes takes values from
         ! the other side of the thin wall in order to interpolate.
         ! UBAR = AFILL2(MFF%U,II-1,JJY,KKZ,XI-II+1,YJ-JJY+.5_EB,ZK-KKZ+.5_EB)
         ! VBAR = AFILL2(MFF%V,IIX,JJ-1,KKZ,XI-IIX+.5_EB,YJ-JJ+1,ZK-KKZ+.5_EB)
         UBAR = (1.0_EB-(XI-II+1))*MFF%U(II-1,JJ,1) + (XI-II+1)*MFF%U(II,JJ,1)
         VBAR = (1.0_EB-(YJ-JJ+1))*MFF%V(II,JJ-1,1) + (YJ-JJ+1)*MFF%V(II,JJ,1)
      End If

      EVEL = Sqrt(UBAR**2 + VBAR**2)  ! (UBAR,VBAR) is an unit vector
      If (EVEL > 0.0_EB) Then
         UBAR = UBAR/EVEL
         VBAR = VBAR/EVEL
      Else
         ! No v0 found for the current location of the agent, use previous value
         UBAR = HR%UBAR
         VBAR = HR%VBAR
      End If
      If (L_Dead) Then
         UBAR = 0.0_EB
         VBAR = 0.0_EB
      End If

      j = Max(0,HR%GROUP_ID)
      If (j == 0 ) Then
         Group_List(0)%Tpre = HR%Tpre + HR%Tdet
         Tpre = HR%Tpre + HR%Tdet
      Else
         Tpre = HR%Tdet
      End If
      If (L_Dead) Tpre = Huge(Tpre)

      ! Direction to the centre of the group (if any)
      If (Group_List(j)%GROUP_SIZE >= 2) Then
         HR%UBAR_Center = (Group_List(j)%GROUP_X - HR%X)
         HR%VBAR_Center = (Group_List(j)%GROUP_Y - HR%Y)
         EVEL = Sqrt(HR%UBAR_Center**2 + HR%VBAR_Center**2)
         If ( EVEL > 0.0_EB .And. .Not. L_Dead ) Then
            HR%UBAR_Center = HR%UBAR_Center / EVEL
            HR%VBAR_Center = HR%VBAR_Center / EVEL
         Else
            HR%UBAR_Center = 0.0_EB
            HR%VBAR_Center = 0.0_EB
         End If
      Else
         HR%UBAR_Center = 0.0_EB ! only one person in the group
         HR%VBAR_Center = 0.0_EB ! only one person in the group
      End If

      If ( j > 0 ) Then
         ! Group is already gathered together, but not yet moving towards the door
         If (Group_List(j)%COMPLETE == 1 .And. T <= Group_List(j)%Tpre+Group_List(j)%Tdoor) Then
            UBAR = 0.0_EB
            VBAR = 0.0_EB
         End If

         EVEL = UBAR**2 + VBAR**2
         If (Group_List(j)%COMPLETE == 1 .And. EVEL > 0.0_EB) Then
            ! The group is already gathered together
            UBAR = (1-GROUP_EFF)*UBAR + GROUP_EFF*HR%UBAR_Center/ &
                 Sqrt(((1-GROUP_EFF)*UBAR + GROUP_EFF*HR%UBAR_Center)**2+ &
                 ((1-GROUP_EFF)*VBAR + GROUP_EFF*HR%VBAR_Center)**2)
            VBAR = (1-GROUP_EFF)*VBAR + GROUP_EFF*HR%VBAR_Center/ &
                 Sqrt(((1-GROUP_EFF)*UBAR + GROUP_EFF*HR%UBAR_Center)**2+ &
                 ((1-GROUP_EFF)*VBAR + GROUP_EFF*HR%VBAR_Center)**2)
         Else
            ! The group is still in the gathering phase
            UBAR = HR%UBAR_Center
            VBAR = HR%VBAR_Center
         End If
      End If

      If ( T <= Tpre ) Then
         ! No movement yet
         UBAR = 0.0_EB
         VBAR = 0.0_EB
         hr_tau = Max(0.1_EB,HR%Tau/10.0_EB)
      End If
      If ( T <= T_BEGIN ) Then
         ! Initialization phase, i.e., flow field calculation
         UBAR = 0.0_EB
         VBAR = 0.0_EB
         hr_tau = HR%Tau
      End If
      Return
    End Subroutine Find_Prefered_Direction


    Subroutine GetStairSpeedAndZ(speed_xm,speed_xp, speed_ym,speed_yp,SP,HP)
      Implicit None
      !
      ! Passed variables
      Real(EB), Intent(OUT) :: speed_xm, speed_xp, speed_ym, speed_yp
      Type (EVAC_STRS_TYPE), Pointer ::  SP
      Type (HUMAN_TYPE), Pointer :: HP
      !
      ! Local variables
      Real(EB) cos_x, cos_y
      Integer J1, J2, J

      If (HP%STR_SUB_INDX > 0) Then
         J1 = Max(HP%STR_SUB_INDX-1,1)
         J2 = Min(HP%STR_SUB_INDX+1,STRP%N_NODES)
      Else
         J1 = 1
         J2 = STRP%N_NODES
      End If

      LoopSubIndx: Do J = J1,J2
         If (HP%X < SP%XB_NODE(J,1)) Cycle                
         If (HP%X > SP%XB_NODE(J,2)) Cycle                
         If (HP%Y < SP%XB_NODE(J,3)) Cycle                
         If (HP%Y > SP%XB_NODE(J,4)) Cycle                
         If (HP%Z < SP%XB_NODE(J,5)) Cycle                
         If (HP%Z > SP%XB_NODE(J,6)) Cycle
         ! Ok, human is inside the subnode J
         HP%STR_SUB_INDX = J      
         If (SP%NODE_TYPE(J)==STRS_STAIR_TYPE) Then
            cos_x = SP%XB_NODE(J,7)
            cos_y = SP%XB_NODE(J,8)
            Select Case (SP%NODE_IOR(J))
            Case(-1)
               HP%Z = SP%XB_NODE(J,5) + (SP%XB_NODE(J,6)-SP%XB_NODE(J,5))* &
                    Abs(SP%XB_NODE(J,2)-HP%X)/Abs(SP%XB_NODE(J,2)-SP%XB_NODE(J,1))
            Case(+1)
               HP%Z = SP%XB_NODE(J,5) + (SP%XB_NODE(J,6)-SP%XB_NODE(J,5))* &
                    Abs(SP%XB_NODE(J,1)-HP%X)/Abs(SP%XB_NODE(J,2)-SP%XB_NODE(J,1))
            Case(-2)
               HP%Z = SP%XB_NODE(J,5) + (SP%XB_NODE(J,6)-SP%XB_NODE(J,5))* &
                    Abs(SP%XB_NODE(J,4)-HP%Y)/Abs(SP%XB_NODE(J,4)-SP%XB_NODE(J,3))
            Case(+2)
               HP%Z = SP%XB_NODE(J,5) + (SP%XB_NODE(J,6)-SP%XB_NODE(J,5))* &
                    Abs(SP%XB_NODE(J,3)-HP%Y)/Abs(SP%XB_NODE(J,4)-SP%XB_NODE(J,3))
            End Select
            !            Select Case (SP%NODE_IOR(J)*HP%STRS_direction)
            Select Case (SP%NODE_IOR(J))
            Case(-1)
               speed_xm = cos_x* HP%Speed* SP%FAC_V0_UP
               speed_xp = cos_x* HP%Speed* SP%FAC_V0_DOWN
               speed_ym =        HP%Speed* SP%FAC_V0_HORI
               speed_yp =        HP%Speed* SP%FAC_V0_HORI
            Case(+1)
               speed_xm = cos_x* HP%Speed* SP%FAC_V0_DOWN
               speed_xp = cos_x* HP%Speed* SP%FAC_V0_UP
               speed_ym =        HP%Speed* SP%FAC_V0_HORI
               speed_yp =        HP%Speed* SP%FAC_V0_HORI
            Case(-2)
               speed_xm =        HP%Speed* SP%FAC_V0_HORI
               speed_xp =        HP%Speed* SP%FAC_V0_HORI
               speed_ym = cos_y* HP%Speed* SP%FAC_V0_UP
               speed_yp = cos_y* HP%Speed* SP%FAC_V0_DOWN
            Case(+2)
               speed_xm =        HP%Speed* SP%FAC_V0_HORI
               speed_xp =        HP%Speed* SP%FAC_V0_HORI
               speed_ym = cos_y* HP%Speed* SP%FAC_V0_DOWN
               speed_yp = cos_y* HP%Speed* SP%FAC_V0_UP
            End Select
         Else
            HP%Z = 0.5_EB*(SP%XB_NODE(J,5)+SP%XB_NODE(J,6))      
         End If
         Exit LoopSubIndx
      End Do LoopSubIndx
    End Subroutine GetStairSpeedAndZ

    Subroutine Find_Target_Node_In_Strs(SP,HP)
      Implicit None
      !
      ! Passed variables
      Type (EVAC_STRS_TYPE), Pointer :: SP
      Type (HUMAN_TYPE), Pointer :: HP
      !
      ! Local variables
      Logical IsKnownDoor,FinalTargetFound
      Integer :: I_Target = 0, I, Id, Final_node, IG, IN, Inode
      Real(EB) z_node, z_final, dz_node, dz_final, z_final_unknown,dz_tmp1, dz_tmp2, dz_node_actual
      Type (EVAC_ENTR_TYPE), Pointer :: PNX =>NULL()

      FinalTargetFound = .FALSE.
      IG = ABS(HP%GROUP_ID)

      z_final_unknown = 0._EB
      dz_node_actual = 0._EB

      Do I = 1, n_nodes
         IsKnownDoor = .FALSE.
         If (IG == 0 ) Then
            ! Human came from entry
            PNX => EVAC_ENTRYS(ABS(HR%IEL))
            If (PNX%N_VENT_FFIELDS == 0) Then
               IsKnownDoor = .True.
            Else
               Do Id = 1, PNX%N_VENT_FFIELDS
                  If (I == EVAC_Node_List(PNX%I_DOOR_NODES(II))%Node_Index) Then
                     IsKnownDoor = .True.
                     Exit
                  End If
               End Do
            Endif
         Else            
            ! Check the group know doors
            If (Human_Known_Doors(IG)%N_nodes == 0) Then
               IsKnownDoor = .True.
            Else
               IN = EVAC_Node_List(I)%Node_index
               Do Id = 1, Human_Known_Doors(IG)%N_nodes
                  If (I == Human_Known_Doors(IG)%I_nodes(Id)) Then
                     IsKnownDoor = .True.
                     Exit
                  End If
               End Do
            End If
         End If
         If (EVAC_NODE_List(I)%Node_type=='Exit') Then
            z_final_unknown = EVAC_EXITS(EVAC_NODE_List(I)%Node_index)%Z1
            If (IsKnownDoor) Then
               Final_node = I
               z_final = z_final_unknown
               FinalTargetFound = .TRUE.
               Exit
            End If
         End If
      End Do
      If (.NOT.FinalTargetFound) z_final = z_final_unknown
      dz_final = z_final - HP%Z

      ! Find the target among the nodes leading out of the STRS (doors and exits)
      FinalTargetFound = .FALSE.
      dz_tmp2 = Huge(1.0_EB)
      FindTargetNodeLoop: Do I = 1,SP%N_NODES_OUT
         Inode = SP%NODES_OUT(I)
         dz_node = -1._EB * SIGN(1.0_EB,dz_final) ! initialize dz_node to different direction than final target
         Select Case(EVAC_NODE_List(Inode)%Node_type)
         Case('Door')
            z_node = EVAC_DOORS(EVAC_Node_List(Inode)%Node_index)%Z1
         Case('Exit')
            z_node = EVAC_EXITS(EVAC_Node_List(Inode)%Node_index)%Z1
         Case default
            write(LU_ERR,*) 'ERROR (debug): unknown node type in Find_Target_Node_In_Strs'
         End Select
         dz_node = z_node - HP%Z
         If (SIGN(1.0_EB,dz_node)==SIGN(1.0_EB,dz_final)) Then
            dz_tmp1 = ABS(dz_final)-ABS(dz_node)
            If ( dz_tmp1 < dz_tmp2) Then
               FinalTargetFound = .TRUE.
               HP%I_Target = EVAC_Node_List(SP%NODES_OUT(I))%Node_index
               dz_tmp2 = dz_tmp1
               dz_node_actual = dz_node
            End If
         End If
      End Do FindTargetNodeLoop

      ! Find the target among the nodes leading in to the STRS, i.e. DOORS
      If (.NOT.FinalTargetFound) Then
      FindTargetNodeLoop2: Do I = 1,SP%N_NODES_IN
         Inode = SP%NODES_IN(I)
         dz_node = -1._EB * SIGN(1.0_EB,dz_final) ! initialize dz_node to different direction than final target
         Select Case(EVAC_NODE_List(Inode)%Node_type)
         Case('Door')
            z_node = EVAC_DOORS(EVAC_Node_List(Inode)%Node_index)%Z1
         Case('Exit')
!            z_node = EVAC_EXITS(EVAC_Node_List(SP%NODES_IN(I))%Node_index)%Z1
         Case('Entry')
!            z_node = EVAC_ENTRYS(EVAC_Node_list(SP%NODES_IN(I))%Node_index)%Z1
         Case default
            Write(LU_ERR,*) 'ERROR (debug): unknown node type in Find_Target_Node_In_Strs'
         End Select
         dz_node = z_node - HP%Z
         If (SIGN(1.0_EB,dz_node)==SIGN(1.0_EB,dz_final)) Then
            dz_tmp1 = ABS(dz_final-dz_node)
            If ( dz_tmp1 < dz_tmp2) Then
               FinalTargetFound = .TRUE.
               HP%I_Target = EVAC_Node_List(Inode)%Node_index
               dz_tmp2 = dz_tmp1
               dz_node_actual = dz_node
            End If
         End If
      End Do FindTargetNodeLoop2
      End If
      if (dz_node_actual >= 0._EB) Then
         HP%STRS_Direction = 1
      Else
         HP%STRS_Direction = -1   
      End If
    End Subroutine Find_Target_Node_In_Strs

    Subroutine STRS_U_AND_V(STRP,I_NODE,X,Y,Direction,UBAR,VBAR)
      Implicit None
      !
      ! Get preferred direction in STRS
      !
      ! Passed variables
      Type (EVAC_STRS_TYPE), Pointer::  STRP
      Integer, Intent(IN) :: I_NODE, Direction
      REAL(EB), Intent(IN) :: X, Y
      REAL(EB), Intent(OUT) :: UBAR, VBAR
      !
      ! Local variables
      Integer I_CORE

      I_CORE = STRP%I_CORE(I_NODE)

      UBAR = 0._EB
      VBAR = 0._EB

      If (STRP%RIGHT_HANDED) Then
         If (Direction > 0) Then
            IF ((X <= STRP%XB_CORE(I_CORE,2)) .AND. (Y <= STRP%XB_CORE(I_CORE,3))) UBAR =  1._EB
            IF ((X >= STRP%XB_CORE(I_CORE,2)) .AND. (Y <= STRP%XB_CORE(I_CORE,4))) VBAR =  1._EB
            IF ((X >= STRP%XB_CORE(I_CORE,1)) .AND. (Y >= STRP%XB_CORE(I_CORE,4))) UBAR = -1._EB
            IF ((X <= STRP%XB_CORE(I_CORE,1)) .AND. (Y >= STRP%XB_CORE(I_CORE,3))) VBAR = -1._EB
         Else
            IF ((X <= STRP%XB_CORE(I_CORE,1)) .AND. (Y <= STRP%XB_CORE(I_CORE,4))) VBAR =  1._EB
            IF ((X <= STRP%XB_CORE(I_CORE,2)) .AND. (Y >= STRP%XB_CORE(I_CORE,4))) UBAR =  1._EB
            IF ((X >= STRP%XB_CORE(I_CORE,2)) .AND. (Y >= STRP%XB_CORE(I_CORE,3))) VBAR = -1._EB
            IF ((X >= STRP%XB_CORE(I_CORE,1)) .AND. (Y <= STRP%XB_CORE(I_CORE,3))) UBAR = -1._EB
         Endif
      Else
         If (Direction < 0) Then
            IF ((X <= STRP%XB_CORE(I_CORE,2)) .AND. (Y <= STRP%XB_CORE(I_CORE,3))) UBAR =  1._EB
            IF ((X >= STRP%XB_CORE(I_CORE,2)) .AND. (Y <= STRP%XB_CORE(I_CORE,4))) VBAR =  1._EB
            IF ((X >= STRP%XB_CORE(I_CORE,1)) .AND. (Y >= STRP%XB_CORE(I_CORE,4))) UBAR = -1._EB
            IF ((X <= STRP%XB_CORE(I_CORE,1)) .AND. (Y >= STRP%XB_CORE(I_CORE,3))) VBAR = -1._EB
         Else
            IF ((X <= STRP%XB_CORE(I_CORE,1)) .AND. (Y <= STRP%XB_CORE(I_CORE,4))) VBAR =  1._EB
            IF ((X <= STRP%XB_CORE(I_CORE,2)) .AND. (Y >= STRP%XB_CORE(I_CORE,4))) UBAR =  1._EB
            IF ((X >= STRP%XB_CORE(I_CORE,2)) .AND. (Y >= STRP%XB_CORE(I_CORE,3))) VBAR = -1._EB
            IF ((X >= STRP%XB_CORE(I_CORE,1)) .AND. (Y <= STRP%XB_CORE(I_CORE,3))) UBAR = -1._EB
         Endif
      Endif
    End Subroutine STRS_U_AND_V

    Subroutine GET_IW(IIin,JJin,KKin,IOR,IW)
      Implicit None
      !
      ! Passed variables
      Integer, Intent(IN) :: IIin, JJin, KKin, IOR
      Integer, Intent(OUT) :: IW
      !
      ! Local variables
      Integer :: ii, jj, kk, ic, I_OBST
      !
      ii = IIin
      jj = JJin
      kk = KKin

      IC  = CELL_INDEX(II,JJ,KK)
      I_OBST = OBST_INDEX_C(IC)

      If (SOLID(IC) .AND. .NOT.OBSTRUCTION(I_OBST)%HIDDEN) Then
         Select Case(IOR)
         Case(-1)
            If (II>0)      II = II-1
         Case( 1)
            If (II<IBAR+1) II = II+1
         Case(-2)
            If (JJ>0)      JJ = JJ-1
         Case( 2)
            If (JJ<JBAR+1) JJ = JJ+1
         Case(-3)
            If (KK>0)      KK = KK-1
         Case( 3)
            If (KK<KBAR+1) KK = KK+1
         End Select
      End If

      IC  = CELL_INDEX(II,JJ,KK)
      IW  = WALL_INDEX(IC,-IOR)

      If (IW<=0) Then
         Select Case(IOR)
         Case(-1)
            If (II>0)      IC = CELL_INDEX(II-1,JJ,KK)
         Case( 1)
            If (II<IBAR+1) IC = CELL_INDEX(II+1,JJ,KK)
         Case(-2)
            If (JJ>0)      IC = CELL_INDEX(II,JJ-1,KK)
         Case( 2)
            If (JJ<JBAR+1) IC = CELL_INDEX(II,JJ+1,KK)
         Case(-3)
            If (KK>0)      IC = CELL_INDEX(II,JJ,KK-1)
         Case( 3)
            If (KK<KBAR+1) IC = CELL_INDEX(II,JJ,KK+1)
         End Select
         IW = WALL_INDEX(IC,-IOR)
      End If

    End Subroutine GET_IW

    !
    Subroutine CHECK_EXITS(T,NM)
      Implicit None
      !
      ! Remove persons if they are found at an exit.
      !
      ! Passed variables
      Real(EB), Intent(IN) :: T
      Integer, Intent(IN) :: NM
      !
      ! Local variables
      Real(EB) x_old, y_old, pexx1, pexx2, pexy1, pexy2
      Integer :: ie,i,n_tmp
      Type (EVAC_EXIT_Type), Pointer :: PEX =>NULL()
      Type (HUMAN_TYPE), Pointer :: HR =>NULL()
      !
      HUMAN(:)%IOR = 0
      PexLoop: Do ie = 1, n_exits
         PEX=>EVAC_EXITS(ie)
         If (PEX%IMESH /= NM ) Cycle PexLoop
         Select Case (PEX%IOR)
         Case (-1,+1)
            pexx1 = PEX%X1
            pexx2 = PEX%X2
            pexy1 = 0.5_EB*(PEX%Y1+PEX%Y2)-0.5_EB*PEX%Width
            pexy2 = 0.5_EB*(PEX%Y1+PEX%Y2)+0.5_EB*PEX%Width
         Case (-2,+2)
            pexx1 = 0.5_EB*(PEX%X1+PEX%X2)-0.5_EB*PEX%Width
            pexx2 = 0.5_EB*(PEX%X1+PEX%X2)+0.5_EB*PEX%Width
            pexy1 = PEX%Y1
            pexy2 = PEX%Y2
         End Select
         HumLoop: Do i = 1, n_humans
            HR=>HUMAN(I)
            If ( HR%IOR > 0 ) Cycle HumLoop
            x_old = HR%X_old
            y_old = HR%Y_old
            Select Case (PEX%IOR)
            Case (+1)
               If ((HR%X >= pex%x1 .And. x_old < pex%x1) .And. (HR%Y >= pex%y1 .And. HR%Y <= pex%y2) ) Then
                  If (PEX%COUNT_ONLY) HR%IOR = 2
                  If (.Not. PEX%COUNT_ONLY) HR%IOR = 1
                  PEX%T_last=T
                  PEX%ICOUNT = PEX%ICOUNT + 1
                  If (PEX%T_first <= T_BEGIN) PEX%T_first = T
               End If
            Case (-1)
               If ((HR%X <= pex%x2 .And. x_old > pex%x2) .And. (HR%Y >= pex%y1 .And. HR%Y <= pex%y2) ) Then
                  If (PEX%COUNT_ONLY) HR%IOR = 2
                  If (.Not. PEX%COUNT_ONLY) HR%IOR = 1
                  PEX%T_last=T
                  PEX%ICOUNT = PEX%ICOUNT + 1
                  If (PEX%T_first <= T_BEGIN) PEX%T_first = T
               End If
            Case (+2)
               If ((HR%Y >= pex%y1 .And. y_old < pex%y1) .And. (HR%X >= pex%x1 .And. HR%X <= pex%x2) ) Then
                  If (PEX%COUNT_ONLY) HR%IOR = 2
                  If (.Not. PEX%COUNT_ONLY) HR%IOR = 1
                  PEX%T_last=T
                  PEX%ICOUNT = PEX%ICOUNT + 1
                  If (PEX%T_first <= T_BEGIN) PEX%T_first = T
               End If
            Case (-2)
               If ((HR%Y <= pex%y2 .And. y_old > pex%y2) .And. (HR%X >= pex%x1 .And. HR%X <= pex%x2) ) Then
                  If (PEX%COUNT_ONLY) HR%IOR = 2
                  If (.Not. PEX%COUNT_ONLY) HR%IOR = 1
                  PEX%T_last=T
                  PEX%ICOUNT = PEX%ICOUNT + 1
                  If (PEX%T_first <= T_BEGIN) PEX%T_first = T
               End If
            End Select
            If (HR%IOR > 0) Then
               If (.Not. PEX%COUNT_ONLY) Write (LU_EVACOUT,fmt='(a,i6,a,f8.2,a,a,a,f8.4,a,i4)') &
                    ' Agent n:o', HR%ILABEL, ' out at ', T, ' s, exit ', Trim(PEX%ID), &
                    ', FED=', HR%IntDose, ', color_i=', HR%COLOR_INDEX
               If (      PEX%COUNT_ONLY) Write (LU_EVACOUT,fmt='(a,i6,a,f8.2,a,a,a,f8.4,a,i4)') &
                    ' Agent n:o', HR%ILABEL, ' counted at ', T, ' s, counter ', Trim(PEX%ID), &
                    ', FED=', HR%IntDose, ', color_i=', HR%COLOR_INDEX
               If (HR%IOR == 2) HR%IOR = HUMAN_NO_TARGET
            End If
         End Do HumLoop
      End Do PexLoop
      n_tmp = n_humans
      Remove_loop: Do i = n_tmp, 1, -1
         HR=>HUMAN(I)
         ! Remove person, if wanted.
         If (HR%IOR >= 1) Then
            HR%IOR = HUMAN_NO_TARGET
            Call Remove_Person(I)
         End If
      End Do Remove_loop
      !
    End Subroutine CHECK_EXITS
    !
    Subroutine CHECK_DOORS(T,NM)
      Implicit None
      !
      ! Replace persons if they are found at a door. 
      !
      ! Passed variables
      Real(EB), Intent(IN) :: T
      Integer, Intent(IN) :: NM
      !
      ! Local variables
      Real(EB) x_old, y_old, xx, yy, zz, pdxx1, pdxx2, pdxy1, pdxy2, v, angle
      Integer :: ie,i,n_tmp, istat, ior_new, inode2, imesh2, n, ior
      Integer :: new_ffield_i, color_index, i_target, inode, STR_INDX, STR_SUB_INDX
      Character(60) :: TO_NODE
      Character(26) :: new_ffield_name
      Logical :: keep_xy, upstream
      Type (EVAC_DOOR_Type), Pointer :: PDX =>NULL()
      Type (HUMAN_TYPE), Pointer :: HR =>NULL()
      !
      keep_xy = .False.
      HUMAN(:)%IOR = HUMAN_NO_TARGET
      PdxLoop: Do ie = 1, n_doors
         PDX=>EVAC_DOORS(ie)
         ! Note: IMESH2 is not good for CORR targets
         If (PDX%IMESH /= NM .AND. PDX%IMESH2 /= NM) Cycle PdxLoop
         If (PDX%IMESH /= NM .And. .Not.NM_STRS_MESH) Cycle PdxLoop
         keep_xy = PDX%KEEP_XY
         Select Case (PDX%IOR)
         Case (-1,+1)
            pdxx1 = PDX%X1
            pdxx2 = PDX%X2
            pdxy1 = 0.5_EB*(PDX%Y1+PDX%Y2)-0.5_EB*PDX%Width
            pdxy2 = 0.5_EB*(PDX%Y1+PDX%Y2)+0.5_EB*PDX%Width
         Case (-2,+2)
            pdxx1 = 0.5_EB*(PDX%X1+PDX%X2)-0.5_EB*PDX%Width
            pdxx2 = 0.5_EB*(PDX%X1+PDX%X2)+0.5_EB*PDX%Width
            pdxy1 = PDX%Y1
            pdxy2 = PDX%Y2
         End Select
         HumLoop: Do i = 1, n_humans
            HR=>HUMAN(I)
            new_ffield_name = Trim(HR%FFIELD_NAME)
            new_ffield_i = HR%I_FFIELD
            ! Check only one door... (Yes, this can happen...)
            If ( HR%IOR /= HUMAN_NO_TARGET) Cycle HumLoop
            ! Coming out upstream?
            upstream = .FALSE.
            If (PDX%IMESH2 == NM .And. NM_STRS_MESH) Then
               ! here, check for correct height
               If ((HR%Z > PDX%Z2) .OR. (HR%Z < PDX%Z1)) Cycle HumLoop
               upstream = .TRUE.
            End If
            x_old = HR%X_old
            y_old = HR%Y_old
            Select Case (ABS(PDX%IOR))
            Case (1)
               If (SIGN(1.0_EB,HR%X-pdx%x1)/=SIGN(1.0_EB,x_old-pdx%x1) .And. (HR%Y >= pdx%y1 .And. HR%Y <= pdx%y2) ) Then
                  HR%IOR = HUMAN_TARGET_UNSPECIFIED
               End If
            Case (2)
               If (SIGN(1.0_EB,HR%Y-pdx%y1)/=SIGN(1.0_EB,y_old-pdx%y1) .And. (HR%X >= pdx%x1 .And. HR%X <= pdx%x2) ) Then
                  HR%IOR = HUMAN_TARGET_UNSPECIFIED
               End If
            End Select
            If ( HR%IOR == HUMAN_TARGET_UNSPECIFIED ) Then
               istat = 0
               If (upstream) Then
                  inode  = PDX%INODE2
                  inode2 = PDX%INODE
               Else
                  inode  = PDX%INODE
                  inode2 = PDX%INODE2
               End If
               Call Check_Target_Node(inode,inode2,istat,xx,yy,zz,ior_new, &
                    imesh2,T,new_ffield_name,new_ffield_i,color_index,&
                    i_target, keep_xy, angle, STR_INDX, STR_SUB_INDX, HR)
               If (istat == 0 ) Then
                  ! Put person to a new node, i.e., target node has empty space
                  HR%X = xx
                  HR%Y = yy
                  HR%Z = zz
                  HR%Angle = angle
                  If (STR_INDX>0) HR%STR_SUB_INDX = STR_SUB_INDX
                  If (keep_xy .And. ior_new == HUMAN_TARGET_UNSPECIFIED) Then
                     If (EVAC_Node_List(INODE2)%Node_Type == 'Door') Then
                        xx = 0.5_EB*(EVAC_DOORS(EVAC_Node_List(INODE2)%Node_Index)%X1 + &
                             EVAC_DOORS(EVAC_Node_List(INODE2)%Node_Index)%X2 - (PDX%X1+PDX%X2))
                        yy = 0.5_EB*(EVAC_DOORS(EVAC_Node_List(INODE2)%Node_Index)%Y1 + &
                             EVAC_DOORS(EVAC_Node_List(INODE2)%Node_Index)%Y2 - (PDX%Y1+PDX%Y2))
                     Else
                        xx = 0.5_EB*(EVAC_ENTRYS(EVAC_Node_List(INODE2)%Node_Index)%X1 + &
                             EVAC_ENTRYS(EVAC_Node_List(INODE2)%Node_Index)%X2 - (PDX%X1+PDX%X2))
                        yy = 0.5_EB*(EVAC_ENTRYS(EVAC_Node_List(INODE2)%Node_Index)%Y1 + &
                             EVAC_ENTRYS(EVAC_Node_List(INODE2)%Node_Index)%Y2 - (PDX%Y1+PDX%Y2))
                     End If
                     If (Abs(PDX%ior) == 2) HR%X = HR%X + xx  
                     If (Abs(PDX%ior) == 1) HR%Y = HR%Y + yy 
                  End If
                  v = Sqrt( HR%U**2 + HR%V**2 )

                  HR%IOR = ior_new
                  If (ior_new == HUMAN_TARGET_UNSPECIFIED .Or. ior_new == HUMAN_STRS_TARGET) Then
                     ! door/entry or STRS-mesh
                     If (HR%IMESH /= imesh2 ) Then
                        ! Put the person to a new floor
                        HR%IOR         = HUMAN_ANOTHER_MESH_TARGET
                        HR%INODE = 0
                        Do n = 1, n_egrids
                           If (EVAC_Node_List(n)%IMESH == imesh2) HR%INODE = n
                        End Do
                        HR%NODE_NAME   = Trim(MESH_NAME(imesh2))
                        HR%FFIELD_NAME = Trim(new_ffield_name)
                        HR%I_FFIELD = new_ffield_i
                        HR%SKIP_WALL_FORCE_IOR = 0
                     Else
                        HR%IOR = -2   ! same floor (door/entry)
                        HR%FFIELD_NAME = Trim(new_ffield_name)
                        HR%I_FFIELD = new_ffield_i
                     End If
                     HR%I_Target = I_Target
                     HR%X_old = HR%X
                     HR%Y_old = HR%Y
                     ! ior is the direction where the human is ejected.
                     If (EVAC_Node_List(INODE2)%Node_Type == 'Door') Then
                        ior = - EVAC_DOORS(EVAC_Node_List(INODE2)%Node_Index)%IOR
                     End If
                     If (EVAC_Node_List(INODE2)%Node_Type == 'Entry') Then
                        ior = EVAC_ENTRYS(EVAC_Node_List(INODE2)%Node_Index)%IOR
                     End If
                     If (EVAC_Node_List(INODE2)%Node_Type == 'Floor') Then
                        ! For STRS the first node should be DOOR
                        ior = EVAC_DOORS(EVAC_Node_List(INODE)%Node_Index)%IOR
                     End If
                     If ( Abs(ior) == 1 ) Then
                        HR%U = v*ior
                        HR%V = 0.0_EB
                        If (.Not.keep_xy) Then
                           ! 180 or 0 degrees, i.e., pi or 0 radians
                           HR%Angle = (0.5_EB-(ior/2.0_EB))*Pi
                        End If
                     End If
                     If ( Abs(ior) == 2 ) Then
                        HR%U = 0.0_EB
                        HR%V = 0.5_EB*v*ior
                        If (.Not.keep_xy) Then
                           ! 270 or 90 degrees, i.e., 3pi/2 or pi/2 radians
                           HR%Angle = (1.0_EB-(ior/4.0_EB))*Pi
                        End If
                     End If

                     If (HR%IMESH /= imesh2 ) Then
                        HR%IMESH       = imesh2
                        If (MESHES(imesh2)%N_HUMANS+1 > MESHES(imesh2)%N_HUMANS_DIM) Then
                           ! Re-allocation is not yet checked.
                           Call SHUTDOWN('ERROR: Humans: no re-allocation yet')
                           Call RE_ALLOCATE_HUMANS(1,imesh2)
                        End If
                        MESHES(imesh2)%N_HUMANS = MESHES(imesh2)%N_HUMANS + 1
                        MESHES(imesh2)%HUMAN(MESHES(imesh2)%N_HUMANS) = HUMAN(I)
                        MESHES(imesh2)%HUMAN(MESHES(imesh2)%N_HUMANS)%IOR = 0
                     End If

                  End If            ! target is door or entry

                  PDX%T_last=T
                  PDX%ICOUNT = PDX%ICOUNT + 1
                  If (PDX%T_first <= T_BEGIN) PDX%T_first = T

                  Write (LU_EVACOUT,fmt='(a,i6,a,f8.2,a,a,a,f8.4)') ' EVAC: Person n:o', HR%ILABEL, ' out at ', T, &
                       ' s, door ', Trim(PDX%ID), ' fed ', HR%IntDose

               Else    ! istat = 1 ==> do not move to node
                  ! Can not move to the next node, so do not allow to move inside
                  ! the door ==> keep old position and put velocity to zero.
                  HR%X = HR%X_old
                  HR%Y = HR%Y_old
                  HR%U = 0.0_EB
                  HR%V = 0.0_EB
                  HR%IOR = -1  ! can not move to target node
               End If
            End If
         End Do HumLoop
      End Do PdxLoop

      ! Remove humans from this floor and put them some other floor
      n_tmp = n_humans
      Remove_loop: Do i = n_tmp, 1, -1
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
         If (HR%IOR >= 1) Then
            HR%IOR = HUMAN_NO_TARGET
            Call Remove_Person(I)
         End If
      End Do Remove_loop
      !
    End Subroutine CHECK_DOORS
    !
    Subroutine CHECK_CORRS(T,NM,DTSP)
      Implicit None
      !
      ! Entry persons from the corridors.
      !
      ! Passed variables
      Real(EB), Intent(IN) :: T,DTSP
      Integer, Intent(IN) :: NM
      !
      ! Local variables
      Real(EB) x_old, y_old, xx, yy, zz, pcxx1, pcxx2, pcxy1, pcxy2, v, x_int, angle
      Integer :: ie,i,n_tmp, istat, ior_new, inode2, imesh2, n, ior
      Integer :: new_ffield_i, color_index, i_target, inode, STR_INDX, STR_SUB_INDX
      Character(60) :: TO_NODE
      Character(26) :: new_ffield_name
      Logical :: keep_xy
      Type (EVAC_CORR_Type),  Pointer :: PCX=>NULL()
      Type (CORR_LL_TYPE), Pointer :: Now_LL=>NULL(), Tmp_LL=>NULL(), Next_LL=>NULL(), Prev_LL=>NULL()
      Type (HUMAN_TYPE), Pointer :: HR =>NULL()
      !
      keep_xy = .False.
      PcxLoop: Do ie = 1, n_corrs
         PCX=>EVAC_CORRS(ie)
         If (PCX%IMESH2 /= NM) Cycle PcxLoop
         If (PCX%n_inside <= 0) Cycle PcxLoop
         Nullify(Prev_LL)
         Nullify(Now_LL)
         Nullify(Next_LL)
         If (Associated(PCX%First)) Then
            Now_LL => PCX%First
            If (Associated(PCX%First%Next)) Next_LL => PCX%First%Next
         End If
         HumLoop: Do
            If (.Not. Associated(Now_LL)) Exit HumLoop
            istat = 0
            inode = PCX%INODE
            inode2 = PCX%INODE2
            HR => Now_LL%HUMAN
            If ( HR%IntDose >= 1.0_EB  ) Then
               If (HR%Tpre /= Huge(HR%Tpre)) Then
                  n_dead = n_dead+1
                  HR%Tpre = Huge(HR%Tpre)
                  HR%Tdet = Huge(HR%Tdet)
                  HR%Tau  = HR%Tau
                  HR%Mass = HR%Mass
               End If
            Else
               fed_max_alive = Max(fed_max_alive,HR%IntDose)
            End If
            fed_max = Max(fed_max,HR%IntDose)

            ! Calculate Purser's fractional effective dose (FED)
            If (T > T_BEGIN) Then
               If ( PCX%FED_MESH2 > 0 ) Then
                  x_int = Min(1.0_EB,Max(0.0_EB,(Now_LL%T_out-T)) / (Now_LL%T_out - Now_LL%T_in))
               Else
                  x_int = 1.0_EB
               End If
               HR%IntDose = DTSP*( (1.0_EB-x_int)*PCX%FED_CO_CO2_O2(2) + x_int*PCX%FED_CO_CO2_O2(1) ) + HR%IntDose
            End If

            If ( (Now_LL%T_out) <= T) Then
               Call Check_Target_Node(inode,inode2,istat,xx,yy,zz,ior_new, &
                    imesh2,T,new_ffield_name,new_ffield_i,color_index, &
                    i_target,keep_xy,angle,STR_INDX,STR_SUB_INDX, HR)

               If (istat == 0 ) Then
                  ! Put person to a new node, i.e., target node has empty space
                  HR%X = xx
                  HR%Y = yy
                  HR%Z = zz
                  HR%Angle = angle

                  v = Sqrt( HR%U**2 + HR%V**2 )
                  If (v > 0.01_EB) v = PCX%Fac_Speed*HR%v0_fac*HR%Speed
                  HR%IOR = 1        ! remove from corridor
                  ! Put the person to a floor
                  If (ior_new == 1) Then ! door or entry
                     HR%IMESH = imesh2
                     HR%INODE = 0
                     Do n = 1, n_egrids
                        If (EVAC_Node_List(n)%IMESH == imesh2) HR%INODE = n
                     End Do
                     HR%NODE_NAME   = Trim(MESH_NAME(imesh2))
                     HR%FFIELD_NAME = Trim(new_ffield_name)
                     HR%I_FFIELD = new_ffield_i
                     HR%I_Target = I_Target

                     If (EVAC_Node_List(INODE2)%Node_Type == 'Door') Then
                        ior = - EVAC_DOORS( EVAC_Node_List(INODE2)%Node_Index)%IOR
                     End If
                     If (EVAC_Node_List(INODE2)%Node_Type == 'Entry') Then
                        ior = EVAC_ENTRYS( EVAC_Node_List(INODE2)%Node_Index)%IOR
                     End If
                     If ( Abs(ior) == 1 ) Then
                        HR%U = v*ior
                        HR%V = 0.0_EB
                     End If
                     If ( Abs(ior) == 2 ) Then
                        HR%U = 0.0_EB
                        HR%V = 0.5_EB*v*ior
                     End If
                     HR%X_old = HR%X
                     HR%Y_old = HR%Y

                     If (MESHES(imesh2)%N_HUMANS+1 > MESHES(imesh2)%N_HUMANS_DIM) Then
                        ! Re-allocation is not yet checked.
                        Call SHUTDOWN('ERROR: Humans: no re-allocation yet')
                        Call RE_ALLOCATE_HUMANS(1,imesh2)
                     End If
                     MESHES(imesh2)%N_HUMANS = MESHES(imesh2)%N_HUMANS + 1
                     MESHES(imesh2)%HUMAN(MESHES(imesh2)%N_HUMANS) = HR
                     Write (LU_EVACOUT,fmt='(a,i6,a,f8.2,a,a,a,f8.4)') ' EVAC: Person n:o', &
                          HR%ILABEL, ' out at ', T, ' s, corr ', Trim(PCX%ID), ' fed ', HR%IntDose
                  End If            ! target is door or entry

                  If (ior_new == 3) Then ! corridor
                     Write (LU_EVACOUT,fmt='(a,i6,a,f8.2,a,a,f8.4)') ' EVAC: Person n:o', &
                          HR%ILABEL, ' change corr ', T, ' s, corr ', Trim(PCX%ID), HR%IntDose
                  End If

                  If (ior_new == 5) Then ! exit
                     Write (LU_EVACOUT,fmt='(a,i6,a,f8.2,a,a,f8.4)') ' EVAC: Person n:o', &
                          HR%ILABEL, ' exits ', T, ' s, corr ', Trim(PCX%ID), HR%IntDose
                  End If
               Else
                  ! Can not move to the next node, so do not allow to move inside
                  ! the door ==> keep old position and put velocity to zero.
                  HR%X = HR%X_old
                  HR%Y = HR%Y_old
                  HR%I_Target = 0
                  HR%U = 0.0_EB
                  HR%V = 0.0_EB
                  HR%IOR = HUMAN_NO_TARGET
               End If         ! istat == 0

            Else             ! T_out > T, i.e., still in corridor
               HR%X = HR%X_old
               HR%Y = HR%Y_old
               HR%IOR = HUMAN_NO_TARGET
               HR%I_Target = 0
            End If

            ! Here person is removed from the corridor (linked list)
            If (HR%IOR > 0) Then
               PCX%n_inside = PCX%n_inside - 1
               If (Associated(Now_LL%Next)) Then
                  Tmp_LL => Now_LL%Next 
                  Nullify(HR)       ! HR points to Now_LL%Human
                  Nullify(Now_LL%Next) ! remove the pointer
                  Deallocate(Now_LL) ! free memory
                  Nullify(Now_LL)   ! remove the pointer
                  Now_LL => Tmp_LL  ! Jump to the next element in LL
                  If (.Not. Associated(Prev_LL)) Then
                     PCX%First => Now_LL   ! first element removed
                  Else
                     Prev_LL%Next => Now_LL ! deleted element skipped
                  End If
               Else                ! remove the last element
                  Nullify(HR)       ! HR points to Now_LL%Human
                  Nullify(Now_LL%Next) ! end list
                  Deallocate(Now_LL) ! free memory
                  Nullify(Now_LL)   ! end loop
                  If (Associated(Prev_LL)) Then
                     Nullify(Prev_LL%Next) ! end the list
                  Else
                     Nullify(PCX%First)   ! empty list
                  End If
               End If
            Else
               Prev_LL => Now_LL  ! save the previous element
               If (Associated(Now_LL%Next)) Then
                  Now_LL => Now_LL%Next ! Jump to the next element in LL
               Else
                  Nullify(Now_LL)  ! end the loop
               End If
            End If
         End Do HumLoop
      End Do PcxLoop
      !
    End Subroutine CHECK_CORRS
    !
    Subroutine Check_Target_Node(inode,inode2,istat,xx,yy,zz,ior_new, &
         imesh2,T,new_ffield_name,new_ffield_i,color_index,i_target,keep_xy,angle,&
         STR_INDX, STR_SUB_INDX, HR)
      Implicit None
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
      Integer, Intent(IN) :: inode, inode2
      Integer, Intent(OUT) :: istat, ior_new, imesh2, color_index, i_target
      Real(EB), Intent(OUT) :: xx, yy, zz, angle
      Integer, Intent(OUT) :: STR_INDX, STR_SUB_INDX
      Real(EB), Intent(IN) :: T
      Integer, Intent(INOUT) :: new_ffield_i
      Logical, Intent(IN) :: keep_xy
      Character(26), Intent(INOUT) :: new_ffield_name
      Type (HUMAN_TYPE), Pointer :: HR
      !
      ! Local variables
      Real(EB) RN, x1, x2, y1, y2, z1, z2, d_max, dist, Width, &
           xx1,yy1, max_fed, ave_K
      Integer  II, JJ, KK, ior, irnmax, irn, ie, izero, j1
      Real(EB), Dimension(6) :: r_tmp, x_tmp, y_tmp
      Integer :: i_tmp, i_tim, iii, jjj, I_OBST
      Logical :: PP_see_door, keep_xy2, NM_STRS_MESH

      Type (CORR_LL_TYPE), Pointer :: TmpCurrent =>NULL(), TmpLoop =>NULL()
      Type (EVAC_STRS_TYPE), Pointer :: STRP =>NULL()
      Type (EVAC_DOOR_Type), Pointer :: PDX2 =>NULL()
      Type (MESH_TYPE), Pointer :: MMF =>NULL()
      Type (EVACUATION_Type), Pointer :: HPT =>NULL()
      Type (HUMAN_TYPE), Pointer :: HRE =>NULL()
      Type (EVAC_ENTR_Type), Pointer :: PNX =>NULL(), PNX2 =>NULL()
      Type (EVAC_CORR_Type), Pointer :: PCX2 =>NULL()
      !
      xx = 0.0_EB ; yy = 0.0_EB ; zz = 0.0_EB 
      I_Target = 0
      STR_INDX = 0
      STR_SUB_INDX = 0

      ! Where is the person going to?
      SelectTargetType: Select Case (EVAC_Node_List(INODE2)%Node_Type)
      Case ('Door','Entry') SelectTargetType
         ior_new = 1
         If (EVAC_Node_List(INODE2)%Node_Type == 'Door') Then
            PDX2 => EVAC_DOORS(EVAC_Node_List(INODE2)%Node_Index)
            imesh2  = PDX2%IMESH
            new_ffield_name = Trim(PDX2%GRID_NAME)
            Irn_Loop1: Do irn = 1, nmeshes
               If ( evacuation_only(irn) .And. Trim(PDX2%GRID_NAME) == Trim(MESH_NAME(irn)) ) Then
                  new_ffield_i = irn
                  Exit Irn_Loop1
               End If
            End Do Irn_Loop1
            X1  = PDX2%X1
            X2  = PDX2%X2
            Y1  = PDX2%Y1
            Y2  = PDX2%Y2
            Z1  = PDX2%Z1
            Z2  = PDX2%Z2
            ior = -PDX2%IOR        ! now entering the room
            Width = PDX2%Width
         Else
            PNX2 => EVAC_ENTRYS(EVAC_Node_List(INODE2)%Node_Index) 
            imesh2  = PNX2%IMESH
            new_ffield_name = Trim(PNX2%GRID_NAME)
            Irn_Loop2: Do irn = 1, nmeshes
               If ( evacuation_only(irn) .And. Trim(PNX2%GRID_NAME) == Trim(MESH_NAME(irn)) ) Then
                  new_ffield_i = irn
                  Exit Irn_Loop2
               End If
            End Do Irn_Loop2
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
            If (STR_INDX>0) Then
               Z1 = EVAC_STRS(STR_INDX)%XB_NODE(STR_SUB_INDX,5)
               Z2 = EVAC_STRS(STR_INDX)%XB_NODE(STR_SUB_INDX,6)
            End If
         End If
         If (Abs(ior) == 1) irnmax = Int(Width*4.0_EB)
         If (Abs(ior) == 2) irnmax = Int(Width*4.0_EB)
         If (Abs(ior) == 3 .Or. ior == 0) irnmax = Int((x2-x1)*4.0_EB)*Int((y2-y1)*4.0_EB)
         irnmax = Max(irnmax,5)
         istat = 1
         MFF=>Meshes(imesh2)
         !
         irn = 0
         angle = 0.0_EB
         CheckPPForce: Do While (irn < irnmax)
            Select Case (ior)
            Case(-1,1)
               ! 180 or 0 degrees, i.e., pi or 0 radians
               angle = (0.5_EB-(ior/2.0_EB))*Pi
               Call Random_number(rn)
               If (keep_xy) Then
                  yy = HR%Y + (irn/Max(irn,1))*(rn-0.5_EB)*0.25_EB*HR%Radius
                  xx = x1 + ior*(1.0_EB*HR%B + HR%Radius)
                  angle = HR%Angle
               Else
                  yy = 0.5_EB*(y1+y2) + (rn-0.5_EB)*Max(0.0_EB,Width-2.0_EB*HR%Radius-2.0_EB*HR%B)
                  xx = x1 + ior*(1.0_EB*HR%B + HR%r_torso)
               End If
            Case(-2,2)
               ! 270 or90 degrees, i.e., 3pi/2 or pi/2 radians
               angle = (1.0_EB-(ior/4.0_EB))*Pi
               Call Random_number(rn)
               If (keep_xy) Then
                  xx = HR%X + (irn/Max(irn,1))*(rn-0.5_EB)*0.25_EB*HR%Radius
                  yy = y1 + (ior/Abs(ior))*(1.0_EB*HR%B + HR%Radius)
                  angle = HR%Angle
               Else
                  xx = 0.5_EB*(x1+x2) + (rn-0.5_EB)*Max(0.0_EB,Width-2.0_EB*HR%Radius-2.0_EB*HR%B)
                  yy = y1 + (ior/Abs(ior))*(1.0_EB*HR%B + HR%r_torso)
               End If
            Case(0,3)
               Call Random_number(rn)
               yy = 0.5_EB*(y1+y2) + (rn-0.5_EB)*Max(0.0_EB,(y2-y1)-2.0_EB*HR%Radius-2.0_EB*HR%B)
               Call Random_number(rn)
               xx = 0.5_EB*(x1+x2) + (rn-0.5_EB)*Max(0.0_EB,(x2-x1)-2.0_EB*HR%Radius-2.0_EB*HR%B)
            End Select
            zz = 0.5_EB*(z1+z2)
            II = Floor( MFF%CELLSI(Floor((xx-MFF%XS)*MFF%RDXINT)) + 1.0_EB )
            JJ = Floor( MFF%CELLSJ(Floor((yy-MFF%YS)*MFF%RDYINT)) + 1.0_EB )
            KK = 1
            I_OBST = MFF%OBST_INDEX_C(MFF%CELL_INDEX(II,JJ,KK))

            irn = irn + 1

            If (MFF%SOLID(MFF%CELL_INDEX(II,JJ,KK)) .AND. .NOT.MFF%OBSTRUCTION(I_OBST)%HIDDEN) Cycle CheckPPForce
            If ( Abs(ior) == 2 .And. .Not. keep_xy ) Then
               xx1 = xx - HR%Radius - 1.0_EB*HR%B
               II = Floor(MFF%CELLSI(Floor((xx1-MFF%XS)*MFF%RDXINT))+1.0_EB)
               I_OBST = MFF%OBST_INDEX_C(MFF%CELL_INDEX(II,JJ,KK))
               If (MFF%SOLID(MFF%CELL_INDEX(II,JJ,KK)) .AND. .NOT.MFF%OBSTRUCTION(I_OBST)%HIDDEN) Cycle CheckPPForce
               xx1 = xx + HR%Radius + 1.0_EB*HR%B
               II = Floor(MFF%CELLSI(Floor((xx1-MFF%XS)*MFF%RDXINT))+1.0_EB)
               I_OBST = MFF%OBST_INDEX_C(MFF%CELL_INDEX(II,JJ,KK))
               If (MFF%SOLID(MFF%CELL_INDEX(II,JJ,KK)) .AND. .NOT.MFF%OBSTRUCTION(I_OBST)%HIDDEN) Cycle CheckPPForce
            End If
            If ( Abs(ior) == 1 .And. .Not. keep_xy ) Then
               yy1 = yy - HR%Radius - 1.0_EB*HR%B
               JJ = Floor(MFF%CELLSJ(Floor((yy1-MFF%YS)*MFF%RDYINT))+1.0_EB)
               I_OBST = MFF%OBST_INDEX_C(MFF%CELL_INDEX(II,JJ,KK))
               If (MFF%SOLID(MFF%CELL_INDEX(II,JJ,KK)) .AND. .NOT.MFF%OBSTRUCTION(I_OBST)%HIDDEN) Cycle CheckPPForce
               yy1 = yy + HR%Radius + 1.0_EB*HR%B
               JJ = Floor(MFF%CELLSJ(Floor((yy1-MFF%YS)*MFF%RDYINT))+1.0_EB)
               I_OBST = MFF%OBST_INDEX_C(MFF%CELL_INDEX(II,JJ,KK))
               If (MFF%SOLID(MFF%CELL_INDEX(II,JJ,KK)) .AND. .NOT.MFF%OBSTRUCTION(I_OBST)%HIDDEN) Cycle CheckPPForce
            End If

            !
            d_max =  1.0_EB*HR%B
            r_tmp(1) = HR%r_shoulder ! right circle
            r_tmp(2) = HR%r_torso     ! center circle
            r_tmp(3) = HR%r_shoulder ! left circle
            y_tmp(1) = yy - Cos(angle)*HR%d_shoulder ! right
            x_tmp(1) = xx + Sin(angle)*HR%d_shoulder
            y_tmp(2) = yy ! torso
            x_tmp(2) = xx
            y_tmp(3) = yy + Cos(angle)*HR%d_shoulder ! left
            x_tmp(3) = xx - Sin(angle)*HR%d_shoulder
            P2PLoop: Do ie = 1, MFF%N_HUMANS
               HRE=>MFF%HUMAN(IE)
               If (STR_SUB_INDX /= HRE%STR_SUB_INDX) Cycle P2PLoop
               r_tmp(4) = HRE%r_shoulder ! right circle
               r_tmp(5) = HRE%r_torso     ! center circle
               r_tmp(6) = HRE%r_shoulder ! left circle
               y_tmp(4) = HRE%Y - Cos(HRE%angle)*HRE%d_shoulder ! right circle
               x_tmp(4) = HRE%X + Sin(HRE%angle)*HRE%d_shoulder
               y_tmp(5) = HRE%Y ! center circle
               x_tmp(5) = HRE%X
               y_tmp(6) = HRE%Y + Cos(HRE%angle)*HRE%d_shoulder ! left circle
               x_tmp(6) = HRE%X - Sin(HRE%angle)*HRE%d_shoulder
               Do iii = 1, 3
                  Do jjj = 4, 6
                     DIST = Sqrt((x_tmp(jjj)-x_tmp(iii))**2 + (y_tmp(jjj)-y_tmp(iii))**2) - (r_tmp(jjj)+r_tmp(iii))
                     If ( DIST < d_max) Cycle CheckPPForce
                  End Do
               End Do
            End Do P2PLoop
            istat = 0    ! target is free
            ior_new = 1  ! target is a door or an entry
            Exit CheckPPForce
         End Do CheckPPForce

         If (istat == 0) Then
            j  =  Max(0,HR%GROUP_ID)
            j1 = -Min(0,HR%GROUP_ID)
            If (j > 0) Then
               Group_Known_Doors(j)%I_Target = 0
               P2PLoop2: Do ie = 1, MFF%N_HUMANS
                  HRE=>MFF%HUMAN(IE)
                  If (HRE%GROUP_ID == j .And. (HRE%Ilabel /= HR%Ilabel)) Then
                     Group_Known_Doors(j)%I_Target = HRE%I_Target
                     Exit P2PLoop2
                  End If
               End Do P2PLoop2
            End If
            If (HR%IEL > 0 ) Then
               HPT => EVACUATION(HR%IEL)
            Else
               PNX => EVAC_ENTRYS(Abs(HR%IEL))
            End If

            i_tmp = 0
            Do ie = 1, n_egrids
               If (EVAC_Node_List(ie)%IMESH == imesh2 ) Then
                  i_tmp = ie
               End If
            End Do
            If (i_tmp == 0 .Or. i_tmp > n_egrids) Then
               Write(MESSAGE,'(A,I6)') 'ERROR: Check_Target_Node, no imesh2 found ',imesh2
               Call SHUTDOWN(MESSAGE)
            End If

            If (j > 0 .And. Group_List(j)%GROUP_I_FFIELDS(i_tmp) > 0) Then
               ! There are already group members on this floor, use the same field
               new_ffield_i = Group_List(j)%GROUP_I_FFIELDS(i_tmp)
               new_ffield_name = Trim(MESH_NAME(new_ffield_i))
               If (j > 0) I_Target = Group_Known_Doors(j)%I_Target
            Else

               Call Change_Target_Door(imesh2, imesh2, 1, j, j1, 0, 2, HR%X, HR%Y, I_Target, color_index, new_ffield_i, HR)

               new_ffield_name = Trim(MESH_NAME(new_ffield_i))
               If ( j > 0 ) Then
                  Group_List(j)%GROUP_I_FFIELDS(i_tmp) = new_ffield_i
               End If
            End If  ! first member of a group?

         End If ! istat=0, i.e., put human to a new node

      Case ('Corr') SelectTargetType
         PCX2 => EVAC_CORRS(EVAC_Node_List(INODE2)%Node_Index)
         If (PCX2%n_inside < PCX2%MAX_HUMANS_INSIDE) Then
            xx = 0.0_EB ; yy = 0.0_EB ; zz = 0.0_EB ! not used for corridors
            ior_new = 3     ! target is a corridor
            imesh2 = HR%IMESH     ! do not put to an another floor yet
            HR%INODE=evac_corrs(EVAC_Node_List(INODE2)%Node_Index)%INODE
            istat = 0       ! target is free
            PCX2%n_inside = PCX2%n_inside + 1
            Allocate(TmpCurrent,STAT=IZERO)
            Call ChkMemErr('LINK','Check_Target_Node',IZERO) 
            TmpCurrent%HUMAN  = HR    ! Save the human data
            TmpCurrent%T_in   = T
            TmpCurrent%T_out  = T + PCX2%Eff_Length/(PCX2%Fac_Speed*HR%v0_fac*HR%Speed)
            TmpCurrent%From1_To2 = .True.
            TmpCurrent%Index     = PCX2%n_inside
            TmpCurrent%Next      => PCX2%First
            PCX2%First            => TmpCurrent
         Else
            ior_new = 3     ! target is a corridor
            istat = 1       ! do not enter the corridor
            imesh2 = HR%IMESH
         End If
      Case ('Floor') SelectTargetType
         ! Next is used for meshes
         !
         NM_STRS_MESH = .FALSE.
         ior_new = 4       ! target is a mesh
         istat = 1         ! Default: no space in the mesh
         Select Case (EVAC_NODE_List(INODE)%Node_Type)
         Case ('Door')
            PDX2 => EVAC_DOORS(EVAC_Node_List(INODE)%Node_Index)
            imesh2 = EVAC_Node_List(INODE2)%IMESH
            MFF=>Meshes(imesh2)
            Do N = 1, N_STRS
               IF (EVAC_STRS(N)%IMESH == imesh2) Then
                  STRS_Indx = N
                  STRP=>EVAC_STRS(N)
                  NM_STRS_MESH = .TRUE.
                  Exit
               End If
            End Do
            X1  = Min(Max(MFF%XS,PDX2%X1),MFF%XF)
            X2  = Min(Max(MFF%XS,PDX2%X2),MFF%XF)
            Y1  = Min(Max(MFF%YS,PDX2%Y1),MFF%YF)
            Y2  = Min(Max(MFF%YS,PDX2%Y2),MFF%YF)
            Z1  = HR%Z   ! use agent's z
            Z2  = HR%Z   ! use agent's z
            zz  = HR%Z   ! use agent's z
            ior = PDX2%IOR        ! direction is mesh1==>door1==>mesh2
            Width = PDX2%Width
            keep_xy2 = .True. ! do not move horizontally
            If (PDX2%STR_INDX > 0) Then
               STR_INDX = PDX2%STR_INDX            
               STR_SUB_INDX = PDX2%STR_SUB_INDX
            End If
         Case ('Entry')
            ! Todo:  Case ('Entry') and geometry part from CHECK_ENTRY
            PNX2 => EVAC_ENTRYS(EVAC_Node_List(INODE2)%Node_Index) 
            If (PNX2%STR_INDX > 0) Then
               STR_INDX = PNX2%STR_INDX            
               STR_SUB_INDX = PNX2%STR_SUB_INDX
            End If
         Case Default
            Write(MESSAGE,'(A)') 'ERROR (debug): Check_Target_Node and STRS: Node 1 is not a door.'
         End Select
         If (Abs(ior) == 1) irnmax = Int(Width*4.0_EB)
         If (Abs(ior) == 2) irnmax = Int(Width*4.0_EB)
         If (Abs(ior) == 3 .Or. ior == 0) irnmax = Int((x2-x1)*4.0_EB)*Int((y2-y1)*4.0_EB)
         irnmax = Max(irnmax,5)
         !
         irn = 0
         angle = 0.0_EB
         CheckPPForceSTRS: Do While (irn < irnmax)
            Select Case (ior)
            Case(-1,1)
               ! 180 or 0 degrees, i.e., pi or 0 radians
               angle = (0.5_EB-(ior/2.0_EB))*Pi
               Call Random_number(rn)
               If (keep_xy2) Then
                  yy = HR%Y + (irn/Max(irn,1))*(rn-0.5_EB)*0.25_EB*HR%Radius
                  xx = x1 + ior*(1.0_EB*HR%B + HR%Radius)
                  angle = HR%Angle
               Else
                  yy = 0.5_EB*(y1+y2) + (rn-0.5_EB)* Max(0.0_EB,Width-2.0_EB*HR%Radius-2.0_EB*HR%B)
                  xx = x1 + ior*(1.0_EB*HR%B + HR%r_torso)
               End If
            Case(-2,2)
               ! 270 or 90 degrees, i.e., 3pi/2 or pi/2 radians
               angle = (1.0_EB-(ior/4.0_EB))*Pi
               Call Random_number(rn)
               If (keep_xy2) Then
                  xx = HR%X + (irn/Max(irn,1))*(rn-0.5_EB)*0.25_EB*HR%Radius
                  yy = y1 + (ior/Abs(ior))*(1.0_EB*HR%B + HR%Radius)
                  angle = HR%Angle
               Else
                  xx = 0.5_EB*(x1+x2) + (rn-0.5_EB)* Max(0.0_EB,Width-2.0_EB*HR%Radius-2.0_EB*HR%B)
                  yy = y1 + (ior/Abs(ior))*(1.0_EB*HR%B + HR%r_torso)
               End If
            Case(0,3)
               Call Random_number(rn)
               yy = 0.5_EB*(y1+y2) + (rn-0.5_EB)* Max(0.0_EB,(y2-y1)-2.0_EB*HR%Radius-2.0_EB*HR%B)
               Call Random_number(rn)
               xx = 0.5_EB*(x1+x2) + (rn-0.5_EB)* Max(0.0_EB,(x2-x1)-2.0_EB*HR%Radius-2.0_EB*HR%B)
            End Select
            II = Floor( MFF%CELLSI(Floor((xx-MFF%XS)*MFF%RDXINT)) + 1.0_EB )
            JJ = Floor( MFF%CELLSJ(Floor((yy-MFF%YS)*MFF%RDYINT)) + 1.0_EB )
            KK = 1

            irn = irn + 1

            I_OBST = MFF%OBST_INDEX_C(MFF%CELL_INDEX(II,JJ,KK))
            If (MFF%SOLID(MFF%CELL_INDEX(II,JJ,KK)) .AND. .NOT.MFF%OBSTRUCTION(I_OBST)%HIDDEN) Cycle CheckPPForceSTRS

            !
            d_max =  1.0_EB*HR%B
            r_tmp(1) = HR%r_shoulder ! right circle
            r_tmp(2) = HR%r_torso     ! center circle
            r_tmp(3) = HR%r_shoulder ! left circle
            y_tmp(1) = yy - Cos(angle)*HR%d_shoulder ! right
            x_tmp(1) = xx + Sin(angle)*HR%d_shoulder
            y_tmp(2) = yy ! torso
            x_tmp(2) = xx
            y_tmp(3) = yy + Cos(angle)*HR%d_shoulder ! left
            x_tmp(3) = xx - Sin(angle)*HR%d_shoulder
            P2PLoopSTRS: Do ie = 1, MFF%N_HUMANS
               HRE=>MFF%HUMAN(IE)
               If (ABS(STR_SUB_INDX-HRE%STR_SUB_INDX)>1) Cycle P2PLoopSTRS
               r_tmp(4) = HRE%r_shoulder ! right circle
               r_tmp(5) = HRE%r_torso     ! center circle
               r_tmp(6) = HRE%r_shoulder ! left circle
               y_tmp(4) = HRE%Y - Cos(HRE%angle)*HRE%d_shoulder ! right circle
               x_tmp(4) = HRE%X + Sin(HRE%angle)*HRE%d_shoulder
               y_tmp(5) = HRE%Y ! center circle
               x_tmp(5) = HRE%X
               y_tmp(6) = HRE%Y + Cos(HRE%angle)*HRE%d_shoulder ! left circle
               x_tmp(6) = HRE%X - Sin(HRE%angle)*HRE%d_shoulder
               Do iii = 1, 3
                  Do jjj = 4, 6
                     DIST = Sqrt((x_tmp(jjj)-x_tmp(iii))**2 + (y_tmp(jjj)-y_tmp(iii))**2) - (r_tmp(jjj)+r_tmp(iii))
                     If ( DIST < d_max) Cycle CheckPPForceSTRS
                  End Do
               End Do
            End Do P2PLoopSTRS
            istat = 0    ! target is free
            Exit CheckPPForceSTRS
         End Do CheckPPForceSTRS
         If (istat == 0) Then
            If (NM_STRS_MESH) Call Find_Target_Node_In_Strs(STRP,HR)
            I_target = HR%I_Target
            new_ffield_i    = imesh2
            new_ffield_name = Trim(MESH_NAME(new_ffield_i))
         End If ! istat=0, i.e., put human to a new node
      Case ('Exit') SelectTargetType
         ior_new = 5  ! target is an exit
         istat = 0    ! remove from the floor (exit is always free)
         imesh2 = HR%IMESH
      Case Default SelectTargetType
         ior_new = 6  ! target is not defined
         istat = 1    ! do not remove from the floor
         imesh2 = HR%IMESH
      End Select SelectTargetType
      !
    End Subroutine Check_Target_Node
    !
    Subroutine REMOVE_PERSON(I)
      Implicit None
      !
      ! Remove a person
      !
      ! Passed variables
      Integer, Intent(IN) :: I
      !
      ! Local variables
      !
      HUMAN(I) = HUMAN(N_HUMANS)
      N_HUMANS = N_HUMANS - 1
      !
    End Subroutine REMOVE_PERSON
    !
    Subroutine REMOVE_OUT_OF_GRIDS(T,NM)
      Implicit None
      !
      ! Remove humans that do not lie in any mesh
      !
      ! Passed variables
      Integer, Intent(IN) :: NM
      Real(EB), Intent(IN) :: T
      !
      ! Local variables
      Integer :: IKILL, I
      !
      IKILL = 0
      DROP_LOOP: Do I=1,N_HUMANS
         !
         HR=>HUMAN(I)
         If (I > N_HUMANS-IKILL) Exit DROP_LOOP
         If (HR%X > XS .And. HR%X < XF .And. HR%Y > YS .And. HR%Y < YF) Cycle DROP_LOOP
         write(LU_ERR,'(A, 3F8.2, I3)') 'WARNING: Person removed at coord: ', HR%X,HR%Y,HR%Z,HR%SKIP_WALL_FORCE_IOR
         HUMAN(I) = HUMAN(N_HUMANS-IKILL)
         IKILL = IKILL + 1
      End Do DROP_LOOP
      N_HUMANS = N_HUMANS - IKILL
      !
    End Subroutine REMOVE_OUT_OF_GRIDS
    !
    Subroutine ENTRY_HUMAN(I_entry, Tin, NM, istat)
      Implicit None
      !
      ! Insert humans into the domain every 1/Flow seconds.
      !
      ! Passed variables
      Real(EB), Intent(IN) :: Tin
      Integer,  Intent(IN) :: I_entry, NM
      Integer, Intent(OUT) :: istat
      !
      ! Local variables
      Real(EB) RN, x1, x2, y1, y2, z1, z2, d_max, dist, xx, yy, zz, xx1, yy1
      Integer  II, JJ, KK, ior, irnmax, irn, ie, NR
      Real(EB), Dimension(6) ::y_tmp, x_tmp, r_tmp

      Type (EVAC_ENTR_Type), Pointer :: PNX =>NULL()
      Type (MESH_TYPE), Pointer :: MFF =>NULL()
      Type (EVAC_PERS_Type), Pointer :: PCP =>NULL()
      Type (HUMAN_TYPE), Pointer :: HR=>NULL(), HRE =>NULL()
      !
      istat = 1
      PNX => EVAC_ENTRYS(I_entry)
      If (PNX%IMESH /= NM ) Return
      If (PNX%Flow <= 0.0_EB ) Return
      If (PNX%T_Start > Tin) Return
      If (PNX%T_Stop < Tin) Return
      If (PNX%Max_Humans < 0) Then
         NR     = -PNX%Max_Humans  
         If (PNX%ICOUNT >= INT(EVALUATE_RAMP(Tin,0._EB,NR))) Return
      Else
         If (PNX%ICOUNT >= PNX%Max_Humans) Return
      Endif
      MFF => MESHES(NM)
      If ( (Tin-PNX%T_last) < (1.0_EB/PNX%Flow) ) Return
      X1  = PNX%X1
      X2  = PNX%X2
      Y1  = PNX%Y1
      Y2  = PNX%Y2
      Z1  = PNX%Z1
      Z2  = PNX%Z2
      ior = PNX%IOR
      !
      If (N_HUMANS+1 > N_HUMANS_DIM) Then
         ! Re-allocation is not yet checked.
         Call SHUTDOWN('ERROR: Insert Humans: no re-allocation yet')
         Call RE_ALLOCATE_HUMANS(1,NM)
         HUMAN=>MESHES(NM)%HUMAN
      End If
      !
      PCP => EVAC_PERSON_CLASSES(PNX%IPC)
      HR  => MESHES(NM)%HUMAN(N_HUMANS+1)
      Call CLASS_PROPERTIES(HR,PCP)
      HR%Tpre = 0.0_EB
      HR%Tdet = T_BEGIN
      HR%IPC  = PNX%IPC
      HR%IEL  = -I_entry
      HR%GROUP_ID = 0
      HR%I_Target = 0
      !
      If (Abs(ior) == 1) irnmax = Int(PNX%Width*4.0_EB)
      If (Abs(ior) == 2) irnmax = Int(PNX%Width*4.0_EB)
      If (Abs(ior) == 3 .Or. ior == 0) irnmax = Int((x2-x1)*4.0_EB)*Int((y2-y1)*4.0_EB)
      irnmax = Max(irnmax,5)
      !
      irn = 0
      CheckPPForce: Do While (irn < irnmax)
         Select Case (ior)
         Case(-1,1)
            Call Random_number(rn)
            yy = 0.5_EB*(y1+y2) + (rn-0.5_EB)*Max(0.0_EB,PNX%Width-2.0_EB*HR%Radius-2.0_EB*HR%B)
            xx = x1 + ior*5.0_EB*HR%B
            HR%Angle = (1-ior)*Pi/2.0_EB  ! ior=1: 0,  ior=-1: pi
         Case(-2,2)
            Call Random_number(rn)
            xx = 0.5_EB*(x1+x2) + (rn-0.5_EB)*Max(0.0_EB,PNX%Width-2.0_EB*HR%Radius-2.0_EB*HR%B)
            yy = y1 + (ior/Abs(ior))*5.0_EB*HR%B
            HR%Angle = Pi/2.0_EB + (2-ior)*Pi/4.0_EB  ! ior=2: (3/2)pi,  ior=-2: pi/2
         Case(0,3)
            Call Random_number(rn)
            yy = 0.5_EB*(y1+y2) + (rn-0.5_EB)*Max(0.0_EB,(y2-y1)-2.0_EB*HR%Radius-2.0_EB*HR%B)
            Call Random_number(rn)
            xx = 0.5_EB*(x1+x2) + (rn-0.5_EB)*Max(0.0_EB,(x2-x1)-2.0_EB*HR%Radius-2.0_EB*HR%B)
            Call Random_number(rn)
            HR%Angle = 2.0_EB*Pi*rn
         End Select
         Do While (HR%Angle >= 2.0_EB*Pi)
            HR%Angle = HR%Angle - 2.0_EB*Pi
         End Do
         Do While (HR%Angle < 0.0_EB)
            HR%Angle = HR%Angle + 2.0_EB*Pi
         End Do
         zz  = 0.5_EB*(z1+z2)
         II = Floor( CELLSI(Floor((xx-XS)*RDXINT)) + 1.0_EB )
         JJ = Floor( CELLSJ(Floor((yy-YS)*RDYINT)) + 1.0_EB )
         KK = 1

         irn = irn + 1

         I_OBST = OBST_INDEX_C(CELL_INDEX(II,JJ,KK))
         If (SOLID(CELL_INDEX(II,JJ,KK)) .AND. .NOT.OBSTRUCTION(I_OBST)%HIDDEN) Cycle CheckPPForce
         If ( Abs(ior) == 2 ) Then
            xx1 = xx - HR%Radius - 1.0_EB*HR%B
            II = Floor(CELLSI(Floor((xx1-XS)*RDXINT))+1.0_EB)
            I_OBST = OBST_INDEX_C(CELL_INDEX(II,JJ,KK))
            If (SOLID(CELL_INDEX(II,JJ,KK)) .AND. .NOT.OBSTRUCTION(I_OBST)%HIDDEN) Cycle CheckPPForce
            xx1 = xx + HR%Radius + 1.0_EB*HR%B
            II = Floor(CELLSI(Floor((xx1-XS)*RDXINT))+1.0_EB)
            I_OBST = OBST_INDEX_C(CELL_INDEX(II,JJ,KK))
            If (SOLID(CELL_INDEX(II,JJ,KK)) .AND. .NOT.OBSTRUCTION(I_OBST)%HIDDEN) Cycle CheckPPForce
         End If
         If ( Abs(ior) == 1 ) Then
            yy1 = yy - HR%Radius - 1.0_EB*HR%B
            JJ = Floor(CELLSJ(Floor((yy1-YS)*RDYINT))+1.0_EB)
            I_OBST = OBST_INDEX_C(CELL_INDEX(II,JJ,KK))
            If (SOLID(CELL_INDEX(II,JJ,KK)) .AND. .NOT.OBSTRUCTION(I_OBST)%HIDDEN) Cycle CheckPPForce
            yy1 = yy + HR%Radius + 1.0_EB*HR%B
            JJ = Floor(CELLSJ(Floor((yy1-YS)*RDYINT))+1.0_EB)
            I_OBST = OBST_INDEX_C(CELL_INDEX(II,JJ,KK))
            If (SOLID(CELL_INDEX(II,JJ,KK)) .AND. .NOT.OBSTRUCTION(I_OBST)%HIDDEN) Cycle CheckPPForce
         End If
         !
         d_max = 1.0_EB*HR%B
         r_tmp(1) = HR%r_shoulder ! right circle
         r_tmp(2) = HR%r_torso     ! center circle
         r_tmp(3) = HR%r_shoulder ! left circle
         y_tmp(1) = yy - Cos(HR%angle)*HR%d_shoulder ! right
         x_tmp(1) = xx + Sin(HR%angle)*HR%d_shoulder
         y_tmp(2) = yy ! torso
         x_tmp(2) = xx
         y_tmp(3) = yy + Cos(HR%angle)*HR%d_shoulder ! left
         x_tmp(3) = xx - Sin(HR%angle)*HR%d_shoulder
         P2PLoop: Do ie = 1, n_humans
            HRE=>HUMAN(IE)
            r_tmp(4) = HRE%r_shoulder ! right circle
            r_tmp(5) = HRE%r_torso     ! center circle
            r_tmp(6) = HRE%r_shoulder ! left circle
            y_tmp(4) = HRE%Y - Cos(HRE%angle)*HRE%d_shoulder ! right circle
            x_tmp(4) = HRE%X + Sin(HRE%angle)*HRE%d_shoulder
            y_tmp(5) = HRE%Y ! center circle
            x_tmp(5) = HRE%X
            y_tmp(6) = HRE%Y + Cos(HRE%angle)*HRE%d_shoulder ! left circle
            x_tmp(6) = HRE%X - Sin(HRE%angle)*HRE%d_shoulder
            Do iii = 1, 3
               Do jjj = 4, 6
                  DIST = Sqrt((x_tmp(jjj)-x_tmp(iii))**2 + (y_tmp(jjj)-y_tmp(iii))**2) - (r_tmp(jjj)+r_tmp(iii))
                  If (DIST < d_max) Cycle CheckPPForce
               End Do
            End Do
         End Do P2PLoop
         istat = 0
         Exit CheckPPForce
      End Do CheckPPForce

      If (istat == 0 ) Then
         N_HUMANS = N_HUMANS + 1
         PNX%T_last = Tin
         PNX%ICOUNT = PNX%ICOUNT + 1
         If (PNX%T_first <= T_BEGIN) PNX%T_first = Tin
         HR%X = xx
         HR%Y = yy
         HR%Z = zz
         ! 
         ILABEL_last = ILABEL_last + 1
         HR%ILABEL = ILABEL_last
         HR%SHOW = .True.    
         HR%COLOR_INDEX = 1
         Select Case (COLOR_METHOD)
         Case (-1)
            HR%COLOR_INDEX = 1
         Case (0)
            HR%COLOR_INDEX = PNX%Avatar_Color_Index
         Case (1,2)
            HR%COLOR_INDEX = 1    ! lonely human
         Case (3)
            HR%COLOR_INDEX = evac_person_classes(PNX%IPC)%Avatar_Color_Index
         Case (4)
            HR%COLOR_INDEX = 1
         Case (5)
            ! Correct color is put, where the flow fields are chosen.
            HR%COLOR_INDEX = 1
         Case Default
            Write(MESSAGE,'(A,I3,A)') 'ERROR: ENTRY_HUMAN COLOR METHOD',COLOR_METHOD, ' is not defined'
            Call SHUTDOWN(MESSAGE)
         End Select
         HR%FFIELD_NAME = Trim(PNX%GRID_NAME)
         HR%I_FFIELD    = 0
         Mesh2Loop: Do i = 1, nmeshes
            If ( evacuation_only(i) .And. Trim(HR%FFIELD_NAME) == Trim(MESH_NAME(i)) ) Then
               HR%I_FFIELD = i
               Exit Mesh2Loop
            End If
         End Do Mesh2Loop
         If ( HR%I_FFIELD == 0 ) Then
            Write(MESSAGE,'(A,A,A,A)') 'ERROR: ENTR line ',Trim(PNX%ID), ' problem with flow field name, ', &
                 Trim(PNX%GRID_NAME),' not found'
            Call SHUTDOWN(MESSAGE)
         End If
         HR%IMESH       = PNX%IMESH
         HR%INODE       = PNX%TO_INODE
         HR%NODE_NAME   = Trim(PNX%TO_NODE)
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
         HR%NewRnd    = .True.
         HR%STR_SUB_INDX = PNX%STR_SUB_INDX
         HR%SKIP_WALL_FORCE_IOR = 0
         Write (LU_EVACOUT,fmt='(a,i6,a,f8.2,a,3a)') ' EVAC: Person n:o', HR%ILABEL, ' inserted ', Tin, &
              ' s, entry ', Trim(PNX%ID),' ffield ', Trim(HR%FFIELD_NAME)
         Write (LU_EVACOUT,fmt='(a,a)') ' person     x       y       z    Tpre    Tdet  ', &
         ' dia    v0   tau   i_gr i_ff'
         Write (LU_EVACOUT,fmt='(i6,5f8.2,3f6.2,i6,i4,i4)') HR%ILABEL, &
              HR%X, HR%Y, HR%Z, HR%Tpre, HR%Tdet,2.0_EB*HR%Radius, &
              HR%Speed, HR%Tau, HR%GROUP_ID, HR%i_ffield, HR%COLOR_INDEX
      End If
      ! 
    End Subroutine ENTRY_HUMAN

    Subroutine Corner_Forces(x1, y1, x11, y11, p2p_dist_max, P2P_U, P2P_V, Social_F, &
         Contact_F, P2P_Torque, d_walls, x_tmp, y_tmp, r_tmp, u_tmp, v_tmp, istat)
      Implicit None

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
      Real(EB), Intent(IN) :: x1, y1, x11, y11, p2p_dist_max
      Real(EB), Dimension(6), Intent(IN) :: x_tmp, y_tmp, r_tmp, u_tmp, v_tmp
      Real(EB), Intent(INOUT) :: P2P_U, P2P_V, Social_F, Contact_F, P2P_Torque, d_walls
      Integer, Intent(OUT) :: istat
      !
      ! Local variables
      Real(EB) :: dist, CosPhiFac, u, v, Fc_x, Fc_y, FricFac, k_fric, evel
      Integer :: iii

      istat = 0
      If (I_Fric_sw >= 1 ) Then
         FricFac = 0.0_EB
         k_fric = HR%Kappa
      Else
         FricFac = 1.0_EB
         k_fric = HR%Gamma
      End If

      dist = Sqrt((y11-y1)**2 + (x11-x1)**2)
      If (dist-HR%Radius > p2p_dist_max) Then
         istat = 1
         Return
      End If
      u = u_tmp(2) ; v = v_tmp(2)

      If ( (u**2+v**2) > 0.0_EB ) Then
         CosPhiFac = ((x11-X1)*u + (y11-Y1)*v) / (dist*Sqrt(u**2+v**2))
         CosPhiFac = LambdaW + 0.5_EB*(1.0_EB-LambdaW)*(1.0_EB+CosPhiFac)
      Else
         CosPhiFac = 1.0_EB
      End If
      P2P_U = P2P_U + (X1-x11)*A_Wall*CosPhiFac*Exp(-(dist-HR%Radius)/B_Wall) / dist
      P2P_V = P2P_V + (Y1-y11)*A_Wall*CosPhiFac*Exp(-(dist-HR%Radius)/B_Wall) / dist
      Social_F = Social_F + Abs(A_Wall*CosPhiFac*Exp(-(dist-HR%Radius)/B_Wall))

      Do iii = 1, 3
         dist = Sqrt( (y11-y_tmp(iii))**2 + (x11-x_tmp(iii))**2 )

         ! Next is |vector1|*|vector2|
         evel = Sqrt((x11-x_tmp(iii))**2+(y11-y_tmp(iii))**2)*Sqrt(u_tmp(iii)**2+v_tmp(iii)**2)
         If (evel > 0.0_EB) evel = ((x11-x_tmp(iii))*u_tmp(iii) + &
              (y11-y_tmp(iii))*v_tmp(iii)) / evel   ! cos theta (scal_prod/(lenght1*length2)
         If (evel > 0.01_EB) Then
            d_walls = Min( (dist-r_tmp(iii))/evel, d_walls)
         Else
            d_walls = Min( (dist-r_tmp(iii))/0.01_EB, d_walls)
         End If

         ! d_walls = Min(dist-r_tmp(iii),d_walls)

         If (dist <= r_tmp(iii)) Then
            Fc_x = (x_tmp(iii)-x11)* 2.0_EB*HR%C_Young*(r_tmp(iii)-dist)/dist
            Fc_y = (y_tmp(iii)-y11)* 2.0_EB*HR%C_Young*(r_tmp(iii)-dist)/dist
            !Only radial contact forces, i.e., pressure calculation
            Contact_F = Contact_F + Sqrt(Fc_x**2 + Fc_y**2)
            ! Tangential contact force:
               Fc_x = Fc_x - k_fric*( (1.0_EB-FricFac)*(r_tmp(iii)-dist)+FricFac ) *(y_tmp(iii)-y11)* &
                    ( (y_tmp(iii)-y11)*u_tmp(iii) - (x_tmp(iii)-x11)*v_tmp(iii) )/dist**2
               Fc_y = Fc_y + k_fric*( (1.0_EB-FricFac)*(r_tmp(iii)-dist)+FricFac ) *(x_tmp(iii)-x11)* &
                    ( (y_tmp(iii)-y11)*u_tmp(iii) - (x_tmp(iii)-x11)*v_tmp(iii) )/dist**2
            P2P_Torque = P2P_Torque + Fc_y*(x11-x_tmp(iii)) - Fc_x*(y11-y_tmp(iii))
            P2P_U = P2P_U + Fc_x
            P2P_V = P2P_V + Fc_y
         End If
      End Do

    End Subroutine Corner_Forces

    Subroutine Door_Forces(nm, x_tmp, y_tmp, r_tmp, u_tmp, v_tmp, p2p_dist_max, d_xy,&
         P2P_U, P2P_V, Social_F, Contact_F, P2P_Torque, FoundWall_xy)
      Implicit None
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
      Integer, Intent(IN) :: nm
      Real(EB), Intent(IN) :: p2p_dist_max
      Real(EB), Dimension(4), Intent(IN) :: d_xy
      Logical, Dimension(4), Intent(IN) :: FoundWall_xy
      Real(EB), Dimension(6), Intent(IN) :: x_tmp, y_tmp, r_tmp, u_tmp, v_tmp
      Real(EB), Intent(INOUT) :: P2P_U, P2P_V, Social_F, Contact_F, P2P_Torque
      !
      ! Local variables
      Integer :: is, idir, iin, jjn, istat
      Real(EB) :: CosPhiFac, dist, dist1, dist2

      
      ! Check if there are doors (vents with vel >0)
      Do ii = 1, N_VENT
         If (Abs(VENTS(ii)%IOR)>2 .Or. SURFACE(VENTS(ii)%IBC)%VEL<=0) Cycle
         dist1 = Sqrt((VENTS(ii)%x1-x_tmp(2))**2 + (VENTS(ii)%y1-y_tmp(2))**2) ! door - agent centre distance
         dist2 = Sqrt((VENTS(ii)%x2-x_tmp(2))**2 + (VENTS(ii)%y2-y_tmp(2))**2)
         If ( (Min(dist1,dist2)-HR%Radius) > P2P_DIST_MAX) Cycle

         
         If (VENTS(ii)%IOR== -1 .And. FoundWall_xy(2)) Cycle
         If (VENTS(ii)%IOR== +1 .And. FoundWall_xy(1)) Cycle
         If (VENTS(ii)%IOR== -2 .And. FoundWall_xy(4)) Cycle
         If (VENTS(ii)%IOR== +2 .And. FoundWall_xy(3)) Cycle
         Select Case(VENTS(ii)%IOR)
         Case(-1)  ! wall at +x direction
            If (.Not.FoundWall_xy(2) .And. Abs(VENTS(ii)%x1-d_xy(2))<0.01_EB .And. &
                 (VENTS(ii)%y1<y_tmp(2).And.VENTS(ii)%y2>y_tmp(2))) Then
               ! There is a outflow vent here
               x11 = VENTS(ii)%x1
               y11 = VENTS(ii)%y1
               Call Corner_Forces(x1, y1, x11, y11, p2p_dist_max, P2P_U, P2P_V, Social_F, &
                    Contact_F, P2P_Torque, d_walls, x_tmp, y_tmp, r_tmp, u_tmp, v_tmp, istat)
               x11 = VENTS(ii)%x2
               y11 = VENTS(ii)%y2
               Call Corner_Forces(x1, y1, x11, y11, p2p_dist_max, P2P_U, P2P_V, Social_F, &
                    Contact_F, P2P_Torque, d_walls, x_tmp, y_tmp, r_tmp, u_tmp, v_tmp, istat)
            End If
         Case(+1)  ! wall at -x direction
            If (.Not.FoundWall_xy(1) .And. Abs(VENTS(ii)%x1-d_xy(1))<0.01_EB .And. &
                 (VENTS(ii)%y1<y_tmp(2).And.VENTS(ii)%y2>y_tmp(2))) Then
               ! There is a outflow vent here
               x11 = VENTS(ii)%x1
               y11 = VENTS(ii)%y1
               Call Corner_Forces(x1, y1, x11, y11, p2p_dist_max, P2P_U, P2P_V, Social_F, &
                    Contact_F, P2P_Torque, d_walls, x_tmp, y_tmp, r_tmp, u_tmp, v_tmp, istat)
               x11 = VENTS(ii)%x2
               y11 = VENTS(ii)%y2
               Call Corner_Forces(x1, y1, x11, y11, p2p_dist_max, P2P_U, P2P_V, Social_F, &
                    Contact_F, P2P_Torque, d_walls, x_tmp, y_tmp, r_tmp, u_tmp, v_tmp, istat)
            End If
         Case(-2)  ! wall at +y direction
            If (.Not.FoundWall_xy(4) .And. Abs(VENTS(ii)%y1-d_xy(4))<0.01_EB .And. &
                 (VENTS(ii)%x1<x_tmp(2).And.VENTS(ii)%x2>x_tmp(2))) Then
               ! There is a outflow vent here
               x11 = VENTS(ii)%x1
               y11 = VENTS(ii)%y1
               Call Corner_Forces(x1, y1, x11, y11, p2p_dist_max, P2P_U, P2P_V, Social_F, &
                    Contact_F, P2P_Torque, d_walls, x_tmp, y_tmp, r_tmp, u_tmp, v_tmp, istat)
               x11 = VENTS(ii)%x2
               y11 = VENTS(ii)%y2
               Call Corner_Forces(x1, y1, x11, y11, p2p_dist_max, P2P_U, P2P_V, Social_F, &
                    Contact_F, P2P_Torque, d_walls, x_tmp, y_tmp, r_tmp, u_tmp, v_tmp, istat)
            End If
         Case(+2)  ! wall at -y direction
            If (.Not.FoundWall_xy(3) .And. Abs(VENTS(ii)%y1-d_xy(3))<0.01_EB .And. &
                 (VENTS(ii)%x1<x_tmp(2).And.VENTS(ii)%x2>x_tmp(2))) Then
               ! There is a outflow vent here
               x11 = VENTS(ii)%x1
               y11 = VENTS(ii)%y1
               Call Corner_Forces(x1, y1, x11, y11, p2p_dist_max, P2P_U, P2P_V, Social_F, &
                    Contact_F, P2P_Torque, d_walls, x_tmp, y_tmp, r_tmp, u_tmp, v_tmp, istat)
               x11 = VENTS(ii)%x2
               y11 = VENTS(ii)%y2
               Call Corner_Forces(x1, y1, x11, y11, p2p_dist_max, P2P_U, P2P_V, Social_F, &
                    Contact_F, P2P_Torque, d_walls, x_tmp, y_tmp, r_tmp, u_tmp, v_tmp, istat)
            End If
         End Select
      End Do

    End Subroutine Door_Forces

    Subroutine Wall_SocialForces(nm, x_tmp, y_tmp, r_tmp, p2p_dist_max, d_xy, P2P_U, P2P_V, Social_F, FoundWall_xy)
      Implicit None
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
      Integer, Intent(IN) :: nm
      Real(EB), Intent(IN) :: p2p_dist_max
      Real(EB), Dimension(6), Intent(IN) :: x_tmp, y_tmp, r_tmp
      Real(EB), Dimension(4), Intent(IN) :: d_xy
      Logical, Dimension(4), Intent(IN) :: FoundWall_xy
      Real(EB), Intent(INOUT) :: P2P_U, P2P_V, Social_F
      !
      ! Local variables
      Integer :: is, idir, iii
      Real(EB) :: CosPhiFac, dist, F_soc, F_tmp

      ! -x direction
      is   = -1
      idir =  1
      F_soc = 0.0_EB
      Do iii = 1,3
         dist = Abs(d_xy(idir) - x_tmp(iii)) ! wall - agent centre distance
         If (dist-r_tmp(iii) <= P2P_DIST_MAX .And. FoundWall_xy(idir)) Then
            If ( (HR%U**2+HR%V**2) > 0.0_EB ) Then
               CosPhiFac = (is*HR%U)/Sqrt(HR%U**2+HR%V**2)
               CosPhiFac = LambdaW + 0.5_EB*(1.0_EB-LambdaW)*(1.0_EB+CosPhiFac)
            Else
               CosPhiFac = 1.0_EB
            End If
            F_tmp = - is*A_Wall*CosPhiFac*Exp( -(dist-r_tmp(iii))/B_Wall )
            If (Abs(F_tmp) > Abs(F_soc)) F_soc = F_tmp
         End If
      End Do
      P2P_U = P2P_U + F_soc
      Social_F = Social_F + Abs(F_soc)

      is   = +1
      idir =  2
      F_soc = 0.0_EB
      Do iii = 1,3
         dist = Abs(d_xy(idir) - x_tmp(iii)) ! wall - agent centre distance
         If (dist-r_tmp(iii) <= P2P_DIST_MAX .And. FoundWall_xy(idir)) Then
            If ( (HR%U**2+HR%V**2) > 0.0_EB ) Then
               CosPhiFac = (is*HR%U)/Sqrt(HR%U**2+HR%V**2)
               CosPhiFac = LambdaW + 0.5_EB*(1.0_EB-LambdaW)*(1.0_EB+CosPhiFac)
            Else
               CosPhiFac = 1.0_EB
            End If
            F_tmp = - is*A_Wall*CosPhiFac*Exp( -(dist-r_tmp(iii))/B_Wall )
            If (Abs(F_tmp) > Abs(F_soc)) F_soc = F_tmp
         End If
      End Do
      P2P_U = P2P_U + F_soc
      Social_F = Social_F + Abs(F_soc)

      is   = -1
      idir =  3
      F_soc = 0.0_EB
      Do iii = 1,3
         dist = Abs(d_xy(idir) - y_tmp(iii)) ! wall - agent centre distance
         If (dist-r_tmp(iii) <= P2P_DIST_MAX .And. FoundWall_xy(idir)) Then
            If ( (HR%U**2+HR%V**2) > 0.0_EB ) Then
               CosPhiFac = (is*HR%V)/Sqrt(HR%U**2+HR%V**2)
               CosPhiFac = LambdaW + 0.5_EB*(1.0_EB-LambdaW)*(1.0_EB+CosPhiFac)
            Else
               CosPhiFac = 1.0_EB
            End If
            F_tmp = - is*A_Wall*CosPhiFac*Exp( -(dist-r_tmp(iii))/B_Wall )
            If (Abs(F_tmp) > Abs(F_soc)) F_soc = F_tmp
         End If
      End Do
      P2P_V = P2P_V + F_soc
      Social_F = Social_F + Abs(F_soc)

      is   = +1
      idir =  4
      F_soc = 0.0_EB
      Do iii = 1,3
         dist = Abs(d_xy(idir) - y_tmp(iii)) ! wall - agent centre distance
         If (dist-r_tmp(iii) <= P2P_DIST_MAX .And. FoundWall_xy(idir)) Then
            If ( (HR%U**2+HR%V**2) > 0.0_EB ) Then
               CosPhiFac = (is*HR%V)/Sqrt(HR%U**2+HR%V**2)
               CosPhiFac = LambdaW + 0.5_EB*(1.0_EB-LambdaW)*(1.0_EB+CosPhiFac)
            Else
               CosPhiFac = 1.0_EB
            End If
            F_tmp = - is*A_Wall*CosPhiFac*Exp( -(dist-r_tmp(iii))/B_Wall )
            If (Abs(F_tmp) > Abs(F_soc)) F_soc = F_tmp
         End If
      End Do
      P2P_V = P2P_V + F_soc
      Social_F = Social_F + Abs(F_soc)

    End Subroutine Wall_SocialForces
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

    Subroutine Wall_ContactForces(nm, x_tmp, y_tmp, r_tmp, u_tmp, v_tmp, d_xy,&
         P2P_U, P2P_V, P2P_Torque, Contact_F, d_walls, FoundWall_xy)
      Implicit None
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
      Integer, Intent(IN) :: nm
      Real(EB), Intent(IN) :: x_tmp, y_tmp, r_tmp, u_tmp, v_tmp
      Real(EB), Dimension(4), Intent(IN) :: d_xy
      Logical, Dimension(4), Intent(IN) :: FoundWall_xy
      Real(EB), Intent(INOUT) :: P2P_U, P2P_V, P2P_Torque, Contact_F, d_walls
      !
      ! Local variables
      Integer :: is, idir
      Real(EB) :: Fc_y, Fc_x, dist, evel

      ! -x direction
      is   = -1
      idir =  1
      dist = Abs(d_xy(idir) - x_tmp)

      If (FoundWall_xy(idir)) Then
         ! Next is |vector1|*|vector2|
         evel = Abs(d_xy(idir)-x_tmp)*Sqrt(u_tmp**2+v_tmp**2)
         ! Next is cos(theta) = (scal_prod/(lenght1*length2)
         If (evel > 0.0_EB) evel = ((d_xy(idir)-x_tmp)*u_tmp) / evel
         If (evel > 0.01_EB) Then
            d_walls = Min( (dist-r_tmp)/evel, d_walls)
         Else
            d_walls = Min( (dist-r_tmp)/0.01_EB, d_walls)
         End If
         ! d_walls = Min(dist-r_tmp,d_walls) ! wall - circle centre distance
      End If

      If (dist-r_tmp <= 0.0_EB .And. FoundWall_xy(idir)) Then
         Fc_x = -is*2.0_EB*HR%C_Young*(r_tmp-dist)
         Fc_y = 0.0_EB
         Contact_F = Contact_F + Abs(Fc_x)
         If (I_Fric_sw >= 1 ) Then
            Fc_y = Fc_y - HR%Kappa*(r_tmp-dist)*v_tmp
         Else
            Fc_y = Fc_y - HR%Gamma*v_tmp
         End If
         !write(lu_evacout,*)'*** -x Fc',Fc_x,Fc_y,x_tmp,y_tmp,u_tmp,v_tmp,icyc
         P2P_Torque = P2P_Torque + Fc_y*(d_xy(idir)-HR%X) - Fc_x*(y_tmp-HR%Y)
         P2P_U = P2P_U + Fc_x
         P2P_V = P2P_V + Fc_y
      End If

      ! +x direction
      is   = +1
      idir =  2
      dist = Abs(d_xy(idir) - x_tmp)

      If (FoundWall_xy(idir)) Then
         ! Next is |vector1|*|vector2|
         evel = Abs(d_xy(idir)-x_tmp)*Sqrt(u_tmp**2+v_tmp**2)
         ! Next is cos(theta) = (scal_prod/(lenght1*length2)
         If (evel > 0.0_EB) evel = ((d_xy(idir)-x_tmp)*u_tmp) / evel
         If (evel > 0.01_EB) Then
            d_walls = Min( (dist-r_tmp)/evel, d_walls)
         Else
            d_walls = Min( (dist-r_tmp)/0.01_EB, d_walls)
         End If
         ! d_walls = Min(dist-r_tmp,d_walls) ! wall - circle centre distance
      End If

      If (dist-r_tmp <= 0.0_EB .And. FoundWall_xy(idir)) Then
         Fc_x = -is*2.0_EB*HR%C_Young*(r_tmp-dist)
         Fc_y = 0.0_EB
         Contact_F = Contact_F + Abs(Fc_x)
         If (I_Fric_sw >= 1 ) Then
            Fc_y = Fc_y - HR%Kappa*(r_tmp-dist)*v_tmp
         Else
            Fc_y = Fc_y - HR%Gamma*v_tmp
         End If
         !write(lu_evacout,*)'*** +x Fc',Fc_x,Fc_y,x_tmp,y_tmp,u_tmp,v_tmp,icyc
         P2P_Torque = P2P_Torque + Fc_y*(d_xy(idir)-HR%X) - Fc_x*(y_tmp-HR%Y)
         P2P_U = P2P_U + Fc_x
         P2P_V = P2P_V + Fc_y
      End If

      ! -y direction
      is   = -1
      idir =  3
      dist = Abs(d_xy(idir) - y_tmp)

      If (FoundWall_xy(idir)) Then
         ! Next is |vector1|*|vector2|
         evel = Abs(d_xy(idir)-y_tmp)*Sqrt(u_tmp**2+v_tmp**2)
         ! Next is cos(theta) = (scal_prod/(lenght1*length2)
         If (evel > 0.0_EB) evel = ((d_xy(idir)-y_tmp)*v_tmp) / evel
         If (evel > 0.01_EB) Then
            d_walls = Min( (dist-r_tmp)/evel, d_walls)
         Else
            d_walls = Min( (dist-r_tmp)/0.01_EB, d_walls)
         End If
         !d_walls = Min(dist-r_tmp,d_walls) ! wall - circle centre distance
      End If

      If (dist-r_tmp <= 0.0_EB .And. FoundWall_xy(idir)) Then
         Fc_y = -is*2.0_EB*HR%C_Young*(r_tmp-dist)
         Fc_x = 0.0_EB
         Contact_F = Contact_F + Abs(Fc_y)
         If (I_Fric_sw >= 1 ) Then
            Fc_x = Fc_x - HR%Kappa*(r_tmp-dist)*u_tmp
         Else
            Fc_x = Fc_x - HR%Gamma*u_tmp
         End If
         !write(lu_evacout,*)'*** -y Fc',Fc_x,Fc_y,x_tmp,y_tmp,u_tmp,v_tmp,icyc
         P2P_Torque = P2P_Torque + Fc_y*(x_tmp-HR%X) - Fc_x*(d_xy(idir)-HR%Y)
         P2P_U = P2P_U + Fc_x
         P2P_V = P2P_V + Fc_y
      End If

      ! +y direction
      is   = +1
      idir =  4
      dist = Abs(d_xy(idir) - y_tmp)

      If (FoundWall_xy(idir)) Then
         ! Next is |vector1|*|vector2|
         evel = Abs(d_xy(idir)-y_tmp)*Sqrt(u_tmp**2+v_tmp**2)
         ! Next is cos(theta) = (scal_prod/(lenght1*length2)
         If (evel > 0.0_EB) evel = ((d_xy(idir)-y_tmp)*v_tmp) / evel
         If (evel > 0.01_EB) Then
            d_walls = Min( (dist-r_tmp)/evel, d_walls)
         Else
            d_walls = Min( (dist-r_tmp)/0.01_EB, d_walls)
         End If
         !d_walls = Min(dist-r_tmp,d_walls) ! wall - circle centre distance
      End If

      If (dist-r_tmp <= 0.0_EB .And. FoundWall_xy(idir)) Then
         Fc_y = -is*2.0_EB*HR%C_Young*(r_tmp-dist)
         Fc_x = 0.0_EB
         Contact_F = Contact_F + Abs(Fc_y)
         If (I_Fric_sw >= 1 ) Then
            Fc_x = Fc_x - HR%Kappa*(r_tmp-dist)*u_tmp
         Else
            Fc_x = Fc_x - HR%Gamma*u_tmp
         End If
         !write(lu_evacout,*)'*** +y Fc',Fc_x,Fc_y,x_tmp,y_tmp,u_tmp,v_tmp,icyc
         P2P_Torque = P2P_Torque + Fc_y*(x_tmp-HR%X) - Fc_x*(d_xy(idir)-HR%Y)
         P2P_U = P2P_U + Fc_x
         P2P_V = P2P_V + Fc_y
      End If

    End Subroutine Wall_ContactForces
    !
  End Subroutine EVACUATE_HUMANS
! ============================================================
! EVACUATE_HUMANS ENDS HERE.
! ============================================================
!
! ============================================================
! NEXT ARE EVAC MODULE SUBPROGRMAS
! ============================================================
!
  Subroutine CLASS_PROPERTIES(HR,PCP)
    Implicit None
    !
    ! Passed variables
    Type (HUMAN_TYPE), Pointer:: HR
    Type (EVAC_PERS_Type), Pointer:: PCP
    !
    ! Local variables
    ! How many rnd numbers per one call to the rnd routine
    Integer, Parameter  :: n_rnd=1, n_max_par=4
    Integer  :: n_par, RandomType
    ! No more than 4 numbers needed to specify the distributions
    Real(EB) :: RandomPara(n_max_par)
    Real(EB) :: rnd_vec(n_rnd)
    ! 1: uniform (TESTED: OK)
    ! 2: Truncated normal (TESTED: OK)
    ! 3: gamma  (TESTED: OK)
    ! 4: normal (TESTED: OK)
    ! 5: lognormal (TESTED: OK)
    ! 6: beta (TESTED: OK)
    ! 7: Triangular (TESTED: OK)
    ! 8: Weibull (TESTED: OK) (alpha=1: Exponential)
    ! 9: Gumbel (TESTED: OK)

    Select Case(PCP%I_VEL_DIST)
    Case(-1)
       Call SHUTDOWN('ERROR: Class_Properties: -1')
    Case(0)
       HR%Speed  = PCP%V_mean
    Case(1)   ! Uniform
       ! Parameters: (ave,min,max) ave not used
       n_par = 3
       Randomtype = 1
       RandomPara(1) = 0.5_EB*(PCP%V_high+PCP%V_low)
       RandomPara(2) = PCP%V_low
       RandomPara(3) = PCP%V_high
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Speed = rnd_vec(1)
    Case(2)   ! Truncated Normal
       ! Parameters: (ave,sigma,min,max)
       n_par = 4
       Randomtype = 2
       RandomPara(1) = PCP%V_mean
       RandomPara(2) = PCP%V_para
       RandomPara(3) = PCP%V_low
       RandomPara(4) = PCP%V_high
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Speed = rnd_vec(1)
    Case(3)   ! Gamma
       ! Parameters: (ave,alpha,beta) ave not used
       n_par = 3
       Randomtype = 3
       RandomPara(1) = PCP%V_mean
       RandomPara(2) = PCP%V_para
       RandomPara(3) = PCP%V_para2
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Speed = rnd_vec(1)
    Case(4)   ! Normal
       ! Parameters: (ave,sigma)
       n_par = 2
       Randomtype = 4
       RandomPara(1) = PCP%V_mean
       RandomPara(2) = PCP%V_para
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Speed = rnd_vec(1)
    Case(5)   ! LogNormal
       ! mean and variance of log(x) should be given
       ! Parameters: (ave,sigma) of ln(x)
       n_par = 4
       Randomtype = 5
       RandomPara(1) = PCP%V_mean  ! mu of ln(x)
       RandomPara(2) = PCP%V_para  ! sigma of ln(x)
       RandomPara(3) = PCP%V_high  ! high end cutoff of x
       RandomPara(4) = PCP%V_para2 ! shift of x
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Speed = rnd_vec(1)
    Case(6)   ! Beta
       ! Parameters: (ave,a,b) ave not used
       n_par = 3
       Randomtype = 6
       RandomPara(1) = PCP%V_mean
       RandomPara(2) = PCP%V_para
       RandomPara(3) = PCP%V_para2
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Speed = rnd_vec(1)
    Case(7)   ! Triangular
       ! Parameters: (peak,min,max)
       n_par = 3
       Randomtype = 7
       RandomPara(1) = PCP%V_mean
       RandomPara(2) = PCP%V_low
       RandomPara(3) = PCP%V_high
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Speed = rnd_vec(1)
    Case(8)   ! Weibull  (alpha=1: Exponential)
       ! Parameters: (ave,alpha,lambda) ave not used
       n_par = 3
       Randomtype = 8
       RandomPara(1) = PCP%V_mean
       RandomPara(2) = PCP%V_para
       RandomPara(3) = PCP%V_para2
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Speed = rnd_vec(1)
    Case(9)   ! Gumbel
       ! Parameters: (ave,alpha) ave not used
       n_par = 2
       Randomtype = 9
       RandomPara(1) = PCP%V_mean
       RandomPara(2) = PCP%V_para
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Speed = rnd_vec(1)
    Case Default
       Call SHUTDOWN('ERROR: Class_Properties I_VEL_DIST')
    End Select

    Select Case(PCP%I_DIA_DIST)
    Case(-1)
       Call SHUTDOWN('ERROR: Class_Properties: -1')
    Case(0)
       HR%Radius  = 0.5_EB*PCP%D_mean
    Case(1)   ! Uniform
       ! Parameters: (ave,min,max) ave not used
       n_par = 3
       Randomtype = 1
       RandomPara(1) = PCP%D_mean
       RandomPara(2) = PCP%D_low
       RandomPara(3) = PCP%D_high
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Radius = 0.5_EB*rnd_vec(1)
    Case(2)   ! Truncated Normal
       n_par = 4
       Randomtype = 2
       RandomPara(1) = PCP%D_mean
       RandomPara(2) = PCP%D_para
       RandomPara(3) = PCP%D_low
       RandomPara(4) = PCP%D_high
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Radius = 0.5_EB*rnd_vec(1)
    Case(3)   ! Gamma
       ! Parameters: (ave,alpha,beta) ave not used
       n_par = 3
       Randomtype = 3
       RandomPara(1) = PCP%D_mean
       RandomPara(2) = PCP%D_para
       RandomPara(3) = PCP%D_para2
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Radius = 0.5_EB*rnd_vec(1)
    Case(4)   ! Normal
       ! Parameters: (ave,sigma)
       n_par = 2
       Randomtype = 4
       RandomPara(1) = PCP%D_mean
       RandomPara(2) = PCP%D_para
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Radius = 0.5_EB*rnd_vec(1)
    Case(5)   ! LogNormal
       ! mean and variance of log(x) should be given
       ! Parameters: (ave,sigma) of ln(x)
       n_par = 4
       Randomtype = 5
       RandomPara(1) = PCP%D_mean  ! mu of ln(x)
       RandomPara(2) = PCP%D_para  ! sigma of ln(x)
       RandomPara(3) = PCP%D_high  ! high end cutoff
       RandomPara(4) = PCP%D_para2 ! shift
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Radius = 0.5_EB*rnd_vec(1)
    Case(6)   ! Beta
       ! Parameters: (ave,a,b) ave not used
       n_par = 3
       Randomtype = 6
       RandomPara(1) = PCP%D_mean
       RandomPara(2) = PCP%D_para
       RandomPara(3) = PCP%D_para2
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Radius = 0.5_EB*rnd_vec(1)
    Case(7)   ! Triangular
       ! Parameters: (peak,min,max)
       n_par = 3
       Randomtype = 7
       RandomPara(1) = PCP%D_mean
       RandomPara(2) = PCP%D_low
       RandomPara(3) = PCP%D_high
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Radius = 0.5_EB*rnd_vec(1)
    Case(8)   ! Weibull  (alpha=1: Exponential)
       ! Parameters: (ave,alpha,lambda)
       n_par = 3
       Randomtype = 8
       RandomPara(1) = PCP%D_mean
       RandomPara(2) = PCP%D_para
       RandomPara(3) = PCP%D_para2
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Radius = 0.5_EB*rnd_vec(1)
    Case(9)   ! Gumbel
       ! Parameters: (ave,alpha)
       n_par = 2
       Randomtype = 9
       RandomPara(1) = PCP%D_mean
       RandomPara(2) = PCP%D_para
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Radius = 0.5_EB*rnd_vec(1)
    Case Default
       Call SHUTDOWN('ERROR: Class_Properties I_DIA_DIST')
    End Select
    HR%Mass   = 80.0_EB*(HR%Radius/0.27_EB)**2 

    Select Case(PCP%I_TAU_DIST)
    Case(-1)
       Call SHUTDOWN('ERROR: Class_Properties: -1')
    Case(0)
       HR%Tau  = PCP%Tau_mean
    Case(1)   ! Uniform
       ! Parameters: (ave,min,max) ave not used
       n_par = 3
       Randomtype = 1
       RandomPara(1) = 0.5_EB*(PCP%Tau_high+PCP%Tau_low)
       PCP%Tau_mean = 0.5_EB*(PCP%Tau_high+PCP%Tau_low)
       RandomPara(2) = PCP%Tau_low
       RandomPara(3) = PCP%Tau_high
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tau = rnd_vec(1)
    Case(2)   ! Truncated Normal
       ! Parameters: (ave,sigma,min,max)
       n_par = 4
       Randomtype = 2
       RandomPara(1) = PCP%Tau_mean
       RandomPara(2) = PCP%Tau_para
       RandomPara(3) = PCP%Tau_low
       RandomPara(4) = PCP%Tau_high
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tau = rnd_vec(1)
    Case(3)   ! Gamma
       ! Parameters: (ave,alpha,beta) ave not used
       n_par = 3
       Randomtype = 3
       RandomPara(1) = PCP%Tau_mean
       RandomPara(2) = PCP%Tau_para
       RandomPara(3) = PCP%Tau_para2
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tau = rnd_vec(1)
    Case(4)   ! Normal
       ! Parameters: (ave,sigma)
       n_par = 2
       Randomtype = 4
       RandomPara(1) = PCP%Tau_mean
       RandomPara(2) = PCP%Tau_para
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tau = rnd_vec(1)
    Case(5)   ! LogNormal
       ! mean and variance of log(x) should be given
       ! Parameters: (ave,sigma) of ln(x)
       n_par = 4
       Randomtype = 5
       RandomPara(1) = PCP%Tau_mean
       RandomPara(2) = PCP%Tau_para
       RandomPara(3) = PCP%Tau_high  ! high end cutoff
       RandomPara(4) = PCP%Tau_para2 ! shift
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tau = rnd_vec(1)
    Case(6)   ! Beta
       ! Parameters: (ave,a,b) ave not used
       n_par = 3
       Randomtype = 6
       RandomPara(1) = PCP%Tau_mean
       RandomPara(2) = PCP%Tau_para
       RandomPara(3) = PCP%Tau_para2
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tau = rnd_vec(1)
    Case(7)   ! Triangular
       ! Parameters: (peak,min,max)
       n_par = 3
       Randomtype = 7
       RandomPara(1) = PCP%Tau_mean
       RandomPara(2) = PCP%Tau_low
       RandomPara(3) = PCP%Tau_high
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tau = rnd_vec(1)
    Case(8)   ! Weibull  (alpha=1: Exponential)
       ! Parameters: (ave,alpha,lambda)
       n_par = 3
       Randomtype = 8
       RandomPara(1) = PCP%Tau_mean
       RandomPara(2) = PCP%Tau_para
       RandomPara(3) = PCP%Tau_para2
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tau = rnd_vec(1)
    Case(9)   ! Gumbel
       ! Parameters: (ave,alpha)
       n_par = 2
       Randomtype = 9
       RandomPara(1) = PCP%Tau_mean
       RandomPara(2) = PCP%Tau_para
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tau = rnd_vec(1)
    Case Default
       Call SHUTDOWN('ERROR: Class_Properties I_TAU_DIST')
    End Select

    Select Case(PCP%I_DET_DIST)
    Case(-1)
       Call SHUTDOWN('ERROR: Class_Properties: -1')
    Case(0)
       HR%Tdet  = PCP%Tdet_mean
    Case(1)   ! Uniform
       ! Parameters: (ave,min,max) ave not used
       n_par = 3
       Randomtype = 1
       RandomPara(1) = 0.5_EB*(PCP%Tdet_high+PCP%Tdet_low)
       RandomPara(2) = PCP%Tdet_low
       RandomPara(3) = PCP%Tdet_high
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tdet = rnd_vec(1)
    Case(2)   ! Truncated Normal
       ! Parameters: (ave,sigma,min,max)
       n_par = 4
       Randomtype = 2
       RandomPara(1) = PCP%Tdet_mean
       RandomPara(2) = PCP%Tdet_para
       RandomPara(3) = PCP%Tdet_low
       RandomPara(4) = PCP%Tdet_high
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tdet = rnd_vec(1)
    Case(3)   ! Gamma
       ! Parameters: (ave,alpha,beta) ave not used
       n_par = 3
       Randomtype = 3
       RandomPara(1) = PCP%Tdet_mean
       RandomPara(2) = PCP%Tdet_para
       RandomPara(3) = PCP%Tdet_para2
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tdet = rnd_vec(1)
    Case(4)   ! Normal
       ! Parameters: (ave,sigma)
       n_par = 2
       Randomtype = 4
       RandomPara(1) = PCP%Tdet_mean
       RandomPara(2) = PCP%Tdet_para
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tdet = rnd_vec(1)
    Case(5)   ! LogNormal
       ! mean and variance of log(x) should be given
       ! Parameters: (ave,sigma) of ln(x)
       n_par = 4
       Randomtype = 5
       RandomPara(1) = PCP%Tdet_mean
       RandomPara(2) = PCP%Tdet_para
       RandomPara(3) = PCP%Tdet_high  ! high end cutoff
       RandomPara(4) = PCP%Tdet_para2 ! shift
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tdet = rnd_vec(1)
    Case(6)   ! Beta
       ! Parameters: (ave,a,b) ave not used
       n_par = 3
       Randomtype = 6
       RandomPara(1) = PCP%Tdet_mean
       RandomPara(2) = PCP%Tdet_para
       RandomPara(3) = PCP%Tdet_para2
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tdet = rnd_vec(1)
    Case(7)   ! Triangular
       ! Parameters: (peak,min,max)
       n_par = 3
       Randomtype = 7
       RandomPara(1) = PCP%Tdet_mean
       RandomPara(2) = PCP%Tdet_low
       RandomPara(3) = PCP%Tdet_high
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tdet = rnd_vec(1)
    Case(8)   ! Weibull  (alpha=1: Exponential)
       ! Parameters: (ave,alpha,lambda)
       n_par = 3
       Randomtype = 8
       RandomPara(1) = PCP%Tdet_mean
       RandomPara(2) = PCP%Tdet_para
       RandomPara(3) = PCP%Tdet_para2
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tdet = rnd_vec(1)
    Case(9)   ! Gumbel
       ! Parameters: (ave,alpha)
       n_par = 2
       Randomtype = 9
       RandomPara(1) = PCP%Tdet_mean
       RandomPara(2) = PCP%Tdet_para
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tdet = rnd_vec(1)
    Case Default
       Call SHUTDOWN('ERROR: Class_Properties I_DET_DIST')
    End Select

    Select Case(PCP%I_PRE_DIST)
    Case(-1)
       Call SHUTDOWN('ERROR: Class_Properties: -1')
    Case(0)
       HR%Tpre  = Max(0._EB,PCP%Tpre_mean)
    Case(1)   ! Uniform
       ! Parameters: (ave,min,max) ave not used
       n_par = 3
       Randomtype = 1
       RandomPara(1) = 0.5_EB*(PCP%Tpre_high+PCP%Tpre_low)
       RandomPara(2) = PCP%Tpre_low
       RandomPara(3) = PCP%Tpre_high
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tpre = Max(0._EB,rnd_vec(1))
    Case(2)   ! Truncated Normal
       ! Parameters: (ave,sigma,min,max)
       n_par = 4
       Randomtype = 2
       RandomPara(1) = PCP%Tpre_mean
       RandomPara(2) = PCP%Tpre_para
       RandomPara(3) = PCP%Tpre_low
       RandomPara(4) = PCP%Tpre_high
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tpre = Max(0._EB,rnd_vec(1))
    Case(3)   ! Gamma
       ! Parameters: (ave,alpha,beta) ave not used
       n_par = 3
       Randomtype = 3
       RandomPara(1) = PCP%Tpre_mean
       RandomPara(2) = PCP%Tpre_para
       RandomPara(3) = PCP%Tpre_para2
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tpre = Max(0._EB,rnd_vec(1))
    Case(4)   ! Normal
       ! Parameters: (ave,sigma)
       n_par = 2
       Randomtype = 4
       RandomPara(1) = PCP%Tpre_mean
       RandomPara(2) = PCP%Tpre_para
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tpre = Max(0._EB,rnd_vec(1))
    Case(5)   ! LogNormal
       ! mean and variance of log(x) should be given
       ! Parameters: (ave,sigma) of ln(x)
       n_par = 4
       Randomtype = 5
       RandomPara(1) = PCP%Tpre_mean
       RandomPara(2) = PCP%Tpre_para
       RandomPara(3) = PCP%Tpre_high  ! high end cutoff
       RandomPara(4) = PCP%Tpre_para2 ! shift
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tpre = Max(0._EB,rnd_vec(1))
    Case(6)   ! Beta
       ! Parameters: (ave,a,b) ave not used
       n_par = 3
       Randomtype = 6
       RandomPara(1) = PCP%Tpre_mean
       RandomPara(2) = PCP%Tpre_para
       RandomPara(3) = PCP%Tpre_para2
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tpre = Max(0._EB,rnd_vec(1))
    Case(7)   ! Triangular
       ! Parameters: (peak,min,max)
       n_par = 3
       Randomtype = 7
       RandomPara(1) = PCP%Tpre_mean
       RandomPara(2) = PCP%Tpre_low
       RandomPara(3) = PCP%Tpre_high
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tpre = Max(0._EB,rnd_vec(1))
    Case(8)   ! Weibull  (alpha=1: Exponential)
       ! Parameters: (ave,alpha,lambda)
       n_par = 3
       Randomtype = 8
       RandomPara(1) = PCP%Tpre_mean
       RandomPara(2) = PCP%Tpre_para
       RandomPara(3) = PCP%Tpre_para2
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tpre = Max(0._EB,rnd_vec(1))
    Case(9)   ! Gumbel
       ! Parameters: (ave,alpha)
       n_par = 2
       Randomtype = 9
       RandomPara(1) = PCP%Tpre_mean
       RandomPara(2) = PCP%Tpre_para
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tpre = Max(0._EB,rnd_vec(1))
    Case Default
       Call SHUTDOWN('ERROR: Class_Properties I_PRE_DIST')
    End Select
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
  End Subroutine CLASS_PROPERTIES
!
  Subroutine RE_ALLOCATE_HUMANS(CODE,NM)
    Implicit None
    !
    ! Passed variables
    Integer, Intent(IN) :: CODE,NM
    !
    ! Local variables
    Integer IZERO
    Type (HUMAN_TYPE), Allocatable, Dimension(:) :: DUMMY
    Type (MESH_TYPE), Pointer :: M =>NULL()
    !
    If (.Not.Any(EVACUATION_GRID)) Return
    If ( .Not.(EVACUATION_ONLY(NM) .And. EVACUATION_GRID(NM)) ) Return

    Select Case(CODE)
       !
    Case(1)
       !
       M=>MESHES(NM)
       Allocate(DUMMY(1:M%N_HUMANS_DIM),STAT=IZERO)
       DUMMY = M%HUMAN
       !
       Deallocate(M%HUMAN)
       Allocate(M%HUMAN(M%N_HUMANS_DIM+1000),STAT=IZERO)
       M%HUMAN(1:M%N_HUMANS_DIM) = DUMMY(1:M%N_HUMANS_DIM)
       M%N_HUMANS_DIM = M%N_HUMANS_DIM+1000
       !
    End Select
    !
    Deallocate(DUMMY)
    !
  End Subroutine RE_ALLOCATE_HUMANS
!
  Subroutine DUMP_EVAC(T,NM)
    Implicit None
    !
    ! Passed variables
    Integer, Intent(IN) :: NM
    Real(EB), Intent(IN) :: T
    !
    ! Local variables
    Integer :: NPP,NPLIM,i,izero,nn,n
    Real(EB) :: TNOW, EVEL, angle_hr
    Real(FB), Allocatable, Dimension(:) :: XP,YP,ZP
    Real(FB), Allocatable, Dimension(:,:) :: QP, AP
    Integer, Allocatable, Dimension(:) :: TA
    Type (HUMAN_TYPE), Pointer :: HR =>NULL()
    !
    TNOW=SECOND() 
    !
    If (.Not.Any(EVACUATION_GRID)) Return
    If (.Not.(EVACUATION_ONLY(NM) .And. EVACUATION_GRID(NM))) Return
    Call POINT_TO_MESH(NM)
    !       Call POINT_TO_EVAC_MESH(NM)

    ! Write the current time to the prt5 file, then start
    ! looping through the particle classes
    Write(LU_PART(NM)) Real(T,FB)

    HUMAN_CLASS_LOOP: Do N=1,N_EVAC
       ! Count the number of humans to dump out
       !
       NPLIM = 0
       Do I=1,N_HUMANS
          HR=>HUMAN(I)
          If (HR%COLOR_INDEX < 1) HR%COLOR_INDEX = 1
          If (HR%COLOR_INDEX > EVAC_AVATAR_NCOLOR) HR%COLOR_INDEX = EVAC_AVATAR_NCOLOR
          If (HR%SHOW .And. N_EVAC == 1) NPLIM = NPLIM + 1
       End Do
       NPLIM = Min(NPPS,NPLIM)
       !
       Allocate(TA(NPLIM),STAT=IZERO)
       Call ChkMemErr('DUMP_EVAC','TA',IZERO)
       Allocate(XP(NPLIM),STAT=IZERO)
       Call ChkMemErr('DUMP_EVAC','XP',IZERO)
       Allocate(YP(NPLIM),STAT=IZERO)
       Call ChkMemErr('DUMP_EVAC','YP',IZERO)
       Allocate(ZP(NPLIM),STAT=IZERO)
       Call ChkMemErr('DUMP_EVAC','ZP',IZERO)
       ! body angle, semi major axis, semi minor axis
       Allocate(AP(NPLIM,4),STAT=IZERO)
       Call ChkMemErr('DUMP_EVAC','AP',IZERO)
       If (EVAC_N_QUANTITIES > 0) Then
          Allocate(QP(NPLIM,EVAC_N_QUANTITIES),STAT=IZERO)
          Call ChkMemErr('DUMP_EVAC','QP',IZERO)
       End If
       !
       ! Load human coordinates into single precision array
       !
       NPP = 0
       PLOOP: Do I=1,N_HUMANS
          HR=>HUMAN(I)
          If (.Not. HR%SHOW ) Cycle PLOOP
          NPP = NPP + 1
          TA(NPP) = HR%ILABEL
          XP(NPP) = Real(HR%X,FB)
          YP(NPP) = Real(HR%Y,FB)
          ZP(NPP) = Min( Max(Real(HR%Z,FB), EVAC_Z_MIN), EVAC_Z_MAX)

          AP(NPP,1) = 180.0_FB*Real(HR%Angle/Pi,FB)
          AP(NPP,2) =   2.0_FB*Real(HR%Radius,FB)
          AP(NPP,3) =   2.0_FB*Real(HR%r_torso,FB)
          ! Height of a human scaled by radius, default male 1.80 m
          AP(NPP,4) =  1.80_FB*Real(HR%Radius/0.27_EB,FB)

          Do NN=1,EVAC_N_QUANTITIES
             Select Case(EVAC_QUANTITIES_INDEX(NN))
             Case(240)  ! MOTIVE_ACCELERATION, Unimpeded walking Speed / tau
                EVEL = Sqrt(HR%UBAR**2 + HR%VBAR**2)
                If (EVEL > 0.0_EB) Then
                   QP(NPP,NN) = Real(HR%Speed/HR%Tau,FB)
                Else
                   QP(NPP,NN) = Real(EVEL,FB)
                End If
             Case(241)  ! FED_DOSE, Fractional Effective Dose
                QP(NPP,NN) = Real(HR%IntDose,FB)
             Case(242)  ! SPEED, Human speed
                QP(NPP,NN) = Real(Sqrt(HR%U**2 + HR%V**2),FB)
             Case(243)  ! ANGULAR_SPEED, Human Angular Velocity
                QP(NPP,NN) = Real(HR%Omega,FB)
             Case(244)  ! ACCELERATION, Human acc.
                QP(NPP,NN) = Real(Sqrt( ((HR%UBAR-HR%U)/HR%Tau + (HR%F_X/HR%Mass))**2 + &
                     ((HR%VBAR-HR%V)/HR%Tau + (HR%F_Y/HR%Mass))**2 ),FB)
             Case(245)  ! CONTACT_LINEFORCE, Human Pressure: contact forces
                QP(NPP,NN) = Real(HR%SumForces ,FB)
             Case(246)  ! TOTAL_LINEFORCE, Human Pressure2: contact + social
                QP(NPP,NN) = Real(HR%SumForces2 ,FB)
             Case(247)  ! COLOR, Human color index
                QP(NPP,NN) = Real(HR%COLOR_INDEX - 1,FB)
             Case(248)  ! MOTIVE_ANGLE, 
                EVEL = Sqrt(HR%UBAR**2 + HR%VBAR**2)
                If (EVEL > 0.0_EB) Then
                   If (HR%VBAR >= 0.0_EB) Then
                      angle_hr = Acos(HR%UBAR/EVEL)
                   Else
                      angle_hr = 2.0_EB*Pi - Acos(HR%UBAR/EVEL)
                   End If
                   If (angle_hr == 2.0_EB*Pi) angle_hr = 0.0_EB  ! agent HR angle is [0,2Pi)
                Else
                   angle_hr = 0.0_EB
                End If
                QP(NPP,NN) = Real(angle_hr*180.0_EB/Pi,FB)
             End Select
          End Do

          If (NPP>=NPPS) Exit PLOOP
       End Do PLOOP
       !
       ! Dump human data into the .prt5 file
       !
       Write(LU_PART(NM)) NPLIM
       Write(LU_PART(NM)) (XP(I),I=1,NPLIM),(YP(I),I=1,NPLIM),(ZP(I),I=1,NPLIM), &
            (AP(I,1),I=1,NPLIM),(AP(I,2),I=1,NPLIM),(AP(I,3),I=1,NPLIM),(AP(I,4),I=1,NPLIM)
       Write(LU_PART(NM)) (TA(I),I=1,NPLIM)
       If (EVAC_N_QUANTITIES > 0) Then
          Write(LU_PART(NM)) ((QP(I,NN),I=1,NPLIM),NN=1,EVAC_N_QUANTITIES)
       End If
       !
       If (EVAC_N_QUANTITIES > 0) Then
          Deallocate(QP)
       End If
       Deallocate(AP)
       Deallocate(ZP)
       Deallocate(YP)
       Deallocate(XP)
       Deallocate(TA)

    End Do HUMAN_CLASS_LOOP

    !
    TUSED(7,NM) = TUSED(7,NM) + SECOND() - TNOW
  End Subroutine DUMP_EVAC
!
  Function GaussRand( gmean, gtheta, gcutmult )
    Implicit None
    !
    ! Random numbers from the Gaussian distribution
    !
    Real(EB) GaussRand
    !  generates a random number (x) with
    !  P(x) = exp[- (x-gmean)^2 / (2*gtheta)], if x is in 
    !         [gmean - gcutmult*sqrt(gtheta), gmean + gcutmult*sqrt(gtheta)]
    !       = 0                              , if not 
    !  
    ! Passed variables
    Real(EB) gmean, gtheta, gcutmult

    ! Local variables
    Real(EB) v1,v2,rsq,fac, rn

    If ( (GaussFlag == 1)  .And. (Abs(GaussSet2-gmean) <= gcutmult*Sqrt(gtheta)) ) Then
       GaussFlag = 0
       GaussRand = GaussSet2
    Else
       GaussFlag = 0
       Do 
          Call Random_number(rn)
          v1 = 1.0_EB - 2.0_EB*rn
          Call Random_number(rn)
          v2 = 1.0_EB - 2.0_EB*rn
          rsq = v1*v1 + v2*v2
          Do While (rsq >= 1.0 )
             Call Random_number(rn)
             v1 = 1.0_EB - 2.0_EB*rn
             Call Random_number(rn)
             v2 = 1.0_EB - 2.0_EB*rn
             rsq=v1*v1+v2*v2
          End Do
          fac = Sqrt(-2.0_EB*gtheta*Log(rsq)/rsq)
          GaussSet1 = v1*fac + gmean
          GaussSet2 = v2*fac + gmean

          If ( (Abs(GaussSet1-gmean) <= gcutmult*Sqrt(gtheta)) .Or. (Abs(GaussSet2-gmean) <= gcutmult*Sqrt(gtheta)) ) Then
             Exit
          End If

       End Do

       If (Abs(GaussSet1-gmean) <= gcutmult*Sqrt(gtheta)) Then
          GaussFlag = 1
          GaussRand = GaussSet1
       Else
          GaussFlag = 0
          GaussRand = GaussSet2
       End If
    End If

  End Function GaussRand
!
  Function GaussTrun( gmean, gsigma, glow, ghigh )
    Implicit None
    !
    ! Random numbers from the Gaussian distribution
    !
    Real(EB) GaussTrun
    !  generates a random number (x) with
    !  P(x) = exp[- (x-gmean)^2 / (2*gsigma^2)], if x is in 
    !         [gmean - gcutmult*gsigma, gmean + gcutmult*gsigma]
    !       = 0                              , if not 
    !  
    ! Passed variables
    Real(EB) gmean, gsigma, glow, ghigh

    ! Local variables
    Real(EB) v1, v2, rsq, fac, rn

    If ( (GTrunFlag == 1)  .And. (GTrunSet2 >= glow) .And. (GTrunSet2 <= ghigh) ) Then
       GTrunFlag = 0
       GaussTrun = GTrunSet2
    Else
       GTrunFlag = 0
       Do 
          Call Random_number(rn)
          v1 = 1.0_EB - 2.0_EB*rn
          Call Random_number(rn)
          v2 = 1.0_EB - 2.0_EB*rn
          rsq = v1*v1 + v2*v2
          Do While (rsq >= 1.0 )
             Call Random_number(rn)
             v1 = 1.0_EB - 2.0_EB*rn
             Call Random_number(rn)
             v2 = 1.0_EB - 2.0_EB*rn
             rsq = v1*v1 + v2*v2
          End Do
          fac = Sqrt(-2.0_EB*(gsigma**2)*Log(rsq)/rsq)
          GTrunSet1 = v1*fac + gmean
          GTrunSet2 = v2*fac + gmean

          If ( ((GTrunSet1 >= glow) .And. (GTrunSet1 <= ghigh)) .Or. ((GTrunSet2 >= glow) .And. (GTrunSet2 <= ghigh)) ) Then
             Exit
          End If

       End Do

       If ( (GTrunSet1 >= glow) .And. (GTrunSet1 <= ghigh) ) Then
          GTrunFlag = 1
          GaussTrun = GTrunSet1
       Else
          GTrunFlag = 0
          GaussTrun = GTrunSet2
       End If
    End If

  End Function GaussTrun
!
  Subroutine DUMP_EVAC_CSV(Tin)
    Implicit None
    !
    ! Dump human data to CHID_evac.csv
    !
    ! Passed variables
    Real(EB), Intent(IN) :: Tin
    !
    ! Local variables
    Character(50) tcform
    Integer n_cols, n_tot_humans, i, ii, izero
    Integer, Allocatable, Dimension(:) :: ITEMP
    !
    If (.Not.Any(EVACUATION_GRID)) Return
    !
    Allocate(ITEMP(Max(1,n_exits-n_co_exits+n_doors)), STAT = IZERO)
    Call ChkMemErr('DUMP_EVAC_CSV','ITEMP', IZERO)
    !
    ! Output first the floors then the corridors
    n_cols = n_egrids + n_corrs + n_exits + n_doors + 1 + n_exits - n_co_exits + n_doors
    n_tot_humans = 0
    Do i = 1, n_egrids
       n_tot_humans = n_tot_humans + MESHES(EVAC_Node_List(i)%IMESH)%N_HUMANS
    End Do
    Do i = 1, n_corrs
       n_tot_humans = n_tot_humans + EVAC_CORRS(i)%n_inside
    End Do
    !
    ii = 0
    Do i = 1, n_exits
       If (.Not.EVAC_EXITS(i)%COUNT_ONLY) Then
          ii = ii + 1
          ITEMP(ii) = EVAC_EXITS(i)%NTARGET(50)
       End If
    End Do
    Do i = 1, n_doors
       ii = ii + 1
       ITEMP(ii) = EVAC_DOORS(i)%NTARGET(50)
    End Do
    !
    If (n_dead >= 0) Then
       ! Write the 'fed' columns
       Write(tcform,'(a,i4.4,a,a)') "(ES13.5E3,",n_cols+1, &
            "(',',i8)", ",',',ES13.5E3,',',ES13.5E3)"
       Write (LU_EVACCSV,fmt=tcform) Tin, n_tot_humans, &
            (MESHES(EVAC_Node_List(i)%IMESH)%N_HUMANS, i=1,n_egrids), &
            (EVAC_CORRS(i)%n_inside, i = 1,n_corrs), &
            (EVAC_EXITS(i)%ICOUNT, i = 1,n_exits), &
            (EVAC_DOORS(i)%ICOUNT, i = 1,n_doors), &
            (ITEMP(i), i = 1,n_exits-n_co_exits+n_doors), &
            n_dead, fed_max, fed_max_alive
    Else
       ! Do not write the 'fed' columns
       Write(tcform,'(a,i4.4,a)') "(ES13.5E3,",n_cols, "(',',i8),i8)"
       Write (LU_EVACCSV,fmt=tcform) Tin, n_tot_humans, &
            (MESHES(EVAC_Node_List(i)%IMESH)%N_HUMANS, i=1,n_egrids), &
            (EVAC_CORRS(i)%n_inside, i = 1,n_corrs), &
            (EVAC_EXITS(i)%ICOUNT, i = 1,n_exits), &
            (EVAC_DOORS(i)%ICOUNT, i = 1,n_doors), &
            (ITEMP(i), i = 1,n_exits-n_co_exits+n_doors)
    End If
    Deallocate(ITEMP)
    !
  End Subroutine DUMP_EVAC_CSV

  Logical Function See_each_other(nm, r1_x, r1_y, r2_x, r2_y)
    ! This function returns true, if the two points have a line-of-sight.
    ! This function does not use smoke information, i.e., it just sees if
    ! there are obstacles between the two points.
    ! Inputs:  nm: mesh index, r1 an r2 should belong to the same mesh
    !          (r1_x,r1_y): co-ordinates of the first agent
    !          (r2_x,r2_y): co-ordinates of the second agent
    ! NOTE: This works for thick and thin OBSTs.
    !
    ! Passed variables
    Integer, Intent(IN) :: nm
    Real(EB), Intent(IN) :: r1_x, r1_y, r2_x, r2_y
    !
    ! Local variables
    Integer :: i, j, isx, isy, i_r1, i_r2, j_r1, j_r2
    Integer :: i_old, j_old, ic, ic2, iw, iw1, iw2
    Real(EB) :: x, y
    Type (MESH_TYPE), Pointer :: M =>NULL()

    M => MESHES(NM)
    See_each_other = .True.  ! Default

    isx = Int(Sign(1.0_EB,r2_x - r1_x)) ! loop increment +1 or -1
    isy = Int(Sign(1.0_EB,r2_y - r1_y)) ! loop increment +1 or -1
    i_r1 = Floor(M%CELLSI(Floor((r1_x-M%XS)*M%RDXINT)) + 1.0_EB) ! II start
    i_r2 = Floor(M%CELLSI(Floor((r2_x-M%XS)*M%RDXINT)) + 1.0_EB) ! II end
    j_r1 = Floor(M%CELLSJ(Floor((r1_y-M%YS)*M%RDYINT)) + 1.0_EB) ! JJ start
    j_r2 = Floor(M%CELLSJ(Floor((r2_y-M%YS)*M%RDYINT)) + 1.0_EB) ! JJ end
    i_r1 = Max(1,Min(i_r1,M%IBAR)) ! To be sure that indices are always ok.
    i_r2 = Max(1,Min(i_r2,M%IBAR))
    j_r1 = Max(1,Min(j_r1,M%JBAR))
    j_r2 = Max(1,Min(j_r2,M%JBAR))

    ! Same cell: sees always each other
    If (Abs(i_r2-i_r1)+Abs(j_r2-j_r1) .Lt. 1) Return

    ! Choose the main direction:
    If (Abs(i_r2-i_r1) .Lt. Abs(j_r2-j_r1)) Then
       ! Now y is the main direction
       i_old = i_r1 ; j_old = j_r1
       y = 0.0_EB
       MainLoopY: Do j = j_r1+isy, j_r2, isy
          y = y + isy*M%DY(j)
          x = Max(M%XS,Min(M%XF,r1_x + y*(r2_x - r1_x)/(r2_y - r1_y)))
          i = Floor(M%CELLSI(Floor((x-M%XS)*M%RDXINT)) + 1.0_EB)
          i = isx*Min(isx*i_r2,isx*i) ! i in interval i_r1...i_r2
          ic  = M%CELL_INDEX(i_old,j_old,1)
          ic2 = M%CELL_INDEX(i    ,j_old,1) ! side cell
          iw  = M%WALL_INDEX(ic, isy*2) ! main direction
          iw1 = M%WALL_INDEX(ic ,isx*1) ! sideways
          iw2 = M%WALL_INDEX(ic2,isy*2) ! side + main direction
          ! iw is zero, if there is no solid boundary
          ! from (i,j)==>(i,jnew):    iw and iw2 are zero, iw1 does not matter
          !                           ic=ic2 ==> iw=iw2
          ! from (i,j)==>(inew,jnew): iw1 and iw2 are zero, iw does not matter
          If ((i==i_old .And. iw/=0) .Or. (i/=i_old .And. (iw1/=0 .Or. iw2/=0))) Then
             See_each_other = .False.
             Exit MainLoopY
          End If
          i_old = i ; j_old = j
       End Do MainLoopY
    Else
       ! Now x is the main direction
       i_old = i_r1 ; j_old = j_r1
       x = 0.0_EB 
       MainLoopX: Do i = i_r1+isx, i_r2, isx
          x = x + isx*M%DX(i)
          y = Max(M%YS,Min(M%YF,r1_y + x*(r2_y - r1_y)/(r2_x - r1_x)))
          j = Floor(M%CELLSJ(Floor((y-M%YS)*M%RDYINT)) + 1.0_EB)
          j = isy*Min(isy*j_r2,isy*j) ! j in interval j_r1...j_r2
          ic  = M%CELL_INDEX(i_old,j_old,1)
          ic2 = M%CELL_INDEX(i_old,j    ,1) ! side cell
          iw  = M%WALL_INDEX(ic, isx*1) ! main direction
          iw1 = M%WALL_INDEX(ic ,isy*2) ! sideways
          iw2 = M%WALL_INDEX(ic2,isx*1) ! side + main direction
          ! iw is zero, if there is no solid boundary
          ! from (i,j)==>(inew,j):    iw and iw2 are zero, iw1 does not matter
          !                           ic=ic2 ==> iw=iw2
          ! from (i,j)==>(inew,jnew): iw1 and iw2 are zero, iw does not matter
          If ((j==j_old .And. iw/=0) .Or. (j/=j_old .And. (iw1/=0 .Or. iw2/=0))) Then
             See_each_other = .False.
             Exit MainLoopX
          End If
          i_old = i ; j_old = j
       End Do MainLoopX
    End If
  End Function See_each_other
  !
  Logical Function See_door(nm, idoor, itarget, r1_x, r1_y, r2_x, r2_y, ave_K, max_fed)
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
    Integer, Intent(IN) :: nm, idoor, itarget
    Real(EB), Intent(IN) :: r1_x, r1_y, r2_x, r2_y
    Real(EB), Intent(OUT) :: ave_K, max_fed
    !
    ! Local variables
    Integer :: i, j, isx, isy, i_r1, i_r2, j_r1, j_r2
    Integer :: i_old, j_old, ic, ic2, iw, iw1, iw2
    Real(EB) :: x, y
    Type (MESH_TYPE), Pointer :: M =>NULL()

    M => MESHES(NM)
    See_door = .True.  ! Default
    ave_K = 0.0_EB
    max_fed = 0.0_EB

    isx = Int(Sign(1.0_EB,r2_x - r1_x)) ! loop increment +1 or -1
    isy = Int(Sign(1.0_EB,r2_y - r1_y)) ! loop increment +1 or -1
    i_r1 = Floor(M%CELLSI(Floor((r1_x-M%XS)*M%RDXINT)) + 1.0_EB) ! II start
    i_r2 = Floor(M%CELLSI(Floor((r2_x-M%XS)*M%RDXINT)) + 1.0_EB) ! II end
    j_r1 = Floor(M%CELLSJ(Floor((r1_y-M%YS)*M%RDYINT)) + 1.0_EB) ! JJ start
    j_r2 = Floor(M%CELLSJ(Floor((r2_y-M%YS)*M%RDYINT)) + 1.0_EB) ! JJ end
    i_r1 = Max(1,Min(i_r1,M%IBAR)) ! To be sure that indices are always ok.
    i_r2 = Max(1,Min(i_r2,M%IBAR))
    j_r1 = Max(1,Min(j_r1,M%JBAR))
    j_r2 = Max(1,Min(j_r2,M%JBAR))

    ! Same cell: sees always each other
    If (Abs(i_r2-i_r1)+Abs(j_r2-j_r1) .Lt. 1) Then
       ave_K = MASS_EXTINCTION_COEFFICIENT*1.0E-6_EB*M%HUMAN_GRID(i_r1,j_r1)%SOOT_DENS
       max_fed = M%HUMAN_GRID(i_r1,j_r1)%FED_CO_CO2_O2
       Return
    End If

    ! Choose the main direction:
    If (Abs(i_r2-i_r1) .Lt. Abs(j_r2-j_r1)) Then
       ! Now y is the main direction

       ! Add the last cell (r1) to the average
       ave_K = Mass_extinction_coefficient*1.0E-6_EB*M%HUMAN_GRID(i_r1,j_r1)%SOOT_DENS/(Abs(j_r1-j_r2)+1)
       max_fed =  M%HUMAN_GRID(i_r1,j_r1)%FED_CO_CO2_O2

       i_old = i_r1 ; j_old = j_r1
       y = 0.0_EB
       MainLoopY: Do j = j_r1+isy, j_r2, isy
          y = y + isy*M%DY(j)
          x = Max(M%XS,Min(M%XF,r1_x + y*(r2_x - r1_x)/(r2_y - r1_y)))
          i = Floor(M%CELLSI(Floor((x-M%XS)*M%RDXINT)) + 1.0_EB)
          i = isx*Min(isx*i_r2,isx*i) ! i in interval j_r1...j_r2
          ic  = M%CELL_INDEX(i_old,j_old,1)
          ic2 = M%CELL_INDEX(i    ,j_old,1) ! side cell
          iw  = M%WALL_INDEX(ic, isy*2) ! main direction
          iw1 = M%WALL_INDEX(ic ,isx*1) ! sideways
          iw2 = M%WALL_INDEX(ic2,isy*2) ! side + main direction
          ! iw is zero, if there is no solid boundary
          ! from (i,j)==>(i,jnew):    iw and iw2 are zero, iw1 does not matter
          !                           ic=ic2 ==> iw=iw2
          ! from (i,j)==>(inew,jnew): iw1 and iw2 are zero, iw does not matter
          If ((i==i_old .And. iw/=0) .Or. (i/=i_old .And. (iw1/=0 .Or. iw2/=0)) .And. (idoor .Ne. itarget)) Then
             See_door = .False.
             ave_K = MASS_EXTINCTION_COEFFICIENT*1.0E-6_EB*M%HUMAN_GRID(i_r1,j_r1)%SOOT_DENS
             max_fed = M%HUMAN_GRID(i_r1,j_r1)%FED_CO_CO2_O2
             Exit MainLoopY
          End If
          ave_K = ave_K + MASS_EXTINCTION_COEFFICIENT*1.0E-6_EB*M%HUMAN_GRID(i,j)%SOOT_DENS/(Abs(j_r1-j_r2)+1)
          max_fed = Max(max_fed, M%HUMAN_GRID(i,j)%FED_CO_CO2_O2)
          i_old = i ; j_old = j
       End Do MainLoopY

    Else
       ! Now x is the main direction
       ! Add the first cell (r1) to the average
       ave_K = Mass_extinction_coefficient*1.0E-6_EB*M%HUMAN_GRID(i_r1,j_r1)%SOOT_DENS/(Abs(i_r1-i_r2)+1)
       max_fed =  M%HUMAN_GRID(i_r1,j_r1)%FED_CO_CO2_O2
       i_old = i_r1 ; j_old = j_r1
       x = 0.0_EB 
       MainLoopX: Do i = i_r1+isx, i_r2, isx
          x = x + isx*M%DX(i)
          y = Max(M%YS,Min(M%YF,r1_y + x*(r2_y - r1_y)/(r2_x - r1_x)))
          j = Floor(M%CELLSJ(Floor((y-M%YS)*M%RDYINT)) + 1.0_EB)
          j = isy*Min(isy*j_r2,isy*j) ! j in interval j_r1...j_r2
          ic  = M%CELL_INDEX(i_old,j_old,1)
          ic2 = M%CELL_INDEX(i_old,j    ,1) ! side cell
          iw  = M%WALL_INDEX(ic, isx*1) ! main direction
          iw1 = M%WALL_INDEX(ic ,isy*2) ! sideways
          iw2 = M%WALL_INDEX(ic2,isx*1) ! side + main direction
          ! iw is zero, if there is no solid boundary
          ! from (i,j)==>(inew,j):    iw and iw2 are zero, iw1 does not matter
          !                           ic=ic2 ==> iw=iw2
          ! from (i,j)==>(inew,jnew): iw1 and iw2 are zero, iw does not matter
          If ((i==i_old .And. iw/=0) .Or. (i/=i_old .And. (iw1/=0 .Or. iw2/=0)) .And. (idoor .Ne. itarget)) Then
             See_door = .False.
             ave_K = MASS_EXTINCTION_COEFFICIENT*1.0E-6_EB*M%HUMAN_GRID(i_r1,j_r1)%SOOT_DENS
             max_fed = M%HUMAN_GRID(i_r1,j_r1)%FED_CO_CO2_O2
             Exit MainLoopX
          End If
          ave_K = ave_K + MASS_EXTINCTION_COEFFICIENT*1.0E-6_EB*M%HUMAN_GRID(i,j)%SOOT_DENS/(Abs(i_r1-i_r2)+1)
          max_fed = Max(max_fed, M%HUMAN_GRID(i,j)%FED_CO_CO2_O2)
          i_old = i ; j_old = j
       End Do MainLoopX
    End If
  End Function See_door

  Subroutine Find_walls(nm, r1_x, r1_y, r_circle, d_cutoff, Skip_Wall_Force_Ior, d_xy, FoundWall_xy, istat)
    Implicit None
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
    Integer, Intent(IN) :: nm, Skip_Wall_Force_Ior
    Real(EB), Intent(IN) :: r1_x, r1_y, r_circle, d_cutoff
    Real(EB), Dimension(4), Intent(OUT) :: d_xy
    Logical, Dimension(4), Intent(OUT) :: FoundWall_xy
    Integer, Intent(OUT) :: istat
    !
    ! Local variables
    Integer :: ii, jj, iw, ic, ibc, is, i_end, iin, jjn, kkn, I_OBST
    Real(EB) :: dx, dy, d_mx, d_px, d_my, d_py
    Type (MESH_TYPE), Pointer :: M =>NULL()

    M => MESHES(NM)
    istat = 0            ! Default
    FoundWall_xy = .False.
    d_mx = M%XS
    d_px = M%XF
    d_my = M%YS
    d_py = M%YF

    ! Where is the circle
    IIN = Floor(M%CELLSI(Floor((r1_x-M%XS)*M%RDXINT))+1.0_EB)
    JJN = Floor(M%CELLSJ(Floor((r1_y-M%YS)*M%RDYINT))+1.0_EB)
    KKN = 1

    ! Find the closest wall at the -x direction

    is = -1    ! minus or plus direction
    If (Skip_Wall_Force_Ior==-1) Then
       d_mx = d_mx + is*(2.0_EB*d_cutoff)
    Else   
       dx = M%DX(iin) ! no stretched meshes for evac
       i_end = Int((d_cutoff+r_circle)/dx) + 1
       i_end = iin + is*i_end
       i_end = Min(M%IBAR,Max(1,i_end))
       ii = iin
       ic = M%cell_index(ii,jjn,kkn)  ! cell index
       I_OBST = M%OBST_INDEX_C(IC)
       If (M%Solid(ic) .AND. .NOT.M%OBSTRUCTION(I_OBST)%HIDDEN) Then
          istat = -1
          Return
       End If
       iw = M%wall_index(ic, is*1)      ! wall index
       Do While (iw==0 .And. ii/=i_end)
          ii = ii + is
          ic = M%cell_index(ii,jjn,kkn)  ! cell index
          iw = M%wall_index(ic, is*1)      ! wall index
       End Do

       If (iw /= 0) Then
          FoundWall_xy(1) = .True.
          ibc = M%IJKW(5,iw)         ! Boundary condition index
          ! There is a 'door', i.e., outflow-boundary (or open boundary)
          ! so no wall forces ==> exit this loop
          d_mx = M%xw(iw)
          I_OBST = M%OBST_INDEX_C(IC)
          If (M%Solid(ic) .AND. .NOT.M%OBSTRUCTION(I_OBST)%HIDDEN) Then
             Write(MESSAGE,'(A,I4,2I6)') 'ERROR: Find_Walls ',nm, ii,jjn
             Call SHUTDOWN(MESSAGE)
          End If
          If (SURFACE(ibc)%VEL> 0.0_EB .Or. M%BOUNDARY_TYPE(iw)==OPEN_BOUNDARY) Then
             !d_mx = d_mx + is*(2.0_EB*d_cutoff)
             FoundWall_xy(1) = .False.
          End If
          If (Abs(d_mx-r1_x) < r_circle) istat = 1
       End If
    End If

    ! Find the closest wall at the +x direction

    is = +1    ! minus or plus direction
    If (Skip_Wall_Force_Ior==1) Then
       d_px = d_px + is*(2.0_EB*d_cutoff)
    Else   
       dx = M%DX(iin) ! no stretched meshes for evac
       i_end = Int((d_cutoff+r_circle)/dx) + 1
       i_end = iin + is*i_end
       i_end = Min(M%IBAR,Max(1,i_end))
       ii = iin
       ic = M%cell_index(ii,jjn,kkn)  ! cell index
       I_OBST = M%OBST_INDEX_C(IC)
       If (M%Solid(ic) .AND. .NOT.M%OBSTRUCTION(I_OBST)%HIDDEN) Then
          istat = -1
          Return
       End If
       iw = M%wall_index(ic, is*1)      ! wall index
       Do While (iw==0 .And. ii/=i_end)
          ii = ii + is
          ic = M%cell_index(ii,jjn,kkn)  ! cell index
          iw = M%wall_index(ic, is*1)      ! wall index
       End Do

       If (iw /= 0) Then
          FoundWall_xy(2) = .True.
          ibc = M%IJKW(5,iw)         ! Boundary condition index
          ! There is a 'door', i.e., outflow-boundary (or open boundary)
          ! so no wall forces ==> exit this loop
          d_px = M%xw(iw)
          I_OBST = M%OBST_INDEX_C(IC)
          If (M%Solid(ic) .AND. .NOT.M%OBSTRUCTION(I_OBST)%HIDDEN) Then
             Write(MESSAGE,'(A,I4,2I6)') 'ERROR: Find_Walls ',nm, ii,jjn
             Call SHUTDOWN(MESSAGE)
          End If
          If (SURFACE(ibc)%VEL> 0.0_EB .Or. M%BOUNDARY_TYPE(iw)==OPEN_BOUNDARY) Then
             !d_px = d_px + is*(2.0_EB*d_cutoff)
             FoundWall_xy(2) = .False.
          End If
          If (Abs(d_px-r1_x) < r_circle) istat = 1
       End If
    End If

    ! Find the closest wall at the -y direction

    is = -1    ! minus or plus direction
    If (Skip_Wall_Force_Ior==-2) Then
       d_my = d_my + is*(2.0_EB*d_cutoff)
    Else   
       dy = M%DY(jjn) ! no stretched meshes for evac
       i_end = Int((d_cutoff+r_circle)/dy) + 1
       i_end = jjn + is*i_end
       i_end = Min(M%JBAR,Max(1,i_end))
       jj = jjn
       ic = M%cell_index(iin,jj,kkn)  ! cell index
       I_OBST = M%OBST_INDEX_C(IC)
       If (M%Solid(ic) .AND. .NOT.M%OBSTRUCTION(I_OBST)%HIDDEN) Then
          istat = -1
          Return
       End If
       iw = M%wall_index(ic, is*2)      ! wall index
       Do While (iw==0 .And. jj/=i_end)
          jj = jj + is
          ic = M%cell_index(iin,jj,kkn)  ! cell index
          iw = M%wall_index(ic, is*2)      ! wall index
       End Do

       If (iw /= 0) Then
          FoundWall_xy(3) = .True.
          ibc = M%IJKW(5,iw)         ! Boundary condition index
          ! There is a 'door', i.e., outflow-boundary (or open boundary)
          ! so no wall forces ==> exit this loop
          d_my = M%yw(iw)
          I_OBST = M%OBST_INDEX_C(IC)
          If (M%Solid(ic) .AND. .NOT.M%OBSTRUCTION(I_OBST)%HIDDEN) Then
             Write(MESSAGE,'(A,I4,2I6)') 'ERROR: Find_Walls ',nm, ii,jjn
             Call SHUTDOWN(MESSAGE)
          End If
          If (SURFACE(ibc)%VEL> 0.0_EB .Or. M%BOUNDARY_TYPE(iw)==OPEN_BOUNDARY) Then
             !d_my = d_my + is*(2.0_EB*d_cutoff)
             FoundWall_xy(3) = .False.
          End If
          If (Abs(d_my-r1_y) < r_circle) istat = 1
       End If
    End If

    ! Find the closest wall at the +y direction

    is = +1    ! minus or plus direction
    If (Skip_Wall_Force_Ior==2) Then
       d_py = d_py + is*(2.0_EB*d_cutoff)
    Else   
       dy = M%DY(jjn) ! no stretched meshes for evac
       i_end = Int((d_cutoff+r_circle)/dy) + 1
       i_end = jjn + is*i_end
       i_end = Min(M%JBAR,Max(1,i_end))
       jj = jjn
       ic = M%cell_index(iin,jj,kkn)  ! cell index
       I_OBST = M%OBST_INDEX_C(IC)
       If (M%Solid(ic) .AND. .NOT.M%OBSTRUCTION(I_OBST)%HIDDEN) Then
          istat = -1
          Return
       End If
       iw = M%wall_index(ic, is*2)      ! wall index
       Do While (iw==0 .And. jj/=i_end)
          jj = jj + is
          ic = M%cell_index(iin,jj,kkn)  ! cell index
          iw = M%wall_index(ic, is*2)      ! wall index
       End Do

       If (iw /= 0) Then
          FoundWall_xy(4) = .True.
          ibc = M%IJKW(5,iw)         ! Boundary condition index
          ! There is a 'door', i.e., outflow-boundary (or open boundary)
          ! so no wall forces ==> exit this loop
          d_py = M%yw(iw)
          I_OBST = M%OBST_INDEX_C(IC)
          If (M%Solid(ic) .AND. .NOT.M%OBSTRUCTION(I_OBST)%HIDDEN) Then
             Write(MESSAGE,'(A,I4,2I6)') 'ERROR: Find_Walls ',nm, ii,jjn
             Call SHUTDOWN(MESSAGE)
          End If
          If (SURFACE(ibc)%VEL> 0.0_EB .Or. M%BOUNDARY_TYPE(iw)==OPEN_BOUNDARY) Then
             !d_py = d_py + is*(2.0_EB*d_cutoff)
             FoundWall_xy(4) = .False.
          End If
          If (Abs(d_py-r1_y) < r_circle) istat = 1
       End If
    End If
    d_xy(1) = d_mx
    d_xy(2) = d_px
    d_xy(3) = d_my
    d_xy(4) = d_py
    Return
  End Subroutine Find_walls

  SUBROUTINE GET_FIRE_CONDITIONS(NOM,I,J,K,fed_indx,soot_dens,gas_temp,rad_flux, YY_GET)
    Implicit None
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
    Call GET_MASS_FRACTION(YY_GET,SOOT_INDEX,Y_MF_INT)
    soot_dens = Y_MF_INT*MESHES(nom)%RHO(I,J,K)*1.E6_EB
    ! Calculate Purser's fractional effective dose (FED)
    fed_indx = FED(YY_GET,MESHES(nom)%RSUM(I,J,K))
    ! Gas temperature, ind=5, C
    gas_temp  = MESHES(nom)%TMP(I,J,K)
    ! Rad flux, ind=18, kW/m2 (no -sigma*Tamb^4 term)
    rad_flux = Max(MESHES(nom)%UII(I,J,K)/4.0_EB,SIGMA*TMPA4)

  END SUBROUTINE GET_FIRE_CONDITIONS

  Subroutine Change_Target_Door(nm, nm2, ie, j, j1, i_egrid, imode, xx, yy, I_Target, I_Color, I_Field, HR)
    Implicit None
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
    Integer, Intent(IN) :: nm, nm2, ie, j, j1, i_egrid, imode
    Real(EB), Intent(IN) :: xx, yy
    Integer, Intent(OUT)  :: I_Target, I_Color, I_Field
    Type (HUMAN_TYPE), Pointer :: HR
    !
    ! Local variables
    Real(EB) :: L2_min, max_fed, ave_K, L2_tmp, rn
    Real(EB) :: x1_old, y1_old, Speed, X11, Y11, x_o, y_o
    Integer :: i_old_ffield, i_tmp, i_new_ffield, IEL, color_index
    Integer :: i, i_o, izero, nm_tmp
    Character(26) :: name_old_ffield, name_new_ffield
    Logical :: PP_see_door
    Real(EB) :: T_tmp, T_tmp1, Width
    Integer :: N_queue, ii
    Type (EVACUATION_Type), Pointer :: HPT =>NULL()
    Type (EVAC_ENTR_Type),  Pointer :: PNX =>NULL()


    If (imode < 2) Then
       nm_tmp = nm
    Else
       ! This is for check_target_node
       nm_tmp = nm2
    End If

    i_old_ffield = HR%I_FFIELD
    If (i_old_ffield > 0) Then
       name_old_ffield = Trim(MESH_NAME(i_old_ffield))
    Else
       name_old_ffield = Trim(MESH_NAME(nm_tmp))
    End If
    If (HR%IEL > 0 ) Then
       ! Agent HR originates from an evac line
       HPT => EVACUATION(HR%IEL)
    Else
       ! Agent HR originates from an entr line (and is a lonely agent)
       PNX => EVAC_ENTRYS(Abs(HR%IEL))
    End If

    ! Only those doors are possible which are in the same main evac mesh.
    K_ave_Door(:)      = 0.0_EB
    FED_max_Door(:)    = 0.0_EB
    Is_Known_Door(:)   = .False.
    Is_Visible_Door(:) = .False.
    ! Group index=0: the agent is from an entry line (no evac line)
    If ( HR%GROUP_ID /= 0 ) Then
       If ( HR%GROUP_ID < 0 ) Then
          ! A lonely soul
          x1_old = xx
          y1_old = yy
          IEL    = HR%IEL
          Speed  = HR%Speed
          Do i = 1, n_doors
             If ( EVAC_DOORS(i)%IMESH == nm_tmp) Then
                Is_Visible_Door(i) = .True.
                If (imode == 0) Cycle   ! Initialization call
                Do i_tmp = 1, Human_Known_Doors(j1)%N_nodes
                   If (EVAC_DOORS(i)%INODE == Human_Known_Doors(j1)%I_nodes(i_tmp)) Is_Known_Door(i) = .True.
                End Do
             End If
          End Do
          Do i = 1, n_exits
             If ( EVAC_EXITS(i)%IMESH == nm_tmp .And. .Not. EVAC_EXITS(i)%COUNT_ONLY ) Then
                Is_Visible_Door(n_doors+i) = .True.
                If (imode == 0) Cycle   ! Initialization call
                Do i_tmp = 1, Human_Known_Doors(j1)%N_nodes
                   If (EVAC_EXITS(i)%INODE == Human_Known_Doors(j1)%I_nodes(i_tmp)) Is_Known_Door(n_doors+i) = .True.
                End Do
             End If
          End Do
       Else
          ! A member of a group
          x1_old = Group_List(j)%GROUP_X
          y1_old = Group_List(j)%GROUP_Y
          IEL    = Group_List(j)%IEL
          Speed  = Group_List(j)%Speed
          Do i = 1, n_doors
             If ( EVAC_DOORS(i)%IMESH == nm_tmp) Then
                Is_Visible_Door(i) = .True.
                If (imode == 0) Cycle   ! Initialization call
                Do i_tmp = 1, Group_Known_Doors(j)%N_nodes
                   If (EVAC_DOORS(i)%INODE == Group_Known_Doors(j)%I_nodes(i_tmp)) Is_Known_Door(i) = .True.
                End Do
             End If
          End Do
          Do i = 1, n_exits
             If ( EVAC_EXITS(i)%IMESH == nm_tmp .And. .Not. EVAC_EXITS(i)%COUNT_ONLY ) Then
                Is_Visible_Door(n_doors+i) = .True.
                If (imode == 0) Cycle   ! Initialization call
                Do i_tmp = 1, Group_Known_Doors(j)%N_nodes
                   If (EVAC_EXITS(i)%INODE == Group_Known_Doors(j)%I_nodes(i_tmp)) Is_Known_Door(n_doors+i) = .True.
                End Do
             End If
          End Do
       End If

       If (imode == 0) Then
          ! Initialization phase: Draw the known doors/exits.
          Do i = 1, HPT%N_VENT_FFIELDS 
             i_tmp = 1
             If (Trim(EVAC_Node_List(HPT%I_DOOR_NODES(i))%Node_Type) == 'Door') Then
                i_tmp = EVAC_Node_List(HPT%I_DOOR_NODES(i))%Node_Index 
             End If
             If (Trim(EVAC_Node_List(HPT%I_DOOR_NODES(i))%Node_Type) == 'Exit' ) Then
                i_tmp = n_doors + EVAC_Node_List(HPT%I_DOOR_NODES(i))%Node_Index 
             End If
             If (HPT%P_VENT_FFIELDS(i) < 1.0_EB) Then
                Call Random_number(RN)
                If ( RN < HPT%P_VENT_FFIELDS(i) ) Then
                   Is_Known_Door(i_tmp) = .True.
                Else
                   Is_Known_Door(i_tmp) = .False.
                End If
             Else
                Is_Known_Door(i_tmp) = .True.
             End If
          End Do

          ! Save the random known door information
          i_tmp = Count(Is_Known_Door)
          If (j > 0) Then    ! The agent is a member of a group
             Group_Known_Doors(j)%N_nodes = i_tmp
             If (Count(Is_Known_Door) > 0) Then
                Allocate(Group_Known_Doors(j)%I_nodes(i_tmp), STAT=IZERO)
                Call ChkMemErr('Change_Target_Door', 'Group_Known_Doors',IZERO) 
                i_tmp = 0
                Do i = 1, n_doors
                   If (Is_Known_Door(i)) Then
                      i_tmp = i_tmp + 1
                      Group_Known_Doors(j)%I_nodes(i_tmp) = EVAC_DOORS(i)%INODE  
                   End If
                End Do
                Do i = 1, n_exits
                   If (Is_Known_Door(i+n_doors)) Then
                      i_tmp = i_tmp + 1
                      Group_Known_Doors(j)%I_nodes(i_tmp) = EVAC_EXITS(i)%INODE  
                   End If
                End Do
             End If        ! there are known doors for this group
          Else    ! The agent is a lonely one
             Human_Known_Doors(j1)%N_nodes = i_tmp
             If (Count(Is_Known_Door) > 0) Then
                Allocate(Human_Known_Doors(j1)%I_nodes(i_tmp),STAT=IZERO)
                Call ChkMemErr('Change_Target_Door', 'Human_Known_Doors',IZERO) 
                i_tmp = 0
                Do i = 1, n_doors
                   If (Is_Known_Door(i)) Then
                      i_tmp = i_tmp + 1
                      Human_Known_Doors(j1)%I_nodes(i_tmp) = EVAC_DOORS(i)%INODE  
                   End If
                End Do
                Do i = 1, n_exits
                   If (Is_Known_Door(i+n_doors)) Then
                      i_tmp = i_tmp + 1
                      Human_Known_Doors(j1)%I_nodes(i_tmp) = EVAC_EXITS(i)%INODE  
                   End If
                End Do
             End If   ! Is there any known doors for this group?
          End If    ! Is the agent a member of a group or not?
          
          ! Check that the door is in correct main evac mesh.
          Do i = 1, n_doors
             If ( EVAC_DOORS(i)%IMESH /= HPT%imesh ) Then
                Is_Known_Door(i) = .False.
             End If
          End Do
          Do i = 1, n_exits
             If ( EVAC_EXITS(i)%IMESH /= HPT%imesh .Or. EVAC_EXITS(i)%COUNT_ONLY ) Then
                Is_Known_Door(n_doors+i) = .False.
             End If
          End Do
       End If  ! imode=0, Initialization call

    Else
       ! The agent is from an entry line. P_door= 0.0 or 1.0, known doors
       x1_old = xx
       y1_old = yy
       IEL    = Abs(HR%IEL)
       Speed  = HR%Speed
       Do i = 1, n_doors
          If ( EVAC_DOORS(i)%IMESH == nm_tmp) Then
             Is_Visible_Door(i) = .True.
          End If
       End Do
       Do i = 1, n_exits
          If ( EVAC_EXITS(i)%IMESH == nm_tmp .And. .Not. EVAC_EXITS(i)%COUNT_ONLY ) Then
             Is_Visible_Door(n_doors+i) = .True.
          End If
       End Do
       Do i = 1, PNX%N_VENT_FFIELDS 
          ! Check that the door/exit is in the correct main evac grid.
          i_tmp = 1
          If (Trim(EVAC_Node_List(PNX%I_DOOR_NODES(i))%Node_Type) == 'Door') Then
             i_tmp = EVAC_Node_List(PNX%I_DOOR_NODES(i))%Node_Index 
          End If
          If (Trim(EVAC_Node_List(PNX%I_DOOR_NODES(i))%Node_Type) == 'Exit' ) Then
             i_tmp = n_doors + EVAC_Node_List(PNX%I_DOOR_NODES(i))%Node_Index 
          End If
          If ( Is_Visible_Door(i_tmp) ) Then
             ! Door/exit is on this main evac mesh.
             If (PNX%P_VENT_FFIELDS(i) < 0.5_EB) Then
                Is_Known_Door(i_tmp) = .False.
             Else
                Is_Known_Door(i_tmp) = .True.
             End If
          End If
       End Do
    End If   ! Is the agent from an entr or from an evac line

    ! Find the visible doors.
    Do i = 1, n_doors + n_exits
       If ( Is_Visible_Door(i) ) Then
          If (EVAC_Node_List(n_egrids+n_entrys+i)%Node_Type == 'Door' ) Then
             X11 = EVAC_DOORS(i)%X 
             Y11 = EVAC_DOORS(i)%Y
          Else
             X11 = EVAC_EXITS(i-n_doors)%X 
             Y11 = EVAC_EXITS(i-n_doors)%Y
          End If
          ! Groups: the first member (x1_old,y1_old) of the group is used.
          PP_see_door = See_door(nm_tmp, i, HR%I_Target, x1_old, y1_old, x11, y11, ave_K, max_fed)
          FED_max_Door(i) = max_fed
          K_ave_Door(i) = ave_K 

          ! Note: a DOOR is not counted as visible door, if it does not have an
          ! EXIT_SIGN, unless it is already been a target door for this agent/group.
          If (PP_see_door) Then
             If (EVAC_Node_List(n_egrids+n_entrys+i)%Node_Type == 'Door') Then
                If (.Not. EVAC_DOORS(i)%EXIT_SIGN .And. .Not. HR%I_Target == i) Then
                   Is_Visible_Door(i) = .False.
                End If
             End If
          Else
             Is_Visible_Door(i) = .False.
          End If
       End If ! correct main evac mesh
    End Do ! all doors and exits


    ! Note: I_Target < 0: not visible, >0: visible
    If (Any(Is_Visible_Door) .And. imode == 1) Then
       Do i = 1, n_doors + n_exits
          If (Abs(HR%I_Target) == i .And. Is_Visible_Door(i)) Is_Known_Door(i) = .True.
       End Do
    End If

    Do i = 1, n_doors
       ! If ( EVAC_DOORS(i)%TIME_OPEN > T .Or. EVAC_DOORS(i)%TIME_CLOSE < T) Then
       If ( Abs(EVAC_DOORS(i)%IMODE)==2) Then
          Is_Visible_Door(i) = .False.
          Is_Known_Door(i) = .False.
       End If
    End Do
    Do i = 1, n_exits
       ! If ( (EVAC_EXITS(i)%TIME_OPEN > T .Or. EVAC_EXITS(i)%TIME_CLOSE < T) .And. .Not. EVAC_EXITS(i)%COUNT_ONLY ) Then
       If ( (Abs(EVAC_EXITS(i)%IMODE)==2) .And. .Not. EVAC_EXITS(i)%COUNT_ONLY ) Then
          Is_Visible_Door(n_doors+i) = .False.
          Is_Known_Door(n_doors+i) = .False.
       End If
    End Do

    If (Any(Is_Known_Door) .Or. Any(Is_Visible_Door)) Then
       i_tmp   = 0
       L2_min = Huge(L2_min)
       Do i = 1, n_doors + n_exits
          If ( Is_Known_Door(i) .And. Is_Visible_Door(i) ) Then
             x_o = 0.0_EB
             y_o = 0.0_EB
             N_queue = 0
             If (Trim(EVAC_Node_List(n_egrids+n_entrys+i)%Node_Type) == 'Door' ) Then
                x_o = EVAC_DOORS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%X
                y_o = EVAC_DOORS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%Y
                i_o = EVAC_DOORS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%I_VENT_FFIELD
                T_tmp1 = 50.0_EB*Sqrt((x1_old-x_o)**2 + (y1_old-y_o)**2)/ &
                     EVAC_DOORS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%R_NTARGET + 1.0_EB
                ii = Min(50,Max(1,Int(T_tmp1)-1))
                Width = EVAC_DOORS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%Width
                N_queue = EVAC_DOORS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%NTARGET(ii)
             Else      ! 'Exit'
                x_o = EVAC_EXITS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%X
                y_o = EVAC_EXITS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%Y
                i_o = EVAC_EXITS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%I_VENT_FFIELD
                T_tmp1 = 50.0_EB*Sqrt((x1_old-x_o)**2 + (y1_old-y_o)**2)/ &
                     EVAC_EXITS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%R_NTARGET + 1.0_EB
                ii = Min(50,Max(1,Int(T_tmp1)-1))
                Width = EVAC_EXITS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%Width
                N_queue = EVAC_EXITS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%NTARGET(ii)
             End If
             If (FED_DOOR_CRIT > 0.0_EB) Then
                L2_tmp = FED_max_Door(i) * Sqrt((x1_old-x_o)**2 + (y1_old-y_o)**2)/Speed
             Else
                L2_tmp = K_ave_Door(i)
             End If
             If (i_o == i_old_ffield) L2_tmp = 0.1_EB*L2_tmp
             If (FAC_DOOR_QUEUE > 0.001_EB) Then
                T_tmp  = Sqrt((x_o-x1_old)**2 + (y_o-y1_old)**2)
                T_tmp1 = Min(1.5_EB*Pi*T_tmp**2/(FAC_DOOR_QUEUE*Width), Real(N_queue,EB)/(FAC_DOOR_QUEUE*Width))
                T_tmp = (T_tmp/Speed) +  T_tmp1
                If (i_o == i_old_ffield) T_tmp = T_tmp*FAC_DOOR_WAIT
                If ( T_tmp < L2_min .And. L2_tmp < Abs(FED_DOOR_CRIT) ) Then
                   L2_min = Max(0.0_EB,T_tmp)
                   i_tmp = i
                End If
             Else
                T_tmp  = Sqrt((x_o-x1_old)**2 + (y_o-y1_old)**2)
                If (i_o == i_old_ffield) T_tmp = T_tmp*FAC_DOOR_WAIT
                If (T_tmp < L2_min .And. L2_tmp < Abs(FED_DOOR_CRIT)) Then
                   L2_min = Max(0.0_EB,T_tmp)
                   i_tmp = i
                End If
             End If
          End If
       End Do
       If (i_tmp > 0 ) Then
          ! Known and visible door, no smoke
          If (EVAC_Node_List(n_egrids+n_entrys+i_tmp)%Node_Type == 'Door' ) Then
             name_new_ffield = Trim(EVAC_DOORS(i_tmp)%VENT_FFIELD)
             i_new_ffield = EVAC_DOORS(i_tmp)%I_VENT_FFIELD
          Else        ! 'Exit'
             name_new_ffield = Trim(EVAC_EXITS(i_tmp-n_doors)%VENT_FFIELD)
             i_new_ffield = EVAC_EXITS(i_tmp-n_doors)%I_VENT_FFIELD
          End If
          color_index = 1
       Else
          ! No visible known doors available, try non-visible known doors
          i_tmp   = 0
          L2_min = Huge(L2_min)
          Do i = 1, n_doors + n_exits
             If ( Is_Known_Door(i) .And. .Not. Is_Visible_Door(i) ) Then
                x_o = 0.0_EB
                y_o = 0.0_EB
                If (EVAC_Node_List(n_egrids+n_entrys+i)%Node_Type == 'Door' ) Then
                   x_o = EVAC_DOORS( EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%X
                   y_o = EVAC_DOORS( EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%Y
                   i_o = EVAC_DOORS( EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%I_VENT_FFIELD
                Else    ! 'Exit'
                   x_o = EVAC_EXITS( EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%X
                   y_o = EVAC_EXITS( EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%Y
                   i_o = EVAC_EXITS( EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%I_VENT_FFIELD
                End If
                If (FED_DOOR_CRIT > 0.0_EB) Then
                   L2_tmp = FED_max_Door(i)*Sqrt((x1_old-x_o)**2 + (y1_old-y_o)**2)/Speed
                Else
                   l2_tmp = K_ave_Door(i)
                End If
                If (i_o == i_old_ffield) L2_tmp = 0.1_EB*L2_tmp
                T_tmp  = Sqrt((x_o-x1_old)**2 + (y_o-y1_old)**2)
                If (i_o == i_old_ffield) T_tmp = T_tmp*FAC_DOOR_WAIT
                If (T_tmp < L2_min .And. L2_tmp < Abs(FED_DOOR_CRIT)) Then
                   L2_min = Max(0.0_EB,T_tmp)
                   i_tmp = i
                End If
             End If
          End Do
          If (i_tmp > 0 ) Then
             !  Non-visible known door, no smoke is found
             If (EVAC_Node_List( n_egrids+n_entrys+i_tmp)%Node_Type == 'Door' ) Then
                name_new_ffield = Trim(EVAC_DOORS(i_tmp)%VENT_FFIELD)
                i_new_ffield = EVAC_DOORS(i_tmp)%I_VENT_FFIELD
             Else      ! 'Exit'
                name_new_ffield = Trim(EVAC_EXITS(i_tmp-n_doors)%VENT_FFIELD)
                i_new_ffield = EVAC_EXITS(i_tmp-n_doors)%I_VENT_FFIELD
             End If
             color_index = 2
          Else
             ! Known doors with no smoke have not been found.
             ! Try visible, not known door with no smoke.
             i_tmp   = 0
             L2_min = Huge(L2_min)
             Do i = 1, n_doors + n_exits
                If (Is_Visible_Door(i)) Then
                   x_o = 0.0_EB
                   y_o = 0.0_EB
                   If (EVAC_Node_List(n_egrids+n_entrys+i)%Node_Type == 'Door' ) Then
                      x_o = EVAC_DOORS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%X
                      y_o = EVAC_DOORS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%Y
                      i_o = EVAC_DOORS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%I_VENT_FFIELD
                   Else  ! 'Exit'
                      x_o = EVAC_EXITS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%X
                      y_o = EVAC_EXITS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%Y
                      i_o = EVAC_EXITS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%I_VENT_FFIELD
                   End If
                   If (FED_DOOR_CRIT > 0.0_EB) Then
                      L2_tmp = FED_max_Door(i)*Sqrt((x1_old-x_o)**2 + (y1_old-y_o)**2)/Speed
                   Else
                      l2_tmp = K_ave_Door(i)
                   End If
                   If (i_o == i_old_ffield) L2_tmp = 0.1_EB*L2_tmp
                   T_tmp  = Sqrt((x_o-x1_old)**2 + (y_o-y1_old)**2)
                   If (i_o == i_old_ffield) T_tmp = T_tmp*FAC_DOOR_WAIT
                   If (T_tmp < L2_min .And. L2_tmp < Abs(FED_DOOR_CRIT)) Then
                      L2_min = Max(0.0_EB,T_tmp)
                      i_tmp = i
                   End If
                End If
             End Do
             If (i_tmp > 0 ) Then
                ! No smoke, visible door (not known) is found
                If (EVAC_Node_List( n_egrids+n_entrys+i_tmp)%Node_Type == 'Door' ) Then
                   name_new_ffield = Trim(EVAC_DOORS(i_tmp)%VENT_FFIELD)
                   i_new_ffield = EVAC_DOORS(i_tmp)%I_VENT_FFIELD
                Else    ! 'Exit'
                   name_new_ffield = Trim( EVAC_EXITS(i_tmp-n_doors)%VENT_FFIELD)
                   i_new_ffield = EVAC_EXITS(i_tmp-n_doors)%I_VENT_FFIELD
                End If
                color_index = 3
             Else
                ! Now we have some smoke and some visible or known doors
                i_tmp   = 0
                L2_min = Huge(L2_min)
                Do i = 1, n_doors + n_exits
                   If ( Is_Known_Door(i) .Or. Is_Visible_Door(i) ) Then
                      x_o = 0.0_EB
                      y_o = 0.0_EB
                      If (EVAC_Node_List(n_egrids+n_entrys+i)%Node_Type == 'Door' ) Then
                         x_o = EVAC_DOORS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%X
                         y_o = EVAC_DOORS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%Y
                         i_o = EVAC_DOORS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%I_VENT_FFIELD
                      Else ! 'Exit'
                         x_o = EVAC_EXITS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%X
                         y_o = EVAC_EXITS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%Y
                         i_o = EVAC_EXITS(EVAC_Node_List(i+n_egrids+n_entrys)%Node_Index)%I_VENT_FFIELD
                      End If
                      If (FED_DOOR_CRIT > 0.0_EB) Then
                         If (j > 0) Then
                            L2_tmp = Group_List(j)%IntDose + FED_max_Door(i) * Sqrt((x1_old-x_o)**2+(y1_old-y_o)**2)/Speed
                         Else
                            L2_tmp = HR%IntDose + FED_max_Door(i) * Sqrt((x1_old-x_o)**2+(y1_old-y_o)**2)/Speed
                         End If
                      Else
                         l2_tmp = (Sqrt((x1_old-x_o)**2+(y1_old-y_o)**2)*0.5_EB)/(3.0_EB/K_ave_Door(i))
                      End If
                      If (i_o == i_old_ffield) L2_tmp = 0.9_EB*L2_tmp
                      If (L2_tmp < L2_min) Then
                         L2_min = L2_tmp
                         i_tmp = i
                      End If
                   End If
                End Do

                If (i_tmp > 0 .And. L2_min < 1.0_EB) Then
                   ! Not too much smoke, (known and/or visible doors)
                   If (EVAC_Node_List(n_egrids+n_entrys+i_tmp)%Node_Type == 'Door' ) Then
                      name_new_ffield = Trim(EVAC_DOORS(i_tmp)%VENT_FFIELD)
                      i_new_ffield = EVAC_DOORS(i_tmp)%I_VENT_FFIELD
                   Else  ! 'Exit'
                      name_new_ffield = Trim(EVAC_EXITS(i_tmp-n_doors)%VENT_FFIELD)
                      i_new_ffield = EVAC_EXITS(i_tmp-n_doors)%I_VENT_FFIELD
                   End If
                   If (Is_Known_Door(i_tmp) .And. Is_Visible_Door(i_tmp)) color_index = 4
                   If (Is_Known_Door(i_tmp) .And. .Not. Is_Visible_Door(i_tmp)) color_index = 5
                   If (.Not. Is_Known_Door(i_tmp) .And. Is_Visible_Door(i_tmp)) color_index = 6
                Else    ! no match 
                   ! No door found (or too much smoke), use the main evac grid ffield
                   ! Note: This should be changed in later versions of the program.
                   i_tmp = 0
                   name_new_ffield = Trim(MESH_NAME(nm_tmp))
                   i_new_ffield    = nm_tmp
                   If (imode == 0) Then  ! Initialization call, use evac line info
                      name_new_ffield = Trim(HPT%GRID_NAME)
                      i_new_ffield    = HPT%I_VENT_FFIELDS(0)
                   End If
                   color_index = EVAC_AVATAR_NCOLOR ! default, cyan
                End If  ! case 4
             End If    ! case 3
          End If      ! case 2
       End If        ! case 1
       If (Color_Method == 4 ) Then
          color_index = EVAC_AVATAR_NCOLOR ! default, cyan
          If (i_tmp > 0 .And. i_tmp <= n_doors ) color_index = EVAC_DOORS(i_tmp)%Avatar_Color_Index
          If (i_tmp > n_doors .And. i_tmp <= n_doors + n_exits) &
               color_index = EVAC_EXITS(i_tmp-n_doors)%Avatar_Color_Index
       End If
       If (imode == 2) Then   ! check_target_node calls
          I_Target = i_tmp
          If (i_tmp > 0 .And. .Not. Is_Visible_Door(Max(1,i_tmp)) ) Then
             ! I_Target >0: visible, <0: not visible
             I_Target = -i_tmp
          End If
          I_Color  = color_index
          I_Field  = i_new_ffield
          Return
       End If

    Else ! No known/visible door
       i_tmp = 0 ! no door found
       color_index = EVAC_AVATAR_NCOLOR ! default, cyan
       If (imode == 2) Then   ! check_target_node calls
          I_Target = 0
          I_Color  = color_index
          I_Field  = nm_tmp
          Return
       End If
       If (HR%IEL > 0 ) Then  
          ! The agent is from some evac line
          If (HPT%IMESH == nm_tmp) Then
             i_new_ffield    = HPT%I_VENT_FFIELDS(0)
             name_new_ffield = Trim(Mesh_Name(i_new_ffield))
          Else
             name_new_ffield = Trim(MESH_NAME(nm_tmp))
             i_new_ffield    = nm_tmp
          End If
          If (Color_Method == 4) color_index = EVAC_AVATAR_NCOLOR
       Else
          ! The agent is from some entr line
          If (PNX%IMESH == nm_tmp) Then
             i_new_ffield    = PNX%I_VENT_FFIELDS(0)
             name_new_ffield = Trim(Mesh_Name(i_new_ffield))
          Else
             name_new_ffield = Trim(MESH_NAME(nm_tmp))
             i_new_ffield    = nm_tmp
          End If
          If (Color_Method == 4) color_index = EVAC_AVATAR_NCOLOR
       End If
    End If ! Any known or visible door
    HR%I_Target = i_tmp
    If (imode < 2) I_Target = i_tmp
    If (i_tmp > 0 .And. .Not. Is_Visible_Door(Max(1,i_tmp))) Then
       ! I_Target > 0: visible, < 0: not visible
       HR%I_Target = -i_tmp
       If (imode < 2) I_Target = -i_tmp
    End If
    If (j > 0) Group_Known_Doors(j)%I_Target = HR%I_Target
    If ( (i_new_ffield /= i_old_ffield) .Or. (imode == 0) ) Then
       ! Change the field of this group/agent.
       If ( j == 0 ) Then
          ! Group index=0: 'lost souls' or lonely agents
          HR%I_FFIELD    = i_new_ffield
          HR%FFIELD_NAME = Trim(name_new_ffield)
          If (COLOR_METHOD == 5) HR%COLOR_INDEX = color_index
          If (COLOR_METHOD == 4) HR%COLOR_INDEX = color_index
          If (FAC_DOOR_QUEUE > 0.001_EB) Return
          If (imode > 0) Then
             Write (LU_EVACOUT,fmt='(a,i5,a,a,a,a)') ' EVAC: Human ',ie,', new ffield: ', &
                  Trim(name_new_ffield), ', old ffield: ',Trim(name_old_ffield)
          Else
             Write (LU_EVACOUT,fmt='(a,i5,a,a)') ' EVAC: Human ',ie,', flow field: ', Trim(name_new_ffield)
          End If
       Else
          Group_Known_Doors(j)%I_Target = HR%I_Target
          Group_List(j)%GROUP_I_FFIELDS(i_egrid) = i_new_ffield
          HR%I_FFIELD = i_new_ffield
          HR%FFIELD_NAME = Trim(name_new_ffield)
          If (COLOR_METHOD == 5) HR%COLOR_INDEX = color_index
          If (COLOR_METHOD == 4) HR%COLOR_INDEX = color_index
          Color_Tmp(j) = color_index
          If (FAC_DOOR_QUEUE > 0.001_EB) Return
          If (imode > 0) Then
             Write (LU_EVACOUT,fmt='(a,i5,a,a,a,a)') ' EVAC: Group ',j,', new ffield: ', &
                  Trim(name_new_ffield), ', old ffield: ', Trim(name_old_ffield)
          Else
             Write (LU_EVACOUT,fmt='(a,i5,a,a)') ' EVAC: Group ',ie,', flow field: ', Trim(name_new_ffield)
          End If
       End If
    Else ! The new door is the same as the old target door.
       If (COLOR_METHOD == 5) HR%COLOR_INDEX = color_index
       If (COLOR_METHOD == 5 .And. j > 0) Color_Tmp(j) = HR%COLOR_INDEX
       If (COLOR_METHOD == 4) HR%COLOR_INDEX = color_index
       If (COLOR_METHOD == 4 .And. j > 0) Color_Tmp(j) = HR%COLOR_INDEX
    End If

  End Subroutine Change_Target_Door
  !
  Subroutine GET_REV_evac(MODULE_REV,MODULE_DATE)
    !
    ! Passed variables
    Integer,Intent(INOUT) :: MODULE_REV
    Character(255),Intent(INOUT) :: MODULE_DATE
    !
    ! Local variables
    !

    Write(MODULE_DATE,'(A)') evacrev(Index(evacrev,':')+1:Len_trim(evacrev)-2)
    Read (MODULE_DATE,'(I5)') MODULE_REV
    Write(MODULE_DATE,'(A)') evacdate

  End Subroutine GET_REV_evac
  
End Module EVAC
