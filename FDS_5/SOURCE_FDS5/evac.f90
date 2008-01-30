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
! Author: Timo Korhonen, VTT Technical Research Centre of Finland, 2007
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
  Use PHYSICAL_FUNCTIONS, ONLY : GET_MASS_FRACTION2
  !
  Implicit None
  CHARACTER(255), PARAMETER :: evacid='$Id$'
  CHARACTER(255), PARAMETER :: evacrev='$Revision$'
  CHARACTER(255), PARAMETER :: evacdate='$Date$'
  !
  Private
  ! Public subprograms (called from the main program)
  Public EVACUATE_HUMANS, INITIALIZE_EVACUATION, INIT_EVAC_GROUPS
  Public READ_EVAC, DUMP_EVAC, DUMP_EVAC_CSV, PREPARE_TO_EVACUATE
  Public EVAC_MESH_EXCHANGE, INITIALIZE_EVAC_DUMPS, GET_REV_evac
  ! Public variables (needed in the main program):
  !
  Character(255):: EVAC_VERSION = '1.10'
  Character(255) :: EVAC_COMPILE_DATE
  INTEGER :: EVAC_MODULE_REV
  !
  ! This is a group of persons, who are initialized together,
  ! i.e., they have same mass, speed, etc distributions and
  ! they are all put in the given rectangle.
  ! (&EVAC lines)
  Type EVACUATION_Type
     Real(EB) :: X1=0._EB,X2=0._EB,Y1=0._EB,Y2=0._EB,Z1=0._EB,Z2=0._EB,T_START=0._EB, Angle=0._EB
     Character(60) :: CLASS_NAME='null', ID_NAME='null'
     Character(30) :: GRID_NAME='null'
     Logical :: EVACFILE=.FALSE., After_Tpre=.FALSE., No_Persons=.FALSE.
     Integer :: N_INITIAL=0,COLOR_INDEX=0,SAMPLING=0, IPC=0, IMESH=0
     Integer :: GN_MIN=0, GN_MAX=0
     Integer :: N_VENT_FFIELDS=0
     Integer, Pointer, Dimension(:) :: I_DOOR_NODES
     Integer, Pointer, Dimension(:) :: I_VENT_FFIELDS
     Real(EB), Pointer, Dimension(:) :: P_VENT_FFIELDS
  End Type EVACUATION_Type
  !
  ! An evacuatio hole, i.e., a rectangle where humans should
  ! not be put.  This makes the &EVAC lines easier to define.
  ! (&EVHO lines)
  Type EVAC_HOLE_Type
     Real(EB) :: X1=0._EB,X2=0._EB,Y1=0._EB,Y2=0._EB,Z1=0._EB,Z2=0._EB
     Character(60) :: ID_NAME='null', PERS_ID='null', EVAC_ID='null'
     Character(30) :: GRID_NAME='null'
     Integer :: IMESH=0
  End Type EVAC_HOLE_Type
  !
  ! A spectator stand. IOR: which x,y line is the bottom line of the stand.
  ! ior=+1 x=x2, ior=-1 x=x1, ior=+2 y=y2, ior=-2 y=y1
  ! H is the height of the stand, S is the length along the incline.
  ! (&EVSS lines)
  Type EVAC_SSTAND_TYPE
     Real(EB) :: X1=0._EB,X2=0._EB,Y1=0._EB,Y2=0._EB,Z1=0._EB,Z2=0._EB, H=0._EB, H0=0._EB, S=0._EB
     Real(EB) :: Esc_SpeedUp=0._EB, Esc_SpeedDn=0._EB
     Real(EB) :: fac_v0_up=1._EB, fac_v0_down=1._EB, fac_v0_hori=1._EB
     Real(EB) :: cos_x=1._EB, cos_y=1._EB, sin_x=0._EB, sin_y=0._EB
     Character(60) :: ID_NAME='null'
     Character(26) :: GRID_NAME='null'
     Integer :: IMESH=0, IOR=0
  End Type EVAC_SSTAND_TYPE
  !
  ! Humans belong to some small group (1 to about 5 persons).  This type
  ! collects the common properties of the group.
  Type GROUP_TYPE
     Real(EB) :: GROUP_X=0._EB, GROUP_Y=0._EB, MAX_DIST_CENTER=0._EB, LIMIT_COMP=0._EB
     Real(EB) :: GROUP_EFF=0._EB, RADIUS_COMPLETE_0=0._EB, RADIUS_COMPLETE_1=0._EB
     Real(EB) :: Speed=0._EB, IntDose=0._EB, Tpre=0._EB, Tdoor=0._EB, Tdet=0._EB
     Integer :: GROUP_SIZE=0, GROUP_ID=0, COMPLETE=0, IEL=0
     Integer, Pointer, Dimension(:) :: GROUP_I_FFIELDS
  End Type GROUP_TYPE
  
  Type KNOWN_DOOR_TYPE
     Integer :: N_nodes=0, I_Target=0
     Integer, Pointer, Dimension(:) :: I_nodes
  End Type KNOWN_DOOR_TYPE
  !
  ! This defines a class of persons, e.g. soccer fan.
  ! (&PERS lines)
  Type EVAC_PERS_Type
     Real(EB) :: D_mean=0._EB, D_para=0._EB, D_para2=0._EB, D_low=0._EB, D_high=0._EB
     Real(EB) :: V_mean=0._EB, V_para=0._EB, V_para2=0._EB, V_low=0._EB, V_high=0._EB
     Real(EB) :: Tau_mean=0._EB, Tau_para=0._EB, Tau_para2=0._EB, Tau_low=0._EB, Tau_high=0._EB
     Real(EB) :: Tpre_mean=0._EB, Tpre_para=0._EB, Tpre_para2=0._EB, Tpre_low=0._EB, Tpre_high=0._EB
     Real(EB) :: Tdet_mean=0._EB, Tdet_para=0._EB, Tdet_para2=0._EB, Tdet_low=0._EB, Tdet_high=0._EB
     Real(EB) :: A=0._EB,B=0._EB,Lambda=0._EB,C_Young=0._EB,Gamma=0._EB,Kappa=0._EB
     Real(EB) :: r_torso=0._EB,r_shoulder=0._EB,d_shoulder=0._EB,m_iner=0._EB, Tau_iner=0._EB
     Character(60) :: ID_NAME='null'
     Integer :: I_DIA_DIST=0, I_VEL_DIST=0, I_PRE_DIST=0, I_DET_DIST=0, I_TAU_DIST=0
  End Type EVAC_PERS_Type
  !
  ! Exit door type: this just count the number of persons
  ! T_first: first person's exit time (saved for output)
  ! CHECK_FLOW: If true then the flow can not exceed Flow_max
  ! (&EXIT lines)
  Type EVAC_EXIT_Type
     Real(EB) :: T_first=0._EB, T_last=0._EB, Flow_max=0._EB, Width=0._EB
     Real(EB) :: X1=0._EB, X2=0._EB, Y1=0._EB, Y2=0._EB, Z1=0._EB, Z2=0._EB, &
          X=0._EB, Y=0._EB, Z=0._EB, Xsmoke=0._EB, Ysmoke=0._EB, Zsmoke=0._EB, &
          TIME_OPEN=0._EB, TIME_CLOSE=0._EB
     Integer :: IOR=0, ICOUNT=0, IMESH=0, INODE=0, NTARGET=0
     Real(EB) :: FED_CO_CO2_O2=0._EB, SOOT_DENS=0._EB, TMP_G=0._EB, RADINT=0._EB
     Integer :: II=0, JJ=0, KK=0, FED_MESH=0, COLOR_INDEX=0
     Logical :: CHECK_FLOW=.FALSE., COUNT_ONLY=.FALSE.
     Character(60) :: ID_NAME='null'
     Character(60) :: TO_NODE='null'
     Character(30) :: GRID_NAME='null'
     Character(26) :: VENT_FFIELD='null'
     Integer :: I_VENT_FFIELD=0
  End Type EVAC_EXIT_Type
  !
  ! Like exit, but door will always put the persons to some
  ! other node. (Thus no count_only option.)
  ! (&DOOR lines)
  Type EVAC_DOOR_Type
     Real(EB) :: T_first=0._EB, T_last=0._EB, Flow_max=0._EB, Width=0._EB
     Real(EB) :: X1=0._EB, X2=0._EB, Y1=0._EB, Y2=0._EB, Z1=0._EB, Z2=0._EB, &
          X=0._EB, Y=0._EB, Z=0._EB, Xsmoke=0._EB, Ysmoke=0._EB, Zsmoke=0._EB, &
          TIME_OPEN=0._EB, TIME_CLOSE=0._EB
     Integer :: IOR=0, ICOUNT=0, INODE=0, INODE2=0, IMESH=0, IMESH2=0, NTARGET=0
     Real(EB) :: FED_CO_CO2_O2=0._EB, SOOT_DENS=0._EB, TMP_G=0._EB, RADINT=0._EB
     Integer :: II=0, JJ=0, KK=0, FED_MESH=0, COLOR_INDEX=0
     Logical :: CHECK_FLOW=.FALSE., EXIT_SIGN=.FALSE., KEEP_XY=.FALSE.
     Character(60) :: ID_NAME='null'
     Character(60) :: TO_NODE='null'
     Character(30) :: GRID_NAME='null'
     Character(26) :: VENT_FFIELD='null'
     Integer :: I_VENT_FFIELD=0
  End Type EVAC_DOOR_Type
  !
  ! Like door, but corr will model stairs (or corridors). 
  ! The parameters, like velocity as function of density etc.
  ! define if it is corridor or stairway
  ! (&CORR lines)
  Type EVAC_CORR_Type
     Real(EB) :: T_first=0._EB, T_last=0._EB, Flow_max=0._EB, Width1=0._EB, Width2=0._EB
     Real(EB) :: X1=0._EB,X2=0._EB,Y1=0._EB,Y2=0._EB,Z1=0._EB,Z2=0._EB, Width=0._EB
     Real(EB) :: Eff_Width=0._EB, Eff_Length=0._EB, Eff_Area=0._EB, Fac_Speed=0._EB
     ! Note: Corridor may have 2 different points, where smoke etc. data
     ! is saved.
     Real(EB), Dimension(2) :: FED_CO_CO2_O2=0._EB, SOOT_DENS=0._EB, TMP_G=0._EB, RADINT=0._EB
     Integer :: FED_MESH=0, FED_MESH2=0
     Integer, Dimension(2) :: II=0, JJ=0, KK=0
     Integer :: IOR=0, ICOUNT=0, INODE=0, INODE2=0, IMESH=0, IMESH2=0
     Integer :: MAX_HUMANS_INSIDE=0, n_inside=0
     Logical :: CHECK_FLOW=.FALSE.
     Character(60) :: ID_NAME='null'
     Character(60) :: TO_NODE='null'
     Character(30) :: GRID_NAME='null'
     Type (CORR_LL_Type), Pointer :: First
  End Type EVAC_CORR_Type
  !
  ! This produces more humans on the floor specified by the
  ! coordinates. the person type ('soccer_fan' etc) are also
  ! defined here for these persons.
  ! (&ENTR lines)
  Type EVAC_ENTR_Type
     Real(EB) :: T_first=0._EB, T_last=0._EB, Flow=0._EB, Width=0._EB, T_Start=0._EB, T_Stop=0._EB
     Real(EB) :: X1=0._EB,X2=0._EB,Y1=0._EB,Y2=0._EB,Z1=0._EB,Z2=0._EB
     Integer :: IOR=0, ICOUNT=0, COLOR_INDEX=0, IPC=0, IMESH=0, INODE=0, &
          TO_INODE=0, N_Initial=0
     Character(60) :: CLASS_NAME='null', ID_NAME='null'
     Character(60) :: TO_NODE='null'
     Character(30) :: GRID_NAME='null'
     Logical :: After_Tpre=.FALSE., No_Persons=.FALSE.
     Integer :: N_VENT_FFIELDS=0
     Integer, Pointer, Dimension(:) :: I_DOOR_NODES
     Integer, Pointer, Dimension(:) :: I_VENT_FFIELDS
     Real(EB), Pointer, Dimension(:) :: P_VENT_FFIELDS
  End Type EVAC_ENTR_Type
  !
  ! coordinates. the person type ('soccer_fan' etc) are also
  ! defined here for these persons.
  Type EVAC_NODE_Type
     Integer :: Node_Index=0, Mesh_Index=0
     Character(60) :: ID_NAME='null', Node_Type='null'
     Character(30) :: GRID_NAME='null'
  End Type EVAC_NODE_Type
  !
  ! Linked list, needed for the corridors
  Type CORR_LL_Type
     Type (HUMAN_Type) :: HUMAN
     Real(EB) :: T_in=0._EB, T_out=0._EB
     Logical :: From1_To2=.FALSE.
     Integer :: Index=0
     Type (CORR_LL_Type), Pointer :: Next
  End Type CORR_LL_Type
  !
  ! Pointers to the allocatable arrays so one can use these as
  ! shorthands to the array elements.
  Type (KNOWN_DOOR_Type), Pointer :: KDT
  Type (GROUP_Type),      Pointer :: GR
  Type (HUMAN_Type),      Pointer :: HR, HRE
  Type (EVACUATION_Type), Pointer :: HPT, HPE
  Type (EVAC_PERS_Type),  Pointer :: PCP
  Type (EVAC_EXIT_Type),  Pointer :: PEX
  Type (EVAC_DOOR_Type),  Pointer :: PDX, PDX2
  Type (EVAC_ENTR_Type),  Pointer :: PNX, PNX2
  Type (EVAC_CORR_Type),  Pointer :: PCX, PCX2
  Type (EVAC_NODE_Type),  Pointer :: NODE
  Type (EVAC_HOLE_Type),  Pointer :: EHX
  Type (EVAC_SSTAND_Type),Pointer :: ESS

  !
  ! Next holds door information for groups
  Type (KNOWN_DOOR_TYPE), Dimension(:), Allocatable, Target :: &
       Group_Known_Doors
  ! Next holds door information for lonely humans (group_id=0)
  Type (KNOWN_DOOR_TYPE), Dimension(:), Allocatable, Target :: &
       Human_Known_Doors
  Integer :: ilh, ilh_dim
  
  ! Holds the list of the different human groups, i33 is a running index 
  ! for the groups, i33_dim is last index, i.e., the dimension of the array.
  Type (GROUP_TYPE), Dimension(:), Allocatable, Target :: &
       Group_List
  Integer :: i33, i33_dim

  ! Holds the information of the nodes
  Type (EVAC_NODE_Type), Dimension(:), Allocatable, Target ::  &
       Evac_Node_List

  ! Holds the information of the EVAC-lines.
  Type (EVACUATION_Type), Dimension(:), Allocatable, Target ::  &
       EVACUATION

  ! Holds the information of the EVHO-lines.
  Type (EVAC_HOLE_Type), Dimension(:), Allocatable, Target ::  &
       EVAC_HOLES

  ! Holds the information of the EVSS-lines.
  Type (EVAC_SSTAND_TYPE), Dimension(:), Allocatable, Target :: &
       EVAC_SSTANDS

  ! Holds the information of the EXIT-lines.
  Type (EVAC_EXIT_Type), Dimension(:), Allocatable, Target ::  &
       EVAC_EXITS

  ! Holds the information of the DOOR-lines.
  Type (EVAC_DOOR_Type), Dimension(:), Allocatable, Target ::  &
       EVAC_DOORS

  ! Holds the information of the ENTR-lines.
  Type (EVAC_ENTR_Type), Dimension(:), Allocatable, Target ::  &
       EVAC_ENTRYS

  ! Holds the information of the CORR-lines.
  Type (EVAC_CORR_Type), Dimension(:), Allocatable, Target ::  &
       EVAC_CORRS

  ! Holds the information of the PERS-lines.
  Type (EVAC_PERS_Type), Dimension(:), Allocatable, Target ::  &
       EVAC_PERSON_CLASSES
  !
  ! Next are needed for the Gaussian random numbers
  Integer GaussFlag
  Real(EB) GaussSet1, GaussSet2
  Integer GTrunFlag
  Real(EB) GTrunSet1, GTrunSet2
  !
  Integer :: NPC_EVAC, NPC_PERS, N_EXITS, N_DOORS, N_ENTRYS, &
       N_CORRS, N_EGRIDS, N_NODES, N_HOLES, N_SSTANDS
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
       FC_DAMPING, EVAC_DT_MIN, V_MAX, V_ANGULAR_MAX, V_ANGULAR
  !
  Real(EB), Dimension(:), Allocatable :: Tsteps
  !
  Integer :: n_dead, icyc_old
  Real(EB) :: fed_max_alive, fed_max
  !
Contains
  !
  Subroutine READ_EVAC
    Implicit None

    !
    Integer :: NUMBER_INITIAL_PERSONS, &
         SAMPLING_FACTOR, IPC, n_tmp, GN_MIN, GN_MAX
    Real(EB) :: DTSAM
    Logical :: EVACFILE

    Real(EB) :: DUMMY
    Real(EB) :: XB(6), XB1(6), XB2(6)
    Real(EB), Dimension(3) :: XYZ, XYZ_SMOKE
    Integer :: IOS, IZERO, N, I, IOR, j
    Character(30) QUANTITY
    Character(60) FYI,ID,PERS_ID,TO_NODE,EVAC_ID, &
         DEFAULT_PROPERTIES, RGB_REVA
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
    Integer :: MAX_HUMANS_INSIDE, n_max_in_corrs, COLOR_INDEX
    Real(EB) :: MAX_FLOW, WIDTH, T_START, T_STOP, WIDTH1, &
         WIDTH2, EFF_WIDTH, EFF_LENGTH, FAC_SPEED, &
         TIME_OPEN, TIME_CLOSE
    Logical :: CHECK_FLOW, COUNT_ONLY, AFTER_REACTION_TIME, &
         EXIT_SIGN, KEEP_XY
    Logical :: OUTPUT_SPEED, OUTPUT_MOTIVE_FORCE, OUTPUT_FED, OUTPUT_OMEGA,&
         OUTPUT_ANGLE, OUTPUT_CONTACT_FORCE, OUTPUT_TOTAL_FORCE
    Integer, Dimension(3) :: RGB_DEAD
    Character(26) :: VENT_FFIELD, MESH_ID
    Real(EB) :: FAC_V0_UP, FAC_V0_DOWN, FAC_V0_HORI, HEIGHT, HEIGHT0, ESC_SPEED

    Character(26), Dimension(51) :: KNOWN_DOOR_NAMES
    Real(EB), Dimension(51) :: KNOWN_DOOR_PROBS

    Type (MESH_Type), Pointer :: M
    Integer :: ii,jj,kk
    
    Integer :: size_rnd
    Integer, Dimension(8) :: t_rnd
    Integer, Dimension(:), Allocatable :: seed_rnd

    Namelist /EXIT/ ID, XB, IOR, FLOW_FIELD_ID, CHECK_FLOW, &
         MAX_FLOW, FYI, COUNT_ONLY, WIDTH, XYZ, VENT_FFIELD, &
         MESH_ID, COLOR_INDEX, XYZ_SMOKE, &
         TIME_OPEN, TIME_CLOSE
    Namelist /DOOR/ ID, XB, IOR, FLOW_FIELD_ID, CHECK_FLOW, &
         MAX_FLOW, TO_NODE, FYI, WIDTH, XYZ, VENT_FFIELD, &
         EXIT_SIGN, MESH_ID, COLOR_INDEX, XYZ_SMOKE, KEEP_XY, &
         TIME_OPEN, TIME_CLOSE
    Namelist /ENTR/ ID, XB, IOR, FLOW_FIELD_ID, MAX_FLOW, &
         FYI, WIDTH, QUANTITY, PERS_ID, T_START, &
         T_STOP, AFTER_REACTION_TIME, &
         KNOWN_DOOR_NAMES, KNOWN_DOOR_PROBS, &
         MESH_ID, COLOR_INDEX
    Namelist /CORR/ ID, XB, IOR, FLOW_FIELD_ID, CHECK_FLOW, &
         MAX_FLOW, TO_NODE, FYI, WIDTH, WIDTH1, WIDTH2, &
         EFF_WIDTH, EFF_LENGTH, MAX_HUMANS_INSIDE, FAC_SPEED, &
         XB1, XB2
    Namelist /EVAC/ NUMBER_INITIAL_PERSONS, QUANTITY, FYI, &
         ID, DTSAM, XB, FLOW_FIELD_ID, PERS_ID, &
         T_START, T_STOP, IOR, MAX_FLOW, WIDTH, ANGLE, &
         AFTER_REACTION_TIME, GN_MIN, GN_MAX, &
         KNOWN_DOOR_NAMES, KNOWN_DOOR_PROBS, MESH_ID, COLOR_INDEX
    Namelist /EVHO/ FYI, ID, XB, EVAC_ID, PERS_ID, MESH_ID

    Namelist /EVSS/ FYI, ID, XB, MESH_ID, HEIGHT, HEIGHT0, IOR, &
         FAC_V0_UP, FAC_V0_DOWN, FAC_V0_HORI, ESC_SPEED

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
         OUTPUT_ANGLE, OUTPUT_CONTACT_FORCE, OUTPUT_TOTAL_FORCE, &
         COLOR_INDEX, RGB_DEAD, RGB_REVA
    !
    NPPS = 30000 ! Number Persons Per Set (dump to a file)
    !
    EVAC_DT = EVAC_DT_FLOWFIELD     ! Initialize the clock
    EVAC_CLOCK = 0.0_EB    ! clock for the CHID_evac.csv file
    EVAC_N_QUANTITIES = 0

    i33 = 0
    ilh = 0

    Allocate(Tsteps(NMESHES),STAT=IZERO)
    Call ChkMemErr('READ','Tsteps',IZERO) 
    Tsteps(:) = EVAC_DT_FLOWFIELD


    !
    ! I_EVAC: 'binary' index:
    ! These are just initialization. Later it is checked if files
    ! exists. If EFF/FED file does not exists, then it is calculated 
    ! (and saved).
    I_EVAC = 16*0 + 8*0 + 4*0 + 2*0 + 1*0 ! do not save soot,fed files
    If (.Not. All(EVACUATION_ONLY) ) Then
       ! There are fire grids ==> save fed and evac flow fields
       I_EVAC = 16*1 + 8*0 + 4*0 + 2*1 + 1*1
    Else
       ! There are no fire grids ==> try to read fed and evac flow fields
       I_EVAC = 16*0 + 8*1 + 4*1 + 2*0 + 1*0
    End If
    !
    ! Every human has an identification number, ILABEL_last is
    ! the last used number, so next human will have an
    ! identification number, which is ILABEL_last + 1
    ILABEL_last = 0
    !
    !
    ! Determine total number of EVAC lines in the input file
    !
    Open(LU_INPUT,File=FN_INPUT)

    NPC_EVAC = 0
    COUNT_EVAC_LOOP: Do
       Call CHECKREAD('EVAC',LU_INPUT,IOS) 
       If (IOS == 1) Exit COUNT_EVAC_LOOP
       Read(LU_INPUT,NML=EVAC,End=219,ERR=220,IOSTAT=IOS)
       NPC_EVAC = NPC_EVAC + 1
220    If (IOS > 0) Call SHUTDOWN('ERROR: Problem with EVAC line')
    End Do COUNT_EVAC_LOOP
219 Rewind(LU_INPUT)
    !
    ! Determine total number of PERS lines in the input file
    !
    NPC_PERS = 0
    COUNT_PERS_LOOP: Do
       Call CHECKREAD('PERS',LU_INPUT,IOS) 
       If (IOS == 1) Exit COUNT_PERS_LOOP
       Read(LU_INPUT,NML=PERS,End=221,ERR=222,IOSTAT=IOS)
       NPC_PERS = NPC_PERS + 1
222    If (IOS > 0) Call SHUTDOWN('ERROR: Problem with PERS line')
    End Do COUNT_PERS_LOOP
221 Rewind(LU_INPUT)
    !
    ! Determine total number of EXIT lines in the input file
    !
    N_EXITS = 0
    COUNT_EXITS_LOOP: Do
       Call CHECKREAD('EXIT',LU_INPUT,IOS) 
       If (IOS == 1) Exit COUNT_EXITS_LOOP
       Read(LU_INPUT,NML=Exit,End=223,ERR=224,IOSTAT=IOS)
       N_EXITS = N_EXITS + 1
224    If (IOS > 0) Call SHUTDOWN('ERROR: Problem with EXIT line')
    End Do COUNT_EXITS_LOOP
223 Rewind(LU_INPUT)
    !
    ! Determine total number of DOOR lines in the input file
    !
    N_DOORS = 0
    COUNT_DOORS_LOOP: Do
       Call CHECKREAD('DOOR',LU_INPUT,IOS) 
       If (IOS == 1) Exit COUNT_DOORS_LOOP
       Read(LU_INPUT,NML=DOOR,End=225,ERR=226,IOSTAT=IOS)
       N_DOORS = N_DOORS + 1
226    If (IOS > 0) Call SHUTDOWN('ERROR: Problem with DOOR line')
    End Do COUNT_DOORS_LOOP
225 Rewind(LU_INPUT)
    !
    ! Determine total number of ENTR lines in the input file
    !
    N_ENTRYS = 0
    COUNT_ENTRYS_LOOP: Do
       Call CHECKREAD('ENTR',LU_INPUT,IOS) 
       If (IOS == 1) Exit COUNT_ENTRYS_LOOP
       Read(LU_INPUT,NML=ENTR,End=227,ERR=228,IOSTAT=IOS)
       N_ENTRYS = N_ENTRYS + 1
228    If (IOS > 0) Call SHUTDOWN('ERROR: Problem with ENTR line')
    End Do COUNT_ENTRYS_LOOP
227 Rewind(LU_INPUT)
    !
    ! Determine total number of CORR lines in the input file
    !
    N_CORRS = 0
    COUNT_CORRS_LOOP: Do
       Call CHECKREAD('CORR',LU_INPUT,IOS) 
       If (IOS == 1) Exit COUNT_CORRS_LOOP
       Read(LU_INPUT,NML=CORR,End=229,ERR=230,IOSTAT=IOS)
       N_CORRS = N_CORRS + 1
230    If (IOS > 0) Call SHUTDOWN('ERROR: Problem with CORR line')
    End Do COUNT_CORRS_LOOP
229 Rewind(LU_INPUT)
    !
    ! Determine total number of EVHO lines in the input file
    !
    N_HOLES = 0
    COUNT_EVHO_LOOP: Do
       Call CHECKREAD('EVHO',LU_INPUT,IOS) 
       If (IOS == 1) Exit COUNT_EVHO_LOOP
       Read(LU_INPUT,NML=EVHO,End=231,ERR=232,IOSTAT=IOS)
       N_HOLES = N_HOLES + 1
232    If (IOS > 0) Call SHUTDOWN('ERROR: Problem with EVHO line')
    End Do COUNT_EVHO_LOOP
231 Rewind(LU_INPUT)
    !
    ! Determine total number of EVSS lines in the input file
    !
    N_SSTANDS = 0
    COUNT_EVSS_LOOP: Do
       Call CHECKREAD('EVSS',LU_INPUT,IOS) 
       If (IOS == 1) Exit COUNT_EVSS_LOOP
       Read(LU_INPUT,NML=EVSS,End=233,ERR=234,IOSTAT=IOS)
       N_SSTANDS = N_SSTANDS + 1
234    If (IOS > 0) Call SHUTDOWN('ERROR: Problem with EVSS line')
    End Do COUNT_EVSS_LOOP
233 Rewind(LU_INPUT)
    !
    ! Allocate quantities for EVAC, PERS, EXIT types
    !
    If (npc_evac > 0 ) Then
       Allocate(EVACUATION(NPC_EVAC),STAT=IZERO)
       Call ChkMemErr('READ','EVACUATION',IZERO) 
    End If

    If (n_holes > 0 ) Then
       Allocate(EVAC_HOLES(N_HOLES),STAT=IZERO)
       Call ChkMemErr('READ','EVAC_HOLES',IZERO) 
    End If

    If (n_sstands > 0 ) Then
       Allocate(EVAC_SSTANDS(N_SSTANDS),STAT=IZERO)
       Call ChkMemErr('READ','EVAC_SSTANDS',IZERO) 
    End If

    If (n_exits > 0 ) Then
       Allocate(EVAC_EXITS(N_EXITS),STAT=IZERO)
       Call ChkMemErr('READ','EVAC_EXITS',IZERO) 
    End If

    If (n_doors > 0 ) Then
       Allocate(EVAC_DOORS(N_DOORS),STAT=IZERO)
       Call ChkMemErr('READ','EVAC_DOORS',IZERO) 
    End If

    If (n_entrys > 0 ) Then
       Allocate(EVAC_ENTRYS(N_ENTRYS),STAT=IZERO)
       Call ChkMemErr('READ','EVAC_ENTRYS',IZERO) 
    End If

    If (n_corrs > 0 ) Then
       Allocate(EVAC_CORRS(N_CORRS),STAT=IZERO)
       Call ChkMemErr('READ','EVAC_CORRS',IZERO) 
    End If

    Allocate(EVAC_PERSON_CLASSES(0:NPC_PERS),STAT=IZERO)
    Call ChkMemErr('READ','EVAC_PERSON_CLASSES',IZERO) 

    n_egrids = 0
    Do n = 1, nmeshes
       If (evacuation_only(n) .And. evacuation_grid(n) ) Then
          n_egrids = n_egrids + 1
       End If
    End  Do

    n_nodes = n_entrys + n_exits + n_doors + n_corrs + n_egrids
    If (n_nodes > 0 ) Then
       Allocate(EVAC_Node_List(1:n_nodes),STAT=IZERO)
       Call ChkMemErr('READ','EVAC_NODE_LIST',IZERO) 
    End If

    If (npc_evac > 0 ) Then
       EVACUATION(1:NPC_EVAC)%COLOR_INDEX = 0
       EVACUATION(1:NPC_EVAC)%GRID_NAME   = 'null'
       EVACUATION(1:NPC_EVAC)%CLASS_NAME  = 'null'
       EVACUATION(1:NPC_EVAC)%IMESH       = 0
       EVACUATION(1:NPC_EVAC)%ID_NAME     = 'null'
    End If
    If (n_holes > 0 ) Then
       EVAC_HOLES(1:N_HOLES)%GRID_NAME   = 'null'
       EVAC_HOLES(1:N_HOLES)%PERS_ID     = 'null'
       EVAC_HOLES(1:N_HOLES)%EVAC_ID     = 'null'
       EVAC_HOLES(1:N_HOLES)%IMESH       = 0
       EVAC_HOLES(1:N_HOLES)%ID_NAME     = 'null'
    End If

    EVAC_PERSON_CLASSES(0:NPC_PERS)%ID_NAME = 'null'

    If (n_exits > 0 ) Then
       EVAC_EXITS(1:N_EXITS)%ID_NAME   = 'null'
       EVAC_EXITS(1:N_EXITS)%TO_NODE   = 'null'
       EVAC_EXITS(1:N_EXITS)%GRID_NAME = 'null'
       EVAC_EXITS(1:N_EXITS)%IMESH     = 0
       EVAC_EXITS(1:N_EXITS)%COLOR_INDEX = 0
    End If

    If (n_doors > 0 ) Then
       EVAC_DOORS(1:N_DOORS)%ID_NAME   = 'null'
       EVAC_DOORS(1:N_DOORS)%TO_NODE   = 'null'
       EVAC_DOORS(1:N_DOORS)%GRID_NAME = 'null'
       EVAC_DOORS(1:N_DOORS)%IMESH     = 0
       EVAC_DOORS(1:N_DOORS)%IMESH2    = 0
       EVAC_DOORS(1:N_DOORS)%COLOR_INDEX = 0
    End If

    If (n_corrs > 0 ) Then
       EVAC_CORRS(1:N_CORRS)%ID_NAME   = 'null'
       EVAC_CORRS(1:N_CORRS)%TO_NODE   = 'null'
       EVAC_CORRS(1:N_CORRS)%GRID_NAME = 'null'
       EVAC_CORRS(1:N_CORRS)%IMESH     = 0
       EVAC_CORRS(1:N_CORRS)%IMESH2    = 0
    End If

    If (n_entrys > 0 ) Then
       EVAC_ENTRYS(1:N_ENTRYS)%ID_NAME     = 'null'
       EVAC_ENTRYS(1:N_ENTRYS)%TO_NODE     = 'null'
       EVAC_ENTRYS(1:N_ENTRYS)%GRID_NAME   = 'null'
       EVAC_ENTRYS(1:N_ENTRYS)%CLASS_NAME  = 'null'
       EVAC_ENTRYS(1:N_ENTRYS)%IMESH       = 0
       EVAC_ENTRYS(1:N_ENTRYS)%COLOR_INDEX = 0
    End If

    !
    ! NEXT PARAMETERS ARE SAME FOR ALL HUMANS. THE LAST
    ! VALUES READ IN FROM 'PERS' LINES ARE VALID.
    FAC_A_WALL  = 1.0_EB
    FAC_B_WALL  = 0.5_EB
    LAMBDA_WALL = 0.2_EB
    NOISEME     = 0.0_EB
    NOISETH     = 0.01_EB
    NOISECM     = 3.0_EB
    I_FRIC_SW   = 1
    FC_DAMPING        = 500.0_EB
    V_MAX             = 20.0_EB
    V_ANGULAR_MAX     = 8.0_EB  ! rps
    V_ANGULAR         = 2.0_EB  ! rps
    GROUP_EFF         = 0.0_EB
    RADIUS_COMPLETE_0 = 0.2_EB
    RADIUS_COMPLETE_1 = 0.5_EB
    NOT_RANDOM = .False.
    ! Which doors are 'smoke free'
    FED_DOOR_CRIT = 0.000001_EB
    ! How to color humans?
    COLOR_METHOD = -1 ! Default is standard human colors in Smokeview
    ! Smoke is detected, when its density is larger than, e.g. 1 mg/m3
    ! Default is no detection due to smoke.
    TDET_SMOKE_DENS = -999.9_EB  ! 
    DENS_INIT       = 0.0_EB
    EVAC_DT_MAX     = 0.01_EB
    EVAC_DT_MIN     = 0.001_EB
    GROUP_DENS      = 0.0_EB
    OUTPUT_SPEED         = .FALSE.
    OUTPUT_MOTIVE_FORCE  = .FALSE.
    OUTPUT_FED           = .FALSE.
    OUTPUT_OMEGA         = .FALSE.
    OUTPUT_ANGLE         = .FALSE.
    OUTPUT_CONTACT_FORCE = .FALSE.
    OUTPUT_TOTAL_FORCE   = .FALSE.
    ! 
    ! Read the PERS lines (no read for default n=0 case)
    !
    READ_PERS_LOOP: Do N=0,NPC_PERS
       PCP=>EVAC_PERSON_CLASSES(N)
       !
       ID  = 'null'
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
       DET_LOW = 0.0_EB
       TAU_LOW = 0.0_EB
       VEL_HIGH = 999.0_EB
       DIA_HIGH = 999.0_EB
       PRE_HIGH = Huge(PRE_HIGH)
       DET_HIGH = Huge(PRE_HIGH)
       TAU_HIGH = 999.0_EB
       ! Default values for persons
       VEL_MEAN = 1.25_EB
       DIA_MEAN = 0.51_EB
       PRE_MEAN = 10.0_EB
       DET_MEAN = 0.0_EB
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
       !
       ! No read for default values
       If ( N > 0 ) Then
          Call CHECKREAD('PERS',LU_INPUT,IOS)
          If (IOS == 1) Exit READ_PERS_LOOP
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
             Write(MESSAGE,'(A,I4,A)') &
                  'ERROR: PERS',N,' problem with DEFAULT_PROPERTIES'
             Call SHUTDOWN(MESSAGE)
          End Select

       End If

       DIAMETER_DIST = Max(0,DIAMETER_DIST)
       VELOCITY_DIST = Max(0,VELOCITY_DIST)
       TAU_EVAC_DIST = Max(0,TAU_EVAC_DIST)
       !
       !
       ! 
       PCP%ID_NAME = ID
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
          PCP%m_iner = -M_INERTIA*(0.25_EB*PCP%D_mean**2+PCP%r_torso**2)**2 / &
               (0.27_EB**2+0.16_EB**2)**2
       Else
          PCP%m_iner = M_INERTIA  ! kg m2
       End If
       !
    End Do READ_PERS_LOOP
24  Rewind(LU_INPUT)

    If (GROUP_DENS .Le. 0.01_EB) GROUP_DENS = 0.25_EB
    If (GROUP_DENS .Gt. 3.50_EB) GROUP_DENS = 3.50_EB
    DENS_INIT = Max(GROUP_DENS,DENS_INIT)
    If (TDET_SMOKE_DENS < 0.0_EB) TDET_SMOKE_DENS = Huge(TDET_SMOKE_DENS)

    If (.Not. NOT_RANDOM ) Then    ! Initialize the generator randomly
       Call Random_Seed(size=size_rnd)
       Allocate(seed_rnd(size_rnd),STAT=IZERO)
       Call ChkMemErr('READ_EVAC','seed_rnd',IZERO)
       Call Date_and_Time(values = t_rnd)
       seed_rnd = 100.0_EB*t_rnd(7) + t_rnd(8)/10.0_EB
       Call Random_Seed(put=seed_rnd)
       Deallocate(seed_rnd)
    End If

    ! Add the number of the evacuation classes
    N_EVAC = 1

    Allocate(EVAC_CLASS_NAME(N_EVAC),STAT=IZERO)
    Call ChkMemErr('READ_EVAC','EVAC_CLASS_NAME',IZERO)
    Allocate(EVAC_CLASS_RGB(3,N_EVAC),STAT=IZERO)
    Call ChkMemErr('READ_EVAC','EVAC_CLASS_RGB',IZERO)

    EVAC_CLASS_NAME(1) = 'Human'
    Do N = 1, N_EVAC
       EVAC_CLASS_RGB(1:3,N) = (/ 39, 64,139/)  ! ROYAL BLUE4
    End Do

    EVAC_AVATAR_NCOLOR = 7
    Allocate(EVAC_AVATAR_RGB(3,MAX(7,EVAC_AVATAR_NCOLOR)),STAT=IZERO)
    Call ChkMemErr('READ_EVAC','EVAC_AVATAR_RGB',IZERO)
    ! Default values for Avatar color index array
    EVAC_AVATAR_RGB(1:3,1) = (/  0,  0,  0/)  ! black
    EVAC_AVATAR_RGB(1:3,2) = (/255,255,  0/)  ! yellow
    EVAC_AVATAR_RGB(1:3,3) = (/  0,  0,255/)  ! blue
    EVAC_AVATAR_RGB(1:3,4) = (/255,  0,  0/)  ! red
    EVAC_AVATAR_RGB(1:3,5) = (/  0,255,  0/)  ! green
    EVAC_AVATAR_RGB(1:3,6) = (/255,  0,255/)  ! magenta
    EVAC_AVATAR_RGB(1:3,7) = (/  0,255,255/)  ! cyan

    !
    ! Read the EXIT lines
    !
    READ_EXIT_LOOP: Do N = 1, N_EXITS
       PEX=>EVAC_EXITS(N)
       !
       ID            = 'null'
       XB            = 0.0_EB
       IOR           = 0
       FLOW_FIELD_ID = 'null'
       VENT_FFIELD   = 'null'
       MESH_ID     = 'null'
       CHECK_FLOW    = .False.
       COUNT_ONLY    = .False.
       MAX_FLOW      = 0.0_EB
       WIDTH         = 0.0_EB
       TIME_OPEN     = 0.0_EB
       TIME_CLOSE    = Huge(TIME_CLOSE)
       XYZ(:)        = Huge(XYZ)
       XYZ_SMOKE(:)  = Huge(XYZ_SMOKE)
       COLOR_INDEX   = 0
       !
       Call CHECKREAD('EXIT',LU_INPUT,IOS)
       If (IOS == 1) Exit READ_EXIT_LOOP
       Read(LU_INPUT,Exit,End=26,IOSTAT=IOS)
       !
       Do I=1,5,2
          If (XB(I) > XB(I+1)) Then
             DUMMY   = XB(I)
             XB(I)   = XB(I+1)
             XB(I+1) = DUMMY
          End If
       End Do
       !
       PEX%X1 = XB(1)
       PEX%X2 = XB(2)
       PEX%Y1 = XB(3)
       PEX%Y2 = XB(4)
       PEX%Z1 = XB(5)
       PEX%Z2 = XB(6)
       PEX%IOR = IOR
       PEX%ID_NAME    = Trim(ID)
       PEX%GRID_NAME  = Trim(FLOW_FIELD_ID)
       PEX%CHECK_FLOW = CHECK_FLOW
       PEX%VENT_FFIELD= VENT_FFIELD
       PEX%INODE      = 0
       PEX%T_first    = 0.0_EB
       PEX%T_last     = 0.0_EB
       PEX%ICOUNT     = 0
       PEX%Flow_max   = 0.0_EB
       PEX%TIME_OPEN  = TIME_OPEN
       PEX%TIME_CLOSE = TIME_CLOSE
       If (CHECK_FLOW) PEX%Flow_max   = MAX_FLOW
       PEX%COUNT_ONLY = .False.
       If (COUNT_ONLY) PEX%COUNT_ONLY = .True.

       PEX%COLOR_INDEX = Mod(Max(0,COLOR_INDEX),8) ! 0-7 always

       PEX%FED_MESH = 0
       If (XYZ(1) < Huge(XYZ)) Then
          PEX%X = XYZ(1)
          PEX%Y = XYZ(2)
          PEX%Z = XYZ(3)
       Else
          PEX%X = 0.5_EB*(XB(1)+XB(2))
          PEX%Y = 0.5_EB*(XB(3)+XB(4))
          PEX%Z = 0.5_EB*(XB(5)+XB(6))
       End If

       If (XYZ_SMOKE(1) < Huge(XYZ_SMOKE)) Then
          PEX%Xsmoke = XYZ_SMOKE(1)
          PEX%Ysmoke = XYZ_SMOKE(2)
          PEX%Zsmoke = XYZ_SMOKE(3)
       Else
          PEX%Xsmoke = PEX%X
          PEX%Ysmoke = PEX%Y
          PEX%Zsmoke = PEX%Z
       End If

       Select Case (IOR)
       Case (-1,+1)
          If (WIDTH <= 0.0_EB) Then
             PEX%Width = XB(4) - XB(3)
          Else
             PEX%Width = WIDTH
          End If
       Case (-2,+2)
          If (WIDTH <= 0.0_EB) Then
             PEX%Width = XB(2) - XB(1)
          Else
             PEX%Width = WIDTH
          End If
       Case (-3)
          If ( (XB(4)-XB(3)) <= 0.0_EB .Or. (XB(2)-XB(1)) <= 0.0_EB) Then
             Write(MESSAGE,'(A,I4,A)') &
                  'ERROR: EXIT',N,' IOR=-3 but not 3-dim object'
             Call SHUTDOWN(MESSAGE)
          End If
       Case (0)
          If ( (XB(4)-XB(3)) <= 0.0_EB .Or. (XB(2)-XB(1)) <= 0.0_EB) Then
             Write(MESSAGE,'(A,I4,A)') &
                  'ERROR: EXIT',N,' no IOR but not 3-dim object'
             Call SHUTDOWN(MESSAGE)
          End If
       Case Default
          Write(MESSAGE,'(A,I4,A)') &
               'ERROR: EXIT',N,' problem with IOR'
          Call SHUTDOWN(MESSAGE)
       End Select
       ! 
       ! Check which evacuation floor
       ii = 0
       PEX_MeshLoop: Do i = 1, nmeshes
          If (evacuation_only(i) .And. evacuation_grid(i)) Then
             If ( (PEX%Z1 >= Meshes(i)%ZS .And. PEX%Z2 <= Meshes(i)%ZF).And. &
                  (PEX%Y1 >= Meshes(i)%YS .And. PEX%Y2 <= Meshes(i)%YF).And. &
                  (PEX%X1 >= Meshes(i)%XS .And. PEX%X2 <= Meshes(i)%XF)) Then
                If (Trim(MESH_ID) == 'null' .Or. &
                     Trim(MESH_ID) == Trim(MESH_NAME(i))) Then
                   ii = ii + 1
                   PEX%IMESH = i
                   !cc             Exit PEX_MeshLoop
                End If
             End If
          End If
       End Do PEX_MeshLoop
       If (PEX%IMESH == 0) Then
          Write(MESSAGE,'(A,A,A)') &
               'ERROR: EXIT line ',Trim(PEX%ID_NAME), &
               ' problem with IMESH, no mesh found'
          Call SHUTDOWN(MESSAGE)
       End If
       If (ii > 1) Then
          Write(MESSAGE,'(A,A,A)') &
               'ERROR: EXIT line ',Trim(PEX%ID_NAME), &
               ' not an unique mesh found '
          Call SHUTDOWN(MESSAGE)
       End If
       ! 
       ! Check which vent field. If VENT_FFIELD is not found, use
       ! the main evac grid.
       PEX%I_VENT_FFIELD = 0
       PEX_Mesh2Loop: Do i = 1, nmeshes
          If ( evacuation_only(i) .And. &
               (Trim(MESH_NAME(i)) == Trim(PEX%VENT_FFIELD)) ) Then
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
          Write(MESSAGE,'(A,A,A)') &
               'ERROR: EXIT line ',Trim(PEX%ID_NAME), &
               ' problem with XYZ, no mesh found'
          Call SHUTDOWN(MESSAGE)
       End If
       ! 
       ! Check, which fire grid and i,j,k (xyz)
       PEX_SmokeLoop: Do i = 1, nmeshes
          If (.Not. evacuation_only(i)) Then
             If ( (PEX%Zsmoke >= Meshes(i)%ZS .And. &
                  PEX%Zsmoke <= Meshes(i)%ZF).And. &
                  (PEX%Ysmoke >= Meshes(i)%YS .And. &
                  PEX%Ysmoke <= Meshes(i)%YF).And. &
                  (PEX%Xsmoke >= Meshes(i)%XS .And. &
                  PEX%Xsmoke <= Meshes(i)%XF)) Then
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
26  Rewind(LU_INPUT)
    !
    ! Read the DOOR lines
    !
    READ_DOOR_LOOP: Do N = 1, N_DOORS
       PDX=>EVAC_DOORS(N)
       !
       ID            = 'null'
       XB            = 0.0_EB
       IOR           = 0
       FLOW_FIELD_ID = 'null'
       VENT_FFIELD   = 'null'
       MESH_ID     = 'null'
       TO_NODE       = 'null'
       CHECK_FLOW    = .False.
       EXIT_SIGN     = .False.
       MAX_FLOW      = 0.0_EB
       WIDTH         = 0.0_EB
       TIME_OPEN     = 0.0_EB
       TIME_CLOSE    = Huge(TIME_CLOSE)
       XYZ(:)        = Huge(XYZ)
       XYZ_SMOKE(:)  = Huge(XYZ_SMOKE)
       COLOR_INDEX   = 0
       KEEP_XY       = .False.
       !
       Call CHECKREAD('DOOR',LU_INPUT,IOS)
       If (IOS == 1) Exit READ_DOOR_LOOP
       Read(LU_INPUT,DOOR,End=27,IOSTAT=IOS)
       !
       Do I=1,5,2
          If (XB(I) > XB(I+1)) Then
             DUMMY   = XB(I)
             XB(I)   = XB(I+1)
             XB(I+1) = DUMMY
          End If
       End Do
       !
       PDX%X1 = XB(1)
       PDX%X2 = XB(2)
       PDX%Y1 = XB(3)
       PDX%Y2 = XB(4)
       PDX%Z1 = XB(5)
       PDX%Z2 = XB(6)
       PDX%IOR        = IOR
       PDX%ID_NAME    = ID
       PDX%GRID_NAME  = FLOW_FIELD_ID
       PDX%VENT_FFIELD= VENT_FFIELD
       PDX%CHECK_FLOW = CHECK_FLOW
       PDX%EXIT_SIGN  = EXIT_SIGN
       PDX%KEEP_XY    = KEEP_XY
       PDX%TO_NODE    = TO_NODE
       PDX%INODE      = 0
       PDX%INODE2     = 0
       PDX%T_first    = 0.0_EB
       PDX%T_last     = 0.0_EB
       PDX%ICOUNT     = 0
       PDX%Flow_max   = 0.0_EB
       PDX%TIME_OPEN  = TIME_OPEN
       PDX%TIME_CLOSE = TIME_CLOSE
       If (CHECK_FLOW) PDX%Flow_max   = MAX_FLOW

       PDX%COLOR_INDEX = Mod(Max(0,COLOR_INDEX),8) ! 0-7 always

       PDX%FED_MESH = 0
       If (XYZ(1) < Huge(XYZ)) Then
          PDX%X = XYZ(1)
          PDX%Y = XYZ(2)
          PDX%Z = XYZ(3)
       Else
          PDX%X = 0.5_EB*(XB(1)+XB(2))
          PDX%Y = 0.5_EB*(XB(3)+XB(4))
          PDX%Z = 0.5_EB*(XB(5)+XB(6))
       End If

       If (XYZ_SMOKE(1) < Huge(XYZ_SMOKE)) Then
          PDX%Xsmoke = XYZ(1)
          PDX%Ysmoke = XYZ(2)
          PDX%Zsmoke = XYZ(3)
       Else
          PDX%Xsmoke = PDX%X
          PDX%Ysmoke = PDX%Y
          PDX%Zsmoke = PDX%Z
       End If

       Select Case (IOR)
       Case (-1,+1)
          If (WIDTH <= 0.0_EB) Then
             PDX%Width = XB(4) - XB(3)
          Else
             PDX%Width = WIDTH
          End If
       Case (-2,+2)
          If (WIDTH <= 0.0_EB) Then
             PDX%Width = XB(2) - XB(1)
          Else
             PDX%Width = WIDTH
          End If
       Case (-3)
          If ( (XB(4)-XB(3)) <= 0.0_EB .Or. (XB(2)-XB(1)) <= 0.0_EB) Then
             Write(MESSAGE,'(A,I4,A)') &
                  'ERROR: DOOR',N,' IOR=-3 but not 3-dim object'
             Call SHUTDOWN(MESSAGE)
          End If
       Case (0)
          If ( (XB(4)-XB(3)) <= 0.0_EB .Or. (XB(2)-XB(1)) <= 0.0_EB) Then
             Write(MESSAGE,'(A,I4,A)') &
                  'ERROR: DOOR',N,' no IOR but not 3-dim object'
             Call SHUTDOWN(MESSAGE)
          End If
       Case Default
          Write(MESSAGE,'(A,I4,A)') &
               'ERROR: DOOR',N,' problem with IOR'
          Call SHUTDOWN(MESSAGE)
       End Select
       ! 
       ! Check which evacuation floor
       ! Now there may be overlapping meshes.
       ii = 0
       PDX_MeshLoop: Do i = 1, nmeshes
          If (evacuation_only(i) .And. evacuation_grid(i)) Then
             If ( (PDX%Z1 >= Meshes(i)%ZS .And. PDX%Z2 <= Meshes(i)%ZF).And. &
                  (PDX%Y1 >= Meshes(i)%YS .And. PDX%Y2 <= Meshes(i)%YF).And. &
                  (PDX%X1 >= Meshes(i)%XS .And. PDX%X2 <= Meshes(i)%XF)) Then
                If (Trim(MESH_ID) == 'null' .Or. &
                     Trim(MESH_ID) == Trim(MESH_NAME(i))) Then
                   ii = ii + 1
                   PDX%IMESH = i
                End If
             End If
          End If
       End Do PDX_MeshLoop
       If (PDX%IMESH == 0) Then
          Write(MESSAGE,'(A,A,A)') &
               'ERROR: DOOR line ',Trim(PDX%ID_NAME), &
               ' problem with IMESH, no mesh found'
          Call SHUTDOWN(MESSAGE)
       End If
       If (ii > 1) Then
          Write(MESSAGE,'(A,A,A)') &
               'ERROR: DOOR line ',Trim(PDX%ID_NAME), &
               ' not an unique mesh found '
          Call SHUTDOWN(MESSAGE)
       End If
       ! 
       ! Check which vent field. If VENT_FFIELD is not found, use
       ! the main evac grid.
       PDX%I_VENT_FFIELD = 0
       PDX_Mesh2Loop: Do i = 1, nmeshes
          If ( evacuation_only(i) .And. &
               (Trim(MESH_NAME(i)) == Trim(PDX%VENT_FFIELD)) ) Then
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
          Write(MESSAGE,'(A,A,A)') &
               'ERROR: DOOR line ',Trim(PDX%ID_NAME), &
               ' problem with XYZ, no mesh found'
          Call SHUTDOWN(MESSAGE)
       End If
       ! 
       ! Check, which fire grid and i,j,k (xyz)
       PDX_SmokeLoop: Do i = 1, nmeshes
          If (.Not. evacuation_only(i)) Then
             If ( (PDX%Zsmoke >= Meshes(i)%ZS .And. &
                  PDX%Zsmoke <= Meshes(i)%ZF).And. &
                  (PDX%Ysmoke >= Meshes(i)%YS .And. &
                  PDX%Ysmoke <= Meshes(i)%YF).And. &
                  (PDX%Xsmoke >= Meshes(i)%XS .And. &
                  PDX%Xsmoke <= Meshes(i)%XF)) Then
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
27  Rewind(LU_INPUT)
    !
    ! Read the CORR line
    !
    n_max_in_corrs = 0
    READ_CORR_LOOP: Do N = 1, N_CORRS
       PCX=>EVAC_CORRS(N)
       !
       ID            = 'null'
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
       If (IOS == 1) Exit READ_CORR_LOOP
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
       PCX%ID_NAME    = ID
       PCX%GRID_NAME  = FLOW_FIELD_ID
       PCX%CHECK_FLOW = CHECK_FLOW
       PCX%TO_NODE    = TO_NODE
       PCX%INODE      = 0
       PCX%INODE2     = 0
       PCX%T_first    = 0.0_EB
       PCX%T_last     = 0.0_EB
       PCX%ICOUNT     = 0

       PCX%MAX_HUMANS_INSIDE = 0
       If (MAX_HUMANS_INSIDE > 0 ) Then
          PCX%MAX_HUMANS_INSIDE = MAX_HUMANS_INSIDE
       Else
          Write(MESSAGE,'(A,I4,A)') &
               'ERROR: CORR',N,' MAX_HUMANS_INSIDE <= 0'
          Call SHUTDOWN(MESSAGE)
       End If

       If (FAC_SPEED < 0 ) Then
          Write(MESSAGE,'(A,I4,A)') &
               'ERROR: CORR',N,' FAC_SPEED < 0'
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
          Write(MESSAGE,'(A,I4,A)') &
               'ERROR: CORR',N,' EFF_LENGTH <= 0'
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
       !
    End Do READ_CORR_LOOP
29  Rewind(LU_INPUT)


    ! Now exits, doors, and corrs are already read in
    If (n_nodes > 0 ) Then
       n_tmp = 0
       Do n = 1, nmeshes
          If (evacuation_only(n).And.evacuation_grid(n)) Then
             n_tmp = n_tmp + 1
             EVAC_Node_List(n_tmp)%Node_Index = n_tmp
             EVAC_Node_List(n_tmp)%Node_Type  = 'Floor'
             EVAC_Node_List(n_tmp)%ID_NAME    = MESH_NAME(n)
             EVAC_Node_List(n_tmp)%GRID_NAME  = MESH_NAME(n)
             EVAC_Node_List(n_tmp)%Mesh_index = n
          End If
       End Do
       Do n = 1, n_entrys
          n_tmp = n_tmp + 1
          evac_entrys(n)%INODE             = n_tmp 
          EVAC_Node_List(n_tmp)%Node_Index = n
          EVAC_Node_List(n_tmp)%Node_Type  = 'Entry'
       End Do
       Do n = 1, n_doors
          n_tmp = n_tmp + 1
          evac_doors(n)%INODE              = n_tmp 
          EVAC_Node_List(n_tmp)%Node_Index = n
          EVAC_Node_List(n_tmp)%Node_Type  = 'Door'
          EVAC_Node_List(n_tmp)%ID_NAME    = EVAC_DOORS(n)%ID_NAME
          EVAC_Node_List(n_tmp)%Mesh_Index = EVAC_DOORS(n)%IMESH
       End Do
       Do n = 1, n_exits
          n_tmp = n_tmp + 1
          evac_exits(n)%INODE              = n_tmp 
          EVAC_Node_List(n_tmp)%Node_Index = n
          EVAC_Node_List(n_tmp)%Node_Type  = 'Exit'
          EVAC_Node_List(n_tmp)%ID_NAME    = EVAC_EXITS(n)%ID_NAME
          EVAC_Node_List(n_tmp)%Mesh_Index = EVAC_EXITS(n)%IMESH
       End Do
       Do n = 1, n_corrs
          n_tmp = n_tmp + 1
          evac_corrs(n)%INODE              = n_tmp 
          EVAC_Node_List(n_tmp)%Node_Index = n
          EVAC_Node_List(n_tmp)%Node_Type  = 'Corr'
          EVAC_Node_List(n_tmp)%ID_NAME    = EVAC_CORRS(n)%ID_NAME
       End Do
    End If

    !
    ! Read the ENTR lines
    !
    READ_ENTR_LOOP: Do N = 1, N_ENTRYS
       PNX=>EVAC_ENTRYS(N)
       !
       ID            = 'null'
       XB            = 0.0_EB
       IOR           = 0
       FLOW_FIELD_ID = 'null'
       MESH_ID     = 'null'
       TO_NODE       = 'null'
       PERS_ID       = 'null'
       QUANTITY      = 'null'
       MAX_FLOW      = 0.0_EB
       WIDTH         = 0.0_EB
       AFTER_REACTION_TIME      = .False.
       T_START                  = -99.0_EB
       T_STOP                   = -99.0_EB
       KNOWN_DOOR_NAMES         = 'null'
       KNOWN_DOOR_PROBS         = 1.0_EB
       !
       !
       Call CHECKREAD('ENTR',LU_INPUT,IOS)
       If (IOS == 1) Exit READ_ENTR_LOOP
       Read(LU_INPUT,ENTR,End=28,IOSTAT=IOS)

       If (Trim(KNOWN_DOOR_NAMES(51)) /= 'null') Then
          Write(MESSAGE,'(A,A,A)') &
               'ERROR: EVAC line ',Trim(HPT%ID_NAME), &
               ' problem with KNOWN_DOOR_NAMES'
          Call SHUTDOWN(MESSAGE)
       End If
       If (Trim(KNOWN_DOOR_NAMES(1)) == 'null') Then
          i = 0 ! no doors given
       Else
          i = 50 ! known door names given
          Do While ( Trim(KNOWN_DOOR_NAMES(i)) == 'null' .And. &
               i > 0)
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
       Do I=1,5,2
          If (XB(I) > XB(I+1)) Then
             DUMMY   = XB(I)
             XB(I)   = XB(I+1)
             XB(I+1) = DUMMY
          End If
       End Do
       !
       PNX%X1 = XB(1)
       PNX%X2 = XB(2)
       PNX%Y1 = XB(3)
       PNX%Y2 = XB(4)
       PNX%Z1 = XB(5)
       PNX%Z2 = XB(6)
       PNX%IOR = IOR
       ! 
       PNX%COLOR_INDEX = 0
       If (QUANTITY == 'BLACK')    PNX%COLOR_INDEX = 0
       If (QUANTITY == 'YELLOW')   PNX%COLOR_INDEX = 1
       If (QUANTITY == 'BLUE')     PNX%COLOR_INDEX = 2
       If (QUANTITY == 'RED')      PNX%COLOR_INDEX = 3
       If (QUANTITY == 'GREEN')    PNX%COLOR_INDEX = 4
       If (QUANTITY == 'MAGENTA')  PNX%COLOR_INDEX = 5
       If (QUANTITY == 'CYAN')     PNX%COLOR_INDEX = 6
       If (QUANTITY == 'WHITE')    PNX%COLOR_INDEX = 7
       ! 
       PNX%ID_NAME    = ID
       PNX%CLASS_NAME = PERS_ID

       PNX%IPC = 0
       Do ipc= 1, npc_pers
          pcp => evac_person_classes(ipc)
          If ( pcp%id_name == PERS_ID ) PNX%IPC = IPC
       End Do
       PNX%TO_NODE    = TO_NODE
       PNX%T_first    = 0.0_EB
       PNX%T_last     = 0.0_EB
       PNX%ICOUNT     = 0
       PNX%Flow       = MAX_FLOW

       Select Case (IOR)
       Case (-1,+1)
          If (WIDTH <= 0.0_EB) Then
             PNX%Width = XB(4) - XB(3)
          Else
             PNX%Width = WIDTH
          End If
       Case (-2,+2)
          If (WIDTH <= 0.0_EB) Then
             PNX%Width = XB(2) - XB(1)
          Else
             PNX%Width = WIDTH
          End If
       Case (3)
          If ( (XB(4)-XB(3)) <= 0.0_EB .Or. (XB(2)-XB(1)) <= 0.0_EB) Then
             Write(MESSAGE,'(A,I4,A)') &
                  'ERROR: ENTR',N,' IOR=3 but not 3-dim object'
             Call SHUTDOWN(MESSAGE)
          End If
       Case (0)
          If ( (XB(4)-XB(3)) <= 0.0_EB .Or. (XB(2)-XB(1)) <= 0.0_EB) Then
             Write(MESSAGE,'(A,I4,A)') &
                  'ERROR: ENTR',N,' no IOR but not 3-dim object'
             Call SHUTDOWN(MESSAGE)
          End If
       Case Default
          Write(MESSAGE,'(A,I4,A)') &
               'ERROR: ENTR',N,' problem with IOR'
          Call SHUTDOWN(MESSAGE)
       End Select
       ! 
       ! Check which evacuation floor
       ii = 0
       n_tmp = 0
       PNX_MeshLoop: Do i = 1, nmeshes
          If (evacuation_only(i) .And. evacuation_grid(i)) Then
             n_tmp = n_tmp + 1
             If ( (PNX%Z1 >= Meshes(i)%ZS .And. PNX%Z2 <= Meshes(i)%ZF).And. &
                  (PNX%Y1 >= Meshes(i)%YS .And. PNX%Y2 <= Meshes(i)%YF).And. &
                  (PNX%X1 >= Meshes(i)%XS .And. PNX%X2 <= Meshes(i)%XF)) Then
                If (Trim(MESH_ID) == 'null' .Or. &
                     Trim(MESH_ID) == Trim(MESH_NAME(i))) Then
                   ii = ii + 1
                   PNX%IMESH = i
                   PNX%TO_INODE = n_tmp
                   PNX%TO_NODE  = MESH_NAME(i)
                End If
             End If
          End If
       End Do PNX_MeshLoop
       If (PNX%IMESH == 0) Then
          Write(MESSAGE,'(A,A,A)') &
               'ERROR: ENTR line ',Trim(PNX%ID_NAME), &
               ' problem with IMESH, no mesh found'
          Call SHUTDOWN(MESSAGE)
       End If
       If (ii > 1) Then
          Write(MESSAGE,'(A,A,A)') &
               'ERROR: ENTR line ',Trim(PNX%ID_NAME), &
               ' not an unique mesh found '
          Call SHUTDOWN(MESSAGE)
       End If

       ! Use the main_evac_grid flow field if none is given
       If (Trim(FLOW_FIELD_ID) == 'null') Then
          PNX%GRID_NAME  = Trim(PNX%TO_NODE)
       Else
          PNX%GRID_NAME  = FLOW_FIELD_ID
       End If

       Do i = 1, PNX%N_VENT_FFIELDS
          PNX%P_VENT_FFIELDS(i) = KNOWN_DOOR_PROBS(i)
          PNX%I_VENT_FFIELDS(i) = 0
          PNX%I_DOOR_NODES(i) = 0
          Do j = 1, n_exits
             If ( Trim(EVAC_EXITS(j)%ID_NAME) == &
                  Trim(KNOWN_DOOR_NAMES(i)) ) Then
                PNX%I_VENT_FFIELDS(i) = EVAC_EXITS(j)%I_VENT_FFIELD
                PNX%I_DOOR_NODES(i)   = EVAC_EXITS(j)%INODE
             End If
          End Do
          Do j = 1, n_doors
             If ( Trim(EVAC_DOORS(j)%ID_NAME) == &
                  Trim(KNOWN_DOOR_NAMES(i)) ) Then
                PNX%I_VENT_FFIELDS(i) = EVAC_DOORS(j)%I_VENT_FFIELD
                PNX%I_DOOR_NODES(i)   = EVAC_DOORS(j)%INODE
             End If
          End Do
          If ( PNX%I_VENT_FFIELDS(i)*PNX%I_DOOR_NODES(i) == 0 ) Then
             Write(MESSAGE,'(A,A,A,A,A)') &
                  'ERROR: ENTR line ',Trim(PNX%ID_NAME), &
                  ' problem with door/exit names, ', &
                  Trim(KNOWN_DOOR_NAMES(i)),' not found'
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
          If ( evacuation_only(i) .And. Trim(PNX%GRID_NAME) == &
               Trim(MESH_NAME(i)) ) Then
             PNX%I_VENT_FFIELDS(0) = i
             Exit PNX_Mesh2Loop
          End If
       End Do PNX_Mesh2Loop
       If ( PNX%I_VENT_FFIELDS(0) == 0 ) Then
          Write(MESSAGE,'(A,A,A,A,A)') &
               'ERROR: ENTR line ',Trim(PNX%ID_NAME), &
               ' problem with flow field name, ', &
               Trim(PNX%GRID_NAME),' not found'
          Call SHUTDOWN(MESSAGE)
       End If
       !
       ! 
    End Do READ_ENTR_LOOP
28  Rewind(LU_INPUT)

    If (n_nodes > 0 ) Then
       n_tmp = n_egrids
       Do n = 1, n_entrys
          n_tmp = n_tmp + 1
          EVAC_Node_List(n_tmp)%ID_NAME    = EVAC_ENTRYS(n)%ID_NAME
          EVAC_Node_List(n_tmp)%Mesh_Index = EVAC_ENTRYS(n)%IMESH
       End Do
    End If
    ! Check that door/corr/entry/exit have unique names
    If (n_nodes > 0 ) Then
       Do n = 1, n_nodes - 1
          Do i = n + 1, n_nodes
             If (Trim(EVAC_Node_List(n)%ID_NAME) == Trim(EVAC_Node_List(i)%ID_NAME)) Then
                
                Write(MESSAGE,'(8A)') &
                     'ERROR: ', Trim(EVAC_Node_List(n)%Node_Type), ': ', &
                     Trim(EVAC_Node_List(n)%ID_NAME), ' has same ID as ', &
                     Trim(EVAC_Node_List(i)%Node_Type), ': ', &
                     Trim(EVAC_Node_List(i)%ID_NAME)
                Call SHUTDOWN(MESSAGE)
             End If
          End Do
       End Do
    End If

    !
    ! Read the EVAC lines
    ! 
    READ_EVAC_LOOP: Do N=1,NPC_EVAC
       HPT=>EVACUATION(N)
       !
       ID                       = 'null'
       QUANTITY                 = 'null'
       FLOW_FIELD_ID            = 'null'
       MESH_ID                = 'null'
       PERS_ID                  = 'null'
       SAMPLING_FACTOR          = 1      
       NUMBER_INITIAL_PERSONS   = 0
       XB                       = 0.0_EB
       ANGLE                    = -1000.0_EB
       EVACFILE                 = .False.
       AFTER_REACTION_TIME      = .False.
       T_START                  = -99.0_EB
       T_STOP                   = -99.0_EB
       GN_MIN                   = 1
       GN_MAX                   = 1      

       KNOWN_DOOR_NAMES         = 'null'
       KNOWN_DOOR_PROBS         = 1.0_EB
       !
       !
       Call CHECKREAD('EVAC',LU_INPUT,IOS)
       If (IOS == 1) Exit READ_EVAC_LOOP
       Read(LU_INPUT,EVAC,End=25,IOSTAT=IOS)

       If (Trim(KNOWN_DOOR_NAMES(51)) /= 'null') Then
          Write(MESSAGE,'(A,A,A)') &
               'ERROR: EVAC line ',Trim(HPT%ID_NAME), &
               ' problem with KNOWN_DOOR_NAMES'
          Call SHUTDOWN(MESSAGE)
       End If
       If (Trim(KNOWN_DOOR_NAMES(1)) == 'null') Then
          i = 0 ! no doors given
       Else
          i = 50 ! known door names given
          Do While ( Trim(KNOWN_DOOR_NAMES(i)) == 'null' .And. &
               i > 1)
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
       If (NUMBER_INITIAL_PERSONS > 0 .And. RESTART) &
            NUMBER_INITIAL_PERSONS =  0
       ! 
       If (NUMBER_INITIAL_PERSONS > 0) EVACFILE = .True.
       !
       !
       HPT%ID_NAME    = ID
       HPT%CLASS_NAME = PERS_ID
       HPT%T_START    = T_START

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
          If ( Trim(pcp%id_name) == Trim(PERS_ID) ) HPT%IPC = IPC
       End Do
       ! 
       HPT%COLOR_INDEX = 0
       If (QUANTITY == 'BLACK')    HPT%COLOR_INDEX = 0
       If (QUANTITY == 'YELLOW')   HPT%COLOR_INDEX = 1
       If (QUANTITY == 'BLUE')     HPT%COLOR_INDEX = 2
       If (QUANTITY == 'RED')      HPT%COLOR_INDEX = 3
       If (QUANTITY == 'GREEN')    HPT%COLOR_INDEX = 4
       If (QUANTITY == 'MAGENTA')  HPT%COLOR_INDEX = 5
       If (QUANTITY == 'CYAN')     HPT%COLOR_INDEX = 6
       If (QUANTITY == 'WHITE')    HPT%COLOR_INDEX = 7
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

       ! Check which evacuation floor
       ii = 0
       HP_MeshLoop: Do i = 1, nmeshes
          If (evacuation_only(i) .And. evacuation_grid(i)) Then
             If ( (HPT%Z1 >= Meshes(i)%ZS .And. HPT%Z2 <= Meshes(i)%ZF).And. &
                  (HPT%Y1 >= Meshes(i)%YS .And. HPT%Y2 <= Meshes(i)%YF).And. &
                  (HPT%X1 >= Meshes(i)%XS .And. HPT%X2 <= Meshes(i)%XF)) Then
                If (Trim(MESH_ID) == 'null' .Or. &
                     Trim(MESH_ID) == Trim(MESH_NAME(i))) Then
                   ii = ii + 1
                   HPT%IMESH = i
                End If
             End If
          End If
       End Do HP_MeshLoop
       If (HPT%IMESH == 0) Then
          Write(MESSAGE,'(A,A,A)') &
               'ERROR: EVAC line ',Trim(HPT%ID_NAME), &
               ' problem with IMESH, no mesh found'
          Call SHUTDOWN(MESSAGE)
       End If
       If (ii > 1) Then
          Write(MESSAGE,'(A,A,A)') &
               'ERROR: EVAC line ',Trim(HPT%ID_NAME), &
               ' not an unique mesh found '
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
          HPT%P_VENT_FFIELDS(i) = KNOWN_DOOR_PROBS(i)
          HPT%I_VENT_FFIELDS(i) = 0
          HPT%I_DOOR_NODES(i) = 0
          Do j = 1, n_exits
             If ( Trim(EVAC_EXITS(j)%ID_NAME) == &
                  Trim(KNOWN_DOOR_NAMES(i)) ) Then
                HPT%I_VENT_FFIELDS(i) = EVAC_EXITS(j)%I_VENT_FFIELD
                HPT%I_DOOR_NODES(i)   = EVAC_EXITS(j)%INODE
             End If
          End Do
          Do j = 1, n_doors
             If ( Trim(EVAC_DOORS(j)%ID_NAME) == &
                  Trim(KNOWN_DOOR_NAMES(i)) ) Then
                HPT%I_VENT_FFIELDS(i) = EVAC_DOORS(j)%I_VENT_FFIELD
                HPT%I_DOOR_NODES(i)   = EVAC_DOORS(j)%INODE
             End If
          End Do
          If ( HPT%I_VENT_FFIELDS(i)*HPT%I_DOOR_NODES(i) == 0 ) Then
             Write(MESSAGE,'(A,A,A,A,A)') &
                  'ERROR: EVAC line ',Trim(HPT%ID_NAME), &
                  ' problem with door/exit names, ', &
                  Trim(KNOWN_DOOR_NAMES(i)),' not found'
             Call SHUTDOWN(MESSAGE)
          End If
       End Do
       !
       ! No knwon doors given, use the flow_field_id value
       ! 
       HPT%P_VENT_FFIELDS(0) = 1.0_EB
       HPT%I_VENT_FFIELDS(0) = 0
       HPT%I_DOOR_NODES(0) = 0
       HP_Mesh2Loop: Do i = 1, nmeshes
          If ( evacuation_only(i) .And. Trim(HPT%GRID_NAME) == &
               Trim(MESH_NAME(i)) ) Then
             HPT%I_VENT_FFIELDS(0) = i
             Exit HP_Mesh2Loop
          End If
       End Do HP_Mesh2Loop
       If ( HPT%I_VENT_FFIELDS(0) == 0 ) Then
          Write(MESSAGE,'(A,A,A,A,A)') &
               'ERROR: EVAC line ',Trim(HPT%ID_NAME), &
               ' problem with flow field name, ', &
               Trim(HPT%GRID_NAME),' not found'
          Call SHUTDOWN(MESSAGE)
       End If
       !
    End Do READ_EVAC_LOOP
25  Rewind(LU_INPUT)
    !
    ! Read the EVHO lines
    !
    READ_EVHO_LOOP: Do N = 1, N_HOLES
       EHX=>EVAC_HOLES(N)
       !
       ID            = 'null'
       XB            = 0.0_EB
       EVAC_ID       = 'null'
       PERS_ID       = 'null'
       MESH_ID     = 'null'
       !
       Call CHECKREAD('EVHO',LU_INPUT,IOS)
       If (IOS == 1) Exit READ_EVHO_LOOP
       Read(LU_INPUT,EVHO,End=30,IOSTAT=IOS)
       !
       Do I=1,5,2
          If (XB(I) > XB(I+1)) Then
             DUMMY   = XB(I)
             XB(I)   = XB(I+1)
             XB(I+1) = DUMMY
          End If
       End Do
       !
       EHX%X1 = XB(1)
       EHX%X2 = XB(2)
       EHX%Y1 = XB(3)
       EHX%Y2 = XB(4)
       EHX%Z1 = XB(5)
       EHX%Z2 = XB(6)
       EHX%ID_NAME   = ID
       EHX%EVAC_ID   = EVAC_ID
       EHX%PERS_ID   = PERS_ID
       ! 
       ! Check which evacuation floor
       ii = 0
       EHX_MeshLoop: Do i = 1, nmeshes
          If (evacuation_only(i) .And. evacuation_grid(i)) Then
             If ( (EHX%Z1 >= Meshes(i)%ZS .And. EHX%Z2 <= Meshes(i)%ZF).And. &
                  (EHX%Y1 >= Meshes(i)%YS .And. EHX%Y2 <= Meshes(i)%YF).And. &
                  (EHX%X1 >= Meshes(i)%XS .And. EHX%X2 <= Meshes(i)%XF)) Then
                If (Trim(MESH_ID) == 'null' .Or. &
                     Trim(MESH_ID) == Trim(MESH_NAME(i))) Then
                   ii = ii + 1
                   EHX%IMESH = i
                   EHX%GRID_NAME  = MESH_NAME(i)
                End If
             End If
          End If
       End Do EHX_MeshLoop
       If (EHX%IMESH == 0) Then
          Write(MESSAGE,'(A,A,A)') &
               'ERROR: EVHO line ',Trim(EHX%ID_NAME), &
               ' problem with IMESH, no mesh found'
          Call SHUTDOWN(MESSAGE)
       End If
       If (ii > 1) Then
          Write(MESSAGE,'(A,A,A)') &
               'ERROR: EVHO line ',Trim(EHX%ID_NAME), &
               ' not an unique mesh found '
          Call SHUTDOWN(MESSAGE)
       End If
    End Do READ_EVHO_LOOP
30  Rewind(LU_INPUT)
    !
    ! Read the EVSS lines
    !
    READ_EVSS_LOOP: Do N = 1, N_SSTANDS
       ESS => EVAC_SSTANDS(N)
       !
       ID            = 'null'
       XB            = 0.0_EB
       MESH_ID     = 'null'
       IOR           = 0
       HEIGHT        = 0.0_EB
       HEIGHT0       = 0.0_EB
       FAC_V0_UP     = 0.0_EB
       FAC_V0_DOWN   = 0.0_EB
       FAC_V0_HORI   = 0.0_EB
       ESC_SPEED     = 0.0_EB
       !
       Call CHECKREAD('EVSS',LU_INPUT,IOS)
       If (IOS == 1) Exit READ_EVSS_LOOP
       Read(LU_INPUT,EVSS,End=31,IOSTAT=IOS)
       !
       Do I=1,5,2
          If (XB(I) > XB(I+1)) Then
             DUMMY   = XB(I)
             XB(I)   = XB(I+1)
             XB(I+1) = DUMMY
          End If
       End Do
       !
       ESS%X1 = XB(1)
       ESS%X2 = XB(2)
       ESS%Y1 = XB(3)
       ESS%Y2 = XB(4)
       ESS%Z1 = XB(5)
       ESS%Z2 = XB(6)
       ESS%ID_NAME     = ID
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
       ESS%IOR         = IOR
       ! 
       ! Check which evacuation floor
       ii = 0
       ESS_MeshLoop: Do i = 1, nmeshes
          If (evacuation_only(i) .And. evacuation_grid(i)) Then
             If ( (ESS%Z1 >= Meshes(i)%ZS .And. ESS%Z2 <= Meshes(i)%ZF).And. &
                  (ESS%Y1 >= Meshes(i)%YS .And. ESS%Y2 <= Meshes(i)%YF).And. &
                  (ESS%X1 >= Meshes(i)%XS .And. ESS%X2 <= Meshes(i)%XF)) Then
                If (Trim(MESH_ID) == 'null' .Or. &
                     Trim(MESH_ID) == Trim(MESH_NAME(i))) Then
                   ii = ii + 1
                   ESS%IMESH = i
                   ESS%GRID_NAME  = MESH_NAME(i)
                End If
             End If
          End If
       End Do ESS_MeshLoop
       If (ESS%IMESH == 0) Then
          Write(MESSAGE,'(A,A,A)') &
               'ERROR: EVSS line ',Trim(ESS%ID_NAME), &
               ' problem with IMESH, no mesh found'
          Call SHUTDOWN(MESSAGE)
       End If
       If (ii > 1) Then
          Write(MESSAGE,'(A,A,A)') &
               'ERROR: EVSS line ',Trim(ESS%ID_NAME), &
               ' not an unique mesh found '
          Call SHUTDOWN(MESSAGE)
       End If

       Select Case (IOR)
       Case(-1,+1)
          ESS%S = Sqrt((ESS%X2-ESS%X1)**2 + (ESS%H-ESS%H0)**2)
          ESS%COS_X = Abs(ESS%X2-ESS%X1)/ &
               Sqrt((ESS%X2-ESS%X1)**2 + (ESS%H-ESS%H0)**2)
          ESS%COS_Y = 1.0_EB
          ESS%SIN_X = Abs(ESS%H-ESS%H0)/ &
               Sqrt((ESS%X2-ESS%X1)**2 + (ESS%H-ESS%H0)**2)
          ESS%SIN_Y = 0.0_EB
       Case(-2,+2)
          ESS%S = Sqrt((ESS%Y2-ESS%Y1)**2 + (ESS%H-ESS%H0)**2)
          ESS%COS_X = 1.0_EB
          ESS%COS_Y = Abs(ESS%Y2-ESS%Y1)/ &
               Sqrt((ESS%Y2-ESS%Y1)**2 + (ESS%H-ESS%H0)**2)
          ESS%SIN_X = 0.0_EB
          ESS%SIN_Y = Abs(ESS%H-ESS%H0)/ &
               Sqrt((ESS%Y2-ESS%Y1)**2 + (ESS%H-ESS%H0)**2)
       Case Default
          Write(MESSAGE,'(A,I4,A)') &
               'ERROR: EVSS',N,' problem with IOR'
          Call SHUTDOWN(MESSAGE)
       End Select
       ! 
    End Do READ_EVSS_LOOP
31  Rewind(LU_INPUT)
    !
    ! Count the evacuation ramps and save their IDs.



    Close (LU_INPUT)    
    !
    ! Set the IMESH and IMESH2 for corridors
    Do n = 1, n_corrs
       Nodeloop2: Do i = 1, n_nodes
          If (EVAC_Node_List(i)%ID_NAME == EVAC_CORRS(n)%TO_NODE) Then
             EVAC_CORRS(n)%INODE2 = i
             If ( Trim(EVAC_Node_List(i)%Node_Type) == 'Exit' .Or. &
                  Trim(EVAC_Node_List(i)%Node_Type) == 'Door' .Or. &
                  Trim(EVAC_Node_List(i)%Node_Type) == 'Entry' ) Then
                EVAC_CORRS(n)%IMESH  = EVAC_Node_List(i)%Mesh_Index
                EVAC_CORRS(n)%IMESH2 = EVAC_Node_List(i)%Mesh_Index
             Else
                EVAC_CORRS(n)%IMESH  = 0
                EVAC_CORRS(n)%IMESH2 = n_egrids
             End If
             EVAC_Node_List(evac_corrs(n)%INODE)%Mesh_Index = &
                  EVAC_CORRS(n)%IMESH2
             Exit Nodeloop2
          End If
       End Do Nodeloop2
       If (EVAC_CORRS(n)%INODE2 == 0 .Or. &
            EVAC_CORRS(n)%IMESH2 == 0) Then
          Write(MESSAGE,'(A,I4,A)') &
               'ERROR: CORR',n,' problem with TO_NODE, loop2'
          Call SHUTDOWN(MESSAGE)
       End If
    End Do
    !
    Do n = 1, n_doors
       NodeLoop: Do i = 1, n_nodes
          If (EVAC_Node_List(i)%ID_NAME == EVAC_DOORS(n)%TO_NODE) Then
             EVAC_DOORS(n)%INODE2 = i
             EVAC_DOORS(n)%IMESH2 = EVAC_Node_List(i)%Mesh_Index
             Exit NodeLoop
          End If
       End Do NodeLoop
       If (EVAC_DOORS(n)%INODE2 == 0 .Or. &
            EVAC_DOORS(n)%IMESH2 == 0) Then
          Write(MESSAGE,'(A,I4,A)') &
               'ERROR: DOOR',n,' problem with TO_NODE'
          Call SHUTDOWN(MESSAGE)
       End If
    End Do


    Do n = 1, n_doors
       If (EVAC_DOORS(n)%KEEP_XY) Then
          i = EVAC_DOORS(n)%INODE2
          If (Trim(EVAC_Node_List(i)%Node_Type) == 'Door') Then
             PDX => EVAC_DOORS(EVAC_Node_List(i)%Node_Index)
             If ((EVAC_DOORS(n)%IOR /= -PDX%IOR) .Or. &
                  Abs(EVAC_DOORS(n)%Width-PDX%Width) > 0.1_EB ) Then
                Write(MESSAGE,'(A,I4,A)') &
                     'ERROR: DOOR',N,' KEEP_XY Problem'
                Call SHUTDOWN(MESSAGE)
             End If
          End If
          If (Trim(EVAC_Node_List(i)%Node_Type) == 'Entry') Then
             PNX => EVAC_ENTRYS(EVAC_Node_List(i)%Node_Index)
             If ((EVAC_DOORS(n)%IOR /= PNX%IOR) .Or. &
                  Abs(EVAC_DOORS(n)%Width-PNX%Width) > 0.1_EB ) Then
                Write(MESSAGE,'(A,I4,A)') &
                     'ERROR: DOOR',N,' KEEP_XY Problem'
                Call SHUTDOWN(MESSAGE)
             End If
          End If
       End If
    End Do

    If (COLOR_METHOD >= 0) EVAC_N_QUANTITIES = EVAC_N_QUANTITIES + 1
    If (OUTPUT_MOTIVE_FORCE) EVAC_N_QUANTITIES = EVAC_N_QUANTITIES + 1
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

       OUTPUT_QUANTITY(248)%NAME       = 'HUMAN_ANGULAR_ACCELERATION'
       OUTPUT_QUANTITY(248)%UNITS      = 'rad/s2'
       OUTPUT_QUANTITY(248)%SHORT_NAME = 'acc'
       
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

       If ( n-1 .ne. EVAC_N_QUANTITIES ) Then
          Write(MESSAGE,'(A,2I4,A)') &
               'ERROR: Evac output quantities ',EVAC_N_QUANTITIES,n-1, &
               ' Some bug in the program.'
          Call SHUTDOWN(MESSAGE)
       End If
 
    End If



  End Subroutine READ_EVAC

  Subroutine Initialize_Evac_Dumps
    Implicit None
    !
    Character(50) tcform
    Integer n_cols, i, j, nm
    Logical L_fed_read, L_fed_save, L_eff_read, L_eff_save, &
         L_status
    integer(4) n_egrids_tmp, ibar_tmp, jbar_tmp, kbar_tmp, &
         ntmp1, ntmp2, ntmp3, ntmp4, ntmp5, ntmp6, ios
    Real(FB) u_tmp, v_tmp
    !
    Type (MESH_TYPE), Pointer :: MFF
    !
    ! Logical unit numbers
    ! LU_EVACCSV: CHID_evac.csv, number of persons
    ! LU_EVACEFF: CHID_evac.eff, evacflow fields, binary
    ! LU_EVACFED: CHID_evac.fed, FED and soot, time dependent, binary
    !      Format: 1. row: n_egrids >=0  (Old Format, version 1.10)
    !              1a. row: n < 0 (New Format)
    !              1b. row: n_egrids,4,n_corrs=0,4 (New Format, version 1.11)
    !
    WRITE(EVAC_COMPILE_DATE,'(A)') evacrev(INDEX(evacrev,':')+1:LEN_TRIM(evacrev)-2)
    WRITE(LU_ERR,*) evacrev(INDEX(evacrev,':')+1:LEN_TRIM(evacrev)-2)
    READ (EVAC_COMPILE_DATE,'(I5)') EVAC_MODULE_REV
    WRITE(EVAC_COMPILE_DATE,'(A)') evacdate
    Call GET_REV_evac(EVAC_MODULE_REV,EVAC_COMPILE_DATE)
!    WRITE(EVAC_COMPILE_DATE,'(A)') EVAC_COMPILE_DATE(INDEX(EVAC_COMPILE_DATE,'(')+1:INDEX(EVAC_COMPILE_DATE,')')-1)
    !
    Write(LU_ERR,'(A)')          ' FDS+Evac Evacuation Module'
    Write(LU_OUTPUT,'(A)')          ' FDS+Evac Evacuation Module'
    Write(LU_ERR,'(A,A)')        ' FDS+Evac Compilation Date: ', &
         Trim(EVAC_COMPILE_DATE(INDEX(EVAC_COMPILE_DATE,'(')+1:INDEX(EVAC_COMPILE_DATE,')')-1))
    Write(LU_OUTPUT,'(A,A)')        ' FDS+Evac Compilation Date: ', &
         Trim(EVAC_COMPILE_DATE(INDEX(EVAC_COMPILE_DATE,'(')+1:INDEX(EVAC_COMPILE_DATE,')')-1))
    Write(LU_ERR,'(A,A)')  ' FDS+Evac Version         : ', &
         Trim(EVAC_VERSION)
    Write(LU_OUTPUT,'(A,A)')  ' FDS+Evac Version         : ', &
         Trim(EVAC_VERSION)
    Write(LU_ERR,'(A,i0/)')  ' FDS+Evac SVN Revision No.: ', &
         EVAC_MODULE_REV
    Write(LU_OUTPUT,'(A,i0/)')  ' FDS+Evac SVN Revision No.: ', &
         EVAC_MODULE_REV    

    Write(LU_ERR,fmt='(/a,i2)')  ' FDS+Evac Color_Method    :', &
         COLOR_METHOD
    If (Fed_Door_Crit >= 0) Then
       Write(LU_ERR,fmt='(a,f14.8)') &
                              ' FDS+Evac Fed_Door_Crit   :', FED_DOOR_CRIT
    Else
       ! Visibility S = 3/K, K is extinction coeff.
       Write(LU_ERR,fmt='(a,f14.8,a)') &
                              ' FDS+Evac Vis_Door_Crit   :', &
            Abs(FED_DOOR_CRIT), ' m'
       FED_DOOR_CRIT = 3.0_EB/FED_DOOR_CRIT ! Extinction coeff (1/m)
    End If
    If (NOT_RANDOM ) Write(LU_ERR,fmt='(a)') &
                              ' FDS+Evac Random seed is not used.'
    !
    L_fed_read = Btest(I_EVAC,3)
    L_fed_save = Btest(I_EVAC,1)
    L_eff_read = Btest(I_EVAC,2)
    L_eff_save = Btest(I_EVAC,0)

    n_cols = n_egrids + n_corrs + n_exits + n_doors + 1 + n_exits + n_doors
    ! Initialize the FED counters:
    icyc_old = -1
    n_dead = -1
    fed_max_alive = 0.0_EB
    fed_max = 0.0_EB
    !
    If (append) Then

       Open (LU_EVACCSV,file=FN_EVACCSV,form='formatted',status='old', &
            position='append')
       !
       If (L_fed_save) Then
          Open (LU_EVACFED,file=FN_EVACFED,form='unformatted', &
               status='old',position='append')
       End If
       ! 
       If (L_fed_read) Then
          Inquire (file=FN_EVACFED,exist=L_status)
          If (.Not. L_status) Then
             Write (LU_ERR,fmt='(a,a,a)') ' FDS+Evac No FED File: ', &
                  Trim(FN_EVACFED), ', FED and soot not used'
             Write (LU_OUTPUT,fmt='(a,a,a)') ' FDS+Evac No FED File: ', &
                  Trim(FN_EVACFED), ', FED and soot not used'
             l_fed_read = .False.
             l_fed_save = .False.
             I_EVAC = Ibclr(I_EVAC,3)  ! do not read FED
             I_EVAC = Ibclr(I_EVAC,1)  ! do not save FED
          Else
             Call SHUTDOWN('ERROR: Evac Dumps: FED, no restart yet')
             Open (LU_EVACFED,file=FN_EVACFED,form='unformatted', &
                  status='old')
             Read (LU_EVACFED,Iostat=ios) n_egrids_tmp
             If (ios.Ne.0) Then
                Write(MESSAGE,'(A)') 'ERROR: Init Evac Dumps: FED READ ERROR'
                Close (LU_EVACFED)
                Call SHUTDOWN(MESSAGE)
             End If
             If (n_egrids_tmp /= n_egrids) Then
                Write(MESSAGE,'(A,2I4)') &
                     'ERROR: Init Evac Dumps: FED ',n_egrids_tmp, n_egrids
                Call SHUTDOWN(MESSAGE)
             End If
          End If
       End If
       If (L_eff_read) Then
          Inquire (file=FN_EVACEFF,exist=L_status)
          If (L_status) Then
             Write (LU_ERR,fmt='(a,a,a/)') ' FDS+Evac EFF File: ', &
                  Trim(FN_EVACEFF), ' is used'
             Write (LU_OUTPUT,fmt='(a,a,a/)') ' FDS+Evac EFF File: ', &
                  Trim(FN_EVACEFF), ' is used'
             l_eff_save = .False.
             I_EVAC = Ibclr(I_EVAC,0)  ! do not save EFF
             Open (LU_EVACEFF,file=FN_EVACEFF,form='unformatted', &
                  status='old')
             Read (LU_EVACEFF,Iostat=ios) n_egrids_tmp
             If (ios.Ne.0) Then
                Write(MESSAGE,'(A)') 'ERROR: Init Evac Dumps: EFF READ ERROR'
                Close (LU_EVACEFF)
                Call SHUTDOWN(MESSAGE)
             End If
             If (n_egrids_tmp /= Count(EVACUATION_ONLY)) Then
                Write(MESSAGE,'(A,2I4)') &
                     'ERROR: Init Evac Dumps: EFF ',n_egrids_tmp, &
                     Count(EVACUATION_ONLY)
                Close (LU_EVACEFF)
                Call SHUTDOWN(MESSAGE)
             End If
          Else
             Write(MESSAGE,'(A,2I4)') &
                  'ERROR: Init Evac Dumps: EFF, no restart yet'
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
          Open (LU_EVACFED,file=FN_EVACFED,form='unformatted', &
               status='replace')
          ! First line: <0 new format
          !             -1: second line: #mesh #reals #corrs #reals #doors+exits #nreals
          !              (#reals: fed,soot,temp,radint,...)
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
          Write (LU_EVACFED) n_egrids_tmp, ntmp2, ntmp3, ntmp4, &
               ntmp5, ntmp6
          Write (LU_ERR,fmt='(a,a,a)') ' FDS+Evac FED File: ', &
               Trim(FN_EVACFED), ' is calculated and used'
          Write (LU_OUTPUT,fmt='(a,a,a)') ' FDS+Evac FED File: ', &
               Trim(FN_EVACFED), ' is calculated and used'
       End If
       !
       ! Number of evac flow fields is same as the number of all evac grids.
       If (L_eff_save) Then
          l_eff_read = .False.
          I_EVAC = Ibclr(I_EVAC,2)  ! do not read EFF
          Open (LU_EVACEFF,file=FN_EVACEFF,form='unformatted', &
               status='replace')
          n_egrids_tmp = Count(EVACUATION_ONLY)
          Write (LU_EVACEFF) n_egrids_tmp
          Write (LU_ERR,fmt='(a,a,a/)') ' FDS+Evac EFF File: ', &
               Trim(FN_EVACEFF), ' is calculated and used'
          Write (LU_OUTPUT,fmt='(a,a,a/)') ' FDS+Evac EFF File: ', &
               Trim(FN_EVACEFF), ' is calculated and used'
       End If
       ! 
       If (L_fed_read) Then
          Inquire (file=FN_EVACFED,exist=L_status)
          If (.Not. L_status) Then
             Write (LU_ERR,fmt='(a,a,a)') ' FDS+Evac No FED File: ', &
                  Trim(FN_EVACFED), ', FED and soot not used'
             Write (LU_OUTPUT,fmt='(a,a,a)') ' FDS+Evac No FED File: ', &
                  Trim(FN_EVACFED), ', FED and soot not used'
             l_fed_read = .False.
             l_fed_save = .False.
             I_EVAC = Ibclr(I_EVAC,3)  ! do not read FED
             I_EVAC = Ibclr(I_EVAC,1)  ! do not save FED
          Else
             Open (LU_EVACFED,file=FN_EVACFED,form='unformatted', &
                  status='old')
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
                Read (LU_EVACFED,Iostat=ios) n_egrids_tmp, ntmp2, ntmp3, ntmp4, &
                     ntmp5, ntmp6
                If (ios.Ne.0) Then
                   Write(MESSAGE,'(A)') 'ERROR: Init Evac Dumps: FED READ ERROR'
                   Close (LU_EVACFED)
                   Call SHUTDOWN(MESSAGE)
                End If
             End If

             ! Do not read old format. Do not read new format, if there the numbers
             ! are not: n_egrids, 4, n_corrs, 8 
             If ( ntmp2 /= 4 .Or. ntmp3 /= n_corrs .Or. ntmp1 >= 0 &
                  .Or. ntmp4 /= 8  .Or. &
                  ntmp5 /= n_doors+n_exits .Or. ntmp6 /= 4) Then
                Write (LU_ERR,fmt='(a,a,a)') ' FDS+Evac Error in FED File: ', &
                     Trim(FN_EVACFED), ', FED and soot not used'
                Write (LU_OUTPUT,fmt='(a,a,a)') ' FDS+Evac Error in FED File: ', &
                     Trim(FN_EVACFED), ', FED and soot not used'
                l_fed_read = .False.
                l_fed_save = .False.
                I_EVAC = Ibclr(I_EVAC,3) ! do not read FED
                I_EVAC = Ibclr(I_EVAC,1) ! do not save FED
                Close (LU_EVACFED)
             End If

             If ( l_fed_read .Or. l_fed_save ) Then
                Write (LU_ERR,fmt='(a,a,a)') ' FDS+Evac FED File: ', &
                     Trim(FN_EVACFED), ' is used'
                Write (LU_OUTPUT,fmt='(a,a,a)') ' FDS+Evac FED File: ', &
                     Trim(FN_EVACFED), ' is used'
             End If
             If (n_egrids_tmp /= n_egrids) Then
                Write(MESSAGE,'(A,2I4)') &
                     'ERROR: Init Evac Dumps: FED ',n_egrids_tmp, n_egrids
                Close (LU_EVACFED)
                Call SHUTDOWN(MESSAGE)
             End If
          End If
       End If
       ! 
       ! Number of evac flow fields is same as the number of all evac grids.
       If (L_eff_read) Then
          Inquire (file=FN_EVACEFF,exist=L_status)
          If (.Not. L_status) Then
             ! If evac flow fields are not found on the disk then 
             ! recalculate these.
             l_eff_save = .True.
             l_eff_read = .False.
             I_EVAC = Ibclr(I_EVAC,2) ! do not read EFF
             I_EVAC = Ibset(I_EVAC,0) ! save EFF
             Open (LU_EVACEFF,file=FN_EVACEFF,form='unformatted', &
                  status='replace')
             n_egrids_tmp = Count(EVACUATION_ONLY)
             Write (LU_EVACEFF) n_egrids_tmp
             Write (LU_ERR,fmt='(a,a,a/)') ' FDS+Evac EFF File: ', &
                  Trim(FN_EVACEFF), ' is (re)calculated'
             Write (LU_OUTPUT,fmt='(a,a,a/)') ' FDS+Evac EFF File: ', &
                  Trim(FN_EVACEFF), ' is (re)calculated'
          Else
             l_eff_save = .False.
             I_EVAC = Ibclr(I_EVAC,0) ! do not save EFF
             Open (LU_EVACEFF,file=FN_EVACEFF,form='unformatted', &
                  status='old')
             Read (LU_EVACEFF,Iostat=ios) n_egrids_tmp
             If (ios.Ne.0) Then
                Write(MESSAGE,'(A)') 'ERROR: Init Evac Dumps: EFF READ ERROR'
                Close (LU_EVACEFF)
                Call SHUTDOWN(MESSAGE)
             End If
             Write (LU_ERR,fmt='(a,a,a/)') ' FDS+Evac EFF File: ', &
                  Trim(FN_EVACEFF), ' is used'
             Write (LU_OUTPUT,fmt='(a,a,a/)') ' FDS+Evac EFF File: ', &
                  Trim(FN_EVACEFF), ' is used'
             If (n_egrids_tmp /= Count(EVACUATION_ONLY) ) Then
                Write(MESSAGE,'(A,2I4)') &
                     'ERROR: Init Evac Dumps: EFF ',n_egrids_tmp, &
                     Count(EVACUATION_ONLY)
                Close (LU_EVACEFF)
                Call SHUTDOWN(MESSAGE)
             End If
          End If
       End If

       !
       If ( l_fed_read .Or. l_fed_save ) Then
          ! Write the 'fed' columns
          n_dead = 0
          Open (LU_EVACCSV,file=FN_EVACCSV,form='formatted',status='replace')
          Write (LU_EVACCSV,*) n_cols+3
          Write (tcform,'(a,i4.4,a)') "(",n_cols+3,"(a,','),a)"
          Write (LU_EVACCSV,tcform) 'Time','Humans', &
               ('Floor', i=1,n_egrids), &
               ('Corridor', i=1,n_corrs), &
               ('Exit', i=1,n_exits), &
               ('Door', i=1,n_doors), &
               ('Exit', i=1,n_exits), &
               ('Door', i=1,n_doors), &
               'Fed','Fed','Fed'
          Write (LU_EVACCSV,tcform) 'Time','Inside', &
               ('Inside', i=1,n_egrids), &
               ('Inside', i=1,n_corrs), &
               ('Counter', i=1,n_exits), &
               ('Counter', i=1,n_doors), &
               ('Target', i=1,n_exits), &
               ('Target', i=1,n_doors), &
               'Counter','Fed','Fed'
          Write (LU_EVACCSV,tcform) 's','All Nodes', &
               (Trim(EVAC_Node_List(i)%GRID_NAME), i=1,n_egrids), &
               (Trim(EVAC_CORRS(i)%ID_NAME), i=1,n_corrs), &
               (Trim(EVAC_EXITS(i)%ID_NAME), i=1,n_exits), &
               (Trim(EVAC_DOORS(i)%ID_NAME), i=1,n_doors), &
               (Trim(EVAC_EXITS(i)%ID_NAME), i=1,n_exits), &
               (Trim(EVAC_DOORS(i)%ID_NAME), i=1,n_doors), &
               'Deads','FED_max','FED_max_alive'
       Else
          ! Do not write the 'fed' columns
          Open (LU_EVACCSV,file=FN_EVACCSV,form='formatted',status='replace')
          Write (LU_EVACCSV,*) n_cols
          Write (tcform,'(a,i4.4,a)') "(",n_cols,"(a,','),a)"
          Write (LU_EVACCSV,tcform) 'Time','Humans', &
               ('Floor', i=1,n_egrids), &
               ('Corridor', i=1,n_corrs), &
               ('Exit', i=1,n_exits), &
               ('Door', i=1,n_doors), &
               ('Exit', i=1,n_exits), &
               ('Door', i=1,n_doors)
          Write (LU_EVACCSV,tcform) 'Time','Inside', &
               ('Inside', i=1,n_egrids), &
               ('Inside', i=1,n_corrs), &
               ('Counter', i=1,n_exits), &
               ('Counter', i=1,n_doors), &
               ('Target', i=1,n_exits), &
               ('Target', i=1,n_doors)
          Write (LU_EVACCSV,tcform) 's','All Nodes', &
               (Trim(EVAC_Node_List(i)%GRID_NAME), i=1,n_egrids), &
               (Trim(EVAC_CORRS(i)%ID_NAME), i=1,n_corrs), &
               (Trim(EVAC_EXITS(i)%ID_NAME), i=1,n_exits), &
               (Trim(EVAC_DOORS(i)%ID_NAME), i=1,n_doors), &
               (Trim(EVAC_EXITS(i)%ID_NAME), i=1,n_exits), &
               (Trim(EVAC_DOORS(i)%ID_NAME), i=1,n_doors)
       End If

    End If                  ! append
    ! 
    ! Read the evac flow fields from the disk, if they exist.
    If ( L_eff_read ) Then
       Do nm = 1, NMESHES
          If (EVACUATION_ONLY(NM)) Then
             MFF=>MESHES(nm)
             Read (LU_EVACEFF,Iostat=ios) ibar_tmp, jbar_tmp, kbar_tmp
             If (ios.Ne.0) Then
                Write(MESSAGE,'(A)') 'ERROR: Init Evac Dumps: EFF READ ERROR'
                Close (LU_EVACEFF)
                Call SHUTDOWN(MESSAGE)
             End If
             If ( MFF%IBAR /= ibar_tmp .Or. MFF%JBAR /= jbar_tmp .Or. &
                  MFF%KBAR /= kbar_tmp ) Then
                Close (LU_EVACEFF)
                Call SHUTDOWN('ERROR: Problems to read the EFF file')
             End If
             Do  i = 0, MFF%IBAR+1
                Do j= 0, MFF%JBAR+1
                   Read (LU_EVACEFF,Iostat=ios) u_tmp, v_tmp
                   If (ios.Ne.0) Then
                      Write(MESSAGE,'(A)') 'ERROR: Init Evac Dumps: EFF READ ERROR'
                      Close (LU_EVACEFF)
                      Call SHUTDOWN(MESSAGE)
                   End If
                   MFF%U(i,j,:) = u_tmp
                   MFF%V(i,j,:) = v_tmp
                   MFF%W(i,j,:) = 0.0_EB
                End Do
             End Do
          End If
       End Do
    End If

    EVAC_Z_MIN =  Huge(EVAC_Z_MIN)
    EVAC_Z_MAX = -Huge(EVAC_Z_MIN)
    Do nm = 1, NMESHES
       MFF=>MESHES(nm)
       EVAC_Z_MIN = Min(EVAC_Z_MIN,Real(MFF%ZS,FB))
       EVAC_Z_MAX = Max(EVAC_Z_MAX,Real(MFF%ZF,FB))
    End Do

  End Subroutine Initialize_Evac_Dumps
      
!
  Subroutine INITIALIZE_EVACUATION(NM,ISTOP)
    Implicit None
    !
    ! Insert humans into the domain at the start of calculation
    !
    Integer, Intent(IN) :: NM
    Integer, Intent(OUT) :: ISTOP
    Real(EB) RN, RN1, simoDX, simoDY, TNOW
    Real(EB) VOL1, VOL2, X1, X2, Y1, Y2, Z1, Z2, &
         dist, d_max, G_mean, G_sd, G_high, G_low, x1_old, y1_old
    integer i,j,ii,jj,kk,ipc, izero, n_tmp, ie, nom
    integer i11, i22, group_size
    logical pp_see_group, is_solid
    integer iie, jje, iio, jjo, jjj, tim_ic, iii, i44
    real(eb) y_o, x_o, delta_y, delta_x, x_now, y_now, &
         xi, yj, x11, y11, group_x_sum, group_y_sum, &
         group_x_center, group_y_center
    integer :: i_endless_loop
    real(eb), dimension(6) :: y_tmp, x_tmp, r_tmp
    ! 
    type (mesh_type), pointer :: m
    !
    tnow=second()

    If ( .Not.(EVACUATION_ONLY(NM) .And. EVACUATION_GRID(NM)) ) Return
!!$    EVAC_MESH_ONLY: If ( EVACUATION_ONLY(NM) .And. &
!!$         EVACUATION_GRID(NM) ) Then
    !
    GTrunFlag = 0
    GaussFlag = 0
    !
    MESHES(NM)%N_HUMANS = 0
    MESHES(NM)%N_HUMANS_DIM = 10000
    Allocate(MESHES(NM)%HUMAN(MESHES(NM)%N_HUMANS_DIM),STAT=IZERO)
    Call ChkMemErr('INIT_EVACUATION','HUMAN',IZERO)
    !
    Allocate(MESHES(NM)%HUMAN_GRID(MESHES(NM)%IBAR,MESHES(NM)%JBAR), &
         STAT=IZERO)
    Call ChkMemErr('INIT_EVACUATION','HUMAN_GRID',IZERO)
    !
    Call POINT_TO_MESH(NM)
    !       Call POINT_TO_EVAC_MESH(NM)
    !
    ! Initialise HUMAN_GRID
    !
    Do i = 1,IBAR
       FED_INNER_LOOP: Do j= 1,JBAR
          x1 = XS + (i-1)*DXI + 0.5_EB*DXI
          y1 = YS + (j-1)*DETA + 0.5_EB*DETA
          z1 = 0.5_EB*(ZF+ZS)
          HUMAN_GRID(i,j)%X = x1
          HUMAN_GRID(i,j)%Y = y1
          HUMAN_GRID(i,j)%Z = z1
          HUMAN_GRID(i,j)%SOOT_DENS  = 0.0_EB
          HUMAN_GRID(i,j)%FED_CO_CO2_O2  = 0.0_EB
          HUMAN_GRID(i,j)%IMESH = 0
          HUMAN_GRID(i,j)%II = i
          HUMAN_GRID(i,j)%JJ = j
          HUMAN_GRID(i,j)%KK = 1

          ! If there are no fire grids, skip this intialization part
          If (.Not. Btest(I_EVAC,4) ) Cycle FED_INNER_LOOP

          If (.Not. SOLID(CELL_INDEX(i,j,1)) ) Then
             MESH_LOOP: Do NOM=1,NMESHES
                M=>MESHES(NOM)
                If (.Not. EVACUATION_ONLY(NOM)) Then
                   If (X1 >= M%XS .And. X1 <= M%XF .And. &
                        Y1 >= M%YS .And. Y1 <= M%YF .And. &
                        Z1 >= M%ZS .And. Z1 <= M%ZF) Then
                      II = Floor( M%CELLSI(Floor((X1-M%XS)*M%RDXINT)) &
                           + 1.0_EB  )
                      JJ = Floor( M%CELLSJ(Floor((Y1-M%YS)*M%RDYINT)) &
                           + 1.0_EB  )
                      KK = Floor( M%CELLSK(Floor((Z1-M%ZS)*M%RDZINT)) &
                           + 1.0_EB  )
                      If ( M%SOLID(M%CELL_INDEX(II,JJ,KK)) ) Then
                         HUMAN_GRID(i,j)%IMESH = 0
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
          Else
             ! This grid cell is solid ==> No humans in this cell
             HUMAN_GRID(i,j)%IMESH = 0
          End If
       End Do FED_INNER_LOOP
    End Do

    EVAC_CLASS_LOOP: Do IPC=1,NPC_EVAC
       !
       HPT=>EVACUATION(IPC)
       !
       ! If there is an initial number of humans, initialize
       !
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
       ! Check which evacuation floor (first mesh found is used.)
       n_tmp = 1
       HP_MeshLoop: Do i = 1, nmeshes
          If (evacuation_only(i) .And. evacuation_grid(i)) Then
             If ( (X1 >= Meshes(i)%XS .And. X2 <= Meshes(i)%XF) .And. &
                  (Y1 >= Meshes(i)%YS .And. Y2 <= Meshes(i)%YF) .And. &
                  (Z1 >= Meshes(i)%ZS .And. Z2 <= Meshes(i)%ZF) ) Then
                Exit HP_MeshLoop
             End If
             n_tmp = n_tmp + 1
          End If
       End Do HP_MeshLoop
       If (n_tmp > n_egrids) Then
          Write(MESSAGE,'(A,A,A)') &
               'ERROR: INIT_EVAC ',Trim(HPT%ID_NAME), &
               ' problem with IMESH, no mesh found'
          Call SHUTDOWN(MESSAGE)
       End If
       !
       PCP => EVAC_PERSON_CLASSES(HPT%IPC)
       !
       !
       i11 = 0

       INITIALIZATION_LOOP: Do I=1,HPT%N_INITIAL
          !
          Call Random_number(RN)
          group_size = HPT%GN_MIN - 1 +  &
               Int((HPT%GN_MAX-HPT%GN_MIN+1)*RN+0.5_EB)
          group_size = Max(Min(group_size, HPT%GN_MAX),HPT%GN_MIN)

          If ( i11+group_size > HPT%N_INITIAL ) Then
             group_size = HPT%N_INITIAL - i11
          End If

          If ( i11 >= HPT%N_INITIAL ) Exit INITIALIZATION_LOOP

          i22 = 0
          i_endless_loop = 0
          GROUP_SIZE_LOOP: Do 
             i22 = i22 + 1
             If (i22 > group_size) Exit GROUP_SIZE_LOOP
             group_X_sum = 0
             group_Y_sum = 0
             i11 = i11 + 1

             If (i22 == 1) Then
                If (group_size > 1) i33 = i33 + 1
                If (group_size == 1) ilh = ilh + 1
                Call Random_number(RN1)
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
             Call CLASS_PROPERTIES
             HR%IPC = HPT%IPC
             HR%IEL = IPC
             HR%I_Target = 0

             !
             BLK_LOOP:  Do
                If (i22 == 1) Then
                   Call Random_number(RN)
                   HR%X = X1 + RN*(X2-X1)
                   x1_old = HR%X
                   Call Random_number(RN)
                   HR%Y = Y1 + RN*(Y2-Y1)
                   y1_old = HR%Y
                Else
                   G_mean = 0.0_EB
                   G_sd   = 4.0_EB  ! units (m) std.dev.
                   G_high = 6.0_EB  ! units (m) cut-off
                   G_low  = -6.0_EB ! units (m) cut-off
                   ! First the angle, then the radial distance
                   Call Random_Number(rn)
                   simoDX = Sin(2.0_EB*Pi*rn)
                   simoDY = Cos(2.0_EB*Pi*rn)
                   G_mean = (2.0_EB/3.0_EB)*Sqrt(group_size/(Pi*GROUP_DENS))
                   G_sd   = G_mean  ! units (m) std.dev.
                   G_high = Max(3.0_EB,3.0_EB*GROUP_DENS)* G_mean ! units (m) cut-off
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

                If (i_endless_loop >= (8.0_EB*Max(1.0_EB,Log10(2.5_EB*VOL2))) / &
                     Max(1.0_EB,Log10((2.5_EB*VOL2)/(2.5_EB*VOL2-1))) ) Then
                   Write(LU_ERR,'(A,I4,A,I4,A,I6)') &
                        'ERROR: Initialize_Humans, EVAC line ', &
                        IPC, ', Mesh ', NM, ', i_human ', n_humans
                   Write(LU_OUTPUT,'(A,I4,A,I4,A,I6)') &
                        'ERROR: Initialize_Humans, EVAC line ', &
                        IPC, ', Mesh ', NM, ', i_human ', n_humans
                   ISTOP = 3
                   HR%SHOW = .True.    
                   HR%COLOR_INDEX = 6  ! Cyan
                   Exit INITIALIZATION_LOOP
                End If


                !Check, that a person is not put on top of some other person
                If (DENS_INIT > 2.0_EB) Then
                   ! High density is wanted
                   d_max = 0.0_EB
                Else
                   d_max = 1.0_EB*HR%B
                End If

                Is_Solid = .False.
                KK = Floor( CELLSK(Floor((HR%Z-ZS)*RDZINT)) + 1.0_EB )

                r_tmp(1) = HR%r_shoulder ! right circle
                r_tmp(2) = HR%r_torso     ! center circle
                r_tmp(3) = HR%r_shoulder ! left circle
                y_tmp(1) = HR%Y - Cos(HR%angle)*(HR%d_shoulder+HR%r_shoulder) ! right
                x_tmp(1) = HR%X + Sin(HR%angle)*(HR%d_shoulder+HR%r_shoulder)
                y_tmp(3) = HR%Y + Cos(HR%angle)*(HR%d_shoulder+HR%r_shoulder)
                x_tmp(3) = HR%X - Sin(HR%angle)*(HR%d_shoulder+HR%r_shoulder)
                II = Floor( CELLSI(Floor((x_tmp(1)-XS)*RDXINT)) + 1.0_EB )
                JJ = Floor( CELLSJ(Floor((y_tmp(1)-YS)*RDYINT)) + 1.0_EB )
                Is_Solid = (Is_Solid .Or. SOLID(CELL_INDEX(II,JJ,KK)))
                II = Floor( CELLSI(Floor((x_tmp(3)-XS)*RDXINT)) + 1.0_EB )
                JJ = Floor( CELLSJ(Floor((y_tmp(3)-YS)*RDYINT)) + 1.0_EB )
                Is_Solid = (Is_Solid .Or. SOLID(CELL_INDEX(II,JJ,KK)))
                y_tmp(2) = HR%Y + Sin(HR%angle)*(HR%r_torso)  ! torso
                x_tmp(2) = HR%X + Cos(HR%angle)*(HR%r_torso)
                II = Floor( CELLSI(Floor((x_tmp(2)-XS)*RDXINT)) + 1.0_EB )
                JJ = Floor( CELLSJ(Floor((y_tmp(2)-YS)*RDYINT)) + 1.0_EB )
                Is_Solid = (Is_Solid .Or. SOLID(CELL_INDEX(II,JJ,KK)))
                y_tmp(2) = HR%Y - Sin(HR%angle)*(HR%r_torso)  ! torso
                x_tmp(2) = HR%X - Cos(HR%angle)*(HR%r_torso)
                II = Floor( CELLSI(Floor((x_tmp(2)-XS)*RDXINT)) + 1.0_EB )
                JJ = Floor( CELLSJ(Floor((y_tmp(2)-YS)*RDYINT)) + 1.0_EB )
                Is_Solid = (Is_Solid .Or. SOLID(CELL_INDEX(II,JJ,KK)))

                r_tmp(1) = HR%r_shoulder ! right circle
                r_tmp(2) = HR%r_torso     ! center circle
                r_tmp(3) = HR%r_shoulder ! left circle
                y_tmp(1) = HR%Y - Cos(HR%angle)*HR%d_shoulder ! right
                x_tmp(1) = HR%X + Sin(HR%angle)*HR%d_shoulder
                y_tmp(2) = HR%Y ! torso
                x_tmp(2) = HR%X
                y_tmp(3) = HR%Y + Cos(HR%angle)*HR%d_shoulder ! left
                x_tmp(3) = HR%X - Sin(HR%angle)*HR%d_shoulder
                II = Floor( CELLSI(Floor((x_tmp(1)-XS)*RDXINT)) + 1.0_EB )
                JJ = Floor( CELLSJ(Floor((y_tmp(1)-YS)*RDYINT)) + 1.0_EB )
                Is_Solid = (Is_Solid .Or. SOLID(CELL_INDEX(II,JJ,KK)))
                II = Floor( CELLSI(Floor((x_tmp(3)-XS)*RDXINT)) + 1.0_EB )
                JJ = Floor( CELLSJ(Floor((y_tmp(3)-YS)*RDYINT)) + 1.0_EB )
                Is_Solid = (Is_Solid .Or. SOLID(CELL_INDEX(II,JJ,KK)))
                II = Floor( CELLSI(Floor((x_tmp(2)-XS)*RDXINT)) + 1.0_EB )
                JJ = Floor( CELLSJ(Floor((y_tmp(2)-YS)*RDYINT)) + 1.0_EB )
                Is_Solid = (Is_Solid .Or. SOLID(CELL_INDEX(II,JJ,KK)))

                If (.Not.Is_Solid) Then
                   P2PLoop: Do ie = 1, n_humans - 1
                      HRE => HUMAN(ie)
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
                            DIST = Sqrt((x_tmp(jjj)-x_tmp(iii))**2 + &
                                 (y_tmp(jjj)-y_tmp(iii))**2) - &
                                 (r_tmp(jjj)+r_tmp(iii))
                            If ( DIST < d_max) Then
                               i_endless_loop = i_endless_loop + 1
                               Cycle BLK_LOOP
                            End If
                         End Do
                      End Do
                   End Do P2PLoop

                   ! Put here EVHO (evac hole) checks, n_holes
                   EH_Loop: Do ie = 1, n_holes
                      EHX => EVAC_HOLES(ie)

                      If ( Trim(EHX%EVAC_ID) /= 'null' .And. &
                           Trim(EHX%EVAC_ID) /= Trim(HPT%ID_NAME)) Then
                         Cycle EH_Loop
                      End If
                      If ( Trim(EHX%PERS_ID) /= 'null' .And. &
                           Trim(EHX%PERS_ID) /= Trim(HPT%CLASS_NAME)) Then
                         Cycle EH_Loop
                      End If

                      If ( (EHX%IMESH == NM) .And. &
                           (EHX%Y1 <= HR%Y .And. EHX%Y2 >= HR%Y) .And. &
                           (EHX%X1 <= HR%X .And. EHX%X2 >= HR%X) ) Then
                         i_endless_loop = i_endless_loop + 1
                         Cycle BLK_LOOP
                      End If
                   End Do EH_Loop

                   If (i22 > 1) Then
                      ! Chech if the new member will see the first member of the group.
                      X11 = HR%X
                      Y11 = HR%Y

                      XI  = CELLSI(Floor((x1_old-XS)*RDXINT))
                      YJ  = CELLSJ(Floor((y1_old-YS)*RDYINT))
                      IIE = Floor(XI+1.0_EB)
                      JJE = Floor(YJ+1.0_EB)
                      PP_see_group = .True.  ! oletusarvo
                      If (Abs(X1_OLD-X11) >= Abs(Y1_OLD-Y11)) Then
                         If ( iie < ii) Then
                            iio = iie
                            jjo = jje
                            iie = ii
                            jje = jj
                            y_o = Y1_OLD
                            x_o = X1_OLD
                            Delta_y = (Y11 - Y1_OLD)
                         Else
                            Delta_y = (Y1_OLD - Y11)
                            iio = ii
                            jjo = jj
                            y_o = Y11
                            x_o = X11
                         End If
                         Delta_x = Abs(X1_OLD - X11)
                         x_now = 0.0_EB
                         PP_see_loop_x: Do iii = iio+1, iie-1
                            x_now = x_now + DX(iii)
                            y_now = y_o + x_now*(Delta_y/Delta_x)
                            jjj = Floor(CELLSJ(Floor((y_now-YS)*RDYINT))+1.0_EB)
                            tim_ic = CELL_INDEX(iii,jjj,KK)
                            If (SOLID(tim_ic)) Then
                               PP_see_group = .False.
                               Exit PP_see_loop_x
                            End If
                         End Do PP_see_loop_x
                      Else 
                         If ( jje < jj) Then
                            iio = iie
                            jjo = jje
                            iie = ii
                            jje = jj
                            y_o = Y1_OLD
                            x_o = X1_OLD
                            Delta_x = (X11 - X1_OLD)
                         Else
                            Delta_x = (X1_OLD - X11)
                            iio = ii
                            jjo = jj
                            y_o = Y11
                            x_o = X11
                         End If
                         Delta_y = Abs(Y1_OLD - Y11)
                         y_now = 0.0_EB
                         PP_see_loop_y: Do jjj = jjo+1, jje-1
                            y_now = y_now + DY(jjj)
                            x_now = x_o + y_now*(Delta_x/Delta_y)
                            iii = Floor(CELLSI(Floor((x_now-XS)*RDXINT))+1.0_EB)
                            tim_ic = CELL_INDEX(iii,jjj,KK)
                            If (SOLID(tim_ic)) Then
                               PP_see_group = .False.
                               Exit PP_see_loop_y
                            End If
                         End Do PP_see_loop_y
                      End If
                   Else
                      PP_see_group = .True.
                   End If

                   If ( .Not. PP_see_group ) Then 
                      i_endless_loop = i_endless_loop + 1
                      Cycle BLK_LOOP
                   End If

                   ! Coordinates are OK, exit the loop
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
             HR%Commitment = RN1
             HR%SHOW = .True.    


             Select Case (COLOR_METHOD)
             Case (-1)
                HR%COLOR_INDEX = 0
             Case (0)
                HR%COLOR_INDEX = HPT%COLOR_INDEX 
             Case (1)
                HR%COLOR_INDEX = Mod(group_size-1,5)
             Case (2)
                If (HR%GROUP_ID > 0 ) Then
                   HR%COLOR_INDEX = Mod(HR%GROUP_ID,5) + 1
                Else
                   HR%COLOR_INDEX = 0 ! lonely human
                End If
             Case (4)
                HR%COLOR_INDEX = HPT%COLOR_INDEX ! default
             Case (5)
                HR%COLOR_INDEX = 0
             Case Default
                Write(MESSAGE,'(A,I3,A)') &
                     'ERROR: READ_EVAC COLOR METHOD',COLOR_METHOD, &
                     ' is not defined'
                Call SHUTDOWN(MESSAGE)
             End Select

             HR%IMESH       = HPT%IMESH
             HR%INODE       = n_tmp
             HR%NODE_NAME   = Trim(MESH_NAME(n_tmp))
             HR%IOR       = 0  
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



             group_X_sum = group_X_sum + HR%X
             group_Y_sum = group_Y_sum + HR%Y
             If ( group_size == 1 ) Exit  GROUP_SIZE_LOOP
             If ( i22 == group_size .And. group_size > 1 ) Then
                group_X_center = group_X_sum / Max(1,group_size)
                group_Y_center = group_Y_sum / Max(1,group_size)
                !
                II = Floor( CELLSI(Floor((group_X_center-XS)* &
                     RDXINT)) + 1.0_EB )
                JJ = Floor( CELLSJ(Floor((group_Y_center-YS)* &
                     RDYINT)) + 1.0_EB )
                KK = Floor( CELLSK(Floor((HR%Z-ZS)*RDZINT)) + 1.0_EB )
                x1_old = group_X_center
                y1_old = group_Y_center

                PP_see_group = .True.
                Do i44 = 1,group_size
                   HR => HUMAN( N_HUMANS - group_size + i44 )
                   X11 = HR%X
                   Y11 = HR%Y
                   XI  = CELLSI(Floor((x1_old-XS)*RDXINT))
                   YJ  = CELLSJ(Floor((y1_old-YS)*RDYINT))
                   IIE = Floor(XI+1.0_EB)
                   JJE = Floor(YJ+1.0_EB)
                   If (Abs(X1_OLD-X11) >= Abs(Y1_OLD-Y11)) Then
                      If ( iie < ii) Then
                         iio = iie
                         jjo = jje
                         iie = ii
                         jje = jj
                         y_o = Y1_OLD
                         x_o = X1_OLD
                         Delta_y = (Y11 - Y1_OLD)
                      Else
                         Delta_y = (Y1_OLD - Y11)
                         iio = ii
                         jjo = jj
                         y_o = Y11
                         x_o = X11
                      End If
                      Delta_x = Abs(X1_OLD - X11)
                      x_now = 0.0_EB
                      PP_see_loop_x2: Do iii = iio+1, iie-1
                         x_now = x_now + DX(iii)
                         y_now = y_o + x_now*(Delta_y/Delta_x)
                         jjj = Floor(CELLSJ(Floor((y_now-YS)*RDYINT))+1.0_EB)
                         tim_ic = CELL_INDEX(iii,jjj,KK)
                         If (SOLID(tim_ic)) Then
                            PP_see_group = .False.
                            Exit PP_see_loop_x2
                         End If
                      End Do PP_see_loop_x2
                   Else 
                      If ( jje < jj) Then
                         iio = iie
                         jjo = jje
                         iie = ii
                         jje = jj
                         y_o = Y1_OLD
                         x_o = X1_OLD
                         Delta_x = (X11 - X1_OLD)
                      Else
                         Delta_x = (X1_OLD - X11)
                         iio = ii
                         jjo = jj
                         y_o = Y11
                         x_o = X11
                      End If
                      Delta_y = Abs(Y1_OLD - Y11)
                      y_now = 0.0_EB
                      PP_see_loop_y2: Do jjj = jjo+1, jje-1
                         y_now = y_now + DY(jjj)
                         x_now = x_o + y_now*(Delta_x/Delta_y)
                         iii = Floor(CELLSI(Floor((x_now-XS)*RDXINT))+1.0_EB)
                         tim_ic = CELL_INDEX(iii,jjj,KK)
                         If (SOLID(tim_ic)) Then
                            PP_see_group = .False.
                            Exit PP_see_loop_y2
                         End If
                      End Do PP_see_loop_y2
                   End If

                   If ( .Not. PP_see_group ) Then 
                      i_endless_loop = i_endless_loop + 1
                      i22 = 0  ! Start at the beginning of the group
                      i33 = i33 - 1  ! group index
                      ilh = ilh - 1  ! lonely human index
                      i11 = i11 - group_size ! human index (evac-line)
                      N_HUMANS = N_HUMANS - group_size ! total # of humans
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

    Write (LU_ERR,fmt='(a,f8.2,a,i0,a,i0/)') ' EVAC: Time ', &
         0.0_EB,' mesh ',nm,' number of humans ',n_humans
    Write (LU_OUTPUT,fmt='(a,f8.2,a,i0,a,i0/)') ' EVAC: Time ', &
         0.0_EB,' mesh ',nm,' number of humans ',n_humans
    !
!!$    End If EVAC_MESH_ONLY
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
    Real(EB) RN
    Real(EB) x1_old, y1_old
    Integer I,J,II,JJ,KK, IZERO, ie, nom, j1
    Logical PP_see_door
    Integer iie, jje, iio, jjo, jjj, tim_ic, iii
    Real(EB) y_o, x_o, Delta_y, Delta_x, x_now, y_now, &
         xi, yj, x11, y11
    Character(26) :: GROUP_FFIELD_NAME
    Integer :: GROUP_FFIELD_I
    Real(EB), Dimension(:), Allocatable :: FED_max_Door
    Logical, Dimension(:), Allocatable :: Is_Known_Door, Is_Visible_Door
    Integer, Dimension(:), Allocatable :: Color_Tmp
    Integer :: i_tmp, i_egrid, color_index
    Real(EB) :: L2_min, max_fed, L2_tmp
    ! 
    Type (MESH_TYPE), Pointer :: M
    !
    If (.NOT.Any(EVACUATION_GRID)) RETURN

    !! TNOW=SECOND()
    !
!!!If (Any(EVACUATION_GRID)) Then
    !
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

    Write (LU_ERR,fmt='(/a)') ' EVAC: Initial positions of the humans'
    Write (LU_ERR,fmt='(a,a)') &
         ' person     x       y       z    Tpre    Tdet  ', &
         ' dia    v0   tau   i_gr i_ff'

    Allocate(Group_Known_Doors(1:i33_dim),STAT=IZERO)
    Call ChkMemErr('Initialize_Evacuation', &
         'Group_Known_Doors',IZERO) 
    Allocate(Human_Known_Doors(1:ilh_dim),STAT=IZERO)
    Call ChkMemErr('Initialize_Evacuation', &
         'Human_Known_Doors',IZERO) 

    !
    Allocate(Is_Known_Door(Max(1,n_doors+n_exits)),STAT=IZERO)
    Call ChkMemErr('Initialize_Evacuation','Is_Known_Door',IZERO) 
    Allocate(Is_Visible_Door(Max(1,n_doors+n_exits)),STAT=IZERO)
    Call ChkMemErr('Initialize_Evacuation','Is_Visible_Door',IZERO)
    Allocate(FED_max_Door(Max(1,n_doors+n_exits)),STAT=IZERO)
    Call ChkMemErr('Initialize_Evacuation','FED_max_Door',IZERO) 
    Allocate(Color_Tmp(Max(1,i33_dim)),STAT=IZERO)
    Call ChkMemErr('Initialize_Evacuation','Color_Tmp',IZERO) 

    i_egrid = 0
    Do nom = 1, NMESHES
       If ( .Not.(EVACUATION_ONLY(NoM) .And. EVACUATION_GRID(NoM)) ) Cycle
!!$       If ( EVACUATION_ONLY(NoM) .And. EVACUATION_GRID(NoM) ) Then
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
       End Do
       Group_List(1:)%GROUP_X = Group_List(1:)%GROUP_X / &
            Max(1,Group_List(1:)%GROUP_SIZE)
       Group_List(1:)%GROUP_Y = Group_List(1:)%GROUP_Y / &
            Max(1,Group_List(1:)%GROUP_SIZE)
       Group_List(1:)%Speed   = Group_List(1:)%Speed / &
            Max(1,Group_List(1:)%GROUP_SIZE)

       Do i = 1, M%N_HUMANS
          HR => M%HUMAN(I)
          ! group_id > 0: +group_id
          ! group_id < 0: -human_id (lonely humans)
          j  =  Max(0,HR%GROUP_ID)
          j1 = -Min(0,HR%GROUP_ID)
          ! Test, if this group has already a ffield (on this floor)
          ! Lonely humans have j=0 and group_i_ffields is always 0.
          If (Group_List(j)%GROUP_I_FFIELDS(i_egrid) == 0) Then
             If ( j == 0 ) Then
                x1_old = HR%X
                y1_old = HR%Y
             Else
                x1_old = Group_List(j)%GROUP_X
                y1_old = Group_List(j)%GROUP_Y
             End If
             HPT => EVACUATION(HR%IEL)
             FED_max_Door(:)    = 0.0_EB
             Is_Known_Door(:)   = .False.
             Is_Visible_Door(:) = .False.
             Do ie = 1, n_doors
                If ( EVAC_DOORS(ie)%IMESH == HPT%imesh ) Then
                   Is_Visible_Door(ie) = .True.
                End If
             End Do
             Do ie = 1, n_exits
                If ( EVAC_EXITS(ie)%IMESH == HPT%imesh .And. &
                     .Not. EVAC_EXITS(ie)%COUNT_ONLY ) Then
                   Is_Visible_Door(n_doors+ie) = .True.
                End If
             End Do

             Do ie = 1, HPT%N_VENT_FFIELDS 
                i_tmp = 1
                If (Trim(EVAC_Node_List(HPT%I_DOOR_NODES(ie) &
                     )%Node_Type) == 'Door') Then
                   i_tmp = EVAC_Node_List(HPT%I_DOOR_NODES(ie) &
                        )%Node_Index 
                End If
                If (Trim(EVAC_Node_List(HPT%I_DOOR_NODES(ie) &
                     )%Node_Type) == 'Exit' ) Then
                   i_tmp = n_doors + &
                        EVAC_Node_List(HPT%I_DOOR_NODES(ie) &
                        )%Node_Index 
                End If

                If (HPT%P_VENT_FFIELDS(ie) < 1.0_EB) Then
                   Call Random_number(RN)
                   If ( RN < HPT%P_VENT_FFIELDS(ie) ) Then
                      Is_Known_Door(i_tmp) = .True.
                   Else
                      Is_Known_Door(i_tmp) = .False.
                   End If
                Else
                   Is_Known_Door(i_tmp) = .True.
                End If
             End Do

             ! Find the visible doors 
             Do ie = 1, n_doors + n_exits
                ! Check that the door/exit is on the correct grid:
                If ( Is_Visible_Door(ie) ) Then
                   If (EVAC_Node_List(n_egrids+n_entrys+ie)%Node_Type &
                        == 'Door' ) Then
                      X11 = EVAC_DOORS(ie)%X 
                      Y11 = EVAC_DOORS(ie)%Y
                   Else        ! 'Exit'
                      X11 = EVAC_EXITS(ie-n_doors)%X 
                      Y11 = EVAC_EXITS(ie-n_doors)%Y
                   End If
                   max_fed = 0.0_EB
                   II = Floor(M%CELLSI(Floor((X11-M%XS)*M%RDXINT)) + 1.0_EB  )
                   JJ = Floor(M%CELLSJ(Floor((Y11-M%YS)*M%RDYINT)) + 1.0_EB  )
                   KK = Floor(M%CELLSK(Floor(( 0.5_EB*(M%ZF+M%ZS)-M%ZS)*M%RDZINT)) + 1.0_EB  )
                   XI  = M%CELLSI(Floor((x1_old-M%XS)*M%RDXINT))
                   YJ  = M%CELLSJ(Floor((y1_old-M%YS)*M%RDYINT))
                   IIE = Floor(XI+1.0_EB)
                   JJE = Floor(YJ+1.0_EB)
                   PP_see_door = .True. ! oletusarvo
                   If (Abs(X1_OLD-X11) >= Abs(Y1_OLD-Y11)) Then
                      If ( iie < ii) Then
                         iio = iie
                         jjo = jje
                         iie = ii
                         jje = jj
                         y_o = Y1_OLD
                         x_o = X1_OLD
                         Delta_y = (Y11 - Y1_OLD)
                      Else
                         Delta_y = (Y1_OLD - Y11)
                         iio = ii
                         jjo = jj
                         y_o = Y11
                         x_o = X11
                      End If
                      Delta_x = Abs(X1_OLD - X11)
                      x_now = 0.0_EB
                      PP_see_door_x: Do iii = iio+1, iie-1
                         x_now = x_now + M%DX(iii)
                         y_now = y_o + x_now*(Delta_y/Delta_x)
                         jjj=Floor(M%CELLSJ(Floor((y_now-M%YS) &
                              *M%RDYINT)) + 1.0_EB  )
                         If (max_fed < &
                              M%HUMAN_GRID(iii,jjj)%FED_CO_CO2_O2) Then
                            max_fed = M%HUMAN_GRID(iii,jjj)%FED_CO_CO2_O2
                         End If
                         tim_ic = M%CELL_INDEX(iii,jjj,KK)
                         If (M%SOLID(tim_ic)) Then
                            PP_see_door = .False.
                            Exit PP_see_door_x
                         End If
                      End Do PP_see_door_x
                   Else 
                      If ( jje < jj) Then
                         iio = iie
                         jjo = jje
                         iie = ii
                         jje = jj
                         y_o = Y1_OLD
                         x_o = X1_OLD
                         Delta_x = (X11 - X1_OLD)
                      Else
                         Delta_x = (X1_OLD - X11)
                         iio = ii
                         jjo = jj
                         y_o = Y11
                         x_o = X11
                      End If
                      Delta_y = Abs(Y1_OLD - Y11)
                      y_now = 0.0_EB
                      PP_see_door_y: Do jjj = jjo+1, jje-1
                         y_now = y_now + M%DY(jjj)
                         x_now = x_o + y_now*(Delta_x/Delta_y)
                         iii=Floor(M%CELLSI(Floor((x_now-M%XS)* &
                              M%RDXINT)) + 1.0_EB)
                         If (max_fed < &
                              M%HUMAN_GRID(iii,jjj)%FED_CO_CO2_O2) Then
                            max_fed = M%HUMAN_GRID(iii,jjj)%FED_CO_CO2_O2
                         End If
                         tim_ic = M%CELL_INDEX(iii,jjj,KK)
                         If (M%SOLID(tim_ic)) Then
                            PP_see_door = .False.
                            Exit PP_see_door_y
                         End If
                      End Do PP_see_door_y
                   End If

                   If (PP_see_door) Then
                      If (EVAC_Node_List(n_egrids+n_entrys+ie &
                           )%Node_Type == 'Door') Then
                         If (.Not. EVAC_DOORS(ie)%EXIT_SIGN) Then
                            Is_Visible_Door(ie) = .False.
                         End If
                      End If
                      FED_max_Door(ie) = max_fed
                   Else
                      ! Not visible, use fed at position of the human
                      iie = Floor(M%CELLSI(Floor((x1_old-M%XS)*M%RDXINT)) + 1.0_EB)
                      jje = Floor(M%CELLSJ(Floor((y1_old-M%YS)*M%RDYINT)) + 1.0_EB)
                      Is_Visible_Door(ie) = .False.
                      FED_max_Door(ie) = &
                           M%HUMAN_GRID(iie,jje)%FED_CO_CO2_O2 
                   End If
                End If        ! correct floor
             End Do          ! doors and exits

             ! Save the random known door information
             i_tmp = Count(Is_Known_Door)
             If (j > 0 ) Then
                Group_Known_Doors(j)%N_nodes = i_tmp
                If (Count(Is_Known_Door) > 0 ) Then
                   Allocate(Group_Known_Doors(j)%I_nodes(i_tmp), &
                        STAT=IZERO)
                   Call ChkMemErr('Initialize_Evacuation', &
                        'Group_Known_Doors',IZERO) 
                   i_tmp = 0
                   Do ie = 1, n_doors
                      If (Is_Known_Door(ie)) Then
                         i_tmp = i_tmp + 1
                         Group_Known_Doors(j)%I_nodes(i_tmp) = &
                              EVAC_DOORS(ie)%INODE  
                      End If
                   End Do
                   Do ie = 1, n_exits
                      If (Is_Known_Door(ie+n_doors)) Then
                         i_tmp = i_tmp + 1
                         Group_Known_Doors(j)%I_nodes(i_tmp) = &
                              EVAC_EXITS(ie)%INODE  
                      End If
                   End Do
                End If        ! there are known doors for this group
             Else
                If (j1 > 0) Then
                   Human_Known_Doors(j1)%N_nodes = i_tmp
                   If (Count(Is_Known_Door) > 0 ) Then
                      Allocate(Human_Known_Doors(j1)%I_nodes(i_tmp), &
                           STAT=IZERO)
                      Call ChkMemErr('Initialize_Evacuation', &
                           'Human_Known_Doors',IZERO) 
                      i_tmp = 0
                      Do ie = 1, n_doors
                         If (Is_Known_Door(ie)) Then
                            i_tmp = i_tmp + 1
                            Human_Known_Doors(j1)%I_nodes(i_tmp) = &
                                 EVAC_DOORS(ie)%INODE  
                         End If
                      End Do
                      Do ie = 1, n_exits
                         If (Is_Known_Door(ie+n_doors)) Then
                            i_tmp = i_tmp + 1
                            Human_Known_Doors(j1)%I_nodes(i_tmp) = &
                                 EVAC_EXITS(ie)%INODE  
                         End If
                      End Do
                   End If        ! there are known doors for this group
                End If
             End If

             Do ie = 1, n_doors
                If ( EVAC_DOORS(ie)%IMESH /= HPT%imesh ) Then
                   Is_Known_Door(ie) = .False.
                End If
             End Do
             Do ie = 1, n_exits
                If ( EVAC_EXITS(ie)%IMESH /= HPT%imesh .Or. &
                     EVAC_EXITS(ie)%COUNT_ONLY ) Then
                   Is_Known_Door(n_doors+ie) = .False.
                End If
             End Do


             Do ie = 1, n_doors
                If ( EVAC_DOORS(ie)%TIME_OPEN > 0.0_EB ) Then
                   Is_Visible_Door(ie) = .False.
                   Is_Known_Door(ie) = .False.
                End If
             End Do
             Do ie = 1, n_exits
                If ( EVAC_EXITS(ie)%TIME_OPEN > 0.0_EB .And. &
                     .Not. EVAC_EXITS(ie)%COUNT_ONLY ) Then
                   Is_Visible_Door(n_doors+ie) = .False.
                   Is_Known_Door(n_doors+ie) = .False.
                End If
             End Do

             If (Any(Is_Known_Door) .Or. Any(Is_Visible_Door)) Then
                i_tmp   = 0
                L2_min = Huge(L2_min)
                Do ie = 1, n_doors + n_exits
                   If ( Is_Known_Door(ie) .And. &
                        Is_Visible_Door(ie) ) Then
                      x_o = 0.0_EB
                      y_o = 0.0_EB
                      If (Trim(EVAC_Node_List(n_egrids+n_entrys+ie &
                           )%Node_Type) == 'Door' ) Then
                         x_o = EVAC_DOORS( EVAC_Node_List( &
                              ie+n_egrids+n_entrys)%Node_Index )%X
                         y_o = EVAC_DOORS( EVAC_Node_List( &
                              ie+n_egrids+n_entrys)%Node_Index )%Y
                      Else      ! 'Exit'
                         x_o = EVAC_EXITS( EVAC_Node_List( &
                              ie+n_egrids+n_entrys)%Node_Index )%X
                         y_o = EVAC_EXITS( EVAC_Node_List( &
                              ie+n_egrids+n_entrys)%Node_Index )%Y
                      End If
                      L2_tmp = 0.0_EB ! initialization, no smoke yet

                      If ( ((x_o-HR%X)**2 + (y_o-HR%Y)**2) < L2_min &
                           .And. L2_tmp < Abs(FED_DOOR_CRIT)) Then
                         L2_min = (x_o-HR%X)**2 + (y_o-HR%Y)**2 
                         L2_min = Max(0.0_EB,L2_min)
                         i_tmp = ie
                      End If
                   End If
                End Do
                If (i_tmp > 0 ) Then
                   ! Known and visible door, no smoke
                   color_index = 0
                   If (EVAC_Node_List(n_egrids+n_entrys+i_tmp &
                        )%Node_Type == 'Door' ) Then
                      GROUP_FFIELD_NAME = &
                           Trim(EVAC_DOORS(i_tmp)%VENT_FFIELD)
                      GROUP_FFIELD_I = &
                           EVAC_DOORS(i_tmp)%I_VENT_FFIELD
                   Else        ! 'Exit'
                      GROUP_FFIELD_NAME = &
                           Trim(EVAC_EXITS(i_tmp-n_doors)%VENT_FFIELD)
                      GROUP_FFIELD_I = &
                           EVAC_EXITS(i_tmp-n_doors)%I_VENT_FFIELD
                   End If
                Else
                   ! No visible known door available, try non-visible known doors
                   i_tmp   = 0
                   L2_min = Huge(L2_min)
                   Do ie = 1, n_doors + n_exits
                      If ( Is_Known_Door(ie) .And. &
                           .Not. Is_Visible_Door(ie) ) Then
                         x_o = 0.0_EB
                         y_o = 0.0_EB
                         If (EVAC_Node_List(n_egrids+n_entrys+ie &
                              )%Node_Type == 'Door' ) Then
                            x_o = EVAC_DOORS( EVAC_Node_List(ie+ &
                                 n_egrids+n_entrys)%Node_Index )%X
                            y_o = EVAC_DOORS( EVAC_Node_List(ie+ &
                                 n_egrids+n_entrys)%Node_Index )%Y
                         Else    ! 'Exit'
                            x_o = EVAC_EXITS( EVAC_Node_List(ie+ &
                                 n_egrids+n_entrys)%Node_Index )%X
                            y_o = EVAC_EXITS( EVAC_Node_List(ie+ &
                                 n_egrids+n_entrys)%Node_Index )%Y
                         End If
                         L2_tmp = 0.0_EB ! initialization, no smoke yet
                         If ( ((x_o-HR%X)**2 + (y_o-HR%Y)**2) < &
                              L2_min .And. L2_tmp < &
                              Abs(FED_DOOR_CRIT)) Then
                            L2_min = (x_o-HR%X)**2 + (y_o-HR%Y)**2 
                            L2_min = Max(0.0_EB,L2_min)
                            i_tmp = ie
                         End If
                      End If
                   End Do
                   If (i_tmp > 0 ) Then
                      ! Non-visible known door, no smoke
                      color_index = 1
                      If (EVAC_Node_List( &
                           n_egrids+n_entrys+i_tmp)%Node_Type &
                           == 'Door' ) Then
                         GROUP_FFIELD_NAME = &
                              Trim(EVAC_DOORS(i_tmp)%VENT_FFIELD)
                         GROUP_FFIELD_I = &
                              EVAC_DOORS(i_tmp)%I_VENT_FFIELD
                      Else      ! 'Exit'
                         GROUP_FFIELD_NAME = Trim( &
                              EVAC_EXITS(i_tmp-n_doors)%VENT_FFIELD)
                         GROUP_FFIELD_I = &
                              EVAC_EXITS(i_tmp-n_doors)%I_VENT_FFIELD
                      End If
                   Else
                      ! known doors with no smoke have not been found
                      i_tmp   = 0
                      L2_min = Huge(L2_min)
                      Do ie = 1, n_doors + n_exits
                         If (Is_Visible_Door(ie)) Then
                            x_o = 0.0_EB
                            y_o = 0.0_EB
                            If (EVAC_Node_List(n_egrids+n_entrys+ &
                                 ie)%Node_Type == 'Door' ) Then
                               x_o = EVAC_DOORS( EVAC_Node_List(ie+ &
                                    n_egrids+n_entrys)%Node_Index )%X
                               y_o = EVAC_DOORS( EVAC_Node_List(ie+ &
                                    n_egrids+n_entrys)%Node_Index )%Y
                            Else  ! 'Exit'
                               x_o = EVAC_EXITS( EVAC_Node_List(ie+ &
                                    n_egrids+n_entrys)%Node_Index )%X
                               y_o = EVAC_EXITS( EVAC_Node_List(ie+ &
                                    n_egrids+n_entrys)%Node_Index )%Y
                            End If
                            L2_tmp = 0.0_EB ! initialization, no smoke yet
                            If ( ((x_o-HR%X)**2+(y_o-HR%Y)**2) < L2_min &
                                 .And. L2_tmp < &
                                 Abs(FED_DOOR_CRIT)) Then
                               L2_min = (x_o-HR%X)**2 + (y_o-HR%Y)**2 
                               L2_min = Max(0.0_EB,L2_min)
                               i_tmp = ie
                            End If
                         End If
                      End Do
                      If (i_tmp > 0 ) Then
                         ! No smoke, visible door (not known)
                         color_index = 2
                         If (EVAC_Node_List( &
                              n_egrids+n_entrys+i_tmp)%Node_Type &
                              == 'Door' ) Then
                            GROUP_FFIELD_NAME = &
                                 Trim(EVAC_DOORS(i_tmp)%VENT_FFIELD)
                            GROUP_FFIELD_I = &
                                 EVAC_DOORS(i_tmp)%I_VENT_FFIELD
                         Else    ! 'Exit'
                            GROUP_FFIELD_NAME = Trim( &
                                 EVAC_EXITS(i_tmp-n_doors)%VENT_FFIELD)
                            GROUP_FFIELD_I = &
                                 EVAC_EXITS(i_tmp-n_doors)%I_VENT_FFIELD
                         End If
                      Else
                         ! Now we have smoke and some visible or known doors
                         i_tmp   = 0
                         L2_min = Huge(L2_min)
                         Do ie = 1, n_doors + n_exits
                            If (Is_Visible_Door(ie) .Or. &
                                 Is_Known_Door(ie) ) Then
                               x_o = 0.0_EB
                               y_o = 0.0_EB
                               If (EVAC_Node_List(n_egrids+n_entrys+ &
                                    ie)%Node_Type == 'Door' ) Then
                                  x_o = EVAC_DOORS(EVAC_Node_List(ie+ &
                                       n_egrids+n_entrys)%Node_Index)%X
                                  y_o = EVAC_DOORS( EVAC_Node_List(ie+ &
                                       n_egrids+n_entrys)%Node_Index)%Y
                               Else ! 'Exit'
                                  x_o = EVAC_EXITS( EVAC_Node_List(ie+ &
                                       n_egrids+n_entrys)%Node_Index)%X
                                  y_o = EVAC_EXITS( EVAC_Node_List(ie+ &
                                       n_egrids+n_entrys)%Node_Index)%Y
                               End If
                               L2_tmp = HR%IntDose + FED_max_Door(ie) * &
                                    Sqrt((HR%X-x_o)**2 + (HR%Y-y_o)**2)/ &
                                    HR%Speed
                               If (L2_tmp < L2_min) Then
                                  L2_min = L2_tmp
                                  i_tmp = ie
                               End If
                            End If
                         End Do

                         If (i_tmp > 0 .And. L2_min < 1.0_EB) Then
                            ! Not too much smoke, i.e., non-lethal amount (known or visible doors)
                            If (EVAC_Node_List( &
                                 n_egrids+n_entrys+i_tmp)%Node_Type &
                                 == 'Door' ) Then
                               GROUP_FFIELD_NAME = &
                                    Trim(EVAC_DOORS(i_tmp)%VENT_FFIELD)
                               GROUP_FFIELD_I = &
                                    EVAC_DOORS(i_tmp)%I_VENT_FFIELD
                            Else  ! 'Exit'
                               GROUP_FFIELD_NAME = Trim( &
                                    EVAC_EXITS(i_tmp-n_doors)%VENT_FFIELD)
                               GROUP_FFIELD_I = &
                                    EVAC_EXITS(i_tmp-n_doors)%I_VENT_FFIELD
                            End If
                            If (Is_Known_Door(i_tmp) .And. &
                                 Is_Visible_Door(i_tmp)) color_index = 3
                            If (Is_Known_Door(i_tmp) .And. .Not. &
                                 Is_Visible_Door(i_tmp)) color_index = 4
                            If (.Not. Is_Known_Door(i_tmp) .And. &
                                 Is_Visible_Door(i_tmp)) color_index = 5
                         Else
                            ! No door found, use the main evac grid ffield or evac-line 
                            i_tmp = 0 ! no door found
                            GROUP_FFIELD_NAME = Trim(HPT%GRID_NAME)
                            GROUP_FFIELD_I    = HPT%I_VENT_FFIELDS(0)
                            color_index = 6
                         End If  ! case 4
                      End If    ! case 3
                   End If      ! case 2
                End If        ! case 1
                If (Color_Method == 4 ) Then
                   color_index =  6 ! default, cyan
                   If (i_tmp > 0 .And. i_tmp <= n_doors) &
                        color_index = EVAC_DOORS(i_tmp)%COLOR_INDEX
                   If (i_tmp > n_doors .And. i_tmp <=  n_doors + n_exits) &
                        color_index = EVAC_EXITS(i_tmp-n_doors)%COLOR_INDEX
                End If

             Else ! No known/visible door
                i_tmp = 0 ! no door found
                GROUP_FFIELD_NAME = Trim(HPT%GRID_NAME)
                GROUP_FFIELD_I    = HPT%I_VENT_FFIELDS(0)
                color_index = 6
                If (Color_Method == 4) color_index = HPT%COLOR_INDEX
             End If

             HR%I_Target = i_tmp
             If (i_tmp > 0 .And. .Not. Is_Visible_Door(Max(1,i_tmp))) Then
                ! I_Target >0: visible, <0: not visible
                HR%I_Target = -i_tmp
             End If
             If (j > 0) Group_Known_Doors(j)%I_Target = HR%I_Target
             HR%FFIELD_NAME = Trim(GROUP_FFIELD_NAME)
             HR%I_FFIELD    = GROUP_FFIELD_I
             If (COLOR_METHOD == 5 ) HR%COLOR_INDEX = color_index
             If (COLOR_METHOD == 4 ) HR%COLOR_INDEX = color_index
             If (j > 0) Color_Tmp(j) = color_index
             If (j > 0) Group_List(j)%GROUP_I_FFIELDS(i_egrid) &
                  = GROUP_FFIELD_I
          Else              ! this group has already a flow field
             ! This group has already tried to change the field
             If (COLOR_METHOD == 5 .And. j > 0) &
                  HR%COLOR_INDEX = Color_Tmp(j)
             If (COLOR_METHOD == 4 .And. j > 0) &
                  HR%COLOR_INDEX = Color_Tmp(j)
             HR%I_FFIELD    = Group_List(j)%GROUP_I_FFIELDS(i_egrid)
             HR%FFIELD_NAME = Trim(MESH_NAME(HR%I_FFIELD))
             HR%I_Target = Group_Known_Doors(j)%I_Target
          End If            ! first member of a group or a lonely soul
       End Do              ! 1, n_humans
!!$       End If                ! main evac grid
    End Do                  ! 1, nmeshes
    Deallocate(Color_Tmp)
    Deallocate(FED_max_Door)
    Deallocate(Is_Visible_Door)
    Deallocate(Is_Known_Door)

    ! Initialize the group_i_ffields
    i_egrid = 0
    Do nom = 1, NMESHES
       If ( .Not.(EVACUATION_ONLY(NoM) .And. EVACUATION_GRID(NoM)) ) Cycle
!!$       If ( EVACUATION_ONLY(NoM) .And. EVACUATION_GRID(NoM) ) Then
       M => MESHES(NOM)
       i_egrid = i_egrid+1
       Do i = 1, M%N_HUMANS
          HR => M%HUMAN(I)

          ! Check spectator stands, correct the z-coordiante
          SS_Loop: Do j = 1, n_sstands
             ESS => EVAC_SSTANDS(j)
             If (ESS%IMESH == nom .And. &
                  (ESS%X1 <= HR%X .And. ESS%X2 >= HR%X) .And. &
                  (ESS%Y1 <= HR%Y .And. ESS%Y2 >= HR%Y) ) Then
                Select Case (ESS%IOR)
                Case(-1)
                   HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + &
                        (ESS%H-ESS%H0)*Abs(ESS%X1-HR%X)/Abs(ESS%X1-ESS%X2)
                Case(+1)
                   HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + &
                        (ESS%H-ESS%H0)*Abs(ESS%X2-HR%X)/Abs(ESS%X1-ESS%X2)
                Case(-2)
                   HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + &
                        (ESS%H-ESS%H0)*Abs(ESS%Y1-HR%Y)/Abs(ESS%Y1-ESS%Y2)
                Case(+2)
                   HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + &
                        (ESS%H-ESS%H0)*Abs(ESS%Y2-HR%Y)/Abs(ESS%Y1-ESS%Y2)
                End Select
                Exit SS_Loop
             End If
          End Do SS_Loop

          j = Max(0,HR%GROUP_ID)
          Group_List(j)%IEL = HR%IEL

          Write (LU_ERR,fmt='(i6,5f8.2,3f6.2,i6,i3,i2)') HR%ILABEL, &
               HR%X, HR%Y, HR%Z, HR%Tpre, HR%Tdet,2.0_EB*HR%Radius, &
               HR%Speed, HR%Tau, HR%GROUP_ID, HR%i_ffield, &
               Abs(HR%I_Target)+n_egrids+n_entrys
       End Do
!!$          End If
    End Do

    !! End If ! any evacuation grid
    !
  End Subroutine INIT_EVAC_GROUPS
!
  Subroutine EVAC_MESH_EXCHANGE(T,T_Save,I_mode, ICYC)
    Implicit None
    !
    Real(EB), Intent(IN) :: T
    Integer, Intent(IN) :: I_mode, ICYC
    Real(EB) T_Save
    !
    Integer nm, nom, i, j, i1, j1, k1
    Integer ios
    Real(EB) tmp_1, y_extra, Y_MF_INT
    Logical L_use_fed, L_fed_read, L_fed_save
    Real(EB) DT_Save
    integer(4) ibar_tmp, jbar_tmp, kbar_tmp, n_tmp
    Real(FB) tmpout1, tmpout2, tmpout3, tmpout4, t_tmp, dt_tmp
    Real(FB) tmpout5, tmpout6, tmpout7, tmpout8
    Real(EB) Z_1,Z_2,Z_3
    !
    If (ICYC < 1) Return
    If (.NOT.Any(EVACUATION_GRID)) RETURN
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
    If ( T >= 0.0_EB .And. Real(T,FB) >= Real(T_Save,FB) ) Then
       T_Save = T + DT_Save
       L_use_fed = .True.
    End If

    L_use_fed = L_use_fed .And. (L_fed_read .Or. L_fed_save)

    If (L_use_fed) Then
       If (L_fed_save) Then
          Write (LU_EVACFED) Real(T,FB), Real(DT_Save,FB)
       Else
          Read (LU_EVACFED,End=324,Iostat=ios) t_tmp, dt_tmp
          T_Save = t_tmp + dt_tmp ! next time point in the file
       End If

       ! Next loop interpolates fire mesh (soot+fed) into human_grids and
       ! saves it to the disk. Or it reads FED+soot from the disk.
       MESH_LOOP: Do NM=1,NMESHES
          If ( .Not.(EVACUATION_GRID(NM) .And. EVACUATION_ONLY(NM)) ) Cycle
!!$          If ( EVACUATION_GRID(NM) .And. EVACUATION_ONLY(NM) ) Then
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
             If (ibar_tmp /= IBAR .Or. jbar_tmp /= JBAR .Or. &
                  n_tmp < 4 ) Then
                Close (LU_EVACFED)
                Call SHUTDOWN('ERROR: Problems to read the FED file')
             End If

          End If

          Do i = 1, IBAR
             Do j= 1, JBAR

                If (L_fed_save) Then

                   If ( HUMAN_GRID(i,j)%IMESH > 0 ) Then
                      i1 = HUMAN_GRID(i,j)%II
                      j1 = HUMAN_GRID(i,j)%JJ
                      k1 = HUMAN_GRID(i,j)%KK
                      nom = HUMAN_GRID(i,j)%IMESH
                      Z_1 = MESHES(nom)%YY(I1,J1,K1,I_FUEL)
                      If (CO_PRODUCTION) Then
                         Z_2 = MESHES(nom)%YY(I1,J1,K1,I_PROG_CO)
                      Else
                         Z_2 = 0._EB
                      End If
                      Z_3 = MESHES(nom)%YY(I1,J1,K1,I_PROG_F)
                      y_extra = MESHES(nom)%Y_SUM(I1,J1,K1)  ! extra species mass fraction
                      ! Mass fraction array ==> soot density (mg/m3)
                      ! Next is for soot (mg/m3)
                      Call GET_MASS_FRACTION2(Z_1,Z_2,Z_3,SOOT_INDEX,y_extra,Y_MF_INT)
                      tmp_1 = Y_MF_INT*MESHES(nom)%RHO(i1,j1,k1)*1.E6_EB
                      HUMAN_GRID(i,j)%SOOT_DENS = tmp_1
                      ! Calculate Purser's fractional effective dose (FED)
                      ! Note: Purser uses minutes, here dt is in seconds
                      !       fed_dose = fed_lco*fed_vco2 + fed_lo
                      ! CO:  (3.317E-5*RMV*t)/D
                      !      [RMV]=ltr/min, D=30% COHb concentration at incapacitation
                      ! Next is for CO (ppm)
                      Call GET_MASS_FRACTION2(Z_1,Z_2,Z_3,CO_INDEX,y_extra,Y_MF_INT)
                      tmp_1 = RCON_MF(CO_INDEX)*Y_MF_INT*1.E6_EB/MESHES(nom)%RSUM(I1,J1,K1)
                      HUMAN_GRID(i,j)%FED_CO_CO2_O2 = 3.317E-5_EB*25.0_EB* &
                           tmp_1**(1.036_EB)/(30.0_EB*60.0_EB)
                      ! VCO2: CO2-induced hyperventilation
                      !      exp(0.1903*c_CO2(%) + 2.0004)
                      ! Next is for CO2
                      Call GET_MASS_FRACTION2(Z_1,Z_2,Z_3,CO2_INDEX,y_extra,Y_MF_INT)
                      tmp_1 = RCON_MF(CO2_INDEX)*Y_MF_INT/MESHES(nom)%RSUM(I1,J1,K1)
                      HUMAN_GRID(i,j)%FED_CO_CO2_O2 = &
                           HUMAN_GRID(i,j)%FED_CO_CO2_O2*Exp( 0.1903_EB*tmp_1 &
                           *100.0_EB + 2.0004_EB )/7.1_EB
                      ! LO: low oxygen
                      ! Next is for O2
                      Call GET_MASS_FRACTION2(Z_1,Z_2,Z_3,O2_INDEX,y_extra,Y_MF_INT)
                      tmp_1 = RCON_MF(O2_INDEX)*Y_MF_INT/MESHES(nom)%RSUM(I1,J1,K1)
                      If ( tmp_1 < 0.20_EB ) Then
                         HUMAN_GRID(i,j)%FED_CO_CO2_O2 = &
                              HUMAN_GRID(i,j)%FED_CO_CO2_O2 + 1.0_EB  / &
                              (60.0_EB*Exp(8.13_EB-0.54_EB*(20.9_EB-100.0_EB*tmp_1)) )
                      End If
                      ! Gas temperature, ind=5, C
                      HUMAN_GRID(i,j)%TMP_G  = MESHES(nom)%TMP(i1,j1,k1)
                      ! Radiant intensity, ind=18, kW/m2 (no -sigma*Tamb^4 term)
                      HUMAN_GRID(i,j)%RADINT = &
                           Max(MESHES(nom)%UII(i1,j1,k1),4.0_EB*SIGMA*TMPA4)


                   End If ! imesh > 0, i.e. fire grid found

                   ! Save Fed, Soot, Temp(C), and RadInt
                   Write (LU_EVACFED) &
                        Real(HUMAN_GRID(i,j)%FED_CO_CO2_O2,FB), &
                        Real(HUMAN_GRID(i,j)%SOOT_DENS,FB), &
                        Real(HUMAN_GRID(i,j)%TMP_G,FB), &
                        Real(HUMAN_GRID(i,j)%RADINT,FB)

                Else     ! Read FED from a file
                   ! Read Fed, Soot, Temp(C), and RadInt
                   Read (LU_EVACFED,Iostat=ios) tmpout1, tmpout2, tmpout3, tmpout4
                   If (ios.Ne.0) Then
                      Write(MESSAGE,'(A)') 'ERROR: Evac Mesh Exchange: FED READ ERROR'
                      Close (LU_EVACFED)
                      Call SHUTDOWN(MESSAGE)
                   End If
                   HUMAN_GRID(i,j)%FED_CO_CO2_O2 = tmpout1
                   HUMAN_GRID(i,j)%SOOT_DENS = tmpout2
                   HUMAN_GRID(i,j)%TMP_G = tmpout3
                   HUMAN_GRID(i,j)%RADINT = tmpout4

                End If   ! Calculate and save FED

             End Do     ! j=1,JBAR
          End Do       ! i=1,IBAR

!!$          End If          ! Main evac grid

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
                Z_1 = MESHES(nom)%YY(I1,J1,K1,I_FUEL)
                If (CO_PRODUCTION) Then
                   Z_2 = MESHES(nom)%YY(I1,J1,K1,I_PROG_CO)
                Else
                   Z_2 = 0._EB
                End If
                Z_3 = MESHES(nom)%YY(I1,J1,K1,I_PROG_F)
                y_extra = MESHES(nom)%Y_SUM(I1,J1,K1)  ! extra species mass fraction

                ! Mass fraction array ==> soot density (mg/m3)
                ! Next is for soot (mg/m3)
                Call GET_MASS_FRACTION2(Z_1,Z_2,Z_3,SOOT_INDEX,y_extra,Y_MF_INT)
                tmp_1 = Y_MF_INT*MESHES(nom)%RHO(i1,j1,k1)*1.E6_EB
                EVAC_CORRS(i)%SOOT_DENS(1) = tmp_1

                ! Calculate Purser's fractional effective dose (FED)

                ! Next is for CO (ppm)
                Call GET_MASS_FRACTION2(Z_1,Z_2,Z_3,CO_INDEX,y_extra,Y_MF_INT)
                tmp_1 = RCON_MF(CO_INDEX)*Y_MF_INT*1.E6_EB/MESHES(nom)%RSUM(I1,J1,K1)
                EVAC_CORRS(i)%FED_CO_CO2_O2(1) = 3.317E-5_EB*25.0_EB* &
                     tmp_1**(1.036_EB)/(30.0_EB*60.0_EB)

                ! Next is for CO2
                Call GET_MASS_FRACTION2(Z_1,Z_2,Z_3,CO2_INDEX,y_extra,Y_MF_INT)
                tmp_1 = RCON_MF(CO2_INDEX)*Y_MF_INT/MESHES(nom)%RSUM(I1,J1,K1)
                EVAC_CORRS(i)%FED_CO_CO2_O2(1) = &
                     EVAC_CORRS(i)%FED_CO_CO2_O2(1)*Exp( 0.1903_EB*tmp_1 &
                     *100.0_EB + 2.0004_EB )/7.1_EB

                ! Next is for O2
                Call GET_MASS_FRACTION2(Z_1,Z_2,Z_3,O2_INDEX,y_extra,Y_MF_INT)
                tmp_1 = RCON_MF(O2_INDEX)*Y_MF_INT/MESHES(nom)%RSUM(I1,J1,K1)
                If ( tmp_1 < 0.20_EB ) Then
                   EVAC_CORRS(i)%FED_CO_CO2_O2(1) = &
                        EVAC_CORRS(i)%FED_CO_CO2_O2(1) + 1.0_EB  / &
                        (60.0_EB*Exp(8.13_EB-0.54_EB*(20.9_EB-100.0_EB*tmp_1)) )
                End If

                ! Gas temperature, ind=5, C
                EVAC_CORRS(i)%TMP_G(1)  = MESHES(nom)%TMP(i1,j1,k1)

                ! Radiant intensity, ind=18, kW/m2 (no -sigma*Tamb^4 term)
                EVAC_CORRS(i)%RADINT(1) = &
                     Max(MESHES(nom)%UII(i1,j1,k1),4.0_EB*SIGMA*TMPA4)

             Else
                ! No fed_mesh found
                EVAC_CORRS(i)%FED_CO_CO2_O2(1) = 0.0_EB
                EVAC_CORRS(i)%SOOT_DENS(1) = 0.0_EB
                EVAC_CORRS(i)%TMP_G(1) = 0.0_EB
                EVAC_CORRS(i)%RADINT(1) = 0.0_EB
             End If                ! fed_mesh > 0, i.e. fire grid found

             If ( EVAC_CORRS(i)%FED_MESH2 > 0 ) Then
                i1 = EVAC_CORRS(i)%II(2)
                j1 = EVAC_CORRS(i)%JJ(2)
                k1 = EVAC_CORRS(i)%KK(2)
                nom = EVAC_CORRS(i)%FED_MESH2
                Z_1 = MESHES(nom)%YY(I1,J1,K1,I_FUEL)
                If (CO_PRODUCTION) Then
                   Z_2 = MESHES(nom)%YY(I1,J1,K1,I_PROG_CO)
                Else
                   Z_2 = 0._EB
                End If
                Z_3 = MESHES(nom)%YY(I1,J1,K1,I_PROG_F)
                y_extra = MESHES(nom)%Y_SUM(I1,J1,K1)  ! extra species mass fraction

                ! Mass fraction array ==> soot density (mg/m3)
                ! Next is for soot (mg/m3)
                Call GET_MASS_FRACTION2(Z_1,Z_2,Z_3,SOOT_INDEX,y_extra,Y_MF_INT)
                tmp_1 = Y_MF_INT*MESHES(nom)%RHO(i1,j1,k1)*1.E6_EB
                EVAC_CORRS(i)%SOOT_DENS(2) = tmp_1

                ! Calculate Purser's fractional effective dose (FED)

                ! Next is for CO (ppm)
                Call GET_MASS_FRACTION2(Z_1,Z_2,Z_3,CO_INDEX,y_extra,Y_MF_INT)
                tmp_1 = RCON_MF(CO_INDEX)*Y_MF_INT*1.E6_EB/MESHES(nom)%RSUM(I1,J1,K1)
                EVAC_CORRS(i)%FED_CO_CO2_O2(2) = 3.317E-5_EB*25.0_EB* &
                     tmp_1**(1.036_EB)/(30.0_EB*60.0_EB)

                ! Next is for CO2
                Call GET_MASS_FRACTION2(Z_1,Z_2,Z_3,CO2_INDEX,y_extra,Y_MF_INT)
                tmp_1 = RCON_MF(CO2_INDEX)*Y_MF_INT/MESHES(nom)%RSUM(I1,J1,K1)
                EVAC_CORRS(i)%FED_CO_CO2_O2(2) = &
                     EVAC_CORRS(i)%FED_CO_CO2_O2(2)*Exp( 0.1903_EB*tmp_1 &
                     *100.0_EB + 2.0004_EB )/7.1_EB

                ! Next is for O2
                Call GET_MASS_FRACTION2(Z_1,Z_2,Z_3,O2_INDEX,y_extra,Y_MF_INT)
                tmp_1 = RCON_MF(O2_INDEX)*Y_MF_INT/MESHES(nom)%RSUM(I1,J1,K1)
                If ( tmp_1 < 0.20_EB ) Then
                   EVAC_CORRS(i)%FED_CO_CO2_O2(2) = &
                        EVAC_CORRS(i)%FED_CO_CO2_O2(2) + 1.0_EB  / &
                        (60.0_EB*Exp(8.13_EB-0.54_EB*(20.9_EB-100.0_EB*tmp_1)) )
                End If
                ! Gas temperature, ind=5, C
                EVAC_CORRS(i)%TMP_G(2)  = MESHES(nom)%TMP(i1,j1,k1)
                ! Radiant intensity, ind=18, kW/m2 (no -sigma*Tamb^4 term)
                EVAC_CORRS(i)%RADINT(2) = &
                     Max(MESHES(nom)%UII(i1,j1,k1),4.0_EB*SIGMA*TMPA4)

             Else
                ! No fed_mesh2 found
                EVAC_CORRS(i)%FED_CO_CO2_O2(2) = 0.0_EB
                EVAC_CORRS(i)%SOOT_DENS(2) = 0.0_EB
                EVAC_CORRS(i)%TMP_G(2) = 0.0_EB
                EVAC_CORRS(i)%RADINT(2) = 0.0_EB
             End If                ! fed_mesh2 > 0, i.e. fire grid found

             ! Save Fed, Soot, Temp(C), and RadInt
             Write (LU_EVACFED) &
                  Real(EVAC_CORRS(i)%FED_CO_CO2_O2(1),FB), &
                  Real(EVAC_CORRS(i)%SOOT_DENS(1),FB), &
                  Real(EVAC_CORRS(i)%TMP_G(1),FB), &
                  Real(EVAC_CORRS(i)%RADINT(1),FB), &
                  Real(EVAC_CORRS(i)%FED_CO_CO2_O2(2),FB), &
                  Real(EVAC_CORRS(i)%SOOT_DENS(2),FB), &
                  Real(EVAC_CORRS(i)%TMP_G(2),FB), &
                  Real(EVAC_CORRS(i)%RADINT(2),FB)

          Else                    ! Read FED from a file
             ! Read Fed, Soot, Temp(C), and RadInt
             Read (LU_EVACFED,Iostat=ios) tmpout1, tmpout2, tmpout3, tmpout4, &
                  tmpout5, tmpout6, tmpout7, tmpout8
             If (ios.Ne.0) Then
                Write(MESSAGE,'(A)') 'ERROR: Evac Mesh Exchange: FED READ ERROR'
                Close (LU_EVACEFF)
                Call SHUTDOWN(MESSAGE)
             End If
             EVAC_CORRS(i)%FED_CO_CO2_O2(1) = tmpout1
             EVAC_CORRS(i)%SOOT_DENS(1) = tmpout2
             EVAC_CORRS(i)%TMP_G(1) = tmpout3
             EVAC_CORRS(i)%RADINT(1) = tmpout4
             EVAC_CORRS(i)%FED_CO_CO2_O2(2) = tmpout5
             EVAC_CORRS(i)%SOOT_DENS(2) = tmpout6
             EVAC_CORRS(i)%TMP_G(2) = tmpout7
             EVAC_CORRS(i)%RADINT(2) = tmpout8

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
                Z_1 = MESHES(nom)%YY(I1,J1,K1,I_FUEL)
                If (CO_PRODUCTION) Then
                   Z_2 = MESHES(nom)%YY(I1,J1,K1,I_PROG_CO)
                Else
                   Z_2 = 0._EB
                End If
                Z_3 = MESHES(nom)%YY(I1,J1,K1,I_PROG_F)
                y_extra = MESHES(nom)%Y_SUM(I1,J1,K1)  ! extra species mass fraction

                ! Mass fraction array ==> soot density (mg/m3)
                ! Next is for soot (mg/m3)
                Call GET_MASS_FRACTION2(Z_1,Z_2,Z_3,SOOT_INDEX,y_extra,Y_MF_INT)
                tmp_1 = Y_MF_INT*MESHES(nom)%RHO(i1,j1,k1)*1.E6_EB
                EVAC_DOORS(i)%SOOT_DENS = tmp_1

                ! Calculate Purser's fractional effective dose (FED)

                ! Next is for CO (ppm)
                Call GET_MASS_FRACTION2(Z_1,Z_2,Z_3,CO_INDEX,y_extra,Y_MF_INT)
                tmp_1 = RCON_MF(CO_INDEX)*Y_MF_INT*1.E6_EB/MESHES(nom)%RSUM(I1,J1,K1)
                EVAC_DOORS(i)%FED_CO_CO2_O2 = 3.317E-5_EB*25.0_EB* &
                     tmp_1**(1.036_EB)/(30.0_EB*60.0_EB)

                ! Next is for CO2
                Call GET_MASS_FRACTION2(Z_1,Z_2,Z_3,CO2_INDEX,y_extra,Y_MF_INT)
                tmp_1 = RCON_MF(CO2_INDEX)*Y_MF_INT/MESHES(nom)%RSUM(I1,J1,K1)
                EVAC_DOORS(i)%FED_CO_CO2_O2 = &
                     EVAC_DOORS(i)%FED_CO_CO2_O2*Exp( 0.1903_EB*tmp_1 &
                     *100.0_EB + 2.0004_EB )/7.1_EB

                ! Next is for O2
                Call GET_MASS_FRACTION2(Z_1,Z_2,Z_3,O2_INDEX,y_extra,Y_MF_INT)
                tmp_1 = RCON_MF(O2_INDEX)*Y_MF_INT/MESHES(nom)%RSUM(I1,J1,K1)
                If ( tmp_1 < 0.20_EB ) Then
                   EVAC_DOORS(i)%FED_CO_CO2_O2 = &
                        EVAC_DOORS(i)%FED_CO_CO2_O2 + 1.0_EB  / &
                        (60.0_EB*Exp(8.13_EB-0.54_EB*(20.9_EB-100.0_EB*tmp_1)) )
                End If

                ! Gas temperature, ind=5, C
                EVAC_DOORS(i)%TMP_G  = MESHES(nom)%TMP(i1,j1,k1)

                ! Radiant intensity, ind=18, kW/m2 (no -sigma*Tamb^4 term)
                EVAC_DOORS(i)%RADINT = &
                     Max(MESHES(nom)%UII(i1,j1,k1),4.0_EB*SIGMA*TMPA4)

             Else
                ! No fed_mesh found
                EVAC_DOORS(i)%FED_CO_CO2_O2 = 0.0_EB
                EVAC_DOORS(i)%SOOT_DENS = 0.0_EB
                EVAC_DOORS(i)%TMP_G = 0.0_EB
                EVAC_DOORS(i)%RADINT = 0.0_EB
             End If                ! fed_mesh > 0, i.e. fire grid found

             ! Save Fed, Soot, Temp(C), and RadInt
             Write (LU_EVACFED) &
                  Real(EVAC_DOORS(i)%FED_CO_CO2_O2,FB), &
                  Real(EVAC_DOORS(i)%SOOT_DENS,FB), &
                  Real(EVAC_DOORS(i)%TMP_G,FB), &
                  Real(EVAC_DOORS(i)%RADINT,FB)

          Else                    ! Read FED from a file
             ! Read Fed, Soot, Temp(C), and RadInt
             Read (LU_EVACFED,Iostat=ios) tmpout1, tmpout2, tmpout3, tmpout4
             If (ios.Ne.0) Then
                Write(MESSAGE,'(A)') 'ERROR: Evac Mesh Exchange: FED READ ERROR'
                Close (LU_EVACEFF)
                Call SHUTDOWN(MESSAGE)
             End If
             EVAC_DOORS(i)%FED_CO_CO2_O2 = tmpout1
             EVAC_DOORS(i)%SOOT_DENS = tmpout2
             EVAC_DOORS(i)%TMP_G = tmpout3
             EVAC_DOORS(i)%RADINT = tmpout4
          End If                  ! Calculate and save FED
       End Do DOOR_LOOP


       EXIT_LOOP: Do i = 1, N_EXITS
          !
          If (L_fed_save) Then
             If ( EVAC_EXITS(i)%FED_MESH > 0 ) Then
                i1 = EVAC_EXITS(i)%II
                j1 = EVAC_EXITS(i)%JJ
                k1 = EVAC_EXITS(i)%KK
                nom = EVAC_EXITS(i)%FED_MESH
                Z_1 = MESHES(nom)%YY(I1,J1,K1,I_FUEL)
                If (CO_PRODUCTION) Then
                   Z_2 = MESHES(nom)%YY(I1,J1,K1,I_PROG_CO)
                Else
                   Z_2 = 0._EB
                End If
                Z_3 = MESHES(nom)%YY(I1,J1,K1,I_PROG_F)
                y_extra = MESHES(nom)%Y_SUM(I1,J1,K1)  ! extra species mass fraction

                ! Mass fraction array ==> soot density (mg/m3)
                ! Next is for soot (mg/m3)
                Call GET_MASS_FRACTION2(Z_1,Z_2,Z_3,SOOT_INDEX,y_extra,Y_MF_INT)
                tmp_1 = Y_MF_INT*MESHES(nom)%RHO(i1,j1,k1)*1.E6_EB
                EVAC_EXITS(i)%SOOT_DENS = tmp_1

                ! Calculate Purser's fractional effective dose (FED)

                ! Next is for CO (ppm)
                Call GET_MASS_FRACTION2(Z_1,Z_2,Z_3,CO_INDEX,y_extra,Y_MF_INT)
                tmp_1 = RCON_MF(CO_INDEX)*Y_MF_INT*1.E6_EB/MESHES(nom)%RSUM(I1,J1,K1)
                EVAC_EXITS(i)%FED_CO_CO2_O2 = 3.317E-5_EB*25.0_EB* &
                     tmp_1**(1.036_EB)/(30.0_EB*60.0_EB)

                ! Next is for CO2
                Call GET_MASS_FRACTION2(Z_1,Z_2,Z_3,CO2_INDEX,y_extra,Y_MF_INT)
                tmp_1 = RCON_MF(CO2_INDEX)*Y_MF_INT/MESHES(nom)%RSUM(I1,J1,K1)
                EVAC_EXITS(i)%FED_CO_CO2_O2 = &
                     EVAC_EXITS(i)%FED_CO_CO2_O2*Exp( 0.1903_EB*tmp_1 &
                     *100.0_EB + 2.0004_EB )/7.1_EB

                ! Next is for O2
                Call GET_MASS_FRACTION2(Z_1,Z_2,Z_3,O2_INDEX,y_extra,Y_MF_INT)
                tmp_1 = RCON_MF(O2_INDEX)*Y_MF_INT/MESHES(nom)%RSUM(I1,J1,K1)
                If ( tmp_1 < 0.20_EB ) Then
                   EVAC_EXITS(i)%FED_CO_CO2_O2 = &
                        EVAC_EXITS(i)%FED_CO_CO2_O2 + 1.0_EB  / &
                        (60.0_EB*Exp(8.13_EB-0.54_EB*(20.9_EB-100.0_EB*tmp_1)) )
                End If

                ! Gas temperature, ind=5, C
                EVAC_EXITS(i)%TMP_G  = MESHES(nom)%TMP(i1,j1,k1)

                ! Radiant intensity, ind=18, kW/m2 (no -sigma*Tamb^4 term)
                EVAC_EXITS(i)%RADINT = &
                     Max(MESHES(nom)%UII(i1,j1,k1),4.0_EB*SIGMA*TMPA4)

             Else
                ! No fed_mesh found
                EVAC_EXITS(i)%FED_CO_CO2_O2 = 0.0_EB
                EVAC_EXITS(i)%SOOT_DENS = 0.0_EB
                EVAC_EXITS(i)%TMP_G = 0.0_EB
                EVAC_EXITS(i)%RADINT = 0.0_EB
             End If                ! fed_mesh > 0, i.e. fire grid found

             ! Save Fed, Soot, Temp(C), and RadInt
             Write (LU_EVACFED) &
                  Real(EVAC_EXITS(i)%FED_CO_CO2_O2,FB), &
                  Real(EVAC_EXITS(i)%SOOT_DENS,FB), &
                  Real(EVAC_EXITS(i)%TMP_G,FB), &
                  Real(EVAC_EXITS(i)%RADINT,FB)

          Else                    ! Read FED from a file
             ! Read Fed, Soot, Temp(C), and RadInt
             Read (LU_EVACFED,Iostat=ios) tmpout1, tmpout2, tmpout3, tmpout4
             If (ios.Ne.0) Then
                Write(MESSAGE,'(A)') 'ERROR: Evac Mesh Exchange: FED READ ERROR'
                Close (LU_EVACEFF)
                Call SHUTDOWN(MESSAGE)
             End If
             EVAC_EXITS(i)%FED_CO_CO2_O2 = tmpout1
             EVAC_EXITS(i)%SOOT_DENS = tmpout2
             EVAC_EXITS(i)%TMP_G = tmpout3
             EVAC_EXITS(i)%RADINT = tmpout4
          End If                  ! Calculate and save FED
       End Do EXIT_LOOP

    End If                    ! l_use_fed
    !    Goto 325

324 Continue
    If (ios < 0) Then 
       Write (LU_ERR,fmt='(a,f12.4,a)') 'FED file EOF: time ', &
            T_Save-DT_Save, ' not found'
       Write (LU_ERR,fmt='(a)') 'FED file EOF: use previous values'
       T_Save = 1.0E15
    End If

    !325 Continue

  End Subroutine EVAC_MESH_EXCHANGE
!
  Subroutine PREPARE_TO_EVACUATE(ICYC)
    Implicit None
    !
    ! Do the mesh independent initializations for the 
    ! subroutine evacuate_humans.
    !
    Integer, Intent(IN) :: ICYC
    !
    Logical L_eff_read, L_eff_save
    Integer(4) ibar_tmp, jbar_tmp, kbar_tmp
    Integer nm_tim, i, j
    ! 
    Type (MESH_TYPE), Pointer :: MFF

    If (.NOT.Any(EVACUATION_GRID)) RETURN

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
    !
    Real(EB), Intent(IN) :: Tin
    Integer, Intent(IN) :: NM,ICYC
    !
    Real(EB) DTSP,UBAR,VBAR, &
         X1,Y1,XI,YJ,ZK
    Integer ICN,I,J,IIN,JJN,KKN,II,JJ,KK,IIX,JJY,KKZ, &
         IBC, ICX, ICY
    Integer  IE, tim_ic, tim_iw, nm_tim, NM_now, j1, tim_iwx, tim_iwy
    Real(EB) P2P_DIST, P2P_DIST_MAX, P2P_U, P2P_V, &
         EVEL, tim_dist
    Integer istat
    !
    !
    Real(EB)  scal_prod_over_rsqr, t_init_timo, &
         U_new, V_new, Vmax_timo, CosPhiFac, &
         Speed_max, Delta_min, Dt_sum, &
         C_Yeff, t_flow_relax, LambdaW, &
         B_Wall, A_Wall, T, Contact_F, Social_F, &
         smoke_beta, smoke_alpha, smoke_speed_fac
    Integer iie, jje, iio, jjo, iii, jjj, i_egrid, i_tmp, i_o
    Real(EB) y_o, x_o, Delta_x, Delta_y, x_now, y_now, &
         d_humans, d_walls, &
         DTSP_new, &
         fac_tim, dt_group_door, x1_old, y1_old, x11, y11, &
         max_fed, L2_min, L2_tmp, Speed, Tpre, ave_K
    Logical PP_see_each, PP_see_door
    Logical L_eff_read, L_eff_save, L_Dead
    Real(EB) :: cos_x, cos_y, &
         speed_xm, speed_xp, speed_ym, speed_yp, hr_z, hr_a, hr_b, hr_tau
    !
    Real(EB) rn
    Real(EB) GaMe, GaTh, GaCM
    !
    Real(EB), Dimension(:), Allocatable :: K_ave_Door
    Real(EB), Dimension(:), Allocatable :: FED_max_Door
    Integer, Dimension(:), Allocatable :: Color_Tmp
    Logical, Dimension(:), Allocatable :: Is_Known_Door, &
         Is_Visible_Door
    Integer :: i_old_ffield, i_new_ffield, IEL, IZERO, color_index
    Character(26) :: name_old_ffield, name_new_ffield
    !
    !
    Real(EB) :: P2P_Torque, Fc_x, Fc_y, Omega_new, angle, A1, Tc_z, Fc_x1, Fc_y1
    Real(EB) :: Omega_max, Omega_0
    Real(EB), Dimension(6) :: y_tmp, x_tmp, r_tmp, v_tmp, u_tmp
    !
    Integer :: Max_Humans_Cell, i_dx, j_dy, ie_max, bl_max
    Real(EB) :: dx_min, dy_min, evac_dt_min2
    Integer, Dimension(:,:,:), Allocatable :: BLOCK_GRID
    Integer, Dimension(:,:), Allocatable :: BLOCK_GRID_N
    Integer, Dimension(:), Allocatable :: BLOCK_LIST

    !
    Real(EB) d_humans_min, d_walls_min
    Real(EB) TNOW, tnow13, tnow14, tnow15
    ! 
    Type (MESH_TYPE), Pointer :: MFF
    !
    TNOW=SECOND()
    If ( .Not.(EVACUATION_ONLY(NM) .And. EVACUATION_GRID(NM)) ) Return
!!$    EVAC_MESH_ONLY: If ( EVACUATION_ONLY(NM) .And. &
!!$         EVACUATION_GRID(NM) ) Then
    !
    t_flow_relax   = 0.0_EB
    t_init_timo    = 0.0_EB + t_flow_relax

    Vmax_timo      = V_MAX
    Omega_max      = V_ANGULAR_MAX*2.0_EB*Pi ! 8 rounds per second
    !Timo Omega_max      = 8.0_EB*2.0_EB*Pi  ! 8 rounds per second
    Omega_0        = V_ANGULAR*2.0_EB*Pi     ! 2 rounds per second
    Dt_sum         = 0.0_EB
    Delta_min      = Huge(Delta_min)
    dt_group_door  = 1.0_EB
    !
    Call POINT_TO_MESH(NM)
    !       Call POINT_TO_EVAC_MESH(NM)
    !
    dx_min = Minval(DX)
    dy_min = Minval(DY)
    Delta_min = Min( dy_min, dx_min, Delta_min )
    !
    L_eff_read = Btest(I_EVAC,2)
    L_eff_save = Btest(I_EVAC,0)
!!$       If ( ICYC == 0 .And. L_eff_save ) Then
!!$          Do nm_tim = 1, NMESHES
!!$             If (EVACUATION_ONLY(nm_tim)) Then
!!$                MFF=>MESHES(nm_tim)
!!$                ibar_tmp = MFF%IBAR
!!$                jbar_tmp = MFF%JBAR
!!$                kbar_tmp = 1
!!$                Write (LU_EVACEFF) ibar_tmp, jbar_tmp, kbar_tmp
!!$                Do  i = 0, MFF%IBAR+1
!!$                   Do j= 0, MFF%JBAR+1
!!$                      Write (LU_EVACEFF) Real(MFF%U(i,j,1),FB), &
!!$                           Real(MFF%V(i,j,1),FB)
!!$                   End Do
!!$                End Do
!!$             End If
!!$          End Do
!!$          ! Clear the save bit, save is done only once.
!!$          L_eff_save = .False.
!!$          I_EVAC = Ibclr(I_EVAC,0)
!!$       End If
!!$       !
!!$       ! Initialize counters only once for each time step.
!!$       If ( ICYC >= 0 .And. icyc_old < ICYC ) Then
!!$          icyc_old = ICYC
!!$          fed_max_alive = 0.0_EB
!!$          fed_max       = 0.0_EB
!!$       End If

    !
    ! Find the egrid index
    i_egrid = 0
    Do i = 1, n_egrids
       If (EVAC_Node_List(i)%Mesh_index == NM) Then
          i_egrid = EVAC_Node_List(i)%Node_Index
       End If
    End Do
    If (i_egrid == 0) Then
       Write(MESSAGE,'(A,I6)') &
            'ERROR: Evacuate_Humans, no mesh found ',NM
       Call SHUTDOWN(MESSAGE)
    End If
    !
    Allocate(Is_Known_Door(Max(1,n_doors+n_exits)),STAT=IZERO)
    Call ChkMemErr('Evacuate_Humans','Is_Known_Door',IZERO) 
    Allocate(Is_Visible_Door(Max(1,n_doors+n_exits)),STAT=IZERO)
    Call ChkMemErr('Evacuate_Humans','Is_Visible_Door',IZERO) 
    Allocate(FED_max_Door(Max(1,n_doors+n_exits)),STAT=IZERO)
    Call ChkMemErr('Evacuate_Humans','FED_max_Door',IZERO) 
    Allocate(K_ave_Door(Max(1,n_doors+n_exits)),STAT=IZERO)
    Call ChkMemErr('Evacuate_Humans','K_ave_Door',IZERO) 
    Allocate(Color_Tmp(Max(1,i33_dim)),STAT=IZERO)
    Call ChkMemErr('Evacuate_Humans','Color_Tmp',IZERO) 

    Allocate(BLOCK_GRID_N(1:IBAR,1:JBAR),STAT=IZERO)
    Call ChkMemErr('EVACUATE_HUMANS','BLOCK_GRID_N',IZERO)
    ! 
    HUMAN_TIME_LOOP: Do While ( Dt_sum < DT )
       DTSP = Min( (DT-Dt_sum), Tsteps(nm) )
       evac_dt_min2 = EVAC_DT_MAX
       d_humans_min = Huge(d_humans_min)
       d_walls_min = Huge(d_walls_min)

       ! ========================================================
       ! Add up the human time steps.
       ! ========================================================
       Dt_sum = Dt_sum + DTSP

       T = Tin - DT + Dt_sum


       ! ========================================================
       ! Evacuation routine, here the humans are moved to the new
       ! positions using the present velocities.  After this the
       ! forces at the new postitions are calculated and the 
       ! velocities are updated.  (The 'dissipative' self-driving
       ! force contribution to the velocities is updated self-
       ! consistently, but the friction terms are not.)
       ! ========================================================

       Group_List(:)%GROUP_SIZE  = 0
       Group_List(:)%GROUP_X = 0.0_EB
       Group_List(:)%GROUP_Y = 0.0_EB
       Group_List(:)%MAX_DIST_CENTER = 0.0_EB
       Group_List(:)%Speed   = 0.0_EB
       Group_List(:)%IntDose = 0.0_EB
       Group_List(:)%Tpre    = 0.0_EB
       Group_List(:)%Tdet    = Huge(Group_List(:)%Tdet)



       !
       Do i = 1, n_doors
          If ( EVAC_DOORS(i)%IMESH == nm) Then
             EVAC_DOORS(i)%NTARGET = 0
          End If
       End Do
       Do i = 1, n_exits
          If ( EVAC_EXITS(i)%IMESH == nm) Then
             EVAC_EXITS(i)%NTARGET = 0
          End If
       End Do

       Do j = 0, i33_dim
          Group_List(j)%GROUP_I_FFIELDS(i_egrid) = 0
       End Do
       Do i = 1, N_HUMANS
          HR=>HUMAN(I)
          j = Max(0,HR%GROUP_ID)
          Group_List(j)%GROUP_SIZE = Group_List(j)%GROUP_SIZE + 1
          Group_List(j)%GROUP_X    = Group_List(j)%GROUP_X + HR%X
          Group_List(j)%GROUP_Y    = Group_List(j)%GROUP_Y + HR%Y
          Group_List(j)%Speed      = Group_List(j)%Speed + HR%Speed
          Group_List(j)%IntDose    = Group_List(j)%IntDose + HR%IntDose
          Group_List(j)%Tpre       = Max(Group_List(j)%Tpre,HR%Tpre)
          Group_List(j)%Tdet       = Min(Group_List(j)%Tdet,HR%Tdet)
       End Do
       Group_List(1:)%GROUP_X = Group_List(1:)%GROUP_X / &
            Max(1,Group_List(1:)%GROUP_SIZE)
       Group_List(1:)%GROUP_Y = Group_List(1:)%GROUP_Y / &
            Max(1,Group_List(1:)%GROUP_SIZE)
       Group_List(1:)%Speed   = Group_List(1:)%Speed / &
            Max(1,Group_List(1:)%GROUP_SIZE)
       Group_List(1:)%IntDose = Group_List(1:)%IntDose / &
            Max(1,Group_List(1:)%GROUP_SIZE)

       Group_List(:)%MAX_DIST_CENTER = 0.0_EB
       Do i = 1, N_HUMANS
          HR=>HUMAN(I)
          ! group_id > 0: +group_id
          ! group_id < 0: -human_id (lonely humans)
          j  =  Max(0,HR%GROUP_ID)
          j1 = -Min(0,HR%GROUP_ID)
          Group_List(j)%MAX_DIST_CENTER = &
               Max(Group_List(j)%MAX_DIST_CENTER, &
               Sqrt((HR%X - Group_List(j)%GROUP_X)**2 + &
               (HR%Y - Group_List(j)%GROUP_Y)**2))

          i_tmp = HR%I_TARGET
          If (i_tmp > 0 .And. i_tmp <= n_doors ) Then
             EVAC_DOORS(i_tmp)%NTARGET = EVAC_DOORS(i_tmp)%NTARGET + 1
          End If
          i_tmp = i_tmp - n_doors
          If (i_tmp > 0 .And. i_tmp <= n_exits ) Then
             EVAC_EXITS(i_tmp)%NTARGET = EVAC_EXITS(i_tmp)%NTARGET + 1
          End If
       End Do

       ! j=0, i.e., lonely humans
       If (N_HUMANS > 0) Then
          Group_List(0)%GROUP_SIZE = 1
          Group_List(0)%GROUP_X    = 0.5_EB*(XS+XF)
          Group_List(0)%GROUP_Y    = 0.5_EB*(YS+YF)
          Group_List(0)%Speed      = 1.0_EB
          Group_List(0)%IntDose    = 0.0_EB
          Group_List(0)%MAX_DIST_CENTER = 0.0_EB
          Group_List(0)%COMPLETE = 1
       End If

       Do j = 1, i33_dim
          Group_List(j)%LIMIT_COMP = RADIUS_COMPLETE_0 + &
               RADIUS_COMPLETE_1*Group_List(j)%GROUP_SIZE
          If ( ((Group_List(j)%MAX_DIST_CENTER <=  &
               Group_List(j)%LIMIT_COMP) .Or. &
               (Group_List(j)%COMPLETE == 1)) .And. &
               Group_List(j)%GROUP_SIZE > 0) Then
             ! Note: If complete=1 already, it stays at 1.
             If (T > Group_List(j)%Tdet) Then
                If (Group_List(j)%COMPLETE == 0) Then
                   Group_List(j)%Tdoor = Max(T,Group_List(j)%Tdet)
                End If
                Group_List(j)%COMPLETE = 1
             End If
          End If
       End Do

       If (T > 0.0_EB) Then
          Change_Door_Loop: Do ie = 1, N_HUMANS
             HR => HUMAN(ie)
             i_old_ffield = HR%I_FFIELD
             name_old_ffield = Trim(MESH_NAME(i_old_ffield))
             j  =  Max(0,HR%GROUP_ID)
             j1 = -Min(0,HR%GROUP_ID)
             !cc          Do j = 1, i33_dim

             If (Group_List(j)%COMPLETE == 0) Cycle Change_Door_Loop
             If (j == 0 .And. T < HR%Tpre+HR%Tdet) &
                  Cycle Change_Door_Loop
             If (j > 0 .And. T < Group_List(j)%Tpre + &
                  Group_List(j)%Tdoor) Cycle Change_Door_Loop

             If (Group_List(j)%GROUP_I_FFIELDS(i_egrid) == 0) Then
                Call Random_number(rn)
                If (rn > Exp(-DTSP/dt_group_door) ) Then
                   i_old_ffield = HR%I_FFIELD
                   name_old_ffield = Trim(MESH_NAME(i_old_ffield))

                   If (HR%IEL > 0 ) Then
                      ! Human HR originates from an evac line
                      HPT => EVACUATION(HR%IEL)
                   Else
                      ! Human HR originates from an entr line
                      PNX => EVAC_ENTRYS(Abs(HR%IEL))
                   End If

                   K_ave_Door(:)      = 0.0_EB
                   FED_max_Door(:)    = 0.0_EB
                   Is_Known_Door(:)   = .False.
                   Is_Visible_Door(:) = .False.
                   If ( HR%GROUP_ID /= 0 ) Then
                      If ( HR%GROUP_ID < 0 ) Then
                         ! A lonely soul
                         x1_old = HR%X
                         y1_old = HR%Y
                         IEL    = HR%IEL
                         Speed  = HR%Speed
                         Do i = 1, n_doors
                            If ( EVAC_DOORS(i)%IMESH == nm) Then
                               Is_Visible_Door(i) = .True.
                               Do i_tmp = 1, Human_Known_Doors(j1)%N_nodes
                                  If (EVAC_DOORS(i)%INODE == &
                                       Human_Known_Doors(j1)%I_nodes(i_tmp)) &
                                       Is_Known_Door(i) = .True.
                               End Do
                            End If
                         End Do
                         Do i = 1, n_exits
                            If ( EVAC_EXITS(i)%IMESH == nm .And. &
                                 .Not. EVAC_EXITS(i)%COUNT_ONLY ) Then
                               Is_Visible_Door(n_doors+i) = .True.
                               Do i_tmp = 1, Human_Known_Doors(j1)%N_nodes
                                  If (EVAC_EXITS(i)%INODE == &
                                       Human_Known_Doors(j1)%I_nodes(i_tmp)) &
                                       Is_Known_Door(n_doors+i) = .True.
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
                            If ( EVAC_DOORS(i)%IMESH == nm) Then
                               Is_Visible_Door(i) = .True.
                               Do i_tmp = 1, Group_Known_Doors(j)%N_nodes
                                  If (EVAC_DOORS(i)%INODE == &
                                       Group_Known_Doors(j)%I_nodes(i_tmp)) &
                                       Is_Known_Door(i) = .True.
                               End Do
                            End If
                         End Do
                         Do i = 1, n_exits
                            If ( EVAC_EXITS(i)%IMESH == nm .And. &
                                 .Not. EVAC_EXITS(i)%COUNT_ONLY ) Then
                               Is_Visible_Door(n_doors+i) = .True.
                               Do i_tmp = 1, Group_Known_Doors(j)%N_nodes
                                  If (EVAC_EXITS(i)%INODE == &
                                       Group_Known_Doors(j)%I_nodes(i_tmp)) &
                                       Is_Known_Door(n_doors+i) = .True.
                               End Do
                            End If
                         End Do
                      End If
                   Else
                      x1_old = HR%X
                      y1_old = HR%Y
                      IEL    = Abs(HR%IEL)
                      Speed  = HR%Speed
                      Do i = 1, n_doors
                         If ( EVAC_DOORS(i)%IMESH == nm) Then
                            Is_Visible_Door(i) = .True.
                         End If
                      End Do
                      Do i = 1, n_exits
                         If ( EVAC_EXITS(i)%IMESH == nm .And. &
                              .Not. EVAC_EXITS(i)%COUNT_ONLY ) Then
                            Is_Visible_Door(n_doors+i) = .True.
                         End If
                      End Do
                      Do i = 1, PNX%N_VENT_FFIELDS 
                         i_tmp = 1
                         If (Trim(EVAC_Node_List(PNX%I_DOOR_NODES(i) &
                              )%Node_Type) == 'Door') Then
                            i_tmp = EVAC_Node_List(PNX%I_DOOR_NODES(i) &
                                 )%Node_Index 
                         End If
                         If (Trim(EVAC_Node_List(PNX%I_DOOR_NODES(i) &
                              )%Node_Type) == 'Exit' ) Then
                            i_tmp = n_doors + &
                                 EVAC_Node_List(PNX%I_DOOR_NODES(i) &
                                 )%Node_Index 
                         End If
                         If ( Is_Visible_Door(i_tmp) ) Then
                            ! Door/exit is on this floor
                            If (PNX%P_VENT_FFIELDS(i) < 0.5_EB) Then
                               Is_Known_Door(i_tmp) = .False.
                            Else
                               Is_Known_Door(i_tmp) = .True.
                            End If
                         End If
                      End Do
                   End If

                   ! Find the visible doors
                   Do i = 1, n_doors + n_exits
                      If ( Is_Visible_Door(i) ) Then
                         If (EVAC_Node_List(n_egrids+n_entrys+i)%Node_Type &
                              == 'Door' ) Then
                            X11 = EVAC_DOORS(i)%X 
                            Y11 = EVAC_DOORS(i)%Y
                         Else        ! 'Exit'
                            X11 = EVAC_EXITS(i-n_doors)%X 
                            Y11 = EVAC_EXITS(i-n_doors)%Y
                         End If
                         max_fed = 0.0_EB
                         ave_K   = 0.0_EB
                         II = Floor(CELLSI(Floor((X11-XS)*RDXINT)) + 1.0_EB  )
                         JJ = Floor(CELLSJ(Floor((Y11-YS)*RDYINT)) + 1.0_EB  )
                         KK = Floor(CELLSK(Floor(( 0.5_EB*(ZF+ZS)-ZS)*RDZINT)) + 1.0_EB  )
                         XI  = CELLSI(Floor((x1_old-XS)*RDXINT))
                         YJ  = CELLSJ(Floor((y1_old-YS)*RDYINT))
                         IIE = Floor(XI+1.0_EB)
                         JJE = Floor(YJ+1.0_EB)
                         PP_see_door = .True. ! oletusarvo
                         If (Abs(X1_OLD-X11) >= Abs(Y1_OLD-Y11)) Then
                            If ( iie < ii) Then
                               iio = iie
                               jjo = jje
                               iie = ii
                               jje = jj
                               y_o = Y1_OLD
                               x_o = X1_OLD
                               Delta_y = (Y11 - Y1_OLD)
                            Else
                               Delta_y = (Y1_OLD - Y11)
                               iio = ii
                               jjo = jj
                               y_o = Y11
                               x_o = X11
                            End If
                            Delta_x = Abs(X1_OLD - X11)
                            x_now = 0.0_EB
                            PP_see_door_x: Do iii = iio+1, iie-1
                               x_now = x_now + DX(iii)
                               y_now = y_o + x_now*(Delta_y/Delta_x)
                               jjj = Floor( CELLSJ(Floor((y_now-YS)*RDYINT)) &
                                    + 1.0_EB )
                               tim_ic = CELL_INDEX(iii,jjj,KK)
                               If (SOLID(tim_ic) .And. .Not. HR%I_Target == i) Then
                                  PP_see_door = .False.
                                  Exit PP_see_door_x
                               Else
                                  ave_K = ave_K + MASS_EXTINCTION_COEFFICIENT* &
                                       1.0E-6_EB*HUMAN_GRID(iii,jjj)%SOOT_DENS / &
                                       ( iie-1 - (iio+1) + 1)
                                  If (max_fed < &
                                       HUMAN_GRID(iii,jjj)%FED_CO_CO2_O2) Then
                                     max_fed = HUMAN_GRID(iii,jjj)%FED_CO_CO2_O2 
                                  End If
                               End If
                            End Do PP_see_door_x
                         Else 
                            If ( jje < jj) Then
                               iio = iie
                               jjo = jje
                               iie = ii
                               jje = jj
                               y_o = Y1_OLD
                               x_o = X1_OLD
                               Delta_x = (X11 - X1_OLD)
                            Else
                               Delta_x = (X1_OLD - X11)
                               iio = ii
                               jjo = jj
                               y_o = Y11
                               x_o = X11
                            End If
                            Delta_y = Abs(Y1_OLD - Y11)
                            y_now = 0.0_EB
                            PP_see_door_y: Do jjj = jjo+1, jje-1
                               y_now = y_now + DY(jjj)
                               x_now = x_o + y_now*(Delta_x/Delta_y)
                               iii = Floor( CELLSI(Floor((x_now-XS)*RDXINT)) &
                                    + 1.0_EB )
                               tim_ic = CELL_INDEX(iii,jjj,KK)
                               If (SOLID(tim_ic) .And. .Not. HR%I_Target == i) Then
                                  PP_see_door = .False.
                                  Exit PP_see_door_y
                               Else
                                  ave_K = ave_K + MASS_EXTINCTION_COEFFICIENT* &
                                       1.0E-6_EB*HUMAN_GRID(iii,jjj)%SOOT_DENS / &
                                       ( jje-1 - (jjo+1) + 1)
                                  If (max_fed < &
                                       HUMAN_GRID(iii,jjj)%FED_CO_CO2_O2) Then
                                     max_fed = HUMAN_GRID(iii,jjj)%FED_CO_CO2_O2 
                                  End If
                               End If
                            End Do PP_see_door_y
                         End If

                         If (PP_see_door) Then
                            If (EVAC_Node_List(n_egrids+n_entrys+i &
                                 )%Node_Type == 'Door') Then
                               If (.Not. EVAC_DOORS(i)%EXIT_SIGN .And. &
                                    .Not. HR%I_Target == i) Then
                                  Is_Visible_Door(i) = .False.
                               End If
                            End If
                            FED_max_Door(i) = max_fed
                            K_ave_Door(i) = ave_K 
                         Else
                            Is_Visible_Door(i) = .False.
                            iie = Floor(CELLSI(Floor((x1_old-XS)*RDXINT)) + 1.0_EB)
                            jje = Floor(CELLSJ(Floor((y1_old-YS)*RDYINT)) + 1.0_EB)
                            FED_max_Door(i) = &
                                 HUMAN_GRID(iie,jje)%FED_CO_CO2_O2 
                            K_ave_Door(i) = MASS_EXTINCTION_COEFFICIENT* &
                                 1.0E-6_EB*HUMAN_GRID(iie,jje)%SOOT_DENS
                         End If
                      End If ! correct floor
                   End Do ! doors and exits



                   If (Any(Is_Visible_Door)) Then
                      Do i = 1, n_doors + n_exits
                         If (Abs(HR%I_Target) == i .And. Is_Visible_Door(i)) &
                              Is_Known_Door(i) = .True.
                      End Do
                   End If


                   Do i = 1, n_doors
                      If ( EVAC_DOORS(i)%TIME_OPEN > T .Or. EVAC_DOORS(i)%TIME_CLOSE < T) Then
                         Is_Visible_Door(i) = .False.
                         Is_Known_Door(i) = .False.
                      End If
                   End Do
                   Do i = 1, n_exits
                      If ( (EVAC_EXITS(i)%TIME_OPEN > T .Or. EVAC_EXITS(i)%TIME_CLOSE < T) .And. &
                           .Not. EVAC_EXITS(i)%COUNT_ONLY ) Then
                         Is_Visible_Door(n_doors+i) = .False.
                         Is_Known_Door(n_doors+i) = .False.
                      End If
                   End Do

                   If (Any(Is_Known_Door) .Or. Any(Is_Visible_Door)) Then
                      i_tmp   = 0
                      L2_min = Huge(L2_min)
                      Do i = 1, n_doors + n_exits
                         If ( Is_Known_Door(i) .And. &
                              Is_Visible_Door(i) ) Then
                            x_o = 0.0_EB
                            y_o = 0.0_EB
                            If (Trim(EVAC_Node_List(n_egrids+n_entrys+i &
                                 )%Node_Type) == 'Door' ) Then
                               x_o = EVAC_DOORS( EVAC_Node_List( &
                                    i+n_egrids+n_entrys)%Node_Index )%X
                               y_o = EVAC_DOORS( EVAC_Node_List( &
                                    i+n_egrids+n_entrys)%Node_Index )%Y
                               i_o = EVAC_DOORS( EVAC_Node_List(i+n_egrids &
                                    + n_entrys)%Node_Index )%I_VENT_FFIELD
                            Else      ! 'Exit'
                               x_o = EVAC_EXITS( EVAC_Node_List( &
                                    i+n_egrids+n_entrys)%Node_Index )%X
                               y_o = EVAC_EXITS( EVAC_Node_List( &
                                    i+n_egrids+n_entrys)%Node_Index )%Y
                               i_o = EVAC_EXITS( EVAC_Node_List(i+n_egrids &
                                    + n_entrys)%Node_Index )%I_VENT_FFIELD
                            End If
                            If (FED_DOOR_CRIT > 0.0_EB) Then
                               L2_tmp = FED_max_Door(i) * Sqrt((x1_old-x_o)**2 &
                                    + (y1_old-y_o)**2)/Speed
                            Else
                               L2_tmp = K_ave_Door(i)
                            End If
                            If (i_o == i_old_ffield) L2_tmp = 0.1_EB*L2_tmp
                            If ( ((x_o-x1_old)**2 + (y_o-y1_old)**2) < &
                                 L2_min .And. L2_tmp < &
                                 Abs(FED_DOOR_CRIT)) Then
                               L2_min = (x_o-x1_old)**2 + (y_o-y1_old)**2 
                               L2_min = Max(0.0_EB,L2_min)
                               i_tmp = i
                            End If
                         End If
                      End Do
                      If (i_tmp > 0 ) Then
                         ! Known and visible door, no smoke
                         If (EVAC_Node_List(n_egrids+n_entrys+i_tmp &
                              )%Node_Type == 'Door' ) Then
                            name_new_ffield = &
                                 Trim(EVAC_DOORS(i_tmp)%VENT_FFIELD)
                            i_new_ffield = &
                                 EVAC_DOORS(i_tmp)%I_VENT_FFIELD
                         Else        ! 'Exit'
                            name_new_ffield = &
                                 Trim(EVAC_EXITS(i_tmp-n_doors)%VENT_FFIELD)
                            i_new_ffield = &
                                 EVAC_EXITS(i_tmp-n_doors)%I_VENT_FFIELD
                         End If
                         color_index = 0
                      Else
                         ! No visible known door available, try non-visible known doors
                         i_tmp   = 0
                         L2_min = Huge(L2_min)
                         Do i = 1, n_doors + n_exits
                            If ( Is_Known_Door(i) .And. &
                                 .Not. Is_Visible_Door(i) ) Then
                               x_o = 0.0_EB
                               y_o = 0.0_EB
                               If (EVAC_Node_List(n_egrids+n_entrys+i &
                                    )%Node_Type == 'Door' ) Then
                                  x_o = EVAC_DOORS( EVAC_Node_List( &
                                       i+n_egrids+n_entrys)%Node_Index )%X
                                  y_o = EVAC_DOORS( EVAC_Node_List( &
                                       i+n_egrids+n_entrys)%Node_Index )%Y
                                  i_o = EVAC_DOORS( EVAC_Node_List(i+n_egrids &
                                       + n_entrys)%Node_Index )%I_VENT_FFIELD
                               Else    ! 'Exit'
                                  x_o = EVAC_EXITS( EVAC_Node_List( &
                                       i+n_egrids+n_entrys)%Node_Index )%X
                                  y_o = EVAC_EXITS( EVAC_Node_List( &
                                       i+n_egrids+n_entrys)%Node_Index )%Y
                                  i_o = EVAC_EXITS( EVAC_Node_List(i+n_egrids &
                                       + n_entrys)%Node_Index )%I_VENT_FFIELD
                               End If
                               If (FED_DOOR_CRIT > 0.0_EB) Then
                                  L2_tmp = FED_max_Door(i)*Sqrt((x1_old-x_o)**2 &
                                       + (y1_old-y_o)**2)/Speed
                               Else
                                  l2_tmp = K_ave_Door(i)
                               End If
                               If (i_o == i_old_ffield) &
                                    L2_tmp = 0.1_EB*L2_tmp
                               If ( ((x_o-x1_old)**2 + (y_o-y1_old)**2) < &
                                    L2_min .And. L2_tmp < &
                                    Abs(FED_DOOR_CRIT)) Then
                                  L2_min = (x_o-x1_old)**2 + (y_o-y1_old)**2 
                                  L2_min = Max(0.0_EB,L2_min)
                                  i_tmp = i
                               End If
                            End If
                         End Do
                         If (i_tmp > 0 ) Then
                            !    Non-visible known door, no smoke
                            If (EVAC_Node_List( &
                                 n_egrids+n_entrys+i_tmp)%Node_Type &
                                 == 'Door' ) Then
                               name_new_ffield = &
                                    Trim(EVAC_DOORS(i_tmp)%VENT_FFIELD)
                               i_new_ffield = &
                                    EVAC_DOORS(i_tmp)%I_VENT_FFIELD
                            Else      ! 'Exit'
                               name_new_ffield = &
                                    Trim(EVAC_EXITS(i_tmp-n_doors)%VENT_FFIELD)
                               i_new_ffield = &
                                    EVAC_EXITS(i_tmp-n_doors)%I_VENT_FFIELD
                            End If
                            color_index = 1
                         Else
                            ! known doors with no smoke have not been found
                            i_tmp   = 0
                            L2_min = Huge(L2_min)
                            Do i = 1, n_doors + n_exits
                               If (Is_Visible_Door(i)) Then
                                  x_o = 0.0_EB
                                  y_o = 0.0_EB
                                  If (EVAC_Node_List(n_egrids+n_entrys+i &
                                       )%Node_Type == 'Door' ) Then
                                     x_o = EVAC_DOORS( EVAC_Node_List(i+ &
                                          n_egrids+n_entrys)%Node_Index )%X
                                     y_o = EVAC_DOORS( EVAC_Node_List(i+ &
                                          n_egrids+n_entrys)%Node_Index )%Y
                                     i_o = EVAC_DOORS( EVAC_Node_List(i+n_egrids &
                                          + n_entrys)%Node_Index )%I_VENT_FFIELD
                                  Else  ! 'Exit'
                                     x_o = EVAC_EXITS( EVAC_Node_List(i+ &
                                          n_egrids+n_entrys)%Node_Index )%X
                                     y_o = EVAC_EXITS( EVAC_Node_List(i+ &
                                          n_egrids+n_entrys)%Node_Index )%Y
                                     i_o = EVAC_EXITS( EVAC_Node_List(i+n_egrids &
                                          + n_entrys)%Node_Index )%I_VENT_FFIELD
                                  End If
                                  If (FED_DOOR_CRIT > 0.0_EB) Then
                                     L2_tmp = FED_max_Door(i)*Sqrt( &
                                          (x1_old-x_o)**2 + (y1_old-y_o)**2) &
                                          /Speed
                                  Else
                                     l2_tmp = K_ave_Door(i)
                                  End If
                                  If (i_o == i_old_ffield) &
                                       L2_tmp = 0.1_EB*L2_tmp
                                  If ( ((x_o-x1_old)**2 + (y_o-y1_old)**2) < &
                                       L2_min .And. L2_tmp < &
                                       Abs(FED_DOOR_CRIT)) Then
                                     L2_min = (x_o-x1_old)**2 + (y_o-y1_old)**2
                                     L2_min = Max(0.0_EB,L2_min)
                                     i_tmp = i
                                  End If
                               End If
                            End Do
                            If (i_tmp > 0 ) Then
                               ! No smoke, visible door (not known)
                               If (EVAC_Node_List( &
                                    n_egrids+n_entrys+i_tmp)%Node_Type &
                                    == 'Door' ) Then
                                  name_new_ffield = &
                                       Trim(EVAC_DOORS(i_tmp)%VENT_FFIELD)
                                  i_new_ffield = &
                                       EVAC_DOORS(i_tmp)%I_VENT_FFIELD
                               Else    ! 'Exit'
                                  name_new_ffield = Trim( &
                                       EVAC_EXITS(i_tmp-n_doors)%VENT_FFIELD)
                                  i_new_ffield = &
                                       EVAC_EXITS(i_tmp-n_doors)%I_VENT_FFIELD
                               End If
                               color_index = 2
                            Else
                               ! Now we have smoke and some visible or known doors
                               i_tmp   = 0
                               L2_min = Huge(L2_min)
                               Do i = 1, n_doors + n_exits
                                  If ( Is_Known_Door(i) .Or. &
                                       Is_Visible_Door(i) ) Then
                                     x_o = 0.0_EB
                                     y_o = 0.0_EB
                                     If (EVAC_Node_List(n_egrids+n_entrys+ &
                                          i)%Node_Type == 'Door' ) Then
                                        x_o = EVAC_DOORS( EVAC_Node_List(i+ &
                                             n_egrids+n_entrys)%Node_Index)%X
                                        y_o = EVAC_DOORS( EVAC_Node_List(i+ &
                                             n_egrids+n_entrys)%Node_Index)%Y
                                        i_o = EVAC_DOORS( EVAC_Node_List(i+ &
                                             n_egrids+ n_entrys)%Node_Index &
                                             )%I_VENT_FFIELD
                                     Else ! 'Exit'
                                        x_o = EVAC_EXITS( EVAC_Node_List(i+ &
                                             n_egrids+n_entrys)%Node_Index)%X
                                        y_o = EVAC_EXITS( EVAC_Node_List(i+ &
                                             n_egrids+n_entrys)%Node_Index)%Y
                                        i_o = EVAC_EXITS( EVAC_Node_List(i+ &
                                             n_egrids+ n_entrys)%Node_Index &
                                             )%I_VENT_FFIELD
                                     End If
                                     If (FED_DOOR_CRIT > 0.0_EB) Then
                                        If (j > 0) Then
                                           L2_tmp = Group_List(j)%IntDose + &
                                                FED_max_Door(i) * Sqrt( &
                                                (x1_old-x_o)**2+(y1_old-y_o)**2)/ &
                                                Speed
                                        Else
                                           L2_tmp = HR%IntDose + &
                                                FED_max_Door(i) * Sqrt( &
                                                (x1_old-x_o)**2+(y1_old-y_o)**2)/ &
                                                Speed
                                        End If
                                     Else
                                        l2_tmp = (Sqrt((x1_old-x_o)**2 &
                                             + (y1_old-y_o)**2)*0.5_EB) / &
                                             (3.0_EB/K_ave_Door(i))
                                     End If
                                     If (i_o == i_old_ffield) &
                                          L2_tmp = 0.9_EB*L2_tmp
                                     If (L2_tmp < L2_min) Then
                                        L2_min = L2_tmp
                                        i_tmp = i
                                     End If
                                  End If
                               End Do

                               If (i_tmp > 0 .And. L2_min < 1.0_EB) Then
                                  ! Not too much smoke, i.e., non-lethal amount (known or visible doors)
                                  If (EVAC_Node_List( &
                                       n_egrids+n_entrys+i_tmp)%Node_Type &
                                       == 'Door' ) Then
                                     name_new_ffield = &
                                          Trim(EVAC_DOORS(i_tmp)%VENT_FFIELD)
                                     i_new_ffield = &
                                          EVAC_DOORS(i_tmp)%I_VENT_FFIELD
                                  Else  ! 'Exit'
                                     name_new_ffield = Trim( &
                                          EVAC_EXITS(i_tmp-n_doors)%VENT_FFIELD)
                                     i_new_ffield = &
                                          EVAC_EXITS(i_tmp-n_doors)%I_VENT_FFIELD
                                  End If
                                  If (Is_Known_Door(i_tmp) .And. &
                                       Is_Visible_Door(i_tmp)) color_index = 3
                                  If (Is_Known_Door(i_tmp) .And. .Not. &
                                       Is_Visible_Door(i_tmp)) color_index = 4
                                  If (.Not. Is_Known_Door(i_tmp) .And. &
                                       Is_Visible_Door(i_tmp)) color_index = 5
                               Else    ! no match using cases 1-3
                                  ! No door found, use the main evac grid ffield
                                  i_tmp = 0
                                  name_new_ffield = Trim(MESH_NAME(nm))
                                  i_new_ffield    = nm
                                  color_index = 6
                               End If  ! case 4
                            End If    ! case 3
                         End If      ! case 2
                      End If        ! case 1
                      If (Color_Method .Eq. 4 ) Then
                         color_index =  6 ! default, cyan
                         If (i_tmp > 0 .And. i_tmp <= n_doors ) &
                              color_index = EVAC_DOORS(i_tmp)%COLOR_INDEX
                         If (i_tmp > n_doors .And. i_tmp <= n_doors + n_exits) &
                              color_index = EVAC_EXITS(i_tmp-n_doors)%COLOR_INDEX
                      End If

                   Else ! No known/visible door
                      i_tmp = 0 ! no door found
                      color_index = 6
                      If (HR%IEL > 0 ) Then
                         If (HPT%IMESH == nm) Then
                            i_new_ffield    = HPT%I_VENT_FFIELDS(0)
                            name_new_ffield = Trim(Mesh_Name(i_new_ffield))
                         Else
                            name_new_ffield = Trim(MESH_NAME(nm))
                            i_new_ffield    = nm
                         End If
                         If (Color_Method == 4) color_index = HPT%COLOR_INDEX
                      Else
                         If (PNX%IMESH == nm) Then
                            i_new_ffield    = PNX%I_VENT_FFIELDS(0)
                            name_new_ffield = Trim(Mesh_Name(i_new_ffield))
                         Else
                            name_new_ffield = Trim(MESH_NAME(nm))
                            i_new_ffield    = nm
                         End If
                         If (Color_Method == 4) color_index = PNX%COLOR_INDEX
                      End If
                   End If ! Any known door
                   HR%I_Target = i_tmp
                   If (i_tmp > 0 .And. .Not. Is_Visible_Door(Max(1,i_tmp))) Then
                      ! I_Target >0: visible, <0: not visible
                      HR%I_Target = -i_tmp
                   End If
                   If (j > 0) Group_Known_Doors(j)%I_Target = HR%I_Target
                   If (i_new_ffield /= i_old_ffield) Then
                      ! Change the field of this group.
                      If ( j == 0 ) Then
                         ! Group index=0: 'lost souls'
                         HR%I_FFIELD    = i_new_ffield
                         HR%FFIELD_NAME = Trim(name_new_ffield)
                         If (COLOR_METHOD == 5) HR%COLOR_INDEX = color_index
                         If (COLOR_METHOD == 4) HR%COLOR_INDEX = color_index
                         Write (LU_ERR,fmt='(a,i5,a,a,a,a)') &
                              ' EVAC: Human ',ie,', new ffield: ', &
                              Trim(name_new_ffield), ', old ffield: ', &
                              Trim(name_old_ffield)
                      Else
                         Group_List(j)%GROUP_I_FFIELDS(i_egrid) = &
                              i_new_ffield
                         HR%I_FFIELD    = i_new_ffield
                         HR%FFIELD_NAME = Trim(name_new_ffield)
                         If (COLOR_METHOD == 5) HR%COLOR_INDEX = color_index
                         If (COLOR_METHOD == 4) HR%COLOR_INDEX = color_index
                         If (j > 0) Color_Tmp(j) = color_index
                         Write (LU_ERR,fmt='(a,i5,a,a,a,a)') &
                              ' EVAC: Group ',j,', new ffield: ', &
                              Trim(name_new_ffield), ', old ffield: ', &
                              Trim(name_old_ffield)
                      End If

                   Else ! door is the same
                      If (COLOR_METHOD == 5) HR%COLOR_INDEX = color_index
                      If (COLOR_METHOD == 5 .And. j > 0) &
                           Color_Tmp(j) = HR%COLOR_INDEX
                      If (COLOR_METHOD == 4) HR%COLOR_INDEX = color_index
                      If (COLOR_METHOD == 4 .And. j > 0) &
                           Color_Tmp(j) = HR%COLOR_INDEX
                   End If

                Else ! Do not update door flow field at this time step
                   If (j > 0) Group_Known_Doors(j)%I_Target = HR%I_Target
                   If (COLOR_METHOD == 5 .And. j > 0) &
                        Color_Tmp(j) = HR%COLOR_INDEX
                   If (COLOR_METHOD == 4 .And. j > 0) &
                        Color_Tmp(j) = HR%COLOR_INDEX
                   If (j > 0) Group_List(j)%GROUP_I_FFIELDS(i_egrid) &
                        = HR%I_FFIELD 
                End If  ! update door flow field, rn large enough 
             Else
                ! This group has already tried to change the field
                If (COLOR_METHOD == 5 .And. j > 0) &
                     HR%COLOR_INDEX = Color_Tmp(j)
                If (COLOR_METHOD == 4 .And. j > 0) &
                     HR%COLOR_INDEX = Color_Tmp(j)
                HR%I_FFIELD    = Group_List(j)%GROUP_I_FFIELDS(i_egrid)
                HR%FFIELD_NAME = Trim(MESH_NAME(HR%I_FFIELD))
                If (j > 0) HR%I_Target = Group_Known_Doors(j)%I_Target
             End If  ! group_i_field=0
          End Do Change_Door_Loop  ! loop over humans, 
       End If  ! t > 0

       ! ========================================================
       ! MOVE LOOP STARTS HERE
       ! ========================================================
       TNOW15=SECOND()
       EVAC_MOVE_LOOP: Do I=1, N_HUMANS

          HR=>HUMAN(I)

          hr_z = 0.5_EB*(ZS+ZF)
          !
          LambdaW = LAMBDA_WALL
          A_Wall  = FAC_A_WALL*HR%A
          B_Wall  = FAC_B_WALL*HR%B
          GaMe    = NOISEME
          GaTh    = NOISETH
          GaCM    = NOISECM
          EVEL = Sqrt(HR%U**2 + HR%V**2)
          GaTh = Max( GaTh, &
               GaTh*(100.0_EB-110.0_EB*Abs(EVEL/HR%Speed)) )

          L_Dead  = .False.
          If ( HR%IntDose >= 1.0_EB  ) Then
             L_Dead = .True.
             ! No random force for a dead person.
             GaTh = 0.0_EB
             ! No psychological force terms for a dead person.
             A_Wall = 0.0_EB
             If (HR%Tpre /= Huge(HR%Tpre)) Then
                n_dead = n_dead+1
             End If
             HR%Tpre = Huge(HR%Tpre)
             HR%Tdet = Huge(HR%Tdet)
             HR%Tau  = HR%Tau
             HR%Mass = HR%Mass
             HR%COLOR_INDEX = 6
          Else
             fed_max_alive = Max(fed_max_alive,HR%IntDose)
          End If
          fed_max = Max(fed_max,HR%IntDose)
          hr_tau = HR%Tau
          !
          ! Check, if new random force is needed on next time step.
          If ( GaTh > 0.0_EB .And. T > 0.0_EB ) Then
             Call Random_number(rn)
             If ( rn > Exp(-DTSP/(0.2_EB*hr_tau)) ) HR%NewRnd = .True.
             If ( HR%NewRnd ) Then
                HR%NewRnd = .False.
                HR%ksi = (GaussRand(GaMe, GaTh, GaCM))
                Call Random_number(rn)
                HR%eta = 2.0_EB*Pi*rn
                HR%ksi = Abs(HR%ksi)
             End If
          End If
          !
          ! Where is the person
          !
          XI = CELLSI(Floor((HR%X-XS)*RDXINT))
          YJ = CELLSJ(Floor((HR%Y-YS)*RDYINT))
          ZK = CELLSK(Floor((HR_Z-ZS)*RDZINT))
          II  = Floor(XI+1.0_EB)
          JJ  = Floor(YJ+1.0_EB)
          KK  = Floor(ZK+1.0_EB)
          IIX = Floor(XI+0.5_EB)
          JJY = Floor(YJ+0.5_EB)
          KKZ = Floor(ZK+0.5_EB)
          HR%W = 0.0_EB

          ! Check the smoke density for Tdet
          If (HUMAN_GRID(ii,jj)%SOOT_DENS > TDET_SMOKE_DENS) Then
             HR%Tdet = Min(HR%Tdet,T)
          End If

          ! Calculate Purser's fractional effective dose (FED)
          smoke_speed_fac = 1.0_EB
          If (T > 0.0_EB) Then
             HR%IntDose = DTSP*HUMAN_GRID(II,JJ)%FED_CO_CO2_O2 + &
                  HR%IntDose
             smoke_beta  = -0.057_EB
             smoke_alpha = 0.706_EB
             smoke_speed_fac = 1.0_EB + (smoke_beta/smoke_alpha)* &
                  MASS_EXTINCTION_COEFFICIENT*1.0E-6_EB* &
                  HUMAN_GRID(II,JJ)%SOOT_DENS
             smoke_speed_fac = Max(smoke_speed_fac, 0.1_EB)
          End If
          HR%v0_fac = smoke_speed_fac


          ! ========================================================
          ! Calculate persons prefered walking direction
          ! ========================================================
          NM_now = NM
          MESH_ID_LOOP: Do nm_tim = 1, NMESHES
             If (MESH_NAME(nm_tim) == HR%FFIELD_NAME) Then
                NM_now = nm_tim
                Exit MESH_ID_LOOP
             End If
          End Do MESH_ID_LOOP
          ! 
          MFF=>MESHES(NM_now)
          UBAR = AFILL(MFF%U(II-1,JJY,KKZ),MFF%U(II,JJY,KKZ), &
               MFF%U(II-1,JJY+1,KKZ),MFF%U(II,JJY+1,KKZ), &
               MFF%U(II-1,JJY,KKZ+1),MFF%U(II,JJY,KKZ+1), &
               MFF%U(II-1,JJY+1,KKZ+1),MFF%U(II,JJY+1,KKZ+1), &
               XI-II+1,YJ-JJY+0.5_EB,ZK-KKZ+0.5_EB)
          VBAR = AFILL(MFF%V(IIX,JJ-1,KKZ),MFF%V(IIX+1,JJ-1,KKZ), &
               MFF%V(IIX,JJ,KKZ),MFF%V(IIX+1,JJ,KKZ), &
               MFF%V(IIX,JJ-1,KKZ+1),MFF%V(IIX+1,JJ-1,KKZ+1), &
               MFF%V(IIX,JJ,KKZ+1),MFF%V(IIX+1,JJ,KKZ+1), &
               XI-IIX+0.5_EB,YJ-JJ+1,ZK-KKZ+0.5_EB)

          EVEL = Sqrt(UBAR**2 + VBAR**2)
          If (EVEL > 0.0_EB) Then
             UBAR = UBAR/EVEL
             VBAR = VBAR/EVEL
          Else
             EVEL = Sqrt(HR%U**2 + HR%V**2)
             If (EVEL > 0.0_EB) Then
                UBAR = HR%U/EVEL
                VBAR = HR%V/EVEL
             End If
          End If
          ! ========================================================
          !
          ! Check if human is on a spectator stand.
          cos_x = 1.0_EB
          cos_y = 1.0_EB
          speed_xm = HR%Speed
          speed_xp = HR%Speed
          speed_ym = HR%Speed
          speed_yp = HR%Speed
          SS_Loop1: Do j = 1, n_sstands
             ESS => EVAC_SSTANDS(j)
             If (ESS%IMESH == nm .And. &
                  (ESS%X1 <= HR%X .And. ESS%X2 >= HR%X) .And. &
                  (ESS%Y1 <= HR%Y .And. ESS%Y2 >= HR%Y) ) Then
                cos_x = ESS%cos_x
                cos_y = ESS%cos_y
                Select Case (ESS%IOR)
                Case(-1)
                   speed_xm = cos_x*HR%Speed*ESS%FAC_V0_DOWN
                   speed_xp = cos_x*HR%Speed*ESS%FAC_V0_UP
                   speed_ym = HR%Speed*ESS%FAC_V0_HORI
                   speed_yp = HR%Speed*ESS%FAC_V0_HORI
                   HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + &
                        (ESS%H-ESS%H0)*Abs(ESS%X1-HR%X)/Abs(ESS%X1-ESS%X2)
                Case(+1)
                   speed_xm = cos_x*HR%Speed*ESS%FAC_V0_UP
                   speed_xp = cos_x*HR%Speed*ESS%FAC_V0_DOWN
                   speed_ym = HR%Speed*ESS%FAC_V0_HORI
                   speed_yp = HR%Speed*ESS%FAC_V0_HORI
                   HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + &
                        (ESS%H-ESS%H0)*Abs(ESS%X2-HR%X)/Abs(ESS%X1-ESS%X2)
                Case(-2)
                   speed_xm = HR%Speed*ESS%FAC_V0_HORI
                   speed_xp = HR%Speed*ESS%FAC_V0_HORI
                   speed_ym = cos_y*HR%Speed*ESS%FAC_V0_DOWN
                   speed_yp = cos_y*HR%Speed*ESS%FAC_V0_UP
                   HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + &
                        (ESS%H-ESS%H0)*Abs(ESS%Y1-HR%Y)/Abs(ESS%Y1-ESS%Y2)
                Case(+2)
                   speed_xm = HR%Speed*ESS%FAC_V0_HORI
                   speed_xp = HR%Speed*ESS%FAC_V0_HORI
                   speed_ym = cos_y*HR%Speed*ESS%FAC_V0_UP
                   speed_yp = cos_y*HR%Speed*ESS%FAC_V0_DOWN
                   HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + &
                        (ESS%H-ESS%H0)*Abs(ESS%Y2-HR%Y)/Abs(ESS%Y1-ESS%Y2)
                End Select
                Exit SS_Loop1
             End If
          End Do SS_Loop1
          !
          ! SC-VV: The new velocities are calculated using the old forces.
          ! Random forces and self-driving force are treated separately.

          U_new = HR%U + 0.5_EB*HR%F_X*DTSP/HR%Mass
          V_new = HR%V + 0.5_EB*HR%F_Y*DTSP/HR%Mass

          ! Rotational motion:
          Omega_new = HR%Omega + 0.5_EB*DTSP*HR%Torque/HR%M_iner

          ! Add random force term
          If ( GaTh > 0.0_EB .And. T > 0.0_EB ) Then
             U_new = U_new + 0.5_EB*DTSP*HR%v0_fac* &
                  HR%Mass*HR%ksi*Cos(HR%eta)/HR%Mass
             V_new = V_new + 0.5_EB*DTSP*HR%v0_fac* &
                  HR%Mass*HR%ksi*Sin(HR%eta)/HR%Mass
             Omega_new = Omega_new + 0.5_EB*DTSP* &
                  1.0_EB*Sign(HR%ksi,HR%eta-Pi)
          End If

          j = Max(0,HR%GROUP_ID)
          If (j == 0 ) Then
             Group_List(0)%Tpre = HR%Tpre + HR%Tdet
             Tpre = HR%Tpre + HR%Tdet
          Else
             Tpre = HR%Tdet
          End If

          If (Group_List(j)%GROUP_SIZE >= 2) Then
             HR%UBAR_Center = (Group_List(j)%GROUP_X - HR%X)
             HR%VBAR_Center = (Group_List(j)%GROUP_Y - HR%Y)
             EVEL = Sqrt(HR%UBAR_Center**2 + HR%VBAR_Center**2)
             If ( EVEL > 0.0_EB ) Then
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
             If (Group_List(j)%COMPLETE == 1 .And. &
                  T <= Group_List(j)%Tpre+Group_List(j)%Tdoor) Then
                UBAR = 0.0_EB
                VBAR = 0.0_EB
             End If
          End If

          EVEL = Sqrt(UBAR**2 + VBAR**2)
          If (Group_List(j)%COMPLETE == 1 .And. EVEL > 0.0_EB) Then
             UBAR = (1-GROUP_EFF)*UBAR + GROUP_EFF*HR%UBAR_Center/ &
                  Sqrt(((1-GROUP_EFF)*UBAR + GROUP_EFF*HR%UBAR_Center)**2+ &
                  ((1-GROUP_EFF)*VBAR + GROUP_EFF*HR%VBAR_Center)**2)
             VBAR = (1-GROUP_EFF)*VBAR + GROUP_EFF*HR%VBAR_Center/ &
                  Sqrt(((1-GROUP_EFF)*UBAR + GROUP_EFF*HR%UBAR_Center)**2+ &
                  ((1-GROUP_EFF)*VBAR + GROUP_EFF*HR%VBAR_Center)**2)
          Else
             UBAR = HR%UBAR_Center
             VBAR = HR%VBAR_Center
          End If

          If ( T <= Tpre ) Then
             ! Fire is not yet detected ==> no movement
             UBAR = 0.0_EB
             VBAR = 0.0_EB
             hr_tau = Max(0.1_EB,HR%Tau/10.0_EB)
          End If
          If ( T <= 0.0_EB ) Then
             ! Initialization phase
             UBAR = 0.0_EB
             VBAR = 0.0_EB
             hr_tau = HR%Tau
          End If

          EVEL = Sqrt(UBAR**2 + VBAR**2)
          If (EVEL > 0.0_EB) Then
             ! Sstands:  U,V are (x,y) plane projection velocities
             speed = speed_xp*(0.5_EB + Sign(0.5_EB,UBAR)) + &
                  speed_xm*(0.5_EB - Sign(0.5_EB,UBAR)) 
             speed = speed*HR%v0_fac
             U_new = U_new + 0.5_EB*(DTSP/hr_tau)* &
                  (speed*(UBAR/EVEL) - HR%U)
             HR%UBAR = speed*(UBAR/EVEL)

             speed = speed_yp*(0.5_EB + Sign(0.5_EB,VBAR)) + &
                  speed_ym*(0.5_EB - Sign(0.5_EB,VBAR)) 
             speed = speed*HR%v0_fac
             V_new = V_new + 0.5_EB*(DTSP/hr_tau)* &
                  (speed*(VBAR/EVEL) - HR%V)
             HR%VBAR = speed*(VBAR/EVEL)


             If (VBAR >= 0.0_EB) Then
                angle = ACos(UBAR/EVEL)
             Else
                angle = 2.0_EB*Pi - ACos(UBAR/EVEL)
             End If  ! angle is [0,2Pi)
             If (angle == 2.0_EB*Pi) angle = 0.0_EB

             ! J(dw/dt) = (J/t_iner)*( ((angle-angle_0/pi))w_0 - w )
             If (Abs(angle-HR%angle) <= Pi ) Then
                ! zero is not crossed.
                Omega_new = Omega_new + 0.5_EB*(DTSP/HR%tau_iner)* &
                     ( (angle-HR%angle)*(Omega_0/Pi) - HR%Omega)
             Else
                ! zero is crossed
                Omega_new = Omega_new + 0.5_EB*(DTSP/HR%tau_iner)* &
                     ( (2.0_EB*Pi-Abs(angle-HR%angle))*Sign(1.0_EB , HR%angle-angle)* &
                     (Omega_0/Pi) - HR%Omega )
             End If
          Else
             U_new = U_new + 0.5_EB*(DTSP/hr_tau)* &
                  (- HR%U)
             V_new = V_new + 0.5_EB*(DTSP/hr_tau)* &
                  (- HR%V)
             ! Slow rotation down if no direction available, i.e., target
             ! Omega_0 is zero.
             Omega_new = Omega_new + 0.5_EB*(DTSP/HR%tau_iner)*(-HR%Omega)
             HR%UBAR = 0.0_EB
             HR%VBAR = 0.0_EB
          End If

          ! Check that velocities are not too large, i.e.,
          ! unphysical (less than 10 m/s for humans)
          EVEL = Sqrt(U_new**2 + V_new**2)
          If ( EVEL > Vmax_timo ) Then
             U_new = U_new*(Vmax_timo/EVEL)
             V_new = V_new*(Vmax_timo/EVEL)
          End If
          ! Check that angular velocity is not too large
          If ( Abs(Omega_new) > Omega_max ) Then
             Omega_new = Sign(Omega_max,Omega_new)
          End If

          ! New coordinates
          X1 = HR%X + U_new*DTSP
          Y1 = HR%Y + V_new*DTSP
          HR%U = U_new
          HR%V = V_new
          ! New angle, should be on the interval [0,2Pi)
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
             If (ESS%IMESH == nm .And. &
                  (ESS%X1 <= HR%X .And. ESS%X2 >= HR%X) .And. &
                  (ESS%Y1 <= HR%Y .And. ESS%Y2 >= HR%Y) ) Then
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
          ! Where is the person, new coordinates (t + dt)?
          XI  = CELLSI(Floor((X1-XS)*RDXINT))
          YJ  = CELLSJ(Floor((Y1-YS)*RDYINT))
          ZK  = CELLSK(Floor((HR_Z-ZS)*RDZINT))
          IIN = Floor(XI+1.0_EB)
          JJN = Floor(YJ+1.0_EB)
          KKN = Floor(ZK+1.0_EB)
          ICN = CELL_INDEX(IIN,JJN,KKN)
          ICX = CELL_INDEX(IIN,JJ ,KKN)
          ICY = CELL_INDEX(II ,JJN,KKN)

          HR%X_old = HR%X
          HR%Y_old = HR%Y
          HR%Angle_old = HR%Angle

          ! Check, if moves inside a solid object ==> might be a open
          ! vent or a 'sucking vent' used to calculate the flow fields.
          If ( SOLID(ICN) ) Then
             If ( Solid(ICX) .And. .Not. Solid(ICY) ) Then
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
             Else If ( Solid(ICY) .And. .Not. Solid(ICX) ) Then
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
                Write(MESSAGE,'(A,I4,A,2F8.2)') &
                     'ERROR: Evacuate_Humans, Solid ICX and ICY, mesh ', &
                     nm, ' pos ',X1,Y1
             End If

          Else
             ! Target cell is not a solid ==> move
             HR%X = X1
             HR%Y = Y1
          End If
          HR%Angle = A1
          !
       End Do EVAC_MOVE_LOOP
       ! ========================================================
       ! MOVE LOOP ENDS HERE
       ! ========================================================
       ! 
       If (N_HUMANS > 0) Call CHECK_DOORS(T,NM)
       If (N_HUMANS > 0) Call CHECK_EXITS(T,NM)
       Call CHECK_CORRS(T,NM,DTSP)
       ! ========================================================
       ! Add persons (from doors and entrys)
       ! ========================================================
       If (T > 0.0_EB ) Then
          Do i = 1, N_ENTRYS
             Call ENTRY_HUMAN(i,T,NM,istat)
          End Do
       End If
       ! ========================================================
       ! Remove out-of-bounds persons (outside the grid)
       ! ========================================================
       If (N_HUMANS > 0) Call REMOVE_OUT_OF_GRIDS(T,NM)

       ! ========================================================
       ! Step (2) of SC-VV ends here
       ! ========================================================

       If ( ICYC >= 0) Then
          DTSP_new = EVAC_DT_STEADY_STATE
       Else
          DTSP_new = EVAC_DT_FLOWFIELD
       End If

       Speed_max  = 0.0_EB

       TUSED(15,NM)=TUSED(15,NM)+SECOND()-TNOW15


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
       Group_List(1:)%GROUP_X = Group_List(1:)%GROUP_X / &
            Max(1,Group_List(1:)%GROUP_SIZE)
       Group_List(1:)%GROUP_Y = Group_List(1:)%GROUP_Y / &
            Max(1,Group_List(1:)%GROUP_SIZE)

       BLOCK_GRID_N = 0
       Group_List(:)%MAX_DIST_CENTER = 0.0_EB
       Do i = 1, N_HUMANS
          HR => HUMAN(i)
          ! Which cell, new coordinates:
          IIN = Floor( CELLSI(Floor((HR%X-XS)*RDXINT)) + 1.0_EB )
          JJN = Floor( CELLSJ(Floor((HR%Y-YS)*RDYINT)) + 1.0_EB )
          BLOCK_GRID_N(IIN,JJN) = BLOCK_GRID_N(IIN,JJN) + 1
          j = Max(0,HR%GROUP_ID)
          !
          Group_List(j)%MAX_DIST_CENTER = &
               Max(Group_List(j)%MAX_DIST_CENTER, &
               Sqrt((HR%X - Group_List(j)%GROUP_X)**2 + &
               (HR%Y - Group_List(j)%GROUP_Y)**2))
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

       Do j = 1, i33_dim
          Group_List(j)%LIMIT_COMP = RADIUS_COMPLETE_0 + &
               RADIUS_COMPLETE_1*Group_List(j)%GROUP_SIZE
          If ( ((Group_List(j)%MAX_DIST_CENTER <=  &
               Group_List(j)%LIMIT_COMP) .Or. &
               (Group_List(j)%COMPLETE == 1)) .And. &
               Group_List(j)%GROUP_SIZE > 0 ) Then
             If (T > Group_List(j)%Tdet) Then
                If (Group_List(j)%COMPLETE == 0) Then
                   Group_List(j)%Tdoor = Max(T,Group_List(j)%Tdet)
                End If
                Group_List(j)%COMPLETE = 1
             End If
          End If
       End Do

       Allocate(BLOCK_GRID(1:IBAR,1:JBAR,Max_Humans_Cell),STAT=IZERO)
       Call ChkMemErr('EVACUATE_HUMANS','BLOCK_GRID',IZERO)
       BLOCK_GRID = 0

       BLOCK_GRID_N = 0
       Do i = 1, N_HUMANS
          HR => HUMAN(i)
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
       Allocate(BLOCK_LIST(bl_max),STAT=IZERO)
       Call ChkMemErr('EVACUATE_HUMANS','BLOCK_LIST',IZERO)
       ! ========================================================
       ! Step (3) of SC-VV starts here: Calculate new forces
       ! ========================================================
       ! Calculate forces on each person.
       TNOW13=SECOND()
       EVAC_FORCE_LOOP: Do I=1,N_HUMANS  
          HR => HUMAN(i)

          hr_z = 0.5_EB*(ZS+ZF)

          ! Calculate the position of the shoulder cirles
          y_tmp(1) = HR%Y - Cos(HR%angle)*HR%d_shoulder ! right
          x_tmp(1) = HR%X + Sin(HR%angle)*HR%d_shoulder
          y_tmp(2) = HR%Y ! torso
          x_tmp(2) = HR%X
          y_tmp(3) = HR%Y + Cos(HR%angle)*HR%d_shoulder ! left
          x_tmp(3) = HR%X - Sin(HR%angle)*HR%d_shoulder
          r_tmp(1) = HR%r_shoulder
          r_tmp(2) = HR%r_torso
          r_tmp(3) = HR%r_shoulder
          v_tmp(1) = HR%V + Sin(HR%angle)*HR%Omega*HR%d_shoulder
          v_tmp(2) = HR%V
          v_tmp(3) = HR%V - Sin(HR%angle)*HR%Omega*HR%d_shoulder
          u_tmp(1) = HR%U + Cos(HR%angle)*HR%Omega*HR%d_shoulder
          u_tmp(2) = HR%U
          u_tmp(3) = HR%U - Cos(HR%angle)*HR%Omega*HR%d_shoulder
          hr_a = HR%A
          hr_b = HR%B

          Contact_F = 0.0_EB
          Social_F  = 0.0_EB
          !
          LambdaW = LAMBDA_WALL
          A_Wall  = FAC_A_WALL*HR%A
          B_Wall  = FAC_B_WALL*HR%B
          GaMe    = NOISEME
          GaTh    = NOISETH
          GaCM    = NOISECM

          d_humans = Huge(d_humans)
          d_walls  = Huge(d_walls)
          L_Dead  = .False.
          If ( HR%IntDose >= 1.0_EB  ) Then
             L_Dead = .True.
             ! No random force for a dead person.
             GaTh = 0.0_EB
             ! No psychological force terms for a dead person.
             A_Wall = 0.0_EB
             HR%Tpre = Huge(HR%Tpre)
             HR%Tdet = Huge(HR%Tdet)
             HR%Tau  = HR%Tau
             HR%Mass = HR%Mass
             HR%COLOR_INDEX = 6
          End If
          hr_tau = HR%Tau
          !
          ! Psychological force: cut-off when acceleration below 0.0001 m/s**2
          P2P_DIST_MAX = HR%B*Log(HR%A/0.0001_EB)
          P2P_DIST_MAX = Min(P2P_DIST_MAX, 5.0_EB)
          ! If large pressure then short range forces only
          If ( HR%SumForces2 > 0.1_EB ) Then
             P2P_DIST_MAX = Min( P2P_DIST_MAX, &
                  -HR%B*Log(HR%SumForces2/(100.0_EB*HR%A)) )
          End If
          P2P_DIST_MAX = Max(P2P_DIST_MAX, 3.0_EB*HR%B)
          ! Speed up dead person loop, only contact forces are needed.
          If ( L_Dead ) P2P_DIST_MAX = 0.0_EB

          ! Where is the person, new coordinates:
          XI = CELLSI(Floor((HR%X-XS)*RDXINT))
          YJ = CELLSJ(Floor((HR%Y-YS)*RDYINT))
          ZK = CELLSK(Floor((HR_Z-ZS)*RDZINT))
          IIN  = Floor(XI+1.0_EB)
          JJN  = Floor(YJ+1.0_EB)
          KKN  = Floor(ZK+1.0_EB)
          IIX = Floor(XI+0.5_EB)
          JJY = Floor(YJ+0.5_EB)
          KKZ = Floor(ZK+0.5_EB)
          ICN = CELL_INDEX(IIN,JJN,KKN)
          X1 = HR%X 
          Y1 = HR%Y 
          HR%W = 0.0_EB
          !
          ! ========================================================
          ! Calculate persons prefered walking direction
          ! (Now the 'flow field' is used.)
          ! ========================================================
          NM_now = NM
          MESH_ID_LOOP2: Do nm_tim = 1, NMESHES
             If (MESH_NAME(nm_tim) == HR%FFIELD_NAME) Then
                NM_now = nm_tim
                Exit MESH_ID_LOOP2
             End If
          End Do MESH_ID_LOOP2
          ! 
          MFF=>MESHES(NM_now)
          UBAR = AFILL(MFF%U(IIN-1,JJY,KKZ),MFF%U(IIN,JJY,KKZ), &
               MFF%U(IIN-1,JJY+1,KKZ),MFF%U(IIN,JJY+1,KKZ), &
               MFF%U(IIN-1,JJY,KKZ+1),MFF%U(IIN,JJY,KKZ+1), &
               MFF%U(IIN-1,JJY+1,KKZ+1),MFF%U(IIN,JJY+1,KKZ+1), &
               XI-IIN+1,YJ-JJY+0.5_EB,ZK-KKZ+0.5_EB)
          VBAR = AFILL(MFF%V(IIX,JJN-1,KKZ),MFF%V(IIX+1,JJN-1,KKZ), &
               MFF%V(IIX,JJN,KKZ),MFF%V(IIX+1,JJN,KKZ), &
               MFF%V(IIX,JJN-1,KKZ+1),MFF%V(IIX+1,JJN-1,KKZ+1), &
               MFF%V(IIX,JJN,KKZ+1),MFF%V(IIX+1,JJN,KKZ+1), &
               XI-IIX+0.5_EB,YJ-JJN+1,ZK-KKZ+0.5_EB)

          EVEL = Sqrt(UBAR**2 + VBAR**2)
          If (EVEL > 0.0_EB) Then
             UBAR = UBAR/EVEL
             VBAR = VBAR/EVEL
          Else
             EVEL = Sqrt(HR%U**2 + HR%V**2)
             If (EVEL > 0.0_EB) Then
                UBAR = HR%U/EVEL
                VBAR = HR%V/EVEL
             End If
          End If

          ! ========================================================
          ! Update the Block_Grid array search ranges
          BLOCK_LIST = 0
          i_dx = Int((2.0_EB*0.3_EB+P2P_DIST_MAX)/DX(IIN)) + 1
          j_dy = Int((2.0_EB*0.3_EB+P2P_DIST_MAX)/DY(JJN)) + 1
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
             Write(MESSAGE,'(A,2I6)') &
                  'ERROR: Evacuate_Humans, ie_max, bl_max ', ie_max, bl_max
             Call SHUTDOWN(MESSAGE)
          End If

          ! ========================================================
          ! Check if human is on a spectator stand.
          cos_x = 1.0_EB
          cos_y = 1.0_EB
          speed_xm = HR%Speed
          speed_xp = HR%Speed
          speed_ym = HR%Speed
          speed_yp = HR%Speed
          SS_Loop2: Do j = 1, n_sstands
             ESS => EVAC_SSTANDS(j)
             If (ESS%IMESH == nm .And. &
                  (ESS%X1 <= HR%X .And. ESS%X2 >= HR%X) .And. &
                  (ESS%Y1 <= HR%Y .And. ESS%Y2 >= HR%Y) ) Then
                cos_x = ESS%cos_x
                cos_y = ESS%cos_y
                Select Case (ESS%IOR)
                Case(-1)
                   speed_xm = cos_x*HR%Speed*ESS%FAC_V0_DOWN
                   speed_xp = cos_x*HR%Speed*ESS%FAC_V0_UP
                   speed_ym = HR%Speed*ESS%FAC_V0_HORI
                   speed_yp = HR%Speed*ESS%FAC_V0_HORI
                   HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + &
                        (ESS%H-ESS%H0)*Abs(ESS%X1-HR%X)/Abs(ESS%X1-ESS%X2)
                Case(+1)
                   speed_xm = cos_x*HR%Speed*ESS%FAC_V0_UP
                   speed_xp = cos_x*HR%Speed*ESS%FAC_V0_DOWN
                   speed_ym = HR%Speed*ESS%FAC_V0_HORI
                   speed_yp = HR%Speed*ESS%FAC_V0_HORI
                   HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + &
                        (ESS%H-ESS%H0)*Abs(ESS%X2-HR%X)/Abs(ESS%X1-ESS%X2)
                Case(-2)
                   speed_xm = HR%Speed*ESS%FAC_V0_HORI
                   speed_xp = HR%Speed*ESS%FAC_V0_HORI
                   speed_ym = cos_y*HR%Speed*ESS%FAC_V0_DOWN
                   speed_yp = cos_y*HR%Speed*ESS%FAC_V0_UP
                   HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + &
                        (ESS%H-ESS%H0)*Abs(ESS%Y1-HR%Y)/Abs(ESS%Y1-ESS%Y2)
                Case(+2)
                   speed_xm = HR%Speed*ESS%FAC_V0_HORI
                   speed_xp = HR%Speed*ESS%FAC_V0_HORI
                   speed_ym = cos_y*HR%Speed*ESS%FAC_V0_UP
                   speed_yp = cos_y*HR%Speed*ESS%FAC_V0_DOWN
                   HR%Z = 0.5_EB*(ESS%Z1+ESS%Z2) + ESS%H0 + &
                        (ESS%H-ESS%H0)*Abs(ESS%Y2-HR%Y)/Abs(ESS%Y1-ESS%Y2)
                End Select
                Exit SS_Loop2
             End If
          End Do SS_Loop2

          hr_a =  HR%A*Max(0.5_EB,(Sqrt(HR%U**2+HR%V**2)/HR%Speed))
          A_Wall = Min(A_Wall, FAC_A_WALL*hr_a)

          ! ========================================================
          ! PERSON-PERSON INTERACTION FORCES
          ! ========================================================
          ! Look for other persons
          P2P_U      = 0.0_EB
          P2P_V      = 0.0_EB
          P2P_Torque = 0.0_EB
          TNOW14=SECOND()   ! person-person force loop timing
          P2PLOOP: Do IE = 1, ie_max
             If (BLOCK_LIST(ie) == I) Cycle P2PLOOP
             HRE => HUMAN(BLOCK_LIST(ie))

             P2P_DIST = ((HRE%X-X1)**2 + (HRE%Y-Y1)**2)
             If ( P2P_DIST > (P2P_DIST_MAX+HR%Radius+HRE%Radius)**2 ) Cycle P2PLOOP
             P2P_DIST = Sqrt(P2P_DIST)
             !
             ! Check, that the persons are seeing each other, i.e., there
             ! are no walls between.
             ! Where is the other person?
             XI  = CELLSI(Floor((HRE%X-XS)*RDXINT))
             YJ  = CELLSJ(Floor((HRE%Y-YS)*RDYINT))
             IIE = Floor(XI+1.0_EB)
             JJE = Floor(YJ+1.0_EB)
             PP_see_each = .True.
             If (Abs(HRE%X-X1) >= Abs(HRE%Y-Y1)) Then
                If ( iie < iin) Then
                   iio = iie
                   jjo = jje
                   iie = iin
                   jje = jjn
                   y_o = HRE%Y
                   x_o = HRE%X
                   Delta_y = (Y1 - HRE%Y)
                Else
                   Delta_y = (HRE%Y - Y1)
                   iio = iin
                   jjo = jjn
                   y_o = Y1
                   x_o = X1
                End If
                Delta_x = Abs(HRE%X - X1)
                x_now = 0.0_EB
                PP_see_loop_x: Do iii = iio+1, iie-1
                   x_now = x_now + DX(iii)
                   y_now = y_o + x_now*(Delta_y/Delta_x)
                   jjj = Floor(CELLSJ(Floor((y_now-YS)*RDYINT)) + 1.0_EB)
                   tim_ic = CELL_INDEX(iii,jjj,KKN)
                   If (SOLID(tim_ic)) Then
                      PP_see_each = .False.
                      Exit PP_see_loop_x
                   End If
                End Do PP_see_loop_x
             Else 
                If ( jje < jjn) Then
                   iio = iie
                   jjo = jje
                   iie = iin
                   jje = jjn
                   y_o = HRE%Y
                   x_o = HRE%X
                   Delta_x = (X1 - HRE%X)
                Else
                   Delta_x = (HRE%X - X1)
                   iio = iin
                   jjo = jjn
                   y_o = Y1
                   x_o = X1
                End If
                Delta_y = Abs(HRE%Y - Y1)
                y_now = 0.0_EB
                PP_see_loop_y: Do jjj = jjo+1, jje-1
                   y_now = y_now + DY(jjj)
                   x_now = x_o + y_now*(Delta_x/Delta_y)
                   iii = Floor(CELLSI(Floor((x_now-XS)*RDXINT)) + 1.0_EB)
                   tim_ic = CELL_INDEX(iii,jjj,KKN)
                   If (SOLID(tim_ic)) Then
                      PP_see_each = .False.
                      Exit PP_see_loop_y
                   End If
                End Do PP_see_loop_y
             End If
             If (.Not. PP_see_each) Cycle P2PLOOP
             !
             ! 
             ! Calculate the combination of spring constant for humans
             C_Yeff = (2.0_EB*HR%C_Young*2.0_EB*HRE%C_Young)/(2.0_EB*HR%C_Young+2.0_EB*HRE%C_Young)
             !
             ! Directional force:
             If ( (HR%U**2 +HR%V**2) > 0.0_EB ) Then
                CosPhiFac = ( (HRE%X-X1)*HR%U + (HRE%Y-Y1)*HR%V ) &
                     / ( Sqrt((HRE%X-X1)**2 + (HRE%Y-Y1)**2)* &
                     Sqrt(HR%U**2 +HR%V**2) )
                CosPhiFac = HR%Lambda + &
                     0.5_EB*(1.0_EB-HR%Lambda)*(1.0_EB+CosPhiFac)
             Else
                CosPhiFac = 1.0_EB
             End If
             ! Calculate the position of the shoulder cirles
             r_tmp(4) = HRE%r_shoulder ! right circle
             r_tmp(5) = HRE%r_torso     ! center circle
             r_tmp(6) = HRE%r_shoulder ! left circle
             y_tmp(4) = HRE%Y - Cos(HRE%angle)*HRE%d_shoulder ! right circle
             x_tmp(4) = HRE%X + Sin(HRE%angle)*HRE%d_shoulder
             y_tmp(5) = HRE%Y ! center circle
             x_tmp(5) = HRE%X
             y_tmp(6) = HRE%Y + Cos(HRE%angle)*HRE%d_shoulder ! left circle
             x_tmp(6) = HRE%X - Sin(HRE%angle)*HRE%d_shoulder
             ! ========================================================
             ! Add psychological force term
             ! ========================================================
             If (.Not. L_Dead) Then
                Fc_x = 0.0_EB
                Fc_y = 0.0_EB
                Do iii = 1, 3
                   Do jjj = 4, 6
                      tim_dist = Sqrt((x_tmp(iii)-x_tmp(jjj))**2 + &
                           (y_tmp(iii)-y_tmp(jjj))**2)
                      d_humans = Min( tim_dist-(r_tmp(iii)+r_tmp(jjj)) , d_humans )
                      Fc_x1 = (x_tmp(iii)-x_tmp(jjj)) * &
                           hr_a*CosPhiFac*Exp( -(tim_dist- &
                           ( r_tmp(iii)+r_tmp(jjj) ))/hr_b )/tim_dist 
                      Fc_y1 = (y_tmp(iii)-y_tmp(jjj)) * &
                           hr_a*CosPhiFac*Exp( -(tim_dist- &
                           ( r_tmp(iii)+r_tmp(jjj) ))/hr_b )/tim_dist 
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
                Do jjj = 4, 6
                   ! First the right shoulder
                   tim_dist = Sqrt( (x_tmp(jjj)-x_tmp(1))**2 + (y_tmp(jjj)-y_tmp(1))**2 )
                   Fc_x = (x_tmp(1)-x_tmp(jjj)) * &
                        hr_a*CosPhiFac*Exp( -(tim_dist - &
                        (r_tmp(1)+r_tmp(jjj)))/hr_b )/tim_dist
                   Fc_y = (y_tmp(1)-y_tmp(jjj)) * &
                        hr_a*CosPhiFac*Exp( -(tim_dist- &
                        (r_tmp(1)+r_tmp(jjj)))/hr_b )/tim_dist
                   If ( Abs(Fc_y*(x_tmp(1)-HR%X) - Fc_x*(y_tmp(1)-HR%Y)) > Abs(Tc_z) ) Then
                      Tc_z = Fc_y*(x_tmp(1)-HR%X) - Fc_x*(y_tmp(1)-HR%Y)
                   End If
                End Do
                P2P_Torque = P2P_Torque + Tc_z
                Tc_z = 0.0_EB
                Do jjj = 4, 6
                   ! Then the left shoulder
                   tim_dist = Sqrt( (x_tmp(jjj)-x_tmp(3))**2 + (y_tmp(jjj)-y_tmp(3))**2 )
                   Fc_x = (x_tmp(3)-x_tmp(jjj)) * &
                        hr_a*CosPhiFac*Exp( -(tim_dist - &
                        (r_tmp(3)+r_tmp(jjj)))/hr_b )/tim_dist
                   Fc_y = (y_tmp(3)-y_tmp(jjj)) * &
                        hr_a*CosPhiFac*Exp( -(tim_dist- &
                        (r_tmp(3)+r_tmp(jjj)))/hr_b )/tim_dist
                   If ( Abs(Fc_y*(x_tmp(3)-HR%X) - Fc_x*(y_tmp(3)-HR%Y)) > Abs(Tc_z) ) Then
                      Tc_z = Fc_y*(x_tmp(3)-HR%X) - Fc_x*(y_tmp(3)-HR%Y)
                   End If
                End Do
                P2P_Torque = P2P_Torque + Tc_z
             End If
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
                      tim_dist = Sqrt((x_tmp(iii)-x_tmp(jjj))**2 + &
                           (y_tmp(iii)-y_tmp(jjj))**2)
                      d_humans = Min( tim_dist-(r_tmp(iii)+r_tmp(jjj)) , d_humans )
                      If (tim_dist <= r_tmp(iii)+r_tmp(jjj) ) Then
                         ! Circles are touching each others
                         Fc_x =(x_tmp(iii)-x_tmp(jjj)) * &
                              C_Yeff*((r_tmp(iii)+r_tmp(jjj))-tim_dist)/tim_dist
                         Fc_y =(y_tmp(iii)-y_tmp(jjj)) * &
                              C_Yeff*((r_tmp(iii)+r_tmp(jjj))-tim_dist)/tim_dist

                         Fc_x = Fc_x - FC_DAMPING*(u_tmp(iii)-u_tmp(jjj))*&
                              (x_tmp(iii)-x_tmp(jjj))/tim_dist
                         Fc_y = Fc_y - FC_DAMPING*(v_tmp(iii)-v_tmp(jjj))*&
                              (y_tmp(iii)-y_tmp(jjj))/tim_dist

                         Contact_F = Contact_F + Sqrt(Fc_x**2 + Fc_y**2)
                         P2P_U = P2P_U + Fc_x
                         P2P_V = P2P_V + Fc_y
                         P2P_Torque = P2P_Torque + Fc_y*(x_tmp(iii)-HR%X) - Fc_x*(y_tmp(iii)-HR%Y)

                         scal_prod_over_rsqr = ((y_tmp(iii)-y_tmp(jjj))*(u_tmp(iii)-u_tmp(jjj)) - &
                              (x_tmp(iii)-x_tmp(jjj))*(v_tmp(iii)-v_tmp(jjj))) / (tim_dist**2)
                         If (I_Fric_sw >= 1 ) Then
                            Fc_x = - HR%Kappa * &
                                 ((r_tmp(iii)+r_tmp(jjj))-tim_dist)* &
                                 ( (y_tmp(iii)-y_tmp(jjj)) * scal_prod_over_rsqr )
                            Fc_y = - HR%Kappa * &
                                 ((r_tmp(iii)+r_tmp(jjj))-tim_dist)* &
                                 (-(x_tmp(iii)-x_tmp(jjj)) * scal_prod_over_rsqr )
                            P2P_U = P2P_U + Fc_x
                            P2P_V = P2P_V + Fc_y

                            P2P_Torque = P2P_Torque + Fc_y*( (x_tmp(iii) + &
                                 (r_tmp(iii)/r_tmp(jjj))*(x_tmp(jjj)-x_tmp(iii)) ) - HR%X ) 
                            P2P_Torque = P2P_Torque - Fc_x*( (y_tmp(iii) + &
                                 (r_tmp(iii)/r_tmp(jjj))*(y_tmp(jjj)-y_tmp(iii)) ) - HR%Y ) 
                         Else
                            Fc_x = - HR%Gamma * ( (y_tmp(iii)-y_tmp(jjj)) &
                                 * scal_prod_over_rsqr )
                            Fc_y = - HR%Gamma * (-(x_tmp(iii)-x_tmp(jjj)) &
                                 * scal_prod_over_rsqr )
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
             End If
          End Do P2PLOOP
          TUSED(14,NM)=TUSED(14,NM)+SECOND()-TNOW14
          ! ========================================================
          ! PERSON-PERSON INTERACTION FORCES ENDS HERE
          ! ========================================================

          ! ========================================================
          ! The wall forces
          ! ========================================================
          !

          ! ========================================================
          ! WALL SURFACE - PERSON FORCES STARTS
          ! ========================================================
          tim_dist = Huge(tim_dist)
          ! Find the closest wall at the -x direction
          x_now = -DX(iin)
          TIM_MX: Do ii = iin,0,-1
             x_now = x_now + DX(ii)
             If ( x_now-HR%Radius > P2P_DIST_MAX) Exit TIM_MX
             tim_ic = cell_index(ii,jjn,kkn)
             Call Get_iw(ii,jjn,kkn,+1,tim_iw)
             If ( tim_iw == 0 ) Then
                Cycle TIM_MX
             Else
                IBC = IJKW(5,tim_iw)
                If (SURFACE(IBC)%VEL> 0.0_EB .Or. BOUNDARY_TYPE(tim_iw)==OPEN_BOUNDARY) Exit TIM_MX
                tim_dist = x1 - xw(tim_iw)
                If ( (HR%U**2 +HR%V**2) > 0.0_EB ) Then
                   CosPhiFac = (-HR%U)/Sqrt(HR%U**2+HR%V**2)
                   CosPhiFac = LambdaW +  &
                        0.5_EB*(1.0_EB-LambdaW)*(1.0_EB+CosPhiFac)
                Else
                   CosPhiFac = 1.0_EB
                End If
                P2P_U = P2P_U + A_Wall*CosPhiFac* &
                     Exp( -(tim_dist-HR%Radius)/B_Wall )
                Social_F = Social_F + Abs(A_Wall*CosPhiFac* &
                     Exp( -(tim_dist-HR%Radius)/B_Wall ))
                Do iii = 1, 3
                   tim_dist = x_tmp(iii) - xw(tim_iw)
                   d_walls = Min(tim_dist-r_tmp(iii),d_walls)
                   If (tim_dist <= r_tmp(iii)) Then
                      Fc_x = + 2.0_EB*HR%C_Young*(r_tmp(iii)-tim_dist)
                      Fc_y = 0.0_EB
                      Contact_F = Contact_F + Abs(Fc_x)
                      If (I_Fric_sw >= 1 ) Then
                         Fc_y = Fc_y - HR%Kappa*(r_tmp(iii)-tim_dist)*v_tmp(iii)
                      Else
                         Fc_y = Fc_y - HR%Gamma*v_tmp(iii)
                      End If
                      P2P_Torque = P2P_Torque + Fc_y*(xw(tim_iw)-HR%X) - &
                           Fc_x*(y_tmp(iii)-HR%Y)
                      P2P_U = P2P_U + Fc_x
                      P2P_V = P2P_V + Fc_y
                   End If
                End Do

                Exit TIM_MX
             End If
          End Do TIM_MX

          ! Find the closest wall at the +x direction
          x_now = -DX(iin)
          TIM_PX: Do ii = iin, IBAR+1
             x_now = x_now + DX(ii)
             If ( x_now-HR%Radius > P2P_DIST_MAX) Exit TIM_PX
             tim_ic = cell_index(ii,jjn,kkn)
             Call Get_iw(ii,jjn,kkn,-1,tim_iw)
             If ( tim_iw == 0 ) Then
                Cycle TIM_PX
             Else
                IBC = IJKW(5,tim_iw)
                If (SURFACE(IBC)%VEL> 0.0_EB .Or. BOUNDARY_TYPE(tim_iw)==OPEN_BOUNDARY) Exit TIM_PX
                tim_dist = xw(tim_iw) - x1
                If ( (HR%U**2 +HR%V**2) > 0.0_EB ) Then
                   CosPhiFac = (HR%U)/Sqrt(HR%U**2+HR%V**2)
                   CosPhiFac = LambdaW +  &
                        0.5_EB*(1.0_EB-LambdaW)*(1.0_EB+CosPhiFac)
                Else
                   CosPhiFac = 1.0_EB
                End If
                P2P_U = P2P_U - A_Wall*CosPhiFac* &
                     Exp( -(tim_dist-HR%Radius)/B_Wall )
                Social_F = Social_F + Abs(A_Wall*CosPhiFac* &
                     Exp( -(tim_dist-HR%Radius)/B_Wall ))
                Do iii = 1, 3
                   tim_dist = xw(tim_iw) - x_tmp(iii)
                   d_walls = Min(tim_dist-r_tmp(iii),d_walls)
                   If (tim_dist <= r_tmp(iii)) Then
                      Fc_x = - 2.0_EB*HR%C_Young*(r_tmp(iii)-tim_dist)
                      Fc_y = 0.0_EB
                      Contact_F = Contact_F + Abs(Fc_x)
                      If (I_Fric_sw >= 1 ) Then
                         Fc_y = Fc_y - HR%Kappa*(r_tmp(iii)-tim_dist)*v_tmp(iii)
                      Else
                         Fc_y = Fc_y - HR%Gamma*v_tmp(iii)
                      End If
                      P2P_Torque = P2P_Torque + Fc_y*(xw(tim_iw)-HR%X) - &
                           Fc_x*(y_tmp(iii)-HR%Y)
                      P2P_U = P2P_U + Fc_x
                      P2P_V = P2P_V + Fc_y
                   End If
                End Do

                Exit TIM_PX
             End If
          End Do TIM_PX

          ! Find the closest wall at the -y direction
          y_now = -DY(jjn)
          TIM_MY: Do jj = jjn,0,-1
             y_now = y_now + DY(jj)
             If ( y_now-HR%Radius > P2P_DIST_MAX) Exit TIM_MY
             tim_ic = cell_index(iin,jj,kkn)
             Call Get_iw(iin,jj,kkn,+2,tim_iw)
             If ( tim_iw == 0 ) Then
                Cycle TIM_MY
             Else
                IBC = IJKW(5,tim_iw)
                If (SURFACE(IBC)%VEL> 0.0_EB .Or. BOUNDARY_TYPE(tim_iw)==OPEN_BOUNDARY) Exit TIM_MY
                tim_dist = y1 - yw(tim_iw)
                If ( (HR%U**2 +HR%V**2) > 0.0_EB ) Then
                   CosPhiFac = (-HR%V)/Sqrt(HR%U**2+HR%V**2)
                   CosPhiFac = LambdaW +  &
                        0.5_EB*(1.0_EB-LambdaW)*(1.0_EB+CosPhiFac)
                Else
                   CosPhiFac = 1.0_EB
                End If
                P2P_V = P2P_V + A_Wall*CosPhiFac* &
                     Exp( -(tim_dist-HR%Radius)/B_Wall )
                Social_F = Social_F + Abs(A_Wall*CosPhiFac* &
                     Exp( -(tim_dist-HR%Radius)/B_Wall ))
                Do iii = 1, 3
                   tim_dist = y_tmp(iii) - yw(tim_iw)
                   d_walls = Min(tim_dist-r_tmp(iii),d_walls)
                   If (tim_dist <= r_tmp(iii)) Then
                      Fc_y = + 2.0_EB*HR%C_Young*(r_tmp(iii)-tim_dist)
                      Fc_x = 0.0_EB
                      Contact_F = Contact_F + Abs(Fc_y)
                      If (I_Fric_sw >= 1 ) Then
                         Fc_x = Fc_x - HR%Kappa*(r_tmp(iii)-tim_dist)*u_tmp(iii)
                      Else
                         Fc_x = Fc_x - HR%Gamma*u_tmp(iii)
                      End If
                      P2P_Torque = P2P_Torque + Fc_y*(x_tmp(iii)-HR%X) - &
                           Fc_x*(yw(tim_iw)-HR%Y)
                      P2P_U = P2P_U + Fc_x
                      P2P_V = P2P_V + Fc_y
                   End If
                End Do

                Exit TIM_MY
             End If
          End Do TIM_MY

          ! Find the closest wall at the +y direction
          y_now = -DY(jjn)
          TIM_PY: Do jj = jjn, JBAR+1
             y_now = y_now + DY(jj)
             If ( y_now-HR%Radius > P2P_DIST_MAX) Exit TIM_PY
             tim_ic = cell_index(iin,jj,kkn)
             Call Get_iw(iin,jj,kkn,-2,tim_iw)
             If ( tim_iw == 0 ) Then
                Cycle TIM_PY
             Else
                IBC = IJKW(5,tim_iw)
                If (SURFACE(IBC)%VEL> 0.0_EB .Or. BOUNDARY_TYPE(tim_iw)==OPEN_BOUNDARY) Exit TIM_PY
                tim_dist = yw(tim_iw) - y1
                If ( (HR%U**2 +HR%V**2) > 0.0_EB ) Then
                   CosPhiFac = (HR%V)/Sqrt(HR%U**2+HR%V**2)
                   CosPhiFac = LambdaW +  &
                        0.5_EB*(1.0_EB-LambdaW)*(1.0_EB+CosPhiFac)
                Else
                   CosPhiFac = 1.0_EB
                End If
                P2P_V = P2P_V - A_Wall*CosPhiFac* &
                     Exp( -(tim_dist-HR%Radius)/B_Wall )
                Social_F = Social_F + Abs(A_Wall*CosPhiFac* &
                     Exp( -(tim_dist-HR%Radius)/B_Wall ))
                Do iii = 1, 3
                   tim_dist = yw(tim_iw) - y_tmp(iii)
                   d_walls = Min(tim_dist-r_tmp(iii),d_walls)
                   If (tim_dist <= r_tmp(iii)) Then
                      Fc_y = - 2.0_EB*HR%C_Young*(r_tmp(iii)-tim_dist)
                      Fc_x = 0.0_EB
                      Contact_F = Contact_F + Abs(Fc_y)
                      If (I_Fric_sw >= 1 ) Then
                         Fc_x = Fc_x - HR%Kappa*(r_tmp(iii)-tim_dist)*u_tmp(iii)
                      Else
                         Fc_x = Fc_x - HR%Gamma*u_tmp(iii)
                      End If
                      P2P_Torque = P2P_Torque + Fc_y*(x_tmp(iii)-HR%X) - &
                           Fc_x*(yw(tim_iw)-HR%Y)
                      P2P_U = P2P_U + Fc_x
                      P2P_V = P2P_V + Fc_y
                   End If
                End Do

                Exit TIM_PY
             End If
          End Do TIM_PY

          ! ========================================================
          ! WALL SURFACE - PERSON FORCES ENDS
          ! ========================================================

          ! ========================================================
          ! WALL CORNER - PERSON FORCES STARTS
          ! ========================================================

          ! top right corner (x > x_human, y > y_human)
          x_now = -DX(iin+1)
          Loop_px: Do ii = iin+1, IBAR+1
             x_now = x_now + DX(ii)
             If (x_now-HR%Radius > P2P_DIST_MAX) Exit Loop_px
             y_now = -DY(jjn-1)
             Loop_pxpy: Do jj = jjn+1, JBAR+1
                y_now = y_now + DY(jj)
                If (Sqrt(x_now**2 + y_now**2)-HR%Radius &
                     > P2P_DIST_MAX) Exit Loop_pxpy
                tim_ic = cell_index(ii,jj,kkn)
                Call Get_iw(ii,jj,kkn,-1,tim_iwx)
                Call Get_iw(ii,jj,kkn,-2,tim_iwy)
                If ( (tim_iwx /= 0) .And. (tim_iwy /= 0) ) Then
                   tim_dist = Sqrt( (yw(tim_iwy) - y1)**2 + &
                        (xw(tim_iwx) - x1)**2 )
                   If (tim_dist-HR%Radius > P2P_DIST_MAX) &
                        Exit Loop_pxpy
                   If ( (HR%U**2 +HR%V**2) > 0.0_EB ) Then
                      CosPhiFac = ((xw(tim_iwx)-X1)*HR%U + &
                           (yw(tim_iwy)-Y1)*HR%V ) &
                           / ( tim_dist* Sqrt(HR%U**2 +HR%V**2) )
                      CosPhiFac = LambdaW + &
                           0.5_EB*(1.0_EB-LambdaW)*(1.0_EB+CosPhiFac)
                   Else
                      CosPhiFac = 1.0_EB
                   End If
                   P2P_U = P2P_U + (X1-xw(tim_iwx))*A_Wall*CosPhiFac* &
                        Exp( -(tim_dist-HR%Radius)/B_Wall ) / &
                        tim_dist
                   P2P_V = P2P_V + (Y1-yw(tim_iwy))*A_Wall*CosPhiFac* &
                        Exp( -(tim_dist-HR%Radius)/B_Wall ) / &
                        tim_dist
                   Social_F = Social_F + Abs(A_Wall*CosPhiFac* &
                        Exp( -(tim_dist-HR%Radius)/B_Wall ))
                   Do iii = 1, 3
                      tim_dist = Sqrt( (yw(tim_iwy) - y_tmp(iii))**2 + &
                           (xw(tim_iwx) - x_tmp(iii))**2 )
                      d_walls = Min(tim_dist-r_tmp(iii),d_walls)
                      If (tim_dist <= r_tmp(iii)) Then
                         Fc_x = (x_tmp(iii)-xw(tim_iwx))* &
                              2.0_EB*HR%C_Young*(r_tmp(iii)-tim_dist)/tim_dist
                         Fc_y = (y_tmp(iii)-yw(tim_iwy))* &
                              2.0_EB*HR%C_Young*(r_tmp(iii)-tim_dist)/tim_dist
                         Contact_F = Contact_F + Sqrt( Fc_x**2 + Fc_y**2 )
                         ! Tangential contact force:
                         If (I_Fric_sw >= 1 ) Then
                            Fc_x = Fc_x - HR%Kappa* &
                                 (r_tmp(iii)-tim_dist)*(y_tmp(iii)-yw(tim_iwy))* &
                                 ( (y_tmp(iii)-yw(tim_iwy))*u_tmp(iii) - &
                                 (x_tmp(iii)-xw(tim_iwx))*v_tmp(iii) )/tim_dist**2
                            Fc_y = Fc_y + HR%Kappa* &
                                 (r_tmp(iii)-tim_dist)*(x_tmp(iii)-xw(tim_iwx))* &
                                 ( (y_tmp(iii)-yw(tim_iwy))*u_tmp(iii) - &
                                 (x_tmp(iii)-xw(tim_iwx))*v_tmp(iii) )/tim_dist**2
                         Else
                            Fc_x = Fc_x - HR%Gamma*(y_tmp(iii)-yw(tim_iwy))* &
                                 ( (y_tmp(iii)-yw(tim_iwy))*u_tmp(iii) - &
                                 (x_tmp(iii)-xw(tim_iwx))*v_tmp(iii) )/tim_dist**2
                            Fc_y = Fc_y + HR%Gamma*(x_tmp(iii)-xw(tim_iwx))* &
                                 ( (y_tmp(iii)-yw(tim_iwy))*u_tmp(iii) - &
                                 (x_tmp(iii)-xw(tim_iwx))*v_tmp(iii) )/tim_dist**2
                         End If
                         P2P_Torque = P2P_Torque + Fc_y*(xw(tim_iwx)-x_tmp(iii)) - &
                              Fc_x*(yw(tim_iwy)-y_tmp(iii))
                         P2P_U = P2P_U + Fc_x
                         P2P_V = P2P_V + Fc_y
                      End If
                   End Do

                   Exit Loop_pxpy
                End If
                If ( Solid(tim_ic) ) Exit Loop_pxpy
             End Do Loop_pxpy
          End Do Loop_px

          ! top left corner (x < x_human, y > y_human)
          x_now = -DX(iin-1)
          Loop_mx: Do ii = iin-1, 0, -1
             x_now = x_now + DX(ii)
             If (x_now-HR%Radius > P2P_DIST_MAX) Exit Loop_mx
             y_now = -DY(jjn-1)
             Loop_mxpy: Do jj = jjn+1, JBAR+1
                y_now = y_now + DY(jj)
                If (Sqrt(x_now**2 + y_now**2)-HR%Radius &
                     > P2P_DIST_MAX) Exit Loop_mxpy
                tim_ic = cell_index(ii,jj,kkn)
                Call Get_iw(ii,jj,kkn,+1,tim_iwx)
                Call Get_iw(ii,jj,kkn,-2,tim_iwy)
                If ( (tim_iwx /= 0) .And. (tim_iwy /= 0) ) Then
                   tim_dist = Sqrt( (yw(tim_iwy) - y1)**2 + &
                        (xw(tim_iwx) - x1)**2 )
                   If (tim_dist-HR%Radius > P2P_DIST_MAX) &
                        Exit Loop_mxpy
                   If ( (HR%U**2 +HR%V**2) > 0.0_EB ) Then
                      CosPhiFac = ((xw(tim_iwx)-X1)*HR%U + &
                           (yw(tim_iwy)-Y1)*HR%V ) &
                           / ( tim_dist* Sqrt(HR%U**2 +HR%V**2) )
                      CosPhiFac = LambdaW + &
                           0.5_EB*(1.0_EB-LambdaW)*(1.0_EB+CosPhiFac)
                   Else
                      CosPhiFac = 1.0_EB
                   End If
                   P2P_U = P2P_U + (X1-xw(tim_iwx))*A_Wall*CosPhiFac* &
                        Exp( -(tim_dist-HR%Radius)/B_Wall ) / &
                        tim_dist
                   P2P_V = P2P_V + (Y1-yw(tim_iwy))*A_Wall*CosPhiFac* &
                        Exp( -(tim_dist-HR%Radius)/B_Wall ) / &
                        tim_dist
                   Social_F = Social_F + Abs(A_Wall*CosPhiFac* &
                        Exp( -(tim_dist-HR%Radius)/B_Wall ))
                   Do iii = 1, 3
                      tim_dist = Sqrt( (yw(tim_iwy) - y_tmp(iii))**2 + &
                           (xw(tim_iwx) - x_tmp(iii))**2 )
                      d_walls = Min(tim_dist-r_tmp(iii),d_walls)
                      If (tim_dist <= r_tmp(iii)) Then
                         Fc_x = (x_tmp(iii)-xw(tim_iwx))* &
                              2.0_EB*HR%C_Young*(r_tmp(iii)-tim_dist)/tim_dist
                         Fc_y = (y_tmp(iii)-yw(tim_iwy))* &
                              2.0_EB*HR%C_Young*(r_tmp(iii)-tim_dist)/tim_dist
                         Contact_F = Contact_F + Sqrt( Fc_x**2 + Fc_y**2 )
                         ! Tangential contact force:
                         If (I_Fric_sw >= 1 ) Then
                            Fc_x = Fc_x - HR%Kappa* &
                                 (r_tmp(iii)-tim_dist)*(y_tmp(iii)-yw(tim_iwy))* &
                                 ( (y_tmp(iii)-yw(tim_iwy))*u_tmp(iii) - &
                                 (x_tmp(iii)-xw(tim_iwx))*v_tmp(iii) )/tim_dist**2
                            Fc_y = Fc_y + HR%Kappa* &
                                 (r_tmp(iii)-tim_dist)*(x_tmp(iii)-xw(tim_iwx))* &
                                 ( (y_tmp(iii)-yw(tim_iwy))*u_tmp(iii) - &
                                 (x_tmp(iii)-xw(tim_iwx))*v_tmp(iii) )/tim_dist**2
                         Else
                            Fc_x = Fc_x - HR%Gamma*(y_tmp(iii)-yw(tim_iwy))* &
                                 ( (y_tmp(iii)-yw(tim_iwy))*u_tmp(iii) - &
                                 (x_tmp(iii)-xw(tim_iwx))*v_tmp(iii) )/tim_dist**2
                            Fc_y = Fc_y + HR%Gamma*(x_tmp(iii)-xw(tim_iwx))* &
                                 ( (y_tmp(iii)-yw(tim_iwy))*u_tmp(iii) - &
                                 (x_tmp(iii)-xw(tim_iwx))*v_tmp(iii) )/tim_dist**2
                         End If
                         P2P_Torque = P2P_Torque + Fc_y*(xw(tim_iwx)-x_tmp(iii)) - &
                              Fc_x*(yw(tim_iwy)-y_tmp(iii))
                         P2P_U = P2P_U + Fc_x
                         P2P_V = P2P_V + Fc_y
                      End If
                   End Do

                   Exit Loop_mxpy
                End If
                If ( Solid(tim_ic) ) Exit Loop_mxpy
             End Do Loop_mxpy
          End Do Loop_mx

          ! bottom right corner (x > x_human, y < y_human)
          x_now = -DX(iin+1)
          Loop_py: Do ii = iin+1, IBAR+1
             x_now = x_now + DX(ii)
             If (x_now-HR%Radius > P2P_DIST_MAX) Exit Loop_py
             y_now = -DY(jjn-1)
             Loop_pxmy: Do jj = jjn-1, 0, -1
                y_now = y_now + DY(jj)
                If (Sqrt(x_now**2 + y_now**2)-HR%Radius &
                     > P2P_DIST_MAX) Exit Loop_pxmy
                tim_ic = cell_index(ii,jj,kkn)
                Call Get_iw(ii,jj,kkn,-1,tim_iwx)
                Call Get_iw(ii,jj,kkn,+2,tim_iwy)
                If ( (tim_iwx /= 0) .And. (tim_iwy /= 0) ) Then
                   tim_dist = Sqrt( (yw(tim_iwy) - y1)**2 + &
                        (xw(tim_iwx) - x1)**2 )
                   If (tim_dist-HR%Radius > P2P_DIST_MAX) &
                        Exit Loop_pxmy
                   If ( (HR%U**2 +HR%V**2) > 0.0_EB ) Then
                      CosPhiFac = ((xw(tim_iwx)-X1)*HR%U + &
                           (yw(tim_iwy)-Y1)*HR%V ) &
                           / ( tim_dist* Sqrt(HR%U**2 +HR%V**2) )
                      CosPhiFac = LambdaW + &
                           0.5_EB*(1.0_EB-LambdaW)*(1.0_EB+CosPhiFac)
                   Else
                      CosPhiFac = 1.0_EB
                   End If
                   P2P_U = P2P_U + (X1-xw(tim_iwx))*A_Wall*CosPhiFac* &
                        Exp( -(tim_dist-HR%Radius)/B_Wall ) / &
                        tim_dist
                   P2P_V = P2P_V + (Y1-yw(tim_iwy))*A_Wall*CosPhiFac* &
                        Exp( -(tim_dist-HR%Radius)/B_Wall ) / &
                        tim_dist
                   Social_F = Social_F + Abs(A_Wall*CosPhiFac* &
                        Exp( -(tim_dist-HR%Radius)/B_Wall ))
                   Do iii = 1, 3
                      tim_dist = Sqrt( (yw(tim_iwy) - y_tmp(iii))**2 + &
                           (xw(tim_iwx) - x_tmp(iii))**2 )
                      d_walls = Min(tim_dist-r_tmp(iii),d_walls)
                      If (tim_dist <= r_tmp(iii)) Then
                         Fc_x = (x_tmp(iii)-xw(tim_iwx))* &
                              2.0_EB*HR%C_Young*(r_tmp(iii)-tim_dist)/tim_dist
                         Fc_y = (y_tmp(iii)-yw(tim_iwy))* &
                              2.0_EB*HR%C_Young*(r_tmp(iii)-tim_dist)/tim_dist
                         Contact_F = Contact_F + Sqrt( Fc_x**2 + Fc_y**2 )
                         ! Tangential contact force:
                         If (I_Fric_sw >= 1 ) Then
                            Fc_x = Fc_x - HR%Kappa* &
                                 (r_tmp(iii)-tim_dist)*(y_tmp(iii)-yw(tim_iwy))* &
                                 ( (y_tmp(iii)-yw(tim_iwy))*u_tmp(iii) - &
                                 (x_tmp(iii)-xw(tim_iwx))*v_tmp(iii) )/tim_dist**2
                            Fc_y = Fc_y + HR%Kappa* &
                                 (r_tmp(iii)-tim_dist)*(x_tmp(iii)-xw(tim_iwx))* &
                                 ( (y_tmp(iii)-yw(tim_iwy))*u_tmp(iii) - &
                                 (x_tmp(iii)-xw(tim_iwx))*v_tmp(iii) )/tim_dist**2
                         Else
                            Fc_x = Fc_x - HR%Gamma*(y_tmp(iii)-yw(tim_iwy))* &
                                 ( (y_tmp(iii)-yw(tim_iwy))*u_tmp(iii) - &
                                 (x_tmp(iii)-xw(tim_iwx))*v_tmp(iii) )/tim_dist**2
                            Fc_y = Fc_y + HR%Gamma*(x_tmp(iii)-xw(tim_iwx))* &
                                 ( (y_tmp(iii)-yw(tim_iwy))*u_tmp(iii) - &
                                 (x_tmp(iii)-xw(tim_iwx))*v_tmp(iii) )/tim_dist**2
                         End If
                         P2P_Torque = P2P_Torque + Fc_y*(xw(tim_iwx)-x_tmp(iii)) - &
                              Fc_x*(yw(tim_iwy)-y_tmp(iii))
                         P2P_U = P2P_U + Fc_x
                         P2P_V = P2P_V + Fc_y
                      End If
                   End Do

                   Exit Loop_pxmy
                End If
                If ( Solid(tim_ic) ) Exit Loop_pxmy
             End Do Loop_pxmy
          End Do Loop_py

          ! bottom left corner (x < x_human, y < y_human)
          x_now = -DX(iin-1)
          Loop_my: Do ii = iin-1, 0, -1
             x_now = x_now + DX(ii)
             If (x_now-HR%Radius > P2P_DIST_MAX) Exit Loop_my
             y_now = -DY(jjn-1)
             Loop_mxmy: Do jj = jjn-1, 0, -1
                y_now = y_now + DY(jj)
                If (Sqrt(x_now**2 + y_now**2)-HR%Radius &
                     > P2P_DIST_MAX) Exit Loop_mxmy
                tim_ic = cell_index(ii,jj,kkn)
                Call Get_iw(ii,jj,kkn,+1,tim_iwx)
                Call Get_iw(ii,jj,kkn,+2,tim_iwy)
                If ( (tim_iwx /= 0) .And. (tim_iwy /= 0) ) Then
                   tim_dist = Sqrt( (yw(tim_iwy) - y1)**2 + &
                        (xw(tim_iwx) - x1)**2 )
                   If (tim_dist-HR%Radius > P2P_DIST_MAX) &
                        Exit Loop_mxmy
                   If ( (HR%U**2 +HR%V**2) > 0.0_EB ) Then
                      CosPhiFac = ((xw(tim_iwx)-X1)*HR%U + &
                           (yw(tim_iwy)-Y1)*HR%V ) &
                           / ( tim_dist* Sqrt(HR%U**2 +HR%V**2) )
                      CosPhiFac = LambdaW + &
                           0.5_EB*(1.0_EB-LambdaW)*(1.0_EB+CosPhiFac)
                   Else
                      CosPhiFac = 1.0_EB
                   End If
                   P2P_U = P2P_U + (X1-xw(tim_iwx))*A_Wall*CosPhiFac* &
                        Exp( -(tim_dist-HR%Radius)/B_Wall ) / &
                        tim_dist
                   P2P_V = P2P_V + (Y1-yw(tim_iwy))*A_Wall*CosPhiFac* &
                        Exp( -(tim_dist-HR%Radius)/B_Wall ) / &
                        tim_dist
                   Social_F = Social_F + Abs(A_Wall*CosPhiFac* &
                        Exp( -(tim_dist-HR%Radius)/B_Wall ))
                   Do iii = 1, 3
                      tim_dist = Sqrt( (yw(tim_iwy) - y_tmp(iii))**2 + &
                           (xw(tim_iwx) - x_tmp(iii))**2 )
                      d_walls = Min(tim_dist-r_tmp(iii),d_walls)
                      If (tim_dist <= r_tmp(iii)) Then
                         Fc_x = (x_tmp(iii)-xw(tim_iwx))* &
                              2.0_EB*HR%C_Young*(r_tmp(iii)-tim_dist)/tim_dist
                         Fc_y = (y_tmp(iii)-yw(tim_iwy))* &
                              2.0_EB*HR%C_Young*(r_tmp(iii)-tim_dist)/tim_dist
                         Contact_F = Contact_F + Sqrt( Fc_x**2 + Fc_y**2 )
                         ! Tangential contact force:
                         If (I_Fric_sw >= 1 ) Then
                            Fc_x = Fc_x - HR%Kappa* &
                                 (r_tmp(iii)-tim_dist)*(y_tmp(iii)-yw(tim_iwy))* &
                                 ( (y_tmp(iii)-yw(tim_iwy))*u_tmp(iii) - &
                                 (x_tmp(iii)-xw(tim_iwx))*v_tmp(iii) )/tim_dist**2
                            Fc_y = Fc_y + HR%Kappa* &
                                 (r_tmp(iii)-tim_dist)*(x_tmp(iii)-xw(tim_iwx))* &
                                 ( (y_tmp(iii)-yw(tim_iwy))*u_tmp(iii) - &
                                 (x_tmp(iii)-xw(tim_iwx))*v_tmp(iii) )/tim_dist**2
                         Else
                            Fc_x = Fc_x - HR%Gamma*(y_tmp(iii)-yw(tim_iwy))* &
                                 ( (y_tmp(iii)-yw(tim_iwy))*u_tmp(iii) - &
                                 (x_tmp(iii)-xw(tim_iwx))*v_tmp(iii) )/tim_dist**2
                            Fc_y = Fc_y + HR%Gamma*(x_tmp(iii)-xw(tim_iwx))* &
                                 ( (y_tmp(iii)-yw(tim_iwy))*u_tmp(iii) - &
                                 (x_tmp(iii)-xw(tim_iwx))*v_tmp(iii) )/tim_dist**2
                         End If
                         P2P_Torque = P2P_Torque + Fc_y*(xw(tim_iwx)-x_tmp(iii)) - &
                              Fc_x*(yw(tim_iwy)-y_tmp(iii))
                         P2P_U = P2P_U + Fc_x
                         P2P_V = P2P_V + Fc_y
                      End If
                   End Do

                   Exit Loop_mxmy
                End If
                If ( Solid(tim_ic) ) Exit Loop_mxmy
             End Do Loop_mxmy
          End Do Loop_my

          ! ========================================================
          ! WALL CORNER - PERSON FORCES ENDS
          ! ========================================================
          !
          ! ========================================================
          ! WALL FORCES ENDS HERE
          ! ========================================================
          HR%F_X = P2P_U
          HR%F_Y = P2P_V
          HR%SumForces = Contact_F ! 21.9.2006 by T.K.
          HR%SumForces2 = Social_F + Contact_F
          HR%Torque = P2P_Torque
          !


          If ( T <= 0.0_EB ) Then
             If ( Abs(P2P_U)/HR%Mass > 550.0_EB ) &
                  P2P_U =550.0_EB*HR%Mass*P2P_U/Abs(P2P_U)
             If ( Abs(P2P_V)/HR%Mass > 550.0_EB ) &
                  P2P_V =550.0_EB*HR%Mass*P2P_V/Abs(P2P_V)
             HR%F_X = P2P_U
             HR%F_Y = P2P_V
          End If

          U_new = HR%U + 0.5_EB*HR%F_X*DTSP/HR%Mass
          V_new = HR%V + 0.5_EB*HR%F_Y*DTSP/HR%Mass
          Omega_new = HR%Omega + 0.5_EB*DTSP*HR%Torque/HR%M_iner

          ! ========================================================
          ! Some random force here
          ! ========================================================
          If (GaTh > 0.0_EB .And. T > 0.0_EB ) Then
             U_new = U_new + 0.5_EB*DTSP*HR%v0_fac* &
                  HR%Mass*HR%ksi*Cos(HR%eta)/HR%Mass
             V_new = V_new + 0.5_EB*DTSP*HR%v0_fac* &
                  HR%Mass*HR%ksi*Sin(HR%eta)/HR%Mass
             P2P_U = P2P_U + HR%v0_fac*HR%Mass*HR%ksi*Cos(HR%eta)
             P2P_V = P2P_V + HR%v0_fac*HR%Mass*HR%ksi*Sin(HR%eta)
             Omega_new = Omega_new + 0.5_EB*DTSP* &
                  1.0_EB*Sign(HR%ksi,HR%eta-Pi)
             P2P_Torque = P2P_Torque + 1.0_EB*Sign(HR%ksi,HR%eta-Pi)*HR%M_iner
          Else
             HR%ksi = 0.0_EB
             HR%eta = 0.0_EB
          End If


          fac_tim =  1.0_EB + (DTSP/(2.0_EB*hr_tau))
          ! 
          ! ========================================================
          ! Add self-propelling force terms, self-consistent
          ! ========================================================

          j = Max(0,HR%GROUP_ID)
          If (j == 0 ) Then
             Group_List(0)%Tpre = HR%Tpre + HR%Tdet
             Tpre = HR%Tpre + HR%Tdet
          Else
             Tpre = HR%Tdet
          End If


          If (Group_List(j)%GROUP_SIZE >= 2) Then
             HR%UBAR_Center = (Group_List(j)%GROUP_X - HR%X)
             HR%VBAR_Center = (Group_List(j)%GROUP_Y - HR%Y)
             EVEL = Sqrt(HR%UBAR_Center**2 + HR%VBAR_Center**2)
             If ( EVEL > 0.0_EB ) Then
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
             If (Group_List(j)%COMPLETE == 1 .And. T+DTSP_new &
                  <= Group_List(j)%Tpre+Group_List(j)%Tdoor) Then
                UBAR = 0.0_EB
                VBAR = 0.0_EB
             End If
          End If


          EVEL = Sqrt(UBAR**2 + VBAR**2)
          If (Group_List(j)%COMPLETE == 1 .And. EVEL > 0.0_EB) Then
             UBAR = (1-GROUP_EFF)*UBAR + GROUP_EFF*HR%UBAR_Center/ &
                  Sqrt(((1-GROUP_EFF)*UBAR + GROUP_EFF*HR%UBAR_Center)**2+ &
                  ((1-GROUP_EFF)*VBAR + GROUP_EFF*HR%VBAR_Center)**2)
             VBAR = (1-GROUP_EFF)*VBAR + GROUP_EFF*HR%VBAR_Center/ &
                  Sqrt(((1-GROUP_EFF)*UBAR + GROUP_EFF*HR%UBAR_Center)**2+ &
                  ((1-GROUP_EFF)*VBAR + GROUP_EFF*HR%VBAR_Center)**2)
          Else
             UBAR = HR%UBAR_Center
             VBAR = HR%VBAR_Center
          End If

          If (T+DTSP_new <= Tpre .Or. T+DTSP_new <= 0.0_EB) Then
             ! Fire is not yet detected ==> no movement
             UBAR = 0.0_EB
             VBAR = 0.0_EB
          End If
          If ( T <= Tpre ) Then
             If ( (T+DTSP_new) > Tpre) Then
                EVEL = Sqrt(UBAR**2 + VBAR**2)
                If (EVEL > 0.0_EB) Then
                   speed = speed_xp*(0.5_EB + Sign(0.5_EB,UBAR)) + &
                        speed_xm*(0.5_EB - Sign(0.5_EB,UBAR)) 
                   speed = speed*HR%v0_fac
                   P2P_U = P2P_U + (HR%Mass/hr_tau)*speed*(UBAR/EVEL)
                   speed = speed_yp*(0.5_EB + Sign(0.5_EB,VBAR)) + &
                        speed_ym*(0.5_EB - Sign(0.5_EB,VBAR)) 
                   speed = speed*HR%v0_fac
                   P2P_V = P2P_V + (HR%Mass/hr_tau)*speed*(VBAR/EVEL)
                End If
             End If
             UBAR = 0.0_EB
             VBAR = 0.0_EB
          End If

          EVEL = Sqrt(UBAR**2 + VBAR**2)
          If (EVEL > 0.0_EB) Then
             speed = speed_xp*(0.5_EB + Sign(0.5_EB,UBAR)) + &
                  speed_xm*(0.5_EB - Sign(0.5_EB,UBAR)) 
             speed = speed*HR%v0_fac
             U_new = (U_new + 0.5_EB*(DTSP/hr_tau)* &
                  speed*(UBAR/EVEL)) / fac_tim
             HR%UBAR = speed*(UBAR/EVEL)
             P2P_U = P2P_U + (HR%Mass/hr_tau)* &
                  (speed*(UBAR/EVEL) - HR%U)

             speed = speed_yp*(0.5_EB + Sign(0.5_EB,VBAR)) + &
                  speed_ym*(0.5_EB - Sign(0.5_EB,VBAR)) 
             speed = speed*HR%v0_fac
             V_new = (V_new + 0.5_EB*(DTSP/hr_tau)* &
                  speed*(VBAR/EVEL)) / fac_tim
             HR%VBAR = speed*(VBAR/EVEL)
             P2P_V = P2P_V + (HR%Mass/hr_tau)* &
                  (speed*(VBAR/EVEL) - HR%V)

             If (VBAR >= 0.0_EB) Then
                angle = ACos(UBAR/EVEL)
             Else
                angle = 2.0_EB*Pi - ACos(UBAR/EVEL)
             End If  ! angle is [0,2Pi)
             If (angle == 2.0_EB*Pi) angle = 0.0_EB

             ! Next is not working, because the angle-force is harmonic, i.e.,
             ! human has harmonic rotational oscilations, no damping at all

             ! J(dw/dt) = (J/t_iner)*( ((angle-angle_0)/pi)w_0 - w )
             If (Abs(angle-HR%angle) <= Pi ) Then
                ! zero is not crossed.
                Omega_new = Omega_new + 0.5_EB*(DTSP/HR%tau_iner)* &
                     ( (angle-HR%angle)*(Omega_0/Pi) - HR%Omega)
                P2P_Torque = P2P_Torque + (HR%M_iner/HR%tau_iner)* &
                     ( (angle-HR%angle)*(Omega_0/Pi) - HR%Omega)
             Else
                ! zero is crossed
                Omega_new = Omega_new + 0.5_EB*(DTSP/HR%Tau_iner)* &
                     ( (2.0_EB*Pi-Abs(angle-HR%angle))*Sign(1.0_EB , HR%angle-angle)* &
                     (Omega_0/Pi) - HR%Omega )
                P2P_Torque = P2P_Torque + (HR%M_iner/HR%Tau_iner)* &
                     ( (2.0_EB*Pi-Abs(angle-HR%angle))*Sign(1.0_EB , HR%angle-angle)* &
                     (Omega_0/Pi) - HR%Omega )
             End If
          Else
             HR%UBAR = 0.0_EB
             HR%VBAR = 0.0_EB
             U_new = U_new / fac_tim
             V_new = V_new / fac_tim
             Omega_new = Omega_new + 0.5_EB*(DTSP/HR%Tau_iner)*(-HR%Omega)
          End If

          ! ========================================================
          ! SELF PROPELLING FORCE ENDS HERE
          ! ========================================================
          ! Check that velocities are not too large, i.e.,
          ! unphysical (less than 10 m/s for humans)
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

          d_humans = Max(d_humans,0.0001_EB)
          d_walls  = Max(d_walls, 0.0001_EB)
          d_humans_min = Min(d_humans_min,d_humans)
          d_walls_min  = Min(d_walls_min, d_walls)

          If ( T > 0.0_EB ) Then
             ! Time step, do not move too close to other humans or 0.5*grid spacing
             dt_Loop: Do
                u_tmp(2) = HR%U + 0.5_EB*DTSP_new*P2P_U/HR%Mass
                v_tmp(2) = HR%V + 0.5_EB*DTSP_new*P2P_V/HR%Mass
                Omega_new = HR%Omega + 0.5_EB*DTSP*P2P_Torque/HR%M_iner
                u_tmp(1) = u_tmp(2) + Cos(HR%angle)*Omega_new*HR%d_shoulder
                u_tmp(3) = u_tmp(2) - Cos(HR%angle)*Omega_new*HR%d_shoulder
                v_tmp(1) = v_tmp(2) + Sin(HR%angle)*Omega_new*HR%d_shoulder
                v_tmp(3) = v_tmp(2) - Sin(HR%angle)*Omega_new*HR%d_shoulder
                If ( Max(u_tmp(1)**2+v_tmp(1)**2, u_tmp(2)**2+v_tmp(2)**2, &
                     u_tmp(3)**2+v_tmp(3)**2)*DTSP_new**2 > &
                     (Min(0.2_EB*d_humans, 0.2_EB*d_walls, 0.5_EB*Delta_min))**2 ) Then
                   DTSP_new = DTSP_new*0.8_EB
                   Cycle dt_Loop
                End If
                Exit dt_Loop
             End Do dt_Loop

             If (Social_F + Contact_F > 100.0_EB) Then
                EVAC_DT_MIN2 = Max( EVAC_DT_MIN, EVAC_DT_MAX*100.0_EB/(Social_F+Contact_F) )
             End If
             DTSP_new = Max(DTSP_new, EVAC_DT_MIN2)

          Else
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

       ! Maximum time step allowed by the SC-VV algorithm
       Tsteps(NM) = DTSP_new

       ! ========================================================
       ! HUMAN TIME LOOP ENDS HERE
       ! ========================================================
    End Do HUMAN_TIME_LOOP
    Deallocate(BLOCK_GRID_N)
    Deallocate(Color_Tmp)
    Deallocate(K_ave_Door)
    Deallocate(FED_max_Door)
    Deallocate(Is_Visible_Door)
    Deallocate(Is_Known_Door)
    ! ========================================================
    ! Evacuation routine ends here
    ! ========================================================
!!$    End If EVAC_MESH_ONLY
    !
    TUSED(12,NM)=TUSED(12,NM)+SECOND()-TNOW

  Contains



    SUBROUTINE GET_IW(IIin,JJin,KKin,IOR,IW)
      Implicit None

      Integer, Intent(IN) :: IIin, JJin, KKin, IOR
      Integer, Intent(OUT) :: IW
      ! Local variables
      Integer :: ii, jj, kk, ic
      !
      ii = IIin
      jj = JJin
      kk = KKin

      IC  = CELL_INDEX(II,JJ,KK)

      If (SOLID(IC)) Then
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
      Endif

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
      Endif

    End Subroutine GET_IW

    !
    Subroutine CHECK_EXITS(T,NM)
      Implicit None
      !
      ! Remove persons if they are found at an exit.
      !
      Real(EB), Intent(IN) :: T
      Integer, Intent(IN) :: NM
      Real(EB) x_old, y_old, pexx1, pexx2, pexy1, pexy2
      Integer :: ie,i,n_tmp
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
               If ((HR%X >= pex%x1 .And. x_old < pex%x1) .And. &
                    (HR%Y >= pex%y1 .And. HR%Y <= pex%y2) ) Then
                  If (PEX%COUNT_ONLY) HR%IOR = 2
                  If (.Not. PEX%COUNT_ONLY) HR%IOR = 1
                  PEX%T_last=T
                  PEX%ICOUNT = PEX%ICOUNT + 1
                  If (PEX%T_first <= 0.0_EB) PEX%T_first = T
               End If
            Case (-1)
               If ((HR%X <= pex%x2 .And. x_old > pex%x2) .And. &
                    (HR%Y >= pex%y1 .And. HR%Y <= pex%y2) ) Then
                  If (PEX%COUNT_ONLY) HR%IOR = 2
                  If (.Not. PEX%COUNT_ONLY) HR%IOR = 1
                  PEX%T_last=T
                  PEX%ICOUNT = PEX%ICOUNT + 1
                  If (PEX%T_first <= 0.0_EB) PEX%T_first = T
               End If
            Case (+2)
               If ((HR%Y >= pex%y1 .And. y_old < pex%y1) .And. &
                    (HR%X >= pex%x1 .And. HR%X <= pex%x2) ) Then
                  If (PEX%COUNT_ONLY) HR%IOR = 2
                  If (.Not. PEX%COUNT_ONLY) HR%IOR = 1
                  PEX%T_last=T
                  PEX%ICOUNT = PEX%ICOUNT + 1
                  If (PEX%T_first <= 0.0_EB) PEX%T_first = T
               End If
            Case (-2)
               If ((HR%Y <= pex%y2 .And. y_old > pex%y2) .And. &
                    (HR%X >= pex%x1 .And. HR%X <= pex%x2) ) Then
                  If (PEX%COUNT_ONLY) HR%IOR = 2
                  If (.Not. PEX%COUNT_ONLY) HR%IOR = 1
                  PEX%T_last=T
                  PEX%ICOUNT = PEX%ICOUNT + 1
                  If (PEX%T_first <= 0.0_EB) PEX%T_first = T
               End If
            End Select
            If (HR%IOR > 0) Then
               Write (LU_ERR,fmt='(a,i6,a,f8.2,a,a,a,f8.4)') &
                    ' EVAC: Person n:o', &
                    HR%ILABEL, ' out at ', T, &
                    ' s, exit ', Trim(PEX%ID_NAME), ' fed ', HR%IntDose
               If (HR%IOR == 2) HR%IOR = 0
            End If
         End Do HumLoop
      End Do PexLoop
      n_tmp = n_humans
      Remove_loop: Do i = n_tmp, 1, -1
         HR=>HUMAN(I)
         ! Remove person, if wanted.
         If (HR%IOR >= 1) Then
            HR%IOR = 0
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
      Real(EB), Intent(IN) :: T
      Integer, Intent(IN) :: NM
      Real(EB) x_old, y_old, xx, yy, zz, pdxx1, pdxx2, pdxy1, pdxy2, v, angle
      Integer :: ie,i,n_tmp, istat, ior_new, inode2, imesh2, n, ior
      Integer :: new_ffield_i, color_index, i_target, inode
      Character(60) :: TO_NODE
      Character(26) :: new_ffield_name
      Logical :: keep_xy
      !
      keep_xy = .False.
      HUMAN(:)%IOR = 0
      PdxLoop: Do ie = 1, n_doors
         PDX=>EVAC_DOORS(ie)
         If (PDX%IMESH /= NM ) Cycle PdxLoop
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
            If ( HR%IOR /= 0 ) Cycle HumLoop
            x_old = HR%X_old
            y_old = HR%Y_old
            Select Case (PDX%IOR)
            Case (+1)
               If ((HR%X >= pdx%x1 .And. x_old < pdx%x1) .And. &
                    (HR%Y >= pdx%y1 .And. HR%Y <= pdx%y2) ) Then
                  HR%IOR = 1
               End If
            Case (-1)
               If ((HR%X <= pdx%x2 .And. x_old > pdx%x2) .And. &
                    (HR%Y >= pdx%y1 .And. HR%Y <= pdx%y2) ) Then
                  HR%IOR = 1
               End If
            Case (+2)
               If ((HR%Y >= pdx%y1 .And. y_old < pdx%y1) .And. &
                    (HR%X >= pdx%x1 .And. HR%X <= pdx%x2) ) Then
                  HR%IOR = 1
               End If
            Case (-2)
               If ((HR%Y <= pdx%y2 .And. y_old > pdx%y2) .And. &
                    (HR%X >= pdx%x1 .And. HR%X <= pdx%x2) ) Then
                  HR%IOR = 1
               End If
            End Select
            If ( HR%IOR == 1 ) Then
               istat = 0
               inode = PDX%INODE
               inode2 = PDX%INODE2
               Call Check_Target_Node(inode,inode2,istat,xx,yy,zz,ior_new, &
                    imesh2,T,new_ffield_name,new_ffield_i,color_index,&
                    i_target, keep_xy, angle)

               If (istat == 0 ) Then
                  ! Put person to a new node, i.e., target node has empty space
                  HR%X = xx
                  HR%Y = yy
                  HR%Z = zz
                  HR%Angle = angle
                  If (keep_xy .And. ior_new == 1) Then
                     If (EVAC_Node_List(INODE2)%Node_Type == 'Door') Then
                        xx = 0.5_EB*(EVAC_DOORS(EVAC_Node_List(INODE2)%Node_Index)%X1 + &
                             EVAC_DOORS(EVAC_Node_List(INODE2)%Node_Index)%X2 - &
                             (PDX%X1+PDX%X2))
                        yy = 0.5_EB*(EVAC_DOORS(EVAC_Node_List(INODE2)%Node_Index)%Y1 + &
                             EVAC_DOORS(EVAC_Node_List(INODE2)%Node_Index)%Y2 - &
                             (PDX%Y1+PDX%Y2))
                     Else
                        xx = 0.5_EB*(EVAC_ENTRYS(EVAC_Node_List(INODE2)%Node_Index)%X1 + &
                             EVAC_ENTRYS(EVAC_Node_List(INODE2)%Node_Index)%X2 - &
                             (PDX%X1+PDX%X2))
                        yy = 0.5_EB*(EVAC_ENTRYS(EVAC_Node_List(INODE2)%Node_Index)%Y1 + &
                             EVAC_ENTRYS(EVAC_Node_List(INODE2)%Node_Index)%Y2 - &
                             (PDX%Y1+PDX%Y2))
                     End If
                     If (Abs(PDX%ior) == 2) HR%X = HR%X + xx  
                     If (Abs(PDX%ior) == 1) HR%Y = HR%Y + yy 
                  End If
                  v = Sqrt( HR%U**2 + HR%V**2 )

                  HR%IOR = ior_new
                  If (ior_new == 1) Then   ! door or entry
                     If (HR%IMESH /= imesh2 ) Then
                        ! Put the person to a new floor
                        HR%IOR         = 2
                        HR%INODE = 0
                        Do n = 1, n_egrids
                           If (EVAC_Node_List(n)%Mesh_Index == imesh2) &
                                HR%INODE = n
                        End Do
                        HR%NODE_NAME   = Trim(MESH_NAME(imesh2))
                        HR%FFIELD_NAME = Trim(new_ffield_name)
                        HR%I_FFIELD = new_ffield_i
                        HR%I_Target = I_Target
                     Else
                        HR%IOR = -2   ! same floor (door/entry)
                        HR%FFIELD_NAME = Trim(new_ffield_name)
                        HR%I_FFIELD = new_ffield_i
                        HR%I_Target = I_Target
                     End If
                     HR%X_old = HR%X
                     HR%Y_old = HR%Y
                     ! ior is the direction where the human is ejected.
                     If (EVAC_Node_List(INODE2)%Node_Type == 'Door') Then
                        ior = - EVAC_DOORS( &
                             EVAC_Node_List(INODE2)%Node_Index)%IOR
                     End If
                     If (EVAC_Node_List(INODE2)%Node_Type == 'Entry') Then
                        ior = EVAC_ENTRYS( &
                             EVAC_Node_List(INODE2)%Node_Index)%IOR
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
                        If (MESHES(imesh2)%N_HUMANS+1 > &
                             MESHES(imesh2)%N_HUMANS_DIM) Then
                           ! Re-allocation is not yet checked.
                           Call SHUTDOWN( &
                                'ERROR: Humans: no re-allocation yet')
                           Call RE_ALLOCATE_HUMANS(1,imesh2)
                        End If
                        MESHES(imesh2)%N_HUMANS = MESHES(imesh2)%N_HUMANS + 1
                        MESHES(imesh2)%HUMAN(MESHES(imesh2)%N_HUMANS) = HUMAN(I)
                        MESHES(imesh2)%HUMAN(MESHES(imesh2)%N_HUMANS)%IOR = 0
                     End If

                  End If            ! target is door or entry

                  PDX%T_last=T
                  PDX%ICOUNT = PDX%ICOUNT + 1
                  If (PDX%T_first <= 0.0_EB) PDX%T_first = T

                  Write (LU_ERR,fmt='(a,i6,a,f8.2,a,a,a,f8.4)') &
                       ' EVAC: Person n:o', &
                       HR%ILABEL, ' out at ', T, &
                       ' s, door ', Trim(PDX%ID_NAME), ' fed ', HR%IntDose

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
            HR%IOR = 0
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
      Real(EB), Intent(IN) :: T,DTSP
      Integer, Intent(IN) :: NM
      Real(EB) x_old, y_old, xx, yy, zz, pcxx1, pcxx2, pcxy1, pcxy2, &
           v, x_int, angle
      Integer :: ie,i,n_tmp, istat, ior_new, inode2, imesh2, n, ior
      Integer :: new_ffield_i, color_index, i_target, inode
      Character(60) :: TO_NODE
      Character(26) :: new_ffield_name
      Logical :: keep_xy
      Type (CORR_LL_TYPE), Pointer :: Now_LL, Tmp_LL, Next_LL, &
           Prev_LL
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
                  HR%Tau  = HR%Tau
                  HR%Mass = HR%Mass
               End If
            Else
               fed_max_alive = Max(fed_max_alive,HR%IntDose)
            End If
            fed_max = Max(fed_max,HR%IntDose)

            ! Calculate Purser's fractional effective dose (FED)
            If (T > 0.0_EB) Then
               If ( PCX%FED_MESH2 > 0 ) Then
                  x_int = Min(1.0_EB,Max(0.0_EB,(Now_LL%T_out-T)) / &
                       (Now_LL%T_out - Now_LL%T_in))
               Else
                  x_int = 1.0_EB
               End If
               HR%IntDose = DTSP*( (1.0_EB-x_int)*PCX%FED_CO_CO2_O2(2) + &
                    x_int*PCX%FED_CO_CO2_O2(1) ) + HR%IntDose
            End If

            If ( (Now_LL%T_out) <= T) Then
               Call Check_Target_Node(inode,inode2,istat,xx,yy,zz,ior_new, &
                    imesh2,T,new_ffield_name,new_ffield_i,color_index,&
                    i_target,keep_xy,angle)

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
                        If (EVAC_Node_List(n)%Mesh_Index == imesh2) &
                             HR%INODE = n
                     End Do
                     HR%NODE_NAME   = Trim(MESH_NAME(imesh2))
                     HR%FFIELD_NAME = Trim(new_ffield_name)
                     HR%I_FFIELD = new_ffield_i
                     HR%I_Target = I_Target

                     If (EVAC_Node_List(INODE2)%Node_Type == 'Door') Then
                        ior = - EVAC_DOORS( &
                             EVAC_Node_List(INODE2)%Node_Index)%IOR
                     End If
                     If (EVAC_Node_List(INODE2)%Node_Type == 'Entry') Then
                        ior = EVAC_ENTRYS( &
                             EVAC_Node_List(INODE2)%Node_Index)%IOR
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

                     If (MESHES(imesh2)%N_HUMANS+1 > &
                          MESHES(imesh2)%N_HUMANS_DIM) Then
                        ! Re-allocation is not yet checked.
                        Call SHUTDOWN('ERROR: Humans: no re-allocation yet')
                        Call RE_ALLOCATE_HUMANS(1,imesh2)
                     End If
                     MESHES(imesh2)%N_HUMANS = MESHES(imesh2)%N_HUMANS + 1
                     MESHES(imesh2)%HUMAN(MESHES(imesh2)%N_HUMANS) = HR
                     Write (LU_ERR,fmt='(a,i6,a,f8.2,a,a,a,f8.4)') &
                          ' EVAC: Person n:o', &
                          HR%ILABEL, ' out at ', T, &
                          ' s, corr ', Trim(PCX%ID_NAME), ' fed ', HR%IntDose
                  End If            ! target is door or entry

                  If (ior_new == 3) Then ! corridor
                     Write (LU_ERR,fmt='(a,i6,a,f8.2,a,a,f8.4)') &
                          ' EVAC: Person n:o', &
                          HR%ILABEL, ' change corr ', T, &
                          ' s, corr ', Trim(PCX%ID_NAME), HR%IntDose
                  End If

                  If (ior_new == 5) Then ! exit
                     Write (LU_ERR,fmt='(a,i6,a,f8.2,a,a,f8.4)') &
                          ' EVAC: Person n:o', &
                          HR%ILABEL, ' exits ', T, &
                          ' s, corr ', Trim(PCX%ID_NAME), HR%IntDose
                  End If
               Else
                  ! Can not move to the next node, so do not allow to move inside
                  ! the door ==> keep old position and put velocity to zero.
                  HR%X = HR%X_old
                  HR%Y = HR%Y_old
                  HR%I_Target = 0
                  HR%U = 0.0_EB
                  HR%V = 0.0_EB
                  HR%IOR = 0
               End If         ! istat == 0

            Else             ! T_out > T, i.e., still in corridor
               HR%X = HR%X_old
               HR%Y = HR%Y_old
               HR%IOR = 0
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
         imesh2,T,new_ffield_name,new_ffield_i,color_index,i_target,keep_xy,angle)
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
      Integer, Intent(in) :: inode, inode2
      Integer, Intent(out) :: istat, ior_new, imesh2, color_index, i_target
      Real(EB), Intent(out) :: xx,yy,zz, angle
      Real(EB), Intent(in) :: T
      Integer, Intent(inout) :: new_ffield_i
      Logical, Intent(in) :: keep_xy
      Character(26), Intent(inout) :: new_ffield_name
      Real(EB) RN, x1, x2, y1, y2, z1, z2, d_max, dist, Width, &
           xx1,yy1, max_fed, ave_K
      Integer  II, JJ, KK, ior, irnmax, irn, ie, izero, j1
      Real(EB), Dimension(6) :: r_tmp, x_tmp, y_tmp
      Real(EB), Dimension(:), Allocatable :: K_ave_Door
      Real(EB), Dimension(:), Allocatable :: FED_max_Door
      Integer, Dimension(:), Allocatable :: Color_Tmp
      Logical, Dimension(:), Allocatable :: Is_Known_Door, &
           Is_Visible_Door
      Integer :: i_tmp, i_tim, iii, jjj
      Logical :: PP_see_door
      Type (CORR_LL_TYPE), Pointer :: TmpCurrent, TmpLoop
      !
      xx = 0.0_EB ; yy = 0.0_EB ; zz = 0.0_EB 
      I_Target = 0
      Select Case (EVAC_Node_List(INODE2)%Node_Type)
      Case ('Door','Entry')
         ior_new = 1
         If (EVAC_Node_List(INODE2)%Node_Type == 'Door') Then
            PDX2 => EVAC_DOORS(EVAC_Node_List(INODE2)%Node_Index)
            imesh2  = PDX2%IMESH
            new_ffield_name = Trim(PDX2%GRID_NAME)
            Irn_Loop1: Do irn = 1, nmeshes
               If ( evacuation_only(irn) .And. Trim(PDX2%GRID_NAME) == &
                    Trim(MESH_NAME(irn)) ) Then
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
               If ( evacuation_only(irn) .And. Trim(PNX2%GRID_NAME) == &
                    Trim(MESH_NAME(irn)) ) Then
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
         End If
         If (Abs(ior) == 1) irnmax = Int(Width*4.0_EB)
         If (Abs(ior) == 2) irnmax = Int(Width*4.0_EB)
         If (Abs(ior) == 3 .Or. ior == 0) irnmax = &
              Int((x2-x1)*4.0_EB)*Int((y2-y1)*4.0_EB)
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
                  yy = 0.5_EB*(y1+y2) + (rn-0.5_EB)* &
                       Max(0.0_EB,Width-2.0_EB*HR%Radius-2.0_EB*HR%B)
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
                  xx = 0.5_EB*(x1+x2) + (rn-0.5_EB)* &
                       Max(0.0_EB,Width-2.0_EB*HR%Radius-2.0_EB*HR%B)
                  yy = y1 + (ior/Abs(ior))*(1.0_EB*HR%B + HR%r_torso)
               End If
            Case(0,3)
               Call Random_number(rn)
               yy = 0.5_EB*(y1+y2) + (rn-0.5_EB)* &
                    Max(0.0_EB,(y2-y1)-2.0_EB*HR%Radius-2.0_EB*HR%B)
               Call Random_number(rn)
               xx = 0.5_EB*(x1+x2) + (rn-0.5_EB)* &
                    Max(0.0_EB,(x2-x1)-2.0_EB*HR%Radius-2.0_EB*HR%B)
            End Select
            zz = 0.5_EB*(z1+z2)
            II = Floor( MFF%CELLSI(Floor((xx-MFF%XS)*MFF%RDXINT)) + 1.0_EB )
            JJ = Floor( MFF%CELLSJ(Floor((yy-MFF%YS)*MFF%RDYINT)) + 1.0_EB )
            KK = Floor( MFF%CELLSK(Floor((zz-MFF%ZS)*MFF%RDZINT)) + 1.0_EB )

            irn = irn + 1

            If (MFF%SOLID(MFF%CELL_INDEX(II,JJ,KK))) Cycle CheckPPForce
            If ( Abs(ior) == 2 .And. .Not. keep_xy ) Then
               xx1 = xx - HR%Radius - 1.0_EB*HR%B
               II = Floor(MFF%CELLSI(Floor((xx1-MFF%XS)*MFF%RDXINT))+1.0_EB)
               If (MFF%SOLID(MFF%CELL_INDEX(II,JJ,KK))) Cycle CheckPPForce
               xx1 = xx + HR%Radius + 1.0_EB*HR%B
               II = Floor(MFF%CELLSI(Floor((xx1-MFF%XS)*MFF%RDXINT))+1.0_EB)
               If (MFF%SOLID(MFF%CELL_INDEX(II,JJ,KK))) Cycle CheckPPForce
            End If
            If ( Abs(ior) == 1 .And. .Not. keep_xy ) Then
               yy1 = yy - HR%Radius - 1.0_EB*HR%B
               JJ = Floor(MFF%CELLSJ(Floor((yy1-MFF%YS)*MFF%RDYINT))+1.0_EB)
               If (MFF%SOLID(MFF%CELL_INDEX(II,JJ,KK))) Cycle CheckPPForce
               yy1 = yy + HR%Radius + 1.0_EB*HR%B
               JJ = Floor(MFF%CELLSJ(Floor((yy1-MFF%YS)*MFF%RDYINT))+1.0_EB)
               If (MFF%SOLID(MFF%CELL_INDEX(II,JJ,KK))) Cycle CheckPPForce
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
                     DIST = Sqrt((x_tmp(jjj)-x_tmp(iii))**2 + &
                          (y_tmp(jjj)-y_tmp(iii))**2) - &
                          (r_tmp(jjj)+r_tmp(iii))
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
               If (EVAC_Node_List(ie)%Mesh_index == imesh2 ) Then
                  i_tmp = ie
               End If
            End Do
            If (i_tmp == 0 .Or. i_tmp > n_egrids) Then
               Write(MESSAGE,'(A,I6)') &
                    'ERROR: Check_Target_Node, no imesh2 found ',imesh2
               Call SHUTDOWN(MESSAGE)
            End If

            If (j > 0 .And. Group_List(j)%GROUP_I_FFIELDS(i_tmp) > 0) Then
               ! There are already group members on this floor, use the same field
               new_ffield_i = Group_List(j)%GROUP_I_FFIELDS(i_tmp)
               new_ffield_name = Trim(MESH_NAME(new_ffield_i))
               If (j > 0) I_Target = Group_Known_Doors(j)%I_Target
            Else
               Allocate(Is_Known_Door(Max(1,n_doors+n_exits)),STAT=IZERO)
               Call ChkMemErr('Is_Known_Door','Color_Tmp',IZERO) 
               Allocate(Is_Visible_Door(Max(1,n_doors+n_exits)),STAT=IZERO)
               Call ChkMemErr('Is_Visible_Door','Color_Tmp',IZERO) 
               Allocate(FED_max_Door(Max(1,n_doors+n_exits)),STAT=IZERO)
               Call ChkMemErr('FED_max_Door','Color_Tmp',IZERO) 
               Allocate(K_ave_Door(Max(1,n_doors+n_exits)),STAT=IZERO)
               Call ChkMemErr('Evacuate_Humans','K_ave_Door',IZERO) 
               Allocate(Color_Tmp(Max(1,i33_dim)),STAT=IZERO)
               Call ChkMemErr('Check_Target_Node','Color_Tmp',IZERO) 
               K_ave_Door(:)      = 0.0_EB
               FED_max_Door(:)    = 0.0_EB
               Is_Known_Door(:)   = .False.
               Is_Visible_Door(:) = .False.
               If ( HR%GROUP_ID /= 0 ) Then
                  If ( HR%GROUP_ID < 0 ) Then
                     ! A lonely soul
                     Do ie = 1, n_doors
                        If ( EVAC_DOORS(ie)%IMESH == imesh2 ) Then
                           Is_Visible_Door(ie) = .True.
                           Do i_tim = 1, Human_Known_Doors(j1)%N_nodes
                              If (EVAC_DOORS(ie)%INODE == &
                                   Human_Known_Doors(j1)%I_nodes(i_tim)) &
                                   Is_Known_Door(ie) = .True.
                           End Do
                        End If
                     End Do
                     Do ie = 1, n_exits
                        If ( EVAC_EXITS(ie)%IMESH == imesh2 .And. &
                             .Not. EVAC_EXITS(ie)%COUNT_ONLY ) Then
                           Is_Visible_Door(n_doors+ie) = .True.
                           Do i_tim = 1, Human_Known_Doors(j1)%N_nodes
                              If (EVAC_EXITS(ie)%INODE == &
                                   Human_Known_Doors(j1)%I_nodes(i_tim)) &
                                   Is_Known_Door(n_doors+ie) = .True.
                           End Do
                        End If
                     End Do
                  Else
                     ! A member of a group
                     Do ie = 1, n_doors
                        If ( EVAC_DOORS(ie)%IMESH == imesh2 ) Then
                           Is_Visible_Door(ie) = .True.
                           Do i_tim = 1, Group_Known_Doors(j)%N_nodes
                              If (EVAC_DOORS(ie)%INODE == &
                                   Group_Known_Doors(j)%I_nodes(i_tim)) &
                                   Is_Known_Door(ie) = .True.
                           End Do
                        End If
                     End Do
                     Do ie = 1, n_exits
                        If ( EVAC_EXITS(ie)%IMESH == imesh2 .And. &
                             .Not. EVAC_EXITS(ie)%COUNT_ONLY ) Then
                           Is_Visible_Door(n_doors+ie) = .True.
                           Do i_tim = 1, Group_Known_Doors(j)%N_nodes
                              If (EVAC_EXITS(ie)%INODE == &
                                   Group_Known_Doors(j)%I_nodes(i_tim)) &
                                   Is_Known_Door(n_doors+ie) = .True.
                           End Do
                        End If
                     End Do
                  End If
               Else
                  Do ie = 1, n_doors
                     If ( EVAC_DOORS(ie)%IMESH == nm) Then
                        Is_Visible_Door(ie) = .True.
                     End If
                  End Do
                  Do ie = 1, n_exits
                     If ( EVAC_EXITS(ie)%IMESH == nm .And. &
                          .Not. EVAC_EXITS(ie)%COUNT_ONLY ) Then
                        Is_Visible_Door(n_doors+ie) = .True.
                     End If
                  End Do
                  Do i = 1, PNX%N_VENT_FFIELDS 
                     i_tim = 1
                     If (Trim(EVAC_Node_List(PNX%I_DOOR_NODES(ie) &
                          )%Node_Type) == 'Door') Then
                        i_tim = &
                             EVAC_Node_List(PNX%I_DOOR_NODES(ie))%Node_Index 
                     End If
                     If (Trim(EVAC_Node_List(PNX%I_DOOR_NODES(ie) &
                          )%Node_Type) == 'Exit' ) Then
                        i_tim = n_doors + &
                             EVAC_Node_List(PNX%I_DOOR_NODES(ie))%Node_Index 
                     End If
                     If ( Is_Visible_Door(i_tim) ) Then
                        ! Door/exit is on this floor
                        If (PNX%P_VENT_FFIELDS(ie) < 0.5_EB) Then
                           Is_Known_Door(i_tim) = .False.
                        Else
                           Is_Known_Door(i_tim) = .True.
                        End If
                     End If
                  End Do
               End If

               ! Find the visible doors (only the first member of a group)
               Do ie = 1, n_doors + n_exits
                  If ( Is_Visible_Door(ie) ) Then
                     If (EVAC_Node_List(n_egrids+n_entrys+ie)%Node_Type &
                          == 'Door' ) Then
                        XX1 = EVAC_DOORS(ie)%X 
                        YY1 = EVAC_DOORS(ie)%Y
                     Else            ! 'Exit'
                        XX1 = EVAC_EXITS(ie-n_doors)%X 
                        YY1 = EVAC_EXITS(ie-n_doors)%Y
                     End If
                     max_fed = 0.0_EB
                     ave_K   = 0.0_EB
                     II = Floor( MFF%CELLSI(Floor((XX1-MFF%XS)*MFF%RDXINT)) &
                          + 1.0_EB )
                     JJ = Floor( MFF%CELLSJ(Floor((YY1-MFF%YS)*MFF%RDYINT)) &
                          + 1.0_EB )
                     KK = Floor( MFF%CELLSK(Floor(( 0.5_EB*(MFF%ZF+MFF%ZS)-MFF%ZS)*MFF%RDZINT)) &
                          + 1.0_EB  )
                     XI  = MFF%CELLSI(Floor((xx-MFF%XS)*MFF%RDXINT))
                     YJ  = MFF%CELLSJ(Floor((yy-MFF%YS)*MFF%RDYINT))
                     IIE = Floor(XI+1.0_EB)
                     JJE = Floor(YJ+1.0_EB)
                     PP_see_door = .True. ! oletusarvo
                     If (Abs(XX-XX1) >= Abs(YY-YY1)) Then
                        If ( iie < ii) Then
                           iio = iie
                           jjo = jje
                           iie = ii
                           jje = jj
                           y_o = YY
                           x_o = XX
                           Delta_y = (YY1 - YY)
                        Else
                           Delta_y = (YY - YY1)
                           iio = ii
                           jjo = jj
                           y_o = YY1
                           x_o = XX1
                        End If
                        Delta_x = Abs(XX - XX1)
                        x_now = 0.0_EB
                        PP_see_door_x: Do iii = iio+1, iie-1
                           x_now = x_now + MFF%DX(iii)
                           y_now = y_o + x_now*(Delta_y/Delta_x)
                           jjj = Floor(MFF%CELLSJ(Floor((y_now-MFF%YS)* &
                                MFF%RDYINT))+1.0_EB)
                           ave_K = ave_K + MASS_EXTINCTION_COEFFICIENT* &
                                1.0E-6_EB*HUMAN_GRID(iii,jjj)%SOOT_DENS / &
                                ( iie-1 - (iio+1) + 1)
                           If (max_fed < &
                                MFF%HUMAN_GRID(iii,jjj)%FED_CO_CO2_O2) Then
                              max_fed = MFF%HUMAN_GRID(iii,jjj)%FED_CO_CO2_O2 
                           End If
                           tim_ic = MFF%CELL_INDEX(iii,jjj,KK)
                           If (MFF%SOLID(tim_ic)) Then
                              PP_see_door = .False.
                              Exit PP_see_door_x
                           End If
                        End Do PP_see_door_x
                     Else 
                        If ( jje < jj) Then
                           iio = iie
                           jjo = jje
                           iie = ii
                           jje = jj
                           y_o = YY
                           x_o = XX
                           Delta_x = (XX1 - XX)
                        Else
                           Delta_x = (XX - XX1)
                           iio = ii
                           jjo = jj
                           y_o = YY1
                           x_o = XX1
                        End If
                        Delta_y = Abs(YY - YY1)
                        y_now = 0.0_EB
                        PP_see_door_y: Do jjj = jjo+1, jje-1
                           y_now = y_now + MFF%DY(jjj)
                           x_now = x_o + y_now*(Delta_x/Delta_y)
                           iii = Floor(MFF%CELLSI(Floor((x_now-MFF%XS)* &
                                MFF%RDXINT))+1.0_EB)
                           ave_K = ave_K + MASS_EXTINCTION_COEFFICIENT* &
                                1.0E-6_EB*HUMAN_GRID(iii,jjj)%SOOT_DENS / &
                                ( jje-1 - (jjo+1) + 1)
                           If (max_fed < &
                                MFF%HUMAN_GRID(iii,jjj)%FED_CO_CO2_O2) Then
                              max_fed = MFF%HUMAN_GRID(iii,jjj)%FED_CO_CO2_O2 
                           End If
                           tim_ic = MFF%CELL_INDEX(iii,jjj,KK)
                           If (MFF%SOLID(tim_ic)) Then
                              PP_see_door = .False.
                              Exit PP_see_door_y
                           End If
                        End Do PP_see_door_y
                     End If

                     If (PP_see_door) Then
                        If (EVAC_Node_List(n_egrids+n_entrys+ie)%Node_Type &
                             == 'Door') Then
                           If (.Not. EVAC_DOORS(ie)%EXIT_SIGN) Then
                              Is_Visible_Door(ie) = .False.
                           End If
                        End If
                        FED_max_Door(ie) = max_fed
                        K_ave_Door(ie) = ave_K 
                     Else
                        Is_Visible_Door(ie) = .False.
                        iie = Floor(MFF%CELLSI(Floor((xx-MFF%XS)*MFF%RDXINT)) + 1.0_EB)
                        jje = Floor(MFF%CELLSJ(Floor((yy-MFF%YS)*MFF%RDYINT)) + 1.0_EB)
                        FED_max_Door(ie) = &
                             MFF%HUMAN_GRID(iie,jje)%FED_CO_CO2_O2 
                        K_ave_Door(ie) = MASS_EXTINCTION_COEFFICIENT* &
                             1.0E-6_EB*HUMAN_GRID(iie,jje)%SOOT_DENS
                     End If
                  End If            ! correct floor
               End Do              ! doors and exits


               Do ie = 1, n_doors
                  If ( EVAC_DOORS(ie)%TIME_OPEN > T .Or. EVAC_DOORS(ie)%TIME_CLOSE < T) Then
                     Is_Visible_Door(ie) = .False.
                     Is_Known_Door(ie) = .False.
                  End If
               End Do
               Do ie = 1, n_exits
                  If ( (EVAC_EXITS(ie)%TIME_OPEN > T .Or. EVAC_EXITS(ie)%TIME_CLOSE < T) .And. &
                       .Not. EVAC_EXITS(ie)%COUNT_ONLY ) Then
                     Is_Visible_Door(n_doors+ie) = .False.
                     Is_Known_Door(n_doors+ie) = .False.
                  End If
               End Do

               If (Any(Is_Known_Door) .Or. Any(Is_Visible_Door)) Then
                  irn = 0
                  d_max = Huge(d_max)
                  Do ie = 1, n_doors + n_exits
                     If (Is_Known_Door(ie) .And. &
                          Is_Visible_Door(ie) ) Then
                        x_o = 0.0_EB
                        y_o = 0.0_EB
                        If (Trim(EVAC_Node_List(n_egrids+n_entrys+ie &
                             )%Node_Type) == 'Door' ) Then
                           x_o = EVAC_DOORS( EVAC_Node_List( &
                                ie+n_egrids+n_entrys)%Node_Index )%X
                           y_o = EVAC_DOORS( EVAC_Node_List( &
                                ie+n_egrids+n_entrys)%Node_Index )%Y
                        Else          ! 'Exit'
                           x_o = EVAC_EXITS( EVAC_Node_List( &
                                ie+n_egrids+n_entrys)%Node_Index )%X
                           y_o = EVAC_EXITS( EVAC_Node_List( &
                                ie+n_egrids+n_entrys)%Node_Index )%Y
                        End If
                        If (FED_DOOR_CRIT > 0.0_EB) Then
                           dist = FED_max_Door(ie) * Sqrt((xx-x_o)**2 + &
                                (yy-y_o)**2)/HR%Speed
                        Else
                           dist = K_ave_Door(ie)
                        End If
                        If (( (x_o-xx)**2 + (y_o-yy)**2 ) < d_max &
                             .And. dist < Abs(FED_DOOR_CRIT)) Then
                           d_max = (x_o-xx)**2 + (y_o-yy)**2 
                           d_max = Max(0.0_EB,d_max)
                           irn = ie
                        End If
                     End If
                  End Do
                  If (irn > 0 ) Then
                     ! Known and visible door, no smoke
                     If (EVAC_Node_List(n_egrids+n_entrys+irn)%Node_Type &
                          == 'Door' ) Then
                        new_ffield_name = &
                             Trim(EVAC_DOORS(irn)%VENT_FFIELD)
                        new_ffield_i = &
                             EVAC_DOORS(irn)%I_VENT_FFIELD
                     Else            ! 'Exit'
                        new_ffield_name = &
                             Trim(EVAC_EXITS(irn-n_doors)%VENT_FFIELD)
                        new_ffield_i = &
                             EVAC_EXITS(irn-n_doors)%I_VENT_FFIELD
                        color_index = 0
                     End If
                  Else
                     ! No visible known door available, try non-visible known doors
                     irn   = 0
                     d_max = Huge(d_max)
                     Do ie = 1, n_doors + n_exits
                        If ( Is_Known_Door(ie) .And. &
                             .Not. Is_Visible_Door(ie) ) Then
                           x_o = 0.0_EB
                           y_o = 0.0_EB
                           If (EVAC_Node_List(n_egrids+n_entrys+ie)%Node_Type &
                                == 'Door' ) Then
                              x_o = EVAC_DOORS( EVAC_Node_List( &
                                   ie+n_egrids+n_entrys)%Node_Index )%X
                              y_o = EVAC_DOORS( EVAC_Node_List( &
                                   ie+n_egrids+n_entrys)%Node_Index )%Y
                           Else        ! 'Exit'
                              x_o = EVAC_EXITS( EVAC_Node_List( &
                                   ie+n_egrids+n_entrys)%Node_Index )%X
                              y_o = EVAC_EXITS( EVAC_Node_List( &
                                   ie+n_egrids+n_entrys)%Node_Index )%Y
                           End If
                           If (FED_DOOR_CRIT > 0.0_EB) Then
                              dist = FED_max_Door(ie) * Sqrt((xx-x_o)**2 + &
                                   (yy-y_o)**2)/HR%Speed
                           Else
                              dist = K_ave_Door(ie)
                           End If
                           If ( ((x_o-xx)**2 + (y_o-yy)**2) < d_max &
                                .And. dist < Abs(FED_DOOR_CRIT)) Then
                              d_max = (x_o-xx)**2 + (y_o-yy)**2 
                              d_max = Max(0.0_EB,d_max)
                              irn = ie
                           End If
                        End If
                     End Do
                     If (irn > 0 ) Then
                        !    Non-visible known door, no smoke
                        If (EVAC_Node_List( &
                             n_egrids+n_entrys+irn)%Node_Type &
                             == 'Door' ) Then
                           new_ffield_name = &
                                Trim(EVAC_DOORS(irn)%VENT_FFIELD)
                           new_ffield_i = &
                                EVAC_DOORS(irn)%I_VENT_FFIELD
                        Else          ! 'Exit'
                           new_ffield_name = &
                                Trim(EVAC_EXITS(irn-n_doors)%VENT_FFIELD)
                           new_ffield_i = &
                                EVAC_EXITS(irn-n_doors)%I_VENT_FFIELD
                        End If
                        color_index = 1
                     Else
                        ! known doors with no smoke have not been found
                        irn   = 0
                        d_max = Huge(d_max)
                        Do ie = 1, n_doors + n_exits
                           If (Is_Visible_Door(ie)) Then
                              x_o = 0.0_EB
                              y_o = 0.0_EB
                              If (EVAC_Node_List(n_egrids+n_entrys+ie &
                                   )%Node_Type == 'Door' ) Then
                                 x_o = EVAC_DOORS( EVAC_Node_List( &
                                      ie+n_egrids+n_entrys)%Node_Index )%X
                                 y_o = EVAC_DOORS( EVAC_Node_List( &
                                      ie+n_egrids+n_entrys)%Node_Index )%Y
                              Else      ! 'Exit'
                                 x_o = EVAC_EXITS( EVAC_Node_List( &
                                      ie+n_egrids+n_entrys)%Node_Index )%X
                                 y_o = EVAC_EXITS( EVAC_Node_List( &
                                      ie+n_egrids+n_entrys)%Node_Index )%Y
                              End If
                              If (FED_DOOR_CRIT > 0.0_EB) Then
                                 dist = FED_max_Door(ie) * Sqrt((xx-x_o)**2 + &
                                      (yy-y_o)**2)/HR%Speed
                              Else
                                 dist = K_ave_Door(ie)
                              End If
                              If ( ((x_o-xx)**2 + (y_o-yy)**2) < d_max &
                                   .And. dist < Abs(FED_DOOR_CRIT)) Then
                                 d_max = (x_o-xx)**2 + (y_o-yy)**2 
                                 d_max = Max(0.0_EB,d_max)
                                 irn = ie
                              End If
                           End If
                        End Do
                        If (irn > 0 ) Then
                           ! No smoke, visible door (not known)
                           If (EVAC_Node_List( &
                                n_egrids+n_entrys+irn)%Node_Type &
                                == 'Door' ) Then
                              new_ffield_name = &
                                   Trim(EVAC_DOORS(irn)%VENT_FFIELD)
                              new_ffield_i = &
                                   EVAC_DOORS(irn)%I_VENT_FFIELD
                           Else        ! 'Exit'
                              new_ffield_name = &
                                   Trim(EVAC_EXITS(irn-n_doors)%VENT_FFIELD)
                              new_ffield_i = &
                                   EVAC_EXITS(irn-n_doors)%I_VENT_FFIELD
                           End If
                           color_index = 2
                        Else
                           ! Now we have smoke and some visible or known doors
                           irn   = 0
                           d_max = Huge(d_max)
                           Do ie = 1, n_doors + n_exits
                              If ( Is_Known_Door(ie) .Or. &
                                   Is_Visible_Door(ie) ) Then
                                 x_o = 0.0_EB
                                 y_o = 0.0_EB
                                 If (EVAC_Node_List(n_egrids+n_entrys+ &
                                      ie)%Node_Type == 'Door' ) Then
                                    x_o = EVAC_DOORS( EVAC_Node_List(ie+ &
                                         n_egrids+n_entrys)%Node_Index )%X
                                    y_o = EVAC_DOORS( EVAC_Node_List(ie+ &
                                         n_egrids+n_entrys)%Node_Index )%Y
                                 Else    ! 'Exit'
                                    x_o = EVAC_EXITS( EVAC_Node_List(ie+ &
                                         n_egrids+n_entrys)%Node_Index )%X
                                    y_o = EVAC_EXITS( EVAC_Node_List(ie+ &
                                         n_egrids+n_entrys)%Node_Index )%Y
                                 End If
                                 If (FED_DOOR_CRIT > 0.0_EB) Then
                                    dist = HR%IntDose + FED_max_Door(ie) * &
                                         Sqrt((xx-x_o)**2 + (yy-y_o)**2)/ &
                                         HR%Speed
                                 Else
                                    ! Check that visibility > 0.5*distance to the door
                                    dist = (Sqrt((xx-x_o)**2 &
                                         + (yy-y_o)**2)*0.5_EB) / &
                                         (3.0_EB/K_ave_Door(ie))
                                 End If
                                 If (dist < d_max) Then
                                    d_max = dist
                                    irn = ie
                                 End If
                              End If
                           End Do

                           If (irn > 0 .And. d_max < 1.0_EB ) Then
                              ! Not too much smoke, i.e., non-lethal amount (known or visible doors)
                              If (EVAC_Node_List( &
                                   n_egrids+n_entrys+irn)%Node_Type &
                                   == 'Door' ) Then
                                 new_ffield_name = &
                                      Trim(EVAC_DOORS(irn)%VENT_FFIELD)
                                 new_ffield_i = &
                                      EVAC_DOORS(irn)%I_VENT_FFIELD
                              Else      ! 'Exit'
                                 new_ffield_name = &
                                      Trim(EVAC_EXITS(irn-n_doors)%VENT_FFIELD)
                                 new_ffield_i = &
                                      EVAC_EXITS(irn-n_doors)%I_VENT_FFIELD
                              End If
                              If (Is_Known_Door(irn) .And. &
                                   Is_Visible_Door(irn)) color_index = 3
                              If (Is_Known_Door(irn) .And. .Not. &
                                   Is_Visible_Door(irn)) color_index = 4
                              If (.Not. Is_Known_Door(irn) .And. &
                                   Is_Visible_Door(irn)) color_index = 5
                           Else        ! not case 1,2,3, i.e., no non-lethal door found
                              ! No door found, use the main evac grid ffield 
                              irn = 0
                              new_ffield_i    = imesh2
                              new_ffield_name = Trim(MESH_NAME(new_ffield_i))
                              color_index = 6
                           End If      ! case 4
                        End If        ! case 3
                     End If          ! case 2
                  End If            ! case 1
                  If (Color_Method == 4 ) Then
                     color_index =  6 ! default, cyan
                     If (irn > 0 .And. irn <= n_doors ) &
                          color_index = EVAC_DOORS(irn)%COLOR_INDEX
                     If (irn > n_doors .And. irn <= n_doors + n_exits) &
                          color_index = EVAC_EXITS(irn-n_doors)%COLOR_INDEX
                  End If
                  I_Target = irn
                  If (irn > 0 .And. .Not. Is_Visible_Door(Max(1,irn)) ) Then
                     ! I_Target >0: visible, <0: not visible
                     I_Target = -irn
                  End If

               Else  ! no known/visible door
                  I_Target = 0
                  new_ffield_i    = imesh2
                  new_ffield_name = Trim(MESH_NAME(new_ffield_i))
               End If
               If ( j > 0 ) Then
                  Group_List(j)%GROUP_I_FFIELDS(i_tmp) = new_ffield_i
               End If
               Deallocate(Color_Tmp)
               Deallocate(K_ave_Door)
               Deallocate(FED_max_Door)
               Deallocate(Is_Visible_Door)
               Deallocate(Is_Known_Door)
            End If  ! first member of a group?

         End If ! istat=0, i.e., put human to a new node

      Case ('Corr')
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
            TmpCurrent%T_out  = T + &
                 PCX2%Eff_Length/(PCX2%Fac_Speed*HR%v0_fac*HR%Speed)
            TmpCurrent%From1_To2 = .True.
            TmpCurrent%Index     = PCX2%n_inside
            TmpCurrent%Next      => PCX2%First
            PCX2%First            => TmpCurrent
         Else
            ior_new = 3     ! target is a corridor
            istat = 1       ! do not enter the corridor
            imesh2 = HR%IMESH
         End If
      Case ('Floor')
         ior_new = 4  ! target is a floor (NOT POSSIBLE, use entry)
         istat = 1         ! do not remove from the floor
         imesh2 = HR%IMESH 
      Case ('Exit')
         ior_new = 5  ! target is an exit
         istat = 0    ! remove from the floor (exit is always free)
         imesh2 = HR%IMESH
      Case Default
         ior_new = 6  ! target is not defined
         istat = 1    ! do not remove from the floor
         imesh2 = HR%IMESH
      End Select
      !
    End Subroutine Check_Target_Node
    !
    Subroutine REMOVE_PERSON(I)
      Implicit None
      !
      ! Remove a person
      !
      Integer, Intent(IN) :: I
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
      Integer :: IKILL,I,NM
      Real(EB) :: T
      !
      IKILL = 0
      DROP_LOOP: Do I=1,N_HUMANS
         !
         HR=>HUMAN(I)
         If (I > N_HUMANS-IKILL) Exit DROP_LOOP
         If (HR%X > XS .And. HR%X < XF .And. &
              HR%Y > YS .And. HR%Y < YF) &
              Cycle DROP_LOOP
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
      Real(EB), Intent(IN) :: Tin
      Integer,  Intent(IN) :: I_entry, NM
      Integer, Intent(OUT) :: istat
      Real(EB) RN, x1, x2, y1, y2, z1, z2, d_max, dist, &
           xx, yy, zz, xx1, yy1
      Integer  II, JJ, KK, ior, irnmax, irn, ie
      Real(EB), Dimension(6) ::y_tmp, x_tmp, r_tmp
      !
      istat = 1
      PNX => EVAC_ENTRYS(I_entry)
      If (PNX%IMESH /= NM ) Return
      If (PNX%Flow <= 0.0_EB ) Return
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
      Call CLASS_PROPERTIES
      HR%Tpre = 0.0_EB
      HR%IPC  = PNX%IPC
      HR%IEL  = -I_entry
      HR%GROUP_ID = 0
      HR%I_Target = 0
      !
      If (Abs(ior) == 1) irnmax = Int(PNX%Width*4.0_EB)
      If (Abs(ior) == 2) irnmax = Int(PNX%Width*4.0_EB)
      If (Abs(ior) == 3 .Or. ior == 0) irnmax = &
           Int((x2-x1)*4.0_EB)*Int((y2-y1)*4.0_EB)
      irnmax = Max(irnmax,5)
      !
      irn = 0
      CheckPPForce: Do While (irn < irnmax)
         Select Case (ior)
         Case(-1,1)
            Call Random_number(rn)
            yy = 0.5_EB*(y1+y2) + (rn-0.5_EB)* &
                 Max(0.0_EB,PNX%Width-2.0_EB*HR%Radius-2.0_EB*HR%B)
            xx = x1 + ior*5.0_EB*HR%B
            HR%Angle = (1-ior)*Pi/2.0_EB  ! ior=1: 0,  ior=-1: pi
         Case(-2,2)
            Call Random_number(rn)
            xx = 0.5_EB*(x1+x2) + (rn-0.5_EB)* &
                 Max(0.0_EB,PNX%Width-2.0_EB*HR%Radius-2.0_EB*HR%B)
            yy = y1 + (ior/Abs(ior))*5.0_EB*HR%B
            HR%Angle = Pi/2.0_EB + (2-ior)*Pi/4.0_EB  ! ior=2: (3/2)pi,  ior=-2: pi/2
         Case(0,3)
            Call Random_number(rn)
            yy = 0.5_EB*(y1+y2) + (rn-0.5_EB)* &
                 Max(0.0_EB,(y2-y1)-2.0_EB*HR%Radius-2.0_EB*HR%B)
            Call Random_number(rn)
            xx = 0.5_EB*(x1+x2) + (rn-0.5_EB)* &
                 Max(0.0_EB,(x2-x1)-2.0_EB*HR%Radius-2.0_EB*HR%B)
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
         KK = Floor( CELLSK(Floor((zz-ZS)*RDZINT)) + 1.0_EB )

         irn = irn + 1

         If (SOLID(CELL_INDEX(II,JJ,KK))) Cycle CheckPPForce
         If ( Abs(ior) == 2 ) Then
            xx1 = xx - HR%Radius - 1.0_EB*HR%B
            II = Floor(CELLSI(Floor((xx1-XS)*RDXINT))+1.0_EB)
            If (SOLID(CELL_INDEX(II,JJ,KK))) Cycle CheckPPForce
            xx1 = xx + HR%Radius + 1.0_EB*HR%B
            II = Floor(CELLSI(Floor((xx1-XS)*RDXINT))+1.0_EB)
            If (SOLID(CELL_INDEX(II,JJ,KK))) Cycle CheckPPForce
         End If
         If ( Abs(ior) == 1 ) Then
            yy1 = yy - HR%Radius - 1.0_EB*HR%B
            JJ = Floor(CELLSJ(Floor((yy1-YS)*RDYINT))+1.0_EB)
            If (SOLID(CELL_INDEX(II,JJ,KK))) Cycle CheckPPForce
            yy1 = yy + HR%Radius + 1.0_EB*HR%B
            JJ = Floor(CELLSJ(Floor((yy1-YS)*RDYINT))+1.0_EB)
            If (SOLID(CELL_INDEX(II,JJ,KK))) Cycle CheckPPForce
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
                  DIST = Sqrt((x_tmp(jjj)-x_tmp(iii))**2 + &
                       (y_tmp(jjj)-y_tmp(iii))**2) - &
                       (r_tmp(jjj)+r_tmp(iii))
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
         If (PNX%T_first <= 0.0_EB) PNX%T_first = Tin
         HR%X = xx
         HR%Y = yy
         HR%Z = zz
         ! 
         ILABEL_last = ILABEL_last + 1
         HR%ILABEL = ILABEL_last
         HR%SHOW = .True.    
         HR%COLOR_INDEX = 0
         Select Case (COLOR_METHOD)
         Case (-1)
            HR%COLOR_INDEX = 0
         Case (0)
            HR%COLOR_INDEX = PNX%COLOR_INDEX
         Case (1,2)
            HR%COLOR_INDEX = 0    ! lonely human
         Case (4)
            HR%COLOR_INDEX = PNX%COLOR_INDEX
         Case (5)
            ! Correct color is put, where the flow fields are chosen.
            HR%COLOR_INDEX = 0
         Case Default
            Write(MESSAGE,'(A,I3,A)') &
                 'ERROR: ENTRY_HUMAN COLOR METHOD',COLOR_METHOD, &
                 ' is not defined'
            Call SHUTDOWN(MESSAGE)
         End Select
         HR%FFIELD_NAME = Trim(PNX%GRID_NAME)
         HR%I_FFIELD    = 0
         Mesh2Loop: Do i = 1, nmeshes
            If ( evacuation_only(i) .And. Trim(HR%FFIELD_NAME) == &
                 Trim(MESH_NAME(i)) ) Then
               HR%I_FFIELD = i
               Exit Mesh2Loop
            End If
         End Do Mesh2Loop
         If ( HR%I_FFIELD == 0 ) Then
            Write(MESSAGE,'(A,A,A,A)') &
                 'ERROR: ENTR line ',Trim(PNX%ID_NAME), &
                 ' problem with flow field name, ', &
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
         Write (LU_ERR,fmt='(a,i6,a,f8.2,a,3a)') &
              ' EVAC: Person n:o', &
              HR%ILABEL, ' inserted ', Tin, &
              ' s, entry ', Trim(PNX%ID_NAME),' ffield ', &
              Trim(HR%FFIELD_NAME)
      End If
      ! 
    End Subroutine ENTRY_HUMAN
    !
    !
  End Subroutine EVACUATE_HUMANS
! ============================================================
! EVACUATE_HUMANS ENDS HERE.
! ============================================================
!
! ============================================================
! NEXT ARE MODULE SUBPROGRMAS
! ============================================================
!
  Subroutine CLASS_PROPERTIES
    Implicit None
    !
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
       n_par = 2
       Randomtype = 5
       RandomPara(1) = PCP%V_mean
       RandomPara(2) = PCP%V_para
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
       ! Parameters: (ave,alpha,lambda)
       n_par = 3
       Randomtype = 8
       RandomPara(1) = PCP%V_mean
       RandomPara(2) = PCP%V_para
       RandomPara(3) = PCP%V_para2
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Speed = rnd_vec(1)
    Case(9)   ! Gumbel
       ! Parameters: (ave,alpha)
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
       RandomPara(1) = 0.5_EB*(PCP%D_high+PCP%D_low)
       PCP%D_mean    = 0.5_EB*(PCP%D_high+PCP%D_low)
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
       n_par = 2
       Randomtype = 5
       RandomPara(1) = PCP%D_mean
       RandomPara(2) = PCP%D_para
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
       n_par = 2
       Randomtype = 5
       RandomPara(1) = PCP%Tau_mean
       RandomPara(2) = PCP%Tau_para
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
       n_par = 2
       Randomtype = 5
       RandomPara(1) = PCP%Tdet_mean
       RandomPara(2) = PCP%Tdet_para
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
       HR%Tpre  = PCP%Tpre_mean
    Case(1)   ! Uniform
       ! Parameters: (ave,min,max) ave not used
       n_par = 3
       Randomtype = 1
       RandomPara(1) = 0.5_EB*(PCP%Tpre_high+PCP%Tpre_low)
       RandomPara(2) = PCP%Tpre_low
       RandomPara(3) = PCP%Tpre_high
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tpre = rnd_vec(1)
    Case(2)   ! Truncated Normal
       ! Parameters: (ave,sigma,min,max)
       n_par = 4
       Randomtype = 2
       RandomPara(1) = PCP%Tpre_mean
       RandomPara(2) = PCP%Tpre_para
       RandomPara(3) = PCP%Tpre_low
       RandomPara(4) = PCP%Tpre_high
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tpre = rnd_vec(1)
    Case(3)   ! Gamma
       ! Parameters: (ave,alpha,beta) ave not used
       n_par = 3
       Randomtype = 3
       RandomPara(1) = PCP%Tpre_mean
       RandomPara(2) = PCP%Tpre_para
       RandomPara(3) = PCP%Tpre_para2
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tpre = rnd_vec(1)
    Case(4)   ! Normal
       ! Parameters: (ave,sigma)
       n_par = 2
       Randomtype = 4
       RandomPara(1) = PCP%Tpre_mean
       RandomPara(2) = PCP%Tpre_para
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tpre = rnd_vec(1)
    Case(5)   ! LogNormal
       ! mean and variance of log(x) should be given
       ! Parameters: (ave,sigma) of ln(x)
       n_par = 2
       Randomtype = 5
       RandomPara(1) = PCP%Tpre_mean
       RandomPara(2) = PCP%Tpre_para
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tpre = rnd_vec(1)
    Case(6)   ! Beta
       ! Parameters: (ave,a,b) ave not used
       n_par = 3
       Randomtype = 6
       RandomPara(1) = PCP%Tpre_mean
       RandomPara(2) = PCP%Tpre_para
       RandomPara(3) = PCP%Tpre_para2
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tpre = rnd_vec(1)
    Case(7)   ! Triangular
       ! Parameters: (peak,min,max)
       n_par = 3
       Randomtype = 7
       RandomPara(1) = PCP%Tpre_mean
       RandomPara(2) = PCP%Tpre_low
       RandomPara(3) = PCP%Tpre_high
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tpre = rnd_vec(1)
    Case(8)   ! Weibull  (alpha=1: Exponential)
       ! Parameters: (ave,alpha,lambda)
       n_par = 3
       Randomtype = 8
       RandomPara(1) = PCP%Tpre_mean
       RandomPara(2) = PCP%Tpre_para
       RandomPara(3) = PCP%Tpre_para2
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tpre = rnd_vec(1)
    Case(9)   ! Gumbel
       ! Parameters: (ave,alpha)
       n_par = 2
       Randomtype = 9
       RandomPara(1) = PCP%Tpre_mean
       RandomPara(2) = PCP%Tpre_para
       Call RandomNumbers(n_rnd, n_par, RandomType, RandomPara(1:n_par), rnd_vec)
       HR%Tpre = rnd_vec(1)
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
    Type (HUMAN_TYPE), Allocatable, Dimension(:) :: DUMMY
    Integer IZERO
    Integer, Intent(IN) :: CODE,NM
    Type (MESH_TYPE), Pointer :: M
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
    !
    Integer, Intent(IN) :: NM
    Real(EB), Intent(IN) :: T
    Integer :: NPP,NPLIM,i,izero,nn,n
    Real(EB) :: TNOW
    Real(FB), Allocatable, Dimension(:) :: XP,YP,ZP
    Real(FB), Allocatable, Dimension(:,:) :: QP, AP
    Integer, Allocatable, Dimension(:) :: TA
    !
    TNOW=SECOND() 
    !
    If (.Not.Any(EVACUATION_GRID)) Return
    If (.Not.(EVACUATION_ONLY(NM) .And. EVACUATION_GRID(NM))) Return
!!$    If (EVACUATION_ONLY(NM) .And. EVACUATION_GRID(NM)) Then
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
          If (HR%COLOR_INDEX < 0) HR%COLOR_INDEX = 0
          If (HR%COLOR_INDEX > 7) HR%COLOR_INDEX = 7
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

          DO NN=1,EVAC_N_QUANTITIES
             SELECT CASE(EVAC_QUANTITIES_INDEX(NN))
             CASE(240)  ! MOTIVE_ACCELERATION, Unimpeded walking Speed / tau
                QP(NPP,NN) = Real(HR%Speed/HR%Tau,FB)
             CASE(241)  ! FED_DOSE, Fractional Effective Dose
                QP(NPP,NN) = Real(HR%IntDose,FB)
             CASE(242)  ! SPEED, Human peed
                QP(NPP,NN) = Real(Sqrt(HR%U**2 + HR%V**2),FB)
             CASE(243)  ! ANGULAR_SPEED, Human Angular Velocity
                QP(NPP,NN) = Real(HR%Omega,FB)
             CASE(244)  ! ACCELERATION, Human acc.
                QP(NPP,NN) = Real(Sqrt( ((HR%UBAR-HR%U)/HR%Tau + (HR%F_X/HR%Mass))**2 + &
                     ((HR%VBAR-HR%V)/HR%Tau + (HR%F_Y/HR%Mass))**2 ),FB)
             CASE(245)  ! CONTACT_LINEFORCE, Human Pressure: contact forces
                QP(NPP,NN) = Real(HR%SumForces ,FB)
             CASE(246)  ! TOTAL_LINEFORCE, Human Pressure2: contact + social
                QP(NPP,NN) = Real(HR%SumForces2 ,FB)
             CASE(247)  ! COLOR, Human color index
                QP(NPP,NN) = Real(HR%COLOR_INDEX,FB)
             CASE(248)  ! ANGULAR_ACCELERATION, 
                QP(NPP,NN) = Real(HR%Torque/HR%M_iner,FB)
             END SELECT
          ENDDO

          If (NPP>=NPPS) Exit PLOOP
       End Do PLOOP
       !
       ! Dump human data into the .prt5 file
       !
       Write(LU_PART(NM)) NPLIM
       WRITE(LU_PART(NM)) (XP(I),I=1,NPLIM),(YP(I),I=1,NPLIM),(ZP(I),I=1,NPLIM), &
            (AP(I,1),I=1,NPLIM),(AP(I,2),I=1,NPLIM),(AP(I,3),I=1,NPLIM),(AP(I,4),I=1,NPLIM)
       WRITE(LU_PART(NM)) (TA(I),I=1,NPLIM)
       If (EVAC_N_QUANTITIES > 0) Then
          WRITE(LU_PART(NM)) ((QP(I,NN),I=1,NPLIM),NN=1,EVAC_N_QUANTITIES)
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

!!$    End If
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

    If ( (GaussFlag == 1)  .And. &
         (Abs(GaussSet2-gmean) <= gcutmult*Sqrt(gtheta)) ) Then
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

          If ( (Abs(GaussSet1-gmean) <= gcutmult*Sqrt(gtheta)) &
               .Or. (Abs(GaussSet2-gmean) <= &
               gcutmult*Sqrt(gtheta)) ) Then
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

    If ( (GTrunFlag == 1)  .And. (GTrunSet2 >= glow) &
         .And. (GTrunSet2 <= ghigh) ) Then
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

          If ( ((GTrunSet1 >= glow) .And. &
               (GTrunSet1 <= ghigh)) .Or. &
               ((GTrunSet2 >= glow) .And. &
               (GTrunSet2 <= ghigh)) ) Then
             Exit
          End If

       End Do

       If ( (GTrunSet1 >= glow) .And. &
            (GTrunSet1 <= ghigh) ) Then
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
    ! Dump human data to CHID_evac.csv
    !
    Real(EB), Intent(In) :: Tin
    Character(50) tcform
    Integer n_cols, n_tot_humans, i
    !
    If (.Not.Any(EVACUATION_GRID)) Return
    !
    ! Output first the floors then the corridors
    n_cols = n_egrids + n_corrs + n_exits + n_doors + 1
    n_tot_humans = 0
    Do i = 1, n_egrids
       n_tot_humans = n_tot_humans + &
            MESHES(EVAC_Node_List(i)%Mesh_index)%N_HUMANS
    End Do
    Do i = 1, n_corrs
       n_tot_humans = n_tot_humans + &
            EVAC_CORRS(i)%n_inside
    End Do
    !
    If (n_dead >= 0) Then
       ! Write the 'fed' columns
       Write(tcform,'(a,i4.4,a,a)') "(ES13.5E3,",n_cols+1+n_exits+n_doors, &
            "(',',i8)", ",',',ES13.5E3,',',ES13.5E3)"
       Write (LU_EVACCSV,fmt=tcform) Tin, n_tot_humans, &
            (MESHES(EVAC_Node_List(i)%Mesh_index)%N_HUMANS, &
            i=1,n_egrids), &
            (EVAC_CORRS(i)%n_inside, i = 1,n_corrs), &
            (EVAC_EXITS(i)%ICOUNT, i = 1,n_exits), &
            (EVAC_DOORS(i)%ICOUNT, i = 1,n_doors), &
            (EVAC_EXITS(i)%NTARGET, i = 1,n_exits), &
            (EVAC_DOORS(i)%NTARGET, i = 1,n_doors), &
            n_dead, fed_max, fed_max_alive
    Else
       ! Do not write the 'fed' columns
       Write(tcform,'(a,i4.4,a)') "(ES13.5E3,",n_cols+n_exits+n_doors, &
            "(',',i8))"
       Write (LU_EVACCSV,fmt=tcform) Tin, n_tot_humans, &
            (MESHES(EVAC_Node_List(i)%Mesh_index)%N_HUMANS, &
            i=1,n_egrids), &
            (EVAC_CORRS(i)%n_inside, i = 1,n_corrs), &
            (EVAC_EXITS(i)%ICOUNT, i = 1,n_exits), &
            (EVAC_DOORS(i)%ICOUNT, i = 1,n_doors), &
            (EVAC_EXITS(i)%NTARGET, i = 1,n_exits), &
            (EVAC_DOORS(i)%NTARGET, i = 1,n_doors)
    End If
    !
  End Subroutine DUMP_EVAC_CSV
!
SUBROUTINE GET_REV_evac(MODULE_REV,MODULE_DATE)
INTEGER,INTENT(INOUT) :: MODULE_REV
CHARACTER(255),INTENT(INOUT) :: MODULE_DATE

WRITE(MODULE_DATE,'(A)') evacrev(INDEX(evacrev,':')+1:LEN_TRIM(evacrev)-2)
READ (MODULE_DATE,'(I5)') MODULE_REV
WRITE(MODULE_DATE,'(A)') evacdate

END SUBROUTINE GET_REV_evac

End Module EVAC
