MODULE VEGE

USE COMP_FUNCTIONS
USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_POINTERS
USE TRAN
USE PART
USE MEMORY_FUNCTIONS, ONLY:CHKMEMERR
USE TYPES, ONLY: LAGRANGIAN_PARTICLE_TYPE, LAGRANGIAN_PARTICLE_CLASS_TYPE! WALL_TYPE,SURFACE_TYPE
IMPLICIT NONE
PRIVATE
PUBLIC INITIALIZE_LEVEL_SET_FIRESPREAD,LEVEL_SET_FIRESPREAD,LEVEL_SET_FIRE_SPREAD_NEW
!TYPE (WALL_TYPE), POINTER :: WC
!TYPE (SURFACE_TYPE), POINTER :: SF
INTEGER :: IZERO
!For Level Set
INTEGER  :: LIMITER_LS,LU_SLCF_LS,NX_LS,NY_LS
REAL(EB) :: DT_LS,DX_LS,DY_LS,DT_OUTPUT,SUM_T_SLCF,SUMTIME,TIME_LS,TIME_FLANKFIRE_QUENCH
REAL(EB) :: DT_COEF,DYN_SR_MAX,IDX_LS,IDY_LS,T_FINAL,ROS_HEAD1,UMAG,UMF_TMP
REAL(EB) :: CPUTIME,LS_T_BEG,ROS_BACKS,ROS_HEADS
REAL(EB) :: B_ROTH,BETA_OP_ROTH,C_ROTH,E_ROTH

CONTAINS


SUBROUTINE INITIALIZE_LEVEL_SET_FIRESPREAD(NM)

!--------------------------------------------------------------------------------------------------
!
! Level set based modeling of fire spread across terrain. Currently, no computation of the wind field
! is needed. Instead, U0 and V0 which are specified on the MISC line of the input file, are used for
! the wind field direction. Does use the extent of the vegetation as defined in the fds input file.
! Level ground spread rates for the head, flank, and back fires are user defined.
! Spread rate dependence on slope is according to McArthur's rules.
!
! Issues:
! 1) Need to make level set computation mesh dependent so the the LS slice file
!    is created only where fire is expected
! 2) Need to make computation capable of running with multiple processors
!
!
INTEGER, INTENT(IN) :: NM
INTEGER :: I,IM1,IM2,IIG,IP1,IP2,IW,J,JJG,JM1,JP1,KKG,KDUM,KWIND
REAL(EB) :: DZT_DUM,LX,SR_MAX,UMAX_LS,VMAX_LS
REAL(EB) :: G_EAST,G_WEST,G_SOUTH,G_NORTH

!---Vars for Farsite emulation (phi_w calculation)---
REAL(EB) :: PHX,PHY,MAG_PHI
REAL(EB) :: PHI_W_X,PHI_W_Y,MAG_PHI_S,UMF_X,UMF_Y

REAL(EB), ALLOCATABLE, DIMENSION(:) :: X_LS,Y_LS
!
!REAL(FB) :: MAG_SR
CHARACTER(30) :: SMOKEVIEW_LABEL,SMOKEVIEW_BAR_LABEL,UNITS

REAL(EB), POINTER, DIMENSION(:,:) :: ZT => NULL()

TYPE (WALL_TYPE),    POINTER :: WC =>NULL()
TYPE (SURFACE_TYPE), POINTER :: SF =>NULL()

CALL CPU_TIME(CPUTIME)
LS_T_BEG = CPUTIME

CALL POINT_TO_MESH(NM)

ZT => LS_Z_TERRAIN

!WRITE(LU_OUTPUT,*)'level set: z(*)',z
!WRITE(LU_OUTPUT,*)'level set: ls_z_terrain(1,1)',ls_z_terrain(:,100)
!
!-Initialize variables
!
!-- Domain specification (meters) from input file (assumes there's only one mesh)
!
 LX = XF - XS ; NX_LS = IBAR
 DX_LS = LX/REAL(NX_LS,EB) ; IDX_LS = 1.0_EB / DX_LS
 LX = YF - YS ; NY_LS = JBAR
 DY_LS = LX/REAL(NY_LS,EB) ; IDY_LS = 1.0_EB / DY_LS
 T_FINAL = T_END
!
! Define spread rate across domain (including no burn areas)
!
ALLOCATE(HEAD_WIDTH(NX_LS,NY_LS))  ; CALL ChkMemErr('VEGE:LEVEL SET','ROS_FLANK',IZERO)
ALLOCATE(ROS_HEAD(NX_LS,NY_LS))    ; CALL ChkMemErr('VEGE:LEVEL SET','ROS_HEAD',IZERO)
ALLOCATE(ROS_FLANK(NX_LS,NY_LS))   ; CALL ChkMemErr('VEGE:LEVEL SET','ROS_FLANK',IZERO)
ALLOCATE(ROS_BACKU(NX_LS,NY_LS))   ; CALL ChkMemErr('VEGE:LEVEL SET','ROS_BACKU',IZERO)
ALLOCATE(FLANKFIRE_LIFETIME(NX_LS,NY_LS))   ; CALL ChkMemErr('VEGE:LEVEL SET','FLANKFIRE_LIFETIME',IZERO)
ALLOCATE(WIND_EXP(NX_LS,NY_LS))   ; CALL ChkMemErr('VEGE:LEVEL SET','WIND_EXP',IZERO)

FLANKFIRE_LIFETIME = 0.0_EB !handles finite length (lifetime) flankfires. For
!                            quenching flanks with lifetimes > TIME_FLANKFIRE_QUENCH
TIME_FLANKFIRE_QUENCH = 20.0_EB !flankfire lifetime in seconds

! Assign spread rates (i.e., vegetation types) to locations on terrain

SUMTIME = 0.0_EB ! Used for time step output

SUM_T_SLCF = 0._EB
DT_LS = 0.1_EB
DT_COEF = 0.5_EB

ROS_HEAD  = 0.0_EB
ROS_HEAD1 = 0.0_EB !needed when dependence on head width is computed
ROS_FLANK = 0.0_EB
ROS_BACKU = 0.0_EB
WIND_EXP  = 1.0_EB

LSET_ELLIPSE = .FALSE. ! Flag for the elliptical spread model
LSET_TAN2    = .FALSE. ! Flag for ROS proportional to Tan(slope)^2
HEAD_WIDTH   = 1.0_EB

!WRITE(LU_OUTPUT,*)'surface ros',surface%veg_lset_ros_head
!WRITE(LU_OUTPUT,*)'surface wind_exp',surface%veg_lset_wind_exp
!
!C_F = 0.2_EB
!
! -- Flux limiter
!LIMITER_LS = 1 !MINMOD
!LIMITER_LS = 2 !SUPERBEE
!LIMITER_LS = 3 !First order upwinding
!
!
LIMITER_LS = I_FLUX_LIMITER
IF (LIMITER_LS > 3) LIMITER_LS = 1

! -- Output file
DT_OUTPUT  = 0.5_EB
IF (DT_SLCF > 0._EB) DT_OUTPUT = DT_SLCF

TIME_LS    = 0._EB
!LU_SLCF_LS = 9999
LU_SLCF_LS = GET_FILE_NUMBER()
SMOKEVIEW_LABEL = 'phifield'
SMOKEVIEW_BAR_LABEL = 'phifield'
UNITS  = 'C'
OPEN(LU_SLCF_LS,FILE=TRIM(CHID)//'_lsfs.sf',FORM='UNFORMATTED',STATUS='REPLACE')
WRITE(LU_SLCF_LS) SMOKEVIEW_LABEL(1:30)
WRITE(LU_SLCF_LS) SMOKEVIEW_LABEL(1:30)
WRITE(LU_SLCF_LS) UNITS(1:30)
WRITE(LU_SLCF_LS)1,NX_LS,1,NY_LS,1,1
!
!
! =============== end of case specifications ========================
!
!-- Allocate arrays

! Level set values (Phi)

ALLOCATE(PHI_LS(NX_LS,NY_LS)) ; CALL ChkMemErr('VEGE:LEVEL SET','PHI_LS',IZERO)
ALLOCATE(PHI0_LS(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','PHI0_LS',IZERO)
ALLOCATE(PHI1_LS(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','PHI1_LS',IZERO)
ALLOCATE(PHI_TEMP(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','PHI_TEMP',IZERO)
ALLOCATE(PHI_OUT(NX_LS,NY_LS)) ; CALL ChkMemErr('VEGE:LEVEL SET','PHI_OUT',IZERO)
PHI_OUT = 0.0_FB

! Wind speeds at ground cells in domain
ALLOCATE(U_LS(NX_LS,NY_LS)) ; CALL ChkMemErr('VEGE:LEVEL SET','U_LS',IZERO) ; U_LS = 0._EB
ALLOCATE(V_LS(NX_LS,NY_LS)) ; CALL ChkMemErr('VEGE:LEVEL SET','V_LS',IZERO) ; V_LS = 0._EB

ALLOCATE(FLUX0_LS(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','FLUX0_LS',IZERO)
ALLOCATE(FLUX1_LS(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','FLUX1_LS',IZERO)

! Slopes (gradients)
ALLOCATE(DZTDX(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','DZDTX',IZERO)
ALLOCATE(DZTDY(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','DZDTY',IZERO)
ALLOCATE(MAG_ZT(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','MAG_ZT',IZERO)

! |ROS|
ALLOCATE(MAG_SR_OUT(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','MAG_SR_OUT',IZERO)

!-- Computational grid
ALLOCATE(X_LS(NX_LS))   ; CALL ChkMemErr('VEGE:LEVEL SET','X_LS',IZERO)
ALLOCATE(Y_LS(NY_LS+1)) ; CALL ChkMemErr('VEGE:LEVEL SET','Y_LS',IZERO)

!'Time of arrival' grid 02dec11
ALLOCATE(TOA(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','LSTOA',IZERO)
TOA = -1.0_EB

!----------Rothermel 'Phi' factors for effects of Wind and Slope on ROS ----------
! Not to be confused with the level set value (Phi)
ALLOCATE(PHI_WS(NX_LS,NY_LS))   ; CALL ChkMemErr('VEGE:LEVEL SET','PHI_W',IZERO)
ALLOCATE(PHI_S(NX_LS,NY_LS))    ; CALL ChkMemErr('VEGE:LEVEL SET','PHI_S',IZERO)
ALLOCATE(PHI_S_X(NX_LS,NY_LS))  ; CALL ChkMemErr('VEGE:LEVEL SET','PHI_S_X',IZERO)
ALLOCATE(PHI_S_Y(NX_LS,NY_LS))  ; CALL ChkMemErr('VEGE:LEVEL SET','PHI_S_Y',IZERO)
ALLOCATE(PHI_W(NX_LS,NY_LS))    ; CALL ChkMemErr('VEGE:LEVEL SET','PHI_W',IZERO)
PHI_WS    = 0.0_EB

! UMF = wind speed at mean flame heights
ALLOCATE(UMF(NX_LS,NY_LS))    ; CALL ChkMemErr('VEGE:LEVEL SET','UMF',IZERO)
ALLOCATE(THETA_ELPS(NX_LS,NY_LS))    ; CALL ChkMemErr('VEGE:LEVEL SET','THETA_ELPS',IZERO)
THETA_ELPS   = 0.0_EB ! Normal to fireline
!-----------------------------------------------------------------------------------


! ROS in X and Y directions for elliptical model
ALLOCATE(SR_X_LS(NX_LS,NY_LS)) ; CALL ChkMemErr('VEGE:LEVEL SET','SR_X_LS',IZERO)
ALLOCATE(SR_Y_LS(NX_LS,NY_LS)) ; CALL ChkMemErr('VEGE:LEVEL SET','SR_Y_LS',IZERO)

! Aspect of terrain slope for elliptical model
ALLOCATE(ASPECT(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','ASPECT',IZERO)


DO I = 0,NX_LS-1
!X_LS(I+1) = -0.5_EB*LX + 0.5_EB*DX_LS + DX_LS*REAL(I,EB)
 X_LS(I+1) = XS + 0.5_EB*DX_LS + DX_LS*REAL(I,EB)
ENDDO
!
DO J = 0,NY_LS
 Y_LS(J+1) = YS + DY_LS*REAL(J,EB)
ENDDO

!---- Compute components of terrain slope gradient and magnitude of gradient

GRADIENT_ILOOP: DO I = 1,NX_LS
 IM1=I-1 ; IM2=I-2
 IP1=I+1 ; IP2=I+2

 IF (I==1) IM1 = I
 IF (I==NX_LS) IP1 = I

 DO J = 1,NY_LS
   JM1=J-1
   JP1=J+1

   IF (J==1) JM1 = J
   IF (J==NX_LS) JP1 = J

   !GIS-type slope calculation
   !Probably not needed, but left in for experimental purposes
   !G_EAST  = ZT(IP1,JP1) + 2._EB * ZT(IP1,J) + ZT(IP1,JM1)
   !G_WEST  = ZT(IM1,JP1) + 2._EB * ZT(IM1,J) + ZT(IM1,JM1)
   !G_NORTH = ZT(IM1,JP1) + 2._EB * ZT(I,JP1) + ZT(IP1,JP1)
   !G_SOUTH = ZT(IM1,JM1) + 2._EB * ZT(I,JM1) + ZT(IP1,JM1)
   !
   !DZTDX(I,J) = (G_EAST-G_WEST) / (8._EB*DX_LS)
   !DZTDY(I,J) = (G_NORTH-G_SOUTH) / (8._EB*DY_LS)

   G_EAST  = 0.5_EB*( ZT(I,J) + ZT(IP1,J) )
   G_WEST  = 0.5_EB*( ZT(I,J) + ZT(IM1,J) )
   G_NORTH = 0.5_EB*( ZT(I,J) + ZT(I,JP1) )
   G_SOUTH = 0.5_EB*( ZT(I,J) + ZT(I,JM1) )

   DZTDX(I,J) = (G_EAST-G_WEST) * IDX_LS
   DZTDY(I,J) = (G_NORTH-G_SOUTH) * IDY_LS


   MAG_ZT(I,J) = SQRT(DZTDX(I,J)**2 + DZTDY(I,J)**2)

   ASPECT(I,J) = PIO2 - ATAN2(-DZTDY(I,J),-DZTDX(I,J))
   IF (ASPECT(I,J) < 0.0_EB) ASPECT(I,J) = 2._EB*PI + ASPECT(I,J)

 ENDDO

ENDDO GRADIENT_ILOOP

!
!_____________________________________________________________________________
!
! Initialize arrays for head, flank, and back fire spread rates with values
! explicitly declared in the input file or from FARSITE head fire and ellipse
! based flank and back fires.
! Fill arrays for the horizontal component of the velocity arrays.
! Initialize level set scalar array PHI

! PHI_MIN_LS = -1._EB
! PHI_MAX_LS =  1._EB
PHI_LS     = PHI_MIN_LS

LSET_INIT_WALL_CELL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
  WC  => WALL(IW)
  IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE LSET_INIT_WALL_CELL_LOOP
  SF  => SURFACE(WC%SURF_INDEX)

  IIG = WC%ONE_D%IIG
  JJG = WC%ONE_D%JJG
  KKG = WC%ONE_D%KKG

! Ignite landscape at user specified location if ignition is at time zero
  IF (SF%VEG_LSET_IGNITE_T == 0.0_EB) PHI_LS(IIG,JJG) = PHI_MAX_LS

! Wind field
  U_LS(IIG,JJG) = U(IIG,JJG,KKG)
  V_LS(IIG,JJG) = V(IIG,JJG,KKG)
!WRITE(LU_OUTPUT,*)'veg: u_ls(i,j)',u_ls(iig,jjg)

  IF (.NOT. SF%VEG_LSET_SPREAD) CYCLE LSET_INIT_WALL_CELL_LOOP
!WRITE(LU_OUTPUT,*)'x,y,z and U,V',X(IIG),Y(JJG),Z(KKG),U(IIG,JJG,KKG),V(IIG,JJG,KKG)

  !Diagnostics
  !WRITE(LU_OUTPUT,*)'IIG,JJG',iig,jjg
  !WRITE(LU_OUTPUT,*)'ROS_HEAD',SF%VEG_LSET_ROS_HEAD
  !WRITE(LU_OUTPUT,*)'ROS_HEAD,ROS_FLANK,ROS_BACK',SF%VEG_LSET_ROS_HEAD,SF%VEG_LSET_ROS_FLANK,SF%VEG_LSET_ROS_BACK

  HEAD_WIDTH(IIG,JJG)= DX_LS

  UMAG     = SQRT(U_LS(IIG,JJG)**2 + V_LS(IIG,JJG)**2)
!AU grassland ROS for infinite head and 6% moisutre
! ROS_HEAD(IIG,JJG)  = SF%VEG_LSET_ROS_HEAD*(0.165_EB + 0.534_EB*UMAG)*0.523_EB

  ROS_HEAD(IIG,JJG)  = SF%VEG_LSET_ROS_HEAD
  ROS_FLANK(IIG,JJG) = SF%VEG_LSET_ROS_FLANK
  ROS_BACKU(IIG,JJG) = SF%VEG_LSET_ROS_BACK
  WIND_EXP(IIG,JJG)  = SF%VEG_LSET_WIND_EXP

!If any surfaces uses tan^2 function for slope, tan^2 will be used throughout simulation
  IF (SF%VEG_LSET_TAN2) LSET_TAN2=.TRUE.

!Use assumed ellipse shape of fireline as in Farsite
  IF_ELLIPSE_UNCOUPLED:IF (SF%VEG_LSET_ELLIPSE) THEN

    ROS_HEAD(IIG,JJG) = SF%VEG_LSET_ELLIPSE_HEAD
    SF%VEG_LSET_HT = MAX(0.001_EB,SF%VEG_LSET_HT)

    !If any surfaces set to ellipse, then elliptical model used for all surfaces
    IF (.NOT. LSET_ELLIPSE) LSET_ELLIPSE=.TRUE.

!Find wind at ~6.1 m height for Farsite
     KWIND = 0
     DO KDUM = KKG,KBAR
      IF(ZC(KDUM)-ZC(KKG) >= 6.1_EB) THEN
       KWIND = KDUM
!      WRITE(LU_OUTPUT,*)'vege: z(+3),z(kwind),kwind',z(kkg+3),z(kwind),kwind
       KKG=KBAR+1
      ENDIF
     ENDDO
     IF (ZC(KBAR) < 6.1_EB) KWIND=1

     U_LS(IIG,JJG) = U(IIG,JJG,KWIND)
     V_LS(IIG,JJG) = V(IIG,JJG,KWIND)
!WRITE(LU_OUTPUT,*)'veg initialize LS:u_ls(i,j)kwind',u_ls(iig,jjg)

    !Wind at midflame height (UMF).
    !From Andrews 2012, USDA FS Gen Tech Rep. RMRS-GTR-266 (with added SI conversion)
    UMF_TMP = 1.83_EB / LOG((20.0_EB + 1.18_EB * SF%VEG_LSET_HT) /(0.43_EB * SF%VEG_LSET_HT))

    !Factor 60 converts U from m/s to m/min which is used in elliptical model.
    UMF_X = UMF_TMP * U_LS(IIG,JJG) * 60.0_EB
    UMF_Y = UMF_TMP * V_LS(IIG,JJG) * 60.0_EB

    !Variables used in Phi_W formulas below (Rothermel model)
    B_ROTH = 0.15988_EB * (SF%VEG_LSET_SIGMA**0.54_EB)
    C_ROTH = 7.47_EB * EXP(-0.8711_EB * (SF%VEG_LSET_SIGMA**0.55_EB))
    E_ROTH = 0.715_EB * EXP(-0.01094_EB * SF%VEG_LSET_SIGMA)
    BETA_OP_ROTH = 0.20395_EB * (SF%VEG_LSET_SIGMA**(-0.8189_EB))! Optimum packing ratio

    ! Components of wind factor - affects spread rate
    PHI_W_X = C_ROTH * ((3.281_EB * ABS(UMF_X))**B_ROTH) * (SF%VEG_LSET_BETA / BETA_OP_ROTH)**(-E_ROTH)
    PHI_W_X = SIGN(PHI_W_X,UMF_X)

    PHI_W_Y = C_ROTH * ((3.281_EB * ABS(UMF_Y))**B_ROTH) * (SF%VEG_LSET_BETA / BETA_OP_ROTH)**(-E_ROTH)
    PHI_W_Y = SIGN(PHI_W_Y,UMF_Y)

    PHI_W(IIG,JJG) =  SQRT(PHI_W_X**2 + PHI_W_Y**2)

    !Limit effect to slope lte 80 degrees
    !Phi_s_x,y are slope factors
    DZT_DUM = MIN(5.67_EB,ABS(DZTDX(IIG,JJG))) ! 5.67 ~ tan 80 deg
    !DZT_DUM = SIGN(DZT_DUM,DZTDX)
    PHI_S_X(IIG,JJG) = 5.275_EB * ((SF%VEG_LSET_BETA)**(-0.3_EB)) * DZT_DUM**2
    PHI_S_X(IIG,JJG) = SIGN(PHI_S_X(IIG,JJG),DZTDX(IIG,JJG))

    DZT_DUM = MIN(1.73_EB,ABS(DZTDY(IIG,JJG))) ! 1.73 ~ tan 60 deg
    !DZT_DUM = SIGN(DZT_DUM,DZTDY)
    PHI_S_Y(IIG,JJG) = 5.275_EB * ((SF%VEG_LSET_BETA)**(-0.3_EB)) * DZT_DUM**2
    PHI_S_Y(IIG,JJG) = SIGN(PHI_S_Y(IIG,JJG),DZTDY(IIG,JJG))

    MAG_PHI_S = SQRT(PHI_S_X(IIG,JJG)**2 + PHI_S_Y(IIG,JJG)**2)

    PHI_S(IIG,JJG) = MAG_PHI_S  !5.275 * MAG_ZT(I,J)**2 * (SF%VEG_LSET_BETA)**-0.3

    ! Slope factor

    IF (MAG_PHI_S > 0.0_EB) THEN

        PHX = PHI_W_X + PHI_S_X(IIG,JJG)
        PHY = PHI_W_Y + PHI_S_Y(IIG,JJG)
        MAG_PHI = SQRT(PHX**2 + PHY**2)

        !Total phi (phi_w + phi_s) for use in spread rate section
        PHI_WS(IIG,JJG) = MAG_PHI

        !Theta_elps is angle of direction (0 to 2pi) of highest spread rate
        !0<=theta_elps<=2pi as measured clockwise from Y-axis
        THETA_ELPS(IIG,JJG) = ATAN2(PHY,PHX)

        !"Effective midflame windspeed" used in length-to-breadth ratio calculation (spread rate routine)
        ! is the wind + slope effect obtained by solving Phi_w eqs. above for UMF
        ! 8/8/13 - Changed phi_ws to Phi_s below to match Farsite, i.e., instead of adding phi_w and phi_s
        ! and then calculating effective wind speed, phi_s is converted to an effected wind speed and added
        ! to UMF calculated from the wind. Effective U has units of m/min in Wilson formula.
        ! 0.3048 ~= 1/3.281
        !if phi_s < 0 then a complex value (NaN) results. Using abs(phi_s) and sign function to correct.

        UMF_TMP = (((ABS(PHI_S_X(IIG,JJG)) * (SF%VEG_LSET_BETA / BETA_OP_ROTH)**E_ROTH)/C_ROTH)**(1/B_ROTH))*0.3048
        UMF_TMP = SIGN(UMF_TMP,PHI_S_X(IIG,JJG))
        UMF_X = UMF_X + UMF_TMP

        UMF_TMP = (((ABS(PHI_S_Y(IIG,JJG)) * (SF%VEG_LSET_BETA / BETA_OP_ROTH)**E_ROTH)/C_ROTH)**(1/B_ROTH))*0.3048
        UMF_TMP = SIGN(UMF_TMP,PHI_S_Y(IIG,JJG))
        UMF_Y = UMF_Y + UMF_TMP


        !write(*,*) 'phi_s_x',phi_s_x(iig,jjg)
        !write(*,*) 'phi_s_y',phi_s_y(iig,jjg)
        !write(*,*) 'umf_x',umf_x
        !write(*,*) 'umf_y',umf_y

    ELSE

        PHI_WS(IIG,JJG) = SQRT (PHI_W_X**2 + PHI_W_Y**2)
        !IF (PHY == 0._EB) PHY = 1.E-6_EB
        !0<= Theta_elps <=2pi as measured clockwise from Y-axis
        THETA_ELPS(IIG,JJG) = ATAN2(PHI_W_Y,PHI_W_X)

    ENDIF

    UMF(IIG,JJG) = SQRT(UMF_X**2 + UMF_Y**2)

     !The following two lines convert ATAN2 output to compass system (0 to 2 pi CW from +Y-axis)
     THETA_ELPS(IIG,JJG) = PIO2 - THETA_ELPS(IIG,JJG)
     IF (THETA_ELPS(IIG,JJG) < 0.0_EB) THETA_ELPS(IIG,JJG) = 2.0_EB*PI + THETA_ELPS(IIG,JJG)

  ENDIF IF_ELLIPSE_UNCOUPLED
!WRITE(LU_OUTPUT,*)'veg LS init:i,j,ros_head',sf%veg_lset_ellipse,iig,jjg,ros_head(iig,jjg)

ENDDO LSET_INIT_WALL_CELL_LOOP

UMAX_LS  = MAXVAL(ABS(U_LS))
VMAX_LS  = MAXVAL(ABS(V_LS))
UMAG     = SQRT(UMAX_LS**2 + VMAX_LS**2)

!WRITE(LU_OUTPUT,*)'before assign ROS'
WRITE(LU_OUTPUT,*)'ROS_HEAD max',MAXVAL(ROS_HEAD)
ROS_HEAD1 = MAXVAL(ROS_HEAD)
WRITE(LU_OUTPUT,*)'ROS_HEAD1',ROS_HEAD1

SR_MAX   = MAXVAL(ROS_HEAD)
DYN_SR_MAX = 0._EB

!WRITE(LU_OUTPUT,*)'SR_MAX',sr_max
SR_MAX   = MAX(SR_MAX,MAXVAL(ROS_FLANK))

!Flank rate not available when ellipse/farsite model is used
IF (LSET_ELLIPSE) THEN
    SR_MAX   = MAXVAL(ROS_HEAD) * (1._EB + MAXVAL(PHI_S) + MAXVAL(PHI_W))
    WRITE(LU_OUTPUT,*)'Phi_S max',MAXVAL(PHI_S)
    WRITE(LU_OUTPUT,*)'Phi_W max',MAXVAL(PHI_W)
    WRITE(LU_OUTPUT,*)'UMF max',MAXVAL(UMF)
    WRITE(LU_OUTPUT,*)'Mag_zt max',MAXVAL(MAG_ZT)
    WRITE(LU_OUTPUT,*)'SR_MAX',sr_max
ENDIF
IF (.NOT. LSET_ELLIPSE) SR_MAX   = 2._EB*SR_MAX !rough accounting for upslope spread aligned with wind

DT_LS = 0.25_EB*MIN(DX_LS,DY_LS)/SR_MAX

!DT_LS = 0.1603_EB !to make AU F19 ignition sequence work

WRITE(LU_OUTPUT,*)'vege: t_final,dt_ls',t_final,dt_ls
WRITE(LU_OUTPUT,*)'flux limiter= ',LIMITER_LS

END SUBROUTINE INITIALIZE_LEVEL_SET_FIRESPREAD

! ************************************************************************************************
! ************************************************************************************************
SUBROUTINE LEVEL_SET_FIRESPREAD(T_CFD,DT,NM)
! ************************************************************************************************
! ************************************************************************************************
!
INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: T_CFD,DT
INTEGER :: J_FLANK,I,IIG,IW,J,JJG,KKG
INTEGER :: IDUM,JDUM,KDUM,KWIND
LOGICAL :: IGNITION = .FALSE.
REAL(EB) :: HEAD_WIDTH_FCTR,IGNITION_WIDTH_Y,ROS_FLANK1,TIME_LS_LAST
!REAL(EB) :: PHI_CHECK
REAL(FB) :: TIME_LS_OUT

!---Vars for Farsite emulation (phi_w calculation)---
REAL(EB) :: UMF_TMP,PHX,PHY,MAG_PHI
REAL(EB) :: PHI_W_X,PHI_W_Y,UMF_X,UMF_Y

TYPE (WALL_TYPE),    POINTER :: WC =>NULL()
TYPE (SURFACE_TYPE), POINTER :: SF =>NULL()
CALL POINT_TO_MESH(NM)

!--- Initialize variables
HEAD_WIDTH_FCTR  = 1._EB
IGNITION_WIDTH_Y = 1
J_FLANK          = 1
ROS_FLANK1       = 0._EB

! LS_NEW ->
!IF (VEG_LEVEL_SET_COUPLED) THEN
! LS_NEW <-
 DT_LS   = DT
 TIME_LS = T_CFD
 T_FINAL = TIME_LS + DT_LS
! LS_NEW ->
!ENDIF
! LS_NEW <-
!WRITE(LU_OUTPUT,*)'vege: dt_ls,time_ls,t_final',dt_ls,time_ls,t_final
!
!-- Time step solution using second order Runge-Kutta -----------------------
!

! LS_NEW ->
!DO WHILE (TIME_LS < T_FINAL)
! LS_NEW <-

!
!-- Find flank-to-flank distance at base of fire assume symmetry about ymid and
!   define spread rate based on AU head fire width dependence
 IF (.NOT. LSET_ELLIPSE) THEN

!********************* Specific to AU grassland fuel experiments *************************
!  IF (SF%VEG_LSET_HEADWIDTH_DEPENDENCE) THEN
!    ROS_WALL_CELL_LOOP2: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
!     WC  => WALL(IW)
!     IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE ROS_WALL_CELL_LOOP2
!     SF  => SURFACE(WC%SURF_INDEX)
!     IF (.NOT. SF%VEG_LSET_SPREAD) CYCLE ROS_WALL_CELL_LOOP2
!     IF (.NOT. SF%VEG_LSET_HEADWIDTH_DEPENDENCE) CYCLE ROS_WALL_CELL_LOOP2
!     IIG = WC%ONE_D%IIG
!     JJG = WC%ONE_D%JJG
!!Case C064
!     IF(TIME_LS > 0._EB .AND. TIME_LS < 27._EB)  HEAD_WIDTH(I,J) = 2._EB*0.9_EB*TIME_LS !Ignition procdure 0.9 m/s rate
!     IF(TIME_LS >= 27._EB) HEAD_WIDTH(I,J) = HEAD_WIDTH(I,J) + 2._EB*ROS_FLANK(I,J)*(TIME_LS-TIME_LS_LAST)
!!Case F19
!!    IF(TIME_LS > 0._EB .AND. TIME_LS < 57._EB)  HEAD_WIDTH(I,J) = 2._EB*1.54_EB*TIME_LS !Ignition procdure 1.54 m/s rate
!!    IF(TIME_LS >= 57._EB .AND. TIME_LS < 100._EB) &
!!                              HEAD_WIDTH(I,J) = HEAD_WIDTH(I,J) + 2._EB*ROS_FLANK(I,J)*(TIME_LS-TIME_LS_LAST)
!!    IF(TIME_LS >= 100._EB) HEAD_WIDTH(I,J) = 100000._EB
!     HEAD_WIDTH_FCTR = EXP(-(0.859_EB + 2.036_EB*UMAG)/HEAD_WIDTH(I,J))
!     IF(ROS_HEAD(I,J) > 0.0_EB) ROS_HEAD(I,J)=ROS_HEAD1*HEAD_WIDTH_FCTR
!    ENDDO ROS_WALL_CELL_LOOP2
! ENDIF



!    IF(TIME_LS > 0._EB .AND. TIME_LS < 24._EB)  HEAD_WIDTH = 2._EB*1._EB*TIME_LS !Ignition procdure 1 m/s rate
!    IF(TIME_LS >= 24._EB) HEAD_WIDTH = HEAD_WIDTH + 2._EB*ROS_FLANK1*(TIME_LS - TIME_LS_LAST)
!    TIME_LS_LAST = TIME_LS
!    HEAD_WIDTH_FCTR = EXP(-(0.859_EB + 2.036_EB*UMAG)/HEAD_WIDTH)
!     DO J = 1,NY_LS
!      DO I = 1,NX_LS
!       IF(ROS_HEAD(I,J) > 0.0_EB) ROS_HEAD(I,J)=1.48*HEAD_WIDTH_FCTR
!      ENDDO
!     ENDDO

!    IF (HEAD_WIDTH_DEPENDENCE) THEN
!     IGNITION_WIDTH_Y = 3
!     J_FLANK = 0
!     DO JJ = NY_LS/2,NY_LS
!   !  IF(PHI_LS(26,JJ) <= 0.0_EB .AND. J_FLANK==0) J_FLANK = JJ
!      IF(PHI_LS(26,JJ) > 0.0_EB) J_FLANK = J_FLANK + 1
!     ENDDO
!   ! HEAD_WIDTH = 2._EB*(J_FLANK - NY_LS/2)*DY_LS
!     HEAD_WIDTH = 2.0_EB*J_FLANK*DY_LS
!     IF (HEAD_WIDTH < IGNITION_WIDTH_Y) HEAD_WIDTH = IGNITION_WIDTH_Y
!     HEAD_WIDTH_FCTR = EXP(-(0.859_EB + 2.036_EB*UMAG)/HEAD_WIDTH)
!     DO J = 1,NY_LS
!      DO I = 1,NX_LS
!       IF(ROS_HEAD(I,J) > 0.0_EB) ROS_HEAD(I,J)=ROS_HEAD1*HEAD_WIDTH_FCTR
!      ENDDO
!     ENDDO
!    ENDIF
 ENDIF
! ****************************************************************************************

 TIME_LS_LAST = TIME_LS

! LS_NEW ->
!!----------------- Output time steps with increasing step (as in FDS)-------------------------
!IF ( (TIME_LS<=10.0_EB) .OR. (SUMTIME > 100.0_EB) ) THEN
! SUMTIME = 0._EB
! WRITE(LU_OUTPUT,*)'vege:LS:-------------------------------------'
! WRITE(LU_OUTPUT,*)'vege:LS:time_ls',time_ls
! WRITE(LU_OUTPUT,*)'vege:ls:dt',dt_ls
!!WRITE(LU_OUTPUT,*)'vege:LS:HW,ros_h',head_width(nx_ls/2,ny_ls/2),ros_head(nx_ls/2,ny_ls/2)
!!WRITE(LU_OUTPUT,*)'vege:LS:ros_f',ros_flank(nx_ls/2,ny_ls/2)
! WRITE(LU_OUTPUT,*)'vege:LS:max_ros for time stepping',dyn_sr_max
!ENDIF
!!------------------------------------------------------------------------------------------------
! LS_NEW <-

WALL_CELL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
  WC  => WALL(IW)
  IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE WALL_CELL_LOOP
  SF  => SURFACE(WC%SURF_INDEX)
  IF (.NOT. SF%VEG_LSET_SPREAD) CYCLE WALL_CELL_LOOP
  IIG = WC%ONE_D%IIG
  JJG = WC%ONE_D%JJG


! Ignite landscape at user specified location(s) and time(s)
  IF (SF%VEG_LSET_IGNITE_T > 0.0_EB .AND. SF%VEG_LSET_IGNITE_T < DT_LS) PHI_LS(IIG,JJG) = PHI_MAX_LS
  IF (SF%VEG_LSET_IGNITE_T >= TIME_LS .AND. SF%VEG_LSET_IGNITE_T <= TIME_LS + DT_LS) PHI_LS(IIG,JJG) = PHI_MAX_LS
!    WRITE(LU_OUTPUT,*)'vege: sf%lset_ignite_t', sf%veg_lset_ignite_t

! Variable update when level set is coupled to CFD computation
  IF_CFD_COUPLED: IF (VEG_LEVEL_SET_COUPLED) THEN
   KKG = WC%ONE_D%KKG

!  IF (PHI_LS(IIG,JJG) == PHI_MAX_LS) WC%ONE_D%TMP_F = 373._EB
   WC%ONE_D%Q_CON_F = 0._EB
   IF (PHI_LS(IIG,JJG) <= 0.5_EB .AND. PHI_LS(IIG,JJG) >= -0.5_EB) THEN
!   WC%ONE_D%TMP_F = 373._EB
    WC%ONE_D%Q_CON_F = SF%VEG_LSET_QCON
!WRITE(LU_OUTPUT,*)'vege: LS QCONF set',SF%VEG_LSET_QCON
   ENDIF

   U_LS(IIG,JJG) = U(IIG,JJG,KKG)
   V_LS(IIG,JJG) = V(IIG,JJG,KKG)

!AU grassland ROS for infinite head and 6% moisutre
  UMAG     = SQRT(U_LS(IIG,JJG)**2 + V_LS(IIG,JJG)**2)
  ROS_HEAD(IIG,JJG)  = SF%VEG_LSET_ROS_HEAD*(0.165_EB + 0.534_EB*UMAG)*0.523_EB

!Use assumed ellipse shape of fireline as in Farsite
!WRITE(LU_OUTPUT,*)'sf%veg_lset_ellipse',SF%VEG_LSET_ELLIPSE
  IF_ELLIPSE_COUPLED: IF (SF%VEG_LSET_ELLIPSE) THEN
! IF_ELLIPSE_COUPLED: IF (LSET_ELLIPSE) THEN

    ROS_HEAD(IIG,JJG) = SF%VEG_LSET_ELLIPSE_HEAD
!WRITE(LU_OUTPUT,*)'time_ls',time_ls
!WRITE(LU_OUTPUT,*)'veg coupled LS1:i,j,ros_head',iig,jjg,ros_head(iig,jjg)

!Find wind at ~6.1 m height for Farsite
     KWIND = 0
     DO KDUM = KKG,KBAR
      IF(ZC(KDUM)-ZC(KKG) >= 6.1_EB) THEN
       KWIND = KDUM
       KKG=KBAR+1
      ENDIF
     ENDDO

     U_LS(IIG,JJG) = U(IIG,JJG,KWIND)
     V_LS(IIG,JJG) = V(IIG,JJG,KWIND)
!WRITE(LU_OUTPUT,*)'veg coupled LS:i,j,u_ls(i,j)kwind',iig,jjg,u_ls(iig,jjg)

     !Wind at midflame height (UMF).
     !From Andrews 2012, USDA FS Gen Tech Rep. RMRS-GTR-266 (with added SI conversion)
     UMF_TMP = 1.83_EB / LOG((20.0_EB + 1.18_EB * SF%VEG_LSET_HT) /(0.43_EB * SF%VEG_LSET_HT))

     !Factor 60 converts U from m/s to m/min which is used in elliptical model.
     UMF_X = UMF_TMP * U_LS(IIG,JJG) * 60.0_EB
     UMF_Y = UMF_TMP * V_LS(IIG,JJG) * 60.0_EB

     ! Components of wind factor - affects spread rate
     PHI_W_X = C_ROTH * ((3.281_EB * ABS(UMF_X))**B_ROTH) * (SF%VEG_LSET_BETA / BETA_OP_ROTH)**(-E_ROTH)
     PHI_W_X = SIGN(PHI_W_X,UMF_X)

     PHI_W_Y = C_ROTH * ((3.281_EB * ABS(UMF_Y))**B_ROTH) * (SF%VEG_LSET_BETA / BETA_OP_ROTH)**(-E_ROTH)
     PHI_W_Y = SIGN(PHI_W_Y,UMF_Y)

     PHI_W(IIG,JJG) =  SQRT(PHI_W_X**2 + PHI_W_Y**2)
!WRITE(LU_OUTPUT,*)'time_ls',time_ls
!WRITE(LU_OUTPUT,*)'veg coupled LS:i,j,phi_w(i,j)',iig,jjg,phi_w(iig,jjg)

     ! Slope factor
     IF (PHI_S(IIG,JJG) > 0.0_EB) THEN

         PHX = PHI_W_X + PHI_S_X(IIG,JJG)
         PHY = PHI_W_Y + PHI_S_Y(IIG,JJG)
         MAG_PHI = SQRT(PHX**2 + PHY**2)

         !Total phi (phi_w + phi_s) for use in spread rate section
         PHI_WS(IIG,JJG) = MAG_PHI

         !Theta_elps is angle of direction (0 to 2pi) of highest spread rate
         !0<=theta_elps<=2pi as measured clockwise from Y-axis
         THETA_ELPS(IIG,JJG) = ATAN2(PHY,PHX)

         !"Effective midflame windspeed" used in length-to-breadth ratio calculation (spread rate routine)
         ! is the wind + slope effect obtained by solving Phi_w eqs. above for UMF
         ! 8/8/13 - Changed phi_ws to Phi_s below to match Farsite, i.e., instead of adding phi_w and phi_s
         ! and then calculating effective wind speed, phi_s is converted to an effected wind speed and added
         ! to UMF calculated from the wind. Effective U has units of m/min in Wilson formula.
         ! 0.3048 ~= 1/3.281
         !if phi_s < 0 then a complex value (NaN) results. Using abs(phi_s) and sign function to correct.

         UMF_TMP = (((ABS(PHI_S_X(IIG,JJG)) * (SF%VEG_LSET_BETA / BETA_OP_ROTH)**E_ROTH)/C_ROTH)**(1/B_ROTH))*0.3048
         UMF_TMP = SIGN(UMF_TMP,PHI_S_X(IIG,JJG))
         UMF_X = UMF_X + UMF_TMP

         UMF_TMP = (((ABS(PHI_S_Y(IIG,JJG)) * (SF%VEG_LSET_BETA / BETA_OP_ROTH)**E_ROTH)/C_ROTH)**(1/B_ROTH))*0.3048
         UMF_TMP = SIGN(UMF_TMP,PHI_S_Y(IIG,JJG))
         UMF_Y = UMF_Y + UMF_TMP


         !write(*,*) 'phi_s_x',phi_s_x(iig,jjg)
         !write(*,*) 'phi_s_y',phi_s_y(iig,jjg)
         !write(*,*) 'umf_x',umf_x
         !write(*,*) 'umf_y',umf_y

     ELSE

         PHI_WS(IIG,JJG) = SQRT(PHI_W_X**2 + PHI_W_Y**2)
         !IF (PHY == 0._EB) PHY = 1.E-6_EB
         !0<= Theta_elps <=2pi as measured clockwise from Y-axis
         THETA_ELPS(IIG,JJG) = ATAN2(PHI_W_Y,PHI_W_X)

     ENDIF

     UMF(IIG,JJG) = SQRT(UMF_X**2 + UMF_Y**2)

     !The following two lines convert ATAN2 output to compass system (0 to 2 pi CW from +Y-axis)
      THETA_ELPS(IIG,JJG) = PIO2 - THETA_ELPS(IIG,JJG)
      IF (THETA_ELPS(IIG,JJG) < 0.0_EB) THETA_ELPS(IIG,JJG) = 2.0_EB*PI + THETA_ELPS(IIG,JJG)

   ENDIF IF_ELLIPSE_COUPLED

  ENDIF IF_CFD_COUPLED

!WRITE(LU_OUTPUT,*)'veg coupled LS2:i,j,ros_head',iig,jjg,ros_head(iig,jjg)
!WRITE(LU_OUTPUT,*)'----'

ENDDO WALL_CELL_LOOP

IF (ANY(PHI_LS==PHI_MAX_LS)) IGNITION = .TRUE.

PHI_TEMP = PHI_LS
!--- RK Stage 1
 RK2_PREDICTOR_LS = .TRUE.
 CALL LEVEL_SET_SPREAD_RATE
 CALL LEVEL_SET_ADVECT_FLUX
 PHI1_LS = PHI_LS - DT_LS*FLUX0_LS

!--- RK Stage2
 RK2_PREDICTOR_LS = .FALSE.
 MAG_SR_OUT       = 0.0_EB
 CALL LEVEL_SET_SPREAD_RATE
 CALL LEVEL_SET_ADVECT_FLUX
 PHI_LS = PHI_LS - 0.5_EB*DT_LS*(FLUX0_LS + FLUX1_LS)

! LS_NEW ->
!IF (VEG_LEVEL_SET_UNCOUPLED) THEN
!
!!-- Variable Time Step for simulations uncoupled from the CFD computation
!
! PHI_CHECK = MAXVAL(ABS(PHI_LS - PHI_TEMP)) !Find max change in phi
!
! !write(*,*)"Phi check ",phi_check
! !write(*,*)"dt_coef",dt_coef
!
! IF (IGNITION) THEN
!     ! If any phi values change by more than 0.5, or all change less
!     ! than 0.1 (both values are arbitrary), during one time step,
!     ! then adjust time step accordingly.
!
!     IF (PHI_CHECK > 0.5_EB) THEN
!         ! Cut time step in half and cycle the do-while loop
!         DT_COEF = 0.5_EB * DT_COEF
!         DT_LS = DT_COEF * MIN(DX_LS,DY_LS)/DYN_SR_MAX
!         DT_LS = MIN(DT_LS,100._EB)
!         PHI_LS = PHI_TEMP ! Exchange for previous phi and cycle
!         WRITE(*,*)"Halving time step."
!         CYCLE
!     ENDIF
!
!     ! Increase time step by 1/4 if changes are small
!     IF (PHI_CHECK < 0.1_EB) DT_COEF = DT_COEF * 1.25_EB
!
!     ! Dynamic Spread Rate Max
!     DYN_SR_MAX = MAX(DYN_SR_MAX,0.01_EB) ! dyn_sr_max must be g.t. zero
!     DT_LS = DT_COEF * MIN(DX_LS,DY_LS)/DYN_SR_MAX
!     DT_LS = MIN(DT_LS,100._EB)
!
! ENDIF
!ENDIF
! LS_NEW <-

 ! **********************************************************

 !Construct Time of Arrival (TOA) grid matching the horizontal domain grid

 DO IDUM=1,NX_LS
     DO JDUM=1,NY_LS
         !Change TOA to time value if Phi >=0. NOTE: Phi_ls goes from -1 to 1
         IF (PHI_LS(IDUM,JDUM)>=0._EB .AND. TOA(IDUM,JDUM)<=-1.0_EB) TOA(IDUM,JDUM)=REAL(TIME_LS,FB)
     ENDDO
 ENDDO

 TIME_LS = TIME_LS + DT_LS
 SUMTIME = SUMTIME + DT_LS
 SUM_T_SLCF = SUM_T_SLCF + DT_LS

 !--- Output slice file for smokeview
! LS_NEW ->
!IF (SUM_T_SLCF >= DT_OUTPUT) THEN
 IF (SUM_T_SLCF >= DT_OUTPUT .AND. .NOT. RK2_PREDICTOR_LS) THEN
! LS_NEW <-
  SUM_T_SLCF = 0._EB
  PHI_OUT = REAL(PHI_LS,FB)
  TIME_LS_OUT = REAL(TIME_LS,FB)
  WRITE(LU_SLCF_LS) TIME_LS_OUT
!negative for consistency with wall thickness output from wfds
  WRITE(LU_SLCF_LS) ((-PHI_OUT(I,J),I=1,NX_LS),J=1,NY_LS)
 ENDIF
!
! LS_NEW ->
!ENDDO !While loop
! LS_NEW <-

! LS_NEW ->
! ******  Write arrays to file **************
!IF (VEG_LEVEL_SET_UNCOUPLED) THEN
! CALL CPU_TIME(CPUTIME)
! LS_T_END = CPUTIME
! WRITE(LU_OUTPUT,*)'Uncoupled Level Set CPU Time: ',LS_T_END - LS_T_BEG
!!
! LU_TOA_LS = GET_FILE_NUMBER()
! OPEN(LU_TOA_LS,FILE='time_of_arrival.toa',STATUS='REPLACE')
! WRITE(LU_TOA_LS,'(I5)') NX_LS,NY_LS
! WRITE(LU_TOA_LS,'(F7.2)') XS,XF,YS,YF
!!Write across row (TOA(1,1), TOA(1,2), ...) to match Farsite output
! WRITE(LU_TOA_LS,'(F7.2)') ((TOA(IDUM,JDUM),JDUM=1,NY_LS),IDUM=1,NX_LS)
! CLOSE(LU_TOA_LS)
!
!CLOSE(LU_SLCF_LS)
!CLOSE(LU_TOA_LS)
!
!ENDIF
! LS_NEW <-

END SUBROUTINE LEVEL_SET_FIRESPREAD

SUBROUTINE LEVEL_SET_SPREAD_RATE
!
! Compute components of spread rate vector
!
INTEGER :: I,J,IM1,IP1,JM1,JP1
REAL(EB) :: COS_THETA_WIND,COS_THETA_SLOPE,COS_THETA_WIND_H,COS_THETA_WIND_B, &
            COS_THETA_SLOPE_H,COS_THETA_SLOPE_B,DPHIDX,DPHIDY,F_EAST,F_WEST,F_NORTH,F_SOUTH, &
            GRAD_SLOPE_DOT_NORMAL_FIRELINE,MAG_F,MAG_SR,MAG_U,WIND_DOT_NORMAL_FIRELINE,NEXP_WIND
REAL(EB) :: RAD_TO_DEGREE,DEGREES_SLOPE,SLOPE_FACTOR

!Variables for elliptical propagation model

REAL(EB) :: COS_THETA,SIN_THETA,XSF,YSF,UMF_DUM
REAL(EB) :: A_ELPS,A_ELPS2,B_ELPS2,B_ELPS,C_ELPS,DENOM,ROS_TMP,LB,LBD,HB
REAL(EB), DIMENSION(:) :: NORMAL_FIRELINE(2)

RAD_TO_DEGREE = 90._EB/ASIN(1._EB)

!NEXP_WIND = 2

IF (RK2_PREDICTOR_LS) PHI0_LS = PHI_LS
IF (.NOT. RK2_PREDICTOR_LS) PHI0_LS = PHI1_LS
SR_X_LS = 0.0_EB ; SR_Y_LS = 0.0_EB
DYN_SR_MAX = 0.0_EB

FLUX_ILOOP: DO I = 1,NX_LS

  IM1=I-1
  IP1=I+1
  IF (I==1) IM1 = I
  IF (I==NX_LS) IP1 = I

  DO J = 1,NY_LS

    JM1=J-1
    JP1=J+1
    IF (J==1) JM1 = J
    IF (J==NX_LS) JP1 = J

    F_EAST  = 0.5_EB*( PHI0_LS(I,J) + PHI0_LS(IP1,J) )
    F_WEST  = 0.5_EB*( PHI0_LS(I,J) + PHI0_LS(IM1,J) )
    F_NORTH = 0.5_EB*( PHI0_LS(I,J) + PHI0_LS(I,JP1) )
    F_SOUTH = 0.5_EB*( PHI0_LS(I,J) + PHI0_LS(I,JM1) )

   DPHIDX = (F_EAST-F_WEST) * IDX_LS
   DPHIDY = (F_NORTH-F_SOUTH) * IDY_LS

   !------------------------------------------------------------------------
   !    The two 'if blocks' below check for rare cases of symmetrically merging fire lines
   !     where the central difference may be zero but there are forward or
   !     backward differences.
   !------------------------------------------------------------------------
   !IF (ABS(DPHIDX) < EPS) THEN
   !
   !     IF ( PHI0_LS(IP1,J) > PHI0_LS(I,J) ) THEN
   !         DPHIDX = F_EAST * IDX_LS
   !     ENDIF
   !
   !     IF ( PHI0_LS(I,JP1) > PHI0_LS(I,J) ) THEN
   !         DPHIDY = F_NORTH * IDY_LS
   !     ENDIF
   !
   ! ENDIF

   !------------------------------------------------------------------------
   !    The two if statements below check for rare cases of symmetrically merging fire lines
   !     where the central difference may be zero but there are forward or
   !     backward differences
   !------------------------------------------------------------------------
   !IF ((DPHIDX < EPS) .AND. (F_EAST > 0._EB)) THEN
   !     IF ((F_EAST == F_WEST) .AND. (PHI0_LS(I,J) < F_EAST)) THEN
   !         DPHIDX =  (PHI0_LS(IP1,J) - PHI0_LS(I,J))* IDX_LS
   !     ENDIF
   !ENDIF
   !
   !IF ((DPHIDY < EPS) .AND. (F_NORTH > 0._EB)) THEN
   !     IF ((F_NORTH == F_SOUTH) .AND. (PHI0_LS(I,J) < F_NORTH)) THEN
   !         DPHIDY = (PHI0_LS(I,JP1) - PHI0_LS(I,J) )* IDY_LS
   !     ENDIF
   !ENDIF

   MAG_F = SQRT(DPHIDX**2 + DPHIDY**2)
   IF (MAG_F > 0._EB) THEN   !components of unit vector normal to PHI contours
        NORMAL_FIRELINE(1) = -DPHIDX/MAG_F
        NORMAL_FIRELINE(2) = -DPHIDY/MAG_F
        ! Lagrangian normal approximation from Rehm and Mcdermott 2009
        ! For elliptical front calculation
        XSF = (DPHIDY / MAG_F ) !* DX_LS
        YSF = (-DPHIDX / MAG_F) !* DY_LS

        GRAD_SLOPE_DOT_NORMAL_FIRELINE = DZTDX(I,J)*(DPHIDY/MAG_F) + DZTDY(I,J)*(-DPHIDY/MAG_F)

   ELSE
        NORMAL_FIRELINE = 0._EB
        GRAD_SLOPE_DOT_NORMAL_FIRELINE = 0._EB
        XSF=0._EB
        YSF=0._EB
   ENDIF

    ! ???
   COS_THETA_SLOPE = 0.0_EB ; COS_THETA_SLOPE_H = 0.0_EB ; COS_THETA_SLOPE_B = 0.0_EB

   IF (MAG_ZT(I,J) > 0.0_EB) THEN
       COS_THETA_SLOPE = GRAD_SLOPE_DOT_NORMAL_FIRELINE/MAG_ZT(I,J)
       !XSF = XSF * COS_THETA_SLOPE
       !YSF = YSF * COS_THETA_SLOPE
   ENDIF

   DEGREES_SLOPE = ATAN(MAG_ZT(I,J))*RAD_TO_DEGREE

   IF (LSET_ELLIPSE) THEN

       ! Effective wind direction (theta) is clockwise from y-axis (Richards 1990)
       COS_THETA = COS(THETA_ELPS(I,J)) !V_LS(I,J) / MAG_U
       SIN_THETA = SIN(THETA_ELPS(I,J)) !U_LS(I,J) / MAG_U

       !ROS_TMP = ROS_HEAD(I,J) * (1.0_EB + PHI_S(I,J) + PHI_W(I,J))
       ROS_TMP = ROS_HEAD(I,J) * (1.0_EB + PHI_WS(I,J))
!IF (ROS_HEAD(I,J) > 0._EB)WRITE(LU_OUTPUT,*)'vege levelset ros: ros_head',ros_head(i,j)

       !Mag of wind speed at midflame ht must be in units of m/s here
       UMF_DUM = UMF(I,J)/60.0_EB

       !Length to breadth ratio of ellipse based on effective UMF
       LB = 0.936_EB * EXP(0.2566_EB * UMF_DUM) + 0.461_EB * EXP(-0.1548_EB * UMF_DUM) - 0.397_EB

       !Constraint LB max = 8 from Finney 2004
       LB = MAX(1.0_EB,MIN(LB,8.0_EB))

       LBD = SQRT(LB**2 - 1.0_EB)

       !Head to back ratio based on LB
       HB = (LB + LBD) / (LB - LBD)

       !*** A_ELPS and B_ELPS are *opposite* in notation from Farsite and Richards
       A_ELPS =  0.5_EB * (ROS_TMP + ROS_TMP/HB)
       A_ELPS2 = A_ELPS**2
       B_ELPS =  A_ELPS / LB !0.5_EB * (ROS_TMP + ROS_TMP/HB) / LB
       B_ELPS2=  B_ELPS**2
       C_ELPS =  A_ELPS - (ROS_TMP/HB)

       ! Denominator used in spread rate equation from Richards 1990 (also in Farsite)
       DENOM = B_ELPS2 * (YSF * COS_THETA + XSF * SIN_THETA)**2 + &
                A_ELPS2 * (XSF * COS_THETA - YSF * SIN_THETA)**2

       ! Finney's formulation
       !DENOM = B_ELPS2 * (XS * SIN_THETA - YS * COS_THETA)**2 - &
       !A_ELPS2 * (XS * COS_THETA + YS * SIN_THETA)**2

       IF (DENOM > 0._EB) THEN

        DENOM = 1._EB / SQRT(DENOM)
       ELSE
        DENOM = 0._EB
       ENDIF

       !IF (PHI_W(I,J) > 0._EB .OR. MAG_ZT(I,J) > 0._EB) THEN

        SR_X_LS(I,J) = DENOM * (B_ELPS2 * COS_THETA * (XSF * SIN_THETA + YSF * COS_THETA) -&
                        A_ELPS2 * SIN_THETA * (XSF * COS_THETA - YSF * SIN_THETA)) + C_ELPS * SIN_THETA


        SR_Y_LS(I,J) = DENOM * (-B_ELPS2 * SIN_THETA * (XSF * SIN_THETA + YSF * COS_THETA) -&
                        A_ELPS2 * COS_THETA * (XSF * COS_THETA - YSF * SIN_THETA)) + C_ELPS * COS_THETA

       !ELSE

            !For no-wind, no-slope case
        !    SR_X_LS(I,J) = ROS_HEAD(I,J) * NORMAL_FIRELINE(1)
        !    SR_Y_LS(I,J) = ROS_HEAD(I,J) * NORMAL_FIRELINE(2)

       !ENDIF


       ! Project spread rates from slope to horizontal plane

       IF (ABS(DZTDX(I,J)) > 0._EB) SR_X_LS(I,J) = SR_X_LS(I,J) * ABS(COS(ATAN(DZTDX(I,J))))
       IF (ABS(DZTDY(I,J)) > 0._EB) SR_Y_LS(I,J) = SR_Y_LS(I,J) * ABS(COS(ATAN(DZTDY(I,J))))

       MAG_SR = SQRT(SR_X_LS(I,J)**2 + SR_Y_LS(I,J)**2)
!WRITE(LU_OUTPUT,*)'vege levelset ros: i,j,mag_sr',i,j,mag_sr

   ELSE !McArthur Spread Model

           WIND_DOT_NORMAL_FIRELINE = U_LS(I,J)*NORMAL_FIRELINE(1) + V_LS(I,J)*NORMAL_FIRELINE(2)
           MAG_U  = SQRT(U_LS(I,J)**2 + V_LS(I,J)**2)

           COS_THETA_WIND = 0.0_EB ; COS_THETA_WIND_H = 0.0_EB ; COS_THETA_WIND_B = 0.0_EB
           IF(MAG_U > 0.0_EB) COS_THETA_WIND = WIND_DOT_NORMAL_FIRELINE/MAG_U

           GRAD_SLOPE_DOT_NORMAL_FIRELINE = DZTDX(I,J)*NORMAL_FIRELINE(1) + DZTDY(I,J)*NORMAL_FIRELINE(2)
           COS_THETA_SLOPE = 0.0_EB ; COS_THETA_SLOPE_H = 0.0_EB ; COS_THETA_SLOPE_B = 0.0_EB

           IF (MAG_ZT(I,J) > 0.0_EB) THEN
               COS_THETA_SLOPE = GRAD_SLOPE_DOT_NORMAL_FIRELINE/MAG_ZT(I,J)
           ENDIF

           DEGREES_SLOPE = ATAN(MAG_ZT(I,J))*RAD_TO_DEGREE

           SLOPE_FACTOR  = MAG_ZT(I,J)**2
           IF (SLOPE_FACTOR > 3._EB) SLOPE_FACTOR = 3._EB


           ROS_HEADS = 0.33_EB*ROS_HEAD(I,J)
            IF(DEGREES_SLOPE >= 5._EB .AND. DEGREES_SLOPE < 10._EB)  ROS_HEADS = 0.33_EB*ROS_HEAD(I,J)
             IF(DEGREES_SLOPE >= 10._EB .AND. DEGREES_SLOPE < 20._EB) ROS_HEADS =         ROS_HEAD(I,J)
             IF(DEGREES_SLOPE >= 20._EB)                              ROS_HEADS =  3._EB*ROS_HEAD(I,J)

              MAG_SR    = 0.0_EB
              ROS_HEADS = 0.0_EB
              ROS_BACKS = 0.0_EB

             NEXP_WIND = WIND_EXP(I,J)

        ! Spread with the wind and upslope
           IF(COS_THETA_WIND >= 0._EB .AND. COS_THETA_SLOPE >= 0._EB) THEN
            IF (.NOT. LSET_TAN2) THEN
                IF(DEGREES_SLOPE >= 5._EB .AND. DEGREES_SLOPE < 10._EB)  ROS_HEADS = 0.33_EB*ROS_HEAD(I,J)
                IF(DEGREES_SLOPE >= 10._EB .AND. DEGREES_SLOPE < 20._EB) ROS_HEADS =         ROS_HEAD(I,J)
                IF(DEGREES_SLOPE >= 20._EB)                              ROS_HEADS =  3._EB*ROS_HEAD(I,J)
            ELSEIF (DEGREES_SLOPE > 0._EB) THEN
                    ROS_HEADS = ROS_HEAD(I,J) * SLOPE_FACTOR !Dependence on TAN(slope)^2
            ENDIF
            MAG_SR = ROS_FLANK(I,J)*(1._EB + COS_THETA_WIND**NEXP_WIND*COS_THETA_SLOPE) + &
                     (ROS_HEAD(I,J) - ROS_FLANK(I,J))*COS_THETA_WIND**NEXP_WIND + &
                     (ROS_HEADS     - ROS_FLANK(I,J))*COS_THETA_SLOPE  !magnitude of spread rate
           ENDIF
        !  IF(ABS(COS_THETA_WIND) < 0.5_EB .AND. MAG_F > 0._EB) MAG_SR = 0.0_EB
        !  IF(ABS(COS_THETA_WIND) < 0.5_EB .AND. MAG_F > 0._EB) FLANKFIRE_LIFETIME(I,J) = FLANKFIRE_LIFETIME(I,J) + DT_LS
        !  IF(FLANKFIRE_LIFETIME(I,J) > TIME_FLANKFIRE_QUENCH) MAG_SR = 0.0_EB

        ! Spread with the wind and downslope
           IF(COS_THETA_WIND >= 0._EB .AND. COS_THETA_SLOPE < 0._EB) THEN
            IF(DEGREES_SLOPE >= 5._EB .AND. DEGREES_SLOPE < 10._EB)  ROS_HEADS =  0.33_EB*ROS_HEAD(I,J)
            IF(DEGREES_SLOPE >= 10._EB .AND. DEGREES_SLOPE < 20._EB) ROS_HEADS =  0.50_EB*ROS_HEAD(I,J)
            IF(DEGREES_SLOPE >= 20._EB)                              ROS_HEADS =  0.75_EB*ROS_HEAD(I,J)
            MAG_SR = ROS_FLANK(I,J)*(1._EB + COS_THETA_WIND*COS_THETA_SLOPE) + &
                     (ROS_HEAD(I,J) - ROS_FLANK(I,J))*COS_THETA_WIND**NEXP_WIND + &
                     (ROS_HEADS     - ROS_FLANK(I,J))*COS_THETA_SLOPE  !magnitude of spread rate
        !   if(cos_theta_wind == 0._EB) FLANKFIRE_LIFETIME(I,J) = FLANKFIRE_LIFETIME(I,J) + DT_LS
        !   if(flankfire_lifetime(i,j) > time_flankfire_quench) mag_sr = 0.0_EB
           ENDIF

        ! Spread against the wind and upslope
           IF(COS_THETA_WIND <  0._EB .AND. COS_THETA_SLOPE >= 0._EB) THEN
            IF (.NOT. LSET_TAN2) THEN
                IF(DEGREES_SLOPE >= 5._EB .AND. DEGREES_SLOPE < 10._EB)  ROS_BACKS = -0.33_EB*ROS_BACKU(I,J)
                IF(DEGREES_SLOPE >= 10._EB .AND. DEGREES_SLOPE < 20._EB) ROS_BACKS =         -ROS_BACKU(I,J)
                IF(DEGREES_SLOPE >= 20._EB)                              ROS_BACKS = -3.0_EB*ROS_BACKU(I,J)
            ELSEIF (DEGREES_SLOPE > 0._EB) THEN
                ROS_HEADS = ROS_HEAD(I,J) * SLOPE_FACTOR !Dependence on TAN(slope)^2
            ENDIF
            MAG_SR = ROS_FLANK(I,J)*(1._EB - ABS(COS_THETA_WIND)**NEXP_WIND*COS_THETA_SLOPE) + &
                     (ROS_FLANK(I,J) - ROS_BACKU(I,J))*(-ABS(COS_THETA_WIND)**NEXP_WIND) + &
                     (ROS_FLANK(I,J) - ROS_BACKS)*COS_THETA_SLOPE  !magnitude of spread rate
           ENDIF

        ! Spread against the wind and downslope
           IF(COS_THETA_WIND <  0._EB .AND. COS_THETA_SLOPE < 0._EB) THEN
            IF(DEGREES_SLOPE >= 5._EB .AND. DEGREES_SLOPE < 10._EB)  ROS_BACKS = 0.33_EB*ROS_BACKU(I,J)
            IF(DEGREES_SLOPE >= 10._EB .AND. DEGREES_SLOPE < 20._EB) ROS_BACKS = 0.50_EB*ROS_BACKU(I,J)
            IF(DEGREES_SLOPE >= 20._EB)                              ROS_BACKS = 0.75_EB*ROS_BACKU(I,J)
            MAG_SR = ROS_FLANK(I,J)*(1._EB - ABS(COS_THETA_WIND)**NEXP_WIND*COS_THETA_SLOPE) + &
                     (ROS_FLANK(I,J) - ROS_BACKU(I,J))*(-ABS(COS_THETA_WIND)**NEXP_WIND) + &
                     (ROS_FLANK(I,J) - ROS_BACKS)*COS_THETA_SLOPE  !magnitude of spread rate
           ENDIF


        !  MAG_SR = ROS_FLANK(I,J) + ROS_HEAD(I,J)*COS_THETA_WIND**1.5 !magnitude of spread rate
        !  MAG_SR = ROS_FLANK(I,J) + ROS_HEAD(I,J)*MAG_U*COS_THETA_WIND**1.5 !magnitude of spread rate
           SR_X_LS(I,J) = MAG_SR*NORMAL_FIRELINE(1) !spread rate components
           SR_Y_LS(I,J) = MAG_SR*NORMAL_FIRELINE(2)
        !  MAG_SR_OUT(I,J) = MAG_SR

   ENDIF !Ellipse or McArthur Spread

   DYN_SR_MAX = MAX(DYN_SR_MAX,MAG_SR)

  ENDDO

ENDDO FLUX_ILOOP

END SUBROUTINE LEVEL_SET_SPREAD_RATE

!--------------------------------------------------------------------
!
SUBROUTINE LEVEL_SET_ADVECT_FLUX
!
! Use the spread rate [SR_X_LS,SR_Y_LS] to compute the limited scalar gradient
! and take dot product with spread rate vector to get advective flux

INTEGER :: I,IM1,IM2,IP1,IP2,J,JM1,JM2,JP1,JP2
REAL(EB), DIMENSION(:) :: Z(4)
REAL(EB), DIMENSION(:,:) :: FLUX_LS(NX_LS,NY_LS)
REAL(EB) :: DPHIDX,DPHIDY,F_EAST,F_WEST,F_NORTH,F_SOUTH
REAL(EB) :: PHIMAG

IF (RK2_PREDICTOR_LS) PHI0_LS = PHI_LS
IF (.NOT. RK2_PREDICTOR_LS) PHI0_LS = PHI1_LS

ILOOP: DO I = 1,NX_LS

 IM1=I-1; IF (IM1<1) IM1=IM1+NX_LS
 IM2=I-2; IF (IM2<1) IM2=IM2+NX_LS

 IP1=I+1; IF (IP1>NX_LS) IP1=IP1-NX_LS
 IP2=I+2; IF (IP2>NX_LS) IP2=IP2-NX_LS

 DO J = 1,NY_LS

   JM1=J-1; IF (JM1<1) JM1=JM1+NY_LS
   JM2=J-2; IF (JM2<1) JM2=JM2+NY_LS

   JP1=J+1; IF (JP1>NY_LS) JP1=JP1-NY_LS
   JP2=J+2; IF (JP2>NY_LS) JP2=JP2-NY_LS

!-- east face
   Z(1) = PHI0_LS(IM1,J)
   Z(2) = PHI0_LS(I,J)
   Z(3) = PHI0_LS(IP1,J)
   Z(4) = PHI0_LS(IP2,J)
   F_EAST = SCALAR_FACE_VALUE_LS(SR_X_LS(I,J),Z,LIMITER_LS)

!-- west face
   Z(1) = PHI0_LS(IM2,J)
   Z(2) = PHI0_LS(IM1,J)
   Z(3) = PHI0_LS(I,J)
   Z(4) = PHI0_LS(IP1,J)
   F_WEST = SCALAR_FACE_VALUE_LS(SR_X_LS(I,J),Z,LIMITER_LS)

 IF (J<2 .OR. J>(NY_LS-2)) THEN

      IF (J==1) THEN
        !    north face
            Z(1) = PHI_MAX_LS
            Z(2) = PHI0_LS(I,J)
            Z(3) = PHI0_LS(I,JP1)
            Z(4) = PHI0_LS(I,JP2)
            F_NORTH = SCALAR_FACE_VALUE_LS(SR_Y_LS(I,J),Z,LIMITER_LS)

        !    south face
            Z(1) = PHI_MAX_LS
            Z(2) = PHI_MAX_LS
            Z(3) = PHI0_LS(I,J)
            Z(4) = PHI0_LS(I,JP1)
            F_SOUTH = SCALAR_FACE_VALUE_LS(SR_Y_LS(I,J),Z,LIMITER_LS)

       ELSEIF (j==2) THEN
        !    north face
            Z(1) = PHI0_LS(I,JM1)
            Z(2) = PHI0_LS(I,J)
            Z(3) = PHI0_LS(I,JP1)
            Z(4) = PHI0_LS(I,JP2)
            F_NORTH = SCALAR_FACE_VALUE_LS(SR_Y_LS(I,J),Z,LIMITER_LS)

        !    south face
            Z(1) = PHI_MAX_LS
            Z(2) = PHI0_LS(I,JM1)
            Z(3) = PHI0_LS(I,J)
            Z(4) = PHI0_LS(I,JP1)
            F_SOUTH = SCALAR_FACE_VALUE_LS(SR_Y_LS(I,J),Z,LIMITER_LS)


        ELSEIF (J == NY_LS-1) THEN
    !    north face
            Z(1) = PHI0_LS(I,JM1)
            Z(2) = PHI0_LS(I,J)
            Z(3) = PHI0_LS(I,JP1)
            Z(4) = PHI_MIN_LS
            F_NORTH = SCALAR_FACE_VALUE_LS(SR_Y_LS(I,J),Z,LIMITER_LS)

        !    south face
            Z(1) = PHI0_LS(I,JM2)
            Z(2) = PHI0_LS(I,JM1)
            Z(3) = PHI0_LS(I,J)
            Z(4) = PHI0_LS(I,JP1)
            F_SOUTH = SCALAR_FACE_VALUE_LS(SR_Y_LS(I,J),Z,LIMITER_LS)

           ELSE ! must be J == NY_LS
        !    north face
            Z(1) = PHI0_LS(I,JM1)
            Z(2) = PHI0_LS(I,J)
            Z(3) = PHI_MIN_LS
            Z(4) = PHI_MIN_LS
            F_NORTH = SCALAR_FACE_VALUE_LS(SR_Y_LS(I,J),Z,LIMITER_LS)

        !    south face
            Z(1) = PHI0_LS(I,JM2)
            Z(2) = PHI0_LS(I,JM1)
            Z(3) = PHI0_LS(I,J)
            Z(4) = PHI_MIN_LS
            F_SOUTH = SCALAR_FACE_VALUE_LS(SR_Y_LS(I,J),Z,LIMITER_LS)

           ENDIF !IF (J==1)

       ELSE

    !    north face
       Z(1) = PHI0_LS(I,JM1)
       Z(2) = PHI0_LS(I,J)
       Z(3) = PHI0_LS(I,JP1)
       Z(4) = PHI0_LS(I,JP2)
       F_NORTH = SCALAR_FACE_VALUE_LS(SR_Y_LS(I,J),Z,LIMITER_LS)

    !    south face
       Z(1) = PHI0_LS(I,JM2)
       Z(2) = PHI0_LS(I,JM1)
       Z(3) = PHI0_LS(I,J)
       Z(4) = PHI0_LS(I,JP1)
       F_SOUTH = SCALAR_FACE_VALUE_LS(SR_Y_LS(I,J),Z,LIMITER_LS)

   ENDIF !IF (J<2 .OR. J>(NY_LS-2)

   DPHIDX = (F_EAST-F_WEST)* IDX_LS

   DPHIDY = (F_NORTH-F_SOUTH)* IDY_LS

   FLUX_LS(I,J) = SR_X_LS(I,J)*DPHIDX + SR_Y_LS(I,J)*DPHIDY

   PHIMAG          = SQRT(DPHIDX**2 + DPHIDY**2)
   MAG_SR_OUT(I,J) = 0.0_EB
   IF (PHIMAG > 0.0_EB) MAG_SR_OUT(I,J) = REAL(FLUX_LS(I,J)/PHIMAG,FB)

!  fx = (f_east-f_west)/dx
!  fy = (f_north-f_south)/dy
!       phi(i,j) = phi0(i,j) - dt*[Fx(i,j) Fy(i,j)]*[fx fy]
 ENDDO

 FLUX_LS(:,1) = FLUX_LS(:,2)

ENDDO ILOOP

IF (RK2_PREDICTOR_LS) FLUX0_LS = FLUX_LS
IF (.NOT. RK2_PREDICTOR_LS) FLUX1_LS = FLUX_LS

END SUBROUTINE LEVEL_SET_ADVECT_FLUX
!
! ----------------------------------------------------
REAL(EB) FUNCTION SCALAR_FACE_VALUE_LS(SR_XY,Z,LIMITER)
!
! From Randy 7-11-08
! This function computes the scalar value on a face.
! The scalar is denoted Z, and the velocity is denoted U.
! The gradient (computed elsewhere) is a central difference across
! the face subject to a flux limiter.  The flux limiter choices are:
!
! limiter = 1 implements the MINMOD limiter
! limiter = 2 implements the SUPERBEE limiter of Roe
! limiter = 3 implements first-order upwinding (monotone)
!
!
!                    location of face
!
!                            f
!    |     o     |     o     |     o     |     o     |
!                     SRXY        SRXY
!                 (if f_east)  (if f_west)
!         Z(1)        Z(2)        Z(3)        Z(4)
!
INTEGER :: LIMITER
REAL(EB) :: SR_XY
REAL(EB), INTENT(IN), DIMENSION(4) :: Z
REAL(EB) :: B,DZLOC,DZUP,R,ZUP,ZDWN

IF (SR_XY > 0._EB) THEN
!     the flow is left to right
 DZLOC = Z(3)-Z(2)
 DZUP  = Z(2)-Z(1)

 IF (ABS(DZLOC) > 0._EB) THEN
  R = DZUP/DZLOC
 ELSE
  R = 0._EB
 ENDIF
 ZUP  = Z(2)
 ZDWN = Z(3)
ELSE
!     the flow is right to left
 DZLOC = Z(3)-Z(2)
 DZUP  = Z(4)-Z(3)

 IF (ABS(DZLOC) > 0._EB) THEN
  R = DZUP/DZLOC
 ELSE
  R = 0._EB
 ENDIF
  ZUP  = Z(3)
  ZDWN = Z(2)
ENDIF

! flux limiter
IF (LIMITER==1) THEN
!     MINMOD
    B = MAX(0._EB,MIN(1._EB,R))
ELSEIF (limiter==2) THEN
!     SUPERBEE
    B = MAX(0._EB,MIN(2._EB*R,1._EB),MIN(R,2._EB))
ELSEIF (limiter==3) THEN
!     first-order upwinding
    B = 0._EB
ENDIF

SCALAR_FACE_VALUE_LS = ZUP + 0.5_EB * B * ( ZDWN - ZUP )

END FUNCTION SCALAR_FACE_VALUE_LS


SUBROUTINE LEVEL_SET_FIRE_SPREAD_NEW(T,DT,NM)

REAL(EB), INTENT(IN) :: T,DT
INTEGER, INTENT(IN) :: NM
INTEGER :: IW,IIG,JJG
TYPE(WALL_TYPE), POINTER :: WC=>NULL()
TYPE(SURFACE_TYPE), POINTER :: SF =>NULL()
REAL(EB) :: DUMMY

DUMMY = DT ! tmp to avoid unused warning
CALL POINT_TO_MESH(NM)

WALL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
   WC=>WALL(IW)
   IF (WC%BOUNDARY_TYPE/=SOLID_BOUNDARY) CYCLE WALL_LOOP
   SF=>SURFACE(WC%SURF_INDEX)
   IF (.NOT.SF%VEG_LSET_SPREAD) CYCLE WALL_LOOP

   IIG = WC%ONE_D%IIG
   JJG = WC%ONE_D%JJG

   IF (T > WC%ONE_D%T_IGN) WC%ONE_D%PHI_NLS = PHI_MAX_LS

ENDDO WALL_LOOP

END SUBROUTINE LEVEL_SET_FIRE_SPREAD_NEW

END MODULE VEGE
