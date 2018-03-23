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
PUBLIC INITIALIZE_LEVEL_SET_FIREFRONT,LEVEL_SET_FIREFRONT_PROPAGATION,END_LEVEL_SET, &
       BNDRY_VEG_MASS_ENERGY_TRANSFER,LEVEL_SET_BC,LEVEL_SET_DT,READ_BRNR
TYPE (LAGRANGIAN_PARTICLE_TYPE), POINTER :: LP=>NULL()
TYPE (LAGRANGIAN_PARTICLE_CLASS_TYPE), POINTER :: LPC=>NULL()
!TYPE (WALL_TYPE), POINTER :: WC
!TYPE (SURFACE_TYPE), POINTER :: SF 
INTEGER :: IZERO

!For Level Set
INTEGER  :: LIMITER_LS,LU_CRWN_PROB_LS,LU_FLI_LS,LU_ROSX_LS,LU_ROSY_LS,LU_TOA_LS  !&
!           ,LU_SLCF_LS,LU_SLCF_FLI_LS,LU_SLCF_PROBC_LS,LU_SLCF_ROS_LS,LU_SLCF_TOA_LS
REAL(EB) :: DX_LS,DY_LS,TIME_FLANKFIRE_QUENCH
!REAL(EB) :: DT_COEF_LS,DYN_SR_MAX_LS,DT_LS,SUM_T_SLCF_LS,SUMTIME_LS,TIME_LS
REAL(EB) :: IDX_LS,IDY_LS,T_FINAL,ROS_HEAD1,UMAG,UMF_TMP 
REAL(EB) :: CPUTIME,LS_T_BEG,LS_T_END,ROS_BACKS,ROS_HEADS
REAL(EB) :: B_ROTH,BETA_OP_ROTH,C_ROTH,E_ROTH

CONTAINS
 

SUBROUTINE BNDRY_VEG_MASS_ENERGY_TRANSFER(T,NM)
!
! Issues:
! 1. Are SF%VEG_FUEL_FLUX_L and SF%VEG_MOIST_FLUX_L needed in linear degradation model?
USE PHYSICAL_FUNCTIONS, ONLY : DRAG,GET_MASS_FRACTION,GET_SPECIFIC_HEAT,GET_VISCOSITY,GET_CONDUCTIVITY
REAL(EB) :: ZZ_GET(0:N_TRACKED_SPECIES)
REAL(EB) :: DT_BC,RDT_BC
REAL(EB), INTENT(IN) ::T
INTEGER, INTENT(IN) :: NM
INTEGER  ::  IW
INTEGER  ::  I,IIG,JJG,KKG,KGRID,KLOC
REAL(EB) :: CP_GAS,CP_MOIST_AND_VEG,DZVEG_L,ETAVEG_H,H_CONV_L, &
            KAPPA_VEG,K_GAS,MU_GAS,QRADM_INC,QRADP_INC,RHO_GAS, &
            TMP_BOIL,TMP_CHAR_MAX,TMP_FILM,TMP_G,DTMP_L,RE_VEG_PART,U2,V2,RE_D,Y_O2,ZVEG
!REAL(EB) :: H_CONV_FDS_WALL,DTMP_FDS_WALL,QCONF_FDS_WALL,LAMBDA_AIR,TMPG_A
INTEGER  IIVEG_L,IVEG_L,J,LBURN,NVEG_L,I_FUEL
!REAL(EB), ALLOCATABLE, DIMENSION(:) :: VEG_DIV_QRNET_EMISS,VEG_DIV_QRNET_INC,
!         VEG_QRNET_EMISS,VEG_QRNET_INC,VEG_QRM_EMISS,VEG_QRP_EMISS, VEG_QRM_INC,VEG_QRP_INC
REAL(EB) :: VEG_DIV_QRNET_EMISS(50),VEG_DIV_QRNET_INC(50),VEG_QRNET_EMISS(0:50),VEG_QRNET_INC(0:50), &
            VEG_QRM_EMISS(0:50),VEG_QRP_EMISS(0:50), VEG_QRM_INC(0:50),VEG_QRP_INC(0:50)
REAL(EB) :: H_H2O_VEG,A_H2O_VEG,E_H2O_VEG,H_PYR_VEG,A_PYR_VEG,E_PYR_VEG,RH_PYR_VEG,                  &
            H_CHAR_VEG,A_CHAR_VEG,E_CHAR_VEG,BETA_CHAR_VEG,NU_CHAR_VEG,NU_ASH_VEG,NU_O2_CHAR_VEG
REAL(EB) :: CP_ASH,CP_CHAR,CP_H2O,CP_VEG,CP_TOTAL,DTMP_VEG,Q_VEG_CHAR,TMP_VEG,TMP_VEG_NEW, &
            CHAR_ENTHALPY_FRACTION_VEG
REAL(EB) :: CHAR_FCTR,CHAR_FCTR2,MPA_MOIST,MPA_MOIST_LOSS,MPA_MOIST_LOSS_MAX,MPA_MOIST_MIN,DMPA_VEG, &
            MPA_CHAR,MPA_VEG,MPA_CHAR_MIN,MPA_VEG_MIN,MPA_VOLIT,MPA_VOLIT_LOSS_MAX,MPA_CHAR_LOSS,MPA_ASH
REAL(EB) :: DETA_VEG,ETA_H,ETAFM_VEG,ETAFP_VEG
REAL(EB) :: QCONF_L,Q_FOR_DRYING,Q_VEG_MOIST,Q_VEG_VOLIT,QNET_VEG,Q_FOR_VOLIT,Q_VOLIT,Q_UPTO_VOLIT
REAL(EB) :: C_DRAG,CM,CN,NUSS_HILPERT_CYL_FORCEDCONV,NUSS_MORGAN_CYL_FREECONV,HCON_VEG_FORCED,HCON_VEG_FREE, &
            LENGTH_SCALE,RAYLEIGH_NUM,ZGRIDCELL,ZGRIDCELL0
!LOGICAL  :: H_VERT_CYLINDER_LAMINAR,H_CYLINDER_RE

!INTEGER  :: IC,II,IOR,JJ,KK,IW_CELL

TYPE (WALL_TYPE),    POINTER :: WC =>NULL()
TYPE (SURFACE_TYPE), POINTER :: SF =>NULL()

!TYPE (WALL_TYPE),    POINTER :: WC1 =>NULL() !to handle qrad on slopes
!TYPE (SURFACE_TYPE), POINTER :: SF1 =>NULL() !to handle qrad on slopes

CALL POINT_TO_MESH(NM)

IF (VEG_LEVEL_SET_COUPLED .OR. VEG_LEVEL_SET_UNCOUPLED) RETURN

TMP_BOIL     = 373._EB
TMP_CHAR_MAX = 1300._EB
CP_ASH       = 800._EB !J/kg/K specific heat of ash
CP_H2O       = 4190._EB !J/kg/K specific heat of water
DT_BC     = T - VEG_CLOCK_BC
RDT_BC    = 1.0_EB/DT_BC

IF (N_REACTIONS>0) I_FUEL = REACTION(1)%FUEL_SMIX_INDEX

! Loop through vegetation wall cells and burn
!
VEG_WALL_CELL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
  WC  => WALL(IW)
  IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE VEG_WALL_CELL_LOOP

  SF  => SURFACE(WC%SURF_INDEX)
!
  IF (.NOT. SF%VEGETATION) CYCLE VEG_WALL_CELL_LOOP

  H_H2O_VEG = SF%VEG_H_H2O !J/kg
  H_PYR_VEG = SF%VEG_H_PYR !J/kg
  RH_PYR_VEG = 1._EB/H_PYR_VEG
  CHAR_FCTR  = 1._EB - SF%VEG_CHAR_FRACTION
  CHAR_FCTR2 = 1._EB/CHAR_FCTR

!Gas quantities 
  IIG = WC%ONE_D%IIG
  JJG = WC%ONE_D%JJG
  KKG = WC%ONE_D%KKG
  TMP_G = TMP(IIG,JJG,KKG)
  IF(SF%VEG_NO_BURN .OR. T <= DT_BC) WC%VEG_HEIGHT = SF%VEG_HEIGHT
! VEG_DRAG(IIG,JJG) = SF%VEG_DRAG_INI*(SF%VEG_CHARFRAC + CHAR_FCTR*WC%VEG_HEIGHT/SF%VEG_HEIGHT)

!Determine drag constant as a function of veg height (implemented in velo.f90)
! VEG_DRAG(IIG,JJG,:) = 0.0_EB
! VEG_DRAG(IIG,JJG,1) = SF%VEG_DRAG_INI*WC%VEG_HEIGHT/SF%VEG_HEIGHT

!-- Multi-layered approach following older, char fraction dependent, Cd
! VEG_DRAG(IIG,JJG,:) = 0.0_EB
! IF (WC%VEG_HEIGHT > 0.0_EB) THEN

!  IF (Z(1) >= SF%VEG_HEIGHT) THEN
!   VEG_DRAG(IIG,JJG,1) = SF%VEG_DRAG_INI*(SF%VEG_CHARFRAC + CHAR_FCTR*WC%VEG_HEIGHT/SF%VEG_HEIGHT)
!   VEG_DRAG(IIG,JJG,1) = VEG_DRAG(IIG,JJG,1)*WC%VEG_HEIGHT/Z(1)
!  ELSE
!   IF (WC%VEG_HEIGHT > Z(1))
!   DO KGRID=1,8
!    IF (Z(KGRID) < S
!    IF (WC%VEG_HEIGHT > Z(KGRID)) VEG_DRAG(IIG,JJG,KGRID) = SF%VEG_DRAG_INI
!    IF (WC%VEG_HEIGHT > Z(KGRID-1) .AND. WC%VEG_HEIGHT < Z(KGRID)) VEG_DRAG(IIG,JJG,KGRID) = & 
!       ((WC%VEG_HEIGHT-Z(KGRID-1))/(Z(KGRID)-Z(KGRID-1)))*SF%VEG_DRAG_INI
!   ENDDO
!  ENDIF

! ENDIF

!-- Drag varires with height above the terrain according to the fraction of the grid cell occupied by veg
!   Implemented in velo.f90 (KKG is the grid cell in the gas phase bordering the terrain).

  VEG_DRAG(IIG,JJG,:) = 0.0_EB
  IF (WC%VEG_HEIGHT > 0.0_EB) THEN
!if (nm==1 .and. iig==10 .and. jjg==10) print '(A,1x,3ES12.3)','vege: SF%VEG_HEIGHT,WC%VEG_HEIGHT,SF%VEG_DRAG_INI', &
!                                                             SF%VEG_HEIGHT,WC%VEG_HEIGHT,SF%VEG_DRAG_INI
    DO KGRID=0,8
      KLOC = KKG + KGRID
!if (nm==1 .and. iig==10 .and. jjg==10) print '(A,1x,3I3,1ES12.4)','vege: KGRID,KKG,KLOC,Z(KLOC)',kgrid,kkg,kloc,z(kloc)
      IF (Z(KLOC) <= WC%VEG_HEIGHT) THEN !grid cell filled with veg
        TMP_G = TMP(IIG,JJG,KLOC)
        RHO_GAS  = RHO(IIG,JJG,KLOC)
        ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(IIG,JJG,KLOC,1:N_TRACKED_SPECIES)
        CALL GET_VISCOSITY(ZZ_GET,MU_GAS,TMP_G)
        U2 = 0.25*(U(IIG,JJG,KLOC)+U(IIG-1,JJG,KLOC))**2
        V2 = 0.25*(V(IIG,JJG,KLOC)+V(IIG,JJG-1,KLOC))**2
        RE_VEG_PART = 4._EB*RHO_GAS*SQRT(U2 + V2 + W(IIG,JJG,KLOC)**2)/SF%VEG_SV/MU_GAS !for cylinder particle
        C_DRAG = 0.0_EB
        IF (RE_VEG_PART > 0.0_EB) C_DRAG = DRAG(RE_VEG_PART,2) !2 is for cylinder, 1 is for sphere
!if (nm==1 .and. iig==10 .and. jjg==10) print '(A,1ES12.4)','vege:filled cell C_DRAG',c_drag
        VEG_DRAG(IIG,JJG,KLOC)= C_DRAG*SF%VEG_DRAG_INI

      ENDIF

      IF (Z(KLOC) >  WC%VEG_HEIGHT .AND. Z(KLOC-1) < WC%VEG_HEIGHT) THEN !grid cell is partially filled with veg
        TMP_G = TMP(IIG,JJG,KLOC)
        RHO_GAS  = RHO(IIG,JJG,KLOC)
        ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(IIG,JJG,KLOC,1:N_TRACKED_SPECIES)
        CALL GET_VISCOSITY(ZZ_GET,MU_GAS,TMP_G)
        U2 = 0.25*(U(IIG,JJG,KLOC)+U(IIG-1,JJG,KLOC))**2
        V2 = 0.25*(V(IIG,JJG,KLOC)+V(IIG,JJG-1,KLOC))**2
        RE_VEG_PART = 4._EB*RHO_GAS*SQRT(U2 + V2 + W(IIG,JJG,KLOC)**2)/SF%VEG_SV/MU_GAS !for cylinder particle
        C_DRAG = 0.0_EB
        IF (RE_VEG_PART > 0.0_EB) C_DRAG = DRAG(RE_VEG_PART,2) !2 is for cylinder, 1 is for sphere
        VEG_DRAG(IIG,JJG,KLOC)= C_DRAG*SF%VEG_DRAG_INI*(WC%VEG_HEIGHT-Z(KLOC-1))/(Z(KLOC)-Z(KLOC-1))
!if (nm==1 .and. iig==10 .and. jjg==10) print '(A,1x,3ES12.4)','vege:unfilled cell Re,C_DRAG,volfrac',re_veg_part,c_drag, &
!                                                       (WC%VEG_HEIGHT-Z(KLOC-1))/(Z(KLOC)-Z(KLOC-1))
                                      
      ENDIF
!if (nm==1 .and. iig==10 .and. jjg==10) print '(A,1x,1I3,1ES12.4)','vege:K,VEG_DRAG(KLOC)',kgrid,veg_drag(iig,jjg,kloc)

    ENDDO

  ENDIF

  IF(SF%VEG_NO_BURN) CYCLE VEG_WALL_CELL_LOOP

! Initialize quantities
  Q_VEG_MOIST     = 0.0_EB
  Q_VEG_VOLIT     = 0.0_EB
  Q_UPTO_VOLIT    = 0.0_EB
  Q_VOLIT         = 0.0_EB
  Q_VEG_CHAR      = 0.0_EB
  MPA_MOIST_LOSS  = 0.0_EB
  MPA_VOLIT       = 0.0_EB
  MPA_CHAR_LOSS   = 0.0_EB
  SF%VEG_DIVQNET_L          = 0.0_EB
  SF%VEG_MOIST_FLUX_L       = 0.0_EB
  SF%VEG_FUEL_FLUX_L        = 0.0_EB
  WC%ONE_D%MASSFLUX(I_FUEL) = 0.0_EB 
  WC%ONE_D%QCONF            = 0.0_EB
  WC%LSET_FIRE       = .FALSE.
  IF (I_WATER > 0) WC%ONE_D%MASSFLUX(I_WATER) = 0.0_EB

! Vegetation variables and minimum bounds
  NVEG_L = SF%NVEG_L
  LBURN  = 0
!  Mininum bound on dry veg. Older approach, linear pyrolysis and no char
! MPA_VEG_MIN   = SF%VEG_CHARFRAC*SF%VEG_LOAD / REAL(NVEG_L,EB) !kg/m^2
!  Minimum bound on dry veg.Newer, linear or Arrhenius degradation and char
  MPA_VEG_MIN   = 0.001_EB*SF%VEG_LOAD/REAL(NVEG_L,EB) !kg/m^2

  MPA_CHAR_MIN  = SF%VEG_CHAR_FRACTION*MPA_VEG_MIN !kg/m^2
  MPA_MOIST_MIN = 0.0001_EB*SF%VEG_MOISTURE*SF%VEG_LOAD/REAL(NVEG_L,EB) !ks/m^2

  IF (SF%VEG_MOISTURE == 0.0_EB) MPA_MOIST_MIN = MPA_VEG_MIN
  DZVEG_L   = SF%VEG_HEIGHT/REAL(NVEG_L,EB)
  KAPPA_VEG = SF%VEG_KAPPA
  DETA_VEG  = DZVEG_L*KAPPA_VEG

! Find top of vegetation which burns downward from the top
  IF (SF%VEG_CHAR_OXIDATION) THEN
    DO IVEG_L = 1,NVEG_L 
      IF(WC%VEG_CHARMASS_L(IVEG_L) <= MPA_CHAR_MIN .AND. WC%VEG_FUELMASS_L(IVEG_L) <= MPA_VEG_MIN ) LBURN = IVEG_L
    ENDDO
  ELSE
    DO IVEG_L = 1,NVEG_L 
      IF(WC%VEG_FUELMASS_L(IVEG_L) <= MPA_VEG_MIN) LBURN = IVEG_L
    ENDDO
  ENDIF
! LBURN = 0
  WC%VEG_HEIGHT = REAL(NVEG_L-LBURN,EB)*DZVEG_L
! LBURN = 0 !keep charred veg
  !FIRELINE_MLR_MAX = w*R*(1-ChiChar)
  MPA_VOLIT_LOSS_MAX = SF%FIRELINE_MLR_MAX*DT_BC*DZVEG_L 
  MPA_MOIST_LOSS_MAX = MPA_VOLIT_LOSS_MAX

! Determine vertical gas-phase grid cell index for each vegetation layer. 
! This is needed for cases in which the vegetation height is larger than the height of the first grid cell
! The WC% and SF% indices are related by WC% goes from LBURN+1 to NVEG_L as SF% goes from 1 to NVEG_L to LBURN.
! Also, with increasing index value we pass from the top of the vegetation to the bottom

! SF%VEG_KGAS_L(:) = 0
! DO IVEG_L = 0, NVEG_L - LBURN
!  ZVEG = WC%VEG_HEIGHT - REAL(IVEG_L,EB)*DZVEG_L 
!  DO KGRID = 1,8
!    IF (ZVEG > Z(KGRID-1) .AND. ZVEG <= Z(KGRID)) SF%VEG_KGAS_L(IVEG_L)=KGRID
!  ENDDO
! ENDDO
!print*,'vege:kgas',SF%VEG_KGAS_L(:)

  DO IVEG_L = 1, NVEG_L - LBURN
   SF%VEG_KGAS_L(NVEG_L-LBURN-IVEG_L+1) = KKG 
   ZVEG = REAL(IVEG_L,EB)*DZVEG_L 
   ZGRIDCELL0 = 0.0_EB
   DO KGRID = 1,5
     ZGRIDCELL = ZGRIDCELL0 + Z(KKG+KGRID) - Z(KKG+KGRID-1)
     IF (ZVEG > ZGRIDCELL0 .AND. ZVEG <= ZGRIDCELL) SF%VEG_KGAS_L(NVEG_L-LBURN-IVEG_L+1) = KKG + KGRID - 1
     ZGRIDCELL0 = ZGRIDCELL
   ENDDO
  ENDDO


! Factors for computing divergence of incident and self emission radiant fluxes
! in vegetation fuel bed. These need to be recomputed as the height of the
! vegetation surface layer decreases with burning

! Factors for computing decay of +/- incident fluxes
  SF%VEG_FINCM_RADFCT_L(:) =  0.0_EB
  SF%VEG_FINCP_RADFCT_L(:) =  0.0_EB
! ETA_H = KAPPA_VEG*WC%VEG_HEIGHT
  ETA_H = KAPPA_VEG*REAL(NVEG_L-LBURN,EB)*DZVEG_L
  DO IVEG_L = 0,NVEG_L - LBURN
    ETAFM_VEG = REAL(IVEG_L,EB)*DETA_VEG
    ETAFP_VEG = ETA_H - ETAFM_VEG
    SF%VEG_FINCM_RADFCT_L(IVEG_L) = EXP(-ETAFM_VEG)
    SF%VEG_FINCP_RADFCT_L(IVEG_L) = EXP(-ETAFP_VEG)
  ENDDO

!  Integrand for computing +/- self emission fluxes
  SF%VEG_SEMISSP_RADFCT_L(:,:) = 0.0_EB
  SF%VEG_SEMISSM_RADFCT_L(:,:) = 0.0_EB
! q+
  DO IIVEG_L = 0,NVEG_L-LBURN !veg grid coordinate
    DO IVEG_L = IIVEG_L,NVEG_L-1-LBURN !integrand index
!    ETAG_VEG = IIVEG_L*DETA_VEG
!    ETAI_VEG =  IVEG_L*DETA_VEG
!    SF%VEG_SEMISSP_RADFCT_L(IVEG_L,IIVEG_L) = EXP(-(ETAI_VEG-ETAG_VEG))
     ETAFM_VEG = REAL((IVEG_L-IIVEG_L),EB)*DETA_VEG
     ETAFP_VEG = ETAFM_VEG + DETA_VEG
     SF%VEG_SEMISSP_RADFCT_L(IVEG_L,IIVEG_L) = EXP(-ETAFM_VEG) - EXP(-ETAFP_VEG)
    ENDDO
  ENDDO
! q-
  DO IIVEG_L = 0,NVEG_L-LBURN
    DO IVEG_L = 1,IIVEG_L
!    ETAG_VEG = IIVEG_L*DETA_VEG
!    ETAI_VEG =  IVEG_L*DETA_VEG
!    SF%VEG_SEMISSM_RADFCT_L(IVEG_L,IIVEG_L) = EXP(-(ETAG_VEG-ETAI_VEG))
     ETAFM_VEG = REAL((IIVEG_L-IVEG_L),EB)*DETA_VEG
     ETAFP_VEG = ETAFM_VEG + DETA_VEG
     SF%VEG_SEMISSM_RADFCT_L(IVEG_L,IIVEG_L) = EXP(-ETAFM_VEG) - EXP(-ETAFP_VEG)
    ENDDO
  ENDDO
!
! -----------------------------------------------
! compute CONVECTIVE HEAT FLUX on vegetation
! -----------------------------------------------
! Divergence of convective and radiative heat fluxes

  DO I=1,NVEG_L-LBURN
    KKG   = SF%VEG_KGAS_L(I)
    TMP_G = TMP(IIG,JJG,KKG)
    DTMP_L = TMP_G - WC%VEG_TMP_L(I+LBURN)

!Convective heat correlation for laminar flow (Holman see ref above) 
    IF (SF%VEG_HCONV_CYLLAM) H_CONV_L = 1.42_EB*(ABS(DTMP_L)/DZVEG_L)**0.25

!Convective heat correlation that accounts for air flow using forced convection correlation for
!a cylinder in a cross flow, Hilpert Correlation; Incropera & Dewitt Forth Edition p. 370
    IF(SF%VEG_HCONV_CYLRE) THEN 
     RHO_GAS  = RHO(IIG,JJG,KKG)
     TMP_FILM = 0.5_EB*(TMP_G + WC%VEG_TMP_L(I+LBURN))
     ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(IIG,JJG,KKG,1:N_TRACKED_SPECIES)
     CALL GET_VISCOSITY(ZZ_GET,MU_GAS,TMP_FILM)
     CALL GET_CONDUCTIVITY(ZZ_GET,K_GAS,TMP_G) !W/m/K
     U2 = 0.25*(U(IIG,JJG,KKG)+U(IIG-1,JJG,KKG))**2
     V2 = 0.25*(V(IIG,JJG,KKG)+V(IIG,JJG-1,KKG))**2
     RE_VEG_PART = 4._EB*RHO_GAS*SQRT(U2 + V2 + W(IIG,JJG,KKG)**2)/SF%VEG_SV/MU_GAS

     IF(RE_VEG_PART < 4._EB) THEN
       CN = 0.989_EB
       CM = 0.330_EB
     ELSE IF (RE_VEG_PART >= 4._EB .AND. RE_VEG_PART < 40._EB) THEN
       CN = 0.911_EB
       CM = 0.385_EB
     ELSE
       CN = 0.683_EB
       CM = 0.466_EB
     ENDIF
     H_CONV_L = 0.25_EB*SF%VEG_SV*K_GAS*CN*(RE_VEG_PART**CM)*PR_ONTH !W/K/m^2
    ENDIF
!
! Use largest of natural and forced convective heat transfer
   
    IF(SF%VEG_HCONV_CYLMAX) THEN 
      RHO_GAS  = RHO(IIG,JJG,KKG)
      TMP_FILM = 0.5_EB*(TMP_G + WC%VEG_TMP_L(I+LBURN))
      ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(IIG,JJG,KKG,1:N_TRACKED_SPECIES)
      CALL GET_VISCOSITY(ZZ_GET,MU_GAS,TMP_G)
      CALL GET_SPECIFIC_HEAT(ZZ_GET,CP_GAS,TMP_G)
      CALL GET_CONDUCTIVITY(ZZ_GET,K_GAS,TMP_G)
      U2 = 0.25*(U(IIG,JJG,KKG)+U(IIG-1,JJG,KKG))**2
      V2 = 0.25*(V(IIG,JJG,KKG)+V(IIG,JJG-1,KKG))**2
      RE_VEG_PART = 4._EB*RHO_GAS*SQRT(U2 + V2 + W(IIG,JJG,KKG)**2)/SF%VEG_SV/MU_GAS !for cylinder SV

! - Forced convection heat transfer coefficients in a layer
!
! Hilpert Correlation (Incropera & DeWitt Fourth Edition, p. 370) for cylinder in crossflow,
! forced convection
      IF(RE_VEG_PART < 4._EB) THEN
        CN = 0.989_EB
        CM = 0.330_EB
      ELSE IF (RE_VEG_PART >= 4._EB .AND. RE_VEG_PART < 40._EB) THEN
        CN = 0.911_EB
        CM = 0.385_EB
      ELSE
        CN = 0.683_EB
        CM = 0.466_EB
      ENDIF
      NUSS_HILPERT_CYL_FORCEDCONV = CN*(RE_VEG_PART**CM)*PR_ONTH !Nusselt number
      HCON_VEG_FORCED = 0.25_EB*SF%VEG_SV*K_GAS*NUSS_HILPERT_CYL_FORCEDCONV !W/m^2 from Hilpert (cylinder)

! - Free convection heat transfer coefficients
      LENGTH_SCALE = 4._EB/SF%VEG_SV !horizontal cylinder diameter
      RAYLEIGH_NUM = 9.8_EB*ABS(DTMP_L)*LENGTH_SCALE**3*RHO_GAS**2*CP_GAS/(TMP_FILM*MU_GAS*K_GAS)

! Morgan correlation (Incropera & DeWitt, 4th Edition, p. 501-502) for horizontal cylinder, free convection
      IF (RAYLEIGH_NUM < 0.01_EB) THEN
        CN = 0.675_EB
        CM = 0.058_EB
      ELSE IF (RAYLEIGH_NUM >= 0.01_EB .AND. RAYLEIGH_NUM < 100._EB) THEN
        CN = 1.02_EB
        CM = 0.148_EB
      ELSE IF (RAYLEIGH_NUM >= 100._EB .AND. RAYLEIGH_NUM < 10**4._EB) THEN
        CN = 0.85_EB
        CM = 0.188_EB
      ELSE IF (RAYLEIGH_NUM >= 10**4._EB .AND. RAYLEIGH_NUM < 10**7._EB) THEN
        CN = 0.48_EB
        CM = 0.25_EB
      ELSE IF (RAYLEIGH_NUM >= 10**7._EB .AND. RAYLEIGH_NUM < 10**12._EB) THEN
        CN = 0.125_EB
        CM = 0.333_EB
      ENDIF

      NUSS_MORGAN_CYL_FREECONV = CN*RAYLEIGH_NUM**CM
      HCON_VEG_FREE = 0.25_EB*SF%VEG_SV*K_GAS*NUSS_MORGAN_CYL_FREECONV !W/K/m^2

      H_CONV_L = MAX(HCON_VEG_FORCED,HCON_VEG_FREE)
    ENDIF

    QCONF_L  = H_CONV_L*DTMP_L
    SF%VEG_DIVQNET_L(I) = SF%VEG_PACKING*SF%VEG_SV*QCONF_L*DZVEG_L !W/m^2

  ENDDO
!
  WC%ONE_D%QCONF = SUM(SF%VEG_DIVQNET_L) !*RDN(IW)*WC%VEG_HEIGHT
! qconf(iw) = 0.0_EB
!
! -----------------------------------------------
! Compute +/- radiation fluxes and their divergence due to self emission within vegetation
! -----------------------------------------------
  LAYER_RAD_FLUXES: IF (LBURN < NVEG_L) THEN
    VEG_QRP_EMISS   = 0.0_EB ; VEG_QRM_EMISS = 0.0_EB 
    VEG_QRNET_EMISS = 0.0_EB ; VEG_DIV_QRNET_EMISS = 0.0_EB
! qe+
    DO J=0,NVEG_L-LBURN !veg grid coordinate loop
      DO I=J,NVEG_L-LBURN !integrand loop 
         VEG_QRP_EMISS(J) =  VEG_QRP_EMISS(J) + SF%VEG_SEMISSP_RADFCT_L(I,J)*WC%VEG_TMP_L(I+LBURN)**4
      ENDDO
    ENDDO
! qe-
    DO J=0,NVEG_L-LBURN  !veg grid coordinate
      DO I=0,J           !integrand for q-
         VEG_QRM_EMISS(J) = VEG_QRM_EMISS(J) + SF%VEG_SEMISSM_RADFCT_L(I,J)*WC%VEG_TMP_L(I+LBURN)**4
      ENDDO
    ENDDO
    VEG_QRP_EMISS =  VEG_QRP_EMISS*SIGMA
    VEG_QRM_EMISS =  VEG_QRM_EMISS*SIGMA
!
    DO I=0,NVEG_L-LBURN
      VEG_QRNET_EMISS(I) = VEG_QRP_EMISS(I)-VEG_QRM_EMISS(I)
    ENDDO
!    DO I=1,NVEG_L-LBURN
!      VEG_QRNET_EMISS(I)  = VEG_QRNET_EMISS(I) - VEG_QRM_EMISS(I)
!    ENDDO
!
    DO I=1,NVEG_L-LBURN
      VEG_DIV_QRNET_EMISS(I) = VEG_QRNET_EMISS(I-1) - VEG_QRNET_EMISS(I)
    ENDDO
!
! Compute +/- radiation fluxes and their divergence due to incident fluxes on boundaries
    QRADM_INC = WC%ONE_D%QRADIN/SF%EMISSIVITY !sigma*Ta^4 + flame
!   QRADM_INC = QRADIN(IW)/E_WALL(IW) + SIGMA*TMP_F(IW)**4 ! as done in FDS4
!   print*,'vege: QRADIN(IW)',qradin(iw)

! Adjust incident radiant flux to account for sloped terrain
! assumes user put VEG_NO_BURN=.TRUE. for vertical faces
! sets qrad on cell downspread of vertical face = qrad on cell face upspread of vertical face
!   QRADM_INC = QRADM_INC*1.0038_EB !adjustment for horizontal faces assuming 5 degree slope
!   II = WC%II
!   JJ = WC%JJ
!   KK = WC%KK
!   IC = CELL_INDEX(II-1,JJ,KK)
!   IOR = 1
!   IW_CELL = WALL_INDEX(IC,IOR) 
!   WC1 => WALL(IW_CELL)
!   SF1 => SURFACE(WC1%SURF_INDEX)
!print*,'vege: i,j,k,iw,sf',ii,jj,kk,iw,sf1%veg_no_burn
!   IF(SF1%VEG_NO_BURN) THEN
!print*,'vege: in vertical face qrad determination'
!!   QRADM_INC_SLOPE_VERTFACE = QRADM_INC_SLOPE_VERTFACE + WALL(IW_CELL)%RADIN/WALL(IW_CELL)%E_WALL
!!   QRADM_INC_SLOPE_VERTFACE = QRADM_INC_SLOPE_VERTFACE*0.0872_EB !assumes 5 degree slope
!!   QRADM_INC = QRADM_INC + QRADM_INC_SLOPE_VERTFACE !adjustment for adjacent vertical faces

!   IOR = -3
!   IW_CELL = WALL_INDEX(IC,IOR)
!adjustment for horizontal faces downspread of vertical face
!set flux = to max of flux up or downspread 
!print*,'vege: i,j,k,iw,qr',ii,jj,kk,wall(iw_cell)%qradin,wall(iw)%qradin
!   WALL(IW)%QRADIN = MAX(WALL(IW_CELL)%QRADIN,WALL(IW)%QRADIN) 
!   QRADM_INC = 1.0038_EB*WALL(IW)%QRADIN/WALL(IW_CELL)%E_WALL !assumes 5 degree slope!!!
!print*,'vege: qradm_inc,wallqrad',qradm_inc,wall(iw)%qradin
!   ENDIF

    ETAVEG_H  = (NVEG_L - LBURN)*DETA_VEG
    !this QRADP_INC ensures zero net radiant fluxes at bottom of vegetation
    IF(SF%VEG_GROUND_ZERO_RAD) QRADP_INC = QRADM_INC*SF%VEG_FINCM_RADFCT_L(NVEG_L-LBURN) + VEG_QRM_EMISS(NVEG_L-LBURN)
    !this QRADP_INC assumes the ground stays at user specified temperature
    IF(.NOT. SF%VEG_GROUND_ZERO_RAD) QRADP_INC = SIGMA*SF%VEG_GROUND_TEMP**4
!   QRADP_INC = SIGMA*WC%VEG_TMP_L(NVEG_L)**4 
!   IF(.NOT. SF%VEG_GROUND_ZERO_RAD) QRADP_INC = SIGMA*TMP_G**4
!   QRADP_INC = SIGMA*WC%VEG_TMP_L(NVEG_L)**4*EXP(-ETAVEG_H) + VEG_QRM_EMISS(NVEG_L-LBURN) !fds4
    VEG_QRM_INC   = 0.0_EB ; VEG_QRP_INC = 0.0_EB 
    VEG_QRNET_INC = 0.0_EB ; VEG_DIV_QRNET_INC = 0.0_EB
    DO I=0,NVEG_L-LBURN
      VEG_QRM_INC(I)   = QRADM_INC*SF%VEG_FINCM_RADFCT_L(I)
      VEG_QRP_INC(I)   = QRADP_INC*SF%VEG_FINCP_RADFCT_L(I)
      VEG_QRNET_INC(I) = VEG_QRP_INC(I)-VEG_QRM_INC(I)
    ENDDO
    DO I=1,NVEG_L-LBURN
      VEG_DIV_QRNET_INC(I) = VEG_QRNET_INC(I-1) - VEG_QRNET_INC(I)
    ENDDO
  ENDIF LAYER_RAD_FLUXES
!
! Add divergence of net radiation flux to divergence of convection flux
  DO I=1,NVEG_L-LBURN
    SF%VEG_DIVQNET_L(I)= SF%VEG_DIVQNET_L(I) - (VEG_DIV_QRNET_INC(I) + VEG_DIV_QRNET_EMISS(I)) !includes self emiss
!   SF%VEG_DIVQNET_L(I)= SF%VEG_DIVQNET_L(I) - VEG_DIV_QRNET_INC(I) !no self emission contribution
  ENDDO
!
!
!      ************** Boundary Fuel Non-Arrehnius (Linear in temp) Degradation model *************************
! Drying occurs if qnet > 0 with Tveg held at 100 c
! Pyrolysis occurs according to Morvan & Dupuy empirical formula. Linear
! temperature dependence with qnet factor
!

  IF_VEG_DEGRADATION_LINEAR: IF (SF%VEG_DEGRADATION == 'LINEAR') THEN

    LAYER_LOOP1: DO IVEG_L = LBURN+1,NVEG_L
!
! Compute temperature of vegetation
!
      MPA_CHAR    = WC%VEG_CHARMASS_L(IVEG_L)
      MPA_VEG     = WC%VEG_FUELMASS_L(IVEG_L)
      MPA_MOIST   = WC%VEG_MOISTMASS_L(IVEG_L)
      TMP_VEG     = WC%VEG_TMP_L(IVEG_L)
      QNET_VEG    = SF%VEG_DIVQNET_L(IVEG_L-LBURN)
      CP_VEG      = (0.01_EB + 0.0037_EB*TMP_VEG)*1000._EB !J/kg/K
      CP_CHAR     = 420._EB + 2.09_EB*TMP_VEG + 6.85E-4_EB*TMP_VEG**2 !J/kg/K Park etal. C&F 2010 147:481-494
      CP_TOTAL    = CP_H2O*MPA_MOIST +  CP_VEG*MPA_VEG + CP_CHAR*MPA_CHAR
      DTMP_VEG    = DT_BC*QNET_VEG/CP_TOTAL
      TMP_VEG_NEW = TMP_VEG + DTMP_VEG 

      IF_DIVQ_L_GE_0: IF(QNET_VEG > 0._EB) THEN 

! -- drying of veg layer
      IF(MPA_MOIST > MPA_MOIST_MIN .AND. TMP_VEG_NEW >= TMP_BOIL) THEN
        Q_FOR_DRYING   = (TMP_VEG_NEW - TMP_BOIL)/DTMP_VEG * QNET_VEG
        MPA_MOIST_LOSS = MIN(DT_BC*Q_FOR_DRYING/H_H2O_VEG,MPA_MOIST_LOSS_MAX)
        MPA_MOIST_LOSS = MIN(MPA_MOIST_LOSS,MPA_MOIST-MPA_MOIST_MIN)
        TMP_VEG_NEW    = TMP_BOIL
        WC%VEG_MOISTMASS_L(IVEG_L) = MPA_MOIST - MPA_MOIST_LOSS !kg/m^2
        IF( WC%VEG_MOISTMASS_L(IVEG_L) <= MPA_MOIST_MIN ) WC%VEG_MOISTMASS_L(IVEG_L) = 0.0_EB
        IF (I_WATER > 0) WC%ONE_D%MASSFLUX(I_WATER) = WC%ONE_D%MASSFLUX(I_WATER) + RDT_BC*MPA_MOIST_LOSS
!       WC%VEG_TMP_L(IVEG_L) = TMP_VEG_NEW
      ENDIF

! -- pyrolysis multiple layers
      IF_VOLITIZATION: IF (MPA_MOIST <= MPA_MOIST_MIN) THEN

!Newer version, includes char which becomes a heat sink since its not oxidized
        IF(TMP_VEG_NEW >= 400._EB .AND. MPA_VEG > MPA_VEG_MIN) THEN
          Q_UPTO_VOLIT = (CP_VEG*MPA_VEG + CP_CHAR*MPA_CHAR)*(400._EB-TMP_VEG)
!         Q_UPTO_VOLIT = CP_VEG*MPA_VEG*(400._EB-TMP_VEG)
          Q_FOR_VOLIT  = DT_BC*QNET_VEG - Q_UPTO_VOLIT
          Q_VOLIT      = Q_FOR_VOLIT*0.01_EB*(TMP_VEG-400._EB)

          MPA_VOLIT    = CHAR_FCTR*Q_VOLIT*RH_PYR_VEG
          MPA_VOLIT    = MAX(MPA_VOLIT,0._EB)
          MPA_VOLIT    = MIN(MPA_VOLIT,MPA_VOLIT_LOSS_MAX) !user specified max

          DMPA_VEG     = CHAR_FCTR2*MPA_VOLIT
          DMPA_VEG     = MIN(DMPA_VEG,(MPA_VEG-MPA_VEG_MIN))
          MPA_VEG      = MPA_VEG - DMPA_VEG

          MPA_VOLIT    = CHAR_FCTR*DMPA_VEG
          MPA_CHAR     = MPA_CHAR + SF%VEG_CHAR_FRACTION*DMPA_VEG
          Q_VOLIT      = MPA_VOLIT*H_PYR_VEG 

          TMP_VEG_NEW  = TMP_VEG + (Q_FOR_VOLIT-Q_VOLIT)/(MPA_VEG*CP_VEG + MPA_CHAR*CP_CHAR)
!         TMP_VEG_NEW  = TMP_VEG + (Q_FOR_VOLIT-Q_VOLIT)/(MPA_VEG*CP_VEG)
          TMP_VEG_NEW  = MIN(TMP_VEG_NEW,500._EB)
          WC%VEG_CHARMASS_L(IVEG_L) = MPA_CHAR
          WC%VEG_FUELMASS_L(IVEG_L) = MPA_VEG
          IF( WC%VEG_FUELMASS_L(IVEG_L) <= MPA_VEG_MIN ) WC%VEG_FUELMASS_L(IVEG_L) = 0.0_EB !**
          WC%ONE_D%MASSFLUX(I_FUEL)= WC%ONE_D%MASSFLUX(I_FUEL) + RDT_BC*MPA_VOLIT
        ENDIF        

!Older version (used for Tony's Catchpole paper) assumes that no char is present in degradation and
!pyrolysis converts a total of (1-Xchar)*(original dry mass) of mass to fuel vapor. Ignores
!heat sink influence of char.
!       IF(TMP_VEG_NEW >= 400._EB .AND. MPA_VEG > MPA_VEG_MIN) THEN
!         Q_UPTO_VOLIT = CP_VEG*MPA_VEG*(400._EB-TMP_VEG)
!         Q_FOR_VOLIT  = DT_BC*QNET_VEG - Q_UPTO_VOLIT
!         Q_VOLIT      = Q_FOR_VOLIT*0.01_EB*(TMP_VEG-400._EB)
!         MPA_VOLIT    = CHAR_FCTR*Q_VOLIT*0.00000239_EB
!         MPA_VOLIT    = MAX(MPA_VOLIT,0._EB)
!         MPA_VOLIT    = MIN(MPA_VOLIT,MPA_VOLIT_LOSS_MAX)
!         MPA_VOLIT    = MIN(MPA_VOLIT,MPA_VEG-MPA_VEG_MIN)
!         MPA_VEG      = MPA_VEG - MPA_VOLIT
!         Q_VOLIT      = MPA_VOLIT*418000._EB
!         TMP_VEG_NEW  = TMP_VEG + (Q_FOR_VOLIT-Q_VOLIT)/(MPA_VEG*CP_VEG)
!         TMP_VEG_NEW  = MIN(TMP_VEG_NEW,500._EB)
!         WC%VEG_FUELMASS_L(IVEG_L) = MPA_VEG
!         WC%ONE_D%MASSFLUX(I_FUEL)= WC%ONE_D%MASSFLUX(I_FUEL) + RDT_BC*MPA_VOLIT
!       ENDIF        


      ENDIF IF_VOLITIZATION

      ENDIF IF_DIVQ_L_GE_0
      
      IF(MPA_VEG <= MPA_VEG_MIN) TMP_VEG_NEW = TMP_G
      WC%VEG_TMP_L(IVEG_L) = TMP_VEG_NEW

    ENDDO LAYER_LOOP1

!   WC%VEG_TMP_L(LBURN) = WC%VEG_TMP_L(LBURN+1)
    WC%VEG_TMP_L(LBURN) = TMP_G

  ENDIF  IF_VEG_DEGRADATION_LINEAR

!      ************** Boundary Fuel Arrehnius Degradation model *************************
! Drying and pyrolysis occur according to Arrehnius expressions obtained 
! from the literature (Porterie et al., Num. Heat Transfer, 47:571-591, 2005
! Predicting wildland fire behavior and emissions using a fine-scale physical
! model

  IF_VEG_DEGRADATION_ARRHENIUS: IF(SF%VEG_DEGRADATION == 'ARRHENIUS') THEN
!   A_H2O_VEG      = 600000._EB !1/s sqrt(K)
!   E_H2O_VEG      = 5800._EB !K

!   A_PYR_VEG      = 36300._EB !1/s
!   E_PYR_VEG      = 7250._EB !K

!   A_CHAR_VEG     = 430._EB !m/s
!   E_CHAR_VEG     = 9000._EB !K
!   H_CHAR_VEG     = -12.0E+6_EB !J/kg

!   BETA_CHAR_VEG  = 0.2_EB
!   NU_CHAR_VEG    = SF%VEG_CHAR_FRACTION
!   NU_ASH_VEG     = 0.1_EB
!   NU_O2_CHAR_VEG = 1.65_EB
!   CHAR_ENTHALPY_FRACTION_VEG = 0.5_EB
!print*,'-----------------------------'
!print 1115,beta_char_veg,nu_char_veg,nu_ash_veg,nu_o2_char_veg,char_enthalpy_fraction_veg

    A_H2O_VEG      = SF%VEG_A_H2O !1/2 sqrt(K)
    E_H2O_VEG      = SF%VEG_E_H2O !K

    A_PYR_VEG      = SF%VEG_A_PYR !1/s
    E_PYR_VEG      = SF%VEG_E_PYR !K

    A_CHAR_VEG     = SF%VEG_A_CHAR !m/s
    E_CHAR_VEG     = SF%VEG_E_CHAR !K
    H_CHAR_VEG     = SF%VEG_H_CHAR !J/kg

    BETA_CHAR_VEG  = SF%VEG_BETA_CHAR
    NU_CHAR_VEG    = SF%VEG_CHAR_FRACTION
    NU_ASH_VEG     = SF%VEG_ASH_FRACTION/SF%VEG_CHAR_FRACTION !fraction of char that can become ash
!   NU_ASH_VEG     = 0.1_EB !fraction of char that can become ash, Porterie et al. 2005 Num. Heat Transfer
    NU_O2_CHAR_VEG = SF%VEG_NU_O2_CHAR
    CHAR_ENTHALPY_FRACTION_VEG = SF%VEG_CHAR_ENTHALPY_FRACTION
!print 1115,nu_ash_veg,sf%veg_ash_fraction,sf%veg_char_fraction
!1115 format('vege:',2x,3(e15.5))

    LAYER_LOOP2: DO IVEG_L = LBURN+1,NVEG_L

      MPA_MOIST = WC%VEG_MOISTMASS_L(IVEG_L)
      MPA_VEG   = WC%VEG_FUELMASS_L(IVEG_L)
      MPA_CHAR  = WC%VEG_CHARMASS_L(IVEG_L)
      MPA_ASH   = WC%VEG_ASHMASS_L(IVEG_L)
      TMP_VEG   = WC%VEG_TMP_L(IVEG_L)

      TEMP_THRESEHOLD: IF (WC%VEG_TMP_L(IVEG_L) > 323._EB) THEN
              !arbitrary thresehold to prevent low-temp hrr reaction
              !added for drainage runs

! Drying of vegetation (Arrhenius)
      IF_DEHYDRATION_2: IF (MPA_MOIST > MPA_MOIST_MIN) THEN
        MPA_MOIST_LOSS = MIN(DT_BC*MPA_MOIST*A_H2O_VEG*EXP(-E_H2O_VEG/TMP_VEG)/SQRT(TMP_VEG), &
                         MPA_MOIST-MPA_MOIST_MIN)
        MPA_MOIST_LOSS = MIN(MPA_MOIST_LOSS,MPA_MOIST_LOSS_MAX) !user specified max
        MPA_MOIST      = MPA_MOIST - MPA_MOIST_LOSS
        WC%VEG_MOISTMASS_L(IVEG_L) = MPA_MOIST !kg/m^2
        IF (MPA_MOIST <= MPA_MOIST_MIN) WC%VEG_MOISTMASS_L(IVEG_L) = 0.0_EB
!print 1114,iveg_l,iig,jjg,tmp_veg,mpa_moist,mpa_moist_loss,dt_bc
!1114 format('(vege)',1x,3(I3),2x,4(e15.5))
!print*,'wwwwwwwwwwwwwwwwwwwww'
      ENDIF IF_DEHYDRATION_2

! Volitalization of vegetation(Arrhenius)
      IF_VOLITALIZATION_2: IF(MPA_VEG > MPA_VEG_MIN) THEN
        MPA_VOLIT = MAX(CHAR_FCTR*DT_BC*MPA_VEG*A_PYR_VEG*EXP(-E_PYR_VEG/TMP_VEG),0._EB)
        MPA_VOLIT = MIN(MPA_VOLIT,MPA_VOLIT_LOSS_MAX) !user specified max

        DMPA_VEG = CHAR_FCTR2*MPA_VOLIT
        DMPA_VEG = MIN(DMPA_VEG,(MPA_VEG - MPA_VEG_MIN))
        MPA_VEG  = MPA_VEG - DMPA_VEG

        MPA_VOLIT = CHAR_FCTR*DMPA_VEG
        MPA_CHAR  = MPA_CHAR + SF%VEG_CHAR_FRACTION*DMPA_VEG !kg/m^2
!print 1114,iveg_l,iig,jjg,tmp_veg,mpa_veg,mpa_volit,dt_bc
!print*,'vvvvvvvvvvvvvvvvvvvvv'

      ENDIF IF_VOLITALIZATION_2

      WC%VEG_FUELMASS_L(IVEG_L) = MPA_VEG
      WC%VEG_CHARMASS_L(IVEG_L) = MPA_CHAR

      WC%ONE_D%MASSFLUX(I_FUEL)= WC%ONE_D%MASSFLUX(I_FUEL) + MPA_VOLIT*RDT_BC
      IF (I_WATER > 0) WC%ONE_D%MASSFLUX(I_WATER) = WC%ONE_D%MASSFLUX(I_WATER) + MPA_MOIST_LOSS*RDT_BC

!Char oxidation oF Vegetation Layer within the Arrhenius pyrolysis model
!(note that this can be handled only approximately with the conserved
!scalar based gas-phase combustion model - no gas phase oxygen is consumed by
!the char oxidation reaction since it would be inconsistent with the state
!relation for oxygen based on the conserved scalar approach for gas phase
!combustion)
      IF_CHAR_OXIDATION: IF (SF%VEG_CHAR_OXIDATION .AND. MPA_CHAR > 0.0_EB) THEN
         ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(IIG,JJG,KKG,1:N_TRACKED_SPECIES)
         CALL GET_MASS_FRACTION(ZZ_GET,O2_INDEX,Y_O2)
         CALL GET_VISCOSITY(ZZ_GET,MU_GAS,TMP_G)
         RE_D = RHO_GAS*SQRT(U2 + V2 + W(IIG,JJG,1)**2)*4._EB/SF%VEG_SV/MU_GAS 
         MPA_CHAR_LOSS = DT_BC*RHO_GAS*Y_O2*A_CHAR_VEG/NU_O2_CHAR_VEG*SF%VEG_SV*  &
                         SF%VEG_PACKING*EXP(-E_CHAR_VEG/WC%VEG_TMP_L(IVEG_L))*  &
                         (1+BETA_CHAR_VEG*SQRT(RE_D))
         MPA_CHAR_LOSS = MIN(MPA_CHAR,MPA_CHAR_LOSS)
         MPA_CHAR      = MPA_CHAR - MPA_CHAR_LOSS
         MPA_ASH       = MPA_ASH + NU_ASH_VEG*MPA_CHAR_LOSS
!        MPA_CHAR_CO2  = (1._EB + NU_O2_CHAR_VEG - NU_ASH_VEG)*MPA_CHAR_LOSS
         WC%VEG_CHARMASS_L(IVEG_L) = MPA_CHAR !kg/m^3
         WC%VEG_ASHMASS_L(IVEG_L)  = MPA_ASH

         IF (MPA_CHAR <= MPA_CHAR_MIN .AND. MPA_VEG <= MPA_VEG_MIN) WC%VEG_CHARMASS_L(IVEG_L) = 0.0_EB
       ENDIF IF_CHAR_OXIDATION

      ENDIF TEMP_THRESEHOLD

! Vegetation temperature (Arrhenius)
      CP_VEG = (0.01_EB + 0.0037_EB*TMP_VEG)*1000._EB !W/kg/K
      CP_CHAR= 420._EB + 2.09_EB*TMP_VEG + 6.85E-4_EB*TMP_VEG**2 !J/kg/K Park etal. C&F 2010 147:481-494
      Q_VEG_CHAR       = MPA_CHAR_LOSS*H_CHAR_VEG
      CP_MOIST_AND_VEG = CP_H2O*WC%VEG_MOISTMASS_L(IVEG_L) + CP_VEG*WC%VEG_FUELMASS_L(IVEG_L) + &
                         CP_CHAR*WC%VEG_CHARMASS_L(IVEG_L) + CP_ASH*WC%VEG_ASHMASS_L(IVEG_L)

      WC%VEG_TMP_L(IVEG_L) = WC%VEG_TMP_L(IVEG_L) + (DT_BC*SF%VEG_DIVQNET_L(IVEG_L-LBURN) - &
                             (MPA_MOIST_LOSS*H_H2O_VEG + MPA_VOLIT*H_PYR_VEG) + CHAR_ENTHALPY_FRACTION_VEG*Q_VEG_CHAR ) &
                             /CP_MOIST_AND_VEG
      WC%VEG_TMP_L(IVEG_L) = MAX( WC%VEG_TMP_L(IVEG_L), TMPA)
      WC%VEG_TMP_L(IVEG_L) = MIN( WC%VEG_TMP_L(IVEG_L), TMP_CHAR_MAX)

    ENDDO LAYER_LOOP2

  ENDIF IF_VEG_DEGRADATION_ARRHENIUS
  
!  if (wc%veg_tmp_L(lburn+1) > 300._EB) wc%massflux(i_fuel)=0.1
!  if (iig == 14 .or. iig==15) wc%massflux(i_fuel)=0.1
  WC%VEG_TMP_L(LBURN) = MAX(TMP_G,TMPA)
  WC%ONE_D%MASSFLUX_SPEC(I_FUEL) = WC%ONE_D%MASSFLUX(I_FUEL)
  IF (I_WATER > 0) WC%ONE_D%MASSFLUX_SPEC(I_WATER) = WC%ONE_D%MASSFLUX(I_WATER)
 
! Temperature boundary condtions 
! Mass boundary conditions are determine in subroutine SPECIES_BC in wall.f90 for case SPECIFIED_MASS_FLUX
! TMP_F(IW) = WC%VEG_TMP_L(NVEG_L)
! IF (LBURN < NVEG_L)  TMP_F(IW) = WC%VEG_TMP_L(1+LBURN)
  IF (LBURN < NVEG_L) THEN
    WC%ONE_D%TMP_F = WC%VEG_TMP_L(1+LBURN)
if (wc%one_d%tmp_f < 0.0_EB) print '(A,1x,1ES16.8)','vege:wall(iw)%tmp_f',wc%one_d%tmp_f
!   TMP_F(IW) = ((VEG_QRP_INC(0)+VEG_QRP_EMISS(0))/SIGMA)**.25 !as done in FDS4
  ELSE
    WC%ONE_D%TMP_F = MAX(TMP_G,TMPA) !Tveg=Tgas if veg is completely burned
!   TMP_F(IW) = TMPA  !Tveg=Tambient if veg is completely burned
  ENDIF
! TMP_F(IW) = MAX(TMP_F(IW),TMPA)

ENDDO VEG_WALL_CELL_LOOP

VEG_CLOCK_BC = T

END SUBROUTINE BNDRY_VEG_MASS_ENERGY_TRANSFER
!
!************************************************************************************************
!
!\/\/\////\/\/\/\/\/\/\\\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
!\/\/\////\/\/\/\/\/\/\\\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
SUBROUTINE INITIALIZE_LEVEL_SET_FIREFRONT(DT,NM)
!\/\/\////\/\/\/\/\/\/\\\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
!\/\/\////\/\/\/\/\/\/\\\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
!
! Level set based modeling of fire front propatagion across terrain. 
! There are four implementations from the simplest in which the wind is constant
! in directior and magnitude to a CFD coupled implementation with buoynacy generated flow.
! See User Guide on Google Docs
!
! Issues:
! 1) Need to make level set computation mesh dependent so the the LS slice file
!    is created only where fire is expected
!
!
INTEGER, INTENT(IN) :: NM
INTEGER  :: I,IM1,IM2,IIG,IP1,IP2,IW,J,JJG,JM1,JP1,KKG
REAL(EB), INTENT(INOUT) :: DT
REAL(EB) :: LX,SR_MAX,UMAX_LS,VMAX_LS
REAL(EB) :: G_EAST,G_WEST,G_SOUTH,G_NORTH
REAL(EB) :: VERT_CANOPY_EXTENT

REAL(EB), ALLOCATABLE, DIMENSION(:) :: X_LS,Y_LS

CHARACTER(30) :: CFORM,SMOKEVIEW_LABEL,SMOKEVIEW_BAR_LABEL,UNITS

REAL(EB), POINTER, DIMENSION(:,:) :: ZT => NULL()

TYPE (WALL_TYPE),    POINTER :: WC =>NULL()
TYPE (SURFACE_TYPE), POINTER :: SF =>NULL()


CALL CPU_TIME(CPUTIME)
LS_T_BEG = CPUTIME

CALL POINT_TO_MESH(NM)

ZT => LS_Z_TERRAIN

!print*,'vege: in initialize LS'
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
 
!******************* Initialize time stepping and other constants

SUMTIME_LS = 0.0_EB ! Used for time step output

SUM_T_SLCF_LS = 0._EB
DT_COEF_LS = 0.5_EB
TIME_FLANKFIRE_QUENCH = 20.0_EB !flankfire lifetime in seconds

LSET_ELLIPSE = .FALSE. ! Default value of flag for the elliptical spread model
LSET_TAN2    = .FALSE. ! Default value: Flag for ROS proportional to Tan(slope)^2 
!HEAD_WIDTH   = 1.0_EB

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
LIMITER_LS = FLUX_LIMITER
IF (LIMITER_LS > 3) LIMITER_LS = 1

!******************* Open output files and write headers (put this in dump.f90 ASSIGN_FILE_NAMES)
   ! Slice Filenames

!  DO N=1,M%N_SLCF
!     LU_SLCF(N,NM) = GET_FILE_NUMBER()
!     IF (NMESHES>1) THEN
!        IF (M%N_SLCF <100) CFORM = '(A,A,I4.4,A,I2.2,A)'
!        IF (M%N_SLCF>=100) CFORM = '(A,A,I4.4,A,I3.3,A)'
!        WRITE(FN_SLCF(N,NM),CFORM) TRIM(CHID),'_',NM,'_',N,'.sf'
!     ELSE
!        IF (M%N_SLCF <100) CFORM = '(A,A,I2.2,A)'
!        IF (M%N_SLCF>=100) CFORM = '(A,A,I3.3,A)'
!        WRITE(FN_SLCF(N,NM),CFORM) TRIM(CHID),'_',N,'.sf'
!     ENDIF
!  ENDDO


TIME_LS    = T_BEGIN

!******************* Assign filenames, open data files and write headers; Put filenames & case info in smv file 
!                    (put this in dump.f90 ASSIGN_FILE_NAMES)
IF (NMESHES>1) THEN
  CFORM = '(A,A,I4.4,A,A)'
ELSE
  CFORM = '(A,A,A)'
ENDIF

!--Level set field for animation via Smokeview
LU_SLCF_LS(1) = GET_FILE_NUMBER()
SMOKEVIEW_LABEL = 'phifield'
SMOKEVIEW_BAR_LABEL = 'phifield'
UNITS  = 'C'
IF(NMESHES  > 1) WRITE(FN_SLCF_LS(1),CFORM) TRIM(CHID),'_',NM,'_','lsfs.sf'
IF(NMESHES == 1) WRITE(FN_SLCF_LS(1),CFORM) TRIM(CHID),'_','lsfs.sf'
OPEN(LU_SLCF_LS(1),FILE=FN_SLCF_LS(1),FORM='UNFORMATTED',STATUS='REPLACE')
!OPEN(LU_SLCF_LS,FILE=TRIM(CHID)//'_lsfs.sf',FORM='UNFORMATTED',STATUS='REPLACE')
WRITE(LU_SLCF_LS(1)) SMOKEVIEW_LABEL(1:30)
WRITE(LU_SLCF_LS(1)) SMOKEVIEW_LABEL(1:30)
WRITE(LU_SLCF_LS(1)) UNITS(1:30)
WRITE(LU_SLCF_LS(1))0,IBAR,0,JBAR,1,1

IF (NM == 1) THEN !write to smv file
  DO I=1,NMESHES 
    IF (NMESHES == 1) THEN
      WRITE(FN_SLCF_LS(1),CFORM) TRIM(CHID),'_','lsfs.sf'
    ELSE
      WRITE(FN_SLCF_LS(1),CFORM) TRIM(CHID),'_',I,'_','lsfs.sf'
    ENDIF
    WRITE(LU_SMV,'(A,5X,I3,5x,F7.2)') 'SLCT ',I,0.1
    WRITE(LU_SMV,'(A)')FN_SLCF_LS(1)
    WRITE(LU_SMV,'(A)') 'phifield'
    WRITE(LU_SMV,'(A)') 'phifield'
    WRITE(LU_SMV,'(A)') '-'
  ENDDO
ENDIF

!--Time of Arrival for animation via Smokeview
LU_SLCF_LS(2) = GET_FILE_NUMBER()
SMOKEVIEW_LABEL = 'LS TOA'
SMOKEVIEW_BAR_LABEL = 'LS TOA'
UNITS  = 's'
IF(NMESHES  > 1) WRITE(FN_SLCF_LS(2),CFORM) TRIM(CHID),'_',NM,'_','lstoa.sf'
IF(NMESHES == 1) WRITE(FN_SLCF_LS(2),CFORM) TRIM(CHID),'_','lstoa.sf'
OPEN(LU_SLCF_LS(2),FILE=FN_SLCF_LS(2),FORM='UNFORMATTED',STATUS='REPLACE')
!OPEN(LU_SLCF_TOA_LS,FILE=TRIM(CHID)//'_lstoa.sf',FORM='UNFORMATTED',STATUS='REPLACE')
WRITE(LU_SLCF_LS(2)) SMOKEVIEW_LABEL(1:30)
WRITE(LU_SLCF_LS(2)) SMOKEVIEW_LABEL(1:30)
WRITE(LU_SLCF_LS(2)) UNITS(1:30)
WRITE(LU_SLCF_LS(2))0,IBAR,0,JBAR,1,1

IF (NM == 1) THEN !write to smv file
  DO I=1,NMESHES 
    IF (NMESHES == 1) THEN
      WRITE(FN_SLCF_LS(2),CFORM) TRIM(CHID),'_','lstoa.sf'
    ELSE
      WRITE(FN_SLCF_LS(2),CFORM) TRIM(CHID),'_',I,'_','lstoa.sf'
    ENDIF
    WRITE(LU_SMV,'(A,5X,I3,5X,F7.2)') 'SLCT ',I,0.1
    WRITE(LU_SMV,'(A)')FN_SLCF_LS(2)
    WRITE(LU_SMV,'(A)') 'LS TOA'
    WRITE(LU_SMV,'(A)') 'LS TOA'
    WRITE(LU_SMV,'(A)') 's'
  ENDDO
ENDIF

!--ROS magnitude for animation in Smokeview
LU_SLCF_LS(3) = GET_FILE_NUMBER()
SMOKEVIEW_LABEL = 'LS ROS'
SMOKEVIEW_BAR_LABEL = 'LS ROS'
UNITS  = 'm/s'
IF(NMESHES  > 1) WRITE(FN_SLCF_LS(3),CFORM) TRIM(CHID),'_',NM,'_','lsros.sf'
IF(NMESHES == 1) WRITE(FN_SLCF_LS(3),CFORM) TRIM(CHID),'_','lsros.sf'
OPEN(LU_SLCF_LS(3),FILE=FN_SLCF_LS(3),FORM='UNFORMATTED',STATUS='REPLACE')
!OPEN(LU_SLCF_ROS_LS,FILE=TRIM(CHID)//'_lsros.sf',FORM='UNFORMATTED',STATUS='REPLACE')
WRITE(LU_SLCF_LS(3)) SMOKEVIEW_LABEL(1:30)
WRITE(LU_SLCF_LS(3)) SMOKEVIEW_LABEL(1:30)
WRITE(LU_SLCF_LS(3)) UNITS(1:30)
WRITE(LU_SLCF_LS(3))0,IBAR,0,JBAR,1,1

IF (NM == 1) THEN !write to smv file
  DO I=1,NMESHES 
    IF (NMESHES == 1) THEN
      WRITE(FN_SLCF_LS(3),CFORM) TRIM(CHID),'_','lsros.sf'
    ELSE
      WRITE(FN_SLCF_LS(3),CFORM) TRIM(CHID),'_',I,'_','lsros.sf'
    ENDIF
    WRITE(LU_SMV,'(A,5X,I3,5X,F7.2)') 'SLCT ',I,0.1
    WRITE(LU_SMV,'(A)')FN_SLCF_LS(3)
    WRITE(LU_SMV,'(A)') 'LS ROS'
    WRITE(LU_SMV,'(A)') 'LS ROS'
    WRITE(LU_SMV,'(A)') 'm/s'
  ENDDO
ENDIF

!--Fire line intensity at time of fire arrival for animation in Smokeview
LU_SLCF_LS(4) = GET_FILE_NUMBER()
SMOKEVIEW_LABEL = 'LS FLI'
SMOKEVIEW_BAR_LABEL = 'LS FLI'
UNITS  = 'kW/m'
IF(NMESHES  > 1) WRITE(FN_SLCF_LS(4),CFORM) TRIM(CHID),'_',NM,'_','lsfli.sf'
IF(NMESHES == 1) WRITE(FN_SLCF_LS(4),CFORM) TRIM(CHID),'_','lsfli.sf'
OPEN(LU_SLCF_LS(4),FILE=FN_SLCF_LS(4),FORM='UNFORMATTED',STATUS='REPLACE')
!OPEN(LU_SLCF_FLI_LS,FILE=TRIM(CHID)//'_lsfli.sf',FORM='UNFORMATTED',STATUS='REPLACE')
WRITE(LU_SLCF_LS(4)) SMOKEVIEW_LABEL(1:30)
WRITE(LU_SLCF_LS(4)) SMOKEVIEW_LABEL(1:30)
WRITE(LU_SLCF_LS(4)) UNITS(1:30)
WRITE(LU_SLCF_LS(4))0,IBAR,0,JBAR,1,1

IF (NM == 1) THEN !write to smv file
  DO I=1,NMESHES 
    IF (NMESHES == 1) THEN
      WRITE(FN_SLCF_LS(4),CFORM) TRIM(CHID),'_','lsfli.sf'
    ELSE
      WRITE(FN_SLCF_LS(4),CFORM) TRIM(CHID),'_',I,'_','lsfli.sf'
    ENDIF
    WRITE(LU_SMV,'(A,5X,I3,5X,F7.2)') 'SLCT ',I,0.1
    WRITE(LU_SMV,'(A)')FN_SLCF_LS(4)
    WRITE(LU_SMV,'(A)') 'LS FLI'
    WRITE(LU_SMV,'(A)') 'LS FLI'
    WRITE(LU_SMV,'(A)') 'kW/m'
  ENDDO
ENDIF

!--HRRPUA at all output times for animation in Smokeview
IF (VEG_LEVEL_SET_COUPLED) THEN

LU_SLCF_LS(5) = GET_FILE_NUMBER()
SMOKEVIEW_LABEL = 'LS HRRPUA'
SMOKEVIEW_BAR_LABEL = 'LS HRRPUA'
UNITS  = 'kW/m^2'
IF(NMESHES  > 1) WRITE(FN_SLCF_LS(5),CFORM) TRIM(CHID),'_',NM,'_','lshrrpua.sf'
IF(NMESHES == 1) WRITE(FN_SLCF_LS(5),CFORM) TRIM(CHID),'_','lshrrpua.sf'
OPEN(LU_SLCF_LS(5),FILE=FN_SLCF_LS(5),FORM='UNFORMATTED',STATUS='REPLACE')
WRITE(LU_SLCF_LS(5)) SMOKEVIEW_LABEL(1:30)
WRITE(LU_SLCF_LS(5)) SMOKEVIEW_LABEL(1:30)
WRITE(LU_SLCF_LS(5)) UNITS(1:30)
WRITE(LU_SLCF_LS(5))0,IBAR,0,JBAR,1,1

IF (NM == 1) THEN !write to smv file
  DO I=1,NMESHES 
    IF (NMESHES == 1) THEN
      WRITE(FN_SLCF_LS(5),CFORM) TRIM(CHID),'_','lshrrpua.sf'
    ELSE
      WRITE(FN_SLCF_LS(5),CFORM) TRIM(CHID),'_',I,'_','lshrrpua.sf'
    ENDIF
    WRITE(LU_SMV,'(A,5X,I3,5X,F7.2)') 'SLCT ',I,0.1
    WRITE(LU_SMV,'(A)')FN_SLCF_LS(5)
    WRITE(LU_SMV,'(A)') 'LS HRRPUA'
    WRITE(LU_SMV,'(A)') 'LS HRRPUA'
    WRITE(LU_SMV,'(A)') 'kW/m^2'
  ENDDO
ENDIF

ENDIF

!--Crown Fire Probablity (Cruz & Alexander) for animation in Smokeview
LU_SLCF_LS(6) = GET_FILE_NUMBER()
SMOKEVIEW_LABEL = 'LS ProbCrown'
SMOKEVIEW_BAR_LABEL = 'LS ProbCrown'
UNITS  = '-'
IF(NMESHES  > 1) WRITE(FN_SLCF_LS(6),CFORM) TRIM(CHID),'_',NM,'_','lsprobc.sf'
IF(NMESHES == 1) WRITE(FN_SLCF_LS(6),CFORM) TRIM(CHID),'_','lsprobc.sf'
OPEN(LU_SLCF_LS(6),FILE=FN_SLCF_LS(6),FORM='UNFORMATTED',STATUS='REPLACE')
WRITE(LU_SLCF_LS(6)) SMOKEVIEW_LABEL(1:30)
WRITE(LU_SLCF_LS(6)) SMOKEVIEW_LABEL(1:30)
WRITE(LU_SLCF_LS(6)) UNITS(1:30)
WRITE(LU_SLCF_LS(6))0,IBAR,0,JBAR,1,1

IF (NM == 1) THEN !write to smv file
  DO I=1,NMESHES 
    IF (NMESHES == 1) THEN
      WRITE(FN_SLCF_LS(6),CFORM) TRIM(CHID),'_','lsprobc.sf'
    ELSE
      WRITE(FN_SLCF_LS(6),CFORM) TRIM(CHID),'_',I,'_','lsprobc.sf'
    ENDIF
    WRITE(LU_SMV,'(A,5X,I3,5X,F7.2)') 'SLCT ',I,0.1
    WRITE(LU_SMV,'(A)')FN_SLCF_LS(6)
    WRITE(LU_SMV,'(A)') 'LS PROBCROWN'
    WRITE(LU_SMV,'(A)') 'LS PROBCROWN'
    WRITE(LU_SMV,'(A)') '-'
  ENDDO
ENDIF

!******************* Open file contained HRRPUA in order to implement the "burner" method which represents a fireline 
!                    using proxies for the actual HRRPUA that are obtained from remote-sensing, or other, sources. 
! ASSUMES SINGLE MESH FOR NOW

IF (VEG_LEVEL_SET_BURNERS_FOR_FIRELINE) THEN
  LU_SLCF_LS(7) = GET_FILE_NUMBER()
  FN_SLCF_LS(7) = BRNRINFO(BURNER_FILE(NM))%BRNRFILE
  OPEN(LU_SLCF_LS(7),FILE=FN_SLCF_LS(7),FORM='UNFORMATTED',STATUS='OLD')
  READ(LU_SLCF_LS(7)) SMOKEVIEW_LABEL(1:30)
  READ(LU_SLCF_LS(7)) SMOKEVIEW_LABEL(1:30)
  READ(LU_SLCF_LS(7)) UNITS(1:30)
  READ(LU_SLCF_LS(7)) I,I,I,I,I,I
  READ(LU_SLCF_LS(7)) LSET_TIME_HRRPUA_BURNER 
  READ(LU_SLCF_LS(7)) ((HRRPUA_IN(I,J),I=0,IBAR),J=0,JBAR) 
ENDIF

!-- ASCII files of level set quantities

!--Time of arrival binary format
!LU_TOA_LS = GET_FILE_NUMBER()
!OPEN(LU_TOA_LS,FILE='time_of_arrival.txt',FORM='UNFORMATTED',STATUS='REPLACE')
!WRITE(LU_TOA_LS) NX_LS,NY_LS
!WRITE(LU_TOA_LS) XS,XF,YS,YF

!--Time of arrival ASCII format
!LU_TOA_LS = GET_FILE_NUMBER()
!OPEN(LU_TOA_LS,FILE='toa_LS.txt',STATUS='REPLACE')
!WRITE(LU_TOA_LS,'(I5)') NX_LS,NY_LS
!WRITE(LU_TOA_LS,'(F7.2)') XS,XF,YS,YF

!Write across row (TOA(1,1), TOA(1,2), ...) to match Farsite output
!WRITE(LU_TOA_LS,'(F7.2)') ((TOA(IDUM,JDUM),JDUM=1,NY_LS),IDUM=1,NX_LS)
!CLOSE(LU_TOA_LS)

!--Rate of spread ASCII format (
!LU_ROSX_LS = GET_FILE_NUMBER()
!OPEN(LU_ROSX_LS,FILE='rosx_LS.txt',STATUS='REPLACE')
!WRITE(LU_ROSX_LS,'(I5)') NX_LS,NY_LS
!WRITE(LU_ROSX_LS,'(F7.2)') XS,XF,YS,YF
!LU_ROSY_LS = GET_FILE_NUMBER()
!OPEN(LU_ROSY_LS,FILE='rosy_LS.txt',STATUS='REPLACE')
!WRITE(LU_ROSY_LS,'(I5)') NX_LS,NY_LS
!WRITE(LU_ROSY_LS,'(F7.2)') XS,XF,YS,YF

!--Fire line intensity ASCII format
!LU_FLI_LS = GET_FILE_NUMBER()
!OPEN(LU_FLI_LS,FILE='fli_LS.txt',STATUS='REPLACE')
!WRITE(LU_FLI_LS,'(I5)') NX_LS,NY_LS
!WRITE(LU_FLI_LS,'(F7.2)') XS,XF,YS,YF

!--Crown Fire Probability (Cruz & Alexander) ASCII format
!LU_CRWN_PROB_LS = GET_FILE_NUMBER()
!OPEN(LU_CRWN_PROB_LS,FILE='crwn_prob_LS.txt',STATUS='REPLACE')
!WRITE(LU_CRWN_PROB_LS,'(I5)') NX_LS,NY_LS
!WRITE(LU_CRWN_PROB_LS,'(F7.2)') XS,XF,YS,YF

! --------------------------------------------------------------------------------------- 
!-- Allocate and initialize arrays
! --------------------------------------------------------------------------------------- 
!
!ALLOCATE(HEAD_WIDTH(NX_LS,NY_LS))  ; CALL ChkMemErr('VEGE:LEVEL SET','HEAD_WIDTH',IZERO) ; HEAD_WIDTH = 1.0_EB
!ALLOCATE(ROS_HEAD(NX_LS,NY_LS))    ; CALL ChkMemErr('VEGE:LEVEL SET','ROS_HEAD',IZERO) ; ROS_HEAD = 0.0_EB
!ALLOCATE(ROS_FLANK(NX_LS,NY_LS))   ; CALL ChkMemErr('VEGE:LEVEL SET','ROS_FLANK',IZERO) ; ROS_FLANK = 0.0_EB
!ALLOCATE(ROS_BACKU(NX_LS,NY_LS))   ; CALL ChkMemErr('VEGE:LEVEL SET','ROS_BACKU',IZERO) ; ROS_BACKU = 0.0_EB
!ALLOCATE(WIND_EXP(NX_LS,NY_LS))   ; CALL ChkMemErr('VEGE:LEVEL SET','WIND_EXP',IZERO) ; WIND_EXP = 1.0_EB
!ALLOCATE(FLANKFIRE_LIFETIME(NX_LS,NY_LS)) ; CALL ChkMemErr('VEGE:LEVEL SET','FLANKFIRE_LIFETIME',IZERO)
!FLANKFIRE_LIFETIME = 0.0_EB
!ROS_HEAD1 needed when dependence on head width is computed
       
!--Level set values (Phi)

!ALLOCATE(PHI_LS(NX_LS,NY_LS)) ; CALL ChkMemErr('VEGE:LEVEL SET','PHI_LS',IZERO)
!ALLOCATE(PHI0_LS(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','PHI0_LS',IZERO)
!ALLOCATE(PHI1_LS(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','PHI1_LS',IZERO)
!ALLOCATE(PHI_TEMP_LS(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','PHI_TEMP_LS',IZERO)

!--Burntime for heat injection
!ALLOCATE(BURN_TIME_LS(NX_LS,NY_LS)) ; CALL ChkMemErr('VEGE:LEVEL SET','BURN_TIME_LS',IZERO) 
!BURN_TIME_LS = 0._EB

!--Crown fraction burned for FLI and HRRPUA
!ALLOCATE(CFB_LS(NX_LS,NY_LS)) ; CALL ChkMemErr('VEGE:LEVEL SET','CFB_LS',IZERO) 
!CFB_LS = 0._EB

!--Wind speeds at ground cells in domain
!ALLOCATE(U_LS(NX_LS,NY_LS)) ; CALL ChkMemErr('VEGE:LEVEL SET','U_LS',IZERO) ; U_LS = 0._EB
!ALLOCATE(V_LS(NX_LS,NY_LS)) ; CALL ChkMemErr('VEGE:LEVEL SET','V_LS',IZERO) ; V_LS = 0._EB

!ALLOCATE(FLUX0_LS(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','FLUX0_LS',IZERO)
!ALLOCATE(FLUX1_LS(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','FLUX1_LS',IZERO)

!--Slopes (gradients)
!ALLOCATE(DZTDX(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','DZDTX',IZERO)
!ALLOCATE(DZTDY(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','DZDTY',IZERO)
!ALLOCATE(MAG_ZT(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','MAG_ZT',IZERO)

!--Computational grid
ALLOCATE(X_LS(NX_LS))   ; CALL ChkMemErr('VEGE:LEVEL SET','X_LS',IZERO)
ALLOCATE(Y_LS(NY_LS+1)) ; CALL ChkMemErr('VEGE:LEVEL SET','Y_LS',IZERO)

!-- Arrays for outputs
!|ROS|
!ALLOCATE(MAG_SR_OUT(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','MAG_SR_OUT',IZERO)

!ALLOCATE(PHI_OUT(NX_LS,NY_LS)) ; CALL ChkMemErr('VEGE:LEVEL SET','PHI_OUT',IZERO)
!PHI_OUT = 0.0

!--Time of arrival, s
!ALLOCATE(TOA(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','TOA',IZERO)
!TOA = -1.0

!--Rate of spread components, m/s
!ALLOCATE(ROS_X_OUT(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','ROS_X_OUT',IZERO)
!ROS_X_OUT =  0.0
!ALLOCATE(ROS_Y_OUT(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','ROS_Y_OUT',IZERO)
!ROS_Y_OUT =  0.0

!--Fire line intensity, kW/m 
!ALLOCATE(FLI_OUT(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','FLI_OUT',IZERO)
!FLI_OUT = -1.0

!--Cruz Crown fire probablity 
!ALLOCATE(CRUZ_CROWN_PROB(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','CRUZ_CROWN_PROB',IZERO)
!CRUZ_CROWN_PROB = 0.0_EB
!ALLOCATE(CRUZ_CROWN_PROB_OUT(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','CRUZ_CROWN_PROB_OUT',IZERO)
!CRUZ_CROWN_PROB_OUT = 0.0


!----------Rothermel 'Phi' factors for effects of Wind and Slope on ROS ----------
!--Not to be confused with the level set value (Phi)
!ALLOCATE(PHI_WS(NX_LS,NY_LS))  ; CALL ChkMemErr('VEGE:LEVEL SET','PHI_W',IZERO)   ; PHI_WS  = 0.0_EB
!ALLOCATE(PHI_S(NX_LS,NY_LS))   ; CALL ChkMemErr('VEGE:LEVEL SET','PHI_S',IZERO)   ; PHI_S   = 0.0_EB
!ALLOCATE(PHI_S_X(NX_LS,NY_LS)) ; CALL ChkMemErr('VEGE:LEVEL SET','PHI_S_X',IZERO) ; PHI_S_X = 0.0_EB
!ALLOCATE(PHI_S_Y(NX_LS,NY_LS)) ; CALL ChkMemErr('VEGE:LEVEL SET','PHI_S_Y',IZERO) ; PHI_S_Y = 0.0_EB
!ALLOCATE(PHI_W_X(NX_LS,NY_LS)) ; CALL ChkMemErr('VEGE:LEVEL SET','PHI_W_X',IZERO) ; PHI_W_X = 0.0_EB
!ALLOCATE(PHI_W_Y(NX_LS,NY_LS)) ; CALL ChkMemErr('VEGE:LEVEL SET','PHI_W_Y',IZERO) ; PHI_W_Y = 0.0_EB
!ALLOCATE(PHI_W(NX_LS,NY_LS))   ; CALL ChkMemErr('VEGE:LEVEL SET','PHI_W',IZERO)   ; PHI_W   = 0.0_EB
!ALLOCATE(UMF_X(NX_LS,NY_LS))   ; CALL ChkMemErr('VEGE:LEVEL SET','UMF_X',IZERO)   ; UMF_X   = 0.0_EB
!ALLOCATE(UMF_Y(NX_LS,NY_LS))   ; CALL ChkMemErr('VEGE:LEVEL SET','UMF_Y',IZERO)   ; UMF_Y   = 0.0_EB

!--UMF = wind speed at mean flame heights
!ALLOCATE(UMF(NX_LS,NY_LS))    ; CALL ChkMemErr('VEGE:LEVEL SET','UMF',IZERO) ; UMF = 0.0_EB
!ALLOCATE(THETA_ELPS(NX_LS,NY_LS)) ; CALL ChkMemErr('VEGE:LEVEL SET','THETA_ELPS',IZERO)
!THETA_ELPS   = 0.0_EB ! Normal to fireline
!-----------------------------------------------------------------------------------  

    
!--ROS in X and Y directions for elliptical model
!ALLOCATE(SR_X_LS(NX_LS,NY_LS)) ; CALL ChkMemErr('VEGE:LEVEL SET','SR_X_LS',IZERO) ; SR_X_LS =0.0_EB
!ALLOCATE(SR_Y_LS(NX_LS,NY_LS)) ; CALL ChkMemErr('VEGE:LEVEL SET','SR_Y_LS',IZERO) ; SR_Y_LS =0.0_EB
    
!--Aspect of terrain slope for elliptical model (currently not used)
ALLOCATE(ASPECT(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','ASPECT',IZERO) ; ASPECT = 0.0_EB
    
!Location of computation grid-cell faces
DO I = 0,NX_LS-1
!X_LS(I+1) = -0.5_EB*LX + 0.5_EB*DX_LS + DX_LS*REAL(I,EB)
!X_LS(I+1) = XS + 0.5_EB*DX_LS + DX_LS*REAL(I,EB)
 X_LS(I+1) = XS + 0.5_EB*DX(I) + DX(I)*REAL(I,EB)
ENDDO
!
DO J = 0,NY_LS
 Y_LS(J+1) = YS + DY_LS*REAL(J,EB)
ENDDO

!Compute components of terrain slope gradient and magnitude of gradient

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

   DZTDX(I,J) = (G_EAST-G_WEST) * RDX(I) !IDX_LS
   DZTDY(I,J) = (G_NORTH-G_SOUTH) * RDY(J) !IDY_LS


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

PHI_MIN_LS = -1._EB
PHI_MAX_LS =  1._EB
!PHI_LS     = PHI_MIN_LS
!PHI_LS     = LSET_PHI(1:NX_LS,1:NY_LS,1)
PHI_LS     = LSET_PHI(0:IBP1,0:JBP1)

LSET_INIT_WALL_CELL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
  WC  => WALL(IW)
  IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE LSET_INIT_WALL_CELL_LOOP
  SF  => SURFACE(WC%SURF_INDEX)
  SF%VEG_LSET_SURF_HEIGHT = MAX(0.001_EB,SF%VEG_LSET_SURF_HEIGHT)
! WC%VEG_HEIGHT = SF%VEG_LSET_SURF_HEIGHT
  WC%VEG_HEIGHT = 0.0_EB

  IIG = WC%ONE_D%IIG
  JJG = WC%ONE_D%JJG
  KKG = WC%ONE_D%KKG

! Ignite landscape at user specified location if ignition is at time zero
  IF (SF%VEG_LSET_IGNITE_TIME == 0.0_EB) THEN 
!print '(A,ES12.4,1x,3I)','veg:LS lset_ignite_time,iig,jjg ',sf%veg_lset_ignite_time,iig,jjg
    PHI_LS(IIG,JJG) = PHI_MAX_LS 
    BURN_TIME_LS(IIG,JJG) = 99999.0_EB
!   IF (SF%HRRPUA > 0.0_EB) THEN !mimic burner
!     WC%VEG_LSET_SURFACE_HEATFLUX = -SF%HRRPUA
!     IF (VEG_LEVEL_SET_SURFACE_HEATFLUX) WC%ONE_D%QCONF = WC%VEG_LSET_SURFACE_HEATFLUX
!     WC%LSET_FIRE = .TRUE.
!   ENDIF
  ENDIF

! Wind field 
  U_LS(IIG,JJG) = U(IIG,JJG,KKG)
  V_LS(IIG,JJG) = V(IIG,JJG,KKG)
!WRITE(LU_OUTPUT,*)'veg: u,v_ls(i,j)',u_ls(iig,jjg),v_ls(iig,jjg)

  IF (.NOT. SF%VEG_LSET_SPREAD) CYCLE LSET_INIT_WALL_CELL_LOOP
!WRITE(LU_OUTPUT,*)'x,y,z and U,V',X(IIG),Y(JJG),Z(KKG),U(IIG,JJG,KKG),V(IIG,JJG,KKG)

  !Diagnostics
  !WRITE(LU_OUTPUT,*)'IIG,JJG',iig,jjg
  !WRITE(LU_OUTPUT,*)'ROS_HEAD',SF%VEG_LSET_ROS_HEAD
  !WRITE(LU_OUTPUT,*)'ROS_HEAD,ROS_FLANK,ROS_BACK',SF%VEG_LSET_ROS_HEAD,SF%VEG_LSET_ROS_FLANK,SF%VEG_LSET_ROS_BACK


  UMAG     = SQRT(U_LS(IIG,JJG)**2 + V_LS(IIG,JJG)**2)

  HEAD_WIDTH(IIG,JJG)= DX(1) !DX_LS
  ROS_HEAD(IIG,JJG)  = SF%VEG_LSET_ROS_HEAD
  ROS_FLANK(IIG,JJG) = SF%VEG_LSET_ROS_FLANK
  ROS_BACKU(IIG,JJG) = SF%VEG_LSET_ROS_BACK
  WIND_EXP(IIG,JJG)  = SF%VEG_LSET_WIND_EXP
!print '(A,1x,L2,4ES12.4)','vegeLS: LSET_SPREAD,ROS:HEAD,FLANK,BACK,WIND_EXP',sf%veg_lset_spread,sf%veg_lset_ros_head, &
!                           sf%veg_lset_ros_flank,sf%veg_lset_ros_back,sf%veg_lset_wind_exp
  
!If any surfaces uses tan^2 function for slope, tan^2 will be used throughout simulation
  IF (SF%VEG_LSET_TAN2) LSET_TAN2=.TRUE.
  
!Compute head ROS and coefficients when using the assumption of ellipsed shaped fireline
  IF_ELLIPSE:IF (SF%VEG_LSET_ELLIPSE) THEN    
    
!-- If any surfaces set to ellipse, then elliptical model used for all surfaces 
    IF (.NOT. LSET_ELLIPSE) LSET_ELLIPSE=.TRUE.

!-- AU grassland fire, head ROS for infinite head width
    IF (SF%VEG_LSET_SURFACE_FIRE_HEAD_ROS_MODEL=='AU GRASS') CALL AUGRASS_HEADROS(NM,IIG,JJG,KKG,SF%VEG_LSET_SURF_EFFM)
!     UMAG = SQRT(U_LS(IIG,JJG)**2 + V_LS(IIG,JJG)**2)
!     ROS_HEAD(IIG,JJG)  = (0.165_EB + 0.534_EB*UMAG)*EXP(-0.108*SF%VEG_LSET_SURF_EFFM)
!   ENDIF

!-- Rothermel surface veg head fire ROS model
    IF (SF%VEG_LSET_ELLIPSE .AND. SF%VEG_LSET_SURFACE_FIRE_HEAD_ROS_MODEL=='ROTHERMEL') THEN    
!---- Slope factors
      CALL ROTH_SLOPE_COEFF(NM,IIG,JJG,SF%VEG_LSET_BETA)
!---- Wind, combined wind & slope, and midflame windspeed factors
       CALL ROTH_WINDANDSLOPE_COEFF_HEADROS(NM,IIG,JJG,KKG,SF%VEG_LSET_BETA,SF%VEG_LSET_SURF_HEIGHT,            &
          SF%VEG_LSET_CANOPY_HEIGHT,SF%VEG_LSET_SIGMA,SF%VEG_LSET_ROTH_ZEROWINDSLOPE_ROS, &
          SF%VEG_LSET_WAF_UNSHELTERED,SF%VEG_LSET_WAF_SHELTERED)
    ENDIF

!-- Cruz et al. crown fire head fire ROS model (needed to determine time step based on surface and crown fire head ROS)
    IF (SF%VEG_LSET_CROWN_FIRE_HEAD_ROS_MODEL=='CRUZ' .AND. VEG_LEVEL_SET_UNCOUPLED) THEN    
      VERT_CANOPY_EXTENT = SF%VEG_LSET_CANOPY_HEIGHT - SF%VEG_LSET_SURF_HEIGHT - SF%VEG_LSET_FUEL_STRATA_GAP
      CALL CRUZ_CROWN_FIRE_HEADROS(NM,IIG,JJG,KKG,SF%VEG_LSET_CANOPY_BULK_DENSITY,SF%VEG_LSET_SURF_EFFM,            &
           SF%VEG_LSET_FUEL_STRATA_GAP,SF%VEG_LSET_SURF_LOAD,SF%VEG_LSET_CRUZ_PROB_PASSIVE,                         &
           SF%VEG_LSET_CRUZ_PROB_ACTIVE,SF%VEG_LSET_CRUZ_PROB_CROWN,SF%VEG_LSET_SURFACE_FIRE_HEAD_ROS_MODEL,        &
           VERT_CANOPY_EXTENT,SF%VEG_LSET_CANOPY_HEIGHT)
    ENDIF
  
  ENDIF IF_ELLIPSE

ENDDO LSET_INIT_WALL_CELL_LOOP

UMAX_LS  = MAXVAL(ABS(U_LS))
VMAX_LS  = MAXVAL(ABS(V_LS))
UMAG     = SQRT(UMAX_LS**2 + VMAX_LS**2)

!WRITE(LU_OUTPUT,*)'before assign ROS'
!WRITE(LU_OUTPUT,*)'ROS_HEAD max',MAXVAL(ROS_HEAD)
!ROS_HEAD1 = MAXVAL(ROS_HEAD)
!WRITE(LU_OUTPUT,*)'ROS_HEAD1',ROS_HEAD1

SR_MAX   = MAXVAL(ROS_HEAD)
SR_MAX   = MAX(SR_MAX,MAXVAL(ROS_FLANK))
DYN_SR_MAX_LS = 0._EB

! Write diagnostic to standard output
!WRITE(LU_OUTPUT,*)'ROS_HEAD max',MAXVAL(ROS_HEAD)
!ROS_HEAD1 = MAXVAL(ROS_HEAD)
!WRITE(LU_OUTPUT,*)'ROS_HEAD1',ROS_HEAD1

IF (LSET_ELLIPSE) THEN
    WRITE(LU_OUTPUT,*)'Mesh number',NM
    WRITE(LU_OUTPUT,*)'Phi_S max',MAXVAL(PHI_S)
    WRITE(LU_OUTPUT,*)'Phi_W max',MAXVAL(PHI_W)
    WRITE(LU_OUTPUT,*)'UMF max',MAXVAL(UMF)
    WRITE(LU_OUTPUT,*)'Mag_zt max',MAXVAL(MAG_ZT)
    WRITE(LU_OUTPUT,*)'SR_MAX',SR_MAX
ENDIF

IF (.NOT. LSET_ELLIPSE) SR_MAX   = 2._EB*SR_MAX !rough accounting for upslope spread aligned with wind

IF (VEG_LEVEL_SET_UNCOUPLED) THEN
 DT_LS = 0.5_EB*MIN(DX(1),DY(1))/SR_MAX
!DT_LS = MESHES(NM)%DT
 DT = DT_LS
!DT_NEXT = DT_LS
ENDIF


LSET_PHI(0:IBP1,0:JBP1) = PHI_LS

!DT_LS = 0.1603_EB !to make AU F19 ignition sequence work

WRITE(LU_OUTPUT,1113)nm,dt_ls
1113 format('vegelsini nm, dt_ls ',1(i2),2x,1(ES12.4))
!WRITE(LU_OUTPUT,*)'flux limiter= ',LIMITER_LS

END SUBROUTINE INITIALIZE_LEVEL_SET_FIREFRONT

!************************************************************************************************
SUBROUTINE ROTH_SLOPE_COEFF(NM,I,J,VEG_BETA)
!************************************************************************************************
!
! Compute components and magnitude of slope coefficient vector that 
! are used in the Rothermel spread rate formula. These, along with the zero wind and zero slope
! Rothermel ROS (given in the input file) and wind coefficient vector (computed below) 
! are used to obtain the local surface fire spread rate
!
INTEGER, INTENT(IN) :: I,J,NM
REAL(EB), INTENT(IN) :: VEG_BETA
REAL(EB) :: DZT_DUM
!REAL(EB) :: DZT_MAG2

CALL POINT_TO_MESH(NM)

!Limit effect to slope lte 80 degrees
!Phi_s_x,y are slope factors
!DZT_DUM = MIN(5.67_EB,ABS(DZTDX(I,J))) ! 5.67 ~ tan 80 deg, used in LS paper, tests show equiv to 60 deg max
DZT_DUM = MIN(1.73_EB,ABS(DZTDX(I,J))) ! 1.73 ~ tan 60 deg
PHI_S_X(I,J) = 5.275_EB * ((VEG_BETA)**(-0.3_EB)) * DZT_DUM**2
PHI_S_X(I,J) = SIGN(PHI_S_X(I,J),DZTDX(I,J))

DZT_DUM = MIN(1.73_EB,ABS(DZTDY(I,J))) ! 1.73 ~ tan 60 deg, used in LS paper
PHI_S_Y(I,J) = 5.275_EB * ((VEG_BETA)**(-0.3_EB)) * DZT_DUM**2
PHI_S_Y(I,J) = SIGN(PHI_S_Y(I,J),DZTDY(I,J))

PHI_S(I,J) = SQRT(PHI_S_X(I,J)**2 + PHI_S_Y(I,J)**2) !used in LS paper
!DZT_MAG2 = DZTDX(I,J)**2 + DZTDY(I,J)**2
!PHI_S(I,J) = 5.275_EB * ((VEG_BETA)**(-0.3_EB)) * DZT_MAG2

END SUBROUTINE ROTH_SLOPE_COEFF

!************************************************************************************************
SUBROUTINE ROTH_WINDANDSLOPE_COEFF_HEADROS(NM,I,J,K,VEG_BETA,SURF_VEG_HT,CANOPY_VEG_HT,VEG_SIGMA,ZEROWINDSLOPE_ROS, &
                                           WAF_UNSHELTERED,WAF_SHELTERED)
!************************************************************************************************
!
! Compute components and magnitude of the wind coefficient vector and the combined
! wind and slope coefficient vectors and use them along with the user defined zero wind, zero slope
! Rothermel (or Behave) surface fire ROS in the Rothermel surface fire spread rate
! formula to obtain the magnitude of the local surface head fire ROS. Top of vegetation is assumed
! to be at the bottom of the computational doamin.

INTEGER,  INTENT(IN) :: I,J,K,NM
REAL(EB), INTENT(IN) :: CANOPY_VEG_HT,SURF_VEG_HT,VEG_BETA,VEG_SIGMA,WAF_UNSHELTERED,WAF_SHELTERED,  &
                        ZEROWINDSLOPE_ROS
LOGICAL :: UNIFORM_UV
INTEGER :: KDUM,KWIND
REAL(EB) :: CONSFCTR,FCTR1,FCTR2,PHX,PHY,MAG_PHI,U6PH,VEG_HT,V6PH,WAF_6M,WAF_MID,Z6PH,ZWFDS

CALL POINT_TO_MESH(NM)

!print*,'n_csvf',n_csvf
!print*,'crown_veg',crown_veg
!print*,'k,z(k-1),z(k)',k,z(k-1),z(k)

VEG_HT = MAX(SURF_VEG_HT,CANOPY_VEG_HT)
FCTR1 = 0.64_EB*VEG_HT !constant in logrithmic wind profile Albini & Baughman INT-221 1979 or 
!                       !Andrews RMRS-GTR-266 2012 (p. 8, Eq. 4)
FCTR2 = 1.0_EB/(0.13_EB*VEG_HT) !constant in log wind profile
Z6PH  = 6.1_EB + VEG_HT
ZWFDS = ZC(K) - Z(K-1) !Height of velocity in first cell above veg, ZC(K)=cell center, Z(K-1)=height of K cell bottom
UNIFORM_UV = .FALSE.
!
!Find the wind components at 6.1 m above the vegetation for the case of a uniform
!wind field (i.e., equivalent to conventional FARSITE). This was used in 2015 LS & 
!FARSITE paper. N_CSVF = 0 when no initial wind field has been read in from a file.
IF (N_CSVF == 0 .AND. VEG_LEVEL_SET_UNCOUPLED) THEN
  U6PH = U_LS(I,J) 
  V6PH = V_LS(I,J)
  UNIFORM_UV = .TRUE.
ENDIF
!
!Find the wind components at 6.1 m above the vegetation for the case of nonuniform wind field. 
!The wind field can be predefined and read in at code initialization or the level set computation 
!is coupled to the CFD computation
!
!---U,V at 6.1 above the veg height computed from the WAF when vegetation height + 6.1 m is above or 
!   equal to the first u,v location on grid
!print*,'zwfds,z6ph,uniform_uv',zwfds,z6ph,uniform_uv
KWIND = 0
KDUM = 0
IF (ZWFDS <= 6.1_EB .AND. .NOT. UNIFORM_UV) THEN 
!Find k array index for first grid cell that has ZC > 6.1 m 
   KWIND = 0
   KDUM  = K
   DO WHILE (ZWFDS <= 6.1_EB) !this assumes the bottom computational boundary = top of veg
      KWIND = KDUM
      KDUM  = KDUM + 1
      ZWFDS = ZC(KDUM) - Z(K-1)
   ENDDO
   ZWFDS = ZC(KWIND) - Z(K-1)
   WAF_6M = LOG((Z6PH-FCTR1)*FCTR2)/LOG((ZWFDS+VEG_HT-FCTR1)*FCTR2) !wind adjustment factor from log wind profile
   U6PH  = WAF_6M*U(I,J,KWIND)
   V6PH  = WAF_6M*V(I,J,KWIND)
ENDIF
!
!---U,V at 6.1 m above the veg height computed from the WAF when vegetation height + 6.1 m is below
!   first u,v location on grid
IF (ZWFDS > 6.1_EB .AND. .NOT. UNIFORM_UV) THEN 
   WAF_6M = LOG((Z6PH-FCTR1)*FCTR2)/LOG((ZWFDS+VEG_HT-FCTR1)*FCTR2) 
   U6PH  = WAF_6M*U(I,J,K)
   V6PH  = WAF_6M*V(I,J,K)
ENDIF
!
!Obtain mid-flame wind adjustment factor
!Log profile based wind adjustiment for unsheltered or sheltered condtions are from 
!Andrews 2012, USDA FS Gen Tech Rep. RMRS-GTR-266 (with added SI conversion)
!When using Andrews log formula for sheltered wind the crown fill portion, f, is 0.2
IF (CANOPY_VEG_HT == 0.0_EB) THEN
  WAF_MID = WAF_UNSHELTERED !WAF is from input file
  IF (WAF_UNSHELTERED == -99.0_EB) &
      WAF_MID=1.83_EB/LOG((20.0_EB + 1.18_EB*VEG_HT)/(0.43_EB*VEG_HT))!used in LS vs FS paper
!if(x(i)==21 .and. y(j)==2) print '(A,2x,L2,1ES12.4)','----crown_veg, waf_mid =',crown_veg,waf_mid
ELSE
  WAF_MID = WAF_SHELTERED !WAF is from input file
  IF (WAF_SHELTERED == -99.0_EB)   &
      WAF_MID=0.555_EB/(SQRT(0.20_EB*3.28_EB*VEG_HT)*LOG((20.0_EB + 1.18_EB*VEG_HT)/(0.43_EB*VEG_HT)))
!if(x(i)==21 .and. y(j)==2) print '(A,2x,L2,1ES12.4)','++++crown_veg, waf_mid =',crown_veg,waf_mid
ENDIF

!if (i==41 .and. j==41) then
!print 1116,k,kwind,zwfds,zc(kwind),u(i,j,k),u(i,j,kwind),waf_6m,waf_mid
!endif
!1116 format('(vege,rothwind)',1x,2(I3),1x,6(e15.5))
!
!!Factor 60 converts U from m/s to m/min which is used in the Rothermel model.  
UMF_X(I,J) = WAF_MID * U6PH * 60.0_EB
UMF_Y(I,J) = WAF_MID * V6PH * 60.0_EB
  
!Variables used in Phi_W formulas below (Rothermel model)
B_ROTH = 0.15988_EB * (VEG_SIGMA**0.54_EB)
C_ROTH = 7.47_EB * EXP(-0.8711_EB * (VEG_SIGMA**0.55_EB))
E_ROTH = 0.715_EB * EXP(-0.01094_EB * VEG_SIGMA)
BETA_OP_ROTH = 0.20395_EB * (VEG_SIGMA**(-0.8189_EB))! Optimum packing ratio
     
! Find components of wind factor PHI_W_X, and PHI_W_Y
CONSFCTR = C_ROTH * (3.281_EB**B_ROTH) * (VEG_BETA / BETA_OP_ROTH)**(-E_ROTH)

PHI_W_X(I,J) = CONSFCTR*(ABS(UMF_X(I,J)))**B_ROTH
PHI_W_X(I,J) = SIGN(PHI_W_X(I,J),UMF_X(I,J))

PHI_W_Y(I,J) = CONSFCTR*(ABS(UMF_Y(I,J)))**B_ROTH
PHI_W_Y(I,J) = SIGN(PHI_W_Y(I,J),UMF_Y(I,J))

PHI_W(I,J) =  SQRT(PHI_W_X(I,J)**2 + PHI_W_Y(I,J)**2) 
     

! Find combined wind and slope factor PHI_WS and effective midflame windspeed UMF

IF (PHI_S(I,J) > 0.0_EB) THEN      
        
  PHX = PHI_W_X(I,J) + PHI_S_X(I,J)
  PHY = PHI_W_Y(I,J) + PHI_S_Y(I,J)
  MAG_PHI = SQRT(PHX**2 + PHY**2)
        
!Magnitude of total phi (phi_w + phi_s) for use in spread rate section
  PHI_WS(I,J) = MAG_PHI
        
!Theta_elps, after adjustment below, is angle of direction (0 to 2pi) of highest spread rate
!0<=theta_elps<=2pi as measured clockwise from Y-axis. ATAN2(y,x) is the angle, measured in the
!counterclockwise direction, between the positive x-axis and the line through (0,0) and (x,y)
!positive x-axis  
  THETA_ELPS(I,J) = ATAN2(PHY,PHX)
        
!"Effective midflame windspeed" used in length-to-breadth ratio calculation (spread rate routine)
! is the wind + slope effect obtained by solving Phi_w eqs. above for UMF
! 8/8/13 - Changed phi_ws to Phi_s below to match Farsite, i.e., instead of adding phi_w and phi_s
! and then calculating effective wind speed, phi_s is converted to an effected wind speed and added
! to UMF calculated from the wind. Effective U has units of m/min in Wilson formula.
! 0.3048 ~= 1/3.281
!if phi_s < 0 then a complex value (NaN) results. Using abs(phi_s) and sign function to correct.
        
  UMF_TMP = (((ABS(PHI_S_X(I,J)) * (VEG_BETA / BETA_OP_ROTH)**E_ROTH)/C_ROTH)**(1/B_ROTH))*0.3048
  UMF_TMP = SIGN(UMF_TMP,PHI_S_X(I,J)) 
  UMF_X(I,J) = UMF_X(I,J) + UMF_TMP
        
  UMF_TMP = (((ABS(PHI_S_Y(I,J)) * (VEG_BETA / BETA_OP_ROTH)**E_ROTH)/C_ROTH)**(1/B_ROTH))*0.3048
  UMF_TMP = SIGN(UMF_TMP,PHI_S_Y(I,J))
  UMF_Y(I,J) = UMF_Y(I,J) + UMF_TMP

ELSE !zero slope case
     
  PHI_WS(I,J) = SQRT (PHI_W_X(I,J)**2 + PHI_W_Y(I,J)**2)
  !IF (PHY == 0._EB) PHY = 1.E-6_EB
  !0<= Theta_elps <=2pi as measured clockwise from Y-axis 
  THETA_ELPS(I,J) = ATAN2(PHI_W_Y(I,J),PHI_W_X(I,J))    

ENDIF
    
UMF(I,J) = SQRT(UMF_X(I,J)**2 + UMF_Y(I,J)**2) !used in LS vs FS paper
!UMF(I,J) = UMF(I,J) + (((PHI_S(I,J) * (VEG_BETA / BETA_OP_ROTH)**E_ROTH)/C_ROTH)**(1/B_ROTH))*0.3048
       
!The following two lines convert ATAN2 output to compass system (0 to 2 pi CW from +Y-axis)
THETA_ELPS(I,J) = PIO2 - THETA_ELPS(I,J)
IF (THETA_ELPS(I,J) < 0.0_EB) THETA_ELPS(I,J) = 2.0_EB*PI + THETA_ELPS(I,J)

ROS_HEAD(I,J) = ZEROWINDSLOPE_ROS*(1.0_EB + PHI_WS(I,J)) !used in LS vs FS paper

!if (i==41 .and. j==41) then
!print 1117,ros_head(i,j)
!print*,'-------------------------'
!endif
!1117 format('(vege,rothROS)',1x,1(e15.5))

END SUBROUTINE ROTH_WINDANDSLOPE_COEFF_HEADROS
!
!************************************************************************************************
SUBROUTINE AUGRASS_HEADROS(NM,I,J,K,VEG_MOIST)
!************************************************************************************************
!
! Compute the magnitude of the head fire from an empirical AU grass fire formula

INTEGER,  INTENT(IN) :: I,J,K,NM
REAL(EB), INTENT(IN) :: VEG_MOIST
INTEGER :: KDUM,KWIND
REAL(EB) :: U2MH,UMAG,UMF_TMP,V2MH,VEG_HT
LOGICAL  :: UNIFORM_UV

CALL POINT_TO_MESH(NM)

UNIFORM_UV = .FALSE.

!Find the wind components at 2 m above the ground for the case of a uniform
!wind field (i.e., equivalent to conventional FARSITE). 
!N_CSVF = 0 when no initial wind field has been read in from a file.
IF (N_CSVF == 0 .AND. VEG_LEVEL_SET_UNCOUPLED) THEN
  U2MH = U_LS(I,J) 
  V2MH = V_LS(I,J)
  UNIFORM_UV = .TRUE.
ENDIF

IF (UNIFORM_UV) THEN
!Find k array index corresponding to ~6.1 m AGL height to determine midflame wind
  KWIND = 0
  KDUM  = K
  DO WHILE (ZC(KDUM)-ZC(K) <= 6.1_EB)
    KWIND = KDUM
    KDUM  = KDUM + 1
  ENDDO
  IF (ZC(KBAR) < 6.1_EB) KWIND=1

!Factor to obtain the wind at midflame height (UMF) based on the wind at 6.1 m AGL.  
!From Andrews 2012, USDA FS Gen Tech Rep. RMRS-GTR-266 (with added SI conversion)
  UMF_TMP = 1.83_EB / LOG((20.0_EB + 1.18_EB * VEG_HT) /(0.43_EB * VEG_HT))
  UMF_X(I,J) = UMF_TMP * U(I,J,KWIND)
  UMF_Y(I,J) = UMF_TMP * V(I,J,KWIND)
ENDIF

!Theta_elps, after adjustment below, is angle of direction (0 to 2pi) of highest spread rate
!0<=theta_elps<=2pi as measured clockwise from Y-axis. ATAN2(y,x) is the angle, measured in the
!counterclockwise direction, between the positive x-axis and the line through (0,0) and (x,y)
!positive x-axis  

!Note, unlike the Rothermel ROS case, the slope is assumed to be zero at this point.
!THETA_ELPS(I,J) = ATAN2(UMF_Y(I,J),UMF_X(I,J))
THETA_ELPS(I,J) = ATAN2(V2MH,U2MH)
        
!The following two lines convert ATAN2 output to compass system (0 to 2 pi CW from +Y-axis)
THETA_ELPS(I,J) = PIO2 - THETA_ELPS(I,J)
IF (THETA_ELPS(I,J) < 0.0_EB) THETA_ELPS(I,J) = 2.0_EB*PI + THETA_ELPS(I,J)

!AU grassland head ROS for infinite head width; See Mell et al. "A physics-based approach to 
!modeling grassland fires" Intnl. J. Wildland Fire, 16:1-22 (2007)
!UMAG = SQRT(UMF_X(I,J)**2 + UMF_Y(I,J)**2)
UMAG = SQRT(U2MH**2 + V2MH**2)
ROS_HEAD(I,J)  = (0.165_EB + 0.534_EB*UMAG)*EXP(-0.108*VEG_MOIST)

END SUBROUTINE AUGRASS_HEADROS
!
!************************************************************************************************
SUBROUTINE CRUZ_CROWN_FIRE_HEADROS(NM,I,J,K,CBD,EFFM,FSG,SFC,PROB_PASSIVE,PROB_ACTIVE,PROB_CROWN,SURFACE_FIRE_HEAD_ROS_MODEL, &
                                   VERT_CANOPY_EXTENT,CANOPY_HEIGHT)
!************************************************************************************************
!
! Compute the magnitude of the head fire from an empirical formula. See
!(1) Cruz et al. "Modeling the likelihood of crown fire occurrence in conifer forest stands," 
!    Forest Science, 50: 640-657 (2004)
!(2) Cruz et al. "Development and testing of models for predicting crown fire rate of spread
!    in conifer forest stands," 35:1626-1639 (2005)
!
! CBD = Canopy Bulk Density (kg/m^3)
! EFFM = Effective Fine Fuel Moisture (%), moisture content of fine, dead-down surface vegetation
! FSG = Fuel Strata Gap (m), distance from top of surface fuel to lower limit of the raised fuel 
!       (ladder and canopy) that sustain fire propatation
! SFC = total Surface Fuel Consumption (kg/m^2), sum of forest floor and dead-down roundwood fuel 
!       consumed, surrogate for the amount of vegetation consumed during flaming combustion
!

INTEGER,  INTENT(IN) :: I,J,K,NM
REAL(EB), INTENT(IN) :: CBD,EFFM,FSG,SFC,VERT_CANOPY_EXTENT,CANOPY_HEIGHT
CHARACTER(25), INTENT(IN) :: SURFACE_FIRE_HEAD_ROS_MODEL
LOGICAL :: UNIFORM_UV
INTEGER :: KDUM,KWIND
REAL(EB) :: CAC,CROSA,CROSP,CRLOAD,MPM_TO_MPS,MPS_TO_KPH,EXPG,G,FCTR1,FCTR2,GMAX,PROB,PROB_PASSIVE,PROB_ACTIVE, &
            PROB_CROWN,U10PH,V10PH,UMAG,VEG_HT,WAF_10m,Z10PH,ZWFDS

CALL POINT_TO_MESH(NM)

MPM_TO_MPS = 1._EB/60._EB
MPS_TO_KPH = 3600._EB/1000._EB

VEG_HT = CANOPY_HEIGHT
FCTR1 = 0.64_EB*VEG_HT !constant in logrithmic wind profile Albini & Baughman INT-221 1979 or 
!                       !Andrews RMRS-GTR-266 2012 (p. 8, Eq. 4)
FCTR2 = 1.0_EB/(0.13_EB*VEG_HT) !constant in log wind profile
Z10PH  = 10.0_EB + VEG_HT
ZWFDS = ZC(K) - Z(K-1) !Height of velocity in first cell above terrain, ZC(K)=cell center, Z(K-1)=height of K cell bottom
UNIFORM_UV = .FALSE.
!
!Find the wind components at 10 m above the vegetation for the case of a uniform
!wind field (i.e., equivalent to conventional FARSITE). 
IF (N_CSVF == 0 .AND. VEG_LEVEL_SET_UNCOUPLED) THEN
  U10PH = U_LS(I,J) 
  V10PH = V_LS(I,J)
  UNIFORM_UV = .TRUE.
ENDIF
!
!Find the wind components at 10 m above the vegetation for the case of nonuniform wind field. 
!The wind field can be predefined and read in at code initialization or the level set computation 
!is coupled to the CFD computation
!
!---U,V at 10 m above the veg height computed from the WAF when vegetation height + 10 m is above or 
!   equal to the first u,v location on grid
!print*,'zwfds,z6ph,uniform_uv',zwfds,z6ph,uniform_uv
KWIND = 0
KDUM  = 0
IF (ZWFDS <= 10.0_EB .AND. .NOT. UNIFORM_UV) THEN 
!Find k array index for first grid cell that has ZC > 10 m 
   KWIND = 0
   KDUM  = K
   DO WHILE (ZWFDS <= 10.0_EB) !this assumes the bottom computational boundary = top of veg
      KWIND = KDUM
      KDUM  = KDUM + 1
      ZWFDS = ZC(KDUM) - Z(K-1)
   ENDDO
   ZWFDS = ZC(KWIND) - Z(K-1)
   WAF_10M = LOG((Z10PH-FCTR1)*FCTR2)/LOG((ZWFDS+VEG_HT-FCTR1)*FCTR2) !wind adjustment factor from log wind profile
   U10PH  = WAF_10M*U(I,J,KWIND)
   V10PH  = WAF_10M*V(I,J,KWIND)
ENDIF
!
!---U,V at 10 m above the veg height computed from the WAF when vegetation height + 10 m is below the 
!   first u,v location on grid
IF (ZWFDS > 10.0_EB .AND. .NOT. UNIFORM_UV) THEN 
   WAF_10M = LOG((Z10PH-FCTR1)*FCTR2)/LOG((ZWFDS+VEG_HT-FCTR1)*FCTR2) 
   U10PH  = WAF_10M*U(I,J,K)
   V10PH  = WAF_10M*V(I,J,K)
ENDIF

!if (i==41 .and. j==41) then
!print 1116,k,kwind,zwfds,zc(kwind),u(i,j,k),u(i,j,kwind),waf_10m
!endif
!1116 format('(vege,cruzwind)',1x,2(I3),1x,5(e15.5))

UMAG = SQRT(U10PH**2 + V10PH**2)*MPS_TO_KPH !wind magnitude at 10 m above canopy, km/hr

!Theta_elps, after adjustment below, is angle of direction (0 to 2pi) of highest spread rate
!0<=theta_elps<=2pi as measured clockwise from Y-axis. ATAN2(y,x) is the angle, measured in the
!counterclockwise direction, between the positive x-axis and the line through (0,0) and (x,y)
!positive x-axis  

!Note, unlike the Rothermel ROS case, the slope is assumed to be zero at this point.
!THETA_ELPS(I,J) = ATAN2(V10PH,U10PH)
        
!The following two lines convert ATAN2 output to compass system (0 to 2 pi CW from +Y-axis)
!THETA_ELPS(I,J) = PIO2 - THETA_ELPS(I,J)
!IF (THETA_ELPS(I,J) < 0.0_EB) THETA_ELPS(I,J) = 2.0_EB*PI + THETA_ELPS(I,J)

!Probability of crowning
GMAX = 4.236_EB + 0.357_EB*UMAG - 0.71_EB*FSG - 0.331*EFFM 
IF (SFC <= 1.0_EB)                     G = GMAX - 4.613_EB
IF (SFC >  1.0_EB .AND. SFC <= 2.0_EB) G = GMAX - 1.856_EB
IF (SFC >  2.0_EB)                     G = GMAX
EXPG = EXP(G)
PROB = EXPG/(1.0_EB + EXPG)

CROSA = 11.02_EB*(UMAG**0.9_EB)*(CBD**0.19_EB)*EXP(-0.17_EB*EFFM)
CROSP = CROSA*EXP(-0.3333_EB*CROSA*CBD)

MIMIC_CRUZ_METHOD: IF (PROB_CROWN <= 1._EB) THEN
!Compute ROS if crown fire exists
  IF(PROB >= PROB_CROWN) THEN

!Theta_elps, after adjustment below, is angle of direction (0 to 2pi) of highest spread rate
!0<=theta_elps<=2pi as measured clockwise from Y-axis. ATAN2(y,x) is the angle, measured in the
!counterclockwise direction, between the positive x-axis and the line through (0,0) and (x,y)
!positive x-axis  

!Note, unlike the Rothermel ROS case, the slope is assumed to be zero at this point so the direction
!of local spread is dependent only on the wind direction.
!This means there is an inconsistency in the handling of the direction of spread in portions of the 
!fireline where PROB >= PROB_CROWN versus PROB < PROB_CROWN. In the latter case the influence of the
!slope on the local direction of spread is accounted for and wil be based on the Rothermel ROS_HEAD.
THETA_ELPS(I,J) = ATAN2(V10PH,U10PH)
        
!The following two lines convert ATAN2 output to compass system (0 to 2 pi CW from +Y-axis)
THETA_ELPS(I,J) = PIO2 - THETA_ELPS(I,J)
IF (THETA_ELPS(I,J) < 0.0_EB) THETA_ELPS(I,J) = 2.0_EB*PI + THETA_ELPS(I,J)

    CFB_LS(I,J) = 0.0_EB
    CRLOAD = CBD*VERT_CANOPY_EXTENT
    CAC = CROSA*CBD/3._EB
    IF (CAC >= 1._EB) THEN !active crown fire
      ROS_HEAD_CROWN(I,J) = CROSA*MPM_TO_MPS !convert m/min to m/s
      ROS_HEAD(I,J) = ROS_HEAD_CROWN(I,J)
      CFB_LS(I,J) = CRLOAD
    ELSE !passive crown fire
      ROS_HEAD_CROWN(I,J) = CROSP*MPM_TO_MPS
      ROS_HEAD(I,J) = ROS_HEAD_CROWN(I,J)
      CFB_LS(I,J) = CRLOAD*MAX(1.0_EB, (PROB - PROB_CROWN)/(1._EB - PROB_CROWN))
    ENDIF
  ENDIF
ENDIF MIMIC_CRUZ_METHOD

PROB_MIN_MAX_METHOD: IF (PROB_CROWN > 1._EB)  THEN !use 
!Compute head fire rate of spread 
  IF (PROB < PROB_PASSIVE) THEN
    IF(SURFACE_FIRE_HEAD_ROS_MODEL .EQ. 'CRUZ') THEN
      ROS_HEAD_CROWN(I,J) = CROSP*MPM_TO_MPS !else use surface ROS specified in input file
      ROS_HEAD(I,J) = ROS_HEAD_CROWN(I,J)
    ENDIF
  ENDIF
  IF (PROB >= PROB_PASSIVE .AND. PROB < PROB_ACTIVE) THEN
    ROS_HEAD_CROWN(I,J) = CROSP*MPM_TO_MPS 
    ROS_HEAD(I,J) = ROS_HEAD_CROWN(I,J)
  ENDIF
  IF (PROB >= PROB_ACTIVE) THEN
    ROS_HEAD_CROWN(I,J) = CROSA*MPM_TO_MPS
    ROS_HEAD(I,J) = ROS_HEAD_CROWN(I,J)
  ENDIF

!Compute crown fraction burned (kg/m^2) for use in QCONF (heat input to atmosphere)
  CFB_LS(I,J) = 0.0_EB
  IF (PROB >= PROB_PASSIVE) THEN
    CRLOAD = CBD*VERT_CANOPY_EXTENT
    CFB_LS(I,J) = CRLOAD*MAX(1.0_EB, (PROB - PROB_PASSIVE)/(PROB_ACTIVE - PROB_PASSIVE))
  ENDIF
ENDIF PROB_MIN_MAX_METHOD

!if (i==41 .and. j==41) then
!print 1117,prob,ros_head(i,j)
!print*,'========================='
!endif
!1117 format('(vege,cruzROS)',1x,2(e15.5))

!Store crown fire probability values
CRUZ_CROWN_PROB(I,J) = PROB 

END SUBROUTINE CRUZ_CROWN_FIRE_HEADROS


!************************************************************************************************
SUBROUTINE LEVEL_SET_FIREFRONT_PROPAGATION(T_CFD,DT,NM)
!************************************************************************************************
!
! Time step the scaler field PHI_LS. 
!
USE PHYSICAL_FUNCTIONS, ONLY : GET_MASS_FRACTION,GET_SPECIFIC_HEAT
INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: DT,T_CFD
INTEGER :: J_FLANK,I,II,IIG,IOR,IPC,IW,J,JJ,JJG,KK,KKG
INTEGER :: IDUM,JDUM
!LOGICAL :: IGNITION = .FALSE.
REAL(EB) :: BURNTIME,BT,FB_TIME_FCTR,FLI,HEAD_WIDTH_FCTR,GRIDCELL_FRACTION,GRIDCELL_TIME, &
            IGNITION_WIDTH_Y,RFIREBASE_TIME,RGRIDCELL_TIME,ROS_FLANK1,ROS_MAG,SHF,TE_TIME_FACTOR,TIME_LS_LAST, &
            TOTAL_FUEL_LOAD,VERT_CANOPY_EXTENT
REAL(EB) :: COSDPHIU,DPHIDX,DPHIDY,DPHIDOTU,DPHIMAG,XI,YJ,ZK,RCP_GAS,TE_HRRPUV,TE_HRR_TOTAL
!REAL(EB) :: PHI_CHECK
REAL(FB) :: TIME_LS_OUT

TYPE (WALL_TYPE),     POINTER :: WC =>NULL()
TYPE (SURFACE_TYPE),  POINTER :: SF =>NULL()

CALL POINT_TO_MESH(NM)

IF (.NOT. VEG_LEVEL_SET_UNCOUPLED .AND. .NOT. VEG_LEVEL_SET_COUPLED) RETURN

!--- Initialize variables
HEAD_WIDTH_FCTR  = 1._EB
IGNITION_WIDTH_Y = 1
J_FLANK          = 1
ROS_FLANK1       = 0._EB

IF (VEG_LEVEL_SET_COUPLED) THEN 
 Q       = 0.0_EB !HRRPUV array
!DT_LS   = MESHES(NM)%DT
 DT_LS   = DT
 TIME_LS = T_CFD
 T_FINAL = TIME_LS + DT_LS
ENDIF

IF (VEG_LEVEL_SET_UNCOUPLED) THEN
!DT_LS   = MESHES(NM)%DT
 DT_LS   = DT
 TIME_LS = T_CFD
 T_FINAL = TIME_LS + DT_LS 
ENDIF

!IF (NM==1) WRITE(LU_OUTPUT,'(A,1(I2),2x,3(E12.4))')'vege: nm,dt_ls,time_ls,t_final',nm,dt_ls,time_ls,t_final
!
!-- Time step solution using second order Runge-Kutta -----------------------
!
PHI_LS = LSET_PHI(0:IBP1,0:JBP1)

DO WHILE (TIME_LS < T_FINAL)

!
!-- Find flank-to-flank distance at base of fire assume symmetry about ymid and
!   define spread rate based on AU head fire width dependence
 IF (.NOT. LSET_ELLIPSE) THEN

!--------------------- Specific to AU grassland fuel experiments --------------------
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
!-----------------------------------------------------------------------------------------

 TIME_LS_LAST = TIME_LS

!----------------- Output time steps with increasing step (as in FDS)-------------------------
!IF ( (TIME_LS<=10.0_EB) .OR. (SUMTIME_LS > 100.0_EB) ) THEN
! SUMTIME_LS = 0._EB
! WRITE(LU_OUTPUT,*)'vege:LS:-------------------------------------'
! WRITE(LU_OUTPUT,*)'vege:LS:time_ls',time_ls
! WRITE(LU_OUTPUT,*)'vege:ls:dt',dt_ls
! WRITE(LU_OUTPUT,*)'vege:LS:HW,ros_h',head_width(nx_ls/2,ny_ls/2),ros_head(nx_ls/2,ny_ls/2)
! WRITE(LU_OUTPUT,*)'vege:LS:ros_f',ros_flank(nx_ls/2,ny_ls/2)
! WRITE(LU_OUTPUT,*)'vege:LS:max_ros for time stepping',dyn_sr_max
!ENDIF
!------------------------------------------------------------------------------------------------

!print*,'vege: ros_head',ros_head
!print*,'vege: crwn_prob',cruz_crown_prob

 WALL_CELL_LOOP1: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
  WC  => WALL(IW)
  IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE WALL_CELL_LOOP1
  SF  => SURFACE(WC%SURF_INDEX)

  !WC%ONE_D%QCONF    = 0.0_EB
  !WC%LSET_FIRE = .FALSE.
  II  = WC%ONE_D%II 
  JJ  = WC%ONE_D%JJ 
  KK  = WC%ONE_D%KK 
  IIG = WC%ONE_D%IIG
  JJG = WC%ONE_D%JJG
  KKG = WC%ONE_D%KKG
  IOR = WC%ONE_D%IOR
  
!-Ignite landscape at user specified location(s) and time(s) to originate Level Set fire front propagation
  IF (SF%VEG_LSET_IGNITE_TIME > 0.0_EB .AND. SF%VEG_LSET_IGNITE_TIME < DT_LS) THEN
    PHI_LS(IIG,JJG) = PHI_MAX_LS 
    BURN_TIME_LS(IIG,JJG) = 99999.0_EB
  ENDIF
  IF (SF%VEG_LSET_IGNITE_TIME >= TIME_LS .AND. SF%VEG_LSET_IGNITE_TIME <= TIME_LS + DT_LS) THEN 
    PHI_LS(IIG,JJG) = PHI_MAX_LS 
    BURN_TIME_LS(IIG,JJG) = 99999.0_EB
  ENDIF

  IF (.NOT. SF%VEG_LSET_SPREAD) CYCLE WALL_CELL_LOOP1

! --- Compute dot product between normal to fireline and wind direction. If location on fire perimeter is between the flank
!     and backing fires, then skip computation of crown fire ROS and use already computed surface fire ROS

  IF (SF%VEG_LSET_CROWN_FIRE_HEAD_ROS_MODEL=='CRUZ' .AND. VEG_LEVEL_SET_UNCOUPLED) THEN

    CALL ROTH_WINDANDSLOPE_COEFF_HEADROS(NM,IIG,JJG,KKG,SF%VEG_LSET_BETA,SF%VEG_LSET_SURF_HEIGHT,            &
      SF%VEG_LSET_CANOPY_HEIGHT,SF%VEG_LSET_SIGMA,SF%VEG_LSET_ROTH_ZEROWINDSLOPE_ROS, &
      SF%VEG_LSET_WAF_UNSHELTERED,SF%VEG_LSET_WAF_SHELTERED)

    DPHIDX = PHI_LS(IIG,JJG) - PHI_LS(IIG-1,JJG)
    DPHIDY = PHI_LS(IIG,JJG) - PHI_LS(IIG  ,JJG-1)
    DPHIDOTU = DPHIDX*U(IIG,JJG,KKG) + DPHIDY*V(IIG,JJG,KKG)
    IF (DPHIDOTU == 0.0_EB) THEN !account for zero wind speed
      COSDPHIU = 1._EB !use surface ROS
    ELSE
      UMAG     = SQRT(U(IIG,JJG,KKG)**2 + V(IIG,JJG,KKG)**2)
      DPHIMAG  = SQRT(DPHIDX**2 + DPHIDY**2)
      COSDPHIU = -DPHIDOTU/(UMAG*DPHIMAG) !minus sign to account for direction caused by phi=1,-1 in burned,unburned 
    ENDIF
    IF (COSDPHIU >= COS(SF%VEG_LSET_CROWNFIRE_ANGLE*PI/180._EB)) THEN  
      VERT_CANOPY_EXTENT = SF%VEG_LSET_CANOPY_HEIGHT - SF%VEG_LSET_SURF_HEIGHT - SF%VEG_LSET_FUEL_STRATA_GAP
      CALL CRUZ_CROWN_FIRE_HEADROS(NM,IIG,JJG,KKG,SF%VEG_LSET_CANOPY_BULK_DENSITY,SF%VEG_LSET_SURF_EFFM,     &
             SF%VEG_LSET_FUEL_STRATA_GAP,SF%VEG_LSET_SURF_LOAD,SF%VEG_LSET_CRUZ_PROB_PASSIVE,                    &
             SF%VEG_LSET_CRUZ_PROB_ACTIVE,SF%VEG_LSET_CRUZ_PROB_CROWN,SF%VEG_LSET_SURFACE_FIRE_HEAD_ROS_MODEL,   &
             VERT_CANOPY_EXTENT,SF%VEG_LSET_CANOPY_HEIGHT)
    ENDIF

  ENDIF


!
!****** Update quantities used in the spread rate computation if the Level Set and CFD computation are coupled.
!
  IF_CFD_COUPLED: IF (VEG_LEVEL_SET_COUPLED) THEN

    U_LS(IIG,JJG) = U(IIG,JJG,KKG)
    V_LS(IIG,JJG) = V(IIG,JJG,KKG)

    IF_ELLIPSE: IF (SF%VEG_LSET_ELLIPSE) THEN

!---Ellipse assumption with AU grassland head fire ROS for infinite head width
      IF (SF%VEG_LSET_SURFACE_FIRE_HEAD_ROS_MODEL=='AU GRASS' .AND. .NOT. SF%VEG_LSET_BURNER) &
        CALL AUGRASS_HEADROS(NM,IIG,JJG,KKG,SF%VEG_LSET_SURF_EFFM)

!---Ellipse assumption with Rothermel head fire ROS (== FARSITE)
      IF (SF%VEG_LSET_SURFACE_FIRE_HEAD_ROS_MODEL=='ROTHERMEL' .AND. .NOT. SF%VEG_LSET_BURNER) THEN 
        CALL ROTH_WINDANDSLOPE_COEFF_HEADROS(NM,IIG,JJG,KKG,SF%VEG_LSET_BETA,SF%VEG_LSET_SURF_HEIGHT,            &
           SF%VEG_LSET_CANOPY_HEIGHT,SF%VEG_LSET_SIGMA,SF%VEG_LSET_ROTH_ZEROWINDSLOPE_ROS, &
           SF%VEG_LSET_WAF_UNSHELTERED,SF%VEG_LSET_WAF_SHELTERED)
      ENDIF

!--- Cruz et al. crown fire head fire ROS model
      IF (SF%VEG_LSET_CROWN_FIRE_HEAD_ROS_MODEL=='CRUZ' .AND. .NOT. SF%VEG_LSET_BURNER) THEN
! --- Compute dot product between normal to fireline and wind direction. If location on fire perimeter is between the flank
!     and backing fires, then skip computation of crown fire ROS and use already computed surface fire ROS
        DPHIDX = PHI_LS(IIG,JJG) - PHI_LS(IIG-1,JJG)
        DPHIDY = PHI_LS(IIG,JJG) - PHI_LS(IIG  ,JJG-1)
        DPHIDOTU = DPHIDX*U_LS(IIG,JJG) + DPHIDY*V_LS(IIG,JJG)
        IF (DPHIDOTU == 0.0_EB) THEN
          COSDPHIU = 1._EB !use surface ROS
        ELSE
          UMAG     = SQRT(U_LS(IIG,JJG)**2 + V_LS(IIG,JJG)**2)
          DPHIMAG  = SQRT(DPHIDX**2 + DPHIDY**2)
          COSDPHIU = -DPHIDOTU/(UMAG*DPHIMAG) !minus sign to account for direction caused by phi=1,-1 in burned,unburned 
        ENDIF
        IF (COSDPHIU >= COS(SF%VEG_LSET_CROWNFIRE_ANGLE*PI/180._EB)) THEN  
          VERT_CANOPY_EXTENT = SF%VEG_LSET_CANOPY_HEIGHT - SF%VEG_LSET_SURF_HEIGHT - SF%VEG_LSET_FUEL_STRATA_GAP
          CALL CRUZ_CROWN_FIRE_HEADROS(NM,IIG,JJG,KKG,SF%VEG_LSET_CANOPY_BULK_DENSITY,SF%VEG_LSET_SURF_EFFM,     &
             SF%VEG_LSET_FUEL_STRATA_GAP,SF%VEG_LSET_SURF_LOAD,SF%VEG_LSET_CRUZ_PROB_PASSIVE,                    &
             SF%VEG_LSET_CRUZ_PROB_ACTIVE,SF%VEG_LSET_CRUZ_PROB_CROWN,SF%VEG_LSET_SURFACE_FIRE_HEAD_ROS_MODEL,   &
             VERT_CANOPY_EXTENT,SF%VEG_LSET_CANOPY_HEIGHT)
        ENDIF
      ENDIF

    ENDIF IF_ELLIPSE

!--- Compute heat flux into atmosphere
    GRIDCELL_TIME  = 0.0_EB
    RFIREBASE_TIME = 1.0_EB/SF%VEG_LSET_FIREBASE_TIME
    ROS_MAG = SQRT(SR_X_LS(IIG,JJG)**2 + SR_Y_LS(IIG,JJG)**2)
    IF(ROS_MAG > 0.0_EB) THEN
      GRIDCELL_TIME = SQRT(DX(IIG)**2 + DY(JJG)**2)/ROS_MAG
      RGRIDCELL_TIME = 1.0_EB/GRIDCELL_TIME
      GRIDCELL_FRACTION = MIN(1.0_EB,SF%VEG_LSET_FIREBASE_TIME*RGRIDCELL_TIME) !assumes spread direction parallel to grid axes
    ENDIF
    BURNTIME = MAX(SF%VEG_LSET_FIREBASE_TIME,GRIDCELL_TIME) !assumes spread direction parallel to grid axes

    BT  = BURN_TIME_LS(IIG,JJG)
    SHF = 0.0_EB !surface heat flux W/m^2
    WC%LSET_FIRE = .FALSE.

!Determine surface heat flux for fire spread through grid cell. Account for fires with a depth that is smaller
!than the grid cell (GRIDCELL_FRACTION). Also account for partial presence of fire base as fire spreads into 
!and out of the grid cell (FB_TIME_FCTR).

    HRRPUA_OUT(IIG,JJG) = 0.0 !kW/m^2 
    TOTAL_FUEL_LOAD = SF%VEG_LSET_SURF_LOAD + CFB_LS(IIG,JJG)

    IF_FIRELINE_PASSAGE: IF (PHI_LS(IIG,JJG) >= -SF%VEG_LSET_PHIDEPTH .AND. .NOT. SF%VEG_LSET_BURNER .AND. &
                             .NOT. VEG_LEVEL_SET_BURNERS_FOR_FIRELINE) THEN 

      WC%LSET_FIRE = .TRUE.
!     TOTAL_FUEL_LOAD = SF%VEG_LSET_SURF_LOAD + CFB_LS(IIG,JJG)
!max surface heat flux, W/m^2
      SHF = (1.0_EB-SF%VEG_CHAR_FRACTION)*SF%VEG_LSET_HEAT_OF_COMBUSTION*TOTAL_FUEL_LOAD*RFIREBASE_TIME 
!if(iig==31 .and. jjg==49 .and. nm==5) print '(A,2x,5ES12.4)','time, charfrac,Hc,w,rfbt =', & 
!            t_cfd,sf%veg_char_fraction,sf%veg_lset_heat_of_combustion,total_fuel_load,rfirebase_time

!Grid cell > fire depth      
      IF (GRIDCELL_FRACTION < 1.0_EB) THEN
        SHF = SHF*GRIDCELL_FRACTION
        FB_TIME_FCTR = 1.0_EB
!       Fire entering cell
        IF (0.0_EB        <= BT .AND. BT <= SF%VEG_LSET_FIREBASE_TIME) FB_TIME_FCTR = BT*RFIREBASE_TIME
!       Fire exiting cell
        IF (GRIDCELL_TIME <  BT .AND. BT <= GRIDCELL_TIME + SF%VEG_LSET_FIREBASE_TIME) FB_TIME_FCTR = &
          1.0_EB - (BT - GRIDCELL_TIME)*RFIREBASE_TIME
!       Fire has left cell
        IF (BT > GRIDCELL_TIME + SF%VEG_LSET_FIREBASE_TIME) WC%LSET_FIRE = .FALSE.
!       WC%VEG_HEIGHT = SF%VEG_LSET_SURF_HEIGHT*(1._EB - BT/(SF%VEG_LSET_FIREBASE_TIME+GRIDCELL_TIME))
        WC%VEG_HEIGHT = 0.0_EB
        BURN_TIME_LS(IIG,JJG) = BURN_TIME_LS(IIG,JJG) + DT_LS

if(iig==46 .and. jjg==49 .and. nm==5) then 
   print '(A,2x,7ES12.4)','----LStime dx>fd, ros, bt, gct, fbt, fctr, shf =', &
            t_cfd,ros_mag,bt,gridcell_time,sf%veg_lset_firebase_time, &
                                                                           fb_time_fctr,-shf*0.001_EB
   print '(A,2x,1ES12.4,L2)','cell fract, lset_fire =',gridcell_fraction,wc%lset_fire
!   print '(A,2x,4ES12.4,L2)','----time dtdx>dtfb, shf kW/m2, fb_time_fctr, cell fract, lset_fire =', &
!       t_cfd,-shf*0.001_EB,fb_time_fctr,gridcell_fraction,wc%lset_fire
!   print '(A,2x,7ES12.4)','ros, hcomb, fuel load, fbt, cfb, probcruz, probin', &
!       ros_mag,sf%veg_lset_heat_of_combustion,total_fuel_load,firebase_time,cfb_ls(iig,jjg), &
!       CRUZ_CROWN_PROB(IIG,JJG),sf%veg_lset_cruz_prob_crown
endif
        WC%VEG_LSET_SURFACE_HEATFLUX = -SHF*FB_TIME_FCTR
      ENDIF

!Grid cell <= fire depth      
      IF (GRIDCELL_FRACTION >= 1.0_EB) THEN
        FB_TIME_FCTR = 1.0_EB
!       Fire entering cell
        IF (0.0_EB        <= BT .AND. BT <= GRIDCELL_TIME) FB_TIME_FCTR = BT*RGRIDCELL_TIME
!       Fire exiting cell
        IF (SF%VEG_LSET_FIREBASE_TIME <  BT .AND. BT <= GRIDCELL_TIME + SF%VEG_LSET_FIREBASE_TIME) FB_TIME_FCTR = &
          1.0_EB - (BT - SF%VEG_LSET_FIREBASE_TIME)*RGRIDCELL_TIME
!       Fire has left cell
        IF (BT > GRIDCELL_TIME + SF%VEG_LSET_FIREBASE_TIME) WC%LSET_FIRE = .FALSE.
!       WC%VEG_HEIGHT = SF%VEG_LSET_SURF_HEIGHT*(1._EB - BT/(SF%VEG_LSET_FIREBASE_TIME+GRIDCELL_TIME))
        WC%VEG_HEIGHT = 0.0_EB
        BURN_TIME_LS(IIG,JJG) = BURN_TIME_LS(IIG,JJG) + DT_LS

!if(iig==46 .and. jjg==49 .and. nm==5) then 
!   print '(A,2x,7ES12.4)','++++time dx<=fd, ros, bt, gct, fbt, fctr, shf =', &
!         t_cfd,ros_mag,bt,gridcell_time,sf%veg_lset_firebase_time,fb_time_fctr,-shf*0.001_EB
!   print '(A,2x,2ES12.4,L2)','cell fract, lset_fire,W =',gridcell_fraction,total_fuel_load,wc%lset_fire
!   print '(A,2x,7ES12.4)','ros, hcomb, fuel load, fbt, cfb, probcruz, probin', &
!                           ros_mag,sf%veg_lset_heat_of_combustion,total_fuel_load,  &
!                           firebase_time,cfb_ls(iig,jjg),CRUZ_CROWN_PROB(IIG,JJG),sf%veg_lset_cruz_prob_crown
!endif
!print '(A,2x,2ES12.4)','----LS shf, fb_time_fctr =', shf,fb_time_fctr
        WC%VEG_LSET_SURFACE_HEATFLUX = -SHF*FB_TIME_FCTR
      ENDIF

!     IF (WC%LSET_FIRE) HRRPUA_OUT(IIG,JJG) = -WC%VEG_LSET_SURFACE_HEATFLUX*0.001 !kW/m^2 for Smokeview output

    ENDIF IF_FIRELINE_PASSAGE

     
! Stop burning if the fire front residence time is exceeded
    IF (PHI_LS(IIG,JJG) >= -SF%VEG_LSET_PHIDEPTH .AND. .NOT. WC%LSET_FIRE) THEN
        WC%VEG_LSET_SURFACE_HEATFLUX = 0.0_EB
        WC%VEG_HEIGHT = 0.0_EB 
        BURN_TIME_LS(IIG,JJG) = 999999999._EB
    ENDIF

!if(x(iig)==29 .and. y(jjg)==1) then 
!    print '(A,2x,8ES12.4)','time,phi_ls,burn_time_ls,shf,Rx, Ry, cell frac, burntime =',t_cfd,phi_ls(iig,jjg), &
!           burn_time_ls(iig,jjg),wc%veg_lset_surface_heatflux/gridcell_fraction, &
!           sr_x_ls(iig,jjg),sr_y_ls(iig,jjg),gridcell_fraction,burntime
!    print '(A,2x,ES12.4)','phi_ls(x+dx,y)',phi_ls(iig+1,jjg)
!endif

!-- Burner placement as explicitly specified (location,timing, etc.) in the input file
    IF (SF%VEG_LSET_BURNER .AND. .NOT. VEG_LEVEL_SET_BURNERS_FOR_FIRELINE) THEN
      IF (TIME_LS >= SF%VEG_LSET_BURNER_TIME_ON .AND. TIME_LS <= SF%VEG_LSET_BURNER_TIME_OFF) THEN
        WC%VEG_LSET_SURFACE_HEATFLUX = -SF%HRRPUA
        PHI_LS(IIG,JJG) = PHI_MAX_LS 
        WC%LSET_FIRE = .TRUE.
      ELSE
        WC%VEG_LSET_SURFACE_HEATFLUX = 0.0_EB
        WC%LSET_FIRE = .FALSE.
      ENDIF
    ENDIF

!-- Placement of burners as determined from reading a file containing a 2D array with HRRPUA values at times
    IF (VEG_LEVEL_SET_BURNERS_FOR_FIRELINE) THEN
!max surface heat flux, W/m^2
      SHF = (1.0_EB-SF%VEG_CHAR_FRACTION)*SF%VEG_LSET_HEAT_OF_COMBUSTION*TOTAL_FUEL_LOAD*RFIREBASE_TIME 
      WC%LSET_FIRE = .FALSE.
      WC%VEG_LSET_SURFACE_HEATFLUX = 0.0_EB
      IF (REAL(TIME_LS,FB) <=  LSET_TIME_HRRPUA_BURNER) THEN
        IF(HRRPUA_IN(IIG,JJG) > 0.0_FB) THEN
          IF (HRRPUA_IN(IIG,JJG)*1000._EB >= SF%VEG_LSET_HRRPUA_MINIMUM_FRAC*SHF) &
             WC%VEG_LSET_SURFACE_HEATFLUX = -HRRPUA_IN(IIG,JJG)*1000._EB !W/m^2
          PHI_LS(IIG,JJG) = PHI_MAX_LS 
          WC%LSET_FIRE = .TRUE.
        ENDIF
      ELSE
        DO WHILE (REAL(TIME_LS,FB) > LSET_TIME_HRRPUA_BURNER)
          READ(LU_SLCF_LS(7),END=111) LSET_TIME_HRRPUA_BURNER 
          READ(LU_SLCF_LS(7)) ((HRRPUA_IN(IDUM,JDUM),IDUM=0,IBAR),JDUM=0,JBAR) 
        ENDDO
111     IF(HRRPUA_IN(IIG,JJG) > 0.0) THEN
          IF (HRRPUA_IN(IIG,JJG)*1000._EB >= SF%VEG_LSET_HRRPUA_MINIMUM_FRAC*SHF) & 
            WC%VEG_LSET_SURFACE_HEATFLUX = -HRRPUA_IN(IIG,JJG)*1000._EB !W/m^2
          PHI_LS(IIG,JJG) = PHI_MAX_LS 
          WC%LSET_FIRE = .TRUE.
        ENDIF
      ENDIF
!if(iig==46 .and. jjg==49 .and. nm==5) print '(A,2x,3ES12.4,L2)','SHF,HRRPUA_IN,WC%SHF=', &
!                                             shf,hrrpua_in(iig,jjg)*1000._EB,wc%veg_lset_surface_heatflux
    ENDIF

    IF (VEG_LEVEL_SET_SURFACE_HEATFLUX) WC%ONE_D%QCONF = WC%VEG_LSET_SURFACE_HEATFLUX
    IF (VEG_LEVEL_SET_THERMAL_ELEMENTS) SF%DT_INSERT = DT_LS !**** is this correct?
    IF (WC%LSET_FIRE) HRRPUA_OUT(IIG,JJG) = -WC%VEG_LSET_SURFACE_HEATFLUX*0.001 !kW/m^2 for Smokeview output

!---Drag constant can vary with height, if hveg > dzgrid
!print '(A,1x,5I3)','velocity_bc_index,nm,i,j,k',sf%velocity_bc_index,nm,iig,jjg,kkg
    VEG_DRAG(IIG,JJG,:) = 0.0_EB
!   IF (WC%VEG_HEIGHT > 0.0_EB) THEN
!    DO KGRID=1,8
!     IF (Z(KGRID) <= WC%VEG_HEIGHT) VEG_DRAG(IIG,JJG,KGRID)= SF%VEG_DRAG_INI
!     IF (Z(KGRID) >  WC%VEG_HEIGHT .AND. Z(KGRID-1) < WC%VEG_HEIGHT) VEG_DRAG(IIG,JJG,KGRID)= &
!                      SF%VEG_DRAG_INI*(WC%VEG_HEIGHT-Z(KGRID-1))/(Z(KGRID)-Z(KGRID-1))
!    ENDDO
!   ENDIF

!   IF (PHI_LS(IIG,JJG) <= SF%VEG_LSET_PHIDEPTH .AND. PHI_LS(IIG,JJG) >= -SF%VEG_LSET_PHIDEPTH) THEN 
!    WC%ONE_D%TMP_F = 373._EB
!    WC%ONE_D%QCONF = SF%VEG_LSET_QCON
!   ENDIF

  ENDIF IF_CFD_COUPLED

! Save Time of Arrival (TOA), Rate of Spread components, Fireline Intensity, etc.
  IF (PHI_LS(IIG,JJG) >= -SF%VEG_LSET_PHIDEPTH .AND. TOA(IIG,JJG) <= -1.0_EB) THEN 
    TOA(IIG,JJG)=TIME_LS
    ROS_X_OUT(IIG,JJG) = SR_X_LS(IIG,JJG)
    ROS_Y_OUT(IIG,JJG) = SR_Y_LS(IIG,JJG)
    TOTAL_FUEL_LOAD = SF%VEG_LSET_SURF_LOAD + CFB_LS(IIG,JJG)
    FLI = SQRT(SR_Y_LS(IIG,JJG)**2 + SR_X_LS(IIG,JJG)**2)*(SF%VEG_LSET_HEAT_OF_COMBUSTION*0.001_EB)* &
         (1.0_EB-SF%VEG_CHAR_FRACTION)*TOTAL_FUEL_LOAD
    FLI_OUT(IIG,JJG) = FLI !kW/m
    CRUZ_CROWN_PROB_OUT(IIG,JJG) = CRUZ_CROWN_PROB(IIG,JJG)
  ENDIF

 ENDDO WALL_CELL_LOOP1

 IF (ANY(PHI_LS==PHI_MAX_LS)) LSET_IGNITION = .TRUE.

!2nd order Runge-Kutta time steppping of level set equation
 PHI_TEMP_LS = PHI_LS 

!RK Stage 1
 RK2_PREDICTOR_LS = .TRUE.
!CALL LEVEL_SET_PERIMETER_SPREAD_RATE(NM)
 CALL LEVEL_SET_ADVECT_FLUX(NM)
 PHI1_LS = PHI_LS - DT_LS*FLUX0_LS
!print*,'max flux0',rk2_predictor_ls,maxval(abs(flux0_ls))

!RK Stage2
 RK2_PREDICTOR_LS = .FALSE.
 MAG_SR_OUT       = 0.0_EB
 CALL LEVEL_SET_PERIMETER_SPREAD_RATE(NM) 
 CALL LEVEL_SET_ADVECT_FLUX(NM)
 PHI_LS = PHI_LS - 0.5_EB*DT_LS*(FLUX0_LS + FLUX1_LS)

!The following is done here instead of in Stage 1 RK so updated ROS can be used in coupled LS when
!determining if fire residence time is shorter than a time step.
 CALL LEVEL_SET_PERIMETER_SPREAD_RATE(NM)

! Account for heat released by thermal elements (if present). Thermal elements are inserted and
! LP%LSET_HRRPUV is determined in part.f90
IF (VEG_LEVEL_SET_THERMAL_ELEMENTS) THEN
  RCP_GAS    = 0.001_EB !1/(J/kg/K)
  TE_HRR_TOTAL = 0.0_EB
  PARTICLE_LOOP: DO I=1,NLP
    LP  => LAGRANGIAN_PARTICLE(I)
    IPC =  LP%CLASS_INDEX
    LPC  => LAGRANGIAN_PARTICLE_CLASS(IPC)
    IF(.NOT. LP%LSET_THERMAL_ELEMENT) CYCLE PARTICLE_LOOP
    CALL GET_IJK(LP%X,LP%Y,LP%Z,NM,XI,YJ,ZK,II,JJ,KK)
!   LP%LSET_HRRPUV = 0.01_EB !W/m^3 
!   TE_TIME_FACTOR = 1.0_EB - (T_FINAL-LP%T)/LPC%TE_BURNTIME !linear decay with time
    TE_TIME_FACTOR = 1.0_EB !no decay with time
!   TE_TIME_FACTOR = MAX(0.0_EB,TE_TIME_FACTOR)
!   LPC%RGB = (/255,0,0/)
!   IF(TE_TIME_FACTOR == 0.0_EB) LPC%RGB = (/0,0,0/)
!   TE_HRR_TOTAL  = TE_HRR_TOTAL + TE_TIME_FACTOR*LP%LSET_HRRPUV*DX(II)*DY(JJ)*DZ(KK)
    TE_HRRPUV =  TE_TIME_FACTOR*LP%LSET_HRRPUV
    IF (T_CFD - LP%T_INSERT > LPC%TE_BURNTIME) TE_HRRPUV = 0.0_EB
    IF (T_CFD - LP%T_INSERT + DT_LS > LPC%TE_BURNTIME) TE_HRRPUV = TE_HRRPUV + & 
                                      TE_HRRPUV*(LPC%TE_BURNTIME-(T_CFD-LP%T_INSERT))/DT_LS
    Q(II,JJ,KK) = Q(II,JJ,KK) + TE_HRRPUV
!   D_LAGRANGIAN(II,JJ,KK) = D_LAGRANGIAN(II,JJ,KK)  +  &
!                             TE_TIME_FACTOR*LP%LSET_HRRPUV*RCP_GAS/(RHO(II,JJ,KK)*TMP(II,JJ,KK))
  ENDDO PARTICLE_LOOP
!  CALL REMOVE_PARTICLES(T_CFD,NM)
!print '(A,2x,I3,2x,2ES12.4)','Mesh, Time, TE_HRR_TOTAL (kW) = ',NM,T_CFD,TE_HRR_TOTAL*0.001_EB
ENDIF

!print*,'min,max phi_ls',minval(phi_ls),maxval(phi_ls)

!print*,'max flux0,flux1',maxval(abs(flux0_ls)),maxval(abs(flux1_ls))

!Variable Time Step for simulations uncoupled from the CFD computation
!IF (VEG_LEVEL_SET_UNCOUPLED) THEN
! IF (.NOT. VEG_LEVEL_SET) THEN
!
!   PHI_CHECK = MAXVAL(ABS(PHI_LS - PHI_TEMP_LS)) !Find max change in phi
! 
!   IF (LSET_IGNITION) THEN
!     ! If any phi values change by more than 0.5, or all change less
!     ! than 0.1 (both values are arbitrary), during one time step,
!     ! then adjust time step accordingly.
!     
!     IF (PHI_CHECK > 0.5_EB) THEN
!         ! Cut time step in half and cycle the do-while loop
!         DT_COEF_LS = 0.5_EB * DT_COEF_LS 
!!        DT_LS = DT_COEF_LS * MIN(DX_LS,DY_LS)/DYN_SR_MAX_LS
!         DT_LS = DT_COEF_LS * MIN(DX(1),DY(1))/DYN_SR_MAX_LS
!         DT_LS = MIN(DT_LS,100._EB)
!         PHI_LS = PHI_TEMP_LS ! Exchange for previous phi and cycle
!         print '(A,1x,E13.5,1x,A,i3)',"Halving time step, dt=  ",dt_ls,' mesh ',nm
!         CYCLE 
!     ENDIF
! 
!     ! Increase time step by 1/4 if changes are small
!     IF (PHI_CHECK < 0.1_EB) DT_COEF_LS = DT_COEF_LS * 1.25_EB
!     
!     ! Dynamic Spread Rate Max
!     DYN_SR_MAX_LS = MAX(DYN_SR_MAX_LS,0.01_EB) ! dyn_sr_max must be g.t. zero
!!    DT_LS = DT_COEF_LS * MIN(DX_LS,DY_LS)/DYN_SR_MAX_LS
!     DT_LS = DT_COEF_LS * MIN(DX(1),DY(1))/DYN_SR_MAX_LS
!     DT_LS = MIN(DT_LS,100._EB)
!     
!   ENDIF
! ENDIF
 

 TIME_LS = TIME_LS + DT_LS
 SUMTIME_LS = SUMTIME_LS + DT_LS
 SUM_T_SLCF_LS = SUM_T_SLCF_LS + DT_LS

 LSET_PHI(0:IBP1,0:JBP1) = PHI_LS
!MESHES(1)%LSET_PHI(IBP1,:,1) = MESHES(2)%LSET_PHI(   1,:,1)
!MESHES(2)%LSET_PHI(   0,:,1) = MESHES(1)%LSET_PHI(IBAR,:,1)
!if (nm == 2) print 1113,nm,meshes(1)%lset_phi(ibar,14,1),meshes(2)%lset_phi(0,14,1)

!Runtime output of slice files containing level set variables for smokeview animation
 IF (SUM_T_SLCF_LS >= DT_OUTPUT_LS) THEN    
  SUM_T_SLCF_LS = 0._EB
  PHI_OUT(0:IBAR,0:JBAR) = PHI_LS(0:IBAR,0:JBAR)
  TIME_LS_OUT = TIME_LS
!-- PHI field
  WRITE(LU_SLCF_LS(1)) TIME_LS_OUT
!if (nm ==3) print*,'vege:phi',meshes(3)%lset_phi(12,0,1)
!if (nm ==3) then
! print*,'vege:phi',phi_ls(12,0)
! print*,'vege:phiout',phi_out(12,0)
!endif
!negative for consistency with wall thickness output from wfds and viz by Smokeview
  WRITE(LU_SLCF_LS(1)) ((-PHI_OUT(I,J),I=0,IBAR),J=0,JBAR) 
!if (nm==3) print '(A,1x,1i3,1x,2E13.5)','vege:outslf',nm,phi_ls(12,0),meshes(1)%lset_phi(12,jbp1,1)
!if (nm==1) print '(A,1x,1i3,1x,2E13.5)','vege:outslfout',nm,phi_out(12,ibar)
!-- Time of Arrival, s
  WRITE(LU_SLCF_LS(2)) TIME_LS_OUT
  TOA(0,0:JBAR) = TOA(1,0:JBAR) ; TOA(0:IBAR,0) = TOA(0:IBAR,1) !for Smokeview
  WRITE(LU_SLCF_LS(2)) ((TOA(I,J),I=0,IBAR),J=0,JBAR) 
!-- ROS magnitude, m/s
  WRITE(LU_SLCF_LS(3)) TIME_LS_OUT
  ROS_X_OUT(0,0:JBAR) = ROS_X_OUT(1,0:JBAR) ; ROS_X_OUT(0:IBAR,0) = ROS_X_OUT(0:IBAR,1)
  ROS_Y_OUT(0,0:JBAR) = ROS_Y_OUT(1,0:JBAR) ; ROS_Y_OUT(0:IBAR,0) = ROS_Y_OUT(0:IBAR,1)
  WRITE(LU_SLCF_LS(3)) ((SQRT(ROS_X_OUT(I,J)**2 + ROS_Y_OUT(I,J)**2),I=0,IBAR),J=0,JBAR) 
!-- Fireline intensity at time of fire arrival, kW/m^2
  WRITE(LU_SLCF_LS(4)) TIME_LS_OUT
  FLI_OUT(0,0:JBAR) = FLI_OUT(1,0:JBAR) ; FLI_OUT(0:IBAR,0) = FLI_OUT(0:IBAR,1) !for Smokeview
  WRITE(LU_SLCF_LS(4)) ((FLI_OUT(I,J),I=0,IBAR),J=0,JBAR) 
!-- HRRPUA at every slice output time, kW/m^2
  IF (VEG_LEVEL_SET_COUPLED) THEN
    WRITE(LU_SLCF_LS(5)) TIME_LS_OUT
    HRRPUA_OUT(0,0:JBAR) = HRRPUA_OUT(1,0:JBAR) ; HRRPUA_OUT(0:IBAR,0) = HRRPUA_OUT(0:IBAR,1) !for Smokeview
    WRITE(LU_SLCF_LS(5)) ((HRRPUA_OUT(I,J),I=0,IBAR),J=0,JBAR) 
  ENDIF
!-- Crown fire Probability (Cruz & Alexancer)
  WRITE(LU_SLCF_LS(6)) TIME_LS_OUT
  CRUZ_CROWN_PROB_OUT(0,0:JBAR) = CRUZ_CROWN_PROB_OUT(1,0:JBAR) !for Smokeview
  CRUZ_CROWN_PROB_OUT(0:IBAR,0) = CRUZ_CROWN_PROB_OUT(0:IBAR,1) !for Smokeview
  WRITE(LU_SLCF_LS(6)) ((CRUZ_CROWN_PROB_OUT(I,J),I=0,IBAR),J=0,JBAR) 
 ENDIF
!
ENDDO !While loop

!if (nm == 1) print*,'vegprop nm,lset_phi(ibp1,1)',nm,phi_ls(ibp1,1)
!if (nm == 2) print*,'vegprop nm,lset_phi(0,1)',nm,phi_ls(0,1)
!if (nm == 1) print*,'vegprop nm,m(1,2)lset_phi',nm,meshes(1)%lset_phi(ibp1,1,1),meshes(2)%lset_phi(0,1,1)
!if (nm == 2) print*,'vegprop nm,m(1,2)lset_phi',nm,meshes(1)%lset_phi(ibp1,1,1),meshes(2)%lset_phi(0,1,1)
!print 1113,nm,meshes(1)%lset_phi(ibar,25,1),meshes(2)%lset_phi(1,25,1)
1113 format('vegelsprop nm,lset_phi',1(i2),2x,2(E12.4))

!CLOSE(LU_SLCF_LS)

! ******  Write arrays to ascii file **************
!IF (VEG_LEVEL_SET_UNCOUPLED .AND. NM == 1) THEN
! CALL CPU_TIME(CPUTIME)
! LS_T_END = CPUTIME
! WRITE(LU_OUTPUT,*)'Uncoupled Level Set CPU Time: ',LS_T_END - LS_T_BEG
!ENDIF
!
!-- Output time of arrival
!LU_TOA_LS = GET_FILE_NUMBER()
!print*,'veg:toa_ls',lu_toa_ls
!OPEN(LU_TOA_LS,FILE='time_of_arrival.toa',STATUS='REPLACE')
!WRITE(LU_TOA_LS,'(I5)') NX_LS,NY_LS
!WRITE(LU_TOA_LS,'(F7.2)') XS,XF,YS,YF
!Write across row (TOA(1,1), TOA(1,2), ...) to match Farsite output
!IF (TIME_LS >= T_END) THEN
! print*,'veg:toaf_ls',lu_toa_ls
! WRITE(LU_TOA_LS,'(F7.2)') ((TOA(IDUM,JDUM),JDUM=1,NY_LS),IDUM=1,NX_LS)
! CLOSE(LU_TOA_LS)
!ENDIF

! Diagnostics at end of run
!OPEN(9998,FILE='Phi_S.txt',STATUS='REPLACE')
!WRITE(9998,'(I5)') NX_LS,NY_LS
!WRITE(9998,'(F7.2)') XS,XF,YS,YF
!WRITE(9998,'(F7.2)') ((PHI_S(IDUM,JDUM),JDUM=1,NY_LS),IDUM=1,NX_LS)
!CLOSE(9998)

!OPEN(9998,FILE='Phi_W.txt',STATUS='REPLACE')
!WRITE(9998,'(I5)') NX_LS,NY_LS
!WRITE(9998,'(F7.2)') XS,XF,YS,YF
!WRITE(9998,'(F7.2)') ((PHI_W(IDUM,JDUM),JDUM=1,NY_LS),IDUM=1,NX_LS)
!CLOSE(9998)

!OPEN(9998,FILE='alt.txt',STATUS='REPLACE')
!WRITE(9998,'(I5)') NX_LS,NY_LS
!WRITE(9998,'(F7.2)') XS,XF,YS,YF
!WRITE(9998,'(F10.5)') ((ZT(IDUM,JDUM),JDUM=1,NY_LS),IDUM=1,NX_LS)
!CLOSE(9998)

!OPEN(9998,FILE='DZTDX.txt',STATUS='REPLACE')
!WRITE(9998,'(I5)') NX_LS,NY_LS
!WRITE(9998,'(F7.2)') XS,XF,YS,YF
!WRITE(9998,'(F7.2)') ((DZTDX(IDUM,JDUM),JDUM=1,NY_LS),IDUM=1,NX_LS)
!CLOSE(9998)

!OPEN(9998,FILE='DZTDY.txt',STATUS='REPLACE')
!WRITE(9998,'(I5)') NX_LS,NY_LS
!WRITE(9998,'(F7.2)') XS,XF,YS,YF
!WRITE(9998,'(F7.2)') ((DZTDY(IDUM,JDUM),JDUM=1,NY_LS),IDUM=1,NX_LS)
!CLOSE(9998)

!OPEN(9998,FILE='Theta_Ellipse.txt',STATUS='REPLACE')
!WRITE(9998,'(I5)') NX_LS,NY_LS
!WRITE(9998,'(F7.2)') XS,XF,YS,YF
!WRITE(9998,'(F7.2)') ((Theta_Elps(IDUM,JDUM),JDUM=1,NY_LS),IDUM=1,NX_LS)
!CLOSE(9998)

!OPEN(9998,FILE='UMF.txt',STATUS='REPLACE')
!WRITE(9998,'(I5)') NX_LS,NY_LS
!WRITE(9998,'(F7.2)') XS,XF,YS,YF
!WRITE(9998,'(F7.2)') ((UMF(IDUM,JDUM),JDUM=1,NY_LS),IDUM=1,NX_LS)
!CLOSE(9998)

END SUBROUTINE LEVEL_SET_FIREFRONT_PROPAGATION

!************************************************************************************************
SUBROUTINE END_LEVEL_SET
!************************************************************************************************
!
! Output quantities at end of level set simulation
!
!INTEGER :: IDUM,JDUM
! Output time of arrival array
!WRITE(LU_TOA_LS,'(F7.2)') ((TOA(IDUM,JDUM),IDUM=1,NX_LS),JDUM=1,NY_LS)
!CLOSE(LU_TOA_LS)
!!WRITE(LU_ROS_LS,'(F7.2)') ((ROS_X_OUT(IDUM,JDUM),IDUM=1,NX_LS),JDUM=1,NY_LS), &
!!                          ((ROS_Y_OUT(IDUM,JDUM),IDUM=1,NX_LS),JDUM=1,NY_LS)
!WRITE(LU_ROSX_LS,'(F7.2)') ((ROS_X_OUT(IDUM,JDUM),IDUM=1,NX_LS),JDUM=1,NY_LS)
!CLOSE(LU_ROSX_LS)
!WRITE(LU_ROSY_LS,'(F7.2)') ((ROS_Y_OUT(IDUM,JDUM),IDUM=1,NX_LS),JDUM=1,NY_LS)
!CLOSE(LU_ROSY_LS)
!WRITE(LU_FLI_LS,'(F9.2)') ((FLI_OUT(IDUM,JDUM),IDUM=1,NX_LS),JDUM=1,NY_LS)
!CLOSE(LU_FLI_LS)
!WRITE(LU_CRWN_PROB_LS,'(F9.2)') ((CRUZ_CROWN_PROB_OUT(IDUM,JDUM),IDUM=1,NX_LS),JDUM=1,NY_LS)
!CLOSE(LU_CRWN_PROB_LS)

END SUBROUTINE END_LEVEL_SET
!
!************************************************************************************************
SUBROUTINE LEVEL_SET_PERIMETER_SPREAD_RATE(NM)
!************************************************************************************************
!
! Compute components of spread rate vector along fire perimeter
!
INTEGER, INTENT(IN) :: NM
INTEGER :: I,J,IM1,IP1,JM1,JP1
REAL(EB) :: COS_THETA_WIND,COS_THETA_SLOPE,COS_THETA_WIND_H,COS_THETA_WIND_B, &
            COS_THETA_SLOPE_H,COS_THETA_SLOPE_B,DPHIDX,DPHIDY,F_EAST,F_WEST,F_NORTH,F_SOUTH, &
            GRAD_SLOPE_DOT_NORMAL_FIRELINE,MAG_F,MAG_SR,MAG_U,WIND_DOT_NORMAL_FIRELINE,NEXP_WIND
REAL(EB) :: RAD_TO_DEGREE,DEGREES_SLOPE,SLOPE_FACTOR

!Variables for elliptical propagation model

REAL(EB) :: COS_THETA,SIN_THETA,XSF,YSF,UMF_DUM
REAL(EB) :: A_ELPS,A_ELPS2,AROS,BROS,B_ELPS2,B_ELPS,C_ELPS,DENOM,ROS_TMP,LB,LBD,HB
REAL(EB), DIMENSION(:) :: NORMAL_FIRELINE(2)

CALL POINT_TO_MESH(NM)
 
RAD_TO_DEGREE = 90._EB/ASIN(1._EB)

!NEXP_WIND = 2

IF (RK2_PREDICTOR_LS) PHI0_LS = PHI_LS
IF (.NOT. RK2_PREDICTOR_LS) PHI0_LS = PHI1_LS
SR_X_LS = 0.0_EB ; SR_Y_LS = 0.0_EB
DYN_SR_MAX_LS = 0.0_EB

FLUX_ILOOP: DO I = 1,NX_LS
  
  IM1=I-1 
  IP1=I+1 
! IF (I==1) IM1 = I
! IF (I==NX_LS) IP1 = I
  
  DO J = 1,NY_LS
    
   JM1=J-1
   JP1=J+1
!  IF (J==1) JM1 = J
!  IF (J==NX_LS) JP1 = J

   F_EAST  = 0.5_EB*( PHI0_LS(I,J) + PHI0_LS(IP1,J) )
   F_WEST  = 0.5_EB*( PHI0_LS(I,J) + PHI0_LS(IM1,J) )
   F_NORTH = 0.5_EB*( PHI0_LS(I,J) + PHI0_LS(I,JP1) )
   F_SOUTH = 0.5_EB*( PHI0_LS(I,J) + PHI0_LS(I,JM1) )
         
   DPHIDX = (F_EAST-F_WEST) * RDX(I) !IDX_LS
   DPHIDY = (F_NORTH-F_SOUTH) * RDY(J) !IDY_LS
   
   MAG_F = SQRT(DPHIDX**2 + DPHIDY**2)
   IF (MAG_F > 0._EB) THEN   !components of unit vector normal to PHI contours
        NORMAL_FIRELINE(1) = -DPHIDX/MAG_F
        NORMAL_FIRELINE(2) = -DPHIDY/MAG_F
        XSF =  DPHIDY
        YSF = -DPHIDX 
        GRAD_SLOPE_DOT_NORMAL_FIRELINE = DZTDX(I,J)*(DPHIDY/MAG_F) + DZTDY(I,J)*(-DPHIDY/MAG_F)
   ELSE
        NORMAL_FIRELINE = 0._EB
        GRAD_SLOPE_DOT_NORMAL_FIRELINE = 0._EB
        XSF=0._EB
        YSF=0._EB
   ENDIF

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

       ROS_TMP = ROS_HEAD(I,J)
       
       !Mag of wind speed at midflame height must be in units of m/s here   
       UMF_DUM = UMF(I,J)/60.0_EB
       
       !Length to breadth ratio of ellipse based on effective UMF
       LB = 0.936_EB * EXP(0.2566_EB * UMF_DUM) + 0.461_EB * EXP(-0.1548_EB * UMF_DUM) - 0.397_EB 
       
       !Constraint LB max = 8 from Finney 2004
       LB = MAX(1.0_EB,MIN(LB,8.0_EB))
       
       LBD = SQRT(LB**2 - 1.0_EB)
       
       !Head to back ratio based on LB
       HB = (LB + LBD) / (LB - LBD)
       
       ! A_ELPS and B_ELPS notation is consistent with Farsite and Richards 
       B_ELPS =  0.5_EB * (ROS_TMP + ROS_TMP/HB)
       B_ELPS2 = B_ELPS**2
       A_ELPS =  B_ELPS / LB
       A_ELPS2=  A_ELPS**2
       C_ELPS =  B_ELPS - (ROS_TMP/HB)
  
       ! Denominator used in spread rate equation from Richards 1990 in final LS vs FS paper 
       AROS  = XSF*COS_THETA - YSF*SIN_THETA
       BROS  = XSF*SIN_THETA + YSF*COS_THETA
       DENOM = A_ELPS2*BROS**2 + B_ELPS2*AROS**2
           
       ! Finney's formulation
       !DENOM = B_ELPS2 * (XS * SIN_THETA - YS * COS_THETA)**2 - &
       !A_ELPS2 * (XS * COS_THETA + YS * SIN_THETA)**2
             
       IF (DENOM > 0._EB) THEN                 
        DENOM = 1._EB / SQRT(DENOM)        
       ELSE
        DENOM = 0._EB
       ENDIF
       
!  
!This is with A_ELPS2 and B_ELPS2 notation consistent with Finney and Richards and in final LS vs FS paper
        SR_X_LS(I,J) = DENOM * ( A_ELPS2*COS_THETA*BROS - B_ELPS2*SIN_THETA*AROS) + C_ELPS*SIN_THETA
        SR_Y_LS(I,J) = DENOM * (-A_ELPS2*SIN_THETA*BROS - B_ELPS2*COS_THETA*AROS) + C_ELPS*COS_THETA
        
        
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
   
     IF (MAG_ZT(I,J) > 0.0_EB) COS_THETA_SLOPE = GRAD_SLOPE_DOT_NORMAL_FIRELINE/MAG_ZT(I,J)
   
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
!if (abs(normal_fireline(1)) > 0._EB) print*,'rf,rh,rs',ros_flank(i,j),ros_head(i,j),ros_heads
!if (abs(normal_fireline(1)) > 0._EB) print*,'i,j',i,j
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
!if (abs(mag_sr) > 0._EB) print*,'mag_sr,nx,ny',mag_sr,normal_fireline(1),normal_fireline(2)
           SR_X_LS(I,J) = MAG_SR*NORMAL_FIRELINE(1) !spread rate components
           SR_Y_LS(I,J) = MAG_SR*NORMAL_FIRELINE(2) 
        !  MAG_SR_OUT(I,J) = MAG_SR
  
   ENDIF !Ellipse or McArthur Spread 
   
   DYN_SR_MAX_LS = MAX(DYN_SR_MAX_LS,MAG_SR) 

  ENDDO

ENDDO FLUX_ILOOP

END SUBROUTINE LEVEL_SET_PERIMETER_SPREAD_RATE 

!--------------------------------------------------------------------
!
SUBROUTINE LEVEL_SET_ADVECT_FLUX(NM)
!
! Use the spread rate [SR_X_LS,SR_Y_LS] to compute the limited scalar gradient
! and take dot product with spread rate vector to get advective flux

INTEGER, INTENT(IN) :: NM
INTEGER :: I,IM1,IM2,IP1,IP2,J,JM1,JM2,JP1,JP2
REAL(EB), DIMENSION(:) :: Z(4)
!REAL(EB), DIMENSION(:,:) :: FLUX_LS(0:IBP1,0:JBP1)
REAL(EB) :: DPHIDX,DPHIDY,F_EAST,F_WEST,F_NORTH,F_SOUTH
REAL(EB) :: PHIMAG

CALL POINT_TO_MESH(NM)

IF (RK2_PREDICTOR_LS) PHI0_LS = PHI_LS
IF (.NOT. RK2_PREDICTOR_LS) PHI0_LS = PHI1_LS

!if (nm == 2) then
!  print*,'predictor',rk2_predictor_ls
!  print '(A,i3,1x,2E13.5)','vege:advect ',nm,phi0_ls(0,12),phi0_ls(1,12)
!endif

ILOOP: DO I=1,NX_LS
 
 IM1=I-1!; IF (IM1<1) IM1=IM1+NX_LS
!IM2=I-2 ; IF (IM2<0) IM2=IM2+NX_LS
 IM2=I-2 ; IF (IM2<0) IM2=0

 IP1=I+1!; IF (IP1>NX_LS) IP1=IP1-NX_LS
!IP2=I+2 ; IF (IP2>NX_LS+1) IP2=IP2-NX_LS
 IP2=I+2 ; IF (IP2>NX_LS+1) IP2=NX_LS+1

 JLOOP: DO J = 1,NY_LS
   
   JM1=J-1!; IF (JM1<1) JM1=JM1+NY_LS
!  JM2=J-2 ; IF (JM2<0) JM2=JM2+NY_LS
   JM2=J-2 ; IF (JM2<0) JM2=0
   
   JP1=J+1!; IF (JP1>NY_LS) JP1=JP1-NY_LS
!  JP2=J+2 ; IF (JP2>NY_LS+1) JP2=JP2-NY_LS
   JP2=J+2 ; IF (JP2>NY_LS+1) JP2=NY_LS+1

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

!-- north face
   Z(1) = PHI0_LS(I,JM1)
   Z(2) = PHI0_LS(I,J)
   Z(3) = PHI0_LS(I,JP1)
   Z(4) = PHI0_LS(I,JP2)
   F_NORTH = SCALAR_FACE_VALUE_LS(SR_Y_LS(I,J),Z,LIMITER_LS)

!-- south face
   Z(1) = PHI0_LS(I,JM2)
   Z(2) = PHI0_LS(I,JM1)
   Z(3) = PHI0_LS(I,J)
   Z(4) = PHI0_LS(I,JP1)
   F_SOUTH = SCALAR_FACE_VALUE_LS(SR_Y_LS(I,J),Z,LIMITER_LS)
   
! IF (J<2 .OR. J>(NY_LS-2)) THEN  
!   
!      IF (J==1) THEN
!        !    north face
!!           Z(1) = PHI_MAX_LS
!            Z(1) = PHI0_LS(I,0)
!            Z(2) = PHI0_LS(I,J)
!            Z(3) = PHI0_LS(I,JP1)
!            Z(4) = PHI0_LS(I,JP2)
!            F_NORTH = SCALAR_FACE_VALUE_LS(SR_Y_LS(I,J),Z,LIMITER_LS)
!
!        !    south face
!!           Z(1) = PHI_MAX_LS
!!           Z(2) = PHI_MAX_LS
!            Z(1) = PHI0_LS(I,0)
!            Z(2) = PHI0_LS(I,0)
!            Z(3) = PHI0_LS(I,J)
!            Z(4) = PHI0_LS(I,JP1)
!            F_SOUTH = SCALAR_FACE_VALUE_LS(SR_Y_LS(I,J),Z,LIMITER_LS)
!!if (i==25 .and. nm==2) print '(A,1(i3),1x,2(E15.3))','vege:advect',nm,f_north,f_south
!
!      ELSEIF (j==2) THEN
!        !    north face
!            Z(1) = PHI0_LS(I,JM1)
!            Z(2) = PHI0_LS(I,J)
!            Z(3) = PHI0_LS(I,JP1)
!            Z(4) = PHI0_LS(I,JP2)
!            F_NORTH = SCALAR_FACE_VALUE_LS(SR_Y_LS(I,J),Z,LIMITER_LS)
!
!        !    south face
!!           Z(1) = PHI_MAX_LS
!            Z(1) = PHI0_LS(I,0)
!            Z(2) = PHI0_LS(I,JM1)
!            Z(3) = PHI0_LS(I,J)
!            Z(4) = PHI0_LS(I,JP1)
!            F_SOUTH = SCALAR_FACE_VALUE_LS(SR_Y_LS(I,J),Z,LIMITER_LS)
!   
!   
!      ELSEIF (J == NY_LS-1) THEN
!    !    north face
!            Z(1) = PHI0_LS(I,JM1)
!            Z(2) = PHI0_LS(I,J)
!            Z(3) = PHI0_LS(I,JP1)
!            Z(4) = PHI_MIN_LS
!            F_NORTH = SCALAR_FACE_VALUE_LS(SR_Y_LS(I,J),Z,LIMITER_LS)
!
!        !    south face
!            Z(1) = PHI0_LS(I,JM2)
!            Z(2) = PHI0_LS(I,JM1)
!            Z(3) = PHI0_LS(I,J)
!            Z(4) = PHI0_LS(I,JP1)
!            F_SOUTH = SCALAR_FACE_VALUE_LS(SR_Y_LS(I,J),Z,LIMITER_LS)
!
!      ELSEIF (J == NY_LS) THEN ! must be J == NY_LS
!        !    north face
!            Z(1) = PHI0_LS(I,JM1)
!            Z(2) = PHI0_LS(I,J)
!            Z(3) = PHI_MIN_LS
!            Z(4) = PHI_MIN_LS
!            F_NORTH = SCALAR_FACE_VALUE_LS(SR_Y_LS(I,J),Z,LIMITER_LS)
!
!        !    south face
!            Z(1) = PHI0_LS(I,JM2)
!            Z(2) = PHI0_LS(I,JM1)
!            Z(3) = PHI0_LS(I,J)
!            Z(4) = PHI_MIN_LS
!            F_SOUTH = SCALAR_FACE_VALUE_LS(SR_Y_LS(I,J),Z,LIMITER_LS)
!    
!      ENDIF
!  
!      ELSE
!
!    !    north face
!       Z(1) = PHI0_LS(I,JM1)
!       Z(2) = PHI0_LS(I,J)
!       Z(3) = PHI0_LS(I,JP1)
!       Z(4) = PHI0_LS(I,JP2)
!       F_NORTH = SCALAR_FACE_VALUE_LS(SR_Y_LS(I,J),Z,LIMITER_LS)
!
!    !    south face
!       Z(1) = PHI0_LS(I,JM2)
!       Z(2) = PHI0_LS(I,JM1)
!       Z(3) = PHI0_LS(I,J)
!       Z(4) = PHI0_LS(I,JP1)
!       F_SOUTH = SCALAR_FACE_VALUE_LS(SR_Y_LS(I,J),Z,LIMITER_LS)
!   
!   ENDIF !IF (J<2 .OR. J>(NY_LS-2) 
        
   DPHIDX = (F_EAST-F_WEST)* RDX(I) !IDX_LS
    
   DPHIDY = (F_NORTH-F_SOUTH)* RDY(J) !IDY_LS
   
   FLUX_LS(I,J) = SR_X_LS(I,J)*DPHIDX + SR_Y_LS(I,J)*DPHIDY
!if (j== 1 .and. i==25 .and. nm==2) print '(A,1(i3),1x,2(E15.3))','vege:advect--',nm,dphidy,flux_ls(i,j)
   
   PHIMAG          = SQRT(DPHIDX**2 + DPHIDY**2)
   MAG_SR_OUT(I,J) = 0.0_EB
   IF(PHIMAG > 0.0_EB) MAG_SR_OUT(I,J) = FLUX_LS(I,J)/PHIMAG
        
!  fx = (f_east-f_west)/dx
!  fy = (f_north-f_south)/dy
!       phi(i,j) = phi0(i,j) - dt*[Fx(i,j) Fy(i,j)]*[fx fy]

 ENDDO JLOOP

!FLUX_LS(:,1) = FLUX_LS(:,2)

ENDDO ILOOP

!print*,'veg advect_flux:maxflux    ',maxval(abs(flux_ls))
!print*,'veg advect_flux:max srx,sry',maxval(sr_x_ls),maxval(sr_y_ls)

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

!--------------------------------------------------------------------
SUBROUTINE LEVEL_SET_BC(NM)

! This finds the values of the level set function along the mesh boundaries through
! interpolation. Follows what's done for RHO in subroutine THERMAL_BC for INTERPOLATED_BC (in wall.f90)

INTEGER, INTENT(IN) :: NM
INTEGER :: II,IIG,IIO,IOR,IW,JJ,JJO,JJG,KK,KKG,NOM
!INTEGER :: KKO
REAL(EB) :: ARO,LSET_PHI_F,LSET_PHI_V
REAL(EB), POINTER, DIMENSION(:,:) :: OM_LSET_PHI =>NULL()
TYPE (WALL_TYPE),     POINTER :: WC =>NULL()
TYPE (EXTERNAL_WALL_TYPE),     POINTER :: EWC =>NULL()
TYPE (OMESH_TYPE),    POINTER :: OM =>NULL()
TYPE (MESH_TYPE),     POINTER :: MM =>NULL()

IF (EVACUATION_ONLY(NM)) RETURN
IF (.NOT. VEG_LEVEL_SET) RETURN

CALL POINT_TO_MESH(NM)

! Set ghost cell values (based on THERMAL_BC case INTERPOLATED_BC in wall.f90

WALL_CELL_LOOP_BC: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
   WC  => WALL(IW)
   IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE WALL_CELL_LOOP_BC
   II  = WC%ONE_D%II 
   JJ  = WC%ONE_D%JJ 
   KK  = WC%ONE_D%KK 
   IIG = WC%ONE_D%IIG
   JJG = WC%ONE_D%JJG
   KKG = WC%ONE_D%KKG
!  IF (KKG /= 1) CYCLE WALL_CELL_LOOP_BC
   IOR = WC%ONE_D%IOR
   IF (IOR == -3 .OR. IOR == 3) CYCLE WALL_CELL_LOOP_BC
   NOM = EWC%NOM
   IF (NOM == 0) CYCLE WALL_CELL_LOOP_BC
   OM  => OMESH(NOM)
   OM_LSET_PHI => OM%LSET_PHI
   MM  => MESHES(NOM)
   LSET_PHI_V = LSET_PHI(IIG,JJG)
   LSET_PHI_F = LSET_PHI_V  ! Initialize face value of LSET_PHI with LSET_PHI_V

!  DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
     DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
       DO IIO=EWC%KKO_MIN,EWC%KKO_MAX
         SELECT CASE(IOR)
         CASE( 1)
           ARO = MIN(1._EB , RDY(JJ)*MM%DY(JJO)) * 2._EB*DX(II)/(MM%DX(IIO)+DX(II))
         CASE(-1)
           ARO = MIN(1._EB , RDY(JJ)*MM%DY(JJO)) * 2._EB*DX(II)/(MM%DX(IIO)+DX(II))
         CASE( 2)
           ARO = MIN(1._EB , RDX(II)*MM%DX(IIO)) * 2._EB*DY(JJ)/(MM%DY(JJO)+DY(JJ))
         CASE(-2)
           ARO = MIN(1._EB , RDX(II)*MM%DX(IIO)) * 2._EB*DY(JJ)/(MM%DY(JJO)+DY(JJ))
!        CASE( 1)
!          ARO = MIN(1._EB , RDY(JJ)*RDZ(KK)*MM%DY(JJO)*MM%DZ(KKO)) * 2._EB*DX(II)/(MM%DX(IIO)+DX(II))
!        CASE(-1)
!          ARO = MIN(1._EB , RDY(JJ)*RDZ(KK)*MM%DY(JJO)*MM%DZ(KKO)) * 2._EB*DX(II)/(MM%DX(IIO)+DX(II))
!        CASE( 2)
!          ARO = MIN(1._EB , RDX(II)*RDZ(KK)*MM%DX(IIO)*MM%DZ(KKO)) * 2._EB*DY(JJ)/(MM%DY(JJO)+DY(JJ))
!        CASE(-2)
!          ARO = MIN(1._EB , RDX(II)*RDZ(KK)*MM%DX(IIO)*MM%DZ(KKO)) * 2._EB*DY(JJ)/(MM%DY(JJO)+DY(JJ))
!        CASE( 3)
!          ARO = MIN(1._EB , RDX(II)*RDY(JJ)*MM%DX(IIO)*MM%DY(JJO)) * 2._EB*DZ(KK)/(MM%DZ(KKO)+DZ(KK))
!        CASE(-3)
!          ARO = MIN(1._EB , RDX(II)*RDY(JJ)*MM%DX(IIO)*MM%DY(JJO)) * 2._EB*DZ(KK)/(MM%DZ(KKO)+DZ(KK))
         END SELECT
!          LSET_PHI_F =  LSET_PHI_F + 0.5_EB*ARO*(OM_LSET_PHI(IIO,JJO,KKO)-LSET_PHI_V)
           LSET_PHI_F =  LSET_PHI_F + 0.5_EB*ARO*(OM_LSET_PHI(IIO,JJO)-LSET_PHI_V)
       ENDDO
     ENDDO
!  ENDDO
   LSET_PHI(II,JJ) = MIN( 1._EB , MAX( -1._EB , 2._EB*LSET_PHI_F-LSET_PHI_V ))

!if (jj == 0 .and. ii==12 .and. nm == 3) then 
!    print*,'==j0 bc'
!    print '(A,1x,1(i3),1x,1(E12.5))','vege:j0 bc',nm,lset_phi(ii,jj,1)
!    print '(A,1x,1(i3),1x,1(E12.5))','vege:lsbc',nom,om_lset_phi(ii,jbar,1)
!endif

!if (ii == 0 .and. jj==12 .and. nm==2) then
!  print*,'--ibar bc'
!  print '(A,1x,1(i3),1x,1(E12.5))','vege:ibp1 bc',nm,lset_phi(ii,jj,1)
!  print '(A,1x,2(i3),1x,1(E12.5))','vege:',nom,wc%nom_ib(1),om_lset_phi(wc%nom_ib(1),jj,1)
!endif

!if (nm==2) then
! print*,'\/\/\/\ mesh 1 & bounding i,j cells'
! print '(A,1x,3i3)','nm, i,j cells',nm,ii,jj
! print '(A,1x,3i3)','nom, i cells',nom,wc%nom_ib(1),wc%nom_ib(4)
! print '(A,1x,3i3)','nom, j cells',nom,wc%nom_ib(2),wc%nom_ib(5)
!endif
!f (jj == jbp1) lset_phi(ii,jj,1) = om_lset_phi(ii,1,1)
!f (jj ==    0) om_lset_phi(ii,jbar,1) = om_lset_phi(ii,0,1)

ENDDO WALL_CELL_LOOP_BC

END SUBROUTINE LEVEL_SET_BC

!--------------------------------------------------------------------
SUBROUTINE LEVEL_SET_DT(DT,NM)
!Variable Time Step for simulations uncoupled from the CFD computation

INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(INOUT) :: DT
REAL(EB) :: PHI_CHECK_LS
!REAL(EB) :: DT_CHECK_LS

!IF (EVACUATION_ONLY(NM)) RETURN
IF (.NOT. VEG_LEVEL_SET) RETURN
IF (VEG_LEVEL_SET_COUPLED) RETURN

!DO NM=1,NMESHES

  CALL POINT_TO_MESH(NM)

  PHI_CHECK_LS = MAXVAL(ABS(PHI_LS - PHI_TEMP_LS)) !Find max change in phi

!print '(A,1x,E13.5,1x,i3)','*************ls_dt,phi check,nm',phi_check_ls, nm !maxval(phi_Ls),nm
!print*,'ls_dt lset_ignition, nm',lset_ignition,nm

 
! IF (LSET_IGNITION) THEN
CHANGE_TIME_STEP_INDEX(NM) = 0
IF (PHI_CHECK_LS > 0.0_EB) THEN
! If any phi values change by more than 0.5, or all change less
! than 0.1 (both values are arbitrary), during one time step,
! then adjust time step accordingly.

IF (PHI_CHECK_LS < 0.5_EB) THEN
!  DT_NEXT = DT
!  IF (PHI_CHECK_LS <= 0.1_EB) DT_NEXT = MIN(1.1_EB*DT,MIN(DX(1),DY(1))/DYN_SR_MAX_LS)
ELSE
   DT = 0.5_EB*DT
   CHANGE_TIME_STEP_INDEX(NM) = -1
ENDIF

!    CHANGE_TIME_STEP(NM) = .FALSE.
!    IF (PHI_CHECK_LS > 0.5_EB) THEN
!      ! Cut time step in half
!      DT_COEF_LS = 0.5_EB * DT_COEF_LS 
!      DT_CHECK_LS = DT_COEF_LS * MIN(DX(1),DY(1))/DYN_SR_MAX_LS
!      DT_CHECK_LS = MIN(DT_CHECK_LS,100._EB)
!!     DT_LS_UNCOUPLED = DT_CHECK_LS
!      print '(A,1x,E13.5,1x,A,i3)',"ls_dt Halving time step to, dt=  ",dt_check_ls,' mesh ',nm
!!     MESHES(NM)%DT_NEXT = DT_CHECK_LS
!      DT = DT_CHECK_LS
!      CHANGE_TIME_STEP(NM) = .TRUE.
!      RETURN
!    ELSE
!      DT_NEXT = DT
!    ENDIF
! 
!    ! Increase time step by 1/4 if changes are small
!!   IF (PHI_CHECK_LS > 0.01_EB .AND. PHI_CHECK_LS < 0.1_EB) DT_COEF_LS = DT_COEF_LS * 1.25_EB
!     
!    ! Dynamic Spread Rate Max
!      print '(A,1x,E13.5,1x,A,i3)',"ls_dt dyn_sr_max =  ",dyn_sr_max,' mesh ',nm
!    DYN_SR_MAX_LS = MAX(DYN_SR_MAX_LS,0.01_EB) ! dyn_sr_max must be g.t. zero
!!   DT_CHECK_LS = DT_COEF_LS * MIN(DX(1),DY(1))/DYN_SR_MAX_LS
!    DT_CHECK_LS = 0.25_EB*DT_COEF_LS * MIN(DX(1),DY(1))/DYN_SR_MAX_LS
!    DT_CHECK_LS = MIN(DT_CHECK_LS,100._EB)
!    IF (DT_CHECK_LS > DT) THEN
!      print '(A,1x,E13.5,1x,A,i3)',"ls_dt increasing time step to, dt=  ",dt_check_ls,' mesh ',nm
!      DT_NEXT = DT_CHECK_LS 
!    ENDIF
!
!!   IF (MAXVAL(PHI_LS) == -1.0_EB) THEN
!!     DT_NEXT = 10._EB*DT
!!     print '(A,1x,E13.5,1x,A,i3)',"ls_dt dphi=0,increasing time step to, dt=  ",dt_next,' mesh ',nm
!!   ENDIF
!   
!!   MESHES(NM)%DT_NEXT = DT_CHECK_LS
!
!!   DT_LS_UNCOUPLED = MAX(DT_LS,DT_LS_UNCOUPLED)
     
ENDIF

!ENDDO

!print '(A,1x,1E13.5)',"ls_dt:meshes(:)%dt=  ",meshes(1:nmeshes)%dt

END SUBROUTINE LEVEL_SET_DT


!************************************************************************************************
SUBROUTINE READ_BRNR
!************************************************************************************************
! Read in, from the run input file, the names of the file(s) containing HRRPUA in each mesh along 
! the bottom of the domain.
!
USE OUTPUT_DATA

INTEGER :: I,IOS,BURNER_MESH_NUMBER
CHARACTER(256) :: BRNRFILE='null'
NAMELIST /BRNR/ BURNER_MESH_NUMBER,BRNRFILE

N_BRNR=0
REWIND(LU_INPUT) ; INPUT_FILE_LINE_NUMBER = 0
COUNT_BRNR_LOOP: DO
   CALL CHECKREAD('BRNR',LU_INPUT,IOS)
   IF (IOS==1) EXIT COUNT_BRNR_LOOP
   READ(LU_INPUT,NML=BRNR,END=16,ERR=17,IOSTAT=IOS)
   N_BRNR=N_BRNR+1
   16 IF (IOS>0) THEN ; CALL SHUTDOWN('ERROR: problem with BRNR line') ; RETURN ; ENDIF
ENDDO COUNT_BRNR_LOOP
17 REWIND(LU_INPUT) ; INPUT_FILE_LINE_NUMBER = 0

IF (N_BRNR==0) RETURN

ALLOCATE(BRNRINFO(N_BRNR),STAT=IZERO) !array for burner data file names
CALL ChkMemErr('VEGE','BRNRINFO',IZERO)
ALLOCATE(BURNER_FILE(NMESHES),STAT=IZERO) ; BURNER_FILE = -99 !array for mesh number of burner data file
CALL ChkMemErr('VEGE','BURNER_FILE',IZERO)

IF (N_BRNR > NMESHES) THEN
  CALL SHUTDOWN('Problem with BRNR lines: N_BRNR > NMESHES') 
  RETURN
ENDIF
IF (.NOT. VEG_LEVEL_SET_BURNERS_FOR_FIRELINE) THEN
  CALL SHUTDOWN('Problem with BRNR lines: use of burners requires VEG_LEVEL_SET_BURNERS_FOR_FIRELINE=.TRUE.') 
  RETURN
ENDIF

READ_BRNR_LOOP: DO I=1,N_BRNR

   CALL CHECKREAD('BRNR',LU_INPUT,IOS)
   IF (IOS==1) EXIT READ_BRNR_LOOP

   ! Read the BRNR line

   READ(LU_INPUT,BRNR,END=37)

   BRNRINFO(I)%BRNRFILE = TRIM(BRNRFILE)
   BURNER_FILE(BURNER_MESH_NUMBER) = I

!print*,'read_brnr, nm, brnrfile',burner_mesh_number,brnrinfo(i)%brnrfile

ENDDO READ_BRNR_LOOP
37 REWIND(LU_INPUT) ; INPUT_FILE_LINE_NUMBER = 0

END SUBROUTINE READ_BRNR

END MODULE VEGE
