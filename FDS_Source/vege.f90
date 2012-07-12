MODULE VEGE
 
USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_POINTERS
USE TRAN
USE PART
USE MEMORY_FUNCTIONS, ONLY:CHKMEMERR
USE TYPES, ONLY: LAGRANGIAN_PARTICLE_TYPE, LAGRANGIAN_PARTICLE_CLASS_TYPE! WALL_TYPE,SURFACE_TYPE 
IMPLICIT NONE
PRIVATE
PUBLIC LEVEL_SET_FIRESPREAD, GET_REV_vege, BNDRY_VEG_MASS_ENERGY_TRANSFER
TYPE (LAGRANGIAN_PARTICLE_TYPE), POINTER :: LP=>NULL()
TYPE (LAGRANGIAN_PARTICLE_CLASS_TYPE), POINTER :: LPC=>NULL()
!TYPE (WALL_TYPE), POINTER :: WC
!TYPE (SURFACE_TYPE), POINTER :: SF 
CHARACTER(255), PARAMETER :: vegeid='$Id$'
CHARACTER(255), PARAMETER :: vegerev='$Revision$'
CHARACTER(255), PARAMETER :: vegedate='$Date$'
LOGICAL, ALLOCATABLE, DIMENSION(:,:,:) :: VEG_PRESENT_FLAG,CELL_TAKEN_FLAG
INTEGER :: IZERO,NLP_VEG_FUEL,NCONE_TREE,NXB,NYB
REAL(EB) :: RCELL,R_TREE,XCELL,XI,YJ,YCELL,ZCELL,ZK
 
CONTAINS
 

SUBROUTINE BNDRY_VEG_MASS_ENERGY_TRANSFER(T,NM)
!
! Issues:
! 1. Are SF%VEG_FUEL_FLUX_L and SF%VEG_MOIST_FLUX_L needed in linear degradation model?
REAL(EB) :: DT_BC,RDT_BC,T
INTEGER, INTENT(IN) :: NM
INTEGER  ::  IW
INTEGER  ::  I,IIG,JJG,KKG
REAL(EB) :: CP_MOIST_AND_VEG,DZVEG_L,ETAVEG_H,H_CONV_FDS_WALL,H_CONV_L, &
            KAPPA_VEG,LAMBDA_AIR,QRADM_INC,QRADP_INC, &
            TMP_BOIL,TMPG_A,TMP_G,DTMP_L,DTMP_FDS_WALL !,RE_VEG_PART,U2,V2
INTEGER  IIVEG_L,IVEG_L,J,LBURN,NVEG_L,I_FUEL
!REAL(EB), ALLOCATABLE, DIMENSION(:) :: VEG_DIV_QRNET_EMISS,VEG_DIV_QRNET_INC,
!         VEG_QRNET_EMISS,VEG_QRNET_INC,VEG_QRM_EMISS,VEG_QRP_EMISS, VEG_QRM_INC,VEG_QRP_INC
REAL(EB) :: VEG_DIV_QRNET_EMISS(20),VEG_DIV_QRNET_INC(20),VEG_QRNET_EMISS(0:20),VEG_QRNET_INC(0:20), &
            VEG_QRM_EMISS(0:20),VEG_QRP_EMISS(0:20), VEG_QRM_INC(0:20),VEG_QRP_INC(0:20)
REAL(EB) :: A_H2O_VEG,E_H2O_VEG,A_PYR_VEG,E_PYR_VEG,L_PYR_VEG
REAL(EB) :: A_CHAR_VEG,E_CHAR_VEG,BETA_CHAR_VEG,NU_CHAR_VEG,NU_ASH_VEG,NU_O2_CHAR_VEG,H_CHAR_OXID
REAL(EB) :: CP_H2O,CP_VEG,DTMP_VEG,H_VAP_H2O,TMP_VEG,TMP_VEG_NEW
REAL(EB) :: CHAR_FCTR,MPA_MOIST,MPA_MOIST_LOSS,MPA_MOIST_LOSS_MAX,MPA_MOIST_MIN,MPA_VEG,MPA_VEG_MIN, & 
            MPA_VOLIT,MPA_VOLIT_MAX
REAL(EB) :: DETA_VEG,ETA_H,ETAFM_VEG,ETAFP_VEG
REAL(EB) :: QCONF_FDS_WALL,QCONF_L,Q_FOR_DRYING,Q_VEG_MOIST,Q_VEG_VOLIT,QNET_VEG,Q_FOR_VOLIT,Q_VOLIT,Q_UPTO_VOLIT
LOGICAL  :: VEG_DEGRADATION_ARRHENIUS,VEG_DEGRADATION_LINEAR
logical  :: fuel_elem_degrad,fds4_degrad

TYPE (WALL_TYPE),    POINTER :: WC =>NULL()
TYPE (SURFACE_TYPE), POINTER :: SF =>NULL()

CALL POINT_TO_MESH(NM)

TMP_BOIL  = 373._EB
CP_H2O    = 4190._EB !J/kg/K specific heat of water
H_VAP_H2O = 2259._EB*1000._EB !J/kg/K heat of vaporization of water
DT_BC     = T - VEG_CLOCK_BC
RDT_BC    = 1.0_EB/DT_BC

IF (N_REACTIONS>0) I_FUEL = REACTION(1)%FUEL_SMIX_INDEX

! Thermal degradation approach parameters
! VEG_DEGRADATION_LINEAR    = .TRUE.
! VEG_DEGRADATION_ARRHENIUS = .FALSE.
  fuel_elem_degrad = .true.
  fds4_degrad      = .false.
!
! Loop through vegetation wall cells and burn
!
VEG_WALL_CELL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
  WC  => WALL(IW)
  IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE VEG_WALL_CELL_LOOP

    SF  => SURFACE(WC%SURF_INDEX)
!
  IF (.NOT. SF%VEGETATION) CYCLE VEG_WALL_CELL_LOOP

  VEG_DEGRADATION_LINEAR    = SF%VEG_LINEAR_DEGRAD
  VEG_DEGRADATION_ARRHENIUS = SF%VEG_ARRHENIUS_DEGRAD


  IIG = WC%ONE_D%IIG
  JJG = WC%ONE_D%JJG
  KKG = WC%ONE_D%KKG
  TMP_G = TMP(IIG,JJG,KKG)
  CHAR_FCTR = 1._EB - SF%VEG_CHARFRAC
  IF(SF%VEG_NO_BURN) WC%VEG_HEIGHT = SF%VEG_HEIGHT
  VEG_DRAG(IIG,JJG) = SF%VEG_DRAG_INI*(SF%VEG_CHARFRAC + CHAR_FCTR*WC%VEG_HEIGHT/SF%VEG_HEIGHT)

  IF(SF%VEG_NO_BURN) CYCLE VEG_WALL_CELL_LOOP

! Initialize quantities
  Q_VEG_MOIST     = 0.0_EB
  Q_VEG_VOLIT     = 0.0_EB
  Q_UPTO_VOLIT    = 0.0_EB
  Q_VOLIT         = 0.0_EB
  MPA_MOIST_LOSS  = 0.0_EB
  MPA_VOLIT       = 0.0_EB
  SF%VEG_DIVQNET_L = 0.0_EB
  SF%VEG_MOIST_FLUX_L = 0.0_EB
  SF%VEG_FUEL_FLUX_L  = 0.0_EB
  WC%ONE_D%MASSFLUX(I_FUEL) = 0.0_EB 
  WC%ONE_D%QCONF           = 0.0_EB
  IF (I_WATER /= 0) WC%ONE_D%MASSFLUX(I_WATER) = 0.0_EB

! Vegetation variables and minimum bounds
  NVEG_L = SF%NVEG_L
  LBURN  = 0
  MPA_VEG_MIN   = SF%VEG_CHARFRAC*SF%VEG_LOAD / REAL(NVEG_L,EB) !kg/m^2
  MPA_MOIST_MIN = 0.001_EB*SF%VEG_MOISTURE*SF%VEG_LOAD/REAL(NVEG_L,EB) !ks/m^2
  IF (ABS(SF%VEG_MOISTURE) <= ZERO_P) MPA_MOIST_MIN = 1._EB
  DZVEG_L   = SF%VEG_HEIGHT/REAL(NVEG_L,EB)
  KAPPA_VEG = SF%VEG_KAPPA
  DETA_VEG  = DZVEG_L*KAPPA_VEG

! Find top of vegetation which burns downward from the top
  DO IVEG_L = 1,NVEG_L 
    IF(WC%VEG_FUELMASS_L(IVEG_L) <= MPA_VEG_MIN) LBURN = IVEG_L
  ENDDO
  WC%VEG_HEIGHT = REAL(NVEG_L-LBURN,EB)*DZVEG_L
  LBURN = 0._EB !keep charred veg
! MPA_VOLIT_MAX      = DT_BC*0.1_EB/REAL(NVEG_L-LBURN,EB) !upper bound F19 AU exp U = 5 m/s
! MPA_VOLIT_MAX      = DT_BC*0.8_EB/REAL(NVEG_L-LBURN,EB) !upper bound F19 AU exp d=20 m U = 20 m/s
  !FIRELINE_MLR_MAX = w*R*(1-ChiChar)
  MPA_VOLIT_MAX      = SF%FIRELINE_MLR_MAX*DT_BC*DX(IIG)/REAL(NVEG_L-LBURN,EB) 
  MPA_MOIST_LOSS_MAX = MPA_VOLIT_MAX
! MPA_MOIST_LOSS_MAX = DT_BC*0.05_EB/REAL(NVEG_L-LBURN,EB)
! MPA_VOLIT_MAX      = DT_BC*0.05_EB/REAL(NVEG_L-LBURN,EB)
! MPA_MOIST_LOSS_MAX = 9999999._EB
! MPA_VOLIT_MAX      = 9999999._EB

! Factors for computing divergence of incident and self emission radiant fluxes
! in vegetation fuel bed. These need to be recomputed as the height of the
! vegetation surface layer decreases with burning

! Factors for computing decay of +/- incident fluxes
  SF%VEG_FINCM_RADFCT_L(:) =  0.0_EB
  SF%VEG_FINCP_RADFCT_L(:) =  0.0_EB
! ETA_H = KAPPA_VEG*WC%VEG_HEIGHT
  ETA_H = KAPPA_VEG*REAL(NVEG_L-LBURN,EB)*DZVEG_L
  DO IVEG_L = 0,SF%NVEG_L - LBURN
    ETAFM_VEG = IVEG_L*DETA_VEG
    ETAFP_VEG = ETA_H - ETAFM_VEG
    SF%VEG_FINCM_RADFCT_L(IVEG_L) = EXP(-ETAFM_VEG)
    SF%VEG_FINCP_RADFCT_L(IVEG_L) = EXP(-ETAFP_VEG)
  ENDDO

!  Integrand for computing +/- self emission fluxes
  SF%VEG_SEMISSP_RADFCT_L(:,:) = 0.0_EB
  SF%VEG_SEMISSM_RADFCT_L(:,:) = 0.0_EB
! q+
  DO IIVEG_L = 0,SF%NVEG_L-LBURN !veg grid coordinate
    DO IVEG_L = IIVEG_L,SF%NVEG_L-1-LBURN !integrand index
!    ETAG_VEG = IIVEG_L*DETA_VEG
!    ETAI_VEG =  IVEG_L*DETA_VEG
!    SF%VEG_SEMISSP_RADFCT_L(IVEG_L,IIVEG_L) = EXP(-(ETAI_VEG-ETAG_VEG))
     ETAFM_VEG = (IVEG_L-IIVEG_L)*DETA_VEG
     ETAFP_VEG = ETAFM_VEG + DETA_VEG
     SF%VEG_SEMISSP_RADFCT_L(IVEG_L,IIVEG_L) = EXP(-ETAFM_VEG) - EXP(-ETAFP_VEG)
    ENDDO
  ENDDO
! q-
  DO IIVEG_L = 0,SF%NVEG_L-LBURN
    DO IVEG_L = 1,IIVEG_L
!    ETAG_VEG = IIVEG_L*DETA_VEG
!    ETAI_VEG =  IVEG_L*DETA_VEG
!    SF%VEG_SEMISSM_RADFCT_L(IVEG_L,IIVEG_L) = EXP(-(ETAG_VEG-ETAI_VEG))
     ETAFM_VEG = (IIVEG_L-IVEG_L)*DETA_VEG
     ETAFP_VEG = ETAFM_VEG + DETA_VEG
     SF%VEG_SEMISSM_RADFCT_L(IVEG_L,IIVEG_L) = EXP(-ETAFM_VEG) - EXP(-ETAFP_VEG)
    ENDDO
  ENDDO
!
! compute CONVECTIVE HEAT FLUX on vegetation
! cylinder heat transfer coefficient, hc, from Albini CST, assumes
! lambda ~ rho*cp*T^1.5/p where cp (of air) is assumed to be 
! independent of temperature. Flux is from Morvan and Dupuy assuming
! constant physical properties and integrating vertically over fuel
! bed to get a factor of h multiplying their qc'''
! DTMP*BETA*sigma*h*hc*(T-Ts)
! hc = 0.350*(sigma/4)*lambda in Albini CST 1985 assumes quiescent air
! hc = 0.683*(sigma/4)*lambda*Re^0.466 ; Re=|u|r/nu, r=2/sigma
!      used by Porterie, cylinders in air flow
! lambda = lambda0*(rho/rho0)(T/T0)^a; a=1.5 below
  TMPG_A     = (TMP_G*0.0033_EB)**1.5
  LAMBDA_AIR = 0.026_EB*RHO(IIG,JJG,KKG)*0.861_EB*TMPG_A
!Albini assumes quiescent air
! H_CONV_L = 0.35*LAMBDA_AIR*SF%VEG_SVRATIO*0.25
!Holman "Heat Transfer",5th Edition, McGraw-Hill, 1981 p.285 
!assumes vertical cylinder laminar air flow
! H_CONV_L = 1.42*(DTMP/VEG_HEIGHT_S(SURF_INDEX))**0.25 !W/m^2/C
!Porterie allow for air flow
! U2 = 0.25*(U(IIG,JJG,KKG)+U(IIG-1,JJG,KKG))**2
! V2 = 0.25*(V(IIG,JJG,KKG)+V(IIG,JJG-1,KKG))**2
! RE_VEG_PART = SQRT(U2 + V2)*2./SF%VEG_SVRATIO/TMPG_A/15.11E-6
! H_CONV_L = 0.5*LAMBDA_AIR*0.683*RE_VEG_PART**0.466*0.5*SF%VEG_SVRATIO
!
  DTMP_FDS_WALL   = TMP_G - WALL(IW)%ONE_D%TMP_F
  H_CONV_FDS_WALL = 1.42_EB*(ABS(DTMP_FDS_WALL)/DZVEG_L)**0.25
  QCONF_FDS_WALL  = H_CONV_FDS_WALL*DTMP_FDS_WALL
! QCONF(IW)       = QCONF_FDS_WALL !W/m^2
! print*,'dtmp_fds_wall,qconf',dtmp_fds_wall,qconf(iw)
! print*,'tmp_g,tmp_f(iw)',tmp_g,tmp_f(iw)
! SF%VEG_DIVQNET_L(1) = SF%VEG_PACKING*SF%VEG_SVRATIO*QCONF_L*DZVEG_L !W/m^2

  DO I=1,NVEG_L-LBURN
    DTMP_L = TMP_G - WC%VEG_TMP_L(I+LBURN)
!   DTMP_L = 0.5_EB*(TMP_G - WC%VEG_TMP_L(I+LBURN))
!Holman see ref above (needs DTMP so its computation is done here)
    H_CONV_L = 1.42_EB*(ABS(DTMP_L)/DZVEG_L)**0.25
    QCONF_L  = H_CONV_L*DTMP_L
!
    SF%VEG_DIVQNET_L(I) = SF%VEG_PACKING*SF%VEG_SVRATIO*QCONF_L*DZVEG_L !W/m^2
!   QCONF(IW) = QCONF(IW) + QCONF_L !W/m^2
  ENDDO
  WALL(IW)%ONE_D%QCONF = SUM(SF%VEG_DIVQNET_L) !*RDN(IW)*WC%VEG_HEIGHT
! qconf(iw) = 0.0_EB
!
! Compute +/- radiation fluxes and their divergence due to self emission within vegetation
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
    QRADM_INC = WALL(IW)%ONE_D%QRADIN/WALL(IW)%ONE_D%EMISSIVITY !sigma*Ta^4 + flame
!   QRADM_INC = WALL(IW)%ONE_D%QRADIN/WALL(IW)%ONE_D%EMISSIVITY + SIGMA*WALL(IW)%ONE_D%TMP_F**4 ! as done in FDS4
!   print*,'vege: WALL(IW)%ONE_D%QRADIN',WALL(IW)%ONE_D%QRADIN
    ETAVEG_H  = (NVEG_L - LBURN)*DETA_VEG
    !this QRADP_INC ensures zero net radiant fluxes at bottom of vegetation
    IF(SF%VEG_GROUND_ZERO_RAD) QRADP_INC = QRADM_INC*SF%VEG_FINCM_RADFCT_L(NVEG_L-LBURN) + VEG_QRM_EMISS(NVEG_L-LBURN)
    !this QRADP_INC assumes the ground stays at ambient temperature
    IF(.NOT. SF%VEG_GROUND_ZERO_RAD) QRADP_INC = SIGMA*SF%VEG_GROUND_TEMP**4
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
    SF%VEG_DIVQNET_L(I)= SF%VEG_DIVQNET_L(I) - (VEG_DIV_QRNET_INC(I) + VEG_DIV_QRNET_EMISS(I))
!   SF%VEG_DIVQNET_L(I)= SF%VEG_DIVQNET_L(I) - VEG_DIV_QRNET_INC(I)
  ENDDO
!
!
!      ************** Boundary Fuel Non-Arrehnius (Linear in temp) Degradation model *************************
! Drying occurs if qnet > 0 with Tveg held at 100 c
! Pyrolysis occurs according to Morvan & Dupuy empirical formula. Linear
! temperature dependence with qnet factor
!

  IF_VEG_DEGRADATION_LINEAR: IF (VEG_DEGRADATION_LINEAR) THEN

  if_fds4_degrad: if(fds4_degrad) then
!
! compute mass flux of H20 vapor or fuel gas
!
!     VEG_MOIST_FLUX_L   = 0.0
!     VEG_FUEL_FLUX_L    = 0.0
!     WC%ONE_D%MASSFLUX(IFUEL) = 0.0
      DO IVEG_L = LBURN+1,NVEG_L

       IF(SF%VEG_DIVQNET_L(IVEG_L-LBURN) > 0._EB) THEN 

!                                 -- boiling 
        IF(WC%VEG_TMP_L(IVEG_L)>=TMP_BOIL .AND. WC%VEG_MOISTMASS_L(IVEG_L)>MPA_MOIST_MIN) THEN
          SF%VEG_MOIST_FLUX_L(IVEG_L) = SF%VEG_DIVQNET_L(IVEG_L-LBURN)/H_VAP_H2O
          WC%VEG_MOISTMASS_L(IVEG_L) = WC%VEG_MOISTMASS_L(IVEG_L) - DT_BC*SF%VEG_MOIST_FLUX_L(IVEG_L)
!         print*,'&& layer, water mass',
!     .  iveg_l,veg_water_mass_per_area_l(iw,iveg_l)
        ENDIF
!                                 -- pyrolysis multiple layers
        IF (WC%VEG_MOISTMASS_L(IVEG_L)<=MPA_MOIST_MIN .AND. WC%VEG_FUELMASS_L(IVEG_L)>MPA_VEG_MIN) THEN

        IF(WC%VEG_TMP_L(IVEG_L)>= 400._EB .AND. WC%VEG_TMP_L(IVEG_L)<= 500._EB)  &
         SF%VEG_FUEL_FLUX_L(IVEG_L) = SF%VEG_DIVQNET_L(IVEG_L-LBURN)*0.0000025*(WC%VEG_TMP_L(IVEG_L)-400.)*0.01 
!
        WC%VEG_FUELMASS_L(IVEG_L) = WC%VEG_FUELMASS_L(IVEG_L) - DT_BC*SF%VEG_FUEL_FLUX_L(IVEG_L)
        WC%ONE_D%MASSFLUX(I_FUEL)= SF%VEG_FUEL_FLUX_L(IVEG_L) + WC%ONE_D%MASSFLUX(I_FUEL)
        ENDIF !pyrolysis models
        WC%VEG_FUELMASS_L(IVEG_L) = MAX(WC%VEG_FUELMASS_L(IVEG_L),MPA_VEG_MIN)
       ENDIF !qnetflux_l > 0
       ENDDO !boil off and pyrolysis
!
! Compute temperature of vegetation
!
      VEG_TEMP_LOOP: DO IVEG_L = LBURN+1,NVEG_L
!
      IF (WC%VEG_MOISTMASS_L(IVEG_L) < MPA_MOIST_MIN) WC%VEG_MOISTMASS_L(IVEG_L) = 0.0
!
      CP_VEG = (0.01 + 0.0037*WC%VEG_TMP_L(IVEG_L))*1000. !W/kg/K
      CP_MOIST_AND_VEG = CP_H2O*WC%VEG_MOISTMASS_L(IVEG_L) + CP_VEG*WC%VEG_FUELMASS_L(IVEG_L)
      WC%VEG_TMP_L(IVEG_L) = WC%VEG_TMP_L(IVEG_L) + DT_BC*( SF%VEG_DIVQNET_L(IVEG_L-LBURN)  &
                             - SF%VEG_MOIST_FLUX_L(IVEG_L)*H_VAP_H2O &
                             - SF%VEG_FUEL_FLUX_L(IVEG_L)*416000. )/CP_MOIST_AND_VEG

!  Set veg. temp to boiling temp if appropriate
      IF (WC%VEG_MOISTMASS_L(IVEG_L)>=MPA_MOIST_MIN .AND. WC%VEG_TMP_L(IVEG_L)>=TMP_BOIL)  &
           WC%VEG_TMP_L(IVEG_L) = TMP_BOIL

!  Set veg. temp to pyroysis temp if appropriate
!     IF (SF%VEG_FUEL_FLUX_L(IVEG_L)>0._EB .AND. WC%VEG_TMP_L(IVEG_L)>500._EB) WC%VEG_TMP_L(IVEG_L) = 500.
      IF (WC%VEG_FUELMASS_L(IVEG_L)>MPA_VEG_MIN .AND. WC%VEG_TMP_L(IVEG_L)>500._EB) WC%VEG_TMP_L(IVEG_L) = 500._EB

      ENDDO VEG_TEMP_LOOP

    WC%VEG_TMP_L(LBURN) = WC%VEG_TMP_L(LBURN+1)

  endif if_fds4_degrad


  if_fuel_elem_degrad: if (fuel_elem_degrad) then

    LAYER_LOOP1: DO IVEG_L = LBURN+1,NVEG_L
!
! Compute temperature of vegetation
!
      MPA_VEG   = WC%VEG_FUELMASS_L(IVEG_L)
      MPA_MOIST = WC%VEG_MOISTMASS_L(IVEG_L)
      TMP_VEG   = WC%VEG_TMP_L(IVEG_L)
      QNET_VEG  = SF%VEG_DIVQNET_L(IVEG_L-LBURN)
      CP_VEG = (0.01_EB + 0.0037_EB*TMP_VEG)*1000._EB !J/kg/K
      CP_MOIST_AND_VEG = CP_H2O*MPA_MOIST +  CP_VEG*MPA_VEG
      DTMP_VEG = DT_BC*QNET_VEG/CP_MOIST_AND_VEG
      TMP_VEG_NEW = TMP_VEG + DTMP_VEG 
!     IF(TMP_VEG_NEW >= 800._EB .AND. MPA_VEG <= MPA_VEG_MIN) TMP_VEG_NEW = 800._EB

      IF_DIVQ_L_GE_0: IF(QNET_VEG > 0._EB) THEN 

! -- drying of veg layer
      IF(MPA_MOIST > MPA_MOIST_MIN .AND. TMP_VEG_NEW >= TMP_BOIL) THEN
        Q_FOR_DRYING   = (TMP_VEG_NEW - TMP_BOIL)/DTMP_VEG * QNET_VEG
        MPA_MOIST_LOSS = MIN(DT_BC*Q_FOR_DRYING/H_VAP_H2O,MPA_MOIST-MPA_MOIST_MIN)
        MPA_MOIST_LOSS = MIN(MPA_MOIST_LOSS,MPA_MOIST_LOSS_MAX) !use specified max
        TMP_VEG_NEW    = TMP_BOIL
        WC%VEG_MOISTMASS_L(IVEG_L) = MPA_MOIST - MPA_MOIST_LOSS !kg/m^2
        IF( WC%VEG_MOISTMASS_L(IVEG_L) <= MPA_MOIST_MIN ) WC%VEG_MOISTMASS_L(IVEG_L) = 0.0_EB
        IF (I_WATER /= 0) WC%ONE_D%MASSFLUX(I_WATER) = WC%ONE_D%MASSFLUX(I_WATER) + RDT_BC*MPA_MOIST_LOSS
!       WC%VEG_TMP_L(IVEG_L) = TMP_VEG_NEW
      ENDIF

! -- pyrolysis multiple layers
      IF_VOLITIZATION: IF (MPA_MOIST <= MPA_MOIST_MIN) THEN

        IF(TMP_VEG_NEW >= 400._EB .AND. MPA_VEG > MPA_VEG_MIN) THEN
          Q_UPTO_VOLIT = CP_VEG*MPA_VEG*(400._EB-TMP_VEG)
          Q_FOR_VOLIT  = DT_BC*QNET_VEG - Q_UPTO_VOLIT
          Q_VOLIT      = Q_FOR_VOLIT*0.01_EB*(TMP_VEG-400._EB)
          MPA_VOLIT    = CHAR_FCTR*Q_VOLIT*0.00000239_EB
          MPA_VOLIT    = MAX(MPA_VOLIT,0._EB)
          MPA_VOLIT    = MIN(MPA_VOLIT,MPA_VOLIT_MAX)
          MPA_VOLIT    = MIN(MPA_VOLIT,MPA_VEG-MPA_VEG_MIN)
          MPA_VEG      = MPA_VEG - MPA_VOLIT
          Q_VOLIT      = MPA_VOLIT*418000._EB
          TMP_VEG_NEW  = TMP_VEG + (Q_FOR_VOLIT-Q_VOLIT)/(MPA_VEG*CP_VEG)
          TMP_VEG_NEW  = MIN(TMP_VEG_NEW,500._EB)
          WC%VEG_FUELMASS_L(IVEG_L) = MPA_VEG
          WC%ONE_D%MASSFLUX(I_FUEL)= WC%ONE_D%MASSFLUX(I_FUEL) + RDT_BC*MPA_VOLIT
!         WC%VEG_TMP_L(IVEG_L) = TMP_VEG_NEW
        ENDIF        

      ENDIF IF_VOLITIZATION

      ENDIF IF_DIVQ_L_GE_0
      
      IF(MPA_VEG <= MPA_VEG_MIN) TMP_VEG_NEW = TMP_G
      WC%VEG_TMP_L(IVEG_L) = TMP_VEG_NEW

    ENDDO LAYER_LOOP1

!   WC%VEG_TMP_L(LBURN) = WC%VEG_TMP_L(LBURN+1)
    WC%VEG_TMP_L(LBURN) = TMP_G

  endif if_fuel_elem_degrad

  ENDIF  IF_VEG_DEGRADATION_LINEAR

!      ************** Boundary Fuel Arrehnius Degradation model *************************
! Drying and pyrolysis occur according to Arrehnius expressions obtained 
! from the literature (Porterie et al., Num. Heat Transfer, 47:571-591, 2005
! Predicting wildland fire behavior and emissions using a fine-scale physical
! model

  IF_VEG_DEGRADATION_ARRHENIUS: IF(VEG_DEGRADATION_ARRHENIUS) THEN
    A_H2O_VEG      = 600000._EB !1/s sqrt(K)
    E_H2O_VEG      = 5800._EB !K
    A_PYR_VEG      = 36300._EB !1/s
    E_PYR_VEG      = 7250._EB !K
    L_PYR_VEG      = 418._EB !J/kg
    A_CHAR_VEG     = 430._EB !m/s
    E_CHAR_VEG     = 9000._EB !K
    BETA_CHAR_VEG  = 0.2_EB
    NU_CHAR_VEG    = 0.3_EB
    NU_ASH_VEG     = 0.1_EB
    NU_O2_CHAR_VEG = 1.65_EB
    H_CHAR_OXID    = -12.0E+6_EB !J/kg

    LAYER_LOOP2: DO IVEG_L = LBURN+1,NVEG_L

      MPA_MOIST = WC%VEG_MOISTMASS_L(IVEG_L)
      MPA_VEG   = WC%VEG_FUELMASS_L(IVEG_L)
      TMP_VEG   = WC%VEG_TMP_L(IVEG_L)

! Drying of vegetation (Arrhenius)
      IF_DEHYDRATION_2: IF (MPA_MOIST > MPA_MOIST_MIN) THEN
        MPA_MOIST_LOSS = MIN(DT_BC*MPA_MOIST*A_H2O_VEG*EXP(-E_H2O_VEG/TMP_VEG)/SQRT(TMP_VEG), &
                         MPA_MOIST-MPA_MOIST_MIN)
        MPA_MOIST      = MPA_MOIST - MPA_MOIST_LOSS
        WC%VEG_MOISTMASS_L(IVEG_L) = MPA_MOIST !kg/m^2
        IF (MPA_MOIST <= MPA_MOIST_MIN) WC%VEG_MOISTMASS_L(IVEG_L) = 0.0_EB
      ENDIF IF_DEHYDRATION_2

! Volitalization of vegetation(Arrhenius)
      IF_VOLITALIZATION_2: IF(MPA_VEG > MPA_VEG_MIN) THEN
        MPA_VOLIT    = CHAR_FCTR*DT_BC*MPA_VEG*A_PYR_VEG*EXP(-E_PYR_VEG/TMP_VEG)
        MPA_VOLIT    = MIN(MPA_VOLIT,(MPA_VEG-MPA_VEG_MIN))
        MPA_VEG      = MPA_VEG - MPA_VOLIT
        WC%VEG_FUELMASS_L(IVEG_L) = MPA_VEG
      ENDIF IF_VOLITALIZATION_2

      WC%ONE_D%MASSFLUX(I_FUEL)= WC%ONE_D%MASSFLUX(I_FUEL) + MPA_VOLIT*RDT_BC
      IF (I_WATER /= 0) WC%ONE_D%MASSFLUX(I_WATER) = WC%ONE_D%MASSFLUX(I_WATER) + MPA_MOIST*RDT_BC

! Vegetation temperature (Arrhenius)
      CP_VEG = (0.01_EB + 0.0037_EB*TMP_VEG)*1000._EB !W/kg/K
      CP_MOIST_AND_VEG = CP_H2O*WC%VEG_MOISTMASS_L(IVEG_L) +  CP_VEG*WC%VEG_FUELMASS_L(IVEG_L)

      WC%VEG_TMP_L(IVEG_L) = WC%VEG_TMP_L(IVEG_L) + (DT_BC*SF%VEG_DIVQNET_L(IVEG_L-LBURN) - &
                             (MPA_MOIST_LOSS*H_VAP_H2O + MPA_VOLIT*L_PYR_VEG) )/CP_MOIST_AND_VEG
      WC%VEG_TMP_L(IVEG_L) = MAX( WC%VEG_TMP_L(IVEG_L), TMPA)

    ENDDO LAYER_LOOP2

  ENDIF IF_VEG_DEGRADATION_ARRHENIUS
  
  WC%VEG_TMP_L(LBURN) = MAX(TMP_G,TMPA)
  WC%ONE_D%MASSFLUX_ACTUAL(I_FUEL) = WC%ONE_D%MASSFLUX(I_FUEL)
  IF (I_WATER /= 0) WC%ONE_D%MASSFLUX_ACTUAL(I_WATER) = WC%ONE_D%MASSFLUX(I_WATER)
 
! Temperature boundary condtions 
! Mass boundary conditions are determine in subroutine SPECIES_BC in wall.f90 for case SPECIFIED_MASS_FLUX
! TMP_F(IW) = WC%VEG_TMP_L(NVEG_L)
! IF (LBURN < NVEG_L)  TMP_F(IW) = WC%VEG_TMP_L(1+LBURN)
  IF (LBURN < NVEG_L) THEN
    WALL(IW)%ONE_D%TMP_F = WC%VEG_TMP_L(1+LBURN)
!   TMP_F(IW) = ((VEG_QRP_INC(0)+VEG_QRP_EMISS(0))/SIGMA)**.25 !as done in FDS4
  ELSE
    WALL(IW)%ONE_D%TMP_F = MAX(TMP_G,TMPA) !Tveg=Tgas if veg is completely burned
!   TMP_F(IW) = TMPA  !Tveg=Tambient if veg is completely burned
  ENDIF
! TMP_F(IW) = MAX(TMP_F(IW),TMPA)

ENDDO VEG_WALL_CELL_LOOP

VEG_CLOCK_BC = T

END SUBROUTINE BNDRY_VEG_MASS_ENERGY_TRANSFER
!
! ***********************************************************************************************
SUBROUTINE LEVEL_SET_FIRESPREAD(NM)
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
! 2) Need to use multiprocessors
!
!
INTEGER, INTENT(IN) :: NM
INTEGER :: I_FLANK,I,IM1,IM2,II,IIG,IP1,IP2,IW,J,JJG,JM1,JP1,LU_SLCF_LS,N_FINAL,N_FIRES,N_STEPS_OUT,N_TIME
REAL(EB) :: VEG_MOISTURE_LS
REAL(EB) :: LX,SR_MAX,T_FINAL,U_AMBIENT,UMAG,UMAX_LS, V_AMBIENT,VMAX_LS,XMAX_LS,XMIN_LS, &
            YMIN_LS,YMAX_LS
REAL(EB) :: G_EAST,G_WEST,G_SOUTH,G_NORTH
REAL(EB) :: HEAD_WIDTH,HEAD_WIDTH_FCTR,IGNITION_WIDTH_X,IGNITION_WIDTH_Y,T_IGN,X_MID,X_IGN_MIN,X_IGN_MAX
REAL(EB) :: R_FLANK
REAL(EB) :: R_HEAD_U0_INFW,R_HEAD_U_INFW
REAL(EB) :: DT_OUTPUT
REAL(EB), ALLOCATABLE, DIMENSION(:) :: X0_FIRE,Y0_FIRE,X_LS,Y_LS
LOGICAL :: HEAD_WIDTH_DEPENDENCE
!
REAL(FB) :: TIME_LS_OUT
!REAL(FB) :: MAG_SR
REAL(FB), ALLOCATABLE, DIMENSION(:,:) :: PHI_OUT
CHARACTER(30) :: SMOKEVIEW_LABEL,SMOKEVIEW_BAR_LABEL,UNITS

REAL(EB), POINTER, DIMENSION(:,:) :: ZT => NULL()

TYPE (WALL_TYPE),    POINTER :: WC =>NULL()
TYPE (SURFACE_TYPE), POINTER :: SF =>NULL()

CALL POINT_TO_MESH(NM)

ZT => LS_Z_TERRAIN

!print*,'level set: z(*)',z
!print*,'level set: ls_z_terrain(1,1)',ls_z_terrain(:,100)
!
!-Initialize variables
!
HEAD_WIDTH_DEPENDENCE = .FALSE.

!-- Domain specification (meters) from input file (assumes there's only one mesh)
!
 LX = XF - XS ; NX_LS = IBAR ; DX_LS = LX/REAL(NX_LS,EB)
 LX = YF - YS ; NY_LS = JBAR ; DY_LS = LX/REAL(NY_LS,EB)
 T_FINAL = T_END
 U_AMBIENT = U0 ; V_AMBIENT = V0 ; VEG_MOISTURE_LS = 5.8_EB
 YMIN_LS = YS
 YMAX_LS = YF
!
!-----------
! -- line fire ignition
!-----------
 T_IGN = -10._EB !set <0 for instantaneous line fire igntion
!
! -- Number of fires and their location
N_FIRES = 1 ; ALLOCATE(X0_FIRE(N_FIRES)) ; ALLOCATE(Y0_FIRE(N_FIRES))
X0_FIRE(1) = 0.0_EB ; Y0_FIRE(1) = 10._EB
!
! -- Constants in spread rate magnitude formula SR = R_FLANK + R_HEAD_U0_INFW + R_HEAD_U_INFW
!    R_FLANK = spread rate of flank fire, assumed to be wind indepent and
!              equal to the spread rate with zero wind on level ground
!    R_HEAD_U0_INFW = spread rate of head fire when ambient wind is zero
!    R_HEAD_U_INFW  = head fire spread rate coefficient which is mulitplied by the wind speed
!
! Australian grassland fires (infinite head width values)
R_HEAD_U0_INFW  = 0.165_EB*EXP(-0.108*VEG_MOISTURE_LS)/SQRT(U_AMBIENT**2 + V_AMBIENT**2)
R_HEAD_U_INFW   = 0.534_EB*EXP(-0.108*VEG_MOISTURE_LS)
!
! WFDS boundary fuel run AU F19, 8m (point ignitor), U2=7 m/s (matches U2=5m/s AU case)
R_FLANK         = 0.4 !flank ROS from wfds
R_HEAD_U0_INFW  = 0.0
R_HEAD_U_INFW   = 1.1
!need to divide by wind speed because this is the spread rate
!R_HEAD_U_INFW   = 1.1/SQRT(U_AMBIENT**2 + V_AMBIENT**2)

! Spread rates from WFDS boundary fuel run AU F19, 8m (point ignitor), U2=1 m/s
!R_HEAD_U0_INFW  = 0.0
!R_HEAD_U_INFW   = 0.25_EB/V_AMBIENT !need to divide by wind speed because this is the spread rate
!R_FLANK       = 0.195 !flank ROS from wfds
!
! Define spread rate across domain (including no burn areas)
!
ALLOCATE(ROS_HEAD(NX_LS,NY_LS))    ; CALL ChkMemErr('VEGE:LEVEL SET','ROS_HEAD',IZERO) ; ROS_HEAD=0.0_EB
!ALLOCATE(ROS_HEAD_U0_INFW(NX_LS,NY_LS)) ; CALL ChkMemErr('VEGE:LEVEL SET','ROS_HEAD_U0_INFW',IZERO)
!ALLOCATE(ROS_HEAD_U_INFW(NX_LS,NY_LS))  ; CALL ChkMemErr('VEGE:LEVEL SET','ROS_HEAD_U_INFW',IZERO)
ALLOCATE(ROS_FLANK(NX_LS,NY_LS))   ; CALL ChkMemErr('VEGE:LEVEL SET','ROS_FLANK',IZERO)
ALLOCATE(ROS_BACKU(NX_LS,NY_LS))   ; CALL ChkMemErr('VEGE:LEVEL SET','ROS_BACKU',IZERO)
ALLOCATE(FLANKFIRE_LIFETIME(NX_LS,NY_LS))   ; CALL ChkMemErr('VEGE:LEVEL SET','FLANKFIRE_LIFETIME',IZERO)
ALLOCATE(WIND_EXP(NX_LS,NY_LS))   ; CALL ChkMemErr('VEGE:LEVEL SET','WIND_EXP',IZERO)

FLANKFIRE_LIFETIME = 0.0_EB !handles finite lenght (lifetime) flankfires. For
!                            quenching flanks with lifetimes > TIME_FLANKFIRE_QUENCH
TIME_FLANKFIRE_QUENCH = 20.0_EB !flankfire lifetime in seconds

!ROS_HEAD_U0_INFW = R_HEAD_U0_INFW
!ROS_HEAD_U_INFW  = R_HEAD_U_INFW
!--- Bernardo Trails ROS from wfds runs
!ROS_HEAD         = 14._EB  !9.5 from wfds flat terrain run
!ROS_FLANK        = 4.0_EB  !4.2 from wfds flat terrain run   
!ROS_BACKU        = 0.1_EB !wind dependent (in general) back fire ROS need look up table for wind depend
!--- Bernardo Trails ROS from wfds runs
ROS_HEAD         = 3.9_EB  
ROS_FLANK        = 0.4_EB   
ROS_BACKU        = 0.0_EB 
! Currently slope dependence of head and back fires are from MkV Forest Danger Meter (hardcoded below)
!ROS_HEADS        = 2._EB  !slope dependent head fire ROS 
!ROS_BACKS        = 0.1_EB !slope dependent back fire ROS

!ROS_HEAD_U0_INFW(90:110,90:110) =  0.0_EB
!ROS_HEAD_U_INFW(90:110,90:110)  =  0.0_EB
!ROS_FLANK(90:110,90:110)  =  0.0_EB

! Assign spread rates (i.e., vegetation types) to locations on terrain

ROS_HEAD  = 0.0_EB
ROS_FLANK = 0.0_EB
ROS_BACKU = 0.0_EB
WIND_EXP  = 1.0_EB

!print*,'surface ros',surface%veg_lset_ros_head
!print*,'surface wind_exp',surface%veg_lset_wind_exp

print*,'before assign ROS'
ROS_WALL_CELL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
  WC  => WALL(IW)
  IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE ROS_WALL_CELL_LOOP

  SF  => SURFACE(WC%SURF_INDEX)

  IF (.NOT. SF%VEG_LSET_SPREAD) CYCLE ROS_WALL_CELL_LOOP



  IIG = WC%ONE_D%IIG
  JJG = WC%ONE_D%JJG
!print*,'IIG,JJG',iig,jjg
!print*,'ROS_HEAD',SF%VEG_LSET_ROS_HEAD
!print*,'ROS_HEAD,ROS_FLANK,ROS_BACK',SF%VEG_LSET_ROS_HEAD,SF%VEG_LSET_ROS_FLANK,SF%VEG_LSET_ROS_BACK

!print*,'IIG,JJG',iig,jjg
!print*,'ROS_HEAD',SF%VEG_LSET_ROS_HEAD
  ROS_HEAD(IIG,JJG)  = SF%VEG_LSET_ROS_HEAD
  ROS_FLANK(IIG,JJG) = SF%VEG_LSET_ROS_FLANK
  ROS_BACKU(IIG,JJG) = SF%VEG_LSET_ROS_BACK
  WIND_EXP(IIG,JJG)  = SF%VEG_LSET_WIND_EXP

! IF (SF%VEG_NO_BURN) THEN
!  ROS_HEAD_U0_INFW(IIG,JJG) =  0.0_EB
!  ROS_HEAD_U_INFW(IIG,JJG)  =  0.0_EB
!  ROS_HEAD(IIG,JJG)         =  0.0_EB
!  ROS_FLANK(IIG,JJG)        =  0.0_EB
!  ROS_BACKU(IIG,JJG)        =  0.0_EB
! ENDIF


ENDDO ROS_WALL_CELL_LOOP

print*,'ROS_HEAD max',MAXVAL(ROS_HEAD)

!
!C_F = 0.2_EB
!
! -- Flux limiter
LIMITER_LS = 1 !MINMOD
!LIMITER_LS = 2 !SUPERBEE
!LIMITER_LS = 3 !First order upwinding
!
! -- Output file
DT_OUTPUT  = 0.5_EB
TIME_LS    = 0._EB
LU_SLCF_LS = 9999
SMOKEVIEW_LABEL = 'phifield'
SMOKEVIEW_BAR_LABEL = 'phifield'
UNITS  = 'C'
OPEN(LU_SLCF_LS,FILE='lsfs.sf',FORM='UNFORMATTED',STATUS='REPLACE')
WRITE(LU_SLCF_LS) SMOKEVIEW_LABEL(1:30)
WRITE(LU_SLCF_LS) SMOKEVIEW_LABEL(1:30)
WRITE(LU_SLCF_LS) UNITS(1:30)
WRITE(LU_SLCF_LS)1,NX_LS,1,NY_LS,1,1
!
!
! =============== end of case specifications ========================
!
!-- Allocate arrays
ALLOCATE(U_LS(NX_LS,NY_LS)) ; CALL ChkMemErr('VEGE:LEVEL SET','U_LS',IZERO) ; U_LS = 0._EB
ALLOCATE(V_LS(NX_LS,NY_LS)) ; CALL ChkMemErr('VEGE:LEVEL SET','V_LS',IZERO) ; V_LS = 0._EB
ALLOCATE(PHI_LS(NX_LS,NY_LS)) ; CALL ChkMemErr('VEGE:LEVEL SET','PHI_LS',IZERO)
ALLOCATE(PHI0_LS(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','PHI0_LS',IZERO)
ALLOCATE(PHI1_LS(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','PHI1_LS',IZERO)
ALLOCATE(PHI_OUT(NX_LS,NY_LS)) ; CALL ChkMemErr('VEGE:LEVEL SET','PHI_OUT',IZERO) ; PHI_OUT = 0.0
ALLOCATE(SR_X_LS(NX_LS,NY_LS)) ; CALL ChkMemErr('VEGE:LEVEL SET','SR_X_LS',IZERO)
ALLOCATE(SR_Y_LS(NX_LS,NY_LS)) ; CALL ChkMemErr('VEGE:LEVEL SET','SR_Y_LS',IZERO)
ALLOCATE(FLUX0_LS(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','FLUX0_LS',IZERO)
ALLOCATE(FLUX1_LS(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','FLUX1_LS',IZERO)
ALLOCATE(DZTDX(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','DZDTX',IZERO)
ALLOCATE(DZTDY(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','DZDTY',IZERO)
ALLOCATE(MAG_ZT(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','MAG_ZT',IZERO)
ALLOCATE(MAG_SR_OUT(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','MAG_SR_OUT',IZERO)
!
!-- Computational grid
ALLOCATE(X_LS(NX_LS))   ; CALL ChkMemErr('VEGE:LEVEL SET','X_LS',IZERO)
ALLOCATE(Y_LS(NY_LS+1)) ; CALL ChkMemErr('VEGE:LEVEL SET','Y_LS',IZERO)

DO I = 0,NX_LS-1
!X_LS(I+1) = -0.5_EB*LX + 0.5_EB*DX_LS + DX_LS*REAL(I,EB)
 X_LS(I+1) = XS + 0.5_EB*DX_LS + DX_LS*REAL(I,EB)
ENDDO
!
DO J = 0,NY_LS
!Y_LS(J+1) = YMIN_LS + DY_LS*REAL(J,EB)
 Y_LS(J+1) = YS + DY_LS*REAL(J,EB)
ENDDO
XMIN_LS = X_LS(1) ; XMAX_LS = X_LS(NX_LS)
!
!-- Build wind field ------------------------------------------------------
!U_LS = U(1:NX_LS,1:NY_LS,1) ; use this if computed wind is to be used
!V_LS = V(1:NX_LS,1:NY_LS,1)

U_LS = U0
V_LS = V0

!U_LS = U_LS + U_AMBIENT
!V_LS = V_LS + V_AMBIENT
!
! Compute time step
!ROS_HEAD_U0_INFW = R_HEAD_U0_INFW
!ROS_HEAD_U_INFW  = R_HEAD_U_INFW
UMAX_LS  = MAXVAL(ABS(U_LS))
VMAX_LS  = MAXVAL(ABS(V_LS))
!UMAG     = SQRT(UMAX_LS**2 + VMAX_LS**2)
!SR_X_MAX = ( MAXVAL(ROS_HEAD_U0_INFW) + MAXVAL(ROS_HEAD_U_INFW) )*UMAX_LS
!SR_Y_MAX = ( MAXVAL(ROS_HEAD_U0_INFW) + MAXVAL(ROS_HEAD_U_INFW) )*VMAX_LS
!SR_MAX   = MAX(SR_Y_MAX,SR_X_MAX)
!SR_MAX   = MAX(SR_MAX,ROS_HEADS)
!SR_X_MAX =  MAXVAL(ROS_HEAD)
SR_MAX   = MAXVAL(ROS_HEAD)
print*,'SR_MAX',sr_max
SR_MAX   = MAX(SR_MAX,MAXVAL(ROS_FLANK))
print*,'SR_MAX',sr_max
SR_MAX   = 2._EB*SR_MAX !rough accounting for upslope spread aligned with wind
!DT_LS    = MIN(DX_LS,DY_LS)/SR_MAX
DT_LS = 0.25_EB*MIN(DX_LS,DY_LS)/SR_MAX
!DT_LS = 0.1603_EB !to make AU F19 ignition sequence work
N_STEPS_OUT = DT_OUTPUT/DT_LS
N_STEPS_OUT = MAX(N_STEPS_OUT,1)
N_FINAL     = MAX(T_FINAL/DT_LS,1._EB)
print*,'vege: t_final,dt_ls',t_final,dt_ls
!
!
!-- Initialize level set field. Fireline is at PHI_LS = 0 -------------------
!
!--- Gaussian shape
!PHI_MIN_LS = -1._EB
!PHI_MAX_LS = -1._EB
!AMP = 2._EB ; DS_LS = 0.01_EB*LX
!DO I = 1,NX_LS
! DO J = 1,NY_LS
!  PHI_LS(I,J) = AMP*EXP(-(X_LS(I)-X0_FIRE(1))**2/DS_LS**2   &
!                        -(Y_LS(J)-Y0_FIRE(1))**2/DS_LS**2) - 1
! ENDDO
!ENDDO
!
!--- Ignition line across entire domain
!PHI_MIN_LS = -1._EB
!PHI_MAX_LS = 1._EB
!LY0 = 0._EB
!MM = (PHI_MAX_LS-PHI_MIN_LS)/(2*YMIN_LS-LY0);
!BB = PHI_MAX _LS- YMIN_LS*MM
!DO I = 1,NX_LS
! DO J = 1,NY_LS
!  IF (Y_LS(J) < (LY0-YMIN_LS)) THEN
!   PHI_LS(I,J) = MM*Y_LS(J)+BB
! ELSE
!   PHI_LS(I,J) = PHI_MIN_LS
!  ENDIF
! ENDDO
!ENDDO
!
!--- Finite width ignition line
PHI_MIN_LS = -1._EB
PHI_MAX_LS = 1._EB
IGNITION_WIDTH_X = 8._EB
IGNITION_WIDTH_Y = 8._EB
X_MID     = XMIN_LS + 0.5*(XMIN_LS + XMAX_LS)
X_IGN_MIN = X_MID - 0.5*IGNITION_WIDTH_X
X_IGN_MAX = X_IGN_MIN + IGNITION_WIDTH_X
PHI_LS = PHI_MIN_LS
!DO I = 1,NX_LS
! DO J = 1,NY_LS
!! IF (Y_LS(J) <  YMIN_LS + 0.02_EB*(YMAX_LS-YMIN_LS) .AND. &
!!     X_LS(I) >= XMIN_LS + 0.55_EB*(XMAX_LS-XMIN_LS) .AND. &
!!     X_LS(I) <= XMIN_LS + 0.65_EB*(XMAX_LS-XMIN_LS) ) THEN
!!  PHI_LS(I,J) = PHI_MAX_LS
!!ELSE
!!  PHI_LS(I,J) = PHI_MIN_LS
!! ENDIF
!  IF (Y_LS(J) <= IGNITION_WIDTH_Y .AND. X_LS(I) >= X_IGN_MIN .AND. X_LS(I) <= X_IGN_MAX) THEN
!    PHI_LS(I,J) = PHI_MAX_LS
!  ELSE
!    PHI_LS(I,J) = PHI_MIN_LS
!  ENDIF
! ENDDO
!ENDDO
!
!--- Time dependent ignition line
!N_STEPS_IGN = 1._EB/DT_LS !ignitor every 1 s
!DELTAT_IGN  = N_STEPS_IGN*DT_LS
!N_STEPS_IGN = 6 !AU C064
!!N_STEPS_IGN = 4 !AU F19 with dx=1m and DT_LS defined to be 0.1603 s
!PHI_MIN_LS = -1._EB ; PHI_MAX_LS = 1._EB
!PHI_LS = PHI_MIN_LS
!N_IGNITION = 0 !counter for ignitions
!DX_IGN     = 1._EB !m, length of ignition segments in each direction
!XIGNL_MIN  = -DX_IGN ; XIGNL_MAX = 0._EB !bounds of left ignitor
!XIGNR_MIN  =  0._EB  ; XIGNR_MAX = DX_IGN !bounds of right ignitor
!HEAD_WIDTH = 2._EB*DX_IGN

!---- Compute components of slope gradient and magnitude of gradient

GRADIENT_ILOOP: DO I = 1,NX_LS
 IM1=I-1 ; IM2=I-2
 IP1=I+1 ; IP2=I+2

 DO J = 2,NY_LS-1
   JM1=J-1
   JP1=J+1

   G_EAST  = 0.5*( ZT(I,J) + ZT(IP1,J) )
   G_WEST  = 0.5*( ZT(I,J) + ZT(IM1,J) )
   G_NORTH = 0.5*( ZT(I,J) + ZT(I,JP1) )
   G_SOUTH = 0.5*( ZT(I,J) + ZT(I,JM1) )

   DZTDX(I,J) = (G_EAST-G_WEST)/DX_LS
   DZTDY(I,J) = (G_NORTH-G_SOUTH)/DY_LS
   MAG_ZT(I,J) = SQRT(DZTDX(I,J)**2 + DZTDY(I,J)**2)
 ENDDO

ENDDO GRADIENT_ILOOP


!
!-- Time step solution using second order Runge-Kutta -----------------------
!
print*,'vege: n_final',n_final
TIMESTEP: DO N_TIME = 1,N_FINAL
!
!--- Put in time dependent ignitor
!
! IF (TIME_LS <= T_IGN) THEN 
!  IF (MOD(N_TIME,N_STEPS_IGN) < 0.0001_EB .AND. HEAD_WIDTH < HEAD_WIDTH_MAX) THEN
!   DO I = 1,NX_LS
!    DO J = 1,2
!     IF ( X_LS(I) >= XIGNL_MIN .AND. X_LS(I) <= XIGNL_MAX ) THEN
!      PHI_LS(I,J) = PHI_MAX_LS
!      PHI_MIN_LS = PHI_MIN_LS
!     ENDIF
!     IF ( X_LS(I) >= XIGNR_MIN .AND. X_LS(I) <= XIGNR_MAX ) THEN
!      PHI_LS(I,J) = PHI_MAX_LS
!      PHI_MIN_LS = PHI_MIN_LS
!     ENDIF
!    ENDDO
!   ENDDO
!   XIGNL_MIN = XIGNL_MIN - DX_IGN
!   XIGNL_MAX = XIGNL_MAX - DX_IGN
!   XIGNR_MIN = XIGNR_MIN + DX_IGN
!   XIGNR_MAX = XIGNR_MAX + DX_IGN
!   N_IGNITION = N_IGNITION + 1
!!  HEAD_WIDTH_FCTR = EXP(-(0.859_EB + 2.036_EB*UMAG)/HEAD_WIDTH)
!!  IF (.NOT. HEAD_WIDTH_DEPENDENCE) HEAD_WIDTH_FCTR = 1._EB
!!  ROS_HEAD = (R_HEAD_U0_INFW + R_HEAD_U_INFW)*HEAD_WIDTH_FCTR - ROS_FLAMK/UMAG
!   HEAD_WIDTH = HEAD_WIDTH + 2._EB*DX_IGN
!  ENDIF
! ENDIF

!-- Find flank-to-flank distance at base of fire assume symmetry about xmid and
!   define spread rate

 IF (HEAD_WIDTH_DEPENDENCE) THEN
  print*,'NEED TO MODIFY CODE TO HANDLE HEAD WIDTH DEPENDENT SPREAD RATE'
  I_FLANK = 0
  DO II = NX_LS/2,NX_LS
   IF(PHI_LS(II,1) <= 0.0_EB .AND. I_FLANK==0) I_FLANK = II
  ENDDO
  HEAD_WIDTH = 2._EB*(I_FLANK - NX_LS/2)*DX_LS
  IF (ABS(TIME_LS) <= ZERO_P) HEAD_WIDTH = IGNITION_WIDTH_X
  HEAD_WIDTH_FCTR = EXP(-(0.859_EB + 2.036_EB*UMAG)/HEAD_WIDTH)
  ROS_HEAD = (ROS_HEAD_U0_INFW + ROS_HEAD_U_INFW)*HEAD_WIDTH_FCTR - ROS_FLANK/UMAG
 ELSE
! ROS_HEAD = ROS_HEAD_U_INFW
! ROS_HEAD = (ROS_HEAD_U0_INFW + ROS_HEAD_U_INFW) - ROS_FLANK
! ROS_HEAD = (ROS_HEAD_U0_INFW + ROS_HEAD_U_INFW) - ROS_FLANK/UMAG
 ENDIF

! -- Ignite landscape at user specified location(s) and time(s)
! ** change to go into the wall cell loop only if time corrensponds to an
! ** ignition time. Need to create a separate array with sorted ignition times

IGNITOR_WALL_CELL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
  WC  => WALL(IW)
  SF  => SURFACE(WC%SURF_INDEX)


  IIG = WC%ONE_D%IIG
  JJG = WC%ONE_D%JJG

  IF (SF%VEG_LSET_IGNITE_T >= TIME_LS .AND. SF%VEG_LSET_IGNITE_T <= TIME_LS + DT_LS) PHI_LS(IIG,JJG) = PHI_MAX_LS 

ENDDO IGNITOR_WALL_CELL_LOOP


!
!--- RK Stage 1
 RK2_PREDICTOR_LS = .TRUE.
 CALL LEVEL_SET_SPREAD_RATE
 CALL LEVEL_SET_ADVECT_FLUX
 PHI1_LS = PHI_LS - DT_LS*FLUX0_LS
 !PHI1_LS = MAX(PHI1_LS,PHI_MIN_LS)
 !PHI1_LS = MIN(PHI1_LS,PHI_MAX_LS)

!--- RK Stage2
 RK2_PREDICTOR_LS = .FALSE.
 MAG_SR_OUT       = 0.0
 CALL LEVEL_SET_SPREAD_RATE
 CALL LEVEL_SET_ADVECT_FLUX
 PHI_LS = PHI_LS - 0.5_EB*DT_LS*(FLUX0_LS + FLUX1_LS)
 !PHI_LS = MAX(PHI_LS,PHI_MIN_LS)
 !PHI_LS = MIN(PHI_LS,PHI_MAX_LS)
!
 TIME_LS = TIME_LS + DT_LS

!--- Output slice file for smokeview
 IF (MOD(N_TIME,N_STEPS_OUT) < 0.0001_EB) THEN
  PHI_OUT = PHI_LS
  TIME_LS_OUT = TIME_LS
  WRITE(LU_SLCF_LS) TIME_LS_OUT
!negative for consistency with wall thicknes output from wfds
  WRITE(LU_SLCF_LS) ((-PHI_OUT(I,J),I=1,NX_LS),J=1,NY_LS) 

! output magnitude of spread rate
!MAG_SR_OUT = 0.0
!DO I = 1,NX_LS
! DO J = 1,NY_LS
!!   MAG_SR = SQRT(SR_X_LS(I,J)**2 + SR_Y_LS(I,J)**2)
!!   IF(MAG_SR  > PHI_OUT(I,J)) PHI_OUT(I,J) = MAG_SR
!!   PHI_OUT(I,J) = MAG_SR
!  IF(PHI_LS(I,J) <= 0.2_EB .AND. PHI_LS(I,J) >= -0.2) MAG_SR_OUT(I,J) =-FLUX1_LS(I,J)
!  IF(PHI_OUT(I,J) < MAG_SR_OUT(I,J)) PHI_OUT(I,J) = MAG_SR_OUT(I,J)
! ENDDO
!ENDDO
!! WRITE(LU_SLCF_LS) ((PHI_OUT(I,J),I=1,NX_LS),J=1,NY_LS) 
!! PHI_OUT = SR_X_LS + SR_Y_LS
! WRITE(LU_SLCF_LS) ((MAG_SR_OUT(I,J),I=1,NX_LS),J=1,NY_LS) 

 ENDIF
!
ENDDO TIMESTEP
!
CLOSE(LU_SLCF_LS)
!
END SUBROUTINE LEVEL_SET_FIRESPREAD


SUBROUTINE LEVEL_SET_SPREAD_RATE
!
! Compute components of spread rate vector
!
INTEGER :: I,J,IM1,IM2,IP1,IP2,JM1,JP1
REAL(EB) :: COS_THETA_WIND,COS_THETA_SLOPE,COS_THETA_WIND_H,COS_THETA_WIND_B, &
            COS_THETA_SLOPE_H,COS_THETA_SLOPE_B,DPHIDX,DPHIDY,F_EAST,F_WEST,F_NORTH,F_SOUTH, &
            GRAD_SLOPE_DOT_NORMAL_FIRELINE,MAG_F,MAG_SR,MAG_U,WIND_DOT_NORMAL_FIRELINE,NEXP_WIND
REAL(EB) :: RAD_TO_DEGREE,DEGREES_SLOPE
REAL(EB), DIMENSION(:)   :: NORMAL_FIRELINE(2)
 
RAD_TO_DEGREE = 90._EB/ASIN(1._EB)
!NEXP_WIND = 2


IF (RK2_PREDICTOR_LS) PHI0_LS = PHI_LS
IF (.NOT. RK2_PREDICTOR_LS) PHI0_LS = PHI1_LS
SR_X_LS = 0.0_EB ; SR_Y_LS = 0.0_EB

FLUX_ILOOP: DO I = 1,NX_LS
!
 IM1=I-1; IF (IM1<1) IM1=IM1+NX_LS
 IM2=I-2; IF (IM2<1) IM2=IM2+NX_LS

 IP1=I+1; IF (IP1>NX_LS) IP1=IP1-NX_LS
 IP2=I+2; IF (IP2>NX_LS) IP2=IP2-NX_LS

  DO J = 2,NY_LS-1
   JM1=J-1
   JP1=J+1

   F_EAST  = 0.5*( PHI0_LS(I,J) + PHI0_LS(IP1,J) )
   F_WEST  = 0.5*( PHI0_LS(I,J) + PHI0_LS(IM1,J) )
   F_NORTH = 0.5*( PHI0_LS(I,J) + PHI0_LS(I,JP1) )
   F_SOUTH = 0.5*( PHI0_LS(I,J) + PHI0_LS(I,JM1) )

   DPHIDX = (F_EAST-F_WEST)/DX_LS
   DPHIDY = (F_NORTH-F_SOUTH)/DY_LS
   MAG_F = SQRT(DPHIDX**2 + DPHIDY**2)
   IF (MAG_F > 0._EB) THEN   !components of unit vector normal to PHI contours
    NORMAL_FIRELINE(1) = -DPHIDX/MAG_F
    NORMAL_FIRELINE(2) = -DPHIDY/MAG_F
   ELSE
    NORMAL_FIRELINE = 0._EB
   ENDIF
                
   WIND_DOT_NORMAL_FIRELINE = U_LS(I,J)*NORMAL_FIRELINE(1) + V_LS(I,J)*NORMAL_FIRELINE(2)
   MAG_U  = SQRT(U_LS(I,J)**2 + V_LS(I,J)**2)

   COS_THETA_WIND = 0.0_EB ; COS_THETA_WIND_H = 0.0_EB ; COS_THETA_WIND_B = 0.0_EB
   IF(MAG_U > 0.0_EB) COS_THETA_WIND = WIND_DOT_NORMAL_FIRELINE/MAG_U

   GRAD_SLOPE_DOT_NORMAL_FIRELINE = DZTDX(I,J)*NORMAL_FIRELINE(1) + DZTDY(I,J)*NORMAL_FIRELINE(2) 
   COS_THETA_SLOPE = 0.0_EB ; COS_THETA_SLOPE_H = 0.0_EB ; COS_THETA_SLOPE_B = 0.0_EB
   IF(MAG_ZT(I,J) > 0.0_EB) COS_THETA_SLOPE = GRAD_SLOPE_DOT_NORMAL_FIRELINE/MAG_ZT(I,J)

   DEGREES_SLOPE = ATAN(MAG_ZT(I,J))*RAD_TO_DEGREE
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
    IF(DEGREES_SLOPE >= 5._EB .AND. DEGREES_SLOPE < 10._EB)  ROS_HEADS = 0.33_EB*ROS_HEAD(I,J)
    IF(DEGREES_SLOPE >= 10._EB .AND. DEGREES_SLOPE < 20._EB) ROS_HEADS =         ROS_HEAD(I,J)
    IF(DEGREES_SLOPE >= 20._EB)                              ROS_HEADS =  3._EB*ROS_HEAD(I,J)
    MAG_SR = ROS_FLANK(I,J)*(1._EB + COS_THETA_WIND**NEXP_WIND*COS_THETA_SLOPE) + &
             (ROS_HEAD(I,J) - ROS_FLANK(I,J))*COS_THETA_WIND**NEXP_WIND + &
             (ROS_HEADS     - ROS_FLANK(I,J))*COS_THETA_SLOPE  !magnitude of spread rate
   ENDIF
!  IF(ABS(COS_THETA_WIND) < 0.5_EB .AND. MAG_F > 0._EB) MAG_SR = 0.0_EB
!  IF(ABS(COS_THETA_WIND) < 0.5_EB .AND. MAG_F > 0._EB) FLANKFIRE_LIFETIME(I,J) = FLANKFIRE_LIFETIME(I,J) + DT_LS
!  IF(FLANKFIRE_LIFETIME(I,J) > TIME_FLANKFIRE_QUENCH) MAG_SR = 0.0_EB

! Spread with the wind and downsslope
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
    IF(DEGREES_SLOPE >= 5._EB .AND. DEGREES_SLOPE < 10._EB)  ROS_BACKS = -0.33_EB*ROS_BACKU(I,J)
    IF(DEGREES_SLOPE >= 10._EB .AND. DEGREES_SLOPE < 20._EB) ROS_BACKS =         -ROS_BACKU(I,J)
    IF(DEGREES_SLOPE >= 20._EB)                              ROS_BACKS = -3.0_EB*ROS_BACKU(I,J)
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
!  MAG_SR = R_FLANK + ROS_HEAD*WIND_DOT_NORMAL_FIRELINE !magnitude of spread rate
   SR_X_LS(I,J) = MAG_SR*NORMAL_FIRELINE(1) !spread rate components
   SR_Y_LS(I,J) = MAG_SR*NORMAL_FIRELINE(2) 
!  MAG_SR_OUT(I,J) = MAG_SR
     
  ENDDO

!-- Inlet boundary
  DO J = 1,1

   JP1 = J+1
        
   F_EAST  = 0.5*( PHI0_LS(I,J) + PHI0_LS(IP1,J) )
   F_WEST  = 0.5*( PHI0_LS(I,J) + PHI0_LS(IM1,J) )
   F_NORTH = 0.5*( PHI0_LS(I,J) + PHI0_LS(I,JP1) )
   F_SOUTH = 0.5*( PHI0_LS(I,J) + PHI_MIN_LS ) !only for finite width ignition line
!   F_SOUTH = 0.5*( PHI0_LS(I,J) + PHI_MAX_LS )
        
   DPHIDX = (F_EAST-F_WEST)/DX_LS
   DPHIDY = (F_NORTH-F_SOUTH)/DY_LS
   MAG_F = SQRT(DPHIDX**2 +  DPHIDY**2)
   If (MAG_F > 0._EB) THEN
    NORMAL_FIRELINE(1) = -DPHIDX/MAG_F
    NORMAL_FIRELINE(2) = -DPHIDY/MAG_F
   ELSE
    NORMAL_FIRELINE = 0._EB
   ENDIF
        
   WIND_DOT_NORMAL_FIRELINE = U_LS(I,J)*NORMAL_FIRELINE(1) + V_LS(I,J)*NORMAL_FIRELINE(2)
   MAG_U  = SQRT(U_LS(I,J)**2 + V_LS(I,J)**2)

   COS_THETA_WIND = 0.0_EB ; COS_THETA_WIND_H = 0.0_EB ; COS_THETA_WIND_B = 0.0_EB
   IF(MAG_U > 0.0_EB) COS_THETA_WIND = WIND_DOT_NORMAL_FIRELINE/MAG_U

   GRAD_SLOPE_DOT_NORMAL_FIRELINE = DZTDX(I,J)*NORMAL_FIRELINE(1) + DZTDY(I,J)*NORMAL_FIRELINE(2) 
   COS_THETA_SLOPE = 0.0_EB ; COS_THETA_SLOPE_H = 0.0_EB ; COS_THETA_SLOPE_B = 0.0_EB
   IF(MAG_ZT(I,J) > 0.0_EB) COS_THETA_SLOPE = GRAD_SLOPE_DOT_NORMAL_FIRELINE/MAG_ZT(I,J)

   DEGREES_SLOPE = ATAN(MAG_ZT(I,J))*RAD_TO_DEGREE
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
    IF(DEGREES_SLOPE >= 5._EB .AND. DEGREES_SLOPE < 10._EB)  ROS_HEADS = 0.33_EB*ROS_HEAD(I,J)
    IF(DEGREES_SLOPE >= 10._EB .AND. DEGREES_SLOPE < 20._EB) ROS_HEADS =         ROS_HEAD(I,J)
    IF(DEGREES_SLOPE >= 20._EB)                              ROS_HEADS =  3._EB*ROS_HEAD(I,J)
    MAG_SR = ROS_FLANK(I,J)*(1._EB + COS_THETA_WIND**NEXP_WIND*COS_THETA_SLOPE) + &
             (ROS_HEAD(I,J) - ROS_FLANK(I,J))*COS_THETA_WIND**NEXP_WIND + &
             (ROS_HEADS     - ROS_FLANK(I,J))*COS_THETA_SLOPE  !magnitude of spread rate
!   if(cos_theta_wind == 0._EB) FLANKFIRE_LIFETIME(I,J) = FLANKFIRE_LIFETIME(I,J) + DT_LS
!   if(flankfire_lifetime(i,j) > time_flankfire_quench) mag_sr = 0.0_EB
   ENDIF

! Spread with the wind and downsslope
   IF(COS_THETA_WIND >= 0._EB .AND. COS_THETA_SLOPE < 0._EB) THEN
    IF(DEGREES_SLOPE >= 5._EB .AND. DEGREES_SLOPE < 10._EB)  ROS_HEADS =  0.33_EB*ROS_HEAD(I,J)
    IF(DEGREES_SLOPE >= 10._EB .AND. DEGREES_SLOPE < 20._EB) ROS_HEADS =  0.50_EB*ROS_HEAD(I,J)
    IF(DEGREES_SLOPE >= 20._EB)                              ROS_HEADS =  0.75_EB*ROS_HEAD(I,J)
    MAG_SR = ROS_FLANK(I,J)*(1._EB + COS_THETA_WIND**NEXP_WIND*COS_THETA_SLOPE) + &
             (ROS_HEAD(I,J) - ROS_FLANK(I,J))*COS_THETA_WIND**NEXP_WIND + &
             (ROS_HEADS     - ROS_FLANK(I,J))*COS_THETA_SLOPE  !magnitude of spread rate
!   if(cos_theta_wind == 0._EB) FLANKFIRE_LIFETIME(I,J) = FLANKFIRE_LIFETIME(I,J) + DT_LS
!   if(flankfire_lifetime(i,j) > time_flankfire_quench) mag_sr = 0.0_EB
   ENDIF

! Spread against the wind and upslope
   IF(COS_THETA_WIND <  0._EB .AND. COS_THETA_SLOPE >= 0._EB) THEN
    IF(DEGREES_SLOPE >= 5._EB .AND. DEGREES_SLOPE < 10._EB)  ROS_BACKS = -0.33_EB*ROS_BACKU(I,J)
    IF(DEGREES_SLOPE >= 10._EB .AND. DEGREES_SLOPE < 20._EB) ROS_BACKS =         -ROS_BACKU(I,J)
    IF(DEGREES_SLOPE >= 20._EB)                              ROS_BACKS = -3.0_EB*ROS_BACKU(I,J)
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
!  MAG_SR = R_FLANK + ROS_HEAD*WIND_DOT_NORMAL_FIRELINE !magnitude of spread rate
   SR_X_LS(I,J) = MAG_SR*NORMAL_FIRELINE(1) 
   SR_Y_LS(I,J) = MAG_SR*NORMAL_FIRELINE(2) 
!  MAG_SR_OUT(I,J) = MAG_SR
        
  ENDDO

!-- Outlet boundary
  DO J = NY_LS,NY_LS
        
   JM1 = J-1

   F_EAST =  0.5*( PHI0_LS(I,J) + PHI0_LS(IP1,J) )
   F_WEST =  0.5*( PHI0_LS(I,J) + PHI0_LS(IM1,J) )
   F_NORTH = 0.5*( PHI0_LS(I,J) + PHI_MIN_LS )
   F_SOUTH = 0.5*( PHI0_LS(I,J) + PHI0_LS(I,JM1) )
        
   DPHIDX = (F_EAST-F_WEST)/DX_LS
   DPHIDY = (F_NORTH-F_SOUTH)/DY_LS
   MAG_F = SQRT(DPHIDX**2 + DPHIDY**2)
   IF (MAG_F > 0._EB) THEN
    NORMAL_FIRELINE(1) = -DPHIDX/MAG_F
    NORMAL_FIRELINE(2) = -DPHIDY/MAG_F
   ELSE
    NORMAL_FIRELINE = 0._EB
   ENDIF
        
   WIND_DOT_NORMAL_FIRELINE = U_LS(I,J)*NORMAL_FIRELINE(1) + V_LS(I,J)*NORMAL_FIRELINE(2)
   MAG_U  = SQRT(U_LS(I,J)**2 + V_LS(I,J)**2)

   COS_THETA_WIND = 0.0_EB ; COS_THETA_WIND_H = 0.0_EB ; COS_THETA_WIND_B = 0.0_EB
   IF(MAG_U > 0.0_EB) COS_THETA_WIND = WIND_DOT_NORMAL_FIRELINE/MAG_U

   GRAD_SLOPE_DOT_NORMAL_FIRELINE = DZTDX(I,J)*NORMAL_FIRELINE(1) + DZTDY(I,J)*NORMAL_FIRELINE(2) 
   COS_THETA_SLOPE = 0.0_EB ; COS_THETA_SLOPE_H = 0.0_EB ; COS_THETA_SLOPE_B = 0.0_EB
   IF(MAG_ZT(I,J) > 0.0_EB) COS_THETA_SLOPE = GRAD_SLOPE_DOT_NORMAL_FIRELINE/MAG_ZT(I,J)

   DEGREES_SLOPE = ATAN(MAG_ZT(I,J))*RAD_TO_DEGREE
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
    IF(DEGREES_SLOPE >= 5._EB .AND. DEGREES_SLOPE < 10._EB)  ROS_HEADS = 0.33_EB*ROS_HEAD(I,J)
    IF(DEGREES_SLOPE >= 10._EB .AND. DEGREES_SLOPE < 20._EB) ROS_HEADS =         ROS_HEAD(I,J)
    IF(DEGREES_SLOPE >= 20._EB)                              ROS_HEADS =  3._EB*ROS_HEAD(I,J)
    MAG_SR = ROS_FLANK(I,J)*(1._EB + COS_THETA_WIND**NEXP_WIND*COS_THETA_SLOPE) + &
             (ROS_HEAD(I,J) - ROS_FLANK(I,J))*COS_THETA_WIND**NEXP_WIND + &
             (ROS_HEADS     - ROS_FLANK(I,J))*COS_THETA_SLOPE  !magnitude of spread rate
!   if(cos_theta_wind == 0._EB) FLANKFIRE_LIFETIME(I,J) = FLANKFIRE_LIFETIME(I,J) + DT_LS
!   if(flankfire_lifetime(i,j) > time_flankfire_quench) mag_sr = 0.0_EB
   ENDIF

! Spread with the wind and downsslope
   IF(COS_THETA_WIND >= 0._EB .AND. COS_THETA_SLOPE < 0._EB) THEN
    IF(DEGREES_SLOPE >= 5._EB .AND. DEGREES_SLOPE < 10._EB)  ROS_HEADS =  0.33_EB*ROS_HEAD(I,J)
    IF(DEGREES_SLOPE >= 10._EB .AND. DEGREES_SLOPE < 20._EB) ROS_HEADS =  0.50_EB*ROS_HEAD(I,J)
    IF(DEGREES_SLOPE >= 20._EB)                              ROS_HEADS =  0.75_EB*ROS_HEAD(I,J)
    MAG_SR = ROS_FLANK(I,J)*(1._EB + COS_THETA_WIND**NEXP_WIND*COS_THETA_SLOPE) + &
             (ROS_HEAD(I,J) - ROS_FLANK(I,J))*COS_THETA_WIND**NEXP_WIND + &
             (ROS_HEADS     - ROS_FLANK(I,J))*COS_THETA_SLOPE  !magnitude of spread rate
!   if(cos_theta_wind == 0._EB) FLANKFIRE_LIFETIME(I,J) = FLANKFIRE_LIFETIME(I,J) + DT_LS
!   if(flankfire_lifetime(i,j) > time_flankfire_quench) mag_sr = 0.0_EB
   ENDIF

! Spread against the wind and upslope
   IF(COS_THETA_WIND <  0._EB .AND. COS_THETA_SLOPE >= 0._EB) THEN
    IF(DEGREES_SLOPE >= 5._EB .AND. DEGREES_SLOPE < 10._EB)  ROS_BACKS = -0.33_EB*ROS_BACKU(I,J)
    IF(DEGREES_SLOPE >= 10._EB .AND. DEGREES_SLOPE < 20._EB) ROS_BACKS =         -ROS_BACKU(I,J)
    IF(DEGREES_SLOPE >= 20._EB)                              ROS_BACKS = -3.0_EB*ROS_BACKU(I,J)
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
!  MAG_SR = R_FLANK + ROS_HEAD*WIND_DOT_NORMAL_FIRELINE !magnitude of spread rate
   SR_X_LS(I,J) = MAG_SR*NORMAL_FIRELINE(1) 
   SR_Y_LS(I,J) = MAG_SR*NORMAL_FIRELINE(2) 
!  MAG_SR_OUT(I,J) = MAG_SR
        
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

 DO J = 3,NY_LS-2
   JM1=J-1
   JM2=J-2
   JP1=J+1
   JP2=J+2

!-- east face
   Z(1) = PHI0_LS(IM1,J)
   Z(2) = PHI0_LS(I,J)
   Z(3) = PHI0_LS(IP1,J)
   Z(4) = PHI0_LS(IP2,J)
   F_EAST = SCALAR_FACE_VALUE(SR_X_LS(I,J),Z,LIMITER_LS)

!-- west face
   Z(1) = PHI0_LS(IM2,J)
   Z(2) = PHI0_LS(IM1,J)
   Z(3) = PHI0_LS(I,J)
   Z(4) = PHI0_LS(IP1,J)
   F_WEST = SCALAR_FACE_VALUe(SR_X_LS(I,J),Z,LIMITER_LS)

!    north face
   Z(1) = PHI0_LS(I,JM1)
   Z(2) = PHI0_LS(I,J)
   Z(3) = PHI0_LS(I,JP1)
   Z(4) = PHI0_LS(I,JP2)
   F_NORTH = SCALAR_FACE_VALUE(SR_Y_LS(I,J),Z,LIMITER_LS)

!    south face
   Z(1) = PHI0_LS(I,JM2)
   Z(2) = PHI0_LS(I,JM1)
   Z(3) = PHI0_LS(I,J)
   Z(4) = PHI0_LS(I,JP1)
   F_SOUTH = SCALAR_FACE_VALUE(SR_Y_LS(I,J),Z,LIMITER_LS)
        
   DPHIDX = (F_EAST-F_WEST)/DX_LS
   DPHIDY = (F_NORTH-F_SOUTH)/DY_LS
   FLUX_LS(I,J) = SR_X_LS(I,J)*DPHIDX + SR_Y_LS(I,J)*DPHIDY
   
   PHIMAG          = SQRT(DPHIDX**2 + DPHIDY**2)
   MAG_SR_OUT(I,J) = 0.0_EB
   IF(PHIMAG > 0.0_EB) MAG_SR_OUT(I,J) = FLUX_LS(I,J)/PHIMAG
        
!  fx = (f_east-f_west)/dx
!  fy = (f_north-f_south)/dy
!       phi(i,j) = phi0(i,j) - dt*[Fx(i,j) Fy(i,j)]*[fx fy]
 ENDDO

! Inlet boundary
 DO J = 1,2
   JM1 = J-1

   JP1 = J+1
   JP2 = J+2

!    east face
   Z(1) = PHI0_LS(IM1,J)
   Z(2) = PHI0_LS(I,J)
   Z(3) = PHI0_LS(IP1,J)
   Z(4) = PHI0_LS(IP2,J)
   F_EAST = SCALAR_FACE_VALUE(SR_X_LS(I,J),Z,LIMITER_LS)

!    west face
   Z(1) = PHI0_LS(IM2,J)
   Z(2) = PHI0_LS(IM1,J)
   Z(3) = PHI0_LS(I,J)
   Z(4) = PHI0_LS(IP1,J)
   F_WEST = SCALAR_FACE_VALUE(SR_X_LS(I,J),Z,LIMITER_LS)

   IF (J==1) THEN
!    north face
    Z(1) = PHI_MAX_LS
    Z(2) = PHI0_LS(I,J)
    Z(3) = PHI0_LS(I,JP1)
    Z(4) = PHI0_LS(I,JP2)
    F_NORTH = SCALAR_FACE_VALUE(SR_Y_LS(I,J),Z,LIMITER_LS)

!    south face
    Z(1) = PHI_MAX_LS
    Z(2) = PHI_MAX_LS
    Z(3) = PHI0_LS(I,J)
    Z(4) = PHI0_LS(I,JP1)
    F_SOUTH = SCALAR_FACE_VALUE(SR_Y_LS(I,J),Z,LIMITER_LS)

   ELSEIF (j==2) THEN
!    north face
    Z(1) = PHI0_LS(I,JM1)
    Z(2) = PHI0_LS(I,J)
    Z(3) = PHI0_LS(I,JP1)
    Z(4) = PHI0_LS(I,JP2)
    F_NORTH = SCALAR_FACE_VALUE(SR_Y_LS(I,J),Z,LIMITER_LS)

!    south face
    Z(1) = PHI_MAX_LS
    Z(2) = PHI0_LS(I,JM1)
    Z(3) = PHI0_LS(I,J)
    Z(4) = PHI0_LS(I,JP1)
    F_SOUTH = SCALAR_FACE_VALUE(SR_Y_LS(I,J),Z,LIMITER_LS)
   ENDIF
        
   DPHIDX = (F_EAST-F_WEST)/DX_LS
   DPHIDY = (F_NORTH-F_SOUTH)/DY_LS
   FLUX_LS(I,J) = SR_X_LS(I,J)*DPHIDX + SR_Y_LS(I,J)*DPHIDY
   
   PHIMAG          = SQRT(DPHIDX**2 + DPHIDY**2)
   MAG_SR_OUT(I,J) = 0.0_EB
   IF(PHIMAG > 0.0_EB) MAG_SR_OUT(I,J) = FLUX_LS(I,J)/PHIMAG
        
!   phi(i,j) = phi0(i,j) - dt*[Fx(i,j) Fy(i,j)]*[fx fy]
 ENDDO
 FLUX_LS(:,1) = FLUX_LS(:,2)

! outlet boundary
 DO J = NY_LS-1,NY_LS
   JM1 = J-1
   JM2 = J-2

   JP1 = J+1

!    east face
   Z(1) = PHI0_LS(IM1,J)
   Z(2) = PHI0_LS(I,J)
   Z(3) = PHI0_LS(IP1,J)
   Z(4) = PHI0_LS(IP2,J)
   F_EAST = SCALAR_FACE_VALUE(SR_X_LS(I,J),Z,LIMITER_LS)

!    west face
   Z(1) = PHI0_LS(IM2,J)
   Z(2) = PHI0_LS(IM1,J)
   Z(3) = PHI0_LS(I,J)
   Z(4) = PHI0_LS(IP1,J)
   F_WEST = SCALAR_FACE_VALUE(SR_X_LS(I,J),Z,LIMITER_LS)

   IF (J == NY_LS-1) THEN
!    north face
    Z(1) = PHI0_LS(I,JM1)
    Z(2) = PHI0_LS(I,J)
    Z(3) = PHI0_LS(I,JP1)
    Z(4) = PHI_MIN_LS
    F_NORTH = SCALAR_FACE_VALUE(SR_Y_LS(I,J),Z,LIMITER_LS)

!    south face
    Z(1) = PHI0_LS(I,JM2)
    Z(2) = PHI0_LS(I,JM1)
    Z(3) = PHI0_LS(I,J)
    Z(4) = PHI0_LS(I,JP1)
    F_SOUTH = SCALAR_FACE_VALUE(SR_Y_LS(I,J),Z,LIMITER_LS)

   ELSEIF (J == NY_LS) THEN
!    north face
    Z(1) = PHI0_LS(I,JM1)
    Z(2) = PHI0_LS(I,J)
    Z(3) = PHI_MIN_LS
    Z(4) = PHI_MIN_LS
    F_NORTH = SCALAR_FACE_VALUE(SR_Y_LS(I,J),Z,LIMITER_LS)

!    south face
    Z(1) = PHI0_LS(I,JM2)
    Z(2) = PHI0_LS(I,JM1)
    Z(3) = PHI0_LS(I,J)
    Z(4) = PHI_MIN_LS
    F_SOUTH = SCALAR_FACE_VALUE(SR_Y_LS(I,J),Z,LIMITER_LS)
    ENDIF

   DPHIDX = (F_EAST-F_WEST)/DX_LS
   DPHIDY = (F_NORTH-F_SOUTH)/DY_LS
   FLUX_LS(I,J) = SR_X_LS(I,J)*DPHIDX + SR_Y_LS(I,J)*DPHIDY
   
   PHIMAG          = SQRT(DPHIDX**2 + DPHIDY**2)
   MAG_SR_OUT(I,J) = 0.0_EB
   IF(PHIMAG > 0.0_EB) MAG_SR_OUT(I,J) = FLUX_LS(I,J)/PHIMAG

!       fx = (f_east-f_west)/dx
!       fy = (f_north-f_south)/dy
        
!       phi(i,j) = phi0(i,j) - dt*[Fx(i,j) Fy(i,j)]*[fx fy]
 ENDDO

ENDDO ILOOP

IF (RK2_PREDICTOR_LS) FLUX0_LS = FLUX_LS
IF (.NOT. RK2_PREDICTOR_LS) FLUX1_LS = FLUX_LS


END SUBROUTINE LEVEL_SET_ADVECT_FLUX 
!
! ----------------------------------------------------
REAL(EB) FUNCTION SCALAR_FACE_VALUE(SR_XY,Z,LIMITER)
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
  r = 0._EB
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

SCALAR_FACE_VALUE = ZUP + 0.5*B*( ZDWN - ZUP )

END FUNCTION SCALAR_FACE_VALUE
!

SUBROUTINE GET_REV_vege(MODULE_REV,MODULE_DATE)
INTEGER,INTENT(INOUT) :: MODULE_REV
CHARACTER(255),INTENT(INOUT) :: MODULE_DATE

WRITE(MODULE_DATE,'(A)') vegerev(INDEX(vegerev,':')+2:LEN_TRIM(vegerev)-2)
READ (MODULE_DATE,'(I5)') MODULE_REV
WRITE(MODULE_DATE,'(A)') vegedate

END SUBROUTINE GET_REV_vege


END MODULE VEGE
