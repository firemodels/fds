!> \brief Level set model of fire spread across terrain

MODULE VEGE

USE COMP_FUNCTIONS, ONLY: CURRENT_TIME
USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_POINTERS
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
IMPLICIT NONE
PRIVATE
PUBLIC INITIALIZE_LEVEL_SET_FIRESPREAD_1,INITIALIZE_LEVEL_SET_FIRESPREAD_2,LEVEL_SET_FIRESPREAD
INTEGER :: IZERO
INTEGER  :: LIMITER_LS
REAL(EB) :: B_ROTH,BETA_OP_ROTH,C_ROTH,E_ROTH,T_NOW
REAL(EB), POINTER, DIMENSION(:,:) :: PHI_LS_P
REAL(EB), PARAMETER :: PHI_LS_MIN=-1._EB, PHI_LS_MAX=1._EB 

CONTAINS


!> \brief Set up the major level set arrays
!>
!> \param NM Mesh number
!> \details Set up the array for level set value PHI_LS, and determine terrain height on each 2D mesh, Z_LS(I,J).
!> After this routine, go back to main, exchange Z_LS, and return for more initialization.

SUBROUTINE INITIALIZE_LEVEL_SET_FIRESPREAD_1(NM)

USE COMPLEX_GEOMETRY, ONLY : IBM_IDCF
INTEGER, INTENT(IN) :: NM
INTEGER :: ICF,IW,I,J,SURF_INDEX
TYPE (MESH_TYPE),    POINTER :: M
TYPE (WALL_TYPE),    POINTER :: WC
TYPE (CFACE_TYPE),   POINTER :: CFA
TYPE (SURFACE_TYPE), POINTER :: SF

T_NOW = CURRENT_TIME()

CALL POINT_TO_MESH(NM)

M => MESHES(NM)

! Loop through all SURFace types and find level set cases that need a calculated RoS

DO SURF_INDEX=0,N_SURF
   SF => SURFACE(SURF_INDEX)
   IF (SF%VEG_LSET_SPREAD .AND. SF%VEG_LSET_FUEL_INDEX>0) THEN
      SF%VEG_LSET_ROS_00 = ROS_NO_WIND_NO_SLOPE(SF%VEG_LSET_FUEL_INDEX,SURF_INDEX)
   ENDIF
   IF (SF%VEG_LSET_SPREAD .AND. SF%VEG_LSET_FUEL_INDEX==0) THEN
     SF%BURN_DURATION = SF%VEG_LSET_FIREBASE_TIME
     IF(LEVEL_SET_COUPLED_FIRE) SF%MASS_FLUX(REACTION(1)%FUEL_SMIX_INDEX) = &
       (1._EB-SF%VEG_LSET_CHAR_FRACTION)*SF%VEG_LSET_SURF_LOAD/SF%VEG_LSET_FIREBASE_TIME
   ENDIF
ENDDO

! Level set values (Phi). PHI1_LS is the first-order accurate estimate at the next time step.

ALLOCATE(M%PHI_LS(0:IBP1,0:JBP1)) ; CALL ChkMemErr('VEGE:LEVEL SET','PHI_LS',IZERO)   ; PHI_LS  => M%PHI_LS  ; PHI_LS = PHI_LS_MIN
ALLOCATE(M%PHI1_LS(0:IBP1,0:JBP1)); CALL ChkMemErr('VEGE:LEVEL SET','PHI1_LS',IZERO)  ; PHI1_LS => M%PHI1_LS ; PHI1_LS = PHI_LS_MIN

! Wind speed components in the center of the first gas phsae cell above the ground.

ALLOCATE(M%U_LS(0:IBP1,0:JBP1)) ; CALL ChkMemErr('VEGE:LEVEL SET','U_LS',IZERO) ; U_LS => M%U_LS ; U_LS = 0._EB 
ALLOCATE(M%V_LS(0:IBP1,0:JBP1)) ; CALL ChkMemErr('VEGE:LEVEL SET','V_LS',IZERO) ; V_LS => M%V_LS ; V_LS = 0._EB

! Terrain height, Z_LS, and z index of the first gas cell above terrain, K_LS

ALLOCATE(M%Z_LS(0:IBP1,0:JBP1),STAT=IZERO) ; CALL ChkMemErr('READ','Z_LS',IZERO) ; Z_LS => M%Z_LS ; Z_LS = 0._EB
ALLOCATE(M%K_LS(0:IBP1,0:JBP1),STAT=IZERO) ; CALL ChkMemErr('READ','K_LS',IZERO) ; K_LS => M%K_LS ; K_LS = 0
ALLOCATE(M%LS_SURF_INDEX(0:IBP1,0:JBP1),STAT=IZERO) ; CALL ChkMemErr('READ','LS_SURF_INDEX',IZERO)  
LS_SURF_INDEX => M%LS_SURF_INDEX ; LS_SURF_INDEX = 0

IF (CC_IBM) THEN

   ALLOCATE(M%LS_KLO_TERRAIN(0:IBP1,0:JBP1),STAT=IZERO) ; CALL ChkMemErr('READ','LS_KLO_TERRAIN',IZERO)
   LS_KLO_TERRAIN => M%LS_KLO_TERRAIN ; LS_KLO_TERRAIN = 2*KBP1+1 ! Number larger that KBP1.
   DO ICF=1,M%N_CUTFACE_MESH
      IF (CUT_FACE(ICF)%STATUS /= 2) CYCLE ! IBM_INBOUNDARY == 2
      ! Location of CFACE with largest AREA, to define SURF_INDEX:
      IW  = MAXLOC(CUT_FACE(ICF)%AREA(1:CUT_FACE(ICF)%NFACE),DIM=1)
      CFA => CFACE( CUT_FACE(ICF)%CFACE_INDEX(IW) )
      IF (CFA%NVEC(KAXIS)>-TWO_EPSILON_EB .AND. CFA%BOUNDARY_TYPE==SOLID_BOUNDARY) THEN
         ! Area averaged Z height of CFACES within this cut-cell (containing IBM_INBOUNDARY CFACES):
         Z_LS(CFA%ONE_D%IIG,CFA%ONE_D%JJG) = DOT_PRODUCT(CUT_FACE(ICF)%XYZCEN(KAXIS,1:CUT_FACE(ICF)%NFACE), &
                                                                   CUT_FACE(ICF)%  AREA(1:CUT_FACE(ICF)%NFACE))     / &
                                                               SUM(CUT_FACE(ICF)% AREA(1:CUT_FACE(ICF)%NFACE))
         IF (CFA%ONE_D%KKG < LS_KLO_TERRAIN(CFA%ONE_D%IIG,CFA%ONE_D%JJG)) & 
            LS_KLO_TERRAIN(CFA%ONE_D%IIG,CFA%ONE_D%JJG) = CFA%ONE_D%KKG
         IF (CFA%ONE_D%KKG > K_LS(CFA%ONE_D%IIG,CFA%ONE_D%JJG)) K_LS(CFA%ONE_D%IIG,CFA%ONE_D%JJG) = CFA%ONE_D%KKG
         LS_SURF_INDEX(CFA%ONE_D%IIG,CFA%ONE_D%JJG) = CFA%SURF_INDEX
      ENDIF
   ENDDO
   DO J=1,JBAR
      DO I=1,IBAR
         IF (K_LS(I,J)==KBAR .AND. FCVAR(I,J,K_LS(I,J),IBM_IDCF,KAXIS)>0) LS_SURF_INDEX(I,J) = 0
      ENDDO
   ENDDO

ELSE

   DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
      WC => WALL(IW)
      IF (WC%ONE_D%IOR==3 .AND. WC%BOUNDARY_TYPE==SOLID_BOUNDARY) THEN
         Z_LS(WC%ONE_D%IIG,WC%ONE_D%JJG) = Z(WC%ONE_D%KKG-1)
         K_LS(WC%ONE_D%IIG,WC%ONE_D%JJG) = WC%ONE_D%KKG
         LS_SURF_INDEX(WC%ONE_D%IIG,WC%ONE_D%JJG)= WC%SURF_INDEX
      ENDIF
   ENDDO

ENDIF

Z_LS(1:IBAR,   0) = 2._EB*Z_LS(1:IBAR,   1) - Z_LS(1:IBAR,   2)
Z_LS(1:IBAR,JBP1) = 2._EB*Z_LS(1:IBAR,JBAR) - Z_LS(1:IBAR,JBM1)
Z_LS(   0,1:JBAR) = 2._EB*Z_LS(   1,1:JBAR) - Z_LS(   2,1:JBAR)
Z_LS(IBP1,1:JBAR) = 2._EB*Z_LS(IBAR,1:JBAR) - Z_LS(IBM1,1:JBAR)

Z_LS(   0,   0) = Z_LS(   1,   1)
Z_LS(IBP1,   0) = Z_LS(IBAR,   1)
Z_LS(   0,JBP1) = Z_LS(   1,JBAR)
Z_LS(IBP1,JBP1) = Z_LS(IBAR,JBAR)

T_USED(15) = T_USED(15) + CURRENT_TIME() - T_NOW
END SUBROUTINE INITIALIZE_LEVEL_SET_FIRESPREAD_1


!> \brief Continue initialialization of level set routines
!>
!> \param NM Mesh number 
!> \details First, retrieve terrain height, Z_LS, from other meshes. Then do various other set up chores.

SUBROUTINE INITIALIZE_LEVEL_SET_FIRESPREAD_2(NM)

INTEGER, INTENT(IN) :: NM
INTEGER :: I,IM1,IM2,IIG,IP1,IP2,J,JJG,JM1,JP1
REAL(EB) :: DZT_DUM,G_EAST,G_WEST,G_SOUTH,G_NORTH
TYPE (MESH_TYPE),    POINTER :: M
TYPE (SURFACE_TYPE), POINTER :: SF

T_NOW = CURRENT_TIME()

CALL POINT_TO_MESH(NM)

M => MESHES(NM)

! Retrieve terrain height, Z_LS, from other meshes above and below current mesh.

CALL GET_BOUNDARY_VALUES

! Allocate some work arrays

ALLOCATE(M%LS_WORK1(0:IBAR,0:JBAR))    ; CALL ChkMemErr('VEGE:LEVEL SET','LS_WORK1',IZERO)
ALLOCATE(M%LS_WORK2(0:IBAR,0:JBAR))    ; CALL ChkMemErr('VEGE:LEVEL SET','LS_WORK2',IZERO)

! Define spread rate across domain (including no burn areas)

ALLOCATE(M%ROS_HEAD(IBAR,JBAR))    ; CALL ChkMemErr('VEGE:LEVEL SET','ROS_HEAD',IZERO)  ; ROS_HEAD => M%ROS_HEAD
ALLOCATE(M%ROS_FLANK(IBAR,JBAR))   ; CALL ChkMemErr('VEGE:LEVEL SET','ROS_FLANK',IZERO) ; ROS_FLANK => M%ROS_FLANK
ALLOCATE(M%ROS_BACKU(IBAR,JBAR))   ; CALL ChkMemErr('VEGE:LEVEL SET','ROS_BACKU',IZERO) ; ROS_BACKU => M%ROS_BACKU
ALLOCATE(M%WIND_EXP(IBAR,JBAR))    ; CALL ChkMemErr('VEGE:LEVEL SET','WIND_EXP',IZERO)  ; WIND_EXP => M%WIND_EXP

! Assign spread rates (i.e., vegetation types) to locations on terrain

ROS_HEAD  = 0.0_EB
ROS_FLANK = 0.0_EB
ROS_BACKU = 0.0_EB
WIND_EXP  = 1.0_EB

LSET_TAN2    = .FALSE. ! Flag for ROS proportional to Tan(slope)^2

! Flux limiters
! LIMITER_LS=1 MINMOD
! LIMITER_LS=2 SUPERBEE
! LIMITER_LS=3 First order upwinding

LIMITER_LS = I_FLUX_LIMITER
IF (LIMITER_LS > 3) LIMITER_LS = 1

! Flux terms

ALLOCATE(M%FLUX0_LS(IBAR,JBAR)); CALL ChkMemErr('VEGE:LEVEL SET','FLUX0_LS',IZERO) ; FLUX0_LS => M%FLUX0_LS
ALLOCATE(M%FLUX1_LS(IBAR,JBAR)); CALL ChkMemErr('VEGE:LEVEL SET','FLUX1_LS',IZERO) ; FLUX1_LS => M%FLUX1_LS

! Slopes (gradients)

ALLOCATE(M%DZTDX(IBAR,JBAR)); CALL ChkMemErr('VEGE:LEVEL SET','DZDTX',IZERO)   ; DZTDX => M%DZTDX
ALLOCATE(M%DZTDY(IBAR,JBAR)); CALL ChkMemErr('VEGE:LEVEL SET','DZDTY',IZERO)   ; DZTDY => M%DZTDY
ALLOCATE(M%MAG_ZT(IBAR,JBAR)); CALL ChkMemErr('VEGE:LEVEL SET','MAG_ZT',IZERO) ; MAG_ZT => M%MAG_ZT

! Rothermel 'Phi' factors for effects of Wind and Slope on ROS

ALLOCATE(M%PHI_WS(IBAR,JBAR))   ; CALL ChkMemErr('VEGE:LEVEL SET','PHI_W',IZERO)   ; PHI_WS => M%PHI_WS    ; PHI_WS = 0.0_EB
ALLOCATE(M%PHI_S(IBAR,JBAR))    ; CALL ChkMemErr('VEGE:LEVEL SET','PHI_S',IZERO)   ; PHI_S => M%PHI_S
ALLOCATE(M%PHI_S_X(IBAR,JBAR))  ; CALL ChkMemErr('VEGE:LEVEL SET','PHI_S_X',IZERO) ; PHI_S_X => M%PHI_S_X
ALLOCATE(M%PHI_S_Y(IBAR,JBAR))  ; CALL ChkMemErr('VEGE:LEVEL SET','PHI_S_Y',IZERO) ; PHI_S_Y => M%PHI_S_Y

! UMF = wind speed at mid-flame height

ALLOCATE(M%UMF(IBAR,JBAR))    ; CALL ChkMemErr('VEGE:LEVEL SET','UMF',IZERO) ; M%UMF = 0._EB ; UMF => M%UMF
ALLOCATE(M%THETA_ELPS(IBAR,JBAR))    ; CALL ChkMemErr('VEGE:LEVEL SET','THETA_ELPS',IZERO) ; THETA_ELPS => M%THETA_ELPS
THETA_ELPS = 0.0_EB ! Normal to fireline

! ROS in X and Y directions

ALLOCATE(M%SR_X_LS(IBAR,JBAR)) ; CALL ChkMemErr('VEGE:LEVEL SET','SR_X_LS',IZERO) ; SR_X_LS => M%SR_X_LS
ALLOCATE(M%SR_Y_LS(IBAR,JBAR)) ; CALL ChkMemErr('VEGE:LEVEL SET','SR_Y_LS',IZERO) ; SR_Y_LS => M%SR_Y_LS

! Compute components of terrain slope gradient and magnitude of gradient

GRADIENT_ILOOP: DO I = 1,IBAR

   IM1=I-1 ; IM2=I-2
   IP1=I+1 ; IP2=I+2
   IF (I==1) IM1 = I
   IF (I==IBAR) IP1 = I

   DO J = 1,JBAR

      JM1=J-1
      JP1=J+1
      IF (J==1) JM1 = J
      IF (J==IBAR) JP1 = J

      G_EAST  = 0.5_EB*( Z_LS(I,J) + Z_LS(IP1,J) )
      G_WEST  = 0.5_EB*( Z_LS(I,J) + Z_LS(IM1,J) )
      G_NORTH = 0.5_EB*( Z_LS(I,J) + Z_LS(I,JP1) )
      G_SOUTH = 0.5_EB*( Z_LS(I,J) + Z_LS(I,JM1) )

      DZTDX(I,J) = (G_EAST-G_WEST)   * RDX(I)
      DZTDY(I,J) = (G_NORTH-G_SOUTH) * RDY(J)
      MAG_ZT(I,J) = SQRT(DZTDX(I,J)**2 + DZTDY(I,J)**2)

   ENDDO

ENDDO GRADIENT_ILOOP

! Initialize arrays for head, flank, and back fire spread rates with values explicitly declared in the input file or 
! from FARSITE head fire and ellipse based flank and back fires.

DO JJG=1,JBAR
   DO IIG=1,IBAR

      SF => SURFACE(LS_SURF_INDEX(IIG,JJG))

      IF (.NOT. SF%VEG_LSET_SPREAD) CYCLE

      ! Initialize various arrays

      ROS_HEAD(IIG,JJG)  = SF%VEG_LSET_ROS_HEAD
      ROS_FLANK(IIG,JJG) = SF%VEG_LSET_ROS_FLANK
      ROS_BACKU(IIG,JJG) = SF%VEG_LSET_ROS_BACK
      WIND_EXP(IIG,JJG)  = SF%VEG_LSET_WIND_EXP

      ! If any surfaces uses tan^2 function for slope, tan^2 will be used throughout simulation

      IF (SF%VEG_LSET_TAN2) LSET_TAN2=.TRUE.

      ! Use assumed ellipse shape of fireline as in Farsite

      IF_ELLIPSE: IF (LEVEL_SET_ELLIPSE) THEN

         SF%VEG_LSET_HT = MAX(0.001_EB,SF%VEG_LSET_HT)

         ! Variables used in Phi_W slope factor formulas below (Rothermel model)

         B_ROTH = 0.15988_EB * (SF%VEG_LSET_SIGMA**0.54_EB)
         C_ROTH = 7.47_EB * EXP(-0.8711_EB * (SF%VEG_LSET_SIGMA**0.55_EB))
         E_ROTH = 0.715_EB * EXP(-0.01094_EB * SF%VEG_LSET_SIGMA)
         BETA_OP_ROTH = 0.20395_EB * (SF%VEG_LSET_SIGMA**(-0.8189_EB))! Optimum packing ratio

         ! Limit effect to slope lte 80 degrees. Phi_s_x,y are slope factors (Rothermel model)

         DZT_DUM = MIN(5.67_EB,ABS(DZTDX(IIG,JJG))) ! 5.67 ~ tan 80 deg
         PHI_S_X(IIG,JJG) = 5.275_EB * ((SF%VEG_LSET_BETA)**(-0.3_EB)) * DZT_DUM**2
         PHI_S_X(IIG,JJG) = SIGN(PHI_S_X(IIG,JJG),DZTDX(IIG,JJG))
         DZT_DUM = MIN(1.73_EB,ABS(DZTDY(IIG,JJG))) ! 1.73 ~ tan 60 deg
         PHI_S_Y(IIG,JJG) = 5.275_EB * ((SF%VEG_LSET_BETA)**(-0.3_EB)) * DZT_DUM**2
         PHI_S_Y(IIG,JJG) = SIGN(PHI_S_Y(IIG,JJG),DZTDY(IIG,JJG))

         PHI_S(IIG,JJG) = SQRT(PHI_S_X(IIG,JJG)**2 + PHI_S_Y(IIG,JJG)**2)

      ENDIF IF_ELLIPSE

   ENDDO
ENDDO

T_USED(15) = T_USED(15) + CURRENT_TIME() - T_NOW
END SUBROUTINE INITIALIZE_LEVEL_SET_FIRESPREAD_2


!> \brief Advance the level set array one time step
!>
!> \param T Current time (s)
!> \param DT Time step (s)
!> \param NM Mesh number
!> \details Predictor: Estimate PHI_LS at next time step. Estimated value is called PHI1_LS.
!> Corrector: Correct PHI_LS at next time step.

SUBROUTINE LEVEL_SET_FIRESPREAD(T,DT,NM)

USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: T,DT
INTEGER :: IIG,IW,JJG,IC
INTEGER :: KDUM,KWIND,ICF,IKT
REAL(EB) :: UMF_TMP,PHX,PHY,MAG_PHI,PHI_W_X,PHI_W_Y,UMF_X,UMF_Y,UMAG,DUMMY=0._EB,ROS_MAG,UMF_MAG
TYPE (ONE_D_M_AND_E_XFER_TYPE), POINTER :: ONE_D
TYPE (SURFACE_TYPE), POINTER :: SF

T_NOW = CURRENT_TIME()

CALL POINT_TO_MESH(NM)

CALL GET_BOUNDARY_VALUES

! Loop over terrain surface cells and update level set field

DO JJG=1,JBAR
   DO IIG=1,IBAR

      SF => SURFACE(LS_SURF_INDEX(IIG,JJG))

      IF (.NOT. SF%VEG_LSET_SPREAD) CYCLE

      ! Ignite landscape at user specified location and time

      IF (PHI_LS(IIG,JJG)<0._EB .AND. T>SF%VEG_LSET_IGNITE_T) THEN
         PHI_LS(IIG,JJG) = PHI_LS_MAX
         IF (LEVEL_SET_MODE==2) THEN
            FREEZE_VELOCITY = .TRUE.
            SOLID_PHASE_ONLY = .TRUE.
         ENDIF
      ENDIF

      ! Establish the wind field

      IF_CFD_COUPLED: IF (LEVEL_SET_COUPLED_WIND) THEN  ! The wind speed is derived from the CFD computation

         U_LS(IIG,JJG) = 0.5_EB*(U(IIG-1,JJG,K_LS(IIG,JJG))+U(IIG,JJG,K_LS(IIG,JJG)))
         V_LS(IIG,JJG) = 0.5_EB*(V(IIG,JJG-1,K_LS(IIG,JJG))+V(IIG,JJG,K_LS(IIG,JJG)))

      ELSE IF_CFD_COUPLED  ! The wind velocity is specified by the user
 
         U_LS(IIG,JJG) = U0*EVALUATE_RAMP(T,DUMMY,I_RAMP_U0_T)
         V_LS(IIG,JJG) = V0*EVALUATE_RAMP(T,DUMMY,I_RAMP_V0_T)

      ENDIF IF_CFD_COUPLED

      IF_ELLIPSE: IF (LEVEL_SET_ELLIPSE) THEN  ! Use assumed elliptical shape of fireline as in Farsite

         ! Find wind at ~6.1 m height for Farsite

         IF (LEVEL_SET_COUPLED_WIND) THEN

            KWIND = 0
            DO KDUM = K_LS(IIG,JJG),KBAR
               IF (ZC(KDUM)-ZC(K_LS(IIG,JJG))>=6.1_EB) THEN
                  KWIND = KDUM
                  EXIT
               ENDIF
            ENDDO

            U_LS(IIG,JJG) = 0.5_EB*(U(IIG-1,JJG,KWIND)+U(IIG,JJG,KWIND))
            V_LS(IIG,JJG) = 0.5_EB*(V(IIG,JJG-1,KWIND)+V(IIG,JJG,KWIND))

         ENDIF 

         IF (LEVEL_SET_MODE == 5) THEN
           U_LS(IIG,JJG) = U0
           V_LS(IIG,JJG) = V0
         ENDIF

         ! Wind at midflame height (UMF). From Andrews 2012, USDA FS Gen Tech Rep. RMRS-GTR-266 (with added SI conversion)

         UMF_TMP = 1.83_EB / LOG((20.0_EB + 1.18_EB * SF%VEG_LSET_HT) /(0.43_EB * SF%VEG_LSET_HT))

         ! Factor 60 converts U from m/s to m/min which is used in elliptical model.

         UMF_X = UMF_TMP * U_LS(IIG,JJG) * 60.0_EB
         UMF_Y = UMF_TMP * V_LS(IIG,JJG) * 60.0_EB
         UMF_MAG = SQRT(UMF_X**2 + UMF_Y**2)

         ! Components of Rothermel wind factor - affects spread rate

         PHI_W_X = C_ROTH * ((3.281_EB * UMF_MAG)**B_ROTH) * (SF%VEG_LSET_BETA / BETA_OP_ROTH)**(-E_ROTH)
         PHI_W_X = PHI_W_X*UMF_X/UMF_MAG

         PHI_W_Y = C_ROTH * ((3.281_EB * UMF_MAG)**B_ROTH) * (SF%VEG_LSET_BETA / BETA_OP_ROTH)**(-E_ROTH)
         PHI_W_Y = PHI_W_Y*UMF_Y/UMF_MAG

         ! Include Rothermel slope factor

         IF (PHI_S(IIG,JJG) > 0.0_EB) THEN

            PHX = PHI_W_X + PHI_S_X(IIG,JJG)
            PHY = PHI_W_Y + PHI_S_Y(IIG,JJG)
            MAG_PHI = SQRT(PHX**2 + PHY**2)

            ! Total phi (phi_w + phi_s) for use in spread rate section

            PHI_WS(IIG,JJG) = MAG_PHI

            ! Theta_elps is angle of direction (0 to 2pi) of highest spread rate
            ! 0<=theta_elps<=2pi as measured clockwise from Y-axis

            THETA_ELPS(IIG,JJG) = ATAN2(PHY,PHX)

            ! "Effective midflame windspeed" used in length-to-breadth ratio calculation (spread rate routine)
            ! is the wind + slope effect obtained by solving Phi_w eqs. above for UMF
            ! 8/8/13 - Changed phi_ws to Phi_s below to match Farsite, i.e., instead of adding phi_w and phi_s
            ! and then calculating effective wind speed, phi_s is converted to an effected wind speed and added
            ! to UMF calculated from the wind. Effective U has units of m/min in Wilson formula.
            ! 0.3048 ~= 1/3.281
            ! if phi_s < 0 then a complex value (NaN) results. Using abs(phi_s) and sign function to correct.

            UMF_TMP = (((ABS(PHI_S_X(IIG,JJG)) * (SF%VEG_LSET_BETA / BETA_OP_ROTH)**E_ROTH)/C_ROTH)**(1/B_ROTH))*0.3048
            UMF_TMP = SIGN(UMF_TMP,PHI_S_X(IIG,JJG))
            UMF_X = UMF_X + UMF_TMP

            UMF_TMP = (((ABS(PHI_S_Y(IIG,JJG)) * (SF%VEG_LSET_BETA / BETA_OP_ROTH)**E_ROTH)/C_ROTH)**(1/B_ROTH))*0.3048
            UMF_TMP = SIGN(UMF_TMP,PHI_S_Y(IIG,JJG))
            UMF_Y = UMF_Y + UMF_TMP

         ELSE

            PHI_WS(IIG,JJG) = SQRT(PHI_W_X**2 + PHI_W_Y**2)
            THETA_ELPS(IIG,JJG) = ATAN2(PHI_W_Y,PHI_W_X)

         ENDIF

         UMF(IIG,JJG) = SQRT(UMF_X**2 + UMF_Y**2)

         ! The following two lines convert ATAN2 output to compass system (0 to 2 pi CW from +Y-axis)

         THETA_ELPS(IIG,JJG) = PIO2 - THETA_ELPS(IIG,JJG)
         IF (THETA_ELPS(IIG,JJG) < 0.0_EB) THETA_ELPS(IIG,JJG) = 2.0_EB*PI + THETA_ELPS(IIG,JJG)

      ELSE IF_ELLIPSE  ! AU grassland ROS for infinite head and 6% moisutre

         UMAG     = SQRT(U_LS(IIG,JJG)**2 + V_LS(IIG,JJG)**2)
         ROS_HEAD(IIG,JJG)  = SF%VEG_LSET_ROS_HEAD*(0.165_EB + 0.534_EB*UMAG)*0.523_EB

      ENDIF IF_ELLIPSE

      IF (SF%VEG_LSET_ROS_00 > 0._EB) ROS_HEAD(IIG,JJG) = SF%VEG_LSET_ROS_00*(1._EB + PHI_WS(IIG,JJG))

   ENDDO
ENDDO

! Runge-Kutta Scheme

CALL LEVEL_SET_SPREAD_RATE
CALL LEVEL_SET_ADVECT_FLUX

IF (PREDICTOR) THEN
   PHI1_LS(1:IBAR,1:JBAR) = PHI_LS(1:IBAR,1:JBAR) - DT*FLUX0_LS(1:IBAR,1:JBAR)
   PHI1_LS = MAX(PHI_LS_MIN,MIN(PHI_LS_MAX,PHI1_LS))
ELSE
   PHI_LS(1:IBAR,1:JBAR) = PHI_LS(1:IBAR,1:JBAR) - 0.5_EB*DT*(FLUX0_LS(1:IBAR,1:JBAR) + FLUX1_LS(1:IBAR,1:JBAR))
   PHI_LS = MAX(PHI_LS_MIN,MIN(PHI_LS_MAX,PHI_LS))
ENDIF

! Loop over all cells and assign the value of PHI_LS to the appropriate WALL or
! CFACE cells. Also, if PHI_LS increases above 0, set the ignition time T_IGN.

IF (.NOT.PREDICTOR) THEN
   IF (CC_IBM) THEN
      DO JJG=1,JBAR
         DO IIG=1,IBAR
            SF => SURFACE(LS_SURF_INDEX(IIG,JJG))
            IF (.NOT. SF%VEG_LSET_SPREAD) CYCLE
            DO IKT=LS_KLO_TERRAIN(IIG,JJG),K_LS(IIG,JJG)
               ! Loop over all CFACEs corresponding to IIG,JJG and set ONE_D%T_IGN and ONE_D%PHI_LS as below
               ICF = CCVAR(IIG,JJG,IKT,3); IF(ICF<1) CYCLE  ! IBM_IDCF = 3 CUT_FCE container for this cell.
               DO IW=1,CUT_FACE(ICF)%NFACE ! All IBM_INBOUNDARY CFACES on this cell.
                  ONE_D => CFACE( CUT_FACE(ICF)%CFACE_INDEX(IW) ) % ONE_D
                  IF (PHI_LS(IIG,JJG)>=0._EB .AND. ONE_D%T_IGN>1.E5_EB) CALL IGNITE_GRID_CELL
                  ONE_D%PHI_LS = PHI_LS(IIG,JJG)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ELSE
      DO JJG=1,JBAR
         DO IIG=1,IBAR
            SF => SURFACE(LS_SURF_INDEX(IIG,JJG))
            IF (.NOT. SF%VEG_LSET_SPREAD) CYCLE
            IF (K_LS(IIG,JJG)<1) CYCLE
            IC = CELL_INDEX(IIG,JJG,K_LS(IIG,JJG))
            IW = WALL_INDEX(IC,-3)
            ONE_D => WALL(IW)%ONE_D
            IF (PHI_LS(IIG,JJG)>=0._EB .AND. ONE_D%T_IGN>1.E5_EB) CALL IGNITE_GRID_CELL
            ONE_D%PHI_LS = PHI_LS(IIG,JJG)
         ENDDO
      ENDDO
   ENDIF
ENDIF

T_USED(15) = T_USED(15) + CURRENT_TIME() - T_NOW

CONTAINS


!> \brief Set the time, burning duration, and burning rate of a newly ignited grid cell

SUBROUTINE IGNITE_GRID_CELL

REAL(EB) :: CROSSING_DISTANCE

ONE_D%T_IGN = T
ROS_MAG = MAX(0.01_EB,SQRT(SR_X_LS(IIG,JJG)**2 + SR_Y_LS(IIG,JJG)**2))  ! Rate Of Spread magnitude
IF (ABS(SR_X_LS(IIG,JJG))<TWO_EPSILON_EB) THEN
   CROSSING_DISTANCE = DY(JJG)
ELSEIF (ABS(SR_Y_LS(IIG,JJG))<TWO_EPSILON_EB) THEN
   CROSSING_DISTANCE = DX(IIG)
ELSE
   CROSSING_DISTANCE = SQRT( (MIN(DX(IIG),ABS(SR_X_LS(IIG,JJG)/SR_Y_LS(IIG,JJG))*DY(JJG)))**2 + &
                             (MIN(DY(JJG),ABS(SR_Y_LS(IIG,JJG)/SR_X_LS(IIG,JJG))*DX(IIG)))**2 )
ENDIF
ONE_D%BURN_DURATION = SF%BURN_DURATION + CROSSING_DISTANCE/ROS_MAG
ONE_D%AREA_ADJUST = SF%BURN_DURATION/ONE_D%BURN_DURATION

END SUBROUTINE IGNITE_GRID_CELL

END SUBROUTINE LEVEL_SET_FIRESPREAD


!> \brief Retrieve various quantities from neighboring meshes after MPI exchange

SUBROUTINE GET_BOUNDARY_VALUES

INTEGER :: IIG,JJG,II,JJ,IOR

IF (PREDICTOR) THEN
   PHI_LS_P => PHI_LS
ELSE
   PHI_LS_P => PHI1_LS
ENDIF

! Fetch values of various quantities at six faces of current mesh.

DO II=1,IBAR
   IIG=II ; JJ=0 ; JJG=1 ; IOR=-2
   CALL FILL_BOUNDARY_VALUES
   IIG=II ; JJ=JBP1 ; JJG=JBAR ; IOR=2
   CALL FILL_BOUNDARY_VALUES
ENDDO
DO JJ=1,JBAR
   JJG=JJ ; II=0 ; IIG=1 ; IOR=-1
   CALL FILL_BOUNDARY_VALUES
   JJG=JJ ; II=IBP1 ; IIG=IBAR ; IOR=1
   CALL FILL_BOUNDARY_VALUES
ENDDO
DO JJ=1,JBAR
   JJG=JJ
   DO II=1,IBAR
      IIG=II
      IOR = -3
      CALL FILL_BOUNDARY_VALUES
      IOR =  3
      CALL FILL_BOUNDARY_VALUES
   ENDDO
ENDDO

CONTAINS

! \brief Grab boundary values of PHI_LS from MPI storage arrays

SUBROUTINE FILL_BOUNDARY_VALUES

USE COMPLEX_GEOMETRY, ONLY : IBM_CGSC,IBM_SOLID,IBM_CUTCFE
INTEGER :: IW,IIO,JJO,N_INT_CELLS,NOM,IC
REAL(EB) :: PHI_LS_OTHER,U_LS_OTHER,V_LS_OTHER,Z_LS_OTHER
TYPE (EXTERNAL_WALL_TYPE), POINTER :: EWC
LOGICAL :: SOLID_CELL

IF (IOR==3) THEN  ! get the CELL_INDEX of the grid cell adjacent to the exterior boundary of the current mesh
   IC = CELL_INDEX(IIG,JJG,KBAR)
ELSEIF (IOR==-3) THEN
   IC = CELL_INDEX(IIG,JJG,1)
ELSE
   IC = CELL_INDEX(IIG,JJG,1)
ENDIF
IW = WALL_INDEX(IC,IOR)
EWC=>EXTERNAL_WALL(IW)
NOM = EWC%NOM
IF (NOM==0) RETURN  ! there is no other mesh adjacent to the boundary

PHI_LS_OTHER = 0._EB
U_LS_OTHER = 0._EB
V_LS_OTHER = 0._EB
Z_LS_OTHER = 0._EB
DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
   DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
      IF (PREDICTOR) THEN
         PHI_LS_OTHER = PHI_LS_OTHER + OMESH(NOM)%PHI_LS(IIO,JJO)
      ELSE
         PHI_LS_OTHER = PHI_LS_OTHER + OMESH(NOM)%PHI1_LS(IIO,JJO)
      ENDIF
      U_LS_OTHER = U_LS_OTHER + OMESH(NOM)%U_LS(IIO,JJO)
      V_LS_OTHER = V_LS_OTHER + OMESH(NOM)%V_LS(IIO,JJO)
      Z_LS_OTHER = Z_LS_OTHER + OMESH(NOM)%Z_LS(IIO,JJO)
   ENDDO
ENDDO
N_INT_CELLS = (EWC%IIO_MAX-EWC%IIO_MIN+1) * (EWC%JJO_MAX-EWC%JJO_MIN+1)

SELECT CASE(IOR)
   CASE(-2:2) 
      PHI_LS_P(II,JJ) = PHI_LS_OTHER/REAL(N_INT_CELLS,EB)
      U_LS(II,JJ)     = U_LS_OTHER/REAL(N_INT_CELLS,EB)
      V_LS(II,JJ)     = V_LS_OTHER/REAL(N_INT_CELLS,EB)
      Z_LS(II,JJ)     = Z_LS_OTHER/REAL(N_INT_CELLS,EB)
   CASE(3)  ! only grab a PHI_LS value from the other mesh if the (II,JJ) cell of the current mesh has no terrain surface
      SOLID_CELL = .FALSE.
      IF (CC_IBM) THEN
         IF (CCVAR(II,JJ,KBAR,IBM_CGSC)==IBM_SOLID .OR. CCVAR(II,JJ,KBAR,IBM_CGSC)==IBM_CUTCFE) SOLID_CELL = .TRUE.
      ELSE
         IF (SOLID(CELL_INDEX(II,JJ,KBAR))) SOLID_CELL = .TRUE.
      ENDIF
      IF (.NOT.SURFACE(LS_SURF_INDEX(II,JJ))%VEG_LSET_SPREAD .AND. SOLID_CELL) THEN
         PHI_LS_P(II,JJ) = PHI_LS_OTHER/REAL(N_INT_CELLS,EB)
         U_LS(II,JJ) = U_LS_OTHER/REAL(N_INT_CELLS,EB)
         V_LS(II,JJ) = V_LS_OTHER/REAL(N_INT_CELLS,EB)
         Z_LS(II,JJ) = Z_LS_OTHER/REAL(N_INT_CELLS,EB)
      ENDIF
   CASE(-3)  ! only grab a PHI_LS value from the other mesh if the (II,JJ) cell of the current mesh has no terrain surface
      SOLID_CELL = .FALSE.
      IF (CC_IBM) THEN
         IF (CCVAR(II,JJ,1,IBM_CGSC)==IBM_SOLID .OR. CCVAR(II,JJ,1,IBM_CGSC)==IBM_CUTCFE) SOLID_CELL = .TRUE.
      ELSE
         IF (SOLID(CELL_INDEX(II,JJ,1))) SOLID_CELL = .TRUE.
      ENDIF
      IF (.NOT.SURFACE(LS_SURF_INDEX(II,JJ))%VEG_LSET_SPREAD .AND. .NOT.SOLID_CELL) THEN
         PHI_LS_P(II,JJ) = PHI_LS_OTHER/REAL(N_INT_CELLS,EB)
         U_LS(II,JJ) = U_LS_OTHER/REAL(N_INT_CELLS,EB)
         V_LS(II,JJ) = V_LS_OTHER/REAL(N_INT_CELLS,EB)
         Z_LS(II,JJ) = Z_LS_OTHER/REAL(N_INT_CELLS,EB)
      ENDIF
END SELECT

END SUBROUTINE FILL_BOUNDARY_VALUES

END SUBROUTINE GET_BOUNDARY_VALUES


!> \brief Compute components of spread rate vector

SUBROUTINE LEVEL_SET_SPREAD_RATE

INTEGER :: I,J,IM1,IP1,JM1,JP1
REAL(EB) :: COS_THETA_WIND,COS_THETA_SLOPE,COS_THETA_WIND_H,COS_THETA_WIND_B, &
            COS_THETA_SLOPE_H,COS_THETA_SLOPE_B,DPHIDX,DPHIDY,F_EAST,F_WEST,F_NORTH,F_SOUTH, &
            GRAD_SLOPE_DOT_NORMAL_FIRELINE,MAG_F,MAG_SR,MAG_U,WIND_DOT_NORMAL_FIRELINE,NEXP_WIND
REAL(EB) :: ROS_BACKS,ROS_HEADS
REAL(EB) :: RAD_TO_DEGREE,DEGREES_SLOPE,SLOPE_FACTOR
REAL(EB) :: COS_THETA,SIN_THETA,XSF,YSF,UMF_DUM
REAL(EB) :: AROS,A_ELPS,A_ELPS2,BROS,B_ELPS2,B_ELPS,C_ELPS,DENOM,ROS_TMP,LB,LBD,HB
REAL(EB), DIMENSION(:) :: NORMAL_FIRELINE(2)

RAD_TO_DEGREE = 90._EB/ASIN(1._EB)

IF (PREDICTOR) THEN
   PHI_LS_P => PHI_LS
ELSE
   PHI_LS_P => PHI1_LS
ENDIF

SR_X_LS = 0.0_EB ; SR_Y_LS = 0.0_EB

FLUX_ILOOP: DO J=1,JBAR

   JM1 = J-1
   JP1 = J+1

   DO I=1,IBAR

      IM1 = I-1
      IP1 = I+1

      F_EAST  = 0.5_EB*( PHI_LS_P(I,J) + PHI_LS_P(IP1,J) )
      F_WEST  = 0.5_EB*( PHI_LS_P(I,J) + PHI_LS_P(IM1,J) )
      F_NORTH = 0.5_EB*( PHI_LS_P(I,J) + PHI_LS_P(I,JP1) )
      F_SOUTH = 0.5_EB*( PHI_LS_P(I,J) + PHI_LS_P(I,JM1) )

      DPHIDX = (F_EAST-F_WEST)   * RDX(I)
      DPHIDY = (F_NORTH-F_SOUTH) * RDY(J)

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

      IF (MAG_ZT(I,J) > 0.0_EB) COS_THETA_SLOPE = GRAD_SLOPE_DOT_NORMAL_FIRELINE/MAG_ZT(I,J)

      DEGREES_SLOPE = ATAN(MAG_ZT(I,J))*RAD_TO_DEGREE

      IF (LEVEL_SET_ELLIPSE) THEN

         ! Effective wind direction (theta) is clockwise from y-axis (Richards 1990)
         COS_THETA = COS(THETA_ELPS(I,J)) !V_LS(I,J) / MAG_U
         SIN_THETA = SIN(THETA_ELPS(I,J)) !U_LS(I,J) / MAG_U

         ROS_TMP = ROS_HEAD(I,J)

         ! Mag of wind speed at midflame ht must be in units of m/s here
         UMF_DUM = UMF(I,J)/60.0_EB

         ! Length to breadth ratio of ellipse based on effective UMF
         LB = 0.936_EB * EXP(0.2566_EB * UMF_DUM) + 0.461_EB * EXP(-0.1548_EB * UMF_DUM) - 0.397_EB

         ! Constraint LB max = 8 from Finney 2004
         LB = MAX(1.0_EB,MIN(LB,8.0_EB))
         LBD = SQRT(LB**2 - 1.0_EB)

         ! Head to back ratio based on LB
         HB = (LB + LBD) / (LB - LBD)

         ! A_ELPS and B_ELPS notation is consistent with Farsite and Richards 
         B_ELPS =  0.5_EB * (ROS_TMP + ROS_TMP/HB)
         B_ELPS2 = B_ELPS**2
         A_ELPS =  B_ELPS / LB
         A_ELPS2=  A_ELPS**2
         C_ELPS =  B_ELPS - (ROS_TMP/HB)

         ! Denominator used in spread rate equation from Richards, Intnl. J. Num. Methods Eng. 1990 
         ! and in LS vs Farsite paper, Bova et al., Intnl. J. Wildland Fire, 25(2):229-241, 2015  
         AROS  = XSF*COS_THETA - YSF*SIN_THETA
         BROS  = XSF*SIN_THETA + YSF*COS_THETA
         DENOM = A_ELPS2*BROS**2 + B_ELPS2*AROS**2

         IF (DENOM > 0._EB) THEN
            DENOM = 1._EB / SQRT(DENOM)
         ELSE
            DENOM = 0._EB
         ENDIF

!        This is with A_ELPS2 and B_ELPS2 notation consistent with Finney and Richards and in 
!        Bova et al. 2015 IJWF 2015
         SR_X_LS(I,J) = DENOM * ( A_ELPS2*COS_THETA*BROS - B_ELPS2*SIN_THETA*AROS) + C_ELPS*SIN_THETA
         SR_Y_LS(I,J) = DENOM * (-A_ELPS2*SIN_THETA*BROS - B_ELPS2*COS_THETA*AROS) + C_ELPS*COS_THETA

         ! Project spread rates from slope to horizontal plane

         IF (ABS(DZTDX(I,J)) > 0._EB) SR_X_LS(I,J) = SR_X_LS(I,J) * ABS(COS(ATAN(DZTDX(I,J))))
         IF (ABS(DZTDY(I,J)) > 0._EB) SR_Y_LS(I,J) = SR_Y_LS(I,J) * ABS(COS(ATAN(DZTDY(I,J))))

         MAG_SR = SQRT(SR_X_LS(I,J)**2 + SR_Y_LS(I,J)**2)

      ELSE ! McArthur Spread Model

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

         IF (DEGREES_SLOPE >= 5._EB  .AND. DEGREES_SLOPE < 10._EB) ROS_HEADS = 0.33_EB*ROS_HEAD(I,J)
         IF (DEGREES_SLOPE >= 10._EB .AND. DEGREES_SLOPE < 20._EB) ROS_HEADS =         ROS_HEAD(I,J)
         IF (DEGREES_SLOPE >= 20._EB)                              ROS_HEADS = 3._EB*  ROS_HEAD(I,J)

         MAG_SR    = 0.0_EB
         ROS_HEADS = 0.0_EB
         ROS_BACKS = 0.0_EB

         NEXP_WIND = WIND_EXP(I,J)

         ! Spread with the wind and upslope

         IF (COS_THETA_WIND >= 0._EB .AND. COS_THETA_SLOPE >= 0._EB) THEN
            IF (.NOT. LSET_TAN2) THEN
                IF (DEGREES_SLOPE >= 5._EB  .AND. DEGREES_SLOPE < 10._EB) ROS_HEADS = 0.33_EB*ROS_HEAD(I,J)
                IF (DEGREES_SLOPE >= 10._EB .AND. DEGREES_SLOPE < 20._EB) ROS_HEADS =         ROS_HEAD(I,J)
                IF (DEGREES_SLOPE >= 20._EB)                              ROS_HEADS =  3._EB*ROS_HEAD(I,J)
            ELSEIF (DEGREES_SLOPE > 0._EB) THEN
                ROS_HEADS = ROS_HEAD(I,J) * SLOPE_FACTOR !Dependence on TAN(slope)^2
            ENDIF
            MAG_SR = ROS_FLANK(I,J)*(1._EB + COS_THETA_WIND**NEXP_WIND*COS_THETA_SLOPE) + &
                     (ROS_HEAD(I,J) - ROS_FLANK(I,J))*COS_THETA_WIND**NEXP_WIND + &
                     (ROS_HEADS     - ROS_FLANK(I,J))*COS_THETA_SLOPE  !magnitude of spread rate
         ENDIF

         ! Spread with the wind and downslope

         IF (COS_THETA_WIND >= 0._EB .AND. COS_THETA_SLOPE < 0._EB) THEN
            IF (DEGREES_SLOPE >= 5._EB  .AND. DEGREES_SLOPE < 10._EB) ROS_HEADS =  0.33_EB*ROS_HEAD(I,J)
            IF (DEGREES_SLOPE >= 10._EB .AND. DEGREES_SLOPE < 20._EB) ROS_HEADS =  0.50_EB*ROS_HEAD(I,J)
            IF (DEGREES_SLOPE >= 20._EB)                              ROS_HEADS =  0.75_EB*ROS_HEAD(I,J)
            MAG_SR = ROS_FLANK(I,J)*(1._EB + COS_THETA_WIND*COS_THETA_SLOPE) + &
                     (ROS_HEAD(I,J) - ROS_FLANK(I,J))*COS_THETA_WIND**NEXP_WIND + &
                     (ROS_HEADS     - ROS_FLANK(I,J))*COS_THETA_SLOPE  !magnitude of spread rate
         ENDIF

         ! Spread against the wind and upslope

         IF (COS_THETA_WIND <  0._EB .AND. COS_THETA_SLOPE >= 0._EB) THEN
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

         IF (COS_THETA_WIND <  0._EB .AND. COS_THETA_SLOPE < 0._EB) THEN
            IF (DEGREES_SLOPE >= 5._EB .AND. DEGREES_SLOPE < 10._EB)  ROS_BACKS = 0.33_EB*ROS_BACKU(I,J)
            IF (DEGREES_SLOPE >= 10._EB .AND. DEGREES_SLOPE < 20._EB) ROS_BACKS = 0.50_EB*ROS_BACKU(I,J)
            IF (DEGREES_SLOPE >= 20._EB)                              ROS_BACKS = 0.75_EB*ROS_BACKU(I,J)
            MAG_SR = ROS_FLANK(I,J)*(1._EB - ABS(COS_THETA_WIND)**NEXP_WIND*COS_THETA_SLOPE) + &
                     (ROS_FLANK(I,J) - ROS_BACKU(I,J))*(-ABS(COS_THETA_WIND)**NEXP_WIND) + &
                     (ROS_FLANK(I,J) - ROS_BACKS)*COS_THETA_SLOPE  !magnitude of spread rate
         ENDIF

         SR_X_LS(I,J) = MAG_SR*NORMAL_FIRELINE(1) !spread rate components
         SR_Y_LS(I,J) = MAG_SR*NORMAL_FIRELINE(2)

      ENDIF !Ellipse or McArthur Spread

   ENDDO

ENDDO FLUX_ILOOP

END SUBROUTINE LEVEL_SET_SPREAD_RATE


!> \brief Compute the flux terms of the level set equation
!>
!> \details Use the spread rate [SR_X_LS,SR_Y_LS] to compute the limited scalar gradient and take dot product with 
!> spread rate vector to get advective flux

SUBROUTINE LEVEL_SET_ADVECT_FLUX

INTEGER :: I,IM1,IP1,IP2,J,JM1,JP1,JP2
REAL(EB), DIMENSION(:) :: Z(4)
REAL(EB), POINTER, DIMENSION(:,:) :: FLUX_LS_P,F_X,F_Y
REAL(EB) :: DPHIDX,DPHIDY,SR_X_AVG,SR_Y_AVG

F_X => LS_WORK1
F_Y => LS_WORK2

IF (PREDICTOR) THEN
   PHI_LS_P => PHI_LS
   FLUX_LS_P => FLUX0_LS
ELSE
   PHI_LS_P => PHI1_LS
   FLUX_LS_P => FLUX1_LS
ENDIF

DO J=1,JBAR
   DO I=0,IBAR
      IM1 = I-1 ; IF (IM1<0) IM1 = I
      IP1 = I+1
      IP2 = I+2 ; IF (IP2>IBP1) IP2 = IP1
      Z(1) = PHI_LS_P(IM1,J)
      Z(2) = PHI_LS_P(I,J)
      Z(3) = PHI_LS_P(IP1,J)
      Z(4) = PHI_LS_P(IP2,J)
      SR_X_AVG = 0.5_EB*(SR_X_LS(MIN(IP1,IBAR),J)+SR_X_LS(MAX(1,I),J))
      F_X(I,J) = SCALAR_FACE_VALUE_LS(SR_X_AVG,Z,LIMITER_LS)
   ENDDO
ENDDO

DO J=0,JBAR
   DO I=1,IBAR
      JM1 = J-1 ; IF (JM1<0) JM1 = J
      JP1 = J+1
      JP2 = J+2 ; IF (JP2>JBP1) JP2 = JP1
      Z(1) = PHI_LS_P(I,JM1)
      Z(2) = PHI_LS_P(I,J)
      Z(3) = PHI_LS_P(I,JP1)
      Z(4) = PHI_LS_P(I,JP2)
      SR_Y_AVG = 0.5_EB*(SR_Y_LS(I,MIN(JP1,JBAR))+SR_Y_LS(I,MAX(1,J)))
      F_Y(I,J) = SCALAR_FACE_VALUE_LS(SR_Y_AVG,Z,LIMITER_LS)
   ENDDO
ENDDO

DO J=1,JBAR
   DO I=1,IBAR
      DPHIDX = (F_X(I,J)-F_X(I-1,J))*RDX(I)
      DPHIDY = (F_Y(I,J)-F_Y(I,J-1))*RDY(J)
      FLUX_LS_P(I,J) = SR_X_LS(I,J)*DPHIDX + SR_Y_LS(I,J)*DPHIDY
   ENDDO
ENDDO

END SUBROUTINE LEVEL_SET_ADVECT_FLUX


!> \brief Compute the scalar value on the flux at a cell face
!>
!> \param SR_XY If positive, indicates that flow is from left to right
!> \param Z Scalar quantity
!> \param LIMITER Flux limiter (1) MINMOD, (2) SUPERBEE, (3) first-order upwinding (monotone)
!> \details
!                           face
!    |     o     |     o     |     o     |     o     |
!         Z(1)        Z(2)        Z(3)        Z(4)
!
REAL(EB) FUNCTION SCALAR_FACE_VALUE_LS(SR_XY,Z,LIMITER)

INTEGER, INTENT(IN) :: LIMITER
REAL(EB) :: SR_XY
REAL(EB), INTENT(IN), DIMENSION(4) :: Z
REAL(EB) :: B,DZLOC,DZUP,R,ZUP,ZDWN

IF (SR_XY > 0._EB) THEN  ! flow is left to right

   DZLOC = Z(3)-Z(2)
   DZUP  = Z(2)-Z(1)
   IF (ABS(DZLOC) > 0._EB) THEN
      R = DZUP/DZLOC
   ELSE
      R = 0._EB
   ENDIF
   ZUP  = Z(2)
   ZDWN = Z(3)

ELSE  ! flow is right to left

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

IF (LIMITER==1) THEN
   B = MAX(0._EB,MIN(1._EB,R))
ELSEIF (LIMITER==2) THEN
   B = MAX(0._EB,MIN(2._EB*R,1._EB),MIN(R,2._EB))
ELSEIF (LIMITER==3) THEN
   B = 0._EB
ENDIF

SCALAR_FACE_VALUE_LS = ZUP + 0.5_EB * B * ( ZDWN - ZUP )

END FUNCTION SCALAR_FACE_VALUE_LS


!> \brief Calculate the Rothermel no-wind, no-slope rate of spread.
!> 
!>
!> \details The Rothermel model as described in Bachmann's thesis.

REAL(EB) FUNCTION ROS_NO_WIND_NO_SLOPE(ROTHERMEL_FUEL_INDEX,SURF_INDEX)

INTEGER, INTENT(IN) :: ROTHERMEL_FUEL_INDEX,SURF_INDEX
REAL(EB) :: w0d1, w0d2, w0d3, w0lh, w0lw, md1, md2, md3, mlh, mlw, svd1, svd2, svd3, svlh, svlw, depth, rhop, heat, st, se, mx
REAL(EB) :: swd1, swd2, swd3, swlh, swlw, swd, swl, swt, s2wt, sw2d, sw2l, swmd, swml, sigma, rhob, beta, &
            betaOpt, wnd, wnl, hnd1, hnd2, hnd3, hnlh, hnlw, hnd, hnl, bigW, hnmd, mfdead, mxlive, rml, rmd, etaMd, etaMl, etaM, &
            etas, gammaMax, bigA, gamma, bigIr, xi, epsd1, epsd2, epsd3, epslh, epslw, bigQd1, bigQd2, &
            bigQd3, bigQlh, bigQlw, hskz, hsk
TYPE(SURFACE_TYPE), POINTER :: SF

SF => SURFACE(SURF_INDEX)

md1 = SF%VEG_LSET_M1
md2 = SF%VEG_LSET_M10
md3 = SF%VEG_LSET_M100
mlw = SF%VEG_LSET_MLW
mlh = SF%VEG_LSET_MLH
               
SELECT CASE(ROTHERMEL_FUEL_INDEX)
   CASE(1)  ! 'Short Grass'
      w0d1=0.1659     ; w0d2=0.        ; w0d3=0.        ; w0lh=0.        ; w0lw=0.     ! dry mass per unit area (kg/m2)
      svd1=11483.     ; svd2=358.      ; svd3=98.       ; svlh=4921.     ; svlw=4921.  ! surface area to volume (1/m)
      mx=0.12         ; depth=0.3048   ; rhop=512.      ; heat=18607.    ; st=0.0555      ; se=0.01
   CASE(2)  ! 'Timbergrass'
      w0d1=0.448      ; w0d2=0.224     ; w0d3=0.112     ; w0lh=0.112     ; w0lw=0. 
      svd1=9842.      ; svd2=358.      ; svd3=98.       ; svlh=4921.     ; svlw=4921. 
      mx=0.15         ; depth=0.3048   ; rhop=512.      ; heat=18607.    ; st=0.0555      ; se=0.01
   CASE(3)  ! 'Tall Grass'
      w0d1=0.675      ; w0d2=0.        ; w0d3=0.        ; w0lh=0.        ; w0lw=0. 
      svd1=4921.      ; svd2=358.      ; svd3=98.       ; svlh=4921.     ; svlw=4921. 
      mx=0.25         ; depth=0.762    ; rhop=512.      ; heat=18607.    ; st=0.0555      ; se=0.01
   CASE(4)  ! 'Chaparral'
      w0d1=1.123      ; w0d2=0.899     ; w0d3=0.448     ; w0lh=1.123     ; w0lw=0. 
      svd1=6562.      ; svd2=358.      ; svd3=98.       ; svlh=4921.     ; svlw=4921. 
      mx=0.20         ; depth=1.829    ; rhop=512.      ; heat=18607.    ; st=0.0555      ; se=0.01
   CASE(5)  ! 'Brush'
      w0d1=0.224      ; w0d2=0.112     ; w0d3=0.        ; w0lh=0.        ; w0lw=0.448
      svd1=6562.      ; svd2=358.      ; svd3=98.       ; svlh=4921.     ; svlw=4921. 
      mx=0.20         ; depth=0.6096   ; rhop=512.      ; heat=18607.    ; st=0.0555      ; se=0.01
   CASE(6)  ! 'Dormant Brush'
      w0d1=0.336      ; w0d2=0.56      ; w0d3=0.448     ; w0lh=0.        ; w0lw=0.   
      svd1=5741.      ; svd2=358.      ; svd3=98.       ; svlh=4921.     ; svlw=4921. 
      mx=0.25         ; depth=0.762    ; rhop=512.      ; heat=18607.    ; st=0.0555      ; se=0.01
   CASE(7)  ! 'Southern Rough'
      w0d1=0.255      ; w0d2=0.419     ; w0d3=0.336     ; w0lh=0.        ; w0lw=0.083
      svd1=5741.      ; svd2=358.      ; svd3=98.       ; svlh=4921.     ; svlw=4921. 
      mx=0.40         ; depth=0.762    ; rhop=512.      ; heat=18607.    ; st=0.0555      ; se=0.01
   CASE(8)  ! 'Closed Timber Litter'
      w0d1=0.336      ; w0d2=0.224     ; w0d3=0.56      ; w0lh=0.        ; w0lw=0.   
      svd1=6562.      ; svd2=358.      ; svd3=98.       ; svlh=4921.     ; svlw=4921. 
      mx=0.30         ; depth=0.06096  ; rhop=512.      ; heat=18607.    ; st=0.0555      ; se=0.01
   CASE(9)  ! ID='Hardwood Litter'
      w0d1=0.655      ; w0d2=0.092     ; w0d3=0.034     ; w0lh=0.        ; w0lw=0.   
      svd1=8202.      ; svd2=358.      ; svd3=98.       ; svlh=4921.     ; svlw=4921. 
      mx=0.25         ; depth=0.06096  ; rhop=512.      ; heat=18607.    ; st=0.0555      ; se=0.01
   CASE(10)  ! 'Timber'
      w0d1=0.675      ; w0d2=0.448     ; w0d3=1.123     ; w0lh=0.        ; w0lw=0.448
      svd1=6562.      ; svd2=358.      ; svd3=98.       ; svlh=4921.     ; svlw=4921. 
      mx=0.25         ; depth=0.3048   ; rhop=512.      ; heat=18607.    ; st=0.0555      ; se=0.01
   CASE(11)  ! 'Light Slash'
      w0d1=0.336      ; w0d2=1.011     ; w0d3=1.235     ; w0lh=0.        ; w0lw=0.   
      svd1=4921.      ; svd2=358.      ; svd3=98.       ; svlh=4921.     ; svlw=4921. 
      mx=0.15         ; depth=0.3048   ; rhop=512.      ; heat=18607.    ; st=0.0555      ; se=0.01
   CASE(12)  ! ID='Medium Slash'
      w0d1=0.899      ; w0d2=3.145     ; w0d3=3.706     ; w0lh=0.        ; w0lw=0.   
      svd1=4921.      ; svd2=358.      ; svd3=98.       ; svlh=4921.     ; svlw=4921. 
      mx=0.20         ; depth=0.70104  ; rhop=512.      ; heat=18607.    ; st=0.0555      ; se=0.01
   CASE(13)  ! 'Heavy Slash'
      w0d1=1.571      ; w0d2=5.165     ; w0d3=6.288     ; w0lh=0.        ; w0lw=0.   
      svd1=4921.      ; svd2=358.      ; svd3=98.       ; svlh=4921.     ; svlw=4921. 
      mx=0.25         ; depth=0.9144   ; rhop=512.      ; heat=18607.    ; st=0.0555      ; se=0.01
END SELECT

SF%VEG_LSET_HT = depth
      
! Auxiliary functions
   
swd1 = svd1*w0d1  
swd2 = svd2*w0d2 
swd3 = svd3*w0d3  
swlh = svlh*w0lh  
swlw = svlw*w0lw 
swd  = swd1 + swd2 + swd3  
swl  = swlh + swlw  
swt  = swd + swl  
s2wt = svd1**2*w0d1 + svd2**2*w0d2 + svd3**2*w0d3 + svlh**2*w0lh + svlw**2*w0lw  
sw2d = svd1*w0d1**2 + svd2*w0d2**2 + svd3*w0d3**2
sw2l = svlh*w0lh**2 + svlw*w0lw**2  
swmd = swd1*md1 + swd2*md2 + swd3*md3  
swml = swlh*mlh + swlw*mlw
   
! Characteristic surface-to-volume ratio [R(71,72)]

sigma = s2wt/swt
SF%VEG_LSET_SIGMA = sigma*0.01_EB  ! Convert from 1/m to 1/cm
   
! Mean bulk density [R(74)]
   
rhob = (w0d1 + w0d2 + w0d3 + w0lh + w0lw)/depth
  
! Mean packing ratio [R(31,73)]
   
beta = rhob/rhop
SF%VEG_LSET_BETA = beta
   
! Optimal packing ratio [R(37)]
   
betaOpt = 8.8578*sigma**(-0.8189)
   
! Net fuel loading [R(60), adjusted by A.(p.88) and R(59)]
   
if (swd==0._eb) then
   wnd = 0._eb
else
   wnd = (sw2d/swd)*(1. - st)
endif

if (swl==0._eb) then
   wnl = 0._eb
else
   wnl = (sw2l/swl)*(1. - st)
endif
   
! Mineral damping coefficient [R(62)]
   
etas = 0.174*se**(-0.19)
   
! Ratio of "fine" fuel loadings,dead/living [Albini,p.89]
   
hnd1 = 0.20482*w0d1*Exp(-452.76/svd1) 
hnd2 = 0.20482*w0d2*exp(-452.76/svd2) 
hnd3 = 0.20482*w0d3*exp(-452.76/svd3) 
hnlh = 0.20482*w0lh*exp(-1640.42/svlh) 
hnlw = 0.20482*w0lw*exp(-1640.42/svlw) 
hnd = hnd1 + hnd2 + hnd3 
hnl = hnlh + hnlw
if (swl==0._eb) then
   bigW = 0._eb
else
   bigW = hnd/hnl
endif
   
! Moisture content of "fine" dead fuel [Albini,p.89]
   
hnmd   = hnd1*md1 + hnd2*md2 + hnd3*md3
mfdead = hnmd/hnd
   
! Moisture of extinction of living fuel [R(88),Albini,p.89]
   
mxlive = 2.9*bigW*(1.0 - (mfdead/mx)) - 0.226
   
! Moisture ratios [R(65,66)]
   
if (swl==0._eb) then
   rml = 0._eb
else
   rml = swml/(swl*mxlive)
endif

rmd = swmd/(swd*mx)
   
! Moisture damping coefficients [R(64)]
   
etaMd = 1.0 - (2.59*rmd) + (5.11*rmd**2) - (3.52*rmd**3) 
etaMl = 1.0 - (2.59*rml) + (5.11*rml**2) - (3.52*rml**3) 
etaM  = wnd*etaMd + wnl*etaMl
   
! Maximum reaction velocity [R(36,68)]
   
gammaMax = (0.16828*sigma**(1.5))/(29700 + 0.5997*sigma**(1.5))
   
! A [R(70),Albini p.88]
   
bigA = 340.53*sigma**(-0.7913)
   
! Potential reaction velocity [R(38)]
   
gamma = gammaMax*(beta/betaOpt)**(bigA)*exp(bigA*(1.0 - (beta/betaOpt)))
   
! Propagating flux ratio [R(42)]
   
xi = exp((0.792 + 0.37597*sqrt(sigma))*(beta + 0.1))/(192.0 + 0.0791*sigma)
   
! Effective heating number [R(14,77)]
   
epsd1 = exp(-452.76/svd1) 
epsd2 = exp(-452.76/svd2) 
epsd3 = exp(-452.76/svd3)
epslh = exp(-452.76/svlh) 
epslw = exp(-452.76/svlw)
   
! Heat of pre-ignition [R(12,78)]
   
bigQd1 = 581.5 + 2595.7*md1 
bigQd2 = 581.5 + 2595.7*md2 
bigQd3 = 581.5 + 2595.7*md3 
bigQlh = 581.5 + 2595.7*mlh 
bigQlw = 581.5 + 2595.7*mlw
   
! Heat sink [R(77)]
   
hskz = svd1*w0d1*epsd1*bigQd1 + svd2*w0d2*epsd2*bigQd2 + svd3*w0d3*epsd3*bigQd3 + svlh*w0lh*epslh*bigQlh + svlw*w0lw*epslw*bigQlw
hsk  = rhob*hskz/swt
   
! Reaction intensity [R(27,58),Albini,p.89]
   
bigIr = gamma*heat*etas*etaM

IF (LEVEL_SET_COUPLED_FIRE) THEN
   SF%MASS_FLUX(REACTION(1)%FUEL_SMIX_INDEX) = bigIr/heat
   SF%BURN_DURATION = 756._EB/SF%VEG_LSET_SIGMA   ! Albini (Eq. 14)
ENDIF
   
! Rate of spread [R(52)] and the rate of spread in the absence of wind and with no slope.
   
ROS_NO_WIND_NO_SLOPE = (bigIr*xi)/hsk

END FUNCTION ROS_NO_WIND_NO_SLOPE

END MODULE VEGE
