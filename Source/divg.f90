MODULE DIVG

USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_POINTERS

IMPLICIT NONE
PRIVATE

PUBLIC DIVERGENCE_PART_1,DIVERGENCE_PART_2,CHECK_DIVERGENCE

CONTAINS

SUBROUTINE DIVERGENCE_PART_1(T,DT,NM)

USE COMP_FUNCTIONS, ONLY: CURRENT_TIME
USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP,INTERPOLATE1D_UNIFORM
USE PHYSICAL_FUNCTIONS, ONLY: GET_CONDUCTIVITY,GET_SPECIFIC_HEAT,GET_SENSIBLE_ENTHALPY_Z,GET_SENSIBLE_ENTHALPY,&
                              GET_VISCOSITY,GET_MOLECULAR_WEIGHT
USE TURBULENCE, ONLY: TENSOR_DIFFUSIVITY_MODEL,WANNIER_FLOW
USE MASS, ONLY: SCALAR_FACE_VALUE
USE GEOMETRY_FUNCTIONS, ONLY: ASSIGN_PRESSURE_ZONE
USE MANUFACTURED_SOLUTIONS, ONLY: DIFF_MMS,UF_MMS,WF_MMS,VD2D_MMS_Z_SRC !,RHO_0_MMS,RHO_1_MMS
USE EVAC, ONLY: EVAC_EMESH_EXITS_TYPE, EMESH_EXITS, EMESH_NFIELDS, N_EXITS, N_CO_EXITS, N_DOORS
USE COMPLEX_GEOMETRY, ONLY : SET_EXIMDIFFLX_3D,SET_DOMAINDIFFLX_3D,SET_EXIMRHOHSLIM_3D,SET_EXIMRHOZZLIM_3D,&
                             IBM_CGSC, IBM_IDCC, IBM_SOLID, IBM_CUTCFE, &
                             CCREGION_DIVERGENCE_PART_1,CFACE_PREDICT_NORMAL_VELOCITY

! Compute contributions to the divergence term

INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: T,DT
REAL(EB), POINTER, DIMENSION(:,:,:) :: KDTDX,KDTDY,KDTDZ,DP,KP,CP, &
          RHO_D,RHOP,H_RHO_D_DZDX,H_RHO_D_DZDY,H_RHO_D_DZDZ,RTRM, &
          U_DOT_DEL_RHO_H_S,RHO_H_S_P,UU,VV,WW,RHO_Z_P,U_DOT_DEL_RHO_Z,RHO_D_TURB,R_H_G
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZP,RHO_D_DZDX,RHO_D_DZDY,RHO_D_DZDZ
REAL(EB), POINTER, DIMENSION(:,:) :: PBAR_P
REAL(EB) :: DELKDELT,VC,VC1,DTDX,DTDY,DTDZ,TNOW, &
            DZDX,DZDY,DZDZ,RDT,RHO_D_DZDN,TSI,TIME_RAMP_FACTOR,DELTA_P,PRES_RAMP_FACTOR,&
            TMP_G,DIV_DIFF_HEAT_FLUX,H_S,ZZZ(1:4),DU,DU_P,DU_M,UN,PROFILE_FACTOR, &
            XHAT,ZHAT,TT,Q_Z,D_Z_TEMP,D_Z_N(0:5000),RHO_D_DZDN_GET(1:N_TRACKED_SPECIES),JCOR,UN_P,TMP_F_GAS,R_PFCT
REAL(EB), ALLOCATABLE, DIMENSION(:) :: ZZ_GET
TYPE(SURFACE_TYPE), POINTER :: SF
TYPE(SPECIES_MIXTURE_TYPE), POINTER :: SM
INTEGER :: IW,N,IOR,II,JJ,KK,IIG,JJG,KKG,I,J,K,IPZ,IOPZ,N_ZZ_MAX,ICC
TYPE(VENTS_TYPE), POINTER :: VT=>NULL()
TYPE(WALL_TYPE), POINTER :: WC=>NULL()

! Check whether to skip this routine

IF (SOLID_PHASE_ONLY) RETURN

! Start the clock and set the pointers

TNOW=CURRENT_TIME()
CALL POINT_TO_MESH(NM)

RDT = 1._EB/DT

SELECT CASE(PREDICTOR)
   CASE(.TRUE.)
      DP => DS
      PBAR_P => PBAR_S
      RHOP => RHOS
   CASE(.FALSE.)
      DP => DDDT
      PBAR_P => PBAR
      RHOP => RHO
END SELECT

R_PBAR = 1._EB/PBAR_P

RTRM => WORK1

! Zero out divergence to start

DP = 0._EB

! Determine if pressure ZONEs have merged

CALL MERGE_PRESSURE_ZONES

! Compute normal component of velocity at boundaries, U_NORMAL_S in the PREDICTOR step, U_NORMAL in the CORRECTOR.

CALL PREDICT_NORMAL_VELOCITY
IF (CC_IBM) CALL CFACE_PREDICT_NORMAL_VELOCITY(T,DT)

IF (EVACUATION_ONLY(NM)) GOTO 1000 ! Evacuation meshes jump directly to pressure zone loop

! Compute species-related finite difference terms

RHO_D_DZDX => SCALAR_WORK1
RHO_D_DZDY => SCALAR_WORK2
RHO_D_DZDZ => SCALAR_WORK3
SELECT CASE(PREDICTOR)
   CASE(.TRUE.)
      ZZP => ZZS
   CASE(.FALSE.)
      ZZP => ZZ
END SELECT

! Save the largest value of the material and thermal diffusion coefficients for use in Von Neumann stability constraint

IF (CHECK_VN) D_Z_MAX = 0._EB

! Add species diffusion terms to divergence expression and compute diffusion term for species equations

SPECIES_GT_1_IF: IF (N_TOTAL_SCALARS>1) THEN

   DEL_RHO_D_DEL_Z = 0._EB
   RHO_D => WORK4
   IF (.NOT.POTENTIAL_TEMPERATURE_CORRECTION) THEN
      ! default
      IF (SIM_MODE/=DNS_MODE) THEN
         IF (SIM_MODE==LES_MODE) THEN
            RHO_D_TURB => WORK9
            RHO_D_TURB = MAX(0._EB,MU-MU_DNS)*RSC
         ELSE
            RHO_D = MAX(0._EB,MU)*RSC
         ENDIF
      ENDIF
   ELSE
      ! dynamic turbulent Schmidt number (Deardorff, 1980)
      IF (SIM_MODE/=DNS_MODE) THEN
         IF (SIM_MODE==LES_MODE) THEN
            RHO_D_TURB => WORK9
            RHO_D_TURB = MAX(0._EB,MU-MU_DNS)/PR_T
         ELSE
            RHO_D = MAX(0._EB,MU)/PR_T
         ENDIF
      ENDIF
   ENDIF

   DIFFUSIVE_FLUX_LOOP: DO N=1,N_TOTAL_SCALARS

      IF (SIM_MODE==DNS_MODE .OR. SIM_MODE==LES_MODE) THEN
         RHO_D = 0._EB
         D_Z_N = D_Z(:,N)
         DO K=0,KBP1
            DO J=0,JBP1
               DO I=0,IBP1
                  CALL INTERPOLATE1D_UNIFORM(LBOUND(D_Z_N,1),D_Z_N,TMP(I,J,K),D_Z_TEMP)
                  RHO_D(I,J,K) = RHOP(I,J,K)*D_Z_TEMP
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      IF (SIM_MODE==LES_MODE) RHO_D = RHO_D + RHO_D_TURB

      ! Manufactured solution

      IF (PERIODIC_TEST==7) RHO_D = DIFF_MMS

      ! Store max diffusivity for stability check

      IF (CHECK_VN) D_Z_MAX = MAX(D_Z_MAX,RHO_D/RHOP)

      ! Compute rho*D del Z

      !$OMP PARALLEL DO PRIVATE(DZDX, DZDY, DZDZ) SCHEDULE (STATIC)
      DO K=0,KBAR
         DO J=0,JBAR
            DO I=0,IBAR
               DZDX = (ZZP(I+1,J,K,N)-ZZP(I,J,K,N))*RDXN(I)
               RHO_D_DZDX(I,J,K,N) = .5_EB*(RHO_D(I+1,J,K)+RHO_D(I,J,K))*DZDX
               DZDY = (ZZP(I,J+1,K,N)-ZZP(I,J,K,N))*RDYN(J)
               RHO_D_DZDY(I,J,K,N) = .5_EB*(RHO_D(I,J+1,K)+RHO_D(I,J,K))*DZDY
               DZDZ = (ZZP(I,J,K+1,N)-ZZP(I,J,K,N))*RDZN(K)
               RHO_D_DZDZ(I,J,K,N) = .5_EB*(RHO_D(I,J,K+1)+RHO_D(I,J,K))*DZDZ
            ENDDO
         ENDDO
      ENDDO
      !$OMP END PARALLEL DO

      ! Tensor diffusivity model (experimental)

      IF (TENSOR_DIFFUSIVITY .AND. SIM_MODE/=DNS_MODE) CALL TENSOR_DIFFUSIVITY_MODEL(NM,N)

      ! Store rho*D_n grad Z_n at OPEN boundaries, zero out otherwise

      WALL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
         WC => WALL(IW)
         IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY .OR. &
             WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY) CYCLE WALL_LOOP
         IIG = WC%ONE_D%IIG
         JJG = WC%ONE_D%JJG
         KKG = WC%ONE_D%KKG
         IOR = WC%ONE_D%IOR
         IF (WC%BOUNDARY_TYPE==OPEN_BOUNDARY) THEN
            SELECT CASE(IOR)
               CASE( 1); WC%ONE_D%RHO_D_DZDN_F(N) =  RHO_D_DZDX(IIG-1,JJG,KKG,N)
               CASE(-1); WC%ONE_D%RHO_D_DZDN_F(N) = -RHO_D_DZDX(IIG,JJG,KKG,N)
               CASE( 2); WC%ONE_D%RHO_D_DZDN_F(N) =  RHO_D_DZDY(IIG,JJG-1,KKG,N)
               CASE(-2); WC%ONE_D%RHO_D_DZDN_F(N) = -RHO_D_DZDY(IIG,JJG,KKG,N)
               CASE( 3); WC%ONE_D%RHO_D_DZDN_F(N) =  RHO_D_DZDZ(IIG,JJG,KKG-1,N)
               CASE(-3); WC%ONE_D%RHO_D_DZDN_F(N) = -RHO_D_DZDZ(IIG,JJG,KKG,N)
            END SELECT
         ELSE
            SELECT CASE(IOR)
               CASE( 1); RHO_D_DZDX(IIG-1,JJG,KKG,N) = 0._EB
               CASE(-1); RHO_D_DZDX(IIG,JJG,KKG,N)   = 0._EB
               CASE( 2); RHO_D_DZDY(IIG,JJG-1,KKG,N) = 0._EB
               CASE(-2); RHO_D_DZDY(IIG,JJG,KKG,N)   = 0._EB
               CASE( 3); RHO_D_DZDZ(IIG,JJG,KKG-1,N) = 0._EB
               CASE(-3); RHO_D_DZDZ(IIG,JJG,KKG,N)   = 0._EB
            END SELECT
         ENDIF
      ENDDO WALL_LOOP

      IF (CHECK_MASS_CONSERVE) THEN
         ! When CHECK_MASS_CONSERVE make zero diffusive mass fluxes on open boundaries.
         ! Using the open boundary mass fluxes leads to a small difference on the integrals done in this test.
         DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
            WC => WALL(IW)
            IF (WC%BOUNDARY_TYPE==OPEN_BOUNDARY) THEN
               IIG = WC%ONE_D%IIG
               JJG = WC%ONE_D%JJG
               KKG = WC%ONE_D%KKG
               IOR = WC%ONE_D%IOR
               SELECT CASE(IOR)
                  CASE( 1)
                     RHO_D_DZDX(IIG-1,JJG,KKG,N) = 0._EB
                  CASE(-1)
                     RHO_D_DZDX(IIG,JJG,KKG,N)   = 0._EB
                  CASE( 2)
                     RHO_D_DZDY(IIG,JJG-1,KKG,N) = 0._EB
                  CASE(-2)
                     RHO_D_DZDY(IIG,JJG,KKG,N)   = 0._EB
                  CASE( 3)
                     RHO_D_DZDZ(IIG,JJG,KKG-1,N) = 0._EB
                  CASE(-3)
                     RHO_D_DZDZ(IIG,JJG,KKG,N)   = 0._EB
               END SELECT
            ENDIF
         ENDDO
      ENDIF

   ENDDO DIFFUSIVE_FLUX_LOOP

   ! Ensure RHO_D terms sum to zero over all species.  Gather error into largest mass fraction present.

   IF (SIM_MODE==DNS_MODE .OR. SIM_MODE==LES_MODE) THEN

      !$OMP PARALLEL DO PRIVATE(N) SCHEDULE(STATIC)
      DO K=0,KBAR
         DO J=0,JBAR
            DO I=0,IBAR
               N=MAXLOC(ZZP(I,J,K,1:N_TRACKED_SPECIES)+ZZP(I+1,J,K,1:N_TRACKED_SPECIES),1)
               RHO_D_DZDX(I,J,K,N) = -(SUM(RHO_D_DZDX(I,J,K,1:N_TRACKED_SPECIES))-RHO_D_DZDX(I,J,K,N))

               N=MAXLOC(ZZP(I,J,K,1:N_TRACKED_SPECIES)+ZZP(I,J+1,K,1:N_TRACKED_SPECIES),1)
               RHO_D_DZDY(I,J,K,N) = -(SUM(RHO_D_DZDY(I,J,K,1:N_TRACKED_SPECIES))-RHO_D_DZDY(I,J,K,N))

               N=MAXLOC(ZZP(I,J,K,1:N_TRACKED_SPECIES)+ZZP(I,J,K+1,1:N_TRACKED_SPECIES),1)
               RHO_D_DZDZ(I,J,K,N) = -(SUM(RHO_D_DZDZ(I,J,K,1:N_TRACKED_SPECIES))-RHO_D_DZDZ(I,J,K,N))
            ENDDO
         ENDDO
      ENDDO
      !$OMP END PARALLEL DO

   ENDIF

   ! Store diffusive species flux on EXIM boundary faces if present

   IF (CC_IBM) CALL SET_EXIMDIFFLX_3D(NM,RHO_D_DZDX,RHO_D_DZDY,RHO_D_DZDZ)

   ! Store diffusive flux for output

   IF (STORE_SPECIES_FLUX) THEN
      IF (PREDICTOR) THEN
         DO N=1,N_TOTAL_SCALARS
            DIF_FX(:,:,:,N) = 0.5_EB*( DIF_FXS(:,:,:,N) - RHO_D_DZDX(:,:,:,N) )
            DIF_FY(:,:,:,N) = 0.5_EB*( DIF_FYS(:,:,:,N) - RHO_D_DZDY(:,:,:,N) )
            DIF_FZ(:,:,:,N) = 0.5_EB*( DIF_FZS(:,:,:,N) - RHO_D_DZDZ(:,:,:,N) )
         ENDDO
      ELSE
         DO N=1,N_TOTAL_SCALARS
            DIF_FXS(:,:,:,N) = -RHO_D_DZDX(:,:,:,N)
            DIF_FYS(:,:,:,N) = -RHO_D_DZDY(:,:,:,N)
            DIF_FZS(:,:,:,N) = -RHO_D_DZDZ(:,:,:,N)
         ENDDO
      ENDIF
   ENDIF

   ! Diffusive heat flux

   SPECIES_LOOP: DO N=1,N_TOTAL_SCALARS

      ! Compute div h_n*rho*D del Z_n (part of div qdot")

      H_RHO_D_DZDX => WORK5
      H_RHO_D_DZDY => WORK6
      H_RHO_D_DZDZ => WORK7

      !$OMP PARALLEL DO PRIVATE(TMP_G, H_S) SCHEDULE(guided)
      DO K=0,KBAR
         DO J=0,JBAR
            DO I=0,IBAR
               ! H_RHO_D_DZDX
               TMP_G = 0.5_EB*(TMP(I+1,J,K)+TMP(I,J,K))
               CALL GET_SENSIBLE_ENTHALPY_Z(N,TMP_G,H_S)
               H_RHO_D_DZDX(I,J,K) = H_S*RHO_D_DZDX(I,J,K,N)

               ! H_RHO_D_DZDY
               TMP_G = 0.5_EB*(TMP(I,J+1,K)+TMP(I,J,K))
               CALL GET_SENSIBLE_ENTHALPY_Z(N,TMP_G,H_S)
               H_RHO_D_DZDY(I,J,K) = H_S*RHO_D_DZDY(I,J,K,N)

               ! H_RHO_D_DZDZ
               TMP_G = 0.5_EB*(TMP(I,J,K+1)+TMP(I,J,K))
               CALL GET_SENSIBLE_ENTHALPY_Z(N,TMP_G,H_S)
               H_RHO_D_DZDZ(I,J,K) = H_S*RHO_D_DZDZ(I,J,K,N)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(DIV_DIFF_HEAT_FLUX) SCHEDULE(STATIC)
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR

               DIV_DIFF_HEAT_FLUX = (R(I)*H_RHO_D_DZDX(I,J,K)-R(I-1)*H_RHO_D_DZDX(I-1,J,K))*RDX(I)*RRN(I) + &
                                    (     H_RHO_D_DZDY(I,J,K)-       H_RHO_D_DZDY(I,J-1,K))*RDY(J)        + &
                                    (     H_RHO_D_DZDZ(I,J,K)-       H_RHO_D_DZDZ(I,J,K-1))*RDZ(K)

               DP(I,J,K) = DP(I,J,K) + DIV_DIFF_HEAT_FLUX

            ENDDO
         ENDDO
      ENDDO
      !$OMP END PARALLEL DO

      ! Compute div rho*D grad Z_n

      !$OMP PARALLEL DO SCHEDULE(STATIC)
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               DEL_RHO_D_DEL_Z(I,J,K,N) = (R(I)*RHO_D_DZDX(I,J,K,N)-R(I-1)*RHO_D_DZDX(I-1,J,K,N))*RDX(I)*RRN(I) + &
                                          (     RHO_D_DZDY(I,J,K,N)-       RHO_D_DZDY(I,J-1,K,N))*RDY(J)        + &
                                          (     RHO_D_DZDZ(I,J,K,N)-       RHO_D_DZDZ(I,J,K-1,N))*RDZ(K)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END PARALLEL DO

      ! Correct rho*D_n grad Z_n and h_n*rho*D_n grad Z_n at boundaries

      WALL_LOOP_2: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
         WC => WALL(IW)
         IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY .OR. &
             WC%BOUNDARY_TYPE==OPEN_BOUNDARY .OR. &
             WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY) CYCLE WALL_LOOP_2
         IIG = WC%ONE_D%IIG
         JJG = WC%ONE_D%JJG
         KKG = WC%ONE_D%KKG
         IOR = WC%ONE_D%IOR

         N_ZZ_MAX = MAXLOC(WC%ONE_D%ZZ_F(1:N_TRACKED_SPECIES),1)
         RHO_D_DZDN = 2._EB*WC%ONE_D%RHO_D_F(N)*(ZZP(IIG,JJG,KKG,N)-WC%ONE_D%ZZ_F(N))*WC%ONE_D%RDN
         IF (N==N_ZZ_MAX) THEN
            RHO_D_DZDN_GET = 2._EB*WC%ONE_D%RHO_D_F(:)*(ZZP(IIG,JJG,KKG,:)-WC%ONE_D%ZZ_F(:))*WC%ONE_D%RDN
            RHO_D_DZDN = -(SUM(RHO_D_DZDN_GET(:))-RHO_D_DZDN)
         ENDIF
         WC%ONE_D%RHO_D_DZDN_F(N) = RHO_D_DZDN

         IF (PREDICTOR) THEN
            UN_P = WC%ONE_D%U_NORMAL_S
         ELSE
            UN_P = WC%ONE_D%U_NORMAL
         ENDIF
         IF (WC%BOUNDARY_TYPE==SOLID_BOUNDARY .AND. UN_P>0._EB) THEN
            TMP_F_GAS = WC%ONE_D%TMP_G
         ELSE
            TMP_F_GAS = WC%ONE_D%TMP_F
         ENDIF

         CALL GET_SENSIBLE_ENTHALPY_Z(N,TMP_F_GAS,H_S)

         SELECT CASE(IOR)
            CASE( 1)
               JCOR = RHO_D_DZDN*RDX(IIG)*RRN(IIG)*R(IIG-1)
            CASE(-1)
               JCOR = RHO_D_DZDN*RDX(IIG)*RRN(IIG)*R(IIG)
            CASE( 2)
               JCOR = RHO_D_DZDN*RDY(JJG)
            CASE(-2)
               JCOR = RHO_D_DZDN*RDY(JJG)
            CASE( 3)
               JCOR = RHO_D_DZDN*RDZ(KKG)
            CASE(-3)
               JCOR = RHO_D_DZDN*RDZ(KKG)
         END SELECT

         DEL_RHO_D_DEL_Z(IIG,JJG,KKG,N) = DEL_RHO_D_DEL_Z(IIG,JJG,KKG,N) - JCOR
         DP(IIG,JJG,KKG) = DP(IIG,JJG,KKG) - H_S*JCOR

      ENDDO WALL_LOOP_2

   ENDDO SPECIES_LOOP

ENDIF SPECIES_GT_1_IF

! Store diffusive species flux at domain boundaries if true (this should be deprecated)

IF (CHECK_MASS_CONSERVE) CALL SET_DOMAINDIFFLX_3D(ZZP,RHO_D_DZDX,RHO_D_DZDY,RHO_D_DZDZ,.NOT.PREDICTOR)

! Get the specific heat

IF (.NOT.CONSTANT_SPECIFIC_HEAT_RATIO) THEN

   CP => WORK5
   R_H_G => WORK9

   !$OMP PARALLEL PRIVATE(ZZ_GET)
   ALLOCATE(ZZ_GET(1:N_TRACKED_SPECIES))
   !$OMP DO SCHEDULE(STATIC)
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
            ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I,J,K,1:N_TRACKED_SPECIES)
            CALL GET_SPECIFIC_HEAT(ZZ_GET,CP(I,J,K),TMP(I,J,K))
            R_H_G(I,J,K) = 1._EB/(CP(I,J,K)*TMP(I,J,K))
         ENDDO
      ENDDO
   ENDDO
   !$OMP END DO
   DEALLOCATE(ZZ_GET)
   !$OMP END PARALLEL

ENDIF

! Compute del dot k del T

KDTDX => WORK1
KDTDY => WORK2
KDTDZ => WORK3
KP    => WORK4

! Compute thermal conductivity k (KP)

K_DNS_OR_LES: IF (SIM_MODE==DNS_MODE .OR. SIM_MODE==LES_MODE) THEN

   ALLOCATE(ZZ_GET(1:N_TRACKED_SPECIES))
   KP = 0._EB
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
            ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I,J,K,1:N_TRACKED_SPECIES)
            CALL GET_CONDUCTIVITY(ZZ_GET,KP(I,J,K),TMP(I,J,K))
         ENDDO
      ENDDO
   ENDDO
   DEALLOCATE(ZZ_GET)

   IF (SIM_MODE==LES_MODE) THEN
      IF (.NOT.POTENTIAL_TEMPERATURE_CORRECTION) THEN
         ! normal LES mode, constant turbulent Prandtl number
         IF(.NOT.CONSTANT_SPECIFIC_HEAT_RATIO) THEN
            KP = KP + MAX(0._EB,(MU-MU_DNS))*CP*RPR
         ELSE
            KP = KP + MAX(0._EB,(MU-MU_DNS))*CPOPR
         ENDIF
      ELSE
         ! dynamic turbulent Prandtl number (Deardorff, 1980)
         IF(.NOT.CONSTANT_SPECIFIC_HEAT_RATIO) THEN
            KP = KP + MAX(0._EB,(MU-MU_DNS))*CP/PR_T
         ELSE
            KP = KP + MAX(0._EB,(MU-MU_DNS))*CPOPR*PR/PR_T
         ENDIF
      ENDIF
   ENDIF

   BOUNDARY_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS
      WC=>WALL(IW)
      II  = WC%ONE_D%II
      JJ  = WC%ONE_D%JJ
      KK  = WC%ONE_D%KK
      IIG = WC%ONE_D%IIG
      JJG = WC%ONE_D%JJG
      KKG = WC%ONE_D%KKG
      KP(II,JJ,KK) = KP(IIG,JJG,KKG)
   ENDDO BOUNDARY_LOOP

ELSE K_DNS_OR_LES

   ! normal VLES mode
   KP = MU*CPOPR

   ! dynamic turbulent Prandtl number (Deardorff, 1980)
   IF (POTENTIAL_TEMPERATURE_CORRECTION)  KP = KP*PR/PR_T

ENDIF K_DNS_OR_LES

! Store max diffusivity for stability check

IF (CHECK_VN .AND. .NOT.CONSTANT_SPECIFIC_HEAT_RATIO) THEN
   !$OMP PARALLEL DO SCHEDULE(STATIC)
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
            D_Z_MAX(I,J,K) = MAX(D_Z_MAX(I,J,K),KP(I,J,K)/(CP(I,J,K)*RHOP(I,J,K)))
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO
ENDIF

! Compute k*dT/dx, etc

!$OMP PARALLEL DO PRIVATE(DTDX, DTDY, DTDZ) SCHEDULE(STATIC)
DO K=0,KBAR
   DO J=0,JBAR
      DO I=0,IBAR
         DTDX = (TMP(I+1,J,K)-TMP(I,J,K))*RDXN(I)
         KDTDX(I,J,K) = .5_EB*(KP(I+1,J,K)+KP(I,J,K))*DTDX
         DTDY = (TMP(I,J+1,K)-TMP(I,J,K))*RDYN(J)
         KDTDY(I,J,K) = .5_EB*(KP(I,J+1,K)+KP(I,J,K))*DTDY
         DTDZ = (TMP(I,J,K+1)-TMP(I,J,K))*RDZN(K)
         KDTDZ(I,J,K) = .5_EB*(KP(I,J,K+1)+KP(I,J,K))*DTDZ
      ENDDO
   ENDDO
ENDDO
!$OMP END PARALLEL DO

! Correct thermal gradient (k dT/dn) at boundaries

CORRECTION_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
   WC => WALL(IW)
   IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY .OR. WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY) CYCLE CORRECTION_LOOP
   II  = WC%ONE_D%II
   JJ  = WC%ONE_D%JJ
   KK  = WC%ONE_D%KK
   IIG = WC%ONE_D%IIG
   JJG = WC%ONE_D%JJG
   KKG = WC%ONE_D%KKG
   IF (WC%BOUNDARY_TYPE==OPEN_BOUNDARY) THEN
      WC%ONE_D%K_G = 0.5_EB*(KP(IIG,JJG,KKG)+KP(II,JJ,KK))
      CYCLE CORRECTION_LOOP
   ELSE
      WC%ONE_D%K_G = KP(IIG,JJG,KKG)
   ENDIF
   IOR = WC%ONE_D%IOR
   SELECT CASE(IOR)
      CASE( 1)
         KDTDX(II,JJ,KK)   = 0._EB
      CASE(-1)
         KDTDX(II-1,JJ,KK) = 0._EB
      CASE( 2)
         KDTDY(II,JJ,KK)   = 0._EB
      CASE(-2)
         KDTDY(II,JJ-1,KK) = 0._EB
      CASE( 3)
         KDTDZ(II,JJ,KK)   = 0._EB
      CASE(-3)
         KDTDZ(II,JJ,KK-1) = 0._EB
   END SELECT
   ! Q_LEAK accounts for enthalpy moving through leakage paths
   DP(IIG,JJG,KKG) = DP(IIG,JJG,KKG) - ( WC%ONE_D%Q_CON_F*WC%ONE_D%RDN - WC%Q_LEAK )
ENDDO CORRECTION_LOOP

! Compute (q + del dot k del T) and add to the divergence

CYLINDER3: SELECT CASE(CYLINDRICAL)
CASE(.FALSE.) CYLINDER3   ! 3D or 2D Cartesian
   !$OMP PARALLEL DO PRIVATE(DELKDELT) SCHEDULE(STATIC)
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            DELKDELT = (KDTDX(I,J,K)-KDTDX(I-1,J,K))*RDX(I) + &
                       (KDTDY(I,J,K)-KDTDY(I,J-1,K))*RDY(J) + &
                       (KDTDZ(I,J,K)-KDTDZ(I,J,K-1))*RDZ(K)
            DP(I,J,K) = DP(I,J,K) + DELKDELT + Q(I,J,K) + QR(I,J,K)
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO
CASE(.TRUE.) CYLINDER3   ! 2D Cylindrical
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            DELKDELT = &
                 (R(I)*KDTDX(I,J,K)-R(I-1)*KDTDX(I-1,J,K))*RDX(I)*RRN(I) + &
                 (KDTDZ(I,J,K)-            KDTDZ(I,J,K-1))*RDZ(K)
            DP(I,J,K) = DP(I,J,K) + DELKDELT + Q(I,J,K) + QR(I,J,K)
         ENDDO
      ENDDO
   ENDDO
END SELECT CYLINDER3

! Point to the appropriate velocity components

IF (PREDICTOR) THEN
   UU=>U
   VV=>V
   WW=>W
ELSE
   UU=>US
   VV=>VS
   WW=>WS
ENDIF

IF (STORE_DIVERGENCE_CORRECTION) DCOR=0._EB

! Compute U_DOT_DEL_RHO_H_S and add to other enthalpy equation source terms

CONST_GAMMA_IF_1: IF (.NOT.CONSTANT_SPECIFIC_HEAT_RATIO) THEN

   CALL ENTHALPY_ADVECTION ! Compute u dot grad rho h_s

   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
            DP(I,J,K) = DP(I,J,K) - U_DOT_DEL_RHO_H_S(I,J,K)
            IF (STORE_DIVERGENCE_CORRECTION) DCOR(I,J,K) = - U_DOT_DEL_RHO_H_S(I,J,K)
         ENDDO
      ENDDO
   ENDDO

   IF (CC_IBM) CALL SET_EXIMRHOHSLIM_3D(NM) ! WORK2,WORK3,WORK4: Get flux limited \bar{rho Hs} on EXIM faces.

ENDIF CONST_GAMMA_IF_1

! Compute RTRM = 1/(rho*c_p*T) and multiply it by divergence terms already summed up

IF (CONSTANT_SPECIFIC_HEAT_RATIO) THEN

   !$OMP PARALLEL DO SCHEDULE(STATIC)
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
            IPZ = PRESSURE_ZONE(I,J,K)
            RTRM(I,J,K) = GM1OG*R_PBAR(K,IPZ)
            DP(I,J,K)   = RTRM(I,J,K)*DP(I,J,K)
        ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO

ELSE

   !$OMP PARALLEL DO SCHEDULE(STATIC)
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
            RTRM(I,J,K) = R_H_G(I,J,K)/RHOP(I,J,K)
            DP(I,J,K) = RTRM(I,J,K)*DP(I,J,K)
            IF (STORE_DIVERGENCE_CORRECTION) DCOR(I,J,K) = RTRM(I,J,K)*DCOR(I,J,K)
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO

ENDIF

! Compute (1/rho) * Sum( (Wbar/W_alpha-h_s,alpha/cp*T) (del dot rho*D del Z_n - u dot del rho*Z_n)

CONST_GAMMA_IF_2: IF (.NOT.CONSTANT_SPECIFIC_HEAT_RATIO) THEN

   DO N=1,N_TRACKED_SPECIES

      CALL SPECIES_ADVECTION ! Compute u dot grad rho Z_n

      SM  => SPECIES_MIXTURE(N)
      !$OMP PARALLEL DO PRIVATE(H_S) SCHEDULE(guided)

      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
               CALL GET_SENSIBLE_ENTHALPY_Z(N,TMP(I,J,K),H_S)
               DP(I,J,K) = DP(I,J,K) + (SM%RCON/RSUM(I,J,K) - H_S*R_H_G(I,J,K))* &
                    ( DEL_RHO_D_DEL_Z(I,J,K,N) - U_DOT_DEL_RHO_Z(I,J,K) )/RHOP(I,J,K)

               IF (STORE_DIVERGENCE_CORRECTION) THEN
                  DCOR(I,J,K) = DCOR(I,J,K) - (SM%RCON/RSUM(I,J,K) - H_S*R_H_G(I,J,K))*U_DOT_DEL_RHO_Z(I,J,K)/RHOP(I,J,K)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      !$OMP END PARALLEL DO

      IF (CC_IBM) CALL SET_EXIMRHOZZLIM_3D(NM,N) ! WORK2,WORK3,WORK4: flux limited \bar{rho Za} on EXIM faces.

   ENDDO

ENDIF CONST_GAMMA_IF_2

! Add contribution of reactions

IF (N_REACTIONS > 0 .OR. N_LP_ARRAY_INDICES>0 .OR. ANY(SPECIES_MIXTURE%DEPOSITING) .OR. &
    ANY(SPECIES_MIXTURE%CONDENSATION_SMIX_INDEX>0)) THEN
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            DP(I,J,K) = DP(I,J,K) + D_SOURCE(I,J,K)
         ENDDO
      ENDDO
   ENDDO
ENDIF

! Atmospheric stratification term

IF (STRATIFICATION) THEN
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
            DP(I,J,K) = DP(I,J,K) + RTRM(I,J,K)*0.5_EB*(WW(I,J,K)+WW(I,J,K-1))*RHO_0(K)*GVEC(3)
         ENDDO
      ENDDO
   ENDDO
ENDIF

! Manufactured solution

MMS_IF: IF (PERIODIC_TEST==7) THEN
   IF (PREDICTOR) TT=T+DT
   IF (CORRECTOR) TT=T
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            ! this term is similar to D_REACTION from fire
            XHAT = XC(I) - UF_MMS*TT
            ZHAT = ZC(K) - WF_MMS*TT
            DO N=1,N_TRACKED_SPECIES
               SM => SPECIES_MIXTURE(N)
               SELECT CASE(N)
                  CASE(1); Q_Z = -VD2D_MMS_Z_SRC(XHAT,ZHAT,TT)
                  CASE(2); Q_Z =  VD2D_MMS_Z_SRC(XHAT,ZHAT,TT)
               END SELECT
               CALL GET_SENSIBLE_ENTHALPY_Z(N,TMP(I,J,K),H_S)
               DP(I,J,K) = DP(I,J,K) + ( SM%RCON/RSUM(I,J,K) - H_S*R_H_G(I,J,K) )*Q_Z/RHOP(I,J,K)
            ENDDO
            ! debug
            !Q_Z = VD2D_MMS_Z_SRC(XHAT,ZHAT,TT)
            !DP(I,J,K) = (1._EB/RHO_1_MMS - 1._EB/RHO_0_MMS) * ( DEL_RHO_D_DEL_Z(I,J,K,2) + Q_Z )
         ENDDO
      ENDDO
   ENDDO
ENDIF MMS_IF

IF (CC_IBM .AND. .NOT.COMPUTE_CUTCELLS_ONLY) THEN
   T_USED(2)=T_USED(2)+CURRENT_TIME()-TNOW
   CALL CCREGION_DIVERGENCE_PART_1(T,DT,NM)
   TNOW=CURRENT_TIME()
ENDIF

1000 CONTINUE ! Evacuation meshes jump here

! Calculate pressure rise in each of the pressure zones by summing divergence expression over each zone

IF_PRESSURE_ZONES: IF (N_ZONE>0) THEN

   USUM(1:N_ZONE,NM) = 0._EB
   DSUM(1:N_ZONE,NM) = 0._EB
   PSUM(1:N_ZONE,NM) = 0._EB

   IF (EVACUATION_ONLY(NM)) RTRM=1._EB

   R_PFCT = 1._EB
   DO K=1,KBAR
      DO J=1,JBAR
         VC1 = DY(J)*DZ(K)
         DO I=1,IBAR
            IF (INTERPOLATED_MESH(I,J,K)>0) CYCLE
            IPZ = PRESSURE_ZONE(I,J,K)
            IF (IPZ<1) CYCLE
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
            VC = DX(I)*RC(I)*VC1
            IF (CC_IBM) THEN
               R_PFCT = 1._EB
               IF (CCVAR(I,J,K,IBM_CGSC) == IBM_SOLID) THEN
                  CYCLE
               ELSEIF(CCVAR(I,J,K,IBM_CGSC) == IBM_CUTCFE) THEN
                  ICC=CCVAR(I,J,K,IBM_IDCC)
                  R_PFCT = SUM(CUT_CELL(ICC)%VOLUME(1:CUT_CELL(ICC)%NCELL)) / VC
               ENDIF
            ENDIF
            DSUM(IPZ,NM) = DSUM(IPZ,NM) + VC*DP(I,J,K)
            PSUM(IPZ,NM) = PSUM(IPZ,NM) + VC*(R_PBAR(K,IPZ)*R_PFCT-RTRM(I,J,K))
         ENDDO
      ENDDO
   ENDDO

   ! Calculate the volume flux to the boundary of the pressure zone (int u dot dA)

   WALL_LOOP4: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
      WC => WALL(IW)
      IF (INTERPOLATED_MESH(WC%ONE_D%IIG,WC%ONE_D%JJG,WC%ONE_D%KKG)>0) CYCLE
      IPZ = WC%ONE_D%PRESSURE_ZONE
      IF (IPZ<1) CYCLE WALL_LOOP4
      IF (WC%BOUNDARY_TYPE/=SOLID_BOUNDARY) CYCLE WALL_LOOP4
      IF (PREDICTOR) USUM(IPZ,NM) = USUM(IPZ,NM) + WC%ONE_D%U_NORMAL_S*WALL(IW)%ONE_D%AREA
      IF (CORRECTOR) USUM(IPZ,NM) = USUM(IPZ,NM) + WC%ONE_D%U_NORMAL* WALL(IW)%ONE_D%AREA
   ENDDO WALL_LOOP4

ENDIF IF_PRESSURE_ZONES

T_USED(2)=T_USED(2)+CURRENT_TIME()-TNOW

CONTAINS


SUBROUTINE ENTHALPY_ADVECTION

REAL(EB), POINTER, DIMENSION(:,:,:) :: FX_H_S=>NULL(),FY_H_S=>NULL(),FZ_H_S=>NULL()
REAL(EB) :: UN_P,TMP_F_GAS

RHO_H_S_P=>WORK1
FX_H_S=>WORK2
FY_H_S=>WORK3
FZ_H_S=>WORK4
U_DOT_DEL_RHO_H_S=>WORK6
U_DOT_DEL_RHO_H_S=0._EB

IF (.NOT.ENTHALPY_TRANSPORT) RETURN

! Compute and store rho*h_s

!$OMP PARALLEL PRIVATE(ZZ_GET, H_S)
ALLOCATE(ZZ_GET(1:N_TRACKED_SPECIES))
!$OMP DO SCHEDULE(static)
DO K=0,KBP1
   DO J=0,JBP1
      DO I=0,IBP1
         ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I,J,K,1:N_TRACKED_SPECIES)
         CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,TMP(I,J,K))
         RHO_H_S_P(I,J,K) = RHOP(I,J,K)*H_S
      ENDDO
   ENDDO
ENDDO
!$OMP END DO
DEALLOCATE(ZZ_GET)
!$OMP END PARALLEL

! Compute scalar face values

!$OMP PARALLEL PRIVATE(ZZZ)
!$OMP DO SCHEDULE(STATIC)
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBM1
         ZZZ(1:4) = RHO_H_S_P(I-1:I+2,J,K)
         FX_H_S(I,J,K) = SCALAR_FACE_VALUE(UU(I,J,K),ZZZ,I_FLUX_LIMITER)
      ENDDO
   ENDDO
ENDDO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO K=1,KBAR
   DO J=1,JBM1
      DO I=1,IBAR
         ZZZ(1:4) = RHO_H_S_P(I,J-1:J+2,K)
         FY_H_S(I,J,K) = SCALAR_FACE_VALUE(VV(I,J,K),ZZZ,I_FLUX_LIMITER)
      ENDDO
   ENDDO
ENDDO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO K=1,KBM1
   DO J=1,JBAR
      DO I=1,IBAR
         ZZZ(1:4) = RHO_H_S_P(I,J,K-1:K+2)
         FZ_H_S(I,J,K) = SCALAR_FACE_VALUE(WW(I,J,K),ZZZ,I_FLUX_LIMITER)
      ENDDO
   ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

ALLOCATE(ZZ_GET(1:N_TRACKED_SPECIES))

WALL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
   WC=>WALL(IW)
   IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE WALL_LOOP

   II  = WC%ONE_D%II
   JJ  = WC%ONE_D%JJ
   KK  = WC%ONE_D%KK
   IIG = WC%ONE_D%IIG
   JJG = WC%ONE_D%JJG
   KKG = WC%ONE_D%KKG
   IOR = WC%ONE_D%IOR

   ! Handle the case where OBST lives on an external boundary
   IF (IW>N_EXTERNAL_WALL_CELLS) THEN
      SELECT CASE(IOR)
         CASE( 1); IF (IIG>IBAR) CYCLE WALL_LOOP
         CASE(-1); IF (IIG<1)    CYCLE WALL_LOOP
         CASE( 2); IF (JJG>JBAR) CYCLE WALL_LOOP
         CASE(-2); IF (JJG<1)    CYCLE WALL_LOOP
         CASE( 3); IF (KKG>KBAR) CYCLE WALL_LOOP
         CASE(-3); IF (KKG<1)    CYCLE WALL_LOOP
      END SELECT
   ENDIF

   ! Calculate the sensible enthalpy at the boundary. If the boundary is solid
   ! and the gas is flowing out, use the gas temperature for the calculation.

   IF (PREDICTOR) THEN
      UN_P = WC%ONE_D%U_NORMAL_S
   ELSE
      UN_P = WC%ONE_D%U_NORMAL
   ENDIF
   IF (WC%BOUNDARY_TYPE==SOLID_BOUNDARY .AND. UN_P>0._EB) THEN
      TMP_F_GAS = WC%ONE_D%TMP_G
   ELSE
      TMP_F_GAS = WC%ONE_D%TMP_F
   ENDIF

   ZZ_GET(1:N_TRACKED_SPECIES) = WC%ONE_D%ZZ_F(1:N_TRACKED_SPECIES)
   CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,TMP_F_GAS)

   ! overwrite first off-wall advective flux if flow is away from the wall and if the face is not also a wall cell

   IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY .AND. WC%BOUNDARY_TYPE/=OPEN_BOUNDARY) THEN

      OFF_WALL_SELECT_1: SELECT CASE(IOR)
         CASE( 1) OFF_WALL_SELECT_1
            !      ghost          FX/UU(II+1)
            ! ///   II   ///  II+1  |  II+2  | ...
            !                       ^ WALL_INDEX(II+1,+1)
            IF ((UU(II+1,JJ,KK)>0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II+1,JJ,KK),+1)>0)) THEN
               ZZZ(1:3) = (/WC%ONE_D%RHO_F*H_S,RHO_H_S_P(II+1:II+2,JJ,KK)/)
               FX_H_S(II+1,JJ,KK) = SCALAR_FACE_VALUE(UU(II+1,JJ,KK),ZZZ,I_FLUX_LIMITER)
            ENDIF
         CASE(-1) OFF_WALL_SELECT_1
            !            FX/UU(II-2)     ghost
            ! ... |  II-2  |  II-1  ///   II   ///
            !              ^ WALL_INDEX(II-1,-1)
            IF ((UU(II-2,JJ,KK)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II-1,JJ,KK),-1)>0)) THEN
               ZZZ(2:4) = (/RHO_H_S_P(II-2:II-1,JJ,KK),WC%ONE_D%RHO_F*H_S/)
               FX_H_S(II-2,JJ,KK) = SCALAR_FACE_VALUE(UU(II-2,JJ,KK),ZZZ,I_FLUX_LIMITER)
            ENDIF
         CASE( 2) OFF_WALL_SELECT_1
            IF ((VV(II,JJ+1,KK)>0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ+1,KK),+2)>0)) THEN
               ZZZ(1:3) = (/WC%ONE_D%RHO_F*H_S,RHO_H_S_P(II,JJ+1:JJ+2,KK)/)
               FY_H_S(II,JJ+1,KK) = SCALAR_FACE_VALUE(VV(II,JJ+1,KK),ZZZ,I_FLUX_LIMITER)
            ENDIF
         CASE(-2) OFF_WALL_SELECT_1
            IF ((VV(II,JJ-2,KK)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ-1,KK),-2)>0)) THEN
               ZZZ(2:4) = (/RHO_H_S_P(II,JJ-2:JJ-1,KK),WC%ONE_D%RHO_F*H_S/)
               FY_H_S(II,JJ-2,KK) = SCALAR_FACE_VALUE(VV(II,JJ-2,KK),ZZZ,I_FLUX_LIMITER)
            ENDIF
         CASE( 3) OFF_WALL_SELECT_1
            IF ((WW(II,JJ,KK+1)>0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ,KK+1),+3)>0)) THEN
               ZZZ(1:3) = (/WC%ONE_D%RHO_F*H_S,RHO_H_S_P(II,JJ,KK+1:KK+2)/)
               FZ_H_S(II,JJ,KK+1) = SCALAR_FACE_VALUE(WW(II,JJ,KK+1),ZZZ,I_FLUX_LIMITER)
            ENDIF
         CASE(-3) OFF_WALL_SELECT_1
            IF ((WW(II,JJ,KK-2)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ,KK-1),-3)>0)) THEN
               ZZZ(2:4) = (/RHO_H_S_P(II,JJ,KK-2:KK-1),WC%ONE_D%RHO_F*H_S/)
               FZ_H_S(II,JJ,KK-2) = SCALAR_FACE_VALUE(WW(II,JJ,KK-2),ZZZ,I_FLUX_LIMITER)
            ENDIF
      END SELECT OFF_WALL_SELECT_1

   ENDIF

   BOUNDARY_SELECT: SELECT CASE(WC%BOUNDARY_TYPE)
      CASE DEFAULT
         IOR_SELECT: SELECT CASE(IOR)
            CASE( 1); UN = UU(II,JJ,KK)
            CASE(-1); UN = UU(II-1,JJ,KK)
            CASE( 2); UN = VV(II,JJ,KK)
            CASE(-2); UN = VV(II,JJ-1,KK)
            CASE( 3); UN = WW(II,JJ,KK)
            CASE(-3); UN = WW(II,JJ,KK-1)
         END SELECT IOR_SELECT
      CASE(SOLID_BOUNDARY)
         IF (PREDICTOR) UN = -SIGN(1._EB,REAL(IOR,EB))*WC%ONE_D%U_NORMAL_S
         IF (CORRECTOR) UN = -SIGN(1._EB,REAL(IOR,EB))*WC%ONE_D%U_NORMAL
      CASE(INTERPOLATED_BOUNDARY)
         UN = UVW_SAVE(IW)
   END SELECT BOUNDARY_SELECT

   DU = (WC%ONE_D%RHO_F*H_S - RHO_H_S_P(IIG,JJG,KKG))*UN
   U_DOT_DEL_RHO_H_S(IIG,JJG,KKG) = U_DOT_DEL_RHO_H_S(IIG,JJG,KKG) - SIGN(1._EB,REAL(IOR,EB))*DU*WC%ONE_D%RDN

ENDDO WALL_LOOP

DEALLOCATE(ZZ_GET)

DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IF (SOLID(CELL_INDEX(I,J,K))) CYCLE

         DU_P = 0._EB
         DU_M = 0._EB
         IF (WALL_INDEX(CELL_INDEX(I,J,K),+1)==0) DU_P = (FX_H_S(I,J,K)   - RHO_H_S_P(I,J,K))*UU(I,J,K)   ! FDS Tech Guide (B.13)
         IF (WALL_INDEX(CELL_INDEX(I,J,K),-1)==0) DU_M = (FX_H_S(I-1,J,K) - RHO_H_S_P(I,J,K))*UU(I-1,J,K) ! FDS Tech Guide (B.14)
         U_DOT_DEL_RHO_H_S(I,J,K) = U_DOT_DEL_RHO_H_S(I,J,K) + (DU_P-DU_M)*RDX(I)                         ! FDS Tech Guide (B.12)

         DU_P = 0._EB
         DU_M = 0._EB
         IF (WALL_INDEX(CELL_INDEX(I,J,K),+2)==0) DU_P = (FY_H_S(I,J,K)   - RHO_H_S_P(I,J,K))*VV(I,J,K)
         IF (WALL_INDEX(CELL_INDEX(I,J,K),-2)==0) DU_M = (FY_H_S(I,J-1,K) - RHO_H_S_P(I,J,K))*VV(I,J-1,K)
         U_DOT_DEL_RHO_H_S(I,J,K) = U_DOT_DEL_RHO_H_S(I,J,K) + (DU_P-DU_M)*RDY(J)

         DU_P = 0._EB
         DU_M = 0._EB
         IF (WALL_INDEX(CELL_INDEX(I,J,K),+3)==0) DU_P = (FZ_H_S(I,J,K)   - RHO_H_S_P(I,J,K))*WW(I,J,K)
         IF (WALL_INDEX(CELL_INDEX(I,J,K),-3)==0) DU_M = (FZ_H_S(I,J,K-1) - RHO_H_S_P(I,J,K))*WW(I,J,K-1)
         U_DOT_DEL_RHO_H_S(I,J,K) = U_DOT_DEL_RHO_H_S(I,J,K) + (DU_P-DU_M)*RDZ(K)

      ENDDO
   ENDDO
ENDDO

END SUBROUTINE ENTHALPY_ADVECTION


SUBROUTINE SPECIES_ADVECTION

REAL(EB), POINTER, DIMENSION(:,:,:) :: FX_ZZ=>NULL(),FY_ZZ=>NULL(),FZ_ZZ=>NULL()

FX_ZZ=>WORK2
FY_ZZ=>WORK3
FZ_ZZ=>WORK4
RHO_Z_P=>WORK6

U_DOT_DEL_RHO_Z=>WORK7
U_DOT_DEL_RHO_Z=0._EB

IF (.NOT.ENTHALPY_TRANSPORT) RETURN

!$OMP PARALLEL PRIVATE(ZZZ)
!$OMP DO SCHEDULE(static) COLLAPSE(3)
DO K=0,KBP1
   DO J=0,JBP1
      DO I=0,IBP1
         RHO_Z_P(I,J,K) = RHOP(I,J,K)*ZZP(I,J,K,N)
      ENDDO
   ENDDO
ENDDO
!$OMP END DO

! Compute scalar face values

!$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBM1
         ZZZ(1:4) = RHO_Z_P(I-1:I+2,J,K)
         FX_ZZ(I,J,K) = SCALAR_FACE_VALUE(UU(I,J,K),ZZZ,I_FLUX_LIMITER)
      ENDDO
   ENDDO
ENDDO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
DO K=1,KBAR
   DO J=1,JBM1
      DO I=1,IBAR
         ZZZ(1:4) = RHO_Z_P(I,J-1:J+2,K)
         FY_ZZ(I,J,K) = SCALAR_FACE_VALUE(VV(I,J,K),ZZZ,I_FLUX_LIMITER)
      ENDDO
   ENDDO
ENDDO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
DO K=1,KBM1
   DO J=1,JBAR
      DO I=1,IBAR
         ZZZ(1:4) = RHO_Z_P(I,J,K-1:K+2)
         FZ_ZZ(I,J,K) = SCALAR_FACE_VALUE(WW(I,J,K),ZZZ,I_FLUX_LIMITER)
      ENDDO
   ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

WALL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
   WC=>WALL(IW)
   IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE WALL_LOOP

   II  = WC%ONE_D%II
   JJ  = WC%ONE_D%JJ
   KK  = WC%ONE_D%KK
   IIG = WC%ONE_D%IIG
   JJG = WC%ONE_D%JJG
   KKG = WC%ONE_D%KKG
   IOR = WC%ONE_D%IOR

   ! Handle the case where OBST lives on an external boundary
   IF (IW>N_EXTERNAL_WALL_CELLS) THEN
      SELECT CASE(IOR)
         CASE( 1); IF (IIG>IBAR) CYCLE WALL_LOOP
         CASE(-1); IF (IIG<1)    CYCLE WALL_LOOP
         CASE( 2); IF (JJG>JBAR) CYCLE WALL_LOOP
         CASE(-2); IF (JJG<1)    CYCLE WALL_LOOP
         CASE( 3); IF (KKG>KBAR) CYCLE WALL_LOOP
         CASE(-3); IF (KKG<1)    CYCLE WALL_LOOP
      END SELECT
   ENDIF

   ! Overwrite first off-wall advective flux if flow is away from the wall and if the face is not also a wall cell

   OFF_WALL_IF_2: IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY .AND. WC%BOUNDARY_TYPE/=OPEN_BOUNDARY) THEN

      OFF_WALL_SELECT_2: SELECT CASE(IOR)
         CASE( 1) OFF_WALL_SELECT_2
            !      ghost          FX/UU(II+1)
            ! ///   II   ///  II+1  |  II+2  | ...
            !                       ^ WALL_INDEX(II+1,+1)
            IF ((UU(II+1,JJ,KK)>0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II+1,JJ,KK),+1)>0)) THEN
               ZZZ(1:3) = (/WC%ONE_D%RHO_F*WC%ONE_D%ZZ_F(N),RHO_Z_P(II+1:II+2,JJ,KK)/)
               FX_ZZ(II+1,JJ,KK) = SCALAR_FACE_VALUE(UU(II+1,JJ,KK),ZZZ,I_FLUX_LIMITER)
            ENDIF
         CASE(-1) OFF_WALL_SELECT_2
            !            FX/UU(II-2)     ghost
            ! ... |  II-2  |  II-1  ///   II   ///
            !              ^ WALL_INDEX(II-1,-1)
            IF ((UU(II-2,JJ,KK)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II-1,JJ,KK),-1)>0)) THEN
               ZZZ(2:4) = (/RHO_Z_P(II-2:II-1,JJ,KK),WC%ONE_D%RHO_F*WC%ONE_D%ZZ_F(N)/)
               FX_ZZ(II-2,JJ,KK) = SCALAR_FACE_VALUE(UU(II-2,JJ,KK),ZZZ,I_FLUX_LIMITER)
            ENDIF
         CASE( 2) OFF_WALL_SELECT_2
            IF ((VV(II,JJ+1,KK)>0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ+1,KK),+2)>0)) THEN
               ZZZ(1:3) = (/WC%ONE_D%RHO_F*WC%ONE_D%ZZ_F(N),RHO_Z_P(II,JJ+1:JJ+2,KK)/)
               FY_ZZ(II,JJ+1,KK) = SCALAR_FACE_VALUE(VV(II,JJ+1,KK),ZZZ,I_FLUX_LIMITER)
            ENDIF
         CASE(-2) OFF_WALL_SELECT_2
            IF ((VV(II,JJ-2,KK)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ-1,KK),-2)>0)) THEN
               ZZZ(2:4) = (/RHO_Z_P(II,JJ-2:JJ-1,KK),WC%ONE_D%RHO_F*WC%ONE_D%ZZ_F(N)/)
               FY_ZZ(II,JJ-2,KK) = SCALAR_FACE_VALUE(VV(II,JJ-2,KK),ZZZ,I_FLUX_LIMITER)
            ENDIF
         CASE( 3) OFF_WALL_SELECT_2
            IF ((WW(II,JJ,KK+1)>0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ,KK+1),+3)>0)) THEN
               ZZZ(1:3) = (/WC%ONE_D%RHO_F*WC%ONE_D%ZZ_F(N),RHO_Z_P(II,JJ,KK+1:KK+2)/)
               FZ_ZZ(II,JJ,KK+1) = SCALAR_FACE_VALUE(WW(II,JJ,KK+1),ZZZ,I_FLUX_LIMITER)
            ENDIF
         CASE(-3) OFF_WALL_SELECT_2
            IF ((WW(II,JJ,KK-2)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ,KK-1),-3)>0)) THEN
               ZZZ(2:4) = (/RHO_Z_P(II,JJ,KK-2:KK-1),WC%ONE_D%RHO_F*WC%ONE_D%ZZ_F(N)/)
               FZ_ZZ(II,JJ,KK-2) = SCALAR_FACE_VALUE(WW(II,JJ,KK-2),ZZZ,I_FLUX_LIMITER)
            ENDIF
      END SELECT OFF_WALL_SELECT_2

   ENDIF OFF_WALL_IF_2

   ! Correct U_DOT_DEL_RHO_Z at the boundary

   BOUNDARY_SELECT: SELECT CASE(WC%BOUNDARY_TYPE)
      CASE DEFAULT
         IOR_SELECT: SELECT CASE(IOR)
            CASE( 1); UN = UU(II,JJ,KK)
            CASE(-1); UN = UU(II-1,JJ,KK)
            CASE( 2); UN = VV(II,JJ,KK)
            CASE(-2); UN = VV(II,JJ-1,KK)
            CASE( 3); UN = WW(II,JJ,KK)
            CASE(-3); UN = WW(II,JJ,KK-1)
         END SELECT IOR_SELECT
      CASE(SOLID_BOUNDARY)
         IF (PREDICTOR) UN = -SIGN(1._EB,REAL(IOR,EB))*WC%ONE_D%U_NORMAL_S
         IF (CORRECTOR) UN = -SIGN(1._EB,REAL(IOR,EB))*WC%ONE_D%U_NORMAL
      CASE(INTERPOLATED_BOUNDARY)
         UN = UVW_SAVE(IW)
   END SELECT BOUNDARY_SELECT

   DU = (WC%ONE_D%RHO_F*WC%ONE_D%ZZ_F(N) - RHO_Z_P(IIG,JJG,KKG))*UN
   U_DOT_DEL_RHO_Z(IIG,JJG,KKG) = U_DOT_DEL_RHO_Z(IIG,JJG,KKG) - SIGN(1._EB,REAL(IOR,EB))*DU*WC%ONE_D%RDN

ENDDO WALL_LOOP

!$OMP PARALLEL DO PRIVATE(DU_P, DU_M) SCHEDULE(static)
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IF (SOLID(CELL_INDEX(I,J,K))) CYCLE

         DU_P = 0._EB
         DU_M = 0._EB
         IF (WALL_INDEX(CELL_INDEX(I,J,K),+1)==0) DU_P = (FX_ZZ(I,J,K)   - RHO_Z_P(I,J,K))*UU(I,J,K)
         IF (WALL_INDEX(CELL_INDEX(I,J,K),-1)==0) DU_M = (FX_ZZ(I-1,J,K) - RHO_Z_P(I,J,K))*UU(I-1,J,K)
         U_DOT_DEL_RHO_Z(I,J,K) = U_DOT_DEL_RHO_Z(I,J,K) + (DU_P-DU_M)*RDX(I)

         DU_P = 0._EB
         DU_M = 0._EB
         IF (WALL_INDEX(CELL_INDEX(I,J,K),+2)==0) DU_P = (FY_ZZ(I,J,K)   - RHO_Z_P(I,J,K))*VV(I,J,K)
         IF (WALL_INDEX(CELL_INDEX(I,J,K),-2)==0) DU_M = (FY_ZZ(I,J-1,K) - RHO_Z_P(I,J,K))*VV(I,J-1,K)
         U_DOT_DEL_RHO_Z(I,J,K) = U_DOT_DEL_RHO_Z(I,J,K) + (DU_P-DU_M)*RDY(J)

         DU_P = 0._EB
         DU_M = 0._EB
         IF (WALL_INDEX(CELL_INDEX(I,J,K),+3)==0) DU_P = (FZ_ZZ(I,J,K)   - RHO_Z_P(I,J,K))*WW(I,J,K)
         IF (WALL_INDEX(CELL_INDEX(I,J,K),-3)==0) DU_M = (FZ_ZZ(I,J,K-1) - RHO_Z_P(I,J,K))*WW(I,J,K-1)
         U_DOT_DEL_RHO_Z(I,J,K) = U_DOT_DEL_RHO_Z(I,J,K) + (DU_P-DU_M)*RDZ(K)

      ENDDO
   ENDDO
ENDDO
!$OMP END PARALLEL DO

END SUBROUTINE SPECIES_ADVECTION


SUBROUTINE MERGE_PRESSURE_ZONES

CONNECTED_ZONES(:,:,NM) = .FALSE.

IF (EVACUATION_ONLY(NM)) THEN
   CALL EVACUATION_PRESSURE_ZONES
   RETURN
END IF

DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
   WC=>WALL(IW)
   IF (WC%BOUNDARY_TYPE/=NULL_BOUNDARY .AND. WC%BOUNDARY_TYPE/=OPEN_BOUNDARY .AND. &
      WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY) CYCLE
   II  = WC%ONE_D%II
   JJ  = WC%ONE_D%JJ
   KK  = WC%ONE_D%KK
   IIG = WC%ONE_D%IIG
   JJG = WC%ONE_D%JJG
   KKG = WC%ONE_D%KKG
   IF (SOLID(CELL_INDEX(IIG,JJG,KKG))) CYCLE
   IPZ  = PRESSURE_ZONE(IIG,JJG,KKG)
   IOPZ = PRESSURE_ZONE(II,JJ,KK)
   IF (IW>N_EXTERNAL_WALL_CELLS .AND. IPZ/=IOPZ) THEN
      CONNECTED_ZONES(IOPZ,IPZ,NM) = .TRUE.
      CONNECTED_ZONES(IPZ,IOPZ,NM) = .TRUE.
   ENDIF
   IF (WC%BOUNDARY_TYPE==OPEN_BOUNDARY) THEN
      CONNECTED_ZONES(0,IPZ,NM) = .TRUE.
      CONNECTED_ZONES(IPZ,0,NM) = .TRUE.
   ELSEIF (WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY) THEN
      CONNECTED_ZONES(IOPZ,IPZ,NM) = .TRUE.
      CONNECTED_ZONES(IPZ,IOPZ,NM) = .TRUE.
   ENDIF
ENDDO

END SUBROUTINE MERGE_PRESSURE_ZONES


SUBROUTINE PREDICT_NORMAL_VELOCITY

PREDICT_NORMALS: IF (PREDICTOR) THEN

   WALL_LOOP3: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS

      WC => WALL(IW)
      IOR = WC%ONE_D%IOR

      WALL_CELL_TYPE: SELECT CASE (WC%BOUNDARY_TYPE)

         CASE (NULL_BOUNDARY)

            WC%ONE_D%U_NORMAL_S = 0._EB

         CASE (SOLID_BOUNDARY)

            SF => SURFACE(WC%SURF_INDEX)

            IF (SF%SPECIES_BC_INDEX==SPECIFIED_MASS_FLUX .OR. &
                SF%SPECIES_BC_INDEX==INTERPOLATED_BC     .OR. &
                WC%ONE_D%NODE_INDEX > 0                  .OR. &
                ANY(SF%LEAK_PATH>0))                          &
                CYCLE WALL_LOOP3

            IF (ABS(WC%ONE_D%T_IGN-T_BEGIN) < SPACING(WC%ONE_D%T_IGN) .AND. SF%RAMP_INDEX(TIME_VELO)>=1) THEN
               TSI = T + DT
            ELSE
               TSI = T + DT - WC%ONE_D%T_IGN
               IF (TSI<0._EB) THEN
                  WC%ONE_D%U_NORMAL_S = 0._EB
                  CYCLE WALL_LOOP3
               ENDIF
            ENDIF
            TIME_RAMP_FACTOR = EVALUATE_RAMP(TSI,SF%TAU(TIME_VELO),SF%RAMP_INDEX(TIME_VELO))
            KK               = WC%ONE_D%KK
            DELTA_P          = PBAR_P(KK,SF%DUCT_PATH(1)) - PBAR_P(KK,SF%DUCT_PATH(2))
            PRES_RAMP_FACTOR = SIGN(1._EB,SF%MAX_PRESSURE-DELTA_P)*SQRT(ABS((DELTA_P-SF%MAX_PRESSURE)/SF%MAX_PRESSURE))
            SELECT CASE(IOR)
               CASE( 1)
                  WC%ONE_D%U_NORMAL_S =-U0 + TIME_RAMP_FACTOR*(WC%ONE_D%U_NORMAL_0+U0)
               CASE(-1)
                  WC%ONE_D%U_NORMAL_S = U0 + TIME_RAMP_FACTOR*(WC%ONE_D%U_NORMAL_0-U0)
               CASE( 2)
                  WC%ONE_D%U_NORMAL_S =-V0 + TIME_RAMP_FACTOR*(WC%ONE_D%U_NORMAL_0+V0)
               CASE(-2)
                  WC%ONE_D%U_NORMAL_S = V0 + TIME_RAMP_FACTOR*(WC%ONE_D%U_NORMAL_0-V0)
               CASE( 3)
                  WC%ONE_D%U_NORMAL_S =-W0 + TIME_RAMP_FACTOR*(WC%ONE_D%U_NORMAL_0+W0)
               CASE(-3)
                  WC%ONE_D%U_NORMAL_S = W0 + TIME_RAMP_FACTOR*(WC%ONE_D%U_NORMAL_0-W0)
            END SELECT
            ! Special Cases
            NEUMANN_IF: IF (SF%SPECIFIED_NORMAL_GRADIENT) THEN
               IIG = WC%ONE_D%IIG
               JJG = WC%ONE_D%JJG
               KKG = WC%ONE_D%KKG
               SELECT CASE(IOR)
                  CASE( 1)
                     WC%ONE_D%U_NORMAL_S =-(U(IIG,JJG,KKG)   + SF%VEL_GRAD*WC%ONE_D%RDN)
                  CASE(-1)
                     WC%ONE_D%U_NORMAL_S = (U(IIG-1,JJG,KKG) + SF%VEL_GRAD*WC%ONE_D%RDN)
                  CASE( 2)
                     WC%ONE_D%U_NORMAL_S =-(V(IIG,JJG,KKG)   + SF%VEL_GRAD*WC%ONE_D%RDN)
                  CASE(-2)
                     WC%ONE_D%U_NORMAL_S = (V(IIG,JJG-1,KKG) + SF%VEL_GRAD*WC%ONE_D%RDN)
                  CASE( 3)
                     WC%ONE_D%U_NORMAL_S =-(W(IIG,JJG,KKG)   + SF%VEL_GRAD*WC%ONE_D%RDN)
                  CASE(-3)
                     WC%ONE_D%U_NORMAL_S = (W(IIG,JJG,KKG-1) + SF%VEL_GRAD*WC%ONE_D%RDN)
               END SELECT
            ENDIF NEUMANN_IF
            IF (ABS(SURFACE(WC%SURF_INDEX)%MASS_FLUX_TOTAL)>=TWO_EPSILON_EB) WC%ONE_D%U_NORMAL_S = &
                                                                             WC%ONE_D%U_NORMAL_S*RHOA/WC%ONE_D%RHO_F
            VENT_IF: IF (WC%VENT_INDEX>0) THEN
               VT=>VENTS(WC%VENT_INDEX)
               IF (VT%N_EDDY>0) THEN ! Synthetic Eddy Method
                  II = WC%ONE_D%II
                  JJ = WC%ONE_D%JJ
                  KK = WC%ONE_D%KK
                  IF (SF%PROFILE/=0 .AND. ABS(SF%VEL)>TWO_EPSILON_EB) THEN
                     PROFILE_FACTOR = ABS(WC%ONE_D%U_NORMAL_0/SF%VEL)
                  ELSE
                     PROFILE_FACTOR = 1._EB
                  ENDIF
                  SELECT CASE(IOR)
                     CASE( 1)
                        WC%ONE_D%U_NORMAL_S = WC%ONE_D%U_NORMAL_S - TIME_RAMP_FACTOR*VT%U_EDDY(JJ,KK)*PROFILE_FACTOR
                     CASE(-1)
                        WC%ONE_D%U_NORMAL_S = WC%ONE_D%U_NORMAL_S + TIME_RAMP_FACTOR*VT%U_EDDY(JJ,KK)*PROFILE_FACTOR
                     CASE( 2)
                        WC%ONE_D%U_NORMAL_S = WC%ONE_D%U_NORMAL_S - TIME_RAMP_FACTOR*VT%V_EDDY(II,KK)*PROFILE_FACTOR
                     CASE(-2)
                        WC%ONE_D%U_NORMAL_S = WC%ONE_D%U_NORMAL_S + TIME_RAMP_FACTOR*VT%V_EDDY(II,KK)*PROFILE_FACTOR
                     CASE( 3)
                        WC%ONE_D%U_NORMAL_S = WC%ONE_D%U_NORMAL_S - TIME_RAMP_FACTOR*VT%W_EDDY(II,JJ)*PROFILE_FACTOR
                     CASE(-3)
                        WC%ONE_D%U_NORMAL_S = WC%ONE_D%U_NORMAL_S + TIME_RAMP_FACTOR*VT%W_EDDY(II,JJ)*PROFILE_FACTOR
                  END SELECT
               ENDIF
               EVACUATION_BC_IF: IF (EVACUATION_ONLY(NM)) THEN
                  II = EVAC_TIME_ITERATIONS / MAXVAL(EMESH_NFIELDS)
                  IF ((ABS(ICYC)+1) > (WC%VENT_INDEX-1)*II .AND. (ABS(ICYC)+1) <= WC%VENT_INDEX*II) THEN
                     TSI = T + DT - (MAXVAL(EMESH_NFIELDS)-WC%VENT_INDEX)*II*EVAC_DT_FLOWFIELD
                     TIME_RAMP_FACTOR = EVALUATE_RAMP(TSI,SF%TAU(TIME_VELO),SF%RAMP_INDEX(TIME_VELO))
                  ELSE
                     TIME_RAMP_FACTOR = 0.0_EB
                  END IF
                  WC%ONE_D%U_NORMAL_S = TIME_RAMP_FACTOR*WC%ONE_D%U_NORMAL_0
               END IF EVACUATION_BC_IF
               WANNIER_BC: IF (PERIODIC_TEST==5) THEN
                  II = WC%ONE_D%II
                  JJ = WC%ONE_D%JJ
                  KK = WC%ONE_D%KK
                  SELECT CASE(IOR)
                     CASE( 1)
                        WC%ONE_D%U_NORMAL_S = -WANNIER_FLOW(X(II),ZC(KK),1)
                     CASE(-1)
                        WC%ONE_D%U_NORMAL_S =  WANNIER_FLOW(X(II-1),ZC(KK),1)
                     CASE( 3)
                        WC%ONE_D%U_NORMAL_S = -WANNIER_FLOW(XC(II),Z(KK),2)
                     CASE(-3)
                        WC%ONE_D%U_NORMAL_S =  WANNIER_FLOW(XC(II),Z(KK-1),2)
                  END SELECT
               ENDIF WANNIER_BC
            ENDIF VENT_IF

         CASE(OPEN_BOUNDARY,INTERPOLATED_BOUNDARY)

            II = WC%ONE_D%II
            JJ = WC%ONE_D%JJ
            KK = WC%ONE_D%KK
            SELECT CASE(IOR)
               CASE( 1)
                  WC%ONE_D%U_NORMAL_S = -U(II,JJ,KK)
               CASE(-1)
                  WC%ONE_D%U_NORMAL_S =  U(II-1,JJ,KK)
               CASE( 2)
                  WC%ONE_D%U_NORMAL_S = -V(II,JJ,KK)
               CASE(-2)
                  WC%ONE_D%U_NORMAL_S =  V(II,JJ-1,KK)
               CASE( 3)
                  WC%ONE_D%U_NORMAL_S = -W(II,JJ,KK)
               CASE(-3)
                  WC%ONE_D%U_NORMAL_S =  W(II,JJ,KK-1)
            END SELECT

      END SELECT WALL_CELL_TYPE

   ENDDO WALL_LOOP3

   ! Calculate du/dt, dv/dt, dw/dt at external boundaries. DUNDT is only used for Neumann BCs.

   !$OMP PARALLEL DO SCHEDULE(STATIC)
   DO IW=1,N_EXTERNAL_WALL_CELLS
      WALL(IW)%DUNDT = RDT*(WALL(IW)%ONE_D%U_NORMAL_S-WALL(IW)%ONE_D%U_NORMAL)
   ENDDO
   !$OMP END PARALLEL DO

ELSE PREDICT_NORMALS

   ! In the CORRECTOR step, the normal component of velocity, U_NORMAL, is the same as the predicted value, U_NORMAL_S.
   ! However, for species mass fluxes and HVAC, U_NORMAL is computed elsewhere (wall.f90).

   DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
      WC => WALL(IW)
      IF (WC%BOUNDARY_TYPE==SOLID_BOUNDARY) THEN
         SF => SURFACE(WC%SURF_INDEX)
         IF (SF%SPECIES_BC_INDEX==SPECIFIED_MASS_FLUX .OR. &
             SF%SPECIES_BC_INDEX==INTERPOLATED_BC     .OR. &
             WC%ONE_D%NODE_INDEX > 0                  .OR. &
             ANY(SF%LEAK_PATH>0)) CYCLE
      ENDIF
      WC%ONE_D%U_NORMAL = WC%ONE_D%U_NORMAL_S
   ENDDO

   ! Calculate du/dt, dv/dt, dw/dt at external boundaries. DUNDT is only used for Neumann BCs. Note the second-order RK formula.

   !$OMP PARALLEL DO SCHEDULE(STATIC)
   DO IW=1,N_EXTERNAL_WALL_CELLS
      WALL(IW)%DUNDT = WALL(IW)%DUNDT + 2._EB*RDT*(WALL(IW)%ONE_D%U_NORMAL-WALL(IW)%ONE_D%U_NORMAL_S)
   ENDDO
   !$OMP END PARALLEL DO

ENDIF PREDICT_NORMALS

END SUBROUTINE PREDICT_NORMAL_VELOCITY


SUBROUTINE EVACUATION_PRESSURE_ZONES
IMPLICIT NONE
REAL(EB) :: X1,Y1,X2,Y2,Z1,Z2
INTEGER :: IW,N,II,JJ,KK,IIG,JJG,KKG,ITMP,I,J,K,I_VENT,N_OVERLAP

IF (.NOT.EVACUATION_ONLY(NM)) RETURN

! Evacuation flow field calculation: Change the outflow vent pressure zone and initialize everything

EVACUATION_PREDICTOR: IF (PREDICTOR) THEN
   ITMP = EVAC_TIME_ITERATIONS / MAXVAL(EMESH_NFIELDS)
   EVACUATION_NEW_FIELD: IF (MOD(ICYC-1,ITMP) == 0 .AND. ICYC < 0) THEN
      ! New exit/door flow field calculation (new outflow-vent), do the necessary initializaions,
      ! because same arrays are used for exits/doors at a main evacuation mesh.  Each evacuation
      ! flow field calculation has just one exit/door, i.e., outflow vent, so the pressure zone
      ! is defined so that the front of the door is in the pressure zone.  So, pressure zone
      ! is also redefined for this main evacuation mesh.
      ! One pressure zone is defined for each main evacuation mesh.

      I_VENT = 0
      ITMP = (ABS(ICYC)+1)/ITMP
      FIND_EXIT_LOOP: DO I=1,N_EXITS-N_CO_EXITS+N_DOORS
         IF (.NOT.EMESH_EXITS(I)%MAINMESH==NM) CYCLE
         IF (.NOT.EMESH_EXITS(I)%DEFINE_MESH) CYCLE
         I_VENT = I_VENT + 1
         IF (I_VENT==ITMP) THEN
            X1 = EMESH_EXITS(I)%XB(1); X2 = EMESH_EXITS(I)%XB(2)
            Y1 = EMESH_EXITS(I)%XB(3); Y2 = EMESH_EXITS(I)%XB(4)
            Z1 = EMESH_EXITS(I)%XB(5); Z2 = EMESH_EXITS(I)%XB(6)
            SELECT CASE (EMESH_EXITS(I)%IOR)
            CASE(+1)
               X1 = EMESH_EXITS(I)%XB(1) - MESHES(EMESH_EXITS(I)%MAINMESH)%DXI
            CASE(-1)
               X2 = EMESH_EXITS(I)%XB(2) + MESHES(EMESH_EXITS(I)%MAINMESH)%DXI
            CASE(+2)
               Y1 = EMESH_EXITS(I)%XB(3) - MESHES(EMESH_EXITS(I)%MAINMESH)%DETA
            CASE(-2)
               Y2 = EMESH_EXITS(I)%XB(4) + MESHES(EMESH_EXITS(I)%MAINMESH)%DETA
            END SELECT
            EXIT FIND_EXIT_LOOP
         END IF
      END DO FIND_EXIT_LOOP
      N = 0
      DO I = 1,N_ZONE
         IF (.NOT.(P_ZONE(I)%EVACUATION)) CYCLE
         IF (P_ZONE(I)%MESH_INDEX==NM) THEN
            N = I ! The ordinar number of the pressure zone of this main evacuation mesh
            EXIT
         END IF
      END DO
      IF (N==0) THEN
         WRITE(LU_ERR,'(A,A)') 'ERROR FDS+Evac: Zone error, no pressure zone found for mesh ',TRIM(MESH_NAME(NM))
      END IF

      U=0._EB; V=0._EB; W=0._EB; US=0._EB; VS=0._EB; WS=0._EB; FVX=0._EB; FVY=0._EB; FVZ=0._EB
      H=0._EB; HS=0._EB; KRES=0._EB; DDDT=0._EB; D=0._EB; DS=0._EB
      P_0=P_INF; TMP_0=TMPA
      PBAR=P_INF; PBAR_S=P_INF; R_PBAR=0._EB; D_PBAR_DT=0._EB; D_PBAR_DT_S=0._EB
      RHO=RHO_0(1); RHOS=RHO_0(1); TMP=TMPA
      USUM(:,NM) = 0.0_EB ; DSUM(:,NM) = 0.0_EB; PSUM(:,NM) = 0.0_EB
      PRESSURE_ZONE = 0

      DO K=0,KBP1
         DO J=0,JBP1
            DO I=0,IBP1
               IF (PRESSURE_ZONE(I,J,K)==N) CYCLE
               IF (XC(I) - X1 >=0._EB .AND. XC(I) < X2 .AND. &
                    YC(J) - Y1 >=0._EB .AND. YC(J) < Y2 .AND. &
                    ZC(K) - Z1 >=0._EB .AND. ZC(K) < Z2) THEN
                  PRESSURE_ZONE(I,J,K) = N
                  IF (.NOT.SOLID(CELL_INDEX(I,J,K))) CALL ASSIGN_PRESSURE_ZONE(NM,XC(I),YC(J),ZC(K),N,N_OVERLAP)
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
         WC=>WALL(IW)
         WC%ONE_D%PRESSURE_ZONE = 0
         WC%ONE_D%U_NORMAL = 0._EB
         WC%ONE_D%U_TAU = 0._EB
         WC%ONE_D%RHO_F = RHO_0(1)
         WC%ONE_D%Y_PLUS = 1._EB
         WC%ONE_D%RHO_D_F = 0.1_EB ! Do not initialize to zero to avoid divide by zero in the first time step
         II  = WC%ONE_D%II
         JJ  = WC%ONE_D%JJ
         KK  = WC%ONE_D%KK
         IIG = WC%ONE_D%IIG
         JJG = WC%ONE_D%JJG
         KKG = WC%ONE_D%KKG
         IF (KK==1) WC%ONE_D%PRESSURE_ZONE = PRESSURE_ZONE(IIG,JJG,KKG)
      END DO
   END IF EVACUATION_NEW_FIELD
END IF EVACUATION_PREDICTOR

END SUBROUTINE EVACUATION_PRESSURE_ZONES


END SUBROUTINE DIVERGENCE_PART_1


SUBROUTINE DIVERGENCE_PART_2(DT,NM)

! Finish computing the divergence of the flow, D, and then compute its time derivative, DDDT

USE COMP_FUNCTIONS, ONLY: CURRENT_TIME
USE COMPLEX_GEOMETRY, ONLY : IBM_CGSC, IBM_IDCC, IBM_SOLID, IBM_CUTCFE

INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: DT
REAL(EB), POINTER, DIMENSION(:,:,:) :: DP,D_NEW,RTRM,DIV
REAL(EB) :: USUM_ADD(N_ZONE),UN_P
REAL(EB) :: RDT,TNOW,P_EQ,SUM_P_PSUM,SUM_USUM,SUM_DSUM,SUM_PSUM,R_PFCT
LOGICAL :: OPEN_ZONE
REAL(EB), POINTER, DIMENSION(:) :: D_PBAR_DT_P
REAL(EB), POINTER, DIMENSION(:,:) :: PBAR_P
INTEGER :: IW,IOR,II,JJ,KK,IIG,JJG,KKG,IC,I,J,K,IPZ,IOPZ,ICC
TYPE(WALL_TYPE), POINTER :: WC=>NULL()

IF (SOLID_PHASE_ONLY) RETURN
IF (PERIODIC_TEST==3) RETURN
IF (PERIODIC_TEST==4) RETURN

TNOW=CURRENT_TIME()
CALL POINT_TO_MESH(NM)

RDT = 1._EB/DT

SELECT CASE(PREDICTOR)
   CASE(.TRUE.)
      DP => DS
      PBAR_P => PBAR_S
   CASE(.FALSE.)
      DP => DDDT
      PBAR_P => PBAR
END SELECT

R_PBAR = 1._EB/PBAR_P

RTRM => WORK1

! Adjust volume flows (USUM) of pressure ZONEs that are connected to equalize background pressure

USUM_ADD = 0._EB

DO IPZ=1,N_ZONE
   IF (EVACUATION_ONLY(NM) .OR. P_ZONE(IPZ)%EVACUATION .OR. P_ZONE(IPZ)%PERIODIC) CYCLE
   SUM_P_PSUM = PBAR_P(1,IPZ)*PSUM(IPZ,NM)
   OPEN_ZONE  = .FALSE.
   SUM_USUM = USUM(IPZ,NM)
   SUM_DSUM = DSUM(IPZ,NM)
   SUM_PSUM = PSUM(IPZ,NM)
   DO IOPZ=N_ZONE,0,-1
      IF (IOPZ==IPZ) CYCLE
      IF (CONNECTED_ZONES(IPZ,IOPZ,NM)) THEN
         IF (IOPZ==0) THEN
            OPEN_ZONE = .TRUE.
         ELSE
            SUM_P_PSUM = SUM_P_PSUM + PBAR_P(1,IOPZ)*PSUM(IOPZ,NM)
            SUM_USUM = SUM_USUM + USUM(IOPZ,NM)
            SUM_DSUM = SUM_DSUM + DSUM(IOPZ,NM)
            SUM_PSUM = SUM_PSUM + PSUM(IOPZ,NM)
         ENDIF
      ENDIF
   ENDDO
   IF (OPEN_ZONE) THEN
      P_EQ          = P_0(1)
      USUM_ADD(IPZ) = PSUM(IPZ,NM)*(PBAR_P(1,IPZ)-P_EQ)/PRESSURE_RELAX_TIME + DSUM(IPZ,NM) - USUM(IPZ,NM)
   ELSE
      P_EQ          = SUM_P_PSUM/SUM_PSUM
      USUM_ADD(IPZ) = PSUM(IPZ,NM)*(PBAR_P(1,IPZ)-P_EQ)/PRESSURE_RELAX_TIME + DSUM(IPZ,NM) - USUM(IPZ,NM) - &
                      PSUM(IPZ,NM)*(SUM_DSUM-SUM_USUM)/SUM_PSUM
   ENDIF
ENDDO

DO IPZ=1,N_ZONE
   USUM(IPZ,NM) = USUM(IPZ,NM) + USUM_ADD(IPZ)
ENDDO

! Compute dP/dt for each pressure ZONE

IF_PRESSURE_ZONES: IF (N_ZONE>0) THEN

   IF (PREDICTOR) D_PBAR_DT_P => D_PBAR_DT_S
   IF (CORRECTOR) D_PBAR_DT_P => D_PBAR_DT

   ! Compute change in background pressure

   DO IPZ=1,N_ZONE
      IF (ABS(PSUM(IPZ,NM)) > TWO_EPSILON_EB) D_PBAR_DT_P(IPZ) = (DSUM(IPZ,NM) - USUM(IPZ,NM))/PSUM(IPZ,NM)
      IF (CORRECTOR) P_ZONE(IPZ)%DPSTAR =  D_PBAR_DT_P(IPZ)
   ENDDO

   ! Add pressure derivative to divergence

   R_PFCT = 1._EB
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IPZ = PRESSURE_ZONE(I,J,K)
            IF (IPZ<1) CYCLE
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
            IF (CC_IBM) THEN
               R_PFCT = 1._EB
               IF (CCVAR(I,J,K,IBM_CGSC) == IBM_SOLID) THEN
                  CYCLE
               ELSEIF(CCVAR(I,J,K,IBM_CGSC) == IBM_CUTCFE) THEN
                  ICC=CCVAR(I,J,K,IBM_IDCC)
                  R_PFCT = SUM(CUT_CELL(ICC)%VOLUME(1:CUT_CELL(ICC)%NCELL)) / (DX(I)*RC(I)*DY(J)*DZ(K))
               ENDIF
            ENDIF
            DP(I,J,K) = DP(I,J,K) - (R_PBAR(K,IPZ)*R_PFCT-RTRM(I,J,K))*D_PBAR_DT_P(IPZ)
         ENDDO
      ENDDO
   ENDDO

ENDIF IF_PRESSURE_ZONES

! Zero out divergence in solid cells

SOLID_LOOP: DO IC=1,CELL_COUNT(NM)
   IF (.NOT.SOLID(IC)) CYCLE SOLID_LOOP
   I = I_CELL(IC)
   J = J_CELL(IC)
   K = K_CELL(IC)
   DP(I,J,K) = 0._EB
ENDDO SOLID_LOOP

! Zero out CC_IBM solid cells:
IF (CC_IBM) THEN
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (CCVAR(I,J,K,IBM_CGSC) /= IBM_SOLID) CYCLE
            DP(I,J,K) = 0._EB
         ENDDO
      ENDDO
   ENDDO
ENDIF

! Specify divergence in boundary cells to account for volume being generated at the walls

BC_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
   WC => WALL(IW)
   IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE BC_LOOP
   II = WC%ONE_D%II
   JJ = WC%ONE_D%JJ
   KK = WC%ONE_D%KK
   SELECT CASE (WC%BOUNDARY_TYPE)
      CASE (SOLID_BOUNDARY)
         IF (.NOT.SOLID(CELL_INDEX(II,JJ,KK))) CYCLE BC_LOOP
         IOR = WC%ONE_D%IOR
         IF (PREDICTOR) THEN
            UN_P = WC%ONE_D%U_NORMAL_S
         ELSE
            UN_P = WC%ONE_D%U_NORMAL
         ENDIF
         SELECT CASE(IOR)
            CASE( 1)
               DP(II,JJ,KK) = DP(II,JJ,KK) - UN_P*RDX(II)*RRN(II)*R(II)
            CASE(-1)
               DP(II,JJ,KK) = DP(II,JJ,KK) - UN_P*RDX(II)*RRN(II)*R(II-1)
            CASE( 2)
               DP(II,JJ,KK) = DP(II,JJ,KK) - UN_P*RDY(JJ)
            CASE(-2)
               DP(II,JJ,KK) = DP(II,JJ,KK) - UN_P*RDY(JJ)
            CASE( 3)
               DP(II,JJ,KK) = DP(II,JJ,KK) - UN_P*RDZ(KK)
            CASE(-3)
               DP(II,JJ,KK) = DP(II,JJ,KK) - UN_P*RDZ(KK)
         END SELECT
      CASE (OPEN_BOUNDARY,MIRROR_BOUNDARY,INTERPOLATED_BOUNDARY)
         IIG = WC%ONE_D%IIG
         JJG = WC%ONE_D%JJG
         KKG = WC%ONE_D%KKG
         DP(II,JJ,KK) = DP(IIG,JJG,KKG)
   END SELECT
ENDDO BC_LOOP

! Compute time derivative of the divergence, dD/dt

TRUE_PROJECTION: IF (PROJECTION) THEN

   DIV=>WORK1

   IF (PREDICTOR) THEN
      DO K = 1,KBAR
         DO J = 1,JBAR
            DO I = 1,IBAR
               DIV(I,J,K) = (R(I)*U(I,J,K)-R(I-1)*U(I-1,J,K))*RDX(I)*RRN(I) + (V(I,J,K)-V(I,J-1,K))*RDY(J) + &
                            (W(I,J,K)-W(I,J,K-1))*RDZ(K)
            ENDDO
         ENDDO
      ENDDO
      DDDT = (DP-DIV)*RDT
   ELSEIF (CORRECTOR) THEN
      DO K = 1,KBAR
         DO J = 1,JBAR
            DO I = 1,IBAR
               DIV(I,J,K) = (R(I)*U(I,J,K) -R(I-1)*U(I-1,J,K)) *RDX(I)*RRN(I) + (V(I,J,K)- V(I,J-1,K)) *RDY(J) + &
                            (W(I,J,K) -W(I,J,K-1)) *RDZ(K) &
                          + (R(I)*US(I,J,K)-R(I-1)*US(I-1,J,K))*RDX(I)*RRN(I) + (VS(I,J,K)-VS(I,J-1,K))*RDY(J) + &
                            (WS(I,J,K)-WS(I,J,K-1))*RDZ(K)
            ENDDO
         ENDDO
      ENDDO
      D = DDDT
      DDDT = (2._EB*DP-DIV)*RDT
   ENDIF

ELSE TRUE_PROJECTION

   IF (PREDICTOR) THEN
      DDDT = (DS-D)*RDT
   ELSE
      D_NEW => WORK1
      D_NEW = DP
      DDDT  = (2._EB*D_NEW-DS-D)*RDT
      D     = D_NEW
   ENDIF

   ! Adjust dD/dt to correct error in divergence due to velocity matching at interpolated boundaries

   NO_SCARC_IF: IF (PRES_METHOD /= 'SCARC' .AND. PRES_METHOD /= 'USCARC') THEN
      DO IW=1,N_EXTERNAL_WALL_CELLS
         IF (EXTERNAL_WALL(IW)%NOM==0) CYCLE
         IIG = WALL(IW)%ONE_D%IIG
         JJG = WALL(IW)%ONE_D%JJG
         KKG = WALL(IW)%ONE_D%KKG
         IF (PREDICTOR) DDDT(IIG,JJG,KKG) = DDDT(IIG,JJG,KKG) + DS_CORR(IW)*RDT
         IF (CORRECTOR) DDDT(IIG,JJG,KKG) = DDDT(IIG,JJG,KKG) + (2._EB*D_CORR(IW)-DS_CORR(IW))*RDT
      ENDDO
   ENDIF NO_SCARC_IF

ENDIF TRUE_PROJECTION

T_USED(2)=T_USED(2)+CURRENT_TIME()-TNOW
END SUBROUTINE DIVERGENCE_PART_2


SUBROUTINE CHECK_DIVERGENCE(NM)

! Computes maximum velocity divergence

USE COMP_FUNCTIONS, ONLY: CURRENT_TIME
INTEGER, INTENT(IN) :: NM
INTEGER  :: I,J,K
REAL(EB) :: DIV,RES,TNOW
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW,DP

TNOW=CURRENT_TIME()
CALL POINT_TO_MESH(NM)

IF (PREDICTOR) THEN
   UU=>US
   VV=>VS
   WW=>WS
   DP=>DS
ELSEIF (CORRECTOR) THEN
   UU=>U
   VV=>V
   WW=>W
   DP=>D
ENDIF

RESMAX = 0._EB
DIVMX  = -10000._EB
DIVMN  =  10000._EB
IMX    = 0
JMX    = 0
KMX    = 0

DO K=1,KBAR
   DO J=1,JBAR
      LOOP1: DO I=1,IBAR
         IF (SOLID(CELL_INDEX(I,J,K))) CYCLE LOOP1
         SELECT CASE(CYLINDRICAL)
            CASE(.FALSE.)
               DIV = (UU(I,J,K)-UU(I-1,J,K))*RDX(I) + &
                     (VV(I,J,K)-VV(I,J-1,K))*RDY(J) + &
                     (WW(I,J,K)-WW(I,J,K-1))*RDZ(K)
            CASE(.TRUE.)
               DIV = (R(I)*UU(I,J,K)-R(I-1)*UU(I-1,J,K))*RDX(I)*RRN(I) +  &
                     (WW(I,J,K)-WW(I,J,K-1))*RDZ(K)
         END SELECT
         RES = ABS(DIV-DP(I,J,K))
         IF (ABS(RES)>=RESMAX) THEN
            RESMAX = ABS(RES)
            IRM=I
            JRM=J
            KRM=K
         ENDIF
         RESMAX = MAX(RES,RESMAX)
         IF (DIV>=DIVMX) THEN
            DIVMX = DIV
            IMX=I
            JMX=J
            KMX=K
         ENDIF
         IF (DIV<DIVMN) THEN
            DIVMN = DIV
            IMN=I
            JMN=J
            KMN=K
         ENDIF
      ENDDO LOOP1
   ENDDO
ENDDO

T_USED(2)=T_USED(2)+CURRENT_TIME()-TNOW
END SUBROUTINE CHECK_DIVERGENCE


END MODULE DIVG
