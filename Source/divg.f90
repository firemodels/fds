MODULE DIVG

USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_POINTERS

IMPLICIT NONE (TYPE,EXTERNAL)
PRIVATE

PUBLIC DIVERGENCE_PART_1,DIVERGENCE_PART_2,CHECK_DIVERGENCE

REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW,RHOP
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZP

CONTAINS

!> \brief Compute contributions to the divergence term
!> \param T Time (s)
!> \param DT Time step (s)
!> \param NM Mesh number

SUBROUTINE DIVERGENCE_PART_1(T,DT,NM)

USE COMP_FUNCTIONS, ONLY: CURRENT_TIME
USE MATH_FUNCTIONS, ONLY: INTERPOLATE1D_UNIFORM
USE PHYSICAL_FUNCTIONS, ONLY: GET_CONDUCTIVITY,GET_SPECIFIC_HEAT,GET_SENSIBLE_ENTHALPY_Z,GET_SENSIBLE_ENTHALPY,&
                              GET_VISCOSITY,GET_MOLECULAR_WEIGHT
USE TURBULENCE, ONLY: TENSOR_DIFFUSIVITY_MODEL
USE MANUFACTURED_SOLUTIONS, ONLY: DIFF_MMS,UF_MMS,WF_MMS,VD2D_MMS_Z_SRC
USE COMPLEX_GEOMETRY, ONLY : CC_CGSC, CC_UNKZ, CC_SOLID, CC_CUTCFE
USE CC_SCALARS, ONLY : ADD_CUTCELL_PSUM,ADD_LINKEDCELL_PSUM,SET_EXIMDIFFLX_3D,SET_EXIMRHOHSLIM_3D,&
                       SET_EXIMRHOZZLIM_3D,CC_DIVERGENCE_PART_1,CC_VELOCITY_FLUX,CFACE_PREDICT_NORMAL_VELOCITY

INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: T,DT
REAL(EB), POINTER, DIMENSION(:,:,:) :: KDTDX,KDTDY,KDTDZ,DP,KP,CP, &
          RHO_D,H_RHO_D_DZDX,H_RHO_D_DZDY,H_RHO_D_DZDZ,RTRM, &
          U_DOT_DEL_RHO_H_S,U_DOT_DEL_RHO_Z,RHO_D_TURB,R_H_G
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: RHO_D_DZDX,RHO_D_DZDY,RHO_D_DZDZ
REAL(EB), POINTER, DIMENSION(:,:) :: PBAR_P
REAL(EB) :: DELKDELT,VC,VC1,DTDX,DTDY,DTDZ,TNOW,DZDX,DZDY,DZDZ,TMP_G,DIV_DIFF_HEAT_FLUX,H_S,XHAT,ZHAT,TT,Q_Z,&
            D_Z_TEMP,D_Z_N(0:I_MAX_TEMP),RHO_D_DZDN_GET(1:N_TRACKED_SPECIES),UN_P,TMP_F_GAS,R_PFCT,RHO_D_DZDN
INTEGER :: IW,N,I,J,K,IPZ,N_ZZ_MAX,ICF
REAL(EB), ALLOCATABLE, DIMENSION(:) :: ZZ_GET
TYPE(SPECIES_MIXTURE_TYPE), POINTER :: SM
TYPE(WALL_TYPE), POINTER :: WC
TYPE(EXTERNAL_WALL_TYPE), POINTER :: EWC
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC
TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1
TYPE(CFACE_TYPE), POINTER :: CFA

IF (SOLID_PHASE_ONLY) RETURN

! Start the clock and set the pointers

TNOW=CURRENT_TIME()
CALL POINT_TO_MESH(NM)

SELECT CASE(PREDICTOR)
   CASE(.TRUE.)
      DP => DS
      PBAR_P => PBAR_S
      RHOP => RHOS
      UU=>U
      VV=>V
      WW=>W
      ZZP=>ZZS
   CASE(.FALSE.)
      DP => D
      PBAR_P => PBAR
      RHOP => RHO
      UU=>US
      VV=>VS
      WW=>WS
      ZZP=>ZZ
END SELECT

R_PBAR = 1._EB/PBAR_P

RTRM => WORK1

! Zero out divergence to start

DP = 0._EB

! Determine if pressure ZONEs have merged

IF (N_ZONE>0 .AND. OBST_CREATED_OR_REMOVED) CALL MERGE_PRESSURE_ZONES

! Compute normal component of velocity at boundaries, U_NORMAL_S in the PREDICTOR step, U_NORMAL in the CORRECTOR.

!$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(IW)
WALL_LOOP3: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
   CALL PREDICT_NORMAL_VELOCITY(IW,T,DT)
ENDDO WALL_LOOP3
!$OMP END PARALLEL DO

IF (CC_IBM) THEN
   CALL CC_VELOCITY_FLUX(NM,DT,PREDICTOR,RHOP,CORRECT_GRAV=.FALSE.) ! Link F.
   CALL CFACE_PREDICT_NORMAL_VELOCITY(T,DT)
ENDIF

! Compute species-related finite difference terms

RHO_D_DZDX => SWORK1
RHO_D_DZDY => SWORK2
RHO_D_DZDZ => SWORK3

! Save the largest value of the material and thermal diffusion coefficients for use in Von Neumann stability constraint

IF (CHECK_VN) D_Z_MAX = 0._EB

! Add species diffusion terms to divergence expression and compute diffusion term for species equations

SPECIES_GT_1_IF: IF (N_TOTAL_SCALARS>1) THEN

   DEL_RHO_D_DEL_Z = 0._EB
   RHO_D => WORK4
   IF (SIM_MODE/=DNS_MODE) THEN
      IF (SIM_MODE==LES_MODE) THEN
         RHO_D_TURB => WORK9
         RHO_D_TURB = MAX(0._EB,MU-MU_DNS)*RSC_T
      ELSE
         RHO_D = MAX(0._EB,MU)*RSC_T
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

      IF (SIM_MODE==LES_MODE .AND. .NOT.TENSOR_DIFFUSIVITY) THEN
         SM=>SPECIES_MIXTURE(N)
         IF (SM%SC_T_USER>TWENTY_EPSILON_EB) THEN
            RHO_D = RHO_D + RHO_D_TURB*SC_T/SM%SC_T_USER
         ELSE
            RHO_D = RHO_D + RHO_D_TURB
         ENDIF
      ENDIF

      ! Manufactured solution

      IF (PERIODIC_TEST==7) RHO_D = DIFF_MMS

      ! Store max diffusivity for stability check

      IF (CHECK_VN) THEN
         DO K=0,KBP1
            DO J=0,JBP1
               DO I=0,IBP1
                  D_Z_MAX(I,J,K) = MAX(D_Z_MAX(I,J,K),RHO_D(I,J,K)/(RHOP(I,J,K)+TWENTY_EPSILON_EB))
               ENDDO
            ENDDO
         ENDDO
      ENDIF

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

      ! If tensor diffusivity, add turbulent scalar flux to rho*D del Z

      IF (TENSOR_DIFFUSIVITY) CALL TENSOR_DIFFUSIVITY_MODEL(NM,N)

      ! Store rho*D_n grad Z_n at OPEN boundaries, flux match at INTERPOLATED boundaries, zero out otherwise

      !$OMP PARALLEL DO PRIVATE(IW,WC,BC,B1,EWC) SCHEDULE(GUIDED)
      WALL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
         WC => WALL(IW)
         IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE WALL_LOOP
         BC => BOUNDARY_COORD(WC%BC_INDEX)
         IF (WC%THIN .AND. BC%IOR<0) CYCLE WALL_LOOP  ! Avoid OpenMP race condition by processing on one side of thin OBST
         BOUNDARY_TYPE_SELECT: SELECT CASE(WC%BOUNDARY_TYPE)
            CASE DEFAULT
               SELECT CASE(BC%IOR)
                  CASE( 1); RHO_D_DZDX(BC%IIG-1,BC%JJG  ,BC%KKG  ,N) = 0._EB
                  CASE(-1); RHO_D_DZDX(BC%IIG  ,BC%JJG  ,BC%KKG  ,N) = 0._EB
                  CASE( 2); RHO_D_DZDY(BC%IIG  ,BC%JJG-1,BC%KKG  ,N) = 0._EB
                  CASE(-2); RHO_D_DZDY(BC%IIG  ,BC%JJG  ,BC%KKG  ,N) = 0._EB
                  CASE( 3); RHO_D_DZDZ(BC%IIG  ,BC%JJG  ,BC%KKG-1,N) = 0._EB
                  CASE(-3); RHO_D_DZDZ(BC%IIG  ,BC%JJG  ,BC%KKG  ,N) = 0._EB
               END SELECT
            CASE(OPEN_BOUNDARY,INTERPOLATED_BOUNDARY)
               B1 => BOUNDARY_PROP1(WC%B1_INDEX)
               EWC => EXTERNAL_WALL(IW)
               IF (EWC%NIC>1) THEN
                  ! overwrite coarse mesh diffusive flux with fine mesh average (flux matched) computed in wall_bc
                  SELECT CASE(BC%IOR)
                     CASE( 1); RHO_D_DZDX(BC%IIG-1,BC%JJG  ,BC%KKG  ,N) =  B1%RHO_D_DZDN_F(N)
                     CASE(-1); RHO_D_DZDX(BC%IIG  ,BC%JJG  ,BC%KKG  ,N) = -B1%RHO_D_DZDN_F(N)
                     CASE( 2); RHO_D_DZDY(BC%IIG  ,BC%JJG-1,BC%KKG  ,N) =  B1%RHO_D_DZDN_F(N)
                     CASE(-2); RHO_D_DZDY(BC%IIG  ,BC%JJG  ,BC%KKG  ,N) = -B1%RHO_D_DZDN_F(N)
                     CASE( 3); RHO_D_DZDZ(BC%IIG  ,BC%JJG  ,BC%KKG-1,N) =  B1%RHO_D_DZDN_F(N)
                     CASE(-3); RHO_D_DZDZ(BC%IIG  ,BC%JJG  ,BC%KKG  ,N) = -B1%RHO_D_DZDN_F(N)
                  END SELECT
               ELSE
                  ! store computed flux for output
                  SELECT CASE(BC%IOR)
                     CASE( 1); B1%RHO_D_DZDN_F(N) =  RHO_D_DZDX(BC%IIG-1,BC%JJG  ,BC%KKG  ,N)
                     CASE(-1); B1%RHO_D_DZDN_F(N) = -RHO_D_DZDX(BC%IIG  ,BC%JJG  ,BC%KKG  ,N)
                     CASE( 2); B1%RHO_D_DZDN_F(N) =  RHO_D_DZDY(BC%IIG  ,BC%JJG-1,BC%KKG  ,N)
                     CASE(-2); B1%RHO_D_DZDN_F(N) = -RHO_D_DZDY(BC%IIG  ,BC%JJG  ,BC%KKG  ,N)
                     CASE( 3); B1%RHO_D_DZDN_F(N) =  RHO_D_DZDZ(BC%IIG  ,BC%JJG  ,BC%KKG-1,N)
                     CASE(-3); B1%RHO_D_DZDN_F(N) = -RHO_D_DZDZ(BC%IIG  ,BC%JJG  ,BC%KKG  ,N)
                  END SELECT
               ENDIF
         END SELECT BOUNDARY_TYPE_SELECT
      ENDDO WALL_LOOP
      !$OMP END PARALLEL DO

   ENDDO DIFFUSIVE_FLUX_LOOP

   ! Ensure RHO_D terms sum to zero over all species.  Gather error into largest mass fraction present.

   IF (SIM_MODE==DNS_MODE .OR. SIM_MODE==LES_MODE .OR. TENSOR_DIFFUSIVITY) THEN
      ! for VLES and SVLES modes, the diffusivity is the same for all species
      ! so, as long as ZZP is realizable, the sum of diffusive fluxes will be zero
      ! and the flux corrections below are not required

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

      !$OMP PARALLEL

      !$OMP DO PRIVATE(TMP_G,H_S) SCHEDULE(GUIDED)
      DO K=0,KBAR
         DO J=0,JBAR
            DO I=0,IBAR
               TMP_G = 0.5_EB*(TMP(I+1,J,K)+TMP(I,J,K))
               CALL GET_SENSIBLE_ENTHALPY_Z(N,TMP_G,H_S)
               H_RHO_D_DZDX(I,J,K) = H_S*RHO_D_DZDX(I,J,K,N)
               TMP_G = 0.5_EB*(TMP(I,J+1,K)+TMP(I,J,K))
               CALL GET_SENSIBLE_ENTHALPY_Z(N,TMP_G,H_S)
               H_RHO_D_DZDY(I,J,K) = H_S*RHO_D_DZDY(I,J,K,N)
               TMP_G = 0.5_EB*(TMP(I,J,K+1)+TMP(I,J,K))
               CALL GET_SENSIBLE_ENTHALPY_Z(N,TMP_G,H_S)
               H_RHO_D_DZDZ(I,J,K) = H_S*RHO_D_DZDZ(I,J,K,N)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO

      ! Correct rho*D_n grad Z_n and h_n*rho*D_n grad Z_n at boundaries

      !$OMP DO PRIVATE(IW,WC,BC,B1,N_ZZ_MAX,RHO_D_DZDN,RHO_D_DZDN_GET,UN_P,TMP_F_GAS,H_S) SCHEDULE(GUIDED)
      WALL_LOOP_2: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
         WC => WALL(IW)
         IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY .OR. &
             WC%BOUNDARY_TYPE==OPEN_BOUNDARY .OR. &
             WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY) CYCLE WALL_LOOP_2
         BC => BOUNDARY_COORD(WC%BC_INDEX)
         B1 => BOUNDARY_PROP1(WC%B1_INDEX)

         N_ZZ_MAX = MAXLOC(B1%ZZ_F(1:N_TRACKED_SPECIES),1)
         RHO_D_DZDN = 2._EB*B1%RHO_D_F(N)*(ZZP(BC%IIG,BC%JJG,BC%KKG,N)-B1%ZZ_F(N))*B1%RDN
         IF (N==N_ZZ_MAX) THEN
            RHO_D_DZDN_GET = 2._EB*B1%RHO_D_F(:)*(ZZP(BC%IIG,BC%JJG,BC%KKG,:)-B1%ZZ_F(:))*B1%RDN
            RHO_D_DZDN = -(SUM(RHO_D_DZDN_GET(:))-RHO_D_DZDN)
         ENDIF
         B1%RHO_D_DZDN_F(N) = RHO_D_DZDN

         IF (WC%THIN .AND. BC%IOR<0) CYCLE WALL_LOOP_2  ! Avoid OpenMP race condition by processing only one side of thin OBST

         IF (STORE_SPECIES_FLUX) THEN
            IF (CORRECTOR) THEN
               SELECT CASE(BC%IOR)
                  CASE(-1) ; DIF_FXS(BC%IIG  ,BC%JJG  ,BC%KKG  ,N) =  RHO_D_DZDN
                  CASE( 1) ; DIF_FXS(BC%IIG-1,BC%JJG  ,BC%KKG  ,N) = -RHO_D_DZDN
                  CASE(-2) ; DIF_FYS(BC%IIG  ,BC%JJG  ,BC%KKG  ,N) =  RHO_D_DZDN
                  CASE( 2) ; DIF_FYS(BC%IIG  ,BC%JJG-1,BC%KKG  ,N) = -RHO_D_DZDN
                  CASE(-3) ; DIF_FZS(BC%IIG  ,BC%JJG  ,BC%KKG  ,N) =  RHO_D_DZDN
                  CASE( 3) ; DIF_FZS(BC%IIG  ,BC%JJG  ,BC%KKG-1,N) = -RHO_D_DZDN
               END SELECT
            ELSE
               SELECT CASE(BC%IOR)
                  CASE(-1) ; DIF_FX(BC%IIG  ,BC%JJG  ,BC%KKG  ,N) = 0.5_EB*(DIF_FXS(BC%IIG  ,BC%JJG  ,BC%KKG  ,N)+RHO_D_DZDN)
                  CASE( 1) ; DIF_FX(BC%IIG-1,BC%JJG  ,BC%KKG  ,N) = 0.5_EB*(DIF_FXS(BC%IIG-1,BC%JJG  ,BC%KKG  ,N)-RHO_D_DZDN)
                  CASE(-2) ; DIF_FY(BC%IIG  ,BC%JJG  ,BC%KKG  ,N) = 0.5_EB*(DIF_FYS(BC%IIG  ,BC%JJG  ,BC%KKG  ,N)+RHO_D_DZDN)
                  CASE( 2) ; DIF_FY(BC%IIG  ,BC%JJG-1,BC%KKG  ,N) = 0.5_EB*(DIF_FYS(BC%IIG  ,BC%JJG-1,BC%KKG  ,N)-RHO_D_DZDN)
                  CASE(-3) ; DIF_FZ(BC%IIG  ,BC%JJG  ,BC%KKG  ,N) = 0.5_EB*(DIF_FZS(BC%IIG  ,BC%JJG  ,BC%KKG  ,N)+RHO_D_DZDN)
                  CASE( 3) ; DIF_FZ(BC%IIG  ,BC%JJG  ,BC%KKG-1,N) = 0.5_EB*(DIF_FZS(BC%IIG  ,BC%JJG  ,BC%KKG-1,N)-RHO_D_DZDN)
               END SELECT
            ENDIF
         ENDIF

         IF (PREDICTOR) THEN
            UN_P = B1%U_NORMAL_S
         ELSE
            UN_P = B1%U_NORMAL
         ENDIF
         IF (WC%BOUNDARY_TYPE==SOLID_BOUNDARY .AND. UN_P>0._EB) THEN
            TMP_F_GAS = TMP(BC%IIG,BC%JJG,BC%KKG)
         ELSE
            TMP_F_GAS = B1%TMP_F
         ENDIF

         CALL GET_SENSIBLE_ENTHALPY_Z(N,TMP_F_GAS,H_S)

         SELECT CASE(BC%IOR)
            CASE( 1) ; RHO_D_DZDX(BC%IIG-1,BC%JJG  ,BC%KKG  ,N) =  RHO_D_DZDN
                     H_RHO_D_DZDX(BC%IIG-1,BC%JJG  ,BC%KKG    ) =  RHO_D_DZDN*H_S
            CASE(-1) ; RHO_D_DZDX(BC%IIG  ,BC%JJG  ,BC%KKG  ,N) = -RHO_D_DZDN
                     H_RHO_D_DZDX(BC%IIG  ,BC%JJG  ,BC%KKG    ) = -RHO_D_DZDN*H_S
            CASE( 2) ; RHO_D_DZDY(BC%IIG  ,BC%JJG-1,BC%KKG  ,N) =  RHO_D_DZDN
                     H_RHO_D_DZDY(BC%IIG  ,BC%JJG-1,BC%KKG    ) =  RHO_D_DZDN*H_S
            CASE(-2) ; RHO_D_DZDY(BC%IIG  ,BC%JJG  ,BC%KKG  ,N) = -RHO_D_DZDN
                     H_RHO_D_DZDY(BC%IIG  ,BC%JJG  ,BC%KKG    ) = -RHO_D_DZDN*H_S
            CASE( 3) ; RHO_D_DZDZ(BC%IIG  ,BC%JJG  ,BC%KKG-1,N) =  RHO_D_DZDN
                     H_RHO_D_DZDZ(BC%IIG  ,BC%JJG  ,BC%KKG-1  ) =  RHO_D_DZDN*H_S
            CASE(-3) ; RHO_D_DZDZ(BC%IIG  ,BC%JJG  ,BC%KKG  ,N) = -RHO_D_DZDN
                     H_RHO_D_DZDZ(BC%IIG  ,BC%JJG  ,BC%KKG    ) = -RHO_D_DZDN*H_S
         END SELECT

      ENDDO WALL_LOOP_2
      !$OMP END DO

      !$OMP DO PRIVATE(DIV_DIFF_HEAT_FLUX) SCHEDULE(STATIC)
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
      !$OMP END DO

      ! Compute div rho*D grad Z_n

      !$OMP DO SCHEDULE(STATIC)
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               DEL_RHO_D_DEL_Z(I,J,K,N) = (R(I)*RHO_D_DZDX(I,J,K,N)-R(I-1)*RHO_D_DZDX(I-1,J,K,N))*RDX(I)*RRN(I) + &
                                          (     RHO_D_DZDY(I,J,K,N)-       RHO_D_DZDY(I,J-1,K,N))*RDY(J)        + &
                                          (     RHO_D_DZDZ(I,J,K,N)-       RHO_D_DZDZ(I,J,K-1,N))*RDZ(K)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO

      !$OMP END PARALLEL

   ENDDO SPECIES_LOOP

ENDIF SPECIES_GT_1_IF

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
            IF (CELL(CELL_INDEX(I,J,K))%SOLID) CYCLE
            ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I,J,K,1:N_TRACKED_SPECIES)
            CALL GET_CONDUCTIVITY(ZZ_GET,KP(I,J,K),TMP(I,J,K))
         ENDDO
      ENDDO
   ENDDO
   DEALLOCATE(ZZ_GET)

   IF (SIM_MODE==LES_MODE .AND. .NOT.TENSOR_DIFFUSIVITY) THEN
      IF(.NOT.CONSTANT_SPECIFIC_HEAT_RATIO) THEN
         KP = KP + MAX(0._EB,(MU-MU_DNS))*CP*RPR_T
      ELSE
         KP = KP + MAX(0._EB,(MU-MU_DNS))*CPOPR
      ENDIF
   ENDIF

   BOUNDARY_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS
      WC=>WALL(IW)
      BC => BOUNDARY_COORD(WC%BC_INDEX)
      KP(BC%II,BC%JJ,BC%KK) = KP(BC%IIG,BC%JJG,BC%KKG)
   ENDDO BOUNDARY_LOOP

ELSE K_DNS_OR_LES

   ! normal VLES mode
   KP = MU*CPOPR

ENDIF K_DNS_OR_LES

! Store max diffusivity for stability check

IF (CHECK_VN .AND. .NOT.CONSTANT_SPECIFIC_HEAT_RATIO) THEN
   !$OMP PARALLEL DO SCHEDULE(STATIC)
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
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

! If tensor diffusivity, add turbulent thermal scalar flux to k*dT/dx, etc

IF (TENSOR_DIFFUSIVITY) CALL TENSOR_DIFFUSIVITY_MODEL(NM)

! Correct thermal gradient (k dT/dn) at boundaries

CORRECTION_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
   WC => WALL(IW)
   IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY .OR. WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY) CYCLE CORRECTION_LOOP
   BC => BOUNDARY_COORD(WC%BC_INDEX)
   B1 => BOUNDARY_PROP1(WC%B1_INDEX)
   IF (WC%BOUNDARY_TYPE==OPEN_BOUNDARY) THEN
      B1%K_G = 0.5_EB*(KP(BC%IIG,BC%JJG,BC%KKG)+KP(BC%II,BC%JJ,BC%KK))
      CYCLE CORRECTION_LOOP
   ELSE
      B1%K_G = KP(BC%IIG,BC%JJG,BC%KKG)
   ENDIF
   ! Q_LEAK accounts for enthalpy moving through leakage paths
   DP(BC%IIG,BC%JJG,BC%KKG) = DP(BC%IIG,BC%JJG,BC%KKG) - ( B1%AREA_ADJUST*B1%Q_CON_F*B1%RDN - B1%Q_LEAK )
   IF (WC%THIN .AND. BC%IOR<0) CYCLE CORRECTION_LOOP  ! Avoid OpenMP race condition by processing on one side of thin OBST
   SELECT CASE(BC%IOR)
      CASE( 1) ; KDTDX(BC%II  ,BC%JJ  ,BC%KK  ) = 0._EB
      CASE(-1) ; KDTDX(BC%II-1,BC%JJ  ,BC%KK  ) = 0._EB
      CASE( 2) ; KDTDY(BC%II  ,BC%JJ  ,BC%KK  ) = 0._EB
      CASE(-2) ; KDTDY(BC%II  ,BC%JJ-1,BC%KK  ) = 0._EB
      CASE( 3) ; KDTDZ(BC%II  ,BC%JJ  ,BC%KK  ) = 0._EB
      CASE(-3) ; KDTDZ(BC%II  ,BC%JJ  ,BC%KK-1) = 0._EB
   END SELECT
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

! Compute U_DOT_DEL_RHO_H_S and add to other enthalpy equation source terms

CONST_GAMMA_IF_1: IF (.NOT.CONSTANT_SPECIFIC_HEAT_RATIO) THEN

   CALL ENTHALPY_ADVECTION_NEW(U_DOT_DEL_RHO_H_S) ! Compute u dot grad rho h_s

   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            DP(I,J,K) = DP(I,J,K) - U_DOT_DEL_RHO_H_S(I,J,K)
         ENDDO
      ENDDO
   ENDDO

   IF (CC_IBM) CALL SET_EXIMRHOHSLIM_3D(NM) ! WORK2,WORK3,WORK4: Get flux limited \bar{rho Hs} on EXIM faces.

ENDIF CONST_GAMMA_IF_1

! Compute RTRM = 1/(rho*c_p*T) and multiply it by divergence terms already summed up

IF (CONSTANT_SPECIFIC_HEAT_RATIO) THEN

   !$OMP PARALLEL DO PRIVATE(IPZ) SCHEDULE(STATIC)
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
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
            RTRM(I,J,K) = R_H_G(I,J,K)/RHOP(I,J,K)
            DP(I,J,K) = RTRM(I,J,K)*DP(I,J,K)
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO

ENDIF

! Compute (1/rho) * Sum( (Wbar/W_alpha-h_s,alpha/cp*T) (del dot rho*D del Z_n - u dot del rho*Z_n)

CONST_GAMMA_IF_2: IF (.NOT.CONSTANT_SPECIFIC_HEAT_RATIO) THEN

   CALL SPECIES_ADVECTION_PART_1_NEW ! Compute and store face values of (rho Z_n)

   DO N=1,N_TRACKED_SPECIES

      CALL SPECIES_ADVECTION_PART_2(N,U_DOT_DEL_RHO_Z) ! Compute u dot grad rho Z_n

      SM  => SPECIES_MIXTURE(N)
      !$OMP PARALLEL DO PRIVATE(H_S) SCHEDULE(guided)
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (CELL(CELL_INDEX(I,J,K))%SOLID) CYCLE
               CALL GET_SENSIBLE_ENTHALPY_Z(N,TMP(I,J,K),H_S)
               DP(I,J,K) = DP(I,J,K) + (SM%RCON/RSUM(I,J,K) - H_S*R_H_G(I,J,K))* &
                    ( DEL_RHO_D_DEL_Z(I,J,K,N) - U_DOT_DEL_RHO_Z(I,J,K) )/RHOP(I,J,K)
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
   !$OMP PARALLEL DO SCHEDULE(STATIC)
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            DP(I,J,K) = DP(I,J,K) + RTRM(I,J,K)*0.5_EB*(WW(I,J,K)+WW(I,J,K-1))*RHO_0(K)*GVEC(3)
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO
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

IF (CC_IBM) THEN
   T_USED(2)=T_USED(2)+CURRENT_TIME()-TNOW
   CALL CC_DIVERGENCE_PART_1(T,DT,NM)
   TNOW=CURRENT_TIME()
ENDIF

! Calculate pressure rise in each of the pressure zones by summing divergence expression over each zone

IF_PRESSURE_ZONES: IF (N_ZONE>0) THEN

   R_PFCT = 1._EB
   DO K=1,KBAR
      DO J=1,JBAR
         VC1 = DY(J)*DZ(K)
         DO I=1,IBAR
            IF (INTERPOLATED_MESH(I,J,K)>0) CYCLE
            IPZ = PRESSURE_ZONE(I,J,K)
            IF (IPZ<1) CYCLE
            IF (CELL(CELL_INDEX(I,J,K))%SOLID) CYCLE
            VC = DX(I)*RC(I)*VC1
            DSUM(IPZ) = DSUM(IPZ) + VC*DP(I,J,K)
            IF (CC_IBM) THEN
               R_PFCT = 1._EB
               IF (CCVAR(I,J,K,CC_CGSC) == CC_SOLID) THEN
                  CYCLE
               ELSEIF(CCVAR(I,J,K,CC_CGSC) == CC_CUTCFE) THEN
                  CALL ADD_CUTCELL_PSUM(I,J,K,PBAR_P(K,IPZ),PSUM(IPZ)); CYCLE
               ELSEIF(CCVAR(I,J,K,CC_UNKZ) > 0) THEN
                  CALL ADD_LINKEDCELL_PSUM(I,J,K,VC,PBAR_P(K,IPZ),RTRM(I,J,K),PSUM(IPZ)); CYCLE
               ENDIF
            ENDIF
            PSUM(IPZ) = PSUM(IPZ) + VC*(R_PBAR(K,IPZ)*R_PFCT-RTRM(I,J,K))
         ENDDO
      ENDDO
   ENDDO

   ! Calculate the volume flux to the boundary of the pressure zone (int u dot dA)

   WALL_LOOP4: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
      WC => WALL(IW)
      BC => BOUNDARY_COORD(WC%BC_INDEX)
      IF (INTERPOLATED_MESH(BC%IIG,BC%JJG,BC%KKG)>0) CYCLE
      B1 => BOUNDARY_PROP1(WC%B1_INDEX)
      IPZ = B1%PRESSURE_ZONE
      IF (IPZ<1) CYCLE WALL_LOOP4
      IF (WC%BOUNDARY_TYPE/=SOLID_BOUNDARY) CYCLE WALL_LOOP4
      IF (PREDICTOR) USUM(IPZ) = USUM(IPZ) + B1%U_NORMAL_S*B1%AREA
      IF (CORRECTOR) USUM(IPZ) = USUM(IPZ) + B1%U_NORMAL  *B1%AREA
   ENDDO WALL_LOOP4


   CFACE_LOOP: DO ICF=INTERNAL_CFACE_CELLS_LB+1,INTERNAL_CFACE_CELLS_LB+N_INTERNAL_CFACE_CELLS
      CFA => CFACE(ICF)
      BC => BOUNDARY_COORD(CFA%BC_INDEX)
      B1 => BOUNDARY_PROP1(CFA%B1_INDEX)
      IPZ = B1%PRESSURE_ZONE
      IF (IPZ<1) CYCLE CFACE_LOOP
      IF (PREDICTOR) USUM(IPZ) = USUM(IPZ) + B1%U_NORMAL_S*B1%AREA
      IF (CORRECTOR) USUM(IPZ) = USUM(IPZ) + B1%U_NORMAL  *B1%AREA
   ENDDO CFACE_LOOP


ENDIF IF_PRESSURE_ZONES

T_USED(2)=T_USED(2)+CURRENT_TIME()-TNOW

END SUBROUTINE DIVERGENCE_PART_1


SUBROUTINE ENTHALPY_ADVECTION_NEW(U_DOT_DEL_RHO_H_S)

USE PHYSICAL_FUNCTIONS, ONLY: GET_SENSIBLE_ENTHALPY
USE MATH_FUNCTIONS, ONLY: GET_SCALAR_FACE_VALUE
REAL(EB), POINTER, DIMENSION(:,:,:) :: FX_H_S,FY_H_S,FZ_H_S,RHO_H_S_P,U_DOT_DEL_RHO_H_S
REAL(EB) :: UN,UN_P,TMP_F_GAS,DU_P,DU_M,DV_P,DV_M,DW_P,DW_M,DU,H_S
REAL(EB), ALLOCATABLE, DIMENSION(:) :: ZZ_GET
REAL(EB), DIMENSION(0:3,0:3,0:3) :: F_TEMP
REAL(EB), DIMENSION(-1:3,-1:3,-1:3) :: U_TEMP,Z_TEMP
INTEGER :: IC,I,J,K,IW
TYPE(WALL_TYPE), POINTER :: WC
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC
TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1

RHO_H_S_P=>WORK_PAD ; RHO_H_S_P = 0._EB
FX_H_S=>WORK2    ; FX_H_S = 0._EB
FY_H_S=>WORK3    ; FY_H_S = 0._EB
FZ_H_S=>WORK4    ; FZ_H_S = 0._EB
U_DOT_DEL_RHO_H_S=>WORK6 ; U_DOT_DEL_RHO_H_S=0._EB

! Compute and store rho*h_s

!$OMP PARALLEL PRIVATE(ZZ_GET, H_S)
ALLOCATE(ZZ_GET(1:N_TRACKED_SPECIES))
!$OMP DO SCHEDULE(static)
DO K=-1,KBP1+1
   DO J=-1,JBP1+1
      DO I=-1,IBP1+1
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

CALL GET_SCALAR_FACE_VALUE(UU,RHO_H_S_P,FX_H_S,0,IBAR,1,JBAR,1,KBAR,1,I_FLUX_LIMITER)
CALL GET_SCALAR_FACE_VALUE(VV,RHO_H_S_P,FY_H_S,1,IBAR,0,JBAR,1,KBAR,2,I_FLUX_LIMITER)
CALL GET_SCALAR_FACE_VALUE(WW,RHO_H_S_P,FZ_H_S,1,IBAR,1,JBAR,0,KBAR,3,I_FLUX_LIMITER)

ALLOCATE(ZZ_GET(1:N_TRACKED_SPECIES))

WALL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS

   WC=>WALL(IW)
   IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE WALL_LOOP
   BC => BOUNDARY_COORD(WC%BC_INDEX)
   B1 => BOUNDARY_PROP1(WC%B1_INDEX)

   ! Calculate the sensible enthalpy at the boundary. If the boundary is solid
   ! and the gas is flowing out, use the gas temperature for the calculation.

   IF (PREDICTOR) THEN
      UN_P = B1%U_NORMAL_S
   ELSE
      UN_P = B1%U_NORMAL
   ENDIF
   IF (WC%BOUNDARY_TYPE==SOLID_BOUNDARY .AND. UN_P>0._EB) THEN
      TMP_F_GAS = TMP(BC%IIG,BC%JJG,BC%KKG)
   ELSE
      TMP_F_GAS = B1%TMP_F
   ENDIF

   ZZ_GET(1:N_TRACKED_SPECIES) = B1%ZZ_F(1:N_TRACKED_SPECIES)
   CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,TMP_F_GAS)

   ! overwrite first off-wall advective flux if flow is away from the wall and if the face is not also a wall cell

   IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY .AND. WC%BOUNDARY_TYPE/=OPEN_BOUNDARY) THEN

      OFF_WALL_SELECT_1: SELECT CASE(BC%IOR)
         CASE( 1) OFF_WALL_SELECT_1
            !      ghost          FX/UU(II+1)
            ! ///   II   ///  II+1  |  II+2  | ...
            !                       ^ WALL_INDEX(II+1,+1)
            IF ((UU(BC%II+1,BC%JJ,BC%KK)>0._EB) .AND. .NOT.(CELL(CELL_INDEX(BC%II+1,BC%JJ,BC%KK))%WALL_INDEX(+1)>0)) THEN
               Z_TEMP(0:2,1,1) = (/RHO_H_S_P(BC%II+1,BC%JJ,BC%KK),RHO_H_S_P(BC%II+1:BC%II+2,BC%JJ,BC%KK)/)
               U_TEMP(1,1,1) = UU(BC%II+1,BC%JJ,BC%KK)
               CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,1,I_FLUX_LIMITER)
               FX_H_S(BC%II+1,BC%JJ,BC%KK) = F_TEMP(1,1,1)
            ENDIF
         CASE(-1) OFF_WALL_SELECT_1
            !            FX/UU(II-2)     ghost
            ! ... |  II-2  |  II-1  ///   II   ///
            !              ^ WALL_INDEX(II-1,-1)
            IF ((UU(BC%II-2,BC%JJ,BC%KK)<0._EB) .AND. .NOT.(CELL(CELL_INDEX(BC%II-1,BC%JJ,BC%KK))%WALL_INDEX(-1)>0)) THEN
               Z_TEMP(1:3,1,1) = (/RHO_H_S_P(BC%II-2:BC%II-1,BC%JJ,BC%KK),RHO_H_S_P(BC%II-1,BC%JJ,BC%KK)/)
               U_TEMP(1,1,1) = UU(BC%II-2,BC%JJ,BC%KK)
               CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,1,I_FLUX_LIMITER)
               FX_H_S(BC%II-2,BC%JJ,BC%KK) = F_TEMP(1,1,1)
            ENDIF
         CASE( 2) OFF_WALL_SELECT_1
            IF ((VV(BC%II,BC%JJ+1,BC%KK)>0._EB) .AND. .NOT.(CELL(CELL_INDEX(BC%II,BC%JJ+1,BC%KK))%WALL_INDEX(+2)>0)) THEN
               Z_TEMP(1,0:2,1) = (/RHO_H_S_P(BC%II,BC%JJ+1,BC%KK),RHO_H_S_P(BC%II,BC%JJ+1:BC%JJ+2,BC%KK)/)
               U_TEMP(1,1,1) = VV(BC%II,BC%JJ+1,BC%KK)
               CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,2,I_FLUX_LIMITER)
               FY_H_S(BC%II,BC%JJ+1,BC%KK) = F_TEMP(1,1,1)
            ENDIF
         CASE(-2) OFF_WALL_SELECT_1
            IF ((VV(BC%II,BC%JJ-2,BC%KK)<0._EB) .AND. .NOT.(CELL(CELL_INDEX(BC%II,BC%JJ-1,BC%KK))%WALL_INDEX(-2)>0)) THEN
               Z_TEMP(1,1:3,1) = (/RHO_H_S_P(BC%II,BC%JJ-2:BC%JJ-1,BC%KK),RHO_H_S_P(BC%II,BC%JJ-1,BC%KK)/)
               U_TEMP(1,1,1) = VV(BC%II,BC%JJ-2,BC%KK)
               CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,2,I_FLUX_LIMITER)
               FY_H_S(BC%II,BC%JJ-2,BC%KK) = F_TEMP(1,1,1)
            ENDIF
         CASE( 3) OFF_WALL_SELECT_1
            IF ((WW(BC%II,BC%JJ,BC%KK+1)>0._EB) .AND. .NOT.(CELL(CELL_INDEX(BC%II,BC%JJ,BC%KK+1))%WALL_INDEX(+3)>0)) THEN
               Z_TEMP(1,1,0:2) = (/RHO_H_S_P(BC%II,BC%JJ,BC%KK+1),RHO_H_S_P(BC%II,BC%JJ,BC%KK+1:BC%KK+2)/)
               U_TEMP(1,1,1) = WW(BC%II,BC%JJ,BC%KK+1)
               CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,3,I_FLUX_LIMITER)
               FZ_H_S(BC%II,BC%JJ,BC%KK+1) = F_TEMP(1,1,1)
            ENDIF
         CASE(-3) OFF_WALL_SELECT_1
            IF ((WW(BC%II,BC%JJ,BC%KK-2)<0._EB) .AND. .NOT.(CELL(CELL_INDEX(BC%II,BC%JJ,BC%KK-1))%WALL_INDEX(-3)>0)) THEN
               Z_TEMP(1,1,1:3) = (/RHO_H_S_P(BC%II,BC%JJ,BC%KK-2:BC%KK-1),RHO_H_S_P(BC%II,BC%JJ,BC%KK-1)/)
               U_TEMP(1,1,1) = WW(BC%II,BC%JJ,BC%KK-2)
               CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,3,I_FLUX_LIMITER)
               FZ_H_S(BC%II,BC%JJ,BC%KK-2) = F_TEMP(1,1,1)
            ENDIF
      END SELECT OFF_WALL_SELECT_1

   ENDIF

   SELECT CASE(WC%BOUNDARY_TYPE)
      CASE DEFAULT
         IOR_SELECT: SELECT CASE(BC%IOR)
            CASE( 1); UN = UU(BC%II  ,BC%JJ  ,BC%KK  )
            CASE(-1); UN = UU(BC%II-1,BC%JJ  ,BC%KK  )
            CASE( 2); UN = VV(BC%II  ,BC%JJ  ,BC%KK  )
            CASE(-2); UN = VV(BC%II  ,BC%JJ-1,BC%KK  )
            CASE( 3); UN = WW(BC%II  ,BC%JJ  ,BC%KK  )
            CASE(-3); UN = WW(BC%II  ,BC%JJ  ,BC%KK-1)
         END SELECT IOR_SELECT
      CASE(SOLID_BOUNDARY)
         IF (PREDICTOR) UN = -SIGN(1._EB,REAL(BC%IOR,EB))*B1%U_NORMAL_S
         IF (CORRECTOR) UN = -SIGN(1._EB,REAL(BC%IOR,EB))*B1%U_NORMAL
      CASE(INTERPOLATED_BOUNDARY)
         UN = UVW_SAVE(IW)
   END SELECT

   DU = (B1%RHO_F*H_S - RHO_H_S_P(BC%IIG,BC%JJG,BC%KKG))*UN
   U_DOT_DEL_RHO_H_S(BC%IIG,BC%JJG,BC%KKG) = U_DOT_DEL_RHO_H_S(BC%IIG,BC%JJG,BC%KKG) - SIGN(1._EB,REAL(BC%IOR,EB))*DU*B1%RDN

ENDDO WALL_LOOP

DEALLOCATE(ZZ_GET)

! FDS Tech Guide (B.12-B.14)

!$OMP PARALLEL DO PRIVATE(IC,DU_P,DU_M,DV_P,DV_M,DW_P,DW_M)
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IC = CELL_INDEX(I,J,K)
         IF (CELL(IC)%SOLID) CYCLE
         DU_P = 0._EB
         DU_M = 0._EB
         DV_P = 0._EB
         DV_M = 0._EB
         DW_P = 0._EB
         DW_M = 0._EB
         IF (CELL(IC)%WALL_INDEX( 1)==0) DU_P = (FX_H_S(I,J,K)   - RHO_H_S_P(I,J,K))*UU(I,J,K)
         IF (CELL(IC)%WALL_INDEX(-1)==0) DU_M = (FX_H_S(I-1,J,K) - RHO_H_S_P(I,J,K))*UU(I-1,J,K)
         IF (CELL(IC)%WALL_INDEX( 2)==0) DV_P = (FY_H_S(I,J,K)   - RHO_H_S_P(I,J,K))*VV(I,J,K)
         IF (CELL(IC)%WALL_INDEX(-2)==0) DV_M = (FY_H_S(I,J-1,K) - RHO_H_S_P(I,J,K))*VV(I,J-1,K)
         IF (CELL(IC)%WALL_INDEX( 3)==0) DW_P = (FZ_H_S(I,J,K)   - RHO_H_S_P(I,J,K))*WW(I,J,K)
         IF (CELL(IC)%WALL_INDEX(-3)==0) DW_M = (FZ_H_S(I,J,K-1) - RHO_H_S_P(I,J,K))*WW(I,J,K-1)
         U_DOT_DEL_RHO_H_S(I,J,K) = U_DOT_DEL_RHO_H_S(I,J,K) + (DU_P-DU_M)*RDX(I) + (DV_P-DV_M)*RDY(J) + (DW_P-DW_M)*RDZ(K)
      ENDDO
   ENDDO
ENDDO
!$OMP END PARALLEL DO

END SUBROUTINE ENTHALPY_ADVECTION_NEW


SUBROUTINE SPECIES_ADVECTION_PART_1_NEW

USE MATH_FUNCTIONS, ONLY: GET_SCALAR_FACE_VALUE
USE PHYSICAL_FUNCTIONS, ONLY: GET_MOLECULAR_WEIGHT
! INTEGER, INTENT(IN) :: N
REAL(EB), POINTER, DIMENSION(:,:,:) :: RHO_Z_P,RHO_RMW
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: FX_ZZ,FY_ZZ,FZ_ZZ
REAL(EB) :: ZZ_GET(1:N_TRACKED_SPECIES),MW_G
REAL(EB), DIMENSION(0:3,0:3,0:3) :: F_TEMP
REAL(EB), DIMENSION(-1:3,-1:3,-1:3) :: U_TEMP,Z_TEMP
INTEGER :: I,J,K,IW,N !,II,JJ,KK,IOR,IC,IIG,JJG,KKG
TYPE(WALL_TYPE), POINTER :: WC
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC
TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1

FX_ZZ=>SWORK1
FY_ZZ=>SWORK2
FZ_ZZ=>SWORK3
RHO_Z_P=>WORK_PAD

! Species face values

SPECIES_LOOP: DO N=1,N_TOTAL_SCALARS

   !$OMP PARALLEL DO
   DO K=-1,KBP1+1
      DO J=-1,JBP1+1
         DO I=-1,IBP1+1
            RHO_Z_P(I,J,K) = RHOP(I,J,K)*ZZP(I,J,K,N)
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO

   ! Compute scalar face values

   CALL GET_SCALAR_FACE_VALUE(UU,RHO_Z_P,FX_ZZ(:,:,:,N),0,IBAR,1,JBAR,1,KBAR,1,I_FLUX_LIMITER)
   CALL GET_SCALAR_FACE_VALUE(VV,RHO_Z_P,FY_ZZ(:,:,:,N),1,IBAR,0,JBAR,1,KBAR,2,I_FLUX_LIMITER)
   CALL GET_SCALAR_FACE_VALUE(WW,RHO_Z_P,FZ_ZZ(:,:,:,N),1,IBAR,1,JBAR,0,KBAR,3,I_FLUX_LIMITER)

   !$OMP PARALLEL DO PRIVATE(IW,WC,BC,B1,U_TEMP,Z_TEMP,F_TEMP)
   WALL_LOOP_2: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
      WC=>WALL(IW)
      IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE WALL_LOOP_2
      BC=>BOUNDARY_COORD(WC%BC_INDEX)
      B1=>BOUNDARY_PROP1(WC%B1_INDEX)

      ! Overwrite first off-wall advective flux if flow is away from the wall and if the face is not also a wall cell

      OFF_WALL_IF_2: IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY .AND. WC%BOUNDARY_TYPE/=OPEN_BOUNDARY) THEN

         OFF_WALL_SELECT_2: SELECT CASE(BC%IOR)
            CASE( 1) OFF_WALL_SELECT_2
               !      ghost          FX/UU(II+1)
               ! ///   II   ///  II+1  |  II+2  | ...
               !                       ^ WALL_INDEX(II+1,+1)
               IF ((UU(BC%II+1,BC%JJ,BC%KK)>0._EB) .AND. .NOT.(CELL(CELL_INDEX(BC%II+1,BC%JJ,BC%KK))%WALL_INDEX(+1)>0)) THEN
                  Z_TEMP(0:2,1,1) = (/RHO_Z_P(BC%II+1,BC%JJ,BC%KK),RHO_Z_P(BC%II+1:BC%II+2,BC%JJ,BC%KK)/)
                  U_TEMP(1,1,1) = UU(BC%II+1,BC%JJ,BC%KK)
                  CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,1,I_FLUX_LIMITER)
                  FX_ZZ(BC%II+1,BC%JJ,BC%KK,N) = F_TEMP(1,1,1)
               ENDIF
            CASE(-1) OFF_WALL_SELECT_2
               !            FX/UU(II-2)     ghost
               ! ... |  II-2  |  II-1  ///   II   ///
               !              ^ WALL_INDEX(II-1,-1)
               IF ((UU(BC%II-2,BC%JJ,BC%KK)<0._EB) .AND. .NOT.(CELL(CELL_INDEX(BC%II-1,BC%JJ,BC%KK))%WALL_INDEX(-1)>0)) THEN
                  Z_TEMP(1:3,1,1) = (/RHO_Z_P(BC%II-2:BC%II-1,BC%JJ,BC%KK),RHO_Z_P(BC%II-1,BC%JJ,BC%KK)/)
                  U_TEMP(1,1,1) = UU(BC%II-2,BC%JJ,BC%KK)
                  CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,1,I_FLUX_LIMITER)
                  FX_ZZ(BC%II-2,BC%JJ,BC%KK,N) = F_TEMP(1,1,1)
               ENDIF
            CASE( 2) OFF_WALL_SELECT_2
               IF ((VV(BC%II,BC%JJ+1,BC%KK)>0._EB) .AND. .NOT.(CELL(CELL_INDEX(BC%II,BC%JJ+1,BC%KK))%WALL_INDEX(+2)>0)) THEN
                  Z_TEMP(1,0:2,1) = (/RHO_Z_P(BC%II,BC%JJ+1,BC%KK),RHO_Z_P(BC%II,BC%JJ+1:BC%JJ+2,BC%KK)/)
                  U_TEMP(1,1,1) = VV(BC%II,BC%JJ+1,BC%KK)
                  CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,2,I_FLUX_LIMITER)
                  FY_ZZ(BC%II,BC%JJ+1,BC%KK,N) = F_TEMP(1,1,1)
               ENDIF
            CASE(-2) OFF_WALL_SELECT_2
               IF ((VV(BC%II,BC%JJ-2,BC%KK)<0._EB) .AND. .NOT.(CELL(CELL_INDEX(BC%II,BC%JJ-1,BC%KK))%WALL_INDEX(-2)>0)) THEN
                  Z_TEMP(1,1:3,1) = (/RHO_Z_P(BC%II,BC%JJ-2:BC%JJ-1,BC%KK),RHO_Z_P(BC%II,BC%JJ-1,BC%KK)/)
                  U_TEMP(1,1,1) = VV(BC%II,BC%JJ-2,BC%KK)
                  CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,2,I_FLUX_LIMITER)
                  FY_ZZ(BC%II,BC%JJ-2,BC%KK,N) = F_TEMP(1,1,1)
               ENDIF
            CASE( 3) OFF_WALL_SELECT_2
               IF ((WW(BC%II,BC%JJ,BC%KK+1)>0._EB) .AND. .NOT.(CELL(CELL_INDEX(BC%II,BC%JJ,BC%KK+1))%WALL_INDEX(+3)>0)) THEN
                  Z_TEMP(1,1,0:2) = (/RHO_Z_P(BC%II,BC%JJ,BC%KK+1),RHO_Z_P(BC%II,BC%JJ,BC%KK+1:BC%KK+2)/)
                  U_TEMP(1,1,1) = WW(BC%II,BC%JJ,BC%KK+1)
                  CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,3,I_FLUX_LIMITER)
                  FZ_ZZ(BC%II,BC%JJ,BC%KK+1,N) = F_TEMP(1,1,1)
               ENDIF
            CASE(-3) OFF_WALL_SELECT_2
               IF ((WW(BC%II,BC%JJ,BC%KK-2)<0._EB) .AND. .NOT.(CELL(CELL_INDEX(BC%II,BC%JJ,BC%KK-1))%WALL_INDEX(-3)>0)) THEN
                  Z_TEMP(1,1,1:3) = (/RHO_Z_P(BC%II,BC%JJ,BC%KK-2:BC%KK-1),RHO_Z_P(BC%II,BC%JJ,BC%KK-1)/)
                  U_TEMP(1,1,1) = WW(BC%II,BC%JJ,BC%KK-2)
                  CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,3,I_FLUX_LIMITER)
                  FZ_ZZ(BC%II,BC%JJ,BC%KK-2,N) = F_TEMP(1,1,1)
               ENDIF
         END SELECT OFF_WALL_SELECT_2

      ENDIF OFF_WALL_IF_2

   ENDDO WALL_LOOP_2
   !$OMP END PARALLEL DO

ENDDO SPECIES_LOOP

FACE_CORRECTION_IF: IF (FLUX_LIMITER_MW_CORRECTION) THEN

   ! Repeat the above for DENSITY

   RHO_RMW=>WORK_PAD

   !$OMP PARALLEL DO PRIVATE(ZZ_GET,MW_G) SCHEDULE(STATIC)
   DO K=-1,KBP1+1
      DO J=-1,JBP1+1
         DO I=-1,IBP1+1
            ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I,J,K,1:N_TRACKED_SPECIES)
            CALL GET_MOLECULAR_WEIGHT(ZZ_GET,MW_G)
            RHO_RMW(I,J,K) = RHOP(I,J,K)/MW_G
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO

   CALL GET_SCALAR_FACE_VALUE(UU,RHO_RMW,FX_ZZ(:,:,:,0),0,IBAR,1,JBAR,1,KBAR,1,I_FLUX_LIMITER)
   CALL GET_SCALAR_FACE_VALUE(VV,RHO_RMW,FY_ZZ(:,:,:,0),1,IBAR,0,JBAR,1,KBAR,2,I_FLUX_LIMITER)
   CALL GET_SCALAR_FACE_VALUE(WW,RHO_RMW,FZ_ZZ(:,:,:,0),1,IBAR,1,JBAR,0,KBAR,3,I_FLUX_LIMITER)

   !$OMP PARALLEL DO PRIVATE(IW,WC,BC,B1,U_TEMP,Z_TEMP,F_TEMP,ZZ_GET)
   WALL_LOOP_3: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
      WC=>WALL(IW)
      IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE WALL_LOOP_3
      BC=>BOUNDARY_COORD(WC%BC_INDEX)
      B1=>BOUNDARY_PROP1(WC%B1_INDEX)

      ! Overwrite first off-wall advective flux if flow is away from the wall and if the face is not also a wall cell

      OFF_WALL_IF_3: IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY .AND. WC%BOUNDARY_TYPE/=OPEN_BOUNDARY) THEN

         OFF_WALL_SELECT_3: SELECT CASE(BC%IOR)
            CASE( 1) OFF_WALL_SELECT_3
               !      ghost          FX/UU(II+1)
               ! ///   II   ///  II+1  |  II+2  | ...
               !                       ^ WALL_INDEX(II+1,+1)
               IF ((UU(BC%II+1,BC%JJ,BC%KK)>0._EB) .AND. .NOT.(CELL(CELL_INDEX(BC%II+1,BC%JJ,BC%KK))%WALL_INDEX(+1)>0)) THEN
                  Z_TEMP(0:2,1,1) = (/RHO_Z_P(BC%II+1,BC%JJ,BC%KK),RHO_Z_P(BC%II+1:BC%II+2,BC%JJ,BC%KK)/)
                  U_TEMP(1,1,1) = UU(BC%II+1,BC%JJ,BC%KK)
                  CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,1,I_FLUX_LIMITER)
                  FX_ZZ(BC%II+1,BC%JJ,BC%KK,0) = F_TEMP(1,1,1)
               ENDIF
            CASE(-1) OFF_WALL_SELECT_3
               !            FX/UU(II-2)     ghost
               ! ... |  II-2  |  II-1  ///   II   ///
               !              ^ WALL_INDEX(II-1,-1)
               IF ((UU(BC%II-2,BC%JJ,BC%KK)<0._EB) .AND. .NOT.(CELL(CELL_INDEX(BC%II-1,BC%JJ,BC%KK))%WALL_INDEX(-1)>0)) THEN
                  Z_TEMP(1:3,1,1) = (/RHO_Z_P(BC%II-2:BC%II-1,BC%JJ,BC%KK),RHO_Z_P(BC%II-1,BC%JJ,BC%KK)/)
                  U_TEMP(1,1,1) = UU(BC%II-2,BC%JJ,BC%KK)
                  CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,1,I_FLUX_LIMITER)
                  FX_ZZ(BC%II-2,BC%JJ,BC%KK,0) = F_TEMP(1,1,1)
               ENDIF
            CASE( 2) OFF_WALL_SELECT_3
               IF ((VV(BC%II,BC%JJ+1,BC%KK)>0._EB) .AND. .NOT.(CELL(CELL_INDEX(BC%II,BC%JJ+1,BC%KK))%WALL_INDEX(+2)>0)) THEN
                  Z_TEMP(1,0:2,1) = (/RHO_Z_P(BC%II,BC%JJ+1,BC%KK),RHO_Z_P(BC%II,BC%JJ+1:BC%JJ+2,BC%KK)/)
                  U_TEMP(1,1,1) = VV(BC%II,BC%JJ+1,BC%KK)
                  CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,2,I_FLUX_LIMITER)
                  FY_ZZ(BC%II,BC%JJ+1,BC%KK,0) = F_TEMP(1,1,1)
               ENDIF
            CASE(-2) OFF_WALL_SELECT_3
               IF ((VV(BC%II,BC%JJ-2,BC%KK)<0._EB) .AND. .NOT.(CELL(CELL_INDEX(BC%II,BC%JJ-1,BC%KK))%WALL_INDEX(-2)>0)) THEN
                  Z_TEMP(1,1:3,1) = (/RHO_Z_P(BC%II,BC%JJ-2:BC%JJ-1,BC%KK),RHO_Z_P(BC%II,BC%JJ-1,BC%KK)/)
                  U_TEMP(1,1,1) = VV(BC%II,BC%JJ-2,BC%KK)
                  CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,2,I_FLUX_LIMITER)
                  FY_ZZ(BC%II,BC%JJ-2,BC%KK,0) = F_TEMP(1,1,1)
               ENDIF
            CASE( 3) OFF_WALL_SELECT_3
               IF ((WW(BC%II,BC%JJ,BC%KK+1)>0._EB) .AND. .NOT.(CELL(CELL_INDEX(BC%II,BC%JJ,BC%KK+1))%WALL_INDEX(+3)>0)) THEN
                  Z_TEMP(1,1,0:2) = (/RHO_Z_P(BC%II,BC%JJ,BC%KK+1),RHO_Z_P(BC%II,BC%JJ,BC%KK+1:BC%KK+2)/)
                  U_TEMP(1,1,1) = WW(BC%II,BC%JJ,BC%KK+1)
                  CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,3,I_FLUX_LIMITER)
                  FZ_ZZ(BC%II,BC%JJ,BC%KK+1,0) = F_TEMP(1,1,1)
               ENDIF
            CASE(-3) OFF_WALL_SELECT_3
               IF ((WW(BC%II,BC%JJ,BC%KK-2)<0._EB) .AND. .NOT.(CELL(CELL_INDEX(BC%II,BC%JJ,BC%KK-1))%WALL_INDEX(-3)>0)) THEN
                  Z_TEMP(1,1,1:3) = (/RHO_Z_P(BC%II,BC%JJ,BC%KK-2:BC%KK-1),RHO_Z_P(BC%II,BC%JJ,BC%KK-1)/)
                  U_TEMP(1,1,1) = WW(BC%II,BC%JJ,BC%KK-2)
                  CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,3,I_FLUX_LIMITER)
                  FZ_ZZ(BC%II,BC%JJ,BC%KK-2,0) = F_TEMP(1,1,1)
               ENDIF
         END SELECT OFF_WALL_SELECT_3

      ENDIF OFF_WALL_IF_3

   ENDDO WALL_LOOP_3
   !$OMP END PARALLEL DO

   ! Now correct the max face value of (RHO*ZZ) such that SUM(RHO*ZZ/MW)_FACE = RHO_FACE/MW_FACE
   ! (necessary condition to preserve isothermal flow)

   !$OMP PARALLEL DO PRIVATE(N,MW_G) SCHEDULE(STATIC)
   DO K=0,KBAR
      DO J=0,JBAR
         DO I=0,IBAR
            N=MAXLOC(FX_ZZ(I,J,K,1:N_TRACKED_SPECIES),1)
            MW_G = SPECIES_MIXTURE(N)%MW
            FX_ZZ(I,J,K,N) = MW_G*MAX( 0._EB, FX_ZZ(I,J,K,0) &
                                      - SUM(FX_ZZ(I,J,K,1:(N-1))/SPECIES_MIXTURE(1:(N-1))%MW) &
                                      - SUM(FX_ZZ(I,J,K,(N+1):N_TRACKED_SPECIES)/SPECIES_MIXTURE((N+1):N_TRACKED_SPECIES)%MW) )

            N=MAXLOC(FY_ZZ(I,J,K,1:N_TRACKED_SPECIES),1)
            MW_G = SPECIES_MIXTURE(N)%MW
            FY_ZZ(I,J,K,N) = MW_G*MAX( 0._EB, FY_ZZ(I,J,K,0) &
                                      - SUM(FY_ZZ(I,J,K,1:(N-1))/SPECIES_MIXTURE(1:(N-1))%MW) &
                                      - SUM(FY_ZZ(I,J,K,(N+1):N_TRACKED_SPECIES)/SPECIES_MIXTURE((N+1):N_TRACKED_SPECIES)%MW) )

            N=MAXLOC(FZ_ZZ(I,J,K,1:N_TRACKED_SPECIES),1)
            MW_G = SPECIES_MIXTURE(N)%MW
            FZ_ZZ(I,J,K,N) = MW_G*MAX( 0._EB, FZ_ZZ(I,J,K,0) &
                                      - SUM(FZ_ZZ(I,J,K,1:(N-1))/SPECIES_MIXTURE(1:(N-1))%MW) &
                                      - SUM(FZ_ZZ(I,J,K,(N+1):N_TRACKED_SPECIES)/SPECIES_MIXTURE((N+1):N_TRACKED_SPECIES)%MW) )
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO

ENDIF FACE_CORRECTION_IF

END SUBROUTINE SPECIES_ADVECTION_PART_1_NEW


SUBROUTINE SPECIES_ADVECTION_PART_2(N,U_DOT_DEL_RHO_Z)

INTEGER, INTENT(IN) :: N
REAL(EB), POINTER, DIMENSION(:,:,:) :: U_DOT_DEL_RHO_Z
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: FX_ZZ,FY_ZZ,FZ_ZZ
REAL(EB) :: UN,DU_P,DU_M,DV_P,DV_M,DW_P,DW_M,DU
INTEGER :: IC,I,J,K,IW
TYPE(WALL_TYPE), POINTER :: WC
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC
TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1

FX_ZZ=>SWORK1
FY_ZZ=>SWORK2
FZ_ZZ=>SWORK3
U_DOT_DEL_RHO_Z=>WORK7 ; U_DOT_DEL_RHO_Z = 0._EB

WALL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS

   WC=>WALL(IW)
   IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE WALL_LOOP
   B1 => BOUNDARY_PROP1(WC%B1_INDEX)
   BC => BOUNDARY_COORD(WC%BC_INDEX)

   ! Correct flux terms at the boundary

   SELECT CASE(WC%BOUNDARY_TYPE)
      CASE DEFAULT
         IOR_SELECT: SELECT CASE(BC%IOR)
            CASE( 1); UN = UU(BC%II  ,BC%JJ  ,BC%KK  )
            CASE(-1); UN = UU(BC%II-1,BC%JJ  ,BC%KK  )
            CASE( 2); UN = VV(BC%II  ,BC%JJ  ,BC%KK  )
            CASE(-2); UN = VV(BC%II  ,BC%JJ-1,BC%KK  )
            CASE( 3); UN = WW(BC%II  ,BC%JJ  ,BC%KK  )
            CASE(-3); UN = WW(BC%II  ,BC%JJ  ,BC%KK-1)
         END SELECT IOR_SELECT
      CASE(SOLID_BOUNDARY)
         IF (PREDICTOR) UN = -SIGN(1._EB,REAL(BC%IOR,EB))*B1%U_NORMAL_S
         IF (CORRECTOR) UN = -SIGN(1._EB,REAL(BC%IOR,EB))*B1%U_NORMAL
      CASE(INTERPOLATED_BOUNDARY)
         UN = UVW_SAVE(IW)
   END SELECT

   DU = (B1%RHO_F*B1%ZZ_F(N) - RHOP(BC%IIG,BC%JJG,BC%KKG)*ZZP(BC%IIG,BC%JJG,BC%KKG,N))*UN
   U_DOT_DEL_RHO_Z(BC%IIG,BC%JJG,BC%KKG) = U_DOT_DEL_RHO_Z(BC%IIG,BC%JJG,BC%KKG) - SIGN(1._EB,REAL(BC%IOR,EB))*DU*B1%RDN

ENDDO WALL_LOOP

!$OMP PARALLEL DO PRIVATE(IC,DU_P,DU_M,DV_P,DV_M,DW_P,DW_M)
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IC = CELL_INDEX(I,J,K)
         IF (CELL(IC)%SOLID) CYCLE
         DU_P = 0._EB
         DU_M = 0._EB
         DV_P = 0._EB
         DV_M = 0._EB
         DW_P = 0._EB
         DW_M = 0._EB
         IF (CELL(IC)%WALL_INDEX( 1)==0) DU_P = (FX_ZZ(I,J,K,N)   - RHOP(I,J,K)*ZZP(I,J,K,N))*UU(I,J,K)
         IF (CELL(IC)%WALL_INDEX(-1)==0) DU_M = (FX_ZZ(I-1,J,K,N) - RHOP(I,J,K)*ZZP(I,J,K,N))*UU(I-1,J,K)
         IF (CELL(IC)%WALL_INDEX( 2)==0) DV_P = (FY_ZZ(I,J,K,N)   - RHOP(I,J,K)*ZZP(I,J,K,N))*VV(I,J,K)
         IF (CELL(IC)%WALL_INDEX(-2)==0) DV_M = (FY_ZZ(I,J-1,K,N) - RHOP(I,J,K)*ZZP(I,J,K,N))*VV(I,J-1,K)
         IF (CELL(IC)%WALL_INDEX( 3)==0) DW_P = (FZ_ZZ(I,J,K,N)   - RHOP(I,J,K)*ZZP(I,J,K,N))*WW(I,J,K)
         IF (CELL(IC)%WALL_INDEX(-3)==0) DW_M = (FZ_ZZ(I,J,K-1,N) - RHOP(I,J,K)*ZZP(I,J,K,N))*WW(I,J,K-1)
         U_DOT_DEL_RHO_Z(I,J,K) = U_DOT_DEL_RHO_Z(I,J,K) + (DU_P-DU_M)*RDX(I) + (DV_P-DV_M)*RDY(J) + (DW_P-DW_M)*RDZ(K)
      ENDDO
   ENDDO
ENDDO
!$OMP END PARALLEL DO

END SUBROUTINE SPECIES_ADVECTION_PART_2


SUBROUTINE MERGE_PRESSURE_ZONES

INTEGER :: IW,IPZ,IOPZ
TYPE(WALL_TYPE), POINTER :: WC
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC

DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
   WC=>WALL(IW)
   IF (WC%BOUNDARY_TYPE/=NULL_BOUNDARY .AND. WC%BOUNDARY_TYPE/=OPEN_BOUNDARY .AND. WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY) CYCLE
   BC => BOUNDARY_COORD(WC%BC_INDEX)
   IF (CELL(CELL_INDEX(BC%IIG,BC%JJG,BC%KKG))%SOLID) CYCLE
   IPZ  = PRESSURE_ZONE(BC%IIG,BC%JJG,BC%KKG)
   IOPZ = PRESSURE_ZONE(BC%II,BC%JJ,BC%KK)
   IF (IW>N_EXTERNAL_WALL_CELLS .AND. IPZ/=IOPZ) THEN
      CONNECTED_ZONES(IOPZ,IPZ) = 1
      CONNECTED_ZONES(IPZ,IOPZ) = 1
   ENDIF
   IF (WC%BOUNDARY_TYPE==OPEN_BOUNDARY) THEN
      CONNECTED_ZONES(0,IPZ) = 1
      CONNECTED_ZONES(IPZ,0) = 1
   ELSEIF (WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY) THEN
      CONNECTED_ZONES(IOPZ,IPZ) = 1
      CONNECTED_ZONES(IPZ,IOPZ) = 1
   ENDIF
ENDDO

END SUBROUTINE MERGE_PRESSURE_ZONES


SUBROUTINE PREDICT_NORMAL_VELOCITY(IW,T,DT)

USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
INTEGER, INTENT(IN) :: IW
REAL(EB), INTENT(IN) :: T,DT
REAL(EB) :: TSI,TIME_RAMP_FACTOR,PROFILE_FACTOR
TYPE(VENTS_TYPE), POINTER :: VT
TYPE(SURFACE_TYPE), POINTER :: SF
TYPE(WALL_TYPE), POINTER :: WC
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC
TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1

WC => WALL(IW)
BC => BOUNDARY_COORD(WC%BC_INDEX)
B1 => BOUNDARY_PROP1(WC%B1_INDEX)
SF => SURFACE(WC%SURF_INDEX)

PREDICT_NORMALS: IF (PREDICTOR) THEN

   WALL_CELL_TYPE: SELECT CASE (WC%BOUNDARY_TYPE)

      CASE (NULL_BOUNDARY)

         B1%U_NORMAL_S = 0._EB

      CASE (SOLID_BOUNDARY)

         IF (SF%SPECIES_BC_INDEX==SPECIFIED_MASS_FLUX .OR. SF%SPECIES_BC_INDEX==INTERPOLATED_BC .OR. &
             B1%NODE_INDEX > 0 .OR. ANY(SF%LEAK_PATH>0)) EXIT WALL_CELL_TYPE

         IF (ABS(B1%T_IGN-T_BEGIN) < SPACING(B1%T_IGN) .AND. SF%RAMP(TIME_VELO)%INDEX>=1) THEN
            TSI = T + DT
         ELSE
            TSI = T + DT - B1%T_IGN
            IF (TSI<0._EB) THEN
               B1%U_NORMAL_S = 0._EB
               EXIT WALL_CELL_TYPE
            ENDIF
         ENDIF
         TIME_RAMP_FACTOR = EVALUATE_RAMP(TSI,SF%RAMP(TIME_VELO)%INDEX,TAU=SF%RAMP(TIME_VELO)%TAU)
         B1%U_NORMAL_S = TIME_RAMP_FACTOR*B1%U_NORMAL_0

         ! Special Cases
         NEUMANN_IF: IF (SF%SPECIFIED_NORMAL_GRADIENT) THEN
            SELECT CASE(BC%IOR)
               CASE( 1); B1%U_NORMAL_S =-(U(BC%IIG,BC%JJG,BC%KKG)   + SF%VEL_GRAD*B1%RDN)
               CASE(-1); B1%U_NORMAL_S = (U(BC%IIG-1,BC%JJG,BC%KKG) + SF%VEL_GRAD*B1%RDN)
               CASE( 2); B1%U_NORMAL_S =-(V(BC%IIG,BC%JJG,BC%KKG)   + SF%VEL_GRAD*B1%RDN)
               CASE(-2); B1%U_NORMAL_S = (V(BC%IIG,BC%JJG-1,BC%KKG) + SF%VEL_GRAD*B1%RDN)
               CASE( 3); B1%U_NORMAL_S =-(W(BC%IIG,BC%JJG,BC%KKG)   + SF%VEL_GRAD*B1%RDN)
               CASE(-3); B1%U_NORMAL_S = (W(BC%IIG,BC%JJG,BC%KKG-1) + SF%VEL_GRAD*B1%RDN)
            END SELECT
         ENDIF NEUMANN_IF
         IF (ABS(SURFACE(WC%SURF_INDEX)%MASS_FLUX_TOTAL)>=TWENTY_EPSILON_EB) B1%U_NORMAL_S = &
                                                                          B1%U_NORMAL_S*RHOA/B1%RHO_F
         VENT_IF: IF (WC%VENT_INDEX>0) THEN
            VT=>VENTS(WC%VENT_INDEX)
            IF (VT%N_EDDY>0) THEN ! Synthetic Eddy Method
               IF (SF%PROFILE/=0 .AND. ABS(SF%VEL)>TWENTY_EPSILON_EB) THEN
                  PROFILE_FACTOR = ABS(B1%U_NORMAL_0/SF%VEL)
               ELSE
                  PROFILE_FACTOR = 1._EB
               ENDIF
               SELECT CASE(BC%IOR)
                  CASE( 1); B1%U_NORMAL_S = B1%U_NORMAL_S - TIME_RAMP_FACTOR*VT%U_EDDY(BC%JJ,BC%KK)*PROFILE_FACTOR
                  CASE(-1); B1%U_NORMAL_S = B1%U_NORMAL_S + TIME_RAMP_FACTOR*VT%U_EDDY(BC%JJ,BC%KK)*PROFILE_FACTOR
                  CASE( 2); B1%U_NORMAL_S = B1%U_NORMAL_S - TIME_RAMP_FACTOR*VT%V_EDDY(BC%II,BC%KK)*PROFILE_FACTOR
                  CASE(-2); B1%U_NORMAL_S = B1%U_NORMAL_S + TIME_RAMP_FACTOR*VT%V_EDDY(BC%II,BC%KK)*PROFILE_FACTOR
                  CASE( 3); B1%U_NORMAL_S = B1%U_NORMAL_S - TIME_RAMP_FACTOR*VT%W_EDDY(BC%II,BC%JJ)*PROFILE_FACTOR
                  CASE(-3); B1%U_NORMAL_S = B1%U_NORMAL_S + TIME_RAMP_FACTOR*VT%W_EDDY(BC%II,BC%JJ)*PROFILE_FACTOR
               END SELECT
            ENDIF
         ENDIF VENT_IF

      CASE(OPEN_BOUNDARY,INTERPOLATED_BOUNDARY)

         SELECT CASE(BC%IOR)
            CASE( 1); B1%U_NORMAL_S = -U(BC%II,BC%JJ,BC%KK)
            CASE(-1); B1%U_NORMAL_S =  U(BC%II-1,BC%JJ,BC%KK)
            CASE( 2); B1%U_NORMAL_S = -V(BC%II,BC%JJ,BC%KK)
            CASE(-2); B1%U_NORMAL_S =  V(BC%II,BC%JJ-1,BC%KK)
            CASE( 3); B1%U_NORMAL_S = -W(BC%II,BC%JJ,BC%KK)
            CASE(-3); B1%U_NORMAL_S =  W(BC%II,BC%JJ,BC%KK-1)
         END SELECT

   END SELECT WALL_CELL_TYPE

   ! Calculate du/dt, dv/dt, dw/dt at external boundaries. DUNDT is only used for Neumann BCs.

   IF (IW<=N_EXTERNAL_WALL_CELLS) EXTERNAL_WALL(IW)%DUNDT = (B1%U_NORMAL_S-B1%U_NORMAL)/DT

ELSE PREDICT_NORMALS

   ! In the CORRECTOR step, the normal component of velocity, U_NORMAL, is the same as the predicted value, U_NORMAL_S.
   ! However, for species mass fluxes and HVAC, U_NORMAL is computed elsewhere (wall.f90).

   IF (WC%BOUNDARY_TYPE/=SOLID_BOUNDARY .OR. (SF%SPECIES_BC_INDEX/=SPECIFIED_MASS_FLUX .AND. &
       SF%SPECIES_BC_INDEX/=INTERPOLATED_BC .AND. B1%NODE_INDEX<=0 .AND. ALL(SF%LEAK_PATH<=0))) B1%U_NORMAL = B1%U_NORMAL_S

   ! Calculate du/dt, dv/dt, dw/dt at external boundaries. DUNDT is only used for Neumann BCs. Note the second-order RK formula.

   IF (IW<=N_EXTERNAL_WALL_CELLS) EXTERNAL_WALL(IW)%DUNDT = EXTERNAL_WALL(IW)%DUNDT + 2._EB*(B1%U_NORMAL-B1%U_NORMAL_S)/DT

ENDIF PREDICT_NORMALS

END SUBROUTINE PREDICT_NORMAL_VELOCITY


!> \brief Finish computing the divergence of the flow, D, and then compute its time derivative, DDDT
!> \param DT Time step (s)
!> \param NM Mesh number

SUBROUTINE DIVERGENCE_PART_2(DT,NM)

USE COMP_FUNCTIONS, ONLY: CURRENT_TIME
USE COMPLEX_GEOMETRY, ONLY : CC_CGSC, CC_UNKZ, CC_SOLID, CC_CUTCFE
USE CC_SCALARS, ONLY : ADD_CUTCELL_D_PBAR_DT, ADD_LINKEDCELL_D_PBAR_DT, GET_CUTCELL_DDDT, GET_LINKED_VELOCITIES

INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: DT
REAL(EB), POINTER, DIMENSION(:,:,:) :: DP,RTRM,DIV
REAL(EB) :: USUM_ADD(N_ZONE),UN_P
REAL(EB) :: RDT,TNOW,P_EQ,SUM_P_PSUM,SUM_USUM,SUM_DSUM,SUM_PSUM
LOGICAL :: OPEN_ZONE
REAL(EB), POINTER, DIMENSION(:) :: D_PBAR_DT_P
REAL(EB), POINTER, DIMENSION(:,:) :: PBAR_P
INTEGER :: IW,IC,I,J,K,IPZ,IOPZ
TYPE(WALL_TYPE), POINTER :: WC
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC
TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1

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
      DP => D
      PBAR_P => PBAR
END SELECT

R_PBAR = 1._EB/PBAR_P

RTRM => WORK1

! Adjust volume flows (USUM) of pressure ZONEs that are connected to equalize background pressure

USUM_ADD = 0._EB

DO IPZ=1,N_ZONE
   SUM_P_PSUM = PBAR_P(1,IPZ)*PSUM(IPZ)
   OPEN_ZONE  = .FALSE.
   SUM_USUM = USUM(IPZ)
   SUM_DSUM = DSUM(IPZ)
   SUM_PSUM = PSUM(IPZ)
   DO IOPZ=N_ZONE,0,-1
      IF (IOPZ==IPZ) CYCLE
      IF (CONNECTED_ZONES(IPZ,IOPZ)>0) THEN
         IF (IOPZ==0) THEN
            OPEN_ZONE = .TRUE.
         ELSE
            SUM_P_PSUM = SUM_P_PSUM + PBAR_P(1,IOPZ)*PSUM(IOPZ)
            SUM_USUM = SUM_USUM + USUM(IOPZ)
            SUM_DSUM = SUM_DSUM + DSUM(IOPZ)
            SUM_PSUM = SUM_PSUM + PSUM(IOPZ)
         ENDIF
      ENDIF
   ENDDO
   IF (OPEN_ZONE) THEN
      P_EQ          = P_0(1)
      USUM_ADD(IPZ) = PSUM(IPZ)*(PBAR_P(1,IPZ)-P_EQ)/PRESSURE_RELAX_TIME + DSUM(IPZ) - USUM(IPZ)
   ELSE
      P_EQ          = SUM_P_PSUM/SUM_PSUM
      USUM_ADD(IPZ) = PSUM(IPZ)*(PBAR_P(1,IPZ)-P_EQ)/PRESSURE_RELAX_TIME + DSUM(IPZ) - USUM(IPZ) - &
                      PSUM(IPZ)*(SUM_DSUM-SUM_USUM)/SUM_PSUM
   ENDIF
ENDDO

DO IPZ=1,N_ZONE
   USUM(IPZ) = USUM(IPZ) + USUM_ADD(IPZ)
ENDDO

! Compute dP/dt for each pressure ZONE

IF_PRESSURE_ZONES: IF (N_ZONE>0) THEN

   IF (PREDICTOR) D_PBAR_DT_P => D_PBAR_DT_S
   IF (CORRECTOR) D_PBAR_DT_P => D_PBAR_DT

   ! Compute change in background pressure

   DO IPZ=1,N_ZONE
      IF (ABS(PSUM(IPZ)) > TWENTY_EPSILON_EB) D_PBAR_DT_P(IPZ) = (DSUM(IPZ) - USUM(IPZ))/PSUM(IPZ)
      IF (CORRECTOR) P_ZONE(IPZ)%DPSTAR =  D_PBAR_DT_P(IPZ)
   ENDDO

   ! Add pressure derivative to divergence

   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IPZ = PRESSURE_ZONE(I,J,K)
            IF (IPZ<1) CYCLE
            IF (CELL(CELL_INDEX(I,J,K))%SOLID) CYCLE
            IF (CC_IBM) THEN
               IF (CCVAR(I,J,K,CC_CGSC) == CC_SOLID) THEN
                  CYCLE
               ELSEIF(CCVAR(I,J,K,CC_CGSC) == CC_CUTCFE) THEN
                  CALL ADD_CUTCELL_D_PBAR_DT(I,J,K,PBAR_P(K,IPZ),D_PBAR_DT_P(IPZ)); CYCLE
               ELSEIF(CCVAR(I,J,K,CC_UNKZ) > 0) THEN
                  CALL ADD_LINKEDCELL_D_PBAR_DT(I,J,K,PBAR_P(K,IPZ),D_PBAR_DT_P(IPZ),RTRM(I,J,K),DP(I,J,K)); CYCLE
               ENDIF
            ENDIF
            DP(I,J,K) = DP(I,J,K) - (R_PBAR(K,IPZ)-RTRM(I,J,K))*D_PBAR_DT_P(IPZ)
         ENDDO
      ENDDO
   ENDDO

ENDIF IF_PRESSURE_ZONES

! Zero out divergence in solid cells

SOLID_LOOP: DO IC=1,CELL_COUNT(NM)
   IF (.NOT.CELL(IC)%SOLID) CYCLE SOLID_LOOP
   I = CELL(IC)%I
   J = CELL(IC)%J
   K = CELL(IC)%K
   DP(I,J,K) = 0._EB
ENDDO SOLID_LOOP

! Zero out CC_IBM solid cells:
IF (CC_IBM) THEN
   T_USED(2)=T_USED(2)+CURRENT_TIME()-TNOW
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (CCVAR(I,J,K,CC_CGSC) /= CC_SOLID) CYCLE
            DP(I,J,K) = 0._EB
         ENDDO
      ENDDO
   ENDDO
   CALL GET_LINKED_VELOCITIES(NM,CORRECTOR,CMP_FLG=.FALSE.)
   TNOW=CURRENT_TIME()
ENDIF

! Specify divergence in boundary cells to account for volume being generated at the walls

BC_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
   WC => WALL(IW)
   IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE BC_LOOP
   BC=>BOUNDARY_COORD(WC%BC_INDEX)
   SELECT CASE (WC%BOUNDARY_TYPE)
      CASE (SOLID_BOUNDARY)
         IF (.NOT.CELL(CELL_INDEX(BC%II,BC%JJ,BC%KK))%SOLID) CYCLE BC_LOOP
         B1 => BOUNDARY_PROP1(WC%BC_INDEX)
         IF (PREDICTOR) THEN
            UN_P = B1%U_NORMAL_S
         ELSE
            UN_P = B1%U_NORMAL
         ENDIF
         SELECT CASE(BC%IOR)
            CASE( 1)
               DP(BC%II,BC%JJ,BC%KK) = DP(BC%II,BC%JJ,BC%KK) - UN_P*RDX(BC%II)*RRN(BC%II)*R(BC%II)
            CASE(-1)
               DP(BC%II,BC%JJ,BC%KK) = DP(BC%II,BC%JJ,BC%KK) - UN_P*RDX(BC%II)*RRN(BC%II)*R(BC%II-1)
            CASE( 2)
               DP(BC%II,BC%JJ,BC%KK) = DP(BC%II,BC%JJ,BC%KK) - UN_P*RDY(BC%JJ)
            CASE(-2)
               DP(BC%II,BC%JJ,BC%KK) = DP(BC%II,BC%JJ,BC%KK) - UN_P*RDY(BC%JJ)
            CASE( 3)
               DP(BC%II,BC%JJ,BC%KK) = DP(BC%II,BC%JJ,BC%KK) - UN_P*RDZ(BC%KK)
            CASE(-3)
               DP(BC%II,BC%JJ,BC%KK) = DP(BC%II,BC%JJ,BC%KK) - UN_P*RDZ(BC%KK)
         END SELECT
      CASE (OPEN_BOUNDARY,MIRROR_BOUNDARY,INTERPOLATED_BOUNDARY)
         DP(BC%II,BC%JJ,BC%KK) = DP(BC%IIG,BC%JJG,BC%KKG)
   END SELECT
ENDDO BC_LOOP

! Compute time derivative of the divergence, dD/dt

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
   DDDT = (2._EB*DP-DIV)*RDT
ENDIF

T_USED(2)=T_USED(2)+CURRENT_TIME()-TNOW

IF(CC_IBM) CALL GET_CUTCELL_DDDT(DT,NM)

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

IF(STORE_CARTESIAN_DIVERGENCE) CARTVELDIV = DP

RESMAX = 0._EB
DIVMX  = -10000._EB
DIVMN  =  10000._EB
IMX    = 0
JMX    = 0
KMX    = 0

DO K=1,KBAR
   DO J=1,JBAR
      LOOP1: DO I=1,IBAR
         IF (CELL(CELL_INDEX(I,J,K))%SOLID) CYCLE LOOP1
         SELECT CASE(CYLINDRICAL)
            CASE(.FALSE.)
               DIV = (UU(I,J,K)-UU(I-1,J,K))*RDX(I) + &
                     (VV(I,J,K)-VV(I,J-1,K))*RDY(J) + &
                     (WW(I,J,K)-WW(I,J,K-1))*RDZ(K)
            CASE(.TRUE.)
               DIV = (R(I)*UU(I,J,K)-R(I-1)*UU(I-1,J,K))*RDX(I)*RRN(I) +  &
                     (WW(I,J,K)-WW(I,J,K-1))*RDZ(K)
         END SELECT
         IF(STORE_CARTESIAN_DIVERGENCE) CARTVELDIV(I,J,K) = DIV
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
