MODULE DIVG              
 
USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_POINTERS
 
IMPLICIT NONE
PRIVATE
CHARACTER(255), PARAMETER :: divgid='$Id$'
CHARACTER(255), PARAMETER :: divgrev='$Revision$'
CHARACTER(255), PARAMETER :: divgdate='$Date$'

PUBLIC DIVERGENCE_PART_1,DIVERGENCE_PART_2,CHECK_DIVERGENCE,GET_REV_divg
 
CONTAINS
 
SUBROUTINE DIVERGENCE_PART_1(T,NM)

USE COMP_FUNCTIONS, ONLY: SECOND 
USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
USE PHYSICAL_FUNCTIONS, ONLY: GET_CONDUCTIVITY,GET_SPECIFIC_HEAT,GET_SENSIBLE_ENTHALPY_DIFF,GET_SENSIBLE_ENTHALPY,&
                              GET_VISCOSITY
USE TURBULENCE, ONLY: TENSOR_DIFFUSIVITY_MODEL,WANNIER_FLOW
USE MASS, ONLY: SCALAR_FACE_VALUE
USE GEOMETRY_FUNCTIONS, ONLY: ASSIGN_PRESSURE_ZONE
USE MANUFACTURED_SOLUTIONS, ONLY: DIFF_MMS,UF_MMS,WF_MMS,VD2D_MMS_Z_SRC,RHO_0_MMS,RHO_1_MMS

! Compute contributions to the divergence term
 
INTEGER, INTENT(IN) :: NM
REAL(EB), POINTER, DIMENSION(:,:,:) :: KDTDX,KDTDY,KDTDZ,DP,KP,CP, &
          RHO_D_DZDX,RHO_D_DZDY,RHO_D_DZDZ,RHO_D,RHOP,H_RHO_D_DZDX,H_RHO_D_DZDY,H_RHO_D_DZDZ,RTRM, &
          U_DOT_DEL_RHO_H_S,RHO_H_S_P,UU,VV,WW,U_DOT_DEL_RHO,RHO_Z_P,U_DOT_DEL_RHO_Z,RHO_D_TURB,R_H_G
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZP
REAL(EB), POINTER, DIMENSION(:,:) :: PBAR_P            
REAL(EB) :: DELKDELT,VC,VC1,DTDX,DTDY,DTDZ,TNOW, &
            HDIFF,DZDX,DZDY,DZDZ,T,RDT,RHO_D_DZDN,TSI,TIME_RAMP_FACTOR,ZONE_VOLUME,DELTA_P,PRES_RAMP_FACTOR,&
            TMP_G,TMP_WGT,DIV_DIFF_HEAT_FLUX,H_S,ZZZ(1:4),DU_P,DU_M,UN,RCON_DIFF,PROFILE_FACTOR, &
            XHAT,ZHAT,TT,Q_Z
REAL(EB), ALLOCATABLE, DIMENSION(:) :: ZZ_GET
TYPE(SURFACE_TYPE), POINTER :: SF
TYPE(SPECIES_MIXTURE_TYPE), POINTER :: SM,SM0
INTEGER :: IW,N,IOR,II,JJ,KK,IIG,JJG,KKG,ITMP,I,J,K,IPZ,IOPZ
TYPE(VENTS_TYPE), POINTER :: VT=>NULL()
TYPE(WALL_TYPE), POINTER :: WC=>NULL()
REAL(EB), PARAMETER :: ADVECTION_EPS=1.E-6_EB
!TYPE(FACET_TYPE), POINTER :: FC=>NULL()
!TYPE(CUTCELL_LINKED_LIST_TYPE), POINTER :: CL=>NULL()
!INTEGER :: NF,IC
 
! Check whether to skip this routine

IF (SOLID_PHASE_ONLY) RETURN

IF (EVACUATION_ONLY(NM)) THEN
   CALL DIVERGENCE_PART_1_EVAC(T,NM)
   RETURN
ENDIF

! Start the clock and set the pointers

TNOW=SECOND()
CALL POINT_TO_MESH(NM)
 
RDT = 1._EB/DT
ALLOCATE(ZZ_GET(0:N_TRACKED_SPECIES))
 
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

CONNECTED_ZONES(:,:,NM) = .FALSE.

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

! Compute species-related finite difference terms

IF (N_TRACKED_SPECIES > 0) THEN
   RHO_D_DZDX  => WORK1
   RHO_D_DZDY  => WORK2
   RHO_D_DZDZ  => WORK3

   SELECT CASE(PREDICTOR)
      CASE(.TRUE.)  
         ZZP => ZZS 
      CASE(.FALSE.) 
         ZZP => ZZ  
   END SELECT
ENDIF

! Add species diffusion terms to divergence expression and compute diffusion term for species equations

IF (N_TRACKED_SPECIES > 0) THEN
   DEL_RHO_D_DEL_Z = 0._EB
   RHO_D => WORK4
   IF (LES) THEN
      IF (.NOT.RESEARCH_MODE) THEN
         RHO_D = MU*RSC
      ELSE
         RHO_D_TURB => WORK9
         RHO_D_TURB = 0._EB
         RHO_D_TURB = MAX(0._EB,MU-MU_DNS)*RSC
      ENDIF
   ENDIF
   IF (CHECK_VN) D_Z_MAX = 0._EB
ENDIF

SPECIES_LOOP: DO N=1,N_TRACKED_SPECIES

   IF (DNS .OR. RESEARCH_MODE) THEN
      RHO_D = 0._EB
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               ITMP = MIN(4999,INT(TMP(I,J,K)))
               TMP_WGT = TMP(I,J,K) - ITMP
               RHO_D(I,J,K) = RHOP(I,J,K)*(D_Z(ITMP,N)+TMP_WGT*(D_Z(ITMP+1,N)-D_Z(ITMP,N)))
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   IF (LES .AND. RESEARCH_MODE) RHO_D = RHO_D + RHO_D_TURB

   ! Store diffusivity for stability check

   IF (CHECK_VN) D_Z_MAX = MAX(D_Z_MAX,RHO_D/RHOP)

   ! Manufactured solution

   IF (PERIODIC_TEST==7) RHO_D = DIFF_MMS

   ! Compute rho*D del Z

   DO K=0,KBAR
      DO J=0,JBAR
         DO I=0,IBAR
            DZDX = (ZZP(I+1,J,K,N)-ZZP(I,J,K,N))*RDXN(I)
            RHO_D_DZDX(I,J,K) = .5_EB*(RHO_D(I+1,J,K)+RHO_D(I,J,K))*DZDX
            DZDY = (ZZP(I,J+1,K,N)-ZZP(I,J,K,N))*RDYN(J)
            RHO_D_DZDY(I,J,K) = .5_EB*(RHO_D(I,J+1,K)+RHO_D(I,J,K))*DZDY
            DZDZ = (ZZP(I,J,K+1,N)-ZZP(I,J,K,N))*RDZN(K)
            RHO_D_DZDZ(I,J,K) = .5_EB*(RHO_D(I,J,K+1)+RHO_D(I,J,K))*DZDZ
         ENDDO
      ENDDO
   ENDDO

   ! Tensor diffusivity model (experimental)

   IF (TENSOR_DIFFUSIVITY .AND. LES) CALL TENSOR_DIFFUSIVITY_MODEL(NM,N)

   ! Compute del dot h_n*rho*D del Z_n (part of del dot qdot")

   H_RHO_D_DZDX => WORK5
   H_RHO_D_DZDY => WORK6
   H_RHO_D_DZDZ => WORK7

   !$OMP PARALLEL DO PRIVATE(TMP_G, HDIFF) SCHEDULE(guided)
   DO K=0,KBAR
      DO J=0,JBAR
         DO I=0,IBAR
            ! H_RHO_D_DZDX
            TMP_G = 0.5_EB*(TMP(I+1,J,K)+TMP(I,J,K))
            CALL GET_SENSIBLE_ENTHALPY_DIFF(N,TMP_G,HDIFF)
            H_RHO_D_DZDX(I,J,K) = HDIFF*RHO_D_DZDX(I,J,K)

            ! H_RHO_D_DZDY
            TMP_G = 0.5_EB*(TMP(I,J+1,K)+TMP(I,J,K))
            CALL GET_SENSIBLE_ENTHALPY_DIFF(N,TMP_G,HDIFF)
            H_RHO_D_DZDY(I,J,K) = HDIFF*RHO_D_DZDY(I,J,K)

            ! H_RHO_D_DZDZ
            TMP_G = 0.5_EB*(TMP(I,J,K+1)+TMP(I,J,K))               
            CALL GET_SENSIBLE_ENTHALPY_DIFF(N,TMP_G,HDIFF)
            H_RHO_D_DZDZ(I,J,K) = HDIFF*RHO_D_DZDZ(I,J,K)
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO

   ! Correct rho*D del Z and del dot h_n*rho*D del Z_n at boundaries and store rho*D at boundaries

   WALL_LOOP2: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
      WC => WALL(IW)
      IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY .OR. &
          WC%BOUNDARY_TYPE==OPEN_BOUNDARY .OR. WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY) CYCLE WALL_LOOP2
      IIG = WC%ONE_D%IIG
      JJG = WC%ONE_D%JJG
      KKG = WC%ONE_D%KKG
      IOR = WC%ONE_D%IOR
      CALL GET_SENSIBLE_ENTHALPY_DIFF(N,WC%ONE_D%TMP_F,HDIFF)
      RHO_D_DZDN = 2._EB*WC%RHODW(N)*(ZZP(IIG,JJG,KKG,N)-WC%ZZ_F(N))*WC%RDN
      SELECT CASE(IOR)
         CASE( 1) 
            RHO_D_DZDX(IIG-1,JJG,KKG)   =  RHO_D_DZDN
            H_RHO_D_DZDX(IIG-1,JJG,KKG) =  HDIFF*RHO_D_DZDN
         CASE(-1) 
            RHO_D_DZDX(IIG,JJG,KKG)     = -RHO_D_DZDN
            H_RHO_D_DZDX(IIG,JJG,KKG)   = -HDIFF*RHO_D_DZDN
         CASE( 2) 
            RHO_D_DZDY(IIG,JJG-1,KKG)   =  RHO_D_DZDN
            H_RHO_D_DZDY(IIG,JJG-1,KKG) =  HDIFF*RHO_D_DZDN
         CASE(-2) 
            RHO_D_DZDY(IIG,JJG,KKG)     = -RHO_D_DZDN
            H_RHO_D_DZDY(IIG,JJG,KKG)   = -HDIFF*RHO_D_DZDN
         CASE( 3) 
            RHO_D_DZDZ(IIG,JJG,KKG-1)   =  RHO_D_DZDN
            H_RHO_D_DZDZ(IIG,JJG,KKG-1) =  HDIFF*RHO_D_DZDN
         CASE(-3) 
            RHO_D_DZDZ(IIG,JJG,KKG)     = -RHO_D_DZDN
            H_RHO_D_DZDZ(IIG,JJG,KKG)   = -HDIFF*RHO_D_DZDN
      END SELECT
   ENDDO WALL_LOOP2

   CYLINDER: SELECT CASE(CYLINDRICAL)
      CASE(.FALSE.) CYLINDER  ! 3D or 2D Cartesian Coords
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR

                  DIV_DIFF_HEAT_FLUX = (H_RHO_D_DZDX(I,J,K)-H_RHO_D_DZDX(I-1,J,K))*RDX(I) + &
                                       (H_RHO_D_DZDY(I,J,K)-H_RHO_D_DZDY(I,J-1,K))*RDY(J) + &
                                       (H_RHO_D_DZDZ(I,J,K)-H_RHO_D_DZDZ(I,J,K-1))*RDZ(K)

                  DP(I,J,K) = DP(I,J,K) + DIV_DIFF_HEAT_FLUX

               ENDDO
            ENDDO
         ENDDO
      CASE(.TRUE.) CYLINDER  ! 2D Cylindrical Coords
         J = 1
         DO K=1,KBAR
            DO I=1,IBAR

               DIV_DIFF_HEAT_FLUX = (R(I)*H_RHO_D_DZDX(I,J,K)-R(I-1)*H_RHO_D_DZDX(I-1,J,K))*RDX(I)*RRN(I) + &
                                    (     H_RHO_D_DZDZ(I,J,K)-       H_RHO_D_DZDZ(I,J,K-1))*RDZ(K)

               DP(I,J,K) = DP(I,J,K) + DIV_DIFF_HEAT_FLUX

            ENDDO
         ENDDO
   END SELECT CYLINDER

   ! Compute del dot rho*D del Z_n

   CYLINDER2: SELECT CASE(CYLINDRICAL)
      CASE(.FALSE.) CYLINDER2  ! 3D or 2D Cartesian Coords
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  DEL_RHO_D_DEL_Z(I,J,K,N) = (RHO_D_DZDX(I,J,K)-RHO_D_DZDX(I-1,J,K))*RDX(I) + &
                                             (RHO_D_DZDY(I,J,K)-RHO_D_DZDY(I,J-1,K))*RDY(J) + &
                                             (RHO_D_DZDZ(I,J,K)-RHO_D_DZDZ(I,J,K-1))*RDZ(K)
               ENDDO
            ENDDO
         ENDDO
      CASE(.TRUE.) CYLINDER2  ! 2D Cylindrical Coords
         J=1
         DO K=1,KBAR
            DO I=1,IBAR
               DEL_RHO_D_DEL_Z(I,J,K,N) = (R(I)*RHO_D_DZDX(I,J,K)-R(I-1)*RHO_D_DZDX(I-1,J,K))*RDX(I)*RRN(I) + &
                                          (     RHO_D_DZDZ(I,J,K)-       RHO_D_DZDZ(I,J,K-1))*RDZ(K)
            ENDDO
         ENDDO
   END SELECT CYLINDER2

ENDDO SPECIES_LOOP

! Get the specific heat

CP => WORK5
R_H_G=> WORK9

IF (N_TRACKED_SPECIES>0) THEN
   !$ DEALLOCATE(ZZ_GET)
   !$OMP PARALLEL PRIVATE(ZZ_GET, H_S)
   !$ ALLOCATE(ZZ_GET(0:N_TRACKED_SPECIES))
   !$OMP DO SCHEDULE(static)
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
   !$ DEALLOCATE(ZZ_GET)
   !$OMP END PARALLEL
   !$ ALLOCATE(ZZ_GET(0:N_TRACKED_SPECIES))
ELSE
   ZZ_GET = 0._EB
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
            CALL GET_SPECIFIC_HEAT(ZZ_GET,CP(I,J,K),TMP(I,J,K))
            R_H_G(I,J,K) = 1._EB/(CP(I,J,K)*TMP(I,J,K))
         ENDDO
      ENDDO
   ENDDO
ENDIF

! Compute del dot k del T

KDTDX => WORK1
KDTDY => WORK2
KDTDZ => WORK3
KP    => WORK4

! Compute thermal conductivity k (KP)

K_DNS_OR_LES: IF (DNS .OR. RESEARCH_MODE) THEN

   KP = 0._EB
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
            IF (N_TRACKED_SPECIES > 0) ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I,J,K,1:N_TRACKED_SPECIES)
            CALL GET_CONDUCTIVITY(ZZ_GET,KP(I,J,K),TMP(I,J,K)) 
         ENDDO
      ENDDO
   ENDDO

   IF (RESEARCH_MODE) THEN
      KP = KP + MAX(0._EB,(MU-MU_DNS))*CP*RPR
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

   KP = MU*CPOPR

ENDIF K_DNS_OR_LES

! Compute k*dT/dx, etc

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
      WC%KW = 0.5_EB*(KP(IIG,JJG,KKG)+KP(II,JJ,KK))
      CYCLE CORRECTION_LOOP
   ELSE
      WC%KW = KP(IIG,JJG,KKG)
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
   DP(IIG,JJG,KKG) = DP(IIG,JJG,KKG) - WC%ONE_D%QCONF*WC%RDN
ENDDO CORRECTION_LOOP

! Compute (q + del dot k del T) and add to the divergence

CYLINDER3: SELECT CASE(CYLINDRICAL)
CASE(.FALSE.) CYLINDER3   ! 3D or 2D Cartesian
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
CASE(.TRUE.) CYLINDER3   ! 2D Cylindrical
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            DELKDELT = & 
                 (R(I)*KDTDX(I,J,K)-R(I-1)*KDTDX(I-1,J,K))*RDX(I)*RRN(I) + &
                 (KDTDZ(I,J,K)-       KDTDZ(I,J,K-1))*RDZ(K)
            DP(I,J,K) = DP(I,J,K) + DELKDELT + Q(I,J,K) + QR(I,J,K)
         ENDDO
      ENDDO
   ENDDO
END SELECT CYLINDER3

! Compute and store rho*h_s

IF (.NOT.CONSTANT_SPECIFIC_HEAT_RATIO) THEN

   RHO_H_S_P=>WORK1
   IF (N_TRACKED_SPECIES > 0) THEN
      !$ DEALLOCATE(ZZ_GET)
      !$OMP PARALLEL PRIVATE(ZZ_GET, H_S)
      !$ ALLOCATE(ZZ_GET(0:N_TRACKED_SPECIES))
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
      !$ DEALLOCATE(ZZ_GET)
      !$OMP END PARALLEL
      !$ ALLOCATE(ZZ_GET(0:N_TRACKED_SPECIES))
   ELSE
      ZZ_GET = 0._EB
      DO K=0,KBP1
         DO J=0,JBP1
            DO I=0,IBP1
               CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,TMP(I,J,K))
               RHO_H_S_P(I,J,K) = RHOP(I,J,K)*H_S
            ENDDO
         ENDDO
      ENDDO
   ENDIF

ENDIF

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

! Compute U_DOT_DEL_RHO_H_S and add to other enthalpy equation source terms

IF (.NOT.CONSTANT_SPECIFIC_HEAT_RATIO) THEN

   CALL ENTHALPY_ADVECTION 

   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
            DP(I,J,K) = DP(I,J,K) - U_DOT_DEL_RHO_H_S(I,J,K)
         ENDDO
      ENDDO
   ENDDO

ENDIF

! Compute RTRM = 1/(rho*c_p*T) and multiply it by divergence terms already summed up

DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
         RTRM(I,J,K) = R_H_G(I,J,K)/RHOP(I,J,K)
         DP(I,J,K) = RTRM(I,J,K)*DP(I,J,K)
      ENDDO
   ENDDO
ENDDO

! Compute mass flux for density equation, U dot grad rho

CALL DENSITY_ADVECTION

! Compute (1/rho) * Sum( (Wbar/W_alpha-h_s,alpha/cp*T) (del dot rho*D del Z_n - u dot del rho*Z_n)

IF (.NOT.CONSTANT_SPECIFIC_HEAT_RATIO) THEN
   SM0 => SPECIES_MIXTURE(0)
   !$ DEALLOCATE(ZZ_GET)
   !$OMP PARALLEL PRIVATE(ZZ_GET, H_S)
   !$ ALLOCATE(ZZ_GET(0:N_TRACKED_SPECIES))
   !$OMP DO SCHEDULE(static)
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
            ZZ_GET = 0._EB
            CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,TMP(I,J,K))
            DP(I,J,K) = DP(I,J,K) - ( SM0%RCON/RSUM(I,J,K) - H_S*R_H_G(I,J,K) )*U_DOT_DEL_RHO(I,J,K)/RHOP(I,J,K)
         ENDDO
      ENDDO
   ENDDO
   !$OMP END DO
   !$ DEALLOCATE(ZZ_GET)
   !$OMP END PARALLEL
   !$ ALLOCATE(ZZ_GET(0:N_TRACKED_SPECIES))
ENDIF

DO N=1,N_TRACKED_SPECIES
   CALL SPECIES_ADVECTION ! compute mass flux for species transport equation
   IF (CONSTANT_SPECIFIC_HEAT_RATIO) CYCLE
   SM  => SPECIES_MIXTURE(N)
   RCON_DIFF = SM%RCON-SM0%RCON
   !$OMP PARALLEL DO PRIVATE(HDIFF) SCHEDULE(guided)
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
            CALL GET_SENSIBLE_ENTHALPY_DIFF(N,TMP(I,J,K),HDIFF)            
            DP(I,J,K) = DP(I,J,K) + (RCON_DIFF/RSUM(I,J,K) - HDIFF*R_H_G(I,J,K))* &
                 ( DEL_RHO_D_DEL_Z(I,J,K,N) - U_DOT_DEL_RHO_Z(I,J,K) )/RHOP(I,J,K)
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO
ENDDO

! Add contribution of reactions

IF (N_REACTIONS > 0 .AND. .NOT.CONSTANT_SPECIFIC_HEAT_RATIO) THEN
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            DP(I,J,K) = DP(I,J,K) + D_REACTION(I,J,K)
         ENDDO
      ENDDO
   ENDDO
ENDIF

! Add contribution of evaporating particles

IF (CALC_D_LAGRANGIAN) THEN
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            DP(I,J,K) = DP(I,J,K) + D_LAGRANGIAN(I,J,K)
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
            DP(I,J,K) = DP(I,J,K) + RTRM(I,J,K)*0.5_EB*(W(I,J,K)+W(I,J,K-1))*RHO_0(K)*GVEC(3)
         ENDDO
      ENDDO
   ENDDO
ENDIF

! Test manufactured solution

IF (PERIODIC_TEST==7) THEN
   IF (PREDICTOR) TT=T+DT
   IF (CORRECTOR) TT=T
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            ! divergence from EOS
            XHAT = XC(I) - UF_MMS*TT
            ZHAT = ZC(K) - WF_MMS*TT
            Q_Z = VD2D_MMS_Z_SRC(XHAT,ZHAT,TT)
            DP(I,J,K) = (1._EB/RHO_1_MMS - 1._EB/RHO_0_MMS) * ( DEL_RHO_D_DEL_Z(I,J,K,1) + Q_Z )
         ENDDO
      ENDDO
   ENDDO
ENDIF

! Compute normal component of velocity at boundaries, UWS

PREDICT_NORMALS: IF (PREDICTOR) THEN
 
   WALL_LOOP3: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
      WC => WALL(IW)
      IOR = WC%ONE_D%IOR
      
      WALL_CELL_TYPE: SELECT CASE (WC%BOUNDARY_TYPE)
         CASE (NULL_BOUNDARY)
            WC%ONE_D%UWS = 0._EB
         CASE (SOLID_BOUNDARY)
            SF => SURFACE(WC%SURF_INDEX)

            IF (SF%SPECIES_BC_INDEX==SPECIFIED_MASS_FLUX .OR. &
                SF%SPECIES_BC_INDEX==INTERPOLATED_BC     .OR. &
                SF%SPECIES_BC_INDEX==HVAC_BOUNDARY       .OR. &
                ANY(SF%LEAK_PATH>0._EB))                      &
                CYCLE WALL_LOOP3
            
            IF (ABS(WC%ONE_D%T_IGN-T_BEGIN) < SPACING(WC%ONE_D%T_IGN) .AND. SF%RAMP_INDEX(TIME_VELO)>=1) THEN
               TSI = T + DT
            ELSE
               TSI = T + DT - WC%ONE_D%T_IGN
               IF (TSI<0._EB) THEN
                  WC%ONE_D%UWS = 0._EB
                  CYCLE WALL_LOOP3
               ENDIF
            ENDIF
            TIME_RAMP_FACTOR = EVALUATE_RAMP(TSI,SF%TAU(TIME_VELO),SF%RAMP_INDEX(TIME_VELO))
            KK               = WC%ONE_D%KK
            DELTA_P          = PBAR_P(KK,SF%DUCT_PATH(1)) - PBAR_P(KK,SF%DUCT_PATH(2))
            PRES_RAMP_FACTOR = SIGN(1._EB,SF%MAX_PRESSURE-DELTA_P)*SQRT(ABS((DELTA_P-SF%MAX_PRESSURE)/SF%MAX_PRESSURE))
            SELECT CASE(IOR) 
               CASE( 1)
                  WC%ONE_D%UWS =-U0 + TIME_RAMP_FACTOR*(WC%UW0+U0)
               CASE(-1)
                  WC%ONE_D%UWS = U0 + TIME_RAMP_FACTOR*(WC%UW0-U0)
               CASE( 2)
                  WC%ONE_D%UWS =-V0 + TIME_RAMP_FACTOR*(WC%UW0+V0)
               CASE(-2)
                  WC%ONE_D%UWS = V0 + TIME_RAMP_FACTOR*(WC%UW0-V0)
               CASE( 3)
                  WC%ONE_D%UWS =-W0 + TIME_RAMP_FACTOR*(WC%UW0+W0)
               CASE(-3)
                  WC%ONE_D%UWS = W0 + TIME_RAMP_FACTOR*(WC%UW0-W0)
            END SELECT
            ! Special Cases
            NEUMANN_IF: IF (SF%SPECIFIED_NORMAL_GRADIENT) THEN
               IIG = WC%ONE_D%IIG
               JJG = WC%ONE_D%JJG
               KKG = WC%ONE_D%KKG
               SELECT CASE(IOR) 
                  CASE( 1)
                     WC%ONE_D%UWS =-(U(IIG,JJG,KKG)   + SF%VEL_GRAD*WC%RDN)
                  CASE(-1)
                     WC%ONE_D%UWS = (U(IIG-1,JJG,KKG) + SF%VEL_GRAD*WC%RDN)
                  CASE( 2)
                     WC%ONE_D%UWS =-(V(IIG,JJG,KKG)   + SF%VEL_GRAD*WC%RDN)
                  CASE(-2)
                     WC%ONE_D%UWS = (V(IIG,JJG-1,KKG) + SF%VEL_GRAD*WC%RDN)
                  CASE( 3)
                     WC%ONE_D%UWS =-(W(IIG,JJG,KKG)   + SF%VEL_GRAD*WC%RDN)
                  CASE(-3)
                     WC%ONE_D%UWS = (W(IIG,JJG,KKG-1) + SF%VEL_GRAD*WC%RDN)
               END SELECT
            ENDIF NEUMANN_IF
            IF (ABS(SURFACE(WC%SURF_INDEX)%MASS_FLUX_TOTAL)>=TWO_EPSILON_EB) WC%ONE_D%UWS = WC%ONE_D%UWS*RHOA/WC%RHO_F
            VENT_IF: IF (WC%VENT_INDEX>0) THEN 
               VT=>VENTS(WC%VENT_INDEX)
               IF (VT%N_EDDY>0) THEN ! Synthetic Eddy Method
                  II = WC%ONE_D%II
                  JJ = WC%ONE_D%JJ
                  KK = WC%ONE_D%KK
                  IF (SF%PROFILE/=0 .AND. ABS(SF%VEL)>TWO_EPSILON_EB) THEN
                     PROFILE_FACTOR = ABS(WC%ONE_D%UWS/SF%VEL)
                  ELSE
                     PROFILE_FACTOR = 1._EB
                  ENDIF
                  SELECT CASE(IOR)
                     CASE( 1)
                        WC%ONE_D%UWS = WC%ONE_D%UWS - TIME_RAMP_FACTOR*VT%U_EDDY(JJ,KK)*PROFILE_FACTOR
                     CASE(-1)
                        WC%ONE_D%UWS = WC%ONE_D%UWS + TIME_RAMP_FACTOR*VT%U_EDDY(JJ,KK)*PROFILE_FACTOR
                     CASE( 2)
                        WC%ONE_D%UWS = WC%ONE_D%UWS - TIME_RAMP_FACTOR*VT%V_EDDY(II,KK)*PROFILE_FACTOR
                     CASE(-2)
                        WC%ONE_D%UWS = WC%ONE_D%UWS + TIME_RAMP_FACTOR*VT%V_EDDY(II,KK)*PROFILE_FACTOR
                     CASE( 3)
                        WC%ONE_D%UWS = WC%ONE_D%UWS - TIME_RAMP_FACTOR*VT%W_EDDY(II,JJ)*PROFILE_FACTOR
                     CASE(-3)
                        WC%ONE_D%UWS = WC%ONE_D%UWS + TIME_RAMP_FACTOR*VT%W_EDDY(II,JJ)*PROFILE_FACTOR
                  END SELECT
               ENDIF
               WANNIER_BC: IF (PERIODIC_TEST==5) THEN
                  II = WC%ONE_D%II
                  JJ = WC%ONE_D%JJ
                  KK = WC%ONE_D%KK
                  SELECT CASE(IOR)
                     CASE( 1)
                        WC%ONE_D%UWS = -WANNIER_FLOW(X(II),ZC(KK),1)
                     CASE(-1)
                        WC%ONE_D%UWS =  WANNIER_FLOW(X(II-1),ZC(KK),1)
                     CASE( 3)
                        WC%ONE_D%UWS = -WANNIER_FLOW(XC(II),Z(KK),2)
                     CASE(-3)
                        WC%ONE_D%UWS =  WANNIER_FLOW(XC(II),Z(KK-1),2)
                  END SELECT
               ENDIF WANNIER_BC
            ENDIF VENT_IF
         CASE(OPEN_BOUNDARY,INTERPOLATED_BOUNDARY)
            II = WC%ONE_D%II
            JJ = WC%ONE_D%JJ
            KK = WC%ONE_D%KK
            SELECT CASE(IOR)
               CASE( 1)
                  WC%ONE_D%UWS = -U(II,JJ,KK)
               CASE(-1)
                  WC%ONE_D%UWS =  U(II-1,JJ,KK)
               CASE( 2)
                  WC%ONE_D%UWS = -V(II,JJ,KK)
               CASE(-2)
                  WC%ONE_D%UWS =  V(II,JJ-1,KK)
               CASE( 3)
                  WC%ONE_D%UWS = -W(II,JJ,KK)
               CASE(-3)
                  WC%ONE_D%UWS =  W(II,JJ,KK-1)
            END SELECT
      END SELECT WALL_CELL_TYPE
   ENDDO WALL_LOOP3

   DO IW=1,N_EXTERNAL_WALL_CELLS
      WALL(IW)%DUWDT = RDT*(WALL(IW)%ONE_D%UWS-WALL(IW)%ONE_D%UW)
   ENDDO
   
ELSE PREDICT_NORMALS
   
   DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
      WALL(IW)%ONE_D%UW = WALL(IW)%ONE_D%UWS
   ENDDO

ENDIF PREDICT_NORMALS

! Calculate pressure rise in each of the pressure zones by summing divergence expression over each zone

PRESSURE_ZONE_LOOP: DO IPZ=1,N_ZONE

   USUM(IPZ,NM) = 0._EB
   DSUM(IPZ,NM) = 0._EB
   PSUM(IPZ,NM) = 0._EB
   ZONE_VOLUME  = 0._EB

   DO K=1,KBAR
      DO J=1,JBAR
         VC1 = DY(J)*DZ(K)
         DO I=1,IBAR
            IF (PRESSURE_ZONE(I,J,K) /= IPZ) CYCLE
            IF (SOLID(CELL_INDEX(I,J,K)))    CYCLE
            VC = DX(I)*RC(I)*VC1
            ZONE_VOLUME = ZONE_VOLUME + VC
            DSUM(IPZ,NM) = DSUM(IPZ,NM) + VC*DP(I,J,K)
            PSUM(IPZ,NM) = PSUM(IPZ,NM) + VC*(R_PBAR(K,IPZ)-RTRM(I,J,K))
         ENDDO
      ENDDO
   ENDDO

   ! Calculate the volume flux to the boundary of the pressure zone (int u dot dA)

   WALL_LOOP4: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
      IF (WALL(IW)%PRESSURE_ZONE/=IPZ) CYCLE WALL_LOOP4
      IF (WALL(IW)%BOUNDARY_TYPE/=SOLID_BOUNDARY .AND. WALL(IW)%BOUNDARY_TYPE/=HVAC_BOUNDARY) CYCLE WALL_LOOP4
      USUM(IPZ,NM) = USUM(IPZ,NM) + WALL(IW)%ONE_D%UWS*WALL(IW)%AW
   ENDDO WALL_LOOP4
   
ENDDO PRESSURE_ZONE_LOOP

TUSED(2,NM)=TUSED(2,NM)+SECOND()-TNOW

CONTAINS


SUBROUTINE ENTHALPY_ADVECTION

REAL(EB), POINTER, DIMENSION(:,:,:) :: HSX=>NULL(),HSY=>NULL(),HSZ=>NULL(),DV=>NULL()
REAL(EB) :: DR,B

HSX=>WORK2!; HSX=1.E20_EB
HSY=>WORK3!; HSY=1.E20_EB
HSZ=>WORK4!; HSZ=1.E20_EB
U_DOT_DEL_RHO_H_S=>WORK6; U_DOT_DEL_RHO_H_S=0._EB

IF (.NOT.ENTHALPY_TRANSPORT) RETURN

LIMITER_SELECT: SELECT CASE (FLUX_LIMITER)

   CASE (SUPERBEE_LIMITER) LIMITER_SELECT

      DV=>WORK7

      ! compute data variation and face value x

      !DV = 1.E20_EB
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=0,IBAR
               DV(I,J,K) = RHO_H_S_P(I+1,J,K) - RHO_H_S_P(I,J,K)
            ENDDO
         ENDDO
      ENDDO

      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBM1

               IF (ABS(DV(I,J,K))>ADVECTION_EPS) THEN
                  IF (UU(I,J,K)>0._EB) THEN
                     DR = DV(I-1,J,K)/DV(I,J,K)
                     B = MAX(0._EB,MIN(2._EB*DR,1._EB),MIN(DR,2._EB))
                     HSX(I,J,K) = RHO_H_S_P(I,J,K)   + 0.5_EB*B*DV(I,J,K)
                  ELSE
                     DR = DV(I+1,J,K)/DV(I,J,K)
                     B = MAX(0._EB,MIN(2._EB*DR,1._EB),MIN(DR,2._EB))
                     HSX(I,J,K) = RHO_H_S_P(I+1,J,K) - 0.5_EB*B*DV(I,J,K)
                  ENDIF
               ELSE
                  HSX(I,J,K) = 0.5_EB*(RHO_H_S_P(I,J,K) + RHO_H_S_P(I+1,J,K))
               ENDIF

            ENDDO
         ENDDO
      ENDDO

      ! compute data variation and face value in y

      !DV = 1.E20_EB
      DO K=1,KBAR
         DO J=0,JBAR
            DO I=1,IBAR
               DV(I,J,K) = RHO_H_S_P(I,J+1,K) - RHO_H_S_P(I,J,K)
            ENDDO
         ENDDO
      ENDDO

      DO K=1,KBAR
         DO J=1,JBM1
            DO I=1,IBAR

               IF (ABS(DV(I,J,K))>ADVECTION_EPS) THEN
                  IF (VV(I,J,K)>0._EB) THEN
                     DR = DV(I,J-1,K)/DV(I,J,K)
                     B = MAX(0._EB,MIN(2._EB*DR,1._EB),MIN(DR,2._EB))
                     HSY(I,J,K) = RHO_H_S_P(I,J,K)   + 0.5_EB*B*DV(I,J,K)
                  ELSE
                     DR = DV(I,J+1,K)/DV(I,J,K)
                     B = MAX(0._EB,MIN(2._EB*DR,1._EB),MIN(DR,2._EB))
                     HSY(I,J,K) = RHO_H_S_P(I,J+1,K) - 0.5_EB*B*DV(I,J,K)
                  ENDIF
               ELSE
                  HSY(I,J,K) = 0.5_EB*(RHO_H_S_P(I,J,K) + RHO_H_S_P(I,J+1,K))
               ENDIF

            ENDDO
         ENDDO
      ENDDO

      ! compute data variation and face value in z

      !DV = 1.E20_EB
      DO K=0,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               DV(I,J,K) = RHO_H_S_P(I,J,K+1) - RHO_H_S_P(I,J,K)
            ENDDO
         ENDDO
      ENDDO

      DO K=1,KBM1   
         DO J=1,JBAR
            DO I=1,IBAR

               IF (ABS(DV(I,J,K))>ADVECTION_EPS) THEN
                  IF (WW(I,J,K)>0._EB) THEN
                     DR = DV(I,J,K-1)/DV(I,J,K)
                     B = MAX(0._EB,MIN(2._EB*DR,1._EB),MIN(DR,2._EB))
                     HSZ(I,J,K) = RHO_H_S_P(I,J,K)   + 0.5_EB*B*DV(I,J,K)
                  ELSE
                     DR = DV(I,J,K+1)/DV(I,J,K)
                     B = MAX(0._EB,MIN(2._EB*DR,1._EB),MIN(DR,2._EB))
                     HSZ(I,J,K) = RHO_H_S_P(I,J,K+1) - 0.5_EB*B*DV(I,J,K)
                  ENDIF
               ELSE
                  HSZ(I,J,K) = 0.5_EB*(RHO_H_S_P(I,J,K) + RHO_H_S_P(I,J,K+1))
               ENDIF

            ENDDO
         ENDDO
      ENDDO

   CASE (CENTRAL_LIMITER) LIMITER_SELECT

      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBM1
               HSX(I,J,K) = 0.5_EB*(RHO_H_S_P(I,J,K) + RHO_H_S_P(I+1,J,K))
            ENDDO
         ENDDO
      ENDDO

      DO K=1,KBAR
         DO J=1,JBM1
            DO I=1,IBAR
               HSY(I,J,K) = 0.5_EB*(RHO_H_S_P(I,J,K) + RHO_H_S_P(I,J+1,K))
            ENDDO
         ENDDO
      ENDDO

      DO K=1,KBM1
         DO J=1,JBAR
            DO I=1,IBAR
               HSZ(I,J,K) = 0.5_EB*(RHO_H_S_P(I,J,K) + RHO_H_S_P(I,J,K+1))
            ENDDO
         ENDDO
      ENDDO

   CASE DEFAULT LIMITER_SELECT

      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBM1
               ZZZ(1:4) = RHO_H_S_P(I-1:I+2,J,K)
               HSX(I,J,K) = SCALAR_FACE_VALUE(UU(I,J,K),ZZZ,FLUX_LIMITER)
            ENDDO
         ENDDO
      ENDDO

      DO K=1,KBAR
         DO J=1,JBM1
            DO I=1,IBAR
               ZZZ(1:4) = RHO_H_S_P(I,J-1:J+2,K)
               HSY(I,J,K) = SCALAR_FACE_VALUE(VV(I,J,K),ZZZ,FLUX_LIMITER)
            ENDDO
         ENDDO
      ENDDO

      DO K=1,KBM1 
         DO J=1,JBAR
            DO I=1,IBAR
               ZZZ(1:4) = RHO_H_S_P(I,J,K-1:K+2)
               HSZ(I,J,K) = SCALAR_FACE_VALUE(WW(I,J,K),ZZZ,FLUX_LIMITER)
            ENDDO
         ENDDO
      ENDDO

END SELECT LIMITER_SELECT

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
   IF (N_TRACKED_SPECIES>0) ZZ_GET(1:N_TRACKED_SPECIES) = WC%ZZ_F(1:N_TRACKED_SPECIES)
   CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,WC%ONE_D%TMP_F)

   ! overwrite first off-wall advective flux if flow is away from the wall and if the face is not also a wall cell

   IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY .AND. WC%BOUNDARY_TYPE/=OPEN_BOUNDARY) THEN

      OFF_WALL_SELECT_1: SELECT CASE(IOR)
         CASE( 1) OFF_WALL_SELECT_1
            !      ghost          FX/UU(II+1)
            ! ///   II   ///  II+1  |  II+2  | ...
            !                       ^ WALL_INDEX(II+1,+1)
            IF ((UU(II+1,JJ,KK)>0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II+1,JJ,KK),+1)>0)) THEN
               ZZZ(1:3) = (/WC%RHO_F*H_S,RHO_H_S_P(II+1:II+2,JJ,KK)/)
               HSX(II+1,JJ,KK) = SCALAR_FACE_VALUE(UU(II+1,JJ,KK),ZZZ,FLUX_LIMITER)
            ENDIF
         CASE(-1) OFF_WALL_SELECT_1
            !            FX/UU(II-2)     ghost
            ! ... |  II-2  |  II-1  ///   II   ///
            !              ^ WALL_INDEX(II-1,-1)
            IF ((UU(II-2,JJ,KK)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II-1,JJ,KK),-1)>0)) THEN
               ZZZ(2:4) = (/RHO_H_S_P(II-2:II-1,JJ,KK),WC%RHO_F*H_S/)
               HSX(II-2,JJ,KK) = SCALAR_FACE_VALUE(UU(II-2,JJ,KK),ZZZ,FLUX_LIMITER)
            ENDIF
         CASE( 2) OFF_WALL_SELECT_1
            IF ((VV(II,JJ+1,KK)>0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ+1,KK),+2)>0)) THEN
               ZZZ(1:3) = (/WC%RHO_F*H_S,RHO_H_S_P(II,JJ+1:JJ+2,KK)/)
               HSY(II,JJ+1,KK) = SCALAR_FACE_VALUE(VV(II,JJ+1,KK),ZZZ,FLUX_LIMITER)
            ENDIF
         CASE(-2) OFF_WALL_SELECT_1
            IF ((VV(II,JJ-2,KK)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ-1,KK),-2)>0)) THEN
               ZZZ(2:4) = (/RHO_H_S_P(II,JJ-2:JJ-1,KK),WC%RHO_F*H_S/)
               HSY(II,JJ-2,KK) = SCALAR_FACE_VALUE(VV(II,JJ-2,KK),ZZZ,FLUX_LIMITER)
            ENDIF
         CASE( 3) OFF_WALL_SELECT_1
            IF ((WW(II,JJ,KK+1)>0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ,KK+1),+3)>0)) THEN
               ZZZ(1:3) = (/WC%RHO_F*H_S,RHO_H_S_P(II,JJ,KK+1:KK+2)/)
               HSZ(II,JJ,KK+1) = SCALAR_FACE_VALUE(WW(II,JJ,KK+1),ZZZ,FLUX_LIMITER)
            ENDIF
         CASE(-3) OFF_WALL_SELECT_1
            IF ((WW(II,JJ,KK-2)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ,KK-1),-3)>0)) THEN
               ZZZ(2:4) = (/RHO_H_S_P(II,JJ,KK-2:KK-1),WC%RHO_F*H_S/)
               HSZ(II,JJ,KK-2) = SCALAR_FACE_VALUE(WW(II,JJ,KK-2),ZZZ,FLUX_LIMITER)
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
      CASE(SOLID_BOUNDARY,HVAC_BOUNDARY)
         IF (PREDICTOR) UN = -SIGN(1._EB,REAL(IOR,EB))*WC%ONE_D%UWS
         IF (CORRECTOR) UN = -SIGN(1._EB,REAL(IOR,EB))*WC%ONE_D%UW
      CASE(INTERPOLATED_BOUNDARY)
         UN = UVW_SAVE(IW)
   END SELECT BOUNDARY_SELECT

   SELECT CASE(IOR)
      CASE( 1)
         HSX(II,JJ,KK)   = RHO_H_S_P(IIG,JJG,KKG) ! zero out DU at wall
         DU_M = (WC%RHO_F*H_S - RHO_H_S_P(IIG,JJG,KKG))*UN
         U_DOT_DEL_RHO_H_S(IIG,JJG,KKG) = U_DOT_DEL_RHO_H_S(IIG,JJG,KKG) - DU_M*WC%RDN
      CASE(-1)
         HSX(II-1,JJ,KK) = RHO_H_S_P(IIG,JJG,KKG)
         DU_P = (WC%RHO_F*H_S - RHO_H_S_P(IIG,JJG,KKG))*UN
         U_DOT_DEL_RHO_H_S(IIG,JJG,KKG) = U_DOT_DEL_RHO_H_S(IIG,JJG,KKG) + DU_P*WC%RDN
      CASE( 2)
         HSY(II,JJ,KK)   = RHO_H_S_P(IIG,JJG,KKG)
         DU_M = (WC%RHO_F*H_S - RHO_H_S_P(IIG,JJG,KKG))*UN
         U_DOT_DEL_RHO_H_S(IIG,JJG,KKG) = U_DOT_DEL_RHO_H_S(IIG,JJG,KKG) - DU_M*WC%RDN
      CASE(-2)
         HSY(II,JJ-1,KK) = RHO_H_S_P(IIG,JJG,KKG)
         DU_P = (WC%RHO_F*H_S - RHO_H_S_P(IIG,JJG,KKG))*UN
         U_DOT_DEL_RHO_H_S(IIG,JJG,KKG) = U_DOT_DEL_RHO_H_S(IIG,JJG,KKG) + DU_P*WC%RDN
      CASE( 3)
         HSZ(II,JJ,KK)   = RHO_H_S_P(IIG,JJG,KKG)
         DU_M = (WC%RHO_F*H_S - RHO_H_S_P(IIG,JJG,KKG))*UN
         U_DOT_DEL_RHO_H_S(IIG,JJG,KKG) = U_DOT_DEL_RHO_H_S(IIG,JJG,KKG) - DU_M*WC%RDN
      CASE(-3)
         HSZ(II,JJ,KK-1) = RHO_H_S_P(IIG,JJG,KKG)
         DU_P = (WC%RHO_F*H_S - RHO_H_S_P(IIG,JJG,KKG))*UN
         U_DOT_DEL_RHO_H_S(IIG,JJG,KKG) = U_DOT_DEL_RHO_H_S(IIG,JJG,KKG) + DU_P*WC%RDN
   END SELECT

ENDDO WALL_LOOP

DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IF (SOLID(CELL_INDEX(I,J,K))) CYCLE

         DU_P = (HSX(I,J,K)   - RHO_H_S_P(I,J,K))*UU(I,J,K)                       ! FDS Tech Guide (B.13)
         DU_M = (HSX(I-1,J,K) - RHO_H_S_P(I,J,K))*UU(I-1,J,K)                     ! FDS Tech Guide (B.14)
         U_DOT_DEL_RHO_H_S(I,J,K) = U_DOT_DEL_RHO_H_S(I,J,K) + (DU_P-DU_M)*RDX(I) ! FDS Tech Guide (B.12)

         DU_P = (HSY(I,J,K)   - RHO_H_S_P(I,J,K))*VV(I,J,K)
         DU_M = (HSY(I,J-1,K) - RHO_H_S_P(I,J,K))*VV(I,J-1,K)
         U_DOT_DEL_RHO_H_S(I,J,K) = U_DOT_DEL_RHO_H_S(I,J,K) + (DU_P-DU_M)*RDY(J)

         DU_P = (HSZ(I,J,K)   - RHO_H_S_P(I,J,K))*WW(I,J,K)
         DU_M = (HSZ(I,J,K-1) - RHO_H_S_P(I,J,K))*WW(I,J,K-1)
         U_DOT_DEL_RHO_H_S(I,J,K) = U_DOT_DEL_RHO_H_S(I,J,K) + (DU_P-DU_M)*RDZ(K)

      ENDDO
   ENDDO 
ENDDO

END SUBROUTINE ENTHALPY_ADVECTION


SUBROUTINE DENSITY_ADVECTION

USE MANUFACTURED_SOLUTIONS, ONLY: VD2D_MMS_RHO

REAL(EB), POINTER, DIMENSION(:,:,:) :: DV=>NULL()
REAL(EB) :: DR,B

!FX(:,:,:,0)=1.E20_EB
!FY(:,:,:,0)=1.E20_EB
!FZ(:,:,:,0)=1.E20_EB

IF (.NOT.CONSTANT_SPECIFIC_HEAT_RATIO) THEN
   U_DOT_DEL_RHO=>WORK8 
   U_DOT_DEL_RHO=0._EB
ENDIF

LIMITER_SELECT: SELECT CASE (FLUX_LIMITER)

   CASE (SUPERBEE_LIMITER) LIMITER_SELECT

      DV=>WORK2

      ! compute data variation and face value x

      !DV = 1.E20_EB
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=0,IBAR
               DV(I,J,K) = RHOP(I+1,J,K) - RHOP(I,J,K)
            ENDDO
         ENDDO
      ENDDO

      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBM1

               IF (ABS(DV(I,J,K))>ADVECTION_EPS) THEN
                  IF (UU(I,J,K)>0._EB) THEN
                     DR = DV(I-1,J,K)/DV(I,J,K)
                     B = MAX(0._EB,MIN(2._EB*DR,1._EB),MIN(DR,2._EB))
                     FX(I,J,K,0) = RHOP(I,J,K)   + 0.5_EB*B*DV(I,J,K)
                  ELSE
                     DR = DV(I+1,J,K)/DV(I,J,K)
                     B = MAX(0._EB,MIN(2._EB*DR,1._EB),MIN(DR,2._EB))
                     FX(I,J,K,0) = RHOP(I+1,J,K) - 0.5_EB*B*DV(I,J,K)
                  ENDIF
               ELSE
                  FX(I,J,K,0) = 0.5_EB*(RHOP(I,J,K) + RHOP(I+1,J,K))
               ENDIF

            ENDDO
         ENDDO
      ENDDO

      ! compute data variation and face value y

      !DV = 1.E20_EB
      DO K=1,KBAR
         DO J=0,JBAR
            DO I=1,IBAR
               DV(I,J,K) = RHOP(I,J+1,K) - RHOP(I,J,K)
            ENDDO
         ENDDO
      ENDDO

      DO K=1,KBAR
         DO J=1,JBM1
            DO I=1,IBAR

               IF (ABS(DV(I,J,K))>ADVECTION_EPS) THEN
                  IF (VV(I,J,K)>0._EB) THEN
                     DR = DV(I,J-1,K)/DV(I,J,K)
                     B = MAX(0._EB,MIN(2._EB*DR,1._EB),MIN(DR,2._EB))
                     FY(I,J,K,0) = RHOP(I,J,K)   + 0.5_EB*B*DV(I,J,K)
                  ELSE
                     DR = DV(I,J+1,K)/DV(I,J,K)
                     B = MAX(0._EB,MIN(2._EB*DR,1._EB),MIN(DR,2._EB))
                     FY(I,J,K,0) = RHOP(I,J+1,K) - 0.5_EB*B*DV(I,J,K)
                  ENDIF
               ELSE
                  FY(I,J,K,0) = 0.5_EB*(RHOP(I,J,K) + RHOP(I,J+1,K))
               ENDIF

            ENDDO
         ENDDO
      ENDDO

      ! compute data variation and face value z

      !DV = 1.E20_EB
      DO K=0,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               DV(I,J,K) = RHOP(I,J,K+1) - RHOP(I,J,K)
            ENDDO
         ENDDO
      ENDDO

      DO K=1,KBM1   
         DO J=1,JBAR
            DO I=1,IBAR

               IF (ABS(DV(I,J,K))>ADVECTION_EPS) THEN
                  IF (WW(I,J,K)>0._EB) THEN
                     DR = DV(I,J,K-1)/DV(I,J,K)
                     B = MAX(0._EB,MIN(2._EB*DR,1._EB),MIN(DR,2._EB))
                     FZ(I,J,K,0) = RHOP(I,J,K)   + 0.5_EB*B*DV(I,J,K)
                  ELSE
                     DR = DV(I,J,K+1)/DV(I,J,K)
                     B = MAX(0._EB,MIN(2._EB*DR,1._EB),MIN(DR,2._EB))
                     FZ(I,J,K,0) = RHOP(I,J,K+1) - 0.5_EB*B*DV(I,J,K)
                  ENDIF
               ELSE
                  FZ(I,J,K,0) = 0.5_EB*(RHOP(I,J,K) + RHOP(I,J,K+1))
               ENDIF

            ENDDO
         ENDDO
      ENDDO

   CASE (CENTRAL_LIMITER) LIMITER_SELECT

      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBM1
               FX(I,J,K,0) = 0.5_EB*(RHOP(I,J,K) + RHOP(I+1,J,K))
            ENDDO
         ENDDO
      ENDDO

      DO K=1,KBAR
         DO J=1,JBM1
            DO I=1,IBAR
               FY(I,J,K,0) = 0.5_EB*(RHOP(I,J,K) + RHOP(I,J+1,K))
            ENDDO
         ENDDO
      ENDDO

      DO K=1,KBM1
         DO J=1,JBAR
            DO I=1,IBAR
               FZ(I,J,K,0) = 0.5_EB*(RHOP(I,J,K) + RHOP(I,J,K+1))
            ENDDO
         ENDDO
      ENDDO

   CASE DEFAULT LIMITER_SELECT

      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBM1
               ZZZ(1:4) = RHOP(I-1:I+2,J,K)
               FX(I,J,K,0) = SCALAR_FACE_VALUE(UU(I,J,K),ZZZ,FLUX_LIMITER)
            ENDDO
         ENDDO
      ENDDO

      DO K=1,KBAR
         DO J=1,JBM1
            DO I=1,IBAR
               ZZZ(1:4) = RHOP(I,J-1:J+2,K)
               FY(I,J,K,0) = SCALAR_FACE_VALUE(VV(I,J,K),ZZZ,FLUX_LIMITER)
            ENDDO
         ENDDO
      ENDDO

      DO K=1,KBM1
         DO J=1,JBAR
            DO I=1,IBAR
               ZZZ(1:4) = RHOP(I,J,K-1:K+2)
               FZ(I,J,K,0) = SCALAR_FACE_VALUE(WW(I,J,K),ZZZ,FLUX_LIMITER)
            ENDDO
         ENDDO
      ENDDO

END SELECT LIMITER_SELECT

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

   ! overwrite first off-wall advective flux if flow is away from the wall and if the face is not also a wall cell

   IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY .AND. WC%BOUNDARY_TYPE/=OPEN_BOUNDARY) THEN

      OFF_WALL_SELECT_2: SELECT CASE(IOR)
         CASE( 1) OFF_WALL_SELECT_2
            !      ghost          FX/UU(II+1)
            ! ///   II   ///  II+1  |  II+2  | ...
            !                       ^ WALL_INDEX(II+1,+1)
            IF ((UU(II+1,JJ,KK)>0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II+1,JJ,KK),+1)>0)) THEN
               ZZZ(1:3) = (/WC%RHO_F,RHOP(II+1:II+2,JJ,KK)/)
               FX(II+1,JJ,KK,0) = SCALAR_FACE_VALUE(UU(II+1,JJ,KK),ZZZ,FLUX_LIMITER)
            ENDIF
         CASE(-1) OFF_WALL_SELECT_2
            !            FX/UU(II-2)     ghost
            ! ... |  II-2  |  II-1  ///   II   ///
            !              ^ WALL_INDEX(II-1,-1)
            IF ((UU(II-2,JJ,KK)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II-1,JJ,KK),-1)>0)) THEN
               ZZZ(2:4) = (/RHOP(II-2:II-1,JJ,KK),WC%RHO_F/)
               FX(II-2,JJ,KK,0) = SCALAR_FACE_VALUE(UU(II-2,JJ,KK),ZZZ,FLUX_LIMITER)
            ENDIF
         CASE( 2) OFF_WALL_SELECT_2
            IF ((VV(II,JJ+1,KK)>0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ+1,KK),+2)>0)) THEN
               ZZZ(1:3) = (/WC%RHO_F,RHOP(II,JJ+1:JJ+2,KK)/)
               FY(II,JJ+1,KK,0) = SCALAR_FACE_VALUE(VV(II,JJ+1,KK),ZZZ,FLUX_LIMITER)
            ENDIF
         CASE(-2) OFF_WALL_SELECT_2
            IF ((VV(II,JJ-2,KK)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ-1,KK),-2)>0)) THEN
               ZZZ(2:4) = (/RHOP(II,JJ-2:JJ-1,KK),WC%RHO_F/)
               FY(II,JJ-2,KK,0) = SCALAR_FACE_VALUE(VV(II,JJ-2,KK),ZZZ,FLUX_LIMITER)
            ENDIF
         CASE( 3) OFF_WALL_SELECT_2
            IF ((WW(II,JJ,KK+1)>0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ,KK+1),+3)>0)) THEN
               ZZZ(1:3) = (/WC%RHO_F,RHOP(II,JJ,KK+1:KK+2)/)
               FZ(II,JJ,KK+1,0) = SCALAR_FACE_VALUE(WW(II,JJ,KK+1),ZZZ,FLUX_LIMITER)
            ENDIF
         CASE(-3) OFF_WALL_SELECT_2
            IF ((WW(II,JJ,KK-2)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ,KK-1),-3)>0)) THEN
               ZZZ(2:4) = (/RHOP(II,JJ,KK-2:KK-1),WC%RHO_F/)
               FZ(II,JJ,KK-2,0) = SCALAR_FACE_VALUE(WW(II,JJ,KK-2),ZZZ,FLUX_LIMITER)
            ENDIF
      END SELECT OFF_WALL_SELECT_2
   
   ENDIF

   SELECT CASE(IOR)
      CASE( 1)
         FX(II,JJ,KK,0)   = RHOP(IIG,JJG,KKG) ! zero out DU at wall
      CASE(-1)
         FX(II-1,JJ,KK,0) = RHOP(IIG,JJG,KKG)
      CASE( 2)
         FY(II,JJ,KK,0)   = RHOP(IIG,JJG,KKG)
      CASE(-2)
         FY(II,JJ-1,KK,0) = RHOP(IIG,JJG,KKG)
      CASE( 3)
         FZ(II,JJ,KK,0)   = RHOP(IIG,JJG,KKG)
      CASE(-3)
         FZ(II,JJ,KK-1,0) = RHOP(IIG,JJG,KKG)
   END SELECT

   ! Correct U_DOT_DEL_RHO at the boundaries

   IF (.NOT.CONSTANT_SPECIFIC_HEAT_RATIO) THEN

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
         CASE(SOLID_BOUNDARY,HVAC_BOUNDARY)
            IF (PREDICTOR) UN = -SIGN(1._EB,REAL(IOR,EB))*WC%ONE_D%UWS
            IF (CORRECTOR) UN = -SIGN(1._EB,REAL(IOR,EB))*WC%ONE_D%UW
         CASE(INTERPOLATED_BOUNDARY)
            UN = UVW_SAVE(IW)
      END SELECT BOUNDARY_SELECT
   
      SELECT CASE(IOR)
         CASE( 1)
            DU_M = (WC%RHO_F - RHOP(IIG,JJG,KKG))*UN
            U_DOT_DEL_RHO(IIG,JJG,KKG) = U_DOT_DEL_RHO(IIG,JJG,KKG) - DU_M*WC%RDN
         CASE(-1)
            DU_P = (WC%RHO_F - RHOP(IIG,JJG,KKG))*UN
            U_DOT_DEL_RHO(IIG,JJG,KKG) = U_DOT_DEL_RHO(IIG,JJG,KKG) + DU_P*WC%RDN
         CASE( 2)
            DU_M = (WC%RHO_F - RHOP(IIG,JJG,KKG))*UN
            U_DOT_DEL_RHO(IIG,JJG,KKG) = U_DOT_DEL_RHO(IIG,JJG,KKG) - DU_M*WC%RDN
         CASE(-2)
            DU_P = (WC%RHO_F - RHOP(IIG,JJG,KKG))*UN
            U_DOT_DEL_RHO(IIG,JJG,KKG) = U_DOT_DEL_RHO(IIG,JJG,KKG) + DU_P*WC%RDN
         CASE( 3)
            DU_M = (WC%RHO_F - RHOP(IIG,JJG,KKG))*UN
            U_DOT_DEL_RHO(IIG,JJG,KKG) = U_DOT_DEL_RHO(IIG,JJG,KKG) - DU_M*WC%RDN
         CASE(-3)
            DU_P = (WC%RHO_F - RHOP(IIG,JJG,KKG))*UN
            U_DOT_DEL_RHO(IIG,JJG,KKG) = U_DOT_DEL_RHO(IIG,JJG,KKG) + DU_P*WC%RDN
      END SELECT

   ENDIF
      
ENDDO WALL_LOOP

! Compute U_DOT_DEL_RHO at internal cells

IF (.NOT.CONSTANT_SPECIFIC_HEAT_RATIO) THEN

   IF (.NOT.ENTHALPY_TRANSPORT) THEN
      U_DOT_DEL_RHO = 0._EB
      RETURN
   ENDIF

   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
   
            DU_P = (FX(I,J,K,0)   - RHOP(I,J,K))*UU(I,J,K)
            DU_M = (FX(I-1,J,K,0) - RHOP(I,J,K))*UU(I-1,J,K)
            U_DOT_DEL_RHO(I,J,K) = U_DOT_DEL_RHO(I,J,K) + (DU_P-DU_M)*RDX(I)
   
            DU_P = (FY(I,J,K,0)   - RHOP(I,J,K))*VV(I,J,K)
            DU_M = (FY(I,J-1,K,0) - RHOP(I,J,K))*VV(I,J-1,K)
            U_DOT_DEL_RHO(I,J,K) = U_DOT_DEL_RHO(I,J,K) + (DU_P-DU_M)*RDY(J)
   
            DU_P = (FZ(I,J,K,0)   - RHOP(I,J,K))*WW(I,J,K)
            DU_M = (FZ(I,J,K-1,0) - RHOP(I,J,K))*WW(I,J,K-1)
            U_DOT_DEL_RHO(I,J,K) = U_DOT_DEL_RHO(I,J,K) + (DU_P-DU_M)*RDZ(K)
   
         ENDDO
      ENDDO 
   ENDDO

ENDIF

END SUBROUTINE DENSITY_ADVECTION


SUBROUTINE SPECIES_ADVECTION

REAL(EB), POINTER, DIMENSION(:,:,:) :: DV=>NULL()
REAL(EB) :: DR,B

RHO_Z_P=>WORK6
!RHO_Z_P=1.E20_EB

!$OMP PARALLEL DO SCHEDULE(static)
DO K=0,KBP1
   DO J=0,JBP1
      DO I=0,IBP1
         RHO_Z_P(I,J,K) = RHOP(I,J,K)*ZZP(I,J,K,N)
      ENDDO
   ENDDO
ENDDO
!$OMP END PARALLEL DO

!FX(:,:,:,N)=1.E20_EB
!FY(:,:,:,N)=1.E20_EB
!FZ(:,:,:,N)=1.E20_EB

IF (.NOT.CONSTANT_SPECIFIC_HEAT_RATIO) THEN
   U_DOT_DEL_RHO_Z=>WORK7 
   U_DOT_DEL_RHO_Z=0._EB
ENDIF

LIMITER_SELECT: SELECT CASE (FLUX_LIMITER)

   CASE (SUPERBEE_LIMITER) LIMITER_SELECT

      DV=>WORK2

      ! compute data variation and face value x

      !DV = 1.E20_EB
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=0,IBAR
               DV(I,J,K) = RHO_Z_P(I+1,J,K) - RHO_Z_P(I,J,K)
            ENDDO
         ENDDO
      ENDDO

      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBM1

               IF (ABS(DV(I,J,K))>ADVECTION_EPS) THEN
                  IF (UU(I,J,K)>0._EB) THEN
                     DR = DV(I-1,J,K)/DV(I,J,K)
                     B = MAX(0._EB,MIN(2._EB*DR,1._EB),MIN(DR,2._EB))
                     FX(I,J,K,N) = RHO_Z_P(I,J,K)   + 0.5_EB*B*DV(I,J,K)
                  ELSE
                     DR = DV(I+1,J,K)/DV(I,J,K)
                     B = MAX(0._EB,MIN(2._EB*DR,1._EB),MIN(DR,2._EB))
                     FX(I,J,K,N) = RHO_Z_P(I+1,J,K) - 0.5_EB*B*DV(I,J,K)
                  ENDIF
               ELSE
                  FX(I,J,K,N) = 0.5_EB*(RHO_Z_P(I,J,K) + RHO_Z_P(I+1,J,K))
               ENDIF

            ENDDO
         ENDDO
      ENDDO

      ! compute data variation and face value y

      !DV = 1.E20_EB
      DO K=1,KBAR
         DO J=0,JBAR
            DO I=1,IBAR
               DV(I,J,K) = RHO_Z_P(I,J+1,K) - RHO_Z_P(I,J,K)
            ENDDO
         ENDDO
      ENDDO

      DO K=1,KBAR
         DO J=1,JBM1
            DO I=1,IBAR

               IF (ABS(DV(I,J,K))>ADVECTION_EPS) THEN
                  IF (VV(I,J,K)>0._EB) THEN
                     DR = DV(I,J-1,K)/DV(I,J,K)
                     B = MAX(0._EB,MIN(2._EB*DR,1._EB),MIN(DR,2._EB))
                     FY(I,J,K,N) = RHO_Z_P(I,J,K)   + 0.5_EB*B*DV(I,J,K)
                  ELSE
                     DR = DV(I,J+1,K)/DV(I,J,K)
                     B = MAX(0._EB,MIN(2._EB*DR,1._EB),MIN(DR,2._EB))
                     FY(I,J,K,N) = RHO_Z_P(I,J+1,K) - 0.5_EB*B*DV(I,J,K)
                  ENDIF
               ELSE
                  FY(I,J,K,N) = 0.5_EB*(RHO_Z_P(I,J,K) + RHO_Z_P(I,J+1,K))
               ENDIF

            ENDDO
         ENDDO
      ENDDO

      ! compute data variation and face value z

      !DV = 1.E20_EB
      DO K=0,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               DV(I,J,K) = RHO_Z_P(I,J,K+1) - RHO_Z_P(I,J,K)
            ENDDO
         ENDDO
      ENDDO

      DO K=1,KBM1   
         DO J=1,JBAR
            DO I=1,IBAR

               IF (ABS(DV(I,J,K))>ADVECTION_EPS) THEN
                  IF (WW(I,J,K)>0._EB) THEN
                     DR = DV(I,J,K-1)/DV(I,J,K)
                     B = MAX(0._EB,MIN(2._EB*DR,1._EB),MIN(DR,2._EB))
                     FZ(I,J,K,N) = RHO_Z_P(I,J,K)   + 0.5_EB*B*DV(I,J,K)
                  ELSE
                     DR = DV(I,J,K+1)/DV(I,J,K)
                     B = MAX(0._EB,MIN(2._EB*DR,1._EB),MIN(DR,2._EB))
                     FZ(I,J,K,N) = RHO_Z_P(I,J,K+1) - 0.5_EB*B*DV(I,J,K)
                  ENDIF
               ELSE
                  FZ(I,J,K,N) = 0.5_EB*(RHO_Z_P(I,J,K) + RHO_Z_P(I,J,K+1))
               ENDIF

            ENDDO
         ENDDO
      ENDDO

   CASE (CENTRAL_LIMITER) LIMITER_SELECT

      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBM1
               FX(I,J,K,N) = 0.5_EB*(RHO_Z_P(I,J,K) + RHO_Z_P(I+1,J,K))
            ENDDO
         ENDDO
      ENDDO

      DO K=1,KBAR
         DO J=1,JBM1
            DO I=1,IBAR
               FY(I,J,K,N) = 0.5_EB*(RHO_Z_P(I,J,K) + RHO_Z_P(I,J+1,K))
            ENDDO
         ENDDO
      ENDDO

      DO K=1,KBM1
         DO J=1,JBAR
            DO I=1,IBAR
               FZ(I,J,K,N) = 0.5_EB*(RHO_Z_P(I,J,K) + RHO_Z_P(I,J,K+1))
            ENDDO
         ENDDO
      ENDDO

   CASE DEFAULT LIMITER_SELECT

      !$OMP PARALLEL PRIVATE(ZZZ)
      !$OMP DO SCHEDULE(static)
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBM1
               ZZZ(1:4) = RHO_Z_P(I-1:I+2,J,K)
               FX(I,J,K,N) = SCALAR_FACE_VALUE(UU(I,J,K),ZZZ,FLUX_LIMITER)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO NOWAIT

      !$OMP DO SCHEDULE(static)
      DO K=1,KBAR
         DO J=1,JBM1
            DO I=1,IBAR
               ZZZ(1:4) = RHO_Z_P(I,J-1:J+2,K)
               FY(I,J,K,N) = SCALAR_FACE_VALUE(VV(I,J,K),ZZZ,FLUX_LIMITER)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO NOWAIT

      !$OMP DO SCHEDULE(static)
      DO K=1,KBM1
         DO J=1,JBAR
            DO I=1,IBAR
               ZZZ(1:4) = RHO_Z_P(I,J,K-1:K+2)
               FZ(I,J,K,N) = SCALAR_FACE_VALUE(WW(I,J,K),ZZZ,FLUX_LIMITER)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO NOWAIT
      !$OMP END PARALLEL

END SELECT LIMITER_SELECT

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

   ! overwrite first off-wall advective flux if flow is away from the wall and if the face is not also a wall cell

   IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY .AND. WC%BOUNDARY_TYPE/=OPEN_BOUNDARY) THEN

      OFF_WALL_SELECT_3: SELECT CASE(IOR)
         CASE( 1) OFF_WALL_SELECT_3
            !      ghost          FX/UU(II+1)
            ! ///   II   ///  II+1  |  II+2  | ...
            !                       ^ WALL_INDEX(II+1,+1)
            IF ((UU(II+1,JJ,KK)>0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II+1,JJ,KK),+1)>0)) THEN
               ZZZ(1:3) = (/WC%RHO_F*WC%ZZ_F(N),RHO_Z_P(II+1:II+2,JJ,KK)/)
               FX(II+1,JJ,KK,N) = SCALAR_FACE_VALUE(UU(II+1,JJ,KK),ZZZ,FLUX_LIMITER)
            ENDIF
         CASE(-1) OFF_WALL_SELECT_3
            !            FX/UU(II-2)     ghost
            ! ... |  II-2  |  II-1  ///   II   ///
            !              ^ WALL_INDEX(II-1,-1)
            IF ((UU(II-2,JJ,KK)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II-1,JJ,KK),-1)>0)) THEN
               ZZZ(2:4) = (/RHO_Z_P(II-2:II-1,JJ,KK),WC%RHO_F*WC%ZZ_F(N)/)
               FX(II-2,JJ,KK,N) = SCALAR_FACE_VALUE(UU(II-2,JJ,KK),ZZZ,FLUX_LIMITER)
            ENDIF
         CASE( 2) OFF_WALL_SELECT_3
            IF ((VV(II,JJ+1,KK)>0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ+1,KK),+2)>0)) THEN
               ZZZ(1:3) = (/WC%RHO_F*WC%ZZ_F(N),RHO_Z_P(II,JJ+1:JJ+2,KK)/)
               FY(II,JJ+1,KK,N) = SCALAR_FACE_VALUE(VV(II,JJ+1,KK),ZZZ,FLUX_LIMITER)
            ENDIF
         CASE(-2) OFF_WALL_SELECT_3
            IF ((VV(II,JJ-2,KK)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ-1,KK),-2)>0)) THEN
               ZZZ(2:4) = (/RHO_Z_P(II,JJ-2:JJ-1,KK),WC%RHO_F*WC%ZZ_F(N)/)
               FY(II,JJ-2,KK,N) = SCALAR_FACE_VALUE(VV(II,JJ-2,KK),ZZZ,FLUX_LIMITER)
            ENDIF
         CASE( 3) OFF_WALL_SELECT_3
            IF ((WW(II,JJ,KK+1)>0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ,KK+1),+3)>0)) THEN
               ZZZ(1:3) = (/WC%RHO_F*WC%ZZ_F(N),RHO_Z_P(II,JJ,KK+1:KK+2)/)
               FZ(II,JJ,KK+1,N) = SCALAR_FACE_VALUE(WW(II,JJ,KK+1),ZZZ,FLUX_LIMITER)
            ENDIF
         CASE(-3) OFF_WALL_SELECT_3
            IF ((WW(II,JJ,KK-2)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ,KK-1),-3)>0)) THEN
               ZZZ(2:4) = (/RHO_Z_P(II,JJ,KK-2:KK-1),WC%RHO_F*WC%ZZ_F(N)/)
               FZ(II,JJ,KK-2,N) = SCALAR_FACE_VALUE(WW(II,JJ,KK-2),ZZZ,FLUX_LIMITER)
            ENDIF
      END SELECT OFF_WALL_SELECT_3
   
   ENDIF

   SELECT CASE(IOR)
      CASE( 1)
         FX(II,JJ,KK,N)   = RHO_Z_P(IIG,JJG,KKG) ! zero out DU at wall
      CASE(-1)
         FX(II-1,JJ,KK,N) = RHO_Z_P(IIG,JJG,KKG)
      CASE( 2)
         FY(II,JJ,KK,N)   = RHO_Z_P(IIG,JJG,KKG)
      CASE(-2)
         FY(II,JJ-1,KK,N) = RHO_Z_P(IIG,JJG,KKG)
      CASE( 3)
         FZ(II,JJ,KK,N)   = RHO_Z_P(IIG,JJG,KKG)
      CASE(-3)
         FZ(II,JJ,KK-1,N) = RHO_Z_P(IIG,JJG,KKG)
   END SELECT

   ! Correct U_DOT_DEL_RHO_Z at the boundary

   IF (.NOT.CONSTANT_SPECIFIC_HEAT_RATIO) THEN

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
         CASE(SOLID_BOUNDARY,HVAC_BOUNDARY)
            IF (PREDICTOR) UN = -SIGN(1._EB,REAL(IOR,EB))*WC%ONE_D%UWS
            IF (CORRECTOR) UN = -SIGN(1._EB,REAL(IOR,EB))*WC%ONE_D%UW
         CASE(INTERPOLATED_BOUNDARY)
            UN = UVW_SAVE(IW)
      END SELECT BOUNDARY_SELECT
   
      SELECT CASE(IOR)
         CASE( 1)
            DU_M = (WC%RHO_F*WC%ZZ_F(N) - RHO_Z_P(IIG,JJG,KKG))*UN
            U_DOT_DEL_RHO_Z(IIG,JJG,KKG) = U_DOT_DEL_RHO_Z(IIG,JJG,KKG) - DU_M*WC%RDN
         CASE(-1)
            DU_P = (WC%RHO_F*WC%ZZ_F(N) - RHO_Z_P(IIG,JJG,KKG))*UN
            U_DOT_DEL_RHO_Z(IIG,JJG,KKG) = U_DOT_DEL_RHO_Z(IIG,JJG,KKG) + DU_P*WC%RDN
         CASE( 2)
            DU_M = (WC%RHO_F*WC%ZZ_F(N) - RHO_Z_P(IIG,JJG,KKG))*UN
            U_DOT_DEL_RHO_Z(IIG,JJG,KKG) = U_DOT_DEL_RHO_Z(IIG,JJG,KKG) - DU_M*WC%RDN
         CASE(-2)
            DU_P = (WC%RHO_F*WC%ZZ_F(N) - RHO_Z_P(IIG,JJG,KKG))*UN
            U_DOT_DEL_RHO_Z(IIG,JJG,KKG) = U_DOT_DEL_RHO_Z(IIG,JJG,KKG) + DU_P*WC%RDN
         CASE( 3)
            DU_M = (WC%RHO_F*WC%ZZ_F(N) - RHO_Z_P(IIG,JJG,KKG))*UN
            U_DOT_DEL_RHO_Z(IIG,JJG,KKG) = U_DOT_DEL_RHO_Z(IIG,JJG,KKG) - DU_M*WC%RDN
         CASE(-3)
            DU_P = (WC%RHO_F*WC%ZZ_F(N) - RHO_Z_P(IIG,JJG,KKG))*UN
            U_DOT_DEL_RHO_Z(IIG,JJG,KKG) = U_DOT_DEL_RHO_Z(IIG,JJG,KKG) + DU_P*WC%RDN
      END SELECT

   ENDIF
      
ENDDO WALL_LOOP

IF (.NOT.CONSTANT_SPECIFIC_HEAT_RATIO) THEN

   IF (.NOT.ENTHALPY_TRANSPORT) THEN
      U_DOT_DEL_RHO_Z = 0._EB
      RETURN
   ENDIF

   !$OMP PARALLEL DO PRIVATE(DU_P, DU_M) SCHEDULE(static)
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
   
            DU_P = (FX(I,J,K,N)   - RHO_Z_P(I,J,K))*UU(I,J,K)
            DU_M = (FX(I-1,J,K,N) - RHO_Z_P(I,J,K))*UU(I-1,J,K)
            U_DOT_DEL_RHO_Z(I,J,K) = U_DOT_DEL_RHO_Z(I,J,K) + (DU_P-DU_M)*RDX(I)
   
            DU_P = (FY(I,J,K,N)   - RHO_Z_P(I,J,K))*VV(I,J,K)
            DU_M = (FY(I,J-1,K,N) - RHO_Z_P(I,J,K))*VV(I,J-1,K)
            U_DOT_DEL_RHO_Z(I,J,K) = U_DOT_DEL_RHO_Z(I,J,K) + (DU_P-DU_M)*RDY(J)
   
            DU_P = (FZ(I,J,K,N)   - RHO_Z_P(I,J,K))*WW(I,J,K)
            DU_M = (FZ(I,J,K-1,N) - RHO_Z_P(I,J,K))*WW(I,J,K-1)
            U_DOT_DEL_RHO_Z(I,J,K) = U_DOT_DEL_RHO_Z(I,J,K) + (DU_P-DU_M)*RDZ(K)
   
         ENDDO
      ENDDO 
   ENDDO
   !$OMP END PARALLEL DO

ENDIF

END SUBROUTINE SPECIES_ADVECTION

END SUBROUTINE DIVERGENCE_PART_1


SUBROUTINE DIVERGENCE_PART_2(NM)

! Finish computing the divergence of the flow, D, and then compute its time derivative, DDDT

USE COMP_FUNCTIONS, ONLY: SECOND
INTEGER, INTENT(IN) :: NM
REAL(EB), POINTER, DIMENSION(:,:,:) :: DP,D_NEW,RTRM,DIV
REAL(EB) :: USUM_ADD(N_ZONE)
REAL(EB) :: RDT,TNOW,P_EQ,SUM_P_PSUM,SUM_USUM,SUM_DSUM,SUM_PSUM
LOGICAL :: OPEN_ZONE
REAL(EB), POINTER, DIMENSION(:) :: D_PBAR_DT_P
REAL(EB), POINTER, DIMENSION(:,:) :: PBAR_P
INTEGER :: IW,IOR,II,JJ,KK,IIG,JJG,KKG,IC,I,J,K,IPZ,IOPZ
TYPE(WALL_TYPE), POINTER :: WC=>NULL()

IF (SOLID_PHASE_ONLY) RETURN
IF (PERIODIC_TEST==3) RETURN
IF (PERIODIC_TEST==4) RETURN

TNOW=SECOND()
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

PRESSURE_ZONE_LOOP: DO IPZ=1,N_ZONE

   IF (EVACUATION_ONLY(NM) .AND. .NOT.EVACUATION_GRID(NM)) CYCLE PRESSURE_ZONE_LOOP

   IF (PREDICTOR) D_PBAR_DT_P => D_PBAR_DT_S
   IF (CORRECTOR) D_PBAR_DT_P => D_PBAR_DT

   ! Compute change in background pressure
 
   IF (ABS(PSUM(IPZ,NM)) > TWO_EPSILON_EB) THEN
      D_PBAR_DT_P(IPZ) = (DSUM(IPZ,NM) - USUM(IPZ,NM))/PSUM(IPZ,NM)
      IF (CORRECTOR) P_ZONE(IPZ)%DPSTAR = D_PBAR_DT_P(IPZ)
   ENDIF

   ! Add pressure derivative to divergence

   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (PRESSURE_ZONE(I,J,K) /= IPZ) CYCLE 
            IF (SOLID(CELL_INDEX(I,J,K)))    CYCLE
            DP(I,J,K) = DP(I,J,K) - (R_PBAR(K,IPZ)-RTRM(I,J,K))*D_PBAR_DT_P(IPZ)
         ENDDO
      ENDDO
   ENDDO

ENDDO PRESSURE_ZONE_LOOP

! Zero out divergence in solid cells

SOLID_LOOP: DO IC=1,CELL_COUNT
   IF (.NOT.SOLID(IC)) CYCLE SOLID_LOOP
   I = I_CELL(IC)
   J = J_CELL(IC)
   K = K_CELL(IC)
   DP(I,J,K) = 0._EB
ENDDO SOLID_LOOP

! Specify divergence in boundary cells to account for volume being generated at the walls

BC_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
   WC => WALL(IW)
   IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE BC_LOOP
   II = WC%ONE_D%II
   JJ = WC%ONE_D%JJ
   KK = WC%ONE_D%KK
   SELECT CASE (WC%BOUNDARY_TYPE)
      CASE (SOLID_BOUNDARY,HVAC_BOUNDARY)
         IF (.NOT.SOLID(CELL_INDEX(II,JJ,KK))) CYCLE BC_LOOP
         IOR = WC%ONE_D%IOR
         SELECT CASE(IOR)
            CASE( 1)
               DP(II,JJ,KK) = DP(II,JJ,KK) - WC%ONE_D%UWS*RDX(II)*RRN(II)*R(II)
            CASE(-1)
               DP(II,JJ,KK) = DP(II,JJ,KK) - WC%ONE_D%UWS*RDX(II)*RRN(II)*R(II-1)
            CASE( 2)
               DP(II,JJ,KK) = DP(II,JJ,KK) - WC%ONE_D%UWS*RDY(JJ)
            CASE(-2)
               DP(II,JJ,KK) = DP(II,JJ,KK) - WC%ONE_D%UWS*RDY(JJ)
            CASE( 3)
               DP(II,JJ,KK) = DP(II,JJ,KK) - WC%ONE_D%UWS*RDZ(KK)
            CASE(-3)
               DP(II,JJ,KK) = DP(II,JJ,KK) - WC%ONE_D%UWS*RDZ(KK)
         END SELECT
      CASE (OPEN_BOUNDARY,MIRROR_BOUNDARY,INTERPOLATED_BOUNDARY) ! should INTERPOLATED_BOUNDARY be here?
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
   
   NO_SCARC_IF: IF (PRES_METHOD /='SCARC') THEN
      DO IW=1,N_EXTERNAL_WALL_CELLS
         IF (WALL(IW)%NOM==0) CYCLE
         IIG = WALL(IW)%ONE_D%IIG
         JJG = WALL(IW)%ONE_D%JJG
         KKG = WALL(IW)%ONE_D%KKG
         IF (PREDICTOR) DDDT(IIG,JJG,KKG) = DDDT(IIG,JJG,KKG) + DS_CORR(IW)*RDT
         IF (CORRECTOR) DDDT(IIG,JJG,KKG) = DDDT(IIG,JJG,KKG) + (2._EB*D_CORR(IW)-DS_CORR(IW))*RDT
      ENDDO
   ENDIF NO_SCARC_IF

ENDIF TRUE_PROJECTION
   
TUSED(2,NM)=TUSED(2,NM)+SECOND()-TNOW
END SUBROUTINE DIVERGENCE_PART_2
 
 
SUBROUTINE CHECK_DIVERGENCE(NM)
USE COMP_FUNCTIONS, ONLY: SECOND 
! Computes maximum velocity divergence
 
INTEGER, INTENT(IN) :: NM
INTEGER  :: I,J,K
REAL(EB) :: DIV,RES,TNOW
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW,DP
 
TNOW=SECOND()
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
 
TUSED(2,NM)=TUSED(2,NM)+SECOND()-TNOW
END SUBROUTINE CHECK_DIVERGENCE


SUBROUTINE DIVERGENCE_PART_1_EVAC(T,NM)

USE COMP_FUNCTIONS, ONLY: SECOND 
USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
USE EVAC, ONLY: EVAC_EMESH_EXITS_TYPE, EMESH_EXITS, EMESH_NFIELDS, N_EXITS, N_CO_EXITS, N_DOORS
USE GEOMETRY_FUNCTIONS, ONLY: ASSIGN_PRESSURE_ZONE

! Compute contributions to the divergence term for evacuation meshes
 
INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: T
REAL(EB), POINTER, DIMENSION(:,:,:) :: DP,RHOP,RTRM
REAL(EB), POINTER, DIMENSION(:,:) :: PBAR_P            
REAL(EB) :: VC,VC1,TNOW,RDT,TSI,TIME_RAMP_FACTOR,ZONE_VOLUME,DELTA_P,PRES_RAMP_FACTOR,&
            X1,Y1,X2,Y2,Z1,Z2
TYPE(SURFACE_TYPE), POINTER :: SF
INTEGER :: IW,N,IOR,II,JJ,KK,IIG,JJG,KKG,ITMP,I,J,K,IPZ,I_VENT
TYPE(VENTS_TYPE), POINTER :: VT=>NULL()
TYPE(WALL_TYPE), POINTER :: WC=>NULL()
 
IF (SOLID_PHASE_ONLY) RETURN
IF (.NOT.EVACUATION_ONLY(NM)) RETURN

TNOW=SECOND()
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

! Evacuation flow field calculation: Change the outflow vent pressure zone and initialize everything

EVACUATION_ZONE_IF: IF (EVACUATION_ONLY(NM) .AND. EVACUATION_GRID(NM) .AND. PREDICTOR) THEN
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
      FRHO= 0._EB; 
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
                  IF (.NOT.SOLID(CELL_INDEX(I,J,K))) CALL ASSIGN_PRESSURE_ZONE(NM,XC(I),YC(J),ZC(K),N)
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
         WC=>WALL(IW)
         WC%PRESSURE_ZONE = 0
         WC%ONE_D%UW = 0._EB
         WC%U_TAU = 0._EB
         WC%RHO_F = RHO_0(1)
         WC%Y_PLUS = 1._EB
         WC%RHODW = 0.1_EB ! Do not initialize to zero to avoid divide by zero in the first time step
         II  = WC%ONE_D%II
         JJ  = WC%ONE_D%JJ
         KK  = WC%ONE_D%KK
         IIG = WC%ONE_D%IIG
         JJG = WC%ONE_D%JJG
         KKG = WC%ONE_D%KKG
         IF (KK==1) WC%PRESSURE_ZONE = PRESSURE_ZONE(IIG,JJG,KKG)
      END DO
   END IF EVACUATION_NEW_FIELD
END IF EVACUATION_ZONE_IF

! Determine if pressure ZONEs have merged

CONNECTED_ZONES(:,:,NM) = .FALSE.

! Compute normal component of velocity at boundaries, UWS

PREDICT_NORMALS: IF (PREDICTOR) THEN
 
   WALL_LOOP3: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
      WC => WALL(IW)
      IOR = WC%ONE_D%IOR
      
      WALL_CELL_TYPE: SELECT CASE (WC%BOUNDARY_TYPE)         
         CASE (NULL_BOUNDARY)
            WC%ONE_D%UWS = 0._EB
         CASE (SOLID_BOUNDARY)
            SF => SURFACE(WC%SURF_INDEX)
            IF (ABS(WC%ONE_D%T_IGN-T_BEGIN) < SPACING(WC%ONE_D%T_IGN) .AND. SF%RAMP_INDEX(TIME_VELO)>=1) THEN
               TSI = T + DT
            ELSE
               TSI = T + DT - WC%ONE_D%T_IGN
               IF (TSI<0._EB) THEN
                  WC%ONE_D%UWS = 0._EB
                  CYCLE WALL_LOOP3
               ENDIF
            ENDIF
            TIME_RAMP_FACTOR = EVALUATE_RAMP(TSI,SF%TAU(TIME_VELO),SF%RAMP_INDEX(TIME_VELO))
            KK               = WC%ONE_D%KK
            DELTA_P          = PBAR_P(KK,SF%DUCT_PATH(1)) - PBAR_P(KK,SF%DUCT_PATH(2))
            PRES_RAMP_FACTOR = SIGN(1._EB,SF%MAX_PRESSURE-DELTA_P)*SQRT(ABS((DELTA_P-SF%MAX_PRESSURE)/SF%MAX_PRESSURE))
            SELECT CASE(IOR) 
               CASE( 1)
                  WC%ONE_D%UWS =-U0 + TIME_RAMP_FACTOR*(WC%UW0+U0)
               CASE(-1)
                  WC%ONE_D%UWS = U0 + TIME_RAMP_FACTOR*(WC%UW0-U0)
               CASE( 2)
                  WC%ONE_D%UWS =-V0 + TIME_RAMP_FACTOR*(WC%UW0+V0)
               CASE(-2)
                  WC%ONE_D%UWS = V0 + TIME_RAMP_FACTOR*(WC%UW0-V0)
               CASE( 3)
                  WC%ONE_D%UWS =-W0 + TIME_RAMP_FACTOR*(WC%UW0+W0)
               CASE(-3)
                  WC%ONE_D%UWS = W0 + TIME_RAMP_FACTOR*(WC%UW0-W0)
            END SELECT
            ! Special Cases
            NEUMANN_IF: IF (SF%SPECIFIED_NORMAL_GRADIENT) THEN
               IIG = WC%ONE_D%IIG
               JJG = WC%ONE_D%JJG
               KKG = WC%ONE_D%KKG
               SELECT CASE(IOR) 
                  CASE( 1)
                     WC%ONE_D%UWS =-(U(IIG,JJG,KKG)   + SF%VEL_GRAD*WC%RDN)
                  CASE(-1)
                     WC%ONE_D%UWS = (U(IIG-1,JJG,KKG) + SF%VEL_GRAD*WC%RDN)
                  CASE( 2)
                     WC%ONE_D%UWS =-(V(IIG,JJG,KKG)   + SF%VEL_GRAD*WC%RDN)
                  CASE(-2)
                     WC%ONE_D%UWS = (V(IIG,JJG-1,KKG) + SF%VEL_GRAD*WC%RDN)
                  CASE( 3)
                     WC%ONE_D%UWS =-(W(IIG,JJG,KKG)   + SF%VEL_GRAD*WC%RDN)
                  CASE(-3)
                     WC%ONE_D%UWS = (W(IIG,JJG,KKG-1) + SF%VEL_GRAD*WC%RDN)
               END SELECT
            ENDIF NEUMANN_IF
            IF (ABS(SURFACE(WC%SURF_INDEX)%MASS_FLUX_TOTAL)>=TWO_EPSILON_EB) WC%ONE_D%UWS = WC%ONE_D%UWS*RHOA/WC%RHO_F
            IF (WC%VENT_INDEX>0) THEN 
               VT=>VENTS(WC%VENT_INDEX)
               IF (VT%N_EDDY>0) THEN ! Synthetic Eddy Method
                  II = WC%ONE_D%II
                  JJ = WC%ONE_D%JJ
                  KK = WC%ONE_D%KK
                  SELECT CASE(IOR)
                     CASE( 1)
                        WC%ONE_D%UWS = WC%ONE_D%UWS - TIME_RAMP_FACTOR*VT%U_EDDY(JJ,KK)
                     CASE(-1)
                        WC%ONE_D%UWS = WC%ONE_D%UWS + TIME_RAMP_FACTOR*VT%U_EDDY(JJ,KK)
                     CASE( 2)
                        WC%ONE_D%UWS = WC%ONE_D%UWS - TIME_RAMP_FACTOR*VT%V_EDDY(II,KK)
                     CASE(-2)
                        WC%ONE_D%UWS = WC%ONE_D%UWS + TIME_RAMP_FACTOR*VT%V_EDDY(II,KK)
                     CASE( 3)
                        WC%ONE_D%UWS = WC%ONE_D%UWS - TIME_RAMP_FACTOR*VT%W_EDDY(II,JJ)
                     CASE(-3)
                        WC%ONE_D%UWS = WC%ONE_D%UWS + TIME_RAMP_FACTOR*VT%W_EDDY(II,JJ)
                  END SELECT
               ENDIF
               EVACUATION_IF: IF (EVACUATION_ONLY(NM) .AND. EVACUATION_GRID(NM)) THEN
                  II = EVAC_TIME_ITERATIONS / MAXVAL(EMESH_NFIELDS)
                  IF ((ABS(ICYC)+1) > (WC%VENT_INDEX-1)*II .AND. (ABS(ICYC)+1) <= WC%VENT_INDEX*II) THEN
                     TSI = T + DT - (MAXVAL(EMESH_NFIELDS)-WC%VENT_INDEX)*II*EVAC_DT_FLOWFIELD
                     TIME_RAMP_FACTOR = EVALUATE_RAMP(TSI,SF%TAU(TIME_VELO),SF%RAMP_INDEX(TIME_VELO))
                  ELSE
                     TIME_RAMP_FACTOR = 0.0_EB
                  END IF
                  WC%ONE_D%UWS = TIME_RAMP_FACTOR*WC%UW0
               END IF EVACUATION_IF
            ENDIF
         CASE(OPEN_BOUNDARY,INTERPOLATED_BOUNDARY)
            II = WC%ONE_D%II
            JJ = WC%ONE_D%JJ
            KK = WC%ONE_D%KK
            SELECT CASE(IOR)
               CASE( 1)
                  WC%ONE_D%UWS = -U(II,JJ,KK)
               CASE(-1)
                  WC%ONE_D%UWS =  U(II-1,JJ,KK)
               CASE( 2)
                  WC%ONE_D%UWS = -V(II,JJ,KK)
               CASE(-2)
                  WC%ONE_D%UWS =  V(II,JJ-1,KK)
               CASE( 3)
                  WC%ONE_D%UWS = -W(II,JJ,KK)
               CASE(-3)
                  WC%ONE_D%UWS =  W(II,JJ,KK-1)
            END SELECT
      END SELECT WALL_CELL_TYPE
   ENDDO WALL_LOOP3

   DO IW=1,N_EXTERNAL_WALL_CELLS
      WALL(IW)%DUWDT = RDT*(WALL(IW)%ONE_D%UWS-WALL(IW)%ONE_D%UW)
   ENDDO
   
ELSE PREDICT_NORMALS
   
   DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
      WALL(IW)%ONE_D%UW = WALL(IW)%ONE_D%UWS
   ENDDO

ENDIF PREDICT_NORMALS

! Calculate pressure rise in each of the pressure zones by summing divergence expression over each zone

PRESSURE_ZONE_LOOP: DO IPZ=1,N_ZONE

   USUM(IPZ,NM) = 0._EB
   DSUM(IPZ,NM) = 0._EB
   PSUM(IPZ,NM) = 0._EB
   ZONE_VOLUME  = 0._EB

   IF (EVACUATION_ONLY(NM) .AND. .NOT.EVACUATION_GRID(NM)) CYCLE PRESSURE_ZONE_LOOP
   IF (.NOT.P_ZONE(IPZ)%EVACUATION) CYCLE PRESSURE_ZONE_LOOP
   RTRM=1._EB

   DO K=1,KBAR
      DO J=1,JBAR
         VC1 = DY(J)*DZ(K)
         DO I=1,IBAR
            IF (PRESSURE_ZONE(I,J,K) /= IPZ) CYCLE
            IF (SOLID(CELL_INDEX(I,J,K)))    CYCLE
            VC   = DX(I)*RC(I)*VC1
            ZONE_VOLUME = ZONE_VOLUME + VC
            DSUM(IPZ,NM) = DSUM(IPZ,NM) + VC*DP(I,J,K)
            PSUM(IPZ,NM) = PSUM(IPZ,NM) + VC*(R_PBAR(K,IPZ)-RTRM(I,J,K))
         ENDDO
      ENDDO
   ENDDO

   ! Calculate the volume flux to the boundary of the pressure zone (int u dot dA)

   WALL_LOOP4: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
      IF (WALL(IW)%PRESSURE_ZONE/=IPZ)            CYCLE WALL_LOOP4
      IF (WALL(IW)%BOUNDARY_TYPE/=SOLID_BOUNDARY) CYCLE WALL_LOOP4
      USUM(IPZ,NM) = USUM(IPZ,NM) + WALL(IW)%ONE_D%UWS*WALL(IW)%AW
   ENDDO WALL_LOOP4
   
ENDDO PRESSURE_ZONE_LOOP

TUSED(2,NM)=TUSED(2,NM)+SECOND()-TNOW

END SUBROUTINE DIVERGENCE_PART_1_EVAC

SUBROUTINE GET_REV_divg(MODULE_REV,MODULE_DATE)
INTEGER,INTENT(INOUT) :: MODULE_REV
CHARACTER(255),INTENT(INOUT) :: MODULE_DATE

WRITE(MODULE_DATE,'(A)') divgrev(INDEX(divgrev,':')+2:LEN_TRIM(divgrev)-2)
READ (MODULE_DATE,'(I5)') MODULE_REV
WRITE(MODULE_DATE,'(A)') divgdate

END SUBROUTINE GET_REV_divg
 
END MODULE DIVG

