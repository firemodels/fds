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
USE PHYSICAL_FUNCTIONS, ONLY: GET_CONDUCTIVITY,GET_SPECIFIC_HEAT,GET_AVERAGE_SPECIFIC_HEAT,GET_AVERAGE_SPECIFIC_HEAT_DIFF
USE EVAC, ONLY: EVAC_EMESH_EXITS_TYPE, EMESH_EXITS, EMESH_NFIELDS, EVAC_FDS6

! Compute contributions to the divergence term
 
INTEGER, INTENT(IN) :: NM
REAL(EB), POINTER, DIMENSION(:,:,:) :: KDTDX,KDTDY,KDTDZ,DP,KP, &
          RHO_D_DZDX,RHO_D_DZDY,RHO_D_DZDZ,RHO_D,RHOP,H_RHO_D_DZDX,H_RHO_D_DZDY,H_RHO_D_DZDZ,RTRM
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZP
REAL(EB) :: DELKDELT,VC,DTDX,DTDY,DTDZ,TNOW,ZZ_GET(0:N_TRACKED_SPECIES), &
            HDIFF,DZDX,DZDY,DZDZ,T,RDT,RHO_D_DZDN,TSI,TIME_RAMP_FACTOR,ZONE_VOLUME,DELTA_P,PRES_RAMP_FACTOR,&
            CP,CPBAR,TMP_G,TMP_WGT,DIV_DIFF_HEAT_FLUX
TYPE(SURFACE_TYPE), POINTER :: SF
INTEGER :: IW,N,IOR,II,JJ,KK,IIG,JJG,KKG,ITMP,IBC,I,J,K,IPZ,IOPZ
TYPE(VENTS_TYPE), POINTER :: VT
 
IF (SOLID_PHASE_ONLY) RETURN
IF (PERIODIC_TEST==3) RETURN
IF (PERIODIC_TEST==4) RETURN

TNOW=SECOND()
CALL POINT_TO_MESH(NM)
 
RDT = 1._EB/DT
 
SELECT CASE(PREDICTOR)
   CASE(.TRUE.)  
      DP => DS   
      R_PBAR = 1._EB/PBAR
      RHOP => RHOS
   CASE(.FALSE.) 
      DP => DDDT 
      R_PBAR = 1._EB/PBAR_S
      RHOP => RHO
END SELECT

! Determine if pressure ZONEs have merged

CONNECTED_ZONES(:,:,NM) = .FALSE.
ASUM(:,:,NM)    = 0._EB

DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
   IF (BOUNDARY_TYPE(IW)/=NULL_BOUNDARY .AND. BOUNDARY_TYPE(IW)/=OPEN_BOUNDARY .AND. BOUNDARY_TYPE(IW)/=INTERPOLATED_BOUNDARY) CYCLE
   IF (EVACUATION_ONLY(NM)) CYCLE
   II  = IJKW(1,IW)
   JJ  = IJKW(2,IW)
   KK  = IJKW(3,IW)
   IIG = IJKW(6,IW)
   JJG = IJKW(7,IW)
   KKG = IJKW(8,IW)
   IF (SOLID(CELL_INDEX(IIG,JJG,KKG))) CYCLE
   IPZ  = PRESSURE_ZONE(IIG,JJG,KKG)
   IOPZ = PRESSURE_ZONE(II,JJ,KK)
   IF (IW>N_EXTERNAL_WALL_CELLS .AND. IPZ/=IOPZ) THEN
      CONNECTED_ZONES(IOPZ,IPZ,NM) = .TRUE.
      CONNECTED_ZONES(IPZ,IOPZ,NM) = .TRUE.
      IF (IPZ>0) ASUM(IPZ,IOPZ,NM) = ASUM(IPZ,IOPZ,NM) + AW(IW)
   ENDIF
   IF (BOUNDARY_TYPE(IW)==OPEN_BOUNDARY) THEN
      CONNECTED_ZONES(0,IPZ,NM) = .TRUE.
      CONNECTED_ZONES(IPZ,0,NM) = .TRUE.
      IF (IPZ>0) ASUM(IPZ,0,NM) = ASUM(IPZ,0,NM) + AW(IW)
   ENDIF
   IF (BOUNDARY_TYPE(IW)==INTERPOLATED_BOUNDARY) THEN
      CONNECTED_ZONES(IOPZ,IPZ,NM) = .TRUE.
      CONNECTED_ZONES(IPZ,IOPZ,NM) = .TRUE.
      IF (IPZ>0) ASUM(IPZ,IOPZ,NM) = ASUM(IPZ,IOPZ,NM) + AW(IW)
   ENDIF
ENDDO

! Compute species-related finite difference terms

IF (N_TRACKED_SPECIES > 0 .AND. .NOT.EVACUATION_ONLY(NM)) THEN
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

! Zero out divergence to start

DP = 0._EB

IF (N_TRACKED_SPECIES > 0 .AND. .NOT.EVACUATION_ONLY(NM)) THEN
   DEL_RHO_D_DEL_Z = 0._EB
ENDIF

! Add species diffusion terms to divergence expression and compute diffusion term for species equations
 
IF (N_TRACKED_SPECIES > 0) THEN
   RHO_D => WORK4
   IF (LES) THEN
      RHO_D = MU*RSC
   ENDIF
ENDIF

SPECIES_LOOP: DO N=1,N_TRACKED_SPECIES

   IF (EVACUATION_ONLY(NM)) Cycle SPECIES_LOOP

   IF (DNS) THEN
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
   
   ! Numerical diffusion
   
   IF (INCLUDE_NUMERICAL_DIFFUSION) THEN
      RHO_D_DZDX = RHO_D_DZDX - DFX(:,:,:,N)
      RHO_D_DZDY = RHO_D_DZDY - DFY(:,:,:,N)
      RHO_D_DZDZ = RHO_D_DZDZ - DFZ(:,:,:,N)
   ENDIF

   ! Correct rho*D del Z at boundaries and store rho*D at boundaries

   WALL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
      IF (BOUNDARY_TYPE(IW)==NULL_BOUNDARY .OR. &
          BOUNDARY_TYPE(IW)==OPEN_BOUNDARY .OR. BOUNDARY_TYPE(IW)==INTERPOLATED_BOUNDARY) CYCLE WALL_LOOP
      IIG = IJKW(6,IW) 
      JJG = IJKW(7,IW) 
      KKG = IJKW(8,IW) 
      RHO_D_DZDN  = 2._EB*RHODW(IW,N)*(ZZP(IIG,JJG,KKG,N)-ZZ_F(IW,N))*RDN(IW)
      IOR = IJKW(4,IW)
      SELECT CASE(IOR) 
         CASE( 1)
            RHO_D_DZDX(IIG-1,JJG,KKG) =  RHO_D_DZDN
         CASE(-1)
            RHO_D_DZDX(IIG,JJG,KKG)   = -RHO_D_DZDN
         CASE( 2)
            RHO_D_DZDY(IIG,JJG-1,KKG) =  RHO_D_DZDN
         CASE(-2)
            RHO_D_DZDY(IIG,JJG,KKG)   = -RHO_D_DZDN
         CASE( 3)
            RHO_D_DZDZ(IIG,JJG,KKG-1) =  RHO_D_DZDN
         CASE(-3)
            RHO_D_DZDZ(IIG,JJG,KKG)   = -RHO_D_DZDN
      END SELECT
      
   ENDDO WALL_LOOP

   ! Compute del dot h_n*rho*D del Z_n (part of del dot qdot")

   H_RHO_D_DZDX => WORK5
   H_RHO_D_DZDY => WORK6
   H_RHO_D_DZDZ => WORK7

   ZZ_GET    = 0._EB
   ZZ_GET(N) = 1._EB

   DO K=0,KBAR
      DO J=0,JBAR
         DO I=0,IBAR
            ! H_RHO_D_DZDX
            TMP_G = .5_EB*(TMP(I+1,J,K)+TMP(I,J,K))
            CALL GET_AVERAGE_SPECIFIC_HEAT_DIFF(N,TMP_G,CPBAR)
            HDIFF = CPBAR*TMP_G
            H_RHO_D_DZDX(I,J,K) = HDIFF*RHO_D_DZDX(I,J,K)
            
            ! H_RHO_D_DZDY
            TMP_G = .5_EB*(TMP(I,J+1,K)+TMP(I,J,K))
            CALL GET_AVERAGE_SPECIFIC_HEAT_DIFF(N,TMP_G,CPBAR)
            HDIFF = CPBAR*TMP_G
            H_RHO_D_DZDY(I,J,K) = HDIFF*RHO_D_DZDY(I,J,K)
            
            ! H_RHO_D_DZDZ
            TMP_G = .5_EB*(TMP(I,J,K+1)+TMP(I,J,K))               
            CALL GET_AVERAGE_SPECIFIC_HEAT_DIFF(N,TMP_G,CPBAR)
            HDIFF = CPBAR*TMP_G
            H_RHO_D_DZDZ(I,J,K) = HDIFF*RHO_D_DZDZ(I,J,K)
         ENDDO
      ENDDO
   ENDDO

   WALL_LOOP2: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
      IF (BOUNDARY_TYPE(IW)==NULL_BOUNDARY .OR. &
          BOUNDARY_TYPE(IW)==OPEN_BOUNDARY .OR. BOUNDARY_TYPE(IW)==INTERPOLATED_BOUNDARY) CYCLE WALL_LOOP2
      IIG = IJKW(6,IW)
      JJG = IJKW(7,IW)
      KKG = IJKW(8,IW)
      IOR = IJKW(4,IW)
      CALL GET_AVERAGE_SPECIFIC_HEAT_DIFF(N,TMP_F(IW),CPBAR)
      HDIFF = CPBAR*TMP_F(IW)

      RHO_D_DZDN = 2._EB*RHODW(IW,N)*(ZZP(IIG,JJG,KKG,N)-ZZ_F(IW,N))*RDN(IW)
      SELECT CASE(IOR)
         CASE( 1) 
            H_RHO_D_DZDX(IIG-1,JJG,KKG) =  HDIFF*RHO_D_DZDN
         CASE(-1) 
            H_RHO_D_DZDX(IIG,JJG,KKG)   = -HDIFF*RHO_D_DZDN
         CASE( 2) 
            H_RHO_D_DZDY(IIG,JJG-1,KKG) =  HDIFF*RHO_D_DZDN
         CASE(-2) 
            H_RHO_D_DZDY(IIG,JJG,KKG)   = -HDIFF*RHO_D_DZDN
         CASE( 3) 
            H_RHO_D_DZDZ(IIG,JJG,KKG-1) =  HDIFF*RHO_D_DZDN
         CASE(-3) 
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

   ! Compute -Sum h_n del dot rho*D del Z_n
   
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            CALL GET_AVERAGE_SPECIFIC_HEAT_DIFF(N,TMP(I,J,K),CPBAR)
            HDIFF = CPBAR*TMP(I,J,K)
            DP(I,J,K) = DP(I,J,K) - HDIFF*DEL_RHO_D_DEL_Z(I,J,K,N)
         ENDDO
      ENDDO
   ENDDO

ENDDO SPECIES_LOOP

! Compute del dot k del T
 
ENERGY: IF (.NOT.EVACUATION_ONLY(NM)) THEN
 
   KDTDX => WORK1
   KDTDY => WORK2
   KDTDZ => WORK3
   KP    => WORK4
   
   KP = K_Z(NINT(TMPA),0)*SPECIES_MIXTURE(0)%MW
   
   ! Compute thermal conductivity k (KP)
 
   K_DNS_OR_LES: IF (DNS .AND. .NOT.EVACUATION_ONLY(NM)) THEN

      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
               IF (N_TRACKED_SPECIES>0) ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I,J,K,1:N_TRACKED_SPECIES)
               CALL GET_CONDUCTIVITY(ZZ_GET,KP(I,J,K),TMP(I,J,K)) 
            ENDDO
         ENDDO
      ENDDO

      BOUNDARY_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS
         II  = IJKW(1,IW)
         JJ  = IJKW(2,IW)
         KK  = IJKW(3,IW)
         IIG = IJKW(6,IW)
         JJG = IJKW(7,IW)
         KKG = IJKW(8,IW)
         KP(II,JJ,KK) = KP(IIG,JJG,KKG)
      ENDDO BOUNDARY_LOOP

   ELSE K_DNS_OR_LES
    
      CP_FTMP_IF: IF (CP_FTMP) THEN

          DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
                  IF (N_TRACKED_SPECIES > 0) ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I,J,K,1:N_TRACKED_SPECIES)
                  CALL GET_SPECIFIC_HEAT(ZZ_GET,CP,TMP(I,J,K)) 
                  KP(I,J,K) = MU(I,J,K)*CP*RPR  
               ENDDO
            ENDDO
         ENDDO

         BOUNDARY_LOOP2: DO IW=1,N_EXTERNAL_WALL_CELLS
            II  = IJKW(1,IW)
            JJ  = IJKW(2,IW)
            KK  = IJKW(3,IW)
            IIG = IJKW(6,IW)
            JJG = IJKW(7,IW)
            KKG = IJKW(8,IW)
            KP(II,JJ,KK) = KP(IIG,JJG,KKG)
         ENDDO BOUNDARY_LOOP2

      ELSE CP_FTMP_IF

         KP = MU*CPOPR

      ENDIF CP_FTMP_IF
      
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
      IF (BOUNDARY_TYPE(IW)==NULL_BOUNDARY .OR. &
          BOUNDARY_TYPE(IW)==OPEN_BOUNDARY .OR. BOUNDARY_TYPE(IW)==INTERPOLATED_BOUNDARY) CYCLE CORRECTION_LOOP
      II  = IJKW(1,IW) 
      JJ  = IJKW(2,IW) 
      KK  = IJKW(3,IW) 
      IIG = IJKW(6,IW) 
      JJG = IJKW(7,IW) 
      KKG = IJKW(8,IW) 
      KW(IW) = KP(IIG,JJG,KKG)
      IOR = IJKW(4,IW)
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
      DP(IIG,JJG,KKG) = DP(IIG,JJG,KKG) - QCONF(IW)*RDN(IW)
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
 
ENDIF ENERGY

! Compute RTRM = R*sum(Z_i/M_i)/(PBAR*C_P) and multiply it by divergence terms already summed up
 
RTRM => WORK1

IF (N_TRACKED_SPECIES>0 .AND. EVACUATION_ONLY(NM)) ZZ_GET(1:N_TRACKED_SPECIES) = 0._EB
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
         IF (N_TRACKED_SPECIES>0 .AND. .NOT.EVACUATION_ONLY(NM)) ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I,J,K,1:N_TRACKED_SPECIES)
         CALL GET_SPECIFIC_HEAT(ZZ_GET,CP,TMP(I,J,K)) 
         RTRM(I,J,K) = R_PBAR(K,PRESSURE_ZONE(I,J,K))*RSUM(I,J,K)/CP
         DP(I,J,K) = RTRM(I,J,K)*DP(I,J,K)
      ENDDO
   ENDDO 
ENDDO

! Compute (Wbar/rho) Sum (1/W_n) del dot rho*D del Z_n

DO N=1,N_TRACKED_SPECIES
   IF (EVACUATION_ONLY(NM)) CYCLE
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            DP(I,J,K) = DP(I,J,K) + (SPECIES_MIXTURE(N)%RCON-SPECIES_MIXTURE(0)%RCON)/(RSUM(I,J,K)*RHOP(I,J,K))*&
                                    DEL_RHO_D_DEL_Z(I,J,K,N)
         ENDDO
      ENDDO
   ENDDO
ENDDO

! Remove numerical diffusion before using DEL_RHO_D_DEL_Z in scalar transport

NUMERICAL_DIFFUSION_IF: IF (INCLUDE_NUMERICAL_DIFFUSION) THEN
   DO N=1,N_TRACKED_SPECIES
      IF (EVACUATION_ONLY(NM)) CYCLE
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               DEL_RHO_D_DEL_Z(I,J,K,N) = DEL_RHO_D_DEL_Z(I,J,K,N) + (DFX(I,J,K,N)-DFX(I-1,J,K,N))*RDX(I) &
                                                                   + (DFY(I,J,K,N)-DFY(I,J-1,K,N))*RDY(J) &
                                                                   + (DFZ(I,J,K,N)-DFZ(I,J,K-1,N))*RDZ(K)
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDIF NUMERICAL_DIFFUSION_IF

! Add contribution of reactions
 
IF (N_REACTIONS > 0 .AND. .NOT.EVACUATION_ONLY(NM)) THEN
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            DP(I,J,K) = DP(I,J,K) + D_REACTION(I,J,K)
         ENDDO
      ENDDO
   ENDDO
ENDIF


! Add contribution of evaporating PARTICLEs

IF (NLP>0 .AND. N_EVAP_INDICES > 0 .AND. .NOT.EVACUATION_ONLY(NM)) THEN
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            DP(I,J,K) = DP(I,J,K) + D_LAGRANGIAN(I,J,K)
         ENDDO
      ENDDO
   ENDDO
ENDIF
 
! Atmospheric Stratification Term

IF (STRATIFICATION .AND. .NOT.EVACUATION_ONLY(NM)) THEN
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
            DP(I,J,K) = DP(I,J,K) + (RTRM(I,J,K)-R_PBAR(K,PRESSURE_ZONE(I,J,K)))*0.5_EB*(W(I,J,K)+W(I,J,K-1))*GVEC(3)*RHO_0(K)
         ENDDO
      ENDDO
   ENDDO
ENDIF

! Compute normal component of velocity at boundaries, UWS

PREDICT_NORMALS: IF (PREDICTOR) THEN
 
   WALL_LOOP3: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS

      IOR = IJKW(4,IW)
      
      WALL_CELL_TYPE: SELECT CASE (BOUNDARY_TYPE(IW))
         CASE (NULL_BOUNDARY)
            UWS(IW) = 0._EB
         CASE (SOLID_BOUNDARY)
            IBC = IJKW(5,IW)
            SF => SURFACE(IBC)
            EVAC_IF_NOT: IF (.NOT.EVACUATION_ONLY(NM)) THEN
            IF (SF%SPECIES_BC_INDEX==SPECIFIED_MASS_FLUX .OR. SF%SPECIES_BC_INDEX==INTERPOLATED_BC .OR. &
                SF%SPECIES_BC_INDEX==HVAC_BOUNDARY .OR. ANY(SF%LEAK_PATH>0._EB)) CYCLE WALL_LOOP3
            ENDIF EVAC_IF_NOT
            IF (ABS(TW(IW)-T_BEGIN) < SPACING(TW(IW)) .AND. SF%RAMP_INDEX(TIME_VELO)>=1) THEN
               TSI = T + DT
            ELSE
               TSI = T + DT - TW(IW)
               IF (TSI<0._EB) THEN
                  UWS(IW) = 0._EB
                  CYCLE WALL_LOOP3
               ENDIF
            ENDIF
            TIME_RAMP_FACTOR = EVALUATE_RAMP(TSI,SF%TAU(TIME_VELO),SF%RAMP_INDEX(TIME_VELO))
            KK               = IJKW(3,IW)
            DELTA_P          = PBAR(KK,SF%DUCT_PATH(1)) - PBAR(KK,SF%DUCT_PATH(2))
            PRES_RAMP_FACTOR = SIGN(1._EB,SF%MAX_PRESSURE-DELTA_P)*SQRT(ABS((DELTA_P-SF%MAX_PRESSURE)/SF%MAX_PRESSURE))
            SELECT CASE(IOR) 
               CASE( 1)
                  UWS(IW) =-U0 + TIME_RAMP_FACTOR*PRES_RAMP_FACTOR*(UW0(IW)+U0)
               CASE(-1)
                  UWS(IW) = U0 + TIME_RAMP_FACTOR*PRES_RAMP_FACTOR*(UW0(IW)-U0)
               CASE( 2)
                  UWS(IW) =-V0 + TIME_RAMP_FACTOR*PRES_RAMP_FACTOR*(UW0(IW)+V0)
               CASE(-2)
                  UWS(IW) = V0 + TIME_RAMP_FACTOR*PRES_RAMP_FACTOR*(UW0(IW)-V0)
               CASE( 3)
                  UWS(IW) =-W0 + TIME_RAMP_FACTOR*PRES_RAMP_FACTOR*(UW0(IW)+W0)
               CASE(-3)
                  UWS(IW) = W0 + TIME_RAMP_FACTOR*PRES_RAMP_FACTOR*(UW0(IW)-W0)
            END SELECT          
            ! Special Cases
            IF (EVACUATION_ONLY(NM) .AND. .NOT.EVAC_FDS6) UWS(IW) = TIME_RAMP_FACTOR*PRES_RAMP_FACTOR*UW0(IW)
            IF (ABS(SURFACE(IBC)%MASS_FLUX_TOTAL)>=ZERO_P) UWS(IW) = UWS(IW)*RHOA/RHO_F(IW)
            IF (VENT_INDEX(IW)>0) THEN 
               VT=>VENTS(VENT_INDEX(IW))
               IF (VT%N_EDDY>0) THEN ! Synthetic Eddy Method
                  II = IJKW(1,IW)
                  JJ = IJKW(2,IW)
                  KK = IJKW(3,IW)
                  SELECT CASE(IOR)
                     CASE( 1)
                        UWS(IW) = UWS(IW) - TIME_RAMP_FACTOR*PRES_RAMP_FACTOR*VT%U_EDDY(JJ,KK)
                     CASE(-1)
                        UWS(IW) = UWS(IW) + TIME_RAMP_FACTOR*PRES_RAMP_FACTOR*VT%U_EDDY(JJ,KK)
                     CASE( 2)
                        UWS(IW) = UWS(IW) - TIME_RAMP_FACTOR*PRES_RAMP_FACTOR*VT%V_EDDY(II,KK)
                     CASE(-2)
                        UWS(IW) = UWS(IW) + TIME_RAMP_FACTOR*PRES_RAMP_FACTOR*VT%V_EDDY(II,KK)
                     CASE( 3)
                        UWS(IW) = UWS(IW) - TIME_RAMP_FACTOR*PRES_RAMP_FACTOR*VT%W_EDDY(II,JJ)
                     CASE(-3)
                        UWS(IW) = UWS(IW) + TIME_RAMP_FACTOR*PRES_RAMP_FACTOR*VT%W_EDDY(II,JJ)
                  END SELECT
               ENDIF
               EVAC_IF: IF (EVACUATION_ONLY(NM) .AND. EVACUATION_GRID(NM) .AND. EVAC_FDS6) THEN
                  II = EVAC_TIME_ITERATIONS / MAXVAL(EMESH_NFIELDS)
                  IF ((ABS(ICYC)+1) > (VENT_INDEX(IW)-1)*II .AND. (ABS(ICYC)+1) <= VENT_INDEX(IW)*II) THEN
                     TSI = T + DT - (MAXVAL(EMESH_NFIELDS)-VENT_INDEX(IW))*II*EVAC_DT_FLOWFIELD
                     TIME_RAMP_FACTOR = EVALUATE_RAMP(TSI,SF%TAU(TIME_VELO),SF%RAMP_INDEX(TIME_VELO))
                  ELSE
                     TIME_RAMP_FACTOR = 0.0_EB
                  END IF
                  UWS(IW) = TIME_RAMP_FACTOR*PRES_RAMP_FACTOR*UW0(IW)
               END IF EVAC_IF
            ENDIF
         CASE(OPEN_BOUNDARY,INTERPOLATED_BOUNDARY)
            II = IJKW(1,IW)
            JJ = IJKW(2,IW)
            KK = IJKW(3,IW)
            SELECT CASE(IOR)
               CASE( 1)
                  UWS(IW) = -U(II,JJ,KK)
               CASE(-1)
                  UWS(IW) =  U(II-1,JJ,KK)
               CASE( 2)
                  UWS(IW) = -V(II,JJ,KK)
               CASE(-2)
                  UWS(IW) =  V(II,JJ-1,KK)
               CASE( 3)
                  UWS(IW) = -W(II,JJ,KK)
               CASE(-3)
                  UWS(IW) =  W(II,JJ,KK-1)
            END SELECT
      END SELECT WALL_CELL_TYPE
   ENDDO WALL_LOOP3

   DUWDT(1:N_EXTERNAL_WALL_CELLS) = RDT*(UWS(1:N_EXTERNAL_WALL_CELLS)-UW(1:N_EXTERNAL_WALL_CELLS))
   
ELSE PREDICT_NORMALS
   
   UW = UWS

ENDIF PREDICT_NORMALS

! Calculate pressure rise in each of the pressure zones by summing divergence expression over each zone

PRESSURE_ZONE_LOOP: DO IPZ=1,N_ZONE

   USUM(IPZ,NM) = 0._EB
   DSUM(IPZ,NM) = 0._EB
   PSUM(IPZ,NM) = 0._EB
   ZONE_VOLUME  = 0._EB

   IF (EVACUATION_ONLY(NM)) CYCLE PRESSURE_ZONE_LOOP

   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (PRESSURE_ZONE(I,J,K) /= IPZ) CYCLE
            IF (SOLID(CELL_INDEX(I,J,K)))    CYCLE
            VC   = DX(I)*RC(I)*DY(J)*DZ(K)
            ZONE_VOLUME = ZONE_VOLUME + VC
            DSUM(IPZ,NM) = DSUM(IPZ,NM) + VC*DP(I,J,K)
            PSUM(IPZ,NM) = PSUM(IPZ,NM) + VC*(R_PBAR(K,IPZ)-RTRM(I,J,K))
         ENDDO
      ENDDO
   ENDDO

   ! Calculate the volume flux to the boundary of the pressure zone (int u dot dA)

   WALL_LOOP4: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
      IF (PRESSURE_ZONE_WALL(IW)/=IPZ)     CYCLE WALL_LOOP4
      IF (BOUNDARY_TYPE(IW)/=SOLID_BOUNDARY) CYCLE WALL_LOOP4
      USUM(IPZ,NM) = USUM(IPZ,NM) + UWS(IW)*AW(IW)  
   ENDDO WALL_LOOP4
ENDDO PRESSURE_ZONE_LOOP

TUSED(2,NM)=TUSED(2,NM)+SECOND()-TNOW
END SUBROUTINE DIVERGENCE_PART_1



SUBROUTINE DIVERGENCE_PART_2(NM)

! Finish computing the divergence of the flow, D, and then compute its time derivative, DDDT

USE COMP_FUNCTIONS, ONLY: SECOND
INTEGER, INTENT(IN) :: NM
REAL(EB), POINTER, DIMENSION(:,:,:) :: DP,D_NEW,RTRM,DIV
REAL(EB) :: RDT,TNOW,P_EQ,P_EQ_NUM,P_EQ_DEN,RF
REAL(EB), POINTER, DIMENSION(:) :: D_PBAR_DT_P
INTEGER :: IW,IOR,II,JJ,KK,IIG,JJG,KKG,IC,I,J,K,IPZ,IOPZ,IOPZ2
REAL(EB), DIMENSION(N_ZONE) :: USUM_ADD

IF (SOLID_PHASE_ONLY) RETURN
IF (PERIODIC_TEST==3) RETURN
IF (PERIODIC_TEST==4) RETURN

TNOW=SECOND()
CALL POINT_TO_MESH(NM)

RDT = 1._EB/DT

SELECT CASE(PREDICTOR)
   CASE(.TRUE.)
      DP => DS
      R_PBAR = 1._EB/PBAR
   CASE(.FALSE.)
      DP => DDDT
      R_PBAR = 1._EB/PBAR_S
END SELECT

RTRM => WORK1

! Adjust volume flows (USUM) of pressure ZONEs that are connected to equalize background pressure

USUM_ADD = 0._EB

DO IPZ=1,N_ZONE
   IF (EVACUATION_ONLY(NM)) CYCLE
   P_EQ_NUM = PBAR(1,IPZ)*PSUM(IPZ,NM)
   P_EQ_DEN = PSUM(IPZ,NM)
   DO IOPZ=N_ZONE,0,-1
      IF (IOPZ==IPZ) CYCLE
      IF (CONNECTED_ZONES(IPZ,IOPZ,NM)) THEN
         IF (IOPZ==0) THEN
            P_EQ_NUM = P_0(1)
            P_EQ_DEN = 1._EB
         ELSE
            P_EQ_NUM = P_EQ_NUM + PBAR(1,IOPZ)*PSUM(IOPZ,NM)
            P_EQ_DEN = P_EQ_DEN + PSUM(IOPZ,NM)
         ENDIF
      ENDIF
   ENDDO
   P_EQ = P_EQ_NUM/P_EQ_DEN
   USUM_ADD(IPZ) = RDT*PSUM(IPZ,NM)*(PBAR(1,IPZ)-P_EQ)
ENDDO

DO IPZ=1,N_ZONE
   IF (EVACUATION_ONLY(NM)) CYCLE
   RF = PRESSURE_RELAX_FACTOR
   DO IOPZ=1,N_ZONE
      IF (IOPZ==IPZ) CYCLE
      IF (CONNECTED_ZONES(IPZ,IOPZ,NM)) THEN
         DO IOPZ2=N_ZONE,0,-1
            IF (IOPZ==IOPZ2) CYCLE
            IF (CONNECTED_ZONES(IOPZ,IOPZ2,NM)) THEN
               IF (IOPZ2/=0) ASUM(IOPZ,IOPZ2,NM) = MAX(ASUM(IOPZ,IOPZ2,NM),ASUM(IOPZ2,IOPZ,NM))
               IF (ASUM(IOPZ,IOPZ2,NM)>0._EB .AND. ABS(USUM_ADD(IOPZ))>=ZERO_P) &
                  RF=MIN(RF,10._EB*ASUM(IOPZ,IOPZ2,NM)/ABS(USUM_ADD(IOPZ)))
            ENDIF
         ENDDO
      ENDIF
   ENDDO
   USUM(IPZ,NM) = USUM(IPZ,NM) + RF*USUM_ADD(IPZ)
ENDDO

! Compute dP/dt for each pressure ZONE

PRESSURE_ZONE_LOOP: DO IPZ=1,N_ZONE

   IF (EVACUATION_ONLY(NM)) CYCLE PRESSURE_ZONE_LOOP

   IF (PREDICTOR) D_PBAR_DT_P => D_PBAR_S_DT
   IF (CORRECTOR) D_PBAR_DT_P => D_PBAR_DT

   ! Compute change in background pressure
 
   IF (ABS(PSUM(IPZ,NM)) > ZERO_P) THEN
      D_PBAR_DT_P(IPZ) = (DSUM(IPZ,NM) - USUM(IPZ,NM))/PSUM(IPZ,NM)
      IF (CORRECTOR) P_ZONE(IPZ)%DPSTAR = D_PBAR_DT_P(IPZ)
   ENDIF

   ! Add pressure derivative to divergence

   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (PRESSURE_ZONE(I,J,K) /= IPZ) CYCLE 
            IF (SOLID(CELL_INDEX(I,J,K)))    CYCLE
            DP(I,J,K) = DP(I,J,K) + (RTRM(I,J,K)-R_PBAR(K,IPZ))*D_PBAR_DT_P(IPZ)
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
   IF (BOUNDARY_TYPE(IW)==NULL_BOUNDARY) CYCLE BC_LOOP
   II = IJKW(1,IW)
   JJ = IJKW(2,IW)
   KK = IJKW(3,IW)
   SELECT CASE (BOUNDARY_TYPE(IW))
      CASE (SOLID_BOUNDARY)
         IF (.NOT.SOLID(CELL_INDEX(II,JJ,KK))) CYCLE BC_LOOP
         IOR = IJKW(4,IW)
         SELECT CASE(IOR)
            CASE( 1)
               DP(II,JJ,KK) = DP(II,JJ,KK) - UWS(IW)*RDX(II)*RRN(II)*R(II)
            CASE(-1)
               DP(II,JJ,KK) = DP(II,JJ,KK) - UWS(IW)*RDX(II)*RRN(II)*R(II-1)
            CASE( 2)
               DP(II,JJ,KK) = DP(II,JJ,KK) - UWS(IW)*RDY(JJ)
            CASE(-2)
               DP(II,JJ,KK) = DP(II,JJ,KK) - UWS(IW)*RDY(JJ)
            CASE( 3)
               DP(II,JJ,KK) = DP(II,JJ,KK) - UWS(IW)*RDZ(KK)
            CASE(-3)
               DP(II,JJ,KK) = DP(II,JJ,KK) - UWS(IW)*RDZ(KK)
         END SELECT
      CASE (OPEN_BOUNDARY,MIRROR_BOUNDARY,INTERPOLATED_BOUNDARY)
         IIG = IJKW(6,IW)
         JJG = IJKW(7,IW)
         KKG = IJKW(8,IW)
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
         IF (IJKW(9,IW)==0) CYCLE
         IIG = IJKW(6,IW)
         JJG = IJKW(7,IW)
         KKG = IJKW(8,IW)
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

SUBROUTINE GET_REV_divg(MODULE_REV,MODULE_DATE)
INTEGER,INTENT(INOUT) :: MODULE_REV
CHARACTER(255),INTENT(INOUT) :: MODULE_DATE

WRITE(MODULE_DATE,'(A)') divgrev(INDEX(divgrev,':')+1:LEN_TRIM(divgrev)-2)
READ (MODULE_DATE,'(I5)') MODULE_REV
WRITE(MODULE_DATE,'(A)') divgdate

END SUBROUTINE GET_REV_divg
 
END MODULE DIVG

