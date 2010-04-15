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
USE PHYSICAL_FUNCTIONS, ONLY: GET_DIFFUSIVITY,GET_CONDUCTIVITY,GET_SPECIFIC_HEAT,GET_AVERAGE_SPECIFIC_HEAT

! Compute contributions to the divergence term
 
INTEGER, INTENT(IN) :: NM
REAL(EB), POINTER, DIMENSION(:,:,:) :: KDTDX,KDTDY,KDTDZ,DP,KP, &
          RHO_D_DYDX,RHO_D_DYDY,RHO_D_DYDZ,RHO_D,RHOP,H_RHO_D_DYDX,H_RHO_D_DYDY,H_RHO_D_DYDZ,RTRM
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: YYP,FX,FY,FZ
REAL(EB) :: DELKDELT,VC,DTDX,DTDY,DTDZ,TNOW, YSUM,YY_GET(1:N_SPECIES),ZZ_GET(1:N_MIX_SPECIES), &
            HDIFF,DYDX,DYDY,DYDZ,T,RDT,RHO_D_DYDN,TSI,VDOT_LEAK,TIME_RAMP_FACTOR,ZONE_VOLUME,CP_MF,DELTA_P,PRES_RAMP_FACTOR,&
            H_G,H_G_A,TMP_G
TYPE(SURFACE_TYPE), POINTER :: SF
INTEGER :: IW,N,IOR,II,JJ,KK,IIG,JJG,KKG,ITMP,IBC,I,J,K,IPZ,IOPZ
 
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

DO IW=1,NWC
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
   IF (IW>NEWC .AND. IPZ/=IOPZ) THEN
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

IF (N_SPECIES > 0 .AND. .NOT.EVACUATION_ONLY(NM)) THEN
   RHO_D_DYDX  => WORK1
   RHO_D_DYDY  => WORK2
   RHO_D_DYDZ  => WORK3
   
   IF (FLUX_LIMITER>=0) THEN
      FX => SCALAR_SAVE1
      FY => SCALAR_SAVE2
      FZ => SCALAR_SAVE3
      !$OMP PARALLEL WORKSHARE
      FX = 0._EB
      FY = 0._EB
      FZ = 0._EB
      !$OMP END PARALLEL WORKSHARE
   ENDIF
   
   SELECT CASE(PREDICTOR)
      CASE(.TRUE.)  
         YYP => YYS 
      CASE(.FALSE.) 
         YYP => YY  
   END SELECT
ENDIF

! Zero out divergence to start

!$OMP PARALLEL
!$OMP WORKSHARE 
DP  = 0._EB
!$OMP END WORKSHARE
IF (N_SPECIES > 0 .AND. .NOT.EVACUATION_ONLY(NM)) THEN
   !$OMP WORKSHARE
   DEL_RHO_D_DEL_Y = 0._EB
   !$OMP END WORKSHARE
ENDIF
!$OMP END PARALLEL

! Add species diffusion terms to divergence expression and compute diffusion term for species equations
 
SPECIES_LOOP: DO N=1,N_SPECIES

   IF (EVACUATION_ONLY(NM)) Cycle SPECIES_LOOP
 
   ! Compute rho*D
 
   RHO_D => WORK4

   !$OMP PARALLEL SHARED(RHO_D) 
   IF (DNS .AND. SPECIES(N)%MODE/=MIXTURE_FRACTION_SPECIES) THEN
      !$OMP WORKSHARE
      RHO_D = 0._EB
      !$OMP END WORKSHARE
      !$OMP DO COLLAPSE(3) PRIVATE(K,J,I,ITMP)
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               !!$ IF ((K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_DIVG_01'
               ITMP = MIN(5000,NINT(TMP(I,J,K)))
               RHO_D(I,J,K) = RHOP(I,J,K)*SPECIES(N)%D(ITMP)
            ENDDO 
         ENDDO
      ENDDO
      !$OMP END DO
   ENDIF
    
   IF (DNS .AND. SPECIES(N)%MODE==MIXTURE_FRACTION_SPECIES) THEN
      !$OMP WORKSHARE
      RHO_D = 0._EB
      !$OMP END WORKSHARE
      !$OMP DO COLLAPSE(3) PRIVATE(K,J,I,ITMP,YSUM)
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               !!$ IF ((K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_DIVG_02'
               ITMP = MIN(5000,NINT(TMP(I,J,K)))
               YSUM = SUM(YYP(I,J,K,:)) - SUM(YYP(I,J,K,I_Z_MIN:I_Z_MAX))
               ZZ_GET(:) = YYP(I,J,K,I_Z_MIN:I_Z_MAX)
               CALL GET_DIFFUSIVITY(ZZ_GET,YSUM,RHO_D(I,J,K),ITMP)
               RHO_D(I,J,K) = RHOP(I,J,K)*RHO_D(I,J,K)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO
   ENDIF

   IF (LES) THEN
      !$OMP WORKSHARE
      RHO_D = MU*RSC
      !$OMP END WORKSHARE
   ENDIF
   
   ! Compute rho*D del Y

   !$OMP DO COLLAPSE(3) PRIVATE(K,J,I,DYDX,DYDY,DYDZ)
   DO K=0,KBAR
      DO J=0,JBAR
         DO I=0,IBAR
            !!$ IF ((K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_DIVG_03'
            DYDX = (YYP(I+1,J,K,N)-YYP(I,J,K,N))*RDXN(I)
            RHO_D_DYDX(I,J,K) = .5_EB*(RHO_D(I+1,J,K)+RHO_D(I,J,K))*DYDX
            DYDY = (YYP(I,J+1,K,N)-YYP(I,J,K,N))*RDYN(J)
            RHO_D_DYDY(I,J,K) = .5_EB*(RHO_D(I,J+1,K)+RHO_D(I,J,K))*DYDY
            DYDZ = (YYP(I,J,K+1,N)-YYP(I,J,K,N))*RDZN(K)
            RHO_D_DYDZ(I,J,K) = .5_EB*(RHO_D(I,J,K+1)+RHO_D(I,J,K))*DYDZ
         ENDDO
      ENDDO
   ENDDO
   !$OMP END DO
   
   IF (FLUX_LIMITER>=0) THEN
      !$OMP WORKSHARE
      FX(:,:,:,N) = -RHO_D_DYDX
      FY(:,:,:,N) = -RHO_D_DYDY
      FZ(:,:,:,N) = -RHO_D_DYDZ
      !$OMP END WORKSHARE
   ENDIF

   ! Correct rho*D del Y at boundaries and store rho*D at boundaries

   !$OMP DO PRIVATE(IW,IIG,JJG,KKG,RHO_D_DYDN,IOR,IBC,SF)
   WALL_LOOP: DO IW=1,NWC
      IF (BOUNDARY_TYPE(IW)==NULL_BOUNDARY .OR. BOUNDARY_TYPE(IW)==POROUS_BOUNDARY) CYCLE WALL_LOOP
      !!$ IF ((IW == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_DIVG_04'
      IIG = IJKW(6,IW) 
      JJG = IJKW(7,IW) 
      KKG = IJKW(8,IW) 
      RHODW(IW,N) = RHO_D(IIG,JJG,KKG)
      RHO_D_DYDN  = 2._EB*RHODW(IW,N)*(YYP(IIG,JJG,KKG,N)-YY_F(IW,N))*RDN(IW)
      DEL_RHO_D_DEL_Y(IIG,JJG,KKG,N) = DEL_RHO_D_DEL_Y(IIG,JJG,KKG,N) - RHO_D_DYDN*RDN(IW)
      IOR = IJKW(4,IW)
      SELECT CASE(IOR) 
         CASE( 1)
            RHO_D_DYDX(IIG-1,JJG,KKG) = 0._EB
         CASE(-1)
            RHO_D_DYDX(IIG,JJG,KKG)   = 0._EB
         CASE( 2)
            RHO_D_DYDY(IIG,JJG-1,KKG) = 0._EB
         CASE(-2)
            RHO_D_DYDY(IIG,JJG,KKG)   = 0._EB
         CASE( 3)
            RHO_D_DYDZ(IIG,JJG,KKG-1) = 0._EB
         CASE(-3)
            RHO_D_DYDZ(IIG,JJG,KKG)   = 0._EB
      END SELECT
      
      IF (FLUX_LIMITER>=0) FW(IW,N) = -SIGN(1._EB,REAL(IOR,EB))*RHO_D_DYDN
      
   ENDDO WALL_LOOP
   !$OMP END DO

   ! Compute del dot h_n*rho*D del Y_n only for non-mixture fraction cases
 
!   SPECIES_DIFFUSION: IF (.NOT.MIXTURE_FRACTION) THEN
   SPECIES_DIFFUSION: IF (SPECIES(N)%MODE/=MIXTURE_FRACTION_SPECIES) THEN

      !$OMP SINGLE
      H_RHO_D_DYDX => WORK4
      H_RHO_D_DYDY => WORK5
      H_RHO_D_DYDZ => WORK6
      !$OMP END SINGLE

      !$OMP DO COLLAPSE(3) PRIVATE(K,J,I,TMP_G,ITMP,YY_GET,H_G_A,H_G,HDIFF)
      DO K=0,KBAR
         DO J=0,JBAR
            DO I=0,IBAR
               !!$ IF ((K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_DIVG_05'
               TMP_G = .5_EB*(TMP(I+1,J,K)+TMP(I,J,K))
               ITMP = MIN(5000,NINT(TMP_G))
               YY_GET=0._EB
               CALL GET_AVERAGE_SPECIFIC_HEAT(YY_GET,H_G_A,ITMP)
               YY_GET(N) = 1._EB
               CALL GET_AVERAGE_SPECIFIC_HEAT(YY_GET,H_G,ITMP)               
               HDIFF = (H_G-H_G_A)*TMP_G
               H_RHO_D_DYDX(I,J,K) = HDIFF*RHO_D_DYDX(I,J,K)
               TMP_G = .5_EB*(TMP(I,J+1,K)+TMP(I,J,K))
               ITMP = MIN(5000,NINT(TMP_G))
               YY_GET=0._EB
               CALL GET_AVERAGE_SPECIFIC_HEAT(YY_GET,H_G_A,ITMP)
               YY_GET(N) = 1._EB
               CALL GET_AVERAGE_SPECIFIC_HEAT(YY_GET,H_G,ITMP)               
               HDIFF = (H_G-H_G_A)*TMP_G
               H_RHO_D_DYDY(I,J,K) = HDIFF*RHO_D_DYDY(I,J,K)
               TMP_G = .5_EB*(TMP(I,J,K+1)+TMP(I,J,K))               
               ITMP = MIN(5000,NINT(TMP_G))
               YY_GET=0._EB
               CALL GET_AVERAGE_SPECIFIC_HEAT(YY_GET,H_G_A,ITMP)
               YY_GET(N) = 1._EB
               CALL GET_AVERAGE_SPECIFIC_HEAT(YY_GET,H_G,ITMP)               
               HDIFF = (H_G-H_G_A)*TMP_G
               H_RHO_D_DYDZ(I,J,K) = HDIFF*RHO_D_DYDZ(I,J,K)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO

      !$OMP DO PRIVATE(IW,IIG,JJG,KKG,IOR,TMP_G,ITMP,YY_GET,H_G_A,H_G,HDIFF,RHO_D_DYDN)
      WALL_LOOP2: DO IW=1,NWC
         !!$ IF ((IW == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_DIVG_06'
         IF (BOUNDARY_TYPE(IW)==NULL_BOUNDARY .OR. BOUNDARY_TYPE(IW)==POROUS_BOUNDARY) CYCLE WALL_LOOP2
         IIG = IJKW(6,IW)
         JJG = IJKW(7,IW)
         KKG = IJKW(8,IW)
         IOR  = IJKW(4,IW)
         TMP_G = 0.5_EB*(TMP(IIG,JJG,KKG)+TMP_F(IW))
         ITMP = MIN(5000,NINT(TMP_G))         
         YY_GET=0._EB
         CALL GET_AVERAGE_SPECIFIC_HEAT(YY_GET,H_G_A,ITMP)
         YY_GET(N) = 1._EB
         CALL GET_AVERAGE_SPECIFIC_HEAT(YY_GET,H_G,ITMP)               
         HDIFF = (H_G-H_G_A)*TMP_G
         RHO_D_DYDN = 2._EB*RHODW(IW,N)*(YYP(IIG,JJG,KKG,N)-YY_F(IW,N))*RDN(IW)
         SELECT CASE(IOR)
            CASE( 1) 
               H_RHO_D_DYDX(IIG-1,JJG,KKG) = 0._EB
            CASE(-1) 
               H_RHO_D_DYDX(IIG,JJG,KKG)   = 0._EB
            CASE( 2) 
               H_RHO_D_DYDY(IIG,JJG-1,KKG) = 0._EB
            CASE(-2) 
               H_RHO_D_DYDY(IIG,JJG,KKG)   = 0._EB
            CASE( 3) 
               H_RHO_D_DYDZ(IIG,JJG,KKG-1) = 0._EB
            CASE(-3) 
               H_RHO_D_DYDZ(IIG,JJG,KKG)   = 0._EB
         END SELECT
         DP(IIG,JJG,KKG) = DP(IIG,JJG,KKG) - RDN(IW)*HDIFF*RHO_D_DYDN
      ENDDO WALL_LOOP2
      !$OMP END DO
 
      CYLINDER: SELECT CASE(CYLINDRICAL)
         CASE(.FALSE.) CYLINDER  ! 3D or 2D Cartesian Coords
            !$OMP DO COLLAPSE(3) PRIVATE(K,J,I)
            DO K=1,KBAR
               DO J=1,JBAR
                  DO I=1,IBAR
                     !!$ IF ((K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_DIVG_07'
                     DP(I,J,K) = DP(I,J,K) + (H_RHO_D_DYDX(I,J,K)-H_RHO_D_DYDX(I-1,J,K))*RDX(I) + &
                                             (H_RHO_D_DYDY(I,J,K)-H_RHO_D_DYDY(I,J-1,K))*RDY(J) + &
                                             (H_RHO_D_DYDZ(I,J,K)-H_RHO_D_DYDZ(I,J,K-1))*RDZ(K)
                  ENDDO
               ENDDO
            ENDDO
            !$OMP END DO

         CASE(.TRUE.) CYLINDER  ! 2D Cylindrical Coords
            !$OMP SINGLE
            J = 1
            !$OMP END SINGLE
            !$OMP DO COLLAPSE(2) PRIVATE(K,I) FIRSTPRIVATE(J)
            DO K=1,KBAR
               DO I=1,IBAR
                  !!$ IF ((K == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_DIVG_08'
                  DP(I,J,K) = DP(I,J,K) + (R(I)*H_RHO_D_DYDX(I,J,K)-R(I-1)*H_RHO_D_DYDX(I-1,J,K))*RDX(I)*RRN(I) + &
                                          (     H_RHO_D_DYDZ(I,J,K)-       H_RHO_D_DYDZ(I,J,K-1))*RDZ(K)
               ENDDO
            ENDDO
            !$OMP END DO
      END SELECT CYLINDER
   ENDIF SPECIES_DIFFUSION

   ! Compute del dot rho*D del Y_n or del dot rho D del Z
 
   CYLINDER2: SELECT CASE(CYLINDRICAL)
      CASE(.FALSE.) CYLINDER2  ! 3D or 2D Cartesian Coords
         !$OMP DO COLLAPSE(3) PRIVATE(K,J,I)
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  !!$ IF ((K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_DIVG_09'
                  DEL_RHO_D_DEL_Y(I,J,K,N) = DEL_RHO_D_DEL_Y(I,J,K,N) + (RHO_D_DYDX(I,J,K)-RHO_D_DYDX(I-1,J,K))*RDX(I) + &
                                                                        (RHO_D_DYDY(I,J,K)-RHO_D_DYDY(I,J-1,K))*RDY(J) + &
                                                                        (RHO_D_DYDZ(I,J,K)-RHO_D_DYDZ(I,J,K-1))*RDZ(K)
               ENDDO
            ENDDO
         ENDDO
         !$OMP END DO
      CASE(.TRUE.) CYLINDER2  ! 2D Cylindrical Coords
         !$OMP SINGLE
         J=1
         !$OMP END SINGLE
         !$OMP DO COLLAPSE(2) PRIVATE(K,I)
         DO K=1,KBAR
            DO I=1,IBAR
               !!$ IF ((K == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_DIVG_10'
               DEL_RHO_D_DEL_Y(I,J,K,N) = DEL_RHO_D_DEL_Y(I,J,K,N) + &
                                                              (R(I)*RHO_D_DYDX(I,J,K)-R(I-1)*RHO_D_DYDX(I-1,J,K))*RDX(I)*RRN(I) + &
                                                              (     RHO_D_DYDZ(I,J,K)-       RHO_D_DYDZ(I,J,K-1))*RDZ(K)
            ENDDO
         ENDDO
         !$OMP END DO
   END SELECT CYLINDER2
   
   ! Compute -Sum h_n del dot rho*D del Y_n
 
!   SPECIES_DIFFUSION_2: IF (.NOT.MIXTURE_FRACTION) THEN
   SPECIES_DIFFUSION_2: IF (SPECIES(N)%MODE/=MIXTURE_FRACTION_SPECIES) THEN
      !$OMP DO COLLAPSE(3) PRIVATE(K,J,I,ITMP,YY_GET,H_G_A,H_G,HDIFF)
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               !!$ IF ((K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_DIVG_11'               
               ITMP = MIN(5000,INT(TMP(I,J,K)))
               YY_GET=0._EB
               CALL GET_AVERAGE_SPECIFIC_HEAT(YY_GET,H_G_A,ITMP)
               YY_GET(N) = 1._EB
               CALL GET_AVERAGE_SPECIFIC_HEAT(YY_GET,H_G,ITMP)               
               HDIFF = (H_G-H_G_A)*TMP(I,J,K)
               DP(I,J,K) = DP(I,J,K) - HDIFF*DEL_RHO_D_DEL_Y(I,J,K,N)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO NOWAIT
   ENDIF SPECIES_DIFFUSION_2
   !$OMP END PARALLEL
ENDDO SPECIES_LOOP

! Compute del dot k del T
 
ENERGY: IF (.NOT.ISOTHERMAL .AND. .NOT.EVACUATION_ONLY(NM)) THEN
 
   KDTDX => WORK1
   KDTDY => WORK2
   KDTDZ => WORK3
   KP    => WORK4
   
   !$OMP PARALLEL SHARED(KP)
   
   !$OMP WORKSHARE
   KP = Y2K_C(NINT(TMPA))*SPECIES(0)%MW
   !$OMP END WORKSHARE
   
   ! Compute thermal conductivity k (KP)
 
   K_DNS_OR_LES: IF (DNS .AND. .NOT.EVACUATION_ONLY(NM)) THEN
      IF (N_SPECIES > 0 ) THEN
         !$OMP DO COLLAPSE(3) PRIVATE(K,J,I,ITMP,YY_GET) 
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  !!$ IF ((K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_DIVG_12'
                  IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
                  ITMP = MIN(5000,NINT(TMP(I,J,K)))
                  YY_GET(:) = YYP(I,J,K,:)
                  CALL GET_CONDUCTIVITY(YY_GET,KP(I,J,K),ITMP)    
               ENDDO
            ENDDO
         ENDDO
         !$OMP END DO
      ELSE
         !$OMP DO COLLAPSE(3) PRIVATE(K,J,I,ITMP)
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  !!$ IF ((K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_DIVG_13'
                  IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
                  ITMP = MIN(5000,NINT(TMP(I,J,K)))
                  KP(I,J,K) = Y2K_C(ITMP)*SPECIES(0)%MW
               ENDDO
            ENDDO
         ENDDO
         !$OMP END DO
      ENDIF
      !$OMP DO PRIVATE(IW,II,JJ,KK,IIG,JJG,KKG)
      BOUNDARY_LOOP: DO IW=1,NEWC
         II  = IJKW(1,IW)
         JJ  = IJKW(2,IW)
         KK  = IJKW(3,IW)
         IIG = IJKW(6,IW)
         JJG = IJKW(7,IW)
         KKG = IJKW(8,IW)
         KP(II,JJ,KK) = KP(IIG,JJG,KKG)
      ENDDO BOUNDARY_LOOP
      !$OMP END DO

   ELSE K_DNS_OR_LES
   
      !$OMP WORKSHARE
      KP = MU*CPOPR
      !$OMP END WORKSHARE
      
   ENDIF K_DNS_OR_LES

   ! Compute k*dT/dx, etc

   !$OMP DO COLLAPSE(3) PRIVATE(K,J,I,DTDX,DTDY,DTDZ)
   DO K=0,KBAR
      DO J=0,JBAR
         DO I=0,IBAR
            !!$ IF ((K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_DIVG_14'
            DTDX = (TMP(I+1,J,K)-TMP(I,J,K))*RDXN(I)
            KDTDX(I,J,K) = .5_EB*(KP(I+1,J,K)+KP(I,J,K))*DTDX
            DTDY = (TMP(I,J+1,K)-TMP(I,J,K))*RDYN(J)
            KDTDY(I,J,K) = .5_EB*(KP(I,J+1,K)+KP(I,J,K))*DTDY
            DTDZ = (TMP(I,J,K+1)-TMP(I,J,K))*RDZN(K)
            KDTDZ(I,J,K) = .5_EB*(KP(I,J,K+1)+KP(I,J,K))*DTDZ
         ENDDO
      ENDDO
   ENDDO
   !$OMP END DO

   ! Correct thermal gradient (k dT/dn) at boundaries

   !$OMP DO PRIVATE(IW,II,JJ,KK,IIG,JJG,KKG,IOR)
   CORRECTION_LOOP: DO IW=1,NWC
      IF (BOUNDARY_TYPE(IW)==NULL_BOUNDARY .OR. BOUNDARY_TYPE(IW)==POROUS_BOUNDARY) CYCLE CORRECTION_LOOP
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
   !$OMP END DO

   ! Compute (q + del dot k del T) and add to the divergence
 
   CYLINDER3: SELECT CASE(CYLINDRICAL)
      CASE(.FALSE.) CYLINDER3   ! 3D or 2D Cartesian
         !$OMP DO COLLAPSE(3) PRIVATE(K,J,I,DELKDELT)
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  !!$ IF ((K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_DIVG_15'
                  DELKDELT = (KDTDX(I,J,K)-KDTDX(I-1,J,K))*RDX(I) + &
                             (KDTDY(I,J,K)-KDTDY(I,J-1,K))*RDY(J) + &
                             (KDTDZ(I,J,K)-KDTDZ(I,J,K-1))*RDZ(K)
                  DP(I,J,K) = DP(I,J,K) + DELKDELT + Q(I,J,K) + QR(I,J,K)
               ENDDO 
            ENDDO
         ENDDO
         !$OMP END DO
      CASE(.TRUE.) CYLINDER3   ! 2D Cylindrical
         !$OMP DO COLLAPSE(3) PRIVATE(K,J,I,DELKDELT)
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  !!$ IF ((K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_DIVG_16'
                  DELKDELT = & 
                  (R(I)*KDTDX(I,J,K)-R(I-1)*KDTDX(I-1,J,K))*RDX(I)*RRN(I) + &
                       (KDTDZ(I,J,K)-       KDTDZ(I,J,K-1))*RDZ(K)
                  DP(I,J,K) = DP(I,J,K) + DELKDELT + Q(I,J,K) + QR(I,J,K)
               ENDDO 
            ENDDO
         ENDDO
         !$OMP END DO
   END SELECT CYLINDER3
   !$OMP END PARALLEL
 
ENDIF ENERGY

! Compute RTRM = R*sum(Y_i/M_i)/(PBAR*C_P) and multiply it by divergence terms already summed up
 
RTRM => WORK1

IF (N_SPECIES==0 .OR. EVACUATION_ONLY(NM)) THEN
  !$OMP PARALLEL DO PRIVATE(K,J,I,ITMP)
   DO K=1,KBAR
      IF (EVACUATION_ONLY(NM)) CYCLE
      DO J=1,JBAR
         DO I=1,IBAR
            !!$ IF ((K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_DIVG_17'
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
            ITMP = MIN(5000,NINT(TMP(I,J,K)))
            RTRM(I,J,K) = R_PBAR(K,PRESSURE_ZONE(I,J,K))*RSUM0/Y2CP_C(ITMP)
            DP(I,J,K) = RTRM(I,J,K)*DP(I,J,K)
         ENDDO
      ENDDO 
   ENDDO
   !$OMP END PARALLEL DO
ELSE
   !$OMP PARALLEL DO PRIVATE(K,J,I,ITMP,YY_GET,CP_MF)
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            !!$ IF ((K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_DIVG_18'
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
            ITMP = MIN(5000,NINT(TMP(I,J,K)))
            YY_GET(:) = YYP(I,J,K,:)
            CALL GET_SPECIFIC_HEAT(YY_GET,CP_MF,ITMP)
            RTRM(I,J,K) = R_PBAR(K,PRESSURE_ZONE(I,J,K))*RSUM(I,J,K)/CP_MF
            DP(I,J,K) = RTRM(I,J,K)*DP(I,J,K)
         ENDDO
      ENDDO 
   ENDDO
   !$OMP END PARALLEL DO
ENDIF 

! Compute (Wbar/rho) Sum (1/W_n) del dot rho*D del Y_n

!SPECIES_DIFFUSION_3: IF (.NOT.MIXTURE_FRACTION) THEN
   DO N=1,N_SPECIES
      IF (EVACUATION_ONLY(NM)) CYCLE
      IF (SPECIES(N)%MODE == MIXTURE_FRACTION_SPECIES) CYCLE
      !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(K,J,I)
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               !!$ IF ((N == 1) .AND. (K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_DIVG_19'
               DP(I,J,K) = DP(I,J,K) + (SPECIES(N)%RCON-SPECIES(0)%RCON)/(RSUM(I,J,K)*RHOP(I,J,K))*DEL_RHO_D_DEL_Y(I,J,K,N)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END PARALLEL DO
   ENDDO
!ENDIF SPECIES_DIFFUSION_3

! Add contribution of evaporating droplets
 
!$OMP PARALLEL
IF (NLP>0 .AND. N_EVAP_INDICES > 0 .AND. .NOT.EVACUATION_ONLY(NM)) THEN
   !$OMP DO COLLAPSE(3) PRIVATE(K,J,I)
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            !!$ IF ((K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_DIVG_20'
            DP(I,J,K) = DP(I,J,K) + D_LAGRANGIAN(I,J,K)
         ENDDO
      ENDDO
   ENDDO
   !$OMP END DO
ENDIF
 
! Atmospheric Stratification Term

IF (STRATIFICATION .AND. .NOT.EVACUATION_ONLY(NM)) THEN
   !$OMP DO COLLAPSE(3) PRIVATE(K,J,I)
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            !!$ IF ((K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_DIVG_21'
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
            DP(I,J,K) = DP(I,J,K) + (RTRM(I,J,K)-R_PBAR(K,PRESSURE_ZONE(I,J,K)))*0.5_EB*(W(I,J,K)+W(I,J,K-1))*GVEC(3)*RHO_0(K)
         ENDDO
      ENDDO
   ENDDO
   !$OMP END DO
ENDIF

! Compute normal component of velocity at boundaries, UWS

PREDICT_NORMALS: IF (PREDICTOR) THEN
 
   !$OMP WORKSHARE
   FDS_LEAK_AREA(:,:,NM) = 0._EB
   !$OMP END WORKSHARE

   !$OMP DO PRIVATE(IW,IOR,IBC,SF,IPZ,IOPZ,TSI,TIME_RAMP_FACTOR,DELTA_P,PRES_RAMP_FACTOR,IIG,JJG,KKG,II,JJ,KK) 
   WALL_LOOP3: DO IW=1,NWC
      !!$ IF ((IW == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_DIVG_22'
      IOR = IJKW(4,IW)
      WALL_CELL_TYPE: SELECT CASE (BOUNDARY_TYPE(IW))
         CASE (NULL_BOUNDARY)
            UWS(IW) = 0._EB
         CASE (SOLID_BOUNDARY,POROUS_BOUNDARY)
            IBC = IJKW(5,IW)
            SF => SURFACE(IBC)
            EVAC_IF_NOT: IF (.NOT.EVACUATION_ONLY(NM)) THEN
            IF (SF%SPECIES_BC_INDEX==SPECIFIED_MASS_FLUX .OR. SF%SPECIES_BC_INDEX==INTERPOLATED_BC .OR. &
                SF%SPECIES_BC_INDEX==HVAC_BOUNDARY) CYCLE WALL_LOOP3
            IF (SF%LEAK_PATH(1) == PRESSURE_ZONE_WALL(IW) .OR. SF%LEAK_PATH(2)==PRESSURE_ZONE_WALL(IW)) THEN
               IPZ  = PRESSURE_ZONE_WALL(IW)
               IF (IPZ/=SF%LEAK_PATH(1)) THEN
                  IOPZ = SF%LEAK_PATH(1)
               ELSE
                  IOPZ = SF%LEAK_PATH(2)
               ENDIF
               !$OMP ATOMIC
               FDS_LEAK_AREA(IPZ,IOPZ,NM) = FDS_LEAK_AREA(IPZ,IOPZ,NM) + AW(IW)
            ENDIF
            ENDIF EVAC_IF_NOT
            IF (TW(IW)==T_BEGIN .AND. SF%RAMP_INDEX(TIME_VELO)>=1) THEN
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
            IF (BOUNDARY_TYPE(IW)==POROUS_BOUNDARY .AND. IOR>0) UWS(IW) = -UWS(IW)  ! One-way flow through POROUS plate
            IF (EVACUATION_ONLY(NM)) UWS(IW) = TIME_RAMP_FACTOR*PRES_RAMP_FACTOR*UW0(IW)
            IF (SURFACE(IBC)%MASS_FLUX_TOTAL /= 0._EB) THEN
               IIG = IJKW(6,IW) 
               JJG = IJKW(7,IW) 
               KKG = IJKW(8,IW) 
               UWS(IW) = UWS(IW)*RHOA/RHO_F(IW)
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
   !$OMP END DO

   !$OMP WORKSHARE
   DUWDT(1:NEWC) = RDT*(UWS(1:NEWC)-UW(1:NEWC))
   !$OMP END WORKSHARE
   
ELSE PREDICT_NORMALS
   
   !$OMP WORKSHARE
   UW = UWS
   !$OMP END WORKSHARE

ENDIF PREDICT_NORMALS

!$OMP END PARALLEL
 
! Calculate pressure rise in each of the pressure zones by summing divergence expression over each zone

PRESSURE_ZONE_LOOP: DO IPZ=1,N_ZONE

   USUM(IPZ,NM) = 0._EB
   DSUM(IPZ,NM) = 0._EB
   PSUM(IPZ,NM) = 0._EB
   ZONE_VOLUME  = 0._EB

   IF (EVACUATION_ONLY(NM)) CYCLE PRESSURE_ZONE_LOOP
 
   !!$OMP PARALLEL DO COLLAPSE(3) PRIVATE(K,J,I,VC) REDUCTION(+:ZONE_VOLUME,DSUM,PSUM)
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
   !!$OMP END PARALLEL DO

   ! Calculate leakage and fan flows to all other pressure zones

   U_LEAK = 0._EB
   OTHER_ZONE_LOOP: DO IOPZ=0,N_ZONE
      IF (IOPZ == IPZ) CYCLE OTHER_ZONE_LOOP
      DELTA_P = PBAR(1,IPZ) - PBAR(1,IOPZ)
      IF (FDS_LEAK_AREA(IPZ,IOPZ,NM) > 0._EB) THEN
         VDOT_LEAK = SIGN(1._EB,DELTA_P)*MIN(ACTUAL_LEAK_AREA(IPZ,IOPZ)*SQRT(2._EB*ABS(DELTA_P)/RHOA), &
                                             ABS(DELTA_P)*RDT*ZONE_VOLUME/(GAMMA*PBAR(1,IPZ)))
         U_LEAK(IOPZ) = FDS_LEAK_AREA_RATIO(IPZ,IOPZ,NM)*VDOT_LEAK/FDS_LEAK_AREA(IPZ,IOPZ,NM)
      ENDIF
   ENDDO OTHER_ZONE_LOOP

   ! Calculate the volume flux to the boundary of the pressure zone (int u dot dA)

   WALL_LOOP4: DO IW=1,NWC
      IF (PRESSURE_ZONE_WALL(IW) /= IPZ)     CYCLE WALL_LOOP4
      IF (BOUNDARY_TYPE(IW)/=SOLID_BOUNDARY .AND. BOUNDARY_TYPE(IW)/=POROUS_BOUNDARY) CYCLE WALL_LOOP4
      IBC = IJKW(5,IW)
      SF=>SURFACE(IBC)
      IF ((SF%LEAK_PATH(1)==IPZ .OR. SF%LEAK_PATH(2)==IPZ) .AND. PREDICTOR) THEN
         IF (SF%LEAK_PATH(1)==IPZ) THEN
            IOPZ = SF%LEAK_PATH(2)
         ELSE
            IOPZ = SF%LEAK_PATH(1)
         ENDIF
         UWS(IW) = UWS(IW) + U_LEAK(IOPZ)
         IF (IW<=NEWC) DUWDT(IW) = RDT*(UWS(IW)-UW(IW))   ! DUWDT is needed by Poisson solver only at external boundary
      ENDIF
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
REAL(EB) :: RDT,TNOW,P_EQUILIBRIUM_NUM,P_EQUILIBRIUM_DEN,RF
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
   P_EQUILIBRIUM_NUM = PBAR(1,IPZ)*PSUM(IPZ,NM)
   P_EQUILIBRIUM_DEN = PSUM(IPZ,NM)
   DO IOPZ=N_ZONE,0,-1
      IF (IOPZ==IPZ) CYCLE
      IF (CONNECTED_ZONES(IPZ,IOPZ,NM)) THEN
         IF (IOPZ==0) THEN
            P_EQUILIBRIUM_NUM = P_0(1)
            P_EQUILIBRIUM_DEN = 1._EB
         ELSE
            P_EQUILIBRIUM_NUM = P_EQUILIBRIUM_NUM + PBAR(1,IOPZ)*PSUM(IOPZ,NM)
            P_EQUILIBRIUM_DEN = P_EQUILIBRIUM_DEN + PSUM(IOPZ,NM)
         ENDIF
      ENDIF
   ENDDO
   IF (P_EQUILIBRIUM_DEN /=0._EB) USUM_ADD(IPZ) = USUM_ADD(IPZ) + RDT*PSUM(IPZ,NM)*(PBAR(1,IPZ)-P_EQUILIBRIUM_NUM/P_EQUILIBRIUM_DEN)
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
               IF (ASUM(IOPZ,IOPZ2,NM)>0._EB .AND. USUM_ADD(IOPZ)/=0._EB) RF=MIN(RF,10._EB*ASUM(IOPZ,IOPZ2,NM)/ABS(USUM_ADD(IOPZ)))
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
 
   IF (PSUM(IPZ,NM)/=0._EB) THEN
      D_PBAR_DT_P(IPZ) = (DSUM(IPZ,NM) - USUM(IPZ,NM))/PSUM(IPZ,NM)
      P_ZONE(IPZ)%DPSTAR = D_PBAR_DT_P(IPZ)
   ENDIF

   ! Add pressure derivative to divergence

   !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(K,J,I)
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            !!$ IF ((K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_DIVG_25'
            IF (PRESSURE_ZONE(I,J,K) /= IPZ) CYCLE 
            IF (SOLID(CELL_INDEX(I,J,K)))    CYCLE
            DP(I,J,K) = DP(I,J,K) + (RTRM(I,J,K)-R_PBAR(K,IPZ))*D_PBAR_DT_P(IPZ)
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO

ENDDO PRESSURE_ZONE_LOOP

! Zero out divergence in solid cells

!$OMP PARALLEL 
!$OMP DO PRIVATE(IC,I,J,K) 
SOLID_LOOP: DO IC=1,CELL_COUNT
   !!$ IF ((IC == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_DIVG_26'
   IF (.NOT.SOLID(IC)) CYCLE SOLID_LOOP
   I = I_CELL(IC)
   J = J_CELL(IC)
   K = K_CELL(IC)
   DP(I,J,K) = 0._EB
ENDDO SOLID_LOOP
!$OMP END DO

! Specify divergence in boundary cells to account for volume being generated at the walls

!$OMP DO PRIVATE(IW,II,JJ,KK,IOR,IIG,JJG,KKG) 
BC_LOOP: DO IW=1,NWC
   !!$ IF ((IW == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_DIVG_27'
   IF (BOUNDARY_TYPE(IW)==NULL_BOUNDARY .OR. BOUNDARY_TYPE(IW)==POROUS_BOUNDARY) CYCLE BC_LOOP
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
!$OMP END DO


! Compute time derivative of the divergence, dD/dt

TRUE_PROJECTION: IF (PROJECTION) THEN

   !$OMP SINGLE
   DIV=>WORK1
   !$OMP END SINGLE
   IF (PREDICTOR) THEN
      !$OMP DO COLLAPSE(3) PRIVATE(K,J,I)
      DO K = 1,KBAR
         DO J = 1,JBAR
            DO I = 1,IBAR
               !!$ IF ((K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_DIVG_28'
               DIV(I,J,K) = (R(I)*U(I,J,K)-R(I-1)*U(I-1,J,K))*RDX(I)*RRN(I) + (V(I,J,K)-V(I,J-1,K))*RDY(J) + &
                            (W(I,J,K)-W(I,J,K-1))*RDZ(K)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO
      !$OMP WORKSHARE
      DDDT = (DP-DIV)*RDT
      !$OMP END WORKSHARE
   ELSEIF (CORRECTOR) THEN
      !$OMP DO COLLAPSE(3) PRIVATE(K,J,I)
      DO K = 1,KBAR
         DO J = 1,JBAR
            DO I = 1,IBAR
               !!$ IF ((K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_DIVG_29'
               DIV(I,J,K) = (R(I)*U(I,J,K) -R(I-1)*U(I-1,J,K)) *RDX(I)*RRN(I) + (V(I,J,K)- V(I,J-1,K)) *RDY(J) + &
                            (W(I,J,K) -W(I,J,K-1)) *RDZ(K) &
                          + (R(I)*US(I,J,K)-R(I-1)*US(I-1,J,K))*RDX(I)*RRN(I) + (VS(I,J,K)-VS(I,J-1,K))*RDY(J) + &
                            (WS(I,J,K)-WS(I,J,K-1))*RDZ(K)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO
      !$OMP WORKSHARE
      D = DDDT
      DDDT = (2._EB*DP-DIV)*RDT
      !$OMP END WORKSHARE
   ENDIF
   
ELSE TRUE_PROJECTION

   IF (PREDICTOR) THEN
      !$OMP WORKSHARE
      DDDT = (DS-D)*RDT
      !$OMP END WORKSHARE
   ELSE
      !$OMP SINGLE
      D_NEW => WORK1
      !$OMP END SINGLE
      !$OMP WORKSHARE
      D_NEW = DP
      DDDT  = (2._EB*D_NEW-DS-D)*RDT
      D     = D_NEW
      !$OMP END WORKSHARE
   ENDIF
   
   ! Adjust dD/dt to correct error in divergence due to velocity matching at interpolated boundaries
   
   !$OMP DO PRIVATE(IW,IIG,JJG,KKG)
   DO IW=1,NEWC
      IF (IJKW(9,IW)==0) CYCLE
      IIG = IJKW(6,IW)
      JJG = IJKW(7,IW)
      KKG = IJKW(8,IW)
      IF (PREDICTOR) DDDT(IIG,JJG,KKG) = DDDT(IIG,JJG,KKG) + DS_CORR(IW)*RDT
      IF (CORRECTOR) DDDT(IIG,JJG,KKG) = DDDT(IIG,JJG,KKG) + (2._EB*D_CORR(IW)-DS_CORR(IW))*RDT
   ENDDO
   !$OMP END DO
ENDIF TRUE_PROJECTION
!$OMP END PARALLEL
   
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
