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
USE PHYSICAL_FUNCTIONS, ONLY: GET_CP,GET_D,GET_K

! Compute contributions to the divergence term
 
INTEGER, INTENT(IN) :: NM
REAL(EB), POINTER, DIMENSION(:,:,:) :: KDTDX,KDTDY,KDTDZ,DP,KP, &
          RHO_D_DYDX,RHO_D_DYDY,RHO_D_DYDZ,RHO_D,RHOP,H_RHO_D_DYDX,H_RHO_D_DYDY,H_RHO_D_DYDZ,RTRM
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: YYP
REAL(EB) :: DELKDELT,VC,DTDX,DTDY,DTDZ,K_SUM,CP_SUM,TNOW, &
            HDIFF,DYDX,DYDY,DYDZ,T,RDT,RHO_D_DYDN,TSI,VDOT_LEAK,TIME_RAMP_FACTOR,ZONE_VOLUME,CP_MF,Z_2,DELTA_P,PRES_RAMP_FACTOR
TYPE(SURFACE_TYPE), POINTER :: SF
INTEGER :: IW,N,IOR,II,JJ,KK,IIG,JJG,KKG,ITMP,IBC,IZ,I,J,K,IPZ,IOPZ
 
IF (SOLID_PHASE_ONLY) RETURN

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
 
IF (N_SPECIES > 0) THEN
   RHO_D_DYDX  => WORK1
   RHO_D_DYDY  => WORK2
   RHO_D_DYDZ  => WORK3
   SELECT CASE(PREDICTOR)
      CASE(.TRUE.)  
         YYP => YYS 
      CASE(.FALSE.) 
         YYP => YY  
   END SELECT
ENDIF

! Zero out divergence to start
 
DP  = 0._EB
IF (N_SPECIES > 0) DEL_RHO_D_DEL_Y = 0._EB

! Add species diffusion terms to divergence expression and compute diffusion term for species equations
 
SPECIES_LOOP: DO N=1,N_SPECIES
 
   ! Compute rho*D
 
   RHO_D => WORK4
    
   IF (DNS .AND. SPECIES(N)%MODE==GAS_SPECIES) THEN
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               ITMP = 0.1_EB*TMP(I,J,K)
               RHO_D(I,J,K) = RHOP(I,J,K)*SPECIES(N)%D(ITMP)
            ENDDO 
         ENDDO
      ENDDO
   ENDIF
    
   IF (DNS .AND. SPECIES(N)%MODE==MIXTURE_FRACTION_SPECIES) THEN
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               ITMP = 0.1_EB*TMP(I,J,K)
               IZ   = NINT(Z_SUM(I,J,K)*100._EB)
               IZ   = MAX(0,MIN(IZ,100))
               Z_2 = 0._EB
               IF (CO_PRODUCTION) Z_2 = YYP(I,J,K,I_PROG_CO)
               CALL GET_D(YYP(I,J,K,I_FUEL),Z_2,YYP(I,J,K,I_PROG_F),Y_SUM(I,J,K),RHO_D(I,J,K),ITMP)
               RHO_D(I,J,K) = RHOP(I,J,K)*RHO_D(I,J,K)
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   IF (LES) RHO_D = MU*RSC
 
   ! Compute rho*D del Y
 
   DO K=0,KBAR
      DO J=0,JBAR
         DO I=0,IBAR
            DYDX = (YYP(I+1,J,K,N)-YYP(I,J,K,N))*RDXN(I)
            RHO_D_DYDX(I,J,K) = .5_EB*(RHO_D(I+1,J,K)+RHO_D(I,J,K))*DYDX
            DYDY = (YYP(I,J+1,K,N)-YYP(I,J,K,N))*RDYN(J)
            RHO_D_DYDY(I,J,K) = .5_EB*(RHO_D(I,J+1,K)+RHO_D(I,J,K))*DYDY
            DYDZ = (YYP(I,J,K+1,N)-YYP(I,J,K,N))*RDZN(K)
            RHO_D_DYDZ(I,J,K) = .5_EB*(RHO_D(I,J,K+1)+RHO_D(I,J,K))*DYDZ
         ENDDO
      ENDDO
   ENDDO
 
   ! Correct rho*D del Y at boundaries and store rho*D at boundaries
 
   WALL_LOOP: DO IW=1,NWC
      IF (BOUNDARY_TYPE(IW)==NULL_BOUNDARY .OR. BOUNDARY_TYPE(IW)==POROUS_BOUNDARY) CYCLE WALL_LOOP
      IIG = IJKW(6,IW) 
      JJG = IJKW(7,IW) 
      KKG = IJKW(8,IW) 
      RHODW(IW,N) = RHO_D(IIG,JJG,KKG)
      RHO_D_DYDN  = RHODW(IW,N)*(YYP(IIG,JJG,KKG,N)-YY_W(IW,N))*RDN(IW)
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
   ENDDO WALL_LOOP

   ! Compute del dot h_n*rho*D del Y_n only for non-mixture fraction cases
 
   SPECIES_DIFFUSION: IF (.NOT.MIXTURE_FRACTION) THEN

      H_RHO_D_DYDX => WORK4
      H_RHO_D_DYDY => WORK5
      H_RHO_D_DYDZ => WORK6
    
      DO K=0,KBAR
         DO J=0,JBAR
            DO I=0,IBAR
               ITMP = .05_EB*(TMP(I+1,J,K)+TMP(I,J,K))
               HDIFF = SPECIES(N)%H_G(ITMP)-SPECIES(0)%H_G(ITMP)
               H_RHO_D_DYDX(I,J,K) = HDIFF*RHO_D_DYDX(I,J,K)
               ITMP = .05_EB*(TMP(I,J+1,K)+TMP(I,J,K))
               HDIFF = SPECIES(N)%H_G(ITMP)-SPECIES(0)%H_G(ITMP)
               H_RHO_D_DYDY(I,J,K) = HDIFF*RHO_D_DYDY(I,J,K)
               ITMP = .05_EB*(TMP(I,J,K+1)+TMP(I,J,K))
               HDIFF = SPECIES(N)%H_G(ITMP)-SPECIES(0)%H_G(ITMP)
               H_RHO_D_DYDZ(I,J,K) = HDIFF*RHO_D_DYDZ(I,J,K)
            ENDDO
         ENDDO
      ENDDO

      WALL_LOOP2: DO IW=1,NWC
         IF (BOUNDARY_TYPE(IW)==NULL_BOUNDARY .OR. BOUNDARY_TYPE(IW)==POROUS_BOUNDARY) CYCLE WALL_LOOP2
         IIG = IJKW(6,IW)
         JJG = IJKW(7,IW)
         KKG = IJKW(8,IW)
         IOR  = IJKW(4,IW)
         ITMP = 0.1_EB*TMP_F(IW)
         HDIFF = SPECIES(N)%H_G(ITMP)-SPECIES(0)%H_G(ITMP)
         RHO_D_DYDN = RHODW(IW,N)*(YYP(IIG,JJG,KKG,N)-YY_W(IW,N))*RDN(IW)
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
 
      CYLINDER: SELECT CASE(CYLINDRICAL)
         CASE(.FALSE.) CYLINDER  ! 3D or 2D Cartesian Coords
            DO K=1,KBAR
               DO J=1,JBAR
                  DO I=1,IBAR
                     DP(I,J,K) = DP(I,J,K) + (H_RHO_D_DYDX(I,J,K)-H_RHO_D_DYDX(I-1,J,K))*RDX(I) + &
                                             (H_RHO_D_DYDY(I,J,K)-H_RHO_D_DYDY(I,J-1,K))*RDY(J) + &
                                             (H_RHO_D_DYDZ(I,J,K)-H_RHO_D_DYDZ(I,J,K-1))*RDZ(K)
                  ENDDO
               ENDDO
            ENDDO
         CASE(.TRUE.) CYLINDER  ! 2D Cylindrical Coords
            J = 1
            DO K=1,KBAR
               DO I=1,IBAR
                  DP(I,J,K) = DP(I,J,K) + (R(I)*H_RHO_D_DYDX(I,J,K)-R(I-1)*H_RHO_D_DYDX(I-1,J,K))*RDX(I)*RRN(I) + &
                                          (     H_RHO_D_DYDZ(I,J,K)-       H_RHO_D_DYDZ(I,J,K-1))*RDZ(K)
               ENDDO
            ENDDO
      END SELECT CYLINDER
    
   ENDIF SPECIES_DIFFUSION

   ! Compute del dot rho*D del Y_n or del dot rho D del Z
 
   CYLINDER2: SELECT CASE(CYLINDRICAL)
      CASE(.FALSE.) CYLINDER2  ! 3D or 2D Cartesian Coords
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  DEL_RHO_D_DEL_Y(I,J,K,N) = DEL_RHO_D_DEL_Y(I,J,K,N) + (RHO_D_DYDX(I,J,K)-RHO_D_DYDX(I-1,J,K))*RDX(I) + &
                                                                        (RHO_D_DYDY(I,J,K)-RHO_D_DYDY(I,J-1,K))*RDY(J) + &
                                                                        (RHO_D_DYDZ(I,J,K)-RHO_D_DYDZ(I,J,K-1))*RDZ(K)
               ENDDO
            ENDDO
         ENDDO
      CASE(.TRUE.) CYLINDER2  ! 2D Cylindrical Coords
         J=1
         DO K=1,KBAR
            DO I=1,IBAR
               DEL_RHO_D_DEL_Y(I,J,K,N) = DEL_RHO_D_DEL_Y(I,J,K,N) + &
                                                              (R(I)*RHO_D_DYDX(I,J,K)-R(I-1)*RHO_D_DYDX(I-1,J,K))*RDX(I)*RRN(I) + &
                                                              (     RHO_D_DYDZ(I,J,K)-       RHO_D_DYDZ(I,J,K-1))*RDZ(K)
            ENDDO
         ENDDO
   END SELECT CYLINDER2

   ! Compute -Sum h_n del dot rho*D del Y_n
 
   SPECIES_DIFFUSION_2: IF (.NOT.MIXTURE_FRACTION) THEN
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               ITMP = 0.1_EB*TMP(I,J,K)
               HDIFF = SPECIES(N)%H_G(ITMP)-SPECIES(0)%H_G(ITMP)
               DP(I,J,K) = DP(I,J,K) - HDIFF*DEL_RHO_D_DEL_Y(I,J,K,N)
            ENDDO
         ENDDO
      ENDDO
   ENDIF SPECIES_DIFFUSION_2
 
ENDDO SPECIES_LOOP
 

! Compute del dot k del T
 
ENERGY: IF (.NOT.ISOTHERMAL) THEN
 
   KDTDX => WORK1
   KDTDY => WORK2
   KDTDZ => WORK3
   KP    => WORK4
 
   ! Compute thermal conductivity k (KP)
 
   K_DNS_OR_LES: IF (DNS) THEN
 
      IF (.NOT.MIXTURE_FRACTION) THEN
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  ITMP = 0.1_EB*TMP(I,J,K)
                  K_SUM = SPECIES(0)%K(ITMP)
                  DO N=1,N_SPECIES
                     K_SUM = K_SUM + YYP(I,J,K,N)*(SPECIES(N)%K(ITMP)-SPECIES(0)%K(ITMP))
                  END DO
                  KP(I,J,K) = K_SUM
               ENDDO
            ENDDO
         ENDDO
      ELSE   ! MIXTURE_FRACTION TRUE
         Z_2 = 0._EB
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  ITMP = 0.1_EB*TMP(I,J,K)
                  IF(CO_PRODUCTION) THEN
                     CALL GET_CP(YY(I,J,K,I_FUEL),YY(I,J,K,I_PROG_CO),YY(I,J,K,I_PROG_F),Y_SUM(I,J,K),CP_MF,ITMP)
                     IF (N_SPECIES > 3) THEN
                        CP_SUM = 0._EB
                        DO N=1,N_SPECIES
                           IF (SPECIES(N)%MODE/=MIXTURE_FRACTION_SPECIES) &
                           CP_SUM = CP_SUM + YYP(I,J,K,N)*(SPECIES(N)%CP(ITMP)-SPECIES(0)%CP(ITMP))
                        END DO
                        CP_MF = CP_SUM+(1._EB-Y_SUM(I,J,K))*CP_MF
                     ENDIF
                  ELSE
                     CALL GET_CP(YY(I,J,K,I_FUEL),0._EB,YY(I,J,K,I_PROG_F),Y_SUM(I,J,K),CP_MF,ITMP)                  
                     IF (N_SPECIES > 2) THEN
                        CP_SUM = 0._EB
                        DO N=1,N_SPECIES
                           IF (SPECIES(N)%MODE/=MIXTURE_FRACTION_SPECIES) &
                           CP_SUM = CP_SUM + YYP(I,J,K,N)*(SPECIES(N)%CP(ITMP)-SPECIES(0)%CP(ITMP))
                        END DO
                        CP_MF = CP_SUM+(1._EB-Y_SUM(I,J,K))*CP_MF
                     ENDIF
                  ENDIF
                  KP(I,J,K) = MU(I,J,K)*CP_MF*RPR
               ENDDO
            ENDDO
         ENDDO
      ENDIF
      BOUNDARY_LOOP: DO IW=1,NEWC
         II  = IJKW(1,IW)
         JJ  = IJKW(2,IW)
         KK  = IJKW(3,IW)
         IIG = IJKW(6,IW)
         JJG = IJKW(7,IW)
         KKG = IJKW(8,IW)
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

! Compute RTRM = R*sum(Y_i/M_i)/(PBAR*C_P) and multiply it by divergence terms already summed up
 
RTRM => WORK1
 
IF (MIXTURE_FRACTION) THEN
   Z_2 = 0._EB
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            ITMP = 0.1_EB*TMP(I,J,K)
            IF(CO_PRODUCTION) THEN
               CALL GET_CP(YY(I,J,K,I_FUEL),YY(I,J,K,I_PROG_CO),YY(I,J,K,I_PROG_F),Y_SUM(I,J,K),CP_MF,ITMP)
               IF (N_SPECIES > 3) THEN
                  CP_SUM = 0._EB
                  DO N=1,N_SPECIES
                     IF (SPECIES(N)%MODE/=MIXTURE_FRACTION_SPECIES) &
                     CP_SUM = CP_SUM + YYP(I,J,K,N)*(SPECIES(N)%CP(ITMP)-SPECIES(0)%CP(ITMP))
                     CP_MF = CP_SUM+(1._EB-Y_SUM(I,J,K))*CP_MF
                  END DO
               ENDIF
            ELSE
               CALL GET_CP(YY(I,J,K,I_FUEL),Z_2,YY(I,J,K,I_PROG_F),Y_SUM(I,J,K),CP_MF,ITMP)  
               IF (N_SPECIES > 2) THEN
                  CP_SUM = 0._EB
                  DO N=1,N_SPECIES
                     IF (SPECIES(N)%MODE/=MIXTURE_FRACTION_SPECIES) &
                     CP_SUM = CP_SUM + YYP(I,J,K,N)*(SPECIES(N)%CP(ITMP)-SPECIES(0)%CP(ITMP))
                  END DO                                
                  CP_MF = CP_SUM+(1._EB-Y_SUM(I,J,K))*CP_MF
               ENDIF
            ENDIF
            RTRM(I,J,K) = R_PBAR(K,PRESSURE_ZONE(I,J,K))*RSUM(I,J,K)/CP_MF
            DP(I,J,K) = RTRM(I,J,K)*DP(I,J,K)
         ENDDO
      ENDDO 
   ENDDO
ENDIF
IF (.NOT.MIXTURE_FRACTION .AND. N_SPECIES>0) THEN
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            ITMP = 0.1_EB*TMP(I,J,K)
            CP_SUM = SPECIES(0)%CP(ITMP)
            DO N=1,N_SPECIES
               CP_SUM = CP_SUM + YYP(I,J,K,N)*(SPECIES(N)%CP(ITMP)-SPECIES(0)%CP(ITMP))
            END DO
            RTRM(I,J,K) = R_PBAR(K,PRESSURE_ZONE(I,J,K))*RSUM(I,J,K)/CP_SUM
            DP(I,J,K) = RTRM(I,J,K)*DP(I,J,K)
         ENDDO
      ENDDO
   ENDDO
ENDIF
IF (.NOT.MIXTURE_FRACTION .AND. N_SPECIES==0) THEN
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            ITMP = 0.1_EB*TMP(I,J,K)
            CP_SUM = SPECIES(0)%CP(ITMP)
            RTRM(I,J,K) = R_PBAR(K,PRESSURE_ZONE(I,J,K))*SPECIES(0)%RCON/CP_SUM
            DP(I,J,K) = RTRM(I,J,K)*DP(I,J,K)
         ENDDO
      ENDDO
   ENDDO
ENDIF

! Compute (Wbar/rho) Sum (1/W_n) del dot rho*D del Y_n

SPECIES_DIFFUSION_3: IF (.NOT.MIXTURE_FRACTION) THEN
   DO N=1,N_SPECIES
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               DP(I,J,K) = DP(I,J,K) + (SPECIES(N)%RCON-SPECIES(0)%RCON)/(RSUM(I,J,K)*RHOP(I,J,K))*DEL_RHO_D_DEL_Y(I,J,K,N)
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDIF SPECIES_DIFFUSION_3

! Add contribution of evaporating droplets
 
IF (NLP>0 .AND. N_EVAP_INDICIES > 0) THEN
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            DP(I,J,K) = DP(I,J,K) + D_VAP(I,J,K)
         ENDDO
      ENDDO
   ENDDO
ENDIF
 
! Atmospheric Stratification Term

IF (STRATIFICATION) THEN
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            DP(I,J,K) = DP(I,J,K) + (RTRM(I,J,K)-R_PBAR(K,PRESSURE_ZONE(I,J,K)))*0.5_EB*(W(I,J,K)+W(I,J,K-1))*GVEC(3)*RHO_0(K)
         ENDDO
      ENDDO
   ENDDO
ENDIF

! Compute normal component of velocity at boundaries, UWS

PREDICT_NORMALS: IF (PREDICTOR) THEN
 
   FDS_LEAK_AREA = 0._EB
 
   WALL_LOOP3: DO IW=1,NWC
      IOR = IJKW(4,IW)
      WALL_CELL_TYPE: SELECT CASE (BOUNDARY_TYPE(IW))
         CASE (NULL_BOUNDARY)
            UWS(IW) = 0._EB
         CASE (SOLID_BOUNDARY,POROUS_BOUNDARY)
            IBC = IJKW(5,IW)
            SF => SURFACE(IBC)
            IF (SF%SPECIES_BC_INDEX==SPECIFIED_MASS_FLUX .OR. SF%SPECIES_BC_INDEX==INTERPOLATED_BC) CYCLE WALL_LOOP3
            IF (SF%LEAK_PATH(1) == PRESSURE_ZONE_WALL(IW) .OR. SF%LEAK_PATH(2)==PRESSURE_ZONE_WALL(IW)) THEN
               IPZ  = PRESSURE_ZONE_WALL(IW)
               IF (IPZ/=SF%LEAK_PATH(1)) THEN
                  IOPZ = SF%LEAK_PATH(1)
               ELSE
                  IOPZ = SF%LEAK_PATH(2)
               ENDIF
               FDS_LEAK_AREA(IPZ,IOPZ) = FDS_LEAK_AREA(IPZ,IOPZ) + AW(IW)
            ENDIF
            IF (TW(IW)==T_BEGIN .AND. SF%RAMP_INDEX(TIME_VELO)>=1) THEN
               TSI = T            
            ELSE
               TSI = T+DT-TW(IW)
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
!            IF (SURFACE(IBC)%MASS_FLUX_TOTAL /= -999._EB .AND. UW0(IW) > 0._EB ) THEN
            IF (SURFACE(IBC)%MASS_FLUX_TOTAL /= -999._EB) THEN
               IIG = IJKW(6,IW) 
               JJG = IJKW(7,IW) 
               KKG = IJKW(8,IW) 
               UWS(IW) = UWS(IW)*2._EB*RHOA / (RHO_W(IW)+RHOP(IIG,JJG,KKG))
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

   DUWDT(1:NEWC) = RDT*(UWS(1:NEWC)-UW(1:NEWC))

ELSE PREDICT_NORMALS
   
   UW = UWS

ENDIF PREDICT_NORMALS
 
! Calculate pressure rise in each of the pressure zones by summing divergence expression over each zone

PRESSURE_ZONE_LOOP: DO IPZ=1,N_ZONE

   USUM(IPZ,NM) = 0._EB
   DSUM(IPZ,NM) = 0._EB
   PSUM(IPZ,NM) = 0._EB
   ZONE_VOLUME  = 0._EB
 
   DO K=1,KBAR
      DO J=1,JBAR
         INNER: DO I=1,IBAR
            IF (PRESSURE_ZONE(I,J,K) /= IPZ) CYCLE INNER
            IF (SOLID(CELL_INDEX(I,J,K)))           CYCLE INNER
            VC   = DX(I)*RC(I)*DY(J)*DZ(K)
            ZONE_VOLUME = ZONE_VOLUME + VC
            DSUM(IPZ,NM) = DSUM(IPZ,NM) + VC*DP(I,J,K)
            PSUM(IPZ,NM) = PSUM(IPZ,NM) + VC*(R_PBAR(K,IPZ)-RTRM(I,J,K))
         ENDDO INNER
      ENDDO
   ENDDO

   ! Calculate leakage and fan flows to all other pressure zones

   U_LEAK = 0._EB
   OTHER_ZONE_LOOP: DO IOPZ=0,N_ZONE
      IF (IOPZ == IPZ) CYCLE OTHER_ZONE_LOOP
      DELTA_P = PBAR(1,IPZ) - PBAR(1,IOPZ)
      IF (FDS_LEAK_AREA(IPZ,IOPZ) > 0._EB) THEN
         VDOT_LEAK = LEAK_AREA(IPZ,IOPZ)*SIGN(1._EB,DELTA_P)*SQRT(2._EB*ABS(DELTA_P)/RHOA)
         U_LEAK(IOPZ) = VDOT_LEAK/FDS_LEAK_AREA(IPZ,IOPZ)
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
REAL(EB), POINTER, DIMENSION(:,:,:) :: DP,D_OLD,RTRM
REAL(EB) :: RDT,TNOW
REAL(EB), POINTER, DIMENSION(:) :: D_PBAR_DT_P
INTEGER :: IW,IOR,II,JJ,KK,IIG,JJG,KKG,IC,I,J,K,IPZ

IF (SOLID_PHASE_ONLY) RETURN

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

PRESSURE_ZONE_LOOP: DO IPZ=1,N_ZONE

   IF (PREDICTOR) D_PBAR_DT_P => D_PBAR_S_DT
   IF (CORRECTOR) D_PBAR_DT_P => D_PBAR_DT

   ! Compute change in background pressure
 
   D_PBAR_DT_P(IPZ) = (DSUM(IPZ,NM) - USUM(IPZ,NM))/PSUM(IPZ,NM)

   ! Add pressure derivative to divergence

   DO K=1,KBAR
      DO J=1,JBAR
         INNER2: DO I=1,IBAR
            IF (PRESSURE_ZONE(I,J,K) /= IPZ) CYCLE INNER2
            DP(I,J,K) = DP(I,J,K) + (RTRM(I,J,K)-R_PBAR(K,IPZ))*D_PBAR_DT_P(IPZ)
         ENDDO INNER2
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
 
BC_LOOP: DO IW=1,NWC
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

! Compute time derivative of the divergence, dD/dt
 
IF (PREDICTOR) THEN
   DDDT = (DS-D)*RDT
ELSE
   D_OLD => WORK1
   D_OLD = DP
   DDDT  = (2._EB*DP-DS-D)*RDT
   D     = D_OLD
ENDIF
 
! Adjust dD/dt to correct error in divergence due to velocity matching at interpolated boundaries

IF (NMESHES>1) THEN
   DO IW=1,NEWC
      IF (BOUNDARY_TYPE(IW)/=INTERPOLATED_BOUNDARY) CYCLE
      IIG = IJKW(6,IW)
      JJG = IJKW(7,IW)
      KKG = IJKW(8,IW)
      IF (PREDICTOR) DDDT(IIG,JJG,KKG) = DDDT(IIG,JJG,KKG) + DS_CORR(IW)*RDT
      IF (CORRECTOR) DDDT(IIG,JJG,KKG) = DDDT(IIG,JJG,KKG) + (2._EB*D_CORR(IW)-DS_CORR(IW))*RDT
   ENDDO
ENDIF

TUSED(2,NM)=TUSED(2,NM)+SECOND()-TNOW
END SUBROUTINE DIVERGENCE_PART_2
 
 
SUBROUTINE CHECK_DIVERGENCE(NM)
USE COMP_FUNCTIONS, ONLY: SECOND 
! Computes maximum velocity divergence
 
INTEGER, INTENT(IN) :: NM
INTEGER  :: I,J,K
REAL(EB) :: DIV,RES,TNOW
 
TNOW=SECOND()
CALL POINT_TO_MESH(NM)
 
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
               DIV = (U(I,J,K)-U(I-1,J,K))*RDX(I) + &
                     (V(I,J,K)-V(I,J-1,K))*RDY(J) + &
                     (W(I,J,K)-W(I,J,K-1))*RDZ(K)
            CASE(.TRUE.)
               DIV = (R(I)*U(I,J,K)-R(I-1)*U(I-1,J,K))*RDX(I)*RRN(I) +  &
                     (W(I,J,K)-W(I,J,K-1))*RDZ(K)
         END SELECT
         RES = ABS(DIV-D(I,J,K))
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
