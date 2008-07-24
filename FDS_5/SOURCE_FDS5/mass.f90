MODULE MASS
 
! Compute the mass equation differences 
 
USE PRECISION_PARAMETERS
USE MESH_POINTERS

IMPLICIT NONE
PRIVATE
CHARACTER(255), PARAMETER :: massid='$Id$'
CHARACTER(255), PARAMETER :: massrev='$Revision$'
CHARACTER(255), PARAMETER :: massdate='$Date$'

REAL(EB), POINTER, DIMENSION(:,:,:,:) :: YYP
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW,RHOP,DP

PUBLIC MASS_FINITE_DIFFERENCES,DENSITY,GET_REV_mass,DENSITY_TVD
 
 
CONTAINS
 
SUBROUTINE MASS_FINITE_DIFFERENCES(NM)
USE COMP_FUNCTIONS, ONLY: SECOND
USE GLOBAL_CONSTANTS, ONLY: N_SPECIES,ISOTHERMAL,NULL_BOUNDARY,POROUS_BOUNDARY,PREDICTOR,CORRECTOR,EVACUATION_ONLY, &
                            SOLID_PHASE_ONLY,TUSED
INTEGER, INTENT(IN) :: NM
REAL(EB) :: FXYZ,PMDT,UDRHODN,TNOW
INTEGER  :: I,J,K,N,II,JJ,KK,IIG,JJG,KKG,IW,IOR
REAL(EB), POINTER, DIMENSION(:) :: UWP
REAL(EB), POINTER, DIMENSION(:,:,:) :: UDRHODX,VDRHODY,WDRHODZ,EPSX,EPSY,EPSZ
 
IF (EVACUATION_ONLY(NM)) RETURN
IF (SOLID_PHASE_ONLY) RETURN

TNOW=SECOND()
CALL POINT_TO_MESH(NM)
 
IF (PREDICTOR) THEN
   UU => U
   VV => V
   WW => W
   DP => D
   RHOP => RHO
   UWP  => UW
   PMDT = DT
ELSE
   UU => US
   VV => VS
   WW => WS
   DP => DS
   RHOP => RHOS
   UWP  => UWS
   PMDT = -DT
ENDIF

! Define local CFL numbers
 
EPSX => WORK1
EPSY => WORK2
EPSZ => WORK3

DO K=0,KBAR
   DO J=0,JBAR
      DO I=0,IBAR
         EPSX(I,J,K) = PMDT*UU(I,J,K)*RDXN(I)
         EPSY(I,J,K) = PMDT*VV(I,J,K)*RDYN(J)
         EPSZ(I,J,K) = PMDT*WW(I,J,K)*RDZN(K)
      ENDDO
   ENDDO
ENDDO

! Compute spatial differences for density equation
 
NOT_ISOTHERMAL_IF: IF (.NOT.ISOTHERMAL) THEN
 
   UDRHODX => WORK4
   VDRHODY => WORK5
   WDRHODZ => WORK6
   
   DO K=0,KBAR
      DO J=0,JBAR
         DO I=0,IBAR
            UDRHODX(I,J,K) = UU(I,J,K)*(RHOP(I+1,J,K)-RHOP(I,J,K))*RDXN(I)
            VDRHODY(I,J,K) = VV(I,J,K)*(RHOP(I,J+1,K)-RHOP(I,J,K))*RDYN(J)
            WDRHODZ(I,J,K) = WW(I,J,K)*(RHOP(I,J,K+1)-RHOP(I,J,K))*RDZN(K)
         ENDDO
      ENDDO
   ENDDO
   
   WLOOP: DO IW=1,NWC
      IF (BOUNDARY_TYPE(IW)==NULL_BOUNDARY .OR. BOUNDARY_TYPE(IW)==POROUS_BOUNDARY) CYCLE WLOOP
      II  = IJKW(1,IW) 
      IIG = IJKW(6,IW)
      JJ  = IJKW(2,IW) 
      JJG = IJKW(7,IW)
      KK  = IJKW(3,IW) 
      KKG = IJKW(8,IW)
      IOR = IJKW(4,IW)
      UDRHODN = UWP(IW)*(RHO_W(IW)-RHOP(IIG,JJG,KKG))*RDN(IW)
      SELECT CASE(IOR)
         CASE( 1)
            UDRHODX(II,JJ,KK)   = UDRHODN
         CASE(-1) 
            UDRHODX(II-1,JJ,KK) = UDRHODN
         CASE( 2) 
            VDRHODY(II,JJ,KK)   = UDRHODN
         CASE(-2) 
            VDRHODY(II,JJ-1,KK) = UDRHODN
         CASE( 3) 
            WDRHODZ(II,JJ,KK)   = UDRHODN
         CASE(-3) 
            WDRHODZ(II,JJ,KK-1) = UDRHODN
      END SELECT
   ENDDO WLOOP
   
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
            FXYZ   = .5_EB*(UDRHODX(I,J,K)  *(1._EB-EPSX(I,J,K))   +  &
                            UDRHODX(I-1,J,K)*(1._EB+EPSX(I-1,J,K)) +  &
                            VDRHODY(I,J,K)  *(1._EB-EPSY(I,J,K))   +  &
                            VDRHODY(I,J-1,K)*(1._EB+EPSY(I,J-1,K)) +  &
                            WDRHODZ(I,J,K)  *(1._EB-EPSZ(I,J,K))   +  &
                            WDRHODZ(I,J,K-1)*(1._EB+EPSZ(I,J,K-1)) )
            FRHO(I,J,K) = FXYZ + RHOP(I,J,K)*DP(I,J,K)
         ENDDO
      ENDDO
   ENDDO
 
ENDIF NOT_ISOTHERMAL_IF
 
! Compute the species equation differences
 
IF (N_SPECIES > 0) THEN
   IF (PREDICTOR) YYP => YY
   IF (CORRECTOR) YYP => YYS
   UDRHODX => WORK4
   VDRHODY => WORK5
   WDRHODZ => WORK6
ENDIF
 
SPECIES_LOOP: DO N=1,N_SPECIES

   DO K=0,KBAR
      DO J=0,JBAR
         DO I=0,IBAR
            UDRHODX(I,J,K) = UU(I,J,K)*( RHOP(I+1,J,K)*YYP(I+1,J,K,N)-RHOP(I,J,K)*YYP(I,J,K,N) )*RDXN(I)
            VDRHODY(I,J,K) = VV(I,J,K)*( RHOP(I,J+1,K)*YYP(I,J+1,K,N)-RHOP(I,J,K)*YYP(I,J,K,N) )*RDYN(J)
            WDRHODZ(I,J,K) = WW(I,J,K)*( RHOP(I,J,K+1)*YYP(I,J,K+1,N)-RHOP(I,J,K)*YYP(I,J,K,N) )*RDZN(K)
         ENDDO
      ENDDO
   ENDDO
 
   ! Correct U d(RHO*Y)/dx etc. on boundaries
 
   WLOOP2: DO IW=1,NWC
      IF (BOUNDARY_TYPE(IW)==NULL_BOUNDARY .OR. BOUNDARY_TYPE(IW)==POROUS_BOUNDARY) CYCLE WLOOP2
      II  = IJKW(1,IW) 
      IIG = IJKW(6,IW)
      JJ  = IJKW(2,IW) 
      JJG = IJKW(7,IW)
      KK  = IJKW(3,IW) 
      KKG = IJKW(8,IW)
      IOR = IJKW(4,IW)
      UDRHODN = UWP(IW)*( RHO_W(IW)*YY_W(IW,N) - RHOP(IIG,JJG,KKG)*YYP(IIG,JJG,KKG,N) )*RDN(IW)
      SELECT CASE(IOR)
         CASE( 1)
            UDRHODX(II,JJ,KK)   = UDRHODN
         CASE(-1)
            UDRHODX(II-1,JJ,KK) = UDRHODN
         CASE( 2)
            VDRHODY(II,JJ,KK)   = UDRHODN
         CASE(-2) 
            VDRHODY(II,JJ-1,KK) = UDRHODN
         CASE( 3) 
            WDRHODZ(II,JJ,KK)   = UDRHODN
         CASE(-3) 
            WDRHODZ(II,JJ,KK-1) = UDRHODN
      END SELECT
   ENDDO WLOOP2
 
  ! Sum up the convective and diffusive terms in the transport equation and store in DEL_RHO_D_DEL_Y
 
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            FXYZ   = .5_EB*(UDRHODX(I,J,K)  *(1._EB-EPSX(I,J,K))   +  &
                            UDRHODX(I-1,J,K)*(1._EB+EPSX(I-1,J,K)) +  &
                            VDRHODY(I,J,K)  *(1._EB-EPSY(I,J,K))   +  &
                            VDRHODY(I,J-1,K)*(1._EB+EPSY(I,J-1,K)) +  &
                            WDRHODZ(I,J,K)  *(1._EB-EPSZ(I,J,K))   +  &
                            WDRHODZ(I,J,K-1)*(1._EB+EPSZ(I,J,K-1)) ) 
            DEL_RHO_D_DEL_Y(I,J,K,N) = -DEL_RHO_D_DEL_Y(I,J,K,N) + FXYZ + RHOP(I,J,K)*YYP(I,J,K,N)*DP(I,J,K) 
         ENDDO
      ENDDO
   ENDDO
 
ENDDO SPECIES_LOOP
 
TUSED(3,NM)=TUSED(3,NM)+SECOND()-TNOW
END SUBROUTINE MASS_FINITE_DIFFERENCES
 
 
SUBROUTINE DENSITY(NM)

! Update the density and species mass fractions

USE COMP_FUNCTIONS, ONLY: SECOND 
USE PHYSICAL_FUNCTIONS, ONLY : GET_MOLECULAR_WEIGHT
USE GLOBAL_CONSTANTS, ONLY: N_SPECIES,CO_PRODUCTION,I_PROG_F,I_PROG_CO,I_FUEL,TMPMAX,TMPMIN,EVACUATION_ONLY,PREDICTOR,CORRECTOR, &
                            CHANGE_TIME_STEP,ISOTHERMAL,TMPA,N_SPEC_DILUENTS, N_ZONE,MIXTURE_FRACTION_SPECIES, &
                            GAS_SPECIES, MIXTURE_FRACTION,R0,SOLID_PHASE_ONLY,TUSED,BAROCLINIC,RHO_LOWER,RHO_UPPER
 
REAL(EB) :: WFAC,DTRATIO,OMDTRATIO,Z_2,TNOW
INTEGER  :: I,J,K,N
INTEGER, INTENT(IN) :: NM
REAL(EB), POINTER, DIMENSION(:,:,:) :: R_SUM_DILUENTS
 
IF (EVACUATION_ONLY(NM)) RETURN
IF (SOLID_PHASE_ONLY) RETURN

TNOW=SECOND()
CALL POINT_TO_MESH(NM)

PREDICTOR_STEP: SELECT CASE (PREDICTOR)

CASE(.TRUE.) PREDICTOR_STEP

   IF (.NOT.CHANGE_TIME_STEP(NM)) THEN

      DO N=1,N_SPECIES
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
                  YYS(I,J,K,N) = RHO(I,J,K)*YY(I,J,K,N) - DT*DEL_RHO_D_DEL_Y(I,J,K,N)
               ENDDO
            ENDDO
         ENDDO
      ENDDO

   ELSE

      DTRATIO   = DT/DTOLD
      OMDTRATIO = 1._EB - DTRATIO
      DO N=1,N_SPECIES
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
                  YYS(I,J,K,N) = OMDTRATIO*RHO(I,J,K) *YY(I,J,K,N) + DTRATIO*RHOS(I,J,K)*YYS(I,J,K,N)
               ENDDO
           ENDDO
         ENDDO
      ENDDO

   ENDIF

   ! Predict the density at the next time step (RHOS or RHO^*)

   IF (.NOT.ISOTHERMAL) THEN

      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               RHOS(I,J,K) = RHO(I,J,K)-DT*FRHO(I,J,K)
            ENDDO
         ENDDO
      ENDDO

   ELSE

      DO K=0,KBP1
         DO J=0,JBP1
            DO I=0,IBP1
               RHOS(I,J,K) = PBAR_S(K,PRESSURE_ZONE(I,J,K))/(TMPA*SPECIES(0)%RCON)
            ENDDO
         ENDDO
      ENDDO

      DO N=1,N_SPECIES
         WFAC = 1._EB - SPECIES(N)%RCON/SPECIES(0)%RCON
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  RHOS(I,J,K) = RHOS(I,J,K) + WFAC*YYS(I,J,K,N)
               ENDDO
            ENDDO
         ENDDO
      ENDDO

   ENDIF
 
   ! Correct densities above or below clip limits

   CALL CHECK_DENSITY
   
   IF (BAROCLINIC) THEN
      ! Set global min and max for rho
      RHO_LOWER = MIN(RHO_LOWER,MINVAL_GASPHASE(RHOS))
      RHO_UPPER = MAX(RHO_UPPER,MAXVAL_GASPHASE(RHOS))
   ENDIF

   ! Extract mass fraction from RHO * YY

   DO N=1,N_SPECIES
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               YYS(I,J,K,N) = YYS(I,J,K,N)/RHOS(I,J,K)
            ENDDO
         ENDDO
      ENDDO
   ENDDO

   ! Correct mass fractions above or below clip limits

   CALL CHECK_MASS_FRACTION

   ! Predict background pressure at next time step

   DO I=1,N_ZONE
      PBAR_S(:,I) = PBAR(:,I) + D_PBAR_DT(I)*DT
   ENDDO

   ! Compute mixture fraction and diluent sums: Y_SUM=Sum(Y_i), Z_SUM=Sum(Z_i)

   IF (MIXTURE_FRACTION) THEN
      Z_SUM  =  0._EB
      Y_SUM  =  0._EB
      IF (N_SPEC_DILUENTS > 0) THEN
         R_SUM_DILUENTS => WORK4
         R_SUM_DILUENTS =  0._EB
      ENDIF
      DO N=1,N_SPECIES
         IF (SPECIES(N)%MODE==MIXTURE_FRACTION_SPECIES) Z_SUM = Z_SUM + YYS(:,:,:,N)
         IF (SPECIES(N)%MODE==GAS_SPECIES) THEN
            Y_SUM = Y_SUM + YYS(:,:,:,N)
            R_SUM_DILUENTS(:,:,:) = R_SUM_DILUENTS(:,:,:) + SPECIES(N)%RCON*YYS(:,:,:,N)
         ENDIF
      ENDDO
   ENDIF

   ! Compute molecular weight term RSUM=R0*SUM(Y_i/M_i)
 
   IF (N_SPECIES>0 .AND. .NOT.MIXTURE_FRACTION) THEN
      RSUM = SPECIES(0)%RCON
      DO N=1,N_SPECIES
         WFAC = SPECIES(N)%RCON - SPECIES(0)%RCON
         RSUM(:,:,:) = RSUM(:,:,:) + WFAC*YYS(:,:,:,N)
      ENDDO
      IF (ISOTHERMAL) THEN
         DO K=0,KBP1
            DO J=0,JBP1
               DO I=0,IBP1
                  RHOS(I,J,K) = PBAR_S(K,PRESSURE_ZONE(I,J,K))/(TMPA*RSUM(I,J,K))
               ENDDO
            ENDDO
         ENDDO
      ENDIF
   ENDIF

   IF (MIXTURE_FRACTION) THEN
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (CO_PRODUCTION) THEN
                  Z_2 = YYS(I,J,K,I_PROG_CO)
               ElSE
                  Z_2 = 0._EB
               ENDIF
               CALL GET_MOLECULAR_WEIGHT(YYS(I,J,K,I_FUEL),Z_2,YYS(I,J,K,I_PROG_F),Y_SUM(I,J,K),RSUM(I,J,K))
               RSUM(I,J,K) = R0/RSUM(I,J,K)
            ENDDO
         ENDDO
      ENDDO
      IF (N_SPEC_DILUENTS > 0) RSUM = RSUM*(1._EB-Y_SUM) + R_SUM_DILUENTS
   ENDIF

   ! Extract predicted temperature at next time step from Equation of State

   IF (.NOT.ISOTHERMAL) THEN
      IF (N_SPECIES==0) THEN
         DO K=0,KBP1
            DO J=0,JBP1
               DO I=0,IBP1
                  TMP(I,J,K) = PBAR_S(K,PRESSURE_ZONE(I,J,K))/(SPECIES(0)%RCON*RHOS(I,J,K))
               ENDDO
            ENDDO
         ENDDO
      ELSE
         DO K=0,KBP1
            DO J=0,JBP1
               DO I=0,IBP1
                  TMP(I,J,K) = PBAR_S(K,PRESSURE_ZONE(I,J,K))/(RSUM(I,J,K)*RHOS(I,J,K))
               ENDDO
            ENDDO
         ENDDO
      ENDIF
      TMP = MAX(TMPMIN,MIN(TMPMAX,TMP))
   ENDIF

! The CORRECTOR step
   
CASE(.FALSE.) PREDICTOR_STEP

   ! Correct species mass fraction at next time step (YY here actually means YY*RHO)

   DO N=1,N_SPECIES
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
               YY(I,J,K,N) = .5_EB*(RHO(I,J,K)*YY(I,J,K,N) + RHOS(I,J,K)*YYS(I,J,K,N) - DT*DEL_RHO_D_DEL_Y(I,J,K,N) ) 
            ENDDO
         ENDDO
      ENDDO
   ENDDO

   ! Correct density at next time step

   IF (.NOT.ISOTHERMAL) THEN
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               RHO(I,J,K) = .5_EB*(RHO(I,J,K)+RHOS(I,J,K)-DT*FRHO(I,J,K))
            ENDDO
         ENDDO
      ENDDO
   ELSE
      DO K=0,KBP1
         DO J=0,JBP1
            DO I=0,IBP1
               RHO(I,J,K) = PBAR(K,PRESSURE_ZONE(I,J,K))/(SPECIES(0)%RCON*TMPA)
            ENDDO
         ENDDO
      ENDDO

      DO N=1,N_SPECIES
         WFAC = 1._EB - SPECIES(N)%RCON/SPECIES(0)%RCON
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  RHO(I,J,K) = RHO(I,J,K) + WFAC*YY(I,J,K,N)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   ! Correct densities above or below clip limits

   CALL CHECK_DENSITY
   
   IF (BAROCLINIC) THEN
      ! Set global min and max for rho
      RHO_LOWER = MIN(RHO_LOWER,MINVAL_GASPHASE(RHO))
      RHO_UPPER = MAX(RHO_UPPER,MAXVAL_GASPHASE(RHO))
   ENDIF
 
   ! Extract Y_n from rho*Y_n

   DO N=1,N_SPECIES
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               YY(I,J,K,N) = YY(I,J,K,N)/RHO(I,J,K)
            ENDDO
         ENDDO
      ENDDO
   ENDDO

   ! Correct mass fractions above or below clip limits

   CALL CHECK_MASS_FRACTION

   ! Correct background pressure

   DO I=1,N_ZONE
      PBAR(:,I) = .5_EB*(PBAR(:,I) + PBAR_S(:,I) + D_PBAR_S_DT(I)*DT)
   ENDDO
 
   ! Compute mixture fraction and diluent sums: Y_SUM=Sum(Y_i), Z_SUM=Sum(Z_i)

   IF (MIXTURE_FRACTION) THEN
      Z_SUM  =  0._EB
      Y_SUM  =  0._EB
      IF (N_SPEC_DILUENTS > 0) THEN
         R_SUM_DILUENTS => WORK4
         R_SUM_DILUENTS =  0._EB
         ENDIF
      DO N=1,N_SPECIES
         IF (SPECIES(N)%MODE==MIXTURE_FRACTION_SPECIES) Z_SUM = Z_SUM + YY(:,:,:,N)
         IF (SPECIES(N)%MODE==GAS_SPECIES) THEN
            Y_SUM = Y_SUM + YY(:,:,:,N)
            R_SUM_DILUENTS(:,:,:) = R_SUM_DILUENTS(:,:,:) + SPECIES(N)%RCON*YY(:,:,:,N)
         ENDIF
      ENDDO
   ENDIF

   ! Compute molecular weight term RSUM=R0*SUM(Y_i/M_i)
 
   IF (N_SPECIES>0 .AND. .NOT. MIXTURE_FRACTION) THEN
      RSUM = SPECIES(0)%RCON
      DO N=1,N_SPECIES
         WFAC = SPECIES(N)%RCON - SPECIES(0)%RCON
         RSUM(:,:,:) = RSUM(:,:,:) + WFAC*YY(:,:,:,N)
      ENDDO
      IF (ISOTHERMAL) THEN
         DO K=0,KBP1
            DO J=0,JBP1
               DO I=0,IBP1
                  RHO(I,J,K) = PBAR(K,PRESSURE_ZONE(I,J,K))/(TMPA*RSUM(I,J,K))
               ENDDO
            ENDDO
         ENDDO
      ENDIF
   ENDIF

   IF (MIXTURE_FRACTION) THEN
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (CO_PRODUCTION) THEN
                  Z_2 = YY(I,J,K,I_PROG_CO)
               ElSE
                  Z_2 = 0._EB
               ENDIF
               CALL GET_MOLECULAR_WEIGHT(YY(I,J,K,I_FUEL),Z_2,YY(I,J,K,I_PROG_F),Y_SUM(I,J,K),RSUM(I,J,K))
               RSUM(I,J,K) = R0/RSUM(I,J,K)
            ENDDO
         ENDDO
      ENDDO
      IF (N_SPEC_DILUENTS > 0) RSUM = RSUM*(1._EB-Y_SUM) + R_SUM_DILUENTS
   ENDIF

   ! Extract temperature from the Equation of State

   IF (.NOT.ISOTHERMAL) THEN
      IF (N_SPECIES==0) THEN
         DO K=0,KBP1
            DO J=0,JBP1
               DO I=0,IBP1
                  TMP(I,J,K) = PBAR(K,PRESSURE_ZONE(I,J,K))/(SPECIES(0)%RCON*RHO(I,J,K))
               ENDDO
            ENDDO
         ENDDO
      ELSE
         DO K=0,KBP1
            DO J=0,JBP1
               DO I=0,IBP1
                  TMP(I,J,K) = PBAR(K,PRESSURE_ZONE(I,J,K))/(RSUM(I,J,K)*RHO(I,J,K))
               ENDDO
            ENDDO
         ENDDO
      ENDIF
      TMP = MAX(TMPMIN,MIN(TMPMAX,TMP))
   ENDIF

END SELECT PREDICTOR_STEP

TUSED(3,NM)=TUSED(3,NM)+SECOND()-TNOW
 
CONTAINS
 

SUBROUTINE CHECK_DENSITY
 
! Redistribute mass from cells below or above the density cut-off limits

USE GLOBAL_CONSTANTS, ONLY : PREDICTOR, CORRECTOR, N_SPECIES,RHOMIN,RHOMAX 
REAL(EB) :: SUM,CONST,RHOMI,RHOPI,RHOMJ,RHOPJ,RHOMK,RHOPK,RHO00,RMIN,RMAX
INTEGER  :: IC,ISUM
LOGICAL :: LC(-3:3)
REAL(EB), POINTER, DIMENSION(:,:,:) :: RHODELTA

RHODELTA => WORK2

IF (PREDICTOR) THEN
   RHOP=>RHOS
   YYP=>YYS
ELSE
   RHOP=>RHO
   YYP=>YY
ENDIF
 
! Correct undershoots

RHODELTA = 0._EB

DO K=1,KBAR
   DO J=1,JBAR
      CHECK_LOOP: DO I=1,IBAR
         IC = CELL_INDEX(I,J,K)
         IF (SOLID(IC)) CYCLE CHECK_LOOP
         RMIN = RHOMIN
         IF (RHOP(I,J,K)>=RMIN) CYCLE CHECK_LOOP
         SUM   = 0.
         ISUM  = 0
         LC    = .FALSE.
         RHO00 = RHOP(I,J,K)
         RHOMI = RHOP(I-1,J,K)
         RHOPI = RHOP(I+1,J,K)
         RHOMJ = RHOP(I,J-1,K)
         RHOPJ = RHOP(I,J+1,K)
         RHOMK = RHOP(I,J,K-1)
         RHOPK = RHOP(I,J,K+1)
         IF (WALL_INDEX(IC,-1)==0 .AND. RHOMI>RMIN) THEN
            SUM = SUM + RHOMI
            ISUM = ISUM + 1
            LC(-1) = .TRUE.
         ENDIF
         IF (WALL_INDEX(IC, 1)==0 .AND. RHOPI>RMIN) THEN
            SUM = SUM + RHOPI
            ISUM = ISUM + 1
            LC( 1) = .TRUE.
         ENDIF
         IF (WALL_INDEX(IC,-2)==0 .AND. RHOMJ>RMIN) THEN
            SUM = SUM + RHOMJ
            ISUM = ISUM + 1
            LC(-2) = .TRUE.
         ENDIF
         IF (WALL_INDEX(IC, 2)==0 .AND. RHOPJ>RMIN) THEN
            SUM = SUM + RHOPJ
            ISUM = ISUM + 1
            LC( 2) = .TRUE.
         ENDIF
         IF (WALL_INDEX(IC,-3)==0 .AND. RHOMK>RMIN) THEN
            SUM = SUM + RHOMK
            ISUM = ISUM + 1
            LC(-3) = .TRUE.
         ENDIF
         IF (WALL_INDEX(IC, 3)==0 .AND. RHOPK>RMIN) THEN
            SUM = SUM + RHOPK
            ISUM = ISUM + 1
            LC( 3) = .TRUE.
         ENDIF
         IF (ISUM==0) THEN
            RHODELTA(I,J,K) = RMIN - RHOP(I,J,K)
            CYCLE CHECK_LOOP
         ELSE
            CONST = (RHOMIN-RHO00)/(SUM-ISUM*RHO00)
            IF (LC(-1)) RHODELTA(I-1,J,K) = RHODELTA(I-1,J,K) + MAX(RMIN,RHOMI+CONST*(RHO00-RHOMI)) - RHOP(I-1,J,K)
            IF (LC( 1)) RHODELTA(I+1,J,K) = RHODELTA(I+1,J,K) + MAX(RMIN,RHOPI+CONST*(RHO00-RHOPI)) - RHOP(I+1,J,K)
            IF (LC(-2)) RHODELTA(I,J-1,K) = RHODELTA(I,J-1,K) + MAX(RMIN,RHOMJ+CONST*(RHO00-RHOMJ)) - RHOP(I,J-1,K)
            IF (LC( 2)) RHODELTA(I,J+1,K) = RHODELTA(I,J+1,K) + MAX(RMIN,RHOPJ+CONST*(RHO00-RHOPJ)) - RHOP(I,J+1,K)
            IF (LC(-3)) RHODELTA(I,J,K-1) = RHODELTA(I,J,K-1) + MAX(RMIN,RHOMK+CONST*(RHO00-RHOMK)) - RHOP(I,J,K-1)
            IF (LC( 3)) RHODELTA(I,J,K+1) = RHODELTA(I,J,K+1) + MAX(RMIN,RHOPK+CONST*(RHO00-RHOPK)) - RHOP(I,J,K+1)
            RHODELTA(I,J,K) = RHODELTA(I,J,K) + RMIN - RHOP(I,J,K)
         ENDIF
      ENDDO CHECK_LOOP
   ENDDO
ENDDO

RHOP = MAX(RHOMIN,RHOP+RHODELTA)

! Correct overshoots

RHODELTA = 0._EB
DO K=1,KBAR
   DO J=1,JBAR
      CHECK_LOOP2: DO I=1,IBAR
         IC = CELL_INDEX(I,J,K)
         IF (SOLID(IC)) CYCLE CHECK_LOOP2
         RMAX = RHOMAX
         IF (RHOP(I,J,K)<=RMAX) CYCLE CHECK_LOOP2
         SUM   = 0.
         ISUM  = 0
         LC    = .FALSE.
         RHO00 = RHOP(I,J,K)
         RHOMI = RHOP(I-1,J,K)
         RHOPI = RHOP(I+1,J,K)
         RHOMJ = RHOP(I,J-1,K)
         RHOPJ = RHOP(I,J+1,K)
         RHOMK = RHOP(I,J,K-1)
         RHOPK = RHOP(I,J,K+1)
         IF (WALL_INDEX(IC,-1)==0 .AND. RHOMI<RMAX) THEN
            SUM = SUM + RHOMI
            ISUM = ISUM + 1
            LC(-1) = .TRUE.
         ENDIF
         IF (WALL_INDEX(IC, 1)==0 .AND. RHOPI<RMAX) THEN
            SUM = SUM + RHOPI
            ISUM = ISUM + 1
            LC( 1) = .TRUE.
         ENDIF
         IF (WALL_INDEX(IC,-2)==0 .AND. RHOMJ<RMAX) THEN
            SUM = SUM + RHOMJ
            ISUM = ISUM + 1
            LC(-2) = .TRUE.
         ENDIF
         IF (WALL_INDEX(IC, 2)==0 .AND. RHOPJ<RMAX) THEN
            SUM = SUM + RHOPJ
            ISUM = ISUM + 1
            LC( 2) = .TRUE.
         ENDIF
         IF (WALL_INDEX(IC,-3)==0 .AND. RHOMK<RMAX) THEN
            SUM = SUM + RHOMK
            ISUM = ISUM + 1
            LC(-3) = .TRUE.
         ENDIF
         IF (WALL_INDEX(IC, 3)==0 .AND. RHOPK<RMAX) THEN
            SUM = SUM + RHOPK
            ISUM = ISUM + 1
            LC( 3) = .TRUE.
         ENDIF
         IF (ISUM==0) THEN
            RHODELTA(I,J,K) = RMAX - RHOP(I,J,K)
            CYCLE CHECK_LOOP2
         ELSE
            CONST = (RMAX-RHO00)/(SUM-ISUM*RHO00)
            IF (LC(-1)) RHODELTA(I-1,J,K) = RHODELTA(I-1,J,K) + MIN(RMAX,RHOMI+CONST*(RHO00-RHOMI)) - RHOP(I-1,J,K)
            IF (LC( 1)) RHODELTA(I+1,J,K) = RHODELTA(I+1,J,K) + MIN(RMAX,RHOPI+CONST*(RHO00-RHOPI)) - RHOP(I+1,J,K)
            IF (LC(-2)) RHODELTA(I,J-1,K) = RHODELTA(I,J-1,K) + MIN(RMAX,RHOMJ+CONST*(RHO00-RHOMJ)) - RHOP(I,J-1,K)
            IF (LC( 2)) RHODELTA(I,J+1,K) = RHODELTA(I,J+1,K) + MIN(RMAX,RHOPJ+CONST*(RHO00-RHOPJ)) - RHOP(I,J+1,K)
            IF (LC(-3)) RHODELTA(I,J,K-1) = RHODELTA(I,J,K-1) + MIN(RMAX,RHOMK+CONST*(RHO00-RHOMK)) - RHOP(I,J,K-1)
            IF (LC( 3)) RHODELTA(I,J,K+1) = RHODELTA(I,J,K+1) + MIN(RMAX,RHOPK+CONST*(RHO00-RHOPK)) - RHOP(I,J,K+1)
            RHODELTA(I,J,K) = RHODELTA(I,J,K) + RMAX - RHOP(I,J,K)
         ENDIF
      ENDDO CHECK_LOOP2
   ENDDO
ENDDO

RHOP = MIN(RHOMAX,RHOP+RHODELTA)

END SUBROUTINE CHECK_DENSITY
 
 
SUBROUTINE CHECK_MASS_FRACTION

! Redistribute species mass from cells below or above the cut-off limits

USE GLOBAL_CONSTANTS, ONLY : PREDICTOR, CORRECTOR, N_SPECIES,YYMIN,YYMAX
REAL(EB) :: SUM,CONST,RHYMI,RHYPI,RHYMJ,RHYPJ,RHYMK,RHYPK,RHY0,YMI,YPI,YMJ,YPJ,YMK,YPK,Y00,YMIN,YMAX
INTEGER  :: IC,N,ISUM, IW_A(-3:3)
LOGICAL  :: LC(-3:3)
REAL(EB), POINTER, DIMENSION(:,:,:) :: YYDELTA

YYDELTA => WORK1
IF (PREDICTOR) THEN
   RHOP    => RHOS
   YYP     => YYS
ELSE
   RHOP    => RHO
   YYP     => YY
ENDIF

! Search the domain for negative values of Y or Z. Redistribute mass where appropriate.

SPECIESLOOP: DO N=1,N_SPECIES
   YYDELTA = 0._EB

   ! Do undershoots

   DO K=1,KBAR
      DO J=1,JBAR
         CHECK_LOOP: DO I=1,IBAR
            IC = CELL_INDEX(I,J,K)
            IF (SOLID(IC)) CYCLE CHECK_LOOP
            IW_A = WALL_INDEX(IC,:)
            Y00   = YYP(I,J,K,N)
            SUM   = 0._EB
            ISUM  = 0
            LC    = .FALSE.
            YMIN  = YYMAX(N) 
            IF (IW_A(-1) == 0) THEN
               YMI = YYP(I-1,J,K,N)
               LC(-1) = .TRUE.
            ELSE
               YMI = YY_W(IW_A(-1),N)  
            ENDIF          
            IF (IW_A( 1) == 0) THEN
               YPI = YYP(I+1,J,K,N)
               LC( 1) = .TRUE.
            ELSE
               YPI = YY_W(IW_A( 1),N)  
            ENDIF           
            IF (IW_A(-2) == 0) THEN
               YMJ = YYP(I,J-1,K,N)
               LC(-2) = .TRUE.
            ELSE
               YMJ = YY_W(IW_A(-2),N)  
            ENDIF         
            IF (IW_A( 2) == 0) THEN
               YPJ = YYP(I,J+1,K,N)
               LC( 2) = .TRUE.
            ELSE
               YPJ = YY_W(IW_A( 2),N)  
            ENDIF         
            IF (IW_A(-3) == 0) THEN
               YMK = YYP(I,J,K-1,N)
               LC(-3) = .TRUE.
            ELSE
               YMK = YY_W(IW_A(-3),N)  
            ENDIF         
            IF (IW_A( 3) == 0) THEN
               YPK = YYP(I,J,K+1,N)
               LC( 3) = .TRUE.
            ELSE
               YPK = YY_W(IW_A( 3),N)  
            ENDIF           
            YMIN  = MIN(YMI,YPI,YMJ,YPJ,YMK,YPK)
            YMIN = MAX(YMIN,YYMIN(N))
            IF ((DEL_RHO_D_DEL_Y(I,J,K,N) > 0._EB .AND. Y00 < YMIN) .OR. Y00 < YYMIN(N)) THEN
               RHY0  = RHOP(I,J,K)  *(YMIN - Y00)
               IF (LC(-1) .AND. YMI>YMIN) THEN! .AND. DEL_RHO_D_DEL_Y(I-1,J,K,N) < 0._EB) THEN
                  RHYMI = RHOP(I-1,J,K)*(YMI - YMIN)
                  SUM  = SUM + RHYMI 
                  ISUM = ISUM + 1
               ELSE
                  LC(-1) = .FALSE.
               ENDIF
               IF (LC( 1) .AND. YPI>YMIN) THEN! .AND. DEL_RHO_D_DEL_Y(I+1,J,K,N) < 0._EB) THEN
                  RHYPI = RHOP(I+1,J,K)*(YPI - YMIN)
                  SUM  = SUM + RHYPI
                  ISUM = ISUM + 1
               ELSE
                  LC( 1) = .FALSE.
               ENDIF
               IF (LC(-2) .AND. YMJ>YMIN) THEN! .AND. DEL_RHO_D_DEL_Y(I,J-1,K,N) < 0._EB) THEN
                  RHYMJ = RHOP(I,J-1,K)*(YMJ - YMIN)
                  SUM  = SUM + RHYMJ
                  ISUM = ISUM + 1
               ELSE
                  LC(-2) = .FALSE.
               ENDIF
               IF (LC( 2) .AND. YPJ>YMIN) THEN! .AND. DEL_RHO_D_DEL_Y(I,J+1,K,N) < 0._EB) THEN
                  RHYPJ = RHOP(I,J+1,K)*(YPJ - YMIN)
                  SUM  = SUM + RHYPJ
                  ISUM = ISUM + 1
                  LC( 2) = .TRUE.
               ELSE
                  LC( 2) = .FALSE.
               ENDIF
               IF (LC(-3) .AND. YMK>YMIN) THEN! .AND. DEL_RHO_D_DEL_Y(I,J,K-1,N) < 0._EB) THEN
               RHYMK = RHOP(I,J,K-1)*(YMK - YMIN)
                  SUM  = SUM + RHYMK
                  ISUM = ISUM + 1
               ELSE
                  LC(-3) = .FALSE.
               ENDIF
               IF (LC( 3) .AND. YPK>YMIN) THEN! .AND. DEL_RHO_D_DEL_Y(I,J,K+1,N) < 0._EB) THEN
                  RHYPK = RHOP(I,J,K+1)*(YPK - YMIN)
                  SUM  = SUM + RHYPK
                  ISUM = ISUM + 1
               ELSE
                  LC( 3) = .FALSE.
               ENDIF                
               IF (ISUM==0) THEN
                  IF (YMIN <= YYMIN(N)) YYDELTA(I,J,K) = YYDELTA(I,J,K) + YMIN - Y00  
                  CYCLE CHECK_LOOP
               ELSE
                  YYDELTA(I,J,K) = YYDELTA(I,J,K) + YMIN - Y00
               
                  CONST = MIN(1._EB,RHY0/SUM)
                  IF (LC(-1)) YYDELTA(I-1,J,K) = YYDELTA(I-1,J,K) - RHYMI*CONST/RHOP(I-1,J,K)
                  IF (LC( 1)) YYDELTA(I+1,J,K) = YYDELTA(I+1,J,K) - RHYPI*CONST/RHOP(I+1,J,K)
                  IF (LC(-2)) YYDELTA(I,J-1,K) = YYDELTA(I,J-1,K) - RHYMJ*CONST/RHOP(I,J-1,K)
                  IF (LC( 2)) YYDELTA(I,J+1,K) = YYDELTA(I,J+1,K) - RHYPJ*CONST/RHOP(I,J+1,K)
                  IF (LC(-3)) YYDELTA(I,J,K-1) = YYDELTA(I,J,K-1) - RHYMK*CONST/RHOP(I,J,K-1)
                  IF (LC( 3)) YYDELTA(I,J,K+1) = YYDELTA(I,J,K+1) - RHYPK*CONST/RHOP(I,J,K+1)
               ENDIF
            ENDIF
         ENDDO CHECK_LOOP
      ENDDO
   ENDDO
   YYP(:,:,:,N) = YYP(:,:,:,N) + YYDELTA
   YYDELTA=0._EB

   ! Do overshoots

   DO K=1,KBAR
      DO J=1,JBAR
         CHECK_LOOP2: DO I=1,IBAR
            IC = CELL_INDEX(I,J,K)
            IF (SOLID(IC)) CYCLE CHECK_LOOP2
            IW_A  = WALL_INDEX(IC,:)
            Y00   = YYP(I,J,K,N)
            SUM   = 0._EB
            ISUM  = 0
            LC    = .FALSE.
            YMIN  = YYMAX(N) 
            IF (IW_A(-1) == 0) THEN
               YMI = YYP(I-1,J,K,N)
               LC(-1) = .TRUE.
            ELSE
               YMI = YY_W(IW_A(-1),N)  
            ENDIF          
            IF (IW_A( 1) == 0) THEN
               YPI = YYP(I+1,J,K,N)
               LC( 1) = .TRUE.
            ELSE
               YPI = YY_W(IW_A( 1),N)  
            ENDIF           
            IF (IW_A(-2) == 0) THEN
               YMJ = YYP(I,J-1,K,N)
               LC(-2) = .TRUE.
            ELSE
               YMJ = YY_W(IW_A(-2),N)  
            ENDIF         
            IF (IW_A( 2) == 0) THEN
               YPJ = YYP(I,J+1,K,N)
               LC( 2) = .TRUE.
            ELSE
               YPJ = YY_W(IW_A( 2),N)  
            ENDIF         
            IF (IW_A(-3) == 0) THEN
               YMK = YYP(I,J,K-1,N)
               LC(-3) = .TRUE.
            ELSE
               YMK = YY_W(IW_A(-3),N)  
            ENDIF         
            IF (IW_A( 3) == 0) THEN
               YPK = YYP(I,J,K+1,N)
               LC( 3) = .TRUE.
            ELSE
               YPK = YY_W(IW_A( 3),N)  
            ENDIF           
            YMAX  = MAX(YMI,YPI,YMJ,YPJ,YMK,YPK)
            YMAX = MIN(YMAX,YYMAX(N))            
            IF ((DEL_RHO_D_DEL_Y(I,J,K,N) < 0._EB .AND. Y00 > YMAX) .OR. Y00 > YYMAX(N)) THEN
               RHY0  = RHOP(I,J,K)  *(Y00 - YMAX)
               IF (LC(-1) .AND. YMI<YMAX) THEN! .AND. DEL_RHO_D_DEL_Y(I-1,J,K,N) > 0._EB) THEN
                  RHYMI = RHOP(I-1,J,K)*(YMAX - YMI)
                  SUM  = SUM + RHYMI
                  ISUM = ISUM + 1
               ELSE
                  LC(-1) = .FALSE.
               ENDIF
               IF (LC( 1) .AND. YPI<YMAX) THEN! .AND. DEL_RHO_D_DEL_Y(I+1,J,K,N) > 0._EB) THEN
                  RHYPI = RHOP(I+1,J,K)*(YMAX - YPI)
                  SUM  = SUM + RHYPI
                  ISUM = ISUM + 1
               ELSE
                  LC( 1) = .FALSE.
               ENDIF
               IF (LC(-2) .AND. YMJ<YMAX) THEN! .AND. DEL_RHO_D_DEL_Y(I,J-1,K,N) > 0._EB) THEN
                  RHYMJ = RHOP(I,J-1,K)*(YMAX - YMJ)
                  SUM  = SUM + RHYMJ
                  ISUM = ISUM + 1
               ELSE
                  LC(-2) = .FALSE.
               ENDIF
               IF (LC( 2) .AND. YPJ<YMAX) THEN! .AND. DEL_RHO_D_DEL_Y(I,J+1,K,N) > 0._EB) THEN
                  RHYPJ = RHOP(I,J+1,K)*(YMAX - YPJ)
                  SUM  = SUM + RHYPJ
                  ISUM = ISUM + 1
               ELSE
                  LC( 2) = .FALSE.
               ENDIF
               IF (LC(-3) .AND. YMK<YMAX) THEN! .AND. DEL_RHO_D_DEL_Y(I,J,K-1,N) > 0._EB) THEN
                  RHYMK = RHOP(I,J,K-1)*(YMAX - YMK)
                  SUM  = SUM + RHYMK
                  ISUM = ISUM + 1
               ELSE
                  LC(-3) = .FALSE.
               ENDIF
               IF (LC( 3) .AND. YPK<YMAX) THEN! .AND. DEL_RHO_D_DEL_Y(I,J,K+1,N) > 0._EB) THEN
                  RHYPK = RHOP(I,J,K+1)*(YMAX - YPK)
                  SUM  = SUM + RHYPK
                  ISUM = ISUM + 1
               ELSE
                  LC( 3) = .FALSE.
               ENDIF                      
               IF (ISUM==0) THEN
                  IF(YMAX >= YYMAX(N)) YYDELTA(I,J,K) = YYDELTA(I,J,K) + YMAX - Y00
                  CYCLE CHECK_LOOP2
               ELSE
                  YYDELTA(I,J,K) = YYDELTA(I,J,K) + YMAX - Y00               
                  CONST = MIN(1._EB,RHY0/SUM)
                  IF (LC(-1)) YYDELTA(I-1,J,K) = YYDELTA(I-1,J,K) + RHYMI*CONST/RHOP(I-1,J,K)
                  IF (LC( 1)) YYDELTA(I+1,J,K) = YYDELTA(I+1,J,K) + RHYPI*CONST/RHOP(I+1,J,K)
                  IF (LC(-2)) YYDELTA(I,J-1,K) = YYDELTA(I,J-1,K) + RHYMJ*CONST/RHOP(I,J-1,K)
                  IF (LC( 2)) YYDELTA(I,J+1,K) = YYDELTA(I,J+1,K) + RHYPJ*CONST/RHOP(I,J+1,K)
                  IF (LC(-3)) YYDELTA(I,J,K-1) = YYDELTA(I,J,K-1) + RHYMK*CONST/RHOP(I,J,K-1)
                  IF (LC( 3)) YYDELTA(I,J,K+1) = YYDELTA(I,J,K+1) + RHYPK*CONST/RHOP(I,J,K+1)
               ENDIF
            ENDIF
         ENDDO CHECK_LOOP2
      ENDDO
   ENDDO   

   YYP(:,:,:,N) = YYP(:,:,:,N) + YYDELTA   

ENDDO SPECIESLOOP

RETURN

END SUBROUTINE CHECK_MASS_FRACTION
 
END SUBROUTINE DENSITY


!===========================================================================
! The following are experimental scalar transport routines which are invoked
! by setting FLUX_LIMITER = {0,1,2,3,4} on the MISC line.  See the function
! SCALAR_FACE_VALUE below for a description of the FLUX_LIMITER value. ~RJM
!===========================================================================

SUBROUTINE DENSITY_TVD(NM)

! Update the density and species mass fractions

USE COMP_FUNCTIONS, ONLY: SECOND 
USE PHYSICAL_FUNCTIONS, ONLY : GET_MOLECULAR_WEIGHT
USE GLOBAL_CONSTANTS, ONLY: N_SPECIES,CO_PRODUCTION,I_PROG_F,I_PROG_CO,I_FUEL,TMPMAX,TMPMIN,EVACUATION_ONLY,PREDICTOR,CORRECTOR, &
                            CHANGE_TIME_STEP,ISOTHERMAL,TMPA,N_SPEC_DILUENTS, N_ZONE,MIXTURE_FRACTION_SPECIES, LU_ERR, &
                            GAS_SPECIES, MIXTURE_FRACTION,R0,SOLID_PHASE_ONLY,TUSED,FLUX_LIMITER,BTOL,RHO_LOWER,RHO_UPPER,BAROCLINIC
 
REAL(EB) :: WFAC,Z_2,TNOW
!REAL(EB) :: DSMIN,DSMAX
INTEGER  :: I,J,K,N
INTEGER, INTENT(IN) :: NM
REAL(EB), POINTER, DIMENSION(:,:,:) :: R_SUM_DILUENTS

REAL(EB), POINTER, DIMENSION(:,:,:) :: RHON
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: YYN,RHOYYP


 
IF (EVACUATION_ONLY(NM)) RETURN
IF (SOLID_PHASE_ONLY) RETURN

TNOW=SECOND()
CALL POINT_TO_MESH(NM)

CALL SCALARF ! Computes FRHOYY and FRHO and populates SCALAR_SAVE3

RHON => SCALAR_SAVE1
IF (N_SPECIES>0) THEN
   YYN    => SCALAR_SAVE2
   RHOYYP => SCALAR_SAVE3
ENDIF

SELECT_SUBSTEP: IF (PREDICTOR) THEN
   
   IF (N_SPECIES>0) YYN = RHOYYP
   RHON = RHO
   
   ! Update the density
   
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            RHOS(I,J,K) = RHON(I,J,K) - DT*FRHO(I,J,K)
            RHOS(I,J,K) = MAX(0._EB,RHOS(I,J,K))
         ENDDO
      ENDDO
   ENDDO
   
   IF (BAROCLINIC) THEN
      ! Set global min and max for rho
      RHO_LOWER = MIN(RHO_LOWER,MINVAL_GASPHASE(RHOS))
      RHO_UPPER = MAX(RHO_UPPER,MAXVAL_GASPHASE(RHOS))
   ENDIF
   
   ! Update mass fractions
    
   DO N=1,N_SPECIES
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR             
               YYS(I,J,K,N) = YYN(I,J,K,N) - DT*FRHOYY(I,J,K,N)
            ENDDO
         ENDDO
      ENDDO
   ENDDO

  
   ! Extract REALIZABLE YY from REALIZABLE RHO*YY
   
   DO N=1,N_SPECIES
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               YYS(I,J,K,N) = YYS(I,J,K,N)/RHOS(I,J,K)
               YYS(I,J,K,N) = MIN(1._EB,YYS(I,J,K,N))
               YYS(I,J,K,N) = MAX(0._EB,YYS(I,J,K,N))
            ENDDO
         ENDDO
      ENDDO
   ENDDO

   ! Predict background pressure at next time step

   DO I=1,N_ZONE
      PBAR_S(:,I) = PBAR(:,I) + D_PBAR_DT(I)*DT
   ENDDO

   ! Compute mixture fraction and diluent sums: Y_SUM=Sum(Y_i), Z_SUM=Sum(Z_i)

   IF (MIXTURE_FRACTION) THEN
      Z_SUM  =  0._EB
      Y_SUM  =  0._EB
      IF (N_SPEC_DILUENTS > 0) THEN
         R_SUM_DILUENTS => WORK4
         R_SUM_DILUENTS =  0._EB
      ENDIF
      DO N=1,N_SPECIES
         IF (SPECIES(N)%MODE==MIXTURE_FRACTION_SPECIES) Z_SUM = Z_SUM + YYS(:,:,:,N)
         IF (SPECIES(N)%MODE==GAS_SPECIES) THEN
            Y_SUM = Y_SUM + YYS(:,:,:,N)
            R_SUM_DILUENTS(:,:,:) = R_SUM_DILUENTS(:,:,:) + SPECIES(N)%RCON*YYS(:,:,:,N)
         ENDIF
      ENDDO
   ENDIF

   ! Compute molecular weight term RSUM=R0*SUM(Y_i/M_i)
 
   IF (N_SPECIES>0 .AND. .NOT.MIXTURE_FRACTION) THEN
      RSUM = SPECIES(0)%RCON
      DO N=1,N_SPECIES
         WFAC = SPECIES(N)%RCON - SPECIES(0)%RCON
         RSUM(:,:,:) = RSUM(:,:,:) + WFAC*YYS(:,:,:,N)
      ENDDO
      IF (ISOTHERMAL) THEN
         DO K=0,KBP1
            DO J=0,JBP1
               DO I=0,IBP1
                  RHOS(I,J,K) = PBAR_S(K,PRESSURE_ZONE(I,J,K))/(TMPA*RSUM(I,J,K))
               ENDDO
            ENDDO
         ENDDO
      ENDIF
   ENDIF

   IF (MIXTURE_FRACTION) THEN
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (CO_PRODUCTION) THEN
                  Z_2 = YYS(I,J,K,I_PROG_CO)
               ElSE
                  Z_2 = 0._EB
               ENDIF
               CALL GET_MOLECULAR_WEIGHT(YYS(I,J,K,I_FUEL),Z_2,YYS(I,J,K,I_PROG_F),Y_SUM(I,J,K),RSUM(I,J,K))
               RSUM(I,J,K) = R0/RSUM(I,J,K)
            ENDDO
         ENDDO
      ENDDO
      IF (N_SPEC_DILUENTS > 0) RSUM = RSUM*(1._EB-Y_SUM) + R_SUM_DILUENTS
   ENDIF

   ! Extract predicted temperature at next time step from Equation of State

   IF (.NOT.ISOTHERMAL) THEN
      IF (N_SPECIES==0) THEN
         DO K=0,KBP1
            DO J=0,JBP1
               DO I=0,IBP1
                  TMP(I,J,K) = PBAR_S(K,PRESSURE_ZONE(I,J,K))/(SPECIES(0)%RCON*RHOS(I,J,K))
                  TMP(I,J,K) = MAX(TMPMIN,MIN(TMPMAX,TMP(I,J,K)))
               ENDDO
            ENDDO
         ENDDO
      ELSE
         DO K=0,KBP1
            DO J=0,JBP1
               DO I=0,IBP1
                  TMP(I,J,K) = PBAR_S(K,PRESSURE_ZONE(I,J,K))/(RSUM(I,J,K)*RHOS(I,J,K))
                  TMP(I,J,K) = MAX(TMPMIN,MIN(TMPMAX,TMP(I,J,K)))
               ENDDO
            ENDDO
         ENDDO
      ENDIF
      !TMP = MAX(TMPMIN,MIN(TMPMAX,TMP))
   ENDIF

! The CORRECTOR step   
ELSEIF (CORRECTOR) THEN

   ! Update the density
   
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            RHO(I,J,K) = RHOS(I,J,K) - DT*FRHO(I,J,K)
            RHO(I,J,K) = MAX(0._EB,RHO(I,J,K))
         ENDDO
      ENDDO
   ENDDO
   
   ! Update mass fractions
   
   DO N=1,N_SPECIES
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               YY(I,J,K,N) = RHOYYP(I,J,K,N) - DT*FRHOYY(I,J,K,N)
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   
   ! Corrector step
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            RHO(I,J,K) = 0.5_EB*( RHON(I,J,K) + RHO(I,J,K) )
            
            ! IF (CLIP_SCALAR) THEN...
            !RHO(I,J,K) = MAX(RHO(I,J,K),0._EB)
         ENDDO
      ENDDO
   ENDDO
   
   IF (BAROCLINIC) THEN
      ! Set global min and max for rho
      RHO_LOWER = MIN(RHO_LOWER,MINVAL_GASPHASE(RHO))
      RHO_UPPER = MAX(RHO_UPPER,MAXVAL_GASPHASE(RHO))
   ENDIF
   
   DO N=1,N_SPECIES
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               YY(I,J,K,N) = 0.5_EB*( YYN(I,J,K,N) + YY(I,J,K,N) )/RHO(I,J,K)
               
               ! IF (CLIP_SCALAR) THEN...
               YY(I,J,K,N) = MAX(YY(I,J,K,N),0._EB)
               YY(I,J,K,N) = MIN(YY(I,J,K,N),1._EB)
            ENDDO
         ENDDO
      ENDDO
   ENDDO

   ! Correct background pressure

   DO I=1,N_ZONE
      PBAR(:,I) = .5_EB*(PBAR(:,I) + PBAR_S(:,I) + D_PBAR_S_DT(I)*DT)
   ENDDO
 
   ! Compute mixture fraction and diluent sums: Y_SUM=Sum(Y_i), Z_SUM=Sum(Z_i)

   IF (MIXTURE_FRACTION) THEN
      Z_SUM  =  0._EB
      Y_SUM  =  0._EB
      IF (N_SPEC_DILUENTS > 0) THEN
         R_SUM_DILUENTS => WORK4
         R_SUM_DILUENTS =  0._EB
         ENDIF
      DO N=1,N_SPECIES
         IF (SPECIES(N)%MODE==MIXTURE_FRACTION_SPECIES) Z_SUM = Z_SUM + YY(:,:,:,N)
         IF (SPECIES(N)%MODE==GAS_SPECIES) THEN
            Y_SUM = Y_SUM + YY(:,:,:,N)
            R_SUM_DILUENTS(:,:,:) = R_SUM_DILUENTS(:,:,:) + SPECIES(N)%RCON*YY(:,:,:,N)
         ENDIF
      ENDDO
   ENDIF

   ! Compute molecular weight term RSUM=R0*SUM(Y_i/M_i)
 
   IF (N_SPECIES>0 .AND. .NOT. MIXTURE_FRACTION) THEN
      RSUM = SPECIES(0)%RCON
      DO N=1,N_SPECIES
         WFAC = SPECIES(N)%RCON - SPECIES(0)%RCON
         RSUM(:,:,:) = RSUM(:,:,:) + WFAC*YY(:,:,:,N)
      ENDDO
      IF (ISOTHERMAL) THEN
         DO K=0,KBP1
            DO J=0,JBP1
               DO I=0,IBP1
                  RHO(I,J,K) = PBAR(K,PRESSURE_ZONE(I,J,K))/(TMPA*RSUM(I,J,K))
               ENDDO
            ENDDO
         ENDDO
      ENDIF
   ENDIF

   IF (MIXTURE_FRACTION) THEN
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (CO_PRODUCTION) THEN
                  Z_2 = YY(I,J,K,I_PROG_CO)
               ElSE
                  Z_2 = 0._EB
               ENDIF
               CALL GET_MOLECULAR_WEIGHT(YY(I,J,K,I_FUEL),Z_2,YY(I,J,K,I_PROG_F),Y_SUM(I,J,K),RSUM(I,J,K))
               RSUM(I,J,K) = R0/RSUM(I,J,K)
            ENDDO
         ENDDO
      ENDDO
      IF (N_SPEC_DILUENTS > 0) RSUM = RSUM*(1._EB-Y_SUM) + R_SUM_DILUENTS
   ENDIF

   ! Extract temperature from the Equation of State

   IF (.NOT.ISOTHERMAL) THEN
      IF (N_SPECIES==0) THEN
         DO K=0,KBP1
            DO J=0,JBP1
               DO I=0,IBP1
                  TMP(I,J,K) = PBAR(K,PRESSURE_ZONE(I,J,K))/(SPECIES(0)%RCON*RHO(I,J,K))
                  TMP(I,J,K) = MAX(TMPMIN,MIN(TMPMAX,TMP(I,J,K)))
               ENDDO
            ENDDO
         ENDDO
      ELSE
         DO K=0,KBP1
            DO J=0,JBP1
               DO I=0,IBP1
                  TMP(I,J,K) = PBAR(K,PRESSURE_ZONE(I,J,K))/(RSUM(I,J,K)*RHO(I,J,K))
                  TMP(I,J,K) = MAX(TMPMIN,MIN(TMPMAX,TMP(I,J,K)))
               ENDDO
            ENDDO
         ENDDO
      ENDIF
      !TMP = MAX(TMPMIN,MIN(TMPMAX,TMP))
   ENDIF

ENDIF SELECT_SUBSTEP

TUSED(3,NM)=TUSED(3,NM)+SECOND()-TNOW
 
END SUBROUTINE DENSITY_TVD


SUBROUTINE SCALARF

USE GLOBAL_CONSTANTS, ONLY: N_SPECIES,PREDICTOR,CORRECTOR,FLUX_LIMITER,NULL_BOUNDARY,POROUS_BOUNDARY,OPEN_BOUNDARY, &
                            LU_ERR,BTOL,INTERPOLATED_BOUNDARY

! Computes the divergence of the scalar advective flux + diffusion

INTEGER :: I,J,K,N,II,JJ,KK,IOR,IW,IIG,JJG,KKG,ICM,ICP
REAL(EB) :: ZZ(4)
REAL(EB), POINTER, DIMENSION(:,:,:) :: RHOP,UU,VV,WW,FX,FY,FZ
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: YYP, RHOYYP

FX => WORK1
FY => WORK2
FZ => WORK3
FX = 0._EB
FY = 0._EB
FZ = 0._EB
IF (N_SPECIES>0) RHOYYP => SCALAR_SAVE3

IF (PREDICTOR) THEN
   UU => U
   VV => V
   WW => W
   RHOP => RHO
   IF (N_SPECIES > 0) YYP => YY
ELSE
   UU => US
   VV => VS
   WW => WS
   RHOP => RHOS
   IF (N_SPECIES > 0) YYP => YYS
ENDIF

DO N=1,N_SPECIES
   DO K=0,KBP1
      DO J=0,JBP1
         DO I=0,IBP1
            RHOYYP(I,J,K,N) = RHOP(I,J,K)*YYP(I,J,K,N)
         ENDDO
      ENDDO
   ENDDO
ENDDO

! Density flux

DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBM1
         ZZ(1:4)   = RHOP(I-1:I+2,J,K)
         ICM = CELL_INDEX(I,J,K)
         ICP = CELL_INDEX(I+1,J,K)
         IF (WALL_INDEX(ICM,-1)==0 .AND. WALL_INDEX(ICP,1)==0) THEN
            FX(I,J,K) = UU(I,J,K)*SCALAR_FACE_VALUE(UU(I,J,K),ZZ,FLUX_LIMITER)
         ELSE
            FX(I,J,K) = UU(I,J,K)*SCALAR_FACE_VALUE(UU(I,J,K),ZZ,1)
         ENDIF
      ENDDO
   ENDDO
ENDDO

DO K=1,KBAR
   DO J=1,JBM1
      DO I=1,IBAR
         ZZ(1:4)   = RHOP(I,J-1:J+2,K)
         ICM = CELL_INDEX(I,J,K)
         ICP = CELL_INDEX(I,J+1,K)
         IF (WALL_INDEX(ICM,-2)==0 .AND. WALL_INDEX(ICP,2)==0) THEN
            FY(I,J,K) = VV(I,J,K)*SCALAR_FACE_VALUE(VV(I,J,K),ZZ,FLUX_LIMITER)
         ELSE
            FY(I,J,K) = VV(I,J,K)*SCALAR_FACE_VALUE(VV(I,J,K),ZZ,1)
         ENDIF
      ENDDO
   ENDDO
ENDDO

DO K=1,KBM1
   DO J=1,JBAR
      DO I=1,IBAR
         ZZ(1:4)   = RHOP(I,J,K-1:K+2)
         ICM = CELL_INDEX(I,J,K)
         ICP = CELL_INDEX(I,J,K+1)
         IF (WALL_INDEX(ICM,-3)==0 .AND. WALL_INDEX(ICP,3)==0) THEN
            FZ(I,J,K) = WW(I,J,K)*SCALAR_FACE_VALUE(WW(I,J,K),ZZ,FLUX_LIMITER)
         ELSE
            FZ(I,J,K) = WW(I,J,K)*SCALAR_FACE_VALUE(WW(I,J,K),ZZ,1)
         ENDIF
      ENDDO
   ENDDO
ENDDO

! Compute mass fluxes at boundaries

WALL_LOOP: DO IW=1,NWC
   IF (BOUNDARY_TYPE(IW)==NULL_BOUNDARY .OR. BOUNDARY_TYPE(IW)==POROUS_BOUNDARY) CYCLE WALL_LOOP
   II  = IJKW(1,IW)
   JJ  = IJKW(2,IW)
   KK  = IJKW(3,IW)
   IIG = IJKW(6,IW)
   JJG = IJKW(7,IW)
   KKG = IJKW(8,IW)
   IOR = IJKW(4,IW)
   SELECT CASE(IOR)
      CASE( 1)
         ZZ(2) = RHO_W(IW)
         ZZ(3) = RHOP(IIG,JJG,KKG)
         IF (BOUNDARY_TYPE(IW)/=INTERPOLATED_BOUNDARY) THEN
            FX(II,JJ,KK) = UU(II,JJ,KK)*SCALAR_FACE_VALUE(UU(II,JJ,KK),ZZ,1)
         ELSE
            FX(II,JJ,KK) = UVW_SAVE(IW)*SCALAR_FACE_VALUE(UVW_SAVE(IW),ZZ,1)
         ENDIF
      CASE(-1)
         ZZ(2) = RHOP(IIG,JJG,KKG)
         ZZ(3) = RHO_W(IW)
         IF (BOUNDARY_TYPE(IW)/=INTERPOLATED_BOUNDARY) THEN
            FX(II-1,JJ,KK) = UU(II-1,JJ,KK)*SCALAR_FACE_VALUE(UU(II-1,JJ,KK),ZZ,1)
         ELSE
            FX(II-1,JJ,KK) = UVW_SAVE(IW)*SCALAR_FACE_VALUE(UVW_SAVE(IW),ZZ,1)
         ENDIF
      CASE( 2)
         ZZ(2) = RHO_W(IW)
         ZZ(3) = RHOP(IIG,JJG,KKG)
         IF (BOUNDARY_TYPE(IW)/=INTERPOLATED_BOUNDARY) THEN
            FY(II,JJ,KK) = VV(II,JJ,KK)*SCALAR_FACE_VALUE(VV(II,JJ,KK),ZZ,1)
         ELSE
            FY(II,JJ,KK) = UVW_SAVE(IW)*SCALAR_FACE_VALUE(UVW_SAVE(IW),ZZ,1)
         ENDIF
      CASE(-2)
         ZZ(2) = RHOP(IIG,JJG,KKG)
         ZZ(3) = RHO_W(IW)
         IF (BOUNDARY_TYPE(IW)/=INTERPOLATED_BOUNDARY) THEN
            FY(II,JJ-1,KK) = VV(II,JJ-1,KK)*SCALAR_FACE_VALUE(VV(II,JJ-1,KK),ZZ,1)
         ELSE
            FY(II,JJ-1,KK) = UVW_SAVE(IW)*SCALAR_FACE_VALUE(UVW_SAVE(IW),ZZ,1)
         ENDIF
      CASE( 3)
         ZZ(2) = RHO_W(IW)
         ZZ(3) = RHOP(IIG,JJG,KKG)
         IF (BOUNDARY_TYPE(IW)/=INTERPOLATED_BOUNDARY) THEN
            FZ(II,JJ,KK) = WW(II,JJ,KK)*SCALAR_FACE_VALUE(WW(II,JJ,KK),ZZ,1)
         ELSE
            FZ(II,JJ,KK) = UVW_SAVE(IW)*SCALAR_FACE_VALUE(UVW_SAVE(IW),ZZ,1)
         ENDIF
      CASE(-3)
         ZZ(2) = RHOP(IIG,JJG,KKG)
         ZZ(3) = RHO_W(IW)
         IF (BOUNDARY_TYPE(IW)/=INTERPOLATED_BOUNDARY) THEN
            FZ(II,JJ,KK-1) = WW(II,JJ,KK-1)*SCALAR_FACE_VALUE(WW(II,JJ,KK-1),ZZ,1)
         ELSE
            FZ(II,JJ,KK-1) = UVW_SAVE(IW)*SCALAR_FACE_VALUE(UVW_SAVE(IW),ZZ,1)
         ENDIF
   END SELECT
ENDDO WALL_LOOP
         
! Compute divergence of advective flux for density

DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         FRHO(I,J,K) = RDX(I)*(FX(I,J,K)-FX(I-1,J,K)) &
                     + RDY(J)*(FY(I,J,K)-FY(I,J-1,K)) &
                     + RDZ(K)*(FZ(I,J,K)-FZ(I,J,K-1))
      ENDDO
   ENDDO
ENDDO

! Compute species fluxes at interior cell faces

SPECIES_LOOP: DO N=1,N_SPECIES

   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBM1
            ZZ(1:4)   = RHOYYP(I-1:I+2,J,K,N)
            ICM = CELL_INDEX(I,J,K)
            ICP = CELL_INDEX(I+1,J,K)
            IF (WALL_INDEX(ICM,-1)==0 .AND. WALL_INDEX(ICP,1)==0) THEN
               FX(I,J,K) = UU(I,J,K)*SCALAR_FACE_VALUE(UU(I,J,K),ZZ,FLUX_LIMITER)
            ELSE
               FX(I,J,K) = UU(I,J,K)*SCALAR_FACE_VALUE(UU(I,J,K),ZZ,1)
            ENDIF
         ENDDO
      ENDDO
   ENDDO

   DO K=1,KBAR
      DO J=1,JBM1
         DO I=1,IBAR
            ZZ(1:4)   = RHOYYP(I,J-1:J+2,K,N)
            ICM = CELL_INDEX(I,J,K)
            ICP = CELL_INDEX(I,J+1,K)
            IF (WALL_INDEX(ICM,-2)==0 .AND. WALL_INDEX(ICP,2)==0) THEN
               FY(I,J,K) = VV(I,J,K)*SCALAR_FACE_VALUE(VV(I,J,K),ZZ,FLUX_LIMITER)
            ELSE
               FY(I,J,K) = VV(I,J,K)*SCALAR_FACE_VALUE(VV(I,J,K),ZZ,1)
            ENDIF
         ENDDO
      ENDDO
   ENDDO

   DO K=1,KBM1
      DO J=1,JBAR
         DO I=1,IBAR
            ZZ(1:4)   = RHOYYP(I,J,K-1:K+2,N)
            ICM = CELL_INDEX(I,J,K)
            ICP = CELL_INDEX(I,J,K+1)
            IF (WALL_INDEX(ICM,-3)==0 .AND. WALL_INDEX(ICP,3)==0) THEN
               FZ(I,J,K) = WW(I,J,K)*SCALAR_FACE_VALUE(WW(I,J,K),ZZ,FLUX_LIMITER)
            ELSE
               FZ(I,J,K) = WW(I,J,K)*SCALAR_FACE_VALUE(WW(I,J,K),ZZ,1)
            ENDIF
         ENDDO
      ENDDO
   ENDDO

   ! Compute species fluxes at boundaries

   WALL_LOOP2: DO IW=1,NWC
      IF (BOUNDARY_TYPE(IW)==NULL_BOUNDARY .OR. BOUNDARY_TYPE(IW)==POROUS_BOUNDARY) CYCLE WALL_LOOP2
      II  = IJKW(1,IW)
      JJ  = IJKW(2,IW)
      KK  = IJKW(3,IW)
      IIG = IJKW(6,IW)
      JJG = IJKW(7,IW)
      KKG = IJKW(8,IW)
      IOR = IJKW(4,IW)
      SELECT CASE(IOR)
      CASE( 1)
         ZZ(2) = RHO_W(IW)*YY_W(IW,N)
         ZZ(3) = RHOYYP(IIG,JJG,KKG,N)
         IF (BOUNDARY_TYPE(IW)/=INTERPOLATED_BOUNDARY) THEN
            FX(II,JJ,KK) = UU(II,JJ,KK)*SCALAR_FACE_VALUE(UU(II,JJ,KK),ZZ,1)
         ELSE
            FX(II,JJ,KK) = UVW_SAVE(IW)*SCALAR_FACE_VALUE(UVW_SAVE(IW),ZZ,1)
         ENDIF
      CASE(-1)
         ZZ(2) = RHOYYP(IIG,JJG,KKG,N)
         ZZ(3) = RHO_W(IW)*YY_W(IW,N)
         IF (BOUNDARY_TYPE(IW)/=INTERPOLATED_BOUNDARY) THEN
            FX(II-1,JJ,KK) = UU(II-1,JJ,KK)*SCALAR_FACE_VALUE(UU(II-1,JJ,KK),ZZ,1)
         ELSE
            FX(II-1,JJ,KK) = UVW_SAVE(IW)*SCALAR_FACE_VALUE(UVW_SAVE(IW),ZZ,1)
         ENDIF
      CASE( 2)
         ZZ(2) = RHO_W(IW)*YY_W(IW,N)
         ZZ(3) = RHOYYP(IIG,JJG,KKG,N)
         IF (BOUNDARY_TYPE(IW)/=INTERPOLATED_BOUNDARY) THEN
            FY(II,JJ,KK) = VV(II,JJ,KK)*SCALAR_FACE_VALUE(VV(II,JJ,KK),ZZ,1)
         ELSE
            FY(II,JJ,KK) = UVW_SAVE(IW)*SCALAR_FACE_VALUE(UVW_SAVE(IW),ZZ,1)
         ENDIF
      CASE(-2)
         ZZ(2) = RHOYYP(IIG,JJG,KKG,N)
         ZZ(3) = RHO_W(IW)*YY_W(IW,N)
         IF (BOUNDARY_TYPE(IW)/=INTERPOLATED_BOUNDARY) THEN
            FY(II,JJ-1,KK) = VV(II,JJ-1,KK)*SCALAR_FACE_VALUE(VV(II,JJ-1,KK),ZZ,1)
         ELSE
            FY(II,JJ-1,KK) = UVW_SAVE(IW)*SCALAR_FACE_VALUE(UVW_SAVE(IW),ZZ,1)
         ENDIF
      CASE( 3)
         ZZ(2) = RHO_W(IW)*YY_W(IW,N)
         ZZ(3) = RHOYYP(IIG,JJG,KKG,N)
         IF (BOUNDARY_TYPE(IW)/=INTERPOLATED_BOUNDARY) THEN
            FZ(II,JJ,KK) = WW(II,JJ,KK)*SCALAR_FACE_VALUE(WW(II,JJ,KK),ZZ,1)
         ELSE
            FZ(II,JJ,KK) = UVW_SAVE(IW)*SCALAR_FACE_VALUE(UVW_SAVE(IW),ZZ,1)
         ENDIF
      CASE(-3)
         ZZ(2) = RHOYYP(IIG,JJG,KKG,N)
         ZZ(3) = RHO_W(IW)*YY_W(IW,N)
         IF (BOUNDARY_TYPE(IW)/=INTERPOLATED_BOUNDARY) THEN
            FZ(II,JJ,KK-1) = WW(II,JJ,KK-1)*SCALAR_FACE_VALUE(WW(II,JJ,KK-1),ZZ,1)
         ELSE
            FZ(II,JJ,KK-1) = UVW_SAVE(IW)*SCALAR_FACE_VALUE(UVW_SAVE(IW),ZZ,1)
         ENDIF
      END SELECT
   ENDDO WALL_LOOP2

   ! Compute divergence of advective flux for species, then add diffusion term

   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            FRHOYY(I,J,K,N) = RDX(I)*(FX(I,J,K)-FX(I-1,J,K)) &
                            + RDY(J)*(FY(I,J,K)-FY(I,J-1,K)) &
                            + RDZ(K)*(FZ(I,J,K)-FZ(I,J,K-1)) &
                            - DEL_RHO_D_DEL_Y(I,J,K,N)
         ENDDO
      ENDDO
   ENDDO

ENDDO SPECIES_LOOP

END SUBROUTINE SCALARF


REAL(EB) FUNCTION SCALAR_FACE_VALUE(A,U,LIMITER)

REAL(EB), INTENT(IN) :: A,U(4)
INTEGER, INTENT(IN) :: LIMITER

! local
REAL(EB) :: R,B,DU_UP,DU_LOC

! This function computes the scalar value on a face.
! The scalar is denoted U, and the velocity is denoted A.
! The divergence (computed elsewhere) uses a central difference across 
! the cell subject to a flux LIMITER.  The flux LIMITER choices are:
! 
! LIMITER = 0 implements central differencing
! LIMITER = 1 implements first-order upwinding (monotone)
! LIMITER = 2 implements the SUPERBEE (SB) LIMITER of Roe
! LIMITER = 3 implements the MINMOD LIMITER
! LIMITER = 4 implements the CHARM LIMITER
!
!                    location of face
!                            
!                            f
!    |     o     |     o     |     o     |     o     |
!                            A
!         U(1)        U(2)        U(3)        U(4)

IF (A>0._EB) THEN
    
   ! the flow is left to right
   DU_UP  = U(2)-U(1)
   DU_LOC = U(3)-U(2)

   R = 0._EB
   B = 0._EB

   SELECT CASE(LIMITER)
      CASE(0) ! central differencing
         SCALAR_FACE_VALUE = 0.5_EB*(U(2)+U(3))
      CASE(1) ! first-order upwinding
         SCALAR_FACE_VALUE = U(2)
      CASE(2) ! SUPERBEE, Roe (1986)
         IF (ABS(DU_LOC)>0._EB) R = DU_UP/DU_LOC
         B = MAX(0._EB,MIN(2._EB*R,1._EB),MIN(R,2._EB))
         SCALAR_FACE_VALUE = U(2) + 0.5_EB*B*(U(3)-U(2))
      CASE(3) ! MINMOD
         IF (ABS(DU_LOC)>0._EB) R = DU_UP/DU_LOC
         B = MAX(0._EB,MIN(1._EB,R))
         SCALAR_FACE_VALUE = U(2) + 0.5_EB*B*(U(3)-U(2))
      CASE(4) ! CHARM
         IF (ABS(DU_UP)>0._EB) R = DU_LOC/DU_UP
         IF (R>0._EB) B = R*(3._EB*R+1._EB)/((R+1._EB)**2)
         SCALAR_FACE_VALUE = U(2) + 0.5_EB*B*(U(2)-U(1))
   END SELECT
    
ELSE

   ! the flow is right to left
   DU_UP  = U(4)-U(3)
   DU_LOC = U(3)-U(2)

   R = 0._EB
   B = 0._EB

   SELECT CASE(LIMITER)
      CASE(0) ! central differencing
         SCALAR_FACE_VALUE = 0.5_EB*(U(2)+U(3))
      CASE(1) ! first-order upwinding
         SCALAR_FACE_VALUE = U(3)
      CASE(2) ! SUPERBEE, Roe (1986)
         IF (ABS(DU_LOC)>0._EB) R = DU_UP/DU_LOC
         B = MAX(0._EB,MIN(2._EB*R,1._EB),MIN(R,2._EB))
         SCALAR_FACE_VALUE = U(3) + 0.5_EB*B*(U(2)-U(3))
      CASE(3) ! MINMOD
         IF (ABS(DU_LOC)>0._EB) R = DU_UP/DU_LOC
         B = MAX(0._EB,MIN(1._EB,R))
         SCALAR_FACE_VALUE = U(3) + 0.5_EB*B*(U(2)-U(3))
      CASE(4) ! CHARM
         IF (ABS(DU_UP)>0._EB) R = DU_LOC/DU_UP
         IF (R>0._EB) B = R*(3._EB*R+1._EB)/((R+1._EB)**2)
         SCALAR_FACE_VALUE = U(3) + 0.5_EB*B*(U(3)-U(4))
    END SELECT
    
ENDIF

END FUNCTION SCALAR_FACE_VALUE



REAL(EB) FUNCTION MINVAL_GASPHASE(PHI)

REAL(EB), INTENT(IN) :: PHI(0:IBP1,0:JBP1,0:KBP1)

! local
INTEGER :: I,J,K

MINVAL_GASPHASE = HUGE(1._EB)
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IF (.NOT.SOLID(CELL_INDEX(I,J,K))) THEN
            IF (PHI(I,J,K)<MINVAL_GASPHASE) THEN
               MINVAL_GASPHASE = PHI(I,J,K)
            ENDIF
         ENDIF
      ENDDO
   ENDDO
ENDDO

END FUNCTION MINVAL_GASPHASE


REAL(EB) FUNCTION MAXVAL_GASPHASE(PHI)

REAL(EB), INTENT(IN) :: PHI(0:IBP1,0:JBP1,0:KBP1)

! local
INTEGER :: I,J,K

MAXVAL_GASPHASE = -HUGE(1._EB)
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IF (.NOT.SOLID(CELL_INDEX(I,J,K))) THEN
            IF (PHI(I,J,K)>MAXVAL_GASPHASE) THEN
               MAXVAL_GASPHASE = PHI(I,J,K)
            ENDIF
         ENDIF
      ENDDO
   ENDDO
ENDDO

END FUNCTION MAXVAL_GASPHASE


!---------------------------------------------------------------------------

SUBROUTINE GET_REV_mass(MODULE_REV,MODULE_DATE)
INTEGER,INTENT(INOUT) :: MODULE_REV
CHARACTER(255),INTENT(INOUT) :: MODULE_DATE
WRITE(MODULE_DATE,'(A)') massrev(INDEX(massrev,':')+1:LEN_TRIM(massrev)-2)
READ (MODULE_DATE,'(I5)') MODULE_REV
WRITE(MODULE_DATE,'(A)') massdate
END SUBROUTINE GET_REV_mass
 
END MODULE MASS
