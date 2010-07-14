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

PUBLIC MASS_FINITE_DIFFERENCES,DENSITY,GET_REV_mass,SCALARF,DENSITY_TVD
 
 
CONTAINS
 
SUBROUTINE MASS_FINITE_DIFFERENCES(NM)
USE COMP_FUNCTIONS, ONLY: SECOND
USE GLOBAL_CONSTANTS, ONLY: N_SPECIES,ISOTHERMAL,NULL_BOUNDARY,POROUS_BOUNDARY,PREDICTOR,CORRECTOR,EVACUATION_ONLY, &
                            SOLID_PHASE_ONLY,TUSED,DEBUG_OPENMP,NOBIAS
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

IF (NOBIAS) PMDT = 0._EB

! Define local CFL numbers
 
EPSX => WORK1
EPSY => WORK2
EPSZ => WORK3

!$OMP PARALLEL 
!$OMP DO COLLAPSE(3) PRIVATE(K,J,I)
DO K=0,KBAR
   DO J=0,JBAR
      DO I=0,IBAR
         !!$ IF ((K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_MASS_FD_01'
         EPSX(I,J,K) = PMDT*UU(I,J,K)*RDXN(I)
         EPSY(I,J,K) = PMDT*VV(I,J,K)*RDYN(J)
         EPSZ(I,J,K) = PMDT*WW(I,J,K)*RDZN(K)
      ENDDO
   ENDDO
ENDDO
!$OMP END DO

! Compute spatial differences for density equation
 
NOT_ISOTHERMAL_IF: IF (.NOT.ISOTHERMAL) THEN
   
   !$OMP SINGLE
   UDRHODX => WORK4
   VDRHODY => WORK5
   WDRHODZ => WORK6
   !$OMP END SINGLE
   
   !$OMP DO COLLAPSE(3) PRIVATE(K,J,I)
   DO K=0,KBAR
      DO J=0,JBAR
         DO I=0,IBAR
            !!$ IF ((K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_MASS_FD_02'
            UDRHODX(I,J,K) = UU(I,J,K)*(RHOP(I+1,J,K)-RHOP(I,J,K))*RDXN(I)
            VDRHODY(I,J,K) = VV(I,J,K)*(RHOP(I,J+1,K)-RHOP(I,J,K))*RDYN(J)
            WDRHODZ(I,J,K) = WW(I,J,K)*(RHOP(I,J,K+1)-RHOP(I,J,K))*RDZN(K)
         ENDDO
      ENDDO
   ENDDO
   !$OMP END DO

   !$OMP DO PRIVATE(IW,II,JJ,KK,IIG,JJG,KKG,IOR,UDRHODN)
   WLOOP: DO IW=1,NWC
      !!$ IF ((IW == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_MASS_FD_03'
      IF (BOUNDARY_TYPE(IW)==NULL_BOUNDARY .OR. BOUNDARY_TYPE(IW)==POROUS_BOUNDARY) CYCLE WLOOP
      II  = IJKW(1,IW) 
      IIG = IJKW(6,IW)
      JJ  = IJKW(2,IW) 
      JJG = IJKW(7,IW)
      KK  = IJKW(3,IW) 
      KKG = IJKW(8,IW)
      IOR = IJKW(4,IW)
      UDRHODN = 2._EB*UWP(IW)*(RHO_F(IW)-RHOP(IIG,JJG,KKG))*RDN(IW)
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
   !$OMP END DO

   !$OMP WORKSHARE
   FRHO = 0._EB
   !$OMP END WORKSHARE
   
   !$OMP DO COLLAPSE(3) PRIVATE(K,J,I,FXYZ)
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
   !$OMP END DO
ENDIF NOT_ISOTHERMAL_IF
!$OMP END PARALLEL 


! Compute the species equation differences
 
IF (N_SPECIES > 0) THEN
   IF (PREDICTOR) YYP => YY
   IF (CORRECTOR) YYP => YYS
   UDRHODX => WORK4
   VDRHODY => WORK5
   WDRHODZ => WORK6
ENDIF
 
SPECIES_LOOP: DO N=1,N_SPECIES

   !$OMP PARALLEL
   !$OMP DO COLLAPSE(3) PRIVATE(K,J,I)
   DO K=0,KBAR
      DO J=0,JBAR
         DO I=0,IBAR
            !!$ IF ((K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_MASS_FD_04'
            UDRHODX(I,J,K) = UU(I,J,K)*( RHOP(I+1,J,K)*YYP(I+1,J,K,N)-RHOP(I,J,K)*YYP(I,J,K,N) )*RDXN(I)
            VDRHODY(I,J,K) = VV(I,J,K)*( RHOP(I,J+1,K)*YYP(I,J+1,K,N)-RHOP(I,J,K)*YYP(I,J,K,N) )*RDYN(J)
            WDRHODZ(I,J,K) = WW(I,J,K)*( RHOP(I,J,K+1)*YYP(I,J,K+1,N)-RHOP(I,J,K)*YYP(I,J,K,N) )*RDZN(K)
         ENDDO
      ENDDO
   ENDDO
   !$OMP END DO
 
   ! Correct U d(RHO*Y)/dx etc. on boundaries

   !$OMP DO PRIVATE(IW,II,JJ,KK,IIG,JJG,KKG,IOR,UDRHODN) 
   WLOOP2: DO IW=1,NWC
      !!$ IF ((IW == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_MASS_FD_05'
      IF (BOUNDARY_TYPE(IW)==NULL_BOUNDARY .OR. BOUNDARY_TYPE(IW)==POROUS_BOUNDARY) CYCLE WLOOP2
      II  = IJKW(1,IW) 
      IIG = IJKW(6,IW)
      JJ  = IJKW(2,IW) 
      JJG = IJKW(7,IW)
      KK  = IJKW(3,IW) 
      KKG = IJKW(8,IW)
      IOR = IJKW(4,IW)
      UDRHODN = 2._EB*UWP(IW)*( RHO_F(IW)*YY_F(IW,N) - RHOP(IIG,JJG,KKG)*YYP(IIG,JJG,KKG,N) )*RDN(IW)
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
   !$OMP END DO
 
  ! Sum up the convective and diffusive terms in the transport equation and store in DEL_RHO_D_DEL_Y

   !$OMP DO COLLAPSE(3) PRIVATE(K,J,I,FXYZ)
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
   !$OMP END DO
   !$OMP END PARALLEL
 
ENDDO SPECIES_LOOP
 
TUSED(3,NM)=TUSED(3,NM)+SECOND()-TNOW
END SUBROUTINE MASS_FINITE_DIFFERENCES

 
SUBROUTINE DENSITY(NM)

! Update the density and species mass fractions

USE COMP_FUNCTIONS, ONLY: SECOND 
USE PHYSICAL_FUNCTIONS, ONLY : GET_SPECIFIC_GAS_CONSTANT
USE GLOBAL_CONSTANTS, ONLY: N_SPECIES,CO_PRODUCTION,I_PROG_F,I_PROG_CO,I_FUEL,TMPMAX,TMPMIN,EVACUATION_ONLY,PREDICTOR,CORRECTOR, &
                            CHANGE_TIME_STEP,ISOTHERMAL,TMPA,N_ZONE,MIXTURE_FRACTION_SPECIES, &
                            GAS_SPECIES, MIXTURE_FRACTION,R0,SOLID_PHASE_ONLY,TUSED, &
                            RSUM0,DEBUG_OPENMP
REAL(EB) :: WFAC,DTRATIO,OMDTRATIO,TNOW,YY_GET(1:N_SPECIES)
INTEGER  :: I,J,K,N
INTEGER, INTENT(IN) :: NM
 
IF (EVACUATION_ONLY(NM)) RETURN
IF (SOLID_PHASE_ONLY) RETURN

TNOW=SECOND()
CALL POINT_TO_MESH(NM)

PREDICTOR_STEP: SELECT CASE (PREDICTOR)

CASE(.TRUE.) PREDICTOR_STEP

   !$OMP PARALLEL PRIVATE(DTRATIO,OMDTRATIO)

   IF (.NOT.CHANGE_TIME_STEP(NM)) THEN

      !$OMP DO COLLAPSE(4) PRIVATE(N,K,J,I)
      DO N=1,N_SPECIES
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  !!$ IF ((N == 1) .AND. (K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_MASS_DENS_01'
                  IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
                  YYS(I,J,K,N) = RHO(I,J,K)*YY(I,J,K,N) - DT*DEL_RHO_D_DEL_Y(I,J,K,N)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO

   ELSE

      DTRATIO   = DT/DT_PREV
      OMDTRATIO = 1._EB - DTRATIO
      !$OMP DO COLLAPSE(4) PRIVATE(N,K,J,I)
      DO N=1,N_SPECIES
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  !!$ IF ((N == 1) .AND. (K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_MASS_DENS_02'
                  IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
                  YYS(I,J,K,N) = OMDTRATIO*RHO(I,J,K) *YY(I,J,K,N) + DTRATIO*RHOS(I,J,K)*YYS(I,J,K,N)
               ENDDO
           ENDDO
         ENDDO
      ENDDO
      !$OMP END DO

   ENDIF
   !$OMP END PARALLEL

   ! Predict the density at the next time step (RHOS or RHO^*)

   IF (.NOT.ISOTHERMAL) THEN
      !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(K,J,I)
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               !!$ IF ((K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_MASS_DENS_03'
               RHOS(I,J,K) = RHO(I,J,K)-DT*FRHO(I,J,K)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END PARALLEL DO
   ELSE

      !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(K,J,I)
      DO K=0,KBP1
         DO J=0,JBP1
            DO I=0,IBP1
               !!$ IF ((K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_MASS_DENS_04'
               RHOS(I,J,K) = PBAR_S(K,PRESSURE_ZONE(I,J,K))/(TMPA*SPECIES(0)%RCON)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END PARALLEL DO

      DO N=1,N_SPECIES
         WFAC = 1._EB - SPECIES(N)%RCON/SPECIES(0)%RCON
         !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(K,J,I)
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  !!$ IF ((N == 1) .AND. (K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_MASS_DENS_05'
                  RHOS(I,J,K) = RHOS(I,J,K) + WFAC*YYS(I,J,K,N)
               ENDDO
            ENDDO
         ENDDO
         !$OMP END PARALLEL DO
      ENDDO

   ENDIF
 
   ! Correct densities above or below clip limits

   CALL CHECK_DENSITY
   
   ! Extract mass fraction from RHO * YY

   !$OMP PARALLEL DO COLLAPSE(4) PRIVATE(N,K,J,I)
   DO N=1,N_SPECIES
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
               YYS(I,J,K,N) = YYS(I,J,K,N)/RHOS(I,J,K)
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO

   ! Correct mass fractions above or below clip limits

   CALL CHECK_MASS_FRACTION

   ! Predict background pressure at next time step

   DO I=1,N_ZONE
      PBAR_S(:,I) = PBAR(:,I) + D_PBAR_DT(I)*DT
   ENDDO

   ! Compute molecular weight term RSUM=R0*SUM(Y_i/M_i)

   !$OMP PARALLEL SHARED(RSUM)
   IF (N_SPECIES>0) THEN
      !$OMP DO COLLAPSE(3) PRIVATE(K,J,I,YY_GET)
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR   
               !!$ IF ((K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_MASS_DENS_06'
               YY_GET(:) = YYS(I,J,K,:)
               CALL GET_SPECIFIC_GAS_CONSTANT(YY_GET,RSUM(I,J,K))
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO
      IF (ISOTHERMAL) THEN
         !$OMP DO COLLAPSE(3) PRIVATE(K,J,I)
         DO K=0,KBP1
            DO J=0,JBP1
               DO I=0,IBP1
                  !!$ IF ((K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_MASS_DENS_07'
                  RHOS(I,J,K) = PBAR_S(K,PRESSURE_ZONE(I,J,K))/(TMPA*RSUM(I,J,K))
               ENDDO
            ENDDO
         ENDDO
         !$OMP END DO
      ENDIF
   ENDIF

! Extract predicted temperature at next time step from Equation of State
   IF (.NOT. ISOTHERMAL) THEN
      IF (N_SPECIES==0) THEN
         !$OMP DO COLLAPSE(3) PRIVATE(K,J,I)
         DO K=0,KBP1
            DO J=0,JBP1
               DO I=0,IBP1
                  !!$ IF ((K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_MASS_DENS_08'
                  TMP(I,J,K) = PBAR_S(K,PRESSURE_ZONE(I,J,K))/(RSUM0*RHOS(I,J,K))
               ENDDO
            ENDDO
         ENDDO
         !$OMP END DO
      ELSE
         !$OMP DO COLLAPSE(3) PRIVATE(K,J,I)
         DO K=0,KBP1
            DO J=0,JBP1
               DO I=0,IBP1
                  !!$ IF ((K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_MASS_DENS_09'
                  TMP(I,J,K) = PBAR_S(K,PRESSURE_ZONE(I,J,K))/(RSUM(I,J,K)*RHOS(I,J,K))
               ENDDO
            ENDDO
         ENDDO
         !$OMP END DO
      ENDIF
      TMP = MAX(TMPMIN,MIN(TMPMAX,TMP))
   ENDIF
   !$OMP END PARALLEL

! The CORRECTOR step
   
CASE(.FALSE.) PREDICTOR_STEP

   ! Correct species mass fraction at next time step (YY here actually means YY*RHO)

   !$OMP PARALLEL DO COLLAPSE(4) PRIVATE(N,K,J,I)
   DO N=1,N_SPECIES
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               !!$ IF ((N == 1) .AND. (K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_MASS_DENS_10'
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
               YY(I,J,K,N) = .5_EB*(RHO(I,J,K)*YY(I,J,K,N) + RHOS(I,J,K)*YYS(I,J,K,N) - DT*DEL_RHO_D_DEL_Y(I,J,K,N) ) 
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO

   ! Correct density at next time step

   IF (.NOT.ISOTHERMAL) THEN
      !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(K,J,I)
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               !!$ IF ((K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_MASS_DENS_11'
               RHO(I,J,K) = .5_EB*(RHO(I,J,K)+RHOS(I,J,K)-DT*FRHO(I,J,K))
            ENDDO
         ENDDO
      ENDDO
      !$OMP END PARALLEL DO
   ELSE
      !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(K,J,I)
      DO K=0,KBP1
         DO J=0,JBP1
            DO I=0,IBP1
               !!$ IF ((K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_MASS_DENS_12'
               RHO(I,J,K) = PBAR(K,PRESSURE_ZONE(I,J,K))/(SPECIES(0)%RCON*TMPA)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END PARALLEL DO

      DO N=1,N_SPECIES
         WFAC = 1._EB - SPECIES(N)%RCON/SPECIES(0)%RCON
         !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(K,J,I)
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  !!$ IF ((N == 1) .AND. (K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_MASS_DENS_13'
                  RHO(I,J,K) = RHO(I,J,K) + WFAC*YY(I,J,K,N)
               ENDDO
            ENDDO
         ENDDO
         !$OMP END PARALLEL DO
      ENDDO
   ENDIF

   ! Correct densities above or below clip limits

   CALL CHECK_DENSITY
   
   ! Extract Y_n from rho*Y_n

   !$OMP PARALLEL DO COLLAPSE(4) PRIVATE(N,K,J,I)
   DO N=1,N_SPECIES
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               !!$ IF ((N == 1) .AND. (K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_MASS_DENS_14'
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
               YY(I,J,K,N) = YY(I,J,K,N)/RHO(I,J,K)
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO

   ! Correct mass fractions above or below clip limits

   CALL CHECK_MASS_FRACTION

   ! Correct background pressure

   DO I=1,N_ZONE
      PBAR(:,I) = .5_EB*(PBAR(:,I) + PBAR_S(:,I) + D_PBAR_S_DT(I)*DT)
   ENDDO
 
   ! Compute molecular weight term RSUM=R0*SUM(Y_i/M_i)

   !$OMP PARALLEL SHARED(RSUM)
   IF (N_SPECIES>0) THEN
      !$OMP DO COLLAPSE(3) PRIVATE(K,J,I,YY_GET)
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR   
               !!$ IF ((K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_MASS_DENS_15'
               YY_GET(:) = YY(I,J,K,:)
               CALL GET_SPECIFIC_GAS_CONSTANT(YY_GET,RSUM(I,J,K))
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO
      IF (ISOTHERMAL) THEN
         !$OMP DO COLLAPSE(3) PRIVATE(K,J,I)
         DO K=0,KBP1
            DO J=0,JBP1
               DO I=0,IBP1
                  !!$ IF ((K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_MASS_DENS_16'
                  RHO(I,J,K) = PBAR(K,PRESSURE_ZONE(I,J,K))/(TMPA*RSUM(I,J,K))
               ENDDO
            ENDDO
         ENDDO
         !$OMP END DO
      ENDIF
   ENDIF

   ! Extract predicted temperature at next time step from Equation of State

   IF (.NOT. ISOTHERMAL) THEN
      IF (N_SPECIES==0) THEN
         !$OMP DO COLLAPSE(3) PRIVATE(K,J,I)
         DO K=0,KBP1
            DO J=0,JBP1
               DO I=0,IBP1
                  !!$ IF ((K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_MASS_DENS_17'
                  TMP(I,J,K) = PBAR(K,PRESSURE_ZONE(I,J,K))/(RSUM0*RHO(I,J,K))
               ENDDO
            ENDDO
         ENDDO
         !$OMP END DO
      ELSE
         !$OMP DO COLLAPSE(3) PRIVATE(K,J,I)
         DO K=0,KBP1
            DO J=0,JBP1
               DO I=0,IBP1
                  !!$ IF ((K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_MASS_DENS_18'
                  TMP(I,J,K) = PBAR(K,PRESSURE_ZONE(I,J,K))/(RSUM(I,J,K)*RHO(I,J,K))
               ENDDO
            ENDDO
         ENDDO
         !$OMP END DO
      ENDIF
      !$OMP WORKSHARE
      TMP = MAX(TMPMIN,MIN(TMPMAX,TMP))
      !$OMP END WORKSHARE
   ENDIF
   !$OMP END PARALLEL
END SELECT PREDICTOR_STEP

TUSED(3,NM)=TUSED(3,NM)+SECOND()-TNOW
 
END SUBROUTINE DENSITY
 

SUBROUTINE CHECK_DENSITY
 
! Redistribute mass from cells below or above the density cut-off limits

USE GLOBAL_CONSTANTS, ONLY : PREDICTOR, CORRECTOR, N_SPECIES,RHOMIN,RHOMAX,DEBUG_OPENMP,TWO_D
REAL(EB) :: SUM,CONST,CONST2,RHOMI,RHOPI,RHOMJ,RHOPJ,RHOMK,RHOPK,RHO00,RMIN,RMAX
INTEGER  :: IC,ISUM,I,J,K
LOGICAL :: LC(-3:3)
REAL(EB), POINTER, DIMENSION(:,:,:) :: RHODELTA,V_CELL

RHODELTA => WORK2

IF (PREDICTOR) THEN
   RHOP=>RHOS
ELSE
   RHOP=>RHO
ENDIF

V_CELL => WORK3

!$OMP PARALLEL
IF (TWO_D) THEN
   !$OMP DO COLLAPSE(3) PRIVATE(K,J,I)
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            V_CELL(I,J,K) = DX(I)*DZ(K)
         ENDDO
      ENDDO
   ENDDO
   !$OMP END DO
ELSE
   !$OMP DO COLLAPSE(3) PRIVATE(K,J,I)
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            V_CELL(I,J,K) = DX(I)*DY(J)*DZ(K)
         ENDDO
      ENDDO
   ENDDO
   !$OMP END DO
ENDIF

 
! Correct undershoots

!$OMP WORKSHARE
RHODELTA = 0._EB
!$OMP END WORKSHARE

!$OMP DO COLLAPSE(3) PRIVATE(K,J,I,IC,RMIN,SUM,ISUM,LC,RHO00,RHOMI,RHOPI,RHOMJ,RHOPJ,RHOMK,RHOPK,CONST,CONST2)
DO K=1,KBAR
   DO J=1,JBAR
      CHECK_LOOP: DO I=1,IBAR
         !!$ IF ((K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_MASS_CHECK_DENSITY_01'
         IC = CELL_INDEX(I,J,K)
         IF (SOLID(IC)) CYCLE CHECK_LOOP
         RMIN = RHOMIN
         IF (RHOP(I,J,K)>=RMIN) CYCLE CHECK_LOOP
         SUM   = 0._EB
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
            IF(SUM-ISUM*RHO00 /= 0._EB) THEN
               CONST = (RHOMIN-RHO00)/(SUM-ISUM*RHO00)
               IF (LC(-1)) THEN
                  CONST2 = CONST*V_CELL(I,J,K)/V_CELL(I-1,J,K)
                  RHODELTA(I-1,J,K) = RHODELTA(I-1,J,K) + MAX(RMIN,RHOMI+CONST2*(RHO00-RHOMI)) - RHOP(I-1,J,K)
               ENDIF
               IF (LC( 1)) THEN
                  CONST2 = CONST*V_CELL(I,J,K)/V_CELL(I+1,J,K)
                  RHODELTA(I+1,J,K) = RHODELTA(I+1,J,K) + MAX(RMIN,RHOPI+CONST2*(RHO00-RHOPI)) - RHOP(I+1,J,K)
               ENDIF
               IF (LC(-2)) THEN
                  CONST2 = CONST*V_CELL(I,J,K)/V_CELL(I,J-1,K)
                  RHODELTA(I,J-1,K) = RHODELTA(I,J-1,K) + MAX(RMIN,RHOMJ+CONST2*(RHO00-RHOMJ)) - RHOP(I,J-1,K)
               ENDIF
               IF (LC( 2)) THEN
                  CONST2 = CONST*V_CELL(I,J,K)/V_CELL(I,J+1,K)
                  RHODELTA(I,J+1,K) = RHODELTA(I,J+1,K) + MAX(RMIN,RHOPJ+CONST2*(RHO00-RHOPJ)) - RHOP(I,J+1,K)
               ENDIF
               IF (LC(-3)) THEN
                  CONST2 = CONST*V_CELL(I,J,K)/V_CELL(I,J,K-1)
                  RHODELTA(I,J,K-1) = RHODELTA(I,J,K-1) + MAX(RMIN,RHOMK+CONST2*(RHO00-RHOMK)) - RHOP(I,J,K-1)
               ENDIF
               IF (LC( 3)) THEN
                  CONST2 = CONST*V_CELL(I,J,K)/V_CELL(I,J,K+1)
                  RHODELTA(I,J,K+1) = RHODELTA(I,J,K+1) + MAX(RMIN,RHOPK+CONST2*(RHO00-RHOPK)) - RHOP(I,J,K+1)
               ENDIF
               RHODELTA(I,J,K) = RHODELTA(I,J,K) + RMIN - RHOP(I,J,K)
            ENDIF
         ENDIF
      ENDDO CHECK_LOOP
   ENDDO
ENDDO
!$OMP END DO

!$OMP WORKSHARE
RHOP = MAX(RHOMIN,RHOP+RHODELTA)

! Correct overshoots

RHODELTA = 0._EB
!$OMP END WORKSHARE

!$OMP DO COLLAPSE(3) PRIVATE(K,J,I,IC,RMAX,SUM,ISUM,LC,RHO00,RHOMI,RHOPI,RHOMJ,RHOPJ,RHOMK,RHOPK,CONST,CONST2)
DO K=1,KBAR
   DO J=1,JBAR
      CHECK_LOOP2: DO I=1,IBAR
         !!$ IF ((K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_MASS_CHECK_DENSITY_02'
         IC = CELL_INDEX(I,J,K)
         IF (SOLID(IC)) CYCLE CHECK_LOOP2
         RMAX = RHOMAX
         IF (RHOP(I,J,K)<=RMAX) CYCLE CHECK_LOOP2
         SUM   = 0._EB
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
            IF(SUM-ISUM*RHO00 /= 0._EB) THEN         
               CONST = (RMAX-RHO00)/(SUM-ISUM*RHO00)
               IF (LC(-1)) THEN
                  CONST2 = CONST*V_CELL(I,J,K)/V_CELL(I-1,J,K)
                  RHODELTA(I-1,J,K) = RHODELTA(I-1,J,K) + MIN(RMAX,RHOMI+CONST2*(RHO00-RHOMI)) - RHOP(I-1,J,K)
               ENDIF
               IF (LC( 1)) THEN
                  CONST2 = CONST*V_CELL(I,J,K)/V_CELL(I+1,J,K)
                  RHODELTA(I+1,J,K) = RHODELTA(I+1,J,K) + MIN(RMAX,RHOPI+CONST2*(RHO00-RHOPI)) - RHOP(I+1,J,K)
               ENDIF
               IF (LC(-2)) THEN
                  CONST2 = CONST*V_CELL(I,J,K)/V_CELL(I,J-1,K)
                  RHODELTA(I,J-1,K) = RHODELTA(I,J-1,K) + MIN(RMAX,RHOMJ+CONST2*(RHO00-RHOMJ)) - RHOP(I,J-1,K)
               ENDIF
               IF (LC( 2)) THEN
                  CONST2 = CONST*V_CELL(I,J,K)/V_CELL(I,J+1,K)
                  RHODELTA(I,J+1,K) = RHODELTA(I,J+1,K) + MIN(RMAX,RHOPJ+CONST2*(RHO00-RHOPJ)) - RHOP(I,J+1,K)
               ENDIF
               IF (LC(-3)) THEN
                  CONST2 = CONST*V_CELL(I,J,K)/V_CELL(I,J,K-1)
                  RHODELTA(I,J,K-1) = RHODELTA(I,J,K-1) + MIN(RMAX,RHOMK+CONST2*(RHO00-RHOMK)) - RHOP(I,J,K-1)
               ENDIF
               IF (LC( 3)) THEN
                  CONST2 = CONST*V_CELL(I,J,K)/V_CELL(I,J,K+1)
                  RHODELTA(I,J,K+1) = RHODELTA(I,J,K+1) + MIN(RMAX,RHOPK+CONST2*(RHO00-RHOPK)) - RHOP(I,J,K+1)
               ENDIF
               RHODELTA(I,J,K) = RHODELTA(I,J,K) + RMAX - RHOP(I,J,K)
            ENDIF
         ENDIF
      ENDDO CHECK_LOOP2
   ENDDO
ENDDO
!$OMP END DO

!$OMP WORKSHARE
RHOP = MIN(RHOMAX,RHOP+RHODELTA)
!$OMP END WORKSHARE
!$OMP END PARALLEL
END SUBROUTINE CHECK_DENSITY
 
 
SUBROUTINE CHECK_MASS_FRACTION

! Redistribute species mass from cells below or above the cut-off limits

USE GLOBAL_CONSTANTS, ONLY : PREDICTOR, CORRECTOR, N_SPECIES,YYMIN,YYMAX,POROUS_BOUNDARY,DEBUG_OPENMP
REAL(EB) :: SUM,CONST,RHYMI,RHYPI,RHYMJ,RHYPJ,RHYMK,RHYPK,RHY0,YMI,YPI,YMJ,YPJ,YMK,YPK,Y00,YMIN,YMAX
INTEGER  :: IC,N,ISUM, IW_A(-3:3),I,J,K
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

   !$OMP PARALLEL
   !$OMP WORKSHARE
   YYDELTA = 0._EB
   !$OMP END WORKSHARE

   ! Do undershoots

   !$OMP DO COLLAPSE(3) PRIVATE(K,J,I,IC,IW_A,Y00,SUM,ISUM,LC,YMIN,YMI,YPI,YMJ,YPJ,YMK,YPK) &
   !$OMP PRIVATE(RHY0,RHYPI,RHYMI,RHYPJ,RHYMJ,RHYPK,RHYMK,CONST)
   DO K=1,KBAR
      DO J=1,JBAR
         CHECK_LOOP: DO I=1,IBAR
            !!$ IF ((K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_MASS_CHECK_M_FRACTION_01'
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
               IF (BOUNDARY_TYPE(IW_A(-1))==POROUS_BOUNDARY) THEN
                 YMI = YYP(I-1,J,K,N)
                 LC(-1) = .TRUE.
               ELSE
                 YMI = YY_F(IW_A(-1),N)  
               ENDIF
            ENDIF          
            IF (IW_A( 1) == 0) THEN
               YPI = YYP(I+1,J,K,N)
               LC( 1) = .TRUE.
            ELSE
               IF (BOUNDARY_TYPE(IW_A(1))==POROUS_BOUNDARY) THEN
                 YPI = YYP(I+1,J,K,N)
                 LC( 1) = .TRUE.
               ELSE
                 YPI = YY_F(IW_A(1),N)  
               ENDIF
            ENDIF           
            IF (IW_A(-2) == 0) THEN
               YMJ = YYP(I,J-1,K,N)
               LC(-2) = .TRUE.
            ELSE
               IF (BOUNDARY_TYPE(IW_A(-2))==POROUS_BOUNDARY) THEN
                 YMJ = YYP(I,J-1,K,N)
                 LC(-2) = .TRUE.
               ELSE
                 YMJ = YY_F(IW_A(-2),N)  
               ENDIF
            ENDIF         
            IF (IW_A( 2) == 0) THEN
               YPJ = YYP(I,J+1,K,N)
               LC( 2) = .TRUE.
            ELSE
               IF (BOUNDARY_TYPE(IW_A( 2))==POROUS_BOUNDARY) THEN
                 YPJ = YYP(I,J+1,K,N)
                 LC( 2) = .TRUE.
               ELSE
                 YPJ = YY_F(IW_A( 2),N)  
               ENDIF
            ENDIF         
            IF (IW_A(-3) == 0) THEN
               YMK = YYP(I,J,K-1,N)
               LC(-3) = .TRUE.
            ELSE
               IF (BOUNDARY_TYPE(IW_A(-3))==POROUS_BOUNDARY) THEN
                 YMK = YYP(I,J,K-1,N)
                 LC(-3) = .TRUE.
               ELSE
                 YMK = YY_F(IW_A(-3),N)  
               ENDIF
            ENDIF         
            IF (IW_A( 3) == 0) THEN
               YPK = YYP(I,J,K+1,N)
               LC( 3) = .TRUE.
            ELSE
               IF (BOUNDARY_TYPE(IW_A( 3))==POROUS_BOUNDARY) THEN
                 YPK = YYP(I,J,K+1,N)
                 LC( 3) = .TRUE.
               ELSE
                 YPK = YY_F(IW_A( 3),N)  
               ENDIF
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
                  IF (SUM/=0._EB) THEN
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
            ENDIF
         ENDDO CHECK_LOOP
      ENDDO
   ENDDO
   !$OMP END DO
   
   !$OMP WORKSHARE
   YYP(:,:,:,N) = YYP(:,:,:,N) + YYDELTA
   YYDELTA=0._EB
   !$OMP END WORKSHARE

   ! Do overshoots
   !$OMP DO COLLAPSE(3) PRIVATE(K,J,I,IC,IW_A,Y00,SUM,ISUM,LC,YMIN,YMI,YPI,YMK,YPK,YMJ,YPJ,YMAX) &
   !$OMP PRIVATE(RHY0,RHYPI,RHYMI,RHYPJ,RHYMJ,RHYPK,RHYMK,CONST)
   DO K=1,KBAR
      DO J=1,JBAR
         CHECK_LOOP2: DO I=1,IBAR
            !!$ IF ((K == 1) .AND. (J == 1) .AND. (I == 1) .AND. DEBUG_OPENMP) WRITE(*,*) 'OpenMP_MASS_CHECK_M_FRACTION_02'
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
               YMI = YY_F(IW_A(-1),N)  
            ENDIF          
            IF (IW_A( 1) == 0) THEN
               YPI = YYP(I+1,J,K,N)
               LC( 1) = .TRUE.
            ELSE
               YPI = YY_F(IW_A( 1),N)  
            ENDIF           
            IF (IW_A(-2) == 0) THEN
               YMJ = YYP(I,J-1,K,N)
               LC(-2) = .TRUE.
            ELSE
               YMJ = YY_F(IW_A(-2),N)  
            ENDIF         
            IF (IW_A( 2) == 0) THEN
               YPJ = YYP(I,J+1,K,N)
               LC( 2) = .TRUE.
            ELSE
               YPJ = YY_F(IW_A( 2),N)  
            ENDIF         
            IF (IW_A(-3) == 0) THEN
               YMK = YYP(I,J,K-1,N)
               LC(-3) = .TRUE.
            ELSE
               YMK = YY_F(IW_A(-3),N)  
            ENDIF         
            IF (IW_A( 3) == 0) THEN
               YPK = YYP(I,J,K+1,N)
               LC( 3) = .TRUE.
            ELSE
               YPK = YY_F(IW_A( 3),N)  
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
                  IF (SUM/=0._EB) THEN
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
            ENDIF
         ENDDO CHECK_LOOP2
      ENDDO
   ENDDO  
   !$OMP END DO 

   !$OMP WORKSHARE
   YYP(:,:,:,N) = YYP(:,:,:,N) + YYDELTA   
   !$OMP END WORKSHARE
   !$OMP END PARALLEL
ENDDO SPECIESLOOP

RETURN

END SUBROUTINE CHECK_MASS_FRACTION


!===========================================================================
! The following are experimental scalar transport routines which are invoked
! by setting FLUX_LIMITER = {0,1,2,3,4} on the MISC line.  See the function
! SCALAR_FACE_VALUE below for a description of the FLUX_LIMITER value. ~RJM
!===========================================================================

SUBROUTINE DENSITY_TVD(NM)

! Update the density and species mass fractions

USE COMP_FUNCTIONS, ONLY: SECOND 
USE PHYSICAL_FUNCTIONS, ONLY : GET_SPECIFIC_GAS_CONSTANT
USE GLOBAL_CONSTANTS, ONLY: N_SPECIES,CO_PRODUCTION,I_PROG_F,I_PROG_CO,I_FUEL,TMPMAX,TMPMIN,EVACUATION_ONLY,PREDICTOR,CORRECTOR, &
                            CHANGE_TIME_STEP,ISOTHERMAL,TMPA, N_ZONE,MIXTURE_FRACTION_SPECIES, LU_ERR, &
                            GAS_SPECIES, MIXTURE_FRACTION,R0,SOLID_PHASE_ONLY,TUSED,FLUX_LIMITER, &
                            RSUM0,CHECK_KINETIC_ENERGY,CYLINDRICAL,CLIP_MASS_FRACTION
REAL(EB) :: TNOW,YY_GET(1:N_SPECIES)
INTEGER  :: I,J,K,N
INTEGER, INTENT(IN) :: NM
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: RHOYYP,FX,FY,FZ,YYN

 
IF (EVACUATION_ONLY(NM)) RETURN
IF (SOLID_PHASE_ONLY) RETURN

TNOW=SECOND()
CALL POINT_TO_MESH(NM)

FX => SCALAR_SAVE1
FY => SCALAR_SAVE2
FZ => SCALAR_SAVE3
RHOYYP => SCALAR_SAVE4
YYN => SCALAR_SAVE5

SELECT_SUBSTEP: IF (PREDICTOR) THEN
   
   ! Update the density
   
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            RHOS(I,J,K) = RHO(I,J,K) - DT*( RDX(I)*(FX(I,J,K,0)-FX(I-1,J,K,0))*RRN(I) &
                                          + RDY(J)*(FY(I,J,K,0)-FY(I,J-1,K,0)) &
                                          + RDZ(K)*(FZ(I,J,K,0)-FZ(I,J,K-1,0)) )
         ENDDO
      ENDDO
   ENDDO
   
   ! Correct densities above or below clip limits

   CALL CHECK_DENSITY
   
   ! Update mass fractions
   
   DO N=1,N_SPECIES
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               YYN(I,J,K,N) = RHOYYP(I,J,K,N)   
               YYS(I,J,K,N) = RHOYYP(I,J,K,N) - DT*( RDX(I)*(FX(I,J,K,N)-FX(I-1,J,K,N))*RRN(I) &
                                                   + RDY(J)*(FY(I,J,K,N)-FY(I,J-1,K,N)) &
                                                   + RDZ(K)*(FZ(I,J,K,N)-FZ(I,J,K-1,N)) )
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
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   
   IF (CLIP_MASS_FRACTION) THEN
      DO N=1,N_SPECIES
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  YYS(I,J,K,N) = MIN(YYS(I,J,K,N),1._EB)
                  YYS(I,J,K,N) = MAX(YYS(I,J,K,N),0._EB)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ELSE
      CALL CHECK_MASS_FRACTION
   ENDIF

   ! Predict background pressure at next time step

   DO I=1,N_ZONE
      PBAR_S(:,I) = PBAR(:,I) + D_PBAR_DT(I)*DT
   ENDDO

   ! Compute molecular weight term RSUM=R0*SUM(Y_i/M_i)
   
   IF (N_SPECIES>0) THEN
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR       
               YY_GET(:) = YYS(I,J,K,:)
               CALL GET_SPECIFIC_GAS_CONSTANT(YY_GET,RSUM(I,J,K))
            ENDDO
         ENDDO
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

   ! Extract predicted temperature at next time step from Equation of State
   
   IF (.NOT. ISOTHERMAL) THEN
      IF (N_SPECIES==0) THEN
         DO K=0,KBP1
            DO J=0,JBP1
               DO I=0,IBP1
                  TMP(I,J,K) = PBAR_S(K,PRESSURE_ZONE(I,J,K))/(RSUM0*RHOS(I,J,K))
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

ELSEIF (CORRECTOR) THEN

   ! Update the density
   
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            RHOS(I,J,K) = RHOS(I,J,K) - DT*( RDX(I)*(FX(I,J,K,0)-FX(I-1,J,K,0))*RRN(I) &
                                           + RDY(J)*(FY(I,J,K,0)-FY(I,J-1,K,0)) &
                                           + RDZ(K)*(FZ(I,J,K,0)-FZ(I,J,K-1,0)) )
         ENDDO
      ENDDO
   ENDDO
   
   ! Update mass fractions
   
   DO N=1,N_SPECIES
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR             
               YYS(I,J,K,N) = RHOYYP(I,J,K,N) - DT*( RDX(I)*(FX(I,J,K,N)-FX(I-1,J,K,N))*RRN(I) &
                                                   + RDY(J)*(FY(I,J,K,N)-FY(I,J-1,K,N)) &
                                                   + RDZ(K)*(FZ(I,J,K,N)-FZ(I,J,K-1,N)) )
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   
   ! Corrector step
   
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            RHO(I,J,K) = 0.5_EB*( RHO(I,J,K) + RHOS(I,J,K) )
         ENDDO
      ENDDO
   ENDDO
   
   ! Correct densities above or below clip limits

   CALL CHECK_DENSITY
   
   DO N=1,N_SPECIES
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               YY(I,J,K,N) = 0.5_EB*( YYN(I,J,K,N) + YYS(I,J,K,N) )/RHO(I,J,K)
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   
   IF (CLIP_MASS_FRACTION) THEN
      DO N=1,N_SPECIES
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  YY(I,J,K,N) = MAX(YY(I,J,K,N),0._EB)
                  YY(I,J,K,N) = MIN(YY(I,J,K,N),1._EB)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ELSE
      CALL CHECK_MASS_FRACTION
   ENDIF
   
   ! Correct background pressure

   DO I=1,N_ZONE
      PBAR(:,I) = .5_EB*(PBAR(:,I) + PBAR_S(:,I) + D_PBAR_S_DT(I)*DT)
   ENDDO
 
   ! Compute molecular weight term RSUM=R0*SUM(Y_i/M_i)
   
   IF (N_SPECIES>0) THEN
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR   
               YY_GET(:) = YY(I,J,K,:)
               CALL GET_SPECIFIC_GAS_CONSTANT(YY_GET,RSUM(I,J,K))
            ENDDO
         ENDDO
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

   ! Extract predicted temperature at next time step from Equation of State
   
   IF (.NOT. ISOTHERMAL) THEN
      IF (N_SPECIES==0) THEN
         DO K=0,KBP1
            DO J=0,JBP1
               DO I=0,IBP1
                  TMP(I,J,K) = PBAR(K,PRESSURE_ZONE(I,J,K))/(RSUM0*RHO(I,J,K))
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
   
ENDIF SELECT_SUBSTEP

TUSED(3,NM)=TUSED(3,NM)+SECOND()-TNOW
 
END SUBROUTINE DENSITY_TVD


SUBROUTINE SCALARF(NM)
USE COMP_FUNCTIONS, ONLY: SECOND
USE GLOBAL_CONSTANTS, ONLY: N_SPECIES,PREDICTOR,CORRECTOR,FLUX_LIMITER,NULL_BOUNDARY,POROUS_BOUNDARY,OPEN_BOUNDARY, &
                            LU_ERR,INTERPOLATED_BOUNDARY,MYID,NO_MASS_FLUX,SOLID_BOUNDARY,INTERPOLATED_SURF_INDEX,  &
                            OPEN_SURF_INDEX,MIRROR_BOUNDARY,MIRROR_SURF_INDEX,SPECIFIED_MASS_FLUX,INFLOW_OUTFLOW,   &
                            CYLINDRICAL,INTERPOLATED_BC,TUSED,EMB_BC

! Computes the scalar advective and diffusive flux
INTEGER, INTENT(IN) :: NM
INTEGER :: I,J,K,N,II,JJ,KK,IOR,IW,IIG,JJG,KKG,IBC,METHOD_ID,NOM,IIO,JJO,KKO ! ICM,ICP,
REAL(EB) :: ZZ(4),TNOW
REAL(EB), POINTER, DIMENSION(:,:,:) :: RHOP,UU,VV,WW
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: RHOYYP,YYP,FX,FY,FZ
TYPE (SURFACE_TYPE), POINTER :: SF
TYPE (OMESH_TYPE), POINTER :: M2

! Notes:
!
! At the moment (FDS 5.4.2) it is critical that SCALARF stay WITHIN the CHANGE_TIME_STEP loop.  I am writing
! myself this note because everytime I go back and try to debug something I wonder why SCALARF just does not
! sit outside this loop in the same way as MASS_FINITE_DIFFERENCES.  The reason is that DIVERGENCE_PART_1
! initializes the FX (flux) arrays and adds the diffusive flux to this array.  So, if SCALARF is computed
! outside the CHANGE_TIME_STEP loop and we happen to iterate on the loop then we clear the advective flux
! and we are left with only the diffusive flux.

TNOW=SECOND()
CALL POINT_TO_MESH(NM)

FX => SCALAR_SAVE1
FY => SCALAR_SAVE2
FZ => SCALAR_SAVE3
RHOYYP => SCALAR_SAVE4
!! eliminate diffusion for testing
!FX = 0._EB
!FY = 0._EB
!FZ = 0._EB

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

! Density flux

DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBM1
         ZZ(1:4) = RHOP(I-1:I+2,J,K)
         FX(I,J,K,0) = R(I)*UU(I,J,K)*SCALAR_FACE_VALUE(UU(I,J,K),ZZ,FLUX_LIMITER)
      ENDDO
   ENDDO
ENDDO

DO K=1,KBAR
   DO J=1,JBM1
      DO I=1,IBAR
         ZZ(1:4) = RHOP(I,J-1:J+2,K)
         FY(I,J,K,0) = VV(I,J,K)*SCALAR_FACE_VALUE(VV(I,J,K),ZZ,FLUX_LIMITER)
      ENDDO
   ENDDO
ENDDO

DO K=1,KBM1
   DO J=1,JBAR
      DO I=1,IBAR
         ZZ(1:4) = RHOP(I,J,K-1:K+2)
         FZ(I,J,K,0) = WW(I,J,K)*SCALAR_FACE_VALUE(WW(I,J,K),ZZ,FLUX_LIMITER)
      ENDDO
   ENDDO
ENDDO

! Compute mass fluxes at boundaries

WALL_LOOP: DO IW=1,NWC
   IF (BOUNDARY_TYPE(IW)==NULL_BOUNDARY .OR. BOUNDARY_TYPE(IW)==POROUS_BOUNDARY) CYCLE WALL_LOOP
   II  = IJKW(1,IW)
   JJ  = IJKW(2,IW)
   KK  = IJKW(3,IW)
   IBC = IJKW(5,IW)
   IIG = IJKW(6,IW)
   JJG = IJKW(7,IW)
   KKG = IJKW(8,IW)
   IOR = IJKW(4,IW)
   IF (BOUNDARY_TYPE(IW)/=SOLID_BOUNDARY) THEN
      IF (BOUNDARY_TYPE(IW)==INTERPOLATED_BOUNDARY) IBC = INTERPOLATED_SURF_INDEX
      IF (BOUNDARY_TYPE(IW)==MIRROR_BOUNDARY)       IBC = MIRROR_SURF_INDEX
   ENDIF
   SF => SURFACE(IBC)
   METHOD_ID = SF%SPECIES_BC_INDEX
   IF (BOUNDARY_TYPE(IW)==OPEN_BOUNDARY) METHOD_ID = INFLOW_OUTFLOW
   IF (BOUNDARY_TYPE(IW)==INTERPOLATED_BOUNDARY .AND. MESH_LEVEL/=0) METHOD_ID = EMB_BC
   IF (BOUNDARY_TYPE(IW)==SOLID_BOUNDARY .AND. IJKW(9,IW)>0) METHOD_ID = 999
   
   ! Apply the different species boundary conditions
   METHOD_OF_MASS_TRANSFER: SELECT CASE(METHOD_ID)

      CASE (INTERPOLATED_BC) METHOD_OF_MASS_TRANSFER
         NOM = IJKW(9,IW)
         IIO = IJKW(10,IW)
         JJO = IJKW(11,IW)
         KKO = IJKW(12,IW)
         M2  => OMESH(NOM)
      
         SELECT CASE(IOR)
            CASE( 1)
               IF (PREDICTOR) ZZ(1) = M2%RHO(IIO-1,JJO,KKO)
               IF (CORRECTOR) ZZ(1) = M2%RHOS(IIO-1,JJO,KKO)
               IF (PREDICTOR) ZZ(2) = M2%RHO(IIO,JJO,KKO)
               IF (CORRECTOR) ZZ(2) = M2%RHOS(IIO,JJO,KKO)
               ZZ(3) = RHOP(IIG,JJG,KKG)
               ZZ(4) = RHOP(IIG+1,JJG,KKG)
               !debug...
               !ZZ(1) = RHOP(IBM1,JJG,KKG)
               !ZZ(2) = RHOP(IBAR,JJG,KKG)
               !ZZ(3) = RHOP(1,JJG,KKG)
               !ZZ(4) = RHOP(2,JJG,KKG)
               FX(II,JJ,KK,0) = UVW_SAVE(IW)*SCALAR_FACE_VALUE(UVW_SAVE(IW),ZZ,FLUX_LIMITER)*R(II)
            CASE(-1)
               ZZ(1) = RHOP(IIG-1,JJG,KKG)
               ZZ(2) = RHOP(IIG,JJG,KKG)
               IF (PREDICTOR) ZZ(3) = M2%RHO(IIO,JJO,KKO)
               IF (CORRECTOR) ZZ(3) = M2%RHOS(IIO,JJO,KKO)
               IF (PREDICTOR) ZZ(4) = M2%RHO(IIO+1,JJO,KKO)
               IF (CORRECTOR) ZZ(4) = M2%RHOS(IIO+1,JJO,KKO)
               !debug...
               !ZZ(1) = RHOP(IBM1,JJG,KKG)
               !ZZ(2) = RHOP(IBAR,JJG,KKG)
               !ZZ(3) = RHOP(1,JJG,KKG)
               !ZZ(4) = RHOP(2,JJG,KKG)
               FX(II-1,JJ,KK,0) = UVW_SAVE(IW)*SCALAR_FACE_VALUE(UVW_SAVE(IW),ZZ,FLUX_LIMITER)*R(II-1)
            CASE( 2)
               IF (PREDICTOR) ZZ(1) = M2%RHO(IIO,JJO-1,KKO)
               IF (CORRECTOR) ZZ(1) = M2%RHOS(IIO,JJO-1,KKO)
               IF (PREDICTOR) ZZ(2) = M2%RHO(IIO,JJO,KKO)
               IF (CORRECTOR) ZZ(2) = M2%RHOS(IIO,JJO,KKO)
               ZZ(3) = RHOP(IIG,JJG,KKG)
               ZZ(4) = RHOP(IIG,JJG+1,KKG)
               FY(II,JJ,KK,0) = UVW_SAVE(IW)*SCALAR_FACE_VALUE(UVW_SAVE(IW),ZZ,FLUX_LIMITER)
            CASE(-2)
               ZZ(1) = RHOP(IIG,JJG-1,KKG)
               ZZ(2) = RHOP(IIG,JJG,KKG)
               IF (PREDICTOR) ZZ(3) = M2%RHO(IIO,JJO,KKO)
               IF (CORRECTOR) ZZ(3) = M2%RHOS(IIO,JJO,KKO)
               IF (PREDICTOR) ZZ(4) = M2%RHO(IIO,JJO+1,KKO)
               IF (CORRECTOR) ZZ(4) = M2%RHOS(IIO,JJO+1,KKO)
               FY(II,JJ-1,KK,0) = UVW_SAVE(IW)*SCALAR_FACE_VALUE(UVW_SAVE(IW),ZZ,FLUX_LIMITER)
            CASE( 3)
               IF (PREDICTOR) ZZ(1) = M2%RHO(IIO,JJO,KKO-1)
               IF (CORRECTOR) ZZ(1) = M2%RHOS(IIO,JJO,KKO-1)
               IF (PREDICTOR) ZZ(2) = M2%RHO(IIO,JJO,KKO)
               IF (CORRECTOR) ZZ(2) = M2%RHOS(IIO,JJO,KKO)
               ZZ(3) = RHOP(IIG,JJG,KKG)
               ZZ(4) = RHOP(IIG,JJG,KKG+1)
               !degug...
               !ZZ(1) = RHOP(IIG,JJG,KBM1)
               !ZZ(2) = RHOP(IIG,JJG,KBAR)
               !ZZ(3) = RHOP(IIG,JJG,1)
               !ZZ(4) = RHOP(IIG,JJG,2)
               FZ(II,JJ,KK,0) = UVW_SAVE(IW)*SCALAR_FACE_VALUE(UVW_SAVE(IW),ZZ,FLUX_LIMITER)
            CASE(-3)
               ZZ(1) = RHOP(IIG,JJG,KKG-1)
               ZZ(2) = RHOP(IIG,JJG,KKG)
               IF (PREDICTOR) ZZ(3) = M2%RHO(IIO,JJO,KKO)
               IF (CORRECTOR) ZZ(3) = M2%RHOS(IIO,JJO,KKO)
               IF (PREDICTOR) ZZ(4) = M2%RHO(IIO,JJO,KKO+1)
               IF (CORRECTOR) ZZ(4) = M2%RHOS(IIO,JJO,KKO+1)
               !debug...
               !ZZ(1) = RHOP(IIG,JJG,KBM1)
               !ZZ(2) = RHOP(IIG,JJG,KBAR)
               !ZZ(3) = RHOP(IIG,JJG,1)
               !ZZ(4) = RHOP(IIG,JJG,2)
               FZ(II,JJ,KK-1,0) = UVW_SAVE(IW)*SCALAR_FACE_VALUE(UVW_SAVE(IW),ZZ,FLUX_LIMITER)
         END SELECT
         
      CASE (999) METHOD_OF_MASS_TRANSFER

         SELECT CASE(IOR)
            CASE( 1)
               FX(II,JJ,KK,0)   = UVW_SAVE(IW)*RHO_F(IW)*R(II)
            CASE(-1)
               FX(II-1,JJ,KK,0) = UVW_SAVE(IW)*RHO_F(IW)*R(II-1)
            CASE( 2)  
               FY(II,JJ,KK,0)   = UVW_SAVE(IW)*RHO_F(IW)
            CASE(-2)
               FY(II,JJ-1,KK,0) = UVW_SAVE(IW)*RHO_F(IW)
            CASE( 3)
               FZ(II,JJ,KK,0)   = UVW_SAVE(IW)*RHO_F(IW)
            CASE(-3)
               FZ(II,JJ,KK-1,0) = UVW_SAVE(IW)*RHO_F(IW)
         END SELECT
            
      CASE DEFAULT METHOD_OF_MASS_TRANSFER

         SELECT CASE(IOR)
            CASE( 1)
               FX(II,JJ,KK,0)   = UU(II,JJ,KK)*RHO_F(IW)*R(II)
            CASE(-1)
               FX(II-1,JJ,KK,0) = UU(II-1,JJ,KK)*RHO_F(IW)*R(II-1)
            CASE( 2)  
               FY(II,JJ,KK,0)   = VV(II,JJ,KK)*RHO_F(IW)
            CASE(-2)
               FY(II,JJ-1,KK,0) = VV(II,JJ-1,KK)*RHO_F(IW)
            CASE( 3)
               FZ(II,JJ,KK,0)   = WW(II,JJ,KK)*RHO_F(IW)
            CASE(-3)
               FZ(II,JJ,KK-1,0) = WW(II,JJ,KK-1)*RHO_F(IW)
         END SELECT
      
   END SELECT METHOD_OF_MASS_TRANSFER
      
ENDDO WALL_LOOP

! Compute species advective fluxes at interior cell faces and 
! add to diffusive flux which is already stored in FX,FY,FZ.

DO N=1,N_SPECIES
   DO K=0,KBP1
      DO J=0,JBP1
         DO I=0,IBP1
            RHOYYP(I,J,K,N) = RHOP(I,J,K)*YYP(I,J,K,N)
         ENDDO
      ENDDO
   ENDDO
ENDDO

SPECIES_LOOP: DO N=1,N_SPECIES

   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBM1
            ZZ(1:4) = RHOYYP(I-1:I+2,J,K,N)
            FX(I,J,K,N) = R(I)*(FX(I,J,K,N) + UU(I,J,K)*SCALAR_FACE_VALUE(UU(I,J,K),ZZ,FLUX_LIMITER))
         ENDDO
      ENDDO
   ENDDO

   DO K=1,KBAR
      DO J=1,JBM1
         DO I=1,IBAR
            ZZ(1:4) = RHOYYP(I,J-1:J+2,K,N)
            FY(I,J,K,N) = FY(I,J,K,N) + VV(I,J,K)*SCALAR_FACE_VALUE(VV(I,J,K),ZZ,FLUX_LIMITER)
         ENDDO
      ENDDO
   ENDDO

   DO K=1,KBM1
      DO J=1,JBAR
         DO I=1,IBAR
            ZZ(1:4) = RHOYYP(I,J,K-1:K+2,N)
            FZ(I,J,K,N) = FZ(I,J,K,N) + WW(I,J,K)*SCALAR_FACE_VALUE(WW(I,J,K),ZZ,FLUX_LIMITER)
         ENDDO
      ENDDO
   ENDDO

   ! Compute species fluxes at boundaries

   WALL_LOOP2: DO IW=1,NWC
      IF (BOUNDARY_TYPE(IW)==NULL_BOUNDARY .OR. BOUNDARY_TYPE(IW)==POROUS_BOUNDARY) CYCLE WALL_LOOP2
      II  = IJKW(1,IW)
      JJ  = IJKW(2,IW)
      KK  = IJKW(3,IW)
      IOR = IJKW(4,IW)
      IBC = IJKW(5,IW)
      IIG = IJKW(6,IW)
      JJG = IJKW(7,IW)
      KKG = IJKW(8,IW)
      IF (BOUNDARY_TYPE(IW)/=SOLID_BOUNDARY) THEN
         IF (BOUNDARY_TYPE(IW)==INTERPOLATED_BOUNDARY) IBC = INTERPOLATED_SURF_INDEX
         IF (BOUNDARY_TYPE(IW)==MIRROR_BOUNDARY)       IBC = MIRROR_SURF_INDEX
      ENDIF
      SF => SURFACE(IBC)
      METHOD_ID = SF%SPECIES_BC_INDEX
      IF (BOUNDARY_TYPE(IW)==OPEN_BOUNDARY) METHOD_ID = INFLOW_OUTFLOW
      IF (BOUNDARY_TYPE(IW)==INTERPOLATED_BOUNDARY .AND. MESH_LEVEL/=0) METHOD_ID = EMB_BC
      !!IF (BOUNDARY_TYPE(IW)==SOLID_BOUNDARY .AND. IJKW(9,IW)>0) METHOD_ID = 999
      
      ! Apply the different species boundary conditions
      METHOD_OF_MASS_TRANSFER2: SELECT CASE(METHOD_ID)
      
         CASE (NO_MASS_FLUX) METHOD_OF_MASS_TRANSFER2

            SELECT CASE(IOR)
               CASE( 1)
                  FX(II,JJ,KK,N)   = 0._EB
               CASE(-1)
                  FX(II-1,JJ,KK,N) = 0._EB
               CASE( 2)   
                  FY(II,JJ,KK,N)   = 0._EB
               CASE(-2)
                  FY(II,JJ-1,KK,N) = 0._EB
               CASE( 3)
                  FZ(II,JJ,KK,N)   = 0._EB
               CASE(-3)
                  FZ(II,JJ,KK-1,N) = 0._EB
            END SELECT
            
         CASE (SPECIFIED_MASS_FLUX) METHOD_OF_MASS_TRANSFER2
         
            SELECT CASE(IOR)
               CASE( 1)
                  FX(II,JJ,KK,N)   = MASSFLUX(IW,N)*R(II)
               CASE(-1)
                  FX(II-1,JJ,KK,N) = -MASSFLUX(IW,N)*R(II-1)
               CASE( 2)   
                  FY(II,JJ,KK,N)   = MASSFLUX(IW,N)
               CASE(-2)
                  FY(II,JJ-1,KK,N) = -MASSFLUX(IW,N)
               CASE( 3)
                  FZ(II,JJ,KK,N)   = MASSFLUX(IW,N)
               CASE(-3)
                  FZ(II,JJ,KK-1,N) = -MASSFLUX(IW,N)
            END SELECT
            
         CASE (INTERPOLATED_BC) METHOD_OF_MASS_TRANSFER2
            NOM = IJKW(9,IW)
            IIO = IJKW(10,IW)
            JJO = IJKW(11,IW)
            KKO = IJKW(12,IW)
            M2  => OMESH(NOM)
               
            SELECT CASE(IOR)
               CASE( 1)
                  IF (PREDICTOR) ZZ(1) = M2%RHO(IIO-1,JJO,KKO)*M2%YY(IIO-1,JJO,KKO,N)
                  IF (CORRECTOR) ZZ(1) = M2%RHOS(IIO-1,JJO,KKO)*M2%YYS(IIO-1,JJO,KKO,N)
                  IF (PREDICTOR) ZZ(2) = M2%RHO(IIO,JJO,KKO)*M2%YY(IIO,JJO,KKO,N)
                  IF (CORRECTOR) ZZ(2) = M2%RHOS(IIO,JJO,KKO)*M2%YYS(IIO,JJO,KKO,N)                  
                  ZZ(3) = RHOYYP(IIG,JJG,KKG,N)
                  ZZ(4) = RHOYYP(IIG+1,JJG,KKG,N)
                  FX(II,JJ,KK,N) = (FX(II,JJ,KK,N) + UU(II,JJ,KK)*SCALAR_FACE_VALUE(UU(II,JJ,KK),ZZ,FLUX_LIMITER))*R(II)
               CASE(-1)
                  ZZ(1) = RHOYYP(IIG-1,JJG,KKG,N)
                  ZZ(2) = RHOYYP(IIG,JJG,KKG,N)
                  IF (PREDICTOR) ZZ(3) = M2%RHO(IIO,JJO,KKO)*M2%YY(IIO,JJO,KKO,N)
                  IF (CORRECTOR) ZZ(3) = M2%RHOS(IIO,JJO,KKO)*M2%YYS(IIO,JJO,KKO,N)
                  IF (PREDICTOR) ZZ(4) = M2%RHO(IIO+1,JJO,KKO)*M2%YY(IIO+1,JJO,KKO,N)
                  IF (CORRECTOR) ZZ(4) = M2%RHOS(IIO+1,JJO,KKO)*M2%YYS(IIO+1,JJO,KKO,N)
                  FX(II-1,JJ,KK,N) = (FX(II-1,JJ,KK,N) + UU(II-1,JJ,KK)*SCALAR_FACE_VALUE(UU(II-1,JJ,KK),ZZ,FLUX_LIMITER))*R(II-1)
               CASE( 2)
                  IF (PREDICTOR) ZZ(1) = M2%RHO(IIO,JJO-1,KKO)*M2%YY(IIO,JJO-1,KKO,N)
                  IF (CORRECTOR) ZZ(1) = M2%RHOS(IIO,JJO-1,KKO)*M2%YYS(IIO,JJO-1,KKO,N)
                  IF (PREDICTOR) ZZ(2) = M2%RHO(IIO,JJO,KKO)*M2%YY(IIO,JJO,KKO,N)
                  IF (CORRECTOR) ZZ(2) = M2%RHOS(IIO,JJO,KKO)*M2%YYS(IIO,JJO,KKO,N)
                  ZZ(3) = RHOYYP(IIG,JJG,KKG,N)
                  ZZ(4) = RHOYYP(IIG,JJG+1,KKG,N)
                  FY(II,JJ,KK,N) = FY(II,JJ,KK,N) + VV(II,JJ,KK)*SCALAR_FACE_VALUE(VV(II,JJ,KK),ZZ,FLUX_LIMITER)
               CASE(-2)
                  ZZ(1) = RHOYYP(IIG,JJG-1,KKG,N)
                  ZZ(2) = RHOYYP(IIG,JJG,KKG,N)
                  IF (PREDICTOR) ZZ(3) = M2%RHO(IIO,JJO,KKO)*M2%YY(IIO,JJO,KKO,N)
                  IF (CORRECTOR) ZZ(3) = M2%RHOS(IIO,JJO,KKO)*M2%YYS(IIO,JJO,KKO,N)
                  IF (PREDICTOR) ZZ(4) = M2%RHO(IIO,JJO+1,KKO)*M2%YY(IIO,JJO+1,KKO,N)
                  IF (CORRECTOR) ZZ(4) = M2%RHOS(IIO,JJO+1,KKO)*M2%YYS(IIO,JJO+1,KKO,N)
                  FY(II,JJ-1,KK,N) = FY(II,JJ-1,KK,N) + VV(II,JJ-1,KK)*SCALAR_FACE_VALUE(VV(II,JJ-1,KK),ZZ,FLUX_LIMITER)
               CASE( 3)
                  IF (PREDICTOR) ZZ(1) = M2%RHO(IIO,JJO,KKO-1)*M2%YY(IIO,JJO,KKO-1,N)
                  IF (CORRECTOR) ZZ(1) = M2%RHOS(IIO,JJO,KKO-1)*M2%YYS(IIO,JJO,KKO-1,N)
                  IF (PREDICTOR) ZZ(2) = M2%RHO(IIO,JJO,KKO)*M2%YY(IIO,JJO,KKO,N)
                  IF (CORRECTOR) ZZ(2) = M2%RHOS(IIO,JJO,KKO)*M2%YYS(IIO,JJO,KKO,N)
                  ZZ(3) = RHOYYP(IIG,JJG,KKG,N)
                  ZZ(4) = RHOYYP(IIG,JJG,KKG+1,N)
                  FZ(II,JJ,KK,N) = FZ(II,JJ,KK,N) + WW(II,JJ,KK)*SCALAR_FACE_VALUE(WW(II,JJ,KK),ZZ,FLUX_LIMITER)
               CASE(-3)
                  ZZ(1) = RHOYYP(IIG,JJG,KKG-1,N)
                  ZZ(2) = RHOYYP(IIG,JJG,KKG,N)
                  IF (PREDICTOR) ZZ(3) = M2%RHO(IIO,JJO,KKO)*M2%YY(IIO,JJO,KKO,N)
                  IF (CORRECTOR) ZZ(3) = M2%RHOS(IIO,JJO,KKO)*M2%YYS(IIO,JJO,KKO,N)
                  IF (PREDICTOR) ZZ(4) = M2%RHO(IIO,JJO,KKO+1)*M2%YY(IIO,JJO,KKO+1,N)
                  IF (CORRECTOR) ZZ(4) = M2%RHOS(IIO,JJO,KKO+1)*M2%YYS(IIO,JJO,KKO+1,N)
                  FZ(II,JJ,KK-1,N) = FZ(II,JJ,KK-1,N) + WW(II,JJ,KK-1)*SCALAR_FACE_VALUE(WW(II,JJ,KK-1),ZZ,FLUX_LIMITER)
            END SELECT
            
!         CASE (999) METHOD_OF_MASS_TRANSFER2
!
!            SELECT CASE(IOR)
!               CASE( 1)
!                  FX(II,JJ,KK,N)   =(FW(IW,N) + UVW_SAVE(IW)*RHO_F(IW)*YY_F(IW,N))*R(II)
!               CASE(-1)
!                  FX(II-1,JJ,KK,N) =(FW(IW,N) + UVW_SAVE(IW)*RHO_F(IW)*YY_F(IW,N))*R(II-1)
!               CASE( 2)   
!                  FY(II,JJ,KK,N)   = FW(IW,N) + UVW_SAVE(IW)*RHO_F(IW)*YY_F(IW,N)
!               CASE(-2)
!                  FY(II,JJ-1,KK,N) = FW(IW,N) + UVW_SAVE(IW)*RHO_F(IW)*YY_F(IW,N)
!               CASE( 3)
!                  FZ(II,JJ,KK,N)   = FW(IW,N) + UVW_SAVE(IW)*RHO_F(IW)*YY_F(IW,N)
!               CASE(-3)
!                  FZ(II,JJ,KK-1,N) = FW(IW,N) + UVW_SAVE(IW)*RHO_F(IW)*YY_F(IW,N)
!            END SELECT
            
         CASE DEFAULT METHOD_OF_MASS_TRANSFER2

            SELECT CASE(IOR)
               CASE( 1)
                  FX(II,JJ,KK,N)   =(FW(IW,N) + UU(II,JJ,KK)  *RHO_F(IW)*YY_F(IW,N))*R(II)
               CASE(-1)
                  FX(II-1,JJ,KK,N) =(FW(IW,N) + UU(II-1,JJ,KK)*RHO_F(IW)*YY_F(IW,N))*R(II-1)
               CASE( 2)   
                  FY(II,JJ,KK,N)   = FW(IW,N) + VV(II,JJ,KK)  *RHO_F(IW)*YY_F(IW,N)
               CASE(-2)
                  FY(II,JJ-1,KK,N) = FW(IW,N) + VV(II,JJ-1,KK)*RHO_F(IW)*YY_F(IW,N)
               CASE( 3)
                  FZ(II,JJ,KK,N)   = FW(IW,N) + WW(II,JJ,KK)  *RHO_F(IW)*YY_F(IW,N)
               CASE(-3)
                  FZ(II,JJ,KK-1,N) = FW(IW,N) + WW(II,JJ,KK-1)*RHO_F(IW)*YY_F(IW,N)
            END SELECT
      
      END SELECT METHOD_OF_MASS_TRANSFER2
      
   ENDDO WALL_LOOP2

ENDDO SPECIES_LOOP

TUSED(3,NM)=TUSED(3,NM)+SECOND()-TNOW
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
