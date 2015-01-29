MODULE MASS

! Compute the mass equation differences

USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_POINTERS

IMPLICIT NONE
PRIVATE

CHARACTER(255), PARAMETER :: massid='$Id$'
CHARACTER(255), PARAMETER :: massrev='$Revision$'
CHARACTER(255), PARAMETER :: massdate='$Date$'

REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZP
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW,RHOP,DP,UP,VP,WP

PUBLIC MASS_FINITE_DIFFERENCES,DENSITY,GET_REV_mass,SCALAR_FACE_VALUE

CONTAINS


SUBROUTINE MASS_FINITE_DIFFERENCES(NM)

! Compute spatial differences for density equation

USE COMP_FUNCTIONS, ONLY: SECOND
USE GLOBAL_CONSTANTS, ONLY: N_TRACKED_SPECIES,PREDICTOR,EVACUATION_ONLY,SOLID_PHASE_ONLY,TUSED

INTEGER, INTENT(IN) :: NM
REAL(EB) :: TNOW,ZZZ(1:4)
INTEGER  :: I,J,K,N,IOR,IW,IIG,JJG,KKG,II,JJ,KK
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZP=>NULL()
REAL(EB), POINTER, DIMENSION(:,:,:) :: RHO_Z_P=>NULL()
TYPE(WALL_TYPE), POINTER :: WC=>NULL()

IF (EVACUATION_ONLY(NM) .OR. SOLID_PHASE_ONLY) RETURN

TNOW=SECOND()
CALL POINT_TO_MESH(NM)

IF (PREDICTOR) THEN
   UU => U
   VV => V
   WW => W
   RHOP => RHO
   ZZP => ZZ
ELSE
   UU => US
   VV => VS
   WW => WS
   RHOP => RHOS
   ZZP => ZZS
ENDIF

! Species face values

SPECIES_LOOP: DO N=1,N_TRACKED_SPECIES

   RHO_Z_P=>WORK1

   !$OMP PARALLEL PRIVATE(ZZZ)
   !$OMP DO SCHEDULE(STATIC)
   DO K=0,KBP1
      DO J=0,JBP1
         DO I=0,IBP1
            RHO_Z_P(I,J,K) = RHOP(I,J,K)*ZZP(I,J,K,N)
         ENDDO
      ENDDO
   ENDDO
   !$OMP END DO

   ! Compute scalar face values

   !$OMP DO SCHEDULE(STATIC)
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBM1
            ZZZ(1:4) = RHO_Z_P(I-1:I+2,J,K)
            FX(I,J,K,N) = SCALAR_FACE_VALUE(UU(I,J,K),ZZZ,FLUX_LIMITER)
         ENDDO
      ENDDO
   ENDDO
   !$OMP END DO NOWAIT

   !$OMP DO SCHEDULE(STATIC)
   DO K=1,KBAR
      DO J=1,JBM1
         DO I=1,IBAR
            ZZZ(1:4) = RHO_Z_P(I,J-1:J+2,K)
            FY(I,J,K,N) = SCALAR_FACE_VALUE(VV(I,J,K),ZZZ,FLUX_LIMITER)
         ENDDO
      ENDDO
   ENDDO
   !$OMP END DO NOWAIT

   !$OMP DO SCHEDULE(STATIC)
   DO K=1,KBM1
      DO J=1,JBAR
         DO I=1,IBAR
            ZZZ(1:4) = RHO_Z_P(I,J,K-1:K+2)
            FZ(I,J,K,N) = SCALAR_FACE_VALUE(WW(I,J,K),ZZZ,FLUX_LIMITER)
         ENDDO
      ENDDO
   ENDDO
   !$OMP END DO
   !$OMP END PARALLEL

   WALL_LOOP_2: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
      WC=>WALL(IW)
      IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE WALL_LOOP_2

      II  = WC%ONE_D%II 
      JJ  = WC%ONE_D%JJ
      KK  = WC%ONE_D%KK
      IIG = WC%ONE_D%IIG 
      JJG = WC%ONE_D%JJG
      KKG = WC%ONE_D%KKG
      IOR = WC%ONE_D%IOR

      SELECT CASE(IOR)
         CASE( 1)
            FX(IIG-1,JJG,KKG,N) = WC%RHO_F*WC%ZZ_F(N)
         CASE(-1)
            FX(IIG,JJG,KKG,N)   = WC%RHO_F*WC%ZZ_F(N)
         CASE( 2)
            FY(IIG,JJG-1,KKG,N) = WC%RHO_F*WC%ZZ_F(N)
         CASE(-2)
            FY(IIG,JJG,KKG,N)   = WC%RHO_F*WC%ZZ_F(N)
         CASE( 3)
            FZ(IIG,JJG,KKG-1,N) = WC%RHO_F*WC%ZZ_F(N)
         CASE(-3)
            FZ(IIG,JJG,KKG,N)   = WC%RHO_F*WC%ZZ_F(N)
      END SELECT

      ! Overwrite first off-wall advective flux if flow is away from the wall and if the face is not also a wall cell

      OFF_WALL_IF_2: IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY .AND. WC%BOUNDARY_TYPE/=OPEN_BOUNDARY) THEN

         OFF_WALL_SELECT_2: SELECT CASE(IOR)
            CASE( 1) OFF_WALL_SELECT_2
               !      ghost          FX/UU(II+1)
               ! ///   II   ///  II+1  |  II+2  | ...
               !                       ^ WALL_INDEX(II+1,+1)
               IF ((UU(II+1,JJ,KK)>0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II+1,JJ,KK),+1)>0)) THEN
                  ZZZ(1:3) = (/WC%RHO_F*WC%ZZ_F(N),RHO_Z_P(II+1:II+2,JJ,KK)/)
                  FX(II+1,JJ,KK,N) = SCALAR_FACE_VALUE(UU(II+1,JJ,KK),ZZZ,FLUX_LIMITER)
               ENDIF
            CASE(-1) OFF_WALL_SELECT_2
               !            FX/UU(II-2)     ghost
               ! ... |  II-2  |  II-1  ///   II   ///
               !              ^ WALL_INDEX(II-1,-1)
               IF ((UU(II-2,JJ,KK)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II-1,JJ,KK),-1)>0)) THEN
                  ZZZ(2:4) = (/RHO_Z_P(II-2:II-1,JJ,KK),WC%RHO_F*WC%ZZ_F(N)/)
                  FX(II-2,JJ,KK,N) = SCALAR_FACE_VALUE(UU(II-2,JJ,KK),ZZZ,FLUX_LIMITER)
               ENDIF
            CASE( 2) OFF_WALL_SELECT_2
               IF ((VV(II,JJ+1,KK)>0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ+1,KK),+2)>0)) THEN
                  ZZZ(1:3) = (/WC%RHO_F*WC%ZZ_F(N),RHO_Z_P(II,JJ+1:JJ+2,KK)/)
                  FY(II,JJ+1,KK,N) = SCALAR_FACE_VALUE(VV(II,JJ+1,KK),ZZZ,FLUX_LIMITER)
               ENDIF
            CASE(-2) OFF_WALL_SELECT_2
               IF ((VV(II,JJ-2,KK)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ-1,KK),-2)>0)) THEN
                  ZZZ(2:4) = (/RHO_Z_P(II,JJ-2:JJ-1,KK),WC%RHO_F*WC%ZZ_F(N)/)
                  FY(II,JJ-2,KK,N) = SCALAR_FACE_VALUE(VV(II,JJ-2,KK),ZZZ,FLUX_LIMITER)
               ENDIF
            CASE( 3) OFF_WALL_SELECT_2
               IF ((WW(II,JJ,KK+1)>0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ,KK+1),+3)>0)) THEN
                  ZZZ(1:3) = (/WC%RHO_F*WC%ZZ_F(N),RHO_Z_P(II,JJ,KK+1:KK+2)/)
                  FZ(II,JJ,KK+1,N) = SCALAR_FACE_VALUE(WW(II,JJ,KK+1),ZZZ,FLUX_LIMITER)
               ENDIF
            CASE(-3) OFF_WALL_SELECT_2
               IF ((WW(II,JJ,KK-2)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ,KK-1),-3)>0)) THEN
                  ZZZ(2:4) = (/RHO_Z_P(II,JJ,KK-2:KK-1),WC%RHO_F*WC%ZZ_F(N)/)
                  FZ(II,JJ,KK-2,N) = SCALAR_FACE_VALUE(WW(II,JJ,KK-2),ZZZ,FLUX_LIMITER)
               ENDIF
         END SELECT OFF_WALL_SELECT_2
         
      ENDIF OFF_WALL_IF_2

   ENDDO WALL_LOOP_2

ENDDO SPECIES_LOOP

TUSED(3,NM)=TUSED(3,NM)+SECOND()-TNOW

END SUBROUTINE MASS_FINITE_DIFFERENCES


SUBROUTINE DENSITY(T,NM)

! Update the species mass fractions and density

USE COMP_FUNCTIONS, ONLY: SECOND,SHUTDOWN
USE PHYSICAL_FUNCTIONS, ONLY : GET_SPECIFIC_GAS_CONSTANT,GET_SENSIBLE_ENTHALPY,GET_SPECIFIC_HEAT
USE GLOBAL_CONSTANTS, ONLY: N_TRACKED_SPECIES,TMPMAX,TMPMIN,EVACUATION_ONLY, &
                            PREDICTOR,N_ZONE,GAS_SPECIES,R0,SOLID_PHASE_ONLY,TUSED
USE MANUFACTURED_SOLUTIONS, ONLY: VD2D_MMS_Z_OF_RHO,VD2D_MMS_Z_SRC,UF_MMS,WF_MMS,VD2D_MMS_RHO_OF_Z,VD2D_MMS_Z_SRC
USE SOOT_ROUTINES, ONLY: SETTLING_VELOCITY
INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: T
REAL(EB) :: TNOW,ZZ_GET(1:N_TRACKED_SPECIES),RHS,UN,Q_Z,XHAT,ZHAT
INTEGER :: I,J,K,N,IW,IOR,IIG,JJG,KKG
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: RHO_ZZ__0=>NULL(),DEL_RHO_D_DEL_Z__0=>NULL()
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU=>NULL(),VV=>NULL(),WW=>NULL()
TYPE(WALL_TYPE), POINTER :: WC=>NULL()

IF (EVACUATION_ONLY(NM)) RETURN
IF (SOLID_PHASE_ONLY) RETURN

! If the RHS of the continuity equation does not yet satisfy the divergence constraint, return.
! This is typical of the case where an initial velocity field is specified by the user.

IF (PROJECTION .AND. ICYC<=1) RETURN
IF (PERIODIC_TEST==5) RETURN
IF (PERIODIC_TEST==8) RETURN

TNOW=SECOND()
CALL POINT_TO_MESH(NM)
UU=>WORK1
VV=>WORK2
WW=>WORK3
RHO_ZZ__0=>SCALAR_WORK1
DEL_RHO_D_DEL_Z__0=>SCALAR_WORK4

PREDICTOR_STEP: SELECT CASE (PREDICTOR)

CASE(.TRUE.) PREDICTOR_STEP

   IF (.NOT.CHANGE_TIME_STEP(NM)) THEN
      DEL_RHO_D_DEL_Z__0 = DEL_RHO_D_DEL_Z
      IF (ANY(SPECIES_MIXTURE%DEPOSITING) .AND. GRAVITATIONAL_DEPOSITION) CALL SETTLING_VELOCITY(NM)
   ENDIF

   ! Correct boundary velocity at wall cells

   UU=U
   VV=V
   WW=W

   WALL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
      WC=>WALL(IW)
      IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY) CYCLE WALL_LOOP

      IIG = WC%ONE_D%IIG 
      JJG = WC%ONE_D%JJG
      KKG = WC%ONE_D%KKG
      IOR = WC%ONE_D%IOR

      UN = UVW_SAVE(IW)
      SELECT CASE(IOR)
         CASE( 1)
            UU(IIG-1,JJG,KKG) = UN
         CASE(-1)
            UU(IIG,JJG,KKG)   = UN
         CASE( 2)
            VV(IIG,JJG-1,KKG) = UN
         CASE(-2)
            VV(IIG,JJG,KKG)   = UN
         CASE( 3)
            WW(IIG,JJG,KKG-1) = UN
         CASE(-3)
            WW(IIG,JJG,KKG)   = UN
      END SELECT
   ENDDO WALL_LOOP

   ! Predictor step for mass density

   DO N=1,N_TRACKED_SPECIES
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE

               RHS = - DEL_RHO_D_DEL_Z__0(I,J,K,N) &
                   + (FX(I,J,K,N)*UU(I,J,K)*R(I) - FX(I-1,J,K,N)*UU(I-1,J,K)*R(I-1))*RDX(I)*RRN(I) &
                   + (FY(I,J,K,N)*VV(I,J,K)      - FY(I,J-1,K,N)*VV(I,J-1,K)       )*RDY(J)        &
                   + (FZ(I,J,K,N)*WW(I,J,K)      - FZ(I,J,K-1,N)*WW(I,J,K-1)       )*RDZ(K)

               RHO_ZZ__0(I,J,K,N) = RHO(I,J,K)*ZZ(I,J,K,N)

               ZZS(I,J,K,N) = RHO_ZZ__0(I,J,K,N) - DT*RHS
            ENDDO
         ENDDO
      ENDDO
   ENDDO

   ! Manufactured solution

   IF (PERIODIC_TEST==7) THEN
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               ! divergence from EOS
               XHAT = XC(I) - UF_MMS*(T)
               ZHAT = ZC(K) - WF_MMS*(T)
               Q_Z = VD2D_MMS_Z_SRC(XHAT,ZHAT,T)
               ZZS(I,J,K,1) = ZZS(I,J,K,1) - DT*Q_Z
               ZZS(I,J,K,2) = ZZS(I,J,K,2) + DT*Q_Z
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   ! Get rho = sum(rho*Y_alpha)

   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
            RHOS(I,J,K) = SUM(ZZS(I,J,K,1:N_TRACKED_SPECIES))
         ENDDO
      ENDDO
   ENDDO

   ! Check mass density for positivity

   CALL CHECK_MASS_DENSITY

   ! Extract mass fraction from RHO * ZZ

   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
               ZZS(I,J,K,1:N_TRACKED_SPECIES) = ZZS(I,J,K,1:N_TRACKED_SPECIES)/RHOS(I,J,K)
         ENDDO
      ENDDO
   ENDDO

   ! Predict background pressure at next time step

   DO I=1,N_ZONE
      PBAR_S(:,I) = PBAR(:,I) + D_PBAR_DT(I)*DT
   ENDDO

   ! Compute molecular weight term RSUM=R0*SUM(Y_i/W_i)

   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
            ZZ_GET(1:N_TRACKED_SPECIES) = ZZS(I,J,K,1:N_TRACKED_SPECIES)            
            CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,RSUM(I,J,K))
         ENDDO
      ENDDO
   ENDDO

   ! Extract predicted temperature at next time step from Equation of State

   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE               
            TMP(I,J,K) = PBAR_S(K,PRESSURE_ZONE(I,J,K))/(RSUM(I,J,K)*RHOS(I,J,K))
         ENDDO
      ENDDO
   ENDDO

   TMP = MAX(TMPMIN,MIN(TMPMAX,TMP))

! The CORRECTOR step

CASE(.FALSE.) PREDICTOR_STEP

   ! Correct boundary velocity at wall cells

   UU=US
   VV=VS
   WW=WS

   WALL_LOOP_2: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
      WC=>WALL(IW)
      IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY) CYCLE WALL_LOOP_2

      IIG = WC%ONE_D%IIG 
      JJG = WC%ONE_D%JJG
      KKG = WC%ONE_D%KKG
      IOR = WC%ONE_D%IOR
      
      UN = UVW_SAVE(IW)
      SELECT CASE(IOR)
         CASE( 1)
            UU(IIG-1,JJG,KKG) = UN
         CASE(-1)
            UU(IIG,JJG,KKG)   = UN
         CASE( 2)
            VV(IIG,JJG-1,KKG) = UN
         CASE(-2)
            VV(IIG,JJG,KKG)   = UN
         CASE( 3)
            WW(IIG,JJG,KKG-1) = UN
         CASE(-3)
            WW(IIG,JJG,KKG)   = UN
      END SELECT
   ENDDO WALL_LOOP_2

   DEL_RHO_D_DEL_Z__0 = DEL_RHO_D_DEL_Z

   IF (ANY(SPECIES_MIXTURE%DEPOSITING) .AND. GRAVITATIONAL_DEPOSITION) CALL SETTLING_VELOCITY(NM)
   
   ! Compute species mass density at the next time step   
   
   DO N=1,N_TRACKED_SPECIES
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE

               RHS = - DEL_RHO_D_DEL_Z__0(I,J,K,N) &
                   + (FX(I,J,K,N)*UU(I,J,K)*R(I) - FX(I-1,J,K,N)*UU(I-1,J,K)*R(I-1))*RDX(I)*RRN(I) &
                   + (FY(I,J,K,N)*VV(I,J,K)      - FY(I,J-1,K,N)*VV(I,J-1,K)       )*RDY(J)        &
                   + (FZ(I,J,K,N)*WW(I,J,K)      - FZ(I,J,K-1,N)*WW(I,J,K-1)       )*RDZ(K)

               RHO_ZZ__0(I,J,K,N) = .5_EB*(RHO(I,J,K)*ZZ(I,J,K,N) + RHOS(I,J,K)*ZZS(I,J,K,N))

               ZZ(I,J,K,N) = .5_EB*(RHO(I,J,K)*ZZ(I,J,K,N) + RHOS(I,J,K)*ZZS(I,J,K,N)) - .5_EB*DT*RHS
            ENDDO
         ENDDO
      ENDDO
   ENDDO

   ! Manufactured solution

   IF (PERIODIC_TEST==7) THEN
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               ! divergence from EOS
               XHAT = XC(I) - UF_MMS*T
               ZHAT = ZC(K) - WF_MMS*T
               Q_Z = VD2D_MMS_Z_SRC(XHAT,ZHAT,T)
               ZZ(I,J,K,1) = ZZ(I,J,K,1) - .5_EB*DT*Q_Z
               ZZ(I,J,K,2) = ZZ(I,J,K,2) + .5_EB*DT*Q_Z
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   ! Get rho = sum(rho*Y_alpha)

   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
            RHO(I,J,K) = SUM(ZZ(I,J,K,1:N_TRACKED_SPECIES))
         ENDDO
      ENDDO
   ENDDO

   ! Check mass density for positivity

   CALL CHECK_MASS_DENSITY

   ! Extract Y_n from rho*Y_n

   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
            ZZ(I,J,K,1:N_TRACKED_SPECIES) = ZZ(I,J,K,1:N_TRACKED_SPECIES)/RHO(I,J,K)
         ENDDO
      ENDDO
   ENDDO

   ! Correct background pressure

   DO I=1,N_ZONE
      PBAR(:,I) = 0.5_EB*(PBAR(:,I) + PBAR_S(:,I) + D_PBAR_DT_S(I)*DT)
   ENDDO

   ! Compute molecular weight term RSUM=R0*SUM(Y_i/W_i)

   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
            ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(I,J,K,1:N_TRACKED_SPECIES)
            CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,RSUM(I,J,K))
         ENDDO
      ENDDO
   ENDDO

   ! Extract predicted temperature at next time step from Equation of State

   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
            TMP(I,J,K) = PBAR(K,PRESSURE_ZONE(I,J,K))/(RSUM(I,J,K)*RHO(I,J,K))
         ENDDO
      ENDDO
   ENDDO

   TMP = MAX(TMPMIN,MIN(TMPMAX,TMP))

END SELECT PREDICTOR_STEP

TUSED(3,NM)=TUSED(3,NM)+SECOND()-TNOW

END SUBROUTINE DENSITY


SUBROUTINE CHECK_MASS_DENSITY

! Redistribute mass from cells below or above the density cut-off limits
! Do not apply OpenMP to this routine

USE GLOBAL_CONSTANTS, ONLY : PREDICTOR,RHOMIN,RHOMAX
REAL(EB) :: MASS_N(-3:3),CONST,MASS_C,ZZ_CUT,RHO_CUT,VC(-3:3),SIGN_FACTOR,SUM_MASS_N,VC1(-3:3)
INTEGER  :: IC,I,J,K,N
REAL(EB), POINTER, DIMENSION(:,:,:) :: DELTA_RHO=>NULL(),DELTA_ZZ=>NULL(),RHOP=>NULL()
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZP=>NULL()
REAL(EB), PARAMETER :: ZZ_MIN=0._EB, ZZ_MAX=HUGE(1._EB)

DELTA_RHO => WORK1
DELTA_RHO =  0._EB

IF (PREDICTOR) THEN
   ZZP=>ZZS
   RHOP=>RHOS
ELSE
   ZZP=>ZZ
   RHOP=>RHO
ENDIF

! Correct density

DO K=1,KBAR
   DO J=1,JBAR
      VC1( 0)  = DY(J)  *DZ(K)
      VC1(-1)  = VC1( 0)
      VC1( 1)  = VC1( 0)
      VC1(-2)  = DY(J-1)*DZ(K)
      VC1( 2)  = DY(J+1)*DZ(K)
      VC1(-3)  = DY(J)  *DZ(K-1)
      VC1( 3)  = DY(J)  *DZ(K+1)
      DO I=1,IBAR

         IF (RHOP(I,J,K)>=RHOMIN .AND. RHOP(I,J,K)<=RHOMAX) CYCLE
         IC = CELL_INDEX(I,J,K)
         IF (SOLID(IC)) CYCLE
         IF (RHOP(I,J,K)<RHOMIN) THEN
            RHO_CUT = RHOMIN
            SIGN_FACTOR = 1._EB
         ELSE
            RHO_CUT = RHOMAX
            SIGN_FACTOR = -1._EB
         ENDIF
         MASS_N = 0._EB
         VC( 0)  = DX(I)  * VC1( 0)
         VC(-1)  = DX(I-1)* VC1(-1)
         VC( 1)  = DX(I+1)* VC1( 1)
         VC(-2)  = DX(I)  * VC1(-2)
         VC( 2)  = DX(I)  * VC1( 2)
         VC(-3)  = DX(I)  * VC1(-3)
         VC( 3)  = DX(I)  * VC1( 3)

         MASS_C = ABS(RHO_CUT-RHOP(I,J,K))*VC(0)
         IF (WALL_INDEX(IC,-1)==0) MASS_N(-1) = ABS(MIN(RHOMAX,MAX(RHOMIN,RHOP(I-1,J,K)))-RHO_CUT)*VC(-1)
         IF (WALL_INDEX(IC, 1)==0) MASS_N( 1) = ABS(MIN(RHOMAX,MAX(RHOMIN,RHOP(I+1,J,K)))-RHO_CUT)*VC( 1)
         IF (WALL_INDEX(IC,-2)==0) MASS_N(-2) = ABS(MIN(RHOMAX,MAX(RHOMIN,RHOP(I,J-1,K)))-RHO_CUT)*VC(-2)
         IF (WALL_INDEX(IC, 2)==0) MASS_N( 2) = ABS(MIN(RHOMAX,MAX(RHOMIN,RHOP(I,J+1,K)))-RHO_CUT)*VC( 2)
         IF (WALL_INDEX(IC,-3)==0) MASS_N(-3) = ABS(MIN(RHOMAX,MAX(RHOMIN,RHOP(I,J,K-1)))-RHO_CUT)*VC(-3)
         IF (WALL_INDEX(IC, 3)==0) MASS_N( 3) = ABS(MIN(RHOMAX,MAX(RHOMIN,RHOP(I,J,K+1)))-RHO_CUT)*VC( 3)
         SUM_MASS_N = SUM(MASS_N)
         IF (SUM_MASS_N<=TWO_EPSILON_EB) CYCLE
         CONST = SIGN_FACTOR*MIN(1._EB,MASS_C/SUM_MASS_N)
         DELTA_RHO(I,J,K)   = DELTA_RHO(I,J,K)   + CONST*SUM_MASS_N/VC( 0)
         DELTA_RHO(I-1,J,K) = DELTA_RHO(I-1,J,K) - CONST*MASS_N(-1)/VC(-1)
         DELTA_RHO(I+1,J,K) = DELTA_RHO(I+1,J,K) - CONST*MASS_N( 1)/VC( 1)
         DELTA_RHO(I,J-1,K) = DELTA_RHO(I,J-1,K) - CONST*MASS_N(-2)/VC(-2)
         DELTA_RHO(I,J+1,K) = DELTA_RHO(I,J+1,K) - CONST*MASS_N( 2)/VC( 2)
         DELTA_RHO(I,J,K-1) = DELTA_RHO(I,J,K-1) - CONST*MASS_N(-3)/VC(-3)
         DELTA_RHO(I,J,K+1) = DELTA_RHO(I,J,K+1) - CONST*MASS_N( 3)/VC( 3)
      ENDDO
   ENDDO
ENDDO

RHOP(1:IBAR,1:JBAR,1:KBAR) = MIN(RHOMAX,MAX(RHOMIN,RHOP(1:IBAR,1:JBAR,1:KBAR)+DELTA_RHO(1:IBAR,1:JBAR,1:KBAR)))

! Correct species mass density

SPECIES_LOOP: DO N=1,N_TRACKED_SPECIES

DELTA_ZZ => WORK2
DELTA_ZZ = 0._EB

DO K=1,KBAR
   DO J=1,JBAR
      VC1( 0)  = DY(J)  *DZ(K)
      VC1(-1)  = VC1( 0)
      VC1( 1)  = VC1( 0)
      VC1(-2)  = DY(J-1)*DZ(K)
      VC1( 2)  = DY(J+1)*DZ(K)
      VC1(-3)  = DY(J)  *DZ(K-1)
      VC1( 3)  = DY(J)  *DZ(K+1)
      DO I=1,IBAR

         IF (ZZP(I,J,K,N)>=ZZ_MIN .AND. ZZP(I,J,K,N)<=ZZ_MAX) CYCLE
         IC = CELL_INDEX(I,J,K)
         IF (SOLID(IC)) CYCLE
         IF (ZZP(I,J,K,N)<ZZ_MIN) THEN
            ZZ_CUT = ZZ_MIN
            SIGN_FACTOR = 1._EB
         ELSE
            ZZ_CUT = ZZ_MAX
            SIGN_FACTOR = -1._EB
         ENDIF
         MASS_N = 0._EB
         VC( 0)  = DX(I)  * VC1( 0)
         VC(-1)  = DX(I-1)* VC1(-1)
         VC( 1)  = DX(I+1)* VC1( 1)
         VC(-2)  = DX(I)  * VC1(-2)
         VC( 2)  = DX(I)  * VC1( 2)
         VC(-3)  = DX(I)  * VC1(-3)
         VC( 3)  = DX(I)  * VC1( 3)

         MASS_C = ABS(ZZ_CUT-ZZP(I,J,K,N))*VC(0)
         IF (WALL_INDEX(IC,-1)==0) MASS_N(-1) = ABS(MIN(ZZ_MAX,MAX(ZZ_MIN,ZZP(I-1,J,K,N)))-ZZ_CUT)*VC(-1)
         IF (WALL_INDEX(IC, 1)==0) MASS_N( 1) = ABS(MIN(ZZ_MAX,MAX(ZZ_MIN,ZZP(I+1,J,K,N)))-ZZ_CUT)*VC( 1)
         IF (WALL_INDEX(IC,-2)==0) MASS_N(-2) = ABS(MIN(ZZ_MAX,MAX(ZZ_MIN,ZZP(I,J-1,K,N)))-ZZ_CUT)*VC(-2)
         IF (WALL_INDEX(IC, 2)==0) MASS_N( 2) = ABS(MIN(ZZ_MAX,MAX(ZZ_MIN,ZZP(I,J+1,K,N)))-ZZ_CUT)*VC( 2)
         IF (WALL_INDEX(IC,-3)==0) MASS_N(-3) = ABS(MIN(ZZ_MAX,MAX(ZZ_MIN,ZZP(I,J,K-1,N)))-ZZ_CUT)*VC(-3)
         IF (WALL_INDEX(IC, 3)==0) MASS_N( 3) = ABS(MIN(ZZ_MAX,MAX(ZZ_MIN,ZZP(I,J,K+1,N)))-ZZ_CUT)*VC( 3)
         SUM_MASS_N = SUM(MASS_N)
         IF (SUM_MASS_N<=TWO_EPSILON_EB) CYCLE
         CONST = SIGN_FACTOR*MIN(1._EB,MASS_C/SUM_MASS_N)
         DELTA_ZZ(I,J,K)   = DELTA_ZZ(I,J,K)   + CONST*SUM_MASS_N/VC( 0)
         DELTA_ZZ(I-1,J,K) = DELTA_ZZ(I-1,J,K) - CONST*MASS_N(-1)/VC(-1)
         DELTA_ZZ(I+1,J,K) = DELTA_ZZ(I+1,J,K) - CONST*MASS_N( 1)/VC( 1)
         DELTA_ZZ(I,J-1,K) = DELTA_ZZ(I,J-1,K) - CONST*MASS_N(-2)/VC(-2)
         DELTA_ZZ(I,J+1,K) = DELTA_ZZ(I,J+1,K) - CONST*MASS_N( 2)/VC( 2)
         DELTA_ZZ(I,J,K-1) = DELTA_ZZ(I,J,K-1) - CONST*MASS_N(-3)/VC(-3)
         DELTA_ZZ(I,J,K+1) = DELTA_ZZ(I,J,K+1) - CONST*MASS_N( 3)/VC( 3)
      ENDDO
   ENDDO
ENDDO

ZZP(1:IBAR,1:JBAR,1:KBAR,N) = MIN(ZZ_MAX,MAX(ZZ_MIN,ZZP(1:IBAR,1:JBAR,1:KBAR,N)+DELTA_ZZ(1:IBAR,1:JBAR,1:KBAR)))

ENDDO SPECIES_LOOP

! Absorb error in most abundant species

DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
         N=MAXLOC(ZZP(I,J,K,:),1)
         ZZP(I,J,K,N) = RHOP(I,J,K) - ( SUM(ZZP(I,J,K,:)) - ZZP(I,J,K,N) )
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE CHECK_MASS_DENSITY


REAL(EB) FUNCTION SCALAR_FACE_VALUE(A,U,LIMITER)

REAL(EB), INTENT(IN) :: A,U(4)
INTEGER, INTENT(IN) :: LIMITER
REAL(EB) :: R,B,DU_UP,DU_LOC,V(5)

! This function computes the scalar value on a face.
! The scalar is denoted U, and the velocity is denoted A.
! The divergence (computed elsewhere) uses a central difference across
! the cell subject to a flux LIMITER.  The flux LIMITER choices are:
!
! CENTRAL_LIMITER  = 0
! GODUNOV_LIMITER  = 1
! SUPERBEE_LIMITER = 2
! MINMOD_LIMITER   = 3
! CHARM_LIMITER    = 4
! MP5_LIMITER      = 5
!
!                    location of face
!
!                            f
!    |     o     |     o     |     o     |     o     |
!                            A
!         U(1)        U(2)        U(3)        U(4)

WIND_DIRECTION_IF: IF (A>0._EB) THEN

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
         IF (ABS(DU_LOC)>TWO_EPSILON_EB) R = DU_UP/DU_LOC
         B = MAX(0._EB,MIN(2._EB*R,1._EB),MIN(R,2._EB))
         SCALAR_FACE_VALUE = U(2) + 0.5_EB*B*DU_LOC
      CASE(3) ! MINMOD
         IF (ABS(DU_LOC)>TWO_EPSILON_EB) R = DU_UP/DU_LOC
         B = MAX(0._EB,MIN(1._EB,R))
         SCALAR_FACE_VALUE = U(2) + 0.5_EB*B*DU_LOC
      CASE(4) ! CHARM
         IF (ABS(DU_UP)>TWO_EPSILON_EB) R = DU_LOC/DU_UP
         IF (R>0._EB) B = R*(3._EB*R+1._EB)/((R+1._EB)**2)
         SCALAR_FACE_VALUE = U(2) + 0.5_EB*B*DU_UP
      CASE(5) ! MP5, Suresh and Huynh (1997)
         V = (/2._EB*U(1)-U(2),U(1:4)/)
         SCALAR_FACE_VALUE = MP5(V)
   END SELECT

ELSE WIND_DIRECTION_IF

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
         IF (ABS(DU_LOC)>TWO_EPSILON_EB) R = DU_UP/DU_LOC
         B = MAX(0._EB,MIN(2._EB*R,1._EB),MIN(R,2._EB))
         SCALAR_FACE_VALUE = U(3) - 0.5_EB*B*DU_LOC
      CASE(3) ! MINMOD
         IF (ABS(DU_LOC)>TWO_EPSILON_EB) R = DU_UP/DU_LOC
         B = MAX(0._EB,MIN(1._EB,R))
         SCALAR_FACE_VALUE = U(3) - 0.5_EB*B*DU_LOC
      CASE(4) ! CHARM
         IF (ABS(DU_UP)>TWO_EPSILON_EB) R = DU_LOC/DU_UP
         IF (R>0._EB) B = R*(3._EB*R+1._EB)/((R+1._EB)**2)
         SCALAR_FACE_VALUE = U(3) - 0.5_EB*B*DU_UP
      CASE(5) ! MP5, Suresh and Huynh (1997)
         V = (/2._EB*U(4)-U(3),U(4),U(3),U(2),U(1)/)
         SCALAR_FACE_VALUE = MP5(V)
    END SELECT

ENDIF WIND_DIRECTION_IF

END FUNCTION SCALAR_FACE_VALUE


REAL(EB) FUNCTION MP5(V)
USE MATH_FUNCTIONS, ONLY: MINMOD2,MINMOD4
REAL(EB), INTENT(IN) :: V(-2:2)
REAL(EB), PARAMETER :: B1 = 0.016666666666667_EB, B2 = 1.333333333333_EB, ALPHA=4._EB, EPSM=1.E-10_EB
REAL(EB) :: VOR,VMP,DJM1,DJ,DJP1,DM4JPH,DM4JMH,VUL,VAV,VMD,VLC,VMIN,VMAX

! Monotonicity preserving 5th-order scheme (MP5) of Suresh and Huynh, JCP 136, 83-99 (1997)

VOR = B1*(2._EB*V(-2)-13._EB*V(-1)+47._EB*V(0)+27._EB*V(1)-3._EB*V(2))
VMP = V(0) + MINMOD2(V(1)-V(0),ALPHA*(V(0)-V(-1)))
IF ((VOR-V(0))*(VOR-VMP)<EPSM) THEN
   MP5=VOR
ELSE
   DJM1 = V(-2)-2._EB*V(-1)+V(0)
   DJ   = V(-1)-2._EB*V(0) +V(1)
   DJP1 = V(0) -2._EB*V(1) +V(2)
   DM4JPH = MINMOD4(4._EB*DJ-DJP1,4._EB*DJP1-DJ,DJ,DJP1)
   DM4JMH = MINMOD4(4._EB*DJ-DJM1,4._EB*DJM1-DJ,DJ,DJM1)
   VUL = V(0) + ALPHA*(V(0)-V(-1))
   VAV = 0.5_EB*(V(0)+V(1))
   VMD = VAV - 0.5_EB*DM4JPH
   VLC = V(0) + 0.5_EB*(V(0)-V(-1)) + B2*DM4JMH
   VMIN = MAX(MIN(V(0),V(1),VMD),MIN(V(0),VUL,VLC))
   VMAX = MIN(MAX(V(0),V(1),VMD),MAX(V(0),VUL,VLC))
   MP5 = VOR + MINMOD2(VMIN-VOR,VMAX-VOR)
ENDIF

END FUNCTION MP5



SUBROUTINE GET_REV_mass(MODULE_REV,MODULE_DATE)
INTEGER,INTENT(INOUT) :: MODULE_REV
CHARACTER(255),INTENT(INOUT) :: MODULE_DATE
WRITE(MODULE_DATE,'(A)') massrev(INDEX(massrev,':')+2:LEN_TRIM(massrev)-2)
READ (MODULE_DATE,'(I5)') MODULE_REV
WRITE(MODULE_DATE,'(A)') massdate
END SUBROUTINE GET_REV_mass

END MODULE MASS
