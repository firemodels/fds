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

PUBLIC MASS_FINITE_DIFFERENCES,DENSITY,GET_REV_mass,SCALAR_FACE_VALUE,MASS_FINITE_DIFFERENCES_2,DENSITY_2

CONTAINS

SUBROUTINE MASS_FINITE_DIFFERENCES(NM)

! Compute spatial differences for density equation

USE COMP_FUNCTIONS, ONLY: SECOND
USE GLOBAL_CONSTANTS, ONLY: N_TRACKED_SPECIES,PREDICTOR,CORRECTOR,EVACUATION_ONLY,SOLID_PHASE_ONLY,TUSED
USE SOOT_ROUTINES, ONLY: SETTLING_VELOCITY

INTEGER, INTENT(IN) :: NM
REAL(EB) :: TNOW,UN,ZZZ(1:4)
INTEGER  :: I,J,K,N,IOR,IW,IIG,JJG,KKG,II,JJ,KK,WALL_BOUNDARY_TYPE
LOGICAL :: ALLOW_TRANSPORT_CYCLE
REAL(EB), POINTER, DIMENSION(:,:,:) :: RHO_Z_P=>NULL()
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZP=>NULL()
TYPE(WALL_TYPE), POINTER :: WC=>NULL()

IF (EVACUATION_ONLY(NM) .OR. SOLID_PHASE_ONLY) RETURN

TNOW=SECOND()
CALL POINT_TO_MESH(NM)

IF (PREDICTOR) THEN
   UU => U
   VV => V
   WW => W
   UP => US
   VP => VS
   WP => WS
   RHOP => RHO
   IF (N_TRACKED_SPECIES > 0) ZZP => ZZ
ELSE
   UU => US
   VV => VS
   WW => WS
   UP => U
   VP => V
   WP => W
   RHOP => RHOS
   IF (N_TRACKED_SPECIES > 0) ZZP => ZZS
ENDIF

! Compute scalar face values

! Note: FX,FY,FZ are computed in divg based on old values of the wind direction.  If the
! wind direction has not changed, there is no need to recompute the face value, which is
! a relatively expensive operation.

IF (ENTHALPY_TRANSPORT .AND. .NOT.CONSTANT_SPECIFIC_HEAT_RATIO) THEN
   ALLOW_TRANSPORT_CYCLE=.TRUE.  ! fluxes may have been computed already in divg
ELSE
   ALLOW_TRANSPORT_CYCLE=.FALSE. ! because fluxes not yet computed in divg
ENDIF

! Special case for RESTARTs

IF (ICYC-ICYC_RESTART < 2) ALLOW_TRANSPORT_CYCLE = .FALSE.

! Compute mass fluxes at cell faces, FX, FY, FZ

!$OMP PARALLEL PRIVATE(ZZZ)
!$OMP DO SCHEDULE(STATIC)
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBM1
         IF (ALLOW_TRANSPORT_CYCLE .AND. SIGN(1._EB,UU(I,J,K))==SIGN(1._EB,UP(I,J,K))) CYCLE
         ZZZ(1:4) = RHOP(I-1:I+2,J,K)
         FX(I,J,K,0) = SCALAR_FACE_VALUE(UU(I,J,K),ZZZ,FLUX_LIMITER)
      ENDDO
   ENDDO
ENDDO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO K=1,KBAR
   DO J=1,JBM1
      DO I=1,IBAR
         IF (ALLOW_TRANSPORT_CYCLE .AND. SIGN(1._EB,VV(I,J,K))==SIGN(1._EB,VP(I,J,K))) CYCLE
         ZZZ(1:4) = RHOP(I,J-1:J+2,K)
         FY(I,J,K,0) = SCALAR_FACE_VALUE(VV(I,J,K),ZZZ,FLUX_LIMITER)
      ENDDO
   ENDDO
ENDDO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO K=1,KBM1
   DO J=1,JBAR
      DO I=1,IBAR
         IF (ALLOW_TRANSPORT_CYCLE .AND. SIGN(1._EB,WW(I,J,K))==SIGN(1._EB,WP(I,J,K))) CYCLE
         ZZZ(1:4) = RHOP(I,J,K-1:K+2)
         FZ(I,J,K,0) = SCALAR_FACE_VALUE(WW(I,J,K),ZZZ,FLUX_LIMITER)
      ENDDO
   ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

FRHO = 0._EB ! Zero out RHS of continuity equation

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

   ! Overwrite first off-wall advective flux if flow is away from the wall and if the face is not also a wall cell

   OFF_WALL_IF: IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY .AND. WC%BOUNDARY_TYPE/=OPEN_BOUNDARY) THEN

      OFF_WALL_SELECT: SELECT CASE(IOR)
         CASE( 1) OFF_WALL_SELECT
            !      ghost          FX/UU(II+1)
            ! ///   II   ///  II+1  |  II+2  | ...
            !                       ^ WALL_INDEX(II+1,+1)
            IF ((UU(II+1,JJ,KK)>0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II+1,JJ,KK),+1)>0)) THEN
               ZZZ(1:3) = (/WC%RHO_F,RHOP(II+1:II+2,JJ,KK)/)
               FX(II+1,JJ,KK,0) = SCALAR_FACE_VALUE(UU(II+1,JJ,KK),ZZZ,FLUX_LIMITER)
            ENDIF
         CASE(-1) OFF_WALL_SELECT
            !            FX/UU(II-2)     ghost
            ! ... |  II-2  |  II-1  ///   II   ///
            !              ^ WALL_INDEX(II-1,-1)
            IF ((UU(II-2,JJ,KK)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II-1,JJ,KK),-1)>0)) THEN
               ZZZ(2:4) = (/RHOP(II-2:II-1,JJ,KK),WC%RHO_F/)
               FX(II-2,JJ,KK,0) = SCALAR_FACE_VALUE(UU(II-2,JJ,KK),ZZZ,FLUX_LIMITER)
            ENDIF
         CASE( 2) OFF_WALL_SELECT
            IF ((VV(II,JJ+1,KK)>0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ+1,KK),+2)>0)) THEN
               ZZZ(1:3) = (/WC%RHO_F,RHOP(II,JJ+1:JJ+2,KK)/)
               FY(II,JJ+1,KK,0) = SCALAR_FACE_VALUE(VV(II,JJ+1,KK),ZZZ,FLUX_LIMITER)
            ENDIF
         CASE(-2) OFF_WALL_SELECT
            IF ((VV(II,JJ-2,KK)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ-1,KK),-2)>0)) THEN
               ZZZ(2:4) = (/RHOP(II,JJ-2:JJ-1,KK),WC%RHO_F/)
               FY(II,JJ-2,KK,0) = SCALAR_FACE_VALUE(VV(II,JJ-2,KK),ZZZ,FLUX_LIMITER)
            ENDIF
         CASE( 3) OFF_WALL_SELECT
            IF ((WW(II,JJ,KK+1)>0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ,KK+1),+3)>0)) THEN
               ZZZ(1:3) = (/WC%RHO_F,RHOP(II,JJ,KK+1:KK+2)/)
               FZ(II,JJ,KK+1,0) = SCALAR_FACE_VALUE(WW(II,JJ,KK+1),ZZZ,FLUX_LIMITER)
            ENDIF
         CASE(-3) OFF_WALL_SELECT
            IF ((WW(II,JJ,KK-2)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ,KK-1),-3)>0)) THEN
               ZZZ(2:4) = (/RHOP(II,JJ,KK-2:KK-1),WC%RHO_F/)
               FZ(II,JJ,KK-2,0) = SCALAR_FACE_VALUE(WW(II,JJ,KK-2),ZZZ,FLUX_LIMITER)
            ENDIF
      END SELECT OFF_WALL_SELECT
   
   ENDIF OFF_WALL_IF

   WALL_BOUNDARY_TYPE=WC%BOUNDARY_TYPE
   IF (WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .AND. EMBEDDED_MESH) WALL_BOUNDARY_TYPE=0

   BOUNDARY_SELECT: SELECT CASE(WALL_BOUNDARY_TYPE)
      CASE DEFAULT
         SELECT CASE(IOR)
            CASE( 1)
               FX(IIG-1,JJG,KKG,0) = WC%RHO_F
            CASE(-1)
               FX(IIG,JJG,KKG,0)   = WC%RHO_F
            CASE( 2)
               FY(IIG,JJG-1,KKG,0) = WC%RHO_F
            CASE(-2)
               FY(IIG,JJG,KKG,0)   = WC%RHO_F
            CASE( 3)
               FZ(IIG,JJG,KKG-1,0) = WC%RHO_F
            CASE(-3)
               FZ(IIG,JJG,KKG,0)   = WC%RHO_F
         END SELECT
      CASE(INTERPOLATED_BOUNDARY)
         UN = UVW_SAVE(IW)
         SELECT CASE(IOR)
            CASE( 1)
               FX(IIG-1,JJG,KKG,0) = 0._EB
               FRHO(IIG,JJG,KKG)   = FRHO(IIG,JJG,KKG) + WC%RHO_F*UN*RDX(IIG)*R(IIG-1)*RRN(IIG)
            CASE(-1)
               FX(IIG,JJG,KKG,0)   = 0._EB
               FRHO(IIG,JJG,KKG)   = FRHO(IIG,JJG,KKG) - WC%RHO_F*UN*RDX(IIG)*R(IIG)*RRN(IIG)
            CASE( 2)
               FY(IIG,JJG-1,KKG,0) = 0._EB
               FRHO(IIG,JJG,KKG)   = FRHO(IIG,JJG,KKG) + WC%RHO_F*UN*RDY(JJG)
            CASE(-2)
               FY(IIG,JJG,KKG,0)   = 0._EB
               FRHO(IIG,JJG,KKG)   = FRHO(IIG,JJG,KKG) - WC%RHO_F*UN*RDY(JJG)
            CASE( 3)
               FZ(IIG,JJG,KKG-1,0) = 0._EB
               FRHO(IIG,JJG,KKG)   = FRHO(IIG,JJG,KKG) + WC%RHO_F*UN*RDZ(KKG)
            CASE(-3)
               FZ(IIG,JJG,KKG,0)   = 0._EB
               FRHO(IIG,JJG,KKG)   = FRHO(IIG,JJG,KKG) - WC%RHO_F*UN*RDZ(KKG)
         END SELECT
   END SELECT BOUNDARY_SELECT

ENDDO WALL_LOOP

!$OMP PARALLEL DO SCHEDULE(STATIC)
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
         FRHO(I,J,K) = -FRHO(I,J,K)                                                                &
                     + (FX(I,J,K,0)*UU(I,J,K)*R(I)-FX(I-1,J,K,0)*UU(I-1,J,K)*R(I-1))*RDX(I)*RRN(I) &
                     + (FY(I,J,K,0)*VV(I,J,K)     -FY(I,J-1,K,0)*VV(I,J-1,K)       )*RDY(J)        &
                     + (FZ(I,J,K,0)*WW(I,J,K)     -FZ(I,J,K-1,0)*WW(I,J,K-1)       )*RDZ(K)
      ENDDO
   ENDDO
ENDDO
!$OMP END PARALLEL DO

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
            IF (ALLOW_TRANSPORT_CYCLE .AND. SIGN(1._EB,UU(I,J,K))==SIGN(1._EB,UP(I,J,K))) CYCLE
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
            IF (ALLOW_TRANSPORT_CYCLE .AND. SIGN(1._EB,VV(I,J,K))==SIGN(1._EB,VP(I,J,K))) CYCLE
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
            IF (ALLOW_TRANSPORT_CYCLE .AND. SIGN(1._EB,WW(I,J,K))==SIGN(1._EB,WP(I,J,K))) CYCLE
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

      WALL_BOUNDARY_TYPE=WC%BOUNDARY_TYPE
      IF (WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .AND. EMBEDDED_MESH) WALL_BOUNDARY_TYPE=0

      BOUNDARY_SELECT_2: SELECT CASE(WALL_BOUNDARY_TYPE)
         CASE DEFAULT
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
         CASE(INTERPOLATED_BOUNDARY,SOLID_BOUNDARY)
            SELECT CASE(WC%BOUNDARY_TYPE)
               CASE(SOLID_BOUNDARY)
                  IF (PREDICTOR) UN = -SIGN(1._EB,REAL(IOR,EB))*WC%ONE_D%UW
                  IF (CORRECTOR) UN = -SIGN(1._EB,REAL(IOR,EB))*WC%ONE_D%UWS
               CASE(INTERPOLATED_BOUNDARY)
                  UN = UVW_SAVE(IW)
            END SELECT
            SELECT CASE(IOR)
               CASE( 1)
                  FX(IIG-1,JJG,KKG,N) = 0._EB
                  DEL_RHO_D_DEL_Z(IIG,JJG,KKG,N) = DEL_RHO_D_DEL_Z(IIG,JJG,KKG,N) &
                                                 + WC%RHO_F*WC%ZZ_F(N)*UN*RDX(IIG)*R(IIG-1)*RRN(IIG)
               CASE(-1)
                  FX(IIG,JJG,KKG,N)   = 0._EB
                  DEL_RHO_D_DEL_Z(IIG,JJG,KKG,N) = DEL_RHO_D_DEL_Z(IIG,JJG,KKG,N) &
                                                 - WC%RHO_F*WC%ZZ_F(N)*UN*RDX(IIG)*R(IIG)*RRN(IIG)
               CASE( 2)
                  FY(IIG,JJG-1,KKG,N) = 0._EB
                  DEL_RHO_D_DEL_Z(IIG,JJG,KKG,N) = DEL_RHO_D_DEL_Z(IIG,JJG,KKG,N) + WC%RHO_F*WC%ZZ_F(N)*UN*RDY(JJG)
               CASE(-2)
                  FY(IIG,JJG,KKG,N)   = 0._EB
                  DEL_RHO_D_DEL_Z(IIG,JJG,KKG,N) = DEL_RHO_D_DEL_Z(IIG,JJG,KKG,N) - WC%RHO_F*WC%ZZ_F(N)*UN*RDY(JJG)
               CASE( 3)
                  FZ(IIG,JJG,KKG-1,N) = 0._EB
                  DEL_RHO_D_DEL_Z(IIG,JJG,KKG,N) = DEL_RHO_D_DEL_Z(IIG,JJG,KKG,N) + WC%RHO_F*WC%ZZ_F(N)*UN*RDZ(KKG)
               CASE(-3)
                  FZ(IIG,JJG,KKG,N)   = 0._EB
                  DEL_RHO_D_DEL_Z(IIG,JJG,KKG,N) = DEL_RHO_D_DEL_Z(IIG,JJG,KKG,N) - WC%RHO_F*WC%ZZ_F(N)*UN*RDZ(KKG)
            END SELECT
      END SELECT BOUNDARY_SELECT_2

   ENDDO WALL_LOOP_2

   !$OMP PARALLEL DO SCHEDULE(STATIC)
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
            DEL_RHO_D_DEL_Z(I,J,K,N) = -DEL_RHO_D_DEL_Z(I,J,K,N)                                                   &
                                     + (FX(I,J,K,N)*UU(I,J,K)*R(I)-FX(I-1,J,K,N)*UU(I-1,J,K)*R(I-1))*RDX(I)*RRN(I) &
                                     + (FY(I,J,K,N)*VV(I,J,K)     -FY(I,J-1,K,N)*VV(I,J-1,K)       )*RDY(J)        &
                                     + (FZ(I,J,K,N)*WW(I,J,K)     -FZ(I,J,K-1,N)*WW(I,J,K-1)       )*RDZ(K)
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO

   IF (SPECIES_MIXTURE(N)%DEPOSITING .AND. GRAVITATIONAL_DEPOSITION) CALL SETTLING_VELOCITY(NM,N)

ENDDO SPECIES_LOOP

TUSED(3,NM)=TUSED(3,NM)+SECOND()-TNOW

END SUBROUTINE MASS_FINITE_DIFFERENCES


SUBROUTINE DENSITY(NM)

! Update the density and species mass fractions

USE COMP_FUNCTIONS, ONLY: SECOND
USE PHYSICAL_FUNCTIONS, ONLY : GET_SPECIFIC_GAS_CONSTANT,GET_SENSIBLE_ENTHALPY,GET_SPECIFIC_HEAT,GET_SENSIBLE_ENTHALPY_DIFF
USE GLOBAL_CONSTANTS, ONLY: N_TRACKED_SPECIES,TMPMAX,TMPMIN,EVACUATION_ONLY, &
                            PREDICTOR,CHANGE_TIME_STEP,N_ZONE, &
                            GAS_SPECIES, R0,SOLID_PHASE_ONLY,TUSED, &
                            CLIP_MASS_FRACTION
USE MANUFACTURED_SOLUTIONS, ONLY: VD2D_MMS_Z_OF_RHO
REAL(EB) :: DTRATIO,OMDTRATIO,TNOW,ZZ_GET(0:N_TRACKED_SPECIES)
INTEGER  :: I,J,K,N
INTEGER, INTENT(IN) :: NM

IF (EVACUATION_ONLY(NM)) RETURN
IF (SOLID_PHASE_ONLY) RETURN

! If the RHS of the continuity equation does not yet satisfy the divergence constraint, return.
! This is typical of the case where an initial velocity field is specified by the user.

IF (PROJECTION .AND. ICYC<=1) RETURN
IF (PERIODIC_TEST==5) RETURN
IF (PERIODIC_TEST==8) RETURN

TNOW=SECOND()
CALL POINT_TO_MESH(NM)

PREDICTOR_STEP: SELECT CASE (PREDICTOR)

CASE(.TRUE.) PREDICTOR_STEP


   IF (.NOT.CHANGE_TIME_STEP(NM)) THEN

   ! NOTE: This IF statement is required because the source terms for species are zeroed out at
   !       the beginning of DIVERGENCE_PART_1, but the array also stores the divergence of the advective
   !       flux which is computed once in MASS_FINITE_DIFFERNENCES above, outside the CHANGE_TIME_STEP loop.
   !       DIVERGENCE_PART_1 is inside the loop.  The source terms are then applied to the next substep in
   !       MASS_FINITE_DIFFERENCES.

      DO N=1,N_TRACKED_SPECIES
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
                  ZZS(I,J,K,N) = RHO(I,J,K)*ZZ(I,J,K,N) - DT*DEL_RHO_D_DEL_Z(I,J,K,N)
               ENDDO
            ENDDO
         ENDDO
      ENDDO

   ELSE

      DTRATIO   = DT/DT_PREV
      OMDTRATIO = 1._EB - DTRATIO
      DO N=1,N_TRACKED_SPECIES
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
                  ZZS(I,J,K,N) = OMDTRATIO*RHO(I,J,K)*ZZ(I,J,K,N) + DTRATIO*RHOS(I,J,K)*ZZS(I,J,K,N)
               ENDDO
           ENDDO
         ENDDO
      ENDDO

   ENDIF

   ! Predict the density at the next time step (RHOS or RHO^*)

   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
            RHOS(I,J,K) = RHO(I,J,K)-DT*FRHO(I,J,K)
         ENDDO
      ENDDO
   ENDDO

   ! Correct densities above or below clip limits

   IF (.NOT.CLIP_MASS_FRACTION) CALL CHECK_DENSITY

   ! Extract mass fraction from RHO * ZZ

   DO N=1,N_TRACKED_SPECIES
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
               ZZS(I,J,K,N) = ZZS(I,J,K,N)/RHOS(I,J,K)
            ENDDO
         ENDDO
      ENDDO
   ENDDO

   ! Predict background pressure at next time step

   DO I=1,N_ZONE
      PBAR_S(:,I) = PBAR(:,I) + D_PBAR_DT(I)*DT
   ENDDO

   ! Manufactured solution

   IF (PERIODIC_TEST==7) THEN
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               ZZS(I,J,K,1) = VD2D_MMS_Z_OF_RHO(RHOS(I,J,K))
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   ! Correct mass fractions above or below clip limits

   IF (CLIP_MASS_FRACTION .AND. N_TRACKED_SPECIES>0) THEN
      ZZS(1:IBAR,1:JBAR,1:KBAR,1:N_TRACKED_SPECIES) = MAX(0._EB,MIN(1._EB,ZZS(1:IBAR,1:JBAR,1:KBAR,1:N_TRACKED_SPECIES)))
   ELSEIF (N_TRACKED_SPECIES>0) THEN
      IF (USE_NEW_CHECK_MASS_FRACTION) THEN
         CALL NEW_CHECK_MASS_FRACTION
      ELSE
         CALL CHECK_MASS_FRACTION
      ENDIF
   ENDIF

   ! Compute molecular weight term RSUM=R0*SUM(Y_i/W_i)

   IF (N_TRACKED_SPECIES>0) THEN
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
               ZZ_GET(1:N_TRACKED_SPECIES) = ZZS(I,J,K,1:N_TRACKED_SPECIES)
               CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,RSUM(I,J,K))
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   ! Extract predicted temperature at next time step from Equation of State

   ISOTHERMAL_IF_1: IF (ISOTHERMAL) THEN
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
               RHOS(I,J,K) = PBAR_S(K,PRESSURE_ZONE(I,J,K))/(RSUM(I,J,K)*TMP(I,J,K))
            ENDDO
         ENDDO
      ENDDO
   ELSE ISOTHERMAL_IF_1
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
               TMP(I,J,K) = PBAR_S(K,PRESSURE_ZONE(I,J,K))/(RSUM(I,J,K)*RHOS(I,J,K))
            ENDDO
         ENDDO
      ENDDO
   ENDIF ISOTHERMAL_IF_1

   TMP = MAX(TMPMIN,MIN(TMPMAX,TMP))

! The CORRECTOR step

CASE(.FALSE.) PREDICTOR_STEP

   ! Correct species mass fraction at next time step (ZZ here also stores RHO*ZZ)

   DO N=1,N_TRACKED_SPECIES
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
               ZZ(I,J,K,N) = .5_EB*(RHO(I,J,K)*ZZ(I,J,K,N) + RHOS(I,J,K)*ZZS(I,J,K,N) - DT*DEL_RHO_D_DEL_Z(I,J,K,N))
            ENDDO
         ENDDO
      ENDDO
   ENDDO

   ! Correct density at next time step

   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
            RHO(I,J,K) = .5_EB*(RHO(I,J,K)+RHOS(I,J,K)-DT*FRHO(I,J,K))
         ENDDO
      ENDDO
   ENDDO

   ! Correct densities above or below clip limits

   IF (.NOT.CLIP_MASS_FRACTION) CALL CHECK_DENSITY

   ! Extract Y_n from rho*Y_n

   DO N=1,N_TRACKED_SPECIES
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
               ZZ(I,J,K,N) = ZZ(I,J,K,N)/RHO(I,J,K)
            ENDDO
         ENDDO
      ENDDO
   ENDDO

   ! Correct background pressure

   DO I=1,N_ZONE
      PBAR(:,I) = 0.5_EB*(PBAR(:,I) + PBAR_S(:,I) + D_PBAR_DT_S(I)*DT)
   ENDDO

   ! Manufactured solution

   IF (PERIODIC_TEST==7) THEN
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               ZZ(I,J,K,1) = VD2D_MMS_Z_OF_RHO(RHO(I,J,K))
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   ! Correct mass fractions above or below clip limits

   IF (CLIP_MASS_FRACTION .AND. N_TRACKED_SPECIES>0) THEN
      ZZ(1:IBAR,1:JBAR,1:KBAR,1:N_TRACKED_SPECIES) = MAX(0._EB,MIN(1._EB,ZZ(1:IBAR,1:JBAR,1:KBAR,1:N_TRACKED_SPECIES)))
   ELSEIF (N_TRACKED_SPECIES>0) THEN
      IF (USE_NEW_CHECK_MASS_FRACTION) THEN
         CALL NEW_CHECK_MASS_FRACTION
      ELSE
         CALL CHECK_MASS_FRACTION
      ENDIF
   ENDIF

   ! Compute molecular weight term RSUM=R0*SUM(Y_i/W_i)

   IF (N_TRACKED_SPECIES>0) THEN
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
               ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(I,J,K,1:N_TRACKED_SPECIES)
               CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,RSUM(I,J,K))
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   ! Extract predicted temperature at next time step from Equation of State

   ISOTHERMAL_IF_2: IF (ISOTHERMAL) THEN
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
               RHO(I,J,K) = PBAR(K,PRESSURE_ZONE(I,J,K))/(RSUM(I,J,K)*TMP(I,J,K))
            ENDDO
         ENDDO
      ENDDO
   ELSE ISOTHERMAL_IF_2
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
               TMP(I,J,K) = PBAR(K,PRESSURE_ZONE(I,J,K))/(RSUM(I,J,K)*RHO(I,J,K))
            ENDDO
         ENDDO
      ENDDO
   ENDIF ISOTHERMAL_IF_2

   TMP = MAX(TMPMIN,MIN(TMPMAX,TMP))

END SELECT PREDICTOR_STEP

TUSED(3,NM)=TUSED(3,NM)+SECOND()-TNOW

END SUBROUTINE DENSITY


SUBROUTINE CHECK_DENSITY

! Redistribute mass from cells below or above the density cut-off limits
! Do not apply OpenMP to this routine

USE GLOBAL_CONSTANTS, ONLY : PREDICTOR,RHOMIN,RHOMAX
REAL(EB) :: MASS_N(-3:3),CONST,MASS_C,RHO_CUT,VC(-3:3),SIGN_FACTOR,SUM_MASS_N,VC1(-3:3)
INTEGER  :: IC,I,J,K
REAL(EB), POINTER, DIMENSION(:,:,:) :: DELTA_RHO=>NULL()

DELTA_RHO => WORK2
DELTA_RHO =  0._EB

IF (PREDICTOR) THEN
   RHOP=>RHOS
ELSE
   RHOP=>RHO
ENDIF

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

END SUBROUTINE CHECK_DENSITY


SUBROUTINE CHECK_MASS_FRACTION

! Redistribute species mass from cells below or above the cut-off limits
! Do not apply OpenMP to this routine

USE PHYSICAL_FUNCTIONS, ONLY: IS_REALIZABLE
USE GLOBAL_CONSTANTS, ONLY : PREDICTOR,N_TRACKED_SPECIES
REAL(EB) :: SUM,CONST,MASS_C,MASS_N(-3:3),ZZ_CUT,SIGN_FACTOR,SUM_MASS_N,VC(-3:3),VC1(-3:3),ZZ_GET(0:N_TRACKED_SPECIES)
INTEGER  :: IC,N,I,J,K
REAL(EB), POINTER, DIMENSION(:,:,:) :: DELTA_ZZ=>NULL()
LOGICAL :: REALIZABLE

DELTA_ZZ => WORK1

IF (PREDICTOR) THEN
   RHOP    => RHOS
   ZZP     => ZZS
ELSE
   RHOP    => RHO
   ZZP     => ZZ
ENDIF

SPECIES_LOOP: DO N=1,N_TRACKED_SPECIES

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

            IF (ZZP(I,J,K,N)>=0._EB .AND. ZZP(I,J,K,N)<=1._EB) CYCLE
            IC = CELL_INDEX(I,J,K)
            IF (SOLID(IC)) CYCLE
            IF (ZZP(I,J,K,N)<0._EB) THEN
               ZZ_CUT      =  0._EB
               SIGN_FACTOR =  1._EB
            ELSE
               ZZ_CUT      =  1._EB
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

            IF (WALL_INDEX(IC,-1)==0) MASS_N(-1) = RHOP(I-1,J,K)*ABS(MIN(1._EB,MAX(0._EB,ZZP(I-1,J,K,N)))-ZZ_CUT)*VC(-1)
            IF (WALL_INDEX(IC, 1)==0) MASS_N( 1) = RHOP(I+1,J,K)*ABS(MIN(1._EB,MAX(0._EB,ZZP(I+1,J,K,N)))-ZZ_CUT)*VC( 1)
            IF (WALL_INDEX(IC,-2)==0) MASS_N(-2) = RHOP(I,J-1,K)*ABS(MIN(1._EB,MAX(0._EB,ZZP(I,J-1,K,N)))-ZZ_CUT)*VC(-2)
            IF (WALL_INDEX(IC, 2)==0) MASS_N( 2) = RHOP(I,J+1,K)*ABS(MIN(1._EB,MAX(0._EB,ZZP(I,J+1,K,N)))-ZZ_CUT)*VC( 2)
            IF (WALL_INDEX(IC,-3)==0) MASS_N(-3) = RHOP(I,J,K-1)*ABS(MIN(1._EB,MAX(0._EB,ZZP(I,J,K-1,N)))-ZZ_CUT)*VC(-3)
            IF (WALL_INDEX(IC, 3)==0) MASS_N( 3) = RHOP(I,J,K+1)*ABS(MIN(1._EB,MAX(0._EB,ZZP(I,J,K+1,N)))-ZZ_CUT)*VC( 3)

            SUM_MASS_N = SUM(MASS_N)
            IF (SUM_MASS_N<=TWO_EPSILON_EB) CYCLE
            MASS_C = RHOP(I,J,K)*ABS(ZZP(I,J,K,N)-ZZ_CUT)*VC(0)
            CONST  = SIGN_FACTOR*MIN(1._EB,MASS_C/SUM_MASS_N)
            DELTA_ZZ(I,J,K)   = DELTA_ZZ(I,J,K)   + CONST*SUM_MASS_N/(RHOP(I,J,K)  *VC( 0))
            DELTA_ZZ(I-1,J,K) = DELTA_ZZ(I-1,J,K) - CONST*MASS_N(-1)/(RHOP(I-1,J,K)*VC(-1))
            DELTA_ZZ(I+1,J,K) = DELTA_ZZ(I+1,J,K) - CONST*MASS_N( 1)/(RHOP(I+1,J,K)*VC( 1))
            DELTA_ZZ(I,J-1,K) = DELTA_ZZ(I,J-1,K) - CONST*MASS_N(-2)/(RHOP(I,J-1,K)*VC(-2))
            DELTA_ZZ(I,J+1,K) = DELTA_ZZ(I,J+1,K) - CONST*MASS_N( 2)/(RHOP(I,J+1,K)*VC( 2))
            DELTA_ZZ(I,J,K-1) = DELTA_ZZ(I,J,K-1) - CONST*MASS_N(-3)/(RHOP(I,J,K-1)*VC(-3))
            DELTA_ZZ(I,J,K+1) = DELTA_ZZ(I,J,K+1) - CONST*MASS_N( 3)/(RHOP(I,J,K+1)*VC( 3))
         ENDDO
      ENDDO
   ENDDO

   ZZP(1:IBAR,1:JBAR,1:KBAR,N) = ZZP(1:IBAR,1:JBAR,1:KBAR,N) + DELTA_ZZ(1:IBAR,1:JBAR,1:KBAR)

   IF (CHECK_REALIZABILITY) THEN
      REALIZABLE=.TRUE.
      OUTER_LOOP: DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I,J,K,1:N_TRACKED_SPECIES)
               ZZ_GET(0) = 1._EB-SUM(ZZ_GET(1:N_TRACKED_SPECIES))
               REALIZABLE=IS_REALIZABLE(ZZ_GET)
               IF (.NOT.REALIZABLE) THEN
                  WRITE(LU_ERR,*) 'ERROR: Unrealizable mass fracitons in CHECK_MASS_FRACTION'
                  STOP_STATUS=REALIZABILITY_STOP
               ENDIF
            ENDDO
         ENDDO
      ENDDO OUTER_LOOP
   ENDIF

ENDDO SPECIES_LOOP

END SUBROUTINE CHECK_MASS_FRACTION


SUBROUTINE NEW_CHECK_MASS_FRACTION

USE COMP_FUNCTIONS, ONLY: SHUTDOWN
USE PHYSICAL_FUNCTIONS, ONLY: IS_REALIZABLE,GET_REALIZABLE_MF
USE GLOBAL_CONSTANTS, ONLY : PREDICTOR,N_TRACKED_SPECIES
REAL(EB) :: ZZ_GET(0:N_TRACKED_SPECIES)
INTEGER  :: I,J,K

IF (PREDICTOR) THEN
   ZZP => ZZS
ELSE
   ZZP => ZZ
ENDIF

DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR

         IF (SOLID(CELL_INDEX(I,J,K))) CYCLE

         ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I,J,K,1:N_TRACKED_SPECIES)
         ZZ_GET(0) = 1._EB - SUM(ZZ_GET(1:N_TRACKED_SPECIES))
         IF (IS_REALIZABLE(ZZ_GET)) CYCLE

         CALL GET_REALIZABLE_MF(ZZ_GET)

         ! double check realizability
         ZZ_GET(0) = 1._EB - SUM(ZZ_GET(1:N_TRACKED_SPECIES))
         IF (.NOT.IS_REALIZABLE(ZZ_GET)) THEN
            WRITE(LU_ERR,*) ZZ_GET,SUM(ZZ_GET)
            WRITE(LU_ERR,*) 'ERROR: Unrealizable mass fractions in NEW_CHECK_MASS_FRACTION'
            STOP_STATUS=REALIZABILITY_STOP
         ENDIF

         ZZP(I,J,K,1:N_TRACKED_SPECIES) = ZZ_GET(1:N_TRACKED_SPECIES)

      ENDDO
   ENDDO
ENDDO

END SUBROUTINE NEW_CHECK_MASS_FRACTION


SUBROUTINE MASS_FINITE_DIFFERENCES_2(NM)

! Compute spatial differences for density equation

USE COMP_FUNCTIONS, ONLY: SECOND
USE GLOBAL_CONSTANTS, ONLY: N_TRACKED_SPECIES,PREDICTOR,EVACUATION_ONLY,SOLID_PHASE_ONLY,TUSED
USE SOOT_ROUTINES, ONLY: SETTLING_VELOCITY

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
   IF (N_TRACKED_SPECIES > 0) ZZP => ZZ
ELSE
   UU => US
   VV => VS
   WW => WS
   RHOP => RHOS
   IF (N_TRACKED_SPECIES > 0) ZZP => ZZS
ENDIF

DENSITY_TRANSPORT_IF: IF (DENSITY_TRANSPORT) THEN

! Compute scalar values at cell faces, FX, FY, FZ

!$OMP PARALLEL PRIVATE(ZZZ)
!$OMP DO SCHEDULE(STATIC)
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBM1
         ZZZ(1:4) = RHOP(I-1:I+2,J,K)
         FX(I,J,K,0) = SCALAR_FACE_VALUE(UU(I,J,K),ZZZ,FLUX_LIMITER)
      ENDDO
   ENDDO
ENDDO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO K=1,KBAR
   DO J=1,JBM1
      DO I=1,IBAR
         ZZZ(1:4) = RHOP(I,J-1:J+2,K)
         FY(I,J,K,0) = SCALAR_FACE_VALUE(VV(I,J,K),ZZZ,FLUX_LIMITER)
      ENDDO
   ENDDO
ENDDO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO K=1,KBM1
   DO J=1,JBAR
      DO I=1,IBAR
         ZZZ(1:4) = RHOP(I,J,K-1:K+2)
         FZ(I,J,K,0) = SCALAR_FACE_VALUE(WW(I,J,K),ZZZ,FLUX_LIMITER)
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

   SELECT CASE(IOR)
      CASE( 1)
         FX(IIG-1,JJG,KKG,0) = WC%RHO_F
      CASE(-1)
         FX(IIG,JJG,KKG,0)   = WC%RHO_F
      CASE( 2)
         FY(IIG,JJG-1,KKG,0) = WC%RHO_F
      CASE(-2)
         FY(IIG,JJG,KKG,0)   = WC%RHO_F
      CASE( 3)
         FZ(IIG,JJG,KKG-1,0) = WC%RHO_F
      CASE(-3)
         FZ(IIG,JJG,KKG,0)   = WC%RHO_F
   END SELECT

   ! Overwrite first off-wall advective flux if flow is away from the wall and if the face is not also a wall cell

   OFF_WALL_IF: IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY .AND. WC%BOUNDARY_TYPE/=OPEN_BOUNDARY) THEN

      OFF_WALL_SELECT: SELECT CASE(IOR)
         CASE( 1) OFF_WALL_SELECT
            !      ghost          FX/UU(II+1)
            ! ///   II   ///  II+1  |  II+2  | ...
            !                       ^ WALL_INDEX(II+1,+1)
            IF ((UU(II+1,JJ,KK)>0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II+1,JJ,KK),+1)>0)) THEN
               ZZZ(1:3) = (/WC%RHO_F,RHOP(II+1:II+2,JJ,KK)/)
               FX(II+1,JJ,KK,0) = SCALAR_FACE_VALUE(UU(II+1,JJ,KK),ZZZ,FLUX_LIMITER)
            ENDIF
         CASE(-1) OFF_WALL_SELECT
            !            FX/UU(II-2)     ghost
            ! ... |  II-2  |  II-1  ///   II   ///
            !              ^ WALL_INDEX(II-1,-1)
            IF ((UU(II-2,JJ,KK)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II-1,JJ,KK),-1)>0)) THEN
               ZZZ(2:4) = (/RHOP(II-2:II-1,JJ,KK),WC%RHO_F/)
               FX(II-2,JJ,KK,0) = SCALAR_FACE_VALUE(UU(II-2,JJ,KK),ZZZ,FLUX_LIMITER)
            ENDIF
         CASE( 2) OFF_WALL_SELECT
            IF ((VV(II,JJ+1,KK)>0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ+1,KK),+2)>0)) THEN
               ZZZ(1:3) = (/WC%RHO_F,RHOP(II,JJ+1:JJ+2,KK)/)
               FY(II,JJ+1,KK,0) = SCALAR_FACE_VALUE(VV(II,JJ+1,KK),ZZZ,FLUX_LIMITER)
            ENDIF
         CASE(-2) OFF_WALL_SELECT
            IF ((VV(II,JJ-2,KK)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ-1,KK),-2)>0)) THEN
               ZZZ(2:4) = (/RHOP(II,JJ-2:JJ-1,KK),WC%RHO_F/)
               FY(II,JJ-2,KK,0) = SCALAR_FACE_VALUE(VV(II,JJ-2,KK),ZZZ,FLUX_LIMITER)
            ENDIF
         CASE( 3) OFF_WALL_SELECT
            IF ((WW(II,JJ,KK+1)>0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ,KK+1),+3)>0)) THEN
               ZZZ(1:3) = (/WC%RHO_F,RHOP(II,JJ,KK+1:KK+2)/)
               FZ(II,JJ,KK+1,0) = SCALAR_FACE_VALUE(WW(II,JJ,KK+1),ZZZ,FLUX_LIMITER)
            ENDIF
         CASE(-3) OFF_WALL_SELECT
            IF ((WW(II,JJ,KK-2)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ,KK-1),-3)>0)) THEN
               ZZZ(2:4) = (/RHOP(II,JJ,KK-2:KK-1),WC%RHO_F/)
               FZ(II,JJ,KK-2,0) = SCALAR_FACE_VALUE(WW(II,JJ,KK-2),ZZZ,FLUX_LIMITER)
            ENDIF
      END SELECT OFF_WALL_SELECT
   
   ENDIF OFF_WALL_IF

ENDDO WALL_LOOP

ENDIF DENSITY_TRANSPORT_IF

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

END SUBROUTINE MASS_FINITE_DIFFERENCES_2


SUBROUTINE DENSITY_2(T,NM)

! Update the species mass fractions and density

USE COMP_FUNCTIONS, ONLY: SECOND,SHUTDOWN
USE PHYSICAL_FUNCTIONS, ONLY : GET_SPECIFIC_GAS_CONSTANT,GET_SENSIBLE_ENTHALPY,GET_SPECIFIC_HEAT,GET_SENSIBLE_ENTHALPY_DIFF,&
                               GET_REALIZABLE_MF
USE GLOBAL_CONSTANTS, ONLY: N_TRACKED_SPECIES,TMPMAX,TMPMIN,EVACUATION_ONLY, &
                            PREDICTOR,N_ZONE,GAS_SPECIES,R0,SOLID_PHASE_ONLY,TUSED
USE MANUFACTURED_SOLUTIONS, ONLY: VD2D_MMS_Z_OF_RHO,VD2D_MMS_Z_SRC,UF_MMS,WF_MMS,VD2D_MMS_RHO_OF_Z,VD2D_MMS_Z_SRC
INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: T
REAL(EB) :: TNOW,ZZ_GET(0:N_TRACKED_SPECIES),RHS,UN,Q_Z,XHAT,ZHAT
INTEGER :: I,J,K,N,IW,IOR,IIG,JJG,KKG
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: RHO_ZZ__0=>NULL(),JX=>NULL(),JY=>NULL(),JZ=>NULL()
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU=>NULL(),VV=>NULL(),WW=>NULL(),RHO__0=>NULL()
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
RHO__0=>WORK4
RHO_ZZ__0=>SCALAR_WORK1

PREDICTOR_STEP: SELECT CASE (PREDICTOR)

CASE(.TRUE.) PREDICTOR_STEP

   ! Diffusive fluxes

   JX=>SCALAR_SAVE1
   JY=>SCALAR_SAVE2
   JZ=>SCALAR_SAVE3
   IF (.NOT.CHANGE_TIME_STEP(NM)) THEN
      JX=FDX
      JY=FDY
      JZ=FDZ
   ENDIF

   ! Correct boundary velocity at wall cells

   UU=U
   VV=V
   WW=W

   WALL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
      WC=>WALL(IW)
      IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE WALL_LOOP

      IIG = WC%ONE_D%IIG 
      JJG = WC%ONE_D%JJG
      KKG = WC%ONE_D%KKG
      IOR = WC%ONE_D%IOR

      BOUNDARY_SELECT: SELECT CASE(WC%BOUNDARY_TYPE)
         CASE(INTERPOLATED_BOUNDARY,SOLID_BOUNDARY)
            SELECT CASE(WC%BOUNDARY_TYPE)
               CASE(SOLID_BOUNDARY)
                  IF (PREDICTOR) UN = -SIGN(1._EB,REAL(IOR,EB))*WC%ONE_D%UW
                  IF (CORRECTOR) UN = -SIGN(1._EB,REAL(IOR,EB))*WC%ONE_D%UWS
               CASE(INTERPOLATED_BOUNDARY)
                  UN = UVW_SAVE(IW)
            END SELECT
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
      END SELECT BOUNDARY_SELECT
   ENDDO WALL_LOOP

   ! Predictor step for mass density

   DO N=1,N_TRACKED_SPECIES
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE

               RHS = ((FX(I,J,K,N)*UU(I,J,K)-JX(I,J,K,N))*R(I)-(FX(I-1,J,K,N)*UU(I-1,J,K)-JX(I-1,J,K,N))*R(I-1))*RDX(I)*RRN(I) &
                   + ((FY(I,J,K,N)*VV(I,J,K)-JY(I,J,K,N))     -(FY(I,J-1,K,N)*VV(I,J-1,K)-JY(I,J-1,K,N))       )*RDY(J)        &
                   + ((FZ(I,J,K,N)*WW(I,J,K)-JZ(I,J,K,N))     -(FZ(I,J,K-1,N)*WW(I,J,K-1)-JZ(I,J,K-1,N))       )*RDZ(K)

               RHO_ZZ__0(I,J,K,N) = RHO(I,J,K)*ZZ(I,J,K,N)

               ZZS(I,J,K,N) = RHO_ZZ__0(I,J,K,N) - DT*RHS
            ENDDO
         ENDDO
      ENDDO
   ENDDO

   ! Manufactured solution

   IF (PERIODIC_TEST==7 .AND. .NOT.DENSITY_TRANSPORT) THEN
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               ! divergence from EOS
               XHAT = XC(I) - UF_MMS*(T)
               ZHAT = ZC(K) - WF_MMS*(T)
               Q_Z = VD2D_MMS_Z_SRC(XHAT,ZHAT,T)
               ZZS(I,J,K,1) = ZZS(I,J,K,1) + DT*Q_Z
               ZZS(I,J,K,2) = ZZS(I,J,K,2) - DT*Q_Z
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   IF (DENSITY_TRANSPORT) THEN
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE

               RHS = (FX(I,J,K,0)*UU(I,J,K)*R(I)-FX(I-1,J,K,0)*UU(I-1,J,K)*R(I-1))*RDX(I)*RRN(I) &
                   + (FY(I,J,K,0)*VV(I,J,K)     -FY(I,J-1,K,0)*VV(I,J-1,K)       )*RDY(J)        &
                   + (FZ(I,J,K,0)*WW(I,J,K)     -FZ(I,J,K-1,0)*WW(I,J,K-1)       )*RDZ(K)

               RHO__0(I,J,K) = RHO(I,J,K)

               RHOS(I,J,K) = RHO__0(I,J,K) - DT*RHS
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   ! Correct fluxes for positivity

   CALL WEIGHTED_AVERAGE_FLUX_CORRECTION

   ! Get rho = sum(rho*Y_alpha)

   IF (.NOT.DENSITY_TRANSPORT) THEN
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
               RHOS(I,J,K) = SUM(ZZS(I,J,K,1:N_TRACKED_SPECIES))
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   ! Extract mass fraction from RHO * ZZ

   IF (N_TRACKED_SPECIES>0) THEN
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
               ZZ_GET(1:N_TRACKED_SPECIES) = ZZS(I,J,K,1:N_TRACKED_SPECIES)/RHOS(I,J,K)
               IF (DENSITY_TRANSPORT) CALL GET_REALIZABLE_MF(ZZ_GET)
               ZZS(I,J,K,1:N_TRACKED_SPECIES) = ZZ_GET(1:N_TRACKED_SPECIES)
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   ! Predict background pressure at next time step

   DO I=1,N_ZONE
      PBAR_S(:,I) = PBAR(:,I) + D_PBAR_DT(I)*DT
   ENDDO

   ! Manufactured solution

   IF (PERIODIC_TEST==7 .AND. DENSITY_TRANSPORT) THEN
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               ZZS(I,J,K,1) = VD2D_MMS_Z_OF_RHO(RHOS(I,J,K))
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   ! Compute molecular weight term RSUM=R0*SUM(Y_i/W_i)

   IF (N_TRACKED_SPECIES>0) THEN
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
               ZZ_GET(1:N_TRACKED_SPECIES) = ZZS(I,J,K,1:N_TRACKED_SPECIES)
               CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,RSUM(I,J,K))
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   ! Extract predicted temperature at next time step from Equation of State

   ISOTHERMAL_IF_1: IF (ISOTHERMAL) THEN
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
               RHOS(I,J,K) = PBAR_S(K,PRESSURE_ZONE(I,J,K))/(RSUM(I,J,K)*TMP(I,J,K))
            ENDDO
         ENDDO
      ENDDO
   ELSE ISOTHERMAL_IF_1
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
               TMP(I,J,K) = PBAR_S(K,PRESSURE_ZONE(I,J,K))/(RSUM(I,J,K)*RHOS(I,J,K))
            ENDDO
         ENDDO
      ENDDO
   ENDIF ISOTHERMAL_IF_1

   TMP = MAX(TMPMIN,MIN(TMPMAX,TMP))

! The CORRECTOR step

CASE(.FALSE.) PREDICTOR_STEP

   ! Correct boundary velocity at wall cells

   UU=US
   VV=VS
   WW=WS

   WALL_LOOP_2: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
      WC=>WALL(IW)
      IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE WALL_LOOP_2

      IIG = WC%ONE_D%IIG 
      JJG = WC%ONE_D%JJG
      KKG = WC%ONE_D%KKG
      IOR = WC%ONE_D%IOR

      BOUNDARY_SELECT_2: SELECT CASE(WC%BOUNDARY_TYPE)
         CASE(INTERPOLATED_BOUNDARY,SOLID_BOUNDARY)
            SELECT CASE(WC%BOUNDARY_TYPE)
               CASE(SOLID_BOUNDARY)
                  IF (PREDICTOR) UN = -SIGN(1._EB,REAL(IOR,EB))*WC%ONE_D%UW
                  IF (CORRECTOR) UN = -SIGN(1._EB,REAL(IOR,EB))*WC%ONE_D%UWS
               CASE(INTERPOLATED_BOUNDARY)
                  UN = UVW_SAVE(IW)
            END SELECT
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
      END SELECT BOUNDARY_SELECT_2
   ENDDO WALL_LOOP_2

   ! Diffusive fluxes

   JX=>SCALAR_SAVE1
   JY=>SCALAR_SAVE2
   JZ=>SCALAR_SAVE3
   JX=FDX
   JY=FDY
   JZ=FDZ

   ! Compute species mass density at the next time step

   DO N=1,N_TRACKED_SPECIES
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE

               RHS = ((FX(I,J,K,N)*UU(I,J,K)-JX(I,J,K,N))*R(I)-(FX(I-1,J,K,N)*UU(I-1,J,K)-JX(I-1,J,K,N))*R(I-1))*RDX(I)*RRN(I) &
                   + ((FY(I,J,K,N)*VV(I,J,K)-JY(I,J,K,N))     -(FY(I,J-1,K,N)*VV(I,J-1,K)-JY(I,J-1,K,N))       )*RDY(J)        &
                   + ((FZ(I,J,K,N)*WW(I,J,K)-JZ(I,J,K,N))     -(FZ(I,J,K-1,N)*WW(I,J,K-1)-JZ(I,J,K-1,N))       )*RDZ(K)

               RHO_ZZ__0(I,J,K,N) = .5_EB*(RHO(I,J,K)*ZZ(I,J,K,N) + RHOS(I,J,K)*ZZS(I,J,K,N))

               ZZ(I,J,K,N) = RHO_ZZ__0(I,J,K,N) - .5_EB*DT*RHS
            ENDDO
         ENDDO
      ENDDO
   ENDDO

   ! Manufactured solution

   IF (PERIODIC_TEST==7 .AND. .NOT.DENSITY_TRANSPORT) THEN
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               ! divergence from EOS
               XHAT = XC(I) - UF_MMS*T
               ZHAT = ZC(K) - WF_MMS*T
               Q_Z = VD2D_MMS_Z_SRC(XHAT,ZHAT,T)
               ZZ(I,J,K,1) = ZZ(I,J,K,1) + .5_EB*DT*Q_Z
               ZZ(I,J,K,2) = ZZ(I,J,K,2) - .5_EB*DT*Q_Z
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   IF (DENSITY_TRANSPORT) THEN
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE

               RHS = (FX(I,J,K,0)*UU(I,J,K)*R(I)-FX(I-1,J,K,0)*UU(I-1,J,K)*R(I-1))*RDX(I)*RRN(I) &
                   + (FY(I,J,K,0)*VV(I,J,K)     -FY(I,J-1,K,0)*VV(I,J-1,K)       )*RDY(J)        &
                   + (FZ(I,J,K,0)*WW(I,J,K)     -FZ(I,J,K-1,0)*WW(I,J,K-1)       )*RDZ(K)

               RHO__0(I,J,K) = .5_EB*(RHO(I,J,K)+RHOS(I,J,K))

               RHO(I,J,K) = RHO__0(I,J,K) - .5_EB*DT*RHS
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   ! Correct fluxes for positivity

   CALL WEIGHTED_AVERAGE_FLUX_CORRECTION

   ! Get rho = sum(rho*Y_alpha)

   IF (.NOT.DENSITY_TRANSPORT) THEN
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
               RHO(I,J,K) = SUM(ZZ(I,J,K,1:N_TRACKED_SPECIES))
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   ! Extract Y_n from rho*Y_n

   IF (N_TRACKED_SPECIES>0) THEN
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
               ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(I,J,K,1:N_TRACKED_SPECIES)/RHO(I,J,K)
               IF (DENSITY_TRANSPORT) CALL GET_REALIZABLE_MF(ZZ_GET)
               ZZ(I,J,K,1:N_TRACKED_SPECIES) = ZZ_GET(1:N_TRACKED_SPECIES)
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   ! Correct background pressure

   DO I=1,N_ZONE
      PBAR(:,I) = 0.5_EB*(PBAR(:,I) + PBAR_S(:,I) + D_PBAR_DT_S(I)*DT)
   ENDDO

   ! Manufactured solution

   IF (PERIODIC_TEST==7 .AND. DENSITY_TRANSPORT) THEN
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               ZZ(I,J,K,1) = VD2D_MMS_Z_OF_RHO(RHO(I,J,K))
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   ! Compute molecular weight term RSUM=R0*SUM(Y_i/W_i)

   IF (N_TRACKED_SPECIES>0) THEN
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
               ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(I,J,K,1:N_TRACKED_SPECIES)
               CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,RSUM(I,J,K))
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   ! Extract predicted temperature at next time step from Equation of State

   ISOTHERMAL_IF_2: IF (ISOTHERMAL) THEN
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
               RHO(I,J,K) = PBAR(K,PRESSURE_ZONE(I,J,K))/(RSUM(I,J,K)*TMP(I,J,K))
            ENDDO
         ENDDO
      ENDDO
   ELSE ISOTHERMAL_IF_2
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
               TMP(I,J,K) = PBAR(K,PRESSURE_ZONE(I,J,K))/(RSUM(I,J,K)*RHO(I,J,K))
            ENDDO
         ENDDO
      ENDDO
   ENDIF ISOTHERMAL_IF_2

   TMP = MAX(TMPMIN,MIN(TMPMAX,TMP))

END SELECT PREDICTOR_STEP

TUSED(3,NM)=TUSED(3,NM)+SECOND()-TNOW

END SUBROUTINE DENSITY_2


SUBROUTINE WEIGHTED_AVERAGE_FLUX_CORRECTION

USE COMP_FUNCTIONS, ONLY: SHUTDOWN

INTEGER :: I,J,K,N
REAL(EB) :: SVDT,GAMMA,RHO_ZZ_GAMMA,RHO_GAMMA,RHS,DT_LOC
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU=>NULL(),VV=>NULL(),WW=>NULL(),RHOP=>NULL(),RHO__0=>NULL()
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZP=>NULL(),RHO_ZZ__0=>NULL(),JX=>NULL(),JY=>NULL(),JZ=>NULL()
INTEGER, POINTER, DIMENSION(:,:,:) :: IOB=>NULL()

IF (PREDICTOR) THEN
   ZZP=>ZZS
   RHOP=>RHOS
   DT_LOC=DT
ELSE
   ZZP=>ZZ
   RHOP=>RHO
   DT_LOC=.5_EB*DT
ENDIF
UU=>WORK1
VV=>WORK2
WW=>WORK3
RHO__0=>WORK4
RHO_ZZ__0=>SCALAR_WORK1
IOB=>IWORK1
JX=>SCALAR_SAVE1
JY=>SCALAR_SAVE2
JZ=>SCALAR_SAVE3

! Species

SPECIES_LOOP: DO N=1,N_TRACKED_SPECIES

   IOB=0 ! integer tag for out of bounds (0=in bounds, 1=out of bounds)

   ! Compute new fluxes

   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
            IF (ZZP(I,J,K,N)>=0._EB) CYCLE
            
            IOB(I,J,K)=1 ! cell is tagged for correction

            SVDT = DT_LOC * ( MAX(0._EB,UU(I,J,K)) - MIN(0._EB,UU(I-1,J,K)) &
                            + MAX(0._EB,VV(I,J,K)) - MIN(0._EB,VV(I,J-1,K)) &
                            + MAX(0._EB,WW(I,J,K)) - MIN(0._EB,WW(I,J,K-1)) )

            GAMMA = (1._EB-EXP(-SVDT))/(SVDT)

            RHO_ZZ_GAMMA = RHO_ZZ__0(I,J,K,N)*GAMMA

            IF ((FX(I,J,K,N)-JX(I,J,K,N))>0._EB) THEN
               JX(I,J,K,N) = 0._EB
               FX(I,J,K,N)  = RHO_ZZ_GAMMA
               IOB(I+1,J,K) = 1 ! right neighbor is tagged for correction
            ENDIF
            IF ((FY(I,J,K,N)-JY(I,J,K,N))>0._EB) THEN
               JY(I,J,K,N) = 0._EB
               FY(I,J,K,N)  = RHO_ZZ_GAMMA
               IOB(I,J+1,K) = 1
            ENDIF
            IF ((FZ(I,J,K,N)-JZ(I,J,K,N))>0._EB) THEN
               JZ(I,J,K,N) = 0._EB
               FZ(I,J,K,N)  = RHO_ZZ_GAMMA
               IOB(I,J,K+1) = 1
            ENDIF

            IF ((FX(I-1,J,K,N)-JX(I-1,J,K,N))<0._EB) THEN
               JX(I-1,J,K,N) = 0._EB
               FX(I-1,J,K,N) = RHO_ZZ_GAMMA
               IOB(I-1,J,K) = 1
            ENDIF
            IF ((FY(I,J-1,K,N)-JY(I,J-1,K,N))<0._EB) THEN
               JY(I,J-1,K,N) = 0._EB
               FY(I,J-1,K,N) = RHO_ZZ_GAMMA
               IOB(I,J-1,K) = 1
            ENDIF
            IF ((FZ(I,J,K-1,N)-JZ(I,J,K-1,N))<0._EB) THEN
               JZ(I,J,K-1,N) = 0._EB
               FZ(I,J,K-1,N) = RHO_ZZ_GAMMA
               IOB(I,J,K-1) = 1
            ENDIF

         ENDDO
      ENDDO
   ENDDO

   ! Only update cells tagged for correction

   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
            IF (IOB(I,J,K)/=1) CYCLE

            RHS = ((FX(I,J,K,N)*UU(I,J,K)-JX(I,J,K,N))*R(I)-(FX(I-1,J,K,N)*UU(I-1,J,K)-JX(I-1,J,K,N))*R(I-1))*RDX(I)*RRN(I) &
                + ((FY(I,J,K,N)*VV(I,J,K)-JY(I,J,K,N))     -(FY(I,J-1,K,N)*VV(I,J-1,K)-JY(I,J-1,K,N))       )*RDY(J)        &
                + ((FZ(I,J,K,N)*WW(I,J,K)-JZ(I,J,K,N))     -(FZ(I,J,K-1,N)*WW(I,J,K-1)-JZ(I,J,K-1,N))       )*RDZ(K)

            ZZP(I,J,K,N) = RHO_ZZ__0(I,J,K,N) - DT_LOC*RHS

         ENDDO
      ENDDO
   ENDDO

ENDDO SPECIES_LOOP

DENSITY_TRANSPORT_IF: IF (DENSITY_TRANSPORT) THEN

! Density

IOB=0

! Compute new fluxes

DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
         IF (RHOP(I,J,K)>=0._EB) CYCLE
         
         IOB(I,J,K)=1 ! cell is tagged for correction

         SVDT = DT_LOC * ( MAX(0._EB,UU(I,J,K)) - MIN(0._EB,UU(I-1,J,K)) &
                         + MAX(0._EB,VV(I,J,K)) - MIN(0._EB,VV(I,J-1,K)) &
                         + MAX(0._EB,WW(I,J,K)) - MIN(0._EB,WW(I,J,K-1)) )

         GAMMA = (1._EB-EXP(-SVDT))/(SVDT)

         RHO_GAMMA = RHO__0(I,J,K)*GAMMA

         IF (FX(I,J,K,0)>0._EB) THEN
            FX(I,J,K,0)  = RHO_GAMMA
            IOB(I+1,J,K) = 1 ! right neighbor is tagged for correction
         ENDIF
         IF (FY(I,J,K,0)>0._EB) THEN
            FY(I,J,K,0)  = RHO_GAMMA
            IOB(I,J+1,K) = 1
         ENDIF
         IF (FZ(I,J,K,0)>0._EB) THEN
            FZ(I,J,K,0)  = RHO_GAMMA
            IOB(I,J,K+1) = 1
         ENDIF

         IF (FX(I-1,J,K,0)<0._EB) THEN
            FX(I-1,J,K,0) = RHO_GAMMA
            IOB(I-1,J,K) = 1
         ENDIF
         IF (FY(I,J-1,K,0)<0._EB) THEN
            FY(I,J-1,K,0) = RHO_GAMMA
            IOB(I,J-1,K) = 1
         ENDIF
         IF (FZ(I,J,K-1,0)<0._EB) THEN
            FZ(I,J,K-1,0) = RHO_GAMMA
            IOB(I,J,K-1) = 1
         ENDIF

      ENDDO
   ENDDO
ENDDO

! Only update cells tagged for correction

DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
         IF (IOB(I,J,K)/=1) CYCLE

         RHS = (FX(I,J,K,0)*UU(I,J,K)*R(I)-FX(I-1,J,K,0)*UU(I-1,J,K)*R(I-1))*RDX(I)*RRN(I) &
             + (FY(I,J,K,0)*VV(I,J,K)     -FY(I,J-1,K,0)*VV(I,J-1,K)       )*RDY(J)        &
             + (FZ(I,J,K,0)*WW(I,J,K)     -FZ(I,J,K-1,0)*WW(I,J,K-1)       )*RDZ(K)

         RHOP(I,J,K) = RHO__0(I,J,K) - DT_LOC*RHS

      ENDDO
   ENDDO
ENDDO

ENDIF DENSITY_TRANSPORT_IF

END SUBROUTINE WEIGHTED_AVERAGE_FLUX_CORRECTION


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
