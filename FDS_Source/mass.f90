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

PUBLIC MASS_FINITE_DIFFERENCES,DENSITY,GET_REV_mass
 
 
CONTAINS
 

SUBROUTINE MASS_FINITE_DIFFERENCES(NM)

! Compute spatial differences for density equation

USE COMP_FUNCTIONS, ONLY: SECOND
USE GLOBAL_CONSTANTS, ONLY: N_TRACKED_SPECIES,NULL_BOUNDARY,POROUS_BOUNDARY,OPEN_BOUNDARY,INTERPOLATED_BOUNDARY, &
                            PREDICTOR,CORRECTOR,EVACUATION_ONLY,SOLID_PHASE_ONLY,TUSED,DEBUG_OPENMP,SOLID_BOUNDARY, &
                            NO_MASS_FLUX,SPECIFIED_MASS_FLUX,HVAC_BOUNDARY
INTEGER, INTENT(IN) :: NM
REAL(EB) :: TNOW,ZZ(4),UN,RHO_D_DYDN
INTEGER  :: I,J,K,N,II,JJ,KK,IIG,JJG,KKG,IW,IOR,IBC
REAL(EB), POINTER, DIMENSION(:) :: UWP
REAL(EB), POINTER, DIMENSION(:,:,:) :: FX=>NULL(),FY=>NULL(),FZ=>NULL(),MASS_COR=>NULL()
 
IF (EVACUATION_ONLY(NM) .OR. SOLID_PHASE_ONLY) RETURN

TNOW=SECOND()
CALL POINT_TO_MESH(NM)
 
IF (PREDICTOR) THEN
   UU => U
   VV => V
   WW => W
   DP => D
   RHOP => RHO
   UWP  => UW
ELSE
   UU => US
   VV => VS
   WW => WS
   DP => DS
   RHOP => RHOS
   UWP  => UWS
ENDIF

!$OMP PARALLEL 

!$OMP SINGLE
FX=>WORK4
FY=>WORK5
FZ=>WORK6
MASS_COR=>WORK7
MASS_COR=0._EB
!$OMP END SINGLE

!$OMP DO COLLAPSE(3) PRIVATE(K,J,I,ZZ)
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBM1
         ZZ(1:4) = RHOP(I-1:I+2,J,K)
         FX(I,J,K) = UU(I,J,K)*SCALAR_FACE_VALUE(UU(I,J,K),ZZ,FLUX_LIMITER)*R(I)
      ENDDO
   ENDDO
ENDDO
!$OMP END DO

!$OMP DO COLLAPSE(3) PRIVATE(K,J,I,ZZ)
DO K=1,KBAR
   DO J=1,JBM1
      DO I=1,IBAR
         ZZ(1:4) = RHOP(I,J-1:J+2,K)
         FY(I,J,K) = VV(I,J,K)*SCALAR_FACE_VALUE(VV(I,J,K),ZZ,FLUX_LIMITER)
      ENDDO
   ENDDO
ENDDO
!$OMP END DO

!$OMP DO COLLAPSE(3) PRIVATE(K,J,I,ZZ)
DO K=1,KBM1
   DO J=1,JBAR
      DO I=1,IBAR
         ZZ(1:4) = RHOP(I,J,K-1:K+2)
         FZ(I,J,K) = WW(I,J,K)*SCALAR_FACE_VALUE(WW(I,J,K),ZZ,FLUX_LIMITER)
      ENDDO
   ENDDO
ENDDO
!$OMP END DO

!$OMP SINGLE PRIVATE(IW,II,JJ,KK,IOR,IBC,IIG,JJG,KKG,UN,ZZ)
!!$OMP DO PRIVATE(IW,II,JJ,KK,IOR,IBC,IIG,JJG,KKG,UN,ZZ)
WLOOP_FL: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
   
   IF (BOUNDARY_TYPE(IW)==NULL_BOUNDARY .OR. BOUNDARY_TYPE(IW)==POROUS_BOUNDARY) CYCLE WLOOP_FL
       
   II  = IJKW(1,IW) 
   JJ  = IJKW(2,IW)
   KK  = IJKW(3,IW)
   IOR = IJKW(4,IW)
   IBC = IJKW(5,IW)
   IIG = IJKW(6,IW)
   JJG = IJKW(7,IW)
   KKG = IJKW(8,IW)
   
   ! overwrite first off-wall advective flux if flow is away from the wall and if the face is not also a wall cell

   OFF_WALL_SELECT_1: SELECT CASE(IOR)
      CASE( 1) OFF_WALL_SELECT_1
         !      ghost          FX/UU(II+1)
         ! ///   II   ///  II+1  |  II+2  | ...
         !                       ^ WALL_INDEX(II+1,+1)
         IF ((UU(II+1,JJ,KK)>0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II+1,JJ,KK),+1)>0)) THEN
            ZZ(1:3) = (/RHO_F(IW),RHOP(II+1:II+2,JJ,KK)/)
            FX(II+1,JJ,KK) = UU(II+1,JJ,KK)*SCALAR_FACE_VALUE(UU(II+1,JJ,KK),ZZ,FLUX_LIMITER)*R(II+1)
         ENDIF
      CASE(-1) OFF_WALL_SELECT_1
         !            FX/UU(II-2)     ghost
         ! ... |  II-2  |  II-1  ///   II   ///
         !              ^ WALL_INDEX(II-1,-1)
         IF ((UU(II-2,JJ,KK)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II-1,JJ,KK),-1)>0)) THEN
            ZZ(2:4) = (/RHOP(II-2:II-1,JJ,KK),RHO_F(IW)/)
            FX(II-2,JJ,KK) = UU(II-2,JJ,KK)*SCALAR_FACE_VALUE(UU(II-2,JJ,KK),ZZ,FLUX_LIMITER)*R(II-2)
         ENDIF
      CASE( 2) OFF_WALL_SELECT_1
         IF ((VV(II,JJ+1,KK)>0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ+1,KK),+2)>0)) THEN
            ZZ(1:3) = (/RHO_F(IW),RHOP(II,JJ+1:JJ+2,KK)/)
            FY(II,JJ+1,KK) = VV(II,JJ+1,KK)*SCALAR_FACE_VALUE(VV(II,JJ+1,KK),ZZ,FLUX_LIMITER)
         ENDIF
      CASE(-2) OFF_WALL_SELECT_1
         IF ((VV(II,JJ-2,KK)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ-1,KK),-2)>0)) THEN
            ZZ(2:4) = (/RHOP(II,JJ-2:JJ-1,KK),RHO_F(IW)/)
            FY(II,JJ-2,KK) = VV(II,JJ-2,KK)*SCALAR_FACE_VALUE(VV(II,JJ-2,KK),ZZ,FLUX_LIMITER)
         ENDIF
      CASE( 3) OFF_WALL_SELECT_1
         IF ((WW(II,JJ,KK+1)>0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ,KK+1),+3)>0)) THEN
            ZZ(1:3) = (/RHO_F(IW),RHOP(II,JJ,KK+1:KK+2)/)
            FZ(II,JJ,KK+1) = WW(II,JJ,KK+1)*SCALAR_FACE_VALUE(WW(II,JJ,KK+1),ZZ,FLUX_LIMITER)
         ENDIF
      CASE(-3) OFF_WALL_SELECT_1
         IF ((WW(II,JJ,KK-2)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ,KK-1),-3)>0)) THEN
            ZZ(2:4) = (/RHOP(II,JJ,KK-2:KK-1),RHO_F(IW)/)
            FZ(II,JJ,KK-2) = WW(II,JJ,KK-2)*SCALAR_FACE_VALUE(WW(II,JJ,KK-2),ZZ,FLUX_LIMITER)
         ENDIF
   END SELECT OFF_WALL_SELECT_1
   
   SELECT CASE(IOR)
      CASE( 1)
         UN = UU(II,JJ,KK)
      CASE(-1)
         UN = UU(II-1,JJ,KK)
      CASE( 2)
         UN = VV(II,JJ,KK)
      CASE(-2)
         UN = VV(II,JJ-1,KK)
      CASE( 3)
         UN = WW(II,JJ,KK)
      CASE(-3)
         UN = WW(II,JJ,KK-1)
   END SELECT
   IF (BOUNDARY_TYPE(IW)==INTERPOLATED_BOUNDARY) UN = UVW_SAVE(IW)
   
   MASS_COR_IF: IF (SURFACE(IBC)%SPECIES_BC_INDEX==NO_MASS_FLUX) THEN
      UN=0._EB
      SELECT CASE(IOR)
         CASE( 1)
            MASS_COR(IIG,JJG,KKG) = MASS_COR(IIG,JJG,KKG) + RHO_F(IW)*(UU(II,JJ,KK)-UN)*RDX(IIG)
         CASE(-1)
            MASS_COR(IIG,JJG,KKG) = MASS_COR(IIG,JJG,KKG) - RHO_F(IW)*(UU(II-1,JJ,KK)-UN)*RDX(IIG)
         CASE( 2)
            MASS_COR(IIG,JJG,KKG) = MASS_COR(IIG,JJG,KKG) + RHO_F(IW)*(VV(II,JJ,KK)-UN)*RDY(JJG)
         CASE(-2)
            MASS_COR(IIG,JJG,KKG) = MASS_COR(IIG,JJG,KKG) - RHO_F(IW)*(VV(II,JJ-1,KK)-UN)*RDY(JJG)
         CASE( 3)
            MASS_COR(IIG,JJG,KKG) = MASS_COR(IIG,JJG,KKG) + RHO_F(IW)*(WW(II,JJ,KK)-UN)*RDZ(KKG)
         CASE(-3)
            MASS_COR(IIG,JJG,KKG) = MASS_COR(IIG,JJG,KKG) - RHO_F(IW)*(WW(II,JJ,KK-1)-UN)*RDZ(KKG)
      END SELECT
   ENDIF MASS_COR_IF
   
   ! compute flux on the face of the wall cell

   SELECT CASE(IOR)
      CASE( 1)
         FX(II,JJ,KK)   = UN*RHO_F(IW)*R(II)
      CASE(-1)
         FX(II-1,JJ,KK) = UN*RHO_F(IW)*R(II-1)
      CASE( 2)
         FY(II,JJ,KK)   = UN*RHO_F(IW)
      CASE(-2)
         FY(II,JJ-1,KK) = UN*RHO_F(IW)
      CASE( 3)
         FZ(II,JJ,KK)   = UN*RHO_F(IW)
      CASE(-3)
         FZ(II,JJ,KK-1) = UN*RHO_F(IW)
   END SELECT
      
ENDDO WLOOP_FL
!!$OMP END DO
!$OMP END SINGLE

!$OMP WORKSHARE
FRHO = 0._EB
!$OMP END WORKSHARE

!$OMP DO COLLAPSE(3) PRIVATE(K,J,I)
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
         FRHO(I,J,K) = (FX(I,J,K)-FX(I-1,J,K))*RDX(I)*RRN(I) &
                     + (FY(I,J,K)-FY(I,J-1,K))*RDY(J)        &
                     + (FZ(I,J,K)-FZ(I,J,K-1))*RDZ(K)        &
                     - MASS_COR(I,J,K)
      ENDDO
   ENDDO
ENDDO
!$OMP END DO

!$OMP END PARALLEL 


! Compute the species equation differences
 
IF (N_TRACKED_SPECIES > 0) THEN
   IF (PREDICTOR) YYP => YY
   IF (CORRECTOR) YYP => YYS
ENDIF
 
SPECIES_LOOP: DO N=1,N_TRACKED_SPECIES

      FX=>WORK4
      FY=>WORK5
      FZ=>WORK6
   
      !$OMP PARALLEL
      !$OMP DO COLLAPSE(3) PRIVATE(K,J,I,ZZ)
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBM1
               ZZ(1:4) = RHOP(I-1:I+2,J,K)*YYP(I-1:I+2,J,K,N)
               FX(I,J,K) = UU(I,J,K)*SCALAR_FACE_VALUE(UU(I,J,K),ZZ,FLUX_LIMITER)*R(I)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO
      
      !$OMP DO COLLAPSE(3) PRIVATE(K,J,I,ZZ)
      DO K=1,KBAR
         DO J=1,JBM1
            DO I=1,IBAR
               ZZ(1:4) = RHOP(I,J-1:J+2,K)*YYP(I,J-1:J+2,K,N)
               FY(I,J,K) = VV(I,J,K)*SCALAR_FACE_VALUE(VV(I,J,K),ZZ,FLUX_LIMITER)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO
      
      !$OMP DO COLLAPSE(3) PRIVATE(K,J,I,ZZ)
      DO K=1,KBM1
         DO J=1,JBAR
            DO I=1,IBAR
               ZZ(1:4) = RHOP(I,J,K-1:K+2)*YYP(I,J,K-1:K+2,N)
               FZ(I,J,K) = WW(I,J,K)*SCALAR_FACE_VALUE(WW(I,J,K),ZZ,FLUX_LIMITER)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO

      
      !$OMP SINGLE PRIVATE(IW,II,JJ,KK,IOR,IBC,IIG,JJG,KKG,RHO_D_DYDN,UN,ZZ)
      !!$OMP DO PRIVATE(IW,II,JJ,KK,IOR,IBC,IIG,JJG,KKG,RHO_D_DYDN,UN,ZZ)
      WLOOP2_FL: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS

         IF (BOUNDARY_TYPE(IW)==NULL_BOUNDARY .OR. BOUNDARY_TYPE(IW)==POROUS_BOUNDARY) CYCLE WLOOP2_FL
             
         II  = IJKW(1,IW)
         JJ  = IJKW(2,IW)
         KK  = IJKW(3,IW) 
         IOR = IJKW(4,IW)
         IBC = IJKW(5,IW)
         IIG = IJKW(6,IW)
         JJG = IJKW(7,IW)
         KKG = IJKW(8,IW)
         
         ! overwrite first off-wall advective flux if flow is away from the wall and if the face is not also a wall cell

         OFF_WALL_SELECT_2: SELECT CASE(IOR)
            CASE( 1) OFF_WALL_SELECT_2
               !      ghost          FX/UU(II+1)
               ! ///   II   ///  II+1  |  II+2  | ...
               !                       ^ WALL_INDEX(II+1,+1)
               IF ((UU(II+1,JJ,KK)>0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II+1,JJ,KK),+1)>0)) THEN
                  ZZ(1:3) = (/RHO_F(IW),RHOP(II+1:II+2,JJ,KK)/)*(/YY_F(IW,N),YYP(II+1:II+2,JJ,KK,N)/)
                  FX(II+1,JJ,KK) = UU(II+1,JJ,KK)*SCALAR_FACE_VALUE(UU(II+1,JJ,KK),ZZ,FLUX_LIMITER)*R(II+1)
               ENDIF
            CASE(-1) OFF_WALL_SELECT_2
               !            FX/UU(II-2)     ghost
               ! ... |  II-2  |  II-1  ///   II   ///
               !              ^ WALL_INDEX(II-1,-1)
               IF ((UU(II-2,JJ,KK)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II-1,JJ,KK),-1)>0)) THEN
                  ZZ(2:4) = (/RHOP(II-2:II-1,JJ,KK),RHO_F(IW)/)*(/YYP(II-2:II-1,JJ,KK,N),YY_F(IW,N)/)
                  FX(II-2,JJ,KK) = UU(II-2,JJ,KK)*SCALAR_FACE_VALUE(UU(II-2,JJ,KK),ZZ,FLUX_LIMITER)*R(II-2)
               ENDIF
            CASE( 2) OFF_WALL_SELECT_2
               IF ((VV(II,JJ+1,KK)>0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ+1,KK),+2)>0)) THEN
                  ZZ(1:3) = (/RHO_F(IW),RHOP(II,JJ+1:JJ+2,KK)/)*(/YY_F(IW,N),YYP(II,JJ+1:JJ+2,KK,N)/)
                  FY(II,JJ+1,KK) = VV(II,JJ+1,KK)*SCALAR_FACE_VALUE(VV(II,JJ+1,KK),ZZ,FLUX_LIMITER)
               ENDIF
            CASE(-2) OFF_WALL_SELECT_2
               IF ((VV(II,JJ-2,KK)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ-1,KK),-2)>0)) THEN
                  ZZ(2:4) = (/RHOP(II,JJ-2:JJ-1,KK),RHO_F(IW)/)*(/YYP(II,JJ-2:JJ-1,KK,N),YY_F(IW,N)/)
                  FY(II,JJ-2,KK) = VV(II,JJ-2,KK)*SCALAR_FACE_VALUE(VV(II,JJ-2,KK),ZZ,FLUX_LIMITER)
               ENDIF
            CASE( 3) OFF_WALL_SELECT_2
               IF ((WW(II,JJ,KK+1)>0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ,KK+1),+3)>0)) THEN
                  ZZ(1:3) = (/RHO_F(IW),RHOP(II,JJ,KK+1:KK+2)/)*(/YY_F(IW,N),YYP(II,JJ,KK+1:KK+2,N)/)
                  FZ(II,JJ,KK+1) = WW(II,JJ,KK+1)*SCALAR_FACE_VALUE(WW(II,JJ,KK+1),ZZ,FLUX_LIMITER)
               ENDIF
            CASE(-3) OFF_WALL_SELECT_2
               IF ((WW(II,JJ,KK-2)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ,KK-1),-3)>0)) THEN
                  ZZ(2:4) = (/RHOP(II,JJ,KK-2:KK-1),RHO_F(IW)/)*(/YYP(II,JJ,KK-2:KK-1,N),YY_F(IW,N)/)
                  FZ(II,JJ,KK-2) = WW(II,JJ,KK-2)*SCALAR_FACE_VALUE(WW(II,JJ,KK-2),ZZ,FLUX_LIMITER)
               ENDIF
         END SELECT OFF_WALL_SELECT_2

         
         SELECT CASE(IOR)
            CASE( 1)
               UN = UU(II,JJ,KK)
            CASE(-1)
               UN = UU(II-1,JJ,KK)
            CASE( 2)
               UN = VV(II,JJ,KK)
            CASE(-2)
               UN = VV(II,JJ-1,KK)
            CASE( 3)
               UN = WW(II,JJ,KK)
            CASE(-3)
               UN = WW(II,JJ,KK-1)
         END SELECT
         IF (BOUNDARY_TYPE(IW)==INTERPOLATED_BOUNDARY) UN = UVW_SAVE(IW)
         IF ((SURFACE(IBC)%SPECIES_BC_INDEX==SPECIFIED_MASS_FLUX .OR. &
             (SURFACE(IBC)%SPECIES_BC_INDEX==HVAC_BOUNDARY       .OR. &
              ANY(SURFACE(IBC)%LEAK_PATH>0._EB)) .AND. UWS(IW)<0._EB) .AND. YY_F(IW,N)>0._EB) THEN
            ! recreate diffusive flux from divg b/c UWP based on old RHODW
            RHO_D_DYDN = 2._EB*RHODW(IW,N)*(YYP(IIG,JJG,KKG,N)-YY_F(IW,N))*RDN(IW)
            UN = SIGN(1._EB,REAL(IOR,EB))*(MASSFLUX(IW,N) + RHO_D_DYDN)/(RHO_F(IW)*YY_F(IW,N))
         ENDIF
         ! compute flux on the face of the wall cell
         SELECT CASE(IOR)
            CASE( 1)
               FX(II,JJ,KK)   = UN*RHO_F(IW)*YY_F(IW,N)*R(II)
            CASE(-1)
               FX(II-1,JJ,KK) = UN*RHO_F(IW)*YY_F(IW,N)*R(II-1)
            CASE( 2)
               FY(II,JJ,KK)   = UN*RHO_F(IW)*YY_F(IW,N)
            CASE(-2)
               FY(II,JJ-1,KK) = UN*RHO_F(IW)*YY_F(IW,N)
            CASE( 3) 
               FZ(II,JJ,KK)   = UN*RHO_F(IW)*YY_F(IW,N)
            CASE(-3) 
               FZ(II,JJ,KK-1) = UN*RHO_F(IW)*YY_F(IW,N)
         END SELECT

      ENDDO WLOOP2_FL
      !!$OMP END DO
      !$OMP END SINGLE
   
      !$OMP DO COLLAPSE(3) PRIVATE(K,J,I)
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
               DEL_RHO_D_DEL_Y(I,J,K,N) = -DEL_RHO_D_DEL_Y(I,J,K,N)             &
                                        + (FX(I,J,K)-FX(I-1,J,K))*RDX(I)*RRN(I) &
                                        + (FY(I,J,K)-FY(I,J-1,K))*RDY(J)        &
                                        + (FZ(I,J,K)-FZ(I,J,K-1))*RDZ(K)
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
USE GLOBAL_CONSTANTS, ONLY: N_TRACKED_SPECIES,CO_PRODUCTION,I_PROG_F,I_PROG_CO,I_FUEL,TMPMAX,TMPMIN,EVACUATION_ONLY, &
                            PREDICTOR,CORRECTOR,CHANGE_TIME_STEP,TMPA,N_ZONE, &
                            GAS_SPECIES, R0,SOLID_PHASE_ONLY,TUSED, &
                            RSUM0,DEBUG_OPENMP,CLIP_MASS_FRACTION
REAL(EB) :: DTRATIO,OMDTRATIO,TNOW,YY_GET(1:N_TRACKED_SPECIES)
INTEGER  :: I,J,K,N
INTEGER, INTENT(IN) :: NM
 
IF (EVACUATION_ONLY(NM)) RETURN
IF (SOLID_PHASE_ONLY) RETURN

TNOW=SECOND()
CALL POINT_TO_MESH(NM)

PREDICTOR_STEP: SELECT CASE (PREDICTOR)

CASE(.TRUE.) PREDICTOR_STEP

   !$OMP PARALLEL

   IF (.NOT.CHANGE_TIME_STEP(NM)) THEN

      !$OMP DO COLLAPSE(4) PRIVATE(N,K,J,I) SCHEDULE(DYNAMIC)
      DO N=1,N_TRACKED_SPECIES
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
                  YYS(I,J,K,N) = RHO(I,J,K)*YY(I,J,K,N) - DT*DEL_RHO_D_DEL_Y(I,J,K,N)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO

   ELSE

      !$OMP SINGLE
      DTRATIO   = DT/DT_PREV
      OMDTRATIO = 1._EB - DTRATIO
      !$OMP END SINGLE
      !$OMP DO COLLAPSE(4) PRIVATE(N,K,J,I) SCHEDULE(DYNAMIC)
      DO N=1,N_TRACKED_SPECIES
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
                  YYS(I,J,K,N) = OMDTRATIO*RHO(I,J,K)*YY(I,J,K,N) + DTRATIO*RHOS(I,J,K)*YYS(I,J,K,N)
               ENDDO
           ENDDO
         ENDDO
      ENDDO
      !$OMP END DO

   ENDIF
   !$OMP END PARALLEL

   ! Predict the density at the next time step (RHOS or RHO^*)

   !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(K,J,I)
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            RHOS(I,J,K) = RHO(I,J,K)-DT*FRHO(I,J,K)
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO

   ! Correct densities above or below clip limits

   CALL CHECK_DENSITY
   
   ! Extract mass fraction from RHO * YY

   !$OMP PARALLEL DO COLLAPSE(4) PRIVATE(N,K,J,I) SCHEDULE(DYNAMIC)
   DO N=1,N_TRACKED_SPECIES
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

   IF (CLIP_MASS_FRACTION) THEN
      YYS(1:IBAR,1:JBAR,1:KBAR,1:N_TRACKED_SPECIES) = MAX(0._EB,MIN(1._EB,YYS(1:IBAR,1:JBAR,1:KBAR,1:N_TRACKED_SPECIES)))
   ELSE
      CALL CHECK_MASS_FRACTION
   ENDIF

   ! Predict background pressure at next time step

   DO I=1,N_ZONE
      PBAR_S(:,I) = PBAR(:,I) + D_PBAR_DT(I)*DT
   ENDDO

   ! Compute molecular weight term RSUM=R0*SUM(Y_i/M_i)

   !$OMP PARALLEL SHARED(RSUM)
   IF (N_TRACKED_SPECIES>0) THEN
      !$OMP DO COLLAPSE(3) PRIVATE(K,J,I,YY_GET)
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR   
               YY_GET(:) = YYS(I,J,K,:)
               CALL GET_SPECIFIC_GAS_CONSTANT(YY_GET,RSUM(I,J,K)) !INTENT: IN,OUT
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO
   ENDIF

   ! Extract predicted temperature at next time step from Equation of State
   
   IF (N_TRACKED_SPECIES==0) THEN
      !$OMP DO COLLAPSE(3) PRIVATE(K,J,I)
      DO K=0,KBP1
         DO J=0,JBP1
            DO I=0,IBP1
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
               TMP(I,J,K) = PBAR_S(K,PRESSURE_ZONE(I,J,K))/(RSUM(I,J,K)*RHOS(I,J,K))
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO
   ENDIF
   !$OMP WORKSHARE
   TMP = MAX(TMPMIN,MIN(TMPMAX,TMP))
   !$OMP END WORKSHARE
   !$OMP END PARALLEL


! The CORRECTOR step
   
CASE(.FALSE.) PREDICTOR_STEP

   ! Correct species mass fraction at next time step (YY here actually means YY*RHO)

   !$OMP PARALLEL DO COLLAPSE(4) PRIVATE(N,K,J,I) SCHEDULE(DYNAMIC)
   DO N=1,N_TRACKED_SPECIES
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
               YY(I,J,K,N) = .5_EB*(RHO(I,J,K)*YY(I,J,K,N) + RHOS(I,J,K)*YYS(I,J,K,N) - DT*DEL_RHO_D_DEL_Y(I,J,K,N) ) 
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO

   ! Correct density at next time step

   !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(K,J,I)
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            RHO(I,J,K) = .5_EB*(RHO(I,J,K)+RHOS(I,J,K)-DT*FRHO(I,J,K))
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO

   ! Correct densities above or below clip limits

   CALL CHECK_DENSITY
   
   ! Extract Y_n from rho*Y_n

   !$OMP PARALLEL DO COLLAPSE(4) PRIVATE(N,K,J,I) SCHEDULE(DYNAMIC)
   DO N=1,N_TRACKED_SPECIES
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
               YY(I,J,K,N) = YY(I,J,K,N)/RHO(I,J,K)
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO

   ! Correct mass fractions above or below clip limits

   IF (CLIP_MASS_FRACTION) THEN
      YY(1:IBAR,1:JBAR,1:KBAR,1:N_TRACKED_SPECIES) = MAX(0._EB,MIN(1._EB,YY(1:IBAR,1:JBAR,1:KBAR,1:N_TRACKED_SPECIES)))
   ELSE
      CALL CHECK_MASS_FRACTION
   ENDIF

   ! Correct background pressure

   DO I=1,N_ZONE
      PBAR(:,I) = .5_EB*(PBAR(:,I) + PBAR_S(:,I) + D_PBAR_S_DT(I)*DT)
   ENDDO
 
   ! Compute molecular weight term RSUM=R0*SUM(Y_i/M_i)

   !$OMP PARALLEL SHARED(RSUM)
   IF (N_TRACKED_SPECIES>0) THEN
      !$OMP DO COLLAPSE(3) PRIVATE(K,J,I,YY_GET)
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR   
               YY_GET(:) = YY(I,J,K,:)
               CALL GET_SPECIFIC_GAS_CONSTANT(YY_GET,RSUM(I,J,K)) !INTENT: IN,OUT
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO
   ENDIF

   ! Extract predicted temperature at next time step from Equation of State

   IF (N_TRACKED_SPECIES==0) THEN
      !$OMP DO COLLAPSE(3) PRIVATE(K,J,I)
      DO K=0,KBP1
         DO J=0,JBP1
            DO I=0,IBP1
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
               TMP(I,J,K) = PBAR(K,PRESSURE_ZONE(I,J,K))/(RSUM(I,J,K)*RHO(I,J,K))
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO
   ENDIF

   !$OMP WORKSHARE
   TMP = MAX(TMPMIN,MIN(TMPMAX,TMP))
   !$OMP END WORKSHARE
   !$OMP END PARALLEL

END SELECT PREDICTOR_STEP

TUSED(3,NM)=TUSED(3,NM)+SECOND()-TNOW
 
END SUBROUTINE DENSITY
 

SUBROUTINE CHECK_DENSITY
 
! Redistribute mass from cells below or above the density cut-off limits

USE GLOBAL_CONSTANTS, ONLY : PREDICTOR, CORRECTOR, N_TRACKED_SPECIES,RHOMIN,RHOMAX,DEBUG_OPENMP,TWO_D
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
!$OMP END PARALLEL

!!$OMP DO COLLAPSE(3) PRIVATE(K,J,I,IC,RMIN,SUM,ISUM,LC,RHO00,RHOMI,RHOPI,RHOMJ,RHOPJ,RHOMK,RHOPK,CONST,CONST2) SCHEDULE(DYNAMIC)
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
            IF(ABS(SUM-ISUM*RHO00) >=ZERO_P) THEN
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
!!$OMP END DO

!$OMP PARALLEL WORKSHARE

RHOP(1:IBAR,1:JBAR,1:KBAR) = MAX(RHOMIN,RHOP(1:IBAR,1:JBAR,1:KBAR)+RHODELTA(1:IBAR,1:JBAR,1:KBAR))

! Correct overshoots

RHODELTA = 0._EB
!$OMP END PARALLEL WORKSHARE

!!$OMP DO COLLAPSE(3) PRIVATE(K,J,I,IC,RMAX,SUM,ISUM,LC,RHO00,RHOMI,RHOPI,RHOMJ,RHOPJ,RHOMK,RHOPK,CONST,CONST2) SCHEDULE(DYNAMIC)
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
            IF(ABS(SUM-ISUM*RHO00) >=ZERO_P) THEN         
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
!!$OMP END DO

!$OMP PARALLEL WORKSHARE
RHOP(1:IBAR,1:JBAR,1:KBAR) = MIN(RHOMAX,RHOP(1:IBAR,1:JBAR,1:KBAR)+RHODELTA(1:IBAR,1:JBAR,1:KBAR))
!$OMP END PARALLEL WORKSHARE
!!$OMP END PARALLEL
END SUBROUTINE CHECK_DENSITY
 
 
SUBROUTINE CHECK_MASS_FRACTION

! Redistribute species mass from cells below or above the cut-off limits

USE GLOBAL_CONSTANTS, ONLY : PREDICTOR, CORRECTOR, N_TRACKED_SPECIES,POROUS_BOUNDARY,DEBUG_OPENMP
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

SPECIESLOOP: DO N=1,N_TRACKED_SPECIES

   !!$OMP PARALLEL
   !$OMP PARALLEL WORKSHARE
   YYDELTA = 0._EB
   !$OMP END PARALLEL WORKSHARE

   ! Do undershoots

   !!$OMP DO COLLAPSE(3) PRIVATE(K,J,I,IC,IW_A,Y00,SUM,ISUM,LC,YMIN,YMI,YPI,YMJ,YPJ,YMK,YPK) &
   !!$OMP PRIVATE(RHY0,RHYPI,RHYMI,RHYPJ,RHYMJ,RHYPK,RHYMK,CONST) SCHEDULE(DYNAMIC)
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
            YMIN  = 1._EB 
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
            YMIN = MIN(YMI,YPI,YMJ,YPJ,YMK,YPK)
            YMIN = MAX(YMIN,0._EB)
            IF ((DEL_RHO_D_DEL_Y(I,J,K,N) > 0._EB .AND. Y00 < YMIN) .OR. Y00 < 0._EB) THEN
               RHY0  = RHOP(I,J,K)  *(YMIN - Y00)
               IF (LC(-1) .AND. YMI>YMIN) THEN
                  RHYMI = RHOP(I-1,J,K)*(YMI - YMIN)
                  SUM  = SUM + RHYMI 
                  ISUM = ISUM + 1
               ELSE
                  LC(-1) = .FALSE.
               ENDIF
               IF (LC( 1) .AND. YPI>YMIN) THEN
                  RHYPI = RHOP(I+1,J,K)*(YPI - YMIN)
                  SUM  = SUM + RHYPI
                  ISUM = ISUM + 1
               ELSE
                  LC( 1) = .FALSE.
               ENDIF
               IF (LC(-2) .AND. YMJ>YMIN) THEN
                  RHYMJ = RHOP(I,J-1,K)*(YMJ - YMIN)
                  SUM  = SUM + RHYMJ
                  ISUM = ISUM + 1
               ELSE
                  LC(-2) = .FALSE.
               ENDIF
               IF (LC( 2) .AND. YPJ>YMIN) THEN
                  RHYPJ = RHOP(I,J+1,K)*(YPJ - YMIN)
                  SUM  = SUM + RHYPJ
                  ISUM = ISUM + 1
                  LC( 2) = .TRUE.
               ELSE
                  LC( 2) = .FALSE.
               ENDIF
               IF (LC(-3) .AND. YMK>YMIN) THEN
               RHYMK = RHOP(I,J,K-1)*(YMK - YMIN)
                  SUM  = SUM + RHYMK
                  ISUM = ISUM + 1
               ELSE
                  LC(-3) = .FALSE.
               ENDIF
               IF (LC( 3) .AND. YPK>YMIN) THEN
                  RHYPK = RHOP(I,J,K+1)*(YPK - YMIN)
                  SUM  = SUM + RHYPK
                  ISUM = ISUM + 1
               ELSE
                  LC( 3) = .FALSE.
               ENDIF                
               IF (ISUM==0) THEN
                  IF (YMIN <= 0._EB) YYDELTA(I,J,K) = YYDELTA(I,J,K) + YMIN - Y00  
                  CYCLE CHECK_LOOP
               ELSE
                  IF (ABS(SUM)>=ZERO_P) THEN
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
   !!$OMP END DO
   
   !$OMP PARALLEL WORKSHARE
   YYP(1:IBAR,1:JBAR,1:KBAR,N) = YYP(1:IBAR,1:JBAR,1:KBAR,N) + YYDELTA(1:IBAR,1:JBAR,1:KBAR)
   YYDELTA=0._EB
   !$OMP END PARALLEL WORKSHARE

   ! Do overshoots
   !!$OMP DO COLLAPSE(3) PRIVATE(K,J,I,IC,IW_A,Y00,SUM,ISUM,LC,YMIN,YMI,YPI,YMK,YPK,YMJ,YPJ,YMAX) &
   !!$OMP PRIVATE(RHY0,RHYPI,RHYMI,RHYPJ,RHYMJ,RHYPK,RHYMK,CONST) SCHEDULE(DYNAMIC)
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
            YMIN  = 1._EB 
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
            YMAX = MAX(YMI,YPI,YMJ,YPJ,YMK,YPK)
            YMAX = MIN(YMAX,1._EB)            
            IF ((DEL_RHO_D_DEL_Y(I,J,K,N) < 0._EB .AND. Y00 > YMAX) .OR. Y00 > 1._EB) THEN
               RHY0  = RHOP(I,J,K)  *(Y00 - YMAX)
               IF (LC(-1) .AND. YMI<YMAX) THEN
                  RHYMI = RHOP(I-1,J,K)*(YMAX - YMI)
                  SUM  = SUM + RHYMI
                  ISUM = ISUM + 1
               ELSE
                  LC(-1) = .FALSE.
               ENDIF
               IF (LC( 1) .AND. YPI<YMAX) THEN
                  RHYPI = RHOP(I+1,J,K)*(YMAX - YPI)
                  SUM  = SUM + RHYPI
                  ISUM = ISUM + 1
               ELSE
                  LC( 1) = .FALSE.
               ENDIF
               IF (LC(-2) .AND. YMJ<YMAX) THEN
                  RHYMJ = RHOP(I,J-1,K)*(YMAX - YMJ)
                  SUM  = SUM + RHYMJ
                  ISUM = ISUM + 1
               ELSE
                  LC(-2) = .FALSE.
               ENDIF
               IF (LC( 2) .AND. YPJ<YMAX) THEN
                  RHYPJ = RHOP(I,J+1,K)*(YMAX - YPJ)
                  SUM  = SUM + RHYPJ
                  ISUM = ISUM + 1
               ELSE
                  LC( 2) = .FALSE.
               ENDIF
               IF (LC(-3) .AND. YMK<YMAX) THEN
                  RHYMK = RHOP(I,J,K-1)*(YMAX - YMK)
                  SUM  = SUM + RHYMK
                  ISUM = ISUM + 1
               ELSE
                  LC(-3) = .FALSE.
               ENDIF
               IF (LC( 3) .AND. YPK<YMAX) THEN
                  RHYPK = RHOP(I,J,K+1)*(YMAX - YPK)
                  SUM  = SUM + RHYPK
                  ISUM = ISUM + 1
               ELSE
                  LC( 3) = .FALSE.
               ENDIF                      
               IF (ISUM==0) THEN
                  IF(YMAX >= 1._EB) YYDELTA(I,J,K) = YYDELTA(I,J,K) + YMAX - Y00
                  CYCLE CHECK_LOOP2
               ELSE
                  IF (ABS(SUM)>=ZERO_P) THEN
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
   !!$OMP END DO 

   !$OMP PARALLEL WORKSHARE
   YYP(1:IBAR,1:JBAR,1:KBAR,N) = YYP(1:IBAR,1:JBAR,1:KBAR,N) + YYDELTA(1:IBAR,1:JBAR,1:KBAR)
   !$OMP END PARALLEL WORKSHARE
   !!$OMP END PARALLEL
ENDDO SPECIESLOOP

RETURN

END SUBROUTINE CHECK_MASS_FRACTION


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
         IF (ABS(DU_LOC)>0._EB) R = DU_UP/DU_LOC
         B = MAX(0._EB,MIN(2._EB*R,1._EB),MIN(R,2._EB))
         SCALAR_FACE_VALUE = U(2) + 0.5_EB*B*DU_LOC
      CASE(3) ! MINMOD
         IF (ABS(DU_LOC)>0._EB) R = DU_UP/DU_LOC
         B = MAX(0._EB,MIN(1._EB,R))
         SCALAR_FACE_VALUE = U(2) + 0.5_EB*B*DU_LOC
      CASE(4) ! CHARM
         IF (ABS(DU_UP)>0._EB) R = DU_LOC/DU_UP
         IF (R>0._EB) B = R*(3._EB*R+1._EB)/((R+1._EB)**2)
         SCALAR_FACE_VALUE = U(2) + 0.5_EB*B*DU_UP
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
         IF (ABS(DU_LOC)>0._EB) R = DU_UP/DU_LOC
         B = MAX(0._EB,MIN(2._EB*R,1._EB),MIN(R,2._EB))
         SCALAR_FACE_VALUE = U(3) - 0.5_EB*B*DU_LOC
      CASE(3) ! MINMOD
         IF (ABS(DU_LOC)>0._EB) R = DU_UP/DU_LOC
         B = MAX(0._EB,MIN(1._EB,R))
         SCALAR_FACE_VALUE = U(3) - 0.5_EB*B*DU_LOC
      CASE(4) ! CHARM
         IF (ABS(DU_UP)>0._EB) R = DU_LOC/DU_UP
         IF (R>0._EB) B = R*(3._EB*R+1._EB)/((R+1._EB)**2)
         SCALAR_FACE_VALUE = U(3) - 0.5_EB*B*DU_UP
    END SELECT
    
ENDIF WIND_DIRECTION_IF

END FUNCTION SCALAR_FACE_VALUE



SUBROUTINE GET_REV_mass(MODULE_REV,MODULE_DATE)
INTEGER,INTENT(INOUT) :: MODULE_REV
CHARACTER(255),INTENT(INOUT) :: MODULE_DATE
WRITE(MODULE_DATE,'(A)') massrev(INDEX(massrev,':')+1:LEN_TRIM(massrev)-2)
READ (MODULE_DATE,'(I5)') MODULE_REV
WRITE(MODULE_DATE,'(A)') massdate
END SUBROUTINE GET_REV_mass
 
END MODULE MASS
