MODULE MASS
 
! Compute the mass equation differences 
 
USE PRECISION_PARAMETERS
USE MESH_POINTERS

IMPLICIT NONE
PRIVATE
CHARACTER(255), PARAMETER :: massid='$Id$'
CHARACTER(255), PARAMETER :: massrev='$Revision$'
CHARACTER(255), PARAMETER :: massdate='$Date$'

REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZP
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW,RHOP,DP

PUBLIC MASS_FINITE_DIFFERENCES,DENSITY,GET_REV_mass
 
 
CONTAINS
 

SUBROUTINE MASS_FINITE_DIFFERENCES(NM)

! Compute spatial differences for density equation

USE COMP_FUNCTIONS, ONLY: SECOND
USE GLOBAL_CONSTANTS, ONLY: N_TRACKED_SPECIES,NULL_BOUNDARY,OPEN_BOUNDARY,INTERPOLATED_BOUNDARY, &
                            PREDICTOR,CORRECTOR,EVACUATION_ONLY,SOLID_PHASE_ONLY,TUSED,DEBUG_OPENMP,SOLID_BOUNDARY, &
                            NO_MASS_FLUX,SPECIFIED_MASS_FLUX,HVAC_BOUNDARY, &
                            INCLUDE_NUMERICAL_DIFFUSION,FLUX_LIMITER
INTEGER, INTENT(IN) :: NM
REAL(EB) :: TNOW,ZZZ(4),UN,RHO_D_DZDN
INTEGER  :: I,J,K,N,II,JJ,KK,IIG,JJG,KKG,IW,IOR,SURF_INDEX
REAL(EB), POINTER, DIMENSION(:,:,:) :: FX=>NULL(),FY=>NULL(),FZ=>NULL()
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZP=>NULL()
TYPE(WALL_TYPE), POINTER :: WC=>NULL()

IF (EVACUATION_ONLY(NM) .OR. SOLID_PHASE_ONLY) RETURN

TNOW=SECOND()
CALL POINT_TO_MESH(NM)
 
IF (PREDICTOR) THEN
   UU => U
   VV => V
   WW => W
   DP => D
   RHOP => RHO
   IF (N_TRACKED_SPECIES > 0) ZZP => ZZ
ELSE
   UU => US
   VV => VS
   WW => WS
   DP => DS
   RHOP => RHOS
   IF (N_TRACKED_SPECIES > 0) ZZP => ZZS
ENDIF

FX=>WORK4
FY=>WORK5
FZ=>WORK6
FX=0._EB
FY=0._EB
FZ=0._EB

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(KBAR,JBAR,IBAR,KBM1,JBM1,IBM1,RHOP,FX,FY,FZ,UU,VV,WW,FLUX_LIMITER,R, &
!$OMP        N_EXTERNAL_WALL_CELLS,N_INTERNAL_WALL_CELLS,WALL_INDEX,CELL_INDEX, &
!$OMP        RHO_F,UVW_SAVE, &
!$OMP        FRHO,SOLID,RDX,RDY,RDZ,RRN, &
!$OMP        N_TRACKED_SPECIES,ZZP,ZZ_F)

!$OMP DO COLLAPSE(3) SCHEDULE(STATIC) PRIVATE(K,J,I,ZZZ)
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBM1
         ZZZ(1:4) = RHOP(I-1:I+2,J,K)
         FX(I,J,K) = UU(I,J,K)*SCALAR_FACE_VALUE(UU(I,J,K),ZZZ,FLUX_LIMITER)*R(I)
      ENDDO
   ENDDO
ENDDO
!$OMP END DO NOWAIT

!$OMP DO COLLAPSE(3) SCHEDULE(STATIC) PRIVATE(K,J,I,ZZZ)
DO K=1,KBAR
   DO J=1,JBM1
      DO I=1,IBAR
         ZZZ(1:4) = RHOP(I,J-1:J+2,K)
         FY(I,J,K) = VV(I,J,K)*SCALAR_FACE_VALUE(VV(I,J,K),ZZZ,FLUX_LIMITER)
      ENDDO
   ENDDO
ENDDO
!$OMP END DO NOWAIT

!$OMP DO COLLAPSE(3) SCHEDULE(STATIC) PRIVATE(K,J,I,ZZZ)
DO K=1,KBM1
   DO J=1,JBAR
      DO I=1,IBAR
         ZZZ(1:4) = RHOP(I,J,K-1:K+2)
         FZ(I,J,K) = WW(I,J,K)*SCALAR_FACE_VALUE(WW(I,J,K),ZZZ,FLUX_LIMITER)
      ENDDO
   ENDDO
ENDDO
!$OMP END DO

!$OMP DO PRIVATE(IW,II,JJ,KK,IOR,SURF_INDEX,IIG,JJG,KKG,ZZZ,UN)
WLOOP_FL: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
   WC=>WALL(IW)
   IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE WLOOP_FL
       
   II  = WC%II 
   JJ  = WC%JJ
   KK  = WC%KK
   IOR = WC%IOR
   SURF_INDEX = WC%SURF_INDEX
   IIG = WC%IIG
   JJG = WC%JJG
   KKG = WC%KKG
   
   ! overwrite first off-wall advective flux if flow is away from the wall and if the face is not also a wall cell

   IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY .AND. WC%BOUNDARY_TYPE/=OPEN_BOUNDARY) THEN

      OFF_WALL_SELECT_1: SELECT CASE(IOR)
         CASE( 1) OFF_WALL_SELECT_1
            !      ghost          FX/UU(II+1)
            ! ///   II   ///  II+1  |  II+2  | ...
            !                       ^ WALL_INDEX(II+1,+1)
            IF ((UU(II+1,JJ,KK)>0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II+1,JJ,KK),+1)>0)) THEN
               ZZZ(1:3) = (/WC%RHO_F,RHOP(II+1:II+2,JJ,KK)/)
               FX(II+1,JJ,KK) = UU(II+1,JJ,KK)*SCALAR_FACE_VALUE(UU(II+1,JJ,KK),ZZZ,FLUX_LIMITER)*R(II+1)
            ENDIF
         CASE(-1) OFF_WALL_SELECT_1
            !            FX/UU(II-2)     ghost
            ! ... |  II-2  |  II-1  ///   II   ///
            !              ^ WALL_INDEX(II-1,-1)
            IF ((UU(II-2,JJ,KK)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II-1,JJ,KK),-1)>0)) THEN
               ZZZ(2:4) = (/RHOP(II-2:II-1,JJ,KK),WC%RHO_F/)
               FX(II-2,JJ,KK) = UU(II-2,JJ,KK)*SCALAR_FACE_VALUE(UU(II-2,JJ,KK),ZZZ,FLUX_LIMITER)*R(II-2)
            ENDIF
         CASE( 2) OFF_WALL_SELECT_1
            IF ((VV(II,JJ+1,KK)>0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ+1,KK),+2)>0)) THEN
               ZZZ(1:3) = (/WC%RHO_F,RHOP(II,JJ+1:JJ+2,KK)/)
               FY(II,JJ+1,KK) = VV(II,JJ+1,KK)*SCALAR_FACE_VALUE(VV(II,JJ+1,KK),ZZZ,FLUX_LIMITER)
            ENDIF
         CASE(-2) OFF_WALL_SELECT_1
            IF ((VV(II,JJ-2,KK)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ-1,KK),-2)>0)) THEN
               ZZZ(2:4) = (/RHOP(II,JJ-2:JJ-1,KK),WC%RHO_F/)
               FY(II,JJ-2,KK) = VV(II,JJ-2,KK)*SCALAR_FACE_VALUE(VV(II,JJ-2,KK),ZZZ,FLUX_LIMITER)
            ENDIF
         CASE( 3) OFF_WALL_SELECT_1
            IF ((WW(II,JJ,KK+1)>0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ,KK+1),+3)>0)) THEN
               ZZZ(1:3) = (/WC%RHO_F,RHOP(II,JJ,KK+1:KK+2)/)
               FZ(II,JJ,KK+1) = WW(II,JJ,KK+1)*SCALAR_FACE_VALUE(WW(II,JJ,KK+1),ZZZ,FLUX_LIMITER)
            ENDIF
         CASE(-3) OFF_WALL_SELECT_1
            IF ((WW(II,JJ,KK-2)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ,KK-1),-3)>0)) THEN
               ZZZ(2:4) = (/RHOP(II,JJ,KK-2:KK-1),WC%RHO_F/)
               FZ(II,JJ,KK-2) = WW(II,JJ,KK-2)*SCALAR_FACE_VALUE(WW(II,JJ,KK-2),ZZZ,FLUX_LIMITER)
            ENDIF
      END SELECT OFF_WALL_SELECT_1
   
   ENDIF
   
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

   ! In case of interpolated boundary, use the original velocity, not the averaged value

   IF (WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY) UN = UVW_SAVE(IW)
   
   ! Compute flux on the face of the wall cell

   SELECT CASE(IOR)
      CASE( 1)
         FX(II,JJ,KK)   = UN*WC%RHO_F*R(II)
      CASE(-1)
         FX(II-1,JJ,KK) = UN*WC%RHO_F*R(II-1)
      CASE( 2)
         FY(II,JJ,KK)   = UN*WC%RHO_F
      CASE(-2)
         FY(II,JJ-1,KK) = UN*WC%RHO_F
      CASE( 3)
         FZ(II,JJ,KK)   = UN*WC%RHO_F
      CASE(-3)
         FZ(II,JJ,KK-1) = UN*WC%RHO_F
   END SELECT
      
ENDDO WLOOP_FL
!$OMP END DO NOWAIT

!$OMP WORKSHARE
FRHO = 0._EB
!$OMP END WORKSHARE

!$OMP DO COLLAPSE(3) SCHEDULE(STATIC) PRIVATE(K,J,I)
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
         FRHO(I,J,K) = (FX(I,J,K)-FX(I-1,J,K))*RDX(I)*RRN(I) &
                     + (FY(I,J,K)-FY(I,J-1,K))*RDY(J)        &
                     + (FZ(I,J,K)-FZ(I,J,K-1))*RDZ(K) 
      ENDDO
   ENDDO
ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL 


! Compute the species equation differences
 
SPECIES_LOOP: DO N=1,N_TRACKED_SPECIES

      FX=0._EB
      FY=0._EB
      FZ=0._EB
      
      ! numerical diffusion (accounted for in divergence)
      
      IF (INCLUDE_NUMERICAL_DIFFUSION) THEN
         DFX(:,:,:,N)=0._EB
         DFY(:,:,:,N)=0._EB
         DFZ(:,:,:,N)=0._EB
      ENDIF
   
      !$OMP PARALLEL DEFAULT(NONE) &
      !$OMP SHARED(N,KBAR,JBAR,IBAR,KBM1,JBM1,IBM1,RHOP,ZZP,FX,FY,FZ,UU,VV,WW,FLUX_LIMITER,R, &
      !$OMP        N_EXTERNAL_WALL_CELLS,N_INTERNAL_WALL_CELLS, &
      !$OMP        WALL_INDEX,CELL_INDEX,RHO_F,ZZ_F,UVW_SAVE, &
      !$OMP        SURFACE,UWS,RHODW,RDN,MASSFLUX, &
      !$OMP        SOLID,DEL_RHO_D_DEL_Z,RDX,RDY,RDZ,RRN)

      !$OMP DO COLLAPSE(3) SCHEDULE(STATIC) PRIVATE(K,J,I,ZZZ)
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBM1
               ZZZ(1:4) = RHOP(I-1:I+2,J,K)*ZZP(I-1:I+2,J,K,N)
               FX(I,J,K) = UU(I,J,K)*SCALAR_FACE_VALUE(UU(I,J,K),ZZZ,FLUX_LIMITER)*R(I)
               IF (INCLUDE_NUMERICAL_DIFFUSION) DFX(I,J,K,N) = FX(I,J,K) - 0.5_EB*(ZZZ(2)+ZZZ(3))*UU(I,J,K)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO NOWAIT
      
      !$OMP DO COLLAPSE(3) SCHEDULE(STATIC) PRIVATE(K,J,I,ZZZ)
      DO K=1,KBAR
         DO J=1,JBM1
            DO I=1,IBAR
               ZZZ(1:4) = RHOP(I,J-1:J+2,K)*ZZP(I,J-1:J+2,K,N)
               FY(I,J,K) = VV(I,J,K)*SCALAR_FACE_VALUE(VV(I,J,K),ZZZ,FLUX_LIMITER)
               IF (INCLUDE_NUMERICAL_DIFFUSION) DFY(I,J,K,N) = FY(I,J,K) - 0.5_EB*(ZZZ(2)+ZZZ(3))*VV(I,J,K)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO NOWAIT
      
      !$OMP DO COLLAPSE(3) SCHEDULE(STATIC) PRIVATE(K,J,I,ZZZ)
      DO K=1,KBM1
         DO J=1,JBAR
            DO I=1,IBAR
               ZZZ(1:4) = RHOP(I,J,K-1:K+2)*ZZP(I,J,K-1:K+2,N)
               FZ(I,J,K) = WW(I,J,K)*SCALAR_FACE_VALUE(WW(I,J,K),ZZZ,FLUX_LIMITER)
               IF (INCLUDE_NUMERICAL_DIFFUSION) DFZ(I,J,K,N) = FZ(I,J,K) - 0.5_EB*(ZZZ(2)+ZZZ(3))*WW(I,J,K)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO

      !$OMP DO SCHEDULE(STATIC) &
      !$OMP PRIVATE(IW,II,JJ,KK,IOR,SURF_INDEX,IIG,JJG,KKG,ZZZ,UN,RHO_D_DZDN)
      WLOOP2_FL: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
         WC=>WALL(IW)
         IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE WLOOP2_FL
             
         II  = WC%II 
         JJ  = WC%JJ
         KK  = WC%KK
         IOR = WC%IOR
         SURF_INDEX = WC%SURF_INDEX
         IIG = WC%IIG
         JJG = WC%JJG
         KKG = WC%KKG
         
         ! overwrite first off-wall advective flux if flow is away from the wall and if the face is not also a wall cell

        IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY .AND. WC%BOUNDARY_TYPE/=OPEN_BOUNDARY) THEN

            OFF_WALL_SELECT_2: SELECT CASE(IOR)
               CASE( 1) OFF_WALL_SELECT_2
                  !      ghost          FX/UU(II+1)
                  ! ///   II   ///  II+1  |  II+2  | ...
                  !                       ^ WALL_INDEX(II+1,+1)
                  IF ((UU(II+1,JJ,KK)>0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II+1,JJ,KK),+1)>0)) THEN
                     ZZZ(1:3) = (/WC%RHO_F,RHOP(II+1:II+2,JJ,KK)/)*(/WC%ZZ_F(N),ZZP(II+1:II+2,JJ,KK,N)/)
                     FX(II+1,JJ,KK) = UU(II+1,JJ,KK)*SCALAR_FACE_VALUE(UU(II+1,JJ,KK),ZZZ,FLUX_LIMITER)*R(II+1)
                     IF (INCLUDE_NUMERICAL_DIFFUSION) DFX(II+1,JJ,KK,N) = FX(II+1,JJ,KK) - 0.5_EB*(ZZZ(2)+ZZZ(3))*UU(II+1,JJ,KK)
                  ENDIF
               CASE(-1) OFF_WALL_SELECT_2
                  !            FX/UU(II-2)     ghost
                  ! ... |  II-2  |  II-1  ///   II   ///
                  !              ^ WALL_INDEX(II-1,-1)
                  IF ((UU(II-2,JJ,KK)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II-1,JJ,KK),-1)>0)) THEN
                     ZZZ(2:4) = (/RHOP(II-2:II-1,JJ,KK),WC%RHO_F/)*(/ZZP(II-2:II-1,JJ,KK,N),WC%ZZ_F(N)/)
                     FX(II-2,JJ,KK) = UU(II-2,JJ,KK)*SCALAR_FACE_VALUE(UU(II-2,JJ,KK),ZZZ,FLUX_LIMITER)*R(II-2)
                     IF (INCLUDE_NUMERICAL_DIFFUSION) DFX(II-2,JJ,KK,N) = FX(II-2,JJ,KK) - 0.5_EB*(ZZZ(2)+ZZZ(3))*UU(II-2,JJ,KK)
                  ENDIF
               CASE( 2) OFF_WALL_SELECT_2
                  IF ((VV(II,JJ+1,KK)>0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ+1,KK),+2)>0)) THEN
                     ZZZ(1:3) = (/WC%RHO_F,RHOP(II,JJ+1:JJ+2,KK)/)*(/WC%ZZ_F(N),ZZP(II,JJ+1:JJ+2,KK,N)/)
                     FY(II,JJ+1,KK) = VV(II,JJ+1,KK)*SCALAR_FACE_VALUE(VV(II,JJ+1,KK),ZZZ,FLUX_LIMITER)
                     IF (INCLUDE_NUMERICAL_DIFFUSION) DFY(II,JJ+1,KK,N) = FY(II,JJ+1,KK) - 0.5_EB*(ZZZ(2)+ZZZ(3))*VV(II,JJ+1,KK)
                  ENDIF
               CASE(-2) OFF_WALL_SELECT_2
                  IF ((VV(II,JJ-2,KK)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ-1,KK),-2)>0)) THEN
                     ZZZ(2:4) = (/RHOP(II,JJ-2:JJ-1,KK),WC%RHO_F/)*(/ZZP(II,JJ-2:JJ-1,KK,N),WC%ZZ_F(N)/)
                     FY(II,JJ-2,KK) = VV(II,JJ-2,KK)*SCALAR_FACE_VALUE(VV(II,JJ-2,KK),ZZZ,FLUX_LIMITER)
                     IF (INCLUDE_NUMERICAL_DIFFUSION) DFY(II,JJ-2,KK,N) = FY(II,JJ-2,KK) - 0.5_EB*(ZZZ(2)+ZZZ(3))*VV(II,JJ-2,KK)
                  ENDIF
               CASE( 3) OFF_WALL_SELECT_2
                  IF ((WW(II,JJ,KK+1)>0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ,KK+1),+3)>0)) THEN
                     ZZZ(1:3) = (/WC%RHO_F,RHOP(II,JJ,KK+1:KK+2)/)*(/WC%ZZ_F(N),ZZP(II,JJ,KK+1:KK+2,N)/)
                     FZ(II,JJ,KK+1) = WW(II,JJ,KK+1)*SCALAR_FACE_VALUE(WW(II,JJ,KK+1),ZZZ,FLUX_LIMITER)
                     IF (INCLUDE_NUMERICAL_DIFFUSION) DFZ(II,JJ,KK+1,N) = FZ(II,JJ,KK+1) - 0.5_EB*(ZZZ(2)+ZZZ(3))*WW(II,JJ,KK+1)
                  ENDIF
               CASE(-3) OFF_WALL_SELECT_2
                  IF ((WW(II,JJ,KK-2)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ,KK-1),-3)>0)) THEN
                     ZZZ(2:4) = (/RHOP(II,JJ,KK-2:KK-1),WC%RHO_F/)*(/ZZP(II,JJ,KK-2:KK-1,N),WC%ZZ_F(N)/)
                     FZ(II,JJ,KK-2) = WW(II,JJ,KK-2)*SCALAR_FACE_VALUE(WW(II,JJ,KK-2),ZZZ,FLUX_LIMITER)
                     IF (INCLUDE_NUMERICAL_DIFFUSION) DFZ(II,JJ,KK-2,N) = FZ(II,JJ,KK-2) - 0.5_EB*(ZZZ(2)+ZZZ(3))*WW(II,JJ,KK-2)
                  ENDIF
            END SELECT OFF_WALL_SELECT_2

         ENDIF

         ! Get the normal components of velocity at the wall

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

         ! At interpolated boundaries, use the actual normal components of velocity, not the ones that have been forced to match.
         ! This line is temporarily commented out because it may not be necessary.
       
     !!  IF (WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY) UN = UVW_SAVE(IW)
         
         ! Special handling of normal component of velocity at surfaces with a specified mass flux.
         
         IF ((SURFACE(SURF_INDEX)%SPECIES_BC_INDEX==SPECIFIED_MASS_FLUX .OR. &
             (SURFACE(SURF_INDEX)%SPECIES_BC_INDEX==HVAC_BOUNDARY       .OR. &
              ANY(SURFACE(SURF_INDEX)%LEAK_PATH>0._EB)) .AND. WC%UWS<0._EB) .AND. WC%ZZ_F(N)>0._EB) THEN
            ! recreate diffusive flux from divg b/c UWP based on old RHODW
            RHO_D_DZDN = 2._EB*WC%RHODW(N)*(ZZP(IIG,JJG,KKG,N)-WC%ZZ_F(N))*WC%RDN
            UN = SIGN(1._EB,REAL(IOR,EB))*(WC%MASSFLUX(N) + RHO_D_DZDN)/(WC%RHO_F*WC%ZZ_F(N))
         ENDIF

         ! Compute species mass flux on the face of the wall cell

         SELECT CASE(IOR)
            CASE( 1)
               FX(II,JJ,KK)   = UN*WC%RHO_F*WC%ZZ_F(N)*R(II)
               IF (INCLUDE_NUMERICAL_DIFFUSION) DFX(II,JJ,KK,N)  = 0._EB
            CASE(-1)
               FX(II-1,JJ,KK) = UN*WC%RHO_F*WC%ZZ_F(N)*R(II-1)
               IF (INCLUDE_NUMERICAL_DIFFUSION) DFX(II-1,JJ,KK,N)= 0._EB
            CASE( 2)
               FY(II,JJ,KK)   = UN*WC%RHO_F*WC%ZZ_F(N)
               IF (INCLUDE_NUMERICAL_DIFFUSION) DFY(II,JJ,KK,N)  = 0._EB
            CASE(-2)
               FY(II,JJ-1,KK) = UN*WC%RHO_F*WC%ZZ_F(N)
               IF (INCLUDE_NUMERICAL_DIFFUSION) DFY(II,JJ-1,KK,N)= 0._EB
            CASE( 3) 
               FZ(II,JJ,KK)   = UN*WC%RHO_F*WC%ZZ_F(N)
               IF (INCLUDE_NUMERICAL_DIFFUSION) DFZ(II,JJ,KK,N)  = 0._EB
            CASE(-3) 
               FZ(II,JJ,KK-1) = UN*WC%RHO_F*WC%ZZ_F(N)
               IF (INCLUDE_NUMERICAL_DIFFUSION) DFZ(II,JJ,KK-1,N)= 0._EB
         END SELECT

      ENDDO WLOOP2_FL
      !$OMP END DO
   
      !$OMP DO COLLAPSE(3) SCHEDULE(STATIC) PRIVATE(K,J,I)
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
               DEL_RHO_D_DEL_Z(I,J,K,N) = -DEL_RHO_D_DEL_Z(I,J,K,N)             &
                                        + (FX(I,J,K)-FX(I-1,J,K))*RDX(I)*RRN(I) &
                                        + (FY(I,J,K)-FY(I,J-1,K))*RDY(J)        &
                                        + (FZ(I,J,K)-FZ(I,J,K-1))*RDZ(K)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO NOWAIT
      !$OMP END PARALLEL
      
ENDDO SPECIES_LOOP
 
TUSED(3,NM)=TUSED(3,NM)+SECOND()-TNOW
END SUBROUTINE MASS_FINITE_DIFFERENCES

 
SUBROUTINE DENSITY(NM)

! Update the density and species mass fractions

USE COMP_FUNCTIONS, ONLY: SECOND 
USE PHYSICAL_FUNCTIONS, ONLY : GET_SPECIFIC_GAS_CONSTANT,GET_SENSIBLE_ENTHALPY
USE GLOBAL_CONSTANTS, ONLY: N_TRACKED_SPECIES,TMPMAX,TMPMIN,EVACUATION_ONLY, &
                            PREDICTOR,CORRECTOR,CHANGE_TIME_STEP,TMPA,N_ZONE, &
                            GAS_SPECIES, R0,SOLID_PHASE_ONLY,TUSED, &
                            DEBUG_OPENMP,CLIP_MASS_FRACTION,ENTHALPY_TRANSPORT
REAL(EB) :: DTRATIO,OMDTRATIO,TNOW,ZZ_GET(0:N_TRACKED_SPECIES),H_S
INTEGER  :: I,J,K,N
INTEGER, INTENT(IN) :: NM
 
IF (EVACUATION_ONLY(NM)) RETURN
IF (SOLID_PHASE_ONLY) RETURN

TNOW=SECOND()
CALL POINT_TO_MESH(NM)

PREDICTOR_STEP: SELECT CASE (PREDICTOR)

CASE(.TRUE.) PREDICTOR_STEP

   !$OMP PARALLEL DEFAULT(NONE) &
   !$OMP SHARED(CHANGE_TIME_STEP,NM,N_TRACKED_SPECIES,KBAR,JBAR,IBAR,SOLID,CELL_INDEX,ZZS,DT,DEL_RHO_D_DEL_Z, &
   !$OMP        RHO,ZZ, &
   !$OMP        DTRATIO,DT_PREV,OMDTRATIO,RHOS, &
   !$OMP        FRHO,CLIP_MASS_FRACTION,N_ZONE,PBAR_S,PBAR,D_PBAR_DT,KBP1,JBP1,IBP1,RSUM,TMP,PRESSURE_ZONE, &
   !$OMP        TMPMIN,TMPMAX)

   IF (.NOT.CHANGE_TIME_STEP(NM)) THEN
   
   ! NOTE: This IF statement is required because the source terms for species are zeroed out at
   !       the beginning of DIVERGENCE_PART_1, but the array also stores the divergence of the advective
   !       flux which is computed once in MASS_FINITE_DIFFERNENCES above, outside the CHANGE_TIME_STEP loop.
   !       DIVERGENCE_PART_1 is inside the loop.  The source terms are then applied to the next substep in
   !       MASS_FINITE_DIFFERENCES.

      ! Store enthalpy for time derivative in divergence

      IF (ENTHALPY_TRANSPORT) THEN
         !$OMP DO COLLAPSE(3) SCHEDULE(STATIC) PRIVATE(K,J,I,ZZ_GET)
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR   
                  IF (N_TRACKED_SPECIES>0) ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(I,J,K,1:N_TRACKED_SPECIES)
                  CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,TMP(I,J,K))
                  RHO_H_S_0(I,J,K) = RHO(I,J,K)*H_S
               ENDDO
            ENDDO
         ENDDO
         !$OMP END DO
      ENDIF

      !$OMP DO COLLAPSE(4) SCHEDULE(DYNAMIC) PRIVATE(N,K,J,I) 
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
      !$OMP END DO NOWAIT

   ELSE

      !$OMP SINGLE
      DTRATIO   = DT/DT_PREV
      OMDTRATIO = 1._EB - DTRATIO
      !$OMP END SINGLE
      !$OMP DO COLLAPSE(4) SCHEDULE(DYNAMIC) PRIVATE(N,K,J,I) 
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
      !$OMP END DO NOWAIT

   ENDIF

   ! Predict the density at the next time step (RHOS or RHO^*)

   !$OMP DO COLLAPSE(3) SCHEDULE(STATIC) PRIVATE(K,J,I)
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            RHOS(I,J,K) = RHO(I,J,K)-DT*FRHO(I,J,K)
         ENDDO
      ENDDO
   ENDDO
   !$OMP END DO

   ! Correct densities above or below clip limits

   !$OMP SINGLE
   CALL CHECK_DENSITY
   !$OMP END SINGLE
   
   ! Extract mass fraction from RHO * ZZ

   !$OMP DO COLLAPSE(4) SCHEDULE(DYNAMIC) PRIVATE(N,K,J,I) 
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
   !$OMP END DO

   ! Correct mass fractions above or below clip limits

   IF (CLIP_MASS_FRACTION) THEN
      !$OMP WORKSHARE
      ZZS(1:IBAR,1:JBAR,1:KBAR,1:N_TRACKED_SPECIES) = MAX(0._EB,MIN(1._EB,ZZS(1:IBAR,1:JBAR,1:KBAR,1:N_TRACKED_SPECIES)))
      !$OMP END WORKSHARE
   ELSE
      !$OMP SINGLE
      CALL CHECK_MASS_FRACTION
      !$OMP END SINGLE
   ENDIF

   ! Predict background pressure at next time step

   !$OMP DO PRIVATE(I)
   DO I=1,N_ZONE
      PBAR_S(:,I) = PBAR(:,I) + D_PBAR_DT(I)*DT
   ENDDO
   !$OMP END DO NOWAIT

   ! Compute molecular weight term RSUM=R0*SUM(Y_i/M_i)

   IF (N_TRACKED_SPECIES>0) THEN
      !$OMP DO COLLAPSE(3) SCHEDULE(STATIC) PRIVATE(K,J,I,ZZ_GET)
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR   
               ZZ_GET(1:N_TRACKED_SPECIES) = ZZS(I,J,K,1:N_TRACKED_SPECIES)
               CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,RSUM(I,J,K))
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO
   ENDIF

   ! Extract predicted temperature at next time step from Equation of State
   
   !$OMP DO COLLAPSE(3) SCHEDULE(STATIC) PRIVATE(K,J,I)
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            TMP(I,J,K) = PBAR_S(K,PRESSURE_ZONE(I,J,K))/(RSUM(I,J,K)*RHOS(I,J,K))
         ENDDO
      ENDDO
   ENDDO
   !$OMP END DO

   !$OMP WORKSHARE
   TMP = MAX(TMPMIN,MIN(TMPMAX,TMP))
   !$OMP END WORKSHARE NOWAIT
   !$OMP END PARALLEL


! The CORRECTOR step
   
CASE(.FALSE.) PREDICTOR_STEP

   ! Store enthalpy for time derivative in divergence

   IF (ENTHALPY_TRANSPORT) THEN
      !$OMP DO COLLAPSE(3) SCHEDULE(STATIC) PRIVATE(K,J,I,ZZ_GET)
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR   
               IF (N_TRACKED_SPECIES>0) ZZ_GET(1:N_TRACKED_SPECIES) = ZZS(I,J,K,1:N_TRACKED_SPECIES)
               CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,TMP(I,J,K))
               RHO_H_S_0(I,J,K) = 0.5_EB*(RHO_H_S_0(I,J,K) + RHOS(I,J,K)*H_S)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO
   ENDIF

   ! Correct species mass fraction at next time step (ZZ here actually means ZZ*RHO)

   !$OMP PARALLEL DEFAULT(NONE) &
   !$OMP SHARED(N_TRACKED_SPECIES,KBAR,JBAR,IBAR,SOLID,CELL_INDEX,ZZ,RHO,RHOS,ZZS,DT,DEL_RHO_D_DEL_Z, &
   !$OMP        FRHO,CLIP_MASS_FRACTION,N_ZONE,PBAR,PBAR_S,D_PBAR_S_DT,KBP1,JBP1,IBP1,RSUM,TMP,PRESSURE_ZONE, &
   !$OMP        TMPMIN,TMPMAX)
   
   !$OMP DO COLLAPSE(4) SCHEDULE(DYNAMIC) PRIVATE(N,K,J,I) 
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
   !$OMP END DO NOWAIT

   ! Correct density at next time step

   !$OMP DO COLLAPSE(3) SCHEDULE(STATIC) PRIVATE(K,J,I)
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            RHO(I,J,K) = .5_EB*(RHO(I,J,K)+RHOS(I,J,K)-DT*FRHO(I,J,K))
         ENDDO
      ENDDO
   ENDDO
   !$OMP END DO

   ! Correct densities above or below clip limits

   !$OMP SINGLE
   CALL CHECK_DENSITY
   !$OMP END SINGLE
   
   ! Extract Y_n from rho*Y_n

   !$OMP DO COLLAPSE(4) SCHEDULE(DYNAMIC) PRIVATE(N,K,J,I) 
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
   !$OMP END DO

   ! Correct mass fractions above or below clip limits

   IF (CLIP_MASS_FRACTION) THEN
      !$OMP WORKSHARE
      ZZ(1:IBAR,1:JBAR,1:KBAR,1:N_TRACKED_SPECIES) = MAX(0._EB,MIN(1._EB,ZZ(1:IBAR,1:JBAR,1:KBAR,1:N_TRACKED_SPECIES)))
      !$OMP END WORKSHARE NOWAIT
   ELSE
      !$OMP SINGLE
      CALL CHECK_MASS_FRACTION
      !$OMP END SINGLE
   ENDIF

   ! Correct background pressure

   !$OMP DO PRIVATE(I)
   DO I=1,N_ZONE
      PBAR(:,I) = .5_EB*(PBAR(:,I) + PBAR_S(:,I) + D_PBAR_S_DT(I)*DT)
   ENDDO
   !$OMP END DO NOWAIT
 
   ! Compute molecular weight term RSUM=R0*SUM(Y_i/M_i)

   IF (N_TRACKED_SPECIES>0) THEN
      !$OMP DO COLLAPSE(3) SCHEDULE(STATIC) PRIVATE(K,J,I,ZZ_GET)
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR   
               ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(I,J,K,1:N_TRACKED_SPECIES)
               CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,RSUM(I,J,K)) 
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO
   ENDIF

   ! Extract predicted temperature at next time step from Equation of State

   !$OMP DO COLLAPSE(3) SCHEDULE(STATIC) PRIVATE(K,J,I)
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            TMP(I,J,K) = PBAR(K,PRESSURE_ZONE(I,J,K))/(RSUM(I,J,K)*RHO(I,J,K))
         ENDDO
      ENDDO
   ENDDO
   !$OMP END DO

   !$OMP WORKSHARE
   TMP = MAX(TMPMIN,MIN(TMPMAX,TMP))
   !$OMP END WORKSHARE NOWAIT
   !$OMP END PARALLEL

END SELECT PREDICTOR_STEP

TUSED(3,NM)=TUSED(3,NM)+SECOND()-TNOW
 
END SUBROUTINE DENSITY
 

SUBROUTINE CHECK_DENSITY
 
! Redistribute mass from cells below or above the density cut-off limits
! Do not apply OpenMP to this routine

USE GLOBAL_CONSTANTS, ONLY : PREDICTOR,RHOMIN,RHOMAX
REAL(EB) :: MASS_N(-3:3),CONST,MASS_C,RHO_CUT,VC(-3:3),SIGN_FACTOR,SUM_MASS_N
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
         VC( 0)  = DX(I)  *DY(J)  *DZ(K)
         VC(-1)  = DX(I-1)*DY(J)  *DZ(K)
         VC( 1)  = DX(I+1)*DY(J)  *DZ(K)
         VC(-2)  = DX(I)  *DY(J-1)*DZ(K)
         VC( 2)  = DX(I)  *DY(J+1)*DZ(K)
         VC(-3)  = DX(I)  *DY(J)  *DZ(K-1)
         VC( 3)  = DX(I)  *DY(J)  *DZ(K+1)

         MASS_C = ABS(RHO_CUT-RHOP(I,J,K))*VC(0)
         IF (WALL_INDEX(IC,-1)==0) MASS_N(-1) = ABS(MIN(RHOMAX,MAX(RHOMIN,RHOP(I-1,J,K)))-RHO_CUT)*VC(-1)
         IF (WALL_INDEX(IC, 1)==0) MASS_N( 1) = ABS(MIN(RHOMAX,MAX(RHOMIN,RHOP(I+1,J,K)))-RHO_CUT)*VC( 1)
         IF (WALL_INDEX(IC,-2)==0) MASS_N(-2) = ABS(MIN(RHOMAX,MAX(RHOMIN,RHOP(I,J-1,K)))-RHO_CUT)*VC(-2)
         IF (WALL_INDEX(IC, 2)==0) MASS_N( 2) = ABS(MIN(RHOMAX,MAX(RHOMIN,RHOP(I,J+1,K)))-RHO_CUT)*VC( 2)
         IF (WALL_INDEX(IC,-3)==0) MASS_N(-3) = ABS(MIN(RHOMAX,MAX(RHOMIN,RHOP(I,J,K-1)))-RHO_CUT)*VC(-3)
         IF (WALL_INDEX(IC, 3)==0) MASS_N( 3) = ABS(MIN(RHOMAX,MAX(RHOMIN,RHOP(I,J,K+1)))-RHO_CUT)*VC( 3)
         SUM_MASS_N = SUM(MASS_N)
         IF (SUM_MASS_N<=ZERO_P) CYCLE
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

USE GLOBAL_CONSTANTS, ONLY : PREDICTOR,N_TRACKED_SPECIES
REAL(EB) :: SUM,CONST,MASS_C,MASS_N(-3:3),ZZ_CUT,SIGN_FACTOR,SUM_MASS_N,VC(-3:3)
INTEGER  :: IC,N,I,J,K
REAL(EB), POINTER, DIMENSION(:,:,:) :: DELTA_ZZ=>NULL()

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
            VC( 0) = DX(I)  *DY(J)  *DZ(K)
            VC(-1) = DX(I-1)*DY(J)  *DZ(K)
            VC( 1) = DX(I+1)*DY(J)  *DZ(K)
            VC(-2) = DX(I)  *DY(J-1)*DZ(K)
            VC( 2) = DX(I)  *DY(J+1)*DZ(K)
            VC(-3) = DX(I)  *DY(J)  *DZ(K-1)
            VC( 3) = DX(I)  *DY(J)  *DZ(K+1)

            IF (WALL_INDEX(IC,-1)==0) MASS_N(-1) = RHOP(I-1,J,K)*ABS(MIN(1._EB,MAX(0._EB,ZZP(I-1,J,K,N)))-ZZ_CUT)*VC(-1)
            IF (WALL_INDEX(IC, 1)==0) MASS_N( 1) = RHOP(I+1,J,K)*ABS(MIN(1._EB,MAX(0._EB,ZZP(I+1,J,K,N)))-ZZ_CUT)*VC( 1)
            IF (WALL_INDEX(IC,-2)==0) MASS_N(-2) = RHOP(I,J-1,K)*ABS(MIN(1._EB,MAX(0._EB,ZZP(I,J-1,K,N)))-ZZ_CUT)*VC(-2)
            IF (WALL_INDEX(IC, 2)==0) MASS_N( 2) = RHOP(I,J+1,K)*ABS(MIN(1._EB,MAX(0._EB,ZZP(I,J+1,K,N)))-ZZ_CUT)*VC( 2)
            IF (WALL_INDEX(IC,-3)==0) MASS_N(-3) = RHOP(I,J,K-1)*ABS(MIN(1._EB,MAX(0._EB,ZZP(I,J,K-1,N)))-ZZ_CUT)*VC(-3)
            IF (WALL_INDEX(IC, 3)==0) MASS_N( 3) = RHOP(I,J,K+1)*ABS(MIN(1._EB,MAX(0._EB,ZZP(I,J,K+1,N)))-ZZ_CUT)*VC( 3)

            SUM_MASS_N = SUM(MASS_N)
            IF (SUM_MASS_N<=ZERO_P) CYCLE
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

ENDDO SPECIES_LOOP

END SUBROUTINE CHECK_MASS_FRACTION



REAL(EB) FUNCTION SCALAR_FACE_VALUE(A,U,LIMITER)

REAL(EB), INTENT(IN) :: A,U(4)
INTEGER, INTENT(IN) :: LIMITER
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
         IF (ABS(DU_LOC)>ZERO_P) R = DU_UP/DU_LOC
         B = MAX(0._EB,MIN(2._EB*R,1._EB),MIN(R,2._EB))
         SCALAR_FACE_VALUE = U(2) + 0.5_EB*B*DU_LOC
      CASE(3) ! MINMOD
         IF (ABS(DU_LOC)>ZERO_P) R = DU_UP/DU_LOC
         B = MAX(0._EB,MIN(1._EB,R))
         SCALAR_FACE_VALUE = U(2) + 0.5_EB*B*DU_LOC
      CASE(4) ! CHARM
         IF (ABS(DU_UP)>ZERO_P) R = DU_LOC/DU_UP
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
         IF (ABS(DU_LOC)>ZERO_P) R = DU_UP/DU_LOC
         B = MAX(0._EB,MIN(2._EB*R,1._EB),MIN(R,2._EB))
         SCALAR_FACE_VALUE = U(3) - 0.5_EB*B*DU_LOC
      CASE(3) ! MINMOD
         IF (ABS(DU_LOC)>ZERO_P) R = DU_UP/DU_LOC
         B = MAX(0._EB,MIN(1._EB,R))
         SCALAR_FACE_VALUE = U(3) - 0.5_EB*B*DU_LOC
      CASE(4) ! CHARM
         IF (ABS(DU_UP)>ZERO_P) R = DU_LOC/DU_UP
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
