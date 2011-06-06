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
USE PHYSICAL_FUNCTIONS, ONLY: GET_AVERAGE_SPECIFIC_HEAT
USE GLOBAL_CONSTANTS, ONLY: N_TRACKED_SPECIES,NULL_BOUNDARY,POROUS_BOUNDARY,OPEN_BOUNDARY,INTERPOLATED_BOUNDARY, &
                            PREDICTOR,CORRECTOR,EVACUATION_ONLY,SOLID_PHASE_ONLY,TUSED,DEBUG_OPENMP,SOLID_BOUNDARY, &
                            NO_MASS_FLUX,SPECIFIED_MASS_FLUX,HVAC_BOUNDARY,ENTHALPY_TRANSPORT
INTEGER, INTENT(IN) :: NM
REAL(EB) :: TNOW,ZZZ(4),UN,CP,ZZ_GET(0:N_TRACKED_SPECIES),E_F
REAL(EB) :: RHO_D_DZDN
INTEGER  :: I,J,K,N,II,JJ,KK,IIG,JJG,KKG,IW,IOR,IBC
REAL(EB), POINTER, DIMENSION(:) :: UWP
REAL(EB), POINTER, DIMENSION(:,:,:) :: FX=>NULL(),FY=>NULL(),FZ=>NULL(),EE=>NULL()
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZP=>NULL()

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
   IF (N_TRACKED_SPECIES > 0) ZZP => ZZ
   EE => E
ELSE
   UU => US
   VV => VS
   WW => WS
   DP => DS
   RHOP => RHOS
   UWP  => UWS
   IF (N_TRACKED_SPECIES > 0) ZZP => ZZS
   EE => ES
ENDIF

FX=>WORK4
FY=>WORK5
FZ=>WORK6
FX=0._EB
FY=0._EB
FZ=0._EB

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(KBAR,JBAR,IBAR,KBM1,JBM1,IBM1,RHOP,FX,FY,FZ,UU,VV,WW,FLUX_LIMITER,R, &
!$OMP        N_EXTERNAL_WALL_CELLS,N_INTERNAL_WALL_CELLS,BOUNDARY_TYPE,IJKW,WALL_INDEX,CELL_INDEX, &
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

!$OMP DO PRIVATE(IW,II,JJ,KK,IOR,IBC,IIG,JJG,KKG,ZZZ,UN)
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
            ZZZ(1:3) = (/RHO_F(IW),RHOP(II+1:II+2,JJ,KK)/)
            FX(II+1,JJ,KK) = UU(II+1,JJ,KK)*SCALAR_FACE_VALUE(UU(II+1,JJ,KK),ZZZ,FLUX_LIMITER)*R(II+1)
         ENDIF
      CASE(-1) OFF_WALL_SELECT_1
         !            FX/UU(II-2)     ghost
         ! ... |  II-2  |  II-1  ///   II   ///
         !              ^ WALL_INDEX(II-1,-1)
         IF ((UU(II-2,JJ,KK)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II-1,JJ,KK),-1)>0)) THEN
            ZZZ(2:4) = (/RHOP(II-2:II-1,JJ,KK),RHO_F(IW)/)
            FX(II-2,JJ,KK) = UU(II-2,JJ,KK)*SCALAR_FACE_VALUE(UU(II-2,JJ,KK),ZZZ,FLUX_LIMITER)*R(II-2)
         ENDIF
      CASE( 2) OFF_WALL_SELECT_1
         IF ((VV(II,JJ+1,KK)>0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ+1,KK),+2)>0)) THEN
            ZZZ(1:3) = (/RHO_F(IW),RHOP(II,JJ+1:JJ+2,KK)/)
            FY(II,JJ+1,KK) = VV(II,JJ+1,KK)*SCALAR_FACE_VALUE(VV(II,JJ+1,KK),ZZZ,FLUX_LIMITER)
         ENDIF
      CASE(-2) OFF_WALL_SELECT_1
         IF ((VV(II,JJ-2,KK)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ-1,KK),-2)>0)) THEN
            ZZZ(2:4) = (/RHOP(II,JJ-2:JJ-1,KK),RHO_F(IW)/)
            FY(II,JJ-2,KK) = VV(II,JJ-2,KK)*SCALAR_FACE_VALUE(VV(II,JJ-2,KK),ZZZ,FLUX_LIMITER)
         ENDIF
      CASE( 3) OFF_WALL_SELECT_1
         IF ((WW(II,JJ,KK+1)>0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ,KK+1),+3)>0)) THEN
            ZZZ(1:3) = (/RHO_F(IW),RHOP(II,JJ,KK+1:KK+2)/)
            FZ(II,JJ,KK+1) = WW(II,JJ,KK+1)*SCALAR_FACE_VALUE(WW(II,JJ,KK+1),ZZZ,FLUX_LIMITER)
         ENDIF
      CASE(-3) OFF_WALL_SELECT_1
         IF ((WW(II,JJ,KK-2)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ,KK-1),-3)>0)) THEN
            ZZZ(2:4) = (/RHOP(II,JJ,KK-2:KK-1),RHO_F(IW)/)
            FZ(II,JJ,KK-2) = WW(II,JJ,KK-2)*SCALAR_FACE_VALUE(WW(II,JJ,KK-2),ZZZ,FLUX_LIMITER)
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

   ! In case of interpolated boundary, use the original velocity, not the averaged value

   IF (BOUNDARY_TYPE(IW)==INTERPOLATED_BOUNDARY) UN = UVW_SAVE(IW)
   
   ! Compute flux on the face of the wall cell

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

!Pointer settings now at the beginning of this subroutine 
!IF (N_TRACKED_SPECIES > 0) THEN
!   IF (PREDICTOR) ZZP => ZZ
!   IF (CORRECTOR) ZZP => ZZS
!ENDIF
 
SPECIES_LOOP: DO N=1,N_TRACKED_SPECIES

      FX=0._EB
      FY=0._EB
      FZ=0._EB
   
      !$OMP PARALLEL DEFAULT(NONE) &
      !$OMP SHARED(N,KBAR,JBAR,IBAR,KBM1,JBM1,IBM1,RHOP,ZZP,FX,FY,FZ,UU,VV,WW,FLUX_LIMITER,R, &
      !$OMP        N_EXTERNAL_WALL_CELLS,N_INTERNAL_WALL_CELLS,BOUNDARY_TYPE,IJKW, &
      !$OMP        WALL_INDEX,CELL_INDEX,RHO_F,ZZ_F,UVW_SAVE, &
      !$OMP        SURFACE,UWS,RHODW,RDN,MASSFLUX, &
      !$OMP        SOLID,DEL_RHO_D_DEL_Z,RDX,RDY,RDZ,RRN)

      !$OMP DO COLLAPSE(3) SCHEDULE(STATIC) PRIVATE(K,J,I,ZZZ)
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBM1
               ZZZ(1:4) = RHOP(I-1:I+2,J,K)*ZZP(I-1:I+2,J,K,N)
               FX(I,J,K) = UU(I,J,K)*SCALAR_FACE_VALUE(UU(I,J,K),ZZZ,FLUX_LIMITER)*R(I)
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
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO

     !$OMP DO SCHEDULE(STATIC) &
     !$OMP PRIVATE(IW,II,JJ,KK,IOR,IBC,IIG,JJG,KKG,ZZZ,UN,RHO_D_DZDN)
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
                  ZZZ(1:3) = (/RHO_F(IW),RHOP(II+1:II+2,JJ,KK)/)*(/ZZ_F(IW,N),ZZP(II+1:II+2,JJ,KK,N)/)
                  FX(II+1,JJ,KK) = UU(II+1,JJ,KK)*SCALAR_FACE_VALUE(UU(II+1,JJ,KK),ZZZ,FLUX_LIMITER)*R(II+1)
               ENDIF
            CASE(-1) OFF_WALL_SELECT_2
               !            FX/UU(II-2)     ghost
               ! ... |  II-2  |  II-1  ///   II   ///
               !              ^ WALL_INDEX(II-1,-1)
               IF ((UU(II-2,JJ,KK)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II-1,JJ,KK),-1)>0)) THEN
                  ZZZ(2:4) = (/RHOP(II-2:II-1,JJ,KK),RHO_F(IW)/)*(/ZZP(II-2:II-1,JJ,KK,N),ZZ_F(IW,N)/)
                  FX(II-2,JJ,KK) = UU(II-2,JJ,KK)*SCALAR_FACE_VALUE(UU(II-2,JJ,KK),ZZZ,FLUX_LIMITER)*R(II-2)
               ENDIF
            CASE( 2) OFF_WALL_SELECT_2
               IF ((VV(II,JJ+1,KK)>0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ+1,KK),+2)>0)) THEN
                  ZZZ(1:3) = (/RHO_F(IW),RHOP(II,JJ+1:JJ+2,KK)/)*(/ZZ_F(IW,N),ZZP(II,JJ+1:JJ+2,KK,N)/)
                  FY(II,JJ+1,KK) = VV(II,JJ+1,KK)*SCALAR_FACE_VALUE(VV(II,JJ+1,KK),ZZZ,FLUX_LIMITER)
               ENDIF
            CASE(-2) OFF_WALL_SELECT_2
               IF ((VV(II,JJ-2,KK)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ-1,KK),-2)>0)) THEN
                  ZZZ(2:4) = (/RHOP(II,JJ-2:JJ-1,KK),RHO_F(IW)/)*(/ZZP(II,JJ-2:JJ-1,KK,N),ZZ_F(IW,N)/)
                  FY(II,JJ-2,KK) = VV(II,JJ-2,KK)*SCALAR_FACE_VALUE(VV(II,JJ-2,KK),ZZZ,FLUX_LIMITER)
               ENDIF
            CASE( 3) OFF_WALL_SELECT_2
               IF ((WW(II,JJ,KK+1)>0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ,KK+1),+3)>0)) THEN
                  ZZZ(1:3) = (/RHO_F(IW),RHOP(II,JJ,KK+1:KK+2)/)*(/ZZ_F(IW,N),ZZP(II,JJ,KK+1:KK+2,N)/)
                  FZ(II,JJ,KK+1) = WW(II,JJ,KK+1)*SCALAR_FACE_VALUE(WW(II,JJ,KK+1),ZZZ,FLUX_LIMITER)
               ENDIF
            CASE(-3) OFF_WALL_SELECT_2
               IF ((WW(II,JJ,KK-2)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ,KK-1),-3)>0)) THEN
                  ZZZ(2:4) = (/RHOP(II,JJ,KK-2:KK-1),RHO_F(IW)/)*(/ZZP(II,JJ,KK-2:KK-1,N),ZZ_F(IW,N)/)
                  FZ(II,JJ,KK-2) = WW(II,JJ,KK-2)*SCALAR_FACE_VALUE(WW(II,JJ,KK-2),ZZZ,FLUX_LIMITER)
               ENDIF
         END SELECT OFF_WALL_SELECT_2

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

         IF (BOUNDARY_TYPE(IW)==INTERPOLATED_BOUNDARY) UN = UVW_SAVE(IW)

         ! At forced flow boundaries, use the specified normal component of velocity

         IF ((SURFACE(IBC)%SPECIES_BC_INDEX==SPECIFIED_MASS_FLUX .OR. &
             (SURFACE(IBC)%SPECIES_BC_INDEX==HVAC_BOUNDARY       .OR. &
              ANY(SURFACE(IBC)%LEAK_PATH>0._EB)) .AND. UWS(IW)<0._EB) .AND. ZZ_F(IW,N)>ZERO_P) THEN
            ! recreate diffusive flux from divg b/c UWP based on old RHODW
            RHO_D_DZDN = 2._EB*RHODW(IW,N)*(ZZP(IIG,JJG,KKG,N)-ZZ_F(IW,N))*RDN(IW)
            UN = SIGN(1._EB,REAL(IOR,EB))*(MASSFLUX(IW,N) + RHO_D_DZDN)/(RHO_F(IW)*ZZ_F(IW,N))
 !!         UN = -SIGN(1._EB,REAL(IOR,EB))*UWP(IW)
         ENDIF

         ! Compute species mass flux on the face of the wall cell

         SELECT CASE(IOR)
            CASE( 1)
               FX(II,JJ,KK)   = UN*RHO_F(IW)*ZZ_F(IW,N)*R(II)
            CASE(-1)
               FX(II-1,JJ,KK) = UN*RHO_F(IW)*ZZ_F(IW,N)*R(II-1)
            CASE( 2)
               FY(II,JJ,KK)   = UN*RHO_F(IW)*ZZ_F(IW,N)
            CASE(-2)
               FY(II,JJ-1,KK) = UN*RHO_F(IW)*ZZ_F(IW,N)
            CASE( 3) 
               FZ(II,JJ,KK)   = UN*RHO_F(IW)*ZZ_F(IW,N)
            CASE(-3) 
               FZ(II,JJ,KK-1) = UN*RHO_F(IW)*ZZ_F(IW,N)
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


! experimental conservative enthalpy transport

ENTHALPY_IF: IF (ENTHALPY_TRANSPORT) THEN

   FX=0._EB
   FY=0._EB
   FZ=0._EB

   !$OMP PARALLEL DEFAULT(NONE) &
   !$OMP SHARED(KBAR,JBAR,IBAR,KBM1,JBM1,IBM1,EE,FX,FY,FZ,UU,VV,WW,FLUX_LIMITER,R, &
   !$OMP        N_EXTERNAL_WALL_CELLS,N_INTERNAL_WALL_CELLS,BOUNDARY_TYPE,IJKW,WALL_INDEX,CELL_INDEX, &
   !$OMP        RHO_F,UVW_SAVE, &
   !$OMP        ENTHALPY_SOURCE,SOLID,RDX,RDY,RDZ,RRN, &
   !$OMP        N_TRACKED_SPECIES,ZZP,ZZ_F)

   !$OMP DO COLLAPSE(3) SCHEDULE(STATIC) PRIVATE(K,J,I,ZZZ)
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBM1
            ZZZ(1:4) = EE(I-1:I+2,J,K)
            FX(I,J,K) = UU(I,J,K)*SCALAR_FACE_VALUE(UU(I,J,K),ZZZ,FLUX_LIMITER)*R(I)
         ENDDO
      ENDDO
   ENDDO
   !$OMP END DO NOWAIT

   !$OMP DO COLLAPSE(3) SCHEDULE(STATIC) PRIVATE(K,J,I,ZZZ)
   DO K=1,KBAR
      DO J=1,JBM1
         DO I=1,IBAR
            ZZZ(1:4) = EE(I,J-1:J+2,K)
            FY(I,J,K) = VV(I,J,K)*SCALAR_FACE_VALUE(VV(I,J,K),ZZZ,FLUX_LIMITER)
         ENDDO
      ENDDO
   ENDDO
   !$OMP END DO NOWAIT

   !$OMP DO COLLAPSE(3) SCHEDULE(STATIC) PRIVATE(K,J,I,ZZZ)
   DO K=1,KBM1
      DO J=1,JBAR
         DO I=1,IBAR
            ZZZ(1:4) = EE(I,J,K-1:K+2)
            FZ(I,J,K) = WW(I,J,K)*SCALAR_FACE_VALUE(WW(I,J,K),ZZZ,FLUX_LIMITER)
         ENDDO
      ENDDO
   ENDDO
   !$OMP END DO

   !$OMP DO PRIVATE(IW,II,JJ,KK,IOR,IBC,IIG,JJG,KKG,ZZZ,UN)
   WLOOP3_FL: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
   
      IF (BOUNDARY_TYPE(IW)==NULL_BOUNDARY .OR. BOUNDARY_TYPE(IW)==POROUS_BOUNDARY) CYCLE WLOOP3_FL
       
      II  = IJKW(1,IW) 
      JJ  = IJKW(2,IW)
      KK  = IJKW(3,IW)
      IOR = IJKW(4,IW)
      IBC = IJKW(5,IW)
      IIG = IJKW(6,IW)
      JJG = IJKW(7,IW)
      KKG = IJKW(8,IW)
      
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

      IF (BOUNDARY_TYPE(IW)==INTERPOLATED_BOUNDARY) UN = UVW_SAVE(IW)
      
      !! At forced flow boundaries, use the specified normal component of velocity
      !
      !IF (SURFACE(IBC)%SPECIES_BC_INDEX==SPECIFIED_MASS_FLUX .OR. &
      !    SURFACE(IBC)%SPECIES_BC_INDEX==HVAC_BOUNDARY) THEN
      !   UN = -SIGN(1._EB,REAL(IOR,EB))*UWP(IW)
      !ENDIF
      
      ! boundary value of enthalpy
      
      IF (UN*REAL(IOR,EB)>ZERO_P) THEN
         ! inflow
         IF (N_TRACKED_SPECIES>0) ZZ_GET(1:N_TRACKED_SPECIES) = ZZ_F(IW,1:N_TRACKED_SPECIES)
         CALL GET_AVERAGE_SPECIFIC_HEAT(ZZ_GET,CP,TMP_F(IW)) 
         E_F = RHO_F(IW)*CP*TMP_F(IW)
      ELSE
         ! outflow
         E_F = EE(IIG,JJG,KKG)
      ENDIF
   
      ! overwrite first off-wall advective flux if flow is away from the wall and if the face is not also a wall cell
      
      OFF_WALL_SELECT_3: SELECT CASE(IOR)
         CASE( 1) OFF_WALL_SELECT_3
            !      ghost          FX/UU(II+1)
            ! ///   II   ///  II+1  |  II+2  | ...
            !                       ^ WALL_INDEX(II+1,+1)
            IF ((UU(II+1,JJ,KK)>0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II+1,JJ,KK),+1)>0)) THEN
               ZZZ(1:3) = (/E_F,EE(II+1:II+2,JJ,KK)/)
               FX(II+1,JJ,KK) = UU(II+1,JJ,KK)*SCALAR_FACE_VALUE(UU(II+1,JJ,KK),ZZZ,FLUX_LIMITER)*R(II+1)
            ENDIF
         CASE(-1) OFF_WALL_SELECT_3
            !            FX/UU(II-2)     ghost
            ! ... |  II-2  |  II-1  ///   II   ///
            !              ^ WALL_INDEX(II-1,-1)
            IF ((UU(II-2,JJ,KK)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II-1,JJ,KK),-1)>0)) THEN
               ZZZ(2:4) = (/EE(II-2:II-1,JJ,KK),E_F/)
               FX(II-2,JJ,KK) = UU(II-2,JJ,KK)*SCALAR_FACE_VALUE(UU(II-2,JJ,KK),ZZZ,FLUX_LIMITER)*R(II-2)
            ENDIF
         CASE( 2) OFF_WALL_SELECT_3
            IF ((VV(II,JJ+1,KK)>0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ+1,KK),+2)>0)) THEN
               ZZZ(1:3) = (/E_F,EE(II,JJ+1:JJ+2,KK)/)
               FY(II,JJ+1,KK) = VV(II,JJ+1,KK)*SCALAR_FACE_VALUE(VV(II,JJ+1,KK),ZZZ,FLUX_LIMITER)
            ENDIF
         CASE(-2) OFF_WALL_SELECT_3
            IF ((VV(II,JJ-2,KK)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ-1,KK),-2)>0)) THEN
               ZZZ(2:4) = (/EE(II,JJ-2:JJ-1,KK),E_F/)
               FY(II,JJ-2,KK) = VV(II,JJ-2,KK)*SCALAR_FACE_VALUE(VV(II,JJ-2,KK),ZZZ,FLUX_LIMITER)
            ENDIF
         CASE( 3) OFF_WALL_SELECT_3
            IF ((WW(II,JJ,KK+1)>0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ,KK+1),+3)>0)) THEN
               ZZZ(1:3) = (/E_F,EE(II,JJ,KK+1:KK+2)/)
               FZ(II,JJ,KK+1) = WW(II,JJ,KK+1)*SCALAR_FACE_VALUE(WW(II,JJ,KK+1),ZZZ,FLUX_LIMITER)
            ENDIF
         CASE(-3) OFF_WALL_SELECT_3
            IF ((WW(II,JJ,KK-2)<0._EB) .AND. .NOT.(WALL_INDEX(CELL_INDEX(II,JJ,KK-1),-3)>0)) THEN
               ZZZ(2:4) = (/EE(II,JJ,KK-2:KK-1),E_F/)
               FZ(II,JJ,KK-2) = WW(II,JJ,KK-2)*SCALAR_FACE_VALUE(WW(II,JJ,KK-2),ZZZ,FLUX_LIMITER)
            ENDIF
      END SELECT OFF_WALL_SELECT_3

      ! Compute flux on the face of the wall cell

      SELECT CASE(IOR)
         CASE( 1)
            FX(II,JJ,KK)   = UN*E_F*R(II)
         CASE(-1)
            FX(II-1,JJ,KK) = UN*E_F*R(II-1)
         CASE( 2)
            FY(II,JJ,KK)   = UN*E_F
         CASE(-2)
            FY(II,JJ-1,KK) = UN*E_F
         CASE( 3)
            FZ(II,JJ,KK)   = UN*E_F
         CASE(-3)
            FZ(II,JJ,KK-1) = UN*E_F
      END SELECT
      
   ENDDO WLOOP3_FL
   !$OMP END DO NOWAIT

   !$OMP DO COLLAPSE(3) SCHEDULE(STATIC) PRIVATE(K,J,I)
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
            ENTHALPY_SOURCE(I,J,K) = -ENTHALPY_SOURCE(I,J,K) &
                                   + (FX(I,J,K)-FX(I-1,J,K))*RDX(I)*RRN(I) &
                                   + (FY(I,J,K)-FY(I,J-1,K))*RDY(J)        &
                                   + (FZ(I,J,K)-FZ(I,J,K-1))*RDZ(K) 
         ENDDO
      ENDDO
   ENDDO
   !$OMP END DO NOWAIT
   !$OMP END PARALLEL

ENDIF ENTHALPY_IF
 
TUSED(3,NM)=TUSED(3,NM)+SECOND()-TNOW
END SUBROUTINE MASS_FINITE_DIFFERENCES

 
SUBROUTINE DENSITY(NM)

! Update the density and species mass fractions

USE COMP_FUNCTIONS, ONLY: SECOND 
USE PHYSICAL_FUNCTIONS, ONLY : GET_SPECIFIC_GAS_CONSTANT
USE GLOBAL_CONSTANTS, ONLY: N_TRACKED_SPECIES,TMPMAX,TMPMIN,EVACUATION_ONLY, &
                            PREDICTOR,CORRECTOR,CHANGE_TIME_STEP,TMPA,N_ZONE, &
                            GAS_SPECIES, R0,SOLID_PHASE_ONLY,TUSED, &
                            DEBUG_OPENMP,CLIP_MASS_FRACTION,ENTHALPY_TRANSPORT
REAL(EB) :: DTRATIO,OMDTRATIO,TNOW,ZZ_GET(0:N_TRACKED_SPECIES)
INTEGER  :: I,J,K,N
INTEGER, INTENT(IN) :: NM
 
IF (EVACUATION_ONLY(NM)) RETURN
IF (SOLID_PHASE_ONLY) RETURN

TNOW=SECOND()
CALL POINT_TO_MESH(NM)

PREDICTOR_STEP: SELECT CASE (PREDICTOR)

CASE(.TRUE.) PREDICTOR_STEP

   !$OMP PARALLEL DEFAULT(NONE) &
   !$OMP SHARED(CHANGE_TIME_STEP,NM,N_TRACKED_SPECIES,KBAR,JBAR,IBAR,SOLID,CELL_INDEX,ZZS,RHO,ZZ,DT, &
   !$OMP        DEL_RHO_D_DEL_Z,DTRATIO,DT_PREV,OMDTRATIO,RHOS, &
   !$OMP        FRHO,CLIP_MASS_FRACTION,N_ZONE,PBAR_S,PBAR,D_PBAR_DT,KBP1,JBP1,IBP1,RSUM,TMP,PRESSURE_ZONE, &
   !$OMP        TMPMIN,TMPMAX)

   IF (.NOT.CHANGE_TIME_STEP(NM)) THEN
   
   ! NOTE: This IF statement is required because the source terms for species and enthalpy are zeroed out at
   !       the beginning of DIVERGENCE_PART_1, but the array also stores the divergence of the advective
   !       flux which computed once in MASS_FINITE_DIFFERNENCES above (outside) the CHANGE_TIME_STEP loop.
   !       DIVERGENCE_PART_1 is inside the loop.  The source terms are then applied to the next substep in
   !       MASS_FINITE_DIFFERENCES.

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
      
      ENTHALPY_IF_1A: IF (ENTHALPY_TRANSPORT) THEN
         !$OMP DO COLLAPSE(3) SCHEDULE(STATIC) PRIVATE(K,J,I)
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
                  ES(I,J,K) = E(I,J,K)-DT*ENTHALPY_SOURCE(I,J,K)
               ENDDO
            ENDDO
         ENDDO
         !$OMP END DO
      ENDIF ENTHALPY_IF_1A

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
      
      ENTHALPY_IF_1B: IF (ENTHALPY_TRANSPORT) THEN
         !$OMP DO COLLAPSE(3) SCHEDULE(STATIC) PRIVATE(K,J,I)
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
                  ES(I,J,K) = OMDTRATIO*E(I,J,K) + DTRATIO*ES(I,J,K)
               ENDDO
            ENDDO
         ENDDO
         !$OMP END DO
      ENDIF ENTHALPY_IF_1B

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
   DO K=0,KBP1
      DO J=0,JBP1
         DO I=0,IBP1
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
               ZZ(I,J,K,N) = .5_EB*(RHO(I,J,K)*ZZ(I,J,K,N) + RHOS(I,J,K)*ZZS(I,J,K,N) - DT*DEL_RHO_D_DEL_Z(I,J,K,N) ) 
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   !$OMP END DO NOWAIT
   
   ENTHALPY_IF_2: IF (ENTHALPY_TRANSPORT) THEN
      !$OMP DO COLLAPSE(3) SCHEDULE(STATIC) PRIVATE(K,J,I)
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               E(I,J,K) = .5_EB*(E(I,J,K)+ES(I,J,K)-DT*ENTHALPY_SOURCE(I,J,K))
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO
   ENDIF ENTHALPY_IF_2

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
   DO K=0,KBP1
      DO J=0,JBP1
         DO I=0,IBP1
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

USE GLOBAL_CONSTANTS, ONLY : PREDICTOR, CORRECTOR, N_TRACKED_SPECIES,RHOMIN,RHOMAX,DEBUG_OPENMP,TWO_D
REAL(EB) :: SUM,CONST,CONST2,RHOMI,RHOPI,RHOMJ,RHOPJ,RHOMK,RHOPK,RHO00,RMIN,RMAX
INTEGER  :: IC,ISUM,I,J,K
LOGICAL :: LC(-3:3)
REAL(EB), POINTER, DIMENSION(:,:,:) :: RHODELTA=>NULL(),V_CELL=>NULL()

RHODELTA => WORK2

IF (PREDICTOR) THEN
   RHOP=>RHOS
ELSE
   RHOP=>RHO
ENDIF

V_CELL => WORK3

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(TWO_D,KBAR,IBAR,JBAR,V_CELL,DX,DY,DZ, &
!$OMP        RHODELTA,CELL_INDEX,SOLID,RHOMIN,RHOP,WALL_INDEX,RHOMAX)

IF (TWO_D) THEN
   !$OMP DO COLLAPSE(3) SCHEDULE(STATIC) PRIVATE(K,J,I)
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            V_CELL(I,J,K) = DX(I)*DZ(K)
         ENDDO
      ENDDO
   ENDDO
   !$OMP END DO NOWAIT
ELSE
   !$OMP DO COLLAPSE(3) SCHEDULE(STATIC) PRIVATE(K,J,I)
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            V_CELL(I,J,K) = DX(I)*DY(J)*DZ(K)
         ENDDO
      ENDDO
   ENDDO
   !$OMP END DO NOWAIT
ENDIF

 
! Correct undershoots

!$OMP WORKSHARE
RHODELTA = 0._EB
!$OMP END WORKSHARE

!$OMP DO COLLAPSE(3) SCHEDULE(DYNAMIC) &
!$OMP PRIVATE(K,J,I,IC,RMIN,SUM,ISUM,LC,RHO00,RHOMI,RHOPI,RHOMJ,RHOPJ,RHOMK,RHOPK,CONST,CONST2)
DO K=1,KBAR
   DO J=1,JBAR
      CHECK_LOOP: DO I=1,IBAR
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
!$OMP END DO

!$OMP WORKSHARE

RHOP(1:IBAR,1:JBAR,1:KBAR) = MAX(RHOMIN,RHOP(1:IBAR,1:JBAR,1:KBAR)+RHODELTA(1:IBAR,1:JBAR,1:KBAR))

! Correct overshoots

RHODELTA = 0._EB
!$OMP END WORKSHARE

!$OMP DO COLLAPSE(3) SCHEDULE(DYNAMIC) &
!$OMP PRIVATE(K,J,I,IC,RMAX,SUM,ISUM,LC,RHO00,RHOMI,RHOPI,RHOMJ,RHOPJ,RHOMK,RHOPK,CONST,CONST2)
DO K=1,KBAR
   DO J=1,JBAR
      CHECK_LOOP2: DO I=1,IBAR
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
!$OMP END DO

!$OMP WORKSHARE
RHOP(1:IBAR,1:JBAR,1:KBAR) = MIN(RHOMAX,RHOP(1:IBAR,1:JBAR,1:KBAR)+RHODELTA(1:IBAR,1:JBAR,1:KBAR))
!$OMP END WORKSHARE NOWAIT
!$OMP END PARALLEL
END SUBROUTINE CHECK_DENSITY
 
 
SUBROUTINE CHECK_MASS_FRACTION

! Redistribute species mass from cells below or above the cut-off limits

USE GLOBAL_CONSTANTS, ONLY : PREDICTOR, CORRECTOR, N_TRACKED_SPECIES,POROUS_BOUNDARY,DEBUG_OPENMP
REAL(EB) :: SUM,CONST,RHZMI,RHZPI,RHZMJ,RHZPJ,RHZMK,RHZPK,RHY0,ZMI,ZPI,ZMJ,ZPJ,ZMK,ZPK,Y00,ZMIN,ZMAX
INTEGER  :: IC,N,ISUM, IW_A(-3:3),I,J,K
LOGICAL  :: LC(-3:3)
REAL(EB), POINTER, DIMENSION(:,:,:) :: ZZDELTA=>NULL()

ZZDELTA => WORK1
IF (PREDICTOR) THEN
   RHOP    => RHOS
   ZZP     => ZZS
ELSE
   RHOP    => RHO
   ZZP     => ZZ
ENDIF

! Search the domain for negative values of Y or Z. Redistribute mass where appropriate.

SPECIESLOOP: DO N=1,N_TRACKED_SPECIES

   !$OMP PARALLEL DEFAULT(NONE) &
   !$OMP SHARED(ZZDELTA,KBAR,JBAR,IBAR,CELL_INDEX,SOLID,WALL_INDEX,ZZP,ZZ_F,N,BOUNDARY_TYPE, &
   !$OMP        DEL_RHO_D_DEL_Z,RHOP)


   !$OMP WORKSHARE
   ZZDELTA = 0._EB
   !$OMP END WORKSHARE

   ! Do undershoots

   !$OMP DO COLLAPSE(3) SCHEDULE(DYNAMIC) &
   !$OMP PRIVATE(K,J,I,IC,IW_A,Y00,SUM,ISUM,LC,ZMIN,ZMI,ZPI,ZMJ,ZPJ,ZMK,ZPK, &
   !$OMP         RHY0,RHZPI,RHZMI,RHZPJ,RHZMJ,RHZPK,RHZMK,CONST)
   DO K=1,KBAR
      DO J=1,JBAR
         CHECK_LOOP: DO I=1,IBAR
            IC = CELL_INDEX(I,J,K)
            IF (SOLID(IC)) CYCLE CHECK_LOOP
            IW_A = WALL_INDEX(IC,:)
            Y00   = ZZP(I,J,K,N)
            SUM   = 0._EB
            ISUM  = 0
            LC    = .FALSE.
            ZMIN  = 1._EB 
            IF (IW_A(-1) == 0) THEN
               ZMI = ZZP(I-1,J,K,N)
               LC(-1) = .TRUE.
            ELSE
               IF (BOUNDARY_TYPE(IW_A(-1))==POROUS_BOUNDARY) THEN
                 ZMI = ZZP(I-1,J,K,N)
                 LC(-1) = .TRUE.
               ELSE
                 ZMI = ZZ_F(IW_A(-1),N)  
               ENDIF
            ENDIF          
            IF (IW_A( 1) == 0) THEN
               ZPI = ZZP(I+1,J,K,N)
               LC( 1) = .TRUE.
            ELSE
               IF (BOUNDARY_TYPE(IW_A(1))==POROUS_BOUNDARY) THEN
                 ZPI = ZZP(I+1,J,K,N)
                 LC( 1) = .TRUE.
               ELSE
                 ZPI = ZZ_F(IW_A(1),N)  
               ENDIF
            ENDIF           
            IF (IW_A(-2) == 0) THEN
               ZMJ = ZZP(I,J-1,K,N)
               LC(-2) = .TRUE.
            ELSE
               IF (BOUNDARY_TYPE(IW_A(-2))==POROUS_BOUNDARY) THEN
                 ZMJ = ZZP(I,J-1,K,N)
                 LC(-2) = .TRUE.
               ELSE
                 ZMJ = ZZ_F(IW_A(-2),N)  
               ENDIF
            ENDIF         
            IF (IW_A( 2) == 0) THEN
               ZPJ = ZZP(I,J+1,K,N)
               LC( 2) = .TRUE.
            ELSE
               IF (BOUNDARY_TYPE(IW_A( 2))==POROUS_BOUNDARY) THEN
                 ZPJ = ZZP(I,J+1,K,N)
                 LC( 2) = .TRUE.
               ELSE
                 ZPJ = ZZ_F(IW_A( 2),N)  
               ENDIF
            ENDIF         
            IF (IW_A(-3) == 0) THEN
               ZMK = ZZP(I,J,K-1,N)
               LC(-3) = .TRUE.
            ELSE
               IF (BOUNDARY_TYPE(IW_A(-3))==POROUS_BOUNDARY) THEN
                 ZMK = ZZP(I,J,K-1,N)
                 LC(-3) = .TRUE.
               ELSE
                 ZMK = ZZ_F(IW_A(-3),N)  
               ENDIF
            ENDIF         
            IF (IW_A( 3) == 0) THEN
               ZPK = ZZP(I,J,K+1,N)
               LC( 3) = .TRUE.
            ELSE
               IF (BOUNDARY_TYPE(IW_A( 3))==POROUS_BOUNDARY) THEN
                 ZPK = ZZP(I,J,K+1,N)
                 LC( 3) = .TRUE.
               ELSE
                 ZPK = ZZ_F(IW_A( 3),N)  
               ENDIF
            ENDIF           
            ZMIN = MIN(ZMI,ZPI,ZMJ,ZPJ,ZMK,ZPK)
            ZMIN = MAX(ZMIN,0._EB)
            IF ((DEL_RHO_D_DEL_Z(I,J,K,N) > 0._EB .AND. Y00 < ZMIN) .OR. Y00 < 0._EB) THEN
               RHY0  = RHOP(I,J,K)  *(ZMIN - Y00)
               IF (LC(-1) .AND. ZMI>ZMIN) THEN
                  RHZMI = RHOP(I-1,J,K)*(ZMI - ZMIN)
                  SUM  = SUM + RHZMI 
                  ISUM = ISUM + 1
               ELSE
                  LC(-1) = .FALSE.
               ENDIF
               IF (LC( 1) .AND. ZPI>ZMIN) THEN
                  RHZPI = RHOP(I+1,J,K)*(ZPI - ZMIN)
                  SUM  = SUM + RHZPI
                  ISUM = ISUM + 1
               ELSE
                  LC( 1) = .FALSE.
               ENDIF
               IF (LC(-2) .AND. ZMJ>ZMIN) THEN
                  RHZMJ = RHOP(I,J-1,K)*(ZMJ - ZMIN)
                  SUM  = SUM + RHZMJ
                  ISUM = ISUM + 1
               ELSE
                  LC(-2) = .FALSE.
               ENDIF
               IF (LC( 2) .AND. ZPJ>ZMIN) THEN
                  RHZPJ = RHOP(I,J+1,K)*(ZPJ - ZMIN)
                  SUM  = SUM + RHZPJ
                  ISUM = ISUM + 1
                  LC( 2) = .TRUE.
               ELSE
                  LC( 2) = .FALSE.
               ENDIF
               IF (LC(-3) .AND. ZMK>ZMIN) THEN
               RHZMK = RHOP(I,J,K-1)*(ZMK - ZMIN)
                  SUM  = SUM + RHZMK
                  ISUM = ISUM + 1
               ELSE
                  LC(-3) = .FALSE.
               ENDIF
               IF (LC( 3) .AND. ZPK>ZMIN) THEN
                  RHZPK = RHOP(I,J,K+1)*(ZPK - ZMIN)
                  SUM  = SUM + RHZPK
                  ISUM = ISUM + 1
               ELSE
                  LC( 3) = .FALSE.
               ENDIF                
               IF (ISUM==0) THEN
                  IF (ZMIN <= 0._EB) ZZDELTA(I,J,K) = ZZDELTA(I,J,K) + ZMIN - Y00  
                  CYCLE CHECK_LOOP
               ELSE
                  IF (ABS(SUM)>=ZERO_P) THEN
                     ZZDELTA(I,J,K) = ZZDELTA(I,J,K) + ZMIN - Y00
                     CONST = MIN(1._EB,RHY0/SUM)
                     IF (LC(-1)) ZZDELTA(I-1,J,K) = ZZDELTA(I-1,J,K) - RHZMI*CONST/RHOP(I-1,J,K)
                     IF (LC( 1)) ZZDELTA(I+1,J,K) = ZZDELTA(I+1,J,K) - RHZPI*CONST/RHOP(I+1,J,K)
                     IF (LC(-2)) ZZDELTA(I,J-1,K) = ZZDELTA(I,J-1,K) - RHZMJ*CONST/RHOP(I,J-1,K)
                     IF (LC( 2)) ZZDELTA(I,J+1,K) = ZZDELTA(I,J+1,K) - RHZPJ*CONST/RHOP(I,J+1,K)
                     IF (LC(-3)) ZZDELTA(I,J,K-1) = ZZDELTA(I,J,K-1) - RHZMK*CONST/RHOP(I,J,K-1)
                     IF (LC( 3)) ZZDELTA(I,J,K+1) = ZZDELTA(I,J,K+1) - RHZPK*CONST/RHOP(I,J,K+1)
                  ENDIF                  
               ENDIF
            ENDIF
         ENDDO CHECK_LOOP
      ENDDO
   ENDDO
   !$OMP END DO
   
   !$OMP WORKSHARE
   ZZP(1:IBAR,1:JBAR,1:KBAR,N) = ZZP(1:IBAR,1:JBAR,1:KBAR,N) + ZZDELTA(1:IBAR,1:JBAR,1:KBAR)
   ZZDELTA=0._EB
   !$OMP END WORKSHARE

   ! Do overshoots
   !$OMP DO COLLAPSE(3) SCHEDULE(DYNAMIC) &
   !$OMP PRIVATE(K,J,I,IC,IW_A,Y00,SUM,ISUM,LC,ZMIN,ZMI,ZPI,ZMK,ZPK,ZMJ,ZPJ,ZMAX, &
   !$OMP         RHY0,RHZPI,RHZMI,RHZPJ,RHZMJ,RHZPK,RHZMK,CONST)
   DO K=1,KBAR
      DO J=1,JBAR
         CHECK_LOOP2: DO I=1,IBAR
            IC = CELL_INDEX(I,J,K)
            IF (SOLID(IC)) CYCLE CHECK_LOOP2
            IW_A  = WALL_INDEX(IC,:)
            Y00   = ZZP(I,J,K,N)
            SUM   = 0._EB
            ISUM  = 0
            LC    = .FALSE.
            ZMIN  = 1._EB 
            IF (IW_A(-1) == 0) THEN
               ZMI = ZZP(I-1,J,K,N)
               LC(-1) = .TRUE.
            ELSE
               ZMI = ZZ_F(IW_A(-1),N)  
            ENDIF          
            IF (IW_A( 1) == 0) THEN
               ZPI = ZZP(I+1,J,K,N)
               LC( 1) = .TRUE.
            ELSE
               ZPI = ZZ_F(IW_A( 1),N)  
            ENDIF           
            IF (IW_A(-2) == 0) THEN
               ZMJ = ZZP(I,J-1,K,N)
               LC(-2) = .TRUE.
            ELSE
               ZMJ = ZZ_F(IW_A(-2),N)  
            ENDIF         
            IF (IW_A( 2) == 0) THEN
               ZPJ = ZZP(I,J+1,K,N)
               LC( 2) = .TRUE.
            ELSE
               ZPJ = ZZ_F(IW_A( 2),N)  
            ENDIF         
            IF (IW_A(-3) == 0) THEN
               ZMK = ZZP(I,J,K-1,N)
               LC(-3) = .TRUE.
            ELSE
               ZMK = ZZ_F(IW_A(-3),N)  
            ENDIF         
            IF (IW_A( 3) == 0) THEN
               ZPK = ZZP(I,J,K+1,N)
               LC( 3) = .TRUE.
            ELSE
               ZPK = ZZ_F(IW_A( 3),N)  
            ENDIF           
            ZMAX = MAX(ZMI,ZPI,ZMJ,ZPJ,ZMK,ZPK)
            ZMAX = MIN(ZMAX,1._EB)            
            IF ((DEL_RHO_D_DEL_Z(I,J,K,N) < 0._EB .AND. Y00 > ZMAX) .OR. Y00 > 1._EB) THEN
               RHY0  = RHOP(I,J,K)  *(Y00 - ZMAX)
               IF (LC(-1) .AND. ZMI<ZMAX) THEN
                  RHZMI = RHOP(I-1,J,K)*(ZMAX - ZMI)
                  SUM  = SUM + RHZMI
                  ISUM = ISUM + 1
               ELSE
                  LC(-1) = .FALSE.
               ENDIF
               IF (LC( 1) .AND. ZPI<ZMAX) THEN
                  RHZPI = RHOP(I+1,J,K)*(ZMAX - ZPI)
                  SUM  = SUM + RHZPI
                  ISUM = ISUM + 1
               ELSE
                  LC( 1) = .FALSE.
               ENDIF
               IF (LC(-2) .AND. ZMJ<ZMAX) THEN
                  RHZMJ = RHOP(I,J-1,K)*(ZMAX - ZMJ)
                  SUM  = SUM + RHZMJ
                  ISUM = ISUM + 1
               ELSE
                  LC(-2) = .FALSE.
               ENDIF
               IF (LC( 2) .AND. ZPJ<ZMAX) THEN
                  RHZPJ = RHOP(I,J+1,K)*(ZMAX - ZPJ)
                  SUM  = SUM + RHZPJ
                  ISUM = ISUM + 1
               ELSE
                  LC( 2) = .FALSE.
               ENDIF
               IF (LC(-3) .AND. ZMK<ZMAX) THEN
                  RHZMK = RHOP(I,J,K-1)*(ZMAX - ZMK)
                  SUM  = SUM + RHZMK
                  ISUM = ISUM + 1
               ELSE
                  LC(-3) = .FALSE.
               ENDIF
               IF (LC( 3) .AND. ZPK<ZMAX) THEN
                  RHZPK = RHOP(I,J,K+1)*(ZMAX - ZPK)
                  SUM  = SUM + RHZPK
                  ISUM = ISUM + 1
               ELSE
                  LC( 3) = .FALSE.
               ENDIF                      
               IF (ISUM==0) THEN
                  IF(ZMAX >= 1._EB) ZZDELTA(I,J,K) = ZZDELTA(I,J,K) + ZMAX - Y00
                  CYCLE CHECK_LOOP2
               ELSE
                  IF (ABS(SUM)>=ZERO_P) THEN
                     ZZDELTA(I,J,K) = ZZDELTA(I,J,K) + ZMAX - Y00               
                     CONST = MIN(1._EB,RHY0/SUM)
                     IF (LC(-1)) ZZDELTA(I-1,J,K) = ZZDELTA(I-1,J,K) + RHZMI*CONST/RHOP(I-1,J,K)
                     IF (LC( 1)) ZZDELTA(I+1,J,K) = ZZDELTA(I+1,J,K) + RHZPI*CONST/RHOP(I+1,J,K)
                     IF (LC(-2)) ZZDELTA(I,J-1,K) = ZZDELTA(I,J-1,K) + RHZMJ*CONST/RHOP(I,J-1,K)
                     IF (LC( 2)) ZZDELTA(I,J+1,K) = ZZDELTA(I,J+1,K) + RHZPJ*CONST/RHOP(I,J+1,K)
                     IF (LC(-3)) ZZDELTA(I,J,K-1) = ZZDELTA(I,J,K-1) + RHZMK*CONST/RHOP(I,J,K-1)
                     IF (LC( 3)) ZZDELTA(I,J,K+1) = ZZDELTA(I,J,K+1) + RHZPK*CONST/RHOP(I,J,K+1)
                  ENDIF
               ENDIF
            ENDIF
         ENDDO CHECK_LOOP2
      ENDDO
   ENDDO  
   !$OMP END DO 

   !$OMP WORKSHARE
   ZZP(1:IBAR,1:JBAR,1:KBAR,N) = ZZP(1:IBAR,1:JBAR,1:KBAR,N) + ZZDELTA(1:IBAR,1:JBAR,1:KBAR)
   !$OMP END WORKSHARE NOWAIT
   !$OMP END PARALLEL
ENDDO SPECIESLOOP

RETURN

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
