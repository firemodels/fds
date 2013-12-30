
! This module is an experimental implementation of my embedded mesh method (EMB),
! a prelude to adaptive mesh refinement.

MODULE EMBEDDED_MESH_METHOD

USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_VARIABLES
USE MESH_POINTERS

IMPLICIT NONE

CHARACTER(255), PARAMETER :: samrid='$Id$'
CHARACTER(255), PARAMETER :: samrrev='$Revision$'
CHARACTER(255), PARAMETER :: samrdate='$Date$'

PRIVATE
PUBLIC GET_REV_samr,SCALARF_EMB,VELOCITY_EMB,RESTRICT_MASS_EMB,RESTRICT_DIV_EMB,PROJECT_VELOCITY, &
       SORT_MESH_LEVEL,MATCH_VELOCITY_EMB,SCALAR_GHOST_EMB

CONTAINS


SUBROUTINE SCALARF_EMB(NM1,NM2,IERROR)

INTEGER, INTENT(IN) :: NM1,NM2

TYPE(MESH_TYPE), POINTER :: M1,M2
INTEGER :: N,I,J,K,I_LO,I_HI,J_LO,J_HI,K_LO,K_HI,II_0,JJ_0,KK_0,II,JJ,KK, &
          NRX,NRY,NRZ,N2X,N2Y,N2Z,II_LO,JJ_LO,KK_LO,INDEX_LIST(12),IERROR
REAL(EB) :: VOLUME_LIST(3)
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: FX1,FY1,FZ1,FX2,FY2,FZ2

!   Comments:
!
!   Assumes uniform grid in each direction and that M2 lies within M1.
!
!   -------------------------------
!   |         |         |         |
!   |         |         |         |<---MESHES(M1)
!   |         |         |         |
!   |         |         |         |
!   -------------------------------
!   |         |-|-|-|-|-|         |
!   |         |-|-|-|-|-|<-------------MESHES(M2)
!   |         |-|-|-|-|-|         |
!   |         |-|-|-|-|-|         |
!   -------------------------------
!   |         |         |         |
!   |         |         |         |
!   |         |         |         |
!   |         |         |         |
!   -------------------------------

CALL LOCATE_MESH(INDEX_LIST,VOLUME_LIST,NM1,NM2,IERROR)
SELECT CASE (IERROR)
  CASE(0)
     I_LO = INDEX_LIST(1)
     I_HI = INDEX_LIST(2)
     J_LO = INDEX_LIST(3)
     J_HI = INDEX_LIST(4)
     K_LO = INDEX_LIST(5)
     K_HI = INDEX_LIST(6)
     II_LO = INDEX_LIST(7)
     JJ_LO = INDEX_LIST(8)
     KK_LO = INDEX_LIST(9)
     NRX = INDEX_LIST(10)
     NRY = INDEX_LIST(11)
     NRZ = INDEX_LIST(12)
  CASE(1)
     RETURN
END SELECT

M1=>MESHES(NM1) ! coarse mesh
M2=>MESHES(NM2) ! fine mesh

N2X = NRY*NRZ
N2Y = NRX*NRZ
N2Z = NRX*NRY

! Fluxes

FX1=>M1%SCALAR_SAVE1
FY1=>M1%SCALAR_SAVE2
FZ1=>M1%SCALAR_SAVE3

FX2=>M2%SCALAR_SAVE1
FY2=>M2%SCALAR_SAVE2
FZ2=>M2%SCALAR_SAVE3

! Restrict fine mesh to coarse mesh for embedded cells

SPECIES_LOOP: DO N=0,N_TRACKED_SPECIES

  ! x-direction fluxes

  DO K = K_LO,K_HI
     KK_0 = KK_LO + (K-K_LO)*NRZ
     DO J = J_LO,J_HI
        JJ_0 = JJ_LO + (J-J_LO)*NRY
        DO I = I_LO-1,I_HI !! note: this includes fine mesh boundary
           II_0 = II_LO + (I-I_LO+1)*NRX !!
                 
           FX1(I,J,K,N) = 0._EB
           DO KK = KK_0+1,KK_0+NRZ
              DO JJ = JJ_0+1,JJ_0+NRY
                 FX1(I,J,K,N) = FX1(I,J,K,N) + FX2(II_0,JJ,KK,N)
              ENDDO
           ENDDO
           FX1(I,J,K,N) = FX1(I,J,K,N)/N2X
        
        ENDDO
     ENDDO
  ENDDO
  
  ! y-direction fluxes

  DO K = K_LO,K_HI
     KK_0 = KK_LO + (K-K_LO)*NRZ
     DO J = J_LO-1,J_HI !!
        JJ_0 = JJ_LO + (J-J_LO+1)*NRY !!
        DO I = I_LO,I_HI 
           II_0 = II_LO + (I-I_LO)*NRX
                 
           FY1(I,J,K,N) = 0._EB
           DO KK = KK_0+1,KK_0+NRZ
              DO II = II_0+1,II_0+NRX
                 FY1(I,J,K,N) = FY1(I,J,K,N) + FY2(II,JJ_0,KK,N)
              ENDDO
           ENDDO
           FY1(I,J,K,N) = FY1(I,J,K,N)/N2Y
        
        ENDDO
     ENDDO
  ENDDO
  
  ! z-direction fluxes

  DO K = K_LO-1,K_HI !!
     KK_0 = KK_LO + (K-K_LO+1)*NRZ !!
     DO J = J_LO,J_HI
        JJ_0 = JJ_LO + (J-J_LO)*NRY
        DO I = I_LO,I_HI 
           II_0 = II_LO + (I-I_LO)*NRX
                 
           FZ1(I,J,K,N) = 0._EB
           DO JJ = JJ_0+1,JJ_0+NRY
              DO II = II_0+1,II_0+NRX
                 FZ1(I,J,K,N) = FZ1(I,J,K,N) + FZ2(II,JJ,KK_0,N)
              ENDDO
           ENDDO
           FZ1(I,J,K,N) = FZ1(I,J,K,N)/N2Z
        
        ENDDO
     ENDDO
  ENDDO
  
ENDDO SPECIES_LOOP

END SUBROUTINE SCALARF_EMB


SUBROUTINE VELOCITY_EMB(NM1,NM2,IERROR)

INTEGER, INTENT(IN) :: NM1,NM2

TYPE(MESH_TYPE), POINTER :: M1,M2
INTEGER :: I,J,K,I_LO,I_HI,J_LO,J_HI,K_LO,K_HI,II_0,JJ_0,KK_0,II,JJ,KK, &
          NRX,NRY,NRZ,N2X,N2Y,N2Z,II_LO,JJ_LO,KK_LO,INDEX_LIST(12),IERROR
REAL(EB) :: VOLUME_LIST(3)
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU1,VV1,WW1,UU2,VV2,WW2

CALL LOCATE_MESH(INDEX_LIST,VOLUME_LIST,NM1,NM2,IERROR)
SELECT CASE (IERROR)
  CASE(0)
     I_LO = INDEX_LIST(1)
     I_HI = INDEX_LIST(2)
     J_LO = INDEX_LIST(3)
     J_HI = INDEX_LIST(4)
     K_LO = INDEX_LIST(5)
     K_HI = INDEX_LIST(6)
     II_LO = INDEX_LIST(7)
     JJ_LO = INDEX_LIST(8)
     KK_LO = INDEX_LIST(9)
     NRX = INDEX_LIST(10)
     NRY = INDEX_LIST(11)
     NRZ = INDEX_LIST(12)
  CASE(1)
     RETURN
END SELECT

M1=>MESHES(NM1) ! coarse mesh
M2=>MESHES(NM2) ! fine mesh

N2X = NRY*NRZ
N2Y = NRX*NRZ
N2Z = NRX*NRY

IF (PREDICTOR) THEN
  UU1=>M1%U
  VV1=>M1%V
  WW1=>M1%W
  UU2=>M2%U
  VV2=>M2%V
  WW2=>M2%W
ELSEIF (CORRECTOR) THEN
  UU1=>M1%US
  VV1=>M1%VS
  WW1=>M1%WS
  UU2=>M2%US
  VV2=>M2%VS
  WW2=>M2%WS
ENDIF

! Restrict fine mesh to coarse mesh for embedded cells

! U-VELOCITY

DO K = K_LO,K_HI
  KK_0 = KK_LO + (K-K_LO)*NRZ
  DO J = J_LO,J_HI
     JJ_0 = JJ_LO + (J-J_LO)*NRY
     DO I = I_LO,I_HI-1 ! excludes boundary values
        II_0 = II_LO + (I-I_LO+1)*NRX
                 
        UU1(I,J,K) = 0._EB
        DO KK = KK_0+1,KK_0+NRZ
           DO JJ = JJ_0+1,JJ_0+NRY
              UU1(I,J,K) = UU1(I,J,K) + UU2(II_0,JJ,KK)
           ENDDO
        ENDDO
        UU1(I,J,K) = UU1(I,J,K)/N2X
        
     ENDDO
  ENDDO
ENDDO
  
! V-VELOCITY

DO K = K_LO,K_HI
  KK_0 = KK_LO + (K-K_LO)*NRZ
  DO J = J_LO,J_HI-1 ! excludes boundary values
     JJ_0 = JJ_LO + (J-J_LO+1)*NRY
     DO I = I_LO,I_HI 
        II_0 = II_LO + (I-I_LO)*NRX
                 
        VV1(I,J,K) = 0._EB
        DO KK = KK_0+1,KK_0+NRZ
           DO II = II_0+1,II_0+NRX
              VV1(I,J,K) = VV1(I,J,K) + VV2(II,JJ_0,KK)
           ENDDO
        ENDDO
        VV1(I,J,K) = VV1(I,J,K)/N2Y
        
     ENDDO
  ENDDO
ENDDO
  
! W-VELOCITY

DO K = K_LO,K_HI-1 ! excludes boundary values
  KK_0 = KK_LO + (K-K_LO+1)*NRZ
  DO J = J_LO,J_HI
     JJ_0 = JJ_LO + (J-J_LO)*NRY
     DO I = I_LO,I_HI 
        II_0 = II_LO + (I-I_LO)*NRX
                 
        WW1(I,J,K) = 0._EB
        DO JJ = JJ_0+1,JJ_0+NRY
           DO II = II_0+1,II_0+NRX
              WW1(I,J,K) = WW1(I,J,K) + WW2(II,JJ,KK_0)
           ENDDO
        ENDDO
        WW1(I,J,K) = WW1(I,J,K)/N2Z
        
     ENDDO
  ENDDO
ENDDO

END SUBROUTINE VELOCITY_EMB


SUBROUTINE RESTRICT_MASS_EMB(NM1,NM2,IERROR)

INTEGER, INTENT(IN) :: NM1,NM2

TYPE(MESH_TYPE), POINTER :: M1,M2
INTEGER :: N,I,J,K,I_LO,I_HI,J_LO,J_HI,K_LO,K_HI,II_0,JJ_0,KK_0,II,JJ,KK, &
          NRX,NRY,NRZ, II_LO,JJ_LO,KK_LO,INDEX_LIST(12),IERROR
REAL(EB) :: DV1,DV2,DVRAT,VOLUME_LIST(3)
REAL(EB), POINTER, DIMENSION(:,:,:) :: RHO1,RHO2
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZ1,ZZ2

CALL LOCATE_MESH(INDEX_LIST,VOLUME_LIST,NM1,NM2,IERROR)
SELECT CASE (IERROR)
  CASE(0)
     I_LO = INDEX_LIST(1)
     I_HI = INDEX_LIST(2)
     J_LO = INDEX_LIST(3)
     J_HI = INDEX_LIST(4)
     K_LO = INDEX_LIST(5)
     K_HI = INDEX_LIST(6)
     II_LO = INDEX_LIST(7)
     JJ_LO = INDEX_LIST(8)
     KK_LO = INDEX_LIST(9)
     NRX = INDEX_LIST(10)
     NRY = INDEX_LIST(11)
     NRZ = INDEX_LIST(12)
     DV1 = VOLUME_LIST(1)
     DV2 = VOLUME_LIST(2)
     DVRAT = VOLUME_LIST(3)
  CASE(1)
     RETURN
END SELECT

M1=>MESHES(NM1) ! coarse mesh
M2=>MESHES(NM2) ! fine mesh

IF (PREDICTOR) THEN
  RHO1 => M1%RHOS
  RHO2 => M2%RHOS
  IF (N_TRACKED_SPECIES>0) ZZ1  => M1%ZZS
  IF (N_TRACKED_SPECIES>0) ZZ2  => M2%ZZS
ELSEIF (CORRECTOR) THEN
  RHO1 => M1%RHO
  RHO2 => M2%RHO
  IF (N_TRACKED_SPECIES>0) ZZ1  => M1%ZZ
  IF (N_TRACKED_SPECIES>0) ZZ2  => M2%ZZ
ENDIF

DO K = K_LO,K_HI
  KK_0 = KK_LO + (K-K_LO)*NRZ
  DO J = J_LO,J_HI
     JJ_0 = JJ_LO + (J-J_LO)*NRY
     DO I = I_LO,I_HI
        II_0 = II_LO + (I-I_LO)*NRX
           
        RHO1(I,J,K) = 0._EB
        
        DO KK = KK_0+1,KK_0+NRZ
           DO JJ = JJ_0+1,JJ_0+NRY
              DO II = II_0+1,II_0+NRX
                
                 RHO1(I,J,K) = RHO1(I,J,K) + RHO2(II,JJ,KK)*DVRAT
                    
              ENDDO
           ENDDO
        ENDDO
     
     ENDDO
  ENDDO
ENDDO

IF (N_TRACKED_SPECIES>0) THEN

  SPECIES_LOOP: DO N=1,N_TRACKED_SPECIES
  
     DO K = K_LO,K_HI
        KK_0 = KK_LO + (K-K_LO)*NRZ
        DO J = J_LO,J_HI
           JJ_0 = JJ_LO + (J-J_LO)*NRY
           DO I = I_LO,I_HI
              II_0 = II_LO + (I-I_LO)*NRX
           
              ZZ1(I,J,K,N) = 0._EB
        
              DO KK = KK_0+1,KK_0+NRZ
                 DO JJ = JJ_0+1,JJ_0+NRY
                    DO II = II_0+1,II_0+NRX
                
                       ZZ1(I,J,K,N) = ZZ1(I,J,K,N) + RHO2(II,JJ,KK)*ZZ2(II,JJ,KK,N)*DV2
                    
                    ENDDO
                 ENDDO
              ENDDO
              
              ZZ1(I,J,K,N) = ZZ1(I,J,K,N)/(RHO1(I,J,K)*DV1)
     
           ENDDO
        ENDDO
     ENDDO

  ENDDO SPECIES_LOOP

ENDIF

END SUBROUTINE RESTRICT_MASS_EMB


SUBROUTINE RESTRICT_DIV_EMB(NM1,NM2,IERROR)

INTEGER, INTENT(IN) :: NM1,NM2

TYPE(MESH_TYPE), POINTER :: M1,M2
INTEGER :: I,J,K,I_LO,I_HI,J_LO,J_HI,K_LO,K_HI,II_0,JJ_0,KK_0,II,JJ,KK, &
          NRX,NRY,NRZ,II_LO,JJ_LO,KK_LO,INDEX_LIST(12),IERROR
REAL(EB) :: DVRAT,VOLUME_LIST(3)
REAL(EB), POINTER, DIMENSION(:,:,:) :: DP1,DP2

CALL LOCATE_MESH(INDEX_LIST,VOLUME_LIST,NM1,NM2,IERROR)
SELECT CASE (IERROR)
  CASE(0)
     I_LO = INDEX_LIST(1)
     I_HI = INDEX_LIST(2)
     J_LO = INDEX_LIST(3)
     J_HI = INDEX_LIST(4)
     K_LO = INDEX_LIST(5)
     K_HI = INDEX_LIST(6)
     II_LO = INDEX_LIST(7)
     JJ_LO = INDEX_LIST(8)
     KK_LO = INDEX_LIST(9)
     NRX = INDEX_LIST(10)
     NRY = INDEX_LIST(11)
     NRZ = INDEX_LIST(12)
     DVRAT = VOLUME_LIST(3)
  CASE(1)
     RETURN
END SELECT

M1=>MESHES(NM1) ! coarse mesh
M2=>MESHES(NM2) ! fine mesh

IF (PREDICTOR) THEN
  DP1 => M1%DS
  DP2 => M2%DS
ELSEIF (CORRECTOR) THEN
  DP1 => M1%DDDT
  DP2 => M2%DDDT
ENDIF

! Restrict divergence
  
DO K = K_LO,K_HI
  KK_0 = KK_LO + (K-K_LO)*NRZ
  DO J = J_LO,J_HI
     JJ_0 = JJ_LO + (J-J_LO)*NRY
     DO I = I_LO,I_HI
        II_0 = II_LO + (I-I_LO)*NRX
           
        DP1(I,J,K) = 0._EB
        
        DO KK = KK_0+1,KK_0+NRZ
           DO JJ = JJ_0+1,JJ_0+NRY
              DO II = II_0+1,II_0+NRX
                
                 DP1(I,J,K) = DP1(I,J,K) + DP2(II,JJ,KK)*DVRAT
                    
              ENDDO
           ENDDO
        ENDDO
     
     ENDDO
  ENDDO
ENDDO

END SUBROUTINE RESTRICT_DIV_EMB


SUBROUTINE PROJECT_VELOCITY(NM)
USE POIS, ONLY: H3CZSS,H2CZSS

INTEGER, INTENT(IN) :: NM
INTEGER :: I,J,K
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW,DP,PP,PRHS_SAVE
REAL(EB) :: DIV,LHSS,RES,POIS_ERR

CALL POINT_TO_MESH(NM)

IF (PREDICTOR) THEN
  ! note: PROJECT_VELOCITY is called AFTER the predictor update of velocity
  UU=>US
  VV=>VS
  WW=>WS
  DP=>D
ELSEIF (CORRECTOR) THEN
  UU=>U
  VV=>V
  WW=>W
  DP=>DS
ENDIF
PP=>WORK1
PRHS_SAVE=>WORK2

! build source

DO K=1,KBAR
  DO J=1,JBAR
     DO I=1,IBAR
        DIV = (UU(I,J,K)-UU(I-1,J,K))*RDX(I) + (VV(I,J,K)-VV(I,J-1,K))*RDY(J) + (WW(I,J,K)-WW(I,J,K-1))*RDZ(K)
        PRHS(I,J,K) = DIV-DP(I,J,K)
     ENDDO
  ENDDO
ENDDO

! solve Poisson equation

BXS = 0._EB
BXF = 0._EB
BYS = 0._EB
BYF = 0._EB
BZS = 0._EB
BZF = 0._EB

PRHS_SAVE(1:IBAR,1:JBAR,1:KBAR) = PRHS
IF (.NOT.TWO_D) CALL H3CZSS(BXS,BXF,BYS,BYF,BZS,BZF,ITRN,JTRN,PRHS,POIS_PTB,SAVE2,WORK,HX)
IF (TWO_D .AND. .NOT. CYLINDRICAL) CALL H2CZSS(BXS,BXF,BZS,BZF,ITRN,PRHS,POIS_PTB,SAVE2,WORK,HX)
PP(1:IBAR,1:JBAR,1:KBAR) = PRHS

! Apply boundary conditions to PP

DO K=1,KBAR
  DO J=1,JBAR
     PP(0,J,K)    = -PP(1,J,K) ! use minus if Dirichlet, plus if Neumann, see init of SAVE2
     PP(IBP1,J,K) = -PP(IBAR,J,K)
  ENDDO
ENDDO

DO K=1,KBAR
  DO I=1,IBAR
     PP(I,0,K)    = -PP(I,1,K)
     PP(I,JBP1,K) = -PP(I,JBAR,K)
  ENDDO
ENDDO

DO J=1,JBAR
  DO I=1,IBAR
     PP(I,J,0)    = -PP(I,J,1)
     PP(I,J,KBP1) = -PP(I,J,KBAR)
  ENDDO
ENDDO

! ************************* Check the Solution *************************

IF (.FALSE.) THEN

  POIS_ERR = 0._EB
  DO K=1,KBAR
     DO J=1,JBAR
        DO I=1,IBAR
           LHSS = ((PP(I+1,J,K)-PP(I,J,K))*RDXN(I) - (PP(I,J,K)-PP(I-1,J,K))*RDXN(I-1) )*RDX(I) &
                + ((PP(I,J+1,K)-PP(I,J,K))*RDYN(J) - (PP(I,J,K)-PP(I,J-1,K))*RDYN(J-1) )*RDY(J) &
                + ((PP(I,J,K+1)-PP(I,J,K))*RDZN(K) - (PP(I,J,K)-PP(I,J,K-1))*RDZN(K-1) )*RDZ(K)
           RES = ABS(PRHS_SAVE(I,J,K)-LHSS)
           POIS_ERR = MAX(RES,POIS_ERR)
        ENDDO
     ENDDO
  ENDDO
  WRITE(0,*) 'POIS ERROR:',pois_ptb,pois_err

ENDIF

! correct velocities

DO K=1,KBAR
  DO J=1,JBAR
     DO I=0,IBAR
        UU(I,J,K) = UU(I,J,K) - RDXN(I)*(PP(I+1,J,K)-PP(I,J,K))
     ENDDO
  ENDDO
ENDDO

DO K=1,KBAR
  DO J=0,JBAR
     DO I=1,IBAR
        VV(I,J,K) = VV(I,J,K) - RDYN(J)*(PP(I,J+1,K)-PP(I,J,K))
     ENDDO
  ENDDO
ENDDO

DO K=0,KBAR
  DO J=1,JBAR
     DO I=1,IBAR
        WW(I,J,K) = WW(I,J,K) - RDZN(K)*(PP(I,J,K+1)-PP(I,J,K))
     ENDDO
  ENDDO
ENDDO

! check divergence

IF (.FALSE.) THEN
  POIS_ERR = 0._EB
  DO K=1,KBAR
     DO J=1,JBAR
        DO I=1,IBAR
           DIV = (UU(I,J,K)-UU(I-1,J,K))*RDX(I) + (VV(I,J,K)-VV(I,J-1,K))*RDY(J) + (WW(I,J,K)-WW(I,J,K-1))*RDZ(K)
           RES = ABS(DIV-DP(I,J,K))
           POIS_ERR = MAX(RES,POIS_ERR)
        ENDDO
     ENDDO
  ENDDO
  WRITE(0,*) NM,MAXVAL(ABS(DP)),POIS_ERR

ENDIF

END SUBROUTINE PROJECT_VELOCITY


SUBROUTINE SORT_MESH_LEVEL

INTEGER :: IRANK,NM,ML,MLMIN,MLMAX

MESH_LIST_EMB = 0

MLMAX = MAXVAL(MESHES(1:NMESHES)%MESH_LEVEL)
MLMIN = MINVAL(MESHES(1:NMESHES)%MESH_LEVEL)

IRANK=0

DO ML=MLMAX,MLMIN,-1
  DO NM=1,NMESHES
  
     IF (MESHES(NM)%MESH_LEVEL==ML) THEN
        IRANK=IRANK+1
        MESH_LIST_EMB(IRANK)=NM
     ENDIF

  ENDDO
ENDDO

!PRINT *,MLMIN,MLMAX
!DO IRANK=1,NMESHES
!   PRINT *,MESH_LIST_EMB(IRANK)
!ENDDO
!STOP

END SUBROUTINE SORT_MESH_LEVEL


SUBROUTINE MATCH_VELOCITY_EMB(NM1,NM2,IERROR)

INTEGER, INTENT(IN) :: NM1,NM2

TYPE(MESH_TYPE), POINTER :: M1,M2
INTEGER :: I,J,K,I_LO,I_HI,J_LO,J_HI,K_LO,K_HI,II_0,JJ_0,KK_0,II,JJ,KK, &
          NRX,NRY,NRZ,II_LO,JJ_LO,KK_LO,INDEX_LIST(12),IERROR,IW,IOR
REAL(EB) :: VOLUME_LIST(3)
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU1,VV1,WW1,UU2,VV2,WW2
REAL(EB), PARAMETER :: RF=0.5_EB,OMRF=0.5_EB
TYPE(WALL_TYPE), POINTER :: WC=>NULL()

CALL LOCATE_MESH(INDEX_LIST,VOLUME_LIST,NM1,NM2,IERROR)
SELECT CASE (IERROR)
  CASE(0)
     I_LO = INDEX_LIST(1)
     I_HI = INDEX_LIST(2)
     J_LO = INDEX_LIST(3)
     J_HI = INDEX_LIST(4)
     K_LO = INDEX_LIST(5)
     K_HI = INDEX_LIST(6)
     II_LO = INDEX_LIST(7)
     JJ_LO = INDEX_LIST(8)
     KK_LO = INDEX_LIST(9)
     NRX = INDEX_LIST(10)
     NRY = INDEX_LIST(11)
     NRZ = INDEX_LIST(12)
  CASE(1)
     RETURN
END SELECT

M1=>MESHES(NM1) ! coarse mesh
M2=>MESHES(NM2) ! fine mesh

IF (PREDICTOR) THEN
  UU1=>M1%US
  VV1=>M1%VS
  WW1=>M1%WS
  UU2=>M2%US
  VV2=>M2%VS
  WW2=>M2%WS
ELSEIF (CORRECTOR) THEN
  UU1=>M1%U
  VV1=>M1%V
  WW1=>M1%W
  UU2=>M2%U
  VV2=>M2%V
  WW2=>M2%W
ENDIF

! Set fine mesh boundary value to corresponding coarse mesh value

! U-VELOCITY

DO K = K_LO,K_HI
  KK_0 = KK_LO + (K-K_LO)*NRZ
  DO J = J_LO,J_HI
     JJ_0 = JJ_LO + (J-J_LO)*NRY

     ! east face
     I = I_HI
     II_0 = II_LO + (I-I_LO+1)*NRX
     IF (II_0==M2%IBAR) THEN
        DO KK = KK_0+1,KK_0+NRZ
           DO JJ = JJ_0+1,JJ_0+NRY
              UU2(II_0,JJ,KK) = UU1(I,J,K)
           ENDDO
        ENDDO
     ENDIF
        
     ! west face
     I = I_LO-1
     II_0 = II_LO + (I-I_LO+1)*NRX
     IF (II_0==0) THEN
        DO KK = KK_0+1,KK_0+NRZ
           DO JJ = JJ_0+1,JJ_0+NRY
              UU2(II_0,JJ,KK) = UU1(I,J,K)
           ENDDO
        ENDDO
     ENDIF
        
  ENDDO
ENDDO
  
! V-VELOCITY

DO K = K_LO,K_HI
  KK_0 = KK_LO + (K-K_LO)*NRZ
  DO I = I_LO,I_HI
     II_0 = II_LO + (I-I_LO)*NRX

     ! north face
     J = J_HI
     JJ_0 = JJ_LO + (J-J_LO+1)*NRY
     IF (JJ_0==M2%JBAR) THEN
        DO KK = KK_0+1,KK_0+NRZ
           DO II = II_0+1,II_0+NRX
              VV2(II,JJ_0,KK) = VV1(I,J,K)
           ENDDO
        ENDDO
     ENDIF
        
     ! south face
     J = J_LO-1
     JJ_0 = JJ_LO + (J-J_LO+1)*NRY
     IF (JJ_0==0) THEN
        DO KK = KK_0+1,KK_0+NRZ
           DO II = II_0+1,II_0+NRX
              VV2(II,JJ_0,KK) = VV1(I,J,K)
           ENDDO
        ENDDO
     ENDIF
        
  ENDDO
ENDDO
  
! W-VELOCITY

DO J = J_LO,J_HI
  JJ_0 = JJ_LO + (J-J_LO)*NRY
  DO I = I_LO,I_HI
     II_0 = II_LO + (I-I_LO)*NRX

     ! top face
     K = K_HI
     KK_0 = KK_LO + (K-K_LO+1)*NRZ
     IF (KK_0==M2%KBAR) THEN
        DO JJ = JJ_0+1,JJ_0+NRY
           DO II = II_0+1,II_0+NRX
              WW2(II,JJ,KK_0) = WW1(I,J,K)
           ENDDO
        ENDDO
     ENDIF
        
     ! bottom face
     K = K_LO-1
     KK_0 = KK_LO + (K-K_LO+1)*NRZ
     IF (KK_0==0) THEN
        DO JJ = JJ_0+1,JJ_0+NRY
           DO II = II_0+1,II_0+NRX
              WW2(II,JJ,KK_0) = WW1(I,J,K)
           ENDDO
        ENDDO
     ENDIF
        
  ENDDO
ENDDO

! fine mesh boundary loop

FINE_MESH_WALL_LOOP: DO IW=1,M2%N_EXTERNAL_WALL_CELLS
   WC => WALL(IW)
   II  = WC%ONE_D%II
   JJ  = WC%ONE_D%JJ
   KK  = WC%ONE_D%KK
   IOR = WC%ONE_D%IOR
   SELECT CASE (IOR)
      CASE(1)
         M2%UVW_SAVE(IW)=UU2(0,JJ,KK)
      CASE(-1)
         M2%UVW_SAVE(IW)=UU2(M2%IBAR,JJ,KK)
      CASE(2)
         M2%UVW_SAVE(IW)=UU2(II,0,KK)
      CASE(-2)
         M2%UVW_SAVE(IW)=UU2(II,M2%JBAR,KK)
      CASE(3)
         M2%UVW_SAVE(IW)=UU2(II,JJ,0)
      CASE(-3)
         M2%UVW_SAVE(IW)=UU2(II,JJ,M2%KBAR)
   END SELECT
ENDDO FINE_MESH_WALL_LOOP

END SUBROUTINE MATCH_VELOCITY_EMB


SUBROUTINE SCALAR_GHOST_EMB(NM1,NM2,IERROR)

INTEGER, INTENT(IN) :: NM1,NM2

TYPE(MESH_TYPE), POINTER :: M1=>NULL(),M2=>NULL()
INTEGER :: N,I,J,K,I_LO,I_HI,J_LO,J_HI,K_LO,K_HI,II_0,JJ_0,KK_0,II,JJ,KK, &
          NRX,NRY,NRZ,II_LO,JJ_LO,KK_LO,INDEX_LIST(12),IERROR,IW
REAL(EB) :: VOLUME_LIST(3)
REAL(EB), POINTER, DIMENSION(:,:,:) :: RHOP1,RHOP2,TMP1,TMP2,HH1,HH2
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZP1,ZZP2
TYPE(WALL_TYPE), POINTER :: WC=>NULL()

CALL LOCATE_MESH(INDEX_LIST,VOLUME_LIST,NM1,NM2,IERROR)
SELECT CASE (IERROR)
   CASE(0)
      I_LO = INDEX_LIST(1)
      I_HI = INDEX_LIST(2)
      J_LO = INDEX_LIST(3)
      J_HI = INDEX_LIST(4)
      K_LO = INDEX_LIST(5)
      K_HI = INDEX_LIST(6)
      II_LO = INDEX_LIST(7)
      JJ_LO = INDEX_LIST(8)
      KK_LO = INDEX_LIST(9)
      NRX = INDEX_LIST(10)
      NRY = INDEX_LIST(11)
      NRZ = INDEX_LIST(12)
   CASE(1)
      RETURN
END SELECT

M1=>MESHES(NM1) ! coarse mesh
M2=>MESHES(NM2) ! fine mesh

IF (PREDICTOR) THEN
   RHOP1 => M1%RHOS
   ZZP1  => M1%ZZS
   HH1   => M1%H

   RHOP2 => M2%RHOS
   ZZP2  => M2%ZZS
   HH2   => M2%H
ELSEIF (CORRECTOR) THEN
   RHOP1 => M1%RHO
   ZZP1  => M1%ZZ
   HH1   => M1%HS

   RHOP2 => M2%RHO
   ZZP2  => M2%ZZ
   HH2   => M2%HS
ENDIF
TMP1 => M1%TMP
TMP2 => M2%TMP


! Set fine mesh boundary value to corresponding coarse mesh value

SPECIES_LOOP: DO N=1,N_TRACKED_SPECIES

   DO K = K_LO,K_HI
      KK_0 = KK_LO + (K-K_LO)*NRZ
      DO J = J_LO,J_HI
         JJ_0 = JJ_LO + (J-J_LO)*NRY

         ! east face
         I = I_HI+1
         II_0 = II_LO + (I-I_LO)*NRX + 1
         IF (II_0==M2%IBP1 .AND. I_HI/=M1%IBAR) THEN
         ! if I_HI==M1%IBAR then this might be an external boundary and the ghost cell value
         ! is handled in WALL_BC; similar conditions apply below
            DO KK = KK_0+1,KK_0+NRZ
               DO JJ = JJ_0+1,JJ_0+NRY
                  RHOP2(II_0,JJ,KK) = RHOP1(I,J,K)
                  TMP2(II_0,JJ,KK) = TMP1(I,J,K)
                  HH2(II_0,JJ,KK) = HH1(I,J,K)
                  IF (N_TRACKED_SPECIES>0) ZZP2(II_0,JJ,KK,N) = ZZP1(I,J,K,N)
               ENDDO
            ENDDO
         ENDIF

         ! west face
         I = I_LO-1
         II_0 = II_LO + (I-I_LO+1)*NRX
         IF (II_0==0 .AND. I_LO/=1) THEN
            DO KK = KK_0+1,KK_0+NRZ
               DO JJ = JJ_0+1,JJ_0+NRY
                  RHOP2(II_0,JJ,KK) = RHOP1(I,J,K)
                  TMP2(II_0,JJ,KK) = TMP1(I,J,K)
                  HH2(II_0,JJ,KK) = HH1(I,J,K)
                  IF (N_TRACKED_SPECIES>0) ZZP2(II_0,JJ,KK,N) = ZZP1(I,J,K,N)
               ENDDO
            ENDDO
         ENDIF
        
      ENDDO
   ENDDO

   DO K = K_LO,K_HI
      KK_0 = KK_LO + (K-K_LO)*NRZ
      DO I = I_LO,I_HI
         II_0 = II_LO + (I-I_LO)*NRX

         ! north face
         J = J_HI+1
         JJ_0 = JJ_LO + (J-J_LO)*NRY + 1
         IF (JJ_0==M2%JBP1 .AND. J_HI/=M1%JBAR) THEN
            DO KK = KK_0+1,KK_0+NRZ
               DO II = II_0+1,II_0+NRX
                  RHOP2(II,JJ_0,KK) = RHOP1(I,J,K)
                  TMP2(II,JJ_0,KK) = TMP1(I,J,K)
                  HH2(II,JJ_0,KK) = HH1(I,J,K)
                  IF (N_TRACKED_SPECIES>0) ZZP2(II,JJ_0,KK,N) = ZZP1(I,J,K,N)
               ENDDO
            ENDDO
         ENDIF

         ! south face
         J = J_LO-1
         JJ_0 = JJ_LO + (J-J_LO+1)*NRY
         IF (JJ_0==0 .AND. J_LO/=1) THEN
            DO KK = KK_0+1,KK_0+NRZ
               DO II = II_0+1,II_0+NRX
                  RHOP2(II,JJ_0,KK) = RHOP1(I,J,K)
                  TMP2(II,JJ_0,KK) = TMP1(I,J,K)
                  HH2(II,JJ_0,KK) = HH1(I,J,K)
                  IF (N_TRACKED_SPECIES>0) ZZP2(II,JJ_0,KK,N) = ZZP1(I,J,K,N)
               ENDDO
            ENDDO
         ENDIF
        
      ENDDO
   ENDDO

   DO J = J_LO,J_HI
      JJ_0 = JJ_LO + (J-J_LO)*NRY
      DO I = I_LO,I_HI
         II_0 = II_LO + (I-I_LO)*NRX

         ! top face
         K = K_HI+1
         KK_0 = KK_LO + (K-K_LO)*NRZ + 1
         IF (KK_0==M2%KBP1  .AND. K_HI/=M1%KBAR) THEN
            DO JJ = JJ_0+1,JJ_0+NRY
               DO II = II_0+1,II_0+NRX
                  RHOP2(II,JJ,KK_0) = RHOP1(I,J,K)
                  TMP2(II,JJ,KK_0) = TMP1(I,J,K)
                  HH2(II,JJ,KK_0) = HH1(I,J,K)
                  IF (N_TRACKED_SPECIES>0) ZZP2(II,JJ,KK_0,N) = ZZP1(I,J,K,N)
               ENDDO
            ENDDO
         ENDIF

         ! bottom face
         K = K_LO-1
         KK_0 = KK_LO + (K-K_LO+1)*NRZ
         IF (KK_0==0 .AND. K_LO/=1) THEN
            DO JJ = JJ_0+1,JJ_0+NRY
               DO II = II_0+1,II_0+NRX
                  RHOP2(II,JJ,KK_0) = RHOP1(I,J,K)
                  TMP2(II,JJ,KK_0) = TMP1(I,J,K)
                  HH2(II,JJ,KK_0) = HH1(I,J,K)
                  IF (N_TRACKED_SPECIES>0) ZZP2(II,JJ,KK_0,N) = ZZP1(I,J,K,N)
               ENDDO
            ENDDO
         ENDIF

      ENDDO
   ENDDO

   WALL_LOOP: DO IW=1,M2%N_EXTERNAL_WALL_CELLS+M2%N_INTERNAL_WALL_CELLS
      WC => WALL(IW)
      IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY) CYCLE WALL_LOOP
      II = WC%ONE_D%II
      JJ = WC%ONE_D%JJ
      KK = WC%ONE_D%KK
      WC%RHO_F = RHOP2(II,JJ,KK) 
      WC%ZZ_F(N) = ZZP2(II,JJ,KK,N)
   ENDDO WALL_LOOP

ENDDO SPECIES_LOOP

END SUBROUTINE SCALAR_GHOST_EMB


SUBROUTINE LOCATE_MESH(INDEX_LIST,VOLUME_LIST,NM1,NM2,IERROR)

INTEGER, INTENT(IN) :: NM1,NM2

INTEGER, INTENT(OUT) :: IERROR,INDEX_LIST(12)
REAL(EB), INTENT(OUT) :: VOLUME_LIST(3)
TYPE (MESH_TYPE), POINTER :: M1,M2

INTEGER :: I_LO,I_HI,J_LO,J_HI,K_LO,K_HI,II_LO,JJ_LO,KK_LO,NRX,NRY,NRZ
REAL(EB) :: DV1,DV2,DVRAT

IERROR=0
INDEX_LIST=0
VOLUME_LIST=0._EB

M1=>MESHES(NM1) ! coarse mesh
M2=>MESHES(NM2) ! fine mesh

! Locate fine mesh within coarse mesh

I_LO = MAX(1,       NINT((M2%XS-M1%XS)/M1%DX(1))+1 )
I_HI = MIN(M1%IBAR, NINT((M2%XF-M1%XS)/M1%DX(1))   )
IF (I_LO>M1%IBAR .OR. I_HI<1) THEN ! meshes do not overlap
  IERROR=1
  RETURN
ENDIF

J_LO = MAX(1,       NINT((M2%YS-M1%YS)/M1%DY(1))+1 )
J_HI = MIN(M1%JBAR, NINT((M2%YF-M1%YS)/M1%DY(1))   )
IF (J_LO>M1%JBAR .OR. J_HI<1) THEN ! meshes do not overlap
  IERROR=1
  RETURN
ENDIF

K_LO = MAX(1,       NINT((M2%ZS-M1%ZS)/M1%DZ(1))+1 )
K_HI = MIN(M1%KBAR, NINT((M2%ZF-M1%ZS)/M1%DZ(1))   )
IF (K_LO>M1%KBAR .OR. K_HI<1) THEN ! meshes do not overlap
  IERROR=1
  RETURN
ENDIF

! Find fine mesh off-set

II_LO = MAX(0, NINT((M1%XS-M2%XS)/M2%DX(1)) )
JJ_LO = MAX(0, NINT((M1%YS-M2%YS)/M2%DY(1)) )
KK_LO = MAX(0, NINT((M1%ZS-M2%ZS)/M2%DZ(1)) )

! Compute refinment ratio in each direction

NRX = NINT(M1%DX(1)/M2%DX(1))
NRY = NINT(M1%DY(1)/M2%DY(1))
NRZ = NINT(M1%DZ(1)/M2%DZ(1))

! Cell volumes

DV1 = M1%DX(1)*M1%DY(1)*M1%DZ(1)
DV2 = M2%DX(1)*M2%DY(1)*M2%DZ(1)
DVRAT = DV2/DV1

INDEX_LIST = (/I_LO,I_HI,J_LO,J_HI,K_LO,K_HI,II_LO,JJ_LO,KK_LO,NRX,NRY,NRZ/)
VOLUME_LIST = (/DV1,DV2,DVRAT/)

END SUBROUTINE LOCATE_MESH


SUBROUTINE GET_REV_samr(MODULE_REV,MODULE_DATE)
INTEGER,INTENT(INOUT) :: MODULE_REV
CHARACTER(255),INTENT(INOUT) :: MODULE_DATE

WRITE(MODULE_DATE,'(A)') samrrev(INDEX(samrrev,':')+2:LEN_TRIM(samrrev)-2)
READ (MODULE_DATE,'(I5)') MODULE_REV
WRITE(MODULE_DATE,'(A)') samrdate

END SUBROUTINE GET_REV_samr


END MODULE EMBEDDED_MESH_METHOD
