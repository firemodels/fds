MODULE SCARC_SOLVER
 
USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_POINTERS
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
USE COMP_FUNCTIONS, ONLY: SECOND, GET_FILE_NUMBER
USE MPI
 
IMPLICIT NONE

CHARACTER(255), PARAMETER :: scrchid='$Id$'
CHARACTER(255), PARAMETER :: scrcrev='$Revision$'
CHARACTER(255), PARAMETER :: scrcdate='$Date$'

PRIVATE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Public structures
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! Public subprograms (called from main and pres)
PUBLIC GET_REV_scrc
PUBLIC SCARC_TIMINGS, SCARC_INITIALIZE, SCARC_SHOW, SCARC_SHOW0,SCARC_SHOW02
PUBLIC SCARC_UPDATE , SCARC_UPDATE_LEVEL
PUBLIC SCARC_CG2D   , SCARC_CG3D
PUBLIC SCARC_MG2D   , SCARC_MG3D
PUBLIC SCARC_BICG2D , SCARC_BICG3D
 
!!! Public variables (needed in main, read, pres, divg)
PUBLIC SCARC_METHOD   , SCARC_DEBUG     , SCARC_CASE
PUBLIC SCARC_LU       , SCARC_FN
PUBLIC SCARC_CAPPA    , SCARC_RES       , SCARC_NIT
PUBLIC SCARC_EPS_REL  , SCARC_EPS_DIVG  , SCARC_BREL    
PUBLIC SCARC_MG_NIT   , SCARC_CG_NIT    , SCARC_BICG_NIT    , SCARC_SM_NIT    , SCARC_CO_NIT    
PUBLIC SCARC_MG_EPS   , SCARC_CG_EPS    , SCARC_BICG_EPS    , SCARC_SM_EPS    , SCARC_CO_EPS
PUBLIC SCARC_MG_PRECON, SCARC_CG_PRECON , SCARC_BICG_PRECON , SCARC_SM_PRECON , SCARC_CO_PRECON
PUBLIC SCARC_MG_OMEGA , SCARC_CG_OMEGA  , SCARC_BICG_OMEGA  , SCARC_SM_OMEGA  , SCARC_CO_OMEGA
PUBLIC SCARC_MG_NLDIFF


CHARACTER (40) :: SCARC_FN                           ! file name for ScaRC debug messages
CHARACTER (10) :: SCARC_METHOD='null'                ! name of method for the solution of the pressure equation

INTEGER        :: SCARC_DEBUG=0,     &               ! debug level (0: no debug messages)
                  SCARC_CASE=0,      &               ! choose different initial solutions
                  SCARC_COUNT=0                      ! counter  for comparison vectors

INTEGER        :: SCARC_LU                           ! unit number for Scarc debug file

CHARACTER (10) :: SCARC_MG_PRECON  ='SSOR', &        ! smoother for mg-method 
                  SCARC_CG_PRECON  ='SSOR', &        ! preconditioner for cg-method         
                  SCARC_SM_PRECON  ='SSOR', &        ! ...                smoother           
                  SCARC_CO_PRECON  ='SSOR', &        ! ...                coarse grid solver 
                  SCARC_BICG_PRECON='SSOR'           ! ...                bicg-method        

INTEGER        :: SCARC_NIT                          ! number of iterations of selected ScaRC variant

INTEGER        :: SCARC_MG_NLDIFF = 2                ! maximum level difference for mg-method

INTEGER        :: SCARC_MG_NIT   =1000, &            ! max number of iterations for mg
                  SCARC_CG_NIT   =1000, &            ! ...                          cg
                  SCARC_SM_NIT   =1000, &            ! ...                          smoother
                  SCARC_CO_NIT   =1000, &            ! ...                          coarse grid solver
                  SCARC_BICG_NIT =1000               ! ...                          bicg

REAL (EB)      :: SCARC_MG_EPS  =1.E-12_EB, &        ! convergence epsilon for mg
                  SCARC_CG_EPS  =1.E-12_EB, &        ! ...                     cg
                  SCARC_SM_EPS  =1.E-12_EB, &        ! ...                     smoother
                  SCARC_CO_EPS  =1.E-12_EB, &        ! ...                     coarse grid solver
                  SCARC_BICG_EPS=1.E-12_EB           ! ...                     bicg

REAL (EB)      :: SCARC_MG_OMEGA  =1.0E+0_EB, &      ! relaxation parameter for mg
                  SCARC_CG_OMEGA  =1.5E+0_EB, &      ! ...                      cg
                  SCARC_SM_OMEGA  =0.9E+0_EB, &      ! ...                      smoother
                  SCARC_CO_OMEGA  =1.0E+0_EB, &      ! ...                      coarse grid solver
                  SCARC_BICG_OMEGA=1.0E+0_EB         ! ...                      bicg

REAL (EB)      :: SCARC_RES                          ! residual of selected ScaRC variant
REAL (EB)      :: SCARC_CAPPA                        ! convergence rate of selected ScaRC variant
REAL (EB)      :: SCARC_EPS_DIVG =1.E+6_EB           ! divergence epsilon for all methods
REAL (EB)      :: SCARC_EPS_REL  =1.E-2_EB           ! minimum relative accuracy for all methods
LOGICAL        :: SCARC_BREL=.FALSE.                 ! relative accuracy activated ? (valid for all methods)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Private structures
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! Private constants
INTEGER, PARAMETER :: NCOM_TYPE_NONE   =  0,  &      ! pure communication along internal boundaries
                      NCOM_TYPE_MATVEC = -1          ! matrix-vector communication along internal boundaries

INTEGER, PARAMETER :: NTYPE_GLOBAL= 1, &             ! (matrix-)vector operations are performed globally
                      NTYPE_LOCAL = 2                ! (matrix-)vector operations are performed locally

INTEGER, PARAMETER :: NCOM_INIT = 0, &               ! initialize communication
                      NCOM_MATV = 1, &               ! matrix-vector communication (only face-based)
                      NCOM_FACE = 1, &               ! only face communication
                      NCOM_FULL = 4                  ! full face-edge-diag communication

INTEGER, PARAMETER :: NMV_NONE= 0, &                 ! communication vector for matrix-vector-multiplication
                      NMV_Y   = 1, &                 ! ...
                      NMV_G   = 2, &                 ! ...
                      NMV_R   = 3, &                 ! ...
                      NMV_D   = 4, &                 ! ...
                      NMV_X   = 5, &                 ! ...
                      NMV_Z   = 6, &                 ! ...
                      NMV_X2  = 7, &                 ! ...
                      NMV_D2  = 8, &                 ! ...
                      NMV_R2  = 9, &                 ! ...
                      NMV_Y2  =10                    ! ...

!!! Private global variables 
LOGICAL :: BCOM_EXTENDED = .FALSE.
INTEGER :: NNLEVEL = 10                              ! max number of levels currently allowed
INTEGER :: SNODE, RNODE, NDIAG
INTEGER :: NFACE
 
INTEGER :: NREQ_FACE                                 ! protocol information for data exchange
INTEGER, ALLOCATABLE, DIMENSION (:)    :: REQ_FACE
INTEGER, ALLOCATABLE, DIMENSION (:, :) :: TAGS_FACE
INTEGER, ALLOCATABLE, DIMENSION (:, :) :: STAT_FACE
INTEGER, ALLOCATABLE, DIMENSION (:, :) :: NBR_FACE

INTEGER, PARAMETER :: N_TIMERS_SCARC=35              ! time measurements within ScaRC
REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: TUSED_SCARC
REAL(EB):: SCARC_TNOW

INTEGER :: NumMaster, NumSlaves                      ! experimental code for including a master process
INTEGER :: myMasterId, mySlaveId, myself              
INTEGER :: MpiGroupWorld, MpiGroupSlaves               
INTEGER :: myMPI_COMM_WORLD, myMpi_COMM_SLAVES         
INTEGER :: NumExclude, RankExclude (4)                
 
 
!!! Private type declarations
!------------------------------------------------------------------------------------------
! Scarc type for multigrid level on 'own' mesh
! in case of the cg-method, only the finest level is used
!------------------------------------------------------------------------------------------
TYPE SCARC_LEVEL_TYPE
 
! local dimensions  on corresponding level
INTEGER :: IBAR, JBAR, KBAR

! number of local and global grid cells on corresponding level
INTEGER :: NCELLS_LOCAL, NCELLS_GLOBAL
INTEGER :: N_EXTERNAL_WALL_CELLS

! mesh coordinates and step sizes
REAL (EB) :: DX,    DY,    DZ
REAL (EB) :: DXI,   DYI,   DZI
REAL (EB) :: DXI2,  DYI2,  DZI2
REAL (EB) :: DXMIN, DYMIN, DZMIN

REAL (EB) :: DX_MIN, DY_MIN, DZ_MIN
REAL (EB) :: DX_MAX, DY_MAX, DZ_MAX

REAL (EB) :: XX_MIN, YY_MIN, ZZ_MIN
REAL (EB) :: XX_MAX, YY_MAX, ZZ_MAX
 
REAL (EB), POINTER, DIMENSION (:) :: XX, YY, ZZ

! communication vectors 
REAL (EB), POINTER, DIMENSION (:) :: SEND_FACE, RECV_FACE

! neighbourship structures 
INTEGER, POINTER, DIMENSION (:, :) :: NIC, IJKW, IOR_FACE
INTEGER, POINTER, DIMENSION (:, :) :: I_MIN, I_MAX, J_MIN, J_MAX, K_MIN, K_MAX

! boundary information
INTEGER, POINTER, DIMENSION (:) :: PRESSURE_BC_INDEX

! matrices and iteration vectors for global defect correction on corresponding level
REAL (EB), POINTER, DIMENSION (:)       :: MDX, LDX, UDX, MDY, MDZ, DAUX, PERIOD
REAL (EB), POINTER, DIMENSION (:, :)    :: AG, AL
REAL (EB), POINTER, DIMENSION (:, :, :) :: X , F , D , Y , G , R, FFT
REAL (EB), POINTER, DIMENSION (:, :, :) :: X2, F2, D2, Y2, G2, R2
REAL (EB), POINTER, DIMENSION (:, :, :) :: Z, TMP
REAL (EB) :: ASUBX, ASUBY, ASUBZ, ASUB
 
! communication vectors
REAL (EB), POINTER, DIMENSION (:, :, :) :: Y_MATV, R_MATV, G_MATV
REAL (EB), POINTER, DIMENSION (:, :, :) :: Y_FACE, R_FACE, G_FACE, Z_FACE
 
! iteration parameters for solver method on corresponding level
REAL (EB) :: RES_CG, RESIN_CG, OMEGA_CG, CAPPA_CG
REAL (EB) :: RES_MG, RESIN_MG, OMEGA_MG, CAPPA_MG
REAL (EB) :: RES_SM, RESIN_SM, OMEGA_SM, CAPPA_SM
REAL (EB) :: RES_CO, RESIN_CO, OMEGA_CO, CAPPA_CO
REAL (EB) :: RES_BICG, RESIN_BICG, OMEGA_BICG, CAPPA_BICG

END TYPE SCARC_LEVEL_TYPE
 
! Scarc type for multigrid level on 'other' mesh
! in case of the cg-method, only the finest level is used
TYPE OSCARC_LEVEL_TYPE
 
! local dimensions  on corresponding level
INTEGER :: IBAR, JBAR, KBAR

! number of local and global grid cells on corresponding level
INTEGER :: NCELLS_LOCAL, NCELLS_GLOBAL
INTEGER :: N_EXTERNAL_WALL_CELLS

! communication vectors
REAL (EB), POINTER, DIMENSION (:) :: SEND_FACE, RECV_FACE

! neighbourship structures along faces, edges and vertices
INTEGER, POINTER, DIMENSION (:, :) :: NIC, IJKW
INTEGER, POINTER, DIMENSION (:, :) :: I_MIN, I_MAX
INTEGER, POINTER, DIMENSION (:, :) :: J_MIN, J_MAX
INTEGER, POINTER, DIMENSION (:, :) :: K_MIN, K_MAX

! iteration vectors for global defect correction on corresponding level
REAL (EB), POINTER, DIMENSION (:, :, :) :: X, Y, F, D, G, R, Z, TMP

! communication vectors
REAL (EB), POINTER, DIMENSION (:, :, :) :: Y_FACE, R_FACE, G_FACE, Z_FACE
 
END TYPE OSCARC_LEVEL_TYPE
 

!------------------------------------------------------------------------------------------
! Scarc type on 'own' mesh
!------------------------------------------------------------------------------------------
TYPE SCARC_TYPE
 
! global MG-levels and MG-cycle variable
INTEGER :: NLEVEL, NLMIN, NLMAX, NLDIFF
INTEGER :: ICYCLE, ITE

! neighborship description arrays
INTEGER, POINTER, DIMENSION (:)         :: MIBAR, MJBAR, MKBAR
INTEGER, POINTER, DIMENSION (:, :)      :: MLEVEL0
INTEGER, POINTER, DIMENSION (:)         :: MLEVEL
INTEGER, POINTER, DIMENSION (:, : )     :: KCYCLE

! neighboring ScaRC-structures
TYPE (OSCARC_TYPE), POINTER, DIMENSION (:) :: OSCARC

! different mesh levels (in case of global multigrid)
INTEGER :: MIBAR_MIN, MJBAR_MIN, MKBAR_MIN, MESH_MIN

TYPE (SCARC_LEVEL_TYPE), POINTER, DIMENSION (:) :: SLEVEL
 
END TYPE SCARC_TYPE


!------------------------------------------------------------------------------------------
! Scarc type on 'other' mesh   
!------------------------------------------------------------------------------------------
TYPE OSCARC_TYPE
TYPE (OSCARC_LEVEL_TYPE), POINTER, DIMENSION (:) :: SLEVEL
END TYPE OSCARC_TYPE


!------------------------------------------------------------------------------------------
! Define necessary types
!------------------------------------------------------------------------------------------
TYPE (SCARC_TYPE),  SAVE, DIMENSION(:), ALLOCATABLE, TARGET :: SCARC
TYPE (OSCARC_TYPE), SAVE, DIMENSION(:), ALLOCATABLE, TARGET :: OSCARC
 
TYPE (SCARC_TYPE),        POINTER ::   S
TYPE (OSCARC_TYPE),       POINTER ::  OS

TYPE (SCARC_LEVEL_TYPE),  POINTER ::  SL,  SLMAX,  SLMIN,  SLHI,  SLLO
TYPE (OSCARC_LEVEL_TYPE), POINTER :: OSL, OSLMAX, OSLMIN, OSLHI, OSLLO

CONTAINS
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! SCARC_INITIALIZE : Initialize Scarc structures for 2D- or 3D-case
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_INITIALIZE (NM)

INTEGER :: IM, NM, IERR, NFINE
INTEGER :: IMIN, JMIN, KMIN
INTEGER :: IBAR0, JBAR0, KBAR0
REAL(EB):: TNOW_INITIALIZE

TNOW_INITIALIZE = SECOND()

IERR=0

!!!
!!! Open debug-file if requested
!!!
IF (SCARC_DEBUG >= 1) THEN
   WRITE (SCARC_FN, '(A,A,i3.3)') TRIM(CHID),'.scarc',NM
   SCARC_LU = GET_FILE_NUMBER()
   OPEN (SCARC_LU, FILE=SCARC_FN)
ENDIF

IF (PRES_METHOD/='SCARC') GOTO 12345

!!!
!!! Allocate SCARC and OSCARC structures
!!!
ALLOCATE (SCARC(NMESHES), STAT=IERR)
CALL CHKMEMERR ('SCARC_INITIALIZE', 'SCARC', IERR)
 
S => SCARC (NM)

ALLOCATE (S%OSCARC(NMESHES), STAT=IERR)
CALL CHKMEMERR ('SCARC_INITIALIZE', 'OSCARC', IERR)


!!!
!!! Initialize time measurement array
!!!
ALLOCATE(TUSED_SCARC(0:N_TIMERS_SCARC,NMESHES),STAT=IERR)
CALL ChkMemErr('SCARC_INITIALIZE','TUSED_SCARC',IERR)

TUSED_SCARC      = 0._EB
TUSED_SCARC(0,:) = SECOND()

!!!
!!! store mesh dimensions on each direction of each mesh within SCARC-structure;
!!! get minimum step size over all meshes for each direction
!!!
ALLOCATE (S%MIBAR(NMESHES), STAT=IERR)
CALL CHKMEMERR ('SCARC_INITIALIZE', 'S%MIBAR', IERR)
 
ALLOCATE (S%MJBAR(NMESHES), STAT=IERR)
CALL CHKMEMERR ('SCARC_INITIALIZE', 'S%MJBAR', IERR)
 
ALLOCATE (S%MKBAR(NMESHES), STAT=IERR)
CALL CHKMEMERR ('SCARC_INITIALIZE', 'S%MKBAR', IERR)

S%MIBAR_MIN=100000
S%MJBAR_MIN=100000
S%MKBAR_MIN=100000

DO IM=1,NMESHES

S%MIBAR(IM)=MESHES(IM)%IBAR
S%MJBAR(IM)=MESHES(IM)%JBAR
S%MKBAR(IM)=MESHES(IM)%KBAR

IF (S%MIBAR(IM)<=S%MIBAR_MIN) THEN
   S%MIBAR_MIN=S%MIBAR(IM)
   IMIN=IM
ENDIF

IF (.NOT.TWO_D) THEN
   IF (S%MJBAR(IM)<=S%MJBAR_MIN) THEN
      S%MJBAR_MIN=S%MJBAR(IM)
      JMIN=IM
   ENDIF
ENDIF

IF (S%MKBAR(IM)<=S%MKBAR_MIN) THEN
   S%MKBAR_MIN=S%MKBAR(IM)
   KMIN=IM
ENDIF

ENDDO


IF (SCARC_DEBUG>=2) THEN
WRITE(SCARC_LU,*) 'S%MIBAR_MIN=',S%MIBAR_MIN
IF (.NOT.TWO_D) WRITE(SCARC_LU,*) 'S%MJBAR_MIN=',S%MJBAR_MIN
WRITE(SCARC_LU,*) 'S%MKBAR_MIN=',S%MKBAR_MIN
ENDIF

!!!
!!! Determine number of grid levels  (1 for CG-method, NLEVEL for MG-method)
!!! and corresponding numbers of cells in each direction
!!!
ALLOCATE (S%MLEVEL0(NMESHES,3), STAT=IERR)
CALL CHKMEMERR ('SCARC_INITIALIZE', 'S%MLEVEL0', IERR)
S%MLEVEL0=0
 
ALLOCATE (S%MLEVEL(NMESHES), STAT=IERR)
CALL CHKMEMERR ('SCARC_INITIALIZE', 'S%MLEVEL', IERR)
S%MLEVEL=0
 
S%NLEVEL=100000
DO IM=1,NMESHES

IBAR0=S%MIBAR(IM)
IF (.NOT.TWO_D) JBAR0=S%MJBAR(IM)
KBAR0=S%MKBAR(IM)

! number of refinements in x-direction
DO NFINE=1,50
   IF (MOD(IBAR0,2)/=0.AND.SCARC_METHOD=='MG') THEN
      WRITE(*,*) 'IBAR=',IBAR0,' NOT YET ALLOWED FOR SCARC-MULTIGRID !'
      EXIT 
   ELSE
      IBAR0=IBAR0/2
      IF (IBAR0==1) EXIT
      IF (NFINE==SCARC_MG_NLDIFF) EXIT
   ENDIF 
ENDDO
S%MLEVEL0(IM,1) = NFINE
S%MLEVEL(IM)    = S%MLEVEL0(IM,1)

IF (S%MLEVEL(IM)<S%NLEVEL) S%NLEVEL = NFINE

IF (SCARC_DEBUG >= 2) THEN
    WRITE(SCARC_LU,*) '1:SCARC_MG_NLDIFF=',SCARC_MG_NLDIFF
    WRITE(SCARC_LU,*) '1:S%MLEVEL(',IM,')=',S%MLEVEL(IM)
    WRITE(SCARC_LU,*) '1:S%NLEVEL     =',S%NLEVEL
ENDIF

! number of refinements in y-direction
IF (.NOT.TWO_D) THEN

   DO NFINE=1,50
      IF (MOD(JBAR0,2)/=0.AND.SCARC_METHOD=='MG') THEN
         WRITE(*,*) 'JBAR=',JBAR0,' NOT YET ALLOWED FOR SCARC-MULTIGRID !'
         EXIT 
      ELSE
         JBAR0=JBAR0/2
         IF (JBAR0==1) EXIT
         IF (NFINE==SCARC_MG_NLDIFF) EXIT
      ENDIF 
   ENDDO
   S%MLEVEL0(IM,2)=NFINE
   IF (S%MLEVEL0(IM,2)<S%MLEVEL(IM)) S%MLEVEL(IM)=S%MLEVEL0(IM,2)

   IF (S%MLEVEL(IM)<S%NLEVEL) S%NLEVEL = NFINE

   IF (SCARC_DEBUG >= 2) THEN
      WRITE(SCARC_LU,*) '2:S%MLEVEL(',IM,')=',S%MLEVEL(IM)
      WRITE(SCARC_LU,*) '2:S%NLEVEL     =',S%NLEVEL
   ENDIF
ELSE
   S%MLEVEL0(IM,2)=1
ENDIF 


! number of refinements in z-direction
DO NFINE=1,50
   IF (MOD(KBAR0,2)/=0.AND.SCARC_METHOD=='MG') THEN
      WRITE(*,*) 'KBAR=',KBAR0,' NOT YET ALLOWED FOR SCARC-MULTIGRID !'
      EXIT 
   ELSE
      KBAR0=KBAR0/2
      IF (KBAR0==1) EXIT
      IF (NFINE==SCARC_MG_NLDIFF) EXIT
   ENDIF 
ENDDO
S%MLEVEL0(IM,3)=NFINE
IF (S%MLEVEL0(IM,3)<S%MLEVEL(IM)) S%MLEVEL(IM)=S%MLEVEL0(IM,3)

IF (S%MLEVEL(IM)<S%NLEVEL) S%NLEVEL =NFINE


IF (SCARC_DEBUG >= 2) then
   WRITE(SCARC_LU,*) '============ MESH ', IM
   WRITE(SCARC_LU,*) '3:S%MLEVEL(',IM,')=',S%MLEVEL(IM)
   WRITE(SCARC_LU,*) '3:S%NLEVEL     =',S%NLEVEL
   WRITE(SCARC_LU,*) 
   WRITE(SCARC_LU,*) 'S%MLEVEL0(',IM,',1)  =',S%MLEVEL0(IM,1)
   WRITE(SCARC_LU,*) 'S%MLEVEL0(',IM,',2)  =',S%MLEVEL0(IM,2)
   WRITE(SCARC_LU,*) 'S%MLEVEL0(',IM,',3)  =',S%MLEVEL0(IM,3)
   WRITE(SCARC_LU,*) 
   WRITE(SCARC_LU,*) 'S%MLEVEL(',IM,')  =',S%MLEVEL(IM)
   WRITE(SCARC_LU,*) 'S%NLEVEL          =',S%NLEVEL   
ENDIF

ENDDO


!!!
!!! Determine number of levels for global method
!!! depends on the minimum number of levels found on all meshes
!!!
 
S%NLMAX  = S%MLEVEL(NM)

! mg-structures with complete grid hierarchy 
IF (SCARC_METHOD == 'MG' .OR. SCARC_CG_PRECON == 'MG'.OR. SCARC_BICG_PRECON=='MG') THEN       
   S%NLMIN  = S%MLEVEL(NM)-S%NLEVEL+1

! only cg-structures on one grid level
ELSE                                                              
   S%NLEVEL = 1
   S%NLMIN  = S%NLMAX
ENDIF

S%NLDIFF = S%NLMAX - S%NLMIN

IF (SCARC_DEBUG >= 2) then
WRITE(SCARC_LU,*) 'S%NLEVEL     =',S%NLEVEL
WRITE(SCARC_LU,*) 'S%NLMAX =',S%NLMAX
WRITE(SCARC_LU,*) 'S%NLMIN =',S%NLMIN
WRITE(SCARC_LU,*) 'S%NLDIFF=',S%NLDIFF
ENDIF

!!!
!!! Allocate corresponding number of SCARC_LEVEL-structures
!!!
ALLOCATE (S%SLEVEL(S%NLMIN:S%NLMAX), STAT=IERR)
CALL CHKMEMERR ('SCARC_INITIALIZE', 'SLEVEL', IERR)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Initialize neighborship structure, matrices, mesh exchange and global solver method
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (TWO_D) THEN

CALL SCARC_INIT_MESHES2D (NM)
CALL SCARC_INIT_NEIGHBORS2D(NM)
CALL SCARC_INIT_MATRICES2D (NM)
CALL SCARC_INIT_COMMUNICATION2D (NM)
CALL SCARC_INIT_SOLVER2D (NM)

ELSE

CALL SCARC_INIT_MESHES3D (NM)         
CALL SCARC_INIT_NEIGHBORS3D (NM)
CALL SCARC_INIT_MATRICES3D (NM)
CALL SCARC_INIT_COMMUNICATION3D (NM)
CALL SCARC_INIT_SOLVER3D (NM)

ENDIF
 
12345 CONTINUE
IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,*) 'Leaving scarc_initialize'

TUSED_SCARC(24,NM)=TUSED_SCARC(24,NM)+SECOND()-TNOW_INITIALIZE
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_INITIALIZE

END SUBROUTINE SCARC_INITIALIZE
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! SCARC_INIT_MESHES2D : initialize grid structures in 2D 
!!!  -  one fine grid on every MESH for CG-method
!!!  -  complete hierarchy of grids on every MESH for MG-method
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_INIT_MESHES2D (NM)
 
INTEGER :: NM, I, K, ILEVEL
INTEGER :: IBAR0, KBAR0, IERR
TYPE (MESH_TYPE), POINTER :: M
REAL(EB):: TNOW_MESHES2D

TNOW_MESHES2D = SECOND()

M => MESHES(NM)

IERR=0

IBAR0=M%IBAR
KBAR0=M%KBAR

GRID_LEVEL_LOOP2D: DO ILEVEL=S%NLMAX,S%NLMIN,-1

SL => S%SLEVEL(ILEVEL)

! ----------------------------------------------------------------------------------
! numbers of cells in x-, y- and z-direction for level 'ILEVEL'
! ----------------------------------------------------------------------------------
SL%IBAR=IBAR0
SL%JBAR=1
SL%KBAR=KBAR0

! ----------------------------------------------------------------------------------
! step widths in x-, y- and z-direction for level 'ILEVEL'
! ----------------------------------------------------------------------------------
SL%DX=(M%XF-M%XS)/REAL(SL%IBAR,EB)
SL%DY=0.0_EB
SL%DZ=(M%ZF-M%ZS)/REAL(SL%KBAR,EB)

SL%DXI=1.0_EB/SL%DX
SL%DYI=1.0_EB
SL%DZI=1.0_EB/SL%DZ

! ----------------------------------------------------------------------------------
! x-, y-, and z-coordinates for level 'ILEVEL'
! -> compute minimum and maximum coordinates in each direction 
! -> compute minimum and maximum step sizes  in each direction
! ----------------------------------------------------------------------------------
ALLOCATE (SL%XX(0:SL%IBAR), STAT=IERR)
CALL CHKMEMERR ('SCARC_INIT_MESHES2D', 'SL%XX', IERR)
SL%XX_MIN= 1.0E+5_EB
SL%XX_MAX=-1.0E+5_EB
SL%DX_MIN= 1.0E+5_EB
SL%DX_MAX=-1.0E+5_EB
DO I=0,SL%IBAR
   SL%XX(I)=M%XS+I*SL%DX
   SL%XX_MIN=MIN(SL%XX_MIN,SL%XX(I))
   SL%XX_MAX=MAX(SL%XX_MAX,SL%XX(I))
   IF (I>0) THEN
     SL%DX_MIN=MIN(SL%DX_MIN,ABS(SL%XX(I)-SL%XX(I-1)))
     SL%DX_MAX=MAX(SL%DX_MAX,ABS(SL%XX(I)-SL%XX(I-1)))
   ENDIF
ENDDO

ALLOCATE (SL%YY(0:0), STAT=IERR)
CALL CHKMEMERR ('SCARC_INIT_MESHES2D', 'SL%YY', IERR)
SL%YY(0) =0.0_EB
SL%YY_MIN=0.0_EB
SL%YY_MAX=0.0_EB
SL%DY_MIN=0.0_EB
SL%DY_MAX=0.0_EB

ALLOCATE (SL%ZZ(0:SL%KBAR), STAT=IERR)
CALL CHKMEMERR ('SCARC_INIT_MESHES2D', 'SL%ZZ', IERR)
SL%ZZ_MIN= 1.0E+5_EB
SL%ZZ_MAX=-1.0E+5_EB
SL%DZ_MIN= 1.0E+5_EB
SL%DZ_MAX=-1.0E+5_EB
DO K=0,SL%KBAR
   SL%ZZ(K)=M%ZS+K*SL%DZ
   SL%ZZ_MIN=MIN(SL%ZZ_MIN,SL%ZZ(K))
   SL%ZZ_MAX=MAX(SL%ZZ_MAX,SL%ZZ(K))
   IF (K>0) THEN
     SL%DZ_MIN=MIN(SL%DX_MIN,ABS(SL%ZZ(K)-SL%ZZ(K-1)))
     SL%DZ_MAX=MAX(SL%DX_MAX,ABS(SL%ZZ(K)-SL%ZZ(K-1)))
   ENDIF
ENDDO

! ----------------------------------------------------------------------------------
! Get global number of grid cells
! ----------------------------------------------------------------------------------
SL%NCELLS_LOCAL = SL%IBAR * SL%KBAR 
IF (NMESHES>1) THEN
   IF (USE_MPI) THEN
     CALL MPI_ALLREDUCE(SL%NCELLS_LOCAL,SL%NCELLS_GLOBAL,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERR)
   ELSE
      WRITE(*,*) 'Serial version not yet implemented'
      stop
   ENDIF
ELSE
   SL%NCELLS_GLOBAL = SL%NCELLS_LOCAL
ENDIF
 
! ----------------------------------------------------------------------------------
! print debug information if requested
! ----------------------------------------------------------------------------------
IF (SCARC_DEBUG>=8) THEN
   WRITE(SCARC_LU,*) ' ======================= SCARC_INIT_MESHES2D =========='
   WRITE(SCARC_LU,*) '==========> LEVEL ', ILEVEL
   WRITE(SCARC_LU,*) 'SL(',ILEVEL,')%IBAR=',SL%IBAR
   WRITE(SCARC_LU,*) 'SL(',ILEVEL,')%JBAR=',SL%JBAR
   WRITE(SCARC_LU,*) 'SL(',ILEVEL,')%KBAR=',SL%KBAR
   WRITE(SCARC_LU,*) 
   WRITE(SCARC_LU,*) 'SL(',ILEVEL,')%DX=',SL%DX
   WRITE(SCARC_LU,*) 'SL(',ILEVEL,')%DY=',SL%DY
   WRITE(SCARC_LU,*) 'SL(',ILEVEL,')%DZ=',SL%DZ
   WRITE(SCARC_LU,*) 
   WRITE(SCARC_LU,*) 'M%XS=',M%XS
   WRITE(SCARC_LU,*) 'M%YS=',M%YS
   WRITE(SCARC_LU,*) 'M%ZS=',M%ZS
   WRITE(SCARC_LU,*) 
   WRITE(SCARC_LU,*) 'SL%XX_MIN=',SL%XX_MIN
   WRITE(SCARC_LU,*) 'SL%XX_MAX=',SL%XX_MAX
   WRITE(SCARC_LU,*) 'SL%DX_MIN=',SL%DX_MIN
   WRITE(SCARC_LU,*) 'SL%DX_MAX=',SL%DX_MAX
   WRITE(SCARC_LU,*) 
   WRITE(SCARC_LU,*) 'SL%YY_MIN=',SL%YY_MIN
   WRITE(SCARC_LU,*) 'SL%YY_MAX=',SL%YY_MAX
   WRITE(SCARC_LU,*) 'SL%DY_MIN=',SL%DY_MIN
   WRITE(SCARC_LU,*) 'SL%DY_MAY=',SL%DY_MAX
   WRITE(SCARC_LU,*) 
   WRITE(SCARC_LU,*) 'SL%ZZ_MIN=',SL%ZZ_MIN
   WRITE(SCARC_LU,*) 'SL%ZZ_MAX=',SL%ZZ_MAX
   WRITE(SCARC_LU,*) 'SL%DZ_MIN=',SL%DZ_MIN
   WRITE(SCARC_LU,*) 'SL%DZ_MAX=',SL%DZ_MAX
   WRITE(SCARC_LU,*) 
   WRITE(SCARC_LU,'(a,i3,a,16f7.3)') 'SL(',ILEVEL,')%XX=',(SL%XX(I),I=0,SL%IBAR)
   WRITE(SCARC_LU,'(a,i3,a,16f7.3)') 'SL(',ILEVEL,')%ZZ=',(SL%ZZ(I),I=0,SL%KBAR)
   WRITE(SCARC_LU,*) 
   WRITE(SCARC_LU,*) 'NMESHES=', NMESHES
   WRITE(SCARC_LU,*) 'NCELLS_LOCAL=', SL%NCELLS_LOCAL
   WRITE(SCARC_LU,*) 'NCELLS_GLOBAL=', SL%NCELLS_GLOBAL
ENDIF

IBAR0=IBAR0/2
KBAR0=KBAR0/2

ENDDO GRID_LEVEL_LOOP2D
 
TUSED_SCARC(1,NM)=TUSED_SCARC(1,NM)+SECOND()-TNOW_MESHES2D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_MESHES2D

END SUBROUTINE SCARC_INIT_MESHES2D
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! SCARC_INIT_MESHES3D : initialize grid structures in 3D 
!!!  -  one fine grid on every MESH for CG-method
!!!  -  complete hierarchy of grids on every MESH for MG-method
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_INIT_MESHES3D (NM)
 
INTEGER :: NM, I, K, ILEVEL
INTEGER :: IBAR0, JBAR0, KBAR0, IERR
TYPE (MESH_TYPE), POINTER :: M
REAL(EB):: TNOW_MESHES3D

TNOW_MESHES3D = SECOND()
 
M => MESHES(NM)

IERR=0

IBAR0=M%IBAR
JBAR0=M%JBAR
KBAR0=M%KBAR

GRID_LEVEL_LOOP3D: DO ILEVEL=S%NLMAX,S%NLMIN,-1

SL => S%SLEVEL(ILEVEL)

! ----------------------------------------------------------------------------------
! numbers of cells in x-, y- and z-direction for level 'ILEVEL'
! ----------------------------------------------------------------------------------
SL%IBAR=IBAR0
SL%JBAR=JBAR0
SL%KBAR=KBAR0

! ----------------------------------------------------------------------------------
! step widths in x-, y- and z-direction for level 'ILEVEL'
! ----------------------------------------------------------------------------------
SL%DX=(M%XF-M%XS)/REAL(SL%IBAR,EB)
SL%DY=(M%YF-M%YS)/REAL(SL%JBAR,EB)
SL%DZ=(M%ZF-M%ZS)/REAL(SL%KBAR,EB)

SL%DXI=1.0_EB/SL%DX
SL%DYI=1.0_EB/SL%DY
SL%DZI=1.0_EB/SL%DZ

! ----------------------------------------------------------------------------------
! x-, y-, and z-coordinates for level 'ILEVEL'
! -> compute minimum and maximum coordinates in each direction 
! -> compute minimum and maximum step sizes  in each direction
! ----------------------------------------------------------------------------------
ALLOCATE (SL%XX(0:SL%IBAR), STAT=IERR)
CALL CHKMEMERR ('SCARC_INIT_MESHES3D', 'SL%XX', IERR)
SL%XX_MIN= 1.0E+5_EB
SL%XX_MAX=-1.0E+5_EB
SL%DX_MIN= 1.0E+5_EB
SL%DX_MAX=-1.0E+5_EB
DO I=0,SL%IBAR
   SL%XX(I)=M%XS+I*SL%DX
   SL%XX_MIN=MIN(SL%XX_MIN,SL%XX(I))
   SL%XX_MAX=MAX(SL%XX_MAX,SL%XX(I))
   IF (I>0) THEN
     SL%DX_MIN=MIN(SL%DX_MIN,ABS(SL%XX(I)-SL%XX(I-1)))
     SL%DX_MAX=MAX(SL%DX_MAX,ABS(SL%XX(I)-SL%XX(I-1)))
   ENDIF
ENDDO

ALLOCATE (SL%YY(0:SL%KBAR), STAT=IERR)
CALL CHKMEMERR ('SCARC_INIT_MESHES3D', 'SL%YY', IERR)
SL%YY_MIN= 1.0E+5_EB
SL%YY_MAX=-1.0E+5_EB
SL%DY_MIN= 1.0E+5_EB
SL%DY_MAX=-1.0E+5_EB
DO K=0,SL%KBAR
   SL%YY(K)=M%YS+K*SL%DY
   SL%YY_MIN=MIN(SL%YY_MIN,SL%YY(K))
   SL%YY_MAX=MAX(SL%YY_MAX,SL%YY(K))
   IF (K>0) THEN
     SL%DY_MIN=MIN(SL%DX_MIN,ABS(SL%YY(K)-SL%YY(K-1)))
     SL%DY_MAX=MAX(SL%DX_MAX,ABS(SL%YY(K)-SL%YY(K-1)))
   ENDIF
ENDDO

ALLOCATE (SL%ZZ(0:SL%KBAR), STAT=IERR)
CALL CHKMEMERR ('SCARC_INIT_MESHES3D', 'SL%ZZ', IERR)
SL%ZZ_MIN= 1.0E+5_EB
SL%ZZ_MAX=-1.0E+5_EB
SL%DZ_MIN= 1.0E+5_EB
SL%DZ_MAX=-1.0E+5_EB
DO K=0,SL%KBAR
   SL%ZZ(K)=M%ZS+K*SL%DZ
   SL%ZZ_MIN=MIN(SL%ZZ_MIN,SL%ZZ(K))
   SL%ZZ_MAX=MAX(SL%ZZ_MAX,SL%ZZ(K))
   IF (K>0) THEN
     SL%DZ_MIN=MIN(SL%DX_MIN,ABS(SL%ZZ(K)-SL%ZZ(K-1)))
     SL%DZ_MAX=MAX(SL%DX_MAX,ABS(SL%ZZ(K)-SL%ZZ(K-1)))
   ENDIF
ENDDO

! ----------------------------------------------------------------------------------
! Get global number of grid cells
! ----------------------------------------------------------------------------------
SL%NCELLS_LOCAL = SL%IBAR * SL%JBAR * SL%KBAR 
IF (NMESHES>1) THEN
   IF (USE_MPI) THEN
      CALL MPI_ALLREDUCE(SL%NCELLS_LOCAL,SL%NCELLS_GLOBAL,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERR)
   ELSE
      WRITE(*,*) 'Serial version not yet implemented'
      stop
   ENDIF
ELSE
   SL%NCELLS_GLOBAL = SL%NCELLS_LOCAL
ENDIF
 
! ----------------------------------------------------------------------------------
! print debug information if requested
! ----------------------------------------------------------------------------------
IF (SCARC_DEBUG.ge.2) THEN
   WRITE(SCARC_LU,*) ' ======================= SCARC_INIT_MESHES3D =========='
   WRITE(SCARC_LU,*) '==========> LEVEL ', ILEVEL
   WRITE(SCARC_LU,*) 'SL(',ILEVEL,')%IBAR=',SL%IBAR
   WRITE(SCARC_LU,*) 'SL(',ILEVEL,')%JBAR=',SL%JBAR
   WRITE(SCARC_LU,*) 'SL(',ILEVEL,')%KBAR=',SL%KBAR
   WRITE(SCARC_LU,*) 
   WRITE(SCARC_LU,*) 'SL(',ILEVEL,')%DX=',SL%DX
   WRITE(SCARC_LU,*) 'SL(',ILEVEL,')%DY=',SL%DY
   WRITE(SCARC_LU,*) 'SL(',ILEVEL,')%DZ=',SL%DZ
   WRITE(SCARC_LU,*) 
   WRITE(SCARC_LU,*) 'M%XS=',M%XS
   WRITE(SCARC_LU,*) 'M%YS=',M%YS
   WRITE(SCARC_LU,*) 'M%ZS=',M%ZS
   WRITE(SCARC_LU,*) 
   WRITE(SCARC_LU,*) 'SL%XX_MIN=',SL%XX_MIN
   WRITE(SCARC_LU,*) 'SL%XX_MAX=',SL%XX_MAX
   WRITE(SCARC_LU,*) 'SL%DX_MIN=',SL%DX_MIN
   WRITE(SCARC_LU,*) 'SL%DX_MAX=',SL%DX_MAX
   WRITE(SCARC_LU,*) 
   WRITE(SCARC_LU,*) 'SL%YY_MIN=',SL%YY_MIN
   WRITE(SCARC_LU,*) 'SL%YY_MAX=',SL%YY_MAX
   WRITE(SCARC_LU,*) 'SL%DY_MIN=',SL%DY_MIN
   WRITE(SCARC_LU,*) 'SL%DY_MAY=',SL%DY_MAX
   WRITE(SCARC_LU,*) 
   WRITE(SCARC_LU,*) 'SL%ZZ_MIN=',SL%ZZ_MIN
   WRITE(SCARC_LU,*) 'SL%ZZ_MAX=',SL%ZZ_MAX
   WRITE(SCARC_LU,*) 'SL%DZ_MIN=',SL%DZ_MIN
   WRITE(SCARC_LU,*) 'SL%DZ_MAX=',SL%DZ_MAX
   WRITE(SCARC_LU,*) 
   !WRITE(SCARC_LU,'(a,i3,a,16f7.3)') 'SL(',ILEVEL,')%XX=',(SL%XX(I),I=0,SL%IBAR)
   !WRITE(SCARC_LU,'(a,i3,a,16f7.3)') 'SL(',ILEVEL,')%YY=',(SL%YY(I),I=0,0)
   !WRITE(SCARC_LU,'(a,i3,a,16f7.3)') 'SL(',ILEVEL,')%ZZ=',(SL%ZZ(I),I=0,SL%KBAR)
   WRITE(SCARC_LU,*) 
   WRITE(SCARC_LU,*) 'NMESHES=', NMESHES
   WRITE(SCARC_LU,*) 'NCELLS_LOCAL=', SL%NCELLS_LOCAL
   WRITE(SCARC_LU,*) 'NCELLS_GLOBAL=', SL%NCELLS_GLOBAL
ENDIF

IBAR0=IBAR0/2
JBAR0=JBAR0/2
KBAR0=KBAR0/2

ENDDO GRID_LEVEL_LOOP3D
 
TUSED_SCARC(1,NM)=TUSED_SCARC(1,NM)+SECOND()-TNOW_MESHES3D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_MESHES3D

END SUBROUTINE SCARC_INIT_MESHES3D
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! determine neighborship/communication structure for data exchange for 2D-case
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_INIT_NEIGHBORS2D(NM)
 
USE GEOMETRY_FUNCTIONS, ONLY: SEARCH_OTHER_MESHES

INTEGER :: NM, NOM, IERR=0, N
INTEGER :: IREFINE, ILEVEL, ILMAX, IOR
INTEGER :: IMIN, IMAX, JMIN, JMAX, KMIN, KMAX
INTEGER :: IMIN1_HI, IMIN2_HI
INTEGER :: IMAX1_HI, IMAX2_HI
INTEGER :: KMIN1_HI, KMIN2_HI
INTEGER :: KMAX1_HI, KMAX2_HI
INTEGER :: IBAR_NOM, KBAR_NOM
INTEGER :: IW1_HI , IW2_HI 
INTEGER :: BC1_HI , BC2_HI 
INTEGER :: NOM1_HI, NOM2_HI
INTEGER :: IW, IW_LO, IM , MMM, NNN, JM, INEWC, IOFFSET, IDIFF, KDIFF
INTEGER :: I, I_LO, J_LO, K_LO
INTEGER, ALLOCATABLE, DIMENSION(:) :: COUNTS,DISPLS,COUNTS2D,DISPLS2D
REAL(EB):: TNOW_NEIGHBORS2D
TYPE (MESH_TYPE), POINTER :: M
 
TNOW_NEIGHBORS2D = SECOND()

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Determine neighborship structures on finest level
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

M => MESHES(NM)

ILMAX =  S%NLMAX
SLMAX => S%SLEVEL(ILMAX)

SLMAX%N_EXTERNAL_WALL_CELLS  =  M%N_EXTERNAL_WALL_CELLS
SLMAX%IJKW  => M%IJKW

ALLOCATE(SLMAX%NIC(NMESHES,NMESHES))
ALLOCATE(SLMAX%I_MIN(NMESHES,NMESHES))
ALLOCATE(SLMAX%I_MAX(NMESHES,NMESHES))
ALLOCATE(SLMAX%J_MIN(NMESHES,NMESHES))
ALLOCATE(SLMAX%J_MAX(NMESHES,NMESHES))
ALLOCATE(SLMAX%K_MIN(NMESHES,NMESHES))
ALLOCATE(SLMAX%K_MAX(NMESHES,NMESHES))

SLMAX%NIC  =NIC
SLMAX%I_MIN=I_MIN
SLMAX%I_MAX=I_MAX
SLMAX%J_MIN=J_MIN
SLMAX%J_MAX=J_MAX
SLMAX%K_MIN=K_MIN
SLMAX%K_MAX=K_MAX

ALLOCATE(SLMAX%IOR_FACE(NMESHES,NMESHES))
SLMAX%IOR_FACE=0

IBAR=SLMAX%IBAR
JBAR=SLMAX%JBAR
KBAR=SLMAX%KBAR

IF (SCARC_DEBUG>=2) THEN
   WRITE(SCARC_LU,*) 'PRESSURE_BC_INDEX:'
   WRITE(SCARC_LU,'(8i4)') (M%PRESSURE_BC_INDEX(IW),IW=1,M%N_EXTERNAL_WALL_CELLS)
   WRITE(SCARC_LU,*) '----------------------------------'
ENDIF

DO IW = 1, M%N_EXTERNAL_WALL_CELLS
   
   IF (ABS(M%IJKW(4,IW))/=2) THEN
      IF (M%BOUNDARY_TYPE(IW)==OPEN_BOUNDARY) THEN
         M%PRESSURE_BC_INDEX(IW)=DIRICHLET
         !M%PRESSURE_BC_INDEX(IW)=M%PRESSURE_BC_INDEX(IW)
      ELSE IF (M%IJKW(9,IW)/=0) THEN
         M%PRESSURE_BC_INDEX(IW)=INTERNAL
      ELSE IF (M%BOUNDARY_TYPE(IW)==NULL_BOUNDARY) THEN
         M%PRESSURE_BC_INDEX(IW)=DIRICHLET
         !M%PRESSURE_BC_INDEX(IW)=M%PRESSURE_BC_INDEX(IW)
      ELSE
         M%PRESSURE_BC_INDEX(IW)=NEUMANN
         !M%PRESSURE_BC_INDEX(IW)=M%PRESSURE_BC_INDEX(IW)
      ENDIF
   ENDIF
!IF (M%IJKW(9,IW)/=0) M%PRESSURE_BC_INDEX(IW)=INTERNAL

   NOM=M%IJKW(9,IW)
   IF (NOM.NE.0) THEN
      SLMAX%IOR_FACE(NM,NOM)= M%IJKW(4,IW)
      SLMAX%IOR_FACE(NOM,NM)=-M%IJKW(4,IW)
      IF (SCARC_DEBUG>=2) THEN
          write(SCARC_LU,*) 'M%IJKW(9,',IW,')=',M%IJKW(9,IW),M%IJKW(4,IW)
          WRITE(SCARC_LU,*) 'SLMAX%IOR_FACE(',NM,',',NOM,')=',SLMAX%IOR_FACE(NM,NOM)
          WRITE(SCARC_LU,*) 'SLMAX%IOR_FACE(',NOM,',',NM,')=',SLMAX%IOR_FACE(NOM,NM)
      ENDIF
  ENDIF

ENDDO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! print some debug messages (only temporarliy)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (SCARC_DEBUG>=2) THEN
WRITE(SCARC_LU,*) '==========================================================================='
WRITE(SCARC_LU,*) '=== SCARC_INIT_NEIGHBORS2D ',NM
WRITE(SCARC_LU,*) '==========================================================================='
WRITE(SCARC_LU,*) 'NMESHES=', NMESHES
WRITE(SCARC_LU,*) 'S%NLEVEL=', S%NLEVEL
WRITE(SCARC_LU,*) 'HIER SCHAUEN, WELCHES GEBRAUCHT WIRD'
WRITE(SCARC_LU,*) 'N_EXTERNAL_WALL_CELLS=', M%N_EXTERNAL_WALL_CELLS
WRITE(SCARC_LU,*) 'NWC=', M%N_EXTERNAL_WALL_CELLS + M%N_INTERNAL_WALL_CELLS
WRITE(SCARC_LU,*) 'BOUNDARY_TYPE:'
WRITE(SCARC_LU,'(8i4)') (M%BOUNDARY_TYPE(IW),IW=1,M%N_EXTERNAL_WALL_CELLS)
WRITE(SCARC_LU,*) '----------------------------------'
WRITE(SCARC_LU,*) 'PRESSURE_BC_INDEX:'
WRITE(SCARC_LU,'(8i4)') (M%PRESSURE_BC_INDEX(IW),IW=1,M%N_EXTERNAL_WALL_CELLS)
WRITE(SCARC_LU,*) '----------------------------------'

WRITE(SCARC_LU,*) 'NIC:'
DO MMM=1,NMESHES
   WRITE(SCARC_LU,'(15i3)') (NIC(NNN,MMM),NNN=1,NMESHES)
ENDDO
WRITE(SCARC_LU,*) 'I_MIN:'
DO MMM=1,NMESHES
   WRITE(SCARC_LU,'(15i3)') (I_MIN(NNN,MMM),NNN=1,NMESHES)
ENDDO
WRITE(SCARC_LU,*) 'I_MAX:'
DO MMM=1,NMESHES
   WRITE(SCARC_LU,'(15i3)') (I_MAX(NNN,MMM),NNN=1,NMESHES)
ENDDO
WRITE(SCARC_LU,*) 'J_MIN:'
DO MMM=1,NMESHES
   WRITE(SCARC_LU,'(15i3)') (J_MIN(NNN,MMM),NNN=1,NMESHES)
ENDDO
WRITE(SCARC_LU,*) 'J_MAX:'
DO MMM=1,NMESHES
   WRITE(SCARC_LU,'(15i3)') (J_MAX(NNN,MMM),NNN=1,NMESHES)
ENDDO
WRITE(SCARC_LU,*) 'K_MIN:'
DO MMM=1,NMESHES
   WRITE(SCARC_LU,'(15i3)') (K_MIN(NNN,MMM),NNN=1,NMESHES)
ENDDO
WRITE(SCARC_LU,*) 'K_MAX:'
DO MMM=1,NMESHES
   WRITE(SCARC_LU,'(15i3)') (K_MAX(NNN,MMM),NNN=1,NMESHES)
ENDDO
WRITE(SCARC_LU,*) 'LEVEL=',S%NLMAX
WRITE(SCARC_LU,*) 'IBAR=',SLMAX%IBAR
WRITE(SCARC_LU,*) 'JBAR=',SLMAX%JBAR
WRITE(SCARC_LU,*) 'KBAR=',SLMAX%KBAR
WRITE(SCARC_LU,*) 'N_EXTERNAL_WALL_CELLS=',SLMAX%N_EXTERNAL_WALL_CELLS
WRITE(SCARC_LU,*) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
WRITE(SCARC_LU,*) "%%% ARRAY IJKW "
WRITE(SCARC_LU,*) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
DO IW=1,M%N_EXTERNAL_WALL_CELLS
   WRITE(SCARC_LU,'(16i5)') IW,(M%IJKW(I,IW),I=1,15)
ENDDO
ENDIF


! communicate NIC-information from other meshes
ALLOCATE(COUNTS(0:NMESHES-1))
ALLOCATE(COUNTS2D(0:NMESHES-1))
ALLOCATE(DISPLS(0:NMESHES-1))
ALLOCATE(DISPLS2D(0:NMESHES-1))
COUNTS=0
DO N=0,NUMPROCS-1
DO IM=1,NMESHES
   IF (PROCESS(IM)==N) COUNTS(N) = COUNTS(N) + 1
ENDDO
ENDDO
DISPLS(0)    = 0
DO N=1,NUMPROCS-1
DISPLS(N) = COUNTS(N-1) + DISPLS(N-1)
ENDDO
COUNTS2D = COUNTS*NMESHES
DISPLS2D = DISPLS*NMESHES


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Only in case of MG-method:
!!! Determine arrays IJKW for coarser levels
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (SCARC_METHOD == 'MG' .OR. SCARC_CG_PRECON =='MG' .OR. SCARC_BICG_PRECON == 'MG') THEN

IREFINE=1
INIT_NBR_LEVEL2D: DO ILEVEL=S%NLMAX-1,S%NLMIN,-1

   IREFINE=IREFINE*2

   SLHI => S%SLEVEL(ILEVEL+1)
   SLLO => S%SLEVEL(ILEVEL)

   SLLO%N_EXTERNAL_WALL_CELLS = 2*SLLO%IBAR + 2*SLLO%KBAR
   ALLOCATE (SLLO%IJKW(15,SLLO%N_EXTERNAL_WALL_CELLS))
   SLLO%IJKW=0

   ALLOCATE (SLLO%PRESSURE_BC_INDEX(SLLO%N_EXTERNAL_WALL_CELLS))
   SLLO%PRESSURE_BC_INDEX=0

   IF (SCARC_DEBUG>=2) THEN
      WRITE(SCARC_LU,*) '========== ILEVEL ',ILEVEL
      WRITE(SCARC_LU,*) 'IREFINE=',IREFINE
      WRITE(SCARC_LU,*) 'SLLO%N_EXTERNAL_WALL_CELLS=',SLLO%N_EXTERNAL_WALL_CELLS
      WRITE(SCARC_LU,*) 'SLLO%IBAR=',SLLO%IBAR 
      WRITE(SCARC_LU,*) 'SLLO%JBAR=',SLLO%JBAR 
      WRITE(SCARC_LU,*) 'SLLO%KBAR=',SLLO%KBAR 
   ENDIF

   IW_LO=1

   !---------------------------------------------------------------------------
   ! cells along IOR=1
   !---------------------------------------------------------------------------
   DO K_LO=1,SLLO%KBAR

      SLLO%IJKW(1,IW_LO)=0           
      SLLO%IJKW(2,IW_LO)=1           ! ghost cell
      SLLO%IJKW(3,IW_LO)=K_LO 
      SLLO%IJKW(4,IW_LO)=1           ! IOR = 1
      SLLO%IJKW(5,IW_LO)=0           ! currently not used
      SLLO%IJKW(6,IW_LO)=1           
      SLLO%IJKW(7,IW_LO)=1           ! boundary cell
      SLLO%IJKW(8,IW_LO)=K_LO 

      IW1_HI=2*K_LO-1
      IW2_HI=2*K_LO

      NOM1_HI=SLHI%IJKW(9,IW1_HI)
      NOM2_HI=SLHI%IJKW(9,IW2_HI)

      IF (NOM1_HI/=NOM2_HI) THEN
         WRITE(*,*) 'Inconsistent neighbors on IOR=1 not allowed!'
         STOP
      ENDIF

      SLLO%IJKW(9,IW_LO)=NOM1_HI 

      IF (ILEVEL==S%NLMAX-1) THEN
         BC1_HI=M%PRESSURE_BC_INDEX(IW1_HI)
         BC2_HI=M%PRESSURE_BC_INDEX(IW2_HI)
      ELSE
         BC1_HI=SLHI%PRESSURE_BC_INDEX(IW1_HI)
         BC2_HI=SLHI%PRESSURE_BC_INDEX(IW2_HI)
      ENDIF

      IF (BC1_HI==INTERNAL.OR.BC2_HI==INTERNAL) THEN
         SLLO%PRESSURE_BC_INDEX(IW_LO)=INTERNAL
      ELSE IF (BC1_HI==DIRICHLET.OR.BC2_HI==DIRICHLET) THEN
         SLLO%PRESSURE_BC_INDEX(IW_LO)=DIRICHLET
      ELSE
         SLLO%PRESSURE_BC_INDEX(IW_LO)=NEUMANN
      ENDIF


      IF (SCARC_DEBUG.ge.2) THEN
         WRITE(SCARC_LU,*) ' 1: =========== K_LO=',K_LO,': J_LO=',J_LO
         WRITE(SCARC_LU,*) 'IW1_HI=',IW1_HI
         WRITE(SCARC_LU,*) 'IW2_HI=',IW2_HI
         WRITE(SCARC_LU,*) 
         WRITE(SCARC_LU,*) 'NOM1_HI=',NOM1_HI
         WRITE(SCARC_LU,*) 'NOM2_HI=',NOM2_HI
         WRITE(SCARC_LU,*) 
         WRITE(SCARC_LU,*) 'BC1_HI=',BC1_HI
         WRITE(SCARC_LU,*) 'BC2_HI=',BC2_HI
         WRITE(SCARC_LU,*) 
         WRITE(SCARC_LU,*) 'SLLO%PRESSURE_BC_INDEX(',IW_LO,')=',SLLO%PRESSURE_BC_INDEX(IW_LO)
      ENDIF

      IF (NOM1_HI > 0) THEN                        ! internal boundary ?

         IBAR_NOM=S%MIBAR(NOM1_HI)/IREFINE
         
         KMIN1_HI=SLHI%IJKW(12,IW1_HI)
         KMIN2_HI=SLHI%IJKW(12,IW2_HI)

         KMAX1_HI=SLHI%IJKW(15,IW1_HI)
         KMAX2_HI=SLHI%IJKW(15,IW2_HI)
        
         KDIFF=KMIN2_HI-KMIN1_HI
   

         ! same JBAR- and KBAR-resolution at neighbor
         IF (KDIFF==1) THEN
            SLLO%IJKW(10,IW_LO)=IBAR_NOM
            SLLO%IJKW(11,IW_LO)=1
            SLLO%IJKW(12,IW_LO)=KMAX2_HI/2
            SLLO%IJKW(13,IW_LO)=IBAR_NOM
            SLLO%IJKW(14,IW_LO)=1
            SLLO%IJKW(15,IW_LO)=KMAX2_HI/2

         ! finer JBAR- and KBAR-resolution at neighbor
         ELSE IF (KDIFF==2) THEN
 
            SLLO%IJKW(10,IW_LO)=IBAR_NOM
            SLLO%IJKW(11,IW_LO)=1
            SLLO%IJKW(12,IW_LO)=KMAX1_HI/2
            SLLO%IJKW(13,IW_LO)=IBAR_NOM
            SLLO%IJKW(14,IW_LO)=1
            SLLO%IJKW(15,IW_LO)=KMAX2_HI/2

         ! coarser JBAR- and KBAR-resolution at neighbor
         ELSE IF (KDIFF==0) THEN
 
            SLLO%IJKW(10,IW_LO)=IBAR_NOM
            SLLO%IJKW(11,IW_LO)=1
            SLLO%IJKW(12,IW_LO)=(KMAX1_HI+1)/2
            SLLO%IJKW(13,IW_LO)=IBAR_NOM
            SLLO%IJKW(14,IW_LO)=1
            SLLO%IJKW(15,IW_LO)=(KMAX1_HI+1)/2

         ELSE                                         
            WRITE(*,*) 'Error1 in INIT_NBR_LEVEL2D at IOR=1 for ILEVEL=',ILEVEL
            WRITE(SCARC_LU,*) 'Error1 in INIT_NBR_LEVEL2D at IOR=1 for ILEVEL=',ILEVEL
            !STOP
         ENDIF
   
         IF (SCARC_DEBUG.ge.2) THEN
            WRITE(SCARC_LU,*) ' 1: =========== K_LO=',K_LO,': J_LO=',J_LO
            WRITE(SCARC_LU,*) 'IW1_HI=',IW1_HI
            WRITE(SCARC_LU,*) 'IW2_HI=',IW2_HI
            WRITE(SCARC_LU,*) 
            WRITE(SCARC_LU,*) 'NOM1_HI=',NOM1_HI
            WRITE(SCARC_LU,*) 'NOM2_HI=',NOM2_HI
            WRITE(SCARC_LU,*) 
            WRITE(SCARC_LU,*) 'IBAR_NOM=',IBAR_NOM
            WRITE(SCARC_LU,*) 
            WRITE(SCARC_LU,*) 'KDIFF=',KDIFF
            WRITE(SCARC_LU,*) 
            WRITE(SCARC_LU,*) 'KMIN1_HI=',KMIN1_HI
            WRITE(SCARC_LU,*) 'KMAX1_HI=',KMAX1_HI
            WRITE(SCARC_LU,*) 
            WRITE(SCARC_LU,*) 'KMIN2_HI=',KMIN2_HI
            WRITE(SCARC_LU,*) 'KMAX2_HI=',KMAX2_HI
            WRITE(SCARC_LU,*) 
            WRITE(SCARC_LU, '(i5,a,15i5)') IW_LO,':',(SLLO%IJKW(I, IW_LO), I=1,15)
         ENDIF
      ENDIF

      IW_LO=IW_LO+1

   ENDDO

   !---------------------------------------------------------------------------
   ! cells along IOR=-1
   !---------------------------------------------------------------------------
   IOFFSET=SLHI%KBAR

   DO K_LO=1,SLLO%KBAR
     
         SLLO%IJKW(1,IW_LO)=SLLO%IBAR+1         
         SLLO%IJKW(2,IW_LO)=1           ! ghost cell
         SLLO%IJKW(3,IW_LO)=K_LO 
         SLLO%IJKW(4,IW_LO)=-1          ! IOR = -1
         SLLO%IJKW(5,IW_LO)=0           ! currently not used
         SLLO%IJKW(6,IW_LO)=SLLO%IBAR           
         SLLO%IJKW(7,IW_LO)=1           ! boundary cell
         SLLO%IJKW(8,IW_LO)=K_LO 
   
         IW1_HI=IOFFSET+2*K_LO-1
         IW2_HI=IOFFSET+2*K_LO

         NOM1_HI=SLHI%IJKW(9,IW1_HI)
         NOM2_HI=SLHI%IJKW(9,IW2_HI)
   
         IF (NOM1_HI/=NOM2_HI) THEN
            WRITE(*,*) 'Inconsistent neighbors on IOR=1 not allowed!'
            STOP
         ENDIF

         SLLO%IJKW(9,IW_LO)=NOM1_HI 

         IF (ILEVEL==S%NLMAX-1) THEN
            BC1_HI=M%PRESSURE_BC_INDEX(IW1_HI)
            BC2_HI=M%PRESSURE_BC_INDEX(IW2_HI)
         ELSE
            BC1_HI=SLHI%PRESSURE_BC_INDEX(IW1_HI)
            BC2_HI=SLHI%PRESSURE_BC_INDEX(IW2_HI)
         ENDIF

         IF (BC1_HI==INTERNAL.OR.BC2_HI==INTERNAL) THEN
            SLLO%PRESSURE_BC_INDEX(IW_LO)=INTERNAL
         ELSE IF (BC1_HI==DIRICHLET.OR.BC2_HI==DIRICHLET) THEN
            SLLO%PRESSURE_BC_INDEX(IW_LO)=DIRICHLET
         ELSE
            SLLO%PRESSURE_BC_INDEX(IW_LO)=NEUMANN
         ENDIF

         IF (SCARC_DEBUG.ge.2) THEN
            WRITE(SCARC_LU,*) ' 1: =========== K_LO=',K_LO,': J_LO=',J_LO
            WRITE(SCARC_LU,*) 'IW1_HI=',IW1_HI
            WRITE(SCARC_LU,*) 'IW2_HI=',IW2_HI
            WRITE(SCARC_LU,*) 
            WRITE(SCARC_LU,*) 'NOM1_HI=',NOM1_HI
            WRITE(SCARC_LU,*) 'NOM2_HI=',NOM2_HI
            WRITE(SCARC_LU,*) 
            WRITE(SCARC_LU,*) 'BC1_HI=',BC1_HI
            WRITE(SCARC_LU,*) 'BC2_HI=',BC2_HI
            WRITE(SCARC_LU,*) 
            WRITE(SCARC_LU,*) 'SLLO%PRESSURE_BC_INDEX(',IW_LO,')=',SLLO%PRESSURE_BC_INDEX(IW_LO)
         ENDIF

         IF (NOM1_HI > 0) THEN                        ! internal boundary ?

            IBAR_NOM=S%MIBAR(NOM1_HI)/IREFINE
            
            KMIN1_HI=SLHI%IJKW(12,IW1_HI)
            KMIN2_HI=SLHI%IJKW(12,IW2_HI)

            KMAX1_HI=SLHI%IJKW(15,IW1_HI)
            KMAX2_HI=SLHI%IJKW(15,IW2_HI)
           
            KDIFF=KMIN2_HI-KMIN1_HI

            ! same JBAR- and KBAR-resolution at neighbor
            IF (KDIFF==1) THEN

               SLLO%IJKW(10,IW_LO)=1
               SLLO%IJKW(11,IW_LO)=1
               SLLO%IJKW(12,IW_LO)=KMAX2_HI/2
               SLLO%IJKW(13,IW_LO)=1
               SLLO%IJKW(14,IW_LO)=1
               SLLO%IJKW(15,IW_LO)=KMAX2_HI/2
   
            ! finer JBAR- and KBAR-resolution at neighbor
            ELSE IF (KDIFF==2) THEN
 
               SLLO%IJKW(10,IW_LO)=1
               SLLO%IJKW(11,IW_LO)=1
               SLLO%IJKW(12,IW_LO)=KMAX1_HI/2
               SLLO%IJKW(13,IW_LO)=1
               SLLO%IJKW(14,IW_LO)=1
               SLLO%IJKW(15,IW_LO)=KMAX2_HI/2
   
            ! coarser JBAR- and KBAR-resolution at neighbor
            ELSE IF (KDIFF==0) THEN
 
               SLLO%IJKW(10,IW_LO)=1
               SLLO%IJKW(11,IW_LO)=1
               SLLO%IJKW(12,IW_LO)=(KMAX1_HI+1)/2
               SLLO%IJKW(13,IW_LO)=1
               SLLO%IJKW(14,IW_LO)=1
               SLLO%IJKW(15,IW_LO)=(KMAX1_HI+1)/2
   
            ELSE                                         
               WRITE(*,*) 'Different JBAR and KBAR not yet implemented!'
               WRITE(SCARC_LU,*) 'Error2 in INIT_NBR_LEVEL2D at IOR=1 for ILEVEL=',ILEVEL
!                  STOP
            ENDIF
      
            IF (SCARC_DEBUG.ge.2) THEN
               WRITE(SCARC_LU,*) ' 1: =========== K_LO=',K_LO,': J_LO=',J_LO
               WRITE(SCARC_LU,*) 'IW1_HI=',IW1_HI
               WRITE(SCARC_LU,*) 'IW2_HI=',IW2_HI
               WRITE(SCARC_LU,*) 
               WRITE(SCARC_LU,*) 'NOM1_HI=',NOM1_HI
               WRITE(SCARC_LU,*) 'NOM2_HI=',NOM2_HI
               WRITE(SCARC_LU,*) 
               WRITE(SCARC_LU,*) 'IBAR_NOM=',IBAR_NOM
               WRITE(SCARC_LU,*) 
               WRITE(SCARC_LU,*) 'KDIFF=',KDIFF
               WRITE(SCARC_LU,*) 
               WRITE(SCARC_LU,*) 'KMIN1_HI=',KMIN1_HI
               WRITE(SCARC_LU,*) 'KMAX1_HI=',KMAX1_HI
               WRITE(SCARC_LU,*) 
               WRITE(SCARC_LU,*) 'KMIN2_HI=',KMIN2_HI
               WRITE(SCARC_LU,*) 'KMAX2_HI=',KMAX2_HI
               WRITE(SCARC_LU,*) 
               WRITE(SCARC_LU, '(i5,a,15i5)') IW_LO,':',(SLLO%IJKW(I, IW_LO), I=1,15)
            ENDIF
         ENDIF

         IW_LO=IW_LO+1
   
      ENDDO
   
      !---------------------------------------------------------------------------
      ! cells along IOR=3
      !---------------------------------------------------------------------------
      IOFFSET=2*SLHI%KBAR + 2*SLHI%IBAR*SLHI%KBAR  ! skip ghost cells in y-direction

      DO I_LO=1,SLLO%IBAR
   
         SLLO%IJKW(1,IW_LO)=I_LO           
         SLLO%IJKW(2,IW_LO)=1              ! ghost cell
         SLLO%IJKW(3,IW_LO)=0 
         SLLO%IJKW(4,IW_LO)=3              ! IOR = 3
         SLLO%IJKW(5,IW_LO)=0              ! currently not used
         SLLO%IJKW(6,IW_LO)=I_LO           
         SLLO%IJKW(7,IW_LO)=1              ! boundary cell
         SLLO%IJKW(8,IW_LO)=1 
   
         IW1_HI=IOFFSET+2*I_LO-1
         IW2_HI=IOFFSET+2*I_LO

         NOM1_HI=SLHI%IJKW(9,IW1_HI)
         NOM2_HI=SLHI%IJKW(9,IW2_HI)

         IF (NOM1_HI/=NOM2_HI) THEN
            WRITE(*,*) 'Inconsistent neighbors on IOR=1 not allowed!'
            STOP
         ENDIF

         SLLO%IJKW(9,IW_LO)=NOM1_HI 

         IF (ILEVEL==S%NLMAX-1) THEN
            BC1_HI=M%PRESSURE_BC_INDEX(IW1_HI)
            BC2_HI=M%PRESSURE_BC_INDEX(IW2_HI)
         ELSE
            BC1_HI=SLHI%PRESSURE_BC_INDEX(IW1_HI)
            BC2_HI=SLHI%PRESSURE_BC_INDEX(IW2_HI)
         ENDIF

         IF (BC1_HI==INTERNAL.OR.BC2_HI==INTERNAL) THEN
            SLLO%PRESSURE_BC_INDEX(IW_LO)=INTERNAL
         ELSE IF (BC1_HI==DIRICHLET.OR.BC2_HI==DIRICHLET) THEN
            SLLO%PRESSURE_BC_INDEX(IW_LO)=DIRICHLET
         ELSE
            SLLO%PRESSURE_BC_INDEX(IW_LO)=NEUMANN
         ENDIF

         IF (SCARC_DEBUG.ge.2) THEN
            WRITE(SCARC_LU,*) ' 1: =========== K_LO=',K_LO,': J_LO=',J_LO
            WRITE(SCARC_LU,*) 'IW1_HI=',IW1_HI
            WRITE(SCARC_LU,*) 'IW2_HI=',IW2_HI
            WRITE(SCARC_LU,*) 
            WRITE(SCARC_LU,*) 'NOM1_HI=',NOM1_HI
            WRITE(SCARC_LU,*) 'NOM2_HI=',NOM2_HI
            WRITE(SCARC_LU,*) 
            WRITE(SCARC_LU,*) 'BC1_HI=',BC1_HI
            WRITE(SCARC_LU,*) 'BC2_HI=',BC2_HI
            WRITE(SCARC_LU,*) 
            WRITE(SCARC_LU,*) 'SLLO%PRESSURE_BC_INDEX(',IW_LO,')=',SLLO%PRESSURE_BC_INDEX(IW_LO)
         ENDIF

         IF (NOM1_HI > 0) THEN                        ! internal boundary ?

            KBAR_NOM=S%MKBAR(NOM1_HI)/IREFINE
            
            IMIN1_HI=SLHI%IJKW(10,IW1_HI)
            IMIN2_HI=SLHI%IJKW(10,IW2_HI)

            IMAX1_HI=SLHI%IJKW(13,IW1_HI)
            IMAX2_HI=SLHI%IJKW(13,IW2_HI)
           
            IDIFF=IMIN2_HI-IMIN1_HI

            ! same IBAR- and JBAR-resolution at neighbor
            IF (IDIFF==1) THEN
               SLLO%IJKW(10,IW_LO)=IMAX2_HI/2
               SLLO%IJKW(11,IW_LO)=1
               SLLO%IJKW(12,IW_LO)=KBAR_NOM
               SLLO%IJKW(13,IW_LO)=IMAX2_HI/2
               SLLO%IJKW(14,IW_LO)=1
               SLLO%IJKW(15,IW_LO)=KBAR_NOM
   
            ! finer JBAR- and JBAR-resolution at neighbor
            ELSE IF (IDIFF==2) THEN
 
               SLLO%IJKW(10,IW_LO)=IMAX1_HI/2
               SLLO%IJKW(11,IW_LO)=1
               SLLO%IJKW(12,IW_LO)=KBAR_NOM
               SLLO%IJKW(13,IW_LO)=IMAX2_HI/2
               SLLO%IJKW(14,IW_LO)=1
               SLLO%IJKW(15,IW_LO)=KBAR_NOM
   
            ! coarser IBAR- and JBAR-resolution at neighbor
            ELSE IF (IDIFF==0) THEN
 
               SLLO%IJKW(10,IW_LO)=(IMAX1_HI+1)/2
               SLLO%IJKW(11,IW_LO)=1
               SLLO%IJKW(12,IW_LO)=KBAR_NOM
               SLLO%IJKW(13,IW_LO)=(IMAX1_HI+1)/2
               SLLO%IJKW(14,IW_LO)=1
               SLLO%IJKW(15,IW_LO)=KBAR_NOM
   
            ELSE                                         
               WRITE(*,*) 'Error5 in INIT_NBR_LEVEL2D at IOR=1 for ILEVEL=',ILEVEL
               WRITE(SCARC_LU,*) 'Error5 in INIT_NBR_LEVEL2D at IOR=1 for ILEVEL=',ILEVEL
               !STOP
            ENDIF
      
            IF (SCARC_DEBUG.ge.2) THEN
               WRITE(SCARC_LU,*) ' 1: =========== K_LO=',K_LO,': J_LO=',J_LO
               WRITE(SCARC_LU,*) 'IW1_HI=',IW1_HI
               WRITE(SCARC_LU,*) 'IW2_HI=',IW2_HI
               WRITE(SCARC_LU,*) 
               WRITE(SCARC_LU,*) 'NOM1_HI=',NOM1_HI
               WRITE(SCARC_LU,*) 'NOM2_HI=',NOM2_HI
               WRITE(SCARC_LU,*) 
               WRITE(SCARC_LU,*) 'KBAR_NOM=',KBAR_NOM
               WRITE(SCARC_LU,*) 
               WRITE(SCARC_LU,*) 'IDIFF=',IDIFF
               WRITE(SCARC_LU,*) 
               WRITE(SCARC_LU,*) 'IMIN1_HI=',IMIN1_HI
               WRITE(SCARC_LU,*) 'IMAX1_HI=',IMAX1_HI
               WRITE(SCARC_LU,*) 
               WRITE(SCARC_LU,*) 'IMIN2_HI=',IMIN2_HI
               WRITE(SCARC_LU,*) 'IMAX2_HI=',IMAX2_HI
               WRITE(SCARC_LU,*) 
               WRITE(SCARC_LU, '(i5,a,15i5)') IW_LO,':',(SLLO%IJKW(I, IW_LO), I=1,15)
            ENDIF
         ENDIF

         IW_LO=IW_LO+1
      
      ENDDO

      !---------------------------------------------------------------------------
      ! cells along IOR=-3
      !---------------------------------------------------------------------------
      IOFFSET=2*SLHI%KBAR + 2*SLHI%IBAR*SLHI%KBAR + SLHI%IBAR

       DO I_LO=1,SLLO%IBAR
   
         SLLO%IJKW(1,IW_LO)=I_LO           
         SLLO%IJKW(2,IW_LO)=1              ! ghost cell
         SLLO%IJKW(3,IW_LO)=SLLO%KBAR+1 
         SLLO%IJKW(4,IW_LO)=-3             ! IOR = -3
         SLLO%IJKW(5,IW_LO)=0              ! currently not used
         SLLO%IJKW(6,IW_LO)=I_LO           
         SLLO%IJKW(7,IW_LO)=1              ! boundary cell
         SLLO%IJKW(8,IW_LO)=SLLO%KBAR 
   
         IW1_HI=IOFFSET+2*I_LO-1
         IW2_HI=IOFFSET+2*I_LO

         NOM1_HI=SLHI%IJKW(9,IW1_HI)
         NOM2_HI=SLHI%IJKW(9,IW2_HI)
   
         IF (NOM1_HI/=NOM2_HI) THEN
            WRITE(*,*) 'Inconsistent neighbors on IOR=1 not allowed!'
            STOP
         ENDIF

         SLLO%IJKW(9,IW_LO)=NOM1_HI 

         IF (ILEVEL==S%NLMAX-1) THEN
            BC1_HI=M%PRESSURE_BC_INDEX(IW1_HI)
            BC2_HI=M%PRESSURE_BC_INDEX(IW2_HI)
         ELSE
            BC1_HI=SLHI%PRESSURE_BC_INDEX(IW1_HI)
            BC2_HI=SLHI%PRESSURE_BC_INDEX(IW2_HI)
         ENDIF

         IF (BC1_HI==INTERNAL.OR.BC2_HI==INTERNAL) THEN
            SLLO%PRESSURE_BC_INDEX(IW_LO)=INTERNAL
         ELSE IF (BC1_HI==DIRICHLET.OR.BC2_HI==DIRICHLET) THEN
            SLLO%PRESSURE_BC_INDEX(IW_LO)=DIRICHLET
         ELSE
            SLLO%PRESSURE_BC_INDEX(IW_LO)=NEUMANN
         ENDIF

         IF (SCARC_DEBUG.ge.2) THEN
            WRITE(SCARC_LU,*) ' 1: =========== K_LO=',K_LO,': J_LO=',J_LO
            WRITE(SCARC_LU,*) 'IW1_HI=',IW1_HI
            WRITE(SCARC_LU,*) 'IW2_HI=',IW2_HI
            WRITE(SCARC_LU,*) 
            WRITE(SCARC_LU,*) 'NOM1_HI=',NOM1_HI
            WRITE(SCARC_LU,*) 'NOM2_HI=',NOM2_HI
            WRITE(SCARC_LU,*) 
            WRITE(SCARC_LU,*) 'BC1_HI=',BC1_HI
            WRITE(SCARC_LU,*) 'BC2_HI=',BC2_HI
            WRITE(SCARC_LU,*) 
            WRITE(SCARC_LU,*) 'SLLO%PRESSURE_BC_INDEX(',IW_LO,')=',SLLO%PRESSURE_BC_INDEX(IW_LO)
         ENDIF

         IF (NOM1_HI > 0) THEN                        ! internal boundary ?

            KBAR_NOM=S%MKBAR(NOM1_HI)/IREFINE
            
            IMIN1_HI=SLHI%IJKW(10,IW1_HI)
            IMIN2_HI=SLHI%IJKW(10,IW2_HI)

            IMAX1_HI=SLHI%IJKW(13,IW1_HI)
            IMAX2_HI=SLHI%IJKW(13,IW2_HI)
           
            IDIFF=IMIN2_HI-IMIN1_HI

         IF (SCARC_DEBUG.ge.2) THEN
            WRITE(SCARC_LU,*) 
            WRITE(SCARC_LU,*) 'IMIN1_HI=',IMIN1_HI
            WRITE(SCARC_LU,*) 'IMIN2_HI=',IMIN2_HI
            WRITE(SCARC_LU,*) 
            WRITE(SCARC_LU,*) 'IDIFF=',IDIFF
         ENDIF

            ! same IBAR- and JBAR-resolution at neighbor
            IF (IDIFF==1) THEN
               SLLO%IJKW(10,IW_LO)=IMAX2_HI/2
               SLLO%IJKW(11,IW_LO)=1
               SLLO%IJKW(12,IW_LO)=1
               SLLO%IJKW(13,IW_LO)=IMAX2_HI/2
               SLLO%IJKW(14,IW_LO)=1
               SLLO%IJKW(15,IW_LO)=1
   
            ! finer JBAR- and JBAR-resolution at neighbor
            ELSE IF (IDIFF==2) THEN
 
               SLLO%IJKW(10,IW_LO)=IMAX1_HI/2
               SLLO%IJKW(11,IW_LO)=1
               SLLO%IJKW(12,IW_LO)=1
               SLLO%IJKW(13,IW_LO)=IMAX2_HI/2
               SLLO%IJKW(14,IW_LO)=1
               SLLO%IJKW(15,IW_LO)=1
   
            ! coarser IBAR- and JBAR-resolution at neighbor
            ELSE IF (IDIFF==0) THEN
 
               SLLO%IJKW(10,IW_LO)=(IMAX2_HI+1)/2
               SLLO%IJKW(11,IW_LO)=1
               SLLO%IJKW(12,IW_LO)=1
               SLLO%IJKW(13,IW_LO)=(IMAX2_HI+1)/2
               SLLO%IJKW(14,IW_LO)=1
               SLLO%IJKW(15,IW_LO)=1
   
            ELSE                                         
               WRITE(*,*) 'Error6 in INIT_NBR_LEVEL2D at IOR=1 for ILEVEL=',ILEVEL
               WRITE(SCARC_LU,*) 'Error6 in INIT_NBR_LEVEL2D at IOR=1 for ILEVEL=',ILEVEL
               !STOP
            ENDIF
      
            IF (SCARC_DEBUG.ge.2) THEN
               WRITE(SCARC_LU,*) ' 1: =========== K_LO=',K_LO,': J_LO=',J_LO
               WRITE(SCARC_LU,*) 'IW1_HI=',IW1_HI
               WRITE(SCARC_LU,*) 'IW2_HI=',IW2_HI
               WRITE(SCARC_LU,*) 
               WRITE(SCARC_LU,*) 'NOM1_HI=',NOM1_HI
               WRITE(SCARC_LU,*) 'NOM2_HI=',NOM2_HI
               WRITE(SCARC_LU,*) 
               WRITE(SCARC_LU,*) 'KBAR_NOM=',KBAR_NOM
               WRITE(SCARC_LU,*) 
               WRITE(SCARC_LU,*) 'IDIFF=',IDIFF
               WRITE(SCARC_LU,*) 
               WRITE(SCARC_LU,*) 'IMIN1_HI=',IMIN1_HI
               WRITE(SCARC_LU,*) 'IMAX1_HI=',IMAX1_HI
               WRITE(SCARC_LU,*) 
               WRITE(SCARC_LU,*) 'IMIN2_HI=',IMIN2_HI
               WRITE(SCARC_LU,*) 'IMAX2_HI=',IMAX2_HI
               WRITE(SCARC_LU,*) 
               WRITE(SCARC_LU, '(i5,a,15i5)') IW_LO,':',(SLLO%IJKW(I, IW_LO), I=1,15)
            ENDIF
         ENDIF

         IW_LO=IW_LO+1
      
      ENDDO

      IF (SCARC_DEBUG >= 6) THEN
         WRITE(SCARC_LU,*) '=================================================='
         WRITE(SCARC_LU,*) '================== LEVEL =',ILEVEL
         WRITE(SCARC_LU,*) 'N_EXTERNAL_WALL_CELLS=',SLLO%N_EXTERNAL_WALL_CELLS
         WRITE(SCARC_LU,*) 'IJKW:'
         DO INEWC=1,SLLO%N_EXTERNAL_WALL_CELLS
            WRITE(SCARC_LU, '(i5,a,15i5)') INEWC,':',(SLLO%IJKW(I, INEWC), I=1,15)
         ENDDO
         WRITE(SCARC_LU,*) 'PRESSURE_BC_INDEX:'
         WRITE(SCARC_LU, '(4i4)') (SLLO%PRESSURE_BC_INDEX(INEWC), INEWC=1,SLLO%N_EXTERNAL_WALL_CELLS)
      ENDIF

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! compute analogues to NIC, I_MIN, I_MAX, K_MIN, K_MAX on level 'ILEVEL'
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ALLOCATE (SLLO%NIC(NMESHES, NMESHES))
      SLLO%NIC=0

      ALLOCATE (SLLO%I_MIN(NMESHES, NMESHES))
      ALLOCATE (SLLO%I_MAX(NMESHES, NMESHES))
      ALLOCATE (SLLO%J_MIN(NMESHES, NMESHES))
      ALLOCATE (SLLO%J_MAX(NMESHES, NMESHES))
      ALLOCATE (SLLO%K_MIN(NMESHES, NMESHES))
      ALLOCATE (SLLO%K_MAX(NMESHES, NMESHES))
      SLLO%I_MIN = 0
      SLLO%I_MAX = 0
      SLLO%J_MIN = 0
      SLLO%J_MAX = 0
      SLLO%K_MIN = 0
      SLLO%K_MAX = 0


      OTHER_MESH_LOOP2D: DO NOM=1,NMESHES
 
         IF (NIC(NM,NOM)==0.AND.NIC(NOM,NM)==0) CYCLE OTHER_MESH_LOOP2D
 
         !!! ACHTUNG: funktioniert nur für 2er-Potenz-Gitterweiten !!!!!
         IMIN=0
         IMAX=S%MIBAR(NOM)/IREFINE+1
         JMIN=0
         JMAX=S%MJBAR(NOM)/IREFINE+1
         KMIN=0
         KMAX=S%MKBAR(NOM)/IREFINE+1

         IF (SCARC_DEBUG.ge.6) THEN
            WRITE(SCARC_LU,*) '============== NOM=',NOM
            WRITE(SCARC_LU,*) 'NIC:'
            DO JM=1,NMESHES
               WRITE(SCARC_LU, '(10i5)') (NIC(IM,JM),IM=1,NMESHES)
            ENDDO
            WRITE(SCARC_LU,*) 'MIBAR(NOM)=',S%MIBAR(NOM)
            WRITE(SCARC_LU,*) 'MJBAR(NOM)=',S%MJBAR(NOM)
            WRITE(SCARC_LU,*) 'MKBAR(NOM)=',S%MKBAR(NOM)
            WRITE(SCARC_LU,*) 'IREFINE=',IREFINE
            WRITE(SCARC_LU,*) 'IMAX=',IMAX
            WRITE(SCARC_LU,*) 'JMAX=',JMAX
            WRITE(SCARC_LU,*) 'KMAX=',KMAX
         ENDIF

         SLLO%NIC(NOM,NM) = 0

         SEARCH_LOOP: DO INEWC=1,SLLO%N_EXTERNAL_WALL_CELLS

            ! neighborship structure already known from finest level
            IF (SLLO%IJKW(9,INEWC)/=NOM) CYCLE SEARCH_LOOP

            !WRITE(SCARC_LU,*) 'IJKW(9,',INEWC,')=',SLLO%IJKW(9,INEWC)

            SLLO%NIC(NOM,NM) = SLLO%NIC(NOM,NM) + 1
            IOR = SLLO%IJKW(4,INEWC)

            !WRITE(SCARC_LU,*) 'NIC(',NOM,',',NM,')=',SLLO%NIC(NOM,NM)

            SELECT CASE(IOR)
               CASE( 1)
                  IMIN=MAX(IMIN,SLLO%IJKW(10,INEWC)-1)
               CASE(-1) 
                  IMAX=MIN(IMAX,SLLO%IJKW(13,INEWC))
               CASE( 2) 
                  JMIN=MAX(JMIN,SLLO%IJKW(11,INEWC)-1)
               CASE(-2) 
                  JMAX=MIN(JMAX,SLLO%IJKW(14,INEWC))
               CASE( 3) 
                  KMIN=MAX(KMIN,SLLO%IJKW(12,INEWC)-1)
               CASE(-3)
                  KMAX=MIN(KMAX,SLLO%IJKW(15,INEWC))
            END SELECT
         ENDDO SEARCH_LOOP

         SLLO%I_MIN(NOM,NM) = IMIN
         SLLO%I_MAX(NOM,NM) = IMAX
         SLLO%J_MIN(NOM,NM) = JMIN
         SLLO%J_MAX(NOM,NM) = JMAX
         SLLO%K_MIN(NOM,NM) = KMIN
         SLLO%K_MAX(NOM,NM) = KMAX

      ENDDO OTHER_MESH_LOOP2D

      IF (SCARC_DEBUG >= 6) THEN
         WRITE(SCARC_LU,*) 'NUMPROCS=',NUMPROCS
         WRITE(SCARC_LU,*) 'MYID  =',MYID  
         WRITE(SCARC_LU,*) 'COUNTS=',COUNTS
         WRITE(SCARC_LU,*) 'COUNTS2D=',COUNTS2D
         WRITE(SCARC_LU,*) 'DISPLS=',DISPLS
         WRITE(SCARC_LU,*) 'DISPLS2D=',DISPLS2D
         WRITE(SCARC_LU,*) 'DISPLS(MYID)+1=',DISPLS(MYID)+1
         WRITE(SCARC_LU,*) 'COUNTS2D(MYID)=',COUNTS2D(MYID)
         WRITE(SCARC_LU,*) 'I_MIN(1,DISPLS(MYID)+1:DISPLS(MYID)+3)=',SLLO%I_MIN(1,DISPLS(MYID)+1:DISPLS(MYID)+3)
         WRITE(SCARC_LU,*) '=======================A:'
         DO JM=1,NMESHES
            WRITE(SCARC_LU, '(10i5)') (SLLO%I_MIN(IM,JM),IM=1,NMESHES)
         ENDDO
      ENDIF

      IF (USE_MPI) THEN
         CALL MPI_ALLGATHERV(SLLO%I_MIN(1,DISPLS(MYID)+1),COUNTS2D(MYID),MPI_INTEGER,&
                             SLLO%I_MIN,COUNTS2D,DISPLS2D,MPI_INTEGER,MPI_COMM_WORLD,IERR)
         CALL MPI_ALLGATHERV(SLLO%I_MAX(1,DISPLS(MYID)+1),COUNTS2D(MYID),MPI_INTEGER,&
                             SLLO%I_MAX,COUNTS2D,DISPLS2D,MPI_INTEGER,MPI_COMM_WORLD,IERR)
         CALL MPI_ALLGATHERV(SLLO%J_MIN(1,DISPLS(MYID)+1),COUNTS2D(MYID),MPI_INTEGER,&
                             SLLO%J_MIN,COUNTS2D,DISPLS2D,MPI_INTEGER,MPI_COMM_WORLD,IERR)
         CALL MPI_ALLGATHERV(SLLO%J_MAX(1,DISPLS(MYID)+1),COUNTS2D(MYID),MPI_INTEGER,&
                             SLLO%J_MAX,COUNTS2D,DISPLS2D,MPI_INTEGER,MPI_COMM_WORLD,IERR)
         CALL MPI_ALLGATHERV(SLLO%K_MIN(1,DISPLS(MYID)+1),COUNTS2D(MYID),MPI_INTEGER,&
                             SLLO%K_MIN,COUNTS2D,DISPLS2D,MPI_INTEGER,MPI_COMM_WORLD,IERR)
         CALL MPI_ALLGATHERV(SLLO%K_MAX(1,DISPLS(MYID)+1),COUNTS2D(MYID),MPI_INTEGER,&
                             SLLO%K_MAX,COUNTS2D,DISPLS2D,MPI_INTEGER,MPI_COMM_WORLD,IERR)
         CALL MPI_ALLGATHERV(SLLO%NIC(1,DISPLS(MYID)+1),  COUNTS2D(MYID),MPI_INTEGER,&
                             SLLO%NIC,  COUNTS2D,DISPLS2D,MPI_INTEGER,MPI_COMM_WORLD,IERR)
      ELSE
         WRITE(*,*) 'Serial version not yet implemented'
         stop
      ENDIF

      SLLO%I_MIN = TRANSPOSE(SLLO%I_MIN)
      SLLO%I_MAX = TRANSPOSE(SLLO%I_MAX)
      SLLO%J_MIN = TRANSPOSE(SLLO%J_MIN)
      SLLO%J_MAX = TRANSPOSE(SLLO%J_MAX)
      SLLO%K_MIN = TRANSPOSE(SLLO%K_MIN)
      SLLO%K_MAX = TRANSPOSE(SLLO%K_MAX)
      SLLO%NIC   = TRANSPOSE(SLLO%NIC)

      IF (SCARC_DEBUG >= 6) THEN
         WRITE(SCARC_LU,*) '=================================================='
         WRITE(SCARC_LU,*) '================== LEVEL =',ILEVEL
         WRITE(SCARC_LU,*) 'N_EXTERNAL_WALL_CELLS=',SLLO%N_EXTERNAL_WALL_CELLS
         WRITE(SCARC_LU,*) 'IJKW:'
         DO INEWC=1,SLLO%N_EXTERNAL_WALL_CELLS
            WRITE(SCARC_LU, '(i5,a,15i5)') INEWC,':',(SLLO%IJKW(I, INEWC), I=1,15)
         ENDDO
         WRITE(SCARC_LU,*) 'NIC:'
         WRITE(SCARC_LU, '(3i4)') ((SLLO%NIC(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
         WRITE(SCARC_LU,*) 'I_MIN:'
         WRITE(SCARC_LU, '(3i4)') ((SLLO%I_MIN(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
         WRITE(SCARC_LU,*) 'I_MAX:'
         WRITE(SCARC_LU, '(3i4)') ((SLLO%I_MAX(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
         WRITE(SCARC_LU,*) 'J_MIN:'
         WRITE(SCARC_LU, '(3i4)') ((SLLO%J_MIN(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
         WRITE(SCARC_LU,*) 'J_MAX:'
         WRITE(SCARC_LU, '(3i4)') ((SLLO%J_MAX(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
         WRITE(SCARC_LU,*) 'K_MIN:'
         WRITE(SCARC_LU, '(3i4)') ((SLLO%K_MIN(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
         WRITE(SCARC_LU,*) 'K_MAX:'
         WRITE(SCARC_LU, '(3i4)') ((SLLO%K_MAX(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
         DO IM=1,NMESHES
            IF (IM/=NM) THEN
               WRITE(SCARC_LU,*) 'FOR MESH ',IM,' ALLOCATING: (',&
                           SLLO%I_MIN(NM,IM),':',SLLO%I_MAX(NM,IM),' , ',&
                           SLLO%J_MIN(NM,IM),':',SLLO%J_MAX(NM,IM),' , ',&
                           SLLO%K_MIN(NM,IM),':',SLLO%K_MAX(NM,IM),')'
            ENDIF
         ENDDO
      ENDIF
       
   ENDDO INIT_NBR_LEVEL2D
ENDIF
    
 
DEALLOCATE(DISPLS2D)
DEALLOCATE(DISPLS)
DEALLOCATE(COUNTS2D)
DEALLOCATE(COUNTS)

TUSED_SCARC(2,NM)=TUSED_SCARC(2,NM)+SECOND()-TNOW_NEIGHBORS2D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_NEIGHBORS2D

END SUBROUTINE SCARC_INIT_NEIGHBORS2D


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! determine neighborship/communication structure for data exchange for 3D-case
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_INIT_NEIGHBORS3D(NM)
 
USE GEOMETRY_FUNCTIONS, ONLY: SEARCH_OTHER_MESHES

INTEGER :: NM, NOM, IERR=0, N
INTEGER :: IREFINE
INTEGER :: ILEVEL
INTEGER :: IW_LO, IOFFSET
INTEGER :: ILMAX, IOR
INTEGER :: IMIN, IMAX, JMIN, JMAX, KMIN, KMAX
INTEGER :: IMIN1_HI, IMIN2_HI, IMIN3_HI, IMIN4_HI
INTEGER :: IMAX1_HI, IMAX2_HI, IMAX3_HI, IMAX4_HI
INTEGER :: JMIN1_HI, JMIN2_HI, JMIN3_HI, JMIN4_HI
INTEGER :: JMAX1_HI, JMAX2_HI, JMAX3_HI, JMAX4_HI
INTEGER :: KMIN1_HI, KMIN2_HI, KMIN3_HI, KMIN4_HI
INTEGER :: KMAX1_HI, KMAX2_HI, KMAX3_HI, KMAX4_HI
INTEGER :: IBAR_NOM, JBAR_NOM, KBAR_NOM
INTEGER :: IDIFF, JDIFF, KDIFF
INTEGER :: IW1_HI , IW2_HI , IW3_HI , IW4_HI
INTEGER :: BC1_HI , BC2_HI , BC3_HI , BC4_HI
INTEGER :: NOM1_HI, NOM2_HI, NOM3_HI, NOM4_HI
INTEGER :: INEWC
INTEGER :: IW, IM , MMM, NNN, JM
INTEGER :: I, I_LO, J_LO, K_LO
INTEGER, ALLOCATABLE, DIMENSION(:) :: COUNTS,DISPLS,COUNTS3D,DISPLS3D
REAL(EB):: TNOW_NEIGHBORS3D
TYPE (MESH_TYPE), POINTER :: M
 
TNOW_NEIGHBORS3D = SECOND()

    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Determine neighborship structures on finest level
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

M => MESHES(NM)

ILMAX =  S%NLMAX
SLMAX => S%SLEVEL(ILMAX)

SLMAX%N_EXTERNAL_WALL_CELLS  =  M%N_EXTERNAL_WALL_CELLS
SLMAX%IJKW  => M%IJKW

ALLOCATE(SLMAX%NIC(NMESHES,NMESHES))
ALLOCATE(SLMAX%I_MIN(NMESHES,NMESHES))
ALLOCATE(SLMAX%I_MAX(NMESHES,NMESHES))
ALLOCATE(SLMAX%J_MIN(NMESHES,NMESHES))
ALLOCATE(SLMAX%J_MAX(NMESHES,NMESHES))
ALLOCATE(SLMAX%K_MIN(NMESHES,NMESHES))
ALLOCATE(SLMAX%K_MAX(NMESHES,NMESHES))

SLMAX%NIC  =NIC
SLMAX%I_MIN=I_MIN
SLMAX%I_MAX=I_MAX
SLMAX%J_MIN=J_MIN
SLMAX%J_MAX=J_MAX
SLMAX%K_MIN=K_MIN
SLMAX%K_MAX=K_MAX

ALLOCATE(SLMAX%IOR_FACE(NMESHES,NMESHES))
SLMAX%IOR_FACE=0

IBAR=SLMAX%IBAR
JBAR=SLMAX%JBAR
KBAR=SLMAX%KBAR

DO IW = 1, M%N_EXTERNAL_WALL_CELLS
   IF (M%BOUNDARY_TYPE(IW)==OPEN_BOUNDARY) THEN
      !M%PRESSURE_BC_INDEX(IW)=DIRICHLET
      M%PRESSURE_BC_INDEX(IW)=M%PRESSURE_BC_INDEX(IW)
   ELSE IF (M%IJKW(9,IW)/=0) THEN
      M%PRESSURE_BC_INDEX(IW)=INTERNAL
   ELSE IF (M%BOUNDARY_TYPE(IW)==NULL_BOUNDARY) THEN
      !M%PRESSURE_BC_INDEX(IW)=DIRICHLET
      M%PRESSURE_BC_INDEX(IW)=M%PRESSURE_BC_INDEX(IW)
   ELSE
      !M%PRESSURE_BC_INDEX(IW)=NEUMANN
      M%PRESSURE_BC_INDEX(IW)=M%PRESSURE_BC_INDEX(IW)
   ENDIF
   !IF (M%IJKW(9,IW)/=0) M%PRESSURE_BC_INDEX(IW)=INTERNAL

   NOM=M%IJKW(9,IW)
   IF (NOM.NE.0) THEN
      SLMAX%IOR_FACE(NM,NOM)= M%IJKW(4,IW)
      SLMAX%IOR_FACE(NOM,NM)=-M%IJKW(4,IW)
      IF (SCARC_DEBUG>=2) THEN
          write(SCARC_LU,*) 'M%IJKW(9,',IW,')=',M%IJKW(9,IW),M%IJKW(4,IW)
          WRITE(SCARC_LU,*) 'SLMAX%IOR_FACE(',NM,',',NOM,')=',SLMAX%IOR_FACE(NM,NOM)
          WRITE(SCARC_LU,*) 'SLMAX%IOR_FACE(',NOM,',',NM,')=',SLMAX%IOR_FACE(NOM,NM)
      ENDIF
  ENDIF

ENDDO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! print some debug messages (only temporarliy)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (SCARC_DEBUG>=2) THEN
   WRITE(SCARC_LU,*) '==========================================================================='
   WRITE(SCARC_LU,*) '=== SCARC_INIT_NEIGHBORS3D ',NM
   WRITE(SCARC_LU,*) '==========================================================================='
   WRITE(SCARC_LU,*) 'NMESHES=', NMESHES
   WRITE(SCARC_LU,*) 'S%NLEVEL=', S%NLEVEL
   WRITE(SCARC_LU,*) 'HIER SCHAUEN, WELCHES GEBRAUCHT WIRD'
   WRITE(SCARC_LU,*) 'N_EXTERNAL_WALL_CELLS=', M%N_EXTERNAL_WALL_CELLS
   WRITE(SCARC_LU,*) 'NWC=', M%N_EXTERNAL_WALL_CELLS + M%N_INTERNAL_WALL_CELLS
   IF (NMESHES==1) THEN
   WRITE(SCARC_LU,*) 'BOUNDARY_TYPE:'
   WRITE(SCARC_LU,'(8i4)') (M%BOUNDARY_TYPE(IW),IW=1,64)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,'(8i4)') (M%BOUNDARY_TYPE(IW),IW=65,128)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,'(8i4)') (M%BOUNDARY_TYPE(IW),IW=129,192)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,'(8i4)') (M%BOUNDARY_TYPE(IW),IW=193,256)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,'(8i4)') (M%BOUNDARY_TYPE(IW),IW=257,320)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,'(8i4)') (M%BOUNDARY_TYPE(IW),IW=321,384)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,*) 'PRESSURE_BC_INDEX:'
   WRITE(SCARC_LU,'(8i4)') (M%PRESSURE_BC_INDEX(IW),IW=1,64)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,'(8i4)') (M%PRESSURE_BC_INDEX(IW),IW=65,128)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,'(8i4)') (M%PRESSURE_BC_INDEX(IW),IW=129,192)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,'(8i4)') (M%PRESSURE_BC_INDEX(IW),IW=193,256)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,'(8i4)') (M%PRESSURE_BC_INDEX(IW),IW=257,320)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,'(8i4)') (M%PRESSURE_BC_INDEX(IW),IW=321,384)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,*) '----------------------------------'
   ELSE
   WRITE(SCARC_LU,*) 'BOUNDARY_TYPE:'
   WRITE(SCARC_LU,'(8i4)') (M%BOUNDARY_TYPE(IW),IW=1,64)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,'(8i4)') (M%BOUNDARY_TYPE(IW),IW=65,128)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,'(4i4)') (M%BOUNDARY_TYPE(IW),IW=129,160)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,'(4i4)') (M%BOUNDARY_TYPE(IW),IW=161,192)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,'(4i4)') (M%BOUNDARY_TYPE(IW),IW=193,224)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,'(4i4)') (M%BOUNDARY_TYPE(IW),IW=225,256)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,*) 'PRESSURE_BC_INDEX:'
   WRITE(SCARC_LU,'(8i4)') (M%PRESSURE_BC_INDEX(IW),IW=1,64)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,'(8i4)') (M%PRESSURE_BC_INDEX(IW),IW=65,128)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,'(4i4)') (M%PRESSURE_BC_INDEX(IW),IW=129,160)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,'(4i4)') (M%PRESSURE_BC_INDEX(IW),IW=161,192)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,'(4i4)') (M%PRESSURE_BC_INDEX(IW),IW=193,224)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,'(4i4)') (M%PRESSURE_BC_INDEX(IW),IW=225,256)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,*) '----------------------------------'
   ENDIF

   WRITE(SCARC_LU,*) 'NIC:'
   DO MMM=1,NMESHES
      WRITE(SCARC_LU,'(15i3)') (NIC(NNN,MMM),NNN=1,NMESHES)
   ENDDO
   WRITE(SCARC_LU,*) 'I_MIN:'
   DO MMM=1,NMESHES
      WRITE(SCARC_LU,'(15i3)') (I_MIN(NNN,MMM),NNN=1,NMESHES)
   ENDDO
   WRITE(SCARC_LU,*) 'I_MAX:'
   DO MMM=1,NMESHES
      WRITE(SCARC_LU,'(15i3)') (I_MAX(NNN,MMM),NNN=1,NMESHES)
   ENDDO
   WRITE(SCARC_LU,*) 'J_MIN:'
   DO MMM=1,NMESHES
      WRITE(SCARC_LU,'(15i3)') (J_MIN(NNN,MMM),NNN=1,NMESHES)
   ENDDO
   WRITE(SCARC_LU,*) 'J_MAX:'
   DO MMM=1,NMESHES
      WRITE(SCARC_LU,'(15i3)') (J_MAX(NNN,MMM),NNN=1,NMESHES)
   ENDDO
   WRITE(SCARC_LU,*) 'K_MIN:'
   DO MMM=1,NMESHES
      WRITE(SCARC_LU,'(15i3)') (K_MIN(NNN,MMM),NNN=1,NMESHES)
   ENDDO
   WRITE(SCARC_LU,*) 'K_MAX:'
   DO MMM=1,NMESHES
      WRITE(SCARC_LU,'(15i3)') (K_MAX(NNN,MMM),NNN=1,NMESHES)
   ENDDO
   WRITE(SCARC_LU,*) 'LEVEL=',S%NLMAX
   WRITE(SCARC_LU,*) 'IBAR=',SLMAX%IBAR
   WRITE(SCARC_LU,*) 'JBAR=',SLMAX%JBAR
   WRITE(SCARC_LU,*) 'KBAR=',SLMAX%KBAR
   WRITE(SCARC_LU,*) 'N_EXTERNAL_WALL_CELLS=',SLMAX%N_EXTERNAL_WALL_CELLS
   WRITE(SCARC_LU,*) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
   WRITE(SCARC_LU,*) "%%% ARRAY IJKW after all INIT_WALL_CELLS in init.f90"
   WRITE(SCARC_LU,*) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
   WRITE(SCARC_LU,'(16i5)') SLMAX%N_EXTERNAL_WALL_CELLS, 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15
   WRITE(SCARC_LU,*) '-------------------------------------------------------------------'
   !DO IW=1,M%N_EXTERNAL_WALL_CELLS
   IF (NM==0) THEN
   WRITE(SCARC_LU,*) '-------------------------------------------------------------------'
   WRITE(SCARC_LU,*) ' IOR =1'
   WRITE(SCARC_LU,*) '-------------------------------------------------------------------'
   DO IW=1,256
      WRITE(SCARC_LU,'(16i5)') IW,(M%IJKW(I,IW),I=1,15)
   ENDDO
   WRITE(SCARC_LU,*) '-------------------------------------------------------------------'
   WRITE(SCARC_LU,*) ' IOR =-1'
   WRITE(SCARC_LU,*) '-------------------------------------------------------------------'
   DO IW=257,512
      WRITE(SCARC_LU,'(16i5)') IW,(M%IJKW(I,IW),I=1,15)
   ENDDO
   WRITE(SCARC_LU,*) '-------------------------------------------------------------------'
   WRITE(SCARC_LU,*) ' IOR =2'
   WRITE(SCARC_LU,*) '-------------------------------------------------------------------'
   DO IW=513,640
      WRITE(SCARC_LU,'(16i5)') IW,(M%IJKW(I,IW),I=1,15)
   ENDDO
   WRITE(SCARC_LU,*) '-------------------------------------------------------------------'
   WRITE(SCARC_LU,*) ' IOR =-2'
   WRITE(SCARC_LU,*) '-------------------------------------------------------------------'
   DO IW=641,768
      WRITE(SCARC_LU,'(16i5)') IW,(M%IJKW(I,IW),I=1,15)
   ENDDO
   WRITE(SCARC_LU,*) '-------------------------------------------------------------------'
   WRITE(SCARC_LU,*) ' IOR =3'
   WRITE(SCARC_LU,*) '-------------------------------------------------------------------'
   DO IW=769,896
      WRITE(SCARC_LU,'(16i5)') IW,(M%IJKW(I,IW),I=1,15)
   ENDDO
   WRITE(SCARC_LU,*) '-------------------------------------------------------------------'
   WRITE(SCARC_LU,*) ' IOR =-3'
   WRITE(SCARC_LU,*) '-------------------------------------------------------------------'
   DO IW=897,1024
      WRITE(SCARC_LU,'(16i5)') IW,(M%IJKW(I,IW),I=1,15)
   ENDDO
   WRITE(SCARC_LU,*) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
   !ELSE IF (NM==2) THEN
   ELSE 
   WRITE(SCARC_LU,*) '-------------------------------------------------------------------'
   WRITE(SCARC_LU,*) ' IOR =1'
   WRITE(SCARC_LU,*) '-------------------------------------------------------------------'
   DO IW=1,64
      WRITE(SCARC_LU,'(16i5)') IW,(M%IJKW(I,IW),I=1,15)
   ENDDO
   WRITE(SCARC_LU,*) '-------------------------------------------------------------------'
   WRITE(SCARC_LU,*) ' IOR =-1'
   WRITE(SCARC_LU,*) '-------------------------------------------------------------------'
   DO IW=65,128
      WRITE(SCARC_LU,'(16i5)') IW,(M%IJKW(I,IW),I=1,15)
   ENDDO
   WRITE(SCARC_LU,*) '-------------------------------------------------------------------'
   WRITE(SCARC_LU,*) ' IOR =2'
   WRITE(SCARC_LU,*) '-------------------------------------------------------------------'
   DO IW=129,160
      WRITE(SCARC_LU,'(16i5)') IW,(M%IJKW(I,IW),I=1,15)
   ENDDO
   WRITE(SCARC_LU,*) '-------------------------------------------------------------------'
   WRITE(SCARC_LU,*) ' IOR =-2'
   WRITE(SCARC_LU,*) '-------------------------------------------------------------------'
   DO IW=161,192
      WRITE(SCARC_LU,'(16i5)') IW,(M%IJKW(I,IW),I=1,15)
   ENDDO
   WRITE(SCARC_LU,*) '-------------------------------------------------------------------'
   WRITE(SCARC_LU,*) ' IOR =3'
   WRITE(SCARC_LU,*) '-------------------------------------------------------------------'
   DO IW=193,224
      WRITE(SCARC_LU,'(16i5)') IW,(M%IJKW(I,IW),I=1,15)
   ENDDO
   WRITE(SCARC_LU,*) '-------------------------------------------------------------------'
   WRITE(SCARC_LU,*) ' IOR =-3'
   WRITE(SCARC_LU,*) '-------------------------------------------------------------------'
   DO IW=224,256
      WRITE(SCARC_LU,'(16i5)') IW,(M%IJKW(I,IW),I=1,15)
   ENDDO
   ENDIF
ENDIF
   

! communicate NIC-information from other meshes
ALLOCATE(COUNTS(0:NMESHES-1))
ALLOCATE(COUNTS3D(0:NMESHES-1))
ALLOCATE(DISPLS(0:NMESHES-1))
ALLOCATE(DISPLS3D(0:NMESHES-1))
COUNTS=0
DO N=0,NUMPROCS-1
   DO IM=1,NMESHES
      IF (PROCESS(IM)==N) COUNTS(N) = COUNTS(N) + 1
   ENDDO
ENDDO
DISPLS(0)    = 0
DO N=1,NUMPROCS-1
   DISPLS(N) = COUNTS(N-1) + DISPLS(N-1)
ENDDO
COUNTS3D = COUNTS*NMESHES
DISPLS3D = DISPLS*NMESHES


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Only in case of MG-method:
!!! Determine arrays IJKW for coarser levels
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (SCARC_METHOD == 'MG' .OR. SCARC_CG_PRECON =='MG' .OR. SCARC_BICG_PRECON == 'MG') THEN

   IREFINE=1
   INIT_NBR_LEVEL3D: DO ILEVEL=S%NLMAX-1,S%NLMIN,-1
   
      IREFINE=IREFINE*2

      SLHI => S%SLEVEL(ILEVEL+1)
      SLLO => S%SLEVEL(ILEVEL)
   
      SLLO%N_EXTERNAL_WALL_CELLS = 2*SLLO%IBAR*SLLO%JBAR + 2*SLLO%IBAR*SLLO%KBAR + 2*SLLO%JBAR*SLLO%KBAR
      ALLOCATE (SLLO%IJKW(15,SLLO%N_EXTERNAL_WALL_CELLS))
      SLLO%IJKW=0

      ALLOCATE (SLLO%PRESSURE_BC_INDEX(SLLO%N_EXTERNAL_WALL_CELLS))
      SLLO%PRESSURE_BC_INDEX=0
   
      IF (SCARC_DEBUG>=2) THEN
         WRITE(SCARC_LU,*) '========== ILEVEL ',ILEVEL
         WRITE(SCARC_LU,*) 'IREFINE=',IREFINE
         WRITE(SCARC_LU,*) 'SLLO%N_EXTERNAL_WALL_CELLS=',SLLO%N_EXTERNAL_WALL_CELLS
         WRITE(SCARC_LU,*) 'SLLO%IBAR=',SLLO%IBAR 
         WRITE(SCARC_LU,*) 'SLLO%JBAR=',SLLO%JBAR 
         WRITE(SCARC_LU,*) 'SLLO%KBAR=',SLLO%KBAR 
      ENDIF
   
      IW_LO=1
   
      !---------------------------------------------------------------------------
      ! cells along IOR=1
      !---------------------------------------------------------------------------
      DO K_LO=1,SLLO%KBAR
        DO J_LO=1,SLLO%JBAR
   
            SLLO%IJKW(1,IW_LO)=0           
            SLLO%IJKW(2,IW_LO)=J_LO        ! ghost cell
            SLLO%IJKW(3,IW_LO)=K_LO 
            SLLO%IJKW(4,IW_LO)=1           ! IOR = 1
            SLLO%IJKW(5,IW_LO)=0           ! currently not used
            SLLO%IJKW(6,IW_LO)=1           
            SLLO%IJKW(7,IW_LO)=J_LO        ! boundary cell
            SLLO%IJKW(8,IW_LO)=K_LO 
      
            IW1_HI=(2*K_LO-2)*SLHI%JBAR+2*J_LO-1
            IW2_HI=(2*K_LO-2)*SLHI%JBAR+2*J_LO
            IW3_HI=(2*K_LO-1)*SLHI%JBAR+2*J_LO-1
            IW4_HI=(2*K_LO-1)*SLHI%JBAR+2*J_LO

            NOM1_HI=SLHI%IJKW(9,IW1_HI)
            NOM2_HI=SLHI%IJKW(9,IW2_HI)
            NOM3_HI=SLHI%IJKW(9,IW3_HI)
            NOM4_HI=SLHI%IJKW(9,IW4_HI)
      
            IF (NOM1_HI/=NOM2_HI .OR. NOM1_HI/=NOM3_HI .OR. NOM1_HI/=NOM4_HI) THEN
               WRITE(*,*) 'Inconsistent neighbors on IOR=1 not allowed!'
               STOP
            ENDIF

            SLLO%IJKW(9,IW_LO)=NOM1_HI 

            IF (ILEVEL==S%NLMAX-1) THEN
               BC1_HI=M%PRESSURE_BC_INDEX(IW1_HI)
               BC2_HI=M%PRESSURE_BC_INDEX(IW2_HI)
               BC3_HI=M%PRESSURE_BC_INDEX(IW3_HI)
               BC4_HI=M%PRESSURE_BC_INDEX(IW4_HI)
            ELSE
               BC1_HI=SLHI%PRESSURE_BC_INDEX(IW1_HI)
               BC2_HI=SLHI%PRESSURE_BC_INDEX(IW2_HI)
               BC3_HI=SLHI%PRESSURE_BC_INDEX(IW3_HI)
               BC4_HI=SLHI%PRESSURE_BC_INDEX(IW4_HI)
            ENDIF

            IF (BC1_HI==INTERNAL.OR.BC2_HI==INTERNAL.OR.BC3_HI==INTERNAL.OR.BC4_HI==INTERNAL) THEN
               SLLO%PRESSURE_BC_INDEX(IW_LO)=INTERNAL
            ELSE IF (BC1_HI==DIRICHLET.OR.BC2_HI==DIRICHLET.OR.BC3_HI==DIRICHLET.OR.BC4_HI==DIRICHLET) THEN
               SLLO%PRESSURE_BC_INDEX(IW_LO)=DIRICHLET
            ELSE
               SLLO%PRESSURE_BC_INDEX(IW_LO)=NEUMANN
            ENDIF


            IF (SCARC_DEBUG.ge.2) THEN
               WRITE(SCARC_LU,*) ' 1: =========== K_LO=',K_LO,': J_LO=',J_LO
               WRITE(SCARC_LU,*) 'IW1_HI=',IW1_HI
               WRITE(SCARC_LU,*) 'IW2_HI=',IW2_HI
               WRITE(SCARC_LU,*) 'IW3_HI=',IW3_HI
               WRITE(SCARC_LU,*) 'IW4_HI=',IW4_HI
               WRITE(SCARC_LU,*) 
               WRITE(SCARC_LU,*) 'NOM1_HI=',NOM1_HI
               WRITE(SCARC_LU,*) 'NOM2_HI=',NOM2_HI
               WRITE(SCARC_LU,*) 'NOM3_HI=',NOM3_HI
               WRITE(SCARC_LU,*) 'NOM4_HI=',NOM4_HI
               WRITE(SCARC_LU,*) 
               WRITE(SCARC_LU,*) 'BC1_HI=',BC1_HI
               WRITE(SCARC_LU,*) 'BC2_HI=',BC2_HI
               WRITE(SCARC_LU,*) 'BC3_HI=',BC3_HI
               WRITE(SCARC_LU,*) 'BC4_HI=',BC4_HI
               WRITE(SCARC_LU,*) 
               WRITE(SCARC_LU,*) 'SLLO%PRESSURE_BC_INDEX(',IW_LO,')=',SLLO%PRESSURE_BC_INDEX(IW_LO)
            ENDIF

            IF (NOM1_HI > 0) THEN                        ! internal boundary ?

               IBAR_NOM=S%MIBAR(NOM1_HI)/IREFINE
               
               JMIN1_HI=SLHI%IJKW(11,IW1_HI)
               JMIN2_HI=SLHI%IJKW(11,IW2_HI)
               JMIN3_HI=SLHI%IJKW(11,IW3_HI)
               JMIN4_HI=SLHI%IJKW(11,IW4_HI)

               JMAX1_HI=SLHI%IJKW(14,IW1_HI)
               JMAX2_HI=SLHI%IJKW(14,IW2_HI)
               JMAX3_HI=SLHI%IJKW(14,IW3_HI)
               JMAX4_HI=SLHI%IJKW(14,IW4_HI)
              
               KMIN1_HI=SLHI%IJKW(12,IW1_HI)
               KMIN2_HI=SLHI%IJKW(12,IW2_HI)
               KMIN3_HI=SLHI%IJKW(12,IW3_HI)
               KMIN4_HI=SLHI%IJKW(12,IW4_HI)

               KMAX1_HI=SLHI%IJKW(15,IW1_HI)
               KMAX2_HI=SLHI%IJKW(15,IW2_HI)
               KMAX3_HI=SLHI%IJKW(15,IW3_HI)
               KMAX4_HI=SLHI%IJKW(15,IW4_HI)
              
               JDIFF=JMIN2_HI-JMIN1_HI
               KDIFF=KMIN3_HI-KMIN1_HI
         

               ! same JBAR- and KBAR-resolution at neighbor
               IF (JDIFF==1.AND.KDIFF==1) THEN
                  SLLO%IJKW(10,IW_LO)=IBAR_NOM
                  SLLO%IJKW(11,IW_LO)=JMAX2_HI/2
                  SLLO%IJKW(12,IW_LO)=KMAX3_HI/2
                  SLLO%IJKW(13,IW_LO)=IBAR_NOM
                  SLLO%IJKW(14,IW_LO)=JMAX2_HI/2
                  SLLO%IJKW(15,IW_LO)=KMAX3_HI/2
      
               ! finer JBAR- and KBAR-resolution at neighbor
               ELSE IF (JDIFF==2.AND.KDIFF==2) THEN
 
                  SLLO%IJKW(10,IW_LO)=IBAR_NOM
                  SLLO%IJKW(11,IW_LO)=JMAX1_HI/2
                  SLLO%IJKW(12,IW_LO)=KMAX1_HI/2
                  SLLO%IJKW(13,IW_LO)=IBAR_NOM
                  SLLO%IJKW(14,IW_LO)=JMAX2_HI/2
                  SLLO%IJKW(15,IW_LO)=KMAX3_HI/2
      
               ! coarser JBAR- and KBAR-resolution at neighbor
               ELSE IF (JDIFF==0.AND.KDIFF==0) THEN
 
                  SLLO%IJKW(10,IW_LO)=IBAR_NOM
                  SLLO%IJKW(11,IW_LO)=(JMAX1_HI+1)/2
                  SLLO%IJKW(12,IW_LO)=(KMAX1_HI+1)/2
                  SLLO%IJKW(13,IW_LO)=IBAR_NOM
                  SLLO%IJKW(14,IW_LO)=(JMAX1_HI+1)/2
                  SLLO%IJKW(15,IW_LO)=(KMAX1_HI+1)/2
      
               ELSE                                         
                  WRITE(*,*) 'Error1 in INIT_NBR_LEVEL3D at IOR=1 for ILEVEL=',ILEVEL
                  WRITE(SCARC_LU,*) 'Error1 in INIT_NBR_LEVEL3D at IOR=1 for ILEVEL=',ILEVEL
                  !STOP
               ENDIF
         
               IF (SCARC_DEBUG.ge.2) THEN
                  WRITE(SCARC_LU,*) ' 1: =========== IOR= 1'
                  WRITE(SCARC_LU,*) ' 1: =========== K_LO=',K_LO,': J_LO=',J_LO
                  WRITE(SCARC_LU,*) 'IW1_HI=',IW1_HI
                  WRITE(SCARC_LU,*) 'IW2_HI=',IW2_HI
                  WRITE(SCARC_LU,*) 'IW3_HI=',IW3_HI
                  WRITE(SCARC_LU,*) 'IW4_HI=',IW4_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'NOM1_HI=',NOM1_HI
                  WRITE(SCARC_LU,*) 'NOM2_HI=',NOM2_HI
                  WRITE(SCARC_LU,*) 'NOM3_HI=',NOM3_HI
                  WRITE(SCARC_LU,*) 'NOM4_HI=',NOM4_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'IBAR_NOM=',IBAR_NOM
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'JDIFF=',JDIFF
                  WRITE(SCARC_LU,*) 'KDIFF=',KDIFF
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'JMIN1_HI=',JMIN1_HI
                  WRITE(SCARC_LU,*) 'JMAX1_HI=',JMAX1_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'JMIN2_HI=',JMIN2_HI
                  WRITE(SCARC_LU,*) 'JMAX2_HI=',JMAX2_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'JMIN3_HI=',JMIN3_HI
                  WRITE(SCARC_LU,*) 'JMAX3_HI=',JMAX3_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'JMIN4_HI=',JMIN4_HI
                  WRITE(SCARC_LU,*) 'JMAX4_HI=',JMAX4_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'KMIN1_HI=',KMIN1_HI
                  WRITE(SCARC_LU,*) 'KMAX1_HI=',KMAX1_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'KMIN2_HI=',KMIN2_HI
                  WRITE(SCARC_LU,*) 'KMAX2_HI=',KMAX2_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'KMIN3_HI=',KMIN3_HI
                  WRITE(SCARC_LU,*) 'KMAX3_HI=',KMAX3_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'KMIN4_HI=',KMIN4_HI
                  WRITE(SCARC_LU,*) 'KMAX4_HI=',KMAX4_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU, '(i5,a,15i5)') IW_LO,':',(SLLO%IJKW(I, IW_LO), I=1,15)
               ENDIF
            ENDIF
   
            IW_LO=IW_LO+1
      
         ENDDO
      ENDDO

      !---------------------------------------------------------------------------
      ! cells along IOR=-1
      !---------------------------------------------------------------------------
      IOFFSET=SLHI%JBAR*SLHI%KBAR

      DO K_LO=1,SLLO%KBAR
        DO J_LO=1,SLLO%JBAR
        
   
            SLLO%IJKW(1,IW_LO)=SLLO%IBAR+1         
            SLLO%IJKW(2,IW_LO)=J_LO        ! ghost cell
            SLLO%IJKW(3,IW_LO)=K_LO 
            SLLO%IJKW(4,IW_LO)=-1          ! IOR = -1
            SLLO%IJKW(5,IW_LO)=0           ! currently not used
            SLLO%IJKW(6,IW_LO)=SLLO%IBAR           
            SLLO%IJKW(7,IW_LO)=J_LO        ! boundary cell
            SLLO%IJKW(8,IW_LO)=K_LO 
      
            IW1_HI=IOFFSET+(2*K_LO-2)*SLHI%JBAR+2*J_LO-1
            IW2_HI=IOFFSET+(2*K_LO-2)*SLHI%JBAR+2*J_LO
            IW3_HI=IOFFSET+(2*K_LO-1)*SLHI%JBAR+2*J_LO-1
            IW4_HI=IOFFSET+(2*K_LO-1)*SLHI%JBAR+2*J_LO

            NOM1_HI=SLHI%IJKW(9,IW1_HI)
            NOM2_HI=SLHI%IJKW(9,IW2_HI)
            NOM3_HI=SLHI%IJKW(9,IW3_HI)
            NOM4_HI=SLHI%IJKW(9,IW4_HI)
      
            IF (NOM1_HI/=NOM2_HI .OR. NOM1_HI/=NOM3_HI .OR. NOM1_HI/=NOM4_HI) THEN
               WRITE(*,*) 'Inconsistent neighbors on IOR=1 not allowed!'
               STOP
            ENDIF

            SLLO%IJKW(9,IW_LO)=NOM1_HI 

            IF (ILEVEL==S%NLMAX-1) THEN
               BC1_HI=M%PRESSURE_BC_INDEX(IW1_HI)
               BC2_HI=M%PRESSURE_BC_INDEX(IW2_HI)
               BC3_HI=M%PRESSURE_BC_INDEX(IW3_HI)
               BC4_HI=M%PRESSURE_BC_INDEX(IW4_HI)
            ELSE
               BC1_HI=SLHI%PRESSURE_BC_INDEX(IW1_HI)
               BC2_HI=SLHI%PRESSURE_BC_INDEX(IW2_HI)
               BC3_HI=SLHI%PRESSURE_BC_INDEX(IW3_HI)
               BC4_HI=SLHI%PRESSURE_BC_INDEX(IW4_HI)
            ENDIF

            IF (BC1_HI==INTERNAL.OR.BC2_HI==INTERNAL.OR.BC3_HI==INTERNAL.OR.BC4_HI==INTERNAL) THEN
               SLLO%PRESSURE_BC_INDEX(IW_LO)=INTERNAL
            ELSE IF (BC1_HI==DIRICHLET.OR.BC2_HI==DIRICHLET.OR.BC3_HI==DIRICHLET.OR.BC4_HI==DIRICHLET) THEN
               SLLO%PRESSURE_BC_INDEX(IW_LO)=DIRICHLET
            ELSE
               SLLO%PRESSURE_BC_INDEX(IW_LO)=NEUMANN
            ENDIF

            IF (SCARC_DEBUG.ge.2) THEN
               WRITE(SCARC_LU,*) ' 1: =========== IOR= -1'
               WRITE(SCARC_LU,*) ' 1: =========== K_LO=',K_LO,': J_LO=',J_LO
               WRITE(SCARC_LU,*) 'IW1_HI=',IW1_HI
               WRITE(SCARC_LU,*) 'IW2_HI=',IW2_HI
               WRITE(SCARC_LU,*) 'IW3_HI=',IW3_HI
               WRITE(SCARC_LU,*) 'IW4_HI=',IW4_HI
               WRITE(SCARC_LU,*) 
               WRITE(SCARC_LU,*) 'NOM1_HI=',NOM1_HI
               WRITE(SCARC_LU,*) 'NOM2_HI=',NOM2_HI
               WRITE(SCARC_LU,*) 'NOM3_HI=',NOM3_HI
               WRITE(SCARC_LU,*) 'NOM4_HI=',NOM4_HI
               WRITE(SCARC_LU,*) 
               WRITE(SCARC_LU,*) 'BC1_HI=',BC1_HI
               WRITE(SCARC_LU,*) 'BC2_HI=',BC2_HI
               WRITE(SCARC_LU,*) 'BC3_HI=',BC3_HI
               WRITE(SCARC_LU,*) 'BC4_HI=',BC4_HI
               WRITE(SCARC_LU,*) 
               WRITE(SCARC_LU,*) 'SLLO%PRESSURE_BC_INDEX(',IW_LO,')=',SLLO%PRESSURE_BC_INDEX(IW_LO)
            ENDIF

            IF (NOM1_HI > 0) THEN                        ! internal boundary ?

               IBAR_NOM=S%MIBAR(NOM1_HI)/IREFINE
               
               JMIN1_HI=SLHI%IJKW(11,IW1_HI)
               JMIN2_HI=SLHI%IJKW(11,IW2_HI)
               JMIN3_HI=SLHI%IJKW(11,IW3_HI)
               JMIN4_HI=SLHI%IJKW(11,IW4_HI)

               JMAX1_HI=SLHI%IJKW(14,IW1_HI)
               JMAX2_HI=SLHI%IJKW(14,IW2_HI)
               JMAX3_HI=SLHI%IJKW(14,IW3_HI)
               JMAX4_HI=SLHI%IJKW(14,IW4_HI)
              
               KMIN1_HI=SLHI%IJKW(12,IW1_HI)
               KMIN2_HI=SLHI%IJKW(12,IW2_HI)
               KMIN3_HI=SLHI%IJKW(12,IW3_HI)
               KMIN4_HI=SLHI%IJKW(12,IW4_HI)

               KMAX1_HI=SLHI%IJKW(15,IW1_HI)
               KMAX2_HI=SLHI%IJKW(15,IW2_HI)
               KMAX3_HI=SLHI%IJKW(15,IW3_HI)
               KMAX4_HI=SLHI%IJKW(15,IW4_HI)
              
               JDIFF=JMIN2_HI-JMIN1_HI
               KDIFF=KMIN3_HI-KMIN1_HI
         

               ! same JBAR- and KBAR-resolution at neighbor
               IF (JDIFF==1.AND.KDIFF==1) THEN

                  SLLO%IJKW(10,IW_LO)=1
                  SLLO%IJKW(11,IW_LO)=JMAX2_HI/2
                  SLLO%IJKW(12,IW_LO)=KMAX3_HI/2
                  SLLO%IJKW(13,IW_LO)=1
                  SLLO%IJKW(14,IW_LO)=JMAX2_HI/2
                  SLLO%IJKW(15,IW_LO)=KMAX3_HI/2
      
               ! finer JBAR- and KBAR-resolution at neighbor
               ELSE IF (JDIFF==2.AND.KDIFF==2) THEN
 
                  SLLO%IJKW(10,IW_LO)=1
                  SLLO%IJKW(11,IW_LO)=JMAX1_HI/2
                  SLLO%IJKW(12,IW_LO)=KMAX1_HI/2
                  SLLO%IJKW(13,IW_LO)=1
                  SLLO%IJKW(14,IW_LO)=JMAX2_HI/2
                  SLLO%IJKW(15,IW_LO)=KMAX3_HI/2
      
               ! coarser JBAR- and KBAR-resolution at neighbor
               ELSE IF (JDIFF==0.AND.KDIFF==0) THEN
 
                  SLLO%IJKW(10,IW_LO)=1
                  SLLO%IJKW(11,IW_LO)=(JMAX1_HI+1)/2
                  SLLO%IJKW(12,IW_LO)=(KMAX1_HI+1)/2
                  SLLO%IJKW(13,IW_LO)=1
                  SLLO%IJKW(14,IW_LO)=(JMAX1_HI+1)/2
                  SLLO%IJKW(15,IW_LO)=(KMAX1_HI+1)/2
      
               ELSE                                         
                  WRITE(*,*) 'Different JBAR and KBAR not yet implemented!'
                  WRITE(SCARC_LU,*) 'Error2 in INIT_NBR_LEVEL3D at IOR=1 for ILEVEL=',ILEVEL
!                  STOP
               ENDIF
         
               IF (SCARC_DEBUG.ge.2) THEN
                  WRITE(SCARC_LU,*) ' 1: =========== IOR= -1'
                  WRITE(SCARC_LU,*) ' 1: =========== K_LO=',K_LO,': J_LO=',J_LO
                  WRITE(SCARC_LU,*) 'IW1_HI=',IW1_HI
                  WRITE(SCARC_LU,*) 'IW2_HI=',IW2_HI
                  WRITE(SCARC_LU,*) 'IW3_HI=',IW3_HI
                  WRITE(SCARC_LU,*) 'IW4_HI=',IW4_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'NOM1_HI=',NOM1_HI
                  WRITE(SCARC_LU,*) 'NOM2_HI=',NOM2_HI
                  WRITE(SCARC_LU,*) 'NOM3_HI=',NOM3_HI
                  WRITE(SCARC_LU,*) 'NOM4_HI=',NOM4_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'IBAR_NOM=',IBAR_NOM
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'JDIFF=',JDIFF
                  WRITE(SCARC_LU,*) 'KDIFF=',KDIFF
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'JMIN1_HI=',JMIN1_HI
                  WRITE(SCARC_LU,*) 'JMAX1_HI=',JMAX1_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'JMIN2_HI=',JMIN2_HI
                  WRITE(SCARC_LU,*) 'JMAX2_HI=',JMAX2_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'JMIN3_HI=',JMIN3_HI
                  WRITE(SCARC_LU,*) 'JMAX3_HI=',JMAX3_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'JMIN4_HI=',JMIN4_HI
                  WRITE(SCARC_LU,*) 'JMAX4_HI=',JMAX4_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'KMIN1_HI=',KMIN1_HI
                  WRITE(SCARC_LU,*) 'KMAX1_HI=',KMAX1_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'KMIN2_HI=',KMIN2_HI
                  WRITE(SCARC_LU,*) 'KMAX2_HI=',KMAX2_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'KMIN3_HI=',KMIN3_HI
                  WRITE(SCARC_LU,*) 'KMAX3_HI=',KMAX3_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'KMIN4_HI=',KMIN4_HI
                  WRITE(SCARC_LU,*) 'KMAX4_HI=',KMAX4_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU, '(i5,a,15i5)') IW_LO,':',(SLLO%IJKW(I, IW_LO), I=1,15)
               ENDIF
            ENDIF
   
            IW_LO=IW_LO+1
      
         ENDDO
      ENDDO
   
      !---------------------------------------------------------------------------
      ! cells along IOR=2
      !---------------------------------------------------------------------------
      IOFFSET=2*SLHI%JBAR*SLHI%KBAR

      DO K_LO=1,SLLO%KBAR
        DO I_LO=1,SLLO%IBAR
   
            SLLO%IJKW(1,IW_LO)=I_LO           
            SLLO%IJKW(2,IW_LO)=0           ! ghost cell
            SLLO%IJKW(3,IW_LO)=K_LO 
            SLLO%IJKW(4,IW_LO)=2           ! IOR = 2
            SLLO%IJKW(5,IW_LO)=0           ! currently not used
            SLLO%IJKW(6,IW_LO)=I_LO           
            SLLO%IJKW(7,IW_LO)=1           ! boundary cell
            SLLO%IJKW(8,IW_LO)=K_LO 
      
            IW1_HI=IOFFSET+(2*K_LO-2)*SLHI%IBAR+2*I_LO-1
            IW2_HI=IOFFSET+(2*K_LO-2)*SLHI%IBAR+2*I_LO
            IW3_HI=IOFFSET+(2*K_LO-1)*SLHI%IBAR+2*I_LO-1
            IW4_HI=IOFFSET+(2*K_LO-1)*SLHI%IBAR+2*I_LO

            NOM1_HI=SLHI%IJKW(9,IW1_HI)
            NOM2_HI=SLHI%IJKW(9,IW2_HI)
            NOM3_HI=SLHI%IJKW(9,IW3_HI)
            NOM4_HI=SLHI%IJKW(9,IW4_HI)
      
            IF (NOM1_HI/=NOM2_HI .OR. NOM1_HI/=NOM3_HI .OR. NOM1_HI/=NOM4_HI) THEN
               WRITE(*,*) 'Inconsistent neighbors on IOR=1 not allowed!'
               STOP
            ENDIF

            SLLO%IJKW(9,IW_LO)=NOM1_HI 

            IF (ILEVEL==S%NLMAX-1) THEN
               BC1_HI=M%PRESSURE_BC_INDEX(IW1_HI)
               BC2_HI=M%PRESSURE_BC_INDEX(IW2_HI)
               BC3_HI=M%PRESSURE_BC_INDEX(IW3_HI)
               BC4_HI=M%PRESSURE_BC_INDEX(IW4_HI)
            ELSE
               BC1_HI=SLHI%PRESSURE_BC_INDEX(IW1_HI)
               BC2_HI=SLHI%PRESSURE_BC_INDEX(IW2_HI)
               BC3_HI=SLHI%PRESSURE_BC_INDEX(IW3_HI)
               BC4_HI=SLHI%PRESSURE_BC_INDEX(IW4_HI)
            ENDIF

            IF (BC1_HI==INTERNAL.OR.BC2_HI==INTERNAL.OR.BC3_HI==INTERNAL.OR.BC4_HI==INTERNAL) THEN
               SLLO%PRESSURE_BC_INDEX(IW_LO)=INTERNAL
            ELSE IF (BC1_HI==DIRICHLET.OR.BC2_HI==DIRICHLET.OR.BC3_HI==DIRICHLET.OR.BC4_HI==DIRICHLET) THEN
               SLLO%PRESSURE_BC_INDEX(IW_LO)=DIRICHLET
            ELSE
               SLLO%PRESSURE_BC_INDEX(IW_LO)=NEUMANN
            ENDIF

            IF (SCARC_DEBUG.ge.2) THEN
               WRITE(SCARC_LU,*) ' 1: =========== IOR= 2'
               WRITE(SCARC_LU,*) ' 1: =========== K_LO=',K_LO,': J_LO=',J_LO
               WRITE(SCARC_LU,*) 'IW1_HI=',IW1_HI
               WRITE(SCARC_LU,*) 'IW2_HI=',IW2_HI
               WRITE(SCARC_LU,*) 'IW3_HI=',IW3_HI
               WRITE(SCARC_LU,*) 'IW4_HI=',IW4_HI
               WRITE(SCARC_LU,*) 
               WRITE(SCARC_LU,*) 'NOM1_HI=',NOM1_HI
               WRITE(SCARC_LU,*) 'NOM2_HI=',NOM2_HI
               WRITE(SCARC_LU,*) 'NOM3_HI=',NOM3_HI
               WRITE(SCARC_LU,*) 'NOM4_HI=',NOM4_HI
               WRITE(SCARC_LU,*) 
               WRITE(SCARC_LU,*) 'BC1_HI=',BC1_HI
               WRITE(SCARC_LU,*) 'BC2_HI=',BC2_HI
               WRITE(SCARC_LU,*) 'BC3_HI=',BC3_HI
               WRITE(SCARC_LU,*) 'BC4_HI=',BC4_HI
               WRITE(SCARC_LU,*) 
               WRITE(SCARC_LU,*) 'SLLO%PRESSURE_BC_INDEX(',IW_LO,')=',SLLO%PRESSURE_BC_INDEX(IW_LO)
            ENDIF

            IF (NOM1_HI > 0) THEN                        ! internal boundary ?

               JBAR_NOM=S%MJBAR(NOM1_HI)/IREFINE
               
               IMIN1_HI=SLHI%IJKW(10,IW1_HI)
               IMIN2_HI=SLHI%IJKW(10,IW2_HI)
               IMIN3_HI=SLHI%IJKW(10,IW3_HI)
               IMIN4_HI=SLHI%IJKW(10,IW4_HI)

               IMAX1_HI=SLHI%IJKW(13,IW1_HI)
               IMAX2_HI=SLHI%IJKW(13,IW2_HI)
               IMAX3_HI=SLHI%IJKW(13,IW3_HI)
               IMAX4_HI=SLHI%IJKW(13,IW4_HI)
              
               KMIN1_HI=SLHI%IJKW(12,IW1_HI)
               KMIN2_HI=SLHI%IJKW(12,IW2_HI)
               KMIN3_HI=SLHI%IJKW(12,IW3_HI)
               KMIN4_HI=SLHI%IJKW(12,IW4_HI)

               KMAX1_HI=SLHI%IJKW(15,IW1_HI)
               KMAX2_HI=SLHI%IJKW(15,IW2_HI)
               KMAX3_HI=SLHI%IJKW(15,IW3_HI)
               KMAX4_HI=SLHI%IJKW(15,IW4_HI)
              
               IDIFF=IMIN2_HI-IMIN1_HI
               KDIFF=KMIN3_HI-KMIN1_HI
         

               ! same JBAR- and KBAR-resolution at neighbor
               IF (IDIFF==1.AND.KDIFF==1) THEN
                  SLLO%IJKW(10,IW_LO)=IMAX2_HI/2
                  SLLO%IJKW(11,IW_LO)=JBAR_NOM
                  SLLO%IJKW(12,IW_LO)=KMAX3_HI/2
                  SLLO%IJKW(13,IW_LO)=IMAX2_HI/2
                  SLLO%IJKW(14,IW_LO)=JBAR_NOM
                  SLLO%IJKW(15,IW_LO)=KMAX3_HI/2
      
               ! finer JBAR- and KBAR-resolution at neighbor
               ELSE IF (IDIFF==2.AND.KDIFF==2) THEN
 
                  SLLO%IJKW(10,IW_LO)=IMAX1_HI/2
                  SLLO%IJKW(11,IW_LO)=JBAR_NOM
                  SLLO%IJKW(12,IW_LO)=KMAX1_HI/2
                  SLLO%IJKW(13,IW_LO)=IMAX2_HI/2
                  SLLO%IJKW(14,IW_LO)=JBAR_NOM
                  SLLO%IJKW(15,IW_LO)=KMAX3_HI/2
      
               ! coarser JBAR- and KBAR-resolution at neighbor
               ELSE IF (IDIFF==0.AND.KDIFF==0) THEN
 
                  SLLO%IJKW(10,IW_LO)=(IMAX1_HI+1)/2
                  SLLO%IJKW(11,IW_LO)=JBAR_NOM
                  SLLO%IJKW(12,IW_LO)=(KMAX1_HI+1)/2
                  SLLO%IJKW(13,IW_LO)=(IMAX1_HI+1)/2
                  SLLO%IJKW(14,IW_LO)=JBAR_NOM
                  SLLO%IJKW(15,IW_LO)=(KMAX1_HI+1)/2
      
               ELSE                                         
                  WRITE(*,*) 'Error3 in INIT_NBR_LEVEL3D at IOR=1 for ILEVEL=',ILEVEL
                  WRITE(SCARC_LU,*) 'Error3 in INIT_NBR_LEVEL3D at IOR=1 for ILEVEL=',ILEVEL
                  !STOP
               ENDIF
         
               IF (SCARC_DEBUG.ge.2) THEN
                  WRITE(SCARC_LU,*) ' 1: =========== IOR= 2'
                  WRITE(SCARC_LU,*) ' 1: =========== K_LO=',K_LO,': J_LO=',J_LO
                  WRITE(SCARC_LU,*) 'IW1_HI=',IW1_HI
                  WRITE(SCARC_LU,*) 'IW2_HI=',IW2_HI
                  WRITE(SCARC_LU,*) 'IW3_HI=',IW3_HI
                  WRITE(SCARC_LU,*) 'IW4_HI=',IW4_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'NOM1_HI=',NOM1_HI
                  WRITE(SCARC_LU,*) 'NOM2_HI=',NOM2_HI
                  WRITE(SCARC_LU,*) 'NOM3_HI=',NOM3_HI
                  WRITE(SCARC_LU,*) 'NOM4_HI=',NOM4_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'JBAR_NOM=',JBAR_NOM
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'IDIFF=',IDIFF
                  WRITE(SCARC_LU,*) 'KDIFF=',KDIFF
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'IMIN1_HI=',IMIN1_HI
                  WRITE(SCARC_LU,*) 'IMAX1_HI=',IMAX1_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'IMIN2_HI=',IMIN2_HI
                  WRITE(SCARC_LU,*) 'IMAX2_HI=',IMAX2_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'IMIN3_HI=',IMIN3_HI
                  WRITE(SCARC_LU,*) 'IMAX3_HI=',IMAX3_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'IMIN4_HI=',IMIN4_HI
                  WRITE(SCARC_LU,*) 'IMAX4_HI=',IMAX4_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'KMIN1_HI=',KMIN1_HI
                  WRITE(SCARC_LU,*) 'KMAX1_HI=',KMAX1_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'KMIN2_HI=',KMIN2_HI
                  WRITE(SCARC_LU,*) 'KMAX2_HI=',KMAX2_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'KMIN3_HI=',KMIN3_HI
                  WRITE(SCARC_LU,*) 'KMAX3_HI=',KMAX3_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'KMIN4_HI=',KMIN4_HI
                  WRITE(SCARC_LU,*) 'KMAX4_HI=',KMAX4_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU, '(i5,a,15i5)') IW_LO,':',(SLLO%IJKW(I, IW_LO), I=1,15)
               ENDIF
            ENDIF
   
            IW_LO=IW_LO+1
      
         ENDDO
      ENDDO

      !---------------------------------------------------------------------------
      ! cells along IOR=-2
      !---------------------------------------------------------------------------
      IOFFSET=2*SLHI%JBAR*SLHI%KBAR+SLHI%IBAR*SLHI%KBAR

      DO K_LO=1,SLLO%KBAR
        DO I_LO=1,SLLO%IBAR
   
            SLLO%IJKW(1,IW_LO)=I_LO           
            SLLO%IJKW(2,IW_LO)=SLLO%JBAR+1       ! ghost cell
            SLLO%IJKW(3,IW_LO)=K_LO 
            SLLO%IJKW(4,IW_LO)=-2                ! IOR = -2
            SLLO%IJKW(5,IW_LO)=0                 ! currently not used
            SLLO%IJKW(6,IW_LO)=I_LO           
            SLLO%IJKW(7,IW_LO)=SLLO%JBAR         ! boundary cell
            SLLO%IJKW(8,IW_LO)=K_LO 
      
            IW1_HI=IOFFSET+(2*K_LO-2)*SLHI%IBAR+2*I_LO-1
            IW2_HI=IOFFSET+(2*K_LO-2)*SLHI%IBAR+2*I_LO
            IW3_HI=IOFFSET+(2*K_LO-1)*SLHI%IBAR+2*I_LO-1
            IW4_HI=IOFFSET+(2*K_LO-1)*SLHI%IBAR+2*I_LO

            NOM1_HI=SLHI%IJKW(9,IW1_HI)
            NOM2_HI=SLHI%IJKW(9,IW2_HI)
            NOM3_HI=SLHI%IJKW(9,IW3_HI)
            NOM4_HI=SLHI%IJKW(9,IW4_HI)
      
            IF (NOM1_HI/=NOM2_HI .OR. NOM1_HI/=NOM3_HI .OR. NOM1_HI/=NOM4_HI) THEN
               WRITE(*,*) 'Inconsistent neighbors on IOR=1 not allowed!'
               STOP
            ENDIF

            SLLO%IJKW(9,IW_LO)=NOM1_HI 

            IF (ILEVEL==S%NLMAX-1) THEN
               BC1_HI=M%PRESSURE_BC_INDEX(IW1_HI)
               BC2_HI=M%PRESSURE_BC_INDEX(IW2_HI)
               BC3_HI=M%PRESSURE_BC_INDEX(IW3_HI)
               BC4_HI=M%PRESSURE_BC_INDEX(IW4_HI)
            ELSE
               BC1_HI=SLHI%PRESSURE_BC_INDEX(IW1_HI)
               BC2_HI=SLHI%PRESSURE_BC_INDEX(IW2_HI)
               BC3_HI=SLHI%PRESSURE_BC_INDEX(IW3_HI)
               BC4_HI=SLHI%PRESSURE_BC_INDEX(IW4_HI)
            ENDIF

            IF (BC1_HI==INTERNAL.OR.BC2_HI==INTERNAL.OR.BC3_HI==INTERNAL.OR.BC4_HI==INTERNAL) THEN
               SLLO%PRESSURE_BC_INDEX(IW_LO)=INTERNAL
            ELSE IF (BC1_HI==DIRICHLET.OR.BC2_HI==DIRICHLET.OR.BC3_HI==DIRICHLET.OR.BC4_HI==DIRICHLET) THEN
               SLLO%PRESSURE_BC_INDEX(IW_LO)=DIRICHLET
            ELSE
               SLLO%PRESSURE_BC_INDEX(IW_LO)=NEUMANN
            ENDIF

            IF (SCARC_DEBUG.ge.2) THEN
               WRITE(SCARC_LU,*) ' 1: =========== IOR= -2'
               WRITE(SCARC_LU,*) ' 1: =========== K_LO=',K_LO,': J_LO=',J_LO
               WRITE(SCARC_LU,*) 'IW1_HI=',IW1_HI
               WRITE(SCARC_LU,*) 'IW2_HI=',IW2_HI
               WRITE(SCARC_LU,*) 'IW3_HI=',IW3_HI
               WRITE(SCARC_LU,*) 'IW4_HI=',IW4_HI
               WRITE(SCARC_LU,*) 
               WRITE(SCARC_LU,*) 'NOM1_HI=',NOM1_HI
               WRITE(SCARC_LU,*) 'NOM2_HI=',NOM2_HI
               WRITE(SCARC_LU,*) 'NOM3_HI=',NOM3_HI
               WRITE(SCARC_LU,*) 'NOM4_HI=',NOM4_HI
               WRITE(SCARC_LU,*) 
               WRITE(SCARC_LU,*) 'BC1_HI=',BC1_HI
               WRITE(SCARC_LU,*) 'BC2_HI=',BC2_HI
               WRITE(SCARC_LU,*) 'BC3_HI=',BC3_HI
               WRITE(SCARC_LU,*) 'BC4_HI=',BC4_HI
               WRITE(SCARC_LU,*) 
               WRITE(SCARC_LU,*) 'SLLO%PRESSURE_BC_INDEX(',IW_LO,')=',SLLO%PRESSURE_BC_INDEX(IW_LO)
            ENDIF

            IF (NOM1_HI > 0) THEN                        ! internal boundary ?

               JBAR_NOM=S%MJBAR(NOM1_HI)/IREFINE
               
               IMIN1_HI=SLHI%IJKW(10,IW1_HI)
               IMIN2_HI=SLHI%IJKW(10,IW2_HI)
               IMIN3_HI=SLHI%IJKW(10,IW3_HI)
               IMIN4_HI=SLHI%IJKW(10,IW4_HI)

               IMAX1_HI=SLHI%IJKW(13,IW1_HI)
               IMAX2_HI=SLHI%IJKW(13,IW2_HI)
               IMAX3_HI=SLHI%IJKW(13,IW3_HI)
               IMAX4_HI=SLHI%IJKW(13,IW4_HI)
              
               KMIN1_HI=SLHI%IJKW(12,IW1_HI)
               KMIN2_HI=SLHI%IJKW(12,IW2_HI)
               KMIN3_HI=SLHI%IJKW(12,IW3_HI)
               KMIN4_HI=SLHI%IJKW(12,IW4_HI)

               KMAX1_HI=SLHI%IJKW(15,IW1_HI)
               KMAX2_HI=SLHI%IJKW(15,IW2_HI)
               KMAX3_HI=SLHI%IJKW(15,IW3_HI)
               KMAX4_HI=SLHI%IJKW(15,IW4_HI)
              
               IDIFF=IMIN2_HI-IMIN1_HI
               KDIFF=KMIN3_HI-KMIN1_HI
         

               ! same JBAR- and KBAR-resolution at neighbor
               IF (IDIFF==1.AND.KDIFF==1) THEN
                  SLLO%IJKW(10,IW_LO)=IMAX2_HI/2
                  SLLO%IJKW(11,IW_LO)=1
                  SLLO%IJKW(12,IW_LO)=KMAX3_HI/2
                  SLLO%IJKW(13,IW_LO)=IMAX2_HI/2
                  SLLO%IJKW(14,IW_LO)=1
                  SLLO%IJKW(15,IW_LO)=KMAX3_HI/2
      
               ! finer JBAR- and KBAR-resolution at neighbor
               ELSE IF (IDIFF==2.AND.KDIFF==2) THEN
 
                  SLLO%IJKW(10,IW_LO)=IMAX1_HI/2
                  SLLO%IJKW(11,IW_LO)=1
                  SLLO%IJKW(12,IW_LO)=KMAX1_HI/2
                  SLLO%IJKW(13,IW_LO)=IMAX2_HI/2
                  SLLO%IJKW(14,IW_LO)=1
                  SLLO%IJKW(15,IW_LO)=KMAX3_HI/2
      
               ! coarser JBAR- and KBAR-resolution at neighbor
               ELSE IF (IDIFF==0.AND.KDIFF==0) THEN
 
                  SLLO%IJKW(10,IW_LO)=(IMAX1_HI+1)/2
                  SLLO%IJKW(11,IW_LO)=1
                  SLLO%IJKW(12,IW_LO)=(KMAX1_HI+1)/2
                  SLLO%IJKW(13,IW_LO)=(IMAX1_HI+1)/2
                  SLLO%IJKW(14,IW_LO)=1
                  SLLO%IJKW(15,IW_LO)=(KMAX1_HI+1)/2
      
               ELSE                                         
                  WRITE(*,*) 'Error4 in INIT_NBR_LEVEL3D at IOR=1 for ILEVEL=',ILEVEL
                  WRITE(SCARC_LU,*) 'Error4 in INIT_NBR_LEVEL3D at IOR=1 for ILEVEL=',ILEVEL
                  !STOP
               ENDIF
         
               IF (SCARC_DEBUG.ge.2) THEN
                  WRITE(SCARC_LU,*) ' 1: =========== IOR= -2'
                  WRITE(SCARC_LU,*) ' 1: =========== K_LO=',K_LO,': J_LO=',J_LO
                  WRITE(SCARC_LU,*) 'IW1_HI=',IW1_HI
                  WRITE(SCARC_LU,*) 'IW2_HI=',IW2_HI
                  WRITE(SCARC_LU,*) 'IW3_HI=',IW3_HI
                  WRITE(SCARC_LU,*) 'IW4_HI=',IW4_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'NOM1_HI=',NOM1_HI
                  WRITE(SCARC_LU,*) 'NOM2_HI=',NOM2_HI
                  WRITE(SCARC_LU,*) 'NOM3_HI=',NOM3_HI
                  WRITE(SCARC_LU,*) 'NOM4_HI=',NOM4_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'JBAR_NOM=',JBAR_NOM
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'IDIFF=',IDIFF
                  WRITE(SCARC_LU,*) 'KDIFF=',KDIFF
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'IMIN1_HI=',IMIN1_HI
                  WRITE(SCARC_LU,*) 'IMAX1_HI=',IMAX1_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'IMIN2_HI=',IMIN2_HI
                  WRITE(SCARC_LU,*) 'IMAX2_HI=',IMAX2_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'IMIN3_HI=',IMIN3_HI
                  WRITE(SCARC_LU,*) 'IMAX3_HI=',IMAX3_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'IMIN4_HI=',IMIN4_HI
                  WRITE(SCARC_LU,*) 'IMAX4_HI=',IMAX4_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'KMIN1_HI=',KMIN1_HI
                  WRITE(SCARC_LU,*) 'KMAX1_HI=',KMAX1_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'KMIN2_HI=',KMIN2_HI
                  WRITE(SCARC_LU,*) 'KMAX2_HI=',KMAX2_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'KMIN3_HI=',KMIN3_HI
                  WRITE(SCARC_LU,*) 'KMAX3_HI=',KMAX3_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'KMIN4_HI=',KMIN4_HI
                  WRITE(SCARC_LU,*) 'KMAX4_HI=',KMAX4_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU, '(i5,a,15i5)') IW_LO,':',(SLLO%IJKW(I, IW_LO), I=1,15)
               ENDIF
            ENDIF
   
            IW_LO=IW_LO+1
      
         ENDDO
      ENDDO

      !---------------------------------------------------------------------------
      ! cells along IOR=3
      !---------------------------------------------------------------------------
      IOFFSET=2*SLHI%JBAR*SLHI%KBAR+2*SLHI%IBAR*SLHI%KBAR

      DO J_LO=1,SLLO%JBAR
        DO I_LO=1,SLLO%IBAR
   
            SLLO%IJKW(1,IW_LO)=I_LO           
            SLLO%IJKW(2,IW_LO)=J_LO           ! ghost cell
            SLLO%IJKW(3,IW_LO)=0 
            SLLO%IJKW(4,IW_LO)=3              ! IOR = 3
            SLLO%IJKW(5,IW_LO)=0              ! currently not used
            SLLO%IJKW(6,IW_LO)=I_LO           
            SLLO%IJKW(7,IW_LO)=J_LO           ! boundary cell
            SLLO%IJKW(8,IW_LO)=1 
      
            IW1_HI=IOFFSET+(2*J_LO-2)*SLHI%IBAR+2*I_LO-1
            IW2_HI=IOFFSET+(2*J_LO-2)*SLHI%IBAR+2*I_LO
            IW3_HI=IOFFSET+(2*J_LO-1)*SLHI%IBAR+2*I_LO-1
            IW4_HI=IOFFSET+(2*J_LO-1)*SLHI%IBAR+2*I_LO

            NOM1_HI=SLHI%IJKW(9,IW1_HI)
            NOM2_HI=SLHI%IJKW(9,IW2_HI)
            NOM3_HI=SLHI%IJKW(9,IW3_HI)
            NOM4_HI=SLHI%IJKW(9,IW4_HI)

            IF (NOM1_HI/=NOM2_HI .OR. NOM1_HI/=NOM3_HI .OR. NOM1_HI/=NOM4_HI) THEN
               WRITE(*,*) 'Inconsistent neighbors on IOR=1 not allowed!'
               STOP
            ENDIF

            SLLO%IJKW(9,IW_LO)=NOM1_HI 

            IF (ILEVEL==S%NLMAX-1) THEN
               BC1_HI=M%PRESSURE_BC_INDEX(IW1_HI)
               BC2_HI=M%PRESSURE_BC_INDEX(IW2_HI)
               BC3_HI=M%PRESSURE_BC_INDEX(IW3_HI)
               BC4_HI=M%PRESSURE_BC_INDEX(IW4_HI)
            ELSE
               BC1_HI=SLHI%PRESSURE_BC_INDEX(IW1_HI)
               BC2_HI=SLHI%PRESSURE_BC_INDEX(IW2_HI)
               BC3_HI=SLHI%PRESSURE_BC_INDEX(IW3_HI)
               BC4_HI=SLHI%PRESSURE_BC_INDEX(IW4_HI)
            ENDIF

            IF (BC1_HI==INTERNAL.OR.BC2_HI==INTERNAL.OR.BC3_HI==INTERNAL.OR.BC4_HI==INTERNAL) THEN
               SLLO%PRESSURE_BC_INDEX(IW_LO)=INTERNAL
            ELSE IF (BC1_HI==DIRICHLET.OR.BC2_HI==DIRICHLET.OR.BC3_HI==DIRICHLET.OR.BC4_HI==DIRICHLET) THEN
               SLLO%PRESSURE_BC_INDEX(IW_LO)=DIRICHLET
            ELSE
               SLLO%PRESSURE_BC_INDEX(IW_LO)=NEUMANN
            ENDIF

            IF (SCARC_DEBUG.ge.2) THEN
               WRITE(SCARC_LU,*) ' 1: =========== IOR= 3'
               WRITE(SCARC_LU,*) ' 1: =========== K_LO=',K_LO,': J_LO=',J_LO
               WRITE(SCARC_LU,*) 'IW1_HI=',IW1_HI
               WRITE(SCARC_LU,*) 'IW2_HI=',IW2_HI
               WRITE(SCARC_LU,*) 'IW3_HI=',IW3_HI
               WRITE(SCARC_LU,*) 'IW4_HI=',IW4_HI
               WRITE(SCARC_LU,*) 
               WRITE(SCARC_LU,*) 'NOM1_HI=',NOM1_HI
               WRITE(SCARC_LU,*) 'NOM2_HI=',NOM2_HI
               WRITE(SCARC_LU,*) 'NOM3_HI=',NOM3_HI
               WRITE(SCARC_LU,*) 'NOM4_HI=',NOM4_HI
               WRITE(SCARC_LU,*) 
               WRITE(SCARC_LU,*) 'BC1_HI=',BC1_HI
               WRITE(SCARC_LU,*) 'BC2_HI=',BC2_HI
               WRITE(SCARC_LU,*) 'BC3_HI=',BC3_HI
               WRITE(SCARC_LU,*) 'BC4_HI=',BC4_HI
               WRITE(SCARC_LU,*) 
               WRITE(SCARC_LU,*) 'SLLO%PRESSURE_BC_INDEX(',IW_LO,')=',SLLO%PRESSURE_BC_INDEX(IW_LO)
            ENDIF

            IF (NOM1_HI > 0) THEN                        ! internal boundary ?

               KBAR_NOM=S%MKBAR(NOM1_HI)/IREFINE
               
               IMIN1_HI=SLHI%IJKW(10,IW1_HI)
               IMIN2_HI=SLHI%IJKW(10,IW2_HI)
               IMIN3_HI=SLHI%IJKW(10,IW3_HI)
               IMIN4_HI=SLHI%IJKW(10,IW4_HI)

               IMAX1_HI=SLHI%IJKW(13,IW1_HI)
               IMAX2_HI=SLHI%IJKW(13,IW2_HI)
               IMAX3_HI=SLHI%IJKW(13,IW3_HI)
               IMAX4_HI=SLHI%IJKW(13,IW4_HI)
              
               JMIN1_HI=SLHI%IJKW(11,IW1_HI)
               JMIN2_HI=SLHI%IJKW(11,IW2_HI)
               JMIN3_HI=SLHI%IJKW(11,IW3_HI)
               JMIN4_HI=SLHI%IJKW(11,IW4_HI)

               JMAX1_HI=SLHI%IJKW(14,IW1_HI)
               JMAX2_HI=SLHI%IJKW(14,IW2_HI)
               JMAX3_HI=SLHI%IJKW(14,IW3_HI)
               JMAX4_HI=SLHI%IJKW(14,IW4_HI)
              
               IDIFF=IMIN2_HI-IMIN1_HI
               JDIFF=JMIN3_HI-JMIN1_HI

               ! same IBAR- and JBAR-resolution at neighbor
               IF (IDIFF==1.AND.JDIFF==1) THEN
                  SLLO%IJKW(10,IW_LO)=IMAX2_HI/2
                  SLLO%IJKW(11,IW_LO)=JMAX3_HI/2
                  SLLO%IJKW(12,IW_LO)=KBAR_NOM
                  SLLO%IJKW(13,IW_LO)=IMAX2_HI/2
                  SLLO%IJKW(14,IW_LO)=JMAX3_HI/2
                  SLLO%IJKW(15,IW_LO)=KBAR_NOM
      
               ! finer JBAR- and JBAR-resolution at neighbor
               ELSE IF (IDIFF==2.AND.JDIFF==2) THEN
 
                  SLLO%IJKW(10,IW_LO)=IMAX1_HI/2
                  SLLO%IJKW(11,IW_LO)=JMAX1_HI/2
                  SLLO%IJKW(12,IW_LO)=KBAR_NOM
                  SLLO%IJKW(13,IW_LO)=IMAX2_HI/2
                  SLLO%IJKW(14,IW_LO)=JMAX3_HI/2
                  SLLO%IJKW(15,IW_LO)=KBAR_NOM
      
               ! coarser IBAR- and JBAR-resolution at neighbor
               ELSE IF (IDIFF==0.AND.JDIFF==0) THEN
 
                  SLLO%IJKW(10,IW_LO)=(IMAX2_HI+1)/2
                  SLLO%IJKW(11,IW_LO)=(JMAX3_HI+1)/2
                  SLLO%IJKW(12,IW_LO)=KBAR_NOM
                  SLLO%IJKW(13,IW_LO)=(IMAX2_HI+1)/2
                  SLLO%IJKW(14,IW_LO)=(JMAX3_HI+1)/2
                  SLLO%IJKW(15,IW_LO)=KBAR_NOM
      
               ELSE                                         
                  WRITE(*,*) 'Error5 in INIT_NBR_LEVEL3D at IOR=1 for ILEVEL=',ILEVEL
                  WRITE(SCARC_LU,*) 'Error5 in INIT_NBR_LEVEL3D at IOR=1 for ILEVEL=',ILEVEL
                  !STOP
               ENDIF
         
               IF (SCARC_DEBUG.ge.2) THEN
                  WRITE(SCARC_LU,*) ' 1: =========== IOR= -3'
                  WRITE(SCARC_LU,*) ' 1: =========== K_LO=',K_LO,': J_LO=',J_LO
                  WRITE(SCARC_LU,*) 'IW1_HI=',IW1_HI
                  WRITE(SCARC_LU,*) 'IW2_HI=',IW2_HI
                  WRITE(SCARC_LU,*) 'IW3_HI=',IW3_HI
                  WRITE(SCARC_LU,*) 'IW4_HI=',IW4_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'NOM1_HI=',NOM1_HI
                  WRITE(SCARC_LU,*) 'NOM2_HI=',NOM2_HI
                  WRITE(SCARC_LU,*) 'NOM3_HI=',NOM3_HI
                  WRITE(SCARC_LU,*) 'NOM4_HI=',NOM4_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'JBAR_NOM=',JBAR_NOM
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'IDIFF=',IDIFF
                  WRITE(SCARC_LU,*) 'KDIFF=',KDIFF
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'IMIN1_HI=',IMIN1_HI
                  WRITE(SCARC_LU,*) 'IMAX1_HI=',IMAX1_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'IMIN2_HI=',IMIN2_HI
                  WRITE(SCARC_LU,*) 'IMAX2_HI=',IMAX2_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'IMIN3_HI=',IMIN3_HI
                  WRITE(SCARC_LU,*) 'IMAX3_HI=',IMAX3_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'IMIN4_HI=',IMIN4_HI
                  WRITE(SCARC_LU,*) 'IMAX4_HI=',IMAX4_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'KMIN1_HI=',KMIN1_HI
                  WRITE(SCARC_LU,*) 'KMAX1_HI=',KMAX1_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'KMIN2_HI=',KMIN2_HI
                  WRITE(SCARC_LU,*) 'KMAX2_HI=',KMAX2_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'KMIN3_HI=',KMIN3_HI
                  WRITE(SCARC_LU,*) 'KMAX3_HI=',KMAX3_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'KMIN4_HI=',KMIN4_HI
                  WRITE(SCARC_LU,*) 'KMAX4_HI=',KMAX4_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU, '(i5,a,15i5)') IW_LO,':',(SLLO%IJKW(I, IW_LO), I=1,15)
               ENDIF
            ENDIF
   
            IW_LO=IW_LO+1
      
         ENDDO
      ENDDO

      !---------------------------------------------------------------------------
      ! cells along IOR=-3
      !---------------------------------------------------------------------------
      IOFFSET=2*SLHI%JBAR*SLHI%KBAR+2*SLHI%IBAR*SLHI%KBAR+SLHI%IBAR*SLHI%JBAR

      DO J_LO=1,SLLO%JBAR
        DO I_LO=1,SLLO%IBAR
   
            SLLO%IJKW(1,IW_LO)=I_LO           
            SLLO%IJKW(2,IW_LO)=J_LO           ! ghost cell
            SLLO%IJKW(3,IW_LO)=SLLO%KBAR+1 
            SLLO%IJKW(4,IW_LO)=-3             ! IOR = -3
            SLLO%IJKW(5,IW_LO)=0              ! currently not used
            SLLO%IJKW(6,IW_LO)=I_LO           
            SLLO%IJKW(7,IW_LO)=J_LO           ! boundary cell
            SLLO%IJKW(8,IW_LO)=SLLO%KBAR 
      
            IW1_HI=IOFFSET+(2*J_LO-2)*SLHI%IBAR+2*I_LO-1
            IW2_HI=IOFFSET+(2*J_LO-2)*SLHI%IBAR+2*I_LO
            IW3_HI=IOFFSET+(2*J_LO-1)*SLHI%IBAR+2*I_LO-1
            IW4_HI=IOFFSET+(2*J_LO-1)*SLHI%IBAR+2*I_LO

            NOM1_HI=SLHI%IJKW(9,IW1_HI)
            NOM2_HI=SLHI%IJKW(9,IW2_HI)
            NOM3_HI=SLHI%IJKW(9,IW3_HI)
            NOM4_HI=SLHI%IJKW(9,IW4_HI)
      
            IF (NOM1_HI/=NOM2_HI .OR. NOM1_HI/=NOM3_HI .OR. NOM1_HI/=NOM4_HI) THEN
               WRITE(*,*) 'Inconsistent neighbors on IOR=1 not allowed!'
               STOP
            ENDIF

            SLLO%IJKW(9,IW_LO)=NOM1_HI 

            IF (ILEVEL==S%NLMAX-1) THEN
               BC1_HI=M%PRESSURE_BC_INDEX(IW1_HI)
               BC2_HI=M%PRESSURE_BC_INDEX(IW2_HI)
               BC3_HI=M%PRESSURE_BC_INDEX(IW3_HI)
               BC4_HI=M%PRESSURE_BC_INDEX(IW4_HI)
            ELSE
               BC1_HI=SLHI%PRESSURE_BC_INDEX(IW1_HI)
               BC2_HI=SLHI%PRESSURE_BC_INDEX(IW2_HI)
               BC3_HI=SLHI%PRESSURE_BC_INDEX(IW3_HI)
               BC4_HI=SLHI%PRESSURE_BC_INDEX(IW4_HI)
            ENDIF

            IF (BC1_HI==INTERNAL.OR.BC2_HI==INTERNAL.OR.BC3_HI==INTERNAL.OR.BC4_HI==INTERNAL) THEN
               SLLO%PRESSURE_BC_INDEX(IW_LO)=INTERNAL
            ELSE IF (BC1_HI==DIRICHLET.OR.BC2_HI==DIRICHLET.OR.BC3_HI==DIRICHLET.OR.BC4_HI==DIRICHLET) THEN
               SLLO%PRESSURE_BC_INDEX(IW_LO)=DIRICHLET
            ELSE
               SLLO%PRESSURE_BC_INDEX(IW_LO)=NEUMANN
            ENDIF

            IF (SCARC_DEBUG.ge.2) THEN
               WRITE(SCARC_LU,*) ' 1: =========== IOR= 3'
               WRITE(SCARC_LU,*) ' 1: =========== K_LO=',K_LO,': J_LO=',J_LO
               WRITE(SCARC_LU,*) 'IW1_HI=',IW1_HI
               WRITE(SCARC_LU,*) 'IW2_HI=',IW2_HI
               WRITE(SCARC_LU,*) 'IW3_HI=',IW3_HI
               WRITE(SCARC_LU,*) 'IW4_HI=',IW4_HI
               WRITE(SCARC_LU,*) 
               WRITE(SCARC_LU,*) 'NOM1_HI=',NOM1_HI
               WRITE(SCARC_LU,*) 'NOM2_HI=',NOM2_HI
               WRITE(SCARC_LU,*) 'NOM3_HI=',NOM3_HI
               WRITE(SCARC_LU,*) 'NOM4_HI=',NOM4_HI
               WRITE(SCARC_LU,*) 
               WRITE(SCARC_LU,*) 'BC1_HI=',BC1_HI
               WRITE(SCARC_LU,*) 'BC2_HI=',BC2_HI
               WRITE(SCARC_LU,*) 'BC3_HI=',BC3_HI
               WRITE(SCARC_LU,*) 'BC4_HI=',BC4_HI
               WRITE(SCARC_LU,*) 
               WRITE(SCARC_LU,*) 'SLLO%PRESSURE_BC_INDEX(',IW_LO,')=',SLLO%PRESSURE_BC_INDEX(IW_LO)
            ENDIF

            IF (NOM1_HI > 0) THEN                        ! internal boundary ?

               KBAR_NOM=S%MKBAR(NOM1_HI)/IREFINE
               
               IMIN1_HI=SLHI%IJKW(10,IW1_HI)
               IMIN2_HI=SLHI%IJKW(10,IW2_HI)
               IMIN3_HI=SLHI%IJKW(10,IW3_HI)
               IMIN4_HI=SLHI%IJKW(10,IW4_HI)

               IMAX1_HI=SLHI%IJKW(13,IW1_HI)
               IMAX2_HI=SLHI%IJKW(13,IW2_HI)
               IMAX3_HI=SLHI%IJKW(13,IW3_HI)
               IMAX4_HI=SLHI%IJKW(13,IW4_HI)
              
               JMIN1_HI=SLHI%IJKW(11,IW1_HI)
               JMIN2_HI=SLHI%IJKW(11,IW2_HI)
               JMIN3_HI=SLHI%IJKW(11,IW3_HI)
               JMIN4_HI=SLHI%IJKW(11,IW4_HI)

               JMAX1_HI=SLHI%IJKW(14,IW1_HI)
               JMAX2_HI=SLHI%IJKW(14,IW2_HI)
               JMAX3_HI=SLHI%IJKW(14,IW3_HI)
               JMAX4_HI=SLHI%IJKW(14,IW4_HI)
              
               IDIFF=IMIN2_HI-IMIN1_HI
               JDIFF=JMIN3_HI-JMIN1_HI

            IF (SCARC_DEBUG.ge.2) THEN
               WRITE(SCARC_LU,*) 
               WRITE(SCARC_LU,*) 'IMIN1_HI=',IMIN1_HI
               WRITE(SCARC_LU,*) 'IMIN2_HI=',IMIN2_HI
               WRITE(SCARC_LU,*) 'IMIN3_HI=',IMIN3_HI
               WRITE(SCARC_LU,*) 'IMIN4_HI=',IMIN4_HI
               WRITE(SCARC_LU,*) 
               WRITE(SCARC_LU,*) 'JMIN1_HI=',JMIN1_HI
               WRITE(SCARC_LU,*) 'JMIN2_HI=',JMIN2_HI
               WRITE(SCARC_LU,*) 'JMIN3_HI=',JMIN3_HI
               WRITE(SCARC_LU,*) 'JMIN4_HI=',JMIN4_HI
               WRITE(SCARC_LU,*) 
               WRITE(SCARC_LU,*) 'IDIFF=',IDIFF
               WRITE(SCARC_LU,*) 'JDIFF=',JDIFF
            ENDIF

               ! same IBAR- and JBAR-resolution at neighbor
               IF (IDIFF==1.AND.JDIFF==1) THEN
                  SLLO%IJKW(10,IW_LO)=IMAX2_HI/2
                  SLLO%IJKW(11,IW_LO)=JMAX3_HI/2
                  SLLO%IJKW(12,IW_LO)=1
                  SLLO%IJKW(13,IW_LO)=IMAX2_HI/2
                  SLLO%IJKW(14,IW_LO)=JMAX3_HI/2
                  SLLO%IJKW(15,IW_LO)=1
      
               ! finer JBAR- and JBAR-resolution at neighbor
               ELSE IF (IDIFF==2.AND.JDIFF==2) THEN
 
                  SLLO%IJKW(10,IW_LO)=IMAX1_HI/2
                  SLLO%IJKW(11,IW_LO)=JMAX1_HI/2
                  SLLO%IJKW(12,IW_LO)=1
                  SLLO%IJKW(13,IW_LO)=IMAX2_HI/2
                  SLLO%IJKW(14,IW_LO)=JMAX3_HI/2
                  SLLO%IJKW(15,IW_LO)=1
      
               ! coarser IBAR- and JBAR-resolution at neighbor
               ELSE IF (IDIFF==0.AND.JDIFF==0) THEN
 
                  SLLO%IJKW(10,IW_LO)=(IMAX2_HI+1)/2
                  SLLO%IJKW(11,IW_LO)=(JMAX3_HI+1)/2
                  SLLO%IJKW(12,IW_LO)=1
                  SLLO%IJKW(13,IW_LO)=(IMAX2_HI+1)/2
                  SLLO%IJKW(14,IW_LO)=(JMAX3_HI+1)/2
                  SLLO%IJKW(15,IW_LO)=1
      
               ELSE                                         
                  WRITE(*,*) 'Error6 in INIT_NBR_LEVEL3D at IOR=1 for ILEVEL=',ILEVEL
                  WRITE(SCARC_LU,*) 'Error6 in INIT_NBR_LEVEL3D at IOR=1 for ILEVEL=',ILEVEL
                  !STOP
               ENDIF
         
               IF (SCARC_DEBUG.ge.2) THEN
                  WRITE(SCARC_LU,*) ' 1: =========== IOR= -3'
                  WRITE(SCARC_LU,*) ' 1: =========== K_LO=',K_LO,': J_LO=',J_LO
                  WRITE(SCARC_LU,*) 'IW1_HI=',IW1_HI
                  WRITE(SCARC_LU,*) 'IW2_HI=',IW2_HI
                  WRITE(SCARC_LU,*) 'IW3_HI=',IW3_HI
                  WRITE(SCARC_LU,*) 'IW4_HI=',IW4_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'NOM1_HI=',NOM1_HI
                  WRITE(SCARC_LU,*) 'NOM2_HI=',NOM2_HI
                  WRITE(SCARC_LU,*) 'NOM3_HI=',NOM3_HI
                  WRITE(SCARC_LU,*) 'NOM4_HI=',NOM4_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'JBAR_NOM=',JBAR_NOM
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'IDIFF=',IDIFF
                  WRITE(SCARC_LU,*) 'KDIFF=',KDIFF
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'IMIN1_HI=',IMIN1_HI
                  WRITE(SCARC_LU,*) 'IMAX1_HI=',IMAX1_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'IMIN2_HI=',IMIN2_HI
                  WRITE(SCARC_LU,*) 'IMAX2_HI=',IMAX2_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'IMIN3_HI=',IMIN3_HI
                  WRITE(SCARC_LU,*) 'IMAX3_HI=',IMAX3_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'IMIN4_HI=',IMIN4_HI
                  WRITE(SCARC_LU,*) 'IMAX4_HI=',IMAX4_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'KMIN1_HI=',KMIN1_HI
                  WRITE(SCARC_LU,*) 'KMAX1_HI=',KMAX1_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'KMIN2_HI=',KMIN2_HI
                  WRITE(SCARC_LU,*) 'KMAX2_HI=',KMAX2_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'KMIN3_HI=',KMIN3_HI
                  WRITE(SCARC_LU,*) 'KMAX3_HI=',KMAX3_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'KMIN4_HI=',KMIN4_HI
                  WRITE(SCARC_LU,*) 'KMAX4_HI=',KMAX4_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU, '(i5,a,15i5)') IW_LO,':',(SLLO%IJKW(I, IW_LO), I=1,15)
               ENDIF
            ENDIF
   
            IW_LO=IW_LO+1
      
         ENDDO
      ENDDO

      IF (SCARC_DEBUG >= 6) THEN
         WRITE(SCARC_LU,*) '=================================================='
         WRITE(SCARC_LU,*) '================== LEVEL =',ILEVEL
         WRITE(SCARC_LU,*) 'N_EXTERNAL_WALL_CELLS=',SLLO%N_EXTERNAL_WALL_CELLS
         WRITE(SCARC_LU,*) 'IJKW:'
         DO INEWC=1,SLLO%N_EXTERNAL_WALL_CELLS
            WRITE(SCARC_LU, '(i5,a,15i5)') INEWC,':',(SLLO%IJKW(I, INEWC), I=1,15)
         ENDDO
         WRITE(SCARC_LU,*) 'PRESSURE_BC_INDEX:'
         WRITE(SCARC_LU, '(4i4)') (SLLO%PRESSURE_BC_INDEX(INEWC), INEWC=1,SLLO%N_EXTERNAL_WALL_CELLS)
      ENDIF

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! compute analogues to NIC, I_MIN, I_MAX, K_MIN, K_MAX on level 'ILEVEL'
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ALLOCATE (SLLO%NIC(NMESHES, NMESHES))
      SLLO%NIC=0

      ALLOCATE (SLLO%I_MIN(NMESHES, NMESHES))
      ALLOCATE (SLLO%I_MAX(NMESHES, NMESHES))
      ALLOCATE (SLLO%J_MIN(NMESHES, NMESHES))
      ALLOCATE (SLLO%J_MAX(NMESHES, NMESHES))
      ALLOCATE (SLLO%K_MIN(NMESHES, NMESHES))
      ALLOCATE (SLLO%K_MAX(NMESHES, NMESHES))
      SLLO%I_MIN = 0
      SLLO%I_MAX = 0
      SLLO%J_MIN = 0
      SLLO%J_MAX = 0
      SLLO%K_MIN = 0
      SLLO%K_MAX = 0


      OTHER_MESH_LOOP3D: DO NOM=1,NMESHES
 
         IF (NIC(NM,NOM)==0.AND.NIC(NOM,NM)==0) CYCLE OTHER_MESH_LOOP3D
 
         !!! ACHTUNG: funktioniert nur für 2er-Potenz-Gitterweiten !!!!!
         IMIN=0
         IMAX=S%MIBAR(NOM)/IREFINE+1
         JMIN=0
         JMAX=S%MJBAR(NOM)/IREFINE+1
         KMIN=0
         KMAX=S%MKBAR(NOM)/IREFINE+1

         IF (SCARC_DEBUG.ge.6) THEN
            WRITE(SCARC_LU,*) '============== NOM=',NOM
            WRITE(SCARC_LU,*) 'NIC:'
            DO JM=1,NMESHES
               WRITE(SCARC_LU, '(10i5)') (NIC(IM,JM),IM=1,NMESHES)
            ENDDO
            WRITE(SCARC_LU,*) 'MIBAR(NOM)=',S%MIBAR(NOM)
            WRITE(SCARC_LU,*) 'MJBAR(NOM)=',S%MJBAR(NOM)
            WRITE(SCARC_LU,*) 'MKBAR(NOM)=',S%MKBAR(NOM)
            WRITE(SCARC_LU,*) 'IREFINE=',IREFINE
            WRITE(SCARC_LU,*) 'IMAX=',IMAX
            WRITE(SCARC_LU,*) 'JMAX=',JMAX
            WRITE(SCARC_LU,*) 'KMAX=',KMAX
         ENDIF

         SLLO%NIC(NOM,NM) = 0

         SEARCH_LOOP: DO INEWC=1,SLLO%N_EXTERNAL_WALL_CELLS

            ! neighborship structure already known from finest level
            IF (SLLO%IJKW(9,INEWC)/=NOM) CYCLE SEARCH_LOOP

            !WRITE(SCARC_LU,*) 'IJKW(9,',INEWC,')=',SLLO%IJKW(9,INEWC)

            SLLO%NIC(NOM,NM) = SLLO%NIC(NOM,NM) + 1
            IOR = SLLO%IJKW(4,INEWC)

            !WRITE(SCARC_LU,*) 'NIC(',NOM,',',NM,')=',SLLO%NIC(NOM,NM)

            SELECT CASE(IOR)
               CASE( 1)
                  IMIN=MAX(IMIN,SLLO%IJKW(10,INEWC)-1)
               CASE(-1) 
                  IMAX=MIN(IMAX,SLLO%IJKW(13,INEWC))
               CASE( 2) 
                  JMIN=MAX(JMIN,SLLO%IJKW(11,INEWC)-1)
               CASE(-2) 
                  JMAX=MIN(JMAX,SLLO%IJKW(14,INEWC))
               CASE( 3) 
                  KMIN=MAX(KMIN,SLLO%IJKW(12,INEWC)-1)
               CASE(-3)
                  KMAX=MIN(KMAX,SLLO%IJKW(15,INEWC))
            END SELECT
         ENDDO SEARCH_LOOP

         SLLO%I_MIN(NOM,NM) = IMIN
         SLLO%I_MAX(NOM,NM) = IMAX
         SLLO%J_MIN(NOM,NM) = JMIN
         SLLO%J_MAX(NOM,NM) = JMAX
         SLLO%K_MIN(NOM,NM) = KMIN
         SLLO%K_MAX(NOM,NM) = KMAX

      ENDDO OTHER_MESH_LOOP3D

      IF (SCARC_DEBUG >= 6) THEN
         WRITE(SCARC_LU,*) 'NUMPROCS=',NUMPROCS
         WRITE(SCARC_LU,*) 'MYID  =',MYID  
         WRITE(SCARC_LU,*) 'COUNTS=',COUNTS
         WRITE(SCARC_LU,*) 'COUNTS3D=',COUNTS3D
         WRITE(SCARC_LU,*) 'DISPLS=',DISPLS
         WRITE(SCARC_LU,*) 'DISPLS3D=',DISPLS3D
         WRITE(SCARC_LU,*) 'DISPLS(MYID)+1=',DISPLS(MYID)+1
         WRITE(SCARC_LU,*) 'COUNTS3D(MYID)=',COUNTS3D(MYID)
         WRITE(SCARC_LU,*) 'I_MIN(1,DISPLS(MYID)+1:DISPLS(MYID)+3)=',SLLO%I_MIN(1,DISPLS(MYID)+1:DISPLS(MYID)+3)
         WRITE(SCARC_LU,*) '=======================A:'
         DO JM=1,NMESHES
            WRITE(SCARC_LU, '(10i5)') (SLLO%I_MIN(IM,JM),IM=1,NMESHES)
         ENDDO
      ENDIF

      IF (USE_MPI) THEN
         CALL MPI_ALLGATHERV(SLLO%I_MIN(1,DISPLS(MYID)+1),COUNTS3D(MYID),MPI_INTEGER,&
                             SLLO%I_MIN,COUNTS3D,DISPLS3D,MPI_INTEGER,MPI_COMM_WORLD,IERR)
         CALL MPI_ALLGATHERV(SLLO%I_MAX(1,DISPLS(MYID)+1),COUNTS3D(MYID),MPI_INTEGER,&
                             SLLO%I_MAX,COUNTS3D,DISPLS3D,MPI_INTEGER,MPI_COMM_WORLD,IERR)
         CALL MPI_ALLGATHERV(SLLO%J_MIN(1,DISPLS(MYID)+1),COUNTS3D(MYID),MPI_INTEGER,&
                             SLLO%J_MIN,COUNTS3D,DISPLS3D,MPI_INTEGER,MPI_COMM_WORLD,IERR)
         CALL MPI_ALLGATHERV(SLLO%J_MAX(1,DISPLS(MYID)+1),COUNTS3D(MYID),MPI_INTEGER,&
                             SLLO%J_MAX,COUNTS3D,DISPLS3D,MPI_INTEGER,MPI_COMM_WORLD,IERR)
         CALL MPI_ALLGATHERV(SLLO%K_MIN(1,DISPLS(MYID)+1),COUNTS3D(MYID),MPI_INTEGER,&
                             SLLO%K_MIN,COUNTS3D,DISPLS3D,MPI_INTEGER,MPI_COMM_WORLD,IERR)
         CALL MPI_ALLGATHERV(SLLO%K_MAX(1,DISPLS(MYID)+1),COUNTS3D(MYID),MPI_INTEGER,&
                             SLLO%K_MAX,COUNTS3D,DISPLS3D,MPI_INTEGER,MPI_COMM_WORLD,IERR)
         CALL MPI_ALLGATHERV(SLLO%NIC(1,DISPLS(MYID)+1),  COUNTS3D(MYID),MPI_INTEGER,&
                             SLLO%NIC,  COUNTS3D,DISPLS3D,MPI_INTEGER,MPI_COMM_WORLD,IERR)
      ELSE
         WRITE(*,*) 'Serial version not yet implemented'
         stop
      ENDIF

      SLLO%I_MIN = TRANSPOSE(SLLO%I_MIN)
      SLLO%I_MAX = TRANSPOSE(SLLO%I_MAX)
      SLLO%J_MIN = TRANSPOSE(SLLO%J_MIN)
      SLLO%J_MAX = TRANSPOSE(SLLO%J_MAX)
      SLLO%K_MIN = TRANSPOSE(SLLO%K_MIN)
      SLLO%K_MAX = TRANSPOSE(SLLO%K_MAX)
      SLLO%NIC   = TRANSPOSE(SLLO%NIC)

      IF (SCARC_DEBUG >= 6) THEN
         WRITE(SCARC_LU,*) '=================================================='
         WRITE(SCARC_LU,*) '================== LEVEL =',ILEVEL
         WRITE(SCARC_LU,*) 'N_EXTERNAL_WALL_CELLS=',SLLO%N_EXTERNAL_WALL_CELLS
         WRITE(SCARC_LU,*) 'IJKW:'
         DO INEWC=1,SLLO%N_EXTERNAL_WALL_CELLS
            WRITE(SCARC_LU, '(i5,a,15i5)') INEWC,':',(SLLO%IJKW(I, INEWC), I=1,15)
         ENDDO
         WRITE(SCARC_LU,*) 'NIC:'
         WRITE(SCARC_LU, '(3i4)') ((SLLO%NIC(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
         WRITE(SCARC_LU,*) 'I_MIN:'
         WRITE(SCARC_LU, '(3i4)') ((SLLO%I_MIN(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
         WRITE(SCARC_LU,*) 'I_MAX:'
         WRITE(SCARC_LU, '(3i4)') ((SLLO%I_MAX(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
         WRITE(SCARC_LU,*) 'J_MIN:'
         WRITE(SCARC_LU, '(3i4)') ((SLLO%J_MIN(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
         WRITE(SCARC_LU,*) 'J_MAX:'
         WRITE(SCARC_LU, '(3i4)') ((SLLO%J_MAX(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
         WRITE(SCARC_LU,*) 'K_MIN:'
         WRITE(SCARC_LU, '(3i4)') ((SLLO%K_MIN(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
         WRITE(SCARC_LU,*) 'K_MAX:'
         WRITE(SCARC_LU, '(3i4)') ((SLLO%K_MAX(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
         DO IM=1,NMESHES
            IF (IM/=NM) THEN
               WRITE(SCARC_LU,*) 'FOR MESH ',IM,' ALLOCATING: (',&
                           SLLO%I_MIN(NM,IM),':',SLLO%I_MAX(NM,IM),' , ',&
                           SLLO%J_MIN(NM,IM),':',SLLO%J_MAX(NM,IM),' , ',&
                           SLLO%K_MIN(NM,IM),':',SLLO%K_MAX(NM,IM),')'
            ENDIF
         ENDDO
      ENDIF
       
   ENDDO INIT_NBR_LEVEL3D
ENDIF
    
 
DEALLOCATE(DISPLS3D)
DEALLOCATE(DISPLS)
DEALLOCATE(COUNTS3D)
DEALLOCATE(COUNTS)

TUSED_SCARC(2,NM)=TUSED_SCARC(2,NM)+SECOND()-TNOW_NEIGHBORS3D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_NEIGHBORS3D

END SUBROUTINE SCARC_INIT_NEIGHBORS3D


 
!------------------------------------------------------------------------------------
! determine neighborship/communication structure for data exchange in 2D
!------------------------------------------------------------------------------------
SUBROUTINE SCARC_INIT_COMMUNICATION2D (NM)
 
INTEGER :: IM, NM, NOM, ILEVEL, IREFINE
INTEGER :: IMIN, IMAX, JMIN, JMAX, KMIN, KMAX, ILMAX
INTEGER :: TAG_FACE
INTEGER :: IERR
REAL(EB):: TNOW_COMMUNICATION2D

TNOW_COMMUNICATION2D = SECOND()
 
S => SCARC (NM)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Allocate protocol arrays for data exchange
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ALLOCATE (STAT_FACE(MPI_STATUS_SIZE, NMESHES*NMESHES*20))
 
ALLOCATE (REQ_FACE(NMESHES*NMESHES*20))
REQ_FACE = MPI_REQUEST_NULL

ALLOCATE (TAGS_FACE(NMESHES, NMESHES))
TAGS_FACE = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Create unique tags arrays for face, edge and vertex exchanges
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
TAG_FACE = 0
 
DO IM = 1, NMESHES
   DO NOM = IM, NMESHES
      TAG_FACE = TAG_FACE + 1
      TAGS_FACE (IM, NOM) = TAG_FACE
      TAGS_FACE (NOM, IM) = TAG_FACE
   ENDDO
ENDDO 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Initialize level structures on neighboring meshes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IERR=0
OSLEVEL_LOOP2D: DO NOM=1,NMESHES

   ILMAX =  S%NLMAX
   SLMAX => S%SLEVEL(ILMAX)

   IF (SLMAX%NIC(NM,NOM)/=0.AND.SLMAX%NIC(NOM,NM)/=0) THEN

      OS => S%OSCARC(NOM)
      ALLOCATE (OS%SLEVEL(S%NLMIN:S%NLMAX), STAT=IERR)
      CALL CHKMEMERR ('SCARC_INIT_COMMUNICATION2D', 'OS%SLEVEL', IERR)

      OSLMAX => OS%SLEVEL(ILMAX)
      OSLMAX%N_EXTERNAL_WALL_CELLS = 2*S%MIBAR(NOM)              + &
                    2*S%MIBAR(NOM)*S%MKBAR(NOM) + &
                    2*S%MKBAR(NOM)
      !ALLOCATE (OSLMAX%IJKW(15,OSLMAX%N_EXTERNAL_WALL_CELLS), STAT=IERR)
      !CALL CHKMEMERR ('SCARC_INIT_COMMUNICATION2D', 'OSLMAX%IJKW', IERR)
      !OSLMAX%IJKW = 0

      IF (SCARC_DEBUG>=2) THEN
         WRITE(SCARC_LU,*) 'OSLMAX zeigt auf OS%SLEVEL(',ILMAX,')'
         WRITE(SCARC_LU,*) 'OSLMAX%N_EXTERNAL_WALL_CELLS=',OSLMAX%N_EXTERNAL_WALL_CELLS
         !WRITE(SCARC_LU,*) 'SIZE(OSLMAX%IJKW)=',SIZE(OSLMAX%IJKW)
         WRITE(SCARC_LU,*) 'S%NLDIFF=',S%NLDIFF
      ENDIF

      IF (S%NLDIFF/=0) THEN
         IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,*) 'S%NLDIFF=',S%NLDIFF, S%NLMAX, S%NLMIN
         IREFINE=1
         DO ILEVEL=S%NLMAX-1,S%NLMIN,-1

            IREFINE=IREFINE*2

         IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,*) 'ILEVEL=',ILEVEL
         IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,*) 'IREFINE=',IREFINE

            OS%SLEVEL(ILEVEL)%N_EXTERNAL_WALL_CELLS=OS%SLEVEL(ILEVEL+1)%N_EXTERNAL_WALL_CELLS/IREFINE
            ALLOCATE (OS%SLEVEL(ILEVEL)%IJKW(15,OSLMAX%N_EXTERNAL_WALL_CELLS), STAT=IERR)
            CALL CHKMEMERR ('SCARC_INIT_COMMUNICATION2D', 'OS%SLEVEL%IJKW', IERR)
            OS%SLEVEL(ILEVEL)%IJKW = 0

            IF (SCARC_DEBUG>=2) THEN
               WRITE(SCARC_LU,*) 'OSL zeigt auf OS%SLEVEL(',ILEVEL,')'
               WRITE(SCARC_LU,*) 'OS%SLEVEL(',ILEVEL,')%N_EXTERNAL_WALL_CELLS=',OS%SLEVEL(ILEVEL)%N_EXTERNAL_WALL_CELLS
               WRITE(SCARC_LU,*) 'SIZE(OS%SLEVEL(ILEVEL)%IJKW)=',SIZE(OS%SLEVEL(ILEVEL)%IJKW)
            ENDIF

         ENDDO
      ENDIF

   ENDIF

ENDDO OSLEVEL_LOOP2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Initialize communication structures on every level
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
COM_LEVEL_LOOP2D: DO ILEVEL=S%NLMAX,S%NLMIN,-1

   SL => S%SLEVEL(ILEVEL)

   !!! ---------------------------------------------------------------------------
   !!! Communication structures for faces
   !!! ---------------------------------------------------------------------------
   FACE_NBR_LOOP2D: DO NOM = 1, NMESHES

      ! is there a face neighbor ?
      IF (NOM == NM) CYCLE FACE_NBR_LOOP2D
      IF (SL%NIC(NM, NOM) == 0 .AND. SL%NIC(NOM, NM) == 0) CYCLE FACE_NBR_LOOP2D

      ! get length of communication vectors for different orientation of neighbors
      IMIN = SL%I_MIN(NM,NOM)
      IMAX = SL%I_MAX(NM,NOM)
      JMIN = SL%J_MIN(NM,NOM)
      JMAX = SL%J_MAX(NM,NOM)
      KMIN = SL%K_MIN(NM,NOM)
      KMAX = SL%K_MAX(NM,NOM)

      ! allocate communication vectors for data exchange with face neighbor NOM on level ILEVEL
      OSL =>S%OSCARC(NOM)%SLEVEL(ILEVEL)
   
      ALLOCATE (OSL%Y_FACE(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
      OSL%Y_FACE = 0.0_EB
    
      ALLOCATE (OSL%R_FACE(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
      OSL%R_FACE = 0.0_EB
    
      ALLOCATE (OSL%G_FACE(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
      OSL%G_FACE = 0.0_EB
   
      ALLOCATE (OSL%Z_FACE(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
      OSL%Z_FACE = 0.0_EB
   
      ! print information if requested    
      IF (SCARC_DEBUG >= 6) THEN
         WRITE(SCARC_LU,*) '---NM=', NM
         WRITE(SCARC_LU,*) '---NOM=', NOM
         WRITE(SCARC_LU,*) 'ILEVEL',ILEVEL
         WRITE(SCARC_LU,*) 'SL%I_MIN=', SL%I_MIN
         WRITE(SCARC_LU,*) 'SL%I_MAX=', SL%I_MAX
         WRITE(SCARC_LU,*) 'SL%K_MIN=', SL%K_MIN
         WRITE(SCARC_LU,*) 'SL%K_MAX=', SL%K_MAX
         WRITE(SCARC_LU,*) 'IMIN=', IMIN
         WRITE(SCARC_LU,*) 'IMAX=', IMAX
         WRITE(SCARC_LU,*) 'JMIN=', JMIN
         WRITE(SCARC_LU,*) 'JMAX=', JMAX
         WRITE(SCARC_LU,*) 'KMIN=', KMIN
         WRITE(SCARC_LU,*) 'KMAX=', KMAX
         WRITE(SCARC_LU,'(a,i2,a,i2,a,i2,a,a,i2,a,i2,a,i2,a,i2,a,i2,a,i2,a,i2,a,i2,a)') 'LEVEL ',ILEVEL,&
                     ': allocate FACE-vectors for SCARC(', NM, ')%OSCARC(', NOM, ')',&
                     '%VEC(',IMIN,':',IMAX,',',JMIN,':',JMAX,',',KMIN,':',KMAX,')'
      ENDIF
    
   ENDDO FACE_NBR_LOOP2D
 
ENDDO COM_LEVEL_LOOP2D
 
! Allocate communication structures for all levels
DO ILEVEL=S%NLMAX,S%NLMIN,-1
   NREQ_FACE=0
   CALL SCARC_RECEIVE (NCOM_INIT, ILEVEL)    
   CALL SCARC_EXCHANGE(NCOM_INIT, ILEVEL, NMV_NONE, NMV_NONE)
ENDDO
 
TUSED_SCARC(3,NM)=TUSED_SCARC(3,NM)+SECOND()-TNOW_COMMUNICATION2D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_COMMUNICATION2D
END SUBROUTINE SCARC_INIT_COMMUNICATION2D
 

!------------------------------------------------------------------------------------
! determine neighborship/communication structure for data exchange in 3D
!------------------------------------------------------------------------------------
SUBROUTINE SCARC_INIT_COMMUNICATION3D (NM)
 
INTEGER :: IM, NM, NOM, ILEVEL, IREFINE
INTEGER :: IMIN, IMAX, JMIN, JMAX, KMIN, KMAX, ILMAX
INTEGER :: TAG_FACE
INTEGER :: IERR
REAL(EB):: TNOW_COMMUNICATION3D

TNOW_COMMUNICATION3D = SECOND()
 
S => SCARC (NM)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Allocate protocol arrays for data exchange
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ALLOCATE (STAT_FACE(MPI_STATUS_SIZE, NMESHES*NMESHES*20))
 
ALLOCATE (REQ_FACE(NMESHES*NMESHES*20))
REQ_FACE = MPI_REQUEST_NULL

ALLOCATE (TAGS_FACE(NMESHES, NMESHES))
TAGS_FACE = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Create unique tags arrays for face, edge and vertex exchanges
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
TAG_FACE = 0
 
DO IM = 1, NMESHES
   DO NOM = IM, NMESHES
      TAG_FACE = TAG_FACE + 1
      TAGS_FACE (IM, NOM) = TAG_FACE
      TAGS_FACE (NOM, IM) = TAG_FACE
   ENDDO
ENDDO 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Initialize level structures on neighboring meshes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IERR=0
OSLEVEL_LOOP3D: DO NOM=1,NMESHES

   ILMAX =  S%NLMAX
   SLMAX => S%SLEVEL(ILMAX)

   IF (SLMAX%NIC(NM,NOM)/=0.AND.SLMAX%NIC(NOM,NM)/=0) THEN

      OS => S%OSCARC(NOM)
      ALLOCATE (OS%SLEVEL(S%NLMIN:S%NLMAX), STAT=IERR)
      CALL CHKMEMERR ('SCARC_INIT_COMMUNICATION3D', 'OS%SLEVEL', IERR)

      OSLMAX => OS%SLEVEL(ILMAX)
      OSLMAX%N_EXTERNAL_WALL_CELLS = 2*S%MIBAR(NOM)*S%MJBAR(NOM) + &
                    2*S%MIBAR(NOM)*S%MKBAR(NOM) + &
                    2*S%MJBAR(NOM)*S%MKBAR(NOM)  
!      ALLOCATE (OSLMAX%IJKW(15,OSLMAX%N_EXTERNAL_WALL_CELLS), STAT=IERR)
!      CALL CHKMEMERR ('SCARC_INIT_COMMUNICATION3D', 'OSLMAX%IJKW', IERR)
!      OSLMAX%IJKW = 0

      IF (SCARC_DEBUG>=2) THEN
         WRITE(SCARC_LU,*) 'OSLMAX zeigt auf OS%SLEVEL(',ILMAX,')'
         WRITE(SCARC_LU,*) 'OSLMAX%N_EXTERNAL_WALL_CELLS=',OSLMAX%N_EXTERNAL_WALL_CELLS
         WRITE(SCARC_LU,*) 'IBAR=',S%MIBAR(NOM)
         WRITE(SCARC_LU,*) 'JBAR=',S%MJBAR(NOM)
         WRITE(SCARC_LU,*) 'KBAR=',S%MKBAR(NOM)
         !WRITE(SCARC_LU,*) 'SIZE(OSLMAX%IJKW)=',SIZE(OSLMAX%IJKW)
         WRITE(SCARC_LU,*) 'S%NLDIFF=',S%NLDIFF
      ENDIF

      IF (S%NLDIFF/=0) THEN
         IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,*) 'S%NLDIFF=',S%NLDIFF, S%NLMAX, S%NLMIN
         IREFINE=1
         DO ILEVEL=S%NLMAX-1,S%NLMIN,-1

            IREFINE=IREFINE*2

            IF (MOD(OS%SLEVEL(ILEVEL+1)%N_EXTERNAL_WALL_CELLS,IREFINE).NE.0) THEN
               WRITE(*,*) 'N_EXTERNAL_WALL_CELLS cannot be divided by 2!'         ! only temporarily
               STOP
            ENDIF

            OS%SLEVEL(ILEVEL)%N_EXTERNAL_WALL_CELLS=OS%SLEVEL(ILEVEL+1)%N_EXTERNAL_WALL_CELLS/IREFINE
            ALLOCATE (OS%SLEVEL(ILEVEL)%IJKW(15,OSLMAX%N_EXTERNAL_WALL_CELLS), STAT=IERR)
            CALL CHKMEMERR ('SCARC_INIT_COMMUNICATION3D', 'OS%SLEVEL%IJKW', IERR)
            OS%SLEVEL(ILEVEL)%IJKW = 0

            IF (SCARC_DEBUG>=2) THEN
               WRITE(SCARC_LU,*) 'ILEVEL=',ILEVEL
               WRITE(SCARC_LU,*) 'IREFINE=',IREFINE
               WRITE(SCARC_LU,*) 'OSL zeigt auf OS%SLEVEL(',ILEVEL,')'
               WRITE(SCARC_LU,*) 'OS%SLEVEL(',ILEVEL,')%N_EXTERNAL_WALL_CELLS=',OS%SLEVEL(ILEVEL)%N_EXTERNAL_WALL_CELLS
               WRITE(SCARC_LU,*) 'SIZE(OS%SLEVEL(ILEVEL)%IJKW)=',SIZE(OS%SLEVEL(ILEVEL)%IJKW)
            ENDIF

         ENDDO
      ENDIF

   ENDIF

ENDDO OSLEVEL_LOOP3D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Initialize communication structures on every level
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
COM_LEVEL_LOOP3D: DO ILEVEL=S%NLMAX,S%NLMIN,-1

   SL => S%SLEVEL(ILEVEL)

   !!! ---------------------------------------------------------------------------
   !!! Communication structures for faces
   !!! ---------------------------------------------------------------------------
   FACE_NBR_LOOP3D: DO NOM = 1, NMESHES

      ! is there a face neighbor ?
      IF (NOM == NM) CYCLE FACE_NBR_LOOP3D
      IF (SL%NIC(NM, NOM) == 0 .AND. SL%NIC(NOM, NM) == 0) CYCLE FACE_NBR_LOOP3D

      ! get length of communication vectors for different orientation of neighbors
      IMIN = SL%I_MIN(NM,NOM)
      IMAX = SL%I_MAX(NM,NOM)
      JMIN = SL%J_MIN(NM,NOM)
      JMAX = SL%J_MAX(NM,NOM)
      KMIN = SL%K_MIN(NM,NOM)
      KMAX = SL%K_MAX(NM,NOM)

      ! allocate communication vectors for data exchange with face neighbor NOM on level ILEVEL
      OSL =>S%OSCARC(NOM)%SLEVEL(ILEVEL)
   
      ALLOCATE (OSL%Y_FACE(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
      OSL%Y_FACE = 0.0_EB
    
      ALLOCATE (OSL%R_FACE(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
      OSL%R_FACE = 0.0_EB
    
      ALLOCATE (OSL%G_FACE(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
      OSL%G_FACE = 0.0_EB
   
      ALLOCATE (OSL%Z_FACE(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
      OSL%Z_FACE = 0.0_EB
   
      ! print information if requested    
      IF (SCARC_DEBUG >= 6) THEN
         WRITE(SCARC_LU,*) '---NM=', NM
         WRITE(SCARC_LU,*) '---NOM=', NOM
         WRITE(SCARC_LU,*) 'ILEVEL',ILEVEL
         WRITE(SCARC_LU,*) 'SL%I_MIN=', SL%I_MIN
         WRITE(SCARC_LU,*) 'SL%I_MAX=', SL%I_MAX
         WRITE(SCARC_LU,*) 'SL%K_MIN=', SL%K_MIN
         WRITE(SCARC_LU,*) 'SL%K_MAX=', SL%K_MAX
         WRITE(SCARC_LU,*) 'IMIN=', IMIN
         WRITE(SCARC_LU,*) 'IMAX=', IMAX
         WRITE(SCARC_LU,*) 'JMIN=', JMIN
         WRITE(SCARC_LU,*) 'JMAX=', JMAX
         WRITE(SCARC_LU,*) 'KMIN=', KMIN
         WRITE(SCARC_LU,*) 'KMAX=', KMAX
         WRITE(SCARC_LU,'(a,i2,a,i2,a,i2,a,a,i2,a,i2,a,i2,a,i2,a,i2,a,i2,a,i2,a,i2,a)') 'LEVEL ',ILEVEL,&
                     ': allocate FACE-vectors for SCARC(', NM, ')%OSCARC(', NOM, ')',&
                     '%VEC(',IMIN,':',IMAX,',',JMIN,':',JMAX,',',KMIN,':',KMAX,')'
      ENDIF
    
   ENDDO FACE_NBR_LOOP3D
 
ENDDO COM_LEVEL_LOOP3D
 
! Allocate communication structures for all levels
DO ILEVEL=S%NLMAX,S%NLMIN,-1
   NREQ_FACE=0
   CALL SCARC_RECEIVE (NCOM_INIT, ILEVEL)    
   CALL SCARC_EXCHANGE(NCOM_INIT, ILEVEL, NMV_NONE, NMV_NONE)
ENDDO
 
TUSED_SCARC(3,NM)=TUSED_SCARC(3,NM)+SECOND()-TNOW_COMMUNICATION3D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_COMMUNICATION3D
END SUBROUTINE SCARC_INIT_COMMUNICATION3D
 
 
 
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Allocate full matrix for the usual 5-point-stencil (2D)
!!! 'full' means that it corresponds to the 1-process matrix of the serial problem
!!!
!!! Note: The matrix is stored "band wise", namely diagonal for diagonal.
!!! Although the sub- and superdiagonals are shorter than the main diagonal,
!!! all bands are stored in the same length and are simply filled with zeros
!!! at redundant positions. The 'wasting' of this storage space is justified
!!! by the possibility to use a much more efficient matrix-vector-multiplication
!!! (which doesn't have to use an expensive referencing logic)
!!!
!!! 1. band:  ibar-th subdiagonal for the southern '1'-entry
!!! 2. band:  1-th subdiagonal for the western  '1'-entry
!!! 3. band:  main diagonal with the '-4'-entry
!!! 4. band:  1-th superdiagonal for the eastern  '1'-entry
!!! 5. band:  ibar-th superdiagonal for the northern '1'-entry
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  2D-version
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_INIT_MATRICES2D (NM)
 
INTEGER :: NM, NOM
INTEGER :: I, K, IC, IW, IERR, IOR, ILEVEL, BC_INDEX
REAL(EB):: DBC
LOGICAL :: BGSTRIX, BCOARSE
REAL(EB):: TNOW_MATRICES2D
TYPE (MESH_TYPE), POINTER :: M
 
TNOW_MATRICES2D = SECOND()

IERR = 0
NDIAG = 5
 
M => MESHES(NM)
 
!!!
!!! Initialize matrices on each level (only 1 level for CG)
!!!
MATRIX_LEVEL_LOOP2D: DO ILEVEL=S%NLMAX,S%NLMIN,-1

   SL => S%SLEVEL(ILEVEL)

   !!! Allocate matrix corresponding to the band-wise storage technique
   ALLOCATE (SL%AG(1:SL%NCELLS_LOCAL, 1:NDIAG), STAT=IERR)
   CALL CHKMEMERR ('SCARC_INIT_MATRICES2D', 'SL%AG', IERR)
   SL%AG = 0.0_EB
    
   SL%DXI2=1.0_EB/(SL%DX)**2
   SL%DZI2=1.0_EB/(SL%DZ)**2

   IF (SCARC_DEBUG >= 2) THEN
      WRITE(SCARC_LU,*) '========= MATRICES2D: LEVEL=',ILEVEL
      WRITE(SCARC_LU,*) ' IBAR=',SL%IBAR
      WRITE(SCARC_LU,*) ' KBAR=',SL%KBAR
      WRITE(SCARC_LU,*) ' DX  =',SL%DX  
      WRITE(SCARC_LU,*) ' DZ  =',SL%DZ  
      WRITE(SCARC_LU,*) ' DXI2=',SL%DXI2
      WRITE(SCARC_LU,*) ' DZI2=',SL%DZI2
   ENDIF
    
   DO K = 1, SL%KBAR
      DO I = 1, SL%IBAR
    
         IC = (K-1) * SL%IBAR + I
    
         SL%AG (IC, 1) = - 2.0_EB*SL%DXI2 - 2.0_EB*SL%DZI2  ! main diagonal

         IF (K > 1)       SL%AG(IC, 2) = SL%DZI2     ! lower subdiagonal corresponding to z-direction
         IF (I > 1)       SL%AG(IC, 3) = SL%DXI2     ! lower subdiagonal corresponding to x-direction

         IF (I < SL%IBAR) SL%AG(IC, 4) = SL%DXI2     ! upper subdiagonal corresponding to x-direction
         IF (K < SL%KBAR) SL%AG(IC, 5) = SL%DZI2     ! upper subdiagonal corresponding to z-direction
    
         IF (SCARC_DEBUG >= 6) THEN
            WRITE(SCARC_LU,*) 'D1: AG(', IC, ',1)=', SL%AG(IC, 1)
            WRITE(SCARC_LU,*) 'D2: AG(', IC, ',2)=', SL%AG(IC, 2)
            WRITE(SCARC_LU,*) 'D3: AG(', IC, ',3)=', SL%AG(IC, 3)
            WRITE(SCARC_LU,*) 'D4: AG(', IC, ',4)=', SL%AG(IC, 4)
            WRITE(SCARC_LU,*) 'D5: AG(', IC, ',5)=', SL%AG(IC, 5)
            WRITE(SCARC_LU,*) 
         ENDIF
    
      ENDDO
   ENDDO
    
   BGSTRIX= SCARC_CG_PRECON=='GSTRIX'.OR.SCARC_BICG_PRECON=='GSTRIX'.OR.SCARC_SM_PRECON=='GSTRIX'
   BCOARSE= SCARC_METHOD=='MG'.AND.ILEVEL==S%NLMIN

   IF (BGSTRIX.AND..NOT.BCOARSE) CALL SCARC_PRECON_INIT_GSTRIX2D(ILEVEL,NM)

   BC_INDEX_LOOP: DO IW = 1,SL%N_EXTERNAL_WALL_CELLS
   
      IF (ILEVEL==S%NLMAX) THEN
       
         IOR = M%IJKW(4,IW)

         IF (ABS(IOR)/=2) THEN
   
            I   = M%IJKW( 6,IW)
            K   = M%IJKW( 8,IW)
            NOM = M%IJKW( 9,IW)
            
            BC_INDEX=M%PRESSURE_BC_INDEX(IW)

         ELSE

           CYCLE BC_INDEX_LOOP

         ENDIF

      ELSE

         IOR = SL%IJKW( 4,IW)
         I   = SL%IJKW( 6,IW)
         K   = SL%IJKW( 8,IW)
         NOM = SL%IJKW( 9,IW)
   
         BC_INDEX=SL%PRESSURE_BC_INDEX(IW)
         
      ENDIF

      SELECT CASE(IOR)
         CASE (1)
            IC = (K-1) * SL%IBAR + I
            DBC= SL%DXI2
         CASE (-1)
            IC = K * SL%IBAR
            DBC= SL%DXI2
         CASE (3)
            IC = I
            DBC= SL%DZI2
         CASE (-3)
            IC = (SL%KBAR-1) * SL%IBAR + I
            DBC= SL%DZI2
      END SELECT

      SELECT CASE(BC_INDEX)
         CASE (DIRICHLET)                    ! set Dirichlet BC's along open boundary cells
            SL%AG(IC, 1) = SL%AG(IC, 1) - DBC
            IF (SCARC_DEBUG > 4) WRITE(SCARC_LU,1001) IW, IOR, IC, SL%AG(IC,1),M%BOUNDARY_TYPE(IW),NOM
         CASE (INTERNAL)                     ! do nothing along internal boundaries (just only debug message)
            IF (SCARC_DEBUG > 4) WRITE(SCARC_LU,1003) IW, IOR, IC, SL%AG(IC,1),M%BOUNDARY_TYPE(IW),NOM
         CASE (NEUMANN)                      ! set Neumann BC's at all other nodes
            SL%AG(IC, 1) = SL%AG(IC, 1) + DBC
            IF (SCARC_DEBUG > 4) WRITE(SCARC_LU,1002) IW, IOR, IC, SL%AG(IC,1),M%BOUNDARY_TYPE(IW),NOM
      END SELECT

   ENDDO BC_INDEX_LOOP

   SL%ASUBX = SL%DXI2
   SL%ASUBZ = SL%DZI2
 
   IF (SCARC_DEBUG >= 2) THEN
      !WRITE(SCARC_LU,'(16i4)') (M%BOUNDARY_TYPE(ICELL),ICELL=1,M%N_EXTERNAL_WALL_CELLS)
      !WRITE(SCARC_LU,*) 'IJKW'
      !DO IW=1,M%N_EXTERNAL_WALL_CELLS+4
      !   WRITE(SCARC_LU,'(i4,a,15i4)') IW,' : ',(M%IJKW(ICELL,IW),ICELL=1,15)
      !ENDDO
      WRITE(SCARC_LU,*) 'ASUB=', SL%ASUB
      WRITE(SCARC_LU,*)
      WRITE(SCARC_LU,*) 'GLOBAL MATRIX AG:'
      DO IC = 1, SL%NCELLS_LOCAL
         WRITE(SCARC_LU, '(5f12.4)') SL%AG(IC, 2),SL%AG(IC,3),SL%AG(IC,1),SL%AG(IC,4),SL%AG(IC,5)
         IF (Mod(IC, SL%IBAR) == 0) WRITE(SCARC_LU,*) '---------------------------------&
                                                &--------------------------------'
      ENDDO
      WRITE(SCARC_LU,*) 'SL%ASUBX=', SL%ASUBX
      WRITE(SCARC_LU,*) 'LBC=', LBC
      WRITE(SCARC_LU,*) 'MBC=', MBC
      WRITE(SCARC_LU,*) 'NBC=', NBC
   ENDIF

ENDDO MATRIX_LEVEL_LOOP2D
 

TUSED_SCARC(4,NM)=TUSED_SCARC(4,NM)+SECOND()-TNOW_MATRICES2D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_MATRICES2D

1001 FORMAT('IW=',i3,': Dirichlet , IOR=',i2,':AG(',i3,'1)=',f12.6,': BT=',i3,': NOM=',i3)
1002 FORMAT('IW=',i3,': Neumann   , IOR=',i2,':AG(',i3,'1)=',f12.6,': BT=',i3,': NOM=',i3)
1003 FORMAT('IW=',i3,': Nothing   , IOR=',i2,':AG(',i3,'1)=',f12.6,': BT=',i3,': NOM=',i3)
END SUBROUTINE SCARC_INIT_MATRICES2D
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  3D-version
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_INIT_MATRICES3D (NM)
 
INTEGER :: NM, NOM
INTEGER :: I, J, K, IC, IW, IERR, ICELL, IOR, ILEVEL, BC_INDEX
REAL(EB):: DBC
REAL(EB):: TNOW_MATRICES3D
TYPE (MESH_TYPE), POINTER :: M
 
TNOW_MATRICES3D = SECOND()

IERR = 0
NDIAG = 7

M => MESHES(NM)

IF (SCARC_DEBUG>=6) THEN
   WRITE(SCARC_LU,*) '============== Starting SCARC_INIT_MATRICES ============='
   WRITE(SCARC_LU,*) 'M%N_EXTERNAL_WALL_CELLS=',M%N_EXTERNAL_WALL_CELLS
   IF (NMESHES==1) THEN
   WRITE(SCARC_LU,*) 'BOUNDARY_TYPE:'
   WRITE(SCARC_LU,'(8i4)') (M%BOUNDARY_TYPE(IW),IW=1,64)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,'(8i4)') (M%BOUNDARY_TYPE(IW),IW=65,128)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,'(8i4)') (M%BOUNDARY_TYPE(IW),IW=129,192)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,'(8i4)') (M%BOUNDARY_TYPE(IW),IW=193,256)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,'(8i4)') (M%BOUNDARY_TYPE(IW),IW=257,320)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,'(8i4)') (M%BOUNDARY_TYPE(IW),IW=321,384)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,*) 'PRESSURE_BC_INDEX:'
   WRITE(SCARC_LU,'(8i4)') (M%PRESSURE_BC_INDEX(IW),IW=1,64)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,'(8i4)') (M%PRESSURE_BC_INDEX(IW),IW=65,128)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,'(8i4)') (M%PRESSURE_BC_INDEX(IW),IW=129,192)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,'(8i4)') (M%PRESSURE_BC_INDEX(IW),IW=193,256)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,'(8i4)') (M%PRESSURE_BC_INDEX(IW),IW=257,320)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,'(8i4)') (M%PRESSURE_BC_INDEX(IW),IW=321,384)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,*) '----------------------------------'
   ELSE
   WRITE(SCARC_LU,*) 'BOUNDARY_TYPE:'
   WRITE(SCARC_LU,'(8i4)') (M%BOUNDARY_TYPE(IW),IW=1,64)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,'(8i4)') (M%BOUNDARY_TYPE(IW),IW=65,128)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,'(4i4)') (M%BOUNDARY_TYPE(IW),IW=129,160)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,'(4i4)') (M%BOUNDARY_TYPE(IW),IW=161,192)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,'(4i4)') (M%BOUNDARY_TYPE(IW),IW=192,224)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,'(4i4)') (M%BOUNDARY_TYPE(IW),IW=225,256)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,*) 'PRESSURE_BC_INDEX:'
   WRITE(SCARC_LU,'(8i4)') (M%PRESSURE_BC_INDEX(IW),IW=1,64)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,'(8i4)') (M%PRESSURE_BC_INDEX(IW),IW=65,128)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,'(4i4)') (M%PRESSURE_BC_INDEX(IW),IW=129,160)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,'(4i4)') (M%PRESSURE_BC_INDEX(IW),IW=161,192)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,'(4i4)') (M%PRESSURE_BC_INDEX(IW),IW=193,224)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,'(4i4)') (M%PRESSURE_BC_INDEX(IW),IW=225,256)
   WRITE(SCARC_LU,*) '----------------------------------'
   WRITE(SCARC_LU,*) '----------------------------------'
   ENDIF
   DO IW=1,M%N_EXTERNAL_WALL_CELLS+4
      WRITE(SCARC_LU,'(i4,a,15i4)') IW,' : ',(M%IJKW(ICELL,IW),ICELL=1,15)
   ENDDO
ENDIF

 
!!!
!!! Initialize matrices on different grid levels (only 1 level for CG/BICG)
!!!
MATRIX_LEVEL_LOOP3D: DO ILEVEL=S%NLMAX,S%NLMIN,-1

   SL => S%SLEVEL(ILEVEL)

   !!! Allocate matrix corresponding to the band-wise storage technique
   ALLOCATE (SL%AG(1:SL%NCELLS_LOCAL, 1:NDIAG), STAT=IERR)
   CALL CHKMEMERR ('SCARC_INIT_MATRICES3D', 'SL%AG', IERR)
   SL%AG = 0.0_EB
    
   SL%DXI2=1.0_EB/(SL%DX)**2
   SL%DYI2=1.0_EB/(SL%DY)**2
   SL%DZI2=1.0_EB/(SL%DZ)**2

   IF (SCARC_DEBUG >= 2) THEN
      WRITE(SCARC_LU,*) '========= LEVEL=',ILEVEL
      WRITE(SCARC_LU,*) ' IBAR=',SL%IBAR
      WRITE(SCARC_LU,*) ' JBAR=',SL%JBAR
      WRITE(SCARC_LU,*) ' KBAR=',SL%KBAR
      WRITE(SCARC_LU,*) ' DX  =',SL%DX  
      WRITE(SCARC_LU,*) ' DY  =',SL%DY  
      WRITE(SCARC_LU,*) ' DZ  =',SL%DZ  
      WRITE(SCARC_LU,*) ' DXI2=',SL%DXI2
      WRITE(SCARC_LU,*) ' DYI2=',SL%DYI2
      WRITE(SCARC_LU,*) ' DZI2=',SL%DZI2
   ENDIF
    
   DO K = 1, SL%KBAR
      DO J = 1, SL%JBAR
         DO I = 1, SL%IBAR
    
            IC = (K-1) * SL%IBAR * SL%JBAR + (J-1) * SL%IBAR + I
    
            SL%AG (IC, 1) = - 2.0_EB*SL%DXI2 - 2.0*SL%DYI2 - 2.0_EB*SL%DZI2  ! main diagonal

            IF (K > 1)       SL%AG(IC, 2) = SL%DZI2     ! lower subdiagonal corresponding to z-direction
            IF (J > 1)       SL%AG(IC, 3) = SL%DYI2     ! lower subdiagonal corresponding to y-direction
            IF (I > 1)       SL%AG(IC, 4) = SL%DXI2     ! lower subdiagonal corresponding to x-direction

            IF (I < SL%IBAR) SL%AG(IC, 5) = SL%DXI2     ! upper subdiagonal corresponding to x-direction
            IF (J < SL%JBAR) SL%AG(IC, 6) = SL%DYI2     ! upper subdiagonal corresponding to y-direction
            IF (K < SL%KBAR) SL%AG(IC, 7) = SL%DZI2     ! upper subdiagonal corresponding to z-direction
    
            IF (SCARC_DEBUG >= 6) THEN
               WRITE(SCARC_LU,*) 'D1: AG(', IC, ',1)=', SL%AG(IC, 1)
               WRITE(SCARC_LU,*) 'D2: AG(', IC, ',2)=', SL%AG(IC, 2)
               WRITE(SCARC_LU,*) 'D3: AG(', IC, ',3)=', SL%AG(IC, 3)
               WRITE(SCARC_LU,*) 'D4: AG(', IC, ',4)=', SL%AG(IC, 4)
               WRITE(SCARC_LU,*) 'D5: AG(', IC, ',5)=', SL%AG(IC, 5)
               WRITE(SCARC_LU,*) 'D6: AG(', IC, ',6)=', SL%AG(IC, 6)
               WRITE(SCARC_LU,*) 'D7: AG(', IC, ',7)=', SL%AG(IC, 7)
               WRITE(SCARC_LU,*) 
            ENDIF
    
         ENDDO
      ENDDO
   ENDDO
    
   ! has still to be set correctly
   !IF (ILEVEL/=S%NLMIN) CALL SCARC_PRECON_INIT_GSTRIX3D(ILEVEL,NM)


   DO IW = 1,SL%N_EXTERNAL_WALL_CELLS

      IF (ILEVEL==S%NLMAX) THEN
   
         IOR = M%IJKW( 4,IW)
         I   = M%IJKW( 6,IW)
         J   = M%IJKW( 7,IW)
         K   = M%IJKW( 8,IW)
         NOM = M%IJKW( 9,IW)
         
         BC_INDEX=M%PRESSURE_BC_INDEX(IW)

      ELSE

         IOR = SL%IJKW( 4,IW)
         I   = SL%IJKW( 6,IW)
         J   = SL%IJKW( 7,IW)
         K   = SL%IJKW( 8,IW)
         NOM = SL%IJKW( 9,IW)
   
         BC_INDEX=SL%PRESSURE_BC_INDEX(IW)
         
      ENDIF

      SELECT CASE(IOR)
         CASE (1)
            IC = (K-1) * SL%IBAR * SL%JBAR + (J-1) * SL%IBAR + I
            DBC= SL%DXI2
         CASE (-1)
            IC = (K-1) * SL%IBAR * SL%JBAR + J * SL%IBAR 
            DBC= SL%DXI2
         CASE (2)
            IC = (K-1) * SL%IBAR * SL%JBAR + I
            DBC= SL%DYI2
         CASE (-2)
            IC = (K-1) * SL%IBAR * SL%JBAR + (SL%JBAR-1) * SL%IBAR + I
            DBC= SL%DYI2
         CASE (3)
            IC = (J-1) * SL%IBAR + I
            DBC= SL%DZI2
         CASE (-3)
            IC = (SL%KBAR-1) * SL%IBAR * SL%JBAR + (J-1) * SL%IBAR + I
            DBC= SL%DZI2
      END SELECT

      SELECT CASE(BC_INDEX)
         CASE (DIRICHLET)                    ! set Dirichlet BC's at open and null boundary cells
            SL%AG(IC, 1) = SL%AG(IC, 1) - DBC
            IF (SCARC_DEBUG > 4) WRITE(SCARC_LU,1001) IW, IOR, IC, SL%AG(IC,1),M%BOUNDARY_TYPE(IW),NOM
         CASE (INTERNAL)                     ! do nothing along internal boundaries (only debug message)
            IF (SCARC_DEBUG > 4) WRITE(SCARC_LU,1003) IW, IOR, IC, SL%AG(IC,1),M%BOUNDARY_TYPE(IW),NOM
         CASE (NEUMANN)                      ! set Neumann BC's at all other cells
            SL%AG(IC, 1) = SL%AG(IC, 1) + DBC
            IF (SCARC_DEBUG > 4) WRITE(SCARC_LU,1002) IW, IOR, IC, SL%AG(IC,1),M%BOUNDARY_TYPE(IW),NOM
      END SELECT

   ENDDO
   
   SL%ASUBX = SL%DXI2
   SL%ASUBY = SL%DYI2
   SL%ASUBZ = SL%DZI2

   IF (SCARC_DEBUG >= 2) THEN
      WRITE(SCARC_LU,*) 'ASUB=', SL%ASUB
      WRITE(SCARC_LU,*) 'NCELLS_LOCAL=', SL%NCELLS_LOCAL
      WRITE(SCARC_LU,*) 'N_EXTERNAL_WALL_CELLS=', SL%N_EXTERNAL_WALL_CELLS
      WRITE(SCARC_LU,*) 'IBAR=', SL%IBAR
      WRITE(SCARC_LU,*) 'JBAR=', SL%JBAR
      WRITE(SCARC_LU,*) 'KBAR=', SL%KBAR
      WRITE(SCARC_LU,*) 'ILEVEL=', ILEVEL
      WRITE(SCARC_LU,*)
      WRITE(SCARC_LU,*) 'GLOBAL MATRIX AG:'
      DO IC = 1, SL%NCELLS_LOCAL
        WRITE(SCARC_LU, '(7f12.4)') SL%AG(IC, 2),SL%AG(IC,3),SL%AG(IC,4), &
                                    SL%AG(IC,1), &
                                    SL%AG(IC,5),SL%AG(IC,6),SL%AG(IC,7)
        IF (Mod(IC, SL%JBAR*SL%IBAR) == 0) THEN
           WRITE(SCARC_LU,*) '===========================================&
                             &==========================================='
        ELSE IF (Mod(IC, SL%IBAR) == 0) THEN
           WRITE(SCARC_LU,*) '-------------------------------------------&
                             &-------------------------------------------'
        ENDIF
     ENDDO
     WRITE(SCARC_LU,*) 'SL%ASUBX=', SL%ASUBX
     WRITE(SCARC_LU,*) 'LBC=', LBC
     WRITE(SCARC_LU,*) 'MBC=', MBC
     WRITE(SCARC_LU,*) 'NBC=', NBC
   ENDIF

ENDDO MATRIX_LEVEL_LOOP3D
 
TUSED_SCARC(4,NM)=TUSED_SCARC(4,NM)+SECOND()-TNOW_MATRICES3D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_MATRICES3D

IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,*) '============== Leaving SCARC_INIT_MATRICES ============='

1001 FORMAT('IW=',i3,': Dirichlet , IOR=',i2,':AG(',i3,',1)=',f12.6,': BT=',i3,': NOM=',i3)
1002 FORMAT('IW=',i3,': Neumann   , IOR=',i2,':AG(',i3,',1)=',f12.6,': BT=',i3,': NOM=',i3)
1003 FORMAT('IW=',i3,': Nothing   , IOR=',i2,':AG(',i3,',1)=',f12.6,': BT=',i3,': NOM=',i3)
END SUBROUTINE SCARC_INIT_MATRICES3D


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Initialize global 2D-solver methods (cg/mg)
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_INIT_SOLVER2D (NM)
 
INTEGER :: NM, IERR, ILEVEL
INTEGER :: IBP0, KBP0
REAL(EB):: TNOW_SOLVER2D
 
TNOW_SOLVER2D = SECOND()
 
IERR = 0

S => SCARC (NM)

!!!
!!! Initialize SLEVEL-structures for single levels
!!!
LEVEL_LOOP2D: DO ILEVEL=S%NLMAX,S%NLMIN,-1

   SL => S%SLEVEL(ILEVEL)

   IBP0=SL%IBAR+1
   KBP0=SL%KBAR+1

   ! ----------------------------------------------------------------------------
   ! auxiliary vectors for update along internal boundaries
   ! ----------------------------------------------------------------------------
   ALLOCATE (SL%Z(0:IBP0, 0:2, 0:KBP0), STAT=IERR)
   CALL CHKMEMERR ('SCARC', 'Z', IERR)
   SL%Z = 0.0_EB
    
   ! ----------------------------------------------------------------------------
   ! working and auxiliary vectors for different solvers
   ! ----------------------------------------------------------------------------
   SELECT CASE(SCARC_METHOD)
      CASE('CG','BICG')

         ALLOCATE (SL%X(0:IBP0, 1, 0:KBP0), STAT=IERR)
         CALL CHKMEMERR ('SCARC', 'X', IERR)
         SL%X = 0.0_EB
    
         ALLOCATE (SL%F(0:IBP0, 1, 0:KBP0), STAT=IERR)
         CALL CHKMEMERR ('SCARC', 'B', IERR)
         SL%F = 0.0_EB
    
         ALLOCATE (SL%D(0:IBP0, 1, 0:KBP0), STAT=IERR)
         CALL CHKMEMERR ('SCARC', 'D', IERR)
         SL%D = 0.0_EB

         ALLOCATE (SL%R(0:IBP0, 1, 0:KBP0), STAT=IERR)
         CALL CHKMEMERR ('SCARC', 'R', IERR)
         SL%R = 0.0_EB
    
         ALLOCATE (SL%G(0:IBP0, 1, 0:KBP0), STAT=IERR)
         CALL CHKMEMERR ('SCARC', 'G', IERR)
         SL%G = 0.0_EB
    
         ALLOCATE (SL%Y(0:IBP0, 1, 0:KBP0), STAT=IERR)
         CALL CHKMEMERR ('SCARC', 'Y', IERR)
         SL%Y = 0.0_EB

         IF (SCARC_CG_PRECON == 'MG' .OR. SCARC_BICG_PRECON == 'MG') THEN
      
            ALLOCATE (SL%X2(0:IBP0, 1, 0:KBP0), STAT=IERR)
            CALL CHKMEMERR ('SCARC', 'X2', IERR)
            SL%X2 = 0.0_EB
      
            ALLOCATE (SL%F2(0:IBP0, 1, 0:KBP0), STAT=IERR)
            CALL CHKMEMERR ('SCARC', 'F2', IERR)
            SL%F2 = 0.0_EB
      
            ALLOCATE (SL%D2(0:IBP0, 1, 0:KBP0), STAT=IERR)
            CALL CHKMEMERR ('SCARC', 'D2', IERR)
            SL%D2 = 0.0_EB

            IF (ILEVEL==S%NLMIN) THEN
   
               ALLOCATE (SL%R2(0:IBP0, 1, 0:KBP0), STAT=IERR)
               CALL CHKMEMERR ('SCARC', 'R2', IERR)
               SL%R2 = 0.0_EB
       
               ALLOCATE (SL%G2(0:IBP0, 1, 0:KBP0), STAT=IERR)
               CALL CHKMEMERR ('SCARC', 'G2', IERR)
               SL%G2 = 0.0_EB
       
               ALLOCATE (SL%Y2(0:IBP0, 1, 0:KBP0), STAT=IERR)
               CALL CHKMEMERR ('SCARC', 'Y2', IERR)
               SL%Y2 = 0.0_EB
   
            ENDIF
          
         ENDIF

        IF (SCARC_CG_PRECON == 'FFT' .OR. SCARC_BICG_PRECON == 'FFT') THEN
            ALLOCATE (SL%FFT(1:IBP0, 1, 1:KBP0), STAT=IERR)
            CALL CHKMEMERR ('SCARC', 'FFT', IERR)
            SL%FFT = 0.0_EB
         ENDIF

      CASE('MG')

         ALLOCATE (SL%X(0:IBP0, 1, 0:KBP0), STAT=IERR)
         CALL CHKMEMERR ('SCARC', 'X', IERR)
         SL%X = 0.0_EB
    
         ALLOCATE (SL%F(0:IBP0, 1, 0:KBP0), STAT=IERR)
         CALL CHKMEMERR ('SCARC', 'B', IERR)
         SL%F = 0.0_EB
    
         ALLOCATE (SL%D(0:IBP0, 1, 0:KBP0), STAT=IERR)
         CALL CHKMEMERR ('SCARC', 'D', IERR)
         SL%D = 0.0_EB

         IF (ILEVEL==S%NLMIN) THEN

            ALLOCATE (SL%R(0:IBP0, 1, 0:KBP0), STAT=IERR)
            CALL CHKMEMERR ('SCARC', 'R', IERR)
            SL%R = 0.0_EB
    
            ALLOCATE (SL%G(0:IBP0, 1, 0:KBP0), STAT=IERR)
            CALL CHKMEMERR ('SCARC', 'G', IERR)
            SL%G = 0.0_EB
    
            ALLOCATE (SL%Y(0:IBP0, 1, 0:KBP0), STAT=IERR)
            CALL CHKMEMERR ('SCARC', 'Y', IERR)
            SL%Y = 0.0_EB

         ENDIF

         IF (SCARC_SM_PRECON == 'FFT') THEN
            ALLOCATE (SL%FFT(1:IBP0, 1, 1:KBP0), STAT=IERR)
            CALL CHKMEMERR ('SCARC', 'FFT', IERR)
            SL%FFT = 0.0_EB
         ENDIF

   END SELECT
    
ENDDO LEVEL_LOOP2D

    
!!!
!!! Allocate array for the description of the MG-cycle
!!!
IF (SCARC_METHOD == 'MG' .OR. SCARC_CG_PRECON =='MG' .OR. SCARC_BICG_PRECON == 'MG') THEN
   ALLOCATE (S%KCYCLE(2, NNLEVEL), STAT=IERR)
   CALL CHKMEMERR ('SCARC', 'KCYCLE', IERR)
   S%KCYCLE = 0
ENDIF
    
IF (SCARC_DEBUG >= 2) WRITE(SCARC_LU,*) 'INIT_SOLVER2D: FINISHED'

TUSED_SCARC(5,NM)=TUSED_SCARC(5,NM)+SECOND()-TNOW_SOLVER2D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_SOLVER2D
 
END SUBROUTINE SCARC_INIT_SOLVER2D

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Initialize global 3D-solver methods (cg/mg)
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_INIT_SOLVER3D (NM)
 
INTEGER :: NM, IERR, ILEVEL
INTEGER :: IBP0, JBP0, KBP0
REAL(EB):: TNOW_SOLVER3D
 
TNOW_SOLVER3D = SECOND()
 
IERR = 0

S => SCARC (NM)

!!!
!!! Initialize SLEVEL-structures for single levels
!!!
LEVEL_LOOP3D: DO ILEVEL=S%NLMAX,S%NLMIN,-1

   SL => S%SLEVEL(ILEVEL)

   IBP0=SL%IBAR+1
   JBP0=SL%JBAR+1
   KBP0=SL%KBAR+1

   ! ----------------------------------------------------------------------------
   ! auxiliary vectors for update along internal boundaries
   ! ----------------------------------------------------------------------------
   ALLOCATE (SL%Z(0:IBP0, 0:JBP0, 0:KBP0), STAT=IERR)
   CALL CHKMEMERR ('SCARC', 'Z', IERR)
   SL%Z = 0.0_EB
    
   ! ----------------------------------------------------------------------------
   ! working and auxiliary vectors for different solvers
   ! ----------------------------------------------------------------------------
   SELECT CASE(SCARC_METHOD)
      CASE('CG','BICG')

         ALLOCATE (SL%X(0:IBP0, 0:JBP0, 0:KBP0), STAT=IERR)
         CALL CHKMEMERR ('SCARC', 'X', IERR)
         SL%X = 0.0_EB
    
         ALLOCATE (SL%F(0:IBP0, 0:JBP0, 0:KBP0), STAT=IERR)
         CALL CHKMEMERR ('SCARC', 'B', IERR)
         SL%F = 0.0_EB
    
         ALLOCATE (SL%D(0:IBP0, 0:JBP0, 0:KBP0), STAT=IERR)
         CALL CHKMEMERR ('SCARC', 'D', IERR)
         SL%D = 0.0_EB

         ALLOCATE (SL%R(0:IBP0, 0:JBP0, 0:KBP0), STAT=IERR)
         CALL CHKMEMERR ('SCARC', 'R', IERR)
         SL%R = 0.0_EB
    
         ALLOCATE (SL%G(0:IBP0, 0:JBP0, 0:KBP0), STAT=IERR)
         CALL CHKMEMERR ('SCARC', 'G', IERR)
         SL%G = 0.0_EB
    
         ALLOCATE (SL%Y(0:IBP0, 0:JBP0, 0:KBP0), STAT=IERR)
         CALL CHKMEMERR ('SCARC', 'Y', IERR)
         SL%Y = 0.0_EB

         IF (SCARC_CG_PRECON == 'MG' .OR. SCARC_BICG_PRECON == 'MG') THEN
      
            ALLOCATE (SL%X2(0:IBP0, 0:JBP0, 0:KBP0), STAT=IERR)
            CALL CHKMEMERR ('SCARC', 'X2', IERR)
            SL%X2 = 0.0_EB
      
            ALLOCATE (SL%F2(0:IBP0, 0:JBP0, 0:KBP0), STAT=IERR)
            CALL CHKMEMERR ('SCARC', 'F2', IERR)
            SL%F2 = 0.0_EB
      
            ALLOCATE (SL%D2(0:IBP0, 0:JBP0, 0:KBP0), STAT=IERR)
            CALL CHKMEMERR ('SCARC', 'D2', IERR)
            SL%D2 = 0.0_EB

            IF (ILEVEL==S%NLMIN) THEN
   
               ALLOCATE (SL%R2(0:IBP0, 0:JBP0, 0:KBP0), STAT=IERR)
               CALL CHKMEMERR ('SCARC', 'R2', IERR)
               SL%R2 = 0.0_EB
       
               ALLOCATE (SL%G2(0:IBP0, 0:JBP0, 0:KBP0), STAT=IERR)
               CALL CHKMEMERR ('SCARC', 'G2', IERR)
               SL%G2 = 0.0_EB
       
               ALLOCATE (SL%Y2(0:IBP0, 0:JBP0, 0:KBP0), STAT=IERR)
               CALL CHKMEMERR ('SCARC', 'Y2', IERR)
               SL%Y2 = 0.0_EB
   
            ENDIF
          
         ENDIF

         IF (SCARC_CG_PRECON == 'FFT' .OR. SCARC_BICG_PRECON == 'FFT') THEN
            ALLOCATE (SL%FFT(1:IBP0, 1:JBP0, 1:KBP0), STAT=IERR)
            CALL CHKMEMERR ('SCARC', 'FFT', IERR)
            SL%FFT = 0.0_EB
         ENDIF

      CASE('MG')

         ALLOCATE (SL%X(0:IBP0, 0:JBP0, 0:KBP0), STAT=IERR)
         CALL CHKMEMERR ('SCARC', 'X', IERR)
         SL%X = 0.0_EB
    
         ALLOCATE (SL%F(0:IBP0, 0:JBP0, 0:KBP0), STAT=IERR)
         CALL CHKMEMERR ('SCARC', 'B', IERR)
         SL%F = 0.0_EB
    
         ALLOCATE (SL%D(0:IBP0, 0:JBP0, 0:KBP0), STAT=IERR)
         CALL CHKMEMERR ('SCARC', 'D', IERR)
         SL%D = 0.0_EB

         IF (ILEVEL==S%NLMIN) THEN

            ALLOCATE (SL%R(0:IBP0, 0:JBP0, 0:KBP0), STAT=IERR)
            CALL CHKMEMERR ('SCARC', 'R', IERR)
            SL%R = 0.0_EB
    
            ALLOCATE (SL%G(0:IBP0, 0:JBP0, 0:KBP0), STAT=IERR)
            CALL CHKMEMERR ('SCARC', 'G', IERR)
            SL%G = 0.0_EB
    
            ALLOCATE (SL%Y(0:IBP0, 0:JBP0, 0:KBP0), STAT=IERR)
            CALL CHKMEMERR ('SCARC', 'Y', IERR)
            SL%Y = 0.0_EB

         ENDIF

         IF (SCARC_SM_PRECON == 'FFT') THEN
            ALLOCATE (SL%FFT(1:IBP0, 1:JBP0, 1:KBP0), STAT=IERR)
            CALL CHKMEMERR ('SCARC', 'FFT', IERR)
            SL%FFT = 0.0_EB
         ENDIF


   END SELECT
    
ENDDO LEVEL_LOOP3D

    
!!!
!!! Allocate array for the description of the MG-cycle
!!!
IF (SCARC_METHOD == 'MG' .OR. SCARC_CG_PRECON =='MG' .OR. SCARC_BICG_PRECON == 'MG') THEN
   ALLOCATE (S%KCYCLE(2, NNLEVEL), STAT=IERR)
   CALL CHKMEMERR ('SCARC', 'KCYCLE', IERR)
   S%KCYCLE = 0
ENDIF
    
IF (SCARC_DEBUG >= 2) WRITE(SCARC_LU,*) 'INIT_SOLVER: FINISHED'
 
TUSED_SCARC(5,NM)=TUSED_SCARC(5,NM)+SECOND()-TNOW_SOLVER3D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_SOLVER3D

END SUBROUTINE SCARC_INIT_SOLVER3D
 

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Perform global 2D cg-method based on global possion-matrix
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_CG2D (NM, HP, PRHS)
 
INTEGER ::  NM, ILMAX, ITE, ITE0, I, K
REAL (EB), POINTER, DIMENSION (:, :, :) :: HP, PRHS
REAL (EB) :: SIGMA0, SIGMA1, ALPHA, GAMMA
REAL (EB) :: TNOW_CG2D
LOGICAL BCONV, BDIVG
TYPE (MESH_TYPE),  POINTER :: M
 
TNOW_CG2D = SECOND()
BCONV=.FALSE.
BDIVG=.FALSE.

S  => SCARC (NM)
M  => MESHES (NM)

! cg only works in finest grid level
ILMAX =  S%NLMAX                              
SL    => S%SLEVEL(ILMAX)
 
! re-initialize auxiliary vectors
SL%G = 0.0_EB                   
SL%Y = 0.0_EB
SL%R = 0.0_EB
SL%D = 0.0_EB
SL%X = 0.0_EB
SL%F = 0.0_EB

! get right hand side and initial vector
DO K = 1, SL%KBAR                                     ! 
   DO I = 1, SL%IBAR
      SL%F (I, 1, K) = PRHS (I, 1, K)
      SL%X (I, 1, K) = HP (I, 1, K)
   ENDDO
ENDDO

IF (SCARC_DEBUG>=6) THEN
   CALL SCARC_SHOW0(SL%X,'SCARC_CG3D','XCG before all')
   CALL SCARC_SHOW0(SL%F,'SCARC_CG3D','FCG before all')
ENDIF

! set boundary conditions along exterior boundaries
CALL SCARC_SETBDRY2D (ILMAX, NM)
 
!!! calculate initial defect Y = B - A*X and get l2-norm of it
NREQ_FACE=0
CALL SCARC_MATVEC2D (SL%X, SL%R, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILMAX, NMV_X, NMV_R)
SL%R = -SL%F + SL%R
SL%RESIN_CG = SCARC_L2NORM2D (SL%R, ILMAX, NM, NTYPE_GLOBAL)

IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) 0, SL%RESIN_CG
IF (SCARC_DEBUG >= 1.AND.NM==1) write (*,1000) 0, SL%RESIN_CG
 

! initial preconditioning
SELECT CASE(SCARC_CG_PRECON)
   CASE('NONE') 
      SIGMA0 = SL%RESIN_CG * SL%RESIN_CG
   CASE('JACOBI') 
      SL%G=SL%R
      CALL SCARC_PRECON_JACOBI2D (SL%G, ILMAX, NM)
      SIGMA0 = SCARC_SCALPROD2D (SL%R, SL%G, ILMAX, NM)
   CASE('SSOR') 
      SL%G=SL%R
      CALL SCARC_PRECON_SSOR2D (SL%G, ILMAX, NM)
      SIGMA0 = SCARC_SCALPROD2D (SL%R, SL%G, ILMAX, NM)
   CASE('GSTRIX') 
      SL%G=SL%R
      CALL SCARC_PRECON_GSTRIX2D (SL%G, ILMAX, NM)
      SIGMA0 = SCARC_SCALPROD2D (SL%R, SL%G, ILMAX, NM)
   CASE('MG')
      SL%G=SL%R
      CALL SCARC_PRECON_MG2D (SL%G, ILMAX, NM)
      SIGMA0 = SCARC_SCALPROD2D (SL%R, SL%G, ILMAX, NM)
   CASE('FFT')
      CALL SCARC_PRECON_FFT2D (SL%R, SL%G, ILMAX, NM)
      SIGMA0 = SCARC_SCALPROD2D (SL%R, SL%G, ILMAX, NM)
END SELECT
SL%D = -SL%G
 
!
! defect correction loop
!
CG_LOOP2D: DO ITE = 1, SCARC_CG_NIT
 
   ! calculate new defect and get L2-norm of it
   NREQ_FACE=0
   CALL SCARC_MATVEC2D (SL%D, SL%Y, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILMAX, NMV_D, NMV_Y)
   ALPHA = SCARC_SCALPROD2D (SL%D, SL%Y, ILMAX, NM)

   ! get new descent direction
   ALPHA = SIGMA0 / ALPHA
   SL%X = ALPHA * SL%D + SL%X
   SL%R = ALPHA * SL%Y + SL%R
   SL%RES_CG = SCARC_L2NORM2D (SL%R, ILMAX, NM, NTYPE_GLOBAL)
 
IF (SCARC_DEBUG>=6) THEN
   WRITE(SCARC_LU,*) 'ALPHA=',ALPHA
   CALL SCARC_SHOW0(SL%X,'SCARC_CG3D','X in new loop ')
   CALL SCARC_SHOW0(SL%R,'SCARC_CG3D','R in new loop ')
ENDIF

   !  exit loop in case of convergence or divergence
   IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) ITE, SL%RES_CG
   IF (SCARC_DEBUG >= 1.AND.NM==1) write (*,1000) ITE, SL%RES_CG
   IF (SCARC_BREL) THEN
      IF (SL%RES_CG <= SL%RESIN_CG*SCARC_CG_EPS) BCONV = .TRUE.
   ELSE
      IF (SL%RES_CG <= SCARC_CG_EPS .AND. SL%RES_CG <= SL%RESIN_CG*SCARC_EPS_REL) BCONV = .TRUE.
   ENDIF
   IF (SL%RES_CG > SCARC_EPS_DIVG) BDIVG = .TRUE.
   IF (BCONV.OR.BDIVG) EXIT CG_LOOP2D
 
   ! preconditioning
   SELECT CASE(SCARC_CG_PRECON)
      CASE('NONE')
         SIGMA1 = SL%RES_CG * SL%RES_CG
      CASE('JACOBI')
         SL%G = SL%R
         CALL SCARC_PRECON_JACOBI2D (SL%G, ILMAX, NM)
         SIGMA1 = SCARC_SCALPROD2D (SL%R, SL%G, ILMAX, NM)
      CASE('SSOR')
         SL%G = SL%R
         CALL SCARC_PRECON_SSOR2D (SL%G, ILMAX, NM)
         SIGMA1 = SCARC_SCALPROD2D (SL%R, SL%G, ILMAX, NM)
      CASE('GSTRIX')
         SL%G = SL%R
         CALL SCARC_PRECON_GSTRIX2D (SL%G, ILMAX, NM)
         SIGMA1 = SCARC_SCALPROD2D (SL%R, SL%G, ILMAX, NM)
      CASE('MG')
         SL%G = SL%R
         CALL SCARC_PRECON_MG2D (SL%G, ILMAX, NM)
         SIGMA1 = SCARC_SCALPROD2D (SL%R, SL%G, ILMAX, NM)
      CASE('FFT')
         CALL SCARC_PRECON_FFT2D (SL%R, SL%G, ILMAX, NM)
         SIGMA1 = SCARC_SCALPROD2D (SL%R, SL%G, ILMAX, NM)
   END SELECT
 
   GAMMA = SIGMA1 / SIGMA0
   SIGMA0 = SIGMA1
   SL%D = -SL%G + GAMMA * SL%D
 
IF (SCARC_DEBUG>=6) THEN
   WRITE(SCARC_LU,*) 'SIGMA0=',SIGMA0
   WRITE(SCARC_LU,*) 'SIGMA1=',SIGMA1
   WRITE(SCARC_LU,*) 'GAMMA=',GAMMA
   CALL SCARC_SHOW0(SL%R,'SCARC_CG3D','R at end of loop ')
   CALL SCARC_SHOW0(SL%G,'SCARC_CG3D','G at end of loop ')
   CALL SCARC_SHOW0(SL%X,'SCARC_CG3D','X at end of loop ')
ENDIF

ENDDO CG_LOOP2D
 
! divergence or convergence ?
IF (BDIVG) THEN                       
   ITE0 = - 1
   SL%CAPPA_CG = 1.0_EB
ELSE 
   IF (BCONV) THEN
      ITE0=ITE
   ELSE
      ITE0=ITE-1
   ENDIF
   IF (SL%RESIN_CG >= 1.0E-70_EB) THEN
      SL%RESIN_CG = SL%RES_CG / SL%RESIN_CG
      SL%CAPPA_CG = SL%RESIN_CG ** (1.0_EB/ITE0)
   ELSE 
      SL%CAPPA_CG = 0.0_EB
   ENDIF
ENDIF
SCARC_NIT=ITE0
SCARC_RES=SL%RES_CG
SCARC_CAPPA=SL%CAPPA_CG
IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,2000) ITE0, SL%RES_CG, SL%CAPPA_CG 
IF (SCARC_DEBUG >= 1.AND.NM==1) WRITE(*,2000) ITE0, SL%RES_CG, SL%CAPPA_CG 

! shift solution back to vector HP
DO K = 1, SL%KBAR                                 
   DO I = 1, SL%IBAR
      HP(I, 1, K) = SL%X (I, 1, K)
   ENDDO
ENDDO

CALL SCARC_SHOW0(HP,'SCARC_CG2D','HP before GHOSTCELLS')

! set ghost cell values along exterior boundaries 
CALL SCARC_GHOSTCELLS(HP,NM)

CALL SCARC_SHOW0(HP,'SCARC_CG2D','HP after GHOSTCELLS')

! get neighbouring data along interior boundaries
!CALL SCARC_UPDATE_LEVEL(HP, NM, ILMAX, 'HPFIN ')
CALL SCARC_UPDATE(HP, NM, 'HPFIN ')

CALL SCARC_SHOW0(HP,'SCARC_CG2D','HP after UPDATE')

TUSED_SCARC(6,NM)=TUSED_SCARC(6,NM)+SECOND()-TNOW_CG2D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_CG2D

1000 FORMAT ('=====SCARC_CG       : #ite= ',i4,': res=',e16.8)
2000 FORMAT ('=====SCARC_CG       : #ite= ',i4,': res=',e16.8,': rate=',e16.8)
END SUBROUTINE SCARC_CG2D


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Perform global 2D cg-method based on global possion-matrix
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_BICG2D (NM, HP, PRHS)
 
INTEGER ::  NM, ILMAX, ITE, ITE0, I, K
REAL (EB), POINTER, DIMENSION (:, :, :) :: HP, PRHS
REAL (EB) :: RHO0, RHO1, DBETA, DTHETA, ALPHA0, ALPHA1, ALPHA2
REAL (EB) :: TNOW_BICG2D
LOGICAL BCONV, BDIVG
TYPE (MESH_TYPE),  POINTER :: M
 
TNOW_BICG2D = SECOND()
BCONV=.FALSE.
BDIVG=.FALSE.

S  => SCARC (NM)
M  => MESHES (NM)

IF (SCARC_DEBUG >=3) WRITE(SCARC_LU,*) 'STARTING BICG2D'

! cg only works in finest grid level
ILMAX =  S%NLMAX                              
SL    => S%SLEVEL(ILMAX)
 
! get right hand side and initial vector
DO K = 1, SL%KBAR                                     ! 
   DO I = 1, SL%IBAR
      SL%F (I, 1, K) = PRHS (I, 1, K)
      SL%X (I, 1, K) = HP (I, 1, K)
   ENDDO
ENDDO

! preset some variables
RHO0  =1.0_EB
RHO1  =0.0_EB
DBETA =0.0_EB
DTHETA=1.0_EB
ALPHA0=1.0_EB

! re-initialize auxiliary vectors
SL%R = 0.0_EB
SL%G = 0.0_EB                   
SL%Z = 0.0_EB
SL%Y = 0.0_EB
SL%D = 0.0_EB
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! initialization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CALL SCARC_SETBDRY2D (ILMAX, NM)
SL%R=SL%F

SELECT CASE(SCARC_BICG_PRECON)
   CASE('JACOBI') 
      CALL SCARC_PRECON_JACOBI2D (SL%R, ILMAX, NM)
   CASE('SSOR') 
      CALL SCARC_PRECON_SSOR2D (SL%R, ILMAX, NM)
   CASE('GSTRIX') 
      CALL SCARC_PRECON_GSTRIX2D (SL%R, ILMAX, NM)
   CASE('MG') 
      CALL SCARC_PRECON_MG2D (SL%R, ILMAX, NM)
   CASE('FFT') 
      CALL SCARC_PRECON_FFT2D (SL%R, SL%R, ILMAX, NM)
END SELECT

CALL SCARC_MATVEC2D (SL%X, SL%R, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILMAX, NMV_X, NMV_R)
SL%R = SL%F - SL%R

SELECT CASE(SCARC_BICG_PRECON)
   CASE('JACOBI') 
      CALL SCARC_PRECON_JACOBI2D (SL%R, ILMAX, NM)
   CASE('SSOR') 
      CALL SCARC_PRECON_SSOR2D (SL%R, ILMAX, NM)
   CASE('GSTRIX') 
      CALL SCARC_PRECON_GSTRIX2D (SL%R, ILMAX, NM)
   CASE('MG') 
      CALL SCARC_PRECON_MG2D (SL%R, ILMAX, NM)
   CASE('FFT') 
      CALL SCARC_PRECON_FFT2D (SL%R, SL%R, ILMAX, NM)
END SELECT

SL%RESIN_BICG = SCARC_L2NORM2D (SL%R, ILMAX, NM, NTYPE_GLOBAL)
SL%G=SL%R

IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) 0, SL%RESIN_BICG
IF (SCARC_DEBUG >= 1.AND.NM==1) write (*,1000) 0, SL%RESIN_BICG


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! iterative correction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
BICG_LOOP2D: DO ITE = 1, SCARC_BICG_NIT
 
   RHO1 = SCARC_SCALPROD2D (SL%G, SL%R, ILMAX, NM)
   DBETA=(RHO1*DTHETA)/(RHO0*ALPHA0)

   RHO0=RHO1

   SL%Z = SL%R + DBETA * SL%Z
   SL%Z = -DBETA * ALPHA0 * SL%Y + SL%Z

   ! matrix-vector multiplication and preconditioning
   CALL SCARC_MATVEC2D (SL%Z, SL%Y, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILMAX, NMV_Z, NMV_Y)

   SELECT CASE(SCARC_BICG_PRECON)
      CASE('JACOBI') 
         CALL SCARC_PRECON_JACOBI2D (SL%Y, ILMAX, NM)
      CASE('SSOR') 
         CALL SCARC_PRECON_SSOR2D (SL%Y, ILMAX, NM)
      CASE('GSTRIX') 
         CALL SCARC_PRECON_GSTRIX2D (SL%Y, ILMAX, NM)
      CASE('MG') 
         CALL SCARC_PRECON_MG2D (SL%Y, ILMAX, NM)
      CASE('FFT') 
         CALL SCARC_PRECON_FFT2D (SL%Y, SL%Y, ILMAX, NM)
   END SELECT

   DTHETA = SCARC_SCALPROD2D (SL%G, SL%Y, ILMAX, NM)
   DTHETA=RHO1/DTHETA
   SL%R = -DTHETA * SL%Y + SL%R

   ! matrix-vector multiplication and preconditioning
   CALL SCARC_MATVEC2D (SL%R, SL%D, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILMAX, NMV_R, NMV_D)

   SELECT CASE(SCARC_BICG_PRECON)
      CASE('JACOBI') 
         CALL SCARC_PRECON_JACOBI2D (SL%D, ILMAX, NM)
      CASE('SSOR') 
         CALL SCARC_PRECON_SSOR2D (SL%D, ILMAX, NM)
      CASE('GSTRIX') 
         CALL SCARC_PRECON_GSTRIX2D (SL%D, ILMAX, NM)
      CASE('MG') 
         CALL SCARC_PRECON_MG2D (SL%D, ILMAX, NM)
      CASE('FFT') 
         CALL SCARC_PRECON_FFT2D (SL%D, SL%D, ILMAX, NM)
   END SELECT

   ALPHA1 = SCARC_SCALPROD2D (SL%D, SL%R, ILMAX, NM)
   ALPHA2 = SCARC_SCALPROD2D (SL%D, SL%D, ILMAX, NM)

   ALPHA0=ALPHA1/ALPHA2

   SL%X =  DTHETA * SL%Z + SL%X
   SL%X =  ALPHA0 * SL%R + SL%X
   SL%R = -ALPHA0 * SL%D + SL%R

   SL%RES_BICG = SCARC_L2NORM2D (SL%R, ILMAX, NM, NTYPE_GLOBAL)
 
   !  exit loop in case of convergence or divergence
   IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) ITE, SL%RES_BICG
   IF (SCARC_DEBUG >= 1.AND.NM==1) write (*,1000) ITE, SL%RES_BICG
   IF (SCARC_BREL) THEN
      IF (SL%RES_BICG <= SL%RESIN_BICG*SCARC_BICG_EPS) BCONV = .TRUE.
   ELSE
      IF (SL%RES_BICG <= SCARC_BICG_EPS .AND. SL%RES_BICG <= SL%RESIN_BICG*SCARC_EPS_REL) BCONV = .TRUE.
   ENDIF
   IF (SL%RES_BICG > SCARC_EPS_DIVG) BDIVG = .TRUE.
   IF (BCONV.OR.BDIVG) EXIT BICG_LOOP2D
 
ENDDO BICG_LOOP2D
 
! divergence or convergence ?
IF (BDIVG) THEN                       
   ITE0 = - 1
   SL%CAPPA_BICG = 1.0_EB
ELSE 
   IF (BCONV) THEN
      ITE0=ITE
   ELSE
      ITE0=ITE-1
   ENDIF
   IF (SL%RESIN_BICG >= 1.0E-70_EB) THEN
      SL%RESIN_BICG = SL%RES_BICG / SL%RESIN_BICG
      SL%CAPPA_BICG = SL%RESIN_BICG ** (1.0_EB/ITE0)
   ELSE 
      SL%CAPPA_BICG = 0.0_EB
   ENDIF
ENDIF
SCARC_NIT=ITE0
SCARC_RES=SL%RES_BICG
SCARC_CAPPA=SL%CAPPA_BICG
IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,2000) ITE0, SL%RES_BICG, SL%CAPPA_BICG 
IF (SCARC_DEBUG >= 1.AND.NM==1) WRITE(*,2000) ITE0, SL%RES_BICG, SL%CAPPA_BICG 

! shift solution back to vector HP
DO K = 1, SL%KBAR                                 
   DO I = 1, SL%IBAR
      HP(I, 1, K) = SL%X (I, 1, K)
   ENDDO
ENDDO

! set ghost cell values along exterior boundaries 
CALL SCARC_GHOSTCELLS(HP,NM)

! get neighbouring data along interior boundaries
CALL SCARC_UPDATE_LEVEL(HP, NM, ILMAX, 'HPFIN ')

TUSED_SCARC(6,NM)=TUSED_SCARC(6,NM)+SECOND()-TNOW_BICG2D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_BICG2D

1000 FORMAT ('=====SCARC_BICG     : #ite= ',i4,': res=',e16.8)
2000 FORMAT ('=====SCARC_BICG     : #ite= ',i4,': res=',e16.8,': rate=',e16.8)
END SUBROUTINE SCARC_BICG2D



 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Perform global MG-method based on global possion-matrix
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_MG2D (NM, HP, PRHS)
 
INTEGER :: NM, ITE, ITE0, ICYCLE
REAL (EB), POINTER, DIMENSION (:, :, :) :: HP, PRHS
INTEGER :: ILEVEL, ILMAX
INTEGER :: I, K
REAL (EB) :: TNOW_MG2D
LOGICAL BMATVEC, BCONV, BDIVG
!TYPE (SCARC_LEVEL_TYPE), POINTER :: SL, SLHI, SLLO
TYPE (MESH_TYPE),  POINTER :: M
 
TNOW_MG2D = SECOND()

! initialize working vectors
S    => SCARC (NM)
M    => MESHES (NM)

ILMAX =  S%NLMAX
SLHI  => S%SLEVEL(ILMAX)
 
! get right hand side and initial vector
DO K = 1, SLHI%KBAR
   DO I = 1, SLHI%IBAR
      SLHI%F (I, 1, K) = PRHS (I, 1, K)
      SLHI%X (I, 1, K) = HP (I, 1, K)
   ENDDO
ENDDO
 
IF (SCARC_DEBUG >= 2) THEN
   WRITE(SCARC_LU,*) 'Starting SCARC_MG2D'
   WRITE(SCARC_LU,*) 'IBAR=:', SLHI%IBAR
   WRITE(SCARC_LU,*) 'JBAR=:', SLHI%JBAR
   WRITE(SCARC_LU,*) 'KBAR=:', SLHI%KBAR
   WRITE(SCARC_LU,*) 'SCARC_MG_NIT=:', SCARC_MG_NIT
   WRITE(SCARC_LU,*) 'SCARC_MG_EPS=:', SCARC_MG_EPS
   WRITE(SCARC_LU,*) 'SCARC_EPS_REL=:', SCARC_EPS_REL
   CALL SCARC_SHOW(SLHI%X, 'SARC',  'X1    ')
   CALL SCARC_SHOW(SLHI%F, 'SARC',  'F1    ')
ENDIF
 
!!! initialize some method parameters

BCONV =.FALSE.
BDIVG =.FALSE.
ICYCLE=1      ! V-cycle
 
!!! adjust boundary values of right hand side
CALL SCARC_SETBDRY2D (ILMAX, NM)
 
!!! ------------------------------------------------------------------------------------
!!! Only one level: solve problem exactly
!!! ------------------------------------------------------------------------------------
ONLY_ONE_LEVEL_IF: IF (S%NLEVEL==1) THEN   

      CALL SCARC_COARSE2D(ILMAX, NM)

!!! ------------------------------------------------------------------------------------
!!! More than one level: start MG-cycling
!!! ------------------------------------------------------------------------------------
ELSE

   !!! save cycle counts for MG-iteration
   S%KCYCLE(2,S%NLMAX)=1
   DO ILEVEL=S%NLMIN+1,S%NLMAX-1
      IF (ICYCLE==0) THEN
         S%KCYCLE(2,ILEVEL)=2
      ELSE
         S%KCYCLE(2,ILEVEL)=ICYCLE
      ENDIF
   ENDDO
     
   !!! calculate initial defect on finest level, Y = B - A*X, and get l2-norm of it
   ILEVEL=ILMAX
   CALL SCARC_MATVEC2D (SLHI%X, SLHI%D, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, NMV_X, NMV_D)       ! 3 korrekt ???
   SLHI%D = SLHI%F - SLHI%D

   SLHI%RESIN_MG = SCARC_L2NORM2D (SLHI%D , ILEVEL, NM, NTYPE_GLOBAL)
   IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) 0, ILMAX, SL%RESIN_MG
   IF (SCARC_DEBUG >= 1.AND.NM==1) write (*,1000) 0, ILMAX, SL%RESIN_MG

   !CALL SCARC_MASTER_INFO (SLHI%ITE, SLHI%RESIN_MG)

    
   !!! ---------------------------------------------------------------------------
   !!! start MG-iteration
   !!! ---------------------------------------------------------------------------
   MG_LOOP2D: DO ITE = 1, SCARC_MG_NIT
    
      NREQ_FACE = 0        

      !!! set level-information
      DO ILEVEL=S%NLMIN,S%NLMAX
         S%KCYCLE(1,ILEVEL)=S%KCYCLE(2,ILEVEL)
      ENDDO

      ILEVEL=S%NLMAX

      !!! ---------------------------------------------------------------------------
      !!!  Presmoothing
      !!! ---------------------------------------------------------------------------
110   IF (ILEVEL/=S%NLMIN) THEN
    
         SLHI => S%SLEVEL(ILEVEL)
         SLLO => S%SLEVEL(ILEVEL-1)

         !IF (SCARC_METHOD>=2) THEN
         !   WRITE(SCARC_LU,*) 'SLHI zeigt auf level ',ILEVEL
         !   WRITE(SCARC_LU,*) 'SLLO zeigt auf level ',ILEVEL-1
         !ENDIF

         IF (ILEVEL==S%NLMAX) THEN
            BMATVEC=.FALSE.
         ELSE
            BMATVEC=.TRUE.
         ENDIF

         CALL SCARC_SMOOTHER2D(SLHI%X, SLHI%F, SLHI%D, BMATVEC, ILEVEL, NM)

         !!! print level information - still missing
         !CALL SCARC_INFO_LEVEL(ILEVEL)

         !!! perform restriction to coarser grid: F:= rest(D)
         CALL SCARC_RESTRICTION2D(SLLO%F, SLHI%D, ILEVEL-1, NM)

         !!! initialize solution vector and set boundary conditions to residuum - necessary ?
         SLLO%X=0.0_EB
         !CALL SCARC_BDRY_RESIDUUM(NM, SLLO%F, ILEVEL-1)

         !!! decrease level
         ILEVEL=ILEVEL-1

         GOTO 110
      
      ENDIF         

      !!! ------------------------------------------------------------------------
      !!! Coarse grid solver
      !!! ------------------------------------------------------------------------
      ILEVEL=S%NLMIN
      CALL SCARC_COARSE2D(ILEVEL, NM)


      !!! ------------------------------------------------------------------------
      !!! Postsmoothing
      !!! ------------------------------------------------------------------------
120   IF (ILEVEL/=S%NLMAX) THEN
    
         BMATVEC=.TRUE.

         SLLO => S%SLEVEL(ILEVEL)

         ILEVEL=ILEVEL+1
         SLHI => S%SLEVEL(ILEVEL)

         !!! print level information - still missing
         !CALL SCARC_INFO_LEVEL(ILEVEL)

         !!! perform prolongation to finer grid: D:=prol(X)
         CALL SCARC_PROLONGATION2D(SLLO%X, SLHI%D, ILEVEL, NM)

         !!! set exterior boundary data of residuum to zero - necessary ?
         !CALL SCARC_BDRY_RESIDUUM(NM, SLHI%D, ILEVEL)
 
         !!! set new solution
         SLHI%X = SLHI%D + SLHI%X

         !!! select smoother for postsmoothing
         CALL SCARC_SMOOTHER2D(SLHI%X, SLHI%F, SLHI%D, BMATVEC, ILEVEL, NM)


         !!! save cycle counts
         S%KCYCLE(1,ILEVEL)=S%KCYCLE(1,ILEVEL)-1
 
         IF (S%KCYCLE(1,ILEVEL)==0) THEN
            IF (ICYCLE==0) THEN
               S%KCYCLE(1,ILEVEL)=1
            ELSE
               S%KCYCLE(1,ILEVEL)=S%KCYCLE(2,ILEVEL)
            ENDIF
            GOTO 120
         ELSE
            GOTO 110
         ENDIF

      ENDIF

      !!! ------------------------------------------------------------------------
      !!! Defect calculation
      !!! ------------------------------------------------------------------------
      CALL SCARC_MATVEC2D (SLHI%X, SLHI%D, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, NMV_X, NMV_D)
      SLHI%D = SLHI%F - SLHI%D
      SLHI%RES_MG = SCARC_L2NORM2D (SLHI%D, ILEVEL, NM, NTYPE_GLOBAL)

      IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,*) 'SCARC-Multigrid, iteration ',ITE,': residuum=',SLHI%RES_MG

      !!! ------------------------------------------------------------------------
      !!! Send error information to Master  - still missing
      !!! ------------------------------------------------------------------------
      !CALL SCARC_INFO_MASTER 
      !CALL SCARC_INFO_LEVEL(ILEVEL)

      !!! ------------------------------------------------------------------------
      !!! Convergence or divergence ?
      !!! ------------------------------------------------------------------------
      IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) ITE, ILEVEL, SL%RES_MG
      IF (SCARC_DEBUG >= 1.AND.NM==1) write (*,1000) ITE, ILEVEL, SL%RES_MG
      IF (SCARC_BREL) THEN
         IF (SL%RES_MG <= SL%RESIN_MG*SCARC_MG_EPS) BCONV = .TRUE.
      ELSE
         IF (SL%RES_MG <= SCARC_MG_EPS .AND. SL%RES_MG <= SL%RESIN_MG*SCARC_EPS_REL) BCONV = .TRUE.
      ENDIF
      IF (SL%RES_MG > SCARC_EPS_DIVG) BDIVG = .TRUE.
      IF (BCONV.OR.BDIVG) EXIT MG_LOOP2D
 
   ENDDO MG_LOOP2D

   !!! ------------------------------------------------------------------------
   !!! Convergence or divergence ?
   !!! ------------------------------------------------------------------------
   IF (BDIVG) THEN                       
      ITE0 = - 1
      SL%CAPPA_MG = 1.0_EB
   ELSE 
      IF (BCONV) THEN
         ITE0=ITE
      ELSE
         ITE0=ITE-1
      ENDIF
      IF (SL%RESIN_MG >= 1.0E-70_EB) THEN
         SL%RESIN_MG = SL%RES_MG / SL%RESIN_MG
         SL%CAPPA_MG = SL%RESIN_MG ** (1.0_EB/ITE0)
      ELSE 
         SL%CAPPA_MG = 0.0_EB
      ENDIF
   ENDIF
   SCARC_NIT=ITE0
   SCARC_RES=SL%RES_MG
   SCARC_CAPPA=SL%CAPPA_MG
   IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,2000) ITE0, ILEVEL, SL%RES_MG, SL%CAPPA_MG 
   IF (SCARC_DEBUG >= 1.AND.NM==1) WRITE (*,2000) ITE0, ILEVEL, SL%RES_MG, SL%CAPPA_MG 

ENDIF ONLY_ONE_LEVEL_IF
 
! shift solution back to vector HP
DO K = 1, SLHI%KBAR                                 
   DO I = 1, SLHI%IBAR
      HP(I, 1, K) = SLHI%X (I, 1, K)
   ENDDO
ENDDO

! set ghost cell values along exterior boundaries 
CALL SCARC_GHOSTCELLS(HP,NM)

! get neighbouring data along interior boundaries
CALL SCARC_UPDATE_LEVEL(HP, NM, ILMAX, 'HPFIN ')

 
TUSED_SCARC(7,NM)=TUSED_SCARC(7,NM)+SECOND()-TNOW_MG2D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_MG2D

1000 FORMAT ('=====SCARC_MG       : #ite= ',i4,': level=',i4,': res=',e16.8)
2000 FORMAT ('=====SCARC_MG       : #ite= ',i4,': level=',i4,': res=',e16.8,': rate=',e16.8)
END SUBROUTINE SCARC_MG2D



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Perform global 2D cg-method based on global possion-matrix
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_COARSE2D (ILEVEL, NM)
 
INTEGER ::  NM, ILEVEL, ITE, ITE0
REAL (EB) :: SIGMA0, SIGMA1, ALPHA, GAMMA
REAL (EB) :: TNOW_COARSE
LOGICAL BCONV, BDIVG
TYPE (MESH_TYPE),  POINTER :: M
 
TNOW_COARSE = SECOND()
BCONV = .FALSE.
BDIVG = .FALSE.

S  => SCARC (NM)
M  => MESHES (NM)

! cg only works in finest grid level
SL    => S%SLEVEL(ILEVEL)
 
! re-initialize auxiliary vectors
SL%X = 0.0_EB                   
SL%G = 0.0_EB                   
SL%Y = 0.0_EB
SL%R = 0.0_EB
SL%D = 0.0_EB
 
IF (SCARC_DEBUG>2) WRITE(SCARC_LU,*) 'ILEVEL=',ILEVEL

!!! calculate initial defect Y = B - A*X and get l2-norm of it
CALL SCARC_MATVEC2D (SL%X, SL%R, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, NMV_X, NMV_R)
SL%R = -SL%F + SL%R
SL%RESIN_CO = SCARC_L2NORM2D (SL%R, ILEVEL, NM, NTYPE_GLOBAL)
IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) 0, ILEVEL, SL%RESIN_CO
IF (SCARC_DEBUG >= 1.AND.NM==1) write (*,1000) 0, ILEVEL, SL%RESIN_CO
 
! initial preconditioning
SELECT CASE(SCARC_CO_PRECON)
   CASE('NONE') 
      SIGMA0 = SL%RESIN_CO * SL%RESIN_CO
   CASE('JACOBI') 
      SL%G=SL%R
      CALL SCARC_PRECON_JACOBI2D (SL%G, ILEVEL, NM)
      SIGMA0 = SCARC_SCALPROD2D (SL%R, SL%G, ILEVEL, NM)
   CASE('SSOR') 
      SL%G=SL%R
      CALL SCARC_PRECON_SSOR2D (SL%G, ILEVEL, NM)
      SIGMA0 = SCARC_SCALPROD2D (SL%R, SL%G, ILEVEL, NM)
   CASE('GSTRIX') 
      SL%G=SL%R
      CALL SCARC_PRECON_GSTRIX2D (SL%G, ILEVEL, NM)
      SIGMA0 = SCARC_SCALPROD2D (SL%R, SL%G, ILEVEL, NM)
END SELECT
SL%D = -SL%G 
 
! start defect correction loop
!
CO_LOOP2D: DO ITE = 1, SCARC_CO_NIT
 
   ! calculate new defect and get L2-norm of it
   CALL SCARC_MATVEC2D (SL%D, SL%Y, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, NMV_D, NMV_Y)
   ALPHA = SCARC_SCALPROD2D (SL%D, SL%Y, ILEVEL, NM)


   ! get new descent direction
   ALPHA = SIGMA0 / ALPHA
   SL%X = ALPHA * SL%D + SL%X
   SL%R = ALPHA * SL%Y + SL%R
   SL%RES_CO = SCARC_L2NORM2D (SL%R, ILEVEL, NM, NTYPE_GLOBAL)

 
   !  exit loop in case of convergence or divergence
   IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) ITE, ILEVEL, SL%RES_CO
   IF (SCARC_DEBUG >= 1.AND.NM==1) write (*,1000) ITE, ILEVEL, SL%RES_CO
   IF (SCARC_BREL) THEN
      IF (SL%RES_CO <= SL%RESIN_CO*SCARC_CO_EPS) BCONV = .TRUE.
   ELSE
      IF (SL%RES_CO <= SCARC_CO_EPS .AND. SL%RES_CO <= SL%RESIN_CO*SCARC_EPS_REL) BCONV = .TRUE.
   ENDIF
   IF (SL%RES_CO > SCARC_EPS_DIVG) BDIVG = .TRUE.
   IF (BCONV.OR.BDIVG) EXIT CO_LOOP2D
 
   ! preconditioning
   SELECT CASE(SCARC_CO_PRECON)
      CASE('NONE')
         SIGMA1 = SL%RES_CO * SL%RES_CO
      CASE('JACOBI')
         SL%G=SL%R
         CALL SCARC_PRECON_JACOBI2D (SL%G, ILEVEL, NM)
         SIGMA1 = SCARC_SCALPROD2D (SL%R, SL%G, ILEVEL, NM)
      CASE('SSOR')
         SL%G=SL%R
         CALL SCARC_PRECON_SSOR2D (SL%G, ILEVEL, NM)
         SIGMA1 = SCARC_SCALPROD2D (SL%R, SL%G, ILEVEL, NM)
      CASE('GSTRIX')
         SL%G=SL%R
         CALL SCARC_PRECON_GSTRIX2D (SL%G, ILEVEL, NM)
         SIGMA1 = SCARC_SCALPROD2D (SL%R, SL%G, ILEVEL, NM)
   END SELECT
 
   GAMMA = SIGMA1 / SIGMA0
   SIGMA0 = SIGMA1
   SL%D = -SL%G + GAMMA * SL%D
 
ENDDO CO_LOOP2D
 
! divergence or convergence ?
IF (BDIVG) THEN                       
   ITE0 = - 1
   SL%CAPPA_CO = 1.0_EB
ELSE 
   IF (BCONV) THEN
      ITE0=ITE
   ELSE
      ITE0=ITE-1
   ENDIF
   IF (SL%RESIN_CO >= 1.0E-70_EB) THEN
      SL%RESIN_CO = SL%RES_CO / SL%RESIN_CO
      SL%CAPPA_CO = SL%RESIN_CO ** (1.0_EB/ITE0)
   ELSE 
      SL%CAPPA_CO = 0.0_EB
   ENDIF
ENDIF
IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,2000) ITE0, ILEVEL, SL%RES_CO, SL%CAPPA_CO 
IF (SCARC_DEBUG >= 1.AND.NM==1) WRITE (*,2000) ITE0, ILEVEL, SL%RES_CO, SL%CAPPA_CO 

 
TUSED_SCARC(9,NM)=TUSED_SCARC(9,NM)+SECOND()-TNOW_COARSE
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_COARSE

1000 FORMAT ('     SCARC_COARSE   : #ite= ',i4,': level=',i4,': res=',e16.8)
2000 FORMAT ('     SCARC_COARSE   : #ite= ',i4,': level=',i4,': res=',e16.8,': rate=',e16.8)
END SUBROUTINE SCARC_COARSE2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Perform global 2D cg-method based on global possion-matrix
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PRECON_COARSE2D (ILEVEL, NM)
 
INTEGER ::  NM, ILEVEL, ITE, ITE0
REAL (EB) :: SIGMA0, SIGMA1, ALPHA, GAMMA
REAL (EB) :: TNOW_COARSE
LOGICAL BCONV, BDIVG
TYPE (MESH_TYPE),  POINTER :: M
 
TNOW_COARSE = SECOND()
BCONV = .FALSE.
BDIVG = .FALSE.

S  => SCARC (NM)
M  => MESHES (NM)

! cg only works in finest grid level
SL    => S%SLEVEL(ILEVEL)
 
! re-initialize auxiliary vectors
SL%X2 = 0.0_EB                   
SL%G2 = 0.0_EB                   
SL%Y2 = 0.0_EB
SL%R2 = 0.0_EB
SL%D2 = 0.0_EB
 
IF (SCARC_DEBUG>2) WRITE(SCARC_LU,*) 'ILEVEL=',ILEVEL

!!! calculate initial defect Y = B - A*X and get l2-norm of it
CALL SCARC_MATVEC2D (SL%X2, SL%R2, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, NMV_X2, NMV_R2)
SL%R2 = -SL%F2 + SL%R2
SL%RESIN_CO = SCARC_L2NORM2D (SL%R2, ILEVEL, NM, NTYPE_GLOBAL)
IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) 0, ILEVEL, SL%RESIN_CO
IF (SCARC_DEBUG >= 1.AND.NM==1) write (*,1000) 0, ILEVEL, SL%RESIN_CO
 
! initial preconditioning
SELECT CASE(SCARC_CO_PRECON)
   CASE('NONE') 
      SIGMA0 = SL%RESIN_CO * SL%RESIN_CO
   CASE('JACOBI') 
      SL%G2=SL%R2
      CALL SCARC_PRECON_JACOBI2D (SL%G2, ILEVEL, NM)
      SIGMA0 = SCARC_SCALPROD2D (SL%R2, SL%G2, ILEVEL, NM)
   CASE('SSOR') 
      SL%G2=SL%R2
      CALL SCARC_PRECON_SSOR2D (SL%G2, ILEVEL, NM)
      SIGMA0 = SCARC_SCALPROD2D (SL%R2, SL%G2, ILEVEL, NM)
   CASE('GSTRIX') 
      SL%G2=SL%R2
      CALL SCARC_PRECON_GSTRIX2D (SL%G2, ILEVEL, NM)
      SIGMA0 = SCARC_SCALPROD2D (SL%R2, SL%G2, ILEVEL, NM)
END SELECT
SL%D2 = -SL%G2 
 
!
! start defect correction loop
!
CO_LOOP2D: DO ITE = 1, SCARC_CO_NIT
 
   ! calculate new defect and get L2-norm of it
   CALL SCARC_MATVEC2D (SL%D2, SL%Y2, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, NMV_D2, NMV_Y2)
   ALPHA = SCARC_SCALPROD2D (SL%D2, SL%Y2, ILEVEL, NM)

   ! get new descent direction
   ALPHA = SIGMA0 / ALPHA
   SL%X2 = ALPHA * SL%D2 + SL%X2
   SL%R2 = ALPHA * SL%Y2 + SL%R2
   SL%RES_CO = SCARC_L2NORM2D (SL%R2, ILEVEL, NM, NTYPE_GLOBAL)

 
   !  exit loop in case of convergence or divergence
   IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) ITE, ILEVEL, SL%RES_CO
   IF (SCARC_DEBUG >= 1.AND.NM==1) write (*,1000) ITE, ILEVEL, SL%RES_CO
   IF (SCARC_BREL) THEN
      IF (SL%RES_CO <= SL%RESIN_CO*SCARC_CO_EPS) BCONV = .TRUE.
   ELSE
      IF (SL%RES_CO <= SCARC_CO_EPS .AND. SL%RES_CO <= SL%RESIN_CO*SCARC_EPS_REL) BCONV = .TRUE.
   ENDIF
   IF (SL%RES_CO > SCARC_EPS_DIVG) BDIVG = .TRUE.
   IF (BCONV.OR.BDIVG) EXIT CO_LOOP2D
 
   ! preconditioning
   SELECT CASE(SCARC_CO_PRECON)
      CASE('NONE')
         SIGMA1 = SL%RES_CO * SL%RES_CO
      CASE('JACOBI')
         SL%G2=SL%R2
         CALL SCARC_PRECON_JACOBI2D (SL%G2, ILEVEL, NM)
         SIGMA1 = SCARC_SCALPROD2D (SL%R2, SL%G, ILEVEL, NM)
      CASE('SSOR')
         SL%G2=SL%R2
         CALL SCARC_PRECON_SSOR2D (SL%G2, ILEVEL, NM)
         SIGMA1 = SCARC_SCALPROD2D (SL%R2, SL%G2, ILEVEL, NM)
      CASE('GSTRIX')
         SL%G2=SL%R2
         CALL SCARC_PRECON_GSTRIX2D (SL%G2, ILEVEL, NM)
         SIGMA1 = SCARC_SCALPROD2D (SL%R2, SL%G2, ILEVEL, NM)
   END SELECT
 
   GAMMA = SIGMA1 / SIGMA0
   SIGMA0 = SIGMA1
   SL%D2 = -SL%G2 + GAMMA * SL%D2
 
ENDDO CO_LOOP2D
 
! divergence or convergence ?
IF (BDIVG) THEN                       
   ITE0 = - 1
   SL%CAPPA_CO = 1.0_EB
ELSE 
   IF (BCONV) THEN
      ITE0=ITE
   ELSE
      ITE0=ITE-1
   ENDIF
   IF (SL%RESIN_CO >= 1.0E-70_EB) THEN
      SL%RESIN_CO = SL%RES_CO / SL%RESIN_CO
      SL%CAPPA_CO = SL%RESIN_CO ** (1.0_EB/ITE0)
   ELSE 
      SL%CAPPA_CO = 0.0_EB
   ENDIF
ENDIF
IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,2000) ITE0, ILEVEL, SL%RES_CO, SL%CAPPA_CO 
IF (SCARC_DEBUG >= 1.AND.NM==1) WRITE (*,2000) ITE0, ILEVEL, SL%RES_CO, SL%CAPPA_CO 

 
TUSED_SCARC(9,NM)=TUSED_SCARC(9,NM)+SECOND()-TNOW_COARSE
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_COARSE

1000 FORMAT ('     SCARC_PRECON_COARSE2D   : #ite= ',i4,': level=',i4,': res=',e16.8)
2000 FORMAT ('     SCARC_PRECON_COARSE2D   : #ite= ',i4,': level=',i4,': res=',e16.8,': rate=',e16.8)
END SUBROUTINE SCARC_PRECON_COARSE2D


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Perform smoothing in SCARC_MG2D
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SMOOTHER2D(X, F, D, BMATVEC, ILEVEL, NM)

REAL (EB), POINTER, DIMENSION (:, :, :) :: X, F, D
INTEGER :: NM
INTEGER :: ILEVEL, ITE, ITE0
REAL(EB):: TNOW_SMOOTHER2D
LOGICAL :: BMATVEC, BL2NORM=.TRUE.       ! necessary ???
LOGICAL :: BCONV, BDIVG


TNOW_SMOOTHER2D = SECOND()

SL => S%SLEVEL(ILEVEL)
BCONV=.FALSE.
BDIVG=.FALSE.

IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,*) 'BMATVEC=',BMATVEC

SL%RESIN_SM=1.0_EB
IF (BMATVEC) THEN
   CALL SCARC_MATVEC2D (X, D, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, NMV_X, NMV_D)   ! global MATVEC
   D = F - D
   IF (BL2NORM) SL%RESIN_SM = SCARC_L2NORM2D (D, ILEVEL, NM, NTYPE_GLOBAL)
ENDIF

IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) 0, ILEVEL, SL%RESIN_SM
IF (SCARC_DEBUG >= 1.AND.NM==1) WRITE(*,1000) 0, ILEVEL, SL%RESIN_SM

SMOOTH_LOOP2D: DO ITE=1,SCARC_SM_NIT
 
   SELECT CASE(SCARC_SM_PRECON)
      CASE('JACOBI') 
         CALL SCARC_PRECON_JACOBI2D (D, ILEVEL, NM)
      CASE('SSOR') 
         CALL SCARC_PRECON_SSOR2D (D, ILEVEL, NM)
      CASE('GSTRIX') 
         CALL SCARC_PRECON_GSTRIX2D (D, ILEVEL, NM)
      CASE('FFT') 
         CALL SCARC_PRECON_FFT2D (D, D, ILEVEL, NM)
   END SELECT

   X = SCARC_SM_OMEGA * D + X
   CALL SCARC_MATVEC2D (X, D, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, NMV_X, NMV_D)   ! global MATVEC
   D = F - D
   IF (BL2NORM) SL%RES_SM = SCARC_L2NORM2D (D, ILEVEL, NM, NTYPE_GLOBAL)

   IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) ITE, ILEVEL, SL%RES_SM
   IF (SCARC_DEBUG >= 1.AND.NM==1) WRITE(*,1000) ITE, ILEVEL, SL%RES_SM
   IF (SCARC_BREL) THEN
      IF (SL%RES_SM <= SL%RESIN_SM*SCARC_SM_EPS) BCONV = .TRUE.
   ELSE
      IF (SL%RES_SM <= SCARC_SM_EPS .AND. SL%RES_SM <= SL%RESIN_SM*SCARC_EPS_REL) BCONV = .TRUE.
   ENDIF
   IF (SL%RES_SM > SCARC_EPS_DIVG) BDIVG = .TRUE.
   IF (BCONV.OR.BDIVG) EXIT SMOOTH_LOOP2D
 

ENDDO SMOOTH_LOOP2D

! divergence or convergence ?
IF (BDIVG) THEN                       
   ITE0 = - 1
   SL%CAPPA_SM = 1.0_EB
ELSE 
   IF (BCONV) THEN
     ITE0=ITE
   ELSE
     ITE0=ITE-1
   ENDIF
   IF (SL%RESIN_SM >= 1.0E-70_EB) THEN
      SL%RESIN_SM = SL%RES_SM / SL%RESIN_SM
      SL%CAPPA_SM = SL%RESIN_SM ** (1.0_EB/ITE0)
   ELSE 
      SL%CAPPA_SM = 0.0_EB
   ENDIF
ENDIF
IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,2000) ITE0, ILEVEL, SL%RES_SM, SL%CAPPA_SM 
IF (SCARC_DEBUG >= 1.AND.NM==1) WRITE (*,2000) ITE0, ILEVEL, SL%RES_SM, SL%CAPPA_SM 

TUSED_SCARC(8,NM)=TUSED_SCARC(8,NM)+SECOND()-TNOW_SMOOTHER2D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_SMOOTHER2D

1000 FORMAT ('     SCARC_SMOOTHER : #ite= ',i4,': level=',i4,': res=',e16.8)
2000 FORMAT ('     SCARC_SMOOTHER : #ite= ',i4,': level=',i4,': res=',e16.8,': rate=',e16.8)
END SUBROUTINE SCARC_SMOOTHER2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Perform smoothing in SCARC_MG2D
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PRECON_SMOOTHER2D(X, F, D, BMATVEC, ILEVEL, NM)

REAL (EB), POINTER, DIMENSION (:, :, :) :: X, F, D
INTEGER :: NM
INTEGER :: ILEVEL, ITE, ITE0
REAL(EB):: TNOW_PRECON_SMOOTHER2D
LOGICAL :: BMATVEC, BL2NORM=.TRUE.       ! necessary ???
LOGICAL :: BCONV, BDIVG


TNOW_PRECON_SMOOTHER2D = SECOND()

SL => S%SLEVEL(ILEVEL)
BCONV=.FALSE.
BDIVG=.FALSE.

IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,*) 'BMATVEC=',BMATVEC

SL%RESIN_SM=1.0_EB
IF (BMATVEC) THEN
   CALL SCARC_MATVEC2D (X, D, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, NMV_X2, NMV_D2)   ! global MATVEC
   D = F - D
   IF (BL2NORM) SL%RESIN_SM = SCARC_L2NORM2D (D, ILEVEL, NM, NTYPE_GLOBAL)
ENDIF

IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) 0, ILEVEL, SL%RESIN_SM
IF (SCARC_DEBUG >= 1.AND.NM==1) WRITE(*,1000) 0, ILEVEL, SL%RESIN_SM

SMOOTH_LOOP2D: DO ITE=1,SCARC_SM_NIT
 
   SELECT CASE(SCARC_SM_PRECON)
      CASE('JACOBI') 
         CALL SCARC_PRECON_JACOBI2D (D, ILEVEL, NM)
      CASE('SSOR') 
         CALL SCARC_PRECON_SSOR2D (D, ILEVEL, NM)
      CASE('GSTRIX') 
         CALL SCARC_PRECON_GSTRIX2D (D, ILEVEL, NM)
      CASE('FFT') 
         CALL SCARC_PRECON_FFT2D (D, D, ILEVEL, NM)
   END SELECT

   X = SCARC_SM_OMEGA * D + X
   CALL SCARC_MATVEC2D (X, D, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, NMV_X2, NMV_D2)   ! global MATVEC
   D = F - D
   IF (BL2NORM) SL%RES_SM = SCARC_L2NORM2D (D, ILEVEL, NM, NTYPE_GLOBAL)

   IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) ITE, ILEVEL, SL%RES_SM
   IF (SCARC_DEBUG >= 1.AND.NM==1) WRITE(*,1000) ITE, ILEVEL, SL%RES_SM
   IF (SCARC_BREL) THEN
      IF (SL%RES_SM <= SL%RESIN_SM*SCARC_SM_EPS) BCONV = .TRUE.
   ELSE
      IF (SL%RES_SM <= SCARC_SM_EPS .AND. SL%RES_SM <= SL%RESIN_SM*SCARC_EPS_REL) BCONV = .TRUE.
   ENDIF
   IF (SL%RES_SM > SCARC_EPS_DIVG) BDIVG = .TRUE.
   IF (BCONV.OR.BDIVG) EXIT SMOOTH_LOOP2D
 

ENDDO SMOOTH_LOOP2D

! divergence or convergence ?
IF (BDIVG) THEN                       
   ITE0 = - 1
   SL%CAPPA_SM = 1.0_EB
ELSE 
   IF (BCONV) THEN
     ITE0=ITE
   ELSE
     ITE0=ITE-1
   ENDIF
   IF (SL%RESIN_SM >= 1.0E-70_EB) THEN
      SL%RESIN_SM = SL%RES_SM / SL%RESIN_SM
      SL%CAPPA_SM = SL%RESIN_SM ** (1.0_EB/ITE0)
   ELSE 
      SL%CAPPA_SM = 0.0_EB
   ENDIF
ENDIF
IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,2000) ITE0, ILEVEL, SL%RES_SM, SL%CAPPA_SM 
IF (SCARC_DEBUG >= 1.AND.NM==1) WRITE (*,2000) ITE0, ILEVEL, SL%RES_SM, SL%CAPPA_SM 

TUSED_SCARC(8,NM)=TUSED_SCARC(8,NM)+SECOND()-TNOW_PRECON_SMOOTHER2D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_PRECON_SMOOTHER2D

1000 FORMAT ('     SCARC_PRECON_SMOOTHER2D : #ite= ',i4,': level=',i4,': res=',e16.8)
2000 FORMAT ('     SCARC_PRECON_SMOOTHER2D : #ite= ',i4,': level=',i4,': res=',e16.8,': rate=',e16.8)
END SUBROUTINE SCARC_PRECON_SMOOTHER2D

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Perform global 3D cg-method based on global possion-matrix
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_CG3D (NM, HP, PRHS)
 
INTEGER ::  NM, ILMAX, ITE, ITE0, I, J, K
REAL (EB), POINTER, DIMENSION (:, :, :) :: HP, PRHS
REAL (EB) :: SIGMA0, SIGMA1, ALPHA, GAMMA
REAL (EB) :: TNOW_CG3D
LOGICAL BCONV, BDIVG
TYPE (MESH_TYPE),  POINTER :: M
 
TNOW_CG3D = SECOND()
BCONV=.FALSE.
BDIVG=.FALSE.


S  => SCARC (NM)
M  => MESHES (NM)

! cg only works in finest grid level
ILMAX =  S%NLMAX                              
SL    => S%SLEVEL(ILMAX)
 
! get right hand side and initial vector
DO K = 1, SL%KBAR                                     ! 
   DO J = 1, SL%JBAR                                     ! 
      DO I = 1, SL%IBAR
         SL%F (I, J, K) = PRHS (I, J, K)
         SL%X (I, J, K) = HP (I, J, K)
      ENDDO
   ENDDO
ENDDO

! re-initialize auxiliary vectors
SL%G = 0.0_EB                   
SL%Y = 0.0_EB
SL%R = 0.0_EB
SL%D = 0.0_EB
 
IF (SCARC_DEBUG>=6) THEN
   CALL SCARC_SHOW0(SL%X,'SCARC_CG3D','XCG before all')
   CALL SCARC_SHOW0(SL%F,'SCARC_CG3D','FCG before all')
ENDIF

! set boundary conditions along exterior boundaries
CALL SCARC_SETBDRY3D (ILMAX, NM)
 
!!! calculate initial defect Y = B - A*X and get l2-norm of it
CALL SCARC_MATVEC3D (SL%X, SL%R, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILMAX, NMV_X, NMV_R)
SL%R = -SL%F + SL%R
SL%RESIN_CG = SCARC_L2NORM3D (SL%R, ILMAX, NM, NTYPE_GLOBAL)

IF (SCARC_DEBUG>=6) CALL SCARC_SHOW0(SL%R,'SCARC_CG3D','RCG before all')

IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) 0, SL%RESIN_CG
IF (SCARC_DEBUG >= 1.AND.NM==1) write (*,1000) 0, SL%RESIN_CG
 
IF (SCARC_DEBUG>=6) THEN
   WRITE(SCARC_LU,*) '============ STARTING SCARC_CG'
   WRITE(SCARC_LU,*) 'SCARC_CG_PRECON=',SCARC_CG_PRECON
ENDIF

! initial preconditioning
SELECT CASE(SCARC_CG_PRECON)
   CASE('NONE') 
      SIGMA0 = SL%RESIN_CG * SL%RESIN_CG
   CASE('JACOBI') 
      SL%G=SL%R
      CALL SCARC_PRECON_JACOBI3D (SL%G, ILMAX, NM)
      SIGMA0 = SCARC_SCALPROD3D (SL%R, SL%G, ILMAX, NM)
   CASE('SSOR') 
      SL%G=SL%R
      CALL SCARC_PRECON_SSOR3D (SL%G, ILMAX, NM)
      SIGMA0 = SCARC_SCALPROD3D (SL%R, SL%G, ILMAX, NM)
   CASE('GSTRIX') 
      SL%G=SL%R
      CALL SCARC_PRECON_GSTRIX3D (SL%G, ILMAX, NM)
      SIGMA0 = SCARC_SCALPROD3D (SL%R, SL%G, ILMAX, NM)
   CASE('MG') 
      SL%G=SL%R
      CALL SCARC_PRECON_MG3D (SL%G, ILMAX, NM)
      SIGMA0 = SCARC_SCALPROD3D (SL%R, SL%G, ILMAX, NM)
      CASE('FFT')
         CALL SCARC_PRECON_FFT3D (SL%R, SL%G, ILMAX, NM)
         SIGMA0 = SCARC_SCALPROD3D (SL%R, SL%G, ILMAX, NM)
END SELECT
SL%D = -SL%G
IF (SCARC_DEBUG>=6) CALL SCARC_SHOW0(SL%D,'SCARC_CG3D','DCG before all')
 
!
! defect correction loop
!
CG_LOOP3D: DO ITE = 1, SCARC_CG_NIT
 
   ! calculate new defect and get L2-norm of it
   CALL SCARC_MATVEC3D (SL%D, SL%Y, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILMAX, NMV_D, NMV_Y)
   ALPHA = SCARC_SCALPROD3D (SL%D, SL%Y, ILMAX, NM)

IF (SCARC_DEBUG>=6) CALL SCARC_SHOW0(SL%Y,'SCARC_CG3D','YITE in  all')

IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,*) 'Alpha=',ALPHA
   ! get new descent direction
   ALPHA = SIGMA0 / ALPHA
   SL%X = ALPHA * SL%D + SL%X
   SL%R = ALPHA * SL%Y + SL%R
   SL%RES_CG = SCARC_L2NORM3D (SL%R, ILMAX, NM, NTYPE_GLOBAL)
IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,*) 'Res=',SL%RES_CG
 
   !  exit loop in case of convergence or divergence
   IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) ITE, SL%RES_CG
   IF (SCARC_DEBUG >= 1.AND.NM==1) write (*,1000) ITE, SL%RES_CG
   IF (VELOCITY_ERROR_FILE) WRITE(LU_VELOCITY_ERROR,'(2(I5,A),E16.8)') ICYC,', ',ITE,', ',SL%RES_CG
   IF (SCARC_BREL) THEN
      IF (SL%RES_CG <= SL%RESIN_CG*SCARC_CG_EPS) BCONV = .TRUE.
   ELSE
      IF (SL%RES_CG <= SCARC_CG_EPS .AND. SL%RES_CG <= SL%RESIN_CG*SCARC_EPS_REL) BCONV = .TRUE.
   ENDIF
   IF (SL%RES_CG > SCARC_EPS_DIVG) BDIVG = .TRUE.
   IF (BCONV.OR.BDIVG) EXIT CG_LOOP3D
 
   ! preconditioning
   SELECT CASE(SCARC_CG_PRECON)
      CASE('NONE')
         SIGMA1 = SL%RES_CG * SL%RES_CG
      CASE('JACOBI')
         SL%G = SL%R
         CALL SCARC_PRECON_JACOBI3D (SL%G, ILMAX, NM)
         SIGMA1 = SCARC_SCALPROD3D (SL%R, SL%G, ILMAX, NM)
      CASE('SSOR')
         SL%G = SL%R
         CALL SCARC_PRECON_SSOR3D (SL%G, ILMAX, NM)
         SIGMA1 = SCARC_SCALPROD3D (SL%R, SL%G, ILMAX, NM)
      CASE('GSTRIX')
         SL%G = SL%R
         CALL SCARC_PRECON_GSTRIX3D (SL%G, ILMAX, NM)
         SIGMA1 = SCARC_SCALPROD3D (SL%R, SL%G, ILMAX, NM)
      CASE('MG')
         SL%G = SL%R
         CALL SCARC_PRECON_MG3D (SL%G, ILMAX, NM)
         SIGMA1 = SCARC_SCALPROD3D (SL%R, SL%G, ILMAX, NM)
      CASE('FFT')
         CALL SCARC_PRECON_FFT3D (SL%R, SL%G, ILMAX, NM)
         SIGMA1 = SCARC_SCALPROD3D (SL%R, SL%G, ILMAX, NM)
   END SELECT
 
   GAMMA = SIGMA1 / SIGMA0
   SIGMA0 = SIGMA1
   SL%D = -SL%G + GAMMA * SL%D
 
ENDDO CG_LOOP3D
 
! divergence or convergence ?
IF (BDIVG) THEN                       
   ITE0 = - 1
   SL%CAPPA_CG = 1.0_EB
ELSE 
   IF (BCONV) THEN
      ITE0=ITE
   ELSE
      ITE0=ITE-1
   ENDIF
   IF (SL%RESIN_CG >= 1.0E-70_EB) THEN
      SL%RESIN_CG = SL%RES_CG / SL%RESIN_CG
      SL%CAPPA_CG = SL%RESIN_CG ** (1.0_EB/ITE0)
   ELSE 
      SL%CAPPA_CG = 0.0_EB
   ENDIF
ENDIF
SCARC_NIT=ITE0
SCARC_RES=SL%RES_CG
SCARC_CAPPA=SL%CAPPA_CG
IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,2000) ITE0, SL%RES_CG, SL%CAPPA_CG 
IF (SCARC_DEBUG >= 1.AND.NM==1) WRITE(*,2000) ITE0, SL%RES_CG, SL%CAPPA_CG 

! shift solution back to vector HP
DO K = 1, SL%KBAR                                 
   DO J = 1, SL%JBAR                                 
      DO I = 1, SL%IBAR
         HP(I, J, K) = SL%X (I, J, K)
      ENDDO
   ENDDO
ENDDO

CALL SCARC_SHOW0(HP,'SCARC_CG3D','HP before GHOSTCELLS')
! set ghost cell values along exterior boundaries 
CALL SCARC_GHOSTCELLS(HP,NM)

CALL SCARC_SHOW0(HP,'SCARC_CG3D','HP after GHOSTCELLS')
! get neighbouring data along interior boundaries
CALL SCARC_UPDATE_LEVEL(HP, NM, ILMAX, 'HPFIN ')
CALL SCARC_SHOW0(HP,'SCARC_CG3D','HP after UPDATE')


TUSED_SCARC(6,NM)=TUSED_SCARC(6,NM)+SECOND()-TNOW_CG3D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_CG3D

1000 FORMAT ('=====SCARC_CG       : #ite= ',i4,': res=',e16.8)
2000 FORMAT ('=====SCARC_CG       : #ite= ',i4,': res=',e16.8,': rate=',e16.8)
END SUBROUTINE SCARC_CG3D


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Perform global 3D cg-method based on global possion-matrix
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_BICG3D (NM, HP, PRHS)
 
INTEGER ::  NM, ILMAX, ITE, ITE0, I, J, K
REAL (EB), POINTER, DIMENSION (:, :, :) :: HP, PRHS
REAL (EB) :: RHO0, RHO1, DBETA, DTHETA, ALPHA0, ALPHA1, ALPHA2
REAL (EB) :: TNOW_BICG3D
LOGICAL BCONV, BDIVG
TYPE (MESH_TYPE),  POINTER :: M
 
TNOW_BICG3D = SECOND()
BCONV=.FALSE.
BDIVG=.FALSE.

S  => SCARC (NM)
M  => MESHES (NM)

IF (SCARC_DEBUG >=3) WRITE(SCARC_LU,*) 'STARTING BICG3D'

! cg only works in finest grid level
ILMAX =  S%NLMAX                              
SL    => S%SLEVEL(ILMAX)
 
! get right hand side and initial vector
DO K = 1, SL%KBAR                                     ! 
   DO J = 1, SL%JBAR                                     ! 
      DO I = 1, SL%IBAR
         SL%F (I, J, K) = PRHS (I, J, K)
         SL%X (I, J, K) = HP (I, J, K)
      ENDDO
   ENDDO
ENDDO

! preset some variables
RHO0  =1.0_EB
RHO1  =0.0_EB
DBETA =0.0_EB
DTHETA=1.0_EB
ALPHA0=1.0_EB

! re-initialize auxiliary vectors
SL%R = 0.0_EB
SL%G = 0.0_EB                   
SL%Z = 0.0_EB
SL%Y = 0.0_EB
SL%D = 0.0_EB
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! initialization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CALL SCARC_SETBDRY3D (ILMAX, NM)
SL%R=SL%F

SELECT CASE(SCARC_BICG_PRECON)
   CASE('JACOBI') 
      CALL SCARC_PRECON_JACOBI3D (SL%R, ILMAX, NM)
   CASE('SSOR') 
      CALL SCARC_PRECON_SSOR3D (SL%R, ILMAX, NM)
   CASE('GSTRIX') 
      CALL SCARC_PRECON_GSTRIX3D (SL%R, ILMAX, NM)
   CASE('MG') 
      CALL SCARC_PRECON_MG3D (SL%R, ILMAX, NM)
   CASE('FFT') 
      CALL SCARC_PRECON_FFT3D (SL%R, SL%R, ILMAX, NM)
END SELECT

CALL SCARC_MATVEC3D (SL%X, SL%R, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILMAX, NMV_X, NMV_R)
SL%R = SL%F - SL%R

SELECT CASE(SCARC_BICG_PRECON)
   CASE('JACOBI') 
      CALL SCARC_PRECON_JACOBI3D (SL%R, ILMAX, NM)
   CASE('SSOR') 
      CALL SCARC_PRECON_SSOR3D (SL%R, ILMAX, NM)
   CASE('GSTRIX') 
      CALL SCARC_PRECON_GSTRIX3D (SL%R, ILMAX, NM)
   CASE('MG') 
      CALL SCARC_PRECON_MG3D (SL%R, ILMAX, NM)
   CASE('FFT') 
      CALL SCARC_PRECON_FFT3D (SL%R, SL%R, ILMAX, NM)
END SELECT

SL%RESIN_BICG = SCARC_L2NORM3D (SL%R, ILMAX, NM, NTYPE_GLOBAL)
SL%G=SL%R

IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) 0, SL%RESIN_BICG
IF (SCARC_DEBUG >= 1.AND.NM==1) write (*,1000) 0, SL%RESIN_BICG

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! iterative correction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
BICG_LOOP3D: DO ITE = 1, SCARC_BICG_NIT
 
   RHO1 = SCARC_SCALPROD3D (SL%G, SL%R, ILMAX, NM)
   DBETA=(RHO1*DTHETA)/(RHO0*ALPHA0)

   RHO0=RHO1

   SL%Z = SL%R + DBETA * SL%Z
   SL%Z = -DBETA * ALPHA0 * SL%Y + SL%Z

   ! matrix-vector multiplication and preconditioning
   CALL SCARC_MATVEC3D (SL%Z, SL%Y, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILMAX, NMV_Z, NMV_Y)

   SELECT CASE(SCARC_BICG_PRECON)
      CASE('JACOBI') 
         CALL SCARC_PRECON_JACOBI3D (SL%Y, ILMAX, NM)
      CASE('SSOR') 
         CALL SCARC_PRECON_SSOR3D (SL%Y, ILMAX, NM)
      CASE('GSTRIX') 
         CALL SCARC_PRECON_GSTRIX3D (SL%Y, ILMAX, NM)
      CASE('MG') 
         CALL SCARC_PRECON_MG3D (SL%Y, ILMAX, NM)
      CASE('FFT') 
         CALL SCARC_PRECON_FFT3D (SL%Y, SL%Y, ILMAX, NM)
   END SELECT

   DTHETA = SCARC_SCALPROD3D (SL%G, SL%Y, ILMAX, NM)
   DTHETA=RHO1/DTHETA
   SL%R = -DTHETA * SL%Y + SL%R

   ! matrix-vector multiplication and preconditioning
   CALL SCARC_MATVEC3D (SL%R, SL%D, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILMAX, NMV_R, NMV_D)

   SELECT CASE(SCARC_BICG_PRECON)
      CASE('JACOBI') 
         CALL SCARC_PRECON_JACOBI3D (SL%D, ILMAX, NM)
      CASE('SSOR') 
         CALL SCARC_PRECON_SSOR3D (SL%D, ILMAX, NM)
      CASE('GSTRIX') 
         CALL SCARC_PRECON_GSTRIX3D (SL%D, ILMAX, NM)
      CASE('MG') 
         CALL SCARC_PRECON_MG3D (SL%D, ILMAX, NM)
      CASE('FFT') 
         CALL SCARC_PRECON_FFT3D (SL%D, SL%D, ILMAX, NM)
   END SELECT

   ALPHA1 = SCARC_SCALPROD3D (SL%D, SL%R, ILMAX, NM)
   ALPHA2 = SCARC_SCALPROD3D (SL%D, SL%D, ILMAX, NM)

   ALPHA0=ALPHA1/ALPHA2

   SL%X =  DTHETA * SL%Z + SL%X
   SL%X =  ALPHA0 * SL%R + SL%X
   SL%R = -ALPHA0 * SL%D + SL%R

   SL%RES_BICG = SCARC_L2NORM3D (SL%R, ILMAX, NM, NTYPE_GLOBAL)
 
   !  exit loop in case of convergence or divergence
   IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) ITE, SL%RES_BICG
   IF (SCARC_DEBUG >= 1.AND.NM==1) write (*,1000) ITE, SL%RES_BICG
   IF (VELOCITY_ERROR_FILE) WRITE(LU_VELOCITY_ERROR,'(2(I5,A),E16.8)') ICYC,', ',ITE,', ',SL%RES_BICG
   IF (SCARC_BREL) THEN
      IF (SL%RES_BICG <= SL%RESIN_BICG*SCARC_BICG_EPS) BCONV = .TRUE.
   ELSE
      IF (SL%RES_BICG <= SCARC_BICG_EPS .AND. SL%RES_BICG <= SL%RESIN_BICG*SCARC_EPS_REL) BCONV = .TRUE.
   ENDIF
   IF (SL%RES_BICG > SCARC_EPS_DIVG) BDIVG = .TRUE.
   IF (BCONV.OR.BDIVG) EXIT BICG_LOOP3D
 
ENDDO BICG_LOOP3D
 
! divergence or convergence ?
IF (BDIVG) THEN                       
   ITE0 = - 1
   SL%CAPPA_BICG = 1.0_EB
ELSE 
   IF (BCONV) THEN
      ITE0=ITE
   ELSE
      ITE0=ITE-1
   ENDIF
   IF (SL%RESIN_BICG >= 1.0E-70_EB) THEN
      SL%RESIN_BICG = SL%RES_BICG / SL%RESIN_BICG
      SL%CAPPA_BICG = SL%RESIN_BICG ** (1.0_EB/ITE0)
   ELSE 
      SL%CAPPA_BICG = 0.0_EB
   ENDIF
ENDIF
SCARC_NIT=ITE0
SCARC_RES=SL%RES_BICG
SCARC_CAPPA=SL%CAPPA_BICG
IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,2000) ITE0, SL%RES_BICG, SL%CAPPA_BICG 
IF (SCARC_DEBUG >= 1.AND.NM==1) WRITE(*,2000) ITE0, SL%RES_BICG, SL%CAPPA_BICG 

! shift solution back to vector HP
DO K = 1, SL%KBAR                                 
   DO J = 1, SL%JBAR                                 
      DO I = 1, SL%IBAR
         HP(I, J, K) = SL%X (I, J, K)
      ENDDO
   ENDDO
ENDDO

! set ghost cell values along exterior boundaries 
CALL SCARC_GHOSTCELLS(HP,NM)

! get neighbouring data along interior boundaries
CALL SCARC_UPDATE_LEVEL(HP, NM, ILMAX, 'HPFIN ')

TUSED_SCARC(6,NM)=TUSED_SCARC(6,NM)+SECOND()-TNOW_BICG3D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_BICG3D

1000 FORMAT ('=====SCARC_BICG     : #ite= ',i4,': res=',e16.8)
2000 FORMAT ('=====SCARC_BICG     : #ite= ',i4,': res=',e16.8,': rate=',e16.8)
END SUBROUTINE SCARC_BICG3D



 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Perform global MG-method based on global possion-matrix
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_MG3D (NM, HP, PRHS)
 
INTEGER :: NM, ITE, ITE0, ICYCLE
REAL (EB), POINTER, DIMENSION (:, :, :) :: HP, PRHS
INTEGER :: ILEVEL, ILMAX
INTEGER :: I, J, K
REAL (EB) :: TNOW_MG3D
LOGICAL BMATVEC, BCONV, BDIVG
TYPE (MESH_TYPE),  POINTER :: M
 
TNOW_MG3D = SECOND()

! initialize working vectors
S    => SCARC (NM)
M    => MESHES (NM)

ILMAX =  S%NLMAX
SLHI  => S%SLEVEL(ILMAX)
 
! get right hand side and initial vector
DO K = 1, SLHI%KBAR
   DO J = 1, SLHI%JBAR
      DO I = 1, SLHI%IBAR
         SLHI%F (I, J, K) = PRHS (I, J, K)
         SLHI%X (I, J, K) = HP (I, J, K)
      ENDDO
   ENDDO
ENDDO
 
IF (SCARC_DEBUG >= 2) THEN
   WRITE(SCARC_LU,*) 'Starting SCARC_MG3D'
   WRITE(SCARC_LU,*) 'IBAR=:', SLHI%IBAR
   WRITE(SCARC_LU,*) 'JBAR=:', SLHI%JBAR
   WRITE(SCARC_LU,*) 'KBAR=:', SLHI%KBAR
   WRITE(SCARC_LU,*) 'SCARC_MG_NIT=:', SCARC_MG_NIT
   WRITE(SCARC_LU,*) 'SCARC_MG_EPS=:', SCARC_MG_EPS
   WRITE(SCARC_LU,*) 'SCARC_EPS_REL=:', SCARC_EPS_REL
   CALL SCARC_SHOW(SLHI%X, 'SARC',  'X1    ')
   CALL SCARC_SHOW(SLHI%F, 'SARC',  'F1    ')
ENDIF
 
!!! initialize some method parameters
NREQ_FACE = 0        

BCONV =.FALSE.
BDIVG =.FALSE.
ICYCLE=1      ! V-cycle
 
!!! adjust boundary values of right hand side
CALL SCARC_SETBDRY3D (ILMAX, NM)
 
!!! ------------------------------------------------------------------------------------
!!! Only one level: solve problem exactly
!!! ------------------------------------------------------------------------------------
ONLY_ONE_LEVEL_IF: IF (S%NLEVEL==1) THEN   

      CALL SCARC_COARSE3D(ILMAX, NM)

!!! ------------------------------------------------------------------------------------
!!! More than one level: start MG-cycling
!!! ------------------------------------------------------------------------------------
ELSE

   !!! save cycle counts for MG-iteration
   S%KCYCLE(2,S%NLMAX)=1
   DO ILEVEL=S%NLMIN+1,S%NLMAX-1
      IF (ICYCLE==0) THEN
         S%KCYCLE(2,ILEVEL)=2
      ELSE
         S%KCYCLE(2,ILEVEL)=ICYCLE
      ENDIF
   ENDDO
     
   !!! calculate initial defect on finest level, Y = B - A*X, and get l2-norm of it
   ILEVEL=ILMAX
   CALL SCARC_MATVEC3D (SLHI%X, SLHI%D, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, NMV_X, NMV_D)  
   SLHI%D = SLHI%F - SLHI%D

   SLHI%RESIN_MG = SCARC_L2NORM3D (SLHI%D , ILEVEL, NM, NTYPE_GLOBAL)
   IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) 0, ILMAX, SL%RESIN_MG
   IF (SCARC_DEBUG >= 1.AND.NM==1) write (*,1000) 0, ILMAX, SL%RESIN_MG

   !CALL SCARC_MASTER_INFO (SLHI%ITE, SLHI%RESIN_MG)

    
   !!! ---------------------------------------------------------------------------
   !!! start MG-iteration
   !!! ---------------------------------------------------------------------------
   MG_LOOP3D: DO ITE = 1, SCARC_MG_NIT
    
IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,*) '============ HALLO, MG-Iteration ',ITE

      !!! set level-information
      DO ILEVEL=S%NLMIN,S%NLMAX
         S%KCYCLE(1,ILEVEL)=S%KCYCLE(2,ILEVEL)
      ENDDO

      ILEVEL=S%NLMAX

      !!! ---------------------------------------------------------------------------
      !!!  Presmoothing
      !!! ---------------------------------------------------------------------------
110   IF (ILEVEL/=S%NLMIN) THEN
    
         SLHI => S%SLEVEL(ILEVEL)
         SLLO => S%SLEVEL(ILEVEL-1)

         !IF (SCARC_METHOD>=2) THEN
         !   WRITE(SCARC_LU,*) 'SLHI zeigt auf level ',ILEVEL
         !   WRITE(SCARC_LU,*) 'SLLO zeigt auf level ',ILEVEL-1
         !ENDIF

         IF (ILEVEL==S%NLMAX) THEN
            BMATVEC=.FALSE.
         ELSE
            BMATVEC=.TRUE.
         ENDIF

         CALL SCARC_SMOOTHER3D(SLHI%X, SLHI%F, SLHI%D, BMATVEC, ILEVEL, NM)

         !!! print level information - still missing
         !CALL SCARC_INFO_LEVEL(ILEVEL)

         !!! perform restriction to coarser grid: F:= rest(D)
         CALL SCARC_RESTRICTION3D(SLLO%F, SLHI%D, ILEVEL-1, NM)

         !!! initialize solution vector and set boundary conditions to residuum - necessary ?
         SLLO%X=0.0_EB
         !CALL SCARC_BDRY_RESIDUUM(NM, SLLO%F, ILEVEL-1)

         !!! decrease level
         ILEVEL=ILEVEL-1

         GOTO 110
      
      ENDIF         

      !!! ------------------------------------------------------------------------
      !!! Coarse grid solver
      !!! ------------------------------------------------------------------------
      ILEVEL=S%NLMIN
      CALL SCARC_COARSE3D(ILEVEL, NM)


      !!! ------------------------------------------------------------------------
      !!! Postsmoothing
      !!! ------------------------------------------------------------------------
120   IF (ILEVEL/=S%NLMAX) THEN
    
         BMATVEC=.TRUE.

         SLLO => S%SLEVEL(ILEVEL)

         ILEVEL=ILEVEL+1
         SLHI => S%SLEVEL(ILEVEL)

         !!! print level information - still missing
         !CALL SCARC_INFO_LEVEL(ILEVEL)

         !!! perform prolongation to finer grid: D:=prol(X)
         CALL SCARC_PROLONGATION3D(SLLO%X, SLHI%D, ILEVEL, NM)

         !!! set exterior boundary data of residuum to zero - necessary ?
         !CALL SCARC_BDRY_RESIDUUM(NM, SLHI%D, ILEVEL)
 
         !!! set new solution
         SLHI%X = SLHI%D + SLHI%X

         !!! select smoother for postsmoothing
         CALL SCARC_SMOOTHER3D(SLHI%X, SLHI%F, SLHI%D, BMATVEC, ILEVEL, NM)


         !!! save cycle counts
         S%KCYCLE(1,ILEVEL)=S%KCYCLE(1,ILEVEL)-1
 
         IF (S%KCYCLE(1,ILEVEL)==0) THEN
            IF (ICYCLE==0) THEN
               S%KCYCLE(1,ILEVEL)=1
            ELSE
               S%KCYCLE(1,ILEVEL)=S%KCYCLE(2,ILEVEL)
            ENDIF
            GOTO 120
         ELSE
            GOTO 110
         ENDIF

      ENDIF

      !!! ------------------------------------------------------------------------
      !!! Defect calculation
      !!! ------------------------------------------------------------------------
      CALL SCARC_MATVEC3D (SLHI%X, SLHI%D, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, NMV_X, NMV_D)
      SLHI%D = SLHI%F - SLHI%D
      SLHI%RES_MG = SCARC_L2NORM3D (SLHI%D, ILEVEL, NM, NTYPE_GLOBAL)

      IF (SCARC_DEBUG>=1) WRITE(SCARC_LU,*) 'SCARC-Multigrid, iteration ',ITE,': residuum=',SLHI%RES_MG
      IF (SCARC_DEBUG>=2.AND.NM==1) WRITE(*,*) 'SCARC-Multigrid, iteration ',ITE,': residuum=',SLHI%RES_MG

      !!! ------------------------------------------------------------------------
      !!! Send error information to Master  - still missing
      !!! ------------------------------------------------------------------------
      !CALL SCARC_INFO_MASTER 
      !CALL SCARC_INFO_LEVEL(ILEVEL)

      !!! ------------------------------------------------------------------------
      !!! Convergence or divergence ?
      !!! ------------------------------------------------------------------------
      IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) ITE, ILEVEL, SL%RES_MG
      IF (SCARC_DEBUG >= 1.AND.NM==1) write (*,1000) ITE, ILEVEL, SL%RES_MG
      IF (VELOCITY_ERROR_FILE) WRITE(LU_VELOCITY_ERROR,'(2(I5,A),E16.8)') ICYC,', ',ITE,', ',SL%RES_MG
      IF (SCARC_BREL) THEN
         IF (SL%RES_MG <= SL%RESIN_MG*SCARC_MG_EPS) BCONV = .TRUE.
      ELSE
         IF (SL%RES_MG <= SCARC_MG_EPS .AND. SL%RES_MG <= SL%RESIN_MG*SCARC_EPS_REL) BCONV = .TRUE.
      ENDIF
      IF (SL%RES_MG > SCARC_EPS_DIVG) BDIVG = .TRUE.
      IF (BCONV.OR.BDIVG) EXIT MG_LOOP3D
 
   ENDDO MG_LOOP3D

   !!! ------------------------------------------------------------------------
   !!! Convergence or divergence ?
   !!! ------------------------------------------------------------------------
   IF (BDIVG) THEN                       
      ITE0 = - 1
      SL%CAPPA_MG = 1.0_EB
   ELSE 
      IF (BCONV) THEN
         ITE0=ITE
      ELSE
         ITE0=ITE-1
      ENDIF
      IF (SL%RESIN_MG >= 1.0E-70_EB) THEN
         SL%RESIN_MG = SL%RES_MG / SL%RESIN_MG
         SL%CAPPA_MG = SL%RESIN_MG ** (1.0_EB/ITE0)
      ELSE 
         SL%CAPPA_MG = 0.0_EB
      ENDIF
   ENDIF

   SCARC_NIT=ITE0
   SCARC_RES=SL%RES_MG
   SCARC_CAPPA=SL%CAPPA_MG
   IF (SCARC_DEBUG >= 1) WRITE (SCARC_LU,2000) ITE0, ILEVEL, SL%RES_MG, SL%CAPPA_MG 
   IF (SCARC_DEBUG >= 1.AND.NM==1) WRITE (*,2000) ITE0, ILEVEL, SL%RES_MG, SL%CAPPA_MG 

ENDIF ONLY_ONE_LEVEL_IF
 
! shift solution back to vector HP
DO K = 1, SLHI%KBAR                                 
   DO J = 1, SLHI%JBAR
      DO I = 1, SLHI%IBAR
         HP(I, J, K) = SLHI%X (I, J, K)
      ENDDO
   ENDDO
ENDDO

! set ghost cell values along exterior boundaries 
CALL SCARC_GHOSTCELLS(HP,NM)

! get neighbouring data along interior boundaries
CALL SCARC_UPDATE_LEVEL(HP, NM, ILMAX, 'HPFIN ')

 
TUSED_SCARC(7,NM)=TUSED_SCARC(7,NM)+SECOND()-TNOW_MG3D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_MG3D

1000 FORMAT ('=====SCARC_MG       : #ite= ',i4,': level=',i4,': res=',e16.8)
2000 FORMAT ('=====SCARC_MG       : #ite= ',i4,': level=',i4,': res=',e16.8,': rate=',e16.8,/)
END SUBROUTINE SCARC_MG3D



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Perform global 3D cg-method based on global possion-matrix
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_COARSE3D (ILEVEL, NM)
 
INTEGER ::  NM, ILEVEL, ITE, ITE0
REAL (EB) :: SIGMA0, SIGMA1, ALPHA, GAMMA
REAL (EB) :: TNOW_COARSE
LOGICAL BCONV, BDIVG
TYPE (MESH_TYPE),  POINTER :: M
 
TNOW_COARSE = SECOND()
BCONV = .FALSE.
BDIVG = .FALSE.

S  => SCARC (NM)
M  => MESHES (NM)

! cg only works in finest grid level
SL    => S%SLEVEL(ILEVEL)
 
! re-initialize auxiliary vectors
SL%X = 0.0_EB                   
SL%G = 0.0_EB                   
SL%Y = 0.0_EB
SL%R = 0.0_EB
SL%D = 0.0_EB
 
!!! calculate initial defect Y = B - A*X and get l2-norm of it
CALL SCARC_MATVEC3D (SL%X, SL%R, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, NMV_X, NMV_R)
SL%R = -SL%F + SL%R
SL%RESIN_CO = SCARC_L2NORM3D (SL%R, ILEVEL, NM, NTYPE_GLOBAL)
IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) 0, ILEVEL, SL%RESIN_CO
IF (SCARC_DEBUG >= 1.AND.NM==1) write (*,1000) 0, ILEVEL, SL%RESIN_CO
 
! initial preconditioning
SELECT CASE(SCARC_CO_PRECON)
   CASE('NONE') 
      SIGMA0 = SL%RESIN_CO * SL%RESIN_CO
   CASE('JACOBI') 
      SL%G=SL%R
      CALL SCARC_PRECON_JACOBI3D (SL%G, ILEVEL, NM)
      SIGMA0 = SCARC_SCALPROD3D (SL%R, SL%G, ILEVEL, NM)
   CASE('SSOR') 
      SL%G=SL%R
      CALL SCARC_PRECON_SSOR3D (SL%G, ILEVEL, NM)
      SIGMA0 = SCARC_SCALPROD3D (SL%R, SL%G, ILEVEL, NM)
   CASE('GSTRIX') 
      SL%G=SL%R
      CALL SCARC_PRECON_GSTRIX3D (SL%G, ILEVEL, NM)
      SIGMA0 = SCARC_SCALPROD3D (SL%R, SL%G, ILEVEL, NM)
END SELECT
SL%D = -SL%G 
 
!
! start defect correction loop
!
CO_LOOP3D: DO ITE = 1, SCARC_CO_NIT
 
   ! calculate new defect and get L2-norm of it
   CALL SCARC_MATVEC3D (SL%D, SL%Y, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, NMV_D, NMV_Y)
   ALPHA = SCARC_SCALPROD3D (SL%D, SL%Y, ILEVEL, NM)

   ! get new descent direction
   ALPHA = SIGMA0 / ALPHA
   SL%X = ALPHA * SL%D + SL%X
   SL%R = ALPHA * SL%Y + SL%R
   SL%RES_CO = SCARC_L2NORM3D (SL%R, ILEVEL, NM, NTYPE_GLOBAL)
 
   !  exit loop in case of convergence or divergence
   IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) ITE, ILEVEL, SL%RES_CO
   IF (SCARC_DEBUG >= 1.AND.NM==1) write (*,1000) ITE, ILEVEL, SL%RES_CO
   IF (SCARC_BREL) THEN
      IF (SL%RES_CO <= SL%RESIN_CO*SCARC_CO_EPS) BCONV = .TRUE.
   ELSE
      IF (SL%RES_CO <= SCARC_CO_EPS .AND. SL%RES_CO <= SL%RESIN_CO*SCARC_EPS_REL) BCONV = .TRUE.
   ENDIF
   IF (SL%RES_CO > SCARC_EPS_DIVG) BDIVG = .TRUE.
   IF (BCONV.OR.BDIVG) EXIT CO_LOOP3D
 
   ! preconditioning
   SELECT CASE(SCARC_CO_PRECON)
      CASE('NONE')
         SIGMA1 = SL%RES_CO * SL%RES_CO
      CASE('JACOBI')
         SL%G=SL%R
         CALL SCARC_PRECON_JACOBI3D (SL%G, ILEVEL, NM)
         SIGMA1 = SCARC_SCALPROD3D (SL%R, SL%G, ILEVEL, NM)
      CASE('SSOR')
         SL%G=SL%R
         CALL SCARC_PRECON_SSOR3D (SL%G, ILEVEL, NM)
         SIGMA1 = SCARC_SCALPROD3D (SL%R, SL%G, ILEVEL, NM)
      CASE('GSTRIX')
         SL%G=SL%R
         CALL SCARC_PRECON_GSTRIX3D (SL%G, ILEVEL, NM)
         SIGMA1 = SCARC_SCALPROD3D (SL%R, SL%G, ILEVEL, NM)
   END SELECT
 
   GAMMA = SIGMA1 / SIGMA0
   SIGMA0 = SIGMA1
   SL%D = -SL%G + GAMMA * SL%D
 
ENDDO CO_LOOP3D
 
! divergence or convergence ?
IF (BDIVG) THEN                       
   ITE0 = - 1
   SL%CAPPA_CO = 1.0_EB
ELSE 
   IF (BCONV) THEN
      ITE0=ITE
   ELSE
      ITE0=ITE-1
   ENDIF
   IF (SL%RESIN_CO >= 1.0E-70_EB) THEN
      SL%RESIN_CO = SL%RES_CO / SL%RESIN_CO
      SL%CAPPA_CO = SL%RESIN_CO ** (1.0_EB/ITE0)
   ELSE 
      SL%CAPPA_CO = 0.0_EB
   ENDIF
ENDIF
IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,2000) ITE0, ILEVEL, SL%RES_CO, SL%CAPPA_CO 
IF (SCARC_DEBUG >= 1.AND.NM==1) WRITE (*,2000) ITE0, ILEVEL, SL%RES_CO, SL%CAPPA_CO 

 
TUSED_SCARC(9,NM)=TUSED_SCARC(9,NM)+SECOND()-TNOW_COARSE
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_COARSE

1000 FORMAT ('     SCARC_COARSE   : #ite= ',i4,': level=',i4,': res=',e16.8)
2000 FORMAT ('     SCARC_COARSE   : #ite= ',i4,': level=',i4,': res=',e16.8,': rate=',e16.8)
END SUBROUTINE SCARC_COARSE3D



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Perform global 3D cg-method based on global possion-matrix
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PRECON_COARSE3D (ILEVEL, NM)
 
INTEGER ::  NM, ILEVEL, ITE, ITE0
REAL (EB) :: SIGMA0, SIGMA1, ALPHA, GAMMA
REAL (EB) :: TNOW_COARSE
LOGICAL BCONV, BDIVG
TYPE (MESH_TYPE),  POINTER :: M
 
TNOW_COARSE = SECOND()
BCONV = .FALSE.
BDIVG = .FALSE.

S  => SCARC (NM)
M  => MESHES (NM)

! cg only works in finest grid level
SL    => S%SLEVEL(ILEVEL)
 
! re-initialize auxiliary vectors
SL%X2 = 0.0_EB                   
SL%G2 = 0.0_EB                   
SL%Y2 = 0.0_EB
SL%R2 = 0.0_EB
SL%D2 = 0.0_EB
 
!!! calculate initial defect Y = B - A*X and get l2-norm of it
CALL SCARC_MATVEC3D (SL%X2, SL%R2, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, NMV_X2, NMV_R2)
SL%R2 = -SL%F2 + SL%R2
SL%RESIN_CO = SCARC_L2NORM3D (SL%R2, ILEVEL, NM, NTYPE_GLOBAL)
IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) 0, ILEVEL, SL%RESIN_CO
IF (SCARC_DEBUG >= 2.AND.NM==1) write (*,1000) 0, ILEVEL, SL%RESIN_CO
 
! initial preconditioning
SELECT CASE(SCARC_CO_PRECON)
   CASE('NONE') 
      SIGMA0 = SL%RESIN_CO * SL%RESIN_CO
   CASE('JACOBI') 
      SL%G2=SL%R2
      CALL SCARC_PRECON_JACOBI3D (SL%G2, ILEVEL, NM)
      SIGMA0 = SCARC_SCALPROD3D (SL%R2, SL%G2, ILEVEL, NM)
   CASE('SSOR') 
      SL%G2=SL%R2
      CALL SCARC_PRECON_SSOR3D (SL%G2, ILEVEL, NM)
      SIGMA0 = SCARC_SCALPROD3D (SL%R2, SL%G2, ILEVEL, NM)
   CASE('GSTRIX') 
      SL%G2=SL%R2
      CALL SCARC_PRECON_GSTRIX3D (SL%G2, ILEVEL, NM)
      SIGMA0 = SCARC_SCALPROD3D (SL%R2, SL%G2, ILEVEL, NM)
END SELECT
SL%D2 = -SL%G2
 
!
! start defect correction loop
!
CO_LOOP3D: DO ITE = 1, SCARC_CO_NIT
 
   ! calculate new defect and get L2-norm of it
   CALL SCARC_MATVEC3D (SL%D2, SL%Y2, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, NMV_D2, NMV_Y2)
   ALPHA = SCARC_SCALPROD3D (SL%D2, SL%Y2, ILEVEL, NM)

   ! get new descent direction
   ALPHA = SIGMA0 / ALPHA
   SL%X2 = ALPHA * SL%D2 + SL%X2
   SL%R2 = ALPHA * SL%Y2 + SL%R2
   SL%RES_CO = SCARC_L2NORM3D (SL%R2, ILEVEL, NM, NTYPE_GLOBAL)
 
   !  exit loop in case of convergence or divergence
   IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) ITE, ILEVEL, SL%RES_CO
   IF (SCARC_DEBUG >= 2.AND.NM==1) write (*,1000) ITE, ILEVEL, SL%RES_CO
   IF (SCARC_BREL) THEN
      IF (SL%RES_CO <= SL%RESIN_CO*SCARC_CO_EPS) BCONV = .TRUE.
   ELSE
      IF (SL%RES_CO <= SCARC_CO_EPS .AND. SL%RES_CO <= SL%RESIN_CO*SCARC_EPS_REL) BCONV = .TRUE.
   ENDIF
   IF (SL%RES_CO > SCARC_EPS_DIVG) BDIVG = .TRUE.
   IF (BCONV.OR.BDIVG) EXIT CO_LOOP3D
 
   ! preconditioning
   SELECT CASE(SCARC_CO_PRECON)
      CASE('NONE')
         SIGMA1 = SL%RES_CO * SL%RES_CO
      CASE('JACOBI')
         SL%G2=SL%R2
         CALL SCARC_PRECON_JACOBI3D (SL%G2, ILEVEL, NM)
         SIGMA1 = SCARC_SCALPROD3D (SL%R2, SL%G2, ILEVEL, NM)
      CASE('SSOR')
         SL%G2=SL%R2
         CALL SCARC_PRECON_SSOR3D (SL%G2, ILEVEL, NM)
         SIGMA1 = SCARC_SCALPROD3D (SL%R2, SL%G2, ILEVEL, NM)
      CASE('GSTRIX')
         SL%G2=SL%R2
         CALL SCARC_PRECON_GSTRIX3D (SL%G2, ILEVEL, NM)
         SIGMA1 = SCARC_SCALPROD3D (SL%R2, SL%G2, ILEVEL, NM)
   END SELECT
 
   GAMMA = SIGMA1 / SIGMA0
   SIGMA0 = SIGMA1
   SL%D2 = -SL%G2 + GAMMA * SL%D2
 
ENDDO CO_LOOP3D
 
! divergence or convergence ?
IF (BDIVG) THEN                       
   ITE0 = - 1
   SL%CAPPA_CO = 1.0_EB
ELSE 
   IF (BCONV) THEN
      ITE0=ITE
   ELSE
      ITE0=ITE-1
   ENDIF
   IF (SL%RESIN_CO >= 1.0E-70_EB) THEN
      SL%RESIN_CO = SL%RES_CO / SL%RESIN_CO
      SL%CAPPA_CO = SL%RESIN_CO ** (1.0_EB/ITE0)
   ELSE 
      SL%CAPPA_CO = 0.0_EB
   ENDIF
ENDIF
IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,2000) ITE0, ILEVEL, SL%RES_CO, SL%CAPPA_CO 
IF (SCARC_DEBUG >= 2.AND.NM==1) WRITE (*,2000) ITE0, ILEVEL, SL%RES_CO, SL%CAPPA_CO 

 
TUSED_SCARC(9,NM)=TUSED_SCARC(9,NM)+SECOND()-TNOW_COARSE
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_COARSE

1000 FORMAT ('     SCARC_PRECON_COARSE3D   : #ite= ',i4,': level=',i4,': res=',e16.8)
2000 FORMAT ('     SCARC_PRECON_COARSE3D   : #ite= ',i4,': level=',i4,': res=',e16.8,': rate=',e16.8)
END SUBROUTINE SCARC_PRECON_COARSE3D


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Perform smoothing in SCARC_MG3D
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SMOOTHER3D(X, F, D, BMATVEC, ILEVEL, NM)

REAL (EB), POINTER, DIMENSION (:, :, :) :: X, F, D
INTEGER :: NM
INTEGER :: ILEVEL, ITE, ITE0
REAL(EB):: TNOW_SMOOTHER3D
LOGICAL :: BMATVEC, BL2NORM=.TRUE.       ! necessary ???
LOGICAL :: BCONV, BDIVG


TNOW_SMOOTHER3D = SECOND()

SL => S%SLEVEL(ILEVEL)
BCONV=.FALSE.
BDIVG=.FALSE.

SL%RESIN_SM=1.0_EB
IF (BMATVEC) THEN
   CALL SCARC_MATVEC3D (X, D, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, NMV_X, NMV_D)  
   D = F - D
   IF (BL2NORM) SL%RESIN_SM = SCARC_L2NORM3D (D, ILEVEL, NM, NTYPE_GLOBAL)
ENDIF

IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) 0, ILEVEL, SL%RESIN_SM
IF (SCARC_DEBUG >= 3.AND.NM==1) WRITE(*,1000) 0, ILEVEL, SL%RESIN_SM

SMOOTH_LOOP3D: DO ITE=1,SCARC_SM_NIT
 
   SELECT CASE(SCARC_SM_PRECON)
      CASE('JACOBI') 
         CALL SCARC_PRECON_JACOBI3D (D, ILEVEL, NM)
      CASE('SSOR') 
         CALL SCARC_PRECON_SSOR3D (D, ILEVEL, NM)
      CASE('GSTRIX') 
         CALL SCARC_PRECON_GSTRIX3D (D, ILEVEL, NM)
      CASE('FFT') 
         CALL SCARC_PRECON_FFT3D (D, D, ILEVEL, NM)
   END SELECT

   X = SCARC_SM_OMEGA * D + X
   CALL SCARC_MATVEC3D (X, D, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, NMV_X, NMV_D)   ! global MATVEC
   D = F - D
   IF (BL2NORM) SL%RES_SM = SCARC_L2NORM3D (D, ILEVEL, NM, NTYPE_GLOBAL)

   IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) ITE, ILEVEL, SL%RES_SM
   IF (SCARC_DEBUG >= 2.AND.NM==1) WRITE(*,1000) ITE, ILEVEL, SL%RES_SM
   IF (SCARC_BREL) THEN
      IF (SL%RES_SM <= SL%RESIN_SM*SCARC_SM_EPS) BCONV = .TRUE.
   ELSE
      IF (SL%RES_SM <= SCARC_SM_EPS .AND. SL%RES_SM <= SL%RESIN_SM*SCARC_EPS_REL) BCONV = .TRUE.
   ENDIF
   IF (SL%RES_SM > SCARC_EPS_DIVG) BDIVG = .TRUE.
   IF (BCONV.OR.BDIVG) EXIT SMOOTH_LOOP3D
 

ENDDO SMOOTH_LOOP3D

! divergence or convergence ?
IF (BDIVG) THEN                       
   ITE0 = - 1
   SL%CAPPA_SM = 1.0_EB
ELSE 
   IF (BCONV) THEN
     ITE0=ITE
   ELSE
     ITE0=ITE-1
   ENDIF
   IF (SL%RESIN_SM >= 1.0E-70_EB) THEN
      SL%RESIN_SM = SL%RES_SM / SL%RESIN_SM
      SL%CAPPA_SM = SL%RESIN_SM ** (1.0_EB/ITE0)
   ELSE 
      SL%CAPPA_SM = 0.0_EB
   ENDIF
ENDIF
IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,2000) ITE0, ILEVEL, SL%RES_SM, SL%CAPPA_SM 
IF (SCARC_DEBUG >= 3.AND.NM==1) WRITE (*,2000) ITE0, ILEVEL, SL%RES_SM, SL%CAPPA_SM 

TUSED_SCARC(8,NM)=TUSED_SCARC(8,NM)+SECOND()-TNOW_SMOOTHER3D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_SMOOTHER3D

1000 FORMAT ('     SCARC_SMOOTHER : #ite= ',i4,': level=',i4,': res=',e16.8)
2000 FORMAT ('     SCARC_SMOOTHER : #ite= ',i4,': level=',i4,': res=',e16.8,': rate=',e16.8)
END SUBROUTINE SCARC_SMOOTHER3D



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Perform smoothing in SCARC_MG3D
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PRECON_SMOOTHER3D(X, F, D, BMATVEC, ILEVEL, NM)

REAL (EB), POINTER, DIMENSION (:, :, :) :: X, F, D
INTEGER :: NM
INTEGER :: ILEVEL, ITE, ITE0
REAL(EB):: TNOW_PRECON_SMOOTHER3D
LOGICAL :: BMATVEC, BL2NORM=.TRUE.       ! necessary ???
LOGICAL :: BCONV, BDIVG


TNOW_PRECON_SMOOTHER3D = SECOND()

SL => S%SLEVEL(ILEVEL)
BCONV=.FALSE.
BDIVG=.FALSE.

SL%RESIN_SM=1.0_EB
IF (BMATVEC) THEN
   CALL SCARC_MATVEC3D (X, D, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, NMV_X2, NMV_D2)  
   D = F - D
   IF (BL2NORM) SL%RESIN_SM = SCARC_L2NORM3D (D, ILEVEL, NM, NTYPE_GLOBAL)
ENDIF

IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) 0, ILEVEL, SL%RESIN_SM
IF (SCARC_DEBUG >= 2.AND.NM==1) WRITE(*,1000) 0, ILEVEL, SL%RESIN_SM

SMOOTH_LOOP3D: DO ITE=1,SCARC_SM_NIT
 
   SELECT CASE(SCARC_SM_PRECON)
      CASE('JACOBI') 
         CALL SCARC_PRECON_JACOBI3D (D, ILEVEL, NM)
      CASE('SSOR') 
         CALL SCARC_PRECON_SSOR3D (D, ILEVEL, NM)
      CASE('GSTRIX') 
         CALL SCARC_PRECON_GSTRIX3D (D, ILEVEL, NM)
      CASE('FFT') 
         CALL SCARC_PRECON_FFT3D (D, D, ILEVEL, NM)
   END SELECT

   X = SCARC_SM_OMEGA * D + X
   CALL SCARC_MATVEC3D (X, D, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, NMV_X2, NMV_D2)   ! global MATVEC
   D = F - D
   IF (BL2NORM) SL%RES_SM = SCARC_L2NORM3D (D, ILEVEL, NM, NTYPE_GLOBAL)

   IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) ITE, ILEVEL, SL%RES_SM
   IF (SCARC_DEBUG >= 2.AND.NM==1) WRITE(*,1000) ITE, ILEVEL, SL%RES_SM
   IF (SCARC_BREL) THEN
      IF (SL%RES_SM <= SL%RESIN_SM*SCARC_SM_EPS) BCONV = .TRUE.
   ELSE
      IF (SL%RES_SM <= SCARC_SM_EPS .AND. SL%RES_SM <= SL%RESIN_SM*SCARC_EPS_REL) BCONV = .TRUE.
   ENDIF
   IF (SL%RES_SM > SCARC_EPS_DIVG) BDIVG = .TRUE.
   IF (BCONV.OR.BDIVG) EXIT SMOOTH_LOOP3D
 

ENDDO SMOOTH_LOOP3D

! divergence or convergence ?
IF (BDIVG) THEN                       
   ITE0 = - 1
   SL%CAPPA_SM = 1.0_EB
ELSE 
   IF (BCONV) THEN
     ITE0=ITE
   ELSE
     ITE0=ITE-1
   ENDIF
   IF (SL%RESIN_SM >= 1.0E-70_EB) THEN
      SL%RESIN_SM = SL%RES_SM / SL%RESIN_SM
      SL%CAPPA_SM = SL%RESIN_SM ** (1.0_EB/ITE0)
   ELSE 
      SL%CAPPA_SM = 0.0_EB
   ENDIF
ENDIF
IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,2000) ITE0, ILEVEL, SL%RES_SM, SL%CAPPA_SM 
IF (SCARC_DEBUG >= 2.AND.NM==1) WRITE (*,2000) ITE0, ILEVEL, SL%RES_SM, SL%CAPPA_SM 

TUSED_SCARC(8,NM)=TUSED_SCARC(8,NM)+SECOND()-TNOW_PRECON_SMOOTHER3D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_PRECON_SMOOTHER3D

1000 FORMAT ('     SCARC_PRECON_SMOOTHER3D : #ite= ',i4,': level=',i4,': res=',e16.8)
2000 FORMAT ('     SCARC_PRECON_SMOOTHER3D : #ite= ',i4,': level=',i4,': res=',e16.8,': rate=',e16.8)
END SUBROUTINE SCARC_PRECON_SMOOTHER3D



 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Perform matrix-vector multiplication with full matrix
!!! corresponding to band-wise storage technique
!!!
!!!  Y = A1 * S%AGLOB * X + A2 * Y
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_MATVEC2D (X, Y, A1, A2, NM, ITYPE, ILEVEL, IMV1, IMV2)
 
REAL (EB) :: A1, A2
REAL (EB), POINTER, DIMENSION (:, :, :) :: X, Y
REAL (EB):: TNOW_MATVEC2D
INTEGER :: NM, ITYPE
INTEGER :: I, K, IC, IMV1, IMV2, ILEVEL, II
 
TNOW_MATVEC2D = SECOND()
 
SL => S%SLEVEL(ILEVEL)
 
IF (SCARC_DEBUG>=4) CALL SCARC_SHOW_LEVEL (X, 'SARC', 'XMAT0 ',ILEVEL)
IF (SCARC_DEBUG>=4) CALL SCARC_SHOW_LEVEL (Y, 'SARC', 'YMAT0 ',ILEVEL)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform 2D-matrix-vector-multiplication
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DO K = 1, SL%KBAR
   DO I = 1, SL%IBAR
      IC = (K-1) * SL%IBAR + I

      IF (K-1<1      .and.X(I,1,K-1).gt.1.0E-14_EB) then
        WRITE(SCARC_LU,*) 'A:WRONG: X(',I,',1,',K-1,')=',X(I,1,K-1), SL%AG(IC,2)
        !stop
      ENDIF
      IF (K+1>SL%KBAR.and.X(I,1,K+1).gt.1.0E-14_EB) then
        WRITE(SCARC_LU,*) 'B:WRONG: X(',I,',1,',K+1,')=',X(I,1,K+1), SL%AG(IC,5)
        !stop
      ENDIF
      IF (I-1<1      .and.X(I-1,1,K).gt.1.0E-14_EB) then
        WRITE(SCARC_LU,*) 'C:WRONG: X(',I-1,',1,',K,')=',X(I-1,1,K), SL%AG(IC,3)
        !stop
      ENDIF
      IF (I+1>SL%IBAR.and.X(I+1,1,K).gt.1.0E-14_EB) then
        WRITE(SCARC_LU,*) 'D:WRONG: X(',I+1,',1,',K,')=',X(I+1,1,K), SL%AG(IC,4)
        !stop
      ENDIF

      IF (K-1<=0.and.X(I,1,K-1).gt.1.0E-14_EB) stop
      Y (I, 1, K) =    A1 * (  SL%AG(IC, 1)*X(I  , 1, K  ) &
                             + SL%AG(IC, 2)*X(I  , 1, K-1) &
                             + SL%AG(IC, 3)*X(I-1, 1, K  ) &
                             + SL%AG(IC, 4)*X(I+1, 1, K  ) &
                             + SL%AG(IC, 5)*X(I  , 1, K+1))&
                     + A2 * Y (I, 1, K)
      IF (SCARC_DEBUG>=6) &
         WRITE(SCARC_LU,'(2i3,11e11.3)') I,K,Y(I,1,K),(SL%AG(IC,II),II=1,5),&
                                  X(I,1,K),X(I,1,K-1),X(I-1,1,K),X(I+1,1,K),X(I,1,K+1)
   ENDDO
ENDDO
 
IF (SCARC_DEBUG>=4) CALL SCARC_SHOW_LEVEL (X, 'SARC', 'XMAT1 ',ILEVEL)
IF (SCARC_DEBUG>=4) CALL SCARC_SHOW_LEVEL (Y, 'SARC', 'YMAT1 ',ILEVEL)
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform data exchange at interior boundaries
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF ((ITYPE == NTYPE_GLOBAL) .AND. (NMESHES > 1)) THEN
   NREQ_FACE = 0
   CALL SCARC_RECEIVE  (NCOM_MATV,  ILEVEL)   
   CALL SCARC_EXCHANGE (NCOM_MATV,  ILEVEL, IMV1, IMV2)
ENDIF
 
IF (SCARC_DEBUG>=4) CALL SCARC_SHOW_LEVEL (X, 'SARC', 'XMAT2 ',ILEVEL)
IF (SCARC_DEBUG>=4) CALL SCARC_SHOW_LEVEL (Y, 'SARC', 'YMAT2 ',ILEVEL)
 
TUSED_SCARC(10,NM)=TUSED_SCARC(10,NM)+SECOND()-TNOW_MATVEC2D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_MATVEC2D

END SUBROUTINE SCARC_MATVEC2D
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Perform matrix-vector multiplication with full matrix
!!! corresponding to band-wise storage technique
!!!
!!!  Y = A1 * S%AGLOB * X + A2 * Y
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_MATVEC3D (X, Y, A1, A2, NM, ITYPE, ILEVEL, IMV1, IMV2)
 
REAL (EB) :: A1, A2
REAL (EB), POINTER, DIMENSION (:, :, :) :: X, Y
REAL (EB):: TNOW_MATVEC3D
INTEGER :: NM, ITYPE
INTEGER :: I, J, K, IC, IMV1, IMV2, ILEVEL, LL
 
 
TNOW_MATVEC3D = SECOND()
 
SL => S%SLEVEL(ILEVEL)
 
IF (SCARC_DEBUG>=4) CALL SCARC_SHOW_LEVEL (X, 'SARC', 'XMAT0 ',ILEVEL)
IF (SCARC_DEBUG>=4) CALL SCARC_SHOW_LEVEL (Y, 'SARC', 'YMAT0 ',ILEVEL)
IF (SCARC_DEBUG>=4) THEN
   write(SCARC_LU,*) 'IBAR=',SL%IBAR
   write(SCARC_LU,*) 'JBAR=',SL%JBAR
   write(SCARC_LU,*) 'KBAR=',SL%KBAR
ENDIF

    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform 3D-matrix-vector-multiplication
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DO K = 1, SL%KBAR
IF (SCARC_DEBUG .GE. 6) write (SCARC_LU,*) '-------------------------------------'

   DO J = 1, SL%JBAR
      DO I = 1, SL%IBAR
         IC = (K-1) * SL%IBAR * SL%JBAR + (J-1) * SL%IBAR + I
         Y(I, J, K) = A1 * ( SL%AG(IC, 1) *X (I  , J  , K  )   &
                         +   SL%AG(IC, 2) *X (I  , J  , K-1)   &
                         +   SL%AG(IC, 3) *X (I  , J-1, K  )   &
                         +   SL%AG(IC, 4) *X (I-1, J  , K  )   &
                         +   SL%AG(IC, 5) *X (I+1, J  , K  )   &
                         +   SL%AG(IC, 6) *X (I  , J+1, K  )   &
                         +   SL%AG(IC, 7) *X (I  , J  , K+1) ) &
                    + A2 * Y (I, J, K)
IF (SCARC_DEBUG .GE. 6) write (SCARC_LU, '(a,i2,a,i2,a,i2,a,i3,a,7f6.0,a,f22.16)') &
                       '(', I, ',', J, ',', K, '):A(', IC, '):', &
                       (SL%AG(IC, LL), LL=1, 7),  ' -> Y=', Y (I, J, K)

      ENDDO
IF (SCARC_DEBUG .GE. 6) write (SCARC_LU,*) '-------------------------------------'
   ENDDO
ENDDO
 
IF (SCARC_DEBUG>=4) THEN
   CALL SCARC_SHOW_LEVEL (Y, 'SARC', 'YMAT1 ',ILEVEL)
   WRITE(SCARC_LU,*) 'ITYPE=',ITYPE
   WRITE(SCARC_LU,*) 'ILEVEL=',ILEVEL
   WRITE(SCARC_LU,*) 'NMESHE=',NMESHES
ENDIF
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform data exchange at interior boundaries
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF ((ITYPE == NTYPE_GLOBAL) .AND. (NMESHES > 1)) THEN
   IF (SCARC_DEBUG>=4)  WRITE(SCARC_LU,*) 'COMMUNICATING',IMV1,IMV2
   NREQ_FACE = 0
   CALL SCARC_RECEIVE  (NCOM_MATV,  ILEVEL)  
   CALL SCARC_EXCHANGE (NCOM_MATV,  ILEVEL, IMV1, IMV2)
ENDIF
 
IF (SCARC_DEBUG>=4) CALL SCARC_SHOW_LEVEL (X, 'SARC', 'XMAT2 ',ILEVEL)
IF (SCARC_DEBUG>=4) CALL SCARC_SHOW_LEVEL (Y, 'SARC', 'YMAT2 ',ILEVEL)
 
TUSED_SCARC(10,NM)=TUSED_SCARC(10,NM)+SECOND()-TNOW_MATVEC3D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_MATVEC3D

END SUBROUTINE SCARC_MATVEC3D
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Restrict vector X from level ILEVEL to vector Y on level ILEVEL-1
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_RESTRICTION2D (X_LO, X_HI, ILEVEL_LO, NM)
 
REAL (EB), POINTER, DIMENSION (:, :, :) :: X_LO, X_HI
INTEGER :: ILEVEL_LO
INTEGER :: I_LO, K_LO, I1_HI, I2_HI, K1_HI, K2_HI,NM
REAL(EB):: TNOW_RESTRICTION2D

TNOW_RESTRICTION2D = SECOND()

SLLO => S%SLEVEL(ILEVEL_LO)

DO K_LO = 1, SLLO%KBAR
   DO I_LO = 1, SLLO%IBAR
      
      I1_HI=2*I_LO-1
      I2_HI=2*I_LO
      K1_HI=2*K_LO-1
      K2_HI=2*K_LO

      X_LO (I_LO, 1, K_LO) = 0.25_EB * (  X_HI(I1_HI,1,K1_HI)  &
                                        + X_HI(I2_HI,1,K1_HI)  &
                                        + X_HI(I1_HI,1,K2_HI)  &
                                        + X_HI(I2_HI,1,K2_HI) )

      IF (SCARC_DEBUG>=4) WRITE(SCARC_LU,1000) I_LO,1,K_LO,X_LO(I_LO,1,K_LO),I1_HI,I2_HI,K1_HI,K2_HI
   ENDDO
ENDDO

TUSED_SCARC(11,NM)=TUSED_SCARC(11,NM)+SECOND()-TNOW_RESTRICTION2D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_RESTRICTION2D
1000 FORMAT('X_LO(',I3,',',I3,',',I3,')=',f12.6,4i3)
END SUBROUTINE SCARC_RESTRICTION2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Restrict vector X from level ILEVEL to vector Y on level ILEVEL-1
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_RESTRICTION3D (X_LO, X_HI, ILEVEL_LO, NM)
 
REAL (EB), POINTER, DIMENSION (:, :, :) :: X_LO, X_HI
INTEGER :: ILEVEL_LO, NM
INTEGER :: I_LO, I1_HI, I2_HI
INTEGER :: J_LO, J1_HI, J2_HI
INTEGER :: K_LO, K1_HI, K2_HI
REAL(EB):: TNOW_RESTRICTION3D

TNOW_RESTRICTION3D = SECOND()

SLLO => S%SLEVEL(ILEVEL_LO)

DO K_LO = 1, SLLO%KBAR
   DO J_LO = 1, SLLO%JBAR
      DO I_LO = 1, SLLO%IBAR
      
         I1_HI=2*I_LO-1
         I2_HI=2*I_LO
         J1_HI=2*J_LO-1
         J2_HI=2*J_LO
         K1_HI=2*K_LO-1
         K2_HI=2*K_LO

         X_LO (I_LO, J_LO , K_LO) = 0.125_EB * (  X_HI(I1_HI,J1_HI,K1_HI)  &
                                                + X_HI(I2_HI,J1_HI,K1_HI)  &
                                                + X_HI(I1_HI,J2_HI,K1_HI)  &
                                                + X_HI(I2_HI,J2_HI,K1_HI)  &
                                                + X_HI(I1_HI,J1_HI,K2_HI)  &
                                                + X_HI(I2_HI,J1_HI,K2_HI)  &
                                                + X_HI(I1_HI,J2_HI,K2_HI)  &
                                                + X_HI(I2_HI,J2_HI,K2_HI) )
         IF (SCARC_DEBUG>=4) WRITE(SCARC_LU,1000) I_LO,J_LO,K_LO,X_LO(I_LO,J_LO,K_LO),I1_HI,I2_HI,J1_HI,J2_HI,K1_HI,K2_HI

      ENDDO
   ENDDO
ENDDO

TUSED_SCARC(11,NM)=TUSED_SCARC(11,NM)+SECOND()-TNOW_RESTRICTION3D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_RESTRICTION3D
1000 FORMAT('X_LO(',I3,',',I3,',',I3,')=',f25.16,6i3)
END SUBROUTINE SCARC_RESTRICTION3D


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Restrict vector X from level ILEVEL to vector Y on level ILEVEL-1
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PROLONGATION2D (X_LO, X_HI, ILEVEL_HI, NM)
 
REAL (EB), POINTER, DIMENSION (:, :, :) :: X_LO, X_HI
INTEGER :: ILEVEL_HI, NM
INTEGER :: I_LO, K_LO, I1_HI, I2_HI, K1_HI, K2_HI
REAL(EB):: TNOW_PROLONGATION2D

TNOW_PROLONGATION2D = SECOND()

SLLO => S%SLEVEL(ILEVEL_HI-1)
SLHI => S%SLEVEL(ILEVEL_HI)

CALL SCARC_SHOW_LEVEL (X_LO, 'SARC',  'PX_LO ',ILEVEL_HI-1)
DO K_LO = 1, SLLO%KBAR
   DO I_LO = 1, SLLO%IBAR
      
      I1_HI=2*I_LO-1
      I2_HI=2*I_LO
      K1_HI=2*K_LO-1
      K2_HI=2*K_LO

      X_HI (I1_HI, 1, K1_HI) = X_LO (I_LO, 1, K_LO)
      X_HI (I2_HI, 1, K1_HI) = X_LO (I_LO, 1, K_LO)
      X_HI (I1_HI, 1, K2_HI) = X_LO (I_LO, 1, K_LO)
      X_HI (I2_HI, 1, K2_HI) = X_LO (I_LO, 1, K_LO)

      IF (SCARC_DEBUG>=4) THEN
         WRITE(SCARC_LU,1000) I1_HI,1,K1_HI,X_LO(I_LO,1,K_LO)
         WRITE(SCARC_LU,1000) I2_HI,1,K1_HI,X_LO(I_LO,1,K_LO)
         WRITE(SCARC_LU,1000) I1_HI,1,K2_HI,X_LO(I_LO,1,K_LO)
         WRITE(SCARC_LU,1000) I2_HI,1,K2_HI,X_LO(I_LO,1,K_LO)
      ENDIF
   ENDDO
ENDDO
CALL SCARC_SHOW_LEVEL (X_HI, 'SARC',  'PX_LO ',ILEVEL_HI)

TUSED_SCARC(12,NM)=TUSED_SCARC(12,NM)+SECOND()-TNOW_PROLONGATION2D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_PROLONGATION2D
1000 FORMAT('X(',I3,',',I3,',',I3,')=',f25.16)
END SUBROUTINE SCARC_PROLONGATION2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Restrict vector X from level ILEVEL to vector Y on level ILEVEL-1
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PROLONGATION3D (X_LO, X_HI, ILEVEL_HI, NM)
 
REAL (EB), POINTER, DIMENSION (:, :, :) :: X_LO, X_HI
INTEGER :: ILEVEL_HI, NM 
INTEGER :: I_LO, I1_HI, I2_HI
INTEGER :: J_LO, J1_HI, J2_HI
INTEGER :: K_LO, K1_HI, K2_HI
REAL(EB):: TNOW_PROLONGATION3D

TNOW_PROLONGATION3D = SECOND()

SLLO => S%SLEVEL(ILEVEL_HI-1)
SLHI => S%SLEVEL(ILEVEL_HI)

CALL SCARC_SHOW_LEVEL (X_LO, 'SARC',  'PX_LO ',ILEVEL_HI-1)
DO K_LO = 1, SLLO%KBAR
   DO J_LO = 1, SLLO%JBAR
      DO I_LO = 1, SLLO%IBAR
      
         I1_HI=2*I_LO-1
         I2_HI=2*I_LO
         J1_HI=2*J_LO-1
         J2_HI=2*J_LO
         K1_HI=2*K_LO-1
         K2_HI=2*K_LO

         X_HI (I1_HI, J1_HI, K1_HI) = X_LO (I_LO, J_LO, K_LO)
         X_HI (I2_HI, J1_HI, K1_HI) = X_LO (I_LO, J_LO, K_LO)
         X_HI (I1_HI, J2_HI, K1_HI) = X_LO (I_LO, J_LO, K_LO)
         X_HI (I2_HI, J2_HI, K1_HI) = X_LO (I_LO, J_LO, K_LO)
         X_HI (I1_HI, J1_HI, K2_HI) = X_LO (I_LO, J_LO, K_LO)
         X_HI (I2_HI, J1_HI, K2_HI) = X_LO (I_LO, J_LO, K_LO)
         X_HI (I1_HI, J2_HI, K2_HI) = X_LO (I_LO, J_LO, K_LO)
         X_HI (I2_HI, J2_HI, K2_HI) = X_LO (I_LO, J_LO, K_LO)

         IF (SCARC_DEBUG>=4) THEN
            WRITE(SCARC_LU,1000) I1_HI,J1_HI,K1_HI,X_LO(I_LO,J_LO,K_LO)
            WRITE(SCARC_LU,1000) I2_HI,J1_HI,K1_HI,X_LO(I_LO,J_LO,K_LO)
            WRITE(SCARC_LU,1000) I1_HI,J2_HI,K1_HI,X_LO(I_LO,J_LO,K_LO)
            WRITE(SCARC_LU,1000) I2_HI,J2_HI,K1_HI,X_LO(I_LO,J_LO,K_LO)
            WRITE(SCARC_LU,1000) I1_HI,J1_HI,K2_HI,X_LO(I_LO,J_LO,K_LO)
            WRITE(SCARC_LU,1000) I2_HI,J1_HI,K2_HI,X_LO(I_LO,J_LO,K_LO)
            WRITE(SCARC_LU,1000) I1_HI,J2_HI,K2_HI,X_LO(I_LO,J_LO,K_LO)
            WRITE(SCARC_LU,1000) I2_HI,J2_HI,K2_HI,X_LO(I_LO,J_LO,K_LO)
         ENDIF

      ENDDO
   ENDDO
ENDDO
CALL SCARC_SHOW_LEVEL (X_HI, 'SARC',  'PX_LO ',ILEVEL_HI)

TUSED_SCARC(12,NM)=TUSED_SCARC(12,NM)+SECOND()-TNOW_PROLONGATION3D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_PROLONGATION3D
1000 FORMAT('X(',I3,',',I3,',',I3,')=',f12.6)
END SUBROUTINE SCARC_PROLONGATION3D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Text ??
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_BDRY_RESIDUUM (NM, F,  ILEVEL)
 
REAL (EB), POINTER, DIMENSION (:, :, :) :: F
INTEGER :: ILEVEL, INEWC, I, J, K, NM
REAL(EB):: TNOW_BDRY_RESIDUUM


TNOW_BDRY_RESIDUUM = SECOND()
SL => S%SLEVEL(ILEVEL)

DO INEWC=1,SL%N_EXTERNAL_WALL_CELLS

   IF (SL%IJKW(9,INEWC)==0) THEN

      I=SL%IJKW(1,INEWC)
      J=SL%IJKW(2,INEWC)
      K=SL%IJKW(3,INEWC)

      F(I,J,K)=0.0_EB

      if (SCARC_DEBUG>=4) WRITE(SCARC_LU,*) 'F(',I,',',J,',',K,')=',F(I,J,K)

   ENDIF

ENDDO

TUSED_SCARC(13,NM)=TUSED_SCARC(13,NM)+SECOND()-TNOW_BDRY_RESIDUUM
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_BDRY_RESIDUUM
END SUBROUTINE SCARC_BDRY_RESIDUUM



 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Compute scalarproduct of vector X and vector Y in 2D: SP=(X,Y)
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(EB) FUNCTION SCARC_SCALPROD2D (X, Y, ILEVEL, NM)
 
REAL (EB) :: SP, SPG
REAL (EB), POINTER, DIMENSION (:, :, :) :: X, Y
INTEGER :: ILEVEL, I, K,  NM, IERR
REAL(EB):: TNOW_SCALPROD2D
 
TNOW_SCALPROD2D = SECOND()
IERR=0
SL => S%SLEVEL(ILEVEL)
 
! compute local scalar product
SP = 0.0_EB
DO K = 1, SL%KBAR
   DO I = 1, SL%IBAR
      SP = SP + X (I, 1, K) * Y (I, 1, K)
      IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,1000) SP,I,1,K,X(I,1,K),Y(I,1,K)
   ENDDO
ENDDO
 
! compute global scalarproduct in case of multi-mesh computation 
IF (NMESHES > 1) THEN
   IF (USE_MPI) THEN
      CALL MPI_ALLREDUCE (SP, SPG, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
      SP = SPG
   ELSE
      WRITE(*,*) 'Serial version not yet implemented'
      stop
   ENDIF
ENDIF

SCARC_SCALPROD2D = SP
IF (SCARC_DEBUG>=3) WRITE(SCARC_LU,*) 'SCARC_SCALPROD2D=',SP

TUSED_SCARC(16,NM)=TUSED_SCARC(16,NM)+SECOND()-TNOW_SCALPROD2D
TUSED_SCARC( 0,NM)=TUSED_SCARC( 0,NM)+SECOND()-TNOW_SCALPROD2D
1000 FORMAT('SP=',f25.16,': (I,J,K)=',3i3,': X=',f25.16,': Y=',f25.16)
END FUNCTION SCARC_SCALPROD2D
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Compute scalarproduct of vector X and vector Y in 3D: SP=(X,Y)
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(EB) FUNCTION SCARC_SCALPROD3D (X, Y, ILEVEL, NM)
 
REAL (EB) :: SP, SPG
REAL (EB), POINTER, DIMENSION (:, :, :) :: X, Y
INTEGER :: ILEVEL, I, J, K, NM, IERR
REAL(EB):: TNOW_SCALPROD2D
 
TNOW_SCALPROD2D = SECOND()
IERR=0
SL => S%SLEVEL(ILEVEL)
 
! compute local scalar product
SP = 0.0_EB
DO K = 1, SL%KBAR
   DO J = 1, SL%JBAR
      DO I = 1, SL%IBAR
         SP = SP + X (I, J, K) * Y (I, J, K)
         IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,1000) SP,I,J,K,X(I,J,K),Y(I,J,K)
      ENDDO
   ENDDO
ENDDO
 
! compute global scalarproduct in case of multi-mesh computation 
IF (NMESHES > 1) THEN
   IF (USE_MPI) THEN
      CALL MPI_ALLREDUCE (SP, SPG, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
      SP = SPG
   ELSE
      WRITE(*,*) 'Serial version not yet implemented'
      stop
   ENDIF
ENDIF

SCARC_SCALPROD3D = SP
IF (SCARC_DEBUG>=3) WRITE(SCARC_LU,*) 'SCARC_SCALPROD3D=',SP

TUSED_SCARC(16,NM)=TUSED_SCARC(16,NM)+SECOND()-TNOW_SCALPROD2D
TUSED_SCARC( 0,NM)=TUSED_SCARC( 0,NM)+SECOND()-TNOW_SCALPROD2D
1000 FORMAT('SP=',f25.16,': (I,J,K)=',3i3,': X=',f25.16,': Y=',f25.16)
END FUNCTION SCARC_SCALPROD3D

 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Compute L2-norm of vector X:   SP = ||X||
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(EB) FUNCTION SCARC_L2NORM2D (X, ILEVEL, NM, ITYPE)
 
REAL (EB) :: SP, SPG
REAL (EB), POINTER, DIMENSION (:, :, :) :: X
INTEGER :: NM, ITYPE
INTEGER :: ILEVEL, IERR, I, K
REAL(EB):: TNOW_L2NORM2D
 
TNOW_L2NORM2D = SECOND()
IERR=0
SL => S%SLEVEL(ILEVEL)

!!! build local scalar product (x,x)
SP = 0.0_EB
DO K = 1, SL%KBAR
   DO I = 1, SL%IBAR
      SP = SP + X (I, 1, K) * X (I, 1, K)
      IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,1000) SP,I,1,K,X(I,1,K)
   ENDDO
ENDDO
 
!!! scale with number of cells
IF (ITYPE == NTYPE_LOCAL .OR. NMESHES==1) THEN                
   SP = Sqrt (SP/REAL(SL%NCELLS_LOCAL, EB))                   ! scale with local  number of cells
   IF (SCARC_DEBUG>=3) WRITE(SCARC_LU,*) 'Scale with ',SL%NCELLS_LOCAL,': SP=',SP
ELSE IF (ITYPE == NTYPE_GLOBAL) THEN                          
   IF (USE_MPI) THEN
      CALL MPI_ALLREDUCE (SP, SPG, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
      SP = Sqrt (SPG/REAL(SL%NCELLS_GLOBAL, EB))                 ! scale with global number of cells
      IF (SCARC_DEBUG>=3) WRITE(SCARC_LU,*) 'Scale with ',SL%NCELLS_GLOBAL,': SP=',SP, ': SPG=',SPG
   ELSE
      WRITE(*,*) 'Serial version not yet implemented'
      stop
   ENDIF
ELSE
   WRITE (*,*) 'Wrong type for SCARC_L2NORM2D ', ITYPE
ENDIF

SCARC_L2NORM2D = SP
IF (SCARC_DEBUG>=3) WRITE(SCARC_LU,*) 'SCARC_L2NORM2D=',SP

TUSED_SCARC(17,NM)=TUSED_SCARC(17,NM)+SECOND()-TNOW_L2NORM2D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_L2NORM2D
1000 FORMAT('SP=',f25.16,': (I,J,K)=',3i3,': X=',f25.16)
END FUNCTION SCARC_L2NORM2D
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Compute L2-norm of vector X:   SP = ||X||
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(EB) FUNCTION SCARC_L2NORM3D (X, ILEVEL, NM, ITYPE)
 
REAL (EB) :: SP, SPG
REAL (EB), POINTER, DIMENSION (:, :, :) :: X
INTEGER :: NM, ITYPE
INTEGER :: ILEVEL, IERR, I, J, K
REAL(EB):: TNOW_L2NORM3D
 
TNOW_L2NORM3D = SECOND()
IERR=0
SL => S%SLEVEL(ILEVEL)

CALL SCARC_SHOW_LEVEL (X, 'SARC',  'L2-X  ',ILEVEL)

!!! build local scalar product (x,x)
SP = 0.0_EB
DO K = 1, SL%KBAR
   DO J = 1, SL%JBAR
      DO I = 1, SL%IBAR
         SP = SP + X (I, J, K) * X (I, J, K)
         IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,1000) SP,I,J,K,X(I,J,K)
      ENDDO
   ENDDO
ENDDO
 
!!! scale with number of cells
IF (ITYPE == NTYPE_LOCAL .OR. NMESHES==1) THEN                
   SP = Sqrt (SP/REAL(SL%NCELLS_LOCAL, EB))                   ! scale with local number of cells
   IF (SCARC_DEBUG>=3) WRITE(SCARC_LU,*) 'Scale with ',SL%NCELLS_LOCAL,': SP=',SP
ELSE IF (ITYPE == NTYPE_GLOBAL) THEN                          
   IF (USE_MPI) THEN
      CALL MPI_ALLREDUCE (SP, SPG, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
      SP = Sqrt (SPG/REAL(SL%NCELLS_GLOBAL, EB))                 ! scale with global number of cells
      IF (SCARC_DEBUG>=3) WRITE(SCARC_LU,*) 'Scale with ',SL%NCELLS_GLOBAL,': SP=',SP, ': SPG=',SPG
   ELSE
      WRITE(*,*) 'Serial version not yet implemented'
      stop
   ENDIF
ELSE
   WRITE (*,*) 'Wrong type for SCARC_L2NORM3D ', ITYPE
ENDIF

SCARC_L2NORM3D = SP
IF (SCARC_DEBUG>=3) WRITE(SCARC_LU,*) 'SCARC_L2NORM3D=',SP

TUSED_SCARC(17,NM)=TUSED_SCARC(17,NM)+SECOND()-TNOW_L2NORM3D
TUSED_SCARC( 0,NM)=TUSED_SCARC( 0,NM)+SECOND()-TNOW_L2NORM3D
1000 FORMAT('SP=',f25.16,': (I,J,K)=',3i3,': X=',f25.16)
END FUNCTION SCARC_L2NORM3D
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Jacobi preconditioner
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PRECON_JACOBI2D (Y, ILEVEL, NM)
REAL (EB), POINTER, DIMENSION (:, :, :) :: Y
REAL (EB) :: TNOW_JACOBI2D, Y_OLD
INTEGER :: I, K, IC, ILEVEL, NM
 
SL => S%SLEVEL(ILEVEL)

TNOW_JACOBI2D = SECOND()
 
DO K = 1, SL%KBAR
   DO I = 1, SL%IBAR
      IC = (K-1) * SL%IBAR + I
      Y_OLD=Y(I,1,K)
      Y (I, 1, K) = Y (I, 1, K) / SL%AG (IC, 1)
      IF (SCARC_DEBUG>=4.AND.I<=4) WRITE(SCARC_LU,1000) I,1,K,Y(I,1,K),Y_OLD,IC,SL%AG(IC,1)
   ENDDO
ENDDO
 
TUSED_SCARC(18,NM)=TUSED_SCARC(18,NM)+SECOND()-TNOW_JACOBI2D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_JACOBI2D
1000 FORMAT('Y(',I3,',',I3,',',I3,')=',F25.16,':  Y_OLD=',F25.16,': AG(',I3,',1)=',F12.6)
END SUBROUTINE SCARC_PRECON_JACOBI2D
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Jacobi preconditioner
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PRECON_JACOBI3D (Y, ILEVEL, NM)
REAL (EB), POINTER, DIMENSION (:, :, :) :: Y
REAL (EB) :: TNOW_JACOBI3D, Y_OLD
INTEGER :: I, J, K, IC, ILEVEL, NM
 
SL => S%SLEVEL(ILEVEL)

TNOW_JACOBI3D = SECOND()
 
CALL SCARC_SHOW_LEVEL (Y, 'SARC', 'YJAC0 ',ILEVEL)
DO K = 1, SL%KBAR
   DO J = 1, SL%JBAR
      DO I = 1, SL%IBAR
         IC = (K-1) * SL%IBAR * SL%JBAR + (J-1) * SL%IBAR + I
         Y_OLD=Y(I,J,K)
         Y (I, J, K) = Y (I, J, K) / SL%AG (IC, 1)
         IF (SCARC_DEBUG>=4.AND.I<=4) WRITE(SCARC_LU,1000) I,J,K,Y(I,J,K),Y_OLD,IC,SL%AG(IC,1)
      ENDDO
   ENDDO
ENDDO
CALL SCARC_SHOW_LEVEL (Y, 'SARC', 'YJAC1 ',ILEVEL)
 
TUSED_SCARC(18,NM)=TUSED_SCARC(18,NM)+SECOND()-TNOW_JACOBI3D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_JACOBI3D
1000 FORMAT('Y(',I3,',',I3,',',I3,')=',F25.16,':  Y_OLD=',F25.16,': AG(',I3,',1)=',F12.6)
END SUBROUTINE SCARC_PRECON_JACOBI3D
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! SSOR preconditioner
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PRECON_SSOR2D (Y, ILEVEL, NM)
REAL (EB), POINTER, DIMENSION (:, :, :) :: Y
REAL (EB) :: AUX, OMEGA=1.5_EB
INTEGER :: I, K, IC, ILEVEL, NM
REAL (EB) :: TNOW_SSOR2D
 
TNOW_SSOR2D = SECOND()

SL => S%SLEVEL(ILEVEL)
 
DO K = 1, SL%KBAR
   DO I = 1, SL%IBAR
      IC = (K-1) * SL%IBAR + I
      AUX =    SL%AG (IC, 2) * Y (I  , 1, K-1)  &
             + SL%AG (IC, 3) * Y (I-1, 1, K  )
      Y (I, 1, K) = (Y(I, 1, K) - AUX*OMEGA) / SL%AG (IC, 1)
   ENDDO
ENDDO
 
DO K = SL%KBAR, 1, - 1
   DO I = SL%IBAR, 1, - 1
      IC = (K-1) * SL%IBAR + I
      AUX =    SL%AG (IC, 5) * Y (I  , 1, K+1) &
             + SL%AG (IC, 4) * Y (I+1, 1, K  )
      Y (I, 1, K) = Y (I, 1, K) - AUX * OMEGA / SL%AG (IC, 1)
   ENDDO
ENDDO
 
TUSED_SCARC(20,NM)=TUSED_SCARC(20,NM)+SECOND()-TNOW_SSOR2D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_SSOR2D
END SUBROUTINE SCARC_PRECON_SSOR2D
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! SSOR preconditioner
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PRECON_SSOR3D (Y, ILEVEL, NM)
REAL (EB), POINTER, DIMENSION (:, :, :) :: Y
REAL (EB) :: AUX, OMEGA=1.5_EB
INTEGER :: I, J, K, IC, ILEVEL, NM
REAL (EB) :: TNOW_SSOR3D
 
TNOW_SSOR3D = SECOND()

SL => S%SLEVEL(ILEVEL)
 
DO K = 1, SL%KBAR
   DO J = 1, SL%JBAR
      DO I = 1, SL%IBAR
         IC = (K-1) * SL%IBAR * SL%JBAR + (J-1) * SL%IBAR + I
         AUX =    SL%AG (IC, 2) * Y (I  , J  , K-1) &
                + SL%AG (IC, 3) * Y (I  , J-1, K  ) &
                + SL%AG (IC, 4) * Y (I-1, J  , K  )
         Y (I, J, K) = (Y(I, J, K) - AUX * OMEGA) / SL%AG (IC, 1)
      ENDDO
   ENDDO
ENDDO
 
DO K = SL%KBAR, 1, - 1
   DO J = SL%JBAR, 1, - 1
      DO I = SL%IBAR, 1, - 1
         IC = (K-1) * SL%IBAR * SL%JBAR + (J-1) * SL%IBAR + I
         AUX =    SL%AG (IC, 7) * Y (I  , J  , K+1) &
                + SL%AG (IC, 6) * Y (I  , J+1, K  ) &
                + SL%AG (IC, 5) * Y (I+1, J  , K  )
         Y (I, J, K) = Y (I, J, K) - AUX * OMEGA / SL%AG (IC, 1)
      ENDDO
   ENDDO
ENDDO
 
TUSED_SCARC(20,NM)=TUSED_SCARC(20,NM)+SECOND()-TNOW_SSOR3D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_SSOR3D
 
END SUBROUTINE SCARC_PRECON_SSOR3D

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! FFT preconditioner in 2D
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PRECON_FFT2D (R, G, ILEVEL, NM)

USE POIS, ONLY: H2CZSS

REAL (EB), POINTER, DIMENSION (:, :, :) :: R, G
INTEGER :: NM, ILEVEL, I, K
TYPE (MESH_TYPE),  POINTER ::  M


SL => S%SLEVEL(ILEVEL)
M  => MESHES(NM)

DO K = 1, SL%KBAR 
   DO I = 1, SL%IBAR 
      SL%FFT(I, 1, K) = R(I, 1, K)
   ENDDO
ENDDO

CALL H2CZSS (M%BXS, M%BXF, M%BZS, M%BZF, &
             SL%IBAR+1, SL%FFT, M%POIS_PTB, M%SAVE1, M%WORK, M%HX)

DO K = 1, SL%KBAR 
   DO I = 1, SL%IBAR 
      G(I, 1, K) = SL%FFT(I, 1, K)
   ENDDO
ENDDO

END SUBROUTINE SCARC_PRECON_FFT2D
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! FFT preconditioner in 3D
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PRECON_FFT3D (R, G, ILEVEL, NM)

USE POIS, ONLY: H3CZSS

REAL (EB), POINTER, DIMENSION (:, :, :) :: R, G
INTEGER :: NM, ILEVEL, I, J, K
TYPE (MESH_TYPE),  POINTER ::  M

SL => S%SLEVEL(ILEVEL)
M  => MESHES(NM)

DO K = 1, SL%KBAR 
   DO J = 1, SL%JBAR 
      DO I = 1, SL%IBAR
         SL%FFT(I, J, K) = R(I, J, K)
      ENDDO
   ENDDO
ENDDO
!SL%FFT(1:SL%IBAR, 1:SL%JBAR, 1:SL%KBAR) = R(1:SL%IBAR, 1:SL%JBAR, 1:SL%KBAR)

CALL H3CZSS (M%BXS, M%BXF, M%BYS, M%BYF, M%BZS, M%BZF, &
             SL%IBAR+1, SL%JBAR+1, SL%FFT, M%POIS_PTB, M%SAVE1, M%WORK, M%HX)

DO K = 1, SL%KBAR 
   DO J = 1, SL%JBAR 
      DO I = 1, SL%IBAR 
         G(I, J, K) = SL%FFT(I, J, K)
      ENDDO
   ENDDO
ENDDO
!G(1:SL%IBAR, 1:SL%JBAR, 1:SL%KBAR) = SL%FFT(1:SL%IBAR, 1:SL%JBAR, 1:SL%KBAR)

END SUBROUTINE SCARC_PRECON_FFT3D

 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Initialize GSTRIX preconditioner
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PRECON_INIT_GSTRIX2D(ILEVEL, NM)

INTEGER :: NM,ILEVEL, IC, NC, IERR
REAL(EB):: TNOW_GSTRIX2D_INIT
 
TNOW_GSTRIX2D_INIT = SECOND()

IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,*) 'INITIALIZING GSTRIX'
 
IERR=0

SL => S%SLEVEL(ILEVEL)
NC =  SL%NCELLS_LOCAL

ALLOCATE (SL%MDX(1:NC), STAT=IERR)
CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX2D', 'SL%MDX', IERR)

ALLOCATE (SL%UDX(1:NC), STAT=IERR)
CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX2D', 'SL%UDX', IERR)

ALLOCATE (SL%LDX(1:NC), STAT=IERR)
CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX2D', 'SL%LDX', IERR)

ALLOCATE (SL%MDZ(1:NC), STAT=IERR)
CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX2D', 'SL%MDZ', IERR)

ALLOCATE (SL%DAUX(1:NC), STAT=IERR)
CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX2D', 'SL%DAUX', IERR)


 
!!! save different diagonals of matrix ('D' corresponds to main, 'L' to lower, 'U' to upper)
DO IC = 1, NC
   SL%MDX(IC) = SL%AG(IC, 1)         ! main  diagonal in x-direction
   SL%MDZ(IC) = SL%AG(IC, 2)         ! main  diagonal in z-direction
   SL%LDX(IC) = SL%AG(IC, 3)         ! lower diagonal in x-direction
   SL%UDX(IC) = SL%AG(IC, 4)         ! upper diagonal in x-direction
ENDDO
 
!!! perform LU-factorization of matrix AG according to bandwise storage technique
DO IC = 2, NC
   SL%LDX (IC) = SL%LDX(IC) / SL%MDX(IC-1)
   SL%MDX (IC) = SL%MDX(IC) - SL%LDX(IC) * SL%UDX(IC-1)
ENDDO
 
!!! replace diagonal values diag(i) by 1/diag(i) and multiply UDX with 1/diag(i)
DO IC = 1, NC
   SL%MDX (IC) = 1.0_EB / SL%MDX(IC)
   SL%UDX (IC) = SL%MDX(IC) * SL%UDX(IC)
   IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,1000) ILEVEL,IC,SL%UDX(IC),IC,SL%MDX(IC)
ENDDO

IF (SCARC_DEBUG>=2) THEN
   WRITE(SCARC_LU,*) 'MDX(',ILEVEL,'):'
   WRITE(SCARC_LU,'(8f12.6)') (SL%MDX(IC),IC=1,NC)
   WRITE(SCARC_LU,*) 'LDX(',ILEVEL,'):'
   WRITE(SCARC_LU,'(8f12.6)') (SL%LDX(IC),IC=1,NC)
   WRITE(SCARC_LU,*) 'UDX(',ILEVEL,'):'
   WRITE(SCARC_LU,'(8f12.6)') (SL%UDX(IC),IC=1,NC)
   WRITE(SCARC_LU,*) 'MDZ(',ILEVEL,'):'
   WRITE(SCARC_LU,'(8f12.6)') (SL%MDZ(IC),IC=1,NC)
ENDIF

TUSED_SCARC(32,NM)=TUSED_SCARC(32,NM)+SECOND()-TNOW_GSTRIX2D_INIT
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_GSTRIX2D_INIT
1000 FORMAT('GSTRIX2D_INIT:',i3,': SL%UDX(',i3,')=',f12.6,': SL%MDX(',i3,')=',f12.6)
END SUBROUTINE SCARC_PRECON_INIT_GSTRIX2D
 
 
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! GSTRIX2D preconditioner
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PRECON_GSTRIX2D (Y, ILEVEL, NM)
REAL (EB), POINTER, DIMENSION (:, :, :) :: Y
INTEGER :: I, K, IC, ILEVEL, NM
REAL (EB) :: TNOW_GSTRIX2D
 
TNOW_GSTRIX2D = SECOND()

SL => S%SLEVEL(ILEVEL)
 
IF (SCARC_DEBUG>=2) THEN
   WRITE(SCARC_LU,*) '===================== STARTING PRECON_GSTRIX2D'
   WRITE(SCARC_LU,*) 'BEFORE PRECON_GSTRIX: Y'
   CALL SCARC_SHOW(Y, 'SARC', 'Y     ')
ENDIF


! backward elimination of first SL%IBAR unkowns (may be solved by tridiagonal system)
DO I=2,SL%IBAR
   Y(I,1,1) = Y(I,1,1)-SL%LDX(I-1)*Y(I-1,1,1)
ENDDO
DO I=1,SL%IBAR
   Y(I,1,1) = Y(I,1,1)*SL%MDX(I)
ENDDO
DO I=SL%IBAR-1,1,-1
   Y(I,1,1) = Y(I,1,1)-SL%UDX(I)*Y(I+1,1,1)
ENDDO
  

! backward elimination of following unknowns (here the subdiagonal in z-direction must be taken into account)
DO K=2,SL%KBAR

   IC = (K-1)*SL%IBAR + I

   ! eliminate subdiagonal in z-direction into right hand side (related solution component already known)
   DO I=1,SL%IBAR
      Y(I,1,K) = Y(I,1,K) - SL%MDZ(IC)*Y(I,1,K-1) 
   ENDDO

   ! perform elimination of matrix lines corresponding to K
   DO I=2,SL%IBAR
      Y(I,1,K) = Y(I,1,K)-SL%LDX(IC-1)*Y(I-1,1,K)
   ENDDO
   DO I=1,SL%IBAR
      Y(I,1,K) = Y(I,1,K)*SL%MDX(IC)
   ENDDO
   DO I=SL%IBAR-1,1,-1
      Y(I,1,K) = Y(I,1,K)-SL%UDX(IC)*Y(I+1,1,K)
   ENDDO
  
ENDDO

IF (SCARC_DEBUG >= 2) THEN
   WRITE(SCARC_LU,*) 'AFTER  PRECON_GSTRIX: Y'
   CALL SCARC_SHOW(Y, 'SARC',  'Y     ')
ENDIF
 
TUSED_SCARC(31,NM)=TUSED_SCARC(20,NM)+SECOND()-TNOW_GSTRIX2D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_GSTRIX2D
 
END SUBROUTINE SCARC_PRECON_GSTRIX2D


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Initialize GSTRIX preconditioner
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PRECON_INIT_GSTRIX3D(ILEVEL, NM)

INTEGER :: NM,ILEVEL, IC, NC, IERR
REAL(EB):: TNOW_GSTRIX3D_INIT
 
TNOW_GSTRIX3D_INIT = SECOND()
 
IERR=0

SL => S%SLEVEL(ILEVEL)
NC =  SL%NCELLS_LOCAL

ALLOCATE (SL%MDX(1:NC), STAT=IERR)
CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX3D', 'SL%MDX', IERR)

ALLOCATE (SL%UDX(1:NC), STAT=IERR)
CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX3D', 'SL%UDX', IERR)

ALLOCATE (SL%LDX(1:NC), STAT=IERR)
CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX3D', 'SL%LDX', IERR)

ALLOCATE (SL%MDY(1:NC), STAT=IERR)
CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX3D', 'SL%MDY', IERR)

ALLOCATE (SL%MDZ(1:NC), STAT=IERR)
CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX3D', 'SL%MDZ', IERR)

ALLOCATE (SL%DAUX(1:NC), STAT=IERR)
CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX3D', 'SL%DAUX', IERR)


!!! save different diagonals of matrix ('D' corresponds to main, 'L' to lower, 'U' to upper)
DO IC = 1, NC
   SL%MDX(IC) = SL%AG(IC, 1)         ! main  diagonal in x-direction
   SL%MDZ(IC) = SL%AG(IC, 2)         ! main  diagonal in z-direction
   SL%MDY(IC) = SL%AG(IC, 3)         ! main  diagonal in y-direction
   SL%LDX(IC) = SL%AG(IC, 4)         ! lower diagonal in x-direction
   SL%UDX(IC) = SL%AG(IC, 5)         ! upper diagonal in x-direction
ENDDO
 
!!! perform LU-factorization of matrix AG according to bandwise storage technique
DO IC = 2, NC
   SL%LDX (IC) = SL%LDX(IC) / SL%MDX(IC-1)
   SL%MDX (IC) = SL%MDX(IC) - SL%LDX(IC) * SL%UDX(IC-1)
ENDDO
 
!!! replace diagonal values diag(i) by 1/diag(i) and multiply UDX with 1/diag(i)
DO IC = 1, NC
   SL%MDX (IC) = 1.0_EB / SL%MDX(IC)
   SL%UDX (IC) = SL%MDX(IC) * SL%UDX(IC)
   IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,1000) ILEVEL,IC,SL%UDX(IC),IC,SL%MDX(IC)
ENDDO

TUSED_SCARC(32,NM)=TUSED_SCARC(32,NM)+SECOND()-TNOW_GSTRIX3D_INIT
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_GSTRIX3D_INIT
1000 FORMAT('GSTRIX3D_INIT:',i3,': SL%UDX(',i3,')=',f12.6,': SL%MDX(',i3,')=',f12.6)
END SUBROUTINE SCARC_PRECON_INIT_GSTRIX3D
 
 
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! GSTRIX3D preconditioner
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PRECON_GSTRIX3D (Y, ILEVEL, NM)
REAL (EB), POINTER, DIMENSION (:, :, :) :: Y
INTEGER :: I, J, K, IC, ILEVEL, NM
REAL (EB) :: TNOW_GSTRIX3D
 
TNOW_GSTRIX3D = SECOND()

SL => S%SLEVEL(ILEVEL)
 
IF (SCARC_DEBUG>=2) THEN
   WRITE(SCARC_LU,*) '===================== STARTING PRECON_GSTRIX3D'
   WRITE(SCARC_LU,*) 'BEFORE PRECON_GSTRIX: Y'
   CALL SCARC_SHOW(Y, 'SARC', 'Y     ')
ENDIF


!!! NOT WORKING YET, has still to be adopted to 3D !!!!!
DO J=1,SL%JBAR

   ! backward elimination of first SL%IBAR unkowns (may be solved by tridiagonal system)
   DO I=2,SL%IBAR
      Y(I,J,1) = Y(I,J,1)-SL%LDX(I-1)*Y(I-1,J,1)
   ENDDO
   DO I=1,SL%IBAR
      Y(I,J,1) = Y(I,J,1)*SL%MDX(I)
   ENDDO
   DO I=SL%IBAR-1,1,-1
      Y(I,J,1) = Y(I,J,1)-SL%UDX(I)*Y(I+1,J,1)
   ENDDO
     
   
   ! backward elimination of following unknowns 
   ! (here the subdiagonals in y- and z-direction must be taken into account)
   DO K=2,SL%KBAR
   
      IC = (K-1)*SL%IBAR + I
   
      ! bring subdiagonal in z-direction to right hand side (corresponding solution component already known)
      DO I=1,SL%IBAR
         Y(I,J,K) = Y(I,J,K) - SL%MDZ(IC)*Y(I,J,K-1) 
      ENDDO
   
      ! perform elimination of matrix lines corresponding to K
      DO I=2,SL%IBAR
         Y(I,J,K) = Y(I,J,K)-SL%LDX(IC-1)*Y(I-1,J,K)
      ENDDO
      DO I=1,SL%IBAR
         Y(I,J,K) = Y(I,J,K)*SL%MDX(IC)
      ENDDO
      DO I=SL%IBAR-1,1,-1
         Y(I,J,K) = Y(I,J,K)-SL%UDX(IC)*Y(I+1,J,K)
      ENDDO
     
   ENDDO

ENDDO

IF (SCARC_DEBUG >= 2) THEN
   WRITE(SCARC_LU,*) 'AFTER  PRECON_GSTRIX: Y'
   CALL SCARC_SHOW(Y, 'SARC',  'Y     ')
ENDIF
 
TUSED_SCARC(31,NM)=TUSED_SCARC(20,NM)+SECOND()-TNOW_GSTRIX3D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_GSTRIX3D
 
END SUBROUTINE SCARC_PRECON_GSTRIX3D




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Perform global MG-method as preconditioner in 2D
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PRECON_MG2D (Y,ILEVEL,NM)

REAL (EB), POINTER, DIMENSION (:, :, :) :: Y
INTEGER :: NM, ITE, ITE0, ICYCLE, ILEVEL
INTEGER :: I, K, IBAR0, KBAR0
REAL (EB) :: TNOW_PRECON_MG2D
LOGICAL BMATVEC, BCONV, BDIVG
!TYPE (SCARC_LEVEL_TYPE), POINTER :: SL, SLHI, SLLO
TYPE (MESH_TYPE),  POINTER :: M
 
TNOW_PRECON_MG2D = SECOND()

! initialize working vectors
S    => SCARC (NM)
M    => MESHES (NM)

SLHI  => S%SLEVEL(ILEVEL)
 
! initialize working vectors
IBAR0=SLHI%IBAR
KBAR0=SLHI%KBAR

SLHI%X2 = 0.0_EB
 
DO K = 1, KBAR0                                ! most possible redundant --- to be checked (use Y instead)
   DO I = 1, IBAR0
      SLHI%F2 (I, 1, K) = Y (I, 1, K)
   ENDDO
ENDDO
 
IF (SCARC_DEBUG >= 3) THEN
   WRITE(SCARC_LU,*) 'Starting SCARC_PRECON_MG2D'
   WRITE(SCARC_LU,*) 'IBAR=:', SLHI%IBAR
   WRITE(SCARC_LU,*) 'JBAR=:', SLHI%JBAR
   WRITE(SCARC_LU,*) 'KBAR=:', SLHI%KBAR
   WRITE(SCARC_LU,*) 'SCARC_MG_NIT=:', SCARC_MG_NIT
   WRITE(SCARC_LU,*) 'SCARC_MG_EPS=:', SCARC_MG_EPS
   WRITE(SCARC_LU,*) 'SCARC_EPS_REL=:', SCARC_EPS_REL
   CALL SCARC_SHOW(SLHI%X2, 'SARC',  'X2    ')
   CALL SCARC_SHOW(SLHI%F2, 'SARC',  'F2    ')
ENDIF
 
!!!
!!! initialize some method parameters
!!!
NREQ_FACE = 0           ! necessary ?

BCONV =.FALSE.
BDIVG =.FALSE.
ICYCLE=1      ! V-cycle
 

!!! save cycle counts for MG-iteration
S%KCYCLE(2,S%NLMAX)=1
DO ILEVEL=S%NLMIN+1,S%NLMAX-1
   IF (ICYCLE==0) THEN
      S%KCYCLE(2,ILEVEL)=2
   ELSE
      S%KCYCLE(2,ILEVEL)=ICYCLE
   ENDIF
ENDDO
  

!!! calculate initial defect on finest level, Y = B - A*X, and get l2-norm of it
!WRITE(SCARC_LU,*) 'KOMMUNIKATIONSPARAMETER 3 prüfen, anders als bei CG !!!'
CALL SCARC_MATVEC2D (SLHI%X2, SLHI%D2, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, NMV_X2, NMV_D2)
SLHI%D2 = SLHI%F2 - SLHI%D2
SLHI%RESIN_MG = SCARC_L2NORM2D (SLHI%D2, ILEVEL, NM, NTYPE_GLOBAL)
IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) 0, ILEVEL, SL%RESIN_MG
 
!WRITE(*,*) 'Info muss noch an Master geschickt werden'
!CALL SCARC_MASTER_INFO (SLHI%ITE, SLHI%RESIN_MG)

 
!!! ---------------------------------------------------------------------------
!!! start MG-iteration
!!! ---------------------------------------------------------------------------
MG_LOOP2D: DO ITE = 1, SCARC_MG_NIT
 
   !!! set level-information
   DO ILEVEL=S%NLMIN,S%NLMAX
      S%KCYCLE(1,ILEVEL)=S%KCYCLE(2,ILEVEL)
   ENDDO

   ILEVEL=S%NLMAX

   !!! ---------------------------------------------------------------------------
   !!!  Presmoothing
   !!! ---------------------------------------------------------------------------
110   IF (ILEVEL/=S%NLMIN) THEN
 
      SLHI => S%SLEVEL(ILEVEL)
      SLLO => S%SLEVEL(ILEVEL-1)

      !IF (SCARC_METHOD>=2) THEN
      !   WRITE(SCARC_LU,*) 'SLHI zeigt auf level ',ILEVEL
      !   WRITE(SCARC_LU,*) 'SLLO zeigt auf level ',ILEVEL-1
      !ENDIF

      IF (ILEVEL==S%NLMAX) THEN
         BMATVEC=.FALSE.
      ELSE
         BMATVEC=.TRUE.
      ENDIF

      CALL SCARC_PRECON_SMOOTHER2D(SLHI%X2, SLHI%F2, SLHI%D2, BMATVEC, ILEVEL, NM)

      !!! print level information - still missing
      !CALL SCARC_INFO_LEVEL(ILEVEL)

      !!! perform restriction to coarser grid: F:= rest(D)
      CALL SCARC_SHOW_LEVEL (SLHI%D2, 'SARC',  'RESTD ',ILEVEL)
      CALL SCARC_RESTRICTION2D(SLLO%F2, SLHI%D2, ILEVEL-1, NM)
      CALL SCARC_SHOW_LEVEL (SLLO%F2, 'SARC',  'RESTF ',ILEVEL-1)

      !!! initialize solution vector and set boundary conditions to residuum - necessary ?
      SLLO%X2=0.0_EB
      !CALL SCARC_BDRY_RESIDUUM(NM, SLLO%F, ILEVEL-1)

      !!! decrease level
      ILEVEL=ILEVEL-1

      GOTO 110
   
   ENDIF         

   !!! ------------------------------------------------------------------------
   !!! Coarse grid solver
   !!! ------------------------------------------------------------------------
   ILEVEL=S%NLMIN
   CALL SCARC_PRECON_COARSE2D(ILEVEL, NM)


   !!! ------------------------------------------------------------------------
   !!! Postsmoothing
   !!! ------------------------------------------------------------------------
120   IF (ILEVEL/=S%NLMAX) THEN
 
      BMATVEC=.TRUE.

      SLLO => S%SLEVEL(ILEVEL)

      ILEVEL=ILEVEL+1
      SLHI => S%SLEVEL(ILEVEL)

      !!! print level information - still missing
      !CALL SCARC_INFO_LEVEL(ILEVEL)

      !!! perform prolongation to finer grid: D:=prol(X)
      CALL SCARC_PROLONGATION2D(SLLO%X2, SLHI%D2, ILEVEL, NM)

      !!! set exterior boundary data of residuum to zero - necessary ?
      !CALL SCARC_BDRY_RESIDUUM(NM, SLHI%D, ILEVEL)
 
      !!! set new solution
      SLHI%X2 = SLHI%D2 + SLHI%X2

      !!! select smoother for postsmoothing
      CALL SCARC_PRECON_SMOOTHER2D(SLHI%X2, SLHI%F2, SLHI%D2, BMATVEC, ILEVEL, NM)


      !!! save cycle counts
      S%KCYCLE(1,ILEVEL)=S%KCYCLE(1,ILEVEL)-1
 
      IF (S%KCYCLE(1,ILEVEL)==0) THEN
         IF (ICYCLE==0) THEN
            S%KCYCLE(1,ILEVEL)=1
         ELSE
            S%KCYCLE(1,ILEVEL)=S%KCYCLE(2,ILEVEL)
         ENDIF
         GOTO 120
      ELSE
         GOTO 110
      ENDIF

   ENDIF

 
   !!! ------------------------------------------------------------------------
   !!! Defect calculation
   !!! ------------------------------------------------------------------------
   CALL SCARC_MATVEC2D (SLHI%X2, SLHI%D2, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, NMV_X2, NMV_D2)
   SLHI%D2 = SLHI%F2 - SLHI%D2
   SLHI%RES_MG = SCARC_L2NORM2D (SLHI%D2, ILEVEL, NM, NTYPE_GLOBAL)

   IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,*) 'SCARC-Multigrid, iteration ',ITE,': residuum=',SLHI%RES_MG

   !!! ------------------------------------------------------------------------
   !!! Send error information to Master  - still missing
   !!! ------------------------------------------------------------------------
   !CALL SCARC_INFO_MASTER 
   !CALL SCARC_INFO_LEVEL(ILEVEL)

   !!! ------------------------------------------------------------------------
   !!! Convergence or divergence ?
   !!! ------------------------------------------------------------------------
   IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) ITE, ILEVEL, SL%RES_MG
   IF (SCARC_DEBUG >= 1.AND.NM==1) write (*,1000) ITE, ILEVEL, SL%RES_MG
   IF (SCARC_BREL) THEN
      IF (SL%RES_MG <= SL%RESIN_MG*SCARC_MG_EPS) BCONV = .TRUE.
   ELSE
      IF (SL%RES_MG <= SCARC_MG_EPS .AND. SL%RES_MG <= SL%RESIN_MG*SCARC_EPS_REL) BCONV = .TRUE.
   ENDIF
   IF (SL%RES_MG > SCARC_EPS_DIVG) BDIVG = .TRUE.
   IF (BCONV.OR.BDIVG) EXIT MG_LOOP2D
 
ENDDO MG_LOOP2D

!!! ------------------------------------------------------------------------
!!! Convergence or divergence ?
!!! ------------------------------------------------------------------------
IF (BDIVG) THEN                       
   ITE0 = - 1
   SL%CAPPA_MG = 1.0_EB
ELSE 
   IF (SL%RESIN_MG >= 1.0E-70_EB) THEN
      SL%RESIN_MG = SL%RES_MG / SL%RESIN_MG
      SL%CAPPA_MG = SL%RESIN_MG ** (1.0_EB/ITE)
   ELSE 
      SL%CAPPA_MG = 0.0_EB
   ENDIF
   IF (BCONV) THEN
      ITE0=ITE
   ELSE
      ITE0=ITE-1
   ENDIF
ENDIF
IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,2000) ITE0, ILEVEL, SL%RES_MG, SL%CAPPA_MG 
IF (SCARC_DEBUG >= 1.AND.NM==1) WRITE(*,2000) ITE0, ILEVEL, SL%RES_MG, SL%CAPPA_MG 
 
DO K = 1, SLHI%KBAR
   DO I = 1, SLHI%IBAR
      Y (I, 1, K) = SLHI%X2(I, 1, K)
   ENDDO
ENDDO
 
TUSED_SCARC(7,NM)=TUSED_SCARC(7,NM)+SECOND()-TNOW_PRECON_MG2D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_PRECON_MG2D

1000 FORMAT (/,'     SCARC_MG_PRECON: #ite= ',i4,': level=',i4,': res=',e16.8,/)
2000 FORMAT (/,'     SCARC_MG_PRECON: #ite= ',i4,': level=',i4,': res=',e16.8,': rate=',e16.8,/)
END SUBROUTINE SCARC_PRECON_MG2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Perform global MG-method as preconditioner in 3D
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PRECON_MG3D (Y,ILEVEL,NM)

REAL (EB), POINTER, DIMENSION (:, :, :) :: Y
INTEGER :: NM, ITE, ITE0, ICYCLE, ILEVEL
INTEGER :: I, J, K, IBAR0, JBAR0, KBAR0
REAL (EB) :: TNOW_PRECON_MG3D
LOGICAL BMATVEC, BCONV, BDIVG
!TYPE (SCARC_LEVEL_TYPE), POINTER :: SL, SLHI, SLLO
TYPE (MESH_TYPE),  POINTER :: M
 
TNOW_PRECON_MG3D = SECOND()

! initialize working vectors
S    => SCARC (NM)
M    => MESHES (NM)

SLHI  => S%SLEVEL(ILEVEL)
 
! initialize working vectors
SLHI%X2 = 0.0_EB

IBAR0=SLHI%IBAR
JBAR0=SLHI%JBAR
KBAR0=SLHI%KBAR
 
DO K = 1, KBAR0                                ! most possibly redundant --- to be checked (use Y instead)
   DO J = 1, JBAR0
      DO I = 1, IBAR0
         SLHI%F2 (I, J, K) = Y (I, J, K)
      ENDDO
   ENDDO
ENDDO
 
IF (SCARC_DEBUG >= 3) THEN
   WRITE(SCARC_LU,*) 'Starting SCARC_PRECON_MG3D'
   WRITE(SCARC_LU,*) 'IBAR=:', SLHI%IBAR
   WRITE(SCARC_LU,*) 'JBAR=:', SLHI%JBAR
   WRITE(SCARC_LU,*) 'KBAR=:', SLHI%KBAR
   WRITE(SCARC_LU,*) 'SCARC_MG_NIT=:', SCARC_MG_NIT
   WRITE(SCARC_LU,*) 'SCARC_MG_EPS=:', SCARC_MG_EPS
   WRITE(SCARC_LU,*) 'SCARC_EPS_REL=:', SCARC_EPS_REL
   CALL SCARC_SHOW(SLHI%X2, 'SARC',  'X2    ')
   CALL SCARC_SHOW(SLHI%F2, 'SARC',  'F2    ')
ENDIF
 
!!!
!!! initialize some method parameters
!!!
NREQ_FACE = 0           ! necessary ?

BCONV =.FALSE.
BDIVG =.FALSE.
ICYCLE=1      ! V-cycle
 
!!! save cycle counts for MG-iteration
S%KCYCLE(2,S%NLMAX)=1
DO ILEVEL=S%NLMIN+1,S%NLMAX-1
   IF (ICYCLE==0) THEN
      S%KCYCLE(2,ILEVEL)=2
   ELSE
      S%KCYCLE(2,ILEVEL)=ICYCLE
   ENDIF
ENDDO
  


!!! calculate initial defect on finest level, Y = B - A*X, and get l2-norm of it
!WRITE(SCARC_LU,*) 'KOMMUNIKATIONSPARAMETER 3 prüfen, anders als bei CG !!!'
CALL SCARC_MATVEC3D (SLHI%X2, SLHI%D2, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, NMV_X2, NMV_D2)
SLHI%D2 = SLHI%F2 - SLHI%D2
CALL SCARC_SHOW(SLHI%D2, 'SARC',  'D2    ')
SLHI%RESIN_MG = SCARC_L2NORM3D (SLHI%D2, ILEVEL, NM, NTYPE_GLOBAL)
IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) 0, ILEVEL, SL%RESIN_MG
 
!WRITE(*,*) 'Info muss noch an Master geschickt werden'
!CALL SCARC_MASTER_INFO (SLHI%ITE, SLHI%RESIN_MG)

 
!!! ---------------------------------------------------------------------------
!!! start MG-iteration
!!! ---------------------------------------------------------------------------
MG_LOOP3D: DO ITE = 1, SCARC_MG_NIT
 
IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,*) '============ HALLO, MG-Precon-Iteration ',ITE
   !!! set level-information
   DO ILEVEL=S%NLMIN,S%NLMAX
      S%KCYCLE(1,ILEVEL)=S%KCYCLE(2,ILEVEL)
   ENDDO

   ILEVEL=S%NLMAX

   !!! ---------------------------------------------------------------------------
   !!!  Presmoothing
   !!! ---------------------------------------------------------------------------
110   IF (ILEVEL/=S%NLMIN) THEN
 
      SLHI => S%SLEVEL(ILEVEL)
      SLLO => S%SLEVEL(ILEVEL-1)

      !IF (SCARC_METHOD>=2) THEN
      !   WRITE(SCARC_LU,*) 'SLHI zeigt auf level ',ILEVEL
      !   WRITE(SCARC_LU,*) 'SLLO zeigt auf level ',ILEVEL-1
      !ENDIF

      IF (ILEVEL==S%NLMAX) THEN
         BMATVEC=.FALSE.
      ELSE
         BMATVEC=.TRUE.
      ENDIF

      CALL SCARC_PRECON_SMOOTHER3D(SLHI%X2, SLHI%F2, SLHI%D2, BMATVEC, ILEVEL, NM)

      !!! print level information - still missing
      !CALL SCARC_INFO_LEVEL(ILEVEL)

      !!! perform restriction to coarser grid: F:= rest(D)
      CALL SCARC_SHOW_LEVEL (SLHI%D2, 'SARC',  'RESTD ',ILEVEL)
      CALL SCARC_RESTRICTION3D(SLLO%F2, SLHI%D2, ILEVEL-1, NM)
      CALL SCARC_SHOW_LEVEL (SLLO%F2, 'SARC',  'RESTF ',ILEVEL-1)

      !!! initialize solution vector and set boundary conditions to residuum - necessary ?
      SLLO%X2=0.0_EB
      !CALL SCARC_BDRY_RESIDUUM(NM, SLLO%F, ILEVEL-1)

      !!! decrease level
      ILEVEL=ILEVEL-1

      GOTO 110
   
   ENDIF         

   !!! ------------------------------------------------------------------------
   !!! Coarse grid solver
   !!! ------------------------------------------------------------------------
   ILEVEL=S%NLMIN
   CALL SCARC_PRECON_COARSE3D(ILEVEL, NM)


   !!! ------------------------------------------------------------------------
   !!! Postsmoothing
   !!! ------------------------------------------------------------------------
120   IF (ILEVEL/=S%NLMAX) THEN
 
      BMATVEC=.TRUE.

      SLLO => S%SLEVEL(ILEVEL)

      ILEVEL=ILEVEL+1
      SLHI => S%SLEVEL(ILEVEL)

      !!! print level information - still missing
      !CALL SCARC_INFO_LEVEL(ILEVEL)

      !!! perform prolongation to finer grid: D:=prol(X)
      CALL SCARC_PROLONGATION3D(SLLO%X2, SLHI%D2, ILEVEL, NM)
      CALL SCARC_SHOW_LEVEL (SLHI%D2, 'SARC',  'PROLD2',ILEVEL-1)

      !!! set exterior boundary data of residuum to zero - necessary ?
      !CALL SCARC_BDRY_RESIDUUM(NM, SLHI%D, ILEVEL)
 
      !!! set new solution
      SLHI%X2 = SLHI%D2 + SLHI%X2
      CALL SCARC_SHOW_LEVEL (SLHI%X2, 'SARC',  'PROLX2',ILEVEL-1)

      !!! select smoother for postsmoothing
      CALL SCARC_PRECON_SMOOTHER3D(SLHI%X2, SLHI%F2, SLHI%D2, BMATVEC, ILEVEL, NM)
      CALL SCARC_SHOW_LEVEL (SLHI%X2, 'SARC',  'SMOX2 ',ILEVEL-1)


      !!! save cycle counts
      S%KCYCLE(1,ILEVEL)=S%KCYCLE(1,ILEVEL)-1
 
      IF (S%KCYCLE(1,ILEVEL)==0) THEN
         IF (ICYCLE==0) THEN
            S%KCYCLE(1,ILEVEL)=1
         ELSE
            S%KCYCLE(1,ILEVEL)=S%KCYCLE(2,ILEVEL)
         ENDIF
         GOTO 120
      ELSE
         GOTO 110
      ENDIF

   ENDIF

 
   !!! ------------------------------------------------------------------------
   !!! Defect calculation
   !!! ------------------------------------------------------------------------
   CALL SCARC_MATVEC3D (SLHI%X2, SLHI%D2, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, NMV_X2, NMV_D2) 
   SLHI%D2 = SLHI%F2 - SLHI%D2
   SLHI%RES_MG = SCARC_L2NORM3D (SLHI%D2, ILEVEL, NM, NTYPE_GLOBAL)

   IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,*) 'SCARC-Multigrid, iteration ',ITE,': residuum=',SLHI%RES_MG

   !!! ------------------------------------------------------------------------
   !!! Send error information to Master  - still missing
   !!! ------------------------------------------------------------------------
   !CALL SCARC_INFO_MASTER 
   !CALL SCARC_INFO_LEVEL(ILEVEL)

   !!! ------------------------------------------------------------------------
   !!! Convergence or divergence ?
   !!! ------------------------------------------------------------------------
   IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) ITE, ILEVEL, SL%RES_MG
   IF (SCARC_DEBUG >= 2.AND.NM==1) write (*,1000) ITE, ILEVEL, SL%RES_MG
   IF (SCARC_BREL) THEN
      IF (SL%RES_MG <= SL%RESIN_MG*SCARC_MG_EPS) BCONV = .TRUE.
   ELSE
      IF (SL%RES_MG <= SCARC_MG_EPS .AND. SL%RES_MG <= SL%RESIN_MG*SCARC_EPS_REL) BCONV = .TRUE.
   ENDIF
   IF (SL%RES_MG > SCARC_EPS_DIVG) BDIVG = .TRUE.
   IF (BCONV.OR.BDIVG) EXIT MG_LOOP3D
 
ENDDO MG_LOOP3D

!!! ------------------------------------------------------------------------
!!! Convergence or divergence ?
!!! ------------------------------------------------------------------------
IF (BDIVG) THEN                       
   ITE0 = - 1
   SL%CAPPA_MG = 1.0_EB
ELSE 
   IF (SL%RESIN_MG >= 1.0E-70_EB) THEN
      SL%RESIN_MG = SL%RES_MG / SL%RESIN_MG
      SL%CAPPA_MG = SL%RESIN_MG ** (1.0_EB/ITE)
   ELSE 
      SL%CAPPA_MG = 0.0_EB
   ENDIF
   IF (BCONV) THEN
      ITE0=ITE
   ELSE
      ITE0=ITE-1
   ENDIF
ENDIF
IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,2000) ITE0, ILEVEL, SL%RES_MG, SL%CAPPA_MG 
IF (SCARC_DEBUG >= 1.AND.NM==1) WRITE(*,2000) ITE0, ILEVEL, SL%RES_MG, SL%CAPPA_MG 
 
DO K = 1, KBAR0
   DO J = 1, JBAR0
      DO I = 1, IBAR0
         Y (I, J, K) = SLHI%X2(I, J, K)
      ENDDO
   ENDDO
ENDDO
 
TUSED_SCARC(7,NM)=TUSED_SCARC(7,NM)+SECOND()-TNOW_PRECON_MG3D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_PRECON_MG3D

1000 FORMAT (/,'     SCARC_MG_PRECON: #ite= ',i4,': level=',i4,': res=',e16.8,/)
2000 FORMAT (/,'     SCARC_MG_PRECON: #ite= ',i4,': level=',i4,': res=',e16.8,': rate=',e16.8,/)
END SUBROUTINE SCARC_PRECON_MG3D


 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! GSADI preconditioner
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!SUBROUTINE SCARC_PRECON_GSADI2D (A, Y, ITE, ILEVEL, NM)
!REAL (EB), POINTER, DIMENSION (:, :) :: A
!REAL (EB), POINTER, DIMENSION (:, :, :) :: Y
!INTEGER :: I, K, IC, ILEVEL, NM, ITE
!REAL (EB) :: TNOW_GSADI2D
! 
!TNOW_GSADI2D = SECOND()
!
!SL => S%SLEVEL(ILEVEL)
! 
!IF (MOD(ITE,2)==1) THEN
!  CALL SCARC_PRECON_GSTRIX2D(A, Y, ILEVEL, NM)
!ELSE
!  CALL SCARC_PRECON_GSTRI2D_Z(A, Y, ILEVEL, NM)
!ENDIF
!
!TUSED_SCARC(31,NM)=TUSED_SCARC(20,NM)+SECOND()-TNOW_GSADI2D
!TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_GSADI2D
! 
!END SUBROUTINE SCARC_PRECON_GSADI2D
 
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  STILL EXPERIMENTAL !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_RECEIVE (CODE, ILEVEL)

INTEGER :: NM, NOM, CODE, ILEVEL
INTEGER :: IERR, II, IW
INTEGER :: ILEN_FACE, IM, JM
INTEGER :: IMIN, IMAX, JMIN, JMAX, KMIN, KMAX
INTEGER :: TAG_FACE
REAL(EB):: TNOW_RECEIVE

TYPE (SCARC_TYPE),  POINTER ::  SNM
TYPE (OSCARC_TYPE), POINTER :: OSNM
TYPE (SCARC_LEVEL_TYPE),  POINTER ::  SNML
TYPE (OSCARC_LEVEL_TYPE), POINTER :: OSNML

TNOW_RECEIVE = SECOND()

IERR=0
RECEIVE_MESH_LOOP: DO NM=1,NMESHES

   IF (PROCESS(NM)/=MYID) CYCLE RECEIVE_MESH_LOOP
 
   RECEIVE_OMESH_LOOP: DO NOM=1,NMESHES
    
      SNODE = PROCESS(NOM)
      IF (PROCESS(NM)==SNODE) CYCLE RECEIVE_OMESH_LOOP

      SNM  => SCARC(NM)                     ! corresponds to M
      SNML => SCARC(NM)%SLEVEL(ILEVEL)      ! ... for the level 'ILEVEL'


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! FACE communication
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      RECEIVE_FACE_IF: IF (SNML%NIC(NM,NOM)/=0 .AND. SNML%NIC(NOM,NM)/=0) THEN
    
         OSNM  => SCARC(NM)%OSCARC(NOM)                  ! corresponds to M3
         OSNML => SCARC(NM)%OSCARC(NOM)%SLEVEL(ILEVEL)   ! ... for the level 'ILEVEL'


         IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,*) '=========== RECEIVE: NM=',NM,', NOM=',NOM
         IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,*) 'OSNML zeigt auf SCARC(',NM,')%OSCARC(',NOM,')%SLEVEL(',ILEVEL,')'

         TAG_FACE = TAGS_FACE(NM,NOM)

         IMIN=I_MIN(NOM,NM)
         IMAX=I_MAX(NOM,NM)
         JMIN=J_MIN(NOM,NM)
         JMAX=J_MAX(NOM,NM)
         KMIN=K_MIN(NOM,NM)
         KMAX=K_MAX(NOM,NM)
         ILEN_FACE=(IMAX-IMIN+1)*(JMAX-JMIN+1)*(KMAX-KMIN+1)

IF (SCARC_DEBUG>=6) THEN
   WRITE(SCARC_LU,*) 'SCARC_EXCHANGE:   I_MIN:'
   WRITE(SCARC_LU,'(2i4)') ((SNML%I_MIN(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
   WRITE(SCARC_LU,*) 'SCARC_EXCHANGE:   I_MAX:'
   WRITE(SCARC_LU,'(2i4)') ((SNML%I_MAX(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
   WRITE(SCARC_LU,*) 'SCARC_EXCHANGE:   J_MIN:'
   WRITE(SCARC_LU,'(2i4)') ((SNML%J_MIN(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
   WRITE(SCARC_LU,*) 'SCARC_EXCHANGE:   J_MAX:'
   WRITE(SCARC_LU,'(2i4)') ((SNML%J_MAX(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
   WRITE(SCARC_LU,*) 'SCARC_EXCHANGE:   K_MIN:'
   WRITE(SCARC_LU,'(2i4)') ((SNML%K_MIN(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
   WRITE(SCARC_LU,*) 'SCARC_EXCHANGE:   K_MAX:'
   WRITE(SCARC_LU,'(2i4)') ((SNML%K_MAX(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
ENDIF

    

         !!! Initialize the communication structures for receiving face data
         INIT_FACE_IF: IF (CODE==NCOM_INIT) THEN

            IF (ILEVEL/=SNM%NLMAX) THEN
               NREQ_FACE = NREQ_FACE+1
               IF (USE_MPI) CALL MPI_IRECV(OSNML%IJKW(1,1),15*OSNML%N_EXTERNAL_WALL_CELLS,MPI_INTEGER,SNODE, &
                                           TAG_FACE,MPI_COMM_WORLD,REQ_FACE(NREQ_FACE),IERR)
            ELSE
               OSNML%IJKW => MESHES(NM)%OMESH(NOM)%IJKW
            ENDIF

            IF (SCARC_DEBUG>=6) THEN
               WRITE(SCARC_LU,*) 'RECEIVE: INIT: LEVEL=',ILEVEL
               WRITE(SCARC_LU,*) 'RECEIVE: INIT: N_EXTERNAL_WALL_CELLS=',OSNML%N_EXTERNAL_WALL_CELLS
               WRITE(SCARC_LU,*) 'RECEIVE: INIT: LEN(IJKW)=',15*OSNML%N_EXTERNAL_WALL_CELLS
               WRITE(SCARC_LU,'(15i4)') SNML%IJKW
               WRITE(SCARC_LU,*) 'RECEIVE: INIT: SIZE(IJKW):  ',SIZE(OSNML%IJKW)
               WRITE(SCARC_LU,*) 'RECEIVING IJKW(',ILEVEL,'): '
               DO IW= 1,OSNML%N_EXTERNAL_WALL_CELLS
                  WRITE(SCARC_LU,'(15i4)') (MESHES(NM)%OMESH(NOM)%IJKW(II,IW),II=1,15)
               ENDDO
            ENDIF

            IF (SNML%NIC(NM, NOM) > 0) THEN
               !ILEN_FACE=(MAX(SNML%NIC(NM, NOM), SNML%NIC(NOM, NM))+2)*2+1
               !ILEN_FACE=(MAX(SNML%NIC(NM, NOM), SNML%NIC(NOM, NM)))*2+1
               IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,*) 'RECEIVE: ILEN_FACE=',ILEN_FACE
               ALLOCATE (OSNML%RECV_FACE(ILEN_FACE))
               OSNML%RECV_FACE = 0.0_EB
            ENDIF

         ENDIF INIT_FACE_IF
   
         !!! Perform exchange for matrix-vector multiplication (only face exchange needed)
         MATV_FACE_IF: IF (CODE==NCOM_MATV) THEN

            NREQ_FACE = NREQ_FACE+1
            IF (USE_MPI) CALL MPI_IRECV(OSNML%RECV_FACE(1),SIZE(OSNML%RECV_FACE),MPI_DOUBLE_PRECISION,&
                                        SNODE,TAG_FACE,MPI_COMM_WORLD,REQ_FACE(NREQ_FACE),IERR)

            IF (SCARC_DEBUG>=6) THEN
               WRITE(SCARC_LU,*) 'RECEIVE: MATV: receive length',SIZE(OSNML%RECV_FACE)
               WRITE(SCARC_LU,*) 'RECEIVE: MATV: SNODE=',SNODE 
               WRITE(SCARC_LU,*) 'RECEIVE: MATV: TAG_FACE=',TAG_FACE 
               WRITE(SCARC_LU,*) 'RECEIVE: MATV: NREQ_FACE=',NREQ_FACE 
               WRITE(SCARC_LU,*) 'RECEIVE: MATV: REQ_FACE =',REQ_FACE(NREQ_FACE) 
               WRITE(SCARC_LU,*) 'RECEIVE: MATV: RECV_FACE(2)'
            ENDIF

         ENDIF MATV_FACE_IF
   
         !!! Perform full exchange including edge (3D)
         FULL_FACE_IF: IF (CODE==NCOM_FULL) THEN

            NREQ_FACE = NREQ_FACE+1
            IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,*) 'RECEIVE:', SIZE(OSNML%RECV_FACE)
            IF (USE_MPI) CALL MPI_IRECV(OSNML%RECV_FACE(1),SIZE(OSNML%RECV_FACE),MPI_DOUBLE_PRECISION,&
                                        SNODE,TAG_FACE,MPI_COMM_WORLD,REQ_FACE(NREQ_FACE),IERR)

         ENDIF FULL_FACE_IF

      ENDIF RECEIVE_FACE_IF


   ENDDO RECEIVE_OMESH_LOOP
ENDDO RECEIVE_MESH_LOOP


TUSED_SCARC(22,MYID+1)=TUSED_SCARC(22,MYID+1)+SECOND()-TNOW_RECEIVE
TUSED_SCARC(0,MYID+1) =TUSED_SCARC(0,MYID+1) +SECOND()-TNOW_RECEIVE
END SUBROUTINE SCARC_RECEIVE
 
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  STILL EXPERIMENTAL !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_EXCHANGE (CODE, ILEVEL, IMV1, IMV2)

INTEGER :: NM, NOM, ILEVEL, CODE, IMV1, IMV2, ISUM
INTEGER :: I, J, K, LL, IW, IWW, IOR0
INTEGER :: IERR, II, JJ, KK
INTEGER :: ILEN_FACE, IM, JM, NPTR
INTEGER :: IMIN, IMAX, JMIN, JMAX, KMIN, KMAX
INTEGER :: TAG_FACE

REAL(EB):: ASUB, ZSUM,yold
REAL(EB):: TNOW_EXCHANGE


TYPE (SCARC_TYPE),  POINTER ::  SNM,   SNOM
TYPE (OSCARC_TYPE), POINTER :: OSNM, OSNOM
TYPE (SCARC_LEVEL_TYPE),  POINTER ::  SNML,  SNOML
TYPE (OSCARC_LEVEL_TYPE), POINTER :: OSNML, OSNOML

TNOW_EXCHANGE = SECOND()

 
IERR=0
IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,*) 'SCARC_EXCHANGE: ', IMV1, IMV2, CODE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Sample communication data in corresponding SEND-buffers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
EXCHANGE_SEND_MESH_LOOP: DO NM=1,NMESHES

   IF (PROCESS(NM)/=MYID) CYCLE EXCHANGE_SEND_MESH_LOOP
 
   EXCHANGE_RECV_MESH_LOOP: DO NOM=1,NMESHES
    
      SNODE = PROCESS(NOM)
      RNODE = PROCESS(NM)
   
      SNM  => SCARC(NM)                     ! corresponds to M
      SNML => SCARC(NM)%SLEVEL(ILEVEL)      ! ... for the level 'ILEVEL'

      IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,*) '============= EXCHANGE: NM=',NM,': NOM=',NOM
      IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,*) 'SNML zeigt auf SCARC(',NM,')%SLEVEL(',ILEVEL,')'

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! FACE communication
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FACE_IF: IF (SNML%NIC(NOM,NM)/=0 .AND. SNML%NIC(NM,NOM)/=0)  THEN
    
         OSNM  => SCARC(NM)%OSCARC(NOM)                 ! corresponds to M3
         OSNML => SCARC(NM)%OSCARC(NOM)%SLEVEL(ILEVEL)  ! ... for the level 'ILEVEL'

         IF (SCARC_DEBUG>=6) THEN
            WRITE(SCARC_LU,*) 'OSNML zeigt auf SCARC(',NM,')%OSCARC(',NOM,')%SLEVEL(',ILEVEL,')'
            WRITE(SCARC_LU,*) 'SIZE(IJKW(1)):  ',SIZE(SCARC(NM)%OSCARC(NOM)%SLEVEL(1)%IJKW)
            WRITE(SCARC_LU,*) 'SIZE(IJKW(2)):  ',SIZE(SCARC(NM)%OSCARC(NOM)%SLEVEL(2)%IJKW)
         ENDIF

         TAG_FACE = TAGS_FACE(NM,NOM)
    
         IMIN=I_MIN(NOM,NM)
         IMAX=I_MAX(NOM,NM)
         JMIN=J_MIN(NOM,NM)
         JMAX=J_MAX(NOM,NM)
         KMIN=K_MIN(NOM,NM)
         KMAX=K_MAX(NOM,NM)
         ILEN_FACE=(IMAX-IMIN+1)*(JMAX-JMIN+1)*(KMAX-KMIN+1)

IF (SCARC_DEBUG>=6) THEN
   WRITE(SCARC_LU,*) 'SCARC_EXCHANGE:   I_MIN:'
   WRITE(SCARC_LU,'(2i4)') ((SNML%I_MIN(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
   WRITE(SCARC_LU,*) 'SCARC_EXCHANGE:   I_MAX:'
   WRITE(SCARC_LU,'(2i4)') ((SNML%I_MAX(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
   WRITE(SCARC_LU,*) 'SCARC_EXCHANGE:   J_MIN:'
   WRITE(SCARC_LU,'(2i4)') ((SNML%J_MIN(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
   WRITE(SCARC_LU,*) 'SCARC_EXCHANGE:   J_MAX:'
   WRITE(SCARC_LU,'(2i4)') ((SNML%J_MAX(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
   WRITE(SCARC_LU,*) 'SCARC_EXCHANGE:   K_MIN:'
   WRITE(SCARC_LU,'(2i4)') ((SNML%K_MIN(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
   WRITE(SCARC_LU,*) 'SCARC_EXCHANGE:   K_MAX:'
   WRITE(SCARC_LU,'(2i4)') ((SNML%K_MAX(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
ENDIF

         !!! Initialize the communication structures for sending data
         EXCHANGE_INIT_FACE_IF: IF (CODE==NCOM_INIT) THEN
    
            ! initialize communication structures for face exchange
            IF (RNODE/=SNODE) THEN

               IF (ILEVEL/=SNM%NLMAX) THEN
                  NREQ_FACE = NREQ_FACE+1
                  IF (USE_MPI) CALL MPI_ISEND(SNML%IJKW(1,1),15*SNML%N_EXTERNAL_WALL_CELLS,MPI_INTEGER,SNODE, &
                                           TAG_FACE,MPI_COMM_WORLD,REQ_FACE(NREQ_FACE),IERR)
               ENDIF

               IF (SCARC_DEBUG>=6) THEN
                 WRITE(SCARC_LU,*) '============================================================='
                 WRITE(SCARC_LU,*) 'EXCHANGE: INIT: LEVEL=',ILEVEL
                 WRITE(SCARC_LU,*) 'EXCHANGE: INIT: SNML%N_EXTERNAL_WALL_CELLS=',SNML%N_EXTERNAL_WALL_CELLS
                 WRITE(SCARC_LU,*) 'EXCHANGE: INIT: LEN(IJKW)=',15*SNML%N_EXTERNAL_WALL_CELLS
                 !WRITE(SCARC_LU,'(15i4)') SNML%IJKW
                 WRITE(SCARC_LU,*) 'EXCHANGE: INIT: SIZE(IJKW):  ',SIZE(SNML%IJKW)
                 WRITE(SCARC_LU,*) 'SENDING IJKW(',ILEVEL,'): '
                 DO IW= 1,SNML%N_EXTERNAL_WALL_CELLS
                    WRITE(SCARC_LU,'(15i4)') (SNML%IJKW(II,IW),II=1,15)
                 ENDDO
               ENDIF

               !ILEN_FACE=(MAX(SNML%NIC(NM, NOM), SNML%NIC(NOM, NM))+2)*2+1   ! extended
               !ILEN_FACE=(MAX(SNML%NIC(NM, NOM), SNML%NIC(NOM, NM)))*2+1 
               IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,*) 'EXCHANGE: ILEN_FACE=',ILEN_FACE
               ALLOCATE (OSNML%SEND_FACE(ILEN_FACE))
               OSNML%SEND_FACE = 0.0_EB

            ENDIF

         ENDIF EXCHANGE_INIT_FACE_IF
   
      !!! Perform exchange for matrix-vector multiplication (only face exchange needed)
         EXCHANGE_MATV_FACE_IF: IF (CODE==NCOM_MATV) THEN
            IF (RNODE/=SNODE) THEN

            IF (SCARC_DEBUG>=6) THEN
               WRITE(SCARC_LU,*) '============================================================='
               WRITE(SCARC_LU,*) 'EXCHANGE: MATV: LEVEL=',ILEVEL
               WRITE(SCARC_LU,*) 'EXCHANGE: NM=', NM, ': NOM=',NOM
               WRITE(SCARC_LU,*) 'N_EXTERNAL_WALL_CELLS: ',OSNML%N_EXTERNAL_WALL_CELLS
               WRITE(SCARC_LU,*) 'IJKW: '
               DO IW= 1,OSNML%N_EXTERNAL_WALL_CELLS
                  WRITE(SCARC_LU,'(15i4)') (OSNML%IJKW(II,IW),II=1,15)
               ENDDO
            ENDIF
 
            LL = 0
            IWW = 0
            PACK_SEND_FACE0: DO IW=1,OSNML%N_EXTERNAL_WALL_CELLS
               IF (OSNML%IJKW(9,IW)/=NM) CYCLE PACK_SEND_FACE0
               DO KK=OSNML%IJKW(12,IW),OSNML%IJKW(15,IW)
                  DO JJ=OSNML%IJKW(11,IW),OSNML%IJKW(14,IW)
                     DO II=OSNML%IJKW(10,IW),OSNML%IJKW(13,IW)
                        IWW = IWW + 1
                        OSNML%SEND_FACE(LL+1) = REAL(IW,EB)
                        SELECT CASE(IMV1)
                           CASE(NMV_Y)
                              OSNML%SEND_FACE(LL+2) = SNML%Y(II,JJ,KK)
if (SCARC_DEBUG>=6) WRITE(SCARC_LU,*) 'IW=',IW,': X(',II,',',JJ,',',KK,')=',SNML%X(II,JJ,KK), IMV1
                           CASE(NMV_G)
                              OSNML%SEND_FACE(LL+2) = SNML%G(II,JJ,KK)
if (SCARC_DEBUG>=6) WRITE(SCARC_LU,*) 'IW=',IW,': G(',II,',',JJ,',',KK,')=',SNML%G(II,JJ,KK), IMV1
                           CASE(NMV_R)
                              OSNML%SEND_FACE(LL+2) = SNML%R(II,JJ,KK)
if (SCARC_DEBUG>=6) WRITE(SCARC_LU,*) 'IW=',IW,': R(',II,',',JJ,',',KK,')=',SNML%R(II,JJ,KK), IMV1
                           CASE(NMV_D)
                              OSNML%SEND_FACE(LL+2) = SNML%D(II,JJ,KK)
if (SCARC_DEBUG>=6) WRITE(SCARC_LU,*) 'IW=',IW,': D(',II,',',JJ,',',KK,')=',OSNML%SEND_FACE(LL+2)
                           CASE(NMV_X)
                              OSNML%SEND_FACE(LL+2) = SNML%X(II,JJ,KK)
if (SCARC_DEBUG>=6) WRITE(SCARC_LU,*) 'IW=',IW,': X(',II,',',JJ,',',KK,')=',SNML%X(II,JJ,KK), IMV1
                           CASE(NMV_Z)
                              OSNML%SEND_FACE(LL+2) = SNML%Z(II,JJ,KK)
if (SCARC_DEBUG>=6) WRITE(SCARC_LU,*) 'IW=',IW,': Z(',II,',',JJ,',',KK,')=',SNML%Z(II,JJ,KK), IMV1
                           CASE(NMV_X2)
                              OSNML%SEND_FACE(LL+2) = SNML%X2(II,JJ,KK)
if (SCARC_DEBUG>=6) WRITE(SCARC_LU,*) 'IW=',IW,': X2(',II,',',JJ,',',KK,')=',SNML%X2(II,JJ,KK), IMV1
                           CASE(NMV_D2)
                              OSNML%SEND_FACE(LL+2) = SNML%D2(II,JJ,KK)
if (SCARC_DEBUG>=6) WRITE(SCARC_LU,*) 'IW=',IW,': D2(',II,',',JJ,',',KK,')=',SNML%D2(II,JJ,KK), IMV1
                           CASE(NMV_R2)
                              OSNML%SEND_FACE(LL+2) = SNML%R2(II,JJ,KK)
if (SCARC_DEBUG>=6) WRITE(SCARC_LU,*) 'IW=',IW,': R2(',II,',',JJ,',',KK,')=',SNML%R2(II,JJ,KK), IMV1
                           CASE(NMV_Y2)
                              OSNML%SEND_FACE(LL+2) = SNML%Y2(II,JJ,KK)
if (SCARC_DEBUG>=6) WRITE(SCARC_LU,*) 'IW=',IW,': R2(',II,',',JJ,',',KK,')=',SNML%Y2(II,JJ,KK), IMV1
                        END SELECT
                        LL = LL+2
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO PACK_SEND_FACE0
            OSNML%SEND_FACE(IWW*2+1) = -999.0_EB
            NREQ_FACE=NREQ_FACE+1
            IF (USE_MPI) CALL MPI_ISEND(OSNML%SEND_FACE(1),IWW*2+1,MPI_DOUBLE_PRECISION,SNODE, &
                                        TAG_FACE,MPI_COMM_WORLD,REQ_FACE(NREQ_FACE),IERR)

if (SCARC_DEBUG>=8) THEN
   WRITE(SCARC_LU,*) 'EXCHANGE: MATV: sending ', IWW*2+1
   WRITE(SCARC_LU,*) 'EXCHANGE: MATV: SNODE=',SNODE
   WRITE(SCARC_LU,*) 'EXCHANGE: MATV: TAG_FACE=',TAG_FACE
   WRITE(SCARC_LU,*) 'EXCHANGE: MATV: NREQ_FACE=',NREQ_FACE
   WRITE(SCARC_LU,*) 'EXCHANGE: MATV: REQ_FACE=',REQ_FACE(NREQ_FACE)
ENDIF
if (SCARC_DEBUG>=10) THEN
   WRITE(SCARC_LU,*) 'EXCHANGE: MATV: SEND_FACE(1)'
   WRITE(SCARC_LU,'(2e20.8)') SCARC(NM)%OSCARC(NOM)%SLEVEL(1)%SEND_FACE
   WRITE(SCARC_LU,*) 'SIZE(OSNML%SEND_FACE(1))=',SIZE(SCARC(NM)%OSCARC(NOM)%SLEVEL(1)%SEND_FACE)
ENDIF
if (SCARC_DEBUG>=8) THEN
   WRITE(SCARC_LU,*) 'EXCHANGE: MATV: SEND_FACE(2)'
   WRITE(SCARC_LU,'(2e20.8)') SCARC(NM)%OSCARC(NOM)%SLEVEL(2)%SEND_FACE
   WRITE(SCARC_LU,*) 'SIZE(OSNML%SEND_FACE(1))=',SIZE(SCARC(NM)%OSCARC(NOM)%SLEVEL(2)%SEND_FACE)
ENDIF
            ELSE
                write(*,*) 'Not yet implemented'
            ENDIF
         ENDIF  EXCHANGE_MATV_FACE_IF
   
      !!! Perform full exchange including ghost cells
         EXCHANGE_FULL_FACE_IF: IF (CODE==NCOM_FULL) THEN

      IF (SCARC_DEBUG>=6) CALL SCARC_SHOW0(SNML%Z,'RECEIVE','Z before all')
            IF (RNODE/=SNODE) THEN 

               !OSNML%SEND_FACE=0.0_EB
               LL = 0
               DO KK=KMIN,KMAX
                  DO JJ=JMIN,JMAX
                     DO II=IMIN,IMAX
                        OSNML%SEND_FACE(LL+1) = SNML%Z(II,JJ,KK)
if (SCARC_DEBUG>=6) WRITE(SCARC_LU,*) 'SEND_FACE(',LL+1,')=',SNML%Z(II,JJ,KK), II, JJ, KK
                        LL=LL+1
                     ENDDO
                  ENDDO
               ENDDO
            NREQ_FACE=NREQ_FACE+1
if (SCARC_DEBUG>=6) WRITE(SCARC_LU,*) 'EXCHANGE: NREQ_FACE=',NREQ_FACE
if (SCARC_DEBUG>=6) WRITE(SCARC_LU,*) 'EXCHANGE: LENGTH=',SIZE(OSNML%SEND_FACE)
if (SCARC_DEBUG>=6) THEN 
   write(scarc_lu,*) 'VERSCHICKEN'
   write(scarc_lu,'(3e12.4)') (OSNML%SEND_FACE(ii),II=1,SIZE(OSNML%SEND_FACE))
ENDIF
            IF (USE_MPI) CALL MPI_ISEND(OSNML%SEND_FACE(1),SIZE(OSNML%SEND_FACE),MPI_DOUBLE_PRECISION,SNODE, &
                                        TAG_FACE,MPI_COMM_WORLD,REQ_FACE(NREQ_FACE),IERR)



         ELSE
                write(*,*) 'Not yet implemented'
         ENDIF 

         ENDIF EXCHANGE_FULL_FACE_IF
      ENDIF FACE_IF

   ENDDO EXCHANGE_RECV_MESH_LOOP
ENDDO EXCHANGE_SEND_MESH_LOOP


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Information from Mesh NM is RECV'ed by Mesh NOM.  
!!! NOM is the receiver, NM is the sender.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (SCARC_DEBUG>=6) THEN
   WRITE(SCARC_LU,*) 'NREQ_FACE=',NREQ_FACE
   WRITE(SCARC_LU,*) 'REQ_FACE=',REQ_FACE(1:NREQ_FACE)
   WRITE(SCARC_LU,*) 'USE_MPI=',USE_MPI
   !WRITE(SCARC_LU,*) 'REQ_FACE=',REQ_FACE(1:NREQ_FACE)
ENDIF
   
IF (USE_MPI.AND.NREQ_FACE/=0) CALL MPI_WAITALL(NREQ_FACE,REQ_FACE(1:NREQ_FACE),STAT_FACE,IERR)


 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Extract communication data from corresponding RECEIVE-buffers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
EXCHANGE_SEND2_MESH_LOOP: DO NOM=1,NMESHES

   IF (PROCESS(NOM)/=MYID) CYCLE EXCHANGE_SEND2_MESH_LOOP
    
   EXCHANGE_RECV2_MESH_LOOP: DO NM=1,NMESHES
    
      SNODE = PROCESS(NOM)
      RNODE = PROCESS(NM)

      SNOM  => SCARC(NOM)                     ! corresponds to M4
      SNOML => SCARC(NOM)%SLEVEL(ILEVEL)      ! ... for the level 'ILEVEL'
      IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,*) 'SNODE=',SNODE,': RNODE=',RNODE,': NM=',NM,': NOM=',NOM
      IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,*) 'SNOML zeigt auf SCARC(',NOM,')%SLEVEL(',ILEVEL,')'


      RECV_FACE_IF: IF (SNOML%NIC(NOM,NM)/=0 .AND. SNOML%NIC(NM,NOM)/=0) THEN

         OSNOM => SCARC(NOM)%OSCARC(NM)                   ! corresponds to M2
         OSNOML=> SCARC(NOM)%OSCARC(NM)%SLEVEL(ILEVEL)  ! ... for the level 'ILEVEL'

         !!! Extract data from neighbor and complete matrix-vector communication on internal boundary cells
         RECV_FACE_MATV_IF: IF (CODE==NCOM_MATV) THEN
     
            IF (RNODE/=SNODE) THEN
               LL = 0

               UNPACK_RECV_FACE0: DO
                  IW = NINT(OSNOML%RECV_FACE(LL+1))
                  IF (IW==-999) EXIT UNPACK_RECV_FACE0
                  ZSUM=0.0_EB
                  DO KK=SNOML%IJKW(12,IW),SNOML%IJKW(15,IW)
                     DO JJ=SNOML%IJKW(11,IW),SNOML%IJKW(14,IW)
                        DO II=SNOML%IJKW(10,IW),SNOML%IJKW(13,IW)
                           OSNOML%Z_FACE(II,JJ,KK) = OSNOML%RECV_FACE(LL+2)
                           IF (SCARC_DEBUG>=10) WRITE(SCARC_LU,*) &
                              'Z_FACE(',II,',',JJ,',',KK,')=',OSNOML%Z_FACE(II,JJ,KK)
                           ZSUM=ZSUM+OSNOML%RECV_FACE(LL+2)
                           LL = LL+2
                        ENDDO
                     ENDDO
                  ENDDO
                  I=SNOML%IJKW(6,IW)
                  J=SNOML%IJKW(7,IW)
                  K=SNOML%IJKW(8,IW)
                  ISUM = (SNOML%IJKW(13,IW)-SNOML%IJKW(10,IW)+1) * &
                         (SNOML%IJKW(14,IW)-SNOML%IJKW(11,IW)+1) * &
                         (SNOML%IJKW(15,IW)-SNOML%IJKW(12,IW)+1)


                  IF (ABS(SNOML%IJKW(4,IW))==1) THEN
                     ASUB=SNOML%ASUBX
                  ELSE IF (ABS(SNOML%IJKW(4,IW))==2) THEN
                     ASUB=SNOML%ASUBY
                  ELSE IF (ABS(SNOML%IJKW(4,IW))==3) THEN
                     ASUB=SNOML%ASUBZ
                  ENDIF

                  SELECT CASE(IMV2)
                     CASE(NMV_Y)
                        yold=SNOML%Y(I,J,K)              ! only temporarily for debugging
                        SNOML%Y(I, J, K) = SNOML%Y(I, J, K) + ASUB * ZSUM/REAL(ISUM,EB)
IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,2f25.16,a,2f12.6,3i3)') &
               '    Y1(',I,',',J,',',K,')=',SNOML%Y(I,J,K), yold, ' :  ASUB=',ASUB, ZSUM, ISUM, IMV2
                     CASE(NMV_G)
                        yold=SNOML%G(I,J,K)
                        SNOML%G(I, J, K) = SNOML%G(I, J, K) + ASUB * ZSUM/REAL(ISUM,EB)
IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,2f25.16,a,2f12.6,3i3)') &
               '    G1(',I,',',J,',',K,')=',SNOML%G(I,J,K), yold, ' :  ASUB=',ASUB, ZSUM, ISUM, IMV2
                     CASE(NMV_R)
                        yold=SNOML%R(I,J,K)
                        SNOML%R(I, J, K) = SNOML%R(I, J, K) + ASUB * ZSUM/REAL(ISUM,EB)
IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,2f25.16,a,2f12.6,3i3)') &
               '    R1(',I,',',J,',',K,')=',SNOML%R(I,J,K), yold, ' :  ASUB=',ASUB, ZSUM, ISUM, IMV2
                     CASE(NMV_D)
                        yold=SNOML%D(I,J,K)
                        SNOML%D(I, J, K) = SNOML%D(I, J, K) + ASUB * ZSUM/REAL(ISUM,EB)
IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,2f25.16,a,2f12.6,3i3)') &
               '    D1(',I,',',J,',',K,')=',SNOML%D(I,J,K), yold, ' :  ASUB=',ASUB, ZSUM, ISUM, IMV2
                     CASE(NMV_X)
                        yold=SNOML%X(I,J,K)
                        SNOML%X(I, J, K) = SNOML%X(I, J, K) + ASUB * ZSUM/REAL(ISUM,EB)
IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,2f25.16,a,2f12.6,3i3)') &
               '    X1(',I,',',J,',',K,')=',SNOML%X(I,J,K), yold, ' :  ASUB=',ASUB, ZSUM, ISUM, IMV2
                     CASE(NMV_X2)
                        yold=SNOML%X2(I,J,K)
                        SNOML%X2(I, J, K) = SNOML%X2(I, J, K) + ASUB * ZSUM/REAL(ISUM,EB)
IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,2f25.16,a,2f12.6,3i3)') &
               '    X2(',I,',',J,',',K,')=',SNOML%X2(I,J,K), yold, ' :  ASUB=',ASUB, ZSUM, ISUM, IMV2
                     CASE(NMV_D2)
                        yold=SNOML%D2(I,J,K)
                        SNOML%D2(I, J, K) = SNOML%D2(I, J, K) + ASUB * ZSUM/REAL(ISUM,EB)
IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,2f25.16,a,2f12.6,3i3)') &
               '    D2(',I,',',J,',',K,')=',SNOML%D2(I,J,K), yold, ' :  ASUB=',ASUB, ZSUM, ISUM, IMV2
                     CASE(NMV_R2)
                        yold=SNOML%R2(I,J,K)
                        SNOML%R2(I, J, K) = SNOML%R2(I, J, K) + ASUB * ZSUM/REAL(ISUM,EB)
IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,2f25.16,a,2f12.6,3i3)') &
               '    R2(',I,',',J,',',K,')=',SNOML%R2(I,J,K), yold, ' :  ASUB=',ASUB, ZSUM, ISUM, IMV2
                     CASE(NMV_Y2)
                        yold=SNOML%Y2(I,J,K)
                        SNOML%Y2(I, J, K) = SNOML%Y2(I, J, K) + ASUB * ZSUM/REAL(ISUM,EB)
IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,2f25.16,a,2f12.6,3i3)') &
               '    Y2(',I,',',J,',',K,')=',SNOML%Y2(I,J,K), yold, ' :  ASUB=',ASUB, ZSUM, ISUM, IMV2
                  END SELECT
               ENDDO UNPACK_RECV_FACE0
            ENDIF

         ENDIF RECV_FACE_MATV_IF


         !!! Extract data from neighbor on ghost cells
         RECV_FACE_FULL_IF: IF (CODE==NCOM_FULL) THEN
   
            IOR0=SNOML%IOR_FACE(NOM,NM)

IF (SCARC_DEBUG>=6) THEN
   write(scarc_lu,*) 'EINSORTIEREN:'
   write(scarc_lu,*) 'OSNOML%RECV_FACE:'
   write(scarc_lu,'(3e12.4)') (OSNOML%RECV_FACE(ii),II=1,SIZE(OSNOML%RECV_FACE))
   WRITE(SCARC_LU,*) 'SNOML zeigt auf SCARC(',NOM,')%SLEVEL(',ILEVEL,')'
   WRITE(SCARC_LU,*) 'NM=',NM, ': IOR0=',IOR0
ENDIF
         IMIN=I_MIN(NOM,NM)
         IMAX=I_MAX(NOM,NM)
         JMIN=J_MIN(NOM,NM)
         JMAX=J_MAX(NOM,NM)
         KMIN=K_MIN(NOM,NM)
         KMAX=K_MAX(NOM,NM)

         OSNOML%Z_FACE=0.0_EB
            IF (RNODE/=SNODE) THEN

               LL = 0
               DO KK=KMIN,KMAX
                  DO JJ=JMIN,JMAX
                     DO II=IMIN,IMAX
                        OSNOML%Z_FACE(II,JJ,KK) = OSNOML%RECV_FACE(LL+1)
if (SCARC_DEBUG>=6) WRITE(SCARC_LU,*) 'RECV_FACE(',II,',',JJ,',',KK,')=',OSNOML%Z_FACE(II,JJ,KK)
                        LL=LL+1
                     ENDDO
                  ENDDO
               ENDDO

               SELECT CASE(IOR0)
               CASE(1)
                  IMIN=SNOML%IBAR
                  IMAX=SNOML%IBAR
                  IF (TWO_D) THEN
                     JMIN=1
                     JMAX=1
                  ELSE
                     JMIN=SNOML%J_MIN(NOM,NM)
                     JMAX=SNOML%J_MAX(NOM,NM)
                  ENDIF
                  KMIN=SNOML%K_MIN(NOM,NM)
                  KMAX=SNOML%K_MAX(NOM,NM)
                  NPTR=0
               CASE(-1)
                  IMIN=1
                  IMAX=1
                  IF (TWO_D) THEN
                     JMIN=1
                     JMAX=1
                  ELSE
                     JMIN=SNOML%J_MIN(NOM,NM)
                     JMAX=SNOML%J_MAX(NOM,NM)
                  ENDIF
                  KMIN=SNOML%K_MIN(NOM,NM)
                  KMAX=SNOML%K_MAX(NOM,NM)
                  NPTR=SNOML%IBAR+1
               CASE(2)
                  IMIN=SNOML%I_MIN(NOM,NM)
                  IMAX=SNOML%I_MAX(NOM,NM)
                  IF (TWO_D) THEN
                     JMIN=1
                     JMAX=1
                  ELSE
                     JMIN=SNOML%JBAR
                     JMAX=SNOML%JBAR
                  ENDIF
                  KMIN=SNOML%K_MIN(NOM,NM)
                  KMAX=SNOML%K_MAX(NOM,NM)
                  NPTR=0
               CASE(-2)
                  IMIN=SNOML%I_MIN(NOM,NM)
                  IMAX=SNOML%I_MAX(NOM,NM)
                  JMIN=1
                  JMAX=1
                  KMIN=SNOML%K_MIN(NOM,NM)
                  KMAX=SNOML%K_MAX(NOM,NM)
                  NPTR=SNOML%JBAR+1
               CASE(3)
                  IMIN=SNOML%I_MIN(NOM,NM)
                  IMAX=SNOML%I_MAX(NOM,NM)
                  IF (TWO_D) THEN
                     JMIN=1
                     JMAX=1
                  ELSE
                     JMIN=SNOML%J_MIN(NOM,NM)
                     JMAX=SNOML%J_MAX(NOM,NM)
                  ENDIF
                  KMIN=SNOML%KBAR
                  KMAX=SNOML%KBAR
                  NPTR=0
               CASE(-3)
                  IMIN=SNOML%I_MIN(NOM,NM)
                  IMAX=SNOML%I_MAX(NOM,NM)
                  IF (TWO_D) THEN
                     JMIN=1
                     JMAX=1
                  ELSE
                     JMIN=SNOML%J_MIN(NOM,NM)
                     JMAX=SNOML%J_MAX(NOM,NM)
                  ENDIF
                  KMIN=1
                  KMAX=1
                  NPTR=SNOML%KBAR+1
               END SELECT
IF (SCARC_DEBUG>=6) THEN
   write(scarc_lu,*) 'IMIN=',IMIN
   write(scarc_lu,*) 'IMAX=',IMAX
   write(scarc_lu,*) 'JMIN=',JMIN
   write(scarc_lu,*) 'JMAX=',JMAX
   write(scarc_lu,*) 'KMIN=',KMIN
   write(scarc_lu,*) 'KMAX=',KMAX
ENDIF

               DO KK=KMIN,KMAX
                  DO JJ=JMIN,JMAX
                     DO II=IMIN,IMAX
                        SELECT CASE(ABS(IOR0))
                        CASE(1)
                           SNOML%Z(NPTR,JJ,KK) = OSNOML%Z_FACE(II,JJ,KK)
if (SCARC_DEBUG>=6) WRITE(SCARC_LU,*) '1:Z(',NPTR,',',JJ,',',KK,')=Z_FACE(',II,',',JJ,',',KK,')=',OSNOML%Z_FACE(II,JJ,KK)
                        CASE (2)
                           SNOML%Z(II,NPTR,KK) = OSNOML%Z_FACE(II,JJ,KK)
if (SCARC_DEBUG>=6) WRITE(SCARC_LU,*) '2:Z(',II,',',NPTR,',',KK,')=Z_FACE(',II,',',JJ,',',KK,')=',OSNOML%Z_FACE(II,JJ,KK)
                        CASE (3)
                           SNOML%Z(II,JJ,NPTR) = OSNOML%Z_FACE(II,JJ,KK)
if (SCARC_DEBUG>=6) WRITE(SCARC_LU,*) '3:Z(',II,',',JJ,',',NPTR,')=Z_FACE(',II,',',JJ,',',KK,')=',OSNOML%Z_FACE(II,JJ,KK)
                        END SELECT
                     ENDDO
                  ENDDO
               ENDDO

            ENDIF

         ENDIF RECV_FACE_FULL_IF

      ENDIF RECV_FACE_IF

   ENDDO EXCHANGE_RECV2_MESH_LOOP
   
ENDDO EXCHANGE_SEND2_MESH_LOOP


TUSED_SCARC(23,MYID+1)=TUSED_SCARC(23,MYID+1)+SECOND()-TNOW_EXCHANGE
TUSED_SCARC(0,MYID+1) =TUSED_SCARC(0,MYID+1) +SECOND()-TNOW_EXCHANGE
 
END SUBROUTINE SCARC_EXCHANGE
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! 2D: After the pressure solution, set the ghost cells correctly
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_GHOSTCELLS (HP,NM)
REAL (EB), POINTER, DIMENSION (:, :, :) :: HP
INTEGER :: NM, IW, IOR0
INTEGER :: I0, J0, K0, I1, J1, K1
REAL(EB):: TNOW_GHOSTCELLS
TYPE (MESH_TYPE),  POINTER ::  M

TNOW_GHOSTCELLS = SECOND()

M  => MESHES(NM)

!CALL SCARC_SHOW (HP, 'SARC', 'GHOST1')
 
DO IW = 1, M%N_EXTERNAL_WALL_CELLS

   I0 = M%IJKW(1,IW)
   J0 = M%IJKW(2,IW)
   K0 = M%IJKW(3,IW)

   I1 = M%IJKW(6,IW)
   J1 = M%IJKW(7,IW)
   K1 = M%IJKW(8,IW)

   IOR0 = M%IJKW(4,IW)

   IF (IOR0 == 1) THEN
      IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
         HP(I0,J0,K0) = -HP(I1,J1,K1) + 2.0_EB * M%BXS(J1,K1)
      ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
         HP(I0,J0,K0) =  HP(I1,J1,K1) - DXI *M%BXS(J1,K1)
      ENDIF
   ELSE IF (IOR0 == -1) THEN
      IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
         HP(I0,J0,K0) = -HP(I1,J1,K1) + 2.0_EB * M%BXF(J1,K1)
      ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
         HP(I0,J0,K0) =  HP(I1,J1,K1) + DXI *M%BXF(J1,K1)
      ENDIF
   ELSE IF (IOR0 == 2) THEN
      IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
         HP(I0,J0,K0) = -HP(I1,J1,K1) + 2.0_EB * M%BYS(I1,K1)
      ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
         HP(I0,J0,K0) =  HP(I1,J1,K1) - DETA *M%BYS(I1,K1)
      ENDIF
   ELSE IF (IOR0 == -2) THEN
      IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
         HP(I0,J0,K0) = -HP(I1,J1,K1) + 2.0_EB * M%BYF(I1,K1)
      ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
         HP(I0,J0,K0) =  HP(I1,J1,K1) + DETA *M%BYF(I1,K1)
      ENDIF
   ELSE IF (IOR0 == 3) THEN
      IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
         HP(I0,J0,K0) = -HP(I1,J1,K1) + 2.0_EB * M%BZS(I1,J1)
         IF (SCARC_DEBUG >=2) WRITE(SCARC_LU,2001) IOR0, I0,J0,K0,I1,J1,K1,2.0_EB, M%BZS(I1,J1),HP(I0,J0,K0)
      ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
         HP(I0,J0,K0) =  HP(I1,J1,K1) - DZETA *M%BZS(I1,J1)
         IF (SCARC_DEBUG >=2) WRITE(SCARC_LU,2002) IOR0, I0,J0,K0,I1,J1,K1,DZETA  , M%BZS(I1,J1),HP(I0,J0,K0)
      ENDIF
   ELSE IF (IOR0 == -3) THEN
      IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
         HP(I0,J0,K0) = -HP(I1,J1,K1) + 2.0_EB * M%BZF(I1,J1)
         IF (SCARC_DEBUG >=2) WRITE(SCARC_LU,2001) IOR0, I0,J0,K0,I1,J1,K1,2.0_EB, M%BZF(I1,J1),HP(I0,J0,K0)
      ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
         HP(I0,J0,K0) =  HP(I1,J1,K1) + DZETA *M%BZF(I1,J1)
         IF (SCARC_DEBUG >=2) WRITE(SCARC_LU,2002) IOR0, I0,J0,K0,I1,J1,K1,DZETA , M%BZF(I1,J1),HP(I0,J0,K0)
      ENDIF
   ENDIF
   IF (SCARC_DEBUG >=8) WRITE(SCARC_LU,1000) I0,J0,K0,I1,J1,K1,IOR0,HP(I0,J0,K0)

ENDDO
 
!CALL SCARC_SHOW (HP, 'SARC', 'GHOST2')

TUSED_SCARC(25,NM)=TUSED_SCARC(25,NM)+SECOND()-TNOW_GHOSTCELLS
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_GHOSTCELLS
1000 FORMAT(' X(',i2,',',i2,',',i2,')=X(',i2,',',i2,',',i2,') + bdry(',i2,') = ',e16.8)
2001 FORMAT('D: IOR=',i3,': HP(',i2,',',i2,',',i2,')=HP(',i2,',',i2,',',i2,') + ',e16.8,'*',e16.8,') = ',e16.8)
2002 FORMAT('N: IOR=',i3,': HP(',i2,',',i2,',',i2,')=HP(',i2,',',i2,',',i2,') + ',e16.8,'*',e16.8,') = ',e16.8)
END SUBROUTINE SCARC_GHOSTCELLS
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! 2D: Set correct boundary conditions for S%AGLOB at exterior boundaries
!!! This is done corresponding to pois.f90
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETBDRY2D (ILEVEL, NM)
INTEGER :: NM, ILEVEL, IW, IOR0, I, J, K
REAL(EB):: TNOW_SETBDRY2D
 
TYPE (SCARC_LEVEL_TYPE),  POINTER ::  SL
TYPE (MESH_TYPE),  POINTER ::  M

TNOW_SETBDRY2D = SECOND()

SL => S%SLEVEL(ILEVEL)
M  => MESHES(NM)
 
DO IW = 1, M%N_EXTERNAL_WALL_CELLS

   I = M%IJKW(6,IW)
   J = M%IJKW(7,IW)
   K = M%IJKW(8,IW)

   IOR0 = M%IJKW(4,IW)

   SELECT CASE(IOR0)
      CASE(1)
         IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
            SL%F(I,J,K) = SL%F(I,J,K) - 2.0_EB * SL%DXI2 * M%BXS(1,K)            ! Dirichlet
            IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,1000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
         ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
            SL%F(I,J,K) = SL%F(I,J,K) + SL%DXI * M%BXS(1,K)                       ! Neumann
            IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,2000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
         ENDIF
      CASE(-1)
         IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
            SL%F(I,J,K) = SL%F(I,J,K) - 2.0_EB * SL%DXI2 *M%BXF(1,K)            ! Dirichlet
            IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,1000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
         ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
            SL%F(I,J,K) = SL%F(I,J,K) - SL%DXI *M%BXF(1,K)                      ! Neumann
            IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,2000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
         ENDIF
      CASE(3)
         IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
            IF (SCARC_DEBUG>=2) THEN
               WRITE(SCARC_LU,*) 'SL%F(',I,',',J,',',K,')=',SL%F(I,J,K)
               WRITE(SCARC_LU,*) 'SL%DZI2=',SL%DZI2
               WRITE(SCARC_LU,*) 'M%BZS=',M%BZS
            ENDIF
            SL%F(I,J,K) = SL%F(I,J,K) - 2.0_EB * SL%DZI2 * M%BZS(I,1)            ! Dirichlet
            IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,1000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
         ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
            SL%F(I,J,K) = SL%F(I,J,K) + SL%DZI * M%BZS(I,1)                       ! Neumann
            IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,2000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
         ENDIF
      CASE(-3)
         IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
            SL%F(I,J,K) = SL%F(I,J,K) - 2.0_EB * SL%DZI2 * M%BZF(I,1)            ! Dirichlet
            IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,1000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
         ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
            SL%F(I,J,K) = SL%F(I,J,K) - SL%DZI  * M%BZF(I,1)                     ! Neumann
            IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,2000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
         ENDIF
   END SELECT

ENDDO

 
TUSED_SCARC(26,NM)=TUSED_SCARC(26,NM)+SECOND()-TNOW_SETBDRY2D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_SETBDRY2D

1000 FORMAT(i3,': IOR=',i3,': Dirichlet :F(',i2,',',i2,',',i2,')=',e25.16,i5,i5)
2000 FORMAT(i3,': IOR=',i3,': Neumann   :F(',i2,',',i2,',',i2,')=',e25.16,i5,i5)
END SUBROUTINE SCARC_SETBDRY2D
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! 3D: Set correct boundary conditions for S%AGLOB at exterior boundaries
!!! This is done corresponding to pois.f90
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETBDRY3D (ILEVEL, NM)
INTEGER :: NM, ILEVEL, IW, IOR0, I, J, K
REAL(EB):: TNOW_SETBDRY3D
 
TYPE (SCARC_LEVEL_TYPE),  POINTER ::  SL
TYPE (MESH_TYPE),  POINTER ::  M

TNOW_SETBDRY3D = SECOND()

SL => S%SLEVEL(ILEVEL)
M  => MESHES(NM)
 
DO IW = 1, M%N_EXTERNAL_WALL_CELLS

   I = M%IJKW(6,IW)
   J = M%IJKW(7,IW)
   K = M%IJKW(8,IW)

   IOR0 = M%IJKW(4,IW)

   SELECT CASE(IOR0)
      CASE(1)
         IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
            SL%F(I,J,K) = SL%F(I,J,K) - 2.0_EB * SL%DXI2 * M%BXS(J,K)            ! Dirichlet
            IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,1000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW), M%BXS(J,K)
         ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
            SL%F(I,J,K) = SL%F(I,J,K) + SL%DXI * M%BXS(J,K)                       ! Neumann
            IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,2000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW), M%BXS(J,K)
         ENDIF
      CASE(-1)
         IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
            SL%F(I,J,K) = SL%F(I,J,K) - 2.0_EB * SL%DXI2 *M%BXF(J,K)            ! Dirichlet
            IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,1000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW), M%BXF(J,K)
         ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
            SL%F(I,J,K) = SL%F(I,J,K) - SL%DXI *M%BXF(J,K)                      ! Neumann
            IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,2000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW), M%BXF(J,K)
         ENDIF
      CASE(2)
         IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
            SL%F(I,J,K) = SL%F(I,J,K) - 2.0_EB * SL%DYI2 * M%BYS(I,K)            ! Dirichlet
            IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,1000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW), M%BYS(I,K)
         ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
            SL%F(I,J,K) = SL%F(I,J,K) + SL%DYI * M%BYS(I,K)                       ! Neumann
            IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,2000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW), M%BYS(I,K)
         ENDIF
      CASE(-2)
         IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
            SL%F(I,J,K) = SL%F(I,J,K) - 2.0_EB * SL%DYI2 *M%BYF(I,K)            ! Dirichlet
            IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,1000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW), M%BYF(I,K)
         ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
            SL%F(I,J,K) = SL%F(I,J,K) - SL%DYI *M%BYF(I,K)                      ! Neumann
            IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,2000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW), M%BYF(I,K)
         ENDIF
      CASE(3)
         IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
            SL%F(I,J,K) = SL%F(I,J,K) - 2.0_EB * SL%DZI2 * M%BZS(I,J)            ! Dirichlet
            IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,1000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW), M%BZS(I,J)
         ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
            SL%F(I,J,K) = SL%F(I,J,K) + SL%DZI * M%BZS(I,J)                       ! Neumann
            IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,2000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW), M%BZS(I,J)
         ENDIF
      CASE(-3)
         IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
            SL%F(I,J,K) = SL%F(I,J,K) - 2.0_EB * SL%DZI2 * M%BZF(I,J)            ! Dirichlet
            IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,1000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW), M%BZF(I,J)
         ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
            SL%F(I,J,K) = SL%F(I,J,K) - SL%DZI  * M%BZF(I,J)                     ! Neumann
            IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,2000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW), M%BZF(I,J)
         ENDIF
   END SELECT


ENDDO
 
 
TUSED_SCARC(27,NM)=TUSED_SCARC(27,NM)+SECOND()-TNOW_SETBDRY3D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_SETBDRY3D
1000 FORMAT(i3,': IOR=',i3,': Dirichlet :F(',i2,',',i2,',',i2,')=',e25.16,i5,i5,e25.16)
2000 FORMAT(i3,': IOR=',i3,': Neumann   :F(',i2,',',i2,',',i2,')=',e25.16,i5,i5,e25.16)
END SUBROUTINE SCARC_SETBDRY3D
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!  update interior boundary values along vertical directions
!!!  to have consistent values all over the domain
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_UPDATE(Z, NM, CNAME)
REAL (EB), POINTER, DIMENSION (:, :, :) :: Z
INTEGER :: NM, I, J, K, ILMAX
CHARACTER (6) :: CNAME
REAL(EB):: TNOW_UPDATE_QUANTITY
 
TNOW_UPDATE_QUANTITY = SECOND()

ILMAX=SCARC(NM)%NLMAX
SL => SCARC(NM)%SLEVEL(ILMAX)
 
IF (SCARC_DEBUG>=2) THEN
   WRITE(SCARC_LU,*) '========================= SCARC_UPDATE'
   WRITE(SCARC_LU,*) 'ILMAX=',ILMAX
   WRITE(SCARC_LU,*) 'NM=',NM
   WRITE(SCARC_LU,*) 'CNAME=',CNAME
ENDIF

IF (SCARC_DEBUG>=6) CALL SCARC_SHOW0(Z,'SCARC_UPDATE','Z before update')

IF (NMESHES > 1) THEN
 
   DO K = 0, SL%KBAR + 1
      DO J = 0, SL%JBAR + 1
         DO I = 0, SL%IBAR + 1
            SL%Z (I, J, K) = Z (I, J, K)
         ENDDO
      ENDDO
   ENDDO
 
   NREQ_FACE = 0
 
!CALL SCARC_SHOW0(SL%Z,'SCARC_UPDATE','SL%Z before update')

   CALL SCARC_RECEIVE  (NCOM_FULL, ILMAX) 
   CALL SCARC_EXCHANGE (NCOM_FULL, ILMAX, NMV_NONE, NMV_NONE)

!CALL SCARC_SHOW0(SL%Z,'SCARC_UPDATE','SL%Z after update')
 
   DO K = 0, SL%KBAR + 1
      DO J = 0, SL%JBAR + 1
         DO I = 0, SL%IBAR + 1
            Z (I, J, K) = SL%Z(I, J, K)
         ENDDO
      ENDDO
   ENDDO
 
ENDIF
 
IF (SCARC_DEBUG>=6) CALL SCARC_SHOW0(Z,'SCARC_UPDATE','Z after update')

 
TUSED_SCARC(30,NM)=TUSED_SCARC(30,NM)+SECOND()-TNOW_UPDATE_QUANTITY
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_UPDATE_QUANTITY
END SUBROUTINE SCARC_UPDATE

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!  update interior boundary values along vertical directions
!!!  to have consistent values all over the domain
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_UPDATE_LEVEL(Z, NM, ILEVEL, CNAME)
REAL (EB), POINTER, DIMENSION (:, :, :) :: Z
INTEGER :: NM, ILEVEL, I, J, K
CHARACTER (6) :: CNAME
REAL(EB):: TNOW_UPDATE_QUANTITY
 
TNOW_UPDATE_QUANTITY = SECOND()

SL => S%SLEVEL(ILEVEL)
!CALL SCARC_SHOW_LEVEL (Z, 'UPDT', CNAME, ILEVEL)
 
IF (NMESHES > 1) THEN
 
WRITE(SCARC_LU,*) 'SHAPE(Z)=',SHAPE(Z)
WRITE(SCARC_LU,*) 'SHAPE(SL%Z)=',SHAPE(SL%Z)

   DO K = 0, SL%KBAR + 1
      DO J = 0, SL%JBAR + 1
         DO I = 0, SL%IBAR + 1
            SL%Z (I, J, K) = Z (I, J, K)
         ENDDO
      ENDDO
   ENDDO
 
   NREQ_FACE = 0

   CALL SCARC_RECEIVE (NCOM_FULL, ILEVEL)  
   CALL SCARC_EXCHANGE(NCOM_FULL, ILEVEL, NMV_NONE, NMV_NONE)
 
   DO K = 0, SL%KBAR + 1
      DO J = 0, SL%JBAR + 1
         DO I = 0, SL%IBAR + 1
            Z (I, J, K) = SL%Z(I, J, K)
         ENDDO
      ENDDO
   ENDDO
 
 
ENDIF
CALL SCARC_SHOW_LEVEL (Z, 'UPDT', CNAME, ILEVEL)
 
TUSED_SCARC(30,NM)=TUSED_SCARC(30,NM)+SECOND()-TNOW_UPDATE_QUANTITY
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_UPDATE_QUANTITY
END SUBROUTINE SCARC_UPDATE_LEVEL
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!  update interior boundary values along vertical directions
!!!  to have consistent values all over the domain
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_UPDATE_FULL (Z, NM, ILEVEL, CNAME)
REAL (EB), POINTER, DIMENSION (:, :, :) :: Z
INTEGER :: NM, ILEVEL, I, J, K
CHARACTER (4) :: CNAME
REAL(EB):: TNOW_UPDATE_QUANTITY
 
TNOW_UPDATE_QUANTITY = SECOND()

SL => SCARC(NM)%SLEVEL(2)
IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,*) 'UPDATING ', CNAME
 
IF (NMESHES > 1) THEN
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! 2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   IF (TWO_D) THEN
 
      DO K = 0, SL%KBAR + 1
         DO I = 0, SL%IBAR + 1
            SL%Z (I, 1, K) = Z (I, 1, K)
         ENDDO
      ENDDO
CALL SCARC_SHOW (SL%Z, 'UPDT', 'ZIN1  ')
      NREQ_FACE = 0
      CALL SCARC_RECEIVE  (NCOM_FULL, ILEVEL)  
      CALL SCARC_EXCHANGE (NCOM_FULL, ILEVEL, NMV_NONE, NMV_NONE)
CALL SCARC_SHOW (SL%Z, 'UPDT', 'ZOUT1 ')
      DO K = 0, SL%KBAR + 1
         DO I = 0, SL%IBAR + 1
            Z (I, 1, K) = SL%Z(I, 1, K)
         ENDDO
      ENDDO
 
      DO K = 0, SL%KBAR + 1
         DO I = 0, SL%IBAR + 1
            SL%Z (I, 1, K) = Z (I, 0, K)
         ENDDO
      ENDDO
CALL SCARC_SHOW (SL%Z, 'UPDT', 'ZIN0  ')
      NREQ_FACE = 0
      CALL SCARC_RECEIVE  (NCOM_FULL, ILEVEL)  
      CALL SCARC_EXCHANGE (NCOM_FULL, ILEVEL, NMV_NONE, NMV_NONE)
CALL SCARC_SHOW (SL%Z, 'UPDT', 'ZOUT0 ')
      DO K = 0, SL%KBAR + 1
         DO I = 0, SL%IBAR + 1
            Z (I, 0, K) = SL%Z(I, 1, K)
         ENDDO
      ENDDO

      DO K = 0, SL%KBAR + 1
         DO I = 0, SL%IBAR + 1
            SL%Z (I, 1, K) = Z (I, 2, K)
         ENDDO
      ENDDO
CALL SCARC_SHOW (SL%Z, 'UPDT', 'ZIN2  ')
      NREQ_FACE = 0
      CALL SCARC_RECEIVE (NCOM_FULL, ILEVEL)  
      CALL SCARC_EXCHANGE (NCOM_FULL, ILEVEL, NMV_NONE, NMV_NONE)
CALL SCARC_SHOW (SL%Z, 'UPDT', 'ZOUT2 ')
      DO K = 0, SL%KBAR + 1
         DO I = 0, SL%IBAR + 1
            Z (I, 2, K) = SL%Z(I, 1, K)
         ENDDO
      ENDDO
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! 3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ELSE
 
      DO K = 0, KBAR + 1
         DO J = 0, JBAR + 1
            DO I = 0, IBAR + 1
               SL%Z (I, J, K) = Z (I, J, K)
            ENDDO
         ENDDO
      ENDDO
 
      NREQ_FACE = 0
 
      CALL SCARC_RECEIVE (NCOM_FULL, ILEVEL) 
      CALL SCARC_EXCHANGE (NCOM_FULL, ILEVEL, NMV_NONE, NMV_NONE)
 
      DO K = 0, KBAR + 1
         DO J = 0, JBAR + 1
            DO I = 0, IBAR + 1
               Z (I, J, K) = SL%Z(I, J, K)
            ENDDO
         ENDDO
      ENDDO
   ENDIF
ENDIF
 
TUSED_SCARC(30,NM)=TUSED_SCARC(30,NM)+SECOND()-TNOW_UPDATE_QUANTITY
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_UPDATE_QUANTITY
END SUBROUTINE SCARC_UPDATE_FULL
 
 
 
 
LOGICAL FUNCTION MATCH (VAL1, VAL2)
REAL (EB) :: VAL1, VAL2, TOL
TOL = 1.0E-8_EB
MATCH = .FALSE.
IF (Abs(VAL1-VAL2) <= TOL) MATCH = .TRUE.
END FUNCTION MATCH
 
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! only for debugging reasons ... show complete vector including ghost cells
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SHOW (X, CROUTINE, CNAME)
REAL (EB), POINTER, DIMENSION (:, :, :) :: X
REAL (EB):: VALUES(10)
INTEGER :: II, JJ, KK, IBAR9,JBAR9,KBAR9
CHARACTER (4) :: CROUTINE
CHARACTER (6) :: CNAME

IBAR9=MIN(9,IBAR+1)
JBAR9=MIN(9,JBAR+1)
KBAR9=MIN(9,KBAR+1)

IF (SCARC_DEBUG >= 2) THEN
   write(SCARC_LU,*) 'SCARC_LU=',SCARC_LU
   WRITE(SCARC_LU,*) '============================================================='
   WRITE(SCARC_LU,*) '===   COMPARE= ',CNAME,': ROUTINE= ',CROUTINE, IBAR, JBAR, KBAR
   WRITE(SCARC_LU,*) '============================================================='
   IF (TWO_D) THEN
      DO KK = KBAR9, 0, - 1
         WRITE(SCARC_LU, '(a,i3,a,10e12.3)') 'k=', KK, ' : ', (X(II, 1, KK), II=0,IBAR9)
      ENDDO
   ELSE
      !DO KK = KBAR9, 0, - 1
      DO KK = 1, 1, - 1
         WRITE(SCARC_LU,'(a,i3,a)')  '----------------------------------------  K = ', KK,&
                                     '----------------------------------------'
         !DO JJ = JBAR9, 0, - 1
         DO JJ = 5, 0, - 1
            DO II=0,IBAR9
               IF (ABS(X(II,JJ,KK))<1.0E-15_EB) THEN
                  VALUES(II)=0.0_EB
               ELSE
                  VALUES(II)=X(II,JJ,KK)
               ENDIF
            ENDDO
            WRITE(SCARC_LU, '(a,i3,a,10e14.5)') 'J= ',JJ,' : ',(VALUES(II), II=0, IBAR9)
         ENDDO
      ENDDO
      WRITE(SCARC_LU,*)  '------------------------------------------------',&
                         '---------------------------------------------------'
      WRITE(SCARC_LU,1000) ('I = ',II,II=0,IBAR9)
      WRITE(SCARC_LU,*)  '--------------------------------------------------',&
                         '-------------------------------------------------'
      WRITE(SCARC_LU,*)

   ENDIF
ENDIF
1000 FORMAT(13x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x)
END SUBROUTINE SCARC_SHOW


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! only for debugging reasons ... show complete vector including ghost cells
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SHOW0 (X, CROUTINE, CNAME)
REAL (EB), POINTER, DIMENSION (:, :, :) :: X
REAL (EB):: VALUES(10)
INTEGER :: II, JJ, KK, IBAR9,JBAR9,KBAR9
CHARACTER (*) :: CROUTINE
CHARACTER (*) :: CNAME

IBAR9=MIN(9,IBAR+1)
JBAR9=MIN(9,JBAR+1)
KBAR9=MIN(9,KBAR+1)

IF (SCARC_DEBUG >= 1) THEN
   WRITE(SCARC_LU,*) '============================================================='
   WRITE(SCARC_LU,*) '===   ROUTINE = ',TRIM(CROUTINE)
   WRITE(SCARC_LU,*) '===   QUANTITY= ',TRIM(CNAME)
   WRITE(SCARC_LU,*) '============================================================='
   IF (TWO_D) THEN
      DO KK = KBAR9, 0, - 1
         WRITE(SCARC_LU, '(a,i3,a,10e12.3)') 'k=', KK, ' : ', (X(II, 1, KK), II=0,IBAR9)
      ENDDO
      WRITE(SCARC_LU,*)  '------------------------------------------------',&
                         '---------------------------------------------------'
      WRITE(SCARC_LU,1000) ('I = ',II,II=0,IBAR9)
   ELSE
      DO KK = KBAR9, 0, - 1
      !DO KK = 4, 4, - 1
         WRITE(SCARC_LU,'(a,i3,a)')  '----------------------------------------  K = ', KK,&
                                     '----------------------------------------'
         DO JJ = JBAR9, 0, - 1
            DO II=0,IBAR9
               IF (ABS(X(II,JJ,KK))<1.0E-15_EB) THEN
                  VALUES(II)=0.0_EB
               ELSE
                  VALUES(II)=X(II,JJ,KK)
               ENDIF
            ENDDO
            WRITE(SCARC_LU, '(a,i3,a,10e13.5)') 'J= ',JJ,' : ',(VALUES(II), II=0, IBAR9)
         ENDDO
      ENDDO
      WRITE(SCARC_LU,*)  '------------------------------------------------',&
                         '---------------------------------------------------'
      WRITE(SCARC_LU,1000) ('I = ',II,II=0,IBAR9)
      WRITE(SCARC_LU,*)  '--------------------------------------------------',&
                         '-------------------------------------------------'
      WRITE(SCARC_LU,*)

   ENDIF
ENDIF
1000 FORMAT(13x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x)
END SUBROUTINE SCARC_SHOW0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! only for debugging reasons ... show complete vector including ghost cells
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SHOW02 (X, CROUTINE, CNAME)
REAL (EB), POINTER, DIMENSION (:, :, :) :: X
REAL (EB):: VALUES(10)
INTEGER :: II, JJ, KK, IBAR9,JBAR9,KBAR9
CHARACTER (*) :: CROUTINE
CHARACTER (*) :: CNAME

IBAR9=MIN(9,IBAR+1)
JBAR9=MIN(9,JBAR+1)
KBAR9=MIN(9,KBAR+1)

IF (SCARC_DEBUG >= 1) THEN
   WRITE(SCARC_LU,*) '============================================================='
   WRITE(SCARC_LU,*) '===   ROUTINE = ',TRIM(CROUTINE)
   WRITE(SCARC_LU,*) '===   QUANTITY= ',TRIM(CNAME)
   WRITE(SCARC_LU,*) '============================================================='
   IF (TWO_D) THEN
   WRITE(SCARC_LU,*) '========J=0'
      DO KK = KBAR9, 0, - 1
         WRITE(SCARC_LU, '(a,i3,a,10e12.3)') 'k=', KK, ' : ', (X(II, 0, KK), II=0,IBAR9)
      ENDDO
   WRITE(SCARC_LU,*) '========J=1'
      DO KK = KBAR9, 0, - 1
         WRITE(SCARC_LU, '(a,i3,a,10e12.3)') 'k=', KK, ' : ', (X(II, 1, KK), II=0,IBAR9)
      ENDDO
      WRITE(SCARC_LU,*)  '------------------------------------------------',&
                         '---------------------------------------------------'
      WRITE(SCARC_LU,1000) ('I = ',II,II=0,IBAR9)
   WRITE(SCARC_LU,*) '========J=2'
      DO KK = KBAR9, 0, - 1
         WRITE(SCARC_LU, '(a,i3,a,10e12.3)') 'k=', KK, ' : ', (X(II, 2, KK), II=0,IBAR9)
      ENDDO
   ELSE
      DO KK = KBAR9, 0, - 1
      !DO KK = 4, 4, - 1
         WRITE(SCARC_LU,'(a,i3,a)')  '----------------------------------------  K = ', KK,&
                                     '----------------------------------------'
         DO JJ = JBAR9, 0, - 1
            DO II=0,IBAR9
               IF (ABS(X(II,JJ,KK))<1.0E-15_EB) THEN
                  VALUES(II)=0.0_EB
               ELSE
                  VALUES(II)=X(II,JJ,KK)
               ENDIF
            ENDDO
            WRITE(SCARC_LU, '(a,i3,a,10e13.5)') 'J= ',JJ,' : ',(VALUES(II), II=0, IBAR9)
         ENDDO
      ENDDO
      WRITE(SCARC_LU,*)  '------------------------------------------------',&
                         '---------------------------------------------------'
      WRITE(SCARC_LU,1000) ('I = ',II,II=0,IBAR9)
      WRITE(SCARC_LU,*)  '--------------------------------------------------',&
                         '-------------------------------------------------'
      WRITE(SCARC_LU,*)

   ENDIF
ENDIF
1000 FORMAT(13x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x)
END SUBROUTINE SCARC_SHOW02



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! only for debugging reasons ... related to single level
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SHOW_LEVEL (X, CROUTINE, CNAME, ILEVEL)
REAL (EB), POINTER, DIMENSION (:, :, :) :: X
REAL (EB):: VALUES(10)
INTEGER :: ILEVEL
INTEGER :: II, JJ, KK, IBAR9,JBAR9,KBAR9
CHARACTER (4) :: CROUTINE
CHARACTER (6) :: CNAME
 
SL => S%SLEVEL(ILEVEL)

IBAR9=MIN(9,SL%IBAR+1)
JBAR9=MIN(9,SL%JBAR+1)
KBAR9=MIN(9,SL%KBAR+1)

IF (SCARC_DEBUG >= 3) THEN
   WRITE(SCARC_LU,*) '============================================================='
   WRITE(SCARC_LU,*) '===   COMPARE= ',CNAME,': ROUTINE= ',CROUTINE, ': LEVEL=',ILEVEL, &
                     SL%IBAR, SL%JBAR, SL%KBAR
   WRITE(SCARC_LU,*) '============================================================='
   IF (TWO_D) THEN
      DO KK = KBAR9, 0, - 1
         WRITE(SCARC_LU, '(a,i3,a,10e12.3)') 'k=', KK, ' : ', (X(II, 1, KK), II=0,IBAR9)
      ENDDO
   ELSE
      DO KK = KBAR9, 0, - 1
      !DO KK = 1, 1, - 1
         WRITE(SCARC_LU,'(a,i3,a)')  '----------------------------------------  K = ', KK,&
                                     '----------------------------------------'
         DO JJ = JBAR9, 0, - 1
            DO II=0,IBAR9
               IF (ABS(X(II,JJ,KK))<1.0E-15_EB) THEN
                  VALUES(II)=0.0_EB
               ELSE
                  VALUES(II)=X(II,JJ,KK)
               ENDIF
            ENDDO
            WRITE(SCARC_LU, '(a,i3,a,10e14.5)') 'J= ',JJ,' : ',(VALUES(II), II=0, IBAR9)
         ENDDO
      ENDDO
      WRITE(SCARC_LU,*)  '-----------------------------------------------------',&
                         '----------------------------------------------'
      WRITE(SCARC_LU,1000) ('I = ',II,II=0,IBAR9)
      WRITE(SCARC_LU,*)  '------------------------------------------------------',&
                         '---------------------------------------------'
      WRITE(SCARC_LU,*)

   ENDIF
ENDIF
1000 FORMAT(13x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x)
 
END SUBROUTINE SCARC_SHOW_LEVEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! only for debugging reasons ... related to single level
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SHOW_LEVEL_LINE(X, CROUTINE, CNAME, ILEVEL)
REAL (EB), POINTER, DIMENSION (:, :, :) :: X
INTEGER :: ILEVEL
INTEGER :: II, JJ, KK
CHARACTER (4) :: CROUTINE
CHARACTER (6) :: CNAME
 
SL => S%SLEVEL(ILEVEL)

IF (SCARC_DEBUG >= 3) THEN
   WRITE(SCARC_LU,*) '============================================================='
   WRITE(SCARC_LU,*) '===   COMPARE= ',CNAME,': ROUTINE= ',CROUTINE, ': LEVEL=',ILEVEL, SL%IBAR, SL%KBAR
   WRITE(SCARC_LU,*) '============================================================='
   IF (TWO_D) THEN
      DO KK = SL%KBAR+1, 0, - 1
         WRITE(SCARC_LU, '(a,i3,a,2f30.16)') 'k=', KK, ' : ', (X(II, 1, KK), II=4,5)
      ENDDO
   ELSE
      DO KK = 5, 0, - 1
         WRITE(SCARC_LU, '(3f26.16)') ((X(II, JJ, KK), II=3,5), JJ=5,3,-1)
         WRITE(SCARC_LU,*) '----------------------------------------'
      ENDDO
   ENDIF
ENDIF
 
END SUBROUTINE SCARC_SHOW_LEVEL_LINE

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! print out timings for ScaRC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_TIMINGS
INTEGER:: NM, N, IERR, I
INTEGER, ALLOCATABLE, DIMENSION(:) :: COUNTS_SCARC, DISPLS_SCARC, COUNTS_SCARC_TIMERS, DISPLS_SCARC_TIMERS
CHARACTER(25) :: NAME_SCARC(0:N_TIMERS_SCARC)
REAL(EB) :: TPCNT_SCARC(0:N_TIMERS_SCARC)

IERR=0
 
ALLOCATE(COUNTS_SCARC(0:NUMPROCS-1))
ALLOCATE(COUNTS_SCARC_TIMERS(0:NUMPROCS-1))
ALLOCATE(DISPLS_SCARC(0:NUMPROCS-1))
ALLOCATE(DISPLS_SCARC_TIMERS(0:NUMPROCS-1))

COUNTS_SCARC = 0
DO N=0,NUMPROCS-1
   DO NM=1,NMESHES
      IF (PROCESS(NM)==N) COUNTS_SCARC(N) = COUNTS_SCARC(N) + 1
   ENDDO
ENDDO
DISPLS_SCARC(0) = 0
DO N=1,NUMPROCS-1
   DISPLS_SCARC(N) = COUNTS_SCARC(N-1) + DISPLS_SCARC(N-1)
ENDDO
COUNTS_SCARC_TIMERS = COUNTS_SCARC*N_TIMERS_SCARC
DISPLS_SCARC_TIMERS = DISPLS_SCARC*N_TIMERS_SCARC

IF (USE_MPI) CALL MPI_GATHERV(TUSED_SCARC(1,DISPLS_SCARC(MYID)+1),COUNTS_SCARC_TIMERS(MYID),&
                              MPI_DOUBLE_PRECISION,TUSED_SCARC, COUNTS_SCARC_TIMERS,DISPLS_SCARC_TIMERS,&
                              MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

! Printout subroutine timings

IF (MYID==0) THEN
   NAME_SCARC     = 'null'
   NAME_SCARC(0)  = 'SCARC'
   NAME_SCARC(1)  = 'INIT_MESHES'
   NAME_SCARC(2)  = 'INIT_NEIGHBORS'
   NAME_SCARC(3)  = 'INIT_COMMUNICATION'
   NAME_SCARC(4)  = 'INIT_MATRICES'
   NAME_SCARC(5)  = 'INIT_SOLVER'
   NAME_SCARC(6)  = 'CG'
   NAME_SCARC(7)  = 'MG'
   NAME_SCARC(8)  = 'SMOOTHER'
   NAME_SCARC(9)  = 'COARSE'
   NAME_SCARC(10) = 'MATVEC'
   NAME_SCARC(11) = 'RESTRICTION'
   NAME_SCARC(12) = 'PROLONGATION'
   NAME_SCARC(13) = 'BDRY_RESIDUUM'
   NAME_SCARC(16) = 'SCALPROD'
   NAME_SCARC(17) = 'L2NORM'
   NAME_SCARC(18) = 'PRECON_JACOBI'
   NAME_SCARC(19) = 'PRECON_GS'
   NAME_SCARC(20) = 'PRECON_SSOR'
   NAME_SCARC(21) = 'PRECON_GSADI'
   NAME_SCARC(22) = 'RECEIVE'
   NAME_SCARC(23) = 'EXCHANGE'
   NAME_SCARC(24) = 'INITIALIZE'
   NAME_SCARC(25) = 'GHOSTCELLS'
   NAME_SCARC(26) = 'SETBDRY2D'
   NAME_SCARC(27) = 'SETBDRY3D'
   NAME_SCARC(28) = 'UPDATE'
   NAME_SCARC(29) = 'UPDATE_LEVEL'
   NAME_SCARC(30) = 'UPDATE_QUANTITY'
   
   DO NM=1,NMESHES
      DO I=0,N_TIMERS_SCARC
         TPCNT_SCARC(I) = 100._EB*TUSED_SCARC(I,NM)/TUSED(1,NM)
      ENDDO
      WRITE(LU_OUTPUT,443) NM
      WRITE(LU_OUTPUT,444) 'MAIN',TUSED(1,NM),100._EB
      WRITE(LU_OUTPUT,444) (NAME_SCARC(I),TUSED_SCARC(I,NM),TPCNT_SCARC(I),I=0,N_TIMERS_SCARC)
   ENDDO
ENDIF

   IF (SCARC_DEBUG >= 1 ) CLOSE(SCARC_LU)

443 FORMAT(//' ScaRC: CPU Time Usage, Mesh ',I3// &
         '                 CPU (s)        %  '/ &
         '       --------------------------------------------------')
444 FORMAT(7X,A25,2F11.2)

END SUBROUTINE SCARC_TIMINGS

SUBROUTINE GET_REV_scrc(MODULE_REV,MODULE_DATE)
INTEGER,INTENT(INOUT) :: MODULE_REV
CHARACTER(255),INTENT(INOUT) :: MODULE_DATE

WRITE(MODULE_DATE,'(A)') scrcrev(INDEX(scrcrev,':')+1:LEN_TRIM(scrcrev)-2)
READ (MODULE_DATE,'(I5)') MODULE_REV
WRITE(MODULE_DATE,'(A)') scrcdate

END SUBROUTINE GET_REV_scrc

 
END MODULE SCARC_SOLVER
 
