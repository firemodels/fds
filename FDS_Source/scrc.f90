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
PUBLIC SCARC_TIMINGS, SCARC_INITIALIZE
PUBLIC SCARC_SHOW   , SCARC_SHOW_LINE
PUBLIC SCARC_UPDATE , SCARC_UPDATE_LEVEL
PUBLIC SCARC_CG2D   , SCARC_CG3D
PUBLIC SCARC_MG2D   , SCARC_MG3D
PUBLIC SCARC_BICG2D , SCARC_BICG3D
 
!!! Public variables (needed in main, read, pres, divg)
PUBLIC SCARC_METHOD   , SCARC_DEBUG     , SCARC_CASE
PUBLIC SCARC_EPS_REL  , SCARC_EPS_DIVG  , SCARC_BREL    
PUBLIC SCARC_MG_NIT   , SCARC_CG_NIT    , SCARC_BICG_NIT    , SCARC_SM_NIT    , SCARC_CO_NIT    
PUBLIC SCARC_MG_EPS   , SCARC_CG_EPS    , SCARC_BICG_EPS    , SCARC_SM_EPS    , SCARC_CO_EPS
PUBLIC SCARC_MG_PRECON, SCARC_CG_PRECON , SCARC_BICG_PRECON , SCARC_SM_PRECON , SCARC_CO_PRECON
PUBLIC SCARC_MG_OMEGA , SCARC_CG_OMEGA  , SCARC_BICG_OMEGA  , SCARC_SM_OMEGA  , SCARC_CO_OMEGA
PUBLIC SCARC_MG_NLDIFF


CHARACTER (40) :: SCARC_FN='msg/   _scarc'           ! file name for ScaRC debug messages
CHARACTER (10) :: SCARC_METHOD='null'                ! name of method for the solution of the pressure equation

INTEGER        :: SCARC_DEBUG=0,     &               ! debug level (0: no debug messages)
                  SCARC_CASE=0,      &               ! choose different initial solutions
                  SCARC_COUNT=0                      ! counter  for comparison vectors

CHARACTER (10) :: SCARC_MG_PRECON  ='JACOBI', &      ! smoother for mg-method 
                  SCARC_CG_PRECON  ='JACOBI', &      ! preconditioner for cg-method         
                  SCARC_SM_PRECON  ='JACOBI', &      ! ...                smoother           
                  SCARC_CO_PRECON  ='JACOBI', &      ! ...                coarse grid solver 
                  SCARC_BICG_PRECON='JACOBI'         ! ...                bicg-method        

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

REAL (EB)      :: SCARC_MG_OMEGA  =1.E+0_EB, &       ! relaxation parameter for mg
                  SCARC_CG_OMEGA  =1.E+0_EB, &       ! ...                      cg
                  SCARC_SM_OMEGA  =1.E+0_EB, &       ! ...                      smoother
                  SCARC_CO_OMEGA  =1.E+0_EB, &       ! ...                      coarse grid solver
                  SCARC_BICG_OMEGA=1.E+0_EB          ! ...                      bicg

INTEGER        :: SCARC_MG_NLDIFF = 2                ! maximum level difference for mg-method
INTEGER        :: SCARC_LU                           ! unit number for Scarc debug file

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
                      NMV_D2  = 8                    ! ...

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
   INTEGER :: NW_FACE

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

   ! neighbourship structures along faces, edges and vertices
   INTEGER, POINTER, DIMENSION (:, :) :: NIC_FACE
   INTEGER, POINTER, DIMENSION (:, :) :: IJKW_FACE
   INTEGER, POINTER, DIMENSION (:, :) :: I_MIN_FACE, I_MAX_FACE
   INTEGER, POINTER, DIMENSION (:, :) :: J_MIN_FACE, J_MAX_FACE
   INTEGER, POINTER, DIMENSION (:, :) :: K_MIN_FACE, K_MAX_FACE

   ! matrices and iteration vectors for global defect correction on corresponding level
   REAL (EB), POINTER, DIMENSION (:)       :: DD, DL, DU, LD, DAUX, PERIOD
   REAL (EB), POINTER, DIMENSION (:, :)    :: AG, AL
   REAL (EB), POINTER, DIMENSION (:, :, :) :: X, F, D, Y, G, R, Z, X2, F2, D2, TMP
   REAL (EB), POINTER, DIMENSION (:, :)    :: BXS0, BXF0, BYS0, BYF0, BZS0, BZF0
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
   INTEGER :: NW_FACE

   ! communication vectors
   REAL (EB), POINTER, DIMENSION (:) :: SEND_FACE, RECV_FACE

   ! neighbourship structures along faces, edges and vertices
   INTEGER, POINTER, DIMENSION (:, :) :: NIC_FACE
   INTEGER, POINTER, DIMENSION (:, :) :: IJKW_FACE
   INTEGER, POINTER, DIMENSION (:, :) :: I_MIN_FACE, I_MAX_FACE
   INTEGER, POINTER, DIMENSION (:, :) :: J_MIN_FACE, J_MAX_FACE
   INTEGER, POINTER, DIMENSION (:, :) :: K_MIN_FACE, K_MAX_FACE

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
   INTEGER, POINTER, DIMENSION (:, : , : ) :: NDIAG_NBR
   INTEGER, POINTER, DIMENSION (:)         :: MIBAR, MJBAR, MKBAR
   INTEGER, POINTER, DIMENSION (:, :)      :: MLEVEL0
   INTEGER, POINTER, DIMENSION (:)         :: MLEVEL
   INTEGER, POINTER, DIMENSION (:, : )     :: KCYCLE

   INTEGER, POINTER, DIMENSION (:, :) :: IOR_FACE

   ! neighboring ScaRC-structures
   TYPE (OSCARC_TYPE), POINTER, DIMENSION (:) :: OSCARC

   ! different mesh levels (in case of global multigrid)
   INTEGER :: NFACE0, NEDGE0, NDIAG0
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
!!! Allocate SCARC and OSCARC structures
!!!
ALLOCATE (SCARC(NMESHES), STAT=IERR)
CALL CHKMEMERR ('SCARC_INITIALIZE', 'SCARC', IERR)
 
S => SCARC (NM)

ALLOCATE (S%OSCARC(NMESHES), STAT=IERR)
CALL CHKMEMERR ('SCARC_INITIALIZE', 'OSCARC', IERR)

!!!
!!! Open debug-file if requested
!!!
IF (SCARC_DEBUG >= 1) THEN
   WRITE (SCARC_FN(5:7), '(i3.3)') NM
   SCARC_LU = GET_FILE_NUMBER()
   OPEN (SCARC_LU, FILE=SCARC_FN)
ENDIF

IF (SCARC_METHOD=='FFT') GOTO 12345

!!!
!!! Initialize time measurement array
!!!
ALLOCATE(TUSED_SCARC(N_TIMERS_SCARC,NMESHES),STAT=IERR)
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
      IF (MOD(IBAR0,2)/=0) THEN
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
         IF (MOD(JBAR0,2)/=0) THEN
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
      IF (MOD(KBAR0,2)/=0) THEN
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
IF (SCARC_METHOD == 'MG' .OR. SCARC_CG_PRECON == 'MG') THEN       ! mg-structures with complete grid hierarchy 
   S%NLMIN  = S%MLEVEL(NM)-S%NLEVEL+1
ELSE                                                              ! only cg-structures on one grid level
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
CALL CHKMEMERR ('SCARC_NEIGHBORS2D', 'SLEVEL', IERR)


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
   WRITE(*,*) '3D-version not finished yet, do not use ...'
   stop
   !CALL SCARC_INIT_NEIGHBORS3D (NM)
   !CALL SCARC_INIT_MATRICES3D (NM)
   !CALL SCARC_INITIALIZE_MESH_EXCHANGE3D (NM)
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
   SL%DYI=1.0_EB/SL%DY
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

   ALLOCATE (SL%YY(0), STAT=IERR)
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
      WRITE(SCARC_LU,'(a,i3,a,16f7.3)') 'SL(',ILEVEL,')%XX=',(SL%XX(I),I=0,SL%IBAR)
      !WRITE(SCARC_LU,'(a,i3,a,16f7.3)') 'SL(',ILEVEL,')%YY=',(SL%YY(I),I=0,0)
      WRITE(SCARC_LU,'(a,i3,a,16f7.3)') 'SL(',ILEVEL,')%ZZ=',(SL%ZZ(I),I=0,SL%KBAR)
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
INTEGER :: IREFINE
INTEGER :: ILEVEL
INTEGER :: IW_LO
INTEGER :: NOM_HI1, NOM_HI2
INTEGER :: IMIN, IMAX, JMIN, JMAX, KMIN, KMAX
INTEGER :: IMIN1_HI, IMIN2_HI, KMIN1_HI, KMIN2_HI
INTEGER :: IMAX1_HI, IMAX2_HI, KMAX1_HI, KMAX2_HI
INTEGER :: IBAR_NOM, KBAR_NOM
INTEGER :: ILMAX, IOR
INTEGER :: IW1_HI, IW2_HI
INTEGER :: IW_FACE
INTEGER :: IW, IM , MMM, NNN, JM
INTEGER :: I, I_LO, K_LO
INTEGER, ALLOCATABLE, DIMENSION(:) :: COUNTS,DISPLS,COUNTS2D,DISPLS2D
REAL(EB):: TNOW_NEIGHBORS2D
TYPE (MESH_TYPE), POINTER :: M
 
TNOW_NEIGHBORS2D = SECOND()

IF (SCARC_DEBUG >= 2) THEN
   WRITE(SCARC_LU,*) '==========================================================================='
   WRITE(SCARC_LU,*) '=== SCARC_INIT_NEIGHBORS2D ',NM
   WRITE(SCARC_LU,*) '==========================================================================='
   WRITE(SCARC_LU,*) 'NMESHES=', NMESHES
   WRITE(SCARC_LU,*) 'S%NLEVEL=', S%NLEVEL
   WRITE(SCARC_LU,*) 'HIER SCHAUEN, WELCHES GEBRAUCHT WIRD'
   WRITE(SCARC_LU,*) 'NEWC=', M%NEWC
   WRITE(SCARC_LU,*) 'NWC=', M%NWC
ENDIF
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Determine initial values
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
M => MESHES(NM)

IBAR=M%IBAR
KBAR=M%KBAR

S%NFACE0 = 4
S%NDIAG0 = 4

ALLOCATE (S%NDIAG_NBR(2,2,S%NDIAG0))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Determine neighborship structure along vertices
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
NDIAG_NBR_LOOP2D: DO IW=1,M%NEWC

   IF (M%IJKW(2,IW)==1) THEN

      ! lower left vertex
      IF (M%IJKW(1,IW)==0.AND.M%IJKW(3,IW)==1) THEN
         S%NDIAG_NBR(1,1,1)=M%IJKW(9,IW)
         S%NDIAG_NBR(1,2,1)=IW
      ENDIF
      IF (M%IJKW(1,IW)==1.AND.M%IJKW(3,IW)==0) THEN
         S%NDIAG_NBR(2,1,1)=M%IJKW(9,IW)
         S%NDIAG_NBR(2,2,1)=IW
      ENDIF

      ! lower right vertex
      IF (M%IJKW(1,IW)==IBAR  .AND.M%IJKW(3,IW)==0) THEN
         S%NDIAG_NBR(1,1,2)=M%IJKW(9,IW)
         S%NDIAG_NBR(1,2,2)=IW
      ENDIF
      IF (M%IJKW(1,IW)==IBAR+1.AND.M%IJKW(3,IW)==1) THEN
         S%NDIAG_NBR(2,1,2)=M%IJKW(9,IW)
         S%NDIAG_NBR(2,2,2)=IW
      ENDIF

      ! upper right vertex
      IF (M%IJKW(1,IW)==IBAR+1.AND.M%IJKW(3,IW)==KBAR) THEN
         S%NDIAG_NBR(1,1,3)=M%IJKW(9,IW)
         S%NDIAG_NBR(1,2,3)=IW
      ENDIF
      IF (M%IJKW(1,IW)==IBAR  .AND.M%IJKW(3,IW)==KBAR+1) THEN
         S%NDIAG_NBR(2,1,3)=M%IJKW(9,IW)
         S%NDIAG_NBR(2,2,3)=IW
      ENDIF

      ! upper left vertex
      IF (M%IJKW(1,IW)==1.AND.M%IJKW(3,IW)==KBAR+1) THEN
         S%NDIAG_NBR(1,1,4)=M%IJKW(9,IW)
         S%NDIAG_NBR(1,2,4)=IW
      ENDIF
      IF (M%IJKW(1,IW)==0.AND.M%IJKW(3,IW)==KBAR) THEN
         S%NDIAG_NBR(2,1,4)=M%IJKW(9,IW)
         S%NDIAG_NBR(2,2,4)=IW
      ENDIF

   ENDIF

ENDDO NDIAG_NBR_LOOP2D

IF (SCARC_DEBUG>=2) THEN
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
   WRITE(SCARC_LU,*) 'K_MIN:'
   DO MMM=1,NMESHES
      WRITE(SCARC_LU,'(15i3)') (K_MIN(NNN,MMM),NNN=1,NMESHES)
   ENDDO
   WRITE(SCARC_LU,*) 'K_MAX:'
   DO MMM=1,NMESHES
      WRITE(SCARC_LU,'(15i3)') (K_MAX(NNN,MMM),NNN=1,NMESHES)
   ENDDO
   WRITE(SCARC_LU,*)
   WRITE(SCARC_LU,*) 'NDIAG_NBR(1:2,1,1)=',S%NDIAG_NBR(1,1,1), S%NDIAG_NBR(2,1,1)
   WRITE(SCARC_LU,*) 'NDIAG_NBR(1:2,1,2)=',S%NDIAG_NBR(1,1,2), S%NDIAG_NBR(2,1,2)
   WRITE(SCARC_LU,*) 'NDIAG_NBR(1:2,1,3)=',S%NDIAG_NBR(1,1,3), S%NDIAG_NBR(2,1,3)
   WRITE(SCARC_LU,*) 'NDIAG_NBR(1:2,1,4)=',S%NDIAG_NBR(1,1,4), S%NDIAG_NBR(2,1,4)
   WRITE(SCARC_LU,*)
   WRITE(SCARC_LU,*) 'NDIAG_NBR(1:2,2,1)=',S%NDIAG_NBR(1,2,1), S%NDIAG_NBR(2,2,1)
   WRITE(SCARC_LU,*) 'NDIAG_NBR(1:2,2,2)=',S%NDIAG_NBR(1,2,2), S%NDIAG_NBR(2,2,2)
   WRITE(SCARC_LU,*) 'NDIAG_NBR(1:2,2,3)=',S%NDIAG_NBR(1,2,3), S%NDIAG_NBR(2,2,3)
   WRITE(SCARC_LU,*) 'NDIAG_NBR(1:2,2,4)=',S%NDIAG_NBR(1,2,4), S%NDIAG_NBR(2,2,4)
ENDIF
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Determine arrays IJKW_FACE for finest level
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ILMAX =  S%NLMAX
SLMAX => S%SLEVEL(ILMAX)

ALLOCATE (S%IOR_FACE(NMESHES, NMESHES))
S%IOR_FACE = 0


!---------------------------------------------------------------------
! determine neighbors along faces
!---------------------------------------------------------------------
ALLOCATE (SLMAX%NIC_FACE(NMESHES, NMESHES))
SLMAX%NIC_FACE = NIC

ALLOCATE (SLMAX%I_MIN_FACE(NMESHES, NMESHES))
ALLOCATE (SLMAX%I_MAX_FACE(NMESHES, NMESHES))
ALLOCATE (SLMAX%J_MIN_FACE(NMESHES, NMESHES))
ALLOCATE (SLMAX%J_MAX_FACE(NMESHES, NMESHES))
ALLOCATE (SLMAX%K_MIN_FACE(NMESHES, NMESHES))
ALLOCATE (SLMAX%K_MAX_FACE(NMESHES, NMESHES))

SLMAX%I_MIN_FACE = -10
SLMAX%I_MAX_FACE = -10
SLMAX%J_MIN_FACE = -10
SLMAX%J_MAX_FACE = -10
SLMAX%K_MIN_FACE = -10
SLMAX%K_MAX_FACE = -10

IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,*) 'NM=',NM

SLMAX%NW_FACE = 2*SLMAX%IBAR + 2*SLMAX%KBAR
ALLOCATE (SLMAX%IJKW_FACE(17,SLMAX%NW_FACE))
SLMAX%IJKW_FACE = 0


IF (SCARC_DEBUG>=2) THEN
   WRITE(SCARC_LU,*) 'LEVEL=',S%NLMAX
   WRITE(SCARC_LU,*) 'SHAPE(IJKW_FACE)=',SHAPE(SLMAX%IJKW_FACE)
   WRITE(SCARC_LU,*) 'SIZE (IJKW_FACE)=',SIZE (SLMAX%IJKW_FACE)
   WRITE(SCARC_LU,*) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
   WRITE(SCARC_LU,*) "%%% ARRAY IJKW after all INIT_WALL_CELLS in init.f90"
   WRITE(SCARC_LU,*) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
   WRITE(SCARC_LU,'(17i5)') 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15
   WRITE(SCARC_LU,*) '-------------------------------------------------------------------'
   DO IW=1,M%NWC
     WRITE(SCARC_LU,'(17i5)') IW,(M%IJKW(I,IW),I=1,15)
   ENDDO
   WRITE(SCARC_LU,*) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
ENDIF
   
IW_FACE = 1
IJKW_FACE_LOOP2D: DO IW = 1, M%NWC

   IF (M%IJKW(2,IW) == 1) THEN

      ! copy information from M%IJKW
      DO I= 1, 15
         SLMAX%IJKW_FACE(I,IW_FACE) = M%IJKW(I,IW)
      ENDDO

      ! store number of corresponding IW
      SLMAX%IJKW_FACE(16,IW_FACE) = IW
      SLMAX%IJKW_FACE(17,IW_FACE) = IW

      ! if there is a neighbor at IW, store orientation of neighborship
      NOM=M%IJKW(9,IW)
      IF (NOM/=0) THEN
         S%IOR_FACE(NM,NOM)=M%IJKW(4,IW)
         S%IOR_FACE(NOM,NM)=M%IJKW(4,IW)
      ENDIF

      IW_FACE = IW_FACE+1

   ENDIF

   IF (IW_FACE>SLMAX%NW_FACE) EXIT IJKW_FACE_LOOP2D

ENDDO IJKW_FACE_LOOP2D

! communicate NIC-information from other meshes
ALLOCATE(COUNTS(0:NMESHES-1))
ALLOCATE(COUNTS2D(0:NMESHES-1))
ALLOCATE(DISPLS(0:NMESHES-1))
ALLOCATE(DISPLS2D(0:NMESHES-1))

COUNTS=0
DO N=0,NUMPROCS-1
   DO IM=1,NMESHES
      IF (PROCESS(IM)==N) COUNTS(N)    = COUNTS(N)    + 1
   ENDDO
ENDDO
DISPLS(0)    = 0
DO N=1,NUMPROCS-1
   DISPLS(N)    = COUNTS(N-1)    + DISPLS(N-1)
ENDDO
COUNTS2D      = COUNTS*NMESHES
DISPLS2D      = DISPLS*NMESHES


OTHER_MESH_LOOP2D_MAX: DO NOM=1,NMESHES

    IF (NIC(NM,NOM)==0.AND.NIC(NOM,NM)==0) CYCLE OTHER_MESH_LOOP2D_MAX

    IMIN=1
    IMAX=S%MIBAR(NOM)
    JMIN=1
    JMAX=1
    KMIN=1
    KMAX=S%MKBAR(NOM)

    IF (SCARC_DEBUG.ge.2) THEN
      WRITE(SCARC_LU,*) '============== NOM=',NOM
      WRITE(SCARC_LU,*) 'NIC_FACE:'
      DO JM=1,NMESHES
         WRITE(SCARC_LU, '(10i5)') (SLMAX%NIC_FACE(IM,JM),IM=1,NMESHES)
      ENDDO
      WRITE(SCARC_LU,*) 'MIBAR(NOM)=',S%MIBAR(NOM)
      WRITE(SCARC_LU,*) 'MKBAR(NOM)=',S%MKBAR(NOM)
      WRITE(SCARC_LU,*) 'IREFINE=',IREFINE
      WRITE(SCARC_LU,*) 'IMIN=',IMIN
      WRITE(SCARC_LU,*) 'IMAX=',IMAX
      WRITE(SCARC_LU,*) 'KMIN=',KMIN
      WRITE(SCARC_LU,*) 'KMAX=',KMAX
   ENDIF

   SLMAX%NIC_FACE(NOM,NM) = 0

   SEARCH_LOOP_MAX: DO IW_FACE=1,SLMAX%NW_FACE

      ! neighborship structure already known from finest level
      IF (SLMAX%IJKW_FACE(9,IW_FACE)/=NOM) CYCLE SEARCH_LOOP_MAX

      SLMAX%NIC_FACE(NOM,NM) = SLMAX%NIC_FACE(NOM,NM) + 1
      IOR = SLMAX%IJKW_FACE(4,IW_FACE)

      SELECT CASE(IOR)
         CASE( 1)
            IMIN=MAX(IMIN,SLMAX%IJKW_FACE(10,IW_FACE))
         CASE(-1) 
            IMAX=MIN(IMAX,SLMAX%IJKW_FACE(13,IW_FACE))
         CASE( 2) 
            JMIN=MAX(JMIN,SLMAX%IJKW_FACE(11,IW_FACE))
         CASE(-2) 
            JMAX=MIN(JMAX,SLMAX%IJKW_FACE(14,IW_FACE))
         CASE( 3) 
            KMIN=MAX(KMIN,SLMAX%IJKW_FACE(12,IW_FACE))
         CASE(-3)
            KMAX=MIN(KMAX,SLMAX%IJKW_FACE(15,IW_FACE))
      END SELECT
   ENDDO SEARCH_LOOP_MAX

   IF (BCOM_EXTENDED) THEN                  ! extended communication (including ghost cells)
      IF (ABS(S%IOR_FACE(NOM,NM))==3) THEN         
         SLMAX%I_MIN_FACE(NOM,NM) = IMIN-1
         SLMAX%I_MAX_FACE(NOM,NM) = IMAX+1
      ELSE
         SLMAX%I_MIN_FACE(NOM,NM) = IMIN
         SLMAX%I_MAX_FACE(NOM,NM) = IMAX
      ENDIF
      SLMAX%J_MIN_FACE(NOM,NM) = JMIN
      SLMAX%J_MAX_FACE(NOM,NM) = JMAX
      IF (ABS(S%IOR_FACE(NOM,NM))==1) THEN        
         SLMAX%K_MIN_FACE(NOM,NM) = KMIN-1
         SLMAX%K_MAX_FACE(NOM,NM) = KMAX+1
      ELSE
         SLMAX%K_MIN_FACE(NOM,NM) = KMIN
         SLMAX%K_MAX_FACE(NOM,NM) = KMAX
      ENDIF
   ELSE                                     ! usual communication (excluding ghost cells)
      SLMAX%I_MIN_FACE(NOM,NM) = IMIN                     !
      SLMAX%I_MAX_FACE(NOM,NM) = IMAX
      SLMAX%J_MIN_FACE(NOM,NM) = JMIN
      SLMAX%J_MAX_FACE(NOM,NM) = JMAX
      SLMAX%K_MIN_FACE(NOM,NM) = KMIN
      SLMAX%K_MAX_FACE(NOM,NM) = KMAX
   ENDIF

ENDDO OTHER_MESH_LOOP2D_MAX

IF (SCARC_DEBUG >= 4) THEN
   WRITE(SCARC_LU,*) 'NUMPROCS=',NUMPROCS
   WRITE(SCARC_LU,*) 'MYID  =',MYID  
   WRITE(SCARC_LU,*) 'COUNTS=',COUNTS
   WRITE(SCARC_LU,*) 'COUNTS2D=',COUNTS2D
   WRITE(SCARC_LU,*) 'DISPLS=',DISPLS
   WRITE(SCARC_LU,*) 'DISPLS2D=',DISPLS2D
   WRITE(SCARC_LU,*) 'DISPLS(MYID)+1=',DISPLS(MYID)+1
   WRITE(SCARC_LU,*) 'COUNTS2D(MYID)=',COUNTS2D(MYID)
   WRITE(SCARC_LU,*) 'I_MIN_FACE0:'
   WRITE(SCARC_LU,'(3i4)') ((SLMAX%I_MIN_FACE(IM,JM),JM=1,NMESHES),IM=1,NMESHES)
   WRITE(SCARC_LU,*) 'I_MAX0:'
   WRITE(SCARC_LU,'(3i4)') ((SLMAX%I_MAX_FACE(IM,JM),JM=1,NMESHES),IM=1,NMESHES)
   WRITE(SCARC_LU,*) 'J_MIN0:'
   WRITE(SCARC_LU,'(3i4)') ((SLMAX%J_MIN_FACE(IM,JM),JM=1,NMESHES),IM=1,NMESHES)
   WRITE(SCARC_LU,*) 'J_MAX0:'
   WRITE(SCARC_LU,'(3i4)') ((SLMAX%J_MAX_FACE(IM,JM),JM=1,NMESHES),IM=1,NMESHES)
   WRITE(SCARC_LU,*) 'K_MIN0:'
   WRITE(SCARC_LU,'(3i4)') ((SLMAX%K_MIN_FACE(IM,JM),JM=1,NMESHES),IM=1,NMESHES)
   WRITE(SCARC_LU,*) 'K_MAX0:'
   WRITE(SCARC_LU,'(3i4)') ((SLMAX%K_MAX_FACE(IM,JM),JM=1,NMESHES),IM=1,NMESHES)
ENDIF

IF (USE_MPI) THEN
   CALL MPI_ALLGATHERV(SLMAX%I_MIN_FACE(1,DISPLS(MYID)+1),COUNTS2D(MYID),MPI_INTEGER,&
                       SLMAX%I_MIN_FACE,COUNTS2D,DISPLS2D,MPI_INTEGER,MPI_COMM_WORLD,IERR)
   CALL MPI_ALLGATHERV(SLMAX%I_MAX_FACE(1,DISPLS(MYID)+1),COUNTS2D(MYID),MPI_INTEGER,&
                       SLMAX%I_MAX_FACE,COUNTS2D,DISPLS2D,MPI_INTEGER,MPI_COMM_WORLD,IERR)
   CALL MPI_ALLGATHERV(SLMAX%J_MIN_FACE(1,DISPLS(MYID)+1),COUNTS2D(MYID),MPI_INTEGER,&
                       SLMAX%J_MIN_FACE,COUNTS2D,DISPLS2D,MPI_INTEGER,MPI_COMM_WORLD,IERR)
   CALL MPI_ALLGATHERV(SLMAX%J_MAX_FACE(1,DISPLS(MYID)+1),COUNTS2D(MYID),MPI_INTEGER,&
                       SLMAX%J_MAX_FACE,COUNTS2D,DISPLS2D,MPI_INTEGER,MPI_COMM_WORLD,IERR)
   CALL MPI_ALLGATHERV(SLMAX%K_MIN_FACE(1,DISPLS(MYID)+1),COUNTS2D(MYID),MPI_INTEGER,&
                       SLMAX%K_MIN_FACE,COUNTS2D,DISPLS2D,MPI_INTEGER,MPI_COMM_WORLD,IERR)
   CALL MPI_ALLGATHERV(SLMAX%K_MAX_FACE(1,DISPLS(MYID)+1),COUNTS2D(MYID),MPI_INTEGER,&
                       SLMAX%K_MAX_FACE,COUNTS2D,DISPLS2D,MPI_INTEGER,MPI_COMM_WORLD,IERR)
   CALL MPI_ALLGATHERV(SLMAX%NIC_FACE(1,DISPLS(MYID)+1),  COUNTS2D(MYID),MPI_INTEGER,&
                       SLMAX%NIC_FACE,  COUNTS2D,DISPLS2D,MPI_INTEGER,MPI_COMM_WORLD,IERR)
ELSE
   WRITE(*,*) 'Serial version not yet implemented'
   stop
ENDIF

IF (SCARC_DEBUG >= 4) THEN
   WRITE(SCARC_LU,*) 'AFTER ALLGATHERV'
   WRITE(SCARC_LU,*) 'I_MIN_FACE0:'
   WRITE(SCARC_LU,'(3i4)') ((SLMAX%I_MIN_FACE(IM,JM),JM=1,NMESHES),IM=1,NMESHES)
   WRITE(SCARC_LU,*) 'I_MAX0:'
   WRITE(SCARC_LU,'(3i4)') ((SLMAX%I_MAX_FACE(IM,JM),JM=1,NMESHES),IM=1,NMESHES)
   WRITE(SCARC_LU,*) 'J_MIN0:'
   WRITE(SCARC_LU,'(3i4)') ((SLMAX%J_MIN_FACE(IM,JM),JM=1,NMESHES),IM=1,NMESHES)
   WRITE(SCARC_LU,*) 'J_MAX0:'
   WRITE(SCARC_LU,'(3i4)') ((SLMAX%J_MAX_FACE(IM,JM),JM=1,NMESHES),IM=1,NMESHES)
   WRITE(SCARC_LU,*) 'K_MIN0:'
   WRITE(SCARC_LU,'(3i4)') ((SLMAX%K_MIN_FACE(IM,JM),JM=1,NMESHES),IM=1,NMESHES)
   WRITE(SCARC_LU,*) 'K_MAX0:'
   WRITE(SCARC_LU,'(3i4)') ((SLMAX%K_MAX_FACE(IM,JM),JM=1,NMESHES),IM=1,NMESHES)
ENDIF

SLMAX%I_MIN_FACE = TRANSPOSE(SLMAX%I_MIN_FACE)
SLMAX%I_MAX_FACE = TRANSPOSE(SLMAX%I_MAX_FACE)
SLMAX%J_MIN_FACE = TRANSPOSE(SLMAX%J_MIN_FACE)
SLMAX%J_MAX_FACE = TRANSPOSE(SLMAX%J_MAX_FACE)
SLMAX%K_MIN_FACE = TRANSPOSE(SLMAX%K_MIN_FACE)
SLMAX%K_MAX_FACE = TRANSPOSE(SLMAX%K_MAX_FACE)
SLMAX%NIC_FACE   = TRANSPOSE(SLMAX%NIC_FACE)

IF (SCARC_DEBUG>=4) THEN
   WRITE(SCARC_LU,*) '======================================== ILEVEL ',ILMAX
   WRITE(SCARC_LU,*) 'AFTER TRANSPOSE'
   WRITE(SCARC_LU,*) 'NIC_FACE:'
   WRITE(SCARC_LU,'(2i4)') ((SLMAX%NIC_FACE(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
   WRITE(SCARC_LU,*) 'I_MIN_FACE:'
   WRITE(SCARC_LU,'(3i4)') ((SLMAX%I_MIN_FACE(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
   WRITE(SCARC_LU,*) 'I_MAX:'
   WRITE(SCARC_LU,'(3i4)') ((SLMAX%I_MAX_FACE(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
   WRITE(SCARC_LU,*) 'J_MIN:'
   WRITE(SCARC_LU,'(3i4)') ((SLMAX%J_MIN_FACE(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
   WRITE(SCARC_LU,*) 'J_MAX:'
   WRITE(SCARC_LU,'(3i4)') ((SLMAX%J_MAX_FACE(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
   WRITE(SCARC_LU,*) 'K_MIN:'
   WRITE(SCARC_LU,'(3i4)') ((SLMAX%K_MIN_FACE(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
   WRITE(SCARC_LU,*) 'K_MAX:'
   WRITE(SCARC_LU,'(3i4)') ((SLMAX%K_MAX_FACE(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
   WRITE(SCARC_LU,*) 'IOR_FACE:'
   WRITE(SCARC_LU,'(3i4)') ((S%IOR_FACE(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
   WRITE(SCARC_LU,*) 'SLMAX%NW_FACE=',SLMAX%NW_FACE
   WRITE(SCARC_LU,*) 'SLMAX%IBAR   =',SLMAX%IBAR 
   WRITE(SCARC_LU,*) 'SLMAX%KBAR   =',SLMAX%KBAR 
   WRITE(SCARC_LU,*) 'S%NLEVEL    =',S%NLEVEL
   WRITE(SCARC_LU,*) 'IJKW_FACE:'
   DO IW_FACE=1,SLMAX%NW_FACE
      WRITE(SCARC_LU, '(18i5)') IW_FACE,(SLMAX%IJKW_FACE(I, IW_FACE), I=1, 17)
   ENDDO
   WRITE(SCARC_LU,*) 'NIC_FACE:'
   DO JM=1,NMESHES
      WRITE(SCARC_LU, '(10i5)') (SLMAX%NIC_FACE(IM,JM),IM=1,NMESHES)
   ENDDO
   DO IM=1,NMESHES
      IF (IM/=NM) THEN
         WRITE(SCARC_LU,*) 'FOR MESH ',IM,' ALLOCATING: (',&
                     SLMAX%I_MIN_FACE(NM,IM),':',SLMAX%I_MAX_FACE(NM,IM),' , ',&
                     SLMAX%J_MIN_FACE(NM,IM),':',SLMAX%J_MAX_FACE(NM,IM),' , ',&
                     SLMAX%K_MIN_FACE(NM,IM),':',SLMAX%K_MAX_FACE(NM,IM),')'
      ENDIF
   ENDDO
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Only in case of MG-method:
!!! Determine arrays IJKW_FACE for coarser levels
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (SCARC_METHOD == 'MG' .OR. SCARC_CG_PRECON =='MG') THEN

   IREFINE=1
   INIT_NBR_LEVEL2D: DO ILEVEL=S%NLMAX-1,S%NLMIN,-1
   
      IREFINE=IREFINE*2

      SLHI => S%SLEVEL(ILEVEL+1)
      SLLO => S%SLEVEL(ILEVEL)
   
      SLLO%NW_FACE = 2*SLLO%IBAR + 2*SLLO%KBAR
      ALLOCATE (SLLO%IJKW_FACE(17,SLLO%NW_FACE))
      SLLO%IJKW_FACE=0
   
      IF (SCARC_DEBUG>=2) THEN
         WRITE(SCARC_LU,*) '========== ILEVEL ',ILEVEL
         WRITE(SCARC_LU,*) 'SHAPE(IJKW_FACE)=',SHAPE(SLLO%IJKW_FACE)
         WRITE(SCARC_LU,*) 'SIZE (IJKW_FACE)=',SIZE (SLLO%IJKW_FACE)
         WRITE(SCARC_LU,*) 'IREFINE=',IREFINE
         WRITE(SCARC_LU,*) 'SLLO%NW_FACE=',SLLO%NW_FACE
         WRITE(SCARC_LU,*) 'SLLO%IBAR       =',SLLO%IBAR 
         WRITE(SCARC_LU,*) 'SLLO%KBAR       =',SLLO%KBAR 
      ENDIF
   
      IW_LO=1
   
      !---------------------------------------------------------------------------
      ! cells along IOR=1
      !---------------------------------------------------------------------------
      DO K_LO=1,SLLO%KBAR
   
         SLLO%IJKW_FACE(1,IW_LO)=0           
         SLLO%IJKW_FACE(2,IW_LO)=1           ! ghost cell
         SLLO%IJKW_FACE(3,IW_LO)=K_LO 
         SLLO%IJKW_FACE(4,IW_LO)=1           ! IOR = 1
         SLLO%IJKW_FACE(5,IW_LO)=0           ! currently not used
         SLLO%IJKW_FACE(6,IW_LO)=1           
         SLLO%IJKW_FACE(7,IW_LO)=1           ! boundary cell
         SLLO%IJKW_FACE(8,IW_LO)=K_LO 
   
         IW1_HI=2*K_LO-1
         IW2_HI=2*K_LO
   
         SLLO%IJKW_FACE(16,IW_LO)=IW1_HI
         SLLO%IJKW_FACE(17,IW_LO)=IW2_HI
         
         NOM_HI1=SLHI%IJKW_FACE(9,IW1_HI)
         NOM_HI2=SLHI%IJKW_FACE(9,IW2_HI)
   
         IF (NOM_HI1 == NOM_HI2) THEN

            SLLO%IJKW_FACE(9,IW_LO)=NOM_HI1 

            IF (NOM_HI1 > 0) THEN

               IMIN1_HI=SLHI%IJKW_FACE(10,IW1_HI)
               IMAX1_HI=SLHI%IJKW_FACE(13,IW1_HI)
         
               IMIN2_HI=SLHI%IJKW_FACE(10,IW2_HI)
               IMAX2_HI=SLHI%IJKW_FACE(13,IW2_HI)
         
               IF (IMIN1_HI/=IMIN2_HI.OR.IMAX2_HI/=IMAX2_HI) THEN
                  WRITE(*,*) '1: Wrong relation of IMIN?_HI and IMAX?_HI for IOR=1'
                  STOP
               ENDIF
      
               KMIN1_HI=SLHI%IJKW_FACE(12,IW1_HI)
               KMAX1_HI=SLHI%IJKW_FACE(15,IW1_HI)
         
               KMIN2_HI=SLHI%IJKW_FACE(12,IW2_HI)
               KMAX2_HI=SLHI%IJKW_FACE(15,IW2_HI)
         
               IBAR_NOM=S%MIBAR(NOM_HI1)/IREFINE
               
               ! same KBAR-resolution on current mesh and on neighbor
               IF ((KMAX2_HI-KMAX1_HI)==1) THEN
                  SLLO%IJKW_FACE(10,IW_LO)=IBAR_NOM
                  SLLO%IJKW_FACE(11,IW_LO)=1
                  SLLO%IJKW_FACE(12,IW_LO)=K_LO
                  SLLO%IJKW_FACE(13,IW_LO)=IBAR_NOM
                  SLLO%IJKW_FACE(14,IW_LO)=1
                  SLLO%IJKW_FACE(15,IW_LO)=K_LO
      
               ! KBAR-resolution coarser on current mesh than on neighbor
               ELSE IF ((KMAX2_HI-KMAX1_HI)==2) THEN
                  SLLO%IJKW_FACE(10,IW_LO)=IBAR_NOM
                  SLLO%IJKW_FACE(11,IW_LO)=1
                  SLLO%IJKW_FACE(12,IW_LO)=2*K_LO-1
                  SLLO%IJKW_FACE(13,IW_LO)=IBAR_NOM
                  SLLO%IJKW_FACE(14,IW_LO)=1
                  SLLO%IJKW_FACE(15,IW_LO)=2*K_LO
      
               ! KBAR-resolution finer on current mesh than on neighbor
               ELSE IF ((KMAX2_HI-KMAX1_HI)==0) THEN
                  SLLO%IJKW_FACE(10,IW_LO)=IBAR_NOM
                  SLLO%IJKW_FACE(11,IW_LO)=1
                  SLLO%IJKW_FACE(12,IW_LO)=MOD(K_LO+1,2)
                  SLLO%IJKW_FACE(13,IW_LO)=IBAR_NOM
                  SLLO%IJKW_FACE(14,IW_LO)=1
                  SLLO%IJKW_FACE(15,IW_LO)=MOD(K_LO+1,2)
      
               ELSE                                         
                  WRITE(*,*) 'Error in INIT_NBR_LEVEL2D at IOR=1 for ILEVEL=',ILEVEL
                  WRITE(SCARC_LU,*) 'Error in INIT_NBR_LEVEL2D at IOR=1 for ILEVEL=',ILEVEL
                  STOP
               ENDIF
         
               IF (SCARC_DEBUG.ge.2) THEN
                  WRITE(SCARC_LU,*) '==== K_LO=',K_LO,' ============================='
                  WRITE(SCARC_LU,*) 'NOM_HI1 =',NOM_HI1
                  WRITE(SCARC_LU,*) 'NOM_HI2 =',NOM_HI2
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'IW1_HI  =',IW1_HI
                  WRITE(SCARC_LU,*) 'IW2_HI  =',IW2_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'IMIN1_HI=',IMIN1_HI
                  WRITE(SCARC_LU,*) 'IMAX1_HI=',IMAX1_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'IMIN2_HI=',IMIN2_HI
                  WRITE(SCARC_LU,*) 'IMAX2_HI=',IMAX2_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'KMIN1_HI=',KMIN1_HI
                  WRITE(SCARC_LU,*) 'KMAX1_HI=',KMAX1_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'KMIN2_HI=',KMIN2_HI
                  WRITE(SCARC_LU,*) 'KMAX2_HI=',KMAX2_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU, '(i5,a,17i5)') IW_LO,':',(SLLO%IJKW_FACE(I, IW_LO), I=1,17)
               ENDIF
            ENDIF

         ELSE

            WRITE(*,*) 'INCONSISTENT NEIGHBORS ON IOR=1, not yet implemented!'
            STOP

         ENDIF
   
         IW_LO=IW_LO+1
   
      ENDDO

   
      !---------------------------------------------------------------------------
      ! cells along IOR=-1
      !---------------------------------------------------------------------------
      DO K_LO=1,SLLO%KBAR
   
         SLLO%IJKW_FACE(1,IW_LO)=SLLO%IBAR+1 
         SLLO%IJKW_FACE(2,IW_LO)=1 
         SLLO%IJKW_FACE(3,IW_LO)=K_LO
         SLLO%IJKW_FACE(4,IW_LO)=-1 
         SLLO%IJKW_FACE(6,IW_LO)=SLLO%IBAR 
         SLLO%IJKW_FACE(7,IW_LO)=1 
         SLLO%IJKW_FACE(8,IW_LO)=K_LO
   
         IW1_HI=SLHI%KBAR+2*K_LO-1
         IW2_HI=SLHI%KBAR+2*K_LO
   
         SLLO%IJKW_FACE(16,IW_LO)=IW1_HI
         SLLO%IJKW_FACE(17,IW_LO)=IW2_HI
         
         NOM_HI1=SLHI%IJKW_FACE(9,IW1_HI)
         NOM_HI2=SLHI%IJKW_FACE(9,IW2_HI)
   
         IF (NOM_HI1 == NOM_HI2) THEN
   
            SLLO%IJKW_FACE(9,IW_LO)=NOM_HI1 
 
            IF (NOM_HI1 > 0) THEN
   
               IMIN1_HI=SLHI%IJKW_FACE(10,IW1_HI)
               IMAX1_HI=SLHI%IJKW_FACE(13,IW1_HI)
         
               IMIN2_HI=SLHI%IJKW_FACE(10,IW2_HI)
               IMAX2_HI=SLHI%IJKW_FACE(13,IW2_HI)
         
               IF (IMIN1_HI/=IMIN2_HI.OR.IMAX2_HI/=IMAX2_HI) THEN
                  WRITE(*,*) '1: Wrong relation of IMIN?_HI and IMAX?_HI for IOR=-1'
                  STOP
               ENDIF
               
               KMIN1_HI=SLHI%IJKW_FACE(12,IW1_HI)
               KMAX1_HI=SLHI%IJKW_FACE(15,IW1_HI)
         
               KMIN2_HI=SLHI%IJKW_FACE(12,IW2_HI)
               KMAX2_HI=SLHI%IJKW_FACE(15,IW2_HI)
         
               ! same KBAR-resolution on current mesh and on neighbor
               IF ((KMAX2_HI-KMAX1_HI)==1) THEN
                  SLLO%IJKW_FACE(10,IW_LO)=1
                  SLLO%IJKW_FACE(11,IW_LO)=1
                  SLLO%IJKW_FACE(12,IW_LO)=K_LO
                  SLLO%IJKW_FACE(13,IW_LO)=1
                  SLLO%IJKW_FACE(14,IW_LO)=1
                  SLLO%IJKW_FACE(15,IW_LO)=K_LO
      
               ! KBAR-resolution coarser on current mesh than on neighbor
               ELSE IF ((KMAX2_HI-KMAX1_HI)==2) THEN
                  SLLO%IJKW_FACE(10,IW_LO)=1
                  SLLO%IJKW_FACE(11,IW_LO)=1
                  SLLO%IJKW_FACE(12,IW_LO)=2*K_LO-1
                  SLLO%IJKW_FACE(13,IW_LO)=1
                  SLLO%IJKW_FACE(14,IW_LO)=1
                  SLLO%IJKW_FACE(15,IW_LO)=2*K_LO
      
               ! KBAR-resolution finer on current mesh than on neighbor
               ELSE IF ((KMAX2_HI-KMAX1_HI)==0) THEN
                  SLLO%IJKW_FACE(10,IW_LO)=1
                  SLLO%IJKW_FACE(11,IW_LO)=1
                  SLLO%IJKW_FACE(12,IW_LO)=MOD(K_LO+1,2)
                  SLLO%IJKW_FACE(13,IW_LO)=1
                  SLLO%IJKW_FACE(14,IW_LO)=1
                  SLLO%IJKW_FACE(15,IW_LO)=MOD(K_LO+1,2)
      
               ELSE                                         
                  WRITE(*,*) 'Error in INIT_NBR_LEVEL2D at IOR=-1 for ILEVEL=',ILEVEL
                  WRITE(SCARC_LU,*) 'Error in INIT_NBR_LEVEL2D at IOR=-1 for ILEVEL=',ILEVEL
                  STOP
               ENDIF
         
               IF (SCARC_DEBUG.ge.2) THEN
                  WRITE(SCARC_LU,*) '==== K_LO=',K_LO,' ============================='
                  WRITE(SCARC_LU,*) 'NOM_HI1 =',NOM_HI1
                  WRITE(SCARC_LU,*) 'NOM_HI2 =',NOM_HI2
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'IW1_HI  =',IW1_HI
                  WRITE(SCARC_LU,*) 'IW2_HI  =',IW2_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'IMIN1_HI=',IMIN1_HI
                  WRITE(SCARC_LU,*) 'IMAX1_HI=',IMAX1_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'IMIN2_HI=',IMIN2_HI
                  WRITE(SCARC_LU,*) 'IMAX2_HI=',IMAX2_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'KMIN1_HI=',KMIN1_HI
                  WRITE(SCARC_LU,*) 'KMAX1_HI=',KMAX1_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'KMIN2_HI=',KMIN2_HI
                  WRITE(SCARC_LU,*) 'KMAX2_HI=',KMAX2_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU, '(i5,a,17i5)') IW_LO,':',(SLLO%IJKW_FACE(I, IW_LO), I=1,17)
               ENDIF
            ENDIF
         ELSE

            WRITE(*,*) 'INCONSISTENT NEIGHBORS ON IOR=-1, not yet implemented!'
            STOP

         ENDIF
   
         IW_LO=IW_LO+1

      ENDDO

   
      !---------------------------------------------------------------------------
      ! cells along IOR=3
      !---------------------------------------------------------------------------
      DO I_LO=1,SLLO%IBAR
   
         SLLO%IJKW_FACE(1,IW_LO)=I_LO
         SLLO%IJKW_FACE(2,IW_LO)=1 
         SLLO%IJKW_FACE(3,IW_LO)=0 
         SLLO%IJKW_FACE(4,IW_LO)=3 
         SLLO%IJKW_FACE(6,IW_LO)=I_LO
         SLLO%IJKW_FACE(7,IW_LO)=1 
         SLLO%IJKW_FACE(8,IW_LO)=1 
   
         IW1_HI=2*SLHI%KBAR+2*I_LO-1
         IW2_HI=2*SLHI%KBAR+2*I_LO

         SLLO%IJKW_FACE(16,IW_LO)=IW1_HI
         SLLO%IJKW_FACE(17,IW_LO)=IW2_HI
         
         NOM_HI1=SLHI%IJKW_FACE(9,IW1_HI)
         NOM_HI2=SLHI%IJKW_FACE(9,IW2_HI)
   
         IF (NOM_HI1 == NOM_HI2) THEN

            SLLO%IJKW_FACE(9,IW_LO)=NOM_HI1 

            IF (NOM_HI1 > 0) THEN

               IMIN1_HI=SLHI%IJKW_FACE(10,IW1_HI)
               IMAX1_HI=SLHI%IJKW_FACE(13,IW1_HI)
         
               IMIN2_HI=SLHI%IJKW_FACE(10,IW2_HI)
               IMAX2_HI=SLHI%IJKW_FACE(13,IW2_HI)
         
               KMIN1_HI=SLHI%IJKW_FACE(12,IW1_HI)
               KMAX1_HI=SLHI%IJKW_FACE(15,IW1_HI)
         
               KMIN2_HI=SLHI%IJKW_FACE(12,IW2_HI)
               KMAX2_HI=SLHI%IJKW_FACE(15,IW2_HI)
         
               IF (KMIN1_HI/=KMIN2_HI.OR.KMAX2_HI/=KMAX2_HI) THEN
                  WRITE(*,*) '1: Wrong relation of KMIN?_HI and KMAX?_HI for IOR=3'
                  STOP
               ENDIF

               KBAR_NOM=S%MKBAR(NOM_HI1)/IREFINE
               
               ! same KBAR-resolution on current mesh and on neighbor
               IF ((IMAX2_HI-IMAX1_HI)==1) THEN
                  SLLO%IJKW_FACE(10,IW_LO)=I_LO
                  SLLO%IJKW_FACE(11,IW_LO)=1
                  SLLO%IJKW_FACE(12,IW_LO)=KBAR_NOM
                  SLLO%IJKW_FACE(13,IW_LO)=I_LO
                  SLLO%IJKW_FACE(14,IW_LO)=1
                  SLLO%IJKW_FACE(15,IW_LO)=KBAR_NOM
      
               ! KBAR-resolution coarser on current mesh than on neighbor
               ELSE IF ((IMAX2_HI-IMAX1_HI)==2) THEN
                  SLLO%IJKW_FACE(10,IW_LO)=2*I_LO-1
                  SLLO%IJKW_FACE(11,IW_LO)=1
                  SLLO%IJKW_FACE(12,IW_LO)=KBAR_NOM
                  SLLO%IJKW_FACE(13,IW_LO)=2*I_LO
                  SLLO%IJKW_FACE(14,IW_LO)=1
                  SLLO%IJKW_FACE(15,IW_LO)=KBAR_NOM
      
               ! KBAR-resolution finer on current mesh than on neighbor
               ELSE IF ((IMAX2_HI-IMAX1_HI)==0) THEN
                  SLLO%IJKW_FACE(10,IW_LO)=MOD(I_LO+1,2)
                  SLLO%IJKW_FACE(11,IW_LO)=1
                  SLLO%IJKW_FACE(12,IW_LO)=KBAR_NOM
                  SLLO%IJKW_FACE(13,IW_LO)=MOD(I_LO+1,2)
                  SLLO%IJKW_FACE(14,IW_LO)=1
                  SLLO%IJKW_FACE(15,IW_LO)=KBAR_NOM
      
               ELSE                                         
                  WRITE(*,*) 'Error in INIT_NBR_LEVEL2D at IOR=3 for ILEVEL=',ILEVEL
                  WRITE(SCARC_LU,*) 'Error in INIT_NBR_LEVEL2D at IOR=3 for ILEVEL=',ILEVEL
                  STOP
               ENDIF
         
               IF (SCARC_DEBUG.ge.2) THEN
                  WRITE(SCARC_LU,*) '==== I_LO=',I_LO,' ============================='
                  WRITE(SCARC_LU,*) 'NOM_HI1 =',NOM_HI1
                  WRITE(SCARC_LU,*) 'NOM_HI2 =',NOM_HI2
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'IW1_HI  =',IW1_HI
                  WRITE(SCARC_LU,*) 'IW2_HI  =',IW2_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'IMIN1_HI=',IMIN1_HI
                  WRITE(SCARC_LU,*) 'IMAX1_HI=',IMAX1_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'IMIN2_HI=',IMIN2_HI
                  WRITE(SCARC_LU,*) 'IMAX2_HI=',IMAX2_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'KMIN1_HI=',KMIN1_HI
                  WRITE(SCARC_LU,*) 'KMAX1_HI=',KMAX1_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'KMIN2_HI=',KMIN2_HI
                  WRITE(SCARC_LU,*) 'KMAX2_HI=',KMAX2_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU, '(i5,a,17i5)') IW_LO,':',(SLLO%IJKW_FACE(I, IW_LO), I=1,17)
               ENDIF
            ENDIF

         ELSE

            WRITE(*,*) 'INCONSISTENT NEIGHBORS ON IOR=3, not yet implemented!'
            STOP

         ENDIF
   
         IW_LO=IW_LO+1
   
      ENDDO
   
      !---------------------------------------------------------------------------
      ! cells along IOR=-3
      !---------------------------------------------------------------------------
      DO I_LO=1,SLLO%IBAR
   
         SLLO%IJKW_FACE(1,IW_LO)=I_LO
         SLLO%IJKW_FACE(2,IW_LO)=1 
         SLLO%IJKW_FACE(3,IW_LO)=SLLO%KBAR+1 
         SLLO%IJKW_FACE(4,IW_LO)=-3 
         SLLO%IJKW_FACE(6,IW_LO)=I_LO
         SLLO%IJKW_FACE(7,IW_LO)=1 
         SLLO%IJKW_FACE(8,IW_LO)=SLLO%KBAR 
   
         IW1_HI=2*SLHI%KBAR+SLHI%IBAR+2*I_LO-1
         IW2_HI=2*SLHI%KBAR+SLHI%IBAR+2*I_LO

         SLLO%IJKW_FACE(16,IW_LO)=IW1_HI
         SLLO%IJKW_FACE(17,IW_LO)=IW2_HI
         
         NOM_HI1=SLHI%IJKW_FACE(9,IW1_HI)
         NOM_HI2=SLHI%IJKW_FACE(9,IW2_HI)
   
         IF (NOM_HI1 == NOM_HI2) THEN

            SLLO%IJKW_FACE(9,IW_LO)=NOM_HI1 
           
            IF (NOM_HI1 > 0) THEN

               IMIN1_HI=SLHI%IJKW_FACE(10,IW1_HI)
               IMAX1_HI=SLHI%IJKW_FACE(13,IW1_HI)
         
               IMIN2_HI=SLHI%IJKW_FACE(10,IW2_HI)
               IMAX2_HI=SLHI%IJKW_FACE(13,IW2_HI)
      
               KMIN1_HI=SLHI%IJKW_FACE(12,IW1_HI)
               KMAX1_HI=SLHI%IJKW_FACE(15,IW1_HI)
         
               KMIN2_HI=SLHI%IJKW_FACE(12,IW2_HI)
               KMAX2_HI=SLHI%IJKW_FACE(15,IW2_HI)
         
               IF (KMIN1_HI/=KMIN2_HI.OR.KMAX2_HI/=KMAX2_HI) THEN
                  WRITE(*,*) '1: Wrong relation of KMIN?_HI and KMAX?_HI for IOR=3'
                  STOP
               ENDIF
               
               ! same KBAR-resolution on current mesh and on neighbor
               IF ((IMAX2_HI-IMAX1_HI)==1) THEN
                  SLLO%IJKW_FACE(10,IW_LO)=I_LO
                  SLLO%IJKW_FACE(11,IW_LO)=1
                  SLLO%IJKW_FACE(12,IW_LO)=1
                  SLLO%IJKW_FACE(13,IW_LO)=I_LO
                  SLLO%IJKW_FACE(14,IW_LO)=1
                  SLLO%IJKW_FACE(15,IW_LO)=1
      
               ! KBAR-resolution coarser on current mesh than on neighbor
               ELSE IF ((IMAX2_HI-IMAX1_HI)==2) THEN
                  SLLO%IJKW_FACE(10,IW_LO)=2*I_LO-1
                  SLLO%IJKW_FACE(11,IW_LO)=1
                  SLLO%IJKW_FACE(12,IW_LO)=1
                  SLLO%IJKW_FACE(13,IW_LO)=2*I_LO
                  SLLO%IJKW_FACE(14,IW_LO)=1
                  SLLO%IJKW_FACE(15,IW_LO)=1
      
               ! KBAR-resolution finer on current mesh than on neighbor
               ELSE IF ((IMAX2_HI-IMAX1_HI)==0) THEN
                  SLLO%IJKW_FACE(10,IW_LO)=MOD(I_LO+1,2)
                  SLLO%IJKW_FACE(11,IW_LO)=1
                  SLLO%IJKW_FACE(12,IW_LO)=1
                  SLLO%IJKW_FACE(13,IW_LO)=MOD(I_LO+1,2)
                  SLLO%IJKW_FACE(14,IW_LO)=1
                  SLLO%IJKW_FACE(15,IW_LO)=1
      
               ELSE                                         
                  WRITE(*,*) 'Error in INIT_NBR_LEVEL2D at IOR=-3 for ILEVEL=',ILEVEL
                  WRITE(SCARC_LU,*) 'Error in INIT_NBR_LEVEL2D at IOR=-3 for ILEVEL=',ILEVEL
                  STOP
               ENDIF
         
               IF (SCARC_DEBUG.ge.2) THEN
                  WRITE(SCARC_LU,*) '==== I_LO=',I_LO,' ============================='
                  WRITE(SCARC_LU,*) 'NOM_HI1 =',NOM_HI1
                  WRITE(SCARC_LU,*) 'NOM_HI2 =',NOM_HI2
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'IW1_HI  =',IW1_HI
                  WRITE(SCARC_LU,*) 'IW2_HI  =',IW2_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'IMIN1_HI=',IMIN1_HI
                  WRITE(SCARC_LU,*) 'IMAX1_HI=',IMAX1_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'IMIN2_HI=',IMIN2_HI
                  WRITE(SCARC_LU,*) 'IMAX2_HI=',IMAX2_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'KMIN1_HI=',KMIN1_HI
                  WRITE(SCARC_LU,*) 'KMAX1_HI=',KMAX1_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'KMIN2_HI=',KMIN2_HI
                  WRITE(SCARC_LU,*) 'KMAX2_HI=',KMAX2_HI
                  WRITE(SCARC_LU,*) 
                  WRITE(SCARC_LU,*) 'Real-Kmax2_HI=',REAL(KMAX2_HI/2,EB)
                  WRITE(SCARC_LU,*) 'Ceil-Kmax2_HI=',CEILING(REAL(KMAX2_HI,EB)/2)
                  WRITE(SCARC_LU, '(i5,a,17i5)') IW_LO,':',(SLLO%IJKW_FACE(I, IW_LO), I=1,17)
               ENDIF
            ENDIF

         ELSE

            WRITE(*,*) 'INCONSISTENT NEIGHBORS ON IOR=-3, not yet implemented!'
            STOP

         ENDIF
   
         IW_LO=IW_LO+1
   
      ENDDO
   
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! compute analogues to NIC_FACE, I_MIN_FACE, I_MAX, K_MIN, K_MAX on level 'ILEVEL'
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ALLOCATE (SLLO%NIC_FACE(NMESHES, NMESHES))
      SLLO%NIC_FACE=0

      ALLOCATE (SLLO%I_MIN_FACE(NMESHES, NMESHES))
      ALLOCATE (SLLO%I_MAX_FACE(NMESHES, NMESHES))
      ALLOCATE (SLLO%J_MIN_FACE(NMESHES, NMESHES))
      ALLOCATE (SLLO%J_MAX_FACE(NMESHES, NMESHES))
      ALLOCATE (SLLO%K_MIN_FACE(NMESHES, NMESHES))
      ALLOCATE (SLLO%K_MAX_FACE(NMESHES, NMESHES))
      SLLO%I_MIN_FACE = 0
      SLLO%I_MAX_FACE = 0
      SLLO%J_MIN_FACE = 0
      SLLO%J_MAX_FACE = 0
      SLLO%K_MIN_FACE = 0
      SLLO%K_MAX_FACE = 0


      OTHER_MESH_LOOP2D: DO NOM=1,NMESHES
 
         IF (SLMAX%NIC_FACE(NM,NOM)==0.AND.SLMAX%NIC_FACE(NOM,NM)==0) CYCLE OTHER_MESH_LOOP2D
 
         !!! ACHTUNG: funktioniert nur für 2er-Potenz-Gitterweiten !!!!!
         IMIN=0
         IMAX=S%MIBAR(NOM)/IREFINE+1
         JMIN=0
         JMAX=2
         KMIN=0
         KMAX=S%MKBAR(NOM)/IREFINE+1

         IF (SCARC_DEBUG.ge.6) THEN
            WRITE(SCARC_LU,*) '============== NOM=',NOM
            WRITE(SCARC_LU,*) 'NIC_FACE:'
            DO JM=1,NMESHES
               WRITE(SCARC_LU, '(10i5)') (SLMAX%NIC_FACE(IM,JM),IM=1,NMESHES)
            ENDDO
            WRITE(SCARC_LU,*) 'MIBAR(NOM)=',S%MIBAR(NOM)
            WRITE(SCARC_LU,*) 'MKBAR(NOM)=',S%MKBAR(NOM)
            WRITE(SCARC_LU,*) 'IREFINE=',IREFINE
            WRITE(SCARC_LU,*) 'IMAX=',IMAX
            WRITE(SCARC_LU,*) 'KMAX=',KMAX
         ENDIF

         SLLO%NIC_FACE(NOM,NM) = 0

         SEARCH_LOOP: DO IW_FACE=1,SLLO%NW_FACE

            ! neighborship structure already known from finest level
            IF (SLLO%IJKW_FACE(9,IW_FACE)/=NOM) CYCLE SEARCH_LOOP

            !WRITE(SCARC_LU,*) 'IJKW_FACE(9,',IW_FACE,')=',SLLO%IJKW_FACE(9,IW_FACE)

            SLLO%NIC_FACE(NOM,NM) = SLLO%NIC_FACE(NOM,NM) + 1
            IOR = SLLO%IJKW_FACE(4,IW_FACE)

            !WRITE(SCARC_LU,*) 'NIC_FACE(',NOM,',',NM,')=',SLLO%NIC_FACE(NOM,NM)

            SELECT CASE(IOR)
               CASE( 1)
                  IMIN=MAX(IMIN,SLLO%IJKW_FACE(10,IW_FACE)-1)
               CASE(-1) 
                  IMAX=MIN(IMAX,SLLO%IJKW_FACE(13,IW_FACE))
               !CASE( 2) 
               !   JMIN=MAX(JMIN,SLLO%IJKW_FACE(11,IW_FACE)-1)
               !CASE(-2) 
               !   JMAX=MIN(JMAX,SLLO%IJKW_FACE(14,IW_FACE))
               CASE( 3) 
                  KMIN=MAX(KMIN,SLLO%IJKW_FACE(12,IW_FACE)-1)
               CASE(-3)
                  KMAX=MIN(KMAX,SLLO%IJKW_FACE(15,IW_FACE))
            END SELECT
         ENDDO SEARCH_LOOP

         !IF (ABS(S%IOR_FACE(NOM,NM))==3) THEN               ! extended communication
         !    SLMAX%I_MIN_FACE(NOM,NM) = IMIN-1
         !    SLMAX%I_MAX_FACE(NOM,NM) = IMAX+1
         !ELSE
         !   SLMAX%I_MIN_FACE(NOM,NM) = IMIN
         !   SLMAX%I_MAX_FACE(NOM,NM) = IMAX
         !ENDIF
         !SLMAX%J_MIN_FACE(NOM,NM) = JMIN
         !SLMAX%J_MAX_FACE(NOM,NM) = JMAX
         !IF (ABS(S%IOR_FACE(NOM,NM))==1) THEN               ! extended communication
         !   SLMAX%K_MIN_FACE(NOM,NM) = KMIN-1
         !   SLMAX%K_MAX_FACE(NOM,NM) = KMAX+1
         !ELSE
         !   SLMAX%K_MIN_FACE(NOM,NM) = KMIN
         !   SLMAX%K_MAX_FACE(NOM,NM) = KMAX
         !ENDIF

         SLLO%I_MIN_FACE(NOM,NM) = IMIN
         SLLO%I_MAX_FACE(NOM,NM) = IMAX
         SLLO%J_MIN_FACE(NOM,NM) = JMIN
         SLLO%J_MAX_FACE(NOM,NM) = JMAX
         SLLO%K_MIN_FACE(NOM,NM) = KMIN
         SLLO%K_MAX_FACE(NOM,NM) = KMAX

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
         WRITE(SCARC_LU,*) 'I_MIN_FACE(1,DISPLS(MYID)+1:DISPLS(MYID)+3)=',SLLO%I_MIN_FACE(1,DISPLS(MYID)+1:DISPLS(MYID)+3)
         WRITE(SCARC_LU,*) '=======================A:'
         DO JM=1,NMESHES
            WRITE(SCARC_LU, '(10i5)') (SLLO%I_MIN_FACE(IM,JM),IM=1,NMESHES)
         ENDDO
      ENDIF

      IF (USE_MPI) THEN
         CALL MPI_ALLGATHERV(SLLO%I_MIN_FACE(1,DISPLS(MYID)+1),COUNTS2D(MYID),MPI_INTEGER,&
                             SLLO%I_MIN_FACE,COUNTS2D,DISPLS2D,MPI_INTEGER,MPI_COMM_WORLD,IERR)
         CALL MPI_ALLGATHERV(SLLO%I_MAX_FACE(1,DISPLS(MYID)+1),COUNTS2D(MYID),MPI_INTEGER,&
                             SLLO%I_MAX_FACE,COUNTS2D,DISPLS2D,MPI_INTEGER,MPI_COMM_WORLD,IERR)
         CALL MPI_ALLGATHERV(SLLO%J_MIN_FACE(1,DISPLS(MYID)+1),COUNTS2D(MYID),MPI_INTEGER,&
                             SLLO%J_MIN_FACE,COUNTS2D,DISPLS2D,MPI_INTEGER,MPI_COMM_WORLD,IERR)
         CALL MPI_ALLGATHERV(SLLO%J_MAX_FACE(1,DISPLS(MYID)+1),COUNTS2D(MYID),MPI_INTEGER,&
                             SLLO%J_MAX_FACE,COUNTS2D,DISPLS2D,MPI_INTEGER,MPI_COMM_WORLD,IERR)
         CALL MPI_ALLGATHERV(SLLO%K_MIN_FACE(1,DISPLS(MYID)+1),COUNTS2D(MYID),MPI_INTEGER,&
                             SLLO%K_MIN_FACE,COUNTS2D,DISPLS2D,MPI_INTEGER,MPI_COMM_WORLD,IERR)
         CALL MPI_ALLGATHERV(SLLO%K_MAX_FACE(1,DISPLS(MYID)+1),COUNTS2D(MYID),MPI_INTEGER,&
                             SLLO%K_MAX_FACE,COUNTS2D,DISPLS2D,MPI_INTEGER,MPI_COMM_WORLD,IERR)
         CALL MPI_ALLGATHERV(SLLO%NIC_FACE(1,DISPLS(MYID)+1),  COUNTS2D(MYID),MPI_INTEGER,&
                             SLLO%NIC_FACE,  COUNTS2D,DISPLS2D,MPI_INTEGER,MPI_COMM_WORLD,IERR)
      ELSE
         WRITE(*,*) 'Serial version not yet implemented'
         stop
      ENDIF

      SLLO%I_MIN_FACE = TRANSPOSE(SLLO%I_MIN_FACE)
      SLLO%I_MAX_FACE = TRANSPOSE(SLLO%I_MAX_FACE)
      SLLO%J_MIN_FACE = TRANSPOSE(SLLO%J_MIN_FACE)
      SLLO%J_MAX_FACE = TRANSPOSE(SLLO%J_MAX_FACE)
      SLLO%K_MIN_FACE = TRANSPOSE(SLLO%K_MIN_FACE)
      SLLO%K_MAX_FACE = TRANSPOSE(SLLO%K_MAX_FACE)
      SLLO%NIC_FACE   = TRANSPOSE(SLLO%NIC_FACE)

      IF (SCARC_DEBUG >= 6) THEN
         WRITE(SCARC_LU,*) '=================================================='
         WRITE(SCARC_LU,*) '================== LEVEL =',ILEVEL
         WRITE(SCARC_LU,*) 'NFACE      =',S%NFACE0
         WRITE(SCARC_LU,*) 'NW_FACE=',SLLO%NW_FACE
         WRITE(SCARC_LU,*) 'IJKW_FACE:'
         DO IW_FACE=1,SLLO%NW_FACE
            WRITE(SCARC_LU, '(i5,a,17i5)') IW_FACE,':',(SLLO%IJKW_FACE(I, IW_FACE), I=1,17)
         ENDDO
         WRITE(SCARC_LU,*) 'NIC_FACE:'
         WRITE(SCARC_LU, '(3i4)') ((SLLO%NIC_FACE(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
         WRITE(SCARC_LU,*) 'I_MIN_FACE:'
         WRITE(SCARC_LU, '(3i4)') ((SLLO%I_MIN_FACE(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
         WRITE(SCARC_LU,*) 'I_MAX:'
         WRITE(SCARC_LU, '(3i4)') ((SLLO%I_MAX_FACE(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
         WRITE(SCARC_LU,*) 'J_MIN:'
         WRITE(SCARC_LU, '(3i4)') ((SLLO%J_MIN_FACE(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
         WRITE(SCARC_LU,*) 'J_MAX:'
         WRITE(SCARC_LU, '(3i4)') ((SLLO%J_MAX_FACE(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
         WRITE(SCARC_LU,*) 'K_MIN:'
         WRITE(SCARC_LU, '(3i4)') ((SLLO%K_MIN_FACE(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
         WRITE(SCARC_LU,*) 'K_MAX:'
         WRITE(SCARC_LU, '(3i4)') ((SLLO%K_MAX_FACE(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
         DO IM=1,NMESHES
            IF (IM/=NM) THEN
               WRITE(SCARC_LU,*) 'FOR MESH ',IM,' ALLOCATING: (',&
                           SLLO%I_MIN_FACE(NM,IM),':',SLLO%I_MAX_FACE(NM,IM),' , ',&
                           SLLO%J_MIN_FACE(NM,IM),':',SLLO%J_MAX_FACE(NM,IM),' , ',&
                           SLLO%K_MIN_FACE(NM,IM),':',SLLO%K_MAX_FACE(NM,IM),')'
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

   IF (SLMAX%NIC_FACE(NM,NOM)/=0.AND.SLMAX%NIC_FACE(NOM,NM)/=0) THEN

      OS => S%OSCARC(NOM)
      ALLOCATE (OS%SLEVEL(S%NLMIN:S%NLMAX), STAT=IERR)
      CALL CHKMEMERR ('SCARC_INIT_COMMUNICATION2D', 'OS%SLEVEL', IERR)

      OSLMAX => OS%SLEVEL(ILMAX)
      OSLMAX%NW_FACE = 2*S%MIBAR(NOM) + 2*S%MKBAR(NOM)
      ALLOCATE (OSLMAX%IJKW_FACE(17,OSLMAX%NW_FACE), STAT=IERR)
      CALL CHKMEMERR ('SCARC_INIT_COMMUNICATION2D', 'OSLMAX%IJKW_FACE', IERR)
      OSLMAX%IJKW_FACE = 0

      IF (SCARC_DEBUG>=2) THEN
         WRITE(SCARC_LU,*) 'OSLMAX zeigt auf OS%SLEVEL(',ILMAX,')'
         WRITE(SCARC_LU,*) 'OSLMAX%NW_FACE=',OSLMAX%NW_FACE
         WRITE(SCARC_LU,*) 'SIZE(OSLMAX%IJKW_FACE)=',SIZE(OSLMAX%IJKW_FACE)
         WRITE(SCARC_LU,*) 'S%NLDIFF=',S%NLDIFF
      ENDIF

      IF (S%NLDIFF/=0) THEN
         IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,*) 'S%NLDIFF=',S%NLDIFF, S%NLMAX, S%NLMIN
         IREFINE=1
         DO ILEVEL=S%NLMAX-1,S%NLMIN,-1

            IREFINE=IREFINE*2

         IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,*) 'ILEVEL=',ILEVEL
         IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,*) 'IREFINE=',IREFINE

            OS%SLEVEL(ILEVEL)%NW_FACE=OS%SLEVEL(ILEVEL+1)%NW_FACE/IREFINE
            ALLOCATE (OS%SLEVEL(ILEVEL)%IJKW_FACE(17,OSLMAX%NW_FACE), STAT=IERR)
            CALL CHKMEMERR ('SCARC_INIT_COMMUNICATION2D', 'OS%SLEVEL%IJKW_FACE', IERR)
            OS%SLEVEL(ILEVEL)%IJKW_FACE = 0

            IF (SCARC_DEBUG>=2) THEN
               WRITE(SCARC_LU,*) 'OSL zeigt auf OS%SLEVEL(',ILEVEL,')'
               WRITE(SCARC_LU,*) 'OS%SLEVEL(',ILEVEL,')%NW_FACE=',OS%SLEVEL(ILEVEL)%NW_FACE
               WRITE(SCARC_LU,*) 'SIZE(OS%SLEVEL(ILEVEL)%IJKW_FACE)=',SIZE(OS%SLEVEL(ILEVEL)%IJKW_FACE)
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
      IF (SL%NIC_FACE(NM, NOM) == 0 .AND. SL%NIC_FACE(NOM, NM) == 0) CYCLE FACE_NBR_LOOP2D

      ! get length of communication vectors for different orientation of neighbors
      IMIN = SL%I_MIN_FACE(NM,NOM)
      IMAX = SL%I_MAX_FACE(NM,NOM)
      !JMIN = SL%J_MIN_FACE(NM,NOM)
      !JMAX = SL%J_MAX_FACE(NM,NOM)
      JMIN = 1
      JMAX = 1
      KMIN = SL%K_MIN_FACE(NM,NOM)
      KMAX = SL%K_MAX_FACE(NM,NOM)

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
         WRITE(SCARC_LU,*) 'SL%I_MIN_FACE=', SL%I_MIN_FACE
         WRITE(SCARC_LU,*) 'SL%I_MAX_FACE=', SL%I_MAX_FACE
         WRITE(SCARC_LU,*) 'SL%K_MIN_FACE=', SL%K_MIN_FACE
         WRITE(SCARC_LU,*) 'SL%K_MAX_FACE=', SL%K_MAX_FACE
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
INTEGER :: I, K, IC, IW, IERR, ICELL, IOR, ILEVEL
REAL(EB):: DBC
REAL(EB):: TNOW_MATRICES2D
TYPE (MESH_TYPE), POINTER :: M
 
TNOW_MATRICES2D = SECOND()

IERR=0
 
M => MESHES(NM)
 
!!!
!!! Initialize matrices on each level (only 1 level for CG)
!!!
MATRIX_LEVEL_LOOP2D: DO ILEVEL=S%NLMAX,S%NLMIN,-1

   SL => S%SLEVEL(ILEVEL)

    
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Initialize matrix diagonals for global matrix
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   NDIAG = 5
    
   !!! Allocate full matrix corresponding to the band-wise storage technique
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
    
   !!! local matrix AL is momentarily set to global matrix AG
   !SL%AL = SL%AG
    
   IF (ILEVEL/=S%NLMIN) CALL SCARC_PRECON_INIT_GSTRIX2D(ILEVEL,NM)
   IF (ILEVEL==S%NLMAX) THEN

   DO ICELL = 1,SL%NW_FACE

      IOR = SL%IJKW_FACE( 4,ICELL)
      I   = SL%IJKW_FACE( 6,ICELL)
      K   = SL%IJKW_FACE( 8,ICELL)
      NOM = SL%IJKW_FACE( 9,ICELL)
      IW  = SL%IJKW_FACE(16,ICELL)

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

      IF (M%BOUNDARY_TYPE(IW)==OPEN_BOUNDARY) THEN
          M%PRESSURE_BC_INDEX(IW)=DIRICHLET
      ELSE IF (M%IJKW(9,IW)/=0) THEN
          M%PRESSURE_BC_INDEX(IW)=INTERNAL
      ELSE IF (M%BOUNDARY_TYPE(IW)==NULL_BOUNDARY) THEN
          M%PRESSURE_BC_INDEX(IW)=DIRICHLET
      ELSE
          M%PRESSURE_BC_INDEX(IW)=NEUMANN
      ENDIF

      
      SELECT CASE(M%PRESSURE_BC_INDEX(IW))
         CASE (DIRICHLET)                    ! set Dirichlet BC's along open boundary cells
            SL%AG(IC, 1) = SL%AG(IC, 1) - DBC
            IF (SCARC_DEBUG .GE. 4) WRITE(SCARC_LU,1001) IW, IOR, IC, SL%AG(IC,1),M%BOUNDARY_TYPE(IW),NOM
         CASE (INTERNAL)                     ! do nothing along internal boundaries (just only debug message)
            IF (SCARC_DEBUG .GE. 4) WRITE(SCARC_LU,1003) IW, IOR, IC, SL%AG(IC,1),M%BOUNDARY_TYPE(IW),NOM
         CASE (NEUMANN)                      ! set Neumann BC's at all other nodes
            SL%AG(IC, 1) = SL%AG(IC, 1) + DBC
            IF (SCARC_DEBUG .GE. 4) WRITE(SCARC_LU,1002) IW, IOR, IC, SL%AG(IC,1),M%BOUNDARY_TYPE(IW),NOM
      END SELECT

   ENDDO
   ENDIF

   SL%ASUBX = SL%DXI2
   SL%ASUBZ = SL%DZI2
 
   IF (SCARC_DEBUG >= 2) THEN
      !WRITE(SCARC_LU,'(16i4)') (M%BOUNDARY_TYPE(ICELL),ICELL=1,M%NEWC)
      !WRITE(SCARC_LU,*) 'IJKW'
      !DO IW=1,M%NEWC+4
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
INTEGER :: I, J, K, IC, IW, IERR, ICELL, IOR, ILEVEL
REAL(EB):: DBC
REAL(EB):: TNOW_MATRICES3D
TYPE (MESH_TYPE), POINTER :: M
 
TNOW_MATRICES3D = SECOND()

IERR=0

M => MESHES(NM)
 
!!!
!!! Initialize matrices on each level (only 1 level for CG)
!!!
MATRIX_LEVEL_LOOP3D: DO ILEVEL=S%NLMAX,S%NLMIN,-1

   SL => S%SLEVEL(ILEVEL)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Initialize matrix diagonals for global matrix
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   NDIAG = 7
    
   !!! Allocate full matrix corresponding to the band-wise storage technique
   ALLOCATE (SL%AG(1:SL%NCELLS_LOCAL, 1:NDIAG), STAT=IERR)
   CALL CHKMEMERR ('SCARC_INIT_MATRICES3D', 'SL%AG', IERR)
   SL%AG = 0.0_EB
    
   SL%DXI2=1.0_EB/(SL%DX)**2
   SL%DYI2=1.0_EB/(SL%DY)**2
   SL%DZI2=1.0_EB/(SL%DZ)**2

   IF (SCARC_DEBUG >= 2) THEN
      WRITE(SCARC_LU,*) '========= LEVEL=',ILEVEL
      WRITE(SCARC_LU,*) ' IBAR=',SL%IBAR
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
    
   !!! local matrix AL is momentarily set to global matrix AG - not yet needed at all
   !SL%AL = SL%AG
    
   IF (ILEVEL/=S%NLMIN) CALL SCARC_PRECON_INIT_GSTRIX3D(ILEVEL,NM)
   IF (ILEVEL/=S%NLMAX) CYCLE

   IF (SCARC_DEBUG>=2) THEN
      WRITE(SCARC_LU,*) 'BOUNDARY_TYPE'
      WRITE(SCARC_LU,'(16i4)') (M%BOUNDARY_TYPE(ICELL),ICELL=1,M%NEWC)
      WRITE(SCARC_LU,*) 'IJKW'
      DO IW=1,M%NEWC+4
         WRITE(SCARC_LU,'(i4,a,15i4)') IW,' : ',(M%IJKW(ICELL,IW),ICELL=1,15)
      ENDDO
   ENDIF

   DO ICELL = 1,SL%NW_FACE

      IOR = SL%IJKW_FACE( 4,ICELL)
      I   = SL%IJKW_FACE( 6,ICELL)
      J   = SL%IJKW_FACE( 7,ICELL)
      K   = SL%IJKW_FACE( 8,ICELL)
      NOM = SL%IJKW_FACE( 9,ICELL)
      IW  = SL%IJKW_FACE(16,ICELL)

      SELECT CASE(IOR)
         CASE (1)
            IC = (K-1) * SL%IBAR * SL%JBAR + (J-1) * SL%IBAR + I
            DBC= SL%DXI2
         CASE (-1)
            IC = K * SL%IBAR * SL%JBAR + (J-1) * SL%IBAR + I
            DBC= SL%DXI2
         CASE (2)
            IC = (K-1) * SL%IBAR * SL%JBAR + I
            DBC= SL%DYI2
         CASE (-2)
            IC = K * SL%IBAR * SL%JBAR + (SL%JBAR-1) * SL%IBAR + I
            DBC= SL%DYI2
         CASE (3)
            IC = (J-1) * SL%IBAR + I
            DBC= SL%DZI2
         CASE (-3)
            IC = (SL%KBAR-1) * SL%IBAR * SL%JBAR + (J-1) * SL%IBAR + I
            DBC= SL%DZI2
      END SELECT

      IF (M%BOUNDARY_TYPE(IW)==OPEN_BOUNDARY) THEN
          M%PRESSURE_BC_INDEX(IW)=DIRICHLET
      ELSE IF (M%IJKW(9,IW)/=0) THEN
          M%PRESSURE_BC_INDEX(IW)=INTERNAL
      ELSE IF (M%BOUNDARY_TYPE(IW)==NULL_BOUNDARY) THEN
          M%PRESSURE_BC_INDEX(IW)=DIRICHLET
      ELSE
          M%PRESSURE_BC_INDEX(IW)=NEUMANN
      ENDIF

      SELECT CASE(M%PRESSURE_BC_INDEX(IW))
         CASE (DIRICHLET)                    ! set Dirichlet BC's for AG and AL in open boundary cells
            SL%AG(IC, 1) = SL%AG(IC, 1) - DBC
            IF (SCARC_DEBUG .GE. 4) WRITE(SCARC_LU,1001) IW, IOR, IC, SL%AG(IC,1),M%BOUNDARY_TYPE(IW),NOM
         CASE (INTERNAL)                     ! do nothing along internal boundaries (only debug message)
            IF (SCARC_DEBUG .GE. 4) WRITE(SCARC_LU,1003) IW, IOR, IC, SL%AG(IC,1),M%BOUNDARY_TYPE(IW),NOM
         CASE (NEUMANN)                      ! set Neumann BC's for AG and AL at all other nodes
            SL%AG(IC, 1) = SL%AG(IC, 1) + DBC
            IF (SCARC_DEBUG .GE. 4) WRITE(SCARC_LU,1002) IW, IOR, IC, SL%AG(IC,1),M%BOUNDARY_TYPE(IW),NOM
      END SELECT

   ENDDO

   IF (SCARC_DEBUG >= 2) THEN
      WRITE(SCARC_LU,*) 'ASUB=', SL%ASUB
      WRITE(SCARC_LU,*)
      WRITE(SCARC_LU,*) 'GLOBAL MATRIX AG:'
      DO IC = 1, SL%NCELLS_LOCAL
         WRITE(SCARC_LU, '(5f12.4)') SL%AG(IC, 2),SL%AG(IC,3),SL%AG(IC,1),SL%AG(IC,4),SL%AG(IC,5)
         IF (Mod(IC, SL%IBAR) == 0) WRITE(SCARC_LU,*) '---------------------------------------------&
                                                 &--------------------'
      ENDDO
   ENDIF
    
   SL%ASUBX = SL%DXI2
   SL%ASUBY = SL%DYI2
   SL%ASUBZ = SL%DZI2
 
   IF (SCARC_DEBUG >= 2) THEN
     WRITE(SCARC_LU,*) 'SL%ASUBX=', SL%ASUBX
     WRITE(SCARC_LU,*) 'LBC=', LBC
     WRITE(SCARC_LU,*) 'MBC=', MBC
     WRITE(SCARC_LU,*) 'NBC=', NBC
   ENDIF

ENDDO MATRIX_LEVEL_LOOP3D
 
TUSED_SCARC(4,NM)=TUSED_SCARC(4,NM)+SECOND()-TNOW_MATRICES3D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_MATRICES3D

1001 FORMAT('IW=',i3,': Dirichlet , IOR=',i2,':AG(',i3,'1)=',f12.6,': BT=',i3,': NOM=',i3)
1002 FORMAT('IW=',i3,': Neumann   , IOR=',i2,':AG(',i3,'1)=',f12.6,': BT=',i3,': NOM=',i3)
1003 FORMAT('IW=',i3,': Nothing   , IOR=',i2,':AG(',i3,'1)=',f12.6,': BT=',i3,': NOM=',i3)
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
   ! working and auxiliary vectors for global CG/MG-method
   ! ----------------------------------------------------------------------------
   ALLOCATE (SL%X(0:IBP0, 1, 0:KBP0), STAT=IERR)
   CALL CHKMEMERR ('SCARC', 'X', IERR)
   SL%X = 0.0_EB
    
   ALLOCATE (SL%F(0:IBP0, 1, 0:KBP0), STAT=IERR)
   CALL CHKMEMERR ('SCARC', 'B', IERR)
   SL%F = 0.0_EB
    
   ALLOCATE (SL%D(0:IBP0, 1, 0:KBP0), STAT=IERR)
   CALL CHKMEMERR ('SCARC', 'D', IERR)
   SL%D = 0.0_EB
    
   IF ((SCARC_METHOD == 'CG' .OR. SCARC_METHOD == 'BICG').OR.ILEVEL==S%NLMIN) THEN
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
    
   ALLOCATE (SL%Z(0:IBP0, 1, 0:KBP0), STAT=IERR)
   CALL CHKMEMERR ('SCARC', 'Z', IERR)
   SL%Z = 0.0_EB
    
   IF (SCARC_CG_PRECON == 'MG') THEN

      ALLOCATE (SL%X2(0:IBP0, 1, 0:KBP0), STAT=IERR)
      CALL CHKMEMERR ('SCARC', 'X2', IERR)
      SL%X2 = 0.0_EB

      ALLOCATE (SL%F2(0:IBP0, 1, 0:KBP0), STAT=IERR)
      CALL CHKMEMERR ('SCARC', 'F2', IERR)
      SL%F2 = 0.0_EB

      ALLOCATE (SL%D2(0:IBP0, 1, 0:KBP0), STAT=IERR)
      CALL CHKMEMERR ('SCARC', 'D2', IERR)
      SL%D2 = 0.0_EB
    
   ENDIF

   ! ----------------------------------------------------------------------------
   ! vectors for the description of the boundary conditions
   ! ----------------------------------------------------------------------------
   ALLOCATE (SL%BXS0(1, KBP0), STAT=IERR)
   CALL CHKMEMERR ('SCARC', 'BXS0', IERR)
   SL%BXS0 = 0.0_EB
  
   ALLOCATE (SL%BXF0(1, KBP0), STAT=IERR)
   CALL CHKMEMERR ('SCARC', 'BXF0', IERR)
   SL%BXF0 = 0.0_EB

   ALLOCATE (SL%BZS0(IBP0, 1), STAT=IERR)
   CALL CHKMEMERR ('SCARC', 'BZS0', IERR)
   SL%BZS0 = 0.0_EB

   ALLOCATE (SL%BZF0(IBP0, 1), STAT=IERR)
   CALL CHKMEMERR ('SCARC', 'BZF0', IERR)
   SL%BZF0 = 0.0_EB

ENDDO LEVEL_LOOP2D

    
!!!
!!! Allocate array for the description of the MG-cycle
!!!
IF (SCARC_METHOD == 'MG' .OR. SCARC_CG_PRECON =='MG') THEN
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
   ! working and auxiliary vectors for MG-method
   ! ----------------------------------------------------------------------------
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
    
   ALLOCATE (SL%Z(0:IBP0, 0:JBP1, 0:KBP0), STAT=IERR)
   CALL CHKMEMERR ('SCARC', 'Z', IERR)
   SL%Z = 0.0_EB
    

   ! ----------------------------------------------------------------------------
   ! vectors for the description of the boundary conditions
   ! ----------------------------------------------------------------------------
   ALLOCATE (SL%BXS0(JBP0, KBP0), STAT=IERR)
   CALL CHKMEMERR ('SCARC', 'BXS0', IERR)
   SL%BXS0 = 0.0_EB
    
   ALLOCATE (SL%BXF0(JBP0, KBP0), STAT=IERR)
   CALL CHKMEMERR ('SCARC', 'BXF0', IERR)
   SL%BXF0 = 0.0_EB
    
   ALLOCATE (SL%BYS0(IBP0, KBP0), STAT=IERR)
   CALL CHKMEMERR ('SCARC', 'BYS0', IERR)
   SL%BYS0 = 0.0_EB
    
   ALLOCATE (SL%BYF0(IBP0, KBP0), STAT=IERR)
   CALL CHKMEMERR ('SCARC', 'BYF0', IERR)
   SL%BYF0 = 0.0_EB
    
   ALLOCATE (SL%BZS0(IBP0, JBP0), STAT=IERR)
   CALL CHKMEMERR ('SCARC', 'BZS0', IERR)
   SL%BZS0 = 0.0_EB
    
   ALLOCATE (SL%BZF0(IBP0, JBP0), STAT=IERR)
   CALL CHKMEMERR ('SCARC', 'BZF0', IERR)
   SL%BZF0 = 0.0_EB
   
ENDDO LEVEL_LOOP3D

    
!!!
!!! Allocate array for the description of the MG-cycle
!!!
IF (SCARC_METHOD == 'MG' .OR. SCARC_CG_PRECON =='MG') THEN
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

! set boundary conditions along exterior boundaries
CALL SCARC_SETBDRY2D (ILMAX, NM)
 
!!! calculate initial defect Y = B - A*X and get l2-norm of it
NREQ_FACE=0
CALL SCARC_MATVEC2D (SL%X, SL%R, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILMAX, NMV_X, NMV_R)
SL%R = -SL%F + SL%R
SL%RESIN_CG = SCARC_L2NORM2D (SL%R, ILMAX, NM, NTYPE_GLOBAL)

IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) 0, ILMAX, SL%RESIN_CG
IF (SCARC_DEBUG >= 1) write (*,1000) 0, ILMAX, SL%RESIN_CG
 
CALL SCARC_SHOW_LEVEL (SL%X, 'SARC', 'XCG0  ',ILMAX)
CALL SCARC_SHOW_LEVEL (SL%F, 'SARC', 'FCG0  ',ILMAX)
CALL SCARC_SHOW_LEVEL (SL%R, 'SARC', 'RCG0  ',ILMAX)

CALL SCARC_SHOW_LEVEL_LINE(SL%F, 'SARC', 'FLINE0',ILMAX)


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
END SELECT
SL%D = -SL%G
 
CALL SCARC_SHOW_LEVEL_LINE(SL%D, 'SARC', 'DLINE0',ILMAX)
CALL SCARC_SHOW_LEVEL (SL%D, 'SARC', 'DLINEA',ILMAX)
!
! defect correction loop
!
CG_LOOP2D: DO ITE = 1, SCARC_CG_NIT
 
   ! calculate new defect and get L2-norm of it
NREQ_FACE=0
   CALL SCARC_MATVEC2D (SL%D, SL%Y, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILMAX, NMV_D, NMV_Y)
   ALPHA = SCARC_SCALPROD2D (SL%D, SL%Y, ILMAX, NM)
CALL SCARC_SHOW_LEVEL_LINE(SL%Y, 'SARC', 'YLINE0',ILMAX)

   ! get new descent direction
   ALPHA = SIGMA0 / ALPHA
   SL%X = ALPHA * SL%D + SL%X
   SL%R = ALPHA * SL%Y + SL%R
   SL%RES_CG = SCARC_L2NORM2D (SL%R, ILMAX, NM, NTYPE_GLOBAL)
 
CALL SCARC_SHOW_LEVEL_LINE(SL%X, 'SARC', 'XLINE1',ILMAX)
CALL SCARC_SHOW_LEVEL_LINE(SL%R, 'SARC', 'RLINE1',ILMAX)
   !  exit loop in case of convergence or divergence
   IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) ITE, ILMAX, SL%RES_CG
   IF (SCARC_DEBUG >= 1) write (*,1000) ITE, ILMAX, SL%RES_CG
   IF (NM == 1)          write (*,1000) ITE, ILMAX, SL%RES_CG
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
         CALL SCARC_MG_PRECON2D (SL%G, ILMAX, NM)
         SIGMA1 = SCARC_SCALPROD2D (SL%R, SL%G, ILMAX, NM)
   END SELECT
 
   GAMMA = SIGMA1 / SIGMA0
   SIGMA0 = SIGMA1
   SL%D = -SL%G + GAMMA * SL%D
 
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
IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,2000) ITE0, ILMAX, SL%RES_CG, SL%CAPPA_CG 
IF (NM == 1)          WRITE (*,2000) ITE0, ILMAX, SL%RES_CG, SL%CAPPA_CG

! shift solution back to vector HP
DO K = 1, SL%KBAR                                 
   DO I = 1, SL%IBAR
      HP(I, 1, K) = SL%X (I, 1, K)
   ENDDO
ENDDO

! set ghost cell values along exterior boundaries 
CALL SCARC_SHOW_LINE(HP, 'SARC', 'HLINE0')
CALL SCARC_GHOSTCELLS(HP,NM)

CALL SCARC_SHOW_LINE(HP, 'SARC', 'HLINE1')
! get neighbouring data along interior boundaries
CALL SCARC_UPDATE_LEVEL(HP, NM, ILMAX, 'HPFIN ')
CALL SCARC_SHOW_LINE(HP, 'SARC', 'HLINE2')


TUSED_SCARC(6,NM)=TUSED_SCARC(6,NM)+SECOND()-TNOW_CG2D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_CG2D

1000 FORMAT ('=====SCARC_CG       : #ite= ',i4,': level=',i4,': res=',e12.4)
2000 FORMAT ('=====SCARC_CG       : #ite= ',i4,': level=',i4,': res=',e12.4,': rate=',e12.4)
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

IF (SCARC_DEBUG >=2) WRITE(SCARC_LU,*) 'STARTING BICG2D'

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
 
CALL SCARC_SHOW_LEVEL (SL%X, 'SARC', 'X0    ',ILMAX)
CALL SCARC_SHOW_LEVEL (SL%F, 'SARC', 'F0    ',ILMAX)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! initialization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CALL SCARC_SETBDRY2D (ILMAX, NM)
SL%R=SL%F

SELECT CASE(SCARC_BICG_PRECON)
   CASE('JACOBI') 
      CALL SCARC_PRECON_JACOBI2D (SL%R, ILMAX, NM)
CALL SCARC_SHOW_LEVEL (SL%R, 'SARC', 'R01   ',ILMAX)
   CASE('SSOR') 
      CALL SCARC_PRECON_SSOR2D (SL%R, ILMAX, NM)
CALL SCARC_SHOW_LEVEL (SL%R, 'SARC', 'R02   ',ILMAX)
   CASE('GSTRIX') 
      CALL SCARC_PRECON_GSTRIX2D (SL%R, ILMAX, NM)
CALL SCARC_SHOW_LEVEL (SL%R, 'SARC', 'R02   ',ILMAX)
   CASE('MG') 
      CALL SCARC_MG_PRECON2D (SL%R, ILMAX, NM)
CALL SCARC_SHOW_LEVEL (SL%R, 'SARC', 'R04   ',ILMAX)
END SELECT

CALL SCARC_MATVEC2D (SL%X, SL%R, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILMAX, NMV_X, NMV_R)
CALL SCARC_SHOW_LEVEL (SL%X, 'SARC', 'X1    ',ILMAX)
CALL SCARC_SHOW_LEVEL (SL%R, 'SARC', 'R1    ',ILMAX)
SL%R = SL%F - SL%R
CALL SCARC_SHOW_LEVEL (SL%R, 'SARC', 'R2    ',ILMAX)


SELECT CASE(SCARC_BICG_PRECON)
   CASE('JACOBI') 
      CALL SCARC_PRECON_JACOBI2D (SL%R, ILMAX, NM)
   CASE('SSOR') 
      CALL SCARC_PRECON_SSOR2D (SL%R, ILMAX, NM)
   CASE('GSTRIX') 
      CALL SCARC_PRECON_GSTRIX2D (SL%R, ILMAX, NM)
   CASE('MG') 
      CALL SCARC_MG_PRECON2D (SL%R, ILMAX, NM)
END SELECT
CALL SCARC_SHOW_LEVEL (SL%R, 'SARC', 'R3    ',ILMAX)

SL%RESIN_BICG = SCARC_L2NORM2D (SL%R, ILMAX, NM, NTYPE_GLOBAL)
SL%G=SL%R

IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) 0, ILMAX, SL%RESIN_BICG
IF (SCARC_DEBUG >= 1) write (*,1000) 0, ILMAX, SL%RESIN_BICG


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! iterative correction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
BICG_LOOP2D: DO ITE = 1, SCARC_BICG_NIT
 
   RHO1 = SCARC_SCALPROD2D (SL%G, SL%R, ILMAX, NM)
   DBETA=(RHO1*DTHETA)/(RHO0*ALPHA0)

   RHO0=RHO1

   SL%Z = SL%R + DBETA * SL%Z
   SL%Z = -DBETA * ALPHA0 * SL%Y + SL%Z

CALL SCARC_SHOW_LEVEL (SL%Z, 'SARC', 'Z4    ',ILMAX)
CALL SCARC_SHOW_LEVEL (SL%Y, 'SARC', 'Y4    ',ILMAX)

   ! matrix-vector multiplication and preconditioning
   CALL SCARC_MATVEC2D (SL%Z, SL%Y, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILMAX, NMV_Z, NMV_Y)

CALL SCARC_SHOW_LEVEL (SL%Z, 'SARC', 'Z5    ',ILMAX)
CALL SCARC_SHOW_LEVEL (SL%Y, 'SARC', 'Y5    ',ILMAX)

   SELECT CASE(SCARC_BICG_PRECON)
      CASE('JACOBI') 
         CALL SCARC_PRECON_JACOBI2D (SL%Y, ILMAX, NM)
      CASE('SSOR') 
         CALL SCARC_PRECON_SSOR2D (SL%Y, ILMAX, NM)
      CASE('GSTRIX') 
         CALL SCARC_PRECON_GSTRIX2D (SL%Y, ILMAX, NM)
      CASE('MG') 
         CALL SCARC_MG_PRECON2D (SL%Y, ILMAX, NM)
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
         CALL SCARC_MG_PRECON2D (SL%D, ILMAX, NM)
   END SELECT

   ALPHA1 = SCARC_SCALPROD2D (SL%D, SL%R, ILMAX, NM)
   ALPHA2 = SCARC_SCALPROD2D (SL%D, SL%D, ILMAX, NM)

   ALPHA0=ALPHA1/ALPHA2

   SL%X =  DTHETA * SL%Z + SL%X
   SL%X =  ALPHA0 * SL%R + SL%X
   SL%R = -ALPHA0 * SL%D + SL%R

   SL%RES_BICG = SCARC_L2NORM2D (SL%R, ILMAX, NM, NTYPE_GLOBAL)
 
   !  exit loop in case of convergence or divergence
   IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) ITE, ILMAX, SL%RES_BICG
   IF (SCARC_DEBUG >= 1) write (*,1000) ITE, ILMAX, SL%RES_BICG
   !IF (NM == 1)          write (*,1000) ITE, ILMAX, SL%RES_BICG
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
IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,2000) ITE0, ILMAX, SL%RES_BICG, SL%CAPPA_BICG 
IF (NM == 1)          WRITE (*,2000) ITE0, ILMAX, SL%RES_BICG, SL%CAPPA_BICG

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

1000 FORMAT ('=====SCARC_BICG     : #ite= ',i4,': level=',i4,': res=',e12.4)
2000 FORMAT ('=====SCARC_BICG     : #ite= ',i4,': level=',i4,': res=',e12.4,': rate=',e12.4)
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

      CALL SCARC_COARSE2D(1, ILMAX, NM)

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
   IF (SCARC_DEBUG >= 1) write (*,1000) 0, ILMAX, SL%RESIN_MG

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
         CALL SCARC_SHOW_LEVEL (SLHI%D, 'SARC',  'RESTD ',ILEVEL)
         CALL SCARC_RESTRICTION2D(SLLO%F, SLHI%D, ILEVEL-1, NM)
         CALL SCARC_SHOW_LEVEL (SLLO%F, 'SARC',  'RESTF ',ILEVEL-1)

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
      CALL SCARC_COARSE2D(1, ILEVEL, NM)


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
      IF (SCARC_DEBUG >= 1) write (*,1000) ITE, ILEVEL, SL%RES_MG
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
   IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,2000) ITE0, ILEVEL, SL%RES_MG, SL%CAPPA_MG 
   IF (SCARC_DEBUG >= 0) WRITE (*,2000) ITE0, ILEVEL, SL%RES_MG, SL%CAPPA_MG 

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

1000 FORMAT ('=====SCARC_MG       : #ite= ',i4,': level=',i4,': res=',e12.4)
2000 FORMAT ('=====SCARC_MG       : #ite= ',i4,': level=',i4,': res=',e12.4,': rate=',e12.4)
END SUBROUTINE SCARC_MG2D



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Perform global 2D cg-method based on global possion-matrix
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_COARSE2D (ILAYER, ILEVEL, NM)
 
INTEGER ::  NM, ILAYER, ILEVEL, ITE, ITE0
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
IF (ILAYER==2) SL%F = SL%F2
SL%X = 0.0_EB                   
SL%G = 0.0_EB                   
SL%Y = 0.0_EB
SL%R = 0.0_EB
SL%D = 0.0_EB
 
IF (SCARC_DEBUG>2) WRITE(SCARC_LU,*) 'ILEVEL=',ILEVEL
CALL SCARC_SHOW_LEVEL (SL%X, 'SARC', 'XCO0  ',ILEVEL)
CALL SCARC_SHOW_LEVEL (SL%F, 'SARC', 'FCO0  ',ILEVEL)

!!! calculate initial defect Y = B - A*X and get l2-norm of it
CALL SCARC_MATVEC2D (SL%X, SL%R, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, NMV_X, NMV_R)
CALL SCARC_SHOW_LEVEL (SL%R, 'SARC', 'RCO00 ',ILEVEL)
SL%R = -SL%F + SL%R
CALL SCARC_SHOW_LEVEL (SL%R, 'SARC', 'RCO01 ',ILEVEL)
SL%RESIN_CO = SCARC_L2NORM2D (SL%R, ILEVEL, NM, NTYPE_GLOBAL)
IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) 0, ILEVEL, SL%RESIN_CO
IF (SCARC_DEBUG >= 1) write (*,1000) 0, ILEVEL, SL%RESIN_CO
 
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
 
CALL SCARC_SHOW_LEVEL (SL%R, 'SARC', 'RCO02 ',ILEVEL)
CALL SCARC_SHOW_LEVEL (SL%G, 'SARC', 'GCO02 ',ILEVEL)
CALL SCARC_SHOW_LEVEL (SL%D, 'SARC', 'DCO02 ',ILEVEL)
!
! start defect correction loop
!
CO_LOOP2D: DO ITE = 1, SCARC_CO_NIT
 
   ! calculate new defect and get L2-norm of it
   CALL SCARC_MATVEC2D (SL%D, SL%Y, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, NMV_D, NMV_Y)
   ALPHA = SCARC_SCALPROD2D (SL%D, SL%Y, ILEVEL, NM)

CALL SCARC_SHOW_LEVEL (SL%Y, 'SARC', 'YCO1  ',ILEVEL)

   ! get new descent direction
   ALPHA = SIGMA0 / ALPHA
   SL%X = ALPHA * SL%D + SL%X
   SL%R = ALPHA * SL%Y + SL%R
   SL%RES_CO = SCARC_L2NORM2D (SL%R, ILEVEL, NM, NTYPE_GLOBAL)

 
   !  exit loop in case of convergence or divergence
   IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) ITE, ILEVEL, SL%RES_CO
   IF (SCARC_DEBUG >= 1) write (*,1000) ITE, ILEVEL, SL%RES_CO
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
 
CALL SCARC_SHOW_LEVEL (SL%R, 'SARC', 'RCO1  ',ILEVEL)
CALL SCARC_SHOW_LEVEL (SL%G, 'SARC', 'GCO1  ',ILEVEL)

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
IF (SCARC_DEBUG >= 1) WRITE (*,2000) ITE0, ILEVEL, SL%RES_CO, SL%CAPPA_CO 

 
TUSED_SCARC(9,NM)=TUSED_SCARC(9,NM)+SECOND()-TNOW_COARSE
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_COARSE

1000 FORMAT ('     SCARC_COARSE   : #ite= ',i4,': level=',i4,': res=',e12.4)
2000 FORMAT ('     SCARC_COARSE   : #ite= ',i4,': level=',i4,': res=',e12.4,': rate=',e12.4)
END SUBROUTINE SCARC_COARSE2D


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

CALL SCARC_SHOW_LEVEL (X, 'SARC',  'SMOX  ',ILEVEL)
CALL SCARC_SHOW_LEVEL (F, 'SARC',  'SMOF  ',ILEVEL)
IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,*) 'BMATVEC=',BMATVEC

SL%RESIN_SM=1.0_EB
IF (BMATVEC) THEN
   CALL SCARC_MATVEC2D (X, D, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, NMV_X, NMV_D)   ! global MATVEC
   D = F - D
   CALL SCARC_SHOW_LEVEL (D, 'SARC',  'SMOD  ',ILEVEL)
   IF (BL2NORM) SL%RESIN_SM = SCARC_L2NORM2D (D, ILEVEL, NM, NTYPE_GLOBAL)
ENDIF

IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) 0, ILEVEL, SL%RESIN_SM
IF (SCARC_DEBUG >= 1) WRITE(*,1000) 0, ILEVEL, SL%RESIN_SM

SMOOTH_LOOP2D: DO ITE=1,SCARC_SM_NIT
 
   CALL SCARC_SHOW_LEVEL (D, 'SARC',  'SMOD  ',ILEVEL)
   SELECT CASE(SCARC_SM_PRECON)
      CASE('JACOBI') 
         CALL SCARC_PRECON_JACOBI2D (D, ILEVEL, NM)
   CALL SCARC_SHOW_LEVEL (D, 'SARC',  'SMODJ ',ILEVEL)
      CASE('SSOR') 
         CALL SCARC_PRECON_SSOR2D (D, ILEVEL, NM)
   CALL SCARC_SHOW_LEVEL (D, 'SARC',  'SMODS ',ILEVEL)
      CASE('GSTRIX') 
         CALL SCARC_PRECON_GSTRIX2D (D, ILEVEL, NM)
   CALL SCARC_SHOW_LEVEL (D, 'SARC',  'SMODG ',ILEVEL)
   END SELECT

   X = SCARC_SM_OMEGA * D + X
   CALL SCARC_SHOW_LEVEL (X, 'SARC',  'SMOX2 ',ILEVEL)
   CALL SCARC_SHOW_LEVEL (D, 'SARC',  'SMOD2 ',ILEVEL)
   CALL SCARC_MATVEC2D (X, D, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, NMV_X, NMV_D)   ! global MATVEC
   CALL SCARC_SHOW_LEVEL (X, 'SARC',  'SMOX3 ',ILEVEL)
   CALL SCARC_SHOW_LEVEL (D, 'SARC',  'SMOD3 ',ILEVEL)
   D = F - D
   CALL SCARC_SHOW_LEVEL (D, 'SARC',  'SMOD4 ',ILEVEL)
   IF (BL2NORM) SL%RES_SM = SCARC_L2NORM2D (D, ILEVEL, NM, NTYPE_GLOBAL)

   IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) ITE, ILEVEL, SL%RES_SM
   IF (SCARC_DEBUG >= 1) WRITE(*,1000) ITE, ILEVEL, SL%RES_SM
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
IF (SCARC_DEBUG >= 1) WRITE (*,2000) ITE0, ILEVEL, SL%RES_SM, SL%CAPPA_SM 

TUSED_SCARC(8,NM)=TUSED_SCARC(8,NM)+SECOND()-TNOW_SMOOTHER2D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_SMOOTHER2D

1000 FORMAT ('     SCARC_SMOOTHER : #ite= ',i4,': level=',i4,': res=',e12.4)
2000 FORMAT ('     SCARC_SMOOTHER : #ite= ',i4,': level=',i4,': res=',e12.4,': rate=',e12.4)
END SUBROUTINE SCARC_SMOOTHER2D

 
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

write(*,*) 'Still experimental'
stop

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
 
! set boundary conditions along exterior boundaries
CALL SCARC_SETBDRY3D (ILMAX, NM)
 
!!! calculate initial defect Y = B - A*X and get l2-norm of it
CALL SCARC_MATVEC3D (SL%X, SL%R, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILMAX, NMV_X, NMV_R)
SL%R = -SL%F + SL%R
SL%RESIN_CG = SCARC_L2NORM3D (SL%R, ILMAX, NM, NTYPE_GLOBAL)

IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) 0, ILMAX, SL%RESIN_CG
IF (SCARC_DEBUG >= 1) write (*,1000) 0, ILMAX, SL%RESIN_CG
 
!CALL SCARC_SHOW_LEVEL (SL%X, 'SARC', 'X0    ',ILMAX)
!CALL SCARC_SHOW_LEVEL (SL%F, 'SARC', 'F0    ',ILMAX)
!CALL SCARC_SHOW_LEVEL (SL%R, 'SARC', 'R0    ',ILMAX)

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
END SELECT
SL%D = -SL%G
 
!
! defect correction loop
!
CG_LOOP3D: DO ITE = 1, SCARC_CG_NIT
 
   ! calculate new defect and get L2-norm of it
   CALL SCARC_MATVEC3D (SL%D, SL%Y, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILMAX, NMV_D, NMV_Y)
   ALPHA = SCARC_SCALPROD3D (SL%D, SL%Y, ILMAX, NM)

   ! get new descent direction
   ALPHA = SIGMA0 / ALPHA
   SL%X = ALPHA * SL%D + SL%X
   SL%R = ALPHA * SL%Y + SL%R
   SL%RES_CG = SCARC_L2NORM3D (SL%R, ILMAX, NM, NTYPE_GLOBAL)
 
   !  exit loop in case of convergence or divergence
   IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) ITE, ILMAX, SL%RES_CG
   IF (SCARC_DEBUG >= 1) write (*,1000) ITE, ILMAX, SL%RES_CG
   IF (NM == 1)          write (*,1000) ITE, ILMAX, SL%RES_CG
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
         CALL SCARC_MG_PRECON3D (SL%G, ILMAX, NM)
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
IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,2000) ITE0, ILMAX, SL%RES_CG, SL%CAPPA_CG 
IF (NM == 1)          WRITE (*,2000) ITE0, ILMAX, SL%RES_CG, SL%CAPPA_CG

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


TUSED_SCARC(6,NM)=TUSED_SCARC(6,NM)+SECOND()-TNOW_CG3D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_CG3D

1000 FORMAT ('=====SCARC_CG       : #ite= ',i4,': level=',i4,': res=',e12.4)
2000 FORMAT ('=====SCARC_CG       : #ite= ',i4,': level=',i4,': res=',e12.4,': rate=',e12.4)
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

write(*,*) 'Still experimental'
stop

S  => SCARC (NM)
M  => MESHES (NM)

IF (SCARC_DEBUG >=2) WRITE(SCARC_LU,*) 'STARTING BICG3D'

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
      CALL SCARC_MG_PRECON3D (SL%R, ILMAX, NM)
END SELECT

CALL SCARC_MATVEC3D (SL%X, SL%R, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILMAX, NMV_X, NMV_R)
SL%R = SL%F - SL%R

CALL SCARC_SHOW(SL%R, 'SARC', 'R0    ')

SELECT CASE(SCARC_BICG_PRECON)
   CASE('JACOBI') 
      CALL SCARC_PRECON_JACOBI3D (SL%R, ILMAX, NM)
   CASE('SSOR') 
      CALL SCARC_PRECON_SSOR3D (SL%R, ILMAX, NM)
   CASE('GSTRIX') 
      CALL SCARC_PRECON_GSTRIX3D (SL%R, ILMAX, NM)
   CASE('MG') 
      CALL SCARC_MG_PRECON3D (SL%R, ILMAX, NM)
END SELECT

SL%RESIN_BICG = SCARC_L2NORM3D (SL%R, ILMAX, NM, NTYPE_GLOBAL)
SL%G=SL%R

IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) 0, ILMAX, SL%RESIN_BICG
IF (SCARC_DEBUG >= 1) write (*,1000) 0, ILMAX, SL%RESIN_BICG

CALL SCARC_SHOW(SL%G, 'SARC', 'G0    ')

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
         CALL SCARC_MG_PRECON3D (SL%Y, ILMAX, NM)
   END SELECT

   DTHETA = SCARC_SCALPROD3D (SL%G, SL%Y, ILMAX, NM)
   DTHETA=RHO1/DTHETA
   SL%R = -DTHETA * SL%Y + SL%R

CALL SCARC_SHOW(SL%R, 'SARC',  'R1    ')
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
         CALL SCARC_MG_PRECON3D (SL%D, ILMAX, NM)
   END SELECT

   ALPHA1 = SCARC_SCALPROD3D (SL%D, SL%R, ILMAX, NM)
   ALPHA2 = SCARC_SCALPROD3D (SL%D, SL%D, ILMAX, NM)

   ALPHA0=ALPHA1/ALPHA2

   SL%X =  DTHETA * SL%Z + SL%X
   SL%X =  ALPHA0 * SL%R + SL%X
   SL%R = -ALPHA0 * SL%D + SL%R

   SL%RES_BICG = SCARC_L2NORM3D (SL%R, ILMAX, NM, NTYPE_GLOBAL)
 
   !  exit loop in case of convergence or divergence
   IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) ITE, ILMAX, SL%RES_BICG
   IF (SCARC_DEBUG >= 1) write (*,1000) ITE, ILMAX, SL%RES_BICG
   !IF (NM == 1)          write (*,1000) ITE, ILMAX, SL%RES_BICG
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
IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,2000) ITE0, ILMAX, SL%RES_BICG, SL%CAPPA_BICG 
IF (NM == 1)          WRITE (*,2000) ITE0, ILMAX, SL%RES_BICG, SL%CAPPA_BICG

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

1000 FORMAT ('=====SCARC_BICG     : #ite= ',i4,': level=',i4,': res=',e12.4)
2000 FORMAT ('=====SCARC_BICG     : #ite= ',i4,': level=',i4,': res=',e12.4,': rate=',e12.4)
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
!TYPE (SCARC_LEVEL_TYPE), POINTER :: SL, SLHI, SLLO
TYPE (MESH_TYPE),  POINTER :: M
 
TNOW_MG3D = SECOND()

write(*,*) 'Still experimental'
stop

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

      CALL SCARC_COARSE3D(1, ILMAX, NM)

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
   IF (SCARC_DEBUG >= 1) write (*,1000) 0, ILMAX, SL%RESIN_MG

   !CALL SCARC_MASTER_INFO (SLHI%ITE, SLHI%RESIN_MG)

    
   !!! ---------------------------------------------------------------------------
   !!! start MG-iteration
   !!! ---------------------------------------------------------------------------
   MG_LOOP3D: DO ITE = 1, SCARC_MG_NIT
    
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
         CALL SCARC_SHOW_LEVEL (SLHI%D, 'SARC',  'RESTD ',ILEVEL)
         CALL SCARC_RESTRICTION3D(SLLO%F, SLHI%D, ILEVEL-1, NM)
         CALL SCARC_SHOW_LEVEL (SLLO%F, 'SARC',  'RESTF ',ILEVEL-1)

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
      CALL SCARC_COARSE3D(1, ILEVEL, NM)


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
      IF (SCARC_DEBUG >= 1) write (*,1000) ITE, ILEVEL, SL%RES_MG
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
   IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,2000) ITE0, ILEVEL, SL%RES_MG, SL%CAPPA_MG 
   IF (SCARC_DEBUG >= 1) WRITE (*,2000) ITE0, ILEVEL, SL%RES_MG, SL%CAPPA_MG 

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

1000 FORMAT ('=====SCARC_MG       : #ite= ',i4,': level=',i4,': res=',e12.4)
2000 FORMAT ('=====SCARC_MG       : #ite= ',i4,': level=',i4,': res=',e12.4,': rate=',e12.4)
END SUBROUTINE SCARC_MG3D



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Perform global 3D cg-method based on global possion-matrix
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_COARSE3D (ILAYER, ILEVEL, NM)
 
INTEGER ::  NM, ILAYER, ILEVEL, ITE, ITE0
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
IF (ILAYER==2) SL%F = SL%F2
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
IF (SCARC_DEBUG >= 1) write (*,1000) 0, ILEVEL, SL%RESIN_CO
 
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
   IF (SCARC_DEBUG >= 1) write (*,1000) ITE, ILEVEL, SL%RES_CO
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
IF (SCARC_DEBUG >= 1) WRITE (*,2000) ITE0, ILEVEL, SL%RES_CO, SL%CAPPA_CO 

 
TUSED_SCARC(9,NM)=TUSED_SCARC(9,NM)+SECOND()-TNOW_COARSE
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_COARSE

1000 FORMAT ('     SCARC_COARSE   : #ite= ',i4,': level=',i4,': res=',e12.4)
2000 FORMAT ('     SCARC_COARSE   : #ite= ',i4,': level=',i4,': res=',e12.4,': rate=',e12.4)
END SUBROUTINE SCARC_COARSE3D


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

CALL SCARC_SHOW_LEVEL (X, 'SARC',  'SMOX  ',ILEVEL)
CALL SCARC_SHOW_LEVEL (F, 'SARC',  'SMOF  ',ILEVEL)

SL%RESIN_SM=1.0_EB
IF (BMATVEC) THEN
   CALL SCARC_MATVEC3D (X, D, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, NMV_X, NMV_D)  
   D = F - D
   CALL SCARC_SHOW_LEVEL (D, 'SARC',  'SMOD  ',ILEVEL)
   IF (BL2NORM) SL%RESIN_SM = SCARC_L2NORM3D (D, ILEVEL, NM, NTYPE_GLOBAL)
ENDIF

IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) 0, ILEVEL, SL%RESIN_SM
IF (SCARC_DEBUG >= 1) WRITE(*,1000) 0, ILEVEL, SL%RESIN_SM

SMOOTH_LOOP3D: DO ITE=1,SCARC_SM_NIT
 
   SELECT CASE(SCARC_SM_PRECON)
      CASE('JACOBI') 
         CALL SCARC_PRECON_JACOBI3D (D, ILEVEL, NM)
      CASE('SSOR') 
         CALL SCARC_PRECON_SSOR3D (D, ILEVEL, NM)
      CASE('GSTRIX') 
         CALL SCARC_PRECON_GSTRIX3D (D, ILEVEL, NM)
   END SELECT

   X = SCARC_SM_OMEGA * D + X
   CALL SCARC_MATVEC3D (X, D, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, NMV_X, NMV_D)   ! global MATVEC
   D = F - D
   IF (BL2NORM) SL%RES_SM = SCARC_L2NORM3D (D, ILEVEL, NM, NTYPE_GLOBAL)

   IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) ITE, ILEVEL, SL%RES_SM
   IF (SCARC_DEBUG >= 1) WRITE(*,1000) ITE, ILEVEL, SL%RES_SM
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
IF (SCARC_DEBUG >= 1) WRITE (*,2000) ITE0, ILEVEL, SL%RES_SM, SL%CAPPA_SM 

TUSED_SCARC(8,NM)=TUSED_SCARC(8,NM)+SECOND()-TNOW_SMOOTHER3D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_SMOOTHER3D

1000 FORMAT ('     SCARC_SMOOTHER : #ite= ',i4,': level=',i4,': res=',e12.4)
2000 FORMAT ('     SCARC_SMOOTHER : #ite= ',i4,': level=',i4,': res=',e12.4,': rate=',e12.4)
END SUBROUTINE SCARC_SMOOTHER3D



 
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
 
IF (SCARC_DEBUG>=6) CALL SCARC_SHOW_LEVEL (X, 'SARC', 'X0    ',ILEVEL)
IF (SCARC_DEBUG>=6) CALL SCARC_SHOW_LEVEL (Y, 'SARC', 'Y0    ',ILEVEL)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform 2D-matrix-vector-multiplication
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DO K = 1, SL%KBAR
   DO I = 1, SL%IBAR
      IC = (K-1) * SL%IBAR + I
      Y (I, 1, K) =    A1 * (  SL%AG(IC, 1)*X(I  , 1, K  ) &
                             + SL%AG(IC, 2)*X(I  , 1, K-1) &
                             + SL%AG(IC, 3)*X(I-1, 1, K  ) &
                             + SL%AG(IC, 4)*X(I+1, 1, K  ) &
                             + SL%AG(IC, 5)*X(I  , 1, K+1))&
                     + A2 * Y (I, 1, K)
      IF (SCARC_DEBUG>=2) &
         WRITE(SCARC_LU,'(2i3,11e11.3)') I,K,Y(I,1,K),(SL%AG(IC,II),II=1,5),&
                                  X(I,1,K),X(I,1,K-1),X(I-1,1,K),X(I+1,1,K),X(I,1,K+1)
   ENDDO
ENDDO
 
IF (SCARC_DEBUG>=6) CALL SCARC_SHOW_LEVEL (X, 'SARC', 'X1    ',ILEVEL)
IF (SCARC_DEBUG>=6) CALL SCARC_SHOW_LEVEL (Y, 'SARC', 'Y1    ',ILEVEL)
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform data exchange at interior boundaries
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF ((ITYPE == NTYPE_GLOBAL) .AND. (NMESHES > 1)) THEN
   NREQ_FACE = 0
   CALL SCARC_RECEIVE  (NCOM_MATV,  ILEVEL)   
   CALL SCARC_EXCHANGE (NCOM_MATV,  ILEVEL, IMV1, IMV2)
ENDIF
 
IF (SCARC_DEBUG>=6) CALL SCARC_SHOW_LEVEL (X, 'SARC', 'X2    ',ILEVEL)
IF (SCARC_DEBUG>=6) CALL SCARC_SHOW_LEVEL (Y, 'SARC', 'Y2    ',ILEVEL)
 
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
INTEGER :: I, J, K, IC, IMV1, IMV2, ILEVEL
 
 
TNOW_MATVEC3D = SECOND()
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform 3D-matrix-vector-multiplication
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DO K = 1, KBAR
   DO J = 1, JBAR
      DO I = 1, IBAR
         IC = (K-1) * IBAR * JBAR + (J-1) * IBAR + I
         Y(I, J, K) = A1 * ( SL%AG(IC, 1) *X (I  , J  , K-1)   &
                         +   SL%AG(IC, 2) *X (I  , J-1, K  )   &
                         +   SL%AG(IC, 3) *X (I-1, J  , K  )   &
                         +   SL%AG(IC, 4) *X (I  , J  , K  )   &
                         +   SL%AG(IC, 5) *X (I+1, J  , K  )   &
                         +   SL%AG(IC, 6) *X (I  , J+1, K  )   &
                         +   SL%AG(IC, 7) *X (I  , J  , K+1) ) &
                    + A2 * Y (I, J, K)
      ENDDO
   ENDDO
ENDDO
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform data exchange at interior boundaries
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF ((ITYPE == NTYPE_GLOBAL) .AND. (NMESHES > 1)) THEN
   NREQ_FACE = 0
   CALL SCARC_RECEIVE  (NCOM_MATV,  ILEVEL)  
   CALL SCARC_EXCHANGE (NCOM_MATV,  ILEVEL, IMV1, IMV2)
ENDIF
 
 
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
1000 FORMAT('X_LO(',I3,',',I3,',',I3,')=',f12.6,6i3)
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
1000 FORMAT('X(',I3,',',I3,',',I3,')=',f12.6)
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
INTEGER :: ILEVEL, IW_FACE, I, J, K, NM
REAL(EB):: TNOW_BDRY_RESIDUUM


TNOW_BDRY_RESIDUUM = SECOND()
SL => S%SLEVEL(ILEVEL)

CALL SCARC_SHOW_LEVEL (F, 'SARC',  'FRES1 ',ILEVEL)
DO IW_FACE=1,SL%NW_FACE

   IF (SL%IJKW_FACE(9,IW_FACE)==0) THEN

      I=SL%IJKW_FACE(1,IW_FACE)
      J=SL%IJKW_FACE(2,IW_FACE)
      K=SL%IJKW_FACE(3,IW_FACE)

      F(I,J,K)=0.0_EB

      if (SCARC_DEBUG>=4) WRITE(SCARC_LU,*) 'F(',I,',',J,',',K,')=',F(I,J,K)

   ENDIF

ENDDO
   CALL SCARC_SHOW_LEVEL (F, 'SARC',  'FRES2 ',ILEVEL)

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
IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,*) 'SCARC_SCALPROD2D=',SP

TUSED_SCARC(16,NM)=TUSED_SCARC(16,NM)+SECOND()-TNOW_SCALPROD2D
TUSED_SCARC( 0,NM)=TUSED_SCARC( 0,NM)+SECOND()-TNOW_SCALPROD2D
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
IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,*) 'SCARC_SCALPROD3D=',SP

TUSED_SCARC(16,NM)=TUSED_SCARC(16,NM)+SECOND()-TNOW_SCALPROD2D
TUSED_SCARC( 0,NM)=TUSED_SCARC( 0,NM)+SECOND()-TNOW_SCALPROD2D
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
   ENDDO
ENDDO
 
!!! scale with number of cells
IF (ITYPE == NTYPE_LOCAL .OR. NMESHES==1) THEN                
   SP = Sqrt (SP/REAL(SL%NCELLS_LOCAL, EB))                   ! scale with local  number of cells
ELSE IF (ITYPE == NTYPE_GLOBAL) THEN                          
   IF (USE_MPI) THEN
      CALL MPI_ALLREDUCE (SP, SPG, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
      SP = Sqrt (SPG/REAL(SL%NCELLS_GLOBAL, EB))                 ! scale with global number of cells
   ELSE
      WRITE(*,*) 'Serial version not yet implemented'
      stop
   ENDIF
ELSE
   WRITE (*,*) 'Wrong type for SCARC_L2NORM2D ', ITYPE
ENDIF

SCARC_L2NORM2D = SP
IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,*) 'SCARC_L2NORM2D=',SP

TUSED_SCARC(17,NM)=TUSED_SCARC(17,NM)+SECOND()-TNOW_L2NORM2D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_L2NORM2D
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

!!! build local scalar product (x,x)
SP = 0.0_EB
DO K = 1, SL%KBAR
   DO J = 1, SL%JBAR
      DO I = 1, SL%IBAR
         SP = SP + X (I, J, K) * X (I, J, K)
      ENDDO
   ENDDO
ENDDO
 
!!! scale with number of cells
IF (ITYPE == NTYPE_LOCAL .OR. NMESHES==1) THEN                
   SP = Sqrt (SP/REAL(SL%NCELLS_LOCAL, EB))                   ! scale with local number of cells
ELSE IF (ITYPE == NTYPE_GLOBAL) THEN                          
   IF (USE_MPI) THEN
      CALL MPI_ALLREDUCE (SP, SPG, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
      SP = Sqrt (SPG/REAL(SL%NCELLS_GLOBAL, EB))                 ! scale with global number of cells
   ELSE
      WRITE(*,*) 'Serial version not yet implemented'
      stop
   ENDIF
ELSE
   WRITE (*,*) 'Wrong type for SCARC_L2NORM3D ', ITYPE
ENDIF

SCARC_L2NORM3D = SP
IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,*) 'SCARC_L2NORM3D=',SP

TUSED_SCARC(17,NM)=TUSED_SCARC(17,NM)+SECOND()-TNOW_L2NORM3D
TUSED_SCARC( 0,NM)=TUSED_SCARC( 0,NM)+SECOND()-TNOW_L2NORM3D
END FUNCTION SCARC_L2NORM3D
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Jacobi preconditioner
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PRECON_JACOBI2D (Y, ILEVEL, NM)
REAL (EB), POINTER, DIMENSION (:, :, :) :: Y
REAL (EB) :: TNOW_JACOBI2D
INTEGER :: I, K, IC, ILEVEL, NM
 
SL => S%SLEVEL(ILEVEL)

TNOW_JACOBI2D = SECOND()
 
DO K = 1, SL%KBAR
   DO I = 1, SL%IBAR
      IC = (K-1) * SL%IBAR + I
      Y (I, 1, K) = Y (I, 1, K) / SL%AG (IC, 1)
   ENDDO
ENDDO
 
TUSED_SCARC(18,NM)=TUSED_SCARC(18,NM)+SECOND()-TNOW_JACOBI2D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_JACOBI2D
END SUBROUTINE SCARC_PRECON_JACOBI2D
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Jacobi preconditioner
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PRECON_JACOBI3D (Y, ILEVEL, NM)
REAL (EB), POINTER, DIMENSION (:, :, :) :: Y
REAL (EB) :: TNOW_JACOBI3D
INTEGER :: I, J, K, IC, ILEVEL, NM
 
SL => S%SLEVEL(ILEVEL)

TNOW_JACOBI3D = SECOND()
 
DO K = 1, SL%KBAR
   DO J = 1, SL%JBAR
      DO I = 1, SL%IBAR
         IC = (K-1) * SL%IBAR * SL%JBAR + (J-1) * SL%IBAR + I
         Y (I, J, K) = Y (I, J, K) / SL%AG (IC, 1)
      ENDDO
   ENDDO
ENDDO
 
TUSED_SCARC(18,NM)=TUSED_SCARC(18,NM)+SECOND()-TNOW_JACOBI3D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_JACOBI3D
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
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Initialize GSTRIX preconditioner
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PRECON_INIT_GSTRIX2D(ILEVEL, NM)

INTEGER :: NM,ILEVEL, IC, NC, IERR
REAL(EB):: TNOW_GSTRIX2D_INIT
 
TNOW_GSTRIX2D_INIT = SECOND()
 
IERR=0

SL => S%SLEVEL(ILEVEL)

ALLOCATE (SL%DD(1:SL%NCELLS_LOCAL), STAT=IERR)
CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX2D', 'SL%DD', IERR)

ALLOCATE (SL%DU(1:SL%NCELLS_LOCAL), STAT=IERR)
CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX2D', 'SL%DU', IERR)

ALLOCATE (SL%DL(1:SL%NCELLS_LOCAL), STAT=IERR)
CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX2D', 'SL%DL', IERR)

ALLOCATE (SL%LD(1:SL%NCELLS_LOCAL), STAT=IERR)
CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX2D', 'SL%LD', IERR)

ALLOCATE (SL%DAUX(1:SL%NCELLS_LOCAL), STAT=IERR)
CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX2D', 'SL%DAUX', IERR)


NC = SL%NCELLS_LOCAL
 
!!! save different diagonals of matrix ('D' corresponds to main, 'L' to lower, 'U' to upper)
DO IC = 1, NC
   SL%DD(IC) = SL%AG(IC, 1)         ! diagonal
   SL%LD(IC) = SL%AG(IC, 2)         ! lower diagonal in z-direction
   SL%DL(IC) = SL%AG(IC, 3)         ! lower diagonal in x-direction
   SL%DU(IC) = SL%AG(IC, 4)         ! upper diagonal in x-direction
ENDDO
 
!!! perform LU-factorization of matrix AG according to bandwise storage technique
DO IC = 2, NC
   SL%DL (IC) = SL%DL(IC) / SL%DD(IC-1)
   SL%DD (IC) = SL%DD(IC) - SL%DL(IC) * SL%DU(IC-1)
ENDDO
 
!!! replace diagonal values diag(i) by 1/diag(i) and multiply DU with 1/diag(i)
DO IC = 1, NC
   SL%DD (IC) = 1.0_EB / SL%DD(IC)
   SL%DU (IC) = SL%DD(IC) * SL%DU(IC)
   IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,1000) ILEVEL,IC,SL%DU(IC),IC,SL%DD(IC)
ENDDO

TUSED_SCARC(32,NM)=TUSED_SCARC(32,NM)+SECOND()-TNOW_GSTRIX2D_INIT
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_GSTRIX2D_INIT
1000 FORMAT('GSTRIX2D_INIT:',i3,': SL%DU(',i3,')=',f12.6,': SL%DD(',i3,')=',f12.6)
END SUBROUTINE SCARC_PRECON_INIT_GSTRIX2D
 
 
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! GSTRIX2D preconditioner
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PRECON_GSTRIX2D (DY, ILEVEL, NM)
REAL (EB), POINTER, DIMENSION (:, :, :) :: DY
INTEGER :: I, K, IC, ILEVEL, NM, N, M
REAL (EB) :: TNOW_GSTRIX2D
 
TNOW_GSTRIX2D = SECOND()

SL => S%SLEVEL(ILEVEL)
 
IF (SCARC_DEBUG>=2) THEN
   WRITE(SCARC_LU,*) '===================== STARTING PRECON_GSTRIX2D'
   WRITE(SCARC_LU,*) 'N=',N
   WRITE(SCARC_LU,*) 'M=',M
   WRITE(SCARC_LU,*) 'BEFORE PRECON_GSTRIX: Y'
   CALL SCARC_SHOW(DY, 'SARC', 'Y     ')
ENDIF


! backward elimination of first SL%IBAR unkowns (may be solved by tridiagonal system)
DO I=2,SL%IBAR
   DY(I,1,1) = DY(I,1,1)-SL%DU(I-1)*DY(I-1,1,1)
ENDDO
DO I=1,SL%IBAR
   DY(I,1,1) = DY(I,1,1)*SL%DD(I)
ENDDO
DO I=SL%IBAR-1,1,-1
   DY(I,1,1) = DY(I,1,1)-SL%DL(I)*DY(I+1,1,1)
ENDDO
  

! backward elimination of following unknowns (here the subdiagonal in z-direction must be taken into account)
DO K=2,SL%KBAR

   IC = (K-1)*SL%IBAR + I

   ! bring subdiagonal in z-direction to right hand side (corresponding solution component already known)
   DO I=1,SL%IBAR
      DY(I,1,K) = DY(I,1,K) - SL%LD(IC)*DY(I,1,K-1) 
   ENDDO

   ! perform elimination of matrix lines corresponding to K
   DO I=2,SL%IBAR
      DY(I,1,K) = DY(I,1,K)-SL%DU(IC-1)*DY(I-1,1,K)
   ENDDO
   DO I=1,SL%IBAR
      DY(I,1,K) = DY(I,1,K)*SL%DD(IC)
   ENDDO
   DO I=SL%IBAR-1,1,-1
      DY(I,1,K) = DY(I,1,K)-SL%DL(IC)*DY(I+1,1,K)
   ENDDO
  
ENDDO

IF (SCARC_DEBUG >= 2) THEN
   WRITE(SCARC_LU,*) 'AFTER  PRECON_GSTRIX: Y'
   CALL SCARC_SHOW(DY, 'SARC',  'Y     ')
ENDIF
 
TUSED_SCARC(31,NM)=TUSED_SCARC(20,NM)+SECOND()-TNOW_GSTRIX2D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_GSTRIX2D
 
END SUBROUTINE SCARC_PRECON_GSTRIX2D


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Initialize GSTRIX3D preconditioner
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PRECON_INIT_GSTRIX3D(ILEVEL, NM)

INTEGER :: NM,ILEVEL, IC, NC, IERR
REAL(EB):: TNOW_GSTRIX3D_INIT
 
TNOW_GSTRIX3D_INIT = SECOND()
 
IERR=0

write(*,*) 'MUSS NOCH AUF 3D angepasst werden !!!'
WRITE(SCARC_LU,*) 'MUSS NOCH AUF 3D angepasst werden !!!'
stop

SL => S%SLEVEL(ILEVEL)

ALLOCATE (SL%DD(1:SL%NCELLS_LOCAL), STAT=IERR)
CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX3D', 'SL%DD', IERR)

ALLOCATE (SL%DU(1:SL%NCELLS_LOCAL), STAT=IERR)
CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX3D', 'SL%DU', IERR)

ALLOCATE (SL%DL(1:SL%NCELLS_LOCAL), STAT=IERR)
CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX3D', 'SL%DL', IERR)

ALLOCATE (SL%LD(1:SL%NCELLS_LOCAL), STAT=IERR)
CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX3D', 'SL%LD', IERR)

ALLOCATE (SL%DAUX(1:SL%NCELLS_LOCAL), STAT=IERR)
CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX3D', 'SL%DAUX', IERR)


NC = SL%NCELLS_LOCAL
 
!!! save different diagonals of matrix ('D' corresponds to main, 'L' to lower, 'U' to upper)
DO IC = 1, NC
   SL%DD(IC) = SL%AG(IC, 1)         ! diagonal
   SL%LD(IC) = SL%AG(IC, 2)         ! lower diagonal in z-direction
   SL%DL(IC) = SL%AG(IC, 3)         ! lower diagonal in x-direction
   SL%DU(IC) = SL%AG(IC, 4)         ! upper diagonal in x-direction
ENDDO
 
!!! perform LU-factorization of matrix AG according to bandwise storage technique
DO IC = 2, NC
   SL%DL (IC) = SL%DL(IC) / SL%DD(IC-1)
   SL%DD (IC) = SL%DD(IC) - SL%DL(IC) * SL%DU(IC-1)
ENDDO
 
!!! replace diagonal values diag(i) by 1/diag(i) and multiply DU with 1/diag(i)
DO IC = 1, NC
   SL%DD (IC) = 1.0_EB / SL%DD(IC)
   SL%DU (IC) = SL%DD(IC) * SL%DU(IC)
   IF (SCARC_DEBUG>=2) WRITE(SCARC_LU,1000) ILEVEL,IC,SL%DU(IC),IC,SL%DD(IC)
ENDDO

TUSED_SCARC(32,NM)=TUSED_SCARC(32,NM)+SECOND()-TNOW_GSTRIX3D_INIT
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_GSTRIX3D_INIT
1000 FORMAT('GSTRIX3D_INIT:',i3,': SL%DU(',i3,')=',f12.6,': SL%DD(',i3,')=',f12.6)
END SUBROUTINE SCARC_PRECON_INIT_GSTRIX3D
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! GSTRIX3D preconditioner
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PRECON_GSTRIX3D (DY, ILEVEL, NM)
REAL (EB), POINTER, DIMENSION (:, :, :) :: DY
INTEGER :: I, K, IC, ILEVEL, NM, N, M
REAL (EB) :: TNOW_GSTRIX3D
 
TNOW_GSTRIX3D = SECOND()

SL => S%SLEVEL(ILEVEL)
 
IF (SCARC_DEBUG>=2) THEN
   WRITE(SCARC_LU,*) '===================== STARTING PRECON_GSTRIX3D'
   WRITE(SCARC_LU,*) 'N=',N
   WRITE(SCARC_LU,*) 'M=',M
   WRITE(SCARC_LU,*) 'BEFORE PRECON_GSTRIX: Y'
   CALL SCARC_SHOW(DY, 'SARC',  'Y     ')
ENDIF

write(*,*) 'MUSS NOCH AUF 3D angepasst werden !!!'
WRITE(SCARC_LU,*) 'MUSS NOCH AUF 3D angepasst werden !!!'
stop

! backward elimination of first SL%IBAR unkowns (may be solved by tridiagonal system)
DO I=2,SL%IBAR
   DY(I,1,1) = DY(I,1,1)-SL%DU(I-1)*DY(I-1,1,1)
ENDDO
DO I=1,SL%IBAR
   DY(I,1,1) = DY(I,1,1)*SL%DD(I)
ENDDO
DO I=SL%IBAR-1,1,-1
   DY(I,1,1) = DY(I,1,1)-SL%DL(I)*DY(I+1,1,1)
ENDDO
  

! backward elimination of following unknowns (here the subdiagonal in z-direction must be taken into account)
DO K=2,SL%KBAR

   IC = (K-1)*SL%IBAR + I

   ! bring subdiagonal in z-direction to right hand side (corresponding solution component already known)
   DO I=1,SL%IBAR
      DY(I,1,K) = DY(I,1,K) - SL%LD(IC)*DY(I,1,K-1) 
   ENDDO

   ! perform elimination of matrix lines corresponding to K
   DO I=2,SL%IBAR
      DY(I,1,K) = DY(I,1,K)-SL%DU(IC-1)*DY(I-1,1,K)
   ENDDO
   DO I=1,SL%IBAR
      DY(I,1,K) = DY(I,1,K)*SL%DD(IC)
   ENDDO
   DO I=SL%IBAR-1,1,-1
      DY(I,1,K) = DY(I,1,K)-SL%DL(IC)*DY(I+1,1,K)
   ENDDO
  
ENDDO

IF (SCARC_DEBUG >= 2) THEN
   WRITE(SCARC_LU,*) 'AFTER  PRECON_GSTRIX: Y'
   CALL SCARC_SHOW(DY, 'SARC',  'Y     ')
ENDIF
 
TUSED_SCARC(31,NM)=TUSED_SCARC(20,NM)+SECOND()-TNOW_GSTRIX3D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_GSTRIX3D
 
END SUBROUTINE SCARC_PRECON_GSTRIX3D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Perform global MG-method as preconditioner in 2D
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_MG_PRECON2D (Y,ILEVEL,NM)

REAL (EB), POINTER, DIMENSION (:, :, :) :: Y
INTEGER :: NM, ITE, ITE0, ICYCLE, ILEVEL
INTEGER :: I, K, IBAR0, JBAR0, KBAR0
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
SLHI%X2 = 0.0_EB
 
DO K = 1, KBAR                                ! most possible redundant --- to be checked (use Y instead)
   DO I = 1, IBAR
      SLHI%F2 (I, 1, K) = Y (I, 1, K)
   ENDDO
ENDDO
 
IF (SCARC_DEBUG >= 3) THEN
   WRITE(SCARC_LU,*) 'Starting SCARC_MG_PRECON2D'
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
 
IBAR0=SLHI%IBAR
JBAR0=SLHI%JBAR
KBAR0=SLHI%KBAR

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
IF (NM == 1)          WRITE(*,1000) 0, ILEVEL, SL%RESIN_MG
 
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

      CALL SCARC_SMOOTHER2D(SLHI%X2, SLHI%F2, SLHI%D2, BMATVEC, ILEVEL, NM)

      !!! print level information - still missing
      !CALL SCARC_INFO_LEVEL(ILEVEL)

      !!! perform restriction to coarser grid: F:= rest(D)
      CALL SCARC_SHOW_LEVEL (SLHI%D2, 'SARC',  'RESTD ',ILEVEL)
      CALL SCARC_RESTRICTION2D(SLLO%F2, SLHI%D2, ILEVEL-1, NM)
      CALL SCARC_SHOW_LEVEL (SLLO%F2, 'SARC',  'RESTF ',ILEVEL-1)

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
   CALL SCARC_COARSE2D(2, ILEVEL, NM)


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
      CALL SCARC_SMOOTHER2D(SLHI%X2, SLHI%F2, SLHI%D2, BMATVEC, ILEVEL, NM)


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
   IF (SCARC_DEBUG >= 1) write (*,1000) ITE, ILEVEL, SL%RES_MG
   IF (NM == 1)          write (*,1000) ITE, ILEVEL, SL%RES_MG
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
!IF (NM == 1)          WRITE (*,2000) ITE0, ILEVEL, SL%RES_MG, SL%CAPPA_MG
 
DO K = 1, SLHI%KBAR
   DO I = 1, SLHI%IBAR
      Y (I, 1, K) = SLHI%X2(I, 1, K)
   ENDDO
ENDDO
 
TUSED_SCARC(7,NM)=TUSED_SCARC(7,NM)+SECOND()-TNOW_PRECON_MG2D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_PRECON_MG2D

1000 FORMAT (/,'     SCARC_MG_PRECON: #ite= ',i4,': level=',i4,': res=',e12.4,/)
2000 FORMAT (/,'     SCARC_MG_PRECON: #ite= ',i4,': level=',i4,': res=',e12.4,': rate=',e12.4,/)
END SUBROUTINE SCARC_MG_PRECON2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Perform global MG-method as preconditioner in 3D
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_MG_PRECON3D (Y,ILEVEL,NM)

REAL (EB), POINTER, DIMENSION (:, :, :) :: Y
INTEGER :: NM, ITE, ITE0, ICYCLE, ILEVEL
INTEGER :: I, K, IBAR0, JBAR0, KBAR0
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
 
DO K = 1, KBAR                                ! most possible redundant --- to be checked (use Y instead)
   DO I = 1, IBAR
      SLHI%F2 (I, 1, K) = Y (I, 1, K)
   ENDDO
ENDDO
 
IF (SCARC_DEBUG >= 3) THEN
   WRITE(SCARC_LU,*) 'Starting SCARC_MG_PRECON3D'
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
 
IBAR0=SLHI%IBAR
JBAR0=SLHI%JBAR
KBAR0=SLHI%KBAR

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
SLHI%RESIN_MG = SCARC_L2NORM3D (SLHI%D2, ILEVEL, NM, NTYPE_GLOBAL)
IF (SCARC_DEBUG >= 1) WRITE(SCARC_LU,1000) 0, ILEVEL, SL%RESIN_MG
IF (NM == 1)          WRITE(*,1000) 0, ILEVEL, SL%RESIN_MG
 
!WRITE(*,*) 'Info muss noch an Master geschickt werden'
!CALL SCARC_MASTER_INFO (SLHI%ITE, SLHI%RESIN_MG)

 
!!! ---------------------------------------------------------------------------
!!! start MG-iteration
!!! ---------------------------------------------------------------------------
MG_LOOP3D: DO ITE = 1, SCARC_MG_NIT
 
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

      CALL SCARC_SMOOTHER3D(SLHI%X2, SLHI%F2, SLHI%D2, BMATVEC, ILEVEL, NM)

      !!! print level information - still missing
      !CALL SCARC_INFO_LEVEL(ILEVEL)

      !!! perform restriction to coarser grid: F:= rest(D)
      CALL SCARC_SHOW_LEVEL (SLHI%D2, 'SARC',  'RESTD ',ILEVEL)
      CALL SCARC_RESTRICTION3D(SLLO%F2, SLHI%D2, ILEVEL-1, NM)
      CALL SCARC_SHOW_LEVEL (SLLO%F2, 'SARC',  'RESTF ',ILEVEL-1)

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
   CALL SCARC_COARSE3D(2, ILEVEL, NM)


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

      !!! set exterior boundary data of residuum to zero - necessary ?
      !CALL SCARC_BDRY_RESIDUUM(NM, SLHI%D, ILEVEL)
 
      !!! set new solution
      SLHI%X2 = SLHI%D2 + SLHI%X2

      !!! select smoother for postsmoothing
      CALL SCARC_SMOOTHER3D(SLHI%X2, SLHI%F2, SLHI%D2, BMATVEC, ILEVEL, NM)


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
   IF (SCARC_DEBUG >= 1) write (*,1000) ITE, ILEVEL, SL%RES_MG
   IF (NM == 1)          write (*,1000) ITE, ILEVEL, SL%RES_MG
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
!IF (NM == 1)          WRITE (*,2000) ITE0, ILEVEL, SL%RES_MG, SL%CAPPA_MG
 
DO K = 1, SLHI%KBAR
   DO I = 1, SLHI%IBAR
      Y (I, 1, K) = SLHI%X2(I, 1, K)
   ENDDO
ENDDO
 
TUSED_SCARC(7,NM)=TUSED_SCARC(7,NM)+SECOND()-TNOW_PRECON_MG3D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_PRECON_MG3D

1000 FORMAT (/,'     SCARC_MG_PRECON: #ite= ',i4,': level=',i4,': res=',e12.4,/)
2000 FORMAT (/,'     SCARC_MG_PRECON: #ite= ',i4,': level=',i4,': res=',e12.4,': rate=',e12.4,/)
END SUBROUTINE SCARC_MG_PRECON3D


 
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
INTEGER :: ILEN_FACE
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
      RECEIVE_FACE_IF: IF (SNML%NIC_FACE(NM,NOM)/=0 .AND. SNML%NIC_FACE(NOM,NM)/=0) THEN
    
         OSNM  => SCARC(NM)%OSCARC(NOM)                  ! corresponds to M3
         OSNML => SCARC(NM)%OSCARC(NOM)%SLEVEL(ILEVEL)   ! ... for the level 'ILEVEL'

         IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,*) '=========== RECEIVE: NM=',NM,', NOM=',NOM
         IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,*) 'OSNML zeigt auf SCARC(',NM,')%OSCARC(',NOM,')%SLEVEL(',ILEVEL,')'

         TAG_FACE = TAGS_FACE(NM,NOM)
    

         !!! Initialize the communication structures for receiving face data
         INIT_FACE_IF: IF (CODE==NCOM_INIT) THEN

            NREQ_FACE = NREQ_FACE+1
            IF (USE_MPI) CALL MPI_IRECV(OSNML%IJKW_FACE(1,1),17*OSNML%NW_FACE,MPI_INTEGER,SNODE, &
                                        TAG_FACE,MPI_COMM_WORLD,REQ_FACE(NREQ_FACE),IERR)

            IF (SCARC_DEBUG>=6) THEN
               WRITE(SCARC_LU,*) 'RECEIVE: INIT: LEVEL=',ILEVEL
               WRITE(SCARC_LU,*) 'RECEIVE: INIT: NW_FACE=',OSNML%NW_FACE
               WRITE(SCARC_LU,*) 'RECEIVE: INIT: LEN(IJKW_FACE)=',17*OSNML%NW_FACE
               WRITE(SCARC_LU,'(17i4)') SNML%IJKW_FACE
               WRITE(SCARC_LU,*) 'RECEIVE: INIT: SIZE(IJKW_FACE):  ',SIZE(OSNML%IJKW_FACE)
               WRITE(SCARC_LU,*) 'RECEIVING IJKW_FACE(',ILEVEL,'): '
               DO IW= 1,OSNML%NW_FACE
                  WRITE(SCARC_LU,'(17i4)') (OSNML%IJKW_FACE(II,IW),II=1,17)
               ENDDO
            ENDIF

            IF (SNML%NIC_FACE(NM, NOM) > 0) THEN
               !ILEN_FACE=(MAX(SNML%NIC_FACE(NM, NOM), SNML%NIC_FACE(NOM, NM))+2)*2+1
               ILEN_FACE=(MAX(SNML%NIC_FACE(NM, NOM), SNML%NIC_FACE(NOM, NM)))*2+1
               IF (SCARC_DEBUG>=8) WRITE(SCARC_LU,*) 'RECEIVE: ILEN_FACE=',ILEN_FACE
               ALLOCATE (OSNML%RECV_FACE(ILEN_FACE))
               OSNML%RECV_FACE = 0.0_EB
            ENDIF

         ENDIF INIT_FACE_IF
   
         !!! Perform exchange for matrix-vector multiplication (only face exchange needed)
         MATV_FACE_IF: IF (CODE==NCOM_MATV) THEN

            NREQ_FACE = NREQ_FACE+1
            IF (USE_MPI) CALL MPI_IRECV(OSNML%RECV_FACE(1),SIZE(OSNML%RECV_FACE),MPI_DOUBLE_PRECISION,&
                                        SNODE,TAG_FACE,MPI_COMM_WORLD,REQ_FACE(NREQ_FACE),IERR)

            IF (SCARC_DEBUG>=8) THEN
               WRITE(SCARC_LU,*) 'RECEIVE: MATV: receive length',SIZE(OSNML%RECV_FACE)
               WRITE(SCARC_LU,*) 'RECEIVE: MATV: SNODE=',SNODE 
               WRITE(SCARC_LU,*) 'RECEIVE: MATV: TAG_FACE=',TAG_FACE 
               WRITE(SCARC_LU,*) 'RECEIVE: MATV: NREQ_FACE=',NREQ_FACE 
               WRITE(SCARC_LU,*) 'RECEIVE: MATV: REQ_FACE =',REQ_FACE(NREQ_FACE) 
            ENDIF
            IF (SCARC_DEBUG>=6) THEN
               WRITE(SCARC_LU,*) 'RECEIVE: MATV: RECV_FACE(1)'
               WRITE(SCARC_LU,'(2e20.8)') SCARC(NM)%OSCARC(NOM)%SLEVEL(1)%RECV_FACE
               WRITE(SCARC_LU,*) 'SIZE(OSNML%RECV_FACE(1))=',SIZE(SCARC(NM)%OSCARC(NOM)%SLEVEL(1)%RECV_FACE)
               WRITE(SCARC_LU,*) 'RECEIVED IJKW_FACE(',ILEVEL,'): ', OSNML%NW_FACE
               DO IW= 1,OSNML%NW_FACE
                  WRITE(SCARC_LU,'(17i4)') (OSNML%IJKW_FACE(II,IW),II=1,17)
               ENDDO
            ENDIF
            IF (SCARC_DEBUG>=8) THEN
               WRITE(SCARC_LU,*) 'RECEIVE: MATV: RECV_FACE(2)'
               WRITE(SCARC_LU,'(2e20.8)') SCARC(NM)%OSCARC(NOM)%SLEVEL(2)%RECV_FACE
               WRITE(SCARC_LU,*) 'SIZE(OSNML%RECV_FACE(1))=',SIZE(SCARC(NM)%OSCARC(NOM)%SLEVEL(2)%RECV_FACE)
            ENDIF

         ENDIF MATV_FACE_IF
   
         !!! Perform full exchange including edge (3D)
         FULL_FACE_IF: IF (CODE==NCOM_FULL) THEN

            NREQ_FACE = NREQ_FACE+1
            IF (SCARC_DEBUG>=8) WRITE(SCARC_LU,*) 'RECEIVE:', SIZE(OSNML%RECV_FACE)
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
INTEGER :: I, J, K, LL, IW, IWW
INTEGER :: IERR, II, JJ, KK
INTEGER :: ILEN_FACE
INTEGER :: TAG_FACE

REAL(EB):: ASUB, ZSUM,yold
REAL(EB):: TNOW_EXCHANGE


TYPE (SCARC_TYPE),  POINTER ::  SNM,   SNOM
TYPE (OSCARC_TYPE), POINTER :: OSNM, OSNOM
TYPE (SCARC_LEVEL_TYPE),  POINTER ::  SNML,  SNOML
TYPE (OSCARC_LEVEL_TYPE), POINTER :: OSNML, OSNOML

TNOW_EXCHANGE = SECOND()

SNM => SCARC(NM)
 
IERR=0
IF (SCARC_DEBUG>=8) WRITE(SCARC_LU,*) 'SCARC_EXCHANGE: ', IMV1, IMV2

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
      FACE_IF: IF (SNML%NIC_FACE(NOM,NM)/=0 .AND. SNML%NIC_FACE(NM,NOM)/=0)  THEN
    
         OSNM  => SCARC(NM)%OSCARC(NOM)                 ! corresponds to M3
         OSNML => SCARC(NM)%OSCARC(NOM)%SLEVEL(ILEVEL)  ! ... for the level 'ILEVEL'

         IF (SCARC_DEBUG>=6) THEN
            WRITE(SCARC_LU,*) 'OSNML zeigt auf SCARC(',NM,')%OSCARC(',NOM,')%SLEVEL(',ILEVEL,')'
            WRITE(SCARC_LU,*) 'SIZE(IJKW_FACE(1)):  ',SIZE(SCARC(NM)%OSCARC(NOM)%SLEVEL(1)%IJKW_FACE)
            WRITE(SCARC_LU,*) 'SIZE(IJKW_FACE(2)):  ',SIZE(SCARC(NM)%OSCARC(NOM)%SLEVEL(2)%IJKW_FACE)
         ENDIF

         TAG_FACE = TAGS_FACE(NM,NOM)
    
         !!! Initialize the communication structures for sending data
         EXCHANGE_INIT_FACE_IF: IF (CODE==NCOM_INIT) THEN
    
            ! initialize communication structures for face exchange
            IF (RNODE/=SNODE) THEN

               NREQ_FACE = NREQ_FACE+1
               IF (USE_MPI) CALL MPI_ISEND(SNML%IJKW_FACE(1,1),17*SNML%NW_FACE,MPI_INTEGER,SNODE, &
                                           TAG_FACE,MPI_COMM_WORLD,REQ_FACE(NREQ_FACE),IERR)

               IF (SCARC_DEBUG>=6) THEN
                 WRITE(SCARC_LU,*) '============================================================='
                 WRITE(SCARC_LU,*) 'EXCHANGE: INIT: LEVEL=',ILEVEL
                 WRITE(SCARC_LU,*) 'EXCHANGE: INIT: LEN(IJKW_FACE)=',17*SNML%NW_FACE
                 WRITE(SCARC_LU,'(17i4)') SNML%IJKW_FACE
                 WRITE(SCARC_LU,*) 'EXCHANGE: INIT: SIZE(IJKW_FACE):  ',SIZE(SNML%IJKW_FACE)
               WRITE(SCARC_LU,*) 'SENDING IJKW_FACE(',ILEVEL,'): '
               DO IW= 1,SNML%NW_FACE
                  WRITE(SCARC_LU,'(17i4)') (SNML%IJKW_FACE(II,IW),II=1,17)
               ENDDO
               ENDIF

               !ILEN_FACE=(MAX(SNML%NIC_FACE(NM, NOM), SNML%NIC_FACE(NOM, NM))+2)*2+1   ! extended
               ILEN_FACE=(MAX(SNML%NIC_FACE(NM, NOM), SNML%NIC_FACE(NOM, NM)))*2+1 
               IF (SCARC_DEBUG>=8) WRITE(SCARC_LU,*) 'EXCHANGE: ILEN_FACE=',ILEN_FACE
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
               WRITE(SCARC_LU,*) 'NW_FACE: ',OSNML%NW_FACE
               WRITE(SCARC_LU,*) 'SIZE(IJKW_FACE):  ',SIZE(OSNML%IJKW_FACE)
               WRITE(SCARC_LU,*) 'SHAPE(IJKW_FACE): ',SHAPE(OSNML%IJKW_FACE)
               WRITE(SCARC_LU,*) 'IJKW_FACE: '
               DO IW= 1,OSNML%NW_FACE
                  WRITE(SCARC_LU,'(17i4)') (OSNML%IJKW_FACE(II,IW),II=1,17)
               ENDDO
            ENDIF
 
            LL = 0
            IWW = 0
IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,*) OSNML%NW_FACE
            PACK_SEND_FACE0: DO IW=1,OSNML%NW_FACE
IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,'(a,i3,a,17i4)') ' IW=',IW,':',(OSNML%IJKW_FACE(II,IW),II=1,17)
               IF (OSNML%IJKW_FACE(9,IW)/=NM) CYCLE PACK_SEND_FACE0
               DO KK=OSNML%IJKW_FACE(12,IW),OSNML%IJKW_FACE(15,IW)
                  DO JJ=OSNML%IJKW_FACE(11,IW),OSNML%IJKW_FACE(14,IW)
                     DO II=OSNML%IJKW_FACE(10,IW),OSNML%IJKW_FACE(13,IW)
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
                           CASE(NMV_D2)
                              OSNML%SEND_FACE(LL+2) = SNML%D2(II,JJ,KK)
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
if (SCARC_DEBUG>=6) THEN
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
   
      !!! Perform full exchange including edge (3D!) and vertex information
         EXCHANGE_FULL_FACE_IF: IF (CODE==NCOM_FULL) THEN
            IF (RNODE/=SNODE) THEN

            LL = 0
            IWW = 0
            PACK_SEND_FACE: DO IW=1,OSNML%NW_FACE
               IF (OSNML%IJKW_FACE(9,IW)/=NM) CYCLE PACK_SEND_FACE
               DO KK=OSNML%IJKW_FACE(12,IW),OSNML%IJKW_FACE(15,IW)
                  DO JJ=OSNML%IJKW_FACE(11,IW),OSNML%IJKW_FACE(14,IW)
                     DO II=OSNML%IJKW_FACE(10,IW),OSNML%IJKW_FACE(13,IW)
                        IWW = IWW + 1
                        OSNML%SEND_FACE(LL+1) = REAL(IW,EB)
                        OSNML%SEND_FACE(LL+2) = SNML%Z(II,JJ,KK)
                        LL = LL+2
!if (SCARC_DEBUG>=8) WRITE(SCARC_LU,*) 'SEND_FACE(',LL,')=',SNML%Z(II,JJ,KK), II, JJ, KK
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO PACK_SEND_FACE
            OSNML%SEND_FACE(IWW*2+1) = -999.0_EB
            NREQ_FACE=NREQ_FACE+1
            IF (USE_MPI) CALL MPI_ISEND(OSNML%SEND_FACE(1),IWW*2+1,MPI_DOUBLE_PRECISION,SNODE, &
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
IF (SCARC_DEBUG>=8) THEN
   WRITE(SCARC_LU,*) 'NREQ_FACE=',NREQ_FACE
   WRITE(SCARC_LU,*) 'REQ_FACE=',REQ_FACE(1:NREQ_FACE)
ENDIF
   
IF (USE_MPI) CALL MPI_WAITALL(NREQ_FACE,REQ_FACE(1:NREQ_FACE),STAT_FACE,IERR)

 
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


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! FACE communication
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      RECV_FACE_IF: IF (SNOML%NIC_FACE(NOM,NM)/=0 .AND. SNOML%NIC_FACE(NM,NOM)/=0) THEN

         OSNOM => SCARC(NOM)%OSCARC(NM)                   ! corresponds to M2
         OSNOML=> SCARC(NOM)%OSCARC(NM)%SLEVEL(ILEVEL)  ! ... for the level 'ILEVEL'

         !!! Extract data for matrix-vector communication
         RECV_FACE_MATV_IF: IF (CODE==NCOM_MATV) THEN
     
            IF (ABS(SNOM%IOR_FACE(NM,NOM))==1) THEN
               ASUB=SNOML%ASUBX
            ELSE IF (ABS(SNOM%IOR_FACE(NM,NOM))==2) THEN
               ASUB=SNOML%ASUBY
            ELSE IF (ABS(SNOM%IOR_FACE(NM,NOM))==3) THEN
               ASUB=SNOML%ASUBZ
            ENDIF

            IF (RNODE/=SNODE) THEN
               LL = 0

IF (SCARC_DEBUG>=6) THEN
   WRITE(SCARC_LU,*) 'OSNOML zeigt auf SCARC(',NOM,')%OSCARC(',NM,')%SLEVEL(',ILEVEL,')'
   WRITE(SCARC_LU,*) 'EXCHANGE: EXTRACT: SIZE(OSNOML%RECV_FACE)=',SIZE(OSNOML%RECV_FACE) 
   WRITE(SCARC_LU,*) 'EXCHANGE: EXTRACT: MATV: RECV_FACE(1)'
   WRITE(SCARC_LU,'(2e20.8)') SCARC(NOM)%OSCARC(NM)%SLEVEL(1)%RECV_FACE
ENDIF
IF (SCARC_DEBUG>=8) THEN
   WRITE(SCARC_LU,*) 'SIZE(OSNML%RECV_FACE(1))=',SIZE(SCARC(NOM)%OSCARC(NM)%SLEVEL(1)%RECV_FACE)
   WRITE(SCARC_LU,*) 'EXCHANGE: EXTRACT: MATV: RECV_FACE(2)'
   WRITE(SCARC_LU,'(2e20.8)') SCARC(NOM)%OSCARC(NM)%SLEVEL(2)%RECV_FACE
   WRITE(SCARC_LU,*) 'SIZE(OSNML%RECV_FACE(1))=',SIZE(SCARC(NOM)%OSCARC(NM)%SLEVEL(2)%RECV_FACE)
ENDIF

               UNPACK_RECV_FACE0: DO
                  IW = NINT(OSNOML%RECV_FACE(LL+1))
!IF (SCARC_DEBUG>=8) WRITE(SCARC_LU,*) '--------------------------- IW=',IW
                  IF (IW==-999) EXIT UNPACK_RECV_FACE0
                  ZSUM=0.0_EB
                  DO KK=SNOML%IJKW_FACE(12,IW),SNOML%IJKW_FACE(15,IW)
                     DO JJ=SNOML%IJKW_FACE(11,IW),SNOML%IJKW_FACE(14,IW)
                        DO II=SNOML%IJKW_FACE(10,IW),SNOML%IJKW_FACE(13,IW)
                           OSNOML%Z_FACE(II,JJ,KK) = OSNOML%RECV_FACE(LL+2)
IF (SCARC_DEBUG>=8) WRITE(SCARC_LU,*) 'Z_FACE(',II,',',JJ,',',KK,')=',OSNOML%Z_FACE(II,JJ,KK)
                           ZSUM=ZSUM+OSNOML%RECV_FACE(LL+2)
                           LL = LL+2
                        ENDDO
                     ENDDO
                  ENDDO
                  I=SNOML%IJKW_FACE(6,IW)
                  J=SNOML%IJKW_FACE(7,IW)
                  K=SNOML%IJKW_FACE(8,IW)
                  ISUM = (SNOML%IJKW_FACE(13,IW)-SNOML%IJKW_FACE(10,IW)+1) * &
                         (SNOML%IJKW_FACE(14,IW)-SNOML%IJKW_FACE(11,IW)+1) * &
                         (SNOML%IJKW_FACE(15,IW)-SNOML%IJKW_FACE(12,IW)+1)

IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,*) 'ASUB=',ASUB, ZSUM, ISUM, IMV2

                  SELECT CASE(IMV2)
                     CASE(NMV_Y)
                        yold=SNOML%Y(I,J,K)
                        SNOML%Y(I, J, K) = SNOML%Y(I, J, K) + ASUB * ZSUM/REAL(ISUM,EB)
IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,*) '    Y1(',I,',',J,',',K,')=',SNOML%Y(I,J,K), yold
                     CASE(NMV_G)
                        yold=SNOML%G(I,J,K)
                        SNOML%G(I, J, K) = SNOML%G(I, J, K) + ASUB * ZSUM/REAL(ISUM,EB)
IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,*) '    G1(',I,',',J,',',K,')=',SNOML%G(I,J,K), yold
                     CASE(NMV_R)
                        yold=SNOML%R(I,J,K)
                        SNOML%R(I, J, K) = SNOML%R(I, J, K) + ASUB * ZSUM/REAL(ISUM,EB)
IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,*) '    R1(',I,',',J,',',K,')=',SNOML%R(I,J,K), yold
                     CASE(NMV_D)
                        yold=SNOML%D(I,J,K)
                        SNOML%D(I, J, K) = SNOML%D(I, J, K) + ASUB * ZSUM/REAL(ISUM,EB)
IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,*) '    D1(',I,',',J,',',K,')=',SNOML%D(I,J,K), yold
                     CASE(NMV_X)
                        yold=SNOML%X(I,J,K)
                        SNOML%X(I, J, K) = SNOML%X(I, J, K) + ASUB * ZSUM/REAL(ISUM,EB)
IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,*) '    X1(',I,',',J,',',K,')=',SNOML%X(I,J,K), yold
                     CASE(NMV_X2)
                        yold=SNOML%X2(I,J,K)
                        SNOML%X2(I, J, K) = SNOML%X2(I, J, K) + ASUB * ZSUM/REAL(ISUM,EB)
IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,*) '   X21(',I,',',J,',',K,')=',SNOML%X2(I,J,K), yold
                     CASE(NMV_D2)
                        yold=SNOML%D2(I,J,K)
                        SNOML%D2(I, J, K) = SNOML%D2(I, J, K) + ASUB * ZSUM/REAL(ISUM,EB)
IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,*) '   D21(',I,',',J,',',K,')=',SNOML%D2(I,J,K), yold
                  END SELECT
               ENDDO UNPACK_RECV_FACE0
            ENDIF

         ENDIF RECV_FACE_MATV_IF


         !!! Extract data for full communication including edge (3D!) and vertex data
         RECV_FACE_FULL_IF: IF (CODE==NCOM_FULL) THEN
   
            IF (RNODE/=SNODE) THEN
               LL = 0
               UNPACK_RECV_FACE: DO
                  IW = NINT(OSNOML%RECV_FACE(LL+1))
                  IF (IW==-999) EXIT UNPACK_RECV_FACE
                  ZSUM=0.0_EB
                  DO KK=SNOML%IJKW_FACE(12,IW),SNOML%IJKW_FACE(15,IW)
                     DO JJ=SNOML%IJKW_FACE(11,IW),SNOML%IJKW_FACE(14,IW)
                        DO II=SNOML%IJKW_FACE(10,IW),SNOML%IJKW_FACE(13,IW)
                           OSNOML%Z_FACE(II,JJ,KK) = OSNOML%RECV_FACE(LL+2)
                           ZSUM=ZSUM+OSNOML%RECV_FACE(LL+2)
!if (SCARC_DEBUG>=8) WRITE(SCARC_LU,*) 'RECV_FACE(',II,',',JJ,',',KK,')=',OSNOML%Z_FACE(II,JJ,KK)
                           LL = LL+2
                        ENDDO
                     ENDDO
                  ENDDO
!if (SCARC_DEBUG>=8) WRITE(SCARC_LU,*) 'ZSUM=',ZSUM
                  I=SNOML%IJKW_FACE(1,IW)
                  J=SNOML%IJKW_FACE(2,IW)
                  K=SNOML%IJKW_FACE(3,IW)
                  ISUM = (SNOML%IJKW_FACE(13,IW)-SNOML%IJKW_FACE(10,IW)+1) * &
                         (SNOML%IJKW_FACE(14,IW)-SNOML%IJKW_FACE(11,IW)+1) * &
                         (SNOML%IJKW_FACE(15,IW)-SNOML%IJKW_FACE(12,IW)+1)
                  SNOML%Z(I, J, K) = ZSUM/REAL(ISUM,EB)
!if (SCARC_DEBUG>=8) WRITE(SCARC_LU,*) 'Z(',I,',',J,',',K,')=',SNOML%Z(I,J,K)
               ENDDO UNPACK_RECV_FACE

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
 
DO IW = 1, M%NEWC

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
 
DO IW = 1, M%NEWC

   I = M%IJKW(6,IW)
   J = M%IJKW(7,IW)
   K = M%IJKW(8,IW)

   IOR0 = M%IJKW(4,IW)

   SELECT CASE(IOR0)
      CASE(1)
         IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
            SL%F(I,J,K) = SL%F(I,J,K) - 2.0_EB * SL%DXI2 * M%BXS(1,K)            ! Dirichlet
            IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,1000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
         ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
            SL%F(I,J,K) = SL%F(I,J,K) + SL%DXI * M%BXS(1,K)                       ! Neumann
            IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,2000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
         ENDIF
      CASE(-1)
         IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
            SL%F(I,J,K) = SL%F(I,J,K) - 2.0_EB * SL%DXI2 *M%BXF(1,K)            ! Dirichlet
            IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,1000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
         ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
            SL%F(I,J,K) = SL%F(I,J,K) - SL%DXI *M%BXF(1,K)                      ! Neumann
            IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,2000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
         ENDIF
      CASE(3)
         IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
            SL%F(I,J,K) = SL%F(I,J,K) - 2.0_EB * SL%DZI2 * M%BZS(I,1)            ! Dirichlet
            IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,1000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
         ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
            SL%F(I,J,K) = SL%F(I,J,K) + SL%DZI * M%BZS(I,1)                       ! Neumann
            IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,2000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
         ENDIF
      CASE(-3)
         IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
            SL%F(I,J,K) = SL%F(I,J,K) - 2.0_EB * SL%DZI2 * M%BZF(I,1)            ! Dirichlet
            IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,1000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
         ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
            SL%F(I,J,K) = SL%F(I,J,K) - SL%DZI  * M%BZF(I,1)                     ! Neumann
            IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,2000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
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
 
DO IW = 1, M%NEWC

   I = M%IJKW(6,IW)
   J = M%IJKW(7,IW)
   K = M%IJKW(8,IW)

   IOR0 = M%IJKW(4,IW)

   SELECT CASE(IOR0)
      CASE(1)
         IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
            SL%F(I,J,K) = SL%F(I,J,K) - 2.0_EB * SL%DXI2 * M%BXS(1,K)            ! Dirichlet
            IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,1000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
         ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
            SL%F(I,J,K) = SL%F(I,J,K) + SL%DXI * M%BXS(1,K)                       ! Neumann
            IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,2000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
         ENDIF
      CASE(-1)
         IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
            SL%F(I,J,K) = SL%F(I,J,K) - 2.0_EB * SL%DXI2 *M%BXF(1,K)            ! Dirichlet
            IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,1000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
         ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
            SL%F(I,J,K) = SL%F(I,J,K) - SL%DXI *M%BXF(1,K)                      ! Neumann
            IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,2000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
         ENDIF
      CASE(3)
         IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
            SL%F(I,J,K) = SL%F(I,J,K) - 2.0_EB * SL%DZI2 * M%BZS(I,1)            ! Dirichlet
            IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,1000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
         ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
            SL%F(I,J,K) = SL%F(I,J,K) + SL%DZI * M%BZS(I,1)                       ! Neumann
            IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,2000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
         ENDIF
      CASE(-3)
         IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
            SL%F(I,J,K) = SL%F(I,J,K) - 2.0_EB * SL%DZI2 * M%BZF(I,1)            ! Dirichlet
            IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,1000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
         ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
            SL%F(I,J,K) = SL%F(I,J,K) - SL%DZI  * M%BZF(I,1)                     ! Neumann
            IF (SCARC_DEBUG>=6) WRITE(SCARC_LU,2000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
         ENDIF
   END SELECT


ENDDO
 
 
TUSED_SCARC(27,NM)=TUSED_SCARC(27,NM)+SECOND()-TNOW_SETBDRY3D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_SETBDRY3D
1000 FORMAT(i3,': IOR=',i3,': Dirichlet :F(',i2,',',i2,',',i2,')=',e25.16,i5,i5)
2000 FORMAT(i3,': IOR=',i3,': Neumann   :F(',i2,',',i2,',',i2,')=',e25.16,i5,i5)
END SUBROUTINE SCARC_SETBDRY3D
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!  update interior boundary values along vertical directions
!!!  to have consistent values all over the domain
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_UPDATE(Z, NM, CNAME)
REAL (EB), POINTER, DIMENSION (:, :, :) :: Z
INTEGER :: NM, ILEVEL, I, J, K, ILMAX
CHARACTER (6) :: CNAME
REAL(EB):: TNOW_UPDATE_QUANTITY
 
TNOW_UPDATE_QUANTITY = SECOND()

ILMAX=SCARC(NM)%NLMAX
SL => SCARC(NM)%SLEVEL(ILMAX)
 
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
      NREQ_FACE = 0
      CALL SCARC_RECEIVE  (NCOM_FULL, ILMAX)   
      CALL SCARC_EXCHANGE (NCOM_FULL, ILMAX, NMV_NONE, NMV_NONE)
      DO K = 0, SL%KBAR + 1
         DO I = 0, SL%IBAR + 1
            Z (I, 1, K) = SL%Z(I, 1, K)
         ENDDO
      ENDDO
      !IF (NMESHES==2.OR.NMESHES==3) THEN
      !   IF (NM==1) Z(SL%IBAR+1,1,0)=Z(SL%IBAR,1,0)
      !   IF (NM==2) Z(0,1,0)=Z(1,1,0)
      !ELSE IF (NMESHES==4) THEN
      !   IF (NM==1) Z(SL%IBAR+1,1,0)=Z(SL%IBAR,1,0)
      !   IF (NM==3) Z(0,1,0)=Z(1,1,0)
      !ENDIF
 
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
 
      CALL SCARC_RECEIVE  (NCOM_FULL, ILEVEL) 
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
CALL SCARC_SHOW_LEVEL (Z, 'UPDT', CNAME, ILEVEL)
 
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
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! 2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   IF (TWO_D) THEN
 
      DO K = 0, SL%KBAR + 1
         DO I = 0, SL%IBAR + 1
            SL%Z (I, 1, K) = Z (I, 1, K)
         ENDDO
      ENDDO
      NREQ_FACE = 0
      CALL SCARC_RECEIVE  (NCOM_FULL, ILEVEL)   
      CALL SCARC_EXCHANGE (NCOM_FULL, ILEVEL, NMV_NONE, NMV_NONE)
      DO K = 0, SL%KBAR + 1
         DO I = 0, SL%IBAR + 1
            Z (I, 1, K) = SL%Z(I, 1, K)
         ENDDO
      ENDDO
      !IF (NMESHES==2.OR.NMESHES==3) THEN
      !   IF (NM==1) Z(SL%IBAR+1,1,0)=Z(SL%IBAR,1,0)
      !   IF (NM==2) Z(0,1,0)=Z(1,1,0)
      !ELSE IF (NMESHES==4) THEN
      !   IF (NM==1) Z(SL%IBAR+1,1,0)=Z(SL%IBAR,1,0)
      !   IF (NM==3) Z(0,1,0)=Z(1,1,0)
      !ENDIF

 
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
 
      CALL SCARC_RECEIVE (NCOM_FULL, ILEVEL)  !Aufruf 2,0
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
INTEGER :: II, jj, KK
CHARACTER (4) :: CROUTINE
CHARACTER (6) :: CNAME
 
IF (SCARC_DEBUG >= 3) THEN
   WRITE(SCARC_LU,*) '============================================================='
   WRITE(SCARC_LU,*) '===   COMPARE= ',CNAME,': ROUTINE= ',CROUTINE
   WRITE(SCARC_LU,*) '============================================================='
   IF (TWO_D) THEN
      DO KK = KBAR+1, 0, - 1
         WRITE(SCARC_LU, '(a,i3,a,10e12.3)') 'k=', KK, ' : ', (X(II, 1, KK), II=0,MIN(8,IBAR)+1)
      ENDDO
   ELSE
      DO KK = KBAR+1, 0, - 1
         WRITE(SCARC_LU, '(10e12.3)') ((X(II, jj, KK), II=0, MIN(8,IBAR)+1), jj=MIN(8,JBAR)+1, 1,-1)
         WRITE(SCARC_LU,*) '----------------------------------------'
      ENDDO
   ENDIF
ENDIF
END SUBROUTINE SCARC_SHOW

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! only for debugging reasons ... print also ghost values
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SHOW_LINE (X, CROUTINE, CNAME)
REAL (EB), POINTER, DIMENSION (:, :, :) :: X
INTEGER :: II, JJ, KK
CHARACTER (4) :: CROUTINE
CHARACTER (6) :: CNAME
 
IF (SCARC_DEBUG >= 3) THEN
   WRITE(SCARC_LU,*) '============================================================='
   WRITE(SCARC_LU,*) '===   COMPARE= ',CNAME,': ROUTINE= ',CROUTINE
   WRITE(SCARC_LU,*) '============================================================='
   IF (TWO_D) THEN
      DO KK = 5, 0, - 1
         WRITE(SCARC_LU, '(a,i3,a,3f26.16)') 'k=', KK, ' : ', (X(II, 1, KK), II=3,5)
      ENDDO
   ELSE
      DO KK = 5, 0, - 1
         WRITE(SCARC_LU, '(3f26.16)') ((X(II, JJ, KK), II=3,5), JJ=5,3,-1)
         WRITE(SCARC_LU,*) '----------------------------------------'
      ENDDO
   ENDIF
ENDIF

END SUBROUTINE SCARC_SHOW_LINE
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! only for debugging reasons ... related to single level
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SHOW_LEVEL (X, CROUTINE, CNAME, ILEVEL)
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
         WRITE(SCARC_LU, '(a,i3,a,10e12.3)') 'k=', KK, ' : ', (X(II, 1, KK), II=0,MIN(8,SL%IBAR)+1)
      ENDDO
   ELSE
      DO KK = SL%KBAR+1, 0, - 1
         WRITE(SCARC_LU, '(10e12.3)') ((X(II, JJ, KK), II=0, MIN(8,SL%IBAR)+1), JJ=MIN(8,SL%JBAR)+1, 0,-1)
         WRITE(SCARC_LU,*) '----------------------------------------'
      ENDDO
   ENDIF
ENDIF
 
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
!!! only for debugging reasons ...
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SHOW2(X, CROUTINE, CNAME)
REAL (EB), POINTER, DIMENSION (:, :, :) :: X
INTEGER :: II, jj, KK
CHARACTER (4) :: CROUTINE
CHARACTER (6) :: CNAME
 
IF (SCARC_DEBUG >= 3) THEN
   WRITE(SCARC_LU,*) '============================================================='
   WRITE(SCARC_LU,*) '===   COMPARE= ',CNAME,': ROUTINE= ',CROUTINE
   WRITE(SCARC_LU,*) '============================================================='
   IF (TWO_D) THEN
      WRITE(SCARC_LU,*) 
      WRITE(SCARC_LU,*) '===   COMPARE= ',CNAME,': SCARC_COUNT= ',SCARC_COUNT,'J=0'
      IF (NMESHES == 1) THEN
         DO KK = KBAR+1, 0, - 1
            WRITE(SCARC_LU, '(a,i3,a,10e12.3)') 'k=', KK, ' : ', (X(II, 0, KK), II=0,IBAR+1)
         ENDDO
      ELSE
         IF (IBAR==8) THEN
         DO KK = KBAR+1, 0, - 1
            WRITE(SCARC_LU, '(a,i3,a,10e12.3)') 'k=', KK, ' : ', (X(II, 0, KK), II=0, IBAR+1)
         ENDDO
         ELSE
         DO KK = KBAR+1, 0, - 1
            WRITE(SCARC_LU, '(a,i3,a,6e12.3)') 'k=', KK, ' : ', (X(II, 0, KK), II=0, IBAR+1)
         ENDDO
         ENDIF
      ENDIF
      WRITE(SCARC_LU,*) 
      WRITE(SCARC_LU,*) '===   COMPARE= ',CNAME,': SCARC_COUNT= ',SCARC_COUNT,'J=1'
      IF (NMESHES == 1) THEN
         DO KK = KBAR+1, 0, - 1
            WRITE(SCARC_LU, '(a,i3,a,10e12.3)') 'k=', KK, ' : ', (X(II, 1, KK), II=0,IBAR+1)
         ENDDO
      ELSE
         IF (IBAR==8) THEN
         DO KK = KBAR+1, 0, - 1
            WRITE(SCARC_LU, '(a,i3,a,10e12.3)') 'k=', KK, ' : ', (X(II, 1, KK), II=0, IBAR+1)
         ENDDO
         ELSE
         DO KK = KBAR+1, 0, - 1
            WRITE(SCARC_LU, '(a,i3,a,6e12.3)') 'k=', KK, ' : ', (X(II, 1, KK), II=0, IBAR+1)
         ENDDO
         ENDIF
      ENDIF
      WRITE(SCARC_LU,*) 
      WRITE(SCARC_LU,*) '===   COMPARE= ',CNAME,': SCARC_COUNT= ',SCARC_COUNT,'J=2'
      IF (NMESHES == 1) THEN
         DO KK = KBAR+1, 0, - 1
            WRITE(SCARC_LU, '(a,i3,a,10e12.3)') 'k=', KK, ' : ', (X(II, 2, KK), II=0,IBAR+1)
         ENDDO
      ELSE
         IF (IBAR==8) THEN
         DO KK = KBAR+1, 0, - 1
            WRITE(SCARC_LU, '(a,i3,a,10e12.3)') 'k=', KK, ' : ', (X(II, 2, KK), II=0, IBAR+1)
         ENDDO
         ELSE
         DO KK = KBAR+1, 0, - 1
            WRITE(SCARC_LU, '(a,i3,a,6e12.3)') 'k=', KK, ' : ', (X(II, 2, KK), II=0, IBAR+1)
         ENDDO
         ENDIF
      ENDIF
            WRITE(SCARC_LU, *) 
   ELSE
      IF (NMESHES == 1) THEN
         DO KK = KBAR+1, 0, - 1
            WRITE(SCARC_LU, '(10e12.3)') ((X(II, jj, KK), II=0, IBAR+1), jj=JBAR, 1,-1)
            WRITE(SCARC_LU,*) '----------------------------------------'
         ENDDO
      ELSE
         DO KK = KBAR+1, 0, - 1
            WRITE(SCARC_LU, '(6e12.3)') ((X(II, jj, KK), II=0, IBAR+1), jj=JBAR, 1,-1)
            WRITE(SCARC_LU,*) '----------------------------------------'
         ENDDO
      ENDIF
   ENDIF
ENDIF
 
END SUBROUTINE SCARC_SHOW2

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
 
