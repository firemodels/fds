MODULE SCARC_SOLVER
 
USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_POINTERS
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
USE POIS, ONLY: H2CZIS, H3CZSS, H2CZSS, H2CYSS, H3CSSS
 
IMPLICIT NONE
 
!INTEGER :: MPI_COMM_WORLD=0, MPI_DOUBLE_PRECISION=0, MPI_INTEGER=0, MPI_SUM=0
!INTEGER :: MPI_REQUEST_NULL=0, MPI_STATUS_SIZE=0
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Global variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CHARACTER (40) :: SCARC_MSG='msg/   _scarc'              ! name of file for ScaRC debug messages
INTEGER        :: SCARC_METHOD=0                         ! method for pressure solver (0:FFT/1:CG/2:MG)
INTEGER        :: SCARC_DEBUG=0                          ! debug level (0: no debug messages)
INTEGER        :: SCARC_COMPARE=0                        ! print out vectors for comparison in separate program
INTEGER        :: SCARC_PERIODIC=0                       ! periodic boundary conditions (0:no/1:yes)
INTEGER        :: SCARC_CASE=0                           ! choose different initial solution
INTEGER        :: SCARC_COUNT=2                          ! counter for comparison vectors 
INTEGER        :: SCARC_COUNT2                           ! counter2 for comparison vectors
INTEGER        :: SCARC_PRECON=4            
INTEGER        :: SCARC_PRECON_CG=4                      ! preconditioner for cg-method (default: SSOR)
INTEGER        :: SCARC_SMOOTH=3                         ! smoother for mg-method (default: GSTRIX)
INTEGER        :: SCARC_NIT=1000
INTEGER        :: SCARC_NIT_CG=1000                      ! max number of iterations for cg-method
INTEGER        :: SCARC_NIT_MG=1000                      ! max number of iterations for mg-method
INTEGER        :: SCARC_NIT_SM=1000                      ! max number of iterations for smoothing-method
REAL (EB)      :: SCARC_EPS=1.E-6_EB
REAL (EB)      :: SCARC_EPS_REL=1.E-2_EB
REAL (EB)      :: SCARC_EPS_STOP=1.E+6_EB
REAL (EB)      :: SCARC_EPS_CG=1.E-6_EB                  ! convergence accuracy for cg-method
REAL (EB)      :: SCARC_EPS_MG=1.E-6_EB                  ! convergence accuracy for mg-method
REAL (EB)      :: SCARC_EPS_SM=1.E-6_EB                  ! convergence accuracy for smoothing-method
REAL (EB)      :: SCARC_EPS_REL_CG=1.E-2_EB              ! relative accuracy for cg-method
REAL (EB)      :: SCARC_EPS_REL_MG=1.E-2_EB              ! relative accuracy for mg-method
REAL (EB)      :: SCARC_EPS_REL_SM=1.E-2_EB              ! relative accuracy for smoothing-method
REAL (EB)      :: SCARC_EPS_STOP_CG=1.E+6_EB             ! divergence accuracy for cg-method
REAL (EB)      :: SCARC_EPS_STOP_MG=1.E+6_EB             ! divergence accuracy for mg-method
REAL (EB)      :: SCARC_EPS_STOP_SM=1.E+6_EB             ! divergence accuracy for mg-method


LOGICAL :: BEDGE, BDIAG     ! diagonal communication needed ?

INTEGER :: NumMaster, NumSlaves
INTEGER :: myMasterId, mySlaveId, myself
INTEGER :: MpiGroupWorld, MpiGroupSlaves
INTEGER :: myMPI_COMM_WORLD, myMpi_COMM_SLAVES
INTEGER :: NumExclude, RankExclude (4)
 
INTEGER :: NTYPE_NONE = 0,  &
           NTYPE_UPDATE = -1

INTEGER :: NCOM_TYPE_GLOBAL= 1, &
           NCOM_TYPE_LOCAL = 2

INTEGER :: NCOM_INIT = 0, &
           NCOM_MATV = 1, &
           NCOM_FACE = 1, &
           NCOM_EDGE = 2, &
           NCOM_DIAG = 3, &
           NCOM_FULL = 4

INTEGER :: NUPDATE_FULL    = 1, &
           NUPDATE_MEDIUM  = 2, &
           NUPDATE_VEL     = 3, &
           NUPDATE_FV      = 4, &
           NUPDATE_Q       = 5, &
           NUPDATE_H       = 6, &
           NUPDATE_HS      = 7

INTEGER :: NSMOOTH_NONE    = 0, &
           NSMOOTH_JACOBI  = 1, &
           NSMOOTH_GS      = 2, &
           NSMOOTH_SOR     = 3, &
           NSMOOTH_SSOR    = 4, &
           NSMOOTH_ILU     = 5, &
           NSMOOTH_TRI1    = 6, &
           NSMOOTH_TRI2    = 8, &
           NSMOOTH_ADI     = 9, &
           NSMOOTH_GSTRI1  =10, &
           NSMOOTH_GSTRI2  =11, &
           NSMOOTH_GSADI   =12


INTEGER :: NEXCHANGE_EXTENDED = 1
 
INTEGER :: NNLEVEL = 10, NLEVEL1=1

INTEGER :: SNODE, RNODE, NDIAG
INTEGER :: NFACE, NEDGE, NVRTX
 
! protocol information for data exchange
INTEGER :: NREQ_FACE, NREQ_EDGE, NREQ_DIAG
!INTEGER :: TAG_FACE  , TAG_EDGE  , TAG_DIAG
INTEGER, ALLOCATABLE, DIMENSION (:)    :: REQ_FACE,  REQ_EDGE,  REQ_DIAG
INTEGER, ALLOCATABLE, DIMENSION (:, :) :: TAGS_FACE, TAGS_EDGE, TAGS_DIAG
INTEGER, ALLOCATABLE, DIMENSION (:, :) :: STAT_FACE, STAT_EDGE, STAT_DIAG
 
! Alte Strukturen, nochmal checken !!!
!INTEGER :: IJKW2_FACE (6, 12), IJKW2_EDGE (12, 12), IJKW2_DIAG (8, 12)
INTEGER, ALLOCATABLE, DIMENSION (:, :)    :: NBR_FACE, NBR_EDGE, NBR_DIAG

INTEGER, PARAMETER :: N_TIMERS_SCARC=35
REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: TUSED_SCARC
REAL(EB):: SCARC_TNOW


 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Scarc type for multigrid level on 'own' mesh
! in case of the cg-method, only the finest level is used
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
TYPE SCARC_LEVEL_TYPE
 
   ! local dimensions  on corresponding level
   INTEGER :: IBAR, JBAR, KBAR

   ! number of local and global grid cells on corresponding level
   INTEGER :: NCELLS_LOCAL, NCELLS_GLOBAL
   INTEGER :: NW_FACE, NW_EDGE, NW_DIAG

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

   ! neighbourship structures along faces, edges and vertices
   INTEGER, POINTER, DIMENSION (:, :) :: NIC_FACE,  NIC_EDGE,  NIC_DIAG
   INTEGER, POINTER, DIMENSION (:, :) :: IJKW_FACE, IJKW_EDGE, IJKW_DIAG
   INTEGER, POINTER, DIMENSION (:, :) :: I_MIN_FACE, I_MAX_FACE
   INTEGER, POINTER, DIMENSION (:, :) :: J_MIN_FACE, J_MAX_FACE
   INTEGER, POINTER, DIMENSION (:, :) :: K_MIN_FACE, K_MAX_FACE
   INTEGER, POINTER, DIMENSION (:, :) :: I_MIN_EDGE, I_MAX_EDGE
   INTEGER, POINTER, DIMENSION (:, :) :: J_MIN_EDGE, J_MAX_EDGE
   INTEGER, POINTER, DIMENSION (:, :) :: K_MIN_EDGE, K_MAX_EDGE
   INTEGER, POINTER, DIMENSION (:, :) :: I_MIN_DIAG, I_MAX_DIAG
   INTEGER, POINTER, DIMENSION (:, :) :: J_MIN_DIAG, J_MAX_DIAG
   INTEGER, POINTER, DIMENSION (:, :) :: K_MIN_DIAG, K_MAX_DIAG

   ! matrices and iteration vectors for global defect correction on corresponding level
   REAL (EB), POINTER, DIMENSION (:)       :: DD, DL, DU, LD, DAUX, PERIOD
   REAL (EB), POINTER, DIMENSION (:, :)    :: AG, AL
   REAL (EB), POINTER, DIMENSION (:, :, :) :: X, Y, F, D, G, R, Z, TMP
   REAL (EB), POINTER, DIMENSION (:, :)    :: BXS0, BXF0, BYS0, BYF0, BZS0, BZF0
   REAL (EB) :: ASUBX, ASUBY, ASUBZ, ASUB
 
   ! communication vectors
   REAL (EB), POINTER, DIMENSION (:, :, :) :: Y_MATV, R_MATV, G_MATV
   REAL (EB), POINTER, DIMENSION (:, :, :) :: Y_FACE, R_FACE, G_FACE
   REAL (EB), POINTER, DIMENSION (:, :, :) :: Y_EDGE, R_EDGE, G_EDGE
   REAL (EB), POINTER, DIMENSION (:, :, :) :: Y_DIAG, R_DIAG, G_DIAG
   REAL (EB), POINTER, DIMENSION (:, :, :) :: Z_FACE, Z_EDGE, Z_DIAG
 
   ! iteration parameters for solver method on corresponding level
   REAL (EB) :: RES, RESIN, OMEGA, CAPPA
   REAL (EB) :: OMEGA_MG, OMEGA_SMOOTH=0.8E+0_EB
   REAL (EB) :: RES_MG, RES_SMOOTH
   REAL (EB) :: RESIN_MG, RESIN_SMOOTH
   REAL (EB) :: CAPPA_MG, CAPPA_SMOOTH
   REAL (EB) :: EPS, EPS_REL
   INTEGER   :: ITE, NIT
   INTEGER   :: ITE_MG, NIT_MG
   INTEGER   :: ITE_SMOOTH, NIT_SMOOTH

END TYPE SCARC_LEVEL_TYPE
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Scarc type for multigrid level on 'other' mesh
! in case of the cg-method, only the finest level is used
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
TYPE OSCARC_LEVEL_TYPE
 
   ! local dimensions  on corresponding level
   INTEGER :: IBAR, JBAR, KBAR

   ! number of local and global grid cells on corresponding level
   INTEGER :: NCELLS_LOCAL, NCELLS_GLOBAL
   INTEGER :: NW_FACE, NW_EDGE, NW_DIAG

   ! neighbourship structures along faces, edges and vertices
   INTEGER, POINTER, DIMENSION (:, :) :: NIC_FACE,  NIC_EDGE,  NIC_DIAG
   INTEGER, POINTER, DIMENSION (:, :) :: IJKW_FACE, IJKW_EDGE, IJKW_DIAG
   INTEGER, POINTER, DIMENSION (:, :) :: I_MIN_FACE, I_MAX_FACE
   INTEGER, POINTER, DIMENSION (:, :) :: J_MIN_FACE, J_MAX_FACE
   INTEGER, POINTER, DIMENSION (:, :) :: K_MIN_FACE, K_MAX_FACE
   INTEGER, POINTER, DIMENSION (:, :) :: I_MIN_EDGE, I_MAX_EDGE
   INTEGER, POINTER, DIMENSION (:, :) :: J_MIN_EDGE, J_MAX_EDGE
   INTEGER, POINTER, DIMENSION (:, :) :: K_MIN_EDGE, K_MAX_EDGE
   INTEGER, POINTER, DIMENSION (:, :) :: I_MIN_DIAG, I_MAX_DIAG
   INTEGER, POINTER, DIMENSION (:, :) :: J_MIN_DIAG, J_MAX_DIAG
   INTEGER, POINTER, DIMENSION (:, :) :: K_MIN_DIAG, K_MAX_DIAG

   ! iteration vectors for global defect correction on corresponding level
   REAL (EB), POINTER, DIMENSION (:, :, :) :: X, Y, F, D, G, R, Z, TMP

   ! communication vectors
   REAL (EB), POINTER, DIMENSION (:, :, :) :: Y_FACE, R_FACE, G_FACE
   REAL (EB), POINTER, DIMENSION (:, :, :) :: Y_EDGE, R_EDGE, G_EDGE
   REAL (EB), POINTER, DIMENSION (:, :, :) :: Y_DIAG, R_DIAG, G_DIAG
   REAL (EB), POINTER, DIMENSION (:, :, :) :: Z_FACE, Z_EDGE, Z_DIAG
 
END TYPE OSCARC_LEVEL_TYPE
 

 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Scarc type on 'own' mesh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
TYPE SCARC_TYPE
 
   ! global MG-levels and MG-cycle variable
   INTEGER :: NLEVEL, NLEVEL_MIN, NLEVEL_MAX, NLEVEL_DIFF
   INTEGER :: ICYCLE, ITE

   ! communication vectors (used for all levels)
   REAL (EB), POINTER, DIMENSION (:) :: SEND_FACE, RECV_FACE
   REAL (EB), POINTER, DIMENSION (:) :: SEND_EDGE, RECV_EDGE
   REAL (EB), POINTER, DIMENSION (:) :: SEND_DIAG, RECV_DIAG

   ! local and global scalar products (used for all levels)
   REAL (EB) :: SP_LOCAL, SP_GLOBAL
   REAL (EB), POINTER, DIMENSION (:) :: SP_LOCAL0, SP_GLOBAL0
   REAL (EB) :: EPS_CONV, EPS_DIVG, EPS_REL, EPS, RES, RESIN, CAPPA
 

   ! neighborship description arrays
   INTEGER, POINTER, DIMENSION (:, : , : ) :: NDIAG_NBR
   INTEGER, POINTER, DIMENSION (:)         :: MIBAR, MJBAR, MKBAR
   INTEGER, POINTER, DIMENSION (:, :)      :: MLEVEL0
   INTEGER, POINTER, DIMENSION (:)         :: MLEVEL
   INTEGER, POINTER, DIMENSION (:, : )     :: KCYCLE
   INTEGER, POINTER, DIMENSION (:)         :: KSMOOTH_PRE, KSMOOTH_POST

   INTEGER, POINTER, DIMENSION (:, :) :: IOR_FACE

   ! neighboring ScaRC-structures
   TYPE (OSCARC_TYPE), POINTER, DIMENSION (:) :: OSCARC

   ! different mesh levels (in case of global multigrid)
   INTEGER :: NFACE0, NEDGE0, NVRTX0
   INTEGER :: MIBAR_MIN, MJBAR_MIN, MKBAR_MIN, MESH_MIN

   TYPE (SCARC_LEVEL_TYPE), POINTER, DIMENSION (:) :: SLEVEL
 
   ! some iteration parameters
   LOGICAL :: BPRECON, BGSADI, BREL_ERR
 
END TYPE SCARC_TYPE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Scarc type on 'other' mesh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
TYPE OSCARC_TYPE
 
   REAL (EB), POINTER, DIMENSION (:) :: SEND_FACE, RECV_FACE
   REAL (EB), POINTER, DIMENSION (:) :: SEND_EDGE, RECV_EDGE
   REAL (EB), POINTER, DIMENSION (:) :: SEND_DIAG, RECV_DIAG

   TYPE (OSCARC_LEVEL_TYPE), POINTER, DIMENSION (:) :: SLEVEL
 
END TYPE OSCARC_TYPE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Define necessary types
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
TYPE (SCARC_TYPE),  SAVE, DIMENSION(:), ALLOCATABLE, TARGET :: SCARC
TYPE (OSCARC_TYPE), SAVE, DIMENSION(:), ALLOCATABLE, TARGET :: OSCARC
 
TYPE (SCARC_TYPE),        POINTER ::   S
TYPE (OSCARC_TYPE),       POINTER ::  OS

TYPE (SCARC_LEVEL_TYPE),  POINTER ::  SL,  SLMAX,  SLMIN,  SLHI,  SLLO
TYPE (OSCARC_LEVEL_TYPE), POINTER :: OSL, OSLMAX, OSLMIN, OSLHI, OSLLO


!INCLUDE 'mpif.h'
 
CONTAINS
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! SCARC_INITIALIZE : Initialize Scarc structures for 2D- or 3D-case
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_INITIALIZE (NM)

INTEGER :: IM, NM, IZERO, NFINE
INTEGER :: IMIN, JMIN, KMIN
INTEGER :: IBAR0, JBAR0, KBAR0
REAL(EB):: TNOW_INITIALIZE

TNOW_INITIALIZE = SECOND()

USE_MPI=.FALSE.

BEDGE=.FALSE.
BDIAG=.FALSE.

!WRITE(9,*) 'HALLO, AN DAS RICHTIGE SETZEN VON BEDGE IN 3D DENKEN !!'
!!!
!!! Allocate SCARC and OSCARC structures
!!!
ALLOCATE (SCARC(NMESHES), STAT=IZERO)
CALL CHKMEMERR ('SCARC_INITIALIZE', 'SCARC', IZERO)
 
S => SCARC (NM)

ALLOCATE (S%OSCARC(NMESHES), STAT=IZERO)
CALL CHKMEMERR ('SCARC_INITIALIZE', 'OSCARC', IZERO)


!!!
!!! Open debug-file if requested
!!!
IF (SCARC_DEBUG >= 1) THEN
   WRITE (SCARC_MSG(5:7), '(i2.2)') NM
   OPEN (9, FILE=SCARC_MSG)
ENDIF

IF (SCARC_METHOD==0) GOTO 12345

!!!
!!! Initialize time measurement array
!!!
ALLOCATE(TUSED_SCARC(N_TIMERS_SCARC,NMESHES),STAT=IZERO)
CALL ChkMemErr('SCARC_INITIALIZE','TUSED_SCARC',IZERO)

TUSED_SCARC      = 0._EB
TUSED_SCARC(0,:) = SECOND()



!!!
!!! store mesh dimensions on each direction of each mesh within SCARC-structure;
!!! get minimum step size over all meshes for each direction
!!!
ALLOCATE (S%MIBAR(NMESHES), STAT=IZERO)
CALL CHKMEMERR ('SCARC_INITIALIZE', 'S%MIBAR', IZERO)
 
ALLOCATE (S%MJBAR(NMESHES), STAT=IZERO)
CALL CHKMEMERR ('SCARC_INITIALIZE', 'S%MJBAR', IZERO)
 
ALLOCATE (S%MKBAR(NMESHES), STAT=IZERO)
CALL CHKMEMERR ('SCARC_INITIALIZE', 'S%MKBAR', IZERO)

 
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

   IF (TWO_D) THEN
      IF (S%MJBAR(IM)<=S%MJBAR_MIN) THEN
         S%MJBAR_MIN=S%MJBAR(IM)
         JMIN=IM
      ENDIF
   ENDIF

   IF (S%MKBAR(IM)<=S%MKBAR_MIN) THEN
      S%MKBAR_MIN=S%MKBAR(IM)
      KMIN=IM
   ENDIF

   IF (SCARC_DEBUG>=2) THEN
      WRITE(9,*) '------------------ Mesh ',IM,'--------'
      WRITE(9,*) 'S%MIBAR(',IM,')=',S%MIBAR(IM)
      IF (TWO_D) WRITE(9,*) 'S%MJBAR(',IM,')=',S%MJBAR(IM)
      WRITE(9,*) 'S%MKBAR(',IM,')=',S%MKBAR(IM)
   ENDIF

ENDDO


IF (SCARC_DEBUG>=2) THEN
   WRITE(9,*) 'S%MIBAR_MIN=',S%MIBAR_MIN
   IF (TWO_D) WRITE(9,*) 'S%MJBAR_MIN=',S%MJBAR_MIN
   WRITE(9,*) 'S%MKBAR_MIN=',S%MKBAR_MIN
ENDIF


!!!
!!! Determine number of grid levels  (1 for CG-method, NLEVEL for MG-method)
!!! and corresponding numbers of cells in each direction
!!!
ALLOCATE (S%MLEVEL0(NMESHES,3), STAT=IZERO)
CALL CHKMEMERR ('SCARC_INITIALIZE', 'S%MLEVEL0', IZERO)
S%MLEVEL0=0
 
ALLOCATE (S%MLEVEL(NMESHES), STAT=IZERO)
CALL CHKMEMERR ('SCARC_INITIALIZE', 'S%MLEVEL', IZERO)
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
      ENDIF 
      !WRITE(9,*) 'NFINE=',NFINE
   ENDDO
   S%MLEVEL0(IM,1) = NFINE
   S%MLEVEL(IM)    = S%MLEVEL0(IM,1)

   IF (S%MLEVEL(IM)<S%NLEVEL) S%NLEVEL = NFINE

 !WRITE(9,*) '1:S%MLEVEL(',IM,')=',S%MLEVEL(IM)
 !WRITE(9,*) '1:S%NLEVEL     =',S%NLEVEL
   
   ! number of refinements in y-direction
   IF (.NOT.TWO_D) THEN

      DO NFINE=1,50
         IF (MOD(JBAR0,2)/=0) THEN
            WRITE(*,*) 'JBAR=',JBAR0,' NOT YET ALLOWED FOR SCARC-MULTIGRID !'
            EXIT 
         ELSE
            JBAR0=JBAR0/2
            IF (JBAR0==1) EXIT
         ENDIF 
      ENDDO
      S%MLEVEL0(IM,2)=NFINE
      IF (S%MLEVEL0(IM,2)<S%MLEVEL(IM)) S%MLEVEL(IM)=S%MLEVEL0(IM,2)

      IF (S%MLEVEL(IM)<S%NLEVEL) S%NLEVEL = NFINE

      IF (SCARC_DEBUG >= 2) THEN
         WRITE(9,*) '2:S%MLEVEL(',IM,')=',S%MLEVEL(IM)
         WRITE(9,*) '2:S%NLEVEL     =',S%NLEVEL
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
      ENDIF 
   ENDDO
   S%MLEVEL0(IM,3)=NFINE
   IF (S%MLEVEL0(IM,3)<S%MLEVEL(IM)) S%MLEVEL(IM)=S%MLEVEL0(IM,3)

   IF (S%MLEVEL(IM)<S%NLEVEL) S%NLEVEL =NFINE
   

   IF (SCARC_DEBUG >= 2) then
      WRITE (9,*) '============ MESH ', IM
      WRITE(9,*) '3:S%MLEVEL(',IM,')=',S%MLEVEL(IM)
      WRITE(9,*) '3:S%NLEVEL     =',S%NLEVEL
      WRITE (9,*) 
      WRITE (9,*) 'S%MLEVEL0(',IM,',1)  =',S%MLEVEL0(IM,1)
      WRITE (9,*) 'S%MLEVEL0(',IM,',2)  =',S%MLEVEL0(IM,2)
      WRITE (9,*) 'S%MLEVEL0(',IM,',3)  =',S%MLEVEL0(IM,3)
      WRITE (9,*) 
      WRITE (9,*) 'S%MLEVEL(',IM,')  =',S%MLEVEL(IM)
      WRITE (9,*) 'S%NLEVEL          =',S%NLEVEL   
   ENDIF

ENDDO


!!!
!!! Determine number of levels for global method
!!! depends on the minimum number of levels found on all meshes
!!!
S%NLEVEL_MAX  = S%MLEVEL(NM)
IF (SCARC_METHOD==1) THEN
   S%NLEVEL = 1
   S%NLEVEL_MIN  = S%NLEVEL_MAX
ELSE
   S%NLEVEL_MIN  = S%MLEVEL(NM)-S%NLEVEL+1
ENDIF
S%NLEVEL_DIFF = S%NLEVEL_MAX - S%NLEVEL_MIN

IF (SCARC_DEBUG >= 2) then
   WRITE (9,*) 'S%NLEVEL     =',S%NLEVEL
   WRITE (9,*) 'S%NLEVEL_MAX =',S%NLEVEL_MAX
   WRITE (9,*) 'S%NLEVEL_MIN =',S%NLEVEL_MIN
   WRITE (9,*) 'S%NLEVEL_DIFF=',S%NLEVEL_DIFF
ENDIF

!!!
!!! Allocate corresponding number of SCARC_LEVEL-structures
!!!
ALLOCATE (S%SLEVEL(S%NLEVEL_MIN:S%NLEVEL_MAX), STAT=IZERO)
CALL CHKMEMERR ('SCARC_NEIGHBORS2D', 'SLEVEL', IZERO)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Initialize neighborship structure, matrices, mesh exchange and global solver method
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (TWO_D) THEN

   CALL SCARC_INITIALIZE_MESHES2D (NM)
   CALL SCARC_INITIALIZE_NEIGHBORS2D(NM)
   CALL SCARC_INITIALIZE_MATRICES2D (NM)
   CALL SCARC_INITIALIZE_COMMUNICATION2D (NM)
   CALL SCARC_INITIALIZE_SOLVER2D (NM)

ELSE

   CALL SCARC_INITIALIZE_MESHES3D (NM)         
   write(*,*) '3D-ScaRC still under construction, do not use at the moment !!'
   stop
   !CALL SCARC_INITIALIZE_NEIGHBORS3D (NM)
   CALL SCARC_INITIALIZE_MATRICES3D (NM)
   !CALL SCARC_INITIALIZE_MESH_EXCHANGE3D (NM)
   CALL SCARC_INITIALIZE_SOLVER3D (NM)

ENDIF
 
12345 CONTINUE
IF (SCARC_DEBUG>=2) WRITE(9,*) 'Leaving scarc_initialize'

TUSED_SCARC(24,NM)=TUSED_SCARC(24,NM)+SECOND()-TNOW_INITIALIZE
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_INITIALIZE
END SUBROUTINE SCARC_INITIALIZE
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! SCARC_INITIALIZE_MESHES2D : initialize grid structures in 2D 
!!!  -  one fine grid on every MESH for CG-method
!!!  -  complete hierarchy of grids on every MESH for MG-method
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_INITIALIZE_MESHES2D (NM)
 
INTEGER :: NM, I, K, ILEVEL
INTEGER :: IBAR0, JBAR0, KBAR0, IZERO=0, IERR=0
TYPE (MESH_TYPE), POINTER :: M
REAL(EB):: TNOW_MESHES2D
INCLUDE 'mpif.h'
 

TNOW_MESHES2D = SECOND()

M => MESHES(NM)

IBAR0=M%IBAR
JBAR0=M%JBAR
KBAR0=M%KBAR


GRID_LEVEL_LOOP2D: DO ILEVEL=S%NLEVEL_MAX,S%NLEVEL_MIN,-1

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
   ALLOCATE (SL%XX(0:SL%IBAR), STAT=IZERO)
   CALL CHKMEMERR ('SCARC_INITIALIZE_MESHES2D', 'SL%XX', IZERO)
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

   ALLOCATE (SL%YY(0), STAT=IZERO)
   CALL CHKMEMERR ('SCARC_INITIALIZE_MESHES2D', 'SL%YY', IZERO)
   SL%YY(0) =0.0_EB
   SL%YY_MIN=0.0_EB
   SL%YY_MAX=0.0_EB
   SL%DY_MIN=0.0_EB
   SL%DY_MAX=0.0_EB

   ALLOCATE (SL%ZZ(0:SL%KBAR), STAT=IZERO)
   CALL CHKMEMERR ('SCARC_INITIALIZE_MESHES2D', 'SL%ZZ', IZERO)
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
      CALL MPI_ALLREDUCE(SL%NCELLS_LOCAL,SL%NCELLS_GLOBAL,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERR)
   ELSE
      SL%NCELLS_GLOBAL = SL%NCELLS_LOCAL
   ENDIF
    
    

   ! ----------------------------------------------------------------------------------
   ! print debug information if requested
   ! ----------------------------------------------------------------------------------
   IF (SCARC_DEBUG.ge.2) THEN
      WRITE(9,*) ' ======================= SCARC_INITIALIZE_MESHES2D =========='
      WRITE(9,*) '==========> LEVEL ', ILEVEL
      WRITE(9,*) 'SL(',ILEVEL,')%IBAR=',SL%IBAR
      WRITE(9,*) 'SL(',ILEVEL,')%JBAR=',SL%JBAR
      WRITE(9,*) 'SL(',ILEVEL,')%KBAR=',SL%KBAR
      WRITE(9,*) 
      WRITE(9,*) 'SL(',ILEVEL,')%DX=',SL%DX
      WRITE(9,*) 'SL(',ILEVEL,')%DY=',SL%DY
      WRITE(9,*) 'SL(',ILEVEL,')%DZ=',SL%DZ
      WRITE(9,*) 
      WRITE(9,*) 'M%XS=',M%XS
      WRITE(9,*) 'M%YS=',M%YS
      WRITE(9,*) 'M%ZS=',M%ZS
      WRITE(9,*) 
      WRITE(9,*) 'SL%XX_MIN=',SL%XX_MIN
      WRITE(9,*) 'SL%XX_MAX=',SL%XX_MAX
      WRITE(9,*) 'SL%DX_MIN=',SL%DX_MIN
      WRITE(9,*) 'SL%DX_MAX=',SL%DX_MAX
      WRITE(9,*) 
      WRITE(9,*) 'SL%YY_MIN=',SL%YY_MIN
      WRITE(9,*) 'SL%YY_MAX=',SL%YY_MAX
      WRITE(9,*) 'SL%DY_MIN=',SL%DY_MIN
      WRITE(9,*) 'SL%DY_MAY=',SL%DY_MAX
      WRITE(9,*) 
      WRITE(9,*) 'SL%ZZ_MIN=',SL%ZZ_MIN
      WRITE(9,*) 'SL%ZZ_MAX=',SL%ZZ_MAX
      WRITE(9,*) 'SL%DZ_MIN=',SL%DZ_MIN
      WRITE(9,*) 'SL%DZ_MAX=',SL%DZ_MAX
      WRITE(9,*) 
      WRITE(9,'(a,i3,a,16f7.3)') 'SL(',ILEVEL,')%XX=',(SL%XX(I),I=0,SL%IBAR)
      !WRITE(9,'(a,i3,a,16f7.3)') 'SL(',ILEVEL,')%YY=',(SL%YY(I),I=0,0)
      WRITE(9,'(a,i3,a,16f7.3)') 'SL(',ILEVEL,')%ZZ=',(SL%ZZ(I),I=0,SL%KBAR)
      WRITE(9,*) 
      WRITE(9,*) 'NMESHES=', NMESHES
      WRITE(9,*) 'NCELLS_LOCAL=', SL%NCELLS_LOCAL
      WRITE(9,*) 'NCELLS_GLOBAL=', SL%NCELLS_GLOBAL
      CALL flush (9)
   ENDIF

   IBAR0=IBAR0/2
   KBAR0=KBAR0/2

ENDDO GRID_LEVEL_LOOP2D
 
TUSED_SCARC(1,NM)=TUSED_SCARC(1,NM)+SECOND()-TNOW_MESHES2D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_MESHES2D

END SUBROUTINE SCARC_INITIALIZE_MESHES2D
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! SCARC_INITIALIZE_MESHES3D : initialize grid structures in 3D 
!!!  -  one fine grid on every MESH for CG-method
!!!  -  complete hierarchy of grids on every MESH for MG-method
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_INITIALIZE_MESHES3D (NM)
 
INTEGER :: NM, I, K, ILEVEL
INTEGER :: IBAR0, JBAR0, KBAR0, IZERO=0, IERR=0
TYPE (MESH_TYPE), POINTER :: M
REAL(EB):: TNOW_MESHES3D
INCLUDE 'mpif.h'

TNOW_MESHES3D = SECOND()
 
M => MESHES(NM)

IBAR0=M%IBAR
JBAR0=M%JBAR
KBAR0=M%KBAR


GRID_LEVEL_LOOP3D: DO ILEVEL=S%NLEVEL_MAX,S%NLEVEL_MIN,-1

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
   ALLOCATE (SL%XX(0:SL%IBAR), STAT=IZERO)
   CALL CHKMEMERR ('SCARC_INITIALIZE_MESHES3D', 'SL%XX', IZERO)
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

   ALLOCATE (SL%YY(0:SL%KBAR), STAT=IZERO)
   CALL CHKMEMERR ('SCARC_INITIALIZE_MESHES3D', 'SL%YY', IZERO)
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

   ALLOCATE (SL%ZZ(0:SL%KBAR), STAT=IZERO)
   CALL CHKMEMERR ('SCARC_INITIALIZE_MESHES3D', 'SL%ZZ', IZERO)
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
      CALL MPI_ALLREDUCE(SL%NCELLS_LOCAL,SL%NCELLS_GLOBAL,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERR)
   ELSE
      SL%NCELLS_GLOBAL = SL%NCELLS_LOCAL
   ENDIF
    
    

   ! ----------------------------------------------------------------------------------
   ! print debug information if requested
   ! ----------------------------------------------------------------------------------
   IF (SCARC_DEBUG.ge.2) THEN
      WRITE(9,*) ' ======================= SCARC_INITIALIZE_MESHES3D =========='
      WRITE(9,*) '==========> LEVEL ', ILEVEL
      WRITE(9,*) 'SL(',ILEVEL,')%IBAR=',SL%IBAR
      WRITE(9,*) 'SL(',ILEVEL,')%JBAR=',SL%JBAR
      WRITE(9,*) 'SL(',ILEVEL,')%KBAR=',SL%KBAR
      WRITE(9,*) 
      WRITE(9,*) 'SL(',ILEVEL,')%DX=',SL%DX
      WRITE(9,*) 'SL(',ILEVEL,')%DY=',SL%DY
      WRITE(9,*) 'SL(',ILEVEL,')%DZ=',SL%DZ
      WRITE(9,*) 
      WRITE(9,*) 'M%XS=',M%XS
      WRITE(9,*) 'M%YS=',M%YS
      WRITE(9,*) 'M%ZS=',M%ZS
      WRITE(9,*) 
      WRITE(9,*) 'SL%XX_MIN=',SL%XX_MIN
      WRITE(9,*) 'SL%XX_MAX=',SL%XX_MAX
      WRITE(9,*) 'SL%DX_MIN=',SL%DX_MIN
      WRITE(9,*) 'SL%DX_MAX=',SL%DX_MAX
      WRITE(9,*) 
      WRITE(9,*) 'SL%YY_MIN=',SL%YY_MIN
      WRITE(9,*) 'SL%YY_MAX=',SL%YY_MAX
      WRITE(9,*) 'SL%DY_MIN=',SL%DY_MIN
      WRITE(9,*) 'SL%DY_MAY=',SL%DY_MAX
      WRITE(9,*) 
      WRITE(9,*) 'SL%ZZ_MIN=',SL%ZZ_MIN
      WRITE(9,*) 'SL%ZZ_MAX=',SL%ZZ_MAX
      WRITE(9,*) 'SL%DZ_MIN=',SL%DZ_MIN
      WRITE(9,*) 'SL%DZ_MAX=',SL%DZ_MAX
      WRITE(9,*) 
      WRITE(9,'(a,i3,a,16f7.3)') 'SL(',ILEVEL,')%XX=',(SL%XX(I),I=0,SL%IBAR)
      !WRITE(9,'(a,i3,a,16f7.3)') 'SL(',ILEVEL,')%YY=',(SL%YY(I),I=0,0)
      WRITE(9,'(a,i3,a,16f7.3)') 'SL(',ILEVEL,')%ZZ=',(SL%ZZ(I),I=0,SL%KBAR)
      WRITE(9,*) 
      WRITE(9,*) 'NMESHES=', NMESHES
      WRITE(9,*) 'NCELLS_LOCAL=', SL%NCELLS_LOCAL
      WRITE(9,*) 'NCELLS_GLOBAL=', SL%NCELLS_GLOBAL
      CALL flush (9)
   ENDIF

   IBAR0=IBAR0/2
   JBAR0=JBAR0/2
   KBAR0=KBAR0/2

ENDDO GRID_LEVEL_LOOP3D
 
TUSED_SCARC(1,NM)=TUSED_SCARC(1,NM)+SECOND()-TNOW_MESHES3D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_MESHES3D

END SUBROUTINE SCARC_INITIALIZE_MESHES3D
 
 
 
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! determine neighborship/communication structure for data exchange for 2D-case
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_INITIALIZE_NEIGHBORS2D(NM)
 
USE GEOMETRY_FUNCTIONS, ONLY: SEARCH_OTHER_MESHES

INTEGER :: NM, NOM, IERR=0, N
INTEGER :: IREFINE
INTEGER :: IVRTX, ILEVEL
INTEGER :: IW_LO
INTEGER :: NOM_HI1, NOM_HI2
INTEGER :: IMIN, IMAX, JMIN, JMAX, KMIN, KMAX
INTEGER :: IMIN1_HI, IMIN2_HI, KMIN1_HI, KMIN2_HI
INTEGER :: IMAX1_HI, IMAX2_HI, KMAX1_HI, KMAX2_HI
INTEGER :: ILMAX, IOR
INTEGER :: IW1_HI, IW2_HI
INTEGER :: IW_FACE, IW_DIAG
INTEGER :: IW, IM , MMM, NNN, ITER, JM
INTEGER :: I, I_LO, K_LO
INTEGER :: IIO, JJO, KKO, IIO_MIN, IIO_MAX, JJO_MIN, JJO_MAX, KKO_MIN, KKO_MAX, NCHECK(4)
INTEGER, ALLOCATABLE, DIMENSION(:) :: COUNTS,DISPLS,COUNTS2D,DISPLS2D
REAL(EB):: XIN, YIN, ZIN
REAL(EB):: TNOW_NEIGHBORS2D
TYPE (MESH_TYPE), POINTER :: M
 
INCLUDE 'mpif.h'
 
TNOW_NEIGHBORS2D = SECOND()

IF (SCARC_DEBUG >= 2) THEN
   WRITE (9,*) '==========================================================================='
   WRITE (9,*) '=== SCARC_INITIALIZE_NEIGHBORS2D ',NM
   WRITE (9,*) '==========================================================================='
   WRITE (9,*) 'NMESHES=', NMESHES
   WRITE (9,*) 'S%NLEVEL=', S%NLEVEL
   WRITE (9,*) 'HIER SCHAUEN, WELCHES GEBRAUCHT WIRD'
   WRITE (9,*) 'NEWC=', M%NEWC
   WRITE (9,*) 'NWC=', M%NWC
   CALL flush (9)
ENDIF
    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Determine initial values
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
M => MESHES(NM)

IBAR=M%IBAR
KBAR=M%KBAR

S%NFACE0 = 4
S%NVRTX0 = 4

ALLOCATE (S%NDIAG_NBR(2,2,S%NVRTX0))

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
   WRITE(9,*) 'NIC:'
   DO MMM=1,NMESHES
      WRITE(9,'(15i3)') (NIC(NNN,MMM),NNN=1,NMESHES)
   ENDDO
   WRITE(9,*) 'I_MIN:'
   DO MMM=1,NMESHES
      WRITE(9,'(15i3)') (I_MIN(NNN,MMM),NNN=1,NMESHES)
   ENDDO
   WRITE(9,*) 'I_MAX:'
   DO MMM=1,NMESHES
      WRITE(9,'(15i3)') (I_MAX(NNN,MMM),NNN=1,NMESHES)
   ENDDO
   WRITE(9,*) 'K_MIN:'
   DO MMM=1,NMESHES
      WRITE(9,'(15i3)') (K_MIN(NNN,MMM),NNN=1,NMESHES)
   ENDDO
   WRITE(9,*) 'K_MAX:'
   DO MMM=1,NMESHES
      WRITE(9,'(15i3)') (K_MAX(NNN,MMM),NNN=1,NMESHES)
   ENDDO
   WRITE(9,*)
   WRITE(9,*) 'NDIAG_NBR(1:2,1,1)=',S%NDIAG_NBR(1,1,1), S%NDIAG_NBR(2,1,1)
   WRITE(9,*) 'NDIAG_NBR(1:2,1,2)=',S%NDIAG_NBR(1,1,2), S%NDIAG_NBR(2,1,2)
   WRITE(9,*) 'NDIAG_NBR(1:2,1,3)=',S%NDIAG_NBR(1,1,3), S%NDIAG_NBR(2,1,3)
   WRITE(9,*) 'NDIAG_NBR(1:2,1,4)=',S%NDIAG_NBR(1,1,4), S%NDIAG_NBR(2,1,4)
   WRITE(9,*)
   WRITE(9,*) 'NDIAG_NBR(1:2,2,1)=',S%NDIAG_NBR(1,2,1), S%NDIAG_NBR(2,2,1)
   WRITE(9,*) 'NDIAG_NBR(1:2,2,2)=',S%NDIAG_NBR(1,2,2), S%NDIAG_NBR(2,2,2)
   WRITE(9,*) 'NDIAG_NBR(1:2,2,3)=',S%NDIAG_NBR(1,2,3), S%NDIAG_NBR(2,2,3)
   WRITE(9,*) 'NDIAG_NBR(1:2,2,4)=',S%NDIAG_NBR(1,2,4), S%NDIAG_NBR(2,2,4)
   flush(9)
ENDIF
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Determine arrays IJKW_FACE and IJKW_DIAG for finest level
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ILMAX =  S%NLEVEL_MAX
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

IF (SCARC_DEBUG>=2) WRITE(9,*) 'NM=',NM

SLMAX%NW_FACE = 2*SLMAX%IBAR + 2*SLMAX%KBAR
ALLOCATE (SLMAX%IJKW_FACE(17,SLMAX%NW_FACE))
SLMAX%IJKW_FACE = 0
   
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
         !S%IOR_FACE(NM,NOM)=ABS(M%IJKW(4,IW))
         S%IOR_FACE(NM,NOM)=M%IJKW(4,IW)
         !S%IOR_FACE(NOM,NM)=ABS(M%IJKW(4,IW))
         S%IOR_FACE(NOM,NM)=M%IJKW(4,IW)
      ENDIF

      IW_FACE = IW_FACE+1

   ENDIF

   IF (IW_FACE>SLMAX%NW_FACE) EXIT IJKW_FACE_LOOP2D

ENDDO IJKW_FACE_LOOP2D

! communicate NIC_DIAG-information from other meshes
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
      WRITE(9,*) '============== NOM=',NOM
      WRITE (9,*) 'NIC_FACE:'
      DO JM=1,NMESHES
         WRITE (9, '(10i5)') (SLMAX%NIC_FACE(IM,JM),IM=1,NMESHES)
      ENDDO
      WRITE(9,*) 'MIBAR(NOM)=',S%MIBAR(NOM)
      WRITE(9,*) 'MKBAR(NOM)=',S%MKBAR(NOM)
      WRITE(9,*) 'IREFINE=',IREFINE
      WRITE(9,*) 'IMIN=',IMIN
      WRITE(9,*) 'IMAX=',IMAX
      WRITE(9,*) 'KMIN=',KMIN
      WRITE(9,*) 'KMAX=',KMAX
   ENDIF

   SLMAX%NIC_FACE(NOM,NM) = 0

   SEARCH_LOOP_MAX: DO IW_FACE=1,SLMAX%NW_FACE

      ! neighborship structure already known from finest level
      IF (SLMAX%IJKW_FACE(9,IW_FACE)/=NOM) CYCLE SEARCH_LOOP_MAX

      !WRITE(9,*) 'IJKW_FACE(9,',IW_FACE,')=',SLMAX%IJKW_FACE(9,IW_FACE)

      SLMAX%NIC_FACE(NOM,NM) = SLMAX%NIC_FACE(NOM,NM) + 1
      IOR = SLMAX%IJKW_FACE(4,IW_FACE)

      !WRITE(9,*) 'NIC_FACE(',NOM,',',NM,')=',SLMAX%NIC_FACE(NOM,NM)

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

   IF (ABS(S%IOR_FACE(NOM,NM))==3) THEN               ! extended communication
      SLMAX%I_MIN_FACE(NOM,NM) = IMIN-1
      SLMAX%I_MAX_FACE(NOM,NM) = IMAX+1
   ELSE
      SLMAX%I_MIN_FACE(NOM,NM) = IMIN
      SLMAX%I_MAX_FACE(NOM,NM) = IMAX
   ENDIF
   SLMAX%J_MIN_FACE(NOM,NM) = JMIN
   SLMAX%J_MAX_FACE(NOM,NM) = JMAX
   IF (ABS(S%IOR_FACE(NOM,NM))==1) THEN               ! extended communication
      SLMAX%K_MIN_FACE(NOM,NM) = KMIN-1
      SLMAX%K_MAX_FACE(NOM,NM) = KMAX+1
   ELSE
      SLMAX%K_MIN_FACE(NOM,NM) = KMIN
      SLMAX%K_MAX_FACE(NOM,NM) = KMAX
   ENDIF
  
   SLMAX%I_MIN_FACE(NOM,NM) = IMIN
   SLMAX%I_MAX_FACE(NOM,NM) = IMAX
   SLMAX%K_MIN_FACE(NOM,NM) = KMIN
   SLMAX%K_MAX_FACE(NOM,NM) = KMAX

ENDDO OTHER_MESH_LOOP2D_MAX

IF (SCARC_DEBUG >= 2) THEN
   WRITE(9,*) 'NUMPROCS=',NUMPROCS
   WRITE(9,*) 'MYID  =',MYID  
   WRITE(9,*) 'COUNTS=',COUNTS
   WRITE(9,*) 'COUNTS2D=',COUNTS2D
   WRITE(9,*) 'DISPLS=',DISPLS
   WRITE(9,*) 'DISPLS2D=',DISPLS2D
   WRITE(9,*) 'DISPLS(MYID)+1=',DISPLS(MYID)+1
   WRITE(9,*) 'COUNTS2D(MYID)=',COUNTS2D(MYID)
   WRITE(9,*) 'I_MIN_FACE0:'
   WRITE(9,'(4i4)') ((SLMAX%I_MIN_FACE(IM,JM),JM=1,NMESHES),IM=1,NMESHES)
   WRITE(9,*) 'I_MAX0:'
   WRITE(9,'(4i4)') ((SLMAX%I_MAX_FACE(IM,JM),JM=1,NMESHES),IM=1,NMESHES)
   WRITE(9,*) 'J_MIN0:'
   WRITE(9,'(4i4)') ((SLMAX%J_MIN_FACE(IM,JM),JM=1,NMESHES),IM=1,NMESHES)
   WRITE(9,*) 'J_MAX0:'
   WRITE(9,'(4i4)') ((SLMAX%J_MAX_FACE(IM,JM),JM=1,NMESHES),IM=1,NMESHES)
   WRITE(9,*) 'K_MIN0:'
   WRITE(9,'(4i4)') ((SLMAX%K_MIN_FACE(IM,JM),JM=1,NMESHES),IM=1,NMESHES)
   WRITE(9,*) 'K_MAX0:'
   WRITE(9,'(4i4)') ((SLMAX%K_MAX_FACE(IM,JM),JM=1,NMESHES),IM=1,NMESHES)
   flush(9)
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
ENDIF

IF (SCARC_DEBUG >= 2) THEN
   WRITE(9,*) 'AFTER ALLGATHERV'
   WRITE(9,*) 'I_MIN_FACE0:'
   WRITE(9,'(4i4)') ((SLMAX%I_MIN_FACE(IM,JM),JM=1,NMESHES),IM=1,NMESHES)
   WRITE(9,*) 'I_MAX0:'
   WRITE(9,'(4i4)') ((SLMAX%I_MAX_FACE(IM,JM),JM=1,NMESHES),IM=1,NMESHES)
   WRITE(9,*) 'J_MIN0:'
   WRITE(9,'(4i4)') ((SLMAX%J_MIN_FACE(IM,JM),JM=1,NMESHES),IM=1,NMESHES)
   WRITE(9,*) 'J_MAX0:'
   WRITE(9,'(4i4)') ((SLMAX%J_MAX_FACE(IM,JM),JM=1,NMESHES),IM=1,NMESHES)
   WRITE(9,*) 'K_MIN0:'
   WRITE(9,'(4i4)') ((SLMAX%K_MIN_FACE(IM,JM),JM=1,NMESHES),IM=1,NMESHES)
   WRITE(9,*) 'K_MAX0:'
   WRITE(9,'(4i4)') ((SLMAX%K_MAX_FACE(IM,JM),JM=1,NMESHES),IM=1,NMESHES)
   flush(9)
ENDIF

SLMAX%I_MIN_FACE = TRANSPOSE(SLMAX%I_MIN_FACE)
SLMAX%I_MAX_FACE = TRANSPOSE(SLMAX%I_MAX_FACE)
SLMAX%J_MIN_FACE = TRANSPOSE(SLMAX%J_MIN_FACE)
SLMAX%J_MAX_FACE = TRANSPOSE(SLMAX%J_MAX_FACE)
SLMAX%K_MIN_FACE = TRANSPOSE(SLMAX%K_MIN_FACE)
SLMAX%K_MAX_FACE = TRANSPOSE(SLMAX%K_MAX_FACE)
SLMAX%NIC_FACE   = TRANSPOSE(SLMAX%NIC_FACE)

IF (SCARC_DEBUG>=2) THEN
   WRITE(9,*) 'AFTER TRANSPOSE'
   WRITE(9,*) 'NIC_FACE:'
   WRITE(9,'(2i4)') ((SLMAX%NIC_FACE(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
   WRITE(9,*) 'I_MIN_FACE:'
   WRITE(9,'(2i4)') ((SLMAX%I_MIN_FACE(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
   WRITE(9,*) 'I_MAX:'
   WRITE(9,'(2i4)') ((SLMAX%I_MAX_FACE(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
   WRITE(9,*) 'J_MIN:'
   WRITE(9,'(2i4)') ((SLMAX%J_MIN_FACE(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
   WRITE(9,*) 'J_MAX:'
   WRITE(9,'(2i4)') ((SLMAX%J_MAX_FACE(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
   WRITE(9,*) 'K_MIN:'
   WRITE(9,'(2i4)') ((SLMAX%K_MIN_FACE(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
   WRITE(9,*) 'K_MAX:'
   WRITE(9,'(2i4)') ((SLMAX%K_MAX_FACE(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
   WRITE(9,*) 'IOR_FACE:'
   WRITE(9,'(2i4)') ((S%IOR_FACE(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
   call flush(9)
ENDIF

    
!---------------------------------------------------------------------
! determine neighbors along vertices
!---------------------------------------------------------------------
ALLOCATE (SLMAX%NIC_DIAG(NMESHES, NMESHES))
SLMAX%NIC_DIAG = 0

SLMAX%NW_DIAG = 4
ALLOCATE (SLMAX%IJKW_DIAG(17,SLMAX%NW_DIAG))
SLMAX%IJKW_DIAG = 0

ALLOCATE (SLMAX%I_MIN_DIAG(NMESHES, NMESHES))
ALLOCATE (SLMAX%I_MAX_DIAG(NMESHES, NMESHES))
ALLOCATE (SLMAX%J_MIN_DIAG(NMESHES, NMESHES))
ALLOCATE (SLMAX%J_MAX_DIAG(NMESHES, NMESHES))
ALLOCATE (SLMAX%K_MIN_DIAG(NMESHES, NMESHES))
ALLOCATE (SLMAX%K_MAX_DIAG(NMESHES, NMESHES))
SLMAX%I_MIN_DIAG = 0
SLMAX%I_MAX_DIAG = 0
SLMAX%J_MIN_DIAG = 0
SLMAX%J_MAX_DIAG = 0
SLMAX%K_MIN_DIAG = 0
SLMAX%K_MAX_DIAG = 0
   
! lower left diagonal node
SLMAX%IJKW_DIAG( 1,1)=0
SLMAX%IJKW_DIAG( 2,1)=1
SLMAX%IJKW_DIAG( 3,1)=0
SLMAX%IJKW_DIAG( 4,1)=1
SLMAX%IJKW_DIAG( 6,1)=1
SLMAX%IJKW_DIAG( 7,1)=1
SLMAX%IJKW_DIAG( 8,1)=1
SLMAX%IJKW_DIAG(16,1)=S%NDIAG_NBR(1,2,1)
SLMAX%IJKW_DIAG(17,1)=S%NDIAG_NBR(2,2,1)

! lower left diagonal node
SLMAX%IJKW_DIAG( 1,2)=SLMAX%IBAR + 1
SLMAX%IJKW_DIAG( 2,2)=1
SLMAX%IJKW_DIAG( 3,2)=0
SLMAX%IJKW_DIAG( 4,2)=2
SLMAX%IJKW_DIAG( 6,2)=SLMAX%IBAR
SLMAX%IJKW_DIAG( 7,2)=1
SLMAX%IJKW_DIAG( 8,2)=1
SLMAX%IJKW_DIAG(16,2)=S%NDIAG_NBR(1,2,2)
SLMAX%IJKW_DIAG(17,2)=S%NDIAG_NBR(2,2,2)

! lower left diagonal node
SLMAX%IJKW_DIAG( 1,3)=SLMAX%IBAR + 1
SLMAX%IJKW_DIAG( 2,3)=1
SLMAX%IJKW_DIAG( 3,3)=SLMAX%KBAR + 1
SLMAX%IJKW_DIAG( 4,3)=3
SLMAX%IJKW_DIAG( 6,3)=SLMAX%IBAR
SLMAX%IJKW_DIAG( 7,3)=1
SLMAX%IJKW_DIAG( 8,3)=SLMAX%KBAR
SLMAX%IJKW_DIAG(16,3)=S%NDIAG_NBR(1,2,3)
SLMAX%IJKW_DIAG(17,3)=S%NDIAG_NBR(2,2,3)

! lower left diagonal node
SLMAX%IJKW_DIAG( 1,4)=0
SLMAX%IJKW_DIAG( 2,4)=1
SLMAX%IJKW_DIAG( 3,4)=SLMAX%KBAR + 1
SLMAX%IJKW_DIAG( 4,4)=4
SLMAX%IJKW_DIAG( 6,4)=1
SLMAX%IJKW_DIAG( 7,4)=1
SLMAX%IJKW_DIAG( 8,4)=SLMAX%KBAR
SLMAX%IJKW_DIAG(16,4)=S%NDIAG_NBR(1,2,4)
SLMAX%IJKW_DIAG(17,4)=S%NDIAG_NBR(2,2,4)

      
! search for diagonal neighbors
SEARCH_DIAG_NBR2D: DO IVRTX=1,S%NVRTX0

   if (SCARC_DEBUG>=2) THEN
      WRITE(9,*) '============================== IVRTX=',IVRTX
      WRITE(9,*) 'SLMAX%IBAR=',SLMAX%IBAR
      WRITE(9,*) 'SLMAX%KBAR=',SLMAX%KBAR
      WRITE(9,*) 'DX(0)   =',M%DX(0)
      WRITE(9,*) 'DX(IBAR)=',M%DX(SLMAX%IBAR)
      WRITE(9,*) 'DZ(0)   =',M%DZ(0)
      WRITE(9,*) 'DZ(KBAR)=',M%DZ(SLMAX%KBAR)
   endif

   IIO_MIN =  1000000
   IIO_MAX = -1000000
   JJO_MIN =  1000000
   JJO_MAX = -1000000
   KKO_MIN =  1000000
   KKO_MAX = -1000000

   NCHECK = 0
   YIN=0.0_EB

   DO ITER=1,4

      ! Lower left vertex
      IF (IVRTX==1) THEN

         IF (ITER==1) THEN
            XIN = M%X(0) - 0.45_EB*M%DX(0)
            ZIN = M%Z(0) - 0.05_EB*M%DZ(0) 
         ELSE IF (ITER==2) THEN
            XIN = M%X(0) - 0.95_EB*M%DX(0)
            ZIN = M%Z(0) - 0.05_EB*M%DZ(0) 
         ELSE IF (ITER==3) THEN
            XIN = M%X(0) - 0.05_EB*M%DX(0)
            ZIN = M%Z(0) - 0.45_EB*M%DZ(0) 
         ELSE IF (ITER==4) THEN
            XIN = M%X(0) - 0.05_EB*M%DX(0)
            ZIN = M%Z(0) - 0.95_EB*M%DZ(0) 
         ENDIF

      ! Lower right vertex
      ELSE IF (IVRTX==2) THEN

         IF (ITER==1) THEN
            XIN = M%X(SLMAX%IBAR) + 0.05_EB*M%DX(SLMAX%IBAR)
            ZIN = M%Z(0)          - 0.45_EB*M%DZ(0) 
         ELSE IF (ITER==2) THEN
            XIN = M%X(SLMAX%IBAR) + 0.05_EB*M%DX(SLMAX%IBAR)
            ZIN = M%Z(0)          - 0.95_EB*M%DZ(0) 
         ELSE IF (ITER==3) THEN
            XIN = M%X(SLMAX%IBAR) + 0.45_EB*M%DX(SLMAX%IBAR)
            ZIN = M%Z(0)          - 0.05_EB*M%DZ(0) 
         ELSE IF (ITER==4) THEN
            XIN = M%X(SLMAX%IBAR) + 0.95_EB*M%DX(SLMAX%IBAR)
            ZIN = M%Z(0)          - 0.05_EB*M%DZ(0) 
         ENDIF

      ! Upper right vertex
      ELSE IF (IVRTX==3) THEN

         IF (ITER==1) THEN
            XIN = M%X(SLMAX%IBAR) + 0.45_EB*M%DX(SLMAX%IBAR)
            ZIN = M%Z(SLMAX%KBAR) + 0.05_EB*M%DZ(SLMAX%KBAR) 
         ELSE IF (ITER==2) THEN
            XIN = M%X(SLMAX%IBAR) + 0.95_EB*M%DX(SLMAX%IBAR)
            ZIN = M%Z(SLMAX%KBAR) + 0.05_EB*M%DZ(SLMAX%KBAR) 
         ELSE IF (ITER==3) THEN
            XIN = M%X(SLMAX%IBAR) + 0.05_EB*M%DX(SLMAX%IBAR)
            ZIN = M%Z(SLMAX%KBAR) + 0.45_EB*M%DZ(SLMAX%KBAR) 
         ELSE IF (ITER==4) THEN
            XIN = M%X(SLMAX%IBAR) + 0.05_EB*M%DX(SLMAX%IBAR)
            ZIN = M%Z(SLMAX%KBAR) + 0.95_EB*M%DZ(SLMAX%KBAR) 
         ENDIF

      ! Upper left vertex
      ELSE IF (IVRTX==4) THEN

         IF (ITER==1) THEN
            XIN = M%X(0)          - 0.05_EB*M%DX(0)
            ZIN = M%Z(SLMAX%KBAR) + 0.45_EB*M%DZ(SLMAX%KBAR) 
         ELSE IF (ITER==2) THEN
            XIN = M%X(0)          - 0.05_EB*M%DX(0)
            ZIN = M%Z(SLMAX%KBAR) + 0.95_EB*M%DZ(SLMAX%KBAR) 
         ELSE IF (ITER==3) THEN
            XIN = M%X(0)          - 0.45_EB*M%DX(0)
            ZIN = M%Z(SLMAX%KBAR) + 0.05_EB*M%DZ(SLMAX%KBAR) 
         ELSE IF (ITER==4) THEN
            XIN = M%X(0)          - 0.95_EB*M%DX(0)
            ZIN = M%Z(SLMAX%KBAR) + 0.05_EB*M%DZ(SLMAX%KBAR) 
         ENDIF

      ELSE 
          WRITE(*,*) 'WRONG VALUE OF IVRTX=',IVRTX,' IN SEARCH_DIAG_LOOP' 
          WRITE(9,*) 'WRONG VALUE OF IVRTX=',IVRTX,' IN SEARCH_DIAG_LOOP' 
          EXIT
      ENDIF

      ! search for neighboring meshes in direction of corresponding vertex
      CALL SEARCH_OTHER_MESHES(XIN,YIN,ZIN,NOM,IIO,JJO,KKO)
      NCHECK(ITER) = NOM

      if (SCARC_DEBUG>=2) THEN
         WRITE(9,*) '========== ITER=',ITER
         WRITE(9,*) 'XIN=',XIN
         WRITE(9,*) 'ZIN=',ZIN
         WRITE(9,*) 'NOM=',NOM
      ENDIF

      IF (NOM/=0) THEN

         IIO_MIN = MIN(IIO_MIN,IIO)
         IIO_MAX = MAX(IIO_MAX,IIO)
         JJO_MIN = MIN(JJO_MIN,JJO)
         JJO_MAX = MAX(JJO_MAX,JJO)
         KKO_MIN = MIN(KKO_MIN,KKO)
         KKO_MAX = MAX(KKO_MAX,KKO)

         IF (SCARC_DEBUG>=2) THEN
            WRITE(9,*) '===> NOM=',NOM
            WRITE(9,*) 'IIO_MIN=',IIO_MIN
            WRITE(9,*) 'IIO_MAX=',IIO_MAX
            WRITE(9,*) 'KKO_MIN=',KKO_MIN
            WRITE(9,*) 'KKO_MAX=',KKO_MAX
            flush(9)
         ENDIF

      ENDIF

   ENDDO

   ! Check to see if the current interpolated cell face spans more than one other mesh
   IF (NCHECK(1)/=NCHECK(2).OR.NCHECK(1)/=NCHECK(3).OR.NCHECK(1)/=NCHECK(4)) THEN
      WRITE(LU_ERR,'(A,I3,A,I3)') 'ERROR: MESH ',NM,' is out of alignment with MESH ',MAXVAL(NCHECK)
      PROCESS_STOP_STATUS = SETUP_STOP
      IERR = 1
      RETURN
   ENDIF
      
   ! store indices of common cell
   IF (NOM>0) THEN

      IF (SCARC_DEBUG>=2) WRITE(9,*) 'DIAG: NOM=',NOM

      SLMAX%IJKW_DIAG( 9,IVRTX)=NOM

      SLMAX%IJKW_DIAG(10,IVRTX)=IIO_MIN
      SLMAX%IJKW_DIAG(11,IVRTX)=1
      SLMAX%IJKW_DIAG(12,IVRTX)=KKO_MIN

      SLMAX%IJKW_DIAG(13,IVRTX)=IIO_MAX
      SLMAX%IJKW_DIAG(14,IVRTX)=1
      SLMAX%IJKW_DIAG(15,IVRTX)=KKO_MAX

      SLMAX%NIC_DIAG(NOM, NM)=(IIO_MAX-IIO_MIN+1)*(KKO_MAX-KKO_MIN+1)

      IF (IVRTX==1) THEN
          SLMAX%I_MIN_DIAG(NOM,NM)=0
          SLMAX%I_MAX_DIAG(NOM,NM)=0
          SLMAX%J_MIN_DIAG(NOM,NM)=1
          SLMAX%J_MAX_DIAG(NOM,NM)=1
          SLMAX%K_MIN_DIAG(NOM,NM)=0
          SLMAX%K_MAX_DIAG(NOM,NM)=0
      ELSE IF (IVRTX==2) THEN
          SLMAX%I_MIN_DIAG(NOM,NM)=SLMAX%IBAR
          SLMAX%I_MAX_DIAG(NOM,NM)=SLMAX%IBAR
          SLMAX%J_MIN_DIAG(NOM,NM)=1
          SLMAX%J_MAX_DIAG(NOM,NM)=1
          SLMAX%K_MIN_DIAG(NOM,NM)=0
          SLMAX%K_MAX_DIAG(NOM,NM)=0
      ELSE IF (IVRTX==3) THEN
          SLMAX%I_MIN_DIAG(NOM,NM)=SLMAX%IBAR
          SLMAX%I_MAX_DIAG(NOM,NM)=SLMAX%IBAR
          SLMAX%J_MIN_DIAG(NOM,NM)=1
          SLMAX%J_MAX_DIAG(NOM,NM)=1
          SLMAX%K_MIN_DIAG(NOM,NM)=SLMAX%KBAR
          SLMAX%K_MAX_DIAG(NOM,NM)=SLMAX%KBAR
      ELSE IF (IVRTX==4) THEN
          SLMAX%I_MIN_DIAG(NOM,NM)=0
          SLMAX%I_MAX_DIAG(NOM,NM)=0
          SLMAX%J_MIN_DIAG(NOM,NM)=1
          SLMAX%J_MAX_DIAG(NOM,NM)=1
          SLMAX%K_MIN_DIAG(NOM,NM)=SLMAX%KBAR
          SLMAX%K_MAX_DIAG(NOM,NM)=SLMAX%KBAR
      ENDIF
     
   ENDIF

ENDDO SEARCH_DIAG_NBR2D

IF (SCARC_DEBUG >= 2) THEN
   WRITE(9,*) 'NUMPROCS=',NUMPROCS
   WRITE(9,*) 'MYID  =',MYID  
   WRITE(9,*) 'COUNTS=',COUNTS
   WRITE(9,*) 'COUNTS2D=',COUNTS2D
   WRITE(9,*) 'DISPLS=',DISPLS
   WRITE(9,*) 'DISPLS2D=',DISPLS2D
   WRITE(9,*) 'DISPLS(MYID)+1=',DISPLS(MYID)+1
   WRITE(9,*) 'COUNTS2D(MYID)=',COUNTS2D(MYID)
   WRITE(9,*) 'I_MIN_DIAG0:'
   WRITE(9,'(2i4)') ((SLMAX%I_MIN_DIAG(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
   WRITE(9,*) 'I_MAX_DIAG0:'
   WRITE(9,'(2i4)') ((SLMAX%I_MAX_DIAG(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
   WRITE(9,*) 'J_MIN_DIAG0:'
   WRITE(9,'(2i4)') ((SLMAX%J_MIN_DIAG(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
   WRITE(9,*) 'J_MAX_DIAG0:'
   WRITE(9,'(2i4)') ((SLMAX%J_MAX_DIAG(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
   WRITE(9,*) 'K_MIN_DIAG0:'
   WRITE(9,'(2i4)') ((SLMAX%K_MIN_DIAG(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
   WRITE(9,*) 'K_MAX_DIAG0:'
   WRITE(9,'(2i4)') ((SLMAX%K_MAX_DIAG(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
   WRITE (9,*) 'NIC_DIAG(1,DISPLS(MYID)+1:DISPLS(MYID)+3)=',SLMAX%NIC_DIAG(1,DISPLS(MYID)+1:DISPLS(MYID)+3)
   WRITE (9,*) '=======================A1:'
   WRITE (9,*) 'NIC_DIAG1:'
   write (9,*) 'USE_MPI=',USE_MPI
   DO JM=1,NMESHES
      WRITE (9, '(2i4)') (SLMAX%NIC_DIAG(IM,JM),IM=1,NMESHES)
   ENDDO
   write (9,*) 'USE_MPI=',USE_MPI
   flush(9)
ENDIF

IF (USE_MPI) THEN
   write (9,*) 'USE_MPI2=',USE_MPI
   CALL MPI_ALLGATHERV(SLMAX%I_MIN_DIAG(1,DISPLS(MYID)+1),COUNTS2D(MYID),MPI_INTEGER,&
                       SLMAX%I_MIN_DIAG,COUNTS2D,DISPLS2D,MPI_INTEGER,MPI_COMM_WORLD,IERR)
   CALL MPI_ALLGATHERV(SLMAX%I_MAX_DIAG(1,DISPLS(MYID)+1),COUNTS2D(MYID),MPI_INTEGER,&
                       SLMAX%I_MAX_DIAG,COUNTS2D,DISPLS2D,MPI_INTEGER,MPI_COMM_WORLD,IERR)
   CALL MPI_ALLGATHERV(SLMAX%J_MIN_DIAG(1,DISPLS(MYID)+1),COUNTS2D(MYID),MPI_INTEGER,&
                       SLMAX%J_MIN_DIAG,COUNTS2D,DISPLS2D,MPI_INTEGER,MPI_COMM_WORLD,IERR)
   CALL MPI_ALLGATHERV(SLMAX%J_MAX_DIAG(1,DISPLS(MYID)+1),COUNTS2D(MYID),MPI_INTEGER,&
                       SLMAX%J_MAX_DIAG,COUNTS2D,DISPLS2D,MPI_INTEGER,MPI_COMM_WORLD,IERR)
   CALL MPI_ALLGATHERV(SLMAX%K_MIN_DIAG(1,DISPLS(MYID)+1),COUNTS2D(MYID),MPI_INTEGER,&
                       SLMAX%K_MIN_DIAG,COUNTS2D,DISPLS2D,MPI_INTEGER,MPI_COMM_WORLD,IERR)
   CALL MPI_ALLGATHERV(SLMAX%K_MAX_DIAG(1,DISPLS(MYID)+1),COUNTS2D(MYID),MPI_INTEGER,&
                       SLMAX%K_MAX_DIAG,COUNTS2D,DISPLS2D,MPI_INTEGER,MPI_COMM_WORLD,IERR)
   CALL MPI_ALLGATHERV(SLMAX%NIC_DIAG(1,DISPLS(MYID)+1),  COUNTS2D(MYID),MPI_INTEGER,&
                       SLMAX%NIC_DIAG,  COUNTS2D,DISPLS2D,MPI_INTEGER,MPI_COMM_WORLD,IERR)
ENDIF


IF (SCARC_DEBUG >= 2) THEN
   WRITE (9,*) '=======================B:'
   WRITE (9,*) 'NIC_DIAG2.2:'
   DO JM=1,NMESHES
      WRITE (9, '(10i5)') (SLMAX%NIC_DIAG(IM,JM),IM=1,NMESHES)
   ENDDO
   WRITE(9,*) '======================================== ILEVEL ',ILMAX
   WRITE(9,*) 'SLMAX%NW_FACE=',SLMAX%NW_FACE
   WRITE(9,*) 'SLMAX%NW_DIAG=',SLMAX%NW_DIAG
   WRITE(9,*) 'SLMAX%IBAR   =',SLMAX%IBAR 
   WRITE(9,*) 'SLMAX%KBAR   =',SLMAX%KBAR 
   WRITE(9,*) 'S%NLEVEL    =',S%NLEVEL
   WRITE (9,*) 'IJKW_FACE:'
   DO IW_FACE=1,SLMAX%NW_FACE
      WRITE (9, '(18i5)') IW_FACE,(SLMAX%IJKW_FACE(I, IW_FACE), I=1, 17)
   ENDDO
   WRITE (9,*) 'IJKW_DIAG:'
   DO IW_DIAG=1,SLMAX%NW_DIAG
      WRITE (9, '(18i5)') IW_DIAG,(SLMAX%IJKW_DIAG(I, IW_DIAG), I=1, 17)
   ENDDO
   WRITE (9,*) 'NIC_FACE:'
   DO JM=1,NMESHES
      WRITE (9, '(10i5)') (SLMAX%NIC_FACE(IM,JM),IM=1,NMESHES)
   ENDDO
   WRITE (9,*) '=======================C:'
   WRITE (9,*) 'NIC_DIAG:'
   DO JM=1,NMESHES
      WRITE (9, '(2i4)') (SLMAX%NIC_DIAG(IM,JM),IM=1,NMESHES)
   ENDDO
   WRITE(9,*) 'I_MIN_DIAG:'
   WRITE(9,'(2i4)') ((SLMAX%I_MIN_DIAG(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
   WRITE(9,*) 'I_MAX_DIAG:'
   WRITE(9,'(2i4)') ((SLMAX%I_MAX_DIAG(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
   WRITE(9,*) 'J_MIN_DIAG:'
   WRITE(9,'(2i4)') ((SLMAX%J_MIN_DIAG(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
   WRITE(9,*) 'J_MAX_DIAG:'
   WRITE(9,'(2i4)') ((SLMAX%J_MAX_DIAG(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
   WRITE(9,*) 'K_MIN_DIAG:'
   WRITE(9,'(2i4)') ((SLMAX%K_MIN_DIAG(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
   WRITE(9,*) 'K_MAX_DIAG:'
   WRITE(9,'(2i4)') ((SLMAX%K_MAX_DIAG(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
   WRITE (9,*) 'SCARC_METHOD:', SCARC_METHOD
   DO IM=1,NMESHES
      IF (IM/=NM) THEN
         WRITE (9,*) 'FOR MESH ',IM,' ALLOCATING: (',&
                     SLMAX%I_MIN_FACE(NM,IM),':',SLMAX%I_MAX_FACE(NM,IM),' , ',&
                     SLMAX%J_MIN_FACE(NM,IM),':',SLMAX%J_MAX_FACE(NM,IM),' , ',&
                     SLMAX%K_MIN_FACE(NM,IM),':',SLMAX%K_MAX_FACE(NM,IM),')'
      ENDIF
   ENDDO
   flush(9)
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Only in case of MG-method:
!!! Determine arrays IJKW_FACE and IJKW_DIAG for coarser levels
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (SCARC_METHOD==2) THEN

   IREFINE=1
   INIT_NBR_LEVEL2D: DO ILEVEL=S%NLEVEL_MAX-1,S%NLEVEL_MIN,-1
   
      IREFINE=IREFINE*2

      SLHI => S%SLEVEL(ILEVEL+1)
      SLLO => S%SLEVEL(ILEVEL)
   
      SLLO%NW_FACE = 2*SLLO%IBAR + 2*SLLO%KBAR
      ALLOCATE (SLLO%IJKW_FACE(17,SLLO%NW_FACE))
      SLLO%IJKW_FACE=0
   
      IF (SCARC_DEBUG>=2) THEN
         WRITE(9,*) '========== ILEVEL ',ILEVEL
         WRITE(9,*) 'IREFINE=',IREFINE
         WRITE(9,*) 'SLLO%NW_FACE=',SLLO%NW_FACE
         WRITE(9,*) 'SLLO%IBAR       =',SLLO%IBAR 
         WRITE(9,*) 'SLLO%KBAR       =',SLLO%KBAR 
         flush(9)
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
         
               
               ! same KBAR-resolution on current mesh and on neighbor
               IF ((KMAX2_HI-KMAX1_HI)==1) THEN
                  SLLO%IJKW_FACE(10,IW_LO)=CEILING(REAL(IMIN2_HI/2,EB))
                  SLLO%IJKW_FACE(11,IW_LO)=1
                  SLLO%IJKW_FACE(12,IW_LO)=CEILING(REAL(KMIN2_HI/2,EB))
                  SLLO%IJKW_FACE(13,IW_LO)=CEILING(REAL(IMIN2_HI/2,EB))
                  SLLO%IJKW_FACE(14,IW_LO)=1
                  SLLO%IJKW_FACE(15,IW_LO)=CEILING(REAL(KMIN2_HI/2,EB))
      
               ! KBAR-resolution coarser on current mesh than on neighbor
               ELSE IF ((KMAX2_HI-KMAX1_HI)==2) THEN
                  SLLO%IJKW_FACE(10,IW_LO)=CEILING(REAL(IMAX1_HI,EB)/2)
                  SLLO%IJKW_FACE(11,IW_LO)=1
                  SLLO%IJKW_FACE(12,IW_LO)=CEILING(REAL(KMAX1_HI,EB)/2)
                  SLLO%IJKW_FACE(13,IW_LO)=CEILING(REAL(IMAX2_HI,EB)/2)
                  SLLO%IJKW_FACE(14,IW_LO)=1
                  SLLO%IJKW_FACE(15,IW_LO)=CEILING(REAL(KMAX2_HI,EB)/2)
      
               ! KBAR-resolution finer on current mesh than on neighbor
               ELSE IF ((KMAX2_HI-KMAX1_HI)==0) THEN
                  SLLO%IJKW_FACE(10,IW_LO)=CEILING(REAL(IMAX2_HI,EB)/2)
                  SLLO%IJKW_FACE(11,IW_LO)=1
                  SLLO%IJKW_FACE(12,IW_LO)=CEILING(REAL(KMAX2_HI,EB)/2)
                  SLLO%IJKW_FACE(13,IW_LO)=CEILING(REAL(IMAX2_HI,EB)/2)
                  SLLO%IJKW_FACE(14,IW_LO)=1
                  SLLO%IJKW_FACE(15,IW_LO)=CEILING(REAL(KMAX2_HI,EB)/2)
      
               ELSE                                         
                  WRITE(*,*) 'Error in INIT_NBR_LEVEL2D at IOR=1 for ILEVEL=',ILEVEL
                  WRITE(9,*) 'Error in INIT_NBR_LEVEL2D at IOR=1 for ILEVEL=',ILEVEL
                  STOP
               ENDIF
         
               IF (SCARC_DEBUG.ge.2) THEN
                  WRITE(9,*) '==== K_LO=',K_LO,' ============================='
                  WRITE(9,*) 'NOM_HI1 =',NOM_HI1
                  WRITE(9,*) 'NOM_HI2 =',NOM_HI2
                  WRITE(9,*) 
                  WRITE(9,*) 'IW1_HI  =',IW1_HI
                  WRITE(9,*) 'IW2_HI  =',IW2_HI
                  WRITE(9,*) 
                  WRITE(9,*) 'IMIN1_HI=',IMIN1_HI
                  WRITE(9,*) 'IMAX1_HI=',IMAX1_HI
                  WRITE(9,*) 
                  WRITE(9,*) 'IMIN2_HI=',IMIN2_HI
                  WRITE(9,*) 'IMAX2_HI=',IMAX2_HI
                  WRITE(9,*) 
                  WRITE(9,*) 'KMIN1_HI=',KMIN1_HI
                  WRITE(9,*) 'KMAX1_HI=',KMAX1_HI
                  WRITE(9,*) 
                  WRITE(9,*) 'KMIN2_HI=',KMIN2_HI
                  WRITE(9,*) 'KMAX2_HI=',KMAX2_HI
                  WRITE(9,*) 
                  WRITE (9, '(i5,a,17i5)') IW_LO,':',(SLLO%IJKW_FACE(I, IW_LO), I=1,17)
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
                  SLLO%IJKW_FACE(10,IW_LO)=CEILING(REAL(IMIN2_HI/2,EB))
                  SLLO%IJKW_FACE(11,IW_LO)=1
                  SLLO%IJKW_FACE(12,IW_LO)=CEILING(REAL(KMIN2_HI/2,EB))
                  SLLO%IJKW_FACE(13,IW_LO)=CEILING(REAL(IMIN2_HI/2,EB))
                  SLLO%IJKW_FACE(14,IW_LO)=1
                  SLLO%IJKW_FACE(15,IW_LO)=CEILING(REAL(KMIN2_HI/2,EB))
      
               ! KBAR-resolution coarser on current mesh than on neighbor
               ELSE IF ((KMAX2_HI-KMAX1_HI)==2) THEN
                  SLLO%IJKW_FACE(10,IW_LO)=CEILING(REAL(IMAX1_HI,EB)/2)
                  SLLO%IJKW_FACE(11,IW_LO)=1
                  SLLO%IJKW_FACE(12,IW_LO)=CEILING(REAL(KMAX1_HI,EB)/2)
                  SLLO%IJKW_FACE(13,IW_LO)=CEILING(REAL(IMAX2_HI,EB)/2)
                  SLLO%IJKW_FACE(14,IW_LO)=1
                  SLLO%IJKW_FACE(15,IW_LO)=CEILING(REAL(KMAX2_HI,EB)/2)
      
               ! KBAR-resolution finer on current mesh than on neighbor
               ELSE IF ((KMAX2_HI-KMAX1_HI)==0) THEN
                  SLLO%IJKW_FACE(10,IW_LO)=CEILING(REAL(IMAX2_HI,EB)/2)
                  SLLO%IJKW_FACE(11,IW_LO)=1
                  SLLO%IJKW_FACE(12,IW_LO)=CEILING(REAL(KMAX2_HI,EB)/2)
                  SLLO%IJKW_FACE(13,IW_LO)=CEILING(REAL(IMAX2_HI,EB)/2)
                  SLLO%IJKW_FACE(14,IW_LO)=1
                  SLLO%IJKW_FACE(15,IW_LO)=CEILING(REAL(KMAX2_HI,EB)/2)
      
               ELSE                                         
                  WRITE(*,*) 'Error in INIT_NBR_LEVEL2D at IOR=-1 for ILEVEL=',ILEVEL
                  WRITE(9,*) 'Error in INIT_NBR_LEVEL2D at IOR=-1 for ILEVEL=',ILEVEL
                  STOP
               ENDIF
         
               IF (SCARC_DEBUG.ge.2) THEN
                  WRITE(9,*) '==== K_LO=',K_LO,' ============================='
                  WRITE(9,*) 'NOM_HI1 =',NOM_HI1
                  WRITE(9,*) 'NOM_HI2 =',NOM_HI2
                  WRITE(9,*) 
                  WRITE(9,*) 'IW1_HI  =',IW1_HI
                  WRITE(9,*) 'IW2_HI  =',IW2_HI
                  WRITE(9,*) 
                  WRITE(9,*) 'IMIN1_HI=',IMIN1_HI
                  WRITE(9,*) 'IMAX1_HI=',IMAX1_HI
                  WRITE(9,*) 
                  WRITE(9,*) 'IMIN2_HI=',IMIN2_HI
                  WRITE(9,*) 'IMAX2_HI=',IMAX2_HI
                  WRITE(9,*) 
                  WRITE(9,*) 'KMIN1_HI=',KMIN1_HI
                  WRITE(9,*) 'KMAX1_HI=',KMAX1_HI
                  WRITE(9,*) 
                  WRITE(9,*) 'KMIN2_HI=',KMIN2_HI
                  WRITE(9,*) 'KMAX2_HI=',KMAX2_HI
                  WRITE(9,*) 
                  WRITE (9, '(i5,a,17i5)') IW_LO,':',(SLLO%IJKW_FACE(I, IW_LO), I=1,17)
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
               
               ! same KBAR-resolution on current mesh and on neighbor
               IF ((IMAX2_HI-IMAX1_HI)==1) THEN
                  SLLO%IJKW_FACE(10,IW_LO)=CEILING(REAL(IMIN2_HI/2,EB))
                  SLLO%IJKW_FACE(11,IW_LO)=1
                  SLLO%IJKW_FACE(12,IW_LO)=CEILING(REAL(KMIN2_HI/2,EB))
                  SLLO%IJKW_FACE(13,IW_LO)=CEILING(REAL(IMIN2_HI/2,EB))
                  SLLO%IJKW_FACE(14,IW_LO)=1
                  SLLO%IJKW_FACE(15,IW_LO)=CEILING(REAL(KMIN2_HI/2,EB))
      
               ! KBAR-resolution coarser on current mesh than on neighbor
               ELSE IF ((IMAX2_HI-IMAX1_HI)==2) THEN
                  SLLO%IJKW_FACE(10,IW_LO)=CEILING(REAL(IMAX1_HI,EB)/2)
                  SLLO%IJKW_FACE(11,IW_LO)=1
                  SLLO%IJKW_FACE(12,IW_LO)=CEILING(REAL(KMAX1_HI,EB)/2)
                  SLLO%IJKW_FACE(13,IW_LO)=CEILING(REAL(IMAX2_HI,EB)/2)
                  SLLO%IJKW_FACE(14,IW_LO)=1
                  SLLO%IJKW_FACE(15,IW_LO)=CEILING(REAL(KMAX2_HI,EB)/2)
      
               ! KBAR-resolution finer on current mesh than on neighbor
               ELSE IF ((IMAX2_HI-IMAX1_HI)==0) THEN
                  SLLO%IJKW_FACE(10,IW_LO)=CEILING(REAL(IMAX2_HI,EB)/2)
                  SLLO%IJKW_FACE(11,IW_LO)=1
                  SLLO%IJKW_FACE(12,IW_LO)=CEILING(REAL(KMAX2_HI,EB)/2)
                  SLLO%IJKW_FACE(13,IW_LO)=CEILING(REAL(IMAX2_HI,EB)/2)
                  SLLO%IJKW_FACE(14,IW_LO)=1
                  SLLO%IJKW_FACE(15,IW_LO)=CEILING(REAL(KMAX2_HI,EB)/2)
      
               ELSE                                         
                  WRITE(*,*) 'Error in INIT_NBR_LEVEL2D at IOR=3 for ILEVEL=',ILEVEL
                  WRITE(9,*) 'Error in INIT_NBR_LEVEL2D at IOR=3 for ILEVEL=',ILEVEL
                  STOP
               ENDIF
         
               IF (SCARC_DEBUG.ge.2) THEN
                  WRITE(9,*) '==== I_LO=',I_LO,' ============================='
                  WRITE(9,*) 'NOM_HI1 =',NOM_HI1
                  WRITE(9,*) 'NOM_HI2 =',NOM_HI2
                  WRITE(9,*) 
                  WRITE(9,*) 'IW1_HI  =',IW1_HI
                  WRITE(9,*) 'IW2_HI  =',IW2_HI
                  WRITE(9,*) 
                  WRITE(9,*) 'IMIN1_HI=',IMIN1_HI
                  WRITE(9,*) 'IMAX1_HI=',IMAX1_HI
                  WRITE(9,*) 
                  WRITE(9,*) 'IMIN2_HI=',IMIN2_HI
                  WRITE(9,*) 'IMAX2_HI=',IMAX2_HI
                  WRITE(9,*) 
                  WRITE(9,*) 'KMIN1_HI=',KMIN1_HI
                  WRITE(9,*) 'KMAX1_HI=',KMAX1_HI
                  WRITE(9,*) 
                  WRITE(9,*) 'KMIN2_HI=',KMIN2_HI
                  WRITE(9,*) 'KMAX2_HI=',KMAX2_HI
                  WRITE(9,*) 
                  WRITE (9, '(i5,a,17i5)') IW_LO,':',(SLLO%IJKW_FACE(I, IW_LO), I=1,17)
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
                  SLLO%IJKW_FACE(10,IW_LO)=CEILING(REAL(IMAX2_HI,EB)/2)
                  SLLO%IJKW_FACE(11,IW_LO)=1
                  SLLO%IJKW_FACE(12,IW_LO)=CEILING(REAL(KMAX2_HI,EB)/2)
                  SLLO%IJKW_FACE(13,IW_LO)=CEILING(REAL(IMAX2_HI,EB)/2)
                  SLLO%IJKW_FACE(14,IW_LO)=1
                  SLLO%IJKW_FACE(15,IW_LO)=CEILING(REAL(KMAX2_HI,EB)/2)
      
               ! KBAR-resolution coarser on current mesh than on neighbor
               ELSE IF ((IMAX2_HI-IMAX1_HI)==2) THEN
                  SLLO%IJKW_FACE(10,IW_LO)=CEILING(REAL(IMAX1_HI,EB)/2)
                  SLLO%IJKW_FACE(11,IW_LO)=1
                  SLLO%IJKW_FACE(12,IW_LO)=CEILING(REAL(KMAX1_HI,EB)/2)
                  SLLO%IJKW_FACE(13,IW_LO)=CEILING(REAL(IMAX2_HI,EB)/2)
                  SLLO%IJKW_FACE(14,IW_LO)=1
                  SLLO%IJKW_FACE(15,IW_LO)=CEILING(REAL(KMAX2_HI,EB)/2)
      
               ! KBAR-resolution finer on current mesh than on neighbor
               ELSE IF ((IMAX2_HI-IMAX1_HI)==0) THEN
                  SLLO%IJKW_FACE(10,IW_LO)=CEILING(REAL(IMAX2_HI,EB)/2)
                  SLLO%IJKW_FACE(11,IW_LO)=1
                  SLLO%IJKW_FACE(12,IW_LO)=CEILING(REAL(KMAX2_HI,EB)/2)
                  SLLO%IJKW_FACE(13,IW_LO)=CEILING(REAL(IMAX2_HI,EB)/2)
                  SLLO%IJKW_FACE(14,IW_LO)=1
                  SLLO%IJKW_FACE(15,IW_LO)=CEILING(REAL(KMAX2_HI,EB)/2)
      
               ELSE                                         
                  WRITE(*,*) 'Error in INIT_NBR_LEVEL2D at IOR=-3 for ILEVEL=',ILEVEL
                  WRITE(9,*) 'Error in INIT_NBR_LEVEL2D at IOR=-3 for ILEVEL=',ILEVEL
                  STOP
               ENDIF
         
               IF (SCARC_DEBUG.ge.2) THEN
                  WRITE(9,*) '==== I_LO=',I_LO,' ============================='
                  WRITE(9,*) 'NOM_HI1 =',NOM_HI1
                  WRITE(9,*) 'NOM_HI2 =',NOM_HI2
                  WRITE(9,*) 
                  WRITE(9,*) 'IW1_HI  =',IW1_HI
                  WRITE(9,*) 'IW2_HI  =',IW2_HI
                  WRITE(9,*) 
                  WRITE(9,*) 'IMIN1_HI=',IMIN1_HI
                  WRITE(9,*) 'IMAX1_HI=',IMAX1_HI
                  WRITE(9,*) 
                  WRITE(9,*) 'IMIN2_HI=',IMIN2_HI
                  WRITE(9,*) 'IMAX2_HI=',IMAX2_HI
                  WRITE(9,*) 
                  WRITE(9,*) 'KMIN1_HI=',KMIN1_HI
                  WRITE(9,*) 'KMAX1_HI=',KMAX1_HI
                  WRITE(9,*) 
                  WRITE(9,*) 'KMIN2_HI=',KMIN2_HI
                  WRITE(9,*) 'KMAX2_HI=',KMAX2_HI
                  WRITE(9,*) 
                  WRITE(9,*) 'Real-Kmax2_HI=',REAL(KMAX2_HI/2,EB)
                  WRITE(9,*) 'Ceil-Kmax2_HI=',CEILING(REAL(KMAX2_HI,EB)/2)
                  WRITE (9, '(i5,a,17i5)') IW_LO,':',(SLLO%IJKW_FACE(I, IW_LO), I=1,17)
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

         IF (SCARC_DEBUG.ge.2) THEN
            WRITE(9,*) '============== NOM=',NOM
            WRITE (9,*) 'NIC_FACE:'
            DO JM=1,NMESHES
               WRITE (9, '(10i5)') (SLMAX%NIC_FACE(IM,JM),IM=1,NMESHES)
            ENDDO
            WRITE(9,*) 'MIBAR(NOM)=',S%MIBAR(NOM)
            WRITE(9,*) 'MKBAR(NOM)=',S%MKBAR(NOM)
            WRITE(9,*) 'IREFINE=',IREFINE
            WRITE(9,*) 'IMAX=',IMAX
            WRITE(9,*) 'KMAX=',KMAX
         ENDIF

         SLLO%NIC_FACE(NOM,NM) = 0

         SEARCH_LOOP: DO IW_FACE=1,SLLO%NW_FACE

            ! neighborship structure already known from finest level
            IF (SLLO%IJKW_FACE(9,IW_FACE)/=NOM) CYCLE SEARCH_LOOP

            WRITE(9,*) 'IJKW_FACE(9,',IW_FACE,')=',SLLO%IJKW_FACE(9,IW_FACE)

            SLLO%NIC_FACE(NOM,NM) = SLLO%NIC_FACE(NOM,NM) + 1
            IOR = SLLO%IJKW_FACE(4,IW_FACE)

            WRITE(9,*) 'NIC_FACE(',NOM,',',NM,')=',SLLO%NIC_FACE(NOM,NM)

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
 
         SLLO%I_MIN_FACE(NOM,NM) = IMIN
         SLLO%I_MAX_FACE(NOM,NM) = IMAX
         SLLO%J_MIN_FACE(NOM,NM) = JMIN
         SLLO%J_MAX_FACE(NOM,NM) = JMAX
         SLLO%K_MIN_FACE(NOM,NM) = KMIN
         SLLO%K_MAX_FACE(NOM,NM) = KMAX

      ENDDO OTHER_MESH_LOOP2D

      IF (SCARC_DEBUG >= 2) THEN
         WRITE(9,*) 'NUMPROCS=',NUMPROCS
         WRITE(9,*) 'MYID  =',MYID  
         WRITE(9,*) 'COUNTS=',COUNTS
         WRITE(9,*) 'COUNTS2D=',COUNTS2D
         WRITE(9,*) 'DISPLS=',DISPLS
         WRITE(9,*) 'DISPLS2D=',DISPLS2D
         WRITE(9,*) 'DISPLS(MYID)+1=',DISPLS(MYID)+1
         WRITE(9,*) 'COUNTS2D(MYID)=',COUNTS2D(MYID)
         WRITE (9,*) 'I_MIN_FACE(1,DISPLS(MYID)+1:DISPLS(MYID)+3)=',SLLO%I_MIN_FACE(1,DISPLS(MYID)+1:DISPLS(MYID)+3)
         WRITE (9,*) '=======================A:'
         WRITE (9,*) 'NIC_DIAG2: USE_MPI', USE_MPI
         DO JM=1,NMESHES
            WRITE (9, '(10i5)') (SLLO%I_MIN_FACE(IM,JM),IM=1,NMESHES)
         ENDDO
         flush(9)
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
      ENDIF

      SLLO%I_MIN_FACE = TRANSPOSE(SLLO%I_MIN_FACE)
      SLLO%I_MAX_FACE = TRANSPOSE(SLLO%I_MAX_FACE)
      SLLO%J_MIN_FACE = TRANSPOSE(SLLO%J_MIN_FACE)
      SLLO%J_MAX_FACE = TRANSPOSE(SLLO%J_MAX_FACE)
      SLLO%K_MIN_FACE = TRANSPOSE(SLLO%K_MIN_FACE)
      SLLO%K_MAX_FACE = TRANSPOSE(SLLO%K_MAX_FACE)
      SLLO%NIC_FACE   = TRANSPOSE(SLLO%NIC_FACE)

      IF (SCARC_DEBUG >= 2) THEN
         WRITE (9,*) '=================================================='
         WRITE (9,*) '================== LEVEL =',ILEVEL
         WRITE (9,*) 'NFACE      =',S%NFACE0
         WRITE (9,*) 'NW_FACE=',SLLO%NW_FACE
         WRITE (9,*) 'IJKW_FACE:'
         DO IW_FACE=1,SLLO%NW_FACE
            WRITE (9, '(i5,a,17i5)') IW_FACE,':',(SLLO%IJKW_FACE(I, IW_FACE), I=1,17)
         ENDDO
         WRITE (9,*) 'NIC_FACE:'
         WRITE (9, '(3i4)') ((SLLO%NIC_FACE(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
         WRITE (9,*) 'I_MIN_FACE:'
         WRITE (9, '(3i4)') ((SLLO%I_MIN_FACE(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
         WRITE (9,*) 'I_MAX:'
         WRITE (9, '(3i4)') ((SLLO%I_MAX_FACE(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
         WRITE (9,*) 'J_MIN:'
         WRITE (9, '(3i4)') ((SLLO%J_MIN_FACE(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
         WRITE (9,*) 'J_MAX:'
         WRITE (9, '(3i4)') ((SLLO%J_MAX_FACE(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
         WRITE (9,*) 'K_MIN:'
         WRITE (9, '(3i4)') ((SLLO%K_MIN_FACE(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
         WRITE (9,*) 'K_MAX:'
         WRITE (9, '(3i4)') ((SLLO%K_MAX_FACE(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
         DO IM=1,NMESHES
            IF (IM/=NM) THEN
               WRITE (9,*) 'FOR MESH ',IM,' ALLOCATING: (',&
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

END SUBROUTINE SCARC_INITIALIZE_NEIGHBORS2D



 

 
 
 
 
!------------------------------------------------------------------------------------
! determine neighborship/communication structure for data exchange in 2D
!------------------------------------------------------------------------------------
SUBROUTINE SCARC_INITIALIZE_COMMUNICATION2D (NM)
 
INTEGER :: IM, NM, NOM, ILEVEL, IREFINE
INTEGER :: IMIN, IMAX, JMIN, JMAX, KMIN, KMAX, ILMAX
INTEGER :: TAG_FACE, TAG_DIAG
INTEGER :: I1, J1
INTEGER :: IZERO
REAL(EB):: TNOW_COMMUNICATION2D

INCLUDE 'mpif.h'
 
TNOW_COMMUNICATION2D = SECOND()
 
IF (SCARC_DEBUG >= 2) THEN
   WRITE (9,*) '=========================================================================='
   WRITE (9,*) '=== SCARC_INITIALIZE_COMMUNICATION2D ',NM
   WRITE (9,*) '=========================================================================='
   WRITE (9,*) 'NM=', NM
   CALL flush (9)
ENDIF
 

S => SCARC (NM)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Allocate protocol arrays for data exchange
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ALLOCATE (STAT_FACE(MPI_STATUS_SIZE, NMESHES*NMESHES*20))
ALLOCATE (STAT_DIAG(MPI_STATUS_SIZE, NMESHES*NMESHES*20))
 
ALLOCATE (REQ_FACE(NMESHES*NMESHES*20))
REQ_FACE = MPI_REQUEST_NULL

ALLOCATE (REQ_DIAG(NMESHES*NMESHES*20))
REQ_DIAG = MPI_REQUEST_NULL
 
ALLOCATE (TAGS_FACE(NMESHES, NMESHES))
TAGS_FACE = 0

ALLOCATE (TAGS_DIAG(NMESHES, NMESHES))
TAGS_DIAG = 0
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Create unique tags arrays for face, edge and vertex exchanges
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
TAG_FACE = 0
TAG_DIAG = 0
 
DO IM = 1, NMESHES
   DO NOM = IM, NMESHES
 
      TAG_FACE = TAG_FACE + 1
      TAG_DIAG = TAG_DIAG + 1
 
      TAGS_FACE (IM, NOM) = TAG_FACE
      TAGS_FACE (NOM, IM) = TAG_FACE
 
      TAGS_DIAG (IM, NOM) = TAG_DIAG
      TAGS_DIAG (NOM, IM) = TAG_DIAG
 
   ENDDO
ENDDO 
 
IF (SCARC_DEBUG >= 2) THEN
   WRITE (9,*) 'NM=',NM
   WRITE (9,*) 'TAGS_FACE:'
   WRITE (9, '(4i16)') ((TAGS_FACE(I1, J1), I1=1, NMESHES), J1=1, NMESHES)
   WRITE (9,*) 'TAGS_DIAG:'
   WRITE (9, '(4i16)') ((TAGS_DIAG(I1, J1), I1=1, NMESHES), J1=1, NMESHES)
ENDIF
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Initialize level structures on neighboring meshes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IZERO=0
OSLEVEL_LOOP2D: DO NOM=1,NMESHES

   ILMAX =  S%NLEVEL_MAX
   SLMAX => S%SLEVEL(ILMAX)

   IF (SLMAX%NIC_FACE(NM,NOM)/=0.AND.SLMAX%NIC_FACE(NOM,NM)/=0) THEN

      OS => S%OSCARC(NOM)
      ALLOCATE (OS%SLEVEL(S%NLEVEL_MIN:S%NLEVEL_MAX), STAT=IZERO)
      CALL CHKMEMERR ('SCARC_INITIALIZE_COMMUNICATION2D', 'OS%SLEVEL', IZERO)

      OSLMAX => OS%SLEVEL(ILMAX)
      OSLMAX%NW_FACE = 2*S%MIBAR(NOM) + 2*S%MKBAR(NOM)
      ALLOCATE (OSLMAX%IJKW_FACE(17,OSLMAX%NW_FACE), STAT=IZERO)
      CALL CHKMEMERR ('SCARC_INITIALIZE_COMMUNICATION2D', 'OSLMAX%IJKW_FACE', IZERO)
      OSLMAX%IJKW_FACE = 0

      IF (S%NLEVEL_DIFF/=0) THEN
         IREFINE=1
         DO ILEVEL=S%NLEVEL_MAX-1,S%NLEVEL_MIN

            IREFINE=IREFINE*2

            OS%SLEVEL(ILEVEL)%NW_FACE=OS%SLEVEL(ILEVEL+1)%NW_FACE/IREFINE
            ALLOCATE (OS%SLEVEL(ILEVEL)%IJKW_FACE(17,OSLMAX%NW_FACE), STAT=IZERO)
            CALL CHKMEMERR ('SCARC_INITIALIZE_COMMUNICATION2D', 'OS%SLEVEL%IJKW_FACE', IZERO)
            OS%SLEVEL(ILEVEL)%IJKW_FACE = 0

         ENDDO
      ENDIF

      IF (SCARC_DEBUG>=2) THEN
         WRITE(9,*) 'S%NLEVEL=',S%NLEVEL
         WRITE(9,*) 'S%NLEVEL_DIFF=',S%NLEVEL_DIFF
         WRITE(9,*) 'OS%SLEVEL(',ILMAX,')%NW_FACE=', OS%SLEVEL(ILMAX)%NW_FACE
         IF (S%NLEVEL_DIFF/=0) THEN
            DO ILEVEL=S%NLEVEL_MAX-1,S%NLEVEL_MIN
               WRITE(9,*) 'OS%SLEVEL(',ILEVEL,')%NW_FACE=', OS%SLEVEL(ILEVEL)%NW_FACE
            ENDDO
         ENDIF
         call flush(9)
      ENDIF

   ENDIF

   IF (SLMAX%NIC_DIAG(NM,NOM)/=0.AND.SLMAX%NIC_DIAG(NOM,NM)/=0) THEN

      OS => S%OSCARC(NOM)
      ALLOCATE (OS%SLEVEL(S%NLEVEL_MIN:S%NLEVEL_MAX), STAT=IZERO)
      CALL CHKMEMERR ('SCARC_INITIALIZE_COMMUNICATION2D', 'OS%SLEVEL', IZERO)

      OSLMAX => OS%SLEVEL(ILMAX)
      OSLMAX%NW_DIAG = 4
      ALLOCATE (OSLMAX%IJKW_DIAG(17,OSLMAX%NW_DIAG), STAT=IZERO)
      CALL CHKMEMERR ('SCARC_INITIALIZE_COMMUNICATION2D', 'OSLMAX%IJKW_DIAG', IZERO)
      OSLMAX%IJKW_DIAG = 0

      IF (S%NLEVEL_DIFF/=0) THEN
         IREFINE=1
         DO ILEVEL=S%NLEVEL_MAX-1,S%NLEVEL_MIN

            IREFINE=IREFINE*2

            OS%SLEVEL(ILEVEL)%NW_DIAG=OS%SLEVEL(ILEVEL+1)%NW_DIAG/IREFINE
            ALLOCATE (OS%SLEVEL(ILEVEL)%IJKW_DIAG(17,OSLMAX%NW_DIAG), STAT=IZERO)
            CALL CHKMEMERR ('SCARC_INITIALIZE_COMMUNICATION2D', 'OS%SLEVEL%IJKW_DIAG', IZERO)
            OS%SLEVEL(ILEVEL)%IJKW_DIAG = 0

         ENDDO
      ENDIF
      IF (SCARC_DEBUG>=2) THEN
         WRITE(9,*) 'HIERHALLO'
         WRITE(9,*) 'S%NLEVEL=',S%NLEVEL
         WRITE(9,*) 'S%NLEVEL_DIFF=',S%NLEVEL_DIFF
         WRITE(9,*) 'OS%SLEVEL(',ILMAX,')%NW_DIAG=', OS%SLEVEL(ILMAX)%NW_DIAG
         IF (S%NLEVEL_DIFF/=0) THEN
            DO ILEVEL=S%NLEVEL_MAX-1,S%NLEVEL_MIN
               WRITE(9,*) 'OS%SLEVEL(',ILEVEL,')%NW_DIAG=', OS%SLEVEL(ILEVEL)%NW_DIAG
            ENDDO
         ENDIF
         call flush(9)
      ENDIF

   ENDIF

ENDDO OSLEVEL_LOOP2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Initialize communication structures on every level
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
COM_LEVEL_LOOP2D: DO ILEVEL=S%NLEVEL_MAX,S%NLEVEL_MIN,-1

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
      IF (SCARC_DEBUG >= 2) THEN
         WRITE (9,*) '---NM=', NM
         WRITE (9,*) '---NOM=', NOM
         WRITE (9,*) 'ILEVEL',ILEVEL
         WRITE (9,*) 'SL%I_MIN_FACE=', SL%I_MIN_FACE
         WRITE (9,*) 'SL%I_MAX_FACE=', SL%I_MAX_FACE
         WRITE (9,*) 'SL%K_MIN_FACE=', SL%K_MIN_FACE
         WRITE (9,*) 'SL%K_MAX_FACE=', SL%K_MAX_FACE
         WRITE (9,*) 'IMIN=', IMIN
         WRITE (9,*) 'IMAX=', IMAX
         WRITE (9,*) 'JMIN=', JMIN
         WRITE (9,*) 'JMAX=', JMAX
         WRITE (9,*) 'KMIN=', KMIN
         WRITE (9,*) 'KMAX=', KMAX
         WRITE (9,'(a,i2,a,i2,a,i2,a,a,i2,a,i2,a,i2,a,i2,a,i2,a,i2,a,i2,a,i2,a)') 'LEVEL ',ILEVEL,&
                     ': allocate FACE-vectors for SCARC(', NM, ')%OSCARC(', NOM, ')',&
                     '%VEC(',IMIN,':',IMAX,',',JMIN,':',JMAX,',',KMIN,':',KMAX,')'
         CALL flush (9)
      ENDIF
    
   ENDDO FACE_NBR_LOOP2D
 
   !!! ---------------------------------------------------------------------------
   !!! Communication structures for vertices
   !!! ---------------------------------------------------------------------------
   DIAG_NBR_LOOP2D: DO NOM = 1, NMESHES

      ! is there a vertex neighbor ?
      IF (NOM == NM) CYCLE DIAG_NBR_LOOP2D
      IF (SL%NIC_DIAG(NM, NOM) == 0 .AND. SL%NIC_DIAG(NOM, NM) == 0) CYCLE DIAG_NBR_LOOP2D

      BDIAG=.TRUE.

      ! set lengths for communication vectors 
      IMIN=1
      IMAX=1
      JMIN=1
      JMAX=1
      KMIN=1
      KMAX=1
    
      ! allocate communication vectors for data exchange with vertex neighbor NOM on level ILEVEL
      OSL =>S%OSCARC(NOM)%SLEVEL(ILEVEL)
   
      ALLOCATE (OSL%Y_DIAG(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
      OSL%Y_DIAG = 0.0_EB
    
      ALLOCATE (OSL%R_DIAG(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
      OSL%R_DIAG = 0.0_EB
    
      ALLOCATE (OSL%G_DIAG(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
      OSL%G_DIAG = 0.0_EB
   
      ALLOCATE (OSL%Z_DIAG(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
      OSL%Z_DIAG = 0.0_EB
   
      ! print information if requested    
      IF (SCARC_DEBUG >= 2) THEN
         WRITE (9,*) '---NM=', NM
         WRITE (9,*) '---NOM=', NOM
         WRITE (9,*) 'IMIN=', IMIN
         WRITE (9,*) 'IMAX=', IMAX
         WRITE (9,*) 'JMIN=', JMIN
         WRITE (9,*) 'JMAX=', JMAX
         WRITE (9,*) 'KMIN=', KMIN
         WRITE (9,*) 'KMAX=', KMAX
         WRITE (9,'(a,i2,a,i2,a,i2,a,a,i2,a,i2,a,i2,a,i2,a,i2,a,i2,a,i2,a,i2,a)') 'LEVEL ',ILEVEL,&
                     ': allocate VRTX-vectors for SCARC(', NM, ')%OSCARC(', NOM, ')',&
                     '%VEC(',IMIN,':',IMAX,',',JMIN,':',JMAX,',',KMIN,':',KMAX,')'
         CALL flush (9)
       ENDIF

   ENDDO DIAG_NBR_LOOP2D

ENDDO COM_LEVEL_LOOP2D
 
 
IF (SCARC_DEBUG >= 2) write (9,*) 'BEFORE SCARC_RECEIVE'
 
DO ILEVEL=S%NLEVEL_MAX,S%NLEVEL_MIN,-1
   CALL SCARC_RECEIVE (NCOM_INIT, ILEVEL)    
   CALL SCARC_EXCHANGE(NCOM_INIT, ILEVEL, NTYPE_NONE)
ENDDO
 
IF (SCARC_DEBUG >= 2) write (9,*) 'Leaving SCARC_INITIALIZE_COMMUNICATION2D'
 
TUSED_SCARC(3,NM)=TUSED_SCARC(3,NM)+SECOND()-TNOW_COMMUNICATION2D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_COMMUNICATION2D
 
END SUBROUTINE SCARC_INITIALIZE_COMMUNICATION2D
 
 
 
 
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
SUBROUTINE SCARC_INITIALIZE_MATRICES2D (NM)
 
INTEGER :: NM, NOM
INTEGER :: I, K, IC, IW, IZERO, IERR, ICELL, IOR, ILEVEL
REAL(EB):: DBC
REAL(EB):: TNOW_MATRICES2D
TYPE (MESH_TYPE), POINTER :: M
 
INCLUDE 'mpif.h'
 
TNOW_MATRICES2D = SECOND()

IERR = 0
 
M => MESHES(NM)
 

IF (SCARC_DEBUG >= 2) THEN
   WRITE (9,*) '==========================================================================='
   WRITE (9,*) '=== SCARC_INITIALIZE_MATRICES2D ',NM
   WRITE (9,*) '==========================================================================='
   WRITE (9,*) 'NLEVEL_MAX=',S%NLEVEL_MAX
   WRITE (9,*) 'NLEVEL_MIN=',S%NLEVEL_MIN
ENDIF

!!!
!!! Initialize matrices on each level (only 1 level for CG)
!!!
MATRIX_LEVEL_LOOP2D: DO ILEVEL=S%NLEVEL_MAX,S%NLEVEL_MIN,-1

   SL => S%SLEVEL(ILEVEL)

    
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Initialize matrix diagonals for global matrix
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   NDIAG = 5
    
   !!! Allocate full matrix corresponding to the band-wise storage technique
   ALLOCATE (SL%AG(1:SL%NCELLS_LOCAL, 1:NDIAG), STAT=IZERO)
   CALL CHKMEMERR ('SCARC_INITIALIZE_MATRICES2D', 'SL%AG', IZERO)
   SL%AG = 0.0_EB
    
   !!! Allocate local matrix corresponding to the band-wise storage technique
   !ALLOCATE (SL%AL(1:SL%NCELLS_LOCAL, 1:NDIAG), STAT=IZERO)
   !CALL CHKMEMERR ('SCARC_INITIALIZE_MATRICES2D', 'SL%AL', IZERO)
    
    
   SL%DXI2=1.0_EB/(SL%DX)**2
   SL%DZI2=1.0_EB/(SL%DZ)**2

   IF (SCARC_DEBUG >= 2) THEN
      WRITE (9,*) '========= LEVEL=',ILEVEL
      WRITE (9,*) ' IBAR=',SL%IBAR
      WRITE (9,*) ' KBAR=',SL%KBAR
      WRITE (9,*) ' DX  =',SL%DX  
      WRITE (9,*) ' DZ  =',SL%DZ  
      WRITE (9,*) ' DXI2=',SL%DXI2
      WRITE (9,*) ' DZI2=',SL%DZI2
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
            WRITE (9,*) 'D1: AG(', IC, ',1)=', SL%AG(IC, 1)
            WRITE (9,*) 'D2: AG(', IC, ',2)=', SL%AG(IC, 2)
            WRITE (9,*) 'D3: AG(', IC, ',3)=', SL%AG(IC, 3)
            WRITE (9,*) 'D4: AG(', IC, ',4)=', SL%AG(IC, 4)
            WRITE (9,*) 'D5: AG(', IC, ',5)=', SL%AG(IC, 5)
            WRITE (9,*) 
         ENDIF
    
      ENDDO
   ENDDO
    
   !!! local matrix AL is momentarily set to global matrix AG
   !SL%AL = SL%AG
    
   IF (ILEVEL/=S%NLEVEL_MIN) THEN
      IF (SCARC_DEBUG>=2) write(9,*) 'CALLING GSTRIX_INIT on level ', ILEVEL
      CALL SCARC_PRECON_GSTRIX2D_INIT(ILEVEL,NM)
   ENDIF

   IF (ILEVEL/=S%NLEVEL_MAX) CYCLE

   IF (SCARC_DEBUG>=2) THEN
      WRITE(9,*) 'BOUNDARY_TYPE'
      WRITE(9,'(16i4)') (M%BOUNDARY_TYPE(ICELL),ICELL=1,M%NEWC)
      WRITE(9,*) 'IJKW'
      DO IW=1,M%NEWC+4
         WRITE(9,'(i4,a,15i4)') IW,' : ',(M%IJKW(ICELL,IW),ICELL=1,15)
      ENDDO
   ENDIF

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

      if (Scarc_debug.ge.8) then
         WRITE(9,*) 'ICELL=',ICELL,' :I=',I,': K=',K,': IW=',IW,': BT=',M%BOUNDARY_TYPE(IW)
      endif


      !IF (M%BOUNDARY_TYPE(IW)==OPEN_BOUNDARY.OR.(M%BOUNDARY_TYPE(IW)==NULL_BOUNDARY.AND.NOM.NE.0)) THEN
      IF (M%BOUNDARY_TYPE(IW)==OPEN_BOUNDARY) THEN
          M%PRESSURE_BC_INDEX(IW)=DIRICHLET
      !   WRITE(9,*) 'D: PRESSURE_BC_INDEX(',IW,')=',M%PRESSURE_BC_INDEX(IW)
      ELSE IF (M%IJKW(9,IW)/=0) THEN
          M%PRESSURE_BC_INDEX(IW)=INTERNAL
      !   WRITE(9,*) 'I: PRESSURE_BC_INDEX(',IW,')=',M%PRESSURE_BC_INDEX(IW)
      ELSE IF (M%BOUNDARY_TYPE(IW)==NULL_BOUNDARY) THEN
          M%PRESSURE_BC_INDEX(IW)=DIRICHLET
      !   WRITE(9,*) 'D: PRESSURE_BC_INDEX(',IW,')=',M%PRESSURE_BC_INDEX(IW)
      ELSE
          M%PRESSURE_BC_INDEX(IW)=NEUMANN
      !   WRITE(9,*) 'N: PRESSURE_BC_INDEX(',IW,')=',M%PRESSURE_BC_INDEX(IW)
      ENDIF

      !
      ! set Dirichlet BC's for AG and AL in open boundary cells
      !
      !IF (M%BOUNDARY_TYPE(IW) == OPEN_BOUNDARY) THEN    
      IF (M%PRESSURE_BC_INDEX(IW) == DIRICHLET) THEN

         SL%AG(IC, 1) = SL%AG(IC, 1) - DBC
         !SL%AL(IC, 1) = SL%AL(IC, 1) - DBC

         IF (SCARC_DEBUG .GE. 4) &
           write (9,'(a,i3,a,i2,a,i3,a,f12.6,a,i3,a,i3)') 'IW=',IW, &
              ': Dirichlet ',IOR,':AG(',IC,',1)=',SL%AG(IC,1),': BT=',M%BOUNDARY_TYPE(IW),': NOM=',NOM

      !
      ! set Dirichlet BC's for AL (NOT for AG!) along internal boundaries
      !
      !ELSE IF (NOM /= 0)  THEN    
      ELSE IF (M%PRESSURE_BC_INDEX(IW) == INTERNAL) THEN

         !SL%AL(IC, 1) = SL%AL(IC, 1) - DBC
         IF (SCARC_DEBUG .GE. 4) &
           write (9,'(a,i3,a,i2,a,i3,a,f12.6,a,i3,a,i3)') 'IW=',IW, &
              ': Nothing   ',IOR,':AG(',IC,',1)=',SL%AG(IC,1),': BT=',M%BOUNDARY_TYPE(IW),': NOM=',NOM


      !
      ! set Neumann BC's for AG and AL at all other nodes
      !
      !ELSE 
      ELSE IF (M%PRESSURE_BC_INDEX(IW) == NEUMANN) THEN

         SL%AG(IC, 1) = SL%AG(IC, 1) + DBC
         !SL%AL(IC, 1) = SL%AL(IC, 1) + DBC

         IF (SCARC_DEBUG .GE. 4) &
           write (9,'(a,i3,a,i2,a,i3,a,f12.6,a,i3,a,i3)') 'IW=',IW, &
              ': Neumann   ',IOR,':AG(',IC,',1)=',SL%AG(IC,1),': BT=',M%BOUNDARY_TYPE(IW),': NOM=',NOM

      ENDIF

   ENDDO


   IF (SCARC_DEBUG >= 2) THEN
      WRITE (9,*) 'ASUB=', SL%ASUB
      WRITE (9,*)
      WRITE (9,*) 'GLOBAL MATRIX AG:'
      DO IC = 1, SL%NCELLS_LOCAL
         WRITE (9, '(5f12.4)') SL%AG(IC, 2),SL%AG(IC,3),SL%AG(IC,1),SL%AG(IC,4),SL%AG(IC,5)
         IF (Mod(IC, SL%IBAR) == 0) write (9,*) '---------------------------------------------&
                                                &--------------------'
      ENDDO
   ENDIF
    
   SL%ASUBX = SL%DXI2
   SL%ASUBZ = SL%DZI2
 
   IF (SCARC_DEBUG >= 2) THEN
     write (9,*) 'SL%ASUBX=', SL%ASUBX
     WRITE (9,*) 'LBC=', LBC
     WRITE (9,*) 'MBC=', MBC
     WRITE (9,*) 'NBC=', NBC
   ENDIF

ENDDO MATRIX_LEVEL_LOOP2D
 

TUSED_SCARC(4,NM)=TUSED_SCARC(4,NM)+SECOND()-TNOW_MATRICES2D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_MATRICES2D

END SUBROUTINE SCARC_INITIALIZE_MATRICES2D
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  3D-version
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_INITIALIZE_MATRICES3D (NM)
 
INTEGER :: NM, NOM
INTEGER :: I, J, K, IC, IW, IZERO, IERR, ICELL, IOR, ILEVEL
REAL(EB):: DBC
REAL(EB):: TNOW_MATRICES3D
TYPE (MESH_TYPE), POINTER :: M
 
INCLUDE 'mpif.h'
 
TNOW_MATRICES3D = SECOND()

IERR = 0
 
M => MESHES(NM)
 
IF (SCARC_DEBUG >= 2) THEN
   WRITE (9,*) '==========================================================================='
   WRITE (9,*) '=== SCARC_INITIALIZE_MATRICES3D ',NM
   WRITE (9,*) '==========================================================================='
   WRITE (9,*) 'NLEVEL_MAX=',S%NLEVEL_MAX
   WRITE (9,*) 'NLEVEL_MIN=',S%NLEVEL_MIN
ENDIF

!!!
!!! Initialize matrices on each level (only 1 level for CG)
!!!
MATRIX_LEVEL_LOOP3D: DO ILEVEL=S%NLEVEL_MAX,S%NLEVEL_MIN,-1

   SL => S%SLEVEL(ILEVEL)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Initialize matrix diagonals for global matrix
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   NDIAG = 7
    
   !!! Allocate full matrix corresponding to the band-wise storage technique
   ALLOCATE (SL%AG(1:SL%NCELLS_LOCAL, 1:NDIAG), STAT=IZERO)
   CALL CHKMEMERR ('SCARC_INITIALIZE_MATRICES3D', 'SL%AG', IZERO)
   SL%AG = 0.0_EB
    
   !!! Allocate local matrix corresponding to the band-wise storage technique
   !ALLOCATE (SL%AL(1:SL%NCELLS_LOCAL, 1:NDIAG), STAT=IZERO)
   !CALL CHKMEMERR ('SCARC_INITIALIZE_MATRICES3D', 'SL%AL', IZERO)
    
   SL%DXI2=1.0_EB/(SL%DX)**2
   SL%DYI2=1.0_EB/(SL%DY)**2
   SL%DZI2=1.0_EB/(SL%DZ)**2

   IF (SCARC_DEBUG >= 2) THEN
      WRITE (9,*) '========= LEVEL=',ILEVEL
      WRITE (9,*) ' IBAR=',SL%IBAR
      WRITE (9,*) ' KBAR=',SL%KBAR
      WRITE (9,*) ' DX  =',SL%DX  
      WRITE (9,*) ' DY  =',SL%DY  
      WRITE (9,*) ' DZ  =',SL%DZ  
      WRITE (9,*) ' DXI2=',SL%DXI2
      WRITE (9,*) ' DYI2=',SL%DYI2
      WRITE (9,*) ' DZI2=',SL%DZI2
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
               WRITE (9,*) 'D1: AG(', IC, ',1)=', SL%AG(IC, 1)
               WRITE (9,*) 'D2: AG(', IC, ',2)=', SL%AG(IC, 2)
               WRITE (9,*) 'D3: AG(', IC, ',3)=', SL%AG(IC, 3)
               WRITE (9,*) 'D4: AG(', IC, ',4)=', SL%AG(IC, 4)
               WRITE (9,*) 'D5: AG(', IC, ',5)=', SL%AG(IC, 5)
               WRITE (9,*) 'D6: AG(', IC, ',6)=', SL%AG(IC, 6)
               WRITE (9,*) 'D7: AG(', IC, ',7)=', SL%AG(IC, 7)
               WRITE (9,*) 
            ENDIF
    
         ENDDO
      ENDDO
   ENDDO
    
   !!! local matrix AL is momentarily set to global matrix AG - not yet needed at all
   !SL%AL = SL%AG
    
   IF (ILEVEL/=S%NLEVEL_MAX) CYCLE

   IF (SCARC_DEBUG>=2) THEN
      WRITE(9,*) 'BOUNDARY_TYPE'
      WRITE(9,'(16i4)') (M%BOUNDARY_TYPE(ICELL),ICELL=1,M%NEWC)
      WRITE(9,*) 'IJKW'
      DO IW=1,M%NEWC+4
         WRITE(9,'(i4,a,15i4)') IW,' : ',(M%IJKW(ICELL,IW),ICELL=1,15)
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

      if (Scarc_debug.ge.8) then
         WRITE(9,*) 'ICELL=',ICELL,' :I=',I,': K=',K,': IW=',IW,': BT=',M%BOUNDARY_TYPE(IW)
      endif


      !IF (M%BOUNDARY_TYPE(IW)==OPEN_BOUNDARY.OR.(M%BOUNDARY_TYPE(IW)==NULL_BOUNDARY.AND.NOM.NE.0)) THEN
      IF (M%BOUNDARY_TYPE(IW)==OPEN_BOUNDARY) THEN
          M%PRESSURE_BC_INDEX(IW)=DIRICHLET
      !   WRITE(9,*) 'D: PRESSURE_BC_INDEX(',IW,')=',M%PRESSURE_BC_INDEX(IW)
      ELSE IF (M%IJKW(9,IW)/=0) THEN
          M%PRESSURE_BC_INDEX(IW)=INTERNAL
      !   WRITE(9,*) 'I: PRESSURE_BC_INDEX(',IW,')=',M%PRESSURE_BC_INDEX(IW)
      ELSE IF (M%BOUNDARY_TYPE(IW)==NULL_BOUNDARY) THEN
          M%PRESSURE_BC_INDEX(IW)=DIRICHLET
      !   WRITE(9,*) 'D: PRESSURE_BC_INDEX(',IW,')=',M%PRESSURE_BC_INDEX(IW)
      ELSE
          M%PRESSURE_BC_INDEX(IW)=NEUMANN
      !   WRITE(9,*) 'N: PRESSURE_BC_INDEX(',IW,')=',M%PRESSURE_BC_INDEX(IW)
      ENDIF

      !
      ! set Dirichlet BC's for AG and AL in open boundary cells
      !
      !IF (M%BOUNDARY_TYPE(IW) == OPEN_BOUNDARY) THEN    
      IF (M%PRESSURE_BC_INDEX(IW) == DIRICHLET) THEN

         SL%AG(IC, 1) = SL%AG(IC, 1) - DBC
         !SL%AL(IC, 1) = SL%AL(IC, 1) - DBC

         IF (SCARC_DEBUG .GE. 4) &
           write (9,'(a,i3,a,i2,a,i3,a,f12.6,a,i3,a,i3)') 'IW=',IW, &
              ': Dirichlet ',IOR,':AG(',IC,',1)=',SL%AG(IC,1),': BT=',M%BOUNDARY_TYPE(IW),': NOM=',NOM

!
      !
      ! set Dirichlet BC's for AL (NOT for AG!) along internal boundaries
      !
      !ELSE IF (NOM /= 0)  THEN    
      ELSE IF (M%PRESSURE_BC_INDEX(IW) == INTERNAL) THEN

         !SL%AL(IC, 1) = SL%AL(IC, 1) - DBC
         IF (SCARC_DEBUG .GE. 4) &
           write (9,'(a,i3,a,i2,a,i3,a,f12.6,a,i3,a,i3)') 'IW=',IW, &
              ': Nothing   ',IOR,':AG(',IC,',1)=',SL%AG(IC,1),': BT=',M%BOUNDARY_TYPE(IW),': NOM=',NOM


      !
      ! set Neumann BC's for AG and AL at all other nodes
      !
      !ELSE 
      ELSE IF (M%PRESSURE_BC_INDEX(IW) == NEUMANN) THEN

         SL%AG(IC, 1) = SL%AG(IC, 1) + DBC
         !SL%AL(IC, 1) = SL%AL(IC, 1) + DBC

         IF (SCARC_DEBUG .GE. 4) &
           write (9,'(a,i3,a,i2,a,i3,a,f12.6,a,i3,a,i3)') 'IW=',IW, &
              ': Neumann   ',IOR,':AG(',IC,',1)=',SL%AG(IC,1),': BT=',M%BOUNDARY_TYPE(IW),': NOM=',NOM

      ENDIF

   ENDDO

   IF (SCARC_DEBUG >= 2) THEN
      WRITE (9,*) 'ASUB=', SL%ASUB
      WRITE (9,*)
      WRITE (9,*) 'GLOBAL MATRIX AG:'
      DO IC = 1, SL%NCELLS_LOCAL
         WRITE (9, '(5f12.4)') SL%AG(IC, 2),SL%AG(IC,3),SL%AG(IC,1),SL%AG(IC,4),SL%AG(IC,5)
         IF (Mod(IC, SL%IBAR) == 0) write (9,*) '---------------------------------------------&
                                                 &--------------------'
      ENDDO
   ENDIF
    
   SL%ASUBX = SL%DXI2
   SL%ASUBY = SL%DYI2
   SL%ASUBZ = SL%DZI2
 
   IF (SCARC_DEBUG >= 2) THEN
     write (9,*) 'SL%ASUBX=', SL%ASUBX
     WRITE (9,*) 'LBC=', LBC
     WRITE (9,*) 'MBC=', MBC
     WRITE (9,*) 'NBC=', NBC
   ENDIF

ENDDO MATRIX_LEVEL_LOOP3D
 
TUSED_SCARC(4,NM)=TUSED_SCARC(4,NM)+SECOND()-TNOW_MATRICES3D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_MATRICES3D

END SUBROUTINE SCARC_INITIALIZE_MATRICES3D


 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Initialize global 2D-solver methods (cg/mg)
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_INITIALIZE_SOLVER2D (NM)
 
INTEGER :: NM, IZERO, ILEVEL
INTEGER :: IBP0, KBP0
REAL(EB):: TNOW_SOLVER2D
 
TNOW_SOLVER2D = SECOND()
 
IF (SCARC_DEBUG >= 2) THEN
   write (9,*) '=========== STARTING SCARC_INITIALIZE_SOLVER'
   flush (9)
ENDIF

IZERO = 0

!!! CG-parameters for global and local solves
S%BPRECON   = .TRUE.
S%BGSADI    = .FALSE.


S => SCARC (NM)

!!!
!!! Initialize SLEVEL-structures for single levels
!!!

LEVEL_LOOP2D: DO ILEVEL=S%NLEVEL_MAX,S%NLEVEL_MIN,-1

IF (SCARC_DEBUG>=2) WRITE(9,*) 'ILEVEL=',ILEVEL

   SL => S%SLEVEL(ILEVEL)

IF (SCARC_DEBUG>=2) WRITE(9,*)' ALLOCATE SCARC(',NM,')%SLEVEL(',ILEVEL,')'

   IBP0=SL%IBAR+1
   KBP0=SL%KBAR+1

   ! ----------------------------------------------------------------------------
   ! working and auxiliary vectors for MG-method
   ! ----------------------------------------------------------------------------
   ALLOCATE (SL%X(0:IBP0, 1, 0:KBP0), STAT=IZERO)
   CALL CHKMEMERR ('SCARC', 'X', IZERO)
   SL%X = 0.0_EB
    
   ALLOCATE (SL%F(0:IBP0, 1, 0:KBP0), STAT=IZERO)
   CALL CHKMEMERR ('SCARC', 'B', IZERO)
   SL%F = 0.0_EB
    
   ALLOCATE (SL%D(0:IBP0, 1, 0:KBP0), STAT=IZERO)
   CALL CHKMEMERR ('SCARC', 'D', IZERO)
   SL%D = 0.0_EB
    
   ALLOCATE (SL%R(0:IBP0, 1, 0:KBP0), STAT=IZERO)
   CALL CHKMEMERR ('SCARC', 'R', IZERO)
   SL%R = 0.0_EB
    
   ALLOCATE (SL%G(0:IBP0, 1, 0:KBP0), STAT=IZERO)
   CALL CHKMEMERR ('SCARC', 'G', IZERO)
   SL%G = 0.0_EB
    
   ALLOCATE (SL%Y(0:IBP0, 1, 0:KBP0), STAT=IZERO)
   CALL CHKMEMERR ('SCARC', 'Y', IZERO)
   SL%Y = 0.0_EB
    
   ALLOCATE (SL%Z(0:IBP0, 1, 0:KBP0), STAT=IZERO)
   CALL CHKMEMERR ('SCARC', 'Z', IZERO)
   SL%Z = 0.0_EB
    

   ! ----------------------------------------------------------------------------
   ! vectors for the description of the boundary conditions
   ! ----------------------------------------------------------------------------
   ALLOCATE (SL%BXS0(1, KBP0), STAT=IZERO)
   CALL CHKMEMERR ('SCARC', 'BXS0', IZERO)
   SL%BXS0 = 0.0_EB
  
   ALLOCATE (SL%BXF0(1, KBP0), STAT=IZERO)
   CALL CHKMEMERR ('SCARC', 'BXF0', IZERO)
   SL%BXF0 = 0.0_EB

   ALLOCATE (SL%BZS0(IBP0, 1), STAT=IZERO)
   CALL CHKMEMERR ('SCARC', 'BZS0', IZERO)
   SL%BZS0 = 0.0_EB

   ALLOCATE (SL%BZF0(IBP0, 1), STAT=IZERO)
   CALL CHKMEMERR ('SCARC', 'BZF0', IZERO)
   SL%BZF0 = 0.0_EB

ENDDO LEVEL_LOOP2D

    
!!!
!!! Allocate vectors for global scalar products 
!!!
ALLOCATE (S%SP_LOCAL0(1:NMESHES), STAT=IZERO)
CALL CHKMEMERR ('SCARC', 'SP_LOCAL0', IZERO)
S%SP_LOCAL0 = 0.0_EB
 
ALLOCATE (S%SP_GLOBAL0(1:NMESHES), STAT=IZERO)
CALL CHKMEMERR ('SCARC', 'SP_GLOBAL0', IZERO)
S%SP_GLOBAL0 = 0.0_EB


!!!
!!! Allocate array for the description of the MG-cycle
!!!
IF (SCARC_METHOD==2) THEN

   ALLOCATE (S%KCYCLE(2, NNLEVEL), STAT=IZERO)
   CALL CHKMEMERR ('SCARC', 'KCYCLE', IZERO)
   S%KCYCLE = 0
    
   ALLOCATE (S%KSMOOTH_PRE(0:NNLEVEL), STAT=IZERO)
   CALL CHKMEMERR ('SCARC', 'KSMOOTH_PRE', IZERO)
   S%KSMOOTH_PRE = 0
    
   ALLOCATE (S%KSMOOTH_POST(0:NNLEVEL), STAT=IZERO)
   CALL CHKMEMERR ('SCARC', 'KSMOOTH_POST', IZERO)
   S%KSMOOTH_POST = 0

ENDIF
    
IF (SCARC_DEBUG >= 2) write (9,*) 'INITIALIZE_SOLVER2D: FINISHED'

TUSED_SCARC(5,NM)=TUSED_SCARC(5,NM)+SECOND()-TNOW_SOLVER2D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_SOLVER2D
 
END SUBROUTINE SCARC_INITIALIZE_SOLVER2D

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Initialize global 3D-solver methods (cg/mg)
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_INITIALIZE_SOLVER3D (NM)
 
INTEGER :: NM, IZERO, ILEVEL
INTEGER :: IBP0, JBP0, KBP0
REAL(EB):: TNOW_SOLVER3D
 
TNOW_SOLVER3D = SECOND()
 
IF (SCARC_DEBUG >= 2) write (9,*) '=========== STARTING SCARC_INITIALIZE_SOLVER'


!!! CG-parameters for global and local solves
S%BPRECON   = .TRUE.
S%BGSADI    = .FALSE.

IZERO = 0

S => SCARC (NM)

!!!
!!! Initialize SLEVEL-structures for single levels
!!!
LEVEL_LOOP3D: DO ILEVEL=S%NLEVEL_MAX,S%NLEVEL_MIN,-1

   SL => S%SLEVEL(ILEVEL)

   IBP0=SL%IBAR+1
   JBP0=SL%JBAR+1
   KBP0=SL%KBAR+1

   ! ----------------------------------------------------------------------------
   ! working and auxiliary vectors for MG-method
   ! ----------------------------------------------------------------------------
   ALLOCATE (SL%X(0:IBP0, 0:JBP0, 0:KBP0), STAT=IZERO)
   CALL CHKMEMERR ('SCARC', 'X', IZERO)
   SL%X = 0.0_EB
    
   ALLOCATE (SL%F(0:IBP0, 0:JBP0, 0:KBP0), STAT=IZERO)
   CALL CHKMEMERR ('SCARC', 'B', IZERO)
   SL%F = 0.0_EB
    
   ALLOCATE (SL%D(0:IBP0, 0:JBP0, 0:KBP0), STAT=IZERO)
   CALL CHKMEMERR ('SCARC', 'D', IZERO)
   SL%D = 0.0_EB
    
   ALLOCATE (SL%R(0:IBP0, 0:JBP0, 0:KBP0), STAT=IZERO)
   CALL CHKMEMERR ('SCARC', 'R', IZERO)
   SL%R = 0.0_EB
    
   ALLOCATE (SL%G(0:IBP0, 0:JBP0, 0:KBP0), STAT=IZERO)
   CALL CHKMEMERR ('SCARC', 'G', IZERO)
   SL%G = 0.0_EB
    
   ALLOCATE (SL%Y(0:IBP0, 0:JBP0, 0:KBP0), STAT=IZERO)
   CALL CHKMEMERR ('SCARC', 'Y', IZERO)
   SL%Y = 0.0_EB
    
   ALLOCATE (SL%Z(0:IBP0, 0:JBP1, 0:KBP0), STAT=IZERO)
   CALL CHKMEMERR ('SCARC', 'Z', IZERO)
   SL%Z = 0.0_EB
    

   ! ----------------------------------------------------------------------------
   ! vectors for the description of the boundary conditions
   ! ----------------------------------------------------------------------------
   ALLOCATE (SL%BXS0(JBP0, KBP0), STAT=IZERO)
   CALL CHKMEMERR ('SCARC', 'BXS0', IZERO)
   SL%BXS0 = 0.0_EB
    
   ALLOCATE (SL%BXF0(JBP0, KBP0), STAT=IZERO)
   CALL CHKMEMERR ('SCARC', 'BXF0', IZERO)
   SL%BXF0 = 0.0_EB
    
   ALLOCATE (SL%BYS0(IBP0, KBP0), STAT=IZERO)
   CALL CHKMEMERR ('SCARC', 'BYS0', IZERO)
   SL%BYS0 = 0.0_EB
    
   ALLOCATE (SL%BYF0(IBP0, KBP0), STAT=IZERO)
   CALL CHKMEMERR ('SCARC', 'BYF0', IZERO)
   SL%BYF0 = 0.0_EB
    
   ALLOCATE (SL%BZS0(IBP0, JBP0), STAT=IZERO)
   CALL CHKMEMERR ('SCARC', 'BZS0', IZERO)
   SL%BZS0 = 0.0_EB
    
   ALLOCATE (SL%BZF0(IBP0, JBP0), STAT=IZERO)
   CALL CHKMEMERR ('SCARC', 'BZF0', IZERO)
   SL%BZF0 = 0.0_EB
   
ENDDO LEVEL_LOOP3D

    
!!!
!!! Allocate vectors for global scalar products 
!!!
ALLOCATE (S%SP_LOCAL0(1:NMESHES), STAT=IZERO)
CALL CHKMEMERR ('SCARC', 'SP_LOCAL0', IZERO)
S%SP_LOCAL0 = 0.0_EB
 
ALLOCATE (S%SP_GLOBAL0(1:NMESHES), STAT=IZERO)
CALL CHKMEMERR ('SCARC', 'SP_GLOBAL0', IZERO)
S%SP_GLOBAL0 = 0.0_EB


!!!
!!! Allocate array for the description of the MG-cycle
!!!
IF (SCARC_METHOD==2) THEN

   ALLOCATE (S%KCYCLE(2, NNLEVEL), STAT=IZERO)
   CALL CHKMEMERR ('SCARC', 'KCYCLE', IZERO)
   S%KCYCLE = 0
    
   ALLOCATE (S%KSMOOTH_PRE(0:NNLEVEL), STAT=IZERO)
   CALL CHKMEMERR ('SCARC', 'KSMOOTH_PRE', IZERO)
   S%KSMOOTH_PRE = 0
    
   ALLOCATE (S%KSMOOTH_POST(0:NNLEVEL), STAT=IZERO)
   CALL CHKMEMERR ('SCARC', 'KSMOOTH_POST', IZERO)
   S%KSMOOTH_POST = 0

ENDIF
    
IF (SCARC_DEBUG >= 2) write (9,*) 'INITIALIZE_SOLVER: FINISHED'
 
TUSED_SCARC(5,NM)=TUSED_SCARC(5,NM)+SECOND()-TNOW_SOLVER3D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_SOLVER3D

END SUBROUTINE SCARC_INITIALIZE_SOLVER3D
 

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Perform global cg-method based on global possion-matrix
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_CG2D (NM)
 
INTEGER ::  NM, ILMAX, ITE, IREL, ICONV, ITYPE, IERR, IPRECON, I, K
REAL (EB) :: SIGMA0, SIGMA1, ALPHA, GAMMA
REAL (EB) :: TNOW_CG2D
TYPE (MESH_TYPE),  POINTER :: M
 
INCLUDE 'mpif.h'
 
TNOW_CG2D = SECOND()
IERR = 0

S  => SCARC (NM)
M  => MESHES (NM)

ILMAX = S%NLEVEL_MAX
SL => S%SLEVEL(ILMAX)
 
! initialize working vectors
SL%G = 0.0_EB
SL%Y = 0.0_EB
SL%R = 0.0_EB
SL%D = 0.0_EB
SL%X = 0.0_EB
 
DO K = 1, KBAR
   DO I = 1, IBAR
      SL%F (I, 1, K) = M%PRHS (I, 1, K)
   ENDDO
ENDDO
 
IF (SCARC_DEBUG >= 3) THEN
   WRITE (9,*) 'Starting SCARC_CG2D'
   WRITE (9,*) 'NMESHES=:', NMESHES
   WRITE (9,*) 'IBAR=:', IBAR
   WRITE (9,*) 'JBAR=:', JBAR
   WRITE (9,*) 'KBAR=:', KBAR
   WRITE (9,*) 'SCARC_NIT=:', SCARC_NIT
   WRITE (9,*) 'SCARC_EPS=:', SCARC_EPS
   WRITE (9,*) 'SCARC_EPS_REL=:', SCARC_EPS_REL
   WRITE (9,*) 'ncells_global=:', SL%Ncells_global
ENDIF
CALL SCARC_COMPARE_SINGLE (SL%X, 'SARC', 1, 'X1    ',0)
CALL SCARC_COMPARE_SINGLE (SL%F, 'SARC', 1, 'F1    ',0)
 
 
!!!
!!! initialize some parameters
!!!
ITYPE = NCOM_TYPE_GLOBAL
IREL = 0
ICONV = 0
IPRECON = 2
 
!!!
!!! set boundary conditions along exterior boundaries
!!!
CALL SCARC_SETBDRY2D (ILMAX, NM)
 
 
!!!
!!! calculate initial defect Y = B - A*X and get l2-norm of it
!!!
CALL SCARC_MATVEC2D (SL%AG, SL%X, SL%R, 1.0_EB, 0.0_EB, NM, ITYPE, ILMAX, 3)
IF (SCARC_DEBUG>=2) WRITE(9,*) 'SCARC_CG2D: FIRST MATVEC2D'
CALL SCARC_VECADD2D (SL%F, SL%R,-1.0_EB, 1.0_EB, ILMAX, NM)
CALL SCARC_L2NORM2D (SL%R, SL%RES, ILMAX, NM, ITYPE)
SL%RESIN = SL%RES
 
IF (SCARC_DEBUG >= 3) THEN
   WRITE (9,*) 'Initial DEFECT '
   WRITE (9,*) 'RES=', SL%RES
   WRITE (9,*) 'SL%X:'
   WRITE (9, '(a,10f12.6)') '2: S%SP_LOCAL0=', (S%SP_LOCAL0(I), I=1, NMESHES)
ENDIF
CALL SCARC_COMPARE_SINGLE (SL%X, 'SARC', 1, 'X3    ',0)
CALL SCARC_COMPARE_SINGLE (SL%F, 'SARC', 1, 'F3    ',0)
CALL SCARC_COMPARE_SINGLE (SL%R, 'SARC', 1, 'R3    ',0)
 
 
IF (S%BPRECON) THEN
 
   IF (SCARC_PRECON == 1) THEN
      WRITE (*,*) 'IPRECON=1 not yet implemented!'
   ELSE IF (SCARC_PRECON == 2) THEN
      CALL SCARC_VECCOPY2D (SL%R, SL%G, ILMAX, NM)
      CALL SCARC_PRECON_JACOBI2D (SL%AG, SL%G, ILMAX, NM)
   ELSE IF (SCARC_PRECON == 3) THEN
      CALL SCARC_VECCOPY2D (SL%R, SL%G, ILMAX, NM)
      CALL SCARC_PRECON_GS2D (SL%AG, SL%G, ILMAX, NM)
   ELSE IF (SCARC_PRECON == 4) THEN
      CALL SCARC_VECCOPY2D (SL%R, SL%G, ILMAX, NM)
      CALL SCARC_PRECON_SSOR2D (SL%AG, SL%G, ILMAX, NM)
   ENDIF
 
   CALL SCARC_SCALPROD2D (SL%R, SL%G, SIGMA0, ILMAX, NM)
 
ELSE
   SIGMA0 = SL%RESIN * SL%RESIN
ENDIF
 
IF (SCARC_DEBUG >= 3) THEN
   WRITE (9,*) 'After initial Precon '
   WRITE (9,*) 'Sigma0=', SIGMA0
   WRITE (9,*) 'SL%X:'
   WRITE (9, '(a,10f12.6)') '4: S%SP_LOCAL0=', (S%SP_LOCAL0(I), I=1, NMESHES)
   CALL SCARC_COMPARE_SINGLE (SL%R, 'SARC', 1, 'R4    ',0)
   CALL SCARC_COMPARE_SINGLE (SL%G, 'SARC', 1, 'G4    ',0)
   CALL flush (9)
ENDIF
 
CALL SCARC_VECADD2D (SL%G, SL%D,-1.0_EB, 0.0_EB, ILMAX, NM)
 
IF (SCARC_DEBUG >= 3) THEN
   WRITE (9,*) 'After initial vecadd SL%D'
   WRITE (9, '(a,10f12.6)') '5: S%SP_LOCAL0=', (S%SP_LOCAL0(I), I=1, &
    NMESHES)
   CALL flush (9)
ENDIF
CALL SCARC_COMPARE_SINGLE (SL%D, 'SARC', 1, 'D4    ',0)
 
 
!!!
!!! start defect correction loop
!!!
CG_LOOP2D: DO ITE = 1, SCARC_NIT
 
   IF (SCARC_DEBUG >= 3) WRITE (9,*) '=========================== STARTING CG-LOOP ', ITE

 
!!!  calculate new defect and get L2-norm of it
   CALL SCARC_MATVEC2D (SL%AG, SL%D, SL%Y, 1.0_EB, 0.0_EB, NM, ITYPE, ILMAX, 5)
   CALL SCARC_SCALPROD2D (SL%D, SL%Y, ALPHA, ILMAX, NM)

   ALPHA = SIGMA0 / ALPHA
 
   CALL SCARC_VECADD2D (SL%D, SL%X, ALPHA, 1.0_EB, ILMAX, NM)
   CALL SCARC_VECADD2D (SL%Y, SL%R, ALPHA, 1.0_EB, ILMAX, NM)
   CALL SCARC_L2NORM2D (SL%R, SL%RES, ILMAX, NM, ITYPE)

   IF (SCARC_DEBUG >= 2) THEN
      CALL SCARC_COMPARE_SINGLE (SL%Y, 'SARC', 1, 'DNEW   ',0)
      CALL SCARC_COMPARE_SINGLE (SL%R, 'SARC', 1, 'YNEW   ',0)
      CALL SCARC_COMPARE_SINGLE (SL%R, 'SARC', 1, 'XNEW   ',0)
      CALL SCARC_COMPARE_SINGLE (SL%R, 'SARC', 1, 'RNEW   ',0)
   ENDIF
 
   IF (SCARC_DEBUG >= 3) THEN
      WRITE (9,*) 'New DEFECT in iteration ', ITE, ': ALPHA=', ALPHA
      WRITE (9,*) 'New Iterate in iteration ', ITE
      WRITE (9, '(a10,f22.16)') 'ALPHA=', ALPHA
      WRITE (9, '(a10,f22.16)') 'SIGMA0=', SIGMA0
      WRITE (9, '(a10,f22.16)') 'RES=', SL%RES
      WRITE (9, '(a,10f12.6)') '6: S%SP_LOCAL0=', (S%SP_LOCAL0(I), I=1, NMESHES)
      CALL SCARC_COMPARE_SINGLE (SL%Y, 'SARC', 1, 'Y5    ',0)
      CALL SCARC_COMPARE_SINGLE (SL%R, 'SARC', 1, 'R5    ',0)
      CALL flush (9)
   ENDIF
   CALL SCARC_COMPARE_SINGLE (SL%X, 'SARC', 1, 'XITE   ',0)
 
 
   !if (MOD(ITE,1)==0) write(*,'(a,i4,a,f25.16,a,f25.16)') '### SCARC_CG, ite=',ITE,' RES=',SL%RES,'  ####'
   !IF (NM == 1 .AND. SCARC_DEBUG >= 1) write (9, '(a,i4,a,f25.16,a,f25.16)') '&
   IF (SCARC_DEBUG >= 1) write (9, '(a,i4,a,f25.16,a,f25.16)') '### SCARC_CG, ite=', ITE, ' RES=', SL%RES, '  ####'
 
!!!  stop in case of convergence or divergence
   IF (IREL == 1) THEN
      IF (SL%RES <= SL%RESIN*SCARC_EPS) ICONV = 1
   ELSE
      IF (SL%RES <= SCARC_EPS .AND. SL%RES <= SL%RESIN*SCARC_EPS_REL) ICONV = 1
   ENDIF
   IF (SL%RES > 10000.0_EB) ICONV = - 1
   IF (Abs(ICONV) == 1) EXIT CG_LOOP2D
 
 
!!! preconditioning
   IF (S%BPRECON) THEN
 
      IF (SCARC_PRECON == 1) THEN
         CALL SCARC_RECEIVE (NCOM_FULL, ILMAX)   !Aufruf 2,0
         CALL SCARC_EXCHANGE (NCOM_FULL, ILMAX, NTYPE_NONE)
         WRITE (*,*) 'IPRECON=1 not yet implemented!'
!     CALL SCARC_CG0(SL%R,NM)
      ELSE IF (SCARC_PRECON == 2) THEN
         CALL SCARC_VECCOPY2D (SL%R, SL%G, ILMAX, NM)
         CALL SCARC_PRECON_JACOBI2D (SL%AG, SL%G, ILMAX, NM)
      ELSE IF (SCARC_PRECON == 3) THEN
         CALL SCARC_VECCOPY2D (SL%R, SL%G, ILMAX, NM)
         CALL SCARC_PRECON_GS2D (SL%AG, SL%G, ILMAX, NM)
      ELSE IF (SCARC_PRECON == 4) THEN
         CALL SCARC_VECCOPY2D (SL%R, SL%G, ILMAX, NM)
         CALL SCARC_PRECON_SSOR2D (SL%AG, SL%G, ILMAX, NM)
      ENDIF
 
      CALL SCARC_SCALPROD2D (SL%R, SL%G, SIGMA1, ILMAX, NM)
 
   ELSE
      SIGMA1 = SL%RES * SL%RES
   ENDIF
 
   CALL SCARC_COMPARE_SINGLE (SL%G, 'SARC', 1, 'G6    ',0)
 
   GAMMA = SIGMA1 / SIGMA0
   SIGMA0 = SIGMA1
 
   IF (SCARC_DEBUG >= 3) THEN
      WRITE (9,*) 'After Precon in iteration ', ITE
      WRITE (9,*) 'Sigma0=', SIGMA0
      WRITE (9,*) 'Gamma =', GAMMA
      CALL SCARC_COMPARE_SINGLE (SL%R, 'SARC', 1, 'R4    ',0)
      CALL flush (9)
   ENDIF
 
   CALL SCARC_VECADD2D (SL%G, SL%D,-1.0_EB, GAMMA, ILMAX, NM)
 
   IF (SCARC_DEBUG >= 3) CALL SCARC_COMPARE_SINGLE (SL%D, 'SARC', 1, 'D7    ',0)
 
ENDDO CG_LOOP2D
 
 
IF (SCARC_DEBUG >= 2) THEN
   write (9,*) 'AFTER CG_LOOP2D:A'
   write (9,*) 'RESIN=',SL%RESIN
   write (9,*) 'RES  =',SL%RES
   write (9,*) 'CAPPA  =',SL%CAPPA
ENDIF
 
!!! 'bad end'
IF (ICONV ==-1) THEN
   ITE = - 1
   SL%CAPPA = 1.0_EB
 
!!! 'good end'
ELSE
   IF (SL%RESIN >= 1.0E-70_EB) THEN
      SL%RESIN = SL%RES / SL%RESIN
      SL%CAPPA = SL%RESIN ** (1.0_EB/ITE)
   ELSE IF (ITE == 0) THEN
      SL%CAPPA = 0.0_EB
   ENDIF
ENDIF
 
IF (SCARC_DEBUG >= 2) THEN
   write (9,*) 'AFTER CG_LOOP2D'
   write (9,*) 'RESIN=',SL%RESIN
   write (9,*) 'RES  =',SL%RES
   write (9,*) 'CAPPA  =',SL%CAPPA
ENDIF
 
IF (SCARC_DEBUG >= 1) THEN
   WRITE (9, '(a,i4,a,f12.6,a,f12.5)') 'SCARC_CG: Ite=', ITE, ': Res=', SL%RES,' : Kappa=', SL%CAPPA 
ENDIF
IF (NM == 1) write (*, '(a20,i4,a,e12.6,a,e12.6)') 'SCARC_CG: Ite=', ITE, &
             ': Res =', SL%RES, ': Kappa=', SL%CAPPA
 
 
DO K = 1, KBAR
   DO I = 1, IBAR
      M%PRHS (I, 1, K) = SL%X(I, 1, K)
   ENDDO
ENDDO
 
IF (SCARC_DEBUG >= 2) THEN
   WRITE (9,*) 'Leaving SCARC_CG2D'
   CALL flush (9)
ENDIF
CALL SCARC_COMPARE_SINGLE (SL%X, 'SARC', 1, 'XFINAL',0)
 
TUSED_SCARC(6,NM)=TUSED_SCARC(6,NM)+SECOND()-TNOW_CG2D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_CG2D

END SUBROUTINE SCARC_CG2D



 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Perform global MG-method based on global possion-matrix
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_MG2D (NM)
 
INTEGER :: NM, ITE, IREL, ICONV, ITYPE, ICYCLE
INTEGER :: IERR, IPRECON, ILEVEL, ILMAX
INTEGER :: I, K, IBAR0, JBAR0, KBAR0
REAL (EB) :: TNOW_MG2D
LOGICAL BMATVEC, BCONV, BRELERR
!TYPE (SCARC_LEVEL_TYPE), POINTER :: SL, SLHI, SLLO
TYPE (MESH_TYPE),  POINTER :: M
 
INCLUDE 'mpif.h'
 
TNOW_MG2D = SECOND()

 
! initialize working vectors
S    => SCARC (NM)
M    => MESHES (NM)

ILMAX =  S%NLEVEL_MAX
SLHI  => S%SLEVEL(ILMAX)
 

! initialize working vectors
SLHI%G = 0.0_EB
SLHI%Y = 0.0_EB
SLHI%R = 0.0_EB
SLHI%D = 0.0_EB
SLHI%X = 0.0_EB
 
DO K = 1, KBAR
   DO I = 1, IBAR
      SLHI%F (I, 1, K) = M%PRHS (I, 1, K)
   ENDDO
ENDDO
 
IF (SCARC_DEBUG >= 3) THEN
   WRITE (9,*) 'Starting SCARC_MG2D'
   WRITE (9,*) 'IBAR=:', SLHI%IBAR
   WRITE (9,*) 'JBAR=:', SLHI%JBAR
   WRITE (9,*) 'KBAR=:', SLHI%KBAR
   WRITE (9,*) 'SCARC_NIT=:', SCARC_NIT
   WRITE (9,*) 'SCARC_EPS=:', SCARC_EPS
   WRITE (9,*) 'SCARC_EPS_REL=:', SCARC_EPS_REL
   CALL SCARC_COMPARE_SINGLE (SLHI%X, 'SARC', 1, 'X1    ',0)
   CALL SCARC_COMPARE_SINGLE (SLHI%F, 'SARC', 1, 'F1    ',0)
ENDIF
 
!!!
!!! initialize some method parameters
!!!
ITYPE = NCOM_TYPE_GLOBAL
IREL = 0
ICONV = 0
IPRECON = 2
 
NREQ_FACE = 0
NREQ_EDGE = 0
NREQ_DIAG = 0

IERR = 0
BCONV  =.FALSE.
BRELERR=.TRUE.
ICYCLE=1      ! V-cycle
 
S%EPS_CONV=1.0E-10_EB
S%EPS_DIVG=1.0E+6_EB
S%EPS_REL=1.0E-2_EB
 
IBAR0=SLHI%IBAR
JBAR0=SLHI%JBAR
KBAR0=SLHI%KBAR

 
!!! adjust boundary values of right hand side
CALL SCARC_SETBDRY2D (ILMAX, NM)
 

!!! ------------------------------------------------------------------------------------
!!! Only one level: solve problem exactly
!!! ------------------------------------------------------------------------------------
IF (S%NLEVEL==1) THEN

      CALL SCARC_COARSE2D(ILMAX, NM)

!!! ------------------------------------------------------------------------------------
!!! More than one level: start MG-cycling
!!! ------------------------------------------------------------------------------------
ELSE

   !!!
   !!! save cycle counts for MG-iteration
   !!!
   S%KCYCLE(2,S%NLEVEL_MAX)=1
   DO ILEVEL=S%NLEVEL_MIN+1,S%NLEVEL_MAX-1
      IF (ICYCLE==0) THEN
         S%KCYCLE(2,ILEVEL)=2
      ELSE
         S%KCYCLE(2,ILEVEL)=ICYCLE
      ENDIF
   ENDDO
     
    
   !!!
   !!! calculate initial defect on finest level, Y = B - A*X, and get l2-norm of it
   !!!
!WRITE(9,*) 'KOMMUNIKATIONSPARAMETER 3 prüfen, anders als bei CG !!!'
!write(*,*) 'KOMMUNIKATIONSPARAMETER 3 prüfen, anders als bei CG !!!'
   ILEVEL=ILMAX
   CALL SCARC_MATVEC2D (SLHI%AG, SLHI%X, SLHI%D, 1.0_EB, 0.0_EB, NM, ITYPE, ILEVEL, 3)

   CALL SCARC_COMPARE_SINGLE (SLHI%D, 'SARC', 1, 'D1    ',0)

   CALL SCARC_VECADD2D (SLHI%F, SLHI%D,1.0_EB, -1.0_EB, ILEVEL, NM)
   CALL SCARC_L2NORM2D (SLHI%D, SLHI%RES, ILEVEL, NM, ITYPE)

   S%ITE   = 0
   S%RESIN = SLHI%RES
    
   !WRITE(*,*) 'Info muss noch an Master geschickt werden'
   !WRITE(9,*) 'Info muss noch an Master geschickt werden'

   !CALL SCARC_MASTER_INFO (SLHI%ITE, SLHI%RESIN)

   !!!
   !!! Debug values if requested
   !!!
   IF (SCARC_DEBUG >= 3) THEN
      WRITE (9,*) 'RESIN=', SLHI%RESIN
      CALL flush (9)
   ENDIF
   CALL SCARC_COMPARE_SINGLE (SLHI%X, 'SARC', 1, 'X3    ',0)
   CALL SCARC_COMPARE_SINGLE (SLHI%F, 'SARC', 1, 'F3    ',0)
   CALL SCARC_COMPARE_SINGLE (SLHI%R, 'SARC', 1, 'R3    ',0)
    
   IF (SLHI%RES<S%EPS_CONV) THEN
      WRITE(*,*) 'Good end - muss noch gestoppt werden'
      WRITE(9,*) 'Good end - muss noch gestoppt werden'
      goto 200
   ENDIF
    
   IF (SLHI%RES>S%EPS_DIVG) THEN
      WRITE(*,*) 'Bad end - muss noch gestoppt werden'
      WRITE(9,*) 'Bad end - muss noch gestoppt werden'
      goto 210
   ENDIF
    
    
   !!!
   !!! start MG-iteration
   !!!
   !WRITE(*,*) 'HALLO HIER BESSEREN WERT FUER SCARC_NIT!!'
   !WRITE(9,*) 'HALLO HIER BESSEREN WERT FUER SCARC_NIT!!'
   MG_LOOP2D: DO ITE = 1, SCARC_NIT
    
      S%ITE=ITE

      IF (SCARC_DEBUG >= 3) THEN
         WRITE (9,*) '=========================== STARTING MG-LOOP ', ITE
      ENDIF

      !!! set level-information
      DO ILEVEL=S%NLEVEL_MIN,S%NLEVEL_MAX
         S%KCYCLE(1,ILEVEL)=S%KCYCLE(2,ILEVEL)
      ENDDO

      ILEVEL=S%NLEVEL_MAX

    
      !!! ---------------------------------------------------------------------------
      !!!  Presmoothing
      !!! ---------------------------------------------------------------------------
110   IF (ILEVEL/=S%NLEVEL_MIN) THEN
    
         SLHI => S%SLEVEL(ILEVEL)
         SLLO => S%SLEVEL(ILEVEL-1)

         !IF (SCARC_METHOD>=2) THEN
         !   WRITE(9,*) 'SLHI zeigt auf level ',ILEVEL
         !   WRITE(9,*) 'SLLO zeigt auf level ',ILEVEL-1
         !ENDIF

         IF (ILEVEL==S%NLEVEL_MAX) THEN
            BMATVEC=.FALSE.
         ELSE
            BMATVEC=.TRUE.
         ENDIF
         BMATVEC=.TRUE.

         CALL SCARC_SMOOTHER2D(SLHI%AG, SLHI%X, SLHI%F, SLHI%D, BMATVEC, ILEVEL, NM)

         !!! select smoother for presmoothing
         !!! print presmoothin information on current level
         !CALL SCARC_INFO_LEVEL(ILEVEL)

         !!! perform restriction to coarser grid: F:= rest(D)
         CALL SCARC_COMPARE_SINGLE_LEV (SLHI%D, 'SARC', 1, 'RESTD ',0,NM,ILEVEL)
         CALL SCARC_RESTRICTION2D(NM, SLLO%F, SLHI%D, SLLO%IBAR, SLLO%KBAR)
         CALL SCARC_COMPARE_SINGLE_LEV (SLLO%F, 'SARC', 1, 'RESTF ',0,NM,ILEVEL-1)

         !!! initialize solution vector and set boundary conditions to residuum
         SLLO%X=0.0_EB
         !CALL SCARC_BDRY_RESIDUUM(NM, SLLO%F, ILEVEL-1)

         !!! decrease level
         ILEVEL=ILEVEL-1

         GOTO 110
      
      ENDIF         

      !!! ------------------------------------------------------------------------
      !!! Coarse grid solver
      !!! ------------------------------------------------------------------------
      ILEVEL=S%NLEVEL_MIN
      CALL SCARC_COARSE2D(ILEVEL, NM)



      !!! ------------------------------------------------------------------------
      !!! Postsmoothing
      !!! ------------------------------------------------------------------------
120   IF (ILEVEL/=S%NLEVEL_MAX) THEN
    
         BMATVEC=.TRUE.

         SLLO => S%SLEVEL(ILEVEL)

         ILEVEL=ILEVEL+1
         SLHI => S%SLEVEL(ILEVEL)

         !!! print information
         !CALL SCARC_INFO_LEVEL(ILEVEL)

         !!! perform prolongation to finer grid
         CALL SCARC_PROLONGATION2D(NM, SLLO%X, SLHI%D, ILEVEL)

         !!! set exterior boundary data of residuum to zero
         !CALL SCARC_BDRY_RESIDUUM(NM, SLHI%D, ILEVEL)
 
         !!! set new solution
         CALL SCARC_VECADD2D(SLHI%D, SLHI%X, 1.0_EB, 1.0_EB, ILEVEL, NM)

         !!! select smoother for postsmoothing
         CALL SCARC_SMOOTHER2D(SLHI%AG, SLHI%X, SLHI%F, SLHI%D, BMATVEC, ILEVEL, NM)


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
      CALL SCARC_MATVEC2D (SLHI%AG, SLHI%X, SLHI%D, 1.0_EB, 0.0_EB, NM, ITYPE, ILEVEL, 3)
      CALL SCARC_VECADD2D (SLHI%F, SLHI%D,1.0_EB, -1.0_EB, ILEVEL, NM)
      CALL SCARC_L2NORM2D (SLHI%D, S%RES, ILEVEL, NM, ITYPE)

      IF (SCARC_DEBUG>=2) WRITE(9,*) 'SCARC-Multigrid, iteration ',S%ITE,': residuum=',S%RES
      !IF (NM==1) write(*,*) 'SCARC-MG2D      :  Iteration ',S%ITE,': Residuum=',S%RES

      !!! ------------------------------------------------------------------------
      !!! Send error information to Master
      !!! ------------------------------------------------------------------------
      !CALL SCARC_INFO_MASTER
      !CALL SCARC_INFO_LEVEL(ILEVEL)


      !!! ------------------------------------------------------------------------
      !!! Convergence or divergence ?
      !!! ------------------------------------------------------------------------
      IF (SL%RES>S%EPS_DIVG) THEN
         WRITE(*,*) 'Bad end - muss noch gestoppt werden'
         WRITE(9,*) 'Bad end - muss noch gestoppt werden'
      ENDIF

      IF (SCARC_DEBUG >=2) THEN
         WRITE(9,*) 'S%EPS_CONV=',S%EPS_CONV
         WRITE(9,*) 'S%RES=',S%RES
         WRITE(9,*) 'S%RESIN=',S%RESIN
      ENDIF

      IF (S%BREL_ERR) THEN
         IF (S%RES <= S%RESIN*S%EPS_CONV) BCONV=.TRUE.
      ELSE
         IF (S%RES <= S%EPS_CONV .AND. S%RES <= S%RESIN*S%EPS_CONV) BCONV=.TRUE.
      ENDIF
      IF (BCONV) GOTO 200
    
   ENDDO MG_LOOP2D


   !!! ------------------------------------------------------------------------
   !!! Convergence or divergence ?
   !!! ------------------------------------------------------------------------
200 CONTINUE
   IF (SCARC_DEBUG >=2) THEN
      WRITE(9,*) 'BEI 200'
      WRITE(9,*) 'S%RES  =',S%RES
      WRITE(9,*) 'S%RESIN=',S%RESIN
      WRITE(9,*) 'S%ITE  =',S%ITE
   ENDIF
   IF (S%RESIN>1.0E-70_EB) S%RESIN=S%RES/S%RESIN
   IF (S%ITE==0) THEN
      S%CAPPA=0.0_EB
   ELSE
      S%CAPPA=S%RESIN**(1.0_EB/REAL(S%ITE,EB))
   ENDIF
   IF (SCARC_DEBUG>=2)    WRITE(9,*) 'S%CAPPA=',S%CAPPA

   GOTO 250
    
   !!! ------------------------------------------------------------------------
   !!! Bad end
   !!! ------------------------------------------------------------------------
210 S%ITE=SCARC_NIT
    S%CAPPA=1.0_EB

   !!! ------------------------------------------------------------------------
   !!! Print convergence information
   !!! ------------------------------------------------------------------------
250 CONTINUE
!250 CALL SCARC_INFO_CAPPA(S%NLEVEL)

IF (SCARC_DEBUG>=2) WRITE(9,*) 'SCARC-Multigrid: convergence rate=',S%CAPPA
!IF (NM==1) write(*,*) 'SCARC_MG2D      : Convergence rate=',S%CAPPA, ', Number of iterations ',S%ITE
  
ENDIF
 
DO K = 1, SLHI%KBAR
   DO I = 1, SLHI%IBAR
      M%PRHS (I, 1, K) = SLHI%X(I, 1, K)
   ENDDO
ENDDO
 
TUSED_SCARC(7,NM)=TUSED_SCARC(7,NM)+SECOND()-TNOW_MG2D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_MG2D

END SUBROUTINE SCARC_MG2D




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Perform smoothing in SCARC_MG2D
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!SUBROUTINE SCARC_MASTER_INFO(ITE,RES)
!
!REAL (EB):: RES
!INTEGER  :: ITE
!
!! if master process: receive from slaves
!IF (MYID==0) THEN
!
!    WRITE(9,*) 'Hallo, ich bin der Master'
!! if slave process: send to master
!ELSE
!    WRITE(9,*) 'Hallo, ich bin ein Slave'
!ENDIF
!
!END SUBROUTINE SCARC_MASTER_INFO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Perform smoothing in SCARC_MG2D
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SMOOTHER2D(AG, X, F, D, BMATVEC, ILEVEL, NM)

REAL (EB), POINTER, DIMENSION (:, :)    :: AG
REAL (EB), POINTER, DIMENSION (:, :, :) :: X, F, D
REAL (EB):: RES, OMEGA_SMOOTH, RES_SMOOTH, RESIN_SMOOTH, CAPPA_SMOOTH, EPS_SMOOTH, EPS_REL
INTEGER :: NM, NIT_SMOOTH, ITE_SMOOTH
INTEGER :: ILEVEL, ITYPE
LOGICAL :: BMATVEC, BL2NORM=.TRUE., BRELERR=.TRUE.
REAL(EB):: TNOW_SMOOTHER2D


TNOW_SMOOTHER2D = SECOND()

SL => S%SLEVEL(ILEVEL)

CALL SCARC_COMPARE_SINGLE_LEV (X, 'SARC', 1, 'SMOX  ',0,NM,ILEVEL)
CALL SCARC_COMPARE_SINGLE_LEV (F, 'SARC', 1, 'SMOF  ',0,NM,ILEVEL)

ITYPE=NCOM_TYPE_GLOBAL  ! ??????????????
IF (BMATVEC) THEN
   CALL SCARC_MATVEC2D (AG, X, D, 1.0_EB, 0.0_EB, NM, ITYPE, ILEVEL, 3)   ! global MATVEC
   CALL SCARC_VECADD2D (F, D, 1.0_EB, -1.0_EB, ILEVEL, NM)
   CALL SCARC_COMPARE_SINGLE_LEV (D, 'SARC', 1, 'SMOD  ',0,NM,ILEVEL)
   IF (BL2NORM) CALL SCARC_L2NORM2D (D, RES, ILEVEL, NM, ITYPE)
ENDIF

! still experimental - must be possible to read in those parameters ...
ITE_SMOOTH=0
NIT_SMOOTH=5
OMEGA_SMOOTH=0.8E+0_EB
EPS_SMOOTH=1.0E-2_EB
EPS_REL=1.0E+0_EB
RESIN_SMOOTH=RES
!WRITE(*,1000) ITE_SMOOTH, RESIN_SMOOTH

IF (RESIN_SMOOTH>1.0E+4_EB) GOTO 210 

SMOOTH_LOOP2D: DO ITE_SMOOTH=1,NIT_SMOOTH
 
   !Y=0.0_EB    ! nötig ??
   !WRITE(9,*) 'SMOOTHER, ITERATION ',ITE_SMOOTH
   !WRITE(*,*) 'HALLO: MATRIX-SHIFT ODER NICHT ??*'

   !CALL SCARC_MATRIX_SHIFT(ILEVEL)

   ! hier muss noch ein Auswahlmechanismus hin ...
   !CALL SCARC_PRECON_JACOBI2D(AG, D, ILEVEL, NM)
   !CALL SCARC_PRECON_SSOR2D(AG, D, ILEVEL, NM)
   CALL SCARC_PRECON_GSTRIX2D (SL%DD, SL%DL, SL%DU, SL%LD, SL%DAUX, D, ILEVEL, NM)

   !CALL SCARC_UPDATE_BDRY

   CALL SCARC_VECADD2D(D, X, OMEGA_SMOOTH, 1.0_EB, ILEVEL, NM)

   CALL SCARC_MATVEC2D (AG, X, D, 1.0_EB, 0.0_EB, NM, ITYPE, ILEVEL, 3)   ! global MATVEC
   CALL SCARC_VECADD2D (F, D, 1.0_EB, -1.0_EB, ILEVEL, NM)
   IF (BL2NORM) CALL SCARC_L2NORM2D (D, RES_SMOOTH, ILEVEL, NM, ITYPE)

   !WRITE(*,1000) ITE_SMOOTH, RES_SMOOTH
   !WRITE(9,1000) ITE_SMOOTH, RES_SMOOTH

   IF (BL2NORM) THEN
      IF (BRELERR) THEN
         IF (RES_SMOOTH<=RESIN_SMOOTH*EPS_SMOOTH) GOTO 200
      ELSE
         IF (RES_SMOOTH<=EPS_SMOOTH.AND.RES_SMOOTH<=RESIN_SMOOTH*EPS_REL) GOTO 200
      ENDIF

      IF (RES_SMOOTH>1.0E+6_EB) goto 210
   ENDIF
ENDDO SMOOTH_LOOP2D

200 IF (RES_SMOOTH>1.0E-70_EB) THEN
       CAPPA_SMOOTH=RES_SMOOTH/RESIN_SMOOTH
       CAPPA_SMOOTH=CAPPA_SMOOTH**(1.0_EB/REAL(ITE_SMOOTH,EB))
    ELSE
       CAPPA_SMOOTH=0.0_EB
    ENDIF

    GOTO 250

210 ITE_SMOOTH=NIT_SMOOTH
    CAPPA_SMOOTH=1.0_EB

    GOTO 250

250 CONTINUE
    !WRITE(*,*) 'SCARC_SMOOTHER2D: Convergence rate ', CAPPA_SMOOTH, ' Iterations ',ITE_SMOOTH

!WRITE(9,*) 'TRARA3, LEAVING SMOOTHER2D, CAPPA_SMOOTH=',CAPPA_SMOOTH

TUSED_SCARC(8,NM)=TUSED_SCARC(8,NM)+SECOND()-TNOW_SMOOTHER2D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_SMOOTHER2D

END SUBROUTINE SCARC_SMOOTHER2D


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Solve global coarse grid problem by cg-method
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_COARSE2D (ILEVEL, NM)
 
INTEGER :: NM, ILEVEL
INTEGER :: ITE, IREL, ICONV, ITYPE, IERR, IPRECON
REAL (EB) :: SIGMA0, SIGMA1, ALPHA, GAMMA
REAL (EB) :: TNOW_COARSE2D
 
 
TNOW_COARSE2D = SECOND()

SL => S%SLEVEL(ILEVEL)
 
IERR = 0
 
! initialize working vectors
SL%G = 0.0_EB
SL%Y = 0.0_EB
SL%R = 0.0_EB
SL%D = 0.0_EB
SL%X = 0.0_EB
 
!!!
!!! initialize some method parameters
!!!
ITYPE = NCOM_TYPE_GLOBAL
IREL = 0
ICONV = 0
IPRECON = 2
 
!CALL SCARC_SETBDRY2D (ILEVEL, NM)
 
 
!!!
!!! calculate initial defect Y = B - A*X and get l2-norm of it
!!!
CALL SCARC_MATVEC2D (SL%AG, SL%X, SL%R, 1.0_EB, 0.0_EB, NM, ITYPE, ILEVEL, 3)
CALL SCARC_VECADD2D (SL%F, SL%R,-1.0_EB, 1.0_EB, ILEVEL, NM)
CALL SCARC_L2NORM2D (SL%R, SL%RES, ILEVEL, NM, ITYPE)
SL%RESIN = SL%RES
 
IF (S%BPRECON) THEN
   IF (SCARC_PRECON == 1) THEN
      WRITE (*,*) 'IPRECON=1 not yet implemented!'
   ELSE IF (SCARC_PRECON == 2) THEN
      CALL SCARC_VECCOPY2D (SL%R, SL%G, ILEVEL, NM)
      CALL SCARC_PRECON_JACOBI2D (SL%AG, SL%G, ILEVEL, NM)
   ELSE IF (SCARC_PRECON == 3) THEN
      CALL SCARC_VECCOPY2D (SL%R, SL%G, ILEVEL, NM)
      CALL SCARC_PRECON_GS2D (SL%AG, SL%G, ILEVEL, NM)
   ELSE IF (SCARC_PRECON == 4) THEN
      CALL SCARC_VECCOPY2D (SL%R, SL%G, ILEVEL, NM)
      CALL SCARC_PRECON_SSOR2D (SL%AG, SL%G, ILEVEL, NM)
   ENDIF
   CALL SCARC_SCALPROD2D (SL%R, SL%G, SIGMA0, ILEVEL, NM)
ELSE
   SIGMA0 = SL%RESIN * SL%RESIN
ENDIF
CALL SCARC_VECADD2D (SL%G, SL%D,-1.0_EB, 0.0_EB, ILEVEL, NM)
 
!!!
!!! start defect correction loop
!!!
CG_LOOP2D: DO ITE = 1, SCARC_NIT
 
!!!  calculate new defect and get L2-norm of it
   CALL SCARC_MATVEC2D (SL%AG, SL%D, SL%Y, 1.0_EB, 0.0_EB, NM, ITYPE, ILEVEL, 5)
   CALL SCARC_SCALPROD2D (SL%D, SL%Y, ALPHA, ILEVEL, NM)

   ALPHA = SIGMA0 / ALPHA
 
   CALL SCARC_VECADD2D (SL%D, SL%X, ALPHA, 1.0_EB, ILEVEL, NM)
   CALL SCARC_VECADD2D (SL%Y, SL%R, ALPHA, 1.0_EB, ILEVEL, NM)
   CALL SCARC_L2NORM2D (SL%R, SL%RES, ILEVEL, NM, ITYPE)

   IF (SCARC_DEBUG >= 1) write (9, '(a,i4,a,f25.16,a,f25.16)') &
                        '### SCARC_COARSE2D, ite=', ITE, ' RES=', SL%RES, '  ####'
 
!!!  stop in case of convergence or divergence
   IF (IREL == 1) THEN
      IF (SL%RES <= SL%RESIN*SCARC_EPS) ICONV = 1
   ELSE
      IF (SL%RES <= SCARC_EPS .AND. SL%RES <= SL%RESIN*SCARC_EPS_REL) ICONV = 1
   ENDIF
   IF (SL%RES > 10000.0_EB) ICONV = - 1
   IF (Abs(ICONV) == 1) EXIT CG_LOOP2D
 
!!! preconditioning
   IF (S%BPRECON) THEN
      IF (SCARC_PRECON == 1) THEN
         WRITE (*,*) 'IPRECON=1 not yet implemented!'
      ELSE IF (SCARC_PRECON == 2) THEN
         CALL SCARC_VECCOPY2D (SL%R, SL%G, ILEVEL, NM)
         CALL SCARC_PRECON_JACOBI2D (SL%AG, SL%G, ILEVEL, NM)
      ELSE IF (SCARC_PRECON == 3) THEN
         CALL SCARC_VECCOPY2D (SL%R, SL%G, ILEVEL, NM)
         CALL SCARC_PRECON_GS2D (SL%AG, SL%G, ILEVEL, NM)
      ELSE IF (SCARC_PRECON == 4) THEN
         CALL SCARC_VECCOPY2D (SL%R, SL%G, ILEVEL, NM)
         CALL SCARC_PRECON_SSOR2D (SL%AG, SL%G, ILEVEL, NM)
      ENDIF
      CALL SCARC_SCALPROD2D (SL%R, SL%G, SIGMA1, ILEVEL, NM)
   ELSE
      SIGMA1 = SL%RES * SL%RES
   ENDIF
 
   GAMMA = SIGMA1 / SIGMA0
   SIGMA0 = SIGMA1
   CALL SCARC_VECADD2D (SL%G, SL%D,-1.0_EB, GAMMA, ILEVEL, NM)
 
ENDDO CG_LOOP2D
 
!!! 'bad end'
IF (ICONV ==-1) THEN
   ITE = - 1
   SL%CAPPA = 1.0_EB
 
!!! 'good end'
ELSE
   IF (SL%RESIN >= 1.0E-70_EB) THEN
      SL%RESIN = SL%RES / SL%RESIN
      SL%CAPPA = SL%RESIN ** (1.0_EB/ITE)
   ELSE IF (ITE == 0) THEN
      SL%CAPPA = 0.0_EB
   ENDIF
ENDIF
 
IF (SCARC_DEBUG >= 1) THEN
   WRITE (9, '(a,i4,a,f12.6,a,f12.5)') 'SCARC_COARSE: Ite=', ITE, ': Res=', SL%RES,'  Kappa=', SL%CAPPA
ENDIF
 
TUSED_SCARC(9,NM)=TUSED_SCARC(9,NM)+SECOND()-TNOW_COARSE2D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_COARSE2D

END SUBROUTINE SCARC_COARSE2D


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Perform smoothing in SCARC_MG2D
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!SUBROUTINE SCARC_DEFCOR(AG, X, F, ILEVEL, NM)

!REAL (EB), POINTER, DIMENSION (:, :)    :: AG
!REAL (EB), POINTER, DIMENSION (:, :, :) :: X, F
!INTEGER :: ILEVEL, NM


!WRITE(*,*) 'MUSS NOCH PROGRAMMIERT WERDEN'
 
!END SUBROUTINE SCARC_DEFCOR

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Perform global cg-method based on global possion-matrix in 3D
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_CG3D (NM)
 
INTEGER :: NM, ILEVEL
INTEGER :: ITE, IREL, ICONV, ITYPE, I, J, K, IERR, IPRECON
REAL (EB):: SIGMA0, SIGMA1, ALPHA, GAMMA
REAL (EB) :: TNOW_CG3D
 
TYPE (MESH_TYPE), POINTER :: M
 
INCLUDE 'mpif.h'
 
TNOW_CG3D = SECOND()

S => SCARC (NM)
M => MESHES (NM)
 
NREQ_FACE = 0
NREQ_EDGE = 0
NREQ_DIAG = 0

IERR = 0
 
ILEVEL=S%NLEVEL_MAX
 
! initialize working vectors
SL%G = 0.0_EB
SL%Y = 0.0_EB
SL%R = 0.0_EB
SL%D = 0.0_EB
 
 
DO K = 0, KBAR + 1
   DO J = 0, JBAR + 1
      DO I = 0, IBAR + 1
         SL%X (I, J, K) = M%H(I, J, K)
      ENDDO
   ENDDO
ENDDO
 
SL%F = 0.0_EB
DO K = 1, KBAR + 1
   DO J = 1, JBAR + 1
      DO I = 1, IBAR + 1
         SL%F (I, J, K) = M%PRHS(I, J, K)
      ENDDO
   ENDDO
ENDDO
 
 
IF (SCARC_DEBUG >= 2) THEN
   WRITE (9,*) 'Starting SCARC_CG3D'
   WRITE (9,*) 'IBAR=:', IBAR
   WRITE (9,*) 'JBAR=:', JBAR
   WRITE (9,*) 'KBAR=:', KBAR
   CALL SCARC_COMPARE_SINGLE (M%H, 'SARC', 1, 'H1    ',0)
   CALL SCARC_COMPARE_SINGLE (SL%X, 'SARC', 1, 'X1    ',0)
   CALL SCARC_COMPARE_SINGLE (SL%F, 'SARC', 1, 'F1    ',0)
   CALL flush (9)
ENDIF
 
!!!
!!! initialize some method parameters
!!!
ITYPE = NCOM_TYPE_GLOBAL
IREL = 0
ICONV = 0
IPRECON = 2
 
!SL%NIT    =5000
!SL%EPS    =1.0E-16_EB
!SL%EPS_REL=1.0E-2_EB
 
 
!!!
!!! initialize boundary values
!!!
CALL SCARC_SETBDRY3D (ILEVEL, NM)
 
 
!!!
!!! calculate initial defect Y = B - A*X and get l2-norm of it
!!!
CALL SCARC_MATVEC3D (SL%AG, SL%X, SL%R, 1.0_EB, 0.0_EB, NM, ITYPE, ILEVEL, 3)
CALL SCARC_VECADD3D (SL%F, SL%R,-1.0_EB, 1.0_EB, ILEVEL, NM)
CALL SCARC_L2NORM3D (SL%R, SL%RES, ILEVEL, NM, ITYPE)
SL%RESIN = SL%RES
 
IF (SCARC_DEBUG >= 2) THEN
   WRITE (9,*) 'Initial DEFECT '
   WRITE (9,*) 'RES=', SL%RES
   WRITE (9,*) 'SL%X:'
   CALL SCARC_COMPARE_SINGLE (SL%X, 'SARC', 1, 'X3    ',0)
   CALL SCARC_COMPARE_SINGLE (SL%F, 'SARC', 1, 'F3    ',0)
   CALL SCARC_COMPARE_SINGLE (SL%R, 'SARC', 1, 'R3    ',0)
   CALL flush (9)
ENDIF
 
 
IF (S%BPRECON) THEN
 
   IF (SCARC_PRECON == 1) THEN
      WRITE (*,*) 'IPRECON=1 not yet implemented!'
   ELSE IF (SCARC_PRECON == 2) THEN
      CALL SCARC_VECCOPY3D (SL%R, SL%G, ILEVEL, NM)
      CALL SCARC_PRECON_JACOBI3D (SL%AG, SL%G, ILEVEL, NM)
   ELSE IF (SCARC_PRECON == 3) THEN
      CALL SCARC_VECCOPY3D (SL%R, SL%G, ILEVEL, NM)
      CALL SCARC_PRECON_GS3D (SL%AG, SL%G, ILEVEL, NM)
   ELSE IF (SCARC_PRECON == 4) THEN
      CALL SCARC_VECCOPY3D (SL%R, SL%G, ILEVEL, NM)
      CALL SCARC_PRECON_SSOR3D (SL%AG, SL%G, ILEVEL, NM)
   ENDIF
 
   CALL SCARC_SCALPROD3D (SL%R, SL%G, SIGMA0, ILEVEL, NM)
 
ELSE
   SIGMA0 = SL%RESIN * SL%RESIN
ENDIF
 
CALL SCARC_VECADD3D (SL%G, SL%D,-1.0_EB, 0.0_EB, ILEVEL, NM)
 
IF (SCARC_DEBUG >= 2) THEN
   WRITE (9,*) 'After initial Precon '
   WRITE (9,*) 'Sigma0=', SIGMA0
   WRITE (9,*) 'SL%X:'
   CALL SCARC_COMPARE_SINGLE (SL%R, 'SARC', 1, 'R4    ',0)
   CALL SCARC_COMPARE_SINGLE (SL%G, 'SARC', 1, 'G4    ',0)
   CALL SCARC_COMPARE_SINGLE (SL%D, 'SARC', 1, 'D4    ',0)
   CALL flush (9)
ENDIF
 
!!!
!!! start defect correction loop
!!!
CG_LOOP3D: DO ITE = 1, SCARC_NIT
 
   IF (SCARC_DEBUG >= 2) write (9,*) '=========================== STARTING CG-&
  &LOOP ', ITE
 
!!!  calculate new defect and get L2-norm of it
   CALL SCARC_MATVEC3D (SL%AG, SL%D, SL%Y, 1.0_EB, 0.0_EB, NM, ITYPE, ILEVEL, 5)
 
   IF (SCARC_DEBUG >= 2) CALL SCARC_COMPARE_SINGLE (SL%Y, 'SARC', 1, 'Y4    ',0)
 
   CALL SCARC_SCALPROD3D (SL%D, SL%Y, ALPHA, NLEVEL1, NM)
   ALPHA = SIGMA0 / ALPHA
 
   IF (SCARC_DEBUG >= 2) CALL SCARC_COMPARE_SINGLE (SL%Y, 'SARC', 1, 'Y5    ',0)
 
   CALL SCARC_VECADD3D (SL%D, SL%X, ALPHA, 1.0_EB, ILEVEL, NM)
   CALL SCARC_VECADD3D (SL%Y, SL%R, ALPHA, 1.0_EB, ILEVEL, NM)
 
   CALL SCARC_L2NORM3D (SL%R, SL%RES, ILEVEL, NM, ITYPE)
 
   IF (SCARC_DEBUG >= 2) THEN
      WRITE (9,*) 'New DEFECT in iteration ', ITE, ': ALPHA=', ALPHA
      WRITE (9,*) 'New Iterate in iteration ', ITE
      WRITE (9,*) 'RES=', SL%RES
      CALL flush (9)
   ENDIF
   CALL SCARC_COMPARE_SINGLE (SL%X, 'SARC', 1, 'X5    ',0)
   CALL SCARC_COMPARE_SINGLE (SL%R, 'SARC', 1, 'R5    ',0)
 
   IF (NM == 1 .AND. SCARC_DEBUG >= 2) write (9, '(a,i4,a,f25.16,a,f25.16)') '&
  &### SCARC_CG, ite=', ITE, ' RES=', SL%RES, '  ####'
   IF (Mod(ITE, 200) == 0 .AND. SCARC_DEBUG == 1) write (9, '(a,i4,a,f25.16,a,f2&
  &5.16)') '### SCARC_CG, ite=', ITE, ' RES=', SL%RES, '  ####'
 
!!!  stop in case of convergence or divergence
   IF (IREL == 1) THEN
      IF (SL%RES <= SL%RESIN*SCARC_EPS) ICONV = 1
   ELSE
      IF (SL%RES <= SCARC_EPS .AND. SL%RES <= SL%RESIN*SCARC_EPS_REL) ICONV = 1
   ENDIF
   IF (SL%RES > 100000.0_EB) ICONV = - 1
   IF (Abs(ICONV) == 1) EXIT CG_LOOP3D
 
!!! preconditioning
   IF (S%BPRECON) THEN
 
      IF (SCARC_PRECON == 1) THEN
         WRITE(*,*) 'HALLO, HIER ALLE KOMMUNIKATIONSTYPEN EINZELN AUFLISTEN !!!'
         CALL SCARC_RECEIVE (NCOM_FULL, ILEVEL)  !Aufruf 2,0
         CALL SCARC_EXCHANGE (NCOM_FULL, ILEVEL, NTYPE_NONE)
         WRITE (*,*) 'IPRECON=1 not yet implemented!'
!     CALL SCARC_CG0(SL%R,NM)
      ELSE IF (SCARC_PRECON == 2) THEN
         WRITE(*,*) 'HALLO, HIER ALLE KOMMUNIKATIONSTYPEN EINZELN AUFLISTEN !!!'
         CALL SCARC_VECCOPY3D (SL%R, SL%G, ILEVEL, NM)
         CALL SCARC_PRECON_JACOBI3D (SL%AG, SL%G, ILEVEL, NM)
      ELSE IF (SCARC_PRECON == 3) THEN
         WRITE(*,*) 'HALLO, HIER ALLE KOMMUNIKATIONSTYPEN EINZELN AUFLISTEN !!!'
         CALL SCARC_VECCOPY3D (SL%R, SL%G, ILEVEL, NM)
         CALL SCARC_PRECON_GS3D (SL%AG, SL%G, ILEVEL, NM)
      ELSE IF (SCARC_PRECON == 4) THEN
         WRITE(*,*) 'HALLO, HIER ALLE KOMMUNIKATIONSTYPEN EINZELN AUFLISTEN !!!'
         CALL SCARC_VECCOPY3D (SL%R, SL%G, ILEVEL, NM)
         CALL SCARC_PRECON_SSOR3D (SL%AG, SL%G, ILEVEL, NM)
      ENDIF
 
      CALL SCARC_SCALPROD3D (SL%R, SL%G, SIGMA1, NLEVEL1, NM)
 
   ELSE
      SIGMA1 = SL%RES * SL%RES
   ENDIF
 
   IF (SCARC_DEBUG >= 3) CALL SCARC_COMPARE_SINGLE (SL%G, 'SARC', 1, 'G6    ',0)
 
   GAMMA = SIGMA1 / SIGMA0
   SIGMA0 = SIGMA1
 
   CALL SCARC_VECADD3D (SL%G, SL%D,-1.0_EB, GAMMA, ILEVEL, NM)
 
   IF (SCARC_DEBUG >= 2) THEN
      WRITE (9,*) 'After Precon in iteration ', ITE
      WRITE (9,*) 'Sigma0=', SIGMA0
      WRITE (9,*) 'Gamma =', GAMMA
      CALL SCARC_COMPARE_SINGLE (SL%R, 'SARC', 1, 'R7    ',0)
      CALL SCARC_COMPARE_SINGLE (SL%D, 'SARC', 1, 'D7    ',0)
      CALL flush (9)
   ENDIF
 
ENDDO CG_LOOP3D
 
IF (SCARC_DEBUG >= 2) write (9,*) 'AFTER CG_LOOP3D'
 
 
!!! 'bad end'
IF (ICONV ==-1) THEN
   ITE = - 1
   SL%CAPPA = 1.0_EB
 
!!! 'good end'
ELSE
   IF (SL%RESIN >= 1.0E-70_EB) THEN
      SL%RESIN = SL%RES / SL%RESIN
      SL%CAPPA = SL%RESIN ** (1.0_EB/ITE)
   ELSE IF (ITE == 0) THEN
      SL%CAPPA = 0.0_EB
   ENDIF
ENDIF
 
 
IF (NM == 1) write (*, '(a20,i4,a,e12.6,a,e12.6)') 'SCARC_CG: Ite=', ITE, ': Res&
  &=', SL%RES, ': Kappa=', SL%CAPPA
IF (SCARC_DEBUG >= 2) THEN
   WRITE (9, '(a,i4,a,f12.6,a,f12.5)') 'SCARC_CG: Ite=', ITE, ': Res=', SL%RES, '&
  &: Kappa=', SL%CAPPA
   CALL flush (9)
ENDIF
CALL flush (6)
 
 
DO K = 0, KBAR + 1
   DO J = 0, JBAR + 1
      DO I = 0, IBAR + 1
         M%H (I, J, K) = SL%X(I, J, K)
      ENDDO
   ENDDO
ENDDO
 
DO K = 1, KBAR + 1
   DO J = 1, JBAR + 1
      DO I = 1, IBAR + 1
         M%PRHS (I, J, K) = SL%X(I, J, K)
      ENDDO
   ENDDO
ENDDO
 
 
IF (SCARC_DEBUG >= 2) THEN
   WRITE (9,*) '-------SOL: SL%X:'
   CALL SCARC_COMPARE_SINGLE (SL%X, 'SARC', 1, 'XFINAL',0)
   CALL flush (9)
ENDIF

TUSED_SCARC(6,NM)=TUSED_SCARC(6,NM)+SECOND()-TNOW_CG3D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_CG3D
 
END SUBROUTINE SCARC_CG3D
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Perform matrix-vector multiplication with full matrix
!!! corresponding to band-wise storage technique
!!!
!!!  Y = A1 * S%AGLOB * X + A2 * Y
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_MATVEC2D (A, X, Y, A1, A2, NM, ITYPE, ILEVEL, IMV)
 
REAL (EB) :: A1, A2
REAL (EB), POINTER, DIMENSION (:, :) :: A
REAL (EB), POINTER, DIMENSION (:, :, :) :: X, Y
REAL (EB):: TNOW_MATVEC2D
INTEGER :: NM, ITYPE
INTEGER :: I, K, IC, IMV, ILEVEL
 
TNOW_MATVEC2D = SECOND()
 
SL => S%SLEVEL(ILEVEL)
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! print debug messages !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (SCARC_DEBUG > 2) THEN
   WRITE (9,*) 'BEFORE MATVEC: '
   WRITE (9,*) 'TWO_D=', TWO_D
   WRITE (9,*) 'KBAR=', KBAR
   WRITE (9,*) 'JBAR=', JBAR
   WRITE (9,*) 'IBAR=', IBAR
   WRITE (9,*) 'A1=', A1
   WRITE (9,*) 'A2=', A2
   CALL SCARC_COMPARE_SINGLE_LEV (X, 'SARC', 1, 'matx12',0, NM, ILEVEL)
   CALL SCARC_COMPARE_SINGLE_LEV (Y, 'SARC', 1, 'maty12',0, NM, ILEVEL)
ENDIF
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform 2D-matrix-vector-multiplication
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DO K = 1, SL%KBAR
IF (SCARC_DEBUG>=6) WRITE(9,*) '================================== K=',K
   DO I = 1, SL%IBAR
      IC = (K-1) * SL%IBAR + I
      Y (I, 1, K) =    A1 * (  A(IC, 1)*X(I  , 1, K  ) &
                             + A(IC, 2)*X(I  , 1, K-1) &
                             + A(IC, 3)*X(I-1, 1, K  ) &
                             + A(IC, 4)*X(I+1, 1, K  ) &
                             + A(IC, 5)*X(I  , 1, K+1))&
                     + A2 * Y (I, 1, K)
   IF (SCARC_DEBUG>=6) WRITE(9,'(a,5f9.2,a,5f9.2,a,f9.2)') 'A: ',A(IC,1),A(IC,2),A(IC,3),A(IC,4),A(IC,5),&
                       '   X: ',X(I,1,K),X(I,1,K-1),X(I-1,1,K),X(I+1,1,K),X(I,1,K+1),'   ===> Y: ',Y(I,1,K)
   ENDDO
ENDDO
IF (SCARC_DEBUG>=2) WRITE(9,*) 'UFF8'
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! print debug messages !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (SCARC_DEBUG > 2) THEN
   CALL SCARC_COMPARE_SINGLE_LEV (Y, 'SARC', 1, 'maty13',0,NM, ILEVEL)
   CALL flush (9)
ENDIF
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform data exchange at interior boundaries
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF ((ITYPE == NCOM_TYPE_GLOBAL) .AND. (NMESHES > 1)) THEN
   NREQ_FACE = 0
   NREQ_EDGE = 0
   NREQ_DIAG = 0
   CALL SCARC_RECEIVE (NCOM_MATV,  ILEVEL)   
   CALL SCARC_EXCHANGE (NCOM_MATV,  ILEVEL, IMV)
   CALL SCARC_COMPARE_SINGLE0 (Y, 'SARC', 1, 'maty14',0)
ENDIF
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! print debug messages !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (SCARC_DEBUG > 2) THEN
   WRITE (9,*) '------->> HALLO AFTER SCARC_EXCHANGE: Y:'
   CALL flush (9)
   !CALL SCARC_COMPARE_SINGLE0 (Y, 'SARC', 1, 'maty1 ',0)
ENDIF
 
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
SUBROUTINE SCARC_MATVEC3D (A, X, Y, A1, A2, NM, ITYPE, ILEVEL, IMV)
 
REAL (EB) :: A1, A2
REAL (EB), POINTER, DIMENSION (:, :) :: A
REAL (EB), POINTER, DIMENSION (:, :, :) :: X, Y
REAL (EB):: TNOW_MATVEC3D
INTEGER :: NM, ITYPE
INTEGER :: I, J, K, LL, IC, IMV, ILEVEL
 
 
TNOW_MATVEC3D = SECOND()
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! print debug messages !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (SCARC_DEBUG > 2) THEN
   WRITE (9,*) 'BEFORE MATVEC: '
!  CALL SCARC_COMPARE_SINGLE(X,'SARC',1,'matX  ',0)
!  CALL SCARC_COMPARE_SINGLE(Y,'SARC',1,'maty  ',0)
   WRITE (9,*) 'TWO_D=', TWO_D
   WRITE (9,*) 'KBAR=', KBAR
   WRITE (9,*) 'JBAR=', JBAR
   WRITE (9,*) 'IBAR=', IBAR
ENDIF
 
!X=REAL(NM,EB)
!Y=0.0_EB
!X=1.0_EB
!Y=0.0_EB
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform 3D-matrix-vector-multiplication
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DO K = 1, KBAR
   IF (SCARC_DEBUG > 2) write (9,*) 'K=', K, ' ==============================&
  &======='
   DO J = 1, JBAR
      DO I = 1, IBAR
         IC = (K-1) * IBAR * JBAR + (J-1) * IBAR + I
!WRITE(9,'(8f22.16)') y(i,j,k),(A(IC,kk),kk=1,7)
         Y (I, J, K) = A1 * (A(IC, 1)*X(I, J, K-1)+A(IC, 2)*X(I, J-1, K)+A(IC, &
          3)*X(I-1, J, K)+A(IC, 4)*X(I, J, K)+A(IC, 5)*X(I+1, J, K)+A(IC, &
          6)*X(I, J+1, K)+A(IC, 7)*X(I, J, K+1)) + A2 * Y (I, J, K)
         IF (SCARC_DEBUG > 2) write (9, '(a,i2,a,i2,a,i2,a,i3,a,7f6.0,a,f22.1&
        &6)') '(', I, ',', J, ',', K, '):A(', IC, '):', (A(IC, LL), LL=1, 7), '-&
        &-> Y=', Y (I, J, K)
      ENDDO
      IF (SCARC_DEBUG > 2) write (9,*) '-------------------------------------&
     &'
   ENDDO
ENDDO
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! print debug messages !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (SCARC_DEBUG > 2) THEN
   WRITE (9,*) '------->> HALLO BEFORE SCARC_EXCHANGE: Y:'
!  CALL SCARC_COMPARE_SINGLE(Y,'SARC',1,'maty2 ',0)
   CALL flush (9)
ENDIF
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform data exchange at interior boundaries
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF ((ITYPE == NCOM_TYPE_GLOBAL) .AND. (NMESHES > 1)) THEN
   NREQ_FACE = 0
   NREQ_EDGE = 0
   NREQ_DIAG = 0
         WRITE(*,*) 'HALLO, HIER ALLE KOMMUNIKATIONSTYPEN EINZELN AUFLISTEN !!!'
   CALL SCARC_RECEIVE (NCOM_MATV,  ILEVEL)  !Aufruf 1, IMV
   CALL SCARC_EXCHANGE (NCOM_MATV,  ILEVEL, IMV)
ENDIF
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! print debug messages !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (SCARC_DEBUG > 2) THEN
   WRITE (9,*) '------->> HALLO AFTER SCARC_EXCHANGE: Y:'
!  CALL SCARC_COMPARE_SINGLE(Y,'SARC',1,'maty3 ',0)
   CALL flush (9)
ENDIF
 
TUSED_SCARC(10,NM)=TUSED_SCARC(10,NM)+SECOND()-TNOW_MATVEC3D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_MATVEC3D

END SUBROUTINE SCARC_MATVEC3D
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Restrict vector X from level ILEVEL to vector Y on level ILEVEL-1
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_RESTRICTION2D (NM, X_LO, X_HI, IBAR_LO, KBAR_LO)
 
REAL (EB), POINTER, DIMENSION (:, :, :) :: X_LO, X_HI
INTEGER :: IBAR_LO, KBAR_LO
INTEGER :: I_LO, K_LO, I1_HI, I2_HI, K1_HI, K2_HI,NM
REAL(EB):: TNOW_RESTRICTION2D


TNOW_RESTRICTION2D = SECOND()
DO K_LO = 1, IBAR_LO
   DO I_LO = 1, KBAR_LO
      
      I1_HI=2*I_LO-1
      I2_HI=2*I_LO
      K1_HI=2*K_LO-1
      K2_HI=2*K_LO

      X_LO (I_LO, 1, K_LO) = 0.25_EB * (X_HI(I1_HI,1,K1_HI) + X_HI(I1_HI,1,K2_HI) + &
                                        X_HI(I2_HI,1,K2_HI) + X_HI(I2_HI,1,K2_HI) )

      IF (SCARC_DEBUG>=4) THEN
          WRITE(9,*) 'X_LO(',I_LO,',',1,',',K_LO,')=',X_LO(I_LO,1,K_LO), I1_HI,I2_HI,K1_HI,K2_HI
      ENDIF
   ENDDO
ENDDO

TUSED_SCARC(11,NM)=TUSED_SCARC(11,NM)+SECOND()-TNOW_RESTRICTION2D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_RESTRICTION2D
END SUBROUTINE SCARC_RESTRICTION2D


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Restrict vector X from level ILEVEL to vector Y on level ILEVEL-1
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PROLONGATION2D (NM,X_LO, X_HI, ILEV_HI)
 
REAL (EB), POINTER, DIMENSION (:, :, :) :: X_LO, X_HI
INTEGER :: ILEV_HI, I_LO, K_LO, I1_HI, I2_HI, K1_HI, K2_HI, NM
REAL(EB):: TNOW_PROLONGATION2D

TNOW_PROLONGATION2D = SECOND()

SLLO => S%SLEVEL(ILEV_HI-1)
SLHI => S%SLEVEL(ILEV_HI)

CALL SCARC_COMPARE_SINGLE_LEV (X_LO, 'SARC', 1, 'PX_LO ',0,NM,ILEV_HI-1)
DO K_LO = 1, SLLO%KBAR
   DO I_LO = 1, SLLO%IBAR
      
      I1_HI=2*I_LO-1
      I2_HI=2*I_LO
      K1_HI=2*K_LO-1
      K2_HI=2*K_LO

      X_HI (I1_HI, 1, K1_HI) = X_LO (I_LO, 1, K_LO)
      X_HI (I1_HI, 1, K2_HI) = X_LO (I_LO, 1, K_LO)
      X_HI (I2_HI, 1, K1_HI) = X_LO (I_LO, 1, K_LO)
      X_HI (I2_HI, 1, K2_HI) = X_LO (I_LO, 1, K_LO)

      IF (SCARC_DEBUG>=4) THEN
          WRITE(9,*) 'X(',I1_HI,',',1,',',K1_HI,')=',X_LO(I_LO,1,K_LO)
          WRITE(9,*) 'X(',I1_HI,',',1,',',K2_HI,')=',X_LO(I_LO,1,K_LO)
          WRITE(9,*) 'X(',I2_HI,',',1,',',K1_HI,')=',X_LO(I_LO,1,K_LO)
          WRITE(9,*) 'X(',I2_HI,',',1,',',K2_HI,')=',X_LO(I_LO,1,K_LO)
      ENDIF
   ENDDO
ENDDO
CALL SCARC_COMPARE_SINGLE_LEV (X_HI, 'SARC', 1, 'PX_LO ',0,NM,ILEV_HI)

TUSED_SCARC(12,NM)=TUSED_SCARC(12,NM)+SECOND()-TNOW_PROLONGATION2D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_PROLONGATION2D
END SUBROUTINE SCARC_PROLONGATION2D


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

   CALL SCARC_COMPARE_SINGLE_LEV (F, 'SARC', 1, 'FRES1 ',0,NM,ILEVEL)
DO IW_FACE=1,SL%NW_FACE

   IF (SL%IJKW_FACE(9,IW_FACE)==0) THEN

      I=SL%IJKW_FACE(1,IW_FACE)
      J=SL%IJKW_FACE(2,IW_FACE)
      K=SL%IJKW_FACE(3,IW_FACE)

      F(I,J,K)=0.0_EB

      if (SCARC_DEBUG>=4) WRITE(9,*) 'F(',I,',',J,',',K,')=',F(I,J,K)

   ENDIF

ENDDO
   CALL SCARC_COMPARE_SINGLE_LEV (F, 'SARC', 1, 'FRES2 ',0,NM,ILEVEL)

TUSED_SCARC(13,NM)=TUSED_SCARC(13,NM)+SECOND()-TNOW_BDRY_RESIDUUM
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_BDRY_RESIDUUM
END SUBROUTINE SCARC_BDRY_RESIDUUM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! add two vectors:
!!!
!!! Y = A1 * X + A2 * Y
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_VECADD2D (X, Y, A1, A2, ILEVEL, NM)
 
REAL (EB) :: A1, A2
REAL (EB), POINTER, DIMENSION (:, :, :) :: X, Y
INTEGER :: I, K, ILEVEL, NM
REAL(EB):: TNOW_VECADD2D
 

TNOW_VECADD2D = SECOND()
SL => S%SLEVEL(ILEVEL)
 
IF (SCARC_DEBUG >= 6) THEN
   WRITE (9,*) 'BEFORE VECADD'
   WRITE (9,*) 'A1=', A1
   WRITE (9,*) 'A2=', A2
   WRITE (9,*) ' X:'
   DO K = SL%KBAR + 1, 0, - 1
      WRITE (9, '(5f10.5)') (X(I, 1, K), I=0, SL%IBAR+1)
   ENDDO
   WRITE (9,*) ' Y:'
   DO K = SL%KBAR + 1, 0, - 1
      WRITE (9, '(5f10.5)') (Y(I, 1, K), I=0, SL%IBAR+1)
   ENDDO
ENDIF
 
 
IF (A2 == 0.0_EB) THEN
   IF (A1 /= 1.0_EB) THEN
      DO K = 1, SL%KBAR
         DO I = 1, SL%IBAR
            Y (I, 1, K) = A1 * X (I, 1, K)
 IF (SCARC_DEBUG>=6) WRITE(9,*) 'Y(',I,',1,',K,')=',Y(I,1,K)
         ENDDO
      ENDDO
   ELSE
      DO K = 1, SL%KBAR
         DO I = 1, SL%IBAR
            Y (I, 1, K) = X (I, 1, K)
 IF (SCARC_DEBUG>=6) WRITE(9,*) 'Y(',I,',1,',K,')=',Y(I,1,K)
         ENDDO
      ENDDO
   ENDIF
ELSE
   IF (A1 /= 1.0_EB) THEN
      DO K = 1, SL%KBAR
         DO I = 1, SL%IBAR
            Y (I, 1, K) = A1 * X (I, 1, K) + A2 * Y (I, 1, K)
 IF (SCARC_DEBUG>=6) WRITE(9,*) 'Y(',I,',1,',K,')=',Y(I,1,K)
         ENDDO
      ENDDO
   ELSE
      DO K = 1, SL%KBAR
         DO I = 1, SL%IBAR
            Y (I, 1, K) = X (I, 1, K) + A2 * Y (I, 1, K)
 IF (SCARC_DEBUG>=6) WRITE(9,*) 'Y(',I,',1,',K,')=',Y(I,1,K)
         ENDDO
      ENDDO
   ENDIF
ENDIF
 
 
IF (SCARC_DEBUG >= 6) THEN
   WRITE (9,*) 'AFTER VECADD'
   WRITE (9,*) ' Y:'
   DO K = SL%KBAR + 1, 0, - 1
      WRITE (9, '(5f10.5)') (Y(I, 1, K), I=0, SL%IBAR+1)
   ENDDO
ENDIF
 
TUSED_SCARC(14,NM)=TUSED_SCARC(14,NM)+SECOND()-TNOW_VECADD2D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_VECADD2D
END SUBROUTINE SCARC_VECADD2D
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! add two vectors:
!!!
!!! Y = A1 * X + A2 * Y
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_VECADD3D (X, Y, A1, A2, ILEVEL, NM)
 
REAL (EB) :: A1, A2
REAL (EB), POINTER, DIMENSION (:, :, :) :: X, Y
INTEGER :: I, J, K, ILEVEL, NM
REAL(EB):: TNOW_VECADD3D
 
TNOW_VECADD3D = SECOND()
 
SL => S%SLEVEL(ILEVEL)

IF (SCARC_DEBUG >= 6) THEN
   WRITE (9,*) 'BEFORE VECADD'
   WRITE (9,*) 'A1=', A1
   WRITE (9,*) 'A2=', A2
   WRITE (9,*) ' X:'
   DO K = SL%KBAR + 1, 0, - 1
      WRITE (9, '(5f10.5)') ((X(I, J, K), I=0, SL%IBAR+1), J=SL%JBAR+1, 1,-1)
   ENDDO
   WRITE (9,*) ' Y:'
   DO K = SL%KBAR + 1, 0, - 1
      WRITE (9, '(5f10.5)') ((Y(I, J, K), I=0, SL%IBAR+1), J=SL%JBAR+1, 1,-1)
   ENDDO
ENDIF
 
 
IF (A2 == 0.0_EB) THEN
   IF (A1 /= 1.0_EB) THEN
      DO K = 1, SL%KBAR
         DO J = 1, SL%JBAR
            DO I = 1, SL%IBAR
               Y (I, J, K) = A1 * X (I, J, K)
!WRITE(9,*) 'Y(',I,',J,',K,')=',Y(I,J,K)
            ENDDO
         ENDDO
      ENDDO
   ELSE
      DO K = 1, SL%KBAR
         DO J = 1, SL%JBAR
            DO I = 1, SL%IBAR
               Y (I, J, K) = X (I, J, K)
!WRITE(9,*) 'Y(',I,',J,',K,')=',Y(I,J,K)
            ENDDO
         ENDDO
      ENDDO
   ENDIF
ELSE
   IF (A1 /= 1.0_EB) THEN
      DO K = 1, SL%KBAR
         DO J = 1, SL%JBAR
            DO I = 1, SL%IBAR
               Y (I, J, K) = A1 * X (I, J, K) + A2 * Y (I, J, K)
!WRITE(9,*) 'Y(',I,',J,',K,')=',Y(I,J,K)
            ENDDO
         ENDDO
      ENDDO
   ELSE
      DO K = 1, SL%KBAR
         DO J = 1, SL%KBAR
            DO I = 1, SL%IBAR
               Y (I, J, K) = X (I, J, K) + A2 * Y (I, J, K)
!WRITE(9,*) 'Y(',I,',J,',K,')=',Y(I,J,K)
            ENDDO
         ENDDO
      ENDDO
   ENDIF
ENDIF
 
 
IF (SCARC_DEBUG >= 2) THEN
   WRITE (9,*) 'AFTER VECADD'
   WRITE (9,*) ' Y:'
   DO K = SL%KBAR + 1, 1, - 1
      WRITE (9, '(5f10.5)') ((Y(I, J, K), I=1, SL%IBAR+1), J=SL%JBAR+1, 1,-1)
   ENDDO
ENDIF
 
 
TUSED_SCARC(14,NM)=TUSED_SCARC(14,NM)+SECOND()-TNOW_VECADD3D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_VECADD3D
END SUBROUTINE SCARC_VECADD3D
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! copy vector
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_VECCOPY2D (X, Y, ILEVEL, NM)
 
REAL (EB), POINTER, DIMENSION (:, :, :) :: X, Y
INTEGER :: I, K, ILEVEL, NM
REAL(EB):: TNOW_VECCOPY2D
 
TNOW_VECCOPY2D = SECOND()
 
SL => S%SLEVEL(ILEVEL)
 
IF (SCARC_DEBUG >= 2) THEN
   WRITE (9,*) 'BEFORE VECCOPY2D'
   WRITE (9,*) ' X:'
   DO K = SL%KBAR, 1, - 1
      WRITE (9, '(4f10.4)') (X(I, 1, K), I=1, SL%IBAR)
   ENDDO
   WRITE (9,*) ' Y:'
   DO K = SL%KBAR, 1, - 1
      WRITE (9, '(4f10.4)') (Y(I, 1, K), I=1, SL%IBAR)
   ENDDO
ENDIF
 
 
DO K = 1, SL%KBAR
   DO I = 1, SL%IBAR
      Y (I, 1, K) = X (I, 1, K)
!    WRITE(9,*) 'Y(',I,',1,',K,')=',Y(I,1,K)
   ENDDO
ENDDO
 
 
IF (SCARC_DEBUG >= 2) THEN
   WRITE (9,*) 'AFTER VECCOPY'
   WRITE (9,*) ' Y:'
   DO K = SL%KBAR, 1, - 1
      WRITE (9, '(4f10.4)') (Y(I, 1, K), I=1, SL%IBAR)
   ENDDO
ENDIF
 
TUSED_SCARC(15,NM)=TUSED_SCARC(15,NM)+SECOND()-TNOW_VECCOPY2D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_VECCOPY2D
END SUBROUTINE SCARC_VECCOPY2D
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! copy vector
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_VECCOPY3D (X, Y, ILEVEL, NM)
 
REAL (EB), POINTER, DIMENSION (:, :, :) :: X, Y
INTEGER :: I, J, K, ILEVEL, NM
REAL(EB):: TNOW_VECCOPY3D
 
TNOW_VECCOPY3D = SECOND()
 
SL => S%SLEVEL(ILEVEL)
 
IF (SCARC_DEBUG >= 2) THEN
   WRITE (9,*) 'BEFORE VECCOPY'
   WRITE (9,*) ' X:'
   DO K = SL%KBAR, 1, - 1
      WRITE (9, '(4f10.4)') ((X(I, J, K), I=1, SL%IBAR), J=SL%JBAR, 1,-1)
   ENDDO
   WRITE (9,*) ' Y:'
   DO K = SL%KBAR, 1, - 1
      WRITE (9, '(4f10.4)') ((Y(I, J, K), I=1, SL%IBAR), J=SL%JBAR, 1,-1)
   ENDDO
ENDIF
 
 
DO K = 1, SL%KBAR
   DO J = 1, SL%JBAR
      DO I = 1, SL%IBAR
         Y (I, J, K) = X (I, J, K)
!      WRITE(9,*) 'Y(',I,',J,',K,')=',Y(I,J,K)
      ENDDO
   ENDDO
ENDDO
 
 
IF (SCARC_DEBUG >= 2) THEN
   WRITE (9,*) 'AFTER VECCOPY3D'
   WRITE (9,*) ' Y:'
   DO K = SL%KBAR, 1, - 1
      WRITE (9, '(4f10.4)') ((Y(I, J, K), I=1, SL%IBAR), J=SL%JBAR, 1,-1)
   ENDDO
ENDIF
 
TUSED_SCARC(15,NM)=TUSED_SCARC(15,NM)+SECOND()-TNOW_VECCOPY3D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_VECCOPY3D
END SUBROUTINE SCARC_VECCOPY3D
 
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Compute scalarproduct of vector X and vector Y
!!!
!!! SP = (X,Y)
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SCALPROD2D (X, Y, SP, ILEVEL, NM)
 
REAL (EB) :: SP
REAL (EB), POINTER, DIMENSION (:, :, :) :: X, Y
INTEGER :: ILEVEL, I, K, IM, NM, IERR
REAL(EB):: TNOW_SCALPROD2D
 
INCLUDE 'mpif.h'
 
TNOW_SCALPROD2D = SECOND()
 
SL => S%SLEVEL(ILEVEL)
 
IF (SCARC_DEBUG >= 2) THEN
   WRITE (9,*) 'SCAL_PROD2D'
   WRITE (9,*) 'NMESHES=', NMESHES
   WRITE (9,*) 'IBAR=', SL%IBAR
   WRITE (9,*) 'KBAR=', SL%KBAR
   CALL SCARC_COMPARE_SINGLE (X, 'SARC', 1, 'X      ',0)
   CALL SCARC_COMPARE_SINGLE (Y, 'SARC', 1, 'Y      ',0)
ENDIF
 
 
SP = 0.0_EB
DO K = 1, SL%KBAR
   DO I = 1, SL%IBAR
      SP = SP + X (I, 1, K) * Y (I, 1, K)
      IF (SCARC_DEBUG >= 2) write (9, 1001) I, K, X (I, 1, K), Y (I, 1, K), SP
   ENDDO
ENDDO
 
IF (SCARC_DEBUG >= 4) write (9,*) 'SP=', SP
 
IF (NMESHES > 1) THEN
 
   S%SP_LOCAL0 = 0.0_EB
   S%SP_LOCAL0 (NM) = SP
 
   S%SP_GLOBAL0 = 0.0_EB
   CALL MPI_ALLREDUCE (S%SP_LOCAL0(1), S%SP_GLOBAL0(1), NMESHES, &
    MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
 
   IF (SCARC_DEBUG >= 2) write (9,*) 'S%SP_LOCAL0='
   IF (SCARC_DEBUG >= 2) write (9, '(f22.16)') (S%SP_LOCAL0(I), I=1, &
    NMESHES)
   IF (SCARC_DEBUG >= 2) write (9,*) 'S%SP_GLOBAL0='
   IF (SCARC_DEBUG >= 2) write (9, '(f22.16)') (S%SP_GLOBAL0(I), I=1, &
    NMESHES)
 
   SP = 0.0_EB
   DO IM = 1, NMESHES
      SP = SP + S%SP_GLOBAL0(IM)
   ENDDO
 
ENDIF
IF (SCARC_DEBUG >= 2) write (9,*) 'SCALPROD2D: SP=', SP

TUSED_SCARC(16,NM)=TUSED_SCARC(16,NM)+SECOND()-TNOW_SCALPROD2D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_SCALPROD2D
 
1001     FORMAT ('sp=sp+[', i4, ',1,', i4, ']:[', f22.16, '*', f22.16, ']:', f22.16)
END SUBROUTINE SCARC_SCALPROD2D
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Compute scalarproduct of vector X and vector Y
!!!
!!! SP = (X,Y)
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SCALPROD3D (X, Y, SP, ILEVEL, NM)
 
REAL (EB) :: SP
REAL (EB), POINTER, DIMENSION (:, :, :) :: X, Y
INTEGER :: ILEVEL, I, J, K, IM, NM, IERR
REAL(EB):: TNOW_SCALPROD3D
 
INCLUDE 'mpif.h'
 
TNOW_SCALPROD3D = SECOND()
 
SL => S%SLEVEL(ILEVEL)

IF (SCARC_DEBUG >= 2) THEN
   WRITE (9,*) 'SCAL_PROD3D'
   CALL SCARC_COMPARE_SINGLE (X, 'SARC', 1, 'X      ',0)
   CALL SCARC_COMPARE_SINGLE (Y, 'SARC', 1, 'Y      ',0)
ENDIF
 
 
SP = 0.0_EB
DO K = 1, SL%KBAR
   DO J = 1, SL%JBAR
      DO I = 1, SL%IBAR
         SP = SP + X (I, J, K) * Y (I, J, K)
         IF (SCARC_DEBUG >= 2) write (9, 1002) I, J, K, X (I, J, K) * Y (I, J, &
          K), SP
      ENDDO
   ENDDO
ENDDO
 
IF (SCARC_DEBUG > 2) write (9,*) 'SP=', SP
 
IF (NMESHES > 1) THEN
 
   S%SP_LOCAL0 = 0.0_EB
   S%SP_LOCAL0 (NM) = SP
 
   S%SP_GLOBAL0 = 0.0_EB
   CALL MPI_ALLREDUCE (S%SP_LOCAL0(1), S%SP_GLOBAL0(1), NMESHES, &
    MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
 
   IF (SCARC_DEBUG >= 2) write (9,*) 'S%SP_LOCAL0='
   IF (SCARC_DEBUG >= 2) write (9, '(f22.16)') (S%SP_LOCAL0(I), I=1, &
    NMESHES)
   IF (SCARC_DEBUG >= 2) write (9,*) 'S%SP_GLOBAL0='
   IF (SCARC_DEBUG >= 2) write (9, '(f22.16)') (S%SP_GLOBAL0(I), I=1, &
    NMESHES)
 
   SP = 0.0_EB
   DO IM = 1, NMESHES
      SP = SP + S%SP_GLOBAL0(IM)
   ENDDO
   IF (SCARC_DEBUG >= 2) write (9,*) 'SP=', SP
 
ENDIF

TUSED_SCARC(16,NM)=TUSED_SCARC(16,NM)+SECOND()-TNOW_SCALPROD3D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_SCALPROD3D
 
1002 FORMAT ('sp=sp+[', i4, ',', i4, ',', i4, ']:[', f12.6, ']:', f24.12)
 
END SUBROUTINE SCARC_SCALPROD3D
 
 
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Compute L2-norm of vector X:   SP = ||X||
!!! depending on ITYPE:
!!! locally: only local cells are considered
!!! globally: all cells (globally) are considered
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_L2NORM2D (X, SP, ILEVEL, NM, ITYPE)
 
REAL (EB) :: SP
REAL (EB), POINTER, DIMENSION (:, :, :) :: X
INTEGER :: NM, ITYPE
INTEGER :: ILEVEL, IM, IERR, I, K
REAL(EB):: TNOW_L2NORM2D
 
INCLUDE 'mpif.h'
 
TNOW_L2NORM2D = SECOND()
 
SL => S%SLEVEL(ILEVEL)

IF (SCARC_DEBUG >= 6) THEN
   WRITE (9,*) '========= L2NORM , ITYPE=', ITYPE
   WRITE (9,*) 'NMESHES=', NMESHES
   WRITE (9,*) 'NM=', NM
   CALL SCARC_COMPARE_SINGLE_LEV (X, 'SARC', 1, 'X     ',0, NM, ILEVEL)
ENDIF
 
 
!!! build local scalar product (x,x)
S%SP_LOCAL = 0.0_EB
DO K = 1, SL%KBAR
   DO I = 1, SL%IBAR
      S%SP_LOCAL = S%SP_LOCAL + X (I, 1, K) * X (I, 1, K)
      IF (SCARC_DEBUG >= 6) write (9, 1002) I, K, X (I, 1, K), X (I, 1, K), &
       S%SP_LOCAL
   ENDDO
ENDDO
 
IF (SCARC_DEBUG >= 2) write (9,*) 'S%SP_LOCAL=', S%SP_LOCAL
 
!!! scale with number of local cells
IF (ITYPE == NCOM_TYPE_LOCAL) THEN
 
   SP = Sqrt (S%SP_LOCAL/REAL(SL%NCELLS_LOCAL, EB))
 
   IF (SCARC_DEBUG >= 2) THEN
      WRITE (9,*) 'SL%NCELLS_LOCAL=', SL%NCELLS_LOCAL
      WRITE (9,*) 'S%SP_LOCAL=', S%SP_LOCAL
      WRITE (9,*) 'NMESHES=', NMESHES
   ENDIF
 
!!! sum up globally and scale with global number of nodes
ELSE IF (ITYPE == NCOM_TYPE_GLOBAL) THEN
 
   IF (NMESHES > 1) THEN
 
      S%SP_LOCAL0 = 0.0_EB
      S%SP_LOCAL0 (NM) = S%SP_LOCAL
 
      S%SP_GLOBAL0 = 0.0_EB
      CALL MPI_ALLREDUCE (S%SP_LOCAL0(1), S%SP_GLOBAL0(1), NMESHES, &
                          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
 
      IF (SCARC_DEBUG >= 2) write (9,*) 'S%SP_LOCALA=',S%SP_LOCAL/SL%NCELLS_LOCAL
      IF (SCARC_DEBUG >= 2) write (9,*) 'S%SP_LOCAL0='
      IF (SCARC_DEBUG >= 2) write (9, '(f22.16)') (S%SP_LOCAL0(I), I=1, NMESHES)
      IF (SCARC_DEBUG >= 2) write (9,*) 'S%SP_GLOBAL0='
      IF (SCARC_DEBUG >= 2) write (9, '(f22.16)') (S%SP_GLOBAL0(I), I=1, NMESHES)
 
      S%SP_GLOBAL = 0.0_EB
      DO IM = 1, NMESHES
         S%SP_GLOBAL = S%SP_GLOBAL + S%SP_GLOBAL0(IM)
         IF (SCARC_DEBUG >= 4) write (9,*) 'S%SP_GLOBAL=', S%SP_GLOBAL
      ENDDO
   ELSE
      S%SP_GLOBAL = S%SP_LOCAL
   ENDIF
 
   SP = Sqrt (S%SP_GLOBAL/REAL(SL%NCELLS_GLOBAL, EB))
 
   IF (SCARC_DEBUG >= 2) write (9,*) 'SL%NCELLS_GLOBAL=', SL%NCELLS_GLOBAL
   IF (SCARC_DEBUG >= 2) write (9,*) 'L2NORM2D  : SP=', SP
 
ELSE
 
   WRITE (*,*) 'Wrong type for SCARC_L2NORM2D ', ITYPE
 
ENDIF
TUSED_SCARC(17,NM)=TUSED_SCARC(17,NM)+SECOND()-TNOW_L2NORM2D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_L2NORM2D
 
1002     FORMAT ('sp=sp+[', i4, ',', i4, ']:[', f12.6, '*', f12.6, ']:', f12.6)
END SUBROUTINE SCARC_L2NORM2D
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Compute L2-norm of vector X:   SP = ||X||
!!! depending on ITYPE:
!!! locally: only local cells are considered
!!! globally: all cells (globally) are considered
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_L2NORM3D (X, SP, ILEVEL, NM, ITYPE)
 
REAL (EB) :: SP
REAL (EB), POINTER, DIMENSION (:, :, :) :: X
INTEGER :: NM, ITYPE
INTEGER :: ILEVEL, IM, IERR, I, J, K
REAL(EB):: TNOW_L2NORM3D
 
INCLUDE 'mpif.h'
 
TNOW_L2NORM3D = SECOND()
 
IF (SCARC_DEBUG >= 2) write (9,*) '========= L2NORM , ITYPE=', ITYPE
IF (SCARC_DEBUG >= 2) write (9,*) 'NMESHES=', NMESHES
 
SL => S%SLEVEL(ILEVEL)
 
!!! build local scalar product (x,x)
S%SP_LOCAL = 0.0_EB
DO K = 1, SL%KBAR
   DO J = 1, SL%JBAR
      DO I = 1, SL%IBAR
         S%SP_LOCAL = S%SP_LOCAL + X (I, J, K) * X (I, J, K)
!if (SCARC_DEBUG.ge.2) WRITE(9,1002) I,J,K,X(I,J,K)*X(I,J,K),SP
      ENDDO
   ENDDO
ENDDO
 
IF (SCARC_DEBUG >= 2) write (9,*) 'S%SP_LOCAL=', S%SP_LOCAL
 
!!! scale with number of local cells
IF (ITYPE == NCOM_TYPE_LOCAL) THEN
 
   SP = Sqrt (S%SP_LOCAL/REAL(SL%NCELLS_LOCAL, EB))
 
   IF (SCARC_DEBUG >= 2) write (9,*) 'SL%NCELLS_LOCAL=', SL%NCELLS_LOCAL
   IF (SCARC_DEBUG >= 2) write (9,*) 'S%SP_LOCAL=', S%SP_LOCAL
   IF (SCARC_DEBUG >= 2) write (9,*) 'NMESHES=', NMESHES
 
!!! sum up globally and scale with global number of nodes
ELSE IF (ITYPE == NCOM_TYPE_GLOBAL) THEN
 
   IF (NMESHES > 1) THEN
 
      S%SP_LOCAL0 = 0.0_EB
      S%SP_LOCAL0 (NM) = S%SP_LOCAL
 
      S%SP_GLOBAL0 = 0.0_EB
      CALL MPI_ALLREDUCE (S%SP_LOCAL0(1), S%SP_GLOBAL0(1), NMESHES, &
       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
 
      IF (SCARC_DEBUG >= 2) write (9,*) 'S%SP_LOCAL0='
      IF (SCARC_DEBUG >= 2) write (9, '(f22.16)') (S%SP_LOCAL0(I), I=1, &
       NMESHES)
      IF (SCARC_DEBUG >= 2) write (9,*) 'S%SP_GLOBAL0='
      IF (SCARC_DEBUG >= 2) write (9, '(f22.16)') (S%SP_GLOBAL0(I), I=1, &
       NMESHES)
 
      S%SP_GLOBAL = 0.0_EB
      DO IM = 1, NMESHES
         S%SP_GLOBAL = S%SP_GLOBAL + S%SP_GLOBAL0(IM)
         IF (SCARC_DEBUG >= 2) write (9,*) 'S%SP_GLOBAL=', S%SP_GLOBAL
      ENDDO
   ELSE
      S%SP_GLOBAL = S%SP_LOCAL
   ENDIF
 
   SP = Sqrt (S%SP_GLOBAL/REAL(SL%NCELLS_GLOBAL, EB))
 
   IF (SCARC_DEBUG >= 2) write (9,*) 'SL%NCELLS_GLOBAL=', SL%NCELLS_GLOBAL
   IF (SCARC_DEBUG >= 2) write (9,*) 'SP=', SP
 
ELSE
 
   WRITE (*,*) 'Wrong type for SCARC_L2NORM ', ITYPE
 
ENDIF
TUSED_SCARC(17,NM)=TUSED_SCARC(17,NM)+SECOND()-TNOW_L2NORM3D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_L2NORM3D
 
END SUBROUTINE SCARC_L2NORM3D
 
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Jacobi preconditioner
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PRECON_JACOBI2D (A, Y, ILEVEL, NM)
REAL (EB), POINTER, DIMENSION (:, :) :: A
REAL (EB), POINTER, DIMENSION (:, :, :) :: Y
REAL (EB) :: TNOW_JACOBI2D
INTEGER :: I, K, IC, ILEVEL, NM
 
SL => S%SLEVEL(ILEVEL)

TNOW_JACOBI2D = SECOND()
 
IF (SCARC_DEBUG >= 2) THEN
   WRITE (9,*) 'BEFORE PRECON_JACOBI2D: Y'
   WRITE (9,*) 'Y:'
   !CALL SCARC_COMPARE_SINGLE_LEV (Y, 'SARC', 1, 'Y     ',0, NM, ILEVEL)
ENDIF
 
 

DO K = 1, SL%KBAR
   DO I = 1, SL%IBAR
      IC = (K-1) * SL%IBAR + I
      Y (I, 1, K) = Y (I, 1, K) / A (IC, 1)
      IF (SCARC_DEBUG >= 2) write (9, '(a,i3,a,f8.2,a,f18.12)') 'A(', IC, ')=',&
        A (IC, 1), ': Y=', Y (I, 1, K)
   ENDDO
   IF (SCARC_DEBUG >= 2) write (9,*) '========================================&
  &======'
ENDDO
 
 
IF (SCARC_DEBUG >= 4) THEN
   WRITE (9,*)
   WRITE (9,*) 'AFTER  PRECON_JACOBI: Y'
   !CALL SCARC_COMPARE_SINGLE (Y, 'SARC', 1, 'Y     ',0)
ENDIF
TUSED_SCARC(18,NM)=TUSED_SCARC(18,NM)+SECOND()-TNOW_JACOBI2D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_JACOBI2D
 
END SUBROUTINE SCARC_PRECON_JACOBI2D
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Jacobi preconditioner
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PRECON_JACOBI3D (A, Y, ILEVEL, NM)
REAL (EB), POINTER, DIMENSION (:, :) :: A
REAL (EB), POINTER, DIMENSION (:, :, :) :: Y
REAL (EB) :: YOLD
REAL (EB) :: TNOW_JACOBI3D
INTEGER :: I, J, K, IC, ILEVEL, NM
 
SL => S%SLEVEL(ILEVEL)

TNOW_JACOBI3D = SECOND()
 
IF (SCARC_DEBUG >= 2) THEN
   WRITE (9,*) 'BEFORE PRECON_JACOBI3D: Y'
   !CALL SCARC_COMPARE_SINGLE (Y, 'SARC', 1, 'Y     ',0)
ENDIF
 
 
DO K = 1, SL%KBAR
   DO J = 1, SL%JBAR
      DO I = 1, SL%IBAR
         IC = (K-1) * SL%IBAR * SL%JBAR + (J-1) * SL%IBAR + I
         YOLD = Y (I, J, K)
         IF (TWO_D) THEN
            Y (I, J, K) = Y (I, J, K) / A (IC, 3)
         ELSE
            Y (I, J, K) = Y (I, J, K) / A (IC, 4)
            IF (SCARC_DEBUG >= 2) write (9, '(a,i3,a,f8.2,a,f18.12)') 'A(', &
             IC, ')=', A (IC, 4), ': Y=', Y (I, J, K)
         ENDIF
      ENDDO
      IF (SCARC_DEBUG >= 2) write (9,*) '-------------'
   ENDDO
   IF (SCARC_DEBUG >= 2) write (9,*) '========================================&
  &======'
ENDDO
 
IF (SCARC_DEBUG >= 2) THEN
   WRITE (9,*) 'AFTER  PRECON_JACOBI3D: Y'
   !CALL SCARC_COMPARE_SINGLE (Y, 'SARC', 1, 'Y     ',0)
ENDIF
TUSED_SCARC(18,NM)=TUSED_SCARC(18,NM)+SECOND()-TNOW_JACOBI3D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_JACOBI3D
 
END SUBROUTINE SCARC_PRECON_JACOBI3D
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! GS preconditioner
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PRECON_GS2D (A, Y, ILEVEL, NM)
REAL (EB), POINTER, DIMENSION (:, :) :: A
REAL (EB), POINTER, DIMENSION (:, :, :) :: Y
REAL (EB) :: AUX, OMEGA
INTEGER :: I, J, K, IC, ILEVEL, NM
REAL(EB):: TNOW_GS2D

TNOW_GS2D = SECOND()

SL => S%SLEVEL(ILEVEL)

IF (SCARC_DEBUG >= 2) THEN
   WRITE (9,*) 'BEFORE PRECON_GS2D: Y'
   !CALL SCARC_COMPARE_SINGLE (Y, 'SARC', 1, 'Y     ',0)
ENDIF
 
 
OMEGA = 1.0_EB
 
DO K = 1, SL%KBAR
   DO I = 1, SL%IBAR
      IC = (K-1) * SL%IBAR + I
 
      AUX = AUX + A (IC, 1) * Y (I, 1, K-1) + A (IC, 2) * Y (I-1, 1, K)
      Y (I, 1, K) = (Y(I, 1, K)-AUX*OMEGA) / A (IC, 3)
      IF (SCARC_DEBUG >= 2) write (9, '(a,i3,a,f8.2,a,f18.12)') 'A(', IC, ')=',&
        A (IC, 3), ': Y=', Y (I, J, K)
   ENDDO
   IF (SCARC_DEBUG >= 2) write (9,*) '========================================&
  &======'
ENDDO
 
 
IF (SCARC_DEBUG >= 2) THEN
   WRITE (9,*) 'AFTER  PRECON_GS: Y'
   !CALL SCARC_COMPARE_SINGLE (Y, 'SARC', 1, 'Y     ',0)
ENDIF
TUSED_SCARC(19,NM)=TUSED_SCARC(19,NM)+SECOND()-TNOW_GS2D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_GS2D
 
END SUBROUTINE SCARC_PRECON_GS2D
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! GS preconditioner
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PRECON_GS3D (A, Y, ILEVEL, NM)
REAL (EB), POINTER, DIMENSION (:, :) :: A
REAL (EB), POINTER, DIMENSION (:, :, :) :: Y
REAL (EB) :: AUX, OMEGA
INTEGER :: I, J, K, IC, ILEVEL, NM
REAL (EB) :: TNOW_GS3D

 
TNOW_GS3D = SECOND()
 
SL => S%SLEVEL(ILEVEL)

IF (SCARC_DEBUG >= 2) THEN
   WRITE (9,*) 'BEFORE PRECON_GS3D: Y'
   !CALL SCARC_COMPARE_SINGLE (Y, 'SARC', 1, 'Y     ',0)
ENDIF
 
OMEGA = 1.0_EB
 
DO K = 1, SL%KBAR
   DO J = 1, SL%JBAR
      DO I = 1, SL%IBAR
         IC = (K-1) * SL%IBAR * SL%JBAR + (J-1) * SL%IBAR + I
         AUX = 0.0_EB
         AUX = A (IC, 1) * Y (I, J, K-1) + A (IC, 2) * Y (I, J-1, K) + A (IC, 3) &
          * Y (I-1, J, K)
         Y (I, J, K) = (Y(I, J, K)-OMEGA*AUX) / A (IC, 4)
         IF (SCARC_DEBUG >= 2) write (9, '(a,i3,a,f8.2,a,f18.12,a,f18.12)') 'A&
        &(', IC, ')=', A (IC, 4), ': AUX=', AUX, ': Y=', Y (I, J, K)
      ENDDO
      IF (SCARC_DEBUG >= 2) write (9,*) '-------------'
   ENDDO
ENDDO
 
IF (SCARC_DEBUG >= 2) THEN
   WRITE (9,*) 'AFTER  PRECON_GS3D: Y'
   !CALL SCARC_COMPARE_SINGLE (Y, 'SARC', 1, 'Y     ',0)
ENDIF
 
TUSED_SCARC(19,NM)=TUSED_SCARC(19,NM)+SECOND()-TNOW_GS3D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_GS3D
END SUBROUTINE SCARC_PRECON_GS3D
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! SSOR preconditioner
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PRECON_SSOR2D (A, Y, ILEVEL, NM)
REAL (EB), POINTER, DIMENSION (:, :) :: A
REAL (EB), POINTER, DIMENSION (:, :, :) :: Y
REAL (EB) :: AUX, OMEGA
INTEGER :: I, K, IC, ILEVEL, NM
REAL (EB) :: TNOW_SSOR2D
 

TNOW_SSOR2D = SECOND()

SL => S%SLEVEL(ILEVEL)
 
IF (SCARC_DEBUG >= 4) THEN
   WRITE (9,*) 'BEFORE PRECON_SSOR2D: Y'
   CALL SCARC_COMPARE_SINGLE (Y, 'SARC', 1, 'Y     ',0)
ENDIF
 
OMEGA = 1.5_EB
 
DO K = 1, SL%KBAR
   DO I = 1, SL%IBAR
      IC = (K-1) * SL%IBAR + I
      AUX = A (IC, 2) * Y (I, 1, K-1) + A (IC, 3) * Y (I-1, 1, K)
      Y (I, 1, K) = (Y(I, 1, K)-AUX*OMEGA) / A (IC, 1)
      IF (SCARC_DEBUG >= 2) write (9, '(a,i3,a,f8.2,a,f18.12,a,f18.12)') 'A(', &
       IC, ')=', A (IC, 1), ': AUX=', AUX, ': Y=', Y (I, 1, K)
   ENDDO
   IF (SCARC_DEBUG >= 2) write (9,*) '========================================&
  &======'
ENDDO
 
DO K = SL%KBAR, 1, - 1
   DO I = SL%IBAR, 1, - 1
      IC = (K-1) * SL%IBAR + I
      AUX = A (IC, 5) * Y (I, 1, K+1) + A (IC, 4) * Y (I+1, 1, K)
      Y (I, 1, K) = Y (I, 1, K) - AUX * OMEGA / A (IC, 1)
      IF (SCARC_DEBUG >= 2) write (9, '(a,i3,a,f8.2,a,f18.12,a,f18.12)') 'A(', &
       IC, ')=', A (IC, 1), ': AUX=', AUX, ': Y=', Y (I, 1, K)
   ENDDO
   IF (SCARC_DEBUG >= 2) write (9,*) '========================================&
  &======'
ENDDO
 
IF (SCARC_DEBUG >= 4) THEN
   WRITE (9,*) 'AFTER  PRECON_SSOR: Y'
   CALL SCARC_COMPARE_SINGLE (Y, 'SARC', 1, 'Y     ',0)
ENDIF
TUSED_SCARC(20,NM)=TUSED_SCARC(20,NM)+SECOND()-TNOW_SSOR2D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_SSOR2D
 
END SUBROUTINE SCARC_PRECON_SSOR2D
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! SSOR preconditioner
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PRECON_SSOR3D (A, Y, ILEVEL, NM)
REAL (EB), POINTER, DIMENSION (:, :) :: A
REAL (EB), POINTER, DIMENSION (:, :, :) :: Y
REAL (EB) :: AUX, OMEGA
INTEGER :: I, J, K, IC, ILEVEL, NM
REAL (EB) :: TNOW_SSOR3D
 
TNOW_SSOR3D = SECOND()

SL => S%SLEVEL(ILEVEL)
 
IF (SCARC_DEBUG >= 2) THEN
   WRITE (9,*) 'BEFORE PRECON_SSOR: Y'
   CALL SCARC_COMPARE_SINGLE (Y, 'SARC', 1, 'Y     ',0)
ENDIF
 
 
OMEGA = 1.6_EB
 
DO K = 1, SL%KBAR
   DO J = 1, SL%JBAR
      DO I = 1, SL%IBAR
         IC = (K-1) * SL%IBAR * SL%JBAR + (J-1) * SL%IBAR + I
 
         AUX = A (IC, 1) * Y (I, J, K-1) + A (IC, 2) * Y (I, J-1, K) + A (IC, 3) &
          * Y (I-1, J, K)
         Y (I, J, K) = (Y(I, J, K)-AUX*OMEGA) / A (IC, 4)
         IF (SCARC_DEBUG >= 2) write (9, '(a,i3,a,f8.2,a,f18.12,a,f18.12)') 'A&
        &(', IC, ')=', A (IC, 4), ': AUX=', AUX, ': Y=', Y (I, J, K)
      ENDDO
      IF (SCARC_DEBUG >= 2) write (9,*) '-------------'
   ENDDO
   IF (SCARC_DEBUG >= 2) write (9,*) '========================================&
  &======'
ENDDO
 
DO K = SL%KBAR, 1, - 1
   DO J = SL%JBAR, 1, - 1
      DO I = SL%IBAR, 1, - 1
 
         IC = (K-1) * SL%IBAR * SL%JBAR + (J-1) * SL%IBAR + I
 
         AUX = A (IC, 7) * Y (I, J, K+1) + A (IC, 6) * Y (I, J+1, K) + A (IC, 5) &
          * Y (I+1, J, K)
 
         Y (I, J, K) = Y (I, J, K) - AUX * OMEGA / A (IC, 4)
         IF (SCARC_DEBUG >= 2) write (9, '(a,i3,a,f8.2,a,f18.12,a,f18.12)') 'A&
        &(', IC, ')=', A (IC, 4), ': AUX=', AUX, ': Y=', Y (I, J, K)
      ENDDO
      IF (SCARC_DEBUG >= 2) write (9,*) '-------------'
   ENDDO
   IF (SCARC_DEBUG >= 2) write (9,*) '========================================&
  &======'
ENDDO
 
IF (SCARC_DEBUG >= 2) THEN
   WRITE (9,*) 'AFTER  PRECON_SSOR3D: Y'
   CALL SCARC_COMPARE_SINGLE (Y, 'SARC', 1, 'Y     ',0)
ENDIF
 
TUSED_SCARC(20,NM)=TUSED_SCARC(20,NM)+SECOND()-TNOW_SSOR3D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_SSOR3D
 
END SUBROUTINE SCARC_PRECON_SSOR3D
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Initialize GSTRIX preconditioner
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PRECON_GSTRIX2D_INIT(ILEVEL, NM)

INTEGER :: NM,ILEVEL, IC, NC, IZERO
REAL(EB):: TNOW_GSTRIX2D_INIT
 
TNOW_GSTRIX2D_INIT = SECOND()
 
SL => S%SLEVEL(ILEVEL)

ALLOCATE (SL%DD(1:SL%NCELLS_LOCAL), STAT=IZERO)
CALL CHKMEMERR ('SCARC_PRECON_GSTRIX2D_INIT', 'SL%DD', IZERO)

ALLOCATE (SL%DU(1:SL%NCELLS_LOCAL), STAT=IZERO)
CALL CHKMEMERR ('SCARC_PRECON_GSTRIX2D_INIT', 'SL%DU', IZERO)

ALLOCATE (SL%DL(1:SL%NCELLS_LOCAL), STAT=IZERO)
CALL CHKMEMERR ('SCARC_PRECON_GSTRIX2D_INIT', 'SL%DL', IZERO)

ALLOCATE (SL%LD(1:SL%NCELLS_LOCAL), STAT=IZERO)
CALL CHKMEMERR ('SCARC_PRECON_GSTRIX2D_INIT', 'SL%LD', IZERO)

ALLOCATE (SL%DAUX(1:SL%NCELLS_LOCAL), STAT=IZERO)
CALL CHKMEMERR ('SCARC_PRECON_GSTRIX2D_INIT', 'SL%DAUX', IZERO)


NC = SL%NCELLS_LOCAL
 
!!! copy diagonal
DO IC = 1, NC
   SL%DD(IC) = SL%AG(IC, 1)         ! diagonal
   SL%LD(IC) = SL%AG(IC, 2)         ! lower diagonal in z-direction
   SL%DL(IC) = SL%AG(IC, 3)         ! lower diagonal in x-direction
   SL%DU(IC) = SL%AG(IC, 4)         ! upper diagonal in x-direction
ENDDO
 
 
!!! perform LU-factorization of matrix AG
DO IC = 2, NC - 1
   SL%DL (IC) = SL%DL(IC) / SL%DD(IC-1)
   SL%DD (IC) = SL%DD(IC) - SL%DL(IC) * SL%DU(IC-1)
ENDDO
 
SL%DL (NC) = SL%DL(NC) / SL%DD(NC-1)
SL%DD (NC) = SL%DD(NC) - SL%DL(NC) * SL%DU(NC-1)
 
!!! replace diagonal values diag(i) by 1/diag(i)
!!! multiply DU with 1/diag(i)
DO IC = 1, NC
   SL%DD (IC) = 1.0_EB / SL%DD(IC)
   SL%DU (IC) = SL%DD(IC) * SL%DU(IC)
   IF (SCARC_DEBUG>=2) THEN
      write(9,*) 'GSTRIX_INIT:',ILEVEL,': SL%DU(',IC,')=',SL%DU(IC),': SL%DD(',IC,')=',SL%DD(IC)
   ENDIF
ENDDO

TUSED_SCARC(32,NM)=TUSED_SCARC(32,NM)+SECOND()-TNOW_GSTRIX2D_INIT
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_GSTRIX2D_INIT
END SUBROUTINE SCARC_PRECON_GSTRIX2D_INIT
 
 
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! GSTRIX2D preconditioner
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PRECON_GSTRIX2D (DD, DL, DU, LD, DAUX, DY, ILEVEL, NM)
REAL (EB), POINTER, DIMENSION (:) :: DD, DL, DU, LD, DAUX
REAL (EB), POINTER, DIMENSION (:, :, :) :: DY
INTEGER :: I, K, IC, ILEVEL, NM, N, M
REAL (EB) :: TNOW_GSTRIX2D
 
TNOW_GSTRIX2D = SECOND()

SL => S%SLEVEL(ILEVEL)
 
N =SL%NCELLS_LOCAL
M =INT(SQRT(REAL(N,EB)))

IF (SCARC_DEBUG>=2) THEN
   write(9,*) '===================== STARTING PRECON_GSTRIX2D'
   write(9,*) 'N=',N
   write(9,*) 'M=',M
   WRITE (9,*) 'BEFORE PRECON_GSTRIX: Y'
   CALL SCARC_COMPARE_SINGLE (DY, 'SARC', 1, 'Y     ',0)
ENDIF

DO K=1,SL%KBAR
   DO I=1,SL%IBAR
      IC = (K-1)*SL%IBAR + I
      DAUX(IC)=DY(I,1,K)
write(9,*) 'DAUX(',IC,')=',DAUX(IC)
   ENDDO
ENDDO

DO I=N-1,N-M+1,-1
   DAUX(I) = DAUX(I)-DU(I)*DAUX(I+1)
ENDDO

DO I=N,N-M+1,-1
   DAUX(I) = DAUX(I)*DD(I)
ENDDO

DO I=N-M+2,N
   DAUX(I) = DAUX(I)-DL(I)*DAUX(I-1)
ENDDO
  

DO NM=M-2,0,-1
   DO I=M,1,-1
      DAUX(NM*M+I) = DAUX(NM*M+I) - LD(NM*M+I)*DAUX((NM+1)*M+I) 
   ENDDO

   DO I=M-1,1,-1
      DAUX(NM*M+I) = DAUX(NM*M+I)-DU(NM*M+I)*DAUX(NM*M+I+1)
   ENDDO

   DO I=M,1,-1
      DAUX(NM*M+I) = DAUX(NM*M+I)*DD(NM*M+I)
   ENDDO
 
   DO I=2,M
      DAUX(NM*M+I) = DAUX(NM*M+I)-DL(NM*M+I)*DAUX(NM*M+I-1)
   ENDDO
  
ENDDO

DO K=1,SL%KBAR
   DO I=1,SL%IBAR
      IC = (K-1)*SL%IBAR + I
      DY(I,1,K)=DAUX(IC)
   ENDDO
ENDDO

IF (SCARC_DEBUG >= 2) THEN
   WRITE (9,*) 'AFTER  PRECON_GSTRIX: Y'
   CALL SCARC_COMPARE_SINGLE (DY, 'SARC', 1, 'Y     ',0)
ENDIF
 
TUSED_SCARC(31,NM)=TUSED_SCARC(20,NM)+SECOND()-TNOW_GSTRIX2D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_GSTRIX2D
 
END SUBROUTINE SCARC_PRECON_GSTRIX2D

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! GSADI preconditioner
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!SUBROUTINE SCARC_PRECON_GSTRIZ2D (A, Y, ITE, ILEVEL, NM)
!REAL (EB), POINTER, DIMENSION (:, :) :: A
!REAL (EB), POINTER, DIMENSION (:, :, :) :: Y
!INTEGER :: I, K, IC, ILEVEL, NM, ITE
!REAL (EB) :: TNOW_GSTRIZ2D
! 
!TNOW_GSTRIZ2D = SECOND()
!
!! NOT FINISHED 
!
!TUSED_SCARC(31,NM)=TUSED_SCARC(20,NM)+SECOND()-TNOW_GSTRIZ2D
!TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_GSTRIZ2D
! 
!END SUBROUTINE SCARC_PRECON_GSTRIZ2D
 
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
!IF (SCARC_DEBUG >= 2) THEN
!   WRITE (9,*) 'BEFORE PRECON_GSADI: Y'
!   CALL SCARC_COMPARE_SINGLE (Y, 'SARC', 1, 'Y     ',0)
!ENDIF
! 
!!IF (MOD(ITE,2)==1) THEN
!!  CALL SCARC_PRECON_GSTRIX2D(A, Y, ILEVEL, NM)
!!ELSE
!!  CALL SCARC_PRECON_GSTRIX2D(A, Y, ILEVEL, NM)
!!ENDIF
!
!IF (SCARC_DEBUG >= 2) THEN
!   WRITE (9,*) 'AFTER  PRECON_GSADI: Y'
!   CALL SCARC_COMPARE_SINGLE (Y, 'SARC', 1, 'Y     ',0)
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
INTEGER :: IMIN, IMAX, JMIN, JMAX, KMIN, KMAX
INTEGER :: IERR
INTEGER :: ILEN_FACE, ILEN_EDGE, ILEN_DIAG
INTEGER :: TAG_FACE, TAG_EDGE, TAG_DIAG
REAL(EB):: TNOW_RECEIVE

TYPE (SCARC_TYPE),  POINTER ::  SNM
TYPE (OSCARC_TYPE), POINTER :: OSNM
TYPE (SCARC_LEVEL_TYPE),  POINTER ::  SNML
TYPE (OSCARC_LEVEL_TYPE), POINTER :: OSNML

INCLUDE 'mpif.h'

TNOW_RECEIVE = SECOND()

RECEIVE_MESH_LOOP: DO NM=1,NMESHES

IF (SCARC_DEBUG>=6) THEN
   WRITE(9,*) '======================================================='
   WRITE(9,*) 'SCARC_RECEIVE:NREQ_FACE=',NREQ_FACE
   WRITE(9,*) 'SCARC_RECEIVE:NREQ_EDGE=',NREQ_EDGE
   WRITE(9,*) 'SCARC_RECEIVE:NREQ_DIAG=',NREQ_DIAG
   WRITE(9,*) 'SCARC_RECEIVE: NM=',NM
   WRITE(9,*) 'SCARC_RECEIVE: NMESHES=',NMESHES
   WRITE(9,*) 'SCARC_RECEIVE: ILEVEL=',ILEVEL
   WRITE(9,*) 'SCARC_RECEIVE: CODE=',CODE
   WRITE(9,*) 'SCARC_RECEIVE: MYID=',MYID
   WRITE(9,*) 'SCARC_RECEIVE: PROCESS(NM)=',PROCESS(NM)
ENDIF

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
    
IF (SCARC_DEBUG>=6) WRITE(9,*) 'IN RECEIVE_FACE_IF'

         OSNM  => SCARC(NM)%OSCARC(NOM)                  ! corresponds to M3
         OSNML => SCARC(NM)%OSCARC(NOM)%SLEVEL(ILEVEL)   ! ... for the level 'ILEVEL'

         TAG_FACE = TAGS_FACE(NM,NOM)
    

         !!! Initialize the communication structures for receiving face data
         INIT_FACE_IF: IF (CODE==NCOM_INIT) THEN

            NREQ_FACE = NREQ_FACE+1
IF (SCARC_DEBUG>=6) THEN
   WRITE(9,*) 'SCARC_RECEIVE_FACE_INIT: NREQ_FACE=',NREQ_FACE
   WRITE(9,*) 'SCARC_RECEIVE_FACE_INIT: NW_FACE  =',17*OSNML%NW_FACE
ENDIF
            CALL MPI_IRECV(OSNML%IJKW_FACE(1,1),17*OSNML%NW_FACE,MPI_INTEGER,SNODE, &
                           TAG_FACE,MPI_COMM_WORLD,REQ_FACE(NREQ_FACE),IERR)

            IF (SNML%NIC_FACE(NM, NOM) > 0) THEN
               ILEN_FACE=(MAX(SNML%NIC_FACE(NM, NOM), SNML%NIC_FACE(NOM, NM))+2)*2+1
               ALLOCATE (OSNM%RECV_FACE(ILEN_FACE))
               OSNM%RECV_FACE = 0.0_EB
            ENDIF
IF (SCARC_DEBUG>=6) WRITE(9,*) 'SCARC_RECEIVE_FACE_INIT: ilen_face  =',ilen_face

         ENDIF INIT_FACE_IF
   
         !!! Perform exchange for matrix-vector multiplication (only face exchange needed)
         MATV_FACE_IF: IF (CODE==NCOM_MATV) THEN

            NREQ_FACE = NREQ_FACE+1
IF (SCARC_DEBUG>=6) WRITE(9,*) 'SCARC_RECEIVE_FACE_MATV: LENGTH  =',OSNM%RECV_FACE
            CALL MPI_IRECV(OSNM%RECV_FACE(1),SIZE(OSNM%RECV_FACE),MPI_DOUBLE_PRECISION,&
                        SNODE,TAG_FACE,MPI_COMM_WORLD,REQ_FACE(NREQ_FACE),IERR)

         ENDIF MATV_FACE_IF
   
         !!! Perform full exchange including edge (3D)
         FULL_FACE_IF: IF (CODE==NCOM_FULL) THEN

            NREQ_FACE = NREQ_FACE+1
IF (SCARC_DEBUG>=6) WRITE(9,*) 'SCARC_RECEIVE_FACE_FULL: LENGTH  =',SIZE(OSNM%RECV_FACE)
            CALL MPI_IRECV(OSNM%RECV_FACE(1),SIZE(OSNM%RECV_FACE),MPI_DOUBLE_PRECISION,&
                           SNODE,TAG_FACE,MPI_COMM_WORLD,REQ_FACE(NREQ_FACE),IERR)

         ENDIF FULL_FACE_IF

      ENDIF RECEIVE_FACE_IF

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! EDGE communication
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      RECEIVE_EDGE_IF: IF (BEDGE.AND.SNML%NIC_EDGE(NM,NOM)/=0 .AND. SNML%NIC_EDGE(NOM,NM)/=0) THEN
    
         OSNM  => SCARC(NM)%OSCARC(NOM)                  ! corresponds to M3
         OSNML => SCARC(NM)%OSCARC(NOM)%SLEVEL(ILEVEL)   ! ... for the level 'ILEVEL'

         TAG_EDGE = TAGS_EDGE(NM,NOM)
    
         IMIN=SNML%I_MIN_EDGE(NM,NOM)
         IMAX=SNML%I_MAX_EDGE(NM,NOM)
         JMIN=SNML%J_MIN_EDGE(NM,NOM)
         JMAX=SNML%J_MAX_EDGE(NM,NOM)
         KMIN=SNML%K_MIN_EDGE(NM,NOM)
         KMAX=SNML%K_MAX_EDGE(NM,NOM)

         ILEN_EDGE=(IMAX-IMIN+1)*(JMAX-JMIN+1)*(KMAX-KMIN+1)

      !!! Initialize the communication structures for receiving face data
         INIT_EDGE_IF: IF (CODE==NCOM_INIT) THEN

            NREQ_EDGE = NREQ_EDGE+1
            CALL MPI_IRECV(OSNML%IJKW_EDGE(1,1),17*OSNML%NW_EDGE,MPI_INTEGER,SNODE, &
                           TAG_EDGE,MPI_COMM_WORLD,REQ_EDGE(NREQ_EDGE),IERR)

            IF (SNML%NIC_EDGE(NM, NOM) > 0) THEN
               ALLOCATE (OSNM%RECV_EDGE(ILEN_EDGE))
               OSNM%RECV_EDGE = 0.0_EB
            ENDIF

         ENDIF INIT_EDGE_IF
   
   
      !!! Perform full exchange including edge (3D)
         FULL_EDGE_IF: IF (CODE==NCOM_FULL) THEN

            NREQ_EDGE = NREQ_EDGE+1
            CALL MPI_IRECV(OSNM%RECV_EDGE(1),SIZE(OSNM%RECV_EDGE),MPI_DOUBLE_PRECISION,&
                           SNODE,TAG_EDGE,MPI_COMM_WORLD,REQ_EDGE(NREQ_EDGE),IERR)

         ENDIF FULL_EDGE_IF

      ENDIF RECEIVE_EDGE_IF

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! DIAG communication
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      RECEIVE_DIAG_IF: IF (BDIAG.AND.SNML%NIC_DIAG(NM,NOM)/=0 .AND. SNML%NIC_DIAG(NOM,NM)/=0) THEN
    
IF (SCARC_DEBUG>=6) WRITE(9,*) 'IN RECEIVE_DIAG_IF'

         OSNM  => SCARC(NM)%OSCARC(NOM)                  ! corresponds to M3
         OSNML => SCARC(NM)%OSCARC(NOM)%SLEVEL(ILEVEL)   ! ... for the level 'ILEVEL'

         TAG_DIAG = TAGS_DIAG(NM,NOM)
    
         IMIN=SNML%I_MIN_DIAG(NM,NOM)
         IMAX=SNML%I_MAX_DIAG(NM,NOM)
         JMIN=SNML%J_MIN_DIAG(NM,NOM)
         JMAX=SNML%J_MAX_DIAG(NM,NOM)
         KMIN=SNML%K_MIN_DIAG(NM,NOM)
         KMAX=SNML%K_MAX_DIAG(NM,NOM)

         ILEN_DIAG=(IMAX-IMIN+1)*(JMAX-JMIN+1)*(KMAX-KMIN+1)

IF (SCARC_DEBUG>=6) THEN
   WRITE(9,*) 'SCARC_RECEIVE_DIAG: IMIN=',IMIN
   WRITE(9,*) 'SCARC_RECEIVE_DIAG: IMAX=',IMAX
   WRITE(9,*) 'SCARC_RECEIVE_DIAG: JMIN=',JMIN
   WRITE(9,*) 'SCARC_RECEIVE_DIAG: JMAX=',JMAX
   WRITE(9,*) 'SCARC_RECEIVE_DIAG: KMIN=',KMIN
   WRITE(9,*) 'SCARC_RECEIVE_DIAG: KMAX=',KMAX
   WRITE(9,*) 'SCARC_RECEIVE_DIAG: ILEN_DIAG=',ILEN_DIAG
ENDIF

      !!! Initialize the communication structures for receiving face data
         INIT_DIAG_IF: IF (CODE==NCOM_INIT) THEN

            NREQ_DIAG = NREQ_DIAG+1
IF (SCARC_DEBUG>=6) THEN
   WRITE(9,*) 'SCARC_RECEIVE_DIAG_INIT: NREQ_diag=',NREQ_diag
   WRITE(9,*) 'SCARC_RECEIVE_DIAG_INIT: NW_FACE  =',17*OSNML%NW_diag
ENDIF
            CALL MPI_IRECV(OSNML%IJKW_DIAG(1,1),17*OSNML%NW_DIAG,MPI_INTEGER,SNODE, &
                           TAG_DIAG,MPI_COMM_WORLD,REQ_DIAG(NREQ_DIAG),IERR)

            IF (SNML%NIC_DIAG(NM, NOM) > 0) THEN
               ILEN_DIAG=(MAX(SNML%NIC_DIAG(NM, NOM), SNML%NIC_DIAG(NOM, NM))+2)*2+1   ! extended
IF (SCARC_DEBUG>=6) WRITE(9,*) 'SCARC_RECEIVE_DIAG_INIT: ilen_diag  =',ilen_diag
               ALLOCATE (OSNM%RECV_DIAG(ILEN_DIAG))
               OSNM%RECV_DIAG = 0.0_EB
            ENDIF

         ENDIF INIT_DIAG_IF
   
   
      !!! Perform full exchange including edge (3D)
         FULL_DIAG_IF: IF (CODE==NCOM_FULL) THEN

            NREQ_DIAG = NREQ_DIAG+1
IF (SCARC_DEBUG>=6) WRITE(9,*) 'SCARC_RECEIVE_DIAG_INIT: LENGTH  =',SIZE(OSNM%RECV_DIAG)
            CALL MPI_IRECV(OSNM%RECV_DIAG(1),SIZE(OSNM%RECV_DIAG),MPI_DOUBLE_PRECISION,&
                           SNODE,TAG_DIAG,MPI_COMM_WORLD,REQ_DIAG(NREQ_DIAG),IERR)

         ENDIF FULL_DIAG_IF

      ENDIF RECEIVE_DIAG_IF

   ENDDO RECEIVE_OMESH_LOOP
ENDDO RECEIVE_MESH_LOOP

IF (SCARC_DEBUG>=6) WRITE(9,*) 'SCARC_RECEIVE: READY'

TUSED_SCARC(22,MYID+1)=TUSED_SCARC(22,MYID+1)+SECOND()-TNOW_RECEIVE
TUSED_SCARC(0,MYID+1) =TUSED_SCARC(0,MYID+1) +SECOND()-TNOW_RECEIVE

END SUBROUTINE SCARC_RECEIVE
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  STILL EXPERIMENTAL !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_EXCHANGE (CODE, ILEVEL, IMV)

INTEGER :: NM, NOM, ILEVEL, CODE, IMV, ISUM
INTEGER :: IMIN, IMAX, JMIN, JMAX, KMIN, KMAX
INTEGER :: I, J, K, LL, IW, IWW
INTEGER :: IERR, II, JJ, KK
INTEGER :: ILEN_FACE, ILEN_EDGE, ILEN_DIAG
INTEGER :: TAG_FACE, TAG_EDGE, TAG_DIAG

REAL(EB):: ASUB, ZSUM,yold
REAL(EB):: TNOW_EXCHANGE


TYPE (SCARC_TYPE),  POINTER ::  SNM,   SNOM
TYPE (OSCARC_TYPE), POINTER :: OSNM, OSNOM
TYPE (SCARC_LEVEL_TYPE),  POINTER ::  SNML,  SNOML
TYPE (OSCARC_LEVEL_TYPE), POINTER :: OSNML, OSNOML

INCLUDE 'mpif.h'

TNOW_EXCHANGE = SECOND()

SNM => SCARC(NM)
IF (SCARC_DEBUG>=6) THEN
   WRITE(9,*) '======================== SCARC_EXCHANGE:'
   WRITE(9,*) 'SCARC_EXCHANGE   :NREQ_FACE=',NREQ_FACE
   WRITE(9,*) 'SCARC_EXCHANGE   :NREQ_EDGE=',NREQ_EDGE
   WRITE(9,*) 'SCARC_EXCHANGE   :NREQ_DIAG=',NREQ_DIAG
   WRITE(9,*) 'SCARC_EXCHANGE   :ILEVEL=',ILEVEL
   WRITE(9,*) 'SCARC_EXCHANGE   :CODE=',CODE
   WRITE(9,*) 'SCARC_EXCHANGE   :NREQ_FACE=',NREQ_FACE
   WRITE(9,*) 'SCARC_EXCHANGE   :NREQ_DIAG=',NREQ_DIAG
ENDIF
 
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

IF (SCARC_DEBUG>=6) WRITE(9,*) 'SCARC_EXCHANGE:   SNML zeigt auf SCARC(',NM,')%SLEVEL(',ILEVEL,')'

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! FACE communication
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FACE_IF: IF (SNML%NIC_FACE(NOM,NM)/=0 .AND. SNML%NIC_FACE(NM,NOM)/=0)  THEN
    
         OSNM  => SCARC(NM)%OSCARC(NOM)                 ! corresponds to M3
         OSNML => SCARC(NM)%OSCARC(NOM)%SLEVEL(ILEVEL)  ! ... for the level 'ILEVEL'

         TAG_FACE = TAGS_FACE(NM,NOM)
    
!         IMIN=SNML%I_MIN_FACE(NOM,NM)
!         IMAX=SNML%I_MAX_FACE(NOM,NM)
!         JMIN=SNML%J_MIN_FACE(NOM,NM)
!         JMAX=SNML%J_MAX_FACE(NOM,NM)
!         KMIN=SNML%K_MIN_FACE(NOM,NM)
!         KMAX=SNML%K_MAX_FACE(NOM,NM)
!
         !ILEN_FACE=(IMAX-IMIN+1)*(JMAX-JMIN+1)*(KMAX-KMIN+1)

!WRITE(9,*) 'SCARC_EXCHANGE:   I_MIN_FACE:'
!WRITE(9,'(2i4)') ((SNML%I_MIN_FACE(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
!WRITE(9,*) 'SCARC_EXCHANGE:   I_MAX:'
!WRITE(9,'(2i4)') ((SNML%I_MAX_FACE(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
!WRITE(9,*) 'SCARC_EXCHANGE:   J_MIN:'
!WRITE(9,'(2i4)') ((SNML%J_MIN_FACE(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
!WRITE(9,*) 'SCARC_EXCHANGE:   J_MAX:'
!WRITE(9,'(2i4)') ((SNML%J_MAX_FACE(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
!WRITE(9,*) 'SCARC_EXCHANGE:   K_MIN:'
!WRITE(9,'(2i4)') ((SNML%K_MIN_FACE(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
!WRITE(9,*) 'SCARC_EXCHANGE:   K_MAX:'
!WRITE(9,'(2i4)') ((SNML%K_MAX_FACE(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
!WRITE(9,*) 'SCARC_EXCHANGE:   IOR_FACE:'
!WRITE(9,'(2i4)') ((SNM%IOR_FACE(IM,JM),IM=1,NMESHES),JM=1,NMESHES)
!
!WRITE(9,*) 'SCARC_EXCHANGE: IMIN=',IMIN
!WRITE(9,*) 'SCARC_EXCHANGE: IMAX=',IMAX
!WRITE(9,*) 'SCARC_EXCHANGE: JMIN=',JMIN
!WRITE(9,*) 'SCARC_EXCHANGE: JMAX=',JMAX
!WRITE(9,*) 'SCARC_EXCHANGE: KMIN=',KMIN
!WRITE(9,*) 'SCARC_EXCHANGE: KMAX=',KMAX
!WRITE(9,*) 'SCARC_EXCHANGE: ILEN_FACE=',ILEN_FACE
!WRITE(9,*) 'SCARC_EXCHANGE :CODE=',CODE
!WRITE(9,*) 'SCARC_EXCHANGE :NM =',NM
!WRITE(9,*) 'SCARC_EXCHANGE :NOM=',NOM

         !!! Initialize the communication structures for sending data
         EXCHANGE_INIT_FACE_IF: IF (CODE==NCOM_INIT) THEN
    
            ! initialize communication structures for face exchange
            IF (RNODE/=SNODE) THEN
              NREQ_FACE = NREQ_FACE+1
IF (SCARC_DEBUG>=6) WRITE(9,*) 'SCARC_EXCHANGE_FACE_INIT: NREQ_FACE=',NREQ_FACE
IF (SCARC_DEBUG>=6) WRITE(9,*) 'SCARC_EXCHANGE_FACE_INIT: NW_FACE  =',17*SNML%NW_FACE
              CALL MPI_ISEND(SNML%IJKW_FACE(1,1),17*SNML%NW_FACE,MPI_INTEGER,SNODE, &
                             TAG_FACE,MPI_COMM_WORLD,REQ_FACE(NREQ_FACE),IERR)
               ILEN_FACE=(MAX(SNML%NIC_FACE(NM, NOM), SNML%NIC_FACE(NOM, NM))+2)*2+1   ! extended
IF (SCARC_DEBUG>=6) WRITE(9,*) 'SCARC_EXCHANGE_FACE_INIT: ilen_face  =',ilen_Face
               ALLOCATE (OSNM%SEND_FACE(ILEN_FACE))
               OSNM%SEND_FACE = 0.0_EB
            ENDIF

         ENDIF EXCHANGE_INIT_FACE_IF
   
      !!! Perform exchange for matrix-vector multiplication (only face exchange needed)
         EXCHANGE_MATV_FACE_IF: IF (CODE==NCOM_MATV) THEN
            IF (RNODE/=SNODE) THEN

   WRITE(9,*) 'IMV=',IMV
            LL = 0
            IWW = 0
            PACK_SEND_FACE0: DO IW=1,OSNML%NW_FACE
               IF (OSNML%IJKW_FACE(9,IW)/=NM) CYCLE PACK_SEND_FACE0
               DO KK=OSNML%IJKW_FACE(12,IW),OSNML%IJKW_FACE(15,IW)
                  DO JJ=OSNML%IJKW_FACE(11,IW),OSNML%IJKW_FACE(14,IW)
                     DO II=OSNML%IJKW_FACE(10,IW),OSNML%IJKW_FACE(13,IW)
                        IWW = IWW + 1
                        OSNM%SEND_FACE(LL+1) = REAL(IW,EB)
                        IF (IMV>=5) THEN
                           OSNM%SEND_FACE(LL+2) = SNML%D(II,JJ,KK)
                        ELSE IF (IMV==3) THEN
                           OSNM%SEND_FACE(LL+2) = SNML%X(II,JJ,KK)
                        ELSE
                           WRITE(*,*) 'WRONG MATVEC-TYPE'
                           stop
                        ENDIF
IF (SCARC_DEBUG>=6) WRITE(9,'(i4,a,i4,a,2e16.6,3i4)') IW,': SCARC_EXCHANGE_FACE_MATV: SENDING SEND_FACE(',LL+1,')=', &
                                OSNM%SEND_FACE(LL+1), OSNM%SEND_FACE(LL+2), II,JJ,KK
                        LL = LL+2
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO PACK_SEND_FACE0
            OSNM%SEND_FACE(IWW*2+1) = -999.0_EB
            NREQ_FACE=NREQ_FACE+1
IF (SCARC_DEBUG>=6) WRITE(9,*) 'SCARC_EXCHANGE_FACE_MATV: LENGTH=',IWW*2+1
            CALL MPI_ISEND(OSNM%SEND_FACE(1),IWW*2+1,MPI_DOUBLE_PRECISION,SNODE, &
                           TAG_FACE,MPI_COMM_WORLD,REQ_FACE(NREQ_FACE),IERR)

IF (SCARC_DEBUG>=6) WRITE(9,*) 'SEND_FACE=',OSNM%SEND_FACE(1)



            ELSE

!              IMIN = SNML%I_MIN_FACE(NOM,NM)
!              IMAX = SNML%I_MAX_FACE(NOM,NM)
!              JMIN = SNML%J_MIN_FACE(NOM,NM)
!              JMAX = SNML%J_MAX_FACE(NOM,NM)
!              KMIN = SNML%K_MIN_FACE(NOM,NM)
!              KMAX = SNML%K_MAX_FACE(NOM,NM)
!     
!              OSNOM=>SCARC(NOM)%OSCARC(NM)
!              OSNOM%D(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)= SNM%D(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
!
            ENDIF
         ENDIF  EXCHANGE_MATV_FACE_IF
   
      !!! Perform full exchange including edge (3D!) and vertex information
         EXCHANGE_FULL_FACE_IF: IF (CODE==NCOM_FULL) THEN
            IF (RNODE/=SNODE) THEN


            LL = 0
            IWW = 0
            PACK_SEND_FACE: DO IW=1,OSNML%NW_FACE
!WRITE(9,'(i4,a,17i4)') IW,': ',OSNML%IJKW_FACE(1:17,IW)
               IF (OSNML%IJKW_FACE(9,IW)/=NM) CYCLE PACK_SEND_FACE
               DO KK=OSNML%IJKW_FACE(12,IW),OSNML%IJKW_FACE(15,IW)
                  DO JJ=OSNML%IJKW_FACE(11,IW),OSNML%IJKW_FACE(14,IW)
                     DO II=OSNML%IJKW_FACE(10,IW),OSNML%IJKW_FACE(13,IW)
                        IWW = IWW + 1
                        OSNM%SEND_FACE(LL+1) = REAL(IW,EB)
                        OSNM%SEND_FACE(LL+2) = SNML%Z(II,JJ,KK)
IF (SCARC_DEBUG>=6) WRITE(9,'(i4,a,i4,a,2e16.6,3i4)') IW,':SCARC_EXCHANGE_FACE_FULL: SENDING SEND_FACE(',LL+1,')=', &
                                OSNM%SEND_FACE(LL+1), OSNM%SEND_FACE(LL+2), II,JJ,KK
!IF (SCARC_DEBUG>=6) THEN
!WRITE(9,*) 'SNML%Z=',LL+1,IW,SNML%Z(II,JJ,KK),OSNM%SEND_FACE(LL+2)
!ENDIF
                        LL = LL+2
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO PACK_SEND_FACE
            OSNM%SEND_FACE(IWW*2+1) = -999.0_EB
            NREQ_FACE=NREQ_FACE+1
IF (SCARC_DEBUG>=6) WRITE(9,*) 'SCARC_EXCHANGE_FACE_FULL: LENGTH=',IWW*2+1
            CALL MPI_ISEND(OSNM%SEND_FACE(1),IWW*2+1,MPI_DOUBLE_PRECISION,SNODE, &
                           TAG_FACE,MPI_COMM_WORLD,REQ_FACE(NREQ_FACE),IERR)


         ELSE
!            OSNOM=>SCARC(NOM)%OSCARC(NM)
!            WRITE(*,*)  'ACHTUNG: UNTERSCHEIDUNG face, edge vertex, IMIN ???'
!            WRITE(9,*)  'ACHTUNG: UNTERSCHEIDUNG face, edge vertex, IMIN ???'
!            OSNOM%D(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)= SNM%D(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
         ENDIF 

         ENDIF EXCHANGE_FULL_FACE_IF
      ENDIF FACE_IF

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! EDGE communication
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      EDGE_IF: IF (BEDGE.AND.SNML%NIC_EDGE(NOM,NM)/=0 .AND. SNML%NIC_EDGE(NM,NOM)/=0)  THEN
    
         OSNM  => SCARC(NM)%OSCARC(NOM)                 ! corresponds to M3
         OSNML => SCARC(NM)%OSCARC(NOM)%SLEVEL(ILEVEL)  ! ... for the level 'ILEVEL'

         TAG_EDGE = TAGS_EDGE(NM,NOM)
    
         !!! Initialize the communication structures for sending data
         EXCHANGE_INIT_EDGE_IF: IF (CODE==NCOM_INIT) THEN
    
            ! initialize communication structures for face exchange
            IF (RNODE/=SNODE) THEN
              NREQ_EDGE = NREQ_EDGE+1
              CALL MPI_ISEND(SNML%IJKW_EDGE(1,1),17*SNML%NW_EDGE,MPI_INTEGER,SNODE, &
                             TAG_EDGE,MPI_COMM_WORLD,REQ_EDGE(NREQ_EDGE),IERR)
            ENDIF
    
            IF (SNML%NIC_EDGE(NM, NOM) > 0 .AND. RNODE/=SNODE) THEN
               !ILEN_EDGE=(MAX(SNML%NIC_EDGE(NM, NOM), SNML%NIC_EDGE(NOM, NM))+2)*2+1   ! extended
               ALLOCATE (OSNM%SEND_EDGE(ILEN_EDGE))
               OSNM%SEND_EDGE = 0.0_EB
            ENDIF

         ENDIF EXCHANGE_INIT_EDGE_IF
   
      !!! Perform full exchange including edge (3D!) and vertex information
         EXCHANGE_FULL_EDGE_IF: IF (CODE==NCOM_FULL) THEN
            IF (RNODE/=SNODE) THEN

               LL = 0
               DO KK=KMIN,KMAX
                  DO JJ=JMIN,JMAX
                     DO II=IMIN,IMAX
                        OSNM%SEND_EDGE(LL+1) = SNML%Z(II,JJ,KK)
IF (SCARC_DEBUG>=6) WRITE(9,*) 'SCARC_EXCHANGE_FULL: SENDING SEND_EDGE(',LL+1,')=',OSNM%SEND_EDGE(LL+1), II, JJ, KK
                        LL = LL+1
                     ENDDO
                  ENDDO
               ENDDO
            NREQ_EDGE=NREQ_EDGE+1
            CALL MPI_ISEND(OSNM%SEND_EDGE(1),ILEN_EDGE,MPI_DOUBLE_PRECISION,SNODE, &
                           TAG_EDGE,MPI_COMM_WORLD,REQ_EDGE(NREQ_EDGE),IERR)
            ENDIF

         ELSE
!            OSNOM=>SCARC(NOM)%OSCARC(NM)
!            WRITE(*,*)  'ACHTUNG: UNTERSCHEIDUNG face, edge vertex, IMIN ???'
!            WRITE(9,*)  'ACHTUNG: UNTERSCHEIDUNG face, edge vertex, IMIN ???'
!            OSNOM%D(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)= SNM%D(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
         ENDIF EXCHANGE_FULL_EDGE_IF

      ENDIF EDGE_IF

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! DIAG
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DIAG_IF: IF (BDIAG.AND.SNML%NIC_DIAG(NOM,NM)/=0 .AND. SNML%NIC_DIAG(NM,NOM)/=0)  THEN
    
         OSNM  => SCARC(NM)%OSCARC(NOM)                 ! corresponds to M3
         OSNML => SCARC(NM)%OSCARC(NOM)%SLEVEL(ILEVEL)  ! ... for the level 'ILEVEL'

         TAG_DIAG = TAGS_DIAG(NM,NOM)
    
         !!! Initialize the communication structures for sending data
         EXCHANGE_INIT_DIAG_IF: IF (CODE==NCOM_INIT) THEN
    
            ! initialize communication structures for face exchange
            IF (RNODE/=SNODE) THEN
              NREQ_DIAG = NREQ_DIAG+1
              CALL MPI_ISEND(SNML%IJKW_DIAG(1,1),17*SNML%NW_DIAG,MPI_INTEGER,SNODE, &
                             TAG_DIAG,MPI_COMM_WORLD,REQ_DIAG(NREQ_DIAG),IERR)
            ENDIF
    
            IF (SNML%NIC_DIAG(NM, NOM) > 0 .AND. RNODE/=SNODE) THEN
               ILEN_DIAG=(MAX(SNML%NIC_DIAG(NM, NOM), SNML%NIC_DIAG(NOM, NM))+2)*2+1   ! extended
IF (SCARC_DEBUG>=6) WRITE(9,*) 'SCARC_EXCHANGE_DIAG_INIT: ILEN_DIAG=',ILEN_DIAG
               ALLOCATE (OSNM%SEND_DIAG(ILEN_DIAG))
               OSNM%SEND_DIAG = 0.0_EB
            ENDIF

         ENDIF EXCHANGE_INIT_DIAG_IF
   
      !!! Perform full exchange including edge (3D!) and vertex information
         EXCHANGE_FULL_DIAG_IF: IF (CODE==NCOM_FULL) THEN
            IF (RNODE/=SNODE) THEN

IF (SCARC_DEBUG>=6) WRITE(9,*) 'INIT_DIAG: NW-DIAG=',OSNML%NW_DIAG
               LL = 0
               IWW = 0
               PACK_SEND_DIAG: DO IW=1,OSNML%NW_DIAG
                  IF (OSNML%IJKW_DIAG(9,IW)/=NM) CYCLE PACK_SEND_DIAG
                  DO KK=OSNML%IJKW_DIAG(12,IW),OSNML%IJKW_DIAG(15,IW)
                     DO JJ=OSNML%IJKW_DIAG(11,IW),OSNML%IJKW_DIAG(14,IW)
                        DO II=OSNML%IJKW_DIAG(10,IW),OSNML%IJKW_DIAG(13,IW)
                           IWW = IWW + 1
                           OSNM%SEND_DIAG(LL+1) = REAL(IW,EB)
                           OSNM%SEND_DIAG(LL+2) = SNML%Z(II,JJ,KK)
IF (SCARC_DEBUG>=6) WRITE(9,*) 'SCARC_EXCHANGE_DIAG_FULL: SENDING SEND_DIAG(',LL+2,')=',OSNM%SEND_DIAG(LL+2), II, JJ, KK
                           LL = LL+2
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO PACK_SEND_DIAG
               OSNM%SEND_DIAG(IWW*2+1) = -999.0_EB
               NREQ_DIAG=NREQ_DIAG+1
IF (SCARC_DEBUG>=6) WRITE(9,*) 'SCARC_EXCHANGE_DIAG_FULL: LENGTH=',IWW*2+1
IF (SCARC_DEBUG>=6) WRITE(9,*) 'SCARC_EXCHANGE_DIAG_FULL: SIZE=',SIZE(OSNM%SEND_DIAG)
               CALL MPI_ISEND(OSNM%SEND_DIAG(1),IWW*2+1,MPI_DOUBLE_PRECISION,SNODE, &
                              TAG_DIAG,MPI_COMM_WORLD,REQ_DIAG(NREQ_DIAG),IERR)

            ENDIF

         ELSE
!            OSNOM=>SCARC(NOM)%OSCARC(NM)
!            WRITE(*,*)  'ACHTUNG: UNTERSCHEIDUNG face, edge vertex, IMIN ???'
!            WRITE(9,*)  'ACHTUNG: UNTERSCHEIDUNG face, edge vertex, IMIN ???'
!            OSNOM%D(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)= SNM%D(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
         ENDIF EXCHANGE_FULL_DIAG_IF

      ENDIF DIAG_IF
   
   ENDDO EXCHANGE_RECV_MESH_LOOP
ENDDO EXCHANGE_SEND_MESH_LOOP

IF (SCARC_DEBUG>=6) WRITE(9,*) 'SCARC_EXCHANGE: AFTER FIRST LOOP'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Information from Mesh NM is RECV'ed by Mesh NOM.  
!!! NOM is the receiver, NM is the sender.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (SCARC_DEBUG>=6) THEN
   WRITE(9,*) 'SCARC_EXCHANGE: NREQ_FACE=',NREQ_FACE
   WRITE(9,*) 'SCARC_EXCHANGE: NREQ_EDGE=',NREQ_EDGE
   WRITE(9,*) 'SCARC_EXCHANGE: NREQ_DIAG=',NREQ_DIAG
   WRITE(9,*) 'SCARC_EXCHANGE: BEDGE=',BEDGE
   WRITE(9,*) 'SCARC_EXCHANGE: BDIAG=',BDIAG
   WRITE(9,*) 'SCARC_EXCHANGE: IERR=',IERR
   WRITE(9,*) 'SCARC_EXCHANGE: REQ_FACE(1)=',REQ_FACE(1)
   WRITE(9,*) 'SCARC_EXCHANGE: REQ_FACE(2)=',REQ_FACE(2)
ENDIF

CALL MPI_WAITALL(NREQ_FACE,REQ_FACE(1:NREQ_FACE),STAT_FACE,IERR)
IF (BEDGE) CALL MPI_WAITALL(NREQ_EDGE,REQ_EDGE(1:NREQ_EDGE),STAT_EDGE,IERR)
IF (BDIAG) THEN
IF (SCARC_DEBUG>=6) WRITE(9,*) 'WAITING FOR DIAG COMM, NREQ_DIAG=',NREQ_DIAG
   CALL MPI_WAITALL(NREQ_DIAG,REQ_DIAG(1:NREQ_DIAG),STAT_DIAG,IERR)
ENDIF

IF (SCARC_DEBUG>=6) THEN
   WRITE(9,*) 'IERR=',IERR
   WRITE(9,*) 'SCARC_EXCHANGE:CODE=',CODE
   WRITE(9,*) 'AFTER WAIT '
ENDIF

 
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

IF (SCARC_DEBUG>=6) THEN
   WRITE(9,*) 'SNOM   zeigt auf SCARC(',NOM,')'
   WRITE(9,*) 'SNOML  zeigt auf SCARC(',NOM,')%SLEVEL(',ILEVEL,')'
ENDIF

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! FACE communication
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      RECV_FACE_IF: IF (SNOML%NIC_FACE(NOM,NM)/=0 .AND. SNOML%NIC_FACE(NM,NOM)/=0) THEN

         OSNOM => SCARC(NOM)%OSCARC(NM)                   ! corresponds to M2
         OSNOML=> SCARC(NOM)%OSCARC(NM)%SLEVEL(ILEVEL)  ! ... for the level 'ILEVEL'

IF (SCARC_DEBUG>=6) THEN
   WRITE(9,*) 'FACE: OSNOM  zeigt auf SCARC(',NOM,')%OSCARC(',NM,')'
   WRITE(9,*) 'FACE: OSNOML zeigt auf SCARC(',NOM,')%OSCARC(',NM,')%SLEVEL(',ILEVEL,')'
   WRITE(9,*) 'FACE: CODE=',CODE
ENDIF
    
         !!! Extract data for matrix-vector communication
         RECV_FACE_MATV_IF: IF (CODE==NCOM_MATV) THEN
     
            IF (ABS(SNOM%IOR_FACE(NM,NOM))==1) THEN
               ASUB=SNOML%ASUBX
            ELSE IF (ABS(SNOM%IOR_FACE(NM,NOM))==2) THEN
               ASUB=SNOML%ASUBY
            ELSE IF (ABS(SNOM%IOR_FACE(NM,NOM))==3) THEN
               ASUB=SNOML%ASUBZ
            ENDIF
IF (SCARC_DEBUG>=6) WRITE(9,*) 'ASUB=',ASUB

            IF (RNODE/=SNODE) THEN
               LL = 0

IF (SCARC_DEBUG>=6) WRITE(9,*) 'RECV_FACE:',OSNOM%RECV_FACE
               UNPACK_RECV_FACE0: DO
                  IW = NINT(OSNOM%RECV_FACE(LL+1))
IF (SCARC_DEBUG>=6) WRITE(9,*) '====================== IW=', IW
                  IF (IW==-999) EXIT UNPACK_RECV_FACE0
IF (SCARC_DEBUG>=6) WRITE(9,'(a,17i4)') 'IJKW_FACE=',SNOML%IJKW_FACE(1:17,IW)
                  ZSUM=0.0_EB
                  DO KK=SNOML%IJKW_FACE(12,IW),SNOML%IJKW_FACE(15,IW)
                     DO JJ=SNOML%IJKW_FACE(11,IW),SNOML%IJKW_FACE(14,IW)
                        DO II=SNOML%IJKW_FACE(10,IW),SNOML%IJKW_FACE(13,IW)
                           OSNOML%Z_FACE(II,JJ,KK) = OSNOM%RECV_FACE(LL+2)
                           ZSUM=ZSUM+OSNOM%RECV_FACE(LL+2)
IF (SCARC_DEBUG>=6) WRITE(9,'(i3,a,i3,a,e12.4,3i3,e12.4)') IW,': SCARC_EXCHANGE_RECV_FACE_FULL:RECV_FACE(',LL+2,')=',&
            OSNOM%RECV_FACE(LL+2), II, JJ, KK, ZSUM
                           LL = LL+2
!WRITE(9,*)  LL,II,JJ,KK
!IF (LL>14) stop
                        ENDDO
                     ENDDO
                  ENDDO
                  I=SNOML%IJKW_FACE(6,IW)
                  J=SNOML%IJKW_FACE(7,IW)
                  K=SNOML%IJKW_FACE(8,IW)
IF (SCARC_DEBUG>=6) WRITE(9,*) 'Uff1', I, J, K, IW
                  ISUM = (SNOML%IJKW_FACE(13,IW)-SNOML%IJKW_FACE(10,IW)+1) * &
                         (SNOML%IJKW_FACE(14,IW)-SNOML%IJKW_FACE(11,IW)+1) * &
                         (SNOML%IJKW_FACE(15,IW)-SNOML%IJKW_FACE(12,IW)+1)
IF (SCARC_DEBUG>=6) WRITE(9,*) 'ISUM=',ISUM, ': ZSUM_MEAN=',ZSUM/REAL(ISUM,EB),': ASUB*ZSUM_MEAN=',ASUB*ZSUM/REAL(ISUM,EB)

                  IF (IMV==1 .OR. IMV==5) THEN
                     yold=SNOML%Y(I,J,K)
                     SNOML%Y(I, J, K) = SNOML%Y(I, J, K) + ASUB * ZSUM/REAL(ISUM,EB)
IF (SCARC_DEBUG>=6) WRITE(9,'(i3,a,i2,a,i2,a,i2,a,2e12.4)') IW,&
      ':1SCARC_EXCHANGE_RECV_FACE_FULL: ====> Y(',I,',',J,',',K,')=',SNOML%Y(I,J,K),yold
                  ELSE IF (IMV==2 .OR. IMV==6) THEN
                     yold=SNOML%G(I,J,K)
                     SNOML%G(I, J, K) = SNOML%G(I, J, K) + ASUB * ZSUM/REAL(ISUM,EB)
IF (SCARC_DEBUG>=6) WRITE(9,'(i3,a,i2,a,i2,a,i2,a,2e12.4)') IW,&
      ':2SCARC_EXCHANGE_RECV_FACE_FULL: ====> G(',I,',',J,',',K,')=',SNOML%G(I,J,K),yold
                  ELSE IF (IMV==3) THEN
                     yold=SNOML%R(I,J,K)
                     SNOML%R(I, J, K) = SNOML%R(I, J, K) + ASUB * ZSUM/REAL(ISUM,EB)
IF (SCARC_DEBUG>=6) WRITE(9,'(i3,a,i2,a,i2,a,i2,a,2e12.4)') IW,&
      ':3SCARC_EXCHANGE_RECV_FACE_FULL: ====> R(',I,',',J,',',K,')=',SNOML%R(I,J,K),yold
                  ENDIF
               ENDDO UNPACK_RECV_FACE0
            ENDIF

         ENDIF RECV_FACE_MATV_IF


         !!! Extract data for full communication including edge (3D!) and vertex data
         RECV_FACE_FULL_IF: IF (CODE==NCOM_FULL) THEN
   
            IF (RNODE/=SNODE) THEN
               LL = 0
               UNPACK_RECV_FACE: DO
                  IW = NINT(OSNOM%RECV_FACE(LL+1))
IF (SCARC_DEBUG>=6) WRITE(9,*) '====================== IW=', IW
                  IF (IW==-999) EXIT UNPACK_RECV_FACE
                  ZSUM=0.0_EB
                  DO KK=SNOML%IJKW_FACE(12,IW),SNOML%IJKW_FACE(15,IW)
                     DO JJ=SNOML%IJKW_FACE(11,IW),SNOML%IJKW_FACE(14,IW)
                        DO II=SNOML%IJKW_FACE(10,IW),SNOML%IJKW_FACE(13,IW)
                           OSNOML%Z_FACE(II,JJ,KK) = OSNOM%RECV_FACE(LL+2)
                           ZSUM=ZSUM+OSNOM%RECV_FACE(LL+2)
                           LL = LL+2
                        ENDDO
                     ENDDO
                  ENDDO
                  I=SNOML%IJKW_FACE(1,IW)
                  J=SNOML%IJKW_FACE(2,IW)
                  K=SNOML%IJKW_FACE(3,IW)
                  ISUM = (SNOML%IJKW_FACE(13,IW)-SNOML%IJKW_FACE(10,IW)+1) * &
                         (SNOML%IJKW_FACE(14,IW)-SNOML%IJKW_FACE(11,IW)+1) * &
                         (SNOML%IJKW_FACE(15,IW)-SNOML%IJKW_FACE(12,IW)+1)
                  SNOML%Z(I, J, K) = ZSUM/REAL(ISUM,EB)
IF (SCARC_DEBUG>=6) WRITE(9,*) IW,': SCARC_EXCHANGE_RECV_FACE_FULL: ====> Z(',I,',',J,',',K,')=',SNOML%Z(I,J,K)
               ENDDO UNPACK_RECV_FACE

            ENDIF
         ENDIF RECV_FACE_FULL_IF

      ENDIF RECV_FACE_IF


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! EDGE communication
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      RECV_EDGE_IF: IF (BEDGE.AND.SNOML%NIC_EDGE(NOM,NM)/=0 .AND. SNOML%NIC_EDGE(NM,NOM)/=0) THEN

         OSNOM => SCARC(NOM)%OSCARC(NM)                   ! corresponds to M2
         OSNOML=> SCARC(NOM)%OSCARC(NM)%SLEVEL(ILEVEL)  ! ... for the level 'ILEVEL'

    
         !!! Extract data for full communication including edge (3D!) and vertex data
         RECV_EDGE_FULL_IF: IF (CODE==NCOM_FULL) THEN
   
            IF (RNODE/=SNODE) THEN
               LL = 0
write(*,*) 'ACHTUNG: EDGE: STIMMT NOCH NICHT !!!!!!!!!!!!!!!!!!'
                  IMIN=SNOML%I_MIN_EDGE(NM,NOM)-1
                  IMAX=SNOML%I_MAX_EDGE(NM,NOM)-1
                  JMIN=SNOML%J_MIN_EDGE(NM,NOM)
                  JMAX=SNOML%J_MAX_EDGE(NM,NOM)
                  KMIN=SNOML%K_MIN_EDGE(NM,NOM)
                  KMAX=SNOML%K_MAX_EDGE(NM,NOM)
IF (SCARC_DEBUG>=6) THEN
   WRITE(9,*) 'NCOM_FULL: SCARC_EXCHANGE: IMIN=',IMIN
   WRITE(9,*) 'NCOM_FULL: SCARC_EXCHANGE: IMAX=',IMAX
   WRITE(9,*) 'NCOM_FULL: SCARC_EXCHANGE: JMIN=',JMIN
   WRITE(9,*) 'NCOM_FULL: SCARC_EXCHANGE: JMAX=',JMAX
   WRITE(9,*) 'NCOM_FULL: SCARC_EXCHANGE: KMIN=',KMIN
   WRITE(9,*) 'NCOM_FULL: SCARC_EXCHANGE: KMAX=',KMAX
ENDIF
               LL = 0
               DO KK=KMIN,KMAX
                  DO JJ=JMIN,JMAX
                     DO II=IMIN,IMAX
                        SNOML%Z(II, JJ, KK) = OSNOM%RECV_EDGE(LL+1)
IF (SCARC_DEBUG>=6) WRITE(9,*) NOM,': RECEIVING from ',NM,': Z(',II,',',JJ,',',KK,')=',OSNOM%RECV_EDGE(LL+1)
                        LL = LL+1
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF
         ENDIF RECV_EDGE_FULL_IF

      ENDIF RECV_EDGE_IF

IF (SCARC_DEBUG>=6) THEN
   WRITE(9,*) 'BDIAG=',BDIAG
   WRITE(9,*) 'NIC_DIAG=',SNOML%NIC_DIAG(NOM,NM)
   WRITE(9,*) 'NIC_DIAG=',SNOML%NIC_DIAG(NM,NOM)
ENDIF
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! DIAG communication
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      RECV_DIAG_IF: IF (BDIAG.AND.SNOML%NIC_DIAG(NOM,NM)/=0 .AND. SNOML%NIC_DIAG(NM,NOM)/=0) THEN

         OSNOM => SCARC(NOM)%OSCARC(NM)                   ! corresponds to M2
         OSNOML=> SCARC(NOM)%OSCARC(NM)%SLEVEL(ILEVEL)  ! ... for the level 'ILEVEL'

IF (SCARC_DEBUG>=6) THEN
   WRITE(9,*) 'DIAG: OSNOM  zeigt auf SCARC(',NOM,')%OSCARC(',NM,')'
   WRITE(9,*) 'DIAG: OSNOML zeigt auf SCARC(',NOM,')%OSCARC(',NM,')%SLEVEL(',ILEVEL,')'
   WRITE(9,*) 'DIAG: CODE=',CODE
ENDIF
    
         !!! Extract data for full communication including edge (3D!) and vertex data
         RECV_DIAG_FULL_IF: IF (CODE==NCOM_FULL) THEN
   
            IF (RNODE/=SNODE) THEN
               LL = 0
IF (SCARC_DEBUG>=6) WRITE(9,'(a,20f8.2)') 'RECV_DIAG=',OSNOM%RECV_DIAG
               UNPACK_RECV_DIAG: DO
                  IW = NINT(OSNOM%RECV_DIAG(LL+1))
IF (SCARC_DEBUG>=6) WRITE(9,*) '====================== IW=', IW
                  IF (IW==-999) EXIT UNPACK_RECV_DIAG
                  ZSUM=0.0_EB
                  DO KK=SNOML%IJKW_DIAG(12,IW),SNOML%IJKW_DIAG(15,IW)
                     DO JJ=SNOML%IJKW_DIAG(11,IW),SNOML%IJKW_DIAG(14,IW)
                        DO II=SNOML%IJKW_DIAG(10,IW),SNOML%IJKW_DIAG(13,IW)
                           OSNOML%Z_DIAG(II,JJ,KK) = OSNOM%RECV_DIAG(LL+2)
                           ZSUM=ZSUM+OSNOM%RECV_DIAG(LL+2)
                           LL = LL+2
                        ENDDO
                     ENDDO
                  ENDDO
                  I=SNOML%IJKW_DIAG(1,IW)
                  J=SNOML%IJKW_DIAG(2,IW)
                  K=SNOML%IJKW_DIAG(3,IW)
                  ISUM = (SNOML%IJKW_DIAG(13,IW)-SNOML%IJKW_DIAG(10,IW)+1) * &
                         (SNOML%IJKW_DIAG(14,IW)-SNOML%IJKW_DIAG(11,IW)+1) * &
                         (SNOML%IJKW_DIAG(15,IW)-SNOML%IJKW_DIAG(12,IW)+1)
                  SNOML%Z(I, J, K) = ZSUM/REAL(ISUM,EB)
IF (SCARC_DEBUG>=6) WRITE(9,*) IW,': SCARC_EXCHANGE_RECV_DIAG: ====> Z(',I,',',J,',',K,')=',SNOML%Z(I,J,K)
               ENDDO UNPACK_RECV_DIAG

            ENDIF

         ENDIF RECV_DIAG_FULL_IF

      ENDIF RECV_DIAG_IF

   
   ENDDO EXCHANGE_RECV2_MESH_LOOP
   
ENDDO EXCHANGE_SEND2_MESH_LOOP
IF (SCARC_DEBUG>=6) WRITE(9,*) 'AFTER SEND_LOOP '

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
INTEGER :: I0, J0, K0, I1, J1, K1, I, K
REAL(EB):: TNOW_GHOSTCELLS
 
TYPE (SCARC_TYPE),  POINTER ::  S
TYPE (SCARC_LEVEL_TYPE),  POINTER ::  SL
TYPE (MESH_TYPE),  POINTER ::  M

TNOW_GHOSTCELLS = SECOND()

S => SCARC(NM)
SL => S%SLEVEL(S%NLEVEL_MAX)
M  => MESHES(NM)
 
IF (SCARC_DEBUG >= 2) THEN
   WRITE (9,*) '========================Starting SCARC_PRESSURE_GHOSTCELLS2D'
   WRITE (9,*)
   WRITE (9,*) 'NM   :',NM
   WRITE (9,*) 'NLEVEL_MAX:',S%NLEVEL_MAX
   WRITE (9,*) 'IBAR0:',M%IBAR
   WRITE (9,*) 'KBAR0:',M%KBAR
   WRITE (9,*) 'NW_FACE:',SL%NW_FACE
   WRITE (9,*) 'NEWC:',M%NEWC
   WRITE (9,*)
   WRITE (9,*) 'DXI :',DXI
   WRITE (9,*) 'DYI :',DETA
   WRITE (9,*) 'DZI :',DZETA
   WRITE (9,*)
   WRITE (9,*) 'IJKW:'
   WRITE (9,'(15i4)') M%IJKW
   WRITE (9,*)
   WRITE (9,*) 'BTYPE:'
   WRITE (9,'(8i4)') M%PRESSURE_BC_INDEX(1:M%NEWC)
   WRITE (9,*)
   WRITE (9,*) 'M%BXS:'
   WRITE (9,'(8f16.6)') (M%BXS(1, K), K=1, M%KBAR)
   WRITE (9,*)
   WRITE (9,*) 'M%BXF:'
   WRITE (9,'(8f16.6)') (M%BXF(M%IBAR, K), K=1, M%KBAR)
   WRITE (9,*)
   WRITE (9,*) 'M%BZS:'
   WRITE (9,'(8f16.6)') (M%BZS(I, 1), I=1, M%IBAR)
   WRITE (9,*)
   WRITE (9,*) 'M%BZF:'
   WRITE (9,'(8f16.6)') (M%BZF(I, M%JBAR), I=1, M%IBAR)
   WRITE (9,*)
   WRITE (9,*)

ENDIF

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
      ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
         HP(I0,J0,K0) =  HP(I1,J1,K1) - DZETA *M%BZS(I1,J1)
      ENDIF
   ELSE IF (IOR0 == -3) THEN
      IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
         HP(I0,J0,K0) = -HP(I1,J1,K1) + 2.0_EB * M%BZF(I1,J1)
      ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
         HP(I0,J0,K0) =  HP(I1,J1,K1) + DZETA *M%BZF(I1,J1)
      ENDIF
   ENDIF

ENDDO
 
TUSED_SCARC(25,NM)=TUSED_SCARC(25,NM)+SECOND()-TNOW_GHOSTCELLS
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_GHOSTCELLS

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
 
IF (SCARC_DEBUG >= 2) THEN
   WRITE (9,*) '========================Starting SCARC_SETBDRY2D'
   WRITE (9,*)
   WRITE (9,*) 'NM   :',NM
   WRITE (9,*) 'ILEVEL:',ILEVEL
   WRITE (9,*) 'IBAR0:',SL%IBAR
   WRITE (9,*) 'KBAR0:',SL%KBAR
   WRITE (9,*) 'NW_FACE:',SL%NW_FACE
   WRITE (9,*)
   WRITE (9,*) 'DXI :',SL%DXI
   WRITE (9,*) 'DZI :',SL%DZI
   WRITE (9,*)
   WRITE (9,*) 'DXI2:',SL%DXI2
   WRITE (9,*) 'DZI2:',SL%DZI2
   WRITE (9,*)
   WRITE (9,*) 'IJKW0:'
   WRITE (9,'(17i4)') SL%IJKW_FACE
   WRITE (9,*)
   WRITE (9,*) 'BTYPE:'
   WRITE (9,'(8i4)') M%PRESSURE_BC_INDEX(1:SL%NW_FACE)
   WRITE (9,*)
   WRITE (9,*) 'M%BXS:'
   WRITE (9,*) (M%BXS(1, K), K=1, M%KBAR)
   WRITE (9,*)
   WRITE (9,*) 'M%BXF:'
   WRITE (9,*) (M%BXF(1, K), K=1, M%KBAR)
   WRITE (9,*)
   WRITE (9,*) 'M%BZS:'
   WRITE (9,*) (M%BZS(I, 1), I=1, M%IBAR)
   WRITE (9,*)
   WRITE (9,*) 'M%BZF:'
   WRITE (9,*) (M%BZF(I, 1), I=1, M%IBAR)
   WRITE (9,*)
   WRITE (9,*) 'F(.,.):INIT'
   WRITE (9, '(8e16.8)') ((SL%F(I, 1, K), I=1, SL%IBAR), K=SL%KBAR, 1,-1)
   WRITE (9,*)
   WRITE (9,*)

ENDIF

DO IW = 1, M%NEWC

   I = M%IJKW(6,IW)
   J = M%IJKW(7,IW)
   K = M%IJKW(8,IW)

   IOR0 = M%IJKW(4,IW)

   IF (IOR0 == 1) THEN
      IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
         SL%F(I,J,K) = SL%F(I,J,K) - 2.0_EB * SL%DXI2 * M%BXS(1,K)            ! Dirichlet
         IF (SCARC_DEBUG>=2) WRITE(9,1000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
      ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
         SL%F(I,J,K) = SL%F(I,J,K) + SL%DXI * M%BXS(1,K)                       ! Neumann
         IF (SCARC_DEBUG>=2) WRITE(9,2000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
      ENDIF
   ELSE IF (IOR0 == -1) THEN
      IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
         SL%F(I,J,K) = SL%F(I,J,K) - 2.0_EB * SL%DXI2 *M%BXF(1,K)            ! Dirichlet
         IF (SCARC_DEBUG>=2) WRITE(9,1000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
      ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
         SL%F(I,J,K) = SL%F(I,J,K) - SL%DXI *M%BXF(1,K)                      ! Neumann
         IF (SCARC_DEBUG>=2) WRITE(9,2000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
      ENDIF
   ELSE IF (IOR0 == 3) THEN
      IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
         SL%F(I,J,K) = SL%F(I,J,K) - 2.0_EB * SL%DZI2 * M%BZS(I,1)            ! Dirichlet
         IF (SCARC_DEBUG>=2) WRITE(9,1000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
      ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
         SL%F(I,J,K) = SL%F(I,J,K) + SL%DZI * M%BZS(I,1)                       ! Neumann
         IF (SCARC_DEBUG>=2) WRITE(9,2000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
      ENDIF
   ELSE IF (IOR0 == -3) THEN
      IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
         SL%F(I,J,K) = SL%F(I,J,K) - 2.0_EB * SL%DZI2 * M%BZF(I,1)            ! Dirichlet
         IF (SCARC_DEBUG>=2) WRITE(9,1000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
      ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
         SL%F(I,J,K) = SL%F(I,J,K) - SL%DZI  * M%BZF(I,1)                     ! Neumann
         IF (SCARC_DEBUG>=2) WRITE(9,2000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
      ENDIF
   ENDIF

ENDDO

IF (SCARC_DEBUG >= 2) THEN
   WRITE (9,*)
   WRITE (9,*) 'AFTER CHANGE: F(.,.)'
   WRITE (9,*) 'IBAR=', IBAR
   WRITE (9,*) 'KBAR=', KBAR
   WRITE (9, '(8e16.8)') ((SL%F(I, 1, K), I=1, SL%IBAR), K=SL%KBAR, 1,-1)
ENDIF
 
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
 
IF (SCARC_DEBUG >= 2) THEN
   WRITE (9,*) '========================Starting SCARC_SETBDRY3D'
   WRITE (9,*)
   WRITE (9,*) 'NM   :',NM
   WRITE (9,*) 'ILEVEL:',ILEVEL
   WRITE (9,*) 'IBAR0:',SL%IBAR
   WRITE (9,*) 'KBAR0:',SL%KBAR
   WRITE (9,*) 'NW_FACE:',SL%NW_FACE
   WRITE (9,*)
   WRITE (9,*) 'DXI :',SL%DXI
   WRITE (9,*) 'DZI :',SL%DZI
   WRITE (9,*)
   WRITE (9,*) 'DXI2:',SL%DXI2
   WRITE (9,*) 'DZI2:',SL%DZI2
   WRITE (9,*)
   WRITE (9,*) 'IJKW0:'
   WRITE (9,'(17i4)') SL%IJKW_FACE
   WRITE (9,*)
   WRITE (9,*) 'BTYPE:'
   WRITE (9,'(8i4)') M%PRESSURE_BC_INDEX(1:SL%NW_FACE)
   WRITE (9,*)
   WRITE (9,*) 'M%BXS:'
   WRITE (9,*) (M%BXS(1, K), K=1, M%KBAR)
   WRITE (9,*)
   WRITE (9,*) 'M%BXF:'
   WRITE (9,*) (M%BXF(1, K), K=1, M%KBAR)
   WRITE (9,*)
   WRITE (9,*) 'M%BYS:'                               ! to be checked
   WRITE (9,*) (M%BYS(J, 1), J=1, M%JBAR)
   WRITE (9,*)
   WRITE (9,*) 'M%BYF:'
   WRITE (9,*) (M%BYF(J, 1), J=1, M%JBAR)
   WRITE (9,*)
   WRITE (9,*) 'M%BZS:'
   WRITE (9,*) (M%BZS(I, 1), I=1, M%IBAR)
   WRITE (9,*)
   WRITE (9,*) 'M%BZF:'
   WRITE (9,*) (M%BZF(I, 1), I=1, M%IBAR)
   WRITE (9,*)
   WRITE (9,*) 'F(.,.):INIT'
   WRITE (9, '(8e16.8)') ((SL%F(I, 1, K), I=1, SL%IBAR), K=SL%KBAR, 1,-1)
   WRITE (9,*)
   WRITE (9,*)

ENDIF

DO IW = 1, M%NEWC

   I = M%IJKW(6,IW)
   J = M%IJKW(7,IW)
   K = M%IJKW(8,IW)

   IOR0 = M%IJKW(4,IW)

   IF (IOR0 == 1) THEN
      IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
         SL%F(I,J,K) = SL%F(I,J,K) - 2.0_EB * SL%DXI2 * M%BXS(1,K)            ! Dirichlet
         IF (SCARC_DEBUG>=2) WRITE(9,1000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
      ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
         SL%F(I,J,K) = SL%F(I,J,K) + SL%DXI * M%BXS(1,K)                       ! Neumann
         IF (SCARC_DEBUG>=2) WRITE(9,2000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
      ENDIF
   ELSE IF (IOR0 == -1) THEN
      IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
         SL%F(I,J,K) = SL%F(I,J,K) - 2.0_EB * SL%DXI2 *M%BXF(1,K)            ! Dirichlet
         IF (SCARC_DEBUG>=2) WRITE(9,1000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
      ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
         SL%F(I,J,K) = SL%F(I,J,K) - SL%DXI *M%BXF(1,K)                      ! Neumann
         IF (SCARC_DEBUG>=2) WRITE(9,2000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
      ENDIF
   ELSE IF (IOR0 == 3) THEN
      IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
         SL%F(I,J,K) = SL%F(I,J,K) - 2.0_EB * SL%DZI2 * M%BZS(I,1)            ! Dirichlet
         IF (SCARC_DEBUG>=2) WRITE(9,1000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
      ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
         SL%F(I,J,K) = SL%F(I,J,K) + SL%DZI * M%BZS(I,1)                       ! Neumann
         IF (SCARC_DEBUG>=2) WRITE(9,2000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
      ENDIF
   ELSE IF (IOR0 == -3) THEN
      IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
         SL%F(I,J,K) = SL%F(I,J,K) - 2.0_EB * SL%DZI2 * M%BZF(I,1)            ! Dirichlet
         IF (SCARC_DEBUG>=2) WRITE(9,1000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
      ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
         SL%F(I,J,K) = SL%F(I,J,K) - SL%DZI  * M%BZF(I,1)                     ! Neumann
         IF (SCARC_DEBUG>=2) WRITE(9,2000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
      ENDIF
   ENDIF

ENDDO
 
IF (SCARC_DEBUG >= 2) THEN
   WRITE (9,*)
   WRITE (9,*) 'AFTER CHANGE: F(.,.)'
   WRITE (9,*) 'IBAR=', IBAR
   WRITE (9,*) 'KBAR=', KBAR
   WRITE (9, '(8e16.8)') ((SL%F(I, 1, K), I=1, SL%IBAR), K=SL%KBAR, 1,-1)
ENDIF
 
TUSED_SCARC(27,NM)=TUSED_SCARC(27,NM)+SECOND()-TNOW_SETBDRY3D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_SETBDRY3D
1000 FORMAT(i3,': IOR=',i3,': Dirichlet :F(',i2,',',i2,',',i2,')=',e25.16,i5,i5)
2000 FORMAT(i3,': IOR=',i3,': Neumann   :F(',i2,',',i2,',',i2,')=',e25.16,i5,i5)
END SUBROUTINE SCARC_SETBDRY3D
 
 

 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Update different quantities according to given type
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_UPDATE (MYID, ITYPE)
INTEGER :: MYID, ILEVEL, ITYPE, NM
REAL(EB):: TNOW_UPDATE

TNOW_UPDATE=SECOND()

NM = MYID + 1
S => SCARC(NM)
 
UPDATE_LEVEL_LOOP: DO ILEVEL=S%NLEVEL_MAX,S%NLEVEL_MIN,-1

   IF (ITYPE == NUPDATE_FULL) THEN
      CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%U, NM, ILEVEL, 'U     ')
      IF ( .NOT. TWO_D) CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%V, NM, ILEVEL, 'V     ')
      CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%W, NM, ILEVEL, 'W     ')
      CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%H, NM, ILEVEL, 'H     ')
      CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%D, NM, ILEVEL, 'D     ')
      CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%RHO, NM, ILEVEL, 'RHO   ')
      CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%MU, NM, ILEVEL, 'MU    ')
      CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%US, NM, ILEVEL, 'US    ')
      IF ( .NOT. TWO_D) CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%VS, NM, ILEVEL, 'VS    ')
      CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%WS, NM, ILEVEL, 'WS    ')
      CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%DS, NM, ILEVEL, 'DS    ')
      CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%RHOS, NM, ILEVEL, 'RHOS  ')
      CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%DS, NM, ILEVEL, 'DS    ')
      CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%DDDT, NM, ILEVEL, 'DDDT  ')
      CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%Q, NM, ILEVEL, 'Q     ')
      CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%QR, NM, ILEVEL, 'QR    ')
      CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%TMP, NM, ILEVEL, 'TMP   ')
   ELSE IF (ITYPE == NUPDATE_MEDIUM) THEN
      CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%U, NM, ILEVEL, 'U     ')
      IF ( .NOT. TWO_D) CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%V, NM, ILEVEL, 'V     ')
      CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%W, NM, ILEVEL, 'W     ')
      CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%H, NM, ILEVEL, 'H     ')
      CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%D, NM, ILEVEL, 'D     ')
      CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%US, NM, ILEVEL, 'US    ')
      IF ( .NOT. TWO_D) CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%VS, NM, ILEVEL, 'VS    ')
      CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%WS, NM, ILEVEL, 'WS    ')
      CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%DS, NM, ILEVEL, 'DS    ')
      CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%DS, NM, ILEVEL, 'DS   ')
      CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%DDDT, NM, ILEVEL, 'DDDT   ')
      CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%TMP, NM, ILEVEL, 'TMP    ')
   ELSE IF (ITYPE == NUPDATE_VEL) THEN
      CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%U, NM, ILEVEL, 'U     ')
      IF ( .NOT. TWO_D) CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%V, NM, ILEVEL, 'V     ')
      CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%W, NM, ILEVEL, 'W     ')
   ELSE IF (ITYPE == NUPDATE_FV) THEN
      CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%FVX, NM, ILEVEL, 'FVX   ')
      IF ( .NOT. TWO_D) CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%FVY, NM, ILEVEL, 'FVY   ')
      CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%FVZ, NM, ILEVEL, 'FVZ   ')
   ELSE IF (ITYPE == NUPDATE_Q) THEN
      CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%Q, NM, ILEVEL, 'Q     ')
      CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%QR, NM, ILEVEL, 'QR    ')
   ELSE IF (ITYPE == NUPDATE_H) THEN
      CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%H, NM, ILEVEL, 'H     ')
   ELSE IF (ITYPE == NUPDATE_HS) THEN
      CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%HS, NM, ILEVEL, 'HS    ')
   ENDIF

ENDDO UPDATE_LEVEL_LOOP
 
TUSED_SCARC(28,NM)=TUSED_SCARC(28,NM)+SECOND()-TNOW_UPDATE
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_UPDATE
 
END SUBROUTINE SCARC_UPDATE



SUBROUTINE SCARC_UPDATE_LEVEL (MYID, ILEVEL, ITYPE)
REAL(EB):: TNOW_UPDATE_LEVEL
INTEGER :: MYID, ILEVEL, ITYPE, NM

TNOW_UPDATE_LEVEL=SECOND()
 
NM = MYID + 1
 
IF (ITYPE == NUPDATE_FULL) THEN
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%U, NM, ILEVEL, 'U     ')
   IF ( .NOT. TWO_D) CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%V, NM, ILEVEL, 'V     ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%W, NM, ILEVEL, 'W     ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%H, NM, ILEVEL, 'H     ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%D, NM, ILEVEL, 'D     ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%RHO, NM, ILEVEL, 'RHO   ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%MU, NM, ILEVEL, 'MU    ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%US, NM, ILEVEL, 'US    ')
   IF ( .NOT. TWO_D) CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%VS, NM, ILEVEL, 'VS    ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%WS, NM, ILEVEL, 'WS    ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%DS, NM, ILEVEL, 'DS    ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%RHOS, NM, ILEVEL, 'RHOS  ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%DS, NM, ILEVEL, 'DS    ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%DDDT, NM, ILEVEL, 'DDDT  ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%Q, NM, ILEVEL, 'Q     ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%QR, NM, ILEVEL, 'QR    ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%TMP, NM, ILEVEL, 'TMP   ')
ELSE IF (ITYPE == NUPDATE_MEDIUM) THEN
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%U, NM, ILEVEL, 'U     ')
   IF ( .NOT. TWO_D) CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%V, NM, ILEVEL, 'V     ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%W, NM, ILEVEL, 'W     ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%H, NM, ILEVEL, 'H     ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%D, NM, ILEVEL, 'D     ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%US, NM, ILEVEL, 'US    ')
   IF ( .NOT. TWO_D) CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%VS, NM, ILEVEL, 'VS    ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%WS, NM, ILEVEL, 'WS    ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%DS, NM, ILEVEL, 'DS    ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%DS, NM, ILEVEL, 'DS   ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%DDDT, NM, ILEVEL, 'DDDT   ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%TMP, NM, ILEVEL, 'TMP    ')
ELSE IF (ITYPE == NUPDATE_VEL) THEN
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%U, NM, ILEVEL, 'U     ')
   IF ( .NOT. TWO_D) CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%V, NM, ILEVEL, 'V     ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%W, NM, ILEVEL, 'W     ')
ELSE IF (ITYPE == NUPDATE_FV) THEN
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%FVX, NM, ILEVEL, 'FVX   ')
   IF ( .NOT. TWO_D) CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%FVY, NM, ILEVEL, 'FVY   ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%FVZ, NM, ILEVEL, 'FVZ   ')
ELSE IF (ITYPE == NUPDATE_Q) THEN
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%Q, NM, ILEVEL, 'Q     ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%QR, NM, ILEVEL, 'QR    ')
ELSE IF (ITYPE == NUPDATE_H) THEN
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%H, NM, ILEVEL, 'H     ')
ENDIF
 
TUSED_SCARC(29,NM)=TUSED_SCARC(29,NM)+SECOND()-TNOW_UPDATE_LEVEL
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_UPDATE_LEVEL
 
END SUBROUTINE SCARC_UPDATE_LEVEL
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!  update interior boundary values along vertical directions
!!!  to have consistent values all over the domain
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_UPDATE_QUANTITY (Z, NM, ILEVEL, CNAME)
REAL (EB), POINTER, DIMENSION (:, :, :) :: Z
INTEGER :: NM, ILEVEL, I, J, K, III, JJJ, NVRTX
CHARACTER (4) :: CNAME
REAL(EB):: TNOW_UPDATE_QUANTITY
 
INCLUDE 'mpif.h'
 

TNOW_UPDATE_QUANTITY = SECOND()

SL => S%SLEVEL(ILEVEL)
 
IF (NMESHES > 1) THEN
 
   if (SCARC_DEBUG>=2) WRITE(9,*) 'UPDATING ',CNAME
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! 2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   IF (TWO_D) THEN
 
      DO K = 0, KBAR + 1
         DO I = 0, IBAR + 1
            SL%Z (I, 1, K) = Z (I, 1, K)
         ENDDO
      ENDDO
      IF (SCARC_DEBUG >= 6) CALL SCARC_COMPARE_SINGLE0 (Z, 'SARC', 1, 'ZIN   ',0)
      NREQ_FACE = 0
      NREQ_DIAG = 0
      CALL SCARC_RECEIVE (NCOM_FULL, ILEVEL)   !Aufruf 2,0
      CALL SCARC_EXCHANGE (NCOM_FULL, ILEVEL, NTYPE_UPDATE)
      DO K = 0, KBAR + 1
         DO I = 0, IBAR + 1
            Z (I, 1, K) = SL%Z(I, 1, K)
         ENDDO
      ENDDO
      IF (NM==1) THEN
         Z(SL%IBAR+1,1,0)=Z(SL%IBAR+1,1,1)
         Z(SL%IBAR+1,1,SL%KBAR+1)=Z(SL%IBAR+1,1,SL%KBAR)
      ELSE IF (NM==2) THEN
         Z(0,1,0)=Z(0,1,1)
         Z(0,1,SL%KBAR+1)=Z(0,1,SL%KBAR)
      ENDIF
      IF (SCARC_DEBUG >= 6) CALL SCARC_COMPARE_SINGLE0 (Z, 'SARC', 1, 'ZOUT  ',0)


 
 
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
 
      IF (SCARC_DEBUG >= 6) THEN
         WRITE (9,*) 'SCARC_UPDATE_QUANTITY: ==================================='
         WRITE (9,*) 'SCARC_UPDATE_QUANTITY: NMESHES=', NMESHES
         WRITE (9,*) 'SCARC_UPDATE_QUANTITY: NM=', NM
         WRITE (9,*) 'SCARC_UPDATE_QUANTITY: CNAME=', CNAME
         WRITE (9,*) 'SCARC_UPDATE_QUANTITY: TWO_D=', TWO_D
         WRITE (9,*) 'NBR_DIAG:'
         WRITE (9, '( 8i4)') ((NBR_DIAG(JJJ, III), III=1, NVRTX), JJJ=1, NMESHES)
         WRITE (9,*) 'NBR_DIAG(', NM, '):'
         WRITE (9, '( 8i4)') (NBR_DIAG(NM, III), III=1, NVRTX)
         WRITE (9,*) 'NBR_EDGE:'
         WRITE (9, '( 12i4)') ((NBR_EDGE(JJJ, III), III=1, 12), JJJ=1, NMESHES)
         WRITE (9,*) 'NBR_EDGE(', NM, '):'
         WRITE (9, '( 12i4)') (NBR_EDGE(NM, III), III=1, 12)
         WRITE (9,*)
         WRITE (9,*) '-------SCARC_COPYBDRY : BEFORE EXCHANGE Z'
         DO K = KBAR + 1, 0, - 1
            WRITE (9, '(6f10.5)') ((Z(I, J, K), I=0, IBAR+1), J=JBAR+1, 0,-1)
            WRITE (9,*) '---------------------------------------'
         ENDDO
         WRITE (9,*) '-------SCARC_COPYBDRY : BEFORE EXCHANGE SL%Z'
         DO K = KBAR + 1, 0, - 1
            WRITE (9, '(6f10.5)') ((SL%Z(I, J, K), I=0, IBAR+1), J=JBAR+1, 0,-1)
            WRITE (9,*) '---------------------------------------'
         ENDDO
         CALL flush (9)
      ENDIF
 
 
      NREQ_FACE = 0
      NREQ_EDGE = 0
      NREQ_DIAG = 0
 
      CALL SCARC_RECEIVE (NCOM_FULL, ILEVEL)  !Aufruf 2,0
      CALL SCARC_EXCHANGE (NCOM_FULL, ILEVEL, NTYPE_UPDATE)
 
 
      DO K = 0, KBAR + 1
         DO J = 0, JBAR + 1
            DO I = 0, IBAR + 1
               Z (I, J, K) = SL%Z(I, J, K)
            ENDDO
         ENDDO
      ENDDO
 
      IF (SCARC_DEBUG >= 6) THEN
         WRITE (9,*) '-------SCARC_COPYBDRY : AFTER EXCHANGE SL%Y'
         IF (NMESHES == 1) THEN
            DO K = KBAR + 1, 0, - 1
               WRITE (9, '(10f10.5)') ((SL%Z(I, J, K), I=0, IBAR+1), J=JBAR+1, 0,-1)
               WRITE (9,*) '---------------------------------------'
            ENDDO
         ELSE
            DO K = KBAR + 1, 0, - 1
               WRITE (9, '(6f10.5)') ((SL%Z(I, J, K), I=0, IBAR+1), J=JBAR+1, 0,-1)
               WRITE (9,*) '---------------------------------------'
            ENDDO
         ENDIF
         CALL flush (9)
      ENDIF
 
   ENDIF
 
ENDIF
 
TUSED_SCARC(30,NM)=TUSED_SCARC(30,NM)+SECOND()-TNOW_UPDATE_QUANTITY
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_UPDATE_QUANTITY
END SUBROUTINE SCARC_UPDATE_QUANTITY
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!  update interior boundary values along vertical directions
!!!  to have consistent values all over the domain
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_UPDATE_QUANTITY2 (Z, NM, ILEVEL, CNAME)
REAL (EB), POINTER, DIMENSION (:, :, :) :: Z
INTEGER :: NM, ILEVEL, I, J, K, III, JJJ, NVRTX
CHARACTER (4) :: CNAME
REAL(EB):: TNOW_UPDATE_QUANTITY
 
INCLUDE 'mpif.h'
 
TNOW_UPDATE_QUANTITY = SECOND()
SL => S%SLEVEL(ILEVEL)
 
IF (NMESHES > 1) THEN
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! 2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   IF (TWO_D) THEN
 
      DO K = 0, KBAR + 1
         DO I = 0, IBAR + 1
            SL%Z (I, 1, K) = Z (I, 1, K)
         ENDDO
      ENDDO
      IF (SCARC_DEBUG >= 6) CALL SCARC_COMPARE_SINGLE0 (SL%Z, 'SARC', 1, 'ZIN   ',0)
      NREQ_FACE = 0
      NREQ_DIAG = 0
      CALL SCARC_RECEIVE (NCOM_FULL, ILEVEL)   !Aufruf 2,0
      CALL SCARC_EXCHANGE (NCOM_FULL, ILEVEL, NTYPE_UPDATE)
      DO K = 0, KBAR + 1
         DO I = 0, IBAR + 1
            Z (I, 1, K) = SL%Z(I, 1, K)
         ENDDO
      ENDDO
      IF (SCARC_DEBUG >= 6) CALL SCARC_COMPARE_SINGLE0 (SL%Z, 'SARC', 1, 'ZOUT  ',0)


 
      DO K = 0, KBAR + 1
         DO I = 0, IBAR + 1
            SL%Z (I, 1, K) = Z (I, 0, K)
         ENDDO
      ENDDO
      IF (SCARC_DEBUG >= 6) CALL SCARC_COMPARE_SINGLE0 (SL%Z, 'SARC', 1, 'ZIN   ',0)
      NREQ_FACE = 0
      NREQ_DIAG = 0
      CALL SCARC_RECEIVE (NCOM_FULL, ILEVEL)   !Aufruf 2,0
      CALL SCARC_EXCHANGE (NCOM_FULL, ILEVEL, NTYPE_UPDATE)
      DO K = 0, KBAR + 1
         DO I = 0, IBAR + 1
            Z (I, 0, K) = SL%Z(I, 1, K)
         ENDDO
      ENDDO
      IF (SCARC_DEBUG >= 6) CALL SCARC_COMPARE_SINGLE0 (SL%Z, 'SARC', 1, 'ZOUT  ',0)
 

      DO K = 0, KBAR + 1
         DO I = 0, IBAR + 1
            SL%Z (I, 1, K) = Z (I, 2, K)
         ENDDO
      ENDDO
      IF (SCARC_DEBUG >= 6) CALL SCARC_COMPARE_SINGLE0 (SL%Z, 'SARC', 1, 'ZIN   ',0)
      NREQ_FACE = 0
      NREQ_DIAG = 0
      CALL SCARC_RECEIVE (NCOM_FULL, ILEVEL)   !Aufruf 2,0
      CALL SCARC_EXCHANGE (NCOM_FULL, ILEVEL, NTYPE_UPDATE)
      DO K = 0, KBAR + 1
         DO I = 0, IBAR + 1
            Z (I, 2, K) = SL%Z(I, 1, K)
         ENDDO
      ENDDO
      IF (SCARC_DEBUG >= 6) CALL SCARC_COMPARE_SINGLE0 (SL%Z, 'SARC', 1, 'ZOUT  ',0)
 
 
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
 
      IF (SCARC_DEBUG >= 6) THEN
         WRITE (9,*) 'SCARC_UPDATE_QUANTITY: ==================================='
         WRITE (9,*) 'SCARC_UPDATE_QUANTITY: NMESHES=', NMESHES
         WRITE (9,*) 'SCARC_UPDATE_QUANTITY: NM=', NM
         WRITE (9,*) 'SCARC_UPDATE_QUANTITY: CNAME=', CNAME
         WRITE (9,*) 'SCARC_UPDATE_QUANTITY: TWO_D=', TWO_D
         WRITE (9,*) 'NBR_DIAG:'
         WRITE (9, '( 8i4)') ((NBR_DIAG(JJJ, III), III=1, NVRTX), JJJ=1, NMESHES)
         WRITE (9,*) 'NBR_DIAG(', NM, '):'
         WRITE (9, '( 8i4)') (NBR_DIAG(NM, III), III=1, NVRTX)
         WRITE (9,*) 'NBR_EDGE:'
         WRITE (9, '( 12i4)') ((NBR_EDGE(JJJ, III), III=1, 12), JJJ=1, NMESHES)
         WRITE (9,*) 'NBR_EDGE(', NM, '):'
         WRITE (9, '( 12i4)') (NBR_EDGE(NM, III), III=1, 12)
         WRITE (9,*)
         WRITE (9,*) '-------SCARC_COPYBDRY : BEFORE EXCHANGE Z'
         DO K = KBAR + 1, 0, - 1
            WRITE (9, '(6f10.5)') ((Z(I, J, K), I=0, IBAR+1), J=JBAR+1, 0,-1)
            WRITE (9,*) '---------------------------------------'
         ENDDO
         WRITE (9,*) '-------SCARC_COPYBDRY : BEFORE EXCHANGE SL%Z'
         DO K = KBAR + 1, 0, - 1
            WRITE (9, '(6f10.5)') ((SL%Z(I, J, K), I=0, IBAR+1), J=JBAR+1, 0,-1)
            WRITE (9,*) '---------------------------------------'
         ENDDO
         CALL flush (9)
      ENDIF
 
 
      NREQ_FACE = 0
      NREQ_EDGE = 0
      NREQ_DIAG = 0
 
      CALL SCARC_RECEIVE (NCOM_FULL, ILEVEL)  !Aufruf 2,0
      CALL SCARC_EXCHANGE (NCOM_FULL, ILEVEL, NTYPE_UPDATE)
 
 
      DO K = 0, KBAR + 1
         DO J = 0, JBAR + 1
            DO I = 0, IBAR + 1
               Z (I, J, K) = SL%Z(I, J, K)
            ENDDO
         ENDDO
      ENDDO
 
      IF (SCARC_DEBUG >= 6) THEN
         WRITE (9,*) '-------SCARC_COPYBDRY : AFTER EXCHANGE SL%Y'
         IF (NMESHES == 1) THEN
            DO K = KBAR + 1, 0, - 1
               WRITE (9, '(10f10.5)') ((SL%Z(I, J, K), I=0, IBAR+1), J=JBAR+1, 0,-1)
               WRITE (9,*) '---------------------------------------'
            ENDDO
         ELSE
            DO K = KBAR + 1, 0, - 1
               WRITE (9, '(6f10.5)') ((SL%Z(I, J, K), I=0, IBAR+1), J=JBAR+1, 0,-1)
               WRITE (9,*) '---------------------------------------'
            ENDDO
         ENDIF
         CALL flush (9)
      ENDIF
 
   ENDIF
 
ENDIF
 
TUSED_SCARC(30,NM)=TUSED_SCARC(30,NM)+SECOND()-TNOW_UPDATE_QUANTITY
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_UPDATE_QUANTITY
END SUBROUTINE SCARC_UPDATE_QUANTITY2
 
 
 
 
 
LOGICAL FUNCTION MATCH (VAL1, VAL2)
REAL (EB) :: VAL1, VAL2, TOL
TOL = 1.0E-8_EB
MATCH = .FALSE.
IF (Abs(VAL1-VAL2) <= TOL) MATCH = .TRUE.
END FUNCTION MATCH
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! print out vectors for comparison with own separate check program
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_COMPARE_ALL (MYID, CROUTINE, INUM, ILINE)
INTEGER :: MYID, INUM, ILINE
CHARACTER (4) :: CROUTINE
 
call flush(6)
call flush(9)

CALL SCARC_COMPARE_SINGLE (MESHES(MYID+1)%H, CROUTINE, INUM, 'H     ',ILINE)
CALL SCARC_COMPARE_SINGLE (MESHES(MYID+1)%MU, CROUTINE, INUM, 'MU    ',ILINE)
CALL SCARC_COMPARE_SINGLE (MESHES(MYID+1)%DDDT, CROUTINE, INUM, 'DDDT  ',ILINE)
CALL SCARC_COMPARE_SINGLE (MESHES(MYID+1)%D, CROUTINE, INUM, 'D     ',ILINE)
CALL SCARC_COMPARE_SINGLE (MESHES(MYID+1)%DS, CROUTINE, INUM, 'DS    ',ILINE)

CALL SCARC_COMPARE_SINGLE (MESHES(MYID+1)%U, CROUTINE, INUM, 'U      ',ILINE)
IF ( .NOT. TWO_D) CALL SCARC_COMPARE_SINGLE (MESHES(MYID+1)%V, CROUTINE, INUM, 'V      ',ILINE)
CALL SCARC_COMPARE_SINGLE (MESHES(MYID+1)%W, CROUTINE, INUM, 'W      ',ILINE)
 
CALL SCARC_COMPARE_SINGLE (MESHES(MYID+1)%US, CROUTINE, INUM, 'US     ',ILINE)
IF ( .NOT. TWO_D) CALL SCARC_COMPARE_SINGLE (MESHES(MYID+1)%VS, CROUTINE, INUM, 'VS     ',ILINE)
CALL SCARC_COMPARE_SINGLE (MESHES(MYID+1)%WS, CROUTINE, INUM, 'WS     ',ILINE)
 
CALL SCARC_COMPARE_SINGLE (MESHES(MYID+1)%FVX, CROUTINE, INUM, 'FVX    ',ILINE)
IF ( .NOT. TWO_D) CALL SCARC_COMPARE_SINGLE (MESHES(MYID+1)%FVY, CROUTINE, INUM, 'FVY    ',ILINE)
CALL SCARC_COMPARE_SINGLE (MESHES(MYID+1)%FVZ, CROUTINE, INUM, 'FVZ    ',ILINE)
 
END SUBROUTINE SCARC_COMPARE_ALL
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! only for debugging reasons ... print also given vector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_COMPARE_SINGLE (X, CROUTINE, INUM, CNAME, ILINE)
REAL (EB), POINTER, DIMENSION (:, :, :) :: X
INTEGER :: INUM, ILINE
INTEGER :: II, jj, KK
CHARACTER (4) :: CROUTINE
CHARACTER (6) :: CNAME
 
IF (SCARC_COMPARE >= 1 .AND. CROUTINE(1:4) /= "SARC") THEN

   SCARC_MSG = "msg/  /    _      _      "
   WRITE (SCARC_MSG(5:6), '(i2.2)') MYID + 1
   WRITE (SCARC_MSG(8:11), '(i4.4)') SCARC_COUNT2
   SCARC_MSG (13:16) = CROUTINE (1:4)
   WRITE (SCARC_MSG(17:18), '(i2.2)') INUM
   SCARC_MSG (20:25) = CNAME (1:6)
   SCARC_COUNT2 = SCARC_COUNT2 + 1
   OPEN (1000, FILE=SCARC_MSG)
 
   IF (TWO_D) THEN
      WRITE (1000,*) '=========================================================='
      WRITE (1000,'(a20,a4)') 'ROUTINE          : ',CROUTINE(1:4)
      WRITE (1000,'(a20,a6)') 'VECTOR           : ',CNAME(1:6)
      WRITE (1000,'(a20,i6)') 'CODE-LINE        : ',ILINE
      WRITE (1000,'(a20,i6)') 'INTERNAL POSITION: ',INUM
      WRITE (1000,'(a20,i6)') 'INTERNAL NUMBER  : ',SCARC_COUNT2-1
      WRITE (1000,*) '=========================================================='
      IF (NMESHES == 1) THEN
         DO KK = KBAR, 1, - 1
            WRITE (1000, '(8f22.16)') (X(II, 1, KK), II=1, IBAR)
         ENDDO
      ELSE
         DO KK = KBAR, 1, - 1
            WRITE (1000, '(4f22.16)') (X(II, 1, KK), II=1, IBAR)
         ENDDO
      ENDIF
   ELSE
      WRITE (1000,*) '=========================================================='
      WRITE (1000,*) '===        SCARC: ', CNAME
      WRITE (1000,*) '=========================================================='
      IF (NMESHES == 1) THEN
         DO KK = KBAR, 1, - 1
            WRITE (1000, '(8f22.16)') ((X(II, jj, KK), II=1, IBAR), jj=JBAR,1,-1)
            WRITE (1000,*) '----------------------------------------'
         ENDDO
      ELSE
         DO KK = KBAR, 1, - 1
            WRITE (1000, '(4f22.16)') ((X(II, jj, KK), II=1, IBAR), jj=JBAR,1,-1)
            WRITE (1000,*) '----------------------------------------'
         ENDDO
      ENDIF
   ENDIF

   CALL flush (1000)
   CLOSE (1000)
ENDIF
 
IF (SCARC_DEBUG >= 2) THEN
   IF (TWO_D) THEN
      WRITE (9,*) '============================================================='
      WRITE (9,*) '===   COMPARE= ',CNAME,': ROUTINE= ',CROUTINE
      WRITE (9,*) '============================================================='
      IF (NMESHES == 1) THEN
         DO KK = KBAR, 1, - 1
            WRITE (9, '(a,i3,a,8e12.3)') 'k=', KK, ' : ', (X(II, 1, KK), II=1,IBAR)
         ENDDO
      ELSE
         IF (IBAR==8) THEN
         DO KK = KBAR+1, 0, - 1
            WRITE (9, '(a,i3,a,8e12.3)') 'k=', KK, ' : ', (X(II, 1, KK), II=1, IBAR)
         ENDDO
         ELSE
         DO KK = KBAR, 1, - 1
            WRITE (9, '(a,i3,a,4e12.3)') 'k=', KK, ' : ', (X(II, 1, KK), II=1, IBAR)
         ENDDO
         ENDIF
      ENDIF
   ELSE
      WRITE (9,*) '============================================================='
      WRITE (9,*) '===   COMPARE= ',CNAME,': ROUTINE= ',CROUTINE
      WRITE (9,*) '============================================================='
      IF (NMESHES == 1) THEN
         DO KK = KBAR, 1, - 1
            WRITE (9, '(8e12.3)') ((X(II, jj, KK), II=1, IBAR), jj=JBAR, 1,-1)
            WRITE (9,*) '----------------------------------------'
         ENDDO
      ELSE
         DO KK = KBAR, 1, - 1
            WRITE (9, '(4e12.3)') ((X(II, jj, KK), II=1, IBAR), jj=JBAR, 1,-1)
            WRITE (9,*) '----------------------------------------'
         ENDDO
      ENDIF
   ENDIF
   CALL flush (9)
ENDIF
 
END SUBROUTINE SCARC_COMPARE_SINGLE
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! only for debugging reasons ... print also ghost values
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_COMPARE_SINGLE0 (X, CROUTINE, INUM, CNAME, ILINE)
REAL (EB), POINTER, DIMENSION (:, :, :) :: X
INTEGER :: INUM, ILINE
INTEGER :: II, jj, KK
CHARACTER (4) :: CROUTINE
CHARACTER (6) :: CNAME
 
IF (SCARC_DEBUG >= 2) THEN
   IF (TWO_D) THEN
      WRITE (9,*) '============================================================='
      WRITE (9,*) '===   COMPARE= ',CNAME,': SCARC_COUNT2= ',SCARC_COUNT2-1, ':ILINE=',ILINE,': INUM=',INUM
      WRITE (9,*) '============================================================='
      IF (NMESHES == 1) THEN
         DO KK = KBAR+1, 0, - 1
            WRITE (9, '(a,i3,a,10e12.3)') 'k=', KK, ' : ', (X(II, 1, KK), II=0,IBAR+1)
         ENDDO
      ELSE
         IF (IBAR==8) THEN
         DO KK = KBAR+1, 0, - 1
            WRITE (9, '(a,i3,a,10e12.3)') 'k=', KK, ' : ', (X(II, 1, KK), II=0, IBAR+1)
         ENDDO
         ELSE
         DO KK = KBAR+1, 0, - 1
            WRITE (9, '(a,i3,a,6e12.3)') 'k=', KK, ' : ', (X(II, 1, KK), II=0, IBAR+1)
         ENDDO
         ENDIF
      ENDIF
   ELSE
      WRITE (9,*) '============================================================='
      WRITE (9,*) '===   COMPARE= ',CNAME,': ROUTINE= ',CROUTINE
      WRITE (9,*) '============================================================='
      IF (NMESHES == 1) THEN
         DO KK = KBAR+1, 0, - 1
            WRITE (9, '(10e12.3)') ((X(II, jj, KK), II=0, IBAR+1), jj=JBAR, 1,-1)
            WRITE (9,*) '----------------------------------------'
         ENDDO
      ELSE
         DO KK = KBAR+1, 0, - 1
            WRITE (9, '(6e12.3)') ((X(II, jj, KK), II=0, IBAR+1), jj=JBAR, 1,-1)
            WRITE (9,*) '----------------------------------------'
         ENDDO
      ENDIF
   ENDIF
   CALL flush (9)
ENDIF
END SUBROUTINE SCARC_COMPARE_SINGLE0
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! only for debugging reasons ... related to single level
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_COMPARE_SINGLE_LEV (X, CROUTINE, INUM, CNAME, ILINE, NM, ILEVEL)
REAL (EB), POINTER, DIMENSION (:, :, :) :: X
INTEGER :: INUM, ILINE, ILEVEL, NM
INTEGER :: II, jj, KK
CHARACTER (4) :: CROUTINE
CHARACTER (6) :: CNAME
 
SL => S%SLEVEL(ILEVEL)

IF (SCARC_DEBUG >= 2) THEN
   IF (TWO_D) THEN
      WRITE (9,*) '============================================================='
      WRITE (9,*) '===   COMPARE= ',CNAME,': SCARC_COUNT2= ',SCARC_COUNT2-1,':ILINE=',ILINE,':INUM=',INUM
      WRITE (9,*) '============================================================='
      IF (NMESHES == 1) THEN
         DO KK = SL%KBAR+1, 0, - 1
            WRITE (9, '(a,i3,a,10e12.3)') 'k=', KK, ' : ', (X(II, 1, KK), II=0,SL%IBAR+1)
         ENDDO
      ELSE
         IF (IBAR==8) THEN
         DO KK = SL%KBAR+1, 0, - 1
            WRITE (9, '(a,i3,a,10e12.3)') 'k=', KK, ' : ', (X(II, 1, KK), II=0, SL%IBAR+1)
         ENDDO
         ELSE
         DO KK = SL%KBAR+1, 0, - 1
            WRITE (9, '(a,i3,a,6e12.3)') 'k=', KK, ' : ', (X(II, 1, KK), II=0, SL%IBAR+1)
         ENDDO
         ENDIF
      ENDIF
   ELSE
      WRITE (9,*) '============================================================='
      WRITE (9,*) '===   COMPARE= ',CNAME,': ROUTINE= ',CROUTINE,': NM=',NM
      WRITE (9,*) '============================================================='
      IF (NMESHES == 1) THEN
         DO KK = SL%KBAR+1, 0, - 1
            WRITE (9, '(10e12.3)') ((X(II, jj, KK), II=0, SL%IBAR+1), jj=SL%JBAR, 1,-1)
            WRITE (9,*) '----------------------------------------'
         ENDDO
      ELSE
         DO KK = SL%KBAR+1, 0, - 1
            WRITE (9, '(6e12.3)') ((X(II, jj, KK), II=0, SL%IBAR+1), jj=SL%JBAR, 1,-1)
            WRITE (9,*) '----------------------------------------'
         ENDDO
      ENDIF
   ENDIF
   CALL flush (9)
ENDIF
 
 
 
END SUBROUTINE SCARC_COMPARE_SINGLE_LEV
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! only for debugging reasons ...
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_COMPARE_SINGLE02(X, CROUTINE, INUM, CNAME, ILINE)
REAL (EB), POINTER, DIMENSION (:, :, :) :: X
INTEGER :: INUM, ILINE
INTEGER :: II, jj, KK
CHARACTER (4) :: CROUTINE
CHARACTER (6) :: CNAME
 
IF (SCARC_DEBUG >= 2) THEN
   IF (TWO_D) THEN
      WRITE (9,*) '============================================================='
      WRITE (9,*) '===   COMPARE= ',CNAME,': SCARC_COUNT2= ',SCARC_COUNT2-1,':ILINE=',ILINE,':INUM=',INUM
      WRITE (9,*) '============================================================='
      IF (NMESHES == 1) THEN
         DO KK = KBAR+1, 0, - 1
            WRITE (9, '(a,i3,a,10e12.3)') 'k=', KK, ' : ', (X(II, 0, KK), II=0,IBAR+1)
         ENDDO
      ELSE
         IF (IBAR==8) THEN
         DO KK = KBAR+1, 0, - 1
            WRITE (9, '(a,i3,a,10e12.3)') 'k=', KK, ' : ', (X(II, 0, KK), II=0, IBAR+1)
         ENDDO
         ELSE
         DO KK = KBAR+1, 0, - 1
            WRITE (9, '(a,i3,a,6e12.3)') 'k=', KK, ' : ', (X(II, 0, KK), II=0, IBAR+1)
         ENDDO
         ENDIF
      ENDIF
      WRITE (9,*) 
      WRITE (9,*) '===   COMPARE= ',CNAME,': SCARC_COUNT2= ',SCARC_COUNT2-1,'J=1'
      IF (NMESHES == 1) THEN
         DO KK = KBAR+1, 0, - 1
            WRITE (9, '(a,i3,a,10e12.3)') 'k=', KK, ' : ', (X(II, 1, KK), II=0,IBAR+1)
         ENDDO
      ELSE
         IF (IBAR==8) THEN
         DO KK = KBAR+1, 0, - 1
            WRITE (9, '(a,i3,a,10e12.3)') 'k=', KK, ' : ', (X(II, 1, KK), II=0, IBAR+1)
         ENDDO
         ELSE
         DO KK = KBAR+1, 0, - 1
            WRITE (9, '(a,i3,a,6e12.3)') 'k=', KK, ' : ', (X(II, 1, KK), II=0, IBAR+1)
         ENDDO
         ENDIF
      ENDIF
      WRITE (9,*) 
      WRITE (9,*) '===   COMPARE= ',CNAME,': SCARC_COUNT2= ',SCARC_COUNT2-1,'J=2'
      IF (NMESHES == 1) THEN
         DO KK = KBAR+1, 0, - 1
            WRITE (9, '(a,i3,a,10e12.3)') 'k=', KK, ' : ', (X(II, 2, KK), II=0,IBAR+1)
         ENDDO
      ELSE
         IF (IBAR==8) THEN
         DO KK = KBAR+1, 0, - 1
            WRITE (9, '(a,i3,a,10e12.3)') 'k=', KK, ' : ', (X(II, 2, KK), II=0, IBAR+1)
         ENDDO
         ELSE
         DO KK = KBAR+1, 0, - 1
            WRITE (9, '(a,i3,a,6e12.3)') 'k=', KK, ' : ', (X(II, 2, KK), II=0, IBAR+1)
         ENDDO
         ENDIF
      ENDIF
            WRITE (9, *) 
   ELSE
      WRITE (9,*) '============================================================='
      WRITE (9,*) '===   COMPARE= ',CNAME,': ROUTINE= ',CROUTINE
      WRITE (9,*) '============================================================='
      IF (NMESHES == 1) THEN
         DO KK = KBAR+1, 0, - 1
            WRITE (9, '(10e12.3)') ((X(II, jj, KK), II=0, IBAR+1), jj=JBAR, 1,-1)
            WRITE (9,*) '----------------------------------------'
         ENDDO
      ELSE
         DO KK = KBAR+1, 0, - 1
            WRITE (9, '(6e12.3)') ((X(II, jj, KK), II=0, IBAR+1), jj=JBAR, 1,-1)
            WRITE (9,*) '----------------------------------------'
         ENDDO
      ENDIF
   ENDIF
   CALL flush (9)
ENDIF
 
END SUBROUTINE SCARC_COMPARE_SINGLE02

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! print out timings for ScaRC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_TIMINGS
INTEGER:: NM, N, IERR, I
INTEGER, ALLOCATABLE, DIMENSION(:) :: COUNTS_SCARC, DISPLS_SCARC, COUNTS_SCARC_TIMERS, DISPLS_SCARC_TIMERS
CHARACTER(25) :: NAME_SCARC(0:N_TIMERS_SCARC)
REAL(EB) :: TPCNT_SCARC(0:N_TIMERS_SCARC)

INCLUDE 'mpif.h'

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
   NAME_SCARC(1)  = 'INITIALIZE_MESHES'
   NAME_SCARC(2)  = 'INITIALIZE_NEIGHBORS'
   NAME_SCARC(3)  = 'INITIALIZE_COMMUNICATION'
   NAME_SCARC(4)  = 'INITIALIZE_MATRICES'
   NAME_SCARC(5)  = 'INITIALIZE_SOLVER'
   NAME_SCARC(6)  = 'CG'
   NAME_SCARC(7)  = 'MG'
   NAME_SCARC(8)  = 'SMOOTHER'
   NAME_SCARC(9)  = 'COARSE'
   NAME_SCARC(10) = 'MATVEC'
   NAME_SCARC(11) = 'RESTRICTION'
   NAME_SCARC(12) = 'PROLONGATION'
   NAME_SCARC(13) = 'BDRY_RESIDUUM'
   NAME_SCARC(14) = 'VECADD'
   NAME_SCARC(15) = 'VECCOPY'
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

443 FORMAT(//' ScaRC: CPU Time Usage, Mesh ',I3// &
         '                 CPU (s)        %  '/ &
         '       --------------------------------------------------')
444 FORMAT(7X,A25,2F11.2)

END SUBROUTINE SCARC_TIMINGS
 
END MODULE SCARC_SOLVER
 
