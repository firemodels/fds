MODULE SCARC_SOLVER
 
USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_POINTERS
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
USE COMP_FUNCTIONS, ONLY: SECOND
USE MPI
 
IMPLICIT NONE

CHARACTER(255), PARAMETER :: scrchid='$Id: mesh.f90 6197 2010-05-11 13:20:05Z mcgratta $'
CHARACTER(255), PARAMETER :: scrcrev='$Revision: 6197 $'
CHARACTER(255), PARAMETER :: scrcdate='$Date: 2010-05-11 09:20:05 -0400 (Tue, 11 May 2010) $'

PRIVATE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Public structures
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! Public subprograms (called from main and pres)
PUBLIC SCARC_INITIALIZE, SCARC_TIMINGS
PUBLIC SCARC_CG2D, SCARC_CG3D, SCARC_MG2D, SCARC_BICG2D, SCARC_GHOSTCELLS, SCARC_UPDATE, GET_REV_scrc
 
!!! Public variables (needed in main, read, pres, divg)
PUBLIC SCARC_METHOD  , SCARC_DEBUG      , SCARC_COMPARE , SCARC_CASE
PUBLIC SCARC_EPS_REL , SCARC_EPS_DIVG   , SCARC_BREL
PUBLIC SCARC_NIT_MG  , SCARC_PRECON_MG  , SCARC_EPS_MG  , SCARC_OMEGA_MG
PUBLIC SCARC_NIT_CG  , SCARC_PRECON_CG  , SCARC_EPS_CG  , SCARC_OMEGA_CG
PUBLIC SCARC_NIT_SM  , SCARC_PRECON_SM  , SCARC_EPS_SM  , SCARC_OMEGA_SM
PUBLIC SCARC_NIT_CO  , SCARC_PRECON_CO  , SCARC_EPS_CO  , SCARC_OMEGA_CO
PUBLIC SCARC_NIT_BICG, SCARC_PRECON_BICG, SCARC_EPS_BICG, SCARC_OMEGA_BICG


CHARACTER (40) :: SCARC_MSG='msg/   _scarc'          ! name of file for ScaRC debug messages
CHARACTER (10) :: SCARC_METHOD='FFT'                 ! name of method for the solution of the pressure equation

INTEGER        :: SCARC_DEBUG=0,     &               ! debug level (0: no debug messages)
                  SCARC_COMPARE=0,   &               ! print out vectors for comparison in separate program
                  SCARC_CASE=0,      &               ! choose different initial solutions
                  SCARC_COUNT=0                      ! counter  for comparison vectors

CHARACTER (10) :: SCARC_PRECON_MG  ='GSTRIX', &      ! smoother for mg-method (default: GSTRIX)
                  SCARC_PRECON_CG  ='SSOR',   &      ! preconditioner for cg-method          (default: SSOR)
                  SCARC_PRECON_SM  ='SSOR',   &      ! ...                smoother           (default: SSOR)
                  SCARC_PRECON_CO  ='SSOR',   &      ! ...                coarse grid solver (default: SSOR)
                  SCARC_PRECON_BICG='SSOR'           ! ...                bicg-method        (default: SSOR)

INTEGER        :: SCARC_NIT_MG   =1000, &            ! max number of iterations for mg
                  SCARC_NIT_CG   =1000, &            ! ...                          cg
                  SCARC_NIT_SM   =1000, &            ! ...                          smoother
                  SCARC_NIT_CO   =1000, &            ! ...                          coarse grid solver
                  SCARC_NIT_BICG =1000               ! ...                          bicg

REAL (EB)      :: SCARC_EPS_MG  =1.E-12_EB, &        ! convergence epsilon for mg
                  SCARC_EPS_CG  =1.E-12_EB, &        ! ...                     cg
                  SCARC_EPS_SM  =1.E-12_EB, &        ! ...                     smoother
                  SCARC_EPS_CO  =1.E-12_EB, &        ! ...                     coarse grid solver
                  SCARC_EPS_BICG=1.E-12_EB           ! ...                     bicg

REAL (EB)      :: SCARC_OMEGA_MG  =1.E+0_EB, &       ! relaxation parameter for mg
                  SCARC_OMEGA_CG  =1.E+0_EB, &       ! ...                      cg
                  SCARC_OMEGA_SM  =1.E+0_EB, &       ! ...                      smoother
                  SCARC_OMEGA_CO  =1.E+0_EB, &       ! ...                      coarse grid solver
                  SCARC_OMEGA_BICG=1.E+0_EB          ! ...                      bicg

REAL (EB)      :: SCARC_EPS_DIVG =1.E+6_EB           ! divergence epsilon for all methods
REAL (EB)      :: SCARC_EPS_REL  =1.E-2_EB           ! minimum relative accuracy for all methods
LOGICAL        :: SCARC_BREL=.FALSE.                  ! relative accuracy activated ? (valid for all methods)


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
                      NCOM_EDGE = 2, &               ! only edge  communication
                      NCOM_DIAG = 3, &               ! only diag communication
                      NCOM_FULL = 4                  ! full face-edge-diag communication

!!! Private global variables 
INTEGER :: NEXCHANGE_EXTENDED = 1
INTEGER :: NNLEVEL = 10                              ! max number of levels currently allowed
INTEGER :: SNODE, RNODE, NDIAG
INTEGER :: NFACE, NEDGE, NVRTX
LOGICAL :: BEDGE, BDIAG                              ! edge and diagonal communication really needed ?
 
INTEGER :: NREQ_FACE, NREQ_EDGE, NREQ_DIAG           ! protocol information for data exchange
INTEGER, ALLOCATABLE, DIMENSION (:)    :: REQ_FACE,  REQ_EDGE,  REQ_DIAG
INTEGER, ALLOCATABLE, DIMENSION (:, :) :: TAGS_FACE, TAGS_EDGE, TAGS_DIAG
INTEGER, ALLOCATABLE, DIMENSION (:, :) :: STAT_FACE, STAT_EDGE, STAT_DIAG
INTEGER, ALLOCATABLE, DIMENSION (:, :) :: NBR_FACE, NBR_EDGE, NBR_DIAG

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
   REAL (EB), POINTER, DIMENSION (:, :, :) :: X, F, D, Y, G, R, Z, X2, F2, D2, TMP
   REAL (EB), POINTER, DIMENSION (:, :)    :: BXS0, BXF0, BYS0, BYF0, BZS0, BZF0
   REAL (EB) :: ASUBX, ASUBY, ASUBZ, ASUB
 
   ! communication vectors
   REAL (EB), POINTER, DIMENSION (:, :, :) :: Y_MATV, R_MATV, G_MATV
   REAL (EB), POINTER, DIMENSION (:, :, :) :: Y_FACE, R_FACE, G_FACE
   REAL (EB), POINTER, DIMENSION (:, :, :) :: Y_EDGE, R_EDGE, G_EDGE
   REAL (EB), POINTER, DIMENSION (:, :, :) :: Y_DIAG, R_DIAG, G_DIAG
   REAL (EB), POINTER, DIMENSION (:, :, :) :: Z_FACE, Z_EDGE, Z_DIAG
 
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
 

!------------------------------------------------------------------------------------------
! Scarc type on 'own' mesh
!------------------------------------------------------------------------------------------
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
   INTEGER :: NFACE0, NEDGE0, NVRTX0
   INTEGER :: MIBAR_MIN, MJBAR_MIN, MKBAR_MIN, MESH_MIN

   TYPE (SCARC_LEVEL_TYPE), POINTER, DIMENSION (:) :: SLEVEL
 
END TYPE SCARC_TYPE


!------------------------------------------------------------------------------------------
! Scarc type on 'other' mesh
!------------------------------------------------------------------------------------------
TYPE OSCARC_TYPE
 
   REAL (EB), POINTER, DIMENSION (:) :: SEND_FACE, RECV_FACE
   REAL (EB), POINTER, DIMENSION (:) :: SEND_EDGE, RECV_EDGE
   REAL (EB), POINTER, DIMENSION (:) :: SEND_DIAG, RECV_DIAG

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
IF (SCARC_DEBUG >= 2) THEN
   WRITE (SCARC_MSG(5:7), '(i3.3)') NM
   OPEN (9, FILE=SCARC_MSG)
ENDIF

IF (SCARC_METHOD=='FFT') GOTO 12345

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
IF (SCARC_METHOD == 'MG' .OR. SCARC_PRECON_CG == 'MG') THEN       ! mg-structures with complete grid hierarchy 
   S%NLEVEL_MIN  = S%MLEVEL(NM)-S%NLEVEL+1
ELSE                                                              ! only cg-structures on one grid level
   S%NLEVEL = 1
   S%NLEVEL_MIN  = S%NLEVEL_MAX
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

   CALL SCARC_INIT_MESHES2D (NM)
   CALL SCARC_INIT_NEIGHBORS2D(NM)
   CALL SCARC_INIT_MATRICES2D (NM)
   CALL SCARC_INIT_COMMUNICATION2D (NM)
   CALL SCARC_INIT_SOLVER2D (NM)

ELSE

   CALL SCARC_INIT_MESHES3D (NM)         
   write(*,*) '3D-ScaRC still under construction, do not use at the moment !!'
   stop
   !CALL SCARC_INIT_NEIGHBORS3D (NM)
   CALL SCARC_INIT_MATRICES3D (NM)
   !CALL SCARC_INITIALIZE_MESH_EXCHANGE3D (NM)
   CALL SCARC_INIT_SOLVER3D (NM)

ENDIF
 
12345 CONTINUE
IF (SCARC_DEBUG>=2) WRITE(9,*) 'Leaving scarc_initialize'

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
INTEGER :: IBAR0, JBAR0, KBAR0, IZERO=0, IERR=0
TYPE (MESH_TYPE), POINTER :: M
REAL(EB):: TNOW_MESHES2D

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
   CALL CHKMEMERR ('SCARC_INIT_MESHES2D', 'SL%XX', IZERO)
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
   CALL CHKMEMERR ('SCARC_INIT_MESHES2D', 'SL%YY', IZERO)
   SL%YY(0) =0.0_EB
   SL%YY_MIN=0.0_EB
   SL%YY_MAX=0.0_EB
   SL%DY_MIN=0.0_EB
   SL%DY_MAX=0.0_EB

   ALLOCATE (SL%ZZ(0:SL%KBAR), STAT=IZERO)
   CALL CHKMEMERR ('SCARC_INIT_MESHES2D', 'SL%ZZ', IZERO)
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
      WRITE(9,*) ' ======================= SCARC_INIT_MESHES2D =========='
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
      !CALL flush (9)
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
INTEGER :: IBAR0, JBAR0, KBAR0, IZERO=0, IERR=0
TYPE (MESH_TYPE), POINTER :: M
REAL(EB):: TNOW_MESHES3D

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
   CALL CHKMEMERR ('SCARC_INIT_MESHES3D', 'SL%XX', IZERO)
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
   CALL CHKMEMERR ('SCARC_INIT_MESHES3D', 'SL%YY', IZERO)
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
   CALL CHKMEMERR ('SCARC_INIT_MESHES3D', 'SL%ZZ', IZERO)
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
      WRITE(9,*) ' ======================= SCARC_INIT_MESHES3D =========='
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
      !CALL flush (9)
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
 
TNOW_NEIGHBORS2D = SECOND()

IF (SCARC_DEBUG >= 2) THEN
   WRITE (9,*) '==========================================================================='
   WRITE (9,*) '=== SCARC_INIT_NEIGHBORS2D ',NM
   WRITE (9,*) '==========================================================================='
   WRITE (9,*) 'NMESHES=', NMESHES
   WRITE (9,*) 'S%NLEVEL=', S%NLEVEL
   WRITE (9,*) 'HIER SCHAUEN, WELCHES GEBRAUCHT WIRD'
   WRITE (9,*) 'NEWC=', M%NEWC
   WRITE (9,*) 'NWC=', M%NWC
   !CALL flush (9)
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
   !flush(9)
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
   !flush(9)
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
   !flush(9)
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
   !call flush(9)
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
            !flush(9)
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
   !flush(9)
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
   !flush(9)
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Only in case of MG-method:
!!! Determine arrays IJKW_FACE and IJKW_DIAG for coarser levels
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (SCARC_METHOD == 'MG' .OR. SCARC_PRECON_CG =='MG') THEN

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
         !flush(9)
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
         !flush(9)
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

END SUBROUTINE SCARC_INIT_NEIGHBORS2D

 
!------------------------------------------------------------------------------------
! determine neighborship/communication structure for data exchange in 2D
!------------------------------------------------------------------------------------
SUBROUTINE SCARC_INIT_COMMUNICATION2D (NM)
 
INTEGER :: IM, NM, NOM, ILEVEL, IREFINE
INTEGER :: IMIN, IMAX, JMIN, JMAX, KMIN, KMAX, ILMAX
INTEGER :: TAG_FACE, TAG_DIAG
INTEGER :: IZERO
REAL(EB):: TNOW_COMMUNICATION2D

TNOW_COMMUNICATION2D = SECOND()
 
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
      CALL CHKMEMERR ('SCARC_INIT_COMMUNICATION2D', 'OS%SLEVEL', IZERO)

      OSLMAX => OS%SLEVEL(ILMAX)
      OSLMAX%NW_FACE = 2*S%MIBAR(NOM) + 2*S%MKBAR(NOM)
      ALLOCATE (OSLMAX%IJKW_FACE(17,OSLMAX%NW_FACE), STAT=IZERO)
      CALL CHKMEMERR ('SCARC_INIT_COMMUNICATION2D', 'OSLMAX%IJKW_FACE', IZERO)
      OSLMAX%IJKW_FACE = 0

      IF (S%NLEVEL_DIFF/=0) THEN
         IREFINE=1
         DO ILEVEL=S%NLEVEL_MAX-1,S%NLEVEL_MIN

            IREFINE=IREFINE*2

            OS%SLEVEL(ILEVEL)%NW_FACE=OS%SLEVEL(ILEVEL+1)%NW_FACE/IREFINE
            ALLOCATE (OS%SLEVEL(ILEVEL)%IJKW_FACE(17,OSLMAX%NW_FACE), STAT=IZERO)
            CALL CHKMEMERR ('SCARC_INIT_COMMUNICATION2D', 'OS%SLEVEL%IJKW_FACE', IZERO)
            OS%SLEVEL(ILEVEL)%IJKW_FACE = 0

         ENDDO
      ENDIF

   ENDIF

   IF (SLMAX%NIC_DIAG(NM,NOM)/=0.AND.SLMAX%NIC_DIAG(NOM,NM)/=0) THEN

      OS => S%OSCARC(NOM)
      ALLOCATE (OS%SLEVEL(S%NLEVEL_MIN:S%NLEVEL_MAX), STAT=IZERO)
      CALL CHKMEMERR ('SCARC_INIT_COMMUNICATION2D', 'OS%SLEVEL', IZERO)

      OSLMAX => OS%SLEVEL(ILMAX)
      OSLMAX%NW_DIAG = 4
      ALLOCATE (OSLMAX%IJKW_DIAG(17,OSLMAX%NW_DIAG), STAT=IZERO)
      CALL CHKMEMERR ('SCARC_INIT_COMMUNICATION2D', 'OSLMAX%IJKW_DIAG', IZERO)
      OSLMAX%IJKW_DIAG = 0

      IF (S%NLEVEL_DIFF/=0) THEN
         IREFINE=1
         DO ILEVEL=S%NLEVEL_MAX-1,S%NLEVEL_MIN

            IREFINE=IREFINE*2

            OS%SLEVEL(ILEVEL)%NW_DIAG=OS%SLEVEL(ILEVEL+1)%NW_DIAG/IREFINE
            ALLOCATE (OS%SLEVEL(ILEVEL)%IJKW_DIAG(17,OSLMAX%NW_DIAG), STAT=IZERO)
            CALL CHKMEMERR ('SCARC_INIT_COMMUNICATION2D', 'OS%SLEVEL%IJKW_DIAG', IZERO)
            OS%SLEVEL(ILEVEL)%IJKW_DIAG = 0

         ENDDO
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
         !CALL flush (9)
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
       ENDIF

   ENDDO DIAG_NBR_LOOP2D

ENDDO COM_LEVEL_LOOP2D
 
DO ILEVEL=S%NLEVEL_MAX,S%NLEVEL_MIN,-1
   CALL SCARC_RECEIVE (NCOM_INIT, ILEVEL)    
   CALL SCARC_EXCHANGE(NCOM_INIT, ILEVEL, NCOM_TYPE_NONE)
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
INTEGER :: I, K, IC, IW, IZERO, IERR, ICELL, IOR, ILEVEL
REAL(EB):: DBC
REAL(EB):: TNOW_MATRICES2D
TYPE (MESH_TYPE), POINTER :: M
 
TNOW_MATRICES2D = SECOND()

IERR = 0
 
M => MESHES(NM)
 
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
   CALL CHKMEMERR ('SCARC_INIT_MATRICES2D', 'SL%AG', IZERO)
   SL%AG = 0.0_EB
    
   SL%DXI2=1.0_EB/(SL%DX)**2
   SL%DZI2=1.0_EB/(SL%DZ)**2

   IF (SCARC_DEBUG >= 2) THEN
      WRITE (9,*) '========= MATRICES2D: LEVEL=',ILEVEL
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
    
   IF (ILEVEL/=S%NLEVEL_MIN) CALL SCARC_PRECON_INIT_GSTRIX2D(ILEVEL,NM)
   IF (ILEVEL==S%NLEVEL_MAX) THEN

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
            IF (SCARC_DEBUG .GE. 4) write(9,1001) IW, IOR, IC, SL%AG(IC,1),M%BOUNDARY_TYPE(IW),NOM
         CASE (INTERNAL)                     ! do nothing along internal boundaries (just only debug message)
            IF (SCARC_DEBUG .GE. 4) write(9,1003) IW, IOR, IC, SL%AG(IC,1),M%BOUNDARY_TYPE(IW),NOM
         CASE (NEUMANN)                      ! set Neumann BC's at all other nodes
            SL%AG(IC, 1) = SL%AG(IC, 1) + DBC
            IF (SCARC_DEBUG .GE. 4) write(9,1002) IW, IOR, IC, SL%AG(IC,1),M%BOUNDARY_TYPE(IW),NOM
      END SELECT

   ENDDO
   ENDIF

   SL%ASUBX = SL%DXI2
   SL%ASUBZ = SL%DZI2
 
   IF (SCARC_DEBUG >= 2) THEN
      !WRITE(9,'(16i4)') (M%BOUNDARY_TYPE(ICELL),ICELL=1,M%NEWC)
      !WRITE(9,*) 'IJKW'
      !DO IW=1,M%NEWC+4
      !   WRITE(9,'(i4,a,15i4)') IW,' : ',(M%IJKW(ICELL,IW),ICELL=1,15)
      !ENDDO
      WRITE (9,*) 'ASUB=', SL%ASUB
      WRITE (9,*)
      WRITE (9,*) 'GLOBAL MATRIX AG:'
      DO IC = 1, SL%NCELLS_LOCAL
         WRITE (9, '(5f12.4)') SL%AG(IC, 2),SL%AG(IC,3),SL%AG(IC,1),SL%AG(IC,4),SL%AG(IC,5)
         IF (Mod(IC, SL%IBAR) == 0) write (9,*) '---------------------------------&
                                                &--------------------------------'
      ENDDO
      write (9,*) 'SL%ASUBX=', SL%ASUBX
      WRITE (9,*) 'LBC=', LBC
      WRITE (9,*) 'MBC=', MBC
      WRITE (9,*) 'NBC=', NBC
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
INTEGER :: I, J, K, IC, IW, IZERO, IERR, ICELL, IOR, ILEVEL
REAL(EB):: DBC
REAL(EB):: TNOW_MATRICES3D
TYPE (MESH_TYPE), POINTER :: M
 
TNOW_MATRICES3D = SECOND()

IERR = 0
 
M => MESHES(NM)
 
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
   CALL CHKMEMERR ('SCARC_INIT_MATRICES3D', 'SL%AG', IZERO)
   SL%AG = 0.0_EB
    
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
    
   IF (ILEVEL/=S%NLEVEL_MIN) CALL SCARC_PRECON_INIT_GSTRIX3D(ILEVEL,NM)
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
            IF (SCARC_DEBUG .GE. 4) write(9,1001) IW, IOR, IC, SL%AG(IC,1),M%BOUNDARY_TYPE(IW),NOM
         CASE (INTERNAL)                     ! do nothing along internal boundaries (only debug message)
            IF (SCARC_DEBUG .GE. 4) write(9,1003) IW, IOR, IC, SL%AG(IC,1),M%BOUNDARY_TYPE(IW),NOM
         CASE (NEUMANN)                      ! set Neumann BC's for AG and AL at all other nodes
            SL%AG(IC, 1) = SL%AG(IC, 1) + DBC
            IF (SCARC_DEBUG .GE. 4) write(9,1002) IW, IOR, IC, SL%AG(IC,1),M%BOUNDARY_TYPE(IW),NOM
      END SELECT

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
 
INTEGER :: NM, IZERO, ILEVEL
INTEGER :: IBP0, KBP0
REAL(EB):: TNOW_SOLVER2D
 
TNOW_SOLVER2D = SECOND()
 
IZERO = 0

S => SCARC (NM)

!!!
!!! Initialize SLEVEL-structures for single levels
!!!
LEVEL_LOOP2D: DO ILEVEL=S%NLEVEL_MAX,S%NLEVEL_MIN,-1

   SL => S%SLEVEL(ILEVEL)

   IBP0=SL%IBAR+1
   KBP0=SL%KBAR+1

   ! ----------------------------------------------------------------------------
   ! working and auxiliary vectors for global CG/MG-method
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
    
   IF ((SCARC_METHOD == 'CG' .OR. SCARC_METHOD == 'BICG').OR.ILEVEL==S%NLEVEL_MIN) THEN
      ALLOCATE (SL%R(0:IBP0, 1, 0:KBP0), STAT=IZERO)
      CALL CHKMEMERR ('SCARC', 'R', IZERO)
      SL%R = 0.0_EB
    
      ALLOCATE (SL%G(0:IBP0, 1, 0:KBP0), STAT=IZERO)
      CALL CHKMEMERR ('SCARC', 'G', IZERO)
      SL%G = 0.0_EB
    
      ALLOCATE (SL%Y(0:IBP0, 1, 0:KBP0), STAT=IZERO)
      CALL CHKMEMERR ('SCARC', 'Y', IZERO)
      SL%Y = 0.0_EB
   ENDIF
    
   IF (SCARC_METHOD == 'BICG') THEN
      ALLOCATE (SL%Z(0:IBP0, 1, 0:KBP0), STAT=IZERO)
      CALL CHKMEMERR ('SCARC', 'Z', IZERO)
      SL%Z = 0.0_EB
   ENDIF
    
   IF (SCARC_PRECON_CG == 'MG') THEN

      ALLOCATE (SL%X2(0:IBP0, 1, 0:KBP0), STAT=IZERO)
      CALL CHKMEMERR ('SCARC', 'X2', IZERO)
      SL%X2 = 0.0_EB

      ALLOCATE (SL%F2(0:IBP0, 1, 0:KBP0), STAT=IZERO)
      CALL CHKMEMERR ('SCARC', 'F2', IZERO)
      SL%F2 = 0.0_EB

      ALLOCATE (SL%D2(0:IBP0, 1, 0:KBP0), STAT=IZERO)
      CALL CHKMEMERR ('SCARC', 'D2', IZERO)
      SL%D2 = 0.0_EB
    
   ENDIF

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
IF (SCARC_METHOD == 'MG' .OR. SCARC_PRECON_CG =='MG') THEN
   ALLOCATE (S%KCYCLE(2, NNLEVEL), STAT=IZERO)
   CALL CHKMEMERR ('SCARC', 'KCYCLE', IZERO)
   S%KCYCLE = 0
ENDIF
    
IF (SCARC_DEBUG >= 2) write (9,*) 'INIT_SOLVER2D: FINISHED'

TUSED_SCARC(5,NM)=TUSED_SCARC(5,NM)+SECOND()-TNOW_SOLVER2D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_SOLVER2D
 
END SUBROUTINE SCARC_INIT_SOLVER2D

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Initialize global 3D-solver methods (cg/mg)
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_INIT_SOLVER3D (NM)
 
INTEGER :: NM, IZERO, ILEVEL
INTEGER :: IBP0, JBP0, KBP0
REAL(EB):: TNOW_SOLVER3D
 
TNOW_SOLVER3D = SECOND()
 
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
IF (SCARC_METHOD == 'MG' .OR. SCARC_PRECON_CG =='MG') THEN
   ALLOCATE (S%KCYCLE(2, NNLEVEL), STAT=IZERO)
   CALL CHKMEMERR ('SCARC', 'KCYCLE', IZERO)
   S%KCYCLE = 0
ENDIF
    
IF (SCARC_DEBUG >= 2) write (9,*) 'INIT_SOLVER: FINISHED'
 
TUSED_SCARC(5,NM)=TUSED_SCARC(5,NM)+SECOND()-TNOW_SOLVER3D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_SOLVER3D

END SUBROUTINE SCARC_INIT_SOLVER3D
 

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Perform global 2D cg-method based on global possion-matrix
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_CG2D (NM)
 
INTEGER ::  NM, ILMAX, ITE, ITE0, I, K
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
ILMAX =  S%NLEVEL_MAX                              
SL    => S%SLEVEL(ILMAX)
 
! get right hand side 
DO K = 1, KBAR                                     ! 
   DO I = 1, IBAR
      SL%F (I, 1, K) = M%PRHS (I, 1, K)
   ENDDO
ENDDO

! re-initialize auxiliary vectors
SL%G = 0.0_EB                   
SL%Y = 0.0_EB
SL%R = 0.0_EB
SL%D = 0.0_EB
SL%X = 0.0_EB
 
! set boundary conditions along exterior boundaries
CALL SCARC_SETBDRY2D (ILMAX, NM)
 
!!! calculate initial defect Y = B - A*X and get l2-norm of it
CALL SCARC_MATVEC2D (SL%X, SL%R, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILMAX, 3)
CALL SCARC_VECADD2D (SL%F, SL%R,-1.0_EB, 1.0_EB, ILMAX, NM)
CALL SCARC_L2NORM2D (SL%R, SL%RESIN_CG, ILMAX, NM, NTYPE_GLOBAL)
IF (SCARC_DEBUG >= 1) write (9,1000) 0, ILMAX, SL%RESIN_CG
IF (SCARC_DEBUG >= 1) write (*,1000) 0, ILMAX, SL%RESIN_CG
 
CALL SCARC_COMPARE_SINGLE_LEV (SL%X, 'SARC', 'X0    ',ILMAX)
CALL SCARC_COMPARE_SINGLE_LEV (SL%F, 'SARC', 'F0    ',ILMAX)
CALL SCARC_COMPARE_SINGLE_LEV (SL%R, 'SARC', 'R0    ',ILMAX)

! initial preconditioning
SELECT CASE(SCARC_PRECON_CG)
   CASE('NONE') 
      SIGMA0 = SL%RESIN_CG * SL%RESIN_CG
   CASE('JACOBI') 
      CALL SCARC_VECCOPY2D (SL%R, SL%G, ILMAX, NM)
      CALL SCARC_PRECON_JACOBI2D (SL%G, ILMAX, NM)
      CALL SCARC_SCALPROD2D (SL%R, SL%G, SIGMA0, ILMAX, NM)
   CASE('SSOR') 
      CALL SCARC_VECCOPY2D (SL%R, SL%G, ILMAX, NM)
      CALL SCARC_PRECON_SSOR2D (SL%G, ILMAX, NM)
      CALL SCARC_SCALPROD2D (SL%R, SL%G, SIGMA0, ILMAX, NM)
   CASE('GSTRIX') 
      CALL SCARC_VECCOPY2D (SL%R, SL%G, ILMAX, NM)
      CALL SCARC_PRECON_GSTRIX2D (SL%G, ILMAX, NM)
      CALL SCARC_SCALPROD2D (SL%R, SL%G, SIGMA0, ILMAX, NM)
   CASE('MG') 
      CALL SCARC_VECCOPY2D (SL%R, SL%G, ILMAX, NM)
      CALL SCARC_PRECON_MG2D (SL%G, ILMAX, NM)
      CALL SCARC_SCALPROD2D (SL%R, SL%G, SIGMA0, ILMAX, NM)
END SELECT
CALL SCARC_VECADD2D (SL%G, SL%D,-1.0_EB, 0.0_EB, ILMAX, NM)
 
CALL SCARC_COMPARE_SINGLE_LEV (SL%D, 'SARC',  'D0    ',ILMAX)
!
! start defect correction loop
!
CG_LOOP2D: DO ITE = 1, SCARC_NIT_CG
 
   ! calculate new defect and get L2-norm of it
   CALL SCARC_MATVEC2D (SL%D, SL%Y, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILMAX, 5)
   CALL SCARC_SCALPROD2D (SL%D, SL%Y, ALPHA, ILMAX, NM)

   ! get new descent direction
   ALPHA = SIGMA0 / ALPHA
   CALL SCARC_VECADD2D (SL%D, SL%X, ALPHA, 1.0_EB, ILMAX, NM)
   CALL SCARC_VECADD2D (SL%Y, SL%R, ALPHA, 1.0_EB, ILMAX, NM)
   CALL SCARC_L2NORM2D (SL%R, SL%RES_CG, ILMAX, NM, NTYPE_GLOBAL)
 
CALL SCARC_COMPARE_SINGLE_LEV (SL%R, 'SARC',  'R1    ',ILMAX)
   !  exit loop in case of convergence or divergence
   IF (SCARC_DEBUG >= 1) write (9,1000) ITE, ILMAX, SL%RES_CG
   IF (SCARC_DEBUG >= 1) write (*,1000) ITE, ILMAX, SL%RES_CG
   IF (NM == 1)          write (*,1000) ITE, ILMAX, SL%RES_CG
   IF (SCARC_BREL) THEN
      IF (SL%RES_CG <= SL%RESIN_CG*SCARC_EPS_CG) BCONV = .TRUE.
   ELSE
      IF (SL%RES_CG <= SCARC_EPS_CG .AND. SL%RES_CG <= SL%RESIN_CG*SCARC_EPS_REL) BCONV = .TRUE.
   ENDIF
   IF (SL%RES_CG > SCARC_EPS_DIVG) BDIVG = .TRUE.
   IF (BCONV.OR.BDIVG) EXIT CG_LOOP2D
 
   ! preconditioning
   SELECT CASE(SCARC_PRECON_CG)
      CASE('NONE')
         SIGMA1 = SL%RES_CG * SL%RES_CG
      CASE('JACOBI')
         CALL SCARC_VECCOPY2D (SL%R, SL%G, ILMAX, NM)
         CALL SCARC_PRECON_JACOBI2D (SL%G, ILMAX, NM)
         CALL SCARC_SCALPROD2D (SL%R, SL%G, SIGMA1, ILMAX, NM)
      CASE('SSOR')
         CALL SCARC_VECCOPY2D (SL%R, SL%G, ILMAX, NM)
         CALL SCARC_PRECON_SSOR2D (SL%G, ILMAX, NM)
         CALL SCARC_SCALPROD2D (SL%R, SL%G, SIGMA1, ILMAX, NM)
      CASE('GSTRIX')
         CALL SCARC_VECCOPY2D (SL%R, SL%G, ILMAX, NM)
         CALL SCARC_PRECON_GSTRIX2D (SL%G, ILMAX, NM)
         CALL SCARC_SCALPROD2D (SL%R, SL%G, SIGMA1, ILMAX, NM)
      CASE('MG')
         CALL SCARC_VECCOPY2D (SL%R, SL%G, ILMAX, NM)
         CALL SCARC_PRECON_MG2D (SL%G, ILMAX, NM)
         CALL SCARC_SCALPROD2D (SL%R, SL%G, SIGMA1, ILMAX, NM)
   END SELECT
 
   GAMMA = SIGMA1 / SIGMA0
   SIGMA0 = SIGMA1
   CALL SCARC_VECADD2D (SL%G, SL%D,-1.0_EB, GAMMA, ILMAX, NM)
 
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
IF (SCARC_DEBUG >= 1) WRITE (9,2000) ITE0, ILMAX, SL%RES_CG, SL%CAPPA_CG 
IF (NM == 1)          WRITE (*,2000) ITE0, ILMAX, SL%RES_CG, SL%CAPPA_CG

! shift solution back to vector PRHS
DO K = 1, KBAR
   DO I = 1, IBAR
      M%PRHS (I, 1, K) = SL%X(I, 1, K)
   ENDDO
ENDDO
 
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
SUBROUTINE SCARC_BICG2D (NM)
 
INTEGER ::  NM, ILMAX, ITE, ITE0, I, K
REAL (EB) :: RHO0, RHO1, DBETA, DTHETA, ALPHA0, ALPHA1, ALPHA2
REAL (EB) :: TNOW_BICG2D
LOGICAL BCONV, BDIVG
TYPE (MESH_TYPE),  POINTER :: M
 
TNOW_BICG2D = SECOND()
BCONV=.FALSE.
BDIVG=.FALSE.

S  => SCARC (NM)
M  => MESHES (NM)

IF (SCARC_DEBUG >=2) WRITE(9,*) 'STARTING BICG2D'

! cg only works in finest grid level
ILMAX =  S%NLEVEL_MAX                              
SL    => S%SLEVEL(ILMAX)
 
! get right hand side 
DO K = 1, KBAR                                     ! 
   DO I = 1, IBAR
      SL%F (I, 1, K) = M%PRHS (I, 1, K)
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
SL%X = 0.0_EB
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! initialization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CALL SCARC_SETBDRY2D (ILMAX, NM)
CALL SCARC_VECCOPY2D (SL%F, SL%R, ILMAX, NM)

SELECT CASE(SCARC_PRECON_BICG)
   CASE('JACOBI') 
      CALL SCARC_PRECON_JACOBI2D (SL%R, ILMAX, NM)
   CASE('SSOR') 
      CALL SCARC_PRECON_SSOR2D (SL%R, ILMAX, NM)
   CASE('GSTRIX') 
      CALL SCARC_PRECON_GSTRIX2D (SL%R, ILMAX, NM)
   CASE('MG') 
      CALL SCARC_PRECON_MG2D (SL%R, ILMAX, NM)
END SELECT

CALL SCARC_MATVEC2D (SL%X, SL%R, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILMAX, 3)
CALL SCARC_VECADD2D (SL%F, SL%R,1.0_EB, -1.0_EB, ILMAX, NM)

CALL SCARC_COMPARE_SINGLE (SL%R, 'SARC', 'R0    ')

SELECT CASE(SCARC_PRECON_BICG)
   CASE('JACOBI') 
      CALL SCARC_PRECON_JACOBI2D (SL%R, ILMAX, NM)
   CASE('SSOR') 
      CALL SCARC_PRECON_SSOR2D (SL%R, ILMAX, NM)
   CASE('GSTRIX') 
      CALL SCARC_PRECON_GSTRIX2D (SL%R, ILMAX, NM)
   CASE('MG') 
      CALL SCARC_PRECON_MG2D (SL%R, ILMAX, NM)
END SELECT

CALL SCARC_L2NORM2D (SL%R, SL%RESIN_BICG, ILMAX, NM, NTYPE_GLOBAL)
CALL SCARC_VECCOPY2D (SL%R, SL%G, ILMAX, NM)

IF (SCARC_DEBUG >= 1) write (9,1000) 0, ILMAX, SL%RESIN_BICG
IF (SCARC_DEBUG >= 1) write (*,1000) 0, ILMAX, SL%RESIN_BICG

CALL SCARC_COMPARE_SINGLE (SL%G, 'SARC', 'G0    ')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! iterative correction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
BICG_LOOP2D: DO ITE = 1, SCARC_NIT_BICG
 
   CALL SCARC_SCALPROD2D (SL%G, SL%R, RHO1, ILMAX, NM)
   DBETA=(RHO1*DTHETA)/(RHO0*ALPHA0)

   RHO0=RHO1

   CALL SCARC_VECADD2D (SL%R, SL%Z, 1.0_EB, DBETA, ILMAX, NM)
   CALL SCARC_VECADD2D (SL%Y, SL%Z, -DBETA*ALPHA0, 1.0_EB, ILMAX, NM)

   ! matrix-vector multiplication and preconditioning
   CALL SCARC_MATVEC2D (SL%Z, SL%Y, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILMAX, 5)

   SELECT CASE(SCARC_PRECON_BICG)
      CASE('JACOBI') 
         CALL SCARC_PRECON_JACOBI2D (SL%Y, ILMAX, NM)
      CASE('SSOR') 
         CALL SCARC_PRECON_SSOR2D (SL%Y, ILMAX, NM)
      CASE('GSTRIX') 
         CALL SCARC_PRECON_GSTRIX2D (SL%Y, ILMAX, NM)
      CASE('MG') 
         CALL SCARC_PRECON_MG2D (SL%Y, ILMAX, NM)
   END SELECT

   CALL SCARC_SCALPROD2D (SL%G, SL%Y, DTHETA, ILMAX, NM)
   DTHETA=RHO1/DTHETA
   CALL SCARC_VECADD2D (SL%Y, SL%R, -DTHETA, 1.0_EB, ILMAX, NM)

CALL SCARC_COMPARE_SINGLE (SL%R, 'SARC',  'R1    ')
   ! matrix-vector multiplication and preconditioning
   CALL SCARC_MATVEC2D (SL%R, SL%D, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILMAX, 5)

   SELECT CASE(SCARC_PRECON_BICG)
      CASE('JACOBI') 
         CALL SCARC_PRECON_JACOBI2D (SL%D, ILMAX, NM)
      CASE('SSOR') 
         CALL SCARC_PRECON_SSOR2D (SL%D, ILMAX, NM)
      CASE('GSTRIX') 
         CALL SCARC_PRECON_GSTRIX2D (SL%D, ILMAX, NM)
      CASE('MG') 
         CALL SCARC_PRECON_MG2D (SL%D, ILMAX, NM)
   END SELECT

   CALL SCARC_SCALPROD2D (SL%D, SL%R, ALPHA1, ILMAX, NM)
   CALL SCARC_SCALPROD2D (SL%D, SL%D, ALPHA2, ILMAX, NM)

   ALPHA0=ALPHA1/ALPHA2

   CALL SCARC_VECADD2D (SL%Z, SL%X,  DTHETA, 1.0_EB, ILMAX, NM)
   CALL SCARC_VECADD2D (SL%R, SL%X,  ALPHA0, 1.0_EB, ILMAX, NM)
   CALL SCARC_VECADD2D (SL%D, SL%R, -ALPHA0, 1.0_EB, ILMAX, NM)

CALL SCARC_COMPARE_SINGLE (SL%R, 'SARC',  'R2    ')
   CALL SCARC_L2NORM2D (SL%R, SL%RES_BICG, ILMAX, NM, NTYPE_GLOBAL)
 
   !  exit loop in case of convergence or divergence
   IF (SCARC_DEBUG >= 1) write (9,1000) ITE, ILMAX, SL%RES_BICG
   IF (SCARC_DEBUG >= 1) write (*,1000) ITE, ILMAX, SL%RES_BICG
   !IF (NM == 1)          write (*,1000) ITE, ILMAX, SL%RES_BICG
   IF (SCARC_BREL) THEN
      IF (SL%RES_BICG <= SL%RESIN_BICG*SCARC_EPS_BICG) BCONV = .TRUE.
   ELSE
      IF (SL%RES_BICG <= SCARC_EPS_BICG .AND. SL%RES_BICG <= SL%RESIN_BICG*SCARC_EPS_REL) BCONV = .TRUE.
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
IF (SCARC_DEBUG >= 1) WRITE (9,2000) ITE0, ILMAX, SL%RES_BICG, SL%CAPPA_BICG 
IF (NM == 1)          WRITE (*,2000) ITE0, ILMAX, SL%RES_BICG, SL%CAPPA_BICG

! shift solution back to vector PRHS
DO K = 1, KBAR
   DO I = 1, IBAR
      M%PRHS (I, 1, K) = SL%X(I, 1, K)
   ENDDO
ENDDO
 
TUSED_SCARC(6,NM)=TUSED_SCARC(6,NM)+SECOND()-TNOW_BICG2D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_BICG2D

1000 FORMAT ('=====SCARC_BICG     : #ite= ',i4,': level=',i4,': res=',e12.4)
2000 FORMAT ('=====SCARC_BICG     : #ite= ',i4,': level=',i4,': res=',e12.4,': rate=',e12.4)
END SUBROUTINE SCARC_BICG2D


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Perform global 3D cg-method based on global possion-matrix
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_CG3D (NM)
 
INTEGER ::  NM, ILMAX, ITE, ITE0, I, K
REAL (EB) :: SIGMA0, SIGMA1, ALPHA, GAMMA
REAL (EB) :: TNOW_CG3D
LOGICAL BCONV, BDIVG
TYPE (MESH_TYPE),  POINTER :: M
 
TNOW_CG3D = SECOND()
BCONV=.FALSE.
BDIVG=.FALSE.

NREQ_FACE=0                   ! check if really necessary ...
NREQ_EDGE=0
NREQ_DIAG=0

S  => SCARC (NM)
M  => MESHES (NM)

! cg only works in finest grid level
ILMAX =  S%NLEVEL_MAX                              
SL    => S%SLEVEL(ILMAX)
 
! get right hand side 
DO K = 1, KBAR                                     ! 
   DO I = 1, IBAR
      SL%F (I, 1, K) = M%PRHS (I, 1, K)
   ENDDO
ENDDO

! re-initialize auxiliary vectors
SL%G = 0.0_EB                   
SL%Y = 0.0_EB
SL%R = 0.0_EB
SL%D = 0.0_EB
SL%X = 0.0_EB
 
! set boundary conditions along exterior boundaries
CALL SCARC_SETBDRY3D (ILMAX, NM)
 
!!! calculate initial defect Y = B - A*X and get l2-norm of it
CALL SCARC_MATVEC3D (SL%X, SL%R, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILMAX, 3)
CALL SCARC_VECADD3D (SL%F, SL%R,-1.0_EB, 1.0_EB, ILMAX, NM)
CALL SCARC_L2NORM3D (SL%R, SL%RESIN_CG, ILMAX, NM, NTYPE_GLOBAL)
IF (SCARC_DEBUG >= 1) write (9,1000) 0, ILMAX, SL%RESIN_CG
IF (SCARC_DEBUG >= 1) write (*,1000) 0, ILMAX, SL%RESIN_CG
 
! initial preconditioning
SELECT CASE(SCARC_PRECON_CG)
   CASE('NONE') 
      SIGMA0 = SL%RESIN_CG * SL%RESIN_CG
   CASE('JACOBI') 
      CALL SCARC_VECCOPY3D (SL%R, SL%G, ILMAX, NM)
      CALL SCARC_PRECON_JACOBI3D (SL%G, ILMAX, NM)
      CALL SCARC_SCALPROD3D (SL%R, SL%G, SIGMA0, ILMAX, NM)
   CASE('SSOR') 
      CALL SCARC_VECCOPY3D (SL%R, SL%G, ILMAX, NM)
      CALL SCARC_PRECON_SSOR3D (SL%G, ILMAX, NM)
      CALL SCARC_SCALPROD3D (SL%R, SL%G, SIGMA0, ILMAX, NM)
   CASE('GSTRIX') 
      CALL SCARC_VECCOPY3D (SL%R, SL%G, ILMAX, NM)
      CALL SCARC_PRECON_GSTRIX3D (SL%G, ILMAX, NM)
      CALL SCARC_SCALPROD3D (SL%R, SL%G, SIGMA0, ILMAX, NM)
END SELECT
CALL SCARC_VECADD3D (SL%G, SL%D,-1.0_EB, 0.0_EB, ILMAX, NM)
 
!
! start defect correction loop
!
CG_LOOP3D: DO ITE = 1, SCARC_NIT_CG
 
   ! calculate new defect and get L2-norm of it
   CALL SCARC_MATVEC3D (SL%D, SL%Y, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILMAX, 5)
   CALL SCARC_SCALPROD3D (SL%D, SL%Y, ALPHA, ILMAX, NM)

   ! get new descent direction
   ALPHA = SIGMA0 / ALPHA
   CALL SCARC_VECADD3D (SL%D, SL%X, ALPHA, 1.0_EB, ILMAX, NM)
   CALL SCARC_VECADD3D (SL%Y, SL%R, ALPHA, 1.0_EB, ILMAX, NM)
   CALL SCARC_L2NORM3D (SL%R, SL%RES_CG, ILMAX, NM, NTYPE_GLOBAL)

   !  exit loop in case of convergence or divergence
   IF (SCARC_DEBUG >= 1) write (9,1000) ITE, ILMAX, SL%RES_CG
   IF (SCARC_DEBUG >= 1) write (*,1000) ITE, ILMAX, SL%RES_CG
   !IF (NM == 1)          write (*,1000) ITE, ILMAX, SL%RES_CG
   IF (SCARC_BREL) THEN
      IF (SL%RES_CG <= SL%RESIN_CG*SCARC_EPS_CG) BCONV = .TRUE.
   ELSE
      IF (SL%RES_CG <= SCARC_EPS_CG .AND. SL%RES_CG <= SL%RESIN_CG*SCARC_EPS_REL) BCONV = .TRUE.
   ENDIF
   IF (SL%RES_CG > SCARC_EPS_DIVG) BDIVG = .TRUE.
   IF (BCONV.OR.BDIVG) EXIT CG_LOOP3D
 
   ! preconditioning
   SELECT CASE(SCARC_PRECON_CG)
      CASE('NONE')
         SIGMA1 = SL%RES_CG * SL%RES_CG
      CASE('JACOBI')
         CALL SCARC_VECCOPY3D (SL%R, SL%G, ILMAX, NM)
         CALL SCARC_PRECON_JACOBI3D (SL%G, ILMAX, NM)
         CALL SCARC_SCALPROD3D (SL%R, SL%G, SIGMA1, ILMAX, NM)
      CASE('SSOR')
         CALL SCARC_VECCOPY3D (SL%R, SL%G, ILMAX, NM)
         CALL SCARC_PRECON_SSOR3D (SL%G, ILMAX, NM)
         CALL SCARC_SCALPROD3D (SL%R, SL%G, SIGMA1, ILMAX, NM)
      CASE('GSTRIX')
         CALL SCARC_VECCOPY3D (SL%R, SL%G, ILMAX, NM)
         CALL SCARC_PRECON_GSTRIX3D (SL%G, ILMAX, NM)
         CALL SCARC_SCALPROD3D (SL%R, SL%G, SIGMA1, ILMAX, NM)
   END SELECT
 
   GAMMA = SIGMA1 / SIGMA0
   SIGMA0 = SIGMA1
   CALL SCARC_VECADD3D (SL%G, SL%D,-1.0_EB, GAMMA, ILMAX, NM)
 
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
IF (SCARC_DEBUG >= 1) WRITE (9,2000) ITE0, ILMAX, ILMAX, SL%RES_CG, SL%CAPPA_CG 
IF (NM == 1)          WRITE (*,2000) ITE0, ILMAX, ILMAX, SL%RES_CG, SL%CAPPA_CG

! shift solution back to vector PRHS
DO K = 1, KBAR
   DO I = 1, IBAR
      M%PRHS (I, 1, K) = SL%X(I, 1, K)
   ENDDO
ENDDO
 
TUSED_SCARC(6,NM)=TUSED_SCARC(6,NM)+SECOND()-TNOW_CG3D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_CG3D

1000 FORMAT ('=====SCARC_CG       : #ite= ',i4,': level=',i4,': res=',e12.4)
2000 FORMAT ('=====SCARC_CG       : #ite= ',i4,': level=',i4,': res=',e12.4,': rate=',e12.4)
END SUBROUTINE SCARC_CG3D

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Perform global MG-method based on global possion-matrix
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_MG2D (NM)
 
INTEGER :: NM, ITE, ITE0, ICYCLE
INTEGER :: ILEVEL, ILMAX
INTEGER :: I, K, IBAR0, JBAR0, KBAR0
REAL (EB) :: TNOW_MG2D
LOGICAL BMATVEC, BCONV, BDIVG
!TYPE (SCARC_LEVEL_TYPE), POINTER :: SL, SLHI, SLLO
TYPE (MESH_TYPE),  POINTER :: M
 
TNOW_MG2D = SECOND()

! initialize working vectors
S    => SCARC (NM)
M    => MESHES (NM)

ILMAX =  S%NLEVEL_MAX
SLHI  => S%SLEVEL(ILMAX)
 
! initialize working vectors
SLHI%X = 0.0_EB               ! necessary ? better without ??
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
   WRITE (9,*) 'SCARC_NIT_MG=:', SCARC_NIT_MG
   WRITE (9,*) 'SCARC_EPS_MG=:', SCARC_EPS_MG
   WRITE (9,*) 'SCARC_EPS_REL=:', SCARC_EPS_REL
   CALL SCARC_COMPARE_SINGLE (SLHI%X, 'SARC',  'X1    ')
   CALL SCARC_COMPARE_SINGLE (SLHI%F, 'SARC',  'F1    ')
ENDIF
 
!!!
!!! initialize some method parameters
!!!
 
NREQ_FACE = 0
NREQ_EDGE = 0
NREQ_DIAG = 0

BCONV =.FALSE.
BDIVG =.FALSE.
ICYCLE=1      ! V-cycle
 
IBAR0=SLHI%IBAR
JBAR0=SLHI%JBAR
KBAR0=SLHI%KBAR

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
   S%KCYCLE(2,S%NLEVEL_MAX)=1
   DO ILEVEL=S%NLEVEL_MIN+1,S%NLEVEL_MAX-1
      IF (ICYCLE==0) THEN
         S%KCYCLE(2,ILEVEL)=2
      ELSE
         S%KCYCLE(2,ILEVEL)=ICYCLE
      ENDIF
   ENDDO
     
   !!! calculate initial defect on finest level, Y = B - A*X, and get l2-norm of it
   ILEVEL=ILMAX
   CALL SCARC_MATVEC2D (SLHI%X, SLHI%D, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, 3)       ! 3 korrekt ???
   CALL SCARC_VECADD2D (SLHI%F , SLHI%D, 1.0_EB, -1.0_EB, ILEVEL, NM)
   CALL SCARC_L2NORM2D (SLHI%D , SLHI%RESIN_MG, ILEVEL, NM, NTYPE_GLOBAL)
   IF (SCARC_DEBUG >= 1) write (9,1000) 0, ILMAX, SL%RESIN_MG
   IF (SCARC_DEBUG >= 1) write (*,1000) 0, ILMAX, SL%RESIN_MG

   !CALL SCARC_MASTER_INFO (SLHI%ITE, SLHI%RESIN_MG)

    
   !!! ---------------------------------------------------------------------------
   !!! start MG-iteration
   !!! ---------------------------------------------------------------------------
   MG_LOOP2D: DO ITE = 1, SCARC_NIT_MG
    
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

         CALL SCARC_SMOOTHER2D(SLHI%X, SLHI%F, SLHI%D, BMATVEC, ILEVEL, NM)

         !!! print level information - still missing
         !CALL SCARC_INFO_LEVEL(ILEVEL)

         !!! perform restriction to coarser grid: F:= rest(D)
         CALL SCARC_COMPARE_SINGLE_LEV (SLHI%D, 'SARC',  'RESTD ',ILEVEL)
         CALL SCARC_RESTRICTION2D(SLLO%F, SLHI%D, ILEVEL-1, NM)
         CALL SCARC_COMPARE_SINGLE_LEV (SLLO%F, 'SARC',  'RESTF ',ILEVEL-1)

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
      ILEVEL=S%NLEVEL_MIN
      CALL SCARC_COARSE2D(1, ILEVEL, NM)


      !!! ------------------------------------------------------------------------
      !!! Postsmoothing
      !!! ------------------------------------------------------------------------
120   IF (ILEVEL/=S%NLEVEL_MAX) THEN
    
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
         CALL SCARC_VECADD2D(SLHI%D, SLHI%X, 1.0_EB, 1.0_EB, ILEVEL, NM)

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
      CALL SCARC_MATVEC2D (SLHI%X, SLHI%D, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, 3)
      CALL SCARC_VECADD2D (SLHI%F, SLHI%D,1.0_EB, -1.0_EB, ILEVEL, NM)
      CALL SCARC_L2NORM2D (SLHI%D, SLHI%RES_MG, ILEVEL, NM, NTYPE_GLOBAL)

      IF (SCARC_DEBUG>=2) WRITE(9,*) 'SCARC-Multigrid, iteration ',ITE,': residuum=',SLHI%RES_MG

      !!! ------------------------------------------------------------------------
      !!! Send error information to Master  - still missing
      !!! ------------------------------------------------------------------------
      !CALL SCARC_INFO_MASTER 
      !CALL SCARC_INFO_LEVEL(ILEVEL)

      !!! ------------------------------------------------------------------------
      !!! Convergence or divergence ?
      !!! ------------------------------------------------------------------------
      IF (SCARC_DEBUG >= 1) write (9,1000) ITE, ILEVEL, SL%RES_MG
      IF (SCARC_DEBUG >= 1) write (*,1000) ITE, ILEVEL, SL%RES_MG
      IF (SCARC_BREL) THEN
         IF (SL%RES_MG <= SL%RESIN_MG*SCARC_EPS_MG) BCONV = .TRUE.
      ELSE
         IF (SL%RES_MG <= SCARC_EPS_MG .AND. SL%RES_MG <= SL%RESIN_MG*SCARC_EPS_REL) BCONV = .TRUE.
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
   IF (SCARC_DEBUG >= 1) WRITE (9,2000) ITE0, ILEVEL, SL%RES_MG, SL%CAPPA_MG 
   IF (SCARC_DEBUG >= 1) WRITE (*,2000) ITE0, ILEVEL, SL%RES_MG, SL%CAPPA_MG 

ENDIF ONLY_ONE_LEVEL_IF
 
DO K = 1, SLHI%KBAR
   DO I = 1, SLHI%IBAR
      M%PRHS (I, 1, K) = SLHI%X(I, 1, K)
   ENDDO
ENDDO
 
TUSED_SCARC(7,NM)=TUSED_SCARC(7,NM)+SECOND()-TNOW_MG2D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_MG2D

1000 FORMAT ('=====SCARC_MG       : #ite= ',i4,': level=',i4,': res=',e12.4)
2000 FORMAT ('=====SCARC_MG       : #ite= ',i4,': level=',i4,': res=',e12.4,': rate=',e12.4)
END SUBROUTINE SCARC_MG2D



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Perform global MG-method based on global possion-matrix
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_MG3D (NM)
 
INTEGER :: NM, ITE, ITE0, ICYCLE
INTEGER :: ILEVEL, ILMAX
INTEGER :: I, K, IBAR0, JBAR0, KBAR0
REAL (EB) :: TNOW_MG3D
LOGICAL BMATVEC, BCONV, BDIVG
!TYPE (SCARC_LEVEL_TYPE), POINTER :: SL, SLHI, SLLO
TYPE (MESH_TYPE),  POINTER :: M
 
TNOW_MG3D = SECOND()

! initialize working vectors
S    => SCARC (NM)
M    => MESHES (NM)

ILMAX =  S%NLEVEL_MAX
SLHI  => S%SLEVEL(ILMAX)
 
! initialize working vectors
SLHI%X = 0.0_EB
DO K = 1, KBAR
   DO I = 1, IBAR
      SLHI%F (I, 1, K) = M%PRHS (I, 1, K)
   ENDDO
ENDDO
 
IF (SCARC_DEBUG >= 3) THEN
   WRITE (9,*) 'Starting SCARC_MG3D'
   WRITE (9,*) 'IBAR=:', SLHI%IBAR
   WRITE (9,*) 'JBAR=:', SLHI%JBAR
   WRITE (9,*) 'KBAR=:', SLHI%KBAR
   WRITE (9,*) 'SCARC_NIT_MG=:', SCARC_NIT_MG
   WRITE (9,*) 'SCARC_EPS_MG=:', SCARC_EPS_MG
   WRITE (9,*) 'SCARC_EPS_REL=:', SCARC_EPS_REL
   CALL SCARC_COMPARE_SINGLE (SLHI%X, 'SARC',  'X1    ')
   CALL SCARC_COMPARE_SINGLE (SLHI%F, 'SARC',  'F1    ')
ENDIF
 
!!!
!!! initialize some method parameters
!!!
NREQ_FACE = 0
NREQ_EDGE = 0
NREQ_DIAG = 0

BCONV =.FALSE.
BDIVG =.FALSE.
ICYCLE=1      ! V-cycle
 
IBAR0=SLHI%IBAR
JBAR0=SLHI%JBAR
KBAR0=SLHI%KBAR

!!! adjust boundary values of right hand side
CALL SCARC_SETBDRY3D (ILMAX, NM)
 
!!! ------------------------------------------------------------------------------------
!!! Only one level: solve problem exactly
!!! ------------------------------------------------------------------------------------
ONLY_ONE_LEVEL_IF: IF (S%NLEVEL==1) THEN   

      CALL SCARC_COARSE3D(1,ILMAX, NM)

!!! ------------------------------------------------------------------------------------
!!! More than one level: start MG-cycling
!!! ------------------------------------------------------------------------------------
ELSE

   !!! save cycle counts for MG-iteration
   S%KCYCLE(2,S%NLEVEL_MAX)=1
   DO ILEVEL=S%NLEVEL_MIN+1,S%NLEVEL_MAX-1
      IF (ICYCLE==0) THEN
         S%KCYCLE(2,ILEVEL)=2
      ELSE
         S%KCYCLE(2,ILEVEL)=ICYCLE
      ENDIF
   ENDDO
     
   !!! calculate initial defect on finest level, Y = B - A*X, and get l2-norm of it
!WRITE(9,*) 'KOMMUNIKATIONSPARAMETER 3 prüfen, anders als bei CG !!!'
   ILEVEL=ILMAX
   CALL SCARC_MATVEC3D (SLHI%X, SLHI%D, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, 3)
   CALL SCARC_VECADD3D (SLHI%F , SLHI%D, 1.0_EB, -1.0_EB, ILEVEL, NM)
   CALL SCARC_L2NORM3D (SLHI%D , SLHI%RES_MG, ILEVEL, NM, NTYPE_GLOBAL)

   ITE   = 0
   SLHI%RESIN_MG = SLHI%RES_MG
    
   !WRITE(*,*) 'Info muss noch an Master geschickt werden'
   !CALL SCARC_MASTER_INFO (SLHI%ITE, SLHI%RESIN_MG)

    
   !!! ---------------------------------------------------------------------------
   !!! start MG-iteration
   !!! ---------------------------------------------------------------------------
   MG_LOOP3D: DO ITE = 1, SCARC_NIT_MG
    
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

         CALL SCARC_SMOOTHER3D(SLHI%X, SLHI%F, SLHI%D, BMATVEC, ILEVEL, NM)

         !!! print level information - still missing
         !CALL SCARC_INFO_LEVEL(ILEVEL)

         !!! perform restriction to coarser grid: F:= rest(D)
         CALL SCARC_COMPARE_SINGLE_LEV (SLHI%D, 'SARC',  'RESTD ',ILEVEL)
         CALL SCARC_RESTRICTION3D(SLLO%F, SLHI%D, ILEVEL-1, NM)
         CALL SCARC_COMPARE_SINGLE_LEV (SLLO%F, 'SARC',  'RESTF ',ILEVEL-1)

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
      ILEVEL=S%NLEVEL_MIN
      CALL SCARC_COARSE3D(1,ILEVEL, NM)


      !!! ------------------------------------------------------------------------
      !!! Postsmoothing
      !!! ------------------------------------------------------------------------
120   IF (ILEVEL/=S%NLEVEL_MAX) THEN
    
         BMATVEC=.TRUE.

         SLLO => S%SLEVEL(ILEVEL)

         ILEVEL=ILEVEL+1
         SLHI => S%SLEVEL(ILEVEL)

         !!! print level information - still missing
         !CALL SCARC_INFO_LEVEL(ILEVEL)

         !!! perform prolongation to finer grid
         CALL SCARC_PROLONGATION3D(SLLO%X, SLHI%D, ILEVEL, NM)

         !!! set exterior boundary data of residuum to zero - necessary ?
         !CALL SCARC_BDRY_RESIDUUM(NM, SLHI%D, ILEVEL)
 
         !!! set new solution
         CALL SCARC_VECADD3D(SLHI%D, SLHI%X, 1.0_EB, 1.0_EB, ILEVEL, NM)

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
      CALL SCARC_MATVEC3D (SLHI%X, SLHI%D, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, 3)
      CALL SCARC_VECADD3D (SLHI%F , SLHI%D,1.0_EB, -1.0_EB, ILEVEL, NM)
      CALL SCARC_L2NORM3D (SLHI%D , SLHI%RES_MG, ILEVEL, NM, NTYPE_GLOBAL)

      IF (SCARC_DEBUG>=2) WRITE(9,*) 'SCARC-Multigrid, iteration ',ITE,': residuum=',SLHI%RES_MG

      !!! ------------------------------------------------------------------------
      !!! Send error information to Master  - still missing
      !!! ------------------------------------------------------------------------
      !CALL SCARC_INFO_MASTER 
      !CALL SCARC_INFO_LEVEL(ILEVEL)

      !!! ------------------------------------------------------------------------
      !!! Convergence or divergence ?
      !!! ------------------------------------------------------------------------
      IF (SCARC_DEBUG >= 2) WRITE (9,1000) ITE, ILEVEL, SL%RES_MG
      !IF (NM == 1)          WRITE (9,1000) ITE, ILEVEL, SL%RES_MG
      IF (SCARC_BREL) THEN
         IF (SL%RES_MG <= SL%RESIN_MG*SCARC_EPS_MG) BCONV = .TRUE.
      ELSE
         IF (SL%RES_MG <= SCARC_EPS_MG .AND. SL%RES_MG <= SL%RESIN_MG*SCARC_EPS_REL) BCONV = .TRUE.
      ENDIF
      IF (SL%RES_MG > SCARC_EPS_DIVG) BDIVG = .TRUE.
      IF (BCONV.OR.BDIVG) EXIT MG_LOOP3D
 
   ENDDO MG_LOOP3D

   !!! ------------------------------------------------------------------------
   !!! Convergence or divergence ?
   !!! ------------------------------------------------------------------------
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
   IF (SCARC_DEBUG >= 1) WRITE (9,2000) ITE0, ILEVEL, SL%RES_CG, SL%CAPPA_CG 
   IF (NM == 1)          WRITE (*,2000) ITE0, ILEVEL, SL%RES_CG, SL%CAPPA_CG

ENDIF ONLY_ONE_LEVEL_IF
 
DO K = 1, SLHI%KBAR
   DO I = 1, SLHI%IBAR
      M%PRHS (I, 1, K) = SLHI%X(I, 1, K)
   ENDDO
ENDDO
 
TUSED_SCARC(7,NM)=TUSED_SCARC(7,NM)+SECOND()-TNOW_MG3D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_MG3D

1000 FORMAT ('=====SCARC_MG       : #ite= ',i4,': level=',i4,': res=',e12.4)
2000 FORMAT ('=====SCARC_MG       : #ite= ',i4,': level=',i4,': res=',e12.4,': rate=',e12.4)
END SUBROUTINE SCARC_MG3D




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
 
!!! calculate initial defect Y = B - A*X and get l2-norm of it
CALL SCARC_MATVEC2D (SL%X, SL%R, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, 3)
CALL SCARC_VECADD2D (SL%F, SL%R,-1.0_EB, 1.0_EB, ILEVEL, NM)
CALL SCARC_L2NORM2D (SL%R, SL%RESIN_CO, ILEVEL, NM, NTYPE_GLOBAL)
IF (SCARC_DEBUG >= 1) write (9,1000) 0, ILEVEL, SL%RESIN_CO
IF (SCARC_DEBUG >= 1) write (*,1000) 0, ILEVEL, SL%RESIN_CO
 
! initial preconditioning
SELECT CASE(SCARC_PRECON_CO)
   CASE('NONE') 
      SIGMA0 = SL%RESIN_CO * SL%RESIN_CO
   CASE('JACOBI') 
      CALL SCARC_VECCOPY2D (SL%R, SL%G, ILEVEL, NM)
      CALL SCARC_PRECON_JACOBI2D (SL%G, ILEVEL, NM)
      CALL SCARC_SCALPROD2D (SL%R, SL%G, SIGMA0, ILEVEL, NM)
   CASE('SSOR') 
      CALL SCARC_VECCOPY2D (SL%R, SL%G, ILEVEL, NM)
      CALL SCARC_PRECON_SSOR2D (SL%G, ILEVEL, NM)
      CALL SCARC_SCALPROD2D (SL%R, SL%G, SIGMA0, ILEVEL, NM)
   CASE('GSTRIX') 
      CALL SCARC_VECCOPY2D (SL%R, SL%G, ILEVEL, NM)
      CALL SCARC_PRECON_GSTRIX2D (SL%G, ILEVEL, NM)
      CALL SCARC_SCALPROD2D (SL%R, SL%G, SIGMA0, ILEVEL, NM)
END SELECT
CALL SCARC_VECADD2D (SL%G, SL%D,-1.0_EB, 0.0_EB, ILEVEL, NM)
 
!
! start defect correction loop
!
CO_LOOP2D: DO ITE = 1, SCARC_NIT_CO
 
   ! calculate new defect and get L2-norm of it
   CALL SCARC_MATVEC2D (SL%D, SL%Y, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, 5)
   CALL SCARC_SCALPROD2D (SL%D, SL%Y, ALPHA, ILEVEL, NM)

   ! get new descent direction
   ALPHA = SIGMA0 / ALPHA
   CALL SCARC_VECADD2D (SL%D, SL%X, ALPHA, 1.0_EB, ILEVEL, NM)
   CALL SCARC_VECADD2D (SL%Y, SL%R, ALPHA, 1.0_EB, ILEVEL, NM)
   CALL SCARC_L2NORM2D (SL%R, SL%RES_CO, ILEVEL, NM, NTYPE_GLOBAL)

 
   !  exit loop in case of convergence or divergence
   IF (SCARC_DEBUG >= 1) write (9,1000) ITE, ILEVEL, SL%RES_CO
   IF (SCARC_DEBUG >= 1) write (*,1000) ITE, ILEVEL, SL%RES_CO
   IF (SCARC_BREL) THEN
      IF (SL%RES_CO <= SL%RESIN_CO*SCARC_EPS_CO) BCONV = .TRUE.
   ELSE
      IF (SL%RES_CO <= SCARC_EPS_CO .AND. SL%RES_CO <= SL%RESIN_CO*SCARC_EPS_REL) BCONV = .TRUE.
   ENDIF
   IF (SL%RES_CO > SCARC_EPS_DIVG) BDIVG = .TRUE.
   IF (BCONV.OR.BDIVG) EXIT CO_LOOP2D
 
   ! preconditioning
   SELECT CASE(SCARC_PRECON_CO)
      CASE('NONE')
         SIGMA1 = SL%RES_CO * SL%RES_CO
      CASE('JACOBI')
         CALL SCARC_VECCOPY2D (SL%R, SL%G, ILEVEL, NM)
         CALL SCARC_PRECON_JACOBI2D (SL%G, ILEVEL, NM)
         CALL SCARC_SCALPROD2D (SL%R, SL%G, SIGMA1, ILEVEL, NM)
      CASE('SSOR')
         CALL SCARC_VECCOPY2D (SL%R, SL%G, ILEVEL, NM)
         CALL SCARC_PRECON_SSOR2D (SL%G, ILEVEL, NM)
         CALL SCARC_SCALPROD2D (SL%R, SL%G, SIGMA1, ILEVEL, NM)
      CASE('GSTRIX')
         CALL SCARC_VECCOPY2D (SL%R, SL%G, ILEVEL, NM)
         CALL SCARC_PRECON_GSTRIX2D (SL%G, ILEVEL, NM)
         CALL SCARC_SCALPROD2D (SL%R, SL%G, SIGMA1, ILEVEL, NM)
   END SELECT
 
   GAMMA = SIGMA1 / SIGMA0
   SIGMA0 = SIGMA1
   CALL SCARC_VECADD2D (SL%G, SL%D,-1.0_EB, GAMMA, ILEVEL, NM)
 
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
IF (SCARC_DEBUG >= 1) WRITE (9,2000) ITE0, ILEVEL, SL%RES_CO, SL%CAPPA_CO 
IF (SCARC_DEBUG >= 1) WRITE (*,2000) ITE0, ILEVEL, SL%RES_CO, SL%CAPPA_CO 

 
TUSED_SCARC(9,NM)=TUSED_SCARC(9,NM)+SECOND()-TNOW_COARSE
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_COARSE

1000 FORMAT ('     SCARC_COARSE   : #ite= ',i4,': level=',i4,': res=',e12.4)
2000 FORMAT ('     SCARC_COARSE   : #ite= ',i4,': level=',i4,': res=',e12.4,': rate=',e12.4)
END SUBROUTINE SCARC_COARSE2D



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Perform global 3D cg-method based on global possion-matrix
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_COARSE3D (LAYER, ILEVEL, NM)
 
INTEGER ::  NM, LAYER, ILEVEL, ITE, ITE0
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
IF (LAYER==2) SL%F = SL%F2
SL%X = 0.0_EB                   
SL%G = 0.0_EB                   
SL%Y = 0.0_EB
SL%R = 0.0_EB
SL%D = 0.0_EB
 
!!! calculate initial defect Y = B - A*X and get l2-norm of it
CALL SCARC_MATVEC3D (SL%X, SL%R, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, 3)
CALL SCARC_VECADD3D (SL%F, SL%R,-1.0_EB, 1.0_EB, ILEVEL, NM)
CALL SCARC_L2NORM3D (SL%R, SL%RESIN_CO, ILEVEL, NM, NTYPE_GLOBAL)
IF (SCARC_DEBUG >= 1) write (9,1000) 0, ILEVEL, SL%RESIN_CO
IF (SCARC_DEBUG >= 1) write (*,1000) 0, ILEVEL, SL%RESIN_CO
 
! initial preconditioning
SELECT CASE(SCARC_PRECON_CO)
   CASE('NONE') 
      SIGMA0 = SL%RESIN_CO * SL%RESIN_CO
   CASE('JACOBI') 
      CALL SCARC_VECCOPY3D (SL%R, SL%G, ILEVEL, NM)
      CALL SCARC_PRECON_JACOBI3D (SL%G, ILEVEL, NM)
      CALL SCARC_SCALPROD3D (SL%R, SL%G, SIGMA0, ILEVEL, NM)
   CASE('SSOR') 
      CALL SCARC_VECCOPY3D (SL%R, SL%G, ILEVEL, NM)
      CALL SCARC_PRECON_SSOR3D (SL%G, ILEVEL, NM)
      CALL SCARC_SCALPROD3D (SL%R, SL%G, SIGMA0, ILEVEL, NM)
   CASE('GSTRIX') 
      CALL SCARC_VECCOPY3D (SL%R, SL%G, ILEVEL, NM)
      CALL SCARC_PRECON_GSTRIX3D (SL%G, ILEVEL, NM)
      CALL SCARC_SCALPROD3D (SL%R, SL%G, SIGMA0, ILEVEL, NM)
END SELECT
CALL SCARC_VECADD3D (SL%G, SL%D,-1.0_EB, 0.0_EB, ILEVEL, NM)
 
!
! start defect correction loop
!
CO_LOOP3D: DO ITE = 1, SCARC_NIT_CO
 
   ! calculate new defect and get L2-norm of it
   CALL SCARC_MATVEC3D (SL%D, SL%Y, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, 5)
   CALL SCARC_SCALPROD3D (SL%D, SL%Y, ALPHA, ILEVEL, NM)

   ! get new descent direction
   ALPHA = SIGMA0 / ALPHA
   CALL SCARC_VECADD3D (SL%D, SL%X, ALPHA, 1.0_EB, ILEVEL, NM)
   CALL SCARC_VECADD3D (SL%Y, SL%R, ALPHA, 1.0_EB, ILEVEL, NM)
   CALL SCARC_L2NORM3D (SL%R, SL%RES_CO, ILEVEL, NM, NTYPE_GLOBAL)

 
   !  exit loop in case of convergence or divergence
   IF (SCARC_DEBUG >= 1) write (9,1000) ITE, ILEVEL, SL%RES_CO
   IF (SCARC_DEBUG >= 1) write (*,1000) ITE, ILEVEL, SL%RES_CO
   IF (SCARC_BREL) THEN
      IF (SL%RES_CO <= SL%RESIN_CO*SCARC_EPS_CO) BCONV = .TRUE.
   ELSE
      IF (SL%RES_CO <= SCARC_EPS_CO .AND. SL%RES_CO <= SL%RESIN_CO*SCARC_EPS_REL) BCONV = .TRUE.
   ENDIF
   IF (SL%RES_CO > SCARC_EPS_DIVG) BDIVG = .TRUE.
   IF (BCONV.OR.BDIVG) EXIT CO_LOOP3D
 
   ! preconditioning
   SELECT CASE(SCARC_PRECON_CO)
      CASE('NONE')
         SIGMA1 = SL%RES_CO * SL%RES_CO
      CASE('JACOBI')
         CALL SCARC_VECCOPY3D (SL%R, SL%G, ILEVEL, NM)
         CALL SCARC_PRECON_JACOBI3D (SL%G, ILEVEL, NM)
         CALL SCARC_SCALPROD3D (SL%R, SL%G, SIGMA1, ILEVEL, NM)
      CASE('SSOR')
         CALL SCARC_VECCOPY3D (SL%R, SL%G, ILEVEL, NM)
         CALL SCARC_PRECON_SSOR3D (SL%G, ILEVEL, NM)
         CALL SCARC_SCALPROD3D (SL%R, SL%G, SIGMA1, ILEVEL, NM)
      CASE('GSTRIX')
         CALL SCARC_VECCOPY3D (SL%R, SL%G, ILEVEL, NM)
         CALL SCARC_PRECON_GSTRIX3D (SL%G, ILEVEL, NM)
         CALL SCARC_SCALPROD3D (SL%R, SL%G, SIGMA1, ILEVEL, NM)
   END SELECT
 
   GAMMA = SIGMA1 / SIGMA0
   SIGMA0 = SIGMA1
   CALL SCARC_VECADD3D (SL%G, SL%D,-1.0_EB, GAMMA, ILEVEL, NM)
 
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
IF (SCARC_DEBUG >= 1) WRITE (9,2000) ITE0, ILEVEL, SL%RES_CO, SL%CAPPA_CO 
IF (NM == 1)          WRITE (*,2000) ITE0, ILEVEL, SL%RES_CO, SL%CAPPA_CO
 
TUSED_SCARC(9,NM)=TUSED_SCARC(9,NM)+SECOND()-TNOW_COARSE
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_COARSE

1000 FORMAT ('     SCARC_COARSE   : #ite= ',i4,': level=',i4,': res=',e12.4)
2000 FORMAT ('     SCARC_COARSE   : #ite= ',i4,': level=',i4,': res=',e12.4,': rate=',e12.4)
END SUBROUTINE SCARC_COARSE3D


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

CALL SCARC_COMPARE_SINGLE_LEV (X, 'SARC',  'SMOX  ',ILEVEL)
CALL SCARC_COMPARE_SINGLE_LEV (F, 'SARC',  'SMOF  ',ILEVEL)

SL%RESIN_SM=1.0_EB
IF (BMATVEC) THEN
   CALL SCARC_MATVEC2D (X, D, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, 3)   ! global MATVEC
   CALL SCARC_VECADD2D (F, D, 1.0_EB, -1.0_EB, ILEVEL, NM)
   CALL SCARC_COMPARE_SINGLE_LEV (D, 'SARC',  'SMOD  ',ILEVEL)
   IF (BL2NORM) CALL SCARC_L2NORM2D (D, SL%RESIN_SM, ILEVEL, NM, NTYPE_GLOBAL)
ENDIF

IF (SCARC_DEBUG >= 1) WRITE(9,1000) 0, ILEVEL, SL%RESIN_SM
IF (SCARC_DEBUG >= 1) WRITE(*,1000) 0, ILEVEL, SL%RESIN_SM

SMOOTH_LOOP2D: DO ITE=1,SCARC_NIT_SM
 
   SELECT CASE(SCARC_PRECON_SM)
      CASE('JACOBI') 
         CALL SCARC_PRECON_JACOBI2D (D, ILEVEL, NM)
      CASE('SSOR') 
         CALL SCARC_PRECON_SSOR2D (D, ILEVEL, NM)
      CASE('GSTRIX') 
         CALL SCARC_PRECON_GSTRIX2D (D, ILEVEL, NM)
   END SELECT

   CALL SCARC_VECADD2D(D, X, SCARC_OMEGA_SM, 1.0_EB, ILEVEL, NM)
   CALL SCARC_MATVEC2D (X, D, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, 3)   ! global MATVEC
   CALL SCARC_VECADD2D (F, D, 1.0_EB, -1.0_EB, ILEVEL, NM)
   IF (BL2NORM) CALL SCARC_L2NORM2D (D, SL%RES_SM, ILEVEL, NM, NTYPE_GLOBAL)

   IF (SCARC_DEBUG >= 1) WRITE(9,1000) ITE, ILEVEL, SL%RES_SM
   IF (SCARC_DEBUG >= 1) WRITE(*,1000) ITE, ILEVEL, SL%RES_SM
   IF (SCARC_BREL) THEN
      IF (SL%RES_SM <= SL%RESIN_SM*SCARC_EPS_SM) BCONV = .TRUE.
   ELSE
      IF (SL%RES_SM <= SCARC_EPS_SM .AND. SL%RES_SM <= SL%RESIN_SM*SCARC_EPS_REL) BCONV = .TRUE.
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
IF (SCARC_DEBUG >= 1) WRITE (9,2000) ITE0, ILEVEL, SL%RES_SM, SL%CAPPA_SM 
IF (SCARC_DEBUG >= 1) WRITE (*,2000) ITE0, ILEVEL, SL%RES_SM, SL%CAPPA_SM 

TUSED_SCARC(8,NM)=TUSED_SCARC(8,NM)+SECOND()-TNOW_SMOOTHER2D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_SMOOTHER2D

1000 FORMAT ('     SCARC_SMOOTHER : #ite= ',i4,': level=',i4,': res=',e12.4)
2000 FORMAT ('     SCARC_SMOOTHER : #ite= ',i4,': level=',i4,': res=',e12.4,': rate=',e12.4)
END SUBROUTINE SCARC_SMOOTHER2D


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

CALL SCARC_COMPARE_SINGLE_LEV (X, 'SARC', 'SMOX  ',ILEVEL)
CALL SCARC_COMPARE_SINGLE_LEV (F, 'SARC', 'SMOF  ',ILEVEL)

SL%RESIN_SM=0.0_EB
IF (BMATVEC) THEN
   CALL SCARC_MATVEC3D (X, D, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, 3)   ! global MATVEC
   CALL SCARC_VECADD3D (F, D, 1.0_EB, -1.0_EB, ILEVEL, NM)
   CALL SCARC_COMPARE_SINGLE_LEV (D, 'SARC',  'SMOD  ',ILEVEL)
   IF (BL2NORM) CALL SCARC_L2NORM3D (D, SL%RESIN_SM, ILEVEL, NM, NTYPE_GLOBAL)
ENDIF

IF (SCARC_DEBUG >= 1) WRITE(9,1000) 0, ILEVEL, SL%RESIN_SM
IF (SCARC_DEBUG >= 1) WRITE(*,1000) 0, ILEVEL, SL%RESIN_SM
!IF (SL%RESIN_SM>SCARC_EPS_DIVG) GOTO 210 

SMOOTH_LOOP3D: DO ITE=1,SCARC_NIT_SM
 
   SELECT CASE(SCARC_PRECON_SM)
      CASE('JACOBI') 
         CALL SCARC_PRECON_JACOBI3D (D, ILEVEL, NM)
      CASE('SSOR') 
         CALL SCARC_PRECON_SSOR3D (D, ILEVEL, NM)
      CASE('GSTRIX') 
         CALL SCARC_PRECON_GSTRIX3D (D, ILEVEL, NM)
   END SELECT

   CALL SCARC_VECADD3D(D, X, SCARC_OMEGA_SM, 1.0_EB, ILEVEL, NM)
   CALL SCARC_MATVEC3D (X, D, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, 3)   ! global MATVEC
   CALL SCARC_VECADD3D (F, D, 1.0_EB, -1.0_EB, ILEVEL, NM)
   IF (BL2NORM) CALL SCARC_L2NORM3D (D, SL%RES_SM, ILEVEL, NM, NTYPE_GLOBAL)

   IF (SCARC_DEBUG >= 1) WRITE(9,1000) ITE, ILEVEL, SL%RES_SM
   IF (SCARC_DEBUG >= 1) WRITE(*,1000) ITE, ILEVEL, SL%RES_SM
   IF (SCARC_BREL) THEN
      IF (SL%RES_SM <= SL%RESIN_SM*SCARC_EPS_SM) BCONV = .TRUE.
   ELSE
      IF (SL%RES_SM <= SCARC_EPS_SM .AND. SL%RES_SM <= SL%RESIN_SM*SCARC_EPS_REL) BCONV = .TRUE.
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
IF (SCARC_DEBUG >= 1) WRITE (9,2000) ITE0, ILEVEL, SL%RES_SM, SL%CAPPA_SM 
IF (NM == 1)          WRITE (*,2000) ITE0, ILEVEL, SL%RES_SM, SL%CAPPA_SM

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
SUBROUTINE SCARC_MATVEC2D (X, Y, A1, A2, NM, ITYPE, ILEVEL, IMV)
 
REAL (EB) :: A1, A2
REAL (EB), POINTER, DIMENSION (:, :, :) :: X, Y
REAL (EB):: TNOW_MATVEC2D
INTEGER :: NM, ITYPE
INTEGER :: I, K, IC, IMV, ILEVEL
 
TNOW_MATVEC2D = SECOND()
 
SL => S%SLEVEL(ILEVEL)
 
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
   ENDDO
ENDDO
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform data exchange at interior boundaries
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF ((ITYPE == NTYPE_GLOBAL) .AND. (NMESHES > 1)) THEN
   NREQ_FACE = 0
   NREQ_EDGE = 0
   NREQ_DIAG = 0
   CALL SCARC_RECEIVE (NCOM_MATV,  ILEVEL)   
   CALL SCARC_EXCHANGE (NCOM_MATV,  ILEVEL, IMV)
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
SUBROUTINE SCARC_MATVEC3D (X, Y, A1, A2, NM, ITYPE, ILEVEL, IMV)
 
REAL (EB) :: A1, A2
REAL (EB), POINTER, DIMENSION (:, :, :) :: X, Y
REAL (EB):: TNOW_MATVEC3D
INTEGER :: NM, ITYPE
INTEGER :: I, J, K, IC, IMV, ILEVEL
 
 
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
   NREQ_EDGE = 0
   NREQ_DIAG = 0
         WRITE(*,*) 'HALLO, HIER ALLE KOMMUNIKATIONSTYPEN EINZELN AUFLISTEN !!!'
   CALL SCARC_RECEIVE (NCOM_MATV,  ILEVEL)  !Aufruf 1, IMV
   CALL SCARC_EXCHANGE (NCOM_MATV,  ILEVEL, IMV)
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

      IF (SCARC_DEBUG>=4) write(9,1000) I_LO,1,K_LO,X_LO(I_LO,1,K_LO),I1_HI,I2_HI,K1_HI,K2_HI
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
         IF (SCARC_DEBUG>=4) write(9,1000) I_LO,J_LO,K_LO,X_LO(I_LO,J_LO,K_LO),I1_HI,I2_HI,J1_HI,J2_HI,K1_HI,K2_HI

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

CALL SCARC_COMPARE_SINGLE_LEV (X_LO, 'SARC',  'PX_LO ',ILEVEL_HI-1)
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
         write(9,1000) I1_HI,1,K1_HI,X_LO(I_LO,1,K_LO)
         write(9,1000) I2_HI,1,K1_HI,X_LO(I_LO,1,K_LO)
         write(9,1000) I1_HI,1,K2_HI,X_LO(I_LO,1,K_LO)
         write(9,1000) I2_HI,1,K2_HI,X_LO(I_LO,1,K_LO)
      ENDIF
   ENDDO
ENDDO
CALL SCARC_COMPARE_SINGLE_LEV (X_HI, 'SARC',  'PX_LO ',ILEVEL_HI)

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

CALL SCARC_COMPARE_SINGLE_LEV (X_LO, 'SARC',  'PX_LO ',ILEVEL_HI-1)
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
            write(9,1000) I1_HI,J1_HI,K1_HI,X_LO(I_LO,J_LO,K_LO)
            write(9,1000) I2_HI,J1_HI,K1_HI,X_LO(I_LO,J_LO,K_LO)
            write(9,1000) I1_HI,J2_HI,K1_HI,X_LO(I_LO,J_LO,K_LO)
            write(9,1000) I2_HI,J2_HI,K1_HI,X_LO(I_LO,J_LO,K_LO)
            write(9,1000) I1_HI,J1_HI,K2_HI,X_LO(I_LO,J_LO,K_LO)
            write(9,1000) I2_HI,J1_HI,K2_HI,X_LO(I_LO,J_LO,K_LO)
            write(9,1000) I1_HI,J2_HI,K2_HI,X_LO(I_LO,J_LO,K_LO)
            write(9,1000) I2_HI,J2_HI,K2_HI,X_LO(I_LO,J_LO,K_LO)
         ENDIF

      ENDDO
   ENDDO
ENDDO
CALL SCARC_COMPARE_SINGLE_LEV (X_HI, 'SARC',  'PX_LO ',ILEVEL_HI)

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

CALL SCARC_COMPARE_SINGLE_LEV (F, 'SARC',  'FRES1 ',ILEVEL)
DO IW_FACE=1,SL%NW_FACE

   IF (SL%IJKW_FACE(9,IW_FACE)==0) THEN

      I=SL%IJKW_FACE(1,IW_FACE)
      J=SL%IJKW_FACE(2,IW_FACE)
      K=SL%IJKW_FACE(3,IW_FACE)

      F(I,J,K)=0.0_EB

      if (SCARC_DEBUG>=4) WRITE(9,*) 'F(',I,',',J,',',K,')=',F(I,J,K)

   ENDIF

ENDDO
   CALL SCARC_COMPARE_SINGLE_LEV (F, 'SARC',  'FRES2 ',ILEVEL)

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
 
IF (A2 == 0.0_EB) THEN
   IF (A1 /= 1.0_EB) THEN
      DO K = 1, SL%KBAR
         DO I = 1, SL%IBAR
            Y (I, 1, K) = A1 * X (I, 1, K)
         ENDDO
      ENDDO
   ELSE
      DO K = 1, SL%KBAR
         DO I = 1, SL%IBAR
            Y (I, 1, K) = X (I, 1, K)
         ENDDO
      ENDDO
   ENDIF
ELSE
   IF (A1 /= 1.0_EB) THEN
      DO K = 1, SL%KBAR
         DO I = 1, SL%IBAR
            Y (I, 1, K) = A1 * X (I, 1, K) + A2 * Y (I, 1, K)
         ENDDO
      ENDDO
   ELSE
      DO K = 1, SL%KBAR
         DO I = 1, SL%IBAR
            Y (I, 1, K) = X (I, 1, K) + A2 * Y (I, 1, K)
         ENDDO
      ENDDO
   ENDIF
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

IF (A2 == 0.0_EB) THEN
   IF (A1 /= 1.0_EB) THEN
      DO K = 1, SL%KBAR
         DO J = 1, SL%JBAR
            DO I = 1, SL%IBAR
               Y (I, J, K) = A1 * X (I, J, K)
            ENDDO
         ENDDO
      ENDDO
   ELSE
      DO K = 1, SL%KBAR
         DO J = 1, SL%JBAR
            DO I = 1, SL%IBAR
               Y (I, J, K) = X (I, J, K)
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
            ENDDO
         ENDDO
      ENDDO
   ELSE
      DO K = 1, SL%KBAR
         DO J = 1, SL%KBAR
            DO I = 1, SL%IBAR
               Y (I, J, K) = X (I, J, K) + A2 * Y (I, J, K)
            ENDDO
         ENDDO
      ENDDO
   ENDIF
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
 
DO K = 1, SL%KBAR
   DO I = 1, SL%IBAR
      Y (I, 1, K) = X (I, 1, K)
   ENDDO
ENDDO
 
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
 
DO K = 1, SL%KBAR
   DO J = 1, SL%JBAR
      DO I = 1, SL%IBAR
         Y (I, J, K) = X (I, J, K)
      ENDDO
   ENDDO
ENDDO
 
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
 
TNOW_SCALPROD2D = SECOND()
 
SL => S%SLEVEL(ILEVEL)
 
SP = 0.0_EB
DO K = 1, SL%KBAR
   DO I = 1, SL%IBAR
      SP = SP + X (I, 1, K) * Y (I, 1, K)
   ENDDO
ENDDO
 
IF (NMESHES > 1) THEN

   S%SP_LOCAL0 = 0.0_EB
   S%SP_LOCAL0 (NM) = SP
   S%SP_GLOBAL0 = 0.0_EB

   CALL MPI_ALLREDUCE (S%SP_LOCAL0(1), S%SP_GLOBAL0(1), NMESHES, &
                       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
   SP = 0.0_EB
   DO IM = 1, NMESHES
      SP = SP + S%SP_GLOBAL0(IM)
   ENDDO
 
ENDIF

TUSED_SCARC(16,NM)=TUSED_SCARC(16,NM)+SECOND()-TNOW_SCALPROD2D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_SCALPROD2D
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
 
TNOW_SCALPROD3D = SECOND()
 
SL => S%SLEVEL(ILEVEL)

SP = 0.0_EB
DO K = 1, SL%KBAR
   DO J = 1, SL%JBAR
      DO I = 1, SL%IBAR
         SP = SP + X (I, J, K) * Y (I, J, K)
      ENDDO
   ENDDO
ENDDO
 
IF (NMESHES > 1) THEN
 
   S%SP_LOCAL0 = 0.0_EB
   S%SP_LOCAL0 (NM) = SP
   S%SP_GLOBAL0 = 0.0_EB

   CALL MPI_ALLREDUCE (S%SP_LOCAL0(1), S%SP_GLOBAL0(1), NMESHES, &
                       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
 
   SP = 0.0_EB
   DO IM = 1, NMESHES
      SP = SP + S%SP_GLOBAL0(IM)
   ENDDO
 
ENDIF

TUSED_SCARC(16,NM)=TUSED_SCARC(16,NM)+SECOND()-TNOW_SCALPROD3D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_SCALPROD3D
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
REAL(EB):: TEST_LOCAL0, TEST_GLOBAL0
REAL(EB):: TNOW_L2NORM2D
 
TNOW_L2NORM2D = SECOND()
 
SL => S%SLEVEL(ILEVEL)

!!! build local scalar product (x,x)
S%SP_LOCAL = 0.0_EB
DO K = 1, SL%KBAR
   DO I = 1, SL%IBAR
      S%SP_LOCAL = S%SP_LOCAL + X (I, 1, K) * X (I, 1, K)
   ENDDO
ENDDO
 
!!! scale with number of local cells
IF (ITYPE == NTYPE_LOCAL) THEN
 
   SP = Sqrt (S%SP_LOCAL/REAL(SL%NCELLS_LOCAL, EB))
 
!!! sum up globally and scale with global number of nodes
ELSE IF (ITYPE == NTYPE_GLOBAL) THEN
 
   IF (NMESHES > 1) THEN
 
      S%SP_LOCAL0 = 0.0_EB
      S%SP_LOCAL0 (NM) = S%SP_LOCAL
 
      S%SP_GLOBAL0 = 0.0_EB
      CALL MPI_ALLREDUCE (S%SP_LOCAL0(1), S%SP_GLOBAL0(1), NMESHES, &
                          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
      CALL MPI_ALLREDUCE (TEST_LOCAL0, TEST_GLOBAL0, NMESHES, &
                          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
 
      S%SP_GLOBAL = 0.0_EB
      DO IM = 1, NMESHES
         S%SP_GLOBAL = S%SP_GLOBAL + S%SP_GLOBAL0(IM)
      ENDDO
   ELSE
      S%SP_GLOBAL = S%SP_LOCAL
   ENDIF
   IF (SCARC_DEBUG>=2) THEN
      WRITE(9,*) '  SP_GLOBAL=',S%SP_GLOBAL
      WRITE(9,*) 'TEST_GLOBAL=',TEST_GLOBAL0
   ENDIF
   SP = Sqrt (S%SP_GLOBAL/REAL(SL%NCELLS_GLOBAL, EB))

ELSE
 
   WRITE (*,*) 'Wrong type for SCARC_L2NORM2D ', ITYPE
 
ENDIF

TUSED_SCARC(17,NM)=TUSED_SCARC(17,NM)+SECOND()-TNOW_L2NORM2D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_L2NORM2D
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
 
TNOW_L2NORM3D = SECOND()
 
SL => S%SLEVEL(ILEVEL)
 
!!! build local scalar product (x,x)
S%SP_LOCAL = 0.0_EB
DO K = 1, SL%KBAR
   DO J = 1, SL%JBAR
      DO I = 1, SL%IBAR
         S%SP_LOCAL = S%SP_LOCAL + X (I, J, K) * X (I, J, K)
      ENDDO
   ENDDO
ENDDO
 
IF (SCARC_DEBUG >= 2) write (9,*) 'S%SP_LOCAL=', S%SP_LOCAL
 
!!! scale with number of local cells
IF (ITYPE == NTYPE_LOCAL) THEN
 
   SP = Sqrt (S%SP_LOCAL/REAL(SL%NCELLS_LOCAL, EB))
 
!!! sum up globally and scale with global number of nodes
ELSE IF (ITYPE == NTYPE_GLOBAL) THEN
 
   IF (NMESHES > 1) THEN
 
      S%SP_LOCAL0 = 0.0_EB
      S%SP_LOCAL0 (NM) = S%SP_LOCAL
 
      S%SP_GLOBAL0 = 0.0_EB
      CALL MPI_ALLREDUCE (S%SP_LOCAL0(1), S%SP_GLOBAL0(1), NMESHES, &
                          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
 
      S%SP_GLOBAL = 0.0_EB
      DO IM = 1, NMESHES
         S%SP_GLOBAL = S%SP_GLOBAL + S%SP_GLOBAL0(IM)
      ENDDO
   ELSE
      S%SP_GLOBAL = S%SP_LOCAL
   ENDIF
 
   SP = Sqrt (S%SP_GLOBAL/REAL(SL%NCELLS_GLOBAL, EB))
 
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

INTEGER :: NM,ILEVEL, IC, NC, IZERO
REAL(EB):: TNOW_GSTRIX2D_INIT
 
TNOW_GSTRIX2D_INIT = SECOND()
 
SL => S%SLEVEL(ILEVEL)

ALLOCATE (SL%DD(1:SL%NCELLS_LOCAL), STAT=IZERO)
CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX2D', 'SL%DD', IZERO)

ALLOCATE (SL%DU(1:SL%NCELLS_LOCAL), STAT=IZERO)
CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX2D', 'SL%DU', IZERO)

ALLOCATE (SL%DL(1:SL%NCELLS_LOCAL), STAT=IZERO)
CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX2D', 'SL%DL', IZERO)

ALLOCATE (SL%LD(1:SL%NCELLS_LOCAL), STAT=IZERO)
CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX2D', 'SL%LD', IZERO)

ALLOCATE (SL%DAUX(1:SL%NCELLS_LOCAL), STAT=IZERO)
CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX2D', 'SL%DAUX', IZERO)


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
   IF (SCARC_DEBUG>=2) write(9,1000) ILEVEL,IC,SL%DU(IC),IC,SL%DD(IC)
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
   write(9,*) '===================== STARTING PRECON_GSTRIX2D'
   write(9,*) 'N=',N
   write(9,*) 'M=',M
   WRITE (9,*) 'BEFORE PRECON_GSTRIX: Y'
   CALL SCARC_COMPARE_SINGLE (DY, 'SARC', 'Y     ')
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
   WRITE (9,*) 'AFTER  PRECON_GSTRIX: Y'
   CALL SCARC_COMPARE_SINGLE (DY, 'SARC',  'Y     ')
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

INTEGER :: NM,ILEVEL, IC, NC, IZERO
REAL(EB):: TNOW_GSTRIX3D_INIT
 
TNOW_GSTRIX3D_INIT = SECOND()
 
write(*,*) 'MUSS NOCH AUF 3D angepasst werden !!!'
write(9,*) 'MUSS NOCH AUF 3D angepasst werden !!!'
stop

SL => S%SLEVEL(ILEVEL)

ALLOCATE (SL%DD(1:SL%NCELLS_LOCAL), STAT=IZERO)
CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX3D', 'SL%DD', IZERO)

ALLOCATE (SL%DU(1:SL%NCELLS_LOCAL), STAT=IZERO)
CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX3D', 'SL%DU', IZERO)

ALLOCATE (SL%DL(1:SL%NCELLS_LOCAL), STAT=IZERO)
CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX3D', 'SL%DL', IZERO)

ALLOCATE (SL%LD(1:SL%NCELLS_LOCAL), STAT=IZERO)
CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX3D', 'SL%LD', IZERO)

ALLOCATE (SL%DAUX(1:SL%NCELLS_LOCAL), STAT=IZERO)
CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX3D', 'SL%DAUX', IZERO)


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
   IF (SCARC_DEBUG>=2) write(9,1000) ILEVEL,IC,SL%DU(IC),IC,SL%DD(IC)
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
   write(9,*) '===================== STARTING PRECON_GSTRIX3D'
   write(9,*) 'N=',N
   write(9,*) 'M=',M
   WRITE (9,*) 'BEFORE PRECON_GSTRIX: Y'
   CALL SCARC_COMPARE_SINGLE (DY, 'SARC',  'Y     ')
ENDIF

write(*,*) 'MUSS NOCH AUF 3D angepasst werden !!!'
write(9,*) 'MUSS NOCH AUF 3D angepasst werden !!!'
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
   WRITE (9,*) 'AFTER  PRECON_GSTRIX: Y'
   CALL SCARC_COMPARE_SINGLE (DY, 'SARC',  'Y     ')
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
   WRITE (9,*) 'Starting SCARC_PRECON_MG2D'
   WRITE (9,*) 'IBAR=:', SLHI%IBAR
   WRITE (9,*) 'JBAR=:', SLHI%JBAR
   WRITE (9,*) 'KBAR=:', SLHI%KBAR
   WRITE (9,*) 'SCARC_NIT_MG=:', SCARC_NIT_MG
   WRITE (9,*) 'SCARC_EPS_MG=:', SCARC_EPS_MG
   WRITE (9,*) 'SCARC_EPS_REL=:', SCARC_EPS_REL
   CALL SCARC_COMPARE_SINGLE (SLHI%X2, 'SARC',  'X2    ')
   CALL SCARC_COMPARE_SINGLE (SLHI%F2, 'SARC',  'F2    ')
ENDIF
 
!!!
!!! initialize some method parameters
!!!
NREQ_FACE = 0           ! necessary ?
NREQ_EDGE = 0
NREQ_DIAG = 0

BCONV =.FALSE.
BDIVG =.FALSE.
ICYCLE=1      ! V-cycle
 
IBAR0=SLHI%IBAR
JBAR0=SLHI%JBAR
KBAR0=SLHI%KBAR

!!! save cycle counts for MG-iteration
S%KCYCLE(2,S%NLEVEL_MAX)=1
DO ILEVEL=S%NLEVEL_MIN+1,S%NLEVEL_MAX-1
   IF (ICYCLE==0) THEN
      S%KCYCLE(2,ILEVEL)=2
   ELSE
      S%KCYCLE(2,ILEVEL)=ICYCLE
   ENDIF
ENDDO
  

!!! calculate initial defect on finest level, Y = B - A*X, and get l2-norm of it
!WRITE(9,*) 'KOMMUNIKATIONSPARAMETER 3 prüfen, anders als bei CG !!!'
CALL SCARC_MATVEC2D (SLHI%X2, SLHI%D2, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, 3)
CALL SCARC_VECADD2D (SLHI%F2, SLHI%D2, 1.0_EB, -1.0_EB, ILEVEL, NM)
CALL SCARC_L2NORM2D (SLHI%D2, SLHI%RESIN_MG, ILEVEL, NM, NTYPE_GLOBAL)
IF (SCARC_DEBUG >= 1) WRITE(9,1000) 0, ILEVEL, SL%RESIN_MG
IF (NM == 1)          WRITE(*,1000) 0, ILEVEL, SL%RESIN_MG
 
!WRITE(*,*) 'Info muss noch an Master geschickt werden'
!CALL SCARC_MASTER_INFO (SLHI%ITE, SLHI%RESIN_MG)

 
!!! ---------------------------------------------------------------------------
!!! start MG-iteration
!!! ---------------------------------------------------------------------------
MG_LOOP2D: DO ITE = 1, SCARC_NIT_MG
 
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

      CALL SCARC_SMOOTHER2D(SLHI%X2, SLHI%F2, SLHI%D2, BMATVEC, ILEVEL, NM)

      !!! print level information - still missing
      !CALL SCARC_INFO_LEVEL(ILEVEL)

      !!! perform restriction to coarser grid: F:= rest(D)
      CALL SCARC_COMPARE_SINGLE_LEV (SLHI%D2, 'SARC',  'RESTD ',ILEVEL)
      CALL SCARC_RESTRICTION2D(SLLO%F2, SLHI%D2, ILEVEL-1, NM)
      CALL SCARC_COMPARE_SINGLE_LEV (SLLO%F2, 'SARC',  'RESTF ',ILEVEL-1)

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
   ILEVEL=S%NLEVEL_MIN
   CALL SCARC_COARSE2D(2, ILEVEL, NM)


   !!! ------------------------------------------------------------------------
   !!! Postsmoothing
   !!! ------------------------------------------------------------------------
120   IF (ILEVEL/=S%NLEVEL_MAX) THEN
 
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
      CALL SCARC_VECADD2D(SLHI%D2, SLHI%X2, 1.0_EB, 1.0_EB, ILEVEL, NM)

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
   CALL SCARC_MATVEC2D (SLHI%X2, SLHI%D2, 1.0_EB, 0.0_EB, NM, NTYPE_GLOBAL, ILEVEL, 3)      ! 3 korrekt ?
   CALL SCARC_VECADD2D (SLHI%F2, SLHI%D2, 1.0_EB, -1.0_EB, ILEVEL, NM)
   CALL SCARC_L2NORM2D (SLHI%D2, SLHI%RES_MG, ILEVEL, NM, NTYPE_GLOBAL)

   IF (SCARC_DEBUG>=2) WRITE(9,*) 'SCARC-Multigrid, iteration ',ITE,': residuum=',SLHI%RES_MG

   !!! ------------------------------------------------------------------------
   !!! Send error information to Master  - still missing
   !!! ------------------------------------------------------------------------
   !CALL SCARC_INFO_MASTER 
   !CALL SCARC_INFO_LEVEL(ILEVEL)

   !!! ------------------------------------------------------------------------
   !!! Convergence or divergence ?
   !!! ------------------------------------------------------------------------
   IF (SCARC_DEBUG >= 1) write (9,1000) ITE, ILEVEL, SL%RES_MG
   IF (SCARC_DEBUG >= 1) write (*,1000) ITE, ILEVEL, SL%RES_MG
   IF (NM == 1)          write (*,1000) ITE, ILEVEL, SL%RES_MG
   IF (SCARC_BREL) THEN
      IF (SL%RES_MG <= SL%RESIN_MG*SCARC_EPS_MG) BCONV = .TRUE.
   ELSE
      IF (SL%RES_MG <= SCARC_EPS_MG .AND. SL%RES_MG <= SL%RESIN_MG*SCARC_EPS_REL) BCONV = .TRUE.
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
IF (SCARC_DEBUG >= 1) WRITE (9,2000) ITE0, ILEVEL, SL%RES_MG, SL%CAPPA_MG 
!IF (NM == 1)          WRITE (*,2000) ITE0, ILEVEL, SL%RES_MG, SL%CAPPA_MG
 
DO K = 1, SLHI%KBAR
   DO I = 1, SLHI%IBAR
      Y (I, 1, K) = SLHI%X2(I, 1, K)
   ENDDO
ENDDO
 
TUSED_SCARC(7,NM)=TUSED_SCARC(7,NM)+SECOND()-TNOW_PRECON_MG2D
TUSED_SCARC(0,NM)=TUSED_SCARC(0,NM)+SECOND()-TNOW_PRECON_MG2D

1000 FORMAT (/,'     SCARC_PRECON_MG: #ite= ',i4,': level=',i4,': res=',e12.4,/)
2000 FORMAT (/,'     SCARC_PRECON_MG: #ite= ',i4,': level=',i4,': res=',e12.4,': rate=',e12.4,/)
END SUBROUTINE SCARC_PRECON_MG2D


 
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
INTEGER :: IMIN, IMAX, JMIN, JMAX, KMIN, KMAX
INTEGER :: IERR
INTEGER :: ILEN_FACE, ILEN_EDGE, ILEN_DIAG
INTEGER :: TAG_FACE, TAG_EDGE, TAG_DIAG
REAL(EB):: TNOW_RECEIVE

TYPE (SCARC_TYPE),  POINTER ::  SNM
TYPE (OSCARC_TYPE), POINTER :: OSNM
TYPE (SCARC_LEVEL_TYPE),  POINTER ::  SNML
TYPE (OSCARC_LEVEL_TYPE), POINTER :: OSNML

TNOW_RECEIVE = SECOND()

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

         TAG_FACE = TAGS_FACE(NM,NOM)
    

         !!! Initialize the communication structures for receiving face data
         INIT_FACE_IF: IF (CODE==NCOM_INIT) THEN

            NREQ_FACE = NREQ_FACE+1
            CALL MPI_IRECV(OSNML%IJKW_FACE(1,1),17*OSNML%NW_FACE,MPI_INTEGER,SNODE, &
                           TAG_FACE,MPI_COMM_WORLD,REQ_FACE(NREQ_FACE),IERR)

            IF (SNML%NIC_FACE(NM, NOM) > 0) THEN
               ILEN_FACE=(MAX(SNML%NIC_FACE(NM, NOM), SNML%NIC_FACE(NOM, NM))+2)*2+1
               ALLOCATE (OSNM%RECV_FACE(ILEN_FACE))
               OSNM%RECV_FACE = 0.0_EB
            ENDIF

         ENDIF INIT_FACE_IF
   
         !!! Perform exchange for matrix-vector multiplication (only face exchange needed)
         MATV_FACE_IF: IF (CODE==NCOM_MATV) THEN

            NREQ_FACE = NREQ_FACE+1
            CALL MPI_IRECV(OSNM%RECV_FACE(1),SIZE(OSNM%RECV_FACE),MPI_DOUBLE_PRECISION,&
                        SNODE,TAG_FACE,MPI_COMM_WORLD,REQ_FACE(NREQ_FACE),IERR)

         ENDIF MATV_FACE_IF
   
         !!! Perform full exchange including edge (3D)
         FULL_FACE_IF: IF (CODE==NCOM_FULL) THEN

            NREQ_FACE = NREQ_FACE+1
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

      !!! Initialize the communication structures for receiving face data
         INIT_DIAG_IF: IF (CODE==NCOM_INIT) THEN

            NREQ_DIAG = NREQ_DIAG+1
            CALL MPI_IRECV(OSNML%IJKW_DIAG(1,1),17*OSNML%NW_DIAG,MPI_INTEGER,SNODE, &
                           TAG_DIAG,MPI_COMM_WORLD,REQ_DIAG(NREQ_DIAG),IERR)

            IF (SNML%NIC_DIAG(NM, NOM) > 0) THEN
               ILEN_DIAG=(MAX(SNML%NIC_DIAG(NM, NOM), SNML%NIC_DIAG(NOM, NM))+2)*2+1   ! extended
               ALLOCATE (OSNM%RECV_DIAG(ILEN_DIAG))
               OSNM%RECV_DIAG = 0.0_EB
            ENDIF

         ENDIF INIT_DIAG_IF
   
   
      !!! Perform full exchange including edge (3D)
         FULL_DIAG_IF: IF (CODE==NCOM_FULL) THEN

            NREQ_DIAG = NREQ_DIAG+1
            CALL MPI_IRECV(OSNM%RECV_DIAG(1),SIZE(OSNM%RECV_DIAG),MPI_DOUBLE_PRECISION,&
                           SNODE,TAG_DIAG,MPI_COMM_WORLD,REQ_DIAG(NREQ_DIAG),IERR)

         ENDIF FULL_DIAG_IF

      ENDIF RECEIVE_DIAG_IF

   ENDDO RECEIVE_OMESH_LOOP
ENDDO RECEIVE_MESH_LOOP


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

TNOW_EXCHANGE = SECOND()

SNM => SCARC(NM)
 
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


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! FACE communication
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FACE_IF: IF (SNML%NIC_FACE(NOM,NM)/=0 .AND. SNML%NIC_FACE(NM,NOM)/=0)  THEN
    
         OSNM  => SCARC(NM)%OSCARC(NOM)                 ! corresponds to M3
         OSNML => SCARC(NM)%OSCARC(NOM)%SLEVEL(ILEVEL)  ! ... for the level 'ILEVEL'

         TAG_FACE = TAGS_FACE(NM,NOM)
    
         !!! Initialize the communication structures for sending data
         EXCHANGE_INIT_FACE_IF: IF (CODE==NCOM_INIT) THEN
    
            ! initialize communication structures for face exchange
            IF (RNODE/=SNODE) THEN
              NREQ_FACE = NREQ_FACE+1
              CALL MPI_ISEND(SNML%IJKW_FACE(1,1),17*SNML%NW_FACE,MPI_INTEGER,SNODE, &
                             TAG_FACE,MPI_COMM_WORLD,REQ_FACE(NREQ_FACE),IERR)
               ILEN_FACE=(MAX(SNML%NIC_FACE(NM, NOM), SNML%NIC_FACE(NOM, NM))+2)*2+1   ! extended
               ALLOCATE (OSNM%SEND_FACE(ILEN_FACE))
               OSNM%SEND_FACE = 0.0_EB
            ENDIF

         ENDIF EXCHANGE_INIT_FACE_IF
   
      !!! Perform exchange for matrix-vector multiplication (only face exchange needed)
         EXCHANGE_MATV_FACE_IF: IF (CODE==NCOM_MATV) THEN
            IF (RNODE/=SNODE) THEN

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
                        LL = LL+2
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO PACK_SEND_FACE0
            OSNM%SEND_FACE(IWW*2+1) = -999.0_EB
            NREQ_FACE=NREQ_FACE+1
            CALL MPI_ISEND(OSNM%SEND_FACE(1),IWW*2+1,MPI_DOUBLE_PRECISION,SNODE, &
                           TAG_FACE,MPI_COMM_WORLD,REQ_FACE(NREQ_FACE),IERR)

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
                        OSNM%SEND_FACE(LL+1) = REAL(IW,EB)
                        OSNM%SEND_FACE(LL+2) = SNML%Z(II,JJ,KK)
                        LL = LL+2
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO PACK_SEND_FACE
            OSNM%SEND_FACE(IWW*2+1) = -999.0_EB
            NREQ_FACE=NREQ_FACE+1
            CALL MPI_ISEND(OSNM%SEND_FACE(1),IWW*2+1,MPI_DOUBLE_PRECISION,SNODE, &
                           TAG_FACE,MPI_COMM_WORLD,REQ_FACE(NREQ_FACE),IERR)


         ELSE
                write(*,*) 'Not yet implemented'
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
                        LL = LL+1
                     ENDDO
                  ENDDO
               ENDDO
            NREQ_EDGE=NREQ_EDGE+1
            CALL MPI_ISEND(OSNM%SEND_EDGE(1),ILEN_EDGE,MPI_DOUBLE_PRECISION,SNODE, &
                           TAG_EDGE,MPI_COMM_WORLD,REQ_EDGE(NREQ_EDGE),IERR)
            ENDIF

         ELSE
                write(*,*) 'Not yet implemented'
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
               ALLOCATE (OSNM%SEND_DIAG(ILEN_DIAG))
               OSNM%SEND_DIAG = 0.0_EB
            ENDIF

         ENDIF EXCHANGE_INIT_DIAG_IF
   
      !!! Perform full exchange including edge (3D!) and vertex information
         EXCHANGE_FULL_DIAG_IF: IF (CODE==NCOM_FULL) THEN
            IF (RNODE/=SNODE) THEN

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
                           LL = LL+2
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO PACK_SEND_DIAG
               OSNM%SEND_DIAG(IWW*2+1) = -999.0_EB
               NREQ_DIAG=NREQ_DIAG+1
               CALL MPI_ISEND(OSNM%SEND_DIAG(1),IWW*2+1,MPI_DOUBLE_PRECISION,SNODE, &
                              TAG_DIAG,MPI_COMM_WORLD,REQ_DIAG(NREQ_DIAG),IERR)

            ENDIF

         ELSE
                write(*,*) 'Not yet implemented'
         ENDIF EXCHANGE_FULL_DIAG_IF

      ENDIF DIAG_IF
   
   ENDDO EXCHANGE_RECV_MESH_LOOP
ENDDO EXCHANGE_SEND_MESH_LOOP


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Information from Mesh NM is RECV'ed by Mesh NOM.  
!!! NOM is the receiver, NM is the sender.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CALL MPI_WAITALL(NREQ_FACE,REQ_FACE(1:NREQ_FACE),STAT_FACE,IERR)
IF (BEDGE) CALL MPI_WAITALL(NREQ_EDGE,REQ_EDGE(1:NREQ_EDGE),STAT_EDGE,IERR)
IF (BDIAG) CALL MPI_WAITALL(NREQ_DIAG,REQ_DIAG(1:NREQ_DIAG),STAT_DIAG,IERR)

 
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

               UNPACK_RECV_FACE0: DO
                  IW = NINT(OSNOM%RECV_FACE(LL+1))
                  IF (IW==-999) EXIT UNPACK_RECV_FACE0
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
                  I=SNOML%IJKW_FACE(6,IW)
                  J=SNOML%IJKW_FACE(7,IW)
                  K=SNOML%IJKW_FACE(8,IW)
                  ISUM = (SNOML%IJKW_FACE(13,IW)-SNOML%IJKW_FACE(10,IW)+1) * &
                         (SNOML%IJKW_FACE(14,IW)-SNOML%IJKW_FACE(11,IW)+1) * &
                         (SNOML%IJKW_FACE(15,IW)-SNOML%IJKW_FACE(12,IW)+1)

                  IF (IMV==1 .OR. IMV==5) THEN
                     yold=SNOML%Y(I,J,K)
                     SNOML%Y(I, J, K) = SNOML%Y(I, J, K) + ASUB * ZSUM/REAL(ISUM,EB)
                  ELSE IF (IMV==2 .OR. IMV==6) THEN
                     yold=SNOML%G(I,J,K)
                     SNOML%G(I, J, K) = SNOML%G(I, J, K) + ASUB * ZSUM/REAL(ISUM,EB)
                  ELSE IF (IMV==3) THEN
                     yold=SNOML%R(I,J,K)
                     SNOML%R(I, J, K) = SNOML%R(I, J, K) + ASUB * ZSUM/REAL(ISUM,EB)
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
               LL = 0
               DO KK=KMIN,KMAX
                  DO JJ=JMIN,JMAX
                     DO II=IMIN,IMAX
                        SNOML%Z(II, JJ, KK) = OSNOM%RECV_EDGE(LL+1)
                        LL = LL+1
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF
         ENDIF RECV_EDGE_FULL_IF

      ENDIF RECV_EDGE_IF

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
               UNPACK_RECV_DIAG: DO
                  IW = NINT(OSNOM%RECV_DIAG(LL+1))
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
               ENDDO UNPACK_RECV_DIAG

            ENDIF

         ENDIF RECV_DIAG_FULL_IF

      ENDIF RECV_DIAG_IF

   
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
 
TYPE (SCARC_TYPE),  POINTER ::  S
TYPE (SCARC_LEVEL_TYPE),  POINTER ::  SL
TYPE (MESH_TYPE),  POINTER ::  M

TNOW_GHOSTCELLS = SECOND()

S => SCARC(NM)
SL => S%SLEVEL(S%NLEVEL_MAX)
M  => MESHES(NM)
 
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
 
DO IW = 1, M%NEWC

   I = M%IJKW(6,IW)
   J = M%IJKW(7,IW)
   K = M%IJKW(8,IW)

   IOR0 = M%IJKW(4,IW)

   SELECT CASE(IOR0)
      CASE(1)
         IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
            SL%F(I,J,K) = SL%F(I,J,K) - 2.0_EB * SL%DXI2 * M%BXS(1,K)            ! Dirichlet
            IF (SCARC_DEBUG>=6) WRITE(9,1000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
         ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
            SL%F(I,J,K) = SL%F(I,J,K) + SL%DXI * M%BXS(1,K)                       ! Neumann
            IF (SCARC_DEBUG>=6) WRITE(9,2000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
         ENDIF
      CASE(-1)
         IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
            SL%F(I,J,K) = SL%F(I,J,K) - 2.0_EB * SL%DXI2 *M%BXF(1,K)            ! Dirichlet
            IF (SCARC_DEBUG>=6) WRITE(9,1000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
         ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
            SL%F(I,J,K) = SL%F(I,J,K) - SL%DXI *M%BXF(1,K)                      ! Neumann
            IF (SCARC_DEBUG>=6) WRITE(9,2000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
         ENDIF
      CASE(3)
         IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
            SL%F(I,J,K) = SL%F(I,J,K) - 2.0_EB * SL%DZI2 * M%BZS(I,1)            ! Dirichlet
            IF (SCARC_DEBUG>=6) WRITE(9,1000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
         ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
            SL%F(I,J,K) = SL%F(I,J,K) + SL%DZI * M%BZS(I,1)                       ! Neumann
            IF (SCARC_DEBUG>=6) WRITE(9,2000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
         ENDIF
      CASE(-3)
         IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
            SL%F(I,J,K) = SL%F(I,J,K) - 2.0_EB * SL%DZI2 * M%BZF(I,1)            ! Dirichlet
            IF (SCARC_DEBUG>=6) WRITE(9,1000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
         ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
            SL%F(I,J,K) = SL%F(I,J,K) - SL%DZI  * M%BZF(I,1)                     ! Neumann
            IF (SCARC_DEBUG>=6) WRITE(9,2000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
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
            IF (SCARC_DEBUG>=6) WRITE(9,1000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
         ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
            SL%F(I,J,K) = SL%F(I,J,K) + SL%DXI * M%BXS(1,K)                       ! Neumann
            IF (SCARC_DEBUG>=6) WRITE(9,2000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
         ENDIF
      CASE(-1)
         IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
            SL%F(I,J,K) = SL%F(I,J,K) - 2.0_EB * SL%DXI2 *M%BXF(1,K)            ! Dirichlet
            IF (SCARC_DEBUG>=6) WRITE(9,1000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
         ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
            SL%F(I,J,K) = SL%F(I,J,K) - SL%DXI *M%BXF(1,K)                      ! Neumann
            IF (SCARC_DEBUG>=6) WRITE(9,2000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
         ENDIF
      CASE(3)
         IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
            SL%F(I,J,K) = SL%F(I,J,K) - 2.0_EB * SL%DZI2 * M%BZS(I,1)            ! Dirichlet
            IF (SCARC_DEBUG>=6) WRITE(9,1000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
         ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
            SL%F(I,J,K) = SL%F(I,J,K) + SL%DZI * M%BZS(I,1)                       ! Neumann
            IF (SCARC_DEBUG>=6) WRITE(9,2000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
         ENDIF
      CASE(-3)
         IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
            SL%F(I,J,K) = SL%F(I,J,K) - 2.0_EB * SL%DZI2 * M%BZF(I,1)            ! Dirichlet
            IF (SCARC_DEBUG>=6) WRITE(9,1000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
         ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
            SL%F(I,J,K) = SL%F(I,J,K) - SL%DZI  * M%BZF(I,1)                     ! Neumann
            IF (SCARC_DEBUG>=6) WRITE(9,2000) IW, IOR0, I, J, K, SL%F(I,J,K), IW, M%PRESSURE_BC_INDEX(IW)
         ENDIF
   END SELECT


ENDDO
 
 
TUSED_SCARC(27,NM)=TUSED_SCARC(27,NM)+SECOND()-TNOW_SETBDRY3D
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_SETBDRY3D
1000 FORMAT(i3,': IOR=',i3,': Dirichlet :F(',i2,',',i2,',',i2,')=',e25.16,i5,i5)
2000 FORMAT(i3,': IOR=',i3,': Neumann   :F(',i2,',',i2,',',i2,')=',e25.16,i5,i5)
END SUBROUTINE SCARC_SETBDRY3D
 
 

 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Update different quantities according to given type
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_UPDATE (MYID, CUPDATE_TYPE)
INTEGER :: MYID, ILEVEL, NM
CHARACTER(4):: CUPDATE_TYPE
REAL(EB):: TNOW_UPDATE

TNOW_UPDATE=SECOND()

NM = MYID + 1
S => SCARC(NM)
 
UPDATE_LEVEL_LOOP: DO ILEVEL=S%NLEVEL_MAX,S%NLEVEL_MIN,-1

   SELECT CASE(CUPDATE_TYPE)
      CASE('FULL')
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
      CASE('VEL ')
         CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%U, NM, ILEVEL, 'U     ')
      IF ( .NOT. TWO_D) CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%V, NM, ILEVEL, 'V     ')
         CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%W, NM, ILEVEL, 'W     ')
      CASE('FV ')
         CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%FVX, NM, ILEVEL, 'FVX   ')
      IF ( .NOT. TWO_D) CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%FVY, NM, ILEVEL, 'FVY   ')
         CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%FVZ, NM, ILEVEL, 'FVZ   ')
      CASE('H   ')
         CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%H, NM, ILEVEL, 'H     ')
      CASE('HS  ')
         CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%HS, NM, ILEVEL, 'HS    ')
   END SELECT

ENDDO UPDATE_LEVEL_LOOP
 
TUSED_SCARC(28,NM)=TUSED_SCARC(28,NM)+SECOND()-TNOW_UPDATE
TUSED_SCARC(0,NM) =TUSED_SCARC(0,NM) +SECOND()-TNOW_UPDATE
 
END SUBROUTINE SCARC_UPDATE



SUBROUTINE SCARC_UPDATE_LEVEL (MYID, ILEVEL, CUPDATE_TYPE)
REAL(EB):: TNOW_UPDATE_LEVEL
INTEGER :: MYID, ILEVEL, NM
CHARACTER(4):: CUPDATE_TYPE

TNOW_UPDATE_LEVEL=SECOND()
 
NM = MYID + 1
 
SELECT CASE(CUPDATE_TYPE)
   CASE('FULL')
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
   CASE('VEL ')
      CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%U, NM, ILEVEL, 'U     ')
      IF ( .NOT. TWO_D) CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%V, NM, ILEVEL, 'V     ')
      CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%W, NM, ILEVEL, 'W     ')
   CASE('FV  ')
      CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%FVX, NM, ILEVEL, 'FVX   ')
      IF ( .NOT. TWO_D) CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%FVY, NM, ILEVEL, 'FVY   ')
      CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%FVZ, NM, ILEVEL, 'FVZ   ')
   CASE('H   ')
      CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%H, NM, ILEVEL, 'H     ')
   CASE('HS  ')
      CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%HS, NM, ILEVEL, 'H     ')
END SELECT
 
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
INTEGER :: NM, ILEVEL, I, J, K
CHARACTER (4) :: CNAME
REAL(EB):: TNOW_UPDATE_QUANTITY
 
TNOW_UPDATE_QUANTITY = SECOND()

SL => S%SLEVEL(ILEVEL)
IF (SCARC_DEBUG>=2) WRITE(9,*) 'UPDATING QUANTITY ', CNAME
 
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
      NREQ_FACE = 0
      NREQ_DIAG = 0
      CALL SCARC_RECEIVE (NCOM_FULL, ILEVEL)   !Aufruf 2,0
      CALL SCARC_EXCHANGE (NCOM_FULL, ILEVEL, NCOM_TYPE_MATVEC)
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
      NREQ_EDGE = 0
      NREQ_DIAG = 0
 
      CALL SCARC_RECEIVE (NCOM_FULL, ILEVEL)  !Aufruf 2,0
      CALL SCARC_EXCHANGE (NCOM_FULL, ILEVEL, NCOM_TYPE_MATVEC)
 
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
END SUBROUTINE SCARC_UPDATE_QUANTITY
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!  update interior boundary values along vertical directions
!!!  to have consistent values all over the domain
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_UPDATE_QUANTITY2 (Z, NM, ILEVEL, CNAME)
REAL (EB), POINTER, DIMENSION (:, :, :) :: Z
INTEGER :: NM, ILEVEL, I, J, K
CHARACTER (4) :: CNAME
REAL(EB):: TNOW_UPDATE_QUANTITY
 
TNOW_UPDATE_QUANTITY = SECOND()

SL => S%SLEVEL(ILEVEL)
IF (SCARC_DEBUG>=2) WRITE(9,*) 'UPDATING QUANTITY ', CNAME
 
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
      NREQ_FACE = 0
      NREQ_DIAG = 0
      CALL SCARC_RECEIVE (NCOM_FULL, ILEVEL)   !Aufruf 2,0
      CALL SCARC_EXCHANGE (NCOM_FULL, ILEVEL, NCOM_TYPE_MATVEC)
      DO K = 0, KBAR + 1
         DO I = 0, IBAR + 1
            Z (I, 1, K) = SL%Z(I, 1, K)
         ENDDO
      ENDDO
 
      DO K = 0, KBAR + 1
         DO I = 0, IBAR + 1
            SL%Z (I, 1, K) = Z (I, 0, K)
         ENDDO
      ENDDO
      NREQ_FACE = 0
      NREQ_DIAG = 0
      CALL SCARC_RECEIVE (NCOM_FULL, ILEVEL)   !Aufruf 2,0
      CALL SCARC_EXCHANGE (NCOM_FULL, ILEVEL, NCOM_TYPE_MATVEC)
      DO K = 0, KBAR + 1
         DO I = 0, IBAR + 1
            Z (I, 0, K) = SL%Z(I, 1, K)
         ENDDO
      ENDDO

      DO K = 0, KBAR + 1
         DO I = 0, IBAR + 1
            SL%Z (I, 1, K) = Z (I, 2, K)
         ENDDO
      ENDDO
      NREQ_FACE = 0
      NREQ_DIAG = 0
      CALL SCARC_RECEIVE (NCOM_FULL, ILEVEL)   !Aufruf 2,0
      CALL SCARC_EXCHANGE (NCOM_FULL, ILEVEL, NCOM_TYPE_MATVEC)
      DO K = 0, KBAR + 1
         DO I = 0, IBAR + 1
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
      NREQ_EDGE = 0
      NREQ_DIAG = 0
 
      CALL SCARC_RECEIVE (NCOM_FULL, ILEVEL)  !Aufruf 2,0
      CALL SCARC_EXCHANGE (NCOM_FULL, ILEVEL, NCOM_TYPE_MATVEC)
 
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
END SUBROUTINE SCARC_UPDATE_QUANTITY2
 
 
 
 
 
LOGICAL FUNCTION MATCH (VAL1, VAL2)
REAL (EB) :: VAL1, VAL2, TOL
TOL = 1.0E-8_EB
MATCH = .FALSE.
IF (Abs(VAL1-VAL2) <= TOL) MATCH = .TRUE.
END FUNCTION MATCH
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! only for debugging reasons ... print also given vector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_COMPARE_SINGLE (X, CROUTINE, CNAME)
REAL (EB), POINTER, DIMENSION (:, :, :) :: X
INTEGER :: II, jj, KK
CHARACTER (4) :: CROUTINE
CHARACTER (6) :: CNAME
 
IF (SCARC_DEBUG >= 2) THEN
   WRITE (9,*) '============================================================='
   WRITE (9,*) '===   COMPARE= ',CNAME,': ROUTINE= ',CROUTINE
   WRITE (9,*) '============================================================='
   IF (TWO_D) THEN
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
   !CALL flush (9)
ENDIF
 
END SUBROUTINE SCARC_COMPARE_SINGLE
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! only for debugging reasons ... print also ghost values
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_COMPARE_SINGLE0 (X, CROUTINE, CNAME)
REAL (EB), POINTER, DIMENSION (:, :, :) :: X
INTEGER :: II, jj, KK
CHARACTER (4) :: CROUTINE
CHARACTER (6) :: CNAME
 
IF (SCARC_DEBUG >= 2) THEN
   WRITE (9,*) '============================================================='
   WRITE (9,*) '===   COMPARE= ',CNAME,': ROUTINE= ',CROUTINE
   WRITE (9,*) '============================================================='
   IF (TWO_D) THEN
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
   !CALL flush (9)
ENDIF
END SUBROUTINE SCARC_COMPARE_SINGLE0
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! only for debugging reasons ... related to single level
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_COMPARE_SINGLE_LEV (X, CROUTINE, CNAME, ILEVEL)
REAL (EB), POINTER, DIMENSION (:, :, :) :: X
INTEGER :: ILEVEL
INTEGER :: II, jj, KK
CHARACTER (4) :: CROUTINE
CHARACTER (6) :: CNAME
 
SL => S%SLEVEL(ILEVEL)

IF (SCARC_DEBUG >= 2) THEN
   WRITE (9,*) '============================================================='
   WRITE (9,*) '===   COMPARE= ',CNAME,': ROUTINE= ',CROUTINE, ': LEVEL=',ILEVEL
   WRITE (9,*) '============================================================='
   IF (TWO_D) THEN
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
   !CALL flush (9)
ENDIF
 
 
 
END SUBROUTINE SCARC_COMPARE_SINGLE_LEV
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! only for debugging reasons ...
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_COMPARE_SINGLE02(X, CROUTINE, CNAME)
REAL (EB), POINTER, DIMENSION (:, :, :) :: X
INTEGER :: II, jj, KK
CHARACTER (4) :: CROUTINE
CHARACTER (6) :: CNAME
 
IF (SCARC_DEBUG >= 2) THEN
   WRITE (9,*) '============================================================='
   WRITE (9,*) '===   COMPARE= ',CNAME,': ROUTINE= ',CROUTINE
   WRITE (9,*) '============================================================='
   IF (TWO_D) THEN
      WRITE (9,*) 
      WRITE (9,*) '===   COMPARE= ',CNAME,': SCARC_COUNT= ',SCARC_COUNT,'J=0'
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
      WRITE (9,*) '===   COMPARE= ',CNAME,': SCARC_COUNT= ',SCARC_COUNT,'J=1'
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
      WRITE (9,*) '===   COMPARE= ',CNAME,': SCARC_COUNT= ',SCARC_COUNT,'J=2'
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
   !CALL flush (9)
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

SUBROUTINE GET_REV_scrc(MODULE_REV,MODULE_DATE)
INTEGER,INTENT(INOUT) :: MODULE_REV
CHARACTER(255),INTENT(INOUT) :: MODULE_DATE

WRITE(MODULE_DATE,'(A)') scrcrev(INDEX(scrcrev,':')+1:LEN_TRIM(scrcrev)-2)
READ (MODULE_DATE,'(I5)') MODULE_REV
WRITE(MODULE_DATE,'(A)') scrcdate

END SUBROUTINE GET_REV_scrc

 
END MODULE SCARC_SOLVER
 
