MODULE SCRC
 
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Public structures (needed in main, read, divg, dump)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!----------------------------------------------------------------------------------------------------
!!! Public subroutines (initialization, solver, time measurement and revisioning)
!!!----------------------------------------------------------------------------------------------------
PUBLIC SCARC_SETUP                   
PUBLIC SCARC_SOLVER                 
PUBLIC SCARC_TIMINGS               
PUBLIC GET_REV_SCRC               

!!!----------------------------------------------------------------------------------------------------
!!! Public variables   (explanations in declaration part below)
!!! Note: For input parameters in character format corresponding INTEGER type-parameters will
!!! be introduced later to simplify inquiries
!!!----------------------------------------------------------------------------------------------------
PUBLIC SCARC_METHOD   
PUBLIC SCARC_SYSTEM  
PUBLIC SCARC_INITIAL  

PUBLIC SCARC_RESIDUAL  
PUBLIC SCARC_ITERATIONS       
PUBLIC SCARC_CAPPA     
PUBLIC SCARC_ACCURACY        
PUBLIC SCARC_ACCURACY_DIVERGENCE 
PUBLIC SCARC_ACCURACY_RELATIVE 

PUBLIC SCARC_MULTIGRID      
PUBLIC SCARC_MULTIGRID_LEVEL
PUBLIC SCARC_MULTIGRID_CYCLE
PUBLIC SCARC_MULTIGRID_COARSENING 
PUBLIC SCARC_MULTIGRID_INTERPOLATION 
PUBLIC SCARC_MULTIGRID_ITERATIONS       
PUBLIC SCARC_MULTIGRID_ACCURACY      

PUBLIC SCARC_KRYLOV      
PUBLIC SCARC_KRYLOV_ITERATIONS
PUBLIC SCARC_KRYLOV_ACCURACY

PUBLIC SCARC_SMOOTH     
PUBLIC SCARC_SMOOTH_ITERATIONS 
PUBLIC SCARC_SMOOTH_ACCURACY 
PUBLIC SCARC_SMOOTH_OMEGA 

PUBLIC SCARC_PRECON      
PUBLIC SCARC_PRECON_ITERATIONS 
PUBLIC SCARC_PRECON_ACCURACY
PUBLIC SCARC_PRECON_OMEGA 

PUBLIC SCARC_COARSE      
PUBLIC SCARC_COARSE_ITERATIONS  
PUBLIC SCARC_COARSE_ACCURACY  
PUBLIC SCARC_COARSE_OMEGA
PUBLIC SCARC_COARSE_PRECON 
 
PUBLIC SCARC_DEBUG 

!!!----------------------------------------------------------------------------------------------------
!!! corresponding declarations (with default settings)
!!!----------------------------------------------------------------------------------------------------

!!! General definitions
CHARACTER(20) :: SCARC_METHOD    = 'null'                     ! requested solver method (KRYLOV/MULTIGRID)
CHARACTER(20) :: SCARC_SYSTEM    = 'null'                     ! matrix storage technique (BANDED/COMPACT)
CHARACTER(20) :: SCARC_INITIAL   = 'null'                     ! initial solution (currently only default is used)

!!! General iteration parameters
REAL (EB)     :: SCARC_RESIDUAL             = -1.0_EB         ! residual of global selected solver
INTEGER       :: SCARC_ITERATIONS           =  0              ! number of iterations of selected ScaRC solver
REAL (EB)     :: SCARC_CAPPA                =  1.0_EB         ! convergence rate of selected ScarC solver
REAL (EB)     :: SCARC_ACCURACY_DIVERGENCE  = 1.E+6_EB        ! divergence epsilon for all solvers
REAL (EB)     :: SCARC_ACCURACY_RELATIVE    = 1.E-2_EB        ! minimum relative accuracy for all solvers
CHARACTER(20) :: SCARC_ACCURACY             = 'ABSOLUTE'      ! accuracy type (ABSOLUTE/RELATIVE)

!!! Parameters for multigrid-type methods
CHARACTER(20) :: SCARC_MULTIGRID               = 'GEOMETRIC'  ! type of MG-method (GEOMETRIC/ALGEBRAIC)
INTEGER       :: SCARC_MULTIGRID_LEVEL         = -1           ! User defined number of MG-levels (optionally)
CHARACTER(1)  :: SCARC_MULTIGRID_CYCLE         = 'V'          ! Cycling type  (F/V/W)
CHARACTER(20) :: SCARC_MULTIGRID_COARSENING    = 'RS3'        ! Coarsening strategy  (RS3/A1/A2/PMIS/FDS...)
CHARACTER(20) :: SCARC_MULTIGRID_INTERPOLATION = 'DIRECT'     ! Interpolation strategy (DIRECT/RS/STANDARD)
INTEGER       :: SCARC_MULTIGRID_ITERATIONS    = 1000         ! max number of iterations
REAL (EB)     :: SCARC_MULTIGRID_ACCURACY      = 1.E-15_EB    ! requested accuracy for convergence

!!! Parameters for Krylov-type methods
CHARACTER(20) :: SCARC_KRYLOV            = 'CG'               ! type of Krylov-method (CG/BICG)
INTEGER       :: SCARC_KRYLOV_ITERATIONS = 1000               ! max number of iterations 
REAL (EB)     :: SCARC_KRYLOV_ACCURACY   = 1.E-15_EB          ! requested accuracy for convergence

!!! Parameters for smoothing method (used in multigrids-methods)
CHARACTER(20) :: SCARC_SMOOTH            = 'JACOBI'           ! smoother for MG (JACOBI/SSOR/GSTRIX)
INTEGER       :: SCARC_SMOOTH_ITERATIONS = 1000               ! max number of iterations 
REAL (EB)     :: SCARC_SMOOTH_ACCURACY   = 1.E-15_EB          ! requested accuracy for convergence
REAL (EB)     :: SCARC_SMOOTH_OMEGA      = 0.90E+0_EB         ! relaxation parameter 

!!! Parameters for preconditioning method (used in Krylov-methods)
CHARACTER(20) :: SCARC_PRECON            = 'JACOBI'           ! preconditioner for CG/BICG (JACOBI/SSOR/GSTRIX/MG)
INTEGER       :: SCARC_PRECON_ITERATIONS = 1000               ! max number of iterations 
REAL (EB)     :: SCARC_PRECON_ACCURACY   = 1.E-15_EB          ! requested accuracy for convergence
REAL (EB)     :: SCARC_PRECON_OMEGA      = 0.90E+0_EB         ! relaxation parameter 

!!! Parameters for coarse grid method
CHARACTER(20) :: SCARC_COARSE            = 'CG'               ! coarse grid solver (CG/Gaussian elimination)
INTEGER       :: SCARC_COARSE_ITERATIONS = 1000               ! max number of iterations 
REAL (EB)     :: SCARC_COARSE_ACCURACY   = 1.E-15_EB          ! requested accuracy for convergence
REAL (EB)     :: SCARC_COARSE_OMEGA      = 1.5E+0_EB          ! relaxation parameter 
CHARACTER(20) :: SCARC_COARSE_PRECON     = 'SSOR'             ! preconditioner
 
!!! debugging parameters
CHARACTER(20) :: SCARC_DEBUG = 'NONE'                         ! debugging level (NONE/LESS/MEDIUM/MUCH)
CHARACTER(40) :: SCARC_FN                                     ! file name for ScaRC debug messages
INTEGER       :: SCARC_LU                                     ! unit number for ScaRC debug file


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Private structures
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!----------------------------------------------------------------------------------------------------
!!! Global constants
!!!----------------------------------------------------------------------------------------------------
INTEGER, PARAMETER :: NSCARC_DIMENSION_NONE         = -1, &
                      NSCARC_DIMENSION_TWO          =  2, &      ! two-dimensional problem
                      NSCARC_DIMENSION_THREE        =  3         ! three-dimensional problem

INTEGER, PARAMETER :: NSCARC_SCOPE_NONE             = -1, &         
                      NSCARC_SCOPE_MAIN             =  1, &      ! method used as main solver
                      NSCARC_SCOPE_SMOOTH           =  2, &      ! method used as smoother
                      NSCARC_SCOPE_PRECON           =  3, &      ! method used as preconditiner
                      NSCARC_SCOPE_COARSE           =  4         ! method used as coarse grid solver

INTEGER, PARAMETER :: NSCARC_METHOD_NONE            = -1, &
                      NSCARC_METHOD_KRYLOV          =  1, &      ! Krylov   -method as global solver
                      NSCARC_METHOD_MULTIGRID       =  2         ! multigrid-method as global solver

INTEGER, PARAMETER :: NSCARC_KRYLOV_NONE            = -1, &
                      NSCARC_KRYLOV_CG              =  1, &      ! CG   as Krylov solver
                      NSCARC_KRYLOV_BICG            =  2         ! BICG as Krylov solver

INTEGER, PARAMETER :: NSCARC_MULTIGRID_NONE         = -1, &
                      NSCARC_MULTIGRID_GEOMETRIC    =  1, &      ! geometric multigrid
                      NSCARC_MULTIGRID_ALGEBRAIC    =  2         ! algebraic multigrid

INTEGER, PARAMETER :: NSCARC_SYSTEM_NONE            = -1, &
                      NSCARC_SYSTEM_BANDED          =  1, &      ! system in banded storage technique
                      NSCARC_SYSTEM_COMPACT         =  2         ! system in compact  storage technique

INTEGER, PARAMETER :: NSCARC_EXCHANGE_NONE          = -1, &
                      NSCARC_EXCHANGE_INIT          =  1, &      ! initialize communication
                      NSCARC_EXCHANGE_VECTOR        =  2, &      ! matrix-vector communication 
                      NSCARC_EXCHANGE_BDRY          =  3, &      ! vector values along internal boundaries
                      NSCARC_EXCHANGE_MATRIX        =  4, &      ! matrix values along internal boundaries
                      NSCARC_EXCHANGE_MEASURE       =  5         ! measure values along internal boundaries

INTEGER, PARAMETER :: NSCARC_SMOOTH_NONE            = -1, &
                      NSCARC_SMOOTH_JACOBI          =  1, &      ! smoothing by JACOBI-method
                      NSCARC_SMOOTH_SSOR            =  2, &      ! smoothing by SSOR-method
                      NSCARC_SMOOTH_GSTRIX          =  3, &      ! smoothing by GSTRIX-method
                      NSCARC_SMOOTH_FFT             =  4         ! smoothing by FFT-method

INTEGER, PARAMETER :: NSCARC_PRECON_NONE            = -1, &
                      NSCARC_PRECON_JACOBI          =  1, &      ! preconditioning by JACOBI-method
                      NSCARC_PRECON_SSOR            =  2, &      ! preconditioning by SSOR-method
                      NSCARC_PRECON_GSTRIX          =  3, &      ! preconditioning by GSTRIX-method
                      NSCARC_PRECON_FFT             =  4, &      ! preconditioning by FFT-method
                      NSCARC_PRECON_MG              =  5         ! preconditioning by MG-method

INTEGER, PARAMETER :: NSCARC_CYCLE_NONE             = -1, &
                      NSCARC_CYCLE_F                =  0, &      ! F-cycle for mg-method
                      NSCARC_CYCLE_V                =  1, &      ! V-cycle for mg-method
                      NSCARC_CYCLE_W                =  2, &      ! W-cycle for mg-method
                      NSCARC_CYCLE_SETUP            =  3, &      ! initialize cycle counts
                      NSCARC_CYCLE_RESET            =  4, &      ! reset cycle counts
                      NSCARC_CYCLE_PROCEED          =  5, &      ! proceed cycle counts
                      NSCARC_CYCLE_PRESMOOTH        =  6, &      ! presmoothing cycle
                      NSCARC_CYCLE_POSTSMOOTH       =  7, &      ! postsmoothing cycle
                      NSCARC_CYCLE_NEXT             =  8, &      ! perform next cycling loop
                      NSCARC_CYCLE_EXIT             =  9         ! exit cycling loop

INTEGER, PARAMETER :: NSCARC_STATE_PROCEED           =  0, &      ! proceed loop
                      NSCARC_STATE_CONV              =  1, &      ! convergence
                      NSCARC_STATE_DIVG              =  2         ! divergence

INTEGER, PARAMETER :: NSCARC_DEBUG_NONE             = -1, &      ! no debugging requested
                      NSCARC_DEBUG_LESS             =  1, &      ! low    level of debugging requested
                      NSCARC_DEBUG_MEDIUM           =  2, &      ! medium level of debugging requested
                      NSCARC_DEBUG_MUCH             =  3         ! strong level of debugging requested

INTEGER, PARAMETER :: NSCARC_COARSENING_NONE        = -1, &
                      NSCARC_COARSENING_RS3         =  1, &      ! parallel Ruge-Stüben 
                      NSCARC_COARSENING_A1          =  2, &      ! aggressive 1 (path=1, length=2)
                      NSCARC_COARSENING_A2          =  3, &      ! aggressive 2 (path=2, length=2)
                      NSCARC_COARSENING_PMIS        =  4, &      ! PMIS 
                      NSCARC_COARSENING_FDSRS3      =  5, &      ! FDSRS3 : ScaRC-FDS variant similar to RS3
                      NSCARC_COARSENING_FDSA1       =  6, &      ! FDSA1  : ScaRC-FDS variant similar to A1
                      NSCARC_COARSENING_FDSA2       =  7, &      ! FDSA2  : ScaRC-FDS variant similar to A2
                      NSCARC_COARSENING_FDSPMIS     =  8         ! FDSPMIS: ScaRC-FDS variant similar to PMIS

INTEGER, PARAMETER :: NSCARC_COARSE_NONE            = -1, &
                      NSCARC_COARSE_CG              =  1, &      ! coarse grid solution by cg-method
                      NSCARC_COARSE_GE              =  2         ! coarse grid solution by Gaussian elimination

INTEGER, PARAMETER :: NSCARC_VECTOR_NONE            = -1, &
                      NSCARC_VECTOR_X               =  1, &      ! selection parameter for vector X
                      NSCARC_VECTOR_F               =  2, &      ! selection parameter for vector F
                      NSCARC_VECTOR_Y               =  3, &      ! selection parameter for vector Y
                      NSCARC_VECTOR_G               =  4, &      ! selection parameter for vector G
                      NSCARC_VECTOR_W               =  5, &      ! selection parameter for vector R
                      NSCARC_VECTOR_D               =  6, &      ! selection parameter for vector D
                      NSCARC_VECTOR_Z               =  7, &      ! selection parameter for vector Z
                      NSCARC_VECTOR_X2              =  8, &      ! selection parameter for vector X2
                      NSCARC_VECTOR_D2              =  9, &      ! selection parameter for vector D2
                      NSCARC_VECTOR_W2              = 10, &      ! selection parameter for vector R2
                      NSCARC_VECTOR_Y2              = 11         ! selection parameter for vector Y2

INTEGER, PARAMETER :: NSCARC_ACCURACY_NONE          = -1, &
                      NSCARC_ACCURACY_ABSOLUTE      =  1, &      ! absolute accuracy must be reached
                      NSCARC_ACCURACY_RELATIVE      =  2         ! relative accuracy must be reached

INTEGER, PARAMETER :: NSCARC_CELLTYPE_NONE          =  0, &
                      NSCARC_CELLTYPE_COARSE        =  1, &      ! coarse-grid cell
                      NSCARC_CELLTYPE_FINE          = -2, &      ! fine-grid cell
                      NSCARC_CELLTYPE_SFINE         = -2, &      ! strongly coupled fine-grid cell
                      NSCARC_CELLTYPE_WFINE         = -3         ! weakly   coupled fine-grid cell

INTEGER, PARAMETER :: NSCARC_TIME_NONE              = -1, &
                      NSCARC_TIME_TOTAL             =  1, &      ! time for complete ScaRC part of FDS
                      NSCARC_TIME_SETUP             =  2, &      ! time for setup phase
                      NSCARC_TIME_SOLVER            =  3, &      ! time for ScaRC solver 
                      NSCARC_TIME_KRYLOV            =  4, &      ! time for Krylov solver
                      NSCARC_TIME_MULTIGRID         =  5, &      ! time for multigrid solver
                      NSCARC_TIME_PRECON            =  6, &      ! time for preconditioner
                      NSCARC_TIME_SMOOTH            =  7, &      ! time for smoother
                      NSCARC_TIME_COARSE            =  8, &      ! time for coarse grid solver
                      NSCARC_TIME_MATVEC            =  9, &      ! time for matrix-vector product
                      NSCARC_TIME_SCALPROD          = 10, &      ! time for scalar product
                      NSCARC_TIME_L2NORM            = 11, &      ! time for l2norm
                      NSCARC_TIME_EXCH_INIT         = 12, &      ! time for exchange initialization
                      NSCARC_TIME_EXCH_VECTOR       = 13, &      ! time for exchange of internal boundary
                      NSCARC_TIME_EXCH_MATRIX       = 14, &      ! time for exchange of internal matrix
                      NSCARC_TIME_EXCH_MEASURE      = 15         ! time for exchange of internal measure

INTEGER, PARAMETER :: NSCARC_LEVEL_NONE             = -1, &      ! no predefined initial solution used
                      NSCARC_LEVEL_MIN              =  0, &      ! minimum multigrid level 
                      NSCARC_LEVEL_MAX              =  8         ! maximum multigrid level 

INTEGER, PARAMETER :: NSCARC_INITIAL_NONE           = -1         ! another initial function ?
 
INTEGER, PARAMETER :: NSCARC_DUMMY                  = -1         ! dummy variable (needed at several places)


!!!----------------------------------------------------------------------------------------------------
!!! Global variables 
!!!----------------------------------------------------------------------------------------------------
!!! use integer types for the user defined input data (based on SCARC_TYPE_... variables)
INTEGER :: TYPE_DIMENSION  = NSCARC_DIMENSION_NONE
INTEGER :: TYPE_SCOPE      = NSCARC_SCOPE_NONE
INTEGER :: TYPE_METHOD     = NSCARC_METHOD_NONE
INTEGER :: TYPE_KRYLOV     = NSCARC_KRYLOV_NONE
INTEGER :: TYPE_MULTIGRID  = NSCARC_MULTIGRID_NONE
INTEGER :: TYPE_SYSTEM     = NSCARC_SYSTEM_NONE
INTEGER :: TYPE_ACCURACY   = NSCARC_ACCURACY_NONE
INTEGER :: TYPE_SMOOTH     = NSCARC_SMOOTH_NONE
INTEGER :: TYPE_PRECON     = NSCARC_PRECON_NONE
INTEGER :: TYPE_CYCLE      = NSCARC_CYCLE_NONE
INTEGER :: TYPE_COARSENING = NSCARC_COARSENING_NONE
INTEGER :: TYPE_COARSE     = NSCARC_COARSE_NONE
INTEGER :: TYPE_DEBUG      = NSCARC_DEBUG_NONE   
INTEGER :: TYPE_INITIAL    = NSCARC_INITIAL_NONE
INTEGER :: TYPE_EXCHANGE   = NSCARC_EXCHANGE_NONE
INTEGER :: TYPE_VECTOR     = NSCARC_VECTOR_NONE

!!! range of meshes which must be processed for MYID
INTEGER :: NMESHES_MIN, NMESHES_MAX                 
 
!!! total, minimum and maximum number of multigrid levels
INTEGER :: NLEVEL, NLEVEL_MAX, NLEVEL_MIN                              

!!! additional arrays for data exchange
INTEGER :: NREQ_SCARC, N_EXCHANGES, TAG_SCARC, SNODE, RNODE, STATUS2_SCARC(MPI_STATUS_SIZE)
INTEGER, ALLOCATABLE, DIMENSION (:)    :: REQ_SCARC

!!! time measurements with ScaRC
INTEGER, PARAMETER :: N_TIMERS_SCARC=18         
REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: TUSED_SCARC

!!! auxiliary variables for global and local scalproducts and number of cells
REAL(EB), ALLOCATABLE, DIMENSION (:) :: SP_LOCAL
INTEGER,  ALLOCATABLE, DIMENSION (:) :: NC_GLOBAL, NC_LOCAL

!!! number of couplings in given matrix stencil and pointer to indices in matrix stencil on given level 
INTEGER :: ID, ILX, ILY, ILZ, IUX, IUY, IUZ

!!! Index numbers of different vector types (used in different ScaRC-solvers)
INTEGER :: VEC_NONE, VEC_X, VEC_F, VEC_Y, VEC_G, VEC_W, VEC_D, VEC_Z, VEC_X2, VEC_D2, VEC_W2, VEC_Y2

 
!!! Private type declarations

!!!----------------------------------------------------------------------------------------------------
!!! SCARC type for mesh information
!!! local numbers of cells (also per direction), wall cells, matrix entries
!!! neighbourship structures, mesh coordinates and step sizes
!!!----------------------------------------------------------------------------------------------------
TYPE SCARC_MESH_TYPE

INTEGER :: IBAR, JBAR, KBAR
INTEGER :: N_CELLS, N_CELLSG, N_CELLS_COARSE, N_CELLS_FINE, N_EXTERNAL_WALL_CELLS
INTEGER, POINTER, DIMENSION (:)    :: PRESSURE_BC_INDEX, ADJACENT_CELL
INTEGER, POINTER, DIMENSION (:, :) :: IJKW
INTEGER   :: SUBDIVISION(3,-3:3)
REAL (EB) :: DX,    DY,    DZ
REAL (EB) :: DXI,   DYI,   DZI
REAL (EB) :: DXI2,  DYI2,  DZI2,  DI2(3)
REAL (EB) :: DXMIN, DYMIN, DZMIN

END TYPE SCARC_MESH_TYPE


!!!----------------------------------------------------------------------------------------------------
!!! OSCARC type for mesh information on other mesh
!!! local numbers of cells (also per direction), wall cells, matrix entries
!!! neighbourship structures 
!!!----------------------------------------------------------------------------------------------------
TYPE OSCARC_MESH_TYPE

INTEGER :: IBAR, JBAR, KBAR
INTEGER :: N_CELLS, N_CELLSG, N_EXTERNAL_WALL_CELLS
INTEGER, POINTER, DIMENSION (:, :) :: IJKW

END TYPE OSCARC_MESH_TYPE


!!!----------------------------------------------------------------------------------------------------
!!! SCARC type for matrix stencil information
!!! offset in stencil, measure of cells, coarse and fine cells, couplings
!!!----------------------------------------------------------------------------------------------------
TYPE SCARC_PRECON_TYPE

REAL (EB), POINTER, DIMENSION (:)       :: MDX, LDX, UDX, MDY, MDZ, DWORK, PERIOD
REAL (EB), POINTER, DIMENSION (:, :, :) :: FFT

END TYPE SCARC_PRECON_TYPE


!!!----------------------------------------------------------------------------------------------------
!!! SCARC types for different matrix storage technique 
!!!   - NX, NY, NZ : correspond to IBAR, JBAR, KBAR
!!!   - NC, NCE    : number of cells and cells_extended (including ghost cells)
!!!   - NA, NP, NR : number of matrix entries in A, P and R
!!!   - NCPL       : number of couplings in matrix stencil on finest level (2D: 5, 3D: 7)
!!!   - NCC, NCF   : number of coarse and fine cells (only used in AMG)
!!!----------------------------------------------------------------------------------------------------
!!! Banded storage technique 
TYPE SCARC_SYSTEM_BANDED_TYPE

INTEGER  , POINTER, DIMENSION (:, :, :) :: CELLTYPE
REAL (EB), POINTER, DIMENSION (:, :)    :: A
REAL (EB), POINTER, DIMENSION (:, :, :) :: X , F , D , Y , G , W, Z
REAL (EB), POINTER, DIMENSION (:, :, :) :: X2, F2, D2, Y2, G2, W2
INTEGER   :: NX, NY, NZ
INTEGER   :: NC, NCE, NCG
INTEGER   :: NA, NCPL

END TYPE SCARC_SYSTEM_BANDED_TYPE


!!! Compact storage technique 
TYPE SCARC_SYSTEM_COMPACT_TYPE

INTEGER,   POINTER, DIMENSION (:) :: CELLTYPE
REAL(EB),  POINTER, DIMENSION (:) :: MEASURE
INTEGER,   POINTER, DIMENSION (:) :: A_ROW, A_COL       ! row and column pointers for system matrix A
INTEGER,   POINTER, DIMENSION (:) :: P_ROW, P_COL       ! row and column pointers for prolongation matrix P
INTEGER,   POINTER, DIMENSION (:) :: R_ROW, R_COL       ! row and column pointers for restriction matrix A
REAL (EB), POINTER, DIMENSION (:) :: A , P , R
REAL (EB), POINTER, DIMENSION (:) :: X , F , D , Y , G , W, Z
REAL (EB), POINTER, DIMENSION (:) :: X2, F2, D2, Y2, G2, W2
INTEGER   :: NX, NY, NZ
INTEGER   :: NC, NCE, NCG
INTEGER   :: NA, NCPL
INTEGER   :: NCC, NCF
INTEGER   :: NP, NR

END TYPE SCARC_SYSTEM_COMPACT_TYPE


!!!----------------------------------------------------------------------------------------------------
!!! General ScaRC type with pointers to different structures on all grid levels
!!!----------------------------------------------------------------------------------------------------
TYPE SCARC_TYPE

INTEGER :: CYCLE_COUNT(2, NSCARC_LEVEL_MAX) = 0

TYPE (SCARC_MESH_TYPE)  , POINTER, DIMENSION(:) :: MESHES
TYPE (SCARC_PRECON_TYPE), POINTER, DIMENSION(:) :: PRECON

TYPE (SCARC_SYSTEM_BANDED_TYPE) , POINTER, DIMENSION(:) :: BANDED
TYPE (SCARC_SYSTEM_COMPACT_TYPE), POINTER, DIMENSION(:) :: COMPACT

TYPE (OSCARC_TYPE), POINTER, DIMENSION(:) :: OSCARC

END TYPE SCARC_TYPE


!!!----------------------------------------------------------------------------------------------------
!!! General OSCARC type on other mesh with mesh and exchange structures
!!!----------------------------------------------------------------------------------------------------
TYPE OSCARC_TYPE

REAL (EB), POINTER, DIMENSION (:) :: SEND_BUF, RECV_BUF
INTEGER :: NICMAX_R=0, NICMAX_S=0
INTEGER :: I_MIN_R=-10,I_MAX_R=-10,J_MIN_R=-10,J_MAX_R=-10,K_MIN_R=-10,K_MAX_R=-10,NIC_R=0, &
           I_MIN_S=-10,I_MAX_S=-10,J_MIN_S=-10,J_MAX_S=-10,K_MIN_S=-10,K_MAX_S=-10,NIC_S=0

TYPE (OSCARC_MESH_TYPE), POINTER, DIMENSION(:) :: MESHES

END TYPE OSCARC_TYPE

!!! globally used types
TYPE ( SCARC_TYPE), SAVE, DIMENSION(:), ALLOCATABLE, TARGET ::  SCARC
TYPE (OSCARC_TYPE), SAVE, DIMENSION(:), ALLOCATABLE, TARGET :: OSCARC


CONTAINS
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! SCARC_SETUP : Initialize ScaRC structures 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP

INTEGER :: IERR
REAL(EB):: TNOW_SETUP

TNOW_SETUP = SECOND()
IERR = 0

!!!----------------------------------------------------------------------------------------------------
!!! Allocate general SCARC structure and time measurement array 
!!!----------------------------------------------------------------------------------------------------
ALLOCATE (SCARC(NMESHES), STAT=IERR)
CALL CHKMEMERR ('SCARC_SETUP', 'SCARC', IERR)
 
ALLOCATE(TUSED_SCARC(0:N_TIMERS_SCARC,NMESHES),STAT=IERR)
CALL ChkMemErr('SCARC_SETUP','TUSED_SCARC',IERR)

TUSED_SCARC = 0._EB
TUSED_SCARC(NSCARC_TIME_TOTAL,:) = SECOND()


!!!----------------------------------------------------------------------------------------------------
!!! Setup different components of ScaRC 
!!!----------------------------------------------------------------------------------------------------
CALL SCARC_SETUP_DIMENSION            ! define dimension of underlying problem
CALL SCARC_SETUP_TYPES                ! store types of input parameters (corresponding to FDS-file)
CALL SCARC_SETUP_DEBUGGING            ! open debug file if requested
CALL SCARC_SETUP_PROCESSES            ! determine set of meshes which must be processed on MYID
CALL SCARC_SETUP_LEVELS               ! define number of necessary grid levels 
CALL SCARC_SETUP_MESHES               ! set mesh information
CALL SCARC_SETUP_WALLCELLS            ! set wall cell information
CALL SCARC_SETUP_EXCHANGE             ! set information for data exchange
CALL SCARC_SETUP_SYSTEM               ! assemble system matrix with boundary conditions and solver vectors
CALL SCARC_SETUP_COARSENING           ! perform coarsening on different grid levels if requested (AMG only)
CALL SCARC_SETUP_VECTORS              ! allocate solution and auxiliary vectors on all needed grid levels
CALL SCARC_SETUP_GLOBALS              ! define some global variables
CALL SCARC_SETUP_NEIGHBORS            ! compute information about abutting neighbors on coarser levels

TUSED_SCARC(NSCARC_TIME_SETUP   ,:)=TUSED_SCARC(NSCARC_TIME_SETUP   ,:)+SECOND()-TNOW_SETUP
TUSED_SCARC(NSCARC_TIME_TOTAL,:)=TUSED_SCARC(NSCARC_TIME_TOTAL,:)+SECOND()-TNOW_SETUP
END SUBROUTINE SCARC_SETUP


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Define dimension of problem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_DIMENSION
IF (TWO_D) THEN
   TYPE_DIMENSION = NSCARC_DIMENSION_TWO
ELSE
   TYPE_DIMENSION = NSCARC_DIMENSION_THREE
ENDIF
END SUBROUTINE SCARC_SETUP_DIMENSION


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Determine types of input parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_TYPES
CHARACTER(300):: CMESSAGE

!!!----------------------------------------------------------------------------------------------------
!!! set type of global solver
!!!----------------------------------------------------------------------------------------------------
SELECT CASE (TRIM(SCARC_METHOD))

   CASE ('KRYLOV')

      TYPE_METHOD = NSCARC_METHOD_KRYLOV

      !!! set type of Krylov-method (CG/BICG)
      SELECT CASE (TRIM(SCARC_KRYLOV))
         CASE ('CG')
            TYPE_KRYLOV = NSCARC_KRYLOV_CG
         CASE ('BICG')
            TYPE_KRYLOV = NSCARC_KRYLOV_BICG
         CASE DEFAULT
            WRITE(CMESSAGE,1002) 'Krylov-solver',TRIM(SCARC_KRYLOV),'ScaRC','CG','BICG'
            CALL SCARC_SHUTDOWN(CMESSAGE)
      END SELECT 

      !!! set type of preconditioner (JACOBI/SSOR/GSTRIX/MG)
      SELECT CASE (TRIM(SCARC_PRECON))
         CASE ('JACOBI')
            TYPE_PRECON = NSCARC_PRECON_JACOBI
         CASE ('SSOR')
            TYPE_PRECON = NSCARC_PRECON_SSOR
         CASE ('GSTRIX')
            TYPE_PRECON = NSCARC_PRECON_GSTRIX
         CASE ('MG')
            TYPE_PRECON = NSCARC_PRECON_MG
         CASE ('FFT')
            TYPE_PRECON = NSCARC_PRECON_FFT
         CASE DEFAULT
            WRITE(CMESSAGE,1005) 'preconditioner',TRIM(SCARC_PRECON),&
                                 'Krylov-method','JACOBI','SSOR','GSTRIX','MG','FFT'
            CALL SCARC_SHUTDOWN(CMESSAGE)
      END SELECT

   CASE ('MULTIGRID')

      TYPE_METHOD = NSCARC_METHOD_MULTIGRID

      !!! set type of multigrid method (GEOMETRIC/ALGEBRAIC)
      SELECT CASE (TRIM(SCARC_MULTIGRID))
         CASE ('GEOMETRIC')
            TYPE_MULTIGRID = NSCARC_MULTIGRID_GEOMETRIC
         CASE ('ALGEBRAIC')
            TYPE_MULTIGRID = NSCARC_MULTIGRID_ALGEBRAIC
         CASE DEFAULT
            WRITE(CMESSAGE,1002) 'multigrid',TRIM(SCARC_MULTIGRID),&
                                 'ScaRC','GEOMETRIC','ALGEBRAIC'
            CALL SCARC_SHUTDOWN(CMESSAGE)
      END SELECT 

      !!! set type of smoother (JACOBI/SSOR/GSTRIX)
      SELECT CASE (TRIM(SCARC_SMOOTH))                          ! use same parameters as for preconditioner
         CASE ('JACOBI')
            TYPE_PRECON = NSCARC_PRECON_JACOBI
         CASE ('SSOR')
            TYPE_PRECON = NSCARC_PRECON_SSOR
         CASE ('GSTRIX')
            TYPE_PRECON = NSCARC_PRECON_GSTRIX
         CASE DEFAULT
            WRITE(CMESSAGE,1003) 'smoother',TRIM(SCARC_SMOOTH),'multigrid','JACOBI','SSOR','GSTRIX'
            CALL SCARC_SHUTDOWN(CMESSAGE)
      END SELECT

   CASE DEFAULT
      WRITE(CMESSAGE,1002) 'solver',TRIM(SCARC_METHOD),'ScaRC','KRYLOV','MULTIGRID'
      CALL SCARC_SHUTDOWN(CMESSAGE)

END SELECT 

!!!----------------------------------------------------------------------------------------------------
!!! if a multigrid solver is used (either as main solver or as preconditioner)
!!! set types for multigrid, coarse grid solver and cycling pattern
!!!----------------------------------------------------------------------------------------------------
IF (TYPE_METHOD == NSCARC_METHOD_MULTIGRID .OR. TYPE_PRECON == NSCARC_PRECON_MG) THEN

   !!! set type of multigrid (GEOMETRIC/ALGEBRAIC with corresponding coarsening strategy)
   SELECT CASE (TRIM(SCARC_MULTIGRID))

      CASE ('GEOMETRIC')
         TYPE_MULTIGRID = NSCARC_MULTIGRID_GEOMETRIC

      CASE ('ALGEBRAIC')
         TYPE_MULTIGRID = NSCARC_MULTIGRID_ALGEBRAIC

         !!! set type of coarsening strategy (STANDARD/AGGRESSIVE)
         SELECT CASE (TRIM(SCARC_MULTIGRID_COARSENING))
            CASE ('RS3')
               TYPE_COARSENING = NSCARC_COARSENING_RS3
            CASE ('A1')
               TYPE_COARSENING = NSCARC_COARSENING_A1
            CASE ('A2')
               TYPE_COARSENING = NSCARC_COARSENING_A2
            CASE ('PMIS')
               TYPE_COARSENING = NSCARC_COARSENING_PMIS
            CASE ('FDSRS3')
               TYPE_COARSENING = NSCARC_COARSENING_FDSRS3
            CASE ('FDSA1')
               TYPE_COARSENING = NSCARC_COARSENING_FDSA1
            CASE ('FDSA2')
               TYPE_COARSENING = NSCARC_COARSENING_FDSA2
            CASE ('FDSPMIS')
               TYPE_COARSENING = NSCARC_COARSENING_FDSPMIS
            CASE DEFAULT
               WRITE(CMESSAGE,1005) 'coarsening',TRIM(SCARC_MULTIGRID_COARSENING),&
                                    'algebraic multigrid','RS3'   ,'A1'   ,'A2'   ,'PMIS',&
                                                          'FDSRS3','FDSA1','FDSA1','FDSPMIS'
               CALL SCARC_SHUTDOWN(CMESSAGE)
         END SELECT

      CASE DEFAULT
         WRITE(CMESSAGE,1002) 'multigrid',TRIM(SCARC_MULTIGRID),&
                              'Krylov-preconditioner','GEOMETRIC','ALGEBRAIC'
         CALL SCARC_SHUTDOWN(CMESSAGE)
   END SELECT

   !!! set type of coarse grid solver (CG/GE)
   SELECT CASE (TRIM(SCARC_COARSE))
      CASE ('CG')
         TYPE_COARSE = NSCARC_COARSE_CG
         TYPE_KRYLOV = NSCARC_KRYLOV_CG
      CASE ('GE')
         TYPE_COARSE = NSCARC_COARSE_GE
      CASE DEFAULT
         WRITE(CMESSAGE,1002) 'coarse grid solver',TRIM(SCARC_COARSE),'multigrid','CG','GE'
         CALL SCARC_SHUTDOWN(CMESSAGE)
   END SELECT

   !!! set type of cycling pattern (F/V/W)
   SELECT CASE (TRIM(SCARC_MULTIGRID_CYCLE))
      CASE ('F')
         TYPE_CYCLE = NSCARC_CYCLE_F
      CASE ('V')
         TYPE_CYCLE = NSCARC_CYCLE_V
      CASE ('W')
         TYPE_CYCLE = NSCARC_CYCLE_W
      CASE DEFAULT
         WRITE(CMESSAGE,1003) 'cycling ',TRIM(SCARC_COARSE),'multigrid','F','V','W'
         CALL SCARC_SHUTDOWN(CMESSAGE)
   END SELECT

ENDIF

!!!----------------------------------------------------------------------------------------------------
!!! set storage type (BANDED/COMPACT)
!!!----------------------------------------------------------------------------------------------------
SELECT CASE (TRIM(SCARC_SYSTEM))
   CASE ('BANDED')
      TYPE_SYSTEM = NSCARC_SYSTEM_BANDED
   CASE ('COMPACT')
      TYPE_SYSTEM = NSCARC_SYSTEM_COMPACT
   CASE DEFAULT
      WRITE(CMESSAGE,1002) 'storage',TRIM(SCARC_SYSTEM),'ScaRC','BANDED','COMPACT'
      CALL SCARC_SHUTDOWN(CMESSAGE)
END SELECT


!!!----------------------------------------------------------------------------------------------------
!!! set type of accuracy (ABSOLUTE/RELATIVE)
!!!----------------------------------------------------------------------------------------------------
SELECT CASE (TRIM(SCARC_ACCURACY))
   CASE ('ABSOLUTE')
      TYPE_ACCURACY = NSCARC_ACCURACY_ABSOLUTE
   CASE ('RELATIVE')
      TYPE_ACCURACY = NSCARC_ACCURACY_RELATIVE
   CASE DEFAULT
      WRITE(CMESSAGE,1002) 'accuracy',TRIM(SCARC_ACCURACY),'ScaRC','ABSOLUTE','RELATIVE'
      CALL SCARC_SHUTDOWN(CMESSAGE)
END SELECT

!!!----------------------------------------------------------------------------------------------------
!!! set level of debugging (NONE/LESS/MEDIUM/MUCH)
!!!----------------------------------------------------------------------------------------------------
SELECT CASE (TRIM(SCARC_DEBUG))
   CASE ('NONE')
      TYPE_DEBUG = NSCARC_DEBUG_NONE
   CASE ('LESS')
      TYPE_DEBUG = NSCARC_DEBUG_LESS
   CASE ('MEDIUM')
      TYPE_DEBUG = NSCARC_DEBUG_MEDIUM
   CASE ('MUCH')
      TYPE_DEBUG = NSCARC_DEBUG_MUCH
   CASE DEFAULT
      WRITE(CMESSAGE,1004) 'debugging',TRIM(SCARC_DEBUG),'ScaRC','NONE','LESS','MEDIUM','MUCH'
      CALL SCARC_SHUTDOWN(CMESSAGE)
END SELECT

!!!----------------------------------------------------------------------------------------------------
!!! set type of initial solution (not yet used, may be used to define own initial vector)
!!!----------------------------------------------------------------------------------------------------
SELECT CASE (TRIM(SCARC_INITIAL))
   CASE ('null')
      TYPE_INITIAL = NSCARC_INITIAL_NONE
   CASE DEFAULT
      TYPE_INITIAL = NSCARC_INITIAL_NONE
END SELECT

1002 FORMAT ('ScaRC: Wrong ',A,' type ',A,' for ',A,/, &
             'ScaRC: Possible choices are ',A,'\/',A,/,&
             'ScaRC: Aborting program ...')
1003 FORMAT ('ScaRC: Wrong ',A,' type ',A,' for ',A,/, &
             'ScaRC: Possible choices are ',A,'\/',A,'\/',A,/,&
             'ScaRC: Aborting program ...')
1004 FORMAT ('ScaRC: Wrong ',A,' type ',A,' for ',A,/, &
             'ScaRC: Possible choices are ',A,'\/',A,'\/',A,'\/',A,/,&
             'ScaRC: Aborting program ...')
1005 FORMAT ('ScaRC: Wrong ',A,' type ',A,' for ',A,/, &
             'ScaRC: Possible choices are ',A,'\/',A,'\/',A,'\/',A,'\/',A,/,&
             'ScaRC: Aborting program ...')
!1006 FORMAT ('ScaRC: Wrong ',A,' type ',A,' for ',A,/, &
!             'ScaRC: Possible choices are ',A,'\/',A,'\/',A,'\/',A,'\/',A,'\/',A,/,&
!             'ScaRC: Aborting program ...')
!1007 FORMAT ('ScaRC: Wrong ',A,' type ',A,' for ',A,/, &
!             'ScaRC: Possible choices are ',A,'\/',A,'\/',A,'\/',A,'\/',A,'\/',A,'\/',A,/,&
!             'ScaRC: Aborting program ...')
!1008 FORMAT ('ScaRC: Wrong ',A,' type ',A,' for ',A,/, &
!             'ScaRC: Possible choices are ',A,'\/',A,'\/',A,'\/',A,'\/',A,'\/',A,'\/',A,'\/',A,/,&
!             'ScaRC: Aborting program ...')
END SUBROUTINE SCARC_SETUP_TYPES


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Setup debug file if requested
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_DEBUGGING
INTEGER:: NM

IF (TYPE_DEBUG /= NSCARC_DEBUG_NONE) THEN
   IF (USE_MPI) THEN
      DO NM=1,NMESHES
         IF (PROCESS(NM)/=MYID) CYCLE
         WRITE (SCARC_FN, '(A,A,i3.3)') TRIM(CHID),'.scarc',MYID+1
         SCARC_LU = GET_FILE_NUMBER()
         OPEN (SCARC_LU, FILE=SCARC_FN)
      ENDDO
   ELSE
      WRITE (SCARC_FN, '(A,A,i3.3)') TRIM(CHID),'.scarc',MYID+1
      SCARC_LU = GET_FILE_NUMBER()
      OPEN (SCARC_LU, FILE=SCARC_FN)
   ENDIF
ENDIF

END SUBROUTINE SCARC_SETUP_DEBUGGING


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Serial or parallel version ? Decide which meshes must be processed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_PROCESSES
IF (USE_MPI) THEN
   NMESHES_MIN = MYID+1
   NMESHES_MAX = MYID+1
ELSE
   NMESHES_MIN = 1
   NMESHES_MAX = NMESHES
ENDIF
END SUBROUTINE SCARC_SETUP_PROCESSES

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Determine number of grid levels  (1 for CG/BICG-method, NLEVEL for MG-method)
!!! Note: NLEVEL_MIN corresponds to finest grid resolution, NLEVEL_MAX to coarsest resolution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_LEVELS
INTEGER :: KLEVEL(3), KLEVEL_MIN, NM

SELECT CASE (TYPE_MULTIGRID)

   !!!------------------------------------------------------------------------------------------------------
   !!! predefined hierarchy of levels in case of geometric multigrid-method
   !!!------------------------------------------------------------------------------------------------------
   CASE (NSCARC_MULTIGRID_GEOMETRIC, NSCARC_MULTIGRID_ALGEBRAIC)

      NLEVEL = NSCARC_LEVEL_MAX
      KLEVEL = NSCARC_LEVEL_MAX

      DO NM=1,NMESHES

         KLEVEL(1)=SCARC_GET_MAXLEVEL(MESHES(NM)%IBAR,1)
         IF (.NOT.TWO_D) KLEVEL(2)=SCARC_GET_MAXLEVEL(MESHES(NM)%JBAR,2)
         KLEVEL(3)=SCARC_GET_MAXLEVEL(MESHES(NM)%KBAR,3)

         KLEVEL_MIN = MINVAL(KLEVEL)
         IF (KLEVEL_MIN<NLEVEL) NLEVEL=KLEVEL_MIN

      ENDDO 
      NLEVEL_MIN  = 1
      IF (SCARC_MULTIGRID_LEVEL /= -1) THEN
         NLEVEL_MAX  = NLEVEL_MIN + SCARC_MULTIGRID_LEVEL - 1
      ELSE
         NLEVEL_MAX  = NLEVEL
      ENDIF

   !!!------------------------------------------------------------------------------------------------------
   !!! first, only finest level is set, further levels are defined during coarsening process
   !!!------------------------------------------------------------------------------------------------------
   !CASE (NSCARC_MULTIGRID_ALGEBRAIC)
!
!      NLEVEL     = NSCARC_LEVEL_MAX
!      NLEVEL_MAX = NSCARC_LEVEL_MAX
!      NLEVEL_MIN = 1

   !!!------------------------------------------------------------------------------------------------------
   !!! no multigrid-hierachy needed in case of a pure Krylov-method: use only one level
   !!!------------------------------------------------------------------------------------------------------
   CASE DEFAULT

      NLEVEL     = 1
      NLEVEL_MIN = 1
      NLEVEL_MAX = 1

END SELECT

END SUBROUTINE SCARC_SETUP_LEVELS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Determine maximum number of possible levels on direction IOR0 of mesh NM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER FUNCTION SCARC_GET_MAXLEVEL(NC, IOR0)
INTEGER, INTENT(IN) :: NC, IOR0
INTEGER :: NC0, NL

!!! In case of the GMG-method, NC must be divisable by 2 at least one time 
IF (MOD(NC,2)/=0 .AND. TYPE_MULTIGRID == NSCARC_MULTIGRID_GEOMETRIC) THEN
   SELECT CASE (IOR0)
      CASE (1)
         WRITE(*,1000) IBAR, NC
      CASE (2)
         WRITE(*,1000) JBAR, NC
      CASE (3)
         WRITE(*,1000) KBAR, NC
   END SELECT
   STOP
ENDIF

!!! divide by 2 as often as possible or till user defined max-level is reached
NC0=NC
DO NL=1,NSCARC_LEVEL_MAX
   NC0=NC0/2
   IF (MOD(NC0,2)/=0) EXIT                  ! NC no longer divisable by two
   IF (NL==SCARC_MULTIGRID_LEVEL) EXIT      ! max number of levels defined by user
   IF (NC0==1) EXIT                         ! NC is power of two, minimum has been reached
ENDDO

SCARC_GET_MAXLEVEL=NL
RETURN
1000 FORMAT(A,'=',I3,' must be divisable by two for ScaRC-Multigrid!')
END FUNCTION SCARC_GET_MAXLEVEL


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Setup geometry information for mesh NM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_MESHES
INTEGER :: IERR, NL, NM
INTEGER :: IBAR0, JBAR0, KBAR0
TYPE (MESH_TYPE)      , POINTER :: M
TYPE (SCARC_TYPE)     , POINTER :: S
TYPE (SCARC_MESH_TYPE), POINTER :: SM

IERR=0

MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX

   M => MESHES(NM)
   S => SCARC(NM)
   
   !!!
   !!!  define hierarchy of meshes depending on the chosen method
   !!!
   ALLOCATE (S%MESHES(NLEVEL_MIN:NLEVEL_MAX), STAT=IERR)
   CALL CHKMEMERR ('SCARC_SETUP_MESH', 'MESH', IERR)
   
   IBAR0=MESHES(NM)%IBAR
   JBAR0=MESHES(NM)%JBAR
   KBAR0=MESHES(NM)%KBAR
   
   LEVEL_LOOP: DO NL = NLEVEL_MIN, NLEVEL_MAX
   
      !!! let SM point to SCARC(NM)%MESHES(NL)
      SM => S%MESHES(NL)
      
      !!!-------------------------------------------------------------------------------------------------
      !!! numbers of cells in x-, y- and z-direction for level 'NL'
      !!!-------------------------------------------------------------------------------------------------
      SM%IBAR = IBAR0
      SM%JBAR = JBAR0
      SM%KBAR = KBAR0
      
      !!!-------------------------------------------------------------------------------------------------
      !!! step widths in x-, y- and z-direction for level 'NL'
      !!!-------------------------------------------------------------------------------------------------
      SM%DX = (M%XF-M%XS)/REAL(SM%IBAR,EB)
      SM%DY = (M%YF-M%YS)/REAL(SM%JBAR,EB)
      SM%DZ = (M%ZF-M%ZS)/REAL(SM%KBAR,EB)
      
      SM%DXI = 1.0_EB/SM%DX
      SM%DYI = 1.0_EB/SM%DY
      SM%DZI = 1.0_EB/SM%DZ
      
      SM%DXI2 = 1.0_EB/(SM%DX)**2
      SM%DYI2 = 1.0_EB/(SM%DY)**2
      SM%DZI2 = 1.0_EB/(SM%DZ)**2
    
      SM%DI2(1) = SM%DXI2
      SM%DI2(2) = SM%DYI2
      SM%DI2(3) = SM%DZI2

      !!!-------------------------------------------------------------------------------------------------
      !!! Get global number of grid cells (internal and including ghost cells)
      !!!-------------------------------------------------------------------------------------------------
      SM%N_CELLS  =  SM%IBAR * SM%JBAR * SM%KBAR 
      SELECT CASE (TYPE_DIMENSION)
         CASE (NSCARC_DIMENSION_TWO)
            SM%N_CELLSG = (SM%IBAR+2) * (SM%KBAR+2) 
         CASE (NSCARC_DIMENSION_THREE)
            SM%N_CELLSG = (SM%IBAR+2) * (SM%JBAR+2) * (SM%KBAR+2) 
      END SELECT
       
      IBAR0=IBAR0/2
      IF (.NOT.TWO_D) JBAR0=JBAR0/2
      KBAR0=KBAR0/2
   
   ENDDO LEVEL_LOOP

ENDDO MESHES_LOOP
 
END SUBROUTINE SCARC_SETUP_MESHES
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Setup communication structure for data exchange 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_EXCHANGE 
INTEGER :: NOM, NM, NL, NLMIN, NLMAX
INTEGER :: IERR, IREFINE
TYPE (MESH_TYPE)       , POINTER :: M
TYPE (SCARC_TYPE)      , POINTER :: S
TYPE (OMESH_TYPE)      , POINTER :: OM
TYPE (OSCARC_TYPE)     , POINTER :: OS
TYPE (OSCARC_MESH_TYPE), POINTER :: OSMF, OSMC

IERR = 0

!!! Initialize communication counter for ScaRC, use same TAG for all communications
N_EXCHANGES = 0
TAG_SCARC  = 99

MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX
   
   M => MESHES(NM)
   S => SCARC(NM)
   
   
   !!!----------------------------------------------------------------------------------------------------
   !!! Initialize level structures on neighboring meshes
   !!! Note that OSMF%IJKW corresponds to MESHES(NM)%OMESH(NOM)%IJKW
   !!!----------------------------------------------------------------------------------------------------
   ALLOCATE (S%OSCARC(NMESHES), STAT=IERR)
   CALL CHKMEMERR ('SCARC_SETUP_EXCHANGE', 'OSCARC', IERR)
      
   OTHER_MESHES_LOOP: DO NOM = 1, NMESHES
      
      IF (NOM == NM) CYCLE OTHER_MESHES_LOOP
   
      !!! Search for this neighbor in IJKW to determine the IOR-value of the boundary
  !    FOUND = .FALSE.
  !    SEARCH_NEIGHBOR_LOOP: DO IW = 1, SMF%N_EXTERNAL_WALL_CELLS
  !       IF (SMF%IJKW(9,IW) == NOM) THEN
  !          FOUND = .TRUE.
  !          EXIT SEARCH_NEIGHBOR_LOOP
  !       ENDIF
  !    ENDDO SEARCH_NEIGHBOR_LOOP
  !
  !    IF (.NOT.FOUND) CYCLE OTHER_MESHES_LOOP

      OM => M%OMESH(NOM)
      OS => S%OSCARC(NOM)
   
      OS%I_MIN_S  = OM%I_MIN_S
      OS%I_MAX_S  = OM%I_MAX_S
      OS%J_MIN_S  = OM%J_MIN_S
      OS%J_MAX_S  = OM%J_MAX_S
      OS%K_MIN_S  = OM%K_MIN_S
      OS%K_MAX_S  = OM%K_MAX_S
      OS%NIC_S    = OM%NIC_S
      OS%NICMAX_S = OM%NIC_S
       
      OS%I_MIN_R  = OM%I_MIN_R
      OS%I_MAX_R  = OM%I_MAX_R
      OS%J_MIN_R  = OM%J_MIN_R
      OS%J_MAX_R  = OM%J_MAX_R
      OS%K_MIN_R  = OM%K_MIN_R
      OS%K_MAX_R  = OM%K_MAX_R
      OS%NIC_R    = OM%NIC_R
      OS%NICMAX_R = OM%NIC_R
   
      IF (OS%NICMAX_S==0 .AND. OS%NICMAX_R==0)  CYCLE OTHER_MESHES_LOOP

      N_EXCHANGES  = N_EXCHANGES+1
   
      !!! Allocate OSCARC grid structure for mesh NM
      ALLOCATE (OS%MESHES(NLEVEL_MIN:NLEVEL_MAX), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_EXCHANGE', 'OS%MESH', IERR)
   
      !!! point to grid structure of OSCARC(NM) on finest level
      OSMF => OS%MESHES(NLEVEL_MIN)
   
      OSMF%IJKW => MESHES(NOM)%IJKW

      OSMF%IBAR =  MESHES(NOM)%IBAR 
      OSMF%JBAR =  MESHES(NOM)%JBAR 
      OSMF%KBAR =  MESHES(NOM)%KBAR 
   
      OSMF%N_EXTERNAL_WALL_CELLS = 2*OSMF%IBAR * OSMF%JBAR + &
                                   2*OSMF%IBAR * OSMF%KBAR + &
                                   2*OSMF%JBAR * OSMF%KBAR  
   
      !!! In case of GMG with a predefined grid hierarchy allocate corresponding level-structures
      IF (NLEVEL_MAX > NLEVEL_MIN.AND.TYPE_MULTIGRID==NSCARC_MULTIGRID_GEOMETRIC) THEN                   
   
         IREFINE=1
         DO NL=NLEVEL_MIN+1,NLEVEL_MAX
   
            IREFINE=IREFINE*2
   
            OSMC => OS%MESHES(NL)                            ! pointer to coarser level
            OSMF => OS%MESHES(NL-1)                          ! pointer to finer level
   
            ! get number of internal cells and external wall cells on neighbor NOM for corresponding level
            OSMC%IBAR=OSMF%IBAR/IREFINE
            SELECT CASE (TYPE_DIMENSION)
               CASE (NSCARC_DIMENSION_TWO)
                  OSMC%JBAR=1
               CASE (NSCARC_DIMENSION_THREE)
                  OSMC%JBAR=OSMF%JBAR/IREFINE
            END SELECT
            OSMC%KBAR=OSMF%KBAR/IREFINE
   
            OSMC%N_EXTERNAL_WALL_CELLS= 2*OSMC%IBAR * OSMC%JBAR + &
                                        2*OSMC%IBAR * OSMC%KBAR + &
                                        2*OSMC%JBAR * OSMC%KBAR  
   
            ! allocate array ijkw for neighbor on corresponding level
            ALLOCATE (OSMC%IJKW(15,OSMC%N_EXTERNAL_WALL_CELLS), STAT=IERR)
            CALL CHKMEMERR ('SCARC_SETUP_EXCHANGE', 'OSMC%IJKW', IERR)
            OSMC%IJKW = 0
   
         ENDDO
      ENDIF
   
   ENDDO OTHER_MESHES_LOOP

ENDDO MESHES_LOOP
   
!!!-------------------------------------------------------------------------------------------------------
!!! Initialize communication structures on finest level (if there is more than 1 mesh) 
!!!-------------------------------------------------------------------------------------------------------
IF (NMESHES>1) THEN

   SELECT CASE (TYPE_MULTIGRID)
      CASE (NSCARC_MULTIGRID_GEOMETRIC)
         NLMIN = NLEVEL_MIN
         NLMAX = NLEVEL_MAX
      CASE DEFAULT
         NLMIN = NLEVEL_MIN
         NLMAX = NLEVEL_MIN
   END SELECT

   DO NL = NLMIN, NLMAX

      TYPE_EXCHANGE = NSCARC_EXCHANGE_INIT
      NREQ_SCARC = 0

      CALL SCARC_RECEIVE  (NL)
      CALL SCARC_EXCHANGE (NL)

   ENDDO

ENDIF

END SUBROUTINE SCARC_SETUP_EXCHANGE

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Setup neighborship structures and boundary conditions on finest level
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_WALLCELLS
USE GEOMETRY_FUNCTIONS, ONLY: SEARCH_OTHER_MESHES
INTEGER :: NM, NL
INTEGER :: IOFFSET, I_LO, J_LO, K_LO, IW_CO, IREFINE, IERR, IOR0, IOR_LAST
INTEGER :: IW, INBR(-3:3)=0
TYPE (MESH_TYPE)       , POINTER :: M
TYPE (SCARC_TYPE)      , POINTER :: S
TYPE (SCARC_MESH_TYPE) , POINTER :: SMF, SMC

IERR=0

MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX

   M => MESHES(NM)
   S => SCARC(NM)
      
   !!!----------------------------------------------------------------------------------------------------
   !!! For all solver: 
   !!! Determine array IJKW and PRESSURE_BC_INDEX on finest level
   !!!----------------------------------------------------------------------------------------------------
   !!! point to ScaRC mesh structure on fine level
   SMF => S%MESHES(NLEVEL_MIN)
   
   SMF%IJKW              => M%IJKW
   SMF%PRESSURE_BC_INDEX => M%PRESSURE_BC_INDEX
   
   SMF%N_EXTERNAL_WALL_CELLS =  M%N_EXTERNAL_WALL_CELLS
   SMF%SUBDIVISION = 0

   ALLOCATE(SMF%ADJACENT_CELL(SMF%N_EXTERNAL_WALL_CELLS), STAT=IERR)
   CALL ChkMemErr('SCARC_SETUP_BOUNDARY','ADJACENT_CELL',IERR)
   SMF%ADJACENT_CELL = 0

   IOR_LAST = 0
   DO IW = 1, M%N_EXTERNAL_WALL_CELLS

      !!! Determine boundary type for each single external wall cell
      IF (M%BOUNDARY_TYPE(IW)==OPEN_BOUNDARY) THEN
         M%PRESSURE_BC_INDEX(IW)=DIRICHLET
      ELSE IF (M%IJKW(9,IW)/=0) THEN
         M%PRESSURE_BC_INDEX(IW)=INTERNAL
      ELSE IF (M%BOUNDARY_TYPE(IW)==NULL_BOUNDARY) THEN
         M%PRESSURE_BC_INDEX(IW)=DIRICHLET
      ELSE
         M%PRESSURE_BC_INDEX(IW)=NEUMANN
      ENDIF

      !!! Store information about single boundary faces
      IOR0 = M%IJKW(4,IW)

      IF (IOR_LAST /= IOR0) SMF%SUBDIVISION(1,IOR0) = IW
      SMF%SUBDIVISION(2,IOR0) = SMF%SUBDIVISION(2,IOR0) + 1
      IF (M%IJKW(9,IW) /= 0 .AND. INBR(IOR0) /= M%IJKW(9,IW)) THEN
         SMF%SUBDIVISION(3,IOR0) = SMF%SUBDIVISION(3,IOR0) + 1
         INBR(IOR0) = M%IJKW(9,IW)
      ENDIF
      
      IOR_LAST = IOR0

   ENDDO
   
   
   !!!----------------------------------------------------------------------------------------------------
   !!! Only in case of MG-method:
   !!! Determine arrays IJKW for coarser levels
   !!!----------------------------------------------------------------------------------------------------
   MG_IF: IF (TYPE_MULTIGRID == NSCARC_MULTIGRID_GEOMETRIC) THEN
   
      IREFINE=1
      INIT_NBR_LEVEL: DO NL = NLEVEL_MIN+1, NLEVEL_MAX
      
         !!! point to SCARC grid structures on coarser and finer level
         SMC => S%MESHES(NL)
         SMF => S%MESHES(NL-1)
   
         IREFINE=IREFINE*2
   
         !!!
         !!! set number of external wall cells on coarser level and allocate corresponding arrays
         !!!
         SMC%N_EXTERNAL_WALL_CELLS = 2 * SMC%IBAR * SMC%JBAR + &
                                     2 * SMC%IBAR * SMC%KBAR + &
                                     2 * SMC%JBAR * SMC%KBAR
   
         ALLOCATE (SMC%IJKW(15,SMC%N_EXTERNAL_WALL_CELLS), STAT=IERR)
         CALL CHKMEMERR ('SCARC_SETUP_WALLCELLS', 'SMC%IJKW', IERR)
         SMC%IJKW=0
      
         ALLOCATE (SMC%PRESSURE_BC_INDEX(SMC%N_EXTERNAL_WALL_CELLS), STAT=IERR)
         CALL CHKMEMERR ('SCARC_SETUP_WALLCELLS', 'SMC%PRESSURE_BC_INDEX', IERR)
         SMC%PRESSURE_BC_INDEX=0
      
         !!! 
         !!! set wall cells for coarser grid and define corresponding IJKW
         !!! 
         IW_CO=1
   
         !!! wall cells along IOR=1
         IOFFSET = 0
         DO K_LO=1,SMC%KBAR
            DO J_LO=1,SMC%JBAR
               CALL SCARC_SETUP_IJKW(SMC%IJKW, SMF%IJKW, SMC%PRESSURE_BC_INDEX, SMF%PRESSURE_BC_INDEX, &
                                     IW_CO,  1, IOFFSET, IREFINE, 0, J_LO, K_LO, SMF%JBAR)
            ENDDO
         ENDDO
   
         !!! wall cells along IOR=-1
        IOFFSET = SMF%JBAR*SMF%KBAR
         DO K_LO=1,SMC%KBAR
            DO J_LO=1,SMC%JBAR
               CALL SCARC_SETUP_IJKW(SMC%IJKW, SMF%IJKW, SMC%PRESSURE_BC_INDEX, SMF%PRESSURE_BC_INDEX, &
                                     IW_CO, -1, IOFFSET, IREFINE, SMC%IBAR+1, J_LO, K_LO, SMF%JBAR)
            ENDDO
         ENDDO
      
         !!! wall cells along IOR=2
         IOFFSET = 2*SMF%JBAR*SMF%KBAR
         DO K_LO=1,SMC%KBAR
            DO I_LO=1,SMC%IBAR
               CALL SCARC_SETUP_IJKW(SMC%IJKW, SMF%IJKW, SMC%PRESSURE_BC_INDEX, SMF%PRESSURE_BC_INDEX, &
                                     IW_CO,  2, IOFFSET, IREFINE, I_LO, 0, K_LO, SMF%IBAR)
            ENDDO
         ENDDO
   
         !!! wall cells along IOR=-2
         IOFFSET = 2*SMF%JBAR*SMF%KBAR + SMF%IBAR*SMF%KBAR
         DO K_LO=1,SMC%KBAR
            DO I_LO=1,SMC%IBAR
               CALL SCARC_SETUP_IJKW(SMC%IJKW, SMF%IJKW, SMC%PRESSURE_BC_INDEX, SMF%PRESSURE_BC_INDEX, &
                                     IW_CO, -2, IOFFSET, IREFINE, I_LO, SMC%JBAR+1, K_LO, SMF%IBAR)
            ENDDO
         ENDDO
   
         !!! wall cells along IOR=3
         IOFFSET = 2*SMF%JBAR*SMF%KBAR + 2*SMF%IBAR*SMF%KBAR
         DO J_LO=1,SMC%JBAR
            DO I_LO=1,SMC%IBAR
               CALL SCARC_SETUP_IJKW(SMC%IJKW, SMF%IJKW, SMC%PRESSURE_BC_INDEX, SMF%PRESSURE_BC_INDEX, &
                                     IW_CO,  3, IOFFSET, IREFINE, I_LO, J_LO, 0, SMF%IBAR)
            ENDDO
         ENDDO
   
         !!! wall cells along IOR=-3
         IOFFSET = 2*SMF%JBAR*SMF%KBAR + 2*SMF%IBAR*SMF%KBAR + SMF%IBAR*SMF%JBAR
         DO J_LO=1,SMC%JBAR
            DO I_LO=1,SMC%IBAR
               CALL SCARC_SETUP_IJKW(SMC%IJKW, SMF%IJKW, SMC%PRESSURE_BC_INDEX, SMF%PRESSURE_BC_INDEX, &
                                     IW_CO, -3, IOFFSET, IREFINE, I_LO, J_LO, SMC%KBAR+1, SMF%IBAR)
            ENDDO
         ENDDO
   
         !!!
         !!! compute analogues to NIC, I_MIN, I_MAX, K_MIN, K_MAX on level 'NL'
         !!!
         CALL SCARC_SETUP_COARSE_DIMENSIONS(IREFINE, NM, NL)
   
         !!! Store information about single boundary faces
         IOR_LAST = 0
         DO IW = 1, SMC%N_EXTERNAL_WALL_CELLS

            IOR0 = SMC%IJKW(4,IW)

            IF (IOR_LAST /= IOR0) SMC%SUBDIVISION(1,IOR0) = IW
            SMC%SUBDIVISION(2,IOR0) = SMC%SUBDIVISION(2,IOR0) + 1
            IF (M%IJKW(9,IW) /= 0 .AND. INBR(IOR0) /= M%IJKW(9,IW)) THEN
               SMC%SUBDIVISION(3,IOR0) = SMC%SUBDIVISION(3,IOR0) + 1
               INBR(IOR0) = SMC%IJKW(9,IW)
            ENDIF
            
            IOR_LAST = IOR0

         ENDDO
   
      ENDDO INIT_NBR_LEVEL
   ENDIF MG_IF

ENDDO MESHES_LOOP
   
END SUBROUTINE SCARC_SETUP_WALLCELLS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Set wall cell information on coarse level
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_IJKW(IJKW_CO, IJKW_HI, PBCI_LO, PBCI_HI, &
                            IW_CO, IOR0, IOFFSET, IREFINE, IP, JP, KP, ILEN_HI)
INTEGER, DIMENSION(:,:), INTENT(OUT) :: IJKW_CO
INTEGER, DIMENSION(:,:), INTENT(IN)  :: IJKW_HI
INTEGER, DIMENSION(:),   INTENT(OUT) :: PBCI_LO
INTEGER, DIMENSION(:),   INTENT(IN)  :: PBCI_HI
INTEGER, INTENT(INOUT) :: IW_CO
INTEGER, INTENT(IN) :: IOR0, IOFFSET, IREFINE
INTEGER, INTENT(IN) :: IP, JP, KP, ILEN_HI
INTEGER :: IW_HI(4) , IBC_HI(4), NOM_HI(4)
INTEGER :: I, I1, I2, J1, J2, K1, K2
INTEGER :: IDIFF, JDIFF, KDIFF

!!!
!!! Set IJKW_CO(1:8, IW_CO)
!!!
SELECT CASE (IOR0)
   CASE (1)
      CALL SCARC_SETUP_IJKW1(IJKW_CO, IW_CO, IOR0, IP, JP, KP, IP+1, JP  , KP  )
   CASE (-1)
      CALL SCARC_SETUP_IJKW1(IJKW_CO, IW_CO, IOR0, IP, JP, KP, IP-1, JP  , KP  )
   CASE (2)
      CALL SCARC_SETUP_IJKW1(IJKW_CO, IW_CO, IOR0, IP, JP, KP, IP  , JP+1, KP  )
   CASE (-2)
      CALL SCARC_SETUP_IJKW1(IJKW_CO, IW_CO, IOR0, IP, JP, KP, IP  , JP-1, KP  )
   CASE (3)
      CALL SCARC_SETUP_IJKW1(IJKW_CO, IW_CO, IOR0, IP, JP, KP, IP  , JP  , KP+1)
   CASE (-3)
      CALL SCARC_SETUP_IJKW1(IJKW_CO, IW_CO, IOR0, IP, JP, KP, IP  , JP  , KP-1)
END SELECT


SELECT CASE (TYPE_DIMENSION)

   !!!----------------------------------------------------------------------------------------------------
   !!! 2D-version
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_DIMENSION_TWO)

      !!! set IJKW(1:8,IW_CO) for coarser grid IW_CO
      SELECT CASE (ABS(IOR0))
         CASE ( 1)
            IW_HI(1) = IOFFSET + 2*KP-1
         CASE ( 2)
            IW_HI(1) = IOFFSET + (2*KP-1)*ILEN_HI + 2*IP - 1
         CASE ( 3)
            IW_HI(1) = IOFFSET + 2*IP-1
      END SELECT
      IW_HI(2) = IW_HI(1)+1
   
      !!! set neighbors IJKW(9,IW_CO) for coarser grid IW_CO
      NOM_HI(1) = IJKW_HI(9,IW_HI(1))
      NOM_HI(2) = IJKW_HI(9,IW_HI(2))
      IF (NOM_HI(1)/=NOM_HI(2)) THEN
         WRITE(*,*) 'SCARC_SETUP_IJKW: Inconsistent neighbors on IOR=', IOR0,' not allowed!'
         STOP
      ENDIF
   
      IJKW_CO(9,IW_CO)=NOM_HI(1) 
   
      !!! set corresponding pressure_bc_index on coarser level
      IBC_HI(1) = PBCI_HI(IW_HI(1))
      IBC_HI(2) = PBCI_HI(IW_HI(2))
      IF (IBC_HI(1)==INTERNAL .OR. IBC_HI(2)==INTERNAL) THEN
         PBCI_LO(IW_CO)=INTERNAL
      ELSE IF (IBC_HI(1)==DIRICHLET .OR. IBC_HI(2)==DIRICHLET) THEN
         PBCI_LO(IW_CO)=DIRICHLET
      ELSE
         PBCI_LO(IW_CO)=NEUMANN
      ENDIF
   
      !!! in case of an internal boundary set IJKW(10:15,IW_CO)
      IF (NOM_HI(1) > 0) THEN   
   
         J1 = 1
         J2 = 1
         SELECT CASE (ABS(IOR0))
            CASE (1)
               KDIFF = IJKW_HI(12, IW_HI(2)) - IJKW_HI(12, IW_HI(1))
               IF (KDIFF == 1) THEN
                  K1 = IJKW_HI(15,IW_HI(2))/2
                  K2 = K1
               ELSE IF (KDIFF == 2) THEN
                  K1 = IJKW_HI(15,IW_HI(1))/2     
                  K2 = IJKW_HI(15,IW_HI(2))/2    
               ELSE IF (KDIFF == 0) THEN
                  K1 = (IJKW_HI(15,IW_HI(1))+1)/2
                  K2 = K1
               ELSE
                  WRITE(*,*) 'WRONG resolutions in SCARC_SETUP_IJKW, IOR0=',IOR0
                  STOP
               ENDIF
            CASE (3)
               IDIFF = IJKW_HI(10, IW_HI(2)) - IJKW_HI(10, IW_HI(1))
               IF (IDIFF == 1) THEN
                  I1 = IJKW_HI(13,IW_HI(2))/2
                  I2 = K1
               ELSE IF (IDIFF == 2) THEN
                  I1 = IJKW_HI(13,IW_HI(1))/2
                  I1 = IJKW_HI(13,IW_HI(2))/2
               ELSE IF (IDIFF == 0) THEN
                  I1 = (IJKW_HI(13,IW_HI(1))+1)/2
                  I2 = K1
               ELSE
                  WRITE(*,*) 'WRONG resolutions in SCARC_SETUP_IJKW, IOR0=',IOR0
                  STOP
               ENDIF
         END SELECT
   
         SELECT CASE (IOR0)
            CASE (1)
               I1 = MESHES(NOM_HI(1))%IBAR/IREFINE
               I2 = I1
            CASE (-1)
               I1 = 1
               I2 = 1
            CASE (3)
               K1 = MESHES(NOM_HI(1))%KBAR/IREFINE
               K2 = K1
            CASE (-3)
               K1 = 1
               K2 = 1
         END SELECT
   
         !!!
         !!! Set IJKW_CO(10:15, IW_CO)
         !!!
         CALL SCARC_SETUP_IJKW2(IJKW_CO, IW_CO, I1, J1, K1, I2, J2, K2)
   
      ENDIF
         

!!!----------------------------------------------------------------------------------------------------
!!! 3D-version
!!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_DIMENSION_THREE)

      !!! set IJKW(1:8,IW_CO) for coarser grid IW_CO
      SELECT CASE (ABS(IOR0))
         CASE (1)
            IW_HI(1) = IOFFSET + (2*KP-2)*ILEN_HI + 2*JP - 1
            IW_HI(3) = IOFFSET + (2*KP-1)*ILEN_HI + 2*JP - 1
         CASE (2)
            IW_HI(1) = IOFFSET + (2*KP-2)*ILEN_HI + 2*IP - 1
            IW_HI(3) = IOFFSET + (2*KP-1)*ILEN_HI + 2*IP - 1
         CASE (3)
            IW_HI(1) = IOFFSET + (2*JP-2)*ILEN_HI + 2*IP - 1
            IW_HI(3) = IOFFSET + (2*JP-1)*ILEN_HI + 2*IP - 1
      END SELECT
      IW_HI(2) = IW_HI(1)+1
      IW_HI(4) = IW_HI(3)+1
   
      !!! set neighbors IJKW(9,IW_CO) for coarser grid IW_CO
      DO I=1,4
         NOM_HI(I) = IJKW_HI(9,IW_HI(I))
      ENDDO
         
      IF (NOM_HI(1)/=NOM_HI(2) .OR. NOM_HI(1)/=NOM_HI(3) .OR. NOM_HI(1)/=NOM_HI(4)) THEN
         WRITE(*,*) 'SCARC_SETUP_IJKW: Inconsistent neighbors on IOR=', IOR0,' not allowed!'
         STOP
      ENDIF
      IJKW_CO(9,IW_CO)=NOM_HI(1) 
   
      !!! set corresponding pressure_bc_index on coarser level
      DO I=1,4
         IBC_HI(I) = PBCI_HI(IW_HI(I))
      ENDDO
      IF (IBC_HI(1)==INTERNAL.OR.IBC_HI(2)==INTERNAL.OR.&
          IBC_HI(3)==INTERNAL.OR.IBC_HI(4)==INTERNAL) THEN
         PBCI_LO(IW_CO)=INTERNAL
      ELSE IF (IBC_HI(1)==DIRICHLET.OR.IBC_HI(2)==DIRICHLET.OR.&
               IBC_HI(3)==DIRICHLET.OR.IBC_HI(4)==DIRICHLET) THEN
         PBCI_LO(IW_CO)=DIRICHLET
      ELSE
         PBCI_LO(IW_CO)=NEUMANN
      ENDIF
   
      !!! in case of an internal boundary set IJKW(10:15,IW_CO)
      IF (NOM_HI(1) > 0) THEN   
   
         SELECT CASE (ABS(IOR0))
            CASE (1)
               JDIFF = IJKW_HI(11, IW_HI(2)) - IJKW_HI(11, IW_HI(1))
               KDIFF = IJKW_HI(12, IW_HI(3)) - IJKW_HI(12, IW_HI(1))
               IF (JDIFF==1 .AND. KDIFF==1) THEN
                  J1 = IJKW_HI(14,IW_HI(2))/2
                  J2 = J1
                  K1 = IJKW_HI(15,IW_HI(3))/2
                  K2 = K1
               ELSE IF (JDIFF==2 .AND. KDIFF==2) THEN
                  J1 = IJKW_HI(14,IW_HI(1))/2
                  J2 = IJKW_HI(14,IW_HI(2))/2
                  K1 = IJKW_HI(15,IW_HI(1))/2
                  K1 = IJKW_HI(15,IW_HI(3))/2
               ELSE IF (JDIFF==0 .AND. KDIFF==0) THEN
                  J1 = IJKW_HI(14,IW_HI(1))/2
                  J2 = J1
                  K1 = IJKW_HI(15,IW_HI(1))/2
                  K2 = K1
               ELSE
                  WRITE(*,*) 'WRONG resolutions in SCARC_SETUP_IJKW, IOR0=',IOR0
                  STOP
               ENDIF
            CASE (2)
               IDIFF = IJKW_HI(10, IW_HI(2)) - IJKW_HI(10, IW_HI(1))
               KDIFF = IJKW_HI(12, IW_HI(3)) - IJKW_HI(12, IW_HI(1))
               IF (IDIFF==1 .AND. KDIFF==1) THEN
                  I1 = IJKW_HI(13,IW_HI(2))/2
                  I2 = I1
                  K1 = IJKW_HI(15,IW_HI(3))/2
                  K2 = K1
               ELSE IF (IDIFF==2 .AND. KDIFF==2) THEN
                  I1 = IJKW_HI(13,IW_HI(1))/2
                  I2 = IJKW_HI(13,IW_HI(2))/2
                  K1 = IJKW_HI(15,IW_HI(1))/2
                  K1 = IJKW_HI(15,IW_HI(3))/2
               ELSE IF (IDIFF==0 .AND. KDIFF==0) THEN
                  I1 = IJKW_HI(13,IW_HI(1))/2
                  I2 = I1
                  K1 = IJKW_HI(15,IW_HI(1))/2
                  K2 = K1
               ELSE
                  WRITE(*,*) 'WRONG resolutions in SCARC_SETUP_IJKW, IOR0=',IOR0
                  STOP
               ENDIF
            CASE (3)
               IDIFF = IJKW_HI(10, IW_HI(2)) - IJKW_HI(10, IW_HI(1))
               JDIFF = IJKW_HI(10, IW_HI(2)) - IJKW_HI(10, IW_HI(1))
               IF (IDIFF==1 .AND. JDIFF==1) THEN
                  I1 = IJKW_HI(13,IW_HI(2))/2
                  I2 = I1
                  J1 = IJKW_HI(14,IW_HI(3))/2
                  J2 = J1
               ELSE IF (IDIFF==2 .AND. JDIFF==2) THEN
                  I1 = IJKW_HI(13,IW_HI(1))/2
                  I2 = IJKW_HI(13,IW_HI(2))/2
                  J1 = IJKW_HI(14,IW_HI(1))/2
                  J1 = IJKW_HI(14,IW_HI(3))/2
               ELSE IF (IDIFF==0 .AND. JDIFF==0) THEN
                  I1 = IJKW_HI(13,IW_HI(2))/2
                  I2 = I1
                  J1 = IJKW_HI(14,IW_HI(3))/2
                  J2 = J1
               ELSE
                  WRITE(*,*) 'WRONG resolutions in SCARC_SETUP_IJKW, IOR0=',IOR0
                  STOP
               ENDIF
         END SELECT
   
         SELECT CASE (IOR0)
            CASE (1)
               I1 = MESHES(NOM_HI(1))%IBAR/IREFINE
               I2 = I1
            CASE (-1)
               I1 = 1
               I2 = I1
            CASE (2)
               J1 = MESHES(NOM_HI(1))%JBAR/IREFINE
               J2 = J1
            CASE (-2)
               J1 = 1
               J2 = J1
            CASE (3)
               K1 = MESHES(NOM_HI(1))%KBAR/IREFINE
               K2 = K1
            CASE (-3)
               K1 = 1
               K2 = K1
         END SELECT
   
         !!!
         !!! Set IJKW_CO(10:15, IW_CO)
         !!!
         CALL SCARC_SETUP_IJKW2(IJKW_CO, IW_CO, I1, J1, K1, I2, J2, K2)
   
      ENDIF
   
END SELECT

IW_CO = IW_CO + 1

END SUBROUTINE SCARC_SETUP_IJKW


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! First part: Set correct IJKW-values in SCARC_INIT_NEIGHBOURS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_IJKW1(IJKW, IW0, IOR0, IMIN, JMIN, KMIN, IMAX, JMAX, KMAX)
INTEGER, DIMENSION(:,:), INTENT(OUT) :: IJKW
INTEGER, INTENT(IN) :: IW0, IOR0
INTEGER, INTENT(IN) :: IMIN, JMIN, KMIN, IMAX, JMAX, KMAX

IJKW(1,IW0)=IMIN           
IJKW(2,IW0)=JMIN    
IJKW(3,IW0)=KMIN 
IJKW(4,IW0)=IOR0     
IJKW(5,IW0)=0        
IJKW(6,IW0)=IMAX           
IJKW(7,IW0)=JMAX 
IJKW(8,IW0)=KMAX 

END SUBROUTINE SCARC_SETUP_IJKW1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Second part: Set correct IJKW-values in SCARC_INIT_NEIGHBOURS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_IJKW2(IJKW,IW0, IMIN, JMIN, KMIN, IMAX, JMAX, KMAX)
INTEGER, DIMENSION(:,:), INTENT(OUT) :: IJKW
INTEGER, INTENT(IN) :: IW0, IMIN, JMIN, KMIN, IMAX, JMAX, KMAX

IJKW(10,IW0)=IMIN
IJKW(11,IW0)=JMIN
IJKW(12,IW0)=KMIN
IJKW(13,IW0)=IMAX
IJKW(14,IW0)=JMAX
IJKW(15,IW0)=KMAX

END SUBROUTINE SCARC_SETUP_IJKW2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Set analogues to I_MIN, IMAX, J_MIN, J_MAX, K_MIN, K_MAX on coarse level
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_COARSE_DIMENSIONS(IREFINE, NM, NL)
INTEGER, INTENT(IN) :: IREFINE, NM, NL
INTEGER :: NOM, IMIN, IMAX, JMIN, JMAX, KMIN, KMAX, IERR, IW
LOGICAL :: FOUND

TYPE (SCARC_TYPE)  , POINTER :: S
TYPE (OSCARC_TYPE) , POINTER :: OS

IERR = 0

S => SCARC(NM)

OTHER_MESH_LOOP: DO NOM = 1, NMESHES

   IF (NOM == NM) CYCLE OTHER_MESH_LOOP

   !!! let OS point to exchange structure on OSCARC(NM)
   OS => S%OSCARC(NOM)
    
   !!! ACHTUNG: funktioniert nur für 2er-Potenz-Gitterweiten !!!!!
   IMIN=0
   IMAX=MESHES(NOM)%IBAR/IREFINE+1
   IF (.NOT.TWO_D) THEN
     JMIN=0
     JMAX=MESHES(NOM)%JBAR/IREFINE+1
   ELSE
     JMIN=0
     JMAX=2
   ENDIF
   KMIN=0
   KMAX=MESHES(NOM)%KBAR/IREFINE+1
   
   OS%NIC_S = 0
   FOUND = .FALSE.

   SEARCH_LOOP: DO IW=1,S%MESHES(NL)%N_EXTERNAL_WALL_CELLS
   
      ! neighborship structure already known from finest level
      IF (MESHES(NM)%IJKW(9,IW)/=NOM) CYCLE SEARCH_LOOP
      OS%NIC_S = OS%NIC_S + 1
      FOUND = .TRUE.
   
      SELECT CASE (IJKW(4,IW))
         CASE ( 1)
            IMIN=MAX(IMIN,IJKW(10,IW)-1)
         CASE (-1) 
            IMAX=MIN(IMAX,IJKW(13,IW))
         CASE ( 2) 
            JMIN=MAX(JMIN,IJKW(11,IW)-1)
         CASE (-2) 
            JMAX=MIN(JMAX,IJKW(14,IW))
         CASE ( 3) 
            KMIN=MAX(KMIN,IJKW(12,IW)-1)
         CASE (-3)
            KMAX=MIN(KMAX,IJKW(15,IW))
      END SELECT
   ENDDO SEARCH_LOOP
   
   IF (.NOT.FOUND) CYCLE OTHER_MESH_LOOP

   N_EXCHANGES = N_EXCHANGES+1

   OS%I_MIN_R = IMIN
   OS%I_MAX_R = IMAX
   OS%J_MIN_R = JMIN
   OS%J_MAX_R = JMAX
   OS%K_MIN_R = KMIN
   OS%K_MAX_R = KMAX
   
ENDDO OTHER_MESH_LOOP

END SUBROUTINE SCARC_SETUP_COARSE_DIMENSIONS

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Allocate several global structures for data exchange 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_GLOBALS
INTEGER :: NM, NL, IERR
INTEGER :: IREFINE

IERR = 0

!!!-------------------------------------------------------------------------------------------------------
!!! Allocate arrays which are used for (global) communciations
!!!-------------------------------------------------------------------------------------------------------
ALLOCATE(NC_GLOBAL(NLEVEL_MIN:NLEVEL_MAX), STAT=IERR)
CALL CHKMEMERR ('SCARC_SETUP', 'NC_GLOBAL', IERR)
NC_GLOBAL = 0

ALLOCATE(NC_LOCAL(NMESHES), STAT=IERR)
CALL CHKMEMERR ('SCARC_SETUP', 'NC_LOCAL', IERR)
NC_LOCAL = 0

ALLOCATE(SP_LOCAL(NMESHES), STAT=IERR)
CALL CHKMEMERR ('SCARC_SETUP_GLOBAL', 'SP_LOCAL', IERR)
SP_LOCAL = 0.0_EB

!!!-------------------------------------------------------------------------------------------------------
!!! Compute global number of cells for all levels
!!!-------------------------------------------------------------------------------------------------------
IREFINE=0
LEVEL_LOOP: DO NL = NLEVEL_MIN, NLEVEL_MAX

   MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX
      SELECT CASE (TYPE_SYSTEM)
         CASE (NSCARC_SYSTEM_BANDED)
            NC_LOCAL(NM)=SCARC(NM)%BANDED(NL)%NC
         CASE (NSCARC_SYSTEM_COMPACT)
            NC_LOCAL(NM)=SCARC(NM)%COMPACT(NL)%NC
      END SELECT
   ENDDO MESHES_LOOP

   !!! Determine global number of cells for all levels 
   IF (NMESHES>1) THEN
      IF (USE_MPI) THEN
         CALL MPI_ALLREDUCE(NC_LOCAL,NC_GLOBAL(NL),NMESHES,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERR)
      ELSE
         NC_GLOBAL(NL)=0
         DO NM=1,NMESHES
            NC_GLOBAL(NL) = NC_GLOBAL(NL) + NC_LOCAL(NM)
         ENDDO
      ENDIF
   ELSE
      NC_GLOBAL(NL) = NC_LOCAL(1)
   ENDIF

ENDDO LEVEL_LOOP

END SUBROUTINE SCARC_SETUP_GLOBALS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Allocate several global structures for data exchange 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_NEIGHBORS
INTEGER :: NM, NOM, NL, IERR
TYPE (OSCARC_TYPE), POINTER :: OS, OSO

IERR = 0

!!!
!!! Exchange information about lengths of abutting faces
!!!
IF (NMESHES>1) THEN
   ALLOCATE (REQ_SCARC(N_EXCHANGES*40))
   CALL CHKMEMERR ('SCARC_SETUP_GLOBAL', 'REQ_SCARC', IERR)
   REQ_SCARC = MPI_REQUEST_NULL
ENDIF

IF (NMESHES < 0) THEN

   LEVEL_LOOP: DO NL = NLEVEL_MIN+1, NLEVEL_MAX
   
      MESHES_LOOP: DO NM=NMESHES_MIN,NMESHES_MAX
         OMESHES_LOOP: DO NOM=1,NMESHES
   
            OS  => SCARC(NM)%OSCARC(NOM)
            OSO => SCARC(NOM)%OSCARC(NM)
   
            IF (USE_MPI) THEN
               CALL MPI_SEND(OS%I_MIN_R,1,MPI_INTEGER,PROCESS(NOM),1,MPI_COMM_WORLD,IERR)
               CALL MPI_SEND(OS%I_MAX_R,1,MPI_INTEGER,PROCESS(NOM),1,MPI_COMM_WORLD,IERR)
               CALL MPI_SEND(OS%J_MIN_R,1,MPI_INTEGER,PROCESS(NOM),1,MPI_COMM_WORLD,IERR)
               CALL MPI_SEND(OS%J_MAX_R,1,MPI_INTEGER,PROCESS(NOM),1,MPI_COMM_WORLD,IERR)
               CALL MPI_SEND(OS%K_MIN_R,1,MPI_INTEGER,PROCESS(NOM),1,MPI_COMM_WORLD,IERR)
               CALL MPI_SEND(OS%K_MAX_R,1,MPI_INTEGER,PROCESS(NOM),1,MPI_COMM_WORLD,IERR)
               CALL MPI_SEND(OS%NIC_S,  1,MPI_INTEGER,PROCESS(NOM),1,MPI_COMM_WORLD,IERR)
            ELSE
               OSO%I_MIN_S = OS%I_MIN_R
               OSO%I_MAX_S = OS%I_MAX_R
               OSO%J_MIN_S = OS%J_MIN_R
               OSO%J_MAX_S = OS%J_MAX_R
               OSO%K_MIN_S = OS%K_MIN_R
               OSO%K_MAX_S = OS%K_MAX_R
               OSO%NIC_R   = OS%NIC_S
            ENDIF

         ENDDO OMESHES_LOOP
      ENDDO MESHES_LOOP
   
      MESHES_LOOP2: DO NM=1,NMESHES
         OMESHES_LOOP2: DO NOM=NMESHES_MIN,NMESHES_MAX
   
            OSO => SCARC(NOM)%OSCARC(NM)
   
            IF (USE_MPI) THEN
               CALL MPI_RECV(OSO%I_MIN_S,1,MPI_INTEGER,PROCESS(NM),1,MPI_COMM_WORLD,STATUS2_SCARC,IERR)
               CALL MPI_RECV(OSO%I_MAX_S,1,MPI_INTEGER,PROCESS(NM),1,MPI_COMM_WORLD,STATUS2_SCARC,IERR)
               CALL MPI_RECV(OSO%J_MIN_S,1,MPI_INTEGER,PROCESS(NM),1,MPI_COMM_WORLD,STATUS2_SCARC,IERR)
               CALL MPI_RECV(OSO%J_MAX_S,1,MPI_INTEGER,PROCESS(NM),1,MPI_COMM_WORLD,STATUS2_SCARC,IERR)
               CALL MPI_RECV(OSO%K_MIN_S,1,MPI_INTEGER,PROCESS(NM),1,MPI_COMM_WORLD,STATUS2_SCARC,IERR)
               CALL MPI_RECV(OSO%K_MAX_S,1,MPI_INTEGER,PROCESS(NM),1,MPI_COMM_WORLD,STATUS2_SCARC,IERR)
               CALL MPI_RECV(OSO%NIC_R,  1,MPI_INTEGER,PROCESS(NM),1,MPI_COMM_WORLD,STATUS2_SCARC,IERR)
            ENDIF
   
         ENDDO OMESHES_LOOP2
      ENDDO MESHES_LOOP2
   
   ENDDO LEVEL_LOOP

ENDIF

END SUBROUTINE SCARC_SETUP_NEIGHBORS

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Setup system of equation:
!!! Define matrix stencils and initialize matrices and boundary conditions on all needed levels
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_SYSTEM
INTEGER :: NM, NL, IERR
TYPE (SCARC_TYPE), POINTER :: S

IERR = 0
 
MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX

   S  => SCARC(NM)

   SELECT_SOLVER: SELECT CASE (TYPE_METHOD)
   
      !!!-------------------------------------------------------------------------------------------------
      !!! Krylov method (CG/BICG) as main solver, different preconditioners possible
      !!!-------------------------------------------------------------------------------------------------
      CASE (NSCARC_METHOD_KRYLOV)
   
         SELECT_PRECON: SELECT CASE (TYPE_PRECON)
   
            !!!-------------------------------------------------------------------------------------------
            !!! in case of multigrid as preconditioner:
            !!!-------------------------------------------------------------------------------------------
            CASE (NSCARC_PRECON_MG)
   
               SELECT_PRECON_MG: SELECT CASE (TYPE_MULTIGRID)
   
                  !!! geometric multigrid:
                  !!!    -  use banded storage technique on all levels unless otherwise specified
                  !!!    -  assemble standard n-point-matrix hierarchy on all levels (finest and all coarser)
                  CASE (NSCARC_MULTIGRID_GEOMETRIC)
                     
                     IF (SCARC_SYSTEM == 'null') TYPE_SYSTEM = NSCARC_SYSTEM_BANDED
   
                     IF (TYPE_SYSTEM == NSCARC_SYSTEM_BANDED) THEN
                        ALLOCATE (S%BANDED(NLEVEL_MIN:NLEVEL_MAX), STAT=IERR)
                        CALL CHKMEMERR ('SCARC_SETUP_SYSTEM', 'SYSTEM_BANDED', IERR)
                     ELSE
                        ALLOCATE (S%COMPACT(NLEVEL_MIN:NLEVEL_MAX), STAT=IERR)
                        CALL CHKMEMERR ('SCARC_SETUP_SYSTEM', 'SYSTEM_COMPACT', IERR)
                     ENDIF
   
                     IF (TYPE_PRECON == NSCARC_PRECON_FFT .OR. TYPE_PRECON == NSCARC_PRECON_GSTRIX) THEN
                        ALLOCATE (S%PRECON(NLEVEL_MIN:NLEVEL_MAX), STAT=IERR)
                        CALL CHKMEMERR ('SCARC_SETUP_SYSTEM', 'PRECON', IERR)
                     ENDIF
   
                     DO NL = NLEVEL_MIN, NLEVEL_MAX
                        CALL SCARC_SETUP_MATRIX  (NM, NL)
                        CALL SCARC_SETUP_BOUNDARY(NM, NL)
                     ENDDO 
   
                  !!! algebraic multigrid:
                  !!!    -  use compact storage technique on all levels (no other choise possible!)
                  !!!    -  assemble standard n-point-matrix only on finest level 
                  !!!    -  construct all coarser levels by requested coarsening strategy
                  CASE (NSCARC_MULTIGRID_ALGEBRAIC)
   
                     TYPE_SYSTEM = NSCARC_SYSTEM_COMPACT
   
                     ALLOCATE (S%COMPACT(NLEVEL_MIN:NLEVEL_MAX), STAT=IERR)
                     CALL CHKMEMERR ('SCARC_SETUP_SYSTEM', 'SYSTEM_COMPACT', IERR)
   
                     IF (TYPE_PRECON == NSCARC_PRECON_FFT .OR. TYPE_PRECON == NSCARC_PRECON_GSTRIX) THEN
                        ALLOCATE (S%PRECON(NLEVEL_MIN:NLEVEL_MIN), STAT=IERR)
                        CALL CHKMEMERR ('SCARC_SETUP_SYSTEM', 'PRECON', IERR)
                     ENDIF
                  
                     CALL SCARC_SETUP_MATRIX  (NM, NLEVEL_MIN)
                     CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)
   
               END SELECT SELECT_PRECON_MG
   
            !!!-------------------------------------------------------------------------------------------
            !!! in case of one-level preconditioners (JACOBI/SSOR/GSTRIX/FFT)
            !!!    -  use banded storage technique on finest level unless otherwise specified 
            !!!    -  assemble standard n-point-matrix on finest level 
            !!!-------------------------------------------------------------------------------------------
            CASE DEFAULT
   
               IF (SCARC_SYSTEM == 'null') TYPE_SYSTEM = NSCARC_SYSTEM_BANDED
   
               IF (TYPE_SYSTEM == NSCARC_SYSTEM_BANDED) THEN
                  ALLOCATE (S%BANDED(NLEVEL_MIN:NLEVEL_MIN), STAT=IERR)
                  CALL CHKMEMERR ('SCARC_SETUP_SYSTEM', 'SYSTEM_BANDED', IERR)
               ELSE
                  ALLOCATE (S%COMPACT(NLEVEL_MIN:NLEVEL_MIN), STAT=IERR)
                  CALL CHKMEMERR ('SCARC_SETUP_SYSTEM', 'SYSTEM_COMPACT', IERR)
               ENDIF
   
               IF (TYPE_PRECON == NSCARC_PRECON_FFT .OR. TYPE_PRECON == NSCARC_PRECON_GSTRIX) THEN
                  ALLOCATE (S%PRECON(NLEVEL_MIN:NLEVEL_MIN), STAT=IERR)
                  CALL CHKMEMERR ('SCARC_SETUP_SYSTEM', 'PRECON', IERR)
               ENDIF
                  
               CALL SCARC_SETUP_MATRIX  (NM, NLEVEL_MIN)
               CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)
   
         END SELECT SELECT_PRECON
   
      !!!-------------------------------------------------------------------------------------------------
      !!! Multigrid as main solver
      !!!-------------------------------------------------------------------------------------------------
      CASE (NSCARC_METHOD_MULTIGRID)
   
         SELECT_MG: SELECT CASE (TYPE_MULTIGRID)
   
            !!!-------------------------------------------------------------------------------------------
            !!! geometric multigrid:
            !!!    -  use banded storage technique on all levels unless otherwise specified 
            !!!    -  assemble standard n-point-matrix hierarchy on all levels 
            !!!-------------------------------------------------------------------------------------------
            CASE (NSCARC_MULTIGRID_GEOMETRIC)
   
               IF (SCARC_SYSTEM == 'null') TYPE_SYSTEM = NSCARC_SYSTEM_BANDED
   
               IF (TYPE_SYSTEM == NSCARC_SYSTEM_BANDED) THEN
                  ALLOCATE (S%BANDED(NLEVEL_MIN:NLEVEL_MAX), STAT=IERR)
                  CALL CHKMEMERR ('SCARC_SETUP_SYSTEM', 'SYSTEM_BANDED', IERR)
               ELSE
                  ALLOCATE (S%COMPACT(NLEVEL_MIN:NLEVEL_MAX), STAT=IERR)
                  CALL CHKMEMERR ('SCARC_SETUP_SYSTEM', 'SYSTEM_COMPACT', IERR)
               ENDIF
   
               DO NL = NLEVEL_MIN, NLEVEL_MAX
                  CALL SCARC_SETUP_MATRIX  (NM, NL)
                  CALL SCARC_SETUP_BOUNDARY(NM, NL)
               ENDDO 
   
            !!!-------------------------------------------------------------------------------------------
            !!! algebraic multigrid:
            !!!    -  use compact storage technique (no other choice possible!)
            !!!    -  assemble standard n-point-matrix only on finest level
            !!!    -  construct all coarser levels later by requested coarsening strategy
            !!!-------------------------------------------------------------------------------------------
            CASE (NSCARC_MULTIGRID_ALGEBRAIC)
   
               TYPE_SYSTEM = NSCARC_SYSTEM_COMPACT
   
               ALLOCATE (S%COMPACT(NLEVEL_MIN:NLEVEL_MAX), STAT=IERR)
               CALL CHKMEMERR ('SCARC_SETUP_SYSTEM', 'SYSTEM_COMPACT', IERR)
   
               CALL SCARC_SETUP_MATRIX  (NM, NLEVEL_MIN)
               CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)

         END SELECT SELECT_MG
   
   END SELECT SELECT_SOLVER
   
ENDDO MESHES_LOOP


!!!-------------------------------------------------------------------------------------------------------
!!! Exchange matrix entries along internal boundaries:
!!!  - in case of GMG as solver or preconditioner: matrices of all levels must be exchanged
!!!  - in all other cases: only matrix of finest level must be exchanged
!!!-------------------------------------------------------------------------------------------------------
IF (NMESHES>1) THEN

   TYPE_EXCHANGE = NSCARC_EXCHANGE_MATRIX

   SELECT CASE (TYPE_MULTIGRID)
      CASE (NSCARC_MULTIGRID_GEOMETRIC)
         DO NL = NLEVEL_MIN, NLEVEL_MAX
            NREQ_SCARC = 0
            CALL SCARC_RECEIVE  (NL)
            CALL SCARC_EXCHANGE (NL)
         ENDDO
      CASE DEFAULT
         NREQ_SCARC = 0
         CALL SCARC_RECEIVE  (NLEVEL_MIN)
         CALL SCARC_EXCHANGE (NLEVEL_MIN)
   END SELECT

ENDIF

END SUBROUTINE SCARC_SETUP_SYSTEM


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Allocate matrix for the usual 5-point-stencil (2D) or 7-point-stencil (3D)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_MATRIX (NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: I, J, K, IC, IP, IW, IW0(-3:3), IL0(-3:3), IERR
TYPE (SCARC_MESH_TYPE)   , POINTER :: SM
TYPE (SCARC_SYSTEM_BANDED_TYPE) , POINTER :: SB
TYPE (SCARC_SYSTEM_COMPACT_TYPE), POINTER :: SC
TYPE (MESH_TYPE), POINTER :: M

M  => MESHES(NM)
SM => SCARC(NM)%MESHES(NL)

SELECT CASE (TYPE_SYSTEM)

!!!----------------------------------------------------------------------------------------------------
!!! Bandwise storage technique:
!!!
!!! The matrix is stored "banded", namely diagonal for diagonal.
!!! Although the sub- and superdiagonals are shorter than the main diagonal,
!!! all bands are stored in the same length and are simply filled with zeros
!!! at redundant positions. The 'wasting' of this storage space is justified
!!! by the possibility to use a much more efficient matrix-vector-multiplication
!!! (which must not use an expensive referencing logic)
!!!
!!! A(.,ILZ)  : lower subdiagonal corresponding to z-coordinate
!!! A(.,ILY)  : lower subdiagonal corresponding to y-coordinate  (only 3D)
!!! A(.,ILX)  : lower subdiagonal corresponding to x-coordinate
!!! A(.,ID )  : main diagonal
!!! A(.,IUX)  : upper subdiagonal corresponding to x-coordinate
!!! A(.,IUY)  : upper subdiagonal corresponding to y-coordinate  (only 3D)
!!! A(.,IUZ)  : upper subdiagonal corresponding to z-coordinate
!!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_SYSTEM_BANDED)

      SB => SCARC(NM)%BANDED(NL)

      SB%NX  = SM%IBAR
      SB%NY  = SM%JBAR
      SB%NZ  = SM%KBAR

      SB%NC  = SM%N_CELLS             
      SB%NCG = SM%N_CELLSG             

      SELECT CASE (TYPE_DIMENSION)

         !!! --------------------- 2D ---------------------
         CASE (NSCARC_DIMENSION_TWO)
      
            SB%NCPL = 5                                              ! number of couplings in matrix stencil
         
            ALLOCATE (SB%A(SB%NC, SB%NCPL), STAT=IERR)
            CALL CHKMEMERR ('SCARC_SETUP_MATRIX', 'SB%A', IERR)
            SB%A = 0.0_EB
      
            ID  = 1                                                  ! pointer to the single subdiagonals
            ILZ = 2
            ILX = 3
            IUX = 4
            IUZ = 5
         
            J = 1
            DO K = 1, SM%KBAR
               DO I = 1, SM%IBAR
              
                  IC = (K-1) * SM%IBAR + I
              
                  !!! main diagonal
                  SB%A(IC,ID) = - 2.0_EB*SM%DXI2 - 2.0_EB*SM%DZI2   
          
                  ! lower subdiagonal in z-direction
                  IF (K > 1) THEN
                     SB%A(IC,ILZ) = SM%DZI2                       
                  ELSE IF (SM%SUBDIVISION(3,3) > 0) THEN
                     IW0(3) = SM%SUBDIVISION(1,3)
                     IL0(3) = SM%SUBDIVISION(2,3)
                     IF (CELL_WITH_NEIGHBOR(I,J,K,IW,IW0(3),IL0(3),SM%IJKW)) SB%A(IC,ILZ) = -IW
                  ENDIF

                  !!! lower subdiagonal in x-direction
                  IF (I > 1) THEN
                     SB%A(IC,ILX) = SM%DXI2                      
                  ELSE IF (SM%SUBDIVISION(3,1) > 0) THEN
                     IW0(1) = SM%SUBDIVISION(1,1)
                     IL0(1) = SM%SUBDIVISION(2,1)
                     IF (CELL_WITH_NEIGHBOR(I,J,K,IW,IW0(1),IL0(1),SM%IJKW)) SB%A(IC,ILX) = -IW
                  ENDIF

                  ! upper subdiagonal in x-direction
                  IF (I < SM%IBAR) THEN
                     SB%A(IC,IUX) = SM%DXI2                     
                  ELSE IF (SM%SUBDIVISION(3,-1) > 0) THEN
                     IW0(-1) = SM%SUBDIVISION(1,-1)
                     IL0(-1) = SM%SUBDIVISION(2,-1)
                     IF (CELL_WITH_NEIGHBOR(I,J,K,IW,IW0(-1),IL0(-1),SM%IJKW)) SB%A(IC,IUX) = -IW
                  ENDIF

                  ! upper subdiagonal in z-direction
                  IF (K < SM%KBAR) THEN
                     SB%A(IC,IUZ) = SM%DZI2                    
                  ELSE  IF (SM%SUBDIVISION(3,-3) > 0) THEN
                     IW0(-3) = SM%SUBDIVISION(1,-3)
                     IL0(-3) = SM%SUBDIVISION(2,-3)
                     IF (CELL_WITH_NEIGHBOR(I,J,K,IW,IW0(-3),IL0(-3),SM%IJKW)) SB%A(IC,IUZ) = -IW
                  ENDIF
              
               ENDDO
            ENDDO
      
         !!! --------------------- 3D ---------------------
         CASE (NSCARC_DIMENSION_THREE)
    
            SB%NCPL = 7                                              ! number of couplings in matrix stencil

            ALLOCATE (SB%A(SB%NC, SB%NCPL), STAT=IERR)
            CALL CHKMEMERR ('SCARC_SETUP_MATRIX', 'SB%A', IERR)
            SB%A = 0.0_EB
      
            ID  = 1                                                  ! pointer to the single subdiagonals
            ILZ = 2
            ILY = 3
            ILX = 4
            IUX = 5
            IUY = 6
            IUZ = 7

            DO K = 1, SM%KBAR
               DO J = 1, SM%JBAR
                  DO I = 1, SM%IBAR
             
                     IC = (K-1) * SM%IBAR * SM%JBAR + (J-1) * SM%IBAR + I
             
                     !!! main diagonal
                     SB%A(IC,ID) = - 2.0_EB*SM%DXI2 - 2.0_EB*SM%DYI2 - 2.0_EB*SM%DZI2   
             
                     ! lower subdiagonal in z-direction
                     IF (K > 1) THEN
                        SB%A(IC,ILZ) = SM%DZI2                       
                     ELSE IF (SM%SUBDIVISION(3,3) > 0) THEN
                        IW0(3) = SM%SUBDIVISION(1,3)
                        IL0(3) = SM%SUBDIVISION(2,3)
                        IF (CELL_WITH_NEIGHBOR(I,J,K,IW,IW0(3),IL0(3),SM%IJKW)) SB%A(IC,ILZ) = -IW
                     ENDIF
   
                     ! lower subdiagonal in y-direction
                     IF (J > 1) THEN
                        SB%A(IC,ILY) = SM%DYI2                       
                     ELSE IF (SM%SUBDIVISION(3,2) > 0) THEN
                        IW0(2) = SM%SUBDIVISION(1,2)
                        IL0(2) = SM%SUBDIVISION(2,2)
                        IF (CELL_WITH_NEIGHBOR(I,J,K,IW,IW0(2),IL0(2),SM%IJKW)) SB%A(IC,ILY) = -IW
                     ENDIF
   
                     !!! lower subdiagonal in x-direction
                     IF (I > 1) THEN
                        SB%A(IC,ILX) = SM%DXI2                      
                     ELSE IF (SM%SUBDIVISION(3,1) > 0) THEN
                        IW0(1) = SM%SUBDIVISION(1,1)
                        IL0(1) = SM%SUBDIVISION(2,1)
                        IF (CELL_WITH_NEIGHBOR(I,J,K,IW,IW0(1),IL0(1),SM%IJKW)) SB%A(IC,ILX) = -IW
                     ENDIF
   
                     ! upper subdiagonal in x-direction
                     IF (I < SM%IBAR) THEN
                        SB%A(IC,IUX) = SM%DXI2                     
                     ELSE IF (SM%SUBDIVISION(3,-1) > 0) THEN
                        IW0(-1) = SM%SUBDIVISION(1,-1)
                        IL0(-1) = SM%SUBDIVISION(2,-1)
                        IF (CELL_WITH_NEIGHBOR(I,J,K,IW,IW0(-1),IL0(-1),SM%IJKW)) SB%A(IC,IUX) = -IW
                     ENDIF
   
                     ! upper subdiagonal in y-direction
                     IF (J < SM%JBAR) THEN
                        SB%A(IC,IUY) = SM%DYI2                     
                     ELSE IF (SM%SUBDIVISION(3,-2) > 0) THEN
                        IW0(-2) = SM%SUBDIVISION(1,-2)
                        IL0(-2) = SM%SUBDIVISION(2,-2)
                        IF (CELL_WITH_NEIGHBOR(I,J,K,IW,IW0(-2),IL0(-2),SM%IJKW)) SB%A(IC,IUY) = -IW
                     ENDIF
   
                     ! upper subdiagonal in z-direction
                     IF (K < SM%KBAR) THEN
                        SB%A(IC,IUZ) = SM%DZI2                    
                     ELSE IF (SM%SUBDIVISION(3,-3) > 0) THEN
                        IW0(-3) = SM%SUBDIVISION(1,-3)
                        IL0(-3) = SM%SUBDIVISION(2,-3)
                        IF (CELL_WITH_NEIGHBOR(I,J,K,IW,IW0(-3),IL0(-3),SM%IJKW)) SB%A(IC,IUZ) = -IW
                     ENDIF
             
                  ENDDO
               ENDDO
            ENDDO
             
      END SELECT
             
      SB%NA   = SB%NC * SB%NCPL                                      ! total number of matrix entries


!!!----------------------------------------------------------------------------------------------------
!!! Compact storage technique:
!!!
!!! Compression technique to store sparse matrices, non-zero entries are stored
!!! in a 1D-vector B(.), row after row, 
!!! each row starts with its diagonal entry followed by the other non-zero entries
!!! In order to identify each element, pointer arrays ROW and COL are needed,
!!! ROW points to the several diagonal entries in vector B(.), 
!!! COL points to the columns which non-zero entries in the matrix stencil
!!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_SYSTEM_COMPACT)

      SC => SCARC(NM)%COMPACT(NL)

      SC%NX  = SM%IBAR
      SC%NY  = SM%JBAR
      SC%NZ  = SM%KBAR

      SC%NC  = SM%N_CELLS
      SC%NCE = SM%N_CELLS
      SC%NCG = SM%N_CELLSG

      ALLOCATE (SC%A_ROW(SC%NCG+1), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_MATRIX', 'SC%A_ROW', IERR)
      SC%A_ROW = 0
      
      SELECT CASE (TYPE_DIMENSION)
      
         !!! --------------------- 2D ---------------------
         CASE (NSCARC_DIMENSION_TWO)

            SC%NCPL = 5                                        ! number of couplings in matrix stencil
            SC%NA   = SC%NCG * SC%NCPL                         ! first guess for number of matrix entries

            ALLOCATE (SC%A(SC%NA), STAT=IERR)
            CALL CHKMEMERR ('SCARC_SETUP_MATRIX', 'SC%A', IERR)
            SC%A = 0.0_EB
            
            ALLOCATE (SC%A_COL(SC%NA), STAT=IERR)
            CALL CHKMEMERR ('SCARC_SETUP_MATRIX', 'SC%A_COL', IERR)
            SC%A_COL = 0.0_EB
      
            J   = 1
            IP  = 1
            IW0 = 0
            DO K = 1, SM%KBAR
               DO I = 1, SM%IBAR
             
                  IC = (K-1) * SM%IBAR + I
  
                  !!! main diagonal
                  SC%A(IP)    = - 2.0_EB*SM%DXI2 - 2.0_EB*SM%DZI2      
                  SC%A_ROW(IC) = IP                                     
                  SC%A_COL(IP) = IC                                     
                  IP = IP + 1
         
                  ! lower subdiagonal in z-direction
                  IF (K > 1) THEN
                     SC%A(IP)    = SM%DZI2           
                     SC%A_COL(IP) = IC - SM%IBAR
                     IP = IP + 1
                  ELSE IF (SM%SUBDIVISION(3,3) > 0) THEN
                     IW0(3) = SM%SUBDIVISION(1,3)
                     IL0(3) = SM%SUBDIVISION(2,3)
                     IF (CELL_WITH_NEIGHBOR(I,J,K,IW,IW0(3),IL0(3),SM%IJKW)) THEN
                        SC%A_COL(IP) = -IW
                        IP = IP + 1
                     ENDIF
                  ENDIF
   
                  !!! lower subdiagonal in x-direction
                  IF (I > 1) THEN
                     SC%A(IP)    = SM%DXI2                            
                     SC%A_COL(IP) = IC - 1
                     IP = IP + 1
                  ELSE IF (SM%SUBDIVISION(3,1) > 0) THEN
                     IW0(1) = SM%SUBDIVISION(1,1)
                     IL0(1) = SM%SUBDIVISION(2,1)
                     IF (CELL_WITH_NEIGHBOR(I,J,K,IW,IW0(1),IL0(1),SM%IJKW)) THEN
                        SC%A_COL(IP) = -IW
                        IP = IP + 1
                     ENDIF
                  ENDIF
   
                  ! upper subdiagonal in x-direction
                  IF (I < SM%IBAR) THEN
                     SC%A(IP)    = SM%DXI2           
                     SC%A_COL(IP) = IC + 1
                     IP = IP + 1
                  ELSE IF (SM%SUBDIVISION(3,-1) > 0) THEN
                     IW0(-1) = SM%SUBDIVISION(1,-1)
                     IL0(-1) = SM%SUBDIVISION(2,-1)
                     IF (CELL_WITH_NEIGHBOR(I,J,K,IW,IW0(-1),IL0(-1),SM%IJKW)) THEN
                        SC%A_COL(IP) = -IW
                        IP = IP + 1
                     ENDIF
                  ENDIF
   
                  ! upper subdiagonal in z-direction
                  IF (K < SM%KBAR) THEN
                     SC%A(IP)    = SM%DZI2           
                     SC%A_COL(IP) = IC + SM%IBAR
                     IP = IP + 1
                  ELSE IF (SM%SUBDIVISION(3,-3) > 0) THEN 
                     IW0(-3) = SM%SUBDIVISION(1,-3)
                     IL0(-3) = SM%SUBDIVISION(2,-3)
                     IF (CELL_WITH_NEIGHBOR(I,J,K,IW,IW0(-3),IL0(-3),SM%IJKW)) THEN
                        SC%A_COL(IP) = -IW
                        IP = IP + 1
                     ENDIF
                  ENDIF
   
               ENDDO
            ENDDO
         
         !!! --------------------- 3D ---------------------
         CASE (NSCARC_DIMENSION_THREE)
         
            SC%NCPL = 7                                        ! number of couplings in matrix stencil
            SC%NA   = SC%NCG * SC%NCPL                         ! first guess for number of matrix entries

            ALLOCATE (SC%A(SC%NA), STAT=IERR)
            CALL CHKMEMERR ('SCARC_SETUP_MATRIX', 'SC%A', IERR)
            SC%A = 0.0_EB
            
            ALLOCATE (SC%A_COL(SC%NA), STAT=IERR)
            CALL CHKMEMERR ('SCARC_SETUP_MATRIX', 'SC%A_COL', IERR)
            SC%A_COL = 0.0_EB
      
            !!! Compute single matrix entries and corresponding row and column pointers
            !!! Along internal boundaries use placeholders for the neighboring matrix entries
            !!! which will be communicated in a following step
            DO K = 1, SM%KBAR
               DO J = 1, SM%JBAR
                  DO I = 1, SM%IBAR
             
                     IC = (K-1) * SM%IBAR * SM%JBAR + (J-1) * SM%IBAR + I
             
                     !!! main diagonal
                     SC%A(IP)    = - 2.0_EB*SM%DXI2 - 2.0*SM%DYI2 - 2.0_EB*SM%DZI2  
                     SC%A_ROW(IC) = IP                                     
                     SC%A_COL(IP) = IC 
                     IP = IP + 1
            
                     ! lower subdiagonal in z-direction
                     IF (K > 1) THEN
                        SC%A(IP)    = SM%DZI2           
                        SC%A_COL(IP) = IC - SM%IBAR * SM%JBAR
                        IP = IP + 1
                     ELSE IF (SM%SUBDIVISION(3,3) > 0) THEN
                        IW0(3) = SM%SUBDIVISION(1,3)
                        IL0(3) = SM%SUBDIVISION(2,3)
                        IF (CELL_WITH_NEIGHBOR(I,J,K,IW,IW0(3),IL0(3),SM%IJKW)) THEN
                           SC%A_COL(IP) = -IW
                           IP = IP + 1
                        ENDIF
                     ENDIF
     
                     !!! lower subdiagonal in y-direction
                     IF (J > 1) THEN
                        SC%A(IP)    =  SM%DYI2                            
                        SC%A_COL(IP) = IC - SM%IBAR
                        IP = IP + 1
                     ELSE IF (SM%SUBDIVISION(3,2) > 0) THEN
                        IW0(2) = SM%SUBDIVISION(1,2)
                        IL0(2) = SM%SUBDIVISION(2,2)
                        IF (CELL_WITH_NEIGHBOR(I,J,K,IW,IW0(2),IL0(2),SM%IJKW)) THEN
                           SC%A_COL(IP) = -IW
                           IP = IP + 1
                        ENDIF
                     ENDIF
     
                     !!! lower subdiagonal in x-direction
                     IF (I > 1) THEN
                        SC%A(IP)    =  SM%DXI2                            
                        SC%A_COL(IP) = IC - 1
                        IP = IP + 1
                     ELSE IF (SM%SUBDIVISION(3,1) > 0) THEN
                        IW0(1) = SM%SUBDIVISION(1,1)
                        IL0(1) = SM%SUBDIVISION(2,1)
                        IF (CELL_WITH_NEIGHBOR(I,J,K,IW,IW0(1),IL0(1),SM%IJKW)) THEN
                           SC%A_COL(IP) = -IW
                           IP = IP + 1
                        ENDIF
                     ENDIF
     
                     ! upper subdiagonal in x-direction
                     IF (I < SM%IBAR) THEN
                        SC%A(IP)    = SM%DXI2           
                        SC%A_COL(IP) = IC + 1
                        IP = IP + 1
                     ELSE IF (SM%SUBDIVISION(3,-1) > 0) THEN
                        IW0(-1) = SM%SUBDIVISION(1,-1)
                        IL0(-1) = SM%SUBDIVISION(2,-1)
                        IF (CELL_WITH_NEIGHBOR(I,J,K,IW,IW0(-1),IL0(-1),SM%IJKW)) THEN
                           SC%A_COL(IP) = -IW
                           IP = IP + 1
                        ENDIF
                     ENDIF
     
                     ! upper subdiagonal in y-direction
                     IF (J < SM%IBAR) THEN
                        SC%A(IP)    = SM%DYI2           
                        SC%A_COL(IP) = IC + SM%IBAR
                        IP = IP + 1
                     ELSE IF (SM%SUBDIVISION(3,-2) > 0) THEN
                        IW0(-2) = SM%SUBDIVISION(1,-2)
                        IL0(-2) = SM%SUBDIVISION(2,-2)
                        IF (CELL_WITH_NEIGHBOR(I,J,K,IW,IW0(-2),IL0(-2),SM%IJKW)) THEN
                           SC%A_COL(IP) = -IW
                           IP = IP + 1
                        ENDIF
                     ENDIF
     
                     ! upper subdiagonal in z-direction
                     IF (K < SM%KBAR) THEN
                        SC%A(IP)    = SM%DZI2           
                        SC%A_COL(IP) = IC + SM%IBAR * SM%JBAR
                        IP = IP + 1
                     ELSE IF (SM%SUBDIVISION(3,-3) > 0) THEN
                        IW0(-3) = SM%SUBDIVISION(1,-3)
                        IL0(-3) = SM%SUBDIVISION(2,-3)
                        IF (CELL_WITH_NEIGHBOR(I,J,K,IW,IW0(-3),IL0(-3),SM%IJKW)) THEN
                           SC%A_COL(IP) = -IW
                           IP = IP + 1
                        ENDIF
                     ENDIF
   
                  ENDDO
               ENDDO
            ENDDO
          
      END SELECT
      SC%A_ROW(SC%NC+1) = IP
      SC%NA            = IP -1                                     ! set correct number of matrix entries

END SELECT


END SUBROUTINE SCARC_SETUP_MATRIX


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Determine if cell (I,J,K) has a neighbor and, if yes, save corresponding IW-value
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
LOGICAL FUNCTION CELL_WITH_NEIGHBOR(I, J, K, IW, IW0, IL0, IJKW)
INTEGER, INTENT(IN)    :: I, J, K, IW0, IL0
INTEGER, INTENT(INOUT) :: IW
INTEGER, DIMENSION(:,:), INTENT(IN) :: IJKW

CELL_WITH_NEIGHBOR = .FALSE.

SEARCH_WALLCELL_LOOP: DO IW = IW0, IW0+IL0-1
  IF (I == IJKW(6,IW) .AND. J == IJKW(7,IW) .AND. K == IJKW(8,IW) .AND. IJKW(9,IW) /= 0) THEN
     CELL_WITH_NEIGHBOR = .TRUE.
     EXIT SEARCH_WALLCELL_LOOP
  ENDIF
ENDDO SEARCH_WALLCELL_LOOP

RETURN
END FUNCTION CELL_WITH_NEIGHBOR


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Set pointer for different structures on level NL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_BOUNDARY (NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: I, J, K, IOR0, IW, IC, NOM, BC_INDEX, IP
REAL(EB) :: DBC
TYPE (SCARC_MESH_TYPE), POINTER :: SM
TYPE (SCARC_SYSTEM_BANDED_TYPE) , POINTER :: SB
TYPE (SCARC_SYSTEM_COMPACT_TYPE), POINTER :: SC

SM  => SCARC(NM)%MESHES(NL)


SELECT CASE (TYPE_DIMENSION)

   !!!----------------------------------------------------------------------------------------------------
   !!! 2D-version
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_DIMENSION_TWO)

      WALL_CELL_LOOP2D: DO IW = 1, SM%N_EXTERNAL_WALL_CELLS
       
         IOR0 = SM%IJKW(4, IW)
         IF (ABS(IOR0) == 2) CYCLE            ! 2D: in case of y-boundary cycle
   
         I    = SM%IJKW(6, IW)
         K    = SM%IJKW(8, IW)
         NOM  = SM%IJKW(9, IW)
   
         BC_INDEX = SM%PRESSURE_BC_INDEX(IW)
   
         SELECT CASE (IOR0)
            CASE (1)
               IC = (K-1) * SM%IBAR + I
               DBC= SM%DXI2
            CASE (-1)
               IC = K * SM%IBAR
               DBC= SM%DXI2
            CASE (3)
               IC = I
               DBC= SM%DZI2
            CASE (-3)
               IC = (SM%KBAR-1) * SM%IBAR + I
               DBC= SM%DZI2
         END SELECT
   
         SELECT_STORAGE2D: SELECT CASE (TYPE_SYSTEM)
            !!!
            !!! banded storage technique
            !!!
            CASE (NSCARC_SYSTEM_BANDED)
   
               SB => SCARC(NM)%BANDED(NL)

               SELECT CASE (BC_INDEX)
                  CASE (DIRICHLET)                        ! set Dirichlet BC's along open boundary cells
                     SB%A(IC,ID) = SB%A(IC,ID) - DBC
                  !CASE (INTERNAL)                        ! do nothing along internal boundaries (only debugging)
                  CASE (NEUMANN)                          ! set Neumann BC's at all other nodes
                     SB%A(IC,ID) = SB%A(IC,ID) + DBC
               END SELECT
   
            !!!
            !!! compact storage technique
            !!!
            CASE (NSCARC_SYSTEM_COMPACT) 

               SC => SCARC(NM)%COMPACT(NL)
   
               IP = SC%A_ROW(IC)
               SELECT CASE (BC_INDEX)
                  CASE (DIRICHLET)                        ! set Dirichlet BC's along open boundary cells
                     SC%A(IP) = SC%A(IP) - DBC
                  !CASE (INTERNAL)                        ! do nothing along internal boundaries (only debugging)
                  CASE (NEUMANN)                          ! set Neumann BC's at all other nodes
                     SC%A(IP) = SC%A(IP) + DBC
               END SELECT
   
         END SELECT SELECT_STORAGE2D
    
      ENDDO WALL_CELL_LOOP2D
   

   !!!----------------------------------------------------------------------------------------------------
   !!! 3D-version
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_DIMENSION_THREE)

      WALL_CELL_LOOP3D: DO IW = 1, SM%N_EXTERNAL_WALL_CELLS
   
         IOR0 = SM%IJKW(4, IW)
         I    = SM%IJKW(6, IW)
         J    = SM%IJKW(7, IW)
         K    = SM%IJKW(8, IW)
         NOM  = SM%IJKW(9, IW)
   
         SM%ADJACENT_CELL(IW) = (K-1)*SM%IBAR*SM%JBAR + (J-1)*SM%IBAR + I

         BC_INDEX = SM%PRESSURE_BC_INDEX(IW)
   
         SELECT CASE (IOR0)
            CASE (1)
               IC = (K-1) * SM%IBAR * SM%JBAR + (J-1) * SM%IBAR + I
               DBC= SM%DXI2
            CASE (-1)
               IC = (K-1) * SM%IBAR * SM%JBAR + J * SM%IBAR 
               DBC= SM%DXI2
            CASE (2)
               IC = (K-1) * SM%IBAR * SM%JBAR + I
               DBC= SM%DYI2
            CASE (-2)
               IC = (K-1) * SM%IBAR * SM%JBAR + (SM%JBAR-1) * SM%IBAR + I
               DBC= SM%DYI2
            CASE (3)
               IC = (J-1) * SM%IBAR + I
               DBC= SM%DZI2
            CASE (-3)
               IC = (SM%KBAR-1) * SM%IBAR * SM%JBAR + (J-1) * SM%IBAR + I
               DBC= SM%DZI2
         END SELECT
   
         SELECT_STORAGE3D: SELECT CASE (TYPE_SYSTEM)
            !!!
            !!! banded storage technique
            !!!
            CASE (NSCARC_SYSTEM_BANDED)
   
               SB => SCARC(NM)%BANDED(NL)
   
               SELECT CASE (BC_INDEX)
                  CASE (DIRICHLET)                        ! set Dirichlet BC's at open and null boundary cells
                     SB%A(IC,ID) = SB%A(IC,ID) - DBC
                  !CASE (INTERNAL)                        ! do nothing along internal boundaries (only debugging)
                  CASE (NEUMANN)                          ! set Neumann BC's at all other cells
                     SB%A(IC,ID) = SB%A(IC,ID) + DBC
               END SELECT
   
            !!!
            !!! compact storage technique
            !!!
            CASE (NSCARC_SYSTEM_COMPACT) 
   
               SC => SCARC(NM)%COMPACT(NL)
   
               IP = SC%A_ROW(IC)
               SELECT CASE (BC_INDEX)
                  CASE (DIRICHLET)                        ! set Dirichlet BC's at open and null boundary cells
                     SC%A(IP) = SC%A(IP) - DBC
                  !CASE (INTERNAL)                        ! do nothing along internal boundaries (only debugging)
                  CASE (NEUMANN)                          ! set Neumann BC's at all other cells
                     SC%A(IP) = SC%A(IP) + DBC
               END SELECT
   
         END SELECT SELECT_STORAGE3D
   
      ENDDO WALL_CELL_LOOP3D
     
END SELECT

END SUBROUTINE SCARC_SETUP_BOUNDARY


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Initialize global solver methods (CG/BICG/GMG/AMG) and corresponding level structures
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_VECTORS
INTEGER :: IERR, IBP1, JBP1, KBP1
TYPE (SCARC_SYSTEM_BANDED_TYPE) , POINTER :: SB
TYPE (SCARC_SYSTEM_COMPACT_TYPE), POINTER :: SC
TYPE (SCARC_PRECON_TYPE) , POINTER :: SP
INTEGER :: NM, NL

IERR = 0

MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX
   LEVEL_LOOP: DO NL = NLEVEL_MIN, NLEVEL_MAX

      SELECT_STORAGE: SELECT CASE (TYPE_SYSTEM) 
      
         !!!-------------------------------------------------------------------------------------------------
         !!! Bandwise storage technique
         !!!-------------------------------------------------------------------------------------------------
         CASE (NSCARC_SYSTEM_BANDED)
            
            !!! let SB point to SYSTEM-BANDED structure
            SB => SCARC(NM)%BANDED(NL)
         
            IBP1 = SB%NX + 1
            JBP1 = SB%NY + 1
            KBP1 = SB%NZ + 1
      
            SELECT_METHOD_BANDED: SELECT CASE (TYPE_METHOD)
         
               !!! working and auxiliary vectors for global CG/BICG-method
               CASE (NSCARC_METHOD_KRYLOV)
         
                  ALLOCATE (SB%X(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'X', IERR)
                  SB%X = 0.0_EB
             
                  ALLOCATE (SB%F(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'B', IERR)
                  SB%F = 0.0_EB
             
                  ALLOCATE (SB%D(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'D', IERR)
                  SB%D = 0.0_EB
         
                  ALLOCATE (SB%W(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'R', IERR)
                  SB%W = 0.0_EB
             
                  ALLOCATE (SB%G(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'G', IERR)
                  SB%G = 0.0_EB
             
                  ALLOCATE (SB%Y(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'Y', IERR)
                  SB%Y = 0.0_EB
         
                  ALLOCATE (SB%Z(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'Z', IERR)
                  SB%Z = 0.0_EB
         
         
                  IF (TYPE_PRECON == NSCARC_PRECON_MG) THEN
               
                     ALLOCATE (SB%X2(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                     CALL CHKMEMERR ('SCARC', 'X2', IERR)
                     SB%X2 = 0.0_EB
               
                     ALLOCATE (SB%F2(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                     CALL CHKMEMERR ('SCARC', 'F2', IERR)
                     SB%F2 = 0.0_EB
               
                     ALLOCATE (SB%D2(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                     CALL CHKMEMERR ('SCARC', 'D2', IERR)
                     SB%D2 = 0.0_EB
         
                     IF (NL==NLEVEL_MAX) THEN
            
                        ALLOCATE (SB%W2(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                        CALL CHKMEMERR ('SCARC', 'R2', IERR)
                        SB%W2 = 0.0_EB
                
                        ALLOCATE (SB%G2(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                        CALL CHKMEMERR ('SCARC', 'G2', IERR)
                        SB%G2 = 0.0_EB
                
                        ALLOCATE (SB%Y2(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                        CALL CHKMEMERR ('SCARC', 'Y2', IERR)
                        SB%Y2 = 0.0_EB
            
                     ENDIF
                   
                  ENDIF
         
                  IF (TYPE_PRECON == NSCARC_PRECON_FFT) THEN 
                     SP => SCARC(NM)%PRECON(NL)
                     ALLOCATE (SP%FFT(1:IBP1, 1:JBP1, 1:KBP1), STAT=IERR)
                     CALL CHKMEMERR ('SCARC', 'FFT', IERR)
                     SP%FFT = 0.0_EB
                  ENDIF
         
               !!! working and auxiliary vectors for global GMG/AMG-method
               CASE (NSCARC_METHOD_MULTIGRID)
         
                  ALLOCATE (SB%X(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'X', IERR)
                  SB%X = 0.0_EB
             
                  ALLOCATE (SB%F(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'B', IERR)
                  SB%F = 0.0_EB
             
                  ALLOCATE (SB%D(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'D', IERR)
                  SB%D = 0.0_EB
         
                  IF (NL==NLEVEL_MAX) THEN
         
                     ALLOCATE (SB%W(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                     CALL CHKMEMERR ('SCARC', 'R', IERR)
                     SB%W = 0.0_EB
             
                     ALLOCATE (SB%G(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                     CALL CHKMEMERR ('SCARC', 'G', IERR)
                     SB%G = 0.0_EB
             
                     ALLOCATE (SB%Y(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                     CALL CHKMEMERR ('SCARC', 'Y', IERR)
                     SB%Y = 0.0_EB
         
                  ENDIF
         
                  IF (TYPE_PRECON == NSCARC_PRECON_FFT) THEN
                     SP => SCARC(NM)%PRECON(NL)
                     ALLOCATE (SP%FFT(1:IBP1, 1:JBP1, 1:KBP1), STAT=IERR)
                     CALL CHKMEMERR ('SCARC', 'FFT', IERR)
                     SP%FFT = 0.0_EB
                  ENDIF
         
            END SELECT SELECT_METHOD_BANDED
            
      
         !!!-------------------------------------------------------------------------------------------------
         !!! Compact storage technique
         !!!-------------------------------------------------------------------------------------------------
         CASE (NSCARC_SYSTEM_COMPACT)
      
            SC => SCARC(NM)%COMPACT(NL)
         
            SELECT_METHOD_COMPACT: SELECT CASE (TRIM(SCARC_METHOD))
         
               !!! working and auxiliary vectors for global CG/BICG-method
               CASE ('KRYLOV')
         
                  ALLOCATE (SC%X(SC%NCE), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'X', IERR)
                  SC%X = 0.0_EB
             
                  ALLOCATE (SC%F(SC%NCE), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'B', IERR)
                  SC%F = 0.0_EB
             
                  ALLOCATE (SC%D(SC%NCE), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'D', IERR)
                  SC%D = 0.0_EB
         
                  ALLOCATE (SC%W(SC%NCE), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'R', IERR)
                  SC%W = 0.0_EB
             
                  ALLOCATE (SC%G(SC%NCE), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'G', IERR)
                  SC%G = 0.0_EB
             
                  ALLOCATE (SC%Y(SC%NCE), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'Y', IERR)
                  SC%Y = 0.0_EB
         
                  ALLOCATE (SC%Z(SC%NCE), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'Z', IERR)
                  SC%Z = 0.0_EB
      
                  IF (TYPE_PRECON == NSCARC_PRECON_MG) THEN
      
                     ALLOCATE (SC%X2(SC%NCE), STAT=IERR)
                     CALL CHKMEMERR ('SCARC', 'X2', IERR)
                     SC%X2 = 0.0_EB
               
                     ALLOCATE (SC%F2(SC%NCE), STAT=IERR)
                     CALL CHKMEMERR ('SCARC', 'F2', IERR)
                     SC%F2 = 0.0_EB
               
                     ALLOCATE (SC%D2(SC%NCE), STAT=IERR)
                     CALL CHKMEMERR ('SCARC', 'D2', IERR)
                     SC%D2 = 0.0_EB
         
                     IF (NL==NLEVEL_MAX) THEN
            
                        ALLOCATE (SC%W2(SC%NCE), STAT=IERR)
                        CALL CHKMEMERR ('SCARC', 'R2', IERR)
                        SC%W2 = 0.0_EB
                
                        ALLOCATE (SC%G2(SC%NCE), STAT=IERR)
                        CALL CHKMEMERR ('SCARC', 'G2', IERR)
                        SC%G2 = 0.0_EB
                
                        ALLOCATE (SC%Y2(SC%NCE), STAT=IERR)
                        CALL CHKMEMERR ('SCARC', 'Y2', IERR)
                        SC%Y2 = 0.0_EB
            
                     ENDIF
                   
                  ENDIF
         
                  IF (TYPE_PRECON == NSCARC_PRECON_FFT) THEN
                     SP => SCARC(NM)%PRECON(NL)
                     ALLOCATE (SP%FFT(1:SC%NX+1, 1:SC%NY+1, 1:SC%NZ+1), STAT=IERR)
                     CALL CHKMEMERR ('SCARC', 'FFT', IERR)
                     SP%FFT = 0.0_EB
                  ENDIF
         
               !!! working and auxiliary vectors for global GMG/AMG-method
               CASE ('MULTIGRID')
         
                  ALLOCATE (SC%X(SC%NCE), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'X', IERR)
                  SC%X = 0.0_EB
             
                  ALLOCATE (SC%F(SC%NCE), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'B', IERR)
                  SC%F = 0.0_EB
             
                  ALLOCATE (SC%D(SC%NCE), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'D', IERR)
                  SC%D = 0.0_EB
         
                  IF (NL==NLEVEL_MAX) THEN
         
                     ALLOCATE (SC%W(SC%NCE), STAT=IERR)
                     CALL CHKMEMERR ('SCARC', 'R', IERR)
                     SC%W = 0.0_EB
             
                     ALLOCATE (SC%G(SC%NCE), STAT=IERR)
                     CALL CHKMEMERR ('SCARC', 'G', IERR)
                     SC%G = 0.0_EB
             
                     ALLOCATE (SC%Y(SC%NCE), STAT=IERR)
                     CALL CHKMEMERR ('SCARC', 'Y', IERR)
                     SC%Y = 0.0_EB
         
                  ENDIF
         
                  IF (TYPE_PRECON == NSCARC_PRECON_FFT) THEN
                     SP => SCARC(NM)%PRECON(NL)
                     ALLOCATE (SP%FFT(1:SC%NX+1, 1:SC%NY+1, 1:SC%NZ+1), STAT=IERR)
                     CALL CHKMEMERR ('SCARC', 'FFT', IERR)
                     SP%FFT = 0.0_EB
                  ENDIF
         
            END SELECT SELECT_METHOD_COMPACT
      
      END SELECT SELECT_STORAGE

   ENDDO LEVEL_LOOP
ENDDO  MESHES_LOOP

END SUBROUTINE SCARC_SETUP_VECTORS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Initialize global 3D-solver methods (cg/mg)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_COARSENING 
INTEGER :: IERR, NL, NM
TYPE (SCARC_SYSTEM_COMPACT_TYPE), POINTER :: SCF, SCC

IERR = 0

IF (TYPE_MULTIGRID /= NSCARC_MULTIGRID_ALGEBRAIC) RETURN

!!! Determine number of multigrid levels
LEVEL_LOOP: DO NL = NLEVEL_MIN, NLEVEL_MAX
   

   !!!-------------------------------------------------------------------------------------------------
   !!! Compute measures of single cells on level NL 
   !!! Perform coarsening with parallel postprocessing: Determine sets of coarse and fine cells
   !!!-------------------------------------------------------------------------------------------------
   DO NM = NMESHES_MIN, NMESHES_MAX

      SCF => SCARC(NM)%COMPACT(NL)                          ! system compact on fine level

      ALLOCATE (SCF%MEASURE(1:SCF%NCE), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'SCF%MEASURE', IERR)
      SCF%MEASURE = 0.0_EB

      ALLOCATE (SCF%CELLTYPE(1:SCF%NCE), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'SCF%CELLTYPE', IERR)
      SCF%CELLTYPE = NSCARC_CELLTYPE_NONE

   ENDDO
    

   SELECT CASE(TYPE_COARSENING)

      !!! RS3 : Original Ruge-Stüben coarsening 
      CASE (NSCARC_COARSENING_RS3)                 

         CALL SCARC_MEASURE_RS3(NL)
         CALL SCARC_COARSENING_RS3(NL)

      !!! A1:  Aggressive1 coarsening (1 path with length 2 corresponding to S^(1,2))
      CASE (NSCARC_COARSENING_A1)                  

         CALL SCARC_MEASURE_A1(NL)
         CALL SCARC_COARSENING_A1(NL)

      !!! A2:  Aggressive2 coarsening (2 paths with length 2 corresponding to S^(2,2))
      CASE (NSCARC_COARSENING_A2)                  

         CALL SCARC_MEASURE_A2(NL)
         CALL SCARC_COARSENING_A2(NL)

      !!! PMIS: Parallel Modified Independent Set coarsening
      CASE (NSCARC_COARSENING_PMIS)                

         CALL SCARC_MEASURE_PMIS(NL)
         CALL SCARC_COARSENING_PMIS(NL)

      !!! FDSRS3: ScaRC-FDS coarsening similar to RS3
      CASE (NSCARC_COARSENING_FDSRS3)                

         CALL SCARC_MEASURE_FDSA1(NL)
         CALL SCARC_COARSENING_PMIS(NL)

      !!! FDSRA1: ScaRC-FDS coarsening similar to A1
      CASE (NSCARC_COARSENING_FDSA1)                

         CALL SCARC_MEASURE_FDSA1(NL)
         CALL SCARC_COARSENING_FDSA1(NL)

   END SELECT


   !!!-------------------------------------------------------------------------------------------------
   !!! Save information of coarse and fine cell information and compute interpolation weights
   !!!-------------------------------------------------------------------------------------------------
   DO NM = NMESHES_MIN, NMESHES_MAX

      SCF => SCARC(NM)%COMPACT(NL)            

      !!! determine sizes of transfer matrices
      CALL SCARC_ASSESS_TRANSFER(NM, NL)

      !!! allocate prolongation matrix including row and column pointers
      ALLOCATE (SCF%P(SCF%NP), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'P', IERR)
      SCF%P = 0.0_EB
   
      ALLOCATE (SCF%P_ROW(SCF%NC+100), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'P_ROW', IERR)
      SCF%P_ROW = 0.0_EB
   
      ALLOCATE (SCF%P_COL(SCF%NP), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'P_COL', IERR)
      SCF%P_COL = 0.0_EB
   
      !!! allocate restriction matrix including row and column pointers
      ALLOCATE (SCF%R(SCF%NR), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'R', IERR)
      SCF%R = 0.0_EB
   
      ALLOCATE (SCF%R_ROW(SCF%NCC+1), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'R_ROW', IERR)
      SCF%R_ROW = 0.0_EB
   
      ALLOCATE (SCF%R_COL(SCF%NR), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'R_COL', IERR)
      SCF%R_COL = 0.0_EB
   
      CALL SCARC_AMG_INTERPOLATION_DIRECT(NM, NL)

      DEALLOCATE(SCF%MEASURE)
      DEALLOCATE(SCF%CELLTYPE)

   ENDDO

   !!!-------------------------------------------------------------------------------------------------
   !!! Allocate coarse grid matrix including pointer arrays
   !!! Note: number of cells on coarse level corresponds to number of c-points on fine level
   !!! Compute coarse grid matrix by multiplication with restriction and prolongation matrix:
   !!!  A_coarse := R * A_fine * P
   !!!-------------------------------------------------------------------------------------------------
   DO NM = NMESHES_MIN, NMESHES_MAX

      SCC => SCARC(NM)%COMPACT(NL+1)               ! Pointer to SYSTEM_COMPACT-structure of coarse level

      SCC%NC = SCF%NCC
   
      ALLOCATE (SCC%A(SCC%NC*9), STAT=IERR)
      CALL CHKMEMERR ('SCARC_AMG_TRANSFER', 'SCC%A', IERR)
      SCC%A = 0.0_EB
   
      ALLOCATE (SCC%A_ROW(SCC%NC+1), STAT=IERR)
      CALL CHKMEMERR ('SCARC_AMG_TRANSFER', 'SCC%A_ROW', IERR)
      SCC%A_ROW = 0
   
      ALLOCATE (SCC%A_COL(SCC%NC*9), STAT=IERR)
      CALL CHKMEMERR ('SCARC_AMG_TRANSFER', 'SCC%A_COL', IERR)
      SCC%A_COL = 0

      CALL SCARC_AMG_TRANSFER (NM, NL)

   ENDDO

   NLEVEL_MAX = NL+1

   EXIT LEVEL_LOOP

ENDDO LEVEL_LOOP

END SUBROUTINE SCARC_SETUP_COARSENING


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Setup environement in every solver call (i.e. set pointers to used vectors)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_ENVIRONMENT(CROUTINE, NSCOPE, NRHS, NL)
CHARACTER(*), INTENT(OUT) :: CROUTINE
INTEGER, INTENT(IN)  :: NSCOPE, NRHS
INTEGER, INTENT(OUT) :: NL
INTEGER:: NMETHOD


SELECT CASE (NSCOPE)
   CASE (NSCARC_SCOPE_COARSE)
      NMETHOD = NSCARC_METHOD_KRYLOV
   CASE DEFAULT
      NMETHOD = TYPE_METHOD
END SELECT

SELECT CASE (NMETHOD)
   
   !!!-------------------------------------------------------------------------------------------------------
   !!! Krylov method
   !!!-------------------------------------------------------------------------------------------------------
   CASE (NSCARC_METHOD_KRYLOV)

      SELECT CASE (TYPE_KRYLOV)
      
         !!! CG-method
         CASE (NSCARC_KRYLOV_CG)
      

            VEC_X = NSCARC_VECTOR_X
            VEC_D = NSCARC_VECTOR_D
            VEC_G = NSCARC_VECTOR_G
            VEC_Y = NSCARC_VECTOR_Y
            VEC_W = NSCARC_VECTOR_W

            VEC_F = NRHS                                             ! set correct right hand side vector
      
            SELECT CASE (NSCOPE)
               CASE (NSCARC_SCOPE_MAIN)
                  CROUTINE = 'SCARC_GLOBAL_CG_BANDED'
                  NL = NLEVEL_MIN
               CASE (NSCARC_SCOPE_COARSE)
                  CROUTINE = 'SCARC_COARSE_CG_BANDED'
                  NL = NLEVEL_MAX
            END SELECT
      
         !!! BICG-method
         CASE (NSCARC_KRYLOV_BICG)
      
            VEC_X = NSCARC_VECTOR_X
            VEC_D = NSCARC_VECTOR_D
            VEC_G = NSCARC_VECTOR_G
            VEC_Y = NSCARC_VECTOR_Y
            VEC_W = NSCARC_VECTOR_W
            VEC_Z = NSCARC_VECTOR_Z
      
            VEC_F = NRHS                                             ! set correct right hand side vector

            SELECT CASE (NSCOPE)
               CASE (NSCARC_SCOPE_MAIN)
                  CROUTINE = 'SCARC_GLOBAL_BICG_BANDED'
                  NL = NLEVEL_MIN
               CASE (NSCARC_SCOPE_COARSE)
                  CROUTINE = 'SCARC_COARSE_BICG_BANDED'
                  NL = NLEVEL_MAX
            END SELECT
      
      END SELECT
   
   !!!-------------------------------------------------------------------------------------------------------
   !!! Multigrid method
   !!!-------------------------------------------------------------------------------------------------------
   CASE (NSCARC_METHOD_MULTIGRID)

      NL = NLEVEL_MIN

      VEC_X = NSCARC_VECTOR_X
      VEC_D = NSCARC_VECTOR_D

      VEC_F = NRHS                                                   ! set correct right hand side vector

      !!! select scope (multigrid as main solver or preconditioner)
      SELECT CASE (NSCOPE)
         CASE (NSCARC_SCOPE_MAIN)
            CROUTINE = 'SCARC_GLOBAL_MULTIGRID_BANDED'
         CASE (NSCARC_SCOPE_PRECON)
            CROUTINE = 'SCARC_PRECON_MULTIGRID_BANDED'
      END SELECT

END SELECT
   
END SUBROUTINE SCARC_SETUP_ENVIRONMENT



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Determine if a cell IC is strongly coupled to another cell JC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
LOGICAL FUNCTION STRONGLY_COUPLED(A, A_ROW, IC, ICOL)
REAL(EB), DIMENSION(:), INTENT(IN) :: A
INTEGER , DIMENSION(:), INTENT(IN) :: A_ROW
INTEGER , INTENT(IN) :: IC, ICOL
INTEGER :: JCOL
REAL(EB) :: AMG_TOL, VAL_MAX

AMG_TOL = 0.25_EB
VAL_MAX = 0.00_EB

DO  JCOL= A_ROW(IC)+1, A_ROW(IC+1)-1       
   IF (JCOL /= ICOL) VAL_MAX = MAX(VAL_MAX, A(JCOL))
ENDDO

IF (A(ICOL) >= AMG_TOL * VAL_MAX) THEN
   STRONGLY_COUPLED=.TRUE.
ELSE
   STRONGLY_COUPLED=.FALSE.
ENDIF
RETURN

END FUNCTION STRONGLY_COUPLED


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Determine the measure of each cell for RS3-coarsening
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_MEASURE_RS3(NL)
INTEGER, INTENT(IN) :: NL
INTEGER :: NM
INTEGER :: IC, ICOL
TYPE (SCARC_MESH_TYPE), POINTER :: SM
TYPE (SCARC_SYSTEM_COMPACT_TYPE), POINTER :: SC

MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX

   SC => SCARC(NM)%COMPACT(NL)            
   SM => SCARC(NM)%MESHES(NL)

   LOCAL_MEASURE_LOOP: DO IC = 1, SC%NC
      DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1
         IF (STRONGLY_COUPLED(SC%A, SC%A_ROW, IC, ICOL)) SC%MEASURE(IC) = SC%MEASURE(IC) + 1.0_EB
      ENDDO
   ENDDO LOCAL_MEASURE_LOOP

ENDDO MESHES_LOOP

END SUBROUTINE SCARC_MEASURE_RS3


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Determine the measure of each cell for A1-coarsening
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_MEASURE_A1(NL)
INTEGER, INTENT(IN) :: NL
INTEGER, PARAMETER  :: NCOUPLINGS = 20
INTEGER :: NM, IC, JC, KC, ICOL, JCOL, IS, IW, JW
INTEGER  :: SCOUPLED(NCOUPLINGS), WCOUPLED(NCOUPLINGS)
TYPE (SCARC_SYSTEM_COMPACT_TYPE), POINTER :: SC

MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX

   SC => SCARC(NM)%COMPACT(NL)            

   !!! First determine the measure of each cell
   MEASURE_LOOP_A1: DO IC = 1, SC%NC

      IS = 1
      IW = 1
      SCOUPLED = 0
      WCOUPLED = 0

      DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1

         !!! set measure for all strongly coupled path-1 neighbors JC of coarse cell IC
         JC = SC%A_COL(ICOL)
         IF (STRONGLY_COUPLED(SC%A, SC%A_ROW, IC, ICOL)) THEN
            SC%MEASURE(IC)  = SC%MEASURE(IC) + 1.0_EB
            SCOUPLED(IS) = JC
            IS = IS + 1
         ENDIF

         !!! Check the neigbors KC of JC (except of IC itself) if they are path-2 neighbors to IC
         DO JCOL = SC%A_ROW(JC)+1, SC%A_ROW(JC+1)-1
            KC = SC%A_COL(JCOL)
            IF (STRONGLY_COUPLED(SC%A, SC%A_ROW, JC, JCOL).AND.SC%A_COL(JCOL)/=IC) THEN
               SC%MEASURE(IC)  = SC%MEASURE(IC) + 1.0_EB
               WCOUPLED(IW) = KC
               IW = IW + 1
            ENDIF
         ENDDO

      ENDDO

      !!! In the above JCOL loop neighbors may have been counted double
      !!! find them out and decrease measure for each of them (i.e. count them only one time)
      DOUBLE_VALUE_LOOP_A1: DO IW = 1, 20
         DO JW = IW+1, 20
            IF (WCOUPLED(JW)>0 .AND. WCOUPLED(IW) == WCOUPLED(JW)) THEN
               SC%MEASURE(IC)  = SC%MEASURE(IC) - 1.0_EB
               CYCLE DOUBLE_VALUE_LOOP_A1
            ENDIF
         ENDDO
      ENDDO DOUBLE_VALUE_LOOP_A1

   ENDDO MEASURE_LOOP_A1

ENDDO MESHES_LOOP

END SUBROUTINE SCARC_MEASURE_A1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Determine the measure of each cell for A2-coarsening
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_MEASURE_A2(NL)
INTEGER, INTENT(IN) :: NL
INTEGER, PARAMETER  :: NCOUPLINGS = 20
INTEGER :: NM, IC, JC, KC, ICOL, JCOL, IS, IW, JW
INTEGER :: SCOUPLED(NCOUPLINGS), WCOUPLED(NCOUPLINGS)
TYPE (SCARC_SYSTEM_COMPACT_TYPE), POINTER :: SC

MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX

   SC => SCARC(NM)%COMPACT(NL)            

   MEASURE_LOOP_A2: DO IC = 1, SC%NC

      IS = 1
      IW = 1
      SCOUPLED = 0
      WCOUPLED = 0

      DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1

         !!! set measure for all strongly coupled fine cell neighbors JC of coarse cell IC
         JC = SC%A_COL(ICOL)
         IF (STRONGLY_COUPLED(SC%A, SC%A_ROW, IC, ICOL)) THEN
            SC%MEASURE(IC)  = SC%MEASURE(IC) + 1.0_EB
            SCOUPLED(IS) = JC
            IS = IS + 1
         ENDIF

         !!! Check the neigbors KC of JC (except of IC itself) if they are strongly coupled to JC
         DO JCOL = SC%A_ROW(JC)+1, SC%A_ROW(JC+1)-1
            KC = SC%A_COL(JCOL)
            IF (STRONGLY_COUPLED(SC%A, SC%A_ROW, JC, JCOL).AND.SC%A_COL(JCOL)/=IC) THEN
               WCOUPLED(IW) = KC
               IW = IW + 1
            ENDIF
         ENDDO

      ENDDO

      !!! Check if some cells are counted double 
      !!! if yes, there is a double path to those cells and the measure for IC must be increased
      !!! Note: This algorithm only works for the 5- and 7-point Laplacian !!
      DOUBLE_VALUE_LOOP_A2: DO IW = 1, 20
         DO JW = IW+1, 20
            IF (WCOUPLED(JW)>0 .AND. WCOUPLED(IW) == WCOUPLED(JW)) THEN
               SC%MEASURE(IC)  =   SC%MEASURE(IC) + 1.0_EB
               CYCLE DOUBLE_VALUE_LOOP_A2
            ENDIF
         ENDDO
      ENDDO DOUBLE_VALUE_LOOP_A2
   ENDDO MEASURE_LOOP_A2

ENDDO MESHES_LOOP

END SUBROUTINE SCARC_MEASURE_A2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Determine the measure of each cell for PMIS-coarsening
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_MEASURE_PMIS(NL)
INTEGER, INTENT(IN) :: NL
INTEGER  :: NM, IC, ICOL
REAL(EB) :: RAND_NUM
TYPE (SCARC_SYSTEM_COMPACT_TYPE), POINTER :: SC

MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX

   SC => SCARC(NM)%COMPACT(NL)            

   !!! First determine the measure of each cell and add a random number between (0,1)
   MEASURE_LOOP_PMIS: DO IC = 1, SC%NC
      DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1
         IF (STRONGLY_COUPLED(SC%A, SC%A_ROW, IC, ICOL)) SC%MEASURE(IC) = SC%MEASURE(IC) + 1.0_EB
      ENDDO
      CALL RANDOM_NUMBER(RAND_NUM)
      !MEASURE(IC) = MEASURE(IC) + REAL(INT(RAND_NUM*10.0_EB),EB)/10.0_EB
      SC%MEASURE(IC) = SC%MEASURE(IC) + RAND_NUM
   ENDDO MEASURE_LOOP_PMIS

ENDDO MESHES_LOOP

IF (NMESHES > 1) THEN
   TYPE_EXCHANGE = NSCARC_EXCHANGE_MEASURE
   NREQ_SCARC = 0
   CALL SCARC_RECEIVE  (NLEVEL_MIN)
   CALL SCARC_EXCHANGE (NLEVEL_MIN)
ENDIF
END SUBROUTINE SCARC_MEASURE_PMIS



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Determine the measure of each cell for PMIS-coarsening
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_MEASURE_BDRY(NL)
INTEGER, INTENT(IN) :: NL
INTEGER  :: NM, IW, IC, ICOL
REAL(EB) :: RAND_NUM
TYPE (SCARC_MESH_TYPE), POINTER :: SM
TYPE (SCARC_SYSTEM_COMPACT_TYPE), POINTER :: SC

MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX

   SC => SCARC(NM)%COMPACT(NL)            
   SM => SCARC(NM)%MESHES(NL)            

   !!! First determine the measure of each cell and add a random number between (0,1)
   MEASURE_LOOP_BDRY: DO IW = 1, SM%N_EXTERNAL_WALL_CELLS
      IC = GET_ADJACENT_CELL(SM%IJKW, IW, SC%NX, SC%NY)
      IF (IC == -1) CYCLE MEASURE_LOOP_BDRY
      SC%MEASURE(IC)=0.0_EB
      DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1
         IF (STRONGLY_COUPLED(SC%A, SC%A_ROW, IC, ICOL)) SC%MEASURE(IC) = SC%MEASURE(IC) + 1.0_EB
      ENDDO
      CALL RANDOM_NUMBER(RAND_NUM)
      SC%MEASURE(IC) = SC%MEASURE(IC) + RAND_NUM
   ENDDO MEASURE_LOOP_BDRY

ENDDO MESHES_LOOP

END SUBROUTINE SCARC_MEASURE_BDRY


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! RS3 :  Original Ruge-Stuben method with parallel postprocessing
!!!      - Produces good C/F splittings but is inherently serial.  
!!!      - May produce AMG hierarchies with relatively high operator complexities.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_COARSENING_RS3(NL)
INTEGER , INTENT(IN)    :: NL
INTEGER  :: NM
INTEGER  :: IC, JC, KC, ICOL, JCOL
REAL(EB) :: MEASURE_MAX, EPS
TYPE (SCARC_SYSTEM_COMPACT_TYPE), POINTER :: SC


MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX

   SC => SCARC(NM)%COMPACT(NL)

   !!! Then determine CELLTYPE of the single cells
   CELLTYPE_LOOP: DO
   
      !!! get maximum (remaining) measure for all cells
      MEASURE_MAX = MAXVAL(SC%MEASURE(1:SC%NC))
      IF (MEASURE_MAX <= EPS) EXIT CELLTYPE_LOOP
   
      CELL_LOOP: DO IC = 1, SC%NC
   
         !!! Take first cell with maximum measure as next coarse cell
         IF (MATCH(MEASURE_MAX, SC%MEASURE(IC))) THEN
   
            SC%MEASURE(IC)  = 0.0_EB
            SC%CELLTYPE(IC) = NSCARC_CELLTYPE_COARSE
   
            !!! Determine set of fine cells 
            COARSE_LOOP: DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1
   
               !!! JC is set to be a fine cell which is no longer measured
               JC = SC%A_COL(ICOL)
   
               SC%MEASURE(JC)  = 0.0_EB
               SC%CELLTYPE(JC) = NSCARC_CELLTYPE_FINE
   
               !!!  increase measures of cells KC adjacent to fine cells JC based on strong couplings
               FINE_LOOP: DO JCOL = SC%A_ROW(JC)+1, SC%A_ROW(JC+1)-1
                  KC = SC%A_COL(JCOL)
                  IF (KC /= IC .AND. (SC%CELLTYPE(KC)==NSCARC_CELLTYPE_NONE)) THEN
                     SC%MEASURE(KC) = SC%MEASURE(KC) + 1.0_EB
                     MEASURE_MAX = MAX(MEASURE_MAX, SC%MEASURE(KC))
                  ENDIF
               ENDDO FINE_LOOP
   
            ENDDO COARSE_LOOP
   
            EXIT CELL_LOOP
         ENDIF
   
      ENDDO CELL_LOOP
   
   ENDDO CELLTYPE_LOOP

ENDDO MESHES_LOOP

END SUBROUTINE SCARC_COARSENING_RS3


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform A1-coarsening, i.e: path l=1, length l=2  ==> S_i^(1,2) with parallel postprocessing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_COARSENING_A1(NL)
INTEGER, INTENT(IN) :: NL
INTEGER  :: NM
INTEGER  :: IC, JC, KC, LC, ICOL, JCOL, KCOL
REAL(EB) :: MEASURE_MAX, EPS
TYPE (SCARC_SYSTEM_COMPACT_TYPE), POINTER :: SC

MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX

   SC => SCARC(NM)%COMPACT(NL)

   !!! Then determine CELLTYPE of the single cells
   CELLTYPE_LOOP_A1: DO
   
      !!! get maximum (remaining) measure for all cells
      MEASURE_MAX = MAXVAL(SC%MEASURE(1:SC%NC))
      IF (MEASURE_MAX <= EPS) EXIT CELLTYPE_LOOP_A1
   
      CELL_LOOP_A1: DO IC = 1, SC%NC
   
         !!! Take first cell with maximum measure as next coarse cell
         IF (MATCH(MEASURE_MAX, SC%MEASURE(IC))) THEN
   
            SC%MEASURE(IC)  = 0.0_EB
            SC%CELLTYPE(IC) = NSCARC_CELLTYPE_COARSE
   
            !!! First, determine set of strongly coupled path-1 neighbors 
            COARSE_LOOP_A1: DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1
   
               !!! JC is set to be a strongly coupled fine cell which is no longer measured
               JC = SC%A_COL(ICOL)
   
               SC%MEASURE(JC)  = 0.0_EB
               SC%CELLTYPE(JC) = NSCARC_CELLTYPE_SFINE
   
               !!! then determine set of weakly coupled path-2 neighbors 
               SFINE_LOOP_A1: DO JCOL = SC%A_ROW(JC)+1, SC%A_ROW(JC+1)-1
   
                  KC = SC%A_COL(JCOL)
   
                  IF (KC/=IC .AND. SC%CELLTYPE(KC)==NSCARC_CELLTYPE_NONE) THEN
   
                     SC%MEASURE(KC)  = 0.0_EB
                     SC%CELLTYPE(KC) = NSCARC_CELLTYPE_WFINE
   
                     WFINE_LOOP_A1: DO KCOL = SC%A_ROW(KC)+1, SC%A_ROW(KC+1)-1
                        LC = SC%A_COL(KCOL)
                        IF (LC /= JC .AND. (SC%CELLTYPE(LC)==NSCARC_CELLTYPE_NONE)) THEN
                           SC%MEASURE(LC) = SC%MEASURE(LC) + 1.0_EB
                           MEASURE_MAX = MAX(MEASURE_MAX, SC%MEASURE(LC))
                        ENDIF
                     ENDDO WFINE_LOOP_A1
   
                  ENDIF
   
               ENDDO SFINE_LOOP_A1
   
            ENDDO COARSE_LOOP_A1
   
            EXIT CELL_LOOP_A1
         ENDIF
   
      ENDDO CELL_LOOP_A1
   
   ENDDO CELLTYPE_LOOP_A1

ENDDO MESHES_LOOP

END SUBROUTINE SCARC_COARSENING_A1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform standard A2-coarsening, i.e: path l=2, length l=2  ==> S_i^(2,2) with parallel postprocessing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_COARSENING_A2(NL)
INTEGER, INTENT(IN) ::  NL
INTEGER, PARAMETER  :: NCOUPLINGS = 20
INTEGER  :: NM
INTEGER  :: IDIR, JDIR, IS, IW, JW, IC, JC, KC, LC, ICOL, JCOL, KCOL
REAL(EB) :: MEASURE_MAX, EPS
INTEGER  :: SCOUPLED(NCOUPLINGS), WCOUPLED(NCOUPLINGS)
TYPE (SCARC_SYSTEM_COMPACT_TYPE), POINTER :: SC

MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX

   SC => SCARC(NM)%COMPACT(NL)

   MEASURE_LOOP_A2: DO IC = 1, SC%NC
   
      IS = 1
      IW = 1
      SCOUPLED = 0
      WCOUPLED = 0
   
      DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1       
         
         !!! set measure for all strongly coupled fine cell neighbors JC of coarse cell IC
         JC = SC%A_COL(ICOL)
         IF (STRONGLY_COUPLED(SC%A, SC%A_ROW, IC, ICOL)) THEN
            SC%MEASURE(IC)  = SC%MEASURE(IC) + 1.0_EB 
            SCOUPLED(IS) = JC
            IS = IS + 1
         ENDIF
   
         !!! Check the neigbors KC of JC (except of IC itself) if they are strongly coupled to JC
         DO JCOL = SC%A_ROW(JC)+1, SC%A_ROW(JC+1)-1
            KC = SC%A_COL(JCOL)
            IF (STRONGLY_COUPLED(SC%A, SC%A_ROW, JC, JCOL).AND.SC%A_COL(JCOL)/=IC) THEN
               WCOUPLED(IW) = KC
               IW = IW + 1
            ENDIF
         ENDDO
   
      ENDDO 
   
      !!! Check if some cells are counted double 
      !!! if yes, there is a double path to those cells and the measure for IC must be increased
      !!! Note: This algorithm only works for the 5- and 7-point Laplacian !!
      DOUBLE_VALUE_LOOP_A2: DO IW = 1, 20
         DO JW = IW+1, 20
            IF (WCOUPLED(JW)>0 .AND. WCOUPLED(IW) == WCOUPLED(JW)) THEN
               SC%MEASURE(IC)  =   SC%MEASURE(IC) + 1.0_EB
               CYCLE DOUBLE_VALUE_LOOP_A2
            ENDIF
         ENDDO
      ENDDO DOUBLE_VALUE_LOOP_A2
   ENDDO MEASURE_LOOP_A2
   
   !!! Then determine CELLTYPE of the single cells
   CELLTYPE_LOOP_A2: DO
   
      !!! get maximum (remaining) measure for all cells
      MEASURE_MAX = MAXVAL(SC%MEASURE(1:SC%NC))
      IF (MEASURE_MAX <= EPS) EXIT CELLTYPE_LOOP_A2
   
      CELL_LOOP_A2: DO IC = 1, SC%NC
   
         !!! Take first cell with maximum measure as next coarse cell
         IF (MATCH(MEASURE_MAX, SC%MEASURE(IC))) THEN
   
            SC%MEASURE(IC)  = 0.0_EB
            SC%CELLTYPE(IC) = NSCARC_CELLTYPE_COARSE
   
            !!! First, determine set of strongly coupled path-1 neighbors
            IS = 1
            SCOUPLED = 0
            COARSE_LOOP_A2: DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1
   
               !!! JC is set to be a strongly coupled fine cell which is no longer measured
               JC = SC%A_COL(ICOL)
   
               SC%MEASURE(JC)  = 0.0_EB
               SCOUPLED(IS) = JC
               IS = IS + 1
   
               SC%CELLTYPE(JC) = NSCARC_CELLTYPE_SFINE
   
               !!! then determine set of weakly coupled path-2 neighbors 
               WFINE_LOOP_A2: DO JCOL = SC%A_ROW(JC)+1, SC%A_ROW(JC+1)-1
   
                  KC = SC%A_COL(JCOL)
   
                  !!! omit neighbors with in same direction as JC
                  IDIR = JC - IC
                  JDIR = KC - JC
   
                  IF (KC/=IC .AND. IDIR/=JDIR .AND. SC%CELLTYPE(KC)==NSCARC_CELLTYPE_NONE) THEN
   
                     SC%MEASURE(KC)  = 0.0_EB
                     SC%CELLTYPE(KC) = NSCARC_CELLTYPE_WFINE
   
                     DO KCOL = SC%A_ROW(KC)+1, SC%A_ROW(KC+1)-1
                        LC = SC%A_COL(KCOL)
                        IF (LC /= JC .AND. (SC%CELLTYPE(LC)==NSCARC_CELLTYPE_NONE)) THEN
                           SC%MEASURE(LC) = SC%MEASURE(LC) + 1.0_EB
                           MEASURE_MAX = MAX(MEASURE_MAX, SC%MEASURE(LC))
                        ENDIF
                     ENDDO 
   
                  ENDIF
   
               ENDDO WFINE_LOOP_A2
            ENDDO COARSE_LOOP_A2
   
            EXIT CELL_LOOP_A2
         ENDIF
   
      ENDDO CELL_LOOP_A2
   ENDDO CELLTYPE_LOOP_A2
   
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_COARSENING_A2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! PMIS: Parallel Modified Independent Set 
!!!     - Very fast construction with low operator complexity.  
!!!     - Convergence can deteriorate with increasing problem size on structured meshes.  
!!!     - Uses method similar to Luby's Maximal Independent Set algorithm.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_COARSENING_PMIS(NL)
INTEGER, INTENT(IN):: NL
INTEGER  :: NM
INTEGER  :: IC, JC, ICOL
REAL(EB) :: MEASURE_MAX, EPS
LOGICAL  :: BIGGEST
TYPE (SCARC_SYSTEM_COMPACT_TYPE), POINTER :: SC

MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX

   SC => SCARC(NM)%COMPACT(NL)

   CELLTYPE_LOOP: DO 
   
      MEASURE_MAX = MAXVAL(SC%MEASURE(1:SC%NCE))
      IF (MEASURE_MAX <= EPS) EXIT CELLTYPE_LOOP
   
      LOCAL_MAX_MEASURE_LOOP: DO IC = 1, SC%NCE
   
         IF (SC%CELLTYPE(IC) /= NSCARC_CELLTYPE_NONE) CYCLE LOCAL_MAX_MEASURE_LOOP

         !!! compare measure of IC with measures of its strongly coupled neighbors
         BIGGEST = .TRUE.
         DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1       
            IF (STRONGLY_COUPLED(SC%A, SC%A_ROW, IC, ICOL)) THEN
               JC = SC%A_COL(ICOL)
               IF (SC%MEASURE(JC) >= SC%MEASURE(IC)) BIGGEST = .FALSE.
            ENDIF
         ENDDO 
      
         !!! if IC has biggest measure set it to be a coarse cell 
         IF (BIGGEST) SC%CELLTYPE(IC) = NSCARC_CELLTYPE_COARSE
      ENDDO LOCAL_MAX_MEASURE_LOOP
   
      INDEPENDENT_SET_LOOP: DO IC = 1, SC%NCE
         !!! set all cells which are strongly coupled to a coarse cell to be a fine cell
         IF (SC%CELLTYPE(IC) == NSCARC_CELLTYPE_COARSE) THEN
            SC%MEASURE(IC) = 0.0_EB
            IF (IC>SC%NC) CYCLE INDEPENDENT_SET_LOOP
            DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1       
               JC = SC%A_COL(ICOL)
               IF (SC%CELLTYPE(JC) == NSCARC_CELLTYPE_NONE) THEN
                  SC%CELLTYPE(JC) = NSCARC_CELLTYPE_SFINE
                  SC%MEASURE(JC) = 0.0_EB
               ENDIF
            ENDDO 
         ENDIF
      ENDDO INDEPENDENT_SET_LOOP
   ENDDO CELLTYPE_LOOP
   
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_COARSENING_PMIS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Determine the measure of each cell for FDSA1-coarsening
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_MEASURE_FDSA1(NL)
INTEGER, INTENT(IN) :: NL
INTEGER  :: NM, IC, ICOL, IW
INTEGER  :: NX1, NX2, NY1, NY2, NZ1, NZ2, IX, IY, IZ
REAL(EB) :: RAND_NUM
LOGICAL  :: BFIVE
TYPE (SCARC_MESH_TYPE), POINTER :: SM
TYPE (SCARC_SYSTEM_COMPACT_TYPE), POINTER :: SC

MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX

   SC => SCARC(NM)%COMPACT(NL)            
   SM => SCARC(NM)%MESHES(NL)            

   NX1 = 2
   NX2 = SC%NX-1
   NZ1 = 2
   NZ2 = SC%NZ-1
   SELECT CASE (TYPE_DIMENSION)
      CASE (NSCARC_DIMENSION_TWO)
         NY1 = 1
         NY2 = 1
      CASE (NSCARC_DIMENSION_THREE)
         NY1 = 2
         NY2 = SC%NY-1
   END SELECT

   !!! First set measures for each cell in the interior of the mesh
   MEASURE_INTERIOR_LOOP: DO IZ = NZ1, NZ2
      DO IY = NY1, NY2
         DO IX = NX1, NX2
            IC = (IZ-1)*SC%NX*SC%NY + (IY-1)*SC%NX + IX
            SELECT CASE(TYPE_DIMENSION) 
               CASE (NSCARC_DIMENSION_TWO)
                  BFIVE = MOD(IX,2)==0.AND.MOD(IZ,2)==0
               CASE (NSCARC_DIMENSION_THREE)
                  BFIVE = MOD(IX,2)==0.AND.MOD(IY,2)==0.AND.MOD(IZ,2)==0
            END SELECT
            IF (BFIVE) THEN
               SC%MEASURE(IC)=5.0_EB
            ELSE
               SC%MEASURE(IC)=4.0_EB
            ENDIF
         ENDDO
      ENDDO
   ENDDO MEASURE_INTERIOR_LOOP

   !!! Then determine the measure of each cell along the internal boundary
   MEASURE_BDRY_LOOP: DO IW = 1, SM%N_EXTERNAL_WALL_CELLS
      IC = GET_ADJACENT_CELL(SM%IJKW, IW, SC%NX, SC%NY)
      IF (IC == -1) CYCLE MEASURE_BDRY_LOOP
      SC%MEASURE(IC)=0.0_EB
      DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1
         IF (STRONGLY_COUPLED(SC%A, SC%A_ROW, IC, ICOL)) SC%MEASURE(IC) = SC%MEASURE(IC) + 1.0_EB
      ENDDO
      CALL RANDOM_NUMBER(RAND_NUM)
      SC%MEASURE(IC) = SC%MEASURE(IC) + RAND_NUM
   ENDDO MEASURE_BDRY_LOOP

ENDDO MESHES_LOOP

IF (NMESHES > 1) THEN
   TYPE_EXCHANGE = NSCARC_EXCHANGE_MEASURE
   NREQ_SCARC = 0
   CALL SCARC_RECEIVE  (NLEVEL_MIN)
   CALL SCARC_EXCHANGE (NLEVEL_MIN)
ENDIF

END SUBROUTINE SCARC_MEASURE_FDSA1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! FDSA1: 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_COARSENING_FDSA1(NL)
INTEGER, INTENT(IN):: NL
INTEGER  :: NM
INTEGER  :: IC, JC, ICOL, IX, IY, IZ
REAL(EB) :: MEASURE_MAX, EPS
LOGICAL :: BIGGEST
TYPE (SCARC_SYSTEM_COMPACT_TYPE), POINTER :: SC

MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX

   SC => SCARC(NM)%COMPACT(NL)

   SELECT CASE (TYPE_DIMENSION)

      CASE (NSCARC_DIMENSION_TWO)

         DO IZ = 2, SC%NZ-1
            DO IX = 2, SC%NX-1

               IC = (IZ-1)*SC%NX + IX

               IF (SC%MEASURE(IC)==5.0_EB) THEN

                  SC%MEASURE(IC)  = 0.0_EB
                  SC%CELLTYPE(IC) = NSCARC_CELLTYPE_COARSE

                  CALL SET_TYPES_A1(SC%MEASURE, SC%CELLTYPE, IC-SC%NX, 1    )
                  CALL SET_TYPES_A1(SC%MEASURE, SC%CELLTYPE, IC-1    , SC%NX)
                  CALL SET_TYPES_A1(SC%MEASURE, SC%CELLTYPE, IC+SC%NX, 1    )
                  CALL SET_TYPES_A1(SC%MEASURE, SC%CELLTYPE, IC+1    , SC%NX)

               ENDIF

            ENDDO
         ENDDO 

      CASE (NSCARC_DIMENSION_THREE)

         DO IZ = 2, SC%NZ-1
            DO IY = 2, SC%NY-1
               DO IX = 2, SC%NX-1

                  IC = (IZ-1)*SC%NX*SC%NY + (IY-1)*SC%NX + IX

                  IF (SC%MEASURE(IC)==5.0_EB) THEN

                     SC%MEASURE(IC)  = 0.0_EB
                     SC%CELLTYPE(IC) = NSCARC_CELLTYPE_COARSE
   
                     CALL SET_TYPES_A1(SC%MEASURE, SC%CELLTYPE, IC-SC%NX*SC%NY, 1          )
                     CALL SET_TYPES_A1(SC%MEASURE, SC%CELLTYPE, IC-SC%NX      , SC%NX      )
                     CALL SET_TYPES_A1(SC%MEASURE, SC%CELLTYPE, IC-1          , SC%NX*SC%NY)
                     CALL SET_TYPES_A1(SC%MEASURE, SC%CELLTYPE, IC+SC%NX*SC%NY, 1          )
                     CALL SET_TYPES_A1(SC%MEASURE, SC%CELLTYPE, IC+SC%NX      , SC%NX      )
                     CALL SET_TYPES_A1(SC%MEASURE, SC%CELLTYPE, IC+1          , SC%NX*SC%NY)

                  ENDIF

               ENDDO
            ENDDO
         ENDDO 

   END SELECT

   CELLTYPE_LOOP: DO 

      MEASURE_MAX = MAXVAL(SC%MEASURE(1:SC%NCE))
      IF (MEASURE_MAX <= EPS) EXIT CELLTYPE_LOOP
   
      LOCAL_MAX_MEASURE_LOOP: DO IC = 1, SC%NCE
   
         IF (SC%CELLTYPE(IC) /= NSCARC_CELLTYPE_NONE) CYCLE LOCAL_MAX_MEASURE_LOOP

         !!! compare measure of IC with measures of its strongly coupled neighbors
         BIGGEST = .TRUE.
         DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1       
            IF (STRONGLY_COUPLED(SC%A, SC%A_ROW, IC, ICOL)) THEN
               JC = SC%A_COL(ICOL)
               IF (SC%MEASURE(JC) >= SC%MEASURE(IC)) BIGGEST = .FALSE.
            ENDIF
         ENDDO 
      
         !!! if IC has biggest measure set it to be a coarse cell 
         IF (BIGGEST) SC%CELLTYPE(IC) = NSCARC_CELLTYPE_COARSE
      ENDDO LOCAL_MAX_MEASURE_LOOP
   
      INDEPENDENT_SET_LOOP: DO IC = 1, SC%NCE
         !!! set all cells which are strongly coupled to a coarse cell to be a fine cell
         IF (SC%CELLTYPE(IC) == NSCARC_CELLTYPE_COARSE) THEN
            SC%MEASURE(IC) = 0.0_EB
            IF (IC>SC%NC) CYCLE INDEPENDENT_SET_LOOP
            DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1       
               JC = SC%A_COL(ICOL)
               IF (SC%CELLTYPE(JC) == NSCARC_CELLTYPE_NONE) THEN
                  SC%CELLTYPE(JC) = NSCARC_CELLTYPE_SFINE
                  SC%MEASURE(JC) = 0.0_EB
               ENDIF
            ENDDO 
         ENDIF
      ENDDO INDEPENDENT_SET_LOOP
   
   ENDDO CELLTYPE_LOOP
   
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_COARSENING_FDSA1


SUBROUTINE SET_TYPES_A1(MEASURE, CELLTYPE, JC, OFFSET)
REAL(EB), DIMENSION(:), INTENT(OUT) :: MEASURE
INTEGER , DIMENSION(:), INTENT(OUT) :: CELLTYPE
INTEGER , INTENT(IN) :: JC, OFFSET

MEASURE(JC)         = 0.0_EB
MEASURE(JC-OFFSET)  = 0.0_EB
MEASURE(JC+OFFSET)  = 0.0_EB

CELLTYPE(JC)        = NSCARC_CELLTYPE_SFINE
CELLTYPE(JC-OFFSET) = NSCARC_CELLTYPE_WFINE
CELLTYPE(JC+OFFSET) = NSCARC_CELLTYPE_WFINE
                     
END SUBROUTINE SET_TYPES_A1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Special coarsening strategy for internal boundaries
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_COARSENING_BDRY(NL)
INTEGER, INTENT(IN):: NL
INTEGER  :: NM
INTEGER  :: IC, JC, ICOL, IW
REAL(EB) :: MEASURE_MAX, EPS
LOGICAL  :: BIGGEST
TYPE (SCARC_MESH_TYPE), POINTER :: SM
TYPE (SCARC_SYSTEM_COMPACT_TYPE), POINTER :: SC

MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX

   SC => SCARC(NM)%COMPACT(NL)

   CELLTYPE_LOOP: DO 
   
      !!! get maximum (remaining) measure for all cells
      MEASURE_MAX = 0.0_EB
      ADJACENT_CELL_LOOP: DO IW = 1, SM%N_EXTERNAL_WALL_CELLS
         IC = GET_ADJACENT_CELL(SM%IJKW, IW, SC%NX, SC%NY)
         IF (IC == -1) CYCLE ADJACENT_CELL_LOOP
         MEASURE_MAX = MAX(MEASURE_MAX, SC%MEASURE(IC))
      ENDDO ADJACENT_CELL_LOOP
      IF (MEASURE_MAX <= EPS) EXIT CELLTYPE_LOOP
   
      MAX_MEASURE_LOOP: DO IW = 1, SM%N_EXTERNAL_WALL_CELLS
   
         IC = GET_ADJACENT_CELL(SM%IJKW, IW, SC%NX, SC%NY)
         IF (IC == -1) CYCLE MAX_MEASURE_LOOP
         IF (SC%CELLTYPE(IC) /= NSCARC_CELLTYPE_NONE) CYCLE MAX_MEASURE_LOOP

         !!! compare measure of IC with measures of its strongly coupled neighbors
         BIGGEST = .TRUE.
         DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1       
            IF (STRONGLY_COUPLED(SC%A, SC%A_ROW, IC, ICOL)) THEN
               JC = SC%A_COL(ICOL)
               IF (SC%MEASURE(JC) >= SC%MEASURE(IC)) BIGGEST = .FALSE.
            ENDIF
         ENDDO 
      
         !!! if IC has biggest measure set it to be a coarse cell 
         IF (BIGGEST) SC%CELLTYPE(IC) = NSCARC_CELLTYPE_COARSE
      ENDDO MAX_MEASURE_LOOP
   
      INDEPENDENT_SET_LOOP: DO IW = 1, SM%N_EXTERNAL_WALL_CELLS
         !!! set all cells which are strongly coupled to a coarse cell to be a fine cell
         IC = GET_ADJACENT_CELL(SM%IJKW, IW, SC%NX, SC%NY)
         IF (IC == -1) CYCLE INDEPENDENT_SET_LOOP
         IF (SC%CELLTYPE(IC) == NSCARC_CELLTYPE_COARSE) THEN
            SC%MEASURE(IC) = 0.0_EB
            DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1       
               JC = SC%A_COL(ICOL)
               IF (SC%CELLTYPE(JC) == NSCARC_CELLTYPE_NONE) THEN
                  SC%CELLTYPE(JC) = NSCARC_CELLTYPE_SFINE
                  SC%MEASURE(JC) = 0.0_EB
               ENDIF
            ENDDO 
         ENDIF
      ENDDO INDEPENDENT_SET_LOOP
   
   ENDDO CELLTYPE_LOOP
   
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_COARSENING_BDRY


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! SFDS: Own ScaRC-optimized coarsening strategy individually fitting the FDS-grid structures 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_COARSENING_FDSRS3(NL)
INTEGER, INTENT(IN) :: NL
INTEGER  :: NM
INTEGER  :: IC, JC, KC, ICOL, JCOL
INTEGER  :: NX1, NX2, NY1, NY2, NZ1, NZ2, IX, IY, IZ
REAL(EB) :: MEASURE_MAX, EPS
TYPE (SCARC_SYSTEM_COMPACT_TYPE), POINTER :: SC


MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX

   SC => SCARC(NM)%COMPACT(NL)

   NX1 = 2
   NX2 = SC%NX-1
   NZ1 = 2
   NZ2 = SC%NZ-1
   SELECT CASE (TYPE_DIMENSION)
      CASE (NSCARC_DIMENSION_TWO)
         NY1 = 1
         NY2 = 1
      CASE (NSCARC_DIMENSION_THREE)
         NY1 = 2
         NY2 = SC%NY-1
   END SELECT

   !!! Then determine CELLTYPE of the single cells
   CELLTYPE_LOOP: DO
   
      !!! get maximum (remaining) measure for all cells
      MEASURE_MAX = INTERNAL_MAX(SC%MEASURE(1:SC%NC),SC%NX,NX1,NX2,SC%NY,NY1,NY2,NZ1,NZ2)
      IF (MEASURE_MAX <= EPS) EXIT CELLTYPE_LOOP
   
      NZ_LOOP: DO IZ = NZ1, NZ2
         NY_LOOP: DO IY = NY1, NY2
            NX_LOOP: DO IX = NX1, NX2

               IC = (IZ-1) * SC%NX * SC%NY + (IY-1) * SC%NX + IX
         
               !!! Take first cell with maximum measure as next coarse cell
               IF (MATCH(MEASURE_MAX, SC%MEASURE(IC))) THEN
         
                  SC%MEASURE(IC)  = 0.0_EB
                  SC%CELLTYPE(IC) = NSCARC_CELLTYPE_COARSE
         
                  !!! Determine set of fine cells 
                  COARSE_LOOP: DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1
         
                     !!! JC is set to be a fine cell which is no longer measured
                     JC = SC%A_COL(ICOL)
         
                     SC%MEASURE(JC)  = 0.0_EB
                     SC%CELLTYPE(JC) = NSCARC_CELLTYPE_FINE
         
                     !!!  increase measures of cells KC adjacent to fine cells JC based on strong couplings
                     FINE_LOOP: DO JCOL = SC%A_ROW(JC)+1, SC%A_ROW(JC+1)-1
                        KC = SC%A_COL(JCOL)
                        IF (KC /= IC .AND. (SC%CELLTYPE(KC)==NSCARC_CELLTYPE_NONE)) THEN
                           SC%MEASURE(KC) = SC%MEASURE(KC) + 1.0_EB
                           MEASURE_MAX = MAX(MEASURE_MAX, SC%MEASURE(KC))
                        ENDIF
                     ENDDO FINE_LOOP
         
                  ENDDO COARSE_LOOP
         
                  EXIT NZ_LOOP
               ENDIF
         
            ENDDO NX_LOOP
         ENDDO NY_LOOP
      ENDDO NZ_LOOP
   
   ENDDO CELLTYPE_LOOP

ENDDO MESHES_LOOP

END SUBROUTINE SCARC_COARSENING_FDSRS3


REAL(EB) FUNCTION INTERNAL_MAX(MEASURE, NX, NX1, NX2, NY, NY1, NY2, NZ1, NZ2)
REAL(EB), DIMENSION(:), INTENT(IN) :: MEASURE
INTEGER , INTENT(IN) :: NX, NX1, NX2, NY, NY1, NY2, NZ1, NZ2
INTEGER  :: IX, IY, IZ, IC
INTERNAL_MAX = 0.0_EB
DO IZ = NZ1, NZ2
   DO IY = NY1, NY2
      DO IX = NX1, NX2
         IC = (IZ-1) * NX * NY + (IY-1) * NX + IX
         INTERNAL_MAX = MAX(INTERNAL_MAX, MEASURE(IC))
      ENDDO
   ENDDO
ENDDO
RETURN
END FUNCTION INTERNAL_MAX

INTEGER FUNCTION GET_ADJACENT_CELL(IJKW, IW, NX, NY)
INTEGER, DIMENSION(:,:), INTENT(IN) :: IJKW
INTEGER, INTENT(IN) :: IW, NX, NY
INTEGER :: IC

SELECT CASE (TYPE_DIMENSION)
   CASE (NSCARC_DIMENSION_TWO)
      IF (ABS(IJKW(4,IW)) == 2) THEN
         IC = -1
      ELSE
         IC = (IJKW(8,IW)-1)*NX + IJKW(6,IW)
      ENDIF
   CASE (NSCARC_DIMENSION_THREE)
      IC = (IJKW(8,IW)-1)*NX*NY + (IJKW(7,IW)-1)*NX + IJKW(6,IW)
END SELECT

GET_ADJACENT_CELL = IC
RETURN

END FUNCTION GET_ADJACENT_CELL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Determine total number of coarse and fine cells and sizes of prolongation and restriction matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_ASSESS_TRANSFER(NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IC, ICP
TYPE (SCARC_SYSTEM_COMPACT_TYPE), POINTER :: SC

SC => SCARC(NM)%COMPACT(NL)

!!! Determine dimensions of restriction and prolongation matrices
SC%NCC = 0
SC%NCF = 0
SC%NP  = 0
SC%NR  = 0
DO IC = 1, SC%NC
   SELECT CASE (SC%CELLTYPE(IC))
      CASE (NSCARC_CELLTYPE_COARSE)
         SC%NCC = SC%NCC + 1
         SC%NP  = SC%NP  + 1
      CASE (NSCARC_CELLTYPE_SFINE, NSCARC_CELLTYPE_WFINE)
         SC%NCF = SC%NCF + 1
         SC%NP  = SC%NP  + SC%A_ROW(IC+1)-SC%A_ROW(IC) - 1
   END SELECT
ENDDO
SC%NR = SC%NP

!!! Determine number of coarse and fine cells and check correctness of computation
IF (SC%NCC + SC%NCF /= SC%NC) THEN
   WRITE(*,*) 'Error in AMG standard coarsening, N_CELLS = ', SC%NCC + SC%NCF, ' differs from N_CELLS = ', SC%NC
!   STOP
ENDIF

!!! Determine new numbering for coarse cells
ICP   = 0
DO IC = 1, SC%NC
   IF (SC%CELLTYPE(IC) == NSCARC_CELLTYPE_COARSE) THEN
      ICP = ICP + 1
      SC%CELLTYPE(IC) = ICP
   ENDIF
ENDDO

END SUBROUTINE SCARC_ASSESS_TRANSFER


SUBROUTINE SCARC_AMG_INTERPOLATION_DIRECT(NM, NL)
INTEGER , INTENT(IN) :: NM, NL
INTEGER  :: ICP, ICP2, IC, ICOL, IDIAG, JDIAG
INTEGER  :: IP, IC2, JC, JCOL, IW, KC, IW0
REAL(EB) :: SUM_COUPLED, SUM_CPOINTS, SCAL
REAL(EB) :: VALUES(20), WEIGHTS(20)
INTEGER  :: NEIGHBOR(20), NWEIGHTS
LOGICAL  :: BFIRST
TYPE (SCARC_SYSTEM_COMPACT_TYPE), POINTER :: SC

SC => SCARC(NM)%COMPACT(NL)

IP = 1
DO IC = 1, SC%NC

   VALUES  = 0.0_EB
   WEIGHTS = 0.0_EB
   NEIGHBOR = 0

   !!!-------------------------------------------------------------------------------------------------
   !!! If IC is a coarse cell, its value is taken
   !!!-------------------------------------------------------------------------------------------------
   IF (SC%CELLTYPE(IC) >= NSCARC_CELLTYPE_COARSE) THEN

      NEIGHBOR(1)= SC%CELLTYPE(IC)
      WEIGHTS(1) = 1.0_EB
      NWEIGHTS   = 1

   !!!-------------------------------------------------------------------------------------------------
   !!! If IC is a fine cell, a mean value must be computed based on surrounding coarse cells
   !!!-------------------------------------------------------------------------------------------------
   ELSE 

      !!! Get main diagonal entry a_ii for that fine cell
      IDIAG = SC%A_ROW(IC)

      !!! Select type of fine cell (weakly/strongly coupled)
      SELECT_FPOINT_TYPE: SELECT CASE(SC%CELLTYPE(IC))

         !!!-------------------------------------------------------------------------------------------
         !!! Strongly coupled fine cell IC
         !!! approximate IC by weighted sum of surrounding strongly coupled coarse cells JC
         !!! Note: N_i: all coupled neighboring cells, C_i: all strongly coupled neighboring cells
         !!!-------------------------------------------------------------------------------------------
         CASE (NSCARC_CELLTYPE_SFINE)

            !!! Compute scaling factor: SCAL = [sum_(k in N_i) a_ik]/[sum_(l in C_i) a_il]
            SUM_COUPLED = 0.0_EB
            SUM_CPOINTS = 0.0_EB

            DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1
               JC = SC%A_COL(ICOL)
               SUM_COUPLED = SUM_COUPLED + SC%A(ICOL)
               IF (SC%CELLTYPE(JC) > 0) SUM_CPOINTS = SUM_CPOINTS + SC%A(ICOL)
            ENDDO

            SCAL = -SUM_COUPLED/SUM_CPOINTS

            !!! for each coupling of IC compute interpolation weight SCAL * a_ij/a_ii
            IW = 1
            DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1
               JC = SC%A_COL(ICOL)
               IF (SC%CELLTYPE(JC) > 0 ) THEN
                  NEIGHBOR(IW) = SC%CELLTYPE(JC)
                  WEIGHTS(IW)  = SCAL * SC%A(ICOL)/SC%A(IDIAG)
                  IW = IW +1
               ENDIF
            ENDDO
            NWEIGHTS = IW - 1

         !!!-------------------------------------------------------------------------------------------
         !!! Weakly coupled fine cell IC:
         !!! Determine strongly coupled fine cells JC surrounding IC and, in turn, replace each of them
         !!! by a mean value of their surrounding strongly coupled coarse cells
         !!!-------------------------------------------------------------------------------------------
         CASE (NSCARC_CELLTYPE_WFINE)


            IW = 1                                       ! weights counter
            DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1           ! loop over couplings of weakly coupled fine cell

               !!! get subdiagonal matrix entry a_ij of weakly coupled fine cell JC to IC
               JC = SC%A_COL(ICOL)

               !!! Find all surrounding (coupled) points of JC and compute scaling factor 
               !!! compute scaling factor SCAL = a_ij/a_ii * [sum_(k in N_j) a_jk]/[sum_(l in C_j) ajl] 
               SUM_COUPLED = 0.0_EB
               SUM_CPOINTS = 0.0_EB

               DO JCOL = SC%A_ROW(JC)+1, SC%A_ROW(JC+1)-1
                  KC = SC%A_COL(JCOL)
                  SUM_COUPLED = SUM_COUPLED + SC%A(JCOL)
                  IF (SC%CELLTYPE(KC) > 0 ) SUM_CPOINTS = SUM_CPOINTS + SC%A(JCOL)
               ENDDO

               SCAL = - SC%A(ICOL)/SC%A(IDIAG) * SUM_COUPLED/SUM_CPOINTS
                   
               !!! Get diagonal matrix a_jj for point JC
               JDIAG = SC%A_ROW(JC)

               !!! Compute interpolation weights for all strong coarse cell couplings KC of JC
               !!! note that a coarse cell KC may be considered several times for different JC's
               DO JCOL = SC%A_ROW(JC)+1, SC%A_ROW(JC+1)-1
                  KC = SC%A_COL(JCOL)
                  IF (SC%CELLTYPE(KC) > 0) THEN
                    NEIGHBOR(IW) = SC%CELLTYPE(KC)
                    WEIGHTS(IW)  = SCAL * SC%A(JCOL)/SC%A(JDIAG)
                    IW = IW +1
                  ENDIF
               ENDDO
                  
            ENDDO
            NWEIGHTS = IW - 1

            !!! make weights unique (add weights for multiple coarse cells)
            DO IW0 = 1, NWEIGHTS
               DO IW = IW0+1, NWEIGHTS
                  IF  (NEIGHBOR(IW0) /= -1 .AND. NEIGHBOR(IW) == NEIGHBOR(IW0)) THEN
                     WEIGHTS(IW0) = WEIGHTS(IW0) + WEIGHTS(IW)
                     WEIGHTS(IW)     = 0.0_EB
                     NEIGHBOR(IW) = -1
                  ENDIF
               ENDDO
            ENDDO

      END SELECT SELECT_FPOINT_TYPE

   ENDIF

   !!! -------------------------------------------------------------------------------------------------
   !!! Define corresponding entry in interpolation matrix P by means of the upper weights
   !!! -------------------------------------------------------------------------------------------------
   SC%P_ROW(IC) = IP
   DO IW = 1, NWEIGHTS
      IF  (NEIGHBOR(IW) /= -1) THEN
         SC%P_COL(IP) = NEIGHBOR(IW)
         SC%P(IP)    = WEIGHTS(IW)
         IP = IP +1
      ENDIF
   ENDDO

ENDDO
SC%P_ROW(SC%NC+1) = IP

!!! ----------------------------------------------------------------------------------------------------
!!! Define restriction matrix R (currently transpose of prolongation matrix P)
!!! In spite of the additinal need for the storing of R, this is done to save computational time
!!! in the later matrix transfer operations 
!!! ----------------------------------------------------------------------------------------------------
IC2 = 1
DO ICP = 1, SC%NCC

   BFIRST = .TRUE.
   DO IC = 1, SC%NC
      
      ROW_LOOP: DO ICP2 = SC%P_ROW(IC),SC%P_ROW(IC+1)-1
         IF (SC%P_COL(ICP2) == ICP) THEN
            SC%R(IC2)    = SC%P(ICP2)
            IF (BFIRST) SC%R_ROW(ICP) = IC2
            SC%R_COL(IC2) = IC
            IC2 = IC2 + 1
            BFIRST = .FALSE.
            EXIT ROW_LOOP
         ENDIF
      ENDDO ROW_LOOP

   ENDDO

ENDDO
SC%R_ROW(SC%NCC+1)=IC2

END SUBROUTINE SCARC_AMG_INTERPOLATION_DIRECT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Compute coarser matrix A_coarse by use of interpolation matrices:
!!!
!!!      A_coarse = I_fine^coarse * A_fine * I_coarse_fine
!!! 
!!! Note the different storage techniques (compact storage technique for 
!!! transfer matrices and coarse matrix, banded for finest matrix)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_AMG_TRANSFER (NM, NL)
INTEGER , INTENT(IN) :: NM, NL
INTEGER :: IC, ICO, ICP, IERR, ICP1, ICP2, IP, ICOL, IDIAG
REAL (EB), ALLOCATABLE, DIMENSION(:) :: VALUE
REAL:: AUX1, AUX2
TYPE (SCARC_SYSTEM_COMPACT_TYPE), POINTER :: SCF, SCC

IERR = 0

SCF => SCARC(NM)%COMPACT(NL)             ! pointer to fine level
SCC => SCARC(NM)%COMPACT(NL+1)           ! pointer to coarse level

ALLOCATE (VALUE(SCC%NC), STAT=IERR)
CALL CHKMEMERR ('SCARC_AMG_TRANSFER', 'VALUE', IERR)
VALUE = 0.0_EB

IP  = 1
ICP = 1
DO ICP1 = 1, SCC%NC
   DO ICP2 = 1, SCC%NC
      AUX2 = 0.0_EB
      ICP = (ICP1-1)*SCC%NC + ICP2
      DO IC = 1, SCF%NC

         IDIAG = SCF%A_ROW(IC)
         AUX1 = SCF%A(IDIAG) * CPOINT_WEIGHT(SCF%P, SCF%P_ROW, SCF%P_COL, IC, ICP2)

         COUPLINGS_LOOP: DO ICOL = SCF%A_ROW(IC)+1, SCF%A_ROW(IC+1)-1
            ICO = SCF%A_COL(ICOL)
            AUX1 = AUX1 +   SCF%A(ICOL)  * CPOINT_WEIGHT(SCF%P, SCF%P_ROW, SCF%P_COL, ICO, ICP2)
         ENDDO COUPLINGS_LOOP

         AUX2 = AUX2 + CPOINT_WEIGHT(SCF%P, SCF%P_ROW, SCF%P_COL, IC, ICP1) * AUX1

      ENDDO
      VALUE(ICP2) = AUX2
      ICP = ICP + 1
   ENDDO


   !!! analyze new matrix line and store it corresponding to compact storage technique:
   !!! (diagonal entry first)
   SCC%A(IP)      = VALUE(ICP1)
   SCC%A_ROW(ICP1) = IP
   SCC%A_COL(IP)   = ICP1

   IP  = IP + 1
   DO ICP2 = 1, SCC%NC
      IF (ICP2 /= ICP1 .AND. ABS(VALUE(ICP2)) >= 1.0E-12_EB) THEN
         SCC%A(IP)    = VALUE(ICP2)
         SCC%A_COL(IP) = ICP2
         IP  = IP + 1
      ENDIF
   ENDDO

ENDDO
SCC%A_ROW(SCC%NC+1) = IP
SCC%NA = IP - 1

DEALLOCATE(VALUE)

END SUBROUTINE SCARC_AMG_TRANSFER


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Determine if corresponding CPOINT contributes a non-zero interpolation weight or not
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(EB) FUNCTION CPOINT_WEIGHT(P, P_ROW, P_COL, IC, ICP)
REAL(EB), DIMENSION(:), INTENT(IN) :: P
INTEGER , DIMENSION(:), INTENT(IN) :: P_ROW, P_COL
INTEGER , INTENT(IN):: IC, ICP
INTEGER :: ICOL

CPOINT_WEIGHT = 0.0_EB
CPOINT_LOOP: DO ICOL = P_ROW(IC), P_ROW(IC+1)-1
   IF (P_COL(ICOL) == ICP) THEN
      CPOINT_WEIGHT =  P(ICOL)
      EXIT CPOINT_LOOP
   ENDIF
ENDDO CPOINT_LOOP

RETURN
END FUNCTION CPOINT_WEIGHT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Interface for the call of ScaRC-solver with requested storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SOLVER
REAL (EB) :: TNOW_SOLVER

TNOW_SOLVER = SECOND()

SELECT_METHOD: SELECT CASE (TYPE_METHOD)

   !!! Krylov method (CG/BICG)
   CASE (NSCARC_METHOD_KRYLOV)

      SELECT_KRYLOV: SELECT CASE (TYPE_KRYLOV)
         CASE (NSCARC_KRYLOV_CG)                           
            CALL SCARC_METHOD_CG  (NSCARC_SCOPE_MAIN, NSCARC_VECTOR_F)
         CASE (NSCARC_KRYLOV_BICG)                         
            CALL SCARC_METHOD_BICG(NSCARC_SCOPE_MAIN, NSCARC_VECTOR_F)
      END SELECT SELECT_KRYLOV

   !!! Multigrid method
   CASE (NSCARC_METHOD_MULTIGRID)

      CALL SCARC_METHOD_MULTIGRID(NSCARC_SCOPE_MAIN, NSCARC_VECTOR_F)

END SELECT SELECT_METHOD

TUSED_SCARC(NSCARC_TIME_SOLVER  ,:)=TUSED_SCARC(NSCARC_TIME_SOLVER  ,:)+SECOND()-TNOW_SOLVER
TUSED_SCARC(NSCARC_TIME_TOTAL,:)=TUSED_SCARC(NSCARC_TIME_TOTAL,:)+SECOND()-TNOW_SOLVER
END SUBROUTINE SCARC_SOLVER


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Set pointer to chosen vector for banded storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION POINT_TO_BANDED_VECTOR(NVECTOR, NM, NL)
REAL(EB), POINTER, DIMENSION(:,:,:) :: POINT_TO_BANDED_VECTOR
INTEGER, INTENT(IN):: NVECTOR, NM, NL

SELECT CASE (NVECTOR)
   CASE (NSCARC_VECTOR_X)
      POINT_TO_BANDED_VECTOR => SCARC(NM)%BANDED(NL)%X
   CASE (NSCARC_VECTOR_F)
      POINT_TO_BANDED_VECTOR => SCARC(NM)%BANDED(NL)%F
   CASE (NSCARC_VECTOR_Y)
      POINT_TO_BANDED_VECTOR => SCARC(NM)%BANDED(NL)%Y
   CASE (NSCARC_VECTOR_G)
      POINT_TO_BANDED_VECTOR => SCARC(NM)%BANDED(NL)%G
   CASE (NSCARC_VECTOR_W)
      POINT_TO_BANDED_VECTOR => SCARC(NM)%BANDED(NL)%W
   CASE (NSCARC_VECTOR_D)
      POINT_TO_BANDED_VECTOR => SCARC(NM)%BANDED(NL)%D
   CASE (NSCARC_VECTOR_Z)
      POINT_TO_BANDED_VECTOR => SCARC(NM)%BANDED(NL)%Z
   CASE (NSCARC_VECTOR_X2)
      POINT_TO_BANDED_VECTOR => SCARC(NM)%BANDED(NL)%X2
   CASE (NSCARC_VECTOR_D2)
      POINT_TO_BANDED_VECTOR => SCARC(NM)%BANDED(NL)%D2
   CASE (NSCARC_VECTOR_W2)
      POINT_TO_BANDED_VECTOR => SCARC(NM)%BANDED(NL)%W2
   CASE (NSCARC_VECTOR_Y2)
      POINT_TO_BANDED_VECTOR => SCARC(NM)%BANDED(NL)%Y2
END SELECT

RETURN
END FUNCTION POINT_TO_BANDED_VECTOR


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Set pointer to chosen vector for banded storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION POINT_TO_COMPACT_VECTOR(NVECTOR, NM, NL)
REAL(EB), POINTER, DIMENSION(:) :: POINT_TO_COMPACT_VECTOR
INTEGER, INTENT(IN):: NVECTOR, NM, NL

SELECT CASE (NVECTOR)
   CASE (NSCARC_VECTOR_X)
      POINT_TO_COMPACT_VECTOR => SCARC(NM)%COMPACT(NL)%X
   CASE (NSCARC_VECTOR_F)
      POINT_TO_COMPACT_VECTOR => SCARC(NM)%COMPACT(NL)%F
   CASE (NSCARC_VECTOR_Y)
      POINT_TO_COMPACT_VECTOR => SCARC(NM)%COMPACT(NL)%Y
   CASE (NSCARC_VECTOR_G)
      POINT_TO_COMPACT_VECTOR => SCARC(NM)%COMPACT(NL)%G
   CASE (NSCARC_VECTOR_W)
      POINT_TO_COMPACT_VECTOR => SCARC(NM)%COMPACT(NL)%W
   CASE (NSCARC_VECTOR_D)
      POINT_TO_COMPACT_VECTOR => SCARC(NM)%COMPACT(NL)%D
   CASE (NSCARC_VECTOR_Z)
      POINT_TO_COMPACT_VECTOR => SCARC(NM)%COMPACT(NL)%Z
   CASE (NSCARC_VECTOR_X2)
      POINT_TO_COMPACT_VECTOR => SCARC(NM)%COMPACT(NL)%X2
   CASE (NSCARC_VECTOR_D2)
      POINT_TO_COMPACT_VECTOR => SCARC(NM)%COMPACT(NL)%D2
   CASE (NSCARC_VECTOR_W2)
      POINT_TO_COMPACT_VECTOR => SCARC(NM)%COMPACT(NL)%W2
   CASE (NSCARC_VECTOR_Y2)
      POINT_TO_COMPACT_VECTOR => SCARC(NM)%COMPACT(NL)%Y2
END SELECT

RETURN
END FUNCTION POINT_TO_COMPACT_VECTOR


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Compute global matrix-vector product (including data exchange along internal boundaries)
!!! for banded storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_MATVEC_PRODUCT(NVECTOR1, NVECTOR2, NL)
INTEGER , INTENT(IN):: NVECTOR1, NVECTOR2, NL
REAL(EB), POINTER, DIMENSION(:,:,:) :: VB1, VB2
REAL(EB), POINTER, DIMENSION(:)     :: VC1, VC2
REAL(EB), POINTER, DIMENSION(:,:)   :: AB
REAL(EB), POINTER, DIMENSION(:)     :: AC
INTEGER , POINTER, DIMENSION(:)     :: AC_ROW, AC_COL
INTEGER , POINTER :: NX, NY, NZ, NC
INTEGER :: NM, I, J, K, IC, JC, ICOL

!!!----------------------------------------------------------------------------------------------------
!!! Exchange internal boundary values of vector1 such that the ghost values contain the corresponding
!!! overlapped values of adjacent neighbor
!!!----------------------------------------------------------------------------------------------------
IF (NMESHES > 1) THEN
   NREQ_SCARC = 0
   TYPE_EXCHANGE = NSCARC_EXCHANGE_VECTOR
   TYPE_VECTOR   = NVECTOR1
   CALL SCARC_RECEIVE  (NL)
   CALL SCARC_EXCHANGE (NL)
ENDIF

!!!----------------------------------------------------------------------------------------------------
!!! Perform global matrix-vector product:
!!! Note: - matrix already contains subdiagonal values from neighbor along internal boundaries
!!!       - if vector1 contains neighboring values, then correct values of global matvec are achieved
!!!----------------------------------------------------------------------------------------------------
SELECT CASE (TYPE_SYSTEM)

   !!! ---------------------------- Banded storage technique ------------------------------------------
   CASE (NSCARC_SYSTEM_BANDED)

      DO NM = NMESHES_MIN, NMESHES_MAX
      
         VB1 => POINT_TO_BANDED_VECTOR (NVECTOR1, NM, NL)
         VB2 => POINT_TO_BANDED_VECTOR (NVECTOR2, NM, NL)
         
         NX => SCARC(NM)%BANDED(NL)%NX                                     
         NY => SCARC(NM)%BANDED(NL)%NY
         NZ => SCARC(NM)%BANDED(NL)%NZ

         AB => SCARC(NM)%BANDED(NL)%A                                        ! system matrix
      
         SELECT CASE (TYPE_DIMENSION)
            
            CASE (NSCARC_DIMENSION_TWO)
      
               DO K = 1, NZ
                  DO I = 1, NX
                     IC = (K-1) * NX + I
                     VB2(I, 1, K) =   (  AB(IC, ID ) * VB1(I  , 1, K  )   &         ! diagonal component
                                       + AB(IC, ILZ) * VB1(I  , 1, K-1)   &         ! lower z-component
                                       + AB(IC, ILX) * VB1(I-1, 1, K  )   &         ! lower x-component
                                       + AB(IC, IUX) * VB1(I+1, 1, K  )   &         ! upper x-component
                                       + AB(IC, IUZ) * VB1(I  , 1, K+1) )           ! upper z-component
                  ENDDO
               ENDDO
         
            CASE (NSCARC_DIMENSION_THREE)
      
               DO K = 1, NZ
                  DO J = 1, NY
                     DO I = 1, NX
                        IC = (K-1) * NX * NY + (J-1) * NX + I
                        VB2(I, J, K) =  ( AB(IC, ID ) * VB1(I  , J  , K  )   &      ! diagonal component
                                      +   AB(IC, ILZ) * VB1(I  , J  , K-1)   &      ! lower z-component
                                      +   AB(IC, ILY) * VB1(I  , J-1, K  )   &      ! lower y-component
                                      +   AB(IC, ILX) * VB1(I-1, J  , K  )   &      ! lower x-component
                                      +   AB(IC, IUX) * VB1(I+1, J  , K  )   &      ! upper x-component
                                      +   AB(IC, IUY) * VB1(I  , J+1, K  )   &      ! upper y-component
                                      +   AB(IC, IUZ) * VB1(I  , J  , K+1) )        ! upper z-component
                     ENDDO
                  ENDDO
               ENDDO
         
         END SELECT
      
      ENDDO 
   
   !!! ---------------------------- Compact storage technique -----------------------------------------
   CASE (NSCARC_SYSTEM_COMPACT)

      DO NM = NMESHES_MIN, NMESHES_MAX
      
         VC1 => POINT_TO_COMPACT_VECTOR (NVECTOR1, NM, NL)
         VC2 => POINT_TO_COMPACT_VECTOR (NVECTOR2, NM, NL)
         
         NC     => SCARC(NM)%COMPACT(NL)%NC                                   ! number of cells

         AC     => SCARC(NM)%COMPACT(NL)%A                                    ! system matrix
         AC_ROW => SCARC(NM)%COMPACT(NL)%A_ROW                                ! row pointer
         AC_COL => SCARC(NM)%COMPACT(NL)%A_COL                                ! column pointer
      
         DO IC = 1, NC
                 
            !!! diagonal entry
            ICOL = AC_ROW(IC)                             
            JC   = AC_COL(ICOL)
      
            VC2 (IC) = AC(ICOL)* VC1(JC)
      
            !!! subdiagonal entries
            DO ICOL = AC_ROW(IC)+1, AC_ROW(IC+1)-1          
               JC = AC_COL(ICOL)
               VC2(IC) =  VC2(IC) + AC(ICOL)* VC1(JC)
            ENDDO
      
         ENDDO
      
      ENDDO 
      
END SELECT

END SUBROUTINE SCARC_MATVEC_PRODUCT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Compute global scalarproductt (including global data exchange) for banded storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(EB) FUNCTION SCARC_SCALAR_PRODUCT(NVECTOR1, NVECTOR2, NL)
INTEGER, INTENT(IN):: NVECTOR1, NVECTOR2, NL
REAL(EB), DIMENSION(:,:,:), POINTER ::  VB1, VB2
REAL(EB), DIMENSION(:)    , POINTER ::  VC1, VC2
INTEGER , POINTER :: NX, NY, NZ, NC
INTEGER  :: NM, IERR, NL0, I, J, K, IC
REAL(EB) :: SP_GLOBAL

!!!----------------------------------------------------------------------------------------------------
!!! Compute local scalarproduct
!!!----------------------------------------------------------------------------------------------------
SELECT CASE (TYPE_SYSTEM)

   !!! ---------------------------- Banded storage technique ------------------------------------------
   CASE (NSCARC_SYSTEM_BANDED)

      DO NM = NMESHES_MIN, NMESHES_MAX
      
         VB1 => POINT_TO_BANDED_VECTOR (NVECTOR1, NM, NL)
         VB2 => POINT_TO_BANDED_VECTOR (NVECTOR2, NM, NL)
         
         NX  => SCARC(NM)%BANDED(NL)%NX
         NY  => SCARC(NM)%BANDED(NL)%NY
         NZ  => SCARC(NM)%BANDED(NL)%NZ
      
         SP_LOCAL(NM) = 0.0_EB
         DO K = 1, NZ
            DO J = 1, NY
               DO I = 1, NX
                  SP_LOCAL(NM) = SP_LOCAL(NM) + VB1 (I, J, K) * VB2 (I, J, K)
               ENDDO
            ENDDO
         ENDDO
      
      ENDDO

   !!! ---------------------------- Compact storage technique -----------------------------------------
   CASE (NSCARC_SYSTEM_COMPACT)

      DO NM = NMESHES_MIN, NMESHES_MAX
      
         VC1 => POINT_TO_COMPACT_VECTOR (NVECTOR1, NM, NL)
         VC2 => POINT_TO_COMPACT_VECTOR (NVECTOR2, NM, NL)
         
         NC  => SCARC(NM)%COMPACT(NL)%NC
      
         SP_LOCAL(NM) = 0.0_EB
         DO IC = 1, NC
            SP_LOCAL(NM) = SP_LOCAL(NM) + VC1(IC) * VC2(IC)
         ENDDO
      
      ENDDO

END SELECT

!!!----------------------------------------------------------------------------------------------------
!!! get global scalarproduct by a global summation of the local values
!!!----------------------------------------------------------------------------------------------------
IERR = 0
NL0  = NL
SP_GLOBAL   = 0.0_EB
IF (NMESHES>1 .AND. USE_MPI) THEN
   CALL MPI_ALLREDUCE (SP_LOCAL, SP_GLOBAL, NMESHES, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
ELSE
   DO NM=1,NMESHES
      SP_GLOBAL = SP_GLOBAL + SP_LOCAL(NM)
   ENDDO
ENDIF

SCARC_SCALAR_PRODUCT = SP_GLOBAL
RETURN

END FUNCTION SCARC_SCALAR_PRODUCT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Compute global L2-norm (including global data exchange) for banded storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(EB) FUNCTION SCARC_L2NORM(NVECTOR1, NL)
INTEGER, INTENT(IN):: NVECTOR1, NL
REAL(EB), DIMENSION(:,:,:), POINTER ::  VB
REAL(EB), DIMENSION(:)    , POINTER ::  VC
INTEGER , POINTER :: NX, NY, NZ, NC
INTEGER  :: NM, IERR, I, J, K, IC
REAL(EB) :: SP_GLOBAL

!!!----------------------------------------------------------------------------------------------------
!!! Compute local scalarproduct
!!!----------------------------------------------------------------------------------------------------
SELECT CASE (TYPE_SYSTEM)

   !!! ---------------------------- Banded storage technique ------------------------------------------
   CASE (NSCARC_SYSTEM_BANDED)

      DO NM = NMESHES_MIN, NMESHES_MAX
      
         VB => POINT_TO_BANDED_VECTOR (NVECTOR1, NM, NL)

         NX => SCARC(NM)%BANDED(NL)%NX
         NY => SCARC(NM)%BANDED(NL)%NY
         NZ => SCARC(NM)%BANDED(NL)%NZ
      
         SP_LOCAL(NM) = 0.0_EB
         DO K = 1, NZ
            DO J = 1, NY
               DO I = 1, NX
                  SP_LOCAL(NM) = SP_LOCAL(NM) + VB (I, J, K) * VB (I, J, K)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      
   !!! ---------------------------- Compact storage technique -----------------------------------------
   CASE (NSCARC_SYSTEM_COMPACT)

      DO NM = NMESHES_MIN, NMESHES_MAX
      
         VC => POINT_TO_COMPACT_VECTOR (NVECTOR1, NM, NL)
      
         NC => SCARC(NM)%COMPACT(NL)%NC
      
         SP_LOCAL(NM) = 0.0_EB
         DO IC = 1, NC
            SP_LOCAL(NM) = SP_LOCAL(NM) + VC(IC) * VC(IC)
         ENDDO
      
      ENDDO

END SELECT

!!!----------------------------------------------------------------------------------------------------
!!! get global scalarproduct by a global summation of the local values
!!!----------------------------------------------------------------------------------------------------
IERR = 0
SP_GLOBAL = 0.0_EB
IF (NMESHES>1 .AND. USE_MPI) THEN
   CALL MPI_ALLREDUCE (SP_LOCAL, SP_GLOBAL, NMESHES, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
ELSE
   DO NM=1,NMESHES
      SP_GLOBAL = SP_GLOBAL + SP_LOCAL(NM)
   ENDDO
ENDIF
SP_GLOBAL = SQRT (SP_GLOBAL/REAL(NC_GLOBAL(NL), EB))   

SCARC_L2NORM = SP_GLOBAL
RETURN

END FUNCTION SCARC_L2NORM



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Compute linear combination of two vectors for banded storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_LINEAR_COMBI(NVECTOR1, NVECTOR2, SCAL1, SCAL2, NL)
INTEGER , INTENT(IN):: NVECTOR1, NVECTOR2, NL
REAL(EB), INTENT(IN):: SCAL1, SCAL2
REAL(EB), DIMENSION(:,:,:), POINTER ::  VB1, VB2
REAL(EB), DIMENSION(:)    , POINTER ::  VC1, VC2
INTEGER  :: NM

SELECT CASE (TYPE_SYSTEM)

   !!! ---------------------------- Banded storage technique ------------------------------------------
   CASE (NSCARC_SYSTEM_BANDED)

      DO NM = NMESHES_MIN, NMESHES_MAX

         VB1 => POINT_TO_BANDED_VECTOR(NVECTOR1, NM, NL)
         VB2 => POINT_TO_BANDED_VECTOR(NVECTOR2, NM, NL)

         VB2 = SCAL1 * VB1 + SCAL2 * VB2

      ENDDO

   !!! ---------------------------- Compact storage technique ------------------------------------------
   CASE (NSCARC_SYSTEM_COMPACT)

      DO NM = NMESHES_MIN, NMESHES_MAX

         VC1 => POINT_TO_COMPACT_VECTOR(NVECTOR1, NM, NL)
         VC2 => POINT_TO_COMPACT_VECTOR(NVECTOR2, NM, NL)

         VC2 = SCAL1 * VC1 + SCAL2 * VC2

      ENDDO

END SELECT

END SUBROUTINE SCARC_LINEAR_COMBI



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Define vector2 to be a scaled copy of vector 1 for banded storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SCALED_COPY(NVECTOR1, NVECTOR2, SCAL1, NL)
INTEGER , INTENT(IN):: NVECTOR1, NVECTOR2, NL
REAL(EB), INTENT(IN):: SCAL1
REAL(EB), DIMENSION(:,:,:), POINTER ::  VB1, VB2
REAL(EB), DIMENSION(:)    , POINTER ::  VC1, VC2
INTEGER  :: NM

SELECT CASE (TYPE_SYSTEM)

   !!! ---------------------------- Banded storage technique ------------------------------------------
   CASE (NSCARC_SYSTEM_BANDED)

      DO NM = NMESHES_MIN, NMESHES_MAX

         VB1 => POINT_TO_BANDED_VECTOR(NVECTOR1, NM, NL)
         VB2 => POINT_TO_BANDED_VECTOR(NVECTOR2, NM, NL)

         VB2 = SCAL1 * VB1 

      ENDDO

   !!! ---------------------------- Compact storage technique ------------------------------------------
   CASE (NSCARC_SYSTEM_COMPACT)

      DO NM = NMESHES_MIN, NMESHES_MAX

         VC1 => POINT_TO_COMPACT_VECTOR(NVECTOR1, NM, NL)
         VC2 => POINT_TO_COMPACT_VECTOR(NVECTOR2, NM, NL)

         VC2 = SCAL1 * VC1 
      
      ENDDO

END SELECT

END SUBROUTINE SCARC_SCALED_COPY


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Clear vector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_CLEAR(NVECTOR, NL)
INTEGER , INTENT(IN):: NVECTOR, NL
REAL(EB), DIMENSION(:,:,:), POINTER ::  VB
REAL(EB), DIMENSION(:)    , POINTER ::  VC
INTEGER  :: NM

SELECT CASE (TYPE_SYSTEM)

   !!! ---------------------------- Banded storage technique ------------------------------------------
   CASE (NSCARC_SYSTEM_BANDED)

      DO NM = NMESHES_MIN, NMESHES_MAX
         VB => POINT_TO_BANDED_VECTOR(NVECTOR, NM, NL)
         VB =  0.0_EB
      ENDDO

   !!! ---------------------------- Compact storage technique ------------------------------------------
   CASE (NSCARC_SYSTEM_COMPACT)

      DO NM = NMESHES_MIN, NMESHES_MAX
         VC => POINT_TO_COMPACT_VECTOR(NVECTOR, NM, NL)
         VC =  0.0_EB
      ENDDO

END SELECT

END SUBROUTINE SCARC_CLEAR


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform preconditioning for banded storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PRECONDITIONER (NVECTOR1, NVECTOR2, NL)
USE POIS, ONLY: H2CZSS, H3CZSS
INTEGER, INTENT(IN):: NVECTOR1, NVECTOR2, NL
REAL(EB), DIMENSION(:,:,:), POINTER ::  VB1, VB2, FFT
REAL(EB), DIMENSION(:)    , POINTER ::  VC1, VC2
REAL(EB), DIMENSION(:,:)  , POINTER ::  AB
REAL(EB), DIMENSION(:)    , POINTER ::  AC
INTEGER , DIMENSION(:)    , POINTER ::  AC_ROW, AC_COL
INTEGER , POINTER:: NX, NY, NZ, NC
INTEGER  :: NM, I, J, K, IC, ICOL
REAL(EB) :: AUX, OMEGA=1.5_EB
TYPE (MESH_TYPE), POINTER :: M

SELECT CASE (TYPE_SYSTEM)

   !!! ---------------------------- Banded storage technique ------------------------------------------
   CASE (NSCARC_SYSTEM_BANDED)

      SELECT CASE (TYPE_PRECON)
      
         !!!--------------------------------------------------------------------------------------------
         !!! Multigrid preconditioner
         !!!--------------------------------------------------------------------------------------------
         CASE (NSCARC_PRECON_MG)
      
            CALL SCARC_METHOD_MULTIGRID (NSCARC_SCOPE_PRECON, NVECTOR2)
      
         !!!--------------------------------------------------------------------------------------------
         !!! Jacobi preconditioner
         !!!--------------------------------------------------------------------------------------------
         CASE (NSCARC_PRECON_JACOBI)
      
            DO NM = NMESHES_MIN, NMESHES_MAX
      
               VB2 => POINT_TO_BANDED_VECTOR(NVECTOR2, NM, NL)
         
               NX => SCARC(NM)%BANDED(NL)%NX
               NY => SCARC(NM)%BANDED(NL)%NY
               NZ => SCARC(NM)%BANDED(NL)%NZ
               AB => SCARC(NM)%BANDED(NL)%A
      
               DO K = 1, NZ
                  DO J = 1, NY
                     DO I = 1, NX
                        IC = (K-1) * NX * NY + (J-1) * NX + I
                        VB2 (I, J, K) = VB2 (I, J, K) / AB(IC, ID)
                     ENDDO
                  ENDDO
               ENDDO
      
            ENDDO
         
         !!!--------------------------------------------------------------------------------------------
         !!! SSOR preconditioner
         !!!--------------------------------------------------------------------------------------------
         CASE (NSCARC_PRECON_SSOR)
         
            DO NM = NMESHES_MIN, NMESHES_MAX
            
               VB2 => POINT_TO_BANDED_VECTOR(NVECTOR2, NM, NL)
         
               NX => SCARC(NM)%BANDED(NL)%NX
               NY => SCARC(NM)%BANDED(NL)%NY
               NZ => SCARC(NM)%BANDED(NL)%NZ
               AB => SCARC(NM)%BANDED(NL)%A
         
               SELECT CASE (TYPE_DIMENSION)
               
                  CASE (NSCARC_DIMENSION_TWO)
               
                     DO K = 1, NZ
                        DO I = 1, NX
                           IC = (K-1) * NX + I
                           AUX =    AB(IC,ILZ) * VB2 (I  , 1, K-1)  &
                                  + AB(IC,ILX) * VB2 (I-1, 1, K  )
                           VB2 (I, 1, K) = (VB2(I, 1, K) - AUX*OMEGA) / AB(IC,ID)
                        ENDDO
                     ENDDO
                     DO K = NZ, 1, - 1
                        DO I = NX, 1, - 1
                           IC = (K-1) * NX + I
                           AUX =    AB(IC,IUZ) * VB2 (I  , 1, K+1) &
                                  + AB(IC,IUX) * VB2 (I+1, 1, K  )
                           VB2 (I, 1, K) = VB2 (I, 1, K) - AUX * OMEGA / AB(IC,ID)
                        ENDDO
                     ENDDO
               
                  CASE (NSCARC_DIMENSION_THREE)
               
                     DO K = 1, NZ
                        DO J = 1, NY
                           DO I = 1, NX
                              IC = (K-1) * NX * NY + (J-1) * NX + I
                              AUX =    AB(IC,ILZ) * VB2 (I  , J  , K-1) &
                                     + AB(IC,ILY) * VB2 (I  , J-1, K  ) &
                                     + AB(IC,ILX) * VB2 (I-1, J  , K  )
                              VB2 (I, J, K) = (VB2(I, J, K) - AUX * OMEGA) / AB(IC,ID)
                           ENDDO
                        ENDDO
                     ENDDO
                     DO K = NZ, 1, - 1
                        DO J = NY, 1, - 1
                           DO I = NX, 1, - 1
                              IC = (K-1) * NX * NY + (J-1) * NX + I
                              AUX =    AB(IC,IUZ) * VB2 (I  , J  , K+1) &
                                     + AB(IC,IUY) * VB2 (I  , J+1, K  ) &
                                     + AB(IC,IUZ) * VB2 (I+1, J  , K  )
                              VB2 (I, J, K) = VB2 (I, J, K) - AUX * OMEGA / AB(IC,ID)
                           ENDDO
                        ENDDO
                     ENDDO
               END SELECT
         
            ENDDO
         
         !!!--------------------------------------------------------------------------------------------
         !!! GSTRIX preconditioner
         !!!--------------------------------------------------------------------------------------------
         !CASE (NSCARC_PRECON_GSTRIX)
         !   CALL SCARC_GSTRIX_BANDED(V2, ...)
            
         !!!--------------------------------------------------------------------------------------------
         !!! FFT preconditioner
         !!!--------------------------------------------------------------------------------------------
         CASE (NSCARC_PRECON_FFT)
            
            DO NM = NMESHES_MIN, NMESHES_MAX
            
               VB1 => POINT_TO_BANDED_VECTOR(NVECTOR1, NM, NL)
               VB2 => POINT_TO_BANDED_VECTOR(NVECTOR2, NM, NL)
         
               M   => MESHES(NM)
               FFT => SCARC(NM)%PRECON(NLEVEL_MIN)%FFT
               
               FFT(1:M%IBAR, 1:M%JBAR, 1:M%KBAR) = VB1(1:M%IBAR, 1:M%JBAR, 1:M%KBAR)
               SELECT CASE (TYPE_DIMENSION)
                  CASE (NSCARC_DIMENSION_TWO)
                     CALL H2CZSS (M%BXS, M%BXF, M%BZS, M%BZF, &
                                  M%IBAR+1, FFT, M%POIS_PTB, M%SAVE1, M%WORK, M%HX)
                  CASE (NSCARC_DIMENSION_THREE)
                     CALL H3CZSS (M%BXS, M%BXF, M%BYS, M%BYF, M%BZS, M%BZF, &
                                  M%IBAR+1, M%JBAR+1, FFT, M%POIS_PTB, M%SAVE1, M%WORK, M%HX)
               END SELECT
               VB2(1:M%IBAR, 1:M%JBAR, 1:M%KBAR) = FFT(1:M%IBAR, 1:M%JBAR, 1:M%KBAR)
               
            ENDDO
                     
      END SELECT
      

   !!! ---------------------------- Compact storage technique ------------------------------------------
   CASE (NSCARC_SYSTEM_COMPACT)

      SELECT CASE (TYPE_PRECON)
      
         !!!--------------------------------------------------------------------------------------------
         !!! Multigrid preconditioner
         !!!--------------------------------------------------------------------------------------------
         CASE (NSCARC_PRECON_MG)
      
            CALL SCARC_METHOD_MULTIGRID (NSCARC_SCOPE_PRECON, NVECTOR2)
      
      
         !!!--------------------------------------------------------------------------------------------
         !!! Jacobi preconditioner
         !!!--------------------------------------------------------------------------------------------
         CASE (NSCARC_PRECON_JACOBI)
      
            DO NM = NMESHES_MIN, NMESHES_MAX
      
               VC2 => POINT_TO_COMPACT_VECTOR(NVECTOR2, NM, NL)
         
               NC     => SCARC(NM)%COMPACT(NL)%NC
               AC     => SCARC(NM)%COMPACT(NL)%A
               AC_ROW => SCARC(NM)%COMPACT(NL)%A_ROW
      
               DO IC = 1, NC
                  VC2 (IC) = VC2 (IC) / AC (AC_ROW(IC))
               ENDDO
      
            ENDDO
         
         !!!--------------------------------------------------------------------------------------------
         !!! SSOR preconditioner
         !!!--------------------------------------------------------------------------------------------
         CASE (NSCARC_PRECON_SSOR)
         
            DO NM = NMESHES_MIN, NMESHES_MAX
            
               VC2 => POINT_TO_COMPACT_VECTOR(NVECTOR2, NM, NL)
         
               NC     => SCARC(NM)%COMPACT(NL)%NC
               AC     => SCARC(NM)%COMPACT(NL)%A
               AC_ROW => SCARC(NM)%COMPACT(NL)%A_ROW
               AC_COL => SCARC(NM)%COMPACT(NL)%A_COL
         
               !!! use only matrix superdiagonals
               FORWARD_CELL_LOOP: DO IC = 1, NC
               
                  AUX = 0.0_EB
                  LOWER_DIAG_LOOP: DO ICOL = AC_ROW(IC)+1, AC_ROW(IC+1)-1
                     IF (AC_COL(ICOL) >= IC) EXIT LOWER_DIAG_LOOP
                     AUX = AUX + AC(ICOL) * VC2(AC_COL(ICOL))
                  ENDDO LOWER_DIAG_LOOP
                
                  VC2(IC) = (VC2(IC) - AUX * OMEGA) / AC(AC_ROW(IC))
               
               ENDDO FORWARD_CELL_LOOP
               
               !!! use only matrix subdiagonals
               BACKWARD_CELL_LOOP: DO IC = NC-1, 1, -1
                  
                  AUX = 0.0_EB
                  UPPER_DIAG_LOOP: DO ICOL = AC_ROW(IC)+1, AC_ROW(IC+1)-1
                     IF (AC_COL(ICOL) <= IC) CYCLE
                     AUX = AUX + AC(ICOL) * VC2(AC_COL(ICOL))
                  ENDDO UPPER_DIAG_LOOP
                  
                  VC2(IC) = VC2(IC) - AUX * OMEGA / AC(AC_ROW(IC))
               
               ENDDO BACKWARD_CELL_LOOP
         
            ENDDO
            
         !!!--------------------------------------------------------------------------------------------
         !!! GSTRIX preconditioner
         !!!--------------------------------------------------------------------------------------------
         !CASE (NSCARC_PRECON_GSTRIX)
         !   CALL SCARC_GSTRIX_BANDED(V2, ...)
            
         !!!--------------------------------------------------------------------------------------------
         !!! FFT preconditioner
         !!!--------------------------------------------------------------------------------------------
         CASE (NSCARC_PRECON_FFT)
            
            DO NM = NMESHES_MIN, NMESHES_MAX
            
               VC1 => POINT_TO_COMPACT_VECTOR(NVECTOR1, NM, NL)
               VC2 => POINT_TO_COMPACT_VECTOR(NVECTOR2, NM, NL)
         
               M   => MESHES(NM)
               FFT => SCARC(NM)%PRECON(NLEVEL_MIN)%FFT
               
               DO K = 1, M%KBAR
                  DO J = 1, M%JBAR
                     DO I = 1, M%IBAR
                        IC = (K-1) * M%IBAR * M%JBAR + (J-1) * M%IBAR + I
                        FFT(I, J, K) = VC1(IC)
                     ENDDO
                  ENDDO
               ENDDO
               SELECT CASE (TYPE_DIMENSION)
                  CASE (NSCARC_DIMENSION_TWO)
                     CALL H2CZSS (M%BXS, M%BXF, M%BZS, M%BZF,&
                                  M%IBAR+1, FFT, M%POIS_PTB, M%SAVE1, M%WORK, M%HX)
                  CASE (NSCARC_DIMENSION_THREE)
                     CALL H3CZSS (M%BXS, M%BXF, M%BYS, M%BYF, M%BZS, M%BZF, &
                                  M%IBAR+1, M%JBAR+1, FFT, M%POIS_PTB, M%SAVE1, M%WORK, M%HX)
               END SELECT
               DO K = 1, M%KBAR
                  DO J = 1, M%JBAR
                     DO I = 1, M%IBAR
                        IC = (K-1) * M%IBAR * M%JBAR + (J-1) * M%IBAR + I
                        VC2(IC) = FFT(I, J, K)
                     ENDDO
                  ENDDO
               ENDDO
      
            ENDDO
      
      END SELECT

END SELECT

END SUBROUTINE SCARC_PRECONDITIONER 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform global CG-method based on global possion-matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_METHOD_CG(NSCOPE, NVECTOR)
INTEGER, INTENT(IN) :: NSCOPE, NVECTOR
INTEGER   :: NL = NSCARC_LEVEL_NONE
INTEGER   :: ITE, NIT, ISTATE
REAL (EB) :: SIGMA0, SIGMA1, ALPHA0, GAMMA0
REAL (EB) :: RES, RESIN, EPS
REAL (EB) :: TNOW_KRYLOV
CHARACTER(30) :: CROUTINE = 'null'

!!!----------------------------------------------------------------------------------------------------
!!! Initialization:
!!!   - Set environment variables and define working level
!!!   - Initialize solution, right hand side vector and auxiliary vectors
!!!   - Define iterations parameters
!!!----------------------------------------------------------------------------------------------------
TNOW_KRYLOV = SECOND()
TYPE_SCOPE  = NSCOPE

EPS = SCARC_KRYLOV_ACCURACY
NIT = SCARC_KRYLOV_ITERATIONS

CALL SCARC_SETUP_ENVIRONMENT(CROUTINE, NSCOPE, NVECTOR, NL)
CALL SCARC_SETUP_SOLVER(NL)

!!!----------------------------------------------------------------------------------------------------
!!! Compute initial residual and perform initial preconditioning
!!!----------------------------------------------------------------------------------------------------
CALL SCARC_MATVEC_PRODUCT (VEC_X, VEC_W, NL)                                !  W := A*X
CALL SCARC_LINEAR_COMBI   (VEC_F, VEC_W, -1.0_EB, 1.0_EB, NL)               !  W := W - F

RESIN = SCARC_L2NORM (VEC_W, NL)                                            !  RESIN := ||W||

CALL SCARC_CONVERGENCE_INFO(RESIN, 0, NL, CROUTINE)
CALL SCARC_SCALED_COPY   (VEC_W, VEC_G, 1.0_EB, NL)                         !  G := W
CALL SCARC_PRECONDITIONER(VEC_W, VEC_G, NL)                                 !  G := PRECON(W)

SIGMA0 = SCARC_SCALAR_PRODUCT(VEC_W, VEC_G, NL)                             !  SIGMA0 := (W,G)

CALL SCARC_SCALED_COPY(VEC_G, VEC_D, -1.0_EB, NL)                           !  D := -G

!!!----------------------------------------------------------------------------------------------------
!!! Perform conjugate gradient looping
!!!----------------------------------------------------------------------------------------------------
CG_LOOP: DO ITE = 1, NIT
 
   CALL SCARC_MATVEC_PRODUCT (VEC_D, VEC_Y, NL)                              !  Y := A*D

   ALPHA0 = SCARC_SCALAR_PRODUCT (VEC_D, VEC_Y, NL)                          !  ALPHA0 := (D,Y)
   ALPHA0 = SIGMA0/ALPHA0                                                 

   CALL SCARC_LINEAR_COMBI (VEC_D, VEC_X, ALPHA0, 1.0_EB, NL)                !  X := ALPHA0*D + X
   CALL SCARC_LINEAR_COMBI (VEC_Y, VEC_W, ALPHA0, 1.0_EB, NL)                !  W := ALPHA0*Y + W

   RES = SCARC_L2NORM (VEC_W, NL)                                            !  RES := ||W||
 
   ISTATE = SCARC_CONVERGENCE_STATE (RESIN, RES, EPS, ITE, NL, CROUTINE)     !  RES < TOL ??
   IF (ISTATE /= NSCARC_STATE_PROCEED) EXIT CG_LOOP
 
   CALL SCARC_SCALED_COPY (VEC_W, VEC_G, 1.0_EB, NL)                         !  G := W
   CALL SCARC_PRECONDITIONER (VEC_W, VEC_G, NL)                              !  G := PRECON(W)

   SIGMA1 = SCARC_SCALAR_PRODUCT (VEC_W, VEC_G, NL)                          !  SIGMA1 := (W,G)
   GAMMA0 = SIGMA1/SIGMA0                                                   
   SIGMA0 = SIGMA1                                                         

   CALL SCARC_LINEAR_COMBI (VEC_G, VEC_D, -1.0_EB, GAMMA0, NL)               !  D := -G + GAMMA0*D

ENDDO CG_LOOP
 
!!!----------------------------------------------------------------------------------------------------
!!! Determine convergence rate and print corresponding information
!!! In case of CG as main solver:
!!!   - Transfer ScaRC solution vector X to FDS pressure vector 
!!!   - Set ghost cell values along external boundaries
!!!   - Exchange values along internal boundaries 
!!!----------------------------------------------------------------------------------------------------
CALL SCARC_CONVERGENCE_RATE(RESIN, RES, ITE, ISTATE, CROUTINE)

IF (TYPE_SCOPE == NSCARC_SCOPE_MAIN) THEN
   CALL SCARC_TERMINATE_SOLVER(NLEVEL_MIN)
   CALL SCARC_UPDATE_GHOSTCELLS(NLEVEL_MIN)
ENDIF

TUSED_SCARC(NSCARC_TIME_KRYLOV,:)=TUSED_SCARC(NSCARC_TIME_KRYLOV,:)+SECOND()-TNOW_KRYLOV
TUSED_SCARC(NSCARC_TIME_TOTAL ,:)=TUSED_SCARC(NSCARC_TIME_TOTAL ,:)+SECOND()-TNOW_KRYLOV
END SUBROUTINE SCARC_METHOD_CG


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform global BICGstab-method based on global possion-matrix - banded storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_METHOD_BICG(NSCOPE, NVECTOR)
INTEGER, INTENT(IN) :: NSCOPE, NVECTOR
INTEGER   :: NL = NSCARC_LEVEL_NONE
INTEGER   :: ITE, NIT, ISTATE
REAL (EB) :: ALPHA0, ALPHA1, ALPHA2, RHO0, RHO1, DTHETA, DBETA
REAL (EB) :: RES, RESIN, EPS
REAL (EB) :: TNOW_KRYLOV
CHARACTER(30) :: CROUTINE = 'null'

!!!----------------------------------------------------------------------------------------------------
!!! Initialization
!!!   - Set environment variables and define working level
!!!   - Initialize solution, right hand side vector and auxiliary vectors
!!!   - Define iterations parameters
!!!----------------------------------------------------------------------------------------------------
TNOW_KRYLOV = SECOND()
TYPE_SCOPE  = NSCOPE

EPS    = SCARC_KRYLOV_ACCURACY
NIT    = SCARC_KRYLOV_ITERATIONS

ALPHA0 = 1.0_EB
RHO0   = 1.0_EB
RHO1   = 0.0_EB
DBETA  = 0.0_EB
DTHETA = 1.0_EB
 
CALL SCARC_SETUP_ENVIRONMENT(CROUTINE, NSCOPE, NVECTOR, NL)
CALL SCARC_SETUP_SOLVER(NL)
   
!!!----------------------------------------------------------------------------------------------------
!!! Compute initial defect and perform (double) initial preconditioning
!!!----------------------------------------------------------------------------------------------------
!CALL SCARC_SCALED_COPY    (VEC_F, VEC_W, 1.0_EB, NL)                         !  W := F
!CALL SCARC_PRECONDITIONER (VEC_W, VEC_W, NL)                                 !  W := PRECON(W)
CALL SCARC_MATVEC_PRODUCT (VEC_X, VEC_W, NL)                                 !  W := A*X
CALL SCARC_LINEAR_COMBI   (VEC_F, VEC_W, 1.0_EB, -1.0_EB, NL)                !  W := F - W
CALL SCARC_PRECONDITIONER (VEC_W, VEC_W, NL)                                 !  W := PRECON(W)

RESIN = SCARC_L2NORM (VEC_W, NL)                                             !  RESIN := ||W||

CALL SCARC_CONVERGENCE_INFO (RESIN, 0, NL, CROUTINE)
CALL SCARC_SCALED_COPY (VEC_W, VEC_G, 1.0_EB, NL)                            !  G := W
   
!!!----------------------------------------------------------------------------------------------------
!!! Perform bi-conjugate gradient looping:
!!!----------------------------------------------------------------------------------------------------
BICG_LOOP: DO ITE = 1, NIT

   RHO1  = SCARC_SCALAR_PRODUCT (VEC_G, VEC_W, NL)                           ! RHO1 := (G,W)
   DBETA = (RHO1*DTHETA)/(RHO0*ALPHA0)                                     
   RHO0  = RHO1

   CALL SCARC_LINEAR_COMBI   (VEC_W, VEC_Z, 1.0_EB       , DBETA , NL)       ! Z := W + DBETA*Z
   CALL SCARC_LINEAR_COMBI   (VEC_Y, VEC_Z, -DBETA*ALPHA0, 1.0_EB, NL)       ! Z := -DBETA*ALPHA0*Y + Z
   CALL SCARC_MATVEC_PRODUCT (VEC_Z, VEC_Y, NL)                              ! Y := A*Z
   CALL SCARC_PRECONDITIONER (VEC_Y, VEC_Y, NL)                              ! Z := PRECON(Z)

   DTHETA = SCARC_SCALAR_PRODUCT (VEC_G, VEC_Y, NL)                          ! DTHETA := (G,Y)
   DTHETA = RHO1/DTHETA

   CALL SCARC_LINEAR_COMBI   (VEC_Y, VEC_W, -DTHETA, 1.0_EB, NL)             ! W := -DTHETA*Y + W
   CALL SCARC_MATVEC_PRODUCT (VEC_W, VEC_D, NL)                              ! D := A*W
   CALL SCARC_PRECONDITIONER (VEC_D, VEC_D, NL)                              ! D := PRECON(D)
   
   ALPHA1 = SCARC_SCALAR_PRODUCT (VEC_D, VEC_W, NL)                          ! ALPHA1 := (D,W)
   ALPHA2 = SCARC_SCALAR_PRODUCT (VEC_D, VEC_D, NL)                          ! ALPHA2 := (D,D)
   ALPHA0 = ALPHA1/ALPHA2

   CALL SCARC_LINEAR_COMBI (VEC_Z, VEC_X,  DTHETA, 1.0_EB, NL)               ! X :=  DTHETA*Z + X
   CALL SCARC_LINEAR_COMBI (VEC_W, VEC_X,  ALPHA0, 1.0_EB, NL)               ! X :=  ALPHA0*W + X
   CALL SCARC_LINEAR_COMBI (VEC_D, VEC_W, -ALPHA0, 1.0_EB, NL)               ! W := -ALPHA0*D + W

   RES = SCARC_L2NORM (VEC_W, NL)                                            ! RES := ||W||

   ISTATE = SCARC_CONVERGENCE_STATE(RESIN, RES, EPS, ITE, NL, CROUTINE)      ! RES < TOL ???
   IF (ISTATE /= NSCARC_STATE_PROCEED) EXIT BICG_LOOP
 
ENDDO BICG_LOOP

!!!----------------------------------------------------------------------------------------------------
!!! Determine convergence rate and print corresponding information
!!! In case of BICG as main solver:
!!!   - Transfer ScaRC solution vector X to FDS pressure vector 
!!!   - Set ghost cell values along external boundaries
!!!   - Exchange values along internal boundaries 
!!!----------------------------------------------------------------------------------------------------
CALL SCARC_CONVERGENCE_RATE(RESIN, RES, ITE, ISTATE, CROUTINE)

IF (TYPE_SCOPE == NSCARC_SCOPE_MAIN) THEN
   CALL SCARC_TERMINATE_SOLVER(NLEVEL_MIN)
   CALL SCARC_UPDATE_GHOSTCELLS(NLEVEL_MIN)
ENDIF

TUSED_SCARC(NSCARC_TIME_KRYLOV,:)=TUSED_SCARC(NSCARC_TIME_KRYLOV,:)+SECOND()-TNOW_KRYLOV
TUSED_SCARC(NSCARC_TIME_TOTAL ,:)=TUSED_SCARC(NSCARC_TIME_TOTAL ,:)+SECOND()-TNOW_KRYLOV
END SUBROUTINE SCARC_METHOD_BICG


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform geometric multigrid method based on global possion-matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_METHOD_MULTIGRID(NSCOPE, NVECTOR)
INTEGER, INTENT(IN) :: NSCOPE, NVECTOR
INTEGER   :: NL = NSCARC_LEVEL_NONE
INTEGER   :: ITE, NIT, ISTATE, ICYCLE
REAL (EB) :: RES, RESIN, EPS
REAL (EB) :: TNOW_MULTIGRID
CHARACTER(30):: CROUTINE='null'

!!!----------------------------------------------------------------------------------------------------
!!! General initialization:
!!!   - Set environment variables and define working level
!!!   - Initialize solution, right hand side vector 
!!!   - Define iterations parameters (NL is set to finest level)
!!!----------------------------------------------------------------------------------------------------
TNOW_MULTIGRID = SECOND()
TYPE_SCOPE     = NSCOPE

EPS = SCARC_MULTIGRID_ACCURACY
NIT = SCARC_MULTIGRID_ITERATIONS

CALL SCARC_SETUP_ENVIRONMENT(CROUTINE, NSCOPE, NVECTOR, NL)
CALL SCARC_SETUP_SOLVER(NL)

!!!----------------------------------------------------------------------------------------------------
!!! Compute initial defect:  RESIN := || F - A*X ||
!!!   - Initialize cycle counts for MG-iteration
!!!   - Perform initial matrix-vector product on finest level
!!!   - calculate norm of initial residual on finest level
!!!----------------------------------------------------------------------------------------------------
CALL SCARC_MATVEC_PRODUCT (VEC_X, VEC_D, NL)                                     !  D := A*X
CALL SCARC_LINEAR_COMBI   (VEC_F, VEC_D, 1.0_EB, -1.0_EB, NL)                    !  D := F - D 

ICYCLE = SCARC_CYCLE_CONTROL(NSCARC_CYCLE_SETUP, NL)
RESIN  = SCARC_L2NORM (VEC_D, NL)                                                !  RESIN := ||D||

CALL SCARC_CONVERGENCE_INFO(RESIN, 0, NL, CROUTINE)

!!!----------------------------------------------------------------------------------------------------
!!! Perform multigrid-looping (start each iteration on finest level) 
!!!----------------------------------------------------------------------------------------------------
MULTIGRID_LOOP: DO ITE = 1, NIT
 
   NL = NLEVEL_MIN
   ICYCLE = SCARC_CYCLE_CONTROL(NSCARC_CYCLE_RESET, NL)

   CYCLE_LOOP: DO WHILE (ICYCLE /= NSCARC_CYCLE_EXIT)

      !!! presmoothing  (smoothing/restriction till coarsest level is reached)
      PRESMOOTH_LOOP: DO WHILE (NL < NLEVEL_MAX)

         CALL SCARC_SMOOTHER (VEC_X, VEC_F, VEC_D, NSCARC_CYCLE_PRESMOOTH, NL)   ! D_fine   := smooth(defect)
         CALL SCARC_RESTRICTION (VEC_D, VEC_F, NL, NL+1)                         ! F_coarse := rest(D_fine)
         NL = NL + 1
         CALL SCARC_CLEAR(VEC_X, NL)

      ENDDO PRESMOOTH_LOOP


      !!! coarse grid solver
      SELECT CASE (TYPE_COARSE)
         CASE (NSCARC_COARSE_CG)
            CALL SCARC_METHOD_CG (NSCARC_SCOPE_COARSE, NSCARC_VECTOR_F)          ! X_coarse := exact_sol(.)
         CASE (NSCARC_COARSE_GE)
            CALL SCARC_METHOD_GE (NL)
      END SELECT
      TYPE_SCOPE = NSCOPE

      !!! postsmoothing (smoothing/restriction till finest level is reached again)
      POSTSMOOTH_LOOP: DO WHILE (NL > NLEVEL_MIN) 

         NL=NL-1
         CALL SCARC_PROLONGATION (VEC_X, VEC_D, NL+1, NL)                         ! D_fine := prol(X_coarse)
         CALL SCARC_LINEAR_COMBI (VEC_D, VEC_X, 1.0_EB, 1.0_EB, NL)               ! X_fine := D_fine + X_fine
        
         CALL SCARC_SMOOTHER (VEC_X, VEC_F, VEC_D, NSCARC_CYCLE_POSTSMOOTH, NL)   ! D_fine := smooth(defect)

         ICYCLE = SCARC_CYCLE_CONTROL(NSCARC_CYCLE_PROCEED, NL)
         IF (ICYCLE /= NSCARC_CYCLE_POSTSMOOTH) CYCLE CYCLE_LOOP

      ENDDO POSTSMOOTH_LOOP

   ENDDO CYCLE_LOOP

   IF (NL /= NLEVEL_MIN) THEN
      WRITE(*,*) 'ERROR in SCARC_MULTIGRID, wrong level ', NL
      STOP
   ENDIF

   !!!-------------------------------------------------------------------------------------------------
   !!! Compute norm of new residual on finest level and  leave loop correspondingly
   !!!-------------------------------------------------------------------------------------------------
   CALL SCARC_MATVEC_PRODUCT (VEC_X, VEC_D, NL)                                   ! D := A*X
   CALL SCARC_LINEAR_COMBI   (VEC_F, VEC_D, 1.0_EB, -1.0_EB, NL)                  ! D := F - D

   RES = SCARC_L2NORM (VEC_D, NL)                                                 ! RES := ||D||

   ISTATE = SCARC_CONVERGENCE_STATE(RESIN, RES, EPS, ITE, NL, CROUTINE)
   IF (ISTATE /= NSCARC_STATE_PROCEED) EXIT MULTIGRID_LOOP
 
ENDDO MULTIGRID_LOOP

!!!----------------------------------------------------------------------------------------------------
!!! Determine convergence rate and print corresponding information
!!! In case of MG as main solver:
!!!   - Transfer ScaRC solution vector X to FDS pressure vector 
!!!   - Set ghost cell values along external boundaries
!!!   - Exchange values along internal boundaries 
!!!----------------------------------------------------------------------------------------------------
CALL SCARC_CONVERGENCE_RATE(RESIN, RES, ITE, ISTATE, CROUTINE)

IF (TYPE_SCOPE == NSCARC_SCOPE_MAIN) THEN
   CALL SCARC_TERMINATE_SOLVER(NLEVEL_MIN)
   CALL SCARC_UPDATE_GHOSTCELLS(NLEVEL_MIN)
ENDIF
 
TUSED_SCARC(NSCARC_TIME_MULTIGRID,:)=TUSED_SCARC(NSCARC_TIME_MULTIGRID,:)+SECOND()-TNOW_MULTIGRID
TUSED_SCARC(NSCARC_TIME_TOTAL    ,:)=TUSED_SCARC(NSCARC_TIME_TOTAL    ,:)+SECOND()-TNOW_MULTIGRID
END SUBROUTINE SCARC_METHOD_MULTIGRID


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Control multigrid cycling (F/V/W)
!!! Note: NLEVEL_MIN corresponds to finest level, NLEVEL_MAX to coarsest level
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER FUNCTION SCARC_CYCLE_CONTROL(NSCOPE, NL)
INTEGER, INTENT(IN) :: NSCOPE, NL
INTEGER :: NM, NL0, ICYCLE

SELECT CASE (NSCOPE)

   !!!-------------------------------------------------------------------------------------------------
   !!! initialize cycle counts at beginning of multigrid method
   !!!-------------------------------------------------------------------------------------------------
   CASE (NSCARC_CYCLE_SETUP)

      DO NM = NMESHES_MIN, NMESHES_MAX
         SCARC(NM)%CYCLE_COUNT(2,NL)=1
         DO NL0 = NLEVEL_MIN+1, NLEVEL_MAX - 1
            IF (TYPE_CYCLE==NSCARC_CYCLE_F) THEN
               SCARC(NM)%CYCLE_COUNT(2,NL0)=2
            ELSE
               SCARC(NM)%CYCLE_COUNT(2,NL0)=TYPE_CYCLE
            ENDIF
         ENDDO
      ENDDO

      ICYCLE = NSCARC_CYCLE_NEXT
      
   !!!-------------------------------------------------------------------------------------------------
   !!! reset cycle counts at beginning of each new multigrid iteration
   !!!-------------------------------------------------------------------------------------------------
   CASE (NSCARC_CYCLE_RESET)

      DO NM = NMESHES_MIN, NMESHES_MAX
         DO NL0 = NLEVEL_MIN, NLEVEL_MAX
            SCARC(NM)%CYCLE_COUNT(1,NL0)=SCARC(NM)%CYCLE_COUNT(2,NL0)
         ENDDO
      ENDDO
      ICYCLE = NSCARC_CYCLE_NEXT

   !!!-------------------------------------------------------------------------------------------------
   !!! determine where to proceed with cycling
   !!!-------------------------------------------------------------------------------------------------
   CASE (NSCARC_CYCLE_PROCEED)

      DO NM = NMESHES_MIN, NMESHES_MAX

         SCARC(NM)%CYCLE_COUNT(1,NL)=SCARC(NM)%CYCLE_COUNT(1,NL)-1

         IF (SCARC(NM)%CYCLE_COUNT(1,NL)==0) THEN
            IF (TYPE_CYCLE==NSCARC_CYCLE_F) THEN
               SCARC(NM)%CYCLE_COUNT(1,NL)=1
            ELSE
               SCARC(NM)%CYCLE_COUNT(1,NL)=SCARC(NM)%CYCLE_COUNT(2,NL)
            ENDIF
            IF (NL == NLEVEL_MIN) THEN
               ICYCLE = NSCARC_CYCLE_EXIT
            ELSE
               ICYCLE = NSCARC_CYCLE_POSTSMOOTH
            ENDIF
         ELSE
            IF (NL == NLEVEL_MIN) THEN
               ICYCLE = NSCARC_CYCLE_EXIT
            ELSE
               ICYCLE = NSCARC_CYCLE_NEXT
            ENDIF
         ENDIF

      ENDDO

END SELECT

SCARC_CYCLE_CONTROL = ICYCLE
RETURN

END FUNCTION SCARC_CYCLE_CONTROL


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform smoothing - banded storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SMOOTHER(VEC0_X, VEC0_F, VEC0_D, NTYPE, NL)
INTEGER , INTENT(IN) :: VEC0_X, VEC0_F, VEC0_D, NTYPE, NL
INTEGER :: ITE, NIT, ISTATE
REAL(EB):: RES, RESIN, EPS, OMEGA
REAL(EB):: TNOW_SMOOTH
LOGICAL :: BMATVEC, BL2NORM
CHARACTER(30) :: CROUTINE

TNOW_SMOOTH = SECOND()

!!!----------------------------------------------------------------------------------------------------
!!! Initialization
!!!----------------------------------------------------------------------------------------------------
SELECT CASE (NTYPE)
   CASE (NSCARC_CYCLE_PRESMOOTH)
      CROUTINE = 'SCARC_PRESMOOTHER'
   CASE (NSCARC_CYCLE_POSTSMOOTH)
      CROUTINE = 'SCARC_POSTSMOOTHER'
END SELECT
BL2NORM  = .TRUE.
BMATVEC  = .TRUE.

NIT      = SCARC_SMOOTH_ITERATIONS
EPS      = SCARC_SMOOTH_ACCURACY
OMEGA    = SCARC_SMOOTH_OMEGA
RESIN    = 1.0_EB
ITE      = 0

IF (NTYPE == NSCARC_CYCLE_PRESMOOTH .AND. NL == NLEVEL_MIN) BMATVEC = .FALSE.


!!!----------------------------------------------------------------------------------------------------
!!! Calculate initial defect on level NL (only if BMATVEC = .TRUE.)
!!! Because initial vector is set to zero, this defect corresponds to F
!!!----------------------------------------------------------------------------------------------------
IF (BMATVEC) THEN

   CALL SCARC_MATVEC_PRODUCT (VEC0_X, VEC0_D, NL)                             !  D := A*X
   CALL SCARC_LINEAR_COMBI   (VEC0_F, VEC0_D, 1.0_EB, -1.0_EB, NL)            !  D := F - D

   IF (BL2NORM) THEN
      RESIN = SCARC_L2NORM (VEC0_D, NL)                                       !  RESIN := ||D||
      CALL SCARC_CONVERGENCE_INFO(RESIN, ITE, NL, CROUTINE)
   ENDIF

ENDIF


!!!----------------------------------------------------------------------------------------------------
!!! Smoothing loop
!!!----------------------------------------------------------------------------------------------------
SMOOTH_LOOP: DO ITE=1, NIT
 
   CALL SCARC_PRECONDITIONER (VEC0_D, VEC0_D, NL)                             !  D := PRECON (D)
   CALL SCARC_LINEAR_COMBI   (VEC0_D, VEC0_X, OMEGA, 1.0_EB, NL)              !  X := OMEGA*D + X
   CALL SCARC_MATVEC_PRODUCT (VEC0_X, VEC0_D, NL)                             !  D := A*X
   CALL SCARC_LINEAR_COMBI   (VEC0_F, VEC0_D, 1.0_EB, -1.0_EB, NL)            !  D := F - D

   IF (BL2NORM) THEN
      RES    = SCARC_L2NORM (VEC0_D, NL)                                      !  RES := ||D||
      ISTATE = SCARC_CONVERGENCE_STATE(RESIN, RES, EPS, ITE, NL, CROUTINE)
      IF (ISTATE /= NSCARC_STATE_PROCEED) EXIT SMOOTH_LOOP                    !  RES < TOL ?
   ENDIF

ENDDO SMOOTH_LOOP

!!!----------------------------------------------------------------------------------------------------
!!! Determine convergence rate and print corresponding information
!!!----------------------------------------------------------------------------------------------------
IF (BL2NORM) CALL SCARC_CONVERGENCE_RATE(RESIN, RES, ITE, ISTATE, CROUTINE)

TUSED_SCARC(NSCARC_TIME_SMOOTH,:)=TUSED_SCARC(NSCARC_TIME_SMOOTH,:)+SECOND()-TNOW_SMOOTH
TUSED_SCARC(NSCARC_TIME_TOTAL ,:)=TUSED_SCARC(NSCARC_TIME_TOTAL ,:)+SECOND()-TNOW_SMOOTH
END SUBROUTINE SCARC_SMOOTHER


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Gaussian elimination for coarse grid solution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_METHOD_GE(NL)
INTEGER, INTENT(IN) :: NL

WRITE(*,*) 'SCARC_METHOD_GE not yet implemented, stopping program', NL
STOP

END SUBROUTINE SCARC_METHOD_GE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Initialize GSTRIX preconditioner/smoother
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_GSTRIX(NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IC, IERR
TYPE (SCARC_PRECON_TYPE), POINTER :: SP
TYPE (SCARC_SYSTEM_BANDED_TYPE), POINTER :: SB

IERR=0

SB => SCARC(NM)%BANDED(NL)
SP => SCARC(NM)%PRECON(NL)

SELECT CASE (TYPE_DIMENSION)

   !!!----------------------------------------------------------------------------------------------------
   !!! 2D-version
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_DIMENSION_TWO)

      ALLOCATE (SP%MDX(1:SB%NC), STAT=IERR)
      CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX', 'SP%MDX', IERR)
      
      ALLOCATE (SP%UDX(1:SB%NC), STAT=IERR)
      CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX', 'SP%UDX', IERR)
      
      ALLOCATE (SP%LDX(1:SB%NC), STAT=IERR)
      CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX', 'SP%LDX', IERR)
      
      ALLOCATE (SP%MDZ(1:SB%NC), STAT=IERR)
      CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX', 'SP%MDZ', IERR)
      
      ALLOCATE (SP%DWORK(1:SB%NC), STAT=IERR)
      CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX', 'SP%DWORK', IERR)
      
      !!! save different diagonals of matrix ('D' corresponds to main, 'L' to lower, 'U' to upper)
      IF (TYPE_SYSTEM == NSCARC_SYSTEM_BANDED) THEN
         DO IC = 1, SB%NC
            SP%MDX(IC) = SB%A(IC,ID)         ! main  diagonal in x-direction
            SP%MDZ(IC) = SB%A(IC,ILZ)        ! main  diagonal in z-direction
            SP%LDX(IC) = SB%A(IC,ILX)        ! lower diagonal in x-direction
            SP%UDX(IC) = SB%A(IC,IUX)        ! upper diagonal in x-direction
         ENDDO
      ELSE
         WRITE(*,*) 'HIER NOCH ANPASSEN !!'
         WRITE(SCARC_LU,*) 'HIER NOCH ANPASSEN !!'
         DO IC = 1, SB%NC
   !          SP%MDX(IC) = SC%A(IC,ID)         ! main  diagonal in x-direction
   !         SP%MDZ(IC) = SC%A(IC,ILZ)        ! main  diagonal in z-direction
   !         SP%LDX(IC) = SC%A(IC,ILX)        ! lower diagonal in x-direction
   !         SP%UDX(IC) = SC%A(IC,IUX)        ! upper diagonal in x-direction
         ENDDO
      ENDIF
       
      !!! perform LU-factorization of matrix AG according to banded storage technique
      DO IC = 2, SB%NC
         SP%LDX (IC) = SP%LDX(IC) / SP%MDX(IC-1)
         SP%MDX (IC) = SP%MDX(IC) - SP%LDX(IC) * SP%UDX(IC-1)
      ENDDO
       
      !!! replace diagonal values diag(i) by 1/diag(i) and multiply UDX with 1/diag(i)
      DO IC = 1, SB%NC
         SP%MDX (IC) = 1.0_EB / SP%MDX(IC)
         SP%UDX (IC) = SP%MDX(IC) * SP%UDX(IC)
      ENDDO
      
   !!!----------------------------------------------------------------------------------------------------
   !!! 3D-version
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_DIMENSION_THREE)
   
      ALLOCATE (SP%MDX(1:SB%NC), STAT=IERR)
      CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX3D', 'SP%MDX', IERR)
      
      ALLOCATE (SP%UDX(1:SB%NC), STAT=IERR)
      CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX3D', 'SP%UDX', IERR)
      
      ALLOCATE (SP%LDX(1:SB%NC), STAT=IERR)
      CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX3D', 'SP%LDX', IERR)
      
      ALLOCATE (SP%MDY(1:SB%NC), STAT=IERR)
      CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX3D', 'SP%MDY', IERR)
      
      ALLOCATE (SP%MDZ(1:SB%NC), STAT=IERR)
      CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX3D', 'SP%MDZ', IERR)
      
      ALLOCATE (SP%DWORK(1:SB%NC), STAT=IERR)
      CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX3D', 'SP%DWORK', IERR)
      
      !!! save different diagonals of matrix ('D' corresponds to main, 'L' to lower, 'U' to upper)
      IF (TYPE_SYSTEM == NSCARC_SYSTEM_BANDED) THEN
         DO IC = 1, SB%NC
            SP%MDX(IC) = SB%A(IC,ID)        ! main  diagonal in x-direction
            SP%MDZ(IC) = SB%A(IC,ILZ)       ! main  diagonal in z-direction
            SP%MDY(IC) = SB%A(IC,ILY)       ! main  diagonal in y-direction
            SP%LDX(IC) = SB%A(IC,ILX)       ! lower diagonal in x-direction
            SP%UDX(IC) = SB%A(IC,IUX)       ! upper diagonal in x-direction
         ENDDO
      ELSE
         WRITE(*,*) 'HIER NOCH ANPASSEN !!'
         WRITE(SCARC_LU,*) 'HIER NOCH ANPASSEN !!'
         DO IC = 1, SB%NC
   !         SP%MDX(IC) = SC%A(IC,ID)        ! main  diagonal in x-direction
   !         SP%MDZ(IC) = SC%A(IC,ILZ)       ! main  diagonal in z-direction
   !         SP%MDY(IC) = SC%A(IC,ILY)       ! main  diagonal in y-direction
   !         SP%LDX(IC) = SC%A(IC,ILX)       ! lower diagonal in x-direction
   !         SP%UDX(IC) = SC%A(IC,IUX)       ! upper diagonal in x-direction
         ENDDO
      ENDIF
       
      !!! perform LU-factorization of matrix AG according to banded storage technique
      DO IC = 2, SB%NC
         SP%LDX (IC) = SP%LDX(IC) / SP%MDX(IC-1)
         SP%MDX (IC) = SP%MDX(IC) - SP%LDX(IC) * SP%UDX(IC-1)
      ENDDO
       
      !!! replace diagonal values diag(i) by 1/diag(i) and multiply UDX with 1/diag(i)
      DO IC = 1, SB%NC
         SP%MDX (IC) = 1.0_EB / SP%MDX(IC)
         SP%UDX (IC) = SP%MDX(IC) * SP%UDX(IC)
      ENDDO
   
END SELECT

END SUBROUTINE SCARC_SETUP_GSTRIX
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! GSTRIX preconditioner/smoother
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_GSTRIX (NVECTOR, NL)
INTEGER, INTENT(IN) :: NVECTOR, NL
REAL(EB), POINTER, DIMENSION(:,:,:) :: VB
REAL(EB), POINTER, DIMENSION(:)     :: VC
INTEGER , POINTER :: NX, NY, NZ
INTEGER :: NM, I, J, K, IC, JC
TYPE (SCARC_PRECON_TYPE) , POINTER :: SP
TYPE (SCARC_SYSTEM_BANDED_TYPE) , POINTER :: SB
TYPE (SCARC_SYSTEM_COMPACT_TYPE), POINTER :: SC


SELECT_SYSTEM: SELECT CASE (TYPE_SYSTEM)
   
   !!!-------------------------------------------------------------------------------------------------
   !!! Banded system
   !!!-------------------------------------------------------------------------------------------------
   CASE (NSCARC_SYSTEM_BANDED)

      DO NM = NMESHES_MIN, NMESHES_MAX
      
         VB => POINT_TO_BANDED_VECTOR (NVECTOR, NM, NL)

         SB => SCARC(NM)%BANDED(NL)
         SP => SCARC(NM)%PRECON(NL)
         
         NX => SB%NX
         NY => SB%NY
         NZ => SB%NZ
         WRITE(SCARC_LU,*) 'ACHTUNG, NOCH ANPASSEN!'
         WRITE(*,*) 'ACHTUNG, NOCH ANPASSEN!'
         
         SELECT_BANDED_DIMENSION: SELECT CASE (TYPE_DIMENSION)
         
            !!!--------------- 2D -----------------
            CASE (NSCARC_DIMENSION_TWO)
         
               ! backward elimination of first NX unkowns (may be solved by tridiagonal system)
               DO I=2,NX
                  VB(I,1,1) = VB(I,1,1)-SP%LDX(I-1)*VB(I-1,1,1)
               ENDDO
               DO I=1,NX
                  VB(I,1,1) = VB(I,1,1)*SP%MDX(I)
               ENDDO
               DO I=NX-1,1,-1
                  VB(I,1,1) = VB(I,1,1)-SP%UDX(I)*VB(I+1,1,1)
               ENDDO
                 
               
               ! backward elimination of following unknowns (here the subdiagonal in z-direction must be taken into account)
               DO K=2,NZ
               
                  IC = (K-1)*NX + I
               
                  ! eliminate subdiagonal in z-direction into right hand side (related solution component already known)
                  DO I=1,NX
                     VB(I,1,K) = VB(I,1,K) - SP%MDZ(IC)*VB(I,1,K-1) 
                  ENDDO
               
                  ! perform elimination of matrix lines corresponding to K
                  DO I=2,NX
                     VB(I,1,K) = VB(I,1,K)-SP%LDX(IC-1)*VB(I-1,1,K)
                  ENDDO
                  DO I=1,NX
                     VB(I,1,K) = VB(I,1,K)*SP%MDX(IC)
                  ENDDO
                  DO I=NX-1,1,-1
                     VB(I,1,K) = VB(I,1,K)-SP%UDX(IC)*VB(I+1,1,K)
                  ENDDO
                 
               ENDDO
               
            !!!--------------- 3D -----------------
            CASE (NSCARC_DIMENSION_THREE)
            
               !!! NOT WORKING YET, has still to be adopted to 3D !!!!!
               DO J=1,NY
               
                  ! backward elimination of first NX unkowns (may be solved by tridiagonal system)
                  DO I=2,NX
                     VB(I,J,1) = VB(I,J,1)-SP%LDX(I-1)*VB(I-1,J,1)
                  ENDDO
                  DO I=1,NX
                     VB(I,J,1) = VB(I,J,1)*SP%MDX(I)
                  ENDDO
                  DO I=NX-1,1,-1
                     VB(I,J,1) = VB(I,J,1)-SP%UDX(I)*VB(I+1,J,1)
                  ENDDO
                    
                  ! backward elimination of following unknowns 
                  ! (here the subdiagonals in y- and z-direction must be taken into account)
                  DO K=2,NZ
                  
                     IC = (K-1)*NX + I
                  
                     ! bring subdiagonal in z-direction to right hand side (corresponding solution component already known)
                     DO I=1,NX
                        VB(I,J,K) = VB(I,J,K) - SP%MDZ(IC)*VB(I,J,K-1) 
                     ENDDO
                  
                     ! perform elimination of matrix lines corresponding to K
                     DO I=2,NX
                        VB(I,J,K) = VB(I,J,K)-SP%LDX(IC-1)*VB(I-1,J,K)
                     ENDDO
                     DO I=1,NX
                        VB(I,J,K) = VB(I,J,K)*SP%MDX(IC)
                     ENDDO
                     DO I=NX-1,1,-1
                        VB(I,J,K) = VB(I,J,K)-SP%UDX(IC)*VB(I+1,J,K)
                     ENDDO
                    
                  ENDDO
               
               ENDDO
               
         END SELECT SELECT_BANDED_DIMENSION

      ENDDO
       

   !!!-------------------------------------------------------------------------------------------------
   !!! Compact system
   !!!-------------------------------------------------------------------------------------------------
   CASE (NSCARC_SYSTEM_COMPACT)

      DO NM = NMESHES_MIN, NMESHES_MAX

         VC => POINT_TO_COMPACT_VECTOR (NVECTOR, NM, NL)

         SC => SCARC(NM)%COMPACT(NL)
         SP => SCARC(NM)%PRECON(NL)
         
         NX => SC%NX
         NY => SC%NY
         NZ => SC%NZ
         
         SELECT_COMPACT_DIMENSION: SELECT CASE (TYPE_DIMENSION)
         
            !!!--------------- 2D -----------------
            CASE (NSCARC_DIMENSION_TWO)
         
               ! backward elimination of first NX unkowns (may be solved by tridiagonal system)
               DO IC=2,NX
                  VC(IC) = VC(IC)-SP%LDX(I-1)*VC(IC-1)
               ENDDO
               DO I=1,NX
                  VC(IC) = VC(IC)*SP%MDX(I)
               ENDDO
               DO I=NX-1,1,-1
                  VC(IC) = VC(IC)-SP%UDX(I)*VC(IC+1)
               ENDDO
                 
               
               ! backward elimination of following unknowns (here the subdiagonal in z-direction must be taken into account)
               DO K=2,NZ
               
                  IC = (K-1)*NX + I
               
                  ! eliminate subdiagonal in z-direction into right hand side (related solution component already known)
                  DO I=1,NX
                     VC(IC) = VC(IC) - SP%MDZ(IC)*VC(IC-1) 
                  ENDDO
               
                  ! perform elimination of matrix lines corresponding to K
                  DO I=2,NX
                     VC(IC) = VC(IC)-SP%LDX(IC-1)*VC(IC-1)
                  ENDDO
                  DO I=1,NX
                     VC(IC) = VC(IC)*SP%MDX(IC)
                  ENDDO
                  DO I=NX-1,1,-1
                     VC(IC) = VC(IC)-SP%UDX(IC)*VC(IC+1)
                  ENDDO
                 
               ENDDO
               
            !!!--------------- 3D -----------------
            CASE (NSCARC_DIMENSION_THREE)
            
               !!! NOT WORKING YET, has still to be adopted to 3D !!!!!
               DO J=1,NY
               
                  JC = J*NX
            
                  ! backward elimination of first NX unkowns (may be solved by tridiagonal system)
                  DO I=2,NX
                     VC(JC+I) = VC(JC+I)-SP%LDX(I-1)*VC(JC+I-1)
                  ENDDO
                  DO I=1,NX
                     VC(JC+I) = VC(JC+I)*SP%MDX(I)
                  ENDDO
                  DO I=NX-1,1,-1
                     VC(JC+I) = VC(JC+I)-SP%UDX(I)*VC(JC+I+1)
                  ENDDO
                    
                  ! backward elimination of following unknowns 
                  ! (here the subdiagonals in y- and z-direction must be taken into account)
                  DO K=2,NZ
                  
                     IC = (K-1)*NX*NY + (J-1)*NX + I
                  
                     ! bring subdiagonal in z-direction to right hand side (corresponding solution component already known)
                     DO I=1,NX
                        VC(IC) = VC(IC) - SP%MDZ(IC)*VC(IC-1) 
                     ENDDO
                  
                     ! perform elimination of matrix lines corresponding to K
                     DO I=2,NX
                        VC(IC) = VC(IC)-SP%LDX(IC-1)*VC(IC-1)
                     ENDDO
                     DO I=1,NX
                        VC(IC) = VC(IC)*SP%MDX(IC)
                     ENDDO
                     DO I=NX-1,1,-1
                        VC(IC) = VC(IC)-SP%UDX(IC)*VC(IC+1)
                     ENDDO
                    
                  ENDDO
               
               ENDDO
            
         END SELECT SELECT_COMPACT_DIMENSION
          
      ENDDO

END SELECT SELECT_SYSTEM

END SUBROUTINE SCARC_GSTRIX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Set initial solution corresponding to boundary data in BXS, BXF, ...
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_SOLVER(NL)
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IW, IOR0, I, J, K, IC
TYPE (      MESH_TYPE), POINTER ::  M
TYPE (SCARC_MESH_TYPE), POINTER :: SM
TYPE (SCARC_SYSTEM_BANDED_TYPE) , POINTER :: SB
TYPE (SCARC_SYSTEM_COMPACT_TYPE), POINTER :: SC

SELECT_SYSTEM: SELECT CASE (TYPE_SYSTEM)

   !!!-------------------------------------------------------------------------------------------------------
   !!! Banded system
   !!!-------------------------------------------------------------------------------------------------------
   CASE (NSCARC_SYSTEM_BANDED)

      SELECT_BANDED_METHOD: SELECT CASE (TYPE_METHOD)

         !!! In case of a Krylov method initialize auxiliary arrays
         CASE (NSCARC_METHOD_KRYLOV)

            DO NM = NMESHES_MIN, NMESHES_MAX
               SCARC(NM)%BANDED(NL)%D = 0.0_EB
               SCARC(NM)%BANDED(NL)%G = 0.0_EB
               SCARC(NM)%BANDED(NL)%Y = 0.0_EB
               SCARC(NM)%BANDED(NL)%W = 0.0_EB
            ENDDO

            IF (TYPE_KRYLOV == NSCARC_KRYLOV_BICG) THEN
               DO NM = NMESHES_MIN, NMESHES_MAX
                  SCARC(NM)%BANDED(NL)%Z = 0.0_EB
               ENDDO
            ENDIF
            
         !!! In case of a multigrid method with coarse grid solution by CG, clear CG-vectors on max level
         CASE (NSCARC_METHOD_MULTIGRID)
            
            DO NM = NMESHES_MIN, NMESHES_MAX
               SCARC(NM)%BANDED(NLEVEL_MAX)%X = 0.0_EB
               SCARC(NM)%BANDED(NLEVEL_MAX)%D = 0.0_EB
               SCARC(NM)%BANDED(NLEVEL_MAX)%G = 0.0_EB
               SCARC(NM)%BANDED(NLEVEL_MAX)%Y = 0.0_EB
               SCARC(NM)%BANDED(NLEVEL_MAX)%W = 0.0_EB
            ENDDO

      END SELECT SELECT_BANDED_METHOD

      IF (TYPE_SCOPE == NSCARC_SCOPE_MAIN) THEN
               
         !!! Initialize solution and right hand side vector corresponding to boundary conditions
         SELECT_BANDED_DIMENSION: SELECT CASE (TYPE_DIMENSION)
         
            !!! -------------- 2D ------------------
            CASE (NSCARC_DIMENSION_TWO)
            
               DO NM = NMESHES_MIN, NMESHES_MAX
         
                  M  => MESHES(NM)
                  SM => SCARC(NM)%MESHES(NL)
                  SB => SCARC(NM)%BANDED(NL)
               
                  ! get right hand side vector
                  SB%F (1:M%IBAR, 1, 1:M%KBAR) = M%PRHS (1:M%IBAR, 1, 1:M%KBAR)
            
                  ! set ghost values of solution vector to zero and use values of H/HS in interior
                  SB%X(1:M%IBAR, 1, 0       ) = 0.0_EB
                  SB%X(1:M%IBAR, 1, M%KBAR+1) = 0.0_EB
                  SB%X(0       , 1, 1:M%KBAR) = 0.0_EB
                  SB%X(M%IBAR+1, 1, 1:M%KBAR) = 0.0_EB
                  IF (PREDICTOR) THEN
                     SB%X (1:M%IBAR, 1, 1:M%KBAR) = M%H (1:M%IBAR, 1, 1:M%KBAR)
                  ELSE
                     SB%X (1:M%IBAR, 1, 1:M%KBAR) = M%HS(1:M%IBAR, 1, 1:M%KBAR)
                  ENDIF
            
                  ! set boundary conditions at exterior boundaries (corresponding to pois.f90)
                  BANDED_WALLCELL_LOOP2D: DO IW = 1, M%N_EXTERNAL_WALL_CELLS
               
                     I = M%IJKW(6,IW)
                     J = M%IJKW(7,IW)
                     K = M%IJKW(8,IW)
            
                     IF (J /= 1) THEN
                        WRITE(*,*) 'Wrong index for J =',J,' in SCARC_SETUP_SOLVER !!!'
                        STOP
                     ENDIF
                  
                     IOR0 = M%IJKW(4,IW)
                  
                     SELECT CASE (IOR0)
                        CASE (1)
                           IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
                              SB%F(I,J,K) = SB%F(I,J,K) - 2.0_EB * SM%DXI2 * M%BXS(1,K)         ! Dirichlet
                           ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
                              SB%F(I,J,K) = SB%F(I,J,K) + SM%DXI * M%BXS(1,K)                   ! Neumann
                           ENDIF
                        CASE (-1)
                           IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
                              SB%F(I,J,K) = SB%F(I,J,K) - 2.0_EB * SM%DXI2 *M%BXF(1,K)          ! Dirichlet
                           ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
                              SB%F(I,J,K) = SB%F(I,J,K) - SM%DXI *M%BXF(1,K)                    ! Neumann
                           ENDIF
                        CASE (3)
                           IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
                              SB%F(I,J,K) = SB%F(I,J,K) - 2.0_EB * SM%DZI2 * M%BZS(I,1)         ! Dirichlet
                           ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
                              SB%F(I,J,K) = SB%F(I,J,K) + SM%DZI * M%BZS(I,1)                   ! Neumann
                           ENDIF
                        CASE (-3)
                           IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
                              SB%F(I,J,K) = SB%F(I,J,K) - 2.0_EB * SM%DZI2 * M%BZF(I,1)         ! Dirichlet
                           ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
                              SB%F(I,J,K) = SB%F(I,J,K) - SM%DZI  * M%BZF(I,1)                  ! Neumann
                           ENDIF
                     END SELECT
                  
                  ENDDO BANDED_WALLCELL_LOOP2D
               
               ENDDO
            
            !!! -------------- 3D ------------------
            CASE (NSCARC_DIMENSION_THREE)
            
               DO NM = NMESHES_MIN, NMESHES_MAX
               
                  M  => MESHES(NM)
                  SM => SCARC(NM)%MESHES(NL)
                  SB => SCARC(NM)%BANDED(NL)
            
                  ! get right hand side vector
                  SB%F (1:M%IBAR, 1:M%JBAR, 1:M%KBAR) = M%PRHS (1:M%IBAR, 1:M%JBAR, 1:M%KBAR)
            
                  ! set ghost values of solution vector to zero and use values of H/HS in interior
                  SB%X(1:M%IBAR, 1:M%JBAR, 0       ) = 0.0_EB
                  SB%X(1:M%IBAR, 1:M%JBAR, M%KBAR+1) = 0.0_EB
                  SB%X(1:M%IBAR, 0       , 1:M%KBAR) = 0.0_EB
                  SB%X(1:M%IBAR, M%JBAR+1, 1:M%KBAR) = 0.0_EB
                  SB%X(0       , 1:M%JBAR, 1:M%KBAR) = 0.0_EB
                  SB%X(M%IBAR+1, 1:M%JBAR, 1:M%KBAR) = 0.0_EB
                  IF (PREDICTOR) THEN
                     SB%X (1:M%IBAR, 1:M%JBAR, 1:M%KBAR) = M%H (1:M%IBAR, 1:M%JBAR, 1:M%KBAR)
                  ELSE
                     SB%X (1:M%IBAR, 1:M%JBAR, 1:M%KBAR) = M%HS(1:M%IBAR, 1:M%JBAR, 1:M%KBAR)
                  ENDIF
            
                  BANDED_WALLCELL_LOOP3D: DO IW = 1, M%N_EXTERNAL_WALL_CELLS
                  
                     I = M%IJKW(6,IW)
                     J = M%IJKW(7,IW)
                     K = M%IJKW(8,IW)
                  
                     IOR0 = M%IJKW(4,IW)
                  
                     SELECT CASE (IOR0)
                        CASE (1)
                           IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
                              SB%F(I,J,K) = SB%F(I,J,K) - 2.0_EB * SM%DXI2 * M%BXS(J,K)         ! Dirichlet
                           ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
                              SB%F(I,J,K) = SB%F(I,J,K) + SM%DXI * M%BXS(J,K)                   ! Neumann
                           ENDIF
                        CASE (-1)
                           IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
                              SB%F(I,J,K) = SB%F(I,J,K) - 2.0_EB * SM%DXI2 *M%BXF(J,K)          ! Dirichlet
                           ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
                              SB%F(I,J,K) = SB%F(I,J,K) - SM%DXI *M%BXF(J,K)                    ! Neumann
                           ENDIF
                        CASE (2)
                           IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
                              SB%F(I,J,K) = SB%F(I,J,K) - 2.0_EB * SM%DYI2 * M%BYS(I,K)         ! Dirichlet
                           ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
                              SB%F(I,J,K) = SB%F(I,J,K) + SM%DYI * M%BYS(I,K)                   ! Neumann
                           ENDIF
                        CASE (-2)
                           IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
                              SB%F(I,J,K) = SB%F(I,J,K) - 2.0_EB * SM%DYI2 *M%BYF(I,K)          ! Dirichlet
                           ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
                              SB%F(I,J,K) = SB%F(I,J,K) - SM%DYI *M%BYF(I,K)                    ! Neumann
                           ENDIF
                        CASE (3)
                           IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
                              SB%F(I,J,K) = SB%F(I,J,K) - 2.0_EB * SM%DZI2 * M%BZS(I,J)         ! Dirichlet
                           ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
                              SB%F(I,J,K) = SB%F(I,J,K) + SM%DZI * M%BZS(I,J)                   ! Neumann
                           ENDIF
                        CASE (-3)
                           IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
                              SB%F(I,J,K) = SB%F(I,J,K) - 2.0_EB * SM%DZI2 * M%BZF(I,J)         ! Dirichlet
                           ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
                              SB%F(I,J,K) = SB%F(I,J,K) - SM%DZI  * M%BZF(I,J)                  ! Neumann
                           ENDIF
                     END SELECT
                  
                  ENDDO BANDED_WALLCELL_LOOP3D
            
               ENDDO
         
         END SELECT SELECT_BANDED_DIMENSION

      ENDIF
   
         
   !!!-------------------------------------------------------------------------------------------------------
   !!! Compact system
   !!!-------------------------------------------------------------------------------------------------------
   CASE (NSCARC_SYSTEM_COMPACT)

      SELECT_COMPACT_METHOD: SELECT CASE (TYPE_METHOD)

         !!! In case of a Krylov method initialize auxiliary arrays
         CASE (NSCARC_METHOD_KRYLOV)

            DO NM = NMESHES_MIN, NMESHES_MAX
               SCARC(NM)%COMPACT(NL)%D = 0.0_EB
               SCARC(NM)%COMPACT(NL)%G = 0.0_EB
               SCARC(NM)%COMPACT(NL)%Y = 0.0_EB
               SCARC(NM)%COMPACT(NL)%W = 0.0_EB
            ENDDO

            IF (TYPE_KRYLOV == NSCARC_KRYLOV_BICG) THEN
               DO NM = NMESHES_MIN, NMESHES_MAX
                  SCARC(NM)%COMPACT(NL)%Z = 0.0_EB
               ENDDO
            ENDIF
            
         !!! In case of a multigrid method with coarse grid solution by CG, clear CG-vectors on max level
         CASE (NSCARC_METHOD_MULTIGRID)
            
            DO NM = NMESHES_MIN, NMESHES_MAX
               SCARC(NM)%COMPACT(NLEVEL_MAX)%X = 0.0_EB
               SCARC(NM)%COMPACT(NLEVEL_MAX)%D = 0.0_EB
               SCARC(NM)%COMPACT(NLEVEL_MAX)%G = 0.0_EB
               SCARC(NM)%COMPACT(NLEVEL_MAX)%Y = 0.0_EB
               SCARC(NM)%COMPACT(NLEVEL_MAX)%W = 0.0_EB
            ENDDO

      END SELECT SELECT_COMPACT_METHOD
      
      IF (TYPE_SCOPE == NSCARC_SCOPE_MAIN) THEN
         
         !!! Initialize solution and right hand side vector corresponding to boundary conditions
         SELECT_COMPACT_DIMENSION: SELECT CASE (TYPE_DIMENSION)
         
            !!! -------------- 2D ------------------
            CASE (NSCARC_DIMENSION_TWO)
         
               DO NM = NMESHES_MIN, NMESHES_MAX
         
                  M  => MESHES(NM)
                  SM => SCARC(NM)%MESHES(NL)
                  SC => SCARC(NM)%COMPACT(NL)
               
                  ! get right hand side and initial vector
                  DO K = 1, M%KBAR
                     DO I = 1, M%IBAR
                        IC = (K-1)*M%IBAR + I
                        SC%F(IC) = M%PRHS (I, 1, K)
                        IF (PREDICTOR) THEN
                           SC%X(IC) = M%H (I, 1, K)
                        ELSE
                          SC% X(IC) = M%HS(I, 1, K)
                        ENDIF
                     ENDDO
                  ENDDO
            
            
                  ! set boundary conditions at exterior boundaries (corresponding to pois.f90)
                  COMPACT_WALLCELL_LOOP2D: DO IW = 1, M%N_EXTERNAL_WALL_CELLS
               
                     I = M%IJKW(6,IW)
                     J = M%IJKW(7,IW)
                     K = M%IJKW(8,IW)
            
                     IF (J /= 1) THEN
                        WRITE(*,*) 'Wrong index for J =',J,' in SCARC_SETUP_SOLVER !!!'
                        STOP
                     ENDIF
                  
                     IOR0 = M%IJKW(4,IW)
                     IC = (K-1)*M%IBAR + I
                     
                     SELECT CASE (IOR0)
                        CASE (1)
                           IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
                              SC%F(IC) = SC%F(IC) - 2.0_EB * SM%DXI2 * M%BXS(1,K)         ! Dirichlet
                           ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
                              SC%F(IC) = SC%F(IC) + SM%DXI * M%BXS(1,K)                   ! Neumann
                           ENDIF
                        CASE (-1)
                           IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
                              SC%F(IC) = SC%F(IC) - 2.0_EB * SM%DXI2 *M%BXF(1,K)          ! Dirichlet
                           ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
                              SC%F(IC) = SC%F(IC) - SM%DXI *M%BXF(1,K)                    ! Neumann
                           ENDIF
                        CASE (3)
                           IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
                              SC%F(IC) = SC%F(IC) - 2.0_EB * SM%DZI2 * M%BZS(I,1)         ! Dirichlet
                           ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
                              SC%F(IC) = SC%F(IC) + SM%DZI * M%BZS(I,1)                   ! Neumann
                           ENDIF
                        CASE (-3)
                           IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
                              SC%F(IC) = SC%F(IC) - 2.0_EB * SM%DZI2 * M%BZF(I,1)         ! Dirichlet
                           ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
                              SC%F(IC) = SC%F(IC) - SM%DZI  * M%BZF(I,1)                  ! Neumann
                           ENDIF
                     END SELECT
                  
                  ENDDO COMPACT_WALLCELL_LOOP2D
         
               ENDDO
            !!! -------------- 3D ------------------
            CASE (NSCARC_DIMENSION_THREE)
            
               DO NM = NMESHES_MIN, NMESHES_MAX
         
                  M  => MESHES(NM)
                  SM => SCARC(NM)%MESHES(NL)
                  SC => SCARC(NM)%COMPACT(NL)
            
                  ! get right hand side and initial vector
                  DO K = 1, M%KBAR
                     DO J = 1, M%JBAR
                        DO I = 1, M%IBAR
                           IC = (K-1) * M%IBAR * M%JBAR + (J-1) * M%IBAR + I
                           SC%F(IC) = M%PRHS (I, J, K)
                           IF (PREDICTOR) THEN
                              SC%X(IC) = M%H (I, J, K)
                           ELSE
                              SC%X(IC) = M%HS(I, J, K)
                           ENDIF
                        ENDDO
                     ENDDO
                  ENDDO
            
                  COMPACT_WALLCELL_LOOP3D: DO IW = 1, M%N_EXTERNAL_WALL_CELLS
                  
                     I = M%IJKW(6,IW)
                     J = M%IJKW(7,IW)
                     K = M%IJKW(8,IW)
                  
                     IOR0 = M%IJKW(4,IW)
                     IC = (K-1) * M%IBAR * M%JBAR + (J-1) * M%IBAR + I
                  
                     SELECT CASE (IOR0)
                        CASE (1)
                           IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
                              SC%F(IC) = SC%F(IC) - 2.0_EB * SM%DXI2 * M%BXS(J,K)         ! Dirichlet
                           ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
                              SC%F(IC) = SC%F(IC) + SM%DXI * M%BXS(J,K)                   ! Neumann
                           ENDIF
                        CASE (-1)
                           IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
                              SC%F(IC) = SC%F(IC) - 2.0_EB * SM%DXI2 *M%BXF(J,K)          ! Dirichlet
                           ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
                              SC%F(IC) = SC%F(IC) - SM%DXI *M%BXF(J,K)                    ! Neumann
                           ENDIF
                        CASE (2)
                           IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
                              SC%F(IC) = SC%F(IC) - 2.0_EB * SM%DYI2 * M%BYS(I,K)         ! Dirichlet
                           ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
                              SC%F(IC) = SC%F(IC) + SM%DYI * M%BYS(I,K)                   ! Neumann
                           ENDIF
                        CASE (-2)
                           IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
                              SC%F(IC) = SC%F(IC) - 2.0_EB * SM%DYI2 *M%BYF(I,K)          ! Dirichlet
                           ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
                              SC%F(IC) = SC%F(IC) - SM%DYI *M%BYF(I,K)                    ! Neumann
                           ENDIF
                        CASE (3)
                           IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
                              SC%F(IC) = SC%F(IC) - 2.0_EB * SM%DZI2 * M%BZS(I,J)         ! Dirichlet
                           ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
                              SC%F(IC) = SC%F(IC) + SM%DZI * M%BZS(I,J)                   ! Neumann
                           ENDIF
                        CASE (-3)
                           IF (M%PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN
                              SC%F(IC) = SC%F(IC) - 2.0_EB * SM%DZI2 * M%BZF(I,J)         ! Dirichlet
                           ELSE IF (M%PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
                              SC%F(IC) = SC%F(IC) - SM%DZI  * M%BZF(I,J)                  ! Neumann
                           ENDIF
                     END SELECT
                  
                  ENDDO COMPACT_WALLCELL_LOOP3D
         
               ENDDO
         
         END SELECT SELECT_COMPACT_DIMENSION

      ENDIF

END SELECT SELECT_SYSTEM

END SUBROUTINE SCARC_SETUP_SOLVER


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Print out residual information for loop ITE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_CONVERGENCE_INFO(RES, ITE, NL, CROUTINE)
INTEGER, INTENT(IN) :: ITE, NL
INTEGER  :: NM
REAL(EB), INTENT(IN) :: RES
CHARACTER(*), INTENT(IN) :: CROUTINE

IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS .AND. MYID==0) write(*,1000) TRIM(CROUTINE), NL, ITE,  RES
DO NM = NMESHES_MIN, NMESHES_MAX
   IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS) WRITE(SCARC_LU,1000) TRIM(CROUTINE), NL, ITE,  RES
ENDDO

1000 FORMAT (5X,A30,': level=',i4,': #ite= ',i4,': res =',e14.6)
END SUBROUTINE SCARC_CONVERGENCE_INFO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Check if solver converges or diverges and print out residual information
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER FUNCTION SCARC_CONVERGENCE_STATE(RESIN, RES, EPS, ITE, NL, CROUTINE)

INTEGER, INTENT(IN) :: ITE, NL
INTEGER :: NM, ISTATE
REAL(EB), INTENT(IN) :: RESIN, RES, EPS
CHARACTER(*), INTENT(IN) :: CROUTINE

ISTATE = NSCARC_STATE_PROCEED

DO NM = NMESHES_MIN, NMESHES_MAX
   IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS)          WRITE(SCARC_LU,1000) TRIM(CROUTINE), NL, ITE, RES
ENDDO
IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS.AND.MYID==0) WRITE(*       ,1000) TRIM(CROUTINE), NL, ITE, RES

SELECT CASE (TYPE_ACCURACY)
   CASE (NSCARC_ACCURACY_RELATIVE)
      IF (RES <= RESIN*EPS)  ISTATE = NSCARC_STATE_CONV
   CASE (NSCARC_ACCURACY_ABSOLUTE)
      IF (RES <= EPS .AND. RES <= RESIN*SCARC_ACCURACY_RELATIVE) ISTATE = NSCARC_STATE_CONV
END SELECT
IF (RES > SCARC_ACCURACY_DIVERGENCE) ISTATE = NSCARC_STATE_DIVG

SCARC_CONVERGENCE_STATE = ISTATE
RETURN

1000 FORMAT (5X,A30,': level=',i4,': #ite= ',i4,': res =',e14.6)
END FUNCTION SCARC_CONVERGENCE_STATE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Compute convergence rate and print out resiual information for final loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_CONVERGENCE_RATE(RESIN, RES, ITE, ISTATE, CROUTINE)

INTEGER, INTENT(IN) :: ITE, ISTATE
INTEGER :: NM
REAL(EB), INTENT(IN) :: RESIN, RES
CHARACTER(*), INTENT(IN) :: CROUTINE

SCARC_RESIDUAL = RES
IF (ISTATE == NSCARC_STATE_DIVG) THEN                       
   SCARC_ITERATIONS = - 1
   SCARC_CAPPA      = 1.0_EB
ELSE 
   IF (ISTATE == NSCARC_STATE_CONV) THEN
     SCARC_ITERATIONS = ITE
   ELSE
     SCARC_ITERATIONS = ITE-1
   ENDIF
   IF (RESIN >= 1.0E-70_EB) THEN
      SCARC_CAPPA = (RES/RESIN) ** (1.0_EB/SCARC_ITERATIONS)
   ELSE 
      SCARC_CAPPA = 0.0_EB
   ENDIF
ENDIF

DO NM = NMESHES_MIN, NMESHES_MAX
   IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS)          WRITE(SCARC_LU,2000) TRIM(CROUTINE), SCARC_CAPPA 
ENDDO
IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS.AND.MYID==0) WRITE(*       ,2000) TRIM(CROUTINE), SCARC_CAPPA 

2000 FORMAT (5X,A30,':',10X,'---> convergence rate =',e14.6,/, &
             5X ,'----------------------------------------------------------------------------------------')
END SUBROUTINE SCARC_CONVERGENCE_RATE

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform restriction from finer to coarser grid in multigrid method
!!!    - 'FI' corresponds to finer   grid
!!!    - 'CO' corresponds to coarser grid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_RESTRICTION (NVECTOR_FI, NVECTOR_CO, NL_FI, NL_CO)
INTEGER, INTENT(IN) :: NVECTOR_FI, NVECTOR_CO, NL_FI, NL_CO
INTEGER :: NM, ICOL, IC
REAL(EB), POINTER, DIMENSION(:,:,:) :: FB_CO, DB_FI
REAL(EB), POINTER, DIMENSION(:)     :: FC_CO, DC_FI, R
INTEGER , POINTER, DIMENSION(:)     :: R_ROW, R_COL
INTEGER , POINTER :: NX_CO, NY_CO, NZ_CO, NC_CO
INTEGER  :: NX_FI, NY_FI, NZ_FI
INTEGER  :: IX_FI, IY_FI, IZ_FI, IC_FI(8)
INTEGER  :: IX_CO, IY_CO, IZ_CO, IC_CO
REAL(EB) :: AUX

SELECT_SYSTEM: SELECT CASE (TYPE_SYSTEM)

   !!!----------------------------------------------------------------------------------------------------
   !!! Banded system
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_SYSTEM_BANDED)

      SELECT_BANDED_MULTIGRID: SELECT CASE (TYPE_MULTIGRID)

         !!!------------------------- Geometric multigrid -----------------------------------------------
         CASE (NSCARC_MULTIGRID_GEOMETRIC)

            DO NM = NMESHES_MIN, NMESHES_MAX
            
               NX_CO => SCARC(NM)%BANDED(NL_CO)%NX
               NY_CO => SCARC(NM)%BANDED(NL_CO)%NY
               NZ_CO => SCARC(NM)%BANDED(NL_CO)%NZ
            
               DB_FI => POINT_TO_BANDED_VECTOR(NVECTOR_FI, NM, NL_FI)
               FB_CO => POINT_TO_BANDED_VECTOR(NVECTOR_CO, NM, NL_CO)
      
               SELECT_BANDED_DIMENSION: SELECT CASE (TYPE_DIMENSION)

                  !!!----------------- 2D -------------------
                  CASE (NSCARC_DIMENSION_TWO)
               
                     DO IZ_CO = 1, NZ_CO
                        DO IX_CO = 1, NX_CO
                 
                           IX_FI = 2*IX_CO
                           IY_FI = 1
                           IZ_FI = 2*IZ_CO
                
                           FB_CO(IX_CO, 1, IZ_CO) = 0.25_EB * (  DB_FI(IX_FI  , IY_FI, IZ_FI-1)  &
                                                               + DB_FI(IX_FI-1, IY_FI, IZ_FI-1)  &
                                                               + DB_FI(IX_FI  , IY_FI, IZ_FI  )  &
                                                               + DB_FI(IX_FI-1, IY_FI, IZ_FI  ) )
                        ENDDO
                     ENDDO
               
                  !!!----------------- 3D -------------------
                  CASE (NSCARC_DIMENSION_THREE)
                     
                     DO IZ_CO = 1, NZ_CO
                        DO IY_CO = 1, NY_CO
                           DO IX_CO = 1, NX_CO
               
                              IX_FI = 2*IX_CO
                              IY_FI = 2*IY_CO
                              IZ_FI = 2*IZ_CO
               
                              FB_CO(IX_CO, IY_CO , IZ_CO) = 0.125_EB * (  DB_FI(IX_FI-1, IY_FI-1, IZ_FI-1)  &
                                                                        + DB_FI(IX_FI  , IY_FI-1, IZ_FI-1)  &
                                                                        + DB_FI(IX_FI-1, IY_FI  , IZ_FI-1)  &
                                                                        + DB_FI(IX_FI  , IY_FI  , IZ_FI-1)  &
                                                                        + DB_FI(IX_FI-1, IY_FI-1, IZ_FI  )  &
                                                                        + DB_FI(IX_FI  , IY_FI-1, IZ_FI  )  &
                                                                        + DB_FI(IX_FI-1, IY_FI  , IZ_FI  )  &
                                                                        + DB_FI(IX_FI  , IY_FI  , IZ_FI  ) )
                           ENDDO
                        ENDDO
                     ENDDO
            
               END SELECT SELECT_BANDED_DIMENSION

            ENDDO
                  
      
         !!!------------------------- Geometric multigrid -----------------------------------------------
         CASE (NSCARC_MULTIGRID_ALGEBRAIC)

            WRITE(*,*) 'No banded restriction for algebraic multigrid available, stopping program!'
            STOP

      END SELECT SELECT_BANDED_MULTIGRID

      
   !!!----------------------------------------------------------------------------------------------------
   !!! Compact system
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_SYSTEM_COMPACT)

      SELECT_COMPACT_MULTIGRID: SELECT CASE (TYPE_MULTIGRID)

         !!!------------------------- Geometric multigrid -----------------------------------------------
         CASE (NSCARC_MULTIGRID_GEOMETRIC)

            DO NM = NMESHES_MIN, NMESHES_MAX
         
               NX_CO => SCARC(NM)%COMPACT(NL_CO)%NX
               NY_CO => SCARC(NM)%COMPACT(NL_CO)%NY
               NZ_CO => SCARC(NM)%COMPACT(NL_CO)%NZ
         
               NX_FI = 2*NX_CO
               NY_FI = 2*NY_CO
               NZ_FI = 2*NZ_CO
   
               DC_FI  => POINT_TO_COMPACT_VECTOR(NVECTOR_FI, NM, NL_FI)
               FC_CO  => POINT_TO_COMPACT_VECTOR(NVECTOR_CO, NM, NL_CO)
   
               SELECT_COMPACT_DIMENSION: SELECT CASE (TYPE_DIMENSION)
               
                  !!!----------------- 2D -------------------
                  CASE (NSCARC_DIMENSION_TWO)
               
                     DO IZ_CO = 1, NZ_CO
                        DO IX_CO = 1, NX_CO
                  
                           IX_FI = 2*IX_CO
                           IZ_FI = 2*IZ_CO
                 
                           IC_CO     = (IZ_CO-1)*NX_CO + IX_CO
                
                           IC_FI(1) = (IZ_FI-2)*NX_FI + IX_FI - 1
                           IC_FI(2) = (IZ_FI-2)*NX_FI + IX_FI   
                           IC_FI(3) = (IZ_FI-1)*NX_FI + IX_FI - 1
                           IC_FI(4) = (IZ_FI-1)*NX_FI + IX_FI    
               
                           FC_CO(IC_CO) = 0.25_EB * (  DC_FI(IC_FI(1)) &
                                                     + DC_FI(IC_FI(2)) &
                                                     + DC_FI(IC_FI(3)) &
                                                     + DC_FI(IC_FI(4)) )
               
                        ENDDO
                     ENDDO
                  
                  !!!----------------- 3D -------------------
                  CASE (NSCARC_DIMENSION_THREE)
                     
                     DO IZ_CO = 1, NZ_CO
                        DO IY_CO = 1, NY_CO
                           DO IX_CO = 1, NX_CO
                  
                              IX_FI = 2*IX_CO
                              IY_FI = 2*IY_CO
                              IZ_FI = 2*IZ_CO
                 
                              IC_CO    = (IZ_CO-1)*NX_CO*NY_CO + (IY_CO-1)*NX_CO + IX_CO
                
                              IC_FI(1) = (IZ_FI-1)*NX_FI*NY_FI + (IY_FI-1)*NX_FI + IX_FI - 1
                              IC_FI(2) = (IZ_FI-1)*NX_FI*NY_FI + (IY_FI-1)*NX_FI + IX_FI    
                              IC_FI(3) = (IZ_FI-1)*NX_FI*NY_FI +  IY_FI   *NX_FI + IX_FI - 1
                              IC_FI(4) = (IZ_FI-1)*NX_FI*NY_FI +  IY_FI   *NX_FI + IX_FI    
                              IC_FI(5) =  IZ_FI   *NX_FI*NY_FI + (IY_FI-1)*NX_FI + IX_FI - 1
                              IC_FI(6) =  IZ_FI   *NX_FI*NY_FI + (IY_FI-1)*NX_FI + IX_FI    
                              IC_FI(7) =  IZ_FI   *NX_FI*NY_FI +  IY_FI   *NX_FI + IX_FI - 1
                              IC_FI(8) =  IZ_FI   *NX_FI*NY_FI +  IY_FI   *NX_FI + IX_FI    
               
                              FC_CO(IC_CO) = 0.125_EB * (  DC_FI(IC_FI(1)) &
                                                         + DC_FI(IC_FI(2)) &
                                                         + DC_FI(IC_FI(3)) &
                                                         + DC_FI(IC_FI(4)) &
                                                         + DC_FI(IC_FI(5)) &
                                                         + DC_FI(IC_FI(6)) &
                                                         + DC_FI(IC_FI(7)) &
                                                         + DC_FI(IC_FI(8)) )
                           ENDDO
                        ENDDO
                     ENDDO
               
               END SELECT SELECT_COMPACT_DIMENSION

            ENDDO

         !!!------------------------- Algebraic multigrid -----------------------------------------------
         CASE (NSCARC_MULTIGRID_ALGEBRAIC)

            DO NM = NMESHES_MIN, NMESHES_MAX

               DC_FI  => POINT_TO_COMPACT_VECTOR(NVECTOR_FI, NM, NL_FI)
               FC_CO  => POINT_TO_COMPACT_VECTOR(NVECTOR_CO, NM, NL_CO)

               NC_CO => SCARC(NM)%COMPACT(NL_CO)%NC

               R     => SCARC(NM)%COMPACT(NL_CO)%R
               R_ROW => SCARC(NM)%COMPACT(NL_CO)%R_ROW
               R_COL => SCARC(NM)%COMPACT(NL_CO)%R_COL
   
               DO IC_CO = 1, NC_CO
                  AUX = 0.0_EB
                  DO ICOL = R_ROW(IC_CO), R_ROW(IC_CO+1)-1
                     IC = R_COL(ICOL)
                     AUX = DC_FI(IC) * R(ICOL)
                  ENDDO
                  FC_CO(IC_CO) = AUX
               ENDDO

            ENDDO

      END SELECT SELECT_COMPACT_MULTIGRID

END SELECT SELECT_SYSTEM

END SUBROUTINE SCARC_RESTRICTION


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform prolongation from coarser to finer grid in multigrid method
!!!    - 'CO' corresponds to coarser grid
!!!    - 'FI' corresponds to finer   grid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PROLONGATION (NVECTOR_CO, NVECTOR_FI, NL_CO, NL_FI)
INTEGER, INTENT(IN) :: NVECTOR_CO, NVECTOR_FI, NL_CO, NL_FI
INTEGER :: NM, ICOL, IC, I
REAL(EB), POINTER, DIMENSION(:,:,:) :: XB_CO, DB_FI
REAL(EB), POINTER, DIMENSION(:)     :: XC_CO, DC_FI, P
INTEGER , POINTER, DIMENSION(:)     :: P_ROW, P_COL
INTEGER , POINTER :: NX_CO, NY_CO, NZ_CO, NC_FI
INTEGER  :: NX_FI, NY_FI, NZ_FI
INTEGER  :: IX_FI, IY_FI, IZ_FI, IC_FI(8)
INTEGER  :: IX_CO, IY_CO, IZ_CO, IC_CO
REAL(EB) :: AUX

SELECT_SYSTEM: SELECT CASE (TYPE_SYSTEM)

   !!!----------------------------------------------------------------------------------------------------
   !!! Banded system
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_SYSTEM_BANDED)

      SELECT_BANDED_MULTIGRID: SELECT CASE (TYPE_MULTIGRID)

         !!!------------------------- Geometric multigrid -----------------------------------------------
         CASE (NSCARC_MULTIGRID_GEOMETRIC)

            DO NM = NMESHES_MIN, NMESHES_MAX
            
               NX_CO => SCARC(NM)%BANDED(NL_CO)%NX
               NY_CO => SCARC(NM)%BANDED(NL_CO)%NY
               NZ_CO => SCARC(NM)%BANDED(NL_CO)%NZ
            
               XB_CO  => POINT_TO_BANDED_VECTOR(NVECTOR_CO, NM, NL_CO)
               DB_FI  => POINT_TO_BANDED_VECTOR(NVECTOR_FI, NM, NL_FI)
      
               SELECT_BANDED_DIMENSION: SELECT CASE (TYPE_DIMENSION)

                  !!!----------------- 2D -------------------
                  CASE (NSCARC_DIMENSION_TWO)
               
                     DO IZ_CO = 1, NZ_CO
                        DO IX_CO = 1, NX_CO
               
                           IX_FI = 2*IX_CO
                           IY_FI = 1
                           IZ_FI = 2*IZ_CO
               
                           DB_FI(IX_FI-1, 1, IZ_FI-1) = XB_CO(IX_CO, 1, IZ_CO)
                           DB_FI(IX_FI  , 1, IZ_FI-1) = XB_CO(IX_CO, 1, IZ_CO)
                           DB_FI(IX_FI-1, 1, IZ_FI  ) = XB_CO(IX_CO, 1, IZ_CO)
                           DB_FI(IX_FI  , 1, IZ_FI  ) = XB_CO(IX_CO, 1, IZ_CO)
               
                        ENDDO
                     ENDDO
               
                  !!!----------------- 3D -------------------
                  CASE (NSCARC_DIMENSION_THREE)
                     
                     DO IZ_CO = 1, NZ_CO
                        DO IY_CO = 1, NY_CO
                           DO IX_CO = 1, NX_CO
               
                              IX_FI = 2*IX_CO
                              IY_FI = 2*IZ_CO
                              IZ_FI = 2*IZ_CO
               
                              DB_FI(IX_FI-1, IY_FI-1, IZ_FI-1) = XB_CO(IX_CO, IY_CO, IZ_CO)
                              DB_FI(IX_FI  , IY_FI-1, IZ_FI-1) = XB_CO(IX_CO, IY_CO, IZ_CO)
                              DB_FI(IX_FI-1, IY_FI  , IZ_FI-1) = XB_CO(IX_CO, IY_CO, IZ_CO)
                              DB_FI(IX_FI  , IY_FI  , IZ_FI-1) = XB_CO(IX_CO, IY_CO, IZ_CO)
                              DB_FI(IX_FI-1, IY_FI-1, IZ_FI  ) = XB_CO(IX_CO, IY_CO, IZ_CO)
                              DB_FI(IX_FI  , IY_FI-1, IZ_FI  ) = XB_CO(IX_CO, IY_CO, IZ_CO)
                              DB_FI(IX_FI-1, IY_FI  , IZ_FI  ) = XB_CO(IX_CO, IY_CO, IZ_CO)
                              DB_FI(IX_FI  , IY_FI  , IZ_FI  ) = XB_CO(IX_CO, IY_CO, IZ_CO)
               
                           ENDDO
                        ENDDO
                     ENDDO
                        
               END SELECT SELECT_BANDED_DIMENSION

            ENDDO
      
         !!!------------------------- Geometric multigrid -----------------------------------------------
         CASE (NSCARC_MULTIGRID_ALGEBRAIC)
            WRITE(*,*) 'No banded restriction for algebraic multigrid available, stopping program!'
            STOP

      END SELECT SELECT_BANDED_MULTIGRID

      
   !!!----------------------------------------------------------------------------------------------------
   !!! Compact system
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_SYSTEM_COMPACT)

      SELECT_COMPACT_MULTIGRID: SELECT CASE (TYPE_MULTIGRID)

         !!!------------------------- Geometric multigrid -----------------------------------------------
         CASE (NSCARC_MULTIGRID_GEOMETRIC)

            DO NM = NMESHES_MIN, NMESHES_MAX
         
               NX_CO => SCARC(NM)%COMPACT(NL_CO)%NX
               NY_CO => SCARC(NM)%COMPACT(NL_CO)%NY
               NZ_CO => SCARC(NM)%COMPACT(NL_CO)%NZ
         
               NX_FI = 2*NX_CO
               NY_FI = 2*NY_CO
               NZ_FI = 2*NZ_CO
   
               XC_CO  => POINT_TO_COMPACT_VECTOR(NVECTOR_CO, NM, NL_CO)
               DC_FI  => POINT_TO_COMPACT_VECTOR(NVECTOR_FI, NM, NL_FI)
   
               SELECT_COMPACT_DIMENSION: SELECT CASE (TYPE_DIMENSION)
               
                  !!!----------------- 2D -------------------
                  CASE (NSCARC_DIMENSION_TWO)
               
                     DO IZ_CO = 1, NZ_CO
                        DO IX_CO = 1, NX_CO
               
                           IX_FI = 2*IX_CO
                           IY_FI = 1
                           IZ_FI = 2*IZ_CO
               
                           IC_CO = (IZ_CO-1)*NX_CO + IX_CO
               
                           IC_FI(1) = (IZ_FI-2)*NX_FI + IX_FI - 1
                           IC_FI(2) = (IZ_FI-2)*NX_FI + IX_FI   
                           IC_FI(3) = (IZ_FI-1)*NX_FI + IX_FI - 1
                           IC_FI(4) = (IZ_FI-1)*NX_FI + IX_FI    
               
                           DO I = 1, 4
                              DC_FI(IC_FI(I)) = XC_CO(IC_CO)
                           ENDDO
               
                        ENDDO
                     ENDDO
               
                              
                  !!!----------------- 3D -------------------
                  CASE (NSCARC_DIMENSION_THREE)
                     
                     DO IZ_CO = 1, NZ_CO
                        DO IY_CO = 1, NY_CO
                           DO IX_CO = 1, NX_CO
               
                              IX_FI = 2*IX_CO
                              IY_FI = 2*IY_CO
                              IZ_FI = 2*IZ_CO
               
                              IC_CO    = (IZ_CO-1)*NX_CO*NY_CO + (IY_CO-1)*NX_CO + IX_CO
               
                              IC_FI(1) = (IZ_FI-2)*NX_FI*NY_FI + (IY_FI-2)*NX_FI + IX_FI - 1
                              IC_FI(2) = (IZ_FI-2)*NX_FI*NY_FI + (IY_FI-2)*NX_FI + IX_FI    
                              IC_FI(3) = (IZ_FI-2)*NX_FI*NY_FI + (IY_FI-1)*NX_FI + IX_FI - 1
                              IC_FI(4) = (IZ_FI-2)*NX_FI*NY_FI + (IY_FI-1)*NX_FI + IX_FI    
                              IC_FI(5) = (IZ_FI-1)*NX_FI*NY_FI + (IY_FI-2)*NX_FI + IX_FI - 1
                              IC_FI(6) = (IZ_FI-1)*NX_FI*NY_FI + (IY_FI-2)*NX_FI + IX_FI    
                              IC_FI(7) = (IZ_FI-1)*NX_FI*NY_FI + (IY_FI-1)*NX_FI + IX_FI - 1
                              IC_FI(8) = (IZ_FI-1)*NX_FI*NY_FI + (IY_FI-1)*NX_FI + IX_FI    
               
                              DO I = 1, 8
                                 DC_FI(IC_FI(I)) = XC_CO(IC_CO)
                              ENDDO
               
                           ENDDO
                        ENDDO
                     ENDDO
   
               END SELECT SELECT_COMPACT_DIMENSION

            ENDDO

         !!!------------------------- Algebraic multigrid -----------------------------------------------
         CASE (NSCARC_MULTIGRID_ALGEBRAIC)

            DO NM = NMESHES_MIN, NMESHES_MAX

               XC_CO  => POINT_TO_COMPACT_VECTOR(NVECTOR_CO, NM, NL_CO)
               DC_FI  => POINT_TO_COMPACT_VECTOR(NVECTOR_FI, NM, NL_FI)

               NC_FI => SCARC(NM)%COMPACT(NL_FI)%NC

               P     => SCARC(NM)%COMPACT(NL_FI)%P
               P_ROW => SCARC(NM)%COMPACT(NL_FI)%P_ROW
               P_COL => SCARC(NM)%COMPACT(NL_FI)%P_COL
   
               DO IC = 1, NC_FI
                  AUX = 0.0_EB
                  DO ICOL = P_ROW(IC), P_ROW(IC+1)-1
                     IC_CO = P_COL(ICOL)
                     AUX = XC_CO(IC_CO) * P(ICOL)
                  ENDDO
                  DC_FI(IC) = AUX
               ENDDO

            ENDDO
            
      END SELECT SELECT_COMPACT_MULTIGRID

END SELECT SELECT_SYSTEM

END SUBROUTINE SCARC_PROLONGATION


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Finalize data - banded storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_TERMINATE_SOLVER(NL)
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, I, J, K, IC
REAL(EB), POINTER, DIMENSION(:,:,:) :: HP
TYPE (MESH_TYPE), POINTER :: M
TYPE (SCARC_SYSTEM_BANDED_TYPE) , POINTER :: SB
TYPE (SCARC_SYSTEM_COMPACT_TYPE), POINTER :: SC

SELECT CASE (TYPE_SYSTEM)
 
   !!!----------------------------------------------------------------------------------------------------
   !!! Banded system
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_SYSTEM_BANDED)

      DO NM = NMESHES_MIN, NMESHES_MAX
      
         M  => MESHES(NM)
         SB => SCARC(NM)%BANDED(NL)
      
         IF (PREDICTOR) THEN
            HP => M%H
         ELSE
            HP => M%HS
         ENDIF
      
         !!! Overwrite internal values of H or HS by corresponding data of X
         HP(1:M%IBAR, 1:M%JBAR, 1:M%KBAR) = SB%X(1:M%IBAR, 1:M%JBAR, 1:M%KBAR)
         
      ENDDO 
      
   !!!----------------------------------------------------------------------------------------------------
   !!! Compact system
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_SYSTEM_COMPACT)

      DO NM = NMESHES_MIN, NMESHES_MAX
      
         M  => MESHES(NM)
         SC => SCARC(NM)%COMPACT(NL)
      
         IF (PREDICTOR) THEN
            HP => M%H
         ELSE
            HP => M%HS
         ENDIF

         DO K = 1, M%KBAR
            DO J = 1, M%JBAR
               DO I = 1, M%IBAR
                  IC = (K-1) * M%IBAR * M%JBAR + (J-1) * M%IBAR + I
                  HP(I, J, K) = SC%X(IC)
               ENDDO
            ENDDO
         ENDDO
      
      ENDDO

END SELECT

END SUBROUTINE SCARC_TERMINATE_SOLVER


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Set correct boundary values at external and internal boundaries
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_UPDATE_GHOSTCELLS(NL)
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IW, IOR0, I0, J0, K0, I1, J1, K1
REAL(EB), POINTER, DIMENSION(:,:,:) :: HP=>NULL()
TYPE (MESH_TYPE), POINTER :: M


!!!---------------------------------------------------------------------------------------------------
!!! Adjust ghost values along external boundaries according to boundary arrays BXS, BXF, ...
!!!---------------------------------------------------------------------------------------------------
DO NM = NMESHES_MIN, NMESHES_MAX

   !!! point to correct pressure vector on mesh 'NM'
   M => MESHES(NM)
   IF (PREDICTOR) THEN
      HP => M%H
   ELSE
      HP => M%HS
   ENDIF
   
   !!! compute ghost cell values
   WALL_CELL_LOOP: DO IW = 1, M%N_EXTERNAL_WALL_CELLS
   
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
   ENDDO WALL_CELL_LOOP
   
ENDDO

!!!---------------------------------------------------------------------------------------------------
!!! Perform data exchange to achieve consistency of ghost values along internal boundaries 
!!!---------------------------------------------------------------------------------------------------
IF (NMESHES > 1) THEN
   NREQ_SCARC = 0
   TYPE_EXCHANGE = NSCARC_EXCHANGE_BDRY
   CALL SCARC_RECEIVE  (NL)
   CALL SCARC_EXCHANGE (NL)
ENDIF

END SUBROUTINE SCARC_UPDATE_GHOSTCELLS

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Receive data from neighbors (corresponds to POST_RECEIVES)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_RECEIVE (NL)

INTEGER, INTENT(IN) :: NL
INTEGER :: NM, NOM, IERR, ILEN
TYPE (SCARC_TYPE)          , POINTER ::  S
TYPE (SCARC_MESH_TYPE)     , POINTER ::  SM
TYPE (OSCARC_TYPE)         , POINTER :: OS
TYPE (OSCARC_MESH_TYPE)    , POINTER :: OSM

IERR=0
RECEIVE_MESH_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX
   RECEIVE_OMESH_LOOP: DO NOM = 1, NMESHES
    
      SNODE = PROCESS(NOM)
      IF (PROCESS(NM)==SNODE) CYCLE RECEIVE_OMESH_LOOP

      S   => SCARC(NM)                            ! corresponds to M
      OS  => SCARC(NM)%OSCARC(NOM)                ! corresponds to M3
      SM  => SCARC(NM)%MESHES(NL)                   ! corresponds to M  for the level 'NL'

      IF (OS%NICMAX_S==0 .AND. OS%NICMAX_R==0) CYCLE RECEIVE_OMESH_LOOP
      OSM => SCARC(NM)%OSCARC(NOM)%MESHES(NL)       ! corresponds to M3 for the level 'NL'

      SELECT_EXCHANGE_TYPE: SELECT CASE (TYPE_EXCHANGE)

         !!!-------------------------------------------------------------------------------------------
         !!! Initialize communication structures for the receiving of data
         !!!-------------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_INIT)

            IF (NL/=NLEVEL_MIN) THEN  ! for lower levels get neighboring IJKW's
               IF (USE_MPI) THEN
                  NREQ_SCARC = NREQ_SCARC+1
                  CALL MPI_IRECV(OSM%IJKW(1,1),15*OSM%N_EXTERNAL_WALL_CELLS,MPI_INTEGER,SNODE, &
                                 TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)
               ENDIF
            ELSE                                ! on maximum level neighboring IJKW already exists
               OSM%IJKW => MESHES(NM)%OMESH(NOM)%IJKW     
            ENDIF
   
            IF (OS%NIC_R > 0) THEN
               ILEN=(MAX(OS%NIC_R, OS%NIC_S)+2)*2+1
               ALLOCATE (OS%RECV_BUF(ILEN))
               OS%RECV_BUF = 0.0_EB
            ENDIF

         !!!-------------------------------------------------------------------------------------------
         !!! Perform exchanges for 
         !!!    - internal values for matrix-vector multiplication 
         !!!    - internal boundariy values
         !!!    - internal subdiagonal matrix values
         !!!-------------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_VECTOR, NSCARC_EXCHANGE_BDRY, NSCARC_EXCHANGE_MATRIX, NSCARC_EXCHANGE_MEASURE)

            NREQ_SCARC = NREQ_SCARC+1
            IF (USE_MPI) CALL MPI_IRECV(OS%RECV_BUF(1),SIZE(OS%RECV_BUF),MPI_DOUBLE_PRECISION,&
                                        SNODE,TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)


      END SELECT SELECT_EXCHANGE_TYPE
   ENDDO RECEIVE_OMESH_LOOP
ENDDO RECEIVE_MESH_LOOP

END SUBROUTINE SCARC_RECEIVE
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Pack send buffer for exchange of internal vector in banded storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PACK_BANDED_VECTOR (X, SEND_BUF, IJKW, NW, NLEN0, NM)
REAL (EB), DIMENSION (0:,0:,0:), INTENT(IN)  :: X
REAL (EB), DIMENSION (:)       , INTENT(OUT) :: SEND_BUF
INTEGER  , DIMENSION (:,:)     , INTENT(IN)  :: IJKW
INTEGER,  INTENT(IN) :: NW
INTEGER, INTENT(OUT) :: NLEN0
INTEGER, INTENT(IN)  :: NM
INTEGER ::  IW, IWW, LL, II, JJ, KK

LL  = 0
IWW = 0
PACK_SEND: DO IW=1,NW
   IF (IJKW(9,IW)/=NM) CYCLE PACK_SEND
   DO KK=IJKW(12,IW),IJKW(15,IW)
      DO JJ=IJKW(11,IW),IJKW(14,IW)
         DO II=IJKW(10,IW),IJKW(13,IW)
            IWW = IWW + 1
            SEND_BUF(LL+1) = REAL(IW,EB)
            SEND_BUF(LL+2) = X(II,JJ,KK)
            LL = LL+2
         ENDDO
      ENDDO
   ENDDO
ENDDO PACK_SEND
NLEN0=2*IWW+1

SEND_BUF(NLEN0) = -999.0_EB

END SUBROUTINE SCARC_PACK_BANDED_VECTOR
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Pack send buffer for exchange of internal vector in compact storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PACK_COMPACT_VECTOR (X, SEND_BUF, IJKW, NW, NLEN0, IBAR, JBAR, NM)
REAL (EB), DIMENSION (:)  , INTENT(IN)  :: X
REAL (EB), DIMENSION (:)  , INTENT(OUT) :: SEND_BUF
INTEGER  , DIMENSION (:,:), INTENT(IN)  :: IJKW
INTEGER, INTENT(IN)  :: NW, IBAR, JBAR
INTEGER, INTENT(OUT) :: NLEN0
INTEGER, INTENT(IN)  :: NM
INTEGER ::  IW, IWW, LL, II, JJ, KK, IC

LL  = 0
IWW = 0
PACK_SEND: DO IW=1,NW
   IF (IJKW(9,IW)/=NM) CYCLE PACK_SEND
   DO KK=IJKW(12,IW),IJKW(15,IW)
      DO JJ=IJKW(11,IW),IJKW(14,IW)
         DO II=IJKW(10,IW),IJKW(13,IW)
            IWW = IWW + 1
            IC = (KK-1)*IBAR*JBAR + (JJ-1)*IBAR + II
            SEND_BUF(LL+1) = REAL(IW,EB)
            SEND_BUF(LL+2) = X(IC)
            LL = LL+2
         ENDDO
      ENDDO
   ENDDO
ENDDO PACK_SEND
NLEN0=2*IWW+1

SEND_BUF(NLEN0) = -999.0_EB

END SUBROUTINE SCARC_PACK_COMPACT_VECTOR
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Pack send buffer for exchange of internal vector in compact storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PACK_BANDED_MATRIX (SEND_BUF, DI2, IJKW, NW, NLEN0, NM)
REAL (EB), DIMENSION (:)  , INTENT(IN)  :: DI2
REAL (EB), DIMENSION (:)  , INTENT(OUT) :: SEND_BUF
INTEGER  , DIMENSION (:,:), INTENT(IN)  :: IJKW
INTEGER  , INTENT(IN)  :: NW, NM
INTEGER  , INTENT(OUT) :: NLEN0
INTEGER ::  IW, IWW, LL, II, JJ, KK, IOR0

LL = 0
IWW = 0
PACK_SEND: DO IW=1,NW

   IF (IJKW(9,IW)/=NM) CYCLE PACK_SEND
   IOR0 = ABS(IJKW(4,IW))

   DO KK=IJKW(12,IW),IJKW(15,IW)
      DO JJ=IJKW(11,IW),IJKW(14,IW)
         DO II=IJKW(10,IW),IJKW(13,IW)
            IWW = IWW + 1
            SEND_BUF(LL+1) = REAL(IW,EB)
            SEND_BUF(LL+2) = DI2(IOR0)
            LL = LL+2
         ENDDO
      ENDDO
   ENDDO

ENDDO PACK_SEND
NLEN0=2*IWW+1

SEND_BUF(NLEN0) = -999.0_EB

END SUBROUTINE SCARC_PACK_BANDED_MATRIX
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Unpack receive buffer for internal update
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_UNPACK_BANDED_VECTOR (X, RECV_BUF, IJKW)
REAL (EB), DIMENSION (0:,0:,0:), INTENT(INOUT) :: X
REAL (EB), DIMENSION (:)    , INTENT(IN)  :: RECV_BUF
INTEGER  , DIMENSION (:,:)  , INTENT(IN)  :: IJKW
REAL (EB):: ZSUM=0.0_EB
INTEGER IW, LL, ISUM, I, J, K, II, JJ, KK

LL = 0
UNPACK_RECV2: DO

   IW = NINT(RECV_BUF(LL+1))
   IF (IW==-999) EXIT UNPACK_RECV2
   ZSUM=0.0_EB
   DO KK=IJKW(12,IW),IJKW(15,IW)
      DO JJ=IJKW(11,IW),IJKW(14,IW)
         DO II=IJKW(10,IW),IJKW(13,IW)
            ZSUM=ZSUM+RECV_BUF(LL+2)
            LL = LL+2
         ENDDO
      ENDDO
   ENDDO

   ISUM = (IJKW(13,IW)-IJKW(10,IW)+1) * &
          (IJKW(14,IW)-IJKW(11,IW)+1) * &
          (IJKW(15,IW)-IJKW(12,IW)+1)

   I=IJKW(1,IW)
   J=IJKW(2,IW)
   K=IJKW(3,IW)
   X(I, J, K) = ZSUM/REAL(ISUM,EB)

ENDDO UNPACK_RECV2

END SUBROUTINE SCARC_UNPACK_BANDED_VECTOR


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Unpack receive buffer for internal update
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_UNPACK_COMPACT_VECTOR (X, RECV_BUF, IJKW, ADJACENT_CELL)
REAL (EB), DIMENSION (:)  , INTENT(OUT):: X
REAL (EB), DIMENSION (:)  , INTENT(IN) :: RECV_BUF
INTEGER  , DIMENSION (:,:), INTENT(IN) :: IJKW
INTEGER  , DIMENSION (:)  , INTENT(IN) :: ADJACENT_CELL
REAL (EB):: ZSUM
INTEGER IW, LL, ISUM, II, JJ, KK, IC

LL = 0
UNPACK_RECV2: DO

   IW = NINT(RECV_BUF(LL+1))
   IF (IW==-999) EXIT UNPACK_RECV2
   ZSUM=0.0_EB
   DO KK=IJKW(12,IW),IJKW(15,IW)
      DO JJ=IJKW(11,IW),IJKW(14,IW)
         DO II=IJKW(10,IW),IJKW(13,IW)
            ZSUM=ZSUM+RECV_BUF(LL+2)
            LL = LL+2
         ENDDO
      ENDDO
   ENDDO

   ISUM = (IJKW(13,IW)-IJKW(10,IW)+1) * &
          (IJKW(14,IW)-IJKW(11,IW)+1) * &
          (IJKW(15,IW)-IJKW(12,IW)+1)

   IC = ADJACENT_CELL(IW)
   X(IC) = ZSUM/REAL(ISUM,EB)

ENDDO UNPACK_RECV2

END SUBROUTINE SCARC_UNPACK_COMPACT_VECTOR


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Unpack receive buffer for matrix-vector multiplication for compact storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_UNPACK_BANDED_MATRIX (RECV_BUF, A, IJKW, ADJACENT_CELL, IBAR, JBAR, NCPL)
REAL (EB), DIMENSION (:,:), INTENT(OUT) :: A
REAL (EB), DIMENSION (:)  , INTENT(IN)  :: RECV_BUF
INTEGER  , DIMENSION (:,:), INTENT(IN)  :: IJKW
INTEGER  , DIMENSION (:)  , INTENT(OUT) :: ADJACENT_CELL
INTEGER  , INTENT(IN) :: IBAR, JBAR, NCPL
REAL (EB):: ZSUM
INTEGER IW, LL, ISUM, I, J, K, II, JJ, KK, IC, ICPL


LL = 0
UNPACK_RECV1: DO

   IW = NINT(RECV_BUF(LL+1))
   IF (IW==-999) EXIT UNPACK_RECV1
   ZSUM=0.0_EB
   DO KK=IJKW(12,IW),IJKW(15,IW)
      DO JJ=IJKW(11,IW),IJKW(14,IW)
         DO II=IJKW(10,IW),IJKW(13,IW)
            ZSUM=ZSUM+RECV_BUF(LL+2)
            LL = LL+2
         ENDDO
      ENDDO
   ENDDO

   ISUM = (IJKW(13,IW)-IJKW(10,IW)+1) * &
          (IJKW(14,IW)-IJKW(11,IW)+1) * &
          (IJKW(15,IW)-IJKW(12,IW)+1)

   I=IJKW(6,IW)
   J=IJKW(7,IW)
   K=IJKW(8,IW)

   IC = (K-1)*IBAR*JBAR + (J-1)*IBAR + I
   DO ICPL = 2, NCPL
      IF (IW == -A(IC,ICPL)) THEN
         A(IC,ICPL) = ZSUM/REAL(ISUM,EB)
         ADJACENT_CELL (IW) = IC
      ENDIF
   ENDDO

ENDDO UNPACK_RECV1

END SUBROUTINE SCARC_UNPACK_BANDED_MATRIX


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Unpack receive buffer for matrix-vector multiplication for compact storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_UNPACK_COMPACT_MATRIX (RECV_BUF, A, A_ROW, A_COL, IJKW, ADJACENT_CELL, IBAR, JBAR, NCE)
REAL (EB), DIMENSION (:)  , INTENT(OUT)   :: A
INTEGER  , DIMENSION (:)  , INTENT(IN)    :: A_ROW
INTEGER  , DIMENSION (:)  , INTENT(INOUT) :: A_COL
REAL (EB), DIMENSION (:)  , INTENT(IN)    :: RECV_BUF
INTEGER  , DIMENSION (:,:), INTENT(IN)    :: IJKW
INTEGER  , DIMENSION (:)  , INTENT(OUT)   :: ADJACENT_CELL
INTEGER  , INTENT(IN)    :: IBAR, JBAR
INTEGER  , INTENT(INOUT) :: NCE
REAL (EB):: ZSUM
INTEGER IW, LL, ISUM, I, J, K, II, JJ, KK, IC, IROW, ICOL

LL = 0
UNPACK_RECV1: DO

   IW = NINT(RECV_BUF(LL+1))
   IF (IW==-999) EXIT UNPACK_RECV1
   ZSUM=0.0_EB
   DO KK=IJKW(12,IW),IJKW(15,IW)
      DO JJ=IJKW(11,IW),IJKW(14,IW)
         DO II=IJKW(10,IW),IJKW(13,IW)
            ZSUM=ZSUM+RECV_BUF(LL+2)
            LL = LL+2
         ENDDO
      ENDDO
   ENDDO

   ISUM = (IJKW(13,IW)-IJKW(10,IW)+1) * &
          (IJKW(14,IW)-IJKW(11,IW)+1) * &
          (IJKW(15,IW)-IJKW(12,IW)+1)

   I=IJKW(6,IW)
   J=IJKW(7,IW)
   K=IJKW(8,IW)

   IC   = (K-1)*IBAR*JBAR + (J-1)*IBAR + I
   IROW = A_ROW(IC)
   DO ICOL = A_ROW(IC)+1, A_ROW(IC+1)-1
      IF (IW == -A_COL(ICOL)) THEN
         A(ICOL)    = ZSUM/REAL(ISUM,EB)
         NCE        = NCE + 1
         A_COL(ICOL) = NCE
         ADJACENT_CELL (IW) = NCE
      ENDIF
   ENDDO

ENDDO UNPACK_RECV1

END SUBROUTINE SCARC_UNPACK_COMPACT_MATRIX


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Send data to neighbors (corresponds to MESH_EXCHANGE)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_EXCHANGE (NL)
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, NOM
INTEGER :: IERR
INTEGER :: ILEN, NLEN
REAL(EB), POINTER, DIMENSION(:)    :: BUFFER0, VECTORC
REAL(EB), POINTER, DIMENSION(:,:,:):: VECTORB
TYPE ( SCARC_TYPE)     , POINTER :: S 
TYPE (OSCARC_TYPE)     , POINTER :: OS, OSO
TYPE ( SCARC_MESH_TYPE), POINTER :: SM, SOM
TYPE (OSCARC_MESH_TYPE), POINTER :: OSM
TYPE (SCARC_SYSTEM_BANDED_TYPE) , POINTER :: SB, SOB
TYPE (SCARC_SYSTEM_COMPACT_TYPE), POINTER :: SC, SOC

IERR = 0

EXCHANGE_SEND_LOOP1: DO NM = NMESHES_MIN, NMESHES_MAX

   IF (PROCESS(NM)/=MYID)  CYCLE EXCHANGE_SEND_LOOP1

   EXCHANGE_RECV_LOOP1: DO NOM = 1, NMESHES
    
      SNODE = PROCESS(NOM)
      RNODE = PROCESS(NM)

      S   => SCARC(NM)                           ! corresponds to M
      OS  => SCARC(NM)%OSCARC(NOM)               ! corresponds to M3

      IF (OS%NICMAX_S == 0 .AND. OS%NICMAX_R == 0) CYCLE EXCHANGE_RECV_LOOP1

      SM  => S%MESHES(NL)                        ! corresponds to M  for the level 'NL'
      OSM => OS%MESHES(NL)                       ! corresponds to M3 for the level 'NL'


      SELECT_EXCHANGE_TYPE1: SELECT CASE(TYPE_EXCHANGE)

         !!!-------------------------------------------------------------------------------------------
         !!! Initialize the communication structures for sending data
         !!!-------------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_INIT)
 
            ! on max level   : take M%IJKW from neighbors (ScaRC pointers already defined, nothing to do ...)
            ! on lover levels: send own IJKW to neighbors
            IF (RNODE /= SNODE) THEN
               IF (NL /= NLEVEL_MIN) THEN
                  NREQ_SCARC = NREQ_SCARC+1
                  IF (USE_MPI) CALL MPI_ISEND(SM%IJKW(1,1),15*SM%N_EXTERNAL_WALL_CELLS,MPI_INTEGER,SNODE, &
                                              TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)
               ENDIF
            ELSE
               SOM      => SCARC(NOM)%MESHES(NL)
               OSM%IJKW => SOM%IJKW(:,1:SOM%N_EXTERNAL_WALL_CELLS)
            ENDIF

            ILEN=(MAX(OS%NIC_R, OS%NIC_S)+2)*2+1   
            ALLOCATE (OS%SEND_BUF(ILEN))
            OS%SEND_BUF = 0.0_EB


         !!!-------------------------------------------------------------------------------------------
         !!! Exchange vector values along internal boundary
         !!!-------------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_VECTOR) 

            SELECT CASE(TYPE_SYSTEM)
   
               !!! banded storage technique
               CASE (NSCARC_SYSTEM_BANDED)
   
                  SB => S%BANDED(NL)
   
                  SELECT CASE (TYPE_VECTOR)
                     CASE (NSCARC_VECTOR_X)
                        VECTORB => SB%X
                     CASE (NSCARC_VECTOR_Y)
                        VECTORB => SB%Y
                     CASE (NSCARC_VECTOR_G)
                        VECTORB => SB%G
                     CASE (NSCARC_VECTOR_W)
                        VECTORB => SB%W
                     CASE (NSCARC_VECTOR_D)
                        VECTORB => SB%D
                     CASE (NSCARC_VECTOR_Z)
                        VECTORB => SB%Z
                     CASE (NSCARC_VECTOR_X2)
                        VECTORB => SB%X2
                     CASE (NSCARC_VECTOR_D2)
                        VECTORB => SB%D2
                     CASE (NSCARC_VECTOR_W2)
                        VECTORB => SB%W2
                     CASE (NSCARC_VECTOR_Y2)
                        VECTORB => SB%Y2
                  END SELECT
   
                  BUFFER0 => OS%SEND_BUF
                  NLEN    =  0
         
                  CALL SCARC_PACK_BANDED_VECTOR(VECTORB, BUFFER0, OSM%IJKW, OSM%N_EXTERNAL_WALL_CELLS, NLEN, NM)
   
   
               !!! compact storage technique
               CASE (NSCARC_SYSTEM_COMPACT)
   
                  SC => S%COMPACT(NL)
   
                  SELECT CASE (TYPE_VECTOR)
                     CASE (NSCARC_VECTOR_X)
                        VECTORC => SC%X
                     CASE (NSCARC_VECTOR_Y)
                        VECTORC => SC%Y
                     CASE (NSCARC_VECTOR_G)
                        VECTORC => SC%G
                     CASE (NSCARC_VECTOR_W)
                        VECTORC => SC%W
                     CASE (NSCARC_VECTOR_D)
                        VECTORC => SC%D
                     CASE (NSCARC_VECTOR_Z)
                        VECTORC => SC%Z
                     CASE (NSCARC_VECTOR_X2)
                        VECTORC => SC%X2
                     CASE (NSCARC_VECTOR_D2)
                        VECTORC => SC%D2
                     CASE (NSCARC_VECTOR_W2)
                        VECTORC => SC%W2
                     CASE (NSCARC_VECTOR_Y2)
                        VECTORC => SC%Y2
                  END SELECT
   
                  BUFFER0 => OS%SEND_BUF
                  NLEN    =  0

                  CALL SCARC_PACK_COMPACT_VECTOR(VECTORC, BUFFER0, OSM%IJKW, OSM%N_EXTERNAL_WALL_CELLS, NLEN, &
                                          SM%IBAR, SM%JBAR, NM)
   
            END SELECT
   
            IF (RNODE/=SNODE) THEN
               NREQ_SCARC=NREQ_SCARC+1
               IF (USE_MPI) CALL MPI_ISEND(BUFFER0, NLEN, MPI_DOUBLE_PRECISION, SNODE, &
                                           TAG_SCARC, MPI_COMM_WORLD, REQ_SCARC(NREQ_SCARC), IERR)
            ENDIF
   
   
         !!!-------------------------------------------------------------------------------------------
         !!! Exchange internal vector values
         !!!-------------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_BDRY) 

            IF (PREDICTOR) THEN
               VECTORB => MESHES(NM)%H
            ELSE
               VECTORB => MESHES(NM)%HS
            ENDIF

            BUFFER0 => OS%SEND_BUF
            NLEN    =  0

            CALL SCARC_PACK_BANDED_VECTOR(VECTORB, BUFFER0, OSM%IJKW, OSM%N_EXTERNAL_WALL_CELLS, NLEN, NM)
   
            IF (RNODE/=SNODE) THEN
               NREQ_SCARC=NREQ_SCARC+1
               IF (USE_MPI) CALL MPI_ISEND(BUFFER0, NLEN, MPI_DOUBLE_PRECISION, SNODE, &
                                           TAG_SCARC, MPI_COMM_WORLD, REQ_SCARC(NREQ_SCARC), IERR)
            ENDIF
   

         !!!-------------------------------------------------------------------------------------------
         !!! Exchange internal matrix values
         !!!-------------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MATRIX) 

            BUFFER0 => OS%SEND_BUF
            NLEN    =  0
          
            CALL SCARC_PACK_BANDED_MATRIX(BUFFER0, SM%DI2, OSM%IJKW, OSM%N_EXTERNAL_WALL_CELLS, NLEN, NM)
   
            IF (RNODE/=SNODE) THEN
               NREQ_SCARC=NREQ_SCARC+1
               IF (USE_MPI) CALL MPI_ISEND(BUFFER0, NLEN, MPI_DOUBLE_PRECISION, SNODE, &
                                           TAG_SCARC, MPI_COMM_WORLD, REQ_SCARC(NREQ_SCARC), IERR)
            ENDIF
   
         !!!-------------------------------------------------------------------------------------------
         !!! Exchange internal measure values
         !!!-------------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MEASURE) 

            VECTORC => S%COMPACT(NL)%MEASURE
            BUFFER0 => OS%SEND_BUF
            NLEN    =  0

            CALL SCARC_PACK_COMPACT_VECTOR(VECTORC, BUFFER0, OSM%IJKW, OSM%N_EXTERNAL_WALL_CELLS, NLEN, &
                                           SM%IBAR, SM%JBAR, NM)
   
            IF (RNODE/=SNODE) THEN
               NREQ_SCARC=NREQ_SCARC+1
               IF (USE_MPI) CALL MPI_ISEND(BUFFER0, NLEN, MPI_DOUBLE_PRECISION, SNODE, &
                                           TAG_SCARC, MPI_COMM_WORLD, REQ_SCARC(NREQ_SCARC), IERR)
            ENDIF
   
      END SELECT SELECT_EXCHANGE_TYPE1

   ENDDO EXCHANGE_RECV_LOOP1

ENDDO EXCHANGE_SEND_LOOP1


!!!----------------------------------------------------------------------------------------------------
!!! Information from Mesh NM is received by Mesh NOM  (NOM receiver, NM sender)
!!!----------------------------------------------------------------------------------------------------
IF (USE_MPI.AND.NREQ_SCARC/=0) CALL MPI_WAITALL(NREQ_SCARC,REQ_SCARC(1:NREQ_SCARC),MPI_STATUS_IGNORE,IERR)


!!!----------------------------------------------------------------------------------------------------
!!! Extract communication data from corresponding RECEIVE-buffers
!!!----------------------------------------------------------------------------------------------------
EXCHANGE_SEND_LOOP2: DO NOM = NMESHES_MIN, NMESHES_MAX

   EXCHANGE_RECV_LOOP2: DO NM=1,NMESHES
    
      SNODE  = PROCESS(NOM)
      RNODE  = PROCESS(NM)

      OS  => SCARC(NM)%OSCARC(NOM)
      OSO => SCARC(NOM)%OSCARC(NM)

      SOM  => SCARC(NOM)%MESHES(NL)

      EXCHANGE_RECV_IF: IF (OSO%NICMAX_S/=0 .AND. OSO%NICMAX_R/=0) THEN

         IF (RNODE/=SNODE) THEN
            BUFFER0 => OSO%RECV_BUF
         ELSE
            BUFFER0 => OS%SEND_BUF
         ENDIF

         SELECT_EXCHANGE_TYPE2: SELECT CASE (TYPE_EXCHANGE)

            !!!----------------------------------------------------------------------------------------
            !!! Extract data for exchange of internal boundary values
            !!!----------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_VECTOR) 

               !!! banded storage technique
               SELECT CASE (TYPE_SYSTEM)
   
                  CASE (NSCARC_SYSTEM_BANDED)
   
                     SOB => SCARC(NOM)%BANDED(NL)
   
                     SELECT CASE (TYPE_VECTOR)
                        CASE (NSCARC_VECTOR_X)
                           VECTORB => SOB%X
                        CASE (NSCARC_VECTOR_Y)
                           VECTORB => SOB%Y
                        CASE (NSCARC_VECTOR_G)
                           VECTORB => SOB%G
                        CASE (NSCARC_VECTOR_W)
                           VECTORB => SOB%W
                        CASE (NSCARC_VECTOR_D)
                           VECTORB => SOB%D
                        CASE (NSCARC_VECTOR_Z)
                           VECTORB => SOB%Z
                        CASE (NSCARC_VECTOR_X2)
                           VECTORB => SOB%X2
                        CASE (NSCARC_VECTOR_D2)
                           VECTORB => SOB%D2
                        CASE (NSCARC_VECTOR_W2)
                           VECTORB => SOB%W2
                        CASE (NSCARC_VECTOR_Y2)
                           VECTORB => SOB%Y2
                     END SELECT
   
                     CALL SCARC_UNPACK_BANDED_VECTOR(VECTORB, BUFFER0, SOM%IJKW)
   
                  !!! compact storage technique
                  CASE (NSCARC_SYSTEM_COMPACT)
   
                     SOC => SCARC(NOM)%COMPACT(NL)
   
                     SELECT CASE (TYPE_VECTOR)
                        CASE (NSCARC_VECTOR_X)
                           VECTORC => SOC%X
                        CASE (NSCARC_VECTOR_Y)
                           VECTORC => SOC%Y
                        CASE (NSCARC_VECTOR_G)
                           VECTORC => SOC%G
                        CASE (NSCARC_VECTOR_W)
                           VECTORC => SOC%W
                        CASE (NSCARC_VECTOR_D)
                           VECTORC => SOC%D
                        CASE (NSCARC_VECTOR_Z)
                           VECTORC => SOC%Z
                        CASE (NSCARC_VECTOR_X2)
                           VECTORC => SOC%X2
                        CASE (NSCARC_VECTOR_D2)
                           VECTORC => SOC%D2
                        CASE (NSCARC_VECTOR_W2)
                           VECTORC => SOC%W2
                        CASE (NSCARC_VECTOR_Y2)
                           VECTORC => SOC%Y2
                     END SELECT
   
                     CALL SCARC_UNPACK_COMPACT_VECTOR(VECTORC, BUFFER0, SOM%IJKW, SOM%ADJACENT_CELL)

               END SELECT
   

            !!!----------------------------------------------------------------------------------------
            !!! Extract data for exchange of internal boundary values
            !!!----------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_BDRY)

               IF (PREDICTOR) THEN
                  VECTORB => MESHES(NOM)%H
               ELSE
                  VECTORB => MESHES(NOM)%HS
               ENDIF

               CALL SCARC_UNPACK_BANDED_VECTOR(VECTORB, BUFFER0, SOM%IJKW)


            !!!----------------------------------------------------------------------------------------
            !!! Extract data for exchange of internal subdiagonal matrix values
            !!!----------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_MATRIX)

               SELECT CASE (TYPE_SYSTEM)
   
                  !!! system in banded storage technique
                  CASE (NSCARC_SYSTEM_BANDED)
   
                     SOB => SCARC(NOM)%BANDED(NL)
                     CALL SCARC_UNPACK_BANDED_MATRIX(BUFFER0, SOB%A, SOM%IJKW, &
                                               SOM%ADJACENT_CELL, SOM%IBAR, SOM%JBAR, SOB%NCPL)

                  !!! system in compact storage technique
                  CASE (NSCARC_SYSTEM_COMPACT)
   
                     SOC => SCARC(NOM)%COMPACT(NL)
                     CALL SCARC_UNPACK_COMPACT_MATRIX(BUFFER0, SOC%A, SOC%A_ROW, SOC%A_COL, SOM%IJKW, &
                                               SOM%ADJACENT_CELL, SOM%IBAR, SOM%JBAR, SOC%NCE)

              END SELECT

            !!!----------------------------------------------------------------------------------------
            !!! Extract data for exchange of internal boundary values
            !!!----------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_MEASURE)

               VECTORC => SCARC(NOM)%COMPACT(NL)%MEASURE
               CALL SCARC_UNPACK_COMPACT_VECTOR(VECTORC, BUFFER0, SOM%IJKW, SOM%ADJACENT_CELL)


          END SELECT SELECT_EXCHANGE_TYPE2

      ENDIF EXCHANGE_RECV_IF
   ENDDO EXCHANGE_RECV_LOOP2
ENDDO EXCHANGE_SEND_LOOP2

END SUBROUTINE SCARC_EXCHANGE
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Check if difference of two values is less than a given tolerance
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
LOGICAL FUNCTION MATCH (VAL1, VAL2)
REAL (EB), INTENT(IN) :: VAL1, VAL2
REAL (EB) :: TOL
TOL = 1.0E-10_EB
MATCH = .FALSE.
IF (Abs(VAL1-VAL2) <= TOL) MATCH = .TRUE.
RETURN
END FUNCTION MATCH
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Only for debugging reasons: print out vector information on level NL 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SHOW (NVECTOR, CROUTINE, CNAME)
INTEGER, INTENT(IN):: NVECTOR
REAL (EB), POINTER, DIMENSION(:,:,:) :: VB
REAL (EB), POINTER, DIMENSION(:)     :: VC
REAL (EB):: VALUES(10)
INTEGER :: NM, II, JJ, KK, IC, NX8, NY8, NZ8
CHARACTER (*), INTENT(IN) :: CROUTINE, CNAME
TYPE (MESH_TYPE), POINTER :: M
 
IF (TYPE_DEBUG < NSCARC_DEBUG_LESS) RETURN

SELECT_SYSTEM: SELECT CASE (TYPE_SYSTEM)

   !!!----------------------------------------------------------------------------------------------------
   !!! Banded system
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_SYSTEM_BANDED)

      DO NM = NMESHES_MIN, NMESHES_MAX

         M  => MESHES(NM)
         VB => POINT_TO_BANDED_VECTOR (NVECTOR, NM, NLEVEL_MIN)

         NX8=MIN(8,M%IBAR)
         NY8=MIN(8,M%JBAR)
         NZ8=MIN(8,M%KBAR)
         
         WRITE(SCARC_LU,*) '=============================================================================='
         WRITE(SCARC_LU,2000) CROUTINE, CNAME
         WRITE(SCARC_LU,*) '=============================================================================='
         SELECT_BANDED_DIMENSION: SELECT CASE (TYPE_DIMENSION)
      
            CASE (NSCARC_DIMENSION_TWO)
      
               DO KK = NZ8, 1, - 1
                  DO II=1,NX8
                     IF (ABS(VB(II,1,KK))<1.0E-14_EB) THEN
                        VALUES(II)=0.0_EB
                     ELSE
                        VALUES(II)=VB(II,1,KK)
                     ENDIF
                  ENDDO
                  WRITE(SCARC_LU, '(a,i3,a,10e13.5)') 'J= ',JJ,' : ',(VALUES(II), II=1, NX8)
               ENDDO
               WRITE(SCARC_LU,*)  '------------------------------------------------',&
                                  '---------------------------------------------------'
               WRITE(SCARC_LU,1000) ('I = ',II,II=1,NX8)
      
            CASE (NSCARC_DIMENSION_THREE)
               DO KK = NZ8, 1, - 1
                  WRITE(SCARC_LU,'(a,i3,a)')  '----------------------------------------  K = ', KK,&
                                              '----------------------------------------'
                  DO JJ = NY8, 1, - 1
                     DO II=1,NX8
                        IF (ABS(VB(II,JJ,KK))<1.0E-14_EB) THEN
                           VALUES(II)=0.0_EB
                        ELSE
                           VALUES(II)=VB(II,JJ,KK)
                        ENDIF
                     ENDDO
                     WRITE(SCARC_LU, '(a,i3,a,10e13.5)') 'J= ',JJ,' : ',(VALUES(II), II=1, NX8)
                  ENDDO
               ENDDO
               WRITE(SCARC_LU,*)  '------------------------------------------------',&
                                  '---------------------------------------------------'
               !WRITE(SCARC_LU,1000) ('I = ',II,II=1,NX8)
               !WRITE(SCARC_LU,*)  '--------------------------------------------------',&
               !                   '-------------------------------------------------'
               !WRITE(SCARC_LU,*)
         END SELECT SELECT_BANDED_DIMENSION
      
      ENDDO
          
   !!!----------------------------------------------------------------------------------------------------
   !!! Compact system
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_SYSTEM_COMPACT)

      DO NM = NMESHES_MIN, NMESHES_MAX
   
         M  => MESHES(NM)
         VC => POINT_TO_COMPACT_VECTOR (NVECTOR, NM, NLEVEL_MIN)
         
         NX8=MIN(8,M%IBAR)
         NY8=MIN(8,M%JBAR)
         NZ8=MIN(8,M%KBAR)
         
         WRITE(SCARC_LU,*) '================================================================================'
         WRITE(SCARC_LU,2000) CROUTINE, CNAME
         WRITE(SCARC_LU,*) '================================================================================'
         SELECT_COMPACT_DIMENSION: SELECT CASE (TYPE_DIMENSION)
         
            CASE (NSCARC_DIMENSION_TWO)
               DO KK = NZ8, 1, - 1
                  DO II=1,NX8
                     IC = (KK-1)*M%IBAR + II
                     IF (ABS(VC(IC))<1.0E-14_EB) THEN
                        VALUES(II)=0.0_EB
                     ELSE
                        VALUES(II)=VC(IC)
                     ENDIF
                  ENDDO
                  WRITE(SCARC_LU, '(a,i3,a,8e13.5)') 'J= ',1,' : ',(VALUES(II), II=1, NX8)
               ENDDO
               WRITE(SCARC_LU,*)  '------------------------------------------------',&
                                  '---------------------------------------------------'
               WRITE(SCARC_LU,1000) ('I = ',II,II=1,NX8)
         
            CASE (NSCARC_DIMENSION_THREE)
               DO KK = NZ8, 1, - 1
                  WRITE(SCARC_LU,'(a,i3,a)')  '----------------------------------------  K = ', KK,&
                                              '----------------------------------------'
                  DO JJ = NY8, 1, - 1
                     DO II=1,NX8
                        IC = (KK-1)*M%IBAR*M%JBAR + (JJ-1)*M%IBAR + II
                        IF (ABS(VC(IC))<1.0E-14_EB) THEN
                           VALUES(II)=0.0_EB
                        ELSE
                           VALUES(II)=VC(IC)
                        ENDIF
                     ENDDO
                     WRITE(SCARC_LU, '(a,i3,a,8e13.5)') 'J= ',JJ,' : ',(VALUES(II), II=1, NX8)
                  ENDDO
               ENDDO
               WRITE(SCARC_LU,*)  '------------------------------------------------',&
                                  '---------------------------------------------------'
               !WRITE(SCARC_LU,1000) ('I = ',II,II=1,NX8)
               !WRITE(SCARC_LU,*)  '--------------------------------------------------',&
               !                   '-------------------------------------------------'
               !WRITE(SCARC_LU,*)
      END SELECT SELECT_COMPACT_DIMENSION

   ENDDO

END SELECT SELECT_SYSTEM
      
1000 FORMAT(13x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x)
2000 FORMAT('=== ',A,' : ',A,' on mesh ',I4,' on level ',I4) 
END SUBROUTINE SCARC_SHOW


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Only for debugging reasons: print out vector information on level NL 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SHOW_LEVEL (NVECTOR, CROUTINE, CNAME, NL)
INTEGER, INTENT(IN):: NVECTOR, NL
REAL (EB), POINTER, DIMENSION(:,:,:) :: VB
REAL (EB), POINTER, DIMENSION(:)     :: VC
REAL (EB):: VALUES(10)
INTEGER :: NM, II, JJ, KK, IC, NX8, NY8, NZ8
CHARACTER (*), INTENT(IN) :: CROUTINE, CNAME
TYPE (SCARC_SYSTEM_BANDED_TYPE) , POINTER :: SB
TYPE (SCARC_SYSTEM_COMPACT_TYPE), POINTER :: SC
 
IF (TYPE_DEBUG < NSCARC_DEBUG_MEDIUM) RETURN

SELECT_SYSTEM: SELECT CASE (TYPE_SYSTEM)

   !!!----------------------------------------------------------------------------------------------------
   !!! Banded system
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_SYSTEM_BANDED)

      DO NM = NMESHES_MIN, NMESHES_MAX

         SB => SCARC(NM)%BANDED(NL)
         VB => POINT_TO_BANDED_VECTOR (NVECTOR, NM, NL)

         NX8=MIN(8,SB%NX)
         NY8=MIN(8,SB%NY)
         NZ8=MIN(8,SB%NZ)
         
         WRITE(SCARC_LU,*) '=============================================================================='
         WRITE(SCARC_LU,2000) CROUTINE, CNAME, NM, NL, NX8, NY8, NZ8
         WRITE(SCARC_LU,*) '=============================================================================='
         SELECT_BANDED_DIMENSION: SELECT CASE (TYPE_DIMENSION)
         
            CASE (NSCARC_DIMENSION_TWO)
               DO KK = NZ8, 1, - 1
                  DO II=1,NX8
                     IF (ABS(VB(II,1,KK))<1.0E-14_EB) THEN
                        VALUES(II)=0.0_EB
                     ELSE
                        VALUES(II)=VB(II,1,KK)
                     ENDIF
                  ENDDO
                  WRITE(SCARC_LU, '(a,i3,a,10e13.5)') 'J= ',1,' : ',(VALUES(II), II=1, NX8)
               ENDDO
               WRITE(SCARC_LU,*)  '------------------------------------------------',&
                                  '---------------------------------------------------'
               !WRITE(SCARC_LU,1000) ('I = ',II,II=1,NX8)
      
            CASE (NSCARC_DIMENSION_THREE)
               DO KK = NZ8, 1, - 1
                  WRITE(SCARC_LU,'(a,i3,a)')  '----------------------------------------  K = ', KK,&
                                              '----------------------------------------'
                  DO JJ = NY8, 1, - 1
                     DO II=1,NX8
                        IF (ABS(VB(II,JJ,KK))<1.0E-14_EB) THEN
                           VALUES(II)=0.0_EB
                        ELSE
                           VALUES(II)=VB(II,JJ,KK)
                        ENDIF
                     ENDDO
                     WRITE(SCARC_LU, '(a,i3,a,10e13.5)') 'J= ',JJ,' : ',(VALUES(II), II=1, NX8)
                  ENDDO
               ENDDO
               WRITE(SCARC_LU,*)  '------------------------------------------------',&
                                  '---------------------------------------------------'
               !WRITE(SCARC_LU,1000) ('I = ',II,II=1,NX8)
               !WRITE(SCARC_LU,*)  '--------------------------------------------------',&
               !                   '-------------------------------------------------'
               !WRITE(SCARC_LU,*)

         END SELECT SELECT_BANDED_DIMENSION

      ENDDO
         
   !!!----------------------------------------------------------------------------------------------------
   !!! Compact system
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_SYSTEM_COMPACT)

      DO NM = NMESHES_MIN, NMESHES_MAX

         SC => SCARC(NM)%COMPACT(NL)
         VC => POINT_TO_COMPACT_VECTOR (NVECTOR, NM, NL)

         NX8=MIN(8,SC%NX)
         NY8=MIN(8,SC%NY)
         NZ8=MIN(8,SC%NZ)
         
         WRITE(SCARC_LU,*) '=============================================================================='
         WRITE(SCARC_LU,2000) CROUTINE, CNAME, NM, NL, NX8, NY8, NZ8
         WRITE(SCARC_LU,*) '=============================================================================='
         SELECT_COMPACT_DIMENSION: SELECT CASE (TYPE_DIMENSION)
            
            CASE (NSCARC_DIMENSION_TWO)
               DO KK = NZ8, 1, - 1
                  DO II=1,NX8
                     IC = (KK-1)*SC%NX + II
                     IF (ABS(VC(IC))<1.0E-14_EB) THEN
                        VALUES(II)=0.0_EB
                     ELSE
                        VALUES(II)=VC(IC)
                     ENDIF
                  ENDDO
                  WRITE(SCARC_LU, '(a,i3,a,8e13.5)') 'J= ',1,' : ',(VALUES(II), II=1, NX8)
               ENDDO
               WRITE(SCARC_LU,*)  '------------------------------------------------',&
                                  '---------------------------------------------------'
               !WRITE(SCARC_LU,1000) ('I = ',II,II=1,NX8)
         
            CASE (NSCARC_DIMENSION_THREE)
               DO KK = NZ8, 1, - 1
                  WRITE(SCARC_LU,'(a,i3,a)')  '----------------------------------------  K = ', KK,&
                                              '----------------------------------------'
                  DO JJ = NY8, 1, - 1
                     DO II=1,NX8
                        IC = (KK-1)*SC%NX*SC%NY + (JJ-1)*SC%NX + II
                        IF (ABS(VC(IC))<1.0E-14_EB) THEN
                           VALUES(II)=0.0_EB
                        ELSE
                           VALUES(II)=VC(IC)
                        ENDIF
                     ENDDO
                     WRITE(SCARC_LU, '(a,i3,a,8e13.5)') 'J= ',JJ,' : ',(VALUES(II), II=1, NX8)
                  ENDDO
               ENDDO
               WRITE(SCARC_LU,*)  '------------------------------------------------',&
                                  '---------------------------------------------------'
               WRITE(SCARC_LU,1000) ('I = ',II,II=1,NX8)
               WRITE(SCARC_LU,*)  '--------------------------------------------------',&
                                  '-------------------------------------------------'
               WRITE(SCARC_LU,*)
            
         END SELECT SELECT_COMPACT_DIMENSION

      ENDDO

END SELECT SELECT_SYSTEM

1000 FORMAT(13x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x,a4,i2,8x)
2000 FORMAT('=== ',A,' : ',A,' on mesh ',I4,' on level ',I4, ': NX, NY, NZ=',3i3)
END SUBROUTINE SCARC_SHOW_LEVEL

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Print out timings for ScaRC - not updated at the moment
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_TIMINGS
INTEGER:: NM, N, IERR, I
INTEGER, ALLOCATABLE, DIMENSION(:) :: COUNTS_SCARC, DISPLS_SCARC, COUNTS_SCARC_TIMERS, DISPLS_SCARC_TIMERS
CHARACTER(40) :: NAME_SCARC(0:N_TIMERS_SCARC)
REAL(EB) :: TPCNT_SCARC(0:N_TIMERS_SCARC)

!!! only temporarily - use routine only in debug mode
IF (TYPE_DEBUG == NSCARC_DEBUG_NONE) RETURN

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

!!! Printout subroutine timings
!!! outdated version, must be revised (not used at the moment)

IF (MYID==0) THEN
   NAME_SCARC                           = 'null'
   NAME_SCARC(NSCARC_TIME_TOTAL)        = 'ScaRC complete'
   NAME_SCARC(NSCARC_TIME_SETUP)        = 'ScaRC setup'
   NAME_SCARC(NSCARC_TIME_SOLVER)       = 'ScaRC solver'
   NAME_SCARC(NSCARC_TIME_KRYLOV)       = 'ScaRC Krylov method'
   NAME_SCARC(NSCARC_TIME_MULTIGRID)    = 'ScaRC multigrid method'
   NAME_SCARC(NSCARC_TIME_PRECON)       = 'ScaRC preconditioner'
   NAME_SCARC(NSCARC_TIME_SMOOTH)       = 'ScaRC smoother'
   NAME_SCARC(NSCARC_TIME_COARSE)       = 'ScaRC coarse grid solver'
   NAME_SCARC(NSCARC_TIME_MATVEC)       = 'ScaRC matrix-vector product'
   NAME_SCARC(NSCARC_TIME_SCALPROD)     = 'ScaRC scalar product'
   NAME_SCARC(NSCARC_TIME_L2NORM)       = 'ScaRC L2-norm'
   NAME_SCARC(NSCARC_TIME_EXCH_INIT)    = 'ScaRC initialization data exchange'
   NAME_SCARC(NSCARC_TIME_EXCH_VECTOR)  = 'ScaRC exchange of vector values'
   NAME_SCARC(NSCARC_TIME_EXCH_MATRIX)  = 'ScaRC exchange of matrix values'
   NAME_SCARC(NSCARC_TIME_EXCH_MEASURE) = 'ScaRC exchange of measure values'
   
   DO NM=1,NMESHES
      DO I=1,N_TIMERS_SCARC
         TPCNT_SCARC(I) = 100._EB*TUSED_SCARC(I,NM)/TUSED(1,NM)
      ENDDO
      WRITE(LU_OUTPUT,443) NM
      !WRITE(LU_OUTPUT,444) 'MAIN',TUSED(1,NM),100._EB
      WRITE(LU_OUTPUT,444) (NAME_SCARC(I),TUSED_SCARC(I,NM),TPCNT_SCARC(I),I=1,N_TIMERS_SCARC)
   ENDDO
ENDIF

   IF (TYPE_DEBUG /= NSCARC_DEBUG_NONE ) CLOSE(SCARC_LU)

443 FORMAT(//' ScaRC: CPU Time Usage, Mesh ',I3// &
         47X,' CPU (s)        %  '/ &
         7X,' -----------------------------------------------------------------')
444 FORMAT(7X,A40,2F11.2)

END SUBROUTINE SCARC_TIMINGS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Stops the code gracefully after writing a message
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SHUTDOWN(CMESSAGE)

CHARACTER(*) CMESSAGE
LOGICAL :: EX

WRITE(LU_ERR,'(/A)') TRIM(CMESSAGE)
INQUIRE(FILE=FN_OUTPUT,EXIST=EX)
IF (.NOT.EX) OPEN(LU_OUTPUT,FILE=TRIM(CHID)//'.out',STATUS='REPLACE',FORM='FORMATTED')
WRITE(LU_OUTPUT,'(/A)') TRIM(CMESSAGE)
STOP

END SUBROUTINE SCARC_SHUTDOWN


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Compute current revision number
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GET_REV_SCRC(MODULE_REV,MODULE_DATE)
INTEGER,INTENT(INOUT) :: MODULE_REV
CHARACTER(255),INTENT(INOUT) :: MODULE_DATE

WRITE(MODULE_DATE,'(A)') SCRCREV(INDEX(SCRCREV,':')+1:LEN_TRIM(SCRCREV)-2)
READ (MODULE_DATE,'(I5)') MODULE_REV
WRITE(MODULE_DATE,'(A)') SCRCDATE

END SUBROUTINE GET_REV_SCRC

END MODULE SCRC
 

